from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.models import SeedConformerRecord
from dsvr.parsing.xtb_outputs import parse_xtb_energy
from dsvr.runners.subprocess_utils import ExternalToolError, run_command


class XtbPrefilterUnavailableError(RuntimeError):
    """Raised when xTB prefilter is enabled but xTB is unavailable."""


@dataclass(frozen=True)
class XtbPrefilterResult:
    seed_id: str
    parent_stereo_id: str | None
    input_molecule_id: str
    energy_kcal_mol: float | None
    workdir: Path
    warnings: list[str]


def run_xtb_prefilter(seed: SeedConformerRecord, config: RunConfig) -> XtbPrefilterResult:
    executable = shutil.which(config.crest.xtb_executable)
    if executable is None:
        raise XtbPrefilterUnavailableError(
            f"xTB prefilter is enabled, but '{config.crest.xtb_executable}' was not found "
            "on PATH. Install xTB or set xtb_prefilter.enabled=false."
        )
    workdir = (
        config.output_dir
        / "xtb_prefilter"
        / seed.input_molecule_id
        / (seed.parent_id or "unknown_stereo")
    )
    workdir.mkdir(parents=True, exist_ok=True)
    input_xyz = workdir / "input.xyz"
    _write_seed_xyz(seed, input_xyz)
    command = [
        executable,
        str(input_xyz.resolve()),
        "--gfn",
        str(config.xtb_prefilter.gfn),
        "--chrg",
        str(seed.formal_charge or 0),
        "-P",
        str(config.xtb_prefilter.nproc),
    ]
    if config.xtb_prefilter.solvent_model != "none":
        command.extend(
            [f"--{config.xtb_prefilter.solvent_model}", config.xtb_prefilter.solvent]
        )
    if config.xtb_prefilter.optimize:
        command.append("--opt")
    warnings: list[str] = []
    try:
        run_command(
            command,
            cwd=workdir,
            timeout_s=config.xtb_prefilter.timeout_seconds_per_variant,
            log_dir=workdir / "logs",
            command_name="xtb_prefilter",
            check=True,
            show_progress=config.logging.tail_subprocess_logs,
        )
    except ExternalToolError as exc:
        warnings.append(f"xTB prefilter failed: {exc}")

    log_path = _first_existing(
        [
            workdir / "logs" / "xtb_prefilter.log",
            workdir / "xtb.out",
            workdir / "output.log",
        ]
    )
    parsed = parse_xtb_energy(log_path.read_text(encoding="utf-8") if log_path else "")
    return XtbPrefilterResult(
        seed_id=seed.id,
        parent_stereo_id=seed.parent_id,
        input_molecule_id=seed.input_molecule_id,
        energy_kcal_mol=getattr(parsed, "electronic_energy_kcal_mol", None),
        workdir=workdir,
        warnings=warnings,
    )


def _write_seed_xyz(seed: SeedConformerRecord, path: Path) -> None:
    if seed.rdkit_mol is None or seed.rdkit_mol.GetNumConformers() == 0:
        raise XtbPrefilterUnavailableError(
            f"Seed {seed.id} has no RDKit conformer coordinates for xTB prefilter."
        )
    molecule = Chem.Mol(seed.rdkit_mol)
    conformer = molecule.GetConformer()
    lines = [str(molecule.GetNumAtoms()), seed.id]
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol()} {position.x:.10f} {position.y:.10f} {position.z:.10f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None
