from __future__ import annotations

import csv
import importlib.util
import json
import re
import shutil
import sys
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord
from dsvr.runners.subprocess_utils import ExternalToolError, run_command
from dsvr.utils.units import hartree_to_kcal_mol

QM_FREE_ENERGY_WARNING = (
    "QM rescoring uses electronic energies only by default; do not interpret these "
    "as high-accuracy free energies unless vibrational/thermal corrections are computed."
)


class Psi4UnavailableError(RuntimeError):
    """Raised when Psi4 rescoring is requested but unavailable."""


class Psi4ExecutionError(RuntimeError):
    """Raised when Psi4 input cannot be prepared."""


def psi4_available() -> bool:
    return shutil.which("psi4") is not None or importlib.util.find_spec("psi4") is not None


def rescore_top_ranked_with_psi4(
    ranked_records: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    if config.refinement.qm_backend != "psi4" and not config.refinement.psi4_enabled:
        return []
    if not psi4_available():
        raise Psi4UnavailableError(
            "Psi4 QM rescoring was requested but Psi4 is unavailable. Install psi4 "
            "or load a module that provides the `psi4` executable, then run `dsvr doctor`."
        )

    selected = _top_n(ranked_records, config)
    rescored: list[RankedVariantRecord] = []
    for candidate in selected:
        workdir = psi4_workdir(candidate, config)
        workdir.mkdir(parents=True, exist_ok=True)
        xyz = _geometry_from_candidate(candidate)
        input_path = workdir / "psi4_input.dat"
        output_path = workdir / "psi4_output.log"
        input_path.write_text(_psi4_input(xyz, candidate, config), encoding="utf-8")
        command = _psi4_command(input_path, output_path, workdir, config)
        warnings = [*candidate.warnings, QM_FREE_ENERGY_WARNING]
        try:
            completed = run_command(
                command,
                cwd=workdir,
                timeout_s=None
                if config.crest.walltime_minutes is None
                else config.crest.walltime_minutes * 60,
                log_dir=workdir / "logs",
                command_name="psi4",
                check=True,
                show_progress=config.logging.tail_subprocess_logs,
            )
            if not output_path.exists():
                output_path.write_text(completed.stdout, encoding="utf-8")
        except ExternalToolError as exc:
            warnings.append(f"Psi4 failed: {exc}")
        energy_h = parse_psi4_energy(output_path)
        rescored.append(_qm_record(candidate, energy_h, "psi4", workdir, command, warnings))

    rescored = rerank_qm_records(rescored)
    write_qm_refined_outputs(rescored, config.output_dir / "qm" / "psi4")
    return rescored


def parse_psi4_energy(logfile: Path) -> float | None:
    if not logfile.exists():
        return None
    text = logfile.read_text(encoding="utf-8", errors="replace")
    patterns = [
        r"Total Energy\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"Final Energy\s*:\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"@\w+.*?Final Energy:\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"energy\s+is\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return float(match.group(1))
    return None


def psi4_workdir(candidate: RankedVariantRecord, config: RunConfig) -> Path:
    return config.output_dir / "qm" / "psi4" / candidate.id


def load_ranked_for_qm(run_dir: Path) -> list[RankedVariantRecord]:
    candidates = [
        run_dir / "censo" / "ranked_variants_refined.json",
        run_dir / "ranking" / "ranked_variants.json",
    ]
    for path in candidates:
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
            return [RankedVariantRecord.model_validate(item) for item in data]
    return []


def write_qm_refined_outputs(records: list[RankedVariantRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_qm_csv(output_dir / "ranked_variants_qm.csv", records)
    (output_dir / "ranked_variants_qm.json").write_text(
        json.dumps([record.model_dump(mode="json") for record in records], indent=2) + "\n",
        encoding="utf-8",
    )


def rerank_qm_records(records: list[RankedVariantRecord]) -> list[RankedVariantRecord]:
    finite = [record.score_kcal_mol for record in records if record.score_kcal_mol is not None]
    minimum = min(finite) if finite else None
    ordered = sorted(
        records,
        key=lambda record: (
            float("inf") if record.score_kcal_mol is None else record.score_kcal_mol,
            record.id,
        ),
    )
    reranked: list[RankedVariantRecord] = []
    for rank, record in enumerate(ordered, start=1):
        data = record.model_dump(mode="python")
        data["rank"] = rank
        if minimum is not None and record.score_kcal_mol is not None:
            data["relative_free_energy_kcal_mol"] = record.score_kcal_mol - minimum
        data["boltzmann_population"] = None
        reranked.append(RankedVariantRecord.model_validate(data))
    return reranked


def _qm_record(
    candidate: RankedVariantRecord,
    energy_hartree: float | None,
    backend: str,
    workdir: Path,
    command: list[str],
    warnings: list[str],
) -> RankedVariantRecord:
    score = None if energy_hartree is None else hartree_to_kcal_mol(energy_hartree)
    data = candidate.model_dump(mode="python")
    data.update(
        {
            "source_software": backend,
            "source_command": " ".join(command),
            "source_python_function": (
                f"dsvr.runners.{backend}_runner.rescore_top_ranked_with_{backend}"
            ),
            "score_kcal_mol": score,
            "warnings": warnings,
            "metadata": candidate.metadata
            | {
                "qm": {
                    "backend": backend,
                    "workdir": str(workdir),
                    "energy_hartree": energy_hartree,
                    "energy_kcal_mol": score,
                    "optimize": False,
                    "single_point": True,
                    "preserves_preliminary_ranking_id": candidate.id,
                }
            },
        }
    )
    return RankedVariantRecord.model_validate(data)


def _psi4_command(
    input_path: Path,
    output_path: Path,
    workdir: Path,
    config: RunConfig,
) -> list[str]:
    executable = shutil.which("psi4")
    if executable is not None:
        return [executable, str(input_path), str(output_path)]
    script = workdir / "run_psi4.py"
    script.write_text(
        "\n".join(
            [
                "import pathlib",
                "import psi4",
                f"psi4.set_output_file({str(output_path)!r}, False)",
                f"psi4.energy({f'{config.refinement.qm_method}/{config.refinement.qm_basis}'!r})",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return [sys.executable, str(script)]


def _psi4_input(xyz: str, candidate: RankedVariantRecord, config: RunConfig) -> str:
    charge = candidate.formal_charge or 0
    multiplicity = 1
    body = _xyz_body(xyz)
    lines = [
        "memory 2 GB",
        "molecule dsvr_candidate {",
        f"{charge} {multiplicity}",
        body,
        "}",
        "",
        "set {",
        f"  basis {config.refinement.qm_basis}",
        "}",
    ]
    if config.refinement.qm_solvent:
        lines.extend(
            [
                "# Solvent requested; exact Psi4 solvent setup is environment/method dependent.",
                f"# qm_solvent = {config.refinement.qm_solvent}",
            ]
        )
    method_basis = f"{config.refinement.qm_method}/{config.refinement.qm_basis}"
    if config.refinement.qm_optimize:
        lines.append(f"optimize('{method_basis}')")
    if config.refinement.qm_single_point:
        lines.append(f"energy('{method_basis}')")
    return "\n".join(lines) + "\n"


def _top_n(
    ranked_records: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    return sorted(
        ranked_records,
        key=lambda record: (
            float("inf")
            if record.relative_free_energy_kcal_mol is None
            else record.relative_free_energy_kcal_mol,
            record.rank,
            record.id,
        ),
    )[: config.refinement.max_candidates_for_refinement]


def _geometry_from_candidate(candidate: RankedVariantRecord) -> str:
    source_workdir = candidate.metadata.get("ranking", {}).get("source_workdir")
    if source_workdir is None:
        source_workdir = candidate.metadata.get("censo", {}).get("workdir")
    if source_workdir is None:
        raise Psi4ExecutionError(f"Candidate {candidate.id} lacks source geometry metadata.")
    for pattern in ("crest_conformer_*.xyz", "*.xyz"):
        matches = sorted(Path(source_workdir).glob(pattern))
        if matches:
            return matches[0].read_text(encoding="utf-8")
    raise Psi4ExecutionError(f"No XYZ geometry found for candidate {candidate.id}.")


def _xyz_body(xyz: str) -> str:
    lines = xyz.splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        return "\n".join(lines[2:])
    return "\n".join(lines)


def _write_qm_csv(path: Path, records: list[RankedVariantRecord]) -> None:
    columns = [
        "rank",
        "id",
        "parent_id",
        "input_molecule_id",
        "molname",
        "score_kcal_mol",
        "relative_free_energy_kcal_mol",
        "boltzmann_population",
        "source_software",
        "warnings",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            row = record.model_dump(mode="json")
            row["warnings"] = " | ".join(record.warnings)
            writer.writerow({column: row.get(column) for column in columns})
