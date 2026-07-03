from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

from dsvr.config import RunConfig
from dsvr.models import RankedVariantRecord
from dsvr.runners.psi4_runner import (
    QM_FREE_ENERGY_WARNING,
    _geometry_from_candidate,
    _qm_record,
    _top_n,
    rerank_qm_records,
    write_qm_refined_outputs,
)
from dsvr.runners.subprocess_utils import ExternalToolError, run_command


class PySCFUnavailableError(RuntimeError):
    """Raised when PySCF rescoring is requested but unavailable."""


def pyscf_available() -> bool:
    return importlib.util.find_spec("pyscf") is not None


def rescore_top_ranked_with_pyscf(
    ranked_records: list[RankedVariantRecord],
    config: RunConfig,
) -> list[RankedVariantRecord]:
    if config.refinement.qm_backend != "pyscf" and not config.refinement.pyscf_enabled:
        return []
    if not pyscf_available():
        raise PySCFUnavailableError(
            "PySCF QM rescoring was requested but PySCF is unavailable. Install pyscf "
            "in the active environment, then run `dsvr doctor`."
        )

    rescored: list[RankedVariantRecord] = []
    for candidate in _top_n(ranked_records, config):
        workdir = pyscf_workdir(candidate, config)
        workdir.mkdir(parents=True, exist_ok=True)
        xyz = _geometry_from_candidate(candidate)
        script = workdir / "pyscf_input.py"
        output_path = workdir / "pyscf_output.log"
        script.write_text(_pyscf_script(xyz, candidate, config, output_path), encoding="utf-8")
        command = [sys.executable, str(script)]
        warnings = [*candidate.warnings, QM_FREE_ENERGY_WARNING]
        try:
            completed = run_command(
                command,
                cwd=workdir,
                timeout_s=None
                if config.crest.walltime_minutes is None
                else config.crest.walltime_minutes * 60,
                log_dir=workdir / "logs",
                command_name="pyscf",
                check=True,
                show_progress=config.logging.tail_subprocess_logs,
            )
            if not output_path.exists():
                output_path.write_text(completed.stdout, encoding="utf-8")
        except ExternalToolError as exc:
            warnings.append(f"PySCF failed: {exc}")
        energy_h = parse_pyscf_energy(output_path)
        rescored.append(_qm_record(candidate, energy_h, "pyscf", workdir, command, warnings))

    rescored = rerank_qm_records(rescored)
    write_qm_refined_outputs(rescored, config.output_dir / "qm" / "pyscf")
    return rescored


def parse_pyscf_energy(logfile: Path) -> float | None:
    if not logfile.exists():
        return None
    text = logfile.read_text(encoding="utf-8", errors="replace")
    patterns = [
        r"DSVR_PYSCF_ENERGY_HARTREE\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"converged\s+SCF\s+energy\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        r"total\s+energy\s*=\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return float(match.group(1))
    return None


def pyscf_workdir(candidate: RankedVariantRecord, config: RunConfig) -> Path:
    return config.output_dir / "qm" / "pyscf" / candidate.id


def _pyscf_script(
    xyz: str,
    candidate: RankedVariantRecord,
    config: RunConfig,
    output_path: Path,
) -> str:
    charge = candidate.formal_charge or 0
    spin = 0
    atom_lines = _xyz_body_for_pyscf(xyz)
    method = config.refinement.qm_method.lower()
    if method not in {"b3lyp", "pbe0", "pbe", "blyp"}:
        method = config.refinement.qm_method
    return (
        "from pathlib import Path\n"
        "from pyscf import dft, gto, scf\n"
        f"mol = gto.M(atom={atom_lines!r}, basis={config.refinement.qm_basis!r}, "
        f"charge={charge}, spin={spin})\n"
        f"mf = dft.RKS(mol) if {method!r} else scf.RHF(mol)\n"
        f"mf.xc = {method!r}\n"
        "energy = mf.kernel()\n"
        f"Path({str(output_path)!r}).write_text(f'DSVR_PYSCF_ENERGY_HARTREE = {{energy}}\\n')\n"
    )


def _xyz_body_for_pyscf(xyz: str) -> str:
    lines = xyz.splitlines()
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
    return "; ".join(line.strip() for line in lines if line.strip())
