from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from dsvr.utils.units import hartree_to_kcal_mol

EV_TO_KCAL_MOL = 23.060548867


@dataclass(frozen=True)
class XtbEnergy:
    electronic_energy_hartree: float | None = None
    electronic_energy_kcal_mol: float | None = None


@dataclass(frozen=True)
class XtbThermo:
    electronic_energy_hartree: float | None = None
    electronic_energy_kcal_mol: float | None = None
    gibbs_free_energy_hartree: float | None = None
    gibbs_free_energy_kcal_mol: float | None = None
    enthalpy_hartree: float | None = None
    enthalpy_kcal_mol: float | None = None
    entropy_cal_mol_k: float | None = None
    raw_values: dict[str, float] | None = None


def parse_xtb_energy(logfile: Path | str) -> XtbEnergy:
    text = _read_text(logfile)
    energy = _first_float_after_patterns(
        text,
        [
            r"TOTAL\s+ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"total\s+energy\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"\|\s*TOTAL ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
    )
    return XtbEnergy(
        electronic_energy_hartree=energy,
        electronic_energy_kcal_mol=None if energy is None else hartree_to_kcal_mol(energy),
    )


def parse_xtb_thermo(logfile: Path | str) -> XtbThermo:
    text = _read_text(logfile)
    energy = parse_xtb_energy(text)
    gibbs_h = _first_float_after_patterns(
        text,
        [
            r"G(?:\(T\))?\s*/?\s*Eh\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"TOTAL\s+FREE\s+ENERGY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"Gibbs\s+free\s+energy\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"free\s+energy\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
        ],
    )
    enthalpy_h = _first_float_after_patterns(
        text,
        [
            r"H(?:\(T\))?\s*/?\s*Eh\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"TOTAL\s+ENTHALPY\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"enthalpy\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+Eh",
        ],
    )
    entropy = _first_float_after_patterns(
        text,
        [
            r"entropy\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)\s+cal",
            r"S(?:\(T\))?\s*/?\s*cal/mol/K\s+([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
    )
    raw = {
        key: value
        for key, value in {
            "electronic_energy_hartree": energy.electronic_energy_hartree,
            "gibbs_free_energy_hartree": gibbs_h,
            "enthalpy_hartree": enthalpy_h,
            "entropy_cal_mol_k": entropy,
        }.items()
        if value is not None
    }
    return XtbThermo(
        electronic_energy_hartree=energy.electronic_energy_hartree,
        electronic_energy_kcal_mol=energy.electronic_energy_kcal_mol,
        gibbs_free_energy_hartree=gibbs_h,
        gibbs_free_energy_kcal_mol=None if gibbs_h is None else hartree_to_kcal_mol(gibbs_h),
        enthalpy_hartree=enthalpy_h,
        enthalpy_kcal_mol=None if enthalpy_h is None else hartree_to_kcal_mol(enthalpy_h),
        entropy_cal_mol_k=entropy,
        raw_values=raw,
    )


def _read_text(logfile: Path | str) -> str:
    if isinstance(logfile, Path):
        return logfile.read_text(encoding="utf-8", errors="replace")
    if "\n" in logfile or "TOTAL" in logfile.upper() or "GIBBS" in logfile.upper():
        return logfile
    path = Path(logfile)
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else logfile


def _first_float_after_patterns(text: str, patterns: list[str]) -> float | None:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return float(match.group(1))
    return None
