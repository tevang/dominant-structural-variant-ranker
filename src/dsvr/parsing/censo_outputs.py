from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from dsvr.utils.units import hartree_to_kcal_mol


@dataclass(frozen=True)
class CensoCandidateResult:
    conformer_index: int
    free_energy_kcal_mol: float | None = None
    relative_free_energy_kcal_mol: float | None = None
    population: float | None = None


@dataclass(frozen=True)
class CensoResult:
    candidates: list[CensoCandidateResult]
    output_path: Path | None
    warnings: list[str]
    raw_excerpt: str


def parse_censo_output(path: Path) -> CensoResult:
    if not path.exists():
        return CensoResult(
            candidates=[],
            output_path=None,
            warnings=[f"CENSO output not found: {path}"],
            raw_excerpt="",
        )
    text = path.read_text(encoding="utf-8", errors="replace")
    candidates = _parse_table_like_results(text)
    if not candidates:
        candidates = _parse_single_result(text)
    warnings = [] if candidates else ["No CENSO refined energies/populations parsed."]
    return CensoResult(
        candidates=candidates,
        output_path=path,
        warnings=warnings,
        raw_excerpt=text[:500],
    )


def parse_censo_summary(text: str) -> dict[str, str]:
    return {"raw_excerpt": text[:200]}


def _parse_table_like_results(text: str) -> list[CensoCandidateResult]:
    results: list[CensoCandidateResult] = []
    for line in text.splitlines():
        if not re.search(r"\d", line):
            continue
        lower = line.lower()
        if not any(token in lower for token in ("conf", "conformer", "candidate")):
            continue
        numbers = [float(match) for match in re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", line)]
        integer = re.search(r"(?:conf(?:ormer)?|candidate)\s*#?\s*(\d+)", line, re.IGNORECASE)
        if not numbers:
            continue
        conformer_index = int(integer.group(1)) if integer else len(results) + 1
        free_energy = _normalize_energy(numbers[0], line)
        relative = _normalize_energy(numbers[1], line) if len(numbers) > 1 else None
        population = numbers[2] if len(numbers) > 2 else _population_from_line(line)
        results.append(
            CensoCandidateResult(
                conformer_index=conformer_index,
                free_energy_kcal_mol=free_energy,
                relative_free_energy_kcal_mol=relative,
                population=population,
            )
        )
    return results


def _parse_single_result(text: str) -> list[CensoCandidateResult]:
    free_energy = _first_energy(
        text,
        [
            r"final\s+gibbs\s+free\s+energy\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"gibbs\s+free\s+energy\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"free\s+energy\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
    )
    relative = _first_energy(
        text,
        [
            r"relative\s+(?:free\s+)?energy\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
            r"delta\s*g\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)",
        ],
    )
    population = _population_from_line(text)
    if free_energy is None and relative is None and population is None:
        return []
    return [
        CensoCandidateResult(
            conformer_index=1,
            free_energy_kcal_mol=free_energy,
            relative_free_energy_kcal_mol=relative,
            population=population,
        )
    ]


def _first_energy(text: str, patterns: list[str]) -> float | None:
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return _normalize_energy(float(match.group(1)), text)
    return None


def _normalize_energy(value: float, context: str) -> float:
    lower = context.lower()
    if "hartree" in lower or " eh" in lower:
        return hartree_to_kcal_mol(value)
    if "ev" in lower:
        return value * 23.060548867
    return value


def _population_from_line(line: str) -> float | None:
    match = re.search(
        r"(?:population|pop)\s*[:=]?\s*([-+]?\d+\.\d+(?:[Ee][-+]?\d+)?)(\s*%)?",
        line,
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    value = float(match.group(1))
    return value / 100.0 if match.group(2) else value
