from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from dsvr.utils.units import hartree_to_kcal_mol


@dataclass(frozen=True)
class ParsedCrestConformer:
    index: int
    xyz: str
    energy_kcal_mol: float | None
    relative_energy_kcal_mol: float | None
    comment: str


@dataclass(frozen=True)
class ParsedCrestOutput:
    conformers: list[ParsedCrestConformer]
    energy_source: Path | None
    warnings: list[str]


def parse_crest_outputs(workdir: Path) -> ParsedCrestOutput:
    warnings: list[str] = []
    xyz_path = _first_existing(
        [
            workdir / "crest_conformers.xyz",
            workdir / "crest_ensemble.xyz",
            workdir / "crest_best.xyz",
        ]
    )
    if xyz_path is None:
        return ParsedCrestOutput(
            conformers=[],
            energy_source=None,
            warnings=["No CREST conformer XYZ file found."],
        )

    energy_source = _first_existing(
        [
            workdir / "crest.energies",
            workdir / "crest_energies.log",
            workdir / "energies.log",
        ]
    )
    energies = parse_crest_energy_file(energy_source) if energy_source is not None else []
    if energy_source is None:
        warnings.append("No CREST energy file found; conformers will not have parsed energies.")

    conformer_blocks = parse_multixyz(xyz_path)
    conformers: list[ParsedCrestConformer] = []
    for index, block in enumerate(conformer_blocks, start=1):
        energy = energies[index - 1] if index - 1 < len(energies) else None
        comment_energy = _energy_from_comment(block.comment)
        if energy is None:
            energy = comment_energy
        conformers.append(
            ParsedCrestConformer(
                index=index,
                xyz=block.xyz,
                energy_kcal_mol=energy,
                relative_energy_kcal_mol=None,
                comment=block.comment,
            )
        )

    absolute = [item.energy_kcal_mol for item in conformers if item.energy_kcal_mol is not None]
    minimum = min(absolute) if absolute else None
    if minimum is not None:
        conformers = [
            ParsedCrestConformer(
                index=item.index,
                xyz=item.xyz,
                energy_kcal_mol=item.energy_kcal_mol,
                relative_energy_kcal_mol=(
                    None if item.energy_kcal_mol is None else item.energy_kcal_mol - minimum
                ),
                comment=item.comment,
            )
            for item in conformers
        ]
    return ParsedCrestOutput(conformers=conformers, energy_source=energy_source, warnings=warnings)


@dataclass(frozen=True)
class _XyzBlock:
    atom_count: int
    comment: str
    xyz: str


def parse_multixyz(path: Path) -> list[_XyzBlock]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    blocks: list[_XyzBlock] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError:
            break
        end = cursor + atom_count + 2
        if end > len(lines):
            break
        block_lines = lines[cursor:end]
        blocks.append(
            _XyzBlock(
                atom_count=atom_count,
                comment=block_lines[1] if len(block_lines) > 1 else "",
                xyz="\n".join(block_lines) + "\n",
            )
        )
        cursor = end
    return blocks


def parse_crest_energy_file(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()
    if "kcal" in lower:
        unit = "kcal/mol"
    elif "ev" in lower:
        unit = "ev"
    else:
        unit = "hartree"

    energies: list[float] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        numbers = [
            float(match)
            for match in re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", stripped)
        ]
        if not numbers:
            continue
        energy = numbers[-1]
        if unit == "hartree":
            energy = hartree_to_kcal_mol(energy)
        elif unit == "ev":
            energy *= 23.060548867
        energies.append(energy)
    return energies


def parse_crest_summary(text: str) -> dict[str, str]:
    return {"raw_excerpt": text[:200]}


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def _energy_from_comment(comment: str) -> float | None:
    lower = comment.lower()
    matches = [float(match) for match in re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", comment)]
    if not matches:
        return None
    value = matches[0]
    if "hartree" in lower or "eh" in lower:
        return hartree_to_kcal_mol(value)
    if "ev" in lower:
        return value * 23.060548867
    if "kcal" in lower:
        return value
    return None
