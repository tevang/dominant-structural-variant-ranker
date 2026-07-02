from __future__ import annotations

from pathlib import Path

from dsvr.models import InputMolecule
from dsvr.utils.hashing import sha256_text


def read_smiles(path: Path) -> list[InputMolecule]:
    molecules: list[InputMolecule] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(maxsplit=1)
            smiles = fields[0]
            name = (
                fields[1].strip()
                if len(fields) == 2 and fields[1].strip()
                else f"mol_{line_number}"
            )
            molecules.append(
                InputMolecule(
                    index=len(molecules),
                    name=name,
                    input_format="smiles",
                    smiles=smiles,
                    source_path=path,
                    input_hash=sha256_text(line),
                    properties={"line_number": str(line_number)},
                )
            )
    return molecules
