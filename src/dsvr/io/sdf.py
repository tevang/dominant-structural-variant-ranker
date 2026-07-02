from __future__ import annotations

from pathlib import Path

from dsvr.models import InputMolecule
from dsvr.utils.hashing import sha256_text


def read_sdf(path: Path) -> list[InputMolecule]:
    text = path.read_text(encoding="utf-8")
    records = [record.strip() for record in text.split("$$$$") if record.strip()]
    molecules: list[InputMolecule] = []
    for index, record in enumerate(records):
        lines = record.splitlines()
        name = lines[0].strip() if lines and lines[0].strip() else f"mol_{index + 1}"
        molecules.append(
            InputMolecule(
                index=index,
                name=name,
                input_format="sdf",
                smiles=None,
                source_path=path,
                input_hash=sha256_text(record),
                properties={"sdf_record": record},
            )
        )
    return molecules

