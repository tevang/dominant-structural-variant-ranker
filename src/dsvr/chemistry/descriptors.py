from __future__ import annotations


def descriptor_placeholder(smiles: str) -> dict[str, int]:
    return {"smiles_length": len(smiles)}

