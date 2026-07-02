from __future__ import annotations

from dsvr.utils.hashing import sha256_text


def stable_molecule_id(text: str) -> str:
    return sha256_text(text)[:16]

