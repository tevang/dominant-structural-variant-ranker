"""Parse the lineage encoded in DSVR variant IDs.

Variant and intermediate IDs encode their provenance in an underscore-separated
token structure such as ``mol_000001_p01_<hash>_t01_<hash>_s01_<hash>_c01_<hash>``
for final variants and ``mol_000001_p01_<hash>_t01_<hash>`` for intermediate
states. Ranked variants append ``_rank0001_<hash>``. Short feature tokens
(``pNN``, ``tNN``, ``sNN``, ``cNN``) are paired with their 8-hex hash token.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class VariantLineage:
    """Decomposed lineage of a variant or intermediate id."""

    molecule: str
    molecule_index: int | None
    protomer: str | None
    tautomer: str | None
    stereoisomer: str | None
    conformer: str | None
    rank: int | None

    def components(self) -> list[tuple[str, str]]:
        parts: list[tuple[str, str]] = [("molecule", self.molecule)]
        for label, value in (
            ("protomer", self.protomer),
            ("tautomer", self.tautomer),
            ("stereoisomer", self.stereoisomer),
            ("conformer", self.conformer),
        ):
            if value:
                parts.append((label, value))
        if self.rank is not None:
            parts.append(("rank", str(self.rank)))
        return parts


_MOL_INDEX_RE = re.compile(r"^\d+$")
_FEATURE_RE = re.compile(r"^(p|t|s|c)(\d+)$")
_RANK_RE = re.compile(r"^rank(\d+)$")
_HASH_RE = re.compile(r"^[0-9a-f]{8,}$")

_KIND_TO_ATTR = {"p": "protomer", "t": "tautomer", "s": "stereoisomer", "c": "conformer"}


def _tokens(variant_id: str) -> list[str]:
    return variant_id.split("_")


def parse_variant_id(variant_id: str) -> VariantLineage:
    """Parse a variant id into its lineage components.

    Unknown or malformed structure is tolerated: the molecule token is always
    captured when present and unrecognised tokens are ignored.
    """
    attrs: dict[str, str | None] = {
        "protomer": None,
        "tautomer": None,
        "stereoisomer": None,
        "conformer": None,
    }
    molecule = ""
    molecule_index: int | None = None
    rank: int | None = None

    tokens = _tokens(variant_id)
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "mol" and index + 1 < len(tokens) and _MOL_INDEX_RE.match(tokens[index + 1]):
            molecule = f"mol_{tokens[index + 1]}"
            molecule_index = int(tokens[index + 1])
            index += 2
            continue
        feature = _FEATURE_RE.match(token)
        if feature:
            kind = feature.group(1)
            label = token
            if index + 1 < len(tokens) and _HASH_RE.match(tokens[index + 1]):
                label = f"{token}_{tokens[index + 1]}"
                index += 2
            else:
                index += 1
            attrs[_KIND_TO_ATTR[kind]] = label
            continue
        rank_match = _RANK_RE.match(token)
        if rank is None and rank_match:
            rank = int(rank_match.group(1))
        index += 1

    return VariantLineage(
        molecule=molecule,
        molecule_index=molecule_index,
        protomer=attrs["protomer"],
        tautomer=attrs["tautomer"],
        stereoisomer=attrs["stereoisomer"],
        conformer=attrs["conformer"],
        rank=rank,
    )
