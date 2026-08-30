#!/usr/bin/env python3
"""Group duplicate structures in an SDF by the duplicate-report criterion.

Two structures are duplicates when they share identical atom types, 3D
coordinates, and bond connectivity — implemented as the V2000 block from the
counts line through ``M END`` — differing only in molname.

Usage: python scripts/check_duplicate_structures.py <final_variants.sdf> [--report]

Prints a summary line ``duplicate_groups: N`` (plus redundant-copy count) and,
with ``--report``, the molnames of every duplicate group. Exit code is 1 when
duplicates are present, 0 otherwise.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path


def _record_blocks(path: Path) -> list[tuple[str, str]]:
    """Return ``(molname, v2000_block)`` per SDF record."""

    text = path.read_text(encoding="utf-8", errors="replace")
    records = [record.strip("\n") for record in text.split("$$$$") if record.strip()]
    result: list[tuple[str, str]] = []
    for record in records:
        lines = record.splitlines()
        molname = lines[0].strip() if lines else ""
        # The counts line is the 4th line of a V2000 molfile (index 3).
        start = None
        for index, line in enumerate(lines[:10]):
            if "V2000" in line:
                start = index
                break
        if start is None:
            continue
        end = start
        for index in range(start, len(lines)):
            end = index
            if lines[index].startswith("M  END"):
                break
        block = "\n".join(lines[start : end + 1])
        result.append((molname, block))
    return result


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    path = Path(argv[0])
    verbose = "--report" in argv[1:]
    records = _record_blocks(path)
    groups: dict[str, list[str]] = defaultdict(list)
    for molname, block in records:
        groups[block].append(molname)
    duplicate_groups = {block: names for block, names in groups.items() if len(names) > 1}
    redundant = sum(len(names) - 1 for names in duplicate_groups.values())
    print(f"total_structures: {len(records)}")
    print(f"unique_structures: {len(groups)}")
    print(f"duplicate_groups: {len(duplicate_groups)}")
    print(f"redundant_copies: {redundant}")
    if verbose and duplicate_groups:
        for index, names in enumerate(duplicate_groups.values(), start=1):
            print(f"\nGROUP {index} ({len(names)} identical copies)")
            for name in names:
                print(f"  {name}")
    return 1 if duplicate_groups else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
