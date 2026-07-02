from __future__ import annotations


def parse_auto3d_summary(text: str) -> dict[str, str]:
    return {"raw_excerpt": text[:200]}

