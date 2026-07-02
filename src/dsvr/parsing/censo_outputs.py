from __future__ import annotations


def parse_censo_summary(text: str) -> dict[str, str]:
    return {"raw_excerpt": text[:200]}

