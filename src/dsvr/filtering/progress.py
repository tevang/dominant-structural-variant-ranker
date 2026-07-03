from __future__ import annotations

from collections import Counter

from dsvr.filtering.selection import FilteringDecision


def decision_counts(decisions: list[FilteringDecision]) -> dict[str, int]:
    counter = Counter("selected" if decision.selected else "filtered" for decision in decisions)
    return dict(counter)
