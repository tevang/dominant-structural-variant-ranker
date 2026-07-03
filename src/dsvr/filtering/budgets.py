from __future__ import annotations

from dsvr.config import RunConfig

MODE_SCALE = {
    "conservative": 0.5,
    "balanced": 1.0,
    "aggressive": 0.5,
    "exhaustive": float("inf"),
}


def budget(config: RunConfig, name: str) -> int | None:
    if not config.variant_filtering.enabled or config.variant_filtering.mode == "exhaustive":
        return None
    value = getattr(config.variant_filtering, name)
    scale = MODE_SCALE.get(config.variant_filtering.mode, 1.0)
    if scale == float("inf"):
        return None
    return max(config.variant_filtering.min_variants_to_keep, int(value * scale))


def seed_budget(config: RunConfig) -> int | None:
    if not config.variant_filtering.enabled or config.variant_filtering.mode == "exhaustive":
        return None
    return config.variant_filtering.max_seeds_per_variant
