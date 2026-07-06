from __future__ import annotations

from collections.abc import Mapping
from typing import Any

DISALLOWED_WITHOUT_EXPLICIT_USER_FLAG: tuple[str, ...] = (
    "patch_code",
    "change_scientific_threshold_defaults",
    "delete_final_outputs",
    "rerun_large_validation_jobs",
    "change_chemistry_assumptions",
)

CONFIG_TWEAK_LIMITS: dict[str, tuple[int | float, int | float]] = {
    "enumeration.max_tautomers_per_protomer": (1, 64),
    "enumeration.max_stereoisomers_per_tautomer": (1, 128),
    "tautomer_filtering.max_rdkit_tautomers_before_auto3d": (1, 128),
    "stereoisomer_filtering.max_stereoisomers_per_tautomer": (1, 64),
    "seeding.auto3d_cpu_workers": (1, 64),
    "final_3d.max_confs": (1, 50),
}


def sanitize_config_tweak(tweak: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if not tweak:
        return None
    key = tweak.get("key")
    value = tweak.get("value")
    if not isinstance(key, str) or key not in CONFIG_TWEAK_LIMITS:
        return None
    low, high = CONFIG_TWEAK_LIMITS[key]
    if not isinstance(value, int | float):
        return None
    clamped = max(low, min(high, value))
    if isinstance(low, int) and isinstance(high, int):
        clamped = int(clamped)
    return {"key": key, "value": clamped}
