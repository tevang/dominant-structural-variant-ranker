"""Automated anomaly detection for a run directory.

Cross-checks counts and flags that a reviewer would otherwise have to find by
hand: irreconcilable input/processed/ranked counts, upstream provenance gaps,
per-stage timeouts and variant rescues.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Anomaly:
    """A single detected discrepancy."""

    level: str
    category: str
    message: str


def _csv_row_count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return max(0, sum(1 for _ in handle) - 1)
    except OSError:
        return 0


def _distinct_ranked_molecules(path: Path) -> int | None:
    if not path.exists():
        return None
    names: set[str] = set()
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                value = row.get("parent_name") or row.get("molname")
                if value:
                    names.add(value)
    except (OSError, csv.Error):
        return None
    return len(names)


def _stage_summary_metrics(path: Path) -> tuple[int | None, int | None, int]:
    accepted: int | None = None
    rejected: int | None = None
    timeout_sum = 0
    if not path.exists():
        return None, None, 0
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if str(row.get("stage", "")).strip() == "Input validation":
                    accepted = _digit_cell(row.get("accepted_count"))
                    rejected = _digit_cell(row.get("rejected_count"))
                timeout_sum += _digit_cell(row.get("timeout_count")) or 0
    except (OSError, csv.Error):
        return None, None, 0
    return accepted, rejected, timeout_sum


def _digit_cell(value: object) -> int | None:
    text = str(value or "").strip()
    return int(text) if text.isdigit() else None


def _rescue_rules(path: Path) -> set[str]:
    rules: set[str] = set()
    if not path.exists():
        return rules
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rule = row.get("rescue_rule")
                if rule:
                    rules.add(rule)
    except (OSError, csv.Error):
        return rules
    return rules


def detect_anomalies(run_dir: Path | str) -> list[Anomaly]:
    """Return the anomalies found for a run directory (empty when none)."""
    run_dir = Path(run_dir)
    anomalies: list[Anomaly] = []
    invalid_count = _csv_row_count(run_dir / "invalid_inputs.csv")
    inputs_count = _csv_row_count(run_dir / "inputs.csv")
    ranked_molecules = _distinct_ranked_molecules(run_dir / "ranked.csv")
    stage_accepted, stage_rejected, timeout_sum = _stage_summary_metrics(
        run_dir / "stage_summary.csv"
    )
    rescue_rules = _rescue_rules(run_dir / "variant_selection.csv")

    if (
        stage_accepted is not None
        and ranked_molecules is not None
        and stage_accepted != ranked_molecules
    ):
        anomalies.append(
            Anomaly(
                level="error",
                category="count_mismatch",
                message=(
                    f"Input validation accepted {stage_accepted} molecules but "
                    f"{ranked_molecules} distinct molecules are ranked in ranked.csv."
                ),
            )
        )

    if ranked_molecules and ranked_molecules > 0 and inputs_count == 0:
        anomalies.append(
            Anomaly(
                level="warning",
                category="provenance_gap",
                message=(
                    f"inputs.csv is empty yet {ranked_molecules} molecules were ranked; "
                    "the valid-input provenance cannot be reconciled."
                ),
            )
        )

    if invalid_count > 0 and stage_rejected is not None and invalid_count != stage_rejected:
        anomalies.append(
            Anomaly(
                level="warning",
                category="invalid_yet_ranked",
                message=(
                    f"invalid_inputs.csv records {invalid_count} failed inputs but the "
                    f"Input Validation stage reports {stage_rejected} rejected; "
                    "the rejection accounting is inconsistent."
                ),
            )
        )

    if timeout_sum > 0:
        anomalies.append(
            Anomaly(
                level="error",
                category="timeouts",
                message=f"Stages report {timeout_sum} timeout(s) across the run.",
            )
        )

    if rescue_rules:
        joined = ", ".join(sorted(rescue_rules))
        anomalies.append(
            Anomaly(
                level="warning",
                category="rescues",
                message=f"Variants were rescued by fallback rule(s): {joined}.",
            )
        )

    return anomalies
