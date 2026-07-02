from __future__ import annotations

from dsvr.models import VariantRecord


def rank_records(records: list[VariantRecord]) -> list[VariantRecord]:
    return sorted(
        records,
        key=lambda record: (
            float("inf")
            if record.relative_energy_kcal_mol is None
            else record.relative_energy_kcal_mol,
            record.variant_id,
        ),
    )
