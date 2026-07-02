from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from dsvr.models import VariantRecord


def write_ranked_csv(path: Path, records: list[VariantRecord]) -> None:
    rows = [
        record.model_dump(mode="json") | {"rank": rank}
        for rank, record in enumerate(records, 1)
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
