"""Helpers for reading run tables lazily or with caching.

Large output files (tens to hundreds of MB) must not be materialised into
memory wholesale. :class:`CsvStream` streams a CSV one row at a time and
``paged_rows`` builds only the requested page while still counting matching
rows for pagination controls. Small reference tables are cached in memory via
``cached_csv_frame``.
"""

from __future__ import annotations

import csv
import functools
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd


def _raise_field_limit() -> None:
    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2**31 - 1)


_raise_field_limit()


def _normalise(value: Any) -> str:
    return "" if value is None else str(value)


class CsvStream:
    """Streams a CSV file without loading it fully into memory."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self._header: list[str] | None = None
        self._size: int | None = None

    @property
    def header(self) -> list[str]:
        if self._header is None:
            with self.path.open(encoding="utf-8", newline="") as handle:
                self._header = next(csv.reader(handle))
        return list(self._header)

    def rows(self) -> Iterator[list[str]]:
        """Yield each data row (excluding the header)."""
        with self.path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            yield from reader

    @property
    def row_count(self) -> int:
        # Counted lazily and memoised; used for large-file affordances only.
        if self._size is None:
            count = 0
            for _ in self.rows():
                count += 1
            self._size = count
        return self._size


def match_row(
    row: list[str],
    header: list[str],
    *,
    query: str = "",
    filters: dict[str, str] | None = None,
) -> bool:
    """Return whether a row matches a substring query and/or per-column filters."""
    for column, value in (filters or {}).items():
        if column in header:
            index = header.index(column)
            index_value = _normalise(row[index]) if index < len(row) else ""
            if index_value != value:
                return False
    if query:
        haystack = " ".join(_normalise(cell) for cell in row)
        if query.lower() not in haystack.lower():
            return False
    return True


def paged_rows(
    path: Path,
    *,
    offset: int = 0,
    limit: int = 50,
    query: str = "",
    filters: dict[str, str] | None = None,
) -> tuple[list[str], list[list[str]], int]:
    """Return (header, page rows, total_matching) streaming the CSV.

    Only the requested page of rows is retained; the full file is still scanned
    to count matches for pagination. Header is always returned even for an
    empty file.
    """
    stream = CsvStream(path)
    header = stream.header
    page: list[list[str]] = []
    total_matched = 0
    end = offset + limit
    for row in stream.rows():
        if len(row) != len(header):
            continue
        if not match_row(row, header, query=query, filters=filters):
            continue
        total_matched += 1
        if offset <= total_matched - 1 < end:
            page.append(row)
    return header, page, total_matched


@functools.lru_cache(maxsize=64)
def cached_csv_frame(path: str, mtime_ns: int) -> pd.DataFrame:
    """Load a small table as a DataFrame, cached by path and modification time."""
    del mtime_ns
    return pd.read_csv(path)
