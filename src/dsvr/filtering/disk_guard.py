from __future__ import annotations

from pathlib import Path

from dsvr.config import RunConfig


class DiskLimitError(RuntimeError):
    pass


def directory_size_gb(path: Path) -> float:
    if not path.exists():
        return 0.0
    total = sum(item.stat().st_size for item in path.rglob("*") if item.is_file())
    return total / (1024**3)


def enforce_run_disk_limit(path: Path, config: RunConfig) -> None:
    size_gb = directory_size_gb(path)
    if size_gb <= config.disk.max_run_dir_gb:
        return
    message = (
        f"Run directory {path} is {size_gb:.2f} GiB, "
        f"above {config.disk.max_run_dir_gb:.2f} GiB"
    )
    if config.disk.fail_on_disk_limit:
        raise DiskLimitError(message)
