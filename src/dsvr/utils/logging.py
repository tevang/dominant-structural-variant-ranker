from __future__ import annotations

import logging
from pathlib import Path

from rich.logging import RichHandler

LOG_FORMAT = "%(asctime)s %(levelname)s [%(name)s] %(message)s"


def configure_logging(
    *,
    level: str = "INFO",
    log_file: Path | None = None,
    rich: bool = True,
) -> None:
    handlers: list[logging.Handler] = []
    if rich:
        handlers.append(RichHandler(markup=True, rich_tracebacks=True))
    else:
        handlers.append(logging.StreamHandler())
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format=LOG_FORMAT,
        handlers=handlers,
        force=True,
    )


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)
