from __future__ import annotations

import signal
from collections.abc import Iterator
from contextlib import contextmanager


class StepTimeoutError(TimeoutError):
    pass


@contextmanager
def time_limit(seconds: int, label: str) -> Iterator[None]:
    if seconds <= 0 or not hasattr(signal, "SIGALRM"):
        yield
        return

    def _handler(signum, frame):
        raise StepTimeoutError(f"{label} exceeded {seconds} seconds")

    previous = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)
