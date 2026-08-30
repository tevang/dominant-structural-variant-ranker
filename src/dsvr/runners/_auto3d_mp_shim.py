"""Multiprocessing shims shared by the generated Auto3D wrapper scripts.

Some sandboxed environments block ``socket.setsockopt`` used by
``multiprocessing.managers.SyncManager``. The wrapper scripts therefore
replace ``multiprocessing.Manager`` with a fake that returns plain queues.

Auto3D communicates with worker processes started in the spawn context. A
queue created in the default fork context cannot be shared with those
workers (``RuntimeError: A SemLock created in a fork context is being
shared with a process in a spawn context``), so the fake manager builds
its queues from the spawn context.
"""

from __future__ import annotations

import multiprocessing as mp
import multiprocessing.context as mp_context


def install_fake_manager() -> None:
    """Replace ``multiprocessing.Manager`` with a spawn-context queue factory."""

    spawn_context = mp.get_context("spawn")

    class _FakeManager:
        def Queue(self, maxsize: int = 0):
            return spawn_context.Queue(maxsize)

    mp.Manager = _FakeManager

    def _fake_context_manager(self):
        return _FakeManager()

    mp_context.BaseContext.Manager = _fake_context_manager
