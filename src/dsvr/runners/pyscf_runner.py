from __future__ import annotations


def pyscf_available() -> bool:
    try:
        import pyscf  # noqa: F401
    except ImportError:
        return False
    return True

