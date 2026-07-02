from __future__ import annotations


def auto3d_available() -> bool:
    try:
        import Auto3D  # noqa: F401
    except ImportError:
        return False
    return True

