"""Streamlit entry script for ``dsvr view``.

This file is executed by `streamlit run`; it reads the run directory and
browser preference from environment variables set by the CLI launcher and
delegates to the shared render function.
"""

from __future__ import annotations

import os

from dsvr.gui.ui.app import render


def main() -> None:
    run_dir = os.environ.get("DSVR_VIEW_RUNDIR", "")
    render(run_dir)


if __name__ == "__main__":
    main()
