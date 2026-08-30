"""Streamlit application entry for the run inspection GUI.

``run_inspection_app`` is invoked by the ``dsvr view`` CLI command. It launches
Streamlit against the committed ``streamlit_entry.py`` script, passing the run
directory and browser preference through environment variables. The actual
Streamlit app body lives in ``render``, which is also called directly by the
entry script.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

_ENV_RUN_DIR = "DSVR_VIEW_RUNDIR"
_ENV_OPEN_BROWSER = "DSVR_VIEW_OPEN_BROWSER"

_VIEWS = [
    ("Overview", "render_overview"),
    ("Molecules", "render_molecules"),
    ("Enumerations", "render_enumerations"),
    ("Ranking", "render_ranking"),
    ("Decisions", "render_decisions"),
    ("Errors & Warnings", "render_errors"),
    ("Logs", "render_logs"),
    ("Artifacts", "render_artifacts"),
]


def _entry_path() -> Path:
    return Path(__file__).with_name("streamlit_entry.py")


def run_inspection_app(run_dir: str, *, open_browser: bool = True) -> None:
    """Launch the Streamlit GUI for a run directory (blocks while running)."""
    from streamlit.web import cli as stcli  # noqa: F401

    env = dict(os.environ)
    env[_ENV_RUN_DIR] = str(run_dir)
    env[_ENV_OPEN_BROWSER] = "1" if open_browser else "0"
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(_entry_path()),
        "--browser.gatherUsageStats",
        "false",
    ]
    if not open_browser:
        command.extend(["--server.headless", "true"])
    subprocess.run(command, env=env, check=False)


def render(run_dir: str) -> None:
    """Render the Streamlit app. Executed as the Streamlit main script body."""
    import streamlit as st

    from dsvr.gui.ui import views  # noqa: F401

    st.set_page_config(page_title="DSVR Run Inspector", layout="wide")
    st.title("DSVR Run Inspector")
    st.caption(f"Run directory: `{run_dir}`")

    path = Path(run_dir)
    with st.sidebar:
        st.markdown("### Views")
        selection = st.radio("Navigate", [label for label, _ in _VIEWS])

    handler_name = dict((label, handler) for label, handler in _VIEWS)[selection]
    handler = getattr(views, handler_name)
    handler(path)
