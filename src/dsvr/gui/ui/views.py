"""View render functions for the run inspection GUI."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

from dsvr.gui.anomalies import detect_anomalies
from dsvr.gui.inventory import load_inventory
from dsvr.gui.tables import cached_csv_frame, paged_rows
from dsvr.gui.ui.depict import smiles_to_svg

_JSONL_COLUMNS = ["timestamp", "event", "stage", "stage_name", "status", "message"]


def _line_count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return sum(1 for _ in handle)
    except OSError:
        return 0


def _read_done_status(run_dir: Path) -> dict:
    path = run_dir / "done.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _read_jsonl(path: Path) -> list[dict]:
    if not path.exists():
        return []
    records: list[dict] = []
    try:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    records.append({"message": line, "level": "unknown"})
    except OSError:
        return []
    return records


def render_paged_table(path: Path, *, key: str) -> None:
    """Render a table with streaming pagination for arbitrarily large files."""
    header, _, total = paged_rows(path, offset=0, limit=1)
    if not header:
        st.info("No rows to display.")
        return
    page_size = st.select_slider(
        f"Page size ({key})",
        options=[25, 50, 100, 200],
        value=50,
    )
    max_page = max(0, (total - 1) // page_size)
    page = st.number_input(
        f"Page (of {max_page + 1})",
        min_value=0,
        max_value=max_page,
        value=0,
        step=1,
        key=f"page_{key}",
    )
    query = st.text_input("Search", key=f"search_{key}").strip()
    _, page_rows, _ = paged_rows(
        path,
        offset=int(page * page_size),
        limit=int(page_size),
        query=query,
        filters={},
    )
    frame = pd.DataFrame(page_rows, columns=header)
    st.caption(f"{total} matching rows")
    st.dataframe(frame, use_container_width=True)


def render_overview(run_dir: Path) -> None:
    st.subheader("Run Overview")
    done = _read_done_status(run_dir)
    if done:
        st.write(f"Status: **{done.get('status', 'unknown')}**")
        if done.get("completed_at"):
            st.write(f"Completed: {done['completed_at']}")

    anomalies = detect_anomalies(run_dir)
    if anomalies:
        st.error(f"{len(anomalies)} potential issue(s) detected")
        for anomaly in anomalies:
            if anomaly.level == "error":
                st.error(f"[{anomaly.category}] {anomaly.message}")
            else:
                st.warning(f"[{anomaly.category}] {anomaly.message}")
    else:
        st.success("No anomalies detected.")

    warning_count = _line_count(run_dir / "warnings.jsonl")
    failure_count = _line_count(run_dir / "failures.jsonl")
    c1, c2 = st.columns(2)
    c1.metric("Warnings", warning_count)
    c2.metric("Failures", failure_count)

    stage_path = run_dir / "stage_summary.csv"
    if stage_path.exists():
        st.subheader("Stage Summary")
        frame = cached_csv_frame(str(stage_path), int(stage_path.stat().st_mtime_ns))
        show = [
            c
            for c in [
                "stage",
                "status",
                "generated_count",
                "accepted_count",
                "selected_count",
                "rejected_count",
                "skipped_count",
                "timeout_count",
                "elapsed_seconds",
                "run_dir_size_mb",
            ]
            if c in frame.columns
        ]
        st.dataframe(frame[show], use_container_width=True)


def render_molecules(run_dir: Path) -> None:
    st.subheader("Molecules")

    invalid_path = run_dir / "invalid_inputs.csv"
    if invalid_path.exists() and _line_count(invalid_path) > 1:
        st.error("Inputs that failed validation:")
        header, rows, _ = paged_rows(invalid_path, offset=0, limit=100)
        if header:
            st.dataframe(pd.DataFrame(rows, columns=header), use_container_width=True)

    inputs_path = run_dir / "inputs.csv"
    if inputs_path.exists() and _line_count(inputs_path) > 1:
        st.subheader("Valid Inputs")
        header, rows, total = paged_rows(inputs_path, offset=0, limit=100)
        st.dataframe(pd.DataFrame(rows, columns=header), use_container_width=True)
    else:
        st.info("inputs.csv is empty or missing.")

    ranked_path = run_dir / "ranked.csv"
    if not ranked_path.exists():
        st.info("No ranked.csv found.")
        return

    rank_frame = cached_csv_frame(str(ranked_path), int(ranked_path.stat().st_mtime_ns))
    name_col = "parent_name" if "parent_name" in rank_frame.columns else "molname"
    if name_col not in rank_frame.columns:
        st.info("ranked.csv has no molecule-name column.")
        return
    molecules = sorted(rank_frame[name_col].dropna().unique().tolist())
    selected = st.selectbox("Molecule", ["-- all --"] + molecules)
    subset = rank_frame if selected == "-- all --" else rank_frame[rank_frame[name_col] == selected]
    view = st.selectbox("Show", ["Table", "Depictions"])
    if view == "Depictions":
        col = "smiles" if "smiles" in subset.columns else None
        if col:
            cols = st.columns(3)
            for i, (_, row) in enumerate(subset.head(30).iterrows()):
                svg = smiles_to_svg(str(row[col]))
                with cols[i % 3]:
                    st.markdown(svg if svg else "_(unparseable SMILES)_", unsafe_allow_html=True)
                    st.caption(str(row.get("variant_id", "")))
        else:
            st.info("No SMILES column available for depictions.")
    else:
        st.dataframe(subset.head(500), use_container_width=True)


def render_enumerations(run_dir: Path) -> None:
    st.subheader("Enumerations")

    for name, title in [
        ("enumeration_counts.csv", "Enumeration Counts"),
        ("variant_counts.csv", "Variant Counts"),
        ("timing_by_stage.csv", "Timing by Stage"),
        ("disk_usage_by_stage.csv", "Disk Usage by Stage"),
    ]:
        path = run_dir / name
        if not path.exists():
            continue
        st.subheader(title)
        header, rows, total = paged_rows(path, offset=0, limit=200)
        st.dataframe(pd.DataFrame(rows, columns=header), use_container_width=True)

    for name, title in [
        ("protomers_rejected.csv", "Rejected Protomers"),
        ("tautomers_rejected.csv", "Rejected Tautomers"),
        ("stereoisomers_rejected.csv", "Rejected Stereoisomers"),
    ]:
        path = run_dir / name
        if not path.exists():
            continue
        st.subheader(title)
        header, rows, total = paged_rows(path, offset=0, limit=100)
        st.dataframe(pd.DataFrame(rows, columns=header), use_container_width=True)


def render_ranking(run_dir: Path) -> None:
    st.subheader("Ranking")
    ranked_path = run_dir / "ranked.csv"
    if not ranked_path.exists():
        st.info("No ranked.csv found.")
        return
    rank_frame = cached_csv_frame(str(ranked_path), int(ranked_path.stat().st_mtime_ns))
    name_col = "parent_name" if "parent_name" in rank_frame.columns else "molname"
    if name_col in rank_frame.columns:
        molecules = ["-- all --"] + sorted(rank_frame[name_col].dropna().unique().tolist())
        selected = st.selectbox("Molecule", molecules)
        if selected != "-- all --":
            rank_frame = rank_frame[rank_frame[name_col] == selected]
    st.dataframe(rank_frame, use_container_width=True)

    energies_path = run_dir / "final_variant_energies.csv"
    if energies_path.exists():
        st.subheader("Final Variant Energies")
        header, rows, _ = paged_rows(energies_path, offset=0, limit=200)
        st.dataframe(pd.DataFrame(rows, columns=header), use_container_width=True)


def render_decisions(run_dir: Path) -> None:
    st.subheader("Decisions")
    selection_path = run_dir / "variant_selection.csv"
    if selection_path.exists():
        selection = cached_csv_frame(str(selection_path), int(selection_path.stat().st_mtime_ns))
        cols = [
            c
            for c in [
                "molname",
                "svp_score_total",
                "accepted_for_3d",
                "rejection_stage",
                "rejection_reason",
                "rescue_rule",
            ]
            if c in selection.columns
        ]
        if "rejection_reason" in selection.columns:
            reasons = ["-- all --"] + sorted(
                selection["rejection_reason"].dropna().unique().tolist()
            )
            reason = st.selectbox("Rejection reason", reasons)
        else:
            reason = "-- all --"
        filtered = selection
        if reason != "-- all --":
            filtered = selection[selection["rejection_reason"] == reason]
        st.dataframe(filtered[cols] if cols else filtered, use_container_width=True)

    decisions_path = run_dir / "variant_decisions.csv"
    if decisions_path.exists():
        st.subheader("Variant Decisions")
        render_paged_table(decisions_path, key="variant_decisions")


def _render_jsonl(title: str, path: Path) -> None:
    records = _read_jsonl(path)
    st.subheader(title)
    if not records:
        st.info("No records.")
        return
    for record in records:
        level = record.get("level", "info")
        message = record.get("message") or record.get("status") or json.dumps(record)
        if len(message) > 500:
            with st.expander(str(message)[:120] + " ..."):
                st.code(message, language=None)
        elif level == "error":
            st.error(str(message))
        elif level == "warning":
            st.warning(str(message))
        else:
            st.write(str(message))


def render_errors(run_dir: Path) -> None:
    st.subheader("Errors & Warnings")
    done = _read_done_status(run_dir)
    if done:
        st.json(done)
    _render_jsonl("Warnings", run_dir / "warnings.jsonl")
    _render_jsonl("Failures", run_dir / "failures.jsonl")


def render_logs(run_dir: Path) -> None:
    st.subheader("Logs")
    logs_dir = run_dir / "logs"
    if not logs_dir.exists():
        st.info("No logs directory found.")
        return
    files = sorted(p.name for p in logs_dir.iterdir() if p.is_file())
    selected = st.selectbox("Log file", files)
    if not selected:
        return
    path = logs_dir / selected
    text = path.read_text(encoding="utf-8", errors="replace")
    st.code(text[-20000:], language=None)
    st.caption(f"Showing last {min(20000, len(text))} of {len(text)} characters.")


def render_artifacts(run_dir: Path) -> None:
    st.subheader("Artifacts")
    inventory = load_inventory(run_dir)
    st.caption("SDF/structural artifacts are shown as metadata only (not rendered).")
    for kind in inventory.kinds():
        with st.expander(f"{kind} ({len(inventory.by_kind(kind))})"):
            artifacts = inventory.by_kind(kind)
            frame = pd.DataFrame(
                [
                    {
                        "artifact": a.name,
                        "exists": a.exists,
                        "size_bytes": a.size_bytes,
                        "record_count": a.record_count,
                        "target": a.target,
                    }
                    for a in artifacts
                ]
            )
            st.dataframe(frame, use_container_width=True)

    st.subheader("Open Artifact")
    names = [
        a.name
        for a in inventory.artifacts
        if a.kind in {"table_csv", "metadata_json", "report_markdown"}
    ]
    selected = st.selectbox("Artifact", names)
    if not selected:
        return
    artifact = inventory.get(selected)
    if artifact is None:
        return
    path = run_dir / artifact.path
    if not path.exists():
        st.warning("Artifact is recorded but missing on disk.")
        return
    if artifact.kind == "table_csv":
        render_paged_table(path, key=f"artifact_{artifact.name}")
    else:
        text = path.read_text(encoding="utf-8", errors="replace")
        st.code(text[:20000], language=None)
