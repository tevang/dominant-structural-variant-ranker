"""Tests for the DSVR run inspection GUI logic and views."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from dsvr.cli import app
from dsvr.gui.anomalies import Anomaly, detect_anomalies
from dsvr.gui.inventory import iter_artifacts, load_inventory
from dsvr.gui.lineage import parse_variant_id
from dsvr.gui.tables import CsvStream, paged_rows


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def make_run_dir(
    root: Path,
    *,
    ranked_molecules: int = 11,
    invalid_inputs: int = 0,
    stage_rejected: int | None = None,
    inputs_empty: bool = True,
    timeouts: int = 0,
    rescue_rule: str | None = "rescue_original_input_state",
) -> Path:
    if stage_rejected is None:
        stage_rejected = invalid_inputs
    run_dir = root / "run"
    run_dir.mkdir(parents=True, exist_ok=True)

    inventory = [
        ["artifact", "path", "kind", "exists", "size_bytes", "record_count", "target"],
        ["stage_summary.csv", "stage_summary.csv", "table_csv", "True", "100", "9", ""],
        ["ranked.csv", "ranked.csv", "table_csv", "True", "100", str(ranked_molecules), ""],
        ["ranked_variants.sdf", "ranked_variants.sdf", "structure_sdf", "True", "100", "76", ""],
        ["summary.md", "summary.md", "report_markdown", "True", "10", "", ""],
    ]
    _write_csv(run_dir / "run_outputs.csv", inventory[0], inventory[1:])

    _write_csv(
        run_dir / "ranked.csv",
        ["variant_id", "parent_name", "smiles", "relative_energy_kcal_mol",
         "approximate_population", "rank"],
        [
            [f"mol_{i:06d}_p01", f"mol{i}", "CCO", "0.0", "1.0", str(i)]
            for i in range(1, ranked_molecules + 1)
        ],
    )

    _write_csv(
        run_dir / "stage_summary.csv",
        ["stage", "status", "accepted_count", "rejected_count", "timeout_count"],
        [
            [
                "Input validation",
                "completed",
                str(ranked_molecules),
                str(stage_rejected),
                str(max(0, timeouts - 1)),
            ],
            ["Ranking", "completed", str(ranked_molecules), "0", str(min(1, timeouts))],
        ],
    )

    if inputs_empty:
        (run_dir / "inputs.csv").write_text("id\n", encoding="utf-8")
    else:
        _write_csv(
            run_dir / "inputs.csv",
            ["id", "molname", "smiles"],
            [[f"mol_{i:06d}", f"mol{i}", "CCO"] for i in range(1, ranked_molecules + 1)],
        )

    invalid_header = [
        "input_id", "source_format", "line_number", "name", "smiles", "raw_record", "error",
    ]
    if invalid_inputs:
        _write_csv(
            run_dir / "invalid_inputs.csv",
            invalid_header,
            [
                [f"mol_{i:06d}", "smiles", str(i), "bad", "not_a_smiles", "x", "failed"]
                for i in range(1, invalid_inputs + 1)
            ],
        )
    else:
        with (run_dir / "invalid_inputs.csv").open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerow(invalid_header)

    if rescue_rule:
        _write_csv(
            run_dir / "variant_selection.csv",
            ["variant_id", "molname", "rescue_rule"],
            [["mol_000001_p01_a_t01_b_c01_c", "m0", rescue_rule]],
        )
    else:
        _write_csv(
            run_dir / "variant_selection.csv",
            ["variant_id", "molname", "rescue_rule"],
            [["mol_000001_p01_a_t01_b_c01_c", "m0", ""]],
        )

    (run_dir / "warnings.jsonl").write_text("", encoding="utf-8")
    (run_dir / "failures.jsonl").write_text("", encoding="utf-8")
    (run_dir / "done.json").write_text(json.dumps({"status": "done"}), encoding="utf-8")
    return run_dir


# --- inventory ---


def test_inventory_parses_artifacts_and_kinds(tmp_path: Path) -> None:
    run_dir = make_run_dir(tmp_path)
    inventory = load_inventory(run_dir)
    assert len(inventory.artifacts) == 4
    by_kind = inventory.by_kind("table_csv")
    assert {a.name for a in by_kind} == {"stage_summary.csv", "ranked.csv"}
    sdf = inventory.get("ranked_variants.sdf")
    assert sdf is not None
    assert sdf.kind == "structure_sdf"
    assert sdf.exists is True
    assert sdf.record_count == 76
    assert inventory.kinds() == ["report_markdown", "structure_sdf", "table_csv"]


def test_inventory_missing_file(tmp_path: Path) -> None:
    empty = tmp_path / "norecord"
    assert iter_artifacts(empty / "run_outputs.csv") == []
    assert load_inventory(empty).exists("anything") is False


# --- lineage ---


def test_parse_variant_id_representative() -> None:
    lineage = parse_variant_id(
        "mol_000001_p01_1b542e9ed2_t01_f2509deb2d_s01_1c8a7a990d_c01_ee65a959b2_rank0001_abc12345"
    )
    assert lineage.molecule == "mol_000001"
    assert lineage.molecule_index == 1
    assert lineage.protomer == "p01_1b542e9ed2"
    assert lineage.tautomer == "t01_f2509deb2d"
    assert lineage.stereoisomer == "s01_1c8a7a990d"
    assert lineage.conformer == "c01_ee65a959b2"
    assert lineage.rank == 1


def test_parse_variant_id_intermediate() -> None:
    lineage = parse_variant_id("mol_000003_p02_be845b0bd4_t01_8f1627e73c")
    assert lineage.molecule == "mol_000003"
    assert lineage.molecule_index == 3
    assert lineage.protomer == "p02_be845b0bd4"
    assert lineage.tautomer == "t01_8f1627e73c"
    assert lineage.stereoisomer is None
    assert lineage.rank is None


def test_parse_variant_id_malformed() -> None:
    lineage = parse_variant_id("totally_not_a_variant")
    assert lineage.molecule == ""
    assert lineage.molecule_index is None
    assert lineage.rank is None


# --- lazy tables ---


def test_paged_rows_returns_page_and_total(tmp_path: Path) -> None:
    path = tmp_path / "big.csv"
    rows = [[f"v{i}", f"mol{i % 3}"] for i in range(1000)]
    _write_csv(path, ["variant", "mol"], rows)
    header, page, total = paged_rows(path, offset=0, limit=50)
    assert header == ["variant", "mol"]
    assert len(page) == 50
    assert page[0] == ["v0", "mol0"]
    assert total == 1000

    header, page2, _ = paged_rows(path, offset=50, limit=50)
    assert page2[0] == ["v50", "mol2"]
    assert page2 != page


def test_stream_is_lazy(tmp_path: Path) -> None:
    path = tmp_path / "big.csv"
    _write_csv(path, ["id"], [[f"v{i}"] for i in range(10000)])
    stream = CsvStream(path)
    assert stream.header == ["id"]
    # rows() is a generator; requesting a page builds only the requested rows.
    iterator = stream.rows()
    assert iter(iterator) is iterator
    first = next(iterator)
    assert first == ["v0"]


# --- anomalies ---


def test_anomalies_clean_run(tmp_path: Path) -> None:
    run_dir = make_run_dir(
        tmp_path,
        ranked_molecules=11,
        invalid_inputs=0,
        inputs_empty=False,
        rescue_rule=None,
    )
    assert detect_anomalies(run_dir) == []


def test_anomalies_invalid_yet_ranked_and_provenance_gap(tmp_path: Path) -> None:
    # empty root inputs.csv plus invalid rows that the stage never accounted
    run_dir = make_run_dir(
        tmp_path,
        ranked_molecules=11,
        invalid_inputs=12,
        stage_rejected=0,
        inputs_empty=True,
    )
    anomalies = detect_anomalies(run_dir)
    categories = {a.category for a in anomalies}
    assert "invalid_yet_ranked" in categories
    assert "provenance_gap" in categories
    assert all(isinstance(a, Anomaly) for a in anomalies)


def test_anomalies_healthy_run_with_rejections(tmp_path: Path) -> None:
    # rejections present but truthfully accounted: no spurious anomalies
    run_dir = make_run_dir(
        tmp_path,
        ranked_molecules=7,
        invalid_inputs=4,
        inputs_empty=False,
        rescue_rule=None,
    )
    assert detect_anomalies(run_dir) == []


def test_anomalies_timeout_and_rescue(tmp_path: Path) -> None:
    run_dir = make_run_dir(tmp_path, ranked_molecules=3, inputs_empty=False, timeouts=2)
    categories = {a.category for a in detect_anomalies(run_dir)}
    assert "timeouts" in categories
    assert "rescues" in categories


# --- CLI ---


def test_view_requires_run_outputs(tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["view", str(tmp_path)])
    assert result.exit_code != 0
    assert "run_outputs.csv" in result.output


def test_view_missing_gui_extra(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    run_dir = make_run_dir(tmp_path)
    real_import = __import__

    def fake_import(name, *args, **kwargs):
        if name == "streamlit" or name.startswith("streamlit."):
            raise ImportError("no streamlit")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", fake_import)
    runner = CliRunner()
    result = runner.invoke(app, ["view", str(run_dir)])
    assert result.exit_code != 0
    assert "gui" in result.output and "extra" in result.output


# --- depiction ---


def test_smiles_to_svg_valid_and_invalid() -> None:
    from dsvr.gui.ui.depict import smiles_to_svg

    svg = smiles_to_svg("CCO")
    assert svg is not None
    assert "<svg" in svg
    assert smiles_to_svg("::::not_a_molecule::::") is None
    assert smiles_to_svg("") is None


# --- streamlit views (skip when streamlit not installed) ---

pytest.importorskip("streamlit", reason="streamlit is not installed")


def test_all_views_render_without_exception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from streamlit.testing.v1 import AppTest

    run_dir = make_run_dir(tmp_path, ranked_molecules=11, invalid_inputs=12, inputs_empty=True)
    entry = (
        Path(__file__).parent.parent / "src" / "dsvr" / "gui" / "ui" / "streamlit_entry.py"
    )
    monkeypatch.setenv("DSVR_VIEW_RUNDIR", str(run_dir))
    at = AppTest.from_file(str(entry))
    at.run()
    assert not at.exception, [str(e) for e in at.exception]
    options = list(at.sidebar.radio[0].options)
    assert "Overview" in options
    at.sidebar.radio[0].set_value("Overview").run()
    assert not at.exception, [str(e) for e in at.exception]
