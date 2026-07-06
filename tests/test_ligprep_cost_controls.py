from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.agent.local_qwen import run_local_diagnostic_agent
from dsvr.chemistry import stereochemistry
from dsvr.chemistry import tautomer_auto3d_filter as tautomer_filter
from dsvr.chemistry.final3d import generate_final_3d_variants
from dsvr.chemistry.protonation import _records_from_candidates
from dsvr.chemistry.stereochemistry import enumerate_stereoisomers
from dsvr.chemistry.tautomer_auto3d_filter import (
    RdkitTautomerFilteringTimeout,
    filter_tautomers_with_auto3d,
)
from dsvr.config import AgentConfig, RunConfig, load_config
from dsvr.models import (
    MoleculeInput,
    ProtomerRecord,
    StereoRecord,
    TautomerRecord,
    make_input_id,
    make_protomer_id,
    make_stereo_id,
    make_tautomer_id,
)
from dsvr.runners.auto3d_runner import Auto3DExecutionError
from dsvr.workflow import engine as engine_module
from dsvr.workflow.engine import run_workflow


def test_ligprep_default_bounds_runaway_enumeration_with_mocked_tools(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"
    config_data = load_config(Path("configs/ligprep_like_default.yaml")).model_dump(mode="python")
    config_data.update({"input_path": input_path, "output_dir": outdir, "overwrite": True, "resume": False})
    config = RunConfig.model_validate(config_data)
    assert config.optional_validation.crest_xtb_enabled is False
    assert config.agent.enabled is False

    def fake_protomers(molecule: MoleculeInput, run_config: RunConfig) -> list[ProtomerRecord]:
        raw = [Chem.MolFromSmiles("C" * count) for count in range(1, 21)]
        assert len(raw) == 20
        output_dir = run_config.output_dir / "enumeration" / "protomers"
        output_dir.mkdir(parents=True, exist_ok=True)
        return _records_from_candidates(
            molecule,
            raw,
            config=run_config,
            source_software="mock-molscrub",
            source_command="mock molscrub emitted 20 protomers",
            output_dir=output_dir,
        )

    monkeypatch.setattr(engine_module, "generate_protomer_candidates", fake_protomers)
    monkeypatch.setattr(engine_module, "filter_tautomers_with_auto3d", _mock_tautomer_filter)
    monkeypatch.setattr(engine_module, "enumerate_stereoisomers", _mock_stereo_enumeration)
    monkeypatch.setattr(engine_module, "filter_stereoisomers_with_auto3d", _mock_stereo_energy_filter)
    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", _mock_final_auto3d)

    def fail_if_crest_runs(*args, **kwargs):
        raise AssertionError("CREST/xTB must not run when optional validation is disabled")

    monkeypatch.setattr(engine_module, "run_crest_for_seed", fail_if_crest_runs)

    result = run_workflow(config)

    assert result.molecule_count == 1
    protomers = _read_csv(outdir / "enumeration" / "protomers" / "protomers_selected.csv")
    tautomers = _read_csv(outdir / "enumeration" / "tautomers" / "tautomers_selected.csv")
    stereos = _read_csv(outdir / "enumeration" / "stereoisomers" / "stereoisomers_selected.csv")
    final_variants = _read_csv(outdir / "final_variants.csv")
    assert len(protomers) <= config.protonation.max_protomers_per_molecule
    assert len(tautomers) <= len(protomers) * config.tautomer_filtering.tauto_k
    assert len(stereos) <= len(tautomers) * config.stereoisomer_filtering.max_stereoisomers_per_tautomer
    assert len(final_variants) <= (
        config.protonation.max_protomers_per_molecule
        * config.tautomer_filtering.tauto_k
        * config.stereoisomer_filtering.max_stereoisomers_per_tautomer
    )

    for protomer in protomers:
        children = [row for row in tautomers if row["parent_protomer_id"] == protomer["protomer_id"]]
        assert len(children) <= config.tautomer_filtering.tauto_k
        relatives = [float(row["relative_energy_kcal_mol"]) for row in children]
        assert all(value <= config.tautomer_filtering.tauto_window_kcal_mol for value in relatives)

    for tautomer in tautomers:
        children = [row for row in stereos if row["parent_tautomer_id"] == tautomer["tautomer_id"]]
        assert len(children) <= config.stereoisomer_filtering.max_stereoisomers_per_tautomer

    assert not list(outdir.rglob("*.xyz"))
    assert not (outdir / "crest_validation.csv").exists()
    manifest = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["optional_validation"]["crest_xtb_enabled"] is False
    assert manifest["filtering"]["final_3d"]["final_conformer_count"] == len(final_variants)
    _assert_rejected_reasons(outdir / "enumeration" / "protomers" / "protomers_rejected.csv", "selection_reason")
    _assert_rejected_reasons(outdir / "enumeration" / "tautomers" / "tautomers_rejected.csv", "reason")
    _assert_rejected_reasons(outdir / "enumeration" / "stereoisomers" / "stereoisomers_rejected.csv", "reason")


def test_rdkit_tautomer_timeout_fallback_is_bounded(tmp_path: Path, monkeypatch) -> None:
    def timeout(*args, **kwargs):
        raise RdkitTautomerFilteringTimeout("mock timeout")

    def fail_auto3d(*args, **kwargs):
        raise AssertionError("single timeout fallback tautomer must not invoke Auto3D")

    monkeypatch.setattr(tautomer_filter, "_enumerate_molblocks_with_timeout", timeout)
    monkeypatch.setattr(tautomer_filter, "run_auto3d", fail_auto3d)
    records = filter_tautomers_with_auto3d([_protomer_record("CCO")], RunConfig(output_dir=tmp_path / "run"))

    assert len(records) == 1
    assert any("TAUTOMER_TIMEOUT_FALLBACK" in warning for warning in records[0].warnings)


def test_stereo_timeout_fallback_is_bounded(tmp_path: Path, monkeypatch) -> None:
    def timeout(*args, **kwargs):
        raise TimeoutError("mock stereo timeout")

    monkeypatch.setattr(stereochemistry, "_enumerate_with_timeout", timeout)
    records = enumerate_stereoisomers(_tautomer_record(_protomer_record("CCO"), 1), RunConfig(output_dir=tmp_path / "run"))

    assert len(records) == 1
    assert any("STEREO_TIMEOUT_FALLBACK" in warning for warning in records[0].warnings)


def test_final_auto3d_gpu_failure_retries_cpu(tmp_path: Path, monkeypatch) -> None:
    calls: list[bool] = []

    def fake_run_auto3d(input_path: Path, output_dir: Path, *, use_gpu: bool, **kwargs):
        calls.append(use_gpu)
        if use_gpu:
            raise Auto3DExecutionError("mock CUDA failure")
        return _write_auto3d_output(input_path, output_dir, energy=-1.0), ["auto3d", "cpu"]

    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", fake_run_auto3d)
    result = generate_final_3d_variants([_stereo_record(_tautomer_record(_protomer_record("CCO"), 1), 1)], RunConfig(output_dir=tmp_path / "run"))

    assert calls[:2] == [True, True]
    assert calls[-1] is False
    assert result.used_fallback is False
    assert result.records[0].source_software == "auto3d"


def test_final_auto3d_batch_failure_retries_smaller_batches(tmp_path: Path, monkeypatch) -> None:
    batch_sizes: list[int] = []

    def fake_run_auto3d(input_path: Path, output_dir: Path, **kwargs):
        size = _sdf_count(input_path)
        batch_sizes.append(size)
        if size > 1:
            raise Auto3DExecutionError("mock batch failure")
        return _write_auto3d_output(input_path, output_dir, energy=-2.0), ["auto3d", "single"]

    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", fake_run_auto3d)
    protomer = _protomer_record("CCO")
    tautomer = _tautomer_record(protomer, 1)
    stereos = [_stereo_record(tautomer, 1), _stereo_record(tautomer, 2)]
    config = RunConfig(output_dir=tmp_path / "run", final_3d={"use_gpu": False})

    result = generate_final_3d_variants(stereos, config)

    assert batch_sizes[:2] == [2, 2]
    assert batch_sizes[-2:] == [1, 1]
    assert len(result.records) == 2
    assert (tmp_path / "run" / "final_3d" / "smaller_batches").exists()


def test_local_agent_default_disabled_and_mocked_enabled_path(monkeypatch) -> None:
    assert AgentConfig().enabled is False
    monkeypatch.setattr("dsvr.agent.local_qwen.shutil.which", lambda _cmd: "/usr/bin/codex")

    def fake_run(command, **kwargs):
        assert command == ["codex", "--oss", "-m", "qwen3.6:35b"]
        assert "Allowed actions:" in kwargs["input"]
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="ACTION: retry_auto3d_cpu\n- mocked diagnostic\n",
            stderr="",
        )

    monkeypatch.setattr("dsvr.agent.local_qwen.subprocess.run", fake_run)
    result = run_local_diagnostic_agent(
        agent=AgentConfig(enabled=True),
        task="classify_failure",
        bug_context="Auto3D GPU failed",
    )

    assert result.available is True
    assert result.decision.action == "retry_auto3d_cpu"
    assert result.decision.valid is True


def _mock_tautomer_filter(protomers: list[ProtomerRecord], config: RunConfig) -> list[TautomerRecord]:
    outdir = config.output_dir / "enumeration" / "tautomers"
    outdir.mkdir(parents=True, exist_ok=True)
    selected_records: list[TautomerRecord] = []
    all_rows = []
    selected_rows = []
    rejected_rows = []
    for protomer in protomers:
        for index in range(1, 201):
            record = _tautomer_record(protomer, index)
            energy = float(index - 1)
            selected = index <= config.tautomer_filtering.tauto_k and energy <= config.tautomer_filtering.tauto_window_kcal_mol
            reason = "selected_by_mock_auto3d_energy" if selected else "rejected_by_mock_auto3d_energy"
            row = {
                "parent_protomer_id": protomer.id,
                "tautomer_id": record.id,
                "input_molecule_id": protomer.input_molecule_id,
                "selected": selected,
                "reason": reason,
                "source": "mock-auto3d",
                "auto3d_rank": index,
                "energy_kcal_mol": energy,
                "relative_energy_kcal_mol": energy,
                "warnings": "mock RDKit emitted 200 tautomers per protomer",
            }
            all_rows.append(row)
            if selected:
                selected_records.append(record)
                selected_rows.append(row)
            else:
                rejected_rows.append(row)
    _write_csv(outdir / "tautomers_all_pre_auto3d.csv", all_rows)
    _write_csv(outdir / "tautomers_auto3d_ranked.csv", selected_rows + rejected_rows)
    _write_csv(outdir / "tautomers_selected.csv", selected_rows)
    _write_csv(outdir / "tautomers_rejected.csv", rejected_rows)
    return selected_records


def _mock_stereo_enumeration(tautomer: TautomerRecord, config: RunConfig) -> list[StereoRecord]:
    outdir = config.output_dir / "enumeration" / "stereoisomers"
    outdir.mkdir(parents=True, exist_ok=True)
    selected_records = []
    all_rows = []
    selected_rows = []
    rejected_rows = []
    for index in range(1, 65):
        record = _stereo_record(tautomer, index)
        selected = index <= config.stereoisomer_filtering.max_stereoisomers_per_tautomer
        reason = "selected_within_max_stereoisomers_per_tautomer" if selected else "rejected_beyond_max_stereoisomers_per_tautomer"
        row = {
            "parent_tautomer_id": tautomer.id,
            "stereo_id": record.id,
            "input_molecule_id": tautomer.input_molecule_id,
            "selected": selected,
            "reason": reason,
            "warnings": "mock RDKit emitted 64 stereoisomers per tautomer",
        }
        all_rows.append(row)
        if selected:
            selected_records.append(record)
            selected_rows.append(row)
        else:
            rejected_rows.append(row)
    _append_csv(outdir / "stereoisomers_all.csv", all_rows)
    _append_csv(outdir / "stereoisomers_selected.csv", selected_rows)
    _append_csv(outdir / "stereoisomers_rejected.csv", rejected_rows)
    return selected_records


def _mock_stereo_energy_filter(records: list[StereoRecord], config: RunConfig):
    return SimpleNamespace(
        all_records=records,
        selected_records=records,
        rejected_records=[],
        decisions=[],
        collapsed_count=0,
        energy_evaluation_count=len(records),
    )


def _mock_final_auto3d(input_path: Path, output_dir: Path, **kwargs):
    return _write_auto3d_output(input_path, output_dir, energy=-5.0), ["auto3d", "mock-final"]


def _protomer_record(smiles: str) -> ProtomerRecord:
    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    input_id = make_input_id("mol", canonical)
    return ProtomerRecord(
        id=make_protomer_id(input_id, 1, canonical, isomeric),
        parent_id=input_id,
        input_molecule_id=input_id,
        molname="mol",
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=0,
        source_software="test",
        protomer_index=1,
        rdkit_mol=mol,
    )


def _tautomer_record(protomer: ProtomerRecord, index: int) -> TautomerRecord:
    mol = Chem.MolFromSmiles("CCO")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"mock_tautomer_index": index}
    return TautomerRecord(
        id=make_tautomer_id(protomer.id, index, canonical, isomeric, metadata),
        parent_id=protomer.id,
        input_molecule_id=protomer.input_molecule_id,
        molname=protomer.molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=0,
        source_software="test",
        tautomer_index=index,
        rdkit_mol=mol,
        metadata=metadata,
    )


def _stereo_record(tautomer: TautomerRecord, index: int) -> StereoRecord:
    mol = Chem.MolFromSmiles("CCO")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"mock_stereo_index": index}
    return StereoRecord(
        id=make_stereo_id(tautomer.id, index, canonical, isomeric, metadata),
        parent_id=tautomer.id,
        input_molecule_id=tautomer.input_molecule_id,
        molname=tautomer.molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=0,
        source_software="test",
        stereo_index=index,
        rdkit_mol=mol,
        metadata=metadata,
    )


def _write_auto3d_output(input_path: Path, output_dir: Path, *, energy: float) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = output_dir / "mock_auto3d.sdf"
    writer = Chem.SDWriter(str(output_sdf))
    for index, mol in enumerate(Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False), start=1):
        if mol is None:
            continue
        mol.SetProp("E_kcal_mol", str(energy - index))
        writer.write(mol)
    writer.close()
    return output_sdf


def _sdf_count(path: Path) -> int:
    return sum(1 for mol in Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False) if mol is not None)


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _append_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def _assert_rejected_reasons(path: Path, column: str) -> None:
    rows = _read_csv(path)
    assert rows, path
    assert all(row.get(column, "").strip() for row in rows)
