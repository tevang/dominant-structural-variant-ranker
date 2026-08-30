import multiprocessing as mp
from pathlib import Path

import pytest
from rdkit import Chem

from dsvr.chemistry import tautomer_auto3d_filter as tautomer_filter
from dsvr.chemistry.tautomer_auto3d_filter import (
    RdkitTautomerFilteringTimeout,
    filter_tautomers_with_auto3d,
)
from dsvr.config import RunConfig
from dsvr.models import ProtomerRecord


def _protomer(smiles: str = "CC(=O)C") -> ProtomerRecord:
    mol = Chem.MolFromSmiles(smiles)
    return ProtomerRecord(
        id="mol_p01",
        parent_id="mol",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        molecular_formula="C3H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        protomer_index=1,
        rdkit_mol=mol,
    )


def test_auto3d_tautomer_filter_selects_ranked_survivors(
    tmp_path: Path,
    monkeypatch,
) -> None:
    def fake_auto3d(input_path: Path, output_dir: Path, **kwargs):
        assert kwargs["internal_tautomer_stereo_enum"] is False
        assert kwargs["model"] == "ANI2xt"
        output_dir.mkdir(parents=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        lines = [
            line.split(maxsplit=1)
            for line in input_path.read_text(encoding="utf-8").splitlines()
        ]
        energies = [-4.0, -1.0, 10.0]
        for index, (smiles, tautomer_id) in enumerate(lines):
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", tautomer_id)
            mol.SetProp("E_kcal_mol", str(energies[min(index, len(energies) - 1)]))
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", str(input_path)]

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake_auto3d)
    config = RunConfig(
        output_dir=tmp_path / "run",
        tautomer_filtering={"tauto_k": 1, "tauto_window_kcal_mol": 20.0},
    )

    records = filter_tautomers_with_auto3d([_protomer()], config)

    assert len(records) == 1
    tautomer_dir = tmp_path / "run" / "enumeration" / "tautomers"
    assert (tautomer_dir / "tautomers_all_pre_auto3d.csv").exists()
    assert (tautomer_dir / "tautomers_auto3d_ranked.csv").exists()
    assert (tautomer_dir / "tautomers_selected.csv").exists()
    assert (tautomer_dir / "tautomers_rejected.csv").exists()
    rejected = (tautomer_dir / "tautomers_rejected.csv").read_text(encoding="utf-8")
    assert "rejected_by_auto3d_energy_filter" in rejected
    selected_ids = {record.id for record in records}
    rejected_ids = {
        row.split(",")[1]
        for row in rejected.splitlines()[1:]
        if row.strip()
    }
    assert selected_ids.isdisjoint(rejected_ids)
    filtering = records[0].metadata["auto3d_tautomer_filtering"]
    assert filtering["score_is_population_estimate"] is False
    assert filtering["scope"] == "fast potential-energy tautomer filter before stereoisomer enumeration"


def test_rdkit_tautomer_timeout_falls_back_to_input(
    tmp_path: Path,
    monkeypatch,
) -> None:
    def timeout(*args, **kwargs):
        raise RdkitTautomerFilteringTimeout("mock timeout")

    def fail_auto3d(*args, **kwargs):
        raise AssertionError("Auto3D should not be needed for a single timeout fallback candidate")

    monkeypatch.setattr(tautomer_filter, "_enumerate_molblocks_with_timeout", timeout)
    monkeypatch.setattr(tautomer_filter, "run_auto3d", fail_auto3d)
    config = RunConfig(output_dir=tmp_path / "run")

    records = filter_tautomers_with_auto3d([_protomer()], config)

    assert len(records) == 1
    selected = (
        tmp_path / "run" / "enumeration" / "tautomers" / "tautomers_selected.csv"
    ).read_text(encoding="utf-8")
    assert "TAUTOMER_TIMEOUT_FALLBACK" in selected
    assert "RDKit tautomer enumeration timeout" in selected


def test_hanging_rdkit_tautomer_worker_is_killed(monkeypatch) -> None:
    def hanging_worker(*args, **kwargs):
        import time

        time.sleep(30)

    monkeypatch.setattr(tautomer_filter, "_tautomer_worker", hanging_worker)
    before = {process.pid for process in mp.active_children()}
    mol = Chem.MolFromSmiles("CC(=O)C")

    with pytest.raises(RdkitTautomerFilteringTimeout):
        tautomer_filter._enumerate_molblocks_with_timeout(
            mol,
            timeout_seconds=1,
            max_tautomers=4,
            max_transforms=8,
            remove_bond_stereo=True,
            remove_sp3_stereo=True,
            reassign_stereo=True,
        )

    for process in mp.active_children():
        process.join(timeout=0.1)
    after = {process.pid for process in mp.active_children()}
    assert after <= before


def test_rdkit_tautomer_cap_warning_is_preserved_with_auto3d(
    tmp_path: Path,
    monkeypatch,
) -> None:
    def capped_enumeration(molecule, **kwargs):
        molblock = Chem.MolToMolBlock(molecule)
        return tautomer_filter._EnumerationResult(
            molblocks=[molblock, molblock],
            warning="RDKit tautomer cap reached; ranked generated subset only",
        )

    def fake_auto3d(input_path: Path, output_dir: Path, **kwargs):
        output_dir.mkdir(parents=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for line in input_path.read_text(encoding="utf-8").splitlines():
            smiles, tautomer_id = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", tautomer_id)
            mol.SetProp("E_kcal_mol", "0.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", str(input_path)]

    monkeypatch.setattr(tautomer_filter, "_enumerate_molblocks_with_timeout", capped_enumeration)
    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake_auto3d)
    config = RunConfig(output_dir=tmp_path / "run", tautomer_filtering={"tauto_k": 1})

    records = filter_tautomers_with_auto3d([_protomer()], config)

    assert records
    assert any("RDKit tautomer cap reached" in warning for warning in records[0].warnings)
    ranked_path = (
        tmp_path
        / "run"
        / "enumeration"
        / "tautomers"
        / "tautomers_auto3d_ranked.csv"
    )
    ranked = ranked_path.read_text(encoding="utf-8")
    assert "RDKit tautomer cap reached" in ranked


def test_energy_from_mol_reads_auto3d_v3_total_energy():
    """Regression test: Auto3D v3 writes E_tot/E_tot(Hartree) in Hartree;
    the tautomer energy reader must convert it so ranking gets usable
    energies instead of silently falling back to RDKit."""
    from rdkit import Chem

    from dsvr.chemistry.tautomer_auto3d_filter import _energy_from_mol
    from dsvr.utils.units import HARTREE_TO_KCAL_MOL

    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("E_tot(Hartree)", "-153.0123")
    assert _energy_from_mol(mol) == pytest.approx(-153.0123 * HARTREE_TO_KCAL_MOL)

    mol2 = Chem.MolFromSmiles("CCO")
    mol2.SetProp("E_tot", "-153.0123")
    mol2.SetProp("E_rel(kcal/mol)", "0.0")
    # absolute total preferred over the near-zero relative conformer energy
    assert _energy_from_mol(mol2) == pytest.approx(-153.0123 * HARTREE_TO_KCAL_MOL)

    mol3 = Chem.MolFromSmiles("CCO")
    mol3.SetProp("E_rel(kcal/mol)", "1.5")
    assert _energy_from_mol(mol3) == pytest.approx(1.5)

    assert _energy_from_mol(Chem.MolFromSmiles("CCO")) is None


def _candidate(protomer: ProtomerRecord, index: int, smiles: str):
    from dsvr.chemistry.tautomer_auto3d_filter import _Candidate
    from dsvr.models import make_tautomer_id

    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    metadata = {"auto3d_tautomer_filtering": True, "tautomer_filtering_stage": "pre_stereo"}
    return _Candidate(
        index=index,
        tautomer_id=make_tautomer_id(protomer.id, index, canonical, isomeric, metadata),
        molecule=mol,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        rdkit_score=0.0,
        is_input_tautomer=index == 1,
        is_canonical_tautomer=index == 1,
    )


def test_mixed_engine_batch_is_split_per_engine(tmp_path, monkeypatch) -> None:
    """Task 2.3: an ANI2xt-configured batch containing an AIMNET-only
    (charged) candidate runs one sub-batch per engine."""

    protomer = _protomer("CC(=O)O")
    neutral = _candidate(protomer, 2, "C=C(O)O")
    charged = _candidate(protomer, 3, "CC(=O)[O-]")
    input_form = _candidate(protomer, 1, "CC(=O)O")
    candidates = [input_form, neutral, charged]
    monkeypatch.setattr(
        tautomer_filter, "_enumerate_candidates", lambda *a, **k: (candidates, None)
    )
    calls: list[tuple[str, list[str]]] = []

    def fake_auto3d(input_path: Path, output_dir: Path, **kwargs):
        lines = input_path.read_text(encoding="utf-8").splitlines()
        calls.append((kwargs["model"], [line.split()[0] for line in lines]))
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for line in lines:
            smiles, tautomer_id = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", tautomer_id)
            mol.SetProp("E_kcal_mol", "0.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run", "--engine", kwargs["model"]]

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake_auto3d)
    config = RunConfig(output_dir=tmp_path / "run")

    records = filter_tautomers_with_auto3d([protomer], config)

    engines_used = [model for model, _smiles in calls]
    assert engines_used == ["ANI2xt", "AIMNET"]
    assert len(calls[0][1]) == 2  # neutral sub-batch
    assert calls[1][1] == ["CC(=O)[O-]"]  # charged sub-batch only
    assert len(records) >= 1
    ranked_csv = (tmp_path / "run" / "enumeration" / "tautomers" / "tautomers_auto3d_ranked.csv").read_text(
        encoding="utf-8"
    )
    assert "cross-engine energy comparison is approximate" in ranked_csv


def test_candidate_unsupported_by_all_engines_uses_rdkit_fallback(tmp_path, monkeypatch) -> None:
    """Task 2.1/2.2: a molecule no configured engine supports is never offered
    to Auto3D and is retained via the recorded RDKit fallback."""

    protomer = _protomer("CC(=O)O")
    neutral = _candidate(protomer, 1, "CC(=O)O")
    # Lithium keeps the structure outside every engine's element set.
    lithium = _candidate(protomer, 2, "[Li+].[O-]C=C")
    candidates = [neutral, lithium]
    monkeypatch.setattr(
        tautomer_filter, "_enumerate_candidates", lambda *a, **k: (candidates, None)
    )
    models_seen: list[str] = []

    def fake_auto3d(input_path: Path, output_dir: Path, **kwargs):
        models_seen.append(kwargs["model"])
        lines = input_path.read_text(encoding="utf-8").splitlines()
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "mock.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for line in lines:
            smiles, tautomer_id = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", tautomer_id)
            mol.SetProp("E_kcal_mol", "0.0")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "run"]

    monkeypatch.setattr(tautomer_filter, "run_auto3d", fake_auto3d)
    config = RunConfig(output_dir=tmp_path / "run")

    records = filter_tautomers_with_auto3d([protomer], config)

    lithium_record = next((r for r in records if "Li" in r.isomeric_smiles), None)
    assert lithium_record is not None, records
    assert lithium_record.source_software == "rdkit_fallback"
    assert any("ENGINE_INCOMPATIBLE" in warning for warning in lithium_record.warnings)
    for model in models_seen:
        for smiles_line in ("[Li+]",):
            assert smiles_line not in model  # lithium never offered to any engine
    # One real invocation for the neutral candidate only
    assert len(models_seen) == 1
