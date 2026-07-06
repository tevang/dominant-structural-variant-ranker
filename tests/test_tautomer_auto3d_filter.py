from pathlib import Path

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
    assert "RDKit tautomer enumeration timeout" in selected
