from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.chemistry import stereo_auto3d_filter
from dsvr.chemistry.stereo_auto3d_filter import filter_stereoisomers_with_auto3d
from dsvr.config import RunConfig
from dsvr.models import StereoRecord


def test_enantiomer_pair_collapses_to_one_auto3d_energy_evaluation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereos = [
        _stereo("F[C@H](Cl)Br", "r", 1),
        _stereo("F[C@@H](Cl)Br", "s", 2),
    ]
    calls: list[list[str]] = []

    def fake_run_auto3d(input_path, output_dir, **kwargs):
        ids = [line.split()[1] for line in input_path.read_text(encoding="utf-8").splitlines()]
        calls.append(ids)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "auto3d_output.sdf"
        _write_energy_sdf(output_sdf, [stereos[0]], {stereos[0].id: 0.0})
        return output_sdf, ["auto3d", "mock"]

    monkeypatch.setattr(stereo_auto3d_filter, "run_auto3d", fake_run_auto3d)
    result = filter_stereoisomers_with_auto3d(
        stereos,
        RunConfig(input_path=tmp_path / "in.smi", output_dir=tmp_path / "run"),
    )

    assert calls == [[stereos[0].id]]
    assert result.energy_evaluation_count == 1
    assert result.collapsed_count == 1
    assert {record.id for record in result.selected_records} == {record.id for record in stereos}
    assert (tmp_path / "run" / "stereo_enantiomer_groups.csv").exists()


def test_diastereomers_are_not_collapsed_and_high_energy_is_rejected(
    tmp_path: Path,
    monkeypatch,
) -> None:
    stereos = [
        _stereo("C[C@H](O)[C@H](F)Cl", "rr", 1),
        _stereo("C[C@H](O)[C@@H](F)Cl", "rs", 2),
    ]

    def fake_run_auto3d(input_path, output_dir, **kwargs):
        ids = [line.split()[1] for line in input_path.read_text(encoding="utf-8").splitlines()]
        assert ids == [stereos[0].id, stereos[1].id]
        output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = output_dir / "auto3d_output.sdf"
        _write_energy_sdf(output_sdf, stereos, {stereos[0].id: 0.0, stereos[1].id: 12.0})
        return output_sdf, ["auto3d", "mock"]

    monkeypatch.setattr(stereo_auto3d_filter, "run_auto3d", fake_run_auto3d)
    result = filter_stereoisomers_with_auto3d(
        stereos,
        RunConfig(
            input_path=tmp_path / "in.smi",
            output_dir=tmp_path / "run",
            stereoisomer_filtering={"stereo_energy_window_kcal_mol": 7.0},
        ),
    )

    assert result.energy_evaluation_count == 2
    assert result.collapsed_count == 0
    assert [record.id for record in result.selected_records] == [stereos[0].id]
    assert [record.id for record in result.rejected_records] == [stereos[1].id]
    assert (tmp_path / "run" / "stereoisomers_all.csv").exists()
    assert (tmp_path / "run" / "stereoisomers_selected.csv").exists()
    assert (tmp_path / "run" / "stereoisomers_rejected.csv").exists()
    assert (tmp_path / "run" / "stereo_energy_ranked.csv").exists()


def _write_energy_sdf(path: Path, records: list[StereoRecord], energies: dict[str, float]) -> None:
    writer = Chem.SDWriter(str(path))
    for record in records:
        mol = Chem.Mol(record.rdkit_mol)
        mol.SetProp("_Name", record.id)
        mol.SetProp("DSVR_STEREO_ID", record.id)
        mol.SetProp("E_kcal_mol", str(energies[record.id]))
        writer.write(mol)
    writer.close()


def _stereo(smiles: str, suffix: str, stereo_index: int) -> StereoRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return StereoRecord(
        id=f"mol_p01_t01_s{stereo_index:02d}_{suffix}",
        parent_id="mol_p01_t01",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=0,
        source_software="test",
        source_python_function="test",
        stereo_index=stereo_index,
        rdkit_mol=mol,
    )
