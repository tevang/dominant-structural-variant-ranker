from pathlib import Path

from rdkit import Chem
from typer.testing import CliRunner

from dsvr.chemistry import conformers_rdkit
from dsvr.chemistry.conformers_rdkit import generate_rdkit_seeds
from dsvr.cli import app
from dsvr.config import RunConfig
from dsvr.models import (
    StereoRecord,
    make_input_id,
    make_protomer_id,
    make_stereo_id,
    make_tautomer_id,
)


def test_generate_rdkit_seeds_for_simple_molecule(tmp_path: Path) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={"rdkit_num_conformers": 3, "rdkit_forcefield": "uff"},
        disk={"keep_raw_xyz": True},
    )

    records = generate_rdkit_seeds(stereo, config)

    assert records
    assert all(record.parent_id == stereo.id for record in records)
    assert all(record.formal_charge == 0 for record in records)
    assert all(record.embedding_status == "success" for record in records)
    assert all(
        record.forcefield_status in {"uff_minimized", "uff_forcefield_unavailable"}
        for record in records
    )
    seed_dir = tmp_path / "run" / "seeding" / "rdkit"
    assert (seed_dir / f"{stereo.id}_seeds.sdf").exists()
    assert (seed_dir / f"{stereo.id}_seeds.csv").exists()
    assert list((seed_dir / "xyz").glob("*.xyz"))


def test_xyz_output_is_basic_xtb_crest_compatible(tmp_path: Path) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={"rdkit_num_conformers": 1},
        disk={"keep_raw_xyz": True},
    )

    generate_rdkit_seeds(stereo, config)

    xyz_path = next((tmp_path / "run" / "seeding" / "rdkit" / "xyz").glob("*.xyz"))
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    atom_count = int(lines[0])
    assert atom_count > 0
    assert len(lines) == atom_count + 2
    assert all(len(line.split()) == 4 for line in lines[2:])


def test_embedding_failure_is_recorded(tmp_path: Path, monkeypatch) -> None:
    stereo = _stereo("ethanol", "CCO")
    config = RunConfig(input_path=tmp_path / "stereo.sdf", output_dir=tmp_path / "run")

    def fail_embed(*args, **kwargs):
        raise RuntimeError("embedding exploded")

    monkeypatch.setattr(conformers_rdkit.AllChem, "EmbedMultipleConfs", fail_embed)
    records = generate_rdkit_seeds(stereo, config)

    assert len(records) == 1
    assert records[0].embedding_status == "failed"
    assert "embedding exploded" in records[0].warnings[0]
    assert (tmp_path / "run" / "seeding" / "rdkit" / f"{stereo.id}_seeds.csv").exists()


def test_etkdg_pool_is_reduced_to_seed_budget_without_xyz_tree(tmp_path: Path) -> None:
    stereo = _stereo("hexane", "CCCCCC")
    config = RunConfig(
        input_path=tmp_path / "stereo.sdf",
        output_dir=tmp_path / "run",
        seeding={"rdkit_num_conformers": 30, "rdkit_forcefield": "uff"},
        variant_filtering={"max_seeds_per_variant": 2},
    )

    records = generate_rdkit_seeds(stereo, config)

    assert len(records) <= 2
    seed_dir = tmp_path / "run" / "seeding" / "rdkit"
    assert (seed_dir / f"{stereo.id}_seed_selection.csv").exists()
    assert not (seed_dir / "xyz").exists()


def test_cli_seed_etkdg_from_stereo_sdf(tmp_path: Path) -> None:
    stereo = _stereo("ethanol", "CCO")
    stereo_sdf = tmp_path / "stereo.sdf"
    writer = Chem.SDWriter(str(stereo_sdf))
    mol = Chem.Mol(stereo.rdkit_mol)
    mol.SetProp("_Name", stereo.id)
    mol.SetProp("DSVR_STEREO_ID", stereo.id)
    mol.SetProp("DSVR_INPUT_ID", stereo.input_molecule_id)
    mol.SetProp("DSVR_PARENT_TAUTOMER_ID", stereo.parent_id or "")
    mol.SetProp("DSVR_MOLNAME", stereo.molname)
    writer.write(mol)
    writer.close()
    outdir = tmp_path / "out"

    result = CliRunner().invoke(
        app,
        [
            "seed-etkdg",
            str(stereo_sdf),
            "--out",
            str(outdir),
            "--num-conformers",
            "2",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (outdir / "seeding" / "rdkit" / "seed_report.json").exists()
    assert list((outdir / "seeding" / "rdkit" / "xyz").glob("*.xyz"))


def _stereo(molname: str, smiles: str) -> StereoRecord:
    molecule = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
    isomeric = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    input_id = make_input_id(molname, canonical)
    protomer_id = make_protomer_id(input_id, 1, canonical, isomeric)
    tautomer_id = make_tautomer_id(protomer_id, 1, canonical, isomeric)
    stereo_id = make_stereo_id(tautomer_id, 1, canonical, isomeric)
    return StereoRecord(
        id=stereo_id,
        parent_id=tautomer_id,
        input_molecule_id=input_id,
        molname=molname,
        canonical_smiles=canonical,
        isomeric_smiles=isomeric,
        molecular_formula="",
        formal_charge=Chem.GetFormalCharge(molecule),
        explicit_proton_count=None,
        source_software="test",
        source_python_function="test",
        stereo_index=1,
        rdkit_mol=molecule,
    )
