from pathlib import Path

import pandas as pd
from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.io.write_outputs import RANKED_VARIANT_COLUMNS, SDF_RANKED_PROPERTIES
from dsvr.workflow.engine import run_workflow


def test_final_ranked_outputs_have_required_columns_and_sdf_properties(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"
    config = RunConfig(
        input_path=input_path,
        output_dir=outdir,
        overwrite=True,
        protonation={"enabled": False},
        tautomer_filtering={"enabled": False},
        stereoisomer_filtering={"enabled": False},
        enumeration={
            "max_protomers_per_molecule": 1,
            "max_tautomers_per_protomer": 1,
            "max_stereoisomers_per_tautomer": 1,
        },
        seeding={"rdkit_num_conformers": 1},
        crest={"enabled": False},
        thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
    )

    run_workflow(config)

    frame = pd.read_csv(outdir / "ranked_variants.csv")
    assert set(RANKED_VARIANT_COLUMNS).issubset(frame.columns)
    assert len(frame) >= 1

    supplier = Chem.SDMolSupplier(str(outdir / "ranked_variants.sdf"), sanitize=True)
    mols = [mol for mol in supplier if mol is not None]
    assert mols
    for prop in SDF_RANKED_PROPERTIES:
        assert mols[0].HasProp(prop), prop


def test_report_generated_even_with_invalid_input_records(tmp_path: Path) -> None:
    input_path = tmp_path / "mixed.smi"
    input_path.write_text("CCO ethanol\nnot_a_smiles bad\n", encoding="utf-8")
    outdir = tmp_path / "run"

    run_workflow(
        RunConfig(
            input_path=input_path,
            output_dir=outdir,
            overwrite=True,
            protonation={"enabled": False},
            tautomer_filtering={"enabled": False},
            stereoisomer_filtering={"enabled": False},
            enumeration={
                "max_protomers_per_molecule": 1,
                "max_tautomers_per_protomer": 1,
                "max_stereoisomers_per_tautomer": 1,
            },
            seeding={"rdkit_num_conformers": 1},
            crest={"enabled": False},
            thermo={"enabled": False, "xtb_hessian": False, "xtb_thermo": False},
        )
    )

    assert (outdir / "invalid_inputs.csv").exists()
    assert (outdir / "report.md").exists()
    report = (outdir / "report.md").read_text(encoding="utf-8")
    assert "Run Settings" in report
    assert "Tool Versions" in report
    assert "Failure Summary" in report
    assert "Output File Locations" in report
    assert "pH 7.0 is used for candidate generation" in report
    assert "Solvent 'water' with solvent model 'alpb'" in report
    assert "Ranking uses approximate final Auto3D conformer energies" in report
    assert "Population scope is 'same_formula'" in report
    assert "micro-pKa/proton chemical-potential corrections" in report



def test_default_ligprep_like_writes_final_auto3d_variants_without_crest(
    tmp_path: Path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"

    def fake_run_auto3d(
        input_path: Path,
        output_dir: Path,
        *,
        k: int,
        model: str,
        internal_tautomer_stereo_enum: bool,
        max_confs: int | None = None,
        patience: int | None = None,
        use_gpu: bool = False,
        timeout_s: int | None = None,
        **kwargs,
    ) -> tuple[Path, list[str]]:
        assert k == 1
        assert max_confs == 10
        assert patience == 100
        assert timeout_s == 1800
        assert internal_tautomer_stereo_enum is False
        supplier = Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
        output_sdf = output_dir / "mock_final_auto3d.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for mol in supplier:
            if mol is None:
                continue
            mol.SetProp("E_kcal_mol", "-12.5")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "mock", "--optimizing_engine", model]

    def fail_if_crest_runs(*args, **kwargs):
        raise AssertionError("CREST must not run in default final_3d mode")

    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", fake_run_auto3d)
    monkeypatch.setattr("dsvr.workflow.engine.run_crest_for_seed", fail_if_crest_runs)

    run_workflow(
        RunConfig(
            input_path=input_path,
            output_dir=outdir,
            overwrite=True,
            protonation={"enabled": False},
            tautomer_filtering={"enabled": False},
            stereoisomer_filtering={"enabled": False},
            enumeration={
                "max_protomers_per_molecule": 1,
                "max_tautomers_per_protomer": 1,
                "max_stereoisomers_per_tautomer": 1,
            },
        )
    )

    final_sdf = outdir / "final_variants.sdf"
    assert final_sdf.exists()
    mols = [mol for mol in Chem.SDMolSupplier(str(final_sdf), sanitize=True, removeHs=False) if mol]
    assert len(mols) == 1
    mol = mols[0]
    assert mol.HasProp("DSVR_FINAL_VARIANT_ID")
    assert mol.HasProp("DSVR_STEREO_ID")
    assert mol.GetProp("DSVR_FINAL_AUTO3D_ENERGY_KCAL_MOL") == "-12.5"
    assert mol.GetProp("DSVR_APPROXIMATE_RANKING") == "True"
    assert "not solvated free energies" in mol.GetProp("DSVR_ENERGY_WARNING")
    assert (outdir / "final_variants.csv").exists()
    assert (outdir / "final_variants.json").exists()
    assert (outdir / "final_variant_energies.csv").exists()
    assert (outdir / "ranked_variants.csv").exists()
    assert not any((outdir / "crest").glob("*/crest_provenance.jsonl"))
