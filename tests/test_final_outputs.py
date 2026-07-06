import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.io.write_outputs import RANKED_VARIANT_COLUMNS, SDF_RANKED_PROPERTIES
from dsvr.reporting.audit import VARIANT_DECISION_COLUMNS
from dsvr.models import CrestConformerRecord
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
    assert (outdir / "variant_decisions.csv").exists()
    for name in (
        "protomers_all.csv",
        "protomers_selected.csv",
        "protomers_rejected.csv",
        "tautomers_all_pre_auto3d.csv",
        "tautomers_auto3d_ranked.csv",
        "tautomers_selected.csv",
        "tautomers_rejected.csv",
        "stereoisomers_all.csv",
        "stereoisomers_selected.csv",
        "stereoisomers_rejected.csv",
    ):
        assert (outdir / name).exists(), name

    decisions = pd.read_csv(outdir / "variant_decisions.csv")
    assert set(VARIANT_DECISION_COLUMNS).issubset(decisions.columns)
    assert "final_variant" in set(decisions["stage"])
    assert decisions["rejection_reason"].fillna("").map(type).eq(str).all()

    report = (outdir / "report.md").read_text(encoding="utf-8")
    assert "Concise Audit Summary" in report
    assert "Molecules read:" in report
    assert "Final variants written:" in report
    assert "Agent interventions enabled:" in report
    assert "Optional validation results enabled:" in report
    assert not any((outdir / "crest").glob("*/crest_provenance.jsonl"))



def test_optional_crest_validation_writes_separate_outputs(tmp_path: Path, monkeypatch) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\n", encoding="utf-8")
    outdir = tmp_path / "run"
    crest_calls = []

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
        supplier = Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
        output_sdf = output_dir / "mock_final_auto3d.sdf"
        writer = Chem.SDWriter(str(output_sdf))
        for mol in supplier:
            if mol is None:
                continue
            mol.SetProp("E_kcal_mol", "-3.25")
            writer.write(mol)
        writer.close()
        return output_sdf, ["auto3d", "mock", "--optimizing_engine", model]

    def fake_run_crest_for_seed(seed, config):
        crest_calls.append((seed, config))
        return [
            CrestConformerRecord(
                id=f"{seed.id}_crest_validation_0001",
                parent_id=seed.id,
                input_molecule_id=seed.input_molecule_id,
                molname=seed.molname,
                canonical_smiles=seed.canonical_smiles,
                isomeric_smiles=seed.isomeric_smiles,
                molecular_formula=seed.molecular_formula,
                formal_charge=seed.formal_charge,
                explicit_proton_count=seed.explicit_proton_count,
                source_software="crest",
                source_python_function="test.fake_run_crest_for_seed",
                warnings=["optional validation mock"],
                metadata={"crest": {"workdir": str(outdir / "crest" / seed.id)}},
                crest_index=1,
                energy_kcal_mol=-4.0,
                relative_energy_kcal_mol=0.0,
            )
        ]

    monkeypatch.setattr("dsvr.chemistry.final3d.run_auto3d", fake_run_auto3d)
    monkeypatch.setattr("dsvr.workflow.engine._tool_available", lambda executable: True)
    monkeypatch.setattr("dsvr.workflow.engine.run_crest_for_seed", fake_run_crest_for_seed)

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
            optional_validation={"crest_xtb_enabled": True},
        )
    )

    assert len(crest_calls) == 1
    validation_csv = outdir / "crest_validation.csv"
    validation_sdf = outdir / "crest_validation.sdf"
    validation_report = outdir / "crest_validation_report.md"
    assert validation_csv.exists()
    assert validation_sdf.exists()
    assert validation_report.exists()
    assert (outdir / "optional_validation" / "selected_final_variants.sdf").exists()

    validation_frame = pd.read_csv(validation_csv)
    assert len(validation_frame) == 1
    assert bool(validation_frame.loc[0, "optional_validation"]) is True
    assert validation_frame.loc[0, "crest_energy_kcal_mol"] == -4.0

    mols = [mol for mol in Chem.SDMolSupplier(str(validation_sdf), sanitize=True, removeHs=False) if mol]
    assert len(mols) == 1
    assert mols[0].GetProp("DSVR_OPTIONAL_VALIDATION") == "CREST/xTB"
    assert mols[0].GetProp("DSVR_OPTIONAL_VALIDATION_DOES_NOT_SET_RANKING") == "True"

    manifest = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["optional_validation"]["crest_xtb_enabled"] is True
    assert manifest["optional_validation"]["selected_count"] == 1
    assert manifest["optional_validation"]["ranking_overwritten"] is False
    assert (outdir / "ranked_variants.csv").exists()
