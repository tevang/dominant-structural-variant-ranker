import json
from pathlib import Path

import pandas as pd
import pytest
from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.io.write_outputs import RANKED_VARIANT_COLUMNS, SDF_RANKED_PROPERTIES
from dsvr.models import CrestConformerRecord
from dsvr.reporting.audit import VARIANT_DECISION_COLUMNS
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
    report_text_pre = (outdir / "report.md").read_text(encoding="utf-8")
    assert "- Molecules read: 1" in report_text_pre
    assert "- Molecules failed: 1" in report_text_pre

    stage_summary = pd.read_csv(outdir / "stage_summary.csv")
    input_row = stage_summary.loc[stage_summary["stage"] == "Input validation"].iloc[0]
    assert input_row["accepted_count"] == 1
    assert input_row["rejected_count"] == 1
    assert input_row["generated_count"] == 2

    warnings = [
        json.loads(line)
        for line in (outdir / "warnings.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    rejections = [
        item
        for item in warnings
        if item["level"] == "warning" and item["stage"] == "Input validation"
    ]
    assert len(rejections) == 1
    assert "\n" not in rejections[0]["message"]
    assert "not_a_smiles" in rejections[0]["message"]

    invalid_frame = pd.read_csv(outdir / "invalid_inputs.csv")
    assert list(invalid_frame.columns) == [
        "input_id", "source_format", "line_number", "name", "smiles", "raw_record", "error",
    ]
    assert len(invalid_frame) == 1
    assert invalid_frame.iloc[0]["smiles"] == "not_a_smiles"
    assert invalid_frame.iloc[0]["name"] == "bad"

    inputs_frame = pd.read_csv(outdir / "inputs.csv")
    assert len(inputs_frame) == 1
    assert inputs_frame.iloc[0]["molname"] == "ethanol"
    assert inputs_frame.iloc[0]["input_molecule_id"] == "mol_000001"
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



def test_clean_run_overwrites_stale_invalid_inputs_with_fresh_header(tmp_path: Path) -> None:
    input_path = tmp_path / "mols.smi"
    input_path.write_text("CCO ethanol\nCCN ethylamine\n", encoding="utf-8")
    outdir = tmp_path / "run"
    outdir.mkdir()
    (outdir / "invalid_inputs.csv").write_text(
        "record_index,raw_record,error\n0,stale,boom\n", encoding="utf-8"
    )

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

    invalid_text = (outdir / "invalid_inputs.csv").read_text(encoding="utf-8")
    assert "stale" not in invalid_text
    assert invalid_text == (
        "input_id,source_format,line_number,name,smiles,raw_record,error\n"
    )

    stage_summary = pd.read_csv(outdir / "stage_summary.csv")
    input_row = stage_summary.loc[stage_summary["stage"] == "Input validation"].iloc[0]
    assert input_row["accepted_count"] == 2
    assert input_row["rejected_count"] == 0

    warnings_text = (outdir / "warnings.jsonl").read_text(encoding="utf-8")
    assert "Rejected input" not in warnings_text

    inputs_frame = pd.read_csv(outdir / "inputs.csv")
    assert len(inputs_frame) == 2
    report = (outdir / "report.md").read_text(encoding="utf-8")
    assert "- Molecules read: 2" in report
    assert "- Molecules failed: 0" in report



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


def test_extract_energy_prefers_total_energy_and_converts_hartree():
    """Regression test: Auto3D v3 writes absolute energy as E_tot/E_tot(Hartree)
    in Hartree, while E_rel(kcal/mol) is 0.0 for the selected best-of-k
    conformer. The extractor must return the converted absolute energy so
    variant ranking keeps discrimination."""
    from rdkit import Chem

    from dsvr.chemistry.final3d import _HARTREE_TO_KCAL_MOL, _extract_energy

    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("E_rel(kcal/mol)", "0.0")
    mol.SetProp("E_tot(Hartree)", "-1360.1065899630466")

    energy, prop = _extract_energy(mol)
    assert prop == "E_tot(Hartree)"
    assert energy == pytest.approx(-1360.1065899630466 * _HARTREE_TO_KCAL_MOL)

    mol2 = Chem.MolFromSmiles("CCO")
    mol2.SetProp("E_tot(kcal/mol)", "-42.5")
    energy2, prop2 = _extract_energy(mol2)
    assert prop2 == "E_tot(kcal/mol)"
    assert energy2 == pytest.approx(-42.5)

    mol3 = Chem.MolFromSmiles("CCO")
    energy3, prop3 = _extract_energy(mol3)
    assert (energy3, prop3) == (None, None)


def test_final3d_fallback_engine_receives_only_supported_molecules(tmp_path, monkeypatch):
    """Escalation-review must-fix: when the primary engine fails and the
    chain advances, the fallback engine must get ONLY the molecules it
    supports, not the full group SDF (this prevents the 'Only AIMNET can
    handle' retry storm)."""

    from dsvr.chemistry import final3d
    from dsvr.chemistry.final3d import _run_final_auto3d_for_engine
    from dsvr.runners.auto3d_runner import Auto3DExecutionError

    primary = "AIMNET"
    fallback_engine = "ANI2xt"
    config = RunConfig(
        output_dir=tmp_path,
        final_3d={
            "optimizing_engine": primary,
            "fallback_optimizing_engine": fallback_engine,
            "use_gpu": False,
        },
    )

    group_dir = tmp_path / "group"
    group_dir.mkdir(parents=True)
    writer = Chem.SDWriter(str(group_dir / "input.sdf"))
    for smiles in ("CCO", "c1ccc([O-])cc1"):
        mol = Chem.MolFromSmiles(smiles)
        writer.write(mol)
    writer.close()

    calls: list[tuple[str, list[str]]] = []

    def fake_run_auto3d(input_path, output_dir, **kwargs):
        size = sum(
            1
            for mol in Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False)
            if mol is not None
        )
        calls.append((kwargs["model"], [input_path.name] * size + [f"n={size}"]))
        if kwargs["model"] == primary:
            raise Auto3DExecutionError("mock primary failure")
        output_sdf = Path(output_dir) / "out.sdf"
        writer2 = Chem.SDWriter(str(output_sdf))
        for mol in Chem.SDMolSupplier(str(input_path), sanitize=True, removeHs=False):
            if mol is not None:
                mol.SetProp("E_kcal_mol", "-1.0")
                writer2.write(mol)
        writer2.close()
        return output_sdf, ["auto3d", "run"]

    monkeypatch.setattr(final3d, "run_auto3d", fake_run_auto3d)

    sdf, _command = _run_final_auto3d_for_engine(group_dir / "input.sdf", group_dir, config, primary)

    assert calls[1][0] == fallback_engine
    # Only the neutral molecule is supported by ANI2xt; the phenolate is not offered.
    assert calls[1][1] == ["final_3d_input_ANI2xt.sdf", "n=1"]
    assert sdf.exists()
