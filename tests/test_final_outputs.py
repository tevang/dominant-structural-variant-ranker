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
    assert "Ranking uses relative CREST/xTB-derived free energies" in report
    assert "Population scope is 'same_formula'" in report
    assert "micro-pKa/proton chemical-potential corrections" in report
