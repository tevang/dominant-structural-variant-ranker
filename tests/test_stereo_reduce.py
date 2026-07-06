from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.config import RunConfig
from dsvr.filtering.stereo_reduce import (
    expand_enantiomer_mapped_crest_records,
    reduce_seeds_for_crest,
)
from dsvr.models import CrestConformerRecord, SeedConformerRecord, StereoRecord
from dsvr.reporting.markdown import write_run_report


def test_enantiomer_pair_is_collapsed_to_one_crest_job() -> None:
    stereos = [
        _stereo("F[C@H](Cl)Br", "r", 1),
        _stereo("F[C@@H](Cl)Br", "s", 2),
    ]
    seeds = [_seed(record) for record in stereos]

    reduction = reduce_seeds_for_crest(seeds, stereos, RunConfig())

    assert len(reduction.selected_seeds) == 1
    assert reduction.jobs_saved == 1
    assert any(
        decision.reason == "crest_skipped_for_enantiomer_equivalent_in_achiral_solvent"
        for decision in reduction.decisions
    )


def test_diastereomers_remain_distinct_for_crest() -> None:
    stereos = [
        _stereo("C[C@H](O)[C@H](F)Cl", "rr", 1),
        _stereo("C[C@H](O)[C@@H](F)Cl", "rs", 2),
    ]
    seeds = [_seed(record) for record in stereos]

    reduction = reduce_seeds_for_crest(seeds, stereos, RunConfig())

    assert len(reduction.selected_seeds) == 2
    assert reduction.jobs_saved == 0
    assert all(decision.selected_for_crest for decision in reduction.decisions)


def test_expanded_crest_records_preserve_all_stereoisomer_identities() -> None:
    stereos = [
        _stereo("F[C@H](Cl)Br", "r", 1),
        _stereo("F[C@@H](Cl)Br", "s", 2),
    ]
    seeds = [_seed(record) for record in stereos]
    reduction = reduce_seeds_for_crest(seeds, stereos, RunConfig())
    representative_seed = reduction.selected_seeds[0]
    crest_record = _crest(representative_seed)

    expanded = expand_enantiomer_mapped_crest_records(
        [crest_record],
        {seed.id: seed for seed in seeds},
        reduction,
        RunConfig(),
    )

    assert len(expanded) == 2
    assert {record.isomeric_smiles for record in expanded} == {
        record.isomeric_smiles for record in seeds
    }
    mapped = [record for record in expanded if record.id != crest_record.id]
    assert mapped
    assert "stereo_reduction" in mapped[0].metadata
    assert any(
        "mapped from an enantiomeric representative" in warning
        for warning in mapped[0].warnings
    )


def test_report_states_stereo_collapse_jobs_saved(tmp_path: Path) -> None:
    report_path = tmp_path / "report.md"

    write_run_report(
        report_path,
        config=RunConfig(output_dir=tmp_path),
        records=[],
        ranked_records=[],
        manifest={
            "filtering": {
                "stereo_energy_filtering": {
                    "enumerated_count": 6,
                    "selected_count": 4,
                    "rejected_count": 2,
                    "collapsed_count": 1,
                    "energy_evaluation_count": 5,
                },
                "stereo_reduction": {
                    "jobs_saved": 3,
                    "enabled": True,
                    "decision_count": 4,
                }
            }
        },
        output_files=[],
    )

    text = report_path.read_text(encoding="utf-8")
    assert "Enumerated stereo states: 6" in text
    assert "Selected stereo states: 4" in text
    assert "Rejected stereo states: 2" in text
    assert "Enantiomer states collapsed for Auto3D energy evaluation: 1" in text
    assert "Auto3D stereo energy evaluations run: 5" in text
    assert "CREST jobs saved by enantiomer collapse: 3" in text


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


def _seed(stereo: StereoRecord) -> SeedConformerRecord:
    return SeedConformerRecord(
        id=f"{stereo.id}_c01",
        parent_id=stereo.id,
        input_molecule_id=stereo.input_molecule_id,
        molname=stereo.molname,
        canonical_smiles=stereo.canonical_smiles,
        isomeric_smiles=stereo.isomeric_smiles,
        molecular_formula=stereo.molecular_formula,
        formal_charge=stereo.formal_charge,
        explicit_proton_count=stereo.explicit_proton_count,
        source_software="test",
        source_python_function="test",
        conformer_index=1,
        energy_kcal_mol=0.0,
        rdkit_mol=stereo.rdkit_mol,
    )


def _crest(seed: SeedConformerRecord) -> CrestConformerRecord:
    return CrestConformerRecord(
        id=f"{seed.id}_crest01",
        parent_id=seed.id,
        input_molecule_id=seed.input_molecule_id,
        molname=seed.molname,
        canonical_smiles=seed.canonical_smiles,
        isomeric_smiles=seed.isomeric_smiles,
        molecular_formula=seed.molecular_formula,
        formal_charge=seed.formal_charge,
        explicit_proton_count=seed.explicit_proton_count,
        source_software="test",
        source_python_function="test",
        crest_index=1,
        energy_kcal_mol=1.0,
    )
