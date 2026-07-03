from __future__ import annotations

from rdkit import Chem

from dsvr.config import RunConfig
from dsvr.filtering.selection import select_seed_records, select_stereo_records
from dsvr.models import SeedConformerRecord, StereoRecord


def test_balanced_filtering_caps_large_variant_counts() -> None:
    records = [_stereo(index) for index in range(1, 41)]
    config = RunConfig(
        variant_filtering={
            "mode": "balanced",
            "max_variants_before_3d_per_molecule": 10,
            "max_variants_after_cheap_score_per_molecule": 6,
            "min_variants_to_keep": 2,
            "absolute_penalty_cutoff": 100.0,
            "relative_penalty_cutoff": 100.0,
            "rescue_rules_enabled": False,
        }
    )

    selected, decisions = select_stereo_records(records, config, "pre_3d")
    selected, decisions_cheap = select_stereo_records(selected, config, "cheap_score")

    assert len(selected) == 6
    assert len([decision for decision in decisions if decision.selected]) == 10
    assert len([decision for decision in decisions_cheap if decision.selected]) == 6
    assert any(decision.reason == "pre_3d_over_budget" for decision in decisions)


def test_seed_filtering_caps_crest_fanout_per_molecule() -> None:
    records = []
    for variant_index in range(1, 11):
        for seed_index in range(1, 5):
            records.append(_seed(variant_index, seed_index))
    config = RunConfig(
        variant_filtering={
            "mode": "balanced",
            "max_variants_for_crest_per_molecule": 3,
            "max_seeds_per_variant": 2,
            "min_variants_to_keep": 1,
        }
    )

    selected, decisions = select_seed_records(records, config)

    assert len(selected) == 6
    assert len({record.parent_id for record in selected}) == 3
    selected_parent_ids = {record.parent_id for record in selected}
    assert all(
        len([record for record in selected if record.parent_id == parent]) == 2
        for parent in selected_parent_ids
    )
    assert any(not decision.selected for decision in decisions)


def test_exhaustive_mode_keeps_all_variants_and_marks_expensive() -> None:
    records = [_stereo(index) for index in range(1, 21)]
    config = RunConfig(variant_filtering={"mode": "exhaustive"})

    selected, decisions = select_stereo_records(records, config, "pre_3d")

    assert selected == records
    assert len(decisions) == len(records)
    assert all(decision.selected for decision in decisions)
    assert {decision.reason for decision in decisions} == {"exhaustive_mode_expensive"}


def test_rescue_selected_variants_include_rescue_reason() -> None:
    records = [_stereo(index) for index in range(1, 4)]
    config = RunConfig(
        variant_filtering={
            "mode": "balanced",
            "max_variants_before_3d_per_molecule": 1,
            "absolute_penalty_cutoff": 0.0,
            "relative_penalty_cutoff": 0.0,
            "keep_original_state": True,
            "keep_best_per_charge_state": False,
            "keep_best_per_formula": False,
            "keep_best_per_protomer": False,
            "keep_best_per_tautomer_family": False,
            "rescue_rules_enabled": True,
        }
    )

    selected, decisions = select_stereo_records(records, config, "pre_3d")

    assert len(selected) == 1
    selected_decision = next(decision for decision in decisions if decision.selected)
    assert selected_decision.rescue_reason == "rescue_original_input_state"
    assert selected_decision.reason == "pre_3d_rescue_rule"


def _stereo(index: int) -> StereoRecord:
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    return StereoRecord(
        id=f"mol_000001_p01_hash_t01_hash_s{index:02d}_hash",
        parent_id="mol_000001_p01_hash_t01_hash",
        input_molecule_id="mol_000001",
        molname="mol",
        canonical_smiles="CCO",
        isomeric_smiles="CCO",
        molecular_formula="C2H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        source_python_function="test",
        stereo_index=index,
        rdkit_mol=mol,
    )


def _seed(variant_index: int, seed_index: int) -> SeedConformerRecord:
    return SeedConformerRecord(
        id=f"seed_{variant_index:02d}_{seed_index:02d}",
        parent_id=f"stereo_{variant_index:02d}",
        input_molecule_id="mol_000001",
        molname="mol",
        canonical_smiles="CCO",
        isomeric_smiles="CCO",
        molecular_formula="C2H6O",
        formal_charge=0,
        explicit_proton_count=6,
        source_software="test",
        source_python_function="test",
        conformer_index=seed_index,
        energy_kcal_mol=float(variant_index * 10 + seed_index),
        embedding_status="success",
    )
