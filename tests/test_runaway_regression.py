from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from dsvr.config import load_config
from dsvr.filtering import selection
from dsvr.filtering.selection import select_stereo_records
from dsvr.filtering.variant_score import PenaltyBreakdown
from dsvr.models import StereoRecord


def test_runaway_variant_count_is_bounded_by_production_balanced(monkeypatch) -> None:
    # 20 protomers * 100 tautomers * 64 stereoisomers = 128,000 structural variants.
    records = [_stereo(index) for index in range(20 * 100 * 64)]
    config = load_config(Path("configs/production_balanced.yaml"))

    def fake_score(record, context):
        index = int(record.metadata["synthetic_index"])
        return PenaltyBreakdown(
            protomer_penalty=0.0,
            tautomer_penalty=0.0,
            stereo_penalty=0.0,
            chemistry_sanity_penalty=0.0,
            complexity_penalty=float(index % 1000) / 1000.0,
            cheap_3d_energy_penalty=0.0,
            total=float(index % 1000) / 1000.0,
            reasons=["synthetic_runaway_score"],
            warnings=[],
        )

    monkeypatch.setattr(selection, "score_variant", fake_score)

    selected_pre3d, decisions_pre3d = select_stereo_records(records, config, "pre_3d")
    selected_cheap, decisions_cheap = select_stereo_records(selected_pre3d, config, "cheap_score")

    assert len(records) == 128_000
    assert (
        len(selected_cheap)
        <= config.variant_filtering.max_variants_after_cheap_score_per_molecule
    )
    assert len([decision for decision in decisions_pre3d if not decision.selected]) > 0
    assert len([decision for decision in decisions_cheap if not decision.selected]) > 0


def _stereo(index: int) -> StereoRecord:
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    return StereoRecord.model_construct(
        id=f"mol_p{index // 6400:02d}_t{index // 64:04d}_s{index % 64:02d}",
        parent_id=f"mol_p{index // 6400:02d}_t{index // 64:04d}",
        input_molecule_id="mol",
        molname="mol",
        canonical_smiles="CCO",
        isomeric_smiles="CCO",
        molecular_formula="C2H6O",
        formal_charge=0,
        explicit_proton_count=6,
        stage_name="stereo",
        source_software="test",
        source_python_function="test",
        warnings=[],
        metadata={"synthetic_index": index},
        output_paths=[],
        stereo_index=index,
        rdkit_mol=mol,
    )
