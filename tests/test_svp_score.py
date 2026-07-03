from __future__ import annotations

import pytest
from rdkit import Chem

from dsvr.filtering.variant_score import (
    RT_LN10_298_KCAL_MOL,
    score_chemistry_sanity,
    score_tautomer,
    score_variant,
)
from dsvr.models import StereoRecord


def test_ph_penalty_uses_mocked_pka_without_claiming_population() -> None:
    record = _stereo("CCN", "amine")

    breakdown = score_variant(record, {"ph": 7.0, "pka_info": {"pka": 9.0}})

    assert breakdown.protomer_penalty == pytest.approx(2.0 * RT_LN10_298_KCAL_MOL)
    assert any("approximate_pH_mismatch" in reason for reason in breakdown.reasons)
    assert any("not final thermodynamics" in warning for warning in breakdown.warnings)


def test_tautomer_score_penalizes_aromaticity_loss_feature() -> None:
    record = _stereo(
        "C1=CC=CC=C1",
        "aromatic_loss",
        metadata={"tautomer": {"aromatic_atom_loss": 2}},
    )

    penalty, reasons = score_tautomer(record)

    assert penalty >= 3.0
    assert any("aromaticity_loss_heuristic" in reason for reason in reasons)
    assert any("not a stability ranking" in reason for reason in reasons)


def test_charge_separation_is_penalized_but_not_silently_rejected() -> None:
    record = _stereo("[NH3+]CC[O-]", "zwitterion")

    penalty, reasons = score_chemistry_sanity(record)

    assert 0.0 < penalty < 1000.0
    assert any("separated_formal_charges" in reason for reason in reasons)


def _stereo(
    smiles: str,
    suffix: str,
    *,
    stereo_index: int = 1,
    metadata: dict[str, object] | None = None,
) -> StereoRecord:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return StereoRecord(
        id=f"mol_000001_p01_hash_t01_hash_s{stereo_index:02d}_{suffix}",
        parent_id="mol_000001_p01_hash_t01_hash",
        input_molecule_id="mol_000001",
        molname="mol",
        canonical_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False),
        isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        molecular_formula="C2H6O",
        formal_charge=Chem.GetFormalCharge(mol),
        explicit_proton_count=0,
        source_software="test",
        source_python_function="test",
        stereo_index=stereo_index,
        rdkit_mol=mol,
        metadata=metadata or {},
    )
