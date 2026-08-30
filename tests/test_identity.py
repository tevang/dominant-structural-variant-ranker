"""Unit tests for the toolkit-local molecular identity policy.

The cases here pin the policy contract from the molecular-identity spec:
equivalent representations merge, while distinct tautomers, distinct
charge-placement/zwitterion states, and distinct stereoisomers stay distinct.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dsvr.chemistry import identity
from dsvr.chemistry.identity import exact_duplicate_key, exact_duplicate_key_with_warning


def _key(smiles: str) -> identity.ExactDuplicateKey:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return exact_duplicate_key(mol)


def test_same_structure_different_kekule_forms_merge() -> None:
    assert _key("c1cc[nH]c1") == _key("C1=CNC=C1")


def test_equivalent_charge_separated_drawings_merge() -> None:
    # Nitro group drawn charge-separated vs hypervalent: same structure.
    assert _key("c1ccc([N+](=O)[O-])cc1") == _key("c1ccc(N(=O)=O)cc1")


def test_distinct_tautomers_stay_distinct() -> None:
    # mol_000007 azole pair: 1H-indazole vs 2H-indazole (proton on different
    # ring nitrogen) are distinct tautomers and must never merge.
    one_h = _key("O=C(NC1CCCCC1)c1cccc(-c2n[nH]c3ccc(-c4nc[nH]n4)cc23)c1")
    two_h = _key("O=C(NC1CCCCC1)c1cccc(-c2nc3ccc(-c4nc[nH]n4)cc3[nH]2)c1")
    assert one_h.formula == two_h.formula
    assert one_h.net_charge == two_h.net_charge
    assert one_h != two_h


def test_charge_placement_variants_stay_distinct() -> None:
    # Zwitterionic vs neutral drawing of the same amino-acid connectivity:
    # same formula, different formal-charge placement.
    zwitterion = _key("C[C@H]([NH3+])C(=O)[O-]")
    neutral = _key("C[C@H](N)C(=O)O")
    assert zwitterion.formula == neutral.formula
    assert zwitterion != neutral
    # A charge-separated heterocycle vs its neutral drawing likewise.
    assert _key("[O-]C(=O)c1cccc[nH+]1") != _key("OC(=O)c1ccccn1")


def test_ionization_states_stay_distinct() -> None:
    assert _key("CC(=O)O") != _key("CC(=O)[O-]")


def test_stereoisomers_stay_distinct() -> None:
    assert _key("C[C@H](N)C(=O)O") != _key("C[C@@H](N)C(=O)O")


def test_key_contents_match_policy_tuple() -> None:
    mol = Chem.MolFromSmiles("CC(=O)O")
    assert mol is not None
    key = exact_duplicate_key(mol)
    assert isinstance(key, tuple)
    assert len(key) == 3
    assert key.formula == rdMolDescriptors.CalcMolFormula(mol)
    assert key.net_charge == 0
    assert key.canonical_isomeric_smiles == Chem.MolToSmiles(
        mol, canonical=True, isomericSmiles=True
    )


def test_cleanup_failure_falls_back_to_uncleaned_molecule(monkeypatch) -> None:
    def raising_cleanup(mol: Chem.Mol) -> Chem.Mol:
        raise RuntimeError("mock cleanup failure")

    monkeypatch.setattr(identity.rdMolStandardize, "Cleanup", raising_cleanup)
    mol = Chem.MolFromSmiles("CC(=O)O")
    assert mol is not None
    key, warning = exact_duplicate_key_with_warning(mol)
    assert warning is not None
    assert "un-cleaned" in warning
    # Fallback key is derived from the un-cleaned molecule.
    assert key == identity.ExactDuplicateKey(
        formula=rdMolDescriptors.CalcMolFormula(mol),
        net_charge=Chem.GetFormalCharge(mol),
        canonical_isomeric_smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
    )


def test_key_derivation_does_not_mutate_input() -> None:
    mol = Chem.MolFromSmiles("c1ccc(N(=O)=O)cc1")
    assert mol is not None
    before = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    exact_duplicate_key(mol)
    after = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    assert before == after
