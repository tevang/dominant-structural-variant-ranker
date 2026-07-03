from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

from dsvr.models import StereoRecord


def molecular_complexity_penalty(record: StereoRecord) -> float:
    mol = record.rdkit_mol
    if mol is None:
        return 10.0
    heavy_atoms = mol.GetNumHeavyAtoms()
    rotors = Lipinski.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    hetero_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6})
    formal_charge = abs(Chem.GetFormalCharge(mol))
    exact_mw = Descriptors.ExactMolWt(mol)
    return (
        0.03 * heavy_atoms
        + 0.18 * rotors
        + 0.08 * hetero_atoms
        + 0.15 * rings
        + 0.5 * formal_charge
        + max(0.0, exact_mw - 550.0) / 250.0
    )
