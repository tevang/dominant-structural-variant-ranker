from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from dsvr.models import StereoRecord

RT_LN10_298_KCAL_MOL = 1.364
REJECTION_PENALTY = 1_000.0
_METALS = {
    3,
    4,
    11,
    12,
    13,
    19,
    20,
    21,
    22,
    23,
    24,
    25,
    26,
    27,
    28,
    29,
    30,
    31,
    37,
    38,
    39,
    40,
    41,
    42,
    43,
    44,
    45,
    46,
    47,
    48,
    49,
    50,
    55,
    56,
    57,
    72,
    73,
    74,
    75,
    76,
    77,
    78,
    79,
    80,
    81,
    82,
    83,
}


@dataclass(frozen=True)
class PenaltyBreakdown:
    """Transparent SVPScore components.

    SVPScore is a triage score for filtering and scheduling expensive
    calculations. It is not a thermodynamic free energy and must not be
    interpreted as a rigorous state population model.
    """

    protomer_penalty: float
    tautomer_penalty: float
    stereo_penalty: float
    chemistry_sanity_penalty: float
    complexity_penalty: float
    cheap_3d_energy_penalty: float
    total: float
    reasons: list[str]
    warnings: list[str]


@dataclass(frozen=True)
class VariantScore:
    record_id: str
    penalty: float
    reason: str


def score_protomer(
    record: Any,
    ph: float,
    pka_info: Mapping[str, Any] | None = None,
) -> tuple[float, list[str]]:
    penalty = 0.0
    reasons = [
        "SVPScore protomer component is a candidate-filtering heuristic, not a pH population model."
    ]
    pka = _first_number(pka_info, ("micro_pka", "site_pka", "pka")) if pka_info else None
    if pka is not None:
        mismatch = abs(ph - pka)
        penalty += RT_LN10_298_KCAL_MOL * mismatch
        reasons.append(
            "approximate_pH_mismatch_penalty_using_RTln10_scale:"
            f"pH={ph:.2f},pKa={pka:.2f},delta={mismatch:.2f}"
        )
    else:
        reasons.append("no_micro_pKa_available_no_pH_population_penalty_invented")

    charge = abs(_formal_charge(record))
    if charge > 1:
        penalty += 2.0 * (charge - 1)
        reasons.append(f"high_formal_charge_magnitude:+{2.0 * (charge - 1):.2f}")

    mol = _mol(record)
    if mol is not None:
        positive, negative = _formal_charge_counts(mol)
        if positive and negative:
            charge_separation = min(6.0, 2.0 + positive + negative)
            if _solvent_from_record(record).lower() not in {"water", "methanol", "ethanol", "dmso"}:
                charge_separation += 2.0
                reasons.append("zwitterion_or_separated_charge_penalized_more_in_nonpolar_solvent")
            penalty += charge_separation
            reasons.append(f"formal_charge_separation:+{charge_separation:.2f}")
    return penalty, reasons


def score_tautomer(
    record: Any,
    parent_best_score: float | None = None,
) -> tuple[float, list[str]]:
    penalty = 0.0
    reasons = [
        "RDKit tautomer features are used only for triage; "
        "canonical tautomer is not a stability ranking."
    ]
    mol = _mol(record)
    metadata = _metadata(record)
    tautomer_metadata = _mapping(metadata.get("tautomer"))
    aromatic_loss = _first_number(tautomer_metadata, ("aromatic_atom_loss", "aromatic_loss"))
    if aromatic_loss is None:
        parent_aromatic = _first_number(metadata, ("parent_aromatic_atom_count",))
        if parent_aromatic is not None and mol is not None:
            aromatic_loss = max(0.0, parent_aromatic - _aromatic_atom_count(mol))
    if aromatic_loss and aromatic_loss > 0:
        value = 1.5 * aromatic_loss
        penalty += value
        reasons.append(f"aromaticity_loss_heuristic:+{value:.2f}")

    transforms = _first_number(tautomer_metadata, ("transform_count", "num_transforms"))
    if transforms and transforms > 2:
        value = 0.5 * (transforms - 2)
        penalty += value
        reasons.append(f"many_tautomer_transforms_from_parent:+{value:.2f}")

    if mol is not None:
        try:
            enumerator = rdMolStandardize.TautomerEnumerator()
            canonical = enumerator.Canonicalize(mol)
            canonical_smiles = Chem.MolToSmiles(canonical, canonical=True, isomericSmiles=True)
            if canonical_smiles != getattr(record, "isomeric_smiles", None):
                penalty += 0.25
                reasons.append("differs_from_rdkit_canonical_tautomer_small_feature_penalty:+0.25")
        except (RuntimeError, ValueError):
            penalty += 2.0
            reasons.append("rdkit_tautomer_canonicalization_failed:+2.00")
        positive, negative = _formal_charge_counts(mol)
        if positive and negative:
            value = min(4.0, 1.0 + positive + negative)
            penalty += value
            reasons.append(f"tautomer_formal_charge_separation:+{value:.2f}")

    if parent_best_score is not None and parent_best_score > 0:
        value = min(3.0, 0.2 * parent_best_score)
        penalty += value
        reasons.append(f"parent_tautomer_priority_offset:+{value:.2f}")
    return penalty, reasons


def score_stereoisomer(
    record: Any,
    achiral_environment: bool = True,
) -> tuple[float, list[str]]:
    penalty = 0.0
    reasons = ["stereo component represents enumeration uncertainty, not intrinsic high energy."]
    stereo_index = int(getattr(record, "stereo_index", 0) or 0)
    if stereo_index > 1:
        value = 0.05 * (stereo_index - 1)
        penalty += value
        reasons.append(f"additional_enumerated_stereoisomer_mild_priority_penalty:+{value:.2f}")
    if achiral_environment and _single_center_chiral(record):
        reasons.append("single-center_enantiomer_candidate_can_be_collapsed_in_achiral_environment")
    metadata = _metadata(record)
    if metadata.get("input_stereo_specified") is True:
        reasons.append("input_stereochemistry_specified_preserved_by_default")
    return penalty, reasons


def score_chemistry_sanity(record: Any) -> tuple[float, list[str]]:
    penalty = 0.0
    reasons: list[str] = []
    mol = _mol(record)
    if mol is None:
        return REJECTION_PENALTY, ["reject_missing_rdkit_mol"]
    try:
        probe = Chem.Mol(mol)
        Chem.SanitizeMol(probe)
    except (RuntimeError, ValueError) as exc:
        return REJECTION_PENALTY, [f"reject_rdkit_sanitization_failure:{exc}"]

    radicals = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if radicals:
        return REJECTION_PENALTY, [f"reject_radicals_not_allowed:{radicals}"]

    positive, negative = _formal_charge_counts(mol)
    charged_atoms = positive + negative
    if charged_atoms >= 3:
        value = 0.75 * charged_atoms
        penalty += value
        reasons.append(f"many_formally_charged_atoms:+{value:.2f}")
    if positive and negative:
        value = min(4.0, 1.0 + charged_atoms)
        penalty += value
        reasons.append(f"separated_formal_charges:+{value:.2f}")

    metals = [atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() in _METALS]
    if metals:
        value = 8.0
        penalty += value
        reasons.append(f"metal_or_unsupported_atom_route_to_warning:{','.join(sorted(set(metals)))}")
    return penalty, reasons or ["chemistry_sanity_checks_passed"]


def score_complexity(record: Any) -> tuple[float, list[str]]:
    mol = _mol(record)
    if mol is None:
        return 10.0, ["missing_mol_complexity_fallback:+10.00"]
    heavy_atoms = mol.GetNumHeavyAtoms()
    rotors = Lipinski.NumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    hetero_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6})
    exact_mw = Descriptors.ExactMolWt(mol)
    stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    penalty = (
        0.03 * heavy_atoms
        + 0.18 * rotors
        + 0.08 * hetero_atoms
        + 0.15 * rings
        + 0.10 * stereo_centers
        + max(0.0, exact_mw - 550.0) / 250.0
    )
    reasons = [
        "complexity_penalty_affects_scheduling_priority_not_physical_stability:"
        f"heavy_atoms={heavy_atoms},rotors={rotors},rings={rings},"
        f"hetero_atoms={hetero_atoms},stereo_centers={stereo_centers}"
    ]
    return penalty, reasons


def score_cheap_3d_energy(
    record: Any,
    best_energy: float | None = None,
) -> tuple[float, list[str]]:
    energy = getattr(record, "energy_kcal_mol", None)
    if energy is None:
        return 0.0, ["no_cheap_3d_energy_available"]
    if best_energy is None:
        return 0.0, ["cheap_3d_energy_available_but_no_group_reference"]
    relative = max(0.0, float(energy) - best_energy)
    return relative, [f"relative_cheap_3d_energy_penalty:+{relative:.2f}"]


def score_variant(record: Any, context: Mapping[str, Any]) -> PenaltyBreakdown:
    ph = float(context.get("ph", 7.0))
    protomer, protomer_reasons = score_protomer(
        record,
        ph,
        _mapping_or_none(context.get("pka_info")),
    )
    tautomer, tautomer_reasons = score_tautomer(
        record,
        _optional_float(context.get("parent_best_score")),
    )
    stereo, stereo_reasons = score_stereoisomer(
        record,
        bool(context.get("achiral_environment", True)),
    )
    sanity, sanity_reasons = score_chemistry_sanity(record)
    complexity, complexity_reasons = score_complexity(record)
    cheap_3d, cheap_3d_reasons = score_cheap_3d_energy(
        record,
        _optional_float(context.get("best_energy")),
    )
    total = protomer + tautomer + stereo + sanity + complexity + cheap_3d
    warnings = [
        "SVPScore is a transparent filtering heuristic, not final thermodynamics.",
        "pH/protomer penalties are approximate unless real micro-pKa corrections are provided.",
    ]
    if sanity >= REJECTION_PENALTY:
        warnings.append("record should be rejected by chemistry sanity checks")
    return PenaltyBreakdown(
        protomer_penalty=protomer,
        tautomer_penalty=tautomer,
        stereo_penalty=stereo,
        chemistry_sanity_penalty=sanity,
        complexity_penalty=complexity,
        cheap_3d_energy_penalty=cheap_3d,
        total=total,
        reasons=[
            *protomer_reasons,
            *tautomer_reasons,
            *stereo_reasons,
            *sanity_reasons,
            *complexity_reasons,
            *cheap_3d_reasons,
        ],
        warnings=warnings,
    )


def cheap_variant_score(record: StereoRecord) -> VariantScore:
    breakdown = score_variant(record, {"ph": 7.0})
    return VariantScore(
        record_id=record.id,
        penalty=breakdown.total,
        reason="SVPScore_total_compatibility_wrapper",
    )


def enantiomer_collapse_key(record: Any) -> str | None:
    """Return a conservative key only for single-center enantiomer pairs.

    A non-isomeric SMILES key is safe for single stereocenter pairs in achiral
    environments. Molecules with multiple stereocenters are not collapsed here
    because they may include diastereomers.
    """

    mol = _mol(record)
    if mol is None or not _single_center_chiral(record):
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


def _mol(record: Any) -> Chem.Mol | None:
    mol = getattr(record, "rdkit_mol", None)
    if mol is not None:
        return mol
    smiles = getattr(record, "isomeric_smiles", None) or getattr(record, "canonical_smiles", None)
    if smiles:
        return Chem.MolFromSmiles(smiles)
    return None


def _metadata(record: Any) -> dict[str, Any]:
    value = getattr(record, "metadata", None)
    return dict(value) if isinstance(value, Mapping) else {}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _mapping_or_none(value: Any) -> Mapping[str, Any] | None:
    return value if isinstance(value, Mapping) else None


def _first_number(mapping: Mapping[str, Any] | None, keys: tuple[str, ...]) -> float | None:
    if not mapping:
        return None
    for key in keys:
        value = mapping.get(key)
        if isinstance(value, int | float):
            return float(value)
    return None


def _optional_float(value: Any) -> float | None:
    if isinstance(value, int | float):
        return float(value)
    return None


def _formal_charge(record: Any) -> int:
    value = getattr(record, "formal_charge", None)
    if isinstance(value, int):
        return value
    mol = _mol(record)
    return Chem.GetFormalCharge(mol) if mol is not None else 0


def _formal_charge_counts(mol: Chem.Mol) -> tuple[int, int]:
    positive = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    negative = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    return positive, negative


def _solvent_from_record(record: Any) -> str:
    metadata = _metadata(record)
    solvent = metadata.get("solvent")
    return str(solvent) if solvent is not None else "water"


def _aromatic_atom_count(mol: Chem.Mol) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())


def _single_center_chiral(record: Any) -> bool:
    mol = _mol(record)
    if mol is None:
        return False
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    return len(centers) == 1
