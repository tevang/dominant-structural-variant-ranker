"""Toolkit-local molecular identity policy for exact-duplicate detection.

This module is the single auditable definition of when the workflow treats two
structures as the *same exact structure*. Every deduplication point in the
pipeline (tautomer selection, stereoisomer selection, final 3D output) derives
keys through :func:`exact_duplicate_key` so the policy stays testable in one
place.

Policy
------
The exact-duplicate key is ``(formula, net_charge, canonical_isomeric_smiles)``
computed from the RDKit-standardized molecule via
:func:`rdkit.Chem.MolStandardize.rdMolStandardize.Cleanup` (ChEMBL-style
functional-group normalization and reionization). Standardization:

* MUST NOT neutralize charges (no ``Uncharger``), strip fragments, remove
  stereochemistry, or collapse tautomers — in particular ``ChargeParent`` and
  ``TautomerParent`` are never applied, so distinct protonation states,
  distinct tautomers, and distinct stereoisomers keep distinct keys.
* merges equivalent drawings of one structure: different Kekulé forms and
  equivalent charge-separated representations (e.g. nitro group drawings)
  of the same connectivity yield the same key.

Deliberately forbidden as dedupe keys:

* **Standard InChI / Standard InChIKey** — the standardization layer of InChI
  deliberately collapses mobile-hydrogen tautomers into one identifier, which
  would over-merge distinct tautomers at exactly the stages that must keep
  them apart.
* **Raw, unstandardized SMILES strings** — no standardization policy, so
  Kekulé or charge-representation differences would under-merge.
* **IUPAC names** — not a structural identity.

Keys are **toolkit-local**: they are valid within this RDKit-based pipeline
under the policy stated here. They are not global identifiers and must not be
treated as interoperable with other toolkits' canonicalization.

If ``Cleanup`` fails for a molecule, the un-cleaned molecule is used as a
fallback and a warning is recorded (returned to the caller for surfacing in
run artifacts); key derivation never mutates the input molecule.
"""

from __future__ import annotations

import logging
from typing import NamedTuple

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger(__name__)


class ExactDuplicateKey(NamedTuple):
    """Identity of one exact structure under the toolkit-local policy."""

    formula: str
    net_charge: int
    canonical_isomeric_smiles: str


def exact_duplicate_key(mol: Chem.Mol) -> ExactDuplicateKey:
    """Return the exact-duplicate key for ``mol`` under the documented policy.

    On standardization failure the un-cleaned molecule is used and a warning
    is logged; use :func:`exact_duplicate_key_with_warning` when the caller
    needs to surface the warning in run artifacts.
    """
    key, warning = exact_duplicate_key_with_warning(mol)
    if warning:
        logger.warning("%s", warning)
    return key


def exact_duplicate_key_with_warning(mol: Chem.Mol) -> tuple[ExactDuplicateKey, str | None]:
    """Return ``(key, warning)``; ``warning`` is set when standardization failed."""

    standardized, warning = _standardize_for_identity_with_warning(mol)
    return (
        ExactDuplicateKey(
            formula=rdMolDescriptors.CalcMolFormula(standardized),
            net_charge=Chem.GetFormalCharge(standardized),
            canonical_isomeric_smiles=Chem.MolToSmiles(
                standardized, canonical=True, isomericSmiles=True
            ),
        ),
        warning,
    )


def _standardize_for_identity(mol: Chem.Mol) -> Chem.Mol:
    """Standardized copy used for key derivation (never mutates ``mol``)."""

    return _standardize_for_identity_with_warning(mol)[0]


def _standardize_for_identity_with_warning(mol: Chem.Mol) -> tuple[Chem.Mol, str | None]:
    """ChEMBL-style ``Cleanup`` with un-cleaned fallback.

    ``Cleanup`` performs functional-group normalization and reionization but
    never uncharges and never returns a tautomer/charge parent; it is applied
    to a copy so the stored structure is never mutated or replaced.
    """

    candidate = Chem.Mol(mol)
    try:
        cleaned = rdMolStandardize.Cleanup(candidate)
    except (RuntimeError, ValueError) as exc:
        warning = (
            "rdMolStandardize.Cleanup failed during exact-duplicate key derivation; "
            f"used the un-cleaned structure instead: {exc}"
        )
        return candidate, warning
    if cleaned is None or cleaned.GetNumAtoms() == 0:
        warning = (
            "rdMolStandardize.Cleanup produced no atoms during exact-duplicate key "
            "derivation; used the un-cleaned structure instead"
        )
        return candidate, warning
    return cleaned, None
