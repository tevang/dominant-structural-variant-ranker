"""Unspecified-stereochemistry policy for Auto3D stages (fix-auto3d-integration §4).

Auto3D warns when molecules with unspecified stereochemistry are submitted
while isomer enumeration is disabled. Before every Auto3D stage call the
pipeline therefore applies one explicit, configured policy
(``config.auto3d.on_unspecified_stereo``):

- ``enumerate`` (default): enumerate unspecified stereochemistry with the
  existing RDKit stereo enumeration up-front; Auto3D only ever receives
  fully specified structures. Energies of enumerated isomers are aggregated
  back to their originating variant (minimum energy).
- ``auto3d_enumerate``: keep the molecules as-is but invoke Auto3D with
  isomer enumeration enabled (``--enumerate-isomer``) for exactly the
  affected sub-batches.

Both policies record their application in variant provenance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from rdkit import Chem

from dsvr.chemistry.stereochemistry import _enumerate_with_timeout
from dsvr.config import RunConfig

ISOMER_SUFFIX_SEPARATOR = "__st"

POLICY_TREATMENT_ENUMERATE = "enumerated_upfront"
POLICY_TREATMENT_AUTO3D_ENUMERATE = "auto3d_enumerate_isomer"

_TREATED_NOTE = (
    "unspecified stereochemistry treated per policy "
    "auto3d.on_unspecified_stereo={policy} ({treatment})"
)


def has_unspecified_stereo(mol: Chem.Mol) -> bool:
    """Return True when the molecule carries unspecified stereo elements."""

    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    if any(label == "?" for _index, label in centers):
        return True
    try:
        bonds = Chem.FindPotentialStereoBonds(mol) or []
    except (AttributeError, RuntimeError):  # pragma: no cover - old RDKit fallback
        return False
    return any(
        getattr(bond, "specified", Chem.StereoSpecified.Specified) == Chem.StereoSpecified.Unspecified
        for bond in bonds
    )


def enumerate_unspecified_isomers(mol: Chem.Mol, config: RunConfig) -> list[Chem.Mol]:
    """Enumerate stereoisomers of an unspecified-stereo molecule.

    Uses the pipeline's existing bounded stereo enumeration; falls back to
    the input molecule when enumeration fails.
    """

    try:
        return _enumerate_with_timeout(
            mol,
            timeout_seconds=config.stereoisomer_filtering.timeout_seconds_per_tautomer,
            try_embedding=False,
            only_unassigned=True,
            unique=config.enumeration.stereo_unique,
            max_isomers=config.stereoisomer_filtering.max_stereoisomers_per_tautomer,
            random_seed=config.enumeration.stereo_random_seed,
        )
    except (TimeoutError, RuntimeError):
        return [Chem.Mol(mol)]


@dataclass(frozen=True)
class StereoPolicyPlan:
    """How an Auto3D batch must be prepared under the configured policy.

    ``expanded`` maps base ids (tautomer / stereo record ids) to the mols
    Auto3D may receive; when a base id expanded to several isomers the
    downstream Auto3D lines carry ``<base_id>__st<n>`` names and energies
    are aggregated back by minimum. ``enumerate_isomer_for`` lists base ids
    that must be invoked with Auto3D isomer enumeration enabled (only used
    by the ``auto3d_enumerate`` policy).
    """

    expanded: list[tuple[str, Chem.Mol, str]]
    unspecified_ids: tuple[str, ...]
    enumerate_isomer_for: frozenset[str]
    treatment: str


def apply_stereo_policy(
    items: list[tuple[str, Chem.Mol]],
    config: RunConfig,
) -> StereoPolicyPlan:
    """Apply ``config.auto3d.on_unspecified_stereo`` to ``(base_id, mol)`` pairs."""

    policy = config.auto3d.on_unspecified_stereo
    unspecified = tuple(
        base_id for base_id, mol in items if has_unspecified_stereo(mol)
    )
    if not unspecified:
        return StereoPolicyPlan(
            expanded=[(base_id, mol, base_id) for base_id, mol in items],
            unspecified_ids=(),
            enumerate_isomer_for=frozenset(),
            treatment="none_needed",
        )
    if policy == "enumerate":
        expanded: list[tuple[str, Chem.Mol, str]] = []
        enumerate_residue: set[str] = set()
        for base_id, mol in items:
            if base_id not in unspecified:
                expanded.append((base_id, mol, base_id))
                continue
            isomers = enumerate_unspecified_isomers(mol, config)
            if len(isomers) <= 1:
                if has_unspecified_stereo(isomers[0]):
                    # RDKit enumeration failed to resolve this molecule — do
                    # not silently send it with enumeration disabled (that is
                    # exactly what the policy promises to prevent). Route it
                    # to an isomer-enumerating Auto3D invocation instead and
                    # record the residue in provenance.
                    enumerate_residue.add(base_id)
                expanded.append((base_id, isomers[0], base_id))
                continue
            for index, isomer in enumerate(isomers, start=1):
                expanded.append(
                    (base_id, isomer, f"{base_id}{ISOMER_SUFFIX_SEPARATOR}{index}")
                )
        return StereoPolicyPlan(
            expanded=expanded,
            unspecified_ids=unspecified,
            enumerate_isomer_for=frozenset(enumerate_residue),
            treatment=POLICY_TREATMENT_ENUMERATE,
        )
    return StereoPolicyPlan(
        expanded=[(base_id, mol, base_id) for base_id, mol in items],
        unspecified_ids=unspecified,
        enumerate_isomer_for=frozenset(unspecified),
        treatment=POLICY_TREATMENT_AUTO3D_ENUMERATE,
    )


def treatment_note(config: RunConfig, treatment: str) -> str:
    """Short provenance note recorded on treated variants."""

    return _TREATED_NOTE.format(
        policy=config.auto3d.on_unspecified_stereo,
        treatment=treatment,
    )


def policy_metadata(config: RunConfig, plan: StereoPolicyPlan) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "policy": config.auto3d.on_unspecified_stereo,
        "treatment": plan.treatment,
        "unspecified_count": len(plan.unspecified_ids),
        "unspecified_ids": list(plan.unspecified_ids),
    }
    if plan.treatment == POLICY_TREATMENT_ENUMERATE and plan.enumerate_isomer_for:
        metadata["auto3d_enumerate_residue"] = sorted(plan.enumerate_isomer_for)
    return metadata
