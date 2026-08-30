## Context

See proposal.md for root cause and evidence. Relevant mechanics: the engine collects `TautomerRecord`s from two paths — `filter_tautomers_with_auto3d` (ligprep_like default) and `enumerate_tautomers` (fallback/non-default) — and resume may also load them from disk (`_load_tautomers`, `_load_existing_tautomer_outputs`). All dedupe today is branch-local: protomer dedupe per input, `_dedupe_candidates` per protomer, stereo per tautomer, `final3d._dedupe_one_conformer_per_variant` per parent. Ranked-but-rejected tautomer candidates exist only transiently inside the filter (summarized to CSV), so refill requires exposing that pool.

## Goals / Non-Goals

**Goals:**
- One shared identity module used by every dedupe point, so the policy is documented and testable in a single place.
- Cross-branch dedupe implemented once in the engine, covering generated and resume-loaded records uniformly.
- Refill to cap without re-running external tools (molscrub/Auto3D) by promoting already-ranked unused candidates.

**Non-Goals:**
- Tautomer-invariant protomer merging at the protonation stage.
- Changing `tauto_k` semantics from per-protomer-branch to per-input-molecule.
- New external dependencies; RDKit `rdMolStandardize` already ships with the project.

## Decisions

1. **Identity keys live in a new `dsvr/chemistry/identity.py`.** `exact_duplicate_key(mol) -> (formula, net_charge, canonical_isomeric_smiles)` computed from `rdMolStandardize.Cleanup(mol)` (ChEMBL-style standardization; never `Uncharger`, `ChargeParent`, or `TautomerParent`, which would over-merge charge/tautomer states). Key derivation never mutates or replaces the stored structure. If `Cleanup` fails, fall back to the un-cleaned molecule and record a warning. Deliberately excluded: Standard InChI/InChIKey (mobile-H tautomer collapsing) and raw SMILES (no policy). Rationale: matches the layered-identity practice in the referenced standardization guidance; one module = one auditable policy. Alternative considered: InChI with FixedH — rejected, adds an interop identifier without addressing the policy problem and duplicates what standardized isomeric SMILES already expresses toolkit-locally.

2. **Cross-branch tautomer dedupe runs in the engine after tautomer collection** (single call site in `engine.py` after the collection/resume-load block), grouping by `(input_molecule_id, exact_duplicate_key)`. Representative selection: lowest `metadata.auto3d_tautomer_filtering.relative_energy_kcal_mol` when present, else RDKit-fallback rank, else deterministic smallest record ID. Merged records contribute their IDs to `metadata["merged_from"]` on the representative; every elimination is appended to `enumeration/tautomers/tautomer_dedupe.csv`. Rationale over embedding the dedupe inside each tautomer path: one implementation covers both generation paths and resume-loaded records, and keeps provenance rewriting in one place.

3. **Both tautomer paths expose a ranked pool for refill.** `filter_tautomers_with_auto3d` returns a small result object (selected records + ranked candidate payloads for rejected candidates, sufficient to materialize `TautomerRecord`s without re-running Auto3D); `enumerate_tautomers` likewise exposes its scored unique pool beyond the selected subset. The engine refills each branch that dropped below `tauto_k` from that branch's own ranked pool, skipping candidates whose exact key is already selected anywhere for that input molecule; shortfalls are recorded in the dedupe CSV and `warnings.jsonl`. Engine call sites adapt to the result object; per-protomer SDF/CSV outputs continue to be written by the paths themselves. Alternative considered (engine re-reads `tautomers_auto3d_ranked.csv`): rejected — couples the engine to a CSV format when the data is already in memory.

4. **Stereoisomer refill by bounded over-enumeration.** `enumerate_stereoisomers` enumerates with `maxIsomers = min(cap * ceiling_multiplier, hard_ceiling)`, dedupes, then keeps the first `cap` unique. New config `stereoisomer_filtering.enumeration_ceiling_multiplier: int = 4` (additive, default preserves fast behavior while enabling refill; `stereo_unique=True` unchanged). Shortfall beyond the ceiling is recorded. Rationale: `maxIsomers` is RDKit's only control; refill "if possible" within that enumeration space.

5. **Cross-tautomer stereo guard lives in `stereo_auto3d_filter`.** It already sees every `StereoRecord` per input molecule before Auto3D ranking, so the exact-key guard dedupes there (audited to the stereoisomer_filtering directory) and each tautomer's selection is refilled from its own ranked unused stereoisomers. Rationale: avoids a second global pass in the engine and reuses existing energy-ranked decisions.

6. **Final safety net in `generate_final_3d_variants`.** After all records are assembled (including RDKit fallbacks), dedupe per input molecule by exact key, keep lowest `energy_kcal_mol` (None sorts last; tie-break by record ID), record merges in metadata and `final_dedupe_audit.csv`. Config `final_3d.dedupe_final_variants: bool = true` (additive). This directly guarantees the behavior in the "Clean final variant files" spec scenario even for resumed/mixed runs.

## Risks / Trade-offs

- `rdMolStandardize.Cleanup` includes reionization, which can relocate protons in zwitterions → keys could over-merge distinct charge-placement variants at the tautomer stage. Mitigation: identity unit tests include zwitterion/charge-placement and azole-tautomer cases; if observed, restrict the identity cleanup to functional-group normalization without reionization (documented in the module).
- Refill promotes candidates that today count as rejected → `tautomers_selected.csv`/stereo selection files gain rows relative to current runs. Mitigation: audit CSVs and `variant_counts.csv`-level summaries keep the deltas explainable.
- Resume behavior: input sets of downstream stages change, so resumed runs re-execute tautomer-and-later stages. Mitigation: existing input-hash resume logic handles this naturally; documented in docs.
- Return-type change of `filter_tautomers_with_auto3d` / `enumerate_tautomers` touches the engine and tests. Mitigation: small result dataclasses; direct callers are few (engine + tests + retry closures), updated mechanically.
- Stereo over-enumeration increases RDKit enumeration work by up to the ceiling multiplier per tautomer. Mitigation: hard ceiling cap and config knob; net runtime still improves because duplicate Auto3D work downstream disappears.

## Migration Plan

Apply code, then rerun affected runs (or let `--resume` re-execute tautomer-and-later stages automatically). No config migration; both new knobs are additive with defaults enabling the fix.

## Open Questions

None blocking.
