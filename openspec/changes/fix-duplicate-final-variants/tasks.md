## 1. Molecular identity module

- [x] 1.1 Create `src/dsvr/chemistry/identity.py` with `exact_duplicate_key(mol)` returning `(formula, net_charge, canonical_isomeric_smiles)` computed from `rdMolStandardize.Cleanup(mol)`, with un-cleaned fallback + warning on cleanup failure; document the policy (no uncharging, no TautomerParent/ChargeParent, no Standard InChI, toolkit-local). Verify: `ruff check src tests` and `mypy src` pass on the new module
- [x] 1.2 Add unit tests `tests/test_identity.py`: same-structure-different-Kekulé merges; distinct tautomers (mol_000007 azole pairs) stay distinct; charge-placement/zwitterion variants stay distinct; stereoisomers stay distinct. Verify: `pytest tests/test_identity.py` passes; if reionization over-merges a charge-placement case, restrict cleanup to functional-group normalization (per design.md risk) and re-verify

## 2. Cross-branch tautomer dedupe and refill

- [x] 2.1 Change `filter_tautomers_with_auto3d` to return a result object carrying selected records plus a ranked pool of rejected candidates (sufficient to materialize `TautomerRecord`s without re-running Auto3D); update engine call site and retry/fallback closures. Verify: existing tautomer-filter tests pass with mechanical updates
- [x] 2.2 Expose the scored unique pool beyond the selected subset from `enumerate_tautomers` (plain RDKit path). Verify: existing `tautomers.py` tests pass with mechanical updates
- [x] 2.3 Implement engine-level `dedupe_and_refill_tautomers` after tautomer collection/resume-load: group per `(input_molecule_id, exact_duplicate_key)`, keep best-ranked representative (lowest tautomer relative energy, else fallback rank, else smallest record ID), record `merged_from` provenance, and refill each branch below `tauto_k` from its own ranked unused unique pool. Verify: unit test on synthetic cross-branch duplicates shows one record per unique structure and branch refilled to `tauto_k`
- [x] 2.4 Write `enumeration/tautomers/tautomer_dedupe.csv` audit rows (retained ID, eliminated IDs, refill promotions, shortfall reason) and surface shortfalls in `warnings.jsonl`. Verify: audit CSV contents in the unit test from 2.3

## 3. Stereoisomer dedupe and refill

- [x] 3.1 Add config `stereoisomer_filtering.enumeration_ceiling_multiplier` (default 4, validated > 0) and apply bounded over-enumeration in `enumerate_stereoisomers` (`maxIsomers = min(cap * multiplier, hard ceiling)`), keeping the first `cap` unique; record shortfalls. Verify: unit test with a molecule producing duplicate-rich enumeration reaches `max_stereoisomers_per_tautomer` unique isomers when available
- [x] 3.2 Add exact-key cross-tautomer guard in `filter_stereoisomers_with_auto3d` before Auto3D ranking, refilling each tautomer's selection from its own ranked unused stereoisomers; audit to the stereoisomer_filtering output directory. Verify: unit test with duplicated stereoisomers across tautomers shows one copy selected and refill from the unused pool

## 4. Final safety-net dedupe

- [x] 4.1 Add config `final_3d.dedupe_final_variants` (default true) and dedupe final records per input molecule by exact key in `generate_final_3d_variants` before `_write_final_outputs`, keeping lowest `energy_kcal_mol` (None last; tie-break by record ID) with `merged_from` metadata. Verify: unit test feeds duplicate-geometry records and asserts one output per unique structure
- [x] 4.2 Write `final_dedupe_audit.csv` alongside final outputs. Verify: audit file exists and lists retained/eliminated IDs in the 4.1 test

## 5. Regression coverage and docs

- [x] 5.1 Add regression test with mol_000007-shaped input (protomers that are tautomers of each other, e.g. multi-azole neutral molecule) running both tautomer paths in RDKit-only/dry-run mode: assert no exact duplicates among selected tautomers, stereoisomers, and final records, and assert refill-to-cap behavior. Verify: `pytest -k duplicate` passes
- [x] 5.2 Verify the real run output: re-check `runs/Vaclav_set1_p8_t10_s24`-style fixture or rerun analysis script and confirm `final_variants.sdf` has no duplicate groups per the report criterion. Verify: duplicate-groups count is 0
- [x] 5.3 Update `docs/plausible_variant_workflow.md` and `docs/limitations.md` with the identity policy (toolkit-local exact key; why Standard InChI is not used; dedupe after every expansion step; fill-to-cap semantics and resume note). Verify: docs render and statements match `identity.py` docstring
- [x] 5.4 Run full test suite and linters. Verify: `pytest`, `ruff check src tests`, `mypy src` all pass
