# Workflow Redesign Report

## 1. Previous Runaway Cause

The previous workflow behaved too much like an exhaustive enumerator. It could
generate many protomer, tautomer, stereoisomer, and ETKDG conformer combinations
and then send too many of those structures to CREST/xTB. In one real run this
created about 58,000 XYZ files in roughly three hours.

The main causes were:

- ETKDG conformers were treated as downstream CREST seeds rather than a local
  sampling pool.
- CREST scheduling was tied too closely to raw seed count.
- Tautomer enumeration had no process-level timeout.
- Enantiomeric duplicates could be sent to CREST independently in achiral
  solvent.
- xTB/CENSO/QM refinement was not sufficiently late-stage and budgeted.
- Disk cleanup and progress/audit reporting were too limited.

## 2. New Funnel Architecture

The redesigned workflow is a staged funnel:

1. Read SMILES/SDF input.
2. Generate pH/protomer candidates with molscrub.
3. Enumerate RDKit tautomers with process timeout and fallback.
4. Enumerate RDKit stereoisomers.
5. Score variants with SVPScore.
6. Apply bounded cheap-score selection.
7. Generate ETKDG/Auto3D seed pools.
8. Select only representative low-energy/diverse seeds.
9. Optionally run xTB quick prefilter.
10. Collapse pure enantiomer CREST jobs in achiral solvent.
11. Run CREST/xTB only on survivors.
12. Run xTB thermo only on capped top conformers/variants.
13. Optionally run CENSO on a small final subset.
14. Optionally run Psi4/PySCF on a tiny final subset.
15. Write ranked outputs, audit tables, progress files, and reports.

## 3. New Default Configs

Recommended production entry point:

- `configs/production_balanced.yaml`

Other profiles:

- `configs/production_conservative.yaml`: lower caps for safer local runs.
- `configs/crest_disk_safe.yaml`: focused disk-safe CREST settings.
- `configs/exhaustive_debug.yaml`: intentionally expensive debug profile.

`README.md` now recommends `production_balanced.yaml` for production-style local
use. Exhaustive mode prints a CLI warning because it may generate very large
variant and XYZ counts.

## 4. New Pruning/Penalty Logic

SVPScore is a transparent heuristic triage score. It is not a thermodynamic
model and is not a rigorous population estimator.

The score includes:

- protomer plausibility terms;
- tautomer heuristic terms;
- stereochemical uncertainty terms;
- chemistry sanity penalties;
- computational complexity penalties;
- optional cheap 3D relative-energy penalties.

Selection is now budgeted:

- `max_variants_after_cheap_score_per_molecule`
- `max_variants_for_xtb_prefilter_per_molecule`
- `max_variants_for_crest_per_molecule`
- `max_seeds_per_variant`

Rescue candidates are prioritized inside the configured budget rather than
silently expanding the budget.

The xTB prefilter is optional and enabled in production configs. It runs quick
single-point/optimization scoring without Hessian or thermo, then keeps
low-energy and rescue variants before CREST.

Pure enantiomeric pairs are collapsed for CREST/xTB scheduling in achiral
solvent, while all stereoisomer identities are preserved by mapping the
representative energy back to equivalent enantiomeric partners.

## 5. New Disk Cleanup Policy

The default production policy is compact:

- seed generation writes compact SDF/CSV outputs, not large XYZ trees;
- temporary XYZ files are generated immediately before CREST/xTB;
- CREST parses a capped number of conformers;
- only top conformer XYZ files are retained;
- intermediate patterns such as `struc*.xyz`, `trial*.xyz`, `coord.*`, and
  `*.tmp` are deleted;
- large raw logs may be gzip-compressed;
- disk guard checks run and molecule directory sizes.

Disk-related audit files include:

- `disk_usage_by_stage.csv`
- `crest_job_plan.csv`
- `progress.json`
- `progress.jsonl`

## 6. New Progress Reporting Behavior

The workflow writes durable progress files:

- `progress.json`
- `progress.jsonl`
- `stage_summary.csv`

`dsvr status RUN_DIR` reports:

- last stage and status;
- active command if known;
- latest log tail;
- disk usage;
- XYZ file count;
- counts by stage.

The completed validation run confirmed:

- `runs/conservative_test/progress.jsonl` was written.
- `runs/balanced_test/progress.jsonl` was written.
- `dsvr status runs/revision1_fast_smoke_test` worked.

## 7. Remaining Scientific Limitations

DSVR remains an orchestrator, not a rigorous pH population engine.

Important limitations:

- molscrub pH handling is candidate generation, not rigorous solution
  speciation.
- Cross-protomer population estimates are approximate without micro-pKa or
  proton chemical-potential corrections.
- RDKit canonical tautomer selection is not tautomer stability ranking.
- SVPScore is a filtering heuristic, not final thermodynamics.
- xTB/CREST rankings depend on conformer search quality, semiempirical method
  limitations, implicit solvent assumptions, and parser correctness.
- Enantiomer collapse is valid only for achiral solvent/environment assumptions;
  disable it for chiral binding pockets, chiral catalysts, or chiral separation
  contexts.

## 8. Exact Production Command

Recommended production-balanced command:

```bash
conda activate dsvr
dsvr run examples/test_molecules.smi \
  --config configs/production_balanced.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/production_balanced_water_pH7 \
  --overwrite
```

Safer conservative command:

```bash
conda activate dsvr
dsvr run examples/test_molecules.smi \
  --config configs/production_conservative.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/production_conservative_water_pH7 \
  --overwrite
```

## Validation Summary

Executed in the `dsvr` conda environment:

```bash
python -m pip install -e .
dsvr doctor --json
dsvr run examples/test_molecules_minimal.smi --config configs/production_conservative.yaml --ph 7.0 --solvent water --out runs/conservative_test --overwrite
dsvr run examples/test_molecules_minimal.smi --config configs/production_balanced.yaml --ph 7.0 --solvent water --out runs/balanced_test --overwrite
pytest
ruff check .
mypy src/dsvr
DSVR_RUN_EXTERNAL=1 pytest -m external
```

Results:

- `pytest`: 108 passed, 2 skipped.
- `ruff check .`: passed.
- `mypy src/dsvr`: passed.
- external tests: 1 passed, 109 deselected.
- required external tools detected by doctor: RDKit, molscrub, xTB, CREST.
- optional tools missing: CENSO, Psi4, PySCF.
