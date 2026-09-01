# Uni-Pka container assets

- `unipka.py` — **patched** vendored copy of the current EasyDock Uni-Pka script
  (https://github.com/ci-lab-cz/easydock `containers/unipka/unipka.py`, BSD-3-Clause;
  based on commit `7a281a95b785278541f2ff08b30c678a29952fea`, 2026-08-21, plus the
  DSVR shared-microstate fix below — it is therefore *not* byte-identical to upstream
  any more).
  DSVR bind-mounts it over `/unipka/unipka.py` at run time because every published
  Zenodo `unipka.sif` build (through record 19627026, 2026-04-19) bakes an older
  script that lacks `-n/--occupancy/--distribution-file`. Override the path with
  `protonation.unipka.script_path` (empty string disables the override).
- `unipka-shared-microstate-fix.patch` — the exact DSVR fork diff (the `---` side is
  upstream file at commit `7a281a9`, the `+++` side is the patched `unipka.py` in this
  directory). Apply with `patch -p1`; re-create after re-syncing upstream with
  `diff -u <upstream>/unipka.py unipka.py > unipka-shared-microstate-fix.patch`.
- `unipka.sif` — local, gitignored. Acquire per `docs/external_tools.md` (Zenodo
  record 19627026, expected size 7713267712 bytes, md5 `64994c54e626ed5eac0dabb4416cc749`)
  or build with `apptainer build unipka.sif` from the EasyDock recipe.

## DSVR fork: shared-microstate (multi-parent) fix

**Upstream bug (Uni-Pka 0.3.2).** The streaming pipeline (`UnipkaStream.process`) keyed
microstates by their plain SMILES and attributed each one to a *single* parent molecule
(`microstate_to_smi[ms] = smi` — last writer wins). When two input molecules of one batch
share microstates by exact SMILES collision — typically a conjugate acid/base pair such as
`aniline.[Cl-]` + `anilinium.[Cl-]`, which share both the neutral form and the cation — the
shared microstates were credited twice to the last-enumerated molecule and zero times to the
first. The first molecule's completion condition (`len(predicted) == total`) could never be
met, so it was silently abandoned and never written to the output. The run exited 0 with
fewer output lines than input lines, and downstream tools (DSVR included) silently fell back
to the unprotonated input state. The same 1:1 assumption also blocked any molecule whose
ensemble lists the same microstate SMILES twice (`total` counts duplicates, a SMILES-keyed
`predicted` dict cannot).

**Fix.** Multi-parent attribution, localised to `UnipkaStream.process` / `UnipkaStream._flush_gpu`:

- the reverse index became many-to-one (`microstate_to_smis: Dict[str, List[str]]`) and every
  currently pending parent of a microstate is credited with its predicted energy;
- each unique microstate SMILES is queued and sent to the GPU exactly once, however many
  molecules share it (prediction dedup preserved);
- predicted energies are kept for the whole run (`predicted_energies`): a molecule that
  registers after a shared microstate was already predicted in an earlier flush is credited
  immediately from that store, and completes on the spot if nothing else is outstanding;
- `total` counts *unique* microstates so duplicate entries within one molecule's ensemble
  cannot block completion;
- a completed molecule's results are remembered (`emitted`): the priority stream does not
  deduplicate, so a repeated input SMILES delivered again after completion emits only its
  not-yet-emitted names instead of repeating all rows (one output row per input row is
  preserved, matching upstream `smi_to_names` semantics);
- molecules still in `pending` when the source is exhausted are now reported with a warning
  naming them, so a truncated output can never be silent again.

**Verified.** `tests/test_unipka_shared_microstate_regression.py` (5 tests) drives the patched
pipeline with a fake predictor and a hand-crafted two-molecule batch sharing both microstates:
conjugate-pair survival, single prediction per shared microstate, intra-molecule duplicate
microstates, late registration across a mid-stream GPU flush, and duplicate-SMILES
re-delivery. All five fail against the pristine upstream copy (the conjugate pair is silently
dropped) and pass against the patched copy. The bug was also reported upstream.

