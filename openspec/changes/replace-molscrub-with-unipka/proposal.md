# Proposal: replace-molscrub-with-unipka

## Why

molscrub generates protomer candidates as structures only — it predicts no protonation thermodynamics, so DSVR's protomer selection relies on a heuristic plausibility score and pH influences candidate filtering in a hand-wavy way ("no rigorous pH population prediction"). Uni-Pka — specifically the error-corrected implementation shipped in the EasyDock project (`ci-lab-cz/easydock` `containers/unipka`) — predicts a free energy for every enumerated microspecies and Boltzmann occupancies at any pH. Making Uni-Pka the default protomer source replaces heuristic selection with prediction-driven, occupancy-based selection and equips every protomer with thermodynamic features usable later for uncertainty prediction alongside features from the other pipeline tools.

## What Changes

- Uni-Pka (EasyDock container implementation, Uni-Mol/Dwars torch model) becomes the **default** `protonation.tool`; molscrub remains selectable (`protonation.tool: molscrub`) with its current heuristic behavior unchanged.
- New Uni-Pka runner executes the tool through a container runtime (Docker or Apptainer; image name or `.sif` path in config), batching all input molecules in one call instead of per-molecule subprocesses.
- Protomer enumeration requests multiple microspecies per molecule (`-n` / `--occupancy`) and selects protomers by predicted occupancy: forms above a configurable occupancy threshold, trimmed by occupancy under `max_protomers_per_molecule`; `keep_best_per_charge` retained.
- Predicted properties stored for downstream uncertainty prediction:
  - per protomer: occupancy at working pH, predicted microstate free energy (dG);
  - per input molecule: top-two occupancy gap, protonation-ensemble entropy at working pH, net-charge population distribution, number of microstates enumerated, distance from working pH to nearest macroscopic pKa transition, isoelectric point.
- The Uni-Pka pH-distribution file (`--distribution-file`, microspecies dG + occupancy over a pH range) becomes a stored run artifact referenced by the run report.
- `dsvr doctor`, workflow step metadata, provenance (`source_software`/`source_command`), stage summaries, reports, docs, configs, and bootstrap scripts updated for the new default tool and its container requirement.
- **BREAKING** for default-config runs: protomer sets and therefore downstream rankings may change versus molscrub-era runs because selection is occupancy-based (provenance records which tool produced each run).
- Verification target: full `prepare-ligands` run on `examples/Vaclavs_heterocycles_set1.smi` (`--max-protomers 8 --tauto-k 10 --max-stereoisomers 24`) completes with zero duplicate final structures.

## Capabilities

### New Capabilities

- `protonation-unipka`: Uni-Pka-backed protomer generation — container execution, batched invocation, occupancy-based enumeration and selection, predicted-property metadata, and pH-distribution artifact.
- `protonation-tool-selection`: dispatch of the protomer generation stage on `protonation.tool` (`unipka` default, `molscrub` alternative), including availability checks and failure behavior per selected tool.

### Modified Capabilities

<!-- No existing capability specs (auto3d-*, run-stabilization) have requirement-level changes. -->

## Impact

- **Code**: new `src/dsvr/runners/unipka_runner.py`; `src/dsvr/chemistry/protonation.py` (dispatch + occupancy selection + metadata); `src/dsvr/config.py` (`protonation.tool` dispatch, new `unipka` options); `src/dsvr/utils/tool_check.py`; `src/dsvr/workflow/steps.py`, `engine.py`, `recovery.py`; `src/dsvr/cli.py` help text; reporting/markdown.
- **Dependencies / runtime**: a container runtime (Docker or Apptainer) plus the Uni-Pka image (`unipka.sif` from Zenodo or built from the EasyDock recipe) becomes a required runtime dependency of the default workflow, replacing the molscrub Python package/CLI requirement (molscrub optional).
- **Configs**: `configs/*.yaml` and `examples/example_config.yaml` gain the Uni-Pka options; `protonation.tool: unipka` becomes the default.
- **Outputs / schema**: `ProtomerRecord.metadata` gains Uni-Pka fields (additive); new per-run distribution artifact; report gains protonation-property sections.
- **Tests / docs**: protonation, doctor, provenance, workflow-smoke, and duplicate-regression tests extended; `docs/` (external_tools, workflow, installation, architecture, limitations, changelog) and `ONBOARDING.md`/`README.md` updated.
- **Out of scope**: removing molscrub, PNG occupancy plots, logP/logD prediction, HPC deployment specifics.
