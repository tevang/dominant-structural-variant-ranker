# Changelog

All notable changes to DSVR are documented here. Format loosely follows
"Keep a Changelog"; dates are ISO.

## Unreleased

### Changed

- **Uni-Pka is the new default protonation tool** (`protonation.tool:
  unipka`). Protomer enumeration requests multiple microspecies per molecule
  and selects by predicted Boltzmann occupancy at the working pH, replacing
  the heuristic plausibility scoring used for molscrub. Default runs are not
  comparable to molscrub-era runs (provenance records the tool).
- Protomers generated with Uni-Pka carry per-form metadata
  (`unipka_occupancy`, `unipka_dg`) and every input molecule gains a
  Uni-Pka summary (top-two occupancy gap, ensemble entropy, charge
  populations, microstate count, nearest pKa transition distance, isoelectric
  point) intended as features for downstream uncertainty analysis.
- The Uni-Pka pH-distribution output is stored per run at
  `enumeration/protomers/unipka_distribution.tsv` and referenced by the run
  report; the report gains a "Protonation Properties (Uni-Pka)" section.
- `dsvr doctor` now checks the Uni-Pka container runtime and image
  (required) and reports molscrub as optional.
- molscrub (`protonation.tool: molscrub`) remains selectable with unchanged
  legacy behavior; `configs/*.yaml` and install docs updated; CI smoke
  configs pin molscrub since no container is available in CI.

### Added

- `protonation.unipka` configuration block: `container`, `runtime`
  (auto/docker/apptainer), `script_path`, `min_occupancy`, `max_forms`,
  `ph_range_*`, `ph_step`, `timeout_seconds`.
- Vendored current EasyDock `unipka.py` at `containers/unipka.py`
  (BSD-3-Clause, see `containers/README.md`). Every published Zenodo
  `unipka.sif` build bakes an outdated script that rejects the occupancy
  flags, so DSVR bind-mounts the vendored copy over `/unipka/unipka.py`
  unless `script_path` is set to another path or to `""`. The default
  `container: unipka` also resolves to a local `containers/unipka.sif`
  when present, so a downloaded image needs no config change.
