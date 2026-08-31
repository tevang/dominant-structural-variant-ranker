# Uni-Pka container assets

- `unipka.py` — vendored copy of the current EasyDock Uni-Pka script
  (https://github.com/ci-lab-cz/easydock `containers/unipka/unipka.py`, BSD-3-Clause;
  synced from commit `7a281a95b785278541f2ff08b30c678a29952fea`, 2026-08-21).
  DSVR bind-mounts it over `/unipka/unipka.py` at run time because every published
  Zenodo `unipka.sif` build (through record 19627026, 2026-04-19) bakes an older
  script that lacks `-n/--occupancy/--distribution-file`. Override the path with
  `protonation.unipka.script_path` (empty string disables the override).
- `unipka.sif` — local, gitignored. Acquire per `docs/external_tools.md` (Zenodo
  record 19627026, expected size 7713267712 bytes, md5 `64994c54e626ed5eac0dabb4416cc749`)
  or build with `apptainer build unipka.sif` from the EasyDock recipe.
