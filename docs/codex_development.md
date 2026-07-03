# Codex Development Guidance

This project contains chemistry-aware orchestration logic and routine engineering tasks. Future Codex work should split model usage by scientific risk and reasoning load.

## Use gpt-5.5

Use `gpt-5.5` for work that changes scientific behavior, workflow assumptions, or chemistry-facing user semantics:

- architecture changes;
- chemistry logic;
- RDKit enumeration logic;
- xTB/CREST/CENSO parser design;
- free-energy and population-ranking logic;
- debugging scientifically suspicious results;
- reviewing workflow assumptions;
- writing user-facing scientific warnings.

These tasks require careful treatment of chemical state definitions, lineage, energy comparability, and documented limitations. Any change in this category should include corresponding tests and documentation updates.

## Use gpt-5.4-mini

Use `gpt-5.4-mini` for low-risk repetitive work after the intended behavior is already specified:

- checking logs for common external-tool failures;
- monitoring long-running subprocesses;
- running pytest/ruff/mypy repeatedly;
- updating routine documentation;
- formatting files;
- checking expected file existence;
- summarizing progress;
- making simple mechanical refactors after gpt-5.5 specifies them.

`gpt-5.4-mini` should not independently introduce new chemistry claims, ranking assumptions, parser semantics, or workflow behavior changes.

## Scientific Claim Discipline

Do not add or strengthen chemistry claims unless the claim is documented in code, tests, and docs. If a method only provides approximate ranking, label it approximate in user-facing outputs and reports.
