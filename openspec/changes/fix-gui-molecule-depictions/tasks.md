# Tasks: fix-gui-molecule-depictions

## 1. Rendering fix

- [x] 1.1 Extract `_render_depictions(subset: pd.DataFrame) -> None` from the
  `"Depictions"` branch of `render_molecules` in `src/dsvr/gui/ui/views.py`,
  preserving: 3-column `st.columns(3)` grid, `subset.head(30)` cap, `variant_id`
  caption per structure, "_(unparseable SMILES)_" placeholder, and the "No SMILES
  column available for depictions." info branch. Verify: `ruff check src` passes;
  `render_molecules` body shows no behavior change outside the extraction
- [x] 1.2 Replace the depiction emission inside `_render_depictions` with
  `st.html(svg)`; remove every `unsafe_allow_html=True` tied to SVG content (keep the
  placeholder/captions on plain `st.caption`/`st.markdown` without HTML). Verify:
  `grep -n "unsafe_allow_html" src/dsvr/gui/` finds no SVG-related occurrence;
  `ruff check src` and `mypy src` pass

## 2. Regression tests

- [x] 2.1 Add an AppTest regression test in `tests/test_gui.py`: open the Molecules
  view on the existing `make_run_dir` fixture (ranked_molecules=11), set the "Show"
  selectbox to "Depictions", rerun, and assert that the number of emitted HTML
  elements containing `<svg` equals `min(30, ranked row count)` and that zero of them
  went through markdown elements. Verify: test passes with the fix; reverting task
  1.2 to `st.markdown(..., unsafe_allow_html=True)` makes it fail
- [x] 2.2 Add a helper-level test: a ranked.csv where one row carries an unparseable
  SMILES string renders the "_(unparseable SMILES)_" placeholder for exactly that
  row while other rows still emit SVG HTML blocks. Verify: `pytest tests/test_gui.py
  -k depict` passes
- [x] 2.3 Run the full GUI test module and linters. Verify: `uv run pytest
  tests/test_gui.py`, `ruff check src tests`, and `mypy src` all pass

## 3. Real-run validation

- [x] 3.1 Headless smoke check via AppTest against real run data:
  `runs/Pavel_set1_p8_t6_s32` (2 molecules, 21 ranked rows) and
  `runs/Vaclav_set2_p8_t10_s24` (12 molecules) — Molecules → Show "Depictions" emits
  one `<svg`-containing HTML block per ranked row on the first page with no
  exceptions. Verify: scripted AppTest run prints zero exceptions and expected block
  counts
- [ ] 3.2 Manual browser check: `dsvr view runs/Pavel_set1_p8_t6_s32
  --no-open-browser` (or with browser) → Molecules → Show "Depictions" displays the
  2D structures visibly in the browser, suitable for the RDKit UGM screenshot.
  Verify: structures are visibly rendered (user confirms)
- [x] 3.3 Run the full test suite. Verify: `uv run pytest` passes end to end
