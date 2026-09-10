# Tasks: fix-gui-molecule-depictions

## 1. Rendering fix

- [x] 1.1 Extract `_render_depictions(subset: pd.DataFrame) -> None` from the
  `"Depictions"` branch of `render_molecules` in `src/dsvr/gui/ui/views.py`,
  preserving: 3-column `st.columns(3)` grid, `subset.head(30)` cap, `variant_id`
  caption per structure, "_(unparseable SMILES)_" placeholder, and the "No SMILES
  column available for depictions." info branch. Verify: `ruff check src` passes;
  `render_molecules` body shows no behavior change outside the extraction
- [x] 1.2 Replace the depiction emission inside `_render_depictions` with
  `st.image(svg)`; remove every `unsafe_allow_html=True` tied to SVG content (keep the
  placeholder/captions on plain `st.caption`/`st.markdown` without HTML). Verify:
  `grep -n "unsafe_allow_html" src/dsvr/gui/` finds no SVG-related occurrence;
  `ruff check src` and `mypy src` pass. Note: an intermediate `st.html(svg)`
  attempt passed AppTest but failed in real browsers — Streamlit's DOMPurify
  `USE_PROFILES: {html: true}` configuration strips `<svg>` entirely (verified in
  real Chromium: 1988-byte SVG sanitized to 2 bytes). `st.image` bypasses the
  sanitizer by emitting a `data:image/svg+xml;base64` `<img>` instead.

## 2. Regression tests

- [x] 2.1 Add an AppTest regression test in `tests/test_gui.py`: open the Molecules
  view on the existing `make_run_dir` fixture (ranked_molecules=11), set the "Show"
  selectbox to "Depictions", rerun, and assert that the number of emitted image
  elements whose URL starts with `data:image/svg+xml;base64,` equals
  `min(30, ranked row count)` and that zero `<svg` strings appear in markdown or
  html element bodies. Verify: test passes with the fix; reverting task
  1.2 to `st.markdown(..., unsafe_allow_html=True)` makes it fail
- [x] 2.2 Add a helper-level test: a ranked.csv where one row carries an unparseable
  SMILES string renders the "_(unparseable SMILES)_" placeholder for exactly that
  row while other rows still emit SVG data-URI images. Verify: `pytest
  tests/test_gui.py -k depict` passes
- [x] 2.3 Run the full GUI test module and linters. Verify: `uv run pytest
  tests/test_gui.py`, `ruff check src tests`, and `mypy src` all pass

## 3. Real-run validation

- [x] 3.1 Headless smoke check via AppTest against real run data:
  `runs/Pavel_set1_p8_t6_s32` (2 molecules, 21 ranked rows) and
  `runs/Vaclav_set2_p8_t10_s24` (12 molecules) — Molecules → Show "Depictions" emits
  one SVG data-URI image per ranked row with no exceptions. Verify: scripted
  AppTest run prints zero exceptions and expected image counts
- [x] 3.2 Browser check via Playwright (headless Chromium): full
  `dsvr view runs/Pavel_set1_p8_t6_s32` → Molecules → Show "Depictions" renders
  exactly 21 `img[src^='data:image/svg+xml']` elements, each with non-zero natural
  size (260x200), and no exception banners. Evidence screenshot:
  `/tmp/gui_check.png`. Verify: scripted Playwright run prints 21 images, 0
  exceptions
- [x] 3.3 Run the full test suite. Verify: `uv run pytest` passes end to end
