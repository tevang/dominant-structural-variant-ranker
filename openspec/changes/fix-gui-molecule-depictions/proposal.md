# Proposal: fix-gui-molecule-depictions

## Why

The run-inspection GUI (`dsvr view`) no longer displays 2D molecular structures. In the
"Molecules" view, selecting Show = "Depictions" produces an empty area in the browser:
no molecule is rendered for any ranked row, on any run directory. The user needs this
view for an RDKit User Group Meeting presentation screenshot.

Root cause (investigated and confirmed before this proposal):

- `src/dsvr/gui/ui/views.py::render_molecules` renders RDKit SVG strings via
  `st.markdown(svg, unsafe_allow_html=True)` (every depicted row, line ~178).
- Streamlit ≥ 1.62 (the version resolved by uv.lock on this machine) sanitizes
  markdown output even when `unsafe_allow_html=True`, stripping the inline `<svg>`
  markup client-side. Code that worked on older Streamlit now renders nothing.
- The SVG generator itself is healthy:
  `src/dsvr/gui/ui/depict.py::smiles_to_svg` produces valid SVG for all rows of real
  run data (verified against `runs/Pavel_set1_p8_t6_s32/ranked.csv` — 21/21 rows
  parse and draw).
- The test suite does not catch the regression because
  `streamlit.testing.v1.AppTest` never exercises the browser-side sanitizer; it only
  asserts that markdown elements exist.
- This is the only affected call site in the codebase — the sole
  `unsafe_allow_html=True` on raw SVG content.

## What Changes

- Render molecule depictions with `st.image(svg)` instead of
  `st.markdown(..., unsafe_allow_html=True)`. `st.image` turns SVG strings into
  base64 `data:image/svg+xml` `<img>` sources, so SVG never passes through a
  client-side sanitizer. (The first attempt, `st.html(svg)`, was disproven by running
  Streamlit's own shipped DOMPurify bundle in a real browser: its html-only profile
  drops `<svg>` entirely — see design.md.)
- Extract a small `_render_depictions(subset)` helper so the depiction branch is
  readable and unit-testable; preserve all existing behavior around it: 3-column grid,
  30-row cap, `variant_id` caption per structure, and the
  "_(unparseable SMILES)_" placeholder for rows whose SMILES fail to parse.
- Add regression coverage that fails if depictions regress to inline-SVG rendering:
  an AppTest-driven test that opens the Molecules view, switches Show to
  "Depictions", and asserts one image element with a `data:image/svg+xml` source per
  depicted ranked row (up to the 30-row cap) and zero inline `<svg` in markdown/html
  elements; plus helper-level tests for the valid/invalid SMILES branches.
- Guard the floor: assert in tests that the `gui` extra keeps `streamlit >= 1.33`
  (the Streamlit floor for the GUI feature set; `st.image` SVG support predates it),
  documented as the requirement in the spec. pyproject.toml already declares
  `streamlit>=1.33` for the `gui` extra, so no dependency change is needed.

## Capabilities

### New Capabilities
- `molecule-depiction-rendering`: how the run-inspection GUI renders 2D molecular
  structures (RDKit SVG via sanitization-safe Streamlit API), including fallback
  rendering for unparseable SMILES and the version floor for the rendering API.

### Modified Capabilities
(none — no archived spec covers the GUI depiction path today; behavior of all other
views is unchanged.)

## Impact

- Code: `src/dsvr/gui/ui/views.py` (rendering call + helper extraction),
  `tests/test_gui.py` (new regression tests).
- No dependency changes: `streamlit>=1.33` already declared; 1.62 is installed.
- Behavior: Depictions view works again in real browsers; all other views unchanged.
  No effect on pipeline outputs, scientific results, or provenance.
