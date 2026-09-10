# Design: fix-gui-molecule-depictions

## Context

`render_molecules` in `src/dsvr/gui/ui/views.py` draws up to 30 ranked variants as RDKit
SVGs in a 3-column grid. Since the environment moved to Streamlit 1.62, the rendering
call `st.markdown(svg, unsafe_allow_html=True)` no longer displays anything: current
Streamlit sanitizes markdown output regardless of the unsafe flag, stripping `<svg>`
in the browser. `AppTest` cannot observe this (no browser, no sanitizer), which is why
`tests/test_gui.py` stayed green while the feature was broken.

## Decision: `st.html(svg)` as the rendering primitive

Options considered:

1. **`st.html(svg)`** (chosen). Available since Streamlit 1.33 — exactly the declared
   floor of the `gui` extra (`streamlit>=1.33` in pyproject.toml), so no dependency
   change is required. Content is DOMPurify-sanitized client-side (the `<svg>`,
   `<line>`, `<ellipse>`, `<path>`, `<rect>`, `<text>` elements RDKit's
   `MolDraw2DSVG` emits are all on the default allow list), injected inline (not
   iframed), and remains zoomable/selectable as vector content.
2. Base64 `<img src="data:image/svg+xml;base64,...">` inside `st.html` or `st.image`.
   Works but rasterizes-by-reference: text in the SVG is no longer selectable, and the
   wrapper adds size/encoding plumbing for no benefit here. Kept as the documented
   fallback only if `st.html` proves to misbehave on some future Streamlit.
3. Pin Streamlit to the last version that tolerated unsanitized markdown SVG.
   Rejected: pinning backwards to preserve a deprecated behavior blocks security
   fixes and is the wrong direction.

## Structure

- `_render_depictions(subset: pd.DataFrame) -> None` extracted from
  `render_molecules`: locates the SMILES column, lays out the 3-column grid, iterates
  `subset.head(30)`, and for each row emits `st.html(svg)` when
  `smiles_to_svg` returns a string, otherwise `st.caption("_(unparseable SMILES)_")`;
  emits `st.caption(str(row.get("variant_id", "")))` after each structure; keeps the
  existing `st.info("No SMILES column available for depictions.")` branch. All
  constants (3 columns, 30-row cap) are preserved verbatim.
- No change to `depict.py::smiles_to_svg` (caching, size defaults, None semantics).
- No change to any other view, to pagination, or to molecule selection.

## Risks and mitigations

- **DOMPurify strips something RDKit emits.** Mitigation: spec requires the emitted
  HTML blocks to contain `<svg` in tests, and the manual browser check (tasks 3.2)
  covers the real render end-to-end. RDKit's output is plain static SVG with no script,
  so allow-listing holds.
- **XML prolog (`<?xml ... encoding='iso-8859-1'?>`) inside inline HTML.** The prolog
  precedes the `<svg>` element in the fragment; browsers/DOMPurify tolerate and drop
  it in fragment parsing. Verified acceptable for the manual check; if it ever
  surfaces as an issue, stripping the prolog in `_render_depictions` is a one-line
  follow-up (not needed now).
- **AppTest regressing silently again.** The new test asserts the *type* of the
  emitted element (HTML block) and the presence of `<svg` in its body, so a future
  reversion to `st.markdown` fails the suite even though AppTest does not sanitize.

## Version floor

`st.html` exists since Streamlit 1.33. The `gui` extra already declares
`streamlit>=1.33`; the spec records this floor as a requirement so a future floor
lowering (or a remove-and-readd of the extra) is caught during spec review.
