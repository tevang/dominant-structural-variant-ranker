# Design: fix-gui-molecule-depictions

## Context

`render_molecules` in `src/dsvr/gui/ui/views.py` draws up to 30 ranked variants as RDKit
SVGs in a 3-column grid. Since the environment moved to Streamlit 1.62, the rendering
call `st.markdown(svg, unsafe_allow_html=True)` no longer displays anything: current
Streamlit sanitizes markdown output regardless of the unsafe flag, stripping `<svg>`
in the browser. `AppTest` cannot observe this (no browser, no sanitizer), which is why
`tests/test_gui.py` stayed green while the feature was broken.

## Decision: `st.image(svg)` (base64 SVG data URI) as the rendering primitive

Original plan (superseded — see verification note below): `st.html(svg)`.

Options considered:

1. **`st.image(svg)`** (chosen; documented fallback of the original plan). Streamlit
   converts SVG strings to a `data:image/svg+xml;base64,...` URL client-side-safe path
   (`image_utils.py` detects the `<?xml`/`<svg` prolog and base64-encodes), so the SVG
   reaches the browser as an `<img>` source and never passes through any HTML
   sanitizer. SVG support in `st.image` long predates the 1.33 floor. Vector rendering
   and browser zoom are preserved (`<img>` of an SVG scales without rasterizing).
2. `st.html(svg)`. **Empirically broken**: Streamlit 1.62's Html component sanitizes
   with DOMPurify `USE_PROFILES:{html:true}, FORCE_BODY:true`, and DOMPurify's
   html-only profile drops `<svg>` and its whole subtree. Running Streamlit's own
   shipped purify bundle in headless Chromium on a real RDKit SVG returned an empty
   string (1988 → 2 bytes). AppTest does not run the sanitizer, so unit tests alone
   cannot distinguish the two — the manual browser check caught this.
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

- **Sanitizer behavior is invisible to AppTest.** Lesson from the st.html attempt:
  only a real browser runs DOMPurify. Mitigation: the rendering path uses `st.image`
  (no sanitizer in the path at all), and the end-to-end check (tasks 3.2) drives a
  real browser and counts SVG `<img>` elements in the DOM — not just AppTest blocks.
- **AppTest regressing silently again.** The tests assert image elements with a
  `data:image/svg+xml` source, and that no inline `<svg` leaks into markdown or html
  elements, so a reversion to any inline-SVG primitive fails the suite.
- **XML prolog (`<?xml ... encoding='iso-8859-1'?>`).** `image_utils._is_svg` accepts
  prolog-prefixed SVGs, so `st.image` handles RDKit output unchanged.

## Version floor

The `gui` extra already declares `streamlit>=1.33`; SVG-string support in `st.image`
predates that floor. The spec records the floor as a requirement so a future floor
lowering (or a remove-and-readd of the extra) is caught during spec review.
