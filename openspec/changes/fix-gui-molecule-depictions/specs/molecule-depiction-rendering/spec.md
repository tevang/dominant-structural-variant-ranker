## Purpose

Defines how the run-inspection GUI (`dsvr view`) renders 2D molecular structures in the
Molecules view so that depictions remain visible across supported Streamlit versions,
with explicit fallback behavior for unparseable SMILES and a declared version floor for
the rendering API.

## ADDED Requirements

### Requirement: Depictions render via a sanitization-safe API

The Molecules view SHALL render each 2D depiction by passing the RDKit SVG string to
`st.html`. Depiction SVG content MUST NOT be emitted through `st.markdown` (with or
without `unsafe_allow_html`), because supported Streamlit versions sanitize markdown
output and strip inline `<svg>` markup.

#### Scenario: Depictions are visible in a real browser

- **WHEN** a user opens the Molecules view for a run whose `ranked.csv` contains a
  `smiles` column and selects Show = "Depictions"
- **THEN** every depicted row emits its SVG through an `st.html` element and no
  depiction SVG is emitted through `st.markdown`

#### Scenario: Regression to markdown rendering is caught by tests

- **WHEN** the automated GUI tests open the Molecules view, switch Show to
  "Depictions", and inspect the emitted elements
- **THEN** the test observes one HTML element whose body contains `<svg` per depicted
  ranked row (up to the 30-row cap), and the test fails if depictions are emitted as
  markdown elements instead

### Requirement: Depiction layout and fallbacks are preserved

The refactor MUST preserve the existing depiction presentation: a 3-column grid, a cap
of 30 depicted rows, and a caption carrying the row's `variant_id` beneath each
structure. Rows whose SMILES fail to parse SHALL render the
"_(unparseable SMILES)_" placeholder in place of the structure, and a run without a
SMILES column SHALL show the existing informational message instead of depictions.

#### Scenario: Unparseable SMILES row

- **WHEN** a depicted row's SMILES cannot be parsed by `smiles_to_svg` (returns None)
- **THEN** the view renders the "_(unparseable SMILES)_" placeholder for that row and
  continues rendering remaining rows

#### Scenario: Missing SMILES column

- **WHEN** the selected subset has no `smiles` column
- **THEN** the view shows "No SMILES column available for depictions." and renders no
  depiction elements

### Requirement: Streamlit version floor for the GUI extra

The rendering path relies on `st.html`, available since Streamlit 1.33. The GUI extra
SHALL keep `streamlit>=1.33` as its lower bound, and capability changes touching the
GUI MUST NOT declare compatibility below that floor without replacing the rendering
primitive.

#### Scenario: Dependency declaration review

- **WHEN** pyproject.toml declares the `gui` extra dependencies
- **THEN** the Streamlit requirement is `>=1.33` or higher, matching the rendering
  API in use
