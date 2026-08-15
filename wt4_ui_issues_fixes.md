# Wt4 UI Issues — Fix Log

This file records how each issue from `wt4_ui_issues.md` was resolved. Only confirmed, working fixes are
recorded here. Do not add entries for partial attempts or approaches that were abandoned. If an item was
fixed as a side effect of fixing another issue, note that in both entries.

Issues are referenced by their number in `wt4_ui_issues.md`.

---

## Issues 1, 10 — Map menu item and Show Map button missing

**Status**: Fixed

**Root cause**: Both `USE_GOOGLE_MAP` and `USE_LEAFLET_MAP` CMake options were set to `OFF` in the Wt4
build configuration. The code is conditionally compiled with `#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )`.

**Fix**: User enabled the Leaflet map option in the CMake configuration. Additionally, the newly-enabled
`LeafletRadMap.cpp` had several Wt3→Wt4 porting issues that were fixed:
- `new WText(msg, parent)` and `new WCheckBox(text, parent)` → `parent->addNew<WText>(msg)` /
  `parent->addNew<WCheckBox>(text)` (Wt4 removed parent-in-constructor).
- `Signal<...>(this)` → `Signal<...>()` (Wt4 signals don't take owner).
- `Wt::Cursor::WhatsThisCursor` → `Wt::Cursor::WhatsThis`.
- `Wt::RenderFlag::RenderFull` → `Wt::RenderFlag::Full` (via `flags.test()`).
- `WObject::sender()` removed in Wt4 — `InterSpec::handleLeafletWarningWindowClose()` and
  `handleLeafletMapClosed()` simplified to use member variables directly.

**Files**: `src/LeafletRadMap.cpp`, `src/InterSpec.cpp`

---

## Issue 2 — Reference Photopeaks thickness input missing

**Status**: Fixed (same root cause and fix as Issues 8 and 9)

See Issue 8/9 entry below.

---

## Issue 3 — Energy Calibration polynomial coefficient rows missing

**Status**: Fixed

**Root cause**: The `CoefDisplay` class constructor accepted a `WContainerWidget *parent` parameter but
never passed it to the `WContainerWidget()` base class. In Wt3, `new Widget(parent)` would automatically
add the widget as a child of the parent. In Wt4, this pattern was removed — constructors no longer accept
a parent parameter. The `CoefDisplay` widgets were created as orphans and never added to `m_coefficients`,
so they never rendered.

**Fix**:
- Removed the unused `parent` parameter from `CoefDisplay` constructor.
- Changed `new CoefDisplay(coefnum, m_coefficients)` to
  `m_coefficients->addNew<CoefDisplay>(coefnum)` so widgets are properly added to the container.

**Files**: `src/EnergyCalTool.cpp` (lines ~890, ~1328)

---

## Issues 4–7 — Gamma XS Calc dialog rendering issues

**Status**: Fixed

**Root cause**: The `GammaXsGui` widget uses a `WGridLayout` with 3 columns (label, value, units). In
Wt4, `WGridLayout` uses absolute positioning for child widgets. Without column stretch configuration and
with insufficient dialog width, the value and units columns overlapped, input fields were too narrow to
display their text, and the material formula was truncated.

**Fix**:
- Added `m_layout->setColumnStretch(1, 1)` so the value column gets available space.
- Set `setMinimumSize(420, WLength::Auto)` on `GammaXsWindow` to ensure adequate dialog width.
- Added `.GammaXsGui { min-width: 350px; }` CSS rule as a safety net.

**Files**: `src/GammaXsGui.cpp` (lines ~108, ~907), `InterSpec_resources/GammaXsGui.css`

---

## Issues 8, 9 — Thickness input missing in Dose Calc and Detection Confidence Tool

**Status**: Fixed

**Root cause**: The `ShieldingSelect` widget's dimension containers (spherical, cylindrical, rectangular)
used `WGridLayout` internally. In Wt4, `WGridLayout` uses absolute positioning, which gives the layout
container zero intrinsic height. When placed inside the `WStackedWidget` (`m_dimensionsStack`), which has
`overflow: hidden`, the zero-height layout container was fully clipped. Only a thin border line was visible.

**Fix**: Replaced `WGridLayout` with CSS flex layout in the dimension containers:
- Changed the `setupDimEdit` lambda to accept a `WContainerWidget*` instead of `WGridLayout*`, creating
  a flex row container (`.DimEditRow`) for each dimension edit (label + input + optional fit checkbox).
- Removed `setLayout(std::make_unique<WGridLayout>())` from spherical, cylindrical, and rectangular divs.
- Added CSS classes `.DimDiv` and `.DimEditRow` with `display: flex` styling.
- Material summary widget is now added to the first `DimEditRow` instead of through the grid layout.

**Files**: `src/ShieldingSelect.cpp` (lines ~2395–2500), `InterSpec_resources/ShieldingSelect.css`

**Notes**: This fix also resolves Issue 2 (Reference Photopeaks thickness) since all three dialogs use
the same `ShieldingSelect` widget.

---

## Issue 11 — Reference Photopeaks "more info" dialog title shows "??title??"

**Status**: Fixed

**Root cause**: In Wt4, `WDialog` creates a caption `WTemplate` in the title bar with a `${title}`
placeholder. If `setWindowTitle()` is never called, the template variable is unbound, and Wt4's
`WTemplate::handleUnresolvedVariable()` outputs `"??" + varName + "??"` → `"??title??"`.
`SimpleDialog::init()` intentionally skipped calling `setWindowTitle()` (to avoid Wt's default `<h4>`
formatting), but this left the template unresolved in Wt4.

**Fix**: Added `setWindowTitle("")` after `setTitleBarEnabled(true)` in `SimpleDialog::init()`. This
binds an empty string to the template variable (rendering an invisible empty `<h4>`), while the custom
`m_title` WText continues to display the actual title.

**Files**: `src/SimpleDialog.cpp` (line ~161)

---

## Issue 12 — Reference Photopeaks "more info" dialog cannot be dismissed

**Status**: Fixed

**Root cause**: The `MoreNuclideInfoWindow` dialog did have a Close button (added via
`SimpleDialog::addButton("Close")`), but the original issue description stated it had none. Testing
confirmed the Close button IS present and functional. What was missing was Escape-key support — the
dialog did not call `rejectWhenEscapePressed()`.

**Fix**: Added `rejectWhenEscapePressed(true)` call in `MoreNuclideInfoWindow` constructor.

**Files**: `src/MoreNuclideInfoDisplay.cpp` (line ~708)

**Notes**: Escape key may not work when focus is inside the scrollable content area, but the Close button
always works.

---

## Issues 13, 14 — Isotopics dialog options sections missing controls

**Status**: Fixed

**Root cause**: The "Spectrum and Peak Options" and "Relative Efficiency Curve Options" sections used
`WGroupBox`. In Wt4, `WGroupBox` always has an internal `WVBoxLayout` (created by
`WGroupBox::setLogicalLayout(nullptr)` override). This causes `createDomChildren()` to take the
`if (layout_)` path and only render layout-managed children. Direct children added via `addNew<>()` are
silently dropped from the DOM.

**Fix**: Replaced `WGroupBox` with plain `WContainerWidget` using `<fieldset>` HTML tag and a `<legend>`
child element for the section title. This is the same workaround used in `ShieldingSourceDisplay.cpp`.
Direct children added via `addNew<>()` render correctly because `WContainerWidget` (without an explicit
layout) renders all direct children.

**Files**: `src/RelActAutoGui.cpp` (lines ~682, ~874)

---

## Issue 15 — Isotopics tab buttons on wrong side

**Status**: Fixed

**Root cause**: In the `RelActAutoGui` constructor, the `WStackedWidget` was added to `upper_div` before
the `WMenu`. In Wt4, DOM order determines visual order, so the stack appeared first (on the left) and the
menu tabs appeared second (on the right). In Wt3, the `WMenu` rendered its tabs on the left regardless of
DOM insertion order.

**Fix**: Reversed the DOM insertion order so the `WMenu` is added to `upper_div` before the
`WStackedWidget`. The `WStackedWidget` is created as a standalone `unique_ptr` first (so it can be passed
to the `WMenu` constructor), then the `WMenu` is added to the container, and finally the stack is added
after via `addWidget()`.

**Files**: `src/RelActAutoGui.cpp` (lines ~519–530)

---

## Issue 16 — Isotopics "add free peaks" link missing

**Status**: Fixed

**Root cause**: This was a side effect of the `WGroupBox` issue (Issues 13, 14). Once the surrounding
containers were corrected and the `WGroupBox` replaced, the Energy Ranges section rendered fully,
including the "add free peaks" anchor at its bottom-right corner.

**Fix**: No additional change needed beyond the `WGroupBox` → `WContainerWidget`/`<fieldset>` fix
applied for Issues 13 and 14 in `src/RelActAutoGui.cpp`.

**Files**: `src/RelActAutoGui.cpp`

---

## Issue 20 — Nuclide Decay Info: Select-Nuclide sub-dialog auto-shows; element names truncated

**Status**: Fixed

**Three sub-problems**, all rooted in Wt4 layout/visibility differences:

1. **Auto-show on parent open**: Wt4's `WDialog::setHidden(true)` invoked from C++ at the end
   of construction does not reliably suppress the first paint — the rendered DOM ends up with
   `visibility:visible` regardless.  In Wt3 the equivalent `hide()` call was honored, so the
   sub-dialog stayed hidden until the user clicked "Add Nuclide...".

2. **Element column truncated**: `DecaySelectNuclide` lays out `m_elementSelection` /
   `m_massSelection` with a `WGridLayout`.  In Wt4 the layout uses absolute positioning sized
   from the dialog's auto-size (which collapses to title-bar width), giving the element list
   ~61 px — only wide enough for "hydrog", "berylliu", etc.

3. **Inline appearance**: a side-effect of (1) — the dialog showed at `inset:0px 0px` because
   `centerWindow()` was called once before being immediately hidden, and Wt4 didn't honor the
   subsequent state changes.

**Fix**:
- Added a new `AuxWindowProperties::StartHidden` flag.  In the `AuxWindow` constructor, when
  this flag is set we still call `WDialog::setHidden(true)` (so C++ state is correct) and
  emit a one-shot `setTimeout(...,0)` JS that sets `el.style.visibility='hidden'` after the
  initial DOM append (the DOM element doesn't exist when `doJavaScript` fires synchronously).
  `AuxWindow::show()` clears that inline visibility so subsequent shows paint normally.  Use
  `visibility` (not `display`) so child layouts still get measured for first render —
  otherwise inner WGridLayout-positioned widgets end up with no width when shown later.
- `DecayActivityDiv::createWidgets()` constructs `m_nuclideSelectDialog` with the new
  `StartHidden` flag, removes the now-redundant `centerWindow()` + `hide()` calls, and moves
  `centerWindow()` into the `Add Nuclide...` click handler so the dialog is centered each
  time it is shown.
- `DecaySelectNuclide::initNuclideMenu()` sets explicit `setMinimumSize` on
  `m_elementSelection` (95 px) and `m_massSelection` (60 px) so the WGridLayout reserves
  enough width for full element names and 4-digit mass numbers.

**Files**: `InterSpec/AuxWindow.h`, `src/AuxWindow.cpp`, `src/DecayActivityDiv.cpp`,
`src/DecaySelectNuclideDiv.cpp`

---

## Issue 23 — Help-icon click causes crash

**Status**: Not reproducible

**Investigation**: Clicking the global help icon in the bottom-left of the page opens the
"InterSpec Help" dialog correctly with the table of contents and content panes populated.
Clicking the help-icon button in dialog footers (e.g. `Activity/Shielding Fit`'s footer help
button) also opens the help dialog cleanly.  No crash, no JS errors in the console — only the
unrelated routine "Failed to find scale factor being adjusted" log lines from
`SpectrumChartD3.js`.  The original issue does not reproduce against the current code base.

**Files**: (none)

---

## Issue 22 — Tool-tab rendering on non-default tabs

**Status**: No longer reproducible

**Investigation**: Comparing tool-tab rendering side-by-side in Wt3 (port 8082) and Wt4
(port 8080) for "Spectrum Files", "Peak Manager", "Reference Photopeaks", "Energy Calibration",
and "Nuclide Search" — including tab transitions — both versions show the same active-tab
styling: bold text inside a tab-shape outline, while inactive tabs are plain bold text on the
chart panel background.  The original description ("looks grey and missing outline") does not
reproduce against the current code; it may have been resolved as a side effect of one of the
earlier WMenu/WTabWidget porting fixes.  Closing without further change.

**Files**: (none)

---

## Issue 19 — HTML Report link not visually disabled in RelActManual

**Status**: Not actually a regression

**Investigation**: The Wt3 and Wt4 source for `RelActManualGui` use the same construction
pattern for `m_downloadHtmlReport` (`WPushButton` / `WAnchor` with a `WLink` to the HTML
resource); neither code path calls `setEnabled(false)` when there is no result. Side-by-side
comparison of port 8082 (Wt3) and port 8080 (Wt4) confirms both render the link in the
active blue style at all times. The original issue was a misobservation. No change needed.

**Files**: (none)

---

## Issue 18 — SimpleDialog title bar: white text on dark background in Wt4

**Status**: Fixed

**Root cause**: Wt4's default theme stylesheet (`themes/default/wt.css`) applies
`.Wt-dialog .titlebar { background:#888888; color:#FFFFFF; padding: 2px 6px 3px; }`.
`SimpleDialog::init()` already calls `removeStyleClass("titlebar")` and adds its own `.title`
class, but Wt4's theme re-adds the `titlebar` class on render — so the dialog's title kept the
dark grey AuxWindow-style appearance instead of the plain heading SimpleDialog wants.

**Fix**: Added a CSS override in `SimpleDialog.css` targeting `.Wt-dialog.simple-dialog .titlebar`
(specificity 0,2,0 — same as Wt's rule, but its later position in the cascade wins) to reset
background, color, and padding so the existing `.simple-dialog .title` rule controls the
heading appearance.

**Files**: `InterSpec_resources/SimpleDialog.css`

---

## Issue 21 — Energy Calibration "More Actions" dialogs have wrong width / clipped content

**Status**: Fixed

**Root cause**: `EnergyCalAddActionsWindow` placed each tool widget inside `AuxWindow::stretcher()`,
which sets a `WGridLayout` on the dialog's contents and `overflow: hidden` on the container.  In
Wt4, `WGridLayout` positions children absolutely with widths derived from the layout container's
own width.  Since the dialog auto-sizes to its title-bar width (the absolutely-positioned children
contribute zero to layout's intrinsic size), the column ends up narrower than the inner table —
which is then clipped by `overflow: hidden`.  Two of the cases (`ConvertToFrf`, `ConvertToPoly`)
worked around this with a hardcoded `setWidth(425)`.

**Fix**: Replaced `stretcher()->addWidget(...)` with `contents()->addNew<...>(...)` for all six
`MoreActionsIndex` cases.  Without a layout, the dialog's content container uses normal block
flow and the dialog auto-sizes to the natural width of its content — which is the desired
behavior.  Removed the hardcoded `setWidth(425)` for the FRF/Poly conversion cases as well; their
content now drives the width.

**Files**: `src/EnergyCalAddActions.cpp` (lines ~194–230)

---

## Issue 17 — Isotopics Close button does not close the dialog

**Status**: Fixed

**Root cause**: Two problems:
1. `RelActAutoGui::createWindow()` connected the `finished()` signal to the wrong handler
   (`InterSpec::closeShieldingSourceFit` instead of `InterSpec::handleRelActAutoClose`). This caused the
   wrong tool's state to be saved and the wrong window to be deleted.
2. `AuxWindow::setHidden()` calls `emitReject()` inline, which emits the `finished()` signal during the
   hide flow. When a `finished()` handler deletes the window (via `deleteAuxWindow` or `removeChild`),
   the stack unwinds through destroyed `AuxWindow` methods — a use-after-free that crashes the server.
   This is a fundamental incompatibility between `AuxWindow::emitReject()` and the standard Wt4
   `done()`/`removeChild()` dialog lifecycle pattern.

**Fix**:
- Removed the wrong `finished()` connection from `createWindow()` (the correct connection is already
  made in `InterSpec::relActAutoWindow()`).
- Changed the close button to directly call `InterSpec::handleRelActAutoClose()` instead of
  `AuxWindow::hide()`, bypassing the `emitReject()` → `finished()` chain.
- `handleRelActAutoClose()` saves state, hides the dialog via `WDialog::setHidden(true)` (bypassing
  `AuxWindow::setHidden`), and clears the observing pointers. The widget tree is hidden but not deleted.

**Files**: `src/RelActAutoGui.cpp` (line ~348, ~384), `src/InterSpec.cpp` (`handleRelActAutoClose`)

**Notes**: The hidden widget tree leaks until the session ends. Proper cleanup requires refactoring
`AuxWindow::setHidden()` to not call `emitReject()` inline, which would enable the standard Wt4 pattern
of `done()` + `removeChild()` from the `finished()` handler. This is a broader AuxWindow issue that
likely affects other dialog close flows in Wt4 as well.

---

## AuxWindow title bar shorter; footer button vertical alignment off (2026-05-07)

**Status**: Fixed

**Root cause**: In Wt3 the WDialog laid out its title-bar / contents / footer with absolute
positioning, so InterSpec's `.titlebar { height: 25px }` and the footer's C++
`setHeight(45)` were honored.  In Wt4 the WDialog wraps those three children in a
WVBoxLayout whose flexbox implementation inlines `height: auto; flex: 1 1 auto;` on each
child, overriding `height` (CSS or inline).  Result:

- Title bar collapsed from 25 px to 24 px (content height).
- Footer collapsed from 45 px to 40 px and lost its `.modal-footer { padding-top: 13px }`,
  because Wt4's theme also adds the `footer` class to the dialog's footer element (even
  after `setStyleClass("modal-footer")`), so `.Wt-dialog .footer { padding: 4px 6px 4px }`
  (specificity 0,2,0, source-ordered after InterSpec.css's `.modal-footer` 0,1,0) won.
- The Close button sat 5 px below the footer top (vs ~15 px in Wt3), so it looked too
  high.

**Fix** (CSS only, in `InterSpec_resources/InterSpec.css`):

- `.titlebar` rule: added `min-height: 25px` so the bar is at least 25 px even when the
  flex layout sets `height: auto`.
- Renamed `.modal-footer` to `.AuxWindow.Wt-dialog .modal-footer` (specificity 0,3,0) to
  outrank `.Wt-dialog .footer`, and added `height: 45px; min-height: 45px;
  box-sizing: border-box;` so the footer keeps its Wt3 dimensions and padding.

After the fix both the title bar and the footer measure the same as in Wt3 (25 px and
45 px), the footer padding is `13px 10px 4px` (Wt3 had `13px 10px 0px` — the Wt4 theme
still injects 4px bottom which is harmless because the buttons are vertical-align middle),
and the Close button sits at the same vertical position as in Wt3.

**Files**: `InterSpec_resources/InterSpec.css` (the `.titlebar` and `.modal-footer` rules)

---

## Tool-tabs / spectrum splitter bar too prominent in Wt4 (2026-05-07)

**Status**: Fixed

**Root cause**: Wt4's default theme added a new rule `.Wt-hrh2 { background-color: #CCC }`
to `themes/default/wt.css` that did not exist in Wt3.  `.Wt-hrh2` is the horizontal
row-resize handle that sits between the chart and the tool-tabs in InterSpec; in Wt3 the
bar was transparent so only the tiny centered `splitter-v.png` glyph was visible, but in
Wt4 the new rule turned the entire 5 px strip light grey, making it look like a deliberate
divider.

**Fix**: Added a same-specificity override in `InterSpec_resources/InterSpec.css`:

```css
.Wt-hrh2 { background-color: transparent; }
```

InterSpec.css is `@import`-ed after the Wt theme, so equal-specificity rules in
InterSpec.css win on source order.  The hover state is unchanged because the theme's
`.Wt-hrh2:hover { background-color: #EEE }` rule is more specific.

**Files**: `InterSpec_resources/InterSpec.css`

---
