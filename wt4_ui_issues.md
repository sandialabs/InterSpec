# Wt4 UI Issues vs Wt3 Reference

## Background

InterSpec is being upgraded from Wt 3.7.1 to Wt 4.x (currently targeting 4.12.6). This is a significant
framework upgrade — Wt 4 changed many layout, CSS, and widget APIs compared to Wt 3. The majority of the
application functionality has been ported, but a number of GUI rendering and layout issues remain. This file
documents the identified remaining visual/functional differences so they can be systematically fixed.

Most issues are expected to be caused by one of:
- **CSS changes**: Wt 4 changed or removed CSS class names, box-model defaults, or stylesheet rules that the
  application relied on. Check `InterSpec_resources/` for app-specific CSS, and compare with how Wt 3 vs Wt 4
  generates HTML/CSS for the relevant widget types.
- **Layout code (C++)**: `WGridLayout`, `WHBoxLayout`, `WVBoxLayout`, and similar container widgets changed
  behavior in Wt 4. In particular, Wt 4 containers in a layout only render children that are explicitly
  managed by the layout — children added to the container directly (not through the layout) may be invisible.
  Look for places where widgets are added with `addWidget()` on the container instead of through the layout
  object, or where `setMinimumSize`/`resize` calls that Wt 3 needed are now incorrect.

**Fix quality**: Well-formed, robust, maintainable fixes are strongly preferred over workarounds. Fixes should
not break other widgets or areas of the application. Avoid hacks that special-case Wt version numbers. Where
the root cause is a layout or CSS pattern used in multiple places, fix it consistently across the codebase
rather than patching each occurrence individually.


**Fix log**: As issues are resolved, record how each was fixed in `wt4_ui_issues_fixes.md`. Only include an
entry once a fix is confirmed working. Do not include partial attempts or abandoned approaches for items that
are not yet fixed.

---

## Build and Run Instructions

### Building the Wt 4 version

The Wt 4 build lives in `build_wt4/` (or `build_vscode/` for the VSCode-integrated build). Build with:

```bash
cd /Users/wcjohns/coding/InterSpec_wt4/build_wt4
cmake --build .
# or equivalently:
ninja
```

### Running the Wt 4 server (port 8085)

```bash
cd /Users/wcjohns/coding/InterSpec_wt4/build_wt4
./InterSpec --docroot . --http-address 127.0.0.1 --http-port 8085 \
    -c ./data/config/wt_config_localweb.xml
```

To run in the background and log output:
```bash
./InterSpec --docroot . --http-address 127.0.0.1 --http-port 8085 \
    -c ./data/config/wt_config_localweb.xml > /tmp/interspec_wt4.log 2>&1 &
```

If the server becomes unresponsive, kill and restart it:
```bash
pkill -f "InterSpec.*8085"
# then relaunch as above
```

### Running the Wt 3 reference server (port 8081)

The Wt 3.7.1 reference build is a separately compiled binary. 
It is currently running in a similar way on port 8081 (its being run with CWD of `/Users/wcjohns/coding/InterSpec_master/build_xcode` with the executable being located in the RelWithDebInfo subdirectory)
It is used as the visual ground-truth for all comparisons in this document.

### Comparing the two versions in Chrome

Claude Code's Chrome integration (`claude --chrome`) can automate visual comparison:

1. Open `http://localhost:8085` in one Chrome tab (Wt 4, the version being fixed).
2. Open `http://localhost:8081` in a second Chrome tab (Wt 3, the reference).
3. Both servers start with a sample spectrum loaded. Navigate to the same tab or dialog in each, take
   screenshots, and compare side by side.
4. Each navigation to `http://localhost:8085` creates a fresh Wt session, which is useful for resetting
   state before reproducing an issue.
5. Browser console logs (`mcp__claude-in-chrome__read_console_messages`) and the network tab can help
   diagnose JavaScript errors or failed Wt widget updates.

---

## Issues

All originally identified issues (1–17) have been confirmed fixed as of 2026-04-08. See
`wt4_ui_issues_fixes.md` for details on each fix.

---

## Dialogs and Areas Checked — No Issues Found

- Units Converter: identical in both versions
- Flux Tool: identical in both versions
- 1/r2 Calculator: identical in both versions
- Energy Range Sum: identical in both versions
- Math/Command Terminal: appears as tab in bottom panel, identical in both versions
- Help menu: identical items in both versions (Welcome..., Help Contents..., Notification Log,
  Options, Language, About InterSpec...)
- View menu submenus (Chart Options, Peak Labels, Detectors): identical items in both versions
- Activity/Shielding Fit: opens and displays correctly in both versions
- Detector Response Select: opens correctly in both versions
- Make Detector Response: opens correctly in both versions
- Peak editor (double-click a fitted peak): opens and displays identically in both versions
- Right-click context menu on spectrum: works correctly in both versions
- Quick MDA (click near a reference line without a fitted peak): works correctly in both versions
- Peak Manager tab: identical in both versions (note: "Peak Fit Opts" collapsible sidebar on the right is new Wt4-only functionality, not a bug)
- Nuclide Search tab: results table and column headers appear identical in both versions
- Energy Calibration tab: layout and controls appear identical in both versions
- File Parameters dialog: group box (fieldset) borders and labels appear identical in both versions
- Dose Calc dialog (Shielding mode): Thickness input field and material selector appear identical in both versions
- Spectrum Files tab: file info metadata layout appears identical in both versions

---

## New Issues Found — 2026-04-08

### Issue 18 — SimpleDialog title bar: white text on dark background in Wt4

**Where:** Reference Photopeaks tab → type "Co60" in Nuclide field → press Enter → click "more info" button

**Wt4 behavior:** The SimpleDialog title bar ("More Info on Co60") displays **bold white text on a dark grey background**, visually identical to an AuxWindow title bar.

**Wt3 behavior:** The SimpleDialog title ("More Info on Co60") is rendered as **bold black text on a plain white/light background** with no background fill — it looks like an in-page heading.

**Likely cause:** In Wt4, the `SimpleDialog` title element is receiving the same CSS styling (`.Wt-dialog .titlebar` or `.dialog-title`) as the `AuxWindow` title bar, rather than its own distinct simpler style. The CSS class used for the SimpleDialog title or its parent container may overlap with the AuxWindow titlebar styles in the Wt4 stylesheet.

**Files to investigate:** `SimpleDialog.cpp/h`, `AuxWindow.cpp/h`, `InterSpec_resources/*.css`

---

### Issue 19 — AuxWindow footer: "HTML Report" link appears active (blue) with no data in RelActManual

**Where:** Tools → Isotopics from peaks → AuxWindow footer

**Wt4 behavior:** The "HTML Report" link/button in the Isotopics from peaks AuxWindow footer is rendered in **blue (active/enabled style)** even when no isotopics data has been computed yet (the tool is in its initial/empty state).

**Wt3 behavior:** The "HTML Report" link is rendered in **grey (disabled style)** when no isotopics data is present.

**Likely cause:** The enabled/disabled CSS state for footer link buttons in AuxWindow is not being correctly applied in Wt4. The Wt4 CSS may not have the same `.disabled` styling for `WAnchor` elements inside the AuxWindow footer, or the C++ code that calls `setEnabled(false)` on the link is not producing the expected visual state.

**Files to investigate:** `RelActManualGui.cpp/h`, `AuxWindow.cpp/h`, `InterSpec_resources/InterSpec.css` (look for `.Wt-dialog .footer` or `.auxwindow-footer` disabled link styles)

---

### Issue 20 — "Select Nuclide To Add" dialog: auto-shows in Wt4, element list too narrow, presented inline

**Where:** Tools menu → Nuclide Decay Info → (dialog opens)

**Wt4 behavior (three related problems):**
1. **Auto-show**: The "Select Nuclide To Add" sub-dialog appears automatically as soon as the Nuclide Decay Information AuxWindow opens, without the user clicking the "Add Nuclide..." button.
2. **Element list column too narrow**: The element name column in the list is truncated — names show as "hydrog", "beryliu", "nitroge", "potassi", "magnes", "aluminu", "phosph" instead of full names.
3. **Inline presentation**: The sub-dialog appears as an inline popup overlapping the parent AuxWindow (no standalone dialog frame/title bar from the OS or Wt), whereas in Wt3 it renders as a separate, clearly bounded dialog with its own distinct title bar.

**Wt3 behavior:**
1. The "Select Nuclide To Add" dialog only appears when the user explicitly clicks "Add Nuclide..." — it does not auto-show.
2. The element name column is wide enough to show full names (hydrogen, beryllium, nitrogen, potassium, etc.).
3. The dialog renders as a proper standalone dialog with a visible dark-grey title bar labelled "Select Nuclide To Add".

**Likely cause:**
- **Auto-show**: The constructor or initialization code for `NuclideDecayInfo` or its "add nuclide" sub-dialog may be calling `show()` unconditionally in Wt4, or there is a signal connection that fires on initial display.
- **Column width**: The `WTable` or list widget containing element names is not being given sufficient width in Wt4's layout, possibly because a `WGridLayout` or fixed-size container is not expanding as expected.
- **Inline presentation**: The sub-dialog may be a `SimpleDialog` or `WDialog` that is positioned relative to the parent AuxWindow in Wt4 rather than being centered/modal over the whole page.

**Files to investigate:** `NuclideDecayInfo.cpp/h`, associated CSS, `SimpleDialog.cpp/h`

## Issue 21 — All of the "More Actions" dialogs from the "Energy Calibration" tab have rendering issues

**Status (2026-05-07)**: Fixed.  See `wt4_ui_issues_fixes.md`.

## Issue 22 — The tool-tabs rendering when clicked on other than default tab looks grey and missing outline

**Status (2026-05-07)**: Not reproducible against current code.  See `wt4_ui_issues_fixes.md`.

## Issue 23 — Clicking help icon on a tab/tool causes crash

**Status (2026-05-07)**: Not reproducible against current code.  Help dialog opens correctly
from both the global help icon and dialog-footer help buttons.  See `wt4_ui_issues_fixes.md`.

---

## Dialogs surveyed 2026-05-07 — no issues found

After resolving Issues 18, 20, 21:

- Reference Photopeaks → "more info" SimpleDialog (Issue 18 fix)
- Energy Calibration → Linearize / Truncate / Combine Channels / To FRF (Issue 21 fix)
- Nuclide Decay Info + "Add Nuclide..." sub-dialog (Issue 20 fix)
- Activity/Shielding Fit (already in earlier "no issues" list)
- Isotopics from peaks (HTML Report link state matches Wt3 — Issue 19)
- Isotopics by nuclides (Relative Act. Isotopics) opens with all sections rendered
- Gamma XS Calc — width-fitted, all rows displayed
- Dose Calc — Dose / Activity / Distance / Shielding modes all open with full content
- Detection Confidence Tool — chart, inputs, gamma-line table all visible
- Spectrum File Query Tool — query rules, columns and footer buttons render
- Help → About InterSpec ("Disclaimers, Licenses, Credit, and Contact") — tabbed sidebar OK
- Help → Options → Color Themes — theme list, color pickers, footer OK
- Help → Options menu — submenu items match Wt3
- Edit menu — items match Wt3
- Tool tabs (Spectrum Files / Peak Manager / Reference Photopeaks / Energy Calibration /
  Nuclide Search) — active-tab outline matches Wt3 across selections



---

## Issues found 2026-08-14 by diffing the Wt 3.7.1 and Wt 4.13.2 sources

These differ in kind from Issues 1-23 above: they were **not** found by looking at the two apps
side by side, but by mechanically diffing the two framework trees
(`/Users/wcjohns/install/wt-3.7.1_src_code/` vs `/Users/wcjohns/install/wt-4.13.2_src_code/`) for
constructors that newly set an **inline style**, and then grepping InterSpec's CSS for rules those
inline styles now silently override. An inline style beats any stylesheet rule, so these are app CSS
rules that have been dead since the migration with no error anywhere.

Each was confirmed by reading both trees plus the InterSpec code/CSS. Items already fixed are listed
at the end for reference. **None of the items in this section are fixed yet.**

### Issue 24 — `WPanel::setTitleBar()` no longer applies the `titlebar` class (LLM card headers unstyled)

**Where:** LLM Assistant → any conversation/turn card header.

**Mechanism:** Wt3 `WPanel.C:119-127` applied `PanelTitleBarRole` from `setTitleBar(true)`. Wt4
`WPanel.C:125-136` does not; only `setTitle()` schedules it (`WPanel.C:98-99` →
`WCssTheme.C:340-342`). `LlmInteractionDisplay` calls `setCollapsible(true)` and populates
`titleBarWidget()` but never `setTitle()`, so the class is never added.

**Symptom:** header loses its background/border/font-weight, and the "…" menu button jumps from the
right edge to immediately after the status text (the `margin-left:auto` rule dies).

**Dead rules:** `LlmToolGui.css:82`, `:93`, `:98`; `InterSpec.css:836`, `:844`.

**Fix:** `wApp->theme()->apply( this, titleBarWidget(), Wt::PanelTitleBar );` after `setCollapsible(true)`.

### Issue 25 — `WSuggestionPopup` lost its inline `z-index: 10000` (suggestions render under dialogs)

**Where:** any nuclide-entry field inside a dialog, e.g. right-click a peak → Peak Editor → type
`u23` in Nuclide.

**Mechanism:** Wt3 `WSuggestionPopup.C:114` wrote `"z-index: 10000; display:none; overflow:auto"`;
Wt4 `:109` dropped the z-index. `calcZIndex()` is unchanged, so stacking order now decides and a
popup created early sits below a dialog opened later.

**Already partly worked around:** `AuxWindow.cpp:162` bumps `.suggestion`, but that class is set only
by `ShieldMaterialSuggestion.cpp:72`. Unprotected: `PeakEdit.cpp:444`, `DetectionLimitTool.cpp:1170`,
`RelActAutoGuiNuclide.cpp:741`, `ColorThemeWidget.cpp:568`, `NuclideSourceEnter.cpp:97`,
`ReferencePhotopeakDisplay.cpp:971`.

**Fix:** widen the `AuxWindow.cpp:162` selector to `'.suggestion, .Wt-suggest'`.

### Issue 26 — `WFormWidget` no longer applies validation styling on first render

**Mechanism:** Wt3 `WFormWidget.C:162-173` ran the validator during `render(RenderFull)`; Wt4
`:141-149` deleted that block, and `validate()` styles only `if (isRendered())` (`:355-361`).

**Symptom:** InterSpec attaches validators in constructors (72 `setValidator` sites), so a field that
starts empty or out-of-range renders plain white and only turns pink after the first interaction.

**Fix:** call `validate()` from the owning widget's `render(RenderFlag::Full)`.

### Issue 27 — `WSpinBox` arrow hit-zone widened to 22px while CSS still reserves 16px

**Where:** worst on phone (`?isphone=1`) → right-click a peak → Quick MDA → click the right edge of
"Num FWHM".

**Mechanism:** Wt3 `js/WSpinBox.js:178` used `offsetWidth - 16`; Wt4 `:207` uses
`bootstrapVersion < 4 && xy.x > offsetWidth - 22`, and `bootstrapVersion` stays `-1` under
`WCssTheme`, so the 22px branch always runs. Theme `padding-right` stayed 16px (`wt.css:624-629`).

**Symptom:** the last 6px of the *text* area behaves as an arrow — ~40% of the 55px phone spin boxes
(`DetectionLimitSimple.css:82`, `:217`).

**Fix:** `input.Wt-spinbox { padding-right: 22px; }` and widen those two rules.

### Issue 28 — Two more `WStackedWidget` overflow exposures

Same root cause as the five fixed on 2026-08-14 (Wt4's ctor sets inline `overflow:hidden`):
- `ExportSpecFile.css:66` — the phone-only `ExportSpecFileTabs` (`ExportSpecFile.cpp:1098-1106`);
  tab content over 300px is clipped **on phones**.
- `D3TimeChart.css:132` — `.D3TimeChartFilters > div > .Wt-stack`.

Also inherited by `WTabWidget`, whose contents stack is a plain `WStackedWidget` and whose
`setOverflow` forwards to the outer wrapper instead (`WTabWidget.C:47-52`, `:229-233`). Seven
instances never touch `contentsStack()`: `InterSpec.cpp:6957` (main tool tabs),
`UseInfoWindow.cpp:212`, `:455`, `MakeDrf.cpp:1731`, `ShieldingSourceDisplay.cpp:3321`,
`ExportSpecFile.cpp:1105`, `CompactFileManager.cpp:163`, `SpecMeasManager.cpp:7178`.

**Fix:** `stack->setOverflow( Overflow::Auto, Orientation::Vertical );`, or
`tabs->contentsStack()->setOverflow(...)`.

### Issue 29 — Box layouts became CSS flexbox, killing stylesheet sizing of dialog bodies

**Mechanism:** Wt4 defaults layouts to `LayoutImplementation::Flex` (`WLayout.C:18`); `WGridLayout`
still forces the old JS impl, but `WDialog` puts title/body/footer in a `WVBoxLayout` in **both**
trees, so every AuxWindow and SimpleDialog now runs `FlexLayoutImpl`. `FlexLayoutImpl.js:65-73`
reads size limits from `from.style.*` (**inline only**) and then writes `height:auto`,
`maxWidth:100%`, `maxHeight:100%` inline on the child — so a limit coming from a stylesheet is
discarded rather than transferred.

**Dead rules** — all `max-height` on `WDialog::contents()` (== `.body` == `.AuxWindow-content`).
**All of these were removed on 2026-08-14** and replaced by the C++ backstop below:
`ExportSpecFile.css`, `SimpleDialog.css:62`, `:78`, `:87`; `BatchGuiWidget.css:4`;
`InjaLogDialog.css:9`; `RefSpectraWidget.css:6`; `ShieldingSourceDisplay.css:8`;
`FitPeaksForNuclidesGui.css:13`.

Two corrections to that list, from the 2026-08-14 survey below: `DrfSelect.css:221` was listed in
error — it is `max-width` on `.DrfFileSelectMainDesc`, an ordinary inner div that flex never touches,
so it is live (the `.DrfFileSelectDialog .body` rule at `:212` sets padding only). And
`DetectionLimitSimple.css:305`/`:313` are dead for a second, unrelated reason: `SimpleMdaBody`
appears nowhere but in that stylesheet (`grep -rn SimpleMdaBody src/ InterSpec/`), so the selector
matches no element. The class was never added on the C++ side when the rule landed in `ad13333d`.
Those two are still there (an AuxWindow, which has its own backstop — see below).

**Dead C++, not just CSS:** `SimpleDialog::setMaximumSize()` used to build a per-instance
`#<id> .body` max-width/max-height rule; its comment claimed it "wins without `!important`", which is
true against other stylesheets but false against an inline style, so it was a no-op on the
scrollable body. Removed 2026-08-14.

**Fix, shipped 2026-08-14: `SimpleDialog::updateBodySizeForWindow()`.** Set the limit on the *widget*
rather than in a stylesheet. `FlexLayoutImpl.js`'s `copySizeLimits()` reads the body's **inline**
`max-height` and copies it onto the flex item that wraps the body, *then* rewrites the body's own to
`100%` — which now resolves against that wrapper, so an inline value survives the round trip while a
stylesheet one is discarded. Measured on a live dialog: setting `contents()->setMaximumSize()` to
629 px puts `max-height: 629px` on the wrapper and gives the body a computed `100%` = 629.15 px.
`FlexLayoutImpl.C` re-asserts the value through `layout.resizeItem(...)` whenever the C++ size
changes, so later updates work too.

Three consequences worth knowing:
- The limit has to be *arithmetic* (`0.95*renderedHeight() - chrome`), because `WLength` cannot
  express `calc(95vh - 90px)` — see the worked example below.
- Because it is arithmetic rather than `95vh`, it goes stale when the browser window changes size,
  so `InterSpec::layoutSizeChanged()` walks `m_trackedDialogs` and re-applies it. Verified live:
  with the Export dialog open, shrinking the window from 837 px to 417 px took the body from 400 px
  to 306 px (= `0.95*417 - 90`) immediately.
- A dialog needing a different allowance says so with `setBodyChromeHeight()` (Export uses 20 px on
  phone, where there is no title bar or footer), and one wanting a specific height rather than
  sizing to content uses `setBodyPreferredHeight()` (Export asks for 400 px). Both are re-clamped on
  resize.

**AuxWindow needs no equivalent — it already had one.** `AuxWindowResizeToFitOnScreen` (at show) and
`AuxWindowOnDomResize` (a `window` resize listener installed once, `AuxWindow.cpp:383-463`) shrink
any AuxWindow bigger than the window via `dialog.wtObj.onresize(...)` and set `overflow-y: auto` on
its body. Verified by injecting 60 lines into an open Gamma XS window in a 757 px viewport: the
window stayed at 753 px, the body capped at 681 px and scrolled, footer at 754 px.

**DECIDED 2026-08-14: do NOT use `WLayout::setDefaultImplementation(LayoutImplementation::JavaScript)`.**
That one-liner would restore the Wt3 layout engine app-wide and retire this whole issue, but the
project's direction is to migrate to flex layout where appropriate rather than pin the framework to
its old engine. So each case here gets fixed on its merits: size from C++ (which flex honours), or
restructure the CSS to be flex-native. Do not re-propose the global switch.

#### Worked example, 2026-08-14: the Export Spectrum File dialog

Symptom: the dialog grew to ~706 px (nearly the whole screen) and hung off the bottom, because the
File Format list never produced a scrollbar. Confirmed the mechanism above end-to-end, and settled
three open questions about the prescribed fixes. `!important` on the existing rules was tried first
and did work (706 px → 428 px), but was taken back the same day in favour of doing the arithmetic in
C++ (`setBodyChromeHeight()` / `setBodyPreferredHeight()` in `ExportSpecFileWindow`'s constructor),
and those CSS rules were deleted. Same result — 428 px with a 400 px body, verified at 837 px, 757 px
and 357 px viewports, the last of which correctly clamps the 400 px request down to 249 px.

**The intended sizing was already written, in CSS, and had simply stopped applying.**
`git blame` puts `height: min(calc(95vh - 90px), 400px)` at commit `72819438`, 2023-09-01. Nothing
about this dialog's sizing has changed since: `ExportSpecFile.css` was last touched before the Wt4
merge, and the migration commits changed 4 and 10 lines of `ExportSpecFile.cpp`, none of them
sizing. The C++ sets width only (`setMaxWidth(95vw)`, `setMinimumSize(650px, Auto)`); height is
`Auto` everywhere. So under Wt3 this dialog's height came from that one stylesheet rule, and under
Wt4 it comes from its content.

**Measured, not inferred.** On the live dialog `.body` carries
`style="flex: 1 1 auto; height: auto; max-width: 100%; max-height: 100%;"`, and its *computed*
`max-height` is `100%` — not the `calc(95vh - 90px)` the stylesheet asks for. A plain stylesheet
`height` on `.body` had no effect at all; the same rule with `!important` took, and survived a
forced reflow (FlexLayoutImpl rewrites its inline style, `!important` still wins).

**Constraint on the "size from C++" fix: `WLength` cannot express `calc()`, `min()` or `clamp()`.**
`WLength::parseCssString` (`WLength.C`) is `strtod` plus a unit suffix, so anything it cannot parse
falls back to `auto`. C++ sizing can therefore only supply a single value+unit — fine for a plain
`400px` or `85vh`, but it cannot reproduce a viewport expression with a fixed offset such as
`calc(95vh - 90px)`, where the offset does not scale with the viewport. Where the intended size
needs that arithmetic, `!important` on the stylesheet rule is the only faithful option.

**The JavaScript engine is not an escape hatch for dialog bodies — tested, does not work.**
This is the per-layout variant, not the global switch ruled out above, and it is *reachable*:
`StdGridLayoutImpl2` still ships in 4.13.2, `WBoxLayout::setImplementation()` builds it whenever
`preferredImplementation() != Flex` (`WBoxLayout::implementationIsFlexLayout()`), and the layout is
one level up from `contents()` — the ancestry is `contents()` → `WContainerWidget 'dialog-layout'`
(holds the `WVBoxLayout`) → `WTemplate` → the dialog — so
`dynamic_cast<WContainerWidget *>( contents()->parent() )->layout()` gets it with no Wt patching.
Both routes render wrong:
- `setPreferredImplementation(JavaScript)` *after* construction: sizing becomes exactly right (the
  inline flex styles vanish, the stylesheet 400 px applies, the list scrolls) but the engine leaves
  `visibility: hidden` on `.body` and never clears it — blank dialog, persists through a forced
  resize.
- Layout built as JS *from the start* (flipping the default around the `SimpleDialog::make` call):
  visible, but the engine writes `width: 15px` on `.body`; the title wraps down a sliver and the
  content is unreachable.

The second failure explains both: the JS engine measures against definite parent dimensions, and in
Wt4 the dialog is a `WTemplate` sized by CSS (`max-width: 50vw`, the C++ `min-width`) rather than by
the explicit pixel sizes Wt3 handed it. Restoring the engine would mean also restoring how the
dialog gets sized — most of the migration, not a switch. Do not retry this per-dialog either.

**Red herring worth not re-chasing:** the dialog also sat ~130 px too low and hung off the bottom.
That was not a separate positioning bug — Wt centres a dialog when it is shown, using the height it
has at that instant, and this one grew afterwards. Fixing the height fixed the position; a
`centerDialog()` re-centre added for it was verified unnecessary and removed.

#### Survey, 2026-08-14: what the other affected dialogs look like today

Every dialog carrying one of the dead rules was opened in the running app and measured (viewport
813 px high, so the intended cap computes to 682.35 px). **None of them is visibly broken today, and
the reason is uniform: each one already computes the very same cap in C++.** `BatchGuiWidget.cpp:109-126`,
`RefSpectraWidget.cpp:154-171` and `FitPeaksForNuclidesGui.cpp:750-766` all set
`m_widget->setHeight( 0.95*renderedHeight() - 90 )` (capped at 650/500/750 px in the landscape
branch); `InjaLogDialog.cpp:228` does `resize( 95%, 95% )`; `BackPeakPreviewDialog` sizes its chart
to `min(0.45*renderedHeight(), 450)`. `0.95*h - 90` *is* `calc(95vh - 90px)`, so the stylesheet rule
was always belt-and-braces — which is why its loss went unnoticed. Measured example: Reference
Spectra's body carries `min-height: 410px; height: 431px` inline, exactly `min(500, 0.95*813) - 90`.
Export Spectrum File was the one dialog that relied on the CSS alone, hence the only one that broke.

**The real exposure is the generic rule, `SimpleDialog.css:62`.** That is the app-wide safety net for
every auto-sized `SimpleDialog` — confirms, warnings, long error messages, `DrfFileSelectDialog`
(`DrfSelect.cpp:3045`, which sets no height at all) — none of which size themselves from C++.
Measured on a live dialog (File → Enter URL, content padded to ~1260 px):

| | body `max-height` | body height | dialog height | result |
|---|---|---|---|---|
| today | `100%` (inline wins) | 1261 px | 772 px (pinned at `max-height:95vh`) | ~490 px of content **and the whole footer** clipped away by `.simple-dialog{overflow:hidden}`; **no scrollbar** |
| rule revived with `!important` | 682.35 px | 682 px | 751 px | body scrolls internally, Cancel/Okay visible |

So the failure mode here is worse than Export's: not merely a too-tall dialog but a modal whose
buttons are unreachable — Escape is the only way out. Verified with the dialog re-centred, to rule
out the stale-centring red herring below.

`AuxWindow` bodies are laid out by the same flex impl (Gamma XS, Dose Calc, DRF Select and Nuclide
Decay Info all read `max-height: 100%` inline), so any future `max-height` written for an
`.AuxWindow-content` in a stylesheet will be discarded too.

**Both follow-ups from this survey, resolved 2026-08-14.** The generic net became
`SimpleDialog::updateBodySizeForWindow()` rather than `!important` (see the fix above); re-run of the
same experiment afterwards, in a 357 px window with 1261 px of content: body capped at 249 px and
scrolling, Cancel/Okay on screen. And the Simple MDA rule needs no `addStyleClass("SimpleMdaBody")`
after all — `DetectionLimitSimpleWindow` is an `AuxWindow`, so it is already covered by the
AuxWindow JS backstop; the two orphan CSS rules are just dead weight.

**Still open, and pre-existing (not caused by any of this):** a `SimpleDialog` is positioned once,
when shown, with an inline `top` in pixels. Shrink the browser window while one is open and it stays
where it was — measured after 837 px → 417 px, a 334 px dialog sat at `top: 206px` and hung 123 px off
the bottom, footer included. AuxWindows re-centre themselves for exactly this reason
(`AuxWindowOnDomResize`); SimpleDialogs have no equivalent. Sizing is now right in that situation;
position is not.

### Issue 30 — Three more tools with a layout-on-self hosted in a layout-less `contents()`

Same shape as the Multi-File Calibration dialog fixed on 2026-08-14: the widget puts a `WGridLayout`
on itself but is added to a parent with no layout and no definite height, so Wt4's JS layout resolves
the height over several passes (visible jitter) and can settle taller than the dialog. Ranked by
whether the layout has a stretch row (the dangerous variant):

| Tool | layout-on-self | mitigation today |
|---|---|---|
| **External RID** (Tools → External RID) | `RemoteRid.cpp:2438`, `setRowStretch(0,1)` | **none** |
| **Nuclide Decay → Add Nuclide…** | `DecaySelectNuclideDiv.cpp:367`, `setRowStretch(0,1)` | only `setMaximumSize`; this is the root cause behind Issue 20 above |
| **Gamma XS Calc** | `GammaXsGui.cpp:107`, no row stretch | `contents()->setOverflow(Auto,Vertical)` at `:905`, so content stays reachable; jitter only |

**Fix:** host in `stretcher()` and give the window a definite size — preferably via the C++
`resize()` rather than `AuxWindow::resizeWindow()`, which is raw JS that does not establish the
layout chain.

### Issue 31 — `WMenu::select()` now auto-selects the parent item

**Mechanism:** Wt4 `WMenu.C:317-322` selects the owning parent item before selecting the child; no
Wt3 equivalent (`WMenu.C:328`). `select()` emits `triggered()` when the index changes
(`WMenu.C:348-360`). Every InterSpec submenu sets `m_parentItem` (`PopupDiv.cpp:1326`).

**Symptom:** picking an entry inside a submenu also fires `triggered()` on the parent item that owns
it, and marks the parent selected. Trigger: right-click a peak → "Change Nuclide" → pick a nuclide.

**Fix:** `parentItem->setSelectable(false)` in `PopupDivMenu::addPopupMenuItem()` — the new Wt4 code
explicitly honours that flag.

### Issue 32 — Smaller confirmed deltas

- **`GammaCountDialog` width lost.** Wt4's default theme adds
  `.Wt-dialog .Wt-fill-width .body { width:100% }` (`wt.css:300-303`), specificity (0,3,0), which
  outranks `GammaCountDialog.css:1-9` (0,2,0). Tools → Energy Range Sum. Fix: `contents()->setWidth(285)`.
- **`.Wt-tooltip { max-width: 280px }` is dead** — Wt4 `js/ToolTip.js:161` writes `maxWidth` inline
  (Wt3 never did). Long help tooltips no longer wrap at 280px. 313 `attachToolTipOn` call sites.
  Since Wt hardcodes it in JS, `!important` is justified here; alternatively cap the inner div.
- **Submenus clamped to a scrolling parent.** Wt4 `Wt.js:1827-1837` makes `fitToWindow` use the
  parent's visible rectangle when the parent is not BODY/domRoot. `.PopupDivMenu` is exactly such a
  scrolling ancestor (`InterSpec.css:1076-1088`). Only `.PopupDivMenu.AppMenu` is protected
  (`InterSpec.css:1094-1099`); plain context menus are not.
- **`Wt-hide-scrollbar` relies on `scrollbar-width:none`** (Wt4 `WTableView.C:94`), unsupported by
  Safari before 18.2 — and per CLAUDE.md the macOS target is Safari 15.4. A stray horizontal
  scrollbar can appear in table-view headers on the native macOS build. Fix:
  `.Wt-hide-scrollbar::-webkit-scrollbar { display:none; }`.
- **`WPanel` collapse icon** is now a CSS-background `Wt-collapse-button` class instead of two
  `<img>`s. Nothing breaks, but dark theme should now add
  `.Wt-collapse-button { filter: invert(100%); }` (Make DRF, Isotopics from peaks, LLM tool).
- **`WText::setPadding` now emits the `padding` shorthand** with unset sides as `0` (Wt3 emitted
  longhands and discarded top/bottom). So `setPadding(x, Side::Left)` on a **WText** zeroes any
  stylesheet-supplied top/right/bottom. Three sites, none currently killing a live rule:
  `EnergyCalAddActions.cpp:1360`, `PeakSearchGuiUtils.cpp:653`, `LlmSubAgentFollowup.cpp:210`.
  `WContainerWidget` is unaffected (no delta). Rule: never mix `WText::setPadding(single side)` with
  CSS padding on the same WText. Stale comment to update: `QrCode.cpp:463`.
- **`WAnchor` with a null link** now emits `href="#"` + `Wt-no-default` instead of omitting `href`,
  so such anchors get the pointer cursor, focus ring and a tab stop. Makes the comment at
  `LlmToolGui.cpp:191-195` wrong. Forward risk: a null-link anchor with no `clicked()` handler would
  navigate to `#` and clobber `window.location.hash`, which InterSpec uses for its internal path.

### Dead-by-attrition CSS (safe to delete; classes emitted by neither Wt tree)

`InterSpec.css:1413-1417` (`br.Wt-tabs-clear`), `LlmToolGui.css:247`
(`.LlmThinkingContent > .Wt-text`), `InterSpec.css:1633` (`.Wt-suggest-item`).
Also dead `#include`s: `CompactFileManager.cpp:36` and `DecayActivityDiv.cpp:44` (`WSlider`, never
instantiated).

### Non-UI items outstanding (tracked here so they are not lost)

Lifetime/crash bugs, not rendering — full write-up with triggers in `/tmp/wt4_followup_prompt.md`:
1. Cat-A tool windows outlive their owner on **File → Clear Session…** and then dereference a freed
   tool (`~EnergyCalTool` is empty; `~InterSpec` misses `m_llmToolWindow`, `m_terminalWindow` and the
   EnergyCal preserve window; `DecayChainChart` has no destructor at all).
2. Decay-tool CSV export: a dropped `wApp->bind()` guard gives a cross-thread use-after-free if the
   dialog is closed mid-stream (`DecayActivityDiv.cpp:1258`, posted at `:1447`). Five more
   dropped-guard sites listed there.
3. Right-click popup menus leak on mobile — `PopupDivMenu::isHidden()` returns true on mobile, so
   Wt4's `WPopupMenu::done()` returns before emitting `aboutToHide()` and the deferred cleanup never
   runs.
4. Dose Calc runs its handler twice per click; "Add peak…" can leave an unclosable modal on a
   degenerate spectrum.

### Already fixed on 2026-08-13/14 (for reference, do not re-report)

`WStackedWidget` inline-overflow in BatchGuiWidget + DrfSelect + RelActAutoGui (×2) + RemoteRid +
EnergyCalPreserveWindow; Multi-File Calibration dialog clipping/jitter; Peak Manager continuum-type
editor throwing `bad_any_cast`; the `delete`-on-parent-owned-widget family; `observing_ptr` captured
across thread boundaries; deferred widget destruction happening on an io-service thread; the
"Keep Previous Calibration?" **No** button doing nothing (it was wired straight to
`AuxWindow::emitReject`, which suppresses its deferred emit unless the dialog is already hidden).


## Issues found 2026-08-15 while verifying the upgrade/Wt4 → upgrade/DetEff merge

All of these are **pre-existing** — each was checked against the pre-merge tips and is not merge
damage. None are fixed. Found by driving the merged app (desktop + phone) and by a DOM/layout audit
(offscreen windows, non-scrolling overflow, children escaping clipping parents).

### Phone mode tears down the session on load — `?isphone=1` is unusable

Reproducible on every load (2/2 on the merged build; the same hang reproduces on a pre-merge binary,
so this is not from the merge). Desktop is unaffected and keeps working after the phone session dies.

Wt emits `Wt4_13_2.$('<id>').setAttribute('title','Update nuclide search.')` for a DOM node that does
not exist, the JS throws `Cannot read properties of null (reading 'setAttribute')`, and the session is
removed a millisecond later.

The string is an undo/redo **step description** (`IsotopeSearchByEnergy.cpp:2346`), routed to an Edit
menu item's tooltip by `UndoRedoManager::m_undoMenuToolTipUpdate` →
`InterSpec.cpp:6844` (`undoMenu->setToolTip(...)`). So the `PopupDivMenuItem` C++ object outlives its
DOM node in the phone menu build and a later `setToolTip()` addresses a stale id. Same family as the
already-fixed `BringAboveDialogs` phone teardown (`PopupDiv.cpp:693-697`) and as item 3 in the
"Non-UI items outstanding" list above: a widget that is alive in C++ but no longer rendered.

### "Modify Detector Response" tabs render as a vertical stack on desktop

`DrfModifyWidget.cpp:89` does `addStyleClass( "VerticalNavMenu HeavyNavMenu HorizontalMenu" )`, but
`ul.HorizontalMenu` is only defined in `InterSpecMobileCommon.css`, which `InterSpecApp` loads only
when `isMobile()`. On desktop nothing counteracts `VerticalNavMenu`, so the four tabs (Name &
Description / Geometry & MC / FWHM / Uncertainty) render as full-width stacked grey bars, each 768 px
wide, instead of a tab strip. It looks correct on a phone/tablet.

`DrfSelect.cpp:2319` gets this right for the same menu style by adding a fourth class,
`DetEditMenuHorizontal` — which is why the Detector Response Select strip one dialog up renders
horizontally. `LicenseAndDisclaimersWindow.cpp:144` uses the same three-class combination and is
worth checking for the same symptom.

### MakeDrf "Peak Fit Prefs" overflows its bordered panel

`.MakeDrfOptions { width: 150px; }` (`MakeDrf.css:6`) is a hard width, but the `PeakFitDetPrefsGui`
block inside it needs ~196 px, so the Det. Type / FWHM Method / Skew Type combos stick ~48 px past the
panel's right border and overlap the chart. Present on both pre-merge tips; `MakeDrf.css` was not
touched by either branch.

### "Create Detector Response Function" opens far taller than its content

The dialog comes up ~810 px tall with ~290 px of content, leaving a large empty band above the footer.
Same family as the dialog-sizing work recorded earlier in this file (Wt 4 sizing a dialog to something
other than its content).

### Non-UI: every stock GADRAS DRF now hashes differently than intended

`DetectorPeakResponse::computeHash()` (`src/DetectorPeakResponse.cpp:722-723`) folds
`m_totalEfficiency` in under a comment saying the optional fields are hashed "only ... when present, so
legacy DRFs keep their existing hash values". Unlike its neighbours (`eff_uncert`, `m_measuredPoints`),
the guard is a bare pointer check with no emptiness test — and `parseEfficiencyCsvFile` populates
`m_totalEfficiency` from the GADRAS `PTOT` column, which every shipped detector has (62 non-zero rows
in `data/GenericGadrasDetectors/HPGe 40%/Efficiency.csv` alone). So the stated intent does not hold for
any stock DRF.

Verified: the guard, the parser path, and the non-zero PTOT data. **Not** verified end-to-end: the
downstream consequence (the hash is the DB key — `DrfSelect.cpp:5279-5283`, `:5447-5455` — so stock
DRFs would re-insert as duplicates and a `UseDrfPref` lookup could throw).
