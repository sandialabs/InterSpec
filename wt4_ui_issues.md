# Wt4 UI Issues vs Wt3 Reference

## Background

InterSpec is being upgraded from Wt 3.7.1 to Wt 4.x (now on 4.13.2). This is a significant
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

Two configured build dirs (CMake, Unix Makefiles). Prefer the Debug one: Release defines `NDEBUG`, which
compiles out Wt's own double-free asserts, so a Wt 4 ownership bug corrupts the heap silently instead of
aborting at the fault.

- `build_wt4_debug/` — Debug (preferred)
- `build_wt4/` — Release

The executable target is `InterSpecExe` (there is no `InterSpec` target; `--target InterSpec` builds
nothing and still exits 0). Keep builds to 4 cores:

```bash
cmake --build /Users/wcjohns/coding/InterSpec_wt4/build_wt4_debug --target InterSpecExe -j4
```

Other sessions share this tree and its build dirs. Rebuilding while someone else's server runs from the
same dir relinks the dylib under it, so check `lsof -nP -iTCP -sTCP:LISTEN | grep InterSpec` first.

### Running the Wt 4 server

Pick a port nobody else is using (8080 and 8082 are often taken):

```bash
cd /Users/wcjohns/coding/InterSpec_wt4/build_wt4_debug
./InterSpec --docroot . --http-address 127.0.0.1 --http-port 8093 \
    -c ./data/config/wt_config_web.xml > /tmp/interspec_wt4.log 2>&1 &
```

To stop it, kill only the process on *your* port. Other InterSpec instances on other ports belong to
other work, so never use `pkill -f InterSpec`:
```bash
PID=$(lsof -nP -iTCP:8093 -sTCP:LISTEN -t); [ -n "$PID" ] && kill $PID
```

### Running the Wt 3 reference server

The Wt 3.7.1 reference is `/Users/wcjohns/coding/InterSpec_master/build_release/InterSpec` (Release,
Wt 3.7.1, last built 2026-08-20). It is the visual ground truth for the comparisons in this document.
It does not start as-is, for two reasons (both checked 2026-09-30):
- its `user_data` is a broken symlink to `../build_xcode/user_data`;
- its `InterSpec_resources` symlinks to the source dir, which lacks the deployed `d3.v3.min.js` /
  `SpectrumChartD3.{js,css}`, so the client dies with "failed loading /InterSpec_resources/d3.v3.min.js".

Run it with a scratch user-data dir and a scratch docroot instead of changing that tree:

```bash
M=/Users/wcjohns/coding/InterSpec_master; D=/tmp/wt3_docroot; mkdir -p $D/InterSpec_resources /tmp/wt3_user_data
for f in $M/InterSpec_resources/*; do ln -sf "$f" $D/InterSpec_resources/; done
for f in d3.v3.min.js SpectrumChartD3.js SpectrumChartD3.css; do ln -sf $M/external_libs/SpecUtils/d3_resources/$f $D/InterSpec_resources/; done
ln -sfn $M/build_release/resources $D/resources; ln -sfn $M/data $D/data
ln -sfn $M/build_release/external_libs $D/external_libs; ln -sfn $M/example_spectra $D/example_spectra
cd $M/build_release && ./InterSpec --docroot $D --http-address 127.0.0.1 --http-port 8094 \
    -c ./data/config/wt_config_web.xml --userdatadir /tmp/wt3_user_data > /tmp/interspec_wt3.log 2>&1 &
```

### Comparing the two versions

**Headless Playwright** needs no browser extension, and it is what the 2026-09-30 re-verification used.
The repo's copy lives in `target/testing/PlaywrightPhoneEmulation/`. Launch it with
`executablePath: ~/Library/Caches/ms-playwright/chromium_headless_shell-1243/chrome-headless-shell-mac-arm64/chrome-headless-shell`,
because the bundled version expects 1223. Things to know:
- Load `?restore=no`, then close the Welcome dialog with its `Close` button (Escape does not close it).
- Menus are `.MenuBar .PopupMenuParentButton`, and the first one is labelled "InterSpec". Samples are
  under InterSpec → Samples → "Ba-133 (16k channel)".
- Phone (`?isphone=1` plus a device descriptor such as `devices['iPhone 14']`) uses `.MobileMenuButton`
  instead, and its first entry is "Spectra".
- Drive the chart from JS: `d=document.querySelector('.D3SpectrumDisplayDiv')`, then
  `d.chart.WtEmit(d.id,{name:'doubleclicked'},energy,0,'',0)` fits a peak and
  `...{name:'rightclicked'},energy,0,pageX,pageY,''` opens the context menu.
- Inject `.Wt-tooltip{display:none !important}`, or tooltips intercept clicks.

**Chrome MCP** (`claude --chrome`) works too: open the Wt 4 and Wt 3 ports in two tabs and compare
screenshots. Each navigation creates a fresh Wt session.

---

## Issues

All originally identified issues (1–17) have been confirmed fixed as of 2026-04-08. See
`wt4_ui_issues_fixes.md` for details on each fix.

### Status at a glance (re-verified 2026-09-30)

Every item below Issue 17 was checked on 2026-09-30 in two ways: against the current source (HEAD
`425c3d7a`), and by driving the app headlessly. The Wt 4 build used was `build_wt4_debug` (built
2026-09-29 22:44); the Wt 3.7.1 reference was `InterSpec_master/build_release`. Views covered were
desktop at 1400×850, plus phone and tablet where noted. Each section below starts with a dated
**Status** line giving the evidence.

| Item | Status |
|---|---|
| 18 SimpleDialog title bar | **Superseded, by design.** SimpleDialog deliberately reuses the AuxWindow chrome (5450eb13). |
| 19 HTML Report link state | Same as Wt 3 |
| 20 Decay "Select Nuclide To Add" | Still fixed |
| 21–23 | Unchanged (fixed / not reproducible) |
| 24 LLM card header lost `titlebar` class | **Open.** Confirmed live. |
| 25 Suggestion popup z-index | **Not a bug.** Wt 4 raises the popup each time it is shown. |
| 26 Validation styling on first render | **Open.** Confirmed live (DRF Select → Formula). |
| 27 Spin-box 22 px arrow zone | **Open.** Confirmed live on desktop and phone. |
| 28 Stacked-widget inline `overflow:hidden` | Latent. Overflow is hidden on both stacks, but nothing is clipped today. |
| 29 Flex dialog-body sizing | Fixed. The "SimpleDialog position on resize" follow-up is now resolved too. Orphan `SimpleMdaBody` CSS remains. |
| 30 Layout-on-self in `contents()` | No jitter seen. Gamma XS opens off-centre (too high). |
| 31 `WMenu::select()` selects the parent | Code unchanged. No observable symptom. |
| 32 Smaller deltas | Tooltip width fixed. Energy Range Sum width still open. Other items are unchanged (see the section). |
| Dead-by-attrition CSS | Still present |
| Non-UI items 1, 4 | Fixed (7760358d) |
| Non-UI item 2 | CSV export fixed. The External RID undo crash and the cross-thread `this` captures in `InterSpec.cpp` were fixed 2026-09-30. |
| Non-UI item 3 (mobile popup leak) | Partly fixed. The leak is now bounded. |
| 2026-08-15: phone teardown | No longer reproduces |
| 2026-08-15: Modify DRF tabs | Fixed (3f1e2e71) |
| 2026-08-15: MakeDrf Peak Fit Prefs overflow | Fixed (10191318 redesign) |
| 2026-08-15: Create DRF dialog tall | **Not a regression.** Wt 3 is identical. |
| 2026-08-15: GADRAS DRF hash | Deliberate, and pinned by tests |

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

**Status (2026-09-30): superseded by design. Do not "fix" this back toward Wt 3.** Commit 5450eb13
(2026-09-22, "Give SimpleDialog the AuxWindow chrome…") deliberately restyled SimpleDialog so that it
reuses AuxWindow's `.titlebar` rule. The goal is consistency: the plain-heading look described below was
the old iOS-alert style, and it made SimpleDialogs read like a different application. The 2026-05-07 CSS
fix in `wt4_ui_issues_fixes.md` has therefore been replaced.

Measured on Edit → Enter URL: the title bar is 25 px, background `#888`, with white text, the same as an
AuxWindow. Footer buttons are ordered by role. The Wt 3 reference still shows the old centred bold
heading with full-width text buttons; that difference is expected now.

**Where:** Reference Photopeaks tab → type "Co60" in Nuclide field → press Enter → click "more info" button

**Wt4 behavior:** The SimpleDialog title bar ("More Info on Co60") displays **bold white text on a dark grey background**, visually identical to an AuxWindow title bar.

**Wt3 behavior:** The SimpleDialog title ("More Info on Co60") is rendered as **bold black text on a plain white/light background** with no background fill — it looks like an in-page heading.

**Likely cause:** In Wt4, the `SimpleDialog` title element is receiving the same CSS styling (`.Wt-dialog .titlebar` or `.dialog-title`) as the `AuxWindow` title bar, rather than its own distinct simpler style. The CSS class used for the SimpleDialog title or its parent container may overlap with the AuxWindow titlebar styles in the Wt4 stylesheet.

**Files to investigate:** `SimpleDialog.cpp/h`, `AuxWindow.cpp/h`, `InterSpec_resources/*.css`

---

### Issue 19 — AuxWindow footer: "HTML Report" link appears active (blue) with no data in RelActManual

**Status (2026-09-30): same as Wt 3.** With a spectrum loaded and no peaks, Tools → Isotopics from peaks
(docked as the "Isotopics" tab) renders the "HTML Report" button the same way in both versions: a
`LinkBtn DownloadBtn` in `rgb(0,0,255)` at opacity 0.7, and not disabled.

**Where:** Tools → Isotopics from peaks → AuxWindow footer

**Wt4 behavior:** The "HTML Report" link/button in the Isotopics from peaks AuxWindow footer is rendered in **blue (active/enabled style)** even when no isotopics data has been computed yet (the tool is in its initial/empty state).

**Wt3 behavior:** The "HTML Report" link is rendered in **grey (disabled style)** when no isotopics data is present.

**Likely cause:** The enabled/disabled CSS state for footer link buttons in AuxWindow is not being correctly applied in Wt4. The Wt4 CSS may not have the same `.disabled` styling for `WAnchor` elements inside the AuxWindow footer, or the C++ code that calls `setEnabled(false)` on the link is not producing the expected visual state.

**Files to investigate:** `RelActManualGui.cpp/h`, `AuxWindow.cpp/h`, `InterSpec_resources/InterSpec.css` (look for `.Wt-dialog .footer` or `.auxwindow-footer` disabled link styles)

---

### Issue 20 — "Select Nuclide To Add" dialog: auto-shows in Wt4, element list too narrow, presented inline

**Status (2026-09-30): still fixed.** The sub-dialog opens only when the add button is clicked, as its
own "Select Nuclide To Add" AuxWindow. Element names show in full (hydrogen … calcium), and its height
settles at 543 px within 120 ms (Wt 3: 539 px). One deliberate difference from Wt 3: the Decay tool's
"Add Nuclide…" and "Remove All" are now icon-only buttons (a ⊕ and a clear icon at the right of the
nuclide strip), and their labels became tooltips (b3539981, `DecayActivityDiv.cpp:1801-1809`).

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

**Re-swept 2026-09-30 (desktop 1400×850, both versions).** Every Tools-menu window was opened and its
height sampled at 0.3, 1.0 and 2.5 s. No window changed height after 300 ms, none extended off-screen,
and none had a `.Wt-invalid` field at open. Sizes match Wt 3 to within 5 px, except for these:
- Energy Range Sum: 407 px wide vs 297 px (Issue 32).
- Gamma XS Calc: 422 px wide vs 304 px, opening at top 38 px vs Wt 3's centred 163 px (Issue 30).
- Activity/Shielding Fit: 1026 px wide vs 1122 px. The redesign changed it, not a regression.
- Isotopics from peaks: opens as a docked tab while tool tabs are shown, as it does in Wt 3.



---

## Issues found 2026-08-14 by diffing the Wt 3.7.1 and Wt 4.13.2 sources

These differ in kind from Issues 1-23 above: they were **not** found by looking at the two apps
side by side, but by mechanically diffing the two framework trees
(`/Users/wcjohns/install/wt-3.7.1_src_code/` vs `/Users/wcjohns/install/wt-4.13.2_src_code/`) for
constructors that newly set an **inline style**, and then grepping InterSpec's CSS for rules those
inline styles now silently override. An inline style beats any stylesheet rule, so these are app CSS
rules that have been dead since the migration with no error anywhere.

Each was confirmed by reading both trees plus the InterSpec code/CSS. Items already fixed are listed
at the end for reference. (As of 2026-09-30 some of these are resolved or turned out not to be
bugs; see each item's **Status** line and the table at the top of "Issues".)

### Issue 24 — `WPanel::setTitleBar()` no longer applies the `titlebar` class (LLM card headers unstyled)

**Status (2026-09-30): open, confirmed live.** The problem is only in the outer per-question card,
`LlmInteractionDisplay`. Its constructor (`LlmInteractionDisplay.cpp:1650`) calls `setCollapsible(true)`
at `:1672` and never calls `setTitle()`. The per-turn cards (`LlmInteractionTurnDisplay` subclasses)
do call `setTitle()` before `setCollapsible()`, so they get the class and are fine.

Measured after sending one question (it fails with a network error, because no API key is set):
- The outer card's title-bar `<div>` has **no class**: transparent background, 0 padding, weight 400,
  `display:block`. Its "…" menu button sits 1071 px from the right edge of a 1368 px bar.
- The error-turn card's bar has class `titlebar`, `display:flex` and 8/10 px padding, and its button is
  10 px from the right edge.

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

**Status (2026-09-30): not a bug. Close it.** The analysis below is right that the constructor's inline
z-index was dropped, but it misses that Wt 4 now raises the popup every time it is shown.
`WSuggestionPopup.js:117-135`: `showPopup()` calls `bringToFront()`, which sets
`zIndex = WT.maxZIndex() + 1`. Creation order therefore does not matter. Measured in both versions,
with `elementFromPoint` landing on the popup in every case:
- Peak Editor, typing `u23`: popup z 4400 over the dialog's 3300.
- Worst case for creation order: the Reference Photopeaks nuclide field used while docked, so its
  popup exists before the window; then View → Hide Tool Tabs → Tools → Reference Photopeaks. The popup
  was z 2200 over the window's 1100.

No `AuxWindow.cpp:162` change is needed.

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

**Status (2026-09-30): open, confirmed live.** A concrete case: Tools → Detector Response Select →
**Formula** tab. The three empty mandatory distance fields (placeholders "2.2 cm", "0 cm" and "1 m";
validator at `DrfSelect.cpp:2299-2302`, `setMandatory(true)`) are `Wt-invalid` (pink) in Wt 3 when the
tab is shown, and plain in Wt 4. No other Tools-menu window had an invalid field at open in either
version.

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

**Status (2026-09-30): open, confirmed live.** Test: open Quick MDA from a right-click near 303 keV with
Ba133 reference lines, then click 19 px in from the right edge of the Num FWHM spin box, which is text
area by the CSS.
- Desktop, 111 px field, `padding-right: 16px`: Wt 4 steps 4 → 5; Wt 3 leaves it at 4.
- iPhone 14 portrait (`?isphone=1`), 61 px field: `padding-right: 0` and no arrow background image, yet
  the same tap steps 4 → 5. So on the phone the right third of the field is an **invisible** stepper.
- An iPhone 14 landscape run did not step. That was not investigated.

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

**Status (2026-09-30): latent. The mechanism is confirmed, but nothing is clipped today.** Both stacks
carry Wt 4's inline `overflow:hidden` (Wt 3's are `visible`), but their content currently fits:
- **Export on phone** (iPhone 14 and iPhone SE, `?isphone=1`, Ba-133): the stack is 300 px. Its File,
  Format and Options tabs hold 283, 290 and 203 px, because the inner lists scroll themselves
  (`overflow-y:auto` at `ExportSpecFile.css:113`, `:119`).
- **Time-chart filters** (Passthrough sample, filter icon): the Int, Filt and Opt tabs are each
  223 px, with `scrollHeight == clientHeight`.

It would bite only if a tab's content grew past its box, for example a file with many detectors or
samples in the Export Options tab. That case was not tried.

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
Those two are still there (an AuxWindow, which has its own backstop — see below). Still present
2026-09-30, at the same lines; they are safe to delete.

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

**Resolved (re-measured 2026-09-30).** The inline pixel `top` is still written, but it no longer
matters. `SimpleDialog.css:24-28` overrides it with `top: 50% !important;
transform: translateY(-50%) !important` on `.Wt-dialog.simple-dialog`, so the dialog is centred by CSS
and follows the window. Measured on Edit → Enter URL, shrinking the window from 837 px to 417 px:
- the inline `top` stayed `305px`;
- the rendered top of the 228 px dialog went from 305 px to 95 px, i.e. (417 − 228)/2, so it stayed
  centred with the footer on screen.

Wt 3, for comparison, re-wrote its inline `top` to 87 px.

#### Related, fixed 2026-09-29: absolutely positioned content lost its dialog-body anchor

Wt3's JS layout set `position:absolute` on every layout item (`StdGridLayoutImpl2.js:432-433`), so the
dialog body was the containing block for anything absolutely positioned inside it. Flex leaves the
body static, so such content now anchors to the whole `.Wt-dialog`. Symptom: the Peak Editor's
prev/next-peak arrows (`bottom:0`) sat on the footer's help/Delete buttons. Fixed with
`position:relative` on `.PeakEdit`. A live scan (every Tools/View/Help item, and every docked tool
opened as a window with tool tabs hidden, desktop and phone) found no other element whose
`offsetParent` is the dialog itself.

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

**Status (2026-09-30): the code shape is unchanged in all three, and only Gamma XS shows anything.**
The three hosting sites are still `RemoteRid.cpp:2402`, `DecayActivityDiv.cpp:1844` and
`GammaXsGui.cpp:907`. Heights were sampled every 120–150 ms after opening (desktop, 1400×850):
- **External RID** (after Continue on the warning): 587×345, constant from the first sample; content
  reachable, no jitter seen.
- **Add Nuclide:** 216×543, constant from the first sample (see Issue 20).
- **Gamma XS:** 422×521, already stable by 300 ms, but it opens with its top at **38 px** rather than
  centred (Wt 3: 163 px, centred). `GammaXsWindow` calls `resizeToFitOnScreen()`,
  `centerWindowHeavyHanded()` and `centerWindow()` exactly as Wt 3 does (`GammaXsGui.cpp:960-963`). So it
  was most likely centred against a taller interim height from the multi-pass layout: 38 px is where an
  ~774 px window would be centred. The extra width is separate and deliberate:
  `setMinimumSize( 420, Auto )` at `GammaXsGui.cpp:903` works around `WGridLayout`'s absolute
  positioning.

Phone was not re-checked for these three.

### Issue 31 — `WMenu::select()` now auto-selects the parent item

**Mechanism:** Wt4 `WMenu.C:317-322` selects the owning parent item before selecting the child; no
Wt3 equivalent (`WMenu.C:328`). `select()` emits `triggered()` when the index changes
(`WMenu.C:348-360`). Every InterSpec submenu sets `m_parentItem` (`PopupDiv.cpp:1326`).

**Symptom:** picking an entry inside a submenu also fires `triggered()` on the parent item that owns
it, and marks the parent selected. Trigger: right-click a peak → "Change Nuclide" → pick a nuclide.

**Fix:** `parentItem->setSelectable(false)` in `PopupDivMenu::addPopupMenuItem()` — the new Wt4 code
explicitly honours that flag.

**Status (2026-09-30): code unchanged, and no observable symptom.** `m_parentItem` is still left
selectable (`PopupDiv.cpp:1435-1439`). Live test on a Ba-133 356 keV peak, right-clicking and then
picking from a submenu:
- Change Nuclide → "U235 356.03 keV": the peak's nuclide became U235.
- Change Skew Type → GaussExp: the skew changed.

In both cases the menu closed and the parent `<li>` was left with only `item submenu` (not selected).

The fix is still worth making cheaply, because a real hazard sits next to it.
`InterSpec::rightClickMenuClosed()` (`InterSpec.cpp:1974`, on `aboutToHide()`) resets
`m_rightClickEnergy`, and the Change Continuum / Change Skew handlers read it (`:2117`, `:2126`).
If the parent-first `done()` path described above ever ran before the child's handler, those two items
would silently do nothing. Change Nuclide is immune, because it captures the peak.

### Issue 32 — Smaller confirmed deltas

- **`GammaCountDialog` width lost.** Wt4's default theme adds
  `.Wt-dialog .Wt-fill-width .body { width:100% }` (`wt.css:300-303`), specificity (0,3,0), which
  outranks `GammaCountDialog.css:1-9` (0,2,0). Tools → Energy Range Sum. Fix: `contents()->setWidth(285)`.
  **Still open, 2026-09-30.** The window is 407 px wide with a 395 px body, against Wt 3's 297 / 285 px.
  The intro sentence and the "Shift-⌘-Drag" hint now sit on one line each instead of wrapping centred.
- **`.Wt-tooltip { max-width: 280px }` is dead** — Wt4 `js/ToolTip.js:161` writes `maxWidth` inline
  (Wt3 never did). Long help tooltips no longer wrap at 280px. 313 `attachToolTipOn` call sites.
  Since Wt hardcodes it in JS, `!important` is justified here; alternatively cap the inner div.
  **Fixed** in 36726c7d (2026-08-22). `InterSpec.css:1634` now has
  `max-width: min( 360px, calc(100vw - 40px) ) !important`, plus a 760 px variant for
  `.Wt-tooltip:has(table)`.
- **Submenus clamped to a scrolling parent.** Wt4 `Wt.js:1827-1837` makes `fitToWindow` use the
  parent's visible rectangle when the parent is not BODY/domRoot. `.PopupDivMenu` is exactly such a
  scrolling ancestor (`InterSpec.css:1076-1088`). Only `.PopupDivMenu.AppMenu` is protected
  (`InterSpec.css:1094-1099`); plain context menus are not.
  **Not observed, 2026-09-30.** Tested with a peak right-click → Change Nuclide in a 1400×400 window.
  The 226 px context menu did not itself scroll. The submenu was parented to `Wt-domRoot` rather than
  re-parented, sat at y = 1 with 362 px height, and scrolled internally (448 px of content). The clamp
  needs the context menu itself to be scrolling, which means a window shorter than about 230 px. The
  AppMenu mitigation is now `.PopupDivMenu.AppMenu:has(.submenu) { overflow-y: hidden }`
  (`InterSpec.css:1318`).
- **`Wt-hide-scrollbar` relies on `scrollbar-width:none`** (Wt4 `WTableView.C:94`), unsupported by
  Safari before 18.2 — and per CLAUDE.md the macOS target is Safari 15.4. A stray horizontal
  scrollbar can appear in table-view headers on the native macOS build. Fix:
  `.Wt-hide-scrollbar::-webkit-scrollbar { display:none; }`.
  **Still open, 2026-09-30.** There is no such rule in `InterSpec_resources/`. This was not testable
  here, since it needs Safari earlier than 18.2 / macOS 12 WKWebView.
- **`WPanel` collapse icon** is now a CSS-background `Wt-collapse-button` class instead of two
  `<img>`s. Nothing breaks, but dark theme should now add
  `.Wt-collapse-button { filter: invert(100%); }` (Make DRF, Isotopics from peaks, LLM tool).
  **Still open, 2026-09-30.** `InterSpec_resources/` has no rule for it (the 2026-09-21 token rework,
  b3539981, did not add one). This was not checked visually in the Dark theme.
- **`WText::setPadding` now emits the `padding` shorthand** with unset sides as `0` (Wt3 emitted
  longhands and discarded top/bottom). So `setPadding(x, Side::Left)` on a **WText** zeroes any
  stylesheet-supplied top/right/bottom. Three sites, none currently killing a live rule:
  `EnergyCalAddActions.cpp:1376`, `PeakSearchGuiUtils.cpp:653`, `LlmSubAgentFollowup.cpp:222`
  (line numbers as of 2026-09-30).
  `WContainerWidget` is unaffected (no delta). Rule: never mix `WText::setPadding(single side)` with
  CSS padding on the same WText. Stale comment to update: `QrCode.cpp:463` (still stale 2026-09-30).
- **`WAnchor` with a null link** now emits `href="#"` + `Wt-no-default` instead of omitting `href`,
  so such anchors get the pointer cursor, focus ring and a tab stop. Makes the comment at
  `LlmToolGui.cpp:346-347` ("Empty WLink => an <a> with no href") wrong; still stale 2026-09-30.
  Forward risk: a null-link anchor with no `clicked()` handler would navigate to `#` and clobber
  `window.location.hash`, which InterSpec uses for its internal path.

### Dead-by-attrition CSS (safe to delete; classes emitted by neither Wt tree)

All still present 2026-09-30, at these line numbers: `InterSpec.css:1660` (`br.Wt-tabs-clear`),
`LlmToolGui.css:318` (`.LlmThinkingContent > .Wt-text`), and `InterSpec.css:1888`
(`.Wt-suggest.nuclide-suggest .Wt-suggest-item`, the second half of a two-selector rule). Also dead
`#include`s: `CompactFileManager.cpp:36` and `DecayActivityDiv.cpp:44` (`WSlider`, never
instantiated). The orphan `SimpleMdaBody` rules (`DetectionLimitSimple.css:305`, `:313`, Issue 29)
belong on this list too.

### Non-UI items outstanding (tracked here so they are not lost)

Lifetime/crash bugs, not rendering. The original write-up, `/tmp/wt4_followup_prompt.md`, no longer
exists. Re-checked against the source on 2026-09-30 by reading code only; none of this was reproduced.
1. Cat-A tool windows outlive their owner on **File → Clear Session…** and then dereference a freed
   tool (`~EnergyCalTool` is empty; `~InterSpec` misses `m_llmToolWindow`, `m_terminalWindow` and the
   EnergyCal preserve window; `DecayChainChart` has no destructor at all).
   **Fixed** (7760358d, 2026-08-14).
   - `~EnergyCalTool` (`EnergyCalTool.cpp:1418-1434`) now tears down its More-Actions and graphical
     recal windows.
   - `~InterSpec` closes the LLM window (`:1394-1397`), the preserve window (`:1413`) and the terminal
     (`:1428-1430`), and then runs the `trackToolDialog` sweep (`:1446-1477`).
   - `~DecayChainChart` exists (`DecayChainChart.cpp:315-318`).
2. Decay-tool CSV export: a dropped `wApp->bind()` guard gives a cross-thread use-after-free if the
   dialog is closed mid-stream (`DecayActivityDiv.cpp:1258`, posted at `:1447`). Five more
   dropped-guard sites listed there.
   **The CSV export is fixed** (7760358d). `DecayCsvResource` takes the update lock
   (`DecayActivityDiv.cpp:1264`) and reaches the widget only through an id-based `WidgetHandle`
   (`:1279-1283`, `:1333`, `:1484`).
   The remaining `InterSpec.cpp` sites found on 2026-09-30 were **fixed the same day** (uncommitted
   when written):
   - **External RID undo step.** It captured the self-deleting `warning` `SimpleDialog`
     (`[warning](){ warning->accept(); }`). Reproduced on the Aug-20 Release build: Tools → External
     RID → Cancel → Edit → Undo gave `SIGSEGV` in `createRemoteRidWindow()::$_3` and killed the
     **whole server**. Undoing while the warning was still up closed it without re-enabling the menu
     item, so External RID could not be reopened for the rest of the session. Now the warning is
     tracked in `m_remoteRidWarning` (an `observing_ptr`, also torn down in `~InterSpec`). The undo
     looks it up there, skips a warning that has already been answered, and otherwise re-enables the
     menu item and `reject()`s it.
   - **Right-click nuclide-suggestion `updater`**, the **DRF-load worker and its post-back**, the
     **External RID auto-call on load**, and the **hint-peak** and **background-peak-recovery**
     threads: each now names the viewer with a `WidgetUtils::WidgetHandle` and resolves it on the
     session thread, instead of capturing `this`.
     - The DRF loader (`loadDetectorResponseFunction`, now `static`) receives the DB session, user id
       and GADRAS search paths as values. To support that, `DrfSelect::getUserPreferredDetector` takes
       a user id instead of a `Dbo::ptr`, and `DrfSelect::initAGadrasDetector` takes the search-path
       string, which the new `DrfSelect::gadrasDrfSearchPaths( InterSpec * )` builds.
     - `doFinishupSetSpectrumWork` is `static`.
     - The parse-warnings worker now posts its messages instead of taking
       `WApplication::UpdateLock` on a raw `WApplication *`.
   - Verified on the fixed build by driving the app: the three undo scenarios (Cancel then Undo; Undo
     while showing; Continue, Undo, Undo, Redo) pass; the default DRF still loads; the Change Nuclide
     candidates still fill in; load-then-Clear-Session ×3 is clean.
   - Still open, low risk: `:5218` and `:12321` post `[this]` from the session thread to run on it
     shortly after.
3. Right-click popup menus leak on mobile — `PopupDivMenu::isHidden()` returns true on mobile, so
   Wt4's `WPopupMenu::done()` returns before emitting `aboutToHide()` and the deferred cleanup never
   runs.
   **Partly fixed.** The `isHidden()` override is gone, replaced by `setHideOnSelect(false)`
   (`PopupDiv.cpp:694`; explained at `PopupDiv.h:236-243`). But a phone item's `triggered()` →
   `mobileHideMenuAndParents` still hides the menu before `WMenu::select` reaches `done()`, and
   `hideOnSelect=false` suppresses `aboutToHide_` anyway, so `aboutToHide()` still does not fire on a
   phone selection. The leak is now bounded: `D3SpectrumDisplayDiv.cpp:3910` and
   `RelActAutoGui.cpp:2320` drop the previous menu on the next open. `InterSpec::rightClickMenuClosed`
   still does not run after a phone pick (a minor issue: stale `m_rightClickEnergy`).
4. Dose Calc runs its handler twice per click; "Add peak…" can leave an unclosable modal on a
   degenerate spectrum.
   **Both fixed** (7760358d).
   - `right_select_item` (`DoseCalcWidget.cpp:88-101`) no longer emits `triggered()` itself.
   - The degenerate-spectrum path of `AddNewPeakDialog` (`:97-122`) now has a Close button and
     `rejectWhenEscapePressed()`.

### Already fixed on 2026-08-13/14 (for reference, do not re-report)

`WStackedWidget` inline-overflow in BatchGuiWidget + DrfSelect + RelActAutoGui (×2) + RemoteRid +
EnergyCalPreserveWindow; Multi-File Calibration dialog clipping/jitter; Peak Manager continuum-type
editor throwing `bad_any_cast`; the `delete`-on-parent-owned-widget family; `observing_ptr` captured
across thread boundaries; deferred widget destruction happening on an io-service thread; the
"Keep Previous Calibration?" **No** button doing nothing (it was wired straight to
`AuxWindow::emitReject`, which suppresses its deferred emit unless the dialog is already hidden).


## Issues found 2026-08-15 while verifying the upgrade/Wt4 → upgrade/DetEff merge

All of these are **pre-existing** — each was checked against the pre-merge tips and is not merge
damage. (As of 2026-09-30, three are fixed and one is not a regression; see each **Status** line.) Found by driving the merged app (desktop + phone) and by a DOM/layout audit
(offscreen windows, non-scrolling overflow, children escaping clipping parents).

### Phone mode tears down the session on load — `?isphone=1` is unusable

**No longer reproduces (checked 2026-09-29):** `?isphone=1` sessions (iPhone 14, iPhone SE, Galaxy S8 and
iPhone 14 landscape via Playwright) load, open the hamburger menu, load a sample, fit peaks and open
tools without the session being torn down. Kept below for the history of the mechanism.
The 2026-09-30 phone runs (iPhone 14, iPhone 14 landscape, iPhone SE) loaded and worked too. The
likely fix is f3a73c99 (2026-09-22, "Stop the phone layout dying on load with a spectrum-chart script
error"), which touched only `D3SpectrumDisplayDiv.cpp`. The undo-menu tooltip wiring described below
is unchanged and unguarded (`InterSpec.cpp:6885`), so either the tooltip diagnosis was wrong or that
stale-id path is still latent.

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

**Status (2026-09-30): fixed.** 3f1e2e71 (2026-08-25) moved `ul.HorizontalMenu` into `InterSpec.css`,
and the dialog now uses a `SideMenu` on desktop and tablet. The tabs (General / Geom & MC / FWHM /
Anchor) are a proper left-hand nav column of 142 px buttons, not full-width stacked bars. It was
checked at 1400×850 and on an iPad Mini (`?istablet=1`). Help → About InterSpec also renders its
Disclaimer/License/Credits/Contact/Data items as a correct side menu.

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

**Status (2026-09-30): fixed by the redesign.** 10191318 (2026-09-13) rebuilt Make Detector Response
around detector geometry. "Peak Fit Prefs" is now its own group box beside "Options", below the
geometry panel, and its Det. Type / FWHM Method / Skew Type combos sit inside its border with no
overlap onto the chart.

`.MakeDrfOptions { width: 150px; }` (`MakeDrf.css:6`) is a hard width, but the `PeakFitDetPrefsGui`
block inside it needs ~196 px, so the Det. Type / FWHM Method / Skew Type combos stick ~48 px past the
panel's right border and overlap the chart. Present on both pre-merge tips; `MakeDrf.css` was not
touched by either branch.

### "Create Detector Response Function" opens far taller than its content

**Status (2026-09-30): not a regression. Wt 3 does the same.** This is Tools → Make Detector Response.
In an 850 px window both versions open it at 902×809 px, with the same large empty area below the
chart. That area is where the per-file peak table goes, and it is empty until peaks are fitted. It is
by design and needs no Wt 4 fix.

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

**Status (2026-09-30): deliberate and pinned by tests.** The code is unchanged
(`DetectorPeakResponse.cpp:778`). The "no emptiness test" point makes no practical difference, because
the parser only sets `m_totalEfficiency` when there are at least two points and at least one is
non-zero. The hash change itself is intended:
- `test_DetectorPeakResponse.cpp:2797` calls it "the accepted lineage break";
- `test_GadrasDetectorDat` (`test_shipped_gadras_drfs_unchanged`) pins the current stock hashes.

The inferred effect is one duplicate "Previous" DRF row per detector after upgrading. `UseDrfPref`
stores the row id, so default-DRF associations should survive. The misleading comment at
`DetectorPeakResponse.cpp:773-774` could still be corrected.
