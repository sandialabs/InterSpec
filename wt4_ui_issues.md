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

**Dead rules** (all target `WDialog::contents()` == `.body` == `.AuxWindow-content`):
`ExportSpecFile.css:10`, `:34`, and knock-on `:39-45`; `SimpleDialog.css:62`, `:78`, `:87`;
`DetectionLimitSimple.css:305-308`, `:313-315`; `BatchGuiWidget.css:4`; `InjaLogDialog.css:9`;
`RefSpectraWidget.css:6`; `ShieldingSourceDisplay.css:8`; `FitPeaksForNuclidesGui.css:13`;
`DrfSelect.css:221`.

**Dead C++, not just CSS:** `SimpleDialog.cpp:236-245` builds a per-instance `#<id> .body`
max-width/max-height rule; its comment claims it "wins without `!important`", which is true against
other stylesheets but false against an inline style. So `SimpleDialog::setMaximumSize()` is a no-op
on the scrollable body.

**Fix:** set these from C++ (`contents()->setHeight(...)` / `setMaximumSize(...)`), which
`FlexLayoutImpl` transfers to the wrapper's flex-basis. A **global escape hatch** also exists:
`WLayout::setDefaultImplementation( LayoutImplementation::JavaScript )` once at startup restores the
Wt3 engine app-wide and retires this whole issue — at the cost of flex's auto-sizing. That is a
judgement call, not an obvious win.

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
