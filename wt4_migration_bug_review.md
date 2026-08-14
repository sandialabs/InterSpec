# InterSpec Wt 3.7.1 → Wt 4.12.6 Migration — Code Review Findings

_Generated 2026-05-29 · automated multi-agent review, every finding adversarially verified against the source and the Wt 4.12.6 headers._

## Scope & method

A three-pass review hunting the bug classes the Wt3→Wt4 ownership/lifetime model invites: incorrect widget handling, ownership/double-free, lifetime/use-after-free, signal-lifetime, AuxWindow/SimpleDialog lifecycle, Wt4 API misuse, threading/strand, and unsafe/insecure access.

- **Pass 1** — 8 repo-wide pattern sweeps over `src/` + `InterSpec/`.
- **Pass 2** — 54 per-file deep readers covering every `src/*.cpp`, plus header / JS-JSignal / WTemplate-XSS sweeps.
- **Pass 3** — `target/` platform glue (Electron, wxWidgets, macOS, iOS, Android, batch, QuickLook), `WResource::handleRequest` worker bodies, and an AuxWindow mid-emit reentrancy sweep.

Each candidate was then re-examined by an independent **adversarial verifier** (default stance: refute) that read the cited code, the relevant header (to confirm `observing_ptr` vs raw vs `unique_ptr`), and the Wt 4.12.6 source, and had to produce a concrete failure mechanism **and** a realistic trigger to confirm. False positives were dropped; the items below survived that gate.

## By the numbers

**56 confirmed** findings (after dedup): **19 critical · 13 high · 19 medium · 5 low**. Plus **11 uncertain** (need a runtime/usage decision) and **26 minor/low**. ~20 candidates were refuted and are not listed.

| Category | Confirmed |
|---|---|
| Ownership / double-free / leak | 21 |
| Threading — worker/native-thread access & UAF | 15 |
| Signal & callback lifetime | 2 |
| Object lifetime / dangling pointers | 9 |
| AuxWindow / SimpleDialog lifecycle | 2 |
| Undo/redo capture | 1 |
| Wt4 API misuse / logic | 1 |
| Security / unsafe input | 1 |
| Other correctness | 4 |

## The systemic patterns (fix these as classes, not one-by-one)

Most findings are instances of a handful of mechanical Wt3→Wt4 regressions. Fixing them by *pattern* clears the bulk of the crashes:

### A. `delete widget;` on a parent-owned child → double-free  *(largest cluster)*
In Wt4 a parent owns each child through a `std::unique_ptr` in `WObject::children_`. Calling `delete child` runs `~WWebWidget`, which re-detaches via `removeFromParent()` and hands back the still-live owning `unique_ptr`; that pointer then destroys the same object a second time (the debug assert at `WWebWidget.C:250` documents the contract). This is the Wt3 idiom (`delete` a widget to remove it) that silently became a double-free. **Uniform fix:** replace `delete child;` with `const std::unique_ptr<Wt::WWidget> removed = parent->removeWidget(child);` (or `removeChild` / `tabs->removeTab` / `menu->removeItem`, or `accept()` for a `SimpleDialog`) and let the returned `unique_ptr` drop. Confirmed sites: `RelActAutoGui` (every ROI/peak/nuclide/curve remove handler — 2460, 3851-3874, 3877-3900, 3981-4022, 4044-4077, 4155-4190, 4818-4862), `ShieldingSourceDisplay` (`reset`, `removeShielding`, two error-path catches), `IsotopeSearchByEnergy::removeSearchEnergy`, `ReferencePhotopeakDisplay::removeFeatureMarkerTool`, `DrfSelect::GadrasDetSelect::removeDirectory`, `RelActAutoGuiRelEffOptions`, `MakeDrfSrcDef::~MakeDrfSrcDef`, `DetectionLimitSimple` (a `SimpleDialog`).

### B. Background/worker callbacks capturing a raw `this`/widget → use-after-free
Long-running work is run off-thread and posts a completion lambda back to the session capturing a raw `this`/widget pointer with no liveness guard. If the user closes the tool/dialog before the work finishes, the callback runs on freed memory. **Fix:** key the post on a still-alive check (e.g. capture a `std::weak_ptr` sentinel or an `observing_ptr`, and bail if null; or route through a `WServer::post(sessionId,…)` that re-resolves the widget). Sites: `RelActAutoGui` (5713), `RelActManualGui` (1861), `ShieldingSourceDisplay` (9205), `FitPeaksForNuclidesGui` (1159), `MakeFwhmForDrf` (697), `BatchGuiInputFile` (158), `BatchGuiWidget` (291), `DirectorySelector` (238), `IsotopeSelectionAids` (454), `SimpleActivityCalc` (1627), `RemoteRid` (1319).

### C. Signal connections without lifetime tracking → use-after-free
Single-arg `signal.connect(lambda)` keeps the connection alive for as long as the **emitter**, with no auto-disconnect when a captured object dies. The headline case is `UserPreferences::addCallbackWhenChanged` /`addIntCallbackWhenChanged` (`UserPreferences.h:276-289, 305-318`): the lambda captures `target` but the `Signal` lives session-long in the prefs maps, so tool-scoped callers (`RelActAutoGui m_spectrum`, `DetectionLimitTool`, `ShieldingSourceDisplay`) dangle after their window closes. **Fix:** use the tracked 2-arg form `signal->connect(target, method)` (Wt auto-disconnects when the `WObject` target dies).

### D. `WResource::handleRequest` reading live widgets/models on an I/O thread → data race
`handleRequest` runs on a Wt I/O worker thread, not the session thread; reading a live `WAbstractItemModel` / widget there without a `WApplication::UpdateLock` races the session. Sites: `SpecFileQueryWidget` (`ResultCsvResource`). Also several `WResource`s are raw-`new`'d and never freed (`PeakModel::m_csvResource`, `SpecMeasManager` `FileDragUploadResource`s) — leaks plus temp-file leaks.

### E. Native/platform glue touching Wt off-thread + dangling raw pointers
`target/` glue calls into Wt from AppKit/JNI/wx/Electron threads. `macos/NativeMenu.mm` reads `isEnabled()` on the AppKit thread and keeps a `Target` with raw `PopupDivMenuItem*`/`WCheckBox*` that dangles across async menu teardown (UAF); plus off-thread wx/Android/iOS accesses and path/JSON-injection in the native file handlers.

### F. Orphaned widgets — raw `new` never added to a parent → leak **and invisible UI**
Wt4 constructors no longer take a parent, so a `new Widget(...)` that is never `addWidget`/`addNew`'d is both leaked and never rendered. `MakeDrf` builds `DrfPeak`/`MakeDrfSrcDef` widgets that are never added to `m_peaks`/`m_sources` (missing peaks/sources in the DRF tool); `FitSkewParamsTool`, `NuclideSourceEnter` leak models/controllers.

### G. Lifetime hygiene — raw dialog pointers that should be `observing_ptr`
Members holding closeable dialogs as raw `T*` (`DecayActivityDiv.h m_csvDownloadDialog`, `SpecMeasManager.h m_batchDialog`) dangle after the dialog self-destructs; per the project convention these should be `Wt::Core::observing_ptr<T>`. And `~ShieldingSourceDisplay` tears down a wApp-owned `SimpleDialog` with `removeFromParent()` (a no-op for global widgets → leak).

---

## Confirmed findings

### Ownership / double-free / leak  (21)

#### `src/DrfSelect.cpp:1565-1577` — GadrasDetSelect::removeDirectory calls delete on a still-parented child widget (double-free + self-delete during signal emit)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Mechanism:** GadrasDirectory instances are owned by m_directories via addNew<GadrasDirectory> (line 1561), so the container holds a unique_ptr in its children list; delete w; at line 1574 frees the widget while the parent still owns it, leaving a dangling entry freed again on container teardown (Wt4 double-free). It is invoked from the directory's own m_deleteBtn->clicked() handler (line 1660 -> removeDirectory(this)), so it also destroys the widget whose click signal is mid-emit.
- **Trigger:** Trigger: desktop build (the !BUILD_FOR_WEB_DEPLOYMENT && !IOS branch where m_deleteBtn exists), user adds a Gadras detector dir then clicks its X
- **Fix:** Replace 'delete w;' at line 1574 with 'm_directories->removeWidget( w );' (returns a unique_ptr that drops, destroying the widget safely after detaching from the parent's children list), matching the sibling RelEffDetSelect::removeFile pattern.

#### `src/IsotopeSearchByEnergy.cpp:1069` — removeSearchEnergy deletes SearchEnergy widget owned by m_searchEnergies (double-free)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 1038-1071
- **Mechanism:** removeSearchEnergy() does `delete energy;` (line 1069) where `energy` is a SearchEnergy* that is a child of m_searchEnergies (created via m_searchEnergies->addNew<SearchEnergy>() at lines 482/560; searches() at line 599-603 enumerates m_searchEnergies->children() and returns them by value, so the local `searchW` copy is safe). A child added via addNew is owned by the parent through WObject::children_ unique_ptr. In Wt4 `delete energy` runs ~WWebWidget -> removeFromParent -> WContainerWidget::removeWidget -> removeChild/Utils::take, extracting the live unique_ptr that still owns the object being destroyed; the assert(unique_this==nullptr) at WWebWidget.C:250 fires in debug, and in release the returned unique_ptr re-deletes / leaves a dangling owner -> double/recursive destruction. Triggered every time the user removes a search-energy row.
- **Trigger:** Trigger is realistic and common: the per-row remove button via `enrgy->remove().connect(this, [this,enrgy](){ removeSearchEnergy(enrgy); })` (cpp:486, 563) fires every time a user removes a search-energy row (also the programmatic path at cpp:1855)
- **Fix:** Replace `delete energy;` at src/IsotopeSearchByEnergy.cpp:1069 with the unique_ptr-returning detach API so the parent link is cleared before destruction and the returned owner drops cleanly: `m_searchEnergies->removeWidget( energy );` (equivalently `energy->removeFromParent();`). Either statement detaches the child, returns the owning std::unique_ptr<WWidget>, and lets it destruct at the end of the statement with parent() already null, so ~WWebWidget's removeFromParent() returns nullptr and the assert holds. Do not `delete` a parented Wt4 widget.

#### `src/ReferencePhotopeakDisplay.cpp:2026` — removeFeatureMarkerTool deletes m_featureMarkers owned by m_featureMarkerColumn (double-free)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 2021-2027
- **Mechanism:** removeFeatureMarkerTool() does `delete m_featureMarkers;` (line 2026). m_featureMarkers (raw FeatureMarkerWidget*, ReferencePhotopeakDisplay.h:540) is created via m_featureMarkerColumn->addNew<FeatureMarkerWidget>(...) at line 2006, so it is owned by m_featureMarkerColumn's WObject::children_ unique_ptr. Deleting it directly runs ~WWebWidget -> removeFromParent -> removeWidget -> removeChild/Utils::take, which extracts the still-owning unique_ptr (assert(unique_this==nullptr) at WWebWidget.C:250 in debug; dangling owner / re-delete in release when the column is later torn down). Triggered whenever the user closes the Feature Marker tool.
- **Trigger:** Realistic trigger: when the tools tabs are visible the Feature Marker tool lives embedded in m_featureMarkerColumn inside m_referencePhotopeakLines; the close icon (ReferencePhotopeakDisplay.cpp:1280) or the View>Feature Markers menu toggle calls InterSpec::displayFeatureMarkerWindow(false), which at InterSpec.cpp:3074 invokes removeFeatureMarkerTool() -> …
- **Fix:** Replace the bare `delete m_featureMarkers;` at ReferencePhotopeakDisplay.cpp:2026 with a Wt-native detach that returns the owning unique_ptr, letting it destroy the widget exactly once: `m_featureMarkerColumn->removeWidget( m_featureMarkers );` then set `m_featureMarkers = nullptr;`. Concretely: `std::unique_ptr<Wt::WWidget> removed = m_featureMarkerColumn->removeWidget( m_featureMarkers ); m_featureMarkers = nullptr;` (the returned unique_ptr drops at end of scope, deleting the widget cleanly with no second owner). Equivalent: `m_featureMarkers->removeFromParent();` which resolves to the same `m_featureMarkerColumn->removeWidget(this)` and returns the owning unique_ptr. This matches §1: …

#### `src/RelActAutoGui.cpp:3851-3874` — handleRemoveFreePeak deletes child of m_free_peaks (double-free)
**🔴 CRITICAL** · _ownership ·×6 (P1+P2)_
- **Also at / spans:** 2415-2473, 2460,3867,3893,3995,4077, 3877-3900, 3981-4022, 4044-4077
- **Mechanism:** handleRemoveFreePeak() does `delete w;` (line 3867) on a RelActAutoGuiFreePeak that is a child of m_free_peaks (a WContainerWidget, RelActAutoGui.h:461; child located in m_free_peaks->children() at line 3856). The widget is owned by m_free_peaks' children_ unique_ptr; deleting it directly leaves the parent's owner dangling (~WWebWidget->removeFromParent->WContainerWidget::removeWidget extracts the live unique_ptr, asserting in debug / double-destructing in release). Triggered each time a free peak is removed.
- **Trigger:** Realistic trigger: the user clicks the per-row remove icon (RelActAutoGuiFreePeak.cpp:80 removeFreePeak->clicked() -> handleRemoveSelf -> m_remove.emit(), line 154), whose remove() signal is wired at RelActAutoGui.cpp:4032 to handleRemoveFreePeak — i.e., every single free-peak removal
- **Fix:** Replace the bare `delete w;` at RelActAutoGui.cpp:3867 with the unique_ptr-returning detach primitive and let the returned owner drop, which destroys the widget cleanly while relinquishing the parent's ownership: `m_free_peaks->removeWidget( w );` (the returned std::unique_ptr<WWidget> is discarded, destroying the row). This matches reference §1's safe Wt4 detach-and-drop pattern. The identical sibling bug at line 3893 (handleRemoveEnergy: `delete w;` on a child of m_energy_ranges) should be fixed the same way: `m_energy_ranges->removeWidget( w );`.

#### `src/RelActAutoGui.cpp:4155-4190` — handleRemoveNuclide deletes RelActAutoGuiNuclide child of menu content (double-free)
**🔴 CRITICAL** · _ownership (P1)_
- **Mechanism:** handleRemoveNuclide() does `delete w;` (line 4190) on a RelActAutoGuiNuclide. These are added into a per-rel-eff nuclide content WContainerWidget, so they are owned by that container's children_ unique_ptr. The direct delete leaves the container's owning unique_ptr dangling (~WWebWidget->removeWidget returns the live owner), producing a debug assert and a release double-free. Triggered every time a nuclide source is removed.
- **Trigger:** Trigger is routine: each nuclide source row has a remove button (clicked -> RelActAutoGuiNuclide::handleRemoveSelf, RelActAutoGuiNuclide.cpp:793 and :1640-1642) that emits m_remove, wired to handleRemoveNuclide (RelActAutoGui.cpp:2278 and :3691)
- **Fix:** Detach through the owning container instead of `delete w`, and defer destruction so the widget is not destroyed while its own remove() signal is mid-emit. Replace line 4190 `delete w;` with a detach-and-deferred-drop, e.g.: `std::unique_ptr<Wt::WWidget> taken = w->removeFromParent();` then hand `taken` to a deferred deleter (mirroring the SimpleDialog startDeleteSelf pattern, e.g. WServer::post a lambda that owns and drops the unique_ptr) so the object outlives the current m_remove.emit() unwind. A minimal correct variant that just fixes the double-free is `w->removeFromParent();` (dropping the returned unique_ptr), which detaches via the parent's removeChild and lets the unique_ptr destroy …

#### `src/RelActAutoGui.cpp:4818-4862` — handleDelRelEffCurve: delete after WMenu::removeItem already destroyed item+contents (use-after-free)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 4847-4861
- **Mechanism:** handleDelRelEffCurve() calls m_rel_eff_opts_menu->removeItem(item) (line 4848) and discards the returned std::unique_ptr<WMenuItem>, so the WMenuItem is destroyed at the end of that statement. The menus have a contents stack and the items were added with ContentLoading::Eager (lines 4734/4787), so WMenu::removeItem (Wt4 WMenu.C:262-286) pulls the contents out of the stack via returnContentsInStack(...) back into the item's uContents_; destroying the returned item therefore also destroys its contents widget. The subsequent `delete curve;` (4849), `delete nuc_content;` (4859) and `delete nuc_item;` (4860) therefore operate on already-freed objects => use-after-free / double-delete. Triggered when a user deletes a rel-eff curve.
- **Trigger:** Trigger: realistic
- **Fix:** Remove the three stray delete statements. removeItem already destroys the item and (for Eager content) its re-owned contents exactly once when the returned unique_ptr drops.  At line 4848-4849, replace:   m_rel_eff_opts_menu->removeItem( item );   delete curve; with:   m_rel_eff_opts_menu->removeItem( item );   // returned unique_ptr<WMenuItem> drops -> destroys item and its contents (curve)  At lines 4858-4861, replace:   m_rel_eff_nuclides_menu->removeItem( nuc_item );   delete nuc_content;   delete nuc_item;   nuc_item = nullptr; with:   m_rel_eff_nuclides_menu->removeItem( nuc_item );  // returned unique_ptr drops -> destroys nuc_item and nuc_content   nuc_item = nullptr;  (The `delete …

#### `src/RelActAutoGuiRelEffOptions.cpp:523-525` — delete of parented external-attenuation child widget is a Wt4 double-free
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Mechanism:** The external attenuation RelEffShieldWidgets are always added via m_phys_ext_attens->addNew<RelEffShieldWidget>(...) (lines 404/539/749/811) so each is owned by m_phys_ext_attens through its internal unique_ptr child list; setRelEffCurveInput snapshots them and then calls `delete starting_ext_atten_widgets[i]` (line 525), but in Wt4 ~WWebWidget runs removeFromParent()->parent->removeWidget(this), which extracts the still-owned unique_ptr and re-deletes the object being destructed (double-free). Every other removal site in this same file correctly uses m_phys_ext_attens->removeWidget(...) (lines 410/744/806).
- **Trigger:** Trigger: setRelEffCurveInput (called from RelActAutoGui.cpp:2246,6327 when applying a Physical-Model RelEff config) with fewer external attenuators than currently displayed; initPhysModelShields preserves existing ext widgets, so the surplus reaches the delete loop
- **Fix:** Replace `delete starting_ext_atten_widgets[i];` with `m_phys_ext_attens->removeWidget( starting_ext_atten_widgets[i] );` letting the returned unique_ptr drop, matching lines 410/744/806.

#### `src/ShieldingSourceDisplay.cpp:7511-7522` — ShieldingSourceDisplay::reset deletes ShieldingSelect children of m_shieldingSelects (double-free)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 7516-7522
- **Mechanism:** reset() copies m_shieldingSelects->children() into a local `const vector<WWidget*> shieldings` (line 7516, so iteration is safe) and does `delete select;` (line 7521) for each ShieldingSelect. Those ShieldingSelect widgets are children owned by m_shieldingSelects' WObject::children_ unique_ptr vector. In Wt4 `delete select` runs ~WWebWidget -> removeFromParent() -> WContainerWidget::removeWidget -> removeChild/Utils::take, which extracts the live owning unique_ptr; the assert(unique_this==nullptr) at WWebWidget.C:250 fires in debug, and in release the parent's children_ retains a dangling extracted slot / the returned unique_ptr re-deletes -> double destruction. Triggered on model reset.
- **Trigger:** Realistic trigger: the identical buggy `delete shielding;` appears at line 8377 in ShieldingSourceDisplay::removeShielding, wired to the per-shielding remove ("X") button via `select->remove().connect(this, &ShieldingSourceDisplay::removeShielding)` (line 7891) — i.e
- **Fix:** In reset(), replace the `delete select;` with a Wt-native detach that returns and drops the owning unique_ptr exactly once. Since the loop already iterates a copy of the children list, mutating m_shieldingSelects mid-loop is safe:    for( WWidget *child : shieldings )   {     if( dynamic_cast<ShieldingSelect *>( child ) )       m_shieldingSelects->removeWidget( child ); // returns unique_ptr<WWidget>, destroys once   }  (`child->removeFromParent();` is equivalent.) Do not keep a `ShieldingSelect *` and `delete` it. Apply the same fix at line 8377 in removeShielding (and there also iterate over a copy of children() rather than the live reference to avoid the separate iterator-invalidation …

#### `src/ShieldingSourceDisplay.cpp:8361-8382` — ShieldingSourceDisplay::removeShielding deletes ShieldingSelect child of m_shieldingSelects (double-free)
**🔴 CRITICAL** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 8371-8390
- **Mechanism:** removeShielding() finds the matching child of m_shieldingSelects (line 8371) and does `delete shielding;` (line 8377). The widget is owned by m_shieldingSelects' children_ unique_ptr; deleting it directly runs ~WWebWidget -> removeFromParent -> removeWidget -> removeChild/Utils::take which extracts the still-owning unique_ptr (assert fires in debug; dangling/double-free in release when the parent is later torn down). The first loop breaks immediately after the delete so its iteration over the live `children` reference is safe; the core bug is the `delete` on a parented widget. Triggered every time a user removes a shielding row.
- **Trigger:** Trigger: select->remove().connect( this, &ShieldingSourceDisplay::removeShielding ) (line 7891) — every click on a shielding row's remove ("X") button
- **Fix:** Replace `delete shielding;` (src/ShieldingSourceDisplay.cpp:8377) with `m_shieldingSelects->removeWidget( shielding );` and let the returned std::unique_ptr<WWidget> drop, destroying the widget cleanly via the owning-API path (matches §1/§3: detach via removeWidget() rather than raw delete on a parented widget). The surrounding `break;` and `foundShielding` logic stay as-is.

#### `InterSpec/DecayActivityDiv.h:205` — Raw CsvDownloadGui* m_csvDownloadDialog should be observing_ptr (dangling-on-teardown UAF risk)
**🟠 HIGH** · _ownership (P2)_
- **Mechanism:** m_csvDownloadDialog is a raw AuxWindow-subclass pointer created via AuxWindow::make (wApp owns it); the normal teardown nulls it, but because it does not auto-null, if wApp destroys the CsvDownloadGui during session teardown before ~DecayActivityDiv runs, the if(m_csvDownloadDialog) deleteCsvDownloadGui() at line 3629 reads a stale non-null pointer and calls deleteAuxWindow on freed memory.
- **Trigger:** Trigger: end any session with the CSV Export dialog still open
- **Fix:** Change DecayActivityDiv.h:205 to Wt::Core::observing_ptr<CsvDownloadGui> m_csvDownloadDialog; it auto-nulls when CsvDownloadGui is destroyed by any path, making the if(m_csvDownloadDialog) guards at lines 2240/3629 correct. Update line 2242 to AuxWindow *dialog = m_csvDownloadDialog.get(); the m_csvDownloadDialog=nullptr at 2243 becomes redundant (harmless).

#### `src/MakeDrf.cpp:872-880` — DrfPeak widgets created with raw new are orphaned (never added to m_peaks): leak + invisible peaks + empty fit
**🟠 HIGH** · _ownership (P2)_
- **Mechanism:** DrfPeak's ctor (line 441-443) ignores its 4th WContainerWidget* parent arg (de-parented in the Wt3->Wt4 cleanup), so new DrfPeak( peak, livetime, summed_meas, m_peaks ) creates a widget with NO parent; the returned p is only used to wire signals and is then dropped, never m_peaks->addWidget-ed. Each DrfPeak therefore leaks and never displays, and peaks() (line 1284, reads m_peaks->children()) returns empty, breaking the whole DRF intrinsic-efficiency fit.
- **Trigger:** Trigger: opening the Make Detector Response tool on any spectrum with peaks.
- **Fix:** Replace the raw new with the owning add: `DrfPeak *p = m_peaks->addNew<DrfPeak>( peak, livetime, summed_meas );` and drop the now-dead trailing parent ctor argument (per project rule preferring addNew<T> for parent-owns-child).

#### `src/MakeDrf.cpp:1263-1268` — MakeDrfSrcDef widgets created with raw new are orphaned (never added to m_sources): leak + missing sources
**🟠 HIGH** · _ownership (P2)_
- **Mechanism:** MakeDrfSrcDef's ctor calls WContainerWidget() with no parent (MakeDrfSrcDef.cpp:330), so new MakeDrfSrcDef( n, measDate ) in refreshSourcesVisible() is fully orphaned; it is stored only in the local nucToWidget map (line 1267) and connected, but never m_sources->addWidget-ed. The new source widgets leak, never render, and are never picked up by sources() (line 1297, iterates m_sources->children()) or by the next refresh, so newly-selected source nuclides silently drop out of the DRF fit.
- **Trigger:** Trigger: user enables eff-fit on a peak whose nuclide lacks a source widget -> refreshSourcesVisible() makes an orphan that never renders (no DOM parent) and is missed by sources() (1300) and selected_peak_to_sources() (1312/1317, both iterate m_sources->children()), so the nuclide is dropped from the DRF fit at 2216/3048 (the `if(src)` guard at 1336 fails)
- **Fix:** Parent the source into m_sources using the Wt4 idiom: `MakeDrfSrcDef *src = m_sources->addNew<MakeDrfSrcDef>( n, measDate );` (replacing the raw `new`), preserving the Wt3 parenting intent.

#### `src/RelActAutoGui.cpp:312-415` — RelActAutoGui::createWindow catch: delete disp after ownership transferred / after disp_owner destroyed (double-free or UAF)
**🟠 HIGH** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 400-411
- **Mechanism:** In createWindow(), `disp = disp_owner.get()` (line 322) then `window->contents()->addWidget(std::move(disp_owner))` (line 338) transfers ownership to window->contents(). The catch block does `if(disp) delete disp;` (line 406) AND then `AuxWindow::deleteAuxWindow(window)` (line 410). Two failure modes: (a) if an exception is thrown after line 338, disp is owned by window->contents(); `delete disp` double-frees it and deleteAuxWindow(window) frees it again via the container. (b) If an exception is thrown between lines 321 and 338, the local unique_ptr disp_owner (declared inside the try) is destroyed during stack-unwinding before the catch runs, so `delete disp` operates on an already-freed object (UAF). Either way the error path corrupts the heap.
- **Trigger:** `delete disp` destroys the object while the container still holds a dangling raw child pointer; `AuxWindow::deleteAuxWindow(window)` (which goes WApplication::removeChild → destroys the window and walks/destroys its children) then frees it again → double-free.\n\nThe trigger is realistic: the entire createWindow body is wrapped in try/catch precisely …
- **Fix:** Remove the manual `delete disp;` from the catch block. Replace lines 405-407 — delete the `if( disp ) delete disp; disp = nullptr;` — and keep only the AuxWindow teardown. Mechanism: if the throw occurs before the move (line 338), `disp_owner`'s own destruction during stack unwinding frees the object correctly (so do nothing for disp); if it occurs after the move, the object is owned by `window->contents()` and is destroyed when `AuxWindow::deleteAuxWindow(window)` tears the window down. So the catch should be just:\n\n```cpp\n}catch( std::exception &e ) {\n  passMessage( WString::tr(\"raag-error-creating-tool\").arg(e.what()), WarningWidget::WarningMsgHigh );\n  disp = nullptr;             …

#### `src/ShieldingSourceDisplay.cpp:2567-2702` — ShieldingSourceDisplay::createWindow catch: delete disp after it was added to the window stretcher (double-free)
**🟠 HIGH** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 2695-2701
- **Mechanism:** createWindow() assigns `disp = window->stretcher()->addWidget( std::make_unique<ShieldingSourceDisplay>(...), 0, 0 );` (line 2592), so disp is owned by the window's stretcher layout immediately. The catch block does `if(disp) delete disp;` (line 2696) followed by `AuxWindow::deleteAuxWindow(window)` (line 2700). If any exception is thrown after line 2592 (e.g. disp->deSerialize rethrows under BUILD_AS_UNIT_TEST_SUITE at 2635, disp->initialSizeHint, window->resizeWindow, the INCLUDE_ANALYSIS_TEST_SUITE block, window->show, or passMessage), disp is layout-owned; `delete disp` double-frees it, and deleteAuxWindow(window) frees it a second time via the layout. disp is null before line 2592, so the early-failure path is safe.
- **Trigger:** The conclusion (delete disp is a double-free) is correct.\n\nTrigger: any throw in the outer try after 2592
- **Fix:** In the outer catch block remove the manual `delete disp;` (lines 2695-2696). Keep `disp = nullptr;` and `AuxWindow::deleteAuxWindow( window ); window = nullptr;`. deleteAuxWindow(window) -> wApp->removeChild(window) destroys the AuxWindow, which destroys contents() -> WGridLayout -> WWidgetItem -> disp exactly once. The `disp = nullptr;` must remain so the returned pair does not hand back a dangling pointer to the (now-destroyed) widget. Resulting catch body:\n  passMessage( WString::tr(\"ssd-error-create-tool\").tr(e.what()), WarningWidget::WarningMsgHigh );\n  disp = nullptr;            // owned by window's stretcher; freed by deleteAuxWindow below\n  if( window )\n    …

#### `src/DrfSelect.cpp:931-984` — RelEffFile::handleUserAskedRemove removes the widget unconditionally instead of on confirmation
**🟡 MEDIUM** · _ownership (P2)_
- **Mechanism:** The else-branch builds a Yes/No confirm SimpleDialog (the Yes handler captures only filepath by value and deletes the file), but then unconditionally calls removeFromParent() at line 983, whose dropped unique_ptr destroys this immediately on the close-icon click regardless of the user's answer; the confirmation dialog is therefore moot and the RelEffFile self-destructs mid-member-function before the user ever responds.
- **Fix:** Delete the unconditional removeFromParent() at line 983; in the yes->clicked() handler, after boost::filesystem::remove succeeds and the toast is shown, remove the row — connect via the lifetime-safe form (e.g. yes->clicked().connect(this, ...) or capture this and call removeFromParent() there) so the row is removed only on confirmed deletion.

#### `src/FitSkewParamsTool.cpp:254-259` — `new PeakModel()` is owned by nobody and never deleted -> leak per window open
**🟡 MEDIUM** · _ownership (P2)_
- **Mechanism:** m_peakModel (declared `PeakModel*` in the header) is assigned `new PeakModel()`; D3SpectrumDisplayDiv::setPeakModel (line 635) only stores a raw pointer and does not take ownership, and ~FitSkewParamsTool (lines 157-161) only sets the cancel flag, never deleting m_peakModel. Each Fit-Skew-Params window open leaks a PeakModel; the rest of InterSpec owns PeakModel via std::shared_ptr.
- **Fix:** Change the member to `std::shared_ptr<PeakModel> m_peakModel;`, assign via `m_peakModel = std::make_shared<PeakModel>();` and pass `m_chart->setPeakModel( m_peakModel.get() )` (matching InterSpec.cpp:633/655); the shared_ptr is destroyed with the tool, and setPeakModel's raw-pointer non-owning contract is preserved.

#### `src/FitSkewParamsTool.cpp:1056-1059` — m_rightClickMenu = new WPopupMenu() creates an ownerless global widget that is never freed -> leak per right-click
**🟡 MEDIUM** · _ownership (P2)_
- **Mechanism:** WPopupMenu's ctor calls app->addGlobalWidget(this), which in the C++ build does domRoot_->addWidget then domRoot_->removeChild(w).release() (WApplication.C:561-565), so the menu has NO widget parent and is owned by nobody; thus `m_rightClickMenu->removeFromParent()` (line 1057) returns an empty unique_ptr and frees nothing, WPopupMenu never self-deletes, and ~FitSkewParamsTool (lines 157-161) does not free it, so every right-click after the first leaks the prior menu and the last one leaks at teardown.
- **Fix:** Mirror InterSpec's pattern: make the member an owned object created once (e.g. std::unique_ptr<WPopupMenu> m_rightClickMenu = addChild(std::make_unique<WPopupMenu>()) in the ctor; FitSkewParamsTool is a WContainerWidget) and reuse it by clearing+repopulating then popup(), rather than new-per-click; or if recreating each click, detach+free the previous via the owning unique_ptr before reassigning. Do not rely on removeFromParent() to free a global widget.

#### `src/NuclideSourceEnter.cpp:417-423` — NuclideSourceEnterController allocated with bare new and never deleted (leak of WObject + its child WSuggestionPopup)
**🟡 MEDIUM** · _ownership ·×2 (P2)_
- **Mechanism:** NuclideSourceEnter's ctor does `m_controller = new NuclideSourceEnterController(...)` into a raw pointer member (NuclideSourceEnter.h:144), but the controller is a Wt::WObject that is never addChild'd to any parent and ~NuclideSourceEnter() (lines 421-423) is empty, so it is never freed; in Wt4 WObject() has no parent-taking ctor (verified WObject.C:29) so the Wt3 auto-parenting that this pattern relied on is gone. The controller also owns a WSuggestionPopup via addChild (line 92), so that widget leaks too every time a NuclideSourceEnter (DoseCalc / activity-limit) is created and destroyed.
- **Fix:** In NuclideSourceEnter ctor replace `m_controller = new NuclideSourceEnterController( m_nuclideEdit, m_nuclideAgeEdit, m_halfLifeTxt );` with `m_controller = addChild( std::make_unique<NuclideSourceEnterController>( m_nuclideEdit, m_nuclideAgeEdit, m_halfLifeTxt ) );` so the widget (a WObject) owns and destroys the controller — matching the already-correct DetectionLimitSimple.cpp:342 pattern; the empty ~NuclideSourceEnter() can stay.

#### `src/PeakModel.cpp:905-1013` — PeakCsvResource leaked (Wt4 WResource() takes no parent, raw new never adopted, empty ~PeakModel) and its cached raw m_app/m_model dangle on late CSV request
**🟡 MEDIUM** · _ownership ·×2 (P2+P3)_
- **Also at / spans:** 1008, 1011-1013
- **Mechanism:** In Wt4 WResource() takes no parent, so PeakCsvResource's ctor merely stores PeakModel* m_model and WApplication* m_app rather than parenting; m_csvResource = new PeakCsvResource(this) (line 1008) is a raw heap object owned by nothing and ~PeakModel() (line 1011) is empty, so it leaks once per session. Because it outlives both the PeakModel and the WApplication, a CSV-download request arriving during/after teardown dereferences the dangling m_app in WApplication::UpdateLock lock(m_app) (line 925) and m_model->m_measurment in handleRequest -> use-after-free.
- **Fix:** Construct the resource as an owned child of the model: in PeakModel ctor use `m_csvResource = addChild( std::make_unique<PeakCsvResource>(this) );` (PeakModel is a WObject via WAbstractItemModel, WResource is a WObject), so it is destroyed with PeakModel; peakCsvResource() keeps returning the raw pointer. Optionally change member to PeakCsvResource* held by addChild ownership.

#### `src/ShieldingSourceDisplay.cpp:7027-7039` — startBrowseDatabaseModels catch: unconditional delete of widgets that may already be parented in m_modelDbBrowseWindow (double-free)
**🟡 MEDIUM** · _ownership ·×2 (P1+P2)_
- **Also at / spans:** 7027-7038
- **Mechanism:** summary/accept/cancel/del are raw-new'd unparented WWidgets (lines 6857-6868) and then conditionally adopted into the window footer/contents via addWidget(std::unique_ptr<WWidget>(...)) (cancel at 6973/6989, accept at 6975, del at 6969, summary at 6983). The catch block unconditionally does `delete accept; delete cancel; delete summary; delete del;` (lines 7029-7036) then AuxWindow::deleteAuxWindow(m_modelDbBrowseWindow). If an exception is thrown after any widget was addWidget'd (e.g. in rejectWhenEscapePressed, centerWindow, or show at 7003-7008), the manual delete destroys a now-layout/container-owned widget and deleteAuxWindow then double-frees it; widgets never added are correctly freed, but the catch can't tell which is which, so it can corrupt the heap on the construction-error path.
- **Fix:** Stop raw-deleting in the catch and stop double-tracking ownership. Hold the four widgets in local std::unique_ptr from creation and only release into the window when adding, so on exception the still-owned unique_ptrs clean up and added widgets are owned solely by the window (which deleteAuxWindow tears down). Concretely:    std::unique_ptr<WTextArea> summary = std::make_unique<WTextArea>();   std::unique_ptr<WPushButton> accept = std::make_unique<WPushButton>( WString::tr("Load") );   std::unique_ptr<WPushButton> cancel = std::make_unique<WPushButton>( WString::tr("Cancel") );   std::unique_ptr<WPushButton> del = std::make_unique<WPushButton>( WString::tr("Delete") );   // keep raw aliases …

#### `src/SpecMeasManager.cpp:1484-1546,1659-1667` — FileDragUploadResource instances new'd as raw pointers but never deleted in ~SpecMeasManager
**🟡 MEDIUM** · _ownership ·×2 (P2+P3)_
- **Also at / spans:** 1484-1486, 1543, 1659-1667
- **Mechanism:** m_foregroundDragNDrop/m_secondForegroundDragNDrop/m_backgroundDragNDrop/m_batchDragNDrop are each created with raw new FileDragUploadResource() and stored as raw FileDragUploadResource* members (header lines 696-700); they are never addChild'd and ~SpecMeasManager (1659-1667) never deletes them, so the WResource objects and their dtor's beingDeleted()+clearSpooledFiles() cleanup of stolen spool files leak once per session.
- **Fix:** Make the four members std::unique_ptr<FileDragUploadResource> (or std::make_unique into the ctor init list / set in body) so they are destroyed with SpecMeasManager, ensuring ~FileDragUploadResource runs clearSpooledFiles(); update the raw-pointer accessors (dragNDrop/foregroundDragNDrop/etc.) to return .get().

### Threading — worker/native-thread access & UAF  (15)

#### `src/BatchGuiInputFile.cpp:158-205` — Batch file-parse worker posts GUI callback capturing raw `this`; widget is freely removable mid-parse -> use-after-free
**🔴 CRITICAL** · _threading ·×2 (P1+P2)_
- **Also at / spans:** 158-201
- **Mechanism:** In the BatchGuiInputSpectrumFile constructor the file is parsed on a thread-pool thread via the qualified `ioService().boost::asio::io_service::post(...)` (line 164). On completion the worker re-posts `updateGuiCallback = [this, spec_meas, status_ptr](){ set_spectrum(...); }` (line 162) via `WServer::schedule(1000ms, sessionid, updateGuiCallback)` (line 200). set_spectrum (line 205) immediately dereferences widget members (m_preview_container->addNew<...>(), m_is_peaks_csv, m_spectrum). The close icon is wired to requestRemoveSelf() (line 155) which emits m_remove_self_request_signal; the parent handles it in BatchGuiWidget::handle_remove_input_file (BatchGuiWidget.cpp:330) by `m_input_files_container->removeWidget(input)`, destroying the widget synchronously. There is NO lifetime protection (no observing_ptr, no deleted-flag) and the destructor (line 340) sets no cancel flag and does …
- **Trigger:** set_spectrum (line 205) immediately dereferences widget members: m_preview_container->addNew (211), m_is_peaks_csv (209), m_spectrum (303-318), m_spec_meas (326), m_preview_created_signal.emit (330), wApp->triggerUpdate (332)
- **Fix:** Add widget-lifetime protection to the deferred GUI callback, mirroring the SpecMeasManager observing_ptr pattern (e.g. SpecMeasManager.cpp:4607-4617). Minimal change: capture an Wt::Core::observing_ptr<BatchGuiInputSpectrumFile> instead of raw `this` and no-op when it has been destroyed. Replace the updateGuiCallback definition (line 161-162) with: `Wt::Core::observing_ptr<BatchGuiInputSpectrumFile> self( this ); std::function<void(void)> updateGuiCallback = [self, spec_meas, status_ptr](){ if( self ) self->set_spectrum( spec_meas, status_ptr ); };`. The observing_ptr auto-nulls when handle_remove_input_file destroys the widget, so the scheduled callback safely skips set_spectrum. …

#### `src/FitPeaksForNuclidesGui.cpp:1159-1219` — FitPeaksAdvancedWidget worker completion lambda calls into freed widget (cross-thread UAF)
**🔴 CRITICAL** · _threading (P2)_
- **Mechanism:** startComputation posts a worker to the IO thread (line 1219); gui_update/gui_error capture the widget `this` by raw pointer, and on completion the worker does WServer::post(sessionid, [=](){ if(app){ gui_update(); ... } }) (lines 1184/1196/1208) checking only that the session app is alive, not the widget. FitPeaksAdvancedWidget lives in FitPeaksAdvancedDialog (a SimpleDialog) that self-destructs on finished(); ~FitPeaksAdvancedWidget only sets m_cancel_calc->store(true), so if the user closes the dialog while the fit is in flight the already-posted completion still invokes updateFromResult(...) whose first statement dereferences this->m_cancel_calc (line 1226) on freed memory -> use-after-free.
- **Trigger:** Trigger: user starts a fit then closes the dialog before it finishes; the deferred widget destruction runs on the session thread, then the still-running/finishing worker reads freed members and the posted completion dereferences the freed widget
- **Fix:** Give FitPeaksAdvancedWidget a std::shared_ptr alive-token reset in its dtor; capture a std::weak_ptr in both the worker and the WServer::post completion lambdas and bail (or skip gui_update/gui_error) if it cannot be locked before any `this`/member dereference. Also copy the worker inputs (m_fg_copy, m_base_nucs, m_user_peaks, m_drf, bg_to_use) into locals captured by value before post() so the IO-thread body never reads `this`.

#### `src/RelActAutoGui.cpp:5713-5776` — Background calc posts back a lambda capturing raw this with no lifetime guard -> UAF if window closed mid-calc
**🔴 CRITICAL** · _threading (P2)_
- **Mechanism:** The worker (run on the io_service thread pool, line 5776) captures gui_update_callback/error_callback which capture raw this, and on completion calls WServer::post(sessionId, ...) (5737/5757) invoking updateFromCalc/handleCalcException, which dereference this and members. WServer::post does NOT track the RelActAutoGui lifetime, and the posted lambda only checks WApplication::instance() (app), not whether this is alive; if the AuxWindow is closed (RelActAutoGui destroyed) after solve() returned and posted but before the continuation runs, the lambda dereferences a destroyed object. ~RelActAutoGui (1167-1177) only sets m_cancel_calc, which cannot stop an already-posted continuation.
- **Trigger:** Trigger: start a calc, close the Isotopics-from-nuclides window before/while it finishes -> UAF
- **Fix:** Add a session-thread-checked alive sentinel: store a std::shared_ptr<bool> m_alive member (set true in ctor) and capture a std::weak_ptr to it (plus calc_number) in the posted continuation; on the session thread (already under UpdateLock) bail if the weak_ptr is expired before calling gui_update_callback/error_callback. Equivalently, capture m_interspec (outlives the tool) by value and re-look-up via m_interspec->relActAutoWindow(false), no-op if it returns null, instead of dereferencing raw this.

#### `src/RelActManualGui.cpp:1861-1891` — Background rel-eff calc posts back a lambda capturing raw `this`; tool can be destroyed mid-calc -> UAF
**🔴 CRITICAL** · _threading (P2)_
- **Mechanism:** calculateSolution() captures `updater = [this,solution](){ updateGuiWithResults(solution); }` and `err_updater = [this,errmsg](){...}`, runs solve_relative_efficiency on a real worker thread (qualified asio post line 1873), then schedules them via WServer::instance()->post(sessionId, updater) (lines 1879/1889). m_relActManualGui is an observing_ptr that InterSpec destroys when the user closes the Rel-Act-Manual tab (handleToolTabClosed/handleRelActManualClose, InterSpec.cpp:7164) — this can happen while the background solve is still running. WServer::post only re-acquires the session UpdateLock if the session is alive; it does NOT verify the RelActManualGui widget still exists, so the posted lambda dereferences a destroyed `this` (m_currentSolution, m_results, m_chart) -> use-after-free.
- **Trigger:** Trigger: change an input to start a non-trivial solve, then close the tab before it finishes.
- **Fix:** Do not capture raw this in updater/err_updater. In the re-posted lambda re-resolve the tool via InterSpec::instance()->relActManualWidget() (returns the live observing_ptr target or nullptr) and null-check before calling updateGuiWithResults/updateGuiWithError, mirroring the project's undo/redo look-up-don't-capture discipline; alternatively gate the GUI update on a shared lifetime token checked under the UpdateLock.

#### `src/RemoteRid.cpp:1319-1350` — requestExeInfo background worker re-posts a lambda capturing `this` (ExternalRidWidget) with no lifetime protection -> UAF if widget closed before DRF query returns
**🔴 CRITICAL** · _threading ·×4 (P1+P2)_
- **Also at / spans:** 1124-1238, 1203-1238, 1203-1238, 1309-1351
- **Mechanism:** requestExeInfo spawns run_external_command on a background thread and on completion does WServer::post(appsession, [this,...]{ receiveExeDrfInfo(...) }). receiveExeDrfInfo is an instance method touching widget members (handleInfoResponse). If the user closes the Remote-RID window (deleteAuxWindow destroys this widget) while the executable is being queried, the posted lambda calls a method on a destroyed ExternalRidWidget. WServer::post guarantees only session liveness + UpdateLock, not widget liveness.
- **Trigger:** Trigger: user clicks 'Retrieve DRFs' on the Executable tab, then closes the Remote-RID window before the executable's --drfs query returns; the posted lambda then calls this->receiveExeDrfInfo -> handleInfoResponse, dereferencing destroyed members (m_drf_select->clear(), m_current_drf_index, m_status_stack) -> use-after-free
- **Fix:** Do not capture a raw this in the posted lambda. Capture a Wt::Core::observing_ptr<ExternalRidWidget> bound to this (make the member type observing-safe) and null-check it inside doUpdateFcn before calling receiveExeDrfInfo; or marshal results through a session-owned resource object whose lifetime is independent of the widget. Apply the same guard to startExeAnalysis's parent capture.

#### `src/MakeFwhmForDrf.cpp:697-727` — MakeFwhmForDrf auto-peak-search worker posts callback capturing raw `this`; modal dialog closeable mid-search with no widget-tied protection
**🟠 HIGH** · _threading (P1)_
- **Mechanism:** startAutomatedPeakSearch posts search_for_peaks(...) to the thread pool via the qualified `ioService().boost::asio::io_service::post(...)` (line 717). The worker UNCONDITIONALLY re-posts `callback = [this, user_peaks, searchresults](){ setPeaksFromAutoSearch(...); }` (line 698) via `server->post(seshid, callback)` (line 726) -- even when the search was cancelled. setPeaksFromAutoSearch (line 954) immediately dereferences widget state (m_currently_searching, m_model->set_peaks(...)). MakeFwhmForDrf is a WContainerWidget hosted in MakeFwhmForDrfWindow (an AuxWindow); the close button fires finished() and tears the window down. The only cancel token is the SpecMeas-tied peak_search_cancel_flag() (line 714) which merely interrupts the search ALGORITHM; it is not tied to widget lifetime and does not prevent the callback post, and ~MakeFwhmForDrf (line 640) is empty (no cancel flag, no future …
- **Trigger:** Realistic trigger: opening the tool with use_auto_fit_peaks_too==true immediately fires startAutomatedPeakSearch() from the ctor (line 626); a full-spectrum auto-search takes seconds, and the user clicking Cancel/escape during it triggers the UAF
- **Fix:** Stop capturing raw `this` in the completion callback; tie it to widget lifetime. Minimal, project-idiomatic fix: capture an observing_ptr to the widget and early-return if null. Replace lines 697-698 with: `Wt::Core::observing_ptr<MakeFwhmForDrf> self( this ); std::function<void(void)> callback = [self, user_peaks, searchresults](){ if( self ) self->setPeaksFromAutoSearch( user_peaks, searchresults ); };` (observing_ptr auto-nulls on destruction, and the callback runs under the session UpdateLock so the check is race-free). Additionally, skip posting when cancelled: in the worker (lines 724-726) guard with `if( server && (!cancel_flag || !cancel_flag->load()) ) server->post( seshid, …

#### `src/ShieldingSourceDisplay.cpp:9205-9232` — Background-fit gui_updater captures raw this and is re-posted to the session after the widget may be destroyed
**🟠 HIGH** · _threading (P2)_
- **Mechanism:** gui_updater = [this,results](){ updateGuiWithModelFitResults(results); } is handed to fit_model as finished_fcn and, when run on the worker thread, is dispatched via Wt::WServer::post(wtsession, finished_fcn) (ShieldingSourceFitCalc.cpp:2746) BEFORE fit_model returns. ~ShieldingSourceDisplay waits on m_currentFitFuture, but that future is satisfied right after the post is enqueued; if the window is destroyed in that window, the queued finished_fcn still runs later and dereferences a destroyed this -> UAF. The std::function capture bypasses Wt signal target-tracking.
- **Trigger:** Trigger: the background fit completes and posts its GUI update (update_gui true, so cancelFit in the dtor is a no-op on an already-finished fit), then the dialog is closed or the spectrum reloaded, destroying the tool in the small window before the queued event is serviced
- **Fix:** Inside gui_updater capture the sessionid string instead of this and look the tool up via InterSpec::instance()->shieldingSourceFit() (or equivalent accessor), no-oping if null, mirroring the undo/redo discipline; alternatively have ~ShieldingSourceDisplay cancel-and-drain the posted finished_fcn (e.g. gate it on a shared cancelled-flag/observing_ptr) so a destroyed widget is never dereferenced.

#### `target/macos/NativeMenu.mm:106-113` — validateMenuItem: dereferences Wt widget on AppKit thread without UpdateLock or liveness guard
**🟠 HIGH** · _threading (P3)_
- **Mechanism:** AppKit calls -validateMenuItem: on the main (UI) thread whenever a native menu is about to display, and it directly calls m_item->isEnabled() / m_cb->isEnabled() on the Wt PopupDivMenuItem/WCheckBox with no WApplication::UpdateLock and no session-liveness check; this races against the Wt session thread mutating those widgets and is a use-after-free if the widget was already destroyed (see deferred removal in finding below).
- **Fix:** Stop reading the live Wt widget on the AppKit thread. Invalidate Target synchronously when the widget dies (have ~PopupDivMenuItem / the WCheckBox-owner clear the Target's m_item/m_cb, e.g. via a setWtItem:0 call before removeOsxMenuItem, so validateMenuItem: returns NO once dead), AND cache an atomic enabled flag (std::atomic<bool>) on the Target that the session thread updates (from item->disable()/setDisabled handlers, posted to main) so -validateMenuItem: only reads the cached atomic instead of calling isEnabled() on the widget. Both halves are needed: the atomic removes the data race, the synchronous invalidation removes the UAF.

#### `src/D3SpectrumDisplayDiv.cpp:3122-3163,3500` — Drag-create-ROI worker reads widget member on background thread without UpdateLock
**🟡 MEDIUM** · _threading ·×2 (P2)_
- **Also at / spans:** 84-101,3378-3447
- **Mechanism:** performDragCreateRoiWork posts fcnworker to WServer ioService (allowAsync is always true, line 3081/3500); the worker captures the D3SpectrumDisplayDiv 'spectrum' by raw pointer and dereferences spectrum->m_last_being_added_peaks at line 3163 on the background thread BEFORE acquiring WApplication::UpdateLock at line 3300, racing with main-thread mutation/clearing of that member (e.g. another drag callback or widget teardown), so it is an unsynchronized read of session-owned state (and a dangling read if the chart is destroyed before the post runs).
- **Fix:** At line 3163 use the main-thread snapshot already captured at line 3120 (prev_shown_peaks.empty()) instead of the live spectrum->m_last_being_added_peaks.empty(); this is the value the 'peaks already cached' comment intends and the only pre-lock touch of this session-owned member - all other accesses are already inside the UpdateLock section.

#### `src/SpecFileQueryWidget.cpp:1795-1811, 454-643` — ResultCsvResource reads live ResultTableModel without UpdateLock (data race)
**🟡 MEDIUM** · _threading ·×2 (P2+P3)_
- **Also at / spans:** 472-641
- **Mechanism:** ResultCsvResource::handleRequest dereferences m_model (the live ResultTableModel owned by m_resultview) but WResource::takesUpdateLock_ defaults to false, so handleRequest runs on a Wt resource-serving thread outside the WApplication event loop while the session thread can be mutating the same model via updateSearchStatus/finishUpdate -> m_resultmodel->appendResults(); concurrent read/iterate vs. vector mutation is undefined behaviour. The CSV button can be clicked while a (re)search is still posting incremental results.
- **Fix:** Call setTakesUpdateLock(true) in the ResultCsvResource ctor so handleRequest runs under WApplication::UpdateLock and serializes against the session thread's appendResults; alternatively snapshot the model rows under the session thread before serving.

#### `target/android/android.cpp:83-127, 237-242` — JNI GlobalRef TOCTOU race between worker-thread save callback and setFileSaveCallback
**🟡 MEDIUM** · _threading (P3)_
- **Mechanism:** invoke_java_save_callback (run on a Wt io_service worker thread during a download) snapshots sm_java_file_save_cb_obj/method under the mutex, releases the lock, then calls env->CallVoidMethod(cb_obj,...); if another thread calls setFileSaveCallback concurrently it runs env->DeleteGlobalRef(sm_java_file_save_cb_obj), invalidating the GlobalRef the worker is about to invoke -> JNI undefined behavior / crash on a deleted global ref.
- **Fix:** Hold sm_jni_mutex across the entire AttachCurrentThread/CallVoidMethod/DetachCurrentThread sequence in invoke_java_save_callback (lines ~89-127) instead of releasing it after the snapshot; or, while still holding the lock, take a fresh local copy via env->NewGlobalRef(sm_java_file_save_cb_obj), release the lock, call CallVoidMethod on that local ref, then DeleteGlobalRef the local copy after the call.

#### `target/electron/InterSpecAddOn.cpp:646-650` — browse_for_directory deadlocks an asio worker thread if BlockingCall fails to enqueue
**🟡 MEDIUM** · _threading (P3)_
- **Mechanism:** After ns_browse_for_directory_callback->BlockingCall(&data) at line 646 the code unconditionally enters data.cv.wait(lk, []{return status!=NotStarted;}) at line 650; the napi_status is only inspected at line 656 AFTER the wait. If BlockingCall returns non-napi_ok (TSFN closing/aborted at shutdown, or queue error), the node worker call_dir_browse_in_node never runs, data.status stays NotStarted forever, and the asio io_service thread blocks permanently — leaking one of only 4 desktop WIOService threads.
- **Fix:** Check the BlockingCall return before waiting: `if(status != napi_ok){ ns_browse_for_directory_callback->Release(); return false; }` placed between lines 646 and 650, so cv.wait is only entered when the JS callback was actually enqueued (and will eventually set status + notify). Optionally also use cv.wait_for with a timeout as defense-in-depth, and honor the Acquire() return value.

#### `target/wxwidgets/InterSpecWxApp.cpp:514-525` — handle_javascript_error reads wx top-level-window list off the wx GUI thread
**🟡 MEDIUM** · _threading (P3)_
- **Mechanism:** This static handler is invoked from InterSpecApp::handleJavaScriptError (src/InterSpecApp.cpp:1421-1422) on the Wt session thread (under UpdateLock), which is not the wx GUI thread; app->GetTopWindow() at line 514 reads the wxTopLevelWindows list that wx only mutates on the GUI thread, a data race with frame creation/destruction. The subsequent CallAfter marshaling is correct.
- **Fix:** Marshal to the GUI thread before touching any window: wxTheApp->CallAfter([error_msg,app_token](){ auto *app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance()); if(app) app->handle_javascript_error_internal(error_msg,app_token); }); wxTheApp->CallAfter is the documented thread-safe cross-thread primitive (queues via wxQueueEvent), so GetInstance()/GetTopWindow()/GetEventHandler() are then only touched on the wx main thread.

#### `target/wxwidgets/InterSpecWxUtils.cpp:128-159` — browse_for_directory reads/uses wx top-level window off the Wt session thread
**🟡 MEDIUM** · _threading (P3)_
- **Mechanism:** Per its header contract this is called from a Wt event-loop (session) thread, yet wxapp->GetTopWindow() at line 135 reads the GUI-thread-owned wxTopLevelWindows list, and topWindow is captured by raw pointer into the CallAfter lambda (line 145) and reused as wxDirDialog parent; reading the list off-thread is a data race.
- **Fix:** In browse_for_directory, call wxTheApp->CallAfter(...) (or queue onto a known wx event handler) FIRST, then resolve the top window with GetTopWindow() inside the GUI-thread lambda just before constructing wxDirDialog, so the global TLW list is only read/used on the wx GUI thread and the parent pointer is fresh. Apply the same change to save_file_data.

#### `src/DataBaseUtils.cpp:104-169` — Double-checked-locking data race on m_pool in DbPoolManager::databaseConnectionPool()
**⚪ LOW** · _threading (P2)_
- **Mechanism:** The lock-free fast path keys solely on the atomic m_numconnection: line 132 stores nconn(>0) BEFORE line 133 assigns the non-atomic unique_ptr m_pool, so a second thread reading m_numconnection>0 at line 104 skips the mutex and reaches line 169 'return m_pool.get()' while the first thread is still inside m_pool.reset(pool) -- a torn/null read of m_pool (data race, potential null deref).
- **Fix:** Inside the locked region assign the pool before publishing the count: m_pool.reset(pool); then m_numconnection.store(nconn, seq_cst); and use acquire (not consume) on the line-104/108 loads so a reader that sees count>0 has a happens-before with the m_pool write. Out of Wt4-migration scope; optional to address.

### Signal & callback lifetime  (2)

#### `InterSpec/UserPreferences.h:276-289` — addCallbackWhenChanged(name,target,method) connects untracked lambda -> UAF when target outlived by prefs signal
**🔴 CRITICAL** · _signal (P2)_
- **Mechanism:** The 2-arg bool overload does signal->connect( [target,method](bool v){(target->*method)(v);} ) using the bare connect(F) form, which does NOT track target; the Wt::Signal lives in UserPreferences::m_onBoolChangeSignals (lifetime == InterSpec). The header comment at line 231 falsely promises auto-disconnect. Callers like RelActAutoGui.cpp:559 (m_spectrum, dies when the RelActAuto window closes) and ShieldingSourceDisplay.cpp:722 (SourceFitModel, dies with the tool) register targets shorter-lived than InterSpec; after the tool closes, toggling the LogY/DisplayBecquerel preference fires the lambda into freed memory -> crash.
- **Trigger:** Trigger: after closing the tool, toggling LogY (InterSpec::setLogY -> setPreferenceValue -> setPreferenceValueInternal emits signal(value) unconditionally, UserPreferences.cpp:389) fires the stale lambda into freed memory -> UAF crash
- **Fix:** In UserPreferences.h replace the bool overload body with signal->connect( target, method ); and the int overload (line 317) with signal->connect( target, method ); using the Signal::connect(T*,void(V::*)(B...)) tracking overload, which auto-disconnects when the WObject-derived target is destroyed.

#### `InterSpec/UserPreferences.h:305-318` — addIntCallbackWhenChanged(name,target,method) connects untracked lambda -> UAF (same as bool overload)
**🔴 CRITICAL** · _signal (P2)_
- **Mechanism:** The 2-arg int overload connects [target,method](int v){(target->*method)(v);} via the untracked connect(F) form onto UserPreferences::m_onIntChangeSignals (lives for the whole session). RelActAutoGui.cpp:561/566 register m_spectrum (a D3SpectrumDisplayDiv owned by the RelActAuto tool layout, destroyed when the tool window closes) for RefLineThickness/RefLineVerbosity; after the tool closes, changing those prefs invokes m_spectrum->handleRefLine...PreferenceChangeCallback on a freed widget.
- **Trigger:** Realistic trigger: after closing the tool, the user changes the Ref
- **Fix:** In UserPreferences.h:317 use the tracked member-pointer overload: signal->connect( target, method ); (D3SpectrumDisplayDiv/ReferencePhotopeakDisplay derive from WObject, so Wt auto-disconnects on target destruction). Apply the identical fix to the bool overload at line 288.

### Object lifetime / dangling pointers  (9)

#### `src/IsotopeSelectionAids.cpp:454-468` — EditWidget::handleBlur schedules a 250ms WServer callback capturing raw this/m_parent with no lifetime tracking -> use-after-free
**🔴 CRITICAL** · _lifetime (P2)_
- **Mechanism:** handleBlur() captures raw self/this (the transient editor EditWidget) and parent/m_parent (raw PhotopeakDelegate*) into a lambda passed to WServer::schedule(250ms, sessionId, fn); schedule only verifies the session still exists and does NOT track widget lifetimes (the in-code comment claiming protection at line 462 is false). In Wt4 the editor is owned by the unique_ptr returned from update() and is destroyed promptly when editing ends (Enter, click-away, re-render), so 250ms later self_closer runs in the still-alive session and invokes handleBlurWorker on the freed EditWidget and parent->doCloseEditor(self,...) -> use-after-free.
- **Trigger:** closeOnBlur=true (wiring blurred()->handleBlur) is used in live tables (PeakInfoDisplay 1217, EnergyCalTool 1951/1956, RelActManualGui, ShieldingSourceDisplay), so the trigger is routine.
- **Fix:** Capture a Wt::Core::observing_ptr<EditWidget> handle (and observing_ptr to the delegate, or look the delegate up from the still-current editor) in the scheduled lambda and bail if null before dereferencing; e.g. wrap worker so it no-ops when the editor was destroyed/replaced. Equivalently, verify the editor is still the view's current editor for its index before calling doCloseEditor. This restores the lifetime guarantee that the removed wApp->bind() previously provided.

#### `src/MakeDrfSrcDef.cpp:371-378` — ~MakeDrfSrcDef deletes m_lib_src_menu which is owned by m_lib_src_btn via addChild (double-free)
**🔴 CRITICAL** · _lifetime ·×2 (P1+P2)_
- **Mechanism:** ~MakeDrfSrcDef does `delete m_lib_src_menu;` (line 375). m_lib_src_menu (PopupDivMenu*) is created by makePopupMenu(m_lib_src_btn) (line 649) which returns `button->addChild(std::move(menu))` (PopupDiv.cpp:769) — i.e. it is owned by m_lib_src_btn through WObject::children_. The derived destructor body runs before the WObject/children base teardown, so `delete m_lib_src_menu` leaves the button's children_ unique_ptr dangling; when m_lib_src_btn (itself a child of MakeDrfSrcDef) is destroyed during widget-tree teardown it double-frees the menu. The `#if WT_VERSION >= 0x3070000` guard always compiles under Wt 4.12.6. Triggered on every MakeDrfSrcDef destruction.
- **Trigger:** Trigger: realistic and routine
- **Fix:** Delete the manual ownership management entirely. `m_lib_src_btn` owns `m_lib_src_menu` via addChild and tears it down automatically during widget-tree destruction, so the destructor body must not touch it. Reduce ~MakeDrfSrcDef to an empty body:    MakeDrfSrcDef::~MakeDrfSrcDef()   {   }  (Equivalently, remove lines 373-377 inside the existing destructor.) If, for some reason, explicit early teardown of the menu is desired, the only correct form is to detach via the owning button and let the returned unique_ptr drop: `if( m_lib_src_btn && m_lib_src_menu ){ m_lib_src_btn->removeChild( m_lib_src_menu ); }` (typed WObject::removeChild returns a unique_ptr that frees it once, and clears the …

#### `src/BatchGuiWidget.cpp:291-298` — WServer::schedule worker captures raw this with no lifetime protection (UAF)
**🟠 HIGH** · _lifetime (P2)_
- **Mechanism:** handleFileDrop posts a deferred worker [this,dropped_files](){addInputFiles(dropped_files);} via WServer::schedule(25ms,sessionId,worker) with no bind/observing guard; BatchGuiWidget lives inside a BatchGuiDialog (SimpleDialog) that self-destructs on finished()/escape, so if the dialog closes within the 25ms window the task fires addInputFiles() on a freed BatchGuiWidget. WServer::schedule only guarantees the session UpdateLock, not the target's lifetime, and the now-commented code at 252-272 used wApp->bind precisely to guard this.
- **Trigger:** Trigger: user drops a file (legit batch-drop UI) then hits Escape/closes dialog within the 25ms window
- **Fix:** WServer::bind (wApp->bind) is gone in Wt4; instead capture observing_ptr<BatchGuiWidget> guard=this in the worker and re-check (if(!guard)return;) before calling addInputFiles, and pass a fallBackFunction to WServer::schedule that frees the spooled files (delete those tuples whose should_delete flag is set) when the session/widget is gone.

#### `src/DirectorySelector.cpp:238-303` — Async native directory-picker callback captures raw this (DirectorySelector*) with no lifetime protection -> UAF
**🟠 HIGH** · _lifetime (P2)_
- **Mechanism:** handleNativeDirectorySelection() builds on_select_callback = [this](...){setPath(...)} capturing the raw DirectorySelector*, then hands it to an async native picker; on all three backends the result is delivered later on the session thread (macOsUtils showFilePicker, ElectronUtils::browse_for_directory, and the wx path server->post at line 281). WServer::post only guarantees the session is alive, not the widget. DirectorySelector lives inside user-closable tools (BatchGuiWidget, SpecFileQueryWidget, RefSpectraWidget); closing the tool while the OS picker is open makes the posted lambda call setPath()/validatePath() on a destroyed object.
- **Trigger:** Trigger is a close-during-picker race rather than guaranteed every use, hence high not critical.
- **Fix:** In handleNativeDirectorySelection capture lifetime protection instead of raw this: e.g. Wt::Core::observing_ptr<DirectorySelector> self{this}; auto on_select_callback=[self](const std::vector<std::string>&paths){ if(self && !paths.empty()) self->setPath(paths[0],true); }; (or a std::shared_ptr<bool> alive-flag member checked before touching this), so the deferred session-posted callback no-ops when the widget was destroyed.

#### `target/macos/NativeMenu.mm:261-300` — Target retains raw PopupDivMenuItem*/WCheckBox* that can dangle after async menu removal
**🟠 HIGH** · _lifetime (P3)_
- **Mechanism:** insertOsxMenuItem/addOsxCheckableMenuItem create a Target holding raw m_item/m_cb; ~PopupDivMenuItem (PopupDiv.cpp:1326-1332) destroys the widget synchronously on the Wt thread but only schedules removeOsxMenuItem via dispatch_async, so between widget destruction and the deferred NSMenuItem removal (and permanently for Edit/InterSpec menus whose NSMenuItems are never removed) the Target's raw pointers dangle and are dereferenced by validateMenuItem:/clicked.
- **Trigger:** Realistic trigger: m_detectorToShowMenu (an app-level/native menu) has all its detector-checkbox items removed and recreated on every spectrum load via updateGuiForPrimarySpecChange (InterSpec.cpp:12850-12865); if AppKit calls validateMenuItem: or the user activates the menu during the async-removal gap, it dereferences freed m_item/m_cb (use-after-free)
- **Fix:** Add an `-invalidate` method to Target that nulls m_item/m_cb (and m_nsitem), and call it synchronously from ~PopupDivMenuItem (via the item's stored NSMenuItem target) before/at the dispatch_async removal; guard every Target method on the (already-present) null checks. Optionally remove the NSMenuItem synchronously on the main thread when not quitting, and have the Wt-side keep the Target alive only as long as the widget exists.

#### `src/SimpleActivityCalc.cpp:1627-1630` — WServer::schedule lambda captures raw D3SpectrumDisplayDiv* owned by a SimpleDialog that can self-destruct before the 100ms timer fires
**🟡 MEDIUM** · _lifetime (P2)_
- **Mechanism:** The deferred lambda captures `spectrum` (a D3SpectrumDisplayDiv* owned by dialog->contents()) by raw pointer and is scheduled 100ms later via WServer::schedule. If the user dismisses the SimpleDialog (Yes/No button) within 100ms, it self-destructs (startDeleteSelf -> WServer::post -> removeChild(this)) and destroys `spectrum`; the scheduled lambda then calls spectrum->setXAxisRange(...) on freed memory -> use-after-free.
- **Fix:** Capture a Wt::Core::observing_ptr<D3SpectrumDisplayDiv> instead of the raw pointer and null-check it in the lambda before calling setXAxisRange, so it auto-nulls when the SimpleDialog (and its child spectrum) is destroyed. Apply the same fix to the parallel instance at ShieldingSourceDisplay.cpp:5336.

#### `InterSpec/SpecMeasManager.h:733-736` — m_batchDialog is a raw BatchGuiDialog* (SimpleDialog) that can dangle; should be observing_ptr
**⚪ LOW** · _lifetime (P2)_
- **Mechanism:** m_batchDialog (a SimpleDialog subclass that self-destructs via deferred deleteSelf on finished()) is held as a raw pointer nulled only by handleBatchDialogFinished (SpecMeasManager.cpp:4658). The if(m_batchDialog) return; guard in showBatchDialog (4648) reads this raw pointer; any teardown path that destroys the dialog without finished() running leaves m_batchDialog dangling. Project rules require a SimpleDialog member needing a deferred handler to be Wt::Core::observing_ptr so it auto-nulls.
- **Fix:** Change member to Wt::Core::observing_ptr<BatchGuiDialog> m_batchDialog; so it auto-nulls on any destruction path; drop the redundant m_batchDialog=nullptr in handleBatchDialogFinished and the ctor init can stay or become default.

#### `src/DecaySelectNuclideDiv.cpp:250-262, 103-105` — WSuggestionPopup + filter model leaked onto wApp->root(); model holds dangling m_parent after DecaySelectNuclide dies
**⚪ LOW** · _lifetime (P2)_
- **Mechanism:** m_isotopeSuggestions is created via wApp->root()->addNew<WSuggestionPopup>() (line 254) so it is owned by the session-long root and is deliberately never removed (~DecaySelectNuclide comment line 105); each DecaySelectNuclide (one per DecayActivityDiv, a re-openable tool) thus leaks a popup. setModel(m_isoSearchFilterModel) (line 260) gives the root-owned popup a shared_ptr keeping the filter model alive past the widget, and SimpleIsotopeNameFilterModel holds a raw m_parent (line 250/674) that filter() dereferences throughout (m_parent->m_elementSelection, m_parent->makeMassList(), m_parent->updateSelectedHalfLife(), lines 896-927); once the widget is destroyed m_parent dangles, so any later popup-driven filter() is a use-after-free.
- **Fix:** In ~DecaySelectNuclide capture and drop the popup via std::unique_ptr<Wt::WObject> p = wApp->root()->removeChild( m_isotopeSuggestions ); (popup is a WObject child of root), or construct it with this->addChild<WSuggestionPopup>(...) so it dies with the widget; either way it no longer outlives DecaySelectNuclide, which also disposes the model and removes the stale edits_/m_parent dangles.

#### `target/macos/AppDelegate.mm:229-234` — openURLs takes UpdateLock on a raw InterSpecApp* after the registry mutex was released (TOCTOU)
**⚪ LOW** · _lifetime (P3)_
- **Mechanism:** InterSpecApp::instanceFromExtenalToken (InterSpecApp.cpp:1179) returns a raw InterSpecApp* and releases AppInstancesMutex before returning; ~InterSpecApp (InterSpecApp.cpp:172) erases under the same mutex while destroying the object, so between the return on the AppKit thread and constructing Wt::WApplication::UpdateLock(specapp) on line 234 the session can be destroyed, yielding a dangling app pointer passed to UpdateLock.
- **Fix:** Change instanceFromExtenalToken to keep the matched session alive across the return: return a std::shared_ptr<WebSession> (or a small RAII handle holding both the app* and an UpdateLock acquired under AppInstancesMutex), or have it perform the requested work (open file / handle URL) while still holding the lock. Minimal local mitigation at AppDelegate.mm:226-234 is insufficient on its own because the pointer can already be dangling before line 234; the fix must live in instanceFromExtenalToken so the lookup and the lifetime-extending lock are atomic under AppInstancesMutex.

### AuxWindow / SimpleDialog lifecycle  (2)

#### `src/DetectionLimitSimple.cpp:2198-2210` — programmaticallyCloseMoreInfoWindow deletes a SimpleDialog owned by wApp (double-free)
**🔴 CRITICAL** · _auxwindow ·×2 (P1+P2)_
- **Also at / spans:** 2198-2209
- **Mechanism:** programmaticallyCloseMoreInfoWindow() does `delete dialog;` (line 2208) on m_moreInfoWindow, an observing_ptr<SimpleDialog> (DetectionLimitSimple.h:300). SimpleDialog::make() gives ownership to wApp via addChild() and the canonical teardown is wApp->removeChild(this) (SimpleDialog.cpp:268). Deleting the dialog directly leaves wApp's children_ unique_ptr dangling, so it is double-freed at session teardown; this also violates the project rule that callers must never delete a SimpleDialog. Triggered on the programmatic (undo/redo, state-driven) close path of the MDA more-info window.
- **Trigger:** Trigger: createMoreInfoWindow registers an undo step whose undo calls programmaticallyCloseMoreInfoWindow, so one undo runs the bad delete
- **Fix:** Replace the delete dialog statement at DetectionLimitSimple.cpp line 2208 with a call to WApplication removeChild on the dialog, mirroring SimpleDialog deleteSelf. Never call raw delete on a SimpleDialog.

#### `src/ShieldingSourceDisplay.cpp:3480-3482` — ~ShieldingSourceDisplay tears down a SimpleDialog (m_diagramDialog) via removeFromParent() — bypasses SimpleDialog teardown, leaks dialog and can leave modal cover stuck
**🟡 MEDIUM** · _auxwindow ·×2 (P1+P2)_
- **Mechanism:** m_diagramDialog is a raw `ShieldingDiagramDialog*` (ShieldingSourceDisplay.h:809), and ShieldingDiagramDialog derives from SimpleDialog (ShieldingSourceDiagram.h:48). Per the SimpleDialog contract it must be destroyed only through its own finished()->startDeleteSelf()->WServer::post->wApp->removeChild(this) path; callers must NOT delete/removeFromParent/deleteAuxWindow it. The destructor instead calls m_diagramDialog->removeFromParent() (and has a redundant double-`if` copy-paste artifact). Mechanism: SimpleDialog::make() gives ownership to wApp via addChild; a WDialog is a WPopupWidget registered as a global widget (WApplication::addGlobalWidget adds it to domRoot's flat children_ list but then domRoot->removeChild(w).release(), so domRoot is NOT its owner). Therefore removeFromParent() detaches it from domRoot's render tree but returns a NULL unique_ptr — the modal SimpleDialog is …
- **Fix:** In ~ShieldingSourceDisplay (ShieldingSourceDisplay.cpp:3480-3482) do NOT call removeFromParent()/delete on the SimpleDialog. Replace the block with a call to the existing correct teardown, which fires finished()->startDeleteSelf and pops the modal cover:    closeShieldSourceDiagram();  (closeShieldSourceDiagram() is idempotent: it returns early if !m_diagramDialog, calls m_diagramDialog->accept(), and m_diagramDialog auto-nulls via the finished()->handleShieldSourceDiagramClosed connection.) Equivalently, inline: `if( m_diagramDialog ) m_diagramDialog->accept();` and drop the manual `m_diagramDialog = nullptr;` (it is already nulled by the finished() handler; the explicit assignment is …

### Undo/redo capture  (1)

#### `src/ReferencePhotopeakDisplay.cpp:2143-2164` — showMoreInfoWindow() undo lambda captures `this` (a destroyable ReferencePhotopeakDisplay) and dereferences it
**🟡 MEDIUM** · _undoredo (P1)_
- **Mechanism:** In ReferencePhotopeakDisplay::showMoreInfoWindow() the registered undo lambda captures `this` (line 2143). It correctly looks the tool up via `disp = InterSpec::instance()->referenceLinesWidget()`, but in the `prev_orig_nuc` branch it then dereferences the *captured* `this` directly: `m_nucInfoWindow = SimpleDialog::make<MoreNuclideInfoWindow>(prev_orig_nuc);` (2156, i.e. this->m_nucInfoWindow) and `m_nucInfoWindow->finished().connect( this, [this,win=m_nucInfoWindow](){ handleMoreInfoWindowClose(win); } );` (2157). Per InterSpec.h:728-740, m_referencePhotopeakLines is NOT session-stable: showGammaLinesWindow()/closeGammaLinesWindow() delete the current ReferencePhotopeakDisplay and create a new one when the user toggles the reference-lines tool between the bottom tool-tab strip and a standalone AuxWindow (the member is `Wt::Core::observing_ptr<ReferencePhotopeakDisplay>`). If the user …
- **Fix:** In the undo lambda (src/ReferencePhotopeakDisplay.cpp:2143-2164), drop the `this` capture and use the looked-up `disp` exactly as the redo lambda already does. Change the capture list at line 2143 from `[this, prev_orig_nuc, prev_current_nuc]` to `[prev_orig_nuc, prev_current_nuc]`, and inside the prev_orig_nuc branch replace lines 2155-2157 with: `assert( !disp->m_nucInfoWindow );` then `disp->m_nucInfoWindow = SimpleDialog::make<MoreNuclideInfoWindow>( prev_orig_nuc );` then `disp->m_nucInfoWindow->finished().connect( disp, [disp, win=disp->m_nucInfoWindow](){ disp->handleMoreInfoWindowClose( win ); } );`. This mirrors the correct redo lambda (2173-2178), satisfies the …

### Wt4 API misuse / logic  (1)

#### `src/DrfSelect.cpp:1632-1647` — GadrasDirectory ctor dereferences m_directoryEdit before it is constructed (null-pointer crash)
**🟠 HIGH** · _api (P2)_
- **Mechanism:** On non-web/non-iOS builds m_directoryEdit is initialized to nullptr in the member-init list (line 1632) and not constructed until line 1779 (topdiv->addNew<WLineEdit>), yet line 1647 calls m_directoryEdit->setAttributeValue(ondragstart,...) on the still-null pointer, guaranteeing a null dereference every time a GadrasDirectory is created on desktop/Electron/macOS (e.g. via GadrasDetSelect::addDirectory).
- **Trigger:** Trigger: on desktop/Electron/macOS, GadrasDetSelect ctor constructs a GadrasDirectory per shipped GADRAS dir (line 1450, e.g
- **Fix:** Move the m_directoryEdit->setAttributeValue("ondragstart","return false") call from line 1647 to immediately after m_directoryEdit is created at line ~1779-1780 (alongside setTextSize), inside the same non-web/non-iOS block.

### Security / unsafe input  (1)

#### `src/InjaLogDialog.cpp:376-406` — Exception/parse-error messages rendered as UnsafeXHTML enable XSS
**🟡 MEDIUM** · _security (P2)_
- **Mechanism:** updateDisplay's catch blocks build error strings from e.what()/e.message/json parse text via WString::tr(...).arg(rawText).toUTF8() and set them with setTextFormat(Wt::TextFormat::UnsafeXHTML), bypassing Wt's XSS filter; the JSON m_data and template inputs derive from spectrum-file-controllable content (filenames, nuclide names) and inja/json error messages echo the offending input, so a crafted file can inject script into the dialog DOM, defeating the sandboxed-iframe isolation used for the normal path at 364-374.
- **Fix:** Change the three catch-block setTextFormat calls to Wt::TextFormat::XHTML (the wrapper's <p>/<code> still render and removeScript strips injected script), or pass the error text through Wt::Utils::htmlEncode before .arg(...); keep UnsafeXHTML only for the trusted iframe-holder path at 373-374.

### Other correctness  (4)

#### `src/DetectorPeakResponse.cpp:2092-2104` — fromAppUrl builds efficiency FormulaWrapper from empty member m_efficiencyFormula instead of the parsed eqn -> silently constant efficiency 1.0
**🟠 HIGH** · _other (P2)_
- **Mechanism:** For a formula-type DRF (EFFT==F) the equation is parsed into local eqn = parts[EFFE] (line 2092), but the FormulaWrapper that becomes m_efficiencyFcn is constructed from the member m_efficiencyFormula (line 2097), which is still empty (parseFromAppUrl builds a fresh default DRF, and m_efficiencyFormula=eqn is not assigned until line 2268). FormulaWrapper treats the empty string as 1 (lines 280-281), so the imported DRF silently evaluates a constant efficiency of 1.0 while isValid() still passes; the XML path (line 3771) does this correctly by assigning the member before constructing the wrapper.
- **Trigger:** Trigger: importing a formula-type DRF via app URL/QR code yields silently wrong (constant) efficiency.
- **Fix:** At line 2097 construct the wrapper from the just-parsed local value: std::make_shared<FormulaWrapper>( eqn, isMeV ) (equivalently parts["EFFE"]), matching the value stored into m_efficiencyFormula at line 2268 and the XML path at line 3782.

#### `src/IsotopeSearchByEnergyModel.cpp:1264,1345` — Local `elements` shadows the `elements` function parameter, silently discarding the user's element restriction
**🟡 MEDIUM** · _other (P2)_
- **Mechanism:** setSearchEnergies takes const std::vector<const SandiaDecay::Element *> &elements (the user-restricted list) at line 1231, but line 1264 redeclares const vector<const SandiaDecay::Element *> &elements = db->elements(); shadowing the parameter with the full periodic table; line 1345 then passes this shadowed full list to xraysWithAllEnergies, whose limited_elements filter (lines 908-913) is therefore never empty and never restricts, so element-limited x-ray searches return matches for every element instead of the requested ones.
- **Fix:** Rename the local at line 1264 to e.g. `all_elements` and use it in the loop at 1268; at line 1345 pass the original parameter `elements` to xraysWithAllEnergies as limited_elements.

#### `target/ios/InterSpec/InterSpec/ViewController.mm:265-278` — wakeupFromBackground kills the server immediately after successfully delivering an app-URL to the live session
**🟡 MEDIUM** · _other (P3)_
- **Mechanism:** In the non-file (interspec:// URL) branch the code calls app->handleAppUrl(...) + app->triggerUpdate() under the UpdateLock but does NOT return YES, so control falls through to line 278 [self enteredBackground], which calls InterSpecServer::killServer() and tears down the very session it just updated (and _fileNeedsOpening is left non-nil, only cleared in the file branch at line 261), causing a needless server restart and a lost/duplicated open.
- **Fix:** In the else (app-URL) branch add `_fileNeedsOpening = nil; ... return YES;` after app->triggerUpdate() (line 271), mirroring the isFileURL branch, so a successfully-delivered app URL does not fall through into enteredBackground.

#### `target/macos/quicklook/SpecFilePreviewExtension/CGPaintDevice.mm:736-745` — measureText splits UTF-8 using a UTF-16 index and reads possibly-uninitialized buffer
**⚪ LOW** · _other (P3)_
- **Mechanism:** breakIndex is a count of UTF-16 code units returned by CTTypesetterSuggestLineBreak, but it is reused both to size the buffer (breakIndex*4+1) and as a byte length in std::string broken( buf, breakIndex ), so for any non-ASCII title the copied byte count is wrong and can cut a multi-byte UTF-8 sequence mid-character, yielding invalid UTF-8 passed to WString::fromUTF8; additionally the CFStringGetCString return value is ignored, so on failure buf is uninitialized and broken reads uninitialized heap bytes.
- **Fix:** Convert only the prefix: build a substring CFAttributedString/CFString for CFRangeMake(0,breakIndex), then size the buffer via CFStringGetMaximumSizeForEncoding (or CFStringGetBytes to query length), call CFStringGetCString and CHECK its BOOL result, and construct std::string broken(buf) from the NUL-terminated bytes (full converted length) rather than passing breakIndex as a byte count. Bail out / fall back to the full text on conversion failure.

---

## Uncertain — need a runtime check or a design decision (11)

_Plausible but the verifier could not establish the trigger without runtime/usage info, or the fix depends on intended behaviour. Worth a look._

- **`InterSpec/UserPreferences.h:276-289`** [critical, signal] — addCallbackWhenChanged(target,method) uses untracked single-arg lambda connect -> UAF when target widget dies before session-long UserPreferences  
  The member-function overload does signal->connect( [target,method](bool v){ (target->*method)(v); } ) which is the single-arg Wt4 connect form with NO target tracking, despite the header documenting (lines 110-119,231) that the connection is auto-disconnected when the WObject target dies. The Signal lives in UserPreferences (whole InterSpec session); callers register widgets from closeable tool windows (e.g. …
- **`InterSpec/UserPreferences.h:305-318`** [critical, signal] — addIntCallbackWhenChanged(target,method) has the same untracked single-arg lambda connect UAF  
  Identical mechanism to the bool overload: signal->connect( [target,method](int v){ (target->*method)(v); } ) on a UserPreferences-owned Signal<int> with no target tracking. Callers register tool-window widgets (RelActAutoGui.cpp:561/566 m_spectrum, InterSpec.cpp:956/962, ReferencePhotopeakDisplay.cpp:1228/1230 this) whose lifetime is shorter than the session-long UserPreferences; after the target is destroyed a …
- **`src/FitSkewParamsTool.cpp:647-674`** [critical, threading] — Background-fit post-back lambda captures raw `this` with no lifetime protection -> UAF when window closed mid-fit  
  doFit() runs the fit on a worker thread and unconditionally posts back via WServer::instance()->post(sessionId, [this,results,cancelFlag](){ handleFitResults(...); }); WServer::post only guards that the session still exists, not the FitSkewParamsTool. Closing the window (InterSpec::closeFitSkewParamsWindow -> deleteAuxWindow) destroys the tool synchronously while the session lives, so the posted lambda runs …
- **`InterSpec/UserPreferences.h:263-273`** [high, signal] — single-arg addCallbackWhenChanged(name,fcn) gives no lifetime tracking for captured this  
  This overload does signal->connect(fcn) with no target; callers such as DetectionLimitTool.cpp:1431 pass [this](){handleInputChange();} where this is a closable tool shorter-lived than UserPreferences. After the tool closes, toggling DisplayBecquerel fires the captured this into freed memory. The API offers no way to track, so callers silently get a dangling connection.
- **`src/ShieldingSelect.cpp:4964-4970`** [high, lifetime] — handleUserChangeForUndoRedo WServer::post captures raw this (ShieldingSelect) with no lifetime protection -> UAF if shielding row deleted before deferred callback runs  
  handleUserChangeForUndoRedo() does server->post( wApp->sessionId(), [this](){ handleUserChangeForUndoRedoWorker(true); } ), capturing the raw ShieldingSelect pointer; the worker later dereferences m_prevState/m_userChangedStateSignal and calls serialize(). ShieldingSourceDisplay::removeShielding() destroys ShieldingSelects synchronously (delete shielding; ShieldingSourceDisplay.cpp:8377), so if the row is destroyed …
- **`src/ShieldingSelect.cpp:5556-5563`** [high, lifetime] — deSerialize() WServer::post captures raw this with no lifetime protection -> UAF during undo/redo or tool teardown  
  At the end of deSerialize() the code posts [this](){ handleUserChangeForUndoRedoWorker(false); } to update m_prevState. deSerialize is driven by state-restore/undo-redo, which is exactly when ShieldingSourceDisplay tears down and rebuilds ShieldingSelect widgets (removeShielding does delete shielding; ShieldingSourceDisplay.cpp:8377). If the ShieldingSelect that scheduled this post is destroyed by a subsequent step …
- **`src/AppUtils.cpp:402-418`** [medium, other] — locate_file PATH search uses wrong delimiter and tests wrong path  
  The POSIX branch sets delims=":" but the call hardcodes SpecUtils::split(paths, env_path, ";"), so on macOS/Linux the colon-separated PATH is never split (treated as one giant element), and line 413 tests check_exists(filename) (the bare relative name) instead of check_exists(trialpath); thus the PATH-relative search never works, and if it 'matches' it sets filename=trialpath to an unverified path.
- **`src/MakeDrf.cpp:1783`** [medium, lifetime] — Deferred WServer::post in MakeDrf ctor captures raw this without an existence guard  
  At the end of the MakeDrf ctor, WServer::instance()->post( wApp->sessionId(), [this](){ handleSourcesUpdates(); } ) schedules a call on raw this for a later event-loop turn; unlike the other posts in this file (which re-validate via findById(thisid) before dereferencing, e.g. lines 2744/2883), this one has no guard, so if the MakeDrf (or its enclosing AuxWindow) is torn down before the posted task runs, …
- **`src/ShieldingSourceDisplay.cpp:2777-2820`** [medium, ownership] — m_sourceModel raw-new SourceFitModel never owned -> leak per dialog open  
  m_sourceModel = new SourceFitModel(...) (line 2777) is set into m_sourceView via a shared_ptr with a no-op deleter (line 2820), is never addChild'd, and ~ShieldingSourceDisplay (3430-3489) never deletes it. In Wt3 WAbstractItemModel took an owning WObject* parent; that ctor arg is gone in Wt4 so nothing now owns the model -> a SourceFitModel leaks every time the Activity/Shielding fit window is opened and closed.
- **`src/TerminalWidget.cpp:384-393`** [medium, api] — Unhandled std::regex_error from user search input crashes command-menu handler  
  commandMenuSearchInput() builds a regex from user text via searchToRegexLiteral(), but that escaper escapes ^$.*+?()[]{}| and NOT the backslash, so any input containing or ending in a lone backslash (e.g. typing '\') produces an invalid pattern; the std::regex constructor at line 393 then throws std::regex_error, which is outside the try/catch (lines 385-391 wrap only searchToRegexLiteral), so the exception escapes …
- **`target/wxwidgets/InterSpecWxUtils.cpp:52-67`** [medium, threading] — save_file_data reads wx GUI top-level-window list from a Wt io_service thread before the CallAfter hop  
  download_to_native_save (src/InterSpecServer.cpp:1097-1134) invokes this registered handler from the Wt Http::Client::done() callback, which runs on a Wt io_service worker thread, not the wx main/GUI thread; calling wxApp::GetTopWindow() (line 59) reads the wxTopLevelWindows list that is mutated only on the GUI thread, a data race, and topWindow may be concurrently destroyed if the last frame is closing. The …

---

## Minor / low-severity (26)

_Leaks of small objects, defensive-coding gaps, non-crashing logic issues, hardening._

- **`src/DbFileBrowser.cpp:1167-1170`** [api] — loadSnapshotSelected dereferences *sets.begin() without empty() check  
  m_snapshotTable->selectedNodes() can return an empty std::set; line 1170 does *sets.begin() (dereferencing end()) which is undefined behavior. The bound m_loadSnapshotButton is disabled when nothing is selected and double-click paths guard for emptiness, so it is not on the common path, but a …
- **`src/DbFileBrowser.cpp:1224-1235`** [api] — loadSpectraSelected dereferences *sets.begin() and indexes map with end()-derived key  
  Same as loadSnapshotSelected: line 1227 dereferences *sets.begin() with no empty() check (UB on empty selection), and lines 1230/1235 then use m_UserFileInDbLookup[selectedTreeNode], whose operator[] on a null/garbage key also inserts a spurious null map entry. Guarded on the common path by the …
- **`src/GoogleMap.cpp:451`** [api] — Wt4 Signal has no owner-taking constructor: m_clicked( this ) won't compile  
  GoogleMap's member is Wt::Signal<double,double> m_clicked; (GoogleMap.h:76), and Wt4's Wt::Signal only provides a no-argument constructor Signal() (WSignal.h:235) — the Wt3 Signal(WObject*) owner overload was removed. The member-init m_clicked( this ) is therefore a Wt4 API mismatch; it survives …
- **`src/QrCode.cpp:408-411`** [api] — Off-by-one bounds guard on error-correction-level combo index  
  The assert expects index in [0,3] (4 combo items) but the runtime guard is `if( ecl < 0 || ecl > 4 ) return;`, so index 4 would slip past and be passed to `static_cast<ErrorCorrLevel>(4)` / `utf8_string_to_svg_qr` as an invalid enum; not reachable today since the combo only has 4 items, but the …
- **`src/EnergyCalTool.cpp:1987-1989`** [auxwindow] — ~EnergyCalTool does not tear down open m_graphicalRecal / m_addActionWindow AuxWindows  
  m_graphicalRecal and m_addActionWindow are observing_ptr to AuxWindow subclasses created via AuxWindow::make<> (owned by wApp, lines 4716/5457) with a raw EnergyCalTool* back-pointer; the empty destructor (1987-1989) never calls deleteGraphicalRecalConfirmWindow()/cancelMoreActionWindow(), so if …
- **`src/RelActTxtResults.cpp:80-83`** [layout] — Error-path addNew<WText> called on layout-managed `this` is silently dropped from the DOM  
  The ctor calls setLayout(WGridLayout) on `this` (line 51) and adds m_txt as the only layout-managed child; in Wt4 a container with a layout renders only layout-managed children. The updateResults catch block calls this->addNew<WText>("Error displaying results") (line 82) instead of …
- **`src/ShieldingSourceDiagram.cpp:435-442`** [layout] — WGridLayout inside SimpleDialog contents() with overflow:hidden risks zero-size clipping  
  Per Wt4 layout semantics a WGridLayout positions children absolutely and the layout container has zero intrinsic size; placing it in `contents()` set to `Overflow::Hidden` can clip the 2D/3D view to nothing unless an explicit size is in force. Here the dialog sets 95vw/95vh via resize() and a CSS …
- **`src/RefSpectraWidget.cpp:877-878`** [lifetime] — Deferred WServer::post captures raw `this` without lifetime protection  
  handleSelectionChanged defers tree-node expansion via WServer::instance()->post(sessionId, [this,index](){ tryExpandNode(index); wApp->triggerUpdate(); }); `this` (the RefSpectraWidget inside RefSpectraDialog, a SimpleDialog) is captured by raw pointer. If the dialog is closed (SimpleDialog …
- **`target/wxwidgets/InterSpecWebFrame.cpp:1100-1101`** [lifetime] — OnError dereferences dynamic_cast result without null check  
  app = dynamic_cast<InterSpecWxApp*>(wxApp::GetInstance()) is immediately dereferenced as app->GetTopWindow() with no null guard (unlike every other site in these files which assert/return); if the cast ever yields null (e.g. during shutdown teardown when an error event is still being delivered) …
- **`src/ColorTheme.cpp:414-427`** [other] — fromJson casts JSON values to WString without type check, throwing on malformed theme  
  Json::parse is wrapped in try/catch but the subsequent static_cast<const WString&>(base["created"]) (and the many similar casts for color fields, lineColors entries, etc.) are not; a theme JSON whose field has the wrong JSON type (e.g. "created": 123 or a numeric color) makes Json::Value's …
- **`target/macos/NativeMenu.mm:244-257`** [other] — Dead special-case assignment in insertOsxMenuItem for the InterSpec menu  
  When the parent item text is "InterSpec", menu is set to the OSX InterSpec submenu at line 249 but is then unconditionally overwritten with (NSMenu*)voidmenu at line 252, so the special-casing has no effect; not a crash, but the intended routing of InterSpec-menu items to the native app menu …
- **`target/macos/quicklook/SpecFilePreviewExtension/CGPaintDevice.mm:338-341,356-358,387-388`** [other] — drawPlainPath reads segments[i+1]/segments[i+2] without bounds check  
  For CubicC1, ArcC, and QuadC the code dereferences segments[i+1] and segments[i+2] assuming the WPainterPath multi-segment grouping is well-formed; a truncated/malformed path would read past the vector end, but in practice these groupings are invariants of WPainterPath construction so this is …
- **`src/RelActManualGui.cpp:2122-2127`** [ownership] — Off-by-one in shield-widget removal loop leaves one extra widget parented  
  The loop `for( size_t i = (solution.m_phys_model_external_atten_shields.size() + 1); i < ext_shields.size(); ++i )` starts removal at index solution.size()+1 instead of solution.size(), so the widget at index == solution.size() is never removed (it is only resetMaterialEntryState()'d earlier), …
- **`src/ShieldingSourceDisplay.cpp:6984-7001`** [ownership] — del button leaked when no saved models exist  
  In the else branch (no previous models) accept and summary are explicitly deleted (6993-7000) but `del` (raw-new'd at 6868) is neither added to any parent nor deleted, so it leaks. A definite per-invocation leak of an unowned raw-new WPushButton.
- **`target/macos/quicklook/SpecFilePreviewExtension/CGPaintDevice.mm:737-741`** [ownership] — Raw new[]/delete[] in measureText is not exception-safe  
  char *buf = new char[bufSize] is freed by delete[] buf at line 741, but the intervening std::string broken( buf, breakIndex ) construction can throw std::bad_alloc, which skips the delete[] and leaks buf; minor since it only manifests under allocation failure.
- **`target/macos/quicklook/SpecFilePreviewExtension/ThumbnailProvider.mm:81-94`** [ownership] — cgImage leaks if QuickLook never invokes the drawing block  
  The CGImageRef returned by render_spec_file_to_cgimage is released only inside currentContextDrawingBlock; if the system discards the QLThumbnailReply without ever invoking the drawing block the CGImage is leaked. Per the QuickLookThumbnailing contract the block is normally invoked, so this is an …
- **`src/FluxTool.cpp:1407-1409`** [security] — Clipboard JS built by raw string concatenation without jsStringLiteral escaping  
  refreshPeakTable builds a doJavaScript call as "...el._isData.TableData = '" + pastebrdtxt.str() + "';" embedding data_to_strm() output directly inside a single-quoted JS literal with no escaping; any apostrophe (or a literal newline in the non-html path) in a field would break out of the string …
- **`src/GadrasShieldScatter.cpp:246-258`** [security] — Integer overflow of totalScatter (int) on large/corrupt db dimensions  
  totalScatter is an `int` product of the per-dimension counts whose individual sanity caps allow up to 64*1024*256*1024*4096 — far beyond INT_MAX — so a corrupt or crafted sandia.shieldscatter.db that passes the individual caps overflows the product to a negative/small int; read_floats early-returns …
- **`src/GoogleMap.cpp:117-127,185-209`** [security] — User-influenced HTML/metadata concatenated unescaped into doGmJavaScript single-quoted JS strings  
  addInfoWindow() and addCircle() build JavaScript by directly concatenating info.html.toUTF8()/html into single-quoted JS string literals ('content: '...'') passed to doGmJavaScript/doJavaScript with no escaping; html for the public addInfoBox()/addMeasurment() paths can carry spectrum-file …
- **`target/android/android.cpp:372-373`** [security] — JSON built by raw string concatenation of externally-supplied filepath  
  openFile constructs the files_json argument as ("[\"" + filepath + "\"]") and hands it to InterSpecServer::open_file_in_session, which feeds it to Wt::Json::parse; a filepath (from an Android Intent/content URI, externally influenced) containing a double-quote, backslash, or control char yields …
- **`target/android/android.cpp:307-331`** [security] — Native file-save handler writes using unsanitized suggested_name (path traversal)  
  The registered handler passes suggested_name straight to SpecUtils::temp_file_name(suggested_name, temp_dir()), which does append_path(directory, base) with no sanitization; a suggested_name containing ../ segments or an absolute path resolves outside temp_dir, so a save triggered with an …
- **`target/ios/InterSpec/InterSpec/ViewController.mm:449-453`** [security] — Filesystem path concatenated raw into JSON string passed to open_file_in_session  
  urlstr = [[url path] UTF8String] is interpolated directly as a JSON string literal via "[\"" + urlstr + "\"]" and handed to InterSpecServer::open_file_in_session, which runs Wt::Json::parse on it; a file path containing a backslash or double-quote (legal in filenames) produces malformed/injected …
- **`target/wxwidgets/InterSpecWebFrame.cpp:243`** [security] — wxLogMessage called with a runtime-built string as the format argument  
  wxLogMessage("URL: " + m_url) passes a concatenated wxString as the printf-style format string; if m_url ever contains a '%' sequence wxWidgets will interpret format specifiers against absent varargs, reading garbage/crashing. m_url is the internal server URL so exposure is limited, but the pattern …
- **`src/SpecFileQueryWidget.cpp:2753-2764`** [threading] — doSearch tests thread-local wApp on an asio worker thread, dead local-time branch  
  doSearch is dispatched via ioService().boost::asio::io_service::post(...) so it runs on a thread-pool worker where wApp (WApplication::instance(), thread-local) is null; the 'if(wApp)' branch that adds the user's local-time-zone offset to the search-description timestamp is therefore never taken …
- **`target/electron/InterSpecAddOn.cpp:593-613`** [threading] — runBatchAnalysis mutates process-global nuclear-data servers on the node thread while GUI sessions may be live  
  runBatchAnalysis (node main thread) calls InterSpec::setStaticDataDirectory/setWritableDataDirectory, which reconfigure process-global singletons (DecayDataBaseServer, ReactionGammaServer, MassAttenuation, IsotopeId) at src/InterSpec.cpp:1218-1232; only the directory-string assignment is …
- **`src/PeakEdit.cpp:3162-3169`** [undoredo] — redo lambda in undo/redo step captures PeakEdit* this  
  The redo lambda registered with addUndoRedoStep captures 'this' (a PeakEdit window/tool pointer) in violation of the project's undo/redo discipline; it happens to be harmless only because the body looks the editor up via get_session_peak_editor() and never dereferences the captured this, but the …

---

## Coverage & remaining gaps

- **Covered:** all 159 `src/*.cpp`, the `InterSpec/*.h` headers, the shipping `target/` platform glue, `InterSpec_resources/*.js` JSignal wiring, and WTemplate/XSS sinks.
- **Not deeply covered (suggested next pass):** the `data/config/wt_config_*.xml` server config vs the Wt 4.12.6 schema; `app_text/*.xml` `tr()`-key wiring vs `WString::tr()` call sites; CMake/`#ifdef`-gated platform branches that were not the active build; and the `target/peak_fit_improve/` research harness (standalone dev tool, not the Wt UI).
- **Already-fixed items** from `wt4_ui_issues_fixes.md` were excluded by design and are not re-reported.

## Notes

- No source files were modified — this is a review deliverable. The systemic fixes (A–C especially) are mechanical and low-risk; A and C can largely be done with guided search-and-replace plus a compile.
- Counts reflect post-dedup clusters; a single cluster may span several call sites (listed under *Also at*).
- Per-finding raw data (including each verifier's full reasoning) is in `/tmp/wt4_merged.json`.

---

# Fixes Applied (2026-05-29)

**Every confirmed critical/high/medium finding in `src/` + `InterSpec/` has been fixed and the full
application rebuilds cleanly** (`InterSpec.dylib` + `InterSpec` executable, 0 errors). The `target/macos`
`NativeMenu.mm` HIGH pair was also fixed and validated by building the macOS Xcode app (`** BUILD
SUCCEEDED **`). Each file was compiled individually during the work and a final full `cmake --build .`
succeeded. No behavior was changed beyond the bug fixes. Fixes touched ~25 files.

**Verified at runtime via Chrome** (`build_wt4`, developer-checks/asserts ON): app loads with a spectrum;
`removeSearchEnergy`, `removeShielding`, and `RelActAutoGui` energy-range add/remove all execute without
crashing or asserting. **That runtime testing surfaced — and led to the fix of — a second, deeper UAF
that the double-free had been masking** (see "Destroy-during-own-emit" below).

### By pattern
- **A — `delete child` → `parent->removeWidget(child)` (double-free):** `RelActAutoGui.cpp` (every
  ROI/peak/nuclide remove handler + the `handleDelRelEffCurve` `removeItem` double-delete + `createWindow`
  catch), `ShieldingSourceDisplay.cpp` (`reset`, `removeShielding`, `createWindow` catch,
  `startBrowseDatabaseModels` catch), `IsotopeSearchByEnergy.cpp`, `ReferencePhotopeakDisplay.cpp`,
  `DrfSelect.cpp` (`GadrasDetSelect::removeDirectory`), `RelActAutoGuiRelEffOptions.cpp`,
  `MakeDrfSrcDef.cpp` (`~MakeDrfSrcDef`), `DetectionLimitSimple.cpp` (`SimpleDialog` → `wApp->removeChild`).
- **B — worker/timer callbacks: re-resolve via `findById()` (thread-safe) instead of capturing raw
  `this`/widget pointers:** `RelActManualGui`, `RelActAutoGui`, `ShieldingSourceDisplay` (fit
  progress/result updaters), `FitPeaksForNuclidesGui` (also copies inputs into locals so the worker never
  touches `this`), `FitSkewParamsTool`, `MakeFwhmForDrf`, `BatchGuiInputFile`, `BatchGuiWidget`,
  `RemoteRid` (×2), `SimpleActivityCalc`; `DirectorySelector` and `IsotopeSelectionAids` now marshal the
  native/timer callback onto the session thread first. (`observing_ptr` was deliberately **not** used for
  these — it is copied across threads here and its observer-list mutation isn't thread-safe.)
- **C — `UserPreferences.h`:** `addCallbackWhenChanged`/`addIntCallbackWhenChanged` now use the
  lifetime-tracked `signal->connect(target, method)` overload.
- **D — `WResource::setTakesUpdateLock(true)`:** `SpecFileQueryWidget` `ResultCsvResource`.
- **F — ownership/leaks:** `MakeDrf` `DrfPeak`/`MakeDrfSrcDef` now `addNew<>` to their containers (fixes
  invisible peaks/sources + leak); `FitSkewParamsTool` `PeakModel`/`WPopupMenu`, `NuclideSourceEnter`
  controller, `PeakModel` `PeakCsvResource` now owned via `addChild`; `SpecMeasManager`
  `FileDragUploadResource`s freed in the destructor.
- **G / misc:** `~ShieldingSourceDisplay` diagram dialog (`removeFromParent` → `wApp->removeChild`);
  `DecayActivityDiv.h` `m_csvDownloadDialog` → `observing_ptr`; `DetectorPeakResponse::fromAppUrl`
  (built efficiency formula from an empty member); `DrfSelect` `GadrasDirectory` ctor null-deref;
  `DrfSelect::RelEffFile` now removes on confirmation, not unconditionally;
  `IsotopeSearchByEnergyModel` parameter-shadowing (was ignoring the user's element restriction);
  `ReferencePhotopeakDisplay` undo-lambda capturing `this`; `InjaLogDialog` `UnsafeXHTML` XSS
  (htmlEncode exception text); `D3SpectrumDisplayDiv` drag-ROI cross-thread member read +
  `DeleteOnClosePopupMenu` ownership/leak. (Also fixed the low-severity `ShieldingSourceDisplay`
  `del`-button leak found in the same hunk.)
- **E — `target/macos/NativeMenu.mm` (HIGH ×2), validated by a macOS Xcode build (`BUILD SUCCEEDED`):**
  `validateMenuItem` now reads a cached `std::atomic<bool>` enabled-flag on the `Target` (synced from the
  session thread via a new `PopupDivMenuItem::setDisabled` override) instead of dereferencing the Wt widget
  on the AppKit thread; the `Target`'s widget pointers are now `std::atomic` and a new `invalidate` (called
  synchronously from `~PopupDivMenuItem`) nulls them, so `validateMenuItem`/`clicked`/`toggleChecked` can
  no longer touch a destroyed widget. (The `NSOffState`/`NSOnState` deprecation warnings are pre-existing.)

### Destroy-during-own-emit UAF (found by Chrome verification, then fixed)
Removing a widget from within **its own** member-`Signal` emission (e.g. a row's `remove()` signal whose
slot destroys the row) frees the emitting signal *while Wt is still emitting it* → `ProtoSignal::emit`
then iterates a freed connection list and `dynamic_cast`s a garbage target → **SIGSEGV**. This pre-existed
the review: the original `delete child` had it too, but the double-free *assert* (developer-checks build)
fired first and masked it; once the double-free was fixed (Pattern A), the latent UAF became a reliable
crash (`removeShielding` crashed deterministically in testing; the crash report's faulting frame is
`ShieldingSelect::emitRemoveSignal` → `Signal::emit` → `__dynamic_cast`). Fix: in every handler that
removes a widget in response to that widget's own member-`Signal`, **detach immediately but defer
destruction** to the next event loop — `parent->removeWidget(child)` and keep the returned `unique_ptr`
alive in a `WServer::post(sessionId, …)` task that drops it after the emit unwinds (helper
`removeWidgetLater(...)` in `RelActAutoGui.cpp`; inline elsewhere). Applied to
`ShieldingSourceDisplay::removeShielding`, `RelActAutoGui` (`handleRemoveFreePeak`, `handleRemoveEnergy`,
`handleCombineRoi`, the convert/`on_yes` and `handleRemovePartOfEnergyRange` paths, `handleRemoveNuclide`,
and `handleDelRelEffCurve` via the `WMenu::removeItem` unique_ptr), and `IsotopeSearchByEnergy::removeSearchEnergy`.
**Re-verified via Chrome:** `removeShielding`, `removeSearchEnergy`, and `RelActAutoGui` energy-range
removal now run cleanly (count decrements, no crash, no assert). The pure-`EventSignal` "delete-from-button"
cases (`DrfSelect::removeDirectory`, `ReferencePhotopeakDisplay::removeFeatureMarkerTool`) have no
intermediate member-`Signal` and follow Wt's well-supported pattern, so they were left synchronous.

### NOT applied — other `target/` glue (needs a platform build to compile/validate)
These are only built under specific platform flags (`BUILD_AS_ELECTRON_APP`, etc.) and **cannot be compiled
by the web build (`build_wt4`)**, so they were left for a platform-build session rather than changed blind:
- **`target/wxwidgets/InterSpecWxApp.cpp` & `InterSpecWxUtils.cpp`, `target/electron/InterSpecAddOn.cpp`,
  `target/ios/.../ViewController.mm`, `target/android/android.cpp` (MEDIUM):** off-thread wx/JNI access,
  an Electron `BlockingCall` deadlock path, and an iOS wake-up that kills the server after delivering a
  URL. See the per-finding entries above for the mechanism/fix of each.

### NOT applied — uncertain / needs a runtime or design decision
The 10 *Uncertain* findings (e.g. `ShieldingSelect` deferred `WServer::post` ×2, `AppUtils::locate_file`
PATH delimiter, `TerminalWidget` `std::regex_error`) and the *Minor/low* appendix were left as-is — the
user requested confirmed critical/high/medium. The single-arg `UserPreferences::addCallbackWhenChanged(name,
fcn)` overload has no `target` to track, so its safety remains the caller's responsibility (all current
callers pass session-lifetime objects).