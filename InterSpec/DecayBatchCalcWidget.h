#ifndef DecayBatchCalcWidget_h
#define DecayBatchCalcWidget_h
/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"

#include <string>
#include <vector>
#include <memory>
#include <utility>

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>
#include <Wt/Core/observing_ptr.hpp>

#include "InterSpec/GroupBox.h"
#include "InterSpec/AuxWindow.h"

// The non-GUI core (decay math + CSV I/O) lives in the DecayBatchCalc namespace.
#include "InterSpec/DecayBatchCalc.h"

class InterSpec;
class FileDragUploadResource;
class DecayBatchCalcWidget;

namespace SandiaDecay
{
  struct Nuclide;
}

namespace Wt
{
  class WLineEdit;
  class WComboBox;
  class WSpinBox;
  class WCheckBox;
  class WAnchor;
  class WPushButton;
  class WSuggestionPopup;
}//namespace Wt


/** A single initial-nuclide input row: nuclide (name-suggest), initial age, initial activity, and a
 remove button.  All three text inputs carry validators.  Child widths are fixed so the row does not
 jitter when the age slot is enabled/disabled for (non-)ageable nuclides.
 */
class DecayBatchCalcNuclide : public Wt::WContainerWidget
{
public:
  DecayBatchCalcNuclide( Wt::WSuggestionPopup *nuclideSuggest );

  /** Detaches this row's nuclide edit from the shared suggestion popup (see #m_nuclideSuggest). */
  virtual ~DecayBatchCalcNuclide();

  /** The currently entered nuclide, or nullptr if the text does not resolve. */
  const SandiaDecay::Nuclide *nuclide() const;

  /** Initial age in SandiaDecay units (seconds); 0 if blank/invalid or not ageable. */
  double age() const;

  /** Initial activity in SandiaDecay units (becquerel); 0 if blank/invalid.  Accepts either a value
   with activity units (e.g. "5 uCi") or a bare scalar (e.g. "5", interpreted as becquerel). */
  double activity() const;

  /** Whether the nuclide and activity inputs both currently resolve.  A zero activity *is* valid: it
   is a legitimate "none detected" measurement (and grouped input depends on those rows surviving). */
  bool isValid() const;

  /** Whether the activity was entered with explicit units (vs. a bare scalar). */
  bool activityHasUnit() const;

  /** The activity's unit token (e.g. "uCi"), or empty if a bare scalar / unparseable. */
  std::string activityUnitStr() const;

  /** Populates the row from resolved values (used when seeding or loading a CSV/URL). */
  void setNuclide( const SandiaDecay::Nuclide *nuc, double activity,
                   bool useCurie, double age, const std::string &activityStr );

  std::string nuclideText() const;
  std::string activityText() const;
  std::string ageText() const;

  /** Sets the activity edit text programmatically and updates the last-good value (used to harmonize
   units across rows).  Does not emit #changed. */
  void setActivityText( const std::string &txt );

  /** Opaque areal/volumetric unit label (e.g. "/m2") parsed from a Product/Value/Unit CSV; carried
   through to the result labels.  Empty for normal input.  Reset by #setNuclide. */
  void setUnitLabel( const std::string &label );
  std::string unitLabel() const;

  /** Emitted when any of the row inputs change. */
  Wt::Signal<> &changed();

  /** Emitted when the user clicks this row's "-" button. */
  Wt::Signal<DecayBatchCalcNuclide *> &removeRequested();

protected:
  void handleNuclideChange();
  void handleAgeChange();
  void handleActivityChange();
  void updateAgeEnabledState();

  /** The shared suggestion popup (owned by the parent DecayBatchCalcWidget).  Held as an
   observing_ptr so that if the popup is torn down before this row, the destructor's removeEdit() is
   safely skipped. */
  Wt::Core::observing_ptr<Wt::WSuggestionPopup> m_nuclideSuggest;

  Wt::WLineEdit *m_nuclideEdit;
  Wt::WContainerWidget *m_ageContainer;
  Wt::WLineEdit *m_ageEdit;
  Wt::WLineEdit *m_activityEdit;
  Wt::WPushButton *m_removeBtn;

  std::string m_unitLabel;

  /** Last accepted age / activity text; an invalid edit reverts to these (the visual feedback). */
  std::string m_prevAgeText;
  std::string m_prevActivityText;

  Wt::Signal<> m_changed;
  Wt::Signal<DecayBatchCalcNuclide *> m_remove;
};//class DecayBatchCalcNuclide


/** One group of input rows: the nuclides measured at a single physical location.

 Grouping only ever comes from an uploaded CSV that names locations (see `BatchNuclide::location`);
 there is no GUI control to create one.  So the ordinary single-set case is a *single* group with an
 empty title, which renders with no legend or border (the "Untitled" style class) and is therefore
 invisible to a user who never uploads such a file.

 Rows live in a child container, so index-based traversal of the rows is never off-by-one from
 `GroupBox`'s legend being child 0.
 */
class DecayBatchCalcLocation : public GroupBox
{
public:
  /** `location` is the physical location name; empty for the untitled/ungrouped case. */
  DecayBatchCalcLocation( const std::string &location );

  /** The location name; empty when this group is the untitled single group. */
  const std::string &location() const;

  /** Sets the location name, updating the legend and the untitled styling. */
  void setLocation( const std::string &location );

  /** The source cells this location's output rows echo verbatim: the text following the nuclide in
   the source "Product" column, the extra columns (e.g. Latitude/Longitude), and the activity unit
   token the values were given in.  See `DecayBatchCalc::BatchNuclide`, which documents them; they
   belong to the location rather than a row, since the core reads one set per group.  All empty for
   the untitled/ungrouped case.
   */
  void setPassthrough( const std::string &product_suffix,
                       const std::vector<std::pair<std::string,std::string>> &extra_columns,
                       const std::string &activity_unit );
  const std::string &productSuffix() const;
  const std::vector<std::pair<std::string,std::string>> &extraColumns() const;
  const std::string &activityUnit() const;

  int numRows() const;
  DecayBatchCalcNuclide *row( int index ) const;

  /** Appends a row; the caller connects its signals. */
  DecayBatchCalcNuclide *addRow( Wt::WSuggestionPopup *nuclideSuggest );

  /** Detaches `row` (destroying it after the current event unwinds - it is normally removed from
   within its own signal).  Returns the number of rows left. */
  int removeRow( DecayBatchCalcNuclide *row );

protected:
  std::string m_location;
  std::string m_productSuffix;
  std::vector<std::pair<std::string,std::string>> m_extraColumns;
  std::string m_activityUnit;
  Wt::WContainerWidget *m_rows;
};//class DecayBatchCalcLocation


/** The embedded "Batch Decay" tool: enter (or drag/drop) many initial nuclides, pick decay options,
 view the result table, and export it to CSV or the clipboard.

 This widget owns its Export/Copy controls (they are children, so they work no matter how the widget
 is hosted).  State is URI-representable (deep link / QR), takes part in undo/redo, and is saved in
 app state (see `InterSpec::createDecayBatchCalcWindow` / `UserState::decayBatchUri`).
 */
class DecayBatchCalcWidget : public Wt::WContainerWidget
{
public:
  DecayBatchCalcWidget( InterSpec *viewer );
  virtual ~DecayBatchCalcWidget();

  /** Adds an initial nuclide input row (signature parallels DecayActivityDiv::addNuclide). */
  void addNuclide( const int z, const int a, const int iso,
                   const double activity, const bool useCurie,
                   const double age, const std::string &activityStr );

  /** Removes all input rows. */
  void clearNuclides();

  /** Encodes tool state to app-url query format (path is "calc"); not url-encoded. */
  std::string encodeStateToUrl() const;

  /** Restores tool state from an app-url path + query (as produced by #encodeStateToUrl). */
  void handleAppUrl( std::string path, std::string query_str );

  /** CSV text of the most recently computed result (used by the download resource). */
  std::string currentResultCsv() const;

protected:
  void addEmptyNuclideRow();
  void handleRemoveRow( DecayBatchCalcNuclide *row );
  void scheduleResultUpdate();
  void updateResult();

  /** The location groups, in display order; always at least one (see #DecayBatchCalcLocation). */
  std::vector<DecayBatchCalcLocation *> locationGroups() const;

  /** Every input row of every group, in display order. */
  std::vector<DecayBatchCalcNuclide *> allRows() const;

  /** The group `row` belongs to, or nullptr if it is not one of ours. */
  DecayBatchCalcLocation *groupOf( DecayBatchCalcNuclide *row ) const;

  /** Removes all groups (and so all rows) and returns a fresh, empty, untitled group.  Used when
   (re)loading from a CSV or URL; the caller is responsible for there ending up being >=1 row. */
  DecayBatchCalcLocation *resetToSingleGroup();

  /** Adds a group with the given location name (empty for untitled). */
  DecayBatchCalcLocation *addLocationGroup( const std::string &location );

  /** Appends a blank, connected row to `group`, seeded with #defaultNewActivityText. */
  DecayBatchCalcNuclide *addRowToGroup( DecayBatchCalcLocation *group );

  /** Connects a freshly added row's signals to this widget. */
  void connectRow( DecayBatchCalcNuclide *row );

  /** Whether any group carries a location name, i.e. grouped (multi-location) input is loaded. */
  bool hasLocationGroups() const;

  /** Reflects #hasLocationGroups in the "Mix inputs" checkbox: grouped input always co-decays within
   a location and never across locations, so the checkbox is forced on and disabled. */
  void updateMixInputState();

  /** Default activity text for a newly added row: "1", or "1 <unit>" if any existing row carries an
   activity unit (so the new row joins the same units style). */
  std::string defaultNewActivityText() const;

  /** Enforces a single activity-units style across rows: if some rows are unit-bearing and others are
   bare scalars, each scalar is given the unit of its nearest unit-bearing sibling (reflected in the
   GUI), keeping cross-nuclide normalization consistent. */
  void harmonizeActivityUnits();

  /** Gathers the current, valid input rows into core BatchNuclide structs. */
  std::vector<DecayBatchCalc::BatchNuclide> gatherInputs() const;

  /** Reads the minimum-activity field.  Returns false when the text is present but not a usable
   non-negative number, in which case `cut` is zero; a blank field is valid and gives a cut of zero
   (i.e. no cut).  #updateResult reports an invalid cut rather than computing without one.
   */
  bool minActivityCut( double &cut ) const;

  /** Reads options widgets into a core BatchDecayOptions struct.  An invalid minimum-activity is read
   as no cut here; #updateResult checks it separately (see #minActivityCut) so it can say so. */
  DecayBatchCalc::BatchDecayOptions gatherOptions() const;

  void handleFileDrop( const std::string &display_name, const std::string &spool_name );
  void loadCsvContents( const std::string &contents );

  void tableCopiedToClipboardCallback( const int copied );
  void updateCopyToClipboardText();

  /** Inserts an undo/redo step if the tool state (nuclides, ages, activities, options) has changed
   since the last recorded state.  Captures the URI-encoded state (see #encodeStateToUrl); the
   undo/redo lambdas re-look-up the tool through InterSpec and restore via #handleAppUrl, so they are
   safe across the tool being closed and reopened.  A no-op while a step is being restored (undo/redo)
   or before the tool has finished constructing. */
  void handleAddUndoPoint();

  InterSpec *m_interspec;

  /** Holds the DecayBatchCalcLocation groups (never the rows directly). */
  Wt::WContainerWidget *m_nuclideRows;
  Wt::WPushButton *m_addRowBtn;
  Wt::WSuggestionPopup *m_nuclideSuggest;

  Wt::WLineEdit *m_timeEdit;
  Wt::WSpinBox *m_stepsSpin;
  Wt::WLineEdit *m_minActEdit;
  Wt::WComboBox *m_actUnitsCombo;
  Wt::WCheckBox *m_mixInput;
  Wt::WCheckBox *m_showProgeny;
  Wt::WCheckBox *m_incActivity;
  Wt::WCheckBox *m_incXrays;
  Wt::WCheckBox *m_incGammas;
  Wt::WCheckBox *m_incAlphas;
  Wt::WCheckBox *m_incBetas;

  Wt::WContainerWidget *m_dropArea;
  std::unique_ptr<FileDragUploadResource> m_uploadResource;

  Wt::WContainerWidget *m_resultContainer;
  Wt::WAnchor *m_csvDownload;
  Wt::WPushButton *m_copyBtn;
  Wt::JSignal<int> m_infoCopied;

  /** The most recently computed result; the CSV resource and clipboard copy read from this. */
  DecayBatchCalc::BatchDecayResult m_lastResult;

  /** When true, #scheduleResultUpdate is a no-op; used to coalesce the many row additions done by
   #loadCsvContents / #handleAppUrl into a single recompute at the end. */
  bool m_suppressUpdate;

  /** The URI-encoded tool state as of the last recorded undo/redo point; #handleAddUndoPoint diffs
   the current state against this to decide whether a change occurred. */
  std::string m_previousStateUri;

  /** False until the constructor finishes, so the initial widget build (and any nuclides seeded right
   after construction) does not generate undo/redo steps. */
  bool m_trackUndo;
};//class DecayBatchCalcWidget


/** Thin AuxWindow host for a DecayBatchCalcWidget; adds help + close + QR footer buttons.  Modeled
 on DecayWindow.  Create via AuxWindow::make<DecayBatchCalcWindow>().
 */
class DecayBatchCalcWindow : public AuxWindow
{
  friend class AuxWindow;

public:
  virtual ~DecayBatchCalcWindow();

  void addNuclide( const int z, const int a, const int iso,
                   const double activity, const bool useCurie,
                   const double age, const std::string &activityStr );
  void clearNuclides();

  void handleAppUrl( const std::string &url );
  void handleAppUrl( const std::string &path, const std::string &query_str );
  std::string encodeStateToUrl() const;

protected:
  DecayBatchCalcWindow( InterSpec *viewer );

  DecayBatchCalcWidget *m_calc;
};//class DecayBatchCalcWindow

#endif //DecayBatchCalcWidget_h
