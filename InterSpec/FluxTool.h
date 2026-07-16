#ifndef FluxTool_h
#define FluxTool_h
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

#include <array>
#include <string>
#include <vector>
#include <memory>

#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

//  ToDo:
//    - Better and more-consistent printing to appropriate number of significant figures.
//    - Havent fully tested that copying to clipboard will work everywhere.
//    - Can maybe improve copying to clipboard using the clipboard API.
//    - Have some capability to automatically fit for a number of pre-defined
//      peaks for easier batch processing of like Co60 density measurements.

//I think copying to the clipboard is working well, but leaving this optional
//  for the moment until I do some more testing.
#define FLUX_USE_COPY_TO_CLIPBOARD 1

class PeakDef;
class InterSpec;
class FluxToolWidget;
class DetectorDisplay;
class RowStretchTreeView;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WCheckBox;
  class WLineEdit;
#if( FLUX_USE_COPY_TO_CLIPBOARD )
  class WPushButton;
#endif
  class WButtonGroup;
}//namespace Wt

namespace FluxToolImp
{
  class FluxModel;
  class FluxCsvResource;
}//namespace FluxToolImp


class FluxToolWindow : public AuxWindow
{
public:
  FluxToolWindow( InterSpec *viewer );
  
  virtual ~FluxToolWindow();
  
  /** See #FluxToolWidget::handleAppUrl */
  void handleAppUrl( const std::string &query_str );
  
  /** See #FluxToolWidget::encodeStateToUrl */
  std::string encodeStateToUrl() const;
  
protected:
  FluxToolWidget *m_fluxTool;
  
  friend class FluxToolWidget;
};//class FluxToolWindow


class FluxToolWidget : public Wt::WContainerWidget
{
public:
  enum class DisplayInfoLevel
  {
    /** Only energy, nuclide, and gammas into 4pi are shown.
     Gammas into 4pi uncertainty gets its own column in CSV, and is given as a percent uncertainty.
     */
    Simple,
    
    /** Nuclide, IntrinsicEff, GeometricEff, FluxOnDet columns are NOT shown.
     Uncertainties get own column in CSV, as actual value (e.g., not percent).
     */
    Normal,
    
    /** All columns are shown.
     Uncertainties are placed in their own column in CSV, as the actual value (e.g., not percent).
     */
    Extended
  };//enum class DisplayInfoLevel
  
  
public:
  FluxToolWidget( InterSpec *viewer,
                  Wt::WContainerWidget *parent = 0 );
  
  
  virtual ~FluxToolWidget();
  
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  enum FluxColumns
  {
    FluxEnergyCol,
    FluxNuclideCol,
    FluxPeakCpsCol,
    FluxIntrinsicEffCol,
    FluxGeometricEffCol,
    FluxFluxOnDetCol,
    FluxFluxPerCm2PerSCol,
    FluxGammasInto4PiCol,
    FluxNumColumns
  };//enum FluxColumns

  /** The model column holding the per-row "include in output" checkbox.

   Deliberately placed after all of #FluxColumns, so the data columns, the CSV/clipboard
   generation, and the sorting comparator (which index #m_data by column) need no changes.
   The model always reports this column - see #FluxToolWidget::setSelectRows for why it is
   hidden rather than removed when rows are not being selected.
   */
  static const int sm_selectColumn = FluxColumns::FluxNumColumns;
  
  Wt::Signal<> &tableUpdated();
  
  DisplayInfoLevel displayInfoLevel() const;
  
  
  /** Handles receiving a "deep-link" url starting with "interspec://flux?dist=1.2m&display=low".

   Example URIs:
   - "interspec://flux?VER=1&dist=1.2m&display=low"
   - "interspec://flux?VER=1.1&dist=1.2m&display=low&selrows=1&sel=661.66,1173.23"

   The optional "SELROWS" and "SEL" fields (added in version 1.1) give the row-selection
   state; both are optional, so URIs written by previous versions still load.  "SELROWS=1"
   means rows are being individually selected.  "SEL" lists the peak means (in keV) of the
   selected rows; each is matched to the nearest peak within 0.5 FWHM of it.  An absent
   "SEL" means every peak is selected, while an empty one ("SEL=") means none are.

   Note that "SEL" is resolved against the peaks loaded when the URI is parsed, since we
   track the *un*-selected rows (so that a peak fit later defaults to being selected).  So
   for a URI applied to a session with no peaks - a QR code scanned before a spectrum is
   loaded, say - the row selection is dropped, and everything gets output.  Restoring a
   saved state is not affected: #InterSpec::loadStateFromDb sets the spectrum, and hence
   the peaks, before it applies this URI.

   @param query_str The query portion of the URI.  So for example, if the URI has a value of
          "interspec://flux?dist=1.2m&display=low", then this string would be "dist=1.2m&display=low".
          This string is in standard URL format of "key1=value1&key2=value2&..." with ordering not mattering.
          Capitalization is not important.
          Assumes the string passed in has already been url-decoded.
          If not a valid query_str, throws exception.
   */
  void handleAppUrl( std::string query_str );

  /** Encodes current tool state to app-url format.  Returned string is just query portion of of URL,
   so will look something like "dist=1.2m&display=low", and it will not be url-encoded.

   Fields that are at their default are left out, to keep the URI (and hence the QR code)
   short - see #handleAppUrl for what an absent field means.
   */
  std::string encodeStateToUrl() const;
  
protected:
  void init();
  void distanceUpdated();
  void setTableNeedsUpdating();
  void refreshPeakTable();
  void handleDrfChange( std::shared_ptr<DetectorPeakResponse> drf );


  void setDisplayInfoLevel( const DisplayInfoLevel disptype, const bool force );

  /** Turns the per-row selection checkbox column on or off, adding an undo/redo step.

   Does not touch which rows are selected, so un-ticking and re-ticking gets the user back
   the rows they had picked.

   @param selectRows Whether rows should be individually selectable.
   */
  void setSelectRows( const bool selectRows );

  /** Handles the user checking/un-checking a single rows checkbox, adding an undo/redo step.

   @param dataRow Index into #m_data (i.e., <em>not</em> the displayed row, which sorting
          may have moved).
   */
  void handleRowSelectionChanged( const size_t dataRow, const bool selected );

  /** Applies a complete selection state, and updates the GUI to match.

   Deliberately does *not* add an undo/redo step - it is what the undo/redo steps call, and
   is also used from #handleAppUrl.  See #setSelectRows and #handleRowSelectionChanged for
   the user-driven entry points that do record steps.

   \p deselected is taken by value since #setSelectRows passes #m_deselectedEnergies itself.
   */
  void applySelectionState( const bool selectRows, std::vector<double> deselected );

  /** Re-derives #m_rowSelected from #m_deselectedEnergies, matching remembered energies up
   to the passed in peaks (which #m_data was just built from).

   Pure - never adds an undo/redo step, since it is called from #refreshPeakTable, which
   runs during rendering.
   */
  void syncRowSelectionFromEnergies( const std::vector<PeakDef> &peaks );

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  void tableCopiedToCliboardCallback( const int copied );

  /** Regenerates the HTML text the "Copy To Clipboard" button places on the clipboard.

   Separate from #refreshPeakTable so a selection change can update the clipboard without
   emitting #m_tableUpdated, which would reset the model and re-render the whole table.
   */
  void updateCopyToClipboardText();
#endif

  InterSpec *m_interspec;
  DetectorDisplay *m_detector;
  
  /** Wether layout out for portrait phones or not. */
  bool m_narrowLayout;
  
  Wt::WText *m_msg;
  Wt::WLineEdit *m_distance;
  Wt::WString m_prevDistance; // For undo/redo
  RowStretchTreeView *m_table;
  FluxToolImp::FluxModel *m_fluxModel;

  /** Lets the user pick which rows get copied to the clipboard, or exported to CSV. */
  Wt::WCheckBox *m_selectRowsCb;

  /** Whether the per-row checkbox column is being shown, and honored by the CSV/clipboard
   output.  When false, all rows are output.
   */
  bool m_selectRows;

  /** The peak means (in keV) of the rows the user has un-checked.

   Kept keyed by energy, rather than by row index, so selections survive the table being
   rebuilt (peak re-fits, peaks added/removed, a transiently invalid distance or DRF, ...);
   #refreshPeakTable matches these back up to peaks within 0.5 FWHM.

   Tracking the *un*-selected rows (rather than the selected ones) is what makes a peak fit
   later default to checked, and keeps this independent of #m_selectRows, so un-ticking
   "select rows" and re-ticking it remembers what the user had.

   Only ever written by a user action (a checkbox, an undo/redo, or an app-URL) - never
   derived from, or pruned against, the current table.  See #refreshPeakTable, which clears
   #m_data and then early-returns on several common conditions; re-deriving this member
   there would silently discard the users selections.
   */
  std::vector<double> m_deselectedEnergies;

  /** Whether each row of #m_data is checked; derived from #m_deselectedEnergies by
   #syncRowSelectionFromEnergies.  Always either empty, or the same size as #m_data, and in
   the same (unsorted) order.
   */
  std::vector<bool> m_rowSelected;

#if( FLUX_USE_COPY_TO_CLIPBOARD )
  Wt::WPushButton *m_copyBtn;
  Wt::JSignal<int> m_infoCopied;
#endif

  /** We will only update dm_data and m_uncertainties from peak model right
      before rendering happens to avoid duplicate work if multiple peaks are
      being added.
   */
  bool m_needsTableRefresh;
  
  /** What columns to show. */
  DisplayInfoLevel m_displayInfoLevel;
  Wt::WButtonGroup *m_displayLevelButtons;
  
  Wt::Signal<> m_tableUpdated;
  
  std::array<Wt::WString,FluxColumns::FluxNumColumns> m_colnames;
  std::array<Wt::WString,FluxColumns::FluxNumColumns> m_colnamesCsv;

  /** Header for #sm_selectColumn.  Kept separate from #m_colnames since that array is
   sized to the data columns, and the selection checkbox never appears in the CSV.
   */
  Wt::WString m_selectColName;
  
  std::vector<std::string> m_nucNames;
  std::vector<std::array<double,FluxColumns::FluxNumColumns>> m_data;
  std::vector<std::array<double,FluxColumns::FluxNumColumns>> m_uncertainties;
  
  friend class FluxToolWindow;
  friend class FluxToolImp::FluxModel;
  friend class FluxToolImp::FluxCsvResource;
};//class FluxToolWidget

#endif //FluxTool_h

