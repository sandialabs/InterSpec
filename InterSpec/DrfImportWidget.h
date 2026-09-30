#ifndef DrfImportWidget_h
#define DrfImportWidget_h
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

#include <memory>
#include <string>
#include <vector>
#include <cstdint>

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/DrfImport.h"

class SimpleDialog;
class EccUncertOptions;
class DetectorPeakResponse;
class FileDragUploadResource;

namespace Wt
{
  class WText;
  class WLineEdit;
  class WComboBox;
  class WPushButton;
}//namespace Wt


/** Imports a detector response function (DRF) from any of the supported file types (see
 `DrfImport::FileKind`): identifies the file, then shows the options that file type has - a drop
 area for its companion file (a Detector.dat for an Efficiency.csv, a DETECTOR.txt for a .par, and
 vice-versa), how to interpret its efficiencies, the detector diameter, and so on.

 While its drop areas are visible, the app's own file drag-n-drop is suppressed (see
 `dropAreaShowing` in src/js_inline/InterSpec.js), so a file dropped anywhere in the window is not
 opened as a spectrum.
 */
class DrfImportWidget : public Wt::WContainerWidget
{
public:
  enum class Host
  {
    /** The "Import" tab of `DrfSelect`: has a drop area (or click to browse) to choose the file.
     DrfSelect applies every rebuilt DRF app-wide (with an undo step), so the diameter, distance
     and name fields rebuild once an edit is committed, not on every keystroke.
     */
    ImportTab,

    /** The dialog for a DRF file dropped onto the app (see #setupDropDialog): the file is already
     chosen, so only its companion file is asked for, and the fields rebuild on every keystroke so
     "Use" is enabled the moment the DRF is complete.
     */
    DropDialog
  };//enum class Host

  explicit DrfImportWidget( const Host host );
  virtual ~DrfImportWidget();

  /** Adds a file.  One half of the current pair (e.g., its Detector.dat) replaces that half; any
   other file replaces the current import.

   Returns false (and shows why) if it is not a DRF file; the current import is then unchanged.
   */
  bool addFile( const std::string &displayName, std::shared_ptr<const std::string> data );
  void addFile( std::shared_ptr<const DrfImport::ParsedFile> file );

  /** The imported DRF, or null if the import is not complete.  A new object after every change. */
  std::shared_ptr<DetectorPeakResponse> candidate() const;

  /** For a geometry-only file (a Detector.dat, an ANGLE .detx), the DRF to characterize with
   Monte Carlo ("Modify..."); null otherwise.
   */
  std::shared_ptr<DetectorPeakResponse> characterizationSeed() const;

  /** Emitted when #candidate or #characterizationSeed may have changed - only ever as a result
   of a file or a user edit.
   */
  Wt::Signal<> &changed();

  /** Emitted when the user clicks "Characterize..." for a geometry-only file. */
  Wt::Signal<> &characterizeRequested();

  /** Fills `dialog` as the dialog shown for a DRF file dropped onto the app: this widget, a chart,
   the default-DRF-for-serial-number/model options, and Cancel / Further options / Use buttons.
   Returns the widget, so further dropped files can be added to it.
   */
  static DrfImportWidget *setupDropDialog( SimpleDialog *dialog,
                                           std::shared_ptr<const DrfImport::ParsedFile> file );

protected:
  void handleUpload( const std::string &displayName, const std::string &spoolName,
                     const bool companionArea );
  void startSource( const size_t record );
  void setSource( std::shared_ptr<const DrfImport::Source> src );
  void handleRecordChanged();
  void updateCandidate();
  DrfImport::Options currentOptions() const;
  void updateFileList();
  void setStatus( const Wt::WString &txt, const bool isError );

  const Host m_host;

  std::unique_ptr<FileDragUploadResource> m_mainUpload;       ///< null for Host::DropDialog
  std::unique_ptr<FileDragUploadResource> m_companionUpload;

  Wt::WContainerWidget *m_mainDrop;                          ///< null for Host::DropDialog
  Wt::WContainerWidget *m_fileList;

  Wt::WContainerWidget *m_companionDrop;
  Wt::WText *m_companionTxt;

  Wt::WContainerWidget *m_recordDiv;
  Wt::WComboBox *m_recordCombo;

  Wt::WContainerWidget *m_interpDiv;
  Wt::WComboBox *m_interpCombo;
  std::vector<DrfImport::Interpretation> m_interpretations;  ///< parallel to m_interpCombo

  Wt::WContainerWidget *m_diameterDiv;
  Wt::WLineEdit *m_diameterEdit;
  Wt::WContainerWidget *m_setbackDiv;
  Wt::WLineEdit *m_setbackEdit;
  Wt::WContainerWidget *m_distanceDiv;
  Wt::WLineEdit *m_distanceEdit;

  Wt::WContainerWidget *m_eccUncertHolder;
  EccUncertOptions *m_eccUncert;

  Wt::WContainerWidget *m_nameDiv;
  Wt::WLineEdit *m_nameEdit;

  Wt::WText *m_status;
  Wt::WContainerWidget *m_notes;
  Wt::WPushButton *m_characterizeBtn;

  std::shared_ptr<const DrfImport::ParsedFile> m_primary;
  std::shared_ptr<const DrfImport::ParsedFile> m_companion;
  std::shared_ptr<const DrfImport::Source> m_source;

  /** Bumped on every file change, so a slow build that finishes after a newer file is dropped
   is discarded rather than shown.
   */
  uint64_t m_generation;

  std::shared_ptr<DetectorPeakResponse> m_candidate;
  std::shared_ptr<DetectorPeakResponse> m_seed;

  Wt::Signal<> m_changed;
  Wt::Signal<> m_characterizeRequested;
};//class DrfImportWidget

#endif //DrfImportWidget_h
