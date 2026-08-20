#ifndef BatchGuiWidget_h
#define BatchGuiWidget_h

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

#include <tuple>
#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget.h>

#include "InterSpec/SimpleDialog.h"
#include "InterSpec/BatchSampleSelect.h"

namespace Wt
{
  class WText;
  class WComboBox;
  class WTreeView;
  class WCheckBox;
  class WPushButton;
  class WSelectionBox;
  class WStackedWidget;
}

class SpecMeas;
class SpectraFileHeader;

namespace SpecUtils{ enum class SpectrumType : int; }

class GroupBox;
class BatchGuiWidget;
class BatchGuiAnaWidget;
class DirectorySelector;
class FileDragUploadResource;
class BatchGuiInputSpectrumFile;


/* TODO items:
 - Add callback to batch act/shield and peak fit ana to display progress to user
 - 
 */


class BatchGuiDialog : public SimpleDialog
{
  friend class SimpleDialog;
public:
  virtual ~BatchGuiDialog();

  BatchGuiWidget *widget() { return m_widget; }

  /** Create a new BatchGuiDialog
   */
  static BatchGuiDialog *createDialog( FileDragUploadResource *uploadResource,
                                       const bool allow_adding_open_files );

protected:
  /** Constructor is protected; use SimpleDialog::make<BatchGuiDialog>() to create.

   @param allow_adding_open_files If true, a "load open file..." link is offered, so the user can
          pull in spectrum files they already have open in InterSpec.  This is wanted when the
          dialog was opened from within the app, but not when it was opened by dropping a bunch
          of files onto InterSpec.
   */
  BatchGuiDialog( FileDragUploadResource *uploadResource,
                  const bool allow_adding_open_files,
                  const Wt::WString &title = "Batch Analysis" );


private:
  BatchGuiWidget *m_widget;
  Wt::WPushButton *m_processBtn;
};//class BatchGuiDialog



class BatchGuiWidget : public Wt::WContainerWidget
{
public:
  BatchGuiWidget( FileDragUploadResource *uploadResource,
                  const bool allow_adding_open_files );
  virtual ~BatchGuiWidget();

  Wt::Signal<bool> &canDoAnalysis();

  void performAnalysis();

  /** Adds spectrum files that are already parsed and resident in memory - i.e., files the user has
   open in InterSpec - to the "Input Files" area.

   The `SpecMeas` objects are held, not copied; the analysis path already makes its own copy of
   every input before running, so the users file is never modified, and any edits they make before
   pressing "Analyze" are picked up.  Note that this does keep the `SpecMeas` in memory until the
   dialog is closed.

   Files whose `SpecMeas` is already in the input area are skipped, so calling this more than once
   with the same file is a no-op.

   The tuple is {display name, path to file on disk (may be empty), measurement}.
   */
  void addInMemoryFiles( const std::vector<std::tuple<std::string,std::string,std::shared_ptr<SpecMeas>>> &files );

  /** The measurements currently in the "Input Files" area. */
  std::vector<std::shared_ptr<const SpecMeas>> currentInputSpecMeas() const;

  /** Sets how input files holding more than one candidate foreground spectrum should be treated,
   on each of the analysis-type panes.
   */
  void setMultiSampleHandling( const BatchSampleSelect::MultiSampleHandling handling );

protected:
  /** Opens the dialog that lets the user pick from the spectrum files open in InterSpec. */
  void handleLoadOpenFileRequest();

  /** Adds the files the user selected in the `handleLoadOpenFileRequest()` dialog. */
  void handleAddOpenFiles( Wt::WSelectionBox *selection,
                           std::shared_ptr<std::vector<std::shared_ptr<SpectraFileHeader>>> headers );

  /** Shows the "load open file..." link only when there is an open file not already added. */
  void updateLoadOpenFileLinkVisibility();

  /** Number of files open in InterSpec that are not already in the input area. */
  size_t numAddableOpenFiles() const;

  void handleFileDrop( const std::string &disp_name, const std::string &spol_name );

  /** Takes ownership of the files, and will eventually delete them if the third argument is true.
   *  The tuple is {display name, path to file, should delete}
   *  Display name may either be the name provided to the http upload, or it may be the
   *  full path to the file for native apps.
   */
  void addInputFiles( const std::vector<std::tuple<std::string,std::string,bool>> &files );

  void handle_remove_input_file( BatchGuiInputSpectrumFile *input );

  void updateCanDoAnalysis();

protected:
  /** The resource used to upload the files.
   * This is _not_ owned by this class, but rather the SpecMeasManager.
  */
  FileDragUploadResource *m_uploadResource;
  
  Wt::WMenu *m_batch_type_menu;
  Wt::WStackedWidget *m_options_stack;
  BatchGuiAnaWidget *m_act_shield_ana_opts;
  BatchGuiAnaWidget *m_peak_fit_opts;
#if( USE_REL_ACT_TOOL )
  BatchGuiAnaWidget *m_iso_from_nucs_opts;
#endif
  BatchGuiAnaWidget *m_file_convert_opts;
  

  /** Wraps the input-files group box, so the "load open file..." link can be positioned against
   it.  The link is a sibling of the group box, rather than a child, both because the group box is
   a click-to-upload target, and because it scrolls.
   */
  Wt::WContainerWidget *m_input_files_holder;
  GroupBox *m_input_files_container;

  /** Link to add a spectrum file already open in InterSpec; null when not allowed. */
  Wt::WPushButton *m_load_open_file_btn;
  const bool m_allow_adding_open_files;

  DirectorySelector *m_output_dir;

  Wt::WText *m_input_status_error;

  bool m_can_do_analysis;

  /** Signal used to indicate if the button for doing the analysis should be enabled/disabled */
  Wt::Signal<bool> m_canDoAnalysis;

  const static int sm_max_spec_file_previews = 18;  //Arbitrarily chosen
};//class BatchGuiWidget

#endif // BatchGuiWidget_h 
