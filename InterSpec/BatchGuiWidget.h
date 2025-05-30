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

#include <Wt/WContainerWidget>

#include "InterSpec/SimpleDialog.h"

namespace Wt 
{
  class WComboBox;
  class WTreeView;
  class WCheckBox;
  class WPushButton;
  class WStackedWidget;
}

namespace SpecUtils{ enum class SpectrumType : int; }

class BatchGuiWidget;
class BatchGuiAnaWidget;
class DirectorySelector;
class FileDragUploadResource;
class BatchGuiInputSpectrumFile;

class BatchGuiDialog : public SimpleDialog
{
public:
  BatchGuiDialog( FileDragUploadResource *uploadResource, const Wt::WString &title = "Batch Analysis" );
  virtual ~BatchGuiDialog();

  BatchGuiWidget *widget() { return m_widget; }

  /** Create a new BatchGuiDialog
   */
  static BatchGuiDialog *createDialog( FileDragUploadResource *uploadResource );

protected:

private:
  BatchGuiWidget *m_widget;
  Wt::WPushButton *m_processBtn;
};//class BatchGuiDialog



class BatchGuiWidget : public Wt::WContainerWidget
{
public:
  BatchGuiWidget( FileDragUploadResource *uploadResource, Wt::WContainerWidget *parent = nullptr );
  virtual ~BatchGuiWidget();

  Wt::Signal<bool> &canDoAnalysis();

  void performAnalysis();

protected:
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
  BatchGuiAnaWidget *m_file_convert_opts;
  

  Wt::WContainerWidget *m_input_files_container;

  DirectorySelector *m_output_dir;

  bool m_can_do_analysis;

  /** Signal used to indicate if the button for doing the analysis should be enabled/disabled */
  Wt::Signal<bool> m_canDoAnalysis;
};//class BatchGuiWidget

#endif // BatchGuiWidget_h 
