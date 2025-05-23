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
class FileDragUploadResource;

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

  void handleClose();
  void handleFileDrop( const std::string &displayName, const std::string &spooledName );
  Wt::Signal<bool> &canDoAnalysis();

  void addInputFiles( const std::vector<std::tuple<std::string,std::string,bool>> &files );

protected:
  void updateCanDoAnalysis();

  
  FileDragUploadResource *m_uploadResource;
  Wt::Signal<bool> m_canDoAnalysis;
  std::vector<std::tuple<std::string,std::string,bool>> m_inputFiles;
};//class BatchGuiWidget

#endif // BatchGuiWidget_h 
