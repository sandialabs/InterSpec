#ifndef RefSpectraWidget_h
#define RefSpectraWidget_h

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

#include <Wt/WContainerWidget>

#include "InterSpec/SimpleDialog.h"

namespace Wt 
{
  class WTreeView;
  class WCheckBox;
  class WPushButton;
}

namespace SpecUtils{ enum class SpectrumType : int; }

class RefSpectraModel;
class RefSpectraWidget;


class RefSpectraDialog : public SimpleDialog
{
public:
  RefSpectraDialog( const Wt::WString &title = "Reference Spectra" );
  virtual ~RefSpectraDialog();

  RefSpectraWidget *widget() { return m_widget; }

  static RefSpectraDialog *createDialog( SpecUtils::SpectrumType type );
private:
  RefSpectraWidget *m_widget;
};//class RefSpectraDialog



class RefSpectraWidget : public Wt::WContainerWidget
{
public:
  RefSpectraWidget( Wt::WContainerWidget *parent = nullptr );
  virtual ~RefSpectraWidget();

  void setLoadSpectrumType( SpecUtils::SpectrumType type );


protected:
  void handleFileSelection( const Wt::WModelIndex &index );
  void addDirectory();
  void removeDirectory();
  void updatePreview();
  void updateTreeView();
  
private:
  RefSpectraModel *m_treeModel;
  Wt::WTreeView *m_treeView;
  Wt::WContainerWidget *m_previewContainer;
  Wt::WCheckBox *m_overlayCurrentDoc;
  Wt::WCheckBox *m_overlayOtherDoc;
  Wt::WContainerWidget *m_addDirButton;
  Wt::WContainerWidget *m_deleteDirButton;
  SpecUtils::SpectrumType m_loadSpectrumType;
  std::vector<std::string> m_baseDirectories;
  
  void setupUI();
  void createTreeModel();
};

#endif // RefSpectraWidget_h 