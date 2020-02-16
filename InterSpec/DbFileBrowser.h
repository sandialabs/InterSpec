#ifndef DbFileBrowser_h
#define DbFileBrowser_h
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

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/SpectraFileModel.h"

namespace Wt
{
  class WPushButton;
  class WTree;
}//namespace Wt

namespace SpecUtils{ enum class SpectrumType : int; }

/**
 Widget to hold all the UI components to show the snapshot/spectra (used in both
 load dialog and welcome dialog).
 
 ToDo:
   - make this class inherit from WContainerWidget
   - add signals when there is an error, or a selection is made
   - Implement deleting a snapshot and renaming it
 
 */
class SnapshotBrowser : public Wt::WContainerWidget
{
public:
  SnapshotBrowser( SpecMeasManager *manager,
                   InterSpec *viewer,
                   const SpecUtils::SpectrumType type,
                   std::shared_ptr<SpectraFileHeader> header,
                   Wt::WContainerWidget *buttonBar = nullptr,
                   Wt::WContainerWidget *parent = nullptr );
  
   Wt::Signal<> &finished();
  
  void loadSnapshotSelected();
  void loadSpectraSelected();

  int numSnaphots() const;
  
protected:
  void selectionChanged();
  void startDeleteSelected();
  void startEditSelected();
  
  void deleteSelected();
  
  std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserState > > m_UserStateLookup;
  std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserFileInDb > > m_UserFileInDbLookup;
  
  void addSpectraNodes(Wt::Dbo::collection< Wt::Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode);
  
  std::shared_ptr<SpecMeas> retrieveMeas( const int dbid );
  
  std::shared_ptr<DataBaseUtils::DbSession> m_session;
  SpecMeasManager  *m_manager;
  InterSpec        *m_viewer;
  Wt::WPushButton  *m_loadSnapshotButton;
  Wt::WPushButton  *m_loadSpectraButton;
  //Wt::WPushButton  *m_deleteButton;
  Wt::WPushButton  *m_renameButton;
  Wt::WButtonGroup *m_buttonGroup;
  Wt::WGroupBox    *m_buttonbox;
  Wt::WTree        *m_snapshotTable;
  Wt::WText        *m_descriptionLabel;
  Wt::WText        *m_timeLabel;
  Wt::WGridLayout  *m_relatedSpectraLayout;
  SpecUtils::SpectrumType     m_type;
  std::shared_ptr<SpectraFileHeader> m_header;
  
  Wt::Signal<> m_finished;
  AuxWindow *m_editWindow;
  
  int m_nrows;
};//class SnapshotBrowser


/**
 Main class to load snapshot/spectra.  Uses SnapshotBrowser to fill in the UI
 */
class DbFileBrowser : public AuxWindow
{
public:
  DbFileBrowser( SpecMeasManager *manager,
                 InterSpec *viewer,
                 SpecUtils::SpectrumType type,
                 std::shared_ptr<SpectraFileHeader> header );
  
protected:
  SnapshotBrowser *m_factory;
};

#endif //DbFileBrowser_h
