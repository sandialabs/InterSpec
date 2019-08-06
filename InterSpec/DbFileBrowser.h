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

/**
 Refactored class holding all the UI components to show the snapshot/spectra (used in both load dialog and welcome dialog
 */
class SnapshotFactory
{
public:
    SnapshotFactory(SpecMeasManager *manager, InterSpec *viewer , std::string uuid,  SpectrumType type, std::shared_ptr<SpectraFileHeader> header, Wt::WGridLayout* layout,AuxWindow* window, Wt::WContainerWidget* footer);
    void loadSnapshotSelected();
    void loadSpectraSelected();
    int size() {return m_size; };
protected:
    void selectionChanged();
    std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserState > > m_UserStateLookup;
    std::map < Wt::WTreeNode *, Wt::Dbo::ptr < UserFileInDb > > m_UserFileInDbLookup;
    
    void addSpectraNodes(Wt::Dbo::collection< Wt::Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode);
    
    std::shared_ptr<DataBaseUtils::DbSession> m_session;
    SpecMeasManager  *m_manager;
    InterSpec   *m_viewer;
    Wt::WPushButton  *m_loadSnapshotButton;
    Wt::WPushButton  *m_loadSpectraButton;
    //  Wt::WPushButton  *m_deleteButton;
    Wt::WButtonGroup *m_buttonGroup;
    Wt::WGroupBox    *m_buttonbox;
    Wt::WTree        *m_snapshotTable;
    Wt::WText       *m_descriptionLabel;
    Wt::WText       *m_timeLabel;
    Wt::WGridLayout* m_relatedSpectraLayout;
    AuxWindow * m_window; //only set if in an auxwindow
    Wt::WContainerWidget* m_footer;
    std::string m_uuid;
    SpectrumType m_type;
    std::shared_ptr<SpectraFileHeader> m_header;
    std::shared_ptr<SpecMeas> retrieveMeas( const int dbid );
    int m_size; //size of snapshots populated
};


/**
 Main class to load snapshot/spectra.  Uses SnapshotFactory to fill in the UI
 */
class DbFileBrowser : public AuxWindow
{
public:
  DbFileBrowser( SpecMeasManager *manager, InterSpec *viewer , std::string uuid,  SpectrumType type, std::shared_ptr<SpectraFileHeader> header);
};

#endif //DbFileBrowser_h
