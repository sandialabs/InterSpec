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

#include <Wt/Dbo/ptr>
#include <Wt/WGroupBox>
#include <Wt/WComboBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WEnvironment>
#include <Wt/Dbo/QueryModel>
#include <Wt/WStringListModel>
#include <Wt/WLabel>
#include <Wt/WTabWidget>
#include <Wt/WMenuItem>
#include <Wt/WStandardItem>
#include <Wt/WHBoxLayout>
#include <Wt/WModelIndex>
#include <Wt/WTree>
#include <Wt/WTree>
#include <Wt/WTreeNode>
#include <Wt/WIconPair>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DbFileBrowser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/UseInfoWindow.h"

using namespace std;
using namespace Wt;


DbFileBrowser::DbFileBrowser( SpecMeasManager *manager, InterSpec *viewer, SpectrumType type, std::shared_ptr<SpectraFileHeader> header)
: AuxWindow( "Previously Stored States",
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                                             | AuxWindowProperties::DisableCollapse
                                             | AuxWindowProperties::EnableResize) ),
  m_factory( nullptr )
{
  WGridLayout *layout = stretcher();
  layout->setContentsMargins(0, 0, 0, 0);
  
  try
  {
    m_factory = new SnapshotBrowser( manager, viewer, type, header, footer(), nullptr );
    layout->addWidget( m_factory, 0, 0 );
    m_factory->finished().connect( this, &AuxWindow::hide );
  }catch( ... )
  {
    //In case something goes wrong inside SnapshotBrowser, we delete this window
    AuxWindow::deleteAuxWindow(this);
    return;
  }
  
  WPushButton *cancel = addCloseButtonToFooter();
  cancel->clicked().connect( this, &AuxWindow::hide );
  
  rejectWhenEscapePressed();
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );

  const int width = std::min( 500, static_cast<int>(0.95*viewer->renderedWidth()) );
  const int height = std::min( 475, static_cast<int>(0.95*viewer->renderedHeight()) );
  resize( WLength(width), WLength(height) );
  
  centerWindow();
  show();
}//DbFileBrowser

/*
SnapshotBrowser is the refactored class to create the UI for loading snapshot/spectra
*/
SnapshotBrowser::SnapshotBrowser( SpecMeasManager *manager,
                                  InterSpec *viewer,
                                  const SpectrumType type,
                                  std::shared_ptr<SpectraFileHeader> header,
                                  Wt::WContainerWidget *buttonBar,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_manager( manager ),
    m_viewer( viewer ),
    m_loadSnapshotButton( nullptr ),
    m_loadSpectraButton( nullptr ),
    m_deleteButton( nullptr ),
    m_renameButton( nullptr ),
    m_buttonGroup( nullptr ),
    m_buttonbox( nullptr ),
    m_snapshotTable( nullptr ),
    m_descriptionLabel( nullptr ),
    m_timeLabel( nullptr ),
    m_relatedSpectraLayout( nullptr ),
    m_type( type ),
    m_header( header ),
    m_size( 0 )
{
  WContainerWidget *footer = buttonBar;
  if( !footer )
    footer = new WContainerWidget();
  
  WGridLayout *layout = new WGridLayout();
  setLayout( layout );
  
  //We have to create a independant Dbo::Session for this class since the
  //  m_viewer->m_user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_session.reset( new DataBaseUtils::DbSession( *m_viewer->sql() ) );
  
  int row=0;
  if( m_header && !m_header->m_uuid.empty() )
  {
    layout->addWidget(new WText("It looks like you have previously loaded and modifed this spectrum file, would you like to resume your previous work?"), 0, 0);
    row++;
  } //!m_uuid.empty()
  
  try
  {
    Dbo::ptr<InterSpecUser> user = m_viewer->m_user;
    WContainerWidget* tablecontainer = new WContainerWidget();
    
    //With m_snapshotTable being a WTree, we cant implement double clicking on
    //  an item to open it right away.  If we use a WTreeView, then we should
    //  be able to do that...
    m_snapshotTable = new WTree();
    WGridLayout* tablelayout = new WGridLayout();
    tablelayout->setContentsMargins(2, 2, 2, 2);
    tablelayout->setRowStretch(0, 1);
    tablelayout->setColumnStretch(0, 1);
    tablelayout->addWidget(m_snapshotTable, 0, 0);
    tablecontainer->setLayout( tablelayout );
    tablecontainer->setMargin(0);
    tablecontainer->setOffsets(0);
    
    //        m_snapshotTable->resize(400,WLength(WLength::Auto));
    
    m_snapshotTable->itemSelectionChanged().connect(boost::bind( &SnapshotBrowser::selectionChanged, this ) );
    
    Wt::WTreeNode *root = new Wt::WTreeNode("Root");
    root->expand();
    
    m_snapshotTable->setTreeRoot(root);
    m_snapshotTable->setSelectionMode(Wt::SingleSelection);
    m_snapshotTable->treeRoot()->setNodeVisible( false ); //makes the tree look like a table! :)
    
    
    Dbo::Transaction transaction( *m_session->session() );
    Dbo::collection< Dbo::ptr<UserState> > query;
    
    if( m_header && !m_header->m_uuid.empty())
    { //uuid filter for snapshots with foreground uuid
      const char *stateQueryTxt = "A.InterSpecUser_id = ? AND A.StateType = ? AND A.SnapshotTagParent_id IS NULL";
      query = m_session->session()->query<Dbo::ptr<UserState> >(
                       "SELECT DISTINCT A FROM UserState A JOIN UserFileInDb B "
                         "ON (A.ForegroundId = B.id OR A.BackgroundId = B.id OR A.SecondForegroundId = B.id) "
                         "AND B.UUID = ? "
                         "LEFT JOIN UserState C "
                         "ON (C.ForegroundId = B.id OR C.BackgroundId = B.id OR C.SecondForegroundId = B.id) "
                         "AND C.SnapshotTagParent_id = A.id")
                      .bind(m_header->m_uuid)
                      .where(stateQueryTxt)
                      .bind( user.id() )
                      .bind(int(UserState::kUserSaved));
    }else
    {
      //No related hash/uuid filtering
      const char *stateQueryTxt = "InterSpecUser_id = ? AND StateType = ? AND SnapshotTagParent_id IS NULL";
      query = m_session->session()->find<UserState>()
      .where( stateQueryTxt )
      .bind( user.id() ).bind(int(UserState::kUserSaved));
    } //uuid filter for snapshots with foreground uuid
    
    m_size = query.size(); //save for future query
    
    if (m_size==0) {
      if( m_header && !m_header->m_uuid.empty() )
      {
        layout->addWidget(new WLabel("No saved application states"), row++, 0);
      }else
      {
        //This prevents the case where the checkIfPreviouslyOpened() is true,
        //  but the previous opened file is a kEndofSessionTemp, so will not
        //  show in the dialog.
        throw "Nothing to display";
      }//m_uuid.empty()
    } //m_size==0
    
    for( Dbo::collection< Dbo::ptr<UserState> >::const_iterator snapshotIterator = query.begin();
        snapshotIterator != query.end(); ++snapshotIterator )
    {
      Wt::WTreeNode *snapshotNode = new Wt::WTreeNode((*snapshotIterator)->name, 0, root);
      
      //m_snapshotTable->doubleClicked
      
      snapshotNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
      snapshotNode->setToolTip((*snapshotIterator)->serializeTime.toString() + (((*snapshotIterator)->description).empty()?"":(" -- " + (*snapshotIterator)->description)));
      
      //Add to lookup table (the HEAD)
      m_UserStateLookup[snapshotNode]=*snapshotIterator;
      
      //If not root, then return how many versions
      typedef Dbo::collection< Dbo::ptr<UserState> > Snapshots;
      typedef Snapshots::iterator SnapshotIter;
      Snapshots snapshots;
      
      if( *snapshotIterator )
      {
        //Create HEAD
        Wt::WTreeNode *versionNode = new Wt::WTreeNode("HEAD", new Wt::WIconPair("InterSpec_resources/images/time.png","InterSpec_resources/images/time.png"), snapshotNode);
        
        versionNode->setToolTip((*snapshotIterator)->serializeTime.toString() + (((*snapshotIterator)->description).empty()?"":(" -- " + (*snapshotIterator)->description)));
        
        versionNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
        
        //Add to lookup table (also the HEAD)
        m_UserStateLookup[versionNode]=*snapshotIterator;
        addSpectraNodes(snapshotIterator, versionNode);
        
        snapshots = m_session->session()->find<UserState>()
        .where( "SnapshotTagParent_id = ? AND StateType = ? " )
        .bind( (*snapshotIterator).id() ).bind(int(UserState::kUserSaved)).orderBy( "id desc" );
        
        for( Dbo::collection< Dbo::ptr<UserState> >::const_iterator versionIterator = snapshots.begin();
            versionIterator != snapshots.end(); ++versionIterator )
        {
          versionNode = new Wt::WTreeNode((*versionIterator)->name, new Wt::WIconPair("InterSpec_resources/images/time.png","InterSpec_resources/images/time.png"), snapshotNode);
          versionNode->setToolTip((*versionIterator)->serializeTime.toString()  + (((*versionIterator)->description).empty()?"":(" -- " + (*versionIterator)->description)));
          
          versionNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
          
          m_UserStateLookup[versionNode]=*versionIterator;
          
          addSpectraNodes(versionIterator, versionNode);
        } //for
      } //*snapshotIterator
      
      if( m_header && !m_header->m_uuid.empty())
        snapshotNode->expand();
    } //for
    
    transaction.commit();
    tablecontainer->setOverflow(Wt::WContainerWidget::OverflowAuto);
    layout->addWidget( tablecontainer, ++row, 0 );
    layout->setRowStretch( row, 1 );
    
    m_timeLabel = new WText();
    layout->addWidget(m_timeLabel, ++row,0);
    
    m_descriptionLabel = new WText();
    layout->addWidget(m_descriptionLabel, ++row,0);
    layout->columnStretch(1);
    
    m_buttonbox = new WGroupBox( "Open Related Spectrum As:" );
    m_buttonbox->setOffsets(10);
    m_buttonGroup = new WButtonGroup( m_buttonbox );
    
    for( int i = 0; i < 3; ++i )
    {
      const char *name = "";
      switch( SpectrumType(i) )
      {
        case kForeground:       name = "Foreground";      break;
        case kSecondForeground: name = "Second Spectrum"; break;
        case kBackground:       name = "Background";      break;
      }//switch( SpectrumType(i) )
      
      WRadioButton *button = new WRadioButton( name, m_buttonbox );
      button->setMargin( 10, Wt::Right );
      button->setInline(false);
      m_buttonGroup->addButton( button, i );
    }//for( int i = 0; i < 3; ++i )
    
    m_buttonGroup->setSelectedButtonIndex( 0 );
    m_buttonbox->setHidden(true);
    layout->addWidget( m_buttonbox, ++row, 0  );
    
    
    m_deleteButton = new WPushButton( "", footer );
    m_deleteButton->setHidden(true);
    m_deleteButton->setIcon( "InterSpec_resources/images/minus_min.png" );
    m_deleteButton->addStyleClass("FloatLeft");
    m_deleteButton->clicked().connect( this, &SnapshotBrowser::startDeleteSelected );
    
    //m_renameButton
    
    //m_renameButton
    m_loadSpectraButton = new WPushButton( "Load Spectrum Only", footer );
    m_loadSpectraButton->clicked().connect( boost::bind( &SnapshotBrowser::loadSpectraSelected, this));
    m_loadSpectraButton->disable();
    
    m_loadSnapshotButton = new WPushButton( "Load App. State", footer);
    m_loadSnapshotButton->clicked().connect( boost::bind(&SnapshotBrowser::loadSnapshotSelected, this));
    m_loadSnapshotButton->setDefault(true);
    //m_loadSnapshotButton->setIcon( "InterSpec_resources/images/time.png" );
    m_loadSnapshotButton->disable();
    
    if( !buttonBar )
      layout->addWidget( footer, layout->rowCount()+1 , 0, AlignRight );
  }catch( std::exception &e )
  {
    if( !buttonBar )
      delete footer;
    
    passMessage( (string("Error creating database spectrum file browser: ") + e.what()).c_str(),
                "", WarningWidget::WarningMsgHigh );
    cerr << "\n\nSnapshotModel caught: " << e.what() << endl << endl;
  }//try / catch
}//SnapshotBrowser


Wt::Signal<> &SnapshotBrowser::finished()
{
  return m_finished;
}


void SnapshotBrowser::addSpectraNodes(Dbo::collection< Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode)
{
    //add in foreground/background/2ndfore
    typedef Dbo::collection< Dbo::ptr<UserFileInDb> > Spectras;
    typedef Spectras::iterator SpectraIter;
   
    for( int i = 0; i < 3; i++ )
    { //loop through each spectrumtype because sometimes there can be one spectra that is in multiple types
        Spectras spectras =  m_session->session()->find< UserFileInDb >().where( " id = ? OR id = ? OR id = ?" ).bind((*versionIterator)->foregroundId).bind((*versionIterator)->backgroundId).bind((*versionIterator)->secondForegroundId);
        
        string pre = "";
        string post = "";
        SpectrumType spectratype = SpectrumType(i);

        for( Dbo::collection< Dbo::ptr<UserFileInDb> >::const_iterator spectraIterator = spectras.begin();
            spectraIterator != spectras.end(); ++spectraIterator )
        {
            Wt::WIconPair *icon = NULL;
            if ((*versionIterator)->foregroundId==spectraIterator->id() && spectratype==kForeground)
            {
                icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_forwards.png","InterSpec_resources/images/shape_move_forwards.png");
                
            } //foreground
            else if ((*versionIterator)->backgroundId==spectraIterator->id() && spectratype==kBackground)
            {
                icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_backwards.png","InterSpec_resources/images/shape_move_backwards.png");
            } //background
            else if ((*versionIterator)->secondForegroundId==spectraIterator->id() && spectratype==kSecondForeground)
            {
                icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_front.png","InterSpec_resources/images/shape_move_front.png");
            } //second foreground
            else
            {
                //nothing found, so just return
                continue;
            }
            if( m_header && !m_header->m_uuid.empty() && (*spectraIterator)->uuid==m_header->m_uuid)
            {
                //bold the spectra that is pulled in
                pre = "<b>";
                post = "</b>";
            } //!m_uuid.empty() && (*spectraIterator)->uuid==m_uuid
            else
            { //default look is not bolded
                pre = "";
                post = "";
            } //default look is not bolded
            
            Wt::WTreeNode *spectraNode = new Wt::WTreeNode(pre + " " + (*spectraIterator)->filename + " <i>(" + descriptionText(spectratype)+")</i>" + post, icon, versionNode);
            
            if( m_header && !m_header->m_uuid.empty() )
            {
                //disable all child spectra tree nodes, so the user can only select the snapshot
                spectraNode->disable();
            } //!m_uuid.empty()
            
            spectraNode->setToolTip((*spectraIterator)->serializeTime.toString()  + (((*spectraIterator)->description).empty()?"":(" -- " + (*spectraIterator)->description)));
            //Add to lookup table (also the HEAD)
            m_UserFileInDbLookup[spectraNode]=*spectraIterator;
            break;
        } //iterate through the spectra contained in this snapshot
    } //for( int i = 0; i < 3; i++ ) -- loop through spectrumtype
    
} //void SnapshotBrowser::addSpectraNodes(Dbo::collection< Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode)

//Updates the buttons when a row is selected
void SnapshotBrowser::selectionChanged()
{
    std::set<Wt::WTreeNode *> sets = m_snapshotTable->selectedNodes();

    //Only one row selected
    Wt::WTreeNode * selectedTreeNode = *sets.begin();
    
    if (m_UserStateLookup[selectedTreeNode])
    { //UserState clicked
        m_buttonbox->disable();
        m_buttonbox->hide();
        m_deleteButton->show();
        m_loadSnapshotButton->enable();
        m_loadSpectraButton->disable();
        m_descriptionLabel->setText("<i>"+m_UserStateLookup[selectedTreeNode]->description+"</i>");
        m_timeLabel->setText("<b>"+m_UserStateLookup[selectedTreeNode]->serializeTime.toString()+"</b>");
    } //UserState clicked
    else if (m_UserFileInDbLookup[selectedTreeNode])
    { //UserFileDB clicked
        m_buttonbox->show();
        m_buttonbox->enable();
//        m_deleteButton->enable();
        m_loadSpectraButton->enable();
      m_deleteButton->hide();
        m_loadSnapshotButton->disable();
        m_descriptionLabel->setText("<i>"+m_UserFileInDbLookup[selectedTreeNode]->description+"</i>");
        m_timeLabel->setText("<b>"+m_UserFileInDbLookup[selectedTreeNode]->serializeTime.toString()+"</b>");
        
    } //UserFileDB clicked
    else
    { //some node not found
        m_buttonbox->hide();
        m_buttonbox->disable();
//        m_deleteButton->disable();
        m_loadSnapshotButton->disable();
      m_deleteButton->hide();
        m_loadSpectraButton->disable();
        m_descriptionLabel->setText("");
        m_timeLabel->setText("");
    } //some node not found
    
    if( (m_header && !m_header->m_uuid.empty()) ||  !m_viewer->measurment( kForeground ) )
    {
      m_buttonbox->hide(); //make sure it is hidden if no uuid selected
    }//!m_uuid.empty()
} //void SnapshotBrowser::selectionChanged()


void SnapshotBrowser::startDeleteSelected()
{
  //deleteSelected();
}//void startDeleteSelected()

std::shared_ptr<SpecMeas> SnapshotBrowser::retrieveMeas( const int dbid )
  {
    std::shared_ptr<SpecMeas> snapshotmeas;
    if( dbid < 0 )
      return snapshotmeas;
    
    try
    {
      Dbo::Transaction transaction( *m_session->session() );
      Dbo::ptr<UserFileInDb> dbsnapshot = m_session->session()->find<UserFileInDb>()
      .where( "id = ?" )
      .bind( dbid );
      
      Dbo::ptr<UserFileInDbData> data;
      if( dbsnapshot && dbsnapshot->filedata.size() )
        data = *(dbsnapshot->filedata.begin());
      transaction.commit();
      if( data )
        snapshotmeas = data->decodeSpectrum();
    }catch( std::exception &e )
    {
      cerr << "retrieveMeas() caught (while trying to load snaphot): "
      << e.what() << endl;
      passMessage( "Sorry, couldnt load requested snapshot",
                  "", WarningWidget::WarningMsgHigh );
    }//try / catch
    
    return snapshotmeas;
  }//std::shared_ptr<SpecMeas> retrieveMeas( const int dbid )
  
  //
  // Delete a saved spectrum
  //
//  void SnapshotBrowser::deleteSelected()
//  {

//currently disabled deleting
      
//    WModelIndexSet indices = m_spectraTable->selectedIndexes();
//    if( !indices.size() )
//    {
//      m_deleteButton->disable();
//      return;
//    }//if( !indices.size() )
//    
//    WModelIndex index = *indices.begin();
//    
//    Dbo::ptr<UserFileInDb> dbfile = m_model->stableResultRow( index.row() );
//    
//    //Now we have to make 'dbfile' be associated with same session of
//    //  m_viewer->m_user.session()
//    if( dbfile.id() >= 0 )
//    {
//      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
//      DataBaseUtils::DbTransaction transaction( *sql );
//      dbfile = sql->session()->find< UserFileInDb >()
//                               .where( "id = ?").bind( dbfile.id() );
//
//      //This will now remove the dbfile from the assigned and loaded states (if it the same)
//      m_manager->removeSpecMeas(dbfile);
//        
//      dbfile.remove();
//      transaction.commit();
//      m_model->reload();
//
//    }//if( dbfile.id() >= 0 )
//  } // deleteSelected()

void SnapshotBrowser::loadSnapshotSelected()
{
      std::set<Wt::WTreeNode *> sets = m_snapshotTable->selectedNodes();
      
      //Only one row selected
      Wt::WTreeNode * selectedTreeNode = *sets.begin();
      
        if (!m_UserStateLookup[selectedTreeNode])
        {
            m_loadSnapshotButton->disable();
          m_deleteButton->hide();
            return;
        }//if( !indices.size() )
        Dbo::ptr<UserState> dbstate = m_UserStateLookup[selectedTreeNode];
        
        //Now we have to make 'dbfile' be associated with same session of
        //  m_viewer->m_user.session()  {I dont know what would happen otherwise)
        if( dbstate.id() >= 0 )
        {
            std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
            DataBaseUtils::DbTransaction transaction( *sql );
            dbstate = sql->session()->find< UserState >()
            .where( "id = ?").bind( dbstate.id() );
            transaction.commit();
        }//if( dbfile.id() >= 0 )
        
        if( !dbstate || dbstate.id() < 0 )
        {
            passMessage( "Error loading from the database",
                        "", WarningWidget::WarningMsgHigh );
            m_finished.emit();
            return;
        }//if( !selected )
        
        try
        {
            m_viewer->loadStateFromDb( dbstate );
            
            WString msg = "Snapshot '";
            msg += dbstate.get()->name;
            msg += "' loaded";
            passMessage( msg.toUTF8(), "", WarningWidget::WarningMsgInfo );
            
        }catch( std::exception &e )
        {
            passMessage( "Failed to load state", "", WarningWidget::WarningMsgHigh );
            cerr << "DbStateBrowser::loadSelected() caught: " << e.what() << endl;
        }//try / catch
  
  m_finished.emit();
}//void loadSnapshotSelected()


void SnapshotBrowser::loadSpectraSelected()
{
    std::set<Wt::WTreeNode *> sets = m_snapshotTable->selectedNodes();
    
    //Only one row selected
    Wt::WTreeNode * selectedTreeNode = *sets.begin();
   
        
        if (!m_UserFileInDbLookup[selectedTreeNode])
        {
            m_loadSpectraButton->disable();
            return;
        }//if( !indices.size() )
        Dbo::ptr<UserFileInDb> dbfile = m_UserFileInDbLookup[selectedTreeNode];

    //Now we have to make 'dbfile' be associated with same session of
    //  m_viewer->m_user.session()
    if( dbfile.id() >= 0 )
    {
        std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
        DataBaseUtils::DbTransaction transaction( *sql );
        dbfile = sql->session()->find< UserFileInDb >()
        .where( "id = ?").bind( dbfile.id() );
        transaction.commit();
    }//if( dbfile.id() >= 0 )
    
    if( !dbfile || dbfile.id() < 0 )
    {
        passMessage( "Error loading from the database",
                    "", WarningWidget::WarningMsgHigh );
        m_finished.emit();
        return;
    }//if( !selected )

    if( !m_header || m_header->m_uuid.empty() )
    {
        
        int snapshot_id = -1;
        
        SpectrumType type = kForeground;
        if( m_buttonbox->isVisible() )
            type = SpectrumType( m_buttonGroup->checkedId() );
        
        int modelrow = -1;
        std::shared_ptr<SpecMeas> measurement;
        std::shared_ptr<SpectraFileHeader> header;
        
        //Now should check if already opened
        SpectraFileModel *specmodel = m_manager->model();
        for( int row = 0; row < specmodel->rowCount(); ++row )
        {
            std::shared_ptr<SpectraFileHeader> thisheader
            = specmodel->fileHeader( row );
            Wt::Dbo::ptr<UserFileInDb> dbentry = thisheader->dbEntry();
            if( dbentry.id() == dbfile.id() )
            {
                modelrow = row;
                dbfile = dbentry;
                header = thisheader;
                measurement = thisheader->parseFile();
                if( measurement == m_viewer->measurment(type) )
                {
                    if( snapshot_id >= 0 )
                    {
                        std::shared_ptr<SpecMeas> meas = retrieveMeas( snapshot_id );
                        if( meas )
                        {
                            thisheader->setMeasurmentInfo( meas );
                            m_manager->displayFile( row, meas, type, false, false, false );
                        }//if( meas )
                    }//if( snapshot_id )
                    
                    passMessage("Already opened","",WarningWidget::WarningMsgInfo);
                    m_finished.emit();
                    return;
                } //measurement == m_viewer->measurment(type)
                break;
            } //if( dbentry.id() == dbfile.id() )
        }//for( int row = 0; row < specmodel->rowCount(); ++row )
        
        if( modelrow < 0 )
            modelrow = m_manager->setDbEntry( dbfile, header, measurement, true );
        
        m_manager->displayFile( modelrow, measurement, type, false, false, false );
    }
    else
    {
        //Copied from 
        //dorevert() is only called from within the application loop
        SpectraFileModel *m_model = m_manager->model();
        if( !dbfile || !m_header || !m_model || !wApp )
            throw runtime_error( "PreviousDbEntry: invalid input or no wApp" );
        
        if( dbfile->userHasModified )
        {
            m_header->setNotACandiateForSavingToDb();
            Wt::WModelIndex index = m_model->index( m_header );
            m_model->removeRows( index.row(), 1 );
            
            std::shared_ptr< SpectraFileHeader > header;
            std::shared_ptr< SpecMeas >  measurement;
            
            //go through and make sure file isnt already open
            for( int row = 0; row < m_model->rowCount(); ++row )
            {
                std::shared_ptr<SpectraFileHeader> header
                = m_model->fileHeader( row );
                Wt::Dbo::ptr<UserFileInDb> entry = header->dbEntry();
                if( entry && entry.id() == dbfile.id() )
                {
                    measurement = header->parseFile();
                    m_manager->displayFile( row, measurement, kForeground, false, false, false );
                    m_finished.emit();
                    return;
                }//if( entry.id() == dbfile.id() )
            }//for( int row = 0; row < m_model->rowCount(); ++row )
            
            try
            {
                const int modelRow = m_manager->setDbEntry( dbfile, header,
                                                           measurement, true );
                m_manager->displayFile( modelRow, measurement, kForeground, false, false, false );
                
            }catch( exception &e )
            {
                cerr << "\n\nSnapshotModel\n\tCaught: " << e.what() << "\n\n";
                passMessage( "Error displaying previous measurment, things may not"
                            " be as expected" , "", WarningWidget::WarningMsgHigh );
            }//try / catch
        }else
        {
            m_header->setDbEntry( dbfile );
        }//if( dbfile->userHasModified )
    }
  
  m_finished.emit();
}//void loadSpectraSelected

