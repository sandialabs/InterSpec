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

#include <Wt/WTree>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/Dbo/ptr>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WGroupBox>
#include <Wt/WTextArea>
#include <Wt/WTreeNode>
#include <Wt/WIconPair>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>


#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DbFileBrowser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"


using namespace std;
using namespace Wt;

namespace
{
  // The time and snapshot descriptions, if totally blank, take up no space, so if you select the
  //  bottom most snapshot, then suddenly the time/description will take up space, shrinking the
  //  table, and causing newly selected state to not be visible.  A way to fix this is to always
  //  have the time/description take up space via a "&nbsp;", but I'm not sure I like all this
  //  wasted space, so we'll leave this kinda wierdness alone for now
  const char * const ns_empty_descrip_label = ""; //"&nbsp;"
}//namespace


namespace
{
/** A class representing a spectrum file in the database; it displays very basic information
 about the spectrum file, and allows the user the ability to load it into InterSpec
 */
class DbSpecFileItem : public WContainerWidget
{
public:
  DbSpecFileItem( SpecUtils::SpectrumType type,
                  SpectraFileModel *model,
                  SpecMeasManager *manager,
                  Dbo::ptr<UserFileInDb> dbentry,
                  std::shared_ptr<SpectraFileHeader> header,
                  WContainerWidget *parent = nullptr )
  : WContainerWidget( parent ),
  m_type( type ),
  m_model( model ),
  m_manager( manager ),
  m_dbentry( dbentry ),
  m_header( header ),
  m_loadedASelection( this )
  {
    if( !m_dbentry || !m_model )
      throw runtime_error( "DbSpecFileItem: invalid input or no wApp" );
    
    addStyleClass( "DbSpecFileItem" );
    
    string msg = "Uploaded: " + dbentry->uploadTime.toString( DATE_TIME_FORMAT_STR ).toUTF8();
    if( dbentry->userHasModified )
      msg += ", was modified";
    
    WText *txt = new WText( msg, this );
    txt->addStyleClass( "DbSpecFileItemTxt" );
    
    WPushButton *button = new WPushButton( "Resume From", this );
    button->addStyleClass( "DbSpecFileItemButton" );
    button->clicked().connect( this, &DbSpecFileItem::dorevert );
    button->setFocus();
  }//DbSpecFileItem(...)
  
  
  void dorevert()
  {
    assert( wApp );
    
    if( !m_header || m_dbentry->userHasModified )
    {
      if( m_header )
      {
        m_header->setNotACandiateForSavingToDb();
        Wt::WModelIndex index = m_model->index( m_header );
        m_model->removeRows( index.row(), 1 );
      }//if( m_header )
      
      std::shared_ptr< SpectraFileHeader > header;
      std::shared_ptr< SpecMeas >  measurement;
      
      //go through and make sure file isnt already open
      for( int row = 0; row < m_model->rowCount(); ++row )
      {
        std::shared_ptr<SpectraFileHeader> header = m_model->fileHeader( row );
        Wt::Dbo::ptr<UserFileInDb> entry = header->dbEntry();
        
        if( entry && entry.id() == m_dbentry.id() )
        {
          measurement = header->parseFile();
          m_manager->displayFile( row, measurement, m_type, false, false, SpecMeasManager::VariantChecksToDo::None );
          m_loadedASelection.emit();
          
          return;
        }//if( entry.id() == m_dbentry.id() )
      }//for( int row = 0; row < m_model->rowCount(); ++row )
      
      try
      {
        const int modelRow = m_manager->setDbEntry( m_dbentry, header, measurement, true );
        m_manager->displayFile( modelRow, measurement, m_type, false, false, SpecMeasManager::VariantChecksToDo::None );
        m_loadedASelection.emit();
      }catch( exception &e )
      {
        cerr << "\n\nDbSpecFileItem::dorevert()\n\tCaught: " << e.what() << "\n\n";
        passMessage( "Error displaying previous measurement, things may not be as expected" ,
                    WarningWidget::WarningMsgHigh );
      }//try / catch
    }else
    {
      m_header->setDbEntry( m_dbentry );
    }//if( m_dbentry->userHasModified )
  }//void dorevert()
  
  ~DbSpecFileItem() noexcept(true)
  {
    try
    {
      m_dbentry.reset();
    }catch(...)
    {
      cerr << "DbSpecFileItem destructor caught exception doing m_dbentry.reset()" << endl;
    }
    
    try
    {
      m_header.reset();
    }catch(...)
    {
      cerr << "DbSpecFileItem destructo caught exception doing m_header.reset()" << endl;
    }
  }//~DbSpecFileItem()
  
  SpecUtils::SpectrumType m_type;
  SpectraFileModel *m_model;
  SpecMeasManager *m_manager;
  Dbo::ptr<UserFileInDb> m_dbentry;
  std::shared_ptr<SpectraFileHeader> m_header;
  Wt::Signal<> m_loadedASelection;
};//class DbSpecFileItem

}//namespace



DbFileBrowser::DbFileBrowser( SpecMeasManager *manager,
                              InterSpec *viewer,
                              std::shared_ptr<SpectraFileHeader> header)
: AuxWindow( "Previously Stored States",
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse)
                                             | AuxWindowProperties::EnableResize) ),
  m_factory( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/DbFileBrowser.css" );
  
  WGridLayout *layout = stretcher();
  layout->setContentsMargins(0, 0, 0, 0);
  
  try
  {
    m_factory = new SnapshotBrowser( manager, viewer, header, footer(), nullptr );
    layout->addWidget( m_factory, 0, 0 );
    
    m_factory->finished().connect( this, &AuxWindow::deleteSelf );
  }catch( std::exception &e )
  {
    m_factory = nullptr;
    WText *tt = new WText( "<b>Error creating SnapshotBrowser</b> - sorry :("
                          "<br />Error: " + string(e.what()) );
    layout->addWidget( tt, 0, 0 );
  }catch( ... )
  {
    WText *tt = new WText( "<b>Unexpected issue creating SnapshotBrowser</b> - sorry :(" );
    layout->addWidget( tt, 0, 0 );
  }
  
  WPushButton *cancel = addCloseButtonToFooter();
  cancel->clicked().connect( this, &AuxWindow::hide );
  
  rejectWhenEscapePressed();
  
  finished().connect( this, &AuxWindow::deleteSelf );

  const int width = std::min( 500, static_cast<int>(0.95*viewer->renderedWidth()) );
  const int height = std::min( 475, static_cast<int>(0.95*viewer->renderedHeight()) );
  resize( WLength(width), WLength(height) );
  
  centerWindow();
  show();
}//DbFileBrowser


int DbFileBrowser::numSnapshots() const
{
  return m_factory ? m_factory->numSnaphots() : 0;
}


Dbo::collection< Dbo::ptr<UserState> >
SnapshotBrowser::get_user_states_collection( Dbo::ptr<InterSpecUser> user,
                                             std::shared_ptr<DataBaseUtils::DbSession> &session,
                                             std::shared_ptr<const SpectraFileHeader> header )
{
  Dbo::collection< Dbo::ptr<UserState> > query;
  
  if( !session || !user )
    return query;
  
  if( header && !header->m_uuid.empty())
  { //uuid filter for snapshots with foreground uuid
    const char *stateQueryTxt = "A.InterSpecUser_id = ? AND A.StateType = ? AND A.SnapshotTagParent_id IS NULL";
    query = session->session()->query<Dbo::ptr<UserState> >(
              "SELECT DISTINCT A FROM UserState A JOIN UserFileInDb B "
              "ON (A.ForegroundId = B.id OR A.BackgroundId = B.id OR A.SecondForegroundId = B.id) "
              "AND B.UUID = ? "
              "LEFT JOIN UserState C "
              "ON (C.ForegroundId = B.id OR C.BackgroundId = B.id OR C.SecondForegroundId = B.id) "
              "AND C.SnapshotTagParent_id = A.id")
            .bind(header->m_uuid)
            .where(stateQueryTxt)
            .bind( user.id() )
            .bind(int(UserState::kUserSaved));
  }else
  {
    //No related hash/uuid filtering
    const char *stateQueryTxt = "InterSpecUser_id = ? AND StateType = ? AND SnapshotTagParent_id IS NULL";
    query = session->session()->find<UserState>()
            .where( stateQueryTxt )
            .bind( user.id() ).bind(int(UserState::kUserSaved));
  }//uuid filter for snapshots with foreground uuid
  
  return query;
}//get_user_states_collection(...)


size_t SnapshotBrowser::num_saved_states( InterSpec *viewer,
                                          shared_ptr<DataBaseUtils::DbSession> session,
                                          shared_ptr<const SpectraFileHeader> header )
{
  if( !session || !viewer )
    return 0;
  
  size_t num_states = 0;
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *session );
    
    Dbo::collection< Dbo::ptr<UserState> > query = SnapshotBrowser::get_user_states_collection( viewer->m_user, session, header );
    
    num_states = query.size();
    
    transaction.commit();
  }catch( std::exception &e )
  {
    cerr << "\n\n\nSnapshotBrowser::num_saved_states: caught exception I wasnt expecting to:\n\t" << e.what() << "\n\n\n" << endl;
  }//
  
  return num_states;
}//size_t num_saved_states(...)


/*
SnapshotBrowser is the refactored class to create the UI for loading snapshot/spectra
*/
SnapshotBrowser::SnapshotBrowser( SpecMeasManager *manager,
                                  InterSpec *viewer,
                                  std::shared_ptr<SpectraFileHeader> header,
                                  Wt::WContainerWidget *buttonBar,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_manager( manager ),
    m_viewer( viewer ),
    m_loadSnapshotButton( nullptr ),
    m_loadSpectraButton( nullptr ),
    //m_deleteButton( nullptr ),
    m_renameButton( nullptr ),
    m_buttonGroup( nullptr ),
    m_buttonbox( nullptr ),
    m_snapshotTable( nullptr ),
    m_descriptionLabel( nullptr ),
    m_timeLabel( nullptr ),
    m_relatedSpectraLayout( nullptr ),
    m_header( header ),
    m_finished( this ),
    m_editWindow( nullptr ),
    m_nrows( 0 )
{
  wApp->useStyleSheet( "InterSpec_resources/DbFileBrowser.css" );
  
  addStyleClass( "SnapshotBrowser" );
  
  WContainerWidget *footer = buttonBar;
  if( !footer )
    footer = new WContainerWidget();
  
  footer->addStyleClass( "SnapshotBrowserFooter" );
  
  WGridLayout *layout = new WGridLayout();
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  setLayout( layout );
  
  //We have to create a independent Dbo::Session for this class since the
  //  m_viewer->m_user.session() could be used in other threads, messing
  //  things up (Dbo::Session is not thread safe).
  m_session.reset( new DataBaseUtils::DbSession( *m_viewer->sql() ) );
  
  int row = 0;
  if( m_header && !m_header->m_uuid.empty() )
  {
    auto txt = new WText("It looks like you have saved this spectrum file as part of a save-state;"
                         " would you like to resume your previous work?");
    txt->addStyleClass( "ResumeWorkTxt" );
    
    layout->addWidget( txt, row++, 0);
  }//if( we want to load a specific state )
  
  try
  {
    Dbo::ptr<InterSpecUser> user = m_viewer->m_user;
    WContainerWidget* tablecontainer = new WContainerWidget();
    
    m_snapshotTable = new WTree();
    m_snapshotTable->addStyleClass( "SnapshotTable" );

    
    WGridLayout *tablelayout = new WGridLayout();
    tablelayout->setContentsMargins( 0, 0, 0, 0 );
    tablelayout->setVerticalSpacing( 0 );
    tablelayout->addWidget(m_snapshotTable, 0, 0);
    tablelayout->setRowStretch(0, 1);
    tablelayout->setColumnStretch(0, 1);
    tablecontainer->setLayout( tablelayout );
    
    m_snapshotTable->itemSelectionChanged().connect( this, &SnapshotBrowser::selectionChanged );
    
    Wt::WTreeNode *root = new Wt::WTreeNode( "Root" );
    root->expand();
    
    m_snapshotTable->setTreeRoot(root);
    m_snapshotTable->setSelectionMode(Wt::SingleSelection);
    m_snapshotTable->treeRoot()->setNodeVisible( false ); //makes the tree look like a table! :)
    
    DataBaseUtils::DbTransaction transaction( *m_session );
    
    Dbo::collection< Dbo::ptr<UserState> > query = SnapshotBrowser::get_user_states_collection( m_viewer->m_user, m_session, m_header );
    
    m_nrows = query.size(); //save for future query
    
    if( m_nrows == 0 )
      layout->addWidget( new WLabel("No saved application states"), row++, 0 );
    
    
    //Some sudo-styling for if we want to style this - or if we implement this as a WTreeView
    //  most of this styling comes along much easier (so I think this should be what is done
    //  sometime in the future).
    //if( !wApp->styleSheet().isDefined("snapshotrow") )
    //{
    //  // We would add these rules to InterSpec.css
    //  wApp->styleSheet().addRule( "li.SnapshotRow > .Wt-item", "height: 30px; line-height: 30px;", "snapshotrow" );
    //  wApp->styleSheet().addRule( ".SnapshotRow", "text-overflow: ellipsis; white-space: nowrap; overflow: hidden;" );
    //  wApp->styleSheet().addRule( ".SnapshotRow .Wt-expand", "height: 30px; padding-top: 5px" );
    //  wApp->styleSheet().addRule( ".SnapshotRow div.Wt-item.Wt-trunk,"
    //                              ".Wt-tree.Wt-trunk,"
    //                              ".Wt-tree .Wt-item.Wt-trunk",
    //                              "background-image: none !important;" );
    //  //Make table alternate row colors using:
    //  wApp->styleSheet().addRule( "ul.Wt-root > li.SnapshotRow", "color: white; background: rgb(41,42,44) !important;" );
    //  wApp->styleSheet().addRule( "ul.Wt-root > li.SnapshotRow:nth-child(odd)", "background: rgb(30,32,34) !important;" );
    //}
    
    for( Dbo::collection< Dbo::ptr<UserState> >::const_iterator snapshot_iter = query.begin();
        snapshot_iter != query.end(); ++snapshot_iter )
    {
      const Dbo::ptr<UserState> &snapshot = *snapshot_iter;
      Wt::WTreeNode *snapshotNode = new Wt::WTreeNode( snapshot->name, 0, root );
    
      //snapshotNode->addStyleClass( "SnapshotRow" );
      
      // We will go-around Wt::WTreeNode API to hack in a way to customize display and options
      //  each row in the table can have by assuming the row is implemented as a WContainerWidget,
      //  which seems to always be true, at least for Wt 3.3.4.
      // \TODO: sub-class WTreeNode in order to allow customization, which is probably a cleaner
      //       solution, or we could use a WTreeView instead of WTree, which is probably best
      //       solution.
      WContainerWidget *rowDiv = nullptr;
      if( snapshotNode->label() )
        rowDiv = dynamic_cast<WContainerWidget *>( snapshotNode->label()->parent() );
      
      if( !rowDiv && snapshotNode->label() )
      {
        //In case Wt changes in new versions
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "nsnapshotNode->label()->parent() Is Not a WContainerWidget" );
#endif
        cerr << "\n\nsnapshotNode->label()->parent() Is Not a WContainerWidget\n\n" << endl;
      }//if( Wt has changed with new version )
      
      // Add ability to double click on a row to load this state
      if( rowDiv )
      {
        auto loader = wApp->bind( boost::bind( &SnapshotBrowser::loadSnapshotSelected, this ) );
        rowDiv->doubleClicked().connect( std::bind([this,loader,snapshotNode](){
          const set<WTreeNode *> sets = m_snapshotTable->selectedNodes();
          if( !sets.empty() && ((*begin(sets)) == snapshotNode) )
            loader();
#if( PERFORM_DEVELOPER_CHECKS )
          else
            log_developer_error( __func__, "Got double click on SnapshotBrowser state but that wasnt what was selected." );
#endif
        }) );
      }//if( rowDiv )
      
      
      // If a 'buttonBar' was specified to this contructor, then this browser is NOT on the
      //  InterSpec introductory screen, but instantiated by the "InterSpec" --> "Previous..." menu.
      if( buttonBar && rowDiv )
      {
        WImage *delBtn = new WImage( "InterSpec_resources/images/minus_min_black.svg", rowDiv );
        delBtn->resize( WLength(16), WLength(16) );
        delBtn->addStyleClass( "DelSnapshotBtn" );
        delBtn->clicked().connect( this, &SnapshotBrowser::startDeleteSelected );
        delBtn->doubleClicked().preventPropagation();
        
        
        WImage *editBtn = new WImage( "InterSpec_resources/images/edit_black.svg", rowDiv );
        editBtn->resize( WLength(16), WLength(16) );
        editBtn->addStyleClass( "DelSnapshotBtn" );
        editBtn->clicked().connect( this, &SnapshotBrowser::startEditSelected );
        editBtn->doubleClicked().preventPropagation();
      }//if( buttonBar && rowDiv )
      
      //snapshotNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
      snapshotNode->setChildCountPolicy(Wt::WTreeNode::Disabled);
      WString tooltip = snapshot->serializeTime.toString()
                        + (snapshot->description.empty() ? "" : " -- ")
                        + snapshot->description;
      snapshotNode->setToolTip( tooltip );
      
      //Add to lookup table (the HEAD)
      m_UserStateLookup[snapshotNode] = snapshot;
      
      //If not root, then return how many versions
      typedef Dbo::collection< Dbo::ptr<UserState> > Snapshots;
      Snapshots snapshots;
      
      //Create HEAD
      auto headIcons = new Wt::WIconPair("InterSpec_resources/images/time.svg","InterSpec_resources/images/time.svg");
      headIcons->icon1()->addStyleClass( "Wt-icon SnapshotIcon" );
      headIcons->icon2()->addStyleClass( "Wt-icon SnapshotIcon" );
      
      Wt::WTreeNode *head_node = new Wt::WTreeNode("HEAD", headIcons, snapshotNode);
      
      tooltip = snapshot->serializeTime.toString()
                + (snapshot->description.empty() ? "" : " -- ")
                + snapshot->description;
      head_node->setToolTip( tooltip );
      
      head_node->setChildCountPolicy(Wt::WTreeNode::Disabled);
      
      //Add to lookup table (also the HEAD)
      m_UserStateLookup[head_node] = snapshot;
      addSpectraNodes( snapshot_iter, head_node );
      
      if( m_header && !m_header->m_uuid.empty())
        snapshotNode->expand();
      
      snapshots = m_session->session()->find<UserState>()
        .where( "SnapshotTagParent_id = ? AND StateType = ? " )
        .bind( snapshot.id() )
        .bind(int(UserState::kUserSaved))
        .orderBy( "id desc" );
      
      for( Dbo::collection<Dbo::ptr<UserState>>::const_iterator version_iter = snapshots.begin();
          version_iter != snapshots.end(); ++version_iter )
      {
        const Dbo::ptr<UserState> &version = *version_iter;
        
        const Wt::WString &name = version->name;
        auto tagIcon = new Wt::WIconPair("InterSpec_resources/images/time.svg","InterSpec_resources/images/time.svg");
        tagIcon->icon1()->addStyleClass( "Wt-icon SnapshotIcon" );
        tagIcon->icon2()->addStyleClass( "Wt-icon SnapshotIcon" );
        
        Wt::WTreeNode *versionNode = new Wt::WTreeNode( name, tagIcon, snapshotNode );
        
        tooltip = version->serializeTime.toString()
                  + (version->description.empty() ? "" : " -- ")
                  + version->description;
        versionNode->setToolTip( tooltip );
        
        versionNode->setChildCountPolicy( Wt::WTreeNode::Disabled );
        
        m_UserStateLookup[versionNode] = *version_iter;
        
        addSpectraNodes( version_iter, versionNode );
      }//for( loop over snapshots )
      
      
      if( m_header && !m_header->m_uuid.empty() )
      {
        snapshotNode->expand();
        head_node->expand();
      }
    }//for( loop over query )
    
    transaction.commit();
    layout->addWidget( tablecontainer, ++row, 0 );
    layout->setRowStretch( row, 1 );
    
    m_timeLabel = new WText( ns_empty_descrip_label );
    m_timeLabel->addStyleClass( "SnapshotTime" );
    layout->addWidget(m_timeLabel, ++row,0);
    
    m_descriptionLabel = new WText( ns_empty_descrip_label );
    m_descriptionLabel->addStyleClass( "SnapshotDesc" );
    layout->addWidget(m_descriptionLabel, ++row,0);
    
    m_buttonbox = new WGroupBox( "Open selected spectrum as:" );
    m_buttonbox->setOffsets(10);
    m_buttonGroup = new WButtonGroup( m_buttonbox );
    
    for( int i = 0; i < 3; ++i )
    {
      const char *name = "";
      switch( SpecUtils::SpectrumType(i) )
      {
        case SpecUtils::SpectrumType::Foreground:       name = "Foreground";      break;
        case SpecUtils::SpectrumType::SecondForeground: name = "Second Spectrum"; break;
        case SpecUtils::SpectrumType::Background:       name = "Background";      break;
      }//switch( SpecUtils::SpectrumType(i) )
      
      WRadioButton *button = new WRadioButton( name, m_buttonbox );
      button->setMargin( 10, Wt::Right );
      button->setInline(false);
      m_buttonGroup->addButton( button, i );
    }//for( int i = 0; i < 3; ++i )
    
    m_buttonGroup->setSelectedButtonIndex( 0 );
    m_buttonbox->setHidden(true);
    layout->addWidget( m_buttonbox, ++row, 0  );
    
    
    //m_deleteButton = new WPushButton( "", footer );
    //m_deleteButton->setHidden(true);
    //m_deleteButton->setIcon( "InterSpec_resources/images/minus_min.png" );
    //m_deleteButton->addStyleClass("FloatLeft");
    //m_deleteButton->clicked().connect( this, &SnapshotBrowser::startDeleteSelected );
    
    //m_renameButton
    
    //m_renameButton
    m_loadSpectraButton = new WPushButton( "Load Spectrum Only", footer );
    m_loadSpectraButton->clicked().connect( boost::bind( &SnapshotBrowser::loadSpectraSelected, this));
    m_loadSpectraButton->disable();
    
    m_loadSnapshotButton = new WPushButton( "Load App State", footer);
    m_loadSnapshotButton->clicked().connect( boost::bind(&SnapshotBrowser::loadSnapshotSelected, this));
    m_loadSnapshotButton->setDefault(true);
    //m_loadSnapshotButton->setIcon( "InterSpec_resources/images/time.svg" );
    m_loadSnapshotButton->disable();
    
    if( !buttonBar )
    {
      layout->addWidget( footer, layout->rowCount()+1 , 0 );
      m_loadSpectraButton->setFloatSide( Wt::Side::Right );
      m_loadSnapshotButton->setFloatSide( Wt::Side::Right );
    }
  }catch( std::exception &e )
  {
    if( !buttonBar )
      delete footer;
    
    // Remove and clear all children widgets so we dont leak memory
    clear();
    
    m_loadSnapshotButton = nullptr;
    m_loadSpectraButton = nullptr;
    m_renameButton = nullptr;
    m_buttonGroup = nullptr;
    m_buttonbox = nullptr;
    m_snapshotTable = nullptr;
    m_descriptionLabel = nullptr;
    m_timeLabel = nullptr;
    m_relatedSpectraLayout = nullptr;
    //m_header = header ),
    //m_finished( this ),
    m_editWindow = nullptr;
    m_nrows = 0;
    
    WText *txt = new WText( "Unexpected error creating database spectrum file browser", this );
    txt->setAttributeValue( "style", "color: red; font-weight: bold; font-size: 22px;" );
    
    passMessage( (string("Error creating database spectrum file browser: ") + e.what()).c_str(),
                WarningWidget::WarningMsgHigh );
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
  
  for( int i = 0; i < 3; i++ )
  {
    //loop through each SpecUtils::SpectrumType because sometimes there can be one spectra that is in multiple types
    Spectras spectras = m_session->session()->find< UserFileInDb >()
                                  .where( " id = ? OR id = ? OR id = ?")
                                  .bind( (*versionIterator)->foregroundId )
                                  .bind( (*versionIterator)->backgroundId )
                                  .bind( (*versionIterator)->secondForegroundId );
    
    string pre = "";
    string post = "";
    SpecUtils::SpectrumType spectratype = SpecUtils::SpectrumType(i);
    
    for( Dbo::collection< Dbo::ptr<UserFileInDb> >::const_iterator spectraIterator = spectras.begin();
        spectraIterator != spectras.end(); ++spectraIterator )
    {
      Wt::WIconPair *icon = nullptr;
      if ((*versionIterator)->foregroundId==spectraIterator->id() && spectratype==SpecUtils::SpectrumType::Foreground)
      {
        //icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_forwards.png","InterSpec_resources/images/shape_move_forwards.png");
      }else if ((*versionIterator)->backgroundId==spectraIterator->id() && spectratype==SpecUtils::SpectrumType::Background)
      {
        //icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_backwards.png","InterSpec_resources/images/shape_move_backwards.png");
      }else if ((*versionIterator)->secondForegroundId==spectraIterator->id() && spectratype==SpecUtils::SpectrumType::SecondForeground)
      {
        //icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_front.png","InterSpec_resources/images/shape_move_front.png");
      }else
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
      
      const auto nodeTxt = pre + " " + (*spectraIterator)->filename
                           + " <i>(" + descriptionText(spectratype)+")</i>" + post;
      Wt::WTreeNode *spectraNode = new Wt::WTreeNode( nodeTxt, icon, versionNode );
      
      // Hack in letting user double click on row to load just the spectrum
      auto label = spectraNode->label();
      WContainerWidget *rowDiv = label ? dynamic_cast<WContainerWidget *>( label->parent() ) : nullptr;
      if( rowDiv )
      {
        rowDiv->doubleClicked().connect( std::bind([this,spectraNode](){
          const set<WTreeNode *> sets = m_snapshotTable->selectedNodes();
          if( !sets.empty() && ((*begin(sets)) == spectraNode) )
            loadSpectraSelected();
#if( PERFORM_DEVELOPER_CHECKS )
          else
            log_developer_error( __func__, "Got double click on SnapshotBrowser spectrum but that wasnt what was selected." );
#endif
        }) );
      }//if( rowDive )
      
      
      
      if( m_header && !m_header->m_uuid.empty() )
      {
        //disable all child spectra tree nodes, so the user can only select the snapshot
        //spectraNode->disable();
      } //!m_uuid.empty()
      
      spectraNode->setToolTip((*spectraIterator)->serializeTime.toString()  + (((*spectraIterator)->description).empty()?"":(" -- " + (*spectraIterator)->description)));
      //Add to lookup table (also the HEAD)
      m_UserFileInDbLookup[spectraNode]=*spectraIterator;
      break;
    } //iterate through the spectra contained in this snapshot
  } //for( int i = 0; i < 3; i++ ) -- loop through SpecUtils::SpectrumType
  
} //void SnapshotBrowser::addSpectraNodes(Dbo::collection< Dbo::ptr<UserState> >::const_iterator versionIterator, Wt::WTreeNode *versionNode)

//Updates the buttons when a row is selected
void SnapshotBrowser::selectionChanged()
{
  std::set<Wt::WTreeNode *> sets = m_snapshotTable->selectedNodes();
  
  Wt::WTreeNode *selectedTreeNode = nullptr;
  if( !sets.empty() ) //Should only select one row at a time
    selectedTreeNode = *sets.begin();
  
  Wt::Dbo::ptr<UserState> userstate;
  Wt::Dbo::ptr<UserFileInDb> dbfile;
  const auto statepos = m_UserStateLookup.find(selectedTreeNode);
  if( statepos != end(m_UserStateLookup) )
    userstate = statepos->second;
  
  const auto filepos = m_UserFileInDbLookup.find(selectedTreeNode);
  if( filepos != end(m_UserFileInDbLookup) )
    dbfile = filepos->second;
  
  if( userstate ) //UserState clicked
  {
    m_buttonbox->disable();
    m_buttonbox->hide();
    //m_deleteButton->show();
    m_loadSnapshotButton->enable();
    m_loadSpectraButton->disable();
    m_descriptionLabel->setText("<i>"+ userstate->description + "</i>");
    m_timeLabel->setText("<b>" + userstate->serializeTime.toString() + "</b>");
    
    if( m_editWindow )
    {
      if( m_editWindow )
        AuxWindow::deleteAuxWindow(m_editWindow);
      m_editWindow = nullptr;
      //else //is edit dialog, in which case just delete the dialog
    }
  }else if( dbfile ) //UserFileDB clicked
  {
    m_buttonbox->show();
    m_buttonbox->enable();
    m_loadSpectraButton->enable();
    //m_deleteButton->hide();
    m_loadSnapshotButton->disable();
    m_descriptionLabel->setText("<i>" + dbfile->description + "</i>");
    m_timeLabel->setText("<b>" + dbfile->serializeTime.toString() + "</b>");
    
    if( m_editWindow )
    {
      AuxWindow::deleteAuxWindow(m_editWindow);
      m_editWindow = nullptr;
    }
  }else //some node not found
  {
    m_buttonbox->hide();
    m_buttonbox->disable();
    m_loadSnapshotButton->disable();
    //m_deleteButton->hide();
    m_loadSpectraButton->disable();
    m_descriptionLabel->setText( ns_empty_descrip_label );
    m_timeLabel->setText( ns_empty_descrip_label );
  } //some node not found
  
  if( (m_header && !m_header->m_uuid.empty()) ||  !m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) )
  {
    m_buttonbox->hide(); //make sure it is hidden if no uuid selected
  }//!m_uuid.empty()
}//void selectionChanged()


void SnapshotBrowser::startDeleteSelected()
{
  const set<WTreeNode *> &selection = m_snapshotTable->selectedNodes();
  if( selection.empty() )
  {
    if( m_editWindow )
      AuxWindow::deleteAuxWindow(m_editWindow);
    m_editWindow = nullptr;
        
    return;
  }//if( selection.size() != 1 )
  
  const char *title = "Confirm Delete";
  
  WTreeNode *node = *begin(selection);
  WText *label = node->label();
  const string name = label ? label->text().toUTF8() : string("");
  
  const string msg =
  "<div>Are you sure you want to delete:</div>"
  "<div style=\"white-space: nowrap; text-align: center; text-overflow: ellipsis; overflow: hidden; font-weight: bold; max-width:98%;\">" + name + "</div>"
  "<div>?</div>";
  
  // We will rely on the SimpleDialog covering everything else to know that the selection didnt change or anything...
  SimpleDialog *dialog = new SimpleDialog( title, msg );
  
  dialog->addButton( "No" );
  WPushButton *yes = dialog->addButton( "Yes" );
  yes->clicked().connect( this, &SnapshotBrowser::deleteSelected );
  
  //Need to update text when selection changes, currently relying on modal underlay to protect against this.
  //m_snapshotTable->itemSelectionChanged().connect(<#WObject *target#>, <#WObject::Method method#>)
}//void startDeleteSelected()


void SnapshotBrowser::startEditSelected()
{
  if( !m_session )
  {
    passMessage( "No active database session", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const set<WTreeNode *> &selection = m_snapshotTable->selectedNodes();
  if( selection.empty() )
  {
    if( m_editWindow )
    {
      AuxWindow::deleteAuxWindow(m_editWindow);
      m_editWindow = nullptr;
    }
    return;
  }//if( selection.size() != 1 )
  
  WTreeNode *node = *begin(selection);
  Wt::Dbo::ptr<UserState> state;
  const auto pos = m_UserStateLookup.find(node);
  if( pos != end(m_UserStateLookup) )
   state = pos->second;
  
  if( !state )  //shouldnt ever happend
  {
    passMessage( "Invalid state selected", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const char *title = "Edit Title/Description";
  if( m_editWindow )
    AuxWindow::deleteAuxWindow( m_editWindow );
  
  m_editWindow = new AuxWindow( title,
                               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                | AuxWindowProperties::DisableCollapse | AuxWindowProperties::PhoneNotFullScreen) );
  
  m_editWindow->setWidth( std::min(425, std::max(m_viewer->renderedWidth(), 250)) );
  m_editWindow->setHeight( std::min(250, std::max(m_viewer->renderedHeight(), 150)) );
  
  WGridLayout *layout = m_editWindow->stretcher();
  
  WLabel *label = new WLabel( "Name" );
  WLineEdit *nameEdit = new WLineEdit();
  nameEdit->setAttributeValue( "ondragstart", "return false" );
  nameEdit->setEmptyText( "(Name to store under)" );
  nameEdit->setText( state->name );
  layout->addWidget( label, 0, 0 );
  layout->addWidget( nameEdit,  0, 1 );
  
  label = new WLabel( "Desc." );
  WTextArea *description = new WTextArea();
  description->setEmptyText( "(Optional description)" );
  description->setText( state->description );
  layout->addWidget( label, 1, 0 );
  layout->addWidget( description, 1, 1 );
  
  layout->setColumnStretch( 1, 1 );
  layout->setRowStretch( 1, 1 );
  
  WContainerWidget *foot = m_editWindow->footer();
  
  WPushButton *cancel = new WPushButton( "Cancel", foot );
  cancel->addStyleClass( "DialogClose" );
  cancel->setFloatSide( Wt::Side::Right );
  
  WPushButton *yes = new WPushButton( "Accept", foot );
  yes->addStyleClass( "DialogClose" );
  yes->setFloatSide( Wt::Side::Right );
  
  cancel->clicked().connect( m_editWindow, &AuxWindow::hide );
  yes->clicked().connect( std::bind([=](){
    if( nameEdit->text().toUTF8().empty() ) //should happen
    {
      passMessage( "You must enter a name", WarningWidget::WarningMsgHigh );
      return;
    }
    
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_session );
      state.modify()->name = nameEdit->text();
      state.modify()->description = description->text();
      transaction.commit();
      if( node->label() )
        node->label()->setText( nameEdit->text() );
      m_descriptionLabel->setText( "<i>" + description->text().toUTF8() + "</i>" );
    }catch( ... )
    {
      passMessage( "Error modifying state - probably not changed", WarningWidget::WarningMsgHigh );
    }//try /catch
    
    m_editWindow->hide();
  }) );
  
  //Need to update text when selection changes, currently relying on modal underlay to protect against this.
  //m_snapshotTable->itemSelectionChanged().connect(<#WObject *target#>, <#WObject::Method method#>)
  
  yes->setDisabled( true );
  auto haschanged = [=](){
    yes->setEnabled( !nameEdit->text().toUTF8().empty() && (nameEdit->text()!=state->name || description->text()!=state->description) );
  };
  
  nameEdit->textInput().connect( std::bind(haschanged) );
  description->textInput().connect( std::bind(haschanged) );
  
  auto deleter = wApp->bind( boost::bind( &AuxWindow::deleteAuxWindow, m_editWindow ) );
  m_editWindow->finished().connect( std::bind( [deleter,this](){ deleter(); m_editWindow = nullptr; } ) );
  
  m_editWindow->show();
  m_editWindow->centerWindow();
}//void startEditSelected()



std::shared_ptr<SpecMeas> SnapshotBrowser::retrieveMeas( const int dbid )
  {
    std::shared_ptr<SpecMeas> snapshotmeas;
    if( dbid < 0 )
      return snapshotmeas;
    
    try
    {
      DataBaseUtils::DbTransaction transaction( *m_session );
      
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
      passMessage( "Sorry, couldnt load requested snapshot", WarningWidget::WarningMsgHigh );
    }//try / catch
    
    return snapshotmeas;
  }//std::shared_ptr<SpecMeas> retrieveMeas( const int dbid )
  

void SnapshotBrowser::deleteSelected()
{
  const set<WTreeNode *> &selection = m_snapshotTable->selectedNodes();
  if( selection.empty() )
  {
    //m_deleteButton->disable();
    passMessage( "No state selected to delete", WarningWidget::WarningMsgMedium );
    return;
  }//if( !indices.size() )
  
  WTreeNode * const node = *begin(selection);
  
  Wt::Dbo::ptr<UserState> state;
  
  const auto statepos = m_UserStateLookup.find(node);
  if( statepos != end(m_UserStateLookup) )
    state = statepos->second;
  
  if( !state )
  {
    passMessage( "Invalid state selected", WarningWidget::WarningMsgMedium );
    return;
  }
  
  try
  {
    if( m_viewer->currentAppStateDbId() == state.id() )
      m_viewer->resetCurrentAppStateDbId();
    
    {
      DataBaseUtils::DbTransaction transaction( *m_session );
      state.reread();  //is this actually necassary?
      state.remove();
      transaction.commit();
    }
    string name;
    if( node->label() )
      name = node->label()->text().toUTF8();
    
    m_UserStateLookup.erase( m_UserStateLookup.find(node) );
    m_snapshotTable->treeRoot()->removeChildNode(node);
    m_nrows -= 1;
    
    passMessage( "Snapshot '" + name + "' removed from database.", WarningWidget::WarningMsgInfo );
  }catch( ... )
  {
    passMessage( "Error removing snapshot from database.", WarningWidget::WarningMsgHigh );
  }
  
  if( m_editWindow )
  {
    AuxWindow::deleteAuxWindow(m_editWindow);
    m_editWindow = nullptr;
  }//if( m_editWindow )
}//void deleteSelected()


void SnapshotBrowser::loadSnapshotSelected()
{
  std::set<Wt::WTreeNode *> sets = m_snapshotTable->selectedNodes();
  
  //Only one row selected
  Wt::WTreeNode * selectedTreeNode = *sets.begin();
  
  Dbo::ptr<UserState> dbstate;
  const auto statepos = m_UserStateLookup.find(selectedTreeNode);
  if( statepos != end(m_UserStateLookup) )
    dbstate = statepos->second;
  
  if( !dbstate )
  {
    m_loadSnapshotButton->disable();
    //m_deleteButton->hide();
    return;
  }//if( !indices.size() )
  
  
  
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
    passMessage( "Error loading from the database", WarningWidget::WarningMsgHigh );
    m_finished.emit();
    return;
  }//if( !selected )
  
  try
  {
    m_viewer->loadStateFromDb( dbstate );
    
    WString msg = "Snapshot '";
    msg += dbstate.get()->name;
    msg += "' loaded";
    
    passMessage( msg.toUTF8(), WarningWidget::WarningMsgInfo );
  }catch( std::exception &e )
  {
    passMessage( "Failed to load state", WarningWidget::WarningMsgHigh );
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
        passMessage( "Error loading from the database", WarningWidget::WarningMsgHigh );
        m_finished.emit();
        return;
    }//if( !selected )

    if( !m_header || m_header->m_uuid.empty() )
    {
        
        int snapshot_id = -1;
        
        SpecUtils::SpectrumType type = SpecUtils::SpectrumType::Foreground;
        if( m_buttonbox->isVisible() )
            type = SpecUtils::SpectrumType( m_buttonGroup->checkedId() );
        
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
                            m_manager->displayFile( row, meas, type, false, false, SpecMeasManager::VariantChecksToDo::None );
                        }//if( meas )
                    }//if( snapshot_id )
                    
                    passMessage("Already opened", WarningWidget::WarningMsgInfo);
                    m_finished.emit();
                    return;
                } //measurement == m_viewer->measurment(type)
                break;
            } //if( dbentry.id() == dbfile.id() )
        }//for( int row = 0; row < specmodel->rowCount(); ++row )
        
        if( modelrow < 0 )
            modelrow = m_manager->setDbEntry( dbfile, header, measurement, true );
        
        m_manager->displayFile( modelrow, measurement, type, false, false, SpecMeasManager::VariantChecksToDo::None );
    }
    else
    {
        //Copied from 
        //dorevert() is only called from within the application loop
        SpectraFileModel *m_model = m_manager->model();
        if( !dbfile || !m_header || !m_model || !wApp )
            throw runtime_error( "DbSpecFileItem: invalid input or no wApp" );
        
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
                    m_manager->displayFile( row, measurement, SpecUtils::SpectrumType::Foreground, false, false, SpecMeasManager::VariantChecksToDo::None );
                    m_finished.emit();
                    return;
                }//if( entry.id() == dbfile.id() )
            }//for( int row = 0; row < m_model->rowCount(); ++row )
            
            try
            {
                const int modelRow = m_manager->setDbEntry( dbfile, header,
                                                           measurement, true );
                m_manager->displayFile( modelRow, measurement, SpecUtils::SpectrumType::Foreground, false, false, SpecMeasManager::VariantChecksToDo::None );
                
            }catch( exception &e )
            {
                cerr << "\n\nSnapshotModel\n\tCaught: " << e.what() << "\n\n";
                passMessage( "Error displaying previous measurment, things may not"
                            " be as expected" , WarningWidget::WarningMsgHigh );
            }//try / catch
        }else
        {
            m_header->setDbEntry( dbfile );
        }//if( dbfile->userHasModified )
    }
  
  m_finished.emit();
}//void loadSpectraSelected


int SnapshotBrowser::numSnaphots() const
{
  return m_nrows;
};


AutosavedSpectrumBrowser::AutosavedSpectrumBrowser(
                           const std::vector<Wt::Dbo::ptr<UserFileInDb>> &modifiedFiles,
                           SpecUtils::SpectrumType type,
                           SpectraFileModel *fileModel,
                           SpecMeasManager *manager,
                           std::shared_ptr<SpectraFileHeader> header,
                           WContainerWidget *parent )
: WContainerWidget( parent ),
  m_loadedASpectrum( this )
{
  wApp->useStyleSheet( "InterSpec_resources/DbFileBrowser.css" );
  
  addStyleClass( "AutosavedSpectrumBrowser" );
  WGridLayout *layout = new WGridLayout( this );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  
  const char *msg = "";
  if( header )
    msg = "It looks like you have previously loaded and modified this spectrum file, would you like"
          " to resume your previous work?";
  else
    msg = "Select your previously used spectrum file to continue from.";
  
  WText *txt = new WText( msg, Wt::XHTMLText);
  txt->addStyleClass( "ResumeWorkTxt" );
  txt->setInline( false );
  
  layout->addWidget(txt,0,0);
  
  WContainerWidget *container = new WContainerWidget();
  container->addStyleClass( "AutosavedSpectra" );
  layout->addWidget( container, 1, 0 );
  layout->setRowStretch( 1, 1 );
  
  for( size_t i = 0; i < modifiedFiles.size(); ++i )
  {
    if( !!modifiedFiles[i] )
    {
      auto entry = new DbSpecFileItem( type, fileModel, manager, modifiedFiles[i], header, container );
      entry->m_loadedASelection.connect( this, &AutosavedSpectrumBrowser::entryWasLoaded );
    }
  }//for( size_t i = 0; i < unModifiedFiles.size(); ++i )
}//AutosavedSpectrumBrowser
 

void AutosavedSpectrumBrowser::entryWasLoaded()
{
  m_loadedASpectrum.emit();
}

Wt::Signal<> &AutosavedSpectrumBrowser::loadedASpectrum()
{
  return m_loadedASpectrum;
}
