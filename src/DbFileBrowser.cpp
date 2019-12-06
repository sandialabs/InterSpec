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
#include <Wt/WTabWidget>
#include <Wt/WModelIndex>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>


#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DbFileBrowser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RowStretchTreeView.h"


using namespace std;
using namespace Wt;


DbFileBrowser::DbFileBrowser( SpecMeasManager *manager,
                              InterSpec *viewer,
                              SpectrumType type,
                              std::shared_ptr<SpectraFileHeader> header)
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
    //m_deleteButton( nullptr ),
    m_renameButton( nullptr ),
    m_buttonGroup( nullptr ),
    m_buttonbox( nullptr ),
    m_snapshotTable( nullptr ),
    m_descriptionLabel( nullptr ),
    m_timeLabel( nullptr ),
    m_relatedSpectraLayout( nullptr ),
    m_type( type ),
    m_header( header ),
    m_finished( this ),
    m_editWindow( nullptr ),
    m_nrows( 0 )
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
  
  int row = 0;
  if( m_header && !m_header->m_uuid.empty() )
  {
    auto txt = new WText("It looks like you have previously loaded and modifed this spectrum file, would you like to resume your previous work?");
    layout->addWidget( txt, row++, 0);
  }//if( we want to load a specific state )
  
  try
  {
    Dbo::ptr<InterSpecUser> user = m_viewer->m_user;
    WContainerWidget* tablecontainer = new WContainerWidget();
    
    //With m_snapshotTable being a WTree, we cant implement double clicking on
    //  an item to open it right away.  If we use a WTreeView, then we should
    //  be able to do that...
    m_snapshotTable = new WTree();
    WGridLayout *tablelayout = new WGridLayout();
    tablelayout->setContentsMargins(2, 2, 2, 2);
    tablelayout->setRowStretch(0, 1);
    tablelayout->setColumnStretch(0, 1);
    tablelayout->addWidget(m_snapshotTable, 0, 0);
    tablecontainer->setLayout( tablelayout );
    tablecontainer->setMargin(0);
    tablecontainer->setOffsets(0);
    
    //m_snapshotTable->resize(400,WLength(WLength::Auto));
    
    m_snapshotTable->itemSelectionChanged().connect(boost::bind( &SnapshotBrowser::selectionChanged, this ) );
    
    Wt::WTreeNode *root = new Wt::WTreeNode( "Root" );
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
    
    m_nrows = query.size(); //save for future query
    
    if( m_nrows == 0 )
    {
      if( m_header && !m_header->m_uuid.empty() )
      {
        layout->addWidget( new WLabel("No saved application states"), row++, 0 );
      }else
      {
        //This prevents the case where the checkIfPreviouslyOpened() is true,
        //  but the previous opened file is a kEndofSessionTemp, so will not
        //  show in the dialog.
        throw std::runtime_error( "Nothing to display" );
      }//m_uuid.empty()
    }//if( m_nrows == 0 )
    
    for( Dbo::collection< Dbo::ptr<UserState> >::const_iterator snapshotIterator = query.begin();
        snapshotIterator != query.end(); ++snapshotIterator )
    {
      Wt::WTreeNode *snapshotNode = new Wt::WTreeNode((*snapshotIterator)->name, 0, root);
      
      //m_snapshotTable->doubleClicked
      
      WContainerWidget *rowDiv = nullptr;
      if( buttonBar && snapshotNode->label() && snapshotNode->label()->parent() )
        rowDiv = dynamic_cast<WContainerWidget *>( snapshotNode->label()->parent() );
      
      if( rowDiv )  //Always seems to be true for Wt 3.1.4
      {
        WImage *delBtn = new WImage( "InterSpec_resources/images/minus_min_black.svg", rowDiv );
        delBtn->resize( WLength(16), WLength(16) );
        delBtn->addStyleClass( "DelSnapshotBtn" );
        delBtn->clicked().connect( this, &SnapshotBrowser::startDeleteSelected );
        
        
        WImage *editBtn = new WImage( "InterSpec_resources/images/edit_black.svg", rowDiv );
        editBtn->resize( WLength(16), WLength(16) );
        editBtn->addStyleClass( "DelSnapshotBtn" );
        editBtn->clicked().connect( this, &SnapshotBrowser::startEditSelected );
        
        if( !wApp->styleSheet().isDefined("snapshotbtn") )
        {
          wApp->styleSheet().addRule( ".DelSnapshotBtn", "float: right; background: none; display: none; cursor: pointer; opacity: 0.4; margin-top: 1px; margin-left: 1px; margin-right: 2px", "snapshotbtn" );
          wApp->styleSheet().addRule( ".DelSnapshotBtn:hover", "opacity: 1;" );
          wApp->styleSheet().addRule( ".Wt-selected .DelSnapshotBtn", "display: inline;" );
        }
      }else
      {
        //In case Wt changes in new versions
        cerr << "\n\nsnapshotNode->label()->parent() Is Not a WContainerWidget\n\n" << endl;
      }//if( rowDiv )
      
      //snapshotNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
      snapshotNode->setChildCountPolicy(Wt::WTreeNode::Disabled);
      snapshotNode->setToolTip((*snapshotIterator)->serializeTime.toString()
                               + (((*snapshotIterator)->description).empty() ? "" : (" -- " + (*snapshotIterator)->description)));
      
      //Add to lookup table (the HEAD)
      m_UserStateLookup[snapshotNode] = *snapshotIterator;
      
      //If not root, then return how many versions
      typedef Dbo::collection< Dbo::ptr<UserState> > Snapshots;
      typedef Snapshots::iterator SnapshotIter;
      Snapshots snapshots;
      
      if( *snapshotIterator )
      {
        //Create HEAD
        Wt::WTreeNode *versionNode = new Wt::WTreeNode("HEAD", new Wt::WIconPair("InterSpec_resources/images/time.png","InterSpec_resources/images/time.png"), snapshotNode);
        
        versionNode->setToolTip((*snapshotIterator)->serializeTime.toString() + (((*snapshotIterator)->description).empty()?"":(" -- " + (*snapshotIterator)->description)));
        
        //versionNode->setChildCountPolicy(Wt::WTreeNode::Enabled);
        versionNode->setChildCountPolicy(Wt::WTreeNode::Disabled);
        
        //Add to lookup table (also the HEAD)
        m_UserStateLookup[versionNode] = *snapshotIterator;
        addSpectraNodes(snapshotIterator, versionNode);
        
        snapshots = m_session->session()->find<UserState>()
                             .where( "SnapshotTagParent_id = ? AND StateType = ? " )
                             .bind( snapshotIterator->id() )
                             .bind(int(UserState::kUserSaved))
                             .orderBy( "id desc" );
        
        for( Dbo::collection< Dbo::ptr<UserState> >::const_iterator versionIterator = snapshots.begin();
            versionIterator != snapshots.end(); ++versionIterator )
        {
          const Wt::WString &name = (*versionIterator)->name;
          auto icon = new Wt::WIconPair("InterSpec_resources/images/time.png","InterSpec_resources/images/time.png");
          versionNode = new Wt::WTreeNode( name, icon, snapshotNode );
          
          auto tooltip = (*versionIterator)->serializeTime.toString()
                         + (((*versionIterator)->description).empty() ? "" : (" -- " + (*versionIterator)->description));
          versionNode->setToolTip( tooltip );
          
          //versionNode->setChildCountPolicy( Wt::WTreeNode::Enabled );
          versionNode->setChildCountPolicy( Wt::WTreeNode::Disabled );
          
          m_UserStateLookup[versionNode] = *versionIterator;
          
          addSpectraNodes( versionIterator, versionNode );
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
  {
    //loop through each spectrumtype because sometimes there can be one spectra that is in multiple types
    Spectras spectras = m_session->session()->find< UserFileInDb >()
                                  .where( " id = ? OR id = ? OR id = ?")
                                  .bind( (*versionIterator)->foregroundId )
                                  .bind( (*versionIterator)->backgroundId )
                                  .bind( (*versionIterator)->secondForegroundId );
    
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
      }else if ((*versionIterator)->backgroundId==spectraIterator->id() && spectratype==kBackground)
      {
        icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_backwards.png","InterSpec_resources/images/shape_move_backwards.png");
      }else if ((*versionIterator)->secondForegroundId==spectraIterator->id() && spectratype==kSecondForeground)
      {
        icon = new Wt::WIconPair("InterSpec_resources/images/shape_move_front.png","InterSpec_resources/images/shape_move_front.png");
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
      bool isDelDialog = (m_editWindow->windowTitle().toUTF8().find("Delete") != string::npos);
      AuxWindow::deleteAuxWindow(m_editWindow);
      m_editWindow = nullptr;
      if( isDelDialog )
        startDeleteSelected();
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
    m_descriptionLabel->setText("");
    m_timeLabel->setText("");
  } //some node not found
  
  if( (m_header && !m_header->m_uuid.empty()) ||  !m_viewer->measurment( kForeground ) )
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
    {
      AuxWindow::deleteAuxWindow(m_editWindow);
      m_editWindow = nullptr;
    }
    return;
  }//if( selection.size() != 1 )
  
  const char *title = "Confirm Delete";
  if( m_editWindow )
    AuxWindow::deleteAuxWindow( m_editWindow );
  
  m_editWindow = new AuxWindow( title,
                                 (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                                  | AuxWindowProperties::DisableCollapse | AuxWindowProperties::PhoneModal) );
  
  WTreeNode *node = *begin(selection);
  WText *label = node->label();
  const string name = label ? label->text().toUTF8() : string("");
  
  WText *txt = new WText( "<div style=\"white-space:nowrap;\">Are you sure you want to delete:</div>"
                         "<div style=\"white-space:nowrap;text-align:center;\"><span style=\"font-weight:bold;\">" + name + "</span>?</div>" );
  m_editWindow->contents()->addWidget( txt );
  m_editWindow->contents()->setPadding( 12, Wt::Side::Left | Wt::Side::Right );
  m_editWindow->contents()->setPadding( 18, Wt::Side::Top | Wt::Side::Bottom );
  
  WContainerWidget *foot = m_editWindow->footer();
  
  WPushButton *cancel = new WPushButton( "No", foot );
  cancel->addStyleClass( "DialogClose" );
  cancel->setFloatSide( Wt::Side::Right );
  
  WPushButton *yes = new WPushButton( "Yes", foot );
  yes->addStyleClass( "DialogClose" );
  yes->setFloatSide( Wt::Side::Right );
  
  cancel->clicked().connect( m_editWindow, &AuxWindow::hide );
  yes->clicked().connect( this, &SnapshotBrowser::deleteSelected );
  
  //Need to update text when selection changes, currently relying on modal underlay to protect against this.
  //m_snapshotTable->itemSelectionChanged().connect(<#WObject *target#>, <#WObject::Method method#>)
  
  auto deleter = wApp->bind( boost::bind( &AuxWindow::deleteAuxWindow, m_editWindow ) );
  m_editWindow->finished().connect( std::bind( [deleter,this](){ deleter(); m_editWindow = nullptr; } ) );
  
  m_editWindow->show();
  m_editWindow->centerWindow();
}//void startDeleteSelected()


void SnapshotBrowser::startEditSelected()
{
  if( !m_session )
  {
    passMessage( "No active database session", "", WarningWidget::WarningMsgHigh );
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
    passMessage( "Invalid state selected", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const char *title = "Edit Title/Description";
  if( m_editWindow )
    AuxWindow::deleteAuxWindow( m_editWindow );
  
  m_editWindow = new AuxWindow( title,
                               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                                | AuxWindowProperties::DisableCollapse | AuxWindowProperties::PhoneModal) );
  
  m_editWindow->setWidth( std::min(425, std::max(m_viewer->renderedWidth(), 250)) );
  m_editWindow->setHeight( std::min(250, std::max(m_viewer->renderedHeight(), 150)) );
  
  WGridLayout *layout = m_editWindow->stretcher();
  
  WLabel *label = new WLabel( "Name" );
  WLineEdit *nameEdit = new WLineEdit();
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
      passMessage( "You must enter a name", "", WarningWidget::WarningMsgHigh );
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
      passMessage( "Error modifying state - probably not changed", "", WarningWidget::WarningMsgHigh );
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
  

void SnapshotBrowser::deleteSelected()
{
  const set<WTreeNode *> &selection = m_snapshotTable->selectedNodes();
  if( selection.empty() )
  {
    //m_deleteButton->disable();
    passMessage( "No state selected to delete", "", WarningWidget::WarningMsgMedium );
    return;
  }//if( !indices.size() )
  
  WTreeNode * const node = *begin(selection);
  
  Wt::Dbo::ptr<UserState> state;
  
  const auto statepos = m_UserStateLookup.find(node);
  if( statepos != end(m_UserStateLookup) )
    state = statepos->second;
  
  if( !state )
  {
    passMessage( "Invalid state selected", "", WarningWidget::WarningMsgMedium );
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
    
    passMessage( "Snapshot '" + name + "' removed from database.", "", WarningWidget::WarningMsgInfo );
  }catch( ... )
  {
    passMessage( "Error removing snapshot from database.", "", WarningWidget::WarningMsgHigh );
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


int SnapshotBrowser::numSnaphots() const
{
  return m_nrows;
};
