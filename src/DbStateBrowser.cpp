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
#include <Wt/WTreeView>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/Dbo/QueryModel>


#include "InterSpec/AuxWindow.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DbStateBrowser.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/RowStretchTreeView.h"


using namespace std;
using namespace Wt;

DbStateBrowser::DbStateBrowser( InterSpec *viewer, bool testStatesOnly )
  : AuxWindow( "Restore Previously Saved Snapshot", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::DisableCollapse | AuxWindowProperties::EnableResize) ),
    m_viewer( viewer ),
    m_model( 0 ),
    m_table( 0 ),
    m_loadButton( 0 )
  {
    int width = 500, height = 375;
    width = std::min( width, static_cast<int>(0.95*viewer->renderedWidth()) );
    height = std::min( height, static_cast<int>(0.95*viewer->renderedHeight()) );
    
    resize( WLength(width), WLength(height) );
    
    //We have to create a independant Dbo::Session for this class since the
    //  m_viewer->m_user.session() could be used in other threads, messing
    //  things up (Dbo::Session is not thread safe).
    m_session.reset( new DataBaseUtils::DbSession( *m_viewer->sql() ) );
    
    try
    {
      WGridLayout *layout = stretcher();
      Dbo::ptr<InterSpecUser> user = m_viewer->m_user;
      
      m_table = new RowStretchTreeView();
      m_table->setRootIsDecorated	(	false); //makes the tree look like a table! :)
      
      m_table->addStyleClass( "DbSpecFileSelectTable" );
      m_model = new Dbo::QueryModel< Dbo::ptr<UserState> >( m_table );
      
      WStringStream query;
      query << "InterSpecUser_id = " << user.id();
      if( testStatesOnly )
        query << " AND StateType = " << int(UserState::kForTest);
      else
        query << " AND StateType <> " << int(UserState::kForTest);
      m_model->setQuery( m_session->session()->find<UserState>().where( query.str() ) );
      m_model->addColumn( "Name" );
      m_model->addColumn( "Description" );
      m_model->addColumn( "SerializeTime" );
      m_model->addColumn( "id" );
      
      //        m_model->setColumnFlags( 0, ItemIsEditable | ItemIsSelectable );
      //        m_model->setColumnFlags( 1, ItemIsEditable | ItemIsSelectable );
      
      m_model->setHeaderData(  0, Horizontal, WString("Snapshot"), DisplayRole );
      m_model->setHeaderData(  1, Horizontal, WString("Description"), DisplayRole );
      m_table->setColumnWidth( 1, 130 );
      m_model->setHeaderData(  2, Horizontal, WString("Save Time"), DisplayRole );
      m_table->setColumnWidth( 2, 130 );
      m_model->setHeaderData(  2, Horizontal, WString("State ID"), DisplayRole );
      m_table->setColumnWidth( 3, 20 );
      
      LocalTimeDelegate *dtDelegate = new LocalTimeDelegate( m_table );
      m_table->setItemDelegateForColumn( 2, dtDelegate );
      
      m_table->setModel( m_model );
      m_table->setAlternatingRowColors( true );
      m_table->sortByColumn( 1, Wt::AscendingOrder );
      m_table->setSelectionMode( SingleSelection );
      
      for( int col = 0; col < m_model->columnCount(); ++col )
      {
        m_table->setColumnHidden( col, false );
        m_table->setSortingEnabled( col, true );
      }//for( int col = 0; col < model->columnCount(); ++col )
      
      m_table->sortByColumn( 2, Wt::DescendingOrder );
      
      m_table->selectionChanged().connect( this, &DbStateBrowser::selectionChanged );
      m_table->doubleClicked().connect( this, &DbStateBrowser::handleDoubleClicked );
      
      layout->addWidget( m_table, 0, 0 );
      layout->setRowStretch( 0, 1 );
      
      
    
//      footer()->resize( WLength::Auto, WLength(50.0) );
      
      m_loadButton = new WPushButton( "Restore", footer() );
      m_loadButton->setStyleClass("DatabaseGoIcon");
      m_loadButton->clicked().connect( this, &DbStateBrowser::loadSelected );
      m_loadButton->setFloatSide(Right);
      m_loadButton->disable();
      
      WPushButton *cancel = addCloseButtonToFooter();
      cancel->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
      
      rejectWhenEscapePressed();
      finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
    }catch( std::exception &e )
    {
      passMessage( "Error creating database spectrum file browser",
                  "", WarningWidget::WarningMsgHigh );
      cerr << "\n\nDbFileBrowser caught: " << e.what() << endl << endl;
    }//try / catch
    
    resizeToFitOnScreen();
    centerWindow();
    show();


  }//DbFileBrowser
  
  DbStateBrowser::~DbStateBrowser()
  {
  }
  
  void DbStateBrowser::selectionChanged()
  {
    if( m_table->selectedIndexes().size() )
      m_loadButton->enable();
  }//void selectionChanged()
  
  void DbStateBrowser::handleDoubleClicked( WModelIndex index, WMouseEvent event )
  {
    if( event.button() != WMouseEvent::LeftButton || !index.isValid() )
      return;
    WModelIndexSet selected;
    selected.insert( index );
    m_table->setSelectedIndexes( selected );
    loadSelected();
  }//void handleDoubleClicked( WModelIndex index, WMouseEvent event )
  
  void DbStateBrowser::loadSelected()
  {
    WModelIndexSet indices = m_table->selectedIndexes();
    if( !indices.size() )
    {
      m_loadButton->disable();
      return;
    }//if( !indices.size() )
    
    WModelIndex index = *indices.begin();
    
    Dbo::ptr<UserState> dbstate = m_model->stableResultRow( index.row() );
    
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
      delete this;
      return;
    }//if( !selected )
    
    try
    {
      m_viewer->loadStateFromDb( dbstate );
    }catch( std::exception &e )
    {
      passMessage( "Failed to load state", "", WarningWidget::WarningMsgHigh );
      cerr << "DbStateBrowser::loadSelected() caught: " << e.what() << endl;
    }//try / catch
    
    delete this;
  }//void loadSelected()
  
//};//class DbStateBrowser
  //#endif //#if( USE_DB_TO_STORE_SPECTRA )
