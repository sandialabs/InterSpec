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

#include <deque>
#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <fstream>
#include <numeric>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// Block out some warnings occurring in boost files.
#pragma warning(disable:4244) // warning C4244: 'initializing' : conversion from 'std::streamoff' to 'size_t', possible loss of data
#pragma warning(disable:4308) // warning C4308: negative integral constant converted to unsigned type
#pragma warning(disable:4355) // warning C4308: negative integral constant converted to unsigned type

#include <boost/any.hpp>
#include <boost/range.hpp>
#include <boost/system/error_code.hpp>

#include <Wt/WLink>
#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WBorder>
#include <Wt/WServer>
#include <Wt/WTextArea>
#include <Wt/WIconPair>
#include <Wt/WIOService>

#include <Wt/WGroupBox>
#include <Wt/WResource>
#include <Wt/WDateTime>
#include <Wt/WGroupBox>
#include <Wt/WComboBox>
#include <Wt/WBoxLayout>

#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WModelIndex>
#include <Wt/WFileUpload>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WEnvironment>
#include <Wt/WApplication>
#include <Wt/WProgressBar>
#include <Wt/WStandardItem>
#include <Wt/WSelectionBox>
#include <Wt/Http/Response>
#include <Wt/WItemDelegate>
#include <Wt/WBorderLayout>
#include <Wt/Dbo/QueryModel>
#include <Wt/WMemoryResource>
#include <Wt/WStringListModel>
#include <Wt/WContainerWidget>
#if( HAS_WTDBOMYSQL )
#include <Wt/Dbo/backend/MySQL>
#endif
#include <Wt/WStandardItemModel>
#include <Wt/WAbstractItemModel>
#include <Wt/WCssDecorationStyle>
#include <Wt/Dbo/backend/Sqlite3>

#if( SUPPORT_ZIPPED_SPECTRUM_FILES )
#include "InterSpec/ZipArchive.h"
#endif


#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeasManager.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/CanvasForDragging.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/RowStretchTreeView.h"
#include "SpecUtils/SpectrumDataStructs.h"

#include "InterSpec/HelpSystem.h"

#if( USE_DB_TO_STORE_SPECTRA )
#include "InterSpec/DbFileBrowser.h"
#endif
  
#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif

using namespace Wt;
using namespace std;

#if( USE_DB_TO_STORE_SPECTRA )
namespace
{
  class PreviousDbEntry : public WContainerWidget
  {
    //a class to display, and then select a previous database entry.  This class
    //  gets displayed in a AuxWindow when the user uploads a spectrum they
    //  have previously used and modified.  This class represents the
    //  previous session with the spectrum
  public:
    PreviousDbEntry( AuxWindow *dialog, WContainerWidget* container, SpectrumType type,
                     SpectraFileModel *model, SpecMeasManager *manager,
                     Dbo::ptr<UserFileInDb> dbentry,
                     std::shared_ptr<SpectraFileHeader> header )
    : WContainerWidget(), m_dialog( dialog ), m_type( type ),
      m_model( model ), m_manager( manager ), m_dbentry( dbentry ),
      m_header( header )
    {
      addStyleClass( "PreviousDbEntry" );
      if( dialog )
        container->addWidget( this );
      string msg = "Uploaded: "
                + dbentry->uploadTime.toString( DATE_TIME_FORMAT_STR ).toUTF8();
      if( dbentry->userHasModified )
        msg += ", was modified";
      WText *txt = new WText( msg, this );
      txt->addStyleClass( "PreviousDbEntryTxt" );
      WPushButton *button = new WPushButton( "Resume From", this );
      button->addStyleClass( "PreviousDbEntryButton" );
      button->clicked().connect( this, &PreviousDbEntry::dorevert );
      button->setFocus();
    }//PreviousDbEntry(...)
    
    void dorevert()
    {
      //dorevert() is only called from within the application loop
      
      if( !m_dbentry || !m_header || !m_model || !wApp )
        throw runtime_error( "PreviousDbEntry: invalid input or no wApp" );
      
      if( m_dbentry->userHasModified )
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
          if( entry && entry.id() == m_dbentry.id() )
          {
            measurement = header->parseFile();
            m_manager->displayFile( row, measurement, m_type, false, false, false );
            m_dialog->hide();
            
            return;
          }//if( entry.id() == m_dbentry.id() )
        }//for( int row = 0; row < m_model->rowCount(); ++row )
        
        try
        {
          const int modelRow = m_manager->setDbEntry( m_dbentry, header,
                                                     measurement, true );
          m_manager->displayFile( modelRow, measurement, m_type, false, false, false );
          m_dialog->hide();
        }catch( exception &e )
        {
          cerr << "\n\n" << SRC_LOCATION <<"\n\tCaught: " << e.what() << "\n\n";
          passMessage( "Error displaying previous measurment, things may not"
                      " be as expected" , "", WarningWidget::WarningMsgHigh );
        }//try / catch
      }else
      {
        m_header->setDbEntry( m_dbentry );
      }//if( m_dbentry->userHasModified )
      
      //      if( m_dialog )
      //        delete m_dialog;
    }//void dorevert()
    
    AuxWindow *m_dialog;
    SpectrumType m_type;
    SpectraFileModel *m_model;
    SpecMeasManager *m_manager;
    Dbo::ptr<UserFileInDb> m_dbentry;
    std::shared_ptr<SpectraFileHeader> m_header;
  };//class PreviousDbEntry
  

  void setHeadersDbEntry( std::shared_ptr<SpectraFileHeader> header, Wt::Dbo::ptr<UserFileInDb> entry )
  {
    header->setDbEntry( entry );
    wApp->triggerUpdate();
  }
}//namespace
#endif //#if( USE_DB_TO_STORE_SPECTRA )

namespace
{
  class FileUploadDialog : public AuxWindow
  {
    //Class used to upload spectrum files.  Specialization of AuxWindow
    //  needed to manage InterSpec::displayedSpectrumChanged() connection.
    
    Wt::Signals::connection m_specChangedConection;
    WFileUpload *m_fileUpload;
    SpecMeasManager *m_manager;
    SpectrumType m_type;
    
  public:
    FileUploadDialog( InterSpec *viewer,
                      SpecMeasManager *manager )
    : AuxWindow( "", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                       | AuxWindowProperties::TabletModal | AuxWindowProperties::DisableCollapse | AuxWindowProperties::PhoneModal ) ),
      m_fileUpload( 0 ),
      m_manager( manager ),
      m_type( kForeground )
    {
      setWindowTitle( "Select File To Open" );
      
      const bool noForeground = !viewer->measurment( kForeground );
      
      string instructions;
      if( noForeground )
      {
        instructions = string(viewer->isMobile() ?"Tap" : "Click")
                       + " below to choost a foreground spectrum file to open.";
      }else
      {
        instructions = "Select how you would like to open the spectrum file, and then ";
        instructions += (viewer->isMobile() ? "tap" : "click");
        instructions += " below to browse for the file.";
      }//if( noForeground )
      
      auto layout = stretcher();
      WText *txt = new WText( instructions );
      
      layout->addWidget( txt, layout->rowCount(), 0 );
      
      
      if( !noForeground )
      {
        WGroupBox *buttons = new WGroupBox( "Open file as:" );
      
        layout->addWidget( buttons, layout->rowCount(), 0 );
        
        WButtonGroup *group = new WButtonGroup( buttons );
        
        WRadioButton *foreground = new WRadioButton( "Forground", buttons );
        foreground->setInline( false );
        foreground->setChecked( true );
        group->addButton( foreground, 0 );
        
        WRadioButton *background = new WRadioButton( "Background", buttons );
        background->setInline( false );
        group->addButton( background, 1 );
        
        WRadioButton *secondary = new WRadioButton( "Secondary", buttons );
        secondary->setInline( false );
        group->addButton( secondary, 2 );
        
        group->checkedChanged().connect( std::bind( [this,group](){
          switch( group->checkedId() ) {
            case 0: m_type = kForeground; break;
            case 1: m_type = kBackground; break;
            case 2: m_type = kSecondForeground; break;
            default:  //shouldnt ever happen
              break;
          }
        } ) );
      }//if( !noForeground )

    
      m_fileUpload = new WFileUpload();
#if( BUILD_FOR_WEB_DEPLOYMENT )
      m_fileUpload->setProgressBar( new WProgressBar() );
#endif
    
      layout->addWidget( m_fileUpload, layout->rowCount(), 0, AlignMiddle | AlignCenter );
      layout->setRowStretch( layout->rowCount()-1, 1 );
      
      const char *msg = "";
      if( !viewer->isMobile() )
      {
        msg = "You can also drag and drop the file directly into the <br />"
              "application window as a quicker alternative.<br />"
              "<br />For more advanced file opening and<br />"
              "manipulation use the <b>File Manager</b>";
      }else
      {
        msg = "You can also directly open spectrum files in <em>InterSpec</em><br />"
              "from your email or file apps.";
      }

      WText *friendlyMsg = new WText( msg );
      friendlyMsg->setStyleClass( "startQuickUploadMsg" );
      
      layout->addWidget( friendlyMsg, layout->rowCount(), 0, AlignMiddle | AlignBottom );
      
      
      WPushButton *cancel = new WPushButton( "Cancel", footer() );
      
      cancel->clicked().connect( this, &FileUploadDialog::userCanceled );
      m_fileUpload->changed().connect( m_fileUpload, &Wt::WFileUpload::upload );
      m_fileUpload->uploaded().connect( this, &FileUploadDialog::finishUpload );
      m_fileUpload->fileTooLarge().connect( boost::bind( &FileUploadDialog::toLarge, this, _1 ) );
      m_specChangedConection = viewer->displayedSpectrumChanged().connect( this, &AuxWindow::emitReject );
      
      finished().connect( this, &FileUploadDialog::userCanceled );
      
      rejectWhenEscapePressed();
      
      show();
      
      if( m_isPhone )
      {
        resizeScaledWindow( 1.0, 1.0 );
        repositionWindow( 0, 0 );
      }else
      {
        resizeToFitOnScreen();
        centerWindow();
      }
    }//FileUploadDialog constructor
    
    void toLarge( const ::int64_t size_tried )
    {
      m_manager->fileTooLarge( size_tried );
      
      if( m_specChangedConection.connected() )
        m_specChangedConection.disconnect();
      delete this;
    }
    
    void userCanceled()
    {
      if( m_specChangedConection.connected() )
        m_specChangedConection.disconnect();
      delete this;
    }
    
    void finishUpload()
    {
      if( m_specChangedConection.connected() )
        m_specChangedConection.disconnect();
      
      m_manager->finishQuickUpload( m_fileUpload, m_type );
      delete this;
    }
    
    virtual ~FileUploadDialog()
    {
      if( m_specChangedConection.connected() )
        m_specChangedConection.disconnect();
    }//~FileUploadDialog()
    
  };//class FileUploadDialog
}//namespace



/* Import Spectrum Files dialog, used to upload foreground, background, 2nd foreground */

class UploadBrowser : public AuxWindow
{
public:
  UploadBrowser( SpecMeasManager *manager )
  : AuxWindow( "Import Spectrum Files",
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                | AuxWindowProperties::DisableCollapse
                | AuxWindowProperties::EnableResize) ),
  m_manager( manager )
  {
    setResizable( true );
  
    
    WGridLayout *layout = new Wt::WGridLayout();
    contents()->setLayout(layout);
    
    WText *uploadText = new WText( "Foreground: " );
    WFileUpload *m_fileUpload = new WFileUpload(  );
    m_fileUpload->changed().connect( m_fileUpload, &Wt::WFileUpload::upload );
    m_fileUpload->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager,m_fileUpload, kForeground));
    m_fileUpload->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge, _1 ) );

    WText *uploadText2 = new WText( "Background: " );
    WFileUpload *m_fileUpload2 = new WFileUpload(  );
    m_fileUpload2->changed().connect( m_fileUpload2, &Wt::WFileUpload::upload );
    m_fileUpload2->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager,m_fileUpload2, kBackground));
    m_fileUpload2->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge, _1 ) );
    
    WText *uploadText3 = new WText( "Secondary Foreground: " );
    WFileUpload *m_fileUpload3 = new WFileUpload(  );
    m_fileUpload3->changed().connect( m_fileUpload3, &Wt::WFileUpload::upload );
    m_fileUpload3->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager,m_fileUpload3, kSecondForeground));
    m_fileUpload3->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge, _1 ) );
    
    layout->addWidget( uploadText, 0, 0 );
    layout->addWidget( m_fileUpload, 0, 1 );
    layout->addWidget( uploadText2, 1, 0 );
    layout->addWidget( m_fileUpload2, 1, 1 );
    layout->addWidget( uploadText3, 2, 0 );
    layout->addWidget( m_fileUpload3, 2, 1 );
    layout->addWidget( new WText(""),3,0);
    
    
    WPushButton *cancel = addCloseButtonToFooter();
    cancel->clicked().connect( boost::bind( &AuxWindow::hide, this ) );

    layout->setRowStretch( 3, 1 );
    
    rejectWhenEscapePressed();
    finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
    
    
    rejectWhenEscapePressed();

    centerWindow();
    show();
  }//UploadBrowser
  
  virtual ~UploadBrowser()
  {
  }
  
//  void selectionChanged()
//  {
//    if( !m_model )
//      return;
//    
//    if( m_snapshotModel->rowCount() )
//      m_snapshotModel->removeRows( 0, m_snapshotModel->rowCount() );
//    
//    const WModelIndexSet indices = m_table->selectedIndexes();
//    
//    if( indices.empty() )
//    {
//      m_snapshot->hide();
//      m_loadButton->disable();
//      return;
//    }//if( !indices.empty() )
//    
//    m_loadButton->enable();
//    const WModelIndex index = *indices.begin();
//    Dbo::ptr<UserFileInDb> dbfile = m_model->stableResultRow( index.row() );
//    
//    try
//    {
//      Dbo::Transaction transaction( *m_session );
//      typedef Dbo::collection< Dbo::ptr<UserFileInDb> > Snapshots;
//      typedef Snapshots::iterator SnapshotIter;
//      Snapshots snapshots = m_session->find<UserFileInDb>()
//      .where( "SnapshotParent_id = ?" )
//      .bind( dbfile.id() );
//      
//      if( snapshots.size() )
//      {
//        m_snapshot->show();
//        m_snapshotModel->addString( "Latest Saved Version" );
//        m_snapshotModel->setData( 0, 0, -1, Wt::UserRole );
//      }else
//      {
//        m_snapshot->hide();
//      }//if( snapshots.size() ) / else
//      
//      int timeoffset = 0;
//      if( wApp )
//        timeoffset = wApp->environment().timeZoneOffset();
//      
//      for( SnapshotIter i = snapshots.begin(); i != snapshots.end(); ++i )
//      {
//        Dbo::ptr<UserFileInDb> snap = *i;
//        string desc = WString::fromUTF8( snap->description ).toUTF8() + " (";
//        desc += snap->serializeTime.addSecs( 60*timeoffset )
//        .toString( DATE_TIME_FORMAT_STR ).toUTF8() + ")";
//        if( desc.size() > 64 )
//          desc = desc.substr( 0, 61 ) + "...";
//        m_snapshotModel->addString( desc );
//        m_snapshotModel->setData( m_snapshotModel->rowCount()-1, 0,
//                                 boost::any(int(snap.id())), Wt::UserRole );
//      }//for( SnapshotIter i = snapshots.begin(); i != snapshots.end(); ++i )
//      
//      transaction.commit();
//    }catch( std::exception &e )
//    {
//      cerr << "UploadBrowser::selectionChanged() caught: "
//      << e.what() << endl;
//    }//try / catch
//  }//void selectionChanged()
  
//  void handleDoubleClicked( WModelIndex index, WMouseEvent event )
//  {
//    if( (event.button() != WMouseEvent::LeftButton) || !index.isValid() )
//      return;
//    
//    WModelIndexSet selected = m_table->selectedIndexes();
//    if( selected.size() && (*selected.begin() != index) )
//    {
//      selected.clear();
//      if( m_snapshotModel->rowCount() )
//        m_snapshotModel->removeRows( 0, m_snapshotModel->rowCount() );
//    }//if( selected.size() && (*selected.begin() != index ) ) / else
//    
//    selected.clear();
//    selected.insert( index );
//    
//    m_table->setSelectedIndexes( selected );
//    loadSelected();
//  }//void handleDoubleClicked( WModelIndex index, WMouseEvent event )
  
//  std::shared_ptr<SpecMeas> retrieveMeas( const int dbid )
//  {
//    std::shared_ptr<SpecMeas> snapshotmeas;
//    if( dbid < 0 )
//      return snapshotmeas;
//    
//    try
//    {
//      Dbo::Transaction transaction( *m_session );
//      Dbo::ptr<UserFileInDb> dbsnapshot = m_session->find<UserFileInDb>()
//      .where( "id = ?" )
//      .bind( dbid );
//      
//      Dbo::ptr<UserFileInDbData> data;
//      if( dbsnapshot && dbsnapshot->filedata.size() )
//        data = *(dbsnapshot->filedata.begin());
//      transaction.commit();
//      if( data )
//        snapshotmeas = data->decodeSpectrum();
//    }catch( std::exception &e )
//    {
//      cerr << "retrieveMeas() caught (while trying to load snaphot): "
//      << e.what() << endl;
//      passMessage( "Sorry, couldnt load requested snapshot",
//                  "", WarningWidget::WarningMsgHigh );
//    }//try / catch
//    
//    return snapshotmeas;
//  }//std::shared_ptr<SpecMeas> retrieveMeas( const int dbid )
  
//  void loadSelected()
//  {
//    if( !m_model )
//      return;
//    
//    WModelIndexSet indices = m_table->selectedIndexes();
//    if( !indices.size() )
//    {
//      m_loadButton->disable();
//      return;
//    }//if( !indices.size() )
//    
//    WModelIndex index = *indices.begin();
//    
//    Dbo::ptr<UserFileInDb> dbfile = m_model->stableResultRow( index.row() );
//    
//    //Now we have to make 'dbfile' be associated with same session of
//    //  m_viewer->m_user.session()  {I dont know what would happen otherwise)
//    if( dbfile.id() >= 0 )
//    {
//      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
//      DataBaseUtils::DbTransaction transaction( *sql );
//      dbfile = sql->session()->find< UserFileInDb >()
//                             .where( "id = ?").bind( dbfile.id() );
//      transaction.commit();
//    }//if( dbfile.id() >= 0 )
//    
//    if( !dbfile || dbfile.id() < 0 )
//    {
//      passMessage( "Error loading from the database",
//                  "", WarningWidget::WarningMsgHigh );
//      delete this;
//      return;
//    }//if( !selected )
//    
//    int snapshot_id = -1;
//    if( m_snapshotModel->rowCount() > 0
//       && m_snapshot->currentIndex() > 0 )
//    {
//      try
//      {
//        const int snapnum = m_snapshot->currentIndex();
//        const WModelIndex index = m_snapshotModel->index( snapnum, 0 );
//        boost::any val = m_snapshotModel->data( index, Wt::UserRole );
//        snapshot_id = boost::any_cast<int>( val );
//      }catch( std::exception &e )
//      {
//        cerr << "loadSelected() caught: " << e.what() << endl;
//      }//try / catch
//    }//if( m_snapshot->rowCount() > 0 )
//    
//    SpectrumType type = kForeground;
//    if( m_buttonGroup )
//      type = SpectrumType( m_buttonGroup->checkedId() );
//    
//    int modelrow = -1;
//    std::shared_ptr<SpecMeas> measurement;
//    std::shared_ptr<SpectraFileHeader> header;
//    
//#if( USE_DB_TO_STORE_SPECTRA )
//    //Now should check if already opened
//    SpectraFileModel *specmodel = m_manager->model();
//    for( int row = 0; row < specmodel->rowCount(); ++row )
//    {
//      std::shared_ptr<SpectraFileHeader> thisheader
//                                            = specmodel->fileHeader( row );
//      Wt::Dbo::ptr<UserFileInDb> dbentry = thisheader->dbEntry();
//      if( dbentry.id() == dbfile.id() )
//      {
//        modelrow = row;
//        dbfile = dbentry;
//        header = thisheader;
//        measurement = thisheader->parseFile();
//        if( measurement == m_viewer->measurment(type) )
//        {
//          if( snapshot_id >= 0 )
//          {
//            std::shared_ptr<SpecMeas> meas = retrieveMeas( snapshot_id );
//            if( meas )
//            {
//              thisheader->setMeasurmentInfo( meas );
//              m_manager->displayFile( row, meas, type, false, true );
//            }//if( meas )
//          }//if( snapshot_id )
//          
//          delete this;
//          return;
//        }
//        break;
//      }
//    }//for( int row = 0; row < specmodel->rowCount(); ++row )
//    
//    if( modelrow < 0 )
//      modelrow = m_manager->setDbEntry( dbfile, header, measurement, true );
//    
//    if( snapshot_id >= 0 && header )
//    {
//      std::shared_ptr<SpecMeas> snapshotmeas = retrieveMeas( snapshot_id );
//      if( snapshotmeas )
//      {
//        header->setMeasurmentInfo( snapshotmeas );
//        measurement = snapshotmeas;
//      }else
//      {
//        cerr << "Could not load spectrum file snapshot with id = "
//        << snapshot_id << endl;
//        passMessage( "Couldnt load requested snapshot",
//                    "", WarningWidget::WarningMsgHigh );
//      }//if( measurment ) / else
//    }//if( snapshot_id >= 0 && header )
//#endif //#if( USE_DB_TO_STORE_SPECTRA )
//    
//    m_manager->displayFile( modelrow, measurement, type, false, true );
//    
//    delete this;
//  }//void loadSelected()
  
protected:
//  std::shared_ptr<Wt::Dbo::Session> m_session;
  SpecMeasManager  *m_manager;
//  InterSpec   *m_viewer;
//  Dbo::QueryModel< Dbo::ptr<UserFileInDb> >    *m_model;
//  RowStretchTreeView       *m_table;
//  WPushButton      *m_loadButton;
//  WButtonGroup     *m_buttonGroup;
//  WComboBox        *m_snapshot;
//  WStringListModel *m_snapshotModel;
};//class UploadBrowser


SpecMeasManager::SpecMeasManager( InterSpec *viewer )
  : WObject(),
    m_treeView( NULL ),
    m_fileModel( NULL ),
    m_fileUpload( NULL ),
    m_viewer( viewer ),
//    m_setPrimaryButton( NULL ),
//    m_setSecondaryButton( NULL ),
//    m_setBackgroundButton( NULL ),
//    m_unDisplayMenuButton( NULL ),
//    m_unDisplayPopup( NULL ),
//    m_unDisplayPrimaryButton( NULL ),
//    m_unDisplaySecondaryButton( NULL ),
//    m_unDisplayBackgroundButton( NULL ),
//    m_makeNewFileButton( NULL ),
//    m_removeFileButton( NULL ),
//    m_saveFileAsButton( NULL ),
//    m_saveAsPopup( NULL ),
    m_setButton ( NULL),
    m_setAsForeground ( NULL),
    m_setAsBackground ( NULL),
    m_setAsSecForeground ( NULL),
    m_combineButton ( NULL),
    m_saveButton ( NULL),
    m_deleteButton ( NULL),
    m_removeForeButton ( NULL),
    m_removeBackButton ( NULL),
    m_removeFore2Button ( NULL)
#if( !ANDROID && !IOS )
    , m_foregroundDragNDrop( new FileDragUploadResource(this) ),
    m_secondForegroundDragNDrop( new FileDragUploadResource(this) ),
    m_backgroundDragNDrop( new FileDragUploadResource(this) )
#endif
    , m_destructMutex( new std::mutex() ),
    m_destructed( new bool(false) )
{
  wApp->useStyleSheet( "InterSpec_resources/SpecMeasManager.css" );
  
  for( SaveSpectrumAsType type = SaveSpectrumAsType(0);
      type < kNumSaveSpectrumAsType;
      type = SaveSpectrumAsType(type+1) )
  {
    m_downloadResources[type] = new DownloadSpectrumResource(type, this, this);
  }
  
//
    m_treeView = new RowStretchTreeView();
    m_fileModel = new SpectraFileModel( m_treeView );
    m_treeView->setModel( m_fileModel );
    m_treeView->selectionChanged().connect( boost::bind( &SpecMeasManager::selectionChanged, this ) );

    startSpectrumManager(); //initializes
    deleteSpectrumManager(); //deletes instance
    
  m_sql = viewer->sql();

  for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
       i < kNumSaveSpectrumAsType;
       i = SaveSpectrumAsType(i+1) )
    m_specificResources[i] = (SpecificSpectrumResource *)0;
  
#if( !ANDROID && !IOS )
  m_foregroundDragNDrop->fileDrop()->connect( boost::bind( &SpecMeasManager::handleFileDrop, this, _1, _2, kForeground ) );
  m_secondForegroundDragNDrop->fileDrop()->connect( boost::bind( &SpecMeasManager::handleFileDrop, this, _1, _2, kSecondForeground ) );
  m_backgroundDragNDrop->fileDrop()->connect( boost::bind( &SpecMeasManager::handleFileDrop, this, _1, _2, kBackground ) );
#endif


}// SpecMeasManager

//Moved what use to be SpecMeasManager, out to a startSpectrumManager() to correct modal issues
void  SpecMeasManager::startSpectrumManager()
{
    const bool showToolTipInstantly = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );

    m_spectrumManagerWindow = new AuxWindow( "Spectrum Manager",
                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                     | AuxWindowProperties::TabletModal) );
    m_spectrumManagerWindow->addStyleClass( "SpecMeasManager" );
    
    WContainerWidget *title = m_spectrumManagerWindow->titleBar();
    title->addStyleClass( "SpectraFileManagerHeader" );
    
    m_spectrumManagerWindow->finished().connect( boost::bind( &SpecMeasManager::displayIsBeingHidden, this ) );
    // collapsed().connect( boost::bind( &SpecMeasManager::displayIsBeingHidden, this ) );
    m_spectrumManagerWindow->expanded().connect( boost::bind( &SpecMeasManager::displayIsBeingShown, this ) );
    
    WContainerWidget *uploadDiv = new WContainerWidget( );
    uploadDiv->setStyleClass( "uploadDiv" );
    
    //WText* spec =
    new WText("Load spectrum from: ", uploadDiv);
    //spec->setIcon( "InterSpec_resources/images/plus_min_white.png" );
    if( m_viewer->isSupportFile() )
    {
      Wt::WPushButton* uploadButton = new Wt::WPushButton("File...",uploadDiv);
      uploadButton->clicked().connect(  this, &SpecMeasManager::uploadSpectrum );
      HelpSystem::attachToolTipOn(uploadButton, "Import spectrum from file", showToolTipInstantly, HelpSystem::Bottom );
      uploadButton->setIcon( "InterSpec_resources/images/file_search.png" );
      uploadButton->setMargin(10,Wt::Left);
    } //isSupportFile()
    
#if( USE_DB_TO_STORE_SPECTRA )
    Wt::WPushButton* importButton = new Wt::WPushButton( "Saved Snapshots...", uploadDiv );
    importButton->clicked().connect( boost::bind(  &SpecMeasManager::browseDatabaseSpectrumFiles, this, "", (SpectrumType)0, std::shared_ptr<SpectraFileHeader>()) );
    HelpSystem::attachToolTipOn(importButton, "Imports previously saved spectrum", showToolTipInstantly , HelpSystem::Bottom);
    importButton->setIcon( "InterSpec_resources/images/db_small_white.png" );
    importButton->setMargin(2,Wt::Left);
    
#endif
    
    Wt::WContainerWidget *treeDiv   = createTreeViewDiv();
    Wt::WContainerWidget *buttonBar = createButtonBar();
    createInfoHandler(); // This sets up the m_infoHandler
    
    WContainerWidget * content = m_spectrumManagerWindow->contents();
  
    WGridLayout *layout = m_spectrumManagerWindow->stretcher();
    layout->addWidget(uploadDiv, 0,0);
    layout->addWidget( treeDiv,       1, 0 );
    //  layout->addWidget( m_infoHandler, 1, 0 );
    layout->addWidget( buttonBar,     2, 0 );
    layout->setRowStretch( 1, 10 );
    
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    layout->setContentsMargins( 5, 5, 5, 5 );
    
    dynamic_cast<WWebWidget *>(layout->parent())->setHiddenKeepsGeometry(true);
    
    content->setOverflow(WContainerWidget::OverflowVisible); //necessary for menu to not be covered by footer
    content->setStyleClass("filemanageroverflow");
    
    
    WPushButton *cancel = m_spectrumManagerWindow->addCloseButtonToFooter();
    cancel->clicked().connect( m_spectrumManagerWindow, &AuxWindow::hide );
    
    
    // Make it so it can't be totally deformed
    //  setMinimumSize( 500, 400 );
    m_spectrumManagerWindow->finished().connect( boost::bind( &SpecMeasManager::deleteSpectrumManager, this ) );
    m_spectrumManagerWindow->rejectWhenEscapePressed();
    m_spectrumManagerWindow->resizeWindow( 800, 550 );
    m_spectrumManagerWindow->centerWindow();
    
    selectionChanged();
} //startSpectrumManager()


SpecMeasManager::~SpecMeasManager()
{
  std::lock_guard<std::mutex> lock( *m_destructMutex );
  
  (*m_destructed) = true;
} // SpecMeasManager::~SpecMeasManager()


#if( !ANDROID && !IOS )
FileDragUploadResource *SpecMeasManager::dragNDrop( SpectrumType type )
{
  switch( type )
  {
    case kForeground:
      return m_foregroundDragNDrop;
    case kSecondForeground:
      return m_secondForegroundDragNDrop;
    case kBackground:
      return m_backgroundDragNDrop;
  }//switch( type )

  throw std::runtime_error( "Seriou problem in SpecMeasManager::dragNDrop(..)" );
  return NULL;
}//FileDragUploadResource *dragNDrop( SpectrumType type )

FileDragUploadResource *SpecMeasManager::foregroundDragNDrop()
{
  return m_foregroundDragNDrop;
}

FileDragUploadResource *SpecMeasManager::secondForegroundDragNDrop()
{
  return m_secondForegroundDragNDrop;
}

FileDragUploadResource *SpecMeasManager::backgroundDragNDrop()
{
  return m_backgroundDragNDrop;
}
#endif  //if( !ANDROID && !IOS )


#if( SUPPORT_ZIPPED_SPECTRUM_FILES )

void SpecMeasManager::extractAndOpenFromZip( const std::string &spoolName,
                                             WButtonGroup *group,
                                             WTreeView *table,
                                             AuxWindow *window,
                                             WModelIndex index )
{
  try
  {
    //const string fileInZip = selection->valueText().toUTF8();
    if( !index.isValid() )
    {
      WModelIndexSet selected = table->selectedIndexes();
      if( selected.size() )
        index = *selected.begin();
    }
    
    if( !index.isValid() )
      throw runtime_error( "No file selected" );
    
    const SpectrumType type = SpectrumType( group->checkedId() );
    
    const string fileInZip = Wt::asString(index.data()).toUTF8();
    
    const string tmppath = UtilityFunctions::temp_dir();
    string tmpfile = UtilityFunctions::temp_file_name( "", tmppath );
    
    ifstream zipfilestrm( spoolName.c_str(), ios::in | ios::binary );
    
    ZipArchive::FilenameToZipHeaderMap headers
                                    = ZipArchive::open_zip_file( zipfilestrm );
    
    if( !headers.count(fileInZip) )
      throw runtime_error( "Couldnt find file in zip" );

    size_t nbytewritten = 0;
    
    {
      ofstream tmpfilestrm( tmpfile.c_str(), ios::out | ios::binary );
      nbytewritten = read_file_from_zip( zipfilestrm, headers[fileInZip], tmpfilestrm );
    }
    
    handleFileDrop( fileInZip, tmpfile, type );
    
    UtilityFunctions::remove_file( tmpfile );
  }catch( std::exception & )
  {
    passMessage( "Error extracting file from zip", "", 2 );
  }//try / catch
  
  delete window;
}//SpecMeasManager::extractAndOpenFromZip(...)


bool SpecMeasManager::handleZippedFile( const std::string &name,
                                        const std::string &spoolName,
                                        const int spectrum_type )
{
  try
  {
    WApplication *app = WApplication::instance();
    std::unique_ptr< WApplication::UpdateLock > lock;
    
    if( !app )
      app = dynamic_cast<WApplication *>( m_viewer->parent() );
    if( app )
      lock.reset( new WApplication::UpdateLock( app ) );
    
    ifstream zipfilestrm( spoolName.c_str(), ios::in | ios::binary );
    
    ZipArchive::FilenameToZipHeaderMap headers = ZipArchive::open_zip_file( zipfilestrm );
    
    
    vector<string> filenames;
    vector<uint32_t> uncompresssize;
    for( const ZipArchive::FilenameToZipHeaderMap::value_type &t : headers )
    {
      filenames.push_back( t.first );
      uncompresssize.push_back( t.second->uncompressed_size );
    }
    
    const bool validtype = ((spectrum_type==kForeground)
                             || (spectrum_type==kSecondForeground)
                             || (spectrum_type==kBackground));

    
    string txt;
    if( name.size() )
      txt += "<b>" + name + "</b>";
    else
      txt += "Uploaded file";
    txt += " is a ZIP file.";
    txt += (m_viewer->isPhone() ? "<br />" : "<br /><br />");
    
    txt += "Select which file in it you'd like to open";
    
    WText *t = new WText( txt );
    //WSelectionBox *selection = new WSelectionBox();
    
    RowStretchTreeView *table = new RowStretchTreeView();
    table->setRootIsDecorated( false );
    table->setAlternatingRowColors( true );
    table->setSelectionMode( Wt::SingleSelection );
    table->addStyleClass( "FilesInZipTable" );
    WStandardItemModel *model = new WStandardItemModel( table );
    table->setModel( model );
    model->insertColumns( 0, 2 );
    model->setHeaderData(  0, Horizontal, WString("Filename"), DisplayRole );
    model->setHeaderData(  1, Horizontal, WString("kb"), DisplayRole );
    WItemDelegate *delegate = new WItemDelegate( table );
    delegate->setTextFormat( "%.1f" );
    table->setItemDelegateForColumn( 1, delegate );
    
    
    WContainerWidget *typecb = new WContainerWidget();
    WGridLayout *cblayout = new WGridLayout( typecb );
    cblayout->setContentsMargins( 2, 0, 0, 0 );
    cblayout->setVerticalSpacing( 0 );
    cblayout->setHorizontalSpacing( 0 );
    WButtonGroup *group = new WButtonGroup(typecb);
    Wt::WRadioButton *button = new WRadioButton( "Foreground" );
    cblayout->addWidget( button, 0, 0, AlignCenter );
    group->addButton(button, kForeground);
    button = new WRadioButton("Background", typecb);
    cblayout->addWidget( button, 0, 1, AlignCenter );
    group->addButton(button, kBackground);
    button = new Wt::WRadioButton("Secondary", typecb);
    cblayout->addWidget( button, 0, 2, AlignCenter );
    group->addButton(button, kSecondForeground);
    
    if( validtype )
    {
      group->setCheckedButton( group->button(spectrum_type) );
      typecb->hide();
    }else if( !m_viewer->displayedHistogram(kForeground) )
    {
      group->setCheckedButton( group->button(kForeground) );
      typecb->hide();
    }else
    {
      group->setCheckedButton( group->button(kForeground) );
    }
    
    AuxWindow *window = new AuxWindow( "Uploaded ZIP File Contents",
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::TabletModal) );
    window->stretcher()->addWidget( t, 0, 0 );
    //window->stretcher()->addWidget( selection, 1, 0 );
    window->stretcher()->addWidget( table, 1, 0 );
    window->stretcher()->addWidget( typecb, 2, 0 );
    window->stretcher()->setRowStretch( 1, 1 );
    window->stretcher()->setContentsMargins( 9, 9, 9, 2 );
    
    size_t maxnumchars = 0;
    for( size_t i = 0; i < filenames.size(); ++i )
    {
      const string &n = filenames[i];
      maxnumchars = std::max( maxnumchars, n.size() );
      
      if( n.size() >= 4
          && UtilityFunctions::iequals(n.substr(n.size()-4), ".zip" ) )
        continue;
      if( n.size() < 1 )
        continue;
      //selection->addItem( n );
      
      WStandardItem *item = new WStandardItem();
      item->setText( n );
      model->setItem( i, 0, item );
      item = new WStandardItem();
      
      const double sizekb = uncompresssize[i]/1024.0;
      
      //item->setText( buff );
      item->setData( sizekb, DisplayRole );
      model->setItem( i, 1, item );
    }//for( size_t i = 0; i < filenames.size(); ++i )
    
    //if( selection->count() == 0 )
    if( model->rowCount() == 0 )
    {
      delete window;
      return false;
    //}if( selection->count() == 1 )
    }else if( model->rowCount() == 1 && validtype )
    {
      extractAndOpenFromZip( spoolName, group, table, window, model->index(0,0) );
      return (uncompresssize[0] > 0);

/*
      //const string fileInZip = selection->itemText(0).toUTF8();
      const string fileInZip = filenames[0]; //Wt::asString(model->data(0, 0)).toUTF8();
      const string tmppath = UtilityFunctions::temp_dir();
      string tmpfile = UtilityFunctions::temp_file_name( fileInZip, tmppath );
      
      ofstream tmpfilestrm( tmpfile.c_str(), ios::out | ios::binary );
      size_t nbytes = 0;
      
      try
      {
        nbytes = ZipArchive::read_file_from_zip( zipfilestrm,
                                        headers.begin()->second, tmpfilestrm );
        handleFileDrop( fileInZip, tmpfile.string<string>(), type );
      }catch( std::exception & )
      {
      }
      
      UtilityFunctions::remove_file( tmpfile.string<string>() );
      
      delete window;
      
      return (nbytes > 0);
*/
    }
    
    //Fix the width of column 1, so it is independant of window width
    table->Wt::WTreeView::setColumnWidth( 1, WLength(5,WLength::FontEm) );
    
    //Set column 0 to have a stretchy width so it will expand to full width.
    table->setColumnWidth( 0, 150 /*WLength(maxnumchars,WLength::FontEx)*/ );
    
    if( !m_viewer->isPhone() )
    {
      //If were not on a phone, lets roughly pre-calculate the window size, so
      //  we can size the window so it will center properly when displayed.
      //  This sizing code has not been tested well at all.
      window->setMaximumSize( m_viewer->renderedWidth(), m_viewer->renderedHeight() - 20 );
      
      //Give the window a concrete width/height so window will center properly.
      const int sh = m_viewer->renderedHeight();
      const int sw = m_viewer->renderedWidth();
      const int nrows = model->rowCount();
      const double tableheight = 10 + nrows*table->rowHeight().toPixels()
                                  + table->headerHeight().toPixels();
      const double height = min( tableheight + 180.0, sh-20.0);
      window->setHeight( height );
      
      const double txtw = 7.9*maxnumchars;
      if( (txtw+60) < 0.85*sw )
        window->setWidth( max(txtw,min(300.0,0.95*(sw-60))) + 60.0 );
      else
        window->setWidth( 0.85*sw );
    }//if( !m_viewer->isPhone() )
    
    window->rejectWhenEscapePressed();
    Wt::WPushButton *closeButton = window->addCloseButtonToFooter("Cancel",true);
    closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
    WPushButton *openButton = new WPushButton( "Display" );
    openButton->disable();
    //selection->activated().connect( openButton, &WPushButton::enable );
    table->clicked().connect( openButton, &WPushButton::enable );
    table->doubleClicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this, spoolName, group, table, window, _1 ) );
    window->footer()->addWidget( openButton );
    
    //openButton->clicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this, spoolName, type, selection, window ) );
    openButton->clicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this, spoolName, group, table, window, WModelIndex() ) );
    
    window->setResizable( true );
    window->centerWindow();
    window->disableCollapse();
    window->show();
    
    if( app )
      app->triggerUpdate();
  }catch( std::exception & )
  {
    return false;
  }
  
  return true;
}//void handleZippedFile(...)
#endif //#if( SUPPORT_ZIPPED_SPECTRUM_FILES )



template<size_t N, size_t M>
bool check_magic_number( const uint8_t (&magic_num)[N], const uint8_t (&data)[M] )
{
  return (0 == memcmp( magic_num, data, (N < M ? N : M)));
}


bool SpecMeasManager::handleNonSpectrumFile( const std::string &displayName,
                                             const std::string &fileLocation )
{
  ifstream infile( fileLocation.c_str(), ios_base::binary|ios_base::in );
  if( !infile )
    return false;
 
  //get the filesize
  infile.seekg(0, ios::end);
  const size_t filesize = infile.tellg();
  infile.seekg(0);
  
  if( filesize <= 1024 )
  {
    passMessage( "File was probably an empty file - not a spectrum file.", "", 2 );
    return true;
  }
  
  AuxWindow *w = new AuxWindow( "Not a spectrum file",
                (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                 | AuxWindowProperties::TabletModal
                 | AuxWindowProperties::SetCloseable
                 | AuxWindowProperties::DisableCollapse) );
  w->centerWindow();
  w->rejectWhenEscapePressed( true );
  WPushButton *b = w->addCloseButtonToFooter();
  b->clicked().connect( w, &AuxWindow::hide );
  w->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, w ) );
  if( !m_viewer->isMobile() && m_viewer->renderedWidth() > 400 && m_viewer->renderedHeight() > 250 )
    w->resize( 400, 250 );
  
  uint8_t data[1024] = { 0x0 };
  
  infile.read( (char *)data, boost::size(data) );
  infile.seekg(0);
  
  // const string filename = UtilityFunctions::filename(displayName);
  
  //Check if ICD2 file
  if( std::find( boost::begin(data), boost::end(data), uint8_t(60)) != boost::end(data) )
  {
    string datastr;
    datastr.resize( boost::size(data) + 1, '\0' );
    memcpy( &(datastr[0]), (const char *)&data[0], boost::size(data) );
    
    if( UtilityFunctions::icontains( datastr, "n42ns:")
       || UtilityFunctions::icontains( datastr, "AnalysisResults" )
       || UtilityFunctions::icontains( datastr, "AlarmInformation" )
       || UtilityFunctions::icontains( datastr, "DNDOARSchema" )
       || UtilityFunctions::icontains( datastr, "DNDOEWSchema" )
       || UtilityFunctions::icontains( datastr, "DNDOARSchema" ) )
    {
      WText *t = new WText( "This looks to be an N42 ICD2 file that contains analysis results rather than raw spectra.<br />"
                            "If you believe this to be a legitimate spectrum file, please email it to <a href=\"mailto:interspec@sandia.gov\" target=\"_blank\">interspec@sandia.gov</a> to support this file type." );
      w->stretcher()->addWidget( t, 0, 0, AlignCenter | AlignMiddle );
      t->setTextAlignment( Wt::AlignCenter );
      w->show();
      
      return true;
    }
  }//if( might be ICD2 )
  
  
  
  const uint8_t zip_mn[]    = { 0x50, 0x4B, 0x03, 0x04 };
  const uint8_t rar4_mn[]   = { 0x52, 0x61, 0x72, 0x21, 0x1A, 0x07, 0x00 };
  const uint8_t rar5_mn[]   = { 0x52, 0x61, 0x72, 0x21, 0x1A, 0x07, 0x01, 0x00 };
  const uint8_t tar_mn[]    = { 0x75, 0x73, 0x74, 0x61, 0x72 };
  const uint8_t zip7_mn[]   = { 0x37, 0x7A, 0xBC, 0xAF, 0x27, 0x1C };
  const uint8_t gz_mn[]     = { 0x1F, 0x8B };

  
  //\sa Wt::Utils::guessImageMimeTypeData(...)
  const uint8_t gifold_mn[] = { 0x47, 0x49, 0x46, 0x38, 0x37, 0x61 };
  const uint8_t gifnew_mn[] = { 0x47, 0x49, 0x46, 0x38, 0x39, 0x61 };
  const uint8_t tiffle_mn[] = { 0x49, 0x49, 0x2A, 0x00 };
  const uint8_t tiffbe_mn[] = { 0x4D, 0x4D, 0x00, 0x2A };
  const uint8_t jpeg1_mn[]  = { 0xFF, 0xD8, 0xFF, 0xDB };
  const uint8_t jpeg2_mn[]  = { 0xFF, 0xD8, 0xFF, 0xE0 };
  const uint8_t jpeg3_mn[]  = { 0x49, 0x46, 0x00, 0x01 };
  const uint8_t jpeg4_mn[]  = { 0xFF, 0xD8, 0xFF, 0xE1 };
  const uint8_t jpeg5_mn[]  = { 0x69, 0x66, 0x00, 0x00 };
  const uint8_t png_mn[]    = { 0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A };
  const uint8_t bmp_mn[]    = { 0x42, 0x4D };
  
  const uint8_t pdf_mn[]    = { 0x25, 0x50, 0x44, 0x46 };
  const uint8_t ps_mn[]     = { 0x25, 0x21, 0x50, 0x53 };
  
  const bool iszip = check_magic_number( data, zip_mn );
  
  const bool israr = (check_magic_number( data, rar4_mn ) || check_magic_number( data, rar5_mn ));
  const bool istar = check_magic_number( data, tar_mn );
  const bool iszip7 = check_magic_number( data, zip7_mn );
  const bool isgz = check_magic_number( data, gz_mn );
  
  const bool ispdf = check_magic_number( data, pdf_mn );
  const bool isps = check_magic_number( data, ps_mn );
  const bool istif = (check_magic_number( data, tiffle_mn ) || check_magic_number( data, tiffbe_mn ));
  
  const bool isgif = (check_magic_number( data, gifold_mn ) || check_magic_number( data, gifnew_mn ));
  const bool isjpg = (check_magic_number( data, jpeg1_mn )
                      || check_magic_number( data, jpeg2_mn )
                      || check_magic_number( data, jpeg3_mn )
                      || check_magic_number( data, jpeg4_mn )
                      || check_magic_number( data, jpeg5_mn ));
  const bool ispng = check_magic_number( data, png_mn );
  const bool isbmp = check_magic_number( data, bmp_mn );
  
  if( iszip )
  {
    const char *msg = NULL;
#if( SUPPORT_ZIPPED_SPECTRUM_FILES )
  //zip (but can be xlsx, pptx, docx, odp, jar, apk) 50 4B 03 04
    msg = "This file appears to be an invalid ZIP file (or xlsx, pptx, <br />"
          "docx, odp, jar, apk, etc), sorry I cant open it :(";
#else
    msg = "This build of InterSpec does not support ZIP files.<br />"
          "Please contact "
          "<a href=\"mailto:wcjohns@sandia.gov\" target=\"_blank\">wcjohns@sandia.gov</a> "
          "for a version that supports ZIP files.";
#endif
    WText *t = new WText( msg );
    w->stretcher()->addWidget( t, 0, 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    w->show();
    
    return true;
  }//if( iszip )
  
  if( israr || istar || iszip7 || isgz )
  {
    const char *msg = "This file looks to be an archive file, not supported by InterSpec yet."
                      "Please contact "
                      "<a href=\"mailto:wcjohns@sandia.gov\" target=\"_blank\">wcjohns@sandia.gov</a> "
                      "if you would like support for this archive type added.";
    
    WText *t = new WText( msg );
    w->stretcher()->addWidget( t, 0, 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    w->show();
    
    return true;
  }//if( israr || istar || iszip7 || isgz )
  
  
  if( ispdf | isps | istif )
  {
    const char *msg = "This file looks to be an document file, and not supported by InterSpec.";
    WText *t = new WText( msg );
    w->stretcher()->addWidget( t, 0, 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    w->show();
    
    return true;
  }//if( ispdf | isps | istif )
  
  
  if( isgif || isjpg || ispng || isbmp )
  {
    const char *msg = "This file looks to be an image, and not a spectrum file.";
    WText *t = new WText( msg );
    t->setTextAlignment( Wt::AlignCenter );
    w->stretcher()->addWidget( t, 0, 0, AlignCenter | AlignMiddle );
   
    const size_t max_disp_size = 5*1024*1024;
    
    if( filesize < max_disp_size )
    {
      const char *mimetype = "";
      if( isgif ) mimetype = "image/gif";
      else if( isjpg ) mimetype = "image/jpeg";
      else if( ispng ) mimetype = "image/png";
      else if( isbmp ) mimetype = "image/bmp";
    
      std::unique_ptr<WImage> image( new WImage() );
      WMemoryResource *resource = new WMemoryResource( mimetype, image.get() );
      vector<uint8_t> totaldata( filesize );
      const bool success = infile.read( (char *)&(totaldata[0]), filesize ).good();
      
      if( !m_viewer->isMobile() )
      {
        w->resize( WLength::Auto, WLength::Auto );
        w->setMaximumSize( 0.6*m_viewer->renderedWidth(), 0.8*m_viewer->renderedHeight() );
      }
      
      if( success )
      {
        resource->setData( totaldata );
        image->setImageLink( WLink(resource) );
        w->stretcher()->addWidget( image.release(), 1, 0, AlignCenter | AlignMiddle );
      }else
      {
        WText *errort = new WText( "Couldnt read uplaoded file." );
        errort->setTextAlignment( Wt::AlignCenter );
        w->stretcher()->addWidget( errort, 1, 0 );
      }
    }else
    {
      WText *errort = new WText( "Uploaded file was to large to try and display." );
      errort->setTextAlignment( Wt::AlignCenter );
      w->stretcher()->addWidget( errort, 1, 0, AlignCenter | AlignMiddle );
    }//if( filesize < max_disp_size ) / else
    
    w->show();
    w->resizeToFitOnScreen();
    w->centerWindowHeavyHanded();
    
    return true;
  }//if( isgif || isjpg || ispng || isbmp )

  delete w;
  
  return false;
}//void handleNonParsableFile(...)


void SpecMeasManager::handleFileDrop( const std::string &name,
                                             const std::string &spoolName,
                                             SpectrumType type )
{
#if( SUPPORT_ZIPPED_SPECTRUM_FILES )
  if( name.length() > 4
     && UtilityFunctions::iequals( name.substr(name.length()-4), ".zip")
     && handleZippedFile( name, spoolName, type ) )
    return;
#endif
  
  std::shared_ptr< SpectraFileHeader > header;
  std::shared_ptr< SpecMeas >   measurement;

  WApplication *app = WApplication::instance();
  std::unique_ptr< WApplication::UpdateLock > lock;
  
  //Should maybe consider posting to the WServer if we arent in the event loop.
  //WServer::instance()->post( appid, boost::bind( &SpecMeasManager::handleFileDrop, this, name, spoolName, type ) );
  
  if( !app )
    app = dynamic_cast<WApplication *>( m_viewer->parent() );
  
  // If the file was uploaded using a POST request to a WResource, then we have to notify
  // WApplication that we have changed things.
  if( app )
    lock.reset( new WApplication::UpdateLock( app ) );
  
  try
  {
    const int modelRow = setFile( name, spoolName, header, measurement );

    displayFile( modelRow, measurement, type, true, true, true );
    
    //It is the responsibility of the caller to clean up the file.
  }catch( exception &e )
  {
    if( !handleNonSpectrumFile( name, spoolName ) )
    {
      displayInvalidFileMsg(name,e.what());
    }
  }
  
  if( app )
    app->triggerUpdate();
}//handleFileDrop(...)

void SpecMeasManager::displayInvalidFileMsg( std::string filename, std::string errormsg )
{
  //make sure we dont display the whole path
  string lastpart = UtilityFunctions::filename(filename);
  if( lastpart.empty() )
    lastpart = filename;
  if( lastpart.size() > 12 )
    lastpart = lastpart.substr(0,9) + "...";
  
  if( errormsg.empty() )
    errormsg = "Unspecified";
  
  lastpart = Wt::Utils::htmlEncode(lastpart);
  errormsg = Wt::Utils::htmlEncode(errormsg);
  
  stringstream msg;
  msg << "<p>Sorry, I couldnt parse the file " << lastpart << ":</p>"
  << "<p>Error: <em>" << errormsg << "</em></p>"
  << "<p>If you think this is a valid spectrum file, please send it to "
  << "<a href=\"mailto:interspec@sandia.gov\" target=\"_blank\">interspec@sandia.gov</a>, and"
  << " we'll try to fix this issue.</p>";
  
  passMessage( msg.str(), "", 2 );
}//void displayInvalidFileMsg( std::string filename, std::string errormsg )

std::set<int> SpecMeasManager::selectedSampleNumbers() const
{
  std::set<int> sample_nums;
  const WModelIndexSet selected = m_treeView->selectedIndexes();

  std::shared_ptr<SpectraFileHeader> header;

  for( const WModelIndex &index : selected )
  {
    const int row = index.row();
    if( row < 0 )
      continue;

    const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
    switch( indexLevel )
    {
      case SpectraFileModel::FileHeaderLevel:
        header = m_fileModel->fileHeader( index.row() );
        if( header )
        {
          for( const SpectraHeader &spectra : header->m_samples )
            sample_nums.insert( spectra.sample_number );
        }//if( header )
      break;

      case SpectraFileModel::SampleLevel:
        header = m_fileModel->fileHeader( index.parent().row() );
        if( header && (row<static_cast<int>(header->m_samples.size()) ) )
          sample_nums.insert( header->m_samples[row].sample_number );
      break;

      case SpectraFileModel::InvalidLevel:
      break;
    } // switch( level(index) )
  } // for( const WModelIndex &index : selected )

  return sample_nums;
} // std::set<int> SpecMeasManager::selectedSampleNumbers()

/**
Returns a vector of all the selected files.  If you just need one, use selectedFile()
 
 Also, only returns FileHeaderLevel files.  No longer returns SampleLevel or InvalidLevels.  This will thus only return the main file, and not say, sample files like in passthrough.
**/
vector<std::shared_ptr<SpectraFileHeader> > SpecMeasManager::getSelectedFiles() const
{
    const WModelIndexSet selected = m_treeView->selectedIndexes();
    vector<std::shared_ptr<SpectraFileHeader> >  ret;
    
    
    if( selected.empty() )
    {
        return ret;
    }
    
    for (WModelIndexSet::iterator iter = selected.begin(); iter!=selected.end(); iter++)
    {
        std::shared_ptr<SpectraFileHeader> header;

        const WModelIndex &index = *(iter);
        
        const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
        switch( indexLevel )
        {
            case SpectraFileModel::FileHeaderLevel:
                header = m_fileModel->fileHeader( index.row() );
                break;
                
            case SpectraFileModel::SampleLevel:
                //Do not return SampleLevel
                continue;
                //header = m_fileModel->fileHeader( index.parent().row() );
                break;
                
            case SpectraFileModel::InvalidLevel:
                continue;
                break;
        } // switch( level(index) )
        ret.push_back(header);
    }
    return ret;
} // std::shared_ptr<SpectraFileHeader> SpecMeasManager::selectedFile() const

//Only returns the first selected file.  use getSelectedFiles() for more than 1
std::shared_ptr<SpectraFileHeader> SpecMeasManager::selectedFile() const
{
  const WModelIndexSet selected = m_treeView->selectedIndexes();

  std::shared_ptr<SpectraFileHeader> header;

  if( selected.empty() )
    return header;

  const WModelIndex &index = *(selected.begin());

  const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
  switch( indexLevel )
  {
    case SpectraFileModel::FileHeaderLevel:
      header = m_fileModel->fileHeader( index.row() );
    break;

    case SpectraFileModel::SampleLevel:
      header = m_fileModel->fileHeader( index.parent().row() );
    break;

    case SpectraFileModel::InvalidLevel:
    break;
  } // switch( level(index) )

  // XXX a logic check that should be removed!
  for( const WModelIndex &index : selected )
  {
    std::shared_ptr<SpectraFileHeader> checkheader;
    if( !index.parent().isValid() )
    {
      checkheader = m_fileModel->fileHeader( index.row() );
    }else checkheader = m_fileModel->fileHeader( index.parent().row() );
    if( checkheader != header )
    {
      passMessage( "Currently cant load samples from multiple files, sorry, will is lazy.", "", 3 );
      return std::shared_ptr<SpectraFileHeader>();
    }
    
  }//for( const WModelIndex &index : selected )

  return header;
} // std::shared_ptr<SpectraFileHeader> SpecMeasManager::selectedFile() const


void SpecMeasManager::unDisplay( SpectrumType type )
{
  std::set<int> displaySampleNums;
  std::shared_ptr<SpecMeas> meas;
  m_viewer->setSpectrum( meas, displaySampleNums, type, false );
  selectionChanged(); // update buttons
} // void SpecMeasManager::unDisplay( SpectrumType type );


// The std::shared_ptr<SpecMeas> dummy_ptr keeps the MeasurementInfo
// object in memory when SpectraFileHeader isnt caching the spectrum, so
// its weak_ptr<> can be used in the call to header->parseFile();
void SpecMeasManager::loadSelected( const SpectrumType type,
                                    std::shared_ptr<SpecMeas> dummy_ptr,
                                    const bool doPreviousEnergyRangeCheck )
{
  dummy_ptr = dummy_ptr; //keep compiler from complaining or optimizing dummy_ptr away
  loadSelected( type, doPreviousEnergyRangeCheck );
} // void SpecMeasManager::loadSelected(...)


void SpecMeasManager::loadSelected( const SpectrumType type,
                                    const bool doPreviousEnergyRangeCheck )
{
  std::shared_ptr<SpectraFileHeader> header = selectedFile();

  if( !header )
    return;

  std::shared_ptr<SpecMeas> meas = header->parseFile();

  if( !meas )
  {
    cerr << SRC_LOCATION << "\n\tSerious error, couldnt parse file" << endl;
    passMessage( "Couldn't parse file!", SRC_LOCATION, 3 );
    return;
  }//if( !meas )

  const set<int> displaySampleNums = selectedSampleNumbers();

  m_viewer->setSpectrum( meas, displaySampleNums, type, doPreviousEnergyRangeCheck );

  addToTempSpectrumInfoCache( meas );

  selectionChanged(); // update buttons
} // void SpecMeasManager::loadSelected(...)


void SpecMeasManager::startQuickUpload()
{
  new FileUploadDialog( m_viewer, this );
}//void startQuickUpload( SpectrumType type )


void SpecMeasManager::finishQuickUpload( Wt::WFileUpload *upload,
                                         const SpectrumType type )
{
  // TODO: The warning messages, and error conditions detected should be greatly improved

  std::shared_ptr<SpecMeas> measement_ptr;
  const int row = dataUploaded( upload, measement_ptr );
  
  if( row < 0 )
  {
    cerr << SRC_LOCATION << "\n\tError uploading file "
         << upload->clientFileName() << endl;
    return;
  } // if( row < 0 )

  displayFile( row, measement_ptr, type, true, true, true );
}//void finishQuickUpload(...)


void SpecMeasManager::selectEnergyBinning( const string binning,
                                  std::shared_ptr<SpectraFileHeader> header,
                                  std::shared_ptr<SpecMeas> meas,
                                  const SpectrumType type,
                                  const bool checkIfPreviouslyOpened,
                                  const bool doPreviousEnergyRangeCheck )
{
  if( binning != "Keep All" )
    meas->keep_energy_cal_variant( binning );
  
  WModelIndex index = m_fileModel->index( header );
  
  if( !index.isValid() )
  {
    passMessage( "Aborting loading of file after selecting energy binning - "
                 "the file is no longer available in memory. please report "
                 "this bug to wcjohns@sandia.gov", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck, false );
}//


bool SpecMeasManager::checkForAndPromptUserForDisplayOptions( std::shared_ptr<SpectraFileHeader> header,
                                            std::shared_ptr<SpecMeas> meas,
                                            const SpectrumType type,
                                            const bool checkIfPreviouslyOpened,
                                            const bool doPreviousEnergyRangeCheck )
{
  if( !header || !meas )
    throw runtime_error( "SpecMeasManager::checkForAndPromptUserForDisplayOptions(): Invalid input" );
  
  const set<string> cals = meas->energy_cal_variants();
  
  if( cals.size() < 2 )
    return false;
  
  
  AuxWindow *dialog = new AuxWindow( "Select Binning",
                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                       | AuxWindowProperties::TabletModal
                       | AuxWindowProperties::DisableCollapse) );
  dialog->rejectWhenEscapePressed( false );
  dialog->setClosable( false );
  dialog->centerWindow();

  const int ww = m_viewer->renderedWidth();
  const int wh = m_viewer->renderedHeight();
  
  if( !m_viewer->isMobile() && ww > 420 && wh > 125 )
  {
    dialog->setWidth( 420 );
    dialog->Wt::WCompositeWidget::setMinimumSize( 420, 125 );
  }
  
  int ncolwide = static_cast<int>( cals.size() + 1 );
  if( ncolwide > 4 )
    ncolwide = 1;
  
  WGridLayout *layout = dialog->stretcher();
  
  WText *msg = new WText( "<div>Multiple energy binnings were found in the spectrum file.</div>"
                          "<div>Please select which one you would like</div>" , XHTMLText );
  layout->addWidget( msg, 0, 0 );
  
  //The buttons dont space properly if they are in the same layout as the text
  //  for some reason.
  WGridLayout *buttonlayout = new WGridLayout();
  layout->addLayout( buttonlayout, 1, 0 );
  
  WPushButton *button = new WPushButton( "Keep All" );
  buttonlayout->addWidget( button, 0, 0 );
  
  button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, dialog ) );
  button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this, string("Keep All"), header, meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
  
  int calnum = 1;
  for( set<string>::const_iterator iter = cals.begin(); iter != cals.end(); ++iter, ++calnum )
  {
    const int row = (calnum / ncolwide);
    const int col = (calnum % ncolwide);
    
    //Make sure the calbration ID isnt too long
    string label = *iter;
    //Not sure if utf8 utilitiities are properly tested, so will skip for now...
    //UtilityFunctions::utf8_limit_str_size( label, 15 );
    if( label.size() > 15 )
      label.resize( 15 );
    
    button = new WPushButton( label );
    buttonlayout->addWidget( button, row, col );

    button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, dialog ) );
    //Note, using *iter instead of label below.
    button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this, *iter, header, meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
  }
  
  //layout->setRowStretch( 0, 1 );
  
  return true;
}//checkForAndPromptUserForDisplayOptions(...)



void SpecMeasManager::displayFile( int row,
                                   std::shared_ptr<SpecMeas> measement_ptr,
                                   const SpectrumType type,
                                   bool checkIfPreviouslyOpened,
                                   const bool doPreviousEnergyRangeCheck,
                                   const bool checkIfAppropriateForViewing )
{
  std::shared_ptr<SpecMeas> old_meas = m_viewer->measurment( type );
  std::shared_ptr<SpecMeas>  old_back;
  if( type == kForeground )
    old_back = m_viewer->measurment( kBackground );

#if( USE_DB_TO_STORE_SPECTRA )
  const bool storeInDb
      = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
#endif
  
  
  if( row < 0 && !measement_ptr )
  {
    m_viewer->setSpectrum( measement_ptr, std::set<int>(), type, true );

    if( old_meas.use_count() == 1 )
      serializeToTempFile( old_meas );

#if( USE_DB_TO_STORE_SPECTRA )
    if( storeInDb )
      saveToDatabase( old_meas );
#endif
    
    return;
  }//if( we specifically wanted to unload the file )

  std::shared_ptr<SpectraFileHeader> header;
  header = m_fileModel->fileHeader( row );

  if( !header
      || (header->measurementIfInMemory() != measement_ptr) )
  {
    const char *msg = "SpecMeasManager::displayFile(...): you must "
                      "pass in the SpectraFileModel row cooresponding to the "
                      "MeasurmentInfo object you pass in";
    cerr << msg << endl;
    throw std::runtime_error( msg );
  }//if( test_meas_ptr != measement_ptr )

  WModelIndex index = m_fileModel->index( row, 0, WModelIndex() );
  
  if( checkIfPreviouslyOpened && wApp )
  {
    //Lets check to see if the user has loaded this same spectrum in this same
    //  session, but hasnt made any changes - if so we'll switch to this
    //  SpectraFileHeader, and essentially delete the passed in
    //  measement_ptr, since if checkIfPreviouslyOpened is true, then its
    //  the case the file was just loaded, and hasnt been used previously.
    for( int testrow = 0; testrow < m_fileModel->rowCount(); ++testrow )
    {
      std::shared_ptr<SpectraFileHeader> oldheader
                                          = m_fileModel->fileHeader( testrow );
      if( oldheader != header
          && oldheader->m_uuid == header->m_uuid
          && oldheader->m_fileDbEntry
          && !oldheader->m_fileDbEntry->userHasModified
          && !oldheader->m_fileDbEntry->isWriteProtected() )
      {
        cerr << "We found a copy of this file in memmorry that hasnt been "
             << "modified, switching to that" << endl;
        m_fileModel->removeRows( row, 1 );
        row = testrow;
        index = m_fileModel->index( row, 0, WModelIndex() );
        measement_ptr = header->parseFile();
        checkIfPreviouslyOpened = false;
        break;
      }
    }//for( int row = 0; row < m_model->rowCount(); ++row )
  }//if( header && checkIfPreviouslyOpened )
  
  
  if( checkIfAppropriateForViewing && wApp )
  {
    if( checkForAndPromptUserForDisplayOptions( header, measement_ptr,
                                type, checkIfPreviouslyOpened,
                               doPreviousEnergyRangeCheck ) )
      return;
  }//if( checkIfAppropriateForViewing )
  
  
  WModelIndexSet selected;
  selected.insert( index );
  const int nrow = m_fileModel->rowCount( index );
  for( int i = 0; i < nrow; ++i )
    selected.insert( index.child(i,0) );

  const size_t ncandiadate_samples = selected.size();

  //backgroundIndexs will only get filled if the a forground or secondforground
  // is desired
  std::set<int> background_sample_numbers;

  WarningWidget::WarningMsgLevel warningmsgLevel = WarningWidget::WarningMsgInfo;
  stringstream warningmsg;

  switch( type )
  {
    case kForeground:
    case kSecondForeground:
    {
      const int nsamples = header->numSamples();
      const bool passthrough = header->passthrough();
      if( !passthrough && (nsamples > 1) )
      {
        int numbackground = 0, numForeground = 0, numIntrinsic = 0;
        int childrow = -1, backrow = -1, intrinsicrow = -1;
        for( size_t i = 0; i < header->m_samples.size(); ++i )
        {
          const SpectraHeader &spechead = header->m_samples[i];
          switch( spechead.spectra_type )
          {
            case Measurement::IntrinsicActivity:
              //lets try to not show intrinsic activity by default
              ++numIntrinsic;
              intrinsicrow = static_cast<int>( i );
            break;
              
            case Measurement::Foreground:
              ++numForeground;
              childrow = static_cast<int>( i );
              //i = header->m_samples.size();  //first foreground, terminate loop
            break;
              
            case Measurement::Background:
              ++numbackground;
              backrow = static_cast<int>( i );
              if( numForeground )
                i = header->m_samples.size();  //We have foreground and background, terminate loop
            break;
            
            case Measurement::Calibration:
              //do nothing
            break;
            
            case Measurement::UnknownSourceType:
              if( childrow < 0 )
                childrow = static_cast<int>( i );
            break;
          }//switch( header->m_samples[i].spectra_type )
        }//for( size_t i = 0; i < header->m_samples.size(); ++i )
        
        if( (!numIntrinsic && !numForeground && !numbackground) || (childrow < 0) )
          childrow = 0;
        
        if( numIntrinsic )
        {
          warningmsg << "The uploaded file contained an intrinsic activity"
                        " spectrum, but currently showing the other spectrum"
                        " in the file.";
        }else
        {
          warningmsg << "The uploaded file contained " << nsamples
                     << " samples<br>";
          
          if( numForeground )
            warningmsg << "The first foreground sample is being shown.<br>";
          else if( childrow == 0 )
            warningmsg << "The first sample is being shown.<br>";
          else
            warningmsg << "Sample " << header->m_samples[childrow].sample_number << " is being shown.<br>";
          
          if( m_viewer->toolTabsVisible() )
            warningmsg << "Use the <b>Spectrum Files</b> tab, or the "
                          "<b>File Manager</b> to select other records";
          else
            warningmsg << "Use the <b>File Manager</b> to select others.";
        }//if( decide how to customize the info message ) / else

        //If we have an unambiguos background, and arent currently displaying a background
        //  from this detector, lets load the unambiguos background
        if( type==kForeground
           && numbackground == 1 && childrow != backrow
           && (childrow>=0) && (childrow<static_cast<int>(header->m_samples.size()))
           && (backrow>=0) && (backrow<static_cast<int>(header->m_samples.size())) //These two conditions should always be true if first condition is true
           && ( !old_back || !measement_ptr
                || (measement_ptr->num_gamma_channels()!=old_back->num_gamma_channels()) || (measement_ptr->instrument_id()!=old_back->instrument_id()) )
           )
          background_sample_numbers.insert( header->m_samples[backrow].sample_number );
        
        selected.clear();
        selected.insert( index.child(childrow,0) );
      } // if( !passthrough && (nsamples > 1) )

      if( passthrough )
      {
        int nspectra_header = static_cast<int>( header->m_samples.size() );

        // A temporary check...
        if( nspectra_header != nsamples )
          cerr << "SpecMeasManager::finishQuickUpload..): nspectra_header != nsamples ("
               << nspectra_header << " != " << nsamples << ")" << endl;
        nspectra_header = min( nspectra_header, nsamples );


        int ncalibration = 0;

        // Remove background and calibration
        for( int sample = 0; sample < nspectra_header; ++sample )
        {
          const SpectraHeader &spectra = header->m_samples[sample];
          const bool back = (spectra.spectra_type == Measurement::Background);
          const bool calib = (spectra.spectra_type == Measurement::Calibration);

          ncalibration += calib;

          if( back && (type!=kSecondForeground) )
            background_sample_numbers.insert( spectra.sample_number );

          if( back || calib )
            selected.erase( index.child(sample,0) );
        }//for( int sample = 0; sample < nsamples; ++sample )

        if( ncalibration > 0 )
        {
          selected.erase( index );

          warningmsgLevel = max( WarningWidget::WarningMsgInfo, warningmsgLevel );
          warningmsg << "The uploaded file contained " << ncalibration
                       << " calibration samples<br>";


          if( ncalibration == 1 ) warningmsg << "It";
          else                    warningmsg << "They";
          warningmsg  << " will not be displayed.<br>"
                      << "Use the <b>File Manager</b> to change this.";
        } // if( ncalibration > 0 )
      }// if( passthrough )

      break;
    } // case kForeground: case kSecondForeground:

    case kBackground:
    {
      int nspectra_header = static_cast<int>( header->m_samples.size() );

      if( nspectra_header > 1 )
      {
        bool foundBackground = false;
        for( int sample = 0; sample < nspectra_header; ++sample )
        {
          const SpectraHeader &spectra = header->m_samples[sample];
          foundBackground = (spectra.spectra_type == Measurement::Background);
          if( foundBackground )
          {
            warningmsgLevel = max( WarningWidget::WarningMsgLow, warningmsgLevel );
            warningmsg << "File has mutliple spectra, so only the background"
                          " spectrum was loaded.";
            selected.clear();
            selected.insert( index.child(sample,0) );
            break;
          }//if( foundBackground )
        } // for( int sample = 0; sample < nsamples; ++sample )

        if( !foundBackground )
        {
          warningmsgLevel = max(WarningWidget::WarningMsgHigh, warningmsgLevel );
          warningmsg << "The uploaded file contained " << nspectra_header << " samples so<br>"
                     << "am displaying the first one.<br>"
                     << "Use the <b>File Manager</b> to change this.";
          selected.clear();
          selected.insert( index.child(0,0) );
        }//if( !foundBackground )
      } // if( !passthrough && (nsamples > 1) )
      break;
    } // case kBackground:
  } // switch( type )

  //Check if we removed any check or background samples, if so do not include
  //  the file "index" in the selected samples
  if( ncandiadate_samples != selected.size() )
    selected.erase( index );
  
  m_treeView->setSelectedIndexes( selected );

  if( warningmsg.str().size() )
  {
    passMessage( warningmsg.str(), "SpecMeasManager", warningmsgLevel );
  }
  
  
  loadSelected( type, doPreviousEnergyRangeCheck );

//if 'old_meas' is the last reference to the SpecMeas object, lets go ahead
//  and try to save it to disk for later access
//    typedef std::deque< std::shared_ptr<const SpecMeas> > queue_type;
//    const queue_type::iterator pos = std::find( m_tempSpectrumInfoCache.begin(),
//                                               m_tempSpectrumInfoCache.end(),
//                                               old_meas );
//    if( pos == m_tempSpectrumInfoCache.end() )
//      serializeToTempFile( old_meas );
  
  if( old_meas.use_count() == 1 )
    serializeToTempFile( old_meas );

#if( USE_DB_TO_STORE_SPECTRA )
  if( storeInDb )
    saveToDatabase( old_meas );
#endif

  if( background_sample_numbers.size() )
  {
    m_viewer->setSpectrum( measement_ptr,
                                   background_sample_numbers,
                                   kBackground, doPreviousEnergyRangeCheck );
  }//if( backgroundIndexs.size() )
  
#if( USE_DB_TO_STORE_SPECTRA )
  WApplication *app = wApp;
  if( checkIfPreviouslyOpened && !!header && app )
  {
    boost::function<void(void)> worker = app->bind(
                      boost::bind( &SpecMeasManager::checkIfPreviouslyOpened,
                                   this, app->sessionId(), header, type, m_destructMutex, m_destructed ) );
//    WServer::instance()->post( app->sessionId(), worker );
    WServer::instance()->ioService().post( worker );
  }//if( checkIfPreviouslyOpened )
#endif //#if( USE_DB_TO_STORE_SPECTRA )
}//void displayFile(...)


void SpecMeasManager::setDisplayedToSelected()
{
  throw runtime_error( "void SpecMeasManager::setDisplayedToSelected() Not implemented" );
}//void SpecMeasManager::setDisplayedToSelected()


void SpecMeasManager::displayQuickSaveAsDialog()
{
  AuxWindow *dialog = new AuxWindow( "Save As...",
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
               | AuxWindowProperties::TabletModal
               | AuxWindowProperties::DisableCollapse) );
  
  dialog->centerWindow();
  dialog->addStyleClass( "QuickSaveAsDialog" );

  std::shared_ptr<const SpecMeas> data, second, background, initial;
  data       = m_viewer->measurment(kForeground);
  second     = m_viewer->measurment(kSecondForeground);
  background = m_viewer->measurment(kBackground);

  initial = data;
  if( !initial )
    initial = second;
  if( !initial )
    initial = background;

  if( !initial )
  {
    const string msgstr = "There are no spectrums loaded into the viewer<br>"
                          "Please use the <b>File Manager</b> for more options.";

    WText *msg = new WText( msgstr, XHTMLUnsafeText, dialog->contents() );
    msg->setInline( false );
    WPushButton *ok = new WPushButton( "Ok", dialog->contents() );
    ok->clicked().connect( dialog, &AuxWindow::hide );
    dialog->show();
    return;
  } // if( !initial )


  WModelIndex fileIndex = m_fileModel->index( initial );
  std::shared_ptr<const SpectraFileHeader> fileHeader
                                  = m_fileModel->fileHeader( fileIndex.row() );
  
  const set<int> detnums = m_viewer->displayedDetectorNumbers();
  set<int> samplenums;
  if( data == initial )
    samplenums = m_viewer->displayedSamples( kForeground );
  
  for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
        i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
  {
    if( !m_specificResources[i] )
      m_specificResources[i] = new SpecificSpectrumResource( i, this );
    
    m_specificResources[i]->setSpectrum( initial, samplenums, detnums );
  }//loop over save as types
  
  const string msgstr = "Please select the file format:";

  WText *msg = new WText( msgstr, XHTMLUnsafeText, dialog->contents() );
  msg->setInline( false );

  int nspecs = 0;
  if( data ) ++nspecs;
  if( second ) ++nspecs;
  if( background ) ++nspecs;

  if( nspecs > 1 )
  {
    WGroupBox *buttons = new WGroupBox( "Which file to save:",
                                        dialog->contents() );
    WButtonGroup *group = new WButtonGroup( dialog->contents() );
    if( data )
    {
      WRadioButton *button = new WRadioButton( "Forground", buttons );
      button->setInline( false );
      button->setChecked( true );
      
      samplenums = m_viewer->displayedSamples( kForeground );
      
      for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
          i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[i], data,
                                       samplenums, detnums ) );
    }//if( data )
    
    if( second )
    {
      WRadioButton *button = new WRadioButton( "Second Forground", buttons );
      button->setInline( false );
      
      samplenums = m_viewer->displayedSamples( kSecondForeground );
      for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
          i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[i], second,
                                       samplenums, detnums ) );
      
      if( !data )
        button->setChecked();
      group->addButton( button, -1 );
    }//if( second )

    if( background )
    {
      WRadioButton *button = new WRadioButton( "Background", buttons );
      button->setInline( false );
      
      samplenums = m_viewer->displayedSamples( kBackground );
      for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
          i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[i], background,
                                       samplenums, detnums ) );
      
      group->addButton( button, -1 );
    } // if( background )
  } // if( nspecs > 1 )

  WContainerWidget *linkDiv = new WContainerWidget( dialog->contents() );
  linkDiv->setList( true, false );

  for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
      i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
  {
    const string linktitle = string("As ") + descriptionText(i) + " File";    
    WAnchor *a = new WAnchor( m_specificResources[i], linktitle, linkDiv );
    a->setTarget( TargetNewWindow );
    a->setInline( false );
    a->setStyleClass( "LoadSpectrumSaveAsLink" );

//Either of the two lines below will work to get rid of the save as dialog (at
//  least it looks like it), but I think looking for the link being clicked is
//  a little safer, since I'm not sure how the connections are managed in the
//  case of using the downloadComplete() slot
    a->clicked().connect( boost::bind( &SpecMeasManager::cleanupQuickSaveAsDialog, this, dialog, wApp ) );
//    m_specificResources[i]->downloadComplete().connect( boost::bind( &SpecMeasManager::cleanupQuickSaveAsDialog, this, dialog, wApp ) );
  }//for( SaveSpectrumAsType i = ... )
  

  WPushButton *cancel = new WPushButton( "Cancel", dialog->contents() );
  cancel->setInline( false );
  cancel->clicked().connect( boost::bind( &SpecMeasManager::cleanupQuickSaveAsDialog, this, dialog, wApp ) );
  dialog->finished().connect( boost::bind( &SpecMeasManager::cleanupQuickSaveAsDialog, this, dialog, wApp ) );
  
  dialog->show();
} // void SpecMeasManager::displayQuickSaveAsDialog();


void SpecMeasManager::removeSpecMeas( std::shared_ptr<const SpecMeas> meas,
                                      const bool undisplay )
{
  //implementation copied from removeSpecMeas(Dbo::ptr<UserFileInDb> remove )
  refreshAdditionalInfo( true );
  
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( row );
    if( header->measurementIfInMemory() == meas )
    {
      std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
      
      if( undisplay )
      {
        //also unassign if it was assigned to foreground/2ndf/background
        if( meas == m_viewer->measurment(kForeground) )
          unDisplay(kForeground);
        if( meas== m_viewer->measurment(kSecondForeground) )
          unDisplay(kSecondForeground);
        if( meas == m_viewer->measurment(kBackground) )
          unDisplay(kBackground);
      }//if( undisplay )
      
      if( header->m_aboutToBeDeletedConnection.connected() )
        header->m_aboutToBeDeletedConnection.disconnect();
      
      //This seems to cause problems, so taking this out for now
      //removeFromSpectrumInfoCache( meas, false );
      m_fileModel->removeRows( row, 1 );
      return;
    }// header->dbEntry()==remove
  }//for( int row = 0 row < m_fileModel->rowCount(); ++row )
}//void removeSpecMeas( std::shared_ptr<const SpecMeas> meas )


#if( USE_DB_TO_STORE_SPECTRA )
// Used to unassign and unload spectrums just deleted from the Spectrum Manager
void SpecMeasManager::removeSpecMeas(Dbo::ptr<UserFileInDb> remove )
{
  // Clear the stuff
  refreshAdditionalInfo( true );
  
  //Goes through each loaded state and checks if it's the same one as
  //the removed UserFileInDb
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader>  header = m_fileModel->fileHeader( row );
    if (header->dbEntry()==remove) {
      std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();

      
      //also unassign if it was assigned to foreground/2ndf/background
      if (meas==m_viewer->measurment(kForeground))
        unDisplay(kForeground);

      if (meas==m_viewer->measurment(kSecondForeground))
        unDisplay(kSecondForeground);
      
      if (meas==m_viewer->measurment(kBackground))
        unDisplay(kBackground);
        
        //This seems to cause problems, so taking this out for now
        //removeFromSpectrumInfoCache( meas, false );
      m_fileModel->removeRows( row, 1 );
    }// header->dbEntry()==remove
  }//for( int row = 0 row < m_fileModel->rowCount(); ++row )
} // void SpecMeasManager::removeSelected()
#endif //#if( USE_DB_TO_STORE_SPECTRA )


/**
 Now supports removing multiple selected spectra
**/
void SpecMeasManager::removeSelected()
{
    const WModelIndexSet selected = m_treeView->selectedIndexes();
    WModelIndexSet selectedFiles;
    
    // Clear the stuff
    refreshAdditionalInfo( true );
    
    for( const WModelIndex &index : selected )
    {
        const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
        if( indexLevel == SpectraFileModel::FileHeaderLevel )
            selectedFiles.insert( index );
    } //for( const WModelIndex &index : selected )
    
    if( selectedFiles.size() < 1 )
    {
        cerr << SRC_LOCATION << "\n\tThere are " << selectedFiles.size()
        << " selected files" << endl;
        return;
    } // if( selectedFiles.size() != 1 )

    std::vector<std::shared_ptr<SpectraFileHeader> >  selectedHeaders = getSelectedFiles();
    
    for( vector<std::shared_ptr<SpectraFileHeader> > ::iterator iter = selectedHeaders.begin(); iter!=selectedHeaders.end();iter++)
    {
        if( *iter )
        {
            std::shared_ptr<SpecMeas> meas = (*iter)->measurementIfInMemory();
            removeFromSpectrumInfoCache( meas, false );
            
            //also unassign if it was assigned to foreground/2ndf/background
            if (meas==m_viewer->measurment(kForeground))
                unDisplay(kForeground);
            
            if (meas==m_viewer->measurment(kSecondForeground))
                unDisplay(kSecondForeground);
            
            if (meas==m_viewer->measurment(kBackground))
                unDisplay(kBackground);
        } // if( selectedHeader )
    } //for( vector<std::shared_ptr<SpectraFileHeader> > ::iterator iter = selectedHeaders.begin(); iter!=selectedHeaders.end();iter++)
 
    m_treeView->setSelectedIndexes( WModelIndexSet() );
    for( WModelIndexSet::reverse_iterator saveiter = selectedFiles.rbegin(); saveiter != selectedFiles.rend(); ++saveiter )
    {
        cout<<saveiter->row()<<endl;
        m_fileModel->removeRows( saveiter->row(), 1 );
    } //for( WModelIndexSet::reverse_iterator saveiter = selectedFiles.rbegin(); saveiter != selectedFiles.rend(); ++saveiter )
} // void SpecMeasManager::removeSelected()


void SpecMeasManager::removeAllFiles()
{
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader>  h = m_fileModel->fileHeader( row );
    std::shared_ptr<SpecMeas> meas = h->measurementIfInMemory();
    removeFromSpectrumInfoCache( meas, false );
  }//for( int row = 0 row < m_fileModel->rowCount(); ++row )
  
  m_fileModel->removeRows( 0, m_fileModel->rowCount() );
}//void removeAllFiles()


void SpecMeasManager::renameSaveAsFile()
{
  //if multiple files
  string origName;
  
  set< std::shared_ptr<SpectraFileHeader> > headers;
  const WModelIndexSet selected = m_treeView->selectedIndexes();
  
  for( const WModelIndex &index : selected )
  {
    std::shared_ptr<SpectraFileHeader> header;
    const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
    
    switch( indexLevel )
    {
      case SpectraFileModel::FileHeaderLevel:
        header = m_fileModel->fileHeader( index.row() );
        break;
      
      case SpectraFileModel::SampleLevel:
        header = m_fileModel->fileHeader( index.parent().row() );
        break;
      
      case SpectraFileModel::InvalidLevel:
        break;
    }// switch( level(index) )
      
    if( header )
      headers.insert( header );
  }//for( const WModelIndex &index : selected )
    
  if( headers.size() > 1 )
  {
    origName = "multiple";
  }else
  {
    std::shared_ptr<SpectraFileHeader> selected = selectedFile();
    if( !selected )
      return;
      
    origName = selected->displayName().toUTF8();
    
    const size_t pos = origName.find_last_of( "." );
    if( pos!=string::npos && ((origName.size()-pos)==4) )
      origName = origName.substr( 0, pos );
      
    std::shared_ptr<SpecMeas> meas = selected->measurementIfInMemory();
    const std::set<int> selectedSamples = selectedSampleNumbers();
    if( !!meas && (selectedSamples != meas->sample_numbers()) && selectedSamples.size() )
    {
      if( selectedSamples.size() == 1 )
        origName += "_sample" + std::to_string( *selectedSamples.begin() );
      else
        origName += "_" + std::to_string( selectedSamples.size() ) + "subsamples";
    }//
  }//if( headers.size() > 1 )
  
  
  
  for( SaveSpectrumAsType type = SaveSpectrumAsType(0);
      type < kNumSaveSpectrumAsType;
      type = SaveSpectrumAsType(type+1) )
  {
    const string name = origName + "." + suggestedNameEnding( type );
    m_downloadResources[type]->suggestFileName( name, WResource::Attachment );
  }
  
//    m_saveAsPopup->jsOpen().exec();
} // void SpecMeasManager::renameSaveAsFile()


std::shared_ptr<SpecMeas>
                     SpecMeasManager::selectedToMeasurementInfo() const
{
  std::shared_ptr<SpecMeas> newspec = std::make_shared<SpecMeas>();
  const WModelIndexSet selected = m_treeView->selectedIndexes();

  bool warnAboutRebinning = false;

  // We'll do a simple caching in case SpectraFileHeader isn't caching the files.
  std::shared_ptr<SpectraFileHeader> lastHeader;
  std::shared_ptr<SpecMeas> measurementinfo;

  for( const WModelIndex &index : selected )
  {
    const SpectraFileModel::Level indexLevel = m_fileModel->level(index);

    bool is_sample = (indexLevel == SpectraFileModel::SampleLevel);
    bool is_measurement = is_sample;

    if( !is_sample )
    {
      std::shared_ptr<SpectraFileHeader> thisHeader;
      thisHeader = m_fileModel->fileHeader( index.row() );
      if( thisHeader )
      {
        //We should have already split up measurments with differnt numbers of
        //  bins, therefore I *think* we can just rreturn the whole SpecMeas
        //  object.  IMplemented this 20130118.
        //20160523: I'm not sure what I was thinking, returning the function
        //          here when multiple files are selected, clearly isnt correct.
        //return thisHeader->parseFile();
//        is_measurement = (thisHeader->numSamples() == 1);
      }//if( thisHeader )
    } // if( !is_sample )


    if( is_measurement )
    {
      // We are relying on the fact that if a user selects a file, then
      // all of its measurements are then selected.
      std::shared_ptr<SpectraFileHeader> thisHeader;
      if( is_sample )
        thisHeader = m_fileModel->fileHeader( index.parent().row() );
      else
        thisHeader = m_fileModel->fileHeader( index.row() );

      if( !thisHeader )
      {
        cout << endl << "!thisHeader continuing" << endl;
        continue;
      }

      if( lastHeader != thisHeader )
      {
        lastHeader = thisHeader;
        measurementinfo = thisHeader->parseFile();
      } // if( lastHeader != thisHeader )

      if( !measurementinfo )
      {
       cout << endl << "!measurementinfo  continuing" << endl;
        continue;
      }

      const set<int> sample_numbers_set = measurementinfo->sample_numbers();
      const vector<int> sample_numbers( sample_numbers_set.begin(), sample_numbers_set.end() );
      const vector<int> &detector_numbers = measurementinfo->detector_numbers();

      vector<MeasurementShrdPtr> meas_to_add;

      const int model_row = (is_sample ? index.row() : 0);


      if( model_row >= (int)sample_numbers.size() )
      {
        cerr << SRC_LOCATION << endl << "\t"
             << "Warning, serious logic error here, please correct" << endl;
        cerr << "sample_numbers.size()=" << sample_numbers.size() << endl;
        cerr << "model_row=" << model_row << endl;
        continue;
      }

      for( int detector : detector_numbers )
      {
        MeasurementConstShrdPtr meas = measurementinfo->measurement( sample_numbers[model_row], detector );
        if( meas )
        {
          MeasurementShrdPtr meas_copy = std::make_shared<Measurement>( *meas );
          
          //Add a slightly meaningful title
          const string &oldtitle = meas_copy->title();
          if( oldtitle.empty() )
          {
            string filename = measurementinfo->filename();
            const string::size_type pos = filename.find_last_of('.');
            if( pos != string::npos )
              filename = filename.substr(0, pos);
            
            char buffer[256];
            snprintf( buffer, sizeof(buffer), "%s (sample %i, detector %i)",
                      filename.c_str(), sample_numbers[model_row], detector );
            meas_copy->set_title( buffer );
          }//if( oldtitle.empty() && !measurementinfo->filename().empty() )
          
          meas_to_add.push_back( meas_copy );
        } // if( meas )
      } // for( int detector : detector_numbers )


      vector<MeasurementConstShrdPtr> newspecMeas = newspec->measurements();
      for( const MeasurementShrdPtr &meas : meas_to_add )
      {
        if( newspecMeas.size() )
        {
          if( newspecMeas[0]->channel_energies()
              && meas->channel_energies()
              && meas->channel_energies()->size() //neutron detectors will have a size of zero
              && newspecMeas[0]->channel_energies()->size()
              && newspecMeas[0]->channel_energies()->size() != meas->channel_energies()->size() )
          {
            warnAboutRebinning = true;
          }
        } // if( newspecMeas.size() )

        newspec->add_measurment( meas, false );
      } // for( MeasurementConstShrdPtr &meas : meas_to_add )
    } // if( is_measurement  )
  } // for( const WModelIndex &index : selected )

  newspec->set_filename( "Combination-" + WDateTime::currentDateTime().toString(DATE_TIME_FORMAT_STR).toUTF8() );
  newspec->cleanup_after_load();

  if( warnAboutRebinning )
  {
    const char *msg = "Not all of the spectrum combined had the same number of "
                      "bins so the results may be a bit wonky, sorry :(";
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  } // if( warnAboutRebinning )

  return newspec;
} // std::shared_ptr<SpecMeas> SpecMeasManager::selectedToMeasurementInfo()


void SpecMeasManager::newFileFromSelection()
{
  try
  {
    std::shared_ptr<SpecMeas> spec = selectedToMeasurementInfo();
    std::shared_ptr<SpectraFileHeader> header
                                     = addFile( spec->filename(), spec );
    
    WModelIndex index = m_fileModel->index( header );
    if( index.isValid() )
    {
      WModelIndexSet selected;
      selected.insert( index );
      m_treeView->setSelectedIndexes( selected );
      selectionChanged();
    }// if( index.isValid() )
  }catch( std::exception &e )
  {
    stringstream msg;
    msg << "Failed in combining files: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
#endif
  }
  
  
  /*
  try
  {
    static int customnum = 0;
    string output( wApp->sessionId()
                                    + "_" + std::to_string(customnum++) );

    output = UtilityFunctions::temp_file_name( output, InterSpecApp::tempDirectory() );

    { // Begin codeblock to write file.
      ofstream ofs( output.generic_string().c_str(), ios::binary|ios::out );
      if( !ofs.is_open() )
      {
        stringstream msg;
        msg << SRC_LOCATION << ":\n\tCould not spectrum file - " << output;
        throw runtime_error( msg.str() );
      } // if( fileptr == 0 )

      newspec->write_2012_N42( ofs );
      cout << "Wrote " << output << endl;
    } // End codeblock to write file.

    std::shared_ptr<SpectraFileHeader> header;
    std::shared_ptr<SpecMeas> measurement;
    try
    {
      setFile( newspec->filename(), output, header, measurement, K2012ICD1Parser );
    }catch( std::exception &e )
    {
      stringstream msg;
      msg << "Failed to parse file created from combining multiple spectra: "
          << e.what();
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
#else
      cerr << msg.str() << endl;
#endif
    }

    if( !UtilityFunctions::remove_file( output ) )
    {
      stringstream msg;
      msg << "Failed to remove temporary file: " << output;
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( BOOST_CURRENT_FUNCTION, msg.str().c_str() );
#else
      cerr << msg.str() << endl;
#endif
    }
    
    WModelIndex index = m_fileModel->index( header );
    if( index.isValid() )
    {
      WModelIndexSet selected;
      selected.insert( index );
      m_treeView->setSelectedIndexes( selected );
      selectionChanged();
    } // if( index.isValid() )
  }
  catch( const std::exception &e )
  {
    string msg = "I'm sorry, I had trouble making the combined spectrum file: ";
    msg += e.what();
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
  } // try / catch
*/
} // void SpecMeasManager::newFileFromSelection()


void SpecMeasManager::selectionChanged()
{
  if (!m_spectrumManagerWindow || m_spectrumManagerWindow->isHidden())
      return;
    
  const WModelIndexSet selected = m_treeView->selectedIndexes();
  WModelIndexSet toSelect = selected;
  WModelIndexSet selectedFiles;

  if( selected.empty() )
  {
      if (m_viewer->isSupportFile())
      {
          m_saveButton->disable();
          m_saveButton->setHidden(true, WAnimation());
      } //isSupportFile()
      m_combineButton->setHidden(true,WAnimation());
      m_deleteButton->setHidden(true,WAnimation());
      m_setButton->setHidden(true,WAnimation());
  } // if( selected.empty() )
  else {

  m_deleteButton->enable();
  m_deleteButton->show();
      
  set<WString> files;
  for( const WModelIndex &index : selected )
  {
    std::shared_ptr<const SpectraFileHeader> header;

    const SpectraFileModel::Level indexLevel = m_fileModel->level(index);

    switch( indexLevel )
    {
      case SpectraFileModel::FileHeaderLevel:
      {
        selectedFiles.insert( index );
        header = m_fileModel->fileHeader( index.row() );

       // Let's set all of the children as selected.
       const int nrow = m_fileModel->rowCount( index );
       for( int i = 0; i < nrow; ++i )
         toSelect.insert( index.child(i,0) );
       break;
      } // case SpectraFileModel::FileHeaderLevel:

      case SpectraFileModel::SampleLevel:
        header = m_fileModel->fileHeader( index.parent().row() );
        break;

      case SpectraFileModel::InvalidLevel:
        break; // no-op
        break;
    } // switch( level( index ) )

    if( header )
      files.insert( header->displayName() );
  } // for( const WModelIndex &index : selected )

  if( files.size() > 1 )
  {
    m_combineButton->enable();
    m_combineButton->show();
    m_setButton->setHidden(true, WAnimation());
    m_setButton->disable();
    m_setAsForeground->disable();
    m_setAsBackground->disable();
    m_setAsSecForeground->disable();
  }else
  {
    m_setAsForeground->enable();
    
    if( m_viewer->measurment(kForeground) )
    {
      m_setAsBackground->enable();
      m_setAsSecForeground->enable();
    }else
    {
      m_setAsBackground->disable();
      m_setAsSecForeground->disable();
    } // if( primary spectrums loaded ) / else

    m_combineButton->disable();
    m_combineButton->hide();
    m_setButton->show();
    m_setButton->enable();
  }// if( multiple files are selected ) / else

  if( (selectedFiles.size()==1) && (files.size()==1) )
  {
      if (m_viewer->isSupportFile())
      {
          m_saveButton->enable();
          m_saveButton->show();
      }//isSupportFile()
    // If only one file is selected, update the additional info.
    refreshAdditionalInfo();
  }else
  {
      const bool enableSave = (files.size()==1);
      if (m_viewer->isSupportFile())
      {
          m_saveButton->setEnabled( enableSave );
          m_saveButton->setHidden( !enableSave );
      }
      
      // Otherwise, set it to n/a
      refreshAdditionalInfo( true );
  } // if( selectedFiles.size() == 1 )

  if( selected.size() != toSelect.size() )
    m_treeView->setSelectedIndexes( toSelect );

  }
  // Disable/hide everything and just show what's needed.
//  m_unDisplayMenuButton->hide();
//  m_unDisplayPrimaryButton->   disable();
//  m_unDisplaySecondaryButton-> disable();
//  m_unDisplayBackgroundButton->disable();

  m_removeForeButton->hide();
  m_removeBackButton->hide();
  m_removeFore2Button->hide();
  
  if( m_viewer->measurment(kForeground)
      || m_viewer->measurment(kSecondForeground)
      || m_viewer->measurment(kBackground) )
  {
//    m_unDisplayMenuButton->show();
  }
  if( m_viewer->measurment(kForeground) ) {
//    m_unDisplayPrimaryButton->enable();
    m_removeForeButton->show();
  }
  if( m_viewer->measurment(kSecondForeground) ) {
//    m_unDisplaySecondaryButton->enable();
    m_removeFore2Button->show();
  }
  if( m_viewer->measurment(kBackground) ) {
//    m_unDisplayBackgroundButton->enable();
    m_removeBackButton->show();
  }
} // void SpecMeasManager::selectionChanged()


WText *messageCombiner( const string& label, const string& value )
{
  WText *boxxed = new WText( "<b>" + label + ":</b> " + value );
  boxxed->setPadding( 2 );

  return boxxed;
} // WText *messageCombiner( const string&, const string& )


void SpecMeasManager::refreshAdditionalInfo( bool clearInfo )
{
  std::shared_ptr< SpecMeas > selected = selectedToMeasurementInfo();

  vector< string > info;
  // If it /isn't/ a clearInfo, there is just one thing selected.
  if( selected && !clearInfo )
  {
    const vector< MeasurementConstShrdPtr > selectedData = selected->measurements();
    // A handful of these assume that the data will be populated for all of the detectors.
    // That is, it will just take the first one out of the detector for the dates, as well
    // as the number of channels.

    // This accounts for potential "mass" samples.
    if( selectedData.size() > 0 )
    {
      info.push_back( boost::posix_time::to_simple_string(
        selectedData.at(0)->start_time() ).substr( 0, 11 ) ); // Date taken
      stringstream liveStream, realStream, sizeStream;
      liveStream << selected->gamma_live_time();
      realStream << selected->gamma_real_time();
      sizeStream << selectedData.at(0)->gamma_counts()->size();
      info.push_back( liveStream.str() ); // Live time
      info.push_back( realStream.str() ); // Real time
      info.push_back( sizeStream.str() ); // Number of channels
      stringstream gammaCounts, neutronCounts, detectorCounts;
      gammaCounts << fixed << setprecision(0) << selected->gamma_count_sum();
      neutronCounts << fixed << setprecision(0) << selected->neutron_counts_sum();
      detectorCounts << selected->detector_names().size();
      info.push_back( gammaCounts.str() ); // Gamma counts
      info.push_back( neutronCounts.str() ); // Neutron counts
      info.push_back( detectorCounts.str() ); // Detector counts
    }
    else
    {
      clearInfo = true;
    }
  }//if( !clearInfo )

  int infoPerRow = 3;

  WGridLayout *boxxyHost = new WGridLayout();

  for( unsigned int i = 0; i < m_labelTags.size(); ++i )
  {
    int cS = 225 - 42 * !( i % 2 );
    WText *boxxy = NULL;

    if( clearInfo || i>=info.size() )
      boxxy = messageCombiner( m_labelTags.at(i), "n/a" );
    else
      boxxy = messageCombiner( m_labelTags.at(i), info.at(i) );

    boxxy->decorationStyle().setBackgroundColor( WColor( cS, cS, cS ) );
    boxxyHost->addWidget( boxxy, (int)( i / infoPerRow ), i % infoPerRow );
  }//for( unsigned int i = 0; i < m_labelTags.size(); ++i )

  boxxyHost->setVerticalSpacing( 8 );
  boxxyHost->setHorizontalSpacing( 10 );
//  boxxyHost->setRowStretch( 0, -1 );
//  boxxyHost->setRowStretch( 1, -1 );
//  boxxyHost->setRowStretch( 2, -1 );

  //Before setting a new layout to m_infoHandler, we have to avoid memorry leaks
  //  by explicitly deleteing the children of the previous layout, as this is
  //  not done by Wt (as noted in documentation for WContainerWidget::setLayout)
  WLayout *prevLayout = m_infoHandler->layout();
  if( prevLayout )
    prevLayout->clear();

  m_infoHandler->setLayout( boxxyHost );
} // void SpecMeasManager::refreshAdditionalInfo()


void SpecMeasManager::createInfoHandler()
{
  m_infoHandler = new WContainerWidget();
  m_infoHandler->setStyleClass( "infoHandler" );
  m_infoHandler->setOverflow( WContainerWidget::OverflowAuto );

  m_labelTags.push_back( "Date taken" );
  m_labelTags.push_back( "Live time (s)" );
  m_labelTags.push_back( "Real time (s)" );
  m_labelTags.push_back( "# of channels" );
  m_labelTags.push_back( "Gamma counts" );
  m_labelTags.push_back( "Neutron counts" );
  m_labelTags.push_back( "# of detectors" );

  refreshAdditionalInfo( true );
}//void createInfoHandler()


void SpecMeasManager::deleteSpectrumManager()
{
    m_spectrumManagertreeDiv->removeWidget(m_treeView);
    AuxWindow::deleteAuxWindow(m_spectrumManagerWindow);
    m_spectrumManagerWindow=NULL;
}//deleteSpectrumManager()


WContainerWidget *SpecMeasManager::createTreeViewDiv()
{
  assert( m_treeView );
  assert( m_fileModel );

  m_spectrumManagertreeDiv = new WContainerWidget();
  m_spectrumManagertreeDiv->setStyleClass( "treeDiv" );
  m_spectrumManagertreeDiv->setOverflow( WContainerWidget::OverflowAuto );

  m_treeView->setHeight( WLength(98.0,WLength::Percentage) );
  m_spectrumManagertreeDiv->addWidget( m_treeView );

  

  m_treeView->setSortingEnabled( true );
  m_treeView->setSelectionMode( ExtendedSelection /*SingleSelection*/ );
  m_treeView->setColumn1Fixed( false );
  m_treeView->setColumnWidth( SpectraFileModel::kDisplayName,     WLength( 17, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kUploadTime,      WLength( 16, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kNumMeasurements, WLength( 10, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kLiveTime,        WLength( 14, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kRealTime,        WLength( 14, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kGammaCounts,     WLength( 15, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kNeutronCounts,   WLength( 15, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kSpectrumTime,    WLength( 18, WLength::FontEx ) );
  m_treeView->setColumnWidth( SpectraFileModel::kNumDetectors,    WLength( 10, WLength::FontEx ) );

  WItemDelegate *delegate = new WItemDelegate( m_treeView );
  delegate->setTextFormat( "%.2f" );
  m_treeView->setItemDelegateForColumn( SpectraFileModel::kLiveTime, delegate );
  m_treeView->setItemDelegateForColumn( SpectraFileModel::kRealTime, delegate );
  m_treeView->setItemDelegateForColumn( SpectraFileModel::kGammaCounts, delegate );

  return m_spectrumManagertreeDiv;
} // WContainerWidget *SpecMeasManager::createTreeViewDiv()


WContainerWidget *SpecMeasManager::createButtonBar()
{
//  assert( !m_unDisplayMenuButton );
//  assert( !m_setPrimaryButton );
  // ...

  WContainerWidget *buttonBar = new WContainerWidget();
  buttonBar->setStyleClass( "SpectraFileManagerButtonDiv" );

  WGridLayout *buttonAlignment = new WGridLayout();
  buttonAlignment->setContentsMargins( 0, 0, 0, 0 );
  WContainerWidget* buttonHost = new WContainerWidget( buttonBar );
  buttonHost->setLayout( buttonAlignment );
  
  // ---- try new bar ----
  WContainerWidget *m_newDiv = new WContainerWidget( );
  buttonAlignment->addWidget( m_newDiv, 1, 0 );
  m_newDiv->setStyleClass( "LoadSpectrumUploadDiv" );
  new WText( "Selected Spectrum: " , m_newDiv );
  
  m_setButton = new Wt::WPushButton("Assign As",m_newDiv);
  m_setButton->setIcon( "InterSpec_resources/images/bullet_arrow_down.png" );
  
  Wt::WPopupMenu *setPopup = new Wt::WPopupMenu();
  m_setAsForeground = setPopup->addItem("Foreground");
  m_setAsForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, kForeground, true) );
  m_setAsBackground = setPopup->addItem("Background");
  m_setAsBackground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, kBackground, true) );
  m_setAsSecForeground = setPopup->addItem("Secondary Foreground");
  m_setAsSecForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, kSecondForeground, true ) );
  m_setButton->setMenu(setPopup);
  m_setButton->setHidden( true, WAnimation() );
//  m_removeFileButton->disable();

  
  m_combineButton = new Wt::WPushButton("Combine",m_newDiv);
  m_combineButton->addStyleClass("JoinIcon");
  m_combineButton->setIcon( "InterSpec_resources/images/arrow_join.png" );
  m_combineButton->clicked().connect( boost::bind( &SpecMeasManager::newFileFromSelection, this ) );
  m_combineButton->setHidden( true, WAnimation() );

  if( m_viewer->isSupportFile() )
  {
    m_saveButton = new Wt::WPushButton("Export",m_newDiv);
    m_saveButton->setIcon( "InterSpec_resources/images/bullet_arrow_down.png" );
    Wt::WPopupMenu *setPopup2 = new Wt::WPopupMenu();
        
    m_saveButton->mouseWentOver().connect( boost::bind( &SpecMeasManager::renameSaveAsFile, this ) );
      
    for( SaveSpectrumAsType type = SaveSpectrumAsType(0);
         type < kNumSaveSpectrumAsType;
         type = SaveSpectrumAsType(type+1) )
    {
      const string descrip = string("As ") + descriptionText(type) + " File";
      WMenuItem *temp = setPopup2->addItem( descrip );
      temp->setLink( m_downloadResources[type] );
      temp->setLinkTarget( TargetNewWindow );
    }
      
    m_saveButton->setMenu( setPopup2 );
    m_saveButton->setHidden( true, WAnimation() );
  }//isSupportFile
  
  m_deleteButton =  new Wt::WPushButton("Unload",m_newDiv);
  m_deleteButton->setIcon( "InterSpec_resources/images/minus_min_white.png" );
  m_deleteButton ->clicked().connect( boost::bind( &SpecMeasManager::removeSelected, this ) );
  m_deleteButton->setHidden( true, WAnimation() );
  
  // ---- try new bar 2 ----
  WContainerWidget *m_newDiv2 = new WContainerWidget( );
  buttonAlignment->addWidget( m_newDiv2, 2, 0);
  m_newDiv2->setStyleClass( "LoadSpectrumUploadDiv" );
  
  //WText *text =
  new WText( "Unassign: ", m_newDiv2 );
  
  m_removeForeButton = new Wt::WPushButton("Foreground",m_newDiv2);
//      m_removeForeButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeForeButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, kForeground) );
  m_removeForeButton  ->setHidden( true, WAnimation() );
//  m_removeForeButton->setHiddenKeepsGeometry( true );
  
  m_removeBackButton = new Wt::WPushButton("Background",m_newDiv2);
//          m_removeBackButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeBackButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, kBackground ) );
  m_removeBackButton  ->setHidden( true, WAnimation() );
//  m_removeBackButton->setHiddenKeepsGeometry( true );
  
  m_removeFore2Button = new Wt::WPushButton("Secondary Foreground",m_newDiv2);
//              m_removeFore2Button->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeFore2Button->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, kSecondForeground ) );
  m_removeFore2Button  ->setHidden( true, WAnimation() );
//  m_removeFore2Button->setHiddenKeepsGeometry( true );
  
//  Wt::WPushButton *m_removeButton = new Wt::WPushButton("Remove",m_newDiv);
  
  return buttonBar;
} // WContainerWidget *SpecMeasManager::createButtonBar()


SpectraFileModel *SpecMeasManager::model()
{
  return m_fileModel;
} // SpectraFileModel *SpecMeasManager::model()


const SpectraFileModel *SpecMeasManager::model() const
{
  return m_fileModel;
} // const SpectraFileModel *SpecMeasManager::model() const


RowStretchTreeView *SpecMeasManager::treeView()
{
  return m_treeView;
} // RowStretchTreeView *SpecMeasManager::treeView()

const RowStretchTreeView *SpecMeasManager::treeView() const
{
  return m_treeView;
} // const RowStretchTreeView *SpecMeasManager::treeView() const


const InterSpec *SpecMeasManager::viewer() const
{
  return m_viewer;
} // const InterSpec *SpecMeasManager::viewer() const

void SpecMeasManager::displayIsBeingShown()
{
//  m_displayShowing = true;
//  clearTempSpectrumInfoCache();
} // void SpecMeasManager::displayIsBeingShown()


void SpecMeasManager::displayIsBeingHidden()
{
//  m_displayShowing = false;
  clearTempSpectrumInfoCache();
} // void SpecMeasManager::displayIsBeingHidden()

#if( USE_DB_TO_STORE_SPECTRA )
void SpecMeasManager::saveToDatabase(
                                std::shared_ptr<const SpecMeas> input ) const
{
  //Make sure we are in the event loop - for thread safety.
  if( !wApp )
    throw runtime_error( "SpecMeasManager::saveToDatabase() must be called from within the event loop." );
  
  std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( input );
  
  if( !!header )
  {
    if( !header->shouldSaveToDb() )
      return;
    
//    boost::function<void(void)> worker
//                      = app->bind( boost::bind( &SpectraFileHeader::saveToDatabaseWorker, meas, header ) );
//    WServer::instance()->post( app->sessionId(), worker );
    
    std::shared_ptr<SpecMeas> meas = header->parseFile();
    boost::function<void(void)> worker = boost::bind( &SpectraFileHeader::saveToDatabaseWorker, meas, header );
    WServer::instance()->ioService().post( worker );
  }//if( headermeas && (headermeas==meas) )
}//void saveToDatabase( std::shared_ptr<const SpecMeas> meas ) const


int SpecMeasManager::setDbEntry( Wt::Dbo::ptr<UserFileInDb> dbfile,
                                 std::shared_ptr<SpectraFileHeader> &header,
                                 std::shared_ptr<SpecMeas> &measurement,
                                 bool enforceUser )
{
  if( !dbfile )
    throw runtime_error( "SpecMeasManager::setDbEntry(...): invalid dbentry" );
  
  if( enforceUser && (dbfile->user != m_viewer->m_user) )
    throw runtime_error( "SpecMeasManager::setDbEntry(...): invalid user" );
  
  int row = -1;
  header.reset();
  measurement.reset();
  header.reset( new SpectraFileHeader( m_viewer->m_user,
                                       false, m_viewer ) );
  measurement = header->resetFromDatabase( dbfile );
  addToTempSpectrumInfoCache( measurement );
  row = m_fileModel->addRow( header );
  return row;
}//int setDbEntry(...)


void SpecMeasManager::userCanceledResumeFromPreviousOpened( AuxWindow *window,
                                  std::shared_ptr<SpectraFileHeader> header )
{
  if( window )
    delete window;

  std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
  
  boost::function<void(void)> f;
  
  if( meas )
    f = boost::bind( &SpectraFileHeader::saveToDatabaseWorker, meas, header );
  else
    f = boost::bind( &SpectraFileHeader::saveToDatabaseFromTempFileWorker,
                     header );
  
//  boost::function<void(void)> worker = wApp->bind( f );
//  WServer::instance()->post( wApp->sessionId(), worker );
  WServer::instance()->ioService().post( f );
}//userCanceledResumeFromPreviousOpened(..)


void SpecMeasManager::createPreviousSpectraDialog( const std::string sessionID,
                                                  std::shared_ptr<SpectraFileHeader> header,
                                                  const SpectrumType type,
                                                  const std::vector< Wt::Dbo::ptr<UserFileInDb> > modifiedFiles,
                                                  const std::vector< Wt::Dbo::ptr<UserFileInDb> > unModifiedFiles )
{
  AuxWindow *message = new AuxWindow( "File Viewed Before",
                          (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                           | AuxWindowProperties::TabletModal
                           | AuxWindowProperties::DisableCollapse) );
  message->setResizable( false );
  message->rejectWhenEscapePressed();
  //      message->setMaximumSize( WLength::Auto,
  //                               WLength( 85.0, WLength::Percentage ) );
  message->contents()->addStyleClass( "PreviousDbDialogContents" );
  WGridLayout *layout = new WGridLayout(message->contents());
  layout->setContentsMargins(5,5,0,0);
  message->contents()->setOverflow(WContainerWidget::OverflowHidden);
  const char *msg = "It looks like you have previously loaded and modified "
  "this spectrum file, would you like to resume "
  "your previous work?";
  WText *txt = new WText( msg, Wt::XHTMLText);
  layout->addWidget(txt,0,0);
  txt->setAttributeValue( "style", "font-size:125%;font-weight:bold;"
                         "font-family:\"Times New Roman\", Times, serif;" );
  txt->setInline( false );

  WContainerWidget* container = new WContainerWidget();
  container->setOverflow(WContainerWidget::OverflowAuto);
  layout->addWidget(container,1,0);
  layout->setRowStretch(1, 1);
  for( size_t i = 0; i < modifiedFiles.size(); ++i )
  {
    if( !!modifiedFiles[i] )
      new PreviousDbEntry( message, container, type, m_fileModel, this,
                           modifiedFiles[i], header );
  }//for( size_t i = 0; i < unModifiedFiles.size(); ++i )
  
  if( unModifiedFiles.size() )
  {
    Dbo::ptr<UserFileInDb> f = unModifiedFiles.front();
    cerr << "Setting file with UUID=" << header->m_uuid << ", and filename "
    << header->m_displayName << " to be connected to DB entry from "
    << "upload at "
    << f->uploadTime.toString(DATE_TIME_FORMAT_STR).toUTF8()
    << " until user determines if they want to resume from a modified"
    << " session" << endl;
    header->setDbEntry( f );
  }//if( unModifiedFiles.size() )
  
 
//  message->footer()->resize( WLength::Auto, WLength(50.0) );
  
  WPushButton *cancel = message->addCloseButtonToFooter();
  
  if( unModifiedFiles.empty() )
  {
    //        int backid = -1;
    Dbo::ptr<UserFileInDb> dpback;
    std::shared_ptr<SpecMeas>  old_back;
    if( type == kForeground )
      old_back = m_viewer->measurment( kBackground );
    //        if( old_back )
    //          backid = m_fileModel->dbEntry( old_back ).id();
    
    cancel->clicked().connect( message, &AuxWindow::hide );
    //          boost::bind( &SpecMeasManager::userCanceledResumeFromPreviousOpened,
    //                       this, message, header ) );
    message->finished().connect(
                                boost::bind( &SpecMeasManager::userCanceledResumeFromPreviousOpened,
                                            this, message, header ) );
  }else
  {
    cancel->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, message ) );
    message->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, message ) );
  }//if( unModifiedFiles.empty() ) / else
  
  //It would be nice to make it so that the created dialog will go away if
  //  another spectrum file is loaded
//  m_viewer
  
  message->show();
  message->resize(WLength(500), WLength(80,WLength::Percentage));
  message->centerWindow();
  message->refresh();
  wApp->triggerUpdate();
}//void createPreviousSpectraDialog(...)

void postErrorMessage( const string msg, const WarningWidget::WarningMsgLevel level )
{
  passMessage( msg, "", level );
  wApp->triggerUpdate();
}//void postErrorMessage( string msg )


void SpecMeasManager::checkIfPreviouslyOpened( const std::string sessionID,
                                 std::shared_ptr<SpectraFileHeader> header,
                                 SpectrumType type,
                                 std::shared_ptr< std::mutex > mutex,
                                 std::shared_ptr<bool> destructed )
{
  std::lock_guard<std::mutex> lock( *mutex );
  
  if( *destructed )
  {
    cerr << "checkIfPreviouslyOpened(): mamnager destructed before I could do ish" << endl;
    return;
  }
  
  try
  {
    const bool storeInDb
      = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );

    if( !storeInDb )
      return;
    
    if( !header )
      throw runtime_error( "Invalid SpectraFileHeader passed in" );

    if( !m_viewer || !m_viewer->m_user )
      throw runtime_error( "Invalid InterSpec or user pointer" );
    
    vector< Wt::Dbo::ptr<UserFileInDb> > modifiedFiles, unModifiedFiles;
    
    {//begin interaction with database
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
      DataBaseUtils::DbTransaction transaction( *sql );
      
      typedef Dbo::collection< Dbo::ptr<UserFileInDb> > UserFileInDbColl;

//      UserFileInDbColl files = m_viewer->m_user->userFiles().find()
//                                          .where( "UUID = ? AND Filename = ? "
//                                                  "AND IsPartOfSaveState = 0" )
//                                          .bind( header->m_uuid )
//                                          .bind( header->m_displayName );
      UserFileInDbColl files = m_viewer->m_user->userFiles().find()
                                 .where( "UUID = ? AND IsPartOfSaveState = 0" )
                                 .bind( header->m_uuid );

      for( UserFileInDbColl::iterator i = files.begin(); i != files.end(); ++i )
      {
        if( (*i)->userHasModified )
          modifiedFiles.push_back( *i );
        else
          unModifiedFiles.push_back( *i );
      }//for( loop over matching files in DB )
   
      m_viewer->m_user.modify()->incrementSpectraFileOpened();
      transaction.commit();
    }//end interaction with database
    

    if( modifiedFiles.empty() && unModifiedFiles.empty() )
    {
      cerr << "File with UUID=" << header->m_uuid << ", and filename "
           << header->m_displayName << " has not been saved to DB before"
           << endl;
      
      try
      {
        if( storeInDb && header->shouldSaveToDb() )
        {
          //see if there is a background measurment
          std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
          if( meas )
            SpectraFileHeader::saveToDatabase( meas, header );
          else
            header->saveToDatabaseFromTempFile();
        }//if( header->shouldSaveToDb() )
      }catch( FileToLargeForDbException &e )
      {
        string msg = "Will not be able to save this file to the database. ";
        msg += e.what();
        cerr << msg << endl;
        passMessage( msg, "", WarningWidget::WarningMsgInfo );
      }
      return;
    }//if( !files.size() )

    if( modifiedFiles.empty() )
    {
      Dbo::ptr<UserFileInDb> dbentry = unModifiedFiles.front();
      cerr << "Setting file with UUID=" << header->m_uuid << ", and filename "
           << header->m_displayName << " to be connected to DB entry from "
           << "upload at "
           << dbentry->uploadTime.toString(DATE_TIME_FORMAT_STR).toUTF8()
           << endl;
      WServer::instance()->post( sessionID, boost::bind( &setHeadersDbEntry, header, dbentry) );
      
      return;
    }else if( !!header )
    {
        
      WServer::instance()->post( sessionID,
                                  boost::bind( &SpecMeasManager::browseDatabaseSpectrumFiles,
                                              this, header->m_uuid, type, header) );

        
      /*
       Old way shows "File Viewed Before" dialog:
       
      WServer::instance()->post( sessionID,
                    boost::bind( &SpecMeasManager::createPreviousSpectraDialog,
                                 this, sessionID, header, type, modifiedFiles,
                                 unModifiedFiles ) );
       */
    }
  }catch( std::exception & )
  {
    WServer::instance()->post( sessionID,
                              boost::bind( &postErrorMessage,
              string("Error checking if this file had been opened previously"),
              WarningWidget::WarningMsgHigh ) );
  }//try / catch
}//void checkIfPreviouslyOpened(...)
#endif  //#if( USE_DB_TO_STORE_SPECTRA )


std::shared_ptr<SpectraFileHeader> SpecMeasManager::addFile( const std::string &displayName,
                              std::shared_ptr<SpecMeas> measurement )
{
  if( !measurement )
    throw runtime_error( "SpecMeasManager::addFile(): invalid input" );
  
  std::shared_ptr<SpectraFileHeader> header
   = std::make_shared<SpectraFileHeader>( m_viewer->m_user, false, m_viewer );
  header->setMeasurmentInfo( measurement );
  
  addToTempSpectrumInfoCache( measurement );
  m_fileModel->addRow( header );
  
  return header;
}//std::shared_ptr<SpectraFileHeader> SpecMeasManager::addFile(...)


int SpecMeasManager::setFile( const std::string &displayName,
                              const std::string &filename,
                              std::shared_ptr<SpectraFileHeader> &header,
                              std::shared_ptr<SpecMeas> &measurement,
                              ParserType parser_type )
{
  //We are relying on SpectraFileHeader::setFile(...) to throw an exception
  //  whenever it isnt successful in loading the file, because we want this
  //  function to throw an exception upon failure
  int row = -1;
  header.reset();
  measurement.reset();
  
  std::shared_ptr<SpectraFileHeader> new_header
     = std::make_shared<SpectraFileHeader>( m_viewer->m_user, false, m_viewer );
  
  std::shared_ptr<SpecMeas> new_measurement
                   = new_header->setFile( displayName, filename, parser_type );
  
  header = new_header;
  measurement = new_measurement;
  
  addToTempSpectrumInfoCache( measurement );
  
  row = m_fileModel->addRow( header );
  
  return row;
} // int SpecMeasManager::setFile(...)


int SpecMeasManager::dataUploaded2( Wt::WFileUpload *upload , SpectrumType type)
{
  std::shared_ptr<SpecMeas> measurement;
  int row= dataUploaded( upload, measurement );
  displayFile( row, measurement, type, true, true, true );
  return row;
} // int SpecMeasManager::dataUploaded( Wt::WFileUpload )


int SpecMeasManager::dataUploaded( Wt::WFileUpload *upload )
{
  std::shared_ptr<SpecMeas> measurement;
  return dataUploaded( upload, measurement );
} // int SpecMeasManager::dataUploaded( Wt::WFileUpload )


bool SpecMeasManager::loadFromFileSystem( const string &name, SpectrumType type,
                                       ParserType parseType )
{
  const string origName = UtilityFunctions::filename( name );
  
  try
  {
    std::shared_ptr<SpecMeas> measurement;
    std::shared_ptr<SpectraFileHeader> header;
    int row = setFile( origName, name, header, measurement, parseType );
    if( row < 0 )
      throw runtime_error( "invalid file" );
     
    WModelIndexSet selected;
    WModelIndex index = m_fileModel->index( row, 0 );
    selected.insert( index );
    m_treeView->setSelectedIndexes( WModelIndexSet() );    
//    passMessage( "Successfully uploaded file.", "", 0 );

    displayFile( row, measurement, type, true, true, true );
  }catch( const std::exception &e )
  {
#if( SUPPORT_ZIPPED_SPECTRUM_FILES )
    {
      ifstream test( name.c_str(), ios::in | ios::binary );
      const bool iszip = (test.get()==0x50 && test.get()==0x4B
                          && test.get()==0x03 && test.get()==0x04);
      test.close();
      
      if( iszip
         /*&& UtilityFunctions::iequals( origName.substr(origName.length()-4), ".zip")*/
         && handleZippedFile( origName, name, type ) )
      return true;
    }
#endif
    
    if( !handleNonSpectrumFile( origName, name ) )
    {
      displayInvalidFileMsg(origName,e.what());
    }
    return false;
  }// try/catch
  
  return true;
}//void loadFromFileSystem( std::string filename )


int SpecMeasManager::dataUploaded( Wt::WFileUpload *upload,
                            std::shared_ptr<SpecMeas> &measurement )
{
  const string fileName = upload->spoolFileName();
  const WString clientFileName = upload->clientFileName();
  const string origName = clientFileName.toUTF8();

  try
  {
    std::shared_ptr<SpectraFileHeader> header;
    int result = setFile( origName, fileName, header, measurement, kAutoParser );
    WModelIndexSet selected;
    WModelIndex index = m_fileModel->index( result, 0 );
    selected.insert( index );
    m_treeView->setSelectedIndexes( WModelIndexSet() );

    passMessage( "Successfully opened file.", "dataUploaded", 0 );

    return result;
  }catch( const std::exception &e )
  {
    displayInvalidFileMsg(origName,e.what());
  }// try/catch

  return -1;
} // int SpecMeasManager::dataUploaded()


void SpecMeasManager::clearTempSpectrumInfoCache()
{
  typedef std::deque< std::shared_ptr<const SpecMeas> > queue_type;

#if( USE_DB_TO_STORE_SPECTRA )
  const bool storeInDb
  = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
#endif

  for( queue_type::iterator iter = m_tempSpectrumInfoCache.begin();
      iter != m_tempSpectrumInfoCache.end(); ++iter )
  {
    serializeToTempFile( *iter );
    
#if( USE_DB_TO_STORE_SPECTRA )
    if( storeInDb )
      saveToDatabase( *iter );
#endif
  }//for( loop over m_tempSpectrumInfoCache to save them to disk )
    
  m_tempSpectrumInfoCache.clear();
} // void SpecMeasManager::clearTempSpectrumInfoCache()


void SpecMeasManager::serializeToTempFile( std::shared_ptr<const SpecMeas> meas ) const
{
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader(row);
    std::shared_ptr<SpecMeas> headermeas = header->measurementIfInMemory();
    if( headermeas && (headermeas==meas) )
    {
      header->saveToFileSystem( headermeas );
      
/*
#if( USE_DB_TO_STORE_SPECTRA )
    if( header->m_app && header->shouldSaveToDb() )
    {
//      boost::function<void(void)> worker
//                        = header->m_app->bind( boost::bind(
//                                      &SpectraFileHeader::saveToDatabaseWorker,
//                                      headermeas, header ) );
//      WServer::instance()->post( header->m_app->sessionId(), worker );
      boost::function<void(void)> worker
                      = boost::bind( &SpectraFileHeader::saveToDatabaseWorker,
                                     headermeas, header );
      WServer::instance()->ioService().post( worker );
    }//if( header->m_app && header->shouldSaveToDb() )
#endif
*/
      
      return;
    }//if( headermeas && (headermeas==meas) )
  }//for( int row = 0; row < m_fileModel->rowCount(); ++row )
}//void serializeToTempFile( std::shared_ptr<const SpecMeas> meas ) const



void SpecMeasManager::removeFromSpectrumInfoCache( std::shared_ptr<const SpecMeas> meas,
                                                   bool saveToDisk ) const
{
  if( !meas )
    return;

  typedef std::deque< std::shared_ptr<const SpecMeas> > queue_type;

  //Only put the ptr in the queue if its not already in there
  const queue_type::iterator pos = std::find( m_tempSpectrumInfoCache.begin(),
                                               m_tempSpectrumInfoCache.end(),
                                               meas );
  if( pos == m_tempSpectrumInfoCache.end() )
    return;

  if( saveToDisk )
    serializeToTempFile( meas );
  
#if( USE_DB_TO_STORE_SPECTRA )
  const bool storeInDb
     = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
  if( saveToDisk && storeInDb )
    saveToDatabase( meas );
#endif
  
  m_tempSpectrumInfoCache.erase( pos );
}//void SpecMeasManager::removeFromSpectrumInfoCache(...) const


void SpecMeasManager::addToTempSpectrumInfoCache( std::shared_ptr<const SpecMeas> meas ) const
{
  if( sm_maxTempCacheSize == 0 )
    return;

  if( !meas )
    return;

  typedef std::deque< std::shared_ptr<const SpecMeas> > queue_type;

  //Only put the ptr in the queue if its not already in there
  const queue_type::iterator pos = std::find( m_tempSpectrumInfoCache.begin(),
                                               m_tempSpectrumInfoCache.end(),
                                               meas );

  if( pos == m_tempSpectrumInfoCache.end() )
  {
    m_tempSpectrumInfoCache.push_back( meas );
  }else
  {
    return; //maybe shouldnt do this because this function isnt real real expensive?
  }
  
#if( USE_DB_TO_STORE_SPECTRA )
  const bool storeInDb
      = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
#endif

  size_t curr_size = 0;
  for( queue_type::reverse_iterator riter = m_tempSpectrumInfoCache.rbegin();
       riter != m_tempSpectrumInfoCache.rend();
       ++riter )
  {
    curr_size += (*riter)->memmorysize();
    if( curr_size > sm_maxTempCacheSize )
    {
      //lets give the SpectraFileHeader class a chance to write the SpecMeas
      //  objects we're about to delete, to disk.  This avoids writing them
      //  to disk in the SpecMeas destructors thread
      for( queue_type::iterator saveiter = m_tempSpectrumInfoCache.begin();
          saveiter != riter.base(); ++saveiter )
      {
        serializeToTempFile( *saveiter );
        
#if( USE_DB_TO_STORE_SPECTRA )
        if( storeInDb )
          saveToDatabase( *saveiter );
#endif
      }//for( loop over measurments about to be removed from cache )
      
      m_tempSpectrumInfoCache.erase( m_tempSpectrumInfoCache.begin(),
                                    riter.base() );
      break;
    } // if( curr_size > sm_maxTempCacheSize )
  } // for( loop over m_tempSpectrumInfoCache )
} // void SpecMeasManager::addToTempSpectrumInfoCache()



void SpecMeasManager::cleanupQuickSaveAsDialog( AuxWindow *dialog, WApplication *app )
{
  if( !dialog || !app )
    return;

  WApplication::UpdateLock lock( app );

  if( lock )
  {
    dialog->accept();
    delete dialog;
    std::shared_ptr<const SpecMeas> dummy;
    
    for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
        i < kNumSaveSpectrumAsType; i = SaveSpectrumAsType(i+1) )
      m_specificResources[i]->setSpectrum( dummy, set<int>(), set<int>() );
    
    app->triggerUpdate();
  }else
    cerr << endl << SRC_LOCATION
         << "\n\tFailed to get app lock - not doing anything!\n\n" << endl;
}//void cleanupQuickSaveAsDialog(...)


void SpecMeasManager::fileTooLarge( const ::int64_t size_tried )
{
  const int max_size = static_cast<int>( WApplication::instance()->maximumRequestSize() );
  stringstream msg;
  msg << "Attempted file is too large; max size is " << (max_size/1012)
      << " you tried to upload " << (size_tried/1012) << " kb";
  passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
} // void SpecMeasManager::fileTooLarge()


void SpecMeasManager::uploadSpectrum() {
  new UploadBrowser(this/*, m_viewer*/);
}


#if( USE_DB_TO_STORE_SPECTRA )
void SpecMeasManager::browseDatabaseSpectrumFiles(std::string uuid,  SpectrumType type, std::shared_ptr<SpectraFileHeader> header)
{
  new DbFileBrowser( this, m_viewer , uuid, type , header);
}//void browseDatabaseSpectrumFiles()


void SpecMeasManager::finishStoreAsSpectrumInDb( Wt::WLineEdit *nameWidget,
                                               Wt::WTextArea *descWidget,
                                        std::shared_ptr<SpecMeas> meas,
                                               AuxWindow *window )
{
  const WString name = nameWidget->text();
  const WString desc = descWidget ? descWidget->text() : "";
  
  if( name.empty() )
  {
    passMessage( "You must specify a name", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( name.empty() )
  
  if( window )
    delete window;
  
  std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( meas );
  if( !header )
  {
    cerr << SRC_LOCATION << "\n\tFailed to save file to DB." << endl;
    passMessage( "Error saving to database", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !header )
  
  header->setDbEntry( Dbo::ptr<UserFileInDb>() );
  
  SpectraFileHeader::saveToDatabase( meas, header );
  
  Wt::Dbo::ptr<UserFileInDb> entry = header->dbEntry();
  

  if( !entry )
  {
    cerr << SRC_LOCATION << "\n\tFailed to save file to DB, apparently" << endl;
    passMessage( "Error saving to database", "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  try
  {
    DataBaseUtils::DbTransaction transaction( *m_sql );
    entry.modify()->filename = name.toUTF8();
    entry.modify()->description = desc.toUTF8();
    
    transaction.commit();
//    passMessage( "Stored " + entry->filename, "", WarningWidget::WarningMsgInfo );
  }catch( std::exception &e )
  {
    cerr << "SpecMeasManager::finishStoreSpectrumInDb() caught: " << e.what()
         << endl;
    passMessage( "Filename or description you entered may not have been used"
                 " when saving the spectrum", "", WarningWidget::WarningMsgHigh );
  }//try / catc
}//void finishStoreAsSpectrumInDb(...)

//Menu: Spectrum->Save Spectrum
void SpecMeasManager::storeSpectraInDb()
{
  for( int i = 0; i < 3; ++i )
  {
    const SpectrumType type = SpectrumType( i );
    std::shared_ptr<SpecMeas> m = m_viewer->measurment( type );

    if( m )
    {
      saveToDatabase( m );  //happens in another thread
      std::shared_ptr<SpectraFileHeader> header= m_fileModel->fileHeader( m );
//      WString msg = "Stored '";
//      msg += header ? header->displayName() : WString("");
//      msg += "'";
//      passMessage( msg, "", WarningWidget::WarningMsgInfo );
    }//if( m )
  }//for( int i = 0; i < 3; ++i )
}//void storeSpectraInDb()

//Menu: Spectrum->Save As
void SpecMeasManager::startStoreSpectraAsInDb()
{
  for( int i = 2; i >= 0; i-- )
  {
    const SpectrumType type = SpectrumType(i);
    std::shared_ptr<SpecMeas> meas = m_viewer->measurment( type );
    if( !meas )
      continue;
    
    const string typestr = descriptionText(type);
    AuxWindow *window = new AuxWindow( "Save " + typestr + " Spectrum As" ,
                                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                                      | AuxWindowProperties::TabletModal
                                      | AuxWindowProperties::SetCloseable) );
    window->rejectWhenEscapePressed();
    window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
    WGridLayout *layout = window->stretcher();
    WLineEdit *edit = new WLineEdit();
    edit->setEmptyText( "(Name to store under)" );
    
    std::shared_ptr<SpectraFileHeader> header
                               = m_fileModel->fileHeader( meas );
    const WString filename = header ? header->displayName() : WString();
    if( !filename.empty() )
      edit->setText( filename );
    
    WText *label = new WText( "Name" );
    layout->addWidget( label, 0, 0 );
    layout->addWidget( edit,  0, 1 );
  
    WTextArea *summary = new WTextArea();
    label = new WText( "Desc." );
    summary->setEmptyText( "(Optional description)" );
    layout->addWidget( label, 1, 0 );
    layout->addWidget( summary, 1, 1 );
    layout->setColumnStretch( 1, 1 );
    layout->setRowStretch( 1, 1 );
  
    Dbo::ptr<UserFileInDb> dbentry = m_fileModel->dbEntry( meas );
    if( dbentry )
    {
      edit->setText( dbentry->filename + " copy" );
      summary->setText( dbentry->description );
    }//if( dbentry )
    
  

    WPushButton *save = new WPushButton( "Save",window->footer() );
    save->setFloatSide(Right);
    save->setIcon( "InterSpec_resources/images/disk2.png" );
    
    WPushButton *cancel = new WPushButton( "Cancel" , window->footer());
    cancel->setFloatSide(Right);
    cancel->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );

    
    save->clicked().connect(
                        boost::bind( &SpecMeasManager::finishStoreAsSpectrumInDb,
                            this, edit, summary, meas, window ) );

    
//    if( m_viewer->toolTabsVisible() )
//    {
//      window->resize( 450, 250 );
//    }else
//    {
//      const int maxHeight = static_cast<int>(0.95*m_viewer->paintedHeight());
//      const int maxWidth = static_cast<int>(0.95*m_viewer->paintedWidth());
//      window->resize( std::min(450,maxWidth), std::min(250,maxHeight) );
//    }//if( m_viewer->toolTabsVisible() )
    
    window->centerWindow();
    window->show();
  }//for( int i = 2; i <= 0; i-- )
}//void startStoreSpectraAsInDb()


void SpecMeasManager::finishSaveSnapshotInDb(
                      const std::vector< std::shared_ptr<SpecMeas> > specs,
                      const std::vector< Wt::Dbo::ptr<UserFileInDb> > dbs,
                      const std::vector< Wt::WLineEdit * > edits,
                      const std::vector< Wt::WCheckBox * > cbs,
                      AuxWindow *window
                      )
{
  assert( specs.size() == dbs.size() );
  assert( specs.size() == edits.size() );
  
  vector<WString> descrips;
  for( const Wt::WLineEdit *edit : edits )
    descrips.push_back( edit->text() );
  delete window;
  
  for( size_t i = 0; i < specs.size(); ++i )
  {
    if( i < cbs.size() && cbs[i] && !cbs[i]->isChecked() )
      continue;
    
    DataBaseUtils::DbTransaction transaction( *m_sql );
    try
    {
      UserFileInDb *info = new UserFileInDb();
      *info = *dbs[i];
      info->snapshotParent = dbs[i];
      info->description = descrips[i].toUTF8();
      info->serializeTime = WDateTime::currentDateTime();
  
      Dbo::ptr<UserFileInDb> dbptr = m_sql->session()->add( info );
      UserFileInDbData *data = new UserFileInDbData();
      data->fileInfo = dbptr;
      try
      {
        data->setFileData( specs[i],
                           UserFileInDbData::sm_defaultSerializationFormat );
      }catch( std::exception & )
      {
        delete data;
        transaction.rollback();
      }//try / catch
    
      Dbo::ptr<UserFileInDbData> dataptr = m_sql->session()->add( data );
      
      UserFileInDb::makeWriteProtected( dbptr );
      
      transaction.commit();
      passMessage( "Saved snapshot of " + dbptr->filename,
                  "", WarningWidget::WarningMsgSave );
    }catch( FileToLargeForDbException &e )
    {
      transaction.rollback();
      passMessage( e.what(), "", WarningWidget::WarningMsgHigh );
    }catch( std::exception &e )
    {
      transaction.rollback();
      cerr << "Error saving snaphshot to database: " << e.what() << endl;
      passMessage( "Error saving snaphshot to database, sorry :(",
                  "", WarningWidget::WarningMsgHigh );
    }//try / catch
  }//for( size_t i = 0; i < specs.size(); ++i )
}//void finishSaveSnapshotInDb(...)


//If the name of the tag is given, then automatically apply,
//otherwise, silently save as a new tag
void SpecMeasManager::storeSpectraSnapshotInDb( const std::string tagname )
{
  vector<WText *> labels;
  vector<Wt::WLineEdit *> edits;
  std::vector< Wt::Dbo::ptr<UserFileInDb> > dbs;
  std::vector< std::shared_ptr<SpecMeas> > specs;
  
  for( int i = 0; i < 3; ++i )
  {
    const SpectrumType type = SpectrumType( i );
    std::shared_ptr<SpecMeas> m = m_viewer->measurment( type );
    std::shared_ptr<SpectraFileHeader> header= m_fileModel->fileHeader( m );
    Dbo::ptr<UserFileInDb> dbentry = m_fileModel->dbEntry( m );
    
    if( !header )
      continue;
    
    if( !dbentry )
    {
      SpectraFileHeader::saveToDatabase( m, header );
      dbentry = header->dbEntry();
//      WString msg = "Stored '";
//      msg += (header ? header->displayName() : WString("")) + "'";
//      passMessage( msg, "", WarningWidget::WarningMsgInfo );
    }//if( !dbentry )
      
    if( !dbentry )
    {
      passMessage( "Couldnt save spectrum to database in order to make a"
                    " snapshot", "", WarningWidget::WarningMsgInfo );
      continue;
    }//if( !dbentry )
    
    edits.push_back( new WLineEdit() );
    labels.push_back( new WText( descriptionText( SpectrumType(i) ) ) );
    dbs.push_back( dbentry );
    specs.push_back( m );
    
  }//for( int i = 0; i < 3; ++i )
  
  if( specs.empty() )
    return;
  
  vector<Wt::WCheckBox *> cbs;
    
  if (tagname.length()==0)
  {
      
      AuxWindow *window = new AuxWindow( "Save Spectrum Version",
                          (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
                           | AuxWindowProperties::IsAlwaysModal) );
      window->setWidth( 250 );
      
      WGridLayout *layout = window->stretcher();
      layout->addWidget( new WText("<b>Spectrum</b>", XHTMLText), 0, 0 );
      layout->addWidget( new WText("<b>Description</b>", XHTMLText), 0, 1 );
      if( specs.size() > 2 )
          layout->addWidget( new WText("<b>Save Snapshot</b>", XHTMLText), 0, 2 );
      
      for( int i = 0; i < static_cast<int>(specs.size()); ++i )
      {
          layout->addWidget( labels[i], i+1, 0, AlignLeft );
          layout->addWidget( edits[i], i+1, 1 );
          
          if( specs.size() > 2 )
          {
              WCheckBox *cb = new WCheckBox();
              cb->setChecked( true );
              layout->addWidget( cb, i+1, 2, AlignLeft );
              cbs.push_back( cb );
          }//if( specs.size() > 2 )
      }//for( int i = 0; i < specs.size(); ++i )
      
      layout->setColumnStretch( 1, 1 );
      
      WPushButton *store = new WPushButton( "Save" , window->footer());
      store->setIcon( "InterSpec_resources/images/disk2.png" );
      store->setFloatSide(Right);
      store->clicked().connect( boost::bind( &SpecMeasManager::finishSaveSnapshotInDb,
                                            this, specs, dbs, edits, cbs, window ) );
      
      WPushButton *cancel = new WPushButton( "Cancel" , window->footer());
      cancel->setFloatSide(Right);
      cancel->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
      window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
      
      if( edits.size() == 1 )
          edits[0]->enterPressed().connect(
                                           boost::bind( &SpecMeasManager::finishSaveSnapshotInDb,
                                                       this, specs, dbs, edits, cbs, window ) );
      
      window->disableCollapse();
      window->rejectWhenEscapePressed();
      window->centerWindow();
      window->show();
  } // tagname.length()==0
  else
  {
        //Given name, save version
        edits.clear();
        for( int i = 0; i < static_cast<int>(specs.size()); ++i )
        {
           //just give it the name
           edits.push_back(new WLineEdit(WString(tagname)));
        }//for( int i = 0; i < specs.size(); ++i )
        finishSaveSnapshotInDb(specs, dbs, edits, cbs, NULL);
      
      WString msg = "Created a new tag '";
      msg += tagname;
      msg += "' for '";
      msg += "'";
      passMessage( msg, "", WarningWidget::WarningMsgInfo );
      
  } // name.length>1
}//void storeSpectraSnapshotInDb()

#endif //#if( USE_DB_TO_STORE_SPECTRA )


