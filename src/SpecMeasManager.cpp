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
#include <boost/scope_exit.hpp>
#include <boost/system/error_code.hpp>

#include <Wt/WLink>
#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WTable>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WString>
#include <Wt/WServer>
#include <Wt/WBorder>
#include <Wt/WServer>
#include <Wt/WTextArea>
#include <Wt/WIconPair>
#include <Wt/WTabWidget>
#include <Wt/WTableCell>
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

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "SpecUtils/DateTime.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/ZipArchive.h"
#include "InterSpec/HelpSystem.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/RowStretchTreeView.h"


#if( USE_DB_TO_STORE_SPECTRA )
#include "InterSpec/DbFileBrowser.h"
#endif
  
#if( !ANDROID && !IOS )
#include "InterSpec/FileDragUploadResource.h"
#endif

using namespace Wt;
using namespace std;


#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif


const int ForeBtnInd = 0;//static_cast<int>(SpecUtils::SpectrumType::Foreground);
const int BackBtnInd = 1;//static_cast<int>(SpecUtils::SpectrumType::Background);
const int SecondBtnInd = 2;//static_cast<int>(SpecUtils::SpectrumType::SecondForeground);

/*
int buttonIndex( const SpecUtils::SpectrumType type )
{
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground: return 0;
    case SpecUtils::SpectrumType::Background: return 1;
    case SpecUtils::SpectrumType::SecondForeground: return 2;
  }
  return -1;
}//buttonIndex
*/

using SpecUtils::SpectrumType;
using SpecUtils::SaveSpectrumAsType;
int toint( const SpectrumType type ){ return static_cast<int>(type); }
int toint( const SaveSpectrumAsType type ){ return static_cast<int>(type); }
SpectrumType typeFromInt( int id ){ return SpectrumType(id); }



namespace
{
void setHeadersDbEntry( std::shared_ptr<SpectraFileHeader> header, Wt::Dbo::ptr<UserFileInDb> entry )
{
  header->setDbEntry( entry );
  wApp->triggerUpdate();
}

  class FileUploadDialog : public AuxWindow
  {
    //Class used to upload spectrum files.  Specialization of AuxWindow
    //  needed to manage InterSpec::displayedSpectrumChanged() connection.
    
    Wt::Signals::connection m_specChangedConection;
    WFileUpload *m_fileUpload;
    SpecMeasManager *m_manager;
    SpecUtils::SpectrumType m_type;
    
  public:
    FileUploadDialog( InterSpec *viewer,
                      SpecMeasManager *manager )
    : AuxWindow( "", (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                       | AuxWindowProperties::PhoneNotFullScreen | AuxWindowProperties::DisableCollapse) ),
      m_fileUpload( 0 ),
      m_manager( manager ),
      m_type( SpectrumType::Foreground )
    {
      setWindowTitle( "Select File To Open" );
      
      const bool noForeground = !viewer->measurment( SpectrumType::Foreground );
      
      string instructions;
      if( !viewer->isPhone() )
      {
        if( noForeground )
        {
          instructions = string(viewer->isMobile() ?"Tap" : "Click")
                         + " below to choose a foreground spectrum file to open.";
        }else
        {
          instructions = "Select how you would like to open the spectrum file, and then ";
          instructions += (viewer->isMobile() ? "tap" : "click");
          instructions += " below to browse for the file.";
        }//if( noForeground )
      }//if( !viewer->isPhone() )
      
      auto layout = stretcher();
      WText *txt = nullptr;
      
      if( !instructions.empty() )
      {
        txt = new WText( instructions );
        layout->addWidget( txt, layout->rowCount(), 0 );
      }
      
      if( !noForeground )
      {
        WGroupBox *buttons = new WGroupBox( "Open file as:" );
      
        layout->addWidget( buttons, layout->rowCount(), 0 );
        
        WButtonGroup *group = new WButtonGroup( buttons );
        
        WRadioButton *foreground = new WRadioButton( "Foreground", buttons );
        foreground->setInline( false );
        foreground->setChecked( true );
        group->addButton( foreground, toint(SpectrumType::Foreground) );
        
        WRadioButton *background = new WRadioButton( "Background", buttons );
        background->setInline( false );
        group->addButton( background, toint(SpectrumType::Background) );
        
        WRadioButton *secondary = new WRadioButton( "Secondary", buttons );
        secondary->setInline( false );
        group->addButton( secondary, toint(SpectrumType::SecondForeground) );
        
        group->checkedChanged().connect( std::bind( [this,group](){
          m_type = typeFromInt( group->checkedId() );
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
        msg = "<span style=\"font-size: 10px;\">"
                "You can also open spectrum files from email or file apps."
              "</span>";
      }

      WText *friendlyMsg = new WText( msg );
      friendlyMsg->setStyleClass( "startQuickUploadMsg" );
      
      layout->addWidget( friendlyMsg, layout->rowCount(), 0, AlignMiddle | AlignBottom );
      
      
      WPushButton *cancel = new WPushButton( "Cancel", footer() );
      
      cancel->clicked().connect( this, &FileUploadDialog::userCanceled );
      m_fileUpload->changed().connect( m_fileUpload, &Wt::WFileUpload::upload );
      m_fileUpload->uploaded().connect( this, &FileUploadDialog::finishUpload );
      m_fileUpload->fileTooLarge().connect( boost::bind( &FileUploadDialog::toLarge, this,
                                                        boost::placeholders::_1 ) );
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
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                | AuxWindowProperties::DisableCollapse
                | AuxWindowProperties::EnableResize) ),
  m_manager( manager )
  {
    WGridLayout *layout = new Wt::WGridLayout();
    contents()->setLayout(layout);
    
    WText *uploadText = new WText( "Foreground: " );
    WFileUpload *m_fileUpload = new WFileUpload(  );
    m_fileUpload->changed().connect( m_fileUpload, &Wt::WFileUpload::upload );
    m_fileUpload->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager, m_fileUpload, SpectrumType::Foreground));
    m_fileUpload->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge,
                                                      boost::placeholders::_1 ) );

    WText *uploadText2 = new WText( "Background: " );
    WFileUpload *m_fileUpload2 = new WFileUpload(  );
    m_fileUpload2->changed().connect( m_fileUpload2, &Wt::WFileUpload::upload );
    m_fileUpload2->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager, m_fileUpload2, SpectrumType::Background));
    m_fileUpload2->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge,
                                                       boost::placeholders::_1 ) );
    
    WText *uploadText3 = new WText( "Secondary Foreground: " );
    WFileUpload *m_fileUpload3 = new WFileUpload(  );
    m_fileUpload3->changed().connect( m_fileUpload3, &Wt::WFileUpload::upload );
    m_fileUpload3->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager, m_fileUpload3, SpectrumType::SecondForeground));
    m_fileUpload3->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge,
                                                       boost::placeholders::_1 ) );
    
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
  
protected:
  SpecMeasManager  *m_manager;
};//class UploadBrowser


SpecMeasManager::SpecMeasManager( InterSpec *viewer )
  : WObject(),
    m_treeView( NULL ),
    m_fileModel( NULL ),
    m_fileUpload( NULL ),
    m_viewer( viewer ),
    m_setButton ( NULL),
    m_setAsForeground ( NULL),
    m_setAsBackground ( NULL),
    m_setAsSecForeground ( NULL),
    m_combineToNewFileButton ( NULL),
    m_subsetOfMeasToNewFileButton( nullptr ),
    m_sumSpectraButton( nullptr ),
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
  
  for( SpecUtils::SaveSpectrumAsType type = SpecUtils::SaveSpectrumAsType(0);
      type < SpecUtils::SaveSpectrumAsType::NumTypes;
      type = SpecUtils::SaveSpectrumAsType(static_cast<int>(type)+1) )
  {
    m_downloadResources[static_cast<int>(type)] = new DownloadSpectrumResource(type, this, this);
  }
  
//
    m_treeView = new RowStretchTreeView();
    m_fileModel = new SpectraFileModel( m_treeView );
    m_treeView->setModel( m_fileModel );
    m_treeView->selectionChanged().connect( boost::bind( &SpecMeasManager::selectionChanged, this ) );

    startSpectrumManager(); //initializes
    deleteSpectrumManager(); //deletes instance
    
  m_sql = viewer->sql();

  for( SpecUtils::SaveSpectrumAsType i = SpecUtils::SaveSpectrumAsType(0);
      i < SpecUtils::SaveSpectrumAsType::NumTypes;
       i = SpecUtils::SaveSpectrumAsType(static_cast<int>(i)+1) )
    m_specificResources[static_cast<int>(i)] = (SpecificSpectrumResource *)0;
  
#if( !ANDROID && !IOS )
  m_foregroundDragNDrop->fileDrop().connect( boost::bind( &SpecMeasManager::handleFileDrop, this,
                                                         boost::placeholders::_1,
                                                         boost::placeholders::_2,
                                                         SpectrumType::Foreground ) );
  m_secondForegroundDragNDrop->fileDrop().connect( boost::bind( &SpecMeasManager::handleFileDrop,
                                                               this, boost::placeholders::_1,
                                                               boost::placeholders::_2,
                                                               SpectrumType::SecondForeground ) );
  m_backgroundDragNDrop->fileDrop().connect( boost::bind( &SpecMeasManager::handleFileDrop, this,
                                                         boost::placeholders::_1,
                                                         boost::placeholders::_2,
                                                         SpectrumType::Background ) );
#endif


}// SpecMeasManager

//Moved what use to be SpecMeasManager, out to a startSpectrumManager() to correct modal issues
void  SpecMeasManager::startSpectrumManager()
{
    const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );

    m_spectrumManagerWindow = new AuxWindow( "Spectrum Manager",
                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                     | AuxWindowProperties::TabletNotFullScreen) );
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
    //spec->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    
    Wt::WPushButton* uploadButton = new Wt::WPushButton("File...",uploadDiv);
    uploadButton->clicked().connect(  this, &SpecMeasManager::uploadSpectrum );
    HelpSystem::attachToolTipOn(uploadButton, "Import spectrum from file", showToolTips, HelpSystem::ToolTipPosition::Bottom );
    uploadButton->setIcon( "InterSpec_resources/images/file_search.png" );
    uploadButton->setMargin(10,Wt::Left);
    
#if( USE_DB_TO_STORE_SPECTRA )
    Wt::WPushButton* importButton = new Wt::WPushButton( "Previous...", uploadDiv );
    importButton->clicked().connect( boost::bind( &SpecMeasManager::browsePrevSpectraAndStatesDb, this ) );
    HelpSystem::attachToolTipOn(importButton, "Imports previously saved spectrum", showToolTips , HelpSystem::ToolTipPosition::Bottom);
    importButton->setIcon( "InterSpec_resources/images/db_small_white.png" );
    importButton->setMargin(2,Wt::Left);
    
#endif
    
    Wt::WContainerWidget *treeDiv   = createTreeViewDiv();
    Wt::WContainerWidget *buttonBar = createButtonBar();
    
    WContainerWidget * content = m_spectrumManagerWindow->contents();
  
    WGridLayout *layout = m_spectrumManagerWindow->stretcher();
    layout->addWidget(uploadDiv, 0,0);
    layout->addWidget( treeDiv,       1, 0 );
    layout->addWidget( buttonBar,     2, 0 );
    layout->setRowStretch( 1, 10 );
    
    layout->setVerticalSpacing( 0 );
    layout->setHorizontalSpacing( 0 );
    layout->setContentsMargins( 5, 5, 5, 5 );
    
    dynamic_cast<WWebWidget *>(layout->parent())->setHiddenKeepsGeometry(true);
    
    content->setOverflow(WContainerWidget::OverflowVisible,Orientation::Vertical); //necessary for menu to not be covered by footer
    //content->setStyleClass("filemanageroverflow");
    
    
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
FileDragUploadResource *SpecMeasManager::dragNDrop( SpecUtils::SpectrumType type )
{
  switch( type )
  {
    case SpectrumType::Foreground:
      return m_foregroundDragNDrop;
    case SpectrumType::SecondForeground:
      return m_secondForegroundDragNDrop;
    case SpectrumType::Background:
      return m_backgroundDragNDrop;
  }//switch( type )

  throw std::runtime_error( "Seriou problem in SpecMeasManager::dragNDrop(..)" );
  return NULL;
}//FileDragUploadResource *dragNDrop( SpecUtils::SpectrumType type )

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
    
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType( group->checkedId() );
    
    const string fileInZip = Wt::asString(index.data()).toUTF8();
    
    const string tmppath = SpecUtils::temp_dir();
    string tmpfile = SpecUtils::temp_file_name( "", tmppath );
    
    ifstream zipfilestrm( spoolName.c_str(), ios::in | ios::binary );
    
    ZipArchive::FilenameToZipHeaderMap headers
                                    = ZipArchive::open_zip_file( zipfilestrm );
    
    if( !headers.count(fileInZip) )
      throw runtime_error( "Couldnt find file in zip" );

    size_t nbytewritten = 0;
    
    {
#ifdef _WIN32
      const std::wstring wtmpfile = SpecUtils::convert_from_utf8_to_utf16(tmpfile);
      ofstream tmpfilestrm( wtmpfile.c_str(), ios::out | ios::binary );
#else
      ofstream tmpfilestrm( tmpfile.c_str(), ios::out | ios::binary );
#endif
      nbytewritten = read_file_from_zip( zipfilestrm, headers[fileInZip], tmpfilestrm );
    }
    
    handleFileDropWorker( fileInZip, tmpfile, type, nullptr, wApp );
    
    SpecUtils::remove_file( tmpfile );
  }catch( std::exception & )
  {
    passMessage( "Error extracting file from zip", "", 2 );
  }//try / catch
  
  delete window;
}//SpecMeasManager::extractAndOpenFromZip(...)


bool SpecMeasManager::handleZippedFile( const std::string &name,
                                        const std::string &spoolName,
                                        const SpecUtils::SpectrumType spectrum_type )
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
    
    const bool validtype = ((spectrum_type==SpectrumType::Foreground)
                             || (spectrum_type==SpectrumType::SecondForeground)
                             || (spectrum_type==SpectrumType::Background));

    
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
    group->addButton(button, toint(SpectrumType::Foreground) );
    button = new WRadioButton("Background", typecb);
    cblayout->addWidget( button, 0, 1, AlignCenter );
    group->addButton(button, toint(SpectrumType::Background) );
    button = new Wt::WRadioButton("Secondary", typecb);
    cblayout->addWidget( button, 0, 2, AlignCenter );
    group->addButton(button, toint(SpectrumType::SecondForeground) );
    
    if( validtype )
    {
      group->setCheckedButton( group->button(toint(spectrum_type)) );
      typecb->hide();
    }else if( !m_viewer->displayedHistogram(SpectrumType::Foreground) )
    {
      group->setCheckedButton( group->button(toint(SpectrumType::Foreground)) );
      typecb->hide();
    }else
    {
      group->setCheckedButton( group->button(toint(SpectrumType::Foreground)) );
    }
    
    AuxWindow *window = new AuxWindow( "Uploaded ZIP File Contents",
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                   | AuxWindowProperties::TabletNotFullScreen
                   | AuxWindowProperties::EnableResize) );
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
          && SpecUtils::iequals_ascii(n.substr(n.size()-4), ".zip" ) )
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
      const string tmppath = SpecUtils::temp_dir();
      string tmpfile = SpecUtils::temp_file_name( fileInZip, tmppath );
      
      ofstream tmpfilestrm( tmpfile.c_str(), ios::out | ios::binary );
      size_t nbytes = 0;
      
      try
      {
        nbytes = ZipArchive::read_file_from_zip( zipfilestrm,
                                        headers.begin()->second, tmpfilestrm );
        handleFileDropWorker( fileInZip, tmpfile.string<string>(), type, nullptr, nullptr );
      }catch( std::exception & )
      {
      }
      
      SpecUtils::remove_file( tmpfile.string<string>() );
      
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
    table->doubleClicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this,
                                                spoolName, group, table, window,
                                                boost::placeholders::_1 ) );
    window->footer()->addWidget( openButton );
    
    //openButton->clicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this, spoolName, type, selection, window ) );
    openButton->clicked().connect( boost::bind( &SpecMeasManager::extractAndOpenFromZip, this, spoolName, group, table, window, WModelIndex() ) );
    
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



template<size_t N, size_t M>
bool check_magic_number( const uint8_t (&magic_num)[N], const uint8_t (&data)[M] )
{
  return (0 == memcmp( magic_num, data, (N < M ? N : M)));
}


bool SpecMeasManager::handleNonSpectrumFile( const std::string &displayName,
                                             const std::string &fileLocation )
{
#ifdef _WIN32
  const std::wstring wpathstr = SpecUtils::convert_from_utf8_to_utf16(fileLocation);
  std::ifstream infile( wpathstr.c_str(), ios::in | ios::binary );
#else
  std::ifstream infile( fileLocation.c_str(), ios::in | ios::binary );
#endif
  
  if( !infile )
    return false;
 
  //get the filesize
  infile.seekg(0, ios::end);
  const size_t filesize = infile.tellg();
  infile.seekg(0);
  
  if( filesize <= 128 )
  {
    passMessage( "File was probably an empty file - not a spectrum file.", "", 2 );
    return true;
  }
  
  uint8_t data[1024] = { 0x0 };
  
  if( !infile.read( (char *)data, std::min(sizeof(data), filesize) ) )
  {
    passMessage( "Failed to read non-spectrum file.", "", 2 );
    return true;
  }
  infile.seekg(0);
  
  //Case insensitive search of 'term' in the header 'data'
  auto position_in_header = [&data]( const std::string &term ) -> int {
    const char * const char_start = (const char *)data;
    const char * const char_end = (const char *)(data + boost::size(data));
    const auto pos = std::search( char_start, char_end, begin(term), end(term),
                                 [](unsigned char a, unsigned char b) -> bool {
      return (a == b);
    } );
    if( pos == char_end )
      return -1;
    return static_cast<int>( pos - char_start );
  };//position_in_header lambda
  
  auto header_contains = [&data]( const std::string &term ) -> bool {
    const char * const char_start = (const char *)data;
    const char * const char_end = (const char *)(data + boost::size(data));
    const auto pos = std::search( char_start, char_end, begin(term), end(term),
                                 [](unsigned char a, unsigned char b) -> bool {
      return (std::tolower(a) == std::tolower(b));
    } );
    return (pos != char_end);
  };//header_contains lambda
  
  
  SimpleDialog *dialog = new SimpleDialog();
  WPushButton *closeButton = dialog->addButton( "Close" );
  WGridLayout *stretcher = new WGridLayout();
  stretcher->setContentsMargins( 0, 0, 0, 0 );
  dialog->contents()->setLayout( stretcher );
  dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowHidden );
  WText *title = new WText( "Not a spectrum file" );
  title->addStyleClass( "title" );
  stretcher->addWidget( title, 0, 0 );
  
  // const string filename = SpecUtils::filename(displayName);
  
  //Check if ICD2 file
  if( std::find( boost::begin(data), boost::end(data), uint8_t(60)) != boost::end(data) )
  {
    string datastr;
    datastr.resize( boost::size(data) + 1, '\0' );
    memcpy( &(datastr[0]), (const char *)&data[0], boost::size(data) );
    
    if( SpecUtils::icontains( datastr, "n42ns:")
       || SpecUtils::icontains( datastr, "AnalysisResults" )
       || SpecUtils::icontains( datastr, "AlarmInformation" )
       || SpecUtils::icontains( datastr, "DNDOARSchema" )
       || SpecUtils::icontains( datastr, "DNDOEWSchema" )
       || SpecUtils::icontains( datastr, "DNDOARSchema" ) )
    {
      WText *t = new WText( "This looks to be an N42 ICD2 file that contains analysis results rather than raw spectra.<br />"
                            "If you believe this to be a legitimate spectrum file, please email it to <a href=\"mailto:interspec@sandia.gov\" target=\"_blank\">interspec@sandia.gov</a> to support this file type." );
      stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
      t->setTextAlignment( Wt::AlignCenter );
      dialog->show();
      
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
  const bool issvg = header_contains( "<svg" );
  
  // Wt::Utils::guessImageMimeTypeData(<#const std::vector<unsigned char> &header#>)
  
  if( iszip )
  {
    const char *msg = NULL;
  //zip (but can be xlsx, pptx, docx, odp, jar, apk) 50 4B 03 04
    msg = "This file appears to be an invalid ZIP file (or xlsx, pptx, <br />"
          "docx, odp, jar, apk, etc), sorry I cant open it :(";
    WText *t = new WText( msg );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    
    return true;
  }//if( iszip )
  
  if( israr || istar || iszip7 || isgz )
  {
    const char *msg = "This file looks to be an archive file, not supported by InterSpec yet."
                      "Please contact "
                      "<a href=\"mailto:wcjohns@sandia.gov\" target=\"_blank\">wcjohns@sandia.gov</a> "
                      "if you would like support for this archive type added.";
    
    WText *t = new WText( msg );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    
    return true;
  }//if( israr || istar || iszip7 || isgz )
  
  
  if( ispdf | isps | istif )
  {
    const char *msg = "This file looks to be a document file, and not supported by InterSpec.";
    WText *t = new WText( msg );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    
    return true;
  }//if( ispdf | isps | istif )
  
  
  if( isgif || isjpg || ispng || isbmp || issvg )
  {
    const char *msg = "This file looks to be an image, and not a spectrum file.";
    WText *t = new WText( msg );
    t->setTextAlignment( Wt::AlignCenter );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
   
    const size_t max_disp_size = 16*1024*1024;
    
    if( filesize < max_disp_size )
    {
      const char *mimetype = "";
      if( isgif ) mimetype = "image/gif";
      else if( isjpg ) mimetype = "image/jpeg";
      else if( ispng ) mimetype = "image/png";
      else if( isbmp ) mimetype = "image/bmp";
      else if( issvg ) mimetype = "image/svg+xml";
    
      std::unique_ptr<WImage> image( new WImage() );
      WMemoryResource *resource = new WMemoryResource( mimetype, image.get() );
      vector<uint8_t> totaldata( filesize );
      const bool success = infile.read( (char *)&(totaldata[0]), filesize ).good();
      
      
      if( success )
      {
        resource->setData( totaldata );
        image->setImageLink( WLink(resource) );
        const int ww = m_viewer->renderedWidth();
        const int wh = m_viewer->renderedHeight();
        if( (ww > 120) && (wh > 120) )
          image->setMaximumSize( WLength(0.45*ww,WLength::Unit::Pixel), WLength(wh - 120, WLength::Unit::Pixel) );
        stretcher->addWidget( image.release(), stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
      }else
      {
        WText *errort = new WText( "Couldn't read uploaded file." );
        errort->setTextAlignment( Wt::AlignCenter );
        stretcher->addWidget( errort, stretcher->rowCount(), 0 );
      }
    }else
    {
      WText *errort = new WText( "Uploaded file was too large to try and display." );
      errort->setTextAlignment( Wt::AlignCenter );
      stretcher->addWidget( errort, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    }//if( filesize < max_disp_size ) / else
        
    return true;
  }//if( isgif || isjpg || ispng || isbmp )

  //Check if CSV giving peak ROIs.
  auto currdata = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  const bool possible_peak_csv = ( currdata
                                  //&& SpecUtils::icontains( SpecUtils::file_extension(displayName), "csv" )
                                  && header_contains("Centroid")
                                  && header_contains("Net_Area")
                                  && header_contains("FWHM") );
  
  const bool possible_gadras_peak_csv = ( currdata
                                  //&& SpecUtils::icontains( SpecUtils::file_extension(displayName), "csv" )
                                  && header_contains("Energy(keV)")
                                  && header_contains("Rate(cps)")
                                  && header_contains("FWHM(keV)")
                                  && header_contains("Centroid"));
  
  
  if( possible_peak_csv || possible_gadras_peak_csv )
  {
    try
    {
      const std::string seessionid = wApp->sessionId();
      const vector<PeakDef> orig_peaks = m_viewer->peakModel()->peakVec();
      
      if( possible_peak_csv )
      {
        const vector<PeakDef> candidate_peaks
                                      = PeakModel::csv_to_candidate_fit_peaks(currdata, infile);
        
        // For peaks from a InterSpec/PeakEasy CSV file, we will re-fit the peaks, as in practice
        //  they might not be from this exact spectrum file.
        Wt::WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
          PeakSearchGuiUtils::fit_template_peaks( m_viewer, currdata, candidate_peaks,
                                                 orig_peaks, PeakSearchGuiUtils::PeakTemplateFitSrc::CsvFile, seessionid );
        } ) );
      }else
      {
        assert( possible_gadras_peak_csv );
        
        const vector<PeakDef> candidate_peaks = PeakModel::gadras_peak_csv_to_peaks(currdata, infile);
        
        Wt::WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=](){
          PeakSearchGuiUtils::prepare_and_add_gadras_peaks( currdata, candidate_peaks,
                                                 orig_peaks, seessionid );
        } ) );
      }//if( possible_peak_csv )
      
      
      delete dialog;

      return true;
    }catch( exception &e )
    {
      WText *errort = new WText( "Uploaded file looked like a Peak CSV file, but was invalid." );
      errort->setTextAlignment( Wt::AlignCenter );
      stretcher->addWidget( errort, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
      
      errort = new WText( string(e.what()) );
      errort->setAttributeValue( "style", "color: red; font-weight: bold; font-family: monospace; " );
      errort->setTextAlignment( Wt::AlignCenter );
      stretcher->addWidget( errort, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
      
      return true;
    }//try / catch get candidate peaks )
  }//if( we could possible care about propagating peaks from a CSV file )
  
  
  // Check if this is an InterSpec exported DRF CSV, or XML file.
  const bool rel_eff_csv_drf = header_contains( "# Detector Response Function" );
  const int xml_drf_pos = position_in_header( "<DetectorPeakResponse" );
  
  if( rel_eff_csv_drf || ((xml_drf_pos >= 0) && (xml_drf_pos <= 20)) )
  {
    shared_ptr<DetectorPeakResponse> det;
    
    if( rel_eff_csv_drf )
    {
      det = DrfSelect::parseRelEffCsvFile( fileLocation );
    }else
    {
      try
      {
        if( filesize > 100*1024 ) //if larger than 100 KB, probably not a DRF
          throw runtime_error( "To large to be XML file" );
          
        rapidxml::file<char> input_file( infile );
        
        rapidxml::xml_document<char> doc;
        doc.parse<rapidxml::parse_default>( input_file.data() );
        auto *node = doc.first_node( "DetectorPeakResponse" );
        if( !node )
          throw runtime_error( "No DetectorPeakResponse node" );
        
        det = make_shared<DetectorPeakResponse>();
        det->fromXml( node );
      }catch( std::exception &e )
      {
        det.reset();
        log("info") << "Failed to parse perspective XML DRF file as DRF: " << e.what();
      }
    }//if( rel_eff_csv_drf ) / else XML DRF
    
    if( det && det->isValid() )
    {
      // TODO: generate a eff plot, and basic info, and display; probably by refactoring DrfSelect::updateChart()
      // TODO: Ask user if they want to use DRF; if so save to `InterSpec::writableDataDirectory() + "UploadedDrfs"`
      // TODO: handle GADRAS style Efficiency.csv files
      // TODO: allow users to rename the DRF.
      
      const string name = Wt::Utils::htmlEncode( det->name() );
      
      string msg = "<p style=\"white-space: nowrap;\">"
        "This file looks to be a Detector Response Function."
      "</p>"
      "<p style=\"text-align: left; white-space: nowrap;\">"
      "Name: " + name +
      "</p>"
      "<p>Would you like to use this DRF?</p>"
      ;
      
      DrfSelect::createChooseDrfDialog( {det}, msg, "" );
      
      delete dialog;
      
      return true;
    }
  }//if( maybe a drf )
  
  
  // Check if this is TSV/CSV file containing multiple DRFs
  if( header_contains( "Relative Eff" )
     && handleMultipleDrfCsv(infile, displayName, fileLocation) )
  {
    delete dialog;
    return true;
  }
  
  
  // Check if this is a PeakEasy CALp file
  if( currdata
     && header_contains( "CALp File" )
     && handleCALpFile(infile, dialog, false) )
  {
    return true;
  }
  
  delete dialog;
  
  return false;
}//void handleNonSpectrumFile(...)


bool SpecMeasManager::handleMultipleDrfCsv( std::istream &input,
                                           const std::string &displayName,
                                           const std::string &fileLocation )
{
  vector<string> credits;
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  
  DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );
  
  if( drfs.empty() )
    return false;
  
  vector<char> fileContents;
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  try
  {
    // We need to copy file contents into memory, because the file may disappear.
    SpecUtils::load_file_data( fileLocation.c_str(), fileContents );
  }catch( std::exception &e )
  {
    fileContents.clear();
    cerr << "handleMultipleDrfCsv: Failed to read spool file '" << fileLocation << "'" << endl;
  }
#endif
  
  std::function<void()> saveDrfFile;
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER )
  saveDrfFile = [displayName,fileContents](){
    if( fileContents.empty() )
      return;
    
    try
    {
      std::string datadir = InterSpec::writableDataDirectory();
      if( datadir.empty() )
        throw runtime_error( "Writable data directory not set." );
      
      datadir = SpecUtils::append_path( datadir, "drfs" );
      
      if( SpecUtils::create_directory(datadir) == 0 ) //-1 means already existed, 1 means created
        throw runtime_error( "Could not create 'drfs' directory in app data directory." );
      
      //displayName
      string filename = SpecUtils::filename( displayName );
      const string orig_extension = SpecUtils::file_extension( filename );
      assert( orig_extension.size() <= filename.size() );
      
      if( orig_extension.size() )
        filename = filename.substr( 0, filename.size() - orig_extension.size() );
      
      const int offset = wApp->environment().timeZoneOffset();
      const boost::posix_time::ptime now = WDateTime::currentDateTime().addSecs(60*offset).toPosixTime();
      string timestr = SpecUtils::to_vax_string(now); //"2014-Sep-19 14:12:01.62"
      const string::size_type pos = timestr.find( ' ' );
      //std::string timestr = SpecUtils::to_extended_iso_string( now ); //"2014-04-14T14:12:01.621543"
      //string::size_type pos = timestr.find( 'T' );
      if( pos != string::npos )
        timestr = timestr.substr(0,pos);
      SpecUtils::ireplace_all( timestr, "-", "_" );
      
      filename += "_" + timestr + orig_extension;
      const string outputname = SpecUtils::append_path( datadir, filename );
      
      
#ifdef _WIN32
      const std::wstring wtmpfile = SpecUtils::convert_from_utf8_to_utf16(outputname);
      ofstream outfilestrm( wtmpfile.c_str(), ios::out | ios::binary );
#else
      ofstream outfilestrm( outputname.c_str(), ios::out | ios::binary );
#endif
      
      if( !outfilestrm )
        throw runtime_error( "Unable to open file '" + outputname + "'" );
      
      if( !outfilestrm.write( &(fileContents[0]), fileContents.size() ) )
      {
        outfilestrm.close();
        SpecUtils::remove_file(outputname);
        
        throw runtime_error( "Failed writing '" + outputname + "'" );
      }//
      
      passMessage( "Saved '" + filename + "' for later use, and will be available in the"
                  " &quot;<em>Rel. Eff.</em>&quot; portion of the"
                  " &quot;<em>Detector Response Function Select</em>&quot; tool.",
                  "", WarningWidget::WarningMsgInfo );
    }catch( std::exception &e )
    {
      cerr << "handleMultipleDrfCsv: error saving multiple DRF file: " << e.what() << endl;
      passMessage( "Error saving DRF file for later use.", "", WarningWidget::WarningMsgHigh );
    }//try / catch to save file
  };//saveDrfFile
#endif
  
  
  string dialogmsg;
  if( drfs.size() == 1 )
    dialogmsg += "This file looks to be a Detector Response Function.";
  else
    dialogmsg += "This file contains multiple Detector Response Functions.";

  string creditsHtml;
  if( credits.size() )
  {
    for( string &s : credits )
      creditsHtml += "<div>" + s + "</div>";
  }//if( credits.size() )
  
  
  DrfSelect::createChooseDrfDialog( drfs, dialogmsg, creditsHtml, saveDrfFile );
  
  return true;
}//bool handleMultipleDrfCsv( std::istream &input, SimpleDialog *dialog )


bool SpecMeasManager::handleCALpFile( std::istream &infile, SimpleDialog *dialog, bool autoApply )
{
  WGridLayout *stretcher = nullptr;
  WPushButton *closeButton = nullptr;
  
  // Make a lamda to clear dialog, if we are going to alter it
  auto clear_dialog = [&](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    closeButton = dialog->addButton( "Close" );
    stretcher = new WGridLayout();
    
    // If we set the contents margins to 0, then scroll-bars may appear.
    //  However doing just the below looks okay, and the scroll bars dont seem to appear
    stretcher->setContentsMargins( 9, 2, 9, 2 );
    //dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowHidden );
    
    dialog->contents()->setLayout( stretcher );
    WText *title = new WText( "Not a spectrum file" );
    title->addStyleClass( "title" );
    stretcher->addWidget( title, 0, 0 );
  };//clear_dialog lamda
  
  
  auto currdata = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if( !currdata )
  {
    clear_dialog();
    assert( stretcher && closeButton );
    
    const string msg = "<p>No currently displayed foreground - skipping CALp file.</p>";
    WText *t = new WText( WString::fromUTF8(msg) );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    
    return true;
  }//if( !currdata )
  
  const size_t nchannel = currdata->num_gamma_channels();
  map<string,shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal;
  
  const std::streampos start_pos = infile.tellg();
  
  while( infile.good() )
  {
    string name;
    const auto cal = EnergyCal::energy_cal_from_CALp_file( infile, nchannel, name );
    
    if( !cal )
    {
      // Display message to user to let them know it was a CALp file, but we couldn't use it.
      //  TODO: improve this error message with details, ex if it was a lower-channel-energy CALp, and number of channels didnt match, we should display this to the user
      if( det_to_cal.empty() /* && SpecUtils::iends_with( displayName, "CALp" ) */ )
      {
        infile.seekg( start_pos, ios::beg );
        
        clear_dialog();
        assert( stretcher && closeButton );
        
        const string msg = "<p style=\"white-space: nowrap;\">"
        "The content of this CALp file was<br />"
        "unreadable or incompatible with the data"
        "</p>";
        WText *t = new WText( WString::fromUTF8(msg) );
        stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
        t->setTextAlignment( Wt::AlignCenter );
        dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowVisible,
                                         Wt::Horizontal | Wt::Vertical );
        return true;
      }//if( we didnt get any calibrations )
      
      break;
    }//if( !cal )
    
    if( !cal->valid() )
    {
      assert( 0 );
      continue;
    }
    
    det_to_cal[name] = cal;
  }//while( true )
  
  set<string> fore_gamma_dets;
  const shared_ptr<SpecMeas> foreground = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  const set<int> &fore_samples = m_viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
  const vector<string> fore_dets = m_viewer->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );
  
  if( det_to_cal.empty() || !foreground )
  {
    infile.seekg( start_pos, ios::beg );
    return false;
  }
  
  clear_dialog();
  assert( stretcher && closeButton );

  
  // We will grab calibrations from our current foreground, to slightly better inform the user
  //  TODO: do a similar thing for background and secondary spectra.
  
  set<shared_ptr<const SpecUtils::EnergyCalibration>> fore_energy_cals;
  //map<string,set<shared_ptr<const SpecUtils::EnergyCalibration>>> fore_meas_cals;
  
  for( const string &det_name : fore_dets )
  {
    for( const int sample : fore_samples )
    {
      const auto m = foreground->measurement( sample, det_name );
      if( m && m->num_gamma_channels() > 3 )
      {
        fore_gamma_dets.insert( m->detector_name() );
        fore_energy_cals.insert( m->energy_calibration() );
        //fore_meas_cals[det_name].insert( m->energy_calibration() );
      }
    }//for( const int sample : fore_samples )
  }//for( const string &det_name : fore_dets )
  
  bool have_cal_for_all_dets = true;
  for( const string &det_name : fore_gamma_dets )
  {
    if( !det_to_cal.count(det_name) )
      have_cal_for_all_dets = false;
  }
  
  // If we only have one calibration for our current data, or only one named detector, lets not
  //  care about matching names up exactly.
  if( (det_to_cal.size() == 1) && ((fore_energy_cals.size() == 1) || (fore_dets.size() == 1)) )
    have_cal_for_all_dets = true;
  
  string msg = "<p style=\"white-space: nowrap;\">"
  "This file looks to contain an energy calibration."
  "</p>";
  
  
  if( (det_to_cal.size() == 1) && det_to_cal.begin()->second )
  {
    shared_ptr<const SpecUtils::EnergyCalibration> cal = det_to_cal.begin()->second;
    assert( cal );
    const SpecUtils::EnergyCalType type = cal->type();
    
    switch( type )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      {
        msg += "<p>";
        if( type == SpecUtils::EnergyCalType::FullRangeFraction )
          msg += "FullRangeFrac:";
        else
          msg += "Polynomial:";
        
        for( size_t i = 0; i < cal->coefficients().size() && i < 4; ++i )
        {
          char buffer[64];
          snprintf( buffer, sizeof(buffer), "%s%.3f", (i ? ", " : " "), cal->coefficients()[i] );
          msg += buffer;
        }
        if( cal->coefficients().size() > 4 )
          msg += "...";
        
        if( cal->deviation_pairs().size() )
          msg += "<br />Plus " + std::to_string(cal->deviation_pairs().size()) + " deviation pairs";
        
        msg += "</p>";
        break;
      }//case polynomial or FRF
        
      case SpecUtils::EnergyCalType::LowerChannelEdge:
        msg += "<p>Lower channel energies.</p>";
        break;
      
      case SpecUtils::EnergyCalType::InvalidEquationType:
        break;
    }//switch( type )
  }else
  {
    msg += "<p>There are calibrations for " + std::to_string(det_to_cal.size()) + " detectors</p>";
  }
  
  if( !have_cal_for_all_dets )
  {
    if( det_to_cal.size() == 1 )
    {
      msg += "<p style=\"text-align: left; white-space: nowrap;\">"
      "Warning: This CALp file contains a single calibration; it will<br />"
      "be applied to the primary displayed energy calibration, and<br />"
      "the relative change will be applied to the other detectors<br />"
      "energy calibrations, which is probably what you want."
      "</p>";
    }else
    {
      msg += "<p style=\"text-align: left; white-space: nowrap;\">"
      "Warning: This energy calibration may not be consistent<br />"
      "with the structure of the current data."
      "</p>";
    }
  }//if( !have_cal_for_all_dets )
  
  
  WText *t = new WText( WString::fromUTF8(msg) );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  WCheckBox *applyOnlyCurrentlyVisible = nullptr;
  WCheckBox *applyForeground = nullptr, *applyBackground = nullptr, *applySecondary = nullptr;
  
  if( fore_samples.size() != foreground->sample_numbers().size() )
  {
    applyOnlyCurrentlyVisible = new WCheckBox( "Apply only to displayed samples" );
    stretcher->addWidget( applyOnlyCurrentlyVisible, stretcher->rowCount(), 0, AlignLeft );
  }//if( not displaying all foreground samples )
  
  const auto back = m_viewer->measurment( SpecUtils::SpectrumType::Background );
  const auto second = m_viewer->measurment( SpecUtils::SpectrumType::SecondForeground );
  
  //Only have
  if( (back && (back != foreground)) || (second && (second != foreground)) )
  {
    applyForeground = new WCheckBox( "Apply to foreground" );
    applyForeground->setChecked( true );
    stretcher->addWidget( applyForeground, stretcher->rowCount(), 0, AlignLeft );
  }
  
  if( back && (back != foreground) )
  {
    applyBackground = new WCheckBox( "Apply to background" );
    applyBackground->setChecked( true );
    stretcher->addWidget( applyBackground, stretcher->rowCount(), 0, AlignLeft );
  }
  
  if( second && (second != foreground) && (second != back) )
  {
    applySecondary = new WCheckBox( "Apply to secondary" );
    applySecondary->setChecked( true );
    stretcher->addWidget( applySecondary, stretcher->rowCount(), 0, AlignLeft );
  }
  
  
  msg = "<p style=\"white-space: nowrap;\">Would you like to use it?</p>";
  t = new WText( WString::fromUTF8(msg) );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  
  dialog->contents()->addStyleClass( "CALp" );
  // TODO: ask if they want to update deviation pairs - maybe?
  
  const auto applyLambda = [=](){
    InterSpec *interspec = InterSpec::instance();
    if( !interspec )
      return;
    
    
    int napplied = 0;
    
    try
    {
      const bool all_detectors = true; // TODO: give user option wether to apply to all detectors or not
      const bool all_samples = (!applyOnlyCurrentlyVisible
                                || !applyOnlyCurrentlyVisible->isChecked());
      
      EnergyCalTool *caltool = interspec->energyCalTool();
      assert( caltool ); //should always be valid
      if( !caltool )
        throw runtime_error( "Invalid EnergyCalTool" );
      
      if( !applyForeground || applyForeground->isChecked() )
      {
        caltool->applyCALpEnergyCal( det_to_cal, SpecUtils::SpectrumType::Foreground, all_detectors, all_samples );
        napplied += 1;
      }
      
      if( applyBackground && applyBackground->isChecked() )
      {
        caltool->applyCALpEnergyCal( det_to_cal, SpecUtils::SpectrumType::Background, all_detectors, all_samples );
        napplied += 1;
      }
      
      if( applySecondary && applySecondary->isChecked() )
      {
        caltool->applyCALpEnergyCal( det_to_cal, SpecUtils::SpectrumType::SecondForeground, all_detectors, all_samples );
        napplied += 1;
      }
    }catch( std::exception &e )
    {
      string msg = "There was an error applying calibration:<br />";
      msg += e.what();
      if( napplied == 1 )
        msg += "<br /><br />Calibration was applied to foreground.";
      else if( napplied == 2 )
        msg += "<br /><br />Calibration was applied to foreground and background.";
      
      passMessage( msg, "", WarningWidget::WarningMsgHigh );
    }//try / catch
  };//applyLambda(...)
  
  
  if( autoApply && !applyOnlyCurrentlyVisible && !applyForeground && !applyBackground && !applySecondary )
  {
    applyLambda();
    dialog->done( Wt::WDialog::DialogCode::Accepted );
    return true;
  }
  
  dialog->addButton( "No" ); //no further action necessary if user clicks no; dialog will close
  closeButton->setText( "Yes" );
  closeButton->clicked().connect( std::bind( applyLambda ) );
  
  return true;
}//void handleCALpFile( std::istream &input )



void SpecMeasManager::handleFileDropWorker( const std::string &name,
                     const std::string &spoolName,
                     SpecUtils::SpectrumType type,
                     SimpleDialog *dialog,
                     Wt::WApplication *app )
{
  // We are outside of the application loop here - we could parse the spectrum file here, instead
  //  of during the loop - but this
  
  if( !app )
    app = WApplication::instance();
  
  WApplication::UpdateLock lock( app );
 
  if( app && !lock )
  {
    cerr << "\n\nFailed to get WApplication::UpdateLock in "
            "SpecMeasManager::handleFileDropWorker(...) - this really shouldnt happen.\n" << endl;
    return;
  }//if( !lock )
 
  assert( WApplication::instance() );
  
  // Make sure we trigger a app update
  BOOST_SCOPE_EXIT(app,dialog){
    
    // TODO: there is a bit of a delay between upload completing, and showing the dialog - should check into that
    // TODO: check that the dialog is actually deleted correctly in all cases.
    if( dialog )
    {
      auto accept = boost::bind(&SimpleDialog::accept, dialog);
      WServer::instance()->post( wApp->sessionId(), std::bind([accept](){
        accept();
        WApplication::instance()->triggerUpdate();
      }) );
      dialog = nullptr;
    }//if( dialog )
    
    WApplication::instance()->triggerUpdate();
  } BOOST_SCOPE_EXIT_END
  
 
  if( name.length() > 4
     && SpecUtils::iequals_ascii( name.substr(name.length()-4), ".zip")
     && handleZippedFile( name, spoolName, type ) )
  {
    return;
  }
  
  
  try
  {
    std::shared_ptr<SpecMeas> measurement;
    std::shared_ptr<SpectraFileHeader> header;
    
    const int modelRow = setFile( name, spoolName, header, measurement );

    displayFile( modelRow, measurement, type, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
    
    //It is the responsibility of the caller to clean up the file.
  }catch( exception &e )
  {
    if( !handleNonSpectrumFile( name, spoolName ) )
    {
      displayInvalidFileMsg( name, e.what() );
    }
  }
  
}//handleFileDropWorker(...)


void SpecMeasManager::handleFileDrop( const std::string &name,
                                             const std::string &spoolName,
                                             SpecUtils::SpectrumType type )
{
  // If file is small, and not csv/txt (these are really slow to parse), dont display the parsing
  //  message.
  if( (SpecUtils::file_size(spoolName) < 512*1024)
     && !SpecUtils::iends_with(name, ".csv") && !SpecUtils::iends_with(name, ".txt") )
  {
    handleFileDropWorker( name, spoolName, type, nullptr, wApp );
    return;
  }
  
  // Its a larger file - display a message letting the user know its being parsed.
  auto dialog = new SimpleDialog( "Parsing File", "This may take a second." );
  
  wApp->triggerUpdate();
  
// When using WServer::instance()->post(...) it seems the "Parsing File" isnt always shown, but
//  posting to the ioService and explicitly taking the WApplication::UpdateLock seems to work a
//  little more reliable - I didnt look into why this is, or how true it is
//  WServer::instance()->post( wApp->sessionId(),
//                             boost::bind( &SpecMeasManager::handleFileDropWorker, this,
//                                          name, spoolName, type, dialog, wApp ) );
  
  WServer::instance()->ioService().boost::asio::io_service::post( boost::bind( &SpecMeasManager::handleFileDropWorker, this,
                                                    name, spoolName, type, dialog, wApp ) );
}//handleFileDrop(...)

void SpecMeasManager::displayInvalidFileMsg( std::string filename, std::string errormsg )
{
  //make sure we dont display the whole path
  string lastpart = SpecUtils::filename(filename);
  if( lastpart.empty() )
    lastpart = filename;
  if( lastpart.size() > 16 )
    lastpart = lastpart.substr(0,9) + "...";
  
  if( errormsg.empty() )
    errormsg = "Unspecified";
  
  lastpart = Wt::Utils::htmlEncode(lastpart);
  errormsg = Wt::Utils::htmlEncode(errormsg);
  
  string msg =
  "<p>Sorry, I couldnt parse the file " + lastpart + ":</p>"
  "<p>Error: <em>" + errormsg + "</em></p>"
  "<p>If you think this is a valid spectrum file, please send it to "
  "<a href=\"mailto:interspec@sandia.gov\" target=\"_blank\">interspec@sandia.gov</a>, and"
  " we'll try to fix this issue.</p>";
  
  //passMessage( msg.str(), "", 2 );
  
  SimpleDialog *dialog = new SimpleDialog( "Could Not Parse File", msg );
  dialog->addButton( "Okay" );
  wApp->triggerUpdate();
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


void SpecMeasManager::unDisplay( SpecUtils::SpectrumType type )
{
  m_viewer->setSpectrum( nullptr, {}, type, 0 );
  selectionChanged(); // update buttons
} // void SpecMeasManager::unDisplay( SpecUtils::SpectrumType type );


// The std::shared_ptr<SpecMeas> dummy_ptr keeps the SpecUtils::SpecFile
// object in memory when SpectraFileHeader isnt caching the spectrum, so
// its weak_ptr<> can be used in the call to header->parseFile();
void SpecMeasManager::loadSelected( const SpecUtils::SpectrumType type,
                                    std::shared_ptr<SpecMeas> dummy_ptr,
                                    const bool doPreviousEnergyRangeCheck )
{
  dummy_ptr = dummy_ptr; //keep compiler from complaining or optimizing dummy_ptr away
  loadSelected( type, doPreviousEnergyRangeCheck );
} // void SpecMeasManager::loadSelected(...)


void SpecMeasManager::loadSelected( const SpecUtils::SpectrumType type,
                                    const bool doPreviousEnergyRangeCheck )
{
  std::shared_ptr<SpectraFileHeader> header = selectedFile();
  std::shared_ptr<SpecMeas> meas = header ? header->parseFile() : nullptr;
  const set<int> displaySampleNums = selectedSampleNumbers();

  WFlags<InterSpec::SetSpectrumOptions> options;
  if( doPreviousEnergyRangeCheck )
  {
    options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
    options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
  }
  
  m_viewer->setSpectrum( meas, displaySampleNums, type, options );

  addToTempSpectrumInfoCache( meas );

  selectionChanged(); // update buttons
} // void SpecMeasManager::loadSelected(...)


void SpecMeasManager::startQuickUpload()
{
  new FileUploadDialog( m_viewer, this );
}//void startQuickUpload( SpecUtils::SpectrumType type )


void SpecMeasManager::finishQuickUpload( Wt::WFileUpload *upload,
                                         const SpecUtils::SpectrumType type )
{
  // TODO: The warning messages, and error conditions detected should be greatly improved

  std::shared_ptr<SpecMeas> measement_ptr;
  const int row = dataUploaded( upload, measement_ptr );
  
  if( row < 0 )
  {
    cerr << "SpecMeasManager::finishQuickUpload(...)\n\tError uploading file "
         << upload->clientFileName() << endl;
    return;
  } // if( row < 0 )

  displayFile( row, measement_ptr, type, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
}//void finishQuickUpload(...)


void SpecMeasManager::selectEnergyBinning( const string binning,
                                  std::shared_ptr<SpectraFileHeader> header,
                                  std::shared_ptr<SpecMeas> meas,
                                  const SpecUtils::SpectrumType type,
                                  const bool checkIfPreviouslyOpened,
                                  const bool doPreviousEnergyRangeCheck )
{
  WModelIndex index = m_fileModel->index( header );
  
  if( !index.isValid() )
  {
    passMessage( "Aborting loading of file after selecting energy binning - "
                "the file is no longer available in memory. please report "
                "this bug to interspec@sandia.gov", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !index.isValid() )
  
  
  if( binning != "Keep All" )
  {
    try
    {
      meas->keep_energy_cal_variant( binning );
    }catch( std::exception &e )
    {
      passMessage( "There was an error separating energy cal type; loading all (error: "
                   + std::string(e.what()) + ")", "", WarningWidget::WarningMsgHigh )
    }
    
    // Trigger a refresh of row info and selected rows in File Manager
    m_fileModel->removeRows( index.row(), 1 );
    header->setMeasurmentInfo( meas );
    m_fileModel->addRow( header );
    index = m_fileModel->index( header );
  }//if( binning != "Keep All" ) / else
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
              SpecMeasManager::VariantChecksToDo::None );
}//


void SpecMeasManager::selectDerivedDataChoice( const SpecMeasManager::DerivedDataToKeep tokeep,
                             std::shared_ptr<SpectraFileHeader> header,
                             std::shared_ptr<SpecMeas> meas,
                             const SpecUtils::SpectrumType type,
                             const bool checkIfPreviouslyOpened,
                             const bool doPreviousEnergyRangeCheck )
{
  WModelIndex index = m_fileModel->index( header );
  
  if( !index.isValid() )
  {
    passMessage( "Aborting loading of file after selecting derived data type to keep - "
                "the file is no longer available in memory. please report "
                "this bug to interspec@sandia.gov", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !index.isValid() )
  
  VariantChecksToDo furtherChecks = VariantChecksToDo::None;
  switch( tokeep )
  {
    case DerivedDataToKeep::All:
      furtherChecks = VariantChecksToDo::MultipleEnergyCal;
      break;
      
    case DerivedDataToKeep::RawOnly:
      furtherChecks = VariantChecksToDo::MultipleEnergyCal;
      break;
      
    case DerivedDataToKeep::DerivedOnly:
      furtherChecks = VariantChecksToDo::None;
      break;
  }//switch( tokeep )
  
  
  switch( tokeep )
  {
    case DerivedDataToKeep::All:
      break;
      
    case DerivedDataToKeep::RawOnly:
    case DerivedDataToKeep::DerivedOnly:
    {
      const auto keeptype = (tokeep==DerivedDataToKeep::DerivedOnly)
                            ? SpecUtils::SpecFile::DerivedVariantToKeep::Derived
                            : SpecUtils::SpecFile::DerivedVariantToKeep::NonDerived;
      
      meas->keep_derived_data_variant( keeptype );
      
      // Trigger a refresh of row info and selected rows in File Manager
      m_fileModel->removeRows( index.row(), 1 );
      header->setMeasurmentInfo( meas );
      m_fileModel->addRow( header );
      index = m_fileModel->index( header );
      break;
    }
  }//switch( tokeep )
  
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
              furtherChecks );
}//void selectDerivedDataChoice(...)



bool SpecMeasManager::checkForAndPromptUserForDisplayOptions( std::shared_ptr<SpectraFileHeader> header,
                                            std::shared_ptr<SpecMeas> meas,
                                            const SpecUtils::SpectrumType type,
                                            const bool checkIfPreviouslyOpened,
                                            const bool doPreviousEnergyRangeCheck,
                                            const VariantChecksToDo viewingChecks )
{
  if( !header || !meas )
    throw runtime_error( "SpecMeasManager::checkForAndPromptUserForDisplayOptions(): Invalid input" );
  
  if( viewingChecks == VariantChecksToDo::None )
    return false;
  
  bool derivedData = false, energyCal = false;
  if( viewingChecks == VariantChecksToDo::DerivedDataAndEnergy )
    derivedData = (meas->contains_derived_data() && meas->contains_non_derived_data());
  
  set<string> cals;
  
  if( !derivedData )
  {
    cals = meas->energy_cal_variants();
    energyCal = (cals.size() > 1);
  }
  
  if( !derivedData && !energyCal )
    return false;

  
  if( derivedData )
  {
    const char *title = "Use Derived Data?";
    const char *msgtxt
    = "<div>This file contained &quot;Derived Data&quot; spectra, which are</div>"
      "<div>usually what is used by the detection system for analysis.</div>"
      "<div style=\"margin-top: 20px; margin-bottom: 5px; white-space: nowrap; font-weight: bold;\">"
        "What data would you like to load?"
      "</div>";
    
    SimpleDialog *dialog = new SimpleDialog( title, msgtxt );
    WPushButton *button = dialog->addButton( "All Data" );
    
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::All, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    button = dialog->addButton( "Raw Data" );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::RawOnly, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    button = dialog->addButton( "Derived Data" );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::DerivedOnly, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    return true;
  }//if( derivedData )

  const char *title = "Select Binning";
  const char *msgtxt
   = "<div style=\"white-space: nowrap;\">Multiple energy binnings were found in the spectrum file.</div>"
     "<div style=\"margin-top: 5px; margin-bottom: 5px;\">Please select which one you would like</div>";
  
  if( cals.size() > 3 )
  {
    SimpleDialog *dialog = new SimpleDialog( title );
    
    int ncolwide = (derivedData ? 3 : static_cast<int>(cals.size() + 1));
    if( ncolwide > 4 )
      ncolwide = 4;
  
    WTable *table = new WTable( dialog->contents() );
  
    WText *msg = new WText( msgtxt , XHTMLText );
    //layout->addWidget( msg, 0, 0, 1, ncolwide );
    WTableCell *cell = table->elementAt( 0, 0 );
    cell->addWidget( msg );
    cell->setColumnSpan( ncolwide );
  
    WPushButton *button = new WPushButton( "Keep All" );
    cell = table->elementAt( 1, 0 );
    cell->addWidget( button );
    button->setWidth( WLength(95.0,WLength::Percentage) );
  
    button->clicked().connect( boost::bind( &SimpleDialog::reject, dialog ) );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this,
                                           string("Keep All"), header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    int calnum = 1;
    for( set<string>::const_iterator iter = cals.begin(); iter != cals.end(); ++iter, ++calnum )
    {
      const int row = 1 + (calnum / ncolwide);
      const int col = (calnum % ncolwide);
      
      //Make sure the calbration ID isnt too long
      string label = *iter;
      SpecUtils::utf8_limit_str_size( label, 15 );
      
      button = new WPushButton( label );
      //layout->addWidget( button, row, col );
      cell = table->elementAt( row, col );
      cell->addWidget( button );
      button->setWidth( WLength(95.0,WLength::Percentage) );
      
      button->clicked().connect( boost::bind( &SimpleDialog::reject, dialog ) );
      //Note, using *iter instead of label below.
      button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this,
                                             *iter, header, meas, type,
                                             checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    }//for( loop over calibrations )
  }else  //if( cals.size() > 3 )
  {
    SimpleDialog *dialog = new SimpleDialog( title );
    
    WText *msg = new WText( msgtxt , XHTMLText, dialog->contents() );
    msg->addStyleClass( "content" );
    
    WPushButton *button = dialog->addButton( "Keep All" );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this,
                                           string("Keep All"), header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    for( set<string>::const_iterator iter = cals.begin(); iter != cals.end(); ++iter )
    {
      //Make sure the calbration ID isnt too long
      string label = *iter;
      SpecUtils::utf8_limit_str_size( label, 15 );
      
      button = dialog->addButton( label );
      button->clicked().connect( boost::bind( &SpecMeasManager::selectEnergyBinning, this,
                                             *iter, header, meas, type,
                                             checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    }//for( loop over calibrations )
  }// if( cals.size() > 3 ) / else
  
  
  return true;
}//checkForAndPromptUserForDisplayOptions(...)



void SpecMeasManager::displayFile( int row,
                                   std::shared_ptr<SpecMeas> measement_ptr,
                                   const SpecUtils::SpectrumType type,
                                   bool checkIfPreviouslyOpened,
                                   const bool doPreviousEnergyRangeCheck,
                                   const SpecMeasManager::VariantChecksToDo viewingChecks )
{
  std::shared_ptr<SpecMeas> old_meas = m_viewer->measurment( type );
  std::shared_ptr<SpecMeas>  old_back;
  if( type == SpecUtils::SpectrumType::Foreground )
    old_back = m_viewer->measurment( SpecUtils::SpectrumType::Background );

#if( USE_DB_TO_STORE_SPECTRA )
  const bool storeInDb
      = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
#endif
  
  
  if( row < 0 && !measement_ptr )
  {
    WFlags<InterSpec::SetSpectrumOptions> options;
    options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
    options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
    
    m_viewer->setSpectrum( measement_ptr, std::set<int>(), type, options );

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

  if( !header || (header->measurementIfInMemory() != measement_ptr) )
  {
    const char *msg = "SpecMeasManager::displayFile(...): you must "
                      "pass in the SpectraFileModel row corresponding to the "
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
  
  
  if( viewingChecks != VariantChecksToDo::None )
  {
    if( checkForAndPromptUserForDisplayOptions( header, measement_ptr,
                                type, checkIfPreviouslyOpened,
                                doPreviousEnergyRangeCheck, viewingChecks ) )
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
    case SpecUtils::SpectrumType::Foreground:
    case SpecUtils::SpectrumType::SecondForeground:
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
            case SpecUtils::SourceType::IntrinsicActivity:
              //lets try to not show intrinsic activity by default
              ++numIntrinsic;
              intrinsicrow = static_cast<int>( i );
            break;
              
            case SpecUtils::SourceType::Foreground:
              ++numForeground;
              childrow = static_cast<int>( i );
              //i = header->m_samples.size();  //first foreground, terminate loop
            break;
              
            case SpecUtils::SourceType::Background:
              ++numbackground;
              backrow = static_cast<int>( i );
              if( numForeground )
                i = header->m_samples.size();  //We have foreground and background, terminate loop
            break;
            
            case SpecUtils::SourceType::Calibration:
              //do nothing
            break;
            
            case SpecUtils::SourceType::Unknown:
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
                     << " samples<br />";
          
          if( numForeground )
            warningmsg << "The first foreground sample is being shown.<br />";
          else if( childrow == 0 )
            warningmsg << "The first sample is being shown.<br />";
          else
            warningmsg << "Sample " << header->m_samples[childrow].sample_number << " is being shown.<br />";
          
          if( m_viewer->toolTabsVisible() )
            warningmsg << "Use the <b>Spectrum Files</b> tab, or the "
                          "<b>File Manager</b> to select other records";
          else
            warningmsg << "Use the <b>File Manager</b> to select others.";
        }//if( decide how to customize the info message ) / else

        //If we have an unambiguos background, and arent currently displaying a background
        //  from this detector, lets load the unambiguos background
        if( type==SpecUtils::SpectrumType::Foreground
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
        const bool hasDerived = measement_ptr ? measement_ptr->contains_derived_data() : false;
        const bool hasNonDerived = measement_ptr ? measement_ptr->contains_non_derived_data() : true;
        
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
          const bool back = (spectra.spectra_type == SpecUtils::SourceType::Background);
          const bool calib = (spectra.spectra_type == SpecUtils::SourceType::Calibration);
          const bool unWantedDerived = (hasDerived && hasNonDerived && spectra.is_derived_data);

          ncalibration += calib;

          if( back && (type != SpecUtils::SpectrumType::SecondForeground) && !unWantedDerived )
            background_sample_numbers.insert( spectra.sample_number );

          if( back || calib || unWantedDerived )
            selected.erase( index.child(sample,0) );
        }//for( int sample = 0; sample < nsamples; ++sample )

        if( ncalibration > 0 )
        {
          selected.erase( index );

          warningmsgLevel = max( WarningWidget::WarningMsgInfo, warningmsgLevel );
          warningmsg << "The uploaded file contained " << ncalibration
                       << " calibration or derived-data samples<br />";


          if( ncalibration == 1 ) warningmsg << "It";
          else                    warningmsg << "They";
          warningmsg  << " will not be displayed.<br />"
                      << "Use the <b>Spectrum Files</b> tab or <b>File Manager</b> to change this.";
        } // if( ncalibration > 0 )
      }// if( passthrough )

      break;
    } // case SpecUtils::SpectrumType::Foreground: case SpecUtils::SpectrumType::SecondForeground:

    case SpecUtils::SpectrumType::Background:
    {
      int nspectra_header = static_cast<int>( header->m_samples.size() );

      if( nspectra_header > 1 )
      {
        bool foundBackground = false;
        for( int sample = 0; sample < nspectra_header; ++sample )
        {
          const SpectraHeader &spectra = header->m_samples[sample];
          foundBackground = (spectra.spectra_type == SpecUtils::SourceType::Background);
          if( foundBackground )
          {
            warningmsgLevel = max( WarningWidget::WarningMsgLow, warningmsgLevel );
            warningmsg << "File has mutliple spectra, so only the background"
                          " spectrum has been displayed.";
            selected.clear();
            selected.insert( index.child(sample,0) );
            break;
          }//if( foundBackground )
        } // for( int sample = 0; sample < nsamples; ++sample )

        if( !foundBackground )
        {
          warningmsgLevel = max(WarningWidget::WarningMsgHigh, warningmsgLevel );
          warningmsg << "The uploaded file contained " << nspectra_header << " samples so<br />"
                     << "am displaying the first one.<br />"
                     << "Use the <b>File Manager</b> to change this.";
          selected.clear();
          selected.insert( index.child(0,0) );
        }//if( !foundBackground )
      } // if( !passthrough && (nsamples > 1) )
      break;
    } // case SpecUtils::SpectrumType::Background:
  } // switch( type )

  //Check if we removed any check or background samples, if so do not include
  //  the file "index" in the selected samples
  if( ncandiadate_samples != selected.size() )
    selected.erase( index );
  
  if( measement_ptr
     && (type == SpecUtils::SpectrumType::Foreground
         || type ==SpecUtils::SpectrumType::SecondForeground)
     && selected.empty() )
  {
    // I think this can happen if its passthrough/search data that are all marked background
    //  We'll just find the first spectrum.
    const int nspectra_header = static_cast<int>( header->m_samples.size() );
    for( int sample = 0; selected.empty() && (sample < nspectra_header); ++sample )
    {
      const SpectraHeader &spectra = header->m_samples[sample];
      
      //Note: we actually are just taking the first sample with non-zero gamma counts, which may
      //      not be a spectrum.  WE could improve this...
      if( spectra.gamma_counts_ > FLT_EPSILON )
        selected.insert( index.child(sample,0) );
    }//for( loop over to find a spectrum, any spectrum )
    
    background_sample_numbers.clear();
  }//if( selected.empty() )
  
  m_treeView->setSelectedIndexes( selected );

  if( warningmsg.str().size() )
    passMessage( warningmsg.str(), "SpecMeasManager", warningmsgLevel );
  
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
    WFlags<InterSpec::SetSpectrumOptions> options;
    if( doPreviousEnergyRangeCheck )
    {
      options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
      options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
    }
    
    m_viewer->setSpectrum( measement_ptr, background_sample_numbers,
                           SpecUtils::SpectrumType::Background, options );
  }//if( backgroundIndexs.size() )
  
#if( USE_DB_TO_STORE_SPECTRA )
  WApplication *app = wApp;
  if( checkIfPreviouslyOpened && !!header && app )
  {
    boost::function<void(void)> worker = app->bind(
                      boost::bind( &SpecMeasManager::checkIfPreviouslyOpened,
                                   this, app->sessionId(), header, type, m_destructMutex, m_destructed ) );
//    WServer::instance()->post( app->sessionId(), worker );
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
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
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
               | AuxWindowProperties::TabletNotFullScreen
               | AuxWindowProperties::DisableCollapse) );
  
  dialog->centerWindow();
  dialog->addStyleClass( "QuickSaveAsDialog" );

  std::shared_ptr<const SpecMeas> data, second, background, initial;
  data       = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  second     = m_viewer->measurment(SpecUtils::SpectrumType::SecondForeground);
  background = m_viewer->measurment(SpecUtils::SpectrumType::Background);

  initial = data;
  if( !initial )
    initial = second;
  if( !initial )
    initial = background;

  if( !initial )
  {
    const string msgstr = "There are no spectrums loaded into the viewer<br />"
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
  
  const vector<string> detnames = m_viewer->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  set<int> samplenums;
  if( data == initial )
    samplenums = m_viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
  
  for( SpecUtils::SaveSpectrumAsType i = SpecUtils::SaveSpectrumAsType(0);
      i < SpecUtils::SaveSpectrumAsType::NumTypes;
      i = SpecUtils::SaveSpectrumAsType(static_cast<int>(i)+1) )
  {
    SpecificSpectrumResource *&resource = m_specificResources[toint(i)];
    if( !resource )
    {
      resource = new SpecificSpectrumResource(i, this);
      resource->setTakesUpdateLock( true );
    }

    resource->setSpectrum( initial, samplenums, detnames );
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
      WRadioButton *button = new WRadioButton( "Foreground", buttons );
      button->setInline( false );
      button->setChecked( true );
      
      samplenums = m_viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
      
      for( SpecUtils::SaveSpectrumAsType i = SpecUtils::SaveSpectrumAsType(0);
           i < SpecUtils::SaveSpectrumAsType::NumTypes;
           i = SpecUtils::SaveSpectrumAsType(static_cast<int>(i)+1) )
      {
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[static_cast<int>(i)], data,
                                       samplenums, detnames ) );
      }
    }//if( data )
    
    if( second )
    {
      WRadioButton *button = new WRadioButton( "Second Foreground", buttons );
      button->setInline( false );
      
      samplenums = m_viewer->displayedSamples( SpectrumType::SecondForeground );
      for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
          i < SaveSpectrumAsType::NumTypes;
          i = SaveSpectrumAsType(static_cast<int>(i)+1) )
      {
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[toint(i)], second,
                                       samplenums, detnames ) );
      }
      
      if( !data )
        button->setChecked();
      group->addButton( button, -1 );
    }//if( second )

    if( background )
    {
      WRadioButton *button = new WRadioButton( "Background", buttons );
      button->setInline( false );
      
      samplenums = m_viewer->displayedSamples( SpecUtils::SpectrumType::Background );
      for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
          i < SaveSpectrumAsType::NumTypes;
          i = SaveSpectrumAsType(toint(i)+1) )
      {
        button->clicked().connect(
                          boost::bind( &SpecificSpectrumResource::setSpectrum,
                                       m_specificResources[toint(i)], background,
                                       samplenums, detnames ) );
      }
      
      group->addButton( button, -1 );
    } // if( background )
  } // if( nspecs > 1 )

  WContainerWidget *linkDiv = new WContainerWidget( dialog->contents() );
  linkDiv->setList( true, false );

  for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
      i < SaveSpectrumAsType::NumTypes;
      i = SaveSpectrumAsType(toint(i)+1) )
  {
    const string linktitle = string("As ") + descriptionText(i) + " File";
    SpecificSpectrumResource *resource = m_specificResources[toint(i)];
    WAnchor *a = new WAnchor( resource, linktitle, linkDiv );
    a->setTarget( TargetNewWindow );
    a->setInline( false );
    a->setStyleClass( "LoadSpectrumSaveAsLink" );
    
#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    a->clicked().connect( std::bind([resource,dialog](){
      android_download_workaround(resource, "spec_file_export");
      dialog->hide();
    }) );
#else
    a->clicked().connect( dialog, &AuxWindow::hide );
#endif
  }//for( SaveSpectrumAsType i = ... )
  

  WPushButton *cancel = new WPushButton( "Cancel", dialog->contents() );
  cancel->setInline( false );
  cancel->clicked().connect( dialog, &AuxWindow::hide );
  dialog->finished().connect( boost::bind( &SpecMeasManager::cleanupQuickSaveAsDialog, this, dialog, wApp ) );
  
  dialog->show();
} // void SpecMeasManager::displayQuickSaveAsDialog();


void SpecMeasManager::removeSpecMeas( std::shared_ptr<const SpecMeas> meas, const bool undisplay )
{
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( row );
    if( header->measurementIfInMemory() == meas )
    {
      std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
      
      if( undisplay )
      {
        //also unassign if it was assigned to foreground/2ndf/background
        if( meas == m_viewer->measurment(SpecUtils::SpectrumType::Foreground) )
          unDisplay(SpecUtils::SpectrumType::Foreground);
        if( meas== m_viewer->measurment(SpecUtils::SpectrumType::SecondForeground) )
          unDisplay(SpecUtils::SpectrumType::SecondForeground);
        if( meas == m_viewer->measurment(SpecUtils::SpectrumType::Background) )
          unDisplay(SpecUtils::SpectrumType::Background);
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
  //Goes through each loaded state and checks if it's the same one as
  //the removed UserFileInDb
  for( int row = 0; row < m_fileModel->rowCount(); ++row )
  {
    std::shared_ptr<SpectraFileHeader>  header = m_fileModel->fileHeader( row );
    if (header->dbEntry()==remove) {
      std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();

      
      //also unassign if it was assigned to foreground/2ndf/background
      if (meas==m_viewer->measurment(SpecUtils::SpectrumType::Foreground))
        unDisplay(SpecUtils::SpectrumType::Foreground);

      if (meas==m_viewer->measurment(SpecUtils::SpectrumType::SecondForeground))
        unDisplay(SpecUtils::SpectrumType::SecondForeground);
      
      if (meas==m_viewer->measurment(SpecUtils::SpectrumType::Background))
        unDisplay(SpecUtils::SpectrumType::Background);
        
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
        
    for( const WModelIndex &index : selected )
    {
        const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
        if( indexLevel == SpectraFileModel::FileHeaderLevel )
            selectedFiles.insert( index );
    } //for( const WModelIndex &index : selected )
    
    if( selectedFiles.size() < 1 )
    {
        cerr << "SpecMeasManager::removeSelected()\n\tThere are " << selectedFiles.size()
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
            if (meas==m_viewer->measurment(SpecUtils::SpectrumType::Foreground))
                unDisplay(SpecUtils::SpectrumType::Foreground);
            
            if (meas==m_viewer->measurment(SpecUtils::SpectrumType::SecondForeground))
                unDisplay(SpecUtils::SpectrumType::SecondForeground);
            
            if (meas==m_viewer->measurment(SpecUtils::SpectrumType::Background))
                unDisplay(SpecUtils::SpectrumType::Background);
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
      type < SaveSpectrumAsType::NumTypes;
      type = SaveSpectrumAsType(static_cast<int>(type)+1) )
  {
    const string name = origName + "." + suggestedNameEnding( type );
    m_downloadResources[toint(type)]->suggestFileName( name, WResource::Attachment );
  }
} // void SpecMeasManager::renameSaveAsFile()


std::shared_ptr<SpecMeas> SpecMeasManager::selectedToSpecMeas() const
{
  const WModelIndexSet selected = m_treeView->selectedIndexes();

  if( selected.empty() )
    throw runtime_error( "No files or samples are selected" );

  // First we will go through and figure out which sample numbers we want from each selected file.
  //
  // Note that currently the Spectrum Manager only shows records at the "sample" level, and
  //  not the {sample,detector} level.
  map<std::shared_ptr<const SpecMeas>,set<int>> files_involved;
  for( const WModelIndex &index : selected )
  {
    const SpectraFileModel::Level index_level = m_fileModel->level(index);
    
    shared_ptr<SpectraFileHeader> header;
    switch( index_level )
    {
      case SpectraFileModel::FileHeaderLevel:
        header = m_fileModel->fileHeader( index.row() );
        break;
        
      case SpectraFileModel::SampleLevel:
        header = m_fileModel->fileHeader( index.parent().row() );
        break;
        
      case SpectraFileModel::InvalidLevel:
        break;
    }//switch( index_level )
    
    if( !header )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      const string msg = "Failed to get header for row " + std::to_string( index.row() );
      log_developer_error( __func__, msg.c_str() );
#endif
      continue;
    }//if( !header )
    
    
    const shared_ptr<const SpecMeas> file = header->parseFile();
    if( !file )
    {
#if( PERFORM_DEVELOPER_CHECKS )
      const string msg = "Failed to parse file '" + header->displayName().toUTF8() + "'";
      log_developer_error( __func__, msg.c_str() );
#endif
      continue;
    }//if( !file )
    
    
    switch( index_level )
    {
      case SpectraFileModel::FileHeaderLevel:
      {
        set<int> &samples = files_involved[file];
        samples = file->sample_numbers();
        
        break;
      }//case SpectraFileModel::FileHeaderLevel
        
      
      case SpectraFileModel::SampleLevel:
      {
        assert( index.row() >= 0 );
        if( index.row() < 0 )
          continue;
        
        const size_t row = static_cast<size_t>( index.row() );
        
        assert( row < header->m_samples.size() );
        if( row >= header->m_samples.size() )
          continue;
        
        set<int> &samples = files_involved[file];
        samples.insert( header->m_samples[row].sample_number );
        break;
      }//case SpectraFileModel::SampleLevel:
        
      case SpectraFileModel::InvalidLevel:
        break;
    }//switch( indexLevel )
  }//for( const WModelIndex &index : selected )
  
  if( files_involved.empty() )
    throw runtime_error( "[logic error] Unable to determine contents to use." );
  
  const string current_time = WDateTime::currentDateTime().toString("yyyyMMdd hh:mm:ss").toUTF8();
  
  // Check if we want samples from a single spectrum file, and if so, we'll clone that SpecMeas
  //  object, and just removed any un-wanted samples.  This helps to preserve all the detector
  //  information and stuff.
  if( files_involved.size() == 1 )
  {
    const shared_ptr<const SpecMeas> &parent_file = files_involved.begin()->first;
    const set<int> &samples_to_use = files_involved.begin()->second;
    assert( parent_file );
    
    shared_ptr<SpecMeas> newspec = make_shared<SpecMeas>();
    newspec->uniqueCopyContents( *parent_file );
    newspec->set_uuid( "" );
    
    if( samples_to_use == parent_file->sample_numbers() )
      return newspec;
    
    vector<shared_ptr<const SpecUtils::Measurement>> meas_to_remove;
    for( const shared_ptr<const SpecUtils::Measurement> m : newspec->measurements() )
    {
      if( m && (samples_to_use.count( m->sample_number() ) == 0) )
        meas_to_remove.push_back( m );
    }
    
    for( const auto m : meas_to_remove )
      newspec->remove_measurement( m, false );
     
    unsigned int cleanup_flags = SpecUtils::SpecFile::CleanupAfterLoadFlags::StandardCleanup
                          | SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples;
    newspec->cleanup_after_load( cleanup_flags );
    
    if( !parent_file->filename().empty() )
      newspec->set_filename( "Subset " + current_time + " " + parent_file->filename() );
    
    newspec->cleanup_orphaned_info();
    newspec->setShieldingSourceModel( nullptr );
    newspec->displayedSpectrumChangedCallback( SpecUtils::SpectrumType::Foreground, newspec, {}, {} );
    
    return newspec;
  }//if( files_involved.size() == 1 )
  
  
  // If we are here, we are combining multiple file.
  //  We could probably to a little better preserving the detection system level "meta-info", but
  //  for the moment it isnt clear what the best thing to do is, in general.
  //  Also, we arent preserving any peaks, which maybe we could do.
  shared_ptr<SpecMeas> newspec = make_shared<SpecMeas>();
  
  
  // We'll create a map from newly created Measurements, to newly created peaks for them; we will
  //  then set the peaks to the SpecMeas object after cleanup (because sample numbers may change
  //  after cleanup, and its sample numbers who peaks belong to)
  map<set<shared_ptr<SpecUtils::Measurement>>,
      shared_ptr<deque<shared_ptr<const PeakDef>>> > peaks_to_keep;
  
  // We'll try to preserve some of the file-level meta information, if it makes sense.
  //  Although using set<string> like below will re-order things, which maybe isnt what we want.
  //  Note: meas->detectors_analysis() is currently not being saved.
  set<int> lane_numbers;
  set<SpecUtils::DetectorType> det_types;
  set<string> locations, inst_types, manufacturers, inst_models, inst_ids, remarks, warnings, insps;
  for( const auto specmeas_samples : files_involved )
  {
    const shared_ptr<const SpecMeas> &meas = specmeas_samples.first;
    const set<int> &samples = specmeas_samples.second;
  
    assert( meas );
    
    det_types.insert( meas->detector_type() );
    locations.insert( meas->measurement_location_name() );
    inst_types.insert( meas->instrument_type() );
    manufacturers.insert( meas->manufacturer() );
    inst_models.insert( meas->instrument_model() );
    inst_ids.insert( meas->instrument_id() );
    insps.insert( meas->inspection() );
    lane_numbers.insert( meas->lane_number() );
    
    const vector<string> &meas_remarks = meas->remarks();
    remarks.insert( begin(meas_remarks), end(meas_remarks) );
    
    const vector<string> &meas_warning = meas->parse_warnings();
    warnings.insert( begin(meas_warning), end(meas_warning) );
    
    
    map< shared_ptr<const SpecUtils::Measurement>, shared_ptr<SpecUtils::Measurement> > old_to_new_meas;
    for( const int sample : samples )
    {
      for( const string &det_name : meas->detector_names() )
      {
        const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det_name );
        
        if( m )
        {
          auto new_meas = make_shared<SpecUtils::Measurement>( *m );
          old_to_new_meas[m] = new_meas;
          newspec->add_measurement( new_meas,  false );
        }
      }//for( const string &det_name : meas->detector_names() )
    }//for( const int sample : samples )
    
    
    for( const set<int> &peakssamples : meas->sampleNumsWithPeaks() )
    {
      shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = meas->peaks(peakssamples);
      if( !peaks || peaks->empty() )
        continue;
      
      bool will_keep_all_meas = true;
      set<shared_ptr<SpecUtils::Measurement>> new_meas_set;
      
      for( const int sample : peakssamples )
      {
        for( const auto m : meas->sample_measurements(sample) )
        {
          if( old_to_new_meas.count(m) )
          {
            new_meas_set.insert( old_to_new_meas[m] );
          }else
          {
            will_keep_all_meas = false;
          }
        }//for( const auto m : meas->sample_measurements(sample) )
      }//for( const int sample : peakssamples )
      
      if( !will_keep_all_meas )
        continue;
      
      
      // We will create completely new instances of all the old peaks
      auto new_peaks = make_shared< deque<shared_ptr<const PeakDef>> >();
      
      // However, PeakDef has a design flaw that shared continuums arent necessarily constant, and
      //  multiple peaks may share a continuum; so we also have to take a little care to clone those
      //  and keep the shared relationship (maybe this has been worked around in some/most places,
      //  but we'll play it safe and so a deep clone here)
      map<shared_ptr<const PeakContinuum>,std::shared_ptr<PeakContinuum>> old_cont_map;
      
      for( shared_ptr<const PeakDef> peak : *peaks )
      {
        assert( peak );
        if( !peak )
          continue;
        
        auto p = make_shared<PeakDef>(*peak);
        shared_ptr<const PeakContinuum> old_cont = peak->continuum();
        if( old_cont_map.count(old_cont) )
        {
          p->setContinuum( old_cont_map[old_cont] );
        }else
        {
          p->makeUniqueNewContinuum();
          old_cont_map[old_cont] = p->continuum();
        }
        
        new_peaks->push_back( p );
      }//for( shared_ptr<const PeakDef> peak : *peaks )
      
      peaks_to_keep[new_meas_set] = new_peaks;
    }//for( const set<int> &peakssamples : meas->sampleNumsWithPeaks() )
  }//for( const auto specmeas_samples : files_involved )
  
  
  if( det_types.size() == 1 )
    newspec->set_detector_type( *begin(det_types) );
  
  if( lane_numbers.size() == 1 )
    newspec->set_lane_number( *begin(lane_numbers) );
  
  newspec->set_remarks( vector<string>( begin(remarks), end(remarks) ) );
  
  newspec->set_parse_warnings( vector<string>( begin(warnings), end(warnings) ) );
  
  const auto strings_to_csv = []( const set<string> &input ) -> string {
    string answer;
    for( const string &val : input )
    {
      if( !val.empty() )
        answer += (answer.empty() ? "" : ", ") + val;
    }
    return answer;
  };//strings_to_csv

  newspec->set_uuid( "" );
  newspec->set_measurement_location_name( strings_to_csv(locations) );
  newspec->set_instrument_id( strings_to_csv(inst_ids) );
  newspec->set_instrument_model( strings_to_csv(inst_models) );
  newspec->set_manufacturer( strings_to_csv(manufacturers) );
  newspec->set_instrument_type( strings_to_csv(inst_types) );
  newspec->set_inspection( strings_to_csv(insps) );
  
  newspec->set_filename( "Combination-" + current_time );
  
  newspec->cleanup_after_load();
  
  for( auto meas_to_peaks : peaks_to_keep )
  {
    const set<shared_ptr<SpecUtils::Measurement>> &meass = meas_to_peaks.first;
    const shared_ptr<deque<shared_ptr<const PeakDef>>> peaks = meas_to_peaks.second;
    
    assert( !meass.empty() );
    assert( peaks && !peaks->empty() );
    if( meass.empty() || !peaks || peaks->empty() ) // Shouldnt ever evaluate to tru, but Jic
    {
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Got an empty set of peaks for new measuremnt..." );
#endif
      continue;
    }
    
    set<int> sample_nums;
    for( auto m : meass )
    {
      assert( m );
      if( m )
        sample_nums.insert( m->sample_number() );
    }
    
    assert( !sample_nums.empty() );
    if( sample_nums.size() )
      newspec->setPeaks( *peaks, sample_nums );
  }//for( auto meas_to_peaks : peaks_to_keep )
  
  
  return newspec;
}//std::shared_ptr<SpecMeas> SpecMeasManager::selectedToSpecMeas()


void SpecMeasManager::newFileFromSelection()
{
  try
  {
    shared_ptr<SpecMeas> spec = selectedToSpecMeas();
    shared_ptr<SpectraFileHeader> header = addFile( spec->filename(), spec );
    
    WModelIndex index = m_fileModel->index( header );
    if( index.isValid() )
    {
      WModelIndexSet selected;
      selected.insert( index );
      m_treeView->setSelectedIndexes( selected );
      selectionChanged();
      
      m_treeView->scrollTo( index, WAbstractItemView::ScrollHint::EnsureVisible );
    }// if( index.isValid() )
  }catch( std::exception &e )
  {
    stringstream msg;
    msg << "Failed in combining files: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.str().c_str() );
#endif
  }
}//void SpecMeasManager::newFileFromSelection()


void SpecMeasManager::sumSelectedSpectra()
{
  try
  {
    shared_ptr<SpecMeas> meas = selectedToSpecMeas();
    assert( meas );
      
    shared_ptr<SpecUtils::Measurement> spec = meas->sum_measurements( meas->sample_numbers(), meas->detector_names(), nullptr );
    spec->set_sample_number( 1 );
    meas->remove_measurements( meas->measurements() );
    meas->add_measurement( spec, true );
    meas->set_filename( "summed" );
    
    shared_ptr<SpectraFileHeader> header = addFile( meas->filename(), meas );
    
    WModelIndex index = m_fileModel->index( header );
    if( index.isValid() )
    {
      WModelIndexSet selected;
      selected.insert( index );
      m_treeView->setSelectedIndexes( selected );
      selectionChanged();
      
      m_treeView->scrollTo( index, WAbstractItemView::ScrollHint::EnsureVisible );
    }// if( index.isValid() )
  }catch( std::exception &e )
  {
    const string msg = "Failed summing spectra: " + string( e.what() );
    passMessage( msg, "", WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
  }//try / catch
}//void SpecMeasManager::sumSelectedSpectra()



void SpecMeasManager::selectionChanged()
{
  if (!m_spectrumManagerWindow || m_spectrumManagerWindow->isHidden())
      return;
    
  const WModelIndexSet selected = m_treeView->selectedIndexes();
  WModelIndexSet toSelect = selected;
  WModelIndexSet selectedFiles;
  
  if( selected.empty() )
  {
    m_saveButton->hide();
    m_combineToNewFileButton->hide();
    m_subsetOfMeasToNewFileButton->hide();
    m_sumSpectraButton->hide();
    m_deleteButton->hide();
    m_setButton->hide();
  }else
  {
    m_deleteButton->show();
    
    set<shared_ptr<const SpectraFileHeader>> files;
    bool fullFileSelected = false;
    for( const WModelIndex &index : selected )
    {
      std::shared_ptr<const SpectraFileHeader> header;
      
      const SpectraFileModel::Level indexLevel = m_fileModel->level(index);
      
      switch( indexLevel )
      {
        case SpectraFileModel::FileHeaderLevel:
        {
          fullFileSelected = true;
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
          break;
      } // switch( level( index ) )
      
      if( header )
        files.insert( header );
    }// for( const WModelIndex &index : selected )
    
    
    if( files.size() > 1 )
    {
      m_combineToNewFileButton->show();
      m_subsetOfMeasToNewFileButton->hide();
      m_sumSpectraButton->show();
      m_setButton->hide();
      m_setAsForeground->disable();
      m_setAsBackground->disable();
      m_setAsSecForeground->disable();
    }else
    {
      m_setAsForeground->enable();
      
      if( m_viewer->measurment(SpecUtils::SpectrumType::Foreground) )
      {
        m_setAsBackground->enable();
        m_setAsSecForeground->enable();
      }else
      {
        m_setAsBackground->disable();
        m_setAsSecForeground->disable();
      }//if( primary spectrums loaded ) / else
      
      m_combineToNewFileButton->hide();
      m_subsetOfMeasToNewFileButton->setHidden( fullFileSelected );
      m_sumSpectraButton->show();
      m_setButton->show();
    }// if( multiple files are selected ) / else
    
    if( (selectedFiles.size()==1) && (files.size()==1) )
    {
      m_saveButton->show();
    }else
    {
      const bool enableSave = (files.size()==1);
      m_saveButton->setHidden( !enableSave );
    } // if( selectedFiles.size() == 1 )
    
    if( selected.size() != toSelect.size() )
      m_treeView->setSelectedIndexes( toSelect );
  }//if( selected.empty() ) / else
  
  // Disable/hide everything and just show what's needed.
  m_removeForeButton->hide();
  m_removeBackButton->hide();
  m_removeFore2Button->hide();
  
  if( m_viewer->measurment(SpecUtils::SpectrumType::Foreground) )
    m_removeForeButton->show();

  if( m_viewer->measurment(SpecUtils::SpectrumType::SecondForeground) )
    m_removeFore2Button->show();
    
  if( m_viewer->measurment(SpecUtils::SpectrumType::Background) )
    m_removeBackButton->show();
}//void selectionChanged()


void SpecMeasManager::deleteSpectrumManager()
{
    m_spectrumManagertreeDiv->removeWidget(m_treeView);
    AuxWindow::deleteAuxWindow(m_spectrumManagerWindow);
    m_spectrumManagerWindow = NULL;
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
  m_setAsForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::Foreground, true) );
  m_setAsBackground = setPopup->addItem("Background");
  m_setAsBackground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::Background, true) );
  m_setAsSecForeground = setPopup->addItem("Secondary Foreground");
  m_setAsSecForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::SecondForeground, true ) );
  m_setButton->setMenu(setPopup);
  m_setButton->hide();

  // Note: the only difference between the m_combineToNewFileButton and
  //  m_subsetOfMeasToNewFileButton buttons are the icons, and a slight difference in text, but the
  //  functionality they trigger is the same either way.
  m_combineToNewFileButton = new Wt::WPushButton( "To New File", m_newDiv );
  m_combineToNewFileButton->setToolTip( "Creates a new in-memory spectrum file from the current selection." );
  m_combineToNewFileButton->addStyleClass("InvertInDark");
  m_combineToNewFileButton->setIcon( "InterSpec_resources/images/arrow_join.svg" );
  m_combineToNewFileButton->clicked().connect( boost::bind( &SpecMeasManager::newFileFromSelection, this ) );
  m_combineToNewFileButton->hide();
  
  m_subsetOfMeasToNewFileButton = new Wt::WPushButton( "As New File", m_newDiv );
  m_subsetOfMeasToNewFileButton->setToolTip( "Creates a new in-memory spectrum file from the current selection." );
  m_subsetOfMeasToNewFileButton->addStyleClass("InvertInDark");
  m_subsetOfMeasToNewFileButton->setIcon( "InterSpec_resources/images/partial.svg" );
  m_subsetOfMeasToNewFileButton->clicked().connect( boost::bind( &SpecMeasManager::newFileFromSelection, this ) );
  m_subsetOfMeasToNewFileButton->hide();
  
  
  m_sumSpectraButton = new Wt::WPushButton( "Sum Spectra", m_newDiv );
  m_sumSpectraButton->setToolTip( "Creates a new in-memory spectrum file from the current selection." );
  m_sumSpectraButton->addStyleClass("InvertInDark");
  m_sumSpectraButton->setIcon( "InterSpec_resources/images/sum_symbol.svg" );
  m_sumSpectraButton->clicked().connect( boost::bind( &SpecMeasManager::sumSelectedSpectra, this ) );
  m_sumSpectraButton->hide();
  
  
  m_saveButton = new Wt::WPushButton("Export",m_newDiv);
  m_saveButton->setIcon( "InterSpec_resources/images/bullet_arrow_down.png" );
  Wt::WPopupMenu *setPopup2 = new Wt::WPopupMenu();
  
  m_saveButton->mouseWentOver().connect( boost::bind( &SpecMeasManager::renameSaveAsFile, this ) );
  
  for( SpecUtils::SaveSpectrumAsType type = SpecUtils::SaveSpectrumAsType(0);
      type < SpecUtils::SaveSpectrumAsType::NumTypes;
      type = SpecUtils::SaveSpectrumAsType(static_cast<int>(type)+1) )
  {
    const string descrip = string("As ") + descriptionText(type) + " File";
    WMenuItem *temp = setPopup2->addItem( descrip );
    DownloadSpectrumResource *resource = m_downloadResources[toint(type)];
    temp->setLink( resource );
    temp->setLinkTarget( TargetNewWindow );
    
#if( ANDROID )
    // Using hacked saving to temporary file in Android, instead of via network download of file.
    temp->clicked().connect( std::bind([resource](){
      android_download_workaround(resource, "spewc_download");
    }) );
#endif //ANDROID
  }//for( loop over save-as spectrum file types )
  
  m_saveButton->setMenu( setPopup2 );
  m_saveButton->hide();
  
  
  m_deleteButton = new Wt::WPushButton("Unload",m_newDiv);
  m_deleteButton->setIcon( "InterSpec_resources/images/minus_min_white.png" );
  m_deleteButton->clicked().connect( boost::bind( &SpecMeasManager::removeSelected, this ) );
  m_deleteButton->hide();
  
  // ---- try new bar 2 ----
  WContainerWidget *m_newDiv2 = new WContainerWidget( );
  buttonAlignment->addWidget( m_newDiv2, 2, 0);
  m_newDiv2->setStyleClass( "LoadSpectrumUploadDiv" );
  
  //WText *text =
  new WText( "Unassign: ", m_newDiv2 );
  
  m_removeForeButton = new Wt::WPushButton("Foreground",m_newDiv2);
//      m_removeForeButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeForeButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::Foreground) );
  m_removeForeButton  ->setHidden( true, WAnimation() );
//  m_removeForeButton->setHiddenKeepsGeometry( true );
  
  m_removeBackButton = new Wt::WPushButton("Background",m_newDiv2);
//          m_removeBackButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeBackButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::Background ) );
  m_removeBackButton  ->setHidden( true, WAnimation() );
//  m_removeBackButton->setHiddenKeepsGeometry( true );
  
  m_removeFore2Button = new Wt::WPushButton("Secondary Foreground",m_newDiv2);
//              m_removeFore2Button->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeFore2Button->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::SecondForeground ) );
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
    WServer::instance()->ioService().boost::asio::io_service::post( worker );
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
  WServer::instance()->ioService().boost::asio::io_service::post( f );
}//userCanceledResumeFromPreviousOpened(..)



void SpecMeasManager::showPreviousSpecFileUsesDialog( std::shared_ptr<SpectraFileHeader> header,
                                    const SpecUtils::SpectrumType type,
                                    const std::vector<Wt::Dbo::ptr<UserFileInDb>> &modifiedFiles,
                                    const std::vector<Wt::Dbo::ptr<UserFileInDb>> &unModifiedFiles,
                                    const std::vector<Wt::Dbo::ptr<UserState>> &userStatesWithFile )
{
  assert( header );
  if( !header )
    return;
  
  // The only place we call this function from has already taken care of the case where
  //  both modifiedFiles and userStatesWithFile are empty.
  assert( !modifiedFiles.empty() || !userStatesWithFile.empty() );
  
  if( modifiedFiles.empty() && unModifiedFiles.empty() && userStatesWithFile.empty() )
    return;

  if( unModifiedFiles.size() )
  {
    try
    {
      Dbo::ptr<UserFileInDb> f = unModifiedFiles.front();
      cerr << "Setting file with UUID=" << header->m_uuid << ", and filename "
      << header->m_displayName << " to be connected to DB entry from "
      << "upload at "
      << f->uploadTime.toString(DATE_TIME_FORMAT_STR).toUTF8()
      << " until user determines if they want to resume from a modified"
      << " session" << endl;
      header->setDbEntry( f );
    }catch( std::exception &e )
    {
      cerr << "Failed to set DB entry to SpectraFileHeader: " << e.what() << endl;
      
      assert( 0 );
    }//try / catch
  }//if( unModifiedFiles.size() )

  AuxWindow *window = new AuxWindow( "Previously Stored States",
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse)
                                     | AuxWindowProperties::EnableResize
                                     | AuxWindowProperties::TabletNotFullScreen) );
  window->rejectWhenEscapePressed();
  window->addStyleClass( "ShowPrevSpecFileUses" );
  
  //bool auto_save_states = false;
  SnapshotBrowser *snapshots = nullptr;
  AutosavedSpectrumBrowser *auto_saved = nullptr;
  
  try
  {
    //auto_save_states = InterSpecUser::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
    
    if( userStatesWithFile.size() )
    {
      // TODO: pass userStatesWithFile into SnapshotBrowser
      snapshots = new SnapshotBrowser( this, m_viewer, header, nullptr, nullptr );
      snapshots->finished().connect( boost::bind( &AuxWindow::deleteSelf, window) );
    }//if( userStatesWithFile.size() )
    
    
    if( !modifiedFiles.empty() )
    {
      auto_saved = new AutosavedSpectrumBrowser( modifiedFiles, type, m_fileModel, this, header );
      auto_saved->loadedASpectrum().connect( window, &AuxWindow::deleteSelf );
      
      WPushButton *cancel = window->addCloseButtonToFooter();
      cancel->clicked().connect( window, &AuxWindow::hide );
      window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
      
      if( unModifiedFiles.empty() )
      {
        window->finished().connect( boost::bind( &SpecMeasManager::userCanceledResumeFromPreviousOpened,
                                                this, nullptr, header ) );
      }else
      {
        window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
      }//if( unModifiedFiles.empty() ) / else
      
    }//if( !modifiedFiles.empty() )
  }catch( std::exception &e )
  {
    if( snapshots )
      delete snapshots;
    if( auto_saved )
      delete auto_saved;
    if( window )
      delete window;
    
    snapshots = nullptr;
    auto_saved = nullptr;
    window = nullptr;
    
    string msg = "Unexpected issue displaying available save-states: " + string( e.what() );
    passMessage( msg, "", WarningWidget::WarningMsgLevel::WarningMsgHigh)
    return;
  }// try / catch
    
  if( !snapshots && !auto_saved )
  {
    assert( 0 );
    delete window;
    return;
  }
  
  WGridLayout *layout = window->stretcher();
  if( snapshots && auto_saved )
  {
    WTabWidget *tabbed = new WTabWidget();
    layout->addWidget( tabbed, 0, 0 );
    
    tabbed->addTab( snapshots, "Saved States", WTabWidget::LoadPolicy::PreLoading );
    tabbed->addTab( auto_saved, "Auto-saved Spectra", WTabWidget::LoadPolicy::PreLoading );
  }else
  {
    if( snapshots )
      layout->addWidget( snapshots, 0, 0 );
    else
      layout->addWidget( auto_saved, 0, 0 );
  }//if( snapshots && auto_saved ) / else
  
  WCheckBox *cb = new WCheckBox( "Automatically store and check for prev. work." );
  cb->addStyleClass( "PrefCb" );
  layout->addWidget( cb, 1, 0 );
  layout->setRowStretch( 0, 1 );
  
  InterSpecUser::associateWidget( m_viewer->m_user, "AutoSaveSpectraToDb", cb, m_viewer );
  
  
  const int width = std::min( 500, static_cast<int>(0.95*m_viewer->renderedWidth()) );
  const int height = std::min( 475, static_cast<int>(0.95*m_viewer->renderedHeight()) );
  window->resize( WLength(width), WLength(height) );
  
  window->centerWindow();
  window->show();
  
  wApp->triggerUpdate();
}//void showPreviousSpecFileUsesDialog(..)


void SpecMeasManager::showDatabaseStatesHavingSpectrumFile( std::shared_ptr<SpectraFileHeader> header )
{
  const size_t num_states = SnapshotBrowser::num_saved_states( m_viewer, m_viewer->sql(), header );
  
  if( !num_states )
    return;
  
  DbFileBrowser *browser = new DbFileBrowser( this, m_viewer, header );
  
  // \TODO: come up with a way to check if there are any states before going through and creating
  //        the widget and everything.
  if( browser->numSnapshots() <= 0 )
  {
    assert( 0 );
    delete browser;
    browser = nullptr;
  }
  
  wApp->triggerUpdate();
}


void postErrorMessage( const string msg, const WarningWidget::WarningMsgLevel level )
{
  passMessage( msg, "", level );
  wApp->triggerUpdate();
}//void postErrorMessage( string msg )


void SpecMeasManager::checkIfPreviouslyOpened( const std::string sessionID,
                                 std::shared_ptr<SpectraFileHeader> header,
                                 SpecUtils::SpectrumType type,
                                 std::shared_ptr< std::mutex > mutex,
                                 std::shared_ptr<bool> destructed )
{
  std::lock_guard<std::mutex> lock( *mutex );
  
  if( *destructed )
  {
    cerr << "checkIfPreviouslyOpened(): manager destructed before I could do ish" << endl;
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
    
    vector<Dbo::ptr<UserState>> userStatesWithFile;
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
   
      
      // Now get user-saved states with this spectrum file in them
      Dbo::collection< Dbo::ptr<UserState> > states_query
                    = SnapshotBrowser::get_user_states_collection( m_viewer->m_user, sql, header );
      
      for( auto iter = states_query.begin(); iter != states_query.end(); ++iter )
        userStatesWithFile.push_back( *iter );
      
      m_viewer->m_user.modify()->incrementSpectraFileOpened();
      
      transaction.commit();
    }//end interaction with database
    
    
    if( modifiedFiles.empty() && unModifiedFiles.empty() && userStatesWithFile.empty() )
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
    }//if( this is a new-to-us file )

    
    if( modifiedFiles.empty() && userStatesWithFile.empty() )
    {
      Dbo::ptr<UserFileInDb> dbentry = unModifiedFiles.front();
      cerr << "Setting file with UUID=" << header->m_uuid << ", and filename "
           << header->m_displayName << " to be connected to DB entry from "
           << "upload at "
           << dbentry->uploadTime.toString(DATE_TIME_FORMAT_STR).toUTF8()
           << endl;
      
      WServer::instance()->post( sessionID, boost::bind( &setHeadersDbEntry, header, dbentry) );
      
      return;
    }//if( user has opened the file, but didnt modify or save it )
    
    
    // If we are here, the user has either modified the file, or has it as part of the save-state.
    //WServer::instance()->post( sessionID,
    //                            boost::bind( &SpecMeasManager::showDatabaseStatesHavingSpectrumFile,
    //                                        this, header) );
      
    
    WServer::instance()->post( sessionID,
                              boost::bind( &SpecMeasManager::showPreviousSpecFileUsesDialog,
                                          this, header, type, modifiedFiles, unModifiedFiles,
                                          userStatesWithFile ) );
    
  }catch( std::exception &e )
  {
    cerr << "Error checking if this file had been opened previously: " << e.what() << endl;
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
                              SpecUtils::ParserType parser_type )
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


int SpecMeasManager::dataUploaded2( Wt::WFileUpload *upload , SpecUtils::SpectrumType type)
{
  std::shared_ptr<SpecMeas> measurement;
  int row= dataUploaded( upload, measurement );
  displayFile( row, measurement, type, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
  return row;
} // int SpecMeasManager::dataUploaded( Wt::WFileUpload )


int SpecMeasManager::dataUploaded( Wt::WFileUpload *upload )
{
  std::shared_ptr<SpecMeas> measurement;
  return dataUploaded( upload, measurement );
} // int SpecMeasManager::dataUploaded( Wt::WFileUpload )


bool SpecMeasManager::loadFromFileSystem( const string &name, SpecUtils::SpectrumType type,
                                         SpecUtils::ParserType parseType )
{
  const string origName = SpecUtils::filename( name );
  
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

    displayFile( row, measurement, type, true, true, SpecMeasManager::VariantChecksToDo::DerivedDataAndEnergy );
  }catch( const std::exception &e )
  {
    {
      ifstream test( name.c_str(), ios::in | ios::binary );
      const bool iszip = (test.get()==0x50 && test.get()==0x4B
                          && test.get()==0x03 && test.get()==0x04);
      test.close();
      
      if( iszip
         /*&& SpecUtils::iequals_ascii( origName.substr(origName.length()-4), ".zip")*/
         && handleZippedFile( origName, name, type ) )
      return true;
    }
    
    if( !handleNonSpectrumFile( origName, name ) )
    {
      displayInvalidFileMsg(origName,e.what());
    }
    return false;
  }// try/catch
  
  return true;
}//void loadFromFileSystem( std::string filename )


int SpecMeasManager::dataUploaded( Wt::WFileUpload *upload, std::shared_ptr<SpecMeas> &measurement )
{
  const string fileName = upload->spoolFileName();
  const WString clientFileName = upload->clientFileName();
  const string origName = clientFileName.toUTF8();

  try
  {
    std::shared_ptr<SpectraFileHeader> header;
    int result = setFile( origName, fileName, header, measurement, SpecUtils::ParserType::Auto );
    WModelIndexSet selected;
    WModelIndex index = m_fileModel->index( result, 0 );
    selected.insert( index );
    m_treeView->setSelectedIndexes( WModelIndexSet() );

    //passMessage( "Successfully opened file.", "dataUploaded", 0 );

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
      WServer::instance()->ioService().boost::asio::io_service::post( worker );
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

  if( !lock )
  {
    //Probably wont ever happen, and if it does, we want to end things anyway.
    cerr << "\nSpecMeasManager::cleanupQuickSaveAsDialog(...)\n\tFailed to get app lock - not doing anything!\n\n" << endl;
    return;
  }
  
  
  //dialog->accept();
  delete dialog;
    
  for( SaveSpectrumAsType i = SaveSpectrumAsType(0);
      i < SaveSpectrumAsType::NumTypes;
      i = SaveSpectrumAsType(static_cast<int>(i)+1) )
  {
    m_specificResources[toint(i)]->setSpectrum( nullptr, {}, {} );
  }
    
  app->triggerUpdate();
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
void SpecMeasManager::browsePrevSpectraAndStatesDb()
{
  // TODO: Make this be the same implementation as SpecMeasManager::showPreviousSpecFileUsesDialog; but to do that, need to make AutosavedSpectrumBrowser be a MVC widget so we dont put like a million elements into the DOM
  new DbFileBrowser( this, m_viewer, nullptr );
}//void browsePrevSpectraAndStatesDb()




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
    cerr << "SpecMeasManager::finishStoreAsSpectrumInDb(...)\n\tFailed to save file to DB." << endl;
    passMessage( "Error saving to database", "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !header )
  
  header->setDbEntry( Dbo::ptr<UserFileInDb>() );
  
  SpectraFileHeader::saveToDatabase( meas, header );
  
  Wt::Dbo::ptr<UserFileInDb> entry = header->dbEntry();
  

  if( !entry )
  {
    cerr << "SpecMeasManager::finishStoreAsSpectrumInDb(...)\n\tFailed to save file to DB, apparently" << endl;
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
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType( i );
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
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(i);
    std::shared_ptr<SpecMeas> meas = m_viewer->measurment( type );
    if( !meas )
      continue;
    
    const string typestr = descriptionText(type);
    AuxWindow *window = new AuxWindow( "Save " + typestr + " Spectrum As" ,
                                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                                      | AuxWindowProperties::TabletNotFullScreen
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
      passMessage( "Saved state of " + dbptr->filename,
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
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType( i );
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
    labels.push_back( new WText( descriptionText( SpecUtils::SpectrumType(i) ) ) );
    dbs.push_back( dbentry );
    specs.push_back( m );
  }//for( int i = 0; i < 3; ++i )
  
  if( specs.empty() )
    return;
  
  vector<Wt::WCheckBox *> cbs;
    
  if (tagname.length()==0)
  {
      
      AuxWindow *window = new AuxWindow( "Save Spectrum Version",
                          (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
                           | AuxWindowProperties::IsModal) );
      window->setWidth( 250 );
      
      WGridLayout *layout = window->stretcher();
      layout->addWidget( new WText("<b>Spectrum</b>", XHTMLText), 0, 0 );
      layout->addWidget( new WText("<b>Description</b>", XHTMLText), 0, 1 );
      if( specs.size() > 2 )
          layout->addWidget( new WText("<b>Save State</b>", XHTMLText), 0, 2 );
      
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


