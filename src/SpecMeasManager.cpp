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
#include <boost/asio/deadline_timer.hpp>

#include <Wt/WLink>
#include <Wt/WText>
#include <Wt/Utils>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WTable>
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
#include <Wt/Chart/WAxis>

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
#include <Wt/WRegExpValidator>
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
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/UriSpectrum.h"


#include "InterSpec/MakeDrf.h"
#include "InterSpec/DrfChart.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/ZipArchive.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/MakeDrfSrcDef.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ExportSpecFile.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/LocalTimeDelegate.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/FileDragUploadResource.h"
#include "InterSpec/ShieldingSourceDisplay.h"

#if( USE_DB_TO_STORE_SPECTRA )
#include "InterSpec/DbFileBrowser.h"
#endif

#if( USE_REL_ACT_TOOL )
#include "InterSpec/RelActAutoGui.h"
#endif

#if( USE_QR_CODES )
#include "rapidxml/rapidxml_print.hpp"

#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/QrCode.h"
#include "InterSpec/QRSpectrum.h"
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
#if( USE_DB_TO_STORE_SPECTRA )
  class PreviousDbEntry : public WContainerWidget
  {
    //a class to display, and then select a previous database entry.  This class
    //  gets displayed in a AuxWindow when the user uploads a spectrum they
    //  have previously used and modified.  This class represents the
    //  previous session with the spectrum
  public:
    PreviousDbEntry( AuxWindow *dialog, WContainerWidget* container, SpecUtils::SpectrumType type,
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
            m_manager->displayFile( row, measurement, m_type, false, false, SpecMeasManager::VariantChecksToDo::None );
            m_dialog->hide();
            
            return;
          }//if( entry.id() == m_dbentry.id() )
        }//for( int row = 0; row < m_model->rowCount(); ++row )
        
        try
        {
          const int modelRow = m_manager->setDbEntry( m_dbentry, header,
                                                     measurement, true );
          m_manager->displayFile( modelRow, measurement, m_type, false, false, SpecMeasManager::VariantChecksToDo::None );
          m_dialog->hide();
        }catch( exception &e )
        {
          cerr << "\n\nPreviousDbEntry::dorevert()\n\tCaught: " << e.what() << "\n\n";
          passMessage( "Error displaying previous measurment, things may not"
                      " be as expected", WarningWidget::WarningMsgHigh );
        }//try / catch
      }else
      {
        m_header->setDbEntry( m_dbentry );
      }//if( m_dbentry->userHasModified )
      
      //      if( m_dialog )
      //        delete m_dialog;
    }//void dorevert()
    
    ~PreviousDbEntry() noexcept(true)
    {
      try
      {
        m_dbentry.reset();
      }catch(...)
      {
        cerr << "PreviousDbEntry destructo caught exception doing m_dbentry.reset()" << endl;
      }
      
      try
      {
        m_header.reset();
      }catch(...)
      {
        cerr << "PreviousDbEntry destructo caught exception doing m_header.reset()" << endl;
      }
      
    }//~PreviousDbEntry()
    
    AuxWindow *m_dialog;
    SpecUtils::SpectrumType m_type;
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
#endif // USE_DB_TO_STORE_SPECTRA

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
                      | AuxWindowProperties::PhoneNotFullScreen
                      | AuxWindowProperties::DisableCollapse
                      | AuxWindowProperties::SetCloseable) ),
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

#if( USE_QR_CODES )
void displayQrDialog( const vector<QRSpectrum::QrCodeEncodedSpec> urls, const size_t index,
                     const SpecUtils::SpectrumType type )
{
  WString seqnum;
  if( urls.size() <= 1 )
    seqnum = WString::tr("smm-qr-sequence-single");
  else
    seqnum = WString::tr("smm-qr-sequence-multi")
              .arg( static_cast<int>(index + 1) )
              .arg( static_cast<int>(urls.size()) );
  
  WString title;
  if( urls.size() > 1 )
    title = WString::tr("smm-qr-title").arg( seqnum );
  
  WString desc = WString::tr("smm-qr-desc")
                  .arg(seqnum)
                  .arg( WString::tr(SpecUtils::descriptionText(type)) );
  
  SimpleDialog *dialog = QrCode::displayTxtAsQrCode( urls[index].m_url, title, desc );
  
  if( (index + 1) < urls.size() )
  {
    WPushButton *btn = dialog->addButton( WString::tr("smm-next-qr") );
    btn->clicked().connect( std::bind([=](){
      displayQrDialog( urls, index + 1, type );
    }) );
  }//if( (index + 1) < urls.size() )
}//void displayQr( vector<QRSpectrum::QrCodeEncodedSpec> urls )

/** This dialog box gets shown when a multi-part QR code is received, but not all parts.
 Once all parts are handed to this dialog, it will load the spectrum, and close.
 
 TODO:
 - When part-1 of a spectrum is received, we could display some of the information/spectrum
 */
class MultiUrlSpectrumDialog : public SimpleDialog
{
  SpecMeasManager *m_manager;
  InterSpec *m_interspec;
  vector<SpecUtils::EncodedSpectraInfo> m_urls;
  
public:
  MultiUrlSpectrumDialog( SpecMeasManager *manager, InterSpec *viewer )
  : SimpleDialog( WString::tr("musd-dialog-title"), "&nbsp;" ),
    m_manager( manager ),
    m_interspec( viewer )
  {
    assert( manager );
    assert( viewer );

#if( IOS || ANDROID )
    WPushButton *btn = addButton( WString::tr("Cancel") );
#else
    addButton( WString::tr("Cancel") );
#endif
    
    finished().connect( m_manager, &SpecMeasManager::multiSpectrumDialogDone );
    
#if( IOS || ANDROID )
    btn->clicked().connect( this, &MultiUrlSpectrumDialog::removePrevUrlsFile );
    
    try
    {
      rapidxml::file<char> input_file( prevUrlFile().c_str() ); //Will throw if doesnt exist
      
      rapidxml::xml_document<char> doc;
      doc.parse<rapidxml::parse_default>( input_file.data() );
      
      XML_FOREACH_CHILD( url_node, &doc, "Url" )
      {
        addUrl( SpecUtils::xml_value_str(url_node) );
      }
    }catch( std::exception & )
    {
      removePrevUrlsFile();
    }//try / catch get previous results
#endif //#if( IOS || ANDROID )
  }//MultiUrlSpectrumDialog
  
  
#if( IOS || ANDROID )
  static std::string prevUrlFile()
  {
    const string datadir = InterSpec::writableDataDirectory();
    return SpecUtils::append_path(datadir, "prev_spec_uri.xml");
  }
  
  void updatePrevUrlsFile()
  {
    removePrevUrlsFile();
    
    if( m_urls.empty() )
      return;
    
    try
    {
      rapidxml::xml_document<char> doc;
      
      for( const SpecUtils::EncodedSpectraInfo &url : m_urls )
      {
        const char *value = doc.allocate_string( url.m_orig_url.c_str(), url.m_orig_url.size() + 1 );
        rapidxml::xml_node<> *node = doc.allocate_node( rapidxml::node_element, "Url", value );
        doc.append_node( node );
      }//for( const QRSpectrum::EncodedSpectraInfo &url : m_urls )
      
      ofstream prev_url_file( prevUrlFile().c_str(), ios::out | ios::binary );
      if( !prev_url_file.is_open() )
        throw runtime_error( "Unable to open file for output ('" + prevUrlFile() + "')" );
      
      rapidxml::print( std::ostream_iterator<char>(prev_url_file), doc, 0 );
    }catch( std::exception &e )
    {
      cerr << "Error writing prev_spec_uri.xml: " << e.what() << endl;
    }//try / catch
  }//void updatePrevUrlsFile()
  
  
  void removePrevUrlsFile()
  {
    try
    {
      const string datafile = prevUrlFile();
      
      if( SpecUtils::is_file(datafile) )
        SpecUtils::remove_file( datafile );
    }catch( std::exception & )
    {
    }
  }//void removePrevUrlsFile()
#endif //#if( IOS || ANDROID )
  
  void addUrl( const string &url_unencoded_url )
  {
    using SpecUtils::EncodedSpectraInfo;
    
    try
    {
      const EncodedSpectraInfo info = SpecUtils::get_spectrum_url_info( url_unencoded_url );
      assert( info.m_number_urls > 1 );
      
      
      // Check if incoming CRC-16 matches previous (i.e., is same spectrum)
      if( !m_urls.empty() && (info.m_crc != m_urls[0].m_crc) )
      {
        passMessage( WString::tr("musd-err-diff-spec"), WarningWidget::WarningMsgMedium );
        m_urls.clear();
      }
      
      // Check to see if we have already received this URL
      for( const SpecUtils::EncodedSpectraInfo &i : m_urls )
      {
        if( i.m_url_num == info.m_url_num )
        {
          passMessage( WString::tr("musd-err-duplicate-url"), WarningWidget::WarningMsgMedium );
          return;
        }
      }//for( const QRSpectrum::EncodedSpectraInfo &i : m_urls )
      
      m_urls.push_back( info );
      
      std::sort( begin(m_urls), end(m_urls),
          []( const EncodedSpectraInfo &lhs, const EncodedSpectraInfo &rhs ) -> bool {
        return lhs.m_url_num < rhs.m_url_num;
      } );
      
      if( m_urls.size() == info.m_number_urls )
      {
        vector<string> urls;
        for( const auto &i : m_urls )
          urls.push_back( i.m_orig_url );
        
        vector<SpecUtils::UrlSpectrum> specs = SpecUtils::decode_spectrum_urls( urls );
        assert( specs.size() == 1 );
        
        
        auto specmeas = make_shared<SpecMeas>();
        {
          shared_ptr<SpecUtils::SpecFile> specfile = SpecUtils::to_spec_file( specs );
          reinterpret_cast<SpecUtils::SpecFile &>( *specmeas ) = *specfile;
        }
        
        m_manager->multiSpectrumDialogDone();
        
        string display_name = WString::tr("musd-from-qr-file-name").toUTF8();
        if( specs.size() && (specs[0].m_title.size() > 2) )
          display_name = specs[0].m_title;
        
#if( IOS || ANDROID )
        removePrevUrlsFile();
#endif
        
        m_interspec->userOpenFile( specmeas, display_name );
        accept();
      }//if( m_urls.size() == info.m_number_urls )
      
      
      // TODO: we could show a partial spectrum, or some of its information here.
      //if( m_urls[0].m_spectrum_number == 0 )
      //  vector<UrlSpectrum> spec = spectrum_decode_first_url( m_urls[0].m_data );
      //std::shared_ptr<SpecUtils::SpecFile> to_spec_file( const std::vector<UrlSpectrum> &meas );
      
      string urls_list;
      for( size_t i = 0; i < m_urls.size(); ++i )
      {
        if( i && (i+1)== m_urls.size() )
          urls_list += " and ";
        else if( i )
          urls_list += ", ";
        urls_list += std::to_string( 1 + m_urls[i].m_url_num );
      }
      
      WString msg = WString::tr("musd-received-urls-list")
                      .arg(urls_list)
                      .arg( static_cast<int>(info.m_number_urls) );
      assert( m_msgContents );
      m_msgContents->setText( msg );
      
#if( IOS || ANDROID )
      updatePrevUrlsFile();
#endif
    }catch( std::exception &e )
    {
      assert( m_title );
      assert( m_msgContents );
      m_msgContents->setText( WString::tr("musd-err-decoding-url").arg(e.what()) );
      
#if( IOS || ANDROID )
      removePrevUrlsFile();
#endif
    }
  }//void addUrl( const string &url )
};//MultiUrlSpectrumDialog
#endif //USE_QR_CODES

WT_DECLARE_WT_MEMBER
(LookForQrCode, Wt::JavaScriptFunction, "LookForQrCode",
async function( sender_id, u8Buffer )
{
  let zxing = await ZXing();
  
  let zxingBuffer = zxing._malloc(u8Buffer.length);
  zxing.HEAPU8.set(u8Buffer, zxingBuffer);
  
  let results = zxing.readBarcodesFromImage(zxingBuffer, u8Buffer.length, true, "QRCode", 0xff);
  zxing._free(zxingBuffer);
  
  const firstRes = (results.size() > 0) ? results.get(0) : null;
  
  if( (results.size() === 0)
     || (firstRes && firstRes.error && (firstRes.error.length !== 0) && (!firstRes.text || firstRes.text.length)) )
  {
    if( firstRes && firstRes.error && (firstRes.error.length !== 0) )
    {
      console.log( "QR-code reading error:", firstRes.error );
      Wt.emit( sender_id, 'QrDecodedFromImg', -1, firstRes.error );
    }else
    {
      // No QR-code found
      console.log( "No QR-code found." );
      Wt.emit( sender_id, 'QrDecodedFromImg', 0, "" );
    }
  }else
  {
    console.log( "Successfully got QR-code results." );
    let uri = "";
    for( let i = 0; i < results.size(); i += 1)
    {
      const { format, text, bytes, error } = results.get(i);
      uri += (uri.length ? "\n" : "") + btoa(text);
      //console.log( "text: ", text, "\nerror:", error );
    }
    
    Wt.emit( sender_id, 'QrDecodedFromImg', results.size(), uri );
  }//if( error ) / else
}
);
  
  
WT_DECLARE_WT_MEMBER
(SearchForQrFromImgData, Wt::JavaScriptFunction, "SearchForQrFromImgData",
function( sender_id, b64_img_data_str )
{
  let binaryString = atob(b64_img_data_str);
  let bytes = new Uint8Array(binaryString.length);
  for( let i = 0; i < binaryString.length; i++)
    bytes[i] = binaryString.charCodeAt(i);
    
  let u8Buffer = new Uint8Array(bytes);
  
  Wt.WT.LookForQrCode(sender_id, u8Buffer).then( function(){
    console.log( "Done calling qr decode" );
  }).catch(function(err){
    console.log( "qr decode error: ", err );
    Wt.emit( sender_id, 'QrDecodedFromImg', -999, "There was an error calling decoding routine." );
  });
}
);
  
  
WT_DECLARE_WT_MEMBER
(SearchForQrUsingCanvas, Wt::JavaScriptFunction, "SearchForQrUsingCanvas",
function( sender_id, img )
{
  //Turn image into PNG, then call Wt.WT.LookForQrCode(sender_id, u8Buffer).then...
  try
  {
    let canvas = document.createElement('canvas');
    canvas.width = Math.max(img.width, img.naturalWidth ? img.naturalWidth : 0);
    canvas.height = Math.max(img.height, img.naturalHeight ? img.naturalHeight : 0);
    
    let ctx = canvas.getContext('2d');
    ctx.drawImage(img, 0, 0);
    
    // Convert the image to PNG format
    function blobToQr(blobData){
      blobData.arrayBuffer().then( function(arrBuf){
        let u8Buffer = new Uint8Array( arrBuf );
        Wt.WT.LookForQrCode(sender_id, u8Buffer).then( function(){
          console.log( "Done calling qr decode for canvas blob" );
        }).catch(function(err){
          console.log( "qr decode error: ", err );
          Wt.emit( sender_id, 'QrDecodedFromImg', -999, "Error calling decoding routine." );
        });
      });
    };//blobToQr
    
    canvas.toBlob( blobToQr, 'image/png' );
  }catch( err )
  {
    console.log( "qr decode error (1): ", err );
    Wt.emit( sender_id, 'QrDecodedFromImg', -999, "There was an error when calling decoding routine." );
  }
}
);
  
  
class UploadedImgDisplay : public WContainerWidget
{
protected:
  InterSpec *m_viewer;
  std::string m_display_name;
  std::string m_mimetype;
  
  SpecUtils::SpectrumType m_upload_type;
  
  /** I like the QR -code auto search running automatically when a image is detected, but
   I'm not quite yet confident enough in the implementation to have things just run.
   
   After some more use/testing (including on the different platforms), we'll maybe always auto search.
   */
  bool m_autoQrCodeSearch;
  
  WImage *m_image;
  WMemoryResource *m_resource;
  
  /** Will be nullptr. if `m_autoQrCodeSearch == true` */
  WPushButton *m_checkForQrCodeBtn;
  
  /** Will be nullptr. if `m_autoQrCodeSearch == false` */
  WText *m_qrCodeStatusTxt;
  
  JSignal<int,std::string> m_qrDecodeSignal;
  boost::function<void()> m_close_parent_dialog;
  
  void close_parent_dialog()
  {
    m_close_parent_dialog();
  }
  
  void apply_uris( const vector<string> &uris )
  {
    assert( !uris.empty() );
    if( uris.empty() )
      return;
    
    // For spectrum files, we could just call `m_viewer->handleAppUrl( uri );`, but then we
    //  cant preserve if the user dropped it in as a foreground, background, or secondary,
    //  and also I'm not sure how multi-qr-code URIs would work out, so we'll decode and
    //  load spectrum file URIs totally here.
    bool is_spec_uris = true;
    for( const string &uri : uris )
    {
      is_spec_uris &= (SpecUtils::istarts_with(uri, "RADDATA://")
                      || SpecUtils::istarts_with(uri, "interspec://G0/"));
    }
    
    if( is_spec_uris )
    {
      try
      {
        bool decoded = false;
        SpecUtils::EncodedSpectraInfo info;
        
        // The URIs may be percent encoded - we'll try up to 3 levels of decoding them
        //  (probably not the case for QR-codes that they would ever be encoded more than
        //  once, but JIC)
        vector<string> unencoded_uris = uris;
        for( size_t i = 0; !decoded && (i < 3); ++i )
        {
          try
          {
            info = SpecUtils::get_spectrum_url_info( unencoded_uris.front() );
            decoded = true;
          }catch( std::exception & )
          {
            // Assume all URIs are percent encoded the same number of times
            for( string &uri : unencoded_uris )
              uri = SpecUtils::url_decode( uri );
          }//
        }//for( try to decode URL as spectrum )
        
        if( decoded )
        {
          shared_ptr<SpecMeas> specmeas;
          
          vector<SpecUtils::UrlSpectrum> spectra;
          if( info.m_number_urls == 1 )
          {
            spectra = SpecUtils::spectrum_decode_first_url( unencoded_uris.front() );
          }else if( (info.m_number_urls > 1) && (info.m_number_urls == unencoded_uris.size()) )
          {
            spectra = SpecUtils::decode_spectrum_urls( unencoded_uris );
          }//if( one URL ) / else
          
          if( spectra.empty() )
            throw runtime_error( "No gamma measurements in URL/QR code" );
          
          shared_ptr<SpecUtils::SpecFile> specfile = SpecUtils::to_spec_file( spectra );
          assert( specfile && specfile->num_measurements() );
          
          if( !specfile || (specfile->num_measurements() < 1) || (specfile->num_gamma_channels() < 7) )
            throw runtime_error( "No gamma measurements in URL/QR code" );
          
          string display_name = m_display_name.empty() ? WString::tr("uid-qr-code-pic").toUTF8() : m_display_name;
          
          specmeas = make_shared<SpecMeas>();
          reinterpret_cast<SpecUtils::SpecFile &>( *specmeas ) = *specfile;
          specmeas->set_filename( display_name );
          
          if( specmeas && (specmeas->num_measurements() >= 1) )
          {
            SpecMeasManager *fileManager = m_viewer->fileManager();
            SpectraFileModel *fileModel = fileManager->model();
            
            auto header = make_shared<SpectraFileHeader>( true, m_viewer );
            header->setFile( display_name, specmeas );
            fileManager->addToTempSpectrumInfoCache( specmeas );
            const int row = fileModel->addRow( header );
            fileManager->displayFile( row, specmeas, m_upload_type, true, true, SpecMeasManager::VariantChecksToDo::None );
            
            m_close_parent_dialog();
            return;
          }//if( specmeas && (specmeas->num_measurements() >= 1) )
        }//if( decoded )
      }catch( std::exception &e )
      {
        cerr << "Failed to decode spectrum URI from image - will let handleAppUrl() deal with things: " << e.what() << endl;
      }
    }//if( a spectrum URI )
    
    // If we are here - its not valid spectrum file URIs, so we'll just handle them as if the OS
    //  had past them off to the app.
    for( const string uri : uris )
      m_viewer->handleAppUrl( uri );
    
    m_close_parent_dialog();
  }//void apply_uris( const vector<string> &uris )
  
  
  void qr_check_result( const int num_qr, const string b64_value )
  {
    vector<string> initial_uris;
    if( num_qr > 1 )
    {
      SpecUtils::split( initial_uris, b64_value, "\n\r" );
    }else
    {
      initial_uris.push_back( b64_value );
    }
    
    if( (num_qr <= 0) || initial_uris.empty() || (initial_uris.size() > 14) )
    {
      WString title, content;
      
      if( (num_qr == 0) && b64_value.empty() )
      {
        if( m_autoQrCodeSearch && m_qrCodeStatusTxt )
        {
          m_qrCodeStatusTxt->setText( WString::tr("uid-no-qr-found") );
          return;
        }
        
        title = WString::tr("uid-no-qr-found-title");
        content = WString::tr("uid-no-qr-content");
      }else
      {
        title = WString::tr("uid-err-qr-search-title");
        content = WString::tr("uid-err-qr-search-content");
        content.arg( num_qr );
        
        if( b64_value.size() < 128 )
          content.arg( Wt::Utils::htmlEncode(b64_value) );
        else
          content.arg( Wt::Utils::htmlEncode(b64_value.substr(0,125) + "...") );
        
        if( m_autoQrCodeSearch && m_qrCodeStatusTxt )
        {
          m_qrCodeStatusTxt->setText( WString("{1} &#9432;").arg(title) );
          passMessage( WString::tr("<div>{1}</div>{2}").arg(title).arg(content), WarningWidget::WarningMsgHigh );
          return;
        }
      }//if( (num_qr == 0) && b64_value.empty() )
      
      SimpleDialog *dialog = new SimpleDialog( title, content );
      dialog->addButton( WString::tr("Okay") );
      
      return;
    }//if( num_qr <= 0 )
    

    vector<string> cleaned_up_uris;
    for( string &uri : initial_uris )
    {
      uri = Wt::Utils::base64Decode(uri);
      
      size_t uri_pos = SpecUtils::ifind_substr_ascii( uri, "interspec://" );
      if( uri_pos == string::npos )
        uri_pos = SpecUtils::ifind_substr_ascii( uri, "raddata://" );
      
      //TODO: maybe also accept text with a "G0/xxx/" or "G0/xxxx/" anywhere (where 'x' is a hex digit).
      
      if( uri_pos == 0 )
        cleaned_up_uris.push_back( uri );
      else if( uri_pos != string::npos )
        cleaned_up_uris.push_back( uri.substr(uri_pos) );
    }//for( string &val : initial_uris )
    
    
    if( (cleaned_up_uris.size() != initial_uris.size()) || cleaned_up_uris.empty() )
    {
      // Right now we'll be conservative, and if any QR code had a invalid URI, we wont
      //  use any of them
      WString content;
      
      if( initial_uris.size() > 1 )
      {
        const size_t num_invalid = initial_uris.size() -  cleaned_up_uris.size();
        content = WString::tr("uid-some-invalid-qrs").arg( static_cast<int>(num_invalid) );
      }else
      {
        content = WString::tr("uid-qr-didnt-have-uri");
      }
 
      if( m_autoQrCodeSearch && m_qrCodeStatusTxt )
      {
        m_qrCodeStatusTxt->setText( WString::tr("uid-invalid-qr") );
        return;
      }
      
      SimpleDialog *dialog = new SimpleDialog( WString::tr("uid-invalid-uri"), content );
      dialog->addButton( WString::tr("Okay") );
      
      return;
    }//if( cleaned_up_uris.size() != initial_uris.size() )
    
    WString title, content;
    if( cleaned_up_uris.size() > 1 )
    {
      title = WString::tr("uid-multi-qr-found-title").arg( static_cast<int>(cleaned_up_uris.size()) );
      content = WString::tr("uid-multi-qr-found-content");
      
      if( m_autoQrCodeSearch && m_qrCodeStatusTxt )
        m_qrCodeStatusTxt->setText(  WString::tr("uid-multi-qr-status").arg( static_cast<int>(cleaned_up_uris.size()) )  );
    }else if( cleaned_up_uris.size() > 0 )
    {
      title = WString::tr("uid-qr-found-title");
      
      assert( !cleaned_up_uris.empty() );
      string short_uri = cleaned_up_uris.front();
      SpecUtils::utf8_limit_str_size(short_uri, 32);
      
      content = WString::tr("uid-qr-found-content").arg( short_uri );
      
      if( m_autoQrCodeSearch && m_qrCodeStatusTxt )
        m_qrCodeStatusTxt->setText( WString::tr("uid-found-qr") );
    }//if( cleaned_up_uris.size() > 1 )
    
    SimpleDialog *dialog = new SimpleDialog( title, content );
    WPushButton *btn = dialog->addButton( WString::tr("Yes") );
    btn->clicked().connect( boost::bind(&UploadedImgDisplay::apply_uris, this, cleaned_up_uris)  );
    btn->clicked().connect( this, &UploadedImgDisplay::close_parent_dialog );
    
    btn = dialog->addButton( WString::tr("No") );
  }//void qr_check_result( const int num_qr, const string b64_value )
  
  
  void check_for_qr_with_raw()
  {
    LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsLookForQrCode);
    LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsSearchForQrFromImgData);
    wApp->require( "InterSpec_resources/assets/js/zxing-cpp-wasm/zxing_reader.js", "zxing_reader.js" );
    
    vector<unsigned char> raw_data = m_resource->data();
    string str_data( raw_data.size(), '\0' );
    memcpy( (void *)&(str_data[0]), (void *)raw_data.data(), raw_data.size() );
    string b64_data = Wt::Utils::base64Encode(str_data, false);
    this->doJavaScript( "Wt.WT.SearchForQrFromImgData('" + this->id() + "','" + b64_data + "');" );
  }//void check_for_qr_with_raw()
  
  
  void check_for_qr_from_canvas()
  {
    LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsLookForQrCode);
    LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsSearchForQrUsingCanvas);
    wApp->require( "InterSpec_resources/assets/js/zxing-cpp-wasm/zxing_reader.js", "zxing_reader.js" );
    
    this->doJavaScript( "Wt.WT.SearchForQrUsingCanvas('" + this->id() + "'," + m_image->jsRef() + ");" );
  }//void check_for_qr_from_canvas()
  
  
  void embed_in_n42()
  {
    if( !m_resource )
      return;
    
    const vector<unsigned char> data = m_resource->data();
    if( data.empty() )
    {
      passMessage( WString::tr("uid-err-getting-data"), WarningWidget::WarningMsgHigh );
      return;
    }
          
    shared_ptr<SpecMeas> meas = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
    if( !meas )
    {
      passMessage( WString::tr("uid-err-no-foreground-embedd"), WarningWidget::WarningMsgHigh );
      return;
    }
          
    SpecUtils::MultimediaData multi;
    multi.remark_ = "Image file embedded using InterSpec.";
    multi.descriptions_= "filename: " + m_display_name;
          
    string data_str;
    data_str.resize( data.size() );
    memcpy( &(data_str[0]), data.data(), data.size() );
    const string base_64_encoded = Wt::Utils::base64Encode( data_str );
    multi.data_.resize( base_64_encoded.size() );
    memcpy( multi.data_.data(), base_64_encoded.data(), base_64_encoded.size() );
    multi.data_encoding_ = SpecUtils::MultimediaData::EncodingType::BinaryBase64;
    multi.capture_start_time_ = SpecUtils::time_point_t{};
    multi.file_uri_ = m_display_name;
    multi.mime_type_ = m_mimetype;
          
    meas->add_multimedia_data( multi );
    m_viewer->checkEnableViewImageMenuItem();
    
    passMessage( WString::tr("uid-image-embed-conf"), WarningWidget::WarningMsgInfo );
          
    //wApp->doJavaScript( "$('#" + dialog->id() + "').hide(); $('.Wt-dialogcover').hide();" );
    wApp->doJavaScript( "$('.Wt-dialogcover').hide();" );
      
    m_close_parent_dialog();
  }//void embed_in_n42()
  
public:
  
  UploadedImgDisplay( InterSpec *interspec,
                     const std::string &display_name, 
                     const char *mimetype,
                     const size_t filesize,
                     std::ifstream &infile,
                     SimpleDialog *dialog,
                     SpecUtils::SpectrumType type )
    : WContainerWidget( dialog->contents() ),
  m_viewer( interspec ),
  m_display_name( display_name ),
  m_mimetype( mimetype ),
  m_upload_type( type ),
  m_autoQrCodeSearch( true ),
  m_image( nullptr ),
  m_resource( nullptr ),
  m_checkForQrCodeBtn( nullptr ),
  m_qrCodeStatusTxt( nullptr ),
  m_qrDecodeSignal( this, "QrDecodedFromImg", false)
  {
    m_close_parent_dialog = wApp->bind( boost::bind( &SimpleDialog::done, dialog, Wt::WDialog::DialogCode::Accepted ) );
    
    WText *t = new WText( WString::tr("uid-image-not-spectrum"), this );
    t->addStyleClass( "NonSpecOtherFile" );
    
    const size_t max_disp_size = 16*1024*1024;
    
    if( filesize > max_disp_size )
    {
      WText *errort = new WText( WString::tr("uid-err-to-large"), this );
      errort->addStyleClass( "NonSpecError" );
      return;
    }//if( filesize > max_disp_size )
      
    vector<uint8_t> totaldata( filesize );
    const bool success = infile.read( (char *)&(totaldata[0]), filesize ).good();
      
    if( !success )
    {
      WText *errort = new WText( WString::tr("uid-err-reading-upload"), this );
      errort->addStyleClass( "NonSpecError" );
      return;
    }//if( !success )
    
    const bool is_gif_jpg_png = (SpecUtils::icontains(mimetype, "gif")
                              || SpecUtils::icontains(mimetype, "jpeg")
                              || SpecUtils::icontains(mimetype, "png"));
    
    m_resource = new WMemoryResource( mimetype, this );
    m_resource->setData( totaldata );
        
    m_image = new WImage();
    m_image->setImageLink( WLink(m_resource) );
    m_image->addStyleClass( "NonSpecImgFile" );
    addWidget( m_image );
        
    WContainerWidget *btn_div = new WContainerWidget( this );
    btn_div->addStyleClass( "ImgFileBtnBar" );
        
    if( m_autoQrCodeSearch )
    {
      m_qrCodeStatusTxt = new WText( btn_div );
      HelpSystem::attachToolTipOn( m_qrCodeStatusTxt, WString::tr("uid-tt-auto-qr"),
                                  true, HelpSystem::ToolTipPosition::Right );
    }else
    {
      m_checkForQrCodeBtn = new WPushButton( WString::tr("uid-check-for-qr-btn"), btn_div );
      m_checkForQrCodeBtn->setStyleClass( "LinkBtn NonSpecQrCodeBtn" );
    }//if( m_autoQrCodeSearch ) / else
    
    
    m_qrDecodeSignal.connect( boost::bind( &UploadedImgDisplay::qr_check_result, this,
                                          boost::placeholders::_1, boost::placeholders::_2 ) );
    
    if( is_gif_jpg_png )
    {
      // We will use the raw image data
      if( m_autoQrCodeSearch )
      {
        assert( m_qrCodeStatusTxt );
        m_qrCodeStatusTxt->setText( WString::tr("uid-looking-for-qr") );
        check_for_qr_with_raw();
      }else
      {
        m_checkForQrCodeBtn->clicked().connect( this, &UploadedImgDisplay::check_for_qr_with_raw );
      }
    }else
    {
      // We will draw the image to a canvas, and use that
      if( m_autoQrCodeSearch )
      {
        LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsLookForQrCode);
        LOAD_JAVASCRIPT(wApp, "SpecMeasManager.cpp", "SpecMeasManager", wtjsSearchForQrUsingCanvas);
        wApp->require( "InterSpec_resources/assets/js/zxing-cpp-wasm/zxing_reader.js", "zxing_reader.js" );
        
        assert( m_qrCodeStatusTxt );
        m_qrCodeStatusTxt->setText( WString::tr("uid-no-qr-found") );
        
        // If `imageLoaded()` is never called, it means the image couldnt be displayed, for example
        //  if image file is invalid, or a HEIC on Windows.
        m_image->imageLoaded().connect( "function(){ Wt.WT.SearchForQrUsingCanvas('" + this->id() + "'," + m_image->jsRef() + "); }" );
        m_image->imageLoaded().connect( boost::bind( &WText::setText, m_qrCodeStatusTxt, WString("Looking for QR-codes.") ) );
      }else
      {
        m_checkForQrCodeBtn->clicked().connect( this, &UploadedImgDisplay::check_for_qr_from_canvas );
        m_checkForQrCodeBtn->disable();
        m_image->imageLoaded().connect( m_checkForQrCodeBtn, &WPushButton::enable );
      }
    }//if( zxing can read the image directly ) / else
    
        
    if( m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) )
    {
      WPushButton *embedbtn = new WPushButton( WString::tr("uid-embed-image-btn"), btn_div );
      embedbtn->setStyleClass( "LinkBtn NonSpecEmbedBtn" );
      embedbtn->clicked().connect( this, &UploadedImgDisplay::embed_in_n42 );
    }//if( m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) )
  }//UploadedImgDisplay constructor
  
};//UploadedImgDisplay
  
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
    
    WText *uploadText = new WText( WString("{1}: ").arg( WString::tr("Foreground") ) );
    WFileUpload *m_fileUpload = new WFileUpload(  );
    m_fileUpload->changed().connect( m_fileUpload, &Wt::WFileUpload::upload );
    m_fileUpload->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager, m_fileUpload, SpectrumType::Foreground));
    m_fileUpload->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge,
                                                      boost::placeholders::_1 ) );

    WText *uploadText2 = new WText( WString("{1}: ").arg( WString::tr("Background") ) );
    WFileUpload *m_fileUpload2 = new WFileUpload(  );
    m_fileUpload2->changed().connect( m_fileUpload2, &Wt::WFileUpload::upload );
    m_fileUpload2->uploaded().connect( boost::bind( &SpecMeasManager::dataUploaded2, m_manager, m_fileUpload2, SpectrumType::Background));
    m_fileUpload2->fileTooLarge().connect( boost::bind( &SpecMeasManager::fileTooLarge,
                                                       boost::placeholders::_1 ) );
    
    WText *uploadText3 = new WText( WString("{1}: ").arg( WString::tr("second-foreground") ) );
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
    m_spectrumManagerWindow( nullptr ),
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
    m_removeFore2Button ( NULL),
    m_foregroundDragNDrop( new FileDragUploadResource(this) ),
    m_secondForegroundDragNDrop( new FileDragUploadResource(this) ),
    m_backgroundDragNDrop( new FileDragUploadResource(this) ),
    m_multiUrlSpectrumDialog( nullptr ),
    m_destructMutex( new std::mutex() ),
    m_destructed( new bool(false) ),
    m_previousStatesDialog( nullptr ),
    m_processingUploadDialog( nullptr ),
    m_nonSpecFileDialog( nullptr ),
    m_processingUploadTimer{}
{
  std::unique_ptr<UndoRedoManager::BlockUndoRedoInserts> undo_blocker;
  if( viewer && viewer->undoRedoManager() )
    undo_blocker.reset( new UndoRedoManager::BlockUndoRedoInserts() );
  
  wApp->useStyleSheet( "InterSpec_resources/SpecMeasManager.css" );
  if( viewer )
    viewer->useMessageResourceBundle( "SpecMeasManager" );
  
  m_treeView = new RowStretchTreeView();
  m_fileModel = new SpectraFileModel( m_treeView );
  m_treeView->setModel( m_fileModel );
  m_treeView->selectionChanged().connect( boost::bind( &SpecMeasManager::selectionChanged, this ) );

  // TODO: (20241028) it doesn't appear necessary to show and then delete the spectrum manager window - but leaving until after v1.0.13 release
  startSpectrumManager(); //initializes
  deleteSpectrumManager(); //deletes instance
    
  m_sql = viewer->sql();

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
  
  m_foregroundDragNDrop->setUploadProgress( true );
  m_foregroundDragNDrop->dataReceived().connect( boost::bind( &SpecMeasManager::handleDataRecievedStatus, this,
                                      boost::placeholders::_1, boost::placeholders::_2, SpectrumType::Foreground ) );
  
  m_secondForegroundDragNDrop->setUploadProgress( true );
  m_secondForegroundDragNDrop->dataReceived().connect( boost::bind( &SpecMeasManager::handleDataRecievedStatus, this,
                                      boost::placeholders::_1, boost::placeholders::_2, SpectrumType::SecondForeground ) );
  
  m_backgroundDragNDrop->setUploadProgress( true );
  m_backgroundDragNDrop->dataReceived().connect( boost::bind( &SpecMeasManager::handleDataRecievedStatus, this,
                                      boost::placeholders::_1, boost::placeholders::_2, SpectrumType::Background ) );
}// SpecMeasManager

//Moved what use to be SpecMeasManager, out to a startSpectrumManager() to correct modal issues
void SpecMeasManager::startSpectrumManager()
{
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );

  if( m_spectrumManagerWindow )
  {
    m_spectrumManagerWindow->show();
    m_spectrumManagerWindow->centerWindow();
    m_spectrumManagerWindow->resizeToFitOnScreen();
    
    return;
  }//if( m_spectrumManagerWindow )
  
    m_spectrumManagerWindow = new AuxWindow( WString::tr("smm-window-title"),
                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsModal)
                     | AuxWindowProperties::TabletNotFullScreen
                     | AuxWindowProperties::DisableCollapse
                     | AuxWindowProperties::SetCloseable
                     ) );
    m_spectrumManagerWindow->addStyleClass( "SpecMeasManager" );
    
    WContainerWidget *title = m_spectrumManagerWindow->titleBar();
    title->addStyleClass( "SpectraFileManagerHeader" );
    
    m_spectrumManagerWindow->finished().connect( boost::bind( &SpecMeasManager::displayIsBeingHidden, this ) );
    // collapsed().connect( boost::bind( &SpecMeasManager::displayIsBeingHidden, this ) );
    m_spectrumManagerWindow->expanded().connect( boost::bind( &SpecMeasManager::displayIsBeingShown, this ) );
    
    WContainerWidget *uploadDiv = new WContainerWidget( );
    uploadDiv->setStyleClass( "uploadDiv" );
    
    //WText* spec =
    new WText( WString::tr("smm-load-spec-from"), uploadDiv);
    //spec->setIcon( "InterSpec_resources/images/plus_min_white.svg" );
    
    Wt::WPushButton* uploadButton = new Wt::WPushButton( WString::tr("app-mi-file-open"),uploadDiv);
    uploadButton->clicked().connect(  this, &SpecMeasManager::uploadSpectrum );
    HelpSystem::attachToolTipOn(uploadButton, WString::tr("smm-tt-upload-file"),
                                showToolTips, HelpSystem::ToolTipPosition::Bottom );
    uploadButton->setIcon( "InterSpec_resources/images/file_search.png" );
    uploadButton->setMargin(10,Wt::Left);
    
#if( USE_DB_TO_STORE_SPECTRA )
    Wt::WPushButton* importButton = new Wt::WPushButton( WString::tr("app-mi-file-prev"), uploadDiv );
    importButton->clicked().connect( boost::bind( &SpecMeasManager::browsePrevSpectraAndStatesDb, this ) );
    HelpSystem::attachToolTipOn(importButton, WString::tr("app-mi-tt-file-prev"), 
                                showToolTips, HelpSystem::ToolTipPosition::Bottom );
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
  
  if( m_viewer && m_viewer->undoRedoManager() && m_viewer->undoRedoManager()->canAddUndoRedoNow() )
    new UndoRedoManager::BlockGuiUndoRedo( m_spectrumManagerWindow ); // BlockGuiUndoRedo is WObject, so this `new` doesnt leak
  
  UndoRedoManager *undoRedo = m_viewer->undoRedoManager();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      InterSpec *viewer = InterSpec::instance();
      SpecMeasManager *manager = viewer->fileManager();
      if( manager )
        manager->deleteSpectrumManager();
    };
    auto redo = [](){
      InterSpec *viewer = InterSpec::instance();
      SpecMeasManager *manager = viewer->fileManager();
      if( manager )
        manager->startSpectrumManager();
    };
    
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Open Spectrum Manager.");
  }//if( undo )
} //startSpectrumManager()


SpecMeasManager::~SpecMeasManager()
{
  std::lock_guard<std::mutex> lock( *m_destructMutex );
  
  (*m_destructed) = true;
  
  if( m_nonSpecFileDialog )
    delete m_nonSpecFileDialog;
} // SpecMeasManager::~SpecMeasManager()


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

  throw std::runtime_error( "Serious problem in SpecMeasManager::dragNDrop(..)" );
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
    passMessage( WString::tr("smm-err-zip"), 2 );
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
    
    if( !app )
      app = dynamic_cast<WApplication *>( m_viewer->parent() );
    
    assert(app);

    WApplication::UpdateLock lock(app);

    if (!lock)
    {
      cerr << "SpecMeasManager::handleZippedFile: failed to get app lock." << endl;
      return false;
    }

    // Make sure we have the CSS we need
    app->useStyleSheet("InterSpec_resources/SpecMeasManager.css");


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

    WString txt = WString::tr( name.size() ? "smm-zip-sel-file" : "smm-zip-sel-file-no-name")
                  .arg( name )
                  .arg( (m_viewer->isPhone() ? "" : "<br />") );
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
    model->setHeaderData(  0, Horizontal, WString::tr("smm-zip-filename-hdr"), DisplayRole );
    model->setHeaderData(  1, Horizontal, WString::tr("smm-zip-kilobyte-hdr"), DisplayRole );
    WItemDelegate *delegate = new WItemDelegate( table );
    delegate->setTextFormat( "%.1f" );
    table->setItemDelegateForColumn( 1, delegate );
    
    
    WContainerWidget *typecb = new WContainerWidget();
    WGridLayout *cblayout = new WGridLayout( typecb );
    cblayout->setContentsMargins( 2, 0, 0, 0 );
    cblayout->setVerticalSpacing( 0 );
    cblayout->setHorizontalSpacing( 0 );
    WButtonGroup *group = new WButtonGroup(typecb);
    Wt::WRadioButton *button = new WRadioButton( WString::tr("Foreground") );
    cblayout->addWidget( button, 0, 0, AlignCenter );
    group->addButton(button, toint(SpectrumType::Foreground) );
    button = new WRadioButton( WString::tr("Background"), typecb);
    cblayout->addWidget( button, 0, 1, AlignCenter );
    group->addButton(button, toint(SpectrumType::Background) );
    button = new Wt::WRadioButton( WString::tr("Secondary"), typecb);
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
    
    AuxWindow *window = new AuxWindow( WString::tr("smm-zip-window-title"),
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
      model->setItem( static_cast<int>(i), 0, item );
      item = new WStandardItem();
      
      const double sizekb = uncompresssize[i]/1024.0;
      
      //item->setText( buff );
      item->setData( sizekb, DisplayRole );
      model->setItem( static_cast<int>(i), 1, item );
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
    Wt::WPushButton *closeButton = window->addCloseButtonToFooter( WString::tr("Cancel"),true);
    closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
    WPushButton *openButton = new WPushButton( WString::tr("smm-zip-display-btn") );
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
                                             const std::string &fileLocation,
                                             SpecUtils::SpectrumType type )
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
    passMessage( WString::tr("smm-empty-file"), 2 );
    return true;
  }
  
  uint8_t data[1024] = { 0x0 };
  
  if( !infile.read( (char *)data, std::min(sizeof(data), filesize) ) )
  {
    passMessage( WString::tr("smm-failed-reading-nonspec-file"), 2 );
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
  
  // SpecMeasManager will keep track of this next dialog, so we can do undo/redo a little better
  SimpleDialog *dialog = new SimpleDialog();
  
  if( m_nonSpecFileDialog )
  {
    delete m_nonSpecFileDialog; //The `destroyed()` signal will set `m_nonSpecFileDialog` to nullptr
    cerr << "m_nonSpecFileDialog was not nullptr" << endl;
  }
  assert( !m_nonSpecFileDialog );
  
  m_nonSpecFileDialog = dialog;
  cout << "Assigning m_nonSpecFileDialog to " << dialog << endl;
  
  // The dialog may get deleted, and never accepted, so we will hook up to the
  //  `destroyed()` signal to keep track of `m_nonSpecFileDialog`
  m_nonSpecFileDialog->destroyed().connect( std::bind([dialog](){
    InterSpec *interspec = InterSpec::instance();
    SpecMeasManager *manager = interspec ? interspec->fileManager() : nullptr;
    assert( manager );
    if( !manager )
      return;
    
    cerr << "deleting m_nonSpecFileDialog of" << manager->m_nonSpecFileDialog << endl;
    
    assert( !manager->m_nonSpecFileDialog || (dialog == manager->m_nonSpecFileDialog) );
    
    manager->m_nonSpecFileDialog = nullptr;
  }) );
  
  dialog->addButton( WString::tr("Close") );
  WContainerWidget *contents = dialog->contents();
  contents->addStyleClass( "NonSpecDialogBody" );
  WText *title = new WText( WString::tr("smm-not-spec-file"), contents );
  title->addStyleClass( "title" );
  
  auto add_undo_redo = [dialog, displayName, filesize, &infile, type](){
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( !undoRedo )
      return;
    
    auto closeDialog = [](){
      InterSpec *interspec = InterSpec::instance();
      SpecMeasManager *manager = interspec ? interspec->fileManager() : nullptr;
      assert( manager );
      if( manager )
        manager->closeNonSpecFileDialog();
    };
    
    std::function<void()> reOpenDialog;
    
    // If input size isnt too horrible large, we'll save it in memory and create a redo point
    if( filesize <= 10*1024*1024 )
    {
      infile.clear();
      infile.seekg(0, std::ios::beg);
        
      auto file_data = make_shared<vector<uint8_t>>( filesize, '\0' );
      infile.read( (char *)(&((*file_data)[0])), static_cast<streamsize>(filesize) );
      infile.clear();
      infile.seekg(0, std::ios::beg);
      
      // If file is more than 64 kb (arbitrarily chosen), we will only keep the data in memory for
      //  ~5 minutes
      std::weak_ptr<vector<uint8_t>> wk_ptr = file_data;
      if( filesize > 64*1024 )
      {
        shared_ptr<vector<uint8_t>> data_ptr = file_data;
        file_data = nullptr;
        
        auto clear_mem = [data_ptr](){
          cout << "Clearing memory of file: " << data_ptr.use_count() << endl;
        };
        
        const int milliSeconds = 5*60*1000;
        WServer::instance()->schedule( milliSeconds, wApp->sessionId(), clear_mem, clear_mem );
      }//if( filesize > 64*1024 )
      
      const std::string dispname = displayName;
      
      reOpenDialog = [wk_ptr, file_data, dispname, type](){
        
        shared_ptr<vector<uint8_t>> data_ptr = wk_ptr.lock();
        if( !data_ptr )
        {
          cout << "Was NOT able to get `file_data` - must have been cleared in memory" << endl;
          return;
        }else
        {
          // Not sure if we need to force the compiler to potentially keep `file_data` around
          cout << "file_data.use_count=" << file_data.use_count() << endl;
        }
        
        InterSpec *interspec = InterSpec::instance();
        SpecMeasManager *manager = interspec ? interspec->fileManager() : nullptr;
        assert( interspec && manager );
        if( !manager )
          return;
        
        const string tmpdir = SpecUtils::temp_dir();
        const string tmpname = SpecUtils::temp_file_name( "non_spec_file", tmpdir );
        assert( !SpecUtils::is_file(tmpname) );
        
        if( SpecUtils::is_file(tmpname) )
        {
          cerr << "Unexpectedly tmp filename is a file: '" << tmpname << "'." << endl;
          return;
        }
        
        bool wrote_tmp_file = false;
          
        {//begin write tmp file
#ifdef _WIN32
          const std::wstring wtmpfile = SpecUtils::convert_from_utf8_to_utf16(tmpname);
          ofstream outfilestrm( wtmpfile.c_str(), ios::out | ios::binary );
#else
          ofstream outfilestrm( tmpname.c_str(), ios::out | ios::binary );
#endif
          wrote_tmp_file = (outfilestrm
                              && outfilestrm.write( (char *)data_ptr->data(), data_ptr->size()) );
        }//end write tmp file
          
        if( wrote_tmp_file )
        {
          try
          {
            manager->handleNonSpectrumFile( dispname, tmpname, type );
          }catch( std::exception &e )
          {
            cerr << "Unexpected exception from `handleNonSpectrumFile(...)`" << endl;
          }
            
          SpecUtils::remove_file( tmpname );
        }else
        {
          cerr << "Failed to write '" << dispname << "' to temporary file." << endl;
        }
      };//redo lambda
    }//if( filesize <= 256*1024 )
    
    if( undoRedo->canAddUndoRedoNow() )
      undoRedo->addUndoRedoStep( closeDialog, reOpenDialog, "Open non-spectrum file dialog." );
    
    // Find close/cancel button
    vector<WWidget *> btns;
    if( dialog->footer() )
      btns = dialog->footer()->children();
    for( WWidget *w : btns )
    {
      WPushButton *btn = dynamic_cast<WPushButton *>( w );
      if( !btn || ((btn->text().key() != "Close") && (btn->text().key() != "Cancel")) )
      {
        // TODO: We could add this undo/redo step for all buttons, so in the case the
        //       user accepted the DRF, or peaks to fit, or whatever, that action should already
        //       be hooked up to undo/redo, but right now if we hooked it up, they would have to do
        //       undo again...  maybe when we aggregate multiple undo/redo steps into a single one
        //       every WApplication event loop, then we can enable this.
        continue;
      }
      
      btn->clicked().connect( std::bind([=](){
        // The `undo` call wont do anything, since it is currently bound to a now deleted dialog,
        //  but leaving in in case we upgrade things a bit in the future to have the dialog tracked
        //  by SpecMeasManager
        undoRedo->addUndoRedoStep( reOpenDialog, closeDialog, "Close non-spectrum file dialog." );
      }) );
    }//for( WWidget *w : btns )
  };//add_undo_redo
  
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
      WText *t = new WText( WString::tr("smm-icd2"), contents );
      t->addStyleClass( "NonSpecOtherFile" );
      add_undo_redo();
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
  const uint8_t heic_mn[]     = { 0x66, 0x74, 0x79, 0x70, 0x68, 0x65, 0x69, 0x63 }; //HEIC (apple) pictures have "ftypheic" at offset 4
  
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
  const bool isheic = (0 == memcmp( heic_mn, data + 4, sizeof(heic_mn) ));
  
  
  // Wt::Utils::guessImageMimeTypeData(<#const std::vector<unsigned char> &header#>)
  
  if( iszip )
  {
    //zip (but can be xlsx, pptx, docx, odp, jar, apk) 50 4B 03 04
    WText *t = new WText( WString::tr("smm-invalid-zip"), contents );
    t->addStyleClass( "NonSpecOtherFile" );
    
    add_undo_redo();
    return true;
  }//if( iszip )
  
  if( israr || istar || iszip7 || isgz )
  {
    WText *t = new WText( WString::tr("smm-unsupported-archive"), contents );
    t->addStyleClass( "NonSpecOtherFile" );
    
    add_undo_redo();
    
    return true;
  }//if( israr || istar || iszip7 || isgz )
  
  
  if( ispdf | isps | istif )
  {
    WText *t = new WText( WString::tr("smm-unsupported-document"), contents );
    t->addStyleClass( "NonSpecOtherFile" );
    
    add_undo_redo();
    
    return true;
  }//if( ispdf | isps | istif )
  
  
  if( isgif || isjpg || ispng || isbmp || issvg || isheic )
  {
    const bool decode_raw = (isgif || isjpg || ispng);
    
    const char *mimetype = "";
    if( isgif ) mimetype = "image/gif";
    else if( isjpg ) mimetype = "image/jpeg";
    else if( ispng ) mimetype = "image/png";
    else if( isbmp ) mimetype = "image/bmp";
    else if( issvg ) mimetype = "image/svg+xml";
    else if( isheic ) mimetype = "image/heic";
    
    new UploadedImgDisplay( m_viewer, displayName, mimetype, filesize, infile, dialog, type );
    
    add_undo_redo();
    
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
        Wt::WServer::instance()->ioService().boost::asio::io_service::post( std::bind( [=,this](){
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
      WText *errort = new WText( WString::tr("smm-invalid-peak-csv"), contents );
      errort->addStyleClass( "NonSpecError" );
      
      errort = new WText( string(e.what()), contents );
      errort->setAttributeValue( "style", "color: red; " );
      
      add_undo_redo();
      
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
      DrfSelect::createChooseDrfDialog( {det}, WString::tr("smm-file-is-drf").arg(name), "" );
      
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
  
  // Check if this is TSV/CSV file containing multiple DRFs
  if( header_contains( "Detector ID" )
     && handleGammaQuantDrfCsv(infile, displayName, fileLocation) )
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
  
#if( USE_REL_ACT_TOOL )
  if( currdata
     && header_contains( "<RelActCalcAuto " )
     && handleRelActAutoXmlFile(infile, dialog) )
  {
    return true;
  }
#endif
  
  // Check if a .ECC file from ISOCS
  if( (header_contains("SGI_template") || header_contains("ISOCS_file_name"))
     && handleEccFile(infile, dialog) )
  {
    add_undo_redo();
    
    return true;
  }//if( a .ECC file from ISOCS )
  
  if( header_contains("<ShieldingSourceFit") && header_contains("<Geometry")
     && (filesize > 128) && (filesize < 1024*1024)
     && handleShieldingSourceFile(infile, dialog) )
  {
    add_undo_redo();
    
    return true;
  }//if( Shielding/Source fit XML file )

  if( m_viewer->makeDrfWindow() )
  {
    //Check if source lib file.  These are text files, with each line defining a source,
    //  and looking like:
    //  "22NA_01551910  5.107E+04  25-Aug-2022 Some Remarks
    string datastr;
    datastr.resize( boost::size(data) + 1, '\0' );
    memcpy( &(datastr[0]), (const char *)&data[0], boost::size(data) );
    
    stringstream strm(datastr);
    if( SrcLibLineInfo::is_candidate_src_lib( strm )
       && handleSourceLibFile(infile, dialog) )
    {
      return true;
    }//if( candidate source lib file )
    
  }//if( m_viewer->makeDrfWindow() )

  
  delete dialog;
  
  return false;
}//void handleNonSpectrumFile(...)


void SpecMeasManager::closeNonSpecFileDialog()
{
  //assert( m_nonSpecFileDialog );
  if( !m_nonSpecFileDialog )
    return;

  const string dialog_id = m_nonSpecFileDialog->id();
  wApp->doJavaScript( "$('#" + dialog_id + "').hide(); $('.Wt-dialogcover').hide();" );
  
  m_nonSpecFileDialog->accept(); //This will delete the dialog, and cause `m_nonSpecFileDialog` to be set to nullptr
  m_nonSpecFileDialog = nullptr;
}//void closeNonSpecFileDialog()


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
      auto now = chrono::time_point_cast<chrono::microseconds>( chrono::system_clock::now() );
      now += chrono::seconds(60*offset);
      
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
      
      passMessage( WString::tr("smm-saved-drf-conf").arg(filename), WarningWidget::WarningMsgInfo );
    }catch( std::exception &e )
    {
      cerr << "handleMultipleDrfCsv: error saving multiple DRF file: " << e.what() << endl;
      passMessage( WString::tr("smm-err-saving-drf"), WarningWidget::WarningMsgHigh );
    }//try / catch to save file
  };//saveDrfFile
#endif
  
  WString dialogmsg = WString::tr( (drfs.size()==1) ? "smm-file-is-drf" : "smm-file-is-multi-drf");

  string creditsHtml;
  if( credits.size() )
  {
    for( const string &s : credits )
      creditsHtml += "<div>" + Wt::Utils::htmlEncode(s) + "</div>";
  }//if( credits.size() )
  
  
  DrfSelect::createChooseDrfDialog( drfs, dialogmsg, creditsHtml, saveDrfFile );
  
  return true;
}//bool handleMultipleDrfCsv( std::istream &input, SimpleDialog *dialog )


bool SpecMeasManager::handleGammaQuantDrfCsv( std::istream &input,
                           const std::string &displayName,
                           const std::string &fileLocation )
{
  vector<string> credits, warnings;
  vector<shared_ptr<DetectorPeakResponse>> drfs;
  
  try
  {
    DetectorPeakResponse::parseGammaQuantRelEffDrfCsv( input, drfs, credits, warnings );
  }catch( std::exception &e )
  {
    return false;
  }
  
  if( drfs.empty() )
    return false;
  
  // Print out as URL - for version of InterSpec pre 20241010
  //cout << "\n\n\n";
  //for( const auto &drf : drfs )
  //  cout << drf->name() << "    UrlEncoded  " << drf->toAppUrl() << endl;
  //cout << "-------- done ---------" << endl;
  
  
  std::function<void()> saveDrfFile;
  WString dialogmsg = WString::tr( (drfs.size()==1) ? "smm-file-is-drf" : "smm-file-is-multi-drf");

  string creditsHtml;
  if( credits.size() )
  {
    for( const string &s : credits )
    {
      if( !s.empty() )
        creditsHtml += "<div>" + Wt::Utils::htmlEncode(s) + "</div>";
    }
  }//if( credits.size() )
  
  if( !warnings.empty() )
  {
    if( !creditsHtml.empty() )
      creditsHtml += "<br />";
    for( const string &s : warnings )
    {
      if( !s.empty() )
        creditsHtml += "<div style=\"color: red\"><b>Warning</b>: " + Wt::Utils::htmlEncode(s) + "</div>";
    }
  }//if( !warnings.empty() )
  
  DrfSelect::createChooseDrfDialog( drfs, dialogmsg, creditsHtml, saveDrfFile );
  
  return true;
}//handleGammaQuantDrfCsv(...)


bool SpecMeasManager::handleCALpFile( std::istream &infile, SimpleDialog *dialog, bool autoApply )
{
  WGridLayout *stretcher = nullptr;
  WPushButton *closeButton = nullptr;
  
  // Make a lamda to clear dialog, if we are going to alter it
  auto clear_dialog = [&](){
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    closeButton = dialog->addButton( WString::tr("Close") );
    stretcher = new WGridLayout();
    
    // If we set the contents margins to 0, then scroll-bars may appear.
    //  However doing just the below looks okay, and the scroll bars dont seem to appear
    stretcher->setContentsMargins( 9, 2, 9, 2 );
    dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowHidden, Wt::Orientation::Horizontal );
    
    dialog->contents()->setLayout( stretcher );
    WText *title = new WText( WString::tr("smm-not-spec-file") );
    title->addStyleClass( "title" );
    stretcher->addWidget( title, 0, 0 );
  };//clear_dialog lamda
  
  
  auto currdata = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if( !currdata )
  {
    clear_dialog();
    assert( stretcher && closeButton );
    
    WText *t = new WText( WString::tr("smm-CALp-no-fore") );
    stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
    t->setTextAlignment( Wt::AlignCenter );
    
    return true;
  }//if( !currdata )
  
  set<string> fore_gamma_dets;
  const shared_ptr<SpecMeas> foreground = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  const set<int> &fore_samples = m_viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
  const vector<string> fore_dets = m_viewer->detectorsToDisplay( SpecUtils::SpectrumType::Foreground );
  
  if( !foreground )
    return false;
  
  const size_t num_display_channel = currdata->num_gamma_channels();
  map<string,shared_ptr<const SpecUtils::EnergyCalibration>> det_to_cal;
  
  const std::streampos start_pos = infile.tellg();
  
  while( infile.good() )
  {
    string name;
    
    // Note that num_display_channel is for the displayed data - the individual detectors may
    //  have different numbers of channels - we'll fix this up after initially loading the cal
    //  (we need to know the detectors name the cal is for in order to fix it up)
    shared_ptr<SpecUtils::EnergyCalibration> cal;
    
    try
    {
      cal = SpecUtils::energy_cal_from_CALp_file( infile, num_display_channel, name );
      assert( cal && cal->valid() );
    }catch( std::exception &e )
    {
      // Display message to user to let them know it was a CALp file, but we couldn't use it.
      //  TODO: improve this error message with details, ex if it was a lower-channel-energy CALp, and number of channels didnt match, we should display this to the user
      
      if( det_to_cal.empty() /* && SpecUtils::iends_with( displayName, "CALp" ) */ )
      {
        infile.seekg( start_pos, ios::beg );
        
        clear_dialog();
        assert( stretcher && closeButton );
        
        WText *t = new WText( WString::tr("smm-CALp-invalid") );
        stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
        t->setTextAlignment( Wt::AlignCenter );
        dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowVisible,
                                         Wt::Horizontal | Wt::Vertical );
        return true;
      }//if( we didnt get any calibrations )
      
      break;
    }//try / catch

    
    if( !cal->valid() )
    {
      assert( 0 );
      continue;
    }
    
    // The calibration may not be for the correct number of channels, for files that have detectors
    //  with different num channels; we'll check for this here and fix the calibration up for this
    //  case.
    //  We could also do this for `name.empty()` case, but we shouldnt need to, I dont think.
    if( !name.empty() )
    {
      for( const int sample_num : fore_samples )
      {
        const auto m = foreground->measurement(sample_num, name);
        const size_t nchan = m ? m->num_gamma_channels() : size_t(0);
        if( nchan > 3 )
        {
          if( nchan != num_display_channel )
          {
            try
            {
              switch( cal->type() )
              {
                case SpecUtils::EnergyCalType::Polynomial:
                case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
                  cal->set_polynomial( nchan, cal->coefficients(), cal->deviation_pairs() );
                  break;
                  
                case SpecUtils::EnergyCalType::FullRangeFraction:
                  cal->set_full_range_fraction( nchan, cal->coefficients(), cal->deviation_pairs() );
                  break;
                  
                case SpecUtils::EnergyCalType::LowerChannelEdge:
                case SpecUtils::EnergyCalType::InvalidEquationType:
                  break;
              }//switch( cal->type() )
            }catch( std::exception &e )
            {
              cerr << "Failed to change number of channels for energy calibration: " << e.what() << endl;
            }
          }//if( m->num_gamma_channels() != num_display_channel )
          
          break;
        }//if( nchan > 3 )
      }//for( const int sample_num : fore_samples )
    }//if( !name.empty() )
    
    det_to_cal[name] = cal;
  }//while( true )
  
  if( det_to_cal.empty() )
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
  
  WString msg = WString::tr("smm-CALp-contains-energy-cal");
  
  
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
          msg += SpecUtils::printCompact( cal->coefficients()[i], 4 );
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
    msg += WString::tr("smm-CALp-multi-dets").arg( static_cast<int>(det_to_cal.size()) );
  }
  
  if( !have_cal_for_all_dets )
    msg += WString::tr( (det_to_cal.size() == 1) ? "smm-warn-single-for-multi" : "smm-warn-multi-for-single" );
  
  WText *t = new WText( msg );
  stretcher->addWidget( t, stretcher->rowCount(), 0, AlignCenter | AlignMiddle );
  t->setTextAlignment( Wt::AlignCenter );
  
  WCheckBox *applyOnlyCurrentlyVisible = nullptr;
  WCheckBox *applyForeground = nullptr, *applyBackground = nullptr, *applySecondary = nullptr;
  
  if( fore_samples.size() != foreground->sample_numbers().size() )
  {
    applyOnlyCurrentlyVisible = new WCheckBox( WString::tr("smm-CALp-cb-disp-samples-only") );
    stretcher->addWidget( applyOnlyCurrentlyVisible, stretcher->rowCount(), 0, AlignLeft );
  }//if( not displaying all foreground samples )
  
  const auto back = m_viewer->measurment( SpecUtils::SpectrumType::Background );
  const auto second = m_viewer->measurment( SpecUtils::SpectrumType::SecondForeground );
  
  //Only have
  if( (back && (back != foreground)) || (second && (second != foreground)) )
  {
    applyForeground = new WCheckBox( WString::tr("smm-CALp-cb-apply-to-fore") );
    applyForeground->setChecked( true );
    stretcher->addWidget( applyForeground, stretcher->rowCount(), 0, AlignLeft );
  }
  
  if( back && (back != foreground) )
  {
    applyBackground = new WCheckBox( WString::tr("smm-CALp-cb-apply-to-back") );
    applyBackground->setChecked( true );
    stretcher->addWidget( applyBackground, stretcher->rowCount(), 0, AlignLeft );
  }
  
  if( second && (second != foreground) && (second != back) )
  {
    applySecondary = new WCheckBox( WString::tr("smm-CALp-cb-apply-to-second") );
    applySecondary->setChecked( true );
    stretcher->addWidget( applySecondary, stretcher->rowCount(), 0, AlignLeft );
  }
  
  
  t = new WText( WString::tr("smm-CALp-like-to-use") );
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
      const char *key = (napplied == 1) ? "smm-CALp-err-applying-single" : "mm-CALp-err-applying-mult";
      passMessage( WString::tr(key).arg( e.what() ), WarningWidget::WarningMsgHigh );
    }//try / catch
  };//applyLambda(...)
  
  
  if( autoApply && !applyOnlyCurrentlyVisible && !applyForeground && !applyBackground && !applySecondary )
  {
    applyLambda();
    dialog->done( Wt::WDialog::DialogCode::Accepted );
    return true;
  }
  
  dialog->addButton( WString::tr("No") ); //no further action necessary if user clicks no; dialog will close
  closeButton->setText( WString::tr("Yes") );
  closeButton->clicked().connect( std::bind( applyLambda ) );
  
  return true;
}//void handleCALpFile( std::istream &input )


#if( USE_REL_ACT_TOOL )
bool SpecMeasManager::handleRelActAutoXmlFile( std::istream &input, SimpleDialog *dialog )
{
  // blah blah blah add undo/redo support
  WString error_msg;
  try
  {
    input.seekg(0, ios::end);
    const size_t filesize = input.tellg();
    input.seekg(0);
    
    if( filesize > 10*1024*1024 )
      throw runtime_error( "Input file larger than expected." );
    
    vector<char> data( filesize + 1, 0x0 );
    
    if( !input.read( &(data[0]), filesize ) )
      throw runtime_error( "Failed to read file data." );
      
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_trim_whitespace>( &(data[0]) );
    
    
    RelActAutoGui *tool = m_viewer->showRelActAutoWindow();
    if( !tool )
      throw runtime_error( "Could not create <em>Isotopics by nuclide</em> tool." );
    
    tool->setGuiStateFromXml( &doc );
    
    dialog->done( Wt::WDialog::DialogCode::Accepted );
    return true;
  }catch( rapidxml::parse_error &e )
  {
    error_msg = WString::tr("smm-err-load-xml-iso-by-nuc").arg( e.what() );
    
    error_msg = "Error parsing config XML: " + string(e.what());
    const char * const position = e.where<char>();
    if( position && *position )
    {
      const char *end_pos = position;
      for( size_t i = 0; (*end_pos) && (i < 80); ++i )
        end_pos += 1;
      error_msg.arg( "<br />&nbsp;&nbsp;At: " + std::string(position, end_pos) );
    }else
    {
      error_msg.arg( "" );
    }//if( position ) / else
  }catch( std::exception &e )
  {
    error_msg = WString::tr("smm-err-load-iso-by-nuc").arg( e.what() ).toUTF8();
  }//try / cat to read the XML
  
  
  dialog->contents()->clear();
  dialog->footer()->clear();
    
  WPushButton *closeButton = dialog->addButton( WString::tr("Close") );
  WGridLayout *stretcher = new WGridLayout();
    
  // If we set the contents margins to 0, then scroll-bars may appear.
  //  However doing just the below looks okay, and the scroll bars dont seem to appear
  stretcher->setContentsMargins( 9, 2, 9, 2 );
  //dialog->contents()->setOverflow( WContainerWidget::Overflow::OverflowHidden );
    
  dialog->contents()->setLayout( stretcher );
  WText *title = new WText( WString::tr("smm-err-iso-by-nuc-window-title") );
  title->addStyleClass( "title" );
  stretcher->addWidget( title, 0, 0 );
  
  WText *content = new WText( error_msg );
  content->addStyleClass( "content" );
  stretcher->addWidget( content, 1, 0 );
  
  return true;
}//bool handleRelActAutoXmlFile( std::istream &input, SimpleDialog *dialog );
#endif


bool SpecMeasManager::handleEccFile( std::istream &input, SimpleDialog *dialog )
{
  const size_t start_pos = input.tellg();
  
  shared_ptr<DetectorPeakResponse> det;
  double source_area = 0.0, source_mass = 0.0;
  try
  {
    tuple<shared_ptr<DetectorPeakResponse>,double,double> det_area_mass
      = DetectorPeakResponse::parseEccFile( input );
    
    det = get<0>(det_area_mass);
    source_area = get<1>(det_area_mass);
    source_mass = get<2>(det_area_mass);
    
    assert( det && det->isValid() );
    if( !det || !det->isValid() )
      throw std::logic_error( "DRF returned from DetectorPeakResponse::parseEccFile() should be valid." );
  }catch( std::exception &e )
  {
    input.seekg( start_pos );
    return false;
  }//try / catch
  
  dialog->addStyleClass( "EccDrfDialog" );
  
  assert( dialog );
  
  dialog->contents()->clear();
  dialog->footer()->clear();
  
  int chartw = 350, charth = 200;
  if( m_viewer->renderedWidth() > 500 )
    chartw = std::min( ((3*m_viewer->renderedWidth()/4) - 50), 500 );
  if( m_viewer->renderedHeight() > 400 )
    charth = std::min( m_viewer->renderedHeight()/4, (4*chartw)/7 );
  chartw = std::max( chartw, 300 );
  charth = std::max( charth, 175 );
  
  WText *title = new WText( WString::tr("smm-ecc-curve"), dialog->contents() );
  title->addStyleClass( "title" );
  title->setInline( false );
  
  DrfChart *chart = new DrfChart( dialog->contents() );
  chart->setMinimumSize( 300, 175 );
  chart->resize( chartw, charth );
  chart->updateChart( det );
  
  auto set_chart_y_range = [=]( shared_ptr<DetectorPeakResponse> drf ){
    // We will override the auto y-axis limits, since it does badly with really small
    //  numbers we might encounter.
    double ymax = -999.0;
    double lower_x = drf->lowerEnergy();
    double upper_x = drf->upperEnergy();
    if( lower_x >= upper_x )
    {
      lower_x = 45;
      upper_x = 3000;
    }
    
    for( double energy = lower_x; energy <= upper_x; energy += 10 )
    {
      const double val = drf->intrinsicEfficiency( energy );
      ymax = std::max( ymax, val );
    }
    if( ymax > 0.0 )
      chart->axis(Chart::Y1Axis).setRange( 0.0, 1.2*ymax );
  };
  
  set_chart_y_range( det );
  
  const string name = Wt::Utils::htmlEncode( det->name() );
  const string desc = Wt::Utils::htmlEncode( det->description() );
    
  string txt_css = "style=\"text-align: left;"
  " max-width: " + std::to_string(chartw-5) + "px;"
  " white-space: nowrap;"
  " text-overflow: ellipsis;"
  " overflow-x: hidden;"
  "\"";
  
  string msg =
    //"<p style=\"white-space: nowrap;\">You can use this .ECC file as a DRF.</p>"
    "<p " + txt_css + ">"
      "Name: " + name +
  "</p>";
  if( !desc.empty() )
    msg += "<p " + txt_css + ">"
        "Desc: " + desc +
    "</p>";
  //msg += "<p>Would you like to use this DRF?</p>";
  
  WText *txt = new WText( msg, TextFormat::XHTMLText, dialog->contents() );

  WContainerWidget *btn_div = new WContainerWidget( dialog->contents() );
  btn_div->addStyleClass( "HowToUseGrp" );
  
  map<int,DetectorPeakResponse::EffGeometryType> index_to_geom;
  
  WLabel *geom_label = new WLabel( WString::tr("smm-ecc-how-to-interpret"), btn_div );
  WComboBox *geom_combo = new WComboBox( btn_div );
  geom_combo->addItem( WString::tr("smm-ecc-far-field") );
  index_to_geom[geom_combo->count() - 1] = DetectorPeakResponse::EffGeometryType::FarField;
  
  geom_combo->addItem( WString::tr("smm-ecc-fix-geom-total-act") );
  index_to_geom[geom_combo->count() - 1] = DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct;
  
  if( source_area > 0.0 )
  {
    geom_combo->addItem( WString::tr("smm-ecc-fix-geom-act-cm2") );
    index_to_geom[geom_combo->count() - 1] = DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2;
    
    geom_combo->addItem( WString::tr("smm-ecc-fix-geom-act-m2") );
    index_to_geom[geom_combo->count() - 1] = DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2;
  }//if( source_area > 0.0 )
  
  if( source_mass > 0 )
  {
    geom_combo->addItem( WString::tr("smm-ecc-fix-geom-act-gram") );
    index_to_geom[geom_combo->count() - 1] = DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram;
  }//if( source_mass > 0 )
  
  geom_combo->setCurrentIndex( 1 );
    
  WTable *far_field_opt = new WTable( dialog->contents() );
  //far_field_opt->setHiddenKeepsGeometry( true );
  far_field_opt->addStyleClass( "FarFieldOptTbl" );
  
  WRegExpValidator *dist_validator = new WRegExpValidator( PhysicalUnits::sm_distanceRegex, this );
  dist_validator->setFlags( Wt::MatchCaseInsensitive );
  dist_validator->setInvalidBlankText( "0.0 cm" );
  dist_validator->setMandatory( true );
    
  WTableCell *cell = far_field_opt->elementAt( 0, 0 );
  WLabel *label = new WLabel( WString::tr("smm-ecc-det-diam"), cell );
  cell = far_field_opt->elementAt( 0, 1 );
  WLineEdit *diameter_edit = new WLineEdit( "", cell );
  label->setBuddy( diameter_edit );
  diameter_edit->setValidator( dist_validator );
  diameter_edit->setEmptyText( "0 cm" );
  
  cell = far_field_opt->elementAt( 1, 0 );
  label = new WLabel( WString::tr("smm-ecc-dist"), cell );
  cell = far_field_opt->elementAt( 1, 1 );
  WLineEdit *distance_edit = new WLineEdit( "", cell );
  label->setBuddy( distance_edit );
  distance_edit->setValidator( dist_validator );
  distance_edit->setEmptyText( "0 cm" );
  
  // TODO: make option to correct for air-attenuation
  
  far_field_opt->hide();
  
  auto fore = InterSpec::instance()->measurment( SpecUtils::SpectrumType::Foreground );
  shared_ptr<DetectorPeakResponse> prev = fore ? fore->detector() : nullptr;
  
  // TODO: make option to make DRF default for detector model, or serial number
  
  dialog->addButton( WString::tr("Cancel") );
  WPushButton *accept = dialog->addButton( WString::tr("smm-ecc-use-drf") );
  
  
  auto try_create_farfield = [=]() -> shared_ptr<DetectorPeakResponse> {
    const double distance = PhysicalUnits::stringToDistance( distance_edit->text().toUTF8() );
    const double diameter = PhysicalUnits::stringToDistance( diameter_edit->text().toUTF8() );
      
    if( distance < 0.0 )
      throw runtime_error( "dist < 0" );
    if( diameter <= 0.0 )
      throw runtime_error( "diam <= 0" );
    
    const bool correct_for_air_atten = true;
    return det->convertFixedGeometryToFarField( diameter, distance, correct_for_air_atten );
  };//try_create_farfield
  
  auto update_state = [=](){
    const int index = geom_combo->currentIndex();
    const auto pos = index_to_geom.find(index);
    assert( pos != end(index_to_geom) );
    if( pos == end(index_to_geom) )
      throw logic_error( "SpecMeasManager::handleEccFile: unexpected index" );
    
    const DetectorPeakResponse::EffGeometryType geom_type = pos->second;
    
    try
    {
      shared_ptr<DetectorPeakResponse> new_drf = det;
      
      far_field_opt->setHidden( (geom_type != DetectorPeakResponse::EffGeometryType::FarField) );
      
      switch( geom_type )
      {
        case DetectorPeakResponse::EffGeometryType::FarField:
          new_drf = try_create_farfield();
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
          new_drf = det->convertFixedGeometryType( source_area, geom_type );
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
          new_drf = det->convertFixedGeometryType( source_mass, geom_type );
          break;
      }//switch( geom_type )
        
      chart->updateChart( new_drf );
      set_chart_y_range( new_drf );
      accept->enable();
    }catch( std::exception & )
    {
      chart->updateChart( nullptr );
      accept->disable();
    }
  };//update_state lambda
  
  geom_combo->activated().connect( std::bind(update_state) );
  distance_edit->textInput().connect( std::bind(update_state) );
  diameter_edit->textInput().connect( std::bind(update_state) );
  
  
  accept->clicked().connect( std::bind( [=](){
    const int index = geom_combo->currentIndex();
    const auto pos = index_to_geom.find(index);
    assert( pos != end(index_to_geom) );
    if( pos == end(index_to_geom) )
      throw logic_error( "SpecMeasManager::handleEccFile: unexpected index" );
    
    const DetectorPeakResponse::EffGeometryType geom_type = pos->second;
    
    auto new_drf = det;
    try
    {
      switch( geom_type )
      {
        case DetectorPeakResponse::EffGeometryType::FarField:
          new_drf = try_create_farfield();
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
          new_drf = det->convertFixedGeometryType( source_area, geom_type );
          break;
          
        case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
          new_drf = det->convertFixedGeometryType( source_mass, geom_type );
          break;
      }//switch( geom_type )
    }catch( std::exception &e )
    {
      passMessage( WString::tr("smm-ecc-error").arg(e.what()), WarningWidget::WarningMsgHigh );
      return;
    }//try / catch
    
    auto interspec = InterSpec::instance();
    if( !new_drf || !interspec )
      return;
      
    shared_ptr<DataBaseUtils::DbSession> sql = interspec->sql();
    const Wt::Dbo::ptr<InterSpecUser> &user = interspec->user();
    DrfSelect::updateLastUsedTimeOrAddToDb( new_drf, user.id(), sql );
    interspec->detectorChanged().emit( new_drf ); //This loads it to the foreground spectrum file
    
    UndoRedoManager *undoManager = InterSpec::instance()->undoRedoManager();
    if( undoManager && undoManager->canAddUndoRedoNow() )
    {
      auto undo = [prev](){
        InterSpec *viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( prev );
      };
      
      auto redo = [new_drf](){
        InterSpec *viewer = InterSpec::instance();
        if( viewer )
          viewer->detectorChanged().emit( new_drf );
      };
       
      // This next undo/redo wont bring up the dialog, but it will at least get us back
      //  to the original detector.
      undoManager->addUndoRedoStep( undo, redo, "Change to ECC DRF" );
    }
  }) );
    
  return true;
}//bool handleEccFile( std::istream &input, SimpleDialog *dialog )


bool SpecMeasManager::handleShieldingSourceFile( std::istream &input, SimpleDialog *dialog )
{
  const size_t start_pos = input.tellg();
  
  try
  {
    //get the filesize
    input.seekg(0, ios::end);
    const size_t end_pos = input.tellg();
    input.seekg(start_pos);
    const size_t file_size = end_pos - start_pos;
    if( (file_size < 128) || (file_size > 1024*1024) )
      throw runtime_error( "invalid size" );
    
    // We need to keep data around as long as the xml_document
    auto data = make_shared<vector<char>>( file_size + 1 );
    
    if( !input.read( (char *)(&((*data)[0])), file_size ) )
      throw runtime_error( "failed to read file" );
    
    (*data)[file_size] = '\0';
    
    auto xml_doc = make_shared<rapidxml::xml_document<char>>();
    const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
    xml_doc->parse<flags>( &((*data)[0]) );
    
  
    MaterialDB *material_db = m_viewer->materialDataBase();
    PeakModel *peak_model = m_viewer->peakModel();
    WSuggestionPopup *shield_suggest = m_viewer->shieldingSuggester();
      
    auto disp = make_unique<ShieldingSourceDisplay>( peak_model, m_viewer, shield_suggest, material_db );
    disp->deSerialize( xml_doc->first_node() );
    
    assert( dialog );
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WText *title = new WText( WString::tr("smm-act-shield-xml"), dialog->contents() );
    title->addStyleClass( "title" );
    title->setInline( false );
    
    WText *content = new WText( WString::tr("smm-act-shield-use"), dialog->contents() );
    content->addStyleClass( "content" );
    content->setInline( false );
    
    dialog->footer()->clear();
    dialog->addButton( WString::tr("Cancel") );
    WPushButton *btn = dialog->addButton( WString::tr("Yes") );
    btn->clicked().connect( std::bind([this,data,xml_doc](){
      InterSpec *viewer = InterSpec::instance();
      if( !viewer || !data || !xml_doc )
        return;
      
      try
      {
        ShieldingSourceDisplay *display = viewer->shieldingSourceFit();
        if( display )
          display->deSerialize( xml_doc->first_node() );
      }catch( std::exception &e )
      {
        passMessage( WString::tr("smm-err-act-shield").arg(e.what()), WarningWidget::WarningMsgHigh );
      }
    }) );
  }catch( std::exception &e )
  {
    input.seekg( start_pos );
    return false;
  }//try / catch
  
  return true;
}//bool handleShieldingSourceFile( std::istream &input, SimpleDialog *dialog );


bool SpecMeasManager::handleSourceLibFile( std::istream &input, SimpleDialog *dialog )
{
  MakeDrfWindow *tool = m_viewer->makeDrfWindow();
  if( !tool )
    return false;
    
  const size_t start_pos = input.tellg();
  
  try
  {
    const vector<SrcLibLineInfo> srcs = SrcLibLineInfo::sources_in_lib( input );
    
    if( srcs.empty() )
      throw runtime_error( "No sources in file." );
    
    vector<shared_ptr<const SrcLibLineInfo>> src_ptrs;
    for( const SrcLibLineInfo &info : srcs )
      src_ptrs.push_back( make_shared<SrcLibLineInfo>(info) );
    
    assert( dialog );
    dialog->contents()->clear();
    dialog->footer()->clear();
    
    WText *title = new WText( WString::tr("smm-source-lib-title"), dialog->contents() );
    title->addStyleClass( "title" );
    title->setInline( false );
    
    WText *content = new WText( WString::tr("smm-source-source-use"), dialog->contents() );
    content->addStyleClass( "content" );
    content->setInline( false );
    
    dialog->footer()->clear();
    
    auto setter = [src_ptrs]( const bool autopopulate ){
      InterSpec *viewer = InterSpec::instance();
      MakeDrfWindow *tool = viewer ? viewer->makeDrfWindow() : nullptr;
      assert( tool );
      if( tool )
        tool->tool()->useSourceLibrary( src_ptrs, autopopulate );
    };
    
    WPushButton *btn = dialog->addButton( WString::tr("Yes") );
    btn->clicked().connect( std::bind([setter](){ setter(true); }) );
    
    btn = dialog->addButton( WString::tr("No") );
    btn->clicked().connect( std::bind([setter](){ setter(false); }) );
    
    dialog->addButton( WString::tr("Cancel") );
  }catch( std::exception &e )
  {
    input.seekg( start_pos );
    return false;
  }//try caltch
  
  return true;
}//bool handleSourceLibFile( std::istream &input, SimpleDialog *dialog )


void SpecMeasManager::handleCancelPreviousStatesDialog( AuxWindow *dialog )
{
  assert( dialog == m_previousStatesDialog );
  if( !dialog )
    return;
  
  if( dialog != m_previousStatesDialog )
  {
    cerr << "SpecMeasManager::handleCancelPreviousStatesDialog: dialog passed in isnt as expected"
    << " - not doing anything." << endl;
    return;
  }
  
  m_previousStatesDialog = nullptr;
  AuxWindow::deleteAuxWindow( dialog );
}//void handleCancelPreviousStatesDialog( AuxWindow *dialog )


void SpecMeasManager::handleClosePreviousStatesDialogAfterSelect( AuxWindow *dialog )
{
  handleCancelPreviousStatesDialog( dialog );
}//void handleClosePreviousStatesDialogAfterSelect( AuxWindow *dialog )


void SpecMeasManager::checkCloseUploadDialog( SimpleDialog *dialog, WApplication *app )
{
  WApplication::UpdateLock lock( app );
  assert( lock );
  if( !lock )
    return;
  
  if( dialog != m_processingUploadDialog )
    return;
  
  if( m_processingUploadTimer )
  {
    boost::system::error_code ec;
    m_processingUploadTimer->cancel( ec );
    if( ec )
      cerr << "SpecMeasManager::checkCloseUploadDialog(): error cancelling timer: " << ec.message() << endl;
    m_processingUploadTimer.reset();
  }//if( m_processingUploadTimer )
  
  if( m_processingUploadDialog )
  {
    m_processingUploadDialog->done(WDialog::DialogCode::Accepted);
    m_processingUploadDialog = nullptr;
  }//if( m_processingUploadDialog )
  
  passMessage( WString::tr("smm-upload-timeout"), WarningWidget::WarningMsgHigh ) ;
}//void checkCloseUploadDialog( SimpleDialog *dialog )


void SpecMeasManager::handleDataRecievedStatus( uint64_t num_bytes_recieved, uint64_t num_bytes_total,
                                               SpecUtils::SpectrumType type )
{
  if( num_bytes_total < sm_minNumBytesShowUploadProgressDialog )
    return;
  
  //cout << "SpecMeasManager::handleDataRecievedStatus: " << num_bytes_recieved << " of " << num_bytes_total << " total." << endl;
  
  auto app = WApplication::instance();
  assert( app );
  if( !app )
    return;
  
  WApplication::UpdateLock lock( app );
  assert( lock );
  if( !lock )
    return;
  
  auto make_timer = [this,app]( SimpleDialog *dialog ){
    Wt::WServer *server = Wt::WServer::instance();
    if( !server )
      return;
    
    if( m_processingUploadTimer )
    {
      boost::system::error_code ec;
      m_processingUploadTimer->cancel( ec );
      if( ec )
        cerr << "SpecMeasManager::handleDataRecievedStatus(): error cancelling timer: " << ec.message() << endl;
      m_processingUploadTimer.reset();
    }//if( m_processingUploadTimer )
    
    auto timeoutfcn = app->bind( boost::bind(&SpecMeasManager::checkCloseUploadDialog, this, dialog, app) );
    m_processingUploadTimer = make_unique<boost::asio::deadline_timer>( server->ioService() );
    m_processingUploadTimer->expires_from_now( boost::posix_time::seconds(120) );
    m_processingUploadTimer->async_wait( [timeoutfcn](const boost::system::error_code &ec){ if(!ec) timeoutfcn(); } );
  };//make_timer lamda
  
  
  if( m_processingUploadDialog )
  {
    assert( m_processingUploadTimer );
    
    // TODO: could occasionally (every few seconds) update dialog text with current status
    
    make_timer( m_processingUploadDialog );
    return;
  }//if( m_processingUploadDialog )
  
  assert( !m_processingUploadTimer );
  
#if( BUILD_AS_LOCAL_SERVER || BUILD_FOR_WEB_DEPLOYMENT )
  const WString title = WString::tr("smm-finishing-upload");
#else
  const WString title = WString::tr("smm-finishing-copying");
#endif
  
  WString msg = WString::tr("smm-finish-up-txt");
  {//begin codeblock to add file size to string
    char filesize_str[64] = { '\0' };
    if( num_bytes_total < 1024*1024 )
      snprintf( filesize_str, sizeof(filesize_str), "%.1f kb", num_bytes_total/1024.0 );
    else
      snprintf( filesize_str, sizeof(filesize_str), "%.1f Mb", num_bytes_total/(1024.0*1024.0) );
    
    msg.arg( filesize_str );
    
    switch( type )
    {
      case SpecUtils::SpectrumType::Foreground:
        msg.arg( WString::tr("foreground") );
        break;
      case SpecUtils::SpectrumType::SecondForeground:
        msg.arg( WString::tr("secondary") );
        break;
      case SpecUtils::SpectrumType::Background:
        msg.arg( WString::tr("background") );
        break;
    }//switch( type )
  }//end codeblock to add file size to string
   
  m_processingUploadDialog = new SimpleDialog( title, msg );
  
  make_timer( m_processingUploadDialog );
  
  app->triggerUpdate();
}//void handleDataRecievedStatus(...)


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
  
 
  if( (name.length() > 4)
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

    displayFile( modelRow, measurement, type, true, true, 
              SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
    
    //It is the responsibility of the caller to clean up the file.
  }catch( exception &e )
  {
    if( !handleNonSpectrumFile( name, spoolName, type ) )
    {
      displayInvalidFileMsg( name, e.what() );
    }
  }
  
}//handleFileDropWorker(...)


void SpecMeasManager::handleFileDrop( const std::string &name,
                                             const std::string &spoolName,
                                             SpecUtils::SpectrumType type )
{
  if( m_processingUploadTimer )
  {
    boost::system::error_code ec;
    m_processingUploadTimer->cancel( ec );
    if( ec )
      cerr << "SpecMeasManager::handleFileDrop(): error cancelling timer: " << ec.message() << endl;
    m_processingUploadTimer.reset();
  }//if( m_processingUploadTimer )
  
  if( m_processingUploadDialog )
  {
    m_processingUploadDialog->done(WDialog::DialogCode::Accepted);
    m_processingUploadDialog = nullptr;
  }//if( m_processingUploadDialog )
  
  if( m_previousStatesDialog )
    handleCancelPreviousStatesDialog( m_previousStatesDialog );
  
  // If file is small, and not csv/txt (these are really slow to parse), dont display the parsing
  //  message.
  if( (SpecUtils::file_size(spoolName) < 512*1024)
     && !SpecUtils::iends_with(name, ".csv") && !SpecUtils::iends_with(name, ".txt") )
  {
    handleFileDropWorker( name, spoolName, type, nullptr, wApp );
    return;
  }
  
  // Its a larger file - display a message letting the user know its being parsed.
  auto dialog = new SimpleDialog( WString::tr("smm-window-title-parsing"),
                                 WString::tr("smm-window-msg-parsing") );
  
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


#if( USE_QR_CODES )
void SpecMeasManager::handleSpectrumUrl( std::string &&unencoded )
{
  try
  {
    if( m_previousStatesDialog )
      handleCancelPreviousStatesDialog( m_previousStatesDialog );
    
    //Remove everything leading up to "RADDATA://G0/"
    const size_t uri_pos = SpecUtils::ifind_substr_ascii( unencoded, "RADDATA://" );
    if( (uri_pos != string::npos) && (uri_pos > 0) )
      unencoded = unencoded.substr( 0, uri_pos );
    
    
    SpecUtils::EncodedSpectraInfo info;
    
    try
    {
      info = SpecUtils::get_spectrum_url_info( unencoded );
    }catch( std::exception &e )
    {
      // We'll try one extra URL decode
      string second_chance = SpecUtils::url_decode( unencoded );
      try
      {
        info = SpecUtils::get_spectrum_url_info( second_chance );
        unencoded = second_chance;
      }catch( std::exception & )
      {
        // No go - give up
        throw runtime_error( e.what() );
      }
    }//try / catch, SpecUtils::get_spectrum_url_info
    
    
    if( info.m_number_urls == 1 )
    {
      if( m_multiUrlSpectrumDialog )
      {
        m_multiUrlSpectrumDialog->accept();
        multiSpectrumDialogDone();
      }
      
      vector<SpecUtils::UrlSpectrum> spectra = SpecUtils::spectrum_decode_first_url( unencoded );
      if( spectra.empty() )
        throw runtime_error( WString::tr("smm-qr-no-meas").toUTF8() );
      
      shared_ptr<SpecUtils::SpecFile> specfile = SpecUtils::to_spec_file( spectra );
      assert( specfile && specfile->num_measurements() );
      
      if( !specfile || (specfile->num_measurements() < 1) || (specfile->num_gamma_channels() < 7) )
        throw runtime_error( WString::tr("smm-qr-no-gamma").toUTF8() );
      
      string display_name = spectra[0].m_title;
      if( display_name.size() < 3 )
        display_name = WString::tr("smm-qr-name-placeholder").toUTF8();

      auto specmeas = make_shared<SpecMeas>();
      reinterpret_cast<SpecUtils::SpecFile &>( *specmeas ) = *specfile;
      
      m_viewer->userOpenFile( specmeas, display_name );
    }else
    {
      auto dialog = dynamic_cast<MultiUrlSpectrumDialog *>( m_multiUrlSpectrumDialog );
      if( !dialog )
        m_multiUrlSpectrumDialog = dialog = new MultiUrlSpectrumDialog( this, m_viewer );
      dialog->addUrl( unencoded );
    }
  }catch( std::exception &e )
  {
    auto dialog = new SimpleDialog( WString::tr("Error"), 
                                   WString::tr("smm-qr-err-decode").arg(e.what()) );
    dialog->addButton( WString::tr("Close") );
  }//try /catch
}//void handleSpectrumUrl( const std::string &url );


void SpecMeasManager::displaySpectrumQrCode( const SpecUtils::SpectrumType type )
{
  const shared_ptr<SpecMeas> meas = m_viewer->measurment( type );
  if( !meas )
  {
    WString msg = WString::tr("smm-no-file-disp")
                    .arg( WString::tr(SpecUtils::descriptionText(type)) );
    passMessage( msg, WarningWidget::WarningMsgHigh ) ;
    return;
  }//if( !spec )
  
  shared_ptr<const SpecUtils::Measurement> spec;
  const set<int> &sample_nums = m_viewer->displayedSamples(type);
  const vector<string> detectors = m_viewer->detectorsToDisplay(type);
  
  // We `m_viewer->displayedHistogram( type )` may have some modifications we'll avoid if we can
  if( (sample_nums.size() == 1) && (detectors.size() == 1 ) )
    spec = meas->measurement( *begin(sample_nums), detectors.front() );
  
  if( !spec )
  {
    try
    {
      spec = meas->sum_measurements( sample_nums, detectors, nullptr );
    }catch( std::exception &e )
    {
      cerr << "SpecMeasManager::displaySpectrumQrCode: failed to sum measurements: "
           << e.what() << endl;
    }
  }//if( !spec )
  
  
  if( !spec )
    spec = m_viewer->displayedHistogram( type );
  
  if( !spec || (spec->num_gamma_channels() < 1) )
  {
    WString msg = WString::tr("smm-no-spec-disp")
                    .arg( WString::tr(SpecUtils::descriptionText(type)) );
    passMessage( msg, WarningWidget::WarningMsgHigh ) ;
    return;
  }//if( !spec )
  
  try
  {
    string model;
    if( meas->detector_type() != SpecUtils::DetectorType::Unknown )
      model = detectorTypeToString( meas->detector_type() );
    if( model.empty() )
      model = meas->instrument_model();
  
    vector<SpecUtils::UrlSpectrum> urlspec = SpecUtils::to_url_spectra( {spec}, model );
    
    const uint8_t encode_options = 0;
    vector<QRSpectrum::QrCodeEncodedSpec> urls;
    QRSpectrum::QrErrorCorrection ecc = QRSpectrum::QrErrorCorrection::High;
    
    // We'll try to only use a single QR-code, but also use highest error-correction we can.
    //  I'm sure this trade-off can be better handled in some way.
    if( spec->num_gamma_channels() < 5000 )
    {
      try
      {
        urls = QRSpectrum::qr_code_encode_spectra( urlspec, ecc, encode_options );
      }catch( std::exception & )
      {
      }//try / catch
    }//if( spec->num_gamma_channels() < 5000 )
    
    if( urls.empty() || (urls.size() > 1) )
    {
      try
      {
        ecc = QRSpectrum::QrErrorCorrection::Medium;
        urls = QRSpectrum::qr_code_encode_spectra( urlspec, ecc, encode_options );
      }catch( std::exception & )
      {
      }//try / catch
    }//if( urls.empty() )
    
    if( urls.empty() || (urls.size() > 1) )
    {
      try
      {
        ecc = QRSpectrum::QrErrorCorrection::Low;
        urls = QRSpectrum::qr_code_encode_spectra( urlspec, ecc, encode_options );
      }catch( std::exception & )
      {
      }//try / catch
    }//if( urls.empty() )
    
    if( urls.empty() )
    {
      auto dialog = new SimpleDialog( WString::tr("Error"), WString::tr("smm-qr-couldnt-encode") );
      dialog->addButton( WString::tr("Close") );
      return;
    }//if( urls.empty() )
    
    displayQrDialog( urls, 0, type );
  }catch( std::exception &e )
  {
    auto dialog = new SimpleDialog( WString::tr("Error"),
                                   WString::tr("smm-qr-encode-fail").arg(e.what()) );
    dialog->addButton( WString::tr("Close") );
  }//try catch
}//void displaySpectrumQrCode( const SpecUtils::SpectrumType type )


void SpecMeasManager::multiSpectrumDialogDone()
{
  m_multiUrlSpectrumDialog = nullptr;
}//void SpecMeasManager::multiSpectrumDialogDone()
#endif


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
  
  WString msg = WString::tr("smm-err-parse-spec").arg(lastpart).arg(errormsg);
  
  SimpleDialog *dialog = new SimpleDialog( WString::tr("smm-err-parse-spec-title"), msg );
  dialog->addButton( WString::tr("Close") );
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
    return ret;
  
  for( WModelIndexSet::iterator iter = selected.begin(); iter!=selected.end(); ++iter )
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
    }// switch( level(index) )
    ret.push_back(header);
  }//for( loop over selected )
  
  return ret;
}//std::shared_ptr<SpectraFileHeader> SpecMeasManager::selectedFile() const


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
      passMessage( "Currently cant load samples from multiple files, sorry, will is lazy.", 3 );
      return std::shared_ptr<SpectraFileHeader>();
    }
    
  }//for( const WModelIndex &index : selected )

  return header;
}//shared_ptr<SpectraFileHeader> SpecMeasManager::selectedFile() const


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
  auto window = new FileUploadDialog( m_viewer, this );
  UndoRedoManager *undoRedo = m_viewer->undoRedoManager();
  
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto closer = wApp->bind( boost::bind(&AuxWindow::hide, window) );
    // We wont use a redo step, because this would be a bit weird to redo.
    undoRedo->addUndoRedoStep( closer, nullptr, "Show file upload dialog" );
  }
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

  displayFile( row, measement_ptr, type, true, true, 
              SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
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
                "this bug to interspec@sandia.gov", WarningWidget::WarningMsgHigh );
    return;
  }//if( !index.isValid() )
  
  
  if( binning != "Keep All" )
  {
    try
    {
      meas->keep_energy_cal_variants( {binning} );
    }catch( std::exception &e )
    {
      // not expected
      passMessage( "There was an error separating energy cal type; loading all (error: "
                   + std::string(e.what()) + ")", WarningWidget::WarningMsgHigh )
    }
    
    // Trigger a refresh of row info and selected rows in File Manager
    m_fileModel->removeRows( index.row(), 1 );
    header->setMeasurmentInfo( meas );
    m_fileModel->addRow( header );
    index = m_fileModel->index( header );
  }//if( binning != "Keep All" ) / else
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
              SpecMeasManager::VariantChecksToDo::MultiVirtualDets );
}//SpecMeasManager::selectEnergyBinning(...)


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
                "this bug to interspec@sandia.gov", WarningWidget::WarningMsgHigh );
    return;
  }//if( !index.isValid() )
  
  VariantChecksToDo furtherChecks = VariantChecksToDo::None;
  switch( tokeep )
  {
    case DerivedDataToKeep::All:
      furtherChecks = VariantChecksToDo::MultiEnergyCalsAndMultiVirtualDets;
      break;
      
    case DerivedDataToKeep::RawOnly:
      furtherChecks = VariantChecksToDo::MultiEnergyCalsAndMultiVirtualDets;
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
      
      if( (meas->detector_type() == SpecUtils::DetectorType::VerifinderNaI)
         || (meas->detector_type() == SpecUtils::DetectorType::VerifinderLaBr)
         || SpecUtils::icontains(meas->manufacturer(), "Symetrica") )
      {
        // Remove the calibration and stabilization measurements; user probably doesnt want these.
        vector<shared_ptr<const SpecUtils::Measurement>> meas_to_remove;
        const vector<shared_ptr<const SpecUtils::Measurement>> orig_meas = meas->measurements();
        for( const shared_ptr<const SpecUtils::Measurement> &m : orig_meas )
        {
          if( SpecUtils::icontains( m->title(), "StabMeas" ) 
             || SpecUtils::icontains( m->title(), "Gamma Cal" )
             || SpecUtils::icontains( m->title(), "CalMeasurement" ) )
          {
            meas_to_remove.push_back( m );
          }
        }//for( loop over remaining measurements )
        
        if( !meas_to_remove.empty() )
        {
          meas->set_uuid( "" );
          meas->remove_measurements( meas_to_remove );
        }
        
        // Also, a lot of times there are multiple energy calibrations, like "ECalVirtual3keV"
        //  and/or something like "ECalGamma-SG_xxxxx-xxxxx", so we'll have detectors named
        //  "DetectorInfoGamma", "DetectorInfoGamma_intercal_ECalVirtual3keV", etc.
        //  So lets also rename these detectors.
        const set<string> energy_cal_variants = meas->energy_cal_variants();
        if( !energy_cal_variants.empty() )
        {
          meas->set_uuid( "" );
          meas->keep_energy_cal_variants( energy_cal_variants );
        }
      }//if( a Symetrica detector )
      
      // Trigger a refresh of row info and selected rows in File Manager
      m_fileModel->removeRows( index.row(), 1 );
      header->setMeasurmentInfo( meas );
      m_fileModel->addRow( header );
      index = m_fileModel->index( header );
      break;
    }//case RawOnly or DerivedOnly
  }//switch( tokeep )
  
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
              furtherChecks );
}//void selectDerivedDataChoice(...)



void SpecMeasManager::selectVirtualDetectorChoice( const std::set<std::string> tokeep,
                         std::shared_ptr<SpectraFileHeader> header,
                         std::shared_ptr<SpecMeas> meas,
                         const SpecUtils::SpectrumType type,
                         const bool checkIfPreviouslyOpened,
                         const bool doPreviousEnergyRangeCheck )
{
  WModelIndex index = m_fileModel->index( header );
  
  if( !index.isValid() || !meas )
  {
    passMessage( "Aborting loading of file after selecting virtual detector type to keep - "
                "the file is no longer available in memory. please report "
                "this bug to interspec@sandia.gov", WarningWidget::WarningMsgHigh );
    return;
  }//if( !index.isValid() )
  
  if( !tokeep.empty() )
  {
    set<string> dets_to_remove;
    for( const string &name : meas->gamma_detector_names() )
    {
      if( SpecUtils::istarts_with(name, "VD") && !tokeep.count(name) )
        dets_to_remove.insert( name );
    }
    
    meas->remove_detectors_data( dets_to_remove );
    
    // Trigger a refresh of row info and selected rows in File Manager
    m_fileModel->removeRows( index.row(), 1 );
    header->setMeasurmentInfo( meas );
    m_fileModel->addRow( header );
    index = m_fileModel->index( header );
  }//if( !tokeep.empty() )
  
  const VariantChecksToDo furtherChecks = VariantChecksToDo::None;
  
  displayFile( index.row(), meas, type, checkIfPreviouslyOpened, doPreviousEnergyRangeCheck,
              furtherChecks );
}//void selectVirtualDetectorChoice(...)


bool SpecMeasManager::checkForAndPromptUserForDisplayOptions( std::shared_ptr<SpectraFileHeader> header,
                                            std::shared_ptr<SpecMeas> meas,
                                            const SpecUtils::SpectrumType type,
                                            const bool checkIfPreviouslyOpened,
                                            const bool doPreviousEnergyRangeCheck,
                                            VariantChecksToDo viewingChecks )
{
  if( !header || !meas )
    throw runtime_error( "SpecMeasManager::checkForAndPromptUserForDisplayOptions(): Invalid input" );
  
  if( viewingChecks == VariantChecksToDo::None )
    return false;
  
  // We will first deal with "Derived Data" being present, then "Multiple Energy Calibrations",
  //  and then finally "Multiple Virtual Detectors"
  
  
  bool derivedData = false, energyCal = false;
  if( viewingChecks == VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets )
    derivedData = (meas->contains_derived_data() && meas->contains_non_derived_data());
  
  set<string> cals;
  
  if( !derivedData )
  {
    cals = meas->energy_cal_variants();
    energyCal = (cals.size() > 1);
  }
  
  if( (!derivedData && !energyCal && (viewingChecks == VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets))
     || (!energyCal && (viewingChecks == VariantChecksToDo::MultiEnergyCalsAndMultiVirtualDets))
     || (viewingChecks == VariantChecksToDo::MultiVirtualDets) )
  {
    vector<string> virtual_detectors;
    for( const string &name : meas->gamma_detector_names() )
    {
      if( SpecUtils::istarts_with( name, "VD" ) )
        virtual_detectors.push_back( name );
    }
    
    if( virtual_detectors.size() > 1 )
    {
      WString msgtxt = WString::tr("smm-vd-load-msg");
      if( virtual_detectors.size() > 5 )
        msgtxt.arg( WString::tr("smm-vd-gt-5dets").arg( static_cast<int>(virtual_detectors.size()) ) );
      else
        msgtxt.arg( "" );
      
      SimpleDialog *dialog = new SimpleDialog( WString::tr("smm-vd-load-title"), msgtxt );
      
      auto add_button = [=,this]( string btn_txt, const set<string> &dets, size_t max_txt_size ){
        if( btn_txt.size() > max_txt_size )
        {
          SpecUtils::utf8_limit_str_size( btn_txt, max_txt_size - 1 );
          btn_txt += "...";
        }
        btn_txt = Wt::Utils::htmlEncode( btn_txt );
          
        WPushButton *button = dialog->addButton( btn_txt );
        button->clicked().connect( boost::bind( &SpecMeasManager::selectVirtualDetectorChoice, this,
                                               dets, header, meas, type,
                                               checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
      };//add_button
      
      add_button( WString::tr("smm-vd-all-btn").toUTF8(), set<string>{}, 6 );
      
      for( size_t i = 0; (i < 5) && (i < virtual_detectors.size()); ++i )
      {
        add_button( virtual_detectors[i], set<string>{virtual_detectors[i]}, 6 );
      }//for( size_t i = 0; (i < 5) && (i < virtual_detectors.size()); ++i )
      
      // If three detectors, let them choose any 1, any 2, or all
      if( virtual_detectors.size() == 3 )
      {
        string btn_txt = virtual_detectors[0] + "+" + virtual_detectors[1];
        add_button( btn_txt, {virtual_detectors[0], virtual_detectors[1]}, 9 );
        
        btn_txt = virtual_detectors[1] + "+" + virtual_detectors[2];
        add_button( btn_txt, {virtual_detectors[1], virtual_detectors[2]}, 9 );
      }//if( virtual_detectors.size() <= 3 )
      
      return true;
    }//if( virtual_detectors.size() > 1 )
  }//if( viewingChecks == VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets )
  
  
  if( !derivedData && !energyCal )
    return false;

  
  if( derivedData )
  {
    SimpleDialog *dialog = new SimpleDialog( WString::tr("smm-derived-window-title"),
                                            WString::tr("smm-derived-window-txt") );
    WPushButton *button = dialog->addButton( WString::tr("smm-derived-all") );
    
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::All, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    button = dialog->addButton( WString::tr("smm-derived-raw") );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::RawOnly, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    button = dialog->addButton( WString::tr("smm-derived-derived") );
    button->clicked().connect( boost::bind( &SpecMeasManager::selectDerivedDataChoice, this,
                                           DerivedDataToKeep::DerivedOnly, header, meas, type,
                                           checkIfPreviouslyOpened, doPreviousEnergyRangeCheck ) );
    
    return true;
  }//if( derivedData )

  // Everything past here is for selecting energy binning - lets check to make sure we are supposed
  //  to check for this
  if( viewingChecks == VariantChecksToDo::MultiVirtualDets )
  {
    return false;
  }
  
  const char *title_key = "smm-multiple-binnings-window-title";
  const char *msgtxt_key = "smm-multiple-binnings-window-txt";
  
  if( cals.size() > 3 )
  {
    SimpleDialog *dialog = new SimpleDialog( WString::tr(title_key) );
    
    int ncolwide = (derivedData ? 3 : static_cast<int>(cals.size() + 1));
    if( ncolwide > 4 )
      ncolwide = 4;
  
    WTable *table = new WTable( dialog->contents() );
  
    WText *msg = new WText( WString::tr(msgtxt_key), XHTMLText );
    //layout->addWidget( msg, 0, 0, 1, ncolwide );
    WTableCell *cell = table->elementAt( 0, 0 );
    cell->addWidget( msg );
    cell->setColumnSpan( ncolwide );
    
    WPushButton *button = new WPushButton( WString::tr("smm-multiple-binnings-keep-all-btn") );
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
    SimpleDialog *dialog = new SimpleDialog( WString::tr(title_key) );
    
    WText *msg = new WText( WString::tr(msgtxt_key), XHTMLText, dialog->contents() );
    msg->addStyleClass( "content" );
    
    WPushButton *button = dialog->addButton( WString::tr("smm-multiple-binnings-keep-all-btn") );
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
      = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
#endif
  
  
  if( row < 0 && !measement_ptr )
  {
    WFlags<InterSpec::SetSpectrumOptions> options;
    options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
    options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
    
    m_viewer->setSpectrum( measement_ptr, std::set<int>(), type, options );

    if( old_meas.use_count() == 1 )
      serializeToTempFile( old_meas );
    
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
          && !oldheader->m_fileDbEntry->userHasModified )
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
    {
      return;
    }
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
          warningmsg << WString::tr("smm-warn-has-intrinsic").toUTF8();
        }else
        {
          warningmsg << WString::tr("smm-warn-upload-mult-spec").arg(nsamples).toUTF8();
          
          if( numForeground )
            warningmsg << WString::tr("smm-warn-upload-mult-showing-fore").toUTF8();
          else if( childrow == 0 )
            warningmsg << WString::tr("smm-warn-upload-mult-showing-first").toUTF8();
          else
            warningmsg << WString::tr("smm-warn-upload-mult-showing-other")
                            .arg(header->m_samples[childrow].sample_number)
                            .toUTF8();
          
          if( m_viewer->toolTabsVisible() )
            warningmsg << WString::tr("smm-warn-upload-use-spec-manager").toUTF8();
          else
            warningmsg << WString::tr("smm-warn-upload-use-file-manager").toUTF8();
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
          warningmsg << WString::tr("smm-multi-cal-present")
            .arg( ncalibration )
            .arg( WString::tr( (ncalibration == 1) ? "smm-multi-cal-it" : "smm-multi-cal-they" ) )
            .toUTF8();
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
            warningmsg << WString::tr("smm-multi-spectra-using-back").toUTF8();
            selected.clear();
            selected.insert( index.child(sample,0) );
            break;
          }//if( foundBackground )
        } // for( int sample = 0; sample < nsamples; ++sample )

        if( !foundBackground )
        {
          warningmsgLevel = max(WarningWidget::WarningMsgHigh, warningmsgLevel );
          warningmsg << WString::tr("smm-multi-spectra-using-first").arg(nspectra_header).toUTF8();
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
    passMessage( warningmsg.str(), warningmsgLevel );
  
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

//#if( USE_DB_TO_STORE_SPECTRA )
//  if( storeInDb )
//    saveToDatabase( old_meas );
//#endif

  if( background_sample_numbers.size() )
  {
    WFlags<InterSpec::SetSpectrumOptions> options;
    if( doPreviousEnergyRangeCheck )
    {
      options |= InterSpec::SetSpectrumOptions::CheckToPreservePreviousEnergyCal;
      options |= InterSpec::SetSpectrumOptions::CheckForRiidResults;
      options |= InterSpec::SetSpectrumOptions::SkipParseWarnings;
#if( USE_REMOTE_RID )
      options |= InterSpec::SetSpectrumOptions::SkipExternalRid;
#endif
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
#if( USE_REL_ACT_TOOL )
    newspec->setRelActManualGuiState( nullptr );
    newspec->setRelActAutoGuiState( nullptr );
#endif
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
    passMessage( WString::tr("smm-failed-combine").arg(e.what()), WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, WString::tr("smm-failed-combine").arg(e.what()).toUTF8().c_str() );
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
    WString msg = WString::tr("smm-failed-sum").arg(e.what());
    passMessage( msg, WarningWidget::WarningMsgHigh );
    
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.toUTF8().c_str() );
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
  if( !m_spectrumManagerWindow )
    return;
  
  m_spectrumManagertreeDiv->removeWidget(m_treeView);
  AuxWindow::deleteAuxWindow(m_spectrumManagerWindow);
  m_spectrumManagerWindow = nullptr;
  
  UndoRedoManager *undoRedo = m_viewer->undoRedoManager();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [](){
      InterSpec *viewer = InterSpec::instance();
      SpecMeasManager *manager = viewer->fileManager();
      if( manager )
        manager->startSpectrumManager();
    };
    auto redo = [](){
      InterSpec *viewer = InterSpec::instance();
      SpecMeasManager *manager = viewer->fileManager();
      if( manager )
        manager->deleteSpectrumManager();
    };
    
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Spectrum Manager.");
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
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
  new WText( WString::tr("smm-selected-spec-label"), m_newDiv );
  
  m_setButton = new WPushButton( WString::tr("smm-disp-as-btn"), m_newDiv);
  m_setButton->setIcon( "InterSpec_resources/images/bullet_arrow_down.png" );
  
  WPopupMenu *setPopup = new WPopupMenu();
  m_setAsForeground = setPopup->addItem( WString::tr("Foreground") );
  m_setAsForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::Foreground, true) );
  m_setAsBackground = setPopup->addItem( WString::tr("Background") );
  m_setAsBackground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::Background, true) );
  m_setAsSecForeground = setPopup->addItem( WString::tr("second-foreground") );
  m_setAsSecForeground->triggered().connect( boost::bind( &SpecMeasManager::loadSelected, this, SpecUtils::SpectrumType::SecondForeground, true ) );
  m_setButton->setMenu(setPopup);
  m_setButton->hide();

  // Note: the only difference between the m_combineToNewFileButton and
  //  m_subsetOfMeasToNewFileButton buttons are the icons, and a slight difference in text, but the
  //  functionality they trigger is the same either way.
  m_combineToNewFileButton = new WPushButton( WString::tr("smm-to-new-file"), m_newDiv );
  m_combineToNewFileButton->setToolTip( WString::tr("smm-tt-to-new-file") );
  m_combineToNewFileButton->addStyleClass("InvertInDark");
  m_combineToNewFileButton->setIcon( "InterSpec_resources/images/arrow_join.svg" );
  m_combineToNewFileButton->clicked().connect( boost::bind( &SpecMeasManager::newFileFromSelection, this ) );
  m_combineToNewFileButton->hide();
  
  m_subsetOfMeasToNewFileButton = new WPushButton( WString::tr("smm-as-new-file"), m_newDiv );
  m_subsetOfMeasToNewFileButton->setToolTip( WString::tr("smm-tt-as-new-file") );
  m_subsetOfMeasToNewFileButton->addStyleClass("InvertInDark");
  m_subsetOfMeasToNewFileButton->setIcon( "InterSpec_resources/images/partial.svg" );
  m_subsetOfMeasToNewFileButton->clicked().connect( boost::bind( &SpecMeasManager::newFileFromSelection, this ) );
  m_subsetOfMeasToNewFileButton->hide();
  
  
  m_sumSpectraButton = new WPushButton( WString::tr("smm-sum-spectra"), m_newDiv );
  m_sumSpectraButton->setToolTip( WString::tr("smm-tt-sum-spectra") );
  m_sumSpectraButton->addStyleClass("InvertInDark");
  m_sumSpectraButton->setIcon( "InterSpec_resources/images/sum_symbol.svg" );
  m_sumSpectraButton->clicked().connect( boost::bind( &SpecMeasManager::sumSelectedSpectra, this ) );
  m_sumSpectraButton->hide();
  
  
  m_saveButton = new WPushButton( WString::tr("smm-export-btn"), m_newDiv );
  m_saveButton->clicked().connect( this, &SpecMeasManager::startSaveSelected );
  m_saveButton->hide();
  
  
  m_deleteButton = new WPushButton( WString::tr("smm-unload-btn"), m_newDiv);
  m_deleteButton->setIcon( "InterSpec_resources/images/minus_min_white.png" );
  m_deleteButton->clicked().connect( boost::bind( &SpecMeasManager::removeSelected, this ) );
  m_deleteButton->hide();
  
  // ---- try new bar 2 ----
  WContainerWidget *m_newDiv2 = new WContainerWidget( );
  buttonAlignment->addWidget( m_newDiv2, 2, 0);
  m_newDiv2->setStyleClass( "LoadSpectrumUploadDiv" );
  
  //WText *text =
  new WText( WString::tr("smm-unassign-label"), m_newDiv2 );
  
  m_removeForeButton = new WPushButton( WString::tr("Foreground"), m_newDiv2 );
//      m_removeForeButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeForeButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::Foreground) );
  m_removeForeButton  ->setHidden( true, WAnimation() );
//  m_removeForeButton->setHiddenKeepsGeometry( true );
  
  m_removeBackButton = new WPushButton( WString::tr("Background"), m_newDiv2 );
//          m_removeBackButton->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeBackButton->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::Background ) );
  m_removeBackButton  ->setHidden( true, WAnimation() );
//  m_removeBackButton->setHiddenKeepsGeometry( true );
  
  m_removeFore2Button = new WPushButton( WString::tr("second-foreground"), m_newDiv2 );
//              m_removeFore2Button->setIcon( "InterSpec_resources/images/minus_min.png" );
  m_removeFore2Button->clicked().connect( boost::bind( &SpecMeasManager::unDisplay, this, SpecUtils::SpectrumType::SecondForeground ) );
  m_removeFore2Button->setHidden( true, WAnimation() );
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
void SpecMeasManager::saveToDatabase( std::shared_ptr<const SpecMeas> input ) const
{
  //Make sure we are in the event loop - for thread safety.
  if( !wApp )
    throw runtime_error( "SpecMeasManager::saveToDatabase() must be called from within the event loop." );
  
  std::shared_ptr<SpectraFileHeader> header = m_fileModel->fileHeader( input );
  
  if( !!header )
  {
    if( !header->shouldSaveToDb() )
      return;
    
    std::shared_ptr<SpecMeas> meas = header->parseFile();
    boost::function<void(void)> worker
                      = wApp->bind( boost::bind( &SpectraFileHeader::saveToDatabaseWorker, meas, header ) );
    
    // Instead of posting immediately, we'll
    //WServer::instance()->post( wApp->sessionId(), worker );
    
    boost::function<void ()> fallback = [](){
      cerr << "Failed to save file to database." << endl;
    };
    
    WServer::instance()->schedule( 50, wApp->sessionId(), worker, fallback );
  }//if( headermeas && (headermeas==meas) )
}//void saveToDatabase( std::shared_ptr<const SpecMeas> meas ) const


int SpecMeasManager::setDbEntry( Wt::Dbo::ptr<UserFileInDb> dbfile,
                                 std::shared_ptr<SpectraFileHeader> &header,
                                 std::shared_ptr<SpecMeas> &measurement,
                                 bool enforceUser )
{
  if( !dbfile )
    throw runtime_error( "SpecMeasManager::setDbEntry(...): invalid dbentry" );
  
  if( enforceUser && (dbfile->user != m_viewer->user()) )
    throw runtime_error( "SpecMeasManager::setDbEntry(...): invalid user" );
  
  int row = -1;
  header.reset();
  measurement.reset();
  header.reset( new SpectraFileHeader( false, m_viewer ) );
  measurement = header->resetFromDatabase( dbfile );
  addToTempSpectrumInfoCache( measurement );
  row = m_fileModel->addRow( header );
  return row;
}//int setDbEntry(...)


void SpecMeasManager::userCanceledResumeFromPreviousOpened( shared_ptr<SpectraFileHeader> header )
{
  std::shared_ptr<SpecMeas> meas = header->measurementIfInMemory();
  
  boost::function<void(void)> f;
  
  if( meas )
    f = boost::bind( &SpectraFileHeader::saveToDatabaseWorker, meas, header );
  else
    f = boost::bind( &SpectraFileHeader::saveToDatabaseFromTempFileWorker, header );
  
  WServer::instance()->post( wApp->sessionId(), f );
  //WServer::instance()->ioService().boost::asio::io_service::post( f );
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
      
      Wt::log("info") << "Setting file with UUID=" << header->m_uuid << ", and filename '"
      << header->m_displayName << "' to be connected to DB entry from "
      << "upload at "
      << f->uploadTime.toString(DATE_TIME_FORMAT_STR).toUTF8()
      << " until user determines if they want to resume from a modified"
      << " session.";
      
      header->setDbEntry( f );
    }catch( std::exception &e )
    {
      cerr << "Failed to set DB entry to SpectraFileHeader: " << e.what() << endl;
      
      assert( 0 );
    }//try / catch
  }//if( unModifiedFiles.size() )

  if( m_previousStatesDialog )
    handleCancelPreviousStatesDialog( m_previousStatesDialog );
  
  AuxWindow *window = new AuxWindow( WString::tr("smm-prev-states-window-title"),
                                    (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::DisableCollapse)
                                     | AuxWindowProperties::EnableResize
                                     | AuxWindowProperties::TabletNotFullScreen
                                     | AuxWindowProperties::SetCloseable) );
  window->rejectWhenEscapePressed();
  window->addStyleClass( "ShowPrevSpecFileUses" );
  WPushButton *cancel = window->addCloseButtonToFooter();
  cancel->clicked().connect( window, &AuxWindow::hide );
  
  //bool auto_save_states = false;
  SnapshotBrowser *snapshots = nullptr;
  AutosavedSpectrumBrowser *auto_saved = nullptr;
  
  try
  {
    //auto_save_states = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
    
    if( userStatesWithFile.size() )
    {
      // TODO: pass userStatesWithFile into SnapshotBrowser
      snapshots = new SnapshotBrowser( this, m_viewer, header, nullptr, nullptr );
      snapshots->finished().connect( boost::bind( &SpecMeasManager::handleClosePreviousStatesDialogAfterSelect,
                                      this, window) );
    }//if( userStatesWithFile.size() )
    
    
    if( !modifiedFiles.empty() )
    {
      auto_saved = new AutosavedSpectrumBrowser( modifiedFiles, type, m_fileModel, this, header );
      auto_saved->loadedASpectrum().connect( boost::bind( &SpecMeasManager::handleClosePreviousStatesDialogAfterSelect,
                                                         this, window) );
      
      if( unModifiedFiles.empty() )
      {
        window->finished().connect( boost::bind( &SpecMeasManager::userCanceledResumeFromPreviousOpened,
                                                this, header ) );
      }
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
    
    WString msg = WString::tr("smm-err-prev-state-unexpected").arg( e.what() );
    passMessage( msg, WarningWidget::WarningMsgLevel::WarningMsgHigh)
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
    
    tabbed->addTab( snapshots, WString::tr("smm-tab-saved-states"), WTabWidget::LoadPolicy::PreLoading );
    tabbed->addTab( auto_saved, WString::tr("smm-tab-auto-saved"), WTabWidget::LoadPolicy::PreLoading );
  }else
  {
    if( snapshots )
      layout->addWidget( snapshots, 0, 0 );
    else
      layout->addWidget( auto_saved, 0, 0 );
  }//if( snapshots && auto_saved ) / else
  
  WCheckBox *cb = new WCheckBox( WString::tr("smm-auto-check-prev-cb") );
  cb->addStyleClass( "PrefCb" );
  layout->addWidget( cb, 1, 0 );
  layout->setRowStretch( 0, 1 );
  
  UserPreferences::associateWidget( "CheckForPrevOnSpecLoad", cb, m_viewer );
  
  
  const int width = std::min( 500, static_cast<int>(0.95*m_viewer->renderedWidth()) );
  const int height = std::min( 475, static_cast<int>(0.95*m_viewer->renderedHeight()) );
  window->resize( WLength(width), WLength(height) );
  
  window->centerWindow();
  window->show();
  
  m_previousStatesDialog = window;
  window->finished().connect( boost::bind( &SpecMeasManager::handleCancelPreviousStatesDialog, this, window ) );
  
  wApp->triggerUpdate();
}//void showPreviousSpecFileUsesDialog(..)


void postErrorMessage( const WString msg, const WarningWidget::WarningMsgLevel level )
{
  passMessage( msg, level );
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
      = UserPreferences::preferenceValue<bool>( "CheckForPrevOnSpecLoad", m_viewer );

    if( !storeInDb )
      return;
    
    if( !header )
      throw runtime_error( "Invalid SpectraFileHeader passed in" );

    if( !m_viewer || !m_viewer->user() )
      throw runtime_error( "Invalid InterSpec or user pointer" );
    
    vector<Dbo::ptr<UserState>> userStatesWithFile;
    vector< Wt::Dbo::ptr<UserFileInDb> > modifiedFiles, unModifiedFiles;
    
    {//begin interaction with database
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_viewer->sql();
      const Wt::Dbo::ptr<InterSpecUser> &user = m_viewer->user();
      DataBaseUtils::DbTransaction transaction( *sql );
      
      typedef Dbo::collection< Dbo::ptr<UserFileInDb> > UserFileInDbColl;

//      UserFileInDbColl files = m_viewer->m_user->userFiles().find()
//                                          .where( "UUID = ? AND Filename = ? "
//                                                  "AND IsPartOfSaveState = 0" )
//                                          .bind( header->m_uuid )
//                                          .bind( header->m_displayName );
      UserFileInDbColl files = user->userFiles().find()
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
                    = SnapshotBrowser::get_user_states_collection( user, sql, header );
      
      for( auto iter = states_query.begin(); iter != states_query.end(); ++iter )
        userStatesWithFile.push_back( *iter );
      
      user.modify()->incrementSpectraFileOpened();
      
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
        WString msg = WString::tr("smm-cant-save").arg( e.message() );
        
        WServer::instance()->post( sessionID,
                  boost::bind( &postErrorMessage, msg, WarningWidget::WarningMsgHigh ) );
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
    WServer::instance()->post( sessionID,
                              boost::bind( &SpecMeasManager::showPreviousSpecFileUsesDialog,
                                          this, header, type, modifiedFiles, unModifiedFiles,
                                          userStatesWithFile ) );
    
  }catch( std::exception &e )
  {
    cerr << "Error checking if this file had been opened previously: " << e.what() << endl;
    WServer::instance()->post( sessionID,
                              boost::bind( &postErrorMessage,
              WString("Error checking if this file had been opened previously"),
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
        = std::make_shared<SpectraFileHeader>( false, m_viewer );
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
     = std::make_shared<SpectraFileHeader>( false, m_viewer );
  
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
  displayFile( row, measurement, type, true, true, 
              SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
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
  if( m_previousStatesDialog )
    handleCancelPreviousStatesDialog( m_previousStatesDialog );
  
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
//    passMessage( "Successfully uploaded file.", 0 );

    displayFile( row, measurement, type, true, true, 
              SpecMeasManager::VariantChecksToDo::DerivedDataAndMultiEnergyAndMultipleVirtualDets );
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
    
    if( !handleNonSpectrumFile( origName, name, type ) )
    {
      displayInvalidFileMsg(origName,e.what());
    }
    return false;
  }// try/catch
  
  return true;
}//void loadFromFileSystem( std::string filename )


int SpecMeasManager::dataUploaded( Wt::WFileUpload *upload, std::shared_ptr<SpecMeas> &measurement )
{
  if( m_previousStatesDialog )
    handleCancelPreviousStatesDialog( m_previousStatesDialog );
  
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

    //passMessage( "Successfully opened file.", 0 );

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

  for( queue_type::iterator iter = m_tempSpectrumInfoCache.begin();
      iter != m_tempSpectrumInfoCache.end(); ++iter )
  {
    serializeToTempFile( *iter );
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
     = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
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
      = UserPreferences::preferenceValue<bool>( "AutoSaveSpectraToDb", m_viewer );
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


void SpecMeasManager::fileTooLarge( const ::int64_t size_tried )
{
  const int max_size = static_cast<int>( WApplication::instance()->maximumRequestSize() );
  
  WString msg = WString::tr("smm-max-file-to-large")
                  .arg( static_cast<int>(max_size/1024) )
                  .arg( static_cast<int>(size_tried/1024) );
  
  passMessage( msg, WarningWidget::WarningMsgHigh );
} // void SpecMeasManager::fileTooLarge()


void SpecMeasManager::uploadSpectrum() {
  new UploadBrowser(this/*, m_viewer*/);
}


void SpecMeasManager::startSaveSelected()
{
  if (!m_spectrumManagerWindow || m_spectrumManagerWindow->isHidden())
    return;
  
  // We shouldnt see these error messages (unless there is a logic error), so we
  //  wont internationalize.
  
  vector<shared_ptr<SpectraFileHeader>> files = getSelectedFiles();
  if( files.empty() )
  {
    passMessage( "No files selected for export.", WarningWidget::WarningMsgHigh );
    return;
  }
  
  if( files.size() > 1 )
  {
    passMessage( "More than one files selected for export; not allowed.", WarningWidget::WarningMsgHigh );
    return;
  }
  
  shared_ptr<SpectraFileHeader> header = files.front();
  
  
  shared_ptr<SpecMeas> spec = header ? header->parseFile() : nullptr;
  if( !spec )
  {
    passMessage( "Programming logic error: couldn't locate spectrum file to export.", WarningWidget::WarningMsgHigh );
    return;
  }
  
  const set<int> sample_numbers = selectedSampleNumbers();
  if( sample_numbers.empty() )
  {
    passMessage( "Error determining sample numbers to export.", WarningWidget::WarningMsgHigh );
    return;
  }
  
  // The undo/redo wont be perfect, because if the user does undo, and then redo, the dialog
  //  wont have the specific spectrum file we want.
  //  We could do things properly, but then we have to trac the created ExportSpecFileWindow
  //  window in this class; for a first go this added complexity isnt worth it for the edge-case.
  ExportSpecFileWindow *window = m_viewer->createExportSpectrumFileDialog();
  if( !window )
  {
    passMessage( "Programming logic error creating export dialog", WarningWidget::WarningMsgHigh );
    return;
  }
  
  window->setSpecificSpectrum( spec, sample_numbers, spec->detector_names(), m_viewer );
}//void startSaveSelected()


#if( USE_DB_TO_STORE_SPECTRA )
void SpecMeasManager::browsePrevSpectraAndStatesDb()
{
  // TODO: Make this be the same implementation as SpecMeasManager::showPreviousSpecFileUsesDialog; but to do that, need to make AutosavedSpectrumBrowser be a MVC widget so we dont put like a million elements into the DOM
  DbFileBrowser *browser = new DbFileBrowser( this, m_viewer, nullptr );
  
  UndoRedoManager *undoRedo = m_viewer->undoRedoManager();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto closer = wApp->bind( boost::bind( &AuxWindow::hide, browser ) );
    undoRedo->addUndoRedoStep( closer, nullptr, "Show browse previous spectra." );
  }//if( add undo )
}//void browsePrevSpectraAndStatesDb()


#endif //#if( USE_DB_TO_STORE_SPECTRA )


