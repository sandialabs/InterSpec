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

#include <regex>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <exception>

#include <Wt/Utils>
#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WServer>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WTabWidget>
#include <Wt/WPushButton>
#include <Wt/Http/Request>
#include <Wt/WApplication>
#include <Wt/Http/Response>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WMessageResourceBundle>

//Following include only needed for lines:
//  answer->setShieldingSourceModel(...);
//  answer->setRelActManualGuiState(...);
#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml.hpp"
#include "external_libs/SpecUtils/3rdparty/rapidxml/rapidxml_print.hpp"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/UriSpectrum.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ExportSpecFile.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

#if( USE_QR_CODES )
#include "InterSpec/QrCode.h"
#include "InterSpec/QRSpectrum.h"
#endif

using namespace std;
using namespace Wt;

#if( ANDROID )
// Defined in target/android/android.cpp
extern void android_download_workaround( Wt::WResource *resource, std::string description );
#endif

namespace
{
  void right_select_item( WMenu *menu, WMenuItem *item )
  {
    menu->select( item );
    item->triggered().emit( item ); //
  }

const std::vector<std::pair<SpecUtils::SaveSpectrumAsType,std::string>> sm_file_type_to_uri_name = {
  { SpecUtils::SaveSpectrumAsType::N42_2012, "N42-2012"},
  { SpecUtils::SaveSpectrumAsType::N42_2006, "N42-2006"},
  { SpecUtils::SaveSpectrumAsType::Chn, "CHN" },
  { SpecUtils::SaveSpectrumAsType::SpeIaea, "IAEA-SPE" },
  { SpecUtils::SaveSpectrumAsType::Csv, "CSV" },
  { SpecUtils::SaveSpectrumAsType::Txt, "TXT" },
  { SpecUtils::SaveSpectrumAsType::Pcf, "PCF" },
  { SpecUtils::SaveSpectrumAsType::Cnf, "CNF" },
  { SpecUtils::SaveSpectrumAsType::SpcBinaryInt, "SPC-int" },
  { SpecUtils::SaveSpectrumAsType::SpcBinaryFloat, "SPC-float" },
  { SpecUtils::SaveSpectrumAsType::SpcAscii, "SPC-ascii" },
  { SpecUtils::SaveSpectrumAsType::Tka, "TKA" },
  { SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0, "GR-130" },
  { SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2, "GR-135" },
#if( SpecUtils_ENABLE_D3_CHART )
  { SpecUtils::SaveSpectrumAsType::HtmlD3, "HTML" },
#endif
#if( USE_QR_CODES )
  { SpecUtils::SaveSpectrumAsType::Uri, "QR-code" },
#endif
};//sm_file_type_to_uri_name
}//namespace

#if( USE_QR_CODES )
namespace
{
  const int ns_min_screen_width_for_wide_layout = 750;
  
std::string clean_uuid( string uuid )
{
  SpecUtils::trim(uuid);
  for( size_t i = 0; i < uuid.size(); ++i )
  {
    if( !std::isalpha( (int)uuid[i] ) && !std::isalnum( (int)uuid[i] ) )
      uuid[i] = ' ';
  }
  
  SpecUtils::erase_any_character( uuid, " $&+,:;=?@'\"<>#%{}|\\^~[]`/" );
  
  for( size_t i = 0; i < uuid.size(); ++i )
  {
    if( !std::isalpha( (int)uuid[i] ) && !std::isalnum( (int)uuid[i] ) )
      uuid[i] = '-';
  }
  
  return uuid;
}//std::string clean_uuid( string uuid )


  void displayQrDialog( const vector<QRSpectrum::QrCodeEncodedSpec> urls,
                       const size_t index,
                       boost::function<void()> successfullyDone,
                       const bool as_emailto,
                       boost::function<void()> toogleEmailVsUri )
  {
    string seqnum = as_emailto ? "QR-code to create email" : "QR-code of spectrum";
    if( urls.size() > 1 )
      seqnum += ", " + to_string(index + 1) + " of " + to_string( urls.size() );
    
    string desc = seqnum;
    string title = seqnum;
    if( urls.size() == 1 )
    {
      title = "";
      
      // If a single URL, and on a phone, save screen space by skipping the description
      InterSpec *viewer = InterSpec::instance();
      if( viewer && viewer->isPhone() )
        desc = "";
    }//if( urls.size() == 1 )
    
    
    SimpleDialog *dialog = QrCode::displayTxtAsQrCode( urls[index].m_url, WString::fromUTF8(title),
                                                      WString::fromUTF8(desc) );
    
    if( (index + 1) < urls.size() )
    {
      WPushButton *btn = dialog->addButton( "Next QR code" );
      btn->clicked().connect( std::bind([=](){
        displayQrDialog( urls, index + 1, successfullyDone, as_emailto, toogleEmailVsUri );
      }) );
    }else if( successfullyDone )
    {
      // We need to find the "Close" button here
      WPushButton *close = nullptr;
      
      WContainerWidget *foot = dialog->footer();
      if( foot )
      {
        for( auto w : foot->children() )
        {
          auto btn = dynamic_cast<WPushButton *>(w);
          close = btn;
          if( close )
            break;
        }//for( auto w : foot->children() )
      }//if( foot )
      assert( close );
      if( close )
        close->clicked().connect( std::bind(successfullyDone) );
    }//if( (index + 1) < urls.size() )
    
    if( toogleEmailVsUri )
    {
      const char *btn_txt = as_emailto ? "URL Link" : "Email Link";
      WPushButton *btn = dialog->addButton( btn_txt );
      dialog->footer()->removeWidget( btn );
      dialog->footer()->insertWidget( 0, btn );
      btn->clicked().connect( std::bind(toogleEmailVsUri) );
    }//if( toogleEmailVsUri )
  }//void displayQr( vector<QRSpectrum::QrCodeEncodedSpec> urls )


void displayQrCode( const vector<SpecUtils::UrlSpectrum> urlspec,
                   shared_ptr<const SpecMeas> spec,
                   boost::function<void()> successfullyDone,
                   const bool as_emailto )
{
  try
  {
    uint8_t encode_options = 0;
    if( as_emailto )
      encode_options = SpecUtils::EncodeOptions::UseUrlSafeBase64
                        | SpecUtils::EncodeOptions::AsMailToUri;
    
    vector<QRSpectrum::QrCodeEncodedSpec> urls;
    QRSpectrum::QrErrorCorrection ecc = QRSpectrum::QrErrorCorrection::High;
    
    auto do_encode = [&](){
      // We'll try to only use a single QR-code, but also use highest error-correction we can.
      //  I'm sure this trade-off can be better handled in some way.
      if( spec->num_gamma_channels() < 2075 )
      {
        try
        {
          urls = QRSpectrum::qr_code_encode_spectra( urlspec, ecc, encode_options );
        }catch( std::exception & )
        {
        }//try / catch
      }//if( spec->num_gamma_channels() < 2075 )
      
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
    };//do_encode lamda
    
    do_encode();
    
    if( urls.empty() )
    {
      string msg = (urlspec.size() == 1) ? "<p>Likely due to requiring more than 9 QR codes<p>"
                                         : "<p>Likely due to not being able to fit multiple spectra"
                                           " into a single QR code.</p>";
      
      auto dialog = new SimpleDialog( "Error",
                                     "<p>Spectrum could not be encoded to a QR code.</p>" + msg );
      dialog->addButton( "Ok" );
      return;
    }//if( urls.empty() )
    
    
    boost::function<void()> toogleEmailVsUri = [=](){
      displayQrCode( urlspec, spec, successfullyDone, !as_emailto );
    };
    
    displayQrDialog( urls, 0, successfullyDone, as_emailto, toogleEmailVsUri );
  }catch( std::exception &e )
  {
    auto dialog = new SimpleDialog( "Error", "Failed to display spectrum as a QR code: " + string(e.what()) );
    dialog->addButton( "Ok" );
  }//try catch
}//void displayQrCode( const vector<SpecUtils::UrlSpectrum> urlspec, const bool as_emailto )
}//namespace
#endif //USE_QR_CODES


namespace ExportSpecFileTool_imp
{
  class DownloadSpectrumResource : public Wt::WResource
  {
    //A simple class to stream a spectrum upon demand for which ever items are
    //  highlighted in the tree view.
  public:
    DownloadSpectrumResource( ExportSpecFileTool *tool );
    virtual ~DownloadSpectrumResource();

    Wt::Signal<> &downloadFinished();
    
  private:
    ExportSpecFileTool *m_tool;
    Wt::WApplication *m_app;
    Wt::Signal<> m_download_finished;

    virtual void handleRequest( const Wt::Http::Request& request,
                                 Wt::Http::Response& response);

  public:
    
    static void write_file( std::ostream &output,
                            const SpecUtils::SaveSpectrumAsType type,
                            const std::shared_ptr<const SpecMeas> Measurement,
                            const std::set<int> &samplenums,
                            const std::vector<std::string> &detectornums,
                            const bool strip_interspec_info );
  };// class DownloadSpectrumResource


  DownloadSpectrumResource::DownloadSpectrumResource( ExportSpecFileTool *tool )
    : WResource( tool ),
      m_tool( tool ),
      m_app( WApplication::instance() ),
      m_download_finished( this )
  {
    assert( m_app );
  }


  DownloadSpectrumResource::~DownloadSpectrumResource()
  {
    beingDeleted();
  }

  Wt::Signal<> &DownloadSpectrumResource::downloadFinished()
  {
    return m_download_finished;
  }

  void DownloadSpectrumResource::handleRequest( const Wt::Http::Request& request,
                                                Wt::Http::Response& response)
  {
    WApplication::UpdateLock lock( m_app );
    
    assert( m_tool );
    if( !m_tool )
      return;
    
    try
    {
      if( !lock )
      {
        log("error") << "Failed to WApplication::UpdateLock in DownloadSpectrumResource.";
#if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, "Failed to WApplication::UpdateLock in DownloadSpectrumResource." );
#endif
        assert( 0 );
        
        throw runtime_error( "Error grabbing application lock to form DownloadSpectrumResource resource; please report to InterSpec@sandia.gov." );
      }//if( !lock )
    
      const SpecUtils::SaveSpectrumAsType save_type = m_tool->currentSaveType();
      const shared_ptr<const SpecMeas> output = m_tool->generateFileToSave();
      const bool strip_interspec_info = m_tool->removeInterSpecInfo();
      
      assert( output );
      if( !output )
        throw runtime_error( "could not generate file." );
      
      switch( save_type )
      {
        case SpecUtils::SaveSpectrumAsType::Txt:
          response.setMimeType( "application/octet-stream" );
    //      response.setMimeType( "text/plain" );
        break;

        case SpecUtils::SaveSpectrumAsType::Csv:
          response.setMimeType( "application/octet-stream" );
    //      response.setMimeType( "text/csv" );
        break;

        case SpecUtils::SaveSpectrumAsType::Pcf:
          response.setMimeType( "application/octet-stream" );
        break;

        case SpecUtils::SaveSpectrumAsType::N42_2006:
          response.setMimeType( "application/octet-stream" );
    //      response.setMimeType( "application/xml" );
        break;

        case SpecUtils::SaveSpectrumAsType::N42_2012:
          response.setMimeType( "application/octet-stream" );
          //      response.setMimeType( "application/xml" );
        break;
          
        case SpecUtils::SaveSpectrumAsType::Chn:
        case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
        case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
        case SpecUtils::SaveSpectrumAsType::SpcAscii:
        case SpecUtils::SaveSpectrumAsType::SpeIaea:
        case SpecUtils::SaveSpectrumAsType::Cnf:
        case SpecUtils::SaveSpectrumAsType::Tka:
          response.setMimeType( "application/octet-stream" );
        break;
          
        case SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0:
        case SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2:
          response.setMimeType( "application/octet-stream" );
        break;
          
    #if( SpecUtils_ENABLE_D3_CHART )
        case SpecUtils::SaveSpectrumAsType::HtmlD3:
          response.setMimeType( "text/html" );
        break;
    #endif //#if( USE_D3_EXPORTING )
          
#if( USE_QR_CODES )
        case SpecUtils::SaveSpectrumAsType::Uri:
        break;
#endif
          
        case SpecUtils::SaveSpectrumAsType::NumTypes:
        break;
      }//switch( type )
      
      
      write_file( response.out(), save_type, output,
                 output->sample_numbers(), output->detector_names(), strip_interspec_info );
      
      m_download_finished.emit();
    }catch( std::exception &e )
    {
      const string msg = "Failed to prepare spectrum file for export: " + string(e.what());
      
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, msg.c_str() );
#endif
      
      passMessage( msg, WarningWidget::WarningMsgHigh );
      
      response.setStatus(500);
    }//try / catch
  }//void DownloadSpectrumResource::handleRequest(...)


  void DownloadSpectrumResource::write_file( std::ostream &output,
                         const SpecUtils::SaveSpectrumAsType type,
                         const std::shared_ptr<const SpecMeas> measurement,
                         const std::set<int> &samplenums,
                         const std::vector<std::string> &detectornames,
                         const bool strip_interspec_info )
  {
    //Convert detector names to detector numbers.
    set<int> detectornums;
    const auto &names = measurement->detector_names();
    const auto &numbers = measurement->detector_numbers();
    assert( names.size() == numbers.size() );
    
    for( const string &n : detectornames )
    {
      auto pos = std::find( begin(names), end(names), n );
      if( pos != end(names) )
        detectornums.insert( numbers[pos - begin(names)] );
    }
    
    
    switch( type )
    {
      case SpecUtils::SaveSpectrumAsType::Txt:
        measurement->write_txt( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::Csv:
        measurement->write_csv( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::Pcf:
        measurement->write_pcf( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::N42_2006:
        if( strip_interspec_info )
          measurement->SpecUtils::SpecFile::write_2006_N42( output );
        else
          measurement->write_2006_N42( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::N42_2012:
      {
        shared_ptr<rapidxml::xml_document<char>> doc = strip_interspec_info
             ? measurement->SpecUtils::SpecFile::create_2012_N42_xml()
             : measurement->create_2012_N42_xml();
        
        if( doc )
        {
          // We currently put in a node like:
          //  <RadInstrumentInformationExtension>
          //    <InterSpec:DetectorType>IdentiFINDER-NG</InterSpec:DetectorType>
          //  </RadInstrumentInformationExtension>
          // Lets get rid of it.
          rapidxml::xml_node<char> *data_node = doc->first_node("RadInstrumentData");
          rapidxml::xml_node<char> *inst_info_node = data_node
                                    ? data_node->first_node( "RadInstrumentInformation" )
                                    : nullptr;
          rapidxml::xml_node<char> *ext_node = inst_info_node
                                    ? inst_info_node->first_node( "RadInstrumentInformationExtension" )
                                    : nullptr;
          if( ext_node )
            inst_info_node->remove_node( ext_node );
            
          rapidxml::print( output, *doc );
        }//if( doc )
        
        break;
      }//case SpecUtils::SaveSpectrumAsType::N42_2012:
        
      case SpecUtils::SaveSpectrumAsType::Chn:
        measurement->write_integer_chn( output, samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
        measurement->write_binary_spc( output,
                                      SpecUtils::SpecFile::IntegerSpcType,
                                      samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
        measurement->write_binary_spc( output,
                                      SpecUtils::SpecFile::FloatSpcType,
                                      samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::SpcAscii:
        measurement->write_ascii_spc( output, samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0:
        measurement->write_binary_exploranium_gr130v0( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2:
        measurement->write_binary_exploranium_gr135v2( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::SpeIaea:
        if( strip_interspec_info )
          measurement->SpecUtils::SpecFile::write_iaea_spe( output, samplenums, detectornums );
        else
          measurement->write_iaea_spe( output, samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::Cnf:
        measurement->write_cnf( output, samplenums, detectornums );
        break;
        
      case SpecUtils::SaveSpectrumAsType::Tka:
        measurement->write_tka( output, samplenums, detectornums );
        break;
        
        
  #if( SpecUtils_ENABLE_D3_CHART )
      case SpecUtils::SaveSpectrumAsType::HtmlD3:
      {
        InterSpec *viewer = InterSpec::instance();
        //For purposes of development, lets cheat and export everything as it is now
        if( viewer )
        {
          const SpecUtils::SpectrumType types[] = { SpecUtils::SpectrumType::Foreground, SpecUtils::SpectrumType::SecondForeground, SpecUtils::SpectrumType::Background };
          
          std::vector< std::pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > measurements;
          
          for( SpecUtils::SpectrumType type : types )
          {
            std::shared_ptr<const SpecUtils::Measurement> histogram = viewer->displayedHistogram( type );
            std::shared_ptr<const SpecMeas> data = viewer->measurment( type );
            
            if( !histogram || !data )
              continue;
            
            //string title = Wt::WWebWidget::escapeText(data->title()).toUTF8();
            //if( title != data->title() )
            //data->set_title( title );  //JIC, proper escaping not implemented in SpecUtils yet.
            
            D3SpectrumExport::D3SpectrumOptions options;
            switch( type )
            {
              case SpecUtils::SpectrumType::Foreground: options.line_color = "black"; break;
              case SpecUtils::SpectrumType::SecondForeground: options.line_color = "steelblue"; break;
              case SpecUtils::SpectrumType::Background: options.line_color = "green"; break;
            }
            options.display_scale_factor = viewer->displayScaleFactor( type );
            options.spectrum_type = type;
            const std::set<int> samples = viewer->displayedSamples( type );
            
            std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > peaks = data->peaks(samples);
            if( peaks )
            {
              vector< std::shared_ptr<const PeakDef> > inpeaks( peaks->begin(), peaks->end() );
              std::shared_ptr<const SpecUtils::Measurement> foreground = viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
              options.peaks_json = PeakDef::peak_json( inpeaks, foreground );
            }
            
            measurements.push_back( pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions>(histogram.get(),options) );
          }//for( SpecUtils::SpectrumType type : types )
          
          write_d3_html( output, measurements, viewer->getD3SpectrumOptions() );
        }//if( viewer )
      }
        break;
  #endif //#if( USE_D3_EXPORTING )
        
#if( USE_QR_CODES )
      case SpecUtils::SaveSpectrumAsType::Uri:
        assert( 0 );
        measurement->write_uri( output, 1, 0x0 );
      break;
#endif
      case SpecUtils::SaveSpectrumAsType::NumTypes:
        break;
    }//switch( type )
  }//DownloadSpectrumResource::write_file(...)

}//namespace ExportSpecFileTool_imp


ExportSpecFileTool::ExportSpecFileTool( InterSpec *viewer, Wt::WContainerWidget *parent )
 : Wt::WContainerWidget( parent ),
  m_interspec ( viewer ),
  m_is_specific_file( false ),
  m_specific_spectrum( nullptr ),
  m_specific_samples{},
  m_specific_detectors{},
  m_current_file(),
  m_done( this ),
  m_fileSelect( nullptr ),
  m_forePlusBack( nullptr ),
  m_fileInfo( nullptr ),
  m_formatMenu( nullptr ),
  m_samplesHolder( nullptr ),
  m_dispForeSamples( nullptr ),
  m_dispBackSamples( nullptr ),
  m_dispSecondSamples( nullptr ),
  m_allSamples( nullptr ),
  m_customSamples( nullptr ),
  m_customSamplesEdit( nullptr ),
  m_filterDetector( nullptr ),
  m_detectorFilterCbs( nullptr ),
  m_optionsHolder( nullptr ),
  m_sumAllToSingleRecord( nullptr ),
  m_sumForeToSingleRecord( nullptr ),
  m_sumBackToSingleRecord( nullptr ),
  m_sumSecoToSingleRecord( nullptr ),
  m_backSubFore( nullptr ),
  m_sumDetsPerSample( nullptr ),
  m_sumSamplesPerDets( nullptr ),
  m_excludeInterSpecInfo( nullptr ),
  m_excludeGpsInfo( nullptr ),
  m_msg( nullptr ),
  m_sampleSelectNotAppTxt( nullptr ),
  m_optionsNotAppTxt( nullptr ),
  m_export_btn( nullptr ),
#if( USE_QR_CODES )
  m_show_qr_btn( nullptr ),
#endif
  m_resource( nullptr ),
  m_last_state_uri{}
{
  init();
}//ExportSpecFileTool()


ExportSpecFileTool::ExportSpecFileTool( const std::shared_ptr<const SpecMeas> &spectrum,
                   const std::set<int> &samples,
                   const std::vector<std::string> &detectors,
                   InterSpec *viewer,
                   Wt::WContainerWidget *parent )
: Wt::WContainerWidget( parent ),
  m_interspec ( viewer ),
  m_is_specific_file( false ),
  m_specific_spectrum( spectrum ),
  m_specific_samples{ samples },
  m_specific_detectors{ detectors },
  m_current_file( spectrum ),
  m_done( this ),
  m_fileSelect( nullptr ),
  m_forePlusBack( nullptr ),
  m_fileInfo( nullptr ),
  m_formatMenu( nullptr ),
  m_samplesHolder( nullptr ),
  m_dispForeSamples( nullptr ),
  m_dispBackSamples( nullptr ),
  m_dispSecondSamples( nullptr ),
  m_allSamples( nullptr ),
  m_customSamples( nullptr ),
  m_customSamplesEdit( nullptr ),
  m_filterDetector( nullptr ),
  m_detectorFilterCbs( nullptr ),
  m_optionsHolder( nullptr ),
  m_sumAllToSingleRecord( nullptr ),
  m_sumForeToSingleRecord( nullptr ),
  m_sumBackToSingleRecord( nullptr ),
  m_sumSecoToSingleRecord( nullptr ),
  m_backSubFore( nullptr ),
  m_sumDetsPerSample( nullptr ),
  m_sumSamplesPerDets( nullptr ),
  m_excludeInterSpecInfo( nullptr ),
  m_excludeGpsInfo( nullptr ),
  m_msg( nullptr ),
  m_sampleSelectNotAppTxt( nullptr ),
  m_optionsNotAppTxt( nullptr ),
  m_export_btn( nullptr ),
#if( USE_QR_CODES )
  m_show_qr_btn( nullptr ),
#endif
  m_resource( nullptr ),
  m_last_state_uri{}
{
  init();
}//
  

void ExportSpecFileTool::init()
{
  if( !m_interspec )
    throw logic_error( "ExportSpecFileTool: Invalid InterSpec pointer." );
  
  addStyleClass( "ExportSpecFileTool" );
  
  // Store Act/Shield fit and Rel Eff fit into file - if these tools are open, so the
  //   tools current states will be available in the foreground N42 files.
  m_interspec->saveShieldingSourceModelToForegroundSpecMeas();
#if( USE_REL_ACT_TOOL )
  m_interspec->saveRelActManualStateToForegroundSpecMeas();
#endif
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_interspec );
  const bool isMobile = m_interspec && m_interspec->isMobile();

  if( isMobile )
    addStyleClass( "ExportSpecFileToolMobile" );
  
  const int screenWidth = m_interspec->renderedWidth();
  const bool isPhone = (screenWidth > 100) ? (screenWidth < ns_min_screen_width_for_wide_layout) : m_interspec->isPhone();
  
  WTabWidget *mobileTabs = nullptr;
  WContainerWidget *body = nullptr;

  if( isPhone )
  {
    addStyleClass( "ExportSpecFileToolPhone" );

    mobileTabs = new WTabWidget( this );
    mobileTabs->addStyleClass( "ExportSpecFileTabs" );
  }else
  {
    body = new WContainerWidget( this );
    body->addStyleClass( "ExportSpecFileBody" );
  }
  
  
  WContainerWidget *fileSelectDiv = new WContainerWidget( body );
  fileSelectDiv->addStyleClass( "ExportSpecSelect" );
  if( isPhone )
    mobileTabs->addTab( fileSelectDiv, "File", Wt::WTabWidget::LoadPolicy::PreLoading );
  
  if( !m_specific_spectrum )
  {
    WText *title = new WText( "File to Export", fileSelectDiv );
    title->addStyleClass( "ExportColTitle" );
    
    m_fileSelect = new WComboBox( fileSelectDiv );
    m_fileSelect->setNoSelectionEnabled( true );
    
    SpecMeasManager * const measManager = m_interspec->fileManager();
    SpectraFileModel * const fileModel = measManager ? measManager->model() : nullptr;
    if( !measManager || !fileModel )
      throw runtime_error( "null SpecMeasManager or SpectraFileModel" );
    
    m_fileSelect->setModel( fileModel );
    m_fileSelect->setModelColumn( SpectraFileModel::DisplayFields::kDisplayName );
    m_fileSelect->activated().connect( this, &ExportSpecFileTool::handleFileSelectionChanged );
    
    
    const auto fore = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    const auto back = m_interspec->measurment(SpecUtils::SpectrumType::Background);
    const auto secondary = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
    
    const bool unique_back = (back && (back != fore));
    const bool unique_seco = (secondary && ((secondary != fore) || (back && (back != secondary))));
                        
    
    if( fore && (unique_back || unique_seco) )
    {
      string title = "For." + string(back ? " + Back." : "") + string(secondary ? " + Secon." : "");
      
      m_forePlusBack = new WCheckBox( title, fileSelectDiv );
      m_forePlusBack->addStyleClass( "ExportForPlusBack" );
      m_forePlusBack->checked().connect( this, &ExportSpecFileTool::handleForePlusBackChanged );
      m_forePlusBack->unChecked().connect( this, &ExportSpecFileTool::handleForePlusBackChanged );
      m_forePlusBack->setToolTip( "Export the foreground and background and/or secondary"
                                 " spectra together into the same file." );
    }
  }//if( !m_specific_spectrum )
  
  m_fileInfo = new WText( fileSelectDiv );
  m_fileInfo->addStyleClass( "ExportSpecInfo" );
  
  // Spectrum format
  WContainerWidget *menuHolder = new WContainerWidget( body );
  menuHolder->addStyleClass( "ExportSpecFormat" );
  if( isPhone )
    mobileTabs->addTab( menuHolder, "Format", Wt::WTabWidget::LoadPolicy::PreLoading );
  
  WText *title = new WText( "File Format", menuHolder );
  title->addStyleClass( "ExportColTitle" );
  
  m_formatMenu = new WMenu( menuHolder );
  m_formatMenu->itemSelected().connect( this, &ExportSpecFileTool::handleFormatChange );
  m_formatMenu->addStyleClass( "SideMenu VerticalNavMenu LightNavMenu ExportSpecFormatMenu" );
  
  Wt::WMessageResourceBundle descrip_bundle;
  if( !isMobile )
  {
    const string docroot = wApp->docRoot();
    const string bundle_file = SpecUtils::append_path(docroot, "InterSpec_resources/static_text/spectrum_file_format_descriptions" );
    descrip_bundle.use(bundle_file,true);
  }//if( !isMobile )
  
  
  auto addFormatItem = [this, &descrip_bundle, isMobile]( const char *label, SpecUtils::SaveSpectrumAsType type ){
    WMenuItem *item = m_formatMenu->addItem( label );
    item->clicked().connect( boost::bind(&right_select_item, m_formatMenu, item) );
    item->setData( reinterpret_cast<void *>(type) );
    
    if( !isMobile )
    {
      string description;
      if( descrip_bundle.resolveKey(label, description) )
      {
        SpecUtils::trim( description );
       
        description = Wt::Utils::htmlEncode( description, Wt::Utils::HtmlEncodingFlag::EncodeNewLines );
        
        WImage *img = new WImage( item );
        img->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
        img->setStyleClass("Wt-icon");
        img->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
        img->setFloatSide( Wt::Side::Right );
        
        HelpSystem::attachToolTipOn( img, description, true,
                                    HelpSystem::ToolTipPosition::Right,
                                    HelpSystem::ToolTipPrefOverride::InstantAlways );
      }//if( we have the description of the file )
    }//if( !isMobile )
  };//addFormatItem lambda
  
  
  addFormatItem( "N42-2012", SpecUtils::SaveSpectrumAsType::N42_2012 );
  addFormatItem( "N42-2006", SpecUtils::SaveSpectrumAsType::N42_2006 );
  addFormatItem( "CHN", SpecUtils::SaveSpectrumAsType::Chn );
  addFormatItem( "IAEA SPE", SpecUtils::SaveSpectrumAsType::SpeIaea );
  addFormatItem( "CSV", SpecUtils::SaveSpectrumAsType::Csv );
  addFormatItem( "TXT", SpecUtils::SaveSpectrumAsType::Txt );
  addFormatItem( "PCF", SpecUtils::SaveSpectrumAsType::Pcf );
  addFormatItem( "CNF", SpecUtils::SaveSpectrumAsType::Cnf );
  addFormatItem( "SPC (int)", SpecUtils::SaveSpectrumAsType::SpcBinaryInt );
  addFormatItem( "SPC (float)", SpecUtils::SaveSpectrumAsType::SpcBinaryFloat );
  addFormatItem( "SPC (ascii)", SpecUtils::SaveSpectrumAsType::SpcAscii );
  addFormatItem( "TKA", SpecUtils::SaveSpectrumAsType::Tka );
  addFormatItem( "GR-130", SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0 );
  addFormatItem( "GR-135", SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2 );
#if( SpecUtils_ENABLE_D3_CHART )
  addFormatItem( "HTML", SpecUtils::SaveSpectrumAsType::HtmlD3 );
#endif
#if( USE_QR_CODES )
  addFormatItem( "QR-code/URL", SpecUtils::SaveSpectrumAsType::Uri );
#endif
  
  // Meas/samples to include
  m_samplesHolder = new WContainerWidget( body );
  m_samplesHolder->addStyleClass( "ExportSpecSamples" );

  if( isPhone )
    mobileTabs->addTab( m_samplesHolder, "Options", Wt::WTabWidget::LoadPolicy::PreLoading );

  title = new WText( "Samples to Include", m_samplesHolder );
  title->addStyleClass( "ExportColTitle" );
  
  m_dispForeSamples   = new WCheckBox( "Disp. Foreground", m_samplesHolder );
  m_dispForeSamples->checked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::Foreground ) );
  m_dispForeSamples->unChecked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::Foreground ) );
  
  m_dispBackSamples   = new WCheckBox( "Disp. Background", m_samplesHolder );
  m_dispBackSamples->checked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::Background ) );
  m_dispBackSamples->unChecked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::Background ) );
  
  m_dispSecondSamples = new WCheckBox( "Disp. Secondary", m_samplesHolder );
  m_dispSecondSamples->checked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::SecondForeground ) );
  m_dispSecondSamples->unChecked().connect( boost::bind( &ExportSpecFileTool::handleDisplaySampleChanged, this, SpecUtils::SpectrumType::SecondForeground ) );
  
  m_allSamples = new WCheckBox( "All Samples", m_samplesHolder );
  m_allSamples->checked().connect( this, &ExportSpecFileTool::handleAllSampleChanged );
  m_allSamples->unChecked().connect( this, &ExportSpecFileTool::handleAllSampleChanged );
  
  m_customSamples = new WCheckBox( "Custom Samples", m_samplesHolder );
  m_customSamples->checked().connect( this, &ExportSpecFileTool::handleCustomSampleChanged );
  m_customSamples->unChecked().connect( this, &ExportSpecFileTool::handleCustomSampleChanged );
  
  
  const char *val_regex = "(\\-?\\d+\\s*((\\-|to|through)\\s*\\d+)?,*\\s*)+";
  WRegExpValidator *validator = new WRegExpValidator( val_regex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  
  m_customSamplesEdit = new WLineEdit( m_samplesHolder );
  m_customSamplesEdit->setAttributeValue( "ondragstart", "return false" );
  m_customSamplesEdit->addStyleClass( "SampleNumInput ExportSpecCustomSamples" );
  m_customSamplesEdit->setValidator( validator );
  m_customSamplesEdit->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_customSamplesEdit->setAttributeValue( "autocorrect", "off" );
  m_customSamplesEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_customSamplesEdit->setTextSize( 10 );
  
  m_customSamplesEdit->enterPressed().connect( this, &ExportSpecFileTool::handleCustomSampleTxtChanged );
  m_customSamplesEdit->blurred().connect( this, &ExportSpecFileTool::handleCustomSampleTxtChanged );
  
  
  const char *tooltip = "Enter the sample number you'de like displayed here"
  ". You may enter a range of sample numbers similar"
  " to '33-39' or '33 to 39', or CSV sample numbers"
  " like '34,39,84'";
  HelpSystem::attachToolTipOn( m_customSamplesEdit, tooltip, showToolTips, HelpSystem::ToolTipPosition::Right );
  
  m_customSamplesEdit->hide();

  m_filterDetector = new WCheckBox( "Filter Detectors", m_samplesHolder );
  m_filterDetector->checked().connect( this, &ExportSpecFileTool::handleFilterDetectorCbChanged );
  m_filterDetector->unChecked().connect( this, &ExportSpecFileTool::handleFilterDetectorCbChanged );
  
  m_detectorFilterCbs = new WContainerWidget( m_samplesHolder );
  m_detectorFilterCbs->addStyleClass( "ExportDetsToFilter" );
  m_filterDetector->hide();
  m_detectorFilterCbs->hide();
  
  // Options
  m_optionsHolder = new WContainerWidget( m_samplesHolder );
  m_optionsHolder->addStyleClass( "ExportSpecOptions" );
  title = new WText( "Options", m_optionsHolder );
  title->addStyleClass( "ExportColTitle" );
  
  
  m_sumAllToSingleRecord = new WCheckBox( "Sum to single record", m_optionsHolder );
  tooltip = "All selected samples and/or detectors will be summed together into a single spectrum.<br />"
            "Caution: if you have both foreground and background selected, they will be summed together.";
  HelpSystem::attachToolTipOn( m_sumAllToSingleRecord, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_sumAllToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumToSingleRecordChanged );
  m_sumAllToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumToSingleRecordChanged );
  
  m_sumForeToSingleRecord = new WCheckBox( "Sum Fore to single record", m_optionsHolder );
  tooltip = "There are either multiple sample numbers, or multiple detectors that comprise the"
            " displayed foreground; checking this option will sum all these to a single spectrum"
            " in the output file.";
  HelpSystem::attachToolTipOn( m_sumForeToSingleRecord, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_sumForeToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumForeToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_sumBackToSingleRecord = new WCheckBox( "Sum Back to single record", m_optionsHolder );
  tooltip = "There are either multiple sample numbers, or multiple detectors that comprise the"
            " displayed background; checking this option will sum all these to a single spectrum"
            " in the output file.";
  HelpSystem::attachToolTipOn( m_sumForeToSingleRecord, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_sumBackToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumBackToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_sumSecoToSingleRecord = new WCheckBox( "Sum Sec to single record", m_optionsHolder );
  m_sumSecoToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumSecoToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_backSubFore = new WCheckBox( "Background subtract", m_optionsHolder );
  tooltip = "The output spectrum will be the foreground, minus the background.";
  HelpSystem::attachToolTipOn( m_backSubFore, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_backSubFore->checked().connect( this, &ExportSpecFileTool::handleBackSubForeChanged );
  m_backSubFore->unChecked().connect( this, &ExportSpecFileTool::handleBackSubForeChanged );
  
  m_sumDetsPerSample = new WCheckBox( "Sum det. per sample", m_optionsHolder );
  tooltip = "For each sample number, sum all detectors for that sample number together.<br />"
            "E.g., if you have 8 detectors in your spectrum file, each with 120 samples"
            " then the results will have 120 spectra, each of which is the sum of the"
            " 8 detectors of the respective samples.";
  HelpSystem::attachToolTipOn( m_sumDetsPerSample, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_sumDetsPerSample->checked().connect( this, &ExportSpecFileTool::handleSumDetPerSampleChanged );
  m_sumDetsPerSample->unChecked().connect( this, &ExportSpecFileTool::handleSumDetPerSampleChanged );
  
  
  m_sumSamplesPerDets = new WCheckBox( "Sum samples per det.", m_optionsHolder );
  tooltip = "For each detector, sums all sample numbers together.<br />"
            "E.g., if you have 8 detectors in your spectrum file, each with 120 samples"
            " then the results will have 8 spectra, each of which is the sum of the"
            " 120 samples of their respective detectors.";
  HelpSystem::attachToolTipOn( m_sumSamplesPerDets, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_sumSamplesPerDets->checked().connect( this, &ExportSpecFileTool::handleSumSamplesPerDetChanged );
  m_sumSamplesPerDets->unChecked().connect( this, &ExportSpecFileTool::handleSumSamplesPerDetChanged );
  
  
  m_excludeInterSpecInfo = new WCheckBox( "Remove InterSpec info", m_optionsHolder );
  tooltip = "Removes detector response, peaks, and other analysis"
            " results you may have added in InterSpec.";
  HelpSystem::attachToolTipOn( m_excludeInterSpecInfo, tooltip, true,
                              HelpSystem::ToolTipPosition::Right,
                              HelpSystem::ToolTipPrefOverride::AlwaysShow );
  m_excludeInterSpecInfo->checked().connect( this, &ExportSpecFileTool::handleIncludeInterSpecInfoChanged );
  m_excludeInterSpecInfo->unChecked().connect( this, &ExportSpecFileTool::handleIncludeInterSpecInfoChanged );
  
  m_excludeGpsInfo = new WCheckBox( "Remove GPS", m_optionsHolder );
  
  m_sampleSelectNotAppTxt = new WText( "Not Applicable" );
  m_sampleSelectNotAppTxt->addStyleClass( "ExportNotAppTxt" );
  m_samplesHolder->insertWidget( 1, m_sampleSelectNotAppTxt ); // Put right below title
  
  m_optionsNotAppTxt = new WText( "None Available" );
  m_optionsNotAppTxt->addStyleClass( "ExportNotAppTxt" );
  m_optionsHolder->insertWidget( 1, m_optionsNotAppTxt );      // Put right below title
  
  WContainerWidget *footer = new WContainerWidget( this );
  footer->addStyleClass( "ExportSpecFileFooter" );
  
  
  m_msg = new WText( "&nbsp;", footer );
  m_msg->addStyleClass( "ExportSpecMsg" );
  
  WContainerWidget *btnsDiv = new WContainerWidget( footer );
  btnsDiv->addStyleClass( "ExportSpecBtns" );
  
  
  WPushButton *cancel_btn = new WPushButton( "Cancel", btnsDiv );
  cancel_btn->addStyleClass( "LightButton" );
  cancel_btn->clicked().connect( boost::bind(&ExportSpecFileTool::emitDone, this, false) );
  
  if( !m_resource )
  {
    m_resource = new ExportSpecFileTool_imp::DownloadSpectrumResource( this );
    m_resource->setTakesUpdateLock( true );
    m_resource->downloadFinished().connect( boost::bind(&ExportSpecFileTool::emitDone, this, true) );
  }//if( !m_resource )
  
  
#if( BUILD_AS_OSX_APP || IOS )
  m_export_btn = new WAnchor( WLink(m_resource), btnsDiv );
  m_export_btn->setTarget( AnchorTarget::TargetNewWindow );
  m_export_btn->setStyleClass( "LinkBtn DownloadLink ExportSpecExportBtn" );
#else
  m_export_btn = new WPushButton( btnsDiv );
  m_export_btn->setLink( WLink(m_resource) );
  m_export_btn->setLinkTarget( AnchorTarget::TargetNewWindow );
  m_export_btn->setStyleClass( "LightButton LinkBtn DownloadBtn ExportSpecExportBtn" );
    
#if( ANDROID )
  // Using hacked saving to temporary file in Android, instead of via network download of file.
  m_export_btn->clicked().connect( std::bind( [this](){
    string filename = m_resource->suggestedFileName().toUTF8();
    if( filename.empty() )
      filename = "spectrum_file";
    android_download_workaround( m_resource, filename );
  } ) );
#endif //ANDROID
#endif
  
  m_export_btn->setText( "Export" );
  m_export_btn->disable();


#if( USE_QR_CODES )
  m_show_qr_btn = new WPushButton( "Show QR-code", btnsDiv );
  m_show_qr_btn->clicked().connect( this, &ExportSpecFileTool::handleGenerateQrCode );
  m_show_qr_btn->disable();
  m_show_qr_btn->hide();
#endif
  
  m_formatMenu->select( 0 );
  
  if( !m_specific_spectrum )
    selectDisplayedFile( SpecUtils::SpectrumType::Foreground );
  
  handleFileSelectionChanged();
  updateExportEnabled();
}//void init()


Wt::Signal<bool> &ExportSpecFileTool::done()
{
  return m_done;
}

void ExportSpecFileTool::emitDone( const bool exported )
{
  m_done.emit(exported);
}


void ExportSpecFileTool::updateExportEnabled()
{
  const bool disable = ((m_fileSelect && (m_fileSelect->currentIndex() < 0))
                        && (!m_forePlusBack || m_forePlusBack->isHidden() || !m_forePlusBack->isChecked()));
  m_export_btn->setDisabled( disable );
#if( USE_QR_CODES )
  m_show_qr_btn->setDisabled( disable );
#endif
  
  scheduleAddingUndoRedo();
}//void updateExportEnabled()


void ExportSpecFileTool::updateInfoAboutSelectedFile()
{
  assert( (!m_specific_spectrum) != (!m_fileSelect) );
  
  shared_ptr<const SpecMeas> meas = currentlySelectedFile();
  
  if( !meas )
  {
    // Put a empty table in - just to keep column width
    const char *dummy_tbl = "<table class=\"ExportSpecInfoTable\">\n"
                            "<tr><th>&nbsp;</th><td>&nbsp;</td></tr>\n"
                            "</table>";
    m_fileInfo->setText( dummy_tbl );
    return;
  }
  
  double min_gamma_cps = 0.0, max_gamma_cps = 0.0, min_neutron_cps = 0.0, max_neutron_cps = 0.0;
  SpecUtils::time_point_t start_time = SpecUtils::time_point_t{} + SpecUtils::time_point_t::duration::max();
  SpecUtils::time_point_t end_time = SpecUtils::time_point_t{} + SpecUtils::time_point_t::duration::min();
  assert( SpecUtils::is_special(start_time) && SpecUtils::is_special(end_time) );
  
  const vector<shared_ptr<const SpecUtils::Measurement>> meass = meas->measurements();
  for( const auto m : meass )
  {
    if( !m )
      continue;
    if( (m->source_type() != SpecUtils::SourceType::IntrinsicActivity)
       && !SpecUtils::is_special(m->start_time()) )
    {
      const SpecUtils::time_point_t begin_time = m->start_time();
      const float meas_dur = std::max( m->live_time(), m->real_time() );
      const int64_t meas_dur_ms = static_cast<int64_t>( std::round(1000*meas_dur) );
      const SpecUtils::time_point_t finish_time = begin_time + std::chrono::milliseconds(meas_dur_ms);
      
      start_time = std::min( start_time, begin_time );
      end_time = std::max( end_time, finish_time );
    }
    
    if( (m->num_gamma_channels() >= 1) && ((m->live_time() > 0.0f) || (m->real_time() > 0.0f)) )
    {
      const double nseconds = (m->live_time() > 0.0f) ? m->live_time() : m->real_time();
      const double ncounts = m->gamma_count_sum();
      const double gamma_cps = ncounts / nseconds;
      
      max_gamma_cps = std::max( max_gamma_cps, gamma_cps );
      min_gamma_cps = std::min( min_gamma_cps, gamma_cps );
    }
    
    if( m->contained_neutron() && ((m->live_time() > 0.0f) || (m->real_time() > 0.0f)) )
    {
      const double nseconds = (m->real_time() > 0.0f) ? m->real_time() : m->live_time();
      const double ncounts = m->neutron_counts_sum();
      const double neutron_cps = ncounts / nseconds;
      
      max_neutron_cps = std::max( max_neutron_cps, neutron_cps );
      min_neutron_cps = std::min( min_neutron_cps, neutron_cps );
    }
  }//for( const auto m : meass )
  
  
  WStringStream tabletxt;
  tabletxt << "<table class=\"ExportSpecInfoTable\">\n";
  
  // If the file is a foreground / background type file - let the user know
  const bool is_fore = (meas && (meas == m_interspec->measurment(SpecUtils::SpectrumType::Foreground)));
  const bool is_back = (meas && (meas == m_interspec->measurment(SpecUtils::SpectrumType::Background)));
  const bool is_sec  = (meas && (meas == m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground)));
  
  if( is_fore || is_back || is_sec )
  {
    // We will only show Displayed Sample Numbers if the spectrum is loaded for only one type
    set<int> samples;
    
    tabletxt << "<tr><th>Displayed As</th><td>";
    if( is_fore && !is_back && !is_sec )
    {
      tabletxt << "Foreground";
      samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    }else if( !is_fore && is_back && !is_sec )
    {
      tabletxt << "Background";
      samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    }else if( !is_fore && !is_back && is_sec )
    {
      tabletxt << "Secondary";
      samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
    }else
    {
      string types;
      if( is_fore )
        types += string(types.empty() ? "" : "/") + "Fore";
      if( is_back )
        types += string(types.empty() ? "" : "/") + "Back";
      if( is_sec )
        types += string(types.empty() ? "" : "/") + "Sec";
      tabletxt << types;
    }
    
    tabletxt << "</td></tr>\n";
    
    
    if( !samples.empty() )
      tabletxt << "<tr><th>Disp.&nbsp;Sample" << string((samples.size() > 1) ? "s" : "")
               << "</th><td>" << SpecUtils::sequencesToBriefString( samples )
               << "</td></tr>\n";
  }//if( is_fore || is_back || is_sec )
  
  
  if( !SpecUtils::is_special(start_time) )
  {
    const string start_str = SpecUtils::to_vax_string( start_time );
    size_t pos = start_str.find( ' ' );
    const string start_date = start_str.substr(0, pos);
    const string start_tod = (pos == string::npos) ? string() : start_str.substr(pos+1);
    
    tabletxt << "<tr><th>Start Time</th><td><div>" << start_date << "</div></td></tr>\n"
                "<tr><th></th><td><div>" << start_tod << "</div></td></tr>\n";
    
    if( (end_time - start_time) > 2*std::chrono::minutes(15) )
    {
      const string end_str = SpecUtils::to_vax_string( end_time );
      pos = end_str.find( ' ' );
      const string end_date = end_str.substr(0, pos);
      const string end_tod = (pos == string::npos) ? string() : end_str.substr(pos+1);
      
      tabletxt << "<tr><th>End Time</th><td><div>" << end_date << "</div></td></tr>\n"
                  "<tr><th></th><td><div>" << end_tod << "</div></td></tr>\n";
    }
  }//if( !SpecUtils::is_special(start_time) )
  
  const set<int> &samples = meas->sample_numbers();
  if( samples.size() <= 1 )
  {
    const int nsample = static_cast<int>(meas->sample_numbers().size());
    tabletxt << "<tr><th>Num. Samples</th><td>" << nsample << "</td></tr>\n";
  }else
  {
    tabletxt << "<tr><th>First Sample</th><td>" << (*samples.begin()) << "</td></tr>\n";
    tabletxt << "<tr><th>Last Sample</th><td>" << (*samples.rbegin()) << "</td></tr>\n";
  }
  
  if( meas->contained_neutron() )
  {
    tabletxt << "<tr><th>Neut. Dets</th><td><div>"
             << static_cast<int>(meas->neutron_detector_names().size()) << "</div></td></tr>\n";
    tabletxt << "<tr><th>Gamma Dets</th><td><div>"
            << static_cast<int>(meas->gamma_detector_names().size()) << "</div></td></tr>\n";
  }else
  {
    tabletxt << "<tr><th>Num Dets</th><td><div>"
             << static_cast<int>(meas->detector_names().size()) << "</div></td></tr>\n";
  }
  
  if( meas->has_gps_info() )
  {
    tabletxt << "<tr><th>Latitude</th><td><div>"
             << SpecUtils::printCompact( meas->mean_latitude(), 7 ) << "</div></td></tr>\n";
    tabletxt << "<tr><th>Longitude</th><td><div>"
             << SpecUtils::printCompact( meas->mean_longitude(), 7 ) << "</div></td></tr>\n";
  }//if( meas->has_gps_info() )
  
  const string total_time = PhysicalUnitsLocalized::printToBestTimeUnits( meas->gamma_real_time() );
  tabletxt << "<tr><th>Total Time</th><td><div>" << total_time << "</div></td></tr>\n";
  
  tabletxt << "<tr><th>Sum Gamma</th><td><div>"
           << SpecUtils::printCompact( meas->gamma_count_sum(), 5) << "</div></td></tr>\n";
  if( meas->contained_neutron() )
    tabletxt << "<tr><th>Sum Neut.</th><td><div>"
             << SpecUtils::printCompact( meas->neutron_counts_sum(), 5) << "</div></td></tr>\n";
  
  //const string &instrument_type = meas->instrument_type();
  const string &manufacturer = meas->manufacturer();
  const string &instrument_model = meas->instrument_model();
  const string &instrument_id = meas->instrument_id();
  
  if( !manufacturer.empty() )
    tabletxt << "<tr><th>Manufacturer</th><td><div>"
             << Wt::Utils::htmlEncode(manufacturer) << "</div></td></tr>\n";
  
  if( meas->detector_type() != SpecUtils::DetectorType::Unknown )
    tabletxt << "<tr><th>Model</th><td><div>"
             << SpecUtils::detectorTypeToString(meas->detector_type()) << "</div></td></tr>\n";
  else if( !instrument_model.empty() )
    tabletxt << "<tr><th>Model</th><td><div>"
             << Wt::Utils::htmlEncode(instrument_model) << "</div></td></tr>\n";
  
  if( !instrument_id.empty() )
    tabletxt << "<tr><th>Serial</th><td><div>" << Wt::Utils::htmlEncode(instrument_id) << "</div></td></tr>\n";
  
  tabletxt << "</table>\n";
  
  m_fileInfo->setText( tabletxt.str() );
}//void updateInfoAboutSelectedFile()


void ExportSpecFileTool::handleFileSelectionChanged()
{
  assert( (!m_specific_spectrum) != (!m_fileSelect) );
  
  shared_ptr<const SpecMeas> prev = m_current_file.lock();
  
  shared_ptr<const SpecMeas> current = currentlySelectedFile();
  
  if( prev != current && (!m_forePlusBack || !m_forePlusBack->isChecked()) )
  {
    m_current_file = current;
    
    if( current && (current == m_interspec->measurment(SpecUtils::SpectrumType::Foreground)) )
      m_dispForeSamples->setChecked( true );
    
    if( current && (current == m_interspec->measurment(SpecUtils::SpectrumType::Background)) )
      m_dispBackSamples->setChecked( true );
    
    if( current && (current == m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground)) )
      m_dispSecondSamples->setChecked( true );
  }//if( prev != current )
  
  updateInfoAboutSelectedFile();
  refreshSampleAndDetectorOptions();
  updateExportEnabled();
}//void ExportSpecFileTool::handleFileSelectionChanged()


void ExportSpecFileTool::selectDisplayedFile( const SpecUtils::SpectrumType type )
{
  assert( (!m_specific_spectrum) != (!m_fileSelect) );
  
  if( m_specific_spectrum )
    return;
  
  const shared_ptr<SpecMeas> spec = m_interspec->measurment(type);
  const set<int> &samples = m_interspec->displayedSamples(type);
  const vector<string> detectors = m_interspec->detectorsToDisplay(type);
  
  const int orig_selection = m_fileSelect->currentIndex();
  
  SpecMeasManager * const measManager = m_interspec->fileManager();
  SpectraFileModel * const fileModel = measManager ? measManager->model() : nullptr;
  
  int updated_selection = -1;
  const int nfiles = fileModel ? fileModel->rowCount() : 0;
  for( int row = 0; (updated_selection < 0) && (row < nfiles); ++row )
  {
    shared_ptr<SpectraFileHeader> header = fileModel->fileHeader(row);
    assert( header );
    shared_ptr<SpecMeas> meas = header ? header->measurementIfInMemory() : nullptr;
    if( meas == spec )
      updated_selection = row;
  }//for( int row = 0; row < nfiles; ++row )
  
  if( updated_selection != orig_selection )
  {
    m_fileSelect->setCurrentIndex( updated_selection );
    handleFileSelectionChanged();
  }//if( updated_selection != orig_selection )
  
  if( updated_selection < 0 )
    return;
  
  // Update sample numbers and detectors here
  refreshSampleAndDetectorOptions();
}//void selectDisplayedFile( const SpecUtils::SpectrumType type )


pair<set<int>,string> ExportSpecFileTool::sampleNumbersFromTxtRange( string fulltxt,
                                                       shared_ptr<const SpecMeas> meas,
                                                       const bool quiet )
{
  set<int> answer;
  
  if( !meas )
    return {answer, ""};
  
  //Lets replace numbers separated by spaces, to be separated by commas.
  //  We cant do a simple replace of spaces to commas, and using a regex would
  //  require using a lookahead or behind
  //  Note that below while loop is a little inefficient, but whatever
  std::smatch mtch;
  std::regex expr( ".*(\\d\\s+\\d).*" );
  while( std::regex_match( fulltxt, mtch, expr ) )
    fulltxt[mtch[1].first - fulltxt.begin() + 1] = ',';
  
//Using a sregex_token_iterator doesnt work for single digit numbers, ex
//  '1 2 3 45 983 193'  will go to '1,2 3,45,983,193'
//  boost::regex expr( "\\d\\s+\\d" );
//  for( boost::sregex_token_iterator iter(fulltxt.begin(),fulltxt.end(),expr,0);
//      iter != boost::sregex_token_iterator(); ++iter )
//    fulltxt[iter->first - fulltxt.begin()+1] = ',';
  
  vector<string> sampleranges;
  SpecUtils::split( sampleranges, fulltxt, "," );

  const set<int> &samples = meas->sample_numbers();

  try
  {
    for( string textStr : sampleranges )
    {
      SpecUtils::trim( textStr );
      if( textStr.empty() )
        continue;
      
      std::smatch matches;
      std::regex range_expression( "(\\d+)\\s*(\\-|to|through)\\s*(\\d+)",
                                  std::regex::ECMAScript | std::regex::icase );
      if( std::regex_match( textStr, matches, range_expression ) )
      {
        string firstStr = string( matches[1].first, matches[1].second );
        string lastStr = string( matches[3].first, matches[3].second );
        
        int first = std::stoi( firstStr );
        int last = std::stoi( lastStr );
        if( last < first )
          std::swap( last, first );
        
        const set<int>::const_iterator firstPos = samples.find(first);
        set<int>::const_iterator lastPos = samples.find(last);
        
        if( firstPos == samples.end() )
          throw runtime_error( "Sample number " + firstStr + " doesnt exist." );
        
        if( lastPos == samples.end() )
          throw runtime_error( "Sample number " + lastStr + " doesnt exist." );
        
        ++lastPos;
        answer.insert( firstPos, lastPos );
      }else
      {
        const int sample = std::stoi( textStr );
        
        if( !samples.count(sample) )
          throw runtime_error( "Sample number " + std::to_string(sample) + " does not exist in the measurement." );
        answer.insert( sample );
      }//if( is sample range ) / else
    }//for( string textStr : sampleranges )
  }catch( std::exception &e )
  {
    if( !quiet )
    {
      passMessage( "Error converting sample numbers: "
                  + string(e.what()), WarningWidget::WarningMsgHigh );
    }
    
    return { set<int>{}, "" };
  }//try / catch
  
  
  const set<int> &total_sample_nums = meas->sample_numbers();

  if( answer.size() == 0 || samples.size() == 0 )
    fulltxt = "";
  else if( answer.size() == 1 )
    fulltxt = std::to_string( *begin(answer) );
  else if( answer.size() == samples.size() )
    fulltxt = to_string( *begin(answer) ) + "-" + to_string(*answer.rbegin());
  else
    fulltxt = SpecUtils::sequencesToBriefString( answer );
  
  return { answer, fulltxt };
}//pair<set<int>,string> sampleNumbersFromTxtRange( )
 

std::shared_ptr<const SpecMeas> ExportSpecFileTool::currentlySelectedFile() const
{
  if( m_specific_spectrum )
    return m_specific_spectrum;
  
  
  if( m_forePlusBack && m_forePlusBack->isChecked() )
  {
    shared_ptr<const SpecMeas> fore = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecMeas> back = m_interspec->measurment(SpecUtils::SpectrumType::Background);
    shared_ptr<const SpecMeas> seco = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
    
    assert( fore );
    if( !fore )
      return nullptr;
    
    shared_ptr<SpecMeas> result = make_shared<SpecMeas>();
    result->uniqueCopyContents( *fore );
    result->remove_measurements( result->measurements() );
    
    for( const set<int> &samples : result->sampleNumsWithPeaks() )
      result->setPeaks( {}, samples );
    
    const set<int> &fore_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    const set<int> &back_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    const set<int> &seco_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
    
    const vector<string> fore_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
    const vector<string> back_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Background);
    const vector<string> seco_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::SecondForeground);
    
    const auto add_meas = [&result]( const shared_ptr<const SpecMeas> &meas,
                                                    const set<int> &samples,
                                                    const vector<string> &dets,
                                                    const SpecUtils::SpectrumType type ){
      if( !meas )
        return;
      
      int current_sample = result->sample_numbers().empty() ? 0 : (*result->sample_numbers().rbegin());
      set<int> used_samples, added_samples;
      
      for( const int sample : samples )
      {
        for( const string &det : dets )
        {
          const shared_ptr<const SpecUtils::Measurement> m = meas->measurement( sample, det );
          if( !m )
            continue;
          
          if( !used_samples.count(m->sample_number()) )
          {
            current_sample += 1;
            used_samples.insert( m->sample_number() );
          }
          
          added_samples.insert( current_sample );
          shared_ptr<SpecUtils::Measurement> new_m = make_shared<SpecUtils::Measurement>( *m );
          new_m->set_sample_number( current_sample );
          
          switch( type )
          {
            case SpecUtils::SpectrumType::Foreground:
              new_m->set_source_type( SpecUtils::SourceType::Foreground );
              break;
              
            case SpecUtils::SpectrumType::SecondForeground:
              new_m->set_source_type( SpecUtils::SourceType::Unknown );
              break;
              
            case SpecUtils::SpectrumType::Background:
              new_m->set_source_type( SpecUtils::SourceType::Background );
              break;
          }//switch( type )
          
          result->add_measurement( new_m, false );
        }//for( const string &det : dets )
      }//for( const int sample : fore_samples )
      
      shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = meas->peaks( samples );
      if( peaks )
        result->setPeaks( *peaks, added_samples );
    };//add_meas lamda
    
    add_meas( fore, fore_samples, fore_dets, SpecUtils::SpectrumType::Foreground );
    add_meas( back, back_samples, back_dets, SpecUtils::SpectrumType::Background );
    add_meas( seco, seco_samples, seco_dets, SpecUtils::SpectrumType::SecondForeground );
    
    result->cleanup_after_load( SpecUtils::SpecFile::CleanupAfterLoadFlags::DontChangeOrReorderSamples );
    
    return result;
  }//if( m_forePlusBack && m_forePlusBack->isChecked() )
  
  
  const int index = m_fileSelect->currentIndex();
  const SpecMeasManager * const measManager = m_interspec->fileManager();
  const SpectraFileModel * const fileModel = measManager ? measManager->model() : nullptr;
  const shared_ptr<const SpectraFileHeader> header = fileModel ? fileModel->fileHeader(index) : nullptr;
  const shared_ptr<const SpecMeas> spec = header ? header->parseFile() : nullptr;
  
  return spec;
}//std::shared_ptr<const SpecMeas> currentlySelectedFile() const


set<int> ExportSpecFileTool::currentlySelectedSamples() const
{
  shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  if( !spec )
    return set<int>{};
  
  if( m_specific_spectrum )
    return m_specific_samples.empty() ? m_specific_spectrum->sample_numbers() : m_specific_samples;
  
  const set<int> &samples = spec->sample_numbers();
  if( (samples.size() == 1) || (m_allSamples->isVisible() && m_allSamples->isChecked()) )
    return samples;
  
  if( m_customSamples->isVisible() && m_customSamples->isChecked() )
  {
    const string txt = m_customSamplesEdit->text().toUTF8();
    pair<set<int>,string> samplenums = sampleNumbersFromTxtRange( txt, spec, true );
    return samplenums.first;
  }
  
  if( m_forePlusBack && m_forePlusBack->isVisible() && m_forePlusBack->isChecked() )
    return spec->sample_numbers();
  
  const bool disp_fore = (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked());
  const bool disp_back = (m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked());
  const bool disp_seco = (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked());
  
  assert( disp_fore || disp_back || disp_seco );
  if( !disp_fore && !disp_back && !disp_seco )
  {
    passMessage( "ExportSpecFileTool::encodeStateToUrl: logic error - will use all sample numbers", WarningWidget::WarningMsgHigh );
    return spec->sample_numbers();
  }
  
  set<int> answer;
  if( disp_fore )
  {
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    answer.insert( begin(samples), end(samples) );
  }
  
  if( disp_back )
  {
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    answer.insert( begin(samples), end(samples) );
  }
  
  if( disp_seco )
  {
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
    answer.insert( begin(samples), end(samples) );
  }
  
  assert( !answer.empty() );
  
  return answer;
}//std::set<int> currentlySelectedSamples() const;


vector<string> ExportSpecFileTool::currentlySelectedDetectors() const
{
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  if( !spec )
    return vector<string>{};
  
  if( m_specific_spectrum )
    return m_specific_detectors.empty() ? m_specific_spectrum->detector_names() : m_specific_detectors;
  
  if( spec->gamma_detector_names().size() <= 1 )
    return spec->detector_names();
  
  if( !m_filterDetector || m_filterDetector->isHidden() || !m_filterDetector->isChecked() )
    return spec->detector_names();
  
  map<string,string> label_to_orig;
  for( const string &name : spec->detector_names() )
  {
    // TODO: use Wt::Utils::htmlEncode, if that is what happens in WCheckBox, instead of creating a WCheckBox
    WCheckBox cb( name );
    label_to_orig[cb.text().toUTF8()] = name;
  }
  
  vector<string> answer;
  const vector<Wt::WWidget *> &detector_cbs = m_detectorFilterCbs->children();
  for( const auto w : detector_cbs )
  {
    WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
    
    if( !cb || !cb->isChecked() )
      continue;
    
    const string cb_label = cb->text().toUTF8();
    const auto pos = label_to_orig.find( cb_label );
    assert( pos != end(label_to_orig) );
    if( pos != end(label_to_orig) )
      answer.push_back( pos->second );
    else
      throw runtime_error( "ExportSpecFileTool::currentlySelectedDetectors:"
                          " Error matching detector names." );
  }//for( const auto w : m_detectorFilterCbs->children() )
  
  return answer;
}//vector<string> currentlySelectedDetectors() const


SpecUtils::SaveSpectrumAsType ExportSpecFileTool::currentSaveType() const
{
  const WMenuItem * const currentFormatItem = m_formatMenu->currentItem();
  assert( currentFormatItem );
  
  if( currentFormatItem )
  {
    const uint64_t data = reinterpret_cast<uint64_t>( currentFormatItem->data() );
    assert( data < static_cast<int>(SpecUtils::SaveSpectrumAsType::NumTypes) );
    return SpecUtils::SaveSpectrumAsType( data );
  }//if( currentFormatItem )
  
  return SpecUtils::SaveSpectrumAsType::N42_2012;
}//SpecUtils::SaveSpectrumAsType currentSaveType() const;


uint16_t ExportSpecFileTool::maxRecordsInCurrentSaveType( shared_ptr<const SpecMeas> spec ) const
{
  const SpecUtils::SaveSpectrumAsType save_format = currentSaveType();
  
  switch( save_format )
  {
    // Spectrum file types that can have many spectra in them
    case SpecUtils::SaveSpectrumAsType::Txt:
    case SpecUtils::SaveSpectrumAsType::Csv:
    case SpecUtils::SaveSpectrumAsType::Pcf:
    case SpecUtils::SaveSpectrumAsType::N42_2006:
    case SpecUtils::SaveSpectrumAsType::N42_2012:
    case SpecUtils::SaveSpectrumAsType::ExploraniumGr130v0:
    case SpecUtils::SaveSpectrumAsType::ExploraniumGr135v2:
      return std::numeric_limits<uint16_t>::max();
      
#if( USE_QR_CODES )
    // Spectrum file types that can have two spectra in them (foreground + background)
    case SpecUtils::SaveSpectrumAsType::Uri:
      if( !spec )
        spec = currentlySelectedFile();
      return (spec && (spec->num_gamma_channels() < 2075)) ? 2 : 1;
#endif
      
#if( SpecUtils_ENABLE_D3_CHART )
    case SpecUtils::SaveSpectrumAsType::HtmlD3:
      return 2;
#endif
      
    // Spectrum file types that can have a single spectra in them
    case SpecUtils::SaveSpectrumAsType::Chn:
    case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
    case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
    case SpecUtils::SaveSpectrumAsType::SpcAscii:
    case SpecUtils::SaveSpectrumAsType::SpeIaea:
    case SpecUtils::SaveSpectrumAsType::Cnf:
    case SpecUtils::SaveSpectrumAsType::Tka:
      return 1;
    
#if( SpecUtils_INJA_TEMPLATES )
    case SpecUtils::SaveSpectrumAsType::Template:
      assert( 0 );
      break;
#endif
    
    case SpecUtils::SaveSpectrumAsType::NumTypes:
      assert( 0 );
      break;
  }//switch( save_format )
  
  assert( 0 );
  
  return 0;
}//uint16_t maxRecordsInCurrentSaveType() const


bool ExportSpecFileTool::removeInterSpecInfo() const
{
  return (m_excludeInterSpecInfo->isVisible() && m_excludeInterSpecInfo->isChecked());
}//bool removeInterSpecInfo() const


void ExportSpecFileTool::refreshSampleAndDetectorOptions()
{
  scheduleAddingUndoRedo();
  
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  const uint16_t max_records = maxRecordsInCurrentSaveType( spec );
  const SpecUtils::SaveSpectrumAsType save_type = currentSaveType();
  
  
  {// Begin update suggested spectrum file name
    string filename = spec ? SpecUtils::filename( spec->filename() ) : string("spectrum");
    const string orig_ext = SpecUtils::file_extension(filename);
    if( (orig_ext.size() > 0) && (orig_ext.size() <= 4) && (orig_ext.size() < filename.size()) )
      filename = filename.substr( 0, filename.size() - orig_ext.size() );
    filename += string(".") + SpecUtils::suggestedNameEnding(save_type);
    m_resource->suggestFileName( filename );
  }// End update suggested spectrum file name
  
  if( !m_dispForeSamples->isChecked()
     && !m_dispBackSamples->isChecked()
     && !m_dispSecondSamples->isChecked()
     && !m_customSamples->isChecked()
     && !m_allSamples->isChecked() )
  {
    m_allSamples->setChecked( true );
  }
  
  
  if( !spec || ((spec->gamma_detector_names().size() <= 1) && (spec->sample_numbers().size() <= 1)) )
  {
    m_filterDetector->setChecked( false );
    m_detectorFilterCbs->clear();
    m_filterDetector->hide();
    m_detectorFilterCbs->hide();
    
    m_sumAllToSingleRecord->hide();
    m_sumForeToSingleRecord->hide();
    m_sumBackToSingleRecord->hide();
    m_sumSecoToSingleRecord->hide();
  }else
  {
    map<string,bool> prev_check;
    for( const auto w : m_detectorFilterCbs->children() )
    {
      WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
      if( cb )
        prev_check[cb->text().toUTF8()] = cb->isChecked();
    }
    
    m_filterDetector->show();
    m_detectorFilterCbs->setHidden( !m_filterDetector->isChecked() );
    
    m_detectorFilterCbs->clear();
    for( const string &name : spec->detector_names() )
    {
      WCheckBox *cb = new WCheckBox( name, m_detectorFilterCbs );
      if( prev_check.count(cb->text().toUTF8()) )
        cb->setChecked( prev_check[cb->text().toUTF8()] );
      else
        cb->setChecked( true );  //Could check if spectrum is displayed, and if so if the det is displayed
      
      // Curiously, if we dont have these next two calls, then the checkboxes wont actually be
      //  registered and not-checked (after the user unchecks them) in
      //  `ExportSpecFileTool::currentlySelectedDetectors()`, if the user immediately clicks
      //  "Export" after un-checking a detector; `ExportSpecFileTool::handleDetectorsToFilterChanged()`
      //  doesnt even have to do anything, just any function to be called.
      cb->checked().connect( this, &ExportSpecFileTool::handleDetectorsToFilterChanged );
      cb->unChecked().connect( this, &ExportSpecFileTool::handleDetectorsToFilterChanged );
    }//
    
    if( (max_records <= 1) || (max_records == 2) || !spec )
    {
      m_sumAllToSingleRecord->hide();
      m_sumForeToSingleRecord->hide();
      m_sumBackToSingleRecord->hide();
      m_sumSecoToSingleRecord->hide();
    }else
    {
      m_sumAllToSingleRecord->show();
      
      const bool forePlusBack = (m_forePlusBack && m_forePlusBack->isVisible() && m_forePlusBack->isChecked());
      const bool useDispForeSamples = (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked());
      const bool useDispBackSamples = (m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked());
      const bool useDispSecondSamples = (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked());
      
      size_t numForeRecords = 0, numBackRecords = 0, numSecRecords = 0;
      if( forePlusBack )
      {
        for( const auto &m : spec->measurements() )
        {
          switch( m->source_type() )
          {
            case SpecUtils::SourceType::IntrinsicActivity: assert( 0 );       break;
            case SpecUtils::SourceType::Calibration:       assert( 0 );       break;
            case SpecUtils::SourceType::Foreground:        ++numForeRecords;  break;
            case SpecUtils::SourceType::Background:        ++numBackRecords;  break;
            case SpecUtils::SourceType::Unknown:           ++numSecRecords;   break;
          }//switch( meass[i]->source_type() )
        }//for( size_t i = 0; (!hasBackground || !hasSecondary) && (i < meass.size()); ++i )
      }else //if( forePlusBack )
      {
        using SpecUtils::SpectrumType;
        for( auto t : {SpectrumType::Foreground, SpectrumType::Background,
                        SpectrumType::SecondForeground} )
        {
          size_t nrecords = 0;
          const vector<string> dets = m_interspec->detectorsToDisplay(t);
          const set<int> &samples = m_interspec->displayedSamples(t);
          for( const int sample : samples )
          {
            for( const string &det : dets )
              nrecords += spec->measurement(sample, det) ? 1 : 0;
          }//for( const int sample : samples )
          
          switch ( t )
          {
            case SpectrumType::Foreground:       numForeRecords = nrecords; break;
            case SpectrumType::SecondForeground: numSecRecords  = nrecords; break;
            case SpectrumType::Background:       numBackRecords = nrecords; break;
          }
        }//for( auto t : {SpectrumType::Foreground, SpectrumType::Background, SpectrumType::SecondForeground} )
      }//if( forePlusBack ) / else
      
      const bool showForeToSingle = ((useDispForeSamples || forePlusBack) && (numForeRecords > 1));
      const bool showBackToSingle = ((useDispBackSamples || forePlusBack) && (numBackRecords > 1));
      const bool showSecondToSingle = ((useDispSecondSamples || forePlusBack) && (numSecRecords > 1));
      
      m_sumForeToSingleRecord->setHidden( !showForeToSingle );
      m_sumBackToSingleRecord->setHidden( !showBackToSingle );
      m_sumSecoToSingleRecord->setHidden( !showSecondToSingle );
    }
  }//if( spec->gamma_detector_names().size() <= 1 ) / else
  
  
  const bool use_fore_disp = (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked());
  const bool use_seco_disp = (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked());
  const bool use_back_disp = (m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked());
  
  if( (use_fore_disp || use_seco_disp) && use_back_disp )
  {
    m_backSubFore->show();
    //if( use_fore_disp && use_seco_disp )
    //  m_backSubFore->setText( "Back. sub. For./Sec." );
    //else if( use_fore_disp )
    //  m_backSubFore->setText( "Back. Sub. For." );
    //else
    //  m_backSubFore->setText( "Back. Sub. Sec." );
  }else
  {
    m_backSubFore->hide();
  }
  
  const bool can_save_interspec_info = ((save_type == SpecUtils::SaveSpectrumAsType::N42_2012)
                                        || (save_type == SpecUtils::SaveSpectrumAsType::N42_2006));
  m_excludeInterSpecInfo->setHidden( !can_save_interspec_info );
  
  
  if( !spec || (max_records <= 2) || (spec->gamma_detector_names().size() <= 1) )
  {
    m_sumDetsPerSample->hide();
  }else
  {
    m_sumDetsPerSample->show();
  }
  
  if( !spec || (max_records <= 2)
     || (spec->sample_numbers().size() <= 1)
     || (spec->gamma_detector_names().size() <= 1) )
  {
    m_sumSamplesPerDets->hide();
  }else
  {
    m_sumSamplesPerDets->show();
  }
  
  
  m_excludeGpsInfo->setHidden( !spec || !spec->has_gps_info() );
  
  if( !spec || (spec->sample_numbers().size() <= 1) )
  {
    m_dispForeSamples->setChecked(false);
    m_dispBackSamples->setChecked(false);
    m_dispSecondSamples->setChecked(false);
    m_allSamples->setChecked(true);
    m_customSamples->setChecked(false);
    m_customSamplesEdit->hide();
    m_customSamplesEdit->setText( "" );
    
    m_dispForeSamples->hide();
    m_dispBackSamples->hide();
    m_dispSecondSamples->hide();
    m_allSamples->hide();
    m_customSamples->hide();
    m_customSamplesEdit->hide();
    
    m_backSubFore->hide();
    m_sumDetsPerSample->hide();
    m_sumSamplesPerDets->hide();
    
    m_sampleSelectNotAppTxt->show();
    m_optionsNotAppTxt->setHidden( m_excludeGpsInfo->isVisible()
                                  || m_excludeInterSpecInfo->isVisible() );
    
    if( spec && (spec->gamma_detector_names().size() > 1) && (max_records < 2) )
    {
      m_msg->setText( "Detectors will be summed together." );
    }else
    {
      m_msg->setText( "&nbsp;" );
    }
    
    return;
  }//if( no reason to show selecting samples )
  
  
  const shared_ptr<const SpecMeas> foreground = m_interspec->measurment(SpecUtils::SpectrumType::Foreground);
  const shared_ptr<const SpecMeas> background = m_interspec->measurment(SpecUtils::SpectrumType::Background);
  const shared_ptr<const SpecMeas> secondary = m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground);
  
  const bool is_for = (spec && (spec == foreground));
  const bool is_back = (spec && (spec == background));
  const bool is_second = (spec && (spec == secondary));
  
  
  const vector<string> fore_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Foreground);
  const vector<string> back_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::Background);
  const vector<string> sec_dets = m_interspec->detectorsToDisplay(SpecUtils::SpectrumType::SecondForeground);
  
  const set<int> &fore_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
  const set<int> &back_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
  const set<int> &sec_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
  
  
  m_dispForeSamples->setHidden( !is_for );
  m_dispForeSamples->setChecked( m_dispForeSamples->isChecked() && is_for );
  
  m_dispBackSamples->setHidden( !is_back );
  m_dispBackSamples->setChecked( m_dispBackSamples->isChecked() && is_back );
  
  m_dispSecondSamples->setHidden( !is_second );
  m_dispSecondSamples->setChecked( m_dispSecondSamples->isChecked() && is_second );

    
  const set<int> &sample = spec->sample_numbers();
  
  if( (sample.size() <= 1) || (m_forePlusBack && m_forePlusBack->isChecked()) )
  {
    m_allSamples->setChecked( true );
    m_allSamples->setHidden( true );
    m_customSamples->setChecked( false );
    m_customSamplesEdit->setHidden( true );
  }else
  {
    m_allSamples->setHidden( false );
    m_customSamples->setHidden( false );
    m_customSamplesEdit->setHidden( !m_customSamples->isChecked() );
  }
  
  if( m_customSamples->isChecked() )
  {
    const string txt = m_customSamplesEdit->text().toUTF8();
    const pair<set<int>,string> samples_txt = sampleNumbersFromTxtRange( txt, spec, true );
    
    m_customSamplesEdit->setText( samples_txt.second );
    
    if( samples_txt.first.empty() )
    {
      m_customSamples->setChecked( false );
      m_customSamplesEdit->hide();
    }
  }//if( m_customSamples->isChecked() )
  
  if( m_dispForeSamples->isChecked()
     || m_dispBackSamples->isChecked()
     || m_dispSecondSamples->isChecked()
     || m_customSamples->isChecked() )
  {
    m_allSamples->setChecked( false );
  }
  
  
  const set<int> samplesToUse = currentlySelectedSamples();
  const vector<string> detsToUse = currentlySelectedDetectors();
  
  if( (max_records >= 2) && ((samplesToUse.size() > 1) || (spec->gamma_detector_names().size() > 1)) )
  {
    m_sumAllToSingleRecord->show();
  }else
  {
    m_sumAllToSingleRecord->hide();
  }
  
  // Loop over sample numbers, and see if there is more than one detector for any sample number,
  //  and use this info to show/hide m_sumDetsPerSample
  bool mult_dets_per_sample = false;
  for( const int sample : samplesToUse )
  {
    size_t num_dets = 0;
    for( const string &det : detsToUse )
      num_dets += !!spec->measurement( sample, det );
    mult_dets_per_sample = (num_dets > 1);
    if( mult_dets_per_sample )
      break;
  }//for( const int sample : samplesToUse )
  
  // Loop over detectors, and see if there is more than one sample number for any detecotr,
  //  and use this info to show/hide m_sumSamplesPerDets
  bool mult_sample_per_dets = false;
  for( const string &det : detsToUse )
  {
    size_t num_samples = 0;
    for( const int sample : samplesToUse )
      num_samples += !!spec->measurement( sample, det );
      
    mult_sample_per_dets = (num_samples > 1);
    if( mult_sample_per_dets )
      break;
  }//for( const string &det : detsToUse )
  
  m_sumDetsPerSample->setHidden( (max_records <= 2) || !mult_dets_per_sample );
  m_sumSamplesPerDets->setHidden( (max_records <= 2) || !mult_sample_per_dets );
  m_filterDetector->setHidden( !mult_dets_per_sample );
  m_detectorFilterCbs->setHidden( !mult_dets_per_sample || !m_filterDetector->isChecked() );
  
  size_t num_sample_showing = 0, num_option_showing = 0;
  for( const auto w : m_samplesHolder->children() )
  {
    auto cb = dynamic_cast<WCheckBox *>( w );
    num_sample_showing += (cb && cb->isVisible());
  }
  
  for( const auto w : m_optionsHolder->children() )
  {
    auto cb = dynamic_cast<WCheckBox *>( w );
    num_option_showing += (cb && cb->isVisible());
  }
  
  m_sampleSelectNotAppTxt->setHidden( num_sample_showing );
  m_optionsNotAppTxt->setHidden( num_option_showing );
  
  
  
  if( (max_records < 2)
     && ( (spec->gamma_detector_names().size() > 1) || (samplesToUse.size() > 1 ) ) )
  {
    if( m_backSubFore->isVisible() && m_backSubFore->isEnabled() && m_backSubFore->isChecked()
       && (use_fore_disp != use_seco_disp) )
    {
      m_msg->setText( "A single spectrum will be produced." );
    }else
    {
      m_msg->setText( "Records will be summed together." );
    }
  }else if( max_records == 2 )
  {
    size_t num_records = 0;
    for( const int sample : samplesToUse )
    {
      for( const string &det : detsToUse )
        num_records += ( spec->measurement(sample, det) ? 1 : 0);
    }//for( const int samples : samplesToUse )
    
    // QR code here
    if( use_fore_disp && use_seco_disp && use_back_disp )
    {
      m_msg->setText( "Will be summed to single spec." );
    }else if( (use_fore_disp || use_seco_disp) && use_back_disp
       && (!m_sumAllToSingleRecord->isVisible() || !m_sumAllToSingleRecord->isChecked()) )
    {
      m_msg->setText( "QR will have 2 spectrum" );
    }else if( m_backSubFore->isVisible() && m_backSubFore->isEnabled() && m_backSubFore->isChecked()
             && (use_fore_disp != use_seco_disp) )
    {
      m_msg->setText( "A single spectrum will be produced." );
    }else if( num_records > 2 )
    {
      m_msg->setText( "Records will be summed together." );
    }else
    {
      m_msg->setText( "&nbsp;" );
    }
  }else
  {
    m_msg->setText( "&nbsp;" );
  }
  
  scheduleAddingUndoRedo();
}//void refreshSampleAndDetectorOptions()


void ExportSpecFileTool::handleAllSampleChanged()
{
  if( m_allSamples->isChecked() )
  {
    m_dispForeSamples->setChecked( false );
    m_dispBackSamples->setChecked( false );
    m_dispSecondSamples->setChecked( false );
    m_customSamples->setChecked( false );
    m_customSamplesEdit->hide();
  }else
  {
    // ...
  }
  
  handleSamplesChanged();
}//void handleAllSampleChanged();


void ExportSpecFileTool::handleDisplaySampleChanged( const SpecUtils::SpectrumType type )
{
  scheduleAddingUndoRedo();
  
  Wt::WCheckBox *cb = nullptr;
  switch( type )
  {
    case SpecUtils::SpectrumType::Foreground:
      cb = m_dispForeSamples;
      break;
      
    case SpecUtils::SpectrumType::SecondForeground:
      cb = m_dispSecondSamples;
      break;
      
    case SpecUtils::SpectrumType::Background:
      cb = m_dispBackSamples;
      break;
  }//switch( type )
  
  if( cb && cb->isChecked() )
  {
    m_allSamples->setChecked( false );
    m_customSamples->setChecked( false );
    m_customSamplesEdit->hide();
  }
  
  handleSamplesChanged();
}//void handleDisplaySampleChanged( const SpecUtils::SpectrumType type );


void ExportSpecFileTool::handleCustomSampleChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_customSamples->isChecked() )
  {
    m_customSamplesEdit->show();
    const string txt = m_customSamplesEdit->text().toUTF8();
    
    const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
    const pair<set<int>,string> samples_txt = sampleNumbersFromTxtRange( txt, spec, true );
    m_customSamplesEdit->setText( samples_txt.second );
    if( samples_txt.second.empty() && spec && !spec->sample_numbers().empty() )
    {
      const set<int> &samples = spec->sample_numbers();
      const string txt = to_string(*samples.begin()) + "-" + to_string(*samples.rbegin());
      m_customSamplesEdit->setText( txt );
    }
    
    m_dispForeSamples->setChecked( false );
    m_dispBackSamples->setChecked( false );
    m_dispSecondSamples->setChecked( false );
    m_allSamples->setChecked( false );
  }else
  {
    m_customSamplesEdit->hide();
    m_dispForeSamples->setChecked( false );
    m_dispBackSamples->setChecked( false );
    m_dispSecondSamples->setChecked( false );
    m_allSamples->setChecked( true );
  }
  
  handleSamplesChanged();
}//void handleCustomSampleChanged()


void ExportSpecFileTool::handleCustomSampleTxtChanged()
{
  scheduleAddingUndoRedo();
  
  const string txt = m_customSamplesEdit->text().toUTF8();
  
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  const pair<set<int>,string> samples_txt = sampleNumbersFromTxtRange( txt, spec, false );
  m_customSamplesEdit->setText( samples_txt.second );
  
  handleSamplesChanged();
}//void handleCustomSampleTxtChanged()


void ExportSpecFileTool::handleSamplesChanged()
{
  scheduleAddingUndoRedo();
  refreshSampleAndDetectorOptions();
}//void handleSamplesChanged()


void ExportSpecFileTool::handleFormatChange()
{
  scheduleAddingUndoRedo();
  
  const SpecUtils::SaveSpectrumAsType save_format = currentSaveType();
  
#if( USE_QR_CODES )
  const bool is_qr = (save_format == SpecUtils::SaveSpectrumAsType::Uri);
  m_export_btn->setHidden( is_qr );
  m_show_qr_btn->setHidden( !is_qr );
#else
  m_export_btn->show();
#endif
  
  refreshSampleAndDetectorOptions();
}//void handleFormatChange();
  

void ExportSpecFileTool::handleForePlusBackChanged()
{
  scheduleAddingUndoRedo();
  
  assert( m_forePlusBack && m_fileSelect );
  if( !m_forePlusBack || !m_fileSelect )
    return;
  
  if( m_forePlusBack->isChecked() )
  {
    m_fileSelect->setCurrentIndex( -1 );
    m_fileSelect->disable();
  }else
  {
    m_fileSelect->enable();
  }
  
  handleFileSelectionChanged();
}//void ExportSpecFileTool::handleForePlusBackChanged()


void ExportSpecFileTool::handleFilterDetectorCbChanged()
{
  scheduleAddingUndoRedo();
  
  assert( m_filterDetector );
  assert( m_detectorFilterCbs );
  if( !m_filterDetector || !m_detectorFilterCbs )
    return;
  
  const bool use_all_det = (m_filterDetector->isHidden() || !m_filterDetector->isChecked());
  m_detectorFilterCbs->setHidden( use_all_det );
  
  if( use_all_det )
  {
    for( const auto w : m_detectorFilterCbs->children() )
    {
      WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
      if( cb )
        cb->setChecked( true );
    }
  }//if( !m_filterDetector->isChecked() )
  
  // TODO: update options
}//void handleFilterDetectorCbChanged()


void ExportSpecFileTool::handleDetectorsToFilterChanged()
{
  scheduleAddingUndoRedo();
}//void handleDetectorsToFilterChanged();


void ExportSpecFileTool::handleSumToSingleRecordChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_sumAllToSingleRecord->isChecked() )
  {
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
    m_sumSamplesPerDets->setChecked( false );
  }//if( m_sumAllToSingleRecord->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleSumToSingleRecordChanged()


void ExportSpecFileTool::handleSumTypeToSingleRecordChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_sumForeToSingleRecord->isChecked()
     || m_sumBackToSingleRecord->isChecked()
     || m_sumSecoToSingleRecord->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
    m_sumSamplesPerDets->setChecked( false );
  }
  
  refreshSampleAndDetectorOptions();
}//void handleSumTypeToSingleRecordChanged()


void ExportSpecFileTool::handleBackSubForeChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_backSubFore->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
    m_sumSamplesPerDets->setChecked( false );
  }//if( m_backSubFore->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleBackSubForeChanged()


void ExportSpecFileTool::handleSumDetPerSampleChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_sumDetsPerSample->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumSamplesPerDets->setChecked( false );
  }//if( m_sumDetsPerSample->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleSumDetPerSampleChanged()


void ExportSpecFileTool::handleSumSamplesPerDetChanged()
{
  scheduleAddingUndoRedo();
  
  if( m_sumSamplesPerDets->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
  }//if( m_sumSamplesPerDets->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleSumSamplesPerDetChanged()


void ExportSpecFileTool::handleIncludeInterSpecInfoChanged()
{
  // Nothing to do here I think
  scheduleAddingUndoRedo();
}//void handleIncludeInterSpecInfoChanged()


std::shared_ptr<const SpecMeas> ExportSpecFileTool::generateFileToSave()
{
  const shared_ptr<const SpecMeas> start_spec = currentlySelectedFile();
  const SpecUtils::SaveSpectrumAsType save_type = currentSaveType();
  const uint16_t max_records = maxRecordsInCurrentSaveType( start_spec );
  
  const bool remove_gps = (m_excludeGpsInfo->isVisible() && m_excludeGpsInfo->isChecked());
  const bool sum_per_sample = (m_sumDetsPerSample->isVisible() && m_sumDetsPerSample->isChecked());
  const bool sum_per_det = (m_sumSamplesPerDets->isVisible() && m_sumSamplesPerDets->isChecked());
  const bool fore_plus_back_files = (m_forePlusBack && m_forePlusBack->isVisible() && m_forePlusBack->isChecked());
  const bool use_disp_fore = (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked());
  const bool use_disp_back = (m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked());
  const bool use_disp_seco = (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked());
  
  
  if( !start_spec )
    throw runtime_error( "No file selected for export." );
  
  // First we'll check for all the cases where we want the whole file
  if( (start_spec->num_measurements() == 1) && !remove_gps )
    return start_spec;
  
  set<int> samples = currentlySelectedSamples();
  vector<string> detectors = currentlySelectedDetectors();
  
  // TODO: we should probably filter out unwanted detectors and samples up front, rather than
  //       waiting for the summing step
  
  const bool backgroundSub = (m_backSubFore && m_backSubFore->isVisible() && m_backSubFore->isChecked());
  const bool sumAll = ((m_sumAllToSingleRecord->isVisible() && m_sumAllToSingleRecord->isChecked())
                       || (max_records < 2));
  const bool foreToSingleRecord = (m_sumForeToSingleRecord->isVisible() && m_sumForeToSingleRecord->isChecked());
  const bool backToSingleRecord = (m_sumBackToSingleRecord->isVisible() && m_sumBackToSingleRecord->isChecked());
  const bool secoToSingleRecord = (m_sumSecoToSingleRecord->isVisible() && m_sumSecoToSingleRecord->isChecked());
  const bool filterDets = (m_filterDetector && m_filterDetector->isVisible() && m_filterDetector->isChecked());
  const bool sumDetectorsPerSample = (sum_per_sample && !sum_per_det && ((max_records >= samples.size()) && (detectors.size() > 1)));
  const bool sumSamplesPerDetector = (sum_per_det && !sum_per_sample && ((max_records >= detectors.size()) && (samples.size() > 1)));
  
  
  shared_ptr<SpecMeas> answer = make_shared<SpecMeas>();
  answer->uniqueCopyContents( *start_spec );

  answer->set_uuid( "" );
  
  if( remove_gps )
  {
    // TODO: the below causes mean lat/lon to be recalculated after each call should just add a SpecMeas::clear_gps_coordinates() function
    for( const auto &m : answer->measurements() )
      answer->set_position( -999.9, -999.9, SpecUtils::time_point_t{}, m );
  }//if( remove_gps )
  
  if( start_spec->num_measurements() == 1 )
    return answer;
  
  set<set<int>> samplesToSum;
  map<set<int>,SpecUtils::SourceType> sampleSourceTypes;
  if( sumDetectorsPerSample )
  {
    for( const int sample : samples )
      samplesToSum.insert( set<int>{sample} );
  }//if( sum detectors per sample )
  

  if( (foreToSingleRecord || (use_disp_fore && (max_records <= 2))) && !fore_plus_back_files )
  {
    assert( use_disp_fore );
    assert( !fore_plus_back_files );
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Foreground);
    samplesToSum.insert( samples );
    sampleSourceTypes[samples] = SpecUtils::SourceType::Foreground;
  }//if( foreground to single record )
  
  
  if( (backToSingleRecord || (use_disp_back && (max_records <= 2))) && !fore_plus_back_files )
  {
    assert( use_disp_back );
    assert( !fore_plus_back_files );
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
    samplesToSum.insert( samples );
    sampleSourceTypes[samples] = SpecUtils::SourceType::Background;
  }//if( background to single record )
  
  
  if( (secoToSingleRecord || (use_disp_seco && (max_records <= 2))) && !fore_plus_back_files )
  {
    assert( use_disp_seco );
    assert( !fore_plus_back_files );
    const set<int> &samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::SecondForeground);
    samplesToSum.insert( samples );
    sampleSourceTypes[samples] = SpecUtils::SourceType::Unknown;
  }//if( secondary to single record )
  
  
  if( fore_plus_back_files )
  {
    vector<SpecUtils::SourceType> src_types;
    if( max_records <= 2 )
    {
      src_types = vector<SpecUtils::SourceType>{ SpecUtils::SourceType::IntrinsicActivity,
        SpecUtils::SourceType::Calibration,
        SpecUtils::SourceType::Foreground,
        SpecUtils::SourceType::Background,
        SpecUtils::SourceType::Unknown
        };//src_types array definition
    }else
    {
      if( foreToSingleRecord )
        src_types.push_back( SpecUtils::SourceType::Foreground );
      if( backToSingleRecord )
        src_types.push_back( SpecUtils::SourceType::Background );
      if( secoToSingleRecord )
        src_types.push_back( SpecUtils::SourceType::Unknown );
    }//if( max_records <= 2 ) / else
    
    for( const SpecUtils::SourceType type : src_types )
    {
      set<int> samples;
      for( const auto &m : answer->measurements() )
      {
        if( m->source_type() == type )
          samples.insert( m->sample_number() );
      }
      
      // We dont expect Intrinsic or Calibration to ever be present when a foreground
      //  and background file has been combined into a single file
      //  See `currentlySelectedFile()` implementation.
      assert( (type != SpecUtils::SourceType::IntrinsicActivity)
             && (type != SpecUtils::SourceType::Calibration)
             || samples.empty() );
      
      if( !samples.empty() )
      {
        samplesToSum.insert( samples );
        sampleSourceTypes[samples] = type;
      }
    }//for( SpecUtils::SourceType type : src_types )
  }//if( fore_plus_back_files )
  
  set<set<int>> peaks_to_remove;
  map<set<int>,shared_ptr<const deque<shared_ptr<const PeakDef>>>> peaks_to_set;
  set<shared_ptr<const SpecUtils::Measurement>> meas_to_remove;
  vector<shared_ptr<SpecUtils::Measurement>> meas_to_add;
  
  if( backgroundSub )
  {
    // We'll check if the foreground and background have the same detectors
    //  for foreground and background, and if so, subtract on a detector by
    //  detector basis; if not we'll sum things, and do that.
    assert( m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked() );
    assert( (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked())
           || (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked()) );
    
    auto get_dets = [this,detectors,answer]( const SpecUtils::SpectrumType type ) -> set<string> {
      set<string> dets;
      
      for( int sample : m_interspec->displayedSamples(type) )
      {
        for( const string &det : detectors )
        {
          const auto m = answer->measurement( sample, det );
          if( m )
            dets.insert( det );
        }
      }//for( int sample : fore_samples )
      
      return dets;
    };//auto get_dets lamda
    
    
    const set<string> seco_dets = get_dets( SpecUtils::SpectrumType::SecondForeground );
    
    
    if( use_disp_fore || use_disp_seco )
    {
      auto make_subtracted = [&]( const SpecUtils::SpectrumType type ){
      
        assert( (type == SpecUtils::SpectrumType::Foreground)
               || (type == SpecUtils::SpectrumType::SecondForeground) );
        
        bool summed_det_by_det = false;
    
        const double fore_sf = m_interspec->displayScaleFactor( type );
        const double back_sf = m_interspec->displayScaleFactor( SpecUtils::SpectrumType::Background );
        
        // In the context of this lamda, we'll call either foreground or secondary spectra, foreground
        const set<string> fore_dets = get_dets( type );
        const set<string> back_dets = get_dets( SpecUtils::SpectrumType::Background );
        const vector<string> fore_dets_vec( begin(fore_dets), end(fore_dets) );
        const vector<string> back_dets_vec( begin(back_dets), end(back_dets) );
        
        const set<int> disp_fore_samples = m_interspec->displayedSamples(type);
        const set<int> disp_back_samples = m_interspec->displayedSamples(SpecUtils::SpectrumType::Background);
        
        if( disp_fore_samples.empty() || disp_back_samples.empty() )
          return;
        
        auto do_sub = [&]( const set<int> &fore_samples, const vector<string> &fore_dets,
                            const set<int> &back_samples, const vector<string> &back_dets ){
          auto foreground = answer->sum_measurements( fore_samples, fore_dets, nullptr );
          auto background = answer->sum_measurements( back_samples, back_dets, nullptr );
          assert( foreground && background );
          if( !foreground || !background )
            throw runtime_error( "Failed to sum foreground or background." );
          
          try
          {
            shared_ptr<const deque<shared_ptr<const PeakDef>>> orig_peaks = answer->peaks(fore_samples);
            shared_ptr<const vector<float>> fore_counts = foreground->gamma_counts();
            shared_ptr<const vector<float>> back_counts = background->gamma_counts();
            
            // Make sure back_counts has the same energy calibration and fore_counts, so we can subtract
            //  on a bin-by-bin basis
            if( background->energy_calibration() != foreground->energy_calibration()
               && (*background->energy_calibration()) != (*foreground->energy_calibration()) )
            {
              auto new_backchan = make_shared<vector<float>>( fore_counts->size(), 0.0f );
              SpecUtils::rebin_by_lower_edge( *background->channel_energies(), *back_counts,
                                             *foreground->channel_energies(), *new_backchan );
              back_counts = new_backchan;
            }
            
            // Create what will be the background subtracted foreground
            auto back_sub_counts = make_shared<vector<float>>( *fore_counts );
            
            //back_counts and fore_counts should always be the same size, but we'll be a bit verbose anyway
            assert( back_counts->size() == fore_counts->size() );
            const size_t nchann = std::min( back_counts->size(), fore_counts->size() );
            
            // Do the actual background subtraction
            for( size_t i = 0; i < nchann; ++i )
            {
              float &val = (*back_sub_counts)[i];
              val *= fore_sf;
              val -= back_sf*(*back_counts)[i];
            }//for( size_t i = 0; i < nchann; ++i )
            
            // Create a new Measurement object, based on the old foreground
            auto newspec = make_shared<SpecUtils::Measurement>( *foreground );
            newspec->set_gamma_counts( back_sub_counts, foreground->live_time(), foreground->real_time() );
            vector<string> remarks = foreground->remarks();
            remarks.push_back( "This spectrum has been background subtracted in InterSpec" );
            newspec->set_remarks( remarks );
            newspec->set_sample_number( *begin(fore_samples) );
            if( fore_dets.size() == 1 )
              newspec->set_detector_name( fore_dets.front() );
            else
              newspec->set_detector_name( "bkg_sub" );
            
            switch( type )
            {
              case SpecUtils::SpectrumType::Foreground:
                newspec->set_source_type( SpecUtils::SourceType::Foreground );
                break;
                
              case SpecUtils::SpectrumType::SecondForeground:
                newspec->set_source_type( SpecUtils::SourceType::Unknown );
                break;
                
              case SpecUtils::SpectrumType::Background:
                assert( 0 );
                break;
            }//switch( type )
            
            meas_to_add.push_back( newspec );
            
            for( const int sample : fore_samples )
            {
              for( const string &det : fore_dets )
              {
                const auto m = answer->measurement( sample, det );
                if( m )
                  meas_to_remove.insert( m );
              }
            }//for( const int sample : fore_samples )
            
            for( const int sample : back_samples )
            {
              for( const string &det : back_dets )
              {
                const auto m = answer->measurement( sample, det );
                if( m )
                  meas_to_remove.insert( m );
              }
            }//for( const int sample : back_samples )
            
            // TODO: refit peaks and then put those into the file - for the moment we'll just leave the original peaks
            if( orig_peaks && orig_peaks->size() )
            {
              peaks_to_set[set<int>{newspec->sample_number()}] = orig_peaks;
              peaks_to_remove.insert( fore_samples );
              peaks_to_remove.insert( back_samples );
            }
            
            /*
             // Re-fit peaks for now being background subtracted
             std::vector<PeakDef> refit_peaks;
             if( orig_peaks && orig_peaks->size() )
             {
             try
             {
             vector<PeakDef> input_peaks;
             for( const auto &i : *orig_peaks )
             input_peaks.push_back( *i );
             
             const double lowE = newspec->gamma_energy_min();
             const double upE = newspec->gamma_energy_max();
             
             refit_peaks = fitPeaksInRange( lowE, upE, 0.0, 0.0, 0.0, input_peaks, newspec, {}, true );
             
             std::deque<std::shared_ptr<const PeakDef> > peakdeque;
             for( const auto &p : refit_peaks )
             peakdeque.push_back( std::make_shared<const PeakDef>(p) );
             
             newmeas->setPeaks( peakdeque, {newspec->sample_number()} );
             }catch( std::exception &e )
             {
             summed_det_by_det = false;
             assert( 0 );
             break;
             }//try / catch to fit peaks
             }//if( we need to refit peaks )
             */
          }catch( std::exception &e )
          {
            cerr << "Error subtracting background from foreground on det-by-det basis: " << e.what() << endl;
            summed_det_by_det = false;
            assert( 0 );
          }//try / catch
         };//auto do_sub lamda
        
        
        if( fore_dets == back_dets )
        {
          summed_det_by_det = true;
          for( const string &det : fore_dets )
          {
            do_sub( disp_fore_samples, {det}, disp_back_samples, {det} );
          }//for( const string &det : fore_dets )
        }//if( fore_dets == back_dets )
        
        
        if( !summed_det_by_det )
        {
          //sum all foreground/background, and then subtract
          do_sub( disp_fore_samples, fore_dets_vec, disp_back_samples, back_dets_vec );
        }//if( !summed_det_by_det )
      };//make_subtracted lamda
      
      
      if( (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked()) )
        make_subtracted( SpecUtils::SpectrumType::Foreground );
      
      if( m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked() )
        make_subtracted( SpecUtils::SpectrumType::SecondForeground );
    }//if( create background subtracted foreground )
  }//if( background subtract )
  
  
  for( const set<int> &sum_samples : samplesToSum )
  {
    bool single_meas = false;
    
    if( sum_samples.size() == 1 )
    {
      // We'll try to just use a copy of the measurement, instead of the sum of a single measurement
      //  (I dont fully trust we copy all meta-info over perfectly when doing a sum).
      size_t ndet = 0;
      shared_ptr<const SpecUtils::Measurement> single_record;
      for( const string &det : detectors )
      {
        auto m = answer->measurement( *begin(sum_samples), det );
        if( m )
        {
          ndet += 1;
          if( ndet > 1 )
            break;
          
          single_record = m;
        }
      }//for( const string &det : detectors )
      
      if( ndet == 1 )
      {
        assert( single_record );
        single_meas = true;
        meas_to_add.push_back( make_shared<SpecUtils::Measurement>(*single_record) );
        
        const auto typePos = sampleSourceTypes.find(sum_samples);
        if( typePos != end(sampleSourceTypes) )
        {
          meas_to_add.back()->set_source_type( typePos->second );
          sampleSourceTypes[samples] = typePos->second;
        }else
        {
          sampleSourceTypes[samples] = SpecUtils::SourceType::Background;
        }
        
        meas_to_remove.insert( single_record );
      }
    }//if( sum_samples.size() == 1 )
    
    
    if( !single_meas )
    {
      assert( !sum_samples.empty() );
      
      shared_ptr<deque<shared_ptr<const PeakDef>>> peaks = answer->peaks(sum_samples);
      
      shared_ptr<SpecUtils::Measurement> m = answer->sum_measurements( sum_samples, detectors, nullptr );
      m->set_sample_number( *begin(sum_samples) );
      const auto typePos = sampleSourceTypes.find(sum_samples);
      if( typePos != end(sampleSourceTypes) )
        m->set_source_type( typePos->second );
      
      meas_to_add.push_back( m );
      
      answer->setPeaks( {}, sum_samples );
      if( peaks && peaks->size() )
        answer->setPeaks( *peaks, {m->sample_number()} );
      else
        answer->setPeaks( {}, {m->sample_number()} );
      
      // We will make a feeble attempt to preserve title, or for portal data
      //  set it to background, for that sample
      bool all_background = true;
      set<string> titles;
      for( const int sample : sum_samples )
      {
        for( const string &det : detectors )
        {
          auto orig = answer->measurement( sample, det );
          if( orig )
          {
            const bool empty_title = orig->title().empty();
            if( !empty_title )
              titles.insert( orig->title() );
            
            all_background &= (SpecUtils::icontains( orig->title(), "Background")
                               || (orig->source_type() == SpecUtils::SourceType::Background)
                               || (answer->passthrough() && (orig->occupied() == SpecUtils::OccupancyStatus::NotOccupied)));
            meas_to_remove.insert( orig );
          }
        }//for( const string &det : detectors )
      }//for( const int sample : sum_samples )
      
      if( titles.size() == 1 )
        m->set_title( *begin(titles) );
      else if( all_background )
        m->set_title( "Background" );
      else
        m->set_title( "" );
    }//if( !single_meas )
  }//for( const set<int> &samples : samplesToSum )
  

  for( const auto &m : answer->measurements() )
  {
    assert( m );
    if( !m )
      continue;
    const bool is_wanted_sample = samples.count( m->sample_number() );
    const auto det_pos = std::find( begin(detectors), end(detectors), m->detector_name() );
    const bool is_wanted_det = (det_pos != end(detectors));
    if( !is_wanted_sample || !is_wanted_det )
      meas_to_remove.insert( m );
  }//for( const auto &m : answer->measurements() )
  
  if( sumSamplesPerDetector )
  {
    assert( meas_to_add.empty() );
    map<string,vector<shared_ptr<const SpecUtils::Measurement>>> det_to_meas;
    for( const auto &m : answer->measurements() )
    {
      if( m && !meas_to_remove.count(m) )
        det_to_meas[m->detector_name()].push_back( m );
    }//for( const auto &m : answer->measurements() )
    
    for( const auto &det_meas : det_to_meas )
    {
      const string &det = det_meas.first;
      const vector<shared_ptr<const SpecUtils::Measurement>> &meass = det_meas.second;
      set<int> samples;
      for( const auto &m : meass )
        samples.insert( m->sample_number() );
      
      if( samples.size() > 1 )
      {
        auto sum_meas = answer->sum_measurements( samples, {det}, nullptr );
        sum_meas->set_sample_number( *begin(samples) );
        sum_meas->set_detector_name( det ); //just to make sure
        
        meas_to_add.push_back( sum_meas );
        for( const auto &m : meass )
          meas_to_remove.insert(m);
      }//if( samples.size() > 1 )
    }//for( const auto &det_meas : det_to_meas )
  }//if( sumSamplesPerDetector )
  
  
  if( !meas_to_remove.empty() )
  {
    vector<shared_ptr<const SpecUtils::Measurement>> to_remove( begin(meas_to_remove), end(meas_to_remove) );
    answer->remove_measurements( to_remove ); //This calls `cleanup_after_load(...)`
  }
  
  if( !meas_to_add.empty() )
  {
    for( const shared_ptr<SpecUtils::Measurement> &i : meas_to_add )
      answer->add_measurement( i, false );
  }
  
  if( !meas_to_remove.empty() || !meas_to_add.empty() )
  {
    answer->set_uuid( "" );
    answer->cleanup_after_load( SpecUtils::SpecFile::DontChangeOrReorderSamples );
    
    // `detectors` and `samples` are maybe no longer valid, since we may have summed them together
    //  - we will fix this up, although it probably makes more sense to just filter out unwanted
    //  detectors and samples, up front (but this maybe creates other problems?)
    const vector<string> &dets_now = answer->detector_names();
    const set<int> &samples_now = answer->sample_numbers();
    vector<string> names_to_remove;
    for( const string &d : detectors )
    {
      if( std::find(begin(dets_now), end(dets_now), d) == end(dets_now) )
        names_to_remove.push_back( d );
    }//for( const string &d : detectors )
    
    for( const string &d : names_to_remove )
    {
      assert( std::find( begin(detectors), end(detectors), d ) != end(detectors) );
      detectors.erase( std::find( begin(detectors), end(detectors), d ) );
    }
    
    set<int> samples_to_remove;
    for( const int sample : samples )
    {
      if( !samples_now.count(sample) )
        samples_to_remove.insert(sample);
    }
    for( const int sample : samples_to_remove )
      samples.erase( sample );
    
    // Now lets map sample numbers to 1 through N
    samples.clear();
    int new_sample_number = 0;
    vector<pair<int,int>> old_to_new_samplenum;
    for( const int old_sample_number : answer->sample_numbers() )
    {
      ++new_sample_number;
      old_to_new_samplenum.emplace_back( old_sample_number, new_sample_number );
      samples.insert( new_sample_number );
    }//for( loop over sample numbers )
    
    answer->change_sample_numbers( old_to_new_samplenum );
  }//if( !meas_to_remove.empty() || !meas_to_add.empty() )
  
  
  if( sumAll )
  {
    const vector<shared_ptr<const SpecUtils::Measurement>> orig_meass = answer->measurements();
    
    // Next call throws exception if invalid sample number, detector name, or cant find energy
    //  binning to use.  And returns nullptr if empty sample numbers or detector names.
    shared_ptr<SpecUtils::Measurement> sum_meas = answer->sum_measurements( samples, detectors, nullptr );
    assert( sum_meas );
    if( !sum_meas )
      throw runtime_error( "Error summing records - perhaps empty sample numbers or detector names." );
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> peaks = answer->peaks(samples);
    answer->setPeaks( {}, samples );
    
    sum_meas->set_sample_number( 1 );
    if( peaks && peaks->size() )
      answer->setPeaks( *peaks, {1} );
    else
      answer->setPeaks( {}, {1} );
    
    answer->remove_measurements( orig_meass );
    answer->add_measurement( sum_meas, true );
  }//if( sumAll )
  
  /*
   // I think we have already summed the detectors, for each sample, together
   //  so we shouldnt do it again
   else if( sumDetectorsPerSample && (answer->detector_names().size() > 1) )
  {
    set<int> problem_samples;
    map<int,vector<shared_ptr<const SpecUtils::Measurement>>> sample_to_meas;
    for( const auto &m : answer->measurements() )
      sample_to_meas[m->sample_number()].push_back( m );
    
    vector<shared_ptr<SpecUtils::Measurement>> meas_to_add;
    vector<shared_ptr<const SpecUtils::Measurement>> meas_to_remove;
    
    for( const auto &sample_meass : sample_to_meas )
    {
      const int sample_num = sample_meass.first;
      const auto &meass = sample_meass.second;
      if( meass.size() < 2 )
        continue;
      
      try
      {
        // `detectors` is no longer valid, since we may have summed them together
        const vector<string> det_names = answer->detector_names();
        shared_ptr<SpecUtils::Measurement> m = answer->sum_measurements( {sample_num}, det_names, nullptr );
        if( !m )
          throw runtime_error( "No gamma spectra found" );
        
        m->set_sample_number( sample_num );
        meas_to_add.push_back( m );
      }catch( exception &e )
      {
        problem_samples.insert( sample_num );
        cerr << "Failed to sum " << meass.size() << " measurements for sample " << sample_num << " to create an output file." << endl;
      }//try / catch
      
      meas_to_remove.insert(end(meas_to_remove), begin(meass), end(meass) );
    }//for( const auto &sample_meass : sample_to_meas )
    
    if( meas_to_add.empty() )
      throw runtime_error( "Error summing individual detectors together to create single spectrum." );
    assert( !meas_to_remove.empty() );
    
    answer->remove_measurements( meas_to_remove );
    
    for( const auto &m : meas_to_add )
      answer->add_measurement( m, false );
    answer->cleanup_after_load();
    
    if( problem_samples.size() )
    {
      const string sample_str = SpecUtils::sequencesToBriefString( problem_samples );
      
      string msg = "There was an issue summing detectors together for ";
      if( sample_str.size() > 20 )
        msg += std::to_string(problem_samples.size()) + " sample numbers.";
      else
        msg += "samples " + sample_str + ".";
      msg += "<br />A possible cause is if the sample number didn't contain any gamma spectra."
             "<br />These samples are not included in the output.";
      passMessage( msg, WarningWidget::WarningMsgHigh );
    }//if( problem_samples.size() )
  }//if( sum to single record ) / else if( sumDetectorsPerSample )
  */
  
  // We will check for this later as well, but we'll remove as much of the InterSpec info
  //  here as well (the displayed sample numbers and wont be removed though).
  if( m_excludeInterSpecInfo->isVisible() && m_excludeInterSpecInfo->isChecked() )
  {
    answer->removeAllPeaks();
    answer->setShieldingSourceModel( std::unique_ptr<rapidxml::xml_document<char>>{} );
#if( USE_REL_ACT_TOOL )
    answer->setRelActManualGuiState( std::unique_ptr<rapidxml::xml_document<char>>{} );
#endif
    answer->setDetector( nullptr );
  }//if( get rid of InterSpec info )
  
  
  return answer;
}//std::shared_ptr<const SpecMeas> generateFileToSave()


#if( USE_QR_CODES )
void ExportSpecFileTool::handleGenerateQrCode()
{
  shared_ptr<const SpecMeas> spec;
  
  try
  {
    spec = ExportSpecFileTool::generateFileToSave();
  }catch( std::exception &e )
  {
    passMessage( "Sorry, there was an issue preparing the file for saving: " + string(e.what()),
                WarningWidget::WarningMsgHigh );
  }//try / catch
  
  assert( spec );
  if( !spec )
    return;
  
  assert( (spec->num_measurements() == 1) || (spec->num_measurements() == 2) );
  if( (spec->num_measurements() != 1) && (spec->num_measurements() != 2) )
  {
    passMessage( "Sorry, there was an internal logic error preparing file for saving.",
                WarningWidget::WarningMsgHigh );
    return;
  }//if( something bad that shouldnt happen )
  
 
  try
  {
    string model;
    if( spec->detector_type() != SpecUtils::DetectorType::Unknown )
      model = detectorTypeToString( spec->detector_type() );
    if( model.empty() )
      model = spec->instrument_model();
  
    const vector<SpecUtils::UrlSpectrum> urlspec = SpecUtils::to_url_spectra( spec->measurements(), model );
    
    auto successfullyDone = wApp->bind( boost::bind( &ExportSpecFileTool::emitDone, this, true ) );
    
    displayQrCode( urlspec, spec, successfullyDone, false );
  }catch( std::exception &e )
  {
    auto dialog = new SimpleDialog( "Error", "Failed to encoded spectrum to a URI: " + string(e.what()) );
    dialog->addButton( "Ok" );
  }//try catch
}//void handleGenerateQrCode()
#endif


void ExportSpecFileTool::scheduleAddingUndoRedo()
{
  // We will actually just mark that the widget will should be rendered,
  //  so this way we will only try to add a undo/redo once per event loop.
  //  (and we will only schedule this if we are not in a undo/redo step)
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
    scheduleRender();
}//void scheduleAddingUndoRedo()


void ExportSpecFileTool::updateUndoRedo()
{
  try
  {
    string appurl = encodeStateToUrl();
    if( m_last_state_uri && ((*m_last_state_uri) == appurl) )
      return;
    
    shared_ptr<const string> curr_state = make_shared<string>( std::move(appurl) );
    shared_ptr<const string> prev_state = m_last_state_uri;
    
    m_last_state_uri = curr_state;
    
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    
    // We shouldnt be in a undo/redo state, if we are here, but we'll check, jic
    if( !prev_state || !undoRedo || undoRedo->isInUndoOrRedo() )
      return;
    
    //cout << "ExportSpecFileTool::updateUndoRedo(): am adding undo/redo step" << endl;
    
    auto undo_redo = [prev_state,curr_state]( const bool undo ){
      try
      {
        InterSpec *interspec = InterSpec::instance();
        ExportSpecFileWindow *export_window = interspec ? interspec->createExportSpectrumFileDialog() : nullptr;
        shared_ptr<const string> state_str = undo ? prev_state : curr_state;
        if( !export_window || !state_str || state_str->empty() )
          throw runtime_error( "No export dialog, or no state avaialable" );
        export_window->handleAppUrl( *prev_state );
      }catch( std::exception &e )
      {
        string msg = "Error executing undo/redo step for Spectrum File Export tool: " + string(e.what());
        passMessage( msg, WarningWidget::WarningMsgHigh );
      }
    };//undo_redo lambda
    
    auto undo = [undo_redo](){ undo_redo(true); };
    auto redo = [undo_redo](){ undo_redo(false); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change spectrum file export tool." );
  }catch( std::exception &e )
  {
    const string msg = "Error adding undo/redo step for Spectrum File Export tool: " + string(e.what());
    log("error") << msg;
    passMessage( msg, WarningWidget::WarningMsgHigh );
  }//try / catch
}//void updateUndoRedo()

  
void ExportSpecFileTool::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  WContainerWidget::render( flags );
  
  // Add the undo/redo step, not from the current event loop, so we dont cause delays
  //  Also, we'll check, JIC, to make sure we arent in an undo/redo step execution,
  //   but this is likely not the case anyway.
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->isInUndoOrRedo() )
  {
    WServer *server = WServer::instance();
    auto worker = wApp->bind( boost::bind(&ExportSpecFileTool::updateUndoRedo, this) );
    server->schedule( 25, wApp->sessionId(), worker );
  }
}//void render( Wt::WFlags<Wt::RenderFlag> flags )

              
void ExportSpecFileTool::handleAppUrl( std::string query_str )
{
  //Example: query_str="V=1&Fore=1&Format=CSV&Samples=1-5&SumAllToSingleRecord=1"
  
  //cout << "query_str=\"" << query_str << "\"" << endl;
  
  map<string,string> parts = AppUtils::query_str_key_values( query_str );
  if( !parts.count("V") || (parts["V"] != "1") )
    throw runtime_error( "fromAppUrl: missing or invalid 'V'" );
  
  
  auto find_spec = [this]( SpecUtils::SpectrumType type ) -> int {
    auto spec = m_interspec->measurment(type);
    if( !spec )
      return -1;
    
    const SpecMeasManager * const measManager = m_interspec->fileManager();
    const SpectraFileModel * const fileModel = measManager ? measManager->model() : nullptr;
    if( !fileModel )
      return -1;
    
    for( int index = 0; index < fileModel->rowCount(); ++index )
    {
      const shared_ptr<const SpectraFileHeader> header = fileModel->fileHeader(index);
      const shared_ptr<const SpecMeas> this_spec = header ? header->parseFile() : nullptr;
      if( this_spec == spec )
        return index;
    }
    
    assert( 0 );
    return -1;
  };
  
  
  if( parts.count("FOREPLUSBACK") )
  {
    if( m_forePlusBack )
    {
      m_forePlusBack->setChecked( true );
      handleForePlusBackChanged();
    }
    
  }else
  {
    if( m_forePlusBack )
    {
      m_forePlusBack->setChecked( false );
      handleForePlusBackChanged();
    }
    
    if( m_fileSelect && m_fileSelect->isHidden() )
      m_fileSelect->show();
    
    int index = -1;
    if( parts.count("FORE") )
      index = find_spec( SpecUtils::SpectrumType::Foreground );
    else if( parts.count("BACK") )
      index = find_spec( SpecUtils::SpectrumType::Background );
    else if( parts.count("SECO") )
      index = find_spec( SpecUtils::SpectrumType::SecondForeground );
    else if( parts.count("FILEUUID") )
    {
      const string wanted_uuid = parts["FILEUUID"];
      
      const SpecMeasManager * const measManager = m_interspec->fileManager();
      const SpectraFileModel * const fileModel = measManager ? measManager->model() : nullptr;
      
      if( fileModel )
      {
        for( int row = 0; row < fileModel->rowCount(); ++row )
        {
          const shared_ptr<const SpectraFileHeader> header = fileModel->fileHeader(row);
          const shared_ptr<const SpecMeas> this_spec = header ? header->parseFile() : nullptr;
          if( this_spec )
          {
            const string this_uuid = clean_uuid( this_spec->uuid() );
            if( SpecUtils::iequals_ascii(this_uuid, wanted_uuid) )
            {
              index = row;
              break;
            }
          }//if( this_spec )
        }//for( int row = 0; row < fileModel->rowCount(); ++row )
      }//if( fileModel )
    }else if( parts.count("FILEINDEX") )
    {
      string str = parts["FILEINDEX"];
      int val;
      if( SpecUtils::parse_int(str.c_str(), str.size(), val) )
        index = val;
    }
    
    if( m_fileSelect && (index < m_fileSelect->count()) )
      m_fileSelect->setCurrentIndex( index );
  }//if( parts.count("FOREPLUSBACK") ) / else
  
  
  handleFileSelectionChanged();
  
  
  SpecUtils::SaveSpectrumAsType file_type = SpecUtils::SaveSpectrumAsType::N42_2012;
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  
  if( parts.count("FORMAT") )
  {
    for( size_t i = 0;  i < sm_file_type_to_uri_name.size(); ++i )
    {
      const string &file_type_str = sm_file_type_to_uri_name[i].second;
      
      if( SpecUtils::iequals_ascii( parts["FORMAT"], file_type_str ) )
      {
        file_type = sm_file_type_to_uri_name[i].first;
        break;
      }
    }//for( size_t i = 0;  i < sm_file_type_to_uri_name.size(); ++i )

    for( WMenuItem *item : m_formatMenu->items() )
    {
      const void *data = item->data();
      const SpecUtils::SaveSpectrumAsType t = reinterpret_cast<SpecUtils::SaveSpectrumAsType&>( data );
      if( t == file_type )
      {
        m_formatMenu->select( item );
        break;
      }
    }
  }//if( parts.count("FORMAT") )
  
  handleFormatChange();
  
  const bool all = parts.count("ALLSAMPLES");
  if( m_dispForeSamples )
    m_dispForeSamples->setChecked( !all && parts.count("DISPFORE") );
  if( m_dispBackSamples )
    m_dispBackSamples->setChecked( !all && parts.count("DISPBACK") );
  if( m_dispSecondSamples )
    m_dispSecondSamples->setChecked( !all && parts.count("DISPSECO") );
  if( m_allSamples )
    m_allSamples->setChecked( all );

  if( m_customSamples && parts.count("SAMPLES") )
  {
    m_customSamples->setChecked( true );
    const string txt = parts["SAMPLES"];
    const pair<set<int>,string> samplenums = sampleNumbersFromTxtRange( txt, spec, true );
    m_customSamplesEdit->setText( samplenums.second );
  }else if( m_customSamples )
  {
    m_customSamples->setChecked( false );
    m_customSamplesEdit->hide();
    m_customSamplesEdit->setText( "" );
  }

  handleSamplesChanged();
  
  vector<string> dets;
  if( m_is_specific_file )
    dets = m_specific_detectors;
  else if( m_filterDetector && m_filterDetector->isVisible() && m_filterDetector->isChecked() )
    dets = currentlySelectedDetectors();
  
  if( parts.count("DETECTORS") )
  {
    if( m_filterDetector )
      m_filterDetector->setChecked( true );
    
    handleFilterDetectorCbChanged();
    
    vector<string> dets_in_url;
    SpecUtils::split( dets_in_url, parts["DETECTORS"], "-" );
    
    map<string,pair<string,bool>> label_to_orig;
    if( spec )
    {
      for( const string &name : spec->detector_names() )
      {
        const string clean_name = clean_uuid( name );
        const bool use = (std::find(begin(dets_in_url), end(dets_in_url), clean_name) != end(dets_in_url));
        
        // TODO: use Wt::Utils::htmlEncode, if that is what happens in WCheckBox, instead of creating a WCheckBox
        WCheckBox cb( name );
        label_to_orig[cb.text().toUTF8()] = {name, use};
      }
    }//if( spec )
    
    vector<string> answer;
    for( const auto w : m_detectorFilterCbs->children() )
    {
      WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
      if( !cb )
        continue;
      
      const string label_txt = cb->text().toUTF8();
      const auto pos = label_to_orig.find( label_txt );
      assert( pos != end(label_to_orig) );
      if( pos != end(label_to_orig) )
        cb->setChecked( pos->second.second );
      else
        cb->setChecked( false );
    }//for( const auto w : m_detectorFilterCbs->children() )
  }else
  {
    if( m_filterDetector )
      m_filterDetector->setChecked( false );
    
    handleFilterDetectorCbChanged();
    
    if( m_detectorFilterCbs )
    {
      for( const auto w : m_detectorFilterCbs->children() )
      {
        WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
        if( cb )
          cb->setChecked( true );
      }//for( const auto w : m_detectorFilterCbs->children() )
    }//if( m_detectorFilterCbs )
  }//if( parts.count("DETECTORS") ) / else
  
  if( m_sumAllToSingleRecord )
    m_sumAllToSingleRecord->setChecked( parts.count("SUMALLTOSINGLERECORD") );
  handleSumToSingleRecordChanged();
  
  if( m_sumForeToSingleRecord )
    m_sumForeToSingleRecord->setChecked( parts.count("SUMFORETOSINGLERECORD") );
  
  if( m_sumBackToSingleRecord )
    m_sumBackToSingleRecord->setChecked( parts.count("SUMBACKTOSINGLERECORD") );
  
  if( m_sumSecoToSingleRecord )
    m_sumSecoToSingleRecord->setChecked( parts.count("SUMSECOTOSINGLERECORD") );

  handleSumTypeToSingleRecordChanged();
  
  if( m_backSubFore )
    m_backSubFore->setChecked( parts.count("BACKSUBFORE") );

  handleBackSubForeChanged();
  
  if( m_sumDetsPerSample )
    m_sumDetsPerSample->setChecked( parts.count("SUMDETSPERSAMPLE") );
  
  if( m_sumSamplesPerDets )
    m_sumSamplesPerDets->setChecked( parts.count("SUMSAMPLEPERDET") );
  
  handleSumDetPerSampleChanged();
  
  if( m_excludeInterSpecInfo )
    m_excludeInterSpecInfo->setChecked( parts.count("NOINTERSPECINFO") );
  
  if( m_excludeGpsInfo )
    m_excludeGpsInfo->setChecked( parts.count("NOGPS") );
}//void handleAppUrl( std::string query_str )


std::string ExportSpecFileTool::encodeStateToUrl() const
{
  // Returns something like "V=1&FORMAT=N42-2012&Samples=1,2&..", and it will not be url-encoded.
  
  string answer = "V=1";
  
  shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  
  if( m_forePlusBack && m_forePlusBack->isVisible() && m_forePlusBack->isChecked() )
  {
    answer += "&ForePlusBack=1";
  }else if( spec )
  {
    if( spec == m_interspec->measurment(SpecUtils::SpectrumType::Foreground) )
      answer += "&Fore=1";
    else if( spec == m_interspec->measurment(SpecUtils::SpectrumType::Background) )
      answer += "&Back=1";
    else if( spec == m_interspec->measurment(SpecUtils::SpectrumType::SecondForeground) )
      answer += "&Seco=1";
    else
    {
      string uuid = clean_uuid( spec->uuid() );
      if( !uuid.empty() )
        answer += "&FileUuid=" + uuid;
      else
        answer += "&FileIndex=" + std::to_string( m_fileSelect->currentIndex() );
    }
  }else
  {
    answer += "&NoFile=1";
  }//if( fore + back ) / else if( spec ) / else
  
  
  const SpecUtils::SaveSpectrumAsType file_type = currentSaveType();
  
  string file_type_str;
  for( size_t i = 0; file_type_str.empty() && i < sm_file_type_to_uri_name.size(); ++i )
  {
    if( sm_file_type_to_uri_name[i].first == file_type )
      file_type_str = sm_file_type_to_uri_name[i].second;
  }
  
  answer += "&Format=" + file_type_str;
  
  if( m_dispForeSamples && m_dispForeSamples->isChecked() )
    answer += "&DispFore=1";
  if( m_dispBackSamples && m_dispBackSamples->isChecked() )
    answer += "&DispBack=1";
  if( m_dispSecondSamples && m_dispSecondSamples->isChecked() )
    answer += "&DispSeco=1";
  
  if( m_allSamples && m_allSamples->isChecked() )
    answer += "&AllSamples=1";
  
  if( m_is_specific_file )
  {
    //m_specific_samples
    //m_specific_detectors
  }else //if( m_is_specific_file )
  {
    if( m_customSamples && m_customSamples->isChecked() )
    {
      const string txt = m_customSamplesEdit->text().toUTF8();
      const pair<set<int>,string> samplenums = sampleNumbersFromTxtRange( txt, spec, true );
      answer += "&Samples=" + samplenums.second;
    }
  }//if( m_is_specific_file ) / else
  
  vector<string> dets;
  if( m_is_specific_file )
    dets = m_specific_detectors;
  else if( m_filterDetector && m_filterDetector->isChecked() )
    dets = currentlySelectedDetectors();
  
  if( dets.size() )
  {
    vector<string> all_dets = spec->detector_names();
    std::sort( begin(all_dets), end(all_dets) );
    std::sort( begin(dets), end(dets) );
    if( all_dets != dets )
    {
      answer += "&Detectors=";
      
      for( size_t i = 0; i < dets.size(); ++i )
      {
        if( i )
          answer += "-";
        answer += clean_uuid( dets[i] );
      }
    }//if( all_dets != dets )
  }//if( dets.size() )
  
  
  if( m_sumAllToSingleRecord && m_sumAllToSingleRecord->isChecked() )
    answer += "&SumAllToSingleRecord=1";
  
  if( m_sumForeToSingleRecord && m_sumForeToSingleRecord->isChecked() )
    answer += "&SumForeToSingleRecord=1";
  
  if( m_sumBackToSingleRecord && m_sumBackToSingleRecord->isChecked() )
    answer += "&SumBackToSingleRecord=1";
  
  if( m_sumSecoToSingleRecord && m_sumSecoToSingleRecord->isChecked() )
    answer += "&SumSecoToSingleRecord=1";
  
  if( m_backSubFore && m_backSubFore->isChecked() )
    answer += "&BackSubFore=1";
  
  if( m_sumDetsPerSample && m_sumDetsPerSample->isChecked() )
    answer += "&SumDetsPerSample=1";
  
  if( m_sumSamplesPerDets && m_sumSamplesPerDets->isChecked() )
    answer += "&SumSamplePerDet=1";
  
  if( m_excludeInterSpecInfo && m_excludeInterSpecInfo->isChecked() )
    answer += "&NoInterSpecInfo=1";
  
  if( m_excludeGpsInfo && m_excludeGpsInfo->isChecked() )
    answer += "&NoGps=1";
  
  //cout << "URL Query string='" << answer << "'" << endl;
  
  return answer;
}//std::string encodeStateToUrl() const


ExportSpecFileWindow::ExportSpecFileWindow( InterSpec *viewer )
  : SimpleDialog( ((viewer && !viewer->isPhone()) ? "Spectrum File Export" : ""), "" ),
  m_tool( nullptr )
{
  // If the CSS isnt available to set widths and stuff when HTML is first formed
  //  then the window wont be sized properly, so we have actually already loaded
  //  ExportSpecFile.css in InterSpecApp
  wApp->useStyleSheet( "InterSpec_resources/ExportSpecFile.css" );
  
  addStyleClass( "export-spec-file" );
  
  const int w = viewer->renderedWidth();
  const bool isPhone = (w > 100) ? (w < ns_min_screen_width_for_wide_layout) : (viewer && viewer->isPhone());

  if( isPhone )
  {
    resize( 320, Wt::WLength::Auto ); //320 px is about the smallest width Android phone to expect
  }else 
  {
    //setMinimumSize( WLength(w > 100 ? std::min(0.95*w, 800.0) : 800.0 ,WLength::Pixel), WLength::Auto );
    setMinimumSize( WLength(w > 100 ? std::min(0.95*w, 650.0) : 650.0, WLength::Pixel), WLength::Auto );
  }
  m_tool = new ExportSpecFileTool( viewer, contents() );
  m_tool->done().connect( boost::bind(&ExportSpecFileWindow::accept, this) );
  
  if( viewer && viewer->isPhone() )
  {
    addStyleClass( "export-spec-phone" );
#if( IOS )
    addStyleClass( "export-spec-iphone" );
#endif
  }
  
  rejectWhenEscapePressed();
}//ExportSpecFileWindow( constructor )


void ExportSpecFileWindow::setSpecificSpectrum( const std::shared_ptr<const SpecMeas> &spectrum,
                         const std::set<int> &samples,
                         const std::vector<std::string> &detectors,
                         InterSpec *viewer )
{
  delete m_tool;
  m_tool = new ExportSpecFileTool( spectrum, samples, detectors, viewer, contents() );
  m_tool->done().connect( boost::bind(&ExportSpecFileWindow::accept, this) );
}//void setSpecificSpectrum(...)


void ExportSpecFileWindow::handleAppUrl( const std::string &query_str )
{
  m_tool->handleAppUrl( query_str );
}//void handleAppUrl( std::string query_str )


std::string ExportSpecFileWindow::encodeStateToUrl() const
{
  return m_tool->encodeStateToUrl();
}//std::string encodeStateToUrl() const


