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
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WMenuItem>
#include <Wt/WResource>
#include <Wt/WPushButton>
#include <Wt/Http/Request>
#include <Wt/WApplication>
#include <Wt/Http/Response>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WMessageResourceBundle>

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ExportSpecFile.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"


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
}//namespace

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
                            const InterSpec *viewer );
    
    //handle_resource_request(): does the actual streaming of the SpecMeas.
    //  'samplenums' and 'detectornums' are only used for file formats where SpecMeas must be
    //  collapsed down into a single spectrum (i.e., Chn, IntegerSpcType, SpcBinaryFloat, SpcAscii,
    //  SpeIaea, Cnf, Tka, will have all the spectra summed into a single spectrum); for other formats
    //  the entire SpecMeas is written out.
    //
    //  Specifying empty detector and sample numbers will default to all samples/detectors.
    static void handle_resource_request(
                                  SpecUtils::SaveSpectrumAsType type,
                                  std::shared_ptr<const SpecMeas> Measurement,
                                  const std::set<int> &samplenums,
                                  const std::vector<std::string> &detectornums,
                                  const InterSpec *viewer,
                                  const Wt::Http::Request& request,
                                  Wt::Http::Response& response );
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
    
    if( !lock )
    {
      log("error") << "Failed to WApplication::UpdateLock in DownloadSpectrumResource.";
      response.out() << "Error grabbing application lock to form DownloadSpectrumResource resource; please report to InterSpec@sandia.gov.";
      response.setStatus(500);
      assert( 0 );
      return;
    }//if( !lock )
    
    assert( m_tool );
    if( !m_tool )
      return;

    /*
    const WModelIndexSet selected = m_display->treeView()->selectedIndexes();

    if( selected.empty() )
      return;
    
    vector<shared_ptr<const SpectraFileHeader> > headers;


    for( const WModelIndex &index : selected )
    {
      if( m_display->model()->level(index) != SpectraFileModel::FileHeaderLevel )
        continue;
      shared_ptr<const SpectraFileHeader> header;
      header = m_display->model()->fileHeader( index.row() );

      if( header )
        headers.push_back( header );
    }//for( const WModelIndex &index : selected )

    const set<int> samplenums;
    const vector<string> detectornames;
    
    try
    {
      shared_ptr<const SpecMeas> measurement;
      if( headers.size() > 1 || headers.empty() )
      {
        measurement = m_display->selectedToSpecMeas();
        handle_resource_request( m_type, measurement, samplenums,
                                detectornames, m_display->viewer(), request, response );
      }else if( headers.size() == 1 )
      {
        shared_ptr<const SpectraFileHeader> header = headers[0];
        measurement = header->parseFile();
        if( measurement )
          handle_resource_request( m_type, measurement, samplenums,
                                  detectornames, m_display->viewer(), request, response );
      }
    }catch( std::exception &e )
    {
      stringstream msg;
      msg << "Failed in creating resource to export: " << e.what();
      passMessage( msg.str(), WarningWidget::WarningMsgHigh );
      
  #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, msg.str().c_str() );
  #endif
    }
     */
    
    suggestFileName( "SomeDumbName.txt" );
    response.out() << "Hello dude.";
    response.setMimeType( "text/plain" );
    
    m_download_finished.emit();
  }//void DownloadSpectrumResource::handleRequest(...)


  void DownloadSpectrumResource::write_file( std::ostream &output,
                         const SpecUtils::SaveSpectrumAsType type,
                         const std::shared_ptr<const SpecMeas> measurement,
                         const std::set<int> &samplenums,
                         const std::vector<std::string> &detectornames,
                         const InterSpec *viewer )
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
        measurement->write_2006_N42( output );
        break;
        
      case SpecUtils::SaveSpectrumAsType::N42_2012:
        measurement->write_2012_N42( output );
        break;
        
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
        
      case SpecUtils::SaveSpectrumAsType::NumTypes:
        break;
    }//switch( type )
  }//DownloadSpectrumResource::write_file(...)


  void DownloadSpectrumResource::handle_resource_request(
                                SpecUtils::SaveSpectrumAsType type,
                                std::shared_ptr<const SpecMeas> measurement,
                                const std::set<int> &samplenums,
                                const std::vector<std::string> &detectornames,
                                const InterSpec *viewer,
                                const Wt::Http::Request& /*request*/,
                                Wt::Http::Response& response )
  {
    switch( type )
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
        
      case SpecUtils::SaveSpectrumAsType::NumTypes:
      break;
    }//switch( type )

    if( !measurement )
      return;
    
    DownloadSpectrumResource::write_file( response.out(), type, measurement, samplenums, detectornames, viewer );
  }//void handle_resource_request(...)
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
  m_includeInterSpecInfo( nullptr ),
  m_msg( nullptr ),
  m_export_btn( nullptr ),
  m_resource( nullptr )
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
  m_includeInterSpecInfo( nullptr ),
  m_msg( nullptr ),
  m_export_btn( nullptr ),
  m_resource( nullptr )
{
  init();
}//
  

void ExportSpecFileTool::init()
{
  if( !m_interspec )
    throw logic_error( "ExportSpecFileTool: Invalid InterSpec pointer." );
  
  addStyleClass( "ExportSpecFileTool" );
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_interspec );
  
  WContainerWidget *body = new WContainerWidget( this );
  body->addStyleClass( "ExportSpecFileBody" );
  
  
  WContainerWidget *fileSelectDiv = new WContainerWidget( body );
  fileSelectDiv->addStyleClass( "ExportSpecSelect" );
  
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
  
  WText *title = new WText( "File Format", menuHolder );
  title->addStyleClass( "ExportColTitle" );
  
  m_formatMenu = new WMenu( menuHolder );
  m_formatMenu->itemSelected().connect( this, &ExportSpecFileTool::handleFormatChange );
  m_formatMenu->addStyleClass( "SideMenu VerticalNavMenu LightNavMenu ExportSpecFormatMenu" );
  
  
  
  
  const string docroot = wApp->docRoot();
  const string bundle_file = SpecUtils::append_path(docroot, "InterSpec_resources/static_text/spectrum_file_format_descriptions" );
  Wt::WMessageResourceBundle descrip_bundle;
  descrip_bundle.use(bundle_file,true);
  
  auto addFormatItem = [this, &descrip_bundle]( const char *label, SpecUtils::SaveSpectrumAsType type ){
    WMenuItem *item = m_formatMenu->addItem( label );
    item->clicked().connect( boost::bind(&right_select_item, m_formatMenu, item) );
    item->setData( reinterpret_cast<void *>(type) );
    
    string description;
    if( descrip_bundle.resolveKey(label, description) )
    {
      SpecUtils::trim( description );
      
      WImage *img = new WImage( item );
      img->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
      img->setStyleClass("Wt-icon");
      img->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
      img->setFloatSide( Wt::Side::Right );
      
      HelpSystem::attachToolTipOn( img, description, true,
                                  HelpSystem::ToolTipPosition::Right,
                                  HelpSystem::ToolTipPrefOverride::AlwaysShow );
    }//if( we have the description of the file )
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
  addFormatItem( "QR-code/URL", SpecUtils::SaveSpectrumAsType::NumTypes );
  
  // Meas/samples to include
  m_samplesHolder = new WContainerWidget( body );
  m_samplesHolder->addStyleClass( "ExportSpecSamples" );
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

  m_filterDetector = new WCheckBox( "Filter Dets.", m_samplesHolder );
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
  m_sumAllToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumToSingleRecordChanged );
  m_sumAllToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumToSingleRecordChanged );
  
  m_sumForeToSingleRecord = new WCheckBox( "Sum Fore to single record", m_optionsHolder );
  m_sumForeToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumForeToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_sumBackToSingleRecord = new WCheckBox( "Sum Back to single record", m_optionsHolder );
  m_sumBackToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumBackToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_sumSecoToSingleRecord = new WCheckBox( "Sum Sec to single record", m_optionsHolder );
  m_sumSecoToSingleRecord->checked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  m_sumSecoToSingleRecord->unChecked().connect( this, &ExportSpecFileTool::handleSumTypeToSingleRecordChanged );
  
  m_backSubFore = new WCheckBox( "Back. sub. Fore/Back", m_optionsHolder );
  m_backSubFore->checked().connect( this, &ExportSpecFileTool::handleBackSubForeChanged );
  m_backSubFore->unChecked().connect( this, &ExportSpecFileTool::handleBackSubForeChanged );
  
  m_sumDetsPerSample = new WCheckBox( "Sum Detectors per sample", m_optionsHolder );
  m_sumDetsPerSample->checked().connect( this, &ExportSpecFileTool::handleSumDetPerSampleChanged );
  m_sumDetsPerSample->unChecked().connect( this, &ExportSpecFileTool::handleSumDetPerSampleChanged );
  
  m_includeInterSpecInfo = new WCheckBox( "Include InterSpec info", m_optionsHolder );
  m_includeInterSpecInfo->checked().connect( this, &ExportSpecFileTool::handleIncludeInterSpecInfoChanged );
  m_includeInterSpecInfo->unChecked().connect( this, &ExportSpecFileTool::handleIncludeInterSpecInfoChanged );
  
  
  
  WContainerWidget *footer = new WContainerWidget( this );
  footer->addStyleClass( "ExportSpecFileFooter" );
  
  
  m_msg = new WText( "Some msg to be updated later", footer );
  m_msg->addStyleClass( "ExportSpecMsg" );
  
  WContainerWidget *btnsDiv = new WContainerWidget( footer );
  btnsDiv->addStyleClass( "ExportSpecBtns" );
  
  
  WPushButton *cancel_btn = new WPushButton( "Cancel", btnsDiv );
  cancel_btn->clicked().connect( boost::bind(&ExportSpecFileTool::emitDone, this, false) );
  
  if( !m_resource )
  {
    m_resource = new ExportSpecFileTool_imp::DownloadSpectrumResource( this );
    m_resource->setTakesUpdateLock( true );
    m_resource->downloadFinished().connect( boost::bind(&ExportSpecFileTool::emitDone, this, true) );
  }//if( !m_resource )
  
  m_export_btn = new WPushButton( "Export", btnsDiv );
  m_export_btn->setLink( WLink(m_resource) );
  m_export_btn->disable();
  
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
  const bool disable = (m_fileSelect && (m_fileSelect->currentIndex() < 0));
  m_export_btn->setDisabled( disable );
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
    
    tabletxt << "<tr><th>Start Time</th><td>" << start_date << "</td></tr>\n"
                "<tr><th></th><td>" << start_tod << "</td></tr>\n";
    
    if( (end_time - start_time) > 2*std::chrono::minutes(15) )
    {
      const string end_str = SpecUtils::to_vax_string( end_time );
      pos = end_str.find( ' ' );
      const string end_date = end_str.substr(0, pos);
      const string end_tod = (pos == string::npos) ? string() : end_str.substr(pos+1);
      
      tabletxt << "<tr><th>End Time</th><td>" << end_date << "</td></tr>\n"
                  "<tr><th></th><td>" << end_tod << "</td></tr>\n";
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
    tabletxt << "<tr><th>Neut. Dets</th><td>"
             << static_cast<int>(meas->neutron_detector_names().size()) << "</td></tr>\n";
    tabletxt << "<tr><th>Gamma Dets</th><td>"
            << static_cast<int>(meas->gamma_detector_names().size()) << "</td></tr>\n";
  }else
  {
    tabletxt << "<tr><th>Num Dets</th><td>"
             << static_cast<int>(meas->detector_names().size()) << "</td></tr>\n";
  }
  
  if( meas->has_gps_info() )
  {
    tabletxt << "<tr><th>Latitude</th><td>"
             << PhysicalUnits::printCompact( meas->mean_latitude(), 7 ) << "</td></tr>\n";
    tabletxt << "<tr><th>Longitude</th><td>"
             << PhysicalUnits::printCompact( meas->mean_longitude(), 7 ) << "</td></tr>\n";
  }//if( meas->has_gps_info() )
  
  const string total_time = PhysicalUnits::printToBestTimeUnits( meas->gamma_real_time() );
  tabletxt << "<tr><th>Total Time</th><td>" << total_time << "</td></tr>\n";
  
  tabletxt << "<tr><th>Sum Gamma</th><td>"
           << PhysicalUnits::printCompact( meas->gamma_count_sum(), 5) << "</td></tr>\n";
  if( meas->contained_neutron() )
    tabletxt << "<tr><th>Sum Neut.</th><td>"
             << PhysicalUnits::printCompact( meas->neutron_counts_sum(), 5) << "</td></tr>\n";
  
  //const string &instrument_type = meas->instrument_type();
  const string &manufacturer = meas->manufacturer();
  const string &instrument_model = meas->instrument_model();
  const string &instrument_id = meas->instrument_id();
  
  if( !manufacturer.empty() )
    tabletxt << "<tr><th>Manufacturer</th><td>"
             << Wt::Utils::htmlEncode(manufacturer) << "</td></tr>\n";
  
  if( meas->detector_type() != SpecUtils::DetectorType::Unknown )
    tabletxt << "<tr><th>Model</th><td>"
             << SpecUtils::detectorTypeToString(meas->detector_type()) << "</td></tr>\n";
  else if( !instrument_model.empty() )
    tabletxt << "<tr><th>Model</th><td>"
             << Wt::Utils::htmlEncode(instrument_model) << "</td></tr>\n";
  
  if( !instrument_id.empty() )
    tabletxt << "<tr><th>Serial</th><td>" << Wt::Utils::htmlEncode(instrument_id) << "</td></tr>\n";
  
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


SpecUtils::SaveSpectrumAsType ExportSpecFileTool::currentSaveType() const
{
  const WMenuItem * const currentFormatItem = m_formatMenu->currentItem();
  assert( currentFormatItem );
  
  if( currentFormatItem )
  {
    const uint64_t data = reinterpret_cast<uint64_t>( currentFormatItem->data() );
    assert( data <= static_cast<int>(SpecUtils::SaveSpectrumAsType::NumTypes) );
    return SpecUtils::SaveSpectrumAsType( data );
  }//if( currentFormatItem )
  
  return SpecUtils::SaveSpectrumAsType::N42_2012;
}//SpecUtils::SaveSpectrumAsType currentSaveType() const;


uint16_t ExportSpecFileTool::maxRecordsInCurrentSaveType() const
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
      break;
    
    // Spectrum file types that can have two spectra in them (foreground + background)
    case SpecUtils::SaveSpectrumAsType::NumTypes:
#if( SpecUtils_ENABLE_D3_CHART )
    case SpecUtils::SaveSpectrumAsType::HtmlD3:
#endif
      return 2;
      break;
      
    // Spectrum file types that can have a single spectra in them
    case SpecUtils::SaveSpectrumAsType::Chn:
    case SpecUtils::SaveSpectrumAsType::SpcBinaryInt:
    case SpecUtils::SaveSpectrumAsType::SpcBinaryFloat:
    case SpecUtils::SaveSpectrumAsType::SpcAscii:
    case SpecUtils::SaveSpectrumAsType::SpeIaea:
    case SpecUtils::SaveSpectrumAsType::Cnf:
    case SpecUtils::SaveSpectrumAsType::Tka:
      return 1;
      break;
    
#if( SpecUtils_INJA_TEMPLATES )
    case SpecUtils::SaveSpectrumAsType::Template:
      assert( 0 );
      break;
#endif
  }//switch( save_format )
  
  assert( 0 );
  
  return 0;
}//uint16_t maxRecordsInCurrentSaveType() const


void ExportSpecFileTool::refreshSampleAndDetectorOptions()
{
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  const uint16_t max_records = maxRecordsInCurrentSaveType();
  const SpecUtils::SaveSpectrumAsType save_type = currentSaveType();
  
  if( !spec || (spec->gamma_detector_names().size() <= 1) )
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
    }//
    
    if( (max_records <= 1) || (max_records == 2) )
    {
      m_sumAllToSingleRecord->hide();
      m_sumForeToSingleRecord->hide();
      m_sumBackToSingleRecord->hide();
      m_sumSecoToSingleRecord->hide();
    }else
    {
      m_sumAllToSingleRecord->show();
      m_sumForeToSingleRecord->setHidden( !(m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked()) );
      m_sumBackToSingleRecord->setHidden( !(m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked()) );
      m_sumSecoToSingleRecord->setHidden( !(m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked()) );
    }
  }//if( spec->gamma_detector_names().size() <= 1 ) / else
  
  
  const bool use_fore_disp = (m_dispForeSamples->isVisible() && m_dispForeSamples->isChecked());
  const bool use_seco_disp = (m_dispSecondSamples->isVisible() && m_dispSecondSamples->isChecked());
  const bool use_back_disp = (m_dispBackSamples->isVisible() && m_dispBackSamples->isChecked());
  
  if( (use_fore_disp || use_seco_disp) && use_back_disp )
  {
    m_backSubFore->show();
    if( use_fore_disp && use_seco_disp )
      m_backSubFore->setText( "Back. Sub. For./Sec." );
    else if( use_fore_disp )
      m_backSubFore->setText( "Back. Sub. For." );
    else
      m_backSubFore->setText( "Back. Sub. Sec." );
  }else
  {
    m_backSubFore->hide();
  }
  
  const bool can_save_interspec_info = ((save_type == SpecUtils::SaveSpectrumAsType::N42_2012)
                                        || (save_type == SpecUtils::SaveSpectrumAsType::N42_2006));
  m_includeInterSpecInfo->setHidden( !can_save_interspec_info );
  
  if( !spec || (max_records <= 2) || (spec->gamma_detector_names().size() <= 1) )
  {
    m_sumDetsPerSample->hide();
  }else
  {
    m_sumDetsPerSample->show();
  }
  
  
  if( !spec || (spec->sample_numbers().size() <= 1) )
  {
    m_samplesHolder->setHidden( (m_detectorFilterCbs->children().size() == 0) );
    
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
    m_includeInterSpecInfo->hide();
    
    handleSamplesChanged();
    return;
  }//if( no reason to show selecting samples )
  
  m_samplesHolder->setHidden(false);
  
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
  
  handleSamplesChanged();
}//void ExportSpecFileTool::refreshSampleAndDetectorOptions()


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
  const string txt = m_customSamplesEdit->text().toUTF8();
  
  const shared_ptr<const SpecMeas> spec = currentlySelectedFile();
  const pair<set<int>,string> samples_txt = sampleNumbersFromTxtRange( txt, spec, false );
  m_customSamplesEdit->setText( samples_txt.second );
  
  handleSamplesChanged();
}//void handleCustomSampleTxtChanged()


void ExportSpecFileTool::handleSamplesChanged()
{
  if( !m_dispForeSamples->isChecked()
     && !m_dispBackSamples->isChecked()
     && !m_dispSecondSamples->isChecked()
     && !m_customSamples->isChecked()
     && !m_allSamples->isChecked() )
  {
    m_allSamples->setChecked( true );
  }
}//void handleSamplesChanged()


void ExportSpecFileTool::handleFormatChange()
{
  refreshSampleAndDetectorOptions();
}//void handleFormatChange();
  

void ExportSpecFileTool::handleForePlusBackChanged()
{
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
  assert( m_filterDetector );
  assert( m_detectorFilterCbs );
  if( !m_filterDetector || !m_detectorFilterCbs )
    return;
  
  if( m_filterDetector->isChecked() )
  {
    m_detectorFilterCbs->show();
  }else
  {
    m_detectorFilterCbs->hide();
    for( const auto w : m_detectorFilterCbs->children() )
    {
      WCheckBox *cb = dynamic_cast<WCheckBox *>( w );
      if( cb )
        cb->setChecked( true );
    }
  }//if( !m_filterDetector->isChecked() )
  
  // TODO: update options
}//void handleFilterDetectorCbChanged()


void ExportSpecFileTool::handleSumToSingleRecordChanged()
{
  if( m_sumAllToSingleRecord->isChecked() )
  {
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
  }//if( m_sumAllToSingleRecord->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleSumToSingleRecordChanged()


void ExportSpecFileTool::handleSumTypeToSingleRecordChanged()
{
  if( m_sumForeToSingleRecord->isChecked()
     || m_sumBackToSingleRecord->isChecked()
     || m_sumSecoToSingleRecord->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
  }
  
  refreshSampleAndDetectorOptions();
}//void handleSumTypeToSingleRecordChanged()


void ExportSpecFileTool::handleBackSubForeChanged()
{
  if( m_backSubFore->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_sumDetsPerSample->setChecked( false );
  }//if( m_backSubFore->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleBackSubForeChanged()


void ExportSpecFileTool::handleSumDetPerSampleChanged()
{
  if( m_sumDetsPerSample->isChecked() )
  {
    m_sumAllToSingleRecord->setChecked( false );
    m_sumForeToSingleRecord->setChecked( false );
    m_sumBackToSingleRecord->setChecked( false );
    m_sumSecoToSingleRecord->setChecked( false );
    m_backSubFore->setChecked( false );
  }//if( m_sumDetsPerSample->isChecked() )
  
  refreshSampleAndDetectorOptions();
}//void handleSumDetPerSampleChanged()


void ExportSpecFileTool::handleIncludeInterSpecInfoChanged()
{
  // Nothing to do here I think
}//void handleIncludeInterSpecInfoChanged()


void ExportSpecFileTool::handleAppUrl( std::string query_str )
{
  passMessage( "ExportSpecFileTool::handleAppUrl not implemented yet", WarningWidget::WarningMsgHigh );
}//void handleAppUrl( std::string query_str )


std::string ExportSpecFileTool::encodeStateToUrl() const
{
  passMessage( "ExportSpecFileTool::encodeStateToUrl not implemented yet", WarningWidget::WarningMsgHigh );
  return "";
}//std::string encodeStateToUrl() const


ExportSpecFileWindow::ExportSpecFileWindow( InterSpec *viewer )
  : SimpleDialog( "Spectrum File Export", "" ),
  m_tool( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/ExportSpecFile.css" );
  
  addStyleClass( "export-spec-file" );
  
  const int w = viewer->renderedWidth();
  setMinimumSize( WLength(w > 100 ? std::min(0.95*w, 800.0) : 800.0 ,WLength::Pixel), WLength::Auto );
  
  m_tool = new ExportSpecFileTool( viewer, contents() );
  m_tool->done().connect( boost::bind(&ExportSpecFileWindow::accept, this) );
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


