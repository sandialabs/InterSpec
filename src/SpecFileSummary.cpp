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

#include <cfloat>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <Wt/WText.h>
#include <Wt/WImage.h>
#include <Wt/WLabel.h>
#include <Wt/WTable.h>
#include <Wt/WServer.h>
#include <Wt/WLineEdit.h>
#include <Wt/WTextArea.h>
#include <Wt/WComboBox.h>
#include <Wt/WGroupBox.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WRadioButton.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WApplication.h>
#include <Wt/WIntValidator.h>
#include <Wt/WSelectionBox.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecFileSummary.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ShowRiidInstrumentsAna.h"

#if( USE_GOOGLE_MAP )
#include "InterSpec/GoogleMap.h"
#elif( USE_LEAFLET_MAP )
#include "InterSpec/LeafletRadMap.h"
#endif

using namespace std;
using namespace Wt;

namespace
{
#if( USE_GOOGLE_MAP )
  void updateCoordText( GoogleMap *googlemap,
                        double lat, double lng,
                        double origLat, double origLng )
  {
    char buffer[128];
    snprintf( buffer, sizeof(buffer), "(%.6f, %.6f)", lat, lng );
    
    googlemap->clearMeasurments();
    googlemap->addMarker( origLat, origLng );
    googlemap->addInfoBox( lat, lng, buffer );
  }//updateCoordText(...)
#endif
  
  double coordVal( Wt::WString instr )
  {
#ifndef WT_NO_STD_WSTRING
    std::wstring input = instr.value();
#else
    std::string input = instr.toUTF8();
#endif
    
    double answer = -999.9f;
    
    boost::algorithm::to_lower( input );
    boost::algorithm::replace_all(input, L"\x00B0", L" ");
    boost::algorithm::replace_all(input, L"'", L" ");
    boost::algorithm::replace_all(input, L"\"", L" ");
    boost::algorithm::replace_all(input, L"\t", L" ");
    boost::algorithm::replace_all(input, L"degree", L" ");
    boost::algorithm::replace_all(input, L"deg.", L" ");
    boost::algorithm::replace_all(input, L"deg", L" ");
    
    double sign = 1.0f;
    
    //Go through and look for n s e w to determine sign
#ifndef WT_NO_STD_WSTRING
    const wchar_t *dirs[] = { L"n", L"s", L"e", L"w" };
#else
    const char *dirs[] = { "n", "s", "e", "w" };
#endif
    
    for( size_t i = 0; i < 4; ++i )
    {
      size_t pos = input.find_first_of( dirs[i] );
      if( pos != wstring::npos )
      {
        input.erase( pos, 1 );
        sign = std::pow( -1.0, double(i) );
        break;
      }
    }//for( size_t i = 0; i < 4; ++i )
    
    //At this point, we should have nothing besides digits, spaces, decimals, and signs
    for( size_t i = 0; i < input.size(); ++i )
      if( !isdigit(input[i]) && input[i]!=' ' && input[i]!= '.' && input[i]!= '-' )
        return answer;
    
    boost::algorithm::trim( input );
    vector<std::wstring> fields;
    boost::algorithm::split( fields, input, boost::is_any_of( L" " ),
                            boost::token_compress_on );
    
    const size_t nfields = fields.size();
    
    if( nfields == 1 )
    {
      try
      {
        answer = std::stod( fields[0] );
        return answer;
      }catch(...){}
    }//if( nfields == 1 )
    
    
    if( nfields==3 )
    {
      try
      {
        const float deg = static_cast<float>( std::stod( fields[0] ) );
        const float min = static_cast<float>( std::stod( fields[1] ) );
        const float sec = static_cast<float>( std::stod( fields[2] ) );
        answer = deg + (min/60.0f) + (sec/3600.0f);
        return sign * answer;
      }catch(...){}
    }//if( nfields==3 )
    
    //If were here, weve failed
    return answer;
  }//float coordVal( std::wstring input )
}//namespace


SpecFileSummary::SpecFileSummary( const SpecUtils::SpectrumType type, InterSpec *specViewer )
  : AuxWindow( WString::tr("window-title-file-parameters"),
              (AuxWindowProperties::IsModal
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse) ),
    m_specViewer( specViewer ),
    m_allowEditGroup( NULL ),
    m_displaySampleDiv( NULL ),
    m_displaySampleNumEdit( NULL ),
    m_displayedPreText( NULL ),
    m_displayedPostText( NULL ),
    m_nextSampleNumButton( NULL ),
    m_prevSampleNumButton( NULL ),
    m_displaySampleNumValidator( NULL ),
    m_gammaCPS( NULL ),
    m_gammaSum( NULL ),
    m_neutronCPS( NULL ),
    m_neutronSum( NULL ),
    m_displayedLiveTime( NULL ),
    m_displayedRealTime( NULL ),
    m_timeStamp( NULL ),
    m_energyRange( NULL ),
    m_numberOfBins( NULL ),
    m_detector( NULL ),
    m_sampleNumber( NULL ),
    m_measurmentRemarks( NULL ),
    m_longitude( NULL ),
    m_latitude( NULL ),
    m_gpsTimeStamp( NULL ),
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
    m_showMapButton( NULL ),
#endif
    m_title( NULL ),
    m_source( NULL ),
    m_spectraGroup( NULL ),
    m_fileRemarks( NULL ),
    m_sizeInMemmory( NULL ),
    m_filename( NULL ),
    m_uuid( NULL ),
    m_laneNumber( NULL ),
    m_measurement_location_name( NULL ),
    m_ana_button( NULL ),
    m_multimedia_button( nullptr ),
    m_inspection( NULL ),
    m_instrument_type( NULL ),
    m_manufacturer( NULL ),
    m_instrument_model( NULL ),
    m_instrument_id( NULL ),
    m_reloadSpectrum( NULL )
{
  init();
  
  if( m_spectraGroup->selectedButtonIndex() != static_cast<int>(type) )
  {
    m_spectraGroup->setSelectedButtonIndex( static_cast<int>(type) );
    handleSpectrumTypeChanged();
  }
}//SpecFileSummary constructor


template <class T>
void addField( T *&edit, WGridLayout *table, const WString &labelstr,
               int row, int col, int rowspan = 1, int colspan = 1 )
{
  WLabel *label = table->addWidget( std::make_unique<WLabel>( labelstr ), row, col,
                                    Wt::AlignmentFlag::Middle );
  label->setStyleClass( "SpectrumFileSummaryLabel" );
  auto editOwner = std::make_unique<T>();
  edit = editOwner.get();
  table->addWidget( std::move(editOwner), row, col+1, rowspan, colspan );
  label->setBuddy( edit );
}//addField(...)


void addField( WText *&edit, WGridLayout *table, const WString &labelstr,
              int row, int col, int rowspan = 1, int colspan = 1 )
{
  WLabel *label = table->addWidget( std::make_unique<WLabel>( labelstr ), row, col,
                                    (rowspan==1 ? Wt::AlignmentFlag::Middle : Wt::AlignmentFlag::Bottom) );
  label->setStyleClass( "SpectrumFileSummaryLabel" );
  auto editOwner = std::make_unique<WText>();
  edit = editOwner.get();
  table->addWidget( std::move(editOwner), row, col+1, rowspan, colspan,
                    Wt::AlignmentFlag::Middle | Wt::AlignmentFlag::Left );
}



void SpecFileSummary::init()
{
  wApp->useStyleSheet( "InterSpec_resources/SpecFileSummary.css" );
  
  if( m_specViewer )
    m_specViewer->useMessageResourceBundle( "SpecFileSummary" );
  
  WContainerWidget *holder = contents()->addNew<WContainerWidget>();

  WGridLayout *overallLayout = holder->setLayout( std::make_unique<WGridLayout>() );
  holder->addStyleClass( "SpecFileSummary" );

  overallLayout->setContentsMargins( 9, 0, 9, 0 );

  WRadioButton *button = nullptr;
  std::unique_ptr<WGroupBox> spectrumGroupBoxOwner = std::make_unique<WGroupBox>( WString::tr("sfs-spectrum") );
  WGroupBox *spectrumGroupBox = spectrumGroupBoxOwner.get();
  spectrumGroupBox->addStyleClass( "SpecSummChoose" );
  // In Wt4, WGroupBox inside WGridLayout requires a layout to render children
  WGridLayout *spectrumBoxLayout = spectrumGroupBox->setLayout( std::make_unique<WGridLayout>() );
  spectrumBoxLayout->setContentsMargins( 0, 0, 0, 0 );
  m_spectraGroup = std::make_shared<WButtonGroup>();
  button = spectrumBoxLayout->addWidget( std::make_unique<WRadioButton>( WString::tr("Foreground") ), 0, 0 );
  button->addStyleClass( "CbNoLineBreak" );
  m_spectraGroup->addButton( button, static_cast<int>(SpecUtils::SpectrumType::Foreground) );
  button = spectrumBoxLayout->addWidget( std::make_unique<WRadioButton>( WString::tr("Secondary") ), 0, 1 );
  button->addStyleClass( "CbNoLineBreak" );
  m_spectraGroup->addButton( button, static_cast<int>(SpecUtils::SpectrumType::SecondForeground) );
  button = spectrumBoxLayout->addWidget( std::make_unique<WRadioButton>( WString::tr("Background") ), 0, 2 );
  button->addStyleClass( "CbNoLineBreak" );
  m_spectraGroup->addButton( button, static_cast<int>(SpecUtils::SpectrumType::Background) );
  m_spectraGroup->setCheckedButton( m_spectraGroup->button( static_cast<int>(SpecUtils::SpectrumType::Foreground) ) );
  m_spectraGroup->checkedChanged().connect( this, &SpecFileSummary::handleSpectrumTypeChanged );


  std::unique_ptr<WGroupBox> editGroupBoxOwner = std::make_unique<WGroupBox>( WString::tr("sfs-alow-edit-cb") );
  WGroupBox *editGroupBox = editGroupBoxOwner.get();
  editGroupBox->addStyleClass( "SpecSummAllowEdit" );
  // In Wt4, WGroupBox inside WGridLayout requires a layout to render children
  WGridLayout *editBoxLayout = editGroupBox->setLayout( std::make_unique<WGridLayout>() );
  editBoxLayout->setContentsMargins( 0, 0, 0, 0 );
  m_allowEditGroup = std::make_shared<WButtonGroup>();
  button = editBoxLayout->addWidget( std::make_unique<WRadioButton>( WString::tr("Yes") ), 0, 0 );
  button->addStyleClass( "CbNoLineBreak" );
  m_allowEditGroup->addButton( button, kAllowModify );
  button = editBoxLayout->addWidget( std::make_unique<WRadioButton>( WString::tr("No") ), 0, 1 );
  button->addStyleClass( "CbNoLineBreak" );
  m_allowEditGroup->addButton( button, kDontAllowModify );
  m_allowEditGroup->setCheckedButton( m_allowEditGroup->button(kDontAllowModify) );
  m_allowEditGroup->checkedChanged().connect( this, &SpecFileSummary::handleAllowModifyStatusChange );

  auto upperdivOwner = std::make_unique<WContainerWidget>();
  WContainerWidget *upperdiv = upperdivOwner.get();
  WGridLayout *upperlayout = upperdiv->setLayout( std::make_unique<WGridLayout>() );
  overallLayout->addWidget( std::move(upperdivOwner), 0, 0 );
  upperlayout->addWidget( std::move(spectrumGroupBoxOwner), 0, 0, 1, 1 );
  upperlayout->addWidget( std::move(editGroupBoxOwner), 0, 1, 1, 1 );
  upperlayout->setColumnStretch( 1, 1 );

  m_reloadSpectrum = editBoxLayout->addWidget( std::make_unique<WPushButton>( WString::tr("sfs-update-display-btn") ), 0, 2 );
  m_reloadSpectrum->setToolTip( WString::tr("sfs-tt-update-display") );
  m_reloadSpectrum->setIcon( WLink("InterSpec_resources/images/arrow_refresh.svg") );
  m_reloadSpectrum->clicked().connect( this, &SpecFileSummary::reloadCurrentSpectrum );
  m_reloadSpectrum->disable();
  

  auto filedivOwner = std::make_unique<WGroupBox>( WString::tr("sfs-file-info-title") );
  WGroupBox *filediv = filedivOwner.get();
  filediv->addStyleClass( "SpecFileInfo" );
  WGridLayout *fileInfoLayout = filediv->setLayout( std::make_unique<WGridLayout>() );
  overallLayout->addWidget( std::move(filedivOwner), 1, 0 );

  addField( m_filename, fileInfoLayout, WString::tr("sfs-file-name-label"), 0, 0, 1, 5 );
  addField( m_sizeInMemmory, fileInfoLayout, WString::tr("sfs-memory-size-label"), 0, 6 );
  addField( m_inspection, fileInfoLayout, WString::tr("sfs-inspection-label"), 1, 0 );
  addField( m_laneNumber, fileInfoLayout, WString::tr("sfs-lane-label"), 1, 2 );
  addField( m_measurement_location_name, fileInfoLayout, WString::tr("sfs-location-label"), 1, 4, 1, 1 );

  auto ana_btn = std::make_unique<WPushButton>( WString::tr("sfs-no-det-id-btn") );
  m_ana_button = ana_btn.get();
  m_ana_button->disable();
  m_ana_button->clicked().connect( this, &SpecFileSummary::showRiidAnalysis );
  fileInfoLayout->addWidget( std::move(ana_btn), 1, 6, 1, 1, Wt::AlignmentFlag::Left );

  auto multimedia_btn = std::make_unique<WPushButton>( WString::tr("sfs-show-pics-btn") );
  m_multimedia_button = multimedia_btn.get();
  m_multimedia_button->hide();
  m_multimedia_button->clicked().connect( this, &SpecFileSummary::showMultimedia );
  fileInfoLayout->addWidget( std::move(multimedia_btn), 1, 7, 1, 1, Wt::AlignmentFlag::Left );

  addField( m_instrument_type, fileInfoLayout, WString::tr("sfs-instrument-type-label"), 2, 0 );
  addField( m_manufacturer, fileInfoLayout, WString::tr("sfs-manufacturer-label"), 2, 2 );
  addField( m_instrument_model, fileInfoLayout, WString::tr("sfs-model-label"), 2, 4, 1, 3 );

  addField( m_instrument_id, fileInfoLayout, WString::tr("sfs-instrument-id-label"), 3, 0, 1, 3 );
  addField( m_uuid, fileInfoLayout, WString::tr("sfs-uuid-label"), 3, 4, 1, 3 );

  auto labelHolderOwner = std::make_unique<WContainerWidget>();
  WContainerWidget *labelHolder = labelHolderOwner.get();
  labelHolder->setAttributeValue( "style", "margin-top: auto; margin-bottom: auto;" );
  WLabel *label = labelHolder->addNew<WLabel>( WString::tr("sfs-file-label") );
  label->setStyleClass( "SpectrumFileSummaryLabel" );
  label->setInline( false );
  label = labelHolder->addNew<WLabel>( WString::tr("sfs-remarks-label") );
  label->setStyleClass( "SpectrumFileSummaryLabel" );
  label->setInline( false );
  fileInfoLayout->addWidget( std::move(labelHolderOwner), 4, 0, Wt::AlignmentFlag::Middle );

  auto remarkHolderOwner = std::make_unique<WContainerWidget>();
  WContainerWidget *remarkHolder = remarkHolderOwner.get();
  fileInfoLayout->addWidget( std::move(remarkHolderOwner), 4, 1, 1, 7 );
  m_fileRemarks = remarkHolder->addNew<WTextArea>();
  m_fileRemarks->setAttributeValue( "style", "resize: none; width: 100%; height: 100%;" );

  fileInfoLayout->setRowStretch( 4, 1 );

  fileInfoLayout->setColumnStretch( 1, 1 );
  fileInfoLayout->setColumnStretch( 3, 1 );
  fileInfoLayout->setColumnStretch( 5, 1 );
  fileInfoLayout->setColumnStretch( 7, 1 );

  
  auto measdivOwner = std::make_unique<WGroupBox>( WString::tr("sfs-measurement-info-label") );
  WGroupBox *measdiv = measdivOwner.get();
  measdiv->addStyleClass( "MeasInfo" );
  WGridLayout *measTable = measdiv->setLayout( std::make_unique<WGridLayout>() );
  overallLayout->addWidget( std::move(measdivOwner), 2, 0 );
  
  
  addField( m_timeStamp, measTable, WString::tr("sfs-date-time-label"), 0, 0, 1, 3 );
  addField( m_displayedLiveTime, measTable, WString::tr("sfs-live-time-label"), 0, 4 );
  addField( m_displayedRealTime, measTable, WString::tr("sfs-real-time-label"), 0, 6 );
  addField( m_detector, measTable, WString::tr("sfs-detector-name-label"), 1, 0 );
  addField( m_sampleNumber, measTable, WString::tr("sfs-sample-number-label"), 1, 2 );
  addField( m_energyRange, measTable, WString::tr("sfs-energy-range-label"), 1, 4 );
  addField( m_numberOfBins, measTable, WString::tr("sfs-number-channels-label"), 1, 6 );
  addField( m_gammaSum, measTable, WString::tr("sfs-sum-gammas-label"), 2, 0 );
  addField( m_gammaCPS, measTable, WString::tr("sfs-gamma-cps-label"), 2, 2 );
  addField( m_neutronSum, measTable, WString::tr("sfs-sum-neutrons-label"), 2, 4 );
  addField( m_neutronCPS, measTable, WString::tr("sfs-neutrons-cps-label"), 2, 6 );
  addField( m_latitude, measTable, WString::tr("sfs-latitude-label"), 3, 0 );
  addField( m_longitude, measTable, WString::tr("sfs-longitude-label"), 3, 2 );
  addField( m_gpsTimeStamp, measTable, WString::tr("sfs-position-time-label"), 3, 4 );
  
  m_latitude->setPlaceholderText( WString::tr("sfs-latitude-empty-text") );
  m_longitude->setPlaceholderText( WString::tr("sfs-longitude-empty-text") );
  
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
  {
    auto mapBtnOwner = std::make_unique<WPushButton>( WString::tr("sfs-show-map-btn") );
    m_showMapButton = mapBtnOwner.get();
    measTable->addWidget( std::move(mapBtnOwner), 3, 6 );
#if( USE_GOOGLE_MAP  )
    m_showMapButton->clicked().connect( this, &SpecFileSummary::showGoogleMap );
#elif( USE_LEAFLET_MAP )
    m_showMapButton->clicked().connect( this, &SpecFileSummary::showLeafletMap );
#endif
    m_showMapButton->disable();
  }
#endif
  
  addField( m_title, measTable, WString::tr("sfs-description-label"), 4, 0, 1, 5 );
  
  addField( m_source, measTable, WString::tr("sfs-source-type-label"), 4, 6 );
  m_source->setNoSelectionEnabled( true );
  
  for( SpecUtils::SourceType i = SpecUtils::SourceType(0);
       i <= SpecUtils::SourceType::Unknown;
       i = SpecUtils::SourceType(static_cast<int>(i)+1) )
  {
    WString val;
    switch( i )
    {
      case SpecUtils::SourceType::IntrinsicActivity: val = WString::tr("intrinsic-activity"); break;
      case SpecUtils::SourceType::Calibration:       val = WString::tr("Calibration");        break;
      case SpecUtils::SourceType::Background:        val = WString::tr("Background");         break;
      case SpecUtils::SourceType::Foreground:        val = WString::tr("Foreground");         break;
      case SpecUtils::SourceType::Unknown: val = WString::tr("Unknown");            break;
      default: break;
    }//switch( i )
    
    m_source->addItem( val );
  }//for( int i = 0; i <= SpecUtils::SourceType::Unknown; ++i )
  

  {
    auto lhOwner = std::make_unique<WContainerWidget>();
    WContainerWidget *lh = lhOwner.get();
    lh->setAttributeValue( "style", "margin-top: auto; margin-bottom: auto;" );
    label = lh->addNew<WLabel>( WString::tr("sfs-spectra-label") );
    label->setStyleClass( "SpectrumFileSummaryLabel" );
    label->setInline( false );
    label = lh->addNew<WLabel>( WString::tr("sfs-remarks-label") );
    label->setStyleClass( "SpectrumFileSummaryLabel" );
    label->setInline( false );
    measTable->addWidget( std::move(lhOwner), 5, 0, Wt::AlignmentFlag::Middle );

    auto rhOwner = std::make_unique<WContainerWidget>();
    WContainerWidget *rh = rhOwner.get();
    measTable->addWidget( std::move(rhOwner), 5, 1, 1, 7 );
    m_measurmentRemarks = rh->addNew<WTextArea>();
    m_measurmentRemarks->setAttributeValue( "style", "resize: none; width: 100%; height: 100%;" );
  }

  measTable->setRowStretch( 5, 1 );

  measTable->setColumnStretch( 1, 1 );
  measTable->setColumnStretch( 3, 1 );
  measTable->setColumnStretch( 5, 1 );
  measTable->setColumnStretch( 7, 1 );


  {
    auto sampleDivOwner = std::make_unique<WContainerWidget>();
    m_displaySampleDiv = sampleDivOwner.get();
    m_displaySampleDiv->addStyleClass( "displaySampleDiv" );
    m_prevSampleNumButton = m_displaySampleDiv->addNew<WImage>( "InterSpec_resources/images/previous_arrow.png" );
    m_prevSampleNumButton->clicked().connect( [this](){ handleUserIncrementSampleNum( false ); } );
    m_displayedPreText = m_displaySampleDiv->addNew<WText>();
    m_displaySampleNumEdit = m_displaySampleDiv->addNew<WLineEdit>();

    m_displaySampleNumEdit->setAutoComplete( false );
    m_displaySampleNumEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_displaySampleNumEdit->setAttributeValue( "autocorrect", "off" );
    m_displaySampleNumEdit->setAttributeValue( "spellcheck", "off" );
#endif

    auto validatorOwner = std::make_shared<WIntValidator>();
    m_displaySampleNumValidator = validatorOwner.get();
    m_displaySampleNumEdit->setValidator( validatorOwner );
    m_displaySampleNumEdit->addStyleClass( "numberValidator" ); //used to detect mobile keyboard
    m_displaySampleNumEdit->setTextSize( 3 );
    m_displayedPostText = m_displaySampleDiv->addNew<WText>();
    m_nextSampleNumButton = m_displaySampleDiv->addNew<WImage>( "InterSpec_resources/images/next_arrow.png" );
    m_nextSampleNumButton->clicked().connect( [this](){ handleUserIncrementSampleNum( true ); } );
    m_displaySampleNumEdit->enterPressed().connect( [this](){ handleUserChangeSampleNum(); } );
    m_displaySampleNumEdit->blurred().connect( [this](){ handleUserChangeSampleNum(); } );
    m_displaySampleDiv->setHiddenKeepsGeometry( false );

    measTable->addWidget( std::move(sampleDivOwner), measTable->rowCount(), 0, 1, measTable->columnCount(),
                          Wt::AlignmentFlag::Bottom );
  }

  overallLayout->setRowStretch( 1, 1 );
  overallLayout->setRowStretch( 2, 1 );

  m_specViewer->displayedSpectrumChanged().connect(
                [this]( SpecUtils::SpectrumType t, std::shared_ptr<SpecMeas> m,
                        std::set<int> s, std::vector<std::string> d ){
                  handleSpectrumChange( t, m, s, d );
                } );

  // Helper lambda to connect all three signals (changed, enterPressed, blurred) to handleFieldUpdate
  auto connectField = [this]( WFormWidget *w, EditableFields f ){
    w->changed().connect( [this, f](){ handleFieldUpdate( f ); } );
    w->enterPressed().connect( [this, f](){ handleFieldUpdate( f ); } );
    w->blurred().connect( [this, f](){ handleFieldUpdate( f ); } );
  };

  connectField( m_displayedLiveTime, kDisplayedLiveTime );
  connectField( m_displayedRealTime, kDisplayedRealTime );
  connectField( m_timeStamp, kTimeStamp );
  connectField( m_longitude, kPosition );
  connectField( m_latitude, kPosition );
  connectField( m_gpsTimeStamp, kPosition );
  connectField( m_title, kDescription );

  m_source->activated().connect( [this]( int ){ handleFieldUpdate( kSourceType ); } );

  connectField( m_measurmentRemarks, kMeasurmentRemarks );
  connectField( m_fileRemarks, kFileRemarks );
  connectField( m_filename, kFilename );
  connectField( m_uuid, kUuid );
  connectField( m_laneNumber, kLaneNumber );
  connectField( m_measurement_location_name, kMeasurement_location_name );
  connectField( m_inspection, kInspection );
  connectField( m_instrument_type, kInstrument_type );
  connectField( m_manufacturer, kManufacturer );
  connectField( m_instrument_model, kInstrument_model );
  connectField( m_instrument_id, kInstrument_id );


  handleSpectrumTypeChanged();
  handleAllowModifyStatusChange();
  finished().connect( [this](){ AuxWindow::deleteAuxWindow( this ); } );

  AuxWindow::addHelpInFooter( footer(), "file-parameters-dialog" );

  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( [this](){ AuxWindow::deleteAuxWindow( this ); } );
  
  holder->setMinimumSize( 780, 500 );
  holder->setWidth( WLength(100.0,WLength::Unit::Percentage) );
  holder->setHeight( WLength(100.0,WLength::Unit::Percentage) );
  WDialog::contents()->setOverflow( Wt::Overflow::Auto );
  
  //set min size so setResizable call before setResizable so Wt/Resizable.js wont cause the initial
  //  size to be the min-size
  setMinimumSize( 640, 480 );
  
//  setMinimumSize( 780, 500 );
  resizeScaledWindow(0.85,0.85);
  refresh();
  
  centerWindow();
  rejectWhenEscapePressed();
  setResizable( true );
  show();
}//void SpecFileSummary::init()


SpecFileSummary::~SpecFileSummary()
{
  //nothing to do here
}


void SpecFileSummary::showRiidAnalysis()
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  
  // The SimpleDialog containing the RIID results will be behind this SpecFileSummary, if we create
  //  it here (I assume do to the mouse click bringing *this to front), so we will create it in the
  //  next event loop iteration, and this way it appears good.
  WServer::instance()->post( wApp->sessionId(), std::bind([=](){
    auto app = WApplication::instance();
    if( app ){
      InterSpec::instance()->showRiidResults( type );
      app->triggerUpdate();
    }
  }) );
}//void showRiidAnalysis()


void SpecFileSummary::showMultimedia()
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  m_specViewer->showMultimedia( type );
}//void showMultimedia();


#if( USE_GOOGLE_MAP )
void SpecFileSummary::showGoogleMap()
{
  if( !m_specViewer )
    return;
  
  double longitude = coordVal( m_longitude->text() );
  double latitude = coordVal( m_latitude->text() );
  if( fabs(longitude) > 180.0 )
    longitude = -999.9;
  if( fabs(latitude) > 90.0 )
    latitude = -999.9;
  
  if( fabs(latitude)>999.0 || fabs(longitude)>999.0 )
    return;
 

  AuxWindow *window = AuxWindow::make( "Google Map" );
  
  const int w = static_cast<int>(0.66*m_specViewer->renderedWidth());
  const int h = static_cast<int>(0.8*m_specViewer->renderedHeight());
  
 
  window->disableCollapse();
  window->finished().connect( [window](){ AuxWindow::deleteAuxWindow( window ); } );
 // window->footer()->setStyleClass( "modal-footer" );
  window->footer()->setHeight(WLength(50,WLength::Unit::Pixel));
  WPushButton *closeButton = window->addCloseButtonToFooter();
  closeButton->clicked().connect( window, &AuxWindow::hide );
  

  WGridLayout *layout = window->stretcher();
  GoogleMap *googlemap = layout->addWidget( std::make_unique<GoogleMap>( false ), 0, 0 );
  googlemap->addMarker( latitude, longitude );
  googlemap->adjustPanAndZoom();
  googlemap->mapClicked().connect( [googlemap, latitude, longitude]( double lat, double lng ){
    updateCoordText( googlemap, lat, lng, latitude, longitude );
  } );
  window->resize( WLength(w), WLength(h) );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->setVerticalSpacing( 0 );
  layout->setHorizontalSpacing( 0 );

  window->show();
  window->centerWindow();
  window->rejectWhenEscapePressed();
  //  window->resizeToFitOnScreen();
  
//  char buffer[128];
//  snprintf( buffer, sizeof(buffer), "(%.6f, %.6f)", latitude, longitude );
//  coordTxt->setText( buffer );
}//void showGoogleMap();
#elif( USE_LEAFLET_MAP ) //#if( USE_GOOGLE_MAP )
void SpecFileSummary::showLeafletMap()
{
  if( !m_specViewer )
    return;
  
  double longitude = coordVal( m_longitude->text() );
  double latitude = coordVal( m_latitude->text() );
  if( fabs(longitude) > 180.0 )
    longitude = -999.9;
  if( fabs(latitude) > 90.0 )
    latitude = -999.9;
  
  if( fabs(latitude)>999.0 || fabs(longitude)>999.0 )
    return;
  
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( type );
  std::shared_ptr<const SpecUtils::Measurement> sample = currentMeasurment();
  
  if( !meas || !sample )
    return;
  
  
  const set<int> sample_numbers{ sample->sample_number() };
  const vector<string> det_names( 1, sample->detector_name() );
  
  const AllowModifyStatus status = AllowModifyStatus(m_allowEditGroup->checkedId());
  if( status == AllowModifyStatus::kDontAllowModify )
  {
    LeafletRadMap::showForMeasurement(meas, sample_numbers, det_names);
    return;
  }
  
  
  try
  {
    shared_ptr<SpecMeas> new_meas = make_shared<SpecMeas>();
    new_meas->uniqueCopyContents( *meas );
    
    auto new_sample = meas->measurement( sample->sample_number(), sample->detector_name() );
    if( !new_sample )
      throw runtime_error( "failed to find measurement from copy of SpecFile." );
    
    new_meas->set_position(longitude, latitude, sample->position_time(), new_sample );
    
    LeafletRadMap::showForMeasurement(new_meas, sample_numbers, det_names);
  }catch( std::exception &e )
  {
    passMessage( WString::tr("sfs-err-show-map").arg(e.what()), WarningWidget::WarningMsgHigh );
  }
}//void showLeafletMap()
#endif


void SpecFileSummary::updateMeasurmentFieldsFromMemory()
{
  char buffer[64];
  try
  {
    std::shared_ptr<const SpecUtils::Measurement> sample = currentMeasurment();

    if( !sample )
      throw runtime_error("");

    const float liveTime = float(sample->live_time() / PhysicalUnits::second);
    const float realTime = float(sample->real_time() / PhysicalUnits::second);

    const double sumGamma = sample->gamma_count_sum();
    const double gammaCps = (liveTime > FLT_EPSILON) ? sumGamma/liveTime : -1.0;
    
    snprintf( buffer, sizeof(buffer), "%.3f", gammaCps );
    if( gammaCps < -FLT_EPSILON )
      m_gammaCPS->setText( "--" );
    else
      m_gammaCPS->setText( buffer );
    
    snprintf( buffer, sizeof(buffer), "%.2f", sumGamma );
    m_gammaSum->setText( buffer );

    if( sample->contained_neutron() )
    {
      const double sumNeutron = sample->neutron_counts_sum();
      const double neutronCps = (realTime > FLT_EPSILON) ? sumNeutron/realTime : -1.0;

      snprintf( buffer, sizeof(buffer), "%.3f", neutronCps );
      if( neutronCps < -FLT_EPSILON )
        m_neutronCPS->setText( "--" );
      else
        m_neutronCPS->setText( buffer );
      
      snprintf( buffer, sizeof(buffer), "%.2f", sumNeutron );
      m_neutronSum->setText( buffer );
    }else
    {
      m_neutronCPS->setText( "N/A" );
      m_neutronSum->setText( "N/A" );
    }

    snprintf( buffer, sizeof(buffer), "%.2f s", liveTime );
    m_displayedLiveTime->setText( buffer );
    snprintf( buffer, sizeof(buffer), "%.2f s", realTime );
    m_displayedRealTime->setText( buffer );

    if( fabs( sample->latitude() ) < 999.0 )
    {
      snprintf( buffer, sizeof(buffer), "%.6f", sample->latitude() );
      m_latitude->setText( buffer );
    }else
      m_latitude->setText( "" );
    if( fabs(sample->longitude()) < 999.0 )
    {
      snprintf( buffer, sizeof(buffer), "%.6f", sample->longitude() );
      m_longitude->setText( buffer );
    }else
      m_longitude->setText( "" );

#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
    bool nomap = (fabs(sample->latitude())>999.0 || fabs(sample->longitude())>999.0);
    m_showMapButton->setDisabled( nomap );
#endif
    
    if( SpecUtils::is_special(sample->position_time()) )
    {
      m_gpsTimeStamp->setText( "" );
    }else
    {
      m_gpsTimeStamp->setText( SpecUtils::to_extended_iso_string(sample->position_time()) );
    }
    
    if( SpecUtils::is_special(sample->start_time()) )
    {
      m_timeStamp->setText( "" );
    }else
    {
      m_timeStamp->setText( SpecUtils::to_extended_iso_string(sample->start_time()) );
    }

    const size_t nchannel = sample->num_gamma_channels();
    m_numberOfBins->setText( std::to_string(nchannel) );
    
    const float min_energy = sample->gamma_energy_min();
    const float max_energy = sample->gamma_energy_max();
    
    if( min_energy != max_energy )
    {
      snprintf( buffer, sizeof(buffer), "%.1f to %.1f", min_energy, max_energy );
      m_energyRange->setText( buffer );
    }else
    {
      m_energyRange->setText( "" );
    }//if( defined energy range ) / else

    m_detector->setText( sample->detector_name() );

    const int specNum = sample->sample_number();
    m_sampleNumber->setText( std::to_string(specNum) );

    string remark;
    const vector<string> &remarks = sample->remarks();
    for( size_t i = 0; i < remarks.size(); ++i )
    {
      if( i )
        remark += '\n';
      remark += remarks[i];
    }//for( size_t i = 0; i < remarks.size(); ++i )

    m_title->setText( WString::fromUTF8(sample->title()) );
    m_source->setCurrentIndex( static_cast<int>(sample->source_type()) );
    
    m_measurmentRemarks->setText( remark );
  }catch(...)
  {
    m_gammaCPS->setText( "-" );
    m_gammaSum->setText( "-" );
    m_neutronCPS->setText( "-" );
    
    m_longitude->setText( "" );
    m_latitude->setText( "" );
    m_gpsTimeStamp->setText( "" );
#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
    m_showMapButton->disable();
#endif
    
    m_neutronSum->setText( "-" );
    m_displayedLiveTime->setText( "" );
    m_displayedRealTime->setText( "" );
    m_timeStamp->setText( "" );

    m_title->setText( "" );
    m_source->setCurrentIndex( -1 );

    m_energyRange->setText( "-" );
    m_numberOfBins->setText( "-" );

    m_detector->setText( "-" );
    m_sampleNumber->setText( "-" );
    m_measurmentRemarks->setText( "-" );
  }//try / catch
}//void updateMeasurmentFieldsFromMemory()


void SpecFileSummary::updateDisplayFromMemory()
{
//  return;
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());

  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( type );
  const size_t nspec = (!!meas ? meas->measurements().size() : 0);

  m_displaySampleDiv->setHidden( nspec < 2 );
  const bool hasForeground = !!m_specViewer->measurment(SpecUtils::SpectrumType::Foreground);
  const bool hasSecondFore = !!m_specViewer->measurment(SpecUtils::SpectrumType::SecondForeground);
  const bool hasBackground = !!m_specViewer->measurment(SpecUtils::SpectrumType::Background);
  m_spectraGroup->button(static_cast<int>(SpecUtils::SpectrumType::Foreground))->setEnabled( hasForeground );
  m_spectraGroup->button(static_cast<int>(SpecUtils::SpectrumType::SecondForeground))->setEnabled( hasSecondFore );
  m_spectraGroup->button(static_cast<int>(SpecUtils::SpectrumType::Background))->setEnabled( hasBackground );
  
  
  m_displayedPostText->setText( "of " + std::to_string(nspec) );
  if( nspec )
  {
    m_displaySampleNumValidator->setRange( 1, static_cast<int>(nspec) );
    
    // Lets set the displayed measurement to be the first one that is actually displayed in the
    //  spectrum
    size_t sample = 0;
    const vector<shared_ptr<const SpecUtils::Measurement>> &measurements = meas->measurements();
    const set<int> &samples = m_specViewer->displayedSamples(type);
    if( samples.size() )
    {
      for( size_t i = 0; i < measurements.size(); ++i )
      {
        if( samples.count(measurements[i]->sample_number()) )
        {
          sample = i;
          break;
        }
      }//for( loop over measurements )
    }//if( samples.size() )
    
    m_displaySampleNumEdit->setText( std::to_string(sample + 1) );
  }else
  {
    m_displaySampleNumValidator->setRange( 0, 0 );
    m_displaySampleNumEdit->setText( "" );
  }//if( nspec ) / else

  if( !!meas )
  {
    double memsize = static_cast<double>( meas->memmorysize() );
    const char *memunit = "b";

    if( memsize > 1024*1024 )
    {
      memunit = "Mb";
      memsize /= (1024.0*1024.0);
    }else if( memsize > 1024 )
    {
      memunit = "kb";
      memsize /= 1024.0;
    }//

    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.1f %s", memsize, memunit );
    m_sizeInMemmory->setText( buffer );
    m_filename->setText( meas->filename() );
    m_uuid->setText( meas->uuid() );

    try
    {
      if( meas->lane_number() < 0 )
        throw std::runtime_error( "" );
      m_laneNumber->setText( std::to_string( meas->lane_number() ) );
    }catch(...)
    {
      m_laneNumber->setText( "" );
    }

    m_measurement_location_name->setText( meas->measurement_location_name() );
    
    if( meas && meas->detectors_analysis() && !m_ana_button->isEnabled() )
    {
      m_ana_button->setText( WString::tr("sfs-show-det-id-btn") );
      m_ana_button->enable();
    }else if( m_ana_button->isEnabled() )
    {
      m_ana_button->setText( WString::tr("sfs-no-det-id-btn") );
      m_ana_button->disable();
    }
    
    m_multimedia_button->setHidden( !meas || meas->multimedia_data().empty() );
    
    m_inspection->setText( meas->inspection() );
    m_instrument_type->setText( meas->instrument_type() );
    m_manufacturer->setText( meas->manufacturer() );
    m_instrument_model->setText( meas->instrument_model() );
    m_instrument_id->setText( meas->instrument_id() );
    
    string remark;
    const vector<string> &remarks = meas->remarks();
    for( size_t i = 0; i < remarks.size(); ++i )
    {
      if( i )
        remark += '\n';
      remark += remarks[i];
    }//for( size_t i = 0; i < remarks.size(); ++i )

    m_fileRemarks->setText( WString::fromUTF8(remark) );
  }else
  {
    m_sizeInMemmory->setText( "" );
    m_filename->setText( "" );
    m_uuid->setText( "" );
    m_laneNumber->setText( "" );
    m_measurement_location_name->setText( "" );
    m_ana_button->setText( "" );
    m_ana_button->disable();
    m_inspection->setText( "" );
    m_instrument_type->setText( "" );
    m_manufacturer->setText( "" );
    m_instrument_model->setText( "" );
    m_instrument_id->setText( "" );
    m_fileRemarks->setText( "" );
  }//if( meas ) / else
  
  updateMeasurmentFieldsFromMemory();
}//void SpecFileSummary::updateDisplayFromMemory()

std::shared_ptr<const SpecUtils::Measurement> SpecFileSummary::currentMeasurment() const
{
  std::shared_ptr<const SpecUtils::Measurement> sample;

  try
  {
    const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
    std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( type );
    if( !meas )
      throw runtime_error( "" );

    const vector<std::shared_ptr<const SpecUtils::Measurement>> &measurements = meas->measurements();
    const string sampleNumStr = m_displaySampleNumEdit->text().toUTF8();
    size_t sampleNum = 0;
    if( m_displaySampleDiv->isEnabled() && measurements.size()>1 )
    {
      sampleNum = std::stoull( sampleNumStr );
      if( sampleNum > 0 )
        sampleNum -= 1;
    }

    if( sampleNum > measurements.size() )
      throw runtime_error("");

    sample = measurements.at( sampleNum );
  }catch(...)
  {
    sample.reset();  //not necessary
  }

  return sample;
}//currentMeasurment() const


void SpecFileSummary::handleFieldUpdate( EditableFields field )
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  std::shared_ptr<SpecMeas> meas = m_specViewer->measurment( type );
  std::shared_ptr<const SpecUtils::Measurement> sample = currentMeasurment();

  if( !sample || !meas )
  {
    passMessage( WString::tr("sfs-no-spectrum-to-update"), WarningWidget::WarningMsgInfo );
    return;
  }//if( !sample )

  m_reloadSpectrum->enable();

  switch( field )
  {
    case kDisplayedLiveTime:
    {
      try
      {
        const string textStr = m_displayedLiveTime->text().toUTF8();
        const float newLiveTime = static_cast<float>(PhysicalUnitsLocalized::stringToTimeDuration( textStr ));
        if( newLiveTime < 0.0 )
          throw runtime_error( "Live time must be zero or greater" );
        
        const float oldLiveTime = sample->live_time();
        if( newLiveTime == oldLiveTime )
          return;
        meas->set_live_time( newLiveTime, sample );
      }catch( exception &e )
      {
        passMessage( e.what(), WarningWidget::WarningMsgHigh );
        
        char text[32];
        snprintf( text, sizeof(text), "%.2f s", sample->live_time() );
        m_displayedLiveTime->setText( text );
      }//try / catch
      break;
    }//case kDisplayedLiveTime:

    case kDisplayedRealTime:
    {
      try
      {
        const string textStr = m_displayedRealTime->text().toUTF8();
        const float newRealTime = static_cast<float>( PhysicalUnitsLocalized::stringToTimeDuration(textStr) );
        if( newRealTime < 0.0 )
          throw runtime_error( "Real time must be zero or greater" );
        
        const float oldRealTime = sample->real_time();
        if( newRealTime == oldRealTime )
          return;
        meas->set_real_time( newRealTime, sample );
      }catch( exception &e )
      {
        passMessage( e.what(), WarningWidget::WarningMsgHigh );
        
        char text[32];
        snprintf( text, sizeof(text), "%.2f", sample->real_time() );
        m_displayedRealTime->setText( text );
      }//try / catch
      break;
    }//case kDisplayedRealTime:

    case kTimeStamp:
    {
      const string newDateStr = m_timeStamp->text().toUTF8();
      SpecUtils::time_point_t newDate = SpecUtils::time_from_string( newDateStr.c_str() );

      if( SpecUtils::is_special(newDate) )
      {
        passMessage( WString::tr("sfs-err-time-format").arg(newDateStr),
                    WarningWidget::WarningMsgHigh );
        m_timeStamp->setText( SpecUtils::to_extended_iso_string( sample->start_time() ) );
        return;
      }//if( newDate.is_special() )

      meas->set_start_time( newDate, sample );
      break;
    }//case kTimeStamp:

    case kPosition:
    {
      double longitude = coordVal( m_longitude->text() );
      double latitude = coordVal( m_latitude->text() );
      
//      cerr << "Long '" << m_longitude->text().toUTF8() << "'--->" << longitude << endl;
//      cerr << "Lat  '" << m_latitude->text().toUTF8() << "'--->" << latitude << endl;
      
      if( fabs(longitude) > 180.0 )
        longitude = -999.9;
      if( fabs(latitude) > 90.0 )
        latitude = -999.9;
      
      bool converror = false;
      if( !m_longitude->text().empty() )
      {
        if( fabs(longitude)>999.0 && fabs(sample->longitude())<999.0 )
        {
          converror = true;
          longitude = sample->longitude();
          char buffer[32];
          snprintf( buffer, sizeof(buffer), "%.6f", longitude );
          m_longitude->setText( buffer );
        }else if( fabs(longitude)>999.0 )
        {
          converror = true;
          m_longitude->setText( "" );
        }
      }//if( !m_longitude->text().empty() )
      
      if( !m_latitude->text().empty() )
      {
        if( fabs(latitude)>999.0 && fabs(sample->latitude())<999.0 )
        {
          converror = true;
          latitude = sample->latitude();
          char buffer[32];
          snprintf( buffer, sizeof(buffer), "%.6f", latitude );
          m_latitude->setText( buffer );
        }else if( fabs(latitude)>999.0 )
        {
          converror = true;
          m_latitude->setText( "" );
        }
      }//if( !m_latitude->text().empty() )
      
      if( converror )
        passMessage( WString::tr("sfs-err-lat-lon"), WarningWidget::WarningMsgHigh );
      
      
      const string newDateStr = m_gpsTimeStamp->text().toUTF8();
      SpecUtils::time_point_t newDate = SpecUtils::time_from_string( newDateStr.c_str() );
      
      if( SpecUtils::is_special(newDate) && newDateStr.size() )
      {
        passMessage( WString::tr("sfs-err-time-format").arg(newDateStr),
                     WarningWidget::WarningMsgHigh );
        newDate = sample->position_time();
        m_timeStamp->setText( SpecUtils::to_extended_iso_string(newDate) );
      }//if( newDate.is_special() )

#if( USE_GOOGLE_MAP || USE_LEAFLET_MAP )
      if( fabs(latitude)<999.0 && fabs(longitude)<999.0 )
        m_showMapButton->enable();
#endif
      
      meas->set_position( longitude, latitude, newDate, sample );
      
      break;
    }//case kPosition:
      
    case kDescription:
    {
      const string newtitle = m_title->text().toUTF8();
      meas->set_title( newtitle, sample );
      break;
    }//case kDescription:
      
    case kSourceType:
    {
      const int index = m_source->currentIndex();
      SpecUtils::SourceType type = SpecUtils::SourceType::Unknown;
      if( index >= 0 )
        type = SpecUtils::SourceType( m_source->currentIndex() );
      
      meas->set_source_type( type, sample );
      break;
    }//case kSourceType:
      
    case kMeasurmentRemarks:
    {
      const string newRemark = m_measurmentRemarks->text().toUTF8();
      vector<string> newRemarks;
      boost::algorithm::split( newRemarks, newRemark,
                               boost::is_any_of( "\r\n" ),
                               boost::token_compress_on );
      meas->set_remarks( newRemarks, sample );
      break;
    }//case kMeasurmentRemarks:


    case kFileRemarks:
    {
      const string newRemark = m_fileRemarks->text().toUTF8();
      vector<string> newRemarks;
      boost::algorithm::split( newRemarks, newRemark,
                               boost::is_any_of( "\r\n" ),
                               boost::token_compress_on );
      meas->set_remarks( newRemarks );
      break;
    }//case kFileRemarks:

    case kFilename:
      meas->set_filename( m_filename->text().toUTF8() );
      break;

    case kUuid:
      meas->set_uuid( m_uuid->text().toUTF8() );
      break;

    case kLaneNumber:
    {
      try
      {
        meas->set_lane_number( std::stoi( m_laneNumber->text().toUTF8() ) );
      }catch(...)
      {
        m_laneNumber->setText( std::to_string( meas->lane_number() ) );
      }
      break;
    }//case kLaneNumber:

    case kMeasurement_location_name:
      meas->set_measurement_location_name( m_measurement_location_name->text().toUTF8() );
      break;

    case kInspection:
      meas->set_inspection( m_inspection->text().toUTF8() );
      break;

    case kInstrument_type:
      meas->set_instrument_type( m_instrument_type->text().toUTF8() );
      break;

    case kManufacturer:
      meas->set_manufacturer( m_manufacturer->text().toUTF8() );
      break;

    case kInstrument_model:
      meas->set_instrument_model( m_instrument_model->text().toUTF8() );
      break;

    case kInstrument_id:
      meas->set_instrument_id( m_instrument_id->text().toUTF8() );
      break;
  }//switch( field )
}//void handleFieldUpdate( EditableFields field )


void SpecFileSummary::handleSpectrumChange( const SpecUtils::SpectrumType type,
                                            const std::shared_ptr<SpecMeas> &meas,
                                            const std::set<int> &displaySample,
                                            const std::vector<std::string> &detectors )
{
  const SpecUtils::SpectrumType display_type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  
  if( type != display_type )
  {
    WRadioButton *button = m_spectraGroup->button( static_cast<int>(type) );
    if( button )
      button->setEnabled( !!meas );
  }else
  {
    m_reloadSpectrum->disable();
    updateDisplayFromMemory();
  }
}//void handleSpectrumChange(...)


void SpecFileSummary::reloadCurrentSpectrum()
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());
  m_specViewer->reloadCurrentSpectrum( type );
  m_reloadSpectrum->disable();
}//void reloadCurrentSpectrum()


void SpecFileSummary::handleSpectrumTypeChanged()
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());

  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( type );
  if( !meas )
  {
    m_displaySampleNumEdit->setText( "" );
    m_displaySampleNumValidator->setRange( 0, 0 );
  }else
  {
    const vector<std::shared_ptr<const SpecUtils::Measurement>> &measurements = meas->measurements();
    
    if( measurements.size() )
    {
      m_displaySampleNumValidator->setRange( 1, static_cast<int>(measurements.size()) );
      m_displaySampleNumEdit->setText( "1" );
    }else
    {
      m_displaySampleNumValidator->setRange( 0, 0 );
      m_displaySampleNumEdit->setText( "" );
    }
  }//if( !meas ) /else

  m_reloadSpectrum->disable();

  updateDisplayFromMemory();
}//void handleSpectrumTypeChanged( Wt::WRadioButton *button )



void SpecFileSummary::handleAllowModifyStatusChange()
{
  const AllowModifyStatus status = AllowModifyStatus(m_allowEditGroup->checkedId());

  bool disable = true;
  switch( status )
  {
    case kAllowModify:
      disable = false;
    break;

    case kDontAllowModify:
      disable = true;
    break;
  }//switch( status )

  m_displayedLiveTime->setDisabled( disable );
  m_displayedRealTime->setDisabled( disable );
  m_timeStamp->setDisabled( disable );
  m_measurmentRemarks->setDisabled( disable );
  m_fileRemarks->setDisabled( disable );

  m_longitude->setDisabled( disable );
  m_latitude->setDisabled( disable );
  m_gpsTimeStamp->setDisabled( disable );
  
  m_title->setDisabled( disable );
  m_source->setDisabled( disable );
  
  m_filename->setDisabled( disable );
  m_uuid->setDisabled( disable );
  m_laneNumber->setDisabled( disable );
  m_measurement_location_name->setDisabled( disable );
  m_inspection->setDisabled( disable );
  m_instrument_type->setDisabled( disable );
  m_manufacturer->setDisabled( disable );
  m_instrument_model->setDisabled( disable );
  m_instrument_id->setDisabled( disable );
}//void handleAllowModifyStatusChange(...)


void SpecFileSummary::handleUserIncrementSampleNum( bool increment )
{
  const SpecUtils::SpectrumType type = SpecUtils::SpectrumType(m_spectraGroup->checkedId());

  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( type );

  try
  {
    if( !meas )
    {
      m_displaySampleNumValidator->setRange( 0, 0 );
      throw runtime_error("");
    }//if( !meas )

    const vector<std::shared_ptr<const SpecUtils::Measurement>> &measurements = meas->measurements();
    if( measurements.size() )
    {
      m_displaySampleNumValidator->setRange( 1, static_cast<int>(measurements.size()) );
    }else
    {
      m_displaySampleNumValidator->setRange( 0, 0 );
      throw runtime_error( "" );
    }

    size_t prevSpecNum;
    try
    {
      prevSpecNum = std::stoull( m_displaySampleNumEdit->text().toUTF8() );
    }catch(...)
    {
      prevSpecNum = 0;
    }

    if( prevSpecNum > 0 )
      prevSpecNum -= 1;

    size_t nextSpecNum;

    if( increment )
    {
      if( prevSpecNum >= (measurements.size()-1) )
        nextSpecNum = 0;
      else
        nextSpecNum = prevSpecNum + 1;
    }else
    {
      if( prevSpecNum == 0 )
        nextSpecNum = measurements.size()-1;
      else
        nextSpecNum = prevSpecNum - 1;
    }//if( increment ) / else

    m_displaySampleNumEdit->setText( std::to_string(nextSpecNum+1) );
  }catch(...)
  {
    m_displaySampleNumEdit->setText( "" );
  }

  updateMeasurmentFieldsFromMemory();
}//void SpecFileSummary::handleUserIncrementSampleNum( bool increment )


void SpecFileSummary::handleUserChangeSampleNum()
{
  updateMeasurmentFieldsFromMemory();
}//void SpecFileSummary::handleUserChangeSampleNum()

