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

#include <string>
#include <vector>
#include <limits>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "Minuit2/MnUserParameters.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/SimpleActivityCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"


#if( USE_QR_CODES )
#include <Wt/Utils>
#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;

SimpleActivityCalcState::SimpleActivityCalcState()
  : peakEnergy( -1 ),
    nuclideName( "" ),
    nuclideAgeStr( "" ),
    distanceStr( "" ),
    geometryType( SimpleActivityGeometryType::Point ),
    shielding{} // Empty optional
{
}

bool SimpleActivityCalcState::operator==( const SimpleActivityCalcState &rhs ) const
{
  if( (peakEnergy != rhs.peakEnergy)
     || (nuclideName != rhs.nuclideName)
     || (nuclideAgeStr != rhs.nuclideAgeStr)
     || (distanceStr != rhs.distanceStr)
     || (geometryType != rhs.geometryType) )
  {
    return false;
  }
  
  // Handle optional shielding comparison
  if( shielding.has_value() != rhs.shielding.has_value() )
    return false;
    
  if( shielding.has_value() )
  {
    try
    {
#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
      ShieldingSourceFitCalc::ShieldingInfo::equalEnough( *shielding, *rhs.shielding );
      return true;
#else
      // For non-test builds, we'll do a simpler comparison since equalEnough throws exceptions
      return (shielding->m_geometry == rhs.shielding->m_geometry)
             && (shielding->m_isGenericMaterial == rhs.shielding->m_isGenericMaterial)
             && (shielding->m_forFitting == rhs.shielding->m_forFitting)
             && (std::abs(shielding->m_dimensions[0] - rhs.shielding->m_dimensions[0]) < 1e-9)
             && (std::abs(shielding->m_dimensions[1] - rhs.shielding->m_dimensions[1]) < 1e-9)
             && (std::abs(shielding->m_dimensions[2] - rhs.shielding->m_dimensions[2]) < 1e-9);
#endif
    }
    catch( const std::exception & )
    {
      return false;
    }
  }
  
  return true; // Both empty optionals
}

bool SimpleActivityCalcState::operator!=( const SimpleActivityCalcState &rhs ) const
{
  return !(*this == rhs);
}

std::string SimpleActivityCalcState::encodeToUrl() const
{
  return "dummy-url-implementation";
}

void SimpleActivityCalcState::decodeFromUrl( const std::string &uri )
{
}

void SimpleActivityCalcState::serialize( rapidxml::xml_node<char> * const parent_node ) const
{
  using namespace rapidxml;
  
  assert( parent_node );
  if( !parent_node || !parent_node->document() )
    throw std::runtime_error( "SimpleActivityCalcState::serialize: invalid parent node" );
    
  xml_document<char> *doc = parent_node->document();
  xml_node<char> *base_node = doc->allocate_node( node_element, "SimpleActivityCalcState" );
  parent_node->append_node( base_node );
  
  XmlUtils::append_version_attrib( base_node, SimpleActivityCalcState::sm_xmlSerializationVersion );
  
  // Serialize basic fields
  XmlUtils::append_float_node( base_node, "PeakEnergy", peakEnergy );
  XmlUtils::append_string_node( base_node, "NuclideName", nuclideName );
  XmlUtils::append_string_node( base_node, "NuclideAgeStr", nuclideAgeStr );
  XmlUtils::append_string_node( base_node, "DistanceStr", distanceStr );
  
  // Serialize geometry type as integer
  XmlUtils::append_int_node( base_node, "GeometryType", static_cast<int>(geometryType) );
  
  // Serialize shielding info only if present
  if( shielding.has_value() )
  {
    shielding->serialize( base_node );
  }
}

void SimpleActivityCalcState::deSerialize( const ::rapidxml::xml_node<char> *src_node, MaterialDB *materialDB )
{
  using namespace rapidxml;
  
  try
  {
    if( !src_node )
      throw std::runtime_error( "SimpleActivityCalcState::deSerialize: invalid input node" );
    
    if( !rapidxml::internal::compare( src_node->name(), src_node->name_size(), "SimpleActivityCalcState", 23, false ) )
      throw std::logic_error( "SimpleActivityCalcState::deSerialize: invalid input node name" );
    
    // Check version compatibility
    static_assert( SimpleActivityCalcState::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    XmlUtils::check_xml_version( src_node, SimpleActivityCalcState::sm_xmlSerializationVersion );
    
    // Deserialize basic fields
    peakEnergy = XmlUtils::get_float_node_value( src_node, "PeakEnergy" );
    
    const rapidxml::xml_node<char> *nuclide_name_node = XML_FIRST_NODE( src_node, "NuclideName" );
    if( nuclide_name_node )
      nuclideName = SpecUtils::xml_value_str( nuclide_name_node );
    
    const rapidxml::xml_node<char> *age_str_node = XML_FIRST_NODE( src_node, "NuclideAgeStr" );
    if( age_str_node )
      nuclideAgeStr = SpecUtils::xml_value_str( age_str_node );
    
    const rapidxml::xml_node<char> *distance_str_node = XML_FIRST_NODE( src_node, "DistanceStr" );
    if( distance_str_node )
      distanceStr = SpecUtils::xml_value_str( distance_str_node );
    
    // Deserialize geometry type from integer
    const rapidxml::xml_node<char> *geom_type_node = XML_FIRST_NODE( src_node, "GeometryType" );
    int geomTypeInt = 0; // Default to Point
    if( geom_type_node )
    {
      const std::string geom_str = SpecUtils::xml_value_str( geom_type_node );
      geomTypeInt = std::stoi( geom_str );
    }
    geometryType = static_cast<SimpleActivityGeometryType>(geomTypeInt);
    
    // Validate geometry type value
    if( geomTypeInt < 0 || geomTypeInt > static_cast<int>(SimpleActivityGeometryType::TraceSrc) )
      throw std::runtime_error( "SimpleActivityCalcState::deSerialize: invalid geometry type value" );
    
    // Deserialize shielding info - find the ShieldingInfo node
    const rapidxml::xml_node<char> *shielding_node = XML_FIRST_NODE( src_node, "Shielding" );
    if( shielding_node )
    {
      shielding = ShieldingSourceFitCalc::ShieldingInfo();
      shielding->deSerialize( shielding_node, materialDB );
    }
    else
    {
      // If no shielding node found, leave as empty optional
      shielding.reset();
    }
  }
  catch( const std::exception &e )
  {
    // Reset to default state on any error
    *this = SimpleActivityCalcState();
    throw std::runtime_error( std::string("SimpleActivityCalcState::deSerialize: ") + e.what() );
  }
}

#if( PERFORM_DEVELOPER_CHECKS || BUILD_AS_UNIT_TEST_SUITE )
void SimpleActivityCalcState::equalEnough( const SimpleActivityCalcState &lhs, const SimpleActivityCalcState &rhs )
{
  // Compare peak energy
  if( std::abs(lhs.peakEnergy - rhs.peakEnergy) > 1e-6 )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: peakEnergy differs - LHS=" 
                             + std::to_string(lhs.peakEnergy) + ", RHS=" + std::to_string(rhs.peakEnergy) );
  }
  
  // Compare nuclide name
  if( lhs.nuclideName != rhs.nuclideName )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: nuclideName differs - LHS='" 
                             + lhs.nuclideName + "', RHS='" + rhs.nuclideName + "'" );
  }
  
  // Compare nuclide age string
  if( lhs.nuclideAgeStr != rhs.nuclideAgeStr )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: nuclideAgeStr differs - LHS='" 
                             + lhs.nuclideAgeStr + "', RHS='" + rhs.nuclideAgeStr + "'" );
  }
  
  // Compare distance string
  if( lhs.distanceStr != rhs.distanceStr )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: distanceStr differs - LHS='" 
                             + lhs.distanceStr + "', RHS='" + rhs.distanceStr + "'" );
  }
  
  // Compare geometry type
  if( lhs.geometryType != rhs.geometryType )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: geometryType differs - LHS=" 
                             + std::to_string(static_cast<int>(lhs.geometryType)) + ", RHS=" + std::to_string(static_cast<int>(rhs.geometryType)) );
  }
  
  // Compare optional shielding
  if( lhs.shielding.has_value() != rhs.shielding.has_value() )
  {
    throw std::runtime_error( "SimpleActivityCalcState::equalEnough: shielding presence differs - LHS has_value=" 
                             + std::string(lhs.shielding.has_value() ? "true" : "false") + ", RHS has_value=" + std::string(rhs.shielding.has_value() ? "true" : "false") );
  }
  
  // If both have shielding, compare using ShieldingInfo::equalEnough
  if( lhs.shielding.has_value() && rhs.shielding.has_value() )
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo::equalEnough( *lhs.shielding, *rhs.shielding );
    }
    catch( const std::exception &e )
    {
      throw std::runtime_error( "SimpleActivityCalcState::equalEnough: shielding differs - " + std::string(e.what()) );
    }
  }
}//void SimpleActivityCalcState::equalEnough( const SimpleActivityCalcState &lhs, const SimpleActivityCalcState &rhs )
#endif

SimpleActivityCalcWindow::SimpleActivityCalcWindow( MaterialDB *materialDB,
                                                  Wt::WSuggestionPopup *materialSuggestion,
                                                  InterSpec* viewer )
: AuxWindow( WString::tr("simple-activity-calc-title"),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse) ),
  m_tool( nullptr )
{
  setModal( false );

  m_tool = new SimpleActivityCalc( materialDB, materialSuggestion, viewer, contents() );
  
  AuxWindow::addHelpInFooter( footer(), "simple-activity-calc" );
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton();
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://quick-activity/" + Wt::Utils::urlEncode(m_tool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("sac-qr-tool-state-title"),
                                 WString::tr("sac-qr-tool-state-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
  if( !viewer->isPhone() )
    footer()->addWidget( qr_btn );
#endif //USE_QR_CODES
  
  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
#if( USE_QR_CODES )
  if( viewer->isPhone() )
    footer()->addWidget( qr_btn );
#endif
  
  centerWindow();
  resizeToFitOnScreen();
  show();
}

SimpleActivityCalcWindow::~SimpleActivityCalcWindow()
{
}

SimpleActivityCalc *SimpleActivityCalcWindow::tool()
{
  return m_tool;
}

SimpleActivityCalc::SimpleActivityCalc( MaterialDB *materialDB,
                                      Wt::WSuggestionPopup *materialSuggestion,
                                      InterSpec *specViewer,
                                      Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_renderFlags( 0x0 ),
  m_viewer( specViewer ),
  m_materialSuggest( materialSuggestion ),
  m_materialDB( materialDB ),
  m_peakSelect( nullptr ),
  m_nuclideInfo( nullptr ),
  m_ageEdit( nullptr ),
  m_distanceEdit( nullptr ),
  m_detectorDisplay( nullptr ),
  m_shieldingSelect( nullptr ),
  m_geometrySelect( nullptr ),
  m_resultText( nullptr ),
  m_errorText( nullptr ),
  m_advancedBtn( nullptr ),
  m_distanceRow( nullptr ),
  m_ageRow( nullptr ),
  m_geometryRow( nullptr ),
  m_currentPeak( nullptr )
{
  if( !m_materialDB )
    throw logic_error( "SimpleActivityCalc requires a valid MaterialDB." );
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance() );
  if( app )
  {
    app->useMessageResourceBundle( "SimpleActivityCalc" );
    app->useStyleSheet( "InterSpec_resources/SimpleActivityCalc.css" );
  }

  addStyleClass( "SimpleActivityCalc" );
  
  init();
}

SimpleActivityCalc::~SimpleActivityCalc()
{
}

void SimpleActivityCalc::init()
{
  // Create the main layout using CSS flexbox
  setLayoutSizeAware( true );
  
  // Peak selection row
  WContainerWidget *peakRow = new WContainerWidget( this );
  peakRow->addStyleClass( "row" );
  WText *peakLabel = new WText( WString::tr("sac-peak-label"), peakRow );
  peakLabel->addStyleClass( "label" );
  m_peakSelect = new WComboBox( peakRow );
  m_peakSelect->addStyleClass( "input" );
  m_peakSelect->changed().connect( this, &SimpleActivityCalc::handlePeakChanged );
  
  // Nuclide info display (shows nuclide, energy, branching ratio, etc.)
  m_nuclideInfo = new WText( this );
  m_nuclideInfo->addStyleClass( "row nuclide-info" );
  
  // Age row (hidden if not appropriate for the nuclide)
  m_ageRow = new WContainerWidget( this );
  m_ageRow->addStyleClass( "row" );
  WText *ageLabel = new WText( WString::tr("sac-age-label"), m_ageRow );
  ageLabel->addStyleClass( "label" );
  m_ageEdit = new WLineEdit( m_ageRow );
  m_ageEdit->addStyleClass( "input" );
  m_ageEdit->changed().connect( this, &SimpleActivityCalc::handleAgeChanged );
  m_ageRow->hide(); // Initially hidden
  
  // Distance row (hidden if fixed geometry)
  m_distanceRow = new WContainerWidget( this );
  m_distanceRow->addStyleClass( "row" );
  WText *distanceLabel = new WText( WString::tr("sac-distance-label"), m_distanceRow );
  distanceLabel->addStyleClass( "label" );
  m_distanceEdit = new WLineEdit( m_distanceRow );
  m_distanceEdit->addStyleClass( "input" );
  m_distanceEdit->setText( "1 m" );
  m_distanceEdit->setValidator( new WRegExpValidator( PhysicalUnits::sm_distanceRegex ) );
  m_distanceEdit->changed().connect( this, &SimpleActivityCalc::handleDistanceChanged );
  
  // Detector row
  WContainerWidget *detectorRow = new WContainerWidget( this );
  detectorRow->addStyleClass( "row" );
  m_detectorDisplay = new DetectorDisplay( m_viewer, m_viewer->fileManager()->model(), detectorRow );
  m_detectorDisplay->addStyleClass( "input" );
  
  // Shielding row
  WContainerWidget *shieldingRow = new WContainerWidget( this );
  shieldingRow->addStyleClass( "row" );
  WText *shieldingLabel = new WText( WString::tr("sac-shielding-label"), shieldingRow );
  shieldingLabel->addStyleClass( "label" );
  m_shieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest, shieldingRow );
  m_shieldingSelect->addStyleClass( "input" );
  m_shieldingSelect->materialChanged().connect( this, &SimpleActivityCalc::handleShieldingChanged );
  m_shieldingSelect->materialModified().connect( this, &SimpleActivityCalc::handleShieldingChanged );
  
  // Geometry row (hidden if fixed geometry detector)
  m_geometryRow = new WContainerWidget( this );
  m_geometryRow->addStyleClass( "row" );
  WText *geometryLabel = new WText( WString::tr("sac-geometry-label"), m_geometryRow );
  geometryLabel->addStyleClass( "label" );
  m_geometrySelect = new WComboBox( m_geometryRow );
  m_geometrySelect->addStyleClass( "input" );
  m_geometrySelect->addItem( WString::tr("sac-geometry-point") );
  m_geometrySelect->changed().connect( this, &SimpleActivityCalc::handleGeometryChanged );
  m_geometryRow->hide(); // Initially hidden, shown if not fixed geometry
  

  // Result display
  m_resultText = new WText( this );
  m_resultText->addStyleClass( "row result" );
  
  // Error display
  m_errorText = new WText( this );
  m_errorText->addStyleClass( "row error" );
  m_errorText->hide();
  
  // Advanced button
  m_advancedBtn = new WPushButton( WString::tr("sac-advanced-btn"), this );
  m_advancedBtn->addStyleClass( "row advanced-btn LinkBtn" );
  m_advancedBtn->clicked().connect( this, &SimpleActivityCalc::handleOpenAdvancedTool );
  
  // Connect to InterSpec signals for detector and spectrum changes
  if( m_viewer )
  {
    m_viewer->detectorChanged().connect( this, &SimpleActivityCalc::handleDetectorChanged );
    m_viewer->detectorModified().connect( this, &SimpleActivityCalc::handleDetectorChanged );
    
    // Connect to PeakModel signals to detect peak changes
    PeakModel *peakModel = m_viewer->peakModel();
    if( peakModel )
    {
      peakModel->dataChanged().connect( this, &SimpleActivityCalc::handlePeaksChanged );
      peakModel->rowsRemoved().connect( this, &SimpleActivityCalc::handlePeaksChanged );
      peakModel->rowsInserted().connect( this, &SimpleActivityCalc::handlePeaksChanged );
      peakModel->layoutChanged().connect( this, &SimpleActivityCalc::handlePeaksChanged );
      peakModel->modelReset().connect( this, &SimpleActivityCalc::handlePeaksChanged );
    }
  }
  
  // Initialize display states
  updatePeakList();
  handleDetectorChanged( m_detectorDisplay->detector() );
  handlePeakChanged();
}

void SimpleActivityCalc::render( WFlags<RenderFlag> flags )
{
  WContainerWidget::render( flags );
}

void SimpleActivityCalc::setPeakFromEnergy( const double energy )
{
  if( !m_peakSelect )
    return;
    
  // Find the peak with the closest energy
  for( size_t i = 0; i < m_peakEnergies.size(); ++i )
  {
    if( std::abs(m_peakEnergies[i] - energy) < 1.0 ) // within 1 keV
    {
      m_peakSelect->setCurrentIndex( static_cast<int>(i) );
      handlePeakChanged();
      break;
    }
  }
}

void SimpleActivityCalc::handleAppUrl( std::string uri )
{
}

std::string SimpleActivityCalc::encodeStateToUrl() const
{
  return "dummy-state-url";
}

SimpleActivityCalcState SimpleActivityCalc::currentState() const
{
  SimpleActivityCalcState state;
  
  // Get current peak and nuclide information
  auto peak = getCurrentPeak();
  if( peak )
  {
    state.peakEnergy = peak->mean();
    const SandiaDecay::Nuclide *nuc = peak->parentNuclide();
    if( nuc )
      state.nuclideName = nuc->symbol;
  }
  
  // Get age string
  if( m_ageEdit )
    state.nuclideAgeStr = m_ageEdit->text().toUTF8();
  
  // Get distance string
  if( m_distanceEdit )
    state.distanceStr = m_distanceEdit->text().toUTF8();
  
  // Get geometry type
  if( m_geometrySelect && m_geometrySelect->currentIndex() >= 0 )
  {
    const std::string key = m_geometrySelect->currentText().toUTF8();
    state.geometryType = geometryTypeFromString( key );
  }
  
  // Get shielding info only if not "no shielding"
  if( m_shieldingSelect && !m_shieldingSelect->isNoShielding() )
  {
    try
    {
      state.shielding = m_shieldingSelect->toShieldingInfo();
    }
    catch( const std::exception &e )
    {
      // If error getting shielding, leave it as empty optional
      state.shielding.reset();
    }
  }
  else
  {
    // No shielding selected
    state.shielding.reset();
  }
  
  return state;
}

void SimpleActivityCalc::setState( const SimpleActivityCalcState &state )
{
  // Set peak based on energy and nuclide name
  if( state.peakEnergy > 0 && !state.nuclideName.empty() )
  {
    setPeakFromEnergy( state.peakEnergy );
  }
  
  // Set age string
  if( m_ageEdit && !state.nuclideAgeStr.empty() )
    m_ageEdit->setText( WString::fromUTF8(state.nuclideAgeStr) );
  
  // Set distance string
  if( m_distanceEdit && !state.distanceStr.empty() )
    m_distanceEdit->setText( WString::fromUTF8(state.distanceStr) );
  
  // Set shielding info first (before geometry, so geometry options can be updated appropriately)
  if( m_shieldingSelect )
  {
    if( state.shielding.has_value() )
    {
      try
      {
        m_shieldingSelect->fromShieldingInfo( state.shielding.value() );
      }
      catch( const std::exception &e )
      {
        // If error setting shielding, clear it
        m_shieldingSelect->setToNoShielding();
      }
    }
    else
    {
      // No shielding in state - set to no shielding
      m_shieldingSelect->setToNoShielding();
    }
  }
  
  // Update geometry options based on current detector, peak, and shielding
  updateGeometryOptions();
  
  // Now try to set geometry type if the desired option is available
  if( m_geometrySelect )
  {
    const std::string geometryStr = geometryTypeToString( state.geometryType );
    for( int i = 0; i < m_geometrySelect->count(); ++i )
    {
      if( m_geometrySelect->itemText(i).toUTF8().find( geometryStr ) != std::string::npos )
      {
        m_geometrySelect->setCurrentIndex( i );
        break;
      }
    }
  }
  
  // Update display
  updateResult();
}

void SimpleActivityCalc::updatePeakList()
{
  if( !m_peakSelect )
    return;
    
  m_peakSelect->clear();
  
  if( !m_viewer )
    return;
    
  std::shared_ptr<const SpecMeas> spec = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = m_viewer->peakModel()->peaks();
  
  if( !spec || !peaks )
    return;
    
  // Store peaks in member variable for later access
  m_peakEnergies.clear();
  
  for( const auto &peak : *peaks )
  {
    if( peak && peak->parentNuclide() )
    {
      const SandiaDecay::Nuclide *nuc = peak->parentNuclide();
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "%.1f keV (%s)", peak->mean(), nuc->symbol.c_str() );
      m_peakSelect->addItem( WString::fromUTF8(buffer) );
      
      // Store the peak energy in member vector
      m_peakEnergies.push_back( peak->mean() );
    }
  }
  
  if( m_peakSelect->count() == 0 )
  {
    m_peakSelect->addItem( "No peaks with assigned nuclides" );
    m_peakSelect->setEnabled( false );
    
    // Clear nuclide info and result when no peaks available
    m_nuclideInfo->setText( "" );
    if( m_resultText )
      m_resultText->setText( "" );
    if( m_errorText )
    {
      m_errorText->setText( WString::tr("sac-error-no-peak") );
      m_errorText->show();
    }
  }
  else
  {
    m_peakSelect->setEnabled( true );
    
    // Hide error message when peaks become available
    if( m_errorText )
      m_errorText->hide();
  }
}

void SimpleActivityCalc::handlePeakChanged()
{
  // Update the current peak tracking
  std::shared_ptr<const PeakDef> old_peak = m_currentPeak;
  m_currentPeak = getCurrentPeak();
  
  if( !old_peak || !m_currentPeak || (old_peak->parentNuclide() != m_currentPeak->parentNuclide()) )
  {
    m_ageEdit->setText( "" );
  }
  
  updateNuclideInfo();
  updateGeometryOptions();
  updateResult();
  scheduleRender( AddUndoRedoStep );
}//void handlePeakChanged()


void SimpleActivityCalc::updateNuclideInfo()
{
  auto peak = getCurrentPeak();
  const SandiaDecay::Nuclide *nuc = peak ? peak->parentNuclide() : nullptr;
  
  if( !nuc || (m_peakSelect->currentIndex() < 0) )
  {
    m_nuclideInfo->setText( "" );
    m_ageRow->hide();
    m_resultText->setText( "" );
    m_resultText->hide();
    m_errorText->setText( WString::tr("sac-error-no-peak") );
    m_errorText->show();
    return;
  }
    
  const double energy = peak->mean();
  const double gamma_energy = peak->gammaParticleEnergy();
  
  const SandiaDecay::Transition *trans = peak ? peak->nuclearTransition() : nullptr;
  
  // Show/hide age row based on whether age fitting is appropriate
  if( PeakDef::ageFitNotAllowed( nuc ) || (trans && (trans->parent == nuc)) )
  {
    const double defaultAge = PeakDef::defaultDecayTime( nuc );
    m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(defaultAge) );
    m_ageRow->hide();
  }else
  {
    m_ageRow->show();
    if( m_ageEdit->text().empty() )
    {
      const double defaultAge = PeakDef::defaultDecayTime( nuc );
      m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(defaultAge) );
    }
  }
  
  double age = 0.0;
  try
  {
    age = PhysicalUnits::stringToTimeDuration( m_ageEdit->text().toUTF8() );
  }catch( std::exception &e )
  {
    age = PeakDef::defaultDecayTime( nuc );
    m_ageEdit->setText( PhysicalUnits::printToBestTimeUnits(age) );
  }
  
  
  
  const double activity = 1.0*PhysicalUnits::bq;
  const double real_time = 1.0;
  const bool decay_correct = false;
  const double cluster_num_sigma = 1.25;
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( nuc, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
  const vector< pair<double,double> > energy_widths( 1, {gamma_energy, peak->sigma()} );
  
  map<double,double> energy_gammas_map;
  GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map, energy_widths,
                mixture, activity, age, cluster_num_sigma, gamma_energy, decay_correct, real_time, nullptr, nullptr );
  
  
  assert( energy_gammas_map.size() == 1 );
  double branchingRatio = -1.0;
  if( energy_gammas_map.size() != 1 )
  {
    m_nuclideInfo->setText( "Error determining BR" );
    return;
  }else
  {
    branchingRatio = energy_gammas_map.begin()->second;
  }
  
  
  char buffer[512];
  snprintf( buffer, sizeof(buffer), "%s, %.1f keV (BR: %.2f%%), T½: %s",
           nuc->symbol.c_str(), gamma_energy, branchingRatio*100.0,
           PhysicalUnits::printToBestTimeUnits(nuc->halfLife).c_str() );
  m_nuclideInfo->setText( WString::fromUTF8(buffer) );
}//void SimpleActivityCalc::updateNuclideInfo()


std::shared_ptr<const PeakDef> SimpleActivityCalc::getCurrentPeak() const
{
  if( !m_viewer || !m_peakSelect || m_peakSelect->currentIndex() < 0 )
    return nullptr;
    
  const int index = m_peakSelect->currentIndex();
  if( index >= static_cast<int>(m_peakEnergies.size()) )
    return nullptr;
    
  const double energy = m_peakEnergies[index];
  
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = m_viewer->peakModel()->peaks();
  if( !peaks )
    return nullptr;
    
  // Find the peak with matching energy
  for( const auto &peak : *peaks )
  {
    if( peak && std::abs(peak->mean() - energy) < 0.05 )
      return peak;
  }
  
  return nullptr;
}//std::shared_ptr<const PeakDef> getCurrentPeak() const


int SimpleActivityCalc::findBestReplacementPeak( std::shared_ptr<const PeakDef> targetPeak ) const
{
  if( !targetPeak || !m_peakSelect || m_peakEnergies.empty() )
    return -1;
    
  if( !m_viewer )
    return -1;
    
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = m_viewer->peakModel()->peaks();
  if( !peaks )
    return -1;
    
  const double targetEnergy = targetPeak->mean();
  const SandiaDecay::Nuclide* targetNuclide = targetPeak->parentNuclide();
  
  // First priority: Try to find the exact same peak object (pointer comparison)
  for( size_t i = 0; i < m_peakEnergies.size(); ++i )
  {
    // Find the corresponding peak in the full peak list
    for( const auto &peak : *peaks )
    {
      if( peak && std::abs(peak->mean() - m_peakEnergies[i]) < 0.1 )
      {
        if( peak == targetPeak ) // Same peak object
        {
          return static_cast<int>(i);
        }
        break;
      }
    }
  }
  
  // Second priority: Find exact same peak (same energy and nuclide)
  for( size_t i = 0; i < m_peakEnergies.size(); ++i )
  {
    if( std::abs(m_peakEnergies[i] - targetEnergy) < 0.1 ) // within 0.1 keV
    {
      // Find the corresponding peak in the full peak list
      for( const auto &peak : *peaks )
      {
        if( peak && std::abs(peak->mean() - m_peakEnergies[i]) < 0.1 )
        {
          if( peak->parentNuclide() == targetNuclide )
          {
            return static_cast<int>(i);
          }
          break;
        }
      }
    }
  }
  
  // Third priority: Find closest peak with the same nuclide
  if( targetNuclide )
  {
    int bestIndex = -1;
    double bestEnergyDiff = std::numeric_limits<double>::max();
    
    for( size_t i = 0; i < m_peakEnergies.size(); ++i )
    {
      // Find the corresponding peak in the full peak list
      for( const auto &peak : *peaks )
      {
        if( peak && std::abs(peak->mean() - m_peakEnergies[i]) < 0.1 )
        {
          if( peak->parentNuclide() == targetNuclide )
          {
            double energyDiff = std::abs(m_peakEnergies[i] - targetEnergy);
            if( energyDiff < bestEnergyDiff )
            {
              bestEnergyDiff = energyDiff;
              bestIndex = static_cast<int>(i);
            }
          }
          break;
        }
      }
    }
    
    if( bestIndex >= 0 )
      return bestIndex;
  }
  
  // Fourth priority: Find closest peak with any nuclide
  int bestIndex = -1;
  double bestEnergyDiff = std::numeric_limits<double>::max();
  
  for( size_t i = 0; i < m_peakEnergies.size(); ++i )
  {
    double energyDiff = std::abs(m_peakEnergies[i] - targetEnergy);
    if( energyDiff < bestEnergyDiff )
    {
      bestEnergyDiff = energyDiff;
      bestIndex = static_cast<int>(i);
    }
  }
  
  return bestIndex;
}

void SimpleActivityCalc::handleAgeChanged()
{
  updateResult();
  scheduleRender( AddUndoRedoStep );
}

void SimpleActivityCalc::handleDistanceChanged()
{
  updateResult();
  scheduleRender( AddUndoRedoStep );
}

void SimpleActivityCalc::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf )
{
  // Check if this is a fixed geometry detector and hide/show distance and geometry rows
  if( new_drf && m_distanceRow && m_geometryRow )
  {
    const bool isFixedGeom = new_drf->isFixedGeometry();
    m_distanceRow->setHidden( isFixedGeom );
    m_geometryRow->setHidden( isFixedGeom );
  }
  
  updateGeometryOptions();
  updateResult();
  scheduleRender( AddUndoRedoStep );
}

void SimpleActivityCalc::handleGeometryChanged()
{
  updateResult();
  scheduleRender( AddUndoRedoStep );
}

void SimpleActivityCalc::handleShieldingChanged()
{
  updateGeometryOptions();
  updateResult();
  scheduleRender( AddUndoRedoStep );
}

void SimpleActivityCalc::handleSpectrumChanged()
{
  updatePeakList();
  updateResult();
}

void SimpleActivityCalc::handlePeaksChanged()
{
  // Store reference to current peak before updating
  std::shared_ptr<const PeakDef> previousPeak = m_currentPeak;
  
  // Update the peak list
  updatePeakList();
  
  // Try to find the best replacement peak if we had a selection before
  if( previousPeak && m_peakSelect && m_peakSelect->count() > 0 && m_peakSelect->isEnabled() )
  {
    int bestIndex = findBestReplacementPeak( previousPeak );
    if( bestIndex >= 0 )
    {
      m_peakSelect->setCurrentIndex( bestIndex );
    }
    else if( m_peakSelect->count() > 0 )
    {
      // Fallback to first available peak
      m_peakSelect->setCurrentIndex( 0 );
    }
  }
  else if( m_peakSelect && m_peakSelect->count() > 0 && m_peakSelect->isEnabled() )
  {
    // If we didn't have a selection before but now have peaks, select the first one
    m_peakSelect->setCurrentIndex( 0 );
  }
  
  // Update the current peak tracking after selection
  m_currentPeak = getCurrentPeak();
  
  // Update displays
  updateNuclideInfo();
  updateResult();
}

void SimpleActivityCalc::scheduleRender( RenderActions action )
{
  m_renderFlags |= action;
  WContainerWidget::scheduleRender();
}

void SimpleActivityCalc::updateResult()
{
  if( !m_resultText || !m_errorText )
    return;
    
  // Clear previous results
  m_resultText->setText( "" );
  m_errorText->setText( "" );
  m_errorText->hide();
  
  // Validate inputs
  auto peak = getCurrentPeak();
  if( !peak )
  {
    m_errorText->setText( WString::tr("sac-error-no-peak") );
    m_errorText->show();
    return;
  }
  
  if( !m_detectorDisplay || !m_detectorDisplay->detector() )
  {
    m_errorText->setText( WString::tr("sac-error-no-detector") );
    m_errorText->show();
    return;
  }
  
  if( m_distanceEdit && !m_distanceEdit->validate() )
  {
    m_errorText->setText( WString::tr("sac-error-invalid-distance") );
    m_errorText->show();
    return;
  }
  
  try
  {
    const SimpleActivityCalcInput input = createCalcInput();
    const SimpleActivityCalcResult result = performCalculation( input );
    
    if( !result.successful )
    {
      m_errorText->setText( WString::fromUTF8("Calculation error: ") + WString::fromUTF8(result.errorMessage) );
      m_errorText->show();
      return;
    }
    
    if( result.activity > 0.0 )
    {
      string resultStr;
      
      switch( input.geometryType )
      {
        case SimpleActivityGeometryType::Point:
        case SimpleActivityGeometryType::TraceSrc:
        case SimpleActivityGeometryType::SelfAttenuating:
        {
          resultStr = "Activity: " + PhysicalUnits::printToBestActivityUnits(result.activity);
          if( result.activityUncertainty > 0.0 )
            resultStr += " ± " + PhysicalUnits::printToBestActivityUnits(result.activityUncertainty);
          
          if( input.geometryType == SimpleActivityGeometryType::SelfAttenuating )
          {
            const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
            const SandiaDecay::Nuclide * const nuc = input.peak ? input.peak->parentNuclide() : nullptr;
            const SandiaDecay::Element *el = (db && nuc) ? db->element( nuc->atomicNumber ) : nullptr;
            assert( el );
            if( el )
              resultStr += "<br/>(used " + el->symbol + " 100% " + nuc->symbol + ")";
            
            if( result.sourceDimensions > 0.0 )
            {
              const std::string dimStr = PhysicalUnits::printToBestLengthUnits(result.sourceDimensions);
              resultStr += "<br/>Source radius: " + dimStr;
            }
            if( result.nuclideMass > 0.0 )
            {
              const std::string srcMassStr = PhysicalUnits::printToBestMassUnits(result.nuclideMass);
              resultStr += "<br/>" + string(nuc ? nuc->symbol : "null") + " mass: " + srcMassStr;
            }
          }//if( input.geometryType == SimpleActivityGeometryType::SelfAttenuating )
          
          break;
        }//case SimpleActivityGeometryType::Point, TraceSrc, SelfAttenuating
          
        case SimpleActivityGeometryType::Plane:
        {
          // We will make the activity per cm2.
          assert( (input.shielding.size() >= 1) && (input.shielding[0].m_traceSources.size() == 1) );
          
          if( (input.shielding.size() >= 1) && (input.shielding[0].m_traceSources.size() > 0) )
          {
            const ShieldingSourceFitCalc::ShieldingInfo &shield = input.shielding[0];
            assert( shield.m_traceSources[0].m_type == GammaInteractionCalc::TraceActivityType::TotalActivity );
            
            const double surface_area = PhysicalUnits::pi * std::pow(shield.m_dimensions[0], 2.0);
            const double surface_area_cm2 = surface_area / PhysicalUnits::cm2;
            const double activity_cm2 = result.activity / surface_area_cm2;
            
            resultStr = "Act/cm<sup>2</sup>: " + PhysicalUnits::printToBestActivityUnits(activity_cm2);
            if( result.activityUncertainty > 0.0 )
            {
              const double uncert_cm2 = result.activityUncertainty / surface_area_cm2;
              resultStr += " ± " + PhysicalUnits::printToBestActivityUnits(uncert_cm2);
            }
            resultStr += " /cm<sup>2</sup>";
          }else
          {
            resultStr = "Unexpected error converting result";
          }
          break;
        }//case SimpleActivityGeometryType::Plane:
      }//switch( input.geometryType )
      
      m_resultText->setText( WString::fromUTF8(resultStr) );
    }
    else
    {
      m_errorText->setText( "Unable to calculate activity" );
      m_errorText->show();
    }
  }
  catch( const std::exception &e )
  {
    m_errorText->setText( WString::fromUTF8("Calculation error: ") + WString::fromUTF8(e.what()) );
    m_errorText->show();
  }
}

SimpleActivityCalcInput SimpleActivityCalc::createCalcInput() const
{
  SimpleActivityCalcInput input;
  
  // Get current peak
  const shared_ptr<const PeakDef> spectrum_peak = getCurrentPeak();
  
  if( !spectrum_peak || !spectrum_peak->parentNuclide() )
    throw std::runtime_error( "No valid peak selected" );
  
  // Make a copy of the peak, so we dont mess the original up, and make sure
  //  it is marked to actually be used for act/shield fit
  shared_ptr<PeakDef> peak = make_shared<PeakDef>( *spectrum_peak );
  peak->setContinuum( make_shared<PeakContinuum>( *peak->continuum() ));
  peak->useForShieldingSourceFit( true );
  input.peak = peak;
  
  const SandiaDecay::Nuclide *nuc = input.peak ? input.peak->parentNuclide() : nullptr;
  const SandiaDecay::Transition *trans = input.peak ? input.peak->nuclearTransition() : nullptr;
  
  
  // Get foreground spectrum
  if( !m_viewer )
    throw std::runtime_error( "No InterSpec viewer available" );
    
  shared_ptr<const SpecMeas> foregroundSpecMeas = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( !foregroundSpecMeas )
    throw std::runtime_error( "No foreground spectrum available" );
  
  // Get the detector
  input.detector = foregroundSpecMeas->detector();
  if( !input.detector )
    throw std::runtime_error( "No detector response available" );
  
  // Get spectra
  input.foreground = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  input.background = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Background);
  
  if( !input.foreground )
    throw std::runtime_error( "Could not get foreground spectrum" );
  
  // Set up geometry and distance
  if( input.detector->isFixedGeometry() )
  {
    input.distance = 0.0;
    input.geometryType = SimpleActivityGeometryType::Point;
  }else
  {
    try
    {
      input.distance = PhysicalUnits::stringToDistance( m_distanceEdit->text().toUTF8() );
    }catch( const std::exception &e )
    {
      throw std::runtime_error( "Invalid distance: " + std::string(e.what()) );
    }
    
    const std::string geom_key = m_geometrySelect->currentText().toUTF8();
    input.geometryType = geometryTypeFromString( geom_key );
    
    switch( input.geometryType )
    {
      case SimpleActivityGeometryType::Point:
      {
        if( !m_shieldingSelect->isNoShielding() )
        {
          ShieldingSourceFitCalc::ShieldingInfo info = m_shieldingSelect->toShieldingInfo();
          input.shielding.push_back( info );
        }//if( !m_shieldingSelect->isNoShielding() )
        
        break;
      }//case SimpleActivityGeometryType::Point:
        
      case SimpleActivityGeometryType::Plane:
      {
        ShieldingSourceFitCalc::TraceSourceInfo trace_info;
        trace_info.m_type = GammaInteractionCalc::TraceActivityType::TotalActivity;
        trace_info.m_fitActivity = true;
        trace_info.m_nuclide = nuc;
        trace_info.m_activity = 1.0E-6*PhysicalUnits::curie; //starting value only
        trace_info.m_relaxationDistance = 0.0; //not applicable
        
        ShieldingSourceFitCalc::ShieldingInfo trace_shield;
        trace_shield.m_geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
        trace_shield.m_isGenericMaterial = false;
        trace_shield.m_forFitting = true;
        const Material *air = m_materialDB->material( "Air" );
        assert( air );
        if( !air )
          throw runtime_error( "No Air material in DB?!?" );
        trace_shield.m_material = make_shared<Material>( *air );
        
        trace_shield.m_dimensions[0] = 50.0*PhysicalUnits::meter; //radius
        trace_shield.m_dimensions[1] = 1.0*PhysicalUnits::mm;     //half-length
        trace_shield.m_dimensions[2] = 0.0; //N/A
        trace_shield.m_fitDimensions[0] = trace_shield.m_fitDimensions[1] = trace_shield.m_fitDimensions[2] = false;
        
        trace_shield.m_traceSources.push_back( trace_info );
        
        input.shielding.insert( begin(input.shielding), trace_shield ); //Make sure goes in the front
        
        if( !m_shieldingSelect->isNoShielding() )
        {
          ShieldingSourceFitCalc::ShieldingInfo info = m_shieldingSelect->toShieldingInfo();
          if( !info.m_isGenericMaterial )
          {
            info.m_geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
            info.m_dimensions[1] = info.m_dimensions[0];
            info.m_dimensions[0] = info.m_dimensions[2] = 0.0;
            info.m_fitDimensions[0] = info.m_fitDimensions[1] = info.m_fitDimensions[2];
          }//if( !info.m_isGenericMaterial )
          
          input.shielding.push_back( info );
        }//if( !m_shieldingSelect->isNoShielding() )
        
        break;
      }//case SimpleActivityGeometryType::Plane:
        
        
      case SimpleActivityGeometryType::SelfAttenuating:
      {
        if( m_shieldingSelect->isNoShielding() || m_shieldingSelect->isGenericMaterial() )
          throw runtime_error( "You must have a material shielding defined for self-attenuation." );
          
        ShieldingSourceFitCalc::ShieldingInfo info = m_shieldingSelect->toShieldingInfo();
        assert( info.m_material );
        if( !info.m_material )
          throw runtime_error( "No valid material defined for self-attenuation." );
        
        info.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
        info.m_isGenericMaterial = false;
        info.m_dimensions[0] = 1.0*PhysicalUnits::cm; //staring radius
        info.m_fitDimensions[0] = true;
        info.m_fitDimensions[1] = info.m_fitDimensions[2] = false;
        
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        if( !db )
          throw runtime_error( "Invalid SandiaDecayDataBase?!?" );
        
        const SandiaDecay::Element * const el = db->element( nuc->atomicNumber );
        assert( el );
        if( !el )
          throw runtime_error( "Invalid element?!?" );
        
        bool has_element = false;
        for( const pair<const SandiaDecay::Element *,float> &el_frac : info.m_material->elements )
        {
          if( el_frac.first != el )
            continue;
          has_element = true;
          break;
        }
        
        for( const pair<const SandiaDecay::Nuclide *,float> &nuc_frac : info.m_material->nuclides )
        {
          if( has_element )
            break;
        }
        
        if( !has_element )
          throw runtime_error( "Shielding material selected does not contain the element for the current source nuclide." );
          
        // TODO: currently we are just putting the isotope in as 100% of the element... need to give the user an option to set this
        
        double massFraction = 1.0;
        info.m_nuclideFractions_[el].push_back( {nuc,massFraction,false} );
        
        input.shielding.push_back( info );
        
        break;
      }//case SimpleActivityGeometryType::SelfAttenuating:
        
        
      case SimpleActivityGeometryType::TraceSrc:
      {
        if( m_shieldingSelect->isNoShielding() || m_shieldingSelect->isGenericMaterial() )
          throw runtime_error( "You must have a material shielding defined for trace-source." );
        
        ShieldingSourceFitCalc::ShieldingInfo info = m_shieldingSelect->toShieldingInfo();
        assert( info.m_material );
        if( !info.m_material )
          throw runtime_error( "No valid material defined for trace-source." );
        
        info.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
        info.m_isGenericMaterial = false;
        if( info.m_dimensions[0] <= 1.0*PhysicalUnits::um )
          throw runtime_error( "Shielding for trace source must have non-zero thickness" );
        
        info.m_fitDimensions[0] = info.m_fitDimensions[1] = info.m_fitDimensions[2] = false;
        
        ShieldingSourceFitCalc::TraceSourceInfo trace_info;
        trace_info.m_type = GammaInteractionCalc::TraceActivityType::TotalActivity;
        trace_info.m_fitActivity = true;
        trace_info.m_nuclide = nuc;
        trace_info.m_activity = 1.0E-6*PhysicalUnits::curie;  //starting value only
        trace_info.m_relaxationDistance = 0.0; //not applicable
        
        info.m_traceSources.push_back( trace_info );
        
        input.shielding.push_back( info );
        break;
      }//case SimpleActivityGeometryType::TraceSrc
    }//switch( input.geometryType )
  }//if( !input.detector->isFixedGeometry() ) / else
  
  
  const bool is_parent_gamma = trans && trans->parent && (trans->parent == nuc);
  
  // Set up age if specified
  input.age = PeakDef::defaultDecayTime( nuc );
  if( !is_parent_gamma && !PeakDef::ageFitNotAllowed(nuc) )
  {
    try
    {
      input.age = PhysicalUnits::stringToTimeDuration( m_ageEdit->text().toUTF8() );
    }catch( const std::exception &e )
    {
      // Use default age on error
    }
  }//if( we want an age )
  
  return input;
}//createCalcInput()


SimpleActivityCalcResult SimpleActivityCalc::performCalculation( const SimpleActivityCalcInput& input )
{
  SimpleActivityCalcResult result;
  
  try
  {
    if( !input.peak )
      throw runtime_error( "No peak provided." );
      
    if( !input.detector || !input.detector->isValid() )
      throw runtime_error( "No detector efficiency provided." );
    
    const SandiaDecay::Nuclide *nuc = input.peak->parentNuclide();
    if( !nuc )
      throw runtime_error( "No nuclide associated with input peak." );
    
    if( input.background_peak )
    {
      if( fabs(input.background_peak->mean() - input.peak->mean()) > input.peak->fwhm() )
        throw runtime_error( "Background peak mean greater than 1 FWHM away from foreground peak." );
    }
    
    assert( !input.detector->isFixedGeometry() || (input.geometryType == SimpleActivityGeometryType::Point) );
    if( input.detector->isFixedGeometry() && (input.geometryType != SimpleActivityGeometryType::Point) )
      throw runtime_error( "Fixed geometry detector efficiencies must be point geometry" );
    
    // Set up the data structures required by the framework
    ShieldingSourceFitCalc::ShieldingSourceFitOptions fit_options;
    fit_options.multiple_nucs_contribute_to_peaks = true;
    fit_options.attenuate_for_air = true;
    fit_options.account_for_decay_during_meas = false;
    fit_options.photopeak_cluster_sigma = 1.25;
    fit_options.background_peak_subtract = !!input.background_peak;
    fit_options.same_age_isotopes = false;
    
    // Set up source definition
    ShieldingSourceFitCalc::SourceFitDef source;
    source.nuclide = nuc;
    source.activity = 1.0 * PhysicalUnits::curie; // will be fit
    source.fitActivity = true;
    source.fitAge = false; // Don't fit age, use user-specified value
    source.age = input.age;
    source.ageDefiningNuc = nullptr;
    
    const vector<ShieldingSourceFitCalc::ShieldingInfo> &shielding = input.shielding;
    GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::Spherical;
    
    switch( input.geometryType )
    {
      case SimpleActivityGeometryType::Point:
      {
        geometry = GammaInteractionCalc::GeometryType::Spherical;
        source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
        break;
      }//case SimpleActivityGeometryType::Point:
        
      case SimpleActivityGeometryType::Plane:
      {
        if( input.detector->isFixedGeometry() )
          throw std::logic_error( "Infinite plane source not allowed for fixed geometry detector" );
        
        geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
        source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Trace;
        
        assert( (shielding.size() == 1) || (shielding.size() == 2) );
        if( (shielding.size() != 1) && (shielding.size() != 2) )
          throw runtime_error( "Invalid inifinite plane input - must have 1 or 2 shieldings." );
        
        if( shielding[0].m_traceSources.size() != 1 )
          throw runtime_error( "Invalid inifinite plane input - must have trace-source defined." );
  
        assert( shielding[0].m_geometry == GammaInteractionCalc::GeometryType::CylinderEndOn );
        if( shielding[0].m_geometry != GammaInteractionCalc::GeometryType::CylinderEndOn )
          throw logic_error( "Unexpected trace source geometery" );
        
        break;
      }//case SimpleActivityGeometryType::Plane:
        
      case SimpleActivityGeometryType::SelfAttenuating:
      {
        if( input.detector->isFixedGeometry() )
          throw std::logic_error( "Self attenuating source not allowed for fixed geometry detector" );
        
        source.fitActivity = false;
        geometry = GammaInteractionCalc::GeometryType::Spherical;
        source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Intrinsic;
        if( shielding.empty() )
          throw runtime_error( "You must have a shielding defined for a self-attenuating source." );
        
        const ShieldingSourceFitCalc::ShieldingInfo &shield = shielding.front();
        if( !shield.m_material || shield.m_isGenericMaterial )
          throw logic_error( "No material defined for self-attenuating source" );
        
        if( shield.m_geometry != GammaInteractionCalc::GeometryType::Spherical )
          throw logic_error( "Geometry for self-attenuating source not spherical" );
        
        if( (shield.m_nuclideFractions_.size() != 1)
           || (begin(shield.m_nuclideFractions_)->second.size() != 1) )
          throw logic_error( "No source nuclides defined for self attenuation" );
        
        if( get<0>(begin(shield.m_nuclideFractions_)->second.front()) != nuc )
          throw logic_error( "Self-attenuating nuclide doesnt match to peak" );
        
        break;
      }//case SimpleActivityGeometryType::SelfAttenuating:
        
      case SimpleActivityGeometryType::TraceSrc:
      {
        if( input.detector->isFixedGeometry() )
          throw std::logic_error( "Trace source not allowed for fixed geometry detector" );
        
        geometry = GammaInteractionCalc::GeometryType::Spherical;
        source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Trace;
        if( shielding.empty() )
          throw runtime_error( "You must have a shielding defined for a trace source." );
        
        break;
      }//case SimpleActivityGeometryType::TraceSrc:
    }//switch( input.geometryType )
    
    if( input.detector->isFixedGeometry() )
      geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  
      
    vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions( 1, source );
    
    
    // Create peak collections for the calculation
    deque<std::shared_ptr<const PeakDef>> foreground_peaks( 1, input.peak );
    shared_ptr<const deque<shared_ptr<const PeakDef>>> background_peaks;
    if( input.background_peak )
    {
      auto bg_peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
      bg_peaks->push_back( input.background_peak );
      background_peaks = bg_peaks;
    }
    
      
    // Create the fitting function
    pair<shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
      GammaInteractionCalc::ShieldingSourceChi2Fcn::create( input.distance, geometry,
                                                           shielding, src_definitions, input.detector,
                                                           input.foreground, input.background, 
                                                           foreground_peaks, 
                                                           background_peaks,
                                                           fit_options );
    
    auto inputPrams = std::make_shared<ROOT::Minuit2::MnUserParameters>();
    *inputPrams = fcn_pars.second;
    
    auto progress = std::make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
    auto fit_results = std::make_shared<ShieldingSourceFitCalc::ModelFitResults>();
    
    auto progress_fcn = [](){
      // No progress updates needed for simple case
    };
    
    bool finished_fit_called = false;
    auto finished_fcn = [&finished_fit_called](){
      finished_fit_called = true;
    };
    
    // Do the fit
    ShieldingSourceFitCalc::fit_model( "", fcn_pars.first, inputPrams, progress, progress_fcn, fit_results, finished_fcn );
    
    if( !finished_fit_called )
      throw std::runtime_error( "Fit was cancelled or did not complete" );
      
    if( fit_results->successful != ShieldingSourceFitCalc::ModelFitResults::FitStatus::Final )
      throw std::runtime_error( "Fit was not successful" );
      
    // Extract the activity result
    if( fit_results->fit_src_info.empty() )
      throw std::runtime_error( "No fit results available" );
      
    const ShieldingSourceFitCalc::IsoFitStruct &fitResult = fit_results->fit_src_info[0];
    
    if( fitResult.nuclide != source.nuclide )
      throw std::runtime_error( "Fit result nuclide mismatch" );
      
    result.successful = true;
    result.activity = fitResult.activity;
    result.activityUncertainty = fitResult.activityUncertainty;
    
    // Calculate nuclide mass
    result.nuclideMass = (result.activity / nuc->activityPerGram()) * PhysicalUnits::gram;
    
    
    // Check if this is self-attenuating
    result.isSelfAttenuating = (input.geometryType == SimpleActivityGeometryType::SelfAttenuating);
    
    if( result.isSelfAttenuating && !input.shielding.empty() )
    {
      // For self-attenuating sources, calculate source dimensions and mass
      const auto& shield = input.shielding[0];
      if( shield.m_geometry == GammaInteractionCalc::GeometryType::Spherical )
      {
        result.sourceDimensions = shield.m_dimensions[0]; // radius
        if( shield.m_material )
        {
          const double volume = (4.0/3.0) * M_PI * pow(result.sourceDimensions, 3);
          result.sourceMass = volume * shield.m_material->density;
        }
      }
    }
  }catch( const std::exception &e )
  {
    result.successful = false;
    result.errorMessage = e.what();
  }//try / catch to compute an answer
  
  return result;
}

void SimpleActivityCalc::updateGeometryOptions()
{
  if( !m_geometrySelect || !m_geometryRow )
    return;
    
  // Get current detector
  auto detector = m_detectorDisplay ? m_detectorDisplay->detector() : nullptr;
  
  // If fixed geometry detector, hide the entire row and set to Point
  if( !detector || detector->isFixedGeometry() )
  {
    m_geometryRow->hide();
    // Set to point source but don't emit signals
    m_geometrySelect->setCurrentIndex( static_cast<int>(SimpleActivityGeometryType::Point) );
    return;
  }
  
  // Show the row for non-fixed geometry detectors
  m_geometryRow->show();
  
  // Get current peak to check nuclide compatibility with shielding
  auto peak = getCurrentPeak();
  
  // Get current shielding info to check for self-attenuation and trace source compatibility
  bool allowSelfAttenuation = false;
  bool allowTraceSrc = false;
  
  if( peak && peak->parentNuclide() && m_shieldingSelect && !m_shieldingSelect->isNoShielding() )
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo shielding = m_shieldingSelect->toShieldingInfo();
      const SandiaDecay::Nuclide *peakNuclide = peak->parentNuclide();
      
      // Check if we should allow trace source option
      // Allow TraceSrc if we have non-generic shielding with non-zero thickness (radius)
      if( !shielding.m_isGenericMaterial && shielding.m_material && shielding.m_dimensions[0] > 1e-6 )
      {
        allowTraceSrc = true;
      }
      
      // Check if the shielding material contains the same element as the peak nuclide
      if( shielding.m_material && peakNuclide )
      {
        const int peakAtomicNumber = peakNuclide->atomicNumber;
        
        // Check material elements by atomic number
        for( const auto &element : shielding.m_material->elements )
        {
          if( element.first && element.first->atomicNumber == peakAtomicNumber && element.second > 1e-6 )
          {
            allowSelfAttenuation = true;
            break;
          }
        }
        
        // Also check nuclide fractions map for exact nuclide match or same element
        for( const auto &elementToNucs : shielding.m_nuclideFractions_ )
        {
          const SandiaDecay::Element *shieldingElement = elementToNucs.first;
          if( shieldingElement && shieldingElement->atomicNumber == peakAtomicNumber )
          {
            // Same element - check if our specific nuclide is present or if there are any nuclides
            allowSelfAttenuation = true;
            
            const auto &nuclideVec = elementToNucs.second;
            for( const auto &nucTuple : nuclideVec )
            {
              const SandiaDecay::Nuclide *nuc = std::get<0>(nucTuple);
              if( nuc == peakNuclide )
              {
                // Exact nuclide match
                allowSelfAttenuation = true;
                break;
              }
            }
            break;
          }
        }
      }
    }
    catch( const std::exception &e )
    {
      // If can't get shielding info, assume no self-attenuation
      allowSelfAttenuation = false;
    }
  }
  
  // Store current selection before clearing
  const int currentIndex = m_geometrySelect->currentIndex();
  Wt::WString currentSelection;
  if( currentIndex >= 0 && currentIndex < m_geometrySelect->count() )
  {
    currentSelection = m_geometrySelect->itemText( currentIndex );
  }
  
  // Clear and rebuild combo box
  m_geometrySelect->clear();
  
  // Always add point and plane options
  m_geometrySelect->addItem( WString::tr("sac-geometry-point") );
  m_geometrySelect->addItem( WString::tr("sac-geometry-plane") );
  
  // Add self-attenuation option only if compatible
  if( allowSelfAttenuation )
  {
    m_geometrySelect->addItem( WString::tr("sac-geometry-self-atten") );
  }
  
  // Add trace source option only if we have appropriate shielding
  if( allowTraceSrc )
  {
    m_geometrySelect->addItem( WString::tr("sac-geometry-trace-src") );
  }
  
  // Try to restore previous selection if still available
  int newIndex = -1;
  if( !currentSelection.empty() )
  {
    for( int i = 0; i < m_geometrySelect->count(); ++i )
    {
      if( m_geometrySelect->itemText(i).key() == currentSelection.key() )
      {
        newIndex = i;
        break;
      }
    }
  }
  
  // If previous selection is no longer available, default to Point
  if( newIndex < 0 )
  {
    newIndex = 0; // Point source
  }
  
  m_geometrySelect->setCurrentIndex( newIndex );
}

void SimpleActivityCalc::handleOpenAdvancedTool()
{
  // TODO: Implement opening ShieldingSourceDisplay with equivalent configuration
  // This will require:
  // 1. Create ShieldingSourceDisplay 
  // 2. Set it to equivalent state based on current SimpleActivityCalc state
  // 3. Add undo/redo point for both tools and peaks
  // 4. Close current dialog
}

SimpleActivityGeometryType SimpleActivityCalc::geometryTypeFromString( const std::string& key )
{
  if( key.find("Point") != std::string::npos || key.find("point") != std::string::npos )
    return SimpleActivityGeometryType::Point;
  else if( key.find("Plane") != std::string::npos || key.find("plane") != std::string::npos )
    return SimpleActivityGeometryType::Plane;
  else if( key.find("Self") != std::string::npos || key.find("self") != std::string::npos )
    return SimpleActivityGeometryType::SelfAttenuating;
  else if( key.find("Trace") != std::string::npos || key.find("trace") != std::string::npos )
    return SimpleActivityGeometryType::TraceSrc;
  else
    return SimpleActivityGeometryType::Point; // Default fallback
}

std::string SimpleActivityCalc::geometryTypeToString( SimpleActivityGeometryType type )
{
  switch( type )
  {
    case SimpleActivityGeometryType::Point: return "Point";
    case SimpleActivityGeometryType::Plane: return "Plane";
    case SimpleActivityGeometryType::SelfAttenuating: return "Self-attenuating";
    case SimpleActivityGeometryType::TraceSrc: return "Trace source";
    default: return "Point";
  }
}
