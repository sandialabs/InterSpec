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
#include <sstream>
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

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/SimpleActivityCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "SandiaDecay/SandiaDecay.h"

#include "Minuit2/MnUserParameters.h"

#if( USE_QR_CODES )
#include <Wt/Utils>
#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;

SimpleActivityCalcState::SimpleActivityCalcState()
  : peakEnergy( -1 ),
    nuclideAgeStr( "" ),
    distanceStr( "" ),
    geometryType( SimpleActivityGeometryType::Point ),
    shielding{}
{
}

bool SimpleActivityCalcState::operator==( const SimpleActivityCalcState &rhs ) const
{
  return (peakEnergy == rhs.peakEnergy)
         && (nuclideAgeStr == rhs.nuclideAgeStr)
         && (distanceStr == rhs.distanceStr)
         && (geometryType == rhs.geometryType)
         // TODO: need to implement equality check for `ShieldingSourceFitCalc::ShieldingInfo` && (shielding == rhs.shielding)
  ;
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

void SimpleActivityCalcState::serialize( std::ostream &out ) const
{
}

void SimpleActivityCalcState::deserialize( std::istream &in )
{
}

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
  return SimpleActivityCalcState();
}

void SimpleActivityCalc::setState( const SimpleActivityCalcState &state )
{
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
  const double cluster_num_sigma = 1.5;
  
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
    const double activity = calculateActivity();
    if( activity > 0.0 )
    {
      char buffer[256];
      const std::string actStr = PhysicalUnits::printToBestActivityUnits(activity);
      snprintf( buffer, sizeof(buffer), "Activity: %s", actStr.c_str() );
      m_resultText->setText( WString::fromUTF8(buffer) );
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

double SimpleActivityCalc::calculateActivity() const
{
  // Get current peak
  auto peak = getCurrentPeak();
  if( !peak || !peak->parentNuclide() )
    throw std::runtime_error( "No valid peak selected" );
    
  // Get detector
  auto detector = m_detectorDisplay ? m_detectorDisplay->detector() : nullptr;
  if( !detector )
    throw std::runtime_error( "No detector response available" );
    
  // Get foreground spectrum
  if( !m_viewer )
    throw std::runtime_error( "No InterSpec viewer available" );
    
  auto foregroundSpecMeas = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( !foregroundSpecMeas )
    throw std::runtime_error( "No foreground spectrum available" );
    
  // Get background spectrum (optional)
  auto backgroundSpecMeas = m_viewer->measurment( SpecUtils::SpectrumType::Background );
  
  // Convert to SpecUtils::Measurement using the displayed sample numbers
  std::set<int> foreground_samples = m_viewer->displayedSamples( SpecUtils::SpectrumType::Foreground );
  std::shared_ptr<const SpecUtils::Measurement> foreground = 
    foregroundSpecMeas->sum_measurements( foreground_samples, foregroundSpecMeas->detector_names(), nullptr );
    
  if( !foreground )
    throw std::runtime_error( "Could not sum foreground measurements" );
    
  std::shared_ptr<const SpecUtils::Measurement> background;
  if( backgroundSpecMeas )
  {
    std::set<int> background_samples = m_viewer->displayedSamples( SpecUtils::SpectrumType::Background );
    background = backgroundSpecMeas->sum_measurements( background_samples, backgroundSpecMeas->detector_names(), nullptr );
  }
  
  // Get current peaks
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> allPeaks = m_viewer->peakModel()->peaks();
  if( !allPeaks )
    throw std::runtime_error( "No peaks available" );
    
  // Set up the data structures required by the framework
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> shield_definitions;
  std::vector<ShieldingSourceFitCalc::SourceFitDef> src_definitions;
  ShieldingSourceFitCalc::ShieldingSourceFitOptions fit_options;
  
  // Set up geometry and distance
  double distance = 1.0 * PhysicalUnits::meter; // default 1 meter
  GammaInteractionCalc::GeometryType geometry = GammaInteractionCalc::GeometryType::Spherical;
  
  if( detector->isFixedGeometry() )
  {
    distance = detector->detectorDiameter(); // Use detector diameter for fixed geometry
  }
  else if( m_distanceEdit && !m_distanceEdit->text().empty() )
  {
    try
    {
      distance = PhysicalUnits::stringToDistance( m_distanceEdit->text().toUTF8() );
    }
    catch( const std::exception &e )
    {
      throw std::runtime_error( "Invalid distance: " + std::string(e.what()) );
    }
  }
  
  // Set up geometry based on user selection
  if( m_geometrySelect && m_geometrySelect->currentIndex() >= 0 && !detector->isFixedGeometry() )
  {
    switch( m_geometrySelect->currentIndex() )
    {
      case 0: geometry = GammaInteractionCalc::GeometryType::Spherical; break;
      case 1: geometry = GammaInteractionCalc::GeometryType::Rectangular; break; // Infinite plane approximation
      case 2: geometry = GammaInteractionCalc::GeometryType::Spherical; break; // Spherical trace source
      case 3: geometry = GammaInteractionCalc::GeometryType::Spherical; break; // Self-attenuating sphere
      default: geometry = GammaInteractionCalc::GeometryType::Spherical; break;
    }
  }
  
  // Set up shielding from ShieldingSelect widget
  if( m_shieldingSelect && !m_shieldingSelect->isNoShielding() )
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo shielding = m_shieldingSelect->toShieldingInfo();
      shield_definitions.push_back( shielding );
    }
    catch( const std::exception &e )
    {
      // If no shielding, that's okay - just use empty vector
    }
  }
  
  const SandiaDecay::Nuclide *nuc = peak ? peak->parentNuclide() : nullptr;
  const SandiaDecay::Transition *trans = peak ? peak->nuclearTransition() : nullptr;
  const bool is_parent_gamma = trans && trans->parent && (trans->parent == nuc);
  
  // Set up source definition
  ShieldingSourceFitCalc::SourceFitDef source;
  source.nuclide = nuc;
  source.activity = 1.0 * PhysicalUnits::curie; // Initial guess - will be fit
  source.fitActivity = true;
  source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
  source.fitAge = false; // Don't fit age, use user-specified value
  source.age = PeakDef::defaultDecayTime( source.nuclide );
  
  // Set up age if specified
  if( !is_parent_gamma && !PeakDef::ageFitNotAllowed(source.nuclide) )
  {
    try
    {
      source.age = PhysicalUnits::stringToTimeDuration( m_ageEdit->text().toUTF8() );
    }catch( const std::exception &e )
    {
      
    }
  }//if( !is_parent_gamma && !PeakDef::ageFitNotAllowed(source.nuclide) )
  
  
  source.ageDefiningNuc = nullptr;
  src_definitions.push_back( source );
  
  // Set up fit options
  fit_options.multiple_nucs_contribute_to_peaks = true; // Incase there are multiple
  //fit_options.photopeak_cluster_sigma = 1.25;
  fit_options.attenuate_for_air = true;
  // We should background subtract only if we have a compatible background peak
  fit_options.background_peak_subtract = (background != nullptr);
  fit_options.same_age_isotopes = false;
  
  // Filter peaks to just those with the selected nuclide
  std::deque<std::shared_ptr<const PeakDef>> foreground_peaks;
  std::deque<std::shared_ptr<const PeakDef>> background_peaks;
  
  // I think we just want to use the one peak the user has selected.
  for( const auto &p : *allPeaks )
  {
    if( p && p->parentNuclide() == source.nuclide )
    {
      foreground_peaks.push_back( p );
    }
  }
  
  if( foreground_peaks.empty() )
    throw std::runtime_error( "No peaks found for selected nuclide" );
    
  // Create the fitting function
  std::pair<std::shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> fcn_pars =
    GammaInteractionCalc::ShieldingSourceChi2Fcn::create( distance, geometry,
                                                         shield_definitions, src_definitions, detector,
                                                         foreground, background, 
                                                         foreground_peaks, 
                                                         std::make_shared<const std::deque<std::shared_ptr<const PeakDef>>>(background_peaks), 
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
    
  return fitResult.activity;
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
  
  // Get current shielding info to check for self-attenuation compatibility
  bool allowSelfAttenuation = false;
  
  if( peak && peak->parentNuclide() && m_shieldingSelect && !m_shieldingSelect->isNoShielding() )
  {
    try
    {
      ShieldingSourceFitCalc::ShieldingInfo shielding = m_shieldingSelect->toShieldingInfo();
      const SandiaDecay::Nuclide *peakNuclide = peak->parentNuclide();
      
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
