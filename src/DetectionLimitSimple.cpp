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


#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WMenuItem>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WTabWidget>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WButtonGroup>
#include <Wt/WRadioButton>
#include <Wt/WStackedWidget>
#include <Wt/WDoubleValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/NuclideSourceEnter.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectionLimitSimple.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;


class CurrieLimitArea : public Wt::WContainerWidget
{
  DetectionLimitSimple *m_calcTool;
  
public:
  CurrieLimitArea( DetectionLimitSimple *calcTool, Wt::WContainerWidget *parent = 0 )
  : Wt::WContainerWidget( parent ),
  m_calcTool( calcTool )
  {
    assert( m_calcTool );
    
    new WText( "Currie Limit Area", this );
    
    // ROI Lower
    // ROI Upper
    // Num Side Channel
    
    // Result area
  }
  
  
};//class CurrieLimitArea


class DeconvolutionLimitArea : public Wt::WContainerWidget
{
  DetectionLimitSimple *m_calcTool;
  
public:
  DeconvolutionLimitArea( DetectionLimitSimple *calcTool, Wt::WContainerWidget *parent = 0 )
  : Wt::WContainerWidget( parent ),
  m_calcTool( calcTool )
  {
    assert( m_calcTool );
    
    new WText( "Deconvolution Limit Area", this );
  }
  

};//class DeconvolutionLimitArea


DetectionLimitSimpleWindow::DetectionLimitSimpleWindow( MaterialDB *materialDB,
                                Wt::WSuggestionPopup *materialSuggestion,
                                InterSpec *viewer )
: AuxWindow( WString::tr("window-title-simple-mda"),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse) )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  rejectWhenEscapePressed( true );
  
  m_tool = new DetectionLimitSimple( materialDB, materialSuggestion, viewer, contents() );
  m_tool->setHeight( WLength(100,WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "simple-mda-dialog" );
  
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://simple-mda/" + Wt::Utils::urlEncode(m_tool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("dlsw-qr-tool-state-title"),
                                 WString::tr("dlsw-qr-tool-state-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES

  
  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  show();
  
  // If we are loading this widget, as we  are creating the InterSpec session,
  //  the screen width and height wont be available, so we'll just assume its
  //  big enough, which it should be.
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  int width = 625, height = 435;
  if( (screenW > 100) && (screenW < width) )
    width = screenW;
  if( (screenH > 100) && (screenH < height) )
    height = screenH;
  
  resizeWindow( width, height );
  
  // But I think this next call should fix things up, even if we do have a tiny screen
  resizeToFitOnScreen();
  
  centerWindowHeavyHanded();
}//DetectionLimitSimpleWindow(...) constructor


DetectionLimitSimpleWindow::~DetectionLimitSimpleWindow()
{
}


DetectionLimitSimple *DetectionLimitSimpleWindow::tool()
{
  return m_tool;
}





DetectionLimitSimple::DetectionLimitSimple( MaterialDB *materialDB,
                                 Wt::WSuggestionPopup *materialSuggestion,
                                 InterSpec *specViewer,
                                 Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_viewer( specViewer ),
  m_materialSuggest( materialSuggestion ),
  m_materialDB( materialDB ),
  m_spectrum( nullptr ),
  m_chartErrMsgStack( nullptr ),
  m_errMsg( nullptr ),
  m_fitFwhmBtn( nullptr ),
  m_nuclideEnter( nullptr ),
  m_photoPeakEnergy( nullptr ),
  m_photoPeakEnergies{},
  m_distance( nullptr ),
  m_confidenceLevel( nullptr ),
  m_detectorDisplay( nullptr ),
  m_methodTabs( nullptr ),
  m_currieLimitArea( nullptr ),
  m_deconLimitArea( nullptr ),
  m_currentNuclide( nullptr ),
  m_currentAge( 0.0 ),
  m_currentEnergy( 0.0 ),
  m_prevDistance{},
  m_stateUri()
{
  init();
}//DoseCalcWidget constructor


void DetectionLimitSimple::init()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  wApp->useStyleSheet( "InterSpec_resources/DetectionLimitSimple.css" );
  m_viewer->useMessageResourceBundle( "DetectionLimitSimple" );
      
  addStyleClass( "DetectionLimitSimple" );
 
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_viewer );
  
  m_chartErrMsgStack = new WStackedWidget( this );
  
  WContainerWidget *errorDiv = new WContainerWidget();
  errorDiv->addStyleClass( "ErrDisplay" );
  m_chartErrMsgStack->addWidget( errorDiv );
  
  m_errMsg = new WText( WString::tr("dls-err-no-input"), errorDiv );
  m_errMsg->addStyleClass( "ErrMsg" );
  
  m_fitFwhmBtn = new WPushButton( "Fit FWHM...", errorDiv );
  m_fitFwhmBtn->addStyleClass( "MdaFitFwhm LightButton" );
  m_fitFwhmBtn->clicked().connect( this, &DetectionLimitSimple::handleFitFwhmRequested );
  m_fitFwhmBtn->hide();
  
  
  m_spectrum = new D3SpectrumDisplayDiv();
  m_chartErrMsgStack->addWidget( m_spectrum );
  m_spectrum->setXAxisTitle( "" );
  m_spectrum->setYAxisTitle( "", "" );
  m_spectrum->setYAxisLog( false );
  m_spectrum->applyColorTheme( m_viewer->getColorTheme() );
  m_viewer->colorThemeChanged().connect( boost::bind( &D3SpectrumDisplayDiv::applyColorTheme, m_spectrum, boost::placeholders::_1 ) );
  m_spectrum->disableLegend();
  m_spectrum->setShowPeakLabel( SpectrumChart::PeakLabels::kShowPeakUserLabel, true );
  
  m_chartErrMsgStack->setCurrentIndex( 0 );
  
  m_viewer->displayedSpectrumChanged().connect( this, &DetectionLimitSimple::handleSpectrumChanged );
  
  //shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  //m_spectrum->setData( hist, true );
  //m_spectrum->setXAxisRange( lower_lower_energy - 0.5*dx, upper_upper_energy + 0.5*dx );
  
  WContainerWidget *generalInput = new WContainerWidget( this );
  generalInput->addStyleClass( "GeneralInput" );
  
  m_nuclideEnter = new NuclideSourceEnter( false, showToolTips, generalInput );
  m_nuclideEnter->changed().connect( this, &DetectionLimitSimple::handleNuclideChanged );
  m_nuclideEnter->addStyleClass( "GridFirstCol GridFirstRow GridSpanTwoCol GridSpanTwoRows" );
  
  WLabel *gammaLabel = new WLabel( WString::tr("Gamma"), generalInput );
  gammaLabel->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
  
  
  m_photoPeakEnergy = new WComboBox( generalInput );
  m_photoPeakEnergy->setWidth( WLength(13.5,WLength::FontEm) );
  m_photoPeakEnergy->activated().connect( this, &DetectionLimitSimple::handleGammaChanged );
  m_photoPeakEnergy->addStyleClass( "GridSecondCol GridThirdRow GridVertCenter" );
  
  m_methodTabs = new WTabWidget( this );
  m_methodTabs->addStyleClass( "MethodTabs" );
  
  // Add Distance input
  WLabel *distanceLabel = new WLabel( WString::tr("Distance"), generalInput );
  distanceLabel->addStyleClass( "GridThirdCol GridFirstRow GridVertCenter" );
  
  m_prevDistance = "100 cm";
  m_distance = new WLineEdit( m_prevDistance, generalInput );
  m_distance->addStyleClass( "GridFourthCol GridFirstRow GridStretchCol" );
  distanceLabel->setBuddy( m_distance );
  
  m_distance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distance->setAttributeValue( "autocorrect", "off" );
  m_distance->setAttributeValue( "spellcheck", "off" );
#endif
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_distance->setValidator( validator );
  //HelpSystem::attachToolTipOn( m_distance, WString::tr("ftw-tt-distance"), showToolTips );
  m_distance->changed().connect( this, &DetectionLimitSimple::handleDistanceChanged );
  m_distance->enterPressed().connect( this, &DetectionLimitSimple::handleDistanceChanged );
  
  
  // Add confidence select
  WLabel *confidenceLabel = new WLabel( "Confidence Level:", generalInput );
  confidenceLabel->addStyleClass( "GridThirdCol GridSecondRow GridVertCenter" );
  m_confidenceLevel = new WComboBox( generalInput );
  m_confidenceLevel->addStyleClass( "GridFourthCol GridSecondRow" );
  
  for( auto cl = ConfidenceLevel(0); cl < NumConfidenceLevel; cl = ConfidenceLevel(cl+1) )
  {
    const char *txt = "";
    
    switch( cl )
    {
      case ConfidenceLevel::OneSigma:   txt = "68%";     break;
      case ConfidenceLevel::TwoSigma:   txt = "95%";     break;
      case ConfidenceLevel::ThreeSigma: txt = "99%";     break;
      case ConfidenceLevel::FourSigma:  txt = "4-sigma"; break;
      case ConfidenceLevel::FiveSigma:  txt = "5-sigma"; break;
      case ConfidenceLevel::NumConfidenceLevel:          break;
    }//switch( cl )
    
    m_confidenceLevel->addItem( txt );
  }//for( loop over confidence levels )
  
  m_confidenceLevel->setCurrentIndex( ConfidenceLevel::TwoSigma );
  m_confidenceLevel->activated().connect(this, &DetectionLimitSimple::handleConfidenceLevelChanged );
  
  
  // Add DRF select
  SpectraFileModel *specFileModel = m_viewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( m_viewer, specFileModel, generalInput );
  m_detectorDisplay->addStyleClass( "GridThirdCol GridThirdRow GridSpanTwoCol GridSpanTwoRows" );
  m_viewer->detectorChanged().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  m_viewer->detectorModified().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  
  
  m_currieLimitArea = new CurrieLimitArea( this );
  m_deconLimitArea = new DeconvolutionLimitArea( this );
  
  m_methodTabs->addTab( m_currieLimitArea, WString::tr("dls-currie-tab-title") );
  m_methodTabs->addTab( m_deconLimitArea, WString::tr("dls-decon-tab-title") );
  
  m_methodTabs->setCurrentIndex( 0 );
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void DetectionLimitSimple::init()


DetectionLimitSimple::~DetectionLimitSimple()
{
  //nothing to do here
}//~DoseCalcWidget()


void DetectionLimitSimple::setNuclide( const SandiaDecay::Nuclide *nuc, 
                                      const double age,
                                      const double energy )
{
  m_nuclideEnter->setNuclideText( nuc ? nuc->symbol : string() );
  if( age >= 0.0 )
  {
    const string agestr = PhysicalUnits::printToBestTimeUnits( age, 4 );
    m_nuclideEnter->setNuclideAgeTxt( agestr );
  }//if( age > 0.0 )
  
  m_currentEnergy = energy;
  
  handleNuclideChanged();
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
  
  // If we dont want to add an undo/redo step, we need to clear this flag, since the
  //  undo/redo step gets added during render, not right now.
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->canAddUndoRedoNow() )
    m_renderFlags.clear( DetectionLimitSimple::RenderActions::AddUndoRedoStep );
}//void setNuclide( const SandiaDecay::Nuclide *nuc, const double age, const double energy )


void DetectionLimitSimple::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  if( m_renderFlags.testFlag(RenderActions::UpdateDisplayedSpectrum) )
  {
    shared_ptr<const SpecUtils::Measurement> hist 
                                = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    m_spectrum->setData( hist, true );
  }//if( update displayed spectrum )
  
  if( m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) )
  {
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    
    if( undoRedo )
    {
      string uri = encodeStateToUrl();
      const bool sameAsPrev = (uri == m_stateUri);
      
      if( !m_stateUri.empty() && undoRedo->canAddUndoRedoNow() && !sameAsPrev )
      {
        const shared_ptr<const string> prev = make_shared<string>( std::move(m_stateUri) );
        const shared_ptr<const string> current = make_shared<string>( uri );
        
        auto undo_redo = [prev, current]( bool is_undo ){
          DetectionLimitSimpleWindow *mdawin = InterSpec::instance()->showSimpleMdaWindow();
          DetectionLimitSimple *tool = mdawin ? mdawin->tool() : nullptr;
          const string &uri = is_undo ? *prev : *current;
          if( tool && !uri.empty() )
            tool->handleAppUrl( uri );
        };//undo_redo
        
        auto undo = [undo_redo](){ undo_redo(true); };
        auto redo = [undo_redo](){ undo_redo(false); };
        
        undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Update Simple MDA values." );
      }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
      
      m_stateUri = std::move(uri);
    }//if( undoRedo )
    
  }//if( m_renderFlags.testFlag(RenderActions::AddUndoRedoStep) )
  
  if( m_renderFlags.testFlag(RenderActions::UpdateLimit) )
    updateResult();
  
  m_renderFlags = 0;
  
  if( m_stateUri.empty() )
    m_stateUri = encodeStateToUrl();
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void DetectionLimitSimple::handleNuclideChanged()
{
  const SandiaDecay::Nuclide *nuc = m_nuclideEnter->nuclide();
  const double age = nuc ? m_nuclideEnter->nuclideAge() : 0.0;
  
  const bool nucChanged = (nuc != m_currentNuclide);
  const bool ageChanged = (m_currentAge != age);
  
  if( !nucChanged && !ageChanged )
  {
    cout << "DetectionLimitSimple::handleNuclideChanged(): Nuclide not actually changed - not doing anything." << endl;
    return;
  }//if( !nucChanged && !ageChanged )
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
  
  m_photoPeakEnergy->clear();
  m_photoPeakEnergies.clear();
  
  m_currentNuclide = nuc;
  m_currentAge = age;
    
  if( !nuc )
    return;
  
  const bool dummy_activity = 0.001*SandiaDecay::curie;
  SandiaDecay::NuclideMixture mixture;
  mixture.addAgedNuclideByActivity( nuc, dummy_activity, age );
  
  const vector<SandiaDecay::EnergyRatePair> xrays = mixture.xrays( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
  const vector<SandiaDecay::EnergyRatePair> gammas = mixture.gammas( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
  
  vector<SandiaDecay::EnergyRatePair> photons;
  photons.insert( end(photons), begin(xrays), end(xrays) );
  photons.insert( end(photons), begin(gammas), end(gammas) );
  
  
  double energyToSelect = m_currentEnergy;
  
  // If we dont currently have an energy selected, pick the largest yield energy
  if( energyToSelect < 10.0 )
  {
    double maxYield = 0.0;
    for( const SandiaDecay::EnergyRatePair &erp : photons )
    {
      // Only consider energies above 10 keV
      if( (erp.energy > 10.0) && (erp.numPerSecond >= maxYield) )
      {
        maxYield = erp.numPerSecond;
        energyToSelect = erp.energy;
      }
    }//for( const SandiaDecay::EnergyRatePair &erp : photons )
  }//if( energyToSelect < 10.0 )
  
  
  shared_ptr<SpecMeas> meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const DetectorPeakResponse> det = meas ? meas->detector() : nullptr;
  
  // Using a positive detector resolution sigma will cause us to consider yield when
  //  selecting the energy to choose.  If we haven't changed nuclide, then we'll use
  //  a negative resolution sigma, which will cause us to select the nearest energy,
  //  without considering yield (i.e. if nuclide is same, then don't change energy).
  const double drfSigma = (det && (energyToSelect > 10.0) && nucChanged)
                          ? det->peakResolutionSigma( static_cast<float>(energyToSelect) )
                          : -1.0;
  
  size_t transition_index = 0;
  const SandiaDecay::Transition *transition = nullptr;
  PeakDef::SourceGammaType sourceGammaType = PeakDef::SourceGammaType::NormalGamma;
  PeakDef::findNearestPhotopeak( nuc, energyToSelect, 4.0*drfSigma, false,
                                transition, transition_index, sourceGammaType );
  if( transition && (transition_index < transition->products.size()) )
  {
    energyToSelect = transition->products[transition_index].energy;
  }else if( sourceGammaType == PeakDef::AnnihilationGamma )
  {
    energyToSelect = 510.998910;
  }else
  {
    assert( 0 );
  }
  
  
  double minDE = DBL_MAX;
  int currentIndex = -1;
   
  for( size_t i = 0; i < photons.size(); ++i )
  {
    double energy = photons[i].energy;
    const double intensity = photons[i].numPerSecond / dummy_activity;
    
    const double dE = fabs( energyToSelect - energy );
    energy = floor(10000.0*energy + 0.5)/10000.0;
    
    if( dE < minDE )
    {
      minDE = dE;
      currentIndex = static_cast<int>( i );
    }
    
    char text[128] = { '\0' };
    if( i < xrays.size() )
    {
      snprintf( text, sizeof(text), "%.4f keV xray I=%.1e", energy, intensity );
    }else
    {
      snprintf( text, sizeof(text), "%.4f keV I=%.1e", energy, intensity );
    }//if( i < xrays.size() )
    
    m_photoPeakEnergies.push_back( energy );
    m_photoPeakEnergy->addItem( text );
  }//for each( const double energy, energies )
  
  
  // we wont change `m_currentEnergy`, since the user might go back to their previous nuclide
  
  m_photoPeakEnergy->setCurrentIndex( currentIndex );
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void handleNuclideChanged()


void DetectionLimitSimple::handleGammaChanged()
{
  if( !m_currentNuclide )
    return;
  
  const int gamma_index = m_photoPeakEnergy->currentIndex();
  
  assert( gamma_index < static_cast<int>(m_photoPeakEnergies.size()) );
  
  if( (gamma_index < 0) || (gamma_index >= static_cast<int>(m_photoPeakEnergies.size())) )
  {
    m_currentEnergy = m_photoPeakEnergies[gamma_index];
  }else
  {
    m_currentEnergy = 0.0;
  }
  
  //const string current_txt = m_photoPeakEnergy->currentText().toUTF8()
  //const bool is_xray = (current_txt.find("xray") != string::npos);
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void handleGammaChanged()


void DetectionLimitSimple::handleDistanceChanged()
{
  WString dist = m_distance->text();
  const WString prev = m_prevDistance;
  
  if( dist == prev )
    return;
  
  try
  {
    if( m_distance->validate() != WValidator::State::Valid )
      throw runtime_error( "Invalid distance" );
    
    PhysicalUnits::stringToDistance( dist.toUTF8() );
    
    m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
    scheduleRender();
  }catch( std::exception & )
  {
    m_distance->setText( prev );
    dist = prev;
  }//try / catch
  
  m_prevDistance = dist;
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleDistanceChanged()


void DetectionLimitSimple::handleConfidenceLevelChanged()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleConfidenceLevelChanged()


void DetectionLimitSimple::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf )
{
  // The DetectorDisplay::setDetector(...) function may not have been called yet because of
  //  the order of signal/slot connections - so we'll set that here to make sure we are up
  //  to date.
  m_detectorDisplay->setDetector( new_drf );
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> new_drf )


void DetectionLimitSimple::handleFitFwhmRequested()
{
  const bool use_auto_fit_peaks_too = true;
  MakeFwhmForDrfWindow *window = m_viewer->fwhmFromForegroundWindow( use_auto_fit_peaks_too );
  if( window )
  {
    // probably nothing to do here
  }
}//void handleFitFwhmRequested()


void DetectionLimitSimple::handleSpectrumChanged()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleSpectrumChanged()


void DetectionLimitSimple::updateResult()
{
  //m_errMsg->setText( WString::tr("dls-err-no-input") );
  m_errMsg->setText( "" );
  
  try
  {
    m_fitFwhmBtn->hide();
    
    if( m_distance->validate() != WValidator::State::Valid )
      throw runtime_error( "Invalid distance" );
    
    std::shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !hist || (hist->num_gamma_channels() < 7) )
      throw runtime_error( "No foreground spectrum loaded." );
    
    const SandiaDecay::Nuclide *nuc = m_nuclideEnter->nuclide();
    if( !nuc )
      throw runtime_error( "Please enter a nuclide." );
    
    const int energyIndex = m_photoPeakEnergy->currentIndex();
    if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergies.size())) )
      throw runtime_error( "Please select gamma energy." );
    
    const double energy = m_photoPeakEnergies[energyIndex];
    
    
    const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
    
    double distance = 0.0;
    
    try
    {
      const string dist_str = m_distance->text().toUTF8();
      distance = PhysicalUnits::stringToDistance( dist_str );
    }catch( std::exception & )
    {
      throw runtime_error( "invalid distance." );
    }
    
    if( distance <= 0.0 )
      throw runtime_error( "Distance cant be zero or negative." );
    
    const int clIndex = m_confidenceLevel->currentIndex();
    if( (clIndex < 0) || (clIndex >= ConfidenceLevel::NumConfidenceLevel) )
      throw runtime_error( "Please select confidence level." );
      
    float confidenceLevel = 0.95;
    const ConfidenceLevel confidence = ConfidenceLevel(clIndex);
    switch( confidence )
    {
      case OneSigma:   confidenceLevel = 0.682689492137086f; break;
      case TwoSigma:   confidenceLevel = 0.954499736103642f; break;
      case ThreeSigma: confidenceLevel = 0.997300203936740f; break;
      case FourSigma:  confidenceLevel = 0.999936657516334f; break;
      case FiveSigma:  confidenceLevel = 0.999999426696856f; break;
      case NumConfidenceLevel: assert(0); break;
    }//switch( confidence )
    
    
    if( m_methodTabs->currentIndex() == 0 )
    {
      throw runtime_error( "Currie-style calculations not implemented yet" );
      
      // We need to calculate Currie-style limit
      DetectionLimitCalc::CurieMdaInput input;
      input.spectrum = hist;
      input.gamma_energy = energy;
      //input.roi_lower_energy = ...;
      //input.roi_upper_energy = ...;
      //input.num_lower_side_channels = ...;
      //input.num_upper_side_channels = ...;
      input.detection_probability = confidenceLevel;
      input.additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
      
      DetectionLimitCalc::CurieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
    }else
    {
      // We need to calculate deconvolution-style limit
      if( !drf )
        throw runtime_error( "Please select a detector efficiency function." );
        
      if( !drf->hasResolutionInfo() )
      {
        m_fitFwhmBtn->show();
        throw runtime_error( "DRF does not have FWHM info - please fit for FWHM, or change DRF." );
      }
      
      throw runtime_error( "Deconvolution calculations not implemented yet" );
      
      DetectionLimitCalc::DeconRoiInfo roiInfo;
      //roiInfo.roi_start = ...; // Will be rounded to nearest channel edge.
      //roiInfo.roi_end = ...;// Will be rounded to nearest channel edge.
      roiInfo.continuum_type = PeakContinuum::OffsetType::Linear;
        
      //Whether to allow the continuum to float in the fit, or to fix the continuum using the peaks bordering the ROI, or use the whole ROI to determine the continuum with the assumption no signal is there.
      //roiInfo.cont_norm_method = DeconContinuumNorm::;
      //roiInfo.num_lower_side_channels = ; //Only used if `cont_norm_method` is `DeconContinuumNorm::FixedByEdges`.
      //roiInfo.num_upper_side_channels = ; //Only used if `cont_norm_method` is `DeconContinuumNorm::FixedByEdges`.
      
      //DetectionLimitCalc::DeconRoiInfo::PeakInfo peakInfo;
      //peakInfo.energy = ;
      //peakInfo.fwhm = ;
      //peakInfo.counts_per_bq_into_4pi = ;//must have effects of shielding already accounted for, but not air atten, or det intrinsic eff
       
      //roiInfo.peak_infos.push_back( peakInfo );
      
      // TODO: put additional peaks...
      
      DetectionLimitCalc::DeconComputeInput input;
      input.distance = distance;
      //dont need to fill out: input.activity;
      input.include_air_attenuation = true;
      input.shielding_thickness = 0.0;
      input.measurement = hist;
      input.drf = drf;
        
      input.roi_info.push_back( roiInfo );
    
      DetectionLimitCalc::DeconComputeResults results = DetectionLimitCalc::decon_compute_peaks( input );
    }//if( calc Currie-style limit ) / else
    
    
    
  }catch( std::exception &e )
  {
    m_chartErrMsgStack->setCurrentIndex( 0 );
    m_errMsg->setText( WString::tr("dls-err-calculating").arg(e.what()) );
  }//try / catch
}//void updateResult()


void DetectionLimitSimple::handleAppUrl( std::string uri )
{
  //blah blah blah handle all this
  /*
#if( PERFORM_DEVELOPER_CHECKS )
  const string expected_uri = path + "?" + query_str;
#endif
  
   m_currentEnergy
   
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  Quantity calcQuantity;
  if( SpecUtils::iequals_ascii(path, "dose") )
    calcQuantity = Quantity::Dose;
  else if( SpecUtils::iequals_ascii(path, "act") )
    calcQuantity = Quantity::Activity;
  else if( SpecUtils::iequals_ascii(path, "dist") )
    calcQuantity = Quantity::Distance;
  else if( SpecUtils::iequals_ascii(path, "shield") )
    calcQuantity = Quantity::Shielding;
  else if( SpecUtils::iequals_ascii(path, "intro") )
  {
    handleQuantityClick( Quantity::NumQuantity );
    return;
  }else
    throw runtime_error( "Dose Calc tool: invalid URI path." );
  
    
  string inshielduri, outshielduri;
  const size_t shieldpos = query_str.find( "&INSHIELD=&" );
  
  if( shieldpos != string::npos )
  {
    inshielduri = query_str.substr(shieldpos + 11);
    query_str = query_str.substr(0, shieldpos);
  }
  
  size_t outshieldpos = inshielduri.find("&OUTSHIELD=&");
  if( outshieldpos != string::npos )
  {
    outshielduri = inshielduri.substr(outshieldpos + 12);
    inshielduri = inshielduri.substr(0, outshieldpos);
  }else if( (outshieldpos = query_str.find("&OUTSHIELD=&")) != string::npos )
  {
    outshielduri = query_str.substr(outshieldpos + 12);
    query_str = query_str.substr(0, outshieldpos);
  }
  
  if( inshielduri.empty() )
    m_enterShieldingSelect->setToNoShielding();
  else
    m_enterShieldingSelect->handleAppUrl( inshielduri );
  
  if( outshielduri.empty() )
    m_answerShieldingSelect->setToNoShielding();
  else
    m_answerShieldingSelect->handleAppUrl( outshielduri );
  
  SpecUtils::ireplace_all( query_str, "%23", "#" );
  SpecUtils::ireplace_all( query_str, "%26", "&" );
  SpecUtils::ireplace_all( query_str, "curries", "curies" ); //fix up me being a bad speller
  
  const map<string,string> parts = AppUtils::query_str_key_values( query_str );
  const auto ver_iter = parts.find( "VER" );
  if( ver_iter == end(parts) )
    Wt::log("warn") << "No 'VER' field in Dose Calc tool URI.";
  else if( ver_iter->second != "1" && !SpecUtils::starts_with(ver_iter->second, "1.") )
    throw runtime_error( "Can not read Dose Calc tool URI version '" + ver_iter->second + "'" );
  
  auto findUnitIndex = [&parts]( const string &key, Wt::WComboBox *combo ) -> int {
    const auto act_unit_iter = parts.find(key);
    if( act_unit_iter == end(parts) )
      return -1;
    
    for( int i = 0; i < combo->count(); ++i )
    {
      if( SpecUtils::iequals_ascii( combo->itemText(i).toUTF8(), act_unit_iter->second ) )
        return i;
    }
    assert( 0 );
    return -1;
  };//findUnitIndex
  
  const int act_in_unit_index = findUnitIndex( "ACTINUNIT", m_activityEnterUnits );
  const int act_out_unit_index = findUnitIndex( "ACTOUTUNIT", m_activityAnswerUnits );
  const int dose_in_unit_index = findUnitIndex( "DOSEINUNIT", m_doseEnterUnits );
  const int dose_out_unit_index = findUnitIndex( "DOSEOUTUNIT", m_doseAnswerUnits );
  
  const auto act_iter = parts.find("ACT");
  if( act_iter != end(parts) && !act_iter->second.empty() && (act_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain activity unit info." );
  
  const auto dose_iter = parts.find("DOSE");
  if( dose_iter != end(parts) && !dose_iter->second.empty() && (dose_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain dose unit info." );
  
  m_menu->select( static_cast<int>(calcQuantity) );
  handleQuantityClick( calcQuantity );
  
  if( act_in_unit_index >= 0 )
    m_activityEnterUnits->setCurrentIndex( act_in_unit_index );
  if( act_out_unit_index >= 0 )
    m_activityAnswerUnits->setCurrentIndex( act_out_unit_index );
  if( dose_in_unit_index >= 0 )
    m_doseEnterUnits->setCurrentIndex( dose_in_unit_index );
  if( dose_out_unit_index >= 0 )
    m_doseAnswerUnits->setCurrentIndex( dose_out_unit_index );
  
  const auto nuc_iter = parts.find("NUC");
  if( nuc_iter != end(parts) )
    m_gammaSource->setNuclideText( nuc_iter->second );
  
  const auto age_iter = parts.find("AGE");
  if( age_iter != end(parts) )
    m_gammaSource->setNuclideAgeTxt( age_iter->second );
  
  if( act_iter != end(parts) )
    m_activityEnter->setText( WString::fromUTF8(act_iter->second) );
  else
    m_activityEnter->setText( "" );
  
  if( dose_iter != end(parts) )
    m_doseEnter->setText( WString::fromUTF8(dose_iter->second) );
  else
    m_doseEnter->setText( "" );
  
  const auto dist_iter = parts.find("DIST");
  if( dist_iter != end(parts) )
    m_distanceEnter->setText( WString::fromUTF8(dist_iter->second) );
  else
    m_distanceEnter->setText( "" );
  
  
   updateResult();
   m_stateUri = encodeStateToUrl();
   m_currentNuclide = m_nuclideEnter->nuclide();
   m_currentEnergy = ;
   m_currentAge = ;
   
  
#if( PERFORM_DEVELOPER_CHECKS )
  if( m_stateUri != expected_uri )
  {
    Wt::log("warn") << "DoseCalcWidget::handleAppUrl: input URI doesnt match current URI.\n\t input: '"
                    << expected_uri.c_str() << "'\n\tresult: '" << m_stateUri.c_str() << "'";
  }
#endif
   
   */
}//void handleAppUrl( std::string uri )


std::string DetectionLimitSimple::encodeStateToUrl() const
{
  /*
  // "interspec://dose/act?nuc=u238&dose=1.1ur/h&dist=100cm&..."
  
  string answer;
  
  switch( m_currentCalcQuantity )
  {
    case Dose:        answer += "dose";   break;
    case Activity:    answer += "act";    break;
    case Distance:    answer += "dist";   break;
    case Shielding:   answer += "shield"; break;
    case NumQuantity: answer += "intro"; break;
  }//switch( m_currentCalcQuantity )
  
  answer += "?VER=1";

  if( m_currentCalcQuantity == NumQuantity )
    return answer;
  
  // We could limit what info we put in the URL, based on current m_currentCalcQuantity,
  //  but we might as well put the full state into the URI.
  
  auto addTxtField = [&answer]( const string &key, Wt::WLineEdit *edit ){
    string txt = edit->text().toUTF8();
    SpecUtils::ireplace_all( txt, "#", "%23" );
    SpecUtils::ireplace_all( txt, "&", "%26" );
    answer += "&" + key + "=" + txt;
  };

  const SandiaDecay::Nuclide *nuc = m_gammaSource->nuclide();
  if( nuc )
    answer += "&NUC=" + nuc->symbol;
  
  if( nuc && !PeakDef::ageFitNotAllowed(nuc) )
    answer += "&AGE=" + m_gammaSource->nuclideAgeStr().toUTF8();
  
  addTxtField( "ACT", m_activityEnter );
  answer += "&ACTINUNIT=" + m_activityEnterUnits->currentText().toUTF8();
  answer += "&ACTOUTUNIT=" + m_activityAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DOSE", m_doseEnter );
  answer += "&DOSEINUNIT=" + m_doseEnterUnits->currentText().toUTF8();
  answer += "&DOSEOUTUNIT=" + m_doseAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DIST", m_distanceEnter );
  
  // We'll mark shielding URL starting with the below, and everything after this is the shielding.
  //  Having an order-independent method would be better, but for the moment...
  switch( m_currentCalcQuantity )
  {
    case Dose:
    case Activity:
    case Distance:
      if( !m_enterShieldingSelect->isNoShielding() )
        answer += "&INSHIELD=&" + m_enterShieldingSelect->encodeStateToUrl();
      break;
      
    case Shielding:
      if( !m_answerShieldingSelect->isNoShielding() )
        answer += "&OUTSHIELD=&" + m_answerShieldingSelect->encodeStateToUrl();
      break;
      
    case NumQuantity:
      break;
  }//switch( m_currentCalcQuantity )
  
  
  
  return answer;
   */
  //blah blah blah handle all this
  return "";
}//std::string encodeStateToUrl() const;


