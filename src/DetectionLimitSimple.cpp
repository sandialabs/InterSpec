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
#include <Wt/WSpinBox>
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
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/NuclideSourceEnter.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectionLimitSimple.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;

namespace
{
  bool use_curie_units()
  {
    InterSpec *interspec = InterSpec::instance();
    if( !interspec )
      return true;
    
    return !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", interspec );
  }//bool use_curie_units()
  
}//namespace


class CurrieLimitArea : public Wt::WContainerWidget
{
public:
  DetectionLimitSimple *m_calcTool;
  
public:
  CurrieLimitArea( DetectionLimitSimple *calcTool, Wt::WContainerWidget *parent = 0 )
  : Wt::WContainerWidget( parent ),
  m_calcTool( calcTool )
  {
    assert( m_calcTool );
    addStyleClass( "CurrieLimitArea" );
    
    
    // Result area

  }//CurrieLimitArea constructor
  
  void handleUserRequestedMoreInfoDialog()
  {
    /*
    const DetectionLimitCalc::CurieMdaInput input = currieInput();
    const DetectionLimitCalc::CurieMdaResult result = DetectionLimitCalc::currie_mda_calc( input );
    
    DetectionLimitTool::createCurrieRoiMoreInfoWindow( const SandiaDecay::Nuclide *const nuclide,
     const DetectionLimitCalc::CurieMdaResult &result,
     std::shared_ptr<const DetectorPeakResponse> drf,
     DetectionLimitTool::LimitType limitType,
               const double distance,
               const bool do_air_attenuation,
               const double branch_ratio,
               const double counts_per_bq_into_4pi );
     */
  }//void handleUserRequestedMoreInfoDialog()
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
  m_peakModel( nullptr ),
  m_chartErrMsgStack( nullptr ),
  m_errMsg( nullptr ),
  m_fitFwhmBtn( nullptr ),
  m_nuclideEdit( nullptr ),
  m_nuclideAgeEdit( nullptr ),
  m_nucEnterController( nullptr ),
  m_photoPeakEnergy( nullptr ),
  m_photoPeakEnergiesAndBr{},
  m_distance( nullptr ),
  m_confidenceLevel( nullptr ),
  m_detectorDisplay( nullptr ),
  m_methodGroup( nullptr ),
  m_methodStack( nullptr ),
  m_currieLimitArea( nullptr ),
  m_deconLimitArea( nullptr ),
  m_numFwhmWide( 2.5f ),
  m_lowerRoi( nullptr ),
  m_upperRoi( nullptr ),
  m_numSideChannelLabel( nullptr ),
  m_numSideChannel( nullptr ),
  m_currentNuclide( nullptr ),
  m_currentAge( 0.0 ),
  m_currentEnergy( 0.0 ),
  m_prevDistance{},
  m_stateUri(),
  m_currentCurrieInput( nullptr ),
  m_currentCurrieResults( nullptr ),
  m_currentDeconResults( nullptr )
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
  
  m_spectrum->existingRoiEdgeDragUpdate().connect( boost::bind( &DetectionLimitSimple::roiDraggedCallback, this, _1, _2, _3, _4, _5, _6 ) );
  
  m_chartErrMsgStack->setCurrentIndex( 0 );
  
  m_viewer->displayedSpectrumChanged().connect( this, &DetectionLimitSimple::handleSpectrumChanged );
  
  //shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  //m_spectrum->setData( hist, true );
  //m_spectrum->setXAxisRange( lower_lower_energy - 0.5*dx, upper_upper_energy + 0.5*dx );
  
  m_peakModel = new PeakModel( m_spectrum );
  m_peakModel->setNoSpecMeasBacking();
  m_spectrum->setPeakModel( m_peakModel );
  
  WContainerWidget *generalInput = new WContainerWidget( this );
  generalInput->addStyleClass( "GeneralInput" );
  
  
  WLabel *nucLabel = new WLabel( WString("{1}:").arg(WString::tr("Nuclide")), generalInput );
  m_nuclideEdit = new WLineEdit( generalInput );
  
  m_nuclideEdit->setMinimumSize( 30, WLength::Auto );
  nucLabel->setBuddy( m_nuclideEdit );
  
  WLabel *ageLabel = new WLabel( WString("{1}:").arg(WString::tr("Age")), generalInput );
  m_nuclideAgeEdit = new WLineEdit( generalInput );
  m_nuclideAgeEdit->setMinimumSize( 30, WLength::Auto );
  m_nuclideAgeEdit->setPlaceholderText( WString::tr("N/A") );
  ageLabel->setBuddy( m_nuclideAgeEdit );
  
  nucLabel->addStyleClass( "GridFirstCol GridFirstRow GridVertCenter" );
  m_nuclideEdit->addStyleClass( "GridSecondCol GridFirstRow" );
  ageLabel->addStyleClass( "GridFirstCol GridSecondRow GridVertCenter" );
  m_nuclideAgeEdit->addStyleClass( "GridSecondCol GridSecondRow" );
  
  m_nucEnterController = new NuclideSourceEnterController( m_nuclideEdit, m_nuclideAgeEdit,
                                                          nullptr, this );
  
  m_nucEnterController->changed().connect( this, &DetectionLimitSimple::handleNuclideChanged );
  
  m_viewer->useMessageResourceBundle( "NuclideSourceEnter" );
  HelpSystem::attachToolTipOn( {nucLabel, m_nuclideEdit},
                              WString::tr("dcw-tt-nuc-edit"), showToolTips );
  
  HelpSystem::attachToolTipOn( {ageLabel, m_nuclideAgeEdit},
                              WString::tr("dcw-tt-age-edit"), showToolTips );
  
  
  WLabel *gammaLabel = new WLabel( WString::tr("Gamma"), generalInput );
  gammaLabel->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
  
  
  m_photoPeakEnergy = new WComboBox( generalInput );
  m_photoPeakEnergy->setWidth( WLength(13.5,WLength::FontEm) );
  m_photoPeakEnergy->activated().connect( this, &DetectionLimitSimple::handleGammaChanged );
  m_photoPeakEnergy->addStyleClass( "GridSecondCol GridThirdRow GridVertCenter" );
  
  // TODO: Add FWHM input.  Add text for DRF default, or the button to fit from data.
  //       when user changes this value - dont change ROI limits, just recalc deconv, and redraw either decon or Currie
  
  
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
  m_detectorDisplay->addStyleClass( "GridThirdCol GridThirdRow GridSpanTwoCol GridSpanTwoRows GridVertCenter GridJustifyCenter" );
  m_viewer->detectorChanged().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  m_viewer->detectorModified().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  
  
  
  WLabel *lowerRoiLabel = new WLabel( "ROI Lower:", generalInput );
  lowerRoiLabel->addStyleClass( "GridFirstCol GridFourthRow GridVertCenter" );
  
  m_lowerRoi = new NativeFloatSpinBox( generalInput );
  m_lowerRoi->setSpinnerHidden();
  lowerRoiLabel->setBuddy( m_lowerRoi );
  m_lowerRoi->addStyleClass( "GridSecondCol GridFourthRow" );
  
  WLabel *upperRoiLabel = new WLabel( "ROI Upper:", generalInput );
  upperRoiLabel->addStyleClass( "GridFirstCol GridFifthRow GridVertCenter" );
  
  m_upperRoi = new NativeFloatSpinBox( generalInput );
  m_upperRoi->setSpinnerHidden();
  upperRoiLabel->setBuddy( m_upperRoi );
  m_upperRoi->addStyleClass( "GridSecondCol GridFifthRow" );
  
  m_lowerRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  m_upperRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  
  // Num Side Channel
  m_numSideChannelLabel = new WLabel( "Num Side Chan.:", generalInput );
  m_numSideChannelLabel->addStyleClass( "GridThirdCol GridFifthRow GridVertCenter" );
  m_numSideChannel = new WSpinBox( generalInput );
  m_numSideChannel->addStyleClass( "GridFourthCol GridFifthRow" );
  m_numSideChannel->setRange( 1, 64 );
  m_numSideChannel->setValue( 4 );
  m_numSideChannelLabel->setBuddy( m_numSideChannel );
  m_numSideChannel->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedNumSideChannel );
  
  
  
  WContainerWidget *container = new WContainerWidget( generalInput );
  container->addStyleClass( "MethodSelect GridFirstCol GridSixthRow GridSpanFourCol" );
  
  WLabel *methodLabel = new WLabel( "Calculation Method:", container);
  
  m_methodGroup = new WButtonGroup( container );
  WRadioButton *currieBtn = new Wt::WRadioButton( WString::tr("dls-currie-tab-title"), container );
  m_methodGroup->addButton(currieBtn, 0);
  
  WRadioButton *deconvBtn = new Wt::WRadioButton( WString::tr("dls-decon-tab-title"), container);
  m_methodGroup->addButton(deconvBtn, 1);
  m_methodGroup->setCheckedButton( currieBtn );
  
  m_methodGroup->checkedChanged().connect( boost::bind(&DetectionLimitSimple::handleMethodChanged, this, boost::placeholders::_1) );
  
  
  m_currieLimitArea = new CurrieLimitArea( this );
  m_deconLimitArea = new DeconvolutionLimitArea( this );
  
  
  
  
  //m_methodTabs = new WTabWidget( this );
  //m_methodTabs->addStyleClass( "MethodTabs" );
  //m_methodTabs->addTab( m_currieLimitArea, WString::tr("dls-currie-tab-title") );
  //m_methodTabs->addTab( m_deconLimitArea, WString::tr("dls-decon-tab-title") );
  //m_methodTabs->setCurrentIndex( 0 );
  
  m_methodStack = new WStackedWidget( this );
  m_methodStack->addStyleClass( "CalcMethodStack" );
  
  m_methodStack->addWidget( m_currieLimitArea );
  m_methodStack->addWidget( m_deconLimitArea );
  
  m_methodStack->setCurrentIndex( 0 );
  
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void DetectionLimitSimple::init()



void DetectionLimitSimple::roiDraggedCallback( double new_roi_lower_energy,
                 double new_roi_upper_energy,
                 double new_roi_lower_px,
                 double new_roi_upper_px,
                 double original_roi_lower_energy,
                 bool is_final_range )
{
  if( !is_final_range )
    return;
  
  if( new_roi_upper_energy < new_roi_lower_energy )
    std::swap( new_roi_upper_energy, new_roi_lower_energy );
  
  if( m_currentNuclide
     && ((m_currentEnergy < new_roi_lower_energy) || (m_currentEnergy > new_roi_upper_energy)) )
  {
    passMessage( "Changing the ROI excluded primary gamma - not changing", WarningWidget::WarningMsgHigh );
    return;
  }
  
  m_lowerRoi->setValue( new_roi_lower_energy );
  m_upperRoi->setValue( new_roi_upper_energy );
}//void roiDraggedCallback(...)


void DetectionLimitSimple::handleUserChangedRoi()
{
  // Round to nearest channel edge, swap values if necessary, and make sure stratles the current mean
  
  bool wasValid = true;
  float lower_val = m_lowerRoi->value();
  float upper_val = m_upperRoi->value();
  
  if( lower_val > upper_val )
  {
    wasValid = false;
    std::swap( lower_val, upper_val );
  }
  
  const float meanEnergy = photopeakEnergy();
  const float fwhm = PeakSearchGuiUtils::estimate_FWHM_of_foreground( meanEnergy );
  
  if( meanEnergy > 10.0f )
  {
    if( lower_val >= meanEnergy )
    {
      wasValid = false;
      lower_val = meanEnergy - 1.25f*std::max(1.0f,fwhm);
    }
    
    if( upper_val <= meanEnergy )
    {
      wasValid = false;
      upper_val = meanEnergy + 1.25f*std::max(1.0f,fwhm);
    }
  }//if( meanEnergy > 10.0f )
  
  if( wasValid && (fwhm > 0.1) && (meanEnergy > 10.0f) )
    m_numFwhmWide = (upper_val - lower_val) / fwhm;
  
  const shared_ptr<const SpecUtils::Measurement> hist
                  = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  
  shared_ptr<const SpecUtils::EnergyCalibration> cal = hist ? hist->energy_calibration() : nullptr;
  
  if( cal && cal->valid() && (cal->num_channels() > 7) )
  {
    try
    {
      const float lower_channel = std::round( cal->channel_for_energy( lower_val ) );
      m_lowerRoi->setValue( cal->energy_for_channel( static_cast<int>(lower_channel) ) );
    }catch( std::exception &e )
    {
      cerr << "Error rounding lower ROI energy: " << e.what() << endl;
      assert( 0 );
    }
    
    try
    {
      const float upper_channel = std::round( cal->channel_for_energy( upper_val ) );
      m_upperRoi->setValue( cal->energy_for_channel( static_cast<int>(upper_channel) ) );
    }catch( std::exception &e )
    {
      cerr << "Error rounding upper ROI energy: " << e.what() << endl;
      assert( 0 );
    }
  }//if( valid energy cal )
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//handleUserChangedRoi()


void DetectionLimitSimple::handleUserChangedNumSideChannel()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleUserChangedNumSideChannel()


DetectionLimitSimple::~DetectionLimitSimple()
{
  //nothing to do here
}//~DoseCalcWidget()


void DetectionLimitSimple::handleMethodChanged( Wt::WRadioButton *btn )
{
  m_methodStack->setCurrentIndex( m_methodGroup->checkedId() );
 
  const bool currieMethod = (m_methodGroup->checkedId() == 0);
  
  m_numSideChannelLabel->setHidden( !currieMethod );
  m_numSideChannel->setHidden( !currieMethod );
  
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleMethodChanged( Wt::WRadioButton *btn )


void DetectionLimitSimple::setNuclide( const SandiaDecay::Nuclide *nuc, 
                                      const double age,
                                      const double energy )
{
  if( energy > 10.0 )
    m_currentEnergy = energy;
  
  m_nucEnterController->setNuclideText( nuc ? nuc->symbol : string() );
  m_currentNuclide = nuc;
  
  if( (age >= 0.0) && !m_nucEnterController->nuclideAgeStr().empty() )
  {
    const string agestr = PhysicalUnits::printToBestTimeUnits( age, 5 );
    m_nucEnterController->setNuclideAgeTxt( agestr );
  }//if( age > 0.0 )
  
  handleGammaChanged();
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
  
  // If we dont want to add an undo/redo step, we need to clear this flag, since the
  //  undo/redo step gets added during render, not right now.
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && !undoRedo->canAddUndoRedoNow() )
    m_renderFlags.clear( DetectionLimitSimple::RenderActions::AddUndoRedoStep );
}//void setNuclide( const SandiaDecay::Nuclide *nuc, const double age, const double energy )


float DetectionLimitSimple::photopeakEnergy() const
{
  if( !m_nucEnterController->nuclide() )
    return 0.0f;
  
  const int energyIndex = m_photoPeakEnergy->currentIndex();
  if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    return 0.0f;
  
  return static_cast<float>( m_photoPeakEnergiesAndBr[energyIndex].first );
}//float energy() const


const SandiaDecay::Nuclide *DetectionLimitSimple::nuclide() const
{
  return m_nucEnterController->nuclide();
}//const SandiaDecay::Nuclide *nuclide()


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
  
  if( m_renderFlags.testFlag(RenderActions::UpdateSpectrumDecorations)
     || m_renderFlags.testFlag(RenderActions::UpdateDisplayedSpectrum) 
     || m_renderFlags.testFlag(RenderActions::UpdateLimit) )
  {
    // Needs to be called after updating results
    updateSpectrumDecorations();
  }//if( update displayed spectrum )
  
  
  m_renderFlags = 0;
  
  if( m_stateUri.empty() )
    m_stateUri = encodeStateToUrl();
  
  WContainerWidget::render( flags );
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


void DetectionLimitSimple::handleNuclideChanged()
{
  const SandiaDecay::Nuclide *nuc = m_nucEnterController->nuclide();
  const double age = nuc ? m_nucEnterController->nuclideAge() : 0.0;
  
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
  m_photoPeakEnergiesAndBr.clear();
  
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
  
  double maxYield = 0.0, maxYieldEnergy = 0.0;
  for( const SandiaDecay::EnergyRatePair &erp : photons )
  {
    // Only consider energies above 10 keV
    if( (erp.energy > 10.0) && (erp.numPerSecond >= maxYield) )
    {
      maxYield = erp.numPerSecond;
      maxYieldEnergy = erp.energy;
      if( energyToSelect < 10.0 )
        energyToSelect = erp.energy;
    }
  }//for( const SandiaDecay::EnergyRatePair &erp : photons )
  
  
  
  shared_ptr<SpecMeas> meas = m_viewer->measurment(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const DetectorPeakResponse> det = meas ? meas->detector() : nullptr;
  
  // Using a positive detector resolution sigma will cause us to consider yield when
  //  selecting the energy to choose.  If we haven't changed nuclide, then we'll use
  //  a negative resolution sigma, which will cause us to select the nearest energy,
  //  without considering yield (i.e. if nuclide is same, then don't change energy).
  const double drfSigma = (det && det->hasResolutionInfo() && (energyToSelect > 10.0) && nucChanged)
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
    energyToSelect = maxYieldEnergy;
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
    
    m_photoPeakEnergiesAndBr.push_back( make_pair(energy, intensity) );
    m_photoPeakEnergy->addItem( text );
  }//for each( const double energy, energies )
  
  
  // we wont change `m_currentEnergy`, since the user might go back to their previous nuclide
  
  m_photoPeakEnergy->setCurrentIndex( currentIndex );
  
  handleGammaChanged();
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void handleNuclideChanged()


void DetectionLimitSimple::handleGammaChanged()
{
  if( m_currentNuclide )
  {
    const int gamma_index = m_photoPeakEnergy->currentIndex();
    assert( gamma_index < static_cast<int>(m_photoPeakEnergiesAndBr.size()) );
    
    if(  (gamma_index >= 0) && (gamma_index < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    {
      m_currentEnergy = m_photoPeakEnergiesAndBr[gamma_index].first;
    }else
    {
      //m_currentEnergy = 0.0;
    }
  }//if( m_currentNuclide )
  
  if( m_currentEnergy > 10.0f )
  {
    const float fwhm = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(m_currentEnergy) );
    
    m_lowerRoi->setValue( m_currentEnergy - 0.5*m_numFwhmWide*fwhm );
    m_upperRoi->setValue( m_currentEnergy + 0.5*m_numFwhmWide*fwhm );
    
    handleUserChangedRoi();
  }//if( energy > 10.0f )
  
  //const string current_txt = m_photoPeakEnergy->currentText().toUTF8()
  //const bool is_xray = (current_txt.find("xray") != string::npos);
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
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
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleSpectrumChanged()



void DetectionLimitSimple::updateSpectrumDecorations()
{
  shared_ptr<const ColorTheme> theme = m_viewer->getColorTheme();
  assert( theme );
  
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  
  m_peakModel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
  m_spectrum->removeAllDecorativeHighlightRegions();
  
  if( m_methodStack->currentIndex() == 0 )
  {
    if( m_currentCurrieInput )
    {
      double gammas_per_bq = -1.0, distance = -1.0, br = -1;
      
      try
      {
        distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
      }catch( std::exception & )
      {
      }
      
      const int energyIndex = m_photoPeakEnergy->currentIndex();
      if( (energyIndex >= 0) && (energyIndex < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
      {
        // TODO: sum nearby gammas - maybe something like within 0.5 FWHM and within the ROI... something to think about
        br = m_photoPeakEnergiesAndBr[energyIndex].second;
      }
        
      
      if( drf && drf->isValid() && (distance >= 0.0) && (br > 0) && m_currentCurrieInput->spectrum )
      {
        const bool fixed_geom = drf->isFixedGeometry();
        const float energy = m_currentCurrieInput->gamma_energy;
        const double det_eff = fixed_geom ? drf->intrinsicEfficiency(energy)
                                          : drf->efficiency(energy, distance);
        
        boost::function<double(float)> att_coef_fcn, air_atten_fcn;
        if( distance > 0.0 )
          air_atten_fcn = boost::bind( &GammaInteractionCalc::transmission_coefficient_air, _1, distance );
        
        const double shield_transmission = att_coef_fcn.empty() ? 1.0 : exp( -1.0*att_coef_fcn(energy) );
        const double air_transmission = air_atten_fcn.empty() ? 1.0 : exp( -1.0*air_atten_fcn(energy) );
        const double counts_per_bq_into_4pi = br * shield_transmission * m_currentCurrieInput->spectrum->live_time();
        const double counts_per_bq_into_4pi_with_air = air_transmission * counts_per_bq_into_4pi;
        const double counts_4pi = fixed_geom ? counts_per_bq_into_4pi : counts_per_bq_into_4pi_with_air;

        gammas_per_bq = counts_4pi * det_eff;
      }//if( drf )
      
      const DetectionLimitTool::LimitType limitType = DetectionLimitTool::LimitType::Activity;
      
      DetectionLimitTool::update_spectrum_for_currie_result( m_spectrum, m_peakModel,
              *m_currentCurrieInput, m_currentCurrieResults.get(), drf, limitType, gammas_per_bq );
    }//if( m_currentCurrieInput )
  }else
  {
    if( !m_currentDeconResults )
    {
      
    }else
    {
      
    }//if( no valid result
  }//if( currently doing Currie-style limit ) / else
  
}//void updateSpectrumDecorations()


void DetectionLimitSimple::updateResult()
{
  //m_errMsg->setText( WString::tr("dls-err-no-input") );
  m_errMsg->setText( "" );
  m_currentDeconResults.reset();
  m_currentCurrieResults.reset();
  
  try
  {
    m_fitFwhmBtn->hide();
    
    std::shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !hist || (hist->num_gamma_channels() < 7) )
      throw runtime_error( "No foreground spectrum loaded." );
    
    double energy = 0.0;
    const int energyIndex = m_photoPeakEnergy->currentIndex();
    if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
      energy = 0.5*(m_lowerRoi->value() + m_upperRoi->value());
    else
      energy = m_photoPeakEnergiesAndBr[energyIndex].first;
    
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
    
    if( m_methodStack->currentIndex() == 0 )
    {
      // We need to calculate Currie-style limit
      auto input = make_shared<DetectionLimitCalc::CurieMdaInput>();
      input->spectrum = hist;
      input->gamma_energy = energy;
      input->roi_lower_energy = m_lowerRoi->value();
      input->roi_upper_energy = m_upperRoi->value();
      input->num_lower_side_channels = static_cast<size_t>( m_numSideChannel->value() );
      input->num_upper_side_channels = input->num_lower_side_channels;
      input->detection_probability = confidenceLevel;
      input->additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
      
      m_currentCurrieInput = input;
      DetectionLimitCalc::CurieMdaResult result = DetectionLimitCalc::currie_mda_calc( *input );
      m_currentCurrieResults = make_shared<DetectionLimitCalc::CurieMdaResult>( result );
    }else
    {
      if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
        throw runtime_error( "Please select gamma energy." );
      
      const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
      
      // We need to calculate deconvolution-style limit
      if( !drf )
        throw runtime_error( "Please select a detector efficiency function." );
        
      if( !drf->hasResolutionInfo() )
      {
        m_fitFwhmBtn->show();
        throw runtime_error( "DRF does not have FWHM info - please fit for FWHM, or change DRF." );
      }
      
      
      if( m_distance->validate() != WValidator::State::Valid )
        throw runtime_error( "Invalid distance" );
      
      double distance = 0.0;
      try
      {
        distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
      }catch( std::exception & )
      {
        throw runtime_error( "invalid distance." );
      }
      
      if( distance < 0.0 )
        throw runtime_error( "Distance can't be negative." );
      
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
    
      m_currentDeconResults = make_shared<DetectionLimitCalc::DeconComputeResults>( results );
    }//if( calc Currie-style limit ) / else
    
    m_chartErrMsgStack->setCurrentIndex( 1 );
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
   m_currentNuclide = m_nucEnterController->nuclide();
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


