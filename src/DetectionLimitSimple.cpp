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
#include "InterSpec/SimpleDialog.h"
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
      QrCode::displayTxtAsQrCode( url, WString::tr("dls-qr-tool-state-title"),
                                 WString::tr("dls-qr-tool-state-txt") );
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
  int width = 525, height = 800;
  if( (screenW > 100) && (screenW < width) )
    width = screenW;
  if( (screenH > 100) && (screenH < height) )
    height = screenH;
  
  //resizeWindow( width, height );
  setWidth( width );
  setMaximumSize( WLength::Auto, height );
  
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
  m_resultTxt( nullptr ),
  m_moreInfoButton( nullptr ),
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
  m_methodDescription( nullptr ),
  m_numFwhmWide( 2.5f ),
  m_lowerRoi( nullptr ),
  m_upperRoi( nullptr ),
  m_numSideChannelLabel( nullptr ),
  m_numSideChannel( nullptr ),
  m_fwhm( nullptr ),
  m_fwhmSuggestTxt( nullptr ),
  m_addFwhmBtn( nullptr ),
  m_selectDetectorBtn( nullptr ),
  m_continuumPriorLabel( nullptr ),
  m_continuumPrior( nullptr ),
  m_continuumTypeLabel( nullptr ),
  m_moreInfoWindow( nullptr ),
  m_continuumType( nullptr ),
  m_currentNuclide( nullptr ),
  m_currentAge( 0.0 ),
  m_currentEnergy( 0.0 ),
  m_allGammasInRoi( true ),
  m_prevDistance{},
  m_stateUri(),
  m_currentCurrieInput( nullptr ),
  m_currentCurrieResults( nullptr ),
  m_currentDeconInput( nullptr ),
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
  
  WContainerWidget *resultsDiv = new WContainerWidget( this );
  resultsDiv->addStyleClass( "ResultsArea" );
  
  m_chartErrMsgStack = new WStackedWidget( resultsDiv );
  
  WContainerWidget *errorDiv = new WContainerWidget();
  errorDiv->addStyleClass( "ErrDisplay" );
  m_chartErrMsgStack->addWidget( errorDiv );
  
  m_errMsg = new WText( WString::tr("dls-err-no-input"), errorDiv );
  m_errMsg->addStyleClass( "ErrMsg" );
  
  m_fitFwhmBtn = new WPushButton( WString::tr("dls-fit-fwhm-btn"), errorDiv );
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
  
  
  m_resultTxt = new WText( "&nbsp;", resultsDiv );
  m_resultTxt->addStyleClass( "ResultsTxtArea" );
  m_resultTxt->setInline( false );
  
  // Now put the "more info..." link below here and to the right
  m_moreInfoButton = new WPushButton( resultsDiv );
  m_moreInfoButton->setText( WString::tr("dls-further-details-link") );
  m_moreInfoButton->setStyleClass( "LinkBtn MdaMoreInfoBtn" );
  m_moreInfoButton->clicked().connect( this, &DetectionLimitSimple::createMoreInfoWindow );
  m_moreInfoButton->setHiddenKeepsGeometry( true );
  m_moreInfoButton->hide();
 
  
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
  
  WText *dummyThirdRow = new WText( "&nbsp;", generalInput );
  dummyThirdRow->addStyleClass( "GridThirdCol GridFirstRow GridStretchCol SpacerColumn" );
  
  
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
  
  
  WLabel *gammaLabel = new WLabel( WString("{1}:").arg(WString::tr("Gamma")), generalInput );
  gammaLabel->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
  
  
  m_photoPeakEnergy = new WComboBox( generalInput );
  m_photoPeakEnergy->activated().connect( this, &DetectionLimitSimple::handleGammaChanged );
  m_photoPeakEnergy->addStyleClass( "GridSecondCol GridThirdRow GridVertCenter PhotopeakComboBox" );
  
  // TODO: Add FWHM input.  Add text for DRF default, or the button to fit from data.
  //       when user changes this value - dont change ROI limits, just recalc deconv, and redraw either decon or Currie
  
  
  // Add Distance input
  WLabel *distanceLabel = new WLabel( WString("{1}:").arg(WString::tr("Distance")), generalInput );
  distanceLabel->addStyleClass( "GridFourthCol GridFirstRow GridVertCenter" );
  
  m_prevDistance = "100 cm";
  m_distance = new WLineEdit( m_prevDistance, generalInput );
  m_distance->addStyleClass( "GridFifthCol GridFirstRow GridStretchCol" );
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
  WLabel *confidenceLabel = new WLabel( WString::tr("dls-conf-level-label"), generalInput );
  confidenceLabel->addStyleClass( "GridFourthCol GridSecondRow GridVertCenter" );
  m_confidenceLevel = new WComboBox( generalInput );
  m_confidenceLevel->addStyleClass( "GridFifthCol GridSecondRow ClComboBox" );
  
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
  m_detectorDisplay->addStyleClass( "DetectorDisplay GridFourthCol GridThirdRow GridSpanTwoCol GridSpanTwoRows GridVertCenter" );
  m_viewer->detectorChanged().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  m_viewer->detectorModified().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  
  
  
  WLabel *lowerRoiLabel = new WLabel( WString::tr("dls-roi-lower-label"), generalInput );
  lowerRoiLabel->addStyleClass( "GridFirstCol GridFourthRow GridVertCenter" );
  
  m_lowerRoi = new NativeFloatSpinBox( generalInput );
  m_lowerRoi->setSpinnerHidden();
  lowerRoiLabel->setBuddy( m_lowerRoi );
  m_lowerRoi->addStyleClass( "GridSecondCol GridFourthRow" );
  
  WLabel *upperRoiLabel = new WLabel( WString::tr("dls-roi-upper-label"), generalInput );
  upperRoiLabel->addStyleClass( "GridFirstCol GridFifthRow GridVertCenter" );
  
  m_upperRoi = new NativeFloatSpinBox( generalInput );
  m_upperRoi->setSpinnerHidden();
  upperRoiLabel->setBuddy( m_upperRoi );
  m_upperRoi->addStyleClass( "GridSecondCol GridFifthRow" );
  
  m_lowerRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  m_upperRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  
  // Num Side Channel
  m_numSideChannelLabel = new WLabel( WString::tr("dls-num-side-channel-label"), generalInput );
  m_numSideChannelLabel->addStyleClass( "GridFourthCol GridFifthRow GridVertCenter" );
  m_numSideChannel = new WSpinBox( generalInput );
  m_numSideChannel->addStyleClass( "GridFifthCol GridFifthRow" );
  m_numSideChannel->setRange( 1, 64 );
  m_numSideChannel->setValue( 4 );
  m_numSideChannelLabel->setBuddy( m_numSideChannel );
  m_numSideChannel->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedNumSideChannel );
  
  m_numSideChannelLabel->setHiddenKeepsGeometry( true );
  m_numSideChannel->setHiddenKeepsGeometry( true );
  
  WLabel *fwhmLabel = new WLabel( WString::tr("dls-fwhm-label"), generalInput );
  fwhmLabel->addStyleClass( "GridFirstCol GridSixthRow" );
  m_fwhm = new NativeFloatSpinBox( generalInput );
  m_fwhm->setRange( 0.05f, 250.0f );
  m_fwhm->setSpinnerHidden();
  m_fwhm->addStyleClass( "GridSecondCol GridSixthRow" );
  fwhmLabel->setBuddy( m_fwhm );
  m_fwhm->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedFwhm );
  
  m_fwhmSuggestTxt = new WText( generalInput );
  m_fwhmSuggestTxt->addStyleClass( "FwhmSuggest GridThirdCol GridSixthRow GridVertCenter GridSpanTwoCol" );
  
  m_addFwhmBtn = new WPushButton( WString::tr("dls-fit-fwhm-btn"), generalInput );
  m_addFwhmBtn->clicked().connect( this, &DetectionLimitSimple::handleFitFwhmRequested );
  m_addFwhmBtn->addStyleClass( "MdaFitFwhm LightButton GridFifthCol GridSixthRow" );
  
  m_selectDetectorBtn = new WPushButton( WString::tr("dls-select-drf-btn"), generalInput );
  m_selectDetectorBtn->clicked().connect( this, &DetectionLimitSimple::handleSelectDetectorRequested );
  m_selectDetectorBtn->addStyleClass( "MdaFitFwhm LightButton GridFifthCol GridSixthRow" );
  
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  m_addFwhmBtn->setHidden( !drf || !drf->isValid() || drf->hasResolutionInfo() );
  m_selectDetectorBtn->setHidden( drf && drf->isValid() );
  
  m_continuumPriorLabel = new WLabel( WString::tr("dls-deon-cont-norm-label"), generalInput );
  m_continuumPriorLabel->addStyleClass( "GridFirstCol GridSeventhRow GridVertCenter" );
  m_continuumPrior = new WComboBox( generalInput );
  m_continuumPrior->addItem( WString::tr("dls-cont-norm-unknown") );
  m_continuumPrior->addItem( WString::tr("dls-cont-norm-not-present") );
  m_continuumPrior->addItem( WString::tr("dls-cont-norm-fixed-by-sides") );
  m_continuumPrior->setCurrentIndex( 0 );
  m_continuumPrior->activated().connect( this, &DetectionLimitSimple::handleDeconPriorChange );
  m_continuumPrior->addStyleClass( "ContTypeCombo GridSecondCol GridSeventhRow" );
  
  m_continuumPriorLabel->setHiddenKeepsGeometry( true );
  m_continuumPrior->setHiddenKeepsGeometry( true );
  
  m_continuumTypeLabel = new WLabel( "Continuum Type:", generalInput );
  m_continuumTypeLabel->addStyleClass( "GridFourthCol GridSeventhRow GridVertCenter" );
  m_continuumType = new WComboBox( generalInput );
  m_continuumType->addItem( WString::tr( PeakContinuum::offset_type_label_tr(PeakContinuum::OffsetType::Linear) ) );
  m_continuumType->addItem( WString::tr( PeakContinuum::offset_type_label_tr(PeakContinuum::OffsetType::Quadratic) ) );
  m_continuumType->setCurrentIndex( 0 );
  m_continuumType->activated().connect( this, &DetectionLimitSimple::handleDeconContinuumTypeChange );
  m_continuumType->addStyleClass( "GridFifthCol GridSeventhRow" );
  
  m_continuumPriorLabel->hide();
  m_continuumPrior->hide();
  m_continuumTypeLabel->hide();
  m_continuumType->hide();
  
  WContainerWidget *container = new WContainerWidget( generalInput );
  container->addStyleClass( "MethodSelect GridFirstCol GridEighthRow GridSpanFiveCol" );
  
  WLabel *methodLabel = new WLabel( WString::tr("dls-calc-method"), container);
  
  m_methodGroup = new WButtonGroup( container );
  WRadioButton *currieBtn = new Wt::WRadioButton( WString::tr("dls-currie-tab-title"), container );
  m_methodGroup->addButton(currieBtn, static_cast<int>(MethodIds::Currie) );
  
  WRadioButton *deconvBtn = new Wt::WRadioButton( WString::tr("dls-decon-tab-title"), container);
  m_methodGroup->addButton(deconvBtn, static_cast<int>(MethodIds::Deconvolution) );
  m_methodGroup->setCheckedButton( currieBtn );
  
  m_methodGroup->checkedChanged().connect( this, &DetectionLimitSimple::handleMethodChanged );
  
  m_methodDescription = new WText( WString::tr("dls-currie-desc"), generalInput );
  m_methodDescription->addStyleClass( "CalcMethodDesc GridSecondCol GridNinthRow GridSpanFourCol" );
  
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
  {
    // TODO: we could implement updating things as the user drags... not quite sure what is required though in PeakModel - need to check what InterSpec does
    return;
  }
  
  if( new_roi_upper_energy < new_roi_lower_energy )
    std::swap( new_roi_upper_energy, new_roi_lower_energy );
  
  if( m_currentNuclide
     && ((m_currentEnergy < new_roi_lower_energy) || (m_currentEnergy > new_roi_upper_energy)) )
  {
    if( is_final_range )
      passMessage( WString::tr("dls-roi-changed-no-gamma"), WarningWidget::WarningMsgHigh );
    return;
  }
  
  m_lowerRoi->setValue( new_roi_lower_energy );
  m_upperRoi->setValue( new_roi_upper_energy );
  
  handleUserChangedRoi();
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
  
  
  // If there isnt a photopeak selected, update predicted FWHM based on center of ROI
  if( meanEnergy <= 10.0f )
    setFwhmFromEstimate();
  
  
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


void DetectionLimitSimple::setFwhmFromEstimate()
{
  float energy = m_currentEnergy;
  if( energy <= 10 )
    energy = 0.5f*(m_lowerRoi->value() + m_upperRoi->value());
  
  float fwhm = 0.1f;
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  
  m_addFwhmBtn->setHidden( !drf || !drf->isValid() || drf->hasResolutionInfo() );
  m_selectDetectorBtn->setHidden( drf && drf->isValid() );
  
  if( drf && drf->hasResolutionInfo() )
  {
    fwhm = drf->peakResolutionFWHM( energy );
  }else
  {
    fwhm = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(energy) );
  }
  
  m_fwhm->setValue( fwhm );
  m_fwhmSuggestTxt->hide();
}//setFwhmFromEstimate();


void DetectionLimitSimple::handleUserChangedFwhm()
{
  float fwhm = m_fwhm->value();
  float energy = m_currentEnergy;
  if( energy <= 10 )
    energy = 0.5f*(m_lowerRoi->value() + m_upperRoi->value());
  
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  
  if( (m_fwhm->validate() != WValidator::State::Valid) || (m_fwhm->value() < 0.1f) )
  {
    // I'm not actually sure if we can ever make it here
    if( drf && drf->hasResolutionInfo() )
      fwhm = drf->peakResolutionFWHM( energy );
    else
      fwhm = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(energy) );
    m_fwhm->setValue( fwhm );
  }//if( invalid FWHM )
  
  m_addFwhmBtn->setHidden( !drf || !drf->isValid() || drf->hasResolutionInfo() );
  m_selectDetectorBtn->setHidden( drf && drf->isValid() );
  
  if( drf && drf->hasResolutionInfo() )
  {
    const double drf_fwhm = drf->peakResolutionFWHM( energy );
    if( fabs(fwhm - drf_fwhm) > 0.1 )
    {
      char text[32] = { '\0' };
      snprintf( text, sizeof(text), "%.2f", drf_fwhm );
      m_fwhmSuggestTxt->setText( WString::tr("dls-suggest-fwhm").arg(text) );
      m_fwhmSuggestTxt->show();
    }else
    {
      m_fwhmSuggestTxt->hide();
    }
  }else
  {
    const float est_fwhm = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(energy) );
    
    char text[32] = { '\0' };
    snprintf( text, sizeof(text), "%.2f", est_fwhm );
    m_fwhmSuggestTxt->setText( WString::tr("dls-rough-est-fwhm").arg(est_fwhm) ); //"No functional FWHM"
    m_fwhmSuggestTxt->show();
  }//if( DRF has FEHM info ) / else
  
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleUserChangedFwhm()


void DetectionLimitSimple::handleDeconPriorChange()
{
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  assert( !currieMethod );
  
  const bool useSideChan = (m_continuumPrior->currentIndex() == 2);
  
  m_numSideChannelLabel->setHidden( !currieMethod && !useSideChan );
  m_numSideChannel->setHidden( !currieMethod && !useSideChan );
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleDeconPriorChange()


void DetectionLimitSimple::handleDeconContinuumTypeChange()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleDeconContinuumTypeChange()

DetectionLimitSimple::~DetectionLimitSimple()
{
  //nothing to do here
}//~DoseCalcWidget()


void DetectionLimitSimple::handleMethodChanged()
{
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  
  m_numSideChannelLabel->setHidden( !currieMethod );
  m_numSideChannel->setHidden( !currieMethod );
  
  const bool useSideChan = (m_continuumPrior->currentIndex() == 2);
  m_numSideChannelLabel->setHidden( !currieMethod && !useSideChan );
  m_numSideChannel->setHidden( !currieMethod && !useSideChan );
  
  m_continuumPriorLabel->setHidden( currieMethod );
  m_continuumPrior->setHidden( currieMethod );
  m_continuumTypeLabel->setHidden( currieMethod );
  m_continuumType->setHidden( currieMethod );
  
  m_methodDescription->setText( WString::tr(currieMethod ? "dls-currie-desc" : "dls-decon-desc") );
  
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleMethodChanged()


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
    updateSpectrumDecorationsAndResultText();
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
  m_photoPeakEnergy->setDisabled( true );
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
    
    if( intensity > std::numeric_limits<double>::epsilon() )
    {
      m_photoPeakEnergiesAndBr.push_back( make_pair(energy, intensity) );
      m_photoPeakEnergy->addItem( text );
    }//if( intensity > 0.0 )
  }//for each( const double energy, energies )
  
  
  // we wont change `m_currentEnergy`, since the user might go back to their previous nuclide
  m_photoPeakEnergy->setDisabled( (m_photoPeakEnergy->count() == 0) );
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
  
  setFwhmFromEstimate();
  
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
  
  m_addFwhmBtn->setHidden( !new_drf || !new_drf->isValid() || new_drf->hasResolutionInfo() );
  m_selectDetectorBtn->setHidden( new_drf && new_drf->isValid() );
  
  handleUserChangedFwhm();
  
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


void DetectionLimitSimple::handleSelectDetectorRequested()
{
  m_detectorDisplay->editDetector();
}//void handleSelectDetectorRequested()


void DetectionLimitSimple::handleSpectrumChanged()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleSpectrumChanged()



void DetectionLimitSimple::updateSpectrumDecorationsAndResultText()
{
  shared_ptr<const ColorTheme> theme = m_viewer->getColorTheme();
  assert( theme );
  
  const bool use_curie = use_curie_units();
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  
  m_resultTxt->setText( "---" );
  m_moreInfoButton->hide();
  m_peakModel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
  m_spectrum->removeAllDecorativeHighlightRegions();
  
  const double confidence_level = m_currentDeconResults ? m_currentDeconResults->confidenceLevel : currentConfidenceLevel();
  
  const string cl_str = ([confidence_level]() -> string {
    char cl_str_buffer[64] = {'\0'};
    if( confidence_level < 0.999 )
      snprintf( cl_str_buffer, sizeof(cl_str_buffer), "%.1f%%", 100.0*confidence_level );
    else
      snprintf( cl_str_buffer, sizeof(cl_str_buffer), "1-%.2G", (1.0-confidence_level) );
    return cl_str_buffer;
  })();
  
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  if( currieMethod )
  {
    // Currie method limit
    if( m_currentCurrieInput )
    {
      double gammas_per_bq = -1.0, distance = -1.0, br = -1;
      vector<DetectionLimitTool::CurrieResultPeak> roi_peaks;
      
      try
      {
        distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
      }catch( std::exception & )
      {
      }
      
      const int energyIndex = m_photoPeakEnergy->currentIndex();
      if( (energyIndex >= 0) && (energyIndex < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
      {
        // TODO: consider being able to draw each peak individually.
        const double energy = m_photoPeakEnergiesAndBr[energyIndex].first;
        assert( fabs(energy - static_cast<float>(m_currentCurrieInput->gamma_energy)) < 0.1 );
        
        if( m_allGammasInRoi )
        {
          //include any gamma within ROI
          br = 0.0;
          
          const double gammaFwhm = m_fwhm->value();
          const bool isDrfFwhm = (drf && drf->hasResolutionInfo()
                                  && (fabs(gammaFwhm - drf->peakResolutionFWHM(energy)) < 0.01));
          
          const double roi_lower = m_lowerRoi->value();
          const double roi_upper = m_upperRoi->value();
          
          for( const pair<double,double> &ppebr : m_photoPeakEnergiesAndBr )
          {
            if( (ppebr.first >= roi_lower) && (ppebr.first <= roi_upper) )
            {
              br += ppebr.second;
              
              DetectionLimitTool::CurrieResultPeak p;
              p.energy = ppebr.first;
              p.fwhm = isDrfFwhm ? drf->peakResolutionFWHM(ppebr.first) : gammaFwhm;
              p.counts_4pi = ppebr.second;
              
              roi_peaks.push_back( std::move(p) );
            }
          }//for( const pair<double,double> &ppebr : m_photoPeakEnergiesAndBr )
          
          assert( br >= m_photoPeakEnergiesAndBr[energyIndex].second );
          br = std::max( br, m_photoPeakEnergiesAndBr[energyIndex].second ); //JIC
        }else
        {
          // Only use the selected gamma
          br = m_photoPeakEnergiesAndBr[energyIndex].second;
          
          DetectionLimitTool::CurrieResultPeak p;
          p.energy = energy;
          p.fwhm = m_fwhm->value();
          p.counts_4pi = br;
          roi_peaks.push_back( std::move(p) );
        }
      }//if( an energy is selected )
        
      
      if( drf && drf->isValid() && (distance >= 0.0) && (br > 0) && m_currentCurrieInput->spectrum )
      {
        const bool fixed_geom = drf->isFixedGeometry();
        
        boost::function<double(float)> att_coef_fcn, air_atten_fcn;
        if( distance > 0.0 )
          air_atten_fcn = boost::bind( &GammaInteractionCalc::transmission_coefficient_air, _1, distance );
        
        {//begin convert `br` to `gammas_per_bq`
          const float energy = m_currentCurrieInput->gamma_energy;
          const double det_eff = fixed_geom ? drf->intrinsicEfficiency(energy)
          : drf->efficiency(energy, distance);
          
          const double shield_transmission = att_coef_fcn.empty() ? 1.0 : exp( -1.0*att_coef_fcn(energy) );
          const double air_transmission = air_atten_fcn.empty() ? 1.0 : exp( -1.0*air_atten_fcn(energy) );
          const double counts_per_bq_into_4pi = br * shield_transmission * m_currentCurrieInput->spectrum->live_time();
          const double counts_per_bq_into_4pi_with_air = air_transmission * counts_per_bq_into_4pi;
          const double counts_4pi = fixed_geom ? counts_per_bq_into_4pi : counts_per_bq_into_4pi_with_air;
          
          gammas_per_bq = counts_4pi * det_eff;
        }//end convert `br` to `gammas_per_bq`
        
        //Now go through and correct get<2>(roi_peaks[i]) and then modify update_spectrum_for_currie_result
        for( DetectionLimitTool::CurrieResultPeak &peak : roi_peaks )
        {
          const double &peak_energy = peak.energy;
          const double peak_br = peak.counts_4pi;
          const double shield_transmission = att_coef_fcn.empty() ? 1.0 : exp( -1.0*att_coef_fcn(peak_energy) );
          const double air_transmission = air_atten_fcn.empty() ? 1.0 : exp( -1.0*air_atten_fcn(peak_energy) );
          const double counts_per_bq_into_4pi = peak_br * shield_transmission * m_currentCurrieInput->spectrum->live_time();
          const double counts_per_bq_into_4pi_with_air = air_transmission * counts_per_bq_into_4pi;
          const double counts_4pi = fixed_geom ? counts_per_bq_into_4pi : counts_per_bq_into_4pi_with_air;
          
          const double det_eff = fixed_geom ? drf->intrinsicEfficiency(peak_energy)
                                            : drf->efficiency(peak_energy, distance);
          
          peak.counts_4pi = counts_4pi * det_eff;
        }//for( DetectionLimitTool::CurrieResultPeak &peak : roi_peaks )
      }//if( drf )
      
      const DetectionLimitTool::LimitType limitType = DetectionLimitTool::LimitType::Activity;
      
      const DetectionLimitCalc::CurrieMdaResult * const result = m_currentCurrieResults.get();
      
      DetectionLimitTool::update_spectrum_for_currie_result( m_spectrum, m_peakModel,
              *m_currentCurrieInput, result, drf, limitType, gammas_per_bq, roi_peaks );
      
      
      WString result_txt;
      const bool use_curie = use_curie_units();
      const DetectorPeakResponse::EffGeometryType det_geom = drf ? drf->geometryType() : DetectorPeakResponse::EffGeometryType::FarField;
      
      if( result->source_counts > result->decision_threshold )
      {
        // There is enough excess counts that we would reliably detect this activity, so we will
        //  give the activity range.
        string lowerstr, upperstr, nomstr;
        
        if( gammas_per_bq > 0.0 )
        {
          const float lower_act = result->lower_limit / gammas_per_bq;
          const float upper_act = result->upper_limit / gammas_per_bq;
          const float nominal_act = result->source_counts / gammas_per_bq;
          
          lowerstr = PhysicalUnits::printToBestActivityUnits( lower_act, 2, use_curie )
          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
          upperstr = PhysicalUnits::printToBestActivityUnits( upper_act, 2, use_curie )
          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
          nomstr = PhysicalUnits::printToBestActivityUnits( nominal_act, 2, use_curie )
          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
          
          result_txt = WString::tr("dls-det-act-with-range").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }else
        {
          lowerstr = SpecUtils::printCompact(result->lower_limit, 4);
          upperstr = SpecUtils::printCompact(result->upper_limit, 4);
          nomstr = SpecUtils::printCompact(result->source_counts, 4);
          
          result_txt = WString::tr("dls-det-counts-with-range").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }
      }else if( result->upper_limit < 0 )
      {
        // This can happen when there are a lot fewer counts in the peak region than predicted
        //  from the sides - since this is non-sensical, we'll just say zero.
        const string unitstr = use_curie ? "Ci" : "Bq";
        
        if( gammas_per_bq > 0.0 )
        {
          result_txt = WString::tr("dls-det-act-less-zero").arg(unitstr);
        }else
        {
          result_txt = WString::tr("dls-det-counts-less-zero").arg(unitstr);
        }
      }else
      {
        // We will provide the upper bound on activity.
        string mdastr;
        if( gammas_per_bq > 0.0 )
        {
          const double simple_mda = result->upper_limit / gammas_per_bq;
          mdastr = PhysicalUnits::printToBestActivityUnits( simple_mda, 2, use_curie )
                  + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
        }else
        {
          mdastr = SpecUtils::printCompact( result->upper_limit, 4 ) + " signal counts";
        }
        
        result_txt = WString::tr("dls-det-upper-bound").arg(mdastr).arg(cl_str);
      }//if( detected ) / else if( ....)
      
      WString full_result_txt( "{1}<br/>{2}" );
      full_result_txt.arg(result_txt);
      
      if( gammas_per_bq > 0.0 )
      {
        const double detection_act = result->detection_limit / gammas_per_bq;
        const string act = PhysicalUnits::printToBestActivityUnits( detection_act, 2, use_curie )
                          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );
        
        full_result_txt.arg( WString::tr("dls-min-detectable-act").arg(act) );
      }else
      {
        const string counts = SpecUtils::printCompact(result->detection_limit, 4);
        full_result_txt.arg( WString::tr("dls-min-detectable-counts").arg( counts ) );
      }
      
      m_resultTxt->setText( full_result_txt );
      m_moreInfoButton->show();
    }//if( m_currentCurrieInput )
  }else
  {
    // Deconvolution method limit
    std::vector<PeakDef> fit_peaks;
    
    if( !m_currentDeconResults )
    {
      m_spectrum->setChartTitle( WString::tr("dls-error-computing-results") );
      m_peakModel->setPeaks( vector<PeakDef>{} );
    }else
    {
      assert( !!m_currentDeconInput );
      assert( m_currentDeconResults->isDistanceLimit == false );
      
      const DetectionLimitCalc::DeconActivityOrDistanceLimitResult &result = *m_currentDeconResults;
      
      //assert( result.foundLowerCl || result.foundUpperCl );
      
      if( !m_currentNuclide )
      {
        // Sanity check that we filled out the limit input how we expect it.
        assert( result.baseInput.roi_info.size() == 1 );
        assert( result.baseInput.roi_info[0].peak_infos.size() == 1 );
        assert( result.baseInput.roi_info[0].peak_infos[0].counts_per_bq_into_4pi == result.baseInput.measurement->live_time() );
      }//if( !m_currentNuclide )
      
      WString chart_title, result_txt;
      double display_activity = 0.0;
      
      if( result.foundLowerCl && result.foundUpperCl )
      {
        display_activity = result.overallBestQuantity;
        
        if( !m_currentNuclide )
        {
          //result.lowerLimit
          assert( result.lowerLimitResults
                 && (result.lowerLimitResults->fit_peaks.size() == 1) );
          assert( result.overallBestResults
                 && (result.overallBestResults->fit_peaks.size() == 1) );
          assert( result.upperLimitResults
                 && (result.upperLimitResults->fit_peaks.size() == 1) );
          
          const double lower_limit_counts = (result.lowerLimitResults
                                  && (result.lowerLimitResults->fit_peaks.size() == 1))
                                        ? result.lowerLimitResults->fit_peaks[0].amplitude() : -1.0;
          const double nominal_counts = (result.overallBestResults
                                  && (result.overallBestResults->fit_peaks.size() == 1))
                                       ? result.overallBestResults->fit_peaks[0].amplitude() : -1.0;
          const double upper_limit_counts = (result.upperLimitResults
                                  && (result.upperLimitResults->fit_peaks.size() == 1))
                                        ? result.upperLimitResults->fit_peaks[0].amplitude() : -1.0;
          
          const string lowerstr = SpecUtils::printCompact(lower_limit_counts, 4);
          const string nomstr = SpecUtils::printCompact(nominal_counts, 4);
          const string upperstr = SpecUtils::printCompact(upper_limit_counts, 4);
          
          chart_title = WString::tr("dls-chart-title-estimated-counts").arg( nomstr );
          result_txt = WString::tr("dls-results-txt-estimated-counts").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }else
        {
          const string nomstr = PhysicalUnits::printToBestActivityUnits(result.overallBestQuantity, 3, use_curie);
          const string lowerstr = PhysicalUnits::printToBestActivityUnits(result.lowerLimit, 3, use_curie);
          const string upperstr = PhysicalUnits::printToBestActivityUnits(result.upperLimit, 3, use_curie);
          
          const string cl_txt = "Estimated activity of " + nomstr + ".";
          
          const string sum_txt = "Detected activity of " + nomstr + "."
          "<br/>"
          "Range: [" + lowerstr + ", " + upperstr + "] @"  + cl_str + " CL";
          
          chart_title = WString::tr("dls-chart-title-estimated-act").arg(nomstr);
          result_txt = WString::tr("dls-results-txt-estimated-act").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }//if( !m_currentNuclide ) / else
        
        assert( result.overallBestResults );
        if( result.overallBestResults )
          fit_peaks = result.overallBestResults->fit_peaks;
      }else if( result.foundLowerCl )
      {
        assert( 0 );
        display_activity = 0.0; //result.foundLowerCl
        const string cl_txt = "Error: Didn't find " + cl_str + " CL activity";
        const string sum_txt = "Error: Didn't find " + cl_str + " CL activity";
        
        chart_title = WString::fromUTF8( cl_txt );
        result_txt = WString::fromUTF8( sum_txt );
        
        //fit_peaks = result.lowerLimitResults.fit_peaks;
      }else if( result.foundUpperCl )
      {
        display_activity = result.foundUpperCl;
        
        if( !m_currentNuclide )
        {
          const double upper_limit_counts = (result.upperLimitResults
                                  && (result.upperLimitResults->fit_peaks.size() == 1))
                                        ? result.upperLimitResults->fit_peaks[0].amplitude() : -1.0;
          const string upperstr = SpecUtils::printCompact(upper_limit_counts, 4);
          
          chart_title = WString::tr("dls-chart-title-upper-bound-counts").arg(upperstr).arg(cl_str);
          result_txt = WString::tr("dls-results-txt-upper-bound-counts").arg(upperstr).arg(cl_str);
        }else
        {
          const string upperstr = PhysicalUnits::printToBestActivityUnits(result.upperLimit, 3, use_curie);
          chart_title = WString::tr("dls-chart-title-upper-bound-act").arg(upperstr).arg(cl_str);
          result_txt = WString::tr("dls-results-txt-upper-bound-act").arg(upperstr).arg(cl_str);
        }//if( !m_currentNuclide ) / else
        
        assert( result.upperLimitResults );
        if( result.upperLimitResults )
          fit_peaks = result.upperLimitResults->fit_peaks;
      }else
      {
        display_activity = 0.0;
        const string cl_txt = "Error: failed upper or lower limits at " + cl_str;
        chart_title = WString::fromUTF8( cl_txt );
        
        const string sum_txt = "Error: failed upper or lower limits at " + cl_str;
        result_txt = WString::fromUTF8( sum_txt );
      }
      
      const float lower_energy = m_lowerRoi->value();
      const float upper_energy = m_upperRoi->value();
      const double dx = upper_energy - lower_energy;
      m_spectrum->setXAxisRange( lower_energy - 0.5*dx, upper_energy + 0.5*dx );
      
      m_spectrum->setChartTitle( chart_title );
      m_peakModel->setPeaks( fit_peaks );
      m_resultTxt->setText( result_txt );
      m_moreInfoButton->show();
      
      //m_currentDeconResults->foundUpperDisplay = false;
      //m_currentDeconResults->upperDisplayRange = 0.0;
      //m_currentDeconResults->foundLowerDisplay = false;
      //m_currentDeconResults->lowerDisplayRange = 0.0;
      //std::vector<std::pair<double,double>> m_currentDeconResults->chi2s;
    }//if( no valid result
  }//if( currently doing Currie-style limit ) / else
}//void updateSpectrumDecorationsAndResultText()


double DetectionLimitSimple::currentConfidenceLevel() const
{
  double confidenceLevel = 0.95;
  
  const int clIndex = m_confidenceLevel->currentIndex();
  const ConfidenceLevel confidence = ConfidenceLevel(clIndex);
  
  switch( confidence )
  {
    case OneSigma:   confidenceLevel = 0.682689492137086; break;
    case TwoSigma:   confidenceLevel = 0.954499736103642; break;
    case ThreeSigma: confidenceLevel = 0.997300203936740; break;
    case FourSigma:  confidenceLevel = 0.999936657516334; break;
    case FiveSigma:  confidenceLevel = 0.999999426696856; break;
    case NumConfidenceLevel: assert(0); break;
  }//switch( confidence )
  
  return confidenceLevel;
}//double currentConfidenceLevel() const


SimpleDialog *DetectionLimitSimple::createDeconvolutionLimitMoreInfo()
{
  if( !m_currentDeconInput || !m_currentDeconResults )
    throw runtime_error( "No solution available." );
  
  const bool use_curie = use_curie_units();
  const DetectionLimitCalc::DeconComputeInput &input = *m_currentDeconInput;
  const DetectionLimitCalc::DeconActivityOrDistanceLimitResult &result = *m_currentDeconResults;
  
  assert( result.baseInput.roi_info.size() == 1 );
  if( !result.baseInput.roi_info.size() )
    throw runtime_error( "No ROI info available" );
  
  assert( !result.isDistanceLimit );
  
  double distance = 0.0;
  try
  {
    distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
  }catch( std::exception & )
  {
    distance = -1.0;
  }
  
  shared_ptr<const SpecUtils::Measurement> measurement = result.baseInput.measurement;
  assert( measurement );
  if( !measurement )
    throw runtime_error( WString::tr("dls-err-no-meas").toUTF8() );
  
  double energy = 0.0;
  const int energyIndex = m_photoPeakEnergy->currentIndex();
  if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    energy = 0.5*(m_lowerRoi->value() + m_upperRoi->value());
  else
    energy = m_photoPeakEnergiesAndBr[energyIndex].first;
  
  const float roi_start = result.baseInput.roi_info[0].roi_start;
  const float roi_end = result.baseInput.roi_info[0].roi_end;
  
  wApp->require( "InterSpec_resources/DetectionLimitTool.js" );
  
  char buffer[256] = { '\0' };
  snprintf( buffer, sizeof(buffer), "%s%.2f keV {1}",
           (m_currentNuclide ? (m_currentNuclide->symbol + " ").c_str() : ""), energy );
  
  SimpleDialog *dialog = new SimpleDialog( WString(buffer).arg(WString::tr("dls-Info")) );
  dialog->addButton( WString::tr("Close") );
  
  WContainerWidget *contents = new WContainerWidget( dialog->contents() );
  contents->addStyleClass( "DeconvMoreInfo" );
  
  // Create a chi2 chart
  WContainerWidget *chi2Chart = new WContainerWidget( contents );
  chi2Chart->addStyleClass( "DeconChi2Chart" );
  
  chi2Chart->setJavaScriptMember( "chart", "new MdaChi2Chart(" + chi2Chart->jsRef() + ", {});");
  const string jsgraph = chi2Chart->jsRef() + ".chart";
  
  chi2Chart->setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + chi2Chart->id() + "') && "
             + chi2Chart->jsRef() + " && " + jsgraph + " )"
          + jsgraph + ".redraw();"
        "}"
      "});"
  );
  chi2Chart->callJavaScriptMember( "resizeObserver.observe", chi2Chart->jsRef() );
  
  const Wt::Json::Object chartJson = DetectionLimitTool::generateChartJson( result, false );
  const string datajson = Wt::Json::serialize(chartJson);
  chi2Chart->doJavaScript( jsgraph + ".setData(" + datajson + ");" );
  
  if( !m_currentNuclide )
  {
    WText *txt = new WText( WString::tr("dls-assumed-br=1"), contents );
    txt->addStyleClass( "AssumedBrNote" );
    txt->setInline( false );
  }//if( !m_currentNuclide )
  
  // Now create rows of text information.
  WTable *table = new WTable( contents );
  table->addStyleClass( "DeconvoMoreInfoTable" );
  
  const auto print_result = [table, use_curie, measurement, roi_start, roi_end]( 
                                        const DetectionLimitCalc::DeconComputeResults &result,
                                        const bool is_best, const WString typestr ){
    WString label = WString("{1} {2}").arg(typestr).arg( WString::tr( (is_best ? "Activity" : "dls-Limit") ) );
    WString value = PhysicalUnits::printToBestActivityUnits( result.input.activity, 3, use_curie );
    
    WTableCell *cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    
    
    label = WString("{1} {2}").arg(typestr).arg( WString::tr( (is_best ? "Counts" : "dls-Limit") ) );
    double counts = 0.0, uncert = 0.0;
    for( const auto peak : result.fit_peaks )
    {
      counts += peak.peakArea();
      uncert += peak.peakAreaUncert() * peak.peakAreaUncert(); //We dont actually hav an uncertainty
    }
    value = SpecUtils::printCompact(counts, 5) + " counts";
    //value = PhysicalUnits::printValueWithUncertainty(counts, sqrt(uncert), 5);
    
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    
    label = WString("{1} &chi;<sup>2</sup>").arg(typestr);
    value = SpecUtils::printCompact(result.chi2, 4);
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    
    if( !is_best )
      return;
    
    label = WString::tr("dls-DOF");
    value = std::to_string( result.num_degree_of_freedom );
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    
    if( result.fit_peaks.size() )
    {
      label = WString::tr("dls-continuum-area");
      const PeakDef &peak = result.fit_peaks.front();
      const double cont_area = peak.continuum()->offset_integral(roi_start, roi_end, measurement);
      value = SpecUtils::printCompact(cont_area, 5);
      
      cell = table->elementAt( table->rowCount(), 0 );
      new WText( label, cell );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      new WText( value, cell );
    }//if( result.overallBestResults->fit_peaks.size() )
  };//const auto print_result
  
  if( result.foundUpperCl && result.upperLimitResults )
    print_result( *result.upperLimitResults, false, WString::tr("dls-Upper") );
  
  if( result.foundLowerCl && result.lowerLimitResults )
    print_result( *result.lowerLimitResults, false, WString::tr("dls-Lower") );
  
  if( result.overallBestResults )
    print_result( *result.overallBestResults, true, WString::tr("dls-Best") );
  
  
  WString label, value;
  
  label = WString::tr("FWHM");
  const double fwhm = m_fwhm->value();
  snprintf( buffer, sizeof(buffer), "%.2f keV", fwhm );
  value = WString::fromUTF8( buffer );
  WTableCell *cell = table->elementAt( table->rowCount(), 0 );
  new WText( label, cell );
  cell = table->elementAt( table->rowCount() - 1, 1 );
  new WText( value, cell );
  
  label = WString::tr("dls-ROI-range-label");
  snprintf( buffer, sizeof(buffer), "[%.2f, %.2f]", roi_start, roi_end );
  value = WString::fromUTF8( buffer );
  cell = table->elementAt( table->rowCount(), 0 );
  new WText( label, cell );
  cell = table->elementAt( table->rowCount() - 1, 1 );
  new WText( value, cell );
  
  label = WString::tr("dls-ROI-Channels-label");
  const size_t lower_chan = measurement->find_gamma_channel( roi_start + 0.0001 );
  const size_t upper_chan = measurement->find_gamma_channel( roi_end - 0.0001 );
  const string channels_str = "[" + std::to_string(lower_chan) + ", " + std::to_string(upper_chan) + "]";
  value = WString::fromUTF8( channels_str );
  cell = table->elementAt( table->rowCount(), 0 );
  new WText( label, cell );
  cell = table->elementAt( table->rowCount() - 1, 1 );
  new WText( value, cell );
  
  label = WString::tr("dls-ROI-Width");
  snprintf( buffer, sizeof(buffer), "%.3f %s", (roi_end - roi_start)/fwhm, WString::tr("FWHM").toUTF8().c_str() );
  value = WString::fromUTF8( buffer );
  cell = table->elementAt( table->rowCount(), 0 );
  new WText( label, cell );
  cell = table->elementAt( table->rowCount() - 1, 1 );
  new WText( value, cell );
  
  shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  
  // Add a blank row
  cell = table->elementAt( table->rowCount(), 0 );
  new WText( "&nbsp;", TextFormat::XHTMLText, cell );
  
  if( drf && drf->isValid() )
  {
    const double intrinsic_eff = drf->intrinsicEfficiency( energy );
    
    label = WString::tr("dls-det-intrinsic-eff");
    value = SpecUtils::printCompact( intrinsic_eff, 5 );
    
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    //addTooltipToRow( "The efficiency for a gamma hitting the detector face,"
    //                " to be detected in the full-energy peak." );
    
    if( distance >= 0.0 )
    {
      const double geom_eff = drf->fractionalSolidAngle( drf->detectorDiameter(), distance );
      
      label = WString::tr("dls-solid-angle-frac");
      value = SpecUtils::printCompact( geom_eff, 5 );
      
      cell = table->elementAt( table->rowCount(), 0 );
      new WText( label, cell );
      cell = table->elementAt( table->rowCount() - 1, 1 );
      new WText( value, cell );
      //addTooltipToRow( "The fraction of the solid angle, the detector face takes up, at the specified distance." );
    }//if( distance >= 0.0 )
  }//if( drf )
  
  
  
  if( (distance > 0.0)
     && (!drf || (drf->geometryType() == DetectorPeakResponse::EffGeometryType::FarField) ) 
     && result.baseInput.include_air_attenuation )
  {
    const double air_atten_coef = GammaInteractionCalc::transmission_coefficient_air( energy, distance );
    const double air_transmission = exp( -1.0 * air_atten_coef );
    
    label = WString::tr("dls-air-trans");
    value = SpecUtils::printCompact( air_transmission, 5 );
    
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    //addTooltipToRow( "The fraction of gammas, at this energy, that will make it through the air (assuming sea level) without interacting." );
  }//if( air_atten )
  
  double branch_ratio = 0.0;
  const float live_time = measurement->live_time();
  
  for( const auto &roi : result.baseInput.roi_info )
  {
    for( const DetectionLimitCalc::DeconRoiInfo::PeakInfo &peak : roi.peak_infos )
      branch_ratio += peak.counts_per_bq_into_4pi / live_time;
  }//for( const auto &roi : result.baseInput.roi_info )
  
  if( branch_ratio > 0.0 )
  {
    label = WString::tr("dls-gamma-intensity");
    value = SpecUtils::printCompact( branch_ratio, 5 );
    
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    //addTooltipToRow( "The number of gamma rays emitted at this energy, from the radioactive"
    //                " source before any shielding, but accounting for nuclide age,"
    //                " per decay of the parent nuclide." );
  }//if( branch_ratio > 0.0 )
  
  return dialog;
}//void createDeconvolutionLimitMoreInfo()


void DetectionLimitSimple::createMoreInfoWindow()
{
  assert( !m_moreInfoWindow );
  m_moreInfoWindow = nullptr; // we really shouldnt need to do this
  
  double distance = 0.0;
  try
  {
    distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
  }catch( std::exception & )
  {
    distance = -1.0;
  }
  
  double energy = 0.0;
  const int energyIndex = m_photoPeakEnergy->currentIndex();
  if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    energy = 0.5*(m_lowerRoi->value() + m_upperRoi->value());
  else
    energy = m_photoPeakEnergiesAndBr[energyIndex].first;
  
  shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  if( drf && !drf->isValid() )
    drf.reset();
  
  try
  {
    const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
    
    if( currieMethod )
    {
      shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
      if( !hist || hist->num_gamma_channels() < 7 )
        throw runtime_error( "No displayed foreground" );
      
      if( !m_currentCurrieInput || !m_currentCurrieResults )
        throw runtime_error( "No current results" );
      
      double branch_ratio = 0.0;
      
      if( m_currentNuclide )
      {
        const double roi_lower = m_lowerRoi->value();
        const double roi_upper = m_upperRoi->value();
        
        for( const pair<double,double> &ppebr : m_photoPeakEnergiesAndBr )
        {
          if( (ppebr.first >= roi_lower) && (ppebr.first <= roi_upper) )
            branch_ratio += ppebr.second;
        }//for( const pair<double,double> &ppebr : m_photoPeakEnergiesAndBr )
        
        assert( branch_ratio >= m_photoPeakEnergiesAndBr[energyIndex].second );
        branch_ratio = std::max( branch_ratio, m_photoPeakEnergiesAndBr[energyIndex].second ); //JIC
      }//if( m_currentNuclide )
      
      const double shield_transmission = 1.0;
      
      const bool do_air_atten = (distance > 0.0);
      m_moreInfoWindow = DetectionLimitTool::createCurrieRoiMoreInfoWindow( m_currentNuclide,
                    *m_currentCurrieResults, drf, DetectionLimitTool::LimitType::Activity,
                    distance, do_air_atten, branch_ratio, shield_transmission );
    }else
    {
      m_moreInfoWindow = createDeconvolutionLimitMoreInfo();
    }//if( currieMethod ) / else
  }catch( std::exception &e )
  {
    assert( !m_moreInfoWindow );
    m_moreInfoWindow = new SimpleDialog( WString::tr("dls-err-more-info-title"),
                                            WString::tr("dls-err-more-info-content").arg(e.what()) );
    m_moreInfoWindow->addButton( WString::tr("Close") );
  }//try / catch
  
  assert( m_moreInfoWindow );
  if( m_moreInfoWindow )
    m_moreInfoWindow->finished().connect( boost::bind(&DetectionLimitSimple::handleMoreInfoWindowClose, this, m_moreInfoWindow) );
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo_redo = []( const bool is_show ){
      DetectionLimitSimpleWindow *mdawin = InterSpec::instance()->showSimpleMdaWindow();
      DetectionLimitSimple *tool = mdawin ? mdawin->tool() : nullptr;
      assert( tool );
      if( tool && is_show )
        tool->createMoreInfoWindow();
      else if( tool )
        tool->programmaticallyCloseMoreInfoWindow();
    };//undo_redo
      
    auto undo = [undo_redo](){ undo_redo(false); };
    auto redo = [undo_redo](){ undo_redo(true); };
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Show Simple MDA more info window." );
  }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
}//void createMoreInfoWindow()


void DetectionLimitSimple::handleMoreInfoWindowClose( SimpleDialog *dialog )
{
  assert( dialog == dynamic_cast<SimpleDialog *>( WObject::sender() ) );
  assert( dialog == m_moreInfoWindow );
  SimpleDialog *current = m_moreInfoWindow;
  m_moreInfoWindow = nullptr;
  
  if( current && (current == dialog) )
  {
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() )
    {
      auto undo_redo = []( const bool is_show ){
        DetectionLimitSimpleWindow *mdawin = InterSpec::instance()->showSimpleMdaWindow();
        DetectionLimitSimple *tool = mdawin ? mdawin->tool() : nullptr;
        assert( tool );
        if( tool && is_show )
          tool->createMoreInfoWindow();
        else if( tool )
          tool->programmaticallyCloseMoreInfoWindow();
      };//undo_redo
        
      auto undo = [undo_redo](){ undo_redo(true); };
      auto redo = [undo_redo](){ undo_redo(false); };
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Close Simple MDA more info window." );
    }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
  }//if( m_moreInfoWindow && (m_moreInfoWindow == dialog) )
}//void handleMoreInfoWindowClose( SimpleDialog *dialog );


void DetectionLimitSimple::programmaticallyCloseMoreInfoWindow()
{
  SimpleDialog *dialog = m_moreInfoWindow;
  m_moreInfoWindow = nullptr;
  
  if( dialog )
  {
    // Note: dialog wont emit the finished() signal
    if( dialog->isModal() )
      dialog->setModal(false);
    delete dialog;
  }//if( dialog )
}//void programmaticallyCloseMoreInfoWindow()


void DetectionLimitSimple::updateResult()
{
  //m_errMsg->setText( WString::tr("dls-err-no-input") );
  m_errMsg->setText( "" );
  
  m_currentDeconInput.reset();
  m_currentDeconResults.reset();
  
  m_currentCurrieInput.reset();
  m_currentCurrieResults.reset();
  
  try
  {
    m_fitFwhmBtn->hide();
    
    std::shared_ptr<const SpecUtils::Measurement> hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !hist || (hist->num_gamma_channels() < 7) )
      throw runtime_error( "No foreground spectrum loaded." );
    
    const float roi_lower_energy = m_lowerRoi->value();
    const float roi_upper_energy = m_upperRoi->value();
    
    double energy = 0.0;
    const int energyIndex = m_photoPeakEnergy->currentIndex();
    if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
      energy = 0.5*(roi_lower_energy + roi_upper_energy);
    else
      energy = m_photoPeakEnergiesAndBr[energyIndex].first;
    
    const int clIndex = m_confidenceLevel->currentIndex();
    if( (clIndex < 0) || (clIndex >= ConfidenceLevel::NumConfidenceLevel) )
      throw runtime_error( "Please select confidence level." );
      
    const double confidenceLevel = currentConfidenceLevel();
    
    // We need to calculate currie-style limit, even if we want the deconvolution-style limit
    auto currie_input = make_shared<DetectionLimitCalc::CurrieMdaInput>();
    currie_input->spectrum = hist;
    currie_input->gamma_energy = static_cast<float>( energy );
    currie_input->roi_lower_energy = m_lowerRoi->value();
    currie_input->roi_upper_energy = m_upperRoi->value();
    currie_input->num_lower_side_channels = static_cast<size_t>( m_numSideChannel->value() );
    currie_input->num_upper_side_channels = currie_input->num_lower_side_channels;
    currie_input->detection_probability = confidenceLevel;
    currie_input->additional_uncertainty = 0.0f;  // TODO: can we get the DRFs contribution to form this?
    
    m_currentCurrieInput = currie_input;
    const DetectionLimitCalc::CurrieMdaResult currie_result = DetectionLimitCalc::currie_mda_calc( *currie_input );
    m_currentCurrieResults = make_shared<DetectionLimitCalc::CurrieMdaResult>( currie_result );
    
    
    // Calculating the deconvolution-style limit is fairly CPU intensive, so we will only computer
    //  it when its what the user actually wants.
    const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
    if( !currieMethod )
    {
      const float live_time = m_currentCurrieInput->spectrum->live_time();
      
      //if( !m_currentNuclide )
      //  throw runtime_error( "Please enter a nuclide." );
      
      //if( (energyIndex < 0) || (energyIndex >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
      //  throw runtime_error( "Please select gamma energy." );
      
      shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
      if( drf && !drf->isValid() )
        drf.reset();
      
      // TODO: we could modify the deconvolution-style computations to not need a DRF for activity limits.
      // We need to calculate deconvolution-style limit
      if( !drf )
        throw runtime_error( WString::tr("dls-err-select-det").toUTF8() );
        
      m_fitFwhmBtn->setHidden( !drf || drf->hasResolutionInfo() );
      if( !drf || !drf->hasResolutionInfo() )
        throw runtime_error( WString::tr("dls-err-no-fwhm-info").toUTF8() );
      
      
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
        throw runtime_error( WString::tr("dls-err-neg-distance").toUTF8() );
      
      DetectionLimitCalc::DeconRoiInfo roiInfo;
      roiInfo.roi_start = m_lowerRoi->value(); // Will be rounded to nearest channel edge.
      roiInfo.roi_end = m_upperRoi->value();// Will be rounded to nearest channel edge.
      
        
      const int continuumTypeIndex = m_continuumType->currentIndex();
      switch( continuumTypeIndex )
      {
        case 0: 
          roiInfo.continuum_type = PeakContinuum::OffsetType::Linear;
          break;
          
        case 1:
          roiInfo.continuum_type = PeakContinuum::OffsetType::Quadratic;
          break;
          
        default:
          assert( 0 );
          throw std::logic_error( "Invalid continuuuum type selected" );
      }//switch( continuumTypeIndex )
      
      
      const int continuumPriorIndex = m_continuumPrior->currentIndex();
      switch( continuumPriorIndex )
      {
        case 0: // "Unknown or Present"
          roiInfo.cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::Floating;
          break;
          
        case 1: // "Not Present"
          roiInfo.cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::FixedByFullRange;
          break;
          
        case 2: // "Cont. from sides"
          roiInfo.cont_norm_method = DetectionLimitCalc::DeconContinuumNorm::FixedByEdges;
          break;
          
        default:
          assert( 0 );
          throw std::logic_error( "Invalid continuuuum prior selected" );
      }//switch( continuumPriorIndex )
      
      if( roiInfo.cont_norm_method == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
      {
        // `roiInfo.num_*_side_channels` only used if `cont_norm_method` is `DeconContinuumNorm::FixedByEdges`.
        roiInfo.num_lower_side_channels = static_cast<int>( m_numSideChannel->value() );
        roiInfo.num_upper_side_channels = roiInfo.num_lower_side_channels;
      }else
      {
        roiInfo.num_lower_side_channels = roiInfo.num_upper_side_channels = 0;
      }
      
      if( !m_currentNuclide )
      {
        DetectionLimitCalc::DeconRoiInfo::PeakInfo peakInfo;
        peakInfo.energy = energy;
        peakInfo.fwhm = m_fwhm->value();
        peakInfo.counts_per_bq_into_4pi = live_time; //Put in BR of 1
        roiInfo.peak_infos.push_back( peakInfo );
      }else
      {
        if( m_allGammasInRoi )
        {
          const double gammaFwhm = m_fwhm->value();
          const bool isDrfFwhm = (fabs(gammaFwhm - drf->peakResolutionFWHM(energy)) < 0.01);
          
          for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
          {
            const double thisEnergy = m_photoPeakEnergiesAndBr[i].first;
            const double thisBr = m_photoPeakEnergiesAndBr[i].second;
            
            if( (thisEnergy >= roi_lower_energy) && (thisEnergy <= roi_upper_energy) )
            {
              DetectionLimitCalc::DeconRoiInfo::PeakInfo peakInfo;
              peakInfo.energy = thisEnergy;
              // if user has left FWHM to detector predicted value, then consult the DRF for each peak
              peakInfo.fwhm = isDrfFwhm ? drf->peakResolutionFWHM(thisEnergy) : gammaFwhm;
              peakInfo.counts_per_bq_into_4pi = live_time * thisBr;//must have effects of shielding already accounted for, but not air atten, or det intrinsic eff
              roiInfo.peak_infos.push_back( peakInfo );
            }
          }//for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
        }else
        {
          const double br = m_photoPeakEnergiesAndBr[energyIndex].second;
          
          DetectionLimitCalc::DeconRoiInfo::PeakInfo peakInfo;
          peakInfo.energy = energy;
          peakInfo.fwhm = m_fwhm->value();
          peakInfo.counts_per_bq_into_4pi = live_time * br;//must have effects of shielding already accounted for, but not air atten, or det intrinsic eff
          roiInfo.peak_infos.push_back( peakInfo );
        }
      }//if( !m_currentNuclide )
      
      auto convo_input = make_shared<DetectionLimitCalc::DeconComputeInput>();
      convo_input->distance = distance;
      convo_input->activity = 0.0; //This will need to be varied
      convo_input->include_air_attenuation = true;
      convo_input->shielding_thickness = 0.0;
      convo_input->measurement = hist;
      convo_input->drf = drf;
      convo_input->roi_info.push_back( roiInfo );
      
      double min_act = 0.0;
      double max_act = 1E3*PhysicalUnits::ci;
      
      {// begin estimate range we should search for deconvolution
        double src_gammas_per_bq = 0.0, gammas_per_bq = 1.0;
        for( const DetectionLimitCalc::DeconRoiInfo::PeakInfo &peak_info : roiInfo.peak_infos )
          src_gammas_per_bq += peak_info.counts_per_bq_into_4pi;
        
        if( drf
           && (distance >= 0.0)
           && (src_gammas_per_bq > 0) 
           && m_currentCurrieInput->spectrum )
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
          const double counts_per_bq_into_4pi = src_gammas_per_bq * shield_transmission;
          const double counts_per_bq_into_4pi_with_air = air_transmission * counts_per_bq_into_4pi;
          const double counts_4pi = fixed_geom ? counts_per_bq_into_4pi : counts_per_bq_into_4pi_with_air;

          gammas_per_bq = counts_4pi * det_eff;
        }//if( drf )
        
        // We want the limits going into `DetectionLimitCalc::get_activity_or_distance_limits(...)`
        //  to definitely cover the entire possible activity range, so we will exaggerate the
        //  expected range from Currie-style limit
        //  The value of 5 is totally arbitrary, and I dont know what is actually a good range yet
        const double diff_multiple = 5.0;
        
        if( currie_result.source_counts > currie_result.decision_threshold )
        {
          // There is enough excess counts to reliably detect this activity
          const double lower_act = currie_result.lower_limit / gammas_per_bq;
          const double upper_act = currie_result.upper_limit / gammas_per_bq;
          const double nominal_act = currie_result.source_counts / gammas_per_bq;
          assert( lower_act <= nominal_act );
          assert( upper_act >= nominal_act );
          
          const double lower_diff = fabs(nominal_act - lower_act);
          const double upper_diff = fabs(upper_act - nominal_act);
          
          min_act = max( 0.0, (nominal_act - diff_multiple*lower_act) );
          max_act = max( 1.0/gammas_per_bq, (nominal_act + diff_multiple*upper_diff) );
        }else if( currie_result.upper_limit < 0 )
        {
          // There are a lot fewer counts in the peak region than predicted from the sides
          // We will just set the activity based on the Poisson uncertainty of the peak region
          min_act = 0.0;
          const double poison_uncert = sqrt(currie_result.peak_region_counts_sum);
          max_act = max( 1.0/gammas_per_bq, diff_multiple*poison_uncert/gammas_per_bq );
        }else
        {
          // No signal was detected, but we can estimate the minimum detectable activity
          const double simple_mda = currie_result.upper_limit / gammas_per_bq;
          min_act = 0.0;
          max_act = diff_multiple*simple_mda;
        }//if( detected signal ) / else / else
      }// end estimate range we should search for deconvolution
      
      m_currentDeconInput = convo_input;
      
      const bool is_dist_limit = false;
      const bool use_curie = use_curie_units();
      cout << "Will search between " << PhysicalUnits::printToBestActivityUnits(min_act, 3, use_curie)
      << " and " << PhysicalUnits::printToBestActivityUnits(max_act, 3, use_curie) << endl;
      
      const DetectionLimitCalc::DeconActivityOrDistanceLimitResult decon_result
                     = DetectionLimitCalc::get_activity_or_distance_limits( confidenceLevel, 
                                                                      convo_input, is_dist_limit,
                                                                      min_act, max_act, use_curie );
    
      m_currentDeconResults = make_shared<DetectionLimitCalc::DeconActivityOrDistanceLimitResult>( decon_result );
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
  // This function is fairly forgiving to the URI, in terms of not requiring
  //  all fields, and trooping along when it encounters invalid values.
  //  We might want to change this at some point.
  //
  // Example input:
  //  uri = "decon?nuc=Cs137&energy=661&dist=100 cm&..."

 
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw std::logic_error( "No SandiaDecayDataBase" );
  
  const bool undoWasSet = m_renderFlags.testFlag(RenderActions::AddUndoRedoStep);
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
 
  string host_str, path_str, query_str, fragment_str;
  AppUtils::split_uri( uri, host_str, path_str, query_str, fragment_str );
  
  const map<string,string> values = AppUtils::query_str_key_values( query_str );
  
  if( path_str.empty() )
    path_str = host_str; //A URI of "CURRIE?VER=1&..." will give a host value, and not a path value
  
  MethodIds methodIndex = MethodIds::Currie;
  if( SpecUtils::istarts_with(path_str, "CUR") )
    methodIndex = MethodIds::Currie;
  else if( SpecUtils::istarts_with(path_str, "DEC") )
    methodIndex = MethodIds::Deconvolution;
  else
    throw runtime_error( "DetectionLimitSimple::handleAppUrl: URI doesnt start with valid path" );
  
  if( m_methodGroup->checkedId() != static_cast<int>(methodIndex) )
  {
    m_methodGroup->setCheckedButton( m_methodGroup->button(static_cast<int>(methodIndex)) );
    handleMethodChanged();
  }//if( m_methodGroup->checkedId() != static_cast<int>(methodIndex) )
  
  auto qpos = values.find( "VER" );
  
  const string ver = ((qpos != end(values)) && !qpos->second.empty()) ? qpos->second : "1";
  if( (ver != "1") && !SpecUtils::istarts_with(ver, "1.") )
    throw runtime_error( "DetectionLimitSimple: invalid URI version" );
  
  qpos = values.find( "NUC" );
  const SandiaDecay::Nuclide *nuc = nullptr;
  if( qpos != end(values) )
    nuc = db->nuclide( qpos->second );
  
  m_nucEnterController->setNuclideText( nuc ? nuc->symbol : string() );
  
  // A sanity check for initial testing
#ifndef NDEBUG
  if( nuc )
  {
    const double age = m_nucEnterController->nuclideAge(); //mmm - using this maybe makes this check not totally independent???
    
    SandiaDecay::NuclideMixture mixture;
    const double dummy_activity = 0.001*SandiaDecay::curie;
    mixture.addAgedNuclideByActivity( nuc, dummy_activity, age );
    
    vector<SandiaDecay::EnergyRatePair> photons = mixture.xrays( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
    const vector<SandiaDecay::EnergyRatePair> gammas = mixture.gammas( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
    photons.insert( end(photons), begin(gammas), end(gammas) );
    size_t num_non_zero = 0;
    for( const auto &erp : photons )
      num_non_zero += ((erp.numPerSecond/dummy_activity) > std::numeric_limits<double>::epsilon());
    assert( m_photoPeakEnergiesAndBr.size() == num_non_zero );
  }else
  {
    assert( m_photoPeakEnergiesAndBr.empty() );
  }
#endif  //#ifndef NDEBUG
  
  
  qpos = values.find( "AGE" );
  if( (qpos != end(values)) && !qpos->second.empty() )
  {
    const string age = qpos->second;
    try
    {
      PhysicalUnits::stringToTimeDurationPossibleHalfLife(age, nuc ? nuc->halfLife : -1.0 );
      m_nucEnterController->setNuclideAgeTxt( age );
    }catch( std::exception &e )
    {
      cerr << "Failed to decode AGE string, '" << age << "', - but trooping along" << endl;
    }
  }//if( qpos != end(values) )
  
  
  
  float energy = 0.0f;
  bool setEnergy = false;
  qpos = values.find( "ENERGY" );
  if( nuc && (qpos != end(values)) && (stringstream(qpos->second) >> energy) )
  {
    // `m_photoPeakEnergiesAndBr` should have been updated
    double bestDiff = energy;
    size_t nearestEnergyIndex = m_photoPeakEnergiesAndBr.size();
    
    for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
    {
      const double diff = fabs( energy - m_photoPeakEnergiesAndBr[i].first );
      if( diff < bestDiff )
      {
        bestDiff = diff;
        nearestEnergyIndex = i;
      }
    }//for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
    
    if( (bestDiff < 0.1) && (nearestEnergyIndex < m_photoPeakEnergiesAndBr.size()) )
    {
      m_photoPeakEnergy->setCurrentIndex( static_cast<int>(nearestEnergyIndex) );
      energy = m_photoPeakEnergiesAndBr[nearestEnergyIndex].first;
      m_currentEnergy = energy;
      setEnergy = true;
      handleGammaChanged();
    }else
    {
      cerr << "Failed to find photopeak for energy=" << energy << endl;
    }
  }//if( URI contains ENERGY )
  
  float lowerRoi = m_lowerRoi->value(), upperRoi = m_upperRoi->value(), dummyFloat;
  qpos = values.find( "LROI" );
  if( (qpos != end(values)) && (stringstream(qpos->second) >> dummyFloat) )
    lowerRoi = dummyFloat;
  
  qpos = values.find( "UROI" );
  if( (qpos != end(values)) && (stringstream(qpos->second) >> dummyFloat) )
    upperRoi = dummyFloat;
  
  if( (energy < 10.0f) || ((lowerRoi < energy) && (upperRoi > energy)) )
  {
    m_lowerRoi->setValue( lowerRoi );
    m_upperRoi->setValue( upperRoi );
  }
  
  if( !setEnergy )
  {
    energy = 0.5*(lowerRoi + upperRoi);
    m_currentEnergy = energy;
  }
  
  qpos = values.find( "DIST" );
  if( qpos != end(values) )
  {
    try 
    {
      const double dist = PhysicalUnits::stringToDistance( qpos->second );
      if( dist >= 0.0 )
        m_distance->setValueText( WString::fromUTF8(qpos->second) );
    }catch( std::exception & )
    {
      cerr << "Failed to convert URI dist, '" << qpos->second << "' to a distance - trooping on" << endl;
    }
  }//if( URI has "DIST" )
  
  qpos = values.find( "CL" );
  if( qpos != end(values) )
  {
    int nsigma;
    if( (stringstream(qpos->second) >> nsigma) && (nsigma >= 1) && (nsigma <= 5) )
    {
      m_confidenceLevel->setCurrentIndex( nsigma - 1 );
    }else
    {
      cerr << "Failed to convert URI CL, '" << qpos->second << "' to a to an int between 1 and 5" << endl;
    }
  }//if( URI has CL )

  int continuumNormIndex = -1;
  qpos = values.find( "CONTNORM" );
  if( qpos != end(values) )
  {
    if( SpecUtils::istarts_with(qpos->second, "UNKN") )
      continuumNormIndex = 0;
    else if( SpecUtils::istarts_with(qpos->second, "NOS") )
      continuumNormIndex = 1;
    else if( SpecUtils::istarts_with(qpos->second, "FIX") )
      continuumNormIndex = 2;
    else
      cerr << "Unexpected 'CONTNORM' value: '" << qpos->second << "'" << endl;
  }//if( URI had continuum norm value )
    
  
  if( (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie))
     || (continuumNormIndex >= 0) )
  {
    qpos = values.find( "NSIDE" );
    if( qpos != end(values) )
    {
      int nside;
      if( (stringstream(qpos->second) >> nside) && (nside >= 1) && (nside <= 64) )
        m_numSideChannel->setValue( nside );
      else
        cerr << "Invalid 'NSIDE' value: '" << qpos->second << "'" << endl;
    }//if( URI has NSIDE )
  }//if( want value of number of side channels )
  
  bool setFwhm = false;
  qpos = values.find( "FWHM" );
  if( qpos != end(values) )
  {
    float fwhm;
    if( stringstream(qpos->second) >> fwhm )
    {
      m_fwhm->setValue( fwhm );
      setFwhm = true;
    }else
    {
      cerr << "Invalid 'FWHM' value: '" << qpos->second << "'" << endl;
    }
  }//if( have FWHM value )
  
  if( !setFwhm )
  {
    m_fwhm->setValue( -1.0f ); //force `handleUserChangedFwhm()` to reset FWHM value
    handleUserChangedFwhm();
  }
  
  if( m_methodGroup->checkedId() == static_cast<int>(MethodIds::Deconvolution) )
  {
    if( continuumNormIndex >= 0 )
      m_continuumPrior->setCurrentIndex( continuumNormIndex );
    
    
    qpos = values.find( "CONTTYPE" );
    if( qpos != end(values) )
    {
      int continuumTypeIndex = -1;
      if( SpecUtils::istarts_with(qpos->second, "LIN") )
        continuumTypeIndex = 0;
      else if( SpecUtils::istarts_with(qpos->second, "QUAD") )
        continuumTypeIndex = 1;
      else
        cerr << "Invalid 'CONTTYPE' value: '" << qpos->second << "'" << endl;
      
      if( continuumTypeIndex >= 0 )
        m_continuumType->setCurrentIndex( continuumTypeIndex );
    }//if( continuum type provided )
  }//if( methodIndex == MethodIds::Deconvolution )
  
  m_allGammasInRoi = true;
  qpos = values.find( "ALLGAMMA" );
  if( qpos != end(values) )
  {
    if( (qpos->second == "0") || SpecUtils::iequals_ascii(qpos->second, "NO") || SpecUtils::iequals_ascii(qpos->second, "FALSE") )
      m_allGammasInRoi = false;
  }//if( URI contains 'ALLGAMMA' )
  
  // Set render flags... JIC
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  
  if( !undoWasSet )
    m_renderFlags.clear(RenderActions::AddUndoRedoStep);
}//void handleAppUrl( std::string uri )


std::string DetectionLimitSimple::encodeStateToUrl() const
{
  string answer;
  
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  answer += currieMethod ? "CURRIE" : "DECON";
  
  answer += "?VER=1";
  
  const SandiaDecay::Nuclide *nuc = m_nucEnterController->nuclide();
  if( nuc )
  {
    answer += "&NUC=" + nuc->symbol;  //SpecUtils::to_upper_ascii(nuc->symbol)
    
    WString age = m_nucEnterController->nuclideAgeStr();
    if( !age.empty() )
    {
      answer += "&AGE=";
      
      // We need to make sure the age is always in English
      const string origAge = age.toUTF8();
      try
      {
        PhysicalUnits::stringToTimeDurationPossibleHalfLife(origAge, nuc->halfLife );
        
        // If here, we could interpret age without current localization - should be in English,
        //  use the string exactly.
        answer += origAge;
      }catch( std::exception & )
      {
        // Count number of orig significant figures
        int num_sig_fig = 0;
        for( const char &c : origAge )
          num_sig_fig += ((c >= '0') && (c <= '9'));
        
        num_sig_fig = std::max( 3, num_sig_fig ); //JIC we messed up, use at least 3 sig figs
        const double age = m_nucEnterController->nuclideAge();
        answer += PhysicalUnits::printToBestTimeUnits( age, num_sig_fig );
      }//try / catch to convert string to
    }//if( !age.empty() )
  }//if( m_nucEnterController->nuclide() )
  
  float energy = DetectionLimitSimple::photopeakEnergy();
  if( energy > 10.0 )
  {
    answer += "&ENERGY=" + SpecUtils::printCompact( energy, 6 );
  }else
  {
    // We wont explicitly put energy here, since tool assumes center of ROI
    assert( !nuc );
    energy = 0.5f*(m_lowerRoi->value() + m_upperRoi->value());
  }//if( have photopeak energy )
  
  answer += "&LROI=" + m_lowerRoi->valueText().toUTF8();
  answer += "&UROI=" + m_upperRoi->valueText().toUTF8();
  
  if( m_distance->validate() == WValidator::State::Valid )
  {
    // We dont (currently?) localize distance, so we should be good to directly
    //  embed the users text to the URL.
    answer += "&DIST=" + m_distance->text().toUTF8();
  }
  
  answer += "&CL=";
  const int clIndex = m_confidenceLevel->currentIndex();
  const ConfidenceLevel confidence = ConfidenceLevel(clIndex);
  
  switch( confidence )
  {
    case OneSigma:   answer += "1"; break;
    case TwoSigma:   answer += "2"; break;
    case ThreeSigma: answer += "3"; break;
    case FourSigma:  answer += "4"; break;
    case FiveSigma:  answer += "5"; break;
    case NumConfidenceLevel: assert(0); break;
  }//switch( confidence )
  
  const bool useSideChan = (m_continuumPrior->currentIndex() == 2);
  if( currieMethod || useSideChan )
    answer += "&NSIDE=" + std::to_string(m_numSideChannel->value());
  
  
  shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  if( drf && (!drf->isValid() || !drf->hasResolutionInfo()) )
    drf.reset();
  
  // Only include FWHM if its set by the user, or there is no DRF
  //  I'm not totally decided how to handle this quantity
  if( !drf || (fabs(m_fwhm->value() - drf->peakResolutionFWHM(energy)) > 0.1) )
    answer += "&FWHM=" + m_fwhm->valueText().toUTF8();
  
  if( !currieMethod )
  {
    answer += "&CONTNORM=";
    switch( m_continuumPrior->currentIndex() )
    {
      case 0: answer += "UNKNOWN"; break;
      case 1: answer += "NOSIG"; break;
      case 2: answer += "FIXED"; break;
        
      default:
        assert( 0 );
        throw logic_error( "Invalid m_continuumPrior" );
    }//switch( m_continuumPrior->currentIndex() )
    
    answer += "&CONTTYPE=";
    
    const int continuumTypeIndex = m_continuumType->currentIndex();
    switch( continuumTypeIndex )
    {
      case 0: answer += "LIN"; break;
      case 1: answer += "QUAD"; break;
      default:
        assert( 0 );
        throw std::logic_error( "Invalid continuuuum type selected" );
    }//switch( continuumTypeIndex )
  }//if( !currieMethod )
  
  if( !m_allGammasInRoi )
    answer += "&ALLGAMMA=0";
  
  return answer;
}//std::string encodeStateToUrl() const;


