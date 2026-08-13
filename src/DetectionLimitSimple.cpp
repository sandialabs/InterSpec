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
#include <thread>
#include <vector>
#include <sstream>


#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WImage>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WSpinBox>
#include <Wt/WCheckBox>
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
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
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
    
    return !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );
  }//bool use_curie_units()
  
}//namespace


DetectionLimitSimpleWindow::DetectionLimitSimpleWindow( Wt::WSuggestionPopup *materialSuggestion,
                                InterSpec *viewer )
: AuxWindow( WString::tr("window-title-simple-mda"),
            (AuxWindowProperties::TabletNotFullScreen
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse) )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;

  rejectWhenEscapePressed( true );

  // The dialog body is what scrolls, so contents that do not fit produce a scrollbar rather than
  //  pushing the footer off the bottom.  The class carries both the overflow and the height cap, in
  //  CSS, so it holds even when Wt's dialog layout has not re-measured - and on phone, where every
  //  C++ sizing call below is a no-op.  \sa DetectionLimitSimple.css
  contents()->addStyleClass( "SimpleMdaBody" );

  m_tool = new DetectionLimitSimple( materialSuggestion, viewer, contents() );
  // Deliberately no `setHeight(100%)`: with a definite height the tool's flex column would resolve
  //  against the squeezed body and its children would be asked to shrink, squashing the 200px chart,
  //  instead of the body scrolling.  Content-sized is also what makes the dialog layout's
  //  preferred-size measurement mean what we want.

  AuxWindow::addHelpInFooter( footer(), "simple-mda-dialog" );

#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton();
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

  // On fullscreen-on-phone, Close (MobileBackBtn) and QR both float left, so the one earlier
  // in the DOM is leftmost.  Put Close first there so it ends up on the far left;
  // for non-fullscreen layouts, QR goes between Help (left) and Close (right).
  if( !isPhone() )
    footer()->addWidget( qr_btn );
#endif //USE_QR_CODES

  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close"), true );
  closeButton->clicked().connect( this, &AuxWindow::hide );

#if( USE_QR_CODES )
  if( isPhone() )
    footer()->addWidget( qr_btn );
#endif
  
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





DetectionLimitSimple::DetectionLimitSimple( Wt::WSuggestionPopup *materialSuggestion,
                                 InterSpec *specViewer,
                                 Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
  m_viewer( specViewer ),
  m_materialSuggest( materialSuggestion ),
  m_spectrum( nullptr ),
  m_peakModel( nullptr ),
  m_resultTxt( nullptr ),
  m_warningTxt( nullptr ),
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
  m_advancedCb( nullptr ),
  m_advancedDiv( nullptr ),
  m_alpha( nullptr ),
  m_beta( nullptr ),
  m_alphaUserSet( false ),
  m_betaUserSet( false ),
  m_distanceUncert( nullptr ),
  m_prevDistanceUncert{},
  m_effUncert( nullptr ),
  m_advancedNote( nullptr ),
  m_systematicNote{},
  m_numFwhmWide( 2.5f ),
  m_lowerRoi( nullptr ),
  m_upperRoi( nullptr ),
  m_numSideChannelLabel( nullptr ),
  m_numSideChannel( nullptr ),
  m_fwhm( nullptr ),
  m_fwhmSuggestTxt( nullptr ),
  m_addFwhmBtn( nullptr ),
  m_selectDetectorBtn( nullptr ),
  m_isBackgroundDiv( nullptr ),
  m_isBackgroundSpectrum( nullptr ),
  m_isBackgroundHelpImg( nullptr ),
  m_planTimeCb( nullptr ),
  m_planTimeEdit( nullptr ),
  m_planTimeDiv( nullptr ),
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
  // Wording shared with the Detection Confidence Tool, so the two describe the same statistical
  //  choices in the same words.
  m_viewer->useMessageResourceBundle( "DetectionLimit" );

  addStyleClass( "DetectionLimitSimple" );

  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );
  const bool isPhone = m_viewer && m_viewer->isPhone();
  
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

  // Notes that qualify a limit that WAS produced - overlapping regions of interest combined, a
  //  profile crossing the threshold more than twice.  A sibling of the result rather than a third
  //  page of `m_chartErrMsgStack`, which is either-or (error page OR chart) and so has nowhere to
  //  show a warning beside a successful answer.  These were dropped entirely before Increment C.
  m_warningTxt = new WText( "", resultsDiv );
  m_warningTxt->addStyleClass( "MdaWarnMsg" );
  m_warningTxt->setInline( false );
  m_warningTxt->hide();
  
  // Now put the "more info..." link below here and to the right
  m_moreInfoButton = new WPushButton( resultsDiv );
  m_moreInfoButton->setText( WString::tr("dls-further-details-link") );
  m_moreInfoButton->setStyleClass( "LinkBtn MdaMoreInfoBtn" );
  m_moreInfoButton->clicked().connect( this, &DetectionLimitSimple::createMoreInfoWindow );
  m_moreInfoButton->setHiddenKeepsGeometry( true );
  m_moreInfoButton->hide();
 
  
  WContainerWidget *generalInput = new WContainerWidget( this );
  generalInput->addStyleClass( "GeneralInput" );
  
  
  WLabel *nucLabel = new WLabel( WString::tr("nuclide-label"), generalInput );
  m_nuclideEdit = new WLineEdit( generalInput );
  
  m_nuclideEdit->setMinimumSize( 30, WLength::Auto );
  nucLabel->setBuddy( m_nuclideEdit );
  
  WLabel *ageLabel = new WLabel( WString::tr("age-label"), generalInput );
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
  
  
  WLabel *gammaLabel = new WLabel( WString::tr("dls-gamma-label"), generalInput );
  gammaLabel->addStyleClass( "GridFirstCol GridThirdRow GridVertCenter" );
  
  
  m_photoPeakEnergy = new WComboBox( generalInput );
  m_photoPeakEnergy->activated().connect( this, &DetectionLimitSimple::handleGammaChanged );
  m_photoPeakEnergy->addStyleClass( "GridSecondCol GridThirdRow GridVertCenter PhotopeakComboBox" );
  
  // TODO: Add FWHM input.  Add text for DRF default, or the button to fit from data.
  //       when user changes this value - dont change ROI limits, just recalc deconv, and redraw either decon or Currie
  
  
  // Add Distance input
  WLabel *distanceLabel = new WLabel( WString::tr("distance-label"), generalInput );
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
  // Simple MDA always reports a one-sided upper limit, so the label says so rather than leaving the
  //  reader to guess which of the two products this confidence applies to.
  WLabel *confidenceLabel = new WLabel( WString::tr(isPhone ? "dl-conf-level-one-sided-short"
                                                            : "dl-conf-level-one-sided"), generalInput );
  confidenceLabel->addStyleClass( "GridFourthCol GridSecondRow GridVertCenter" );
  m_confidenceLevel = new WComboBox( generalInput );
  m_confidenceLevel->addStyleClass( "GridFifthCol GridSecondRow ClComboBox" );
  HelpSystem::attachToolTipOn( {confidenceLabel, m_confidenceLevel},
                              WString::tr("dl-conf-level-tt"), showToolTips );
  
  for( auto cl = ConfidenceLevel(0); cl < NumConfidenceLevel; cl = ConfidenceLevel(cl+1) )
  {
    // The values are unchanged; only the labels are.  The old "1σ (68.2%)" ... "5σ" wording claimed
    //  a sigma multiple these are not - they are central-normal probabilities applied as one-sided
    //  confidence levels - and gave 4σ/5σ no percentage at all.  See `dl-conf-level-tt`.
    const char *key = "";

    switch( cl )
    {
      case ConfidenceLevel::NinetyFivePercent:  key = "dl-conf-95";     break;
      case ConfidenceLevel::NinetyNinePercent:  key = "dl-conf-99";     break;
      case ConfidenceLevel::OneSigma:           key = "dl-conf-1sigma"; break;
      case ConfidenceLevel::TwoSigma:           key = "dl-conf-2sigma"; break;
      case ConfidenceLevel::ThreeSigma:         key = "dl-conf-3sigma"; break;
      case ConfidenceLevel::FourSigma:          key = "dl-conf-4sigma"; break;
      case ConfidenceLevel::FiveSigma:          key = "dl-conf-5sigma"; break;
      case ConfidenceLevel::NumConfidenceLevel: assert( 0 );            break;
    }//switch( cl )

    m_confidenceLevel->addItem( WString::tr(key) );
  }//for( loop over confidence levels )
  
  m_confidenceLevel->setCurrentIndex( ConfidenceLevel::NinetyFivePercent );
  m_confidenceLevel->activated().connect(this, &DetectionLimitSimple::handleConfidenceLevelChanged );
  
  
  // Add DRF select
  SpectraFileModel *specFileModel = m_viewer->fileManager()->model();
  m_detectorDisplay = new DetectorDisplay( m_viewer, specFileModel, generalInput );
  m_detectorDisplay->addStyleClass( "DetectorDisplay GridFourthCol GridThirdRow GridSpanTwoCol GridSpanTwoRows GridVertCenter" );
  m_viewer->detectorChanged().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  m_viewer->detectorModified().connect( boost::bind( &DetectionLimitSimple::handleDetectorChanged, this, boost::placeholders::_1 ) );
  
  
  
  WLabel *lowerRoiLabel = new WLabel( WString::tr(isPhone ? "dls-roi-lower-label-short" : "dls-roi-lower-label"), generalInput );
  lowerRoiLabel->addStyleClass( "GridFirstCol GridFourthRow GridVertCenter" );

  m_lowerRoi = new NativeFloatSpinBox( generalInput );
  m_lowerRoi->setSpinnerHidden();
  lowerRoiLabel->setBuddy( m_lowerRoi );
  m_lowerRoi->addStyleClass( "GridSecondCol GridFourthRow" );

  WLabel *upperRoiLabel = new WLabel( WString::tr(isPhone ? "dls-roi-upper-label-short" : "dls-roi-upper-label"), generalInput );
  upperRoiLabel->addStyleClass( "GridFirstCol GridFifthRow GridVertCenter" );
  
  m_upperRoi = new NativeFloatSpinBox( generalInput );
  m_upperRoi->setSpinnerHidden();
  upperRoiLabel->setBuddy( m_upperRoi );
  m_upperRoi->addStyleClass( "GridSecondCol GridFifthRow" );
  
  m_lowerRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  m_upperRoi->valueChanged().connect( this, &DetectionLimitSimple::handleUserChangedRoi );
  
  // Num Side Channel
  m_numSideChannelLabel = new WLabel( WString::tr(isPhone ? "dls-num-side-channel-label-short" : "dls-num-side-channel-label"), generalInput );
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
  
  m_isBackgroundDiv = new WContainerWidget( generalInput );
  m_isBackgroundDiv->addStyleClass( "BackCbDiv GridFirstCol GridSeventhRow GridVertCenter GridSpanTwoCol" );
  m_isBackgroundSpectrum = new WCheckBox( WString::tr(isPhone ? "dls-is-background-spectrum-cb-short" : "dls-is-background-spectrum-cb"), m_isBackgroundDiv );
  m_isBackgroundSpectrum->addStyleClass( "CbNoLineBreak" );
  m_isBackgroundSpectrum->checked().connect( this, &DetectionLimitSimple::handleNoSignalPresentChanged );
  m_isBackgroundSpectrum->unChecked().connect( this, &DetectionLimitSimple::handleNoSignalPresentChanged );

  {
    // One checkbox, one assertion - "this spectrum has no signal in it" - but it drives different
    //  machinery per method: zero Currie side channels, or the deconvolution's BackgroundReference
    //  measurement model.  `handleMethodChanged()` swaps the label and tooltip to say which.
    m_isBackgroundSpectrum->setWordWrap( false );
    m_isBackgroundHelpImg = new WImage( m_isBackgroundDiv );
    m_isBackgroundHelpImg->setImageLink(Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
    m_isBackgroundHelpImg->resize( 16, 16 );  //setStyleClass("Wt-icon");
    m_isBackgroundHelpImg->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
    HelpSystem::attachToolTipOn( m_isBackgroundHelpImg,
                                WString::tr("dls-is-background-spectrum-tt"), true,
                                HelpSystem::ToolTipPosition::Right,
                                HelpSystem::ToolTipPrefOverride::InstantAlways );
  }

  // Planned-measurement-time control on the right of row 7, beside the background checkbox in
  //  cols 1-2.  The deconvolution continuum controls used to share these cells and overlapped them;
  //  they now live on row 8.

  m_planTimeCb = new WCheckBox( WString::tr(isPhone ? "dl-plan-time-cb-short" : "dl-plan-time-cb"), generalInput );
  m_planTimeCb->addStyleClass( "CbNoLineBreak GridFourthCol GridSeventhRow GridVertCenter" );
  m_planTimeCb->setWordWrap( false );

  // These two hide under the deconvolution method until the spectrum is called a background
  //  reference.  Hiding them outright drops them out of the grid, and since the line edit is taller
  //  than the checkbox sharing the row, the row - and everything below it - jumps as the background
  //  checkbox is toggled.  Keeping the geometry reserves the cells either way, so only the contents
  //  appear and disappear.
  m_planTimeCb->setHiddenKeepsGeometry( true );
  m_planTimeCb->checked().connect( this, &DetectionLimitSimple::handlePlanTimeChanged );
  m_planTimeCb->unChecked().connect( this, &DetectionLimitSimple::handlePlanTimeChanged );

  m_planTimeDiv = new WContainerWidget( generalInput );
  m_planTimeDiv->addStyleClass( "ScaleEntryDiv GridFifthCol GridSeventhRow GridVertCenter" );
  m_planTimeDiv->setHiddenKeepsGeometry( true );   // the taller of the two - see above
  m_planTimeEdit = new WLineEdit( m_planTimeDiv );
  m_planTimeEdit->setEmptyText( WString::tr("dl-plan-time-empty-text") );
  m_planTimeEdit->setDisabled( true );
  {
    WRegExpValidator *scaleTimeValidator
              = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationRegex(), m_planTimeEdit );
    scaleTimeValidator->setFlags( Wt::MatchCaseInsensitive );
    m_planTimeEdit->setValidator( scaleTimeValidator );
  }
  m_planTimeEdit->changed().connect( this, &DetectionLimitSimple::handlePlanTimeChanged );
  m_planTimeEdit->blurred().connect( this, &DetectionLimitSimple::handlePlanTimeChanged );
  m_planTimeEdit->enterPressed().connect( this, &DetectionLimitSimple::handlePlanTimeChanged );

  {
    WImage *img = new WImage( m_planTimeDiv );
    img->setImageLink( Wt::WLink("InterSpec_resources/images/help_minimal.svg") );
    img->resize( 16, 16 );
    img->decorationStyle().setCursor( Wt::Cursor::WhatsThisCursor );
    HelpSystem::attachToolTipOn( img, WString::tr("dl-plan-time-tt"), true,
                                HelpSystem::ToolTipPosition::Left,
                                HelpSystem::ToolTipPrefOverride::InstantAlways );
  }

  // Pre-fill the disabled input with the current foreground's real time.
  {
    const shared_ptr<const SpecUtils::Measurement> hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( hist && (hist->real_time() > 0.0f) )
      m_planTimeEdit->setText( WString::fromUTF8(
          PhysicalUnits::printToBestTimeUnits( hist->real_time(), 3 ) ) );
    else
      m_planTimeCb->setDisabled( true );
  }
  
  m_continuumPriorLabel = new WLabel( WString::tr("dls-decon-cont-norm-label"), generalInput );
  m_continuumPriorLabel->addStyleClass( "GridFirstCol GridEighthRow GridVertCenter" );
  m_continuumPrior = new WComboBox( generalInput );
  // Only the two selectable treatments; the "as no signal" option is deprecated - it measured about
  //  40% coverage where 95% was claimed - and a stored state naming it is migrated to a background
  //  reference with a visible notice.  Order must match `continuum_norm_from_index`.
  m_continuumPrior->addItem( WString::tr("dl-cont-norm-floating") );
  m_continuumPrior->addItem( WString::tr("dl-cont-norm-edges-short") );
  assert( m_continuumPrior->count() == DetectionLimitCalc::num_selectable_continuum_norms() );
  m_continuumPrior->setCurrentIndex( 0 );
  m_continuumPrior->activated().connect( this, &DetectionLimitSimple::handleDeconPriorChange );
  m_continuumPrior->addStyleClass( "ContTypeCombo GridSecondCol GridEighthRow" );
  
  //m_continuumPriorLabel->setHiddenKeepsGeometry( true );
  //m_continuumPrior->setHiddenKeepsGeometry( true );
  m_continuumPriorLabel->hide();
  m_continuumPrior->hide();
  
  m_continuumTypeLabel = new WLabel( "Continuum Type:", generalInput );
  m_continuumTypeLabel->addStyleClass( "GridFourthCol GridEighthRow GridVertCenter" );
  m_continuumType = new WComboBox( generalInput );
  m_continuumType->addItem( WString::tr( PeakContinuum::offset_type_label_tr(PeakContinuum::OffsetType::Linear) ) );
  m_continuumType->addItem( WString::tr( PeakContinuum::offset_type_label_tr(PeakContinuum::OffsetType::Quadratic) ) );
  m_continuumType->setCurrentIndex( 0 );
  m_continuumType->activated().connect( this, &DetectionLimitSimple::handleDeconContinuumTypeChange );
  m_continuumType->addStyleClass( "GridFifthCol GridEighthRow" );
  m_continuumTypeLabel->setHidden( true );
  m_continuumType->setHidden( true );
  
  WContainerWidget *container = new WContainerWidget( generalInput );
  container->addStyleClass( "MethodSelect GridFirstCol GridNinthRow GridSpanFiveCol" );
  
  WLabel *methodLabel = new WLabel( WString::tr("dls-calc-method"), container);
  
  m_methodGroup = new WButtonGroup( container );
  WRadioButton *currieBtn = new Wt::WRadioButton( WString::tr("dls-currie-tab-title"), container );
  m_methodGroup->addButton(currieBtn, static_cast<int>(MethodIds::Currie) );
  
  WRadioButton *deconvBtn = new Wt::WRadioButton( WString::tr("dls-decon-tab-title"), container);
  m_methodGroup->addButton(deconvBtn, static_cast<int>(MethodIds::Deconvolution) );
  m_methodGroup->setCheckedButton( currieBtn );
  
  m_methodGroup->checkedChanged().connect( this, &DetectionLimitSimple::handleMethodChanged );

  // "Advanced" goes on the far right of this row.  An auto margin in the flex row does the pushing,
  //  so `justify-content: flex-start` can stay and the label + radios keep hugging the left edge.
  m_advancedCb = new WCheckBox( WString::tr(isPhone ? "dl-advanced-cb-short" : "dl-advanced-cb"),
                                container );
  m_advancedCb->addStyleClass( "CbNoLineBreak AdvancedCb" );
  m_advancedCb->setWordWrap( false );
  m_advancedCb->checked().connect( this, &DetectionLimitSimple::handleAdvancedToggled );
  m_advancedCb->unChecked().connect( this, &DetectionLimitSimple::handleAdvancedToggled );
  HelpSystem::attachToolTipOn( m_advancedCb, WString::tr("dl-advanced-tt"), showToolTips );

  m_methodDescription = new WText( WString::tr("dls-currie-desc"), generalInput );
  m_methodDescription->addStyleClass( "CalcMethodDesc GridSecondCol GridTenthRow GridSpanFourCol" );

  // The advanced statistical inputs; a sibling of `generalInput` rather than an eleventh grid row -
  //  see `m_advancedDiv`'s doc comment for why.
  m_advancedDiv = new WContainerWidget( this );
  m_advancedDiv->addStyleClass( "AdvancedInput" );
  m_advancedDiv->hide();  //Deliberately NOT setHiddenKeepsGeometry: this must take up no room.

  // Labels and fields are direct grid children, so the columns line up across both rows; wrapping
  //  them in per-pair divs would let each pair size independently and lose that alignment.
  WLabel *alphaLabel = new WLabel( WString::tr(isPhone ? "dl-alpha-label-short" : "dl-alpha-label"),
                                   m_advancedDiv );
  m_alpha = new NativeFloatSpinBox( m_advancedDiv );
  m_alpha->setSpinnerHidden();
  // A probability, not a percent.  The upper bound is 0.5 because at or above it the "threshold"
  //  would sit at or below the expected background, and the arithmetic stops meaning what it says.
  //  The lower bound is 1E-7 because below that the normal approximation the whole Currie
  //  formulation rests on is not worth much - and it is under the smallest value the
  //  confidence-level combo can produce (1 - 0.999999426696856 = 5.7E-7).
  m_alpha->setRange( 1.0E-7f, 0.4999f );
  m_alpha->setValue( static_cast<float>( 1.0 - currentConfidenceLevel() ) );
  alphaLabel->setBuddy( m_alpha );
  m_alpha->valueChanged().connect( this, &DetectionLimitSimple::handleAlphaChanged );
  HelpSystem::attachToolTipOn( {alphaLabel, m_alpha}, WString::tr("dl-alpha-tt"), showToolTips );

  WLabel *betaLabel = new WLabel( WString::tr(isPhone ? "dl-beta-label-short" : "dl-beta-label"),
                                  m_advancedDiv );
  m_beta = new NativeFloatSpinBox( m_advancedDiv );
  m_beta->setSpinnerHidden();
  m_beta->setRange( 1.0E-7f, 0.4999f );  //Same reasoning as alpha, above.
  m_beta->setValue( static_cast<float>( 1.0 - currentConfidenceLevel() ) );
  betaLabel->setBuddy( m_beta );
  m_beta->valueChanged().connect( this, &DetectionLimitSimple::handleBetaChanged );
  HelpSystem::attachToolTipOn( {betaLabel, m_beta}, WString::tr("dl-beta-tt"), showToolTips );

  WLabel *distUncertLabel = new WLabel(
            WString::tr(isPhone ? "dl-dist-uncert-label-short" : "dl-dist-uncert-label"),
            m_advancedDiv );
  m_distanceUncert = new WLineEdit( m_advancedDiv );
  m_distanceUncert->setEmptyText( WString::tr("dl-dist-uncert-empty-text") );
  distUncertLabel->setBuddy( m_distanceUncert );
  m_distanceUncert->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_distanceUncert->setAttributeValue( "autocorrect", "off" );
  m_distanceUncert->setAttributeValue( "spellcheck", "off" );
#endif
  {
    // Same grammar as the distance field above it.  Note `PhysicalUnits::stringToDistance` requires
    //  a unit *unless* the value is exactly "0" - which is what makes "0" the natural spelling of
    //  "none", and stops a bare "1" being silently read as some unit.
    WRegExpValidator *distUncertVal
                = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
    distUncertVal->setFlags( Wt::MatchCaseInsensitive );
    m_distanceUncert->setValidator( distUncertVal );
  }
  m_distanceUncert->changed().connect( this, &DetectionLimitSimple::handleSystematicUncertChanged );
  m_distanceUncert->enterPressed().connect( this, &DetectionLimitSimple::handleSystematicUncertChanged );
  HelpSystem::attachToolTipOn( {distUncertLabel, m_distanceUncert},
                              WString::tr("dl-dist-uncert-tt"), showToolTips );

  WLabel *effUncertLabel = new WLabel(
            WString::tr(isPhone ? "dl-eff-uncert-label-short" : "dl-eff-uncert-label"),
            m_advancedDiv );
  m_effUncert = new WLineEdit( m_advancedDiv );
  m_effUncert->setEmptyText( WString::tr("dl-eff-uncert-empty-text") );
  effUncertLabel->setBuddy( m_effUncert );
  {
    // Not mandatory, so an empty field validates - blank means "none".  The 99.9 cap is only a first
    //  pass; the *combined* systematic is what has to stay under 100%, and
    //  `currentSystematicUncertainty()` is where that is enforced.
    WDoubleValidator *effUncertVal = new WDoubleValidator( 0.0, 99.9, this );
    m_effUncert->setValidator( effUncertVal );
  }
  m_effUncert->changed().connect( this, &DetectionLimitSimple::handleSystematicUncertChanged );
  m_effUncert->enterPressed().connect( this, &DetectionLimitSimple::handleSystematicUncertChanged );
  HelpSystem::attachToolTipOn( {effUncertLabel, m_effUncert},
                              WString::tr("dl-eff-uncert-tt"), showToolTips );

  m_advancedNote = new WText( WString::tr("dls-advanced-decon-note"), m_advancedDiv );
  m_advancedNote->addStyleClass( "AdvancedNote" );
  m_advancedNote->setInline( false );
  m_advancedNote->hide();  //Only shown under the Deconvolution method; \sa handleMethodChanged

  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void DetectionLimitSimple::init()



void DetectionLimitSimple::roiDraggedCallback( double new_roi_lower_energy,
                 double new_roi_upper_energy,
                 double new_roi_px,
                 double original_roi_lower_energy,
                 const std::string &spec_type,
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
  assert( m_methodGroup->checkedId() != static_cast<int>(MethodIds::Currie) );

  // One function owns the whole show/hide table; see `handleMethodChanged()`.
  handleMethodChanged();
}//void handleDeconPriorChange()


void DetectionLimitSimple::handleNoSignalPresentChanged()
{
  // The checkbox now applies to both methods, so the whole show/hide table lives in one function
  //  rather than being duplicated here and drifting.
  handleMethodChanged();
}//void handleNoSignalPresentChanged()


void DetectionLimitSimple::handleDeconContinuumTypeChange()
{
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  scheduleRender();
}//void handleDeconContinuumTypeChange()


void DetectionLimitSimple::handlePlanTimeChanged()
{
  const bool checked = m_planTimeCb->isChecked();
  m_planTimeEdit->setEnabled( checked );

  // Whenever the field is empty (whitespace counts), or the checkbox is off,
  // repopulate it with the current foreground's real time so the displayed
  // value is always meaningful and `currentEffectiveForeground()` won't see an
  // empty string.
  const bool field_empty
              = SpecUtils::trim_copy( m_planTimeEdit->text().toUTF8() ).empty();
  if( !checked || field_empty )
  {
    const shared_ptr<const SpecUtils::Measurement> hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( hist && (hist->real_time() > 0.0f) )
      m_planTimeEdit->setText( WString::fromUTF8(
          PhysicalUnits::printToBestTimeUnits( hist->real_time(), 3 ) ) );
    else if( !checked )
      m_planTimeEdit->setText( "" );
  }

  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;
  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  scheduleRender();
}//void handlePlanTimeChanged()


double DetectionLimitSimple::currentPlanTimeSeconds() const
{
  // Hidden, unchecked, or blank all mean "not asked for"; none of them is an error.
  if( !m_planTimeCb || !m_planTimeEdit || m_planTimeCb->isHidden() || !m_planTimeCb->isChecked() )
    return 0.0;

  const string txt = SpecUtils::trim_copy( m_planTimeEdit->text().toUTF8() );
  if( txt.empty() )
    return 0.0;

  const double t = PhysicalUnits::stringToTimeDuration( txt ); //throws on parse error
  if( t <= 0.0 )
    throw runtime_error( WString::tr("dl-err-bad-plan-time").toUTF8() );

  return t;
}//currentPlanTimeSeconds()


DetectionLimitCalc::DeconMeasurementModel DetectionLimitSimple::currentMeasurementModel() const
{
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));

  // Under Currie the checkbox means "use the peak region itself to estimate the background" (zero
  //  side channels); there is no future-measurement model there.  Under the deconvolution method
  //  the same assertion is what BackgroundReference is.
  if( currieMethod || !m_isBackgroundSpectrum || !m_isBackgroundSpectrum->isChecked() )
    return DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum;

  return DetectionLimitCalc::DeconMeasurementModel::BackgroundReference;
}//currentMeasurementModel()


DetectionLimitCalc::PlannedMeasurement DetectionLimitSimple::currentEffectiveForeground() const
{
  const shared_ptr<const SpecUtils::Measurement> hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);

  return DetectionLimitCalc::plan_measurement( hist, currentPlanTimeSeconds(),
                                               currentMeasurementModel() );
}//currentEffectiveForeground()


DetectionLimitSimple::~DetectionLimitSimple()
{
  //nothing to do here
}//~DoseCalcWidget()


void DetectionLimitSimple::handleMethodChanged()
{
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  const bool isBackground = m_isBackgroundSpectrum->isChecked();

  // The background checkbox applies to BOTH methods now - it is the same assertion either way,
  //  "this spectrum has no signal in it" - but it drives different machinery, so only the wording
  //  changes.  Under Currie it zeroes the side channels; under deconvolution it selects the
  //  BackgroundReference measurement model.
  const bool isPhone = (m_viewer && m_viewer->isPhone());
  m_isBackgroundSpectrum->setText( WString::tr( currieMethod
                    ? (isPhone ? "dls-is-background-spectrum-cb-short" : "dls-is-background-spectrum-cb")
                    : "dl-model-backref-cb" ) );
  if( m_isBackgroundHelpImg )
  {
    // On mobile `attachToolTipOn` connects a NEW clicked() handler rather than replacing a qTip,
    //  and nothing disconnects the old ones - so without this every method/background change would
    //  add another stacked dialog, each showing the text captured when it was attached.
    HelpSystem::removeToolTipOn( m_isBackgroundHelpImg );
    HelpSystem::attachToolTipOn( m_isBackgroundHelpImg,
                                WString::tr( currieMethod ? "dls-is-background-spectrum-tt"
                                                          : "dl-model-backref-tt" ),
                                true, HelpSystem::ToolTipPosition::Right,
                                HelpSystem::ToolTipPrefOverride::InstantAlways );
  }

  m_continuumPriorLabel->setHidden( currieMethod );
  m_continuumPrior->setHidden( currieMethod );

  const bool useSideChan = !currieMethod
        && (DetectionLimitCalc::continuum_norm_from_index( m_continuumPrior->currentIndex() )
            == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges);

  m_numSideChannelLabel->setHidden( currieMethod ? isBackground : !useSideChan );
  m_numSideChannel->setHidden( currieMethod ? isBackground : !useSideChan );
  m_continuumTypeLabel->setHidden( currieMethod || useSideChan );
  m_continuumType->setHidden( currieMethod || useSideChan );

  // T_s means "the dwell I am asking about".  Under Currie that is always answerable.  Under the
  //  deconvolution method it is only answerable once the spectrum has been asserted to be a
  //  background reference - projecting the spectrum in hand and then bounding the signal in it is
  //  circular.  \sa DetectionLimitCalc::DeconMeasurementModel
  const bool planTimeApplies = currieMethod || isBackground;
  m_planTimeCb->setHidden( !planTimeApplies );
  m_planTimeDiv->setHidden( !planTimeApplies );

  m_methodDescription->setText( WString::tr(currieMethod ? "dls-currie-desc" : "dls-decon-desc") );

  // The advanced inputs stay visible under the deconvolution method but go disabled, with a note
  //  saying why - hiding them would make the checkbox look broken, and the values are still part of
  //  the state.  They are Currie-method quantities, and the deconvolution limit does not consume
  //  them.
  // TODO: `DetectionLimitCalc::decon_characteristic_limits()` (on the feature/DeconLdLc branch)
  //       computes L_c/L_d for the deconvolution method from an alpha and a beta.  Wire these two
  //       fields into it, and add a systematic term to `DeconComputeInput`, when that calculation
  //       reaches the GUI.
  m_alpha->setDisabled( !currieMethod );
  m_beta->setDisabled( !currieMethod );
  m_distanceUncert->setDisabled( !currieMethod );
  m_effUncert->setDisabled( !currieMethod );
  m_advancedNote->setHidden( currieMethod );

  // Section visibility follows the checkbox, and is set here as well as in `handleAdvancedToggled()`
  //  because `handleAppUrl()` drives this one function to bring all dependent visibility into line
  //  after decoding.
  m_advancedDiv->setHidden( !m_advancedCb->isChecked() );

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
  assert( m_currentNuclide == nuc );
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
    shared_ptr<const SpecUtils::Measurement> hist;
    try
    {
      // Under a background reference this is the UNSCALED spectrum - the counts the likelihood
      //  actually sees, with the projection carried as an exposure instead.
      hist = currentEffectiveForeground().decon;
    }catch( std::exception & )
    {
      // Bad planned-time string; fall back to raw foreground for display.
      hist = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    }
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
  
  const double dummy_activity = 0.001*SandiaDecay::curie;
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
  
  const float fwhm = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(energyToSelect) );
  
  size_t transition_index = 0;
  const SandiaDecay::Transition *transition = nullptr;
  PeakDef::SourceGammaType sourceGammaType = PeakDef::SourceGammaType::NormalGamma;
  PeakDef::findNearestPhotopeak( nuc, energyToSelect, 4.0*drfSigma, false, -1.0,
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
  
  
  double min_scale_delta_e = DBL_MAX;
  int currentIndex = -1;
   
  for( size_t i = 0; i < photons.size(); ++i )
  {
    double energy = photons[i].energy;
    const double intensity = photons[i].numPerSecond / dummy_activity;
    
    if( intensity > std::numeric_limits<double>::epsilon() )
    {
      const double delta_e = fabs( energyToSelect - energy );
      const double scale_delta_e = (0.1*fwhm + delta_e) / intensity;
      
      energy = floor(10000.0*energy + 0.5)/10000.0;
      
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
      
      if( scale_delta_e < min_scale_delta_e )
      {
        min_scale_delta_e = scale_delta_e;
        currentIndex = static_cast<int>( m_photoPeakEnergiesAndBr.size() - 1 );
      }
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
    
    if( (gamma_index >= 0) && (gamma_index < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
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
  // Until the user edits them, the two error rates follow the confidence level as 1 - CL.  Only the
  //  displayed value is updated here; while un-latched the calculation still receives the sentinel,
  //  so it uses the confidence level itself rather than this rounded-for-display copy of it.
  const float complement = static_cast<float>( 1.0 - currentConfidenceLevel() );
  if( !m_alphaUserSet )
    m_alpha->setValue( complement );
  if( !m_betaUserSet )
    m_beta->setValue( complement );

  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleConfidenceLevelChanged()


void DetectionLimitSimple::handleAdvancedToggled()
{
  m_advancedDiv->setHidden( !m_advancedCb->isChecked() );

  // Wt's dialog layout re-measures only when something schedules an adjust; a widget appearing
  //  inside the body does not.  This is the app's established way of asking for one - and it is
  //  `wApp`-scoped, so the tool does not need to know whether a dialog is around it.
  wApp->doJavaScript( wApp->javaScriptClass() + ".TriggerResizeEvent();" );

  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleAdvancedToggled()


void DetectionLimitSimple::handleAlphaChanged()
{
  // Latching here is what stops the value following the confidence-level combo from now on, and it
  //  is also what makes the ALPHA token appear in the state URI.
  m_alphaUserSet = true;

  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleAlphaChanged()


void DetectionLimitSimple::handleBetaChanged()
{
  m_betaUserSet = true;

  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleBetaChanged()


void DetectionLimitSimple::handleSystematicUncertChanged()
{
  // The distance-uncertainty field takes a length string, so an invalid entry is reverted the same
  //  way `handleDistanceChanged()` reverts the distance itself.  The efficiency field is guarded by
  //  a `WDoubleValidator`, so it needs no equivalent.
  const WString dist_uncert = m_distanceUncert->text();
  if( m_distanceUncert->validate() == WValidator::State::Valid )
    m_prevDistanceUncert = dist_uncert;
  else
    m_distanceUncert->setText( m_prevDistanceUncert );

  m_renderFlags |= DetectionLimitSimple::RenderActions::AddUndoRedoStep;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  scheduleRender();
}//void handleSystematicUncertChanged()


float DetectionLimitSimple::currentSystematicUncertainty( Wt::WString &note ) const
{
  // With "Advanced" off, nothing here reaches the calculation: the tool must give exactly the
  //  answers it gave before this section existed.
  if( !m_advancedCb || !m_advancedCb->isChecked() )
    return 0.0f;

  double dist_uncert = 0.0;
  string txt = m_distanceUncert->text().toUTF8();
  SpecUtils::trim( txt );
  if( !txt.empty() )
  {
    try
    {
      dist_uncert = PhysicalUnits::stringToDistance( txt );
    }catch( std::exception & )
    {
      throw runtime_error( WString::tr("dl-err-bad-dist-uncert").toUTF8() );
    }
  }//if( a distance uncertainty was entered )

  double eff_uncert = 0.0;
  txt = m_effUncert->text().toUTF8();
  SpecUtils::trim( txt );
  if( !txt.empty() )
  {
    if( !(stringstream(txt) >> eff_uncert) || (eff_uncert < 0.0) )
      throw runtime_error( WString::tr("dl-err-bad-eff-uncert").toUTF8() );
    eff_uncert /= 100.0;  //the field is a percent
  }//if( an efficiency uncertainty was entered )

  if( (dist_uncert <= 0.0) && (eff_uncert <= 0.0) )
    return 0.0f;

  // A fixed-geometry response has no distance in it, so 1/r^2 does not hold and a distance
  //  uncertainty simply does not propagate.  Say so, rather than quietly using or ignoring it.
  const shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  const bool inverse_square = !(drf && drf->isValid() && drf->isFixedGeometry());
  if( (dist_uncert > 0.0) && !inverse_square )
    note = WString::tr("dl-warn-dist-uncert-fixed-geom");

  double distance = 0.0;
  if( (dist_uncert > 0.0) && inverse_square )
  {
    try
    {
      distance = PhysicalUnits::stringToDistance( m_distance->text().toUTF8() );
    }catch( std::exception & )
    {
      distance = 0.0;
    }

    // The user explicitly asked for this term, so a missing distance is an error, not a silent drop.
    //  (The Currie limit itself does not need a distance - it is in counts - which is why this is
    //  checked here and not with the rest of the input.)
    if( distance <= 0.0 )
      throw runtime_error( WString::tr("dl-err-dist-uncert-no-distance").toUTF8() );
  }//if( a distance term applies )

  const double u_rel = DetectionLimitCalc::combine_systematic_uncertainty( distance,
                                    inverse_square ? dist_uncert : 0.0, eff_uncert, inverse_square );

  // Clamping would silently change the number the user entered, so this is a hard error.  Note the
  //  *detection limit* gives out earlier than this, at u = 1/k_beta (about 61% at beta = 0.05); that
  //  case still yields a valid upper limit, so it is only a warning - see `updateResult()`.
  if( u_rel >= 1.0 )
    throw runtime_error( WString::tr("dl-err-systematic-too-large")
                        .arg( SpecUtils::printCompact(100.0*u_rel, 3) ).toUTF8() );

  return static_cast<float>( u_rel );
}//float currentSystematicUncertainty( Wt::WString &note ) const


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


void DetectionLimitSimple::handleSpectrumChanged( const SpecUtils::SpectrumType type )
{
  // Background/secondary changes don't affect this tool's foreground-driven state, and
  // re-running the FWHM/ROI re-estimate on those would silently mutate the user's ROI.
  if( type != SpecUtils::SpectrumType::Foreground )
    return;

  // Update Scale-input behavior to track the new foreground.
  if( m_planTimeCb )
  {
    const shared_ptr<const SpecUtils::Measurement> hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    const bool has_real_time = (hist && (hist->real_time() > 0.0f));

    if( !has_real_time )
    {
      // Can't compute a sensible scale ratio - turn the feature off until a usable
      //  spectrum is loaded.
      m_planTimeCb->setChecked( false );
      m_planTimeCb->setDisabled( true );
      m_planTimeEdit->setDisabled( true );
      m_planTimeEdit->setText( "" );
    }else
    {
      m_planTimeCb->setDisabled( false );
      // Only auto-overwrite the field when the user is NOT actively scaling - otherwise
      //  preserve what they typed so they can compare the same target dwell across spectra.
      if( !m_planTimeCb->isChecked() )
      {
        m_planTimeEdit->setText( WString::fromUTF8(
            PhysicalUnits::printToBestTimeUnits( hist->real_time(), 3 ) ) );
      }
    }
  }

  // The new foreground may have a different resolution and/or DRF, so recenter the
  // ROI and re-estimate FWHM the same way handleGammaChanged() does - otherwise the
  // calc runs the new spectrum's counts through the previous spectrum's peak shape
  // and side-channel widths, and the limit comes out wildly wrong.
  if( m_currentEnergy > 10.0f )
  {
    const float fwhm
              = std::max( 0.1f, PeakSearchGuiUtils::estimate_FWHM_of_foreground(m_currentEnergy) );
    m_lowerRoi->setValue( m_currentEnergy - 0.5*m_numFwhmWide*fwhm );
    m_upperRoi->setValue( m_currentEnergy + 0.5*m_numFwhmWide*fwhm );
  }
  setFwhmFromEstimate();

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

  // Anything that qualifies the limit without preventing it - overlapping regions of interest
  //  combined, a profile crossing the threshold more than twice.  Each changes how the number
  //  should be read, so none may be dropped just because the calculation succeeded.  Set here
  //  rather than in `updateResult()`, which runs first and would have this reset wipe it.
  m_warningTxt->setText( "" );
  m_warningTxt->hide();
  {
    WString warning_text;

    if( m_currentDeconResults && !m_currentDeconResults->warnings.empty() )
    {
      string decon_warnings;
      for( const string &warning : m_currentDeconResults->warnings )
        decon_warnings += (decon_warnings.empty() ? "" : "<br />") + warning;
      warning_text = WString::fromUTF8( decon_warnings );
    }

    // A qualification raised while building the Currie input - currently only a distance uncertainty
    //  dropped because the detector response is fixed-geometry.  \sa currentSystematicUncertainty
    if( !m_systematicNote.empty() )
    {
      if( !warning_text.empty() )
        warning_text += WString::fromUTF8( "<br />" );
      warning_text += m_systematicNote;
    }

    // A systematic uncertainty at or above 1/k_beta leaves no finite detection limit - the true
    //  signal cannot be separated from the scale uncertainty however long you count.  The upper
    //  limit is unaffected, so this qualifies the answer rather than replacing it.
    //  \sa DetectionLimitCalc::currie_mda_calc, which returns -999 for exactly this case.
    if( m_currentCurrieResults && (m_currentCurrieResults->detection_limit <= -999.0f) )
    {
      if( !warning_text.empty() )
        warning_text += WString::fromUTF8( "<br />" );
      warning_text += WString::tr("dl-warn-no-detection-limit");
    }

    if( !warning_text.empty() )
    {
      m_warningTxt->setText( warning_text );
      m_warningTxt->show();
    }
  }
  m_moreInfoButton->hide();
  m_peakModel->setPeaks( vector<shared_ptr<const PeakDef>>{} );
  m_spectrum->removeAllDecorativeHighlightRegions();
  
  const double confidence_level = m_currentDeconResults ? m_currentDeconResults->confidenceLevel : currentConfidenceLevel();
  
  const string cl_str = DetectionLimitCalc::confidence_level_pct_str( confidence_level );
  
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));

  // The background checkbox applies to BOTH methods, so it is always visible - it asserts the same
  //  thing either way, and only the machinery behind it differs.  This used to assert it was
  //  visible only under Currie, which is no longer the invariant.  \sa handleMethodChanged
  assert( !m_isBackgroundDiv->isHidden() );

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
      const DetectorPeakResponse::EffGeometryType det_geom = drf ? drf->geometryType() : DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic;
      
      assert( !result
             || (result->input.num_lower_side_channels != 0)
             || (result->input.num_upper_side_channels != 0)
             || (result->input.num_lower_side_channels == result->input.num_upper_side_channels) );
      
      const bool assertedIsBackground = (result
                                         && (result->input.num_lower_side_channels == 0)
                                          && (result->input.num_upper_side_channels == 0));
      
      if( !result )
      {
        result_txt = WString::tr("dls-det-no-result");
      }else if( result->source_counts > result->decision_threshold )
      {
        assert( !assertedIsBackground );
        
        // There is enough excess counts that we would reliably detect this activity, so we will
        //  give the activity range.
        string lowerstr, upperstr, nomstr;
        
        // The quoted value is the ISO 11929-1:2019 best estimate, Formula (44) - the mean of the
        //  Gaussian after truncating at zero - not the primary result `source_counts`.  The decision
        //  above is still made on the primary result, per clause 10; only what is reported changes.
        //  The two agree above 4*u(y) and diverge as the signal weakens (~7% at the decision
        //  threshold for a well populated continuum, more at low counts).  The range is unchanged:
        //  it is already the truncated interval from the same distribution.
        if( gammas_per_bq > 0.0 )
        {
          const float lower_act = result->lower_limit / gammas_per_bq;
          const float upper_act = result->upper_limit / gammas_per_bq;
          const float nominal_act = result->best_estimate / gammas_per_bq;

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
          nomstr = SpecUtils::printCompact(result->best_estimate, 4);

          result_txt = WString::tr("dls-det-counts-with-range").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }
      }else if( result->upper_limit < 0 )
      {
        assert( !assertedIsBackground );
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
      
      // `currie_mda_calc` returns -999 when the systematic uncertainty leaves no finite detection
      //  limit.  Formatting that sentinel gives a nonsense negative activity, so the line is dropped
      //  entirely and `m_warningTxt` says why it is missing.  Only reachable since the "Advanced"
      //  section let a systematic uncertainty be entered at all.
      const bool haveDetLimit = (result && (result->detection_limit > -999.0f));

      WString mda_txt;
      if( haveDetLimit && (gammas_per_bq > 0.0) )
      {
        const double detection_act = result->detection_limit / gammas_per_bq;
        const string act = PhysicalUnits::printToBestActivityUnits( detection_act, 2, use_curie )
                          + DetectorPeakResponse::det_eff_geom_type_postfix( det_geom );

        mda_txt = WString::tr("dls-min-detectable-act").arg(act);
      }else if( haveDetLimit )
      {
        const string counts = SpecUtils::printCompact(result->detection_limit, 4);
        mda_txt = WString::tr("dls-min-detectable-counts").arg( counts );
      }//if( have a detection limit to report )

      // A background reference describes a measurement nobody has taken, so the minimum detectable
      //  line is the whole answer there - unless there isn't one, in which case fall back to the
      //  result line rather than showing nothing.
      WString full_result_txt;
      if( !result )
        full_result_txt = result_txt;
      else if( assertedIsBackground )
        full_result_txt = haveDetLimit ? mda_txt : result_txt;
      else if( haveDetLimit )
        full_result_txt = WString("{1}<br/>{2}").arg(result_txt).arg(mda_txt);
      else
        full_result_txt = result_txt;

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

      // `fit_peaks` come back in the exposure of the spectrum they are drawn over, so under a
      //  background reference their amplitudes are counts in the *reference*.  Every counts figure
      //  quoted here describes the measurement being predicted, so it is projected back up.
      //  Defined out here because a predicted sensitivity always has `foundLowerCl == false`
      //  (`get_activity_or_distance_limits` forces it), so the upper-limit-only branch below is
      //  the one that actually runs for exactly the case where the ratio is not one.
      const auto peak_counts = []( const shared_ptr<const DetectionLimitCalc::DeconComputeResults> &r )
        -> double {
          if( !r || (r->fit_peaks.size() != 1) )
            return -1.0;
          return r->exposure_ratio * r->fit_peaks[0].amplitude();
      };

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
          
          const double lower_limit_counts = peak_counts( result.lowerLimitResults );
          const double nominal_counts = peak_counts( result.overallBestResults );
          const double upper_limit_counts = peak_counts( result.upperLimitResults );
          
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
          
          chart_title = WString::tr("dls-chart-title-estimated-act").arg(nomstr);
          result_txt = WString::tr("dls-results-txt-estimated-act").arg(nomstr).arg(lowerstr).arg(upperstr).arg(cl_str);
        }//if( !m_currentNuclide ) / else
        
        assert( result.overallBestResults );
        if( result.overallBestResults )
          fit_peaks = result.overallBestResults->fit_peaks;
      }else if( result.foundLowerCl )
      {
        display_activity = 0.0; //result.foundLowerCl
        const string cl_txt = "Error: Didn't find " + cl_str + " CL activity";
        const string sum_txt = "Error: Didn't find " + cl_str + " CL activity";
        
        chart_title = WString::fromUTF8( cl_txt );
        result_txt = WString::fromUTF8( sum_txt );
        
        //fit_peaks = result.lowerLimitResults.fit_peaks;
      }else if( result.foundUpperCl )
      {
        display_activity = result.upperLimit;

        // A background-reference scan is a statement about a measurement nobody has taken, so it
        //  must not be worded as a bound on the loaded spectrum - the Detection Confidence Tool
        //  says "Predicted sensitivity: ..." for the same result, and this tool used to say "Less
        //  than X @95% CL", which reads as a bound on the spectrum in hand.
        //  \sa DetectionLimitCalc::DeconMeasurementModel::BackgroundReference
        const bool predicted = result.is_predicted_sensitivity;
        const string dwell = (result.sampleRealTime > 0.0)
              ? PhysicalUnits::printToBestTimeUnits( result.sampleRealTime, 2 )
              : string();

        if( !m_currentNuclide )
        {
          const double upper_limit_counts = peak_counts( result.upperLimitResults );
          const string upperstr = SpecUtils::printCompact(upper_limit_counts, 4);

          if( predicted && !dwell.empty() )
          {
            chart_title = WString::tr("dls-chart-title-predicted-counts").arg(upperstr).arg(cl_str);
            result_txt = WString::tr("dls-results-txt-predicted-counts")
                            .arg(upperstr).arg(cl_str).arg(dwell);
          }else
          {
            chart_title = WString::tr("dls-chart-title-upper-bound-counts").arg(upperstr).arg(cl_str);
            result_txt = WString::tr("dls-results-txt-upper-bound-counts").arg(upperstr).arg(cl_str);
          }
        }else
        {
          const string upperstr = PhysicalUnits::printToBestActivityUnits(result.upperLimit, 3, use_curie);

          if( predicted && !dwell.empty() )
          {
            chart_title = WString::tr("dls-chart-title-predicted-act").arg(upperstr).arg(cl_str);
            result_txt = WString::tr("dls-results-txt-predicted-act")
                            .arg(upperstr).arg(cl_str).arg(dwell);
          }else
          {
            chart_title = WString::tr("dls-chart-title-upper-bound-act").arg(upperstr).arg(cl_str);
            result_txt = WString::tr("dls-results-txt-upper-bound-act").arg(upperstr).arg(cl_str);
          }
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

  // The projected live time used to be appended here as a "(LT=...)" marker.  It cost a line of a
  //  three-line box that only has room for two, for something the user can already see: they typed
  //  the real time into the field beside the result, and the derived live time - with the ratio and
  //  the original dwell - is spelled out in the "further details" dialog.  \sa dls-scale-more-info-note

  // A background-reference result is a prediction even with no measurement time entered - no time
  //  means "another measurement just like this one" - so the band applies there too.
  {
    // A projected limit is a prediction about a measurement nobody has taken.  Quoting only its
    //  middle hides how far the answer can move, and by more than a plain scaling suggests - the
    //  predictive spread grows as sqrt(1+k) in the projection factor, where a scaling's does not
    //  grow at all.  So the band goes on screen beside the number.
    // Given as a multiple of the quoted limit rather than in counts or activity.  Counts and
    //  activity differ by a constant factor, so a fraction of the median is the same number in
    //  either - which keeps this correct without plumbing the conversion into this scope, and
    //  anchors the band on the figure actually on screen.
    string low_multiple, high_multiple;
    if( DetectionLimitCalc::projected_band_endpoints( m_currentProjectedLimit,
                                                      low_multiple, high_multiple ) )
    {
      const WString band = WString::tr("dls-result-projected-band")
          .arg( low_multiple )
          .arg( high_multiple );

      // Deliberately no tooltip here: this runs on every result update, and `attachToolTipOn`
      //  re-runs its qTip setup each call (and on mobile connects another `clicked()` handler),
      //  so attaching from an update path accumulates.  What the band means is explained in the
      //  help page's "Measurement time" section instead.
      m_resultTxt->setText( m_resultTxt->text() + band );
    }//if( a projected band was computed )
  }
}//void updateSpectrumDecorationsAndResultText()


double DetectionLimitSimple::currentConfidenceLevel() const
{
  double confidenceLevel = 0.95;
  
  const int clIndex = m_confidenceLevel->currentIndex();
  const ConfidenceLevel confidence = ConfidenceLevel(clIndex);
  
  switch( confidence )
  {
    case ConfidenceLevel::NinetyFivePercent:  confidenceLevel = 0.95;     break;
    case ConfidenceLevel::NinetyNinePercent:  confidenceLevel = 0.99;     break;
    case ConfidenceLevel::OneSigma:           confidenceLevel = 0.682689492137086; break;
    case ConfidenceLevel::TwoSigma:           confidenceLevel = 0.954499736103642; break;
    case ConfidenceLevel::ThreeSigma:         confidenceLevel = 0.997300203936740; break;
    case ConfidenceLevel::FourSigma:          confidenceLevel = 0.999936657516334; break;
    case ConfidenceLevel::FiveSigma:          confidenceLevel = 0.999999426696856; break;
    case ConfidenceLevel::NumConfidenceLevel: assert(0); break;
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

    // `fit_peaks` are in the reference spectrum's exposure, but the activity row above this one is
    //  the projected activity, so counts are projected to match.  Two rows of the same table in
    //  different exposures is worse than either choice on its own.  The continuum-area row further
    //  down stays in reference counts, and says so.
    counts *= result.exposure_ratio;
    value = SpecUtils::printCompact(counts, 5) + " counts";
    //value = PhysicalUnits::printValueWithUncertainty(counts, sqrt(uncert), 5);
    
    cell = table->elementAt( table->rowCount(), 0 );
    new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    new WText( value, cell );
    
    // The value is a Cash statistic, not a chi-square; `statistic_name` had been set since Stage 1
    //  and read nowhere, so this is the first place that label reaches a user.  The tooltip goes on
    //  both this row and the DOF row, because their ratio is the thing most likely to be misread.
    const WString stat_tt = WString::tr("dl-stat-is-cash-tt")
          .arg( WString::fromUTF8( result.statistic_name.empty() ? string("Cash")
                                                                 : result.statistic_name ) );

    label = WString("{1} &chi;<sup>2</sup>").arg(typestr);
    value = SpecUtils::printCompact(result.chi2, 4);
    cell = table->elementAt( table->rowCount(), 0 );
    WText *stat_label = new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    WText *stat_value = new WText( value, cell );
    HelpSystem::attachToolTipOn( {stat_label, stat_value}, stat_tt, true,
                                HelpSystem::ToolTipPosition::Left,
                                HelpSystem::ToolTipPrefOverride::InstantAlways );

    if( !is_best )
      return;

    label = WString::tr("dls-DOF");
    value = std::to_string( result.num_degree_of_freedom );
    cell = table->elementAt( table->rowCount(), 0 );
    WText *dof_label = new WText( label, cell );
    cell = table->elementAt( table->rowCount() - 1, 1 );
    WText *dof_value = new WText( value, cell );
    HelpSystem::attachToolTipOn( {dof_label, dof_value}, stat_tt, true,
                                HelpSystem::ToolTipPosition::Left,
                                HelpSystem::ToolTipPrefOverride::InstantAlways );

    if( result.fit_peaks.size() )
    {
      label = WString::tr("dls-continuum-area");
      const PeakDef &peak = result.fit_peaks.front();
      // CDF step types wont typically occur in detection limit context; if they do, use single peak
      const PeakDef *peak_ptr = &peak;
      const double cont_area = peak.continuum()->offset_integral( roi_start, roi_end, measurement, &peak_ptr, 1 );

      // Deliberately NOT projected: the continuum belongs to the reference spectrum, which is the
      //  one drawn on the chart beside this dialog.  Saying so keeps it from being read against the
      //  projected counts row above.
      value = SpecUtils::printCompact(cont_area, 5);
      if( result.exposure_ratio != 1.0 )
        value += " (in the reference spectrum)";
      
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
      const double geom_eff = drf->fractionalSolidAngle( drf->detectorDiameter(), distance + drf->detectorSetback() );
      
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
     && (!drf || !drf->isFixedGeometry()) 
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

  // Make the scale factor obvious in the more-info dialog so the user knows the
  //  numbers are for the scaled-to dwell, not the literal measurement live time.
  if( m_moreInfoWindow && m_planTimeCb && m_planTimeCb->isChecked() && !m_planTimeCb->isHidden() )
  {
    const shared_ptr<const SpecUtils::Measurement> raw_hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> scaled;
    try { scaled = currentEffectiveForeground().currie; } catch( std::exception & ){}

    if( raw_hist && scaled && (raw_hist->real_time() > 0.0f) && (scaled->real_time() > 0.0f) )
    {
      const double ratio = static_cast<double>(scaled->real_time()) / static_cast<double>(raw_hist->real_time());
      WText *note = new WText( WString::tr("dls-scale-more-info-note")
                              .arg( PhysicalUnits::printToBestTimeUnits(scaled->real_time(), 3) )
                              .arg( PhysicalUnits::printToBestTimeUnits(scaled->live_time(), 3) )
                              .arg( SpecUtils::printCompact(ratio, 3) )
                              .arg( PhysicalUnits::printToBestTimeUnits(raw_hist->real_time(), 3) ) );
      note->addStyleClass( "ScaleMoreInfoNote" );
      note->setInline( false );
      m_moreInfoWindow->contents()->insertWidget( 0, note );
    }
  }

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
  m_currentProjectedLimit = DetectionLimitCalc::ProjectedLimit{};

  // Parked here, rather than set on `m_warningTxt` directly, because
  //  `updateSpectrumDecorationsAndResultText()` runs after this function and clears that text.
  m_systematicNote = WString();

  try
  {
    m_fitFwhmBtn->hide();
    
    // Route through currentEffectiveForeground() so the empty-string-as-unrequested logic stays in
    //  one place and the calc here matches what the chart shows.  `.currie` is the projected
    //  spectrum, `.decon` is the unscaled reference under a background-reference model.
    DetectionLimitCalc::PlannedMeasurement eff;
    try
    {
      eff = currentEffectiveForeground();
    }catch( std::exception & )
    {
      throw runtime_error( WString::tr("dl-err-bad-plan-time").toUTF8() );
    }
    const std::shared_ptr<const SpecUtils::Measurement> hist = eff.currie;
    if( !hist || (hist->num_gamma_channels() < 7) )
      throw runtime_error( "No foreground spectrum loaded." );

    const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
    
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
    if( currieMethod && m_isBackgroundSpectrum->isChecked() )
    {
      currie_input->num_lower_side_channels = 0;
      currie_input->num_upper_side_channels = 0;
    }
    
    currie_input->detection_probability = confidenceLevel;

    // The two error rates and the systematic uncertainty come from the "Advanced" section, and apply
    //  to the Currie method only - which is what the section's note tells the user.  The gate on the
    //  method matters even though the deconvolution branch reports its own limit: it still uses this
    //  Currie result to seed the activity search range below, so applying these here would move the
    //  deconvolution answer, and a bad entry would abort it with an error about a field the UI says
    //  does not apply to it.
    //  Zero means "not specified", so the confidence level supplies all three roles - which is what
    //  happened before the section existed.
    if( currieMethod )
    {
      currie_input->alpha = m_alphaUserSet ? m_alpha->value() : 0.0;
      currie_input->beta = m_betaUserSet ? m_beta->value() : 0.0;
      currie_input->additional_uncertainty = currentSystematicUncertainty( m_systematicNote );
    }else
    {
      currie_input->alpha = currie_input->beta = 0.0;
      currie_input->additional_uncertainty = 0.0f;
    }

    m_currentCurrieInput = currie_input;
    const DetectionLimitCalc::CurrieMdaResult currie_result = DetectionLimitCalc::currie_mda_calc( *currie_input );
    m_currentCurrieResults = make_shared<DetectionLimitCalc::CurrieMdaResult>( currie_result );

    // How far the answer could move for a measurement nobody has taken yet.  Only meaningful when
    //  a dwell other than this spectrum's is being asked about, and only for the Currie method
    //  here - the deconvolution branch below computes its own, since it is the one whose limit is
    //  a profile scan rather than arithmetic.
    //  The region is drawn from the *unprojected* foreground: the Monte Carlo does the projecting,
    //  and handing it an already-scaled spectrum would project twice.
    const shared_ptr<const SpecUtils::Measurement> unprojected
                      = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );

    if( currieMethod && (eff.planned_real_time > 0.0) && unprojected )
    {
      DetectionLimitCalc::CurrieMdaInput reference_input = *currie_input;
      reference_input.spectrum = unprojected;

      // A few hundred realisations keeps the Monte Carlo's own sampling error well under a percent
      //  of the band, and the Currie limit is closed-form arithmetic, so this costs milliseconds.
      m_currentProjectedLimit = DetectionLimitCalc::currie_projected_limit( reference_input,
                                                          eff.planned_real_time, 1024 );
    }
    
    
    // Calculating the deconvolution-style limit is fairly CPU intensive, so we will only computer
    //  it when its what the user actually wants.
    if( !currieMethod )
    {
      // From the DECON spectrum, not the projected Currie one: `counts_per_bq_into_4pi` carries this
      //  live time, and `decon_compute_peaks` multiplies the trial signal by
      //  `sample_exposure/live_time` on top of it.  Taking the projected live time here would apply
      //  the projection twice - silently, and quadratically in the ratio.
      const float live_time = eff.decon->live_time();
      
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
      
      // Through the shared mapping rather than a local switch: this combo's indices used to run
      //  1 -> FixedByFullRange, 2 -> FixedByEdges, the reverse of the enum, with nothing guarding
      //  it.  `continuum_norm_from_index` is `static_assert`-backed and shared with both URL
      //  decoders, so the order now lives in exactly one place.
      roiInfo.cont_norm_method
          = DetectionLimitCalc::continuum_norm_from_index( m_continuumPrior->currentIndex() );

      if( roiInfo.cont_norm_method == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges )
      {
        // Four-plus-four side channels barely determine an offset and a slope, so a quadratic is
        //  not supportable here; the combo is disabled to match.
        roiInfo.continuum_type = PeakContinuum::OffsetType::Linear;
      }

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
      // The unscaled reference under a background-reference model; the planned dwell rides along as
      //  `sample_exposure`.  \sa DetectionLimitCalc::plan_measurement
      convo_input->measurement = eff.decon;
      convo_input->measurement_model = currentMeasurementModel();
      convo_input->sample_exposure = eff.sample_exposure;
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

          // The Currie result below is for the PROJECTED measurement, while `counts_per_bq_into_4pi`
          //  is at the reference exposure.  Bring them onto the same footing before dividing, or
          //  the seeded search range is wrong by the exposure ratio.  Exactly 1 unless a background
          //  reference is being projected.
          gammas_per_bq = counts_4pi * det_eff * eff.exposure_ratio;
        }//if( drf )

        // We want the limits going into `DetectionLimitCalc::get_activity_or_distance_limits(...)`
        //  to definitely cover the entire possible activity range, so we will exaggerate the
        //  expected range from Currie-style limit
        //  The value of 5 is totally arbitrary, and I dont know what is actually a good range yet
        const double diff_multiple = 50.0;

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
          
          //if( assertNoSignal )
          //{
          //  assert( currie_result.source_counts == 0.0f );
          //  const double niave_max = currie_result.peak_region_counts_sum / gammas_per_bq;
          //  max_act = diff_multiple * std::max( simple_mda, niave_max );
          //}
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

      // A background-reference result is a prediction about a measurement nobody has taken - and it
      //  is one whether or not a measurement time was entered, since no time entered means "another
      //  measurement just like this one".  Reporting only its middle says nothing about how far the
      //  answer could land from it, so the spread goes beside it.
      //
      //  Scored `JointWithReference` because that is what this mode reports: the reference and the
      //  sample against one continuum.  A sample-only prediction would be a different, looser
      //  quantity.  \sa DetectionLimitCalc::ProjectedLimitScoring
      if( (convo_input->measurement_model
             == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference)
         && convo_input->measurement && (convo_input->measurement->real_time() > 0.0f) )
      {
        // Users enter, and results quote, a REAL time.
        const double planned_real = (eff.planned_real_time > 0.0)
                       ? eff.planned_real_time
                       : static_cast<double>( convo_input->measurement->real_time() );

        // 128 realisations costs about half a second across the available cores for this tool's
        //  single region; the deconvolution limit itself is already the slow part of this path.
        //  The Detection Confidence Tool computes one of these per gamma line, so it does not run
        //  them inline.
        const unsigned hardware_threads = std::thread::hardware_concurrency();
        const size_t num_threads = (std::max)( 1u, hardware_threads ? hardware_threads : 4u );

        m_currentProjectedLimit = DetectionLimitCalc::decon_projected_limit( convo_input,
                                      confidenceLevel, planned_real, max_act, use_curie,
                                      DetectionLimitCalc::DeconLimitType::OneSidedUpperLimit,
                                      DetectionLimitCalc::ProjectedLimitScoring::JointWithReference,
                                      128, num_threads );
      }//if( a background-reference prediction )
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
  
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  
  auto qpos = values.find( "VER" );
  
  const string ver = ((qpos != end(values)) && !qpos->second.empty()) ? qpos->second : "1";
  // Accept the major version only, so a later "2.1" adding an optional token stays readable here.
  // Major version only, so a later "2.1" adding an optional token stays readable - but "10"/"27"
  //  are different majors, not v1/v2, so match the separator too.
  const bool known_ver = (ver == "1") || (ver == "2")
                         || SpecUtils::istarts_with(ver, "1.")
                         || SpecUtils::istarts_with(ver, "2.");
  if( !known_ver )
    throw runtime_error( "DetectionLimitSimple: invalid URI version" );

  // Collected while decoding, shown once at the end - a retired option is never applied silently.
  bool migrated_deprecated_norm = false;
  
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
    if( (stringstream(qpos->second) >> nsigma) )
    {
      ConfidenceLevel cl = ConfidenceLevel::NinetyFivePercent;
      
      switch( nsigma )
      {
        case 1:  cl = ConfidenceLevel::OneSigma;          break;
        case 2:  cl = ConfidenceLevel::TwoSigma;          break;
        case 3:  cl = ConfidenceLevel::ThreeSigma;        break;
        case 4:  cl = ConfidenceLevel::FourSigma;         break;
        case 5:  cl = ConfidenceLevel::FiveSigma;         break;
        case 95: cl = ConfidenceLevel::NinetyFivePercent; break;
        case 99: cl = ConfidenceLevel::NinetyNinePercent; break;
        default:
          cerr << "Invalid URI CL, '" << qpos->second << "' to a valid CL." << endl;
          break;
      }//switch( nsigma )
      
      m_confidenceLevel->setCurrentIndex( static_cast<int>(cl) );
    }else
    {
      cerr << "Failed to convert URI CL, '" << qpos->second << "' to a to an integer." << endl;
    }
  }//if( URI has CL )
  
  if( currieMethod )
  {
    bool noSignal = false;
    qpos = values.find( "ISBACK" );
    if( qpos != end(values) )
    {
      if( (qpos->second == "0") || SpecUtils::iequals_ascii(qpos->second, "NO") || SpecUtils::iequals_ascii(qpos->second, "FALSE") )
        noSignal = false;
      else if( (qpos->second == "1") || SpecUtils::iequals_ascii(qpos->second, "YES") || SpecUtils::iequals_ascii(qpos->second, "TRUE") )
        noSignal = true;
      else
        cerr << "Unexpected 'ISBACK' value: '" << qpos->second << "' (non-bool)" << endl;
    }//if( URI had is background value )
    
    m_isBackgroundSpectrum->setChecked( noSignal );
  }//if( currieMethod )

  // NSIDE is emitted for BOTH methods (the deconvolution "from sides" treatment uses it too), so it
  //  has to be decoded for both - it used to be read only under Currie, which silently reset the
  //  side-channel count on every reload of a deconvolution state.  No assert on the background
  //  checkbox here: a legacy VER=1 Currie URI could legitimately carry both, and an assert that
  //  fires in a Debug build kills the whole session.
  qpos = values.find( "NSIDE" );
  if( qpos != end(values) )
  {
    int nside;
    if( (stringstream(qpos->second) >> nside) && (nside >= 1) && (nside <= 64) )
      m_numSideChannel->setValue( nside );
    else
      cerr << "Invalid 'NSIDE' value: '" << qpos->second << "'" << endl;
  }//if( URI has NSIDE )
  
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
  
  // Which measurement the limit describes.  Emitted for both methods since VER=2; a VER=1 URI says
  //  the same thing with ISBACK (Currie) or CONTNORM=NOSIG (deconvolution).
  qpos = values.find( "MODEL" );
  if( qpos != end(values) )
  {
    if( SpecUtils::iequals_ascii(qpos->second, "BACKREF") )
      m_isBackgroundSpectrum->setChecked( true );
    else if( SpecUtils::iequals_ascii(qpos->second, "CUR") )
      m_isBackgroundSpectrum->setChecked( false );
    else
      cerr << "Invalid 'MODEL' value: '" << qpos->second << "'" << endl;
  }

  if( !currieMethod )
  {
    // Decodes both vocabularies - the v1 UNKNOWN/NOSIG/FIXED and the v2 FLOAT/EDGES - and migrates
    //  the retired "no signal" option onto a background-reference measurement.  Note this sets the
    //  ENUM, then converts; the old code decoded straight into a combo index whose order was the
    //  reverse of the enum for two of the three values.
    qpos = values.find( "CONTNORM" );
    if( qpos != end(values) )
    {
      DetectionLimitCalc::DeconContinuumNorm norm = DetectionLimitCalc::DeconContinuumNorm::Floating;
      DetectionLimitCalc::DeconMeasurementModel model
            = DetectionLimitCalc::DeconMeasurementModel::CurrentSpectrum;
      bool migrated = false;

      if( DetectionLimitCalc::decode_continuum_norm_token( qpos->second, norm, model, migrated ) )
      {
        m_continuumPrior->setCurrentIndex( DetectionLimitCalc::index_from_continuum_norm(norm) );
        if( migrated )
        {
          // The retired option asserted the spectrum is signal-free; that assertion is now the
          //  background checkbox, and it must be visibly reinstated rather than dropped.
          m_isBackgroundSpectrum->setChecked(
                model == DetectionLimitCalc::DeconMeasurementModel::BackgroundReference );
          migrated_deprecated_norm = true;
        }
      }else
      {
        cerr << "Invalid 'CONTNORM' value: '" << qpos->second << "'" << endl;
      }
    }//if( continuum prior provided )

    
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

  // Restore "Scale to dwell" state.  Encoded as an absolute target dwell (e.g. "60s",
  //  "5 min") so the targeted dwell is preserved regardless of which foreground is
  //  loaded at decode time.  Presence implies checkbox on; any failure leaves it off.
  qpos = values.find( "SCALE" );
  if( qpos != end(values) && m_planTimeCb )
  {
    bool ok = false;
    try
    {
      const double t = PhysicalUnits::stringToTimeDuration( qpos->second );
      if( t > 0.0 )
      {
        m_planTimeCb->setChecked( true );
        m_planTimeEdit->setEnabled( true );
        m_planTimeEdit->setText( WString::fromUTF8(
            PhysicalUnits::printToBestTimeUnits( t, 4 ) ) );
        ok = true;
      }
    }catch( std::exception & ){ }

    if( !ok )
    {
      cerr << "Invalid 'SCALE' value: '" << qpos->second << "'" << endl;
      m_planTimeCb->setChecked( false );
      m_planTimeEdit->setEnabled( false );
      // handleSpectrumChanged()/init() will refresh the disabled text on the next render.
    }
  }else if( m_planTimeCb )
  {
    // No SCALE in URL - ensure the checkbox is off (and the disabled text gets repopulated).
    m_planTimeCb->setChecked( false );
    m_planTimeEdit->setEnabled( false );
    const shared_ptr<const SpecUtils::Measurement> hist
                          = m_viewer->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( hist && (hist->real_time() > 0.0f) )
      m_planTimeEdit->setText( WString::fromUTF8(
          PhysicalUnits::printToBestTimeUnits( hist->real_time(), 3 ) ) );
  }
  
  // The advanced statistical inputs.  Reset all five to their defaults first, because for these the
  //  *absence* of a token is meaningful: no ALPHA means "still tracking the confidence level", and a
  //  VER=2 or VER=1 URI has none of them at all and must decode to Advanced-off.
  m_advancedCb->setChecked( false );
  m_alphaUserSet = m_betaUserSet = false;
  m_alpha->setValue( static_cast<float>( 1.0 - currentConfidenceLevel() ) );
  m_beta->setValue( static_cast<float>( 1.0 - currentConfidenceLevel() ) );
  m_distanceUncert->setText( "" );
  m_prevDistanceUncert = "";
  m_effUncert->setText( "" );

  qpos = values.find( "ADV" );
  if( qpos != end(values) )
  {
    if( (qpos->second == "1") || SpecUtils::iequals_ascii(qpos->second, "YES")
       || SpecUtils::iequals_ascii(qpos->second, "TRUE") )
      m_advancedCb->setChecked( true );
    else if( (qpos->second != "0") && !SpecUtils::iequals_ascii(qpos->second, "NO")
            && !SpecUtils::iequals_ascii(qpos->second, "FALSE") )
      cerr << "Invalid 'ADV' value: '" << qpos->second << "'" << endl;
  }//if( URI has ADV )

  qpos = values.find( "ALPHA" );
  if( qpos != end(values) )
  {
    float alpha;
    if( (stringstream(qpos->second) >> alpha) && (alpha > 0.0f) && (alpha < 0.5f) )
    {
      m_alpha->setValue( alpha );
      m_alphaUserSet = true;
    }else
    {
      cerr << "Invalid 'ALPHA' value: '" << qpos->second << "'" << endl;
    }
  }//if( URI has ALPHA )

  qpos = values.find( "BETA" );
  if( qpos != end(values) )
  {
    float beta;
    if( (stringstream(qpos->second) >> beta) && (beta > 0.0f) && (beta < 0.5f) )
    {
      m_beta->setValue( beta );
      m_betaUserSet = true;
    }else
    {
      cerr << "Invalid 'BETA' value: '" << qpos->second << "'" << endl;
    }
  }//if( URI has BETA )

  qpos = values.find( "DISTUNCERT" );
  if( qpos != end(values) )
  {
    try
    {
      const double dist_uncert = PhysicalUnits::stringToDistance( qpos->second );
      if( dist_uncert < 0.0 )
        throw runtime_error( "negative" );
      m_distanceUncert->setValueText( WString::fromUTF8(qpos->second) );
      m_prevDistanceUncert = m_distanceUncert->text();
    }catch( std::exception & )
    {
      cerr << "Invalid 'DISTUNCERT' value: '" << qpos->second << "'" << endl;
    }
  }//if( URI has DISTUNCERT )

  qpos = values.find( "EFFUNCERT" );
  if( qpos != end(values) )
  {
    double eff_uncert;
    if( (stringstream(qpos->second) >> eff_uncert) && (eff_uncert >= 0.0) && (eff_uncert < 100.0) )
      m_effUncert->setValueText( WString::fromUTF8(qpos->second) );
    else
      cerr << "Invalid 'EFFUNCERT' value: '" << qpos->second << "'" << endl;
  }//if( URI has EFFUNCERT )

  // Set render flags... JIC
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateLimit;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateDisplayedSpectrum;
  m_renderFlags |= DetectionLimitSimple::RenderActions::UpdateSpectrumDecorations;

  // The decode above only set widget values.  Drive the one function that owns the show/hide table
  //  so the dependent visibility matches, rather than hand-setting it here - two copies of those
  //  rules is exactly how they drifted before.
  //
  // MUST come before the `undoWasSet` clear below: `handleMethodChanged()` unconditionally re-sets
  //  AddUndoRedoStep, so calling it afterwards silently un-does the clear and makes every URL load
  //  - including the one an undo itself performs - insert a fresh undo step.
  handleMethodChanged();

  if( !undoWasSet )
    m_renderFlags.clear(RenderActions::AddUndoRedoStep);

  // A retired option was reinterpreted.  Loading it as something else without saying so would
  //  change what a saved number means without telling anyone, so this is never silent.  Safe inside
  //  the `BlockUndoRedoInserts` scope: `passMessage` goes to `InterSpecApp::svlog`, not to the
  //  undo/redo manager.
  if( migrated_deprecated_norm )
    passMessage( WString::tr("dl-migrate-nosig-notice"), WarningWidget::WarningMsgMedium );
}//void handleAppUrl( std::string uri )


std::string DetectionLimitSimple::encodeStateToUrl() const
{
  string answer;
  
  const bool currieMethod = (m_methodGroup->checkedId() == static_cast<int>(MethodIds::Currie));
  answer += currieMethod ? "CURRIE" : "DECON";
  
  // 2.1, not 3: the advanced tokens are purely additive and optional, and the version gate in
  //  `handleAppUrl()` already accepts a "2."-prefixed minor - it was written so exactly this could
  //  be added without breaking QR codes and stored states held by already-shipped builds.  Those
  //  builds decode a 2.1 URI losing only the advanced fields, which is the correct degradation.
  answer += "?VER=2.1";
  
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
    case ConfidenceLevel::NinetyFivePercent:  answer += "95"; break;
    case ConfidenceLevel::NinetyNinePercent:  answer += "99"; break;
    case ConfidenceLevel::OneSigma:           answer += "1";  break;
    case ConfidenceLevel::TwoSigma:           answer += "2";  break;
    case ConfidenceLevel::ThreeSigma: 	      answer += "3";  break;
    case ConfidenceLevel::FourSigma:          answer += "4";  break;
    case ConfidenceLevel::FiveSigma:          answer += "5";  break;
    case ConfidenceLevel::NumConfidenceLevel: assert(0);      break;
  }//switch( confidence )
  
  const bool useSideChan
        = (DetectionLimitCalc::continuum_norm_from_index( m_continuumPrior->currentIndex() )
           == DetectionLimitCalc::DeconContinuumNorm::FixedByEdges);
  if( (currieMethod && !m_isBackgroundSpectrum->isChecked()) || useSideChan )
    answer += "&NSIDE=" + std::to_string(m_numSideChannel->value());
  
  shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
  if( drf && (!drf->isValid() || !drf->hasResolutionInfo()) )
    drf.reset();
  
  // Only include FWHM if its set by the user, or there is no DRF
  //  I'm not totally decided how to handle this quantity
  if( !drf || (fabs(m_fwhm->value() - drf->peakResolutionFWHM(energy)) > 0.1) )
    answer += "&FWHM=" + m_fwhm->valueText().toUTF8();
  
  // Which measurement the limit describes.  One token for both methods since VER=2 - the checkbox
  //  means the same assertion either way, so ISBACK and CONTNORM=NOSIG would have been two spellings
  //  of one setting, free to contradict each other.  ISBACK is decode-only now.
  answer += "&MODEL=" + string( m_isBackgroundSpectrum->isChecked() ? "BACKREF" : "CUR" );

  if( !currieMethod )
  {
    // Never re-emits the retired value; a state carrying it was migrated on decode.
    answer += "&CONTNORM=" + DetectionLimitCalc::continuum_norm_token(
          DetectionLimitCalc::continuum_norm_from_index( m_continuumPrior->currentIndex() ) );

    if( !useSideChan )
    {
      answer += "&CONTTYPE=";
      switch( m_continuumType->currentIndex() )
      {
        case 0: answer += "LIN"; break;
        case 1: answer += "QUAD"; break;
        default:
          assert( 0 );
          throw std::logic_error( "Invalid continuum type selected" );
      }//switch( continuumTypeIndex )
    }//if( m_continuumPrior->currentIndex() != 2 )
  }//if( currieMethod ) / else
  
  if( !m_allGammasInRoi )
    answer += "&ALLGAMMA=0";

  // Encode scale as the absolute target dwell (in seconds), so reloading the URL against
  //  a different foreground preserves "the dwell I asked about" rather than rescaling it.
  //  Use printCompact (no units) for the seconds-value to keep the URL field free of spaces.
  // Hidden means the calculation ignores it, so it must not be persisted either.
  if( m_planTimeCb && m_planTimeCb->isChecked() && !m_planTimeCb->isHidden() )
  {
    try
    {
      const double t = PhysicalUnits::stringToTimeDuration( m_planTimeEdit->text().toUTF8() );
      if( t > 0.0 )
        answer += "&SCALE=" + SpecUtils::printCompact( t / PhysicalUnits::second, 6 ) + "s";
    }catch( std::exception & )
    {
      // Invalid time string - skip; result is reflected by SCALE absent from URL.
    }
  }

  // The advanced statistical inputs.  Each value is emitted independently of ADV, so unticking
  //  "Advanced" and then undoing gets the typed values back; ADV alone records whether they apply.
  //  The *absence* of ALPHA/BETA is how "still tracking the confidence level" is encoded - that is
  //  what makes a decoded state respond to the confidence-level combo the way a live one does.
  if( m_advancedCb->isChecked() )
    answer += "&ADV=1";

  if( m_alphaUserSet )
    answer += "&ALPHA=" + SpecUtils::printCompact( m_alpha->value(), 6 );

  if( m_betaUserSet )
    answer += "&BETA=" + SpecUtils::printCompact( m_beta->value(), 6 );

  if( m_distanceUncert->validate() == WValidator::State::Valid )
  {
    string dist_uncert = m_distanceUncert->text().toUTF8();
    SpecUtils::trim( dist_uncert );
    // Like DIST, distance is not localized, so the user's text can go straight into the URL.
    if( !dist_uncert.empty() )
      answer += "&DISTUNCERT=" + dist_uncert;
  }

  if( m_effUncert->validate() == WValidator::State::Valid )
  {
    string eff_uncert = m_effUncert->text().toUTF8();
    SpecUtils::trim( eff_uncert );
    if( !eff_uncert.empty() )
      answer += "&EFFUNCERT=" + eff_uncert;
  }

  return answer;
}//std::string encodeStateToUrl() const;

