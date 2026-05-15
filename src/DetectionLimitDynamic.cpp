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

#include <map>
#include <string>
#include <sstream>
#include <algorithm>

#include <Wt/Utils.h>
#include <Wt/WText.h>
#include <Wt/WTable.h>
#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>
#include <Wt/WLabel.h>
#include <Wt/WSpinBox.h>
#include <Wt/WCheckBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WComboBox.h>
#include <Wt/WTabWidget.h>
#include <Wt/WTableView.h>
#include <Wt/WPushButton.h>
#include <Wt/WGridLayout.h>
#include <Wt/WButtonGroup.h>
#include <Wt/WRadioButton.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WSuggestionPopup.h>
#include <Wt/WStandardItemModel.h>
#include <Wt/WStandardItem.h>

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/SimpleDialog.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/DetectionLimitTool.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/NuclideSourceEnter.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/DetectionLimitDynamic.h"
#include "InterSpec/DetectorPeakResponse.h"


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
  }

  /** Number of confidence levels in `ConfidenceLevel` below. */
  enum ConfidenceLevel
  {
    NinetyFivePercent = 0, // 0.95
    NinetyNinePercent,     // 0.99
    OneSigma,              // 0.682689492137086
    TwoSigma,              // 0.954499736103642
    ThreeSigma,            // 0.997300203936740
    NumConfidenceLevel
  };

  double cl_to_double( const int idx )
  {
    switch( idx )
    {
      case NinetyFivePercent:  return 0.95;
      case NinetyNinePercent:  return 0.99;
      case OneSigma:           return 0.682689492137086;
      case TwoSigma:           return 0.954499736103642;
      case ThreeSigma:         return 0.997300203936740;
      default: break;
    }
    return 0.95;
  }

  // Speed regex: a positive decimal optionally followed by one of
  // {m/s, cm/s, km/h, km/hr, mph, mi/h, ft/s}.  Whitespace-tolerant.
  const char *const sm_speedUnitOptionalRegex
      = "^\\s*\\+?\\s*\\d+(\\.\\d*)?([eE][+\\-]?\\d+)?\\s*"
        "(m\\/s|cm\\/s|km\\/h|km\\/hr|mph|mi\\/h|ft\\/s)?\\s*$";

  /** Parses a user-provided speed string and returns the speed in m/s.
   Throws std::runtime_error if it cannot be parsed. */
  double stringToSpeed_m_per_s( const std::string &input )
  {
    std::string s = input;
    // Trim
    size_t start = 0;
    while( start < s.size() && std::isspace(static_cast<unsigned char>(s[start])) ) ++start;
    size_t end = s.size();
    while( end > start && std::isspace(static_cast<unsigned char>(s[end-1])) ) --end;
    s = s.substr( start, end - start );
    if( s.empty() )
      throw std::runtime_error( "empty speed string" );

    // Find first non-numeric / non-{., e, E, +, -} character to split number/unit.
    size_t i = 0;
    while( i < s.size() )
    {
      const char c = s[i];
      const bool numChar = std::isdigit(static_cast<unsigned char>(c))
                        || (c=='.') || (c=='+') || (c=='-')
                        || (c=='e') || (c=='E');
      if( !numChar )
        break;
      ++i;
    }
    if( i == 0 )
      throw std::runtime_error( "no numeric prefix" );

    double value = 0.0;
    try
    {
      value = std::stod( s.substr(0, i) );
    }catch( const std::exception &e )
    {
      throw std::runtime_error( std::string("could not parse numeric prefix: ") + e.what() );
    }

    std::string unit = s.substr(i);
    // Trim
    while( !unit.empty() && std::isspace(static_cast<unsigned char>(unit.front())) )
      unit.erase( unit.begin() );
    while( !unit.empty() && std::isspace(static_cast<unsigned char>(unit.back())) )
      unit.pop_back();
    SpecUtils::to_lower_ascii( unit );

    if( unit.empty() || unit == "m/s" )
      return value;
    if( unit == "cm/s" )
      return value * 0.01;
    if( unit == "km/h" || unit == "km/hr" )
      return value * (1000.0 / 3600.0);
    if( unit == "mph" || unit == "mi/h" )
      return value * 0.44704;
    if( unit == "ft/s" )
      return value * 0.3048;

    throw std::runtime_error( "unknown speed unit '" + unit + "'" );
  }//stringToSpeed_m_per_s
}//namespace


// ===========================================================================
// DetectionLimitDynamicWindow
// ===========================================================================
DetectionLimitDynamicWindow::DetectionLimitDynamicWindow( InterSpec *viewer )
 : AuxWindow( WString::tr("window-title-dynamic-mda"),
             (AuxWindowProperties::TabletNotFullScreen
              | AuxWindowProperties::SetCloseable
              | AuxWindowProperties::DisableCollapse) ),
   m_tool( nullptr )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;

  rejectWhenEscapePressed( true );

  m_tool = contents()->addNew<DetectionLimitDynamic>( viewer );
  m_tool->setHeight( WLength(100, WLength::Unit::Percentage) );

  AuxWindow::addHelpInFooter( footer(), "dynamic-mda-dialog" );

  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );

  show();

  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  int width = 625, height = 800;
  if( (screenW > 100) && (screenW < width) )
    width = screenW;
  if( (screenH > 100) && (screenH < height) )
    height = screenH;

  setWidth( width );
  setMaximumSize( WLength::Auto, height );

  resizeToFitOnScreen();
  centerWindowHeavyHanded();
}//DetectionLimitDynamicWindow()


DetectionLimitDynamicWindow::~DetectionLimitDynamicWindow()
{
}


DetectionLimitDynamic *DetectionLimitDynamicWindow::tool()
{
  return m_tool;
}


// ===========================================================================
// DetectionLimitDynamic
// ===========================================================================
DetectionLimitDynamic::DetectionLimitDynamic( InterSpec *viewer )
 : WContainerWidget(),
   m_viewer( viewer ),
   m_spectrum( nullptr ),
   m_peakModel( nullptr ),
   m_tabs( nullptr ),
   m_nuclideEdit( nullptr ),
   m_nuclideAgeEdit( nullptr ),
   m_nucEnterController( nullptr ),
   m_photoPeakEnergy( nullptr ),
   m_photoPeakEnergiesAndBr{},
   m_confidenceLevel( nullptr ),
   m_fpPerHour( nullptr ),
   m_detectorDisplay( nullptr ),
   m_lowerRoi( nullptr ),
   m_upperRoi( nullptr ),
   m_numSideChannel( nullptr ),
   m_shieldingSelect( nullptr ),
   m_solveForGroup( nullptr ),
   m_speedRow( nullptr ),
   m_speedLabel( nullptr ),
   m_speedInput( nullptr ),
   m_dcaRow( nullptr ),
   m_dcaLabel( nullptr ),
   m_dcaInput( nullptr ),
   m_activityRow( nullptr ),
   m_activityLabel( nullptr ),
   m_activityInput( nullptr ),
   m_windowModeGroup( nullptr ),
   m_windowLengthManual( nullptr ),
   m_strideModeGroup( nullptr ),
   m_strideSeconds( nullptr ),
   m_strideSamples( nullptr ),
   m_strideDerivedTxt( nullptr ),
   m_snapSamplesCb( nullptr ),
   m_trialsAlphaTxt( nullptr ),
   m_resultTxt( nullptr ),
   m_calcErrTxt( nullptr ),
   m_moreInfoBtn( nullptr ),
   m_moreInfoWindow( nullptr ),
   m_searchInfoTxt( nullptr ),
   m_runSearchBtn( nullptr ),
   m_hitsView( nullptr ),
   m_hitsModel( nullptr ),
   m_searchErrTxt( nullptr ),
   m_currentInput( nullptr ),
   m_currentResult( nullptr ),
   m_currentSearch( nullptr ),
   m_stateUri()
{
  init();
}


DetectionLimitDynamic::~DetectionLimitDynamic()
{
}


void DetectionLimitDynamic::init()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  wApp->useStyleSheet( "InterSpec_resources/DetectionLimitDynamic.css" );
  m_viewer->useMessageResourceBundle( "DetectionLimitDynamic" );
  m_viewer->useMessageResourceBundle( "NuclideSourceEnter" );

  addStyleClass( "DetectionLimitDynamic" );

  // ----- chart at top: shows the background spectrum with ROI / side-channel
  // highlights, and a peak Gaussian when calc succeeds.
  {
    WContainerWidget *chartArea = addNew<WContainerWidget>();
    chartArea->addStyleClass( "ChartArea" );

    auto spectrumOwner = std::make_unique<D3SpectrumDisplayDiv>();
    m_spectrum = spectrumOwner.get();
    chartArea->addWidget( std::move(spectrumOwner) );

    m_spectrum->setXAxisTitle( "" );
    m_spectrum->setYAxisTitle( "", "" );
    m_spectrum->setYAxisLog( false );
    m_spectrum->disableLegend();
    m_spectrum->setShowPeakLabel( SpectrumChart::PeakLabels::kShowPeakUserLabel, true );

    m_peakModel = m_spectrum->addChild( std::make_unique<PeakModel>() );
    m_peakModel->setNoSpecMeasBacking();
    m_spectrum->setPeakModel( m_peakModel );

    // Let the user drag the ROI edges on the chart to set the ROI bounds.
    m_spectrum->existingRoiEdgeDragUpdate().connect( this,
      [this]( double new_lower, double new_upper, double /*px*/, double /*orig_lower*/,
              const std::string & /*spec_type*/, bool is_final_range )
      {
        if( !is_final_range )
          return;
        if( new_upper < new_lower )
          std::swap( new_lower, new_upper );

        // If we have a selected gamma, refuse a drag that would exclude it
        // (same behavior as Simple MDA).
        const float gamma_e = photopeakEnergy();
        if( (gamma_e > 0.0f) && ((gamma_e < new_lower) || (gamma_e > new_upper)) )
        {
          passMessage( WString::tr("dlsd-roi-changed-no-gamma"),
                       WarningWidget::WarningMsgHigh );
          return;
        }

        m_lowerRoi->setValue( static_cast<float>(new_lower) );
        m_upperRoi->setValue( static_cast<float>(new_upper) );
        m_renderFlags |= UpdateLimit;
  m_renderFlags |= AddUndoRedoStep;
        scheduleRender();
      } );
  }

  // ----- shared inputs (live on this widget, both tabs see them) -----
  WContainerWidget *generalInput = addNew<WContainerWidget>();
  generalInput->addStyleClass( "GeneralInput" );

  // Helpers: organize the general-input area into explicit lines, each holding
  // one or more label+input field cells.  `field_row` adds a cell to a line;
  // `line_row` starts a new horizontal line under `generalInput`.
  auto line_row = [generalInput]() -> WContainerWidget *
  {
    WContainerWidget *line = generalInput->addNew<WContainerWidget>();
    line->addStyleClass( "DynMdaInputLine" );
    return line;
  };

  auto field_row = []( WContainerWidget *parent ) -> WContainerWidget *
  {
    WContainerWidget *row = parent->addNew<WContainerWidget>();
    row->addStyleClass( "DynMdaFieldRow" );
    return row;
  };

  // Line 1: Nuclide, Age, Gamma
  WContainerWidget *line1 = line_row();
  {
    WContainerWidget *r = field_row( line1 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("nuclide-label") );
    m_nuclideEdit = r->addNew<WLineEdit>();
    m_nuclideEdit->setMinimumSize( 30, WLength::Auto );
    lbl->setBuddy( m_nuclideEdit );
  }

  {
    WContainerWidget *r = field_row( line1 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("age-label") );
    m_nuclideAgeEdit = r->addNew<WLineEdit>();
    m_nuclideAgeEdit->setMinimumSize( 30, WLength::Auto );
    m_nuclideAgeEdit->setPlaceholderText( WString::tr("N/A") );
    lbl->setBuddy( m_nuclideAgeEdit );
  }

  m_nucEnterController = addChild( std::make_unique<NuclideSourceEnterController>(
      m_nuclideEdit, m_nuclideAgeEdit, nullptr ) );
  m_nucEnterController->changed().connect( this, &DetectionLimitDynamic::handleNuclideChanged );

  {
    WContainerWidget *r = field_row( line1 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-gamma-label") );
    m_photoPeakEnergy = r->addNew<WComboBox>();
    m_photoPeakEnergy->activated().connect( this, &DetectionLimitDynamic::handleGammaChanged );
    lbl->setBuddy( m_photoPeakEnergy );
  }

  // Line 2: Confidence Level, Max FP/hour
  WContainerWidget *line2 = line_row();
  {
    WContainerWidget *r = field_row( line2 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-cl-label") );
    m_confidenceLevel = r->addNew<WComboBox>();
    m_confidenceLevel->addItem( "95%" );
    m_confidenceLevel->addItem( "99%" );
    m_confidenceLevel->addItem( WString::tr("dlsd-cl-1sigma") );
    m_confidenceLevel->addItem( WString::tr("dlsd-cl-2sigma") );
    m_confidenceLevel->addItem( WString::tr("dlsd-cl-3sigma") );
    m_confidenceLevel->setCurrentIndex( NinetyFivePercent );
    m_confidenceLevel->activated().connect( this, &DetectionLimitDynamic::handleConfidenceLevelChanged );
    lbl->setBuddy( m_confidenceLevel );
  }

  {
    WContainerWidget *r = field_row( line2 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-fp-per-hour-label") );
    m_fpPerHour = r->addNew<NativeFloatSpinBox>();
    m_fpPerHour->setRange( 1.0e-4f, 1000.0f );
    m_fpPerHour->setValue( 1.0f );
    m_fpPerHour->setSpinnerHidden();
    m_fpPerHour->valueChanged().connect( this, &DetectionLimitDynamic::handleFpPerHourChanged );
    lbl->setBuddy( m_fpPerHour );

    const bool showToolTips
        = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
    HelpSystem::attachToolTipOn( {lbl, m_fpPerHour},
                                 WString::tr("dlsd-tt-fp-per-hour"), showToolTips,
                                 HelpSystem::ToolTipPrefOverride::RespectPreference );
  }

  // Line 3: ROI Lower, ROI Upper, Side Channels
  WContainerWidget *line3 = line_row();
  {
    WContainerWidget *r = field_row( line3 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-roi-lower-label") );
    m_lowerRoi = r->addNew<NativeFloatSpinBox>();
    m_lowerRoi->setSpinnerHidden();
    m_lowerRoi->valueChanged().connect( this, [this](){
      m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
    } );
    lbl->setBuddy( m_lowerRoi );
  }

  {
    WContainerWidget *r = field_row( line3 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-roi-upper-label") );
    m_upperRoi = r->addNew<NativeFloatSpinBox>();
    m_upperRoi->setSpinnerHidden();
    m_upperRoi->valueChanged().connect( this, [this](){
      m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
    } );
    lbl->setBuddy( m_upperRoi );
  }

  {
    WContainerWidget *r = field_row( line3 );
    WLabel *lbl = r->addNew<WLabel>( WString::tr("dlsd-num-side-channel-label") );
    m_numSideChannel = r->addNew<WSpinBox>();
    m_numSideChannel->setRange( 1, 64 );
    m_numSideChannel->setValue( 4 );
    m_numSideChannel->valueChanged().connect( this, [this](){
      m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
    } );
    lbl->setBuddy( m_numSideChannel );
  }

  // Line 4: Detector display + Shielding (no label on shielding to save space)
  WContainerWidget *line4 = line_row();
  {
    WContainerWidget *r = field_row( line4 );
    SpectraFileModel *specFileModel = m_viewer->fileManager()->model();
    m_detectorDisplay = r->addNew<DetectorDisplay>( m_viewer, specFileModel );
  }
  m_viewer->detectorChanged().connect( this, [this]( std::shared_ptr<DetectorPeakResponse> det ){
    handleDetectorChanged( det );
  } );
  m_viewer->detectorModified().connect( this, [this]( std::shared_ptr<DetectorPeakResponse> det ){
    handleDetectorChanged( det );
  } );

  {
    m_viewer->initMaterialDbAndSuggestions();
    WSuggestionPopup * const materialSuggest = m_viewer->shieldingSuggester();
    WContainerWidget *r = field_row( line4 );
    r->addStyleClass( "DynMdaShieldingRow" );
    m_shieldingSelect = r->addNew<ShieldingSelect>( materialSuggest );
    m_shieldingSelect->materialEdit()->setPlaceholderText(
        WString("<{1}>").arg( WString::tr("dlsd-shielding-empty-text") ) );
    m_shieldingSelect->materialChanged().connect( this, &DetectionLimitDynamic::handleShieldingChanged );
    m_shieldingSelect->materialModified().connect( this, &DetectionLimitDynamic::handleShieldingChanged );
  }

  // Tabs
  m_tabs = addNew<WTabWidget>();
  m_tabs->addStyleClass( "DynamicMdaTabs" );

  {
    auto calcOwner = std::make_unique<WContainerWidget>();
    WContainerWidget *calc = calcOwner.get();
    buildCalculatorTab( calc );
    m_tabs->addTab( std::move(calcOwner), WString::tr("dlsd-tab-calculator") );
  }

  {
    auto searchOwner = std::make_unique<WContainerWidget>();
    WContainerWidget *search = searchOwner.get();
    buildSearchTab( search );
    m_tabs->addTab( std::move(searchOwner), WString::tr("dlsd-tab-search") );
  }

  m_viewer->displayedSpectrumChanged().connect( this, &DetectionLimitDynamic::handleSpectrumChanged );

  // Set initial enabled/disabled states for radio-driven inputs.
  handleSolveForChanged();
  handleWindowModeChanged();
  handleStrideModeChanged();

  // Pre-populate nuclide / age (and shielding) from the Reference Photopeaks
  // widget, if the user is currently displaying a nuclide there.
  if( ReferencePhotopeakDisplay * const reflines = m_viewer->referenceLinesWidget() )
  {
    const ReferenceLineInfo &current = reflines->currentlyShowingNuclide();
    if( current.m_nuclide )
    {
      m_nucEnterController->setNuclideText( current.m_nuclide->symbol );
      if( !current.m_input.m_age.empty() )
        m_nucEnterController->setNuclideAgeTxt( current.m_input.m_age );

      const ShieldingSelect * const refShield = reflines->shieldingSelect();
      if( refShield && m_shieldingSelect )
      {
        if( refShield->isGenericMaterial() )
        {
          if( refShield->arealDensity() > 0.0 )
            m_shieldingSelect->setAtomicNumberAndArealDensity( refShield->atomicNumber(),
                                                               refShield->arealDensity() );
        }else if( refShield->materialEdit() && refShield->thicknessEdit() )
        {
          const std::string mat = refShield->materialEdit()->text().toUTF8();
          const std::string thick = refShield->thicknessEdit()->text().toUTF8();
          if( !mat.empty() && !thick.empty() )
            m_shieldingSelect->setMaterialNameAndThickness( mat, thick );
        }
      }
    }
  }

  handleSpectrumChanged();
  m_renderFlags |= UpdateLimit;
  scheduleRender();
}//init()


void DetectionLimitDynamic::buildCalculatorTab( WContainerWidget *container )
{
  container->addStyleClass( "CalculatorTab" );

  // Solve-for radio group
  WContainerWidget *solveFor = container->addNew<WContainerWidget>();
  solveFor->addStyleClass( "SolveForRow" );
  solveFor->addNew<WLabel>( WString::tr("dlsd-solve-for-label") );
  m_solveForGroup = std::make_shared<WButtonGroup>();
  WRadioButton *rActivity = solveFor->addNew<WRadioButton>( WString::tr("dlsd-solve-activity") );
  m_solveForGroup->addButton( rActivity, static_cast<int>(SolveForId::Activity) );
  WRadioButton *rSpeed = solveFor->addNew<WRadioButton>( WString::tr("dlsd-solve-speed") );
  m_solveForGroup->addButton( rSpeed, static_cast<int>(SolveForId::Speed) );
  WRadioButton *rDca = solveFor->addNew<WRadioButton>( WString::tr("dlsd-solve-dca") );
  m_solveForGroup->addButton( rDca, static_cast<int>(SolveForId::Dca) );
  m_solveForGroup->setCheckedButton( rActivity );
  m_solveForGroup->checkedChanged().connect( this, &DetectionLimitDynamic::handleSolveForChanged );

  // Geometry inputs — unit-aware text fields with regex validators.
  // Each (label, input) pair lives in its own flex sub-row so we can hide
  // the row corresponding to "solve for" as a single unit.
  // Layout defined in InterSpec_resources/DetectionLimitDynamic.css.
  auto field_row_calc = []( WContainerWidget *parent ) -> WContainerWidget *
  {
    WContainerWidget *row = parent->addNew<WContainerWidget>();
    row->addStyleClass( "DynMdaFieldRow" );
    return row;
  };

  auto make_validator = []( const char *regex )
  {
    auto v = std::make_shared<WRegExpValidator>( regex );
    v->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );
    return v;
  };

  WContainerWidget *geom = container->addNew<WContainerWidget>();
  geom->addStyleClass( "GeometryRow" );

  // Speed row
  m_speedRow = field_row_calc( geom );
  m_speedLabel = m_speedRow->addNew<WLabel>( WString::tr("dlsd-speed-label") );
  m_speedInput = m_speedRow->addNew<WLineEdit>( "1 m/s" );
  m_speedInput->setValidator( make_validator( sm_speedUnitOptionalRegex ) );
  m_speedInput->setAttributeValue( "ondragstart", "return false" );
  m_speedLabel->setBuddy( m_speedInput );
  m_speedInput->changed().connect( this, &DetectionLimitDynamic::handleSpeedChanged );
  m_speedInput->enterPressed().connect( this, &DetectionLimitDynamic::handleSpeedChanged );

  // DCA row
  m_dcaRow = field_row_calc( geom );
  m_dcaLabel = m_dcaRow->addNew<WLabel>( WString::tr("dlsd-dca-label") );
  m_dcaInput = m_dcaRow->addNew<WLineEdit>( "1 m" );
  m_dcaInput->setValidator( make_validator( PhysicalUnits::sm_distanceUnitOptionalRegex ) );
  m_dcaInput->setAttributeValue( "ondragstart", "return false" );
  m_dcaLabel->setBuddy( m_dcaInput );
  m_dcaInput->changed().connect( this, &DetectionLimitDynamic::handleDcaChanged );
  m_dcaInput->enterPressed().connect( this, &DetectionLimitDynamic::handleDcaChanged );

  // Activity row
  m_activityRow = field_row_calc( geom );
  m_activityLabel = m_activityRow->addNew<WLabel>( WString::tr("dlsd-activity-label") );
  m_activityInput = m_activityRow->addNew<WLineEdit>( "" );
  m_activityInput->setPlaceholderText( WString::tr("dlsd-activity-placeholder") );
  m_activityInput->setValidator( make_validator( PhysicalUnits::sm_activityUnitOptionalRegex ) );
  m_activityInput->setAttributeValue( "ondragstart", "return false" );
  m_activityLabel->setBuddy( m_activityInput );
  m_activityInput->changed().connect( this, &DetectionLimitDynamic::handleActivityChanged );
  m_activityInput->enterPressed().connect( this, &DetectionLimitDynamic::handleActivityChanged );

  // Window controls
  WContainerWidget *wnd = container->addNew<WContainerWidget>();
  wnd->addStyleClass( "WindowRow" );
  wnd->addNew<WLabel>( WString::tr("dlsd-window-len-label") );
  m_windowModeGroup = std::make_shared<WButtonGroup>();
  WRadioButton *wAuto = wnd->addNew<WRadioButton>( WString::tr("dlsd-window-auto") );
  m_windowModeGroup->addButton( wAuto, static_cast<int>(WindowModeId::Auto) );
  WRadioButton *wMan = wnd->addNew<WRadioButton>( WString::tr("dlsd-window-manual") );
  m_windowModeGroup->addButton( wMan, static_cast<int>(WindowModeId::Manual) );
  m_windowModeGroup->setCheckedButton( wAuto );
  m_windowModeGroup->checkedChanged().connect( this, &DetectionLimitDynamic::handleWindowModeChanged );

  m_windowLengthManual = wnd->addNew<WLineEdit>( "1 s" );
  m_windowLengthManual->setValidator( make_validator( PhysicalUnits::sm_timeDurationRegex ) );
  m_windowLengthManual->setAttributeValue( "ondragstart", "return false" );
  m_windowLengthManual->changed().connect( this, &DetectionLimitDynamic::handleWindowLengthChanged );
  m_windowLengthManual->enterPressed().connect( this, &DetectionLimitDynamic::handleWindowLengthChanged );

  m_snapSamplesCb = wnd->addNew<WCheckBox>( WString::tr("dlsd-snap-samples") );
  m_snapSamplesCb->setChecked( true );
  m_snapSamplesCb->checked().connect( this, &DetectionLimitDynamic::handleSnapSamplesChanged );
  m_snapSamplesCb->unChecked().connect( this, &DetectionLimitDynamic::handleSnapSamplesChanged );

  // Stride
  WContainerWidget *strideC = container->addNew<WContainerWidget>();
  strideC->addStyleClass( "StrideRow" );
  strideC->addNew<WLabel>( WString::tr("dlsd-stride-label") );
  m_strideModeGroup = std::make_shared<WButtonGroup>();
  WRadioButton *sSec = strideC->addNew<WRadioButton>( WString::tr("dlsd-stride-seconds") );
  m_strideModeGroup->addButton( sSec, static_cast<int>(StrideModeId::Seconds) );
  WRadioButton *sSamp = strideC->addNew<WRadioButton>( WString::tr("dlsd-stride-samples") );
  m_strideModeGroup->addButton( sSamp, static_cast<int>(StrideModeId::Samples) );
  m_strideModeGroup->setCheckedButton( sSec );
  m_strideModeGroup->checkedChanged().connect( this, &DetectionLimitDynamic::handleStrideModeChanged );

  m_strideSeconds = strideC->addNew<WLineEdit>( "1 s" );
  m_strideSeconds->setValidator( make_validator( PhysicalUnits::sm_timeDurationRegex ) );
  m_strideSeconds->setAttributeValue( "ondragstart", "return false" );
  m_strideSeconds->changed().connect( this, &DetectionLimitDynamic::handleStrideChanged );
  m_strideSeconds->enterPressed().connect( this, &DetectionLimitDynamic::handleStrideChanged );

  m_strideSamples = strideC->addNew<WSpinBox>();
  m_strideSamples->setRange( 1, 1000 );
  m_strideSamples->setValue( 1 );
  m_strideSamples->valueChanged().connect( this, [this](){
    handleStrideChanged();
  } );

  m_strideDerivedTxt = strideC->addNew<WText>( "" );
  m_strideDerivedTxt->addStyleClass( "StrideDerived" );

  m_trialsAlphaTxt = container->addNew<WText>( "" );
  m_trialsAlphaTxt->addStyleClass( "TrialsAlphaTxt" );
  m_trialsAlphaTxt->setInline( false );

  // Result line
  m_calcErrTxt = container->addNew<WText>( "" );
  m_calcErrTxt->addStyleClass( "ErrMsg" );
  m_calcErrTxt->setInline( false );

  m_resultTxt = container->addNew<WText>( "&nbsp;" );
  m_resultTxt->addStyleClass( "ResultsTxtArea" );
  m_resultTxt->setInline( false );

  m_moreInfoBtn = container->addNew<WPushButton>( WString::tr("dlsd-more-info") );
  m_moreInfoBtn->setStyleClass( "LinkBtn" );
  m_moreInfoBtn->clicked().connect( this, &DetectionLimitDynamic::showMoreInfoDialog );
  m_moreInfoBtn->setHiddenKeepsGeometry( true );
  m_moreInfoBtn->hide();
}//buildCalculatorTab(...)


void DetectionLimitDynamic::buildSearchTab( WContainerWidget *container )
{
  container->addStyleClass( "SearchTab" );

  m_searchInfoTxt = container->addNew<WText>( "" );
  m_searchInfoTxt->setInline( false );
  m_searchInfoTxt->addStyleClass( "SearchInfo" );

  m_runSearchBtn = container->addNew<WPushButton>( WString::tr("dlsd-run-search") );
  m_runSearchBtn->clicked().connect( this, &DetectionLimitDynamic::runSearch );

  m_searchErrTxt = container->addNew<WText>( "" );
  m_searchErrTxt->addStyleClass( "ErrMsg" );
  m_searchErrTxt->setInline( false );

  m_hitsModel = std::make_shared<WStandardItemModel>();
  m_hitsModel->insertColumns( 0, 8 );
  m_hitsModel->setHeaderData( 0, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-first-sample"), Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 1, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-last-sample"),  Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 2, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-time"),         Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 3, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-real-time"),    Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 4, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-counts"),       Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 5, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-bg"),           Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 6, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-lc"),           Wt::ItemDataRole::Display );
  m_hitsModel->setHeaderData( 7, Wt::Orientation::Horizontal, WString::tr("dlsd-hit-col-sigma"),        Wt::ItemDataRole::Display );

  m_hitsView = container->addNew<WTableView>();
  m_hitsView->setModel( m_hitsModel );
  m_hitsView->setAlternatingRowColors( true );
  m_hitsView->setSelectionMode( SelectionMode::Single );
  m_hitsView->setColumnResizeEnabled( true );
  m_hitsView->setMinimumSize( WLength::Auto, 200 );
}//buildSearchTab(...)


void DetectionLimitDynamic::render( WFlags<RenderFlag> flags )
{
  WContainerWidget::render( flags );

  if( m_renderFlags.test(AddUndoRedoStep) )
  {
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo )
    {
      std::string uri = encodeStateToUrl();
      const bool sameAsPrev = (uri == m_stateUri);
      if( !m_stateUri.empty() && undoRedo->canAddUndoRedoNow() && !sameAsPrev )
      {
        const std::shared_ptr<const std::string> prev = std::make_shared<std::string>( std::move(m_stateUri) );
        const std::shared_ptr<const std::string> current = std::make_shared<std::string>( uri );
        auto undo_redo = [prev, current]( bool is_undo ){
          DetectionLimitDynamicWindow *win = InterSpec::instance()->showDynamicMdaWindow();
          DetectionLimitDynamic *tool = win ? win->tool() : nullptr;
          const std::string &target = is_undo ? *prev : *current;
          if( tool && !target.empty() )
            tool->handleAppUrl( target );
        };
        auto undo = [undo_redo](){ undo_redo(true); };
        auto redo = [undo_redo](){ undo_redo(false); };
        undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Update Dynamic MDA values" );
      }
      m_stateUri = std::move(uri);
    }
  }

  if( m_renderFlags.test(UpdateLimit) )
    updateCalcResult();

  m_renderFlags = WFlags<RenderActions>{};

  if( m_stateUri.empty() )
    m_stateUri = encodeStateToUrl();
}//render(...)


// ----- input handlers ------------------------------------------------------
void DetectionLimitDynamic::handleSpectrumChanged()
{
  // Invalidate the sample-period cache; foreground may have changed.
  m_samplePeriodCacheKey = nullptr;
  m_samplePeriodCache = 0.0;

  refreshGammasForCurrentNuclide();

  // Push the current background (or foreground if no background) into the chart.
  if( m_spectrum )
  {
    std::shared_ptr<const SpecUtils::Measurement> hist
        = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    if( !hist )
      hist = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    m_spectrum->setData( hist, true );
  }

  // Update search-tab info text
  const std::shared_ptr<const SpecUtils::SpecFile> foreground = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  WString info;
  if( foreground )
  {
    const bool is_pass = foreground->passthrough();
    const double period = currentSamplePeriod();
    char periodBuf[32];
    snprintf( periodBuf, sizeof(periodBuf), "%.3f", period );
    info = WString::tr("dlsd-search-info-template")
             .arg( WString::fromUTF8(foreground->filename()) )
             .arg( WString::tr( is_pass ? "dlsd-yes" : "dlsd-no" ) )
             .arg( WString::fromUTF8(periodBuf) );
  }
  else
  {
    info = WString::tr("dlsd-no-foreground");
  }
  if( m_searchInfoTxt )
    m_searchInfoTxt->setText( info );

  m_renderFlags |= UpdateLimit;
  scheduleRender();
}//handleSpectrumChanged()


void DetectionLimitDynamic::handleNuclideChanged()
{
  refreshGammasForCurrentNuclide();
  m_renderFlags |= UpdateLimit;
  m_renderFlags |= AddUndoRedoStep;
  scheduleRender();
}


void DetectionLimitDynamic::handleGammaChanged()
{
  // When the gamma changes, reset ROI to ~2.5 FWHM if DRF available;
  // otherwise reasonable default.
  const int idx = m_photoPeakEnergy->currentIndex();
  if( (idx >= 0) && (idx < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
  {
    const double energy = m_photoPeakEnergiesAndBr[idx].first;
    const std::shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
    double fwhm = 0.05 * energy;
    if( drf && drf->hasResolutionInfo() )
      fwhm = std::max( 0.1, static_cast<double>(drf->peakResolutionFWHM( energy )) );
    m_lowerRoi->setValue( energy - 1.25 * fwhm );
    m_upperRoi->setValue( energy + 1.25 * fwhm );
  }
  m_renderFlags |= UpdateLimit;
  m_renderFlags |= AddUndoRedoStep;
  scheduleRender();
}


void DetectionLimitDynamic::handleConfidenceLevelChanged()  { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleFpPerHourChanged()        { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleSolveForChanged()
{
  // Hide the row corresponding to the quantity we're solving for, so it's
  // clear to the user which fields they need to fill in.
  const SolveForId id = static_cast<SolveForId>( m_solveForGroup->checkedId() );

  if( m_speedRow )
    m_speedRow->setHidden( id == SolveForId::Speed );
  if( m_dcaRow )
    m_dcaRow->setHidden( id == SolveForId::Dca );
  if( m_activityRow )
    m_activityRow->setHidden( id == SolveForId::Activity );

  m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
}
void DetectionLimitDynamic::handleSpeedChanged()    { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleDcaChanged()      { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleActivityChanged() { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleWindowModeChanged()
{
  // Hide the manual input field when "Auto" is selected so it's clear the
  // user does not need (and cannot) provide a window length.
  const WindowModeId id = static_cast<WindowModeId>( m_windowModeGroup->checkedId() );
  m_windowLengthManual->setHidden( id == WindowModeId::Auto );
  m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
}
void DetectionLimitDynamic::handleWindowLengthChanged() { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleStrideModeChanged()
{
  // Hide the inactive stride input so it's clear which one to fill in.
  const StrideModeId id = static_cast<StrideModeId>( m_strideModeGroup->checkedId() );
  m_strideSeconds->setHidden( id == StrideModeId::Samples );
  m_strideSamples->setHidden( id == StrideModeId::Seconds );
  m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender();
}
void DetectionLimitDynamic::handleStrideChanged()       { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }
void DetectionLimitDynamic::handleSnapSamplesChanged()  { m_renderFlags |= UpdateLimit; m_renderFlags |= AddUndoRedoStep; scheduleRender(); }

void DetectionLimitDynamic::handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> )
{
  m_renderFlags |= UpdateLimit;
  scheduleRender();
}

void DetectionLimitDynamic::handleShieldingChanged()
{
  m_renderFlags |= UpdateLimit;
  m_renderFlags |= AddUndoRedoStep;
  scheduleRender();
}


void DetectionLimitDynamic::refreshGammasForCurrentNuclide()
{
  m_photoPeakEnergiesAndBr.clear();
  if( m_photoPeakEnergy )
    m_photoPeakEnergy->clear();

  const SandiaDecay::Nuclide *nuc = m_nucEnterController->nuclide();
  if( !nuc )
    return;

  const double age = m_nucEnterController->nuclideAge();

  SandiaDecay::NuclideMixture mix;
  const double dummy = 0.001 * SandiaDecay::curie;
  mix.addAgedNuclideByActivity( nuc, dummy, age );

  std::vector<SandiaDecay::EnergyRatePair> photons
      = mix.xrays( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
  const std::vector<SandiaDecay::EnergyRatePair> gammas
      = mix.gammas( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, true );
  photons.insert( end(photons), begin(gammas), end(gammas) );

  for( const SandiaDecay::EnergyRatePair &erp : photons )
  {
    const double br = erp.numPerSecond / dummy;
    if( br <= std::numeric_limits<double>::epsilon() )
      continue;
    if( erp.energy < 10.0 )
      continue;
    char buf[64];
    snprintf( buf, sizeof(buf), "%.2f keV (BR=%.3g)", erp.energy, br );
    m_photoPeakEnergiesAndBr.emplace_back( erp.energy, br );
    if( m_photoPeakEnergy )
      m_photoPeakEnergy->addItem( buf );
  }

  if( !m_photoPeakEnergiesAndBr.empty() && m_photoPeakEnergy )
  {
    // Pick the most intense gamma above 50 keV by default.
    size_t bestIdx = 0;
    double bestBr = 0.0;
    for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
    {
      if( (m_photoPeakEnergiesAndBr[i].first > 50.0)
          && (m_photoPeakEnergiesAndBr[i].second > bestBr) )
      {
        bestBr = m_photoPeakEnergiesAndBr[i].second;
        bestIdx = i;
      }
    }
    m_photoPeakEnergy->setCurrentIndex( static_cast<int>(bestIdx) );
    handleGammaChanged();
  }
}


float DetectionLimitDynamic::photopeakEnergy() const
{
  const int idx = m_photoPeakEnergy->currentIndex();
  if( (idx < 0) || (idx >= static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    return 0.0f;
  return static_cast<float>( m_photoPeakEnergiesAndBr[idx].first );
}

const SandiaDecay::Nuclide *DetectionLimitDynamic::nuclide() const
{
  return m_nucEnterController ? m_nucEnterController->nuclide() : nullptr;
}


double DetectionLimitDynamic::currentConfidenceLevel() const
{
  return cl_to_double( m_confidenceLevel->currentIndex() );
}


double DetectionLimitDynamic::currentSamplePeriod() const
{
  const std::shared_ptr<const SpecUtils::SpecFile> file = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( !file )
  {
    m_samplePeriodCacheKey = nullptr;
    m_samplePeriodCache = 0.0;
    return 0.0;
  }

  if( file.get() == m_samplePeriodCacheKey )
    return m_samplePeriodCache;

  const std::set<int> &samples = file->sample_numbers();
  double result = 0.0;
  if( samples.size() >= 2 )
  {
    // Iterate raw measurements: real_time is a scalar, no allocation needed.
    std::vector<double> rts;
    rts.reserve( samples.size() );
    const std::vector<std::shared_ptr<const SpecUtils::Measurement>> &meas = file->measurements();
    for( const std::shared_ptr<const SpecUtils::Measurement> &m : meas )
    {
      if( m && m->real_time() > 0.0f )
        rts.push_back( m->real_time() );
    }
    if( !rts.empty() )
    {
      std::sort( rts.begin(), rts.end() );
      result = rts[ rts.size()/2 ];
    }
  }

  m_samplePeriodCacheKey = file.get();
  m_samplePeriodCache = result;
  return result;
}


void DetectionLimitDynamic::updateTrialsDisplay()
{
  if( !m_currentResult || !m_trialsAlphaTxt )
    return;
  char nphBuf[32], alphaBuf[32], kaBuf[32], kbBuf[32];
  snprintf( nphBuf,   sizeof(nphBuf),   "%.1f", m_currentResult->trials_per_hour );
  snprintf( alphaBuf, sizeof(alphaBuf), "%.3g", m_currentResult->per_trial_alpha );
  snprintf( kaBuf,    sizeof(kaBuf),    "%.3f", m_currentResult->k_alpha );
  snprintf( kbBuf,    sizeof(kbBuf),    "%.3f", m_currentResult->k_beta );
  m_trialsAlphaTxt->setText( WString::tr("dlsd-trials-line")
      .arg( WString::fromUTF8(nphBuf) )
      .arg( WString::fromUTF8(alphaBuf) )
      .arg( WString::fromUTF8(kaBuf) )
      .arg( WString::fromUTF8(kbBuf) ) );
  m_trialsAlphaTxt->setTextFormat( TextFormat::XHTML );
}


DetectionLimitCalc::DynamicMdaInput DetectionLimitDynamic::buildCalcInput() const
{
  DetectionLimitCalc::DynamicMdaInput in;

  // Background spectrum: prefer the dedicated background; fallback foreground.
  in.background_spectrum = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
  if( !in.background_spectrum )
    in.background_spectrum = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );

  in.drf = m_detectorDisplay->detector();
  in.nuclide = nuclide();
  in.age = m_nucEnterController ? m_nucEnterController->nuclideAge() : 0.0;
  in.gamma_energy = photopeakEnergy();

  const int gIdx = m_photoPeakEnergy->currentIndex();
  if( (gIdx >= 0) && (gIdx < static_cast<int>(m_photoPeakEnergiesAndBr.size())) )
    in.branch_ratio = m_photoPeakEnergiesAndBr[gIdx].second;
  else
    in.branch_ratio = 1.0;

  in.roi_lower_energy = m_lowerRoi->value();
  in.roi_upper_energy = m_upperRoi->value();
  in.num_lower_side_channels = static_cast<size_t>( std::max(0, m_numSideChannel->value()) );
  in.num_upper_side_channels = in.num_lower_side_channels;

  in.detection_probability = currentConfidenceLevel();
  in.max_false_positives_per_hour = std::max( 1.0e-4, static_cast<double>(m_fpPerHour->value()) );
  in.additional_uncertainty = 0.0;

  // Source shielding (constant transmission multiplier at gamma_energy)
  if( m_shieldingSelect )
  {
    in.shielding_generic = m_shieldingSelect->isGenericMaterial();
    if( in.shielding_generic )
    {
      try
      {
        in.shielding_atomic_number = m_shieldingSelect->atomicNumber();
        in.shielding_areal_density = m_shieldingSelect->arealDensity();
      }catch( std::exception & )
      {
        in.shielding_atomic_number = 0.0;
        in.shielding_areal_density = 0.0;
      }
    }else
    {
      in.shielding_material = m_shieldingSelect->material();
      try
      {
        in.shielding_thickness = in.shielding_material ? m_shieldingSelect->thickness() : 0.0;
      }catch( std::exception & )
      {
        in.shielding_thickness = 0.0;
      }
    }
  }

  const SolveForId solveFor = static_cast<SolveForId>( m_solveForGroup->checkedId() );

  // Parse the unit-aware text fields.  The user can enter values like
  // "5 m/s", "100 cm", "1 uCi".  Missing units default to the SI base
  // (m/s, m, Bq).
  if( solveFor == SolveForId::Speed )
  {
    in.speed_m_per_s = 0.0;
  }else
  {
    try
    {
      in.speed_m_per_s = stringToSpeed_m_per_s( m_speedInput->text().toUTF8() );
    }catch( std::exception &e )
    {
      throw std::runtime_error( WString::tr("dlsd-err-invalid-speed").arg( WString::fromUTF8(e.what()) ).toUTF8() );
    }
  }

  if( solveFor == SolveForId::Dca )
  {
    in.dca_m = 0.0;
  }else
  {
    try
    {
      // stringToDistance returns InterSpec internal length units; divide
      // by PhysicalUnits::meter to get meters.
      const double d_internal = PhysicalUnits::stringToDistance( m_dcaInput->text().toUTF8() );
      in.dca_m = d_internal / PhysicalUnits::meter;
    }catch( std::exception &e )
    {
      throw std::runtime_error( WString::tr("dlsd-err-invalid-distance").arg( WString::fromUTF8(e.what()) ).toUTF8() );
    }
  }

  if( solveFor == SolveForId::Activity )
  {
    in.activity_bq = 0.0;
  }else
  {
    const std::string actTxt = m_activityInput->text().toUTF8();
    if( actTxt.empty() )
      throw std::runtime_error( WString::tr("dlsd-err-no-activity").toUTF8() );
    try
    {
      // stringToActivity returns in internal units (becquerel = 1).
      const double a_internal = PhysicalUnits::stringToActivity( actTxt );
      in.activity_bq = a_internal / PhysicalUnits::becquerel;
    }catch( std::exception &e )
    {
      throw std::runtime_error( WString::tr("dlsd-err-invalid-activity").arg( WString::fromUTF8(e.what()) ).toUTF8() );
    }
  }

  // Window length
  const WindowModeId wm = static_cast<WindowModeId>( m_windowModeGroup->checkedId() );
  if( wm == WindowModeId::Auto )
  {
    in.window_length_s = 0.0;
  }else
  {
    const std::string s = m_windowLengthManual->text().toUTF8();
    try
    {
      in.window_length_s = PhysicalUnits::stringToTimeDuration( s ) / PhysicalUnits::second;
    }catch( std::exception &e )
    {
      throw std::runtime_error( WString::tr("dlsd-err-invalid-window-length").arg( WString::fromUTF8(e.what()) ).toUTF8() );
    }
  }
  in.snap_window_to_samples = m_snapSamplesCb->isChecked();
  in.sample_period_s = currentSamplePeriod();

  // Stride
  const StrideModeId sm = static_cast<StrideModeId>( m_strideModeGroup->checkedId() );
  if( sm == StrideModeId::Samples )
  {
    const double period = (in.sample_period_s > 0.0) ? in.sample_period_s : 1.0;
    in.window_stride_s = std::max( 1, m_strideSamples->value() ) * period;
  }
  else
  {
    const std::string s = m_strideSeconds->text().toUTF8();
    try
    {
      const double t = PhysicalUnits::stringToTimeDuration( s ) / PhysicalUnits::second;
      in.window_stride_s = std::max( 1.0e-3, t );
    }catch( std::exception &e )
    {
      throw std::runtime_error( WString::tr("dlsd-err-invalid-stride").arg( WString::fromUTF8(e.what()) ).toUTF8() );
    }
  }

  return in;
}


void DetectionLimitDynamic::updateCalcResult()
{
  m_currentInput.reset();
  m_currentResult.reset();

  if( m_calcErrTxt )  m_calcErrTxt->setText( "" );
  if( m_resultTxt )   m_resultTxt->setText( "&nbsp;" );
  if( m_moreInfoBtn ) m_moreInfoBtn->hide();
  if( m_trialsAlphaTxt ) m_trialsAlphaTxt->setText( "" );
  if( m_strideDerivedTxt ) m_strideDerivedTxt->setText( "" );

  // Build a CurrieMdaInput so we can decorate the chart with ROI / side-channel
  // highlights regardless of whether the dynamic-MDA calc succeeds.
  auto build_curie_input = [this]() -> DetectionLimitCalc::CurrieMdaInput
  {
    DetectionLimitCalc::CurrieMdaInput cin;
    cin.spectrum = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Background );
    if( !cin.spectrum )
      cin.spectrum = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    cin.gamma_energy = photopeakEnergy();
    cin.roi_lower_energy = m_lowerRoi->value();
    cin.roi_upper_energy = m_upperRoi->value();
    cin.num_lower_side_channels = std::max(0, m_numSideChannel->value());
    cin.num_upper_side_channels = cin.num_lower_side_channels;
    cin.detection_probability = currentConfidenceLevel();
    cin.additional_uncertainty = 0.0f;
    return cin;
  };

  // Clear the chart decorations + peak.  We'll repopulate below if anything works.
  if( m_spectrum && m_peakModel )
  {
    m_peakModel->setPeaks( std::vector<std::shared_ptr<const PeakDef>>{} );
    m_spectrum->removeAllDecorativeHighlightRegions();
  }

  try
  {
    DetectionLimitCalc::DynamicMdaInput in = buildCalcInput();
    DetectionLimitCalc::DynamicMdaResult r = DetectionLimitCalc::compute_dynamic_mda( in );

    m_currentInput  = std::make_shared<DetectionLimitCalc::DynamicMdaInput>( std::move(in) );
    m_currentResult = std::make_shared<DetectionLimitCalc::DynamicMdaResult>( std::move(r) );

    const bool use_curie = use_curie_units();

    const SolveForId solveFor = static_cast<SolveForId>( m_solveForGroup->checkedId() );
    const std::string actStr = PhysicalUnits::printToBestActivityUnits( m_currentResult->activity_bq * SandiaDecay::becquerel, 3, use_curie );

    // Pre-format each numeric piece with the desired precision, then feed
    // the formatted strings into the localized template.
    char vBuf[32], dBuf[32], clBuf[32], fpBuf[32];
    snprintf( vBuf,  sizeof(vBuf),  "%.2f", m_currentResult->speed_m_per_s );
    snprintf( dBuf,  sizeof(dBuf),  "%.2f", m_currentResult->dca_m );
    snprintf( clBuf, sizeof(clBuf), "%.1f", 100.0 * m_currentInput->detection_probability );
    snprintf( fpBuf, sizeof(fpBuf), "%.3g", m_currentInput->max_false_positives_per_hour );

    WString resultText;
    switch( solveFor )
    {
      case SolveForId::Activity:
        resultText = WString::tr("dlsd-result-solve-activity")
            .arg( WString::fromUTF8(vBuf) )
            .arg( WString::fromUTF8(dBuf) )
            .arg( WString::fromUTF8(actStr) )
            .arg( WString::fromUTF8(clBuf) )
            .arg( WString::fromUTF8(fpBuf) );
        break;

      case SolveForId::Speed:
      {
        char sBuf[32];
        snprintf( sBuf, sizeof(sBuf), "%.3f", m_currentResult->speed_m_per_s );
        resultText = WString::tr("dlsd-result-solve-speed")
            .arg( WString::fromUTF8(actStr) )
            .arg( WString::fromUTF8(dBuf) )
            .arg( WString::fromUTF8(sBuf) )
            .arg( WString::fromUTF8(clBuf) );
        break;
      }

      case SolveForId::Dca:
      {
        char dBuf2[32];
        snprintf( dBuf2, sizeof(dBuf2), "%.3f", m_currentResult->dca_m );
        resultText = WString::tr("dlsd-result-solve-dca")
            .arg( WString::fromUTF8(actStr) )
            .arg( WString::fromUTF8(vBuf) )
            .arg( WString::fromUTF8(dBuf2) )
            .arg( WString::fromUTF8(clBuf) );
        break;
      }
    }

    if( m_resultTxt )
    {
      m_resultTxt->setText( resultText );
      m_resultTxt->setTextFormat( TextFormat::XHTML );
    }
    if( !m_currentResult->warning.empty() && m_calcErrTxt )
      m_calcErrTxt->setText( WString::fromUTF8( m_currentResult->warning ) );
    else if( m_calcErrTxt )
      m_calcErrTxt->setText( "" );

    // If the calc had to clamp the stride up to T_w, write the realized value
    // back into the widget so the user sees what was actually used.
    if( m_strideModeGroup && m_strideSeconds
        && (static_cast<StrideModeId>(m_strideModeGroup->checkedId()) == StrideModeId::Seconds)
        && (m_currentResult->window_stride_s > m_currentInput->window_stride_s + 1.0e-9) )
    {
      const std::string clamped = PhysicalUnits::printToBestTimeUnits(
          m_currentResult->window_stride_s * PhysicalUnits::second, 3 );
      m_strideSeconds->setValueText( WString::fromUTF8(clamped) );
    }

    if( m_moreInfoBtn )
      m_moreInfoBtn->show();
    updateTrialsDisplay();

    // Derived stride seconds (samples mode)
    if( m_strideDerivedTxt )
    {
      const double period = currentSamplePeriod();
      if( (static_cast<StrideModeId>(m_strideModeGroup->checkedId()) == StrideModeId::Samples)
          && (period > 0.0) )
      {
        char sBuf[32], pBuf[32];
        snprintf( sBuf, sizeof(sBuf), "%.3f", m_currentResult->window_stride_s );
        snprintf( pBuf, sizeof(pBuf), "%.3f", period );
        m_strideDerivedTxt->setText( WString::tr("dlsd-stride-derived")
            .arg( WString::fromUTF8(sBuf) )
            .arg( WString::fromUTF8(pBuf) ) );
      }
    }

    // ----- Update the chart with ROI highlights + a peak Gaussian -----
    //
    // Scale the displayed spectrum to T_w worth of counts (counts_orig ·
    // T_w / real_time_orig) so the chart shows what an actual transit-window
    // observation would look like.  The peak area is then literally
    //   signal_counts_per_bq · activity_bq
    // i.e. the source counts integrated across the transit window centered
    // on closest approach.
    //
    // `update_spectrum_for_currie_result` derives the peak area from the
    // CurrieMdaResult's `source_counts` / `upper_limit` / `detection_limit`;
    // we override those to steer it into its "upper-bound" branch with our
    // own peak area.
    if( m_spectrum && m_peakModel )
    {
      try
      {
        DetectionLimitCalc::CurrieMdaInput cin = build_curie_input();
        if( cin.spectrum && (cin.spectrum->real_time() > 0.0f) )
        {
          const double T_w = std::max( 1.0e-9, m_currentResult->window_length_s );
          const double real_time_orig = cin.spectrum->real_time();
          const double scale = T_w / real_time_orig;

          // Build a scaled copy of the background measurement so the chart
          // visualizes T_w seconds of data.
          std::shared_ptr<SpecUtils::Measurement> scaled
              = std::make_shared<SpecUtils::Measurement>( *cin.spectrum );
          const std::shared_ptr<const std::vector<float>> orig_counts = scaled->gamma_counts();
          if( orig_counts )
          {
            auto new_counts = std::make_shared<std::vector<float>>( *orig_counts );
            for( float &c : *new_counts )
              c *= static_cast<float>( scale );
            const float new_live = static_cast<float>( scaled->live_time() * scale );
            const float new_real = static_cast<float>( real_time_orig * scale );
            scaled->set_gamma_counts( new_counts, new_live, new_real );
          }
          cin.spectrum = scaled;
          m_spectrum->setData( scaled, true );

          DetectionLimitCalc::CurrieMdaResult cres = DetectionLimitCalc::currie_mda_calc( cin );

          // Peak area = source counts deposited in T_w (no further scaling).
          const double peak_counts = m_currentResult->signal_counts_per_bq
                                   * std::max( 0.0, m_currentResult->activity_bq );
          const double gammas_per_bq = m_currentResult->signal_counts_per_bq;

          // Steer the helper into its "upper-bound" branch:
          //   source_counts <= decision_threshold  AND  upper_limit > 0.
          cres.source_counts      = 0.0f;
          cres.decision_threshold = static_cast<float>( std::max( 1.0e-6, 0.5*peak_counts ) );
          cres.upper_limit        = static_cast<float>( peak_counts );
          cres.lower_limit        = 0.0f;
          cres.detection_limit    = static_cast<float>( peak_counts );

          std::vector<DetectionLimitTool::CurrieResultPeak> peaks;
          const std::shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
          DetectionLimitTool::CurrieResultPeak p;
          p.energy = m_currentResult->input.gamma_energy;
          if( drf && drf->hasResolutionInfo() )
            p.fwhm = drf->peakResolutionFWHM( p.energy );
          else
            p.fwhm = 0.05 * p.energy;
          p.counts_4pi = std::max( 1.0e-30, peak_counts );
          peaks.push_back( p );

          DetectionLimitTool::update_spectrum_for_currie_result(
              m_spectrum, m_peakModel, cin, &cres,
              drf, DetectionLimitTool::LimitType::Activity,
              gammas_per_bq, peaks );
        }
      }catch( std::exception & )
      {
        // Chart decoration is best-effort; ignore.
      }
    }
  }
  catch( std::exception &e )
  {
    if( m_calcErrTxt )
      m_calcErrTxt->setText( WString::fromUTF8(e.what()) );

    // Even on failure, still try to draw the ROI / side-channel highlights
    // (without a peak) so the user can see what region the calc was attempting.
    if( m_spectrum && m_peakModel )
    {
      try
      {
        DetectionLimitCalc::CurrieMdaInput cin = build_curie_input();
        if( cin.spectrum )
        {
          const std::shared_ptr<const DetectorPeakResponse> drf = m_detectorDisplay->detector();
          DetectionLimitTool::update_spectrum_for_currie_result(
              m_spectrum, m_peakModel, cin, nullptr,
              drf, DetectionLimitTool::LimitType::Activity,
              0.0, {} );
        }
      }catch( std::exception & ) {}
    }
  }
}//updateCalcResult()


void DetectionLimitDynamic::runSearch()
{
  m_currentSearch.reset();
  if( m_searchErrTxt )
    m_searchErrTxt->setText( "" );
  if( m_hitsModel )
  {
    while( m_hitsModel->rowCount() > 0 )
      m_hitsModel->removeRow( 0 );
  }

  if( !m_currentInput || !m_currentResult )
  {
    if( m_searchErrTxt )
      m_searchErrTxt->setText( WString::tr("dlsd-search-no-mda") );
    return;
  }

  const std::shared_ptr<const SpecUtils::SpecFile> file = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  if( !file )
  {
    if( m_searchErrTxt )
      m_searchErrTxt->setText( WString::tr("dlsd-search-no-foreground") );
    return;
  }

  // Offload the actual search to a worker thread; long passthrough files take
  // a non-trivial amount of time and would otherwise freeze the session.
  if( m_runSearchBtn )
    m_runSearchBtn->disable();
  if( m_searchInfoTxt )
    m_searchInfoTxt->setText( WString::tr("dlsd-search-running") );

  const std::vector<std::string> dets = file->gamma_detector_names();
  const std::set<int> background_samples; // empty ⇒ use mda_result.input.background_spectrum
  const std::shared_ptr<const DetectionLimitCalc::DynamicMdaResult> mda_result = m_currentResult;
  const std::string sessionid = wApp->sessionId();

  // Output goes via shared_ptr so both worker and GUI-side update can share it.
  auto result_holder = std::make_shared<DetectionLimitCalc::DynamicSearchResult>();
  auto error_msg = std::make_shared<std::string>();

  Wt::Core::observing_ptr<DetectionLimitDynamic> weak_this( this );

  auto update_gui = [weak_this, result_holder, error_msg]()
  {
    if( !weak_this )
      return;
    DetectionLimitDynamic *self = weak_this.get();
    if( self->m_runSearchBtn )
      self->m_runSearchBtn->enable();

    if( !error_msg->empty() )
    {
      if( self->m_searchErrTxt )
        self->m_searchErrTxt->setText( WString::fromUTF8(*error_msg) );
      if( self->m_searchInfoTxt )
        self->m_searchInfoTxt->setText( WString() );
      wApp->triggerUpdate();
      return;
    }

    self->m_currentSearch = result_holder;

    if( !result_holder->warning.empty() && self->m_searchErrTxt )
      self->m_searchErrTxt->setText( WString::fromUTF8(result_holder->warning) );

    for( const DetectionLimitCalc::DynamicSearchHit &h : result_holder->hits )
    {
      const int row = self->m_hitsModel->rowCount();
      char ts[32], rt[32], cnt[32], bg[32], lc[32], sig[32];
      snprintf( ts, sizeof(ts),  "%.2f s", h.window_start_s );
      snprintf( rt, sizeof(rt),  "%.2f s", h.window_real_time_s );
      snprintf( cnt, sizeof(cnt),"%.1f",  h.counts_in_roi );
      snprintf( bg, sizeof(bg),  "%.2f",  h.expected_background );
      snprintf( lc, sizeof(lc),  "%.2f",  h.decision_threshold );
      snprintf( sig, sizeof(sig),"%.2f",  h.n_sigma_excess );

      self->m_hitsModel->insertRow( row );
      self->m_hitsModel->setItem( row, 0, std::make_unique<WStandardItem>( std::to_string(h.first_sample) ) );
      self->m_hitsModel->setItem( row, 1, std::make_unique<WStandardItem>( std::to_string(h.last_sample) ) );
      self->m_hitsModel->setItem( row, 2, std::make_unique<WStandardItem>( std::string(ts) ) );
      self->m_hitsModel->setItem( row, 3, std::make_unique<WStandardItem>( std::string(rt) ) );
      self->m_hitsModel->setItem( row, 4, std::make_unique<WStandardItem>( std::string(cnt) ) );
      self->m_hitsModel->setItem( row, 5, std::make_unique<WStandardItem>( std::string(bg) ) );
      self->m_hitsModel->setItem( row, 6, std::make_unique<WStandardItem>( std::string(lc) ) );
      self->m_hitsModel->setItem( row, 7, std::make_unique<WStandardItem>( std::string(sig) ) );
    }

    char twBuf[32], tsBuf[32];
    snprintf( twBuf, sizeof(twBuf), "%.3f", self->m_currentResult ? self->m_currentResult->window_length_s : 0.0 );
    snprintf( tsBuf, sizeof(tsBuf), "%.3f", self->m_currentResult ? self->m_currentResult->window_stride_s : 0.0 );
    if( self->m_searchInfoTxt )
      self->m_searchInfoTxt->setText( WString::tr("dlsd-search-summary")
          .arg( static_cast<int>(result_holder->hits.size()) )
          .arg( static_cast<int>(result_holder->total_windows_evaluated) )
          .arg( WString::fromUTF8(twBuf) )
          .arg( WString::fromUTF8(tsBuf) ) );

    wApp->triggerUpdate();
  };//update_gui

  auto do_work = [file, dets, background_samples, mda_result, result_holder, error_msg, sessionid, update_gui]()
  {
    try
    {
      *result_holder = DetectionLimitCalc::search_passthrough_for_source(
          file, dets, background_samples, *mda_result );
    }catch( std::exception &e )
    {
      *error_msg = e.what();
    }
    Wt::WServer::instance()->post( sessionid, update_gui );
  };//do_work

  Wt::WServer::instance()->ioService().boost::asio::io_service::post( do_work );
}//runSearch()


void DetectionLimitDynamic::clearSearchHighlights()
{
  // Time-chart highlights deferred to a follow-up.  See plan.
}


void DetectionLimitDynamic::showMoreInfoDialog()
{
  if( !m_currentResult || !m_currentInput )
    return;

  if( m_moreInfoWindow )
    return;

  m_moreInfoWindow = SimpleDialog::make( WString::tr("dlsd-more-info-title") );
  m_moreInfoWindow->addButton( WString::tr("Close") );

  WContainerWidget *contents = m_moreInfoWindow->contents()->addNew<WContainerWidget>();
  contents->addStyleClass( "DynamicMdaMoreInfo" );

  WTable *table = contents->addNew<WTable>();
  table->addStyleClass( "DynamicMdaMoreInfoTable" );

  auto addRow = [table]( const WString &label, const std::string &value ){
    const int row = table->rowCount();
    table->elementAt( row, 0 )->addNew<WText>( label );
    table->elementAt( row, 1 )->addNew<WText>( WString::fromUTF8(value) );
  };

  char buf[128];
  snprintf( buf, sizeof(buf), "%.3f s", m_currentResult->window_length_s );    addRow( WString::tr("dlsd-mi-tw"), buf );
  snprintf( buf, sizeof(buf), "%.3f s", m_currentResult->window_stride_s );    addRow( WString::tr("dlsd-mi-tstep"), buf );
  snprintf( buf, sizeof(buf), "%.1f",   m_currentResult->trials_per_hour );    addRow( WString::tr("dlsd-mi-trials-hr"), buf );
  snprintf( buf, sizeof(buf), "%.3g",   m_currentResult->per_trial_alpha );    addRow( WString::tr("dlsd-mi-alpha-trial"), buf );
  snprintf( buf, sizeof(buf), "%.4f",   m_currentResult->k_alpha );            addRow( WString::tr("dlsd-mi-k-alpha"), buf );
  snprintf( buf, sizeof(buf), "%.4f",   m_currentResult->k_beta );             addRow( WString::tr("dlsd-mi-k-beta"), buf );
  snprintf( buf, sizeof(buf), "%.3f",   m_currentResult->background_counts_in_window ); addRow( WString::tr("dlsd-mi-b-counts"), buf );
  snprintf( buf, sizeof(buf), "%.3f",   m_currentResult->thresholds.sigma_B ); addRow( WString::tr("dlsd-mi-sigma-b"), buf );
  snprintf( buf, sizeof(buf), "%.3f",   m_currentResult->thresholds.L_c );     addRow( WString::tr("dlsd-mi-lc"), buf );
  snprintf( buf, sizeof(buf), "%.3f",   m_currentResult->thresholds.L_d );     addRow( WString::tr("dlsd-mi-ld"), buf );
  snprintf( buf, sizeof(buf), "%.4g",   m_currentResult->signal_counts_per_bq ); addRow( WString::tr("dlsd-mi-s-per-bq"), buf );

  const bool use_curie = use_curie_units();
  const std::string actStr = PhysicalUnits::printToBestActivityUnits(
      m_currentResult->activity_bq * SandiaDecay::becquerel, 3, use_curie );
  addRow( WString::tr("dlsd-mi-activity"), actStr );

  // Diagnostic sweep, ASCII'd.
  if( !m_currentResult->sweep.empty() )
  {
    std::ostringstream oss;
    oss << WString::tr("dlsd-mi-x-fx").toUTF8() << "\n";
    for( const auto &p : m_currentResult->sweep )
      oss << p.first << ", " << p.second << "\n";
    WText *swp = contents->addNew<WText>( WString::fromUTF8(oss.str()) );
    swp->setStyleClass( "MdaSweepCsv" );
    swp->setInline( false );
  }

  m_moreInfoWindow->finished().connect( this, [this, win=m_moreInfoWindow.get()](){
    handleMoreInfoClose( win );
  } );
}


void DetectionLimitDynamic::handleMoreInfoClose( SimpleDialog *dialog )
{
  if( m_moreInfoWindow.get() == dialog )
    m_moreInfoWindow = nullptr;
}


// ===========================================================================
// State serialization (URLs)
// ===========================================================================
std::string DetectionLimitDynamic::encodeStateToUrl() const
{
  std::ostringstream oss;
  // Path = tab name
  const int curTab = m_tabs ? m_tabs->currentIndex() : 0;
  oss << ((curTab == 1) ? "search" : "calculator");
  oss << "?VER=1";

  const SandiaDecay::Nuclide *nuc = m_nucEnterController ? m_nucEnterController->nuclide() : nullptr;
  if( nuc )
  {
    oss << "&NUC=" << nuc->symbol;
    const WString ageStr = m_nucEnterController->nuclideAgeStr();
    if( !ageStr.empty() )
      oss << "&AGE=" << ageStr.toUTF8();
  }

  // URL-encode each free-text value so that spaces, slashes, etc. (e.g. "1 m/s")
  // round-trip correctly.  Numeric values don't need encoding but cost ~nothing.
  auto enc = []( const std::string &v ){ return Wt::Utils::urlEncode(v); };

  oss << "&ENERGY=" << photopeakEnergy();
  oss << "&LROI="   << m_lowerRoi->value();
  oss << "&UROI="   << m_upperRoi->value();
  oss << "&NSIDE="  << m_numSideChannel->value();
  oss << "&CL="     << m_confidenceLevel->currentIndex();
  oss << "&FPHR="   << m_fpPerHour->value();
  oss << "&SOLVE="  << m_solveForGroup->checkedId();
  oss << "&V="      << enc( m_speedInput->text().toUTF8() );
  oss << "&D0="     << enc( m_dcaInput->text().toUTF8() );
  oss << "&A="      << enc( m_activityInput->text().toUTF8() );
  oss << "&WMODE="  << m_windowModeGroup->checkedId();
  oss << "&WLEN="   << enc( m_windowLengthManual->text().toUTF8() );
  oss << "&SMODE="  << m_strideModeGroup->checkedId();
  oss << "&SSEC="   << enc( m_strideSeconds->text().toUTF8() );
  oss << "&SSAMP="  << m_strideSamples->value();
  oss << "&SNAP="   << (m_snapSamplesCb->isChecked() ? "1" : "0");

  // Shielding (only emit keys when populated, to keep URLs compact)
  if( m_shieldingSelect )
  {
    if( m_shieldingSelect->isGenericMaterial() )
    {
      const std::string an = m_shieldingSelect->atomicNumberEdit()->text().toUTF8();
      const std::string ad = m_shieldingSelect->arealDensityEdit()->text().toUTF8();
      if( !an.empty() || !ad.empty() )
      {
        oss << "&SHGEN=1";
        oss << "&SHAN=" << enc(an);
        oss << "&SHAD=" << enc(ad);
      }
    }else
    {
      const std::string mat = m_shieldingSelect->materialEdit()->text().toUTF8();
      if( !mat.empty() )
      {
        const std::string th = m_shieldingSelect->thicknessEdit()->text().toUTF8();
        oss << "&SHMAT=" << enc(mat);
        oss << "&SHTH="  << enc(th);
      }
    }
  }

  return oss.str();
}//encodeStateToUrl()


void DetectionLimitDynamic::handleAppUrl( std::string uri )
{
  std::string host, path, query, fragment;
  AppUtils::split_uri( uri, host, path, query, fragment );
  if( path.empty() )
    path = host;

  if( m_tabs )
  {
    if( SpecUtils::istarts_with(path, "SEARCH") )
      m_tabs->setCurrentIndex( 1 );
    else
      m_tabs->setCurrentIndex( 0 );
  }

  const std::map<std::string,std::string> vals = AppUtils::query_str_key_values( query );
  auto get = [&vals]( const char *k ) -> std::string {
    auto it = vals.find( k );
    if( it == vals.end() )
      return std::string();
    return Wt::Utils::urlDecode( it->second );
  };

  // Version check.  Bumping the format requires either a downgrade path here
  // or a refusal — we currently only know version 1.
  {
    const std::string verStr = get("VER");
    if( !verStr.empty() && verStr != "1" )
      throw std::runtime_error( "Dynamic MDA state URI has unsupported VER=" + verStr );
  }

  const std::string nucStr = get("NUC");
  if( !nucStr.empty() && m_nucEnterController )
    m_nucEnterController->setNuclideText( nucStr );

  const std::string ageStr = get("AGE");
  if( !ageStr.empty() && m_nucEnterController )
    m_nucEnterController->setNuclideAgeTxt( ageStr );

  refreshGammasForCurrentNuclide();

  auto parseFloat = [&]( const char *k, float &out ){
    const std::string s = get(k);
    if( !s.empty() )
      out = static_cast<float>( atof( s.c_str() ) );
  };
  auto parseInt = [&]( const char *k, int &out ){
    const std::string s = get(k);
    if( !s.empty() )
      out = atoi( s.c_str() );
  };

  float energy = photopeakEnergy();
  parseFloat( "ENERGY", energy );
  // Select the closest matching gamma in the combo
  if( energy > 10.0f && !m_photoPeakEnergiesAndBr.empty() )
  {
    size_t best = 0;
    double bestDiff = 1.0e9;
    for( size_t i = 0; i < m_photoPeakEnergiesAndBr.size(); ++i )
    {
      const double d = fabs( m_photoPeakEnergiesAndBr[i].first - energy );
      if( d < bestDiff ) { bestDiff = d; best = i; }
    }
    if( m_photoPeakEnergy )
      m_photoPeakEnergy->setCurrentIndex( static_cast<int>(best) );
  }

  float fval;
  fval = m_lowerRoi->value(); parseFloat( "LROI", fval ); m_lowerRoi->setValue( fval );
  fval = m_upperRoi->value(); parseFloat( "UROI", fval ); m_upperRoi->setValue( fval );

  int ival;
  ival = m_numSideChannel->value(); parseInt( "NSIDE", ival ); m_numSideChannel->setValue( std::max(1, std::min(64, ival)) );
  ival = m_confidenceLevel->currentIndex(); parseInt( "CL", ival ); m_confidenceLevel->setCurrentIndex( std::max(0, std::min<int>(NumConfidenceLevel-1, ival)) );

  fval = m_fpPerHour->value(); parseFloat( "FPHR", fval ); m_fpPerHour->setValue( fval );

  ival = m_solveForGroup->checkedId(); parseInt( "SOLVE", ival );
  if( WRadioButton *btn = m_solveForGroup->button( ival ) )
    m_solveForGroup->setCheckedButton( btn );

  {
    const std::string s = get("V");
    if( !s.empty() )
      m_speedInput->setValueText( WString::fromUTF8(s) );
  }
  {
    const std::string s = get("D0");
    if( !s.empty() )
      m_dcaInput->setValueText( WString::fromUTF8(s) );
  }
  {
    const std::string s = get("A");
    if( !s.empty() )
      m_activityInput->setValueText( WString::fromUTF8(s) );
  }

  ival = m_windowModeGroup->checkedId(); parseInt( "WMODE", ival );
  if( WRadioButton *btn = m_windowModeGroup->button( ival ) )
    m_windowModeGroup->setCheckedButton( btn );
  {
    const std::string s = get("WLEN");
    if( !s.empty() )
      m_windowLengthManual->setValueText( WString::fromUTF8(s) );
  }

  ival = m_strideModeGroup->checkedId(); parseInt( "SMODE", ival );
  if( WRadioButton *btn = m_strideModeGroup->button( ival ) )
    m_strideModeGroup->setCheckedButton( btn );
  {
    const std::string s = get("SSEC");
    if( !s.empty() )
      m_strideSeconds->setValueText( WString::fromUTF8(s) );
  }
  ival = m_strideSamples->value(); parseInt( "SSAMP", ival ); m_strideSamples->setValue( std::max(1, ival) );

  const std::string snap = get("SNAP");
  m_snapSamplesCb->setChecked( !snap.empty() && (snap[0] == '1') );

  // Shielding
  if( m_shieldingSelect )
  {
    const std::string shGen = get("SHGEN");
    if( !shGen.empty() && (shGen[0] == '1') )
    {
      const std::string an = get("SHAN");
      const std::string ad = get("SHAD");
      if( !an.empty() || !ad.empty() )
        m_shieldingSelect->setAtomicNumberAndArealDensity( an, ad );
    }else
    {
      const std::string mat = get("SHMAT");
      if( !mat.empty() )
      {
        const std::string th = get("SHTH");
        m_shieldingSelect->setMaterialNameAndThickness( mat, th );
      }
    }
  }

  handleSolveForChanged();
  handleWindowModeChanged();
  handleStrideModeChanged();

  m_renderFlags |= UpdateLimit;
  scheduleRender();
}//handleAppUrl(...)
