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

#include <set>
#include <memory>
#include <cassert>
#include <functional>

#include <boost/asio/io_service.hpp>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WPoint.h>
#include <Wt/WServer.h>
#include <Wt/WMenuItem.h>
#include <Wt/WIOService.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WPopupMenu.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/FitSkewParamsTool.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;


namespace
{
  // Returns the PeakEdit.xml message IDs for label and tooltip of a given skew parameter.
  // Duplicated from PeakFitDetPrefsGui.cpp to keep files self-contained.
  void skew_param_msg_ids( const PeakDef::SkewType skewType, const size_t paramIndex,
                           const char *&labelId, const char *&tooltipId )
  {
    labelId = nullptr;
    tooltipId = nullptr;

    switch( skewType )
    {
      case PeakDef::NoSkew:
      case PeakDef::NumSkewType:
        return;

      case PeakDef::Bortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-bortel-tau"; tooltipId = "pe-tt-skew-bortel-tau"; }
        return;

      case PeakDef::GaussExp:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-gaussexp-k"; tooltipId = "pe-tt-skew-gaussexp-k"; }
        return;

      case PeakDef::CrystalBall:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-crystalball-alpha"; tooltipId = "pe-tt-skew-crystalball-alpha"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-crystalball-n"; tooltipId = "pe-tt-skew-crystalball-n"; }
        return;

      case PeakDef::ExpGaussExp:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-expgaussexp-kl"; tooltipId = "pe-tt-skew-expgaussexp-kl"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-expgaussexp-kh"; tooltipId = "pe-tt-skew-expgaussexp-kh"; }
        return;

      case PeakDef::DoubleSidedCrystalBall:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-dscb-alphalow"; tooltipId = "pe-tt-skew-dscb-alphalow"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-dscb-nlow"; tooltipId = "pe-tt-skew-dscb-nlow"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-dscb-alphahigh"; tooltipId = "pe-tt-skew-dscb-alphahigh"; }
        if( paramIndex == 3 ){ labelId = "pe-label-skew-dscb-nhigh"; tooltipId = "pe-tt-skew-dscb-nhigh"; }
        return;

      case PeakDef::VoigtPlusBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-voigtplusbortel-gamma"; tooltipId = "pe-tt-skew-voigtplusbortel-gamma"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-voigtplusbortel-r"; tooltipId = "pe-tt-skew-voigtplusbortel-r"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-voigtplusbortel-tau"; tooltipId = "pe-tt-skew-voigtplusbortel-tau"; }
        return;

      case PeakDef::GaussPlusBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-gaussplusbortel-r"; tooltipId = "pe-tt-skew-gaussplusbortel-r"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-gaussplusbortel-tau"; tooltipId = "pe-tt-skew-gaussplusbortel-tau"; }
        return;

      case PeakDef::DoubleBortel:
        if( paramIndex == 0 ){ labelId = "pe-label-skew-doublebortel-tau1"; tooltipId = "pe-tt-skew-doublebortel-tau1"; }
        if( paramIndex == 1 ){ labelId = "pe-label-skew-doublebortel-deltatau2"; tooltipId = "pe-tt-skew-doublebortel-deltatau2"; }
        if( paramIndex == 2 ){ labelId = "pe-label-skew-doublebortel-eta"; tooltipId = "pe-tt-skew-doublebortel-eta"; }
        return;
    }//switch( skewType )
  }//void skew_param_msg_ids(...)
}//anonymous namespace


FitSkewParamsTool::FitSkewParamsTool( InterSpec *viewer )
  : WContainerWidget(),
    m_viewer( viewer ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_skewTypeCombo( nullptr ),
    m_paramsDiv( nullptr ),
    m_updatePeaksCb( nullptr ),
    m_fitBtn( nullptr ),
    m_statusText( nullptr ),
    m_isCalculating( false )
{
  assert( m_viewer );

  for( int i = 0; i < 4; ++i )
  {
    m_lowerSpin[i] = nullptr;
    m_upperSpin[i] = nullptr;
    m_fitCb[i] = nullptr;
  }

  m_rightClickMenu = nullptr;
  m_rightClickEnergy = 0.0;

  initWidgets();
}//FitSkewParamsTool constructor


FitSkewParamsTool::~FitSkewParamsTool()
{
  if( m_cancelCalc )
    m_cancelCalc->store( true);
}


Wt::Signal<> &FitSkewParamsTool::resultUpdated()
{
  return m_resultUpdated;
}


bool FitSkewParamsTool::canAccept() const
{
  return !m_isCalculating;
}


Wt::WCheckBox *FitSkewParamsTool::updatePeaksCb()
{
  return m_updatePeaksCb;
}


void FitSkewParamsTool::initWidgets()
{
  addStyleClass( "FitSkewParamsTool");

  InterSpecApp *app = dynamic_cast<InterSpecApp *>( WApplication::instance());
  if( app )
  {
    app->useMessageResourceBundle( "FitSkewParamsTool");
    app->useMessageResourceBundle( "PeakEdit");
  }

  wApp->useStyleSheet( "InterSpec_resources/FitSkewParamsTool.css");

  // Spectrum display on top — flex: 1 fills available space
  m_chart = addNew<D3SpectrumDisplayDiv>();
  m_chart->addStyleClass( "FswChart");
  m_chart->setMinimumSize( 300, 150);
  m_chart->setCompactAxis( true);

  // Controls area below the chart — flex column, centered
  WContainerWidget *controlsDiv = addNew<WContainerWidget>();
  controlsDiv->addStyleClass( "FswControls");

  // Skew type row: label + combo side by side
  WContainerWidget *skewRow = controlsDiv->addNew<WContainerWidget>();
  skewRow->addStyleClass( "FswSkewRow");

  WLabel *skewLabel = skewRow->addNew<WLabel>( WString::tr( "fsw-skew-type-label" ));
  skewLabel->addStyleClass( "FswLabel");

  m_skewTypeCombo = skewRow->addNew<WComboBox>();
  m_skewTypeCombo->addStyleClass( "FswCombo");
  for( int i = 0; i < static_cast<int>( PeakDef::NumSkewType); ++i )
  {
    const PeakDef::SkewType st = static_cast<PeakDef::SkewType>( i);
    m_skewTypeCombo->addItem( WString::fromUTF8( PeakDef::to_label( st ) ));
  }
  m_skewTypeCombo->setCurrentIndex( 0);
  m_skewTypeCombo->activated().connect( std::bind( [this](){
    updateSkewParamRows();
    userEditedSkewValue();
  }));

  // Skew parameter rows container
  m_paramsDiv = controlsDiv->addNew<WContainerWidget>();
  m_paramsDiv->addStyleClass( "FswSkewParams");

  // Fit button and status on one line
  WContainerWidget *fitRow = controlsDiv->addNew<WContainerWidget>();
  fitRow->addStyleClass( "FswFitRow");

  m_fitBtn = fitRow->addNew<WPushButton>( WString::tr( "fsw-fit-btn" ));
  m_fitBtn->addStyleClass( "FswFitBtn");
  m_fitBtn->clicked().connect( this, &FitSkewParamsTool::doFit);

  m_statusText = fitRow->addNew<WText>();
  m_statusText->addStyleClass( "FswStatus");

  // "Update analysis peaks on accept" checkbox - created here but may be reparented by the window
  m_updatePeaksCb = controlsDiv->addNew<WCheckBox>( WString::tr( "fsw-update-peaks-cb" ));
  m_updatePeaksCb->addStyleClass( "FswUpdatePeaksCb");
  m_updatePeaksCb->setChecked( false);

  // Get a copy of the foreground measurement
  shared_ptr<const SpecUtils::Measurement> foreground
    = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground);
  if( foreground )
  {
    m_spectrum = make_shared<SpecUtils::Measurement>( *foreground);
    m_chart->setData( m_spectrum, false);
  }

  // Create standalone PeakModel for the chart
  m_peakModel = new PeakModel();
  m_peakModel->setNoSpecMeasBacking();
  if( m_spectrum )
    m_peakModel->setForeground( m_spectrum);
  m_chart->setPeakModel( m_peakModel);

  // Initialize from current PeakFitDetPrefs
  shared_ptr<const SpecMeas> meas
    = m_viewer->measurment( SpecUtils::SpectrumType::Foreground);
  shared_ptr<const PeakFitDetPrefs> prefs = meas ? meas->peakFitDetPrefs() : nullptr;
  assert( prefs);

  if( prefs )
  {
    const int skewIdx = static_cast<int>( prefs->m_peak_skew_type);
    if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
      m_skewTypeCombo->setCurrentIndex( skewIdx);
  }

  // Build skew param rows
  updateSkewParamRows();

  // If prefs have skew param values, fill them in
  if( prefs )
  {
    const PeakDef::SkewType skewType = prefs->m_peak_skew_type;
    const size_t nparams = PeakDef::num_skew_parameters( skewType);
    for( size_t p = 0; p < nparams; ++p )
    {
      if( m_lowerSpin[p] && prefs->m_lower_energy_skew[p].has_value() )
        m_lowerSpin[p]->setValue( static_cast<float>( prefs->m_lower_energy_skew[p].value() ));
      if( m_upperSpin[p] && prefs->m_upper_energy_skew[p].has_value() )
        m_upperSpin[p]->setValue( static_cast<float>( prefs->m_upper_energy_skew[p].value() ));
    }
  }//if( prefs )

  // Detect peaks and display them
  if( m_spectrum )
  {
    try
    {
      AnalystChecks::DetectedPeaksOptions peakOpts;
      peakOpts.specType = SpecUtils::SpectrumType::Foreground;
      peakOpts.nonBackgroundPeaksOnly = false;
      const AnalystChecks::DetectedPeakStatus detected
        = AnalystChecks::detected_peaks( peakOpts, m_viewer);
      m_detectedPeaks = detected.peaks;
    }catch( std::exception & )
    {
      // Peak detection failed; continue with empty list
    }

    // Apply current skew values to peaks and display
    if( !m_detectedPeaks.empty() )
    {
      const vector<shared_ptr<const PeakDef>> displayPeaks = applySkewToPeaks( m_detectedPeaks);
      m_peakModel->setPeaks( displayPeaks);
    }
  }//if( m_spectrum )

  // Disable fit button if no peaks or NoSkew
  const int skewIdx = m_skewTypeCombo->currentIndex();
  const PeakDef::SkewType st = static_cast<PeakDef::SkewType>( skewIdx);
  const bool canFit = (PeakDef::num_skew_parameters( st ) > 0) && !m_detectedPeaks.empty();
  m_fitBtn->setEnabled( canFit);

  // Connect ROI edge drag signal for user adjustment of ROI bounds
  m_chart->existingRoiEdgeDragUpdate().connect( this, &FitSkewParamsTool::handleRoiDrag);

  // Connect right-click signal for context menu
  m_chart->rightClicked().connect( this, &FitSkewParamsTool::handleRightClick);

  // Connect shift+drag to delete peaks in dragged range
  m_chart->shiftKeyDragged().connect( this, &FitSkewParamsTool::handleShiftKeyDrag);
}//void initWidgets()


void FitSkewParamsTool::updateSkewParamRows()
{
  m_paramsDiv->clear();
  for( int i = 0; i < 4; ++i )
  {
    m_lowerSpin[i] = nullptr;
    m_upperSpin[i] = nullptr;
    m_fitCb[i] = nullptr;
  }

  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx < 0 || skewIdx >= static_cast<int>( PeakDef::NumSkewType ) )
    return;

  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( skewIdx);
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  if( nparams == 0 )
  {
    m_fitBtn->setEnabled( false);
    return;
  }

  m_fitBtn->setEnabled( !m_detectedPeaks.empty());

  // Check if any parameter is energy-dependent
  bool hasEnergyDep = false;
  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType
      = static_cast<PeakDef::CoefficientType>(
          static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));
    if( PeakDef::is_energy_dependent( skewType, coefType ) )
    {
      hasEnergyDep = true;
      break;
    }
  }

  // Grid: name | fit_cb | lower_val | upper_val  (or name | fit_cb | val for non-energy-dep only)
  WContainerWidget *table = m_paramsDiv->addNew<WContainerWidget>();
  if( hasEnergyDep )
    table->addStyleClass( "FswParamTable");
  else
    table->addStyleClass( "FswParamTable FswParamTableSingleCol");

  // Column headers: name | lower | upper | fit_cb  (or name | val | fit_cb)
  if( hasEnergyDep )
  {
    table->addNew<WText>( ""); // name col placeholder
    WText *lowHeader = table->addNew<WText>( WString::tr( "fsw-lower-header" ));
    lowHeader->addStyleClass( "FswColHeader");
    WText *highHeader = table->addNew<WText>( WString::tr( "fsw-upper-header" ));
    highHeader->addStyleClass( "FswColHeader");
    WText *fitHeader = table->addNew<WText>( WString::tr( "fsw-fit-cb-header" ));
    fitHeader->addStyleClass( "FswColHeader");
  }
  else
  {
    table->addNew<WText>( "");
    table->addNew<WText>( ""); // val col placeholder (no header needed for single-col)
    WText *fitHeader = table->addNew<WText>( WString::tr( "fsw-fit-cb-header" ));
    fitHeader->addStyleClass( "FswColHeader");
  }

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType
      = static_cast<PeakDef::CoefficientType>(
          static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));

    const bool energyDep = PeakDef::is_energy_dependent( skewType, coefType);

    double range_lower = 0, range_upper = 0, start_val = 0, step_size = 0;
    PeakDef::skew_parameter_range( skewType, coefType, range_lower, range_upper, start_val, step_size);

    const char *labelMsgId = nullptr;
    const char *tooltipMsgId = nullptr;
    skew_param_msg_ids( skewType, p, labelMsgId, tooltipMsgId);

    // Build tooltip
    WString tooltipText;
    if( tooltipMsgId )
      tooltipText = WString::tr( tooltipMsgId);

    // Parameter label
    WText *paramLabel = table->addNew<WText>(
      labelMsgId ? WString::tr( labelMsgId ) : WString::tr( "fsw-param-label" ).arg( static_cast<int>( p ) ));
    paramLabel->addStyleClass( "FswParamName");

    if( !tooltipText.empty() )
    {
      HelpSystem::attachToolTipOn( paramLabel, tooltipText, true,
                                   HelpSystem::ToolTipPosition::Right);
    }

    if( energyDep )
    {
      // Lower spin
      NativeFloatSpinBox *lowerSpin = table->addNew<NativeFloatSpinBox>();
      lowerSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ));
      lowerSpin->setFormatString( "%.4G");
      lowerSpin->setSpinnerHidden( true);
      lowerSpin->addStyleClass( "FswSpin");
      lowerSpin->setValue( static_cast<float>( start_val ));
      m_lowerSpin[p] = lowerSpin;
      lowerSpin->valueChanged().connect( std::bind( [this]( float ){
        userEditedSkewValue();
      }, std::placeholders::_1 ));

      // Upper spin
      NativeFloatSpinBox *upperSpin = table->addNew<NativeFloatSpinBox>();
      upperSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ));
      upperSpin->setFormatString( "%.4G");
      upperSpin->setSpinnerHidden( true);
      upperSpin->addStyleClass( "FswSpin");
      upperSpin->setValue( static_cast<float>( start_val ));
      m_upperSpin[p] = upperSpin;
      upperSpin->valueChanged().connect( std::bind( [this]( float ){
        userEditedSkewValue();
      }, std::placeholders::_1 ));

      if( !tooltipText.empty() )
      {
        HelpSystem::attachToolTipOn( lowerSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right);
        HelpSystem::attachToolTipOn( upperSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right);
      }
    }
    else
    {
      // Single value spin
      NativeFloatSpinBox *valSpin = table->addNew<NativeFloatSpinBox>();
      valSpin->setRange( static_cast<float>( range_lower ), static_cast<float>( range_upper ));
      valSpin->setFormatString( "%.4G");
      valSpin->setSpinnerHidden( true);
      valSpin->addStyleClass( hasEnergyDep ? "FswSpin FswSpinWide" : "FswSpin");
      valSpin->setValue( static_cast<float>( start_val ));
      m_lowerSpin[p] = valSpin;
      valSpin->valueChanged().connect( std::bind( [this]( float ){
        userEditedSkewValue();
      }, std::placeholders::_1 ));

      if( !tooltipText.empty() )
      {
        HelpSystem::attachToolTipOn( valSpin, tooltipText, true,
                                     HelpSystem::ToolTipPosition::Right);
      }
    }

    // "Fit" checkbox — last column (right side)
    m_fitCb[p] = table->addNew<WCheckBox>();
    m_fitCb[p]->setChecked( true);
  }//for( each skew param )
}//void updateSkewParamRows()


vector<shared_ptr<const PeakDef>> FitSkewParamsTool::applySkewToPeaks(
  const vector<shared_ptr<const PeakDef>> &peaks ) const
{
  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx < 0 || skewIdx >= static_cast<int>( PeakDef::NumSkewType ) )
    return peaks;

  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( skewIdx);
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  if( nparams == 0 )
    return peaks;

  // Determine energy range for interpolation
  double minEnergy = 1e30, maxEnergy = -1e30;
  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    const double mean = peak->mean();
    if( mean < minEnergy ) minEnergy = mean;
    if( mean > maxEnergy ) maxEnergy = mean;
  }

  const double energySpan = maxEnergy - minEnergy;

  vector<shared_ptr<const PeakDef>> result;
  result.reserve( peaks.size());

  for( const shared_ptr<const PeakDef> &peak : peaks )
  {
    shared_ptr<PeakDef> newPeak = make_shared<PeakDef>( *peak);
    newPeak->setSkewType( skewType);

    for( size_t p = 0; p < nparams; ++p )
    {
      const PeakDef::CoefficientType coefType
        = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));

      const bool energyDep = PeakDef::is_energy_dependent( skewType, coefType);

      double value = 0;
      if( energyDep && m_lowerSpin[p] && m_upperSpin[p] && (energySpan > 1.0) )
      {
        // Interpolate between lower and upper based on peak energy
        const double lowerVal = static_cast<double>( m_lowerSpin[p]->value());
        const double upperVal = static_cast<double>( m_upperSpin[p]->value());
        const double frac = (newPeak->mean() - minEnergy) / energySpan;
        value = lowerVal + frac * (upperVal - lowerVal);
      }
      else if( m_lowerSpin[p] )
      {
        value = static_cast<double>( m_lowerSpin[p]->value());
      }

      newPeak->set_coefficient( value, coefType);
    }//for( each param )

    result.push_back( newPeak);
  }//for( each peak )

  return result;
}//vector<shared_ptr<const PeakDef>> applySkewToPeaks(...)


void FitSkewParamsTool::userEditedSkewValue()
{
  if( m_detectedPeaks.empty() )
    return;

  // Use fit peaks if available, otherwise detected peaks
  const vector<shared_ptr<const PeakDef>> &basePeaks
    = m_fitPeaks.empty() ? m_detectedPeaks : m_fitPeaks;

  const vector<shared_ptr<const PeakDef>> displayPeaks = applySkewToPeaks( basePeaks);
  m_peakModel->setPeaks( displayPeaks);

  // Update fit button state
  const int skewIdx = m_skewTypeCombo->currentIndex();
  const PeakDef::SkewType st = static_cast<PeakDef::SkewType>( skewIdx);
  m_fitBtn->setEnabled( PeakDef::num_skew_parameters( st ) > 0 && !m_isCalculating);
}//void userEditedSkewValue()


void FitSkewParamsTool::doFit()
{
  if( m_isCalculating )
    return;

  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx < 0 || skewIdx >= static_cast<int>( PeakDef::NumSkewType ) )
    return;

  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( skewIdx);
  const size_t nparams = PeakDef::num_skew_parameters( skewType);
  if( nparams == 0 )
    return;

  // Use current peaks from the model — reflects user's ROI bound, continuum, and deletion changes
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> modelPeaks = m_peakModel->peaks();
  if( !modelPeaks || modelPeaks->empty() )
    return;

  m_isCalculating = true;
  m_fitBtn->setEnabled( false);
  m_statusText->setText( WString::tr( "fsw-fitting" ));

  // Cancel any previous calculation
  if( m_cancelCalc )
    m_cancelCalc->store( true);
  m_cancelCalc = make_shared<atomic_bool>( false);

  // Apply current GUI skew values to the model's peaks
  const vector<shared_ptr<const PeakDef>> basePeaks( modelPeaks->begin(), modelPeaks->end());
  vector<shared_ptr<const PeakDef>> inputPeaks = applySkewToPeaks( basePeaks);

  // Set fitFor flags based on the "fit" checkboxes
  for( shared_ptr<const PeakDef> &constPeak : inputPeaks )
  {
    // We need mutable access - applySkewToPeaks already made copies
    shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
    for( size_t p = 0; p < nparams; ++p )
    {
      const PeakDef::CoefficientType coefType
        = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));

      const bool shouldFit = m_fitCb[p] && m_fitCb[p]->isChecked();
      peak->setFitFor( coefType, shouldFit);
    }
  }

  // Get detector type
  shared_ptr<const SpecMeas> meas
    = m_viewer->measurment( SpecUtils::SpectrumType::Foreground);
  shared_ptr<const PeakFitDetPrefs> prefs = meas ? meas->peakFitDetPrefs() : nullptr;
  assert( prefs);

  optional<PeakFitUtils::CoarseResolutionType> resType = PeakFitUtils::CoarseResolutionType::Unknown;
  if( prefs )
    resType = prefs->m_det_type;

  // Capture data for background thread
  shared_ptr<const SpecUtils::Measurement> specCopy = m_spectrum;
  shared_ptr<atomic_bool> cancelFlag = m_cancelCalc;
  const string sessionId = wApp->sessionId();

  WServer *server = WServer::instance();
  if( !server )
  {
    m_isCalculating = false;
    m_fitBtn->setEnabled( true);
    m_statusText->setText( "");
    return;
  }

  // Pre-allocate results; WServer::post() safely no-ops if session is gone by the time
  //  the background thread finishes.
  shared_ptr<PeakFitLM::FitPeaksResults> results
    = make_shared<PeakFitLM::FitPeaksResults>();

  // Run fit in background thread
  server->ioService().boost::asio::io_service::post( std::bind( [=](){
    if( cancelFlag->load() )
      return;

    try
    {
      *results = PeakFitLM::fit_peaks_in_spectrum_LM(
        inputPeaks,
        specCopy,
        0.0,  // stat_threshold - keep all peaks
        0.0,  // hypothesis_threshold - keep all peaks
        resType,
        skewType,
        PeakFitLM::SmallRefinementOnly
     );
    }catch( std::exception &e )
    {
      results->status = PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Failure;
      results->error_message = std::string( "Fit failed: " ) + e.what();
    }catch( ... )
    {
      results->status = PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Failure;
      results->error_message = "Fit failed with unknown error";
    }

    // Post results back to GUI thread; safe even if session was destroyed
    WServer::instance()->post( sessionId, [this, results, cancelFlag](){ handleFitResults( results, cancelFlag); });
  }));
}//void doFit()


void FitSkewParamsTool::handleFitResults( const shared_ptr<PeakFitLM::FitPeaksResults> &results,
                                           const shared_ptr<atomic_bool> &cancelFlag )
{
  // If this calculation was cancelled, ignore results
  if( cancelFlag->load() || cancelFlag != m_cancelCalc )
    return;

  m_isCalculating = false;
  m_fitBtn->setEnabled( true);

  if( !results
     || results->status != PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success )
  {
    const string errMsg = results ? results->error_message : "Unknown error";
    m_statusText->setText( WString::tr( "fsw-fit-failed" ).arg( errMsg ));
    return;
  }

  m_fitPeaks = results->fit_peaks;
  m_fitSkewRelation = results->skew_relation;

  // Update spin boxes from SkewRelation if available
  if( m_fitSkewRelation.has_value() )
  {
    const PeakFitLM::FitPeaksResults::SkewRelation &sr = m_fitSkewRelation.value();

    for( size_t p = 0; p < PeakFitLM::FitPeaksResults::SkewRelation::sm_max_num_skew_pars; ++p )
    {
      // Energy-dependent params: pair is {value, uncertainty}
      if( sr.energy_dependent_skew_pars[p].has_value() )
      {
        // The SkewRelation stores values at energy_range lower and upper anchors.
        // energy_dependent_skew_pars[p] = {lower_energy_value, upper_energy_value} ... actually
        // looking at the struct, it's {value, uncertainty}. Let me use fit peaks directly instead.
      }

      // Non-energy-dependent params
      if( sr.non_energy_dependent_skew_pars[p].has_value() )
      {
        if( m_lowerSpin[p] )
          m_lowerSpin[p]->setValue( static_cast<float>( sr.non_energy_dependent_skew_pars[p].value().first ));
      }
    }
  }//if( m_fitSkewRelation )

  // If SkewRelation didn't fully fill in values, extract from fit peaks directly
  // by looking at the lowest and highest energy peaks
  if( !m_fitPeaks.empty() )
  {
    const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( m_skewTypeCombo->currentIndex());
    const size_t nparams = PeakDef::num_skew_parameters( skewType);

    // Find lowest and highest energy peaks
    shared_ptr<const PeakDef> lowestPeak = m_fitPeaks.front();
    shared_ptr<const PeakDef> highestPeak = m_fitPeaks.front();
    for( const shared_ptr<const PeakDef> &pk : m_fitPeaks )
    {
      if( pk->mean() < lowestPeak->mean() )
        lowestPeak = pk;
      if( pk->mean() > highestPeak->mean() )
        highestPeak = pk;
    }

    for( size_t p = 0; p < nparams; ++p )
    {
      const PeakDef::CoefficientType coefType
        = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));
      const bool energyDep = PeakDef::is_energy_dependent( skewType, coefType);

      if( energyDep )
      {
        if( m_lowerSpin[p] )
          m_lowerSpin[p]->setValue( static_cast<float>( lowestPeak->coefficient( coefType ) ));
        if( m_upperSpin[p] )
          m_upperSpin[p]->setValue( static_cast<float>( highestPeak->coefficient( coefType ) ));
      }
      else
      {
        // Non-energy-dep: all peaks should have the same value; use lowest peak's
        if( m_lowerSpin[p] )
          m_lowerSpin[p]->setValue( static_cast<float>( lowestPeak->coefficient( coefType ) ));
      }
    }//for( each param )
  }//if( !m_fitPeaks.empty() )

  // Update spectrum display with fit peaks
  m_peakModel->setPeaks( m_fitPeaks);

  // Compute average chi2/dof across unique ROIs
  double chi2dofSum = 0.0;
  int nRois = 0;
  set<shared_ptr<const PeakContinuum>> seenContinua;
  for( const shared_ptr<const PeakDef> &pk : m_fitPeaks )
  {
    if( pk->chi2Defined() && seenContinua.insert( pk->continuum() ).second )
    {
      chi2dofSum += pk->chi2dof();
      nRois += 1;
    }
  }

  WString statusMsg = WString::tr( "fsw-fit-success" ).arg( static_cast<int>( m_fitPeaks.size() ));
  if( nRois > 0 )
  {
    char chi2Buf[32];
    snprintf( chi2Buf, sizeof( chi2Buf ), "%.1f", chi2dofSum / nRois);
    statusMsg += WString::tr( "fsw-fit-chi2" ).arg( chi2Buf);
  }

  m_statusText->setText( statusMsg);

  m_resultUpdated.emit();

  wApp->triggerUpdate();
}//void handleFitResults(...)


void FitSkewParamsTool::acceptResults()
{
  shared_ptr<SpecMeas> meas = m_viewer
    ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground )
    : nullptr;
  if( !meas )
    return;

  // Build PeakFitDetPrefs from the dialog's current values
  shared_ptr<PeakFitDetPrefs> newPrefs = make_shared<PeakFitDetPrefs>();

  // Keep the detector type from current prefs
  shared_ptr<const PeakFitDetPrefs> oldPrefs = meas->peakFitDetPrefs();
  if( oldPrefs )
    newPrefs->m_det_type = oldPrefs->m_det_type;

  // Skew type from combo
  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
    newPrefs->m_peak_skew_type = static_cast<PeakDef::SkewType>( skewIdx);
  else
    newPrefs->m_peak_skew_type = PeakDef::NoSkew;

  // Skew param values from spin boxes
  const size_t nparams = PeakDef::num_skew_parameters( newPrefs->m_peak_skew_type);
  for( size_t p = 0; p < nparams; ++p )
  {
    if( m_lowerSpin[p] )
      newPrefs->m_lower_energy_skew[p] = static_cast<double>( m_lowerSpin[p]->value());
    if( m_upperSpin[p] )
      newPrefs->m_upper_energy_skew[p] = static_cast<double>( m_upperSpin[p]->value());
  }

  newPrefs->m_roi_independent_skew = false; // This tool is for related skew
  newPrefs->m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;

  // Apply to SpecMeas and notify
  meas->setPeakFitDetPrefs( newPrefs);
  m_viewer->peakFitDetPrefsChanged().emit();

  // Optionally refit analysis peaks with the new skew parameters
  if( m_updatePeaksCb->isChecked() )
  {
    PeakModel *peakModel = m_viewer->peakModel();
    const shared_ptr<const deque<PeakModel::PeakShrdPtr>> currentPeaks
      = peakModel ? peakModel->peaks() : nullptr;
    const shared_ptr<const SpecUtils::Measurement> data
      = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground);

    if( peakModel && currentPeaks && !currentPeaks->empty() && data )
    {
      UndoRedoManager::PeakModelChange peak_undo_creator;

      // Apply fitted skew parameters to copies of the existing analysis peaks
      const vector<shared_ptr<const PeakDef>> existing( currentPeaks->begin(), currentPeaks->end());
      const vector<shared_ptr<const PeakDef>> withSkew = applySkewToPeaks( existing);

      // Lock skew params so refitting only adjusts amplitude, mean, sigma, continuum
      for( const shared_ptr<const PeakDef> &constPeak : withSkew )
      {
        shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
        const size_t nSkew = PeakDef::num_skew_parameters( peak->skewType());
        for( size_t i = 0; i < nSkew; ++i )
        {
          const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( i ));
          peak->setFitFor( ct, false);
        }
      }

      // Group by continuum (ROI) and refit each ROI
      shared_ptr<const DetectorPeakResponse> drf = meas->detector();

      map<shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>> roiGroups;
      for( const shared_ptr<const PeakDef> &pk : withSkew )
        roiGroups[pk->continuum()].push_back( pk);

      shared_ptr<const PeakFitDetPrefs> skewPrefs = meas ? meas->peakFitDetPrefs() : nullptr;
      if( !skewPrefs && drf )
        skewPrefs = drf->peakFitDetPrefs();
      assert( skewPrefs);
      const PeakFitUtils::CoarseResolutionType skewDetType
        = skewPrefs ? skewPrefs->m_det_type : PeakFitUtils::coarse_det_type( data, meas);

      vector<shared_ptr<const PeakDef>> allRefit;
      for( const auto &entry : roiGroups )
      {
        const vector<shared_ptr<const PeakDef>> refit
          = refitPeaksThatShareROI( data, drf, entry.second,
                                    skewDetType, PeakFitLM::SmallRefinementOnly);

        if( refit.size() == entry.second.size() )
          allRefit.insert( allRefit.end(), refit.begin(), refit.end());
        else
          allRefit.insert( allRefit.end(), entry.second.begin(), entry.second.end());
      }//for( each ROI group )

      std::sort( allRefit.begin(), allRefit.end(), &PeakDef::lessThanByMeanShrdPtr);
      peakModel->setPeaks( allRefit);
    }//if( have peaks and data )
  }//if( update peaks )
}//void acceptResults()


void FitSkewParamsTool::handleRoiDrag( double new_lower, double new_upper, double roi_px,
                                       double orig_lower, std::string /*spec_type*/, bool is_final )
{
  if( !m_peakModel )
    return;

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> allPeaks = m_peakModel->peaks();
  if( !allPeaks || allPeaks->empty() )
    return;

  // Find the continuum matching the original lower energy
  double minDe = 999999.9;
  shared_ptr<const PeakContinuum> continuum;
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    const double de = fabs( p->continuum()->lowerEnergy() - orig_lower);
    if( de < minDe )
    {
      minDe = de;
      continuum = p->continuum();
    }
  }

  if( !continuum || minDe > 1.0 )
    return;

  // Preserve the edge that isn't being dragged
  const bool draggingUpper
    = (fabs( new_lower - continuum->lowerEnergy() ) < fabs( new_upper - continuum->upperEnergy() ));
  if( draggingUpper )
    new_lower = continuum->lowerEnergy();
  else
    new_upper = continuum->upperEnergy();

  // Create new continuum with updated range
  shared_ptr<PeakContinuum> newContinuum = make_shared<PeakContinuum>( *continuum);
  newContinuum->setRange( new_lower, new_upper);

  // Build new peaks for this ROI with the new continuum
  vector<shared_ptr<const PeakDef>> roiPeaks, otherPeaks;
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    if( p->continuum() == continuum )
    {
      // Skip peaks whose mean is more than 1 sigma outside the new range
      if( (p->mean() + p->sigma()) < new_lower || (p->mean() - p->sigma()) > new_upper )
        continue;
      shared_ptr<PeakDef> newPeak = make_shared<PeakDef>( *p);
      newPeak->setContinuum( newContinuum);
      roiPeaks.push_back( newPeak);
    }
    else
    {
      otherPeaks.push_back( p);
    }
  }//for( each peak )

  if( roiPeaks.empty() )
    return;

  // Set skew fitFor(false) on all ROI peaks so skew is preserved during refit;
  //  the skew values are already set correctly from the display peaks.
  for( const shared_ptr<const PeakDef> &constPeak : roiPeaks )
  {
    shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
    const size_t nSkew = PeakDef::num_skew_parameters( peak->skewType());
    for( size_t i = 0; i < nSkew; ++i )
    {
      const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>(
        static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( i ));
      peak->setFitFor( ct, false);
    }
  }//for( set skew fitFor to false )

  if( is_final && m_spectrum && (roi_px > 10.0) )
  {
    // Refit the peaks in this ROI (skew locked, mean/sigma/amplitude free)
    shared_ptr<const DetectorPeakResponse> detector;
    const shared_ptr<const SpecMeas> dragMeas = m_viewer
      ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) : nullptr;
    shared_ptr<const PeakFitDetPrefs> dragFitPrefs = dragMeas ? dragMeas->peakFitDetPrefs() : nullptr;
    if( !dragFitPrefs && dragMeas )
      dragFitPrefs = dragMeas->detector() ? dragMeas->detector()->peakFitDetPrefs() : nullptr;
    assert( dragFitPrefs);
    const PeakFitUtils::CoarseResolutionType dragDetType
      = dragFitPrefs ? dragFitPrefs->m_det_type : PeakFitUtils::coarse_det_type( m_spectrum, dragMeas);

    const vector<shared_ptr<const PeakDef>> refitPeaks
      = refitPeaksThatShareROI_LM( m_spectrum, detector, roiPeaks,
                                   dragDetType, PeakFitLM::SmallRefinementOnly);
    if( !refitPeaks.empty() )
      roiPeaks = refitPeaks;
  }//if( is_final && worth refitting )

  if( is_final )
  {
    // Combine refit ROI peaks with other peaks and update model
    vector<shared_ptr<const PeakDef>> combined;
    combined.reserve( otherPeaks.size() + roiPeaks.size());
    combined.insert( combined.end(), otherPeaks.begin(), otherPeaks.end());
    combined.insert( combined.end(), roiPeaks.begin(), roiPeaks.end());
    std::sort( combined.begin(), combined.end(), &PeakDef::lessThanByMeanShrdPtr);
    m_peakModel->setPeaks( combined);

    // Also update m_detectedPeaks to match (remove old ROI peaks, add new ones)
    vector<shared_ptr<const PeakDef>> updatedDetected;
    for( const shared_ptr<const PeakDef> &p : m_detectedPeaks )
    {
      if( p->continuum() != continuum )
        updatedDetected.push_back( p);
    }
    // Add back the refit peaks (without skew, matching the source peaks)
    for( const shared_ptr<const PeakDef> &p : roiPeaks )
      updatedDetected.push_back( p);
    std::sort( updatedDetected.begin(), updatedDetected.end(), &PeakDef::lessThanByMeanShrdPtr);
    m_detectedPeaks = updatedDetected;

    m_fitPeaks.clear();
  }
  else
  {
    // Intermediate drag preview
    m_chart->updateRoiBeingDragged( roiPeaks);
  }
}//void handleRoiDrag(...)


void FitSkewParamsTool::handleRightClick( double energy, double /*counts*/,
                                           int pageX, int pageY, std::string /*ref_line*/ )
{
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> allPeaks = m_peakModel->peaks();
  if( !allPeaks || allPeaks->empty() )
    return;

  // Find the nearest peak whose ROI contains the click energy
  shared_ptr<const PeakDef> nearestPeak;
  double minDist = numeric_limits<double>::max();
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    if( energy >= p->continuum()->lowerEnergy() && energy <= p->continuum()->upperEnergy() )
    {
      const double dist = fabs( p->mean() - energy);
      if( dist < minDist )
      {
        minDist = dist;
        nearestPeak = p;
      }
    }
  }//for( each peak )

  if( !nearestPeak )
    return;

  m_rightClickEnergy = energy;

  // Build context menu
  if( m_rightClickMenu )
    if( m_rightClickMenu ) m_rightClickMenu->removeFromParent();

  m_rightClickMenu = new WPopupMenu();

  // "Change Continuum" submenu
  auto contMenuOwned = std::make_unique<WPopupMenu>();
  WPopupMenu *contMenu = contMenuOwned.get();
  for( int t = static_cast<int>( PeakContinuum::NoOffset);
       t <= static_cast<int>( PeakContinuum::External);
       t = t + 1 )
  {
    const PeakContinuum::OffsetType ot = static_cast<PeakContinuum::OffsetType>( t);
    WMenuItem *item = contMenu->addItem( WString::tr( PeakContinuum::offset_type_label_tr( ot ) ));
    item->triggered().connect( std::bind( [this, t](){
      changeContinuumTypeNearEnergy( m_rightClickEnergy, t);
    }));
  }
  m_rightClickMenu->addMenu( WString::tr( "fsw-rclick-change-cont" ), std::move(contMenuOwned));

  // "Delete Peak" item
  WMenuItem *delItem = m_rightClickMenu->addItem( WString::tr( "fsw-rclick-delete-peak" ));
  delItem->triggered().connect( std::bind( [this](){
    deletePeakNearEnergy( m_rightClickEnergy);
  }));

  m_rightClickMenu->popup( WPoint( pageX, pageY ));
}//void handleRightClick(...)


void FitSkewParamsTool::deletePeakNearEnergy( double energy )
{
  const shared_ptr<const deque<shared_ptr<const PeakDef>>> allPeaks = m_peakModel->peaks();
  if( !allPeaks || allPeaks->empty() )
    return;

  // Find nearest peak
  shared_ptr<const PeakDef> target;
  double minDist = numeric_limits<double>::max();
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    const double dist = fabs( p->mean() - energy);
    if( dist < minDist )
    {
      minDist = dist;
      target = p;
    }
  }

  if( !target )
    return;

  const double targetMean = target->mean();

  // Remove from display model
  m_peakModel->removePeak( target);

  // Remove corresponding peak from m_detectedPeaks (match by mean energy)
  for( auto it = m_detectedPeaks.begin(); it != m_detectedPeaks.end(); ++it )
  {
    if( fabs( (*it)->mean() - targetMean ) < 0.1 )
    {
      m_detectedPeaks.erase( it);
      break;
    }
  }

  // Also remove from m_fitPeaks if present
  for( auto it = m_fitPeaks.begin(); it != m_fitPeaks.end(); ++it )
  {
    if( fabs( (*it)->mean() - targetMean ) < 0.1 )
    {
      m_fitPeaks.erase( it);
      break;
    }
  }
}//void deletePeakNearEnergy(...)


void FitSkewParamsTool::handleShiftKeyDrag( double x0, double x1 )
{
  if( x0 > x1 )
    std::swap( x0, x1);

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> allPeaks = m_peakModel->peaks();
  if( !allPeaks || allPeaks->empty() )
    return;

  // Find peaks whose mean is within the dragged range
  vector<shared_ptr<const PeakDef>> toRemove;
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    if( p->mean() >= x0 && p->mean() <= x1 )
      toRemove.push_back( p);
  }

  if( toRemove.empty() )
    return;

  // Remove from display model
  for( const shared_ptr<const PeakDef> &p : toRemove )
    m_peakModel->removePeak( p);

  // Remove from m_detectedPeaks and m_fitPeaks by matching mean energy
  for( const shared_ptr<const PeakDef> &removed : toRemove )
  {
    const double mean = removed->mean();

    for( auto it = m_detectedPeaks.begin(); it != m_detectedPeaks.end(); ++it )
    {
      if( fabs( (*it)->mean() - mean ) < 0.1 )
      {
        m_detectedPeaks.erase( it);
        break;
      }
    }

    for( auto it = m_fitPeaks.begin(); it != m_fitPeaks.end(); ++it )
    {
      if( fabs( (*it)->mean() - mean ) < 0.1 )
      {
        m_fitPeaks.erase( it);
        break;
      }
    }
  }//for( each removed peak )
}//void handleShiftKeyDrag(...)


void FitSkewParamsTool::changeContinuumTypeNearEnergy( double energy, int continuum_type )
{
  if( continuum_type < 0 || continuum_type > static_cast<int>( PeakContinuum::External ) )
    return;

  const PeakContinuum::OffsetType newType = static_cast<PeakContinuum::OffsetType>( continuum_type);

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> allPeaks = m_peakModel->peaks();
  if( !allPeaks || allPeaks->empty() )
    return;

  // Find the nearest peak whose ROI contains the click energy
  shared_ptr<const PeakDef> nearestPeak;
  double minDist = numeric_limits<double>::max();
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    if( energy >= p->continuum()->lowerEnergy() && energy <= p->continuum()->upperEnergy() )
    {
      const double dist = fabs( p->mean() - energy);
      if( dist < minDist )
      {
        minDist = dist;
        nearestPeak = p;
      }
    }
  }

  if( !nearestPeak )
    return;

  const shared_ptr<const PeakContinuum> oldContinuum = nearestPeak->continuum();

  // Create new continuum with the requested type
  shared_ptr<PeakContinuum> newContinuum = make_shared<PeakContinuum>( *oldContinuum);
  newContinuum->setType( newType);

  // Gather all peaks sharing this ROI and update their continuum
  vector<shared_ptr<const PeakDef>> roiPeaks, otherPeaks;
  for( const shared_ptr<const PeakDef> &p : *allPeaks )
  {
    if( p->continuum() == oldContinuum )
    {
      shared_ptr<PeakDef> newPeak = make_shared<PeakDef>( *p);
      newPeak->setContinuum( newContinuum);
      roiPeaks.push_back( newPeak);
    }
    else
    {
      otherPeaks.push_back( p);
    }
  }

  if( roiPeaks.empty() )
    return;

  // Set skew fitFor(false) on all ROI peaks so skew is preserved during refit
  for( const shared_ptr<const PeakDef> &constPeak : roiPeaks )
  {
    shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
    const size_t nSkew = PeakDef::num_skew_parameters( peak->skewType());
    for( size_t i = 0; i < nSkew; ++i )
    {
      const PeakDef::CoefficientType ct = static_cast<PeakDef::CoefficientType>(
        static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( i ));
      peak->setFitFor( ct, false);
    }
  }//for( set skew fitFor to false )

  // Refit the ROI with the new continuum type (skew locked)
  if( m_spectrum )
  {
    shared_ptr<const DetectorPeakResponse> detector;
    const shared_ptr<const SpecMeas> contMeas = m_viewer
      ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) : nullptr;
    shared_ptr<const PeakFitDetPrefs> contFitPrefs = contMeas ? contMeas->peakFitDetPrefs() : nullptr;
    if( !contFitPrefs && contMeas )
      contFitPrefs = contMeas->detector() ? contMeas->detector()->peakFitDetPrefs() : nullptr;
    assert( contFitPrefs);
    const PeakFitUtils::CoarseResolutionType contDetType
      = contFitPrefs ? contFitPrefs->m_det_type : PeakFitUtils::coarse_det_type( m_spectrum, contMeas);

    const vector<shared_ptr<const PeakDef>> refitPeaks
      = refitPeaksThatShareROI_LM( m_spectrum, detector, roiPeaks,
                                   contDetType, PeakFitLM::SmallRefinementOnly);
    if( !refitPeaks.empty() )
      roiPeaks = refitPeaks;
  }

  // Combine and update model
  vector<shared_ptr<const PeakDef>> combined;
  combined.reserve( otherPeaks.size() + roiPeaks.size());
  combined.insert( combined.end(), otherPeaks.begin(), otherPeaks.end());
  combined.insert( combined.end(), roiPeaks.begin(), roiPeaks.end());
  std::sort( combined.begin(), combined.end(), &PeakDef::lessThanByMeanShrdPtr);
  m_peakModel->setPeaks( combined);

  // Also update m_detectedPeaks to reflect the continuum change
  vector<shared_ptr<const PeakDef>> updatedDetected;
  for( const shared_ptr<const PeakDef> &p : m_detectedPeaks )
  {
    if( p->continuum() != oldContinuum )
      updatedDetected.push_back( p);
  }
  for( const shared_ptr<const PeakDef> &p : roiPeaks )
    updatedDetected.push_back( p);
  std::sort( updatedDetected.begin(), updatedDetected.end(), &PeakDef::lessThanByMeanShrdPtr);
  m_detectedPeaks = updatedDetected;

  m_fitPeaks.clear();
}//void changeContinuumTypeNearEnergy(...)


// ============================================================================
// FitSkewParamsWindow
// ============================================================================

FitSkewParamsWindow::FitSkewParamsWindow( InterSpec *viewer )
  : AuxWindow( WString::tr( "fsw-title" ),
               Wt::WFlags<AuxWindowProperties>( AuxWindowProperties::EnableResize )
               | AuxWindowProperties::DisableCollapse ),
    m_tool( nullptr ),
    m_acceptBtn( nullptr )
{
  rejectWhenEscapePressed();

  m_tool = contents()->addNew<FitSkewParamsTool>( viewer);

  // Move the "update peaks" checkbox to the footer
  WCheckBox *updateCb = m_tool->updatePeaksCb();
  if( updateCb && updateCb->parent() )
  {
    WContainerWidget *oldParent = dynamic_cast<WContainerWidget *>( updateCb->parent());
    if( oldParent )
    {
      std::unique_ptr<Wt::WWidget> cbOwned = oldParent->removeWidget( updateCb);
      footer()->addWidget( std::move(cbOwned));
    }
  }

  // Cancel button - routes through InterSpec for undo/redo tracking
  WPushButton *cancelBtn = addCloseButtonToFooter( WString::tr( "fsw-cancel-btn" ), true);
  cancelBtn->clicked().connect( viewer, &InterSpec::closeFitSkewParamsWindow);

  // Also route the finished() signal (escape key, close button) through InterSpec
  finished().connect( viewer, &InterSpec::closeFitSkewParamsWindow);

  // Accept button - routes through InterSpec for undo/redo tracking
  m_acceptBtn = footer()->addNew<WPushButton>( WString::tr( "fsw-accept-btn" ));
  m_acceptBtn->addStyleClass( "Wt-btn");
  m_acceptBtn->clicked().connect( viewer, &InterSpec::acceptFitSkewParamsWindow);

  m_tool->resultUpdated().connect( std::bind( [this](){
    m_acceptBtn->setEnabled( m_tool->canAccept());
  }));

  resizeScaledWindow( 0.7, 0.7);
  centerWindow();
  show();
}//FitSkewParamsWindow constructor


FitSkewParamsWindow::~FitSkewParamsWindow()
{
}


FitSkewParamsTool *FitSkewParamsWindow::tool()
{
  return m_tool;
}
