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
#include "InterSpec/SkewParamsGrid.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/FitSkewParamsTool.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;


FitSkewParamsTool::FitSkewParamsTool( InterSpec *viewer )
  : WContainerWidget(),
    m_viewer( viewer ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_skewTypeCombo( nullptr ),
    m_paramsGrid( nullptr ),
    m_updatePeaksCb( nullptr ),
    m_fitBtn( nullptr ),
    m_statusText( nullptr ),
    m_isCalculating( false )
{
  assert( m_viewer );

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
  m_skewTypeCombo->activated().connect( this, [this](){
    updateSkewParamRows();
    userEditedSkewValue();
  } );

  // Skew parameter rows container
  m_paramsGrid = controlsDiv->addNew<SkewParamsGrid>( true, false );
  m_paramsGrid->addStyleClass( "FswSkewParams");
  m_paramsGrid->userChanged().connect( this, &FitSkewParamsTool::userEditedSkewValue );

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

  // PeakModel for the chart: owned by this tool's `m_peakModel` (shared_ptr) and co-owned by the chart
  //  via setPeakModel().
  m_peakModel = PeakModel::create();
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

  // Detect peaks (before building the skew rows, whose defaults depend on the peaks)
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
  }//if( m_spectrum )

  // Build skew param rows (their values come from the prefs, where they are for this skew type)
  updateSkewParamRows();

  // Apply current skew values to peaks and display
  if( !m_detectedPeaks.empty() )
  {
    const vector<shared_ptr<const PeakDef>> displayPeaks = applySkewToPeaks( m_detectedPeaks);
    m_peakModel->setPeaks( displayPeaks);
  }

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
  const int skewIdx = m_skewTypeCombo->currentIndex();
  const PeakDef::SkewType skewType = ((skewIdx >= 0) && (skewIdx < static_cast<int>( PeakDef::NumSkewType )))
                                     ? static_cast<PeakDef::SkewType>( skewIdx ) : PeakDef::NoSkew;
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  m_paramsGrid->setSkewType( skewType );
  m_fitBtn->setEnabled( (nparams > 0) && !m_detectedPeaks.empty() );
  if( nparams == 0 )
    return;

  // Starting values come from the spectrum's (or else the detector's) peak-fit prefs, when they are
  //  for this skew type - so, e.g., a GADRAS detector's own tail powers and extents are used.
  const shared_ptr<const SpecMeas> meas = m_viewer->measurment( SpecUtils::SpectrumType::Foreground );
  const shared_ptr<const PeakFitDetPrefs> measPrefs = meas ? meas->peakFitDetPrefs() : nullptr;
  const shared_ptr<const DetectorPeakResponse> drf = meas ? meas->detector() : nullptr;
  const shared_ptr<const PeakFitDetPrefs> drfPrefs = drf ? drf->peakFitDetPrefs() : nullptr;
  const bool measPrefsMatch = measPrefs && (measPrefs->m_peak_skew_type == skewType);
  const bool isGadras = (skewType == PeakDef::SkewType::GadrasGeneric)
                        || (skewType == PeakDef::SkewType::GadrasCZT);

  // The GADRAS powers (their energy dependence) are fit by default only when the peaks can pin them
  //  down: more than one ROI, spanning at least the 100 keV PeakFitLM requires for energy dependence.
  bool peaksSpanEnergy = false;
  if( isGadras && !m_detectedPeaks.empty() )
  {
    set<shared_ptr<const PeakContinuum>> rois;
    double minEnergy = m_detectedPeaks.front()->mean(), maxEnergy = minEnergy;
    for( const shared_ptr<const PeakDef> &peak : m_detectedPeaks )
    {
      rois.insert( peak->continuum() );
      minEnergy = (std::min)( minEnergy, peak->mean() );
      maxEnergy = (std::max)( maxEnergy, peak->mean() );
    }
    peaksSpanEnergy = (rois.size() > 1) && ((maxEnergy - minEnergy) >= 100.0);
  }//if( isGadras && !m_detectedPeaks.empty() )

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType = PeakDef::CoefficientType( PeakDef::SkewPar0 + p );

    if( PeakDef::is_energy_dependent( skewType, coefType ) )
    {
      double range_lower = 0, range_upper = 0, start_val = 0, step_size = 0;
      PeakDef::skew_parameter_range( skewType, coefType, range_lower, range_upper, start_val, step_size );
      const optional<double> prefsLower = measPrefsMatch ? measPrefs->m_lower_energy_skew[p] : optional<double>{};
      const optional<double> prefsUpper = measPrefsMatch ? measPrefs->m_upper_energy_skew[p] : optional<double>{};
      const double lowerVal = prefsLower.value_or( start_val );
      m_paramsGrid->setValue( p, lowerVal, prefsUpper.value_or( lowerVal ) );
    }else
    {
      m_paramsGrid->setValue( p, skew_starting_value( skewType, coefType, measPrefs.get(), drfPrefs.get() ),
                              std::nullopt );
    }

    // Fitting across the whole spectrum is what can pin down the GADRAS powers (their energy
    //  dependence), so they are fit by default here too.
    const bool fitByDefault = PeakDef::skew_parameter_fit_by_default( skewType, coefType )
          || (peaksSpanEnergy && ((coefType == PeakDef::CoefficientType::SkewPar2)
                                  || (coefType == PeakDef::CoefficientType::SkewPar3)));
    m_paramsGrid->setFit( p, fitByDefault );
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

  // Energy-dependent values are linear between the spectrum's lowest and highest energies (the same
  //  anchors PeakFitDetPrefs and PeakFitLM use), evaluated at each ROI's center (as PeakFitLM does).
  const size_t nchannel = m_spectrum ? m_spectrum->num_gamma_channels() : size_t(0);
  const double lowerEnergy = nchannel ? m_spectrum->gamma_channel_lower( 0 ) : 0.0;
  const double upperEnergy = nchannel ? m_spectrum->gamma_channel_upper( nchannel - 1 ) : 0.0;
  const double energySpan = upperEnergy - lowerEnergy;

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

      const double lowerVal = m_paramsGrid->lowerValue( p ).value_or( 0.0 );
      const optional<double> upperVal = m_paramsGrid->upperValue( p );
      double value = lowerVal;
      if( energyDep && upperVal.has_value() && (energySpan > 1.0) )
      {
        const shared_ptr<const PeakContinuum> cont = newPeak->continuum();
        const double energy = cont->energyRangeDefined()
                              ? 0.5*(cont->lowerEnergy() + cont->upperEnergy())
                              : newPeak->mean();
        const double frac = (energy - lowerEnergy) / energySpan;
        value = lowerVal + frac * (upperVal.value() - lowerVal);
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

      peak->setFitFor( coefType, m_paramsGrid->isFit( p ));
    }
  }

  // The refinement fit below keeps the starting skew as given, but a fit could never leave both
  //  GADRAS tail amplitudes at zero (e.g., the shipped generic GADRAS detectors) - so start from the
  //  defaults instead.  GADRAS skew has no energy dependence, so all peaks have the same values.
  if( !inputPeaks.empty() )
  {
    vector<double> startValues( nparams, 0.0 );
    vector<bool> isFit( nparams, false );
    for( size_t p = 0; p < nparams; ++p )
    {
      const PeakDef::CoefficientType coefType
        = static_cast<PeakDef::CoefficientType>(
            static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));
      startValues[p] = inputPeaks.front()->coefficient( coefType );
      isFit[p] = inputPeaks.front()->fitFor( coefType );
    }

    const vector<double> origValues = startValues;
    PeakDef::avoid_stationary_skew_start( skewType, startValues, isFit );
    if( startValues != origValues )
    {
      for( shared_ptr<const PeakDef> &constPeak : inputPeaks )
      {
        shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
        for( size_t p = 0; p < nparams; ++p )
          peak->set_coefficient( startValues[p], static_cast<PeakDef::CoefficientType>(
                                   static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p )) );
      }
    }//if( moved off the stationary point )
  }//if( !inputPeaks.empty() )

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

  // Wt4: capture this widget's id (not a raw `this`) so the GUI post-back can re-resolve it via
  //  findById() - this tool can be closed mid-fit, and WServer::post only guards the session.
  const string thisid = id();

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

    // Post results back to GUI thread; re-resolve the tool via findById() (it may have been closed
    //  mid-fit - WServer::post only guards the session, not this widget) -> avoids use-after-free.
    WServer::instance()->post( sessionId, [thisid, results, cancelFlag](){
      FitSkewParamsTool *self = dynamic_cast<FitSkewParamsTool *>( wApp->domRoot() ? wApp->domRoot()->findById(thisid) : nullptr );
      if( self )
        self->handleFitResults( results, cancelFlag );
    });
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

  // Update the spin boxes from the fit.  The skew relation gives energy-dependent values at the
  //  spectrum's lowest and highest energies (what the spin boxes, and PeakFitDetPrefs, hold); a
  //  single ROI has no relation, so then all its peaks have the same values.
  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( m_skewTypeCombo->currentIndex());
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  shared_ptr<const PeakDef> skewPeak;
  for( size_t i = 0; !skewPeak && (i < m_fitPeaks.size()); ++i )
  {
    if( m_fitPeaks[i]->skewType() == skewType )
      skewPeak = m_fitPeaks[i];
  }

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType
      = static_cast<PeakDef::CoefficientType>(
          static_cast<int>( PeakDef::CoefficientType::SkewPar0 ) + static_cast<int>( p ));

    optional<pair<double,double>> lowerUpper;
    if( m_fitSkewRelation.has_value() && (m_fitSkewRelation->skew_type == skewType) )
    {
      const PeakFitLM::FitPeaksResults::SkewRelation &sr = m_fitSkewRelation.value();
      if( sr.energy_dependent_skew_pars[p].has_value() )
        lowerUpper = sr.energy_dependent_skew_pars[p];  // {lower energy value, upper energy value}
      else if( sr.non_energy_dependent_skew_pars[p].has_value() )  // {value, uncertainty}
        lowerUpper = make_pair( sr.non_energy_dependent_skew_pars[p]->first,
                                sr.non_energy_dependent_skew_pars[p]->first );
    }else if( skewPeak )
    {
      lowerUpper = make_pair( skewPeak->coefficient( coefType ), skewPeak->coefficient( coefType ) );
    }

    if( !lowerUpper.has_value() )
      continue;

    m_paramsGrid->setValue( p, lowerUpper->first, lowerUpper->second );
  }//for( each param )

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

  // Build PeakFitDetPrefs from the current prefs (keeping detector type, FWHM method, etc), with
  //  the skew from the dialog's current values
  const shared_ptr<const PeakFitDetPrefs> oldPrefs = meas->peakFitDetPrefs();
  const shared_ptr<PeakFitDetPrefs> newPrefs = oldPrefs ? make_shared<PeakFitDetPrefs>( *oldPrefs )
                                                        : make_shared<PeakFitDetPrefs>();

  // Skew type from combo
  const int skewIdx = m_skewTypeCombo->currentIndex();
  if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
    newPrefs->m_peak_skew_type = static_cast<PeakDef::SkewType>( skewIdx);
  else
    newPrefs->m_peak_skew_type = PeakDef::NoSkew;

  // Skew param values (the upper value is only given for energy-dependent params)
  for( size_t p = 0; p < std::size( newPrefs->m_lower_energy_skew ); ++p )
  {
    newPrefs->m_lower_energy_skew[p] = m_paramsGrid->lowerValue( p );
    newPrefs->m_upper_energy_skew[p] = m_paramsGrid->upperValue( p );
  }

  newPrefs->m_roi_independent_skew = false; // This tool is for related skew
  newPrefs->m_source = PeakFitDetPrefs::LoadingSource::UserInputInGui;

  // Apply to SpecMeas and notify
  meas->setPeakFitDetPrefs( newPrefs);
  m_viewer->peakFitDetPrefsChanged().emit();

  // Optionally refit analysis peaks with the new skew parameters
  if( m_updatePeaksCb->isChecked() )
  {
    const std::shared_ptr<PeakModel> peakModel = m_viewer->peakModel();
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

  // Build context menu (owned by this tool; Wt4: a bare `new WPopupMenu()` was an ownerless global
  //  widget that leaked on every right-click, and removeFromParent() is a no-op for global widgets).
  if( m_rightClickMenu )
    removeChild( m_rightClickMenu );

  m_rightClickMenu = addChild( std::make_unique<WPopupMenu>() );

  // "Change Continuum" submenu
  auto contMenuOwned = std::make_unique<WPopupMenu>();
  WPopupMenu *contMenu = contMenuOwned.get();
  for( int t = static_cast<int>( PeakContinuum::NoOffset);
       t <= static_cast<int>( PeakContinuum::External);
       t = t + 1 )
  {
    const PeakContinuum::OffsetType ot = static_cast<PeakContinuum::OffsetType>( t);
    WMenuItem *item = contMenu->addItem( WString::tr( PeakContinuum::offset_type_label_tr( ot ) ));
    item->triggered().connect( this, [this, t](){
      changeContinuumTypeNearEnergy( m_rightClickEnergy, t);
    } );
  }
  m_rightClickMenu->addMenu( WString::tr( "fsw-rclick-change-cont" ), std::move(contMenuOwned));

  // "Delete Peak" item
  WMenuItem *delItem = m_rightClickMenu->addItem( WString::tr( "fsw-rclick-delete-peak" ));
  delItem->triggered().connect( this, [this](){
    deletePeakNearEnergy( m_rightClickEnergy);
  } );

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
  WPushButton *cancelBtn = addCloseButtonToFooter( WString::tr( "fsw-cancel-btn" ));
  cancelBtn->clicked().connect( viewer, &InterSpec::closeFitSkewParamsWindow);

  // Also route the finished() signal (escape key, close button) through InterSpec
  finished().connect( viewer, &InterSpec::closeFitSkewParamsWindow);

  // Accept button - routes through InterSpec for undo/redo tracking
  m_acceptBtn = footer()->addNew<WPushButton>( WString::tr( "fsw-accept-btn" ));
  WidgetUtils::applyButtonRole( m_acceptBtn, WidgetUtils::ButtonRole::Affirm );
  m_acceptBtn->addStyleClass( "Wt-btn");
  m_acceptBtn->clicked().connect( viewer, &InterSpec::acceptFitSkewParamsWindow);

  m_tool->resultUpdated().connect( this, [this](){
    m_acceptBtn->setEnabled( m_tool->canAccept());
  } );

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
