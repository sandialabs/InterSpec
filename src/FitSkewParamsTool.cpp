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
#include <set>
#include <cmath>
#include <memory>
#include <cassert>
#include <algorithm>
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
#include "InterSpec/UserPreferences.h"
#include "InterSpec/FitSkewParamsTool.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"

using namespace std;
using namespace Wt;


namespace
{
  /** Refits each ROI of `peaks` with the skew held at its current values; returns, for each input peak
   (in the same order), its refit peak - or nullptr where its ROI could not be refit.
   */
  vector<shared_ptr<const PeakDef>> refit_peaks_with_fixed_skew( const shared_ptr<const SpecUtils::Measurement> &data,
                                                                 const shared_ptr<const DetectorPeakResponse> &drf,
                                                                 const vector<shared_ptr<const PeakDef>> &peaks,
                                                                 const PeakFitUtils::CoarseResolutionType det_type )
  {
    // Copies of the peaks, with the skew not fit, grouped by ROI (as indices into `peaks`)
    vector<shared_ptr<const PeakDef>> fixed_skew( peaks.size() );
    map<shared_ptr<const PeakContinuum>, vector<size_t>> roi_indices;
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      const shared_ptr<PeakDef> peak = make_shared<PeakDef>( *peaks[i] );
      for( size_t k = 0; k < PeakDef::num_skew_parameters( peak->skewType() ); ++k )
        peak->setFitFor( PeakDef::CoefficientType( PeakDef::SkewPar0 + k ), false );
      fixed_skew[i] = peak;
      roi_indices[peak->continuum()].push_back( i );
    }//for( size_t i = 0; i < peaks.size(); ++i )

    vector<shared_ptr<const PeakDef>> answer( peaks.size() );
    if( !data )
      return answer;

    for( pair<const shared_ptr<const PeakContinuum>, vector<size_t>> &roi : roi_indices )
    {
      vector<size_t> &indices = roi.second;
      std::sort( begin(indices), end(indices), [&fixed_skew]( const size_t lhs, const size_t rhs ){
        return fixed_skew[lhs]->mean() < fixed_skew[rhs]->mean();
      } );

      vector<shared_ptr<const PeakDef>> roi_peaks;
      for( const size_t index : indices )
        roi_peaks.push_back( fixed_skew[index] );

      vector<shared_ptr<const PeakDef>> refit
        = refitPeaksThatShareROI( data, drf, roi_peaks, det_type, PeakFitLM::SmallRefinementOnly );
      if( refit.size() != roi_peaks.size() )
        continue;

      std::sort( begin(refit), end(refit), &PeakDef::lessThanByMeanShrdPtr );
      for( size_t i = 0; i < indices.size(); ++i )
        answer[indices[i]] = refit[i];
    }//for( loop over ROIs )

    return answer;
  }//refit_peaks_with_fixed_skew(...)
}//namespace


FitSkewParamsTool::FitSkewParamsTool( InterSpec *viewer, std::vector<FitSkewSpectrumInput> drfSpectra )
  : WContainerWidget(),
    m_viewer( viewer ),
    m_drfInputs( std::move(drfSpectra) ),
    m_currentSpectrum( 0 ),
    m_anchorLowerEnergy( 0.0 ),
    m_anchorUpperEnergy( 0.0 ),
    m_chart( nullptr ),
    m_peakModel( nullptr ),
    m_spectrumCombo( nullptr ),
    m_onlyCharPeaksCb( nullptr ),
    m_onlyDisplayedCb( nullptr ),
    m_skewTypeCombo( nullptr ),
    m_paramsGrid( nullptr ),
    m_updatePeaksCb( nullptr ),
    m_fitBtn( nullptr ),
    m_statusText( nullptr ),
    m_descText( nullptr ),
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


bool FitSkewParamsTool::isDrfMode() const
{
  return !m_drfInputs.empty();
}


const std::map<std::shared_ptr<const PeakDef>,std::shared_ptr<const PeakDef>> &FitSkewParamsTool::drfPeakReplacements() const
{
  return m_drfPeakReplacements;
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

  // The spectra to fit across: the DRF tool's non-background spectra, or else the foreground with the
  //  peaks detected in it.
  // If the DRF tool has any (non-background) peaks selected to use; if not, all peaks are used.
  bool drfAnySelected = false;
  if( isDrfMode() )
  {
    for( const FitSkewSpectrumInput &input : m_drfInputs )
      drfAnySelected = (drfAnySelected || (!input.isBackground && !input.selectedPeaks.empty()));
    const bool anySelected = drfAnySelected;

    for( size_t i = 0; i < m_drfInputs.size(); ++i )
    {
      const FitSkewSpectrumInput &input = m_drfInputs[i];
      if( input.isBackground || !input.spectrum )
        continue;

      if( input.isDisplayedForeground )
        m_currentSpectrum = m_spectra.size();

      SpectrumState state;
      state.spectrum = input.spectrum;
      state.title = input.title;
      state.drfInputIndex = i;
      state.peaks = anySelected ? input.selectedPeaks : input.allPeaks;
      m_spectra.push_back( std::move(state) );
    }//for( loop over DRF tool spectra )

    // With only background spectra there is nothing to fit across, so just use the foreground
    if( m_spectra.empty() )
      m_drfInputs.clear();
  }//if( isDrfMode() )

  // Get a copy of the foreground measurement
  const shared_ptr<const SpecUtils::Measurement> foreground
    = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground);
  if( !isDrfMode() && foreground )
  {
    SpectrumState state;
    state.spectrum = make_shared<SpecUtils::Measurement>( *foreground);

    // Detect peaks
    try
    {
      AnalystChecks::DetectedPeaksOptions peakOpts;
      peakOpts.specType = SpecUtils::SpectrumType::Foreground;
      peakOpts.nonBackgroundPeaksOnly = false;
      const AnalystChecks::DetectedPeakStatus detected
        = AnalystChecks::detected_peaks( peakOpts, m_viewer);
      state.peaks = detected.peaks;
    }catch( std::exception & )
    {
      // Peak detection failed; continue with empty list
    }

    m_spectra.push_back( std::move(state) );
  }//if( !isDrfMode() && foreground )

  // The skew values in the GUI (and PeakFitDetPrefs) are at the foreground's lowest and highest
  //  energies; without a foreground, use the range of the spectra.
  if( foreground && foreground->num_gamma_channels() )
  {
    m_anchorLowerEnergy = foreground->gamma_channel_lower( 0 );
    m_anchorUpperEnergy = foreground->gamma_channel_upper( foreground->num_gamma_channels() - 1 );
  }else
  {
    bool first = true;
    for( const SpectrumState &state : m_spectra )
    {
      const size_t nchannel = state.spectrum ? state.spectrum->num_gamma_channels() : size_t(0);
      if( !nchannel )
        continue;
      const double lower = state.spectrum->gamma_channel_lower( 0 );
      const double upper = state.spectrum->gamma_channel_upper( nchannel - 1 );
      m_anchorLowerEnergy = first ? lower : (std::min)( m_anchorLowerEnergy, lower );
      m_anchorUpperEnergy = first ? upper : (std::max)( m_anchorUpperEnergy, upper );
      first = false;
    }
  }//if( foreground ) / else
  assert( m_currentSpectrum < std::max( m_spectra.size(), size_t(1) ) );

  // Spectrum display on top — flex: 1 fills available space
  m_chart = addNew<D3SpectrumDisplayDiv>();
  m_chart->addStyleClass( "FswChart");
  m_chart->setMinimumSize( 300, 150);
  m_chart->setCompactAxis( true);

  // Controls below the chart: the spectra options on the left, the skew options in the middle, and
  //  the fit status and button on the right - stacked on narrow screens (see FitSkewParamsTool.css).
  WContainerWidget *controlsDiv = addNew<WContainerWidget>();
  controlsDiv->addStyleClass( "FswControls");

  WContainerWidget *spectraOpts = controlsDiv->addNew<WContainerWidget>();
  spectraOpts->addStyleClass( "FswSpectraOpts");
  WContainerWidget *skewOpts = controlsDiv->addNew<WContainerWidget>();
  skewOpts->addStyleClass( "FswSkewOpts");
  WContainerWidget *fitOpts = controlsDiv->addNew<WContainerWidget>();
  fitOpts->addStyleClass( "FswFitOpts");

  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );

  if( isDrfMode() )
  {
    // Spectrum selector, and the choice of peaks and spectra to fit
    if( m_spectra.size() > 1 )
    {
      WLabel *spectrumLabel = spectraOpts->addNew<WLabel>( WString::tr( "fsw-spectrum-label" ));
      spectrumLabel->addStyleClass( "FswLabel");

      m_spectrumCombo = spectraOpts->addNew<WComboBox>();
      m_spectrumCombo->addStyleClass( "FswCombo FswSpectrumCombo");
      spectrumLabel->setBuddy( m_spectrumCombo );
      for( size_t i = 0; i < m_spectra.size(); ++i )
        m_spectrumCombo->addItem( "" );
      updateSpectrumComboText();
      m_spectrumCombo->setCurrentIndex( static_cast<int>( m_currentSpectrum ) );
      m_spectrumCombo->activated().connect( this, [this]( const int index ){
        if( (index >= 0) && (static_cast<size_t>(index) < m_spectra.size()) )
          showSpectrum( static_cast<size_t>(index) );
        updateDescription();
        updateFitButton();
      } );
    }//if( m_spectra.size() > 1 )

    m_onlyCharPeaksCb = spectraOpts->addNew<WCheckBox>( WString::tr( "fsw-only-char-peaks-cb" ));
    m_onlyCharPeaksCb->addStyleClass( "FswOptCb");
    HelpSystem::attachToolTipOn( m_onlyCharPeaksCb, WString::tr( "fsw-only-char-peaks-tt" ), showToolTips );
    m_onlyCharPeaksCb->setChecked( drfAnySelected );
    m_onlyCharPeaksCb->changed().connect( this, &FitSkewParamsTool::handleOnlyCharPeaksToggled );

    if( m_spectra.size() > 1 )
    {
      m_onlyDisplayedCb = spectraOpts->addNew<WCheckBox>( WString::tr( "fsw-only-displayed-cb" ));
      m_onlyDisplayedCb->addStyleClass( "FswOptCb");
      HelpSystem::attachToolTipOn( m_onlyDisplayedCb, WString::tr( "fsw-only-displayed-tt" ), showToolTips );
      m_onlyDisplayedCb->setChecked( false );
      m_onlyDisplayedCb->changed().connect( this, [this](){
        updateDescription();
        updateFitButton();
      } );
    }//if( m_spectra.size() > 1 )
  }//if( isDrfMode() )

  // Skew type row: label + combo side by side
  WContainerWidget *skewRow = skewOpts->addNew<WContainerWidget>();
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
  m_paramsGrid = skewOpts->addNew<SkewParamsGrid>( true, false );
  m_paramsGrid->addStyleClass( "FswSkewParams");
  m_paramsGrid->userChanged().connect( this, &FitSkewParamsTool::userEditedSkewValue );

  // What will be fit, and the fit status, to the left of the "Fit" button (which then stays put)
  WContainerWidget *statusDiv = fitOpts->addNew<WContainerWidget>();
  statusDiv->addStyleClass( "FswStatusDiv");

  m_descText = statusDiv->addNew<WText>();
  m_descText->addStyleClass( "FswDesc");

  m_statusText = statusDiv->addNew<WText>();
  m_statusText->addStyleClass( "FswStatus");

  m_fitBtn = fitOpts->addNew<WPushButton>( WString::tr( "fsw-fit-btn" ));
  m_fitBtn->addStyleClass( "FswFitBtn");
  m_fitBtn->clicked().connect( this, &FitSkewParamsTool::doFit);

  // "Update analysis peaks on accept" checkbox - created here but may be reparented by the window
  m_updatePeaksCb = controlsDiv->addNew<WCheckBox>( WString::tr( isDrfMode() ? "fsw-refit-drf-peaks-cb"
                                                                             : "fsw-update-peaks-cb" ));
  m_updatePeaksCb->addStyleClass( "FswUpdatePeaksCb");
  m_updatePeaksCb->setChecked( false);
  m_updatePeaksCb->changed().connect( this, &FitSkewParamsTool::updateDescription );

  // PeakModel for the chart: owned by this tool's `m_peakModel` (shared_ptr) and co-owned by the chart
  //  via setPeakModel().
  m_peakModel = PeakModel::create();
  m_peakModel->setNoSpecMeasBacking();
  m_chart->setPeakModel( m_peakModel);

  // Initialize from current PeakFitDetPrefs
  shared_ptr<const SpecMeas> meas
    = m_viewer->measurment( SpecUtils::SpectrumType::Foreground);
  shared_ptr<const PeakFitDetPrefs> prefs = meas ? meas->peakFitDetPrefs() : nullptr;
  assert( prefs || isDrfMode() );

  if( prefs )
  {
    const int skewIdx = static_cast<int>( prefs->m_peak_skew_type);
    if( skewIdx >= 0 && skewIdx < static_cast<int>( PeakDef::NumSkewType ) )
      m_skewTypeCombo->setCurrentIndex( skewIdx);
  }

  // Build skew param rows (their values come from the prefs, where they are for this skew type; and
  //  their defaults depend on the peaks)
  updateSkewParamRows();

  // Show the spectrum, with the current skew values applied to its peaks
  if( !m_spectra.empty() )
    showSpectrum( m_currentSpectrum );

  // Disable fit button if no peaks
  updateDescription();
  updateFitButton();

  // Connect ROI edge drag signal for user adjustment of ROI bounds
  m_chart->existingRoiEdgeDragUpdate().connect( this, &FitSkewParamsTool::handleRoiDrag);

  // Connect right-click signal for context menu
  m_chart->rightClicked().connect( this, &FitSkewParamsTool::handleRightClick);

  // Connect shift+drag to delete peaks in dragged range
  m_chart->shiftKeyDragged().connect( this, &FitSkewParamsTool::handleShiftKeyDrag);
}//void initWidgets()


vector<size_t> FitSkewParamsTool::spectraToFit() const
{
  if( m_onlyDisplayedCb && m_onlyDisplayedCb->isChecked() && (m_currentSpectrum < m_spectra.size()) )
    return { m_currentSpectrum };

  vector<size_t> answer( m_spectra.size() );
  for( size_t i = 0; i < m_spectra.size(); ++i )
    answer[i] = i;
  return answer;
}//vector<size_t> spectraToFit() const


size_t FitSkewParamsTool::numPeaksToFit() const
{
  size_t npeaks = 0;
  for( const size_t index : spectraToFit() )
    npeaks += m_spectra[index].peaks.size();
  return npeaks;
}//size_t numPeaksToFit() const


void FitSkewParamsTool::updateFitButton()
{
  // NoSkew is allowed, so the user can compare chi2 and visuals to a skewed fit
  m_fitBtn->setEnabled( (numPeaksToFit() > 0) && !m_isCalculating );
}//void updateFitButton()


void FitSkewParamsTool::updateDescription()
{
  if( !m_descText )
    return;

  const size_t npeaks = numPeaksToFit();

  string desc;
  if( isDrfMode() )
  {
    const bool onlyChar = m_onlyCharPeaksCb && m_onlyCharPeaksCb->isChecked();
    const bool onlyDisplayed = (spectraToFit().size() == 1);
    const char *key = onlyChar ? (onlyDisplayed ? "fsw-desc-char-peaks-displayed" : "fsw-desc-char-peaks")
                               : (onlyDisplayed ? "fsw-desc-all-peaks-displayed" : "fsw-desc-all-peaks");
    desc = WString::tr( key ).arg( static_cast<int>(npeaks) ).arg( static_cast<int>(m_spectra.size()) ).toUTF8();

    const bool anyBackground = std::any_of( begin(m_drfInputs), end(m_drfInputs),
                                            []( const FitSkewSpectrumInput &input ){ return input.isBackground; } );
    if( anyBackground && !onlyDisplayed )
      desc += " " + WString::tr( "fsw-desc-bkg-excluded" ).toUTF8();

    if( onlyChar && !npeaks )
      desc += " " + WString::tr( "fsw-no-selected-peaks" ).toUTF8();

    if( m_updatePeaksCb && m_updatePeaksCb->isChecked() )
      desc += " " + WString::tr( "fsw-desc-refit-drf" ).arg( static_cast<int>(m_drfInputs.size()) ).toUTF8();
  }else
  {
    desc = WString::tr( "fsw-desc-foreground" ).arg( static_cast<int>(npeaks) ).toUTF8();
    if( m_updatePeaksCb && m_updatePeaksCb->isChecked() )
      desc += " " + WString::tr( "fsw-desc-refit-foreground" ).toUTF8();
  }//if( isDrfMode() ) / else

  m_descText->setText( WString::fromUTF8( desc ) );
}//void updateDescription()


void FitSkewParamsTool::peaksChanged()
{
  refreshDisplayedPeaks();
  updateSpectrumComboText();
  updateDescription();
  updateFitButton();
}//void peaksChanged()


void FitSkewParamsTool::setSpectraOptionsEnabled( const bool enabled )
{
  // The peaks, and spectra, to fit cant change while a fit is running
  if( m_onlyCharPeaksCb )
    m_onlyCharPeaksCb->setEnabled( enabled );
  if( m_onlyDisplayedCb )
    m_onlyDisplayedCb->setEnabled( enabled );
}//void setSpectraOptionsEnabled( const bool enabled )


void FitSkewParamsTool::showSpectrum( const size_t index )
{
  assert( index < m_spectra.size() );
  if( index >= m_spectra.size() )
    return;

  m_currentSpectrum = index;
  m_chart->setData( m_spectra[index].spectrum, false );
  m_peakModel->setForeground( m_spectra[index].spectrum );
  refreshDisplayedPeaks();
}//void showSpectrum( const size_t index )


void FitSkewParamsTool::refreshDisplayedPeaks()
{
  if( m_currentSpectrum >= m_spectra.size() )
    return;

  m_peakModel->setPeaks( applySkewToPeaks( m_spectra[m_currentSpectrum].peaks ) );
}//void refreshDisplayedPeaks()


void FitSkewParamsTool::updateSpectrumComboText()
{
  if( !m_spectrumCombo )
    return;

  for( size_t i = 0; i < m_spectra.size(); ++i )
  {
    const int npeaks = static_cast<int>( m_spectra[i].peaks.size() );
    m_spectrumCombo->setItemText( static_cast<int>(i),
                                  WString::tr( "fsw-spectrum-item" ).arg( m_spectra[i].title ).arg( npeaks ) );
  }
}//void updateSpectrumComboText()


void FitSkewParamsTool::handleOnlyCharPeaksToggled()
{
  if( !isDrfMode() || !m_onlyCharPeaksCb || m_isCalculating )
    return;

  // Starts over from the DRF tool's peaks, dropping any edits, and fit results (but not the skew values)
  const bool onlyChar = m_onlyCharPeaksCb->isChecked();
  for( SpectrumState &state : m_spectra )
  {
    const FitSkewSpectrumInput &input = m_drfInputs[state.drfInputIndex];
    state.peaks = onlyChar ? input.selectedPeaks : input.allPeaks;
  }

  m_statusText->setText( WString() );
  peaksChanged();
}//void handleOnlyCharPeaksToggled()


void FitSkewParamsTool::updateSkewParamRows()
{
  const int skewIdx = m_skewTypeCombo->currentIndex();
  const PeakDef::SkewType skewType = ((skewIdx >= 0) && (skewIdx < static_cast<int>( PeakDef::NumSkewType )))
                                     ? static_cast<PeakDef::SkewType>( skewIdx ) : PeakDef::NoSkew;
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  m_paramsGrid->setSkewType( skewType );
  updateFitButton();
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

  // The GADRAS powers (their energy dependence) are fit by default only when the peaks, of all the
  //  spectra, can pin them down: more than one ROI, spanning at least the 100 keV PeakFitLM requires
  //  for energy dependence.
  bool peaksSpanEnergy = false;
  if( isGadras )
  {
    size_t nrois = 0;
    bool havePeak = false;
    double minEnergy = 0.0, maxEnergy = 0.0;
    for( const size_t index : spectraToFit() )
    {
      set<shared_ptr<const PeakContinuum>> rois;
      for( const shared_ptr<const PeakDef> &peak : m_spectra[index].peaks )
      {
        rois.insert( peak->continuum() );
        minEnergy = havePeak ? (std::min)( minEnergy, peak->mean() ) : peak->mean();
        maxEnergy = havePeak ? (std::max)( maxEnergy, peak->mean() ) : peak->mean();
        havePeak = true;
      }
      nrois += rois.size();
    }//for( const size_t index : spectraToFit() )

    peaksSpanEnergy = (nrois > 1) && ((maxEnergy - minEnergy) >= 100.0);
  }//if( isGadras )

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

  // For NoSkew, strip any skew the input peaks (e.g., from the DRF tool) already have
  if( nparams == 0 )
  {
    vector<shared_ptr<const PeakDef>> result = peaks;
    for( shared_ptr<const PeakDef> &peak : result )
    {
      if( peak->skewType() != PeakDef::NoSkew )
      {
        shared_ptr<PeakDef> newPeak = make_shared<PeakDef>( *peak );
        newPeak->setSkewType( PeakDef::NoSkew );
        peak = newPeak;
      }
    }
    return result;
  }//if( nparams == 0 )

  // Energy-dependent values are linear between the anchor energies (the foreground's lowest and
  //  highest energies, as PeakFitDetPrefs uses), evaluated at each ROI's center (as PeakFitLM does).
  const double energySpan = m_anchorUpperEnergy - m_anchorLowerEnergy;

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
        const double frac = (energy - m_anchorLowerEnergy) / energySpan;
        value = lowerVal + frac * (upperVal.value() - lowerVal);

        // A ROI outside the anchors (e.g., of a DRF spectrum with a wider energy range) extrapolates
        double range_lower = 0.0, range_upper = 0.0, start_val = 0.0, step_size = 0.0;
        if( PeakDef::skew_parameter_range( skewType, coefType, range_lower, range_upper, start_val, step_size ) )
          value = std::clamp( value, range_lower, range_upper );
      }

      newPeak->set_coefficient( value, coefType);
    }//for( each param )

    result.push_back( newPeak);
  }//for( each peak )

  return result;
}//vector<shared_ptr<const PeakDef>> applySkewToPeaks(...)


void FitSkewParamsTool::userEditedSkewValue()
{
  refreshDisplayedPeaks();
  updateFitButton();
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

  // The peaks of the spectra to fit - they reflect the user's ROI bound, continuum, and deletion changes
  if( !numPeaksToFit() )
    return;

  m_isCalculating = true;
  m_fitBtn->setEnabled( false);
  setSpectraOptionsEnabled( false );
  m_statusText->setText( WString::tr( "fsw-fitting" ));
  m_resultUpdated.emit(); // So the dialog disables "Accept" until the fit finishes

  // Cancel any previous calculation
  if( m_cancelCalc )
    m_cancelCalc->store( true);
  m_cancelCalc = make_shared<atomic_bool>( false);

  // Apply current GUI skew values to the peaks, and set fitFor flags based on the "fit" checkboxes
  //  (applySkewToPeaks makes copies, so we can modify them).
  m_fitSpectrumIndices = spectraToFit();
  vector<PeakFitLM::SpectrumPeaks> spectra;
  vector<shared_ptr<PeakDef>> allPeaks;
  for( const size_t index : m_fitSpectrumIndices )
  {
    const SpectrumState &state = m_spectra[index];
    PeakFitLM::SpectrumPeaks spectrum;
    spectrum.data = state.spectrum;
    spectrum.peaks = applySkewToPeaks( state.peaks );

    for( const shared_ptr<const PeakDef> &constPeak : spectrum.peaks )
    {
      shared_ptr<PeakDef> peak = const_pointer_cast<PeakDef>( constPeak);
      for( size_t p = 0; p < nparams; ++p )
        peak->setFitFor( PeakDef::CoefficientType( PeakDef::SkewPar0 + p ), m_paramsGrid->isFit( p ));
      allPeaks.push_back( peak );
    }

    spectra.push_back( std::move(spectrum) );
  }//for( const size_t index : m_fitSpectrumIndices )

  // The refinement fit below keeps the starting skew as given, but a fit could never leave both
  //  GADRAS tail amplitudes at zero (e.g., the shipped generic GADRAS detectors) - so start from the
  //  defaults instead.  GADRAS skew has no energy dependence, so all peaks have the same values.
  if( !allPeaks.empty() && (nparams > 0) )
  {
    vector<double> startValues( nparams, 0.0 );
    vector<bool> isFit( nparams, false );
    for( size_t p = 0; p < nparams; ++p )
    {
      const PeakDef::CoefficientType coefType = PeakDef::CoefficientType( PeakDef::SkewPar0 + p );
      startValues[p] = allPeaks.front()->coefficient( coefType );
      isFit[p] = allPeaks.front()->fitFor( coefType );
    }

    const vector<double> origValues = startValues;
    PeakDef::avoid_stationary_skew_start( skewType, startValues, isFit );
    if( startValues != origValues )
    {
      for( const shared_ptr<PeakDef> &peak : allPeaks )
      {
        for( size_t p = 0; p < nparams; ++p )
          peak->set_coefficient( startValues[p], PeakDef::CoefficientType( PeakDef::SkewPar0 + p ) );
      }
    }//if( moved off the stationary point )
  }//if( !allPeaks.empty() )

  // Get detector type
  shared_ptr<const SpecMeas> meas
    = m_viewer->measurment( SpecUtils::SpectrumType::Foreground);
  shared_ptr<const PeakFitDetPrefs> prefs = meas ? meas->peakFitDetPrefs() : nullptr;
  assert( prefs || isDrfMode() );

  optional<PeakFitUtils::CoarseResolutionType> resType = PeakFitUtils::CoarseResolutionType::Unknown;
  if( prefs )
    resType = prefs->m_det_type;

  // Capture data for background thread
  shared_ptr<atomic_bool> cancelFlag = m_cancelCalc;
  const string sessionId = wApp->sessionId();

  WServer *server = WServer::instance();
  if( !server )
  {
    m_isCalculating = false;
    updateFitButton();
    setSpectraOptionsEnabled( true );
    m_statusText->setText( "");
    m_resultUpdated.emit();
    return;
  }

  // Pre-allocate results; WServer::post() safely no-ops if session is gone by the time
  //  the background thread finishes.
  shared_ptr<PeakFitLM::MultiSpectrumFitResults> results
    = make_shared<PeakFitLM::MultiSpectrumFitResults>();

  // Wt4: capture this widget's id (not a raw `this`) so the GUI post-back can re-resolve it via
  //  findById() - this tool can be closed mid-fit, and WServer::post only guards the session.
  const string thisid = id();

  // Run fit in background thread
  server->ioService().boost::asio::io_service::post( std::bind( [=](){
    if( cancelFlag->load() )
      return;

    // fit_peaks_in_spectra_LM catches all exceptions, and reports them in the results
    *results = PeakFitLM::fit_peaks_in_spectra_LM(
        spectra,
        0.0,  // stat_threshold - keep all peaks
        0.0,  // hypothesis_threshold - keep all peaks
        resType,
        skewType,
        PeakFitLM::SmallRefinementOnly,
        {},
        cancelFlag
     );

    // Post results back to GUI thread; re-resolve the tool via findById() (it may have been closed
    //  mid-fit - WServer::post only guards the session, not this widget) -> avoids use-after-free.
    WServer::instance()->post( sessionId, [thisid, results, cancelFlag](){
      FitSkewParamsTool *self = dynamic_cast<FitSkewParamsTool *>( wApp->domRoot() ? wApp->domRoot()->findById(thisid) : nullptr );
      if( self )
        self->handleFitResults( results, cancelFlag );
    });
  }));
}//void doFit()


void FitSkewParamsTool::handleFitResults( const shared_ptr<PeakFitLM::MultiSpectrumFitResults> &results,
                                           const shared_ptr<atomic_bool> &cancelFlag )
{
  // If this calculation was cancelled, ignore results
  if( cancelFlag->load() || cancelFlag != m_cancelCalc )
    return;

  m_isCalculating = false;
  updateFitButton();
  setSpectraOptionsEnabled( true );
  m_resultUpdated.emit(); // Re-enables "Accept", whether or not the fit succeeded

  if( !results
     || (results->status != PeakFitLM::FitPeaksResults::FitPeaksResultsStatus::Success)
     || (results->spectra.size() != m_fitSpectrumIndices.size()) )
  {
    const string errMsg = results ? results->error_message : "Unknown error";
    m_statusText->setText( WString::tr( "fsw-fit-failed" ).arg( errMsg ));
    wApp->triggerUpdate();
    return;
  }

  for( size_t i = 0; i < m_fitSpectrumIndices.size(); ++i )
  {
    assert( m_fitSpectrumIndices[i] < m_spectra.size() );
    m_spectra.at( m_fitSpectrumIndices[i] ).peaks = results->spectra[i].fit_peaks;
  }
  m_fitSkewRelation = results->skew_relation;

  // Update the spin boxes from the fit.  The skew relation gives energy-dependent values linear in
  //  energy; the spin boxes (and PeakFitDetPrefs) hold the values at the anchor energies.  A single ROI
  //  has no relation, so then all its peaks have the same values.
  const PeakDef::SkewType skewType = static_cast<PeakDef::SkewType>( m_skewTypeCombo->currentIndex());
  const size_t nparams = PeakDef::num_skew_parameters( skewType);

  shared_ptr<const PeakDef> skewPeak;
  for( size_t i = 0; !skewPeak && (i < results->spectra.size()); ++i )
  {
    for( const shared_ptr<const PeakDef> &peak : results->spectra[i].fit_peaks )
    {
      if( !skewPeak && (peak->skewType() == skewType) )
        skewPeak = peak;
    }
  }

  for( size_t p = 0; p < nparams; ++p )
  {
    const PeakDef::CoefficientType coefType = PeakDef::CoefficientType( PeakDef::SkewPar0 + p );

    optional<pair<double,double>> lowerUpper;
    if( m_fitSkewRelation.has_value() && (m_fitSkewRelation->skew_type == skewType) )
    {
      const optional<double> lower = m_fitSkewRelation->value_at( p, m_anchorLowerEnergy );
      const optional<double> upper = m_fitSkewRelation->value_at( p, m_anchorUpperEnergy );
      if( lower.has_value() && upper.has_value() )
        lowerUpper = make_pair( lower.value(), upper.value() );
    }else if( skewPeak )
    {
      lowerUpper = make_pair( skewPeak->coefficient( coefType ), skewPeak->coefficient( coefType ) );
    }

    if( !lowerUpper.has_value() )
      continue;

    // Re-expressing the fit at the foreground's energies can extrapolate, if the spectra fit span less
    double range_lower = 0.0, range_upper = 0.0, start_val = 0.0, step_size = 0.0;
    if( PeakDef::skew_parameter_range( skewType, coefType, range_lower, range_upper, start_val, step_size ) )
    {
      lowerUpper->first = std::clamp( lowerUpper->first, range_lower, range_upper );
      lowerUpper->second = std::clamp( lowerUpper->second, range_lower, range_upper );
    }

    m_paramsGrid->setValue( p, lowerUpper->first, lowerUpper->second );
  }//for( each param )

  // Update spectrum display with fit peaks
  peaksChanged();

  // Compute average chi2/dof across unique ROIs, of all the spectra fit
  double chi2dofSum = 0.0;
  int nRois = 0, nPeaks = 0;
  for( const PeakFitLM::FitPeaksResults &spectrumResult : results->spectra )
  {
    set<shared_ptr<const PeakContinuum>> seenContinua;
    for( const shared_ptr<const PeakDef> &pk : spectrumResult.fit_peaks )
    {
      nPeaks += 1;
      if( pk->chi2Defined() && seenContinua.insert( pk->continuum() ).second )
      {
        chi2dofSum += pk->chi2dof();
        nRois += 1;
      }
    }
  }//for( loop over spectra fit )

  const int nSpectraFit = static_cast<int>( results->spectra.size() );
  WString statusMsg = (nSpectraFit > 1)
                       ? WString::tr( "fsw-fit-success-multi" ).arg( nPeaks ).arg( nSpectraFit )
                       : WString::tr( "fsw-fit-success" ).arg( nPeaks );
  if( nRois > 0 )
  {
    char chi2Buf[32];
    snprintf( chi2Buf, sizeof( chi2Buf ), "%.1f", chi2dofSum / nRois);
    statusMsg += WString::tr( "fsw-fit-chi2" ).arg( chi2Buf);
  }

  m_statusText->setText( statusMsg);

  wApp->triggerUpdate();
}//void handleFitResults(...)


void FitSkewParamsTool::acceptResults()
{
  m_drfPeakReplacements.clear();

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

  if( !m_updatePeaksCb->isChecked() )
    return;

  // Refit peaks with the new skew parameters (held fixed), so only amplitude, mean, sigma, and
  //  continuum adjust.
  const shared_ptr<const DetectorPeakResponse> drf = meas->detector();
  shared_ptr<const PeakFitDetPrefs> skewPrefs = meas->peakFitDetPrefs();
  if( !skewPrefs && drf )
    skewPrefs = drf->peakFitDetPrefs();
  assert( skewPrefs );

  if( isDrfMode() )
  {
    // The peaks of all the DRF tool spectra - background ones too, so background subtraction uses
    //  the same peak shape.  A peak whose ROI could not be refit is left as it was.
    for( const FitSkewSpectrumInput &input : m_drfInputs )
    {
      const PeakFitUtils::CoarseResolutionType detType
        = skewPrefs ? skewPrefs->m_det_type : PeakFitUtils::coarse_det_type( input.spectrum, meas );
      const vector<shared_ptr<const PeakDef>> withSkew = applySkewToPeaks( input.allPeaks );
      const vector<shared_ptr<const PeakDef>> refit
        = refit_peaks_with_fixed_skew( input.spectrum, drf, withSkew, detType );
      assert( refit.size() == input.allPeaks.size() );

      for( size_t i = 0; (i < refit.size()) && (i < input.allPeaks.size()); ++i )
      {
        if( refit[i] )
          m_drfPeakReplacements[input.allPeaks[i]] = refit[i];
      }
    }//for( const FitSkewSpectrumInput &input : m_drfInputs )

    return;
  }//if( isDrfMode() )

  // Refit the analysis peaks of the foreground
  const std::shared_ptr<PeakModel> peakModel = m_viewer->peakModel();
  const shared_ptr<const deque<PeakModel::PeakShrdPtr>> currentPeaks
    = peakModel ? peakModel->peaks() : nullptr;
  const shared_ptr<const SpecUtils::Measurement> data
    = m_viewer->displayedHistogram( SpecUtils::SpectrumType::Foreground);

  if( peakModel && currentPeaks && !currentPeaks->empty() && data )
  {
    UndoRedoManager::PeakModelChange peak_undo_creator;

    const PeakFitUtils::CoarseResolutionType skewDetType
      = skewPrefs ? skewPrefs->m_det_type : PeakFitUtils::coarse_det_type( data, meas);

    // Apply fitted skew parameters to copies of the existing analysis peaks, and refit each ROI; a
    //  ROI that could not be refit keeps the new skew, as is.
    const vector<shared_ptr<const PeakDef>> existing( currentPeaks->begin(), currentPeaks->end());
    const vector<shared_ptr<const PeakDef>> withSkew = applySkewToPeaks( existing);
    const vector<shared_ptr<const PeakDef>> refit = refit_peaks_with_fixed_skew( data, drf, withSkew, skewDetType );
    assert( refit.size() == withSkew.size() );

    vector<shared_ptr<const PeakDef>> allRefit;
    for( size_t i = 0; i < withSkew.size(); ++i )
      allRefit.push_back( ((i < refit.size()) && refit[i]) ? refit[i] : withSkew[i] );

    std::sort( allRefit.begin(), allRefit.end(), &PeakDef::lessThanByMeanShrdPtr);
    peakModel->setPeaks( allRefit);
  }//if( have peaks and data )
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

  const shared_ptr<const SpecUtils::Measurement> spectrum
    = (m_currentSpectrum < m_spectra.size()) ? m_spectra[m_currentSpectrum].spectrum : nullptr;

  if( is_final && spectrum && (roi_px > 10.0) )
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
      = dragFitPrefs ? dragFitPrefs->m_det_type : PeakFitUtils::coarse_det_type( spectrum, dragMeas);

    const vector<shared_ptr<const PeakDef>> refitPeaks
      = refitPeaksThatShareROI_LM( spectrum, detector, roiPeaks,
                                   dragDetType, PeakFitLM::SmallRefinementOnly);
    if( !refitPeaks.empty() )
      roiPeaks = refitPeaks;
  }//if( is_final && worth refitting )

  if( is_final )
  {
    // Replace the ROI's peaks with the refit ones (the displayed peaks share continua with ours)
    replaceRoiPeaks( continuum, roiPeaks );
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
  if( m_currentSpectrum >= m_spectra.size() )
    return;

  // Find nearest peak (the displayed peaks are copies, with the same means and continua, of ours)
  vector<shared_ptr<const PeakDef>> &peaks = m_spectra[m_currentSpectrum].peaks;
  vector<shared_ptr<const PeakDef>>::iterator target = end(peaks);
  double minDist = numeric_limits<double>::max();
  for( auto it = begin(peaks); it != end(peaks); ++it )
  {
    const double dist = fabs( (*it)->mean() - energy);
    if( dist < minDist )
    {
      minDist = dist;
      target = it;
    }
  }

  if( target == end(peaks) )
    return;

  peaks.erase( target );

  peaksChanged();
}//void deletePeakNearEnergy(...)


void FitSkewParamsTool::handleShiftKeyDrag( double x0, double x1 )
{
  if( x0 > x1 )
    std::swap( x0, x1);

  if( m_currentSpectrum >= m_spectra.size() )
    return;

  // Remove peaks whose mean is within the dragged range
  vector<shared_ptr<const PeakDef>> &peaks = m_spectra[m_currentSpectrum].peaks;
  const size_t nbefore = peaks.size();
  peaks.erase( std::remove_if( begin(peaks), end(peaks), [x0,x1]( const shared_ptr<const PeakDef> &p ){
    return (p->mean() >= x0) && (p->mean() <= x1);
  } ), end(peaks) );

  if( peaks.size() == nbefore )
    return;

  peaksChanged();
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
  const shared_ptr<const SpecUtils::Measurement> spectrum
    = (m_currentSpectrum < m_spectra.size()) ? m_spectra[m_currentSpectrum].spectrum : nullptr;
  if( spectrum )
  {
    shared_ptr<const DetectorPeakResponse> detector;
    const shared_ptr<const SpecMeas> contMeas = m_viewer
      ? m_viewer->measurment( SpecUtils::SpectrumType::Foreground ) : nullptr;
    shared_ptr<const PeakFitDetPrefs> contFitPrefs = contMeas ? contMeas->peakFitDetPrefs() : nullptr;
    if( !contFitPrefs && contMeas )
      contFitPrefs = contMeas->detector() ? contMeas->detector()->peakFitDetPrefs() : nullptr;
    assert( contFitPrefs);
    const PeakFitUtils::CoarseResolutionType contDetType
      = contFitPrefs ? contFitPrefs->m_det_type : PeakFitUtils::coarse_det_type( spectrum, contMeas);

    const vector<shared_ptr<const PeakDef>> refitPeaks
      = refitPeaksThatShareROI_LM( spectrum, detector, roiPeaks,
                                   contDetType, PeakFitLM::SmallRefinementOnly);
    if( !refitPeaks.empty() )
      roiPeaks = refitPeaks;
  }

  replaceRoiPeaks( oldContinuum, roiPeaks );
}//void changeContinuumTypeNearEnergy(...)


void FitSkewParamsTool::replaceRoiPeaks( const shared_ptr<const PeakContinuum> &continuum,
                                         const vector<shared_ptr<const PeakDef>> &roiPeaks )
{
  if( m_currentSpectrum >= m_spectra.size() )
    return;

  vector<shared_ptr<const PeakDef>> &peaks = m_spectra[m_currentSpectrum].peaks;
  vector<shared_ptr<const PeakDef>> updated;
  for( const shared_ptr<const PeakDef> &p : peaks )
  {
    if( p->continuum() != continuum )
      updated.push_back( p );
  }
  updated.insert( end(updated), begin(roiPeaks), end(roiPeaks) );
  std::sort( begin(updated), end(updated), &PeakDef::lessThanByMeanShrdPtr );
  peaks = std::move( updated );

  peaksChanged();
}//void replaceRoiPeaks(...)


// ============================================================================
// FitSkewParamsWindow
// ============================================================================

FitSkewParamsWindow::FitSkewParamsWindow( InterSpec *viewer, std::vector<FitSkewSpectrumInput> drfSpectra )
  : AuxWindow( WString::tr( "fsw-title" ),
               Wt::WFlags<AuxWindowProperties>( AuxWindowProperties::EnableResize )
               | AuxWindowProperties::DisableCollapse ),
    m_tool( nullptr ),
    m_acceptBtn( nullptr )
{
  rejectWhenEscapePressed();

  m_tool = contents()->addNew<FitSkewParamsTool>( viewer, std::move(drfSpectra) );

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
  m_acceptBtn->setEnabled( m_tool->canAccept() );

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
