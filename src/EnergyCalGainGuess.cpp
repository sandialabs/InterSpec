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
#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <stdexcept>

#include <Wt/WText.h>
#include <Wt/WTable.h>
#include <Wt/WLabel.h>
#include <Wt/WServer.h>
#include <Wt/WTableRow.h>
#include <Wt/WTableCell.h>
#include <Wt/WIOService.h>
#include <Wt/WLineEdit.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WSuggestionPopup.h>
#include <Wt/WContainerWidget.h>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/EnergyCal.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/EnergyCalTool.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/ReferenceLineInfo.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/EnergyCalGainGuess.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/D3SpectrumDisplayDiv.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

using namespace std;

using SpecUtils::Measurement;
using SpecUtils::EnergyCalibration;


namespace
{
  const double sm_annihilation_energy = 510.9989;

  /** Lines below this are too unreliable to match against (huge line density, heavy
   attenuation, and the candidate finder is noisy there).
   */
  const double sm_min_line_energy = 45.0;

  /** A hypothesis must place the top of the spectrum in this range.  Real-world spectra we
   care about span at least ~250 keV full-scale and no more than ~13 MeV, so hypotheses
   outside this are rejected outright.
   */
  const double sm_min_spectrum_span = 250.0;   //keV
  const double sm_max_spectrum_span = 13000.0; //keV

  const double sm_gate_fail_score = -999.0;


  /** Nominal expected FWHM (keV) at an energy - only used for *plausibility* gating, so the
   gates that use it are kept loose (real detectors vary by ~2x from these).
   For HPGe: ~1.05 keV at 122 keV, ~1.9 keV at 1332 keV.
   For NaI-ish: 7.5% at 661.7 keV, scaling as sqrt(E).
   */
  double expected_fwhm_kev( const double energy, const bool is_high_res )
  {
    const double e = std::max( energy, 10.0 );
    if( is_high_res )
      return sqrt( 0.9 + 0.002*e );
    return 1.93 * sqrt( e );
  }


  bool is_cancelled( const shared_ptr<const atomic<bool>> &cancel )
  {
    return cancel && cancel->load();
  }


  /** The compact %.1f formatting used in provenance strings. */
  string energy_str( const double energy )
  {
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.1f", energy );
    return buffer;
  }
}//namespace


namespace GainGuessCalc
{

std::vector<CandidatePeak>
find_candidate_peaks_channelspace( const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                   bool &is_high_res )
{
  is_high_res = true;

  if( !meas )
    throw runtime_error( "find_candidate_peaks_channelspace: null spectrum" );

  const size_t nchan = meas->num_gamma_channels();
  if( nchan < 64 )
    throw runtime_error( "find_candidate_peaks_channelspace: too few channels" );

  // Install an identity calibration on a copy, so candidate means/sigmas come back in channel
  //  units no matter what (possibly nonsense) calibration the file carried.
  auto identity_cal = make_shared<EnergyCalibration>();
  identity_cal->set_default_polynomial( nchan, {0.0f, 1.0f}, {} );

  auto chanmeas = make_shared<Measurement>( *meas );
  chanmeas->set_energy_calibration( identity_cal );

  // The finder is told the resolution class through the prefs (the FWHM-based auto
  //  classification would be meaningless with channel-unit "energies").
  auto run_finder = [&chanmeas]( const PeakFitUtils::CoarseResolutionType type )
        -> vector<tuple<float,float,float>>
  {
    auto prefs = make_shared<PeakFitDetPrefs>();
    prefs->m_det_type = type;
    vector<tuple<float,float,float>> cands; //(mean, sigma, area)
    secondDerivativePeakCanidates( chanmeas, prefs, 0, 0, cands );
    return cands;
  };

  vector<tuple<float,float,float>> raw = run_finder( PeakFitUtils::CoarseResolutionType::High );

  // Classify high vs low/medium resolution from the relative widths of the strongest
  //  candidates: sigma_ch/channel ~= sigma_E/E, which is ~0.1-0.3% for HPGe and ~2-4% for NaI.
  auto median_rel_width = []( vector<tuple<float,float,float>> cands ) -> double {
    std::sort( begin(cands), end(cands),
        []( const tuple<float,float,float> &a, const tuple<float,float,float> &b ){
          return get<2>(a) > get<2>(b);
        } );
    if( cands.size() > 10 )
      cands.resize( 10 );
    vector<double> rel;
    for( const tuple<float,float,float> &c : cands )
      if( (get<0>(c) > 1.0f) && (get<1>(c) > 0.0f) )
        rel.push_back( get<1>(c) / get<0>(c) );
    if( rel.empty() )
      return 0.0;
    std::sort( begin(rel), end(rel) );
    return rel[rel.size()/2];
  };

  const double rel_width = median_rel_width( raw );
  is_high_res = (rel_width > 0.0) && (rel_width < 0.0125);

  if( !is_high_res )
    raw = run_finder( PeakFitUtils::CoarseResolutionType::LowOrMedRes );

  // Build CandidatePeak entries, with a rough significance = area / sqrt(gross counts under
  //  the peak) (the finder does not report an area uncertainty).
  vector<CandidatePeak> peaks;
  for( const tuple<float,float,float> &c : raw )
  {
    CandidatePeak peak;
    peak.channel = get<0>( c );
    peak.sigma_ch = std::max( 0.25, static_cast<double>(get<1>(c)) );
    peak.area = get<2>( c );
    if( (peak.area <= 0.0) || (peak.channel <= 1.0) || (peak.channel >= (nchan - 2)) )
      continue;

    const size_t lo = static_cast<size_t>( std::max( 0.0, peak.channel - 2.0*peak.sigma_ch ) );
    const size_t hi = std::min( nchan - 1,
                        static_cast<size_t>( std::max(0.0, peak.channel + 2.0*peak.sigma_ch) ) );
    const double gross = chanmeas->gamma_channels_sum( lo, hi );
    peak.area_uncert = sqrt( std::max( gross, 1.0 ) );
    peaks.push_back( peak );
  }//for( loop over raw candidates )

  // Drop marginal detections (second-derivative bumps on a Poisson continuum are typically
  //  ~2-4 sigma), keep the ~30 most significant, and mark clearly-real ones "must explain".
  auto significance = []( const CandidatePeak &p ){ return p.area / p.area_uncert; };
  peaks.erase( std::remove_if( begin(peaks), end(peaks),
      [&significance]( const CandidatePeak &p ){ return significance(p) < 4.0; } ),
      end(peaks) );

  std::sort( begin(peaks), end(peaks),
      [&significance]( const CandidatePeak &a, const CandidatePeak &b ){
        return significance(a) > significance(b);
      } );
  if( peaks.size() > 30 )
    peaks.resize( 30 );

  const double max_sig = peaks.empty() ? 0.0 : significance( peaks.front() );
  for( size_t i = 0; (i < peaks.size()) && (i < 15); ++i )
    peaks[i].significant = ( significance(peaks[i]) >= std::max( 8.0, 0.04*max_sig ) );

  std::sort( begin(peaks), end(peaks),
      []( const CandidatePeak &a, const CandidatePeak &b ){ return a.channel < b.channel; } );

  // Flag width outliers (511 keV / x-ray suspects) for high-res spectra: fit sigma^2 vs
  //  channel (linear, one robust re-fit pass), and flag peaks well above the trend.
  if( is_high_res )
  {
    vector<const CandidatePeak *> sig_peaks;
    for( const CandidatePeak &p : peaks )
      if( p.significant )
        sig_peaks.push_back( &p );

    if( sig_peaks.size() >= 6 )
    {
      double a = 0.0, b = 0.0;
      auto fit_line = [&a, &b]( const vector<const CandidatePeak *> &pts ){
        double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
        const double n = static_cast<double>( pts.size() );
        for( const CandidatePeak *p : pts )
        {
          const double x = p->channel, y = p->sigma_ch * p->sigma_ch;
          sx += x; sy += y; sxx += x*x; sxy += x*y;
        }
        const double denom = (n*sxx - sx*sx);
        b = (fabs(denom) > 1.0E-9) ? (n*sxy - sx*sy) / denom : 0.0;
        a = (sy - b*sx) / n;
      };

      fit_line( sig_peaks );

      // One robust pass: refit without points far above the first trend.
      vector<const CandidatePeak *> keep;
      for( const CandidatePeak *p : sig_peaks )
      {
        const double pred = std::max( a + b*p->channel, 0.0625 );
        if( p->sigma_ch*p->sigma_ch < 2.0*pred )
          keep.push_back( p );
      }
      if( keep.size() >= 4 )
        fit_line( keep );

      for( CandidatePeak &p : peaks )
      {
        const double pred = std::max( a + b*p.channel, 0.0625 );
        p.width_outlier = ( p.significant && (p.sigma_ch > 1.3*sqrt(pred)) );
      }
    }//if( enough significant peaks to establish a width trend )
  }//if( is_high_res )

  return peaks;
}//find_candidate_peaks_channelspace(...)


LineLibrary build_line_library( const std::vector<const SandiaDecay::Nuclide *> &user_nuclides )
{
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "build_line_library: decay database not initialized" );

  LineLibrary library;

  // Curated common parents: NORM, common industrial/medical/check sources, and SNM.  Kept
  //  deliberately small - every extra parent raises the accidental-match line density and so
  //  lowers everyones score resolution.
  const char * const common_parents[] = {
    "K40", "Th232", "U238", "U235", "Ra226", "Cs137", "Co60", "Co57", "Ba133", "Eu152",
    "Am241", "Na22", "Mn54", "I131", "Ir192",
    "Np237", "Np239", "Pu238", "Pu239", "Pu240", "Pu241"
  };

  vector<pair<const SandiaDecay::Nuclide *, bool>> parents; //(nuclide, user_supplied)
  for( const char *name : common_parents )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
    if( nuc )
      parents.emplace_back( nuc, false );
  }

  for( const SandiaDecay::Nuclide * const nuc : user_nuclides )
  {
    if( !nuc )
      continue;
    bool already = false;
    for( pair<const SandiaDecay::Nuclide *, bool> &p : parents )
    {
      if( p.first == nuc )
      {
        p.second = true;  //upgrade the common entry to user-supplied (score bonus)
        already = true;
      }
    }
    if( !already )
      parents.emplace_back( nuc, true );
  }//for( user nuclides )

  for( const pair<const SandiaDecay::Nuclide *, bool> &parent : parents )
  {
    const SandiaDecay::Nuclide * const nuc = parent.first;

    const double age = PeakDef::defaultDecayTime( nuc );

    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuc, 1.0E6 * SandiaDecay::becquerel );

    // photons() = gammas + annihilation + decay x-rays; decay x-rays (e.g. the U/Pu x-ray
    //  groups) are often the most prominent feature of an actinide spectrum.
    const vector<SandiaDecay::EnergyRatePair> photons
                  = mix.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );

    double max_rate = 0.0;
    for( const SandiaDecay::EnergyRatePair &p : photons )
      if( p.energy >= sm_min_line_energy )
        max_rate = std::max( max_rate, p.numPerSecond );

    if( max_rate <= 0.0 )
      continue;

    vector<LineSource> nuc_lines;
    for( const SandiaDecay::EnergyRatePair &p : photons )
    {
      if( (p.energy < sm_min_line_energy) || (p.numPerSecond < 0.02*max_rate) )
        continue;

      LineSource line;
      const bool is_annih = ( fabs(p.energy - sm_annihilation_energy) < 0.75 );
      line.type = is_annih ? LineSource::Type::Annihilation : LineSource::Type::NuclideGamma;
      line.energy = p.energy;
      line.rel_intensity = static_cast<float>( p.numPerSecond / max_rate );
      // The parent is kept even for the annihilation line: a beta+ emitter (e.g. Na-22)
      //  producing its 511 is real evidence for that parent in the self-consistency scoring.
      line.nuc = nuc;
      line.user_supplied = parent.second;
      nuc_lines.push_back( line );
    }//for( photons of this parent )

    // Cap the per-parent line count (strongest first) to bound the library density.
    std::sort( begin(nuc_lines), end(nuc_lines),
        []( const LineSource &a, const LineSource &b ){ return a.rel_intensity > b.rel_intensity; } );
    if( nuc_lines.size() > 20 )
      nuc_lines.resize( 20 );

    library.lines.insert( end(library.lines), begin(nuc_lines), end(nuc_lines) );
  }//for( loop over parents )

  // A free-standing annihilation line (e.g. from pair production of hard background gammas),
  //  in case no included parent produced one.
  {
    bool have_annih = false;
    for( const LineSource &l : library.lines )
      have_annih = ( have_annih || (l.type == LineSource::Type::Annihilation) );
    if( !have_annih )
    {
      LineSource line;
      line.type = LineSource::Type::Annihilation;
      line.energy = sm_annihilation_energy;
      line.rel_intensity = 1.0f;
      library.lines.push_back( line );
    }
  }

  // Fluorescence x-rays of the elements that most often dominate shielded/self-attenuating
  //  unknowns; the Ka/Kb spacing is a gain-invariant fingerprint for the pair matcher.
  const int fluor_elements[] = { 74 /*W*/, 82 /*Pb*/, 83 /*Bi*/, 90 /*Th*/, 92 /*U*/ };
  for( const int z : fluor_elements )
  {
    const SandiaDecay::Element * const el = db->element( z );
    if( !el )
      continue;

    vector<SandiaDecay::EnergyIntensityPair> xrays = el->xrays;
    std::sort( begin(xrays), end(xrays),
        []( const SandiaDecay::EnergyIntensityPair &a, const SandiaDecay::EnergyIntensityPair &b ){
          return a.intensity > b.intensity;
        } );

    // Keep the K-series thicket (Ka1/Ka2/Kb1/Kb3/Kb2) - a correct hypothesis must be able to
    //  claim all the x-ray peaks a heavy absorber produces, not just the top three.
    double max_intensity = 0.0;
    size_t nadded = 0;
    for( const SandiaDecay::EnergyIntensityPair &xray : xrays )
    {
      if( (xray.energy < sm_min_line_energy) || (nadded >= 5) )
        continue;
      if( max_intensity <= 0.0 )
        max_intensity = xray.intensity;
      if( xray.intensity < 0.05*max_intensity )
        continue;

      LineSource line;
      line.type = LineSource::Type::FluorXray;
      line.energy = xray.energy;
      line.rel_intensity = static_cast<float>( xray.intensity / max_intensity );
      line.el = el;
      library.lines.push_back( line );
      ++nadded;
    }//for( xrays of this element )
  }//for( fluorescence elements )

  std::sort( begin(library.lines), end(library.lines),
      []( const LineSource &a, const LineSource &b ){ return a.energy < b.energy; } );

  // Pairwise energy-ratio table over the reasonably strong lines.
  vector<uint32_t> strong;
  for( size_t i = 0; i < library.lines.size(); ++i )
    if( library.lines[i].rel_intensity >= 0.1f )
      strong.push_back( static_cast<uint32_t>(i) );

  for( size_t a = 0; a < strong.size(); ++a )
  {
    for( size_t b = a + 1; b < strong.size(); ++b )
    {
      const double e_low = library.lines[strong[a]].energy;
      const double e_high = library.lines[strong[b]].energy;
      const double ratio = e_low / e_high;
      if( (ratio <= 0.02) || (ratio >= 0.995) )
        continue;
      library.ratios.push_back( LineLibrary::RatioEntry{
                    static_cast<float>(ratio), strong[a], strong[b] } );
    }
  }//for( pairs of strong lines )

  std::sort( begin(library.ratios), end(library.ratios),
      []( const LineLibrary::RatioEntry &a, const LineLibrary::RatioEntry &b ){
        return a.ratio < b.ratio;
      } );

  return library;
}//build_line_library(...)


double score_hypothesis( const double offset, const double gain,
                         const std::vector<CandidatePeak> &peaks,
                         const size_t nchannels,
                         const LineLibrary &library,
                         const bool is_high_res,
                         std::vector<PeakAssignment> *assignments,
                         double *explained_frac )
{
  if( assignments )
    assignments->clear();
  if( explained_frac )
    *explained_frac = 0.0;

  // Sanity gates.
  const double e_top = offset + gain * static_cast<double>(nchannels);
  if( (gain <= 0.0) || (e_top < sm_min_spectrum_span) || (e_top > sm_max_spectrum_span) )
    return sm_gate_fail_score;
  if( fabs(offset) > std::min( 100.0, 0.1*gain*static_cast<double>(nchannels) ) )
    return sm_gate_fail_score;

  if( peaks.empty() || library.lines.empty() )
    return 0.0;

  // Weights: sqrt(area), normalized so the significant peaks sum to 1.
  double sig_weight_sum = 0.0;
  for( const CandidatePeak &p : peaks )
    if( p.significant )
      sig_weight_sum += sqrt( std::max( p.area, 1.0 ) );
  if( sig_weight_sum <= 0.0 )
    return 0.0;

  const vector<LineSource> &lines = library.lines;

  // Lines per keV around an energy (for the accidental-match discount).
  auto line_density = [&lines]( const double energy ) -> double {
    const double window = 50.0;
    const LineSource lo_val{}, hi_val{};
    auto lo = std::lower_bound( begin(lines), end(lines), energy - window,
        []( const LineSource &l, const double e ){ return l.energy < e; } );
    auto hi = std::upper_bound( begin(lines), end(lines), energy + window,
        []( const double e, const LineSource &l ){ return e < l.energy; } );
    (void)lo_val; (void)hi_val;
    return static_cast<double>( std::distance( lo, hi ) ) / (2.0*window);
  };

  double score = 0.0;
  double explained_sig_weight = 0.0;

  struct MatchInfo
  {
    bool matched = false;
    size_t line_index = 0;
    double delta = 0.0;
    double quality = 0.0;   //exp( -delta^2 / ... ), in (0,1]
  };
  vector<MatchInfo> matches( peaks.size() );

  // --- Pass 1: match each candidate peak to its best library line -------------------------
  for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
  {
    const CandidatePeak &peak = peaks[peak_index];
    const double e_pred = offset + gain*peak.channel;
    if( e_pred < sm_min_line_energy )
      continue;

    const double tol = std::max( { 3.0*gain*peak.sigma_ch, 0.002*e_pred, 0.5 } );

    auto first = std::lower_bound( begin(lines), end(lines), e_pred - tol,
        []( const LineSource &l, const double e ){ return l.energy < e; } );

    double best_pick = 0.0;
    MatchInfo &match = matches[peak_index];
    for( auto it = first; (it != end(lines)) && (it->energy <= e_pred + tol); ++it )
    {
      const double delta = e_pred - it->energy;
      const double sigma_match = tol / 3.0;
      const double quality = exp( -0.5 * (delta/sigma_match) * (delta/sigma_match) );
      // Prefer stronger lines when several are in-tolerance, but let proximity dominate.
      const double pick = quality * ( 0.25 + 0.75*it->rel_intensity )
                          * ( it->user_supplied ? 1.3 : 1.0 );
      if( pick > best_pick )
      {
        best_pick = pick;
        match.matched = true;
        match.line_index = static_cast<size_t>( std::distance( begin(lines), it ) );
        match.delta = delta;
        match.quality = quality;
      }
    }//for( lines within tolerance )

    if( !match.matched )
      continue;

    const LineSource &line = lines[match.line_index];
    const double w = sqrt( std::max( peak.area, 1.0 ) ) / sig_weight_sum;

    // Accidental-match discount: a match in a line-dense region proves little.  The window
    //  used here is the peaks *position resolution* (not the acceptance tolerance, whose
    //  percent-of-energy slack would unfairly discount precise high-energy matches).
    const double match_window = std::max( 3.0*gain*peak.sigma_ch, 0.5 );
    const double density = std::max( line_density( e_pred ), 0.01 );
    const double rarity = std::max( 0.0, std::min( 3.0, -log( 2.0*match_window*density ) ) );

    double contribution = w * match.quality * rarity * ( line.user_supplied ? 1.3 : 1.0 );

    // FWHM plausibility.  A peak notably *narrower* than the detector resolution is
    //  physically impossible - it is the signature of a degenerate too-small gain squeezing
    //  the spectrum - so the narrow side is gated hard.  The wide side stays soft (summing,
    //  unresolved multiplets, and our FWHM curves are only nominal).
    const double obs_fwhm = 2.3548 * gain * peak.sigma_ch;
    const double exp_fwhm = expected_fwhm_kev( line.energy, is_high_res );
    const double fwhm_ratio = obs_fwhm / exp_fwhm;
    const double narrow_limit = is_high_res ? 0.25 : 0.2;

    double fwhm_factor = 1.0;
    if( fwhm_ratio < narrow_limit )
    {
      fwhm_factor = 0.0;   //treated as unmatched (and penalized) below
    }else if( fwhm_ratio < 2.0*narrow_limit )
    {
      fwhm_factor = 0.2;
    }else if( line.type == LineSource::Type::Annihilation )
    {
      if( is_high_res )
      {
        // 511 keV is Doppler broadened - the extra width is *evidence* for this hypothesis.
        if( (fwhm_ratio > 1.25) && (fwhm_ratio < 3.5) )
          score += 0.25 * w;
        else if( fwhm_ratio < 1.1 )
          score -= 0.15 * w;
      }
    }else if( fwhm_ratio > 2.5 )
    {
      fwhm_factor = 0.5;
    }

    if( fwhm_factor <= 0.0 )
    {
      match.matched = false;
      continue;
    }

    score += contribution * fwhm_factor;
    if( peak.significant )
      explained_sig_weight += w * match.quality * fwhm_factor;

    if( assignments )
      assignments->push_back( PeakAssignment{ peak_index, line, match.delta } );
  }//for( loop over candidate peaks - pass 1 )

  // --- Pass 2: unexplained significant peaks - escape/backscatter credit, else penalty ----
  vector<double> matched_energies;
  for( size_t i = 0; i < peaks.size(); ++i )
    if( matches[i].matched )
      matched_energies.push_back( lines[matches[i].line_index].energy );

  for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
  {
    const CandidatePeak &peak = peaks[peak_index];
    if( !peak.significant || matches[peak_index].matched )
      continue;

    const double e_pred = offset + gain*peak.channel;
    const double w = sqrt( std::max( peak.area, 1.0 ) ) / sig_weight_sum;
    const double tol = std::max( { 3.0*gain*peak.sigma_ch, 0.002*e_pred, 1.0 } );

    LineSource derived;
    bool explained = false;
    for( const double e_match : matched_energies )
    {
      if( is_high_res && (e_match > 1300.0) && (fabs(e_pred - (e_match - 510.99)) < tol) )
      {
        derived.type = LineSource::Type::SingleEscape;
        derived.energy = e_match - 510.99;
        explained = true;
      }else if( is_high_res && (e_match > 1600.0) && (fabs(e_pred - (e_match - 1022.0)) < tol) )
      {
        derived.type = LineSource::Type::DoubleEscape;
        derived.energy = e_match - 1022.0;
        explained = true;
      }else if( e_match > 250.0 )
      {
        // Backscatter peaks are broad and slightly asymmetric - use a wide tolerance.
        const double e_bs = e_match / ( 1.0 + 2.0*e_match/510.99 );
        if( fabs(e_pred - e_bs) < std::max( 2.0*tol, 10.0 ) )
        {
          derived.type = LineSource::Type::Backscatter;
          derived.energy = e_bs;
          explained = true;
        }
      }

      if( explained )
        break;
    }//for( matched primary energies )

    if( !explained && (fabs(e_pred - sm_annihilation_energy) < tol) )
    {
      derived.type = LineSource::Type::Annihilation;
      derived.energy = sm_annihilation_energy;
      explained = true;
    }

    if( explained )
    {
      score += 0.1 * w;
      explained_sig_weight += w;
      if( assignments )
        assignments->push_back( PeakAssignment{ peak_index, derived, e_pred - derived.energy } );
    }else
    {
      score -= 0.5 * w;
    }
  }//for( loop over candidate peaks - pass 2 )

  // --- Nuclide self-consistency: matched parents should show their other strong lines -----
  map<const SandiaDecay::Nuclide *, vector<size_t>> parent_matches; //-> matched line indices
  map<const SandiaDecay::Nuclide *, double> parent_quality_sum;
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    if( matches[i].matched && lines[matches[i].line_index].nuc )
    {
      parent_matches[ lines[matches[i].line_index].nuc ].push_back( matches[i].line_index );
      parent_quality_sum[ lines[matches[i].line_index].nuc ] += matches[i].quality;
    }
  }

  for( const auto &pm : parent_matches )
  {
    const SandiaDecay::Nuclide * const nuc = pm.first;
    const vector<size_t> &matched_lines = pm.second;

    float max_matched_intensity = 0.0f;
    for( const size_t index : matched_lines )
      max_matched_intensity = std::max( max_matched_intensity, lines[index].rel_intensity );

    // Lines of this parent we should have seen: at least ~25% as strong as its strongest
    //  matched line, and inside the implied energy range of the spectrum.
    size_t n_expected = 0, n_matched_expected = 0;
    for( size_t i = 0; i < lines.size(); ++i )
    {
      const LineSource &line = lines[i];
      if( (line.nuc != nuc) || (line.rel_intensity < 0.25f*max_matched_intensity) )
        continue;
      if( (line.energy < std::max(sm_min_line_energy, offset + 0.01*gain*nchannels))
          || (line.energy > 0.98*e_top) )
        continue;
      ++n_expected;
      if( std::find( begin(matched_lines), end(matched_lines), i ) != end(matched_lines) )
        ++n_matched_expected;
    }

    if( !n_expected )
      continue;

    const double fraction = static_cast<double>(n_matched_expected)
                            / static_cast<double>(n_expected);
    const double n_extra = static_cast<double>( std::min( matched_lines.size(), size_t(5) ) - 1 );
    // Scaled by the average match quality, so a parent "matched" by several sloppy
    //  few-keV-off accidentals does not collect the same bonus as clean matches.
    const double avg_quality = parent_quality_sum[nuc] / static_cast<double>(matched_lines.size());
    score += 0.35 * n_extra * fraction * avg_quality;
    if( (n_expected >= 3) && (fraction < 0.5) )
      score -= 0.4 * (1.0 - fraction);
  }//for( matched parents )

  // --- Amplitude plausibility (weak): within a parent, stronger lines should tend to give
  //     bigger peaks.  Spearman rank correlation between peak area and line intensity.
  for( const auto &pm : parent_matches )
  {
    const vector<size_t> &matched_lines = pm.second;
    if( matched_lines.size() < 3 )
      continue;

    // Collect (area, intensity) for this parents matched peaks.
    vector<pair<double,double>> pts;
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      if( !matches[i].matched || (lines[matches[i].line_index].nuc != pm.first) )
        continue;
      pts.emplace_back( peaks[i].area, lines[matches[i].line_index].rel_intensity );
    }
    if( pts.size() < 3 )
      continue;

    auto ranks = []( const vector<pair<double,double>> &v, const bool first ) -> vector<double> {
      vector<size_t> order( v.size() );
      for( size_t i = 0; i < v.size(); ++i )
        order[i] = i;
      std::sort( begin(order), end(order), [&v, first]( const size_t a, const size_t b ){
        return first ? (v[a].first < v[b].first) : (v[a].second < v[b].second);
      } );
      vector<double> r( v.size() );
      for( size_t i = 0; i < order.size(); ++i )
        r[order[i]] = static_cast<double>( i );
      return r;
    };

    const vector<double> ra = ranks( pts, true ), rb = ranks( pts, false );
    double d2 = 0.0;
    for( size_t i = 0; i < ra.size(); ++i )
      d2 += (ra[i] - rb[i]) * (ra[i] - rb[i]);
    const double n = static_cast<double>( ra.size() );
    const double spearman = 1.0 - 6.0*d2 / ( n*(n*n - 1.0) );
    if( spearman > 0.0 )
      score += 0.1 * spearman;
  }//for( matched parents - amplitude plausibility )

  // Completeness: explaining *all* the significant peaks is worth much more than a partial
  //  (possibly accidental) explanation - the quadratic keeps partial credit small.
  const double explained = std::min( explained_sig_weight, 1.0 );
  score += 1.2 * explained * explained;

  if( explained_frac )
    *explained_frac = explained;

  return score;
}//score_hypothesis(...)


std::vector<GuessedCal>
guess_energy_cal( const std::shared_ptr<const SpecUtils::Measurement> &meas,
                  const GuessOptions &options )
{
  vector<GuessedCal> answer;

  try
  {
    bool is_high_res = true;
    const vector<CandidatePeak> peaks = find_candidate_peaks_channelspace( meas, is_high_res );
    if( peaks.empty() || is_cancelled(options.cancel) )
      return answer;

    const size_t nchan = meas->num_gamma_channels();

    const LineLibrary library = build_line_library( options.user_nuclides );
    if( library.lines.empty() || is_cancelled(options.cancel) )
      return answer;

    // The "right edge" of real data, for the endpoint strategy (trailing zero channels would
    //  otherwise make endpoint hypotheses meaningless).
    size_t extent_lower = 0, extent_upper = nchan - 1;
    {
      auto identity_cal = make_shared<EnergyCalibration>();
      identity_cal->set_default_polynomial( nchan, {0.0f, 1.0f}, {} );
      auto chanmeas = make_shared<Measurement>( *meas );
      chanmeas->set_energy_calibration( identity_cal );
      if( !ExperimentalPeakSearch::find_spectroscopic_extent( chanmeas, extent_lower, extent_upper )
          || (extent_upper <= extent_lower) )
      {
        extent_lower = 0;
        extent_upper = nchan - 1;
      }
    }

    // ---- Hypothesis seeds --------------------------------------------------------------
    struct Seed
    {
      double offset, gain;
      SeedStrategy strategy;
      double strategy_energy1, strategy_energy2;
      string provenance;
    };
    vector<Seed> seeds;

    vector<const CandidatePeak *> sig_peaks;
    for( const CandidatePeak &p : peaks )
      if( p.significant )
        sig_peaks.push_back( &p );

    // G1: pairs of significant peaks whose channel ratio matches a known line-energy ratio
    //  (gain-invariant when offset ~ 0).
    for( size_t i = 0; i < sig_peaks.size(); ++i )
    {
      for( size_t j = i + 1; j < sig_peaks.size(); ++j )
      {
        const CandidatePeak *low = sig_peaks[i], *high = sig_peaks[j];
        if( low->channel > high->channel )
          std::swap( low, high );

        const double ratio = low->channel / high->channel;
        const double ratio_tol = ratio * ( low->sigma_ch/low->channel
                                           + high->sigma_ch/high->channel + 0.005 );

        auto first = std::lower_bound( begin(library.ratios), end(library.ratios),
            static_cast<float>(ratio - ratio_tol),
            []( const LineLibrary::RatioEntry &e, const float r ){ return e.ratio < r; } );

        for( auto it = first;
             (it != end(library.ratios)) && (it->ratio <= ratio + ratio_tol); ++it )
        {
          const LineSource &line_low = library.lines[it->low_index];
          const LineSource &line_high = library.lines[it->high_index];
          Seed seed;
          seed.offset = 0.0;
          seed.gain = line_high.energy / high->channel;
          seed.strategy = SeedStrategy::PeakPair;
          seed.strategy_energy1 = line_low.energy;
          seed.strategy_energy2 = line_high.energy;
          seed.provenance = "peak pair " + energy_str(line_low.energy) + "/"
                            + energy_str(line_high.energy) + " keV";
          seeds.push_back( seed );
        }
      }//for( j )

      if( is_cancelled(options.cancel) )
        return answer;
    }//for( i - G1 )

    // G2: single-peak anchors.
    if( !sig_peaks.empty() )
    {
      const CandidatePeak *highest = sig_peaks.front(), *biggest = sig_peaks.front();
      for( const CandidatePeak *p : sig_peaks )
      {
        if( p->channel > highest->channel )
          highest = p;
        if( p->area > biggest->area )
          biggest = p;
      }

      for( const double e : { 2614.5, 1460.8 } )
        seeds.push_back( Seed{ 0.0, e/highest->channel, SeedStrategy::SingleAnchor, e, 0.0,
                               "highest peak as " + energy_str(e) + " keV" } );
      for( const double e : { 1460.8, 661.7, 185.7 } )
        seeds.push_back( Seed{ 0.0, e/biggest->channel, SeedStrategy::SingleAnchor, e, 0.0,
                               "largest peak as " + energy_str(e) + " keV" } );
      for( const CandidatePeak *p : sig_peaks )
        if( p->width_outlier )
          seeds.push_back( Seed{ 0.0, sm_annihilation_energy/p->channel, SeedStrategy::Wide511,
                                 sm_annihilation_energy, 0.0, "wide peak as 511 keV" } );
    }//if( have significant peaks )

    // G3: common endpoint energies at the right edge of the data.
    for( const double e : { 400.0, 800.0, 1500.0, 2000.0, 3000.0, 8000.0, 10000.0, 13000.0 } )
      seeds.push_back( Seed{ 0.0, e/static_cast<double>(extent_upper), SeedStrategy::Endpoint,
                             e, 0.0, "spectrum ends at " + energy_str(e) + " keV" } );

    // G4: common gains.
    for( const double g : { 0.075, 0.1, 0.125, 0.25, 0.3, 0.5, 0.75, 1.0, 1.5, 3.0 } )
      seeds.push_back( Seed{ 0.0, g, SeedStrategy::CommonGain, g, 0.0,
                             "common gain " + energy_str(g) + " keV/channel" } );

    // G5: coarse full-range gain scan (0.5% steps) as a backstop for spectra where only one
    //  genuine line exists among the candidates.
    for( double e_top = sm_min_spectrum_span; e_top <= sm_max_spectrum_span; e_top *= 1.005 )
      seeds.push_back( Seed{ 0.0, e_top/static_cast<double>(nchan), SeedStrategy::GainScan,
                             0.0, 0.0, "gain scan" } );

    // ---- Score all seeds ---------------------------------------------------------------
    struct Scored
    {
      double offset, gain, score, explained_frac;
      SeedStrategy strategy;
      double strategy_energy1, strategy_energy2;
      string provenance;
      vector<PeakAssignment> assignments;
    };
    vector<Scored> scored;
    scored.reserve( seeds.size() );

    for( size_t i = 0; i < seeds.size(); ++i )
    {
      if( ((i % 64) == 0) && is_cancelled(options.cancel) )
        return answer;

      const Seed &seed = seeds[i];
      Scored s;
      s.offset = seed.offset;
      s.gain = seed.gain;
      s.strategy = seed.strategy;
      s.strategy_energy1 = seed.strategy_energy1;
      s.strategy_energy2 = seed.strategy_energy2;
      s.provenance = seed.provenance;
      s.score = score_hypothesis( seed.offset, seed.gain, peaks, nchan, library,
                                  is_high_res, &s.assignments, &s.explained_frac );
      if( s.score > 0.0 )
        scored.push_back( std::move(s) );
    }//for( seeds )

    std::sort( begin(scored), end(scored),
        []( const Scored &a, const Scored &b ){ return a.score > b.score; } );

    // Greedy de-dup: two hypotheses are the same if their predictions differ by less than
    //  about one peak width at the top of the spectrum.
    auto is_duplicate = [nchan, is_high_res]( const Scored &a, const Scored &b ) -> bool {
      const double e_top = a.offset + a.gain*static_cast<double>(nchan);
      const double sep = fabs( a.gain - b.gain )*static_cast<double>(nchan)
                         + fabs( a.offset - b.offset );
      return sep < std::max( 2.0, expected_fwhm_kev( e_top, is_high_res ) );
    };

    vector<Scored> survivors;
    for( Scored &s : scored )
    {
      bool dup = false;
      for( const Scored &kept : survivors )
        dup = ( dup || is_duplicate( kept, s ) );
      if( !dup )
        survivors.push_back( std::move(s) );
      if( survivors.size() >= 20 )
        break;
    }

    // ---- Refine the survivors: a local gain polish (so each seed snaps to the maximum of
    //      its score basin, and near-duplicate seeds converge), then a linear (offset,gain)
    //      LLS fit to the assignments, re-matched / re-fit once more.
    for( Scored &s : survivors )
    {
      if( is_cancelled(options.cancel) )
        return answer;

      for( int step = -20; step <= 20; ++step )
      {
        if( !step )
          continue;
        const double gain = s.gain * ( 1.0 + 0.0005*step );
        vector<PeakAssignment> trial_assignments;
        double trial_explained = 0.0;
        const double trial = score_hypothesis( s.offset, gain, peaks, nchan, library,
                                    is_high_res, &trial_assignments, &trial_explained );
        if( trial > s.score )
        {
          s.gain = gain;
          s.score = trial;
          s.explained_frac = trial_explained;
          s.assignments = std::move( trial_assignments );
        }
      }//for( local gain polish )

      for( int iteration = 0; iteration < 2; ++iteration )
      {
        vector<EnergyCal::RecalPeakInfo> infos;
        for( const PeakAssignment &assign : s.assignments )
        {
          if( (assign.line.type != LineSource::Type::NuclideGamma)
              && (assign.line.type != LineSource::Type::Annihilation)
              && (assign.line.type != LineSource::Type::FluorXray) )
            continue;

          // Only well-placed matches may anchor the fit - a marginal (few-keV-off, likely
          //  accidental) match used as an anchor drags the whole calibration to it.
          const CandidatePeak &peak = peaks[assign.peak_index];
          const double e_pred = s.offset + s.gain*peak.channel;
          const double tol = std::max( { 3.0*s.gain*peak.sigma_ch, 0.002*e_pred, 0.5 } );
          if( fabs(assign.delta_kev) > 0.6*tol )
            continue;

          EnergyCal::RecalPeakInfo info;
          info.peakMean = e_pred;
          info.peakMeanUncert = std::max( s.gain*peak.sigma_ch, 0.01 );
          info.peakMeanBinNumber = peak.channel;
          info.photopeakEnergy = assign.line.energy;
          infos.push_back( info );
        }

        if( infos.size() < 2 )
          break;

        try
        {
          // Freeing the offset with only two anchors lets the fit *exactly* interpolate any
          //  two (possibly wrong) lines, erasing the evidence the deltas carry - so the
          //  offset stays fixed at zero until there are three or more anchors.
          const bool fit_offset = ( infos.size() >= 3 );
          vector<float> coefs{ 0.0f, static_cast<float>(s.gain) }, coefs_uncert;
          EnergyCal::fit_energy_cal_poly( infos, {fit_offset, true}, nchan, {},
                                          coefs, coefs_uncert );
          if( (coefs.size() < 2) || (coefs[1] <= 0.0f) )
            break;

          vector<PeakAssignment> new_assignments;
          double new_explained = 0.0;
          const double new_score = score_hypothesis( coefs[0], coefs[1], peaks, nchan, library,
                                          is_high_res, &new_assignments, &new_explained );
          if( new_score <= s.score )
            break;

          s.offset = coefs[0];
          s.gain = coefs[1];
          s.score = new_score;
          s.explained_frac = new_explained;
          s.assignments = std::move( new_assignments );
        }catch( std::exception & )
        {
          break;
        }
      }//for( refinement iterations )
    }//for( survivors )

    std::sort( begin(survivors), end(survivors),
        []( const Scored &a, const Scored &b ){ return a.score > b.score; } );

    // ---- Final de-dup, truncate, and package -------------------------------------------
    vector<Scored> final_list;
    for( Scored &s : survivors )
    {
      bool dup = false;
      for( const Scored &kept : final_list )
        dup = ( dup || is_duplicate( kept, s ) );
      if( !dup )
        final_list.push_back( std::move(s) );
      if( final_list.size() >= options.max_results )
        break;
    }

    for( Scored &s : final_list )
    {
      GuessedCal result;
      result.offset = s.offset;
      result.gain = s.gain;
      result.score = s.score;
      result.explained_frac = s.explained_frac;
      result.strategy = s.strategy;
      result.strategy_energy1 = s.strategy_energy1;
      result.strategy_energy2 = s.strategy_energy2;
      result.provenance = s.provenance;
      result.assignments = std::move( s.assignments );

      // Implied parents, ordered by their (weight-summed) contribution.  A 511 keV match
      //  alone does not implicate a specific beta+ parent, so annihilation lines dont count.
      map<const SandiaDecay::Nuclide *, double> parent_weight;
      for( const PeakAssignment &assign : result.assignments )
        if( assign.line.nuc && (assign.line.type != LineSource::Type::Annihilation) )
          parent_weight[assign.line.nuc] += sqrt( std::max( peaks[assign.peak_index].area, 1.0 ) );

      vector<pair<const SandiaDecay::Nuclide *, double>> ordered( begin(parent_weight),
                                                                  end(parent_weight) );
      std::sort( begin(ordered), end(ordered),
          []( const pair<const SandiaDecay::Nuclide *, double> &a,
              const pair<const SandiaDecay::Nuclide *, double> &b ){ return a.second > b.second; } );
      for( const pair<const SandiaDecay::Nuclide *, double> &p : ordered )
        result.nuclides.push_back( p.first );

      answer.push_back( std::move(result) );
    }//for( final results )
  }catch( std::exception & )
  {
    answer.clear();
  }

  return answer;
}//guess_energy_cal(...)

}//namespace GainGuessCalc


using SpecUtils::SpectrumType;
using GainGuessCalc::GuessedCal;
using GainGuessCalc::SeedStrategy;

using namespace Wt;


namespace
{
  string format_number( const char * const fmt, const double value )
  {
    char buffer[48];
    snprintf( buffer, sizeof(buffer), fmt, value );
    return buffer;
  }
}//namespace


EnergyCalGainGuess::EnergyCalGainGuess( std::shared_ptr<std::vector<MeasToApplyCoefChangeTo>> measToChange,
                                        EnergyCalTool *cal, AuxWindow *parent )
  : WContainerWidget(),
    m_interspec( InterSpec::instance() ),
    m_calibrator( cal ),
    m_parent( parent ),
    m_measToChange( measToChange ),
    m_foreground( nullptr ),
    m_userNuclides{},
    m_results{},
    m_selectedRow( -1 ),
    m_generation( 0 ),
    m_calcRunning( false ),
    m_cancelCalc( nullptr ),
    m_chart( nullptr ),
    m_nuclideEdit( nullptr ),
    m_nuclideSuggest( nullptr ),
    m_nuclideChips( nullptr ),
    m_resultTable( nullptr ),
    m_resultRows{},
    m_status( nullptr ),
    m_cancel( nullptr ),
    m_use( nullptr )
{
  assert( m_interspec && cal && parent );

  wApp->useStyleSheet( "InterSpec_resources/EnergyCalGainGuess.css" );
  m_interspec->useMessageResourceBundle( "EnergyCalGainGuess" );

  addStyleClass( "EnergyCalGainGuess" );

  m_foreground = m_interspec->displayedHistogram( SpectrumType::Foreground );

  // Chart on top (InterSpec tools put the spectrum first), previewing the selected candidate.
  m_chart = addNew<D3SpectrumDisplayDiv>();
  m_chart->addStyleClass( "GainGuessChart" );
  m_chart->setCompactAxis( true );
  m_chart->setYAxisLog( true );
  m_chart->disableLegend();
  m_chart->showYAxisScalers( false );

  WContainerWidget * const controls = addNew<WContainerWidget>();
  controls->addStyleClass( "GainGuessControls" );

  // Candidate-nuclide hint row.
  WContainerWidget * const hintRow = controls->addNew<WContainerWidget>();
  hintRow->addStyleClass( "GainGuessRow" );

  hintRow->addNew<WLabel>( WString::tr("ecgg-nuclide-hint-label") );
  m_nuclideEdit = hintRow->addNew<WLineEdit>();
  m_nuclideEdit->addStyleClass( "GainGuessNuclideEdit" );
  m_nuclideEdit->setAutoComplete( false );

  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  std::shared_ptr<IsotopeNameFilterModel> isoSuggestModel
                                        = std::make_shared<IsotopeNameFilterModel>();
  m_nuclideSuggest = addNew<WSuggestionPopup>( matcherJs, replacerJs );
  m_nuclideSuggest->addStyleClass( "nuclide-suggest" );
  IsotopeNameFilterModel::setQuickTypeFixHackjs( m_nuclideSuggest );
  IsotopeNameFilterModel::setEnterKeyMatchFixJs( m_nuclideSuggest, m_nuclideEdit );
  isoSuggestModel->filter( "" );
  m_nuclideSuggest->setFilterLength( -1 );
  m_nuclideSuggest->setModel( isoSuggestModel );
  m_nuclideSuggest->filterModel().connect( isoSuggestModel.get(), &IsotopeNameFilterModel::filter );
  m_nuclideSuggest->forEdit( m_nuclideEdit, PopupTrigger::Editing );

  WPushButton * const addBtn = hintRow->addNew<WPushButton>( WString::tr("ecgg-add") );
  addBtn->clicked().connect( this, &EnergyCalGainGuess::handleAddNuclide );
  m_nuclideEdit->enterPressed().connect( this, &EnergyCalGainGuess::handleAddNuclide );

  m_nuclideChips = hintRow->addNew<WContainerWidget>();
  m_nuclideChips->addStyleClass( "GainGuessChips" );

  // Ranked candidate table, inside a fixed-height scroll container.
  WContainerWidget * const tableWrap = controls->addNew<WContainerWidget>();
  tableWrap->addStyleClass( "GainGuessTableWrap" );
  m_resultTable = tableWrap->addNew<WTable>();
  m_resultTable->addStyleClass( "GainGuessTable" );

  m_status = controls->addNew<WText>();
  m_status->addStyleClass( "GainGuessStatus" );

  // Footer buttons.
  WContainerWidget * const footer = parent->footer();
  AuxWindow::addHelpInFooter( footer, "energy-cal-guess-gain" );
  m_cancel = footer->addNew<WPushButton>( WString::tr("Cancel") );
  m_cancel->clicked().connect( this, [this](){ handleFinish( Wt::DialogCode::Rejected ); } );
  m_use = footer->addNew<WPushButton>( WString::tr("Use") );
  m_use->clicked().connect( this, [this](){ handleFinish( Wt::DialogCode::Accepted ); } );
  m_use->disable();

  startCalc();

  parent->resizeScaledWindow( 0.85, 0.85 );
  parent->centerWindow();
}//EnergyCalGainGuess constructor


EnergyCalGainGuess::~EnergyCalGainGuess()
{
  if( m_cancelCalc )
    m_cancelCalc->store( true );
}


void EnergyCalGainGuess::startCalc()
{
  if( m_cancelCalc )
    m_cancelCalc->store( true );
  m_cancelCalc = std::make_shared<std::atomic<bool>>( false );

  const int generation = ++m_generation;
  m_calcRunning = true;
  m_results.clear();
  m_selectedRow = -1;
  m_use->disable();
  rebuildResultTable();
  updatePreview();
  m_status->removeStyleClass( "GainGuessLowConfidence" );
  m_status->setText( WString::tr("ecgg-searching") );

  GainGuessCalc::GuessOptions options;
  options.user_nuclides = m_userNuclides;
  options.cancel = m_cancelCalc;

  const shared_ptr<const Measurement> meas = m_foreground;
  const string widgetId = id();
  const string sessionId = wApp ? wApp->sessionId() : string();

  Wt::WServer * const server = Wt::WServer::instance();
  if( !server || sessionId.empty() )
  {
    // No server (e.g. unusual test context) - fall back to a synchronous computation.
    onCalcDone( generation, GainGuessCalc::guess_energy_cal( meas, options ) );
    return;
  }

  server->ioService().boost::asio::io_service::post( [=](){
    // --- worker thread: the CPU-heavy search on the (immutable) spectrum snapshot ---
    const vector<GuessedCal> results = GainGuessCalc::guess_energy_cal( meas, options );

    // --- back on the session thread: only the inert widget id crosses the boundary; the
    //     widget-tree lookup happens here, under the applications UpdateLock (see CLAUDE.md) ---
    server->post( sessionId, [widgetId, generation, results](){
      WApplication * const app = WApplication::instance();
      if( !app || !app->domRoot() )
        return;
      WWidget * const w = app->domRoot()->findById( widgetId );
      EnergyCalGainGuess * const self = dynamic_cast<EnergyCalGainGuess *>( w );
      if( self )
      {
        self->onCalcDone( generation, results );
        app->triggerUpdate();
      }
    } );
  } );
}//startCalc()


void EnergyCalGainGuess::onCalcDone( int generation, const std::vector<GuessedCal> &results )
{
  if( generation != m_generation )
    return;  //a newer computation (nuclide hint change) superseded this result

  m_calcRunning = false;
  m_results = results;
  rebuildResultTable();

  if( m_results.empty() )
  {
    m_status->setText( WString::tr("ecgg-no-results") );
    updatePreview();
    return;
  }

  if( m_results.front().explained_frac < GainGuessCalc::sm_min_confident_explained_frac )
  {
    m_status->addStyleClass( "GainGuessLowConfidence" );
    m_status->setText( WString::tr("ecgg-low-confidence") );
  }else
  {
    m_status->setText( WString::tr("ecgg-num-results") );
  }

  handleRowSelected( 0 );
}//onCalcDone(...)


void EnergyCalGainGuess::handleAddNuclide()
{
  const string txt = m_nuclideEdit->text().toUTF8();
  if( txt.empty() )
    return;

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide * const nuc = db ? db->nuclide( txt ) : nullptr;
  if( !nuc )
  {
    m_interspec->logMessage( WString::tr("ecgg-unknown-nuclide").arg(txt), 2 );
    return;
  }

  m_nuclideEdit->setText( "" );

  if( std::find( begin(m_userNuclides), end(m_userNuclides), nuc ) != end(m_userNuclides) )
    return;

  m_userNuclides.push_back( nuc );

  WContainerWidget * const chip = m_nuclideChips->addNew<WContainerWidget>();
  chip->addStyleClass( "GainGuessChip" );
  chip->addNew<WText>( WString::fromUTF8(nuc->symbol) );
  WText * const close = chip->addNew<WText>( WString::fromUTF8("×") );
  close->addStyleClass( "GainGuessChipClose" );
  close->clicked().connect( this, [this, nuc, chip](){
    const std::unique_ptr<Wt::WWidget> doomed = chip->removeFromParent();
    removeNuclide( nuc );
  } );

  startCalc();
}//handleAddNuclide()


void EnergyCalGainGuess::removeNuclide( const SandiaDecay::Nuclide *nuc )
{
  const auto pos = std::find( begin(m_userNuclides), end(m_userNuclides), nuc );
  if( pos != end(m_userNuclides) )
    m_userNuclides.erase( pos );
  startCalc();
}//removeNuclide(...)


void EnergyCalGainGuess::rebuildResultTable()
{
  m_resultTable->clear();
  m_resultRows.clear();

  if( m_results.empty() )
    return;

  // Header row.
  const char * const header_keys[] = { "ecgg-col-rank", "ecgg-col-gain", "ecgg-col-offset",
                                       "ecgg-col-explained", "ecgg-col-nuclides",
                                       "ecgg-col-strategy" };
  for( int col = 0; col < 6; ++col )
    m_resultTable->elementAt( 0, col )->addNew<WText>( WString::tr(header_keys[col]) );
  m_resultTable->rowAt( 0 )->addStyleClass( "GainGuessHeaderRow" );

  for( size_t i = 0; i < m_results.size(); ++i )
  {
    const GuessedCal &result = m_results[i];
    const int table_row = static_cast<int>( i ) + 1;

    string nucs;
    for( size_t j = 0; (j < result.nuclides.size()) && (j < 4); ++j )
      nucs += (j ? ", " : "") + result.nuclides[j]->symbol;
    if( result.nuclides.size() > 4 )
      nucs += ", …";

    WString strategy;
    switch( result.strategy )
    {
      case SeedStrategy::PeakPair:
        strategy = WString::tr("ecgg-strategy-pair")
                     .arg( format_number("%.1f", result.strategy_energy1) )
                     .arg( format_number("%.1f", result.strategy_energy2) );
        break;
      case SeedStrategy::SingleAnchor:
        strategy = WString::tr("ecgg-strategy-anchor")
                     .arg( format_number("%.1f", result.strategy_energy1) );
        break;
      case SeedStrategy::Wide511:
        strategy = WString::tr("ecgg-strategy-wide511");
        break;
      case SeedStrategy::Endpoint:
        strategy = WString::tr("ecgg-strategy-endpoint")
                     .arg( format_number("%.0f", result.strategy_energy1) );
        break;
      case SeedStrategy::CommonGain:
        strategy = WString::tr("ecgg-strategy-common-gain");
        break;
      case SeedStrategy::GainScan:
        strategy = WString::tr("ecgg-strategy-scan");
        break;
    }//switch( result.strategy )

    const WString cell_text[6] = {
      WString::fromUTF8( std::to_string(i + 1) ),
      WString::fromUTF8( format_number("%.5g", result.gain) ),
      WString::fromUTF8( format_number("%.2f", result.offset) ),
      WString::fromUTF8( format_number("%.0f%%", 100.0*result.explained_frac) ),
      WString::fromUTF8( nucs ),
      strategy
    };

    const int row_index = static_cast<int>( i );
    for( int col = 0; col < 6; ++col )
    {
      WTableCell * const cell = m_resultTable->elementAt( table_row, col );
      cell->addNew<WText>( cell_text[col] );
      if( col == 4 )
        cell->addStyleClass( "GainGuessNucCell" );
      cell->clicked().connect( this, [this, row_index](){ handleRowSelected( row_index ); } );
    }//for( columns )

    WTableRow * const row = m_resultTable->rowAt( table_row );
    row->addStyleClass( "GainGuessResultRow" );
    m_resultRows.push_back( row );
  }//for( loop over results )
}//rebuildResultTable()


void EnergyCalGainGuess::handleRowSelected( const int row_index )
{
  if( (row_index < 0) || (row_index >= static_cast<int>(m_results.size())) )
    return;

  m_selectedRow = row_index;
  for( size_t i = 0; i < m_resultRows.size(); ++i )
  {
    if( static_cast<int>(i) == row_index )
      m_resultRows[i]->addStyleClass( "selected" );
    else
      m_resultRows[i]->removeStyleClass( "selected" );
  }

  m_use->enable();
  updatePreview();
}//handleRowSelected(...)


void EnergyCalGainGuess::updatePreview()
{
  m_chart->clearAllReferncePhotoPeakLines();

  if( !m_foreground || (m_selectedRow < 0)
      || (m_selectedRow >= static_cast<int>(m_results.size())) )
  {
    m_chart->setData( nullptr, false );
    return;
  }

  const GuessedCal &guess = m_results[m_selectedRow];

  try
  {
    const size_t nchan = m_foreground->num_gamma_channels();
    auto cal = std::make_shared<EnergyCalibration>();
    cal->set_polynomial( nchan, { static_cast<float>(guess.offset),
                                  static_cast<float>(guess.gain) }, {} );
    auto adjusted = std::make_shared<Measurement>( *m_foreground );
    adjusted->set_energy_calibration( cal );
    m_chart->setData( adjusted, false );
  }catch( std::exception & )
  {
    m_chart->setData( nullptr, false );
    return;
  }

  showPreviewRefLines( guess );
}//updatePreview()


void EnergyCalGainGuess::showPreviewRefLines( const GuessedCal &guess )
{
  // Draw the reference lines of the implied nuclides onto the preview chart, so the user can
  //  eyeball how well the candidate calibration lines the peaks up before committing.
  const vector<Wt::WColor> &colors = ReferencePhotopeakDisplay::sm_def_line_colors;

  size_t shown = 0;
  for( size_t i = 0; (i < guess.nuclides.size()) && (shown < 4); ++i )
  {
    const SandiaDecay::Nuclide * const nuc = guess.nuclides[i];
    if( !nuc )
      continue;

    RefLineInput input;
    input.m_input_txt = nuc->symbol;
    input.m_color = colors.empty() ? Wt::WColor("#0000FF") : colors[shown % colors.size()];
    input.m_showGammas = true;
    input.m_showXrays = true;

    try
    {
      const shared_ptr<ReferenceLineInfo> ref = ReferenceLineInfo::generateRefLineInfo( input );
      if( !ref || (ref->m_validity != ReferenceLineInfo::InputValidity::Valid) )
        continue;

      m_chart->setReferncePhotoPeakLines( *ref );
      m_chart->persistCurrentReferncePhotoPeakLines();
      ++shown;
    }catch( std::exception & )
    {
    }
  }//for( implied nuclides )
}//showPreviewRefLines(...)


void EnergyCalGainGuess::applyChanges()
{
  if( (m_selectedRow < 0) || (m_selectedRow >= static_cast<int>(m_results.size()))
      || !m_foreground || !m_measToChange )
    return;

  const GuessedCal &guess = m_results[m_selectedRow];

  try
  {
    const size_t nchan = m_foreground->num_gamma_channels();
    auto newcal = std::make_shared<EnergyCalibration>();
    newcal->set_polynomial( nchan, { static_cast<float>(guess.offset),
                                     static_cast<float>(guess.gain) }, {} );

    const shared_ptr<const EnergyCalibration> orig = m_foreground->energy_calibration();

    // applyCalChange() translates peaks, refreshes the GUI and charts, and registers the
    //  undo/redo step (through its EnergyCalUndoRedoSentry).
    m_calibrator->applyCalChange( orig, newcal, *m_measToChange, false );

    // Turn on reference lines for the implied nuclides (primary + one more), so the user can
    //  immediately eyeball how well the spectrum lines up.
    ReferencePhotopeakDisplay * const reflines = m_interspec->referenceLinesWidget();
    if( reflines && !guess.nuclides.empty() )
    {
      reflines->setIsotope( guess.nuclides[0] );
      if( guess.nuclides.size() > 1 )
      {
        reflines->persistCurentLines();
        reflines->setIsotope( guess.nuclides[1] );
      }

      m_interspec->logMessage( WString::tr("ecgg-applied-toast")
                                 .arg( guess.nuclides[0]->symbol ), 1 );
    }//if( have implied nuclides )
  }catch( std::exception &e )
  {
    m_interspec->logMessage( WString::tr("ecgg-apply-fail").arg( string(e.what()) ), 2 );
  }
}//applyChanges()


void EnergyCalGainGuess::handleFinish( Wt::DialogCode result )
{
  switch( result )
  {
    case Wt::DialogCode::Rejected:
    {
      UndoRedoManager * const undoManager = m_interspec->undoRedoManager();
      if( m_parent && undoManager )
      {
        auto undo = [](){
          InterSpec * const viewer = InterSpec::instance();
          EnergyCalTool * const tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->moreActionBtnClicked( MoreActionsIndex::GuessGain );
        };

        auto redo = [](){
          InterSpec * const viewer = InterSpec::instance();
          EnergyCalTool * const tool = viewer ? viewer->energyCalTool() : nullptr;
          if( tool )
            tool->cancelMoreActionWindow();
        };

        undoManager->addUndoRedoStep( undo, redo, "Cancel energy cal gain guess" );
      }//if( m_parent && undoManager )
      break;
    }//case Rejected

    case Wt::DialogCode::Accepted:
      applyChanges();
      break;
  }//switch( result )

  if( m_parent )
    m_parent->hide();
}//handleFinish(...)
