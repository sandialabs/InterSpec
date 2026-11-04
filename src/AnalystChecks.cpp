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
#include <deque>
#include <future>
#include <string>
#include <limits>
#include <algorithm>

#include "Eigen/Dense"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/AnalystChecks.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/PeakSearchGuiUtils.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/RelActCalcAuto.h"

#include "SpecUtils/SpecUtilsAsync.h"

#include "SandiaDecay/SandiaDecay.h"

#include "cascade/CascadeTypes.h"
#include "cascade/SandiaDecayCascade.h"

#include <Wt/WServer.h>
#include <Wt/WIOService.h>
#include <Wt/WApplication.h>

using namespace std;

namespace
{
  /** Decay series a NORM line belongs to, used for the chain-consistency (secular-equilibrium) rescue
   of borderline peaks in detected_peaks. */
  enum class NormSeries : int { None = 0, K40, Th232, U238Ra, Annih };

  struct NormLine { double energy_kev; NormSeries series; };

  /** Well-known NORM / ubiquitous environmental-background gamma-ray lines (keV), tagged by decay
   series.

   A foreground peak matching one of these is held to a stricter elevation bar before being reported
   as "elevated over background" (see detected_peaks), since these lines are expected to be present -
   and to fluctuate - in any background, so a small foreground excess is usually noise rather than a
   real foreground source.  Covers K-40; the Th-232 chain (Ac-228, Pb-212, Bi-212, Tl-208); the
   U-238/Ra-226 chain (Th-234, Pa-234m, Ra-226, Pb-214, Bi-214); and the 511 keV annihilation line.
   */
  const NormLine sk_norm_lines[] = {
    // U-238 / Ra-226 chain
    {  63.29, NormSeries::U238Ra }, {  92.58, NormSeries::U238Ra }, { 186.21, NormSeries::U238Ra },
    { 240.99, NormSeries::U238Ra }, { 242.00, NormSeries::U238Ra }, { 295.22, NormSeries::U238Ra },
    { 351.93, NormSeries::U238Ra }, { 609.32, NormSeries::U238Ra }, { 768.36, NormSeries::U238Ra },
    { 934.06, NormSeries::U238Ra }, {1001.03, NormSeries::U238Ra }, {1120.29, NormSeries::U238Ra },
    {1155.19, NormSeries::U238Ra }, {1238.11, NormSeries::U238Ra }, {1280.96, NormSeries::U238Ra },
    {1377.67, NormSeries::U238Ra }, {1401.50, NormSeries::U238Ra }, {1407.98, NormSeries::U238Ra },
    {1509.23, NormSeries::U238Ra }, {1661.28, NormSeries::U238Ra }, {1729.60, NormSeries::U238Ra },
    {1764.49, NormSeries::U238Ra }, {1847.42, NormSeries::U238Ra }, {2118.55, NormSeries::U238Ra },
    {2204.21, NormSeries::U238Ra }, {2447.86, NormSeries::U238Ra },
    // Th-232 chain
    { 238.63, NormSeries::Th232 }, { 300.09, NormSeries::Th232 }, { 338.32, NormSeries::Th232 },
    { 463.00, NormSeries::Th232 }, { 510.77, NormSeries::Th232 }, { 583.19, NormSeries::Th232 },
    { 727.33, NormSeries::Th232 }, { 785.37, NormSeries::Th232 }, { 794.95, NormSeries::Th232 },
    { 860.56, NormSeries::Th232 }, { 911.20, NormSeries::Th232 }, { 964.77, NormSeries::Th232 },
    { 968.97, NormSeries::Th232 }, {1588.20, NormSeries::Th232 }, {1620.50, NormSeries::Th232 },
    {1630.63, NormSeries::Th232 }, {2614.51, NormSeries::Th232 },
    // K-40
    {1460.82, NormSeries::K40 },
    // Positron annihilation
    { 511.00, NormSeries::Annih }
  };

  /** Returns the NORM decay series `energy_kev` falls in (NormSeries::None if not a known NORM line),
   choosing the closest table entry within `tolerance_kev`. */
  NormSeries norm_series_for_energy( const double energy_kev, const double tolerance_kev )
  {
    NormSeries best = NormSeries::None;
    double bestDE = tolerance_kev;
    for( const NormLine &l : sk_norm_lines )
    {
      const double de = fabs( energy_kev - l.energy_kev );
      if( de <= bestDE ) { bestDE = de; best = l.series; }
    }
    return best;
  }//norm_series_for_energy(...)

  /** Returns true if `energy_kev` is within `tolerance_kev` of a known NORM/background line. */
  bool is_known_norm_energy( const double energy_kev, const double tolerance_kev )
  {
    return norm_series_for_energy( energy_kev, tolerance_kev ) != NormSeries::None;
  }//is_known_norm_energy(...)

  /** "Seed" gamma-ray lines (keV): Cs-137 (661.66) and Na-22 (511, 1274.53).

   On low-resolution (NaI) field detectors these nuclides are sometimes embedded as a stabilization/
   energy-calibration "seed", so the background spectrum also contains these peaks.  Unlike NORM,
   a seed's intensity is fixed (it does not vary by hours or ~10 ft of distance), so even a small but
   statistically-significant foreground excess indicates a genuine added foreground source.  Such peaks
   therefore get only a significance test - no blanket percent-elevation floor.
   */
  const double sk_seed_energies_kev[] = { 511.00, 661.66, 1274.53 };

  /** Returns true if `energy_kev` is within `tolerance_kev` of a known seed line (Cs-137 / Na-22). */
  bool is_seed_energy( const double energy_kev, const double tolerance_kev )
  {
    for( const double e : sk_seed_energies_kev )
    {
      if( fabs( energy_kev - e ) <= tolerance_kev )
        return true;
    }
    return false;
  }//is_seed_energy(...)

  /** Estimate the net (continuum-subtracted) counts a peak-sized feature would have in `hist` at the
   given `mean`/`fwhm`, using a simple linear side-band continuum.  Used to gauge whether a foreground
   peak at a known NORM energy is genuinely elevated over the background *data* in the case where the
   background auto-search did not fit a peak there (a weak background NORM line may simply have been
   too small to fit).  Returns 0 if the histogram/fwhm are unusable. */
  double estimate_net_counts_in_hist( const std::shared_ptr<const SpecUtils::Measurement> &hist,
                                      const double mean, const double fwhm )
  {
    if( !hist || (fwhm <= 0.0) )
      return 0.0;

    const double half = 1.25 * fwhm;            // peak window half-width
    const double lo = mean - half, hi = mean + half;
    const double side = fwhm;                   // side-band width on each side
    const double gross = hist->gamma_integral( static_cast<float>(lo), static_cast<float>(hi) );
    const double left  = hist->gamma_integral( static_cast<float>(lo - side), static_cast<float>(lo) );
    const double right = hist->gamma_integral( static_cast<float>(hi), static_cast<float>(hi + side) );
    const double cont_per_kev = (left + right) / (2.0 * side);   // continuum counts per keV
    const double continuum = cont_per_kev * (hi - lo);
    return gross - continuum;                   // net counts under the peak window
  }//estimate_net_counts_in_hist(...)


  /** Compares only RGB, as `WColor::operator==` also compares the color name. */
  bool same_rgb( const Wt::WColor &lhs, const Wt::WColor &rhs )
  {
    return (lhs.red() == rhs.red()) && (lhs.green() == rhs.green()) && (lhs.blue() == rhs.blue());
  }


  /** Returns the color for a peak being added with a source assigned, chosen the way showing
   reference lines for the source would, so each source gets its own color (see
   `ReferencePhotopeakDisplay::colorForNewPeak(...)`).  Without a reference-lines widget (e.g., tool
   tabs hidden on phones), picks from the color theme's reference-line palette, avoiding the colors
   of existing peaks and of `pending`.

   @param pending Source to color picks for peaks not yet added to the PeakModel (e.g., earlier
          peaks of the same batch); updated by this call.
   */
  Wt::WColor color_for_new_peak( const PeakDef &peak, InterSpec *viewer,
                                 std::map<std::string,Wt::WColor> &pending )
  {
    if( !viewer )
      return Wt::WColor();

    ReferencePhotopeakDisplay * const ref_lines = viewer->referenceLinesWidget();
    if( ref_lines )
      return ref_lines->colorForNewPeak( peak, pending );

    const string src = peak.sourceName();
    const shared_ptr<const ColorTheme> theme = viewer->getColorTheme();
    if( src.empty() || (theme && !theme->peaksTakeOnReferenceLineColor) )
      return Wt::WColor();

    const auto pending_pos = pending.find( src );
    if( pending_pos != end(pending) )
      return pending_pos->second;

    // Use the color of existing peaks of the source, if any, otherwise avoid all existing colors
    Wt::WColor color;
    vector<Wt::WColor> used;
    const shared_ptr<PeakModel> pmodel = viewer->peakModel();
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = pmodel ? pmodel->peaks() : nullptr;
    for( size_t i = 0; peaks && (i < peaks->size()); ++i )
    {
      const Wt::WColor &c = (*peaks)[i]->lineColor();
      if( c.isDefault() )
        continue;

      used.push_back( c );
      if( color.isDefault() && ((*peaks)[i]->sourceName() == src) )
        color = c;
    }//for( loop over existing peaks )

    if( color.isDefault() && theme )
    {
      const auto specific_pos = theme->referenceLineColorForSources.find( src );
      if( specific_pos != end(theme->referenceLineColorForSources) )
        color = specific_pos->second;
    }

    if( color.isDefault() )
    {
      for( const auto &p : pending )
        used.push_back( p.second );

      const vector<Wt::WColor> &palette = (theme && !theme->referenceLineColor.empty())
                                            ? theme->referenceLineColor
                                            : ReferencePhotopeakDisplay::sm_def_line_colors;
      color = ReferencePhotopeakDisplay::firstUnusedColor( palette, used, pending.size() );
    }//if( color.isDefault() )

    if( !color.isDefault() )
      pending[src] = color;

    return color;
  }//color_for_new_peak(...)


  /** Assigns `source` to `peak`, the way a user entering it into the peak table would, and colors it
   to match other peaks of the same source, or else a new color for the source (see
   `color_for_new_peak(...)`).

   A source containing "unknown" (or starting with "unk") sets the user label to "UNKNOWN" instead,
   which `ToolRegistry::executeGetUnidentifiedDetectedPeaks(...)` recognizes.  An empty source
   does nothing.

   @param other_peaks The other peaks in the spectrum, to match the color of.
   @param viewer The session, to pick new colors from; if null, no new color is picked.
   @param pending Source to color picks for peaks not yet added to the PeakModel; updated by this call.
   @returns false if `source` is not a valid source string.
   */
  bool assign_source_and_color( PeakDef &peak, const std::string &source,
                                const std::vector<std::shared_ptr<const PeakDef>> &other_peaks,
                                InterSpec *viewer, std::map<std::string,Wt::WColor> &pending )
  {
    if( source.empty() )
      return true;

    if( SpecUtils::icontains(source, "unknown") || SpecUtils::istarts_with(source, "unk") )
    {
      peak.setUserLabel( "UNKNOWN" );
      return true;
    }

    const PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( peak, source, 4.0 );
    if( res == PeakModel::SetGammaSource::FailedSourceChange )
      return false;

    vector<Wt::WColor> src_colors;
    size_t num_same_src_peaks = 0;
    for( const shared_ptr<const PeakDef> &srcPeak : other_peaks )
    {
      if( (srcPeak->parentNuclide() || srcPeak->xrayElement() || srcPeak->reaction())
          && srcPeak->parentNuclide()==peak.parentNuclide()
          && srcPeak->xrayElement()==peak.xrayElement()
          && srcPeak->reaction()==peak.reaction() )
      {
        num_same_src_peaks += 1;

        const Wt::WColor &c = srcPeak->lineColor();
        const bool have_color = std::any_of( begin(src_colors), end(src_colors), [&c]( const Wt::WColor &rhs ){
          return same_rgb( c, rhs );
        } );

        if( !c.isDefault() && !have_color )
          src_colors.push_back( c );
      }
    }//for( const auto &srcPeak : other_peaks )

    if( src_colors.size() == 1 )
      peak.setLineColor( src_colors[0] );
    else if( num_same_src_peaks == 0 )
      peak.setLineColor( color_for_new_peak( peak, viewer, pending ) );

    return true;
  }//assign_source_and_color(...)


  /** Assigns the source to the newly fit peak in `peaks_to_add_in` (the one without a source,
   nearest `energy`), for the legacy `fit_user_peak` path.

   @returns The newly fit peak (possibly a copy, if a source was assigned).
   */
  shared_ptr<const PeakDef> assign_source_to_user_fit_peak( vector<shared_ptr<const PeakDef>> &peaks_to_add_in,
                                                const double energy, const string &source,
                                                const vector<shared_ptr<const PeakDef>> &orig_peaks,
                                                InterSpec *viewer )
  {
    // The new peak may share a ROI with previously existing peaks, so we have to figure out
    //  what peak we just fit.
    //  TODO: For the moment, we'll just match up to nearest peak, but in a normal user double click for a peak, we match up peaks in ROI to previous peaks, and then take whats left-over, which is maybe a better way.
    size_t fit_peak_index = peaks_to_add_in.size();
    for( size_t index = 0; index < peaks_to_add_in.size(); ++index )
    {
      const shared_ptr<const PeakDef> &peak = peaks_to_add_in[index];
      if( peak->parentNuclide() || peak->xrayElement() || peak->reaction() )
        continue;  //We know new peak doesnt have a source assigned yet

      const double diff = fabs(peak->mean() - energy);
      if( (fit_peak_index >= peaks_to_add_in.size())
         || (diff < fabs(peaks_to_add_in[fit_peak_index]->mean() - energy)) )
      {
        fit_peak_index = index;
      }
    }//for( loop over peaks_to_add_in )

    assert( fit_peak_index < peaks_to_add_in.size() );
    if( fit_peak_index >= peaks_to_add_in.size() )
      throw runtime_error( "Could not identify newly fit peak." );

    if( source.empty() )
      return peaks_to_add_in[fit_peak_index];

    shared_ptr<PeakDef> new_fit_peak = make_shared<PeakDef>( *peaks_to_add_in[fit_peak_index] );
    map<string,Wt::WColor> pending_colors;
    if( !assign_source_and_color( *new_fit_peak, source, orig_peaks, viewer, pending_colors ) )
      throw runtime_error( "Peak was fit, but source string '" + source + "' was not valid to set source." );

    // Long-standing behavior, kept as-is: an "UNKNOWN" label only goes on the returned peak, not the
    //  peak added to the model (probably unintended).
    const bool unknown_label = (SpecUtils::icontains(source, "unknown") || SpecUtils::istarts_with(source, "unk"));
    if( !unknown_label )
      peaks_to_add_in[fit_peak_index] = new_fit_peak;

    return new_fit_peak;
  }//assign_source_to_user_fit_peak(...)


  /** Fits a Gaussian of `diag.expectedFwhm`, fixed at `energy`, with a quadratic continuum, alongside
   the `nearby_peaks` whose tails reach its window, and records it in `diag`.  Thread-safe. */
  void fixed_position_fit( AnalystChecks::UnfitPeakDiagnostic &diag, const double energy,
                           const shared_ptr<const SpecUtils::Measurement> &data,
                           const vector<shared_ptr<const PeakDef>> &nearby_peaks )
  {
    diag.fixedFitValid = false;
    const double expected_fwhm = diag.expectedFwhm;
    if( !data || !(expected_fwhm > 0.0) )
      return;

    vector<pair<double,double>> neighbours;
    for( const shared_ptr<const PeakDef> &p : nearby_peaks )
    {
      if( !p || !p->gausPeak() )
        continue;
      const double de = fabs( p->mean() - energy );
      if( (de < 0.5*expected_fwhm) || (de > (2.0*expected_fwhm + 1.5*p->fwhm())) )
        continue;
      const bool dup = std::any_of( begin(neighbours), end(neighbours), [&p]( const pair<double,double> &n ){
        return fabs( n.first - p->mean() ) < 0.25*p->fwhm(); } );
      if( !dup )
        neighbours.push_back( { p->mean(), p->fwhm() } );
    }//for( loop over nearby peaks )

    const FitPeaksForNuclides::detail::FixedShapePeakFit fit
                = FitPeaksForNuclides::detail::fixed_shape_peak_fit( data, energy, expected_fwhm, neighbours );
    if( !fit.valid )
      return;

    const size_t num_pars = 3 + 1 + neighbours.size();
    diag.fixedFitValid = true;
    diag.fixedFitArea = fit.area;
    diag.fixedFitAreaUncert = fit.area_uncert;
    diag.fixedFitChi2Dof = (fit.num_channels > num_pars) ? (fit.chi2 / (fit.num_channels - num_pars)) : -1.0;
    diag.fixedFitLowerEnergy = fit.lower_energy;
    diag.fixedFitUpperEnergy = fit.upper_energy;
  }//fixed_position_fit(...)


  /** What the data show at `energy`, where the automated-threshold search accepted no peak: the
   candidate the lenient user-click thresholds find, and a Gaussian of the expected FWHM fit at
   exactly `energy`.  Touches no session state, so may be called from a worker thread.
   */
  AnalystChecks::UnfitPeakDiagnostic diagnose_unfit_peak( const double energy,
                                    const std::string &auto_reject_reason,
                                    const double expected_fwhm,
                                    const std::string &expected_fwhm_source,
                                    const shared_ptr<const SpecUtils::Measurement> &data,
                                    const vector<shared_ptr<const PeakDef>> &orig_peaks,
                                    const shared_ptr<const DetectorPeakResponse> &det,
                                    const shared_ptr<const deque<shared_ptr<const PeakDef>>> &auto_search_peaks,
                                    const shared_ptr<const PeakFitDetPrefs> &fit_prefs )
  {
    AnalystChecks::UnfitPeakDiagnostic diag;
    diag.expectedFwhm = expected_fwhm;
    diag.expectedFwhmSource = expected_fwhm_source;
    diag.reason = auto_reject_reason.empty() ? string("the automated peak search found no peak there")
                                             : ("the automated peak search rejected the candidate: " + auto_reject_reason);

    // The lenient thresholds of a user double-click (same candidate search and ROI as the automated
    //  search, which uses 5 pixels per keV too).
    try
    {
      const pair<vector<shared_ptr<const PeakDef>>, vector<shared_ptr<const PeakDef>>> lenient
              = searchForPeakFromUser( energy, 5.0, data, orig_peaks, det, auto_search_peaks, fit_prefs );
      for( const shared_ptr<const PeakDef> &p : lenient.first )
      {
        if( !p || !p->gausPeak() )
          continue;

        // The ROI's peaks are all new objects, including refits of existing peaks; those are not the
        //  candidate, so they are recognized by energy.
        const bool is_existing = std::any_of( begin(orig_peaks), end(orig_peaks),
                                              [&p]( const shared_ptr<const PeakDef> &orig ){
          return orig && (fabs( orig->mean() - p->mean() ) < 0.25*std::max( orig->fwhm(), p->fwhm() ));
        } );

        if( !is_existing && (!diag.candidate || (fabs( p->mean() - energy ) < fabs( diag.candidate->mean() - energy ))) )
          diag.candidate = p;
      }
      if( diag.candidate )
        diag.candidateRoiPeaks = lenient.first;
    }catch( std::exception & )
    {
    }

    vector<shared_ptr<const PeakDef>> nearby_peaks = orig_peaks;
    if( auto_search_peaks )
      nearby_peaks.insert( end(nearby_peaks), begin(*auto_search_peaks), end(*auto_search_peaks) );
    fixed_position_fit( diag, energy, data, nearby_peaks );

    return diag;
  }//diagnose_unfit_peak(...)


  /** A few sentences summarizing an `UnfitPeakDiagnostic`, e.g., for an error message. */
  std::string unfit_peak_summary( const double energy, const AnalystChecks::UnfitPeakDiagnostic &diag )
  {
    char buffer[512] = { '\0' };
    string answer = diag.reason;
    if( !answer.empty() )
      answer[0] = static_cast<char>( std::toupper( answer[0] ) );
    answer += ".";

    if( diag.candidate )
    {
      const shared_ptr<const PeakDef> &c = diag.candidate;
      const double nsigma = (c->amplitudeUncert() > 0.0) ? (c->amplitude() / c->amplitudeUncert()) : 0.0;
      snprintf( buffer, sizeof(buffer), "  The lower thresholds of a user double-click give a peak at %.2f keV, with"
                " FWHM %.3g keV (expected %.3g keV), and area %.4g +- %.2g (%.2g sigma).",
                c->mean(), c->fwhm(), diag.expectedFwhm, c->amplitude(), c->amplitudeUncert(), nsigma );
      answer += buffer;
    }else
    {
      answer += "  Not even the lower thresholds of a user double-click find a peak there.";
    }

    if( diag.fixedFitValid )
    {
      const double nsigma = (diag.fixedFitAreaUncert > 0.0) ? (diag.fixedFitArea / diag.fixedFitAreaUncert) : 0.0;
      snprintf( buffer, sizeof(buffer), "  A Gaussian with the expected FWHM (%.3g keV), fixed at %.2f keV, fits"
                " %.4g +- %.2g counts (%.2g sigma).", diag.expectedFwhm, energy, diag.fixedFitArea,
                diag.fixedFitAreaUncert, nsigma );
      answer += buffer;
    }

    return answer;
  }//unfit_peak_summary(...)


  /** The error message when a single-energy peak fit adds no peak. */
  std::string unfit_peak_error_message( const double energy, const bool within_existing_roi,
                                        const AnalystChecks::UnfitPeakDiagnostic &diag )
  {
    string msg = "Could not fit peak at " + SpecUtils::printCompact(energy, 4) + " keV";
    if( within_existing_roi )
      msg += " - the energy is within the ROI of an existing peak, so maybe the peak you wanted already exists"
             " in the analysis peak list, or another peak just cant be fit in the same ROI";
    msg += ".  " + unfit_peak_summary( energy, diag )
           + "  To get this measurement as a result instead, set DoNotAddToAnalysisPeaks to true; to add a peak"
             " here regardless of its significance, give the energy as an array, e.g., [" + SpecUtils::printCompact(energy, 5) + "].";
    return msg;
  }//unfit_peak_error_message(...)


  /** The detector class for expected-FWHM estimates: from the peak-fit preferences (`fit_prefs`, else
   the spectrum's, else the DRF's), else classified from the spectrum; an undetermined class is
   taken as HPGe when the spectrum looks high resolution.  GUI thread (`is_likely_high_res` reads the
   session). */
  PeakFitUtils::CoarseResolutionType resolve_det_type( InterSpec *interspec,
                                                       const shared_ptr<const SpecMeas> &meas,
                                                       const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                                       shared_ptr<const PeakFitDetPrefs> fit_prefs = nullptr )
  {
    if( !fit_prefs && meas )
      fit_prefs = meas->peakFitDetPrefs();
    if( !fit_prefs && meas && meas->detector() )
      fit_prefs = meas->detector()->peakFitDetPrefs();

    PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::effective_det_type( fit_prefs, spectrum, meas );
    if( (det_type == PeakFitUtils::CoarseResolutionType::Unknown) && interspec && PeakFitUtils::is_likely_high_res( interspec ) )
      det_type = PeakFitUtils::CoarseResolutionType::High;
    return det_type;
  }//resolve_det_type(...)


  /** The expected-FWHM function for a single-energy peak fit, from the DRF, the analysis peaks, and
   the automated-search peaks (either may be empty).  Touches no session state. */
  std::function<float(double)> fwhm_fcn_for_peak_fit( const PeakFitUtils::CoarseResolutionType det_type,
                                    const shared_ptr<const SpecUtils::Measurement> &data,
                                    const shared_ptr<const DetectorPeakResponse> &det,
                                    const vector<shared_ptr<const PeakDef>> &orig_peaks,
                                    const shared_ptr<const deque<shared_ptr<const PeakDef>>> &auto_search_peaks,
                                    std::string &source )
  {
    auto resolution_peaks = make_shared<deque<shared_ptr<const PeakDef>>>( begin(orig_peaks), end(orig_peaks) );
    if( auto_search_peaks )
      resolution_peaks->insert( end(*resolution_peaks), begin(*auto_search_peaks), end(*auto_search_peaks) );

    return AnalystChecks::make_expected_fwhm_fcn( det, resolution_peaks, data, det_type, &source );
  }//fwhm_fcn_for_peak_fit(...)


  /** If `diag`'s expected FWHM is only a generic estimate (no peaks were at hand), gets it from the
   detected peaks instead - running the automated peak search if it has not been, so GUI thread
   only; used by the synchronous `fit_user_peak` - and redoes the fixed-position fit with it. */
  void refine_unfit_diagnostic_fwhm( AnalystChecks::UnfitPeakDiagnostic &diag, const double energy,
                                     InterSpec *interspec, const SpecUtils::SpectrumType type,
                                     const shared_ptr<const SpecUtils::Measurement> &data,
                                     const shared_ptr<const DetectorPeakResponse> &det,
                                     const shared_ptr<const PeakFitDetPrefs> &fit_prefs )
  {
    if( !interspec || !SpecUtils::istarts_with( diag.expectedFwhmSource, "generic" ) )
      return;

    try
    {
      AnalystChecks::DetectedPeaksOptions opts;
      opts.specType = type;
      const AnalystChecks::DetectedPeakStatus detected = AnalystChecks::detected_peaks( opts, interspec );
      if( detected.peaks.empty() )
        return;

      auto peaks = make_shared<deque<shared_ptr<const PeakDef>>>( begin(detected.peaks), end(detected.peaks) );
      const PeakFitUtils::CoarseResolutionType det_type = resolve_det_type( interspec, nullptr, data, fit_prefs );

      string source;
      const std::function<float(double)> fwhm_fcn = AnalystChecks::make_expected_fwhm_fcn( det, peaks, data, det_type, &source );
      diag.expectedFwhm = fwhm_fcn( energy );
      diag.expectedFwhmSource = source;
      fixed_position_fit( diag, energy, data, detected.peaks );
    }catch( std::exception & )
    {
    }
  }//refine_unfit_diagnostic_fwhm(...)


  /** The decay cascades of a nuclide (and its progeny, grown in over its default age), with x-ray
   and annihilation members, for finding coincident photon pairs; `per_parent_decay` is set to the
   factor that turns a cascade's joint emission probability into one per decay of the nuclide.
   */
  std::vector<ceelo::DecayCascade> nuclide_cascades( const SandiaDecay::Nuclide *nuc, double &per_parent_decay )
  {
    per_parent_decay = 0.0;
    if( !nuc )
      return {};

    const double age = PeakDef::defaultDecayTime( nuc, nullptr ) / PhysicalUnits::second;

    ceelo::cascade_adapter::CascadeOptions copt;
    copt.age_seconds = age;
    copt.prompt_equilibrium = false;
    copt.include_xrays = true;
    copt.vacancy_xray_model = false;  // x-rays linked to the gammas, and K x-rays to each other
    copt.include_annihilation = true;

    // The adapter starts from a unit parent activity, decayed to `age`.
    const double parent_activity = std::exp( -nuc->decayConstant() * age );
    if( !(parent_activity > 0.0) || IsInf(parent_activity) )
      return {};
    per_parent_decay = 1.0 / parent_activity;

    return ceelo::cascade_adapter::build_cascades( nuc, copt );
  }//nuclide_cascades(...)


  /** Probability per decay that two coincident photons, within `tol1` of `e1` and `tol2` of `e2`,
   are emitted (summed over every such pair of lines, including two photons of one line, such as
   two Kalpha1 x-rays from two K vacancies); `type` is set to the kind of pair with the largest
   probability, e.g. "x-ray-gamma" or "x-ray-x-ray".  Zero if there is no such coincidence.
   */
  double cascade_pair_probability( const std::vector<ceelo::DecayCascade> &cascades,
                                   const double per_parent_decay,
                                   const double e1, const double tol1,
                                   const double e2, const double tol2,
                                   std::string &type )
  {
    const auto type_name = []( const ceelo::CascadeParticleType t ) -> const char * {
      switch( t )
      {
        case ceelo::CascadeParticleType::Gamma:    return "gamma";
        case ceelo::CascadeParticleType::Xray:     return "x-ray";
        case ceelo::CascadeParticleType::Annih511: return "annihilation";
      }
      return "photon";
    };

    double total = 0.0, largest = 0.0;
    type.clear();
    for( const ceelo::DecayCascade &dc : cascades )
    {
      std::set<std::pair<size_t,size_t>> seen;  // a link may be listed from either member
      for( size_t a = 0; a < dc.members.size(); ++a )
      {
        const ceelo::CascadeMember &ma = dc.members[a];

        // Two photons of this one line (a link always joins two different members)
        if( (ma.self_pair_prob > 0.0)
           && (fabs(ma.energy_keV - e1) <= tol1) && (fabs(ma.energy_keV - e2) <= tol2) )
        {
          const double prob = per_parent_decay * dc.branch_weight * ma.self_pair_prob;
          total += prob;
          if( prob > largest )
          {
            largest = prob;
            type = string( type_name( ma.type ) ) + "-" + type_name( ma.type );
          }
        }//if( member a pairs with itself )

        for( const ceelo::CoincidenceLink &link : ma.coincident )
        {
          const size_t b = link.partner;
          if( (b >= dc.members.size()) || (b == a) || !seen.insert( {std::min(a,b), std::max(a,b)} ).second )
            continue;

          const ceelo::CascadeMember &mb = dc.members[b];
          const bool forward = (fabs(ma.energy_keV - e1) <= tol1) && (fabs(mb.energy_keV - e2) <= tol2);
          const bool backward = (fabs(ma.energy_keV - e2) <= tol2) && (fabs(mb.energy_keV - e1) <= tol1);
          if( !forward && !backward )
            continue;

          const double prob = per_parent_decay * dc.branch_weight * ma.intensity * link.prob;
          if( !(prob > 0.0) )
            continue;

          total += prob;
          if( prob > largest )
          {
            largest = prob;
            // Gamma last, so the names read "x-ray-gamma", "annihilation-gamma"
            const bool a_first = (ma.type != ceelo::CascadeParticleType::Gamma);
            type = string( type_name( a_first ? ma.type : mb.type ) ) + "-" + type_name( a_first ? mb.type : ma.type );
          }
        }//for( loop over links of member a )
      }//for( loop over members )
    }//for( loop over cascades )

    return total;
  }//cascade_pair_probability(...)


}//anonymous namespace


namespace AnalystChecks
{
  double elevation_test_uncert( const PeakDef &peak, const vector<shared_ptr<const PeakDef>> &roi_peaks,
                                const shared_ptr<const SpecUtils::Measurement> &data )
  {
    if( data && (peak.amplitude() > 0.0) )
    {
      const double z = PeakFitLM::peak_detection_significance( peak, roi_peaks, data, /*chi2_weights=*/ true );
      if( z > 0.0 )
        return peak.amplitude() / z;
    }
    return peak.amplitudeUncert();
  }//elevation_test_uncert(...)


  DetectedPeakStatus detected_peaks( const DetectedPeaksOptions& options, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");

    // Validate nonBackgroundPeaksOnly option
    if( options.nonBackgroundPeaksOnly )
    {
      if( options.specType != SpecUtils::SpectrumType::Foreground )
        throw std::runtime_error("nonBackgroundPeaksOnly may only be specified for Foreground spectrum");

      std::shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
      if( !background )
        throw std::runtime_error("nonBackgroundPeaksOnly may only be specified when a Background spectrum is loaded");
    }

    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for "
          + std::string(SpecUtils::descriptionText(options.specType))
          + " spectrum");
    }

    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for " 
          + std::string(SpecUtils::descriptionText(options.specType)) 
          + " spectrum");
    }

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(options.specType);
    if (!spectrum) {
      throw std::runtime_error("No spectrum displayed for " 
          + std::string(SpecUtils::descriptionText(options.specType)) 
          + " spectrum");
    }
    
    shared_ptr<const deque<shared_ptr<const PeakDef>>> user_peaks = meas->peaks(sample_nums);

    // Whether the detector is high-resolution (HPGe); used for the peak search and recovery below.
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);

    // Peak-fit prefs (detector resolution type, skew, ...) to search/fit with.
    const shared_ptr<const PeakFitDetPrefs> fitPrefs = meas->peakFitDetPrefs();

    // Obtain the automated-search peaks through the shared, de-duplicated accessor - this launches
    // at most one search per (SpecMeas, sample numbers) and coalesces with any concurrent search
    // (GUI hint-peak search, other tools), instead of each caller running its own.  We are on the
    // GUI thread here, so the accessor may snapshot meas->peaks() itself for the search seed.
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks
        = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( meas, sample_nums, spectrum,
                                                                    meas->detector(), fitPrefs ).get();

    vector<shared_ptr<const PeakDef>> all_peaks;
    vector<shared_ptr<const PeakDef>> user_analysis_peaks;

    if (user_peaks)
    {
      all_peaks.insert(all_peaks.end(), user_peaks->begin(), user_peaks->end());
      user_analysis_peaks.insert(user_analysis_peaks.end(), user_peaks->begin(), user_peaks->end());
    }

    // We will add auto-search peaks only if there isnt already a user peak at essentially the same energy.
    // Track which auto-search continuum pointers had peaks suppressed, and which user analysis
    // peak(s) caused each suppression, so we can unify ROI bounds for remaining siblings.
    if (auto_peaks) {
      //Lambda to find the nearest peak, so far, to a given energy.
      const auto nearest_peak = [&all_peaks](const float energy) -> std::shared_ptr<const PeakDef> {
        std::shared_ptr<const PeakDef> nearest;
        double minDE = std::numeric_limits<double>::infinity();

        for (const std::shared_ptr<const PeakDef> &peak : all_peaks) {
          const double dE = fabs(peak->mean() - energy);
          const bool in_roi = (energy > peak->lowerX()) && (energy < peak->upperX());
          const double pk_sigma = peak->gausPeak() ? peak->sigma() : 0.25 * peak->roiWidth();
          const bool within_sigma = (dE < pk_sigma);
          if ((dE < minDE) && (in_roi || within_sigma)) {
            minDE = dE;
            nearest = peak;
          }
        }

        return nearest;
      };

      // Map from suppressed auto-search continuum → user analysis peak(s) that caused suppression.
      // A multimap handles the case where multiple auto-search peaks from the same continuum are
      // each suppressed by different user analysis peaks.
      std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>> suppressed_to_nearpeaks;

      for (const shared_ptr<const PeakDef> &peak : *auto_peaks) {
        std::shared_ptr<const PeakDef> nearpeak = nearest_peak(peak->mean());
        const double peak_sigma = peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth();
        if (!nearpeak || (fabs(nearpeak->mean() - peak->mean()) > peak_sigma))
        {
          all_peaks.push_back(peak);
        }else
        {
          // Auto-search peak suppressed by nearpeak; record for sibling ROI unification
          suppressed_to_nearpeaks[peak->continuum().get()].push_back( nearpeak );
        }
      }

      // Post-processing: unify ROI bounds for auto-search siblings whose peer was suppressed.
      // When an auto-search peak is suppressed, remaining siblings from the same auto-search ROI
      // still reference the original wide ROI bounds, which may overlap with the user analysis
      // peak's ROI. Fix by giving all involved peaks a shared continuum with unified bounds.
      if( !suppressed_to_nearpeaks.empty() )
      {
        // Track which user analysis peaks have already been replaced (to avoid double-replacement
        // when multiple siblings reference the same suppression)
        std::map<const PeakDef *, std::shared_ptr<PeakContinuum>> replaced_user_peaks;

        for( size_t i = 0; i < all_peaks.size(); ++i )
        {
          const std::shared_ptr<const PeakDef> &peak = all_peaks[i];
          const std::map<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>>::const_iterator
            it = suppressed_to_nearpeaks.find( peak->continuum().get() );
          if( it == suppressed_to_nearpeaks.end() )
            continue;

          // This peak is a kept auto-search sibling of a suppressed peak.
          // Check if we've already created a unified continuum for this suppression group.
          const std::map<const PeakDef *, std::shared_ptr<PeakContinuum>>::const_iterator
            replaced_it = replaced_user_peaks.find( it->second.front().get() );
          std::shared_ptr<PeakContinuum> unified_cont;

          if( replaced_it != replaced_user_peaks.end() )
          {
            // Reuse previously created unified continuum and expand bounds if needed
            unified_cont = replaced_it->second;
            const double new_lower = std::min( unified_cont->lowerEnergy(), peak->lowerX() );
            const double new_upper = std::max( unified_cont->upperEnergy(), peak->upperX() );
            unified_cont->setRange( new_lower, new_upper );
          }else
          {
            // Create a new unified continuum; start from the first user analysis peak's continuum
            const std::vector<std::shared_ptr<const PeakDef>> &nearpeaks = it->second;
            unified_cont = std::make_shared<PeakContinuum>( *nearpeaks.front()->continuum() );

            // Compute union of all involved ROI bounds
            double union_lower = std::min( peak->lowerX(), unified_cont->lowerEnergy() );
            double union_upper = std::max( peak->upperX(), unified_cont->upperEnergy() );
            for( const std::shared_ptr<const PeakDef> &np : nearpeaks )
            {
              union_lower = std::min( union_lower, np->lowerX() );
              union_upper = std::max( union_upper, np->upperX() );
            }
            unified_cont->setRange( union_lower, union_upper );

            // Replace the user analysis peaks in all_peaks (and user_analysis_peaks) with copies
            // using the unified continuum, so the copies inherit the analysis-peak status.
            for( const std::shared_ptr<const PeakDef> &np : nearpeaks )
            {
              replaced_user_peaks[np.get()] = unified_cont;
              std::shared_ptr<PeakDef> np_copy = std::make_shared<PeakDef>( *np );
              np_copy->setContinuum( unified_cont );

              for( size_t j = 0; j < all_peaks.size(); ++j )
              {
                if( all_peaks[j].get() == np.get() )
                {
                  all_peaks[j] = np_copy;
                  break;
                }
              }

              // Also replace in user_analysis_peaks so isAnalysisPeak detection
              // (which uses pointer equality) recognizes the copy.
              for( size_t j = 0; j < user_analysis_peaks.size(); ++j )
              {
                if( user_analysis_peaks[j].get() == np.get() )
                {
                  user_analysis_peaks[j] = np_copy;
                  break;
                }
              }
            }
          }

          // Replace this auto-search sibling with a copy using the unified continuum
          std::shared_ptr<PeakDef> sibling_copy = std::make_shared<PeakDef>( *peak );
          sibling_copy->setContinuum( unified_cont );
          all_peaks[i] = sibling_copy;
        }//for( all_peaks )
      }//if( !suppressed_to_nearpeaks.empty() )
    }

    std::sort(begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr);

    // Apply nonBackgroundPeaksOnly filtering if requested
    if( options.nonBackgroundPeaksOnly )
    {
      // Get background spectrum and its peaks
      std::shared_ptr<SpecMeas> background_meas = interspec->measurment(SpecUtils::SpectrumType::Background);
      std::shared_ptr<const SpecUtils::Measurement> background_spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
      
      assert( background_spectrum );
      if( !background_spectrum )
        throw std::logic_error( "Somehow lost background spectrum being loaded?" );
      const set<int> background_sample_nums = interspec->displayedSamples(SpecUtils::SpectrumType::Background);

      // Ensure the background's automated-search peaks include any real NORM peaks recovered from
      // underneath the foreground peaks, so non-elevated background lines that the background
      // auto-search may have missed are not mis-flagged as elevated.  This obtains both foreground
      // and background auto-search peaks via the shared accessor, runs recovery, and persists the
      // augmented background set (guarded to run once per foreground/background pairing).  We are on
      // the GUI thread, where blocking on the search/recovery is acceptable (matches prior behavior).
      PeakSearchGuiUtils::ensure_background_peaks_recovered( meas, sample_nums, spectrum,
                                                            background_meas, background_sample_nums,
                                                            background_spectrum, fitPrefs );

      // Fetch the (now recovered/augmented) background auto peaks via the shared, de-duplicated
      // accessor - a cache hit after the call above.
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> background_auto_peaks
          = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( background_meas,
                          background_sample_nums, background_spectrum,
                          background_meas->detector(), fitPrefs ).get();
      const vector<shared_ptr<const PeakDef>> background_auto_peak_vec = background_auto_peaks
          ? vector<shared_ptr<const PeakDef>>( begin(*background_auto_peaks), end(*background_auto_peaks) )
          : vector<shared_ptr<const PeakDef>>{};

      // Filter foreground peaks based on background
      vector<shared_ptr<const PeakDef>> filtered_peaks;
      vector<shared_ptr<const PeakDef>> filtered_analysis_peaks;
      const double foreground_live_time = (spectrum->live_time() > 0.0f) ? spectrum->live_time() : 1.0f;
      const double background_live_time = (background_spectrum->live_time() > 0.0f ? background_spectrum->live_time() : 1.0f);

#if( PERFORM_DEVELOPER_CHECKS )
      // Debug dump of the foreground candidate peaks and the (recovered) background auto-search peaks,
      // to aid tuning of the elevated-above-background logic.  Per-peak KEEP/drop decisions follow below.
      {
        cerr << "=== detected_peaks: elevated-above-background debug (fg live=" << foreground_live_time
             << "s, bg live=" << background_live_time << "s) ===" << endl;
        cerr << " Background auto-search peaks ("
             << (background_auto_peaks ? background_auto_peaks->size() : size_t(0)) << "):" << endl;
        if( background_auto_peaks )
        {
          for( const shared_ptr<const PeakDef> &bp : *background_auto_peaks )
            cerr << "  bg E=" << bp->mean() << " keV  amp=" << bp->amplitude()
                 << " cps=" << (bp->amplitude()/background_live_time) << endl;
        }
        cerr << " Foreground candidate peaks (" << all_peaks.size() << "):" << endl;
      }
#endif

      // Two-pass elevated-above-background filter.  Pass 1 makes the clear accept/reject decisions and
      // *defers* borderline NORM peaks (a known NORM line that fails the 4-sigma bar but is still above
      // 2 sigma).  Pass 2 rescues a borderline NORM peak only if another line of the SAME decay series
      // is already accepted and this peak's foreground/background ratio to the nearest such line is
      // consistent with them being the same chain (secular equilibrium) - otherwise it stays dropped.
      struct FgPeakDecision {
        shared_ptr<const PeakDef> peak;
        bool include = false;
        bool is_user = false;
        bool is_seed = false;
        bool borderline_norm = false;   // deferred to pass 2
        NormSeries series = NormSeries::None;
        double mean = 0.0;
        double fg_cps = 0.0;
        double bg_cps = 0.0;            // >0 only when a matching background peak was found
      };

      const double norm_min_sigma = 4.0;         // net-significance bar for known NORM lines
      const double norm_borderline_sigma = 2.0;  // below this a NORM peak is simply background
      const double seed_min_sigma = 2.25;
      const double other_min_sigma = 2.25;       // a non-NORM, non-seed peak elevated >20%

      vector<FgPeakDecision> decisions;
      decisions.reserve( all_peaks.size() );

      // --- Pass 1: per-peak decision (deferring borderline NORM peaks) ---
      for( const shared_ptr<const PeakDef> &fg_peak : all_peaks )
      {
        FgPeakDecision dec;
        dec.peak = fg_peak;
        dec.mean = fg_peak->mean();
        dec.fg_cps = fg_peak->amplitude() / foreground_live_time;
        dec.is_user = (std::find(user_analysis_peaks.begin(), user_analysis_peaks.end(), fg_peak) != user_analysis_peaks.end());

        // Known NORM/background lines are held to a stricter bar (they are expected in any background
        // and fluctuate).  Seed nuclides (Cs-137/Na-22) on low-res detectors get only a significance
        // test (stable background intensity).  Tolerance scales with resolution (peak FWHM).
        const double norm_tol = std::max( 2.0, 0.5 * fg_peak->fwhm() );
        dec.series = norm_series_for_energy( fg_peak->mean(), norm_tol );
        const bool is_norm = (dec.series != NormSeries::None);
        dec.is_seed = (!isHPGe && is_seed_energy( fg_peak->mean(), norm_tol ));

        bool include_peak = true;

        if( background_auto_peaks )
        {
          // Find the closest matching background peak
          shared_ptr<const PeakDef> closest_bg_peak;
          double smallest_energy_diff = std::numeric_limits<double>::infinity();
          for( const shared_ptr<const PeakDef> &bg_peak : *background_auto_peaks )
          {
            assert( fg_peak != bg_peak );
            const double energy_diff = fabs(fg_peak->mean() - bg_peak->mean());
            const double avg_fwhm = 0.75 * (fg_peak->fwhm() + bg_peak->fwhm()) / 2.0;
            if( (energy_diff < avg_fwhm) && (energy_diff < smallest_energy_diff) )
            {
              closest_bg_peak = bg_peak;
              smallest_energy_diff = energy_diff;
            }
          }

          if( closest_bg_peak )
          {
            const double fg_cps = dec.fg_cps;
            const double bg_cps = (closest_bg_peak->amplitude() / background_live_time);
            dec.bg_cps = bg_cps;
            const double fg_amp_uncert = elevation_test_uncert( *fg_peak, all_peaks, spectrum );
            const double bg_amp_uncert = elevation_test_uncert( *closest_bg_peak, background_auto_peak_vec,
                                                                background_spectrum );
            const bool have_uncert = (fg_amp_uncert > 0.0) && (bg_amp_uncert > 0.0);

            double sigma_elevation = 0.0;
            if( have_uncert )
            {
              const double fg_cps_uncert = fg_amp_uncert / foreground_live_time;
              const double bg_cps_uncert = bg_amp_uncert / background_live_time;
              const double combined_uncert = sqrt(fg_cps_uncert*fg_cps_uncert + bg_cps_uncert*bg_cps_uncert);
              if( combined_uncert > 0.0 )
                sigma_elevation = ((fg_cps - bg_cps) / combined_uncert);
            }

            if( dec.is_seed )
            {
              // Seed nuclide on a low-res detector: significance test only, no percent floor.
              include_peak = have_uncert ? (sigma_elevation > seed_min_sigma) : (fg_cps > (bg_cps * 1.20));
            }else if( is_norm )
            {
              // NORM: strict 4-sigma bar; 2-4 sigma is deferred to the pass-2 chain-consistency check;
              // below 2 sigma is background.  (Flat 20% ratio dropped - meaningless at low counts.)
              if( have_uncert )
              {
                if( sigma_elevation > norm_min_sigma )
                  include_peak = true;
                else if( sigma_elevation > norm_borderline_sigma )
                { include_peak = false; dec.borderline_norm = true; }
                else
                  include_peak = false;
              }else
              {
                include_peak = (fg_cps > (bg_cps * 1.75));
              }
            }else
            {
              // Non-NORM, non-seed: original behavior - >20% elevated, then >2.25 sigma if available.
              const double elevation_threshold = 1.20;
              if( fg_cps > (bg_cps * elevation_threshold) )
                include_peak = have_uncert ? (sigma_elevation > other_min_sigma) : true;
              else
                include_peak = false;
            }//if( is_seed ) / else if( is_norm ) / else
          }else
          {
            // No matching background peak.  Non-NORM (and seed) energies keep the existing behavior:
            // the background recovery (run above) already fit-tested here and found nothing, so treat
            // the foreground peak as a real source (include_peak stays true).  For a known NORM energy,
            // do NOT auto-pass - require it to be clearly elevated over the background *data* counts in
            // the peak region (a weak background NORM line may simply have been too small to fit).
            if( is_norm && !dec.is_seed )
            {
              const double bg_net = estimate_net_counts_in_hist( background_spectrum, fg_peak->mean(), fg_peak->fwhm() );
              const double bg_net_cps = std::max( 0.0, bg_net ) / background_live_time;
              const double fg_net_cps = fg_peak->amplitude() / foreground_live_time;
              const double fg_cps_uncert = elevation_test_uncert( *fg_peak, all_peaks, spectrum ) / foreground_live_time;
              const double bg_net_uncert_cps = sqrt( std::max( fabs(bg_net), 1.0 ) ) / background_live_time;
              const double combined = sqrt( fg_cps_uncert*fg_cps_uncert + bg_net_uncert_cps*bg_net_uncert_cps );
              const double sigma_vs_data = (combined > 0.0) ? ((fg_net_cps - bg_net_cps) / combined) : 0.0;
              include_peak = (sigma_vs_data > norm_min_sigma);
            }
          }//if( closest_bg_peak )
        }//if( background_auto_peaks )

        dec.include = include_peak;
        decisions.push_back( std::move(dec) );
      }//for( fg_peak : all_peaks ) - pass 1

      // --- Pass 2: rescue borderline NORM peaks consistent with an accepted same-series line ---
      for( FgPeakDecision &dec : decisions )
      {
        if( !dec.borderline_norm || (dec.bg_cps <= 0.0) )
          continue;   // no way to form a ratio => stays rejected

        // Nearest PASS-1-accepted peak of the same series that also has a background counterpart.
        const FgPeakDecision *ref = nullptr;
        double bestDE = std::numeric_limits<double>::infinity();
        for( const FgPeakDecision &other : decisions )
        {
          if( (&other == &dec) || !other.include || other.borderline_norm )
            continue;
          if( (other.series != dec.series) || (other.bg_cps <= 0.0) || (other.fg_cps <= 0.0) )
            continue;
          const double de = fabs( other.mean - dec.mean );
          if( de < bestDE ) { bestDE = de; ref = &other; }
        }

        if( ref )
        {
          // If this peak's foreground/background ratio to the accepted line is consistent with them
          // being the same chain, keep it; if it is under-represented in the foreground relative to
          // that line, it's a background fluctuation and is dropped.
          const double R_fg = dec.fg_cps / ref->fg_cps;
          const double R_bg = dec.bg_cps / ref->bg_cps;
          const double consistency = (R_bg > 0.0) ? (R_fg / R_bg) : 0.0;
          dec.include = (consistency >= 0.5) && (consistency <= 3.0);
        }
        // else: no same-series accepted reference => stays rejected.
      }//for( pass 2 )

      // --- Emit results (+ optional debug dump) ---
      for( const FgPeakDecision &dec : decisions )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        cerr << "  fg E=" << dec.mean << " keV  amp=" << dec.peak->amplitude() << " cps=" << dec.fg_cps
             << (dec.is_seed ? "  [SEED]" : ((dec.series != NormSeries::None) ? "  [NORM]" : "        "))
             << (dec.borderline_norm ? " (borderline)" : "")
             << "  -> " << (dec.include ? "KEEP" : "drop (background)") << endl;
#endif
        if( dec.include )
        {
          filtered_peaks.push_back( dec.peak );
          if( dec.is_user )
            filtered_analysis_peaks.push_back( dec.peak );
        }
      }//for( emit results )

      all_peaks = std::move(filtered_peaks);
      user_analysis_peaks = std::move(filtered_analysis_peaks);
    }//if( options.nonBackgroundPeaksOnly )


    // Apply optional energy filtering
    const bool has_lower = options.lowerEnergy.has_value();
    const bool has_upper = options.upperEnergy.has_value();

    if( has_lower && has_upper && (*options.upperEnergy < *options.lowerEnergy) )
    {
      throw std::runtime_error( "upperEnergy (" + std::to_string(*options.upperEnergy)
                                + " keV) must be greater than or equal to lowerEnergy ("
                                + std::to_string(*options.lowerEnergy) + " keV)" );
    }

    if( has_lower || has_upper )
    {
      const double lower_energy = has_lower ? *options.lowerEnergy : -std::numeric_limits<double>::infinity();
      const double upper_energy = has_upper ? *options.upperEnergy : std::numeric_limits<double>::infinity();

      auto filter_peaks_in_range = [lower_energy, upper_energy]( const vector<shared_ptr<const PeakDef>> &peaks ){
        vector<shared_ptr<const PeakDef>> filtered;
        filtered.reserve( peaks.size() );

        for( const shared_ptr<const PeakDef> &peak : peaks )
        {
          if( peak && (peak->mean() >= lower_energy) && (peak->mean() <= upper_energy) )
            filtered.push_back( peak );
        }

        return filtered;
      };

      all_peaks = filter_peaks_in_range( all_peaks );
      user_analysis_peaks = filter_peaks_in_range( user_analysis_peaks );
    }

    DetectedPeakStatus result;
    result.peaks = std::move(all_peaks);
    result.analysis_peaks = std::move(user_analysis_peaks);

    return result;
  }


  PeakContinuum::OffsetType parse_continuum_type_str( const std::string &str );


  /** Everything needed to fit, and then commit, a `FitPeakOptions::explicitRoi` request; gathered on
   the GUI thread by `prepare_explicit_roi(...)`, so `fit_explicit_roi(...)` can run on any thread.
   */
  struct ExplicitRoiSetup
  {
    SpecUtils::SpectrumType specType = SpecUtils::SpectrumType::Foreground;
    bool doNotAddToAnalysisPeaks = false;

    shared_ptr<const SpecUtils::Measurement> data;
    shared_ptr<const DetectorPeakResponse> drf;
    shared_ptr<const PeakFitDetPrefs> fitPrefs;
    PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::CoarseResolutionType::Unknown;

    /** The first requested energy, as given; the peak nearest it is reported as `FitPeakStatus::fitPeak`. */
    double first_energy = 0.0;

    /** Requested starting means, sorted, with their sources (empty string for none). */
    vector<double> energies;
    vector<string> sources;

    /** Expected FWHM at each of `energies`. */
    vector<double> fwhms;

    double roi_lower = 0.0, roi_upper = 0.0;
    PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::Linear;

    /** Existing peaks to fit together with the requested ones, in the new ROI (as when a user
     double-clicks to add a peak to, or next to, an existing ROI); these keep their source, color,
     etc.  Sorted by mean. */
    vector<shared_ptr<const PeakDef>> cofit_peaks;

    /** Existing peaks removed, as a requested energy is at about the same energy; a requested peak
     without a source keeps the source, color, etc. of the replaced peak. */
    vector<shared_ptr<const PeakDef>> replaced_peaks;
  };//struct ExplicitRoiSetup


  /** Formats an energy for messages to the LLM, to at most two decimals, e.g., "742.81 keV". */
  string kev_str( const double energy )
  {
    char buffer[64] = { '\0' };
    snprintf( buffer, sizeof(buffer), "%.2f", energy );
    string answer = buffer;
    if( answer.find('.') != string::npos )
    {
      while( !answer.empty() && (answer.back() == '0') )
        answer.pop_back();
      if( !answer.empty() && (answer.back() == '.') )
        answer.pop_back();
    }
    return answer + " keV";
  }


  /** Groups `peaks` by ROI (i.e., shared continuum). */
  vector<vector<shared_ptr<const PeakDef>>> peaks_by_roi( const deque<shared_ptr<const PeakDef>> &peaks )
  {
    vector<vector<shared_ptr<const PeakDef>>> rois;
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      auto pos = std::find_if( begin(rois), end(rois),
                               [&p]( const auto &roi ){ return roi.front()->continuum() == p->continuum(); } );
      if( pos == end(rois) )
        rois.push_back( { p } );
      else
        pos->push_back( p );
    }
    return rois;
  }//peaks_by_roi(...)


  /** Works out what to do with the existing peaks of the spectrum, for the requested peaks/ROI of
   `setup`, filling in `setup.cofit_peaks` and `setup.replaced_peaks`, and widening the ROI of
   `setup` (for bounds not given) to include ROIs being combined.

   An existing ROI is combined with the new one (its peaks fit together with the requested ones) if:
   - no ROI bounds were given, and a requested peak is close to one of its peaks, or the ROI the
     data suggests for the requested peaks overlaps it substantially - the same criteria as when a
     user double-clicks to add a peak near an existing one (see `combine_peaks_to_roi(...)`).
   - ROI bounds were given, and at least half the existing ROI is inside them, or all its peaks are.
   An existing ROI that overlaps the new ROI, but isn't combined, is left as-is if none of its peaks
   are inside the new ROI (the ROIs just share an edge); otherwise it's an error, as what is wanted
   is ambiguous.

   An existing peak within 0.5 FWHM of a requested energy is replaced by the requested peak.

   Throws with a message intended for the LLM to act on.
   */
  void resolve_existing_rois( ExplicitRoiSetup &setup, const deque<shared_ptr<const PeakDef>> &existing,
                              const bool lower_given, const bool upper_given, const bool isHPGe )
  {
    setup.cofit_peaks.clear();
    setup.replaced_peaks.clear();

    const vector<vector<shared_ptr<const PeakDef>>> rois = peaks_by_roi( existing );
    vector<bool> combine( rois.size(), false );

    const auto roi_range = []( const vector<shared_ptr<const PeakDef>> &roi ) -> pair<double,double> {
      double lower = std::numeric_limits<double>::infinity(), upper = -lower;
      for( const shared_ptr<const PeakDef> &p : roi )
      {
        lower = std::min( lower, p->lowerX() );
        upper = std::max( upper, p->upperX() );
      }
      return { lower, upper };
    };

    const auto num_inside = [&setup]( const vector<shared_ptr<const PeakDef>> &roi ) -> size_t {
      size_t n = 0;
      for( const shared_ptr<const PeakDef> &p : roi )
        n += ((p->mean() > setup.roi_lower) && (p->mean() < setup.roi_upper));
      return n;
    };

    // Decide which ROIs to combine; with a bound not given, combining a ROI widens the new ROI, which
    //  may then reach further ROIs, so iterate until nothing changes.
    for( bool changed = true; changed; )
    {
      changed = false;
      for( size_t roi_index = 0; roi_index < rois.size(); ++roi_index )
      {
        if( combine[roi_index] )
          continue;

        const vector<shared_ptr<const PeakDef>> &roi = rois[roi_index];
        const pair<double,double> range = roi_range( roi );
        const double overlap = std::min( range.second, setup.roi_upper ) - std::max( range.first, setup.roi_lower );
        if( overlap <= 0.0 )
          continue;

        bool do_combine = false;
        if( lower_given && upper_given )
        {
          do_combine = (num_inside( roi ) == roi.size()) || (overlap >= 0.5*(range.second - range.first));
        }else
        {
          // The double-click criteria (`combine_peaks_to_roi(...)`): near in sigma, or the ROIs
          //  overlap by enough of the new ROI.
          const double overlap_frac = overlap / (setup.roi_upper - setup.roi_lower);
          for( const shared_ptr<const PeakDef> &p : roi )
          {
            for( const double energy : setup.energies )
            {
              const double nsigma = fabs( p->mean() - energy ) / p->sigma();
              if( isHPGe )
                do_combine |= ((((1.0 + fabs( p->mean() - energy )) / p->sigma()) < 8.5)
                               || ((overlap_frac > 0.35) && (nsigma <= 10.0)));
              else
                do_combine |= (((overlap_frac > 0.2) || (nsigma < 5.0)) && (nsigma < 10.0));
            }
          }//for( loop over ROI peaks )

          // A ROI with a peak inside the new ROI must be combined
          do_combine |= (num_inside( roi ) > 0);
        }//if( both bounds given ) / else

        if( !do_combine )
        {
          if( num_inside( roi ) )
          {
            string inside_str;
            for( const shared_ptr<const PeakDef> &p : roi )
              inside_str += (inside_str.empty() ? "" : ", ") + kev_str( p->mean() );

            throw runtime_error( "The ROI (" + kev_str(setup.roi_lower) + " to " + kev_str(setup.roi_upper)
                                 + ") contains part of the existing ROI from " + kev_str(range.first) + " to "
                                 + kev_str(range.second) + " (peaks at " + inside_str + "), but less than half of it;"
                                 " either widen the ROI to cover most of the existing ROI (its peaks are then fit"
                                 " along with yours), or narrow it to not contain any of its peaks." );
          }//if( an existing peak is inside the new ROI )

          continue;
        }//if( !do_combine )

        for( const shared_ptr<const PeakDef> &p : roi )
        {
          if( !p->gausPeak() )
            throw runtime_error( "The ROI (" + kev_str(setup.roi_lower) + " to " + kev_str(setup.roi_upper)
                                 + ") overlaps the data-defined peak covering " + kev_str(p->lowerX()) + " to "
                                 + kev_str(p->upperX()) + "; peaks can not be fit within the ROI of a data-defined peak." );
        }

        combine[roi_index] = true;
        changed = true;

        // Expand any bounds not given, to include the combined ROI
        if( !lower_given && (range.first < setup.roi_lower) )
          setup.roi_lower = range.first;
        if( !upper_given && (range.second > setup.roi_upper) )
          setup.roi_upper = range.second;
      }//for( loop over ROIs )
    }//for( until nothing changes )

    // The peaks of combined ROIs are fit along with the requested peaks, except those at about the
    //  same energy as a requested peak, which are replaced by it.
    for( size_t roi_index = 0; roi_index < rois.size(); ++roi_index )
    {
      if( !combine[roi_index] )
        continue;

      for( const shared_ptr<const PeakDef> &p : rois[roi_index] )
      {
        bool is_replaced = false;
        for( size_t i = 0; i < setup.energies.size(); ++i )
          is_replaced |= (fabs( p->mean() - setup.energies[i] ) < 0.5*std::max( p->fwhm(), setup.fwhms[i] ));

        (is_replaced ? setup.replaced_peaks : setup.cofit_peaks).push_back( p );
      }
    }//for( loop over ROIs )

    std::sort( begin(setup.cofit_peaks), end(setup.cofit_peaks), &PeakDef::lessThanByMeanShrdPtr );

    // Co-fit peaks must be inside the ROI (bounds given may cut a combined ROI)
    for( const shared_ptr<const PeakDef> &p : setup.cofit_peaks )
    {
      if( (p->mean() <= setup.roi_lower) || (p->mean() >= setup.roi_upper) )
        throw runtime_error( "The existing peak at " + kev_str(p->mean()) + " would be fit along with your peak(s),"
                             " but is outside the ROI (" + kev_str(setup.roi_lower) + " to " + kev_str(setup.roi_upper)
                             + "); include it in the ROI, or narrow the ROI to not overlap its ROI." );
    }

    // Co-fit peaks too close to each other, or a requested peak, can't be fit independently
    vector<pair<double,double>> means_fwhms;
    for( size_t i = 0; i < setup.energies.size(); ++i )
      means_fwhms.emplace_back( setup.energies[i], setup.fwhms[i] );
    for( const shared_ptr<const PeakDef> &p : setup.cofit_peaks )
      means_fwhms.emplace_back( p->mean(), p->fwhm() );
    std::sort( begin(means_fwhms), end(means_fwhms) );
    if( means_fwhms.size() > FitPeakOptions::ExplicitRoi::sm_max_peaks + 2 )
      throw runtime_error( "The ROI would have " + std::to_string(means_fwhms.size()) + " peaks, including "
                           + std::to_string(setup.cofit_peaks.size()) + " existing peaks; please use a narrower ROI." );
  }//resolve_existing_rois(...)


  /** Stage A of an explicit-ROI fit (GUI thread): gathers the data, validates the request, chooses
   any default ROI range and continuum type, and works out what to do with existing peaks.

   Throws with a message intended for the LLM to act on, on any problem with the request.
   */
  ExplicitRoiSetup prepare_explicit_roi( const FitPeakOptions &options, InterSpec *interspec )
  {
    assert( options.explicitRoi );
    if( !options.explicitRoi )
      throw logic_error( "prepare_explicit_roi: no explicit ROI specified" );

    const FitPeakOptions::ExplicitRoi &request = *options.explicitRoi;
    const string spec_desc = SpecUtils::descriptionText( options.specType );

    ExplicitRoiSetup setup;
    setup.specType = options.specType;
    setup.doNotAddToAnalysisPeaks = options.doNotAddToAnalysisPeaks;

    const shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
    if( !meas )
      throw runtime_error( "No measurement loaded for " + spec_desc + " spectrum" );

    const set<int> sample_nums = interspec->displayedSamples( options.specType );
    if( sample_nums.empty() )
      throw runtime_error( "No samples displayed for " + spec_desc + " spectrum" );

    setup.data = interspec->displayedHistogram( options.specType );
    const shared_ptr<const SpecUtils::Measurement> &data = setup.data;
    if( !data || (data->num_gamma_channels() < 16) )
      throw runtime_error( "No spectrum displayed for " + spec_desc + " spectrum" );

    setup.drf = meas->detector();
    setup.fitPrefs = meas->peakFitDetPrefs();
    if( !setup.fitPrefs && setup.drf )
      setup.fitPrefs = setup.drf->peakFitDetPrefs();
    if( !setup.fitPrefs )
      setup.fitPrefs = make_shared<PeakFitDetPrefs>();
    // Classify from the data, not the file metadata, like the single-energy search does
    setup.det_type = PeakFitUtils::effective_det_type( setup.fitPrefs, data, nullptr );
    const bool isHPGe = (setup.det_type == PeakFitUtils::CoarseResolutionType::High)
                        || ((setup.det_type == PeakFitUtils::CoarseResolutionType::Unknown)
                            && PeakFitUtils::is_likely_high_res( interspec ));

    // The analysis peaks of this spectrum; we will not wait for the automated search, so its peaks
    //  are used for estimating resolution only if they are already available.
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks( sample_nums );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );

    auto resolution_peaks = make_shared<deque<shared_ptr<const PeakDef>>>();
    for( const shared_ptr<const PeakDef> &p : *meas_peaks )
    {
      if( p->gausPeak() )
        resolution_peaks->push_back( p );
    }
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = meas->automatedSearchPeaks( sample_nums );
    if( auto_peaks )
      resolution_peaks->insert( end(*resolution_peaks), begin(*auto_peaks), end(*auto_peaks) );

    const PeakFitUtils::CoarseResolutionType fwhm_det_type = isHPGe ? PeakFitUtils::CoarseResolutionType::High
                                                                    : setup.det_type;
    const std::function<float(double)> fwhm_fcn = make_expected_fwhm_fcn( setup.drf, resolution_peaks, data, fwhm_det_type );

    // Validate, and sort, the requested energies
    if( request.energies.empty() )
      throw runtime_error( "At least one energy must be specified." );
    if( request.energies.size() > FitPeakOptions::ExplicitRoi::sm_max_peaks )
      throw runtime_error( "At most " + std::to_string(FitPeakOptions::ExplicitRoi::sm_max_peaks)
                           + " peaks may be fit in one ROI; " + std::to_string(request.energies.size())
                           + " energies were given." );
    if( !request.sources.empty() && (request.sources.size() != request.energies.size()) )
      throw runtime_error( "The number of sources (" + std::to_string(request.sources.size())
                           + ") must match the number of energies (" + std::to_string(request.energies.size()) + ")." );

    const double spec_lower = data->gamma_energy_min(), spec_upper = data->gamma_energy_max();

    setup.first_energy = request.energies.front();
    vector<pair<double,string>> energy_sources;
    for( size_t i = 0; i < request.energies.size(); ++i )
    {
      const double energy = request.energies[i];
      if( IsNan(energy) || IsInf(energy) || (energy <= spec_lower) || (energy >= spec_upper) )
        throw runtime_error( "Energy " + kev_str(energy) + " is outside the spectrum energy range ("
                             + kev_str(spec_lower) + " to " + kev_str(spec_upper) + ")." );
      energy_sources.emplace_back( energy, request.sources.empty() ? string() : request.sources[i] );
    }

    std::stable_sort( begin(energy_sources), end(energy_sources),
                      []( const pair<double,string> &lhs, const pair<double,string> &rhs ){ return lhs.first < rhs.first; } );

    for( const pair<double,string> &es : energy_sources )
    {
      double fwhm = fwhm_fcn( es.first );
      if( (fwhm <= 0.0) || IsNan(fwhm) || IsInf(fwhm) )
        throw runtime_error( "Could not estimate the peak width at " + kev_str(es.first) + "." );

      // Keep the starting width within what the fitter will allow
      float min_sigma = 0.0f, max_sigma = 0.0f;
      expected_peak_width_limits( static_cast<float>(es.first), setup.det_type, data, min_sigma, max_sigma );
      if( (min_sigma > 0.0f) && (max_sigma > min_sigma) )
        fwhm = std::clamp( fwhm, PhysicalUnits::fwhm_nsigma*1.05*min_sigma, PhysicalUnits::fwhm_nsigma*0.95*max_sigma );

      setup.energies.push_back( es.first );
      setup.sources.push_back( es.second );
      setup.fwhms.push_back( fwhm );
    }//for( const pair<double,string> &es : energy_sources )

    const size_t npeaks = setup.energies.size();
    for( size_t i = 1; i < npeaks; ++i )
    {
      const double min_sep = 0.5 * 0.5 * (setup.fwhms[i-1] + setup.fwhms[i]);
      if( (setup.energies[i] - setup.energies[i-1]) < min_sep )
        throw runtime_error( "Energies " + kev_str(setup.energies[i-1]) + " and " + kev_str(setup.energies[i])
                             + " are within half a FWHM of each other, which is too close to fit as separate peaks." );
    }

    // Test each source string is valid, before doing any fitting
    for( size_t i = 0; i < npeaks; ++i )
    {
      PeakDef test_peak( setup.energies[i], setup.fwhms[i] / PhysicalUnits::fwhm_nsigma, 100.0 );
      map<string,Wt::WColor> unused_colors;
      if( !assign_source_and_color( test_peak, setup.sources[i], {}, nullptr, unused_colors ) )
        throw runtime_error( "Source '" + setup.sources[i] + "' (for the peak at " + kev_str(setup.energies[i])
                             + ") is not a valid nuclide, x-ray element, or reaction." );
    }

    if( request.continuumType.has_value() )
    {
      try
      {
        setup.continuum_type = parse_continuum_type_str( *request.continuumType );
      }catch( std::exception & )
      {
        throw runtime_error( "Invalid continuumType '" + *request.continuumType + "'; valid values are NoOffset,"
                             " Constant, Linear, Quadratic, Cubic, FlatStep, LinearStep, BiLinearStep, FlatStepCDF,"
                             " LinearStepCDF, and BiLinearStepCDF." );
      }

      if( setup.continuum_type == PeakContinuum::OffsetType::External )
        throw runtime_error( "The External continuumType can not be used when adding a peak." );
    }//if( request.continuumType )

    // ROI bounds not given: the ROI limit found from the data for the outer peaks (like for a
    //  single peak ROI), but at least 1.5 FWHM past them; this is then widened to include any
    //  existing ROIs combined with it.
    setup.roi_lower = request.lowerEnergy.value_or( 0.0 );
    setup.roi_upper = request.upperEnergy.value_or( 0.0 );
    if( !request.lowerEnergy.has_value() )
    {
      const PeakDef low_peak( setup.energies.front(), setup.fwhms.front() / PhysicalUnits::fwhm_nsigma, 100.0 );
      double lower = 0.0, upper = 0.0;
      findROIEnergyLimits( lower, upper, low_peak, data, isHPGe );
      setup.roi_lower = std::min( lower, setup.energies.front() - 1.5*setup.fwhms.front() );
    }

    if( !request.upperEnergy.has_value() )
    {
      const PeakDef high_peak( setup.energies.back(), setup.fwhms.back() / PhysicalUnits::fwhm_nsigma, 100.0 );
      double lower = 0.0, upper = 0.0;
      findROIEnergyLimits( lower, upper, high_peak, data, isHPGe );
      setup.roi_upper = std::max( upper, setup.energies.back() + 1.5*setup.fwhms.back() );
    }

    if( IsNan(setup.roi_lower) || IsInf(setup.roi_lower) || IsNan(setup.roi_upper) || IsInf(setup.roi_upper) )
      throw runtime_error( "Invalid ROI energy range." );

    if( setup.roi_lower >= setup.roi_upper )
      throw runtime_error( "roiLowerEnergy (" + kev_str(setup.roi_lower) + ") must be less than roiUpperEnergy ("
                           + kev_str(setup.roi_upper) + ")." );

    for( const double energy : setup.energies )
    {
      if( (energy <= setup.roi_lower) || (energy >= setup.roi_upper) )
        throw runtime_error( "Energy " + kev_str(energy) + " is not inside the ROI (" + kev_str(setup.roi_lower)
                             + " to " + kev_str(setup.roi_upper) + "); each energy must be within roiLowerEnergy and roiUpperEnergy." );
    }

    resolve_existing_rois( setup, *meas_peaks, request.lowerEnergy.has_value(), request.upperEnergy.has_value(), isHPGe );

    // Defaulted bounds stay within the spectrum
    if( !request.lowerEnergy.has_value() )
      setup.roi_lower = std::max( setup.roi_lower, static_cast<double>( data->gamma_channel_lower(1) ) );
    if( !request.upperEnergy.has_value() )
      setup.roi_upper = std::min( setup.roi_upper, static_cast<double>( data->gamma_channel_upper(data->num_gamma_channels() - 2) ) );

    if( (setup.roi_lower < spec_lower) || (setup.roi_upper > spec_upper) )
      throw runtime_error( "The ROI (" + kev_str(setup.roi_lower) + " to " + kev_str(setup.roi_upper)
                           + ") extends past the spectrum energy range (" + kev_str(spec_lower) + " to " + kev_str(spec_upper) + ")." );

    // Without a requested continuum: the default for the number of peaks, or a combined ROI's
    //  continuum if it is a higher order (as `fit_peak_for_user_click_LM(...)` does)
    const size_t num_roi_peaks = npeaks + setup.cofit_peaks.size();
    if( !request.continuumType.has_value() )
    {
      setup.continuum_type = (num_roi_peaks < (isHPGe ? 3 : 2)) ? PeakContinuum::OffsetType::Linear
                                                                : PeakContinuum::OffsetType::Quadratic;
      for( const vector<shared_ptr<const PeakDef>> *peaks : { &setup.cofit_peaks, &setup.replaced_peaks } )
      {
        for( const shared_ptr<const PeakDef> &p : *peaks )
        {
          const PeakContinuum::OffsetType type = p->continuum()->type();
          if( (type != PeakContinuum::OffsetType::External) && (type > setup.continuum_type) )
            setup.continuum_type = type;
        }
      }
    }//if( !request.continuumType.has_value() )

    const size_t lower_channel = data->find_gamma_channel( static_cast<float>(setup.roi_lower) );
    const size_t upper_channel = data->find_gamma_channel( static_cast<float>(setup.roi_upper) );
    const size_t num_roi_channels = (upper_channel > lower_channel) ? (1 + upper_channel - lower_channel) : 1;
    const size_t num_fit_pars = PeakContinuum::num_parameters( setup.continuum_type ) + 3*num_roi_peaks;
    if( num_roi_channels <= num_fit_pars )
      throw runtime_error( "The ROI (" + kev_str(setup.roi_lower) + " to " + kev_str(setup.roi_upper) + ") only has "
                           + std::to_string(num_roi_channels) + " channels, which is too few to fit "
                           + std::to_string(num_roi_peaks) + " peak(s) with a " + PeakContinuum::offset_type_name(setup.continuum_type)
                           + " continuum; please widen the ROI." );

    return setup;
  }//prepare_explicit_roi(...)


  /** The fit results of an explicit ROI. */
  struct ExplicitRoiFitResult
  {
    /** The fit requested peaks, in the same order as `ExplicitRoiSetup::energies`. */
    vector<shared_ptr<const PeakDef>> requested_peaks;

    /** The fit co-fit (existing) peaks, in the same order as `ExplicitRoiSetup::cofit_peaks`. */
    vector<shared_ptr<const PeakDef>> cofit_peaks;

    vector<string> warnings;
  };//struct ExplicitRoiFitResult


  /** Stage B of an explicit-ROI fit: fits the requested peaks, along with any existing peaks being
   combined into the ROI.

   Does not touch any session state, so may be called from a background thread.
   Throws on failure.
   */
  ExplicitRoiFitResult fit_explicit_roi( const ExplicitRoiSetup &setup )
  {
    const shared_ptr<const SpecUtils::Measurement> &data = setup.data;
    assert( data );
    assert( !setup.energies.empty() );
    assert( setup.energies.size() == setup.fwhms.size() );

    ExplicitRoiFitResult result;

    // One continuum shared by all the peaks; start it as a straight line between the ROI edges.
    //  The starting areas are estimated using this line, before changing the continuum type (a
    //  CDF-step continuum can't be integrated without its peaks).
    auto continuum = make_shared<PeakContinuum>();
    continuum->calc_linear_continuum_eqn( data, setup.roi_lower, setup.roi_lower, setup.roi_upper, 2, 2 );

    vector<shared_ptr<PeakDef>> new_peaks;
    for( size_t i = 0; i < setup.energies.size(); ++i )
    {
      const double mean = setup.energies[i];
      const double fwhm = setup.fwhms[i];

      // Initial area from the counts above the continuum within +-0.5 FWHM (76.1% of a Gaussian);
      //  the amplitudes are solved for linearly by the fit, so this just needs to be reasonable.
      const double x0 = mean - 0.5*fwhm, x1 = mean + 0.5*fwhm;
      const double data_counts = data->gamma_integral( static_cast<float>(x0), static_cast<float>(x1) );
      const double cont_counts = continuum->offset_integral( x0, x1, data, vector<shared_ptr<const PeakDef>>{} );
      const double area = std::max( 1.0, (data_counts - cont_counts) / 0.761 );

      auto peak = make_shared<PeakDef>( mean, fwhm / PhysicalUnits::fwhm_nsigma, area );
      peak->setContinuum( continuum );
      new_peaks.push_back( peak );
    }//for( loop over energies )

    // Skew and FWHM method from the peak-fit preferences, like a user fit uses
    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
    assert( setup.fitPrefs );
    if( setup.fitPrefs )
      apply_fit_prefs_to_peaks( new_peaks, data, setup.drf, *setup.fitPrefs, fit_options );

    // Existing peaks being combined into the ROI keep their parameters as starting values, and their
    //  skew type is made to match (the fit uses one skew type per ROI)
    vector<shared_ptr<PeakDef>> cofit_peaks;
    for( const shared_ptr<const PeakDef> &p : setup.cofit_peaks )
    {
      auto peak = make_shared<PeakDef>( *p );
      peak->setContinuum( continuum );
      if( setup.fitPrefs )
        apply_prefs_skew_to_peak( *peak, *setup.fitPrefs, data,
                                  setup.drf ? setup.drf->peakFitDetPrefs().get() : nullptr );
      cofit_peaks.push_back( peak );
    }

    continuum->setType( setup.continuum_type );

    vector<shared_ptr<const PeakDef>> input_peaks( begin(new_peaks), end(new_peaks) );
    input_peaks.insert( end(input_peaks), begin(cofit_peaks), end(cofit_peaks) );
    std::sort( begin(input_peaks), end(input_peaks), &PeakDef::lessThanByMeanShrdPtr );

    vector<shared_ptr<const PeakDef>> fit_peaks = PeakFitLM::fit_peaks_in_roi_LM( input_peaks, data, setup.det_type, fit_options );

    if( fit_peaks.size() != input_peaks.size() )
      throw runtime_error( "Fitting the peaks in the ROI failed." );

    // The fit returns peaks sorted by mean; match them back to the input peaks by order, which the
    //  fit keeps unless peaks cross.
    std::sort( begin(fit_peaks), end(fit_peaks), &PeakDef::lessThanByMeanShrdPtr );

    for( size_t i = 0; i < fit_peaks.size(); ++i )
    {
      const shared_ptr<const PeakDef> &peak = fit_peaks[i];
      if( !peak->gausPeak() || IsNan(peak->mean()) || IsNan(peak->amplitude()) )
        throw runtime_error( "Fitting the peaks in the ROI gave invalid results." );

      if( i && (fabs( peak->mean() - fit_peaks[i-1]->mean() ) < 0.5*std::min( peak->fwhm(), fit_peaks[i-1]->fwhm() )) )
        throw runtime_error( "The peaks started at " + kev_str(input_peaks[i-1]->mean()) + " and " + kev_str(input_peaks[i]->mean())
                             + " fit to nearly the same energy (" + kev_str(peak->mean()) + "), so the data likely only"
                             " supports one peak there; please remove one of the energies." );

      const auto pos = std::find( begin(new_peaks), end(new_peaks), input_peaks[i] );
      if( pos != end(new_peaks) )
      {
        const size_t index = static_cast<size_t>( pos - begin(new_peaks) );
        if( result.requested_peaks.empty() )
          result.requested_peaks.resize( new_peaks.size() );
        result.requested_peaks[index] = peak;

        if( fabs(peak->mean() - setup.energies[index]) > setup.fwhms[index] )
          result.warnings.push_back( "The peak started at " + kev_str(setup.energies[index]) + " moved to "
                                     + kev_str(peak->mean()) + " (more than one FWHM); check if this is the peak you intended." );

        const double nsigma = (peak->amplitudeUncert() > 0.0) ? (peak->amplitude() / peak->amplitudeUncert()) : 0.0;
        if( peak->amplitude() <= 0.0 )
          result.warnings.push_back( "The peak at " + kev_str(peak->mean()) + " fit to a negative area, so there is"
                                     " likely no peak there; it was marked to not be used for activity or relative"
                                     " efficiency fits - consider deleting it." );
        else if( nsigma < 2.0 )
          result.warnings.push_back( "The peak at " + kev_str(peak->mean()) + " is not statistically significant ("
                                     + SpecUtils::printCompact(nsigma, 2) + " sigma); consider deleting it." );
      }else
      {
        result.cofit_peaks.push_back( peak );
      }
    }//for( loop over fit peaks )

    assert( result.requested_peaks.size() == new_peaks.size() );
    assert( result.cofit_peaks.size() == setup.cofit_peaks.size() );

#if( PERFORM_DEVELOPER_CHECKS )
    {
      const shared_ptr<const PeakContinuum> fit_cont = fit_peaks.front()->continuum();
      for( const shared_ptr<const PeakDef> &p : fit_peaks )
      {
        if( p->continuum() != fit_cont )
          log_developer_error( __func__, "Explicit ROI fit peaks do not share a continuum." );
      }
      if( (fabs(fit_cont->lowerEnergy() - setup.roi_lower) > 0.01)
         || (fabs(fit_cont->upperEnergy() - setup.roi_upper) > 0.01)
         || (fit_cont->type() != setup.continuum_type) )
      {
        log_developer_error( __func__, "Explicit ROI fit did not preserve ROI range or continuum type." );
      }

      for( const shared_ptr<const PeakDef> &p : result.cofit_peaks )
      {
        // The fit carries over sources (etc.) by closest mean (`inheritUserSelectedOptions`); check the
        //  co-fit peaks kept theirs
        const auto orig = std::find_if( begin(setup.cofit_peaks), end(setup.cofit_peaks),
          [&p]( const shared_ptr<const PeakDef> &o ){ return fabs( o->mean() - p->mean() ) < o->fwhm(); } );
        if( (orig != end(setup.cofit_peaks)) && ((*orig)->sourceName() != p->sourceName()) )
          log_developer_error( __func__, ("Co-fit peak at " + kev_str(p->mean()) + " changed source.").c_str() );
      }
    }
#endif

    return result;
  }//fit_explicit_roi(...)


  /** Stage C of an explicit-ROI fit (GUI thread): assigns sources, and unless
   `setup.doNotAddToAnalysisPeaks`, updates the analysis peaks, as one undo/redo step.

   The current analysis peaks are re-checked, so other tool calls (e.g., run in parallel) committed
   since `setup` was made are fine, unless they changed the peaks this fit uses or replaces.

   Throws if the existing peaks this fit depends on changed, or a source can't be assigned.
   */
  FitPeakStatus commit_explicit_roi( const ExplicitRoiSetup &setup, const ExplicitRoiFitResult &fit,
                                     InterSpec *interspec )
  {
    assert( interspec );
    assert( fit.requested_peaks.size() == setup.energies.size() );
    assert( fit.cofit_peaks.size() == setup.cofit_peaks.size() );

    const shared_ptr<SpecMeas> meas = interspec->measurment( setup.specType );
    const set<int> sample_nums = interspec->displayedSamples( setup.specType );
    const shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas ? meas->peaks( sample_nums ) : nullptr;
    if( !meas_peaks )
      throw runtime_error( "The spectrum changed while fitting; please retry." );

    // The existing peaks this fit uses or replaces must still be there, and no other peak may have
    //  been added inside the ROI.
    vector<shared_ptr<const PeakDef>> to_remove = setup.cofit_peaks;
    to_remove.insert( end(to_remove), begin(setup.replaced_peaks), end(setup.replaced_peaks) );
    for( const shared_ptr<const PeakDef> &p : to_remove )
    {
      if( std::find( begin(*meas_peaks), end(*meas_peaks), p ) == end(*meas_peaks) )
        throw runtime_error( "The analysis peak at " + kev_str(p->mean()) + " was changed while fitting; please retry." );
    }
    for( const shared_ptr<const PeakDef> &p : *meas_peaks )
    {
      if( (p->mean() > setup.roi_lower) && (p->mean() < setup.roi_upper)
         && (std::find( begin(to_remove), end(to_remove), p ) == end(to_remove)) )
        throw runtime_error( "A peak at " + kev_str(p->mean()) + " was added inside the ROI while fitting; please retry." );
    }

    // The other peaks of the spectrum, for matching source colors
    vector<shared_ptr<const PeakDef>> other_peaks;
    for( const shared_ptr<const PeakDef> &p : *meas_peaks )
    {
      if( std::find( begin(setup.replaced_peaks), end(setup.replaced_peaks), p ) == end(setup.replaced_peaks) )
        other_peaks.push_back( p );
    }

    map<string,Wt::WColor> pending_colors;  // So each new source of the ROI gets its own color

    vector<shared_ptr<const PeakDef>> new_peaks;
    vector<bool> inherited( setup.replaced_peaks.size(), false );
    for( size_t i = 0; i < fit.requested_peaks.size(); ++i )
    {
      auto peak = make_shared<PeakDef>( *fit.requested_peaks[i] );
      peak->clearSources();
      peak->setUserLabel( "" );

      if( !setup.sources[i].empty() )
      {
        if( !assign_source_and_color( *peak, setup.sources[i], other_peaks, interspec, pending_colors ) )
          throw runtime_error( "Source '" + setup.sources[i] + "' is not valid." );
      }else
      {
        // Without a requested source, keep the source (and color, label, and use-for flags) of a
        //  replaced peak at about the same energy.
        size_t best = setup.replaced_peaks.size();
        for( size_t j = 0; j < setup.replaced_peaks.size(); ++j )
        {
          const shared_ptr<const PeakDef> &old = setup.replaced_peaks[j];
          const double diff = fabs( old->mean() - peak->mean() );
          if( !inherited[j] && (diff < peak->fwhm())
             && ((best >= setup.replaced_peaks.size()) || (diff < fabs(setup.replaced_peaks[best]->mean() - peak->mean()))) )
          {
            best = j;
          }
        }

        if( best < setup.replaced_peaks.size() )
        {
          inherited[best] = true;
          const PeakDef &old = *setup.replaced_peaks[best];
          if( old.parentNuclide() )
            peak->setNuclearTransition( old.parentNuclide(), old.nuclearTransition(),
                                        old.decayParticleIndex(), old.sourceGammaType() );
          else if( old.xrayElement() )
            peak->setXray( old.xrayElement(), old.xrayEnergy() );
          else if( old.reaction() )
            peak->setReaction( old.reaction(), old.reactionEnergy(), old.sourceGammaType() );

          peak->setLineColor( old.lineColor() );
          peak->setUserLabel( old.userLabel() );
          peak->useForEnergyCalibration( old.useForEnergyCalibrationUserPreference() );
          peak->useForShieldingSourceFit( old.useForShieldingSourceFitUserPreference() );
          peak->useForManualRelEff( old.useForManualRelEffUserPreference() );
        }//if( found a replaced peak to inherit from )
      }//if( source requested ) / else

      // A peak with no area shouldn't feed activity or relative efficiency fits
      if( peak->amplitude() <= 0.0 )
      {
        peak->useForShieldingSourceFit( false );
        peak->useForManualRelEff( false );
      }

      new_peaks.push_back( peak );
    }//for( loop over requested peaks )

    // The co-fit peaks keep the source, color, label, and use-for flags they had
    vector<shared_ptr<const PeakDef>> refit_peaks;
    for( size_t i = 0; i < fit.cofit_peaks.size(); ++i )
    {
      auto peak = make_shared<PeakDef>( *fit.cofit_peaks[i] );
      peak->inheritUserSelectedOptions( *setup.cofit_peaks[i], false );
      // `inheritUserSelectedOptions` copies the fit-for flags too, which may have been changed for the fit
      for( PeakDef::CoefficientType t = PeakDef::CoefficientType(0); t < PeakDef::NumCoefficientTypes;
          t = PeakDef::CoefficientType(t+1) )
      {
        peak->setFitFor( t, fit.cofit_peaks[i]->fitFor(t) );
      }
      refit_peaks.push_back( peak );
    }

    vector<shared_ptr<const PeakDef>> to_add = new_peaks;
    to_add.insert( end(to_add), begin(refit_peaks), end(refit_peaks) );

    // All the ROI peaks share a continuum
    assert( std::all_of( begin(to_add), end(to_add),
                         [&to_add]( const auto &p ){ return p->continuum() == to_add.front()->continuum(); } ) );

    if( !setup.doNotAddToAnalysisPeaks )
    {
      const std::shared_ptr<PeakModel> pmodel = interspec->peakModel();
      if( !pmodel )
        throw runtime_error( "No peak model available." );

      UndoRedoManager::PeakModelChange peak_undo_creator;

      if( setup.specType == SpecUtils::SpectrumType::Foreground )
      {
        pmodel->removePeaks( to_remove );
        pmodel->addPeaks( to_add );
      }else
      {
        deque<shared_ptr<const PeakDef>> peak_deque = *meas_peaks;
        for( const shared_ptr<const PeakDef> &old : to_remove )
        {
          const auto pos = std::find( begin(peak_deque), end(peak_deque), old );
          assert( pos != end(peak_deque) );
          if( pos != end(peak_deque) )
            peak_deque.erase( pos );
        }

        peak_deque.insert( end(peak_deque), begin(to_add), end(to_add) );
        std::sort( begin(peak_deque), end(peak_deque), &PeakDef::lessThanByMeanShrdPtr );

        pmodel->setPeaks( peak_deque, setup.specType );
      }//if( foreground ) / else

#if( PERFORM_DEVELOPER_CHECKS )
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> final_peaks = meas->peaks( sample_nums );
      size_t num_in_new_roi = 0;
      for( const shared_ptr<const PeakDef> &p : (final_peaks ? *final_peaks : deque<shared_ptr<const PeakDef>>{}) )
      {
        if( p->continuum() == to_add.front()->continuum() )
        {
          num_in_new_roi += 1;
          continue;
        }
        // Other ROIs may overlap the new ROIs edges, but not have peaks inside it
        if( (p->mean() > setup.roi_lower) && (p->mean() < setup.roi_upper) )
          log_developer_error( __func__, ("Another ROI has a peak inside the explicit ROI, at " + kev_str(p->mean())).c_str() );
      }
      // Foreground peaks are copied as they are added to the PeakModel, but keep their continuum
      if( num_in_new_roi != to_add.size() )
        log_developer_error( __func__, "Not all explicit ROI peaks made it into the analysis peaks." );
#endif
    }//if( !setup.doNotAddToAnalysisPeaks )

    // Report the peaks as they are now in the analysis peaks (the PeakModel copies foreground peaks)
    vector<shared_ptr<const PeakDef>> reported = to_add;
    std::sort( begin(reported), end(reported), &PeakDef::lessThanByMeanShrdPtr );
    if( !setup.doNotAddToAnalysisPeaks )
    {
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> final_peaks = meas->peaks( sample_nums );
      if( final_peaks )
      {
        vector<shared_ptr<const PeakDef>> in_model;
        for( const shared_ptr<const PeakDef> &p : *final_peaks )
        {
          if( p->continuum() == to_add.front()->continuum() )
            in_model.push_back( p );
        }
        if( in_model.size() == to_add.size() )
          reported = in_model;
      }
    }//if( peaks were added )

    FitPeakStatus status;
    status.explicitRoi = true;
    status.peaksInRoi = reported;
    if( !setup.doNotAddToAnalysisPeaks )
      status.removedPeaks = setup.replaced_peaks;
    status.warnings = fit.warnings;

    // `fitPeak` is the peak nearest the first requested energy (as given, before sorting)
    const auto first_pos = std::find( begin(setup.energies), end(setup.energies), setup.first_energy );
    assert( first_pos != end(setup.energies) );
    const shared_ptr<const PeakDef> first_peak = new_peaks[ (first_pos != end(setup.energies)) ? (first_pos - begin(setup.energies)) : 0 ];
    for( const shared_ptr<const PeakDef> &p : reported )
    {
      if( !status.fitPeak || (fabs(p->mean() - first_peak->mean()) < fabs(status.fitPeak->mean() - first_peak->mean())) )
        status.fitPeak = p;
    }

    return status;
  }//commit_explicit_roi(...)


  ExplicitRoiPreparedFit prepare_explicit_roi_fit( const FitPeakOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error( "No InterSpec session available" );

    ExplicitRoiPreparedFit prepared;
    prepared.setup = make_shared<ExplicitRoiSetup>( prepare_explicit_roi( options, interspec ) );
    return prepared;
  }//prepare_explicit_roi_fit(...)


  void run_explicit_roi_fit( ExplicitRoiPreparedFit &prepared )
  {
    assert( prepared.setup );
    if( !prepared.setup )
      throw std::logic_error( "run_explicit_roi_fit: not prepared" );
    prepared.fit = make_shared<ExplicitRoiFitResult>( fit_explicit_roi( *prepared.setup ) );
  }//run_explicit_roi_fit(...)


  FitPeakStatus commit_explicit_roi_fit( const ExplicitRoiPreparedFit &prepared, InterSpec *interspec )
  {
    assert( prepared.setup && prepared.fit );
    if( !prepared.setup || !prepared.fit )
      throw std::logic_error( "commit_explicit_roi_fit: not fit" );
    if( !interspec )
      throw std::runtime_error( "No InterSpec session available" );
    return commit_explicit_roi( *prepared.setup, *prepared.fit, interspec );
  }//commit_explicit_roi_fit(...)


  FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");

    if( options.explicitRoi )
    {
      ExplicitRoiPreparedFit prepared = prepare_explicit_roi_fit( options, interspec );
      run_explicit_roi_fit( prepared );
      return commit_explicit_roi_fit( prepared, interspec );
    }//if( options.explicitRoi )

    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram(options.specType);
    if (!data) {
      throw std::runtime_error("No spectrum displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<const DetectorPeakResponse> det = meas->detector();
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks(sample_nums);
    assert( meas_peaks );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );

    shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_search_peaks = meas->automatedSearchPeaks(sample_nums);

    shared_ptr<const PeakFitDetPrefs> fitPrefs = meas->peakFitDetPrefs();

    vector<shared_ptr<const PeakDef>> origPeaks;

    for( const shared_ptr<const PeakDef> &p : *meas_peaks )
    {
      //Avoid fitting a peak in the same area a data defined peak is.
      if( p->gausPeak() )
        origPeaks.push_back( p );
      else if( (options.energy >= p->lowerX()) && (options.energy <= p->upperX()) )
        throw runtime_error( "Can not fit a peak within the ROI of a data-defined peak." );
    }//if( pmodel->peaks() )

    double pixelPerKev = -1.0; //This triggers an "automed" peak fit, which has higher thresholds for keeping peak.
    pair<vector<shared_ptr<const PeakDef>>, vector<shared_ptr<const PeakDef>>> foundPeaks;
    string reject_reason;
    foundPeaks = searchForPeakFromUser( options.energy, pixelPerKev, data, origPeaks, det, auto_search_peaks, fitPrefs,
                                        nullptr, &reject_reason );
    // The peaks become the user's: fit them like any other (the automated fit returns its chi2 decision fits).
    if( !foundPeaks.first.empty() )
      foundPeaks.first = ExperimentalAutomatedPeakSearch::refit_sparse_rois( foundPeaks.first, data,
                                                PeakFitUtils::effective_det_type( fitPrefs, data, nullptr ), {} );

    vector<shared_ptr<const PeakDef>> &peaks_to_add_in = foundPeaks.first;
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove = foundPeaks.second;
    
    if( peaks_to_add_in.empty() )
    {
      // Check if the requested energy is within the ROI of any existing peak
      bool within_existing_roi = false;
      for( const shared_ptr<const PeakDef> &p : origPeaks )
      {
        if( (options.energy >= p->lowerX()) && (options.energy <= p->upperX()) )
        {
          within_existing_roi = true;
          break;
        }
      }

      // What the data do show there, so the caller gets a measurement, or at least a reason.
      string fwhm_source;
      const std::function<float(double)> fwhm_fcn = fwhm_fcn_for_peak_fit( resolve_det_type( interspec, meas, data, fitPrefs ),
                                                            data, det, origPeaks, auto_search_peaks, fwhm_source );
      UnfitPeakDiagnostic diag = diagnose_unfit_peak( options.energy, reject_reason, fwhm_fcn( options.energy ),
                                       fwhm_source, data, origPeaks, det, auto_search_peaks, fitPrefs );
      refine_unfit_diagnostic_fwhm( diag, options.energy, interspec, options.specType, data, det, fitPrefs );

      if( options.doNotAddToAnalysisPeaks )
      {
        FitPeakStatus result;
        result.significant = false;
        result.fitPeak = diag.candidate;
        result.peaksInRoi = diag.candidateRoiPeaks;
        result.diagnostic = diag;
        return result;
      }//if( only probing the energy )

      throw runtime_error( unfit_peak_error_message( options.energy, within_existing_roi, diag ) );
    }
    
    const string source = options.source.has_value() ? options.source.value() : ""s;
    const shared_ptr<const PeakDef> fit_peak = assign_source_to_user_fit_peak( peaks_to_add_in,
                                                options.energy, source, origPeaks, interspec );

    const std::shared_ptr<PeakModel> pmodel = interspec->peakModel();
    assert( pmodel );
    if( !options.doNotAddToAnalysisPeaks && pmodel )
    {
      if( options.specType == SpecUtils::SpectrumType::Foreground )
      {
        pmodel->removePeaks( peaks_to_remove );
        pmodel->addPeaks( peaks_to_add_in );
      }else
      {
        // It would also be valid to handle the Foreground in this logic block, I think
        shared_ptr<deque<shared_ptr<const PeakDef>>> new_deque = make_shared<deque<shared_ptr<const PeakDef>>>(*meas_peaks);
        deque<shared_ptr<const PeakDef>> &peak_deque = *new_deque;

        for( shared_ptr<const PeakDef> old : peaks_to_remove )
        {
          auto pos = std::find( begin(peak_deque), end(peak_deque), old );
          assert( pos != end(peak_deque) );
          if( pos != end(peak_deque) )
            peak_deque.erase( pos );
        }

        for( shared_ptr<const PeakDef> new_peak : peaks_to_add_in )
          peak_deque.push_back( new_peak );
        std::sort( begin(peak_deque), end(peak_deque), &PeakDef::lessThanByMeanShrdPtr );

        pmodel->setPeaks( peak_deque, options.specType );
      }
    }//if( !options.doNotAddToAnalysisPeaks )
    
    
    
    FitPeakStatus result;

    result.fitPeak = fit_peak;
    result.peaksInRoi = peaks_to_add_in;

    return result;
  }//FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );


  void fit_user_peak_async( const FitPeakOptions &options, InterSpec *interspec,
                             FitPeakCallback callback )
  {
    assert( callback );
    assert( interspec );
    assert( wApp );

    if( !interspec )
    {
      callback( string("No InterSpec session available") );
      return;
    }

#if( BUILD_AS_UNIT_TEST_SUITE )
    // The unit-test build has no GUI event loop to pump the WServer-posted work back through, so
    //  run the fit synchronously and invoke the callback inline.  This lets the LLM tool test
    //  suites exercise the async peak-fitting tools through the normal executeTool() path.
    try
    {
      const FitPeakStatus result = fit_user_peak( options, interspec );
      callback( result );
    }catch( const std::exception &e )
    {
      callback( string("Error in fit_user_peak: ") + e.what() );
    }
    return;
#endif

    Wt::WServer *server = Wt::WServer::instance();
    if( !server )
    {
      callback( string("No WServer instance available") );
      return;
    }

    if( options.explicitRoi )
    {
      // Validate and gather on this (GUI) thread, fit in the background, then commit on the GUI thread
      shared_ptr<const ExplicitRoiSetup> setup;
      try
      {
        setup = make_shared<const ExplicitRoiSetup>( prepare_explicit_roi( options, interspec ) );
      }catch( const exception &e )
      {
        callback( string("Error in fit_user_peak: ") + e.what() );
        return;
      }

      const string sessionid = wApp->sessionId();
      const weak_ptr<const SpecUtils::Measurement> weak_data = setup->data;

      // The background thread fits to its own copy of the spectrum
      auto bg_setup = make_shared<ExplicitRoiSetup>( *setup );
      bg_setup->data = make_shared<SpecUtils::Measurement>( *setup->data );

      server->ioService().boost::asio::io_service::post( [=](){
        const shared_ptr<ExplicitRoiFitResult> fit = make_shared<ExplicitRoiFitResult>();
        const shared_ptr<string> bg_error = make_shared<string>();
        try
        {
          *fit = fit_explicit_roi( *bg_setup );
        }catch( const exception &e )
        {
          *bg_error = string("Error in fit_user_peak: ") + e.what();
        }catch( ... )
        {
          *bg_error = "Error in fit_user_peak: unknown error fitting ROI";
        }

        Wt::WServer::instance()->post( sessionid, [=](){
          InterSpec *viewer = InterSpec::instance();
          if( !wApp || !viewer )
            return;

          if( !bg_error->empty() )
          {
            callback( *bg_error );
            return;
          }

          const shared_ptr<const SpecUtils::Measurement> current_data = weak_data.lock();
          if( !current_data || (viewer->displayedHistogram( setup->specType ) != current_data) )
          {
            callback( string("Spectrum changed while fitting peak - results discarded") );
            return;
          }

          try
          {
            const FitPeakStatus result = commit_explicit_roi( *setup, *fit, viewer );
            callback( result );
          }catch( const exception &e )
          {
            callback( string("Error in fit_user_peak: ") + e.what() );
          }

          wApp->triggerUpdate();
        } );//GUI-thread post
      } );//background post

      return;
    }//if( options.explicitRoi )

    // Stage A (GUI thread): capture all data needed for background computation
    try
    {
      shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
      if( !meas )
        throw runtime_error( "No measurement loaded for "
          + string(SpecUtils::descriptionText(options.specType)) + " spectrum" );

      const set<int> sample_nums = interspec->displayedSamples( options.specType );
      if( sample_nums.empty() )
        throw runtime_error( "No samples displayed for "
          + string(SpecUtils::descriptionText(options.specType)) + " spectrum" );

      shared_ptr<const SpecUtils::Measurement> data = interspec->displayedHistogram( options.specType );
      if( !data )
        throw runtime_error( "No spectrum displayed for "
          + string(SpecUtils::descriptionText(options.specType)) + " spectrum" );

      shared_ptr<const DetectorPeakResponse> det = meas->detector();

      shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks( sample_nums );
      assert( meas_peaks );
      if( !meas_peaks )
        throw runtime_error( "Unexpected error getting existing peak list" );

      shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_search_peaks = meas->automatedSearchPeaks( sample_nums );
      const shared_ptr<const PeakFitDetPrefs> fitPrefs = meas->peakFitDetPrefs();

      vector<shared_ptr<const PeakDef>> origPeaks;
      for( const shared_ptr<const PeakDef> &p : *meas_peaks )
      {
        if( p->gausPeak() )
          origPeaks.push_back( p );
        else if( (options.energy >= p->lowerX()) && (options.energy <= p->upperX()) )
          throw runtime_error( "Can not fit a peak within the ROI of a data-defined peak." );
      }

      // Make a copy of spectrum data for the background thread
      shared_ptr<SpecUtils::Measurement> data_copy = make_shared<SpecUtils::Measurement>( *data );
      weak_ptr<const SpecUtils::Measurement> weak_data = data;

      const string sessionid = wApp->sessionId();
      const double energy = options.energy;
      const string source = options.source.has_value() ? options.source.value() : ""s;
      const bool doNotAddToAnalysisPeaks = options.doNotAddToAnalysisPeaks;
      const SpecUtils::SpectrumType specType = options.specType;

      // Capture meas_peaks by value (shared_ptr copy) for Stage C
      const shared_ptr<deque<shared_ptr<const PeakDef>>> captured_meas_peaks = meas_peaks;

      // For describing what the data show if no peak is accepted: the detector class (which needs the
      //  session), and the automated-search peaks to fit a resolution to - if the search has not
      //  finished, its future, which the background thread can wait on without blocking this one.
      const PeakFitUtils::CoarseResolutionType det_type = resolve_det_type( interspec, meas, data, fitPrefs );
      std::shared_future<shared_ptr<const deque<shared_ptr<const PeakDef>>>> search_future;
      if( !auto_search_peaks )
      {
        try
        {
          search_future = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( meas, sample_nums, data, det, fitPrefs );
        }catch( std::exception & )
        {
        }
      }//if( !auto_search_peaks )

      using PeakPairResult = pair<vector<shared_ptr<const PeakDef>>, vector<shared_ptr<const PeakDef>>>;
      shared_ptr<PeakPairResult> found_peaks_result = make_shared<PeakPairResult>();
      shared_ptr<UnfitPeakDiagnostic> unfit_diag = make_shared<UnfitPeakDiagnostic>();
      shared_ptr<string> bg_error = make_shared<string>();

      // Stage B (background thread): run the expensive peak search
      server->ioService().boost::asio::io_service::post( std::bind( [=](){
        try
        {
          double pixelPerKev = -1.0; // Triggers "automed" peak fit
          string reject_reason;
          *found_peaks_result = searchForPeakFromUser( energy, pixelPerKev, data_copy, origPeaks, det, auto_search_peaks,
                                                       fitPrefs, nullptr, &reject_reason );
          // The peaks become the user's: fit them like any other (see the synchronous version above).
          if( !found_peaks_result->first.empty() )
            found_peaks_result->first = ExperimentalAutomatedPeakSearch::refit_sparse_rois( found_peaks_result->first,
                                          data_copy, PeakFitUtils::effective_det_type( fitPrefs, data_copy, nullptr ), {} );
          else
          {
            shared_ptr<const deque<shared_ptr<const PeakDef>>> search_peaks = auto_search_peaks;
            if( !search_peaks && search_future.valid() )
            {
              try
              {
                search_peaks = search_future.get();
              }catch( std::exception & )
              {
              }
            }//if( need to wait for the automated search )

            string fwhm_source;
            const std::function<float(double)> fwhm_fcn = fwhm_fcn_for_peak_fit( det_type, data_copy, det,
                                                                   origPeaks, search_peaks, fwhm_source );
            *unfit_diag = diagnose_unfit_peak( energy, reject_reason, fwhm_fcn( energy ), fwhm_source, data_copy,
                                               origPeaks, det, search_peaks, fitPrefs );
          }//if( found a peak ) / else
        }catch( const exception &e )
        {
          *bg_error = string("Peak search failed: ") + e.what();
        }catch( ... )
        {
          *bg_error = "Peak search failed with unknown error";
        }

        // Stage C (GUI thread): post-process results and call callback
        Wt::WServer::instance()->post( sessionid, std::bind( [=](){
          InterSpec *viewer = InterSpec::instance();
          if( !wApp || !viewer )
            return;

          // Check for background thread error
          if( !bg_error->empty() )
          {
            callback( *bg_error );
            return;
          }

          // Stale check
          shared_ptr<const SpecUtils::Measurement> current_data = weak_data.lock();
          if( !current_data || viewer->displayedHistogram( specType ) != current_data )
          {
            callback( string("Spectrum changed while fitting peak - results discarded") );
            return;
          }

          vector<shared_ptr<const PeakDef>> peaks_to_add_in = found_peaks_result->first;
          const vector<shared_ptr<const PeakDef>> &peaks_to_remove = found_peaks_result->second;

          if( peaks_to_add_in.empty() )
          {
            bool within_existing_roi = false;
            for( const shared_ptr<const PeakDef> &p : origPeaks )
            {
              if( (energy >= p->lowerX()) && (energy <= p->upperX()) )
              {
                within_existing_roi = true;
                break;
              }
            }

            if( doNotAddToAnalysisPeaks )
            {
              FitPeakStatus result;
              result.significant = false;
              result.fitPeak = unfit_diag->candidate;
              result.peaksInRoi = unfit_diag->candidateRoiPeaks;
              result.diagnostic = *unfit_diag;
              callback( result );
              return;
            }//if( only probing the energy )

            callback( unfit_peak_error_message( energy, within_existing_roi, *unfit_diag ) );
            return;
          }

          // Find the newly fit peak (the one without a source assigned), and assign source if specified
          shared_ptr<const PeakDef> fit_peak;
          try
          {
            fit_peak = assign_source_to_user_fit_peak( peaks_to_add_in, energy, source, origPeaks, viewer );
          }catch( std::exception &e )
          {
            callback( string(e.what()) );
            return;
          }

          // Update PeakModel
          const std::shared_ptr<PeakModel> pmodel = viewer->peakModel();
          assert( pmodel );
          if( !doNotAddToAnalysisPeaks && pmodel )
          {
            if( specType == SpecUtils::SpectrumType::Foreground )
            {
              pmodel->removePeaks( peaks_to_remove );
              pmodel->addPeaks( peaks_to_add_in );
            }else
            {
              shared_ptr<deque<shared_ptr<const PeakDef>>> new_deque
                = make_shared<deque<shared_ptr<const PeakDef>>>( *captured_meas_peaks );
              deque<shared_ptr<const PeakDef>> &peak_deque = *new_deque;

              for( const shared_ptr<const PeakDef> &old : peaks_to_remove )
              {
                auto pos = find( begin(peak_deque), end(peak_deque), old );
                assert( pos != end(peak_deque) );
                if( pos != end(peak_deque) )
                  peak_deque.erase( pos );
              }

              for( const shared_ptr<const PeakDef> &new_peak : peaks_to_add_in )
                peak_deque.push_back( new_peak );
              sort( begin(peak_deque), end(peak_deque), &PeakDef::lessThanByMeanShrdPtr );

              pmodel->setPeaks( peak_deque, specType );
            }
          }//if( !doNotAddToAnalysisPeaks )

          FitPeakStatus result;
          result.fitPeak = fit_peak;
          result.peaksInRoi = peaks_to_add_in;

          callback( result );

          wApp->triggerUpdate();
        } ) );//Stage C post
      } ) );//Stage B post

    }catch( const exception &e )
    {
      // Stage A failed
      callback( string("Error in fit_user_peak: ") + e.what() );
    }
  }//void fit_user_peak_async(...)


  GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    std::shared_ptr<SpecMeas> meas = interspec->measurment(options.specType);
    if (!meas) {
      throw std::runtime_error("No measurement loaded for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    const set<int> sample_nums = interspec->displayedSamples(options.specType);
    if (sample_nums.empty()) {
      throw std::runtime_error("No samples displayed for "
                               + std::string(SpecUtils::descriptionText(options.specType))
                               + " spectrum");
    }
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks(sample_nums);
    assert( meas_peaks );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );


    GetUserPeakStatus status;

    if( meas_peaks )
    {
      // Filter peaks by energy range if specified
      const bool has_lower = options.lowerEnergy.has_value();
      const bool has_upper = options.upperEnergy.has_value();

      // Validate energy range if both bounds are specified
      if( has_lower && has_upper && (*options.upperEnergy < *options.lowerEnergy) )
      {
        throw std::runtime_error( "upperEnergy (" + std::to_string(*options.upperEnergy)
                                  + " keV) must be greater than or equal to lowerEnergy ("
                                  + std::to_string(*options.lowerEnergy) + " keV)" );
      }

      if( !has_lower && !has_upper )
      {
        // No filtering, return all peaks
        status.peaks.insert( end(status.peaks), begin(*meas_peaks), end(*meas_peaks) );
      }else
      {
        const double lower_energy = has_lower ? *options.lowerEnergy : -std::numeric_limits<double>::infinity();
        const double upper_energy = has_upper ? *options.upperEnergy : std::numeric_limits<double>::infinity();

        for( const shared_ptr<const PeakDef> &peak : *meas_peaks )
        {
          if( peak && peak->mean() >= lower_energy && peak->mean() <= upper_energy )
            status.peaks.push_back( peak );
        }
      }
    }

    return status;
  }


  std::vector<std::string> get_identified_sources( const SpecUtils::SpectrumType specType, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error( "No InterSpec session available" );

    std::shared_ptr<SpecMeas> meas = interspec->measurment( specType );
    if( !meas )
    {
      throw std::runtime_error( "No measurement loaded for "
                                + std::string(SpecUtils::descriptionText(specType))
                                + " spectrum" );
    }

    const set<int> sample_nums = interspec->displayedSamples( specType );
    if( sample_nums.empty() )
    {
      throw std::runtime_error( "No samples displayed for "
                                + std::string(SpecUtils::descriptionText(specType))
                                + " spectrum" );
    }

    shared_ptr<const deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks( sample_nums );
    if( !meas_peaks )
      return std::vector<std::string>();

    // Use a set to track unique sources
    std::set<std::string> unique_sources;

    for( const shared_ptr<const PeakDef> &peak : *meas_peaks )
    {
      if( !peak )
        continue;

      // Check for nuclide source
      if( const SandiaDecay::Nuclide * const nuc = peak->parentNuclide() )
      {
        unique_sources.insert( nuc->symbol );
      }
      // Check for x-ray element source
      else if( const SandiaDecay::Element * const el = peak->xrayElement() )
      {
        unique_sources.insert( el->symbol );
      }
      // Check for reaction source
      else if( const ReactionGamma::Reaction * const rctn = peak->reaction() )
      {
        unique_sources.insert( rctn->name() );
      }
      // Check for user label (if it doesn't contain "unknown")
      else
      {
        const std::string &label = peak->userLabel();
        if( !label.empty() && !SpecUtils::icontains(label, "unknown") )
          unique_sources.insert( label );
      }
    }

    // Convert set to vector
    return std::vector<std::string>( unique_sources.begin(), unique_sources.end() );
  }


  vector<tuple<float,float,float>> cacl_estimated_gamma_importance( const vector<tuple<float,float>> &gamma_energies_and_yields )
  {
    // Photons that are never emitted are skipped: SandiaDecay always lists a 511 keV annihilation
    //  entry, with a zero rate for nuclides without positrons, and with no other photons (H3, P32,
    //  a stable nuclide) its importance would be 0/0 = NaN - which no "importance < cut" test removes.
    const auto emitted = []( const tuple<float,float> &gamma ) -> bool {
      return (std::get<0>(gamma) > 0.0f) && (std::get<1>(gamma) > 0.0f);
    };

    float sum_yields_sqrt_energy = 0.0f;
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      if( emitted(gamma) )
        sum_yields_sqrt_energy += std::get<1>(gamma) * sqrt(std::get<0>(gamma));
    }

    vector<tuple<float,float,float>> answer;
    if( !(sum_yields_sqrt_energy > 0.0f) )
      return answer;

    answer.reserve(gamma_energies_and_yields.size());
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      if( !emitted(gamma) )
        continue;
      const float importance = std::get<1>(gamma) * sqrt(std::get<0>(gamma)) / sum_yields_sqrt_energy;
      answer.push_back({std::get<0>(gamma), std::get<1>(gamma), importance});
    }
    return answer;
  }
  
  std::vector<float> get_characteristic_gammas( const std::string &nuclide )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    if( !db )
      throw std::runtime_error("No decay database available");

    const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
    if( !nuc )
    {
      // See if its a element for nuclide fluorescence
      string el_str = nuclide;
      const bool contained_xray = (SpecUtils::icontains(el_str, "x-ray")
                                   || SpecUtils::icontains(el_str, "x ray")
                                   || SpecUtils::icontains(el_str, "xray")
                                   || SpecUtils::icontains(el_str, "element")
                                   || SpecUtils::icontains(el_str, "fluorescence") );
      if( contained_xray )
      {
        SpecUtils::ireplace_all( el_str, "x-ray", "");
        SpecUtils::ireplace_all( el_str, "x ray", "");
        SpecUtils::ireplace_all( el_str, "xray", "");
        SpecUtils::ireplace_all( el_str, "element", "");
        SpecUtils::ireplace_all( el_str, "fluorescence", "");
      }
      SpecUtils::trim( el_str );
      
      const SandiaDecay::Element * const el = db->element( el_str );
      if( !el && !contained_xray )
      {
        if( SpecUtils::istarts_with( nuclide, "Neut" ) )
          return {478.0f, 847.0f, 2223.3f, 2235.3f};

        // See if its a reaction
        try
        {
          const ReactionGamma * const rctnDb = ReactionGammaServer::database();
          vector<ReactionGamma::ReactionPhotopeak> reactions;
          
          if( rctnDb )
            rctnDb->gammas( nuclide, reactions );

          if( reactions.empty() )
            throw std::runtime_error("No reactions found for '" + nuclide + "'");

          vector<tuple<float,float>> energies_and_yields;
          for( const auto &reaction : reactions )
            energies_and_yields.push_back({reaction.energy, reaction.abundance});

          vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(energies_and_yields);
          std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
              return std::get<2>(a) > std::get<2>(b);
          });

          if( importances.empty() )
            return {};

          const float max_importance = std::get<2>(importances[0]);
          vector<float> answer;
          for( const auto &importance : importances )
          {
            if( std::get<2>(importance) < max_importance * 0.25 )
              break;
            answer.push_back( std::get<0>(importance) );
          }
        

          if( answer.size() > 3 )
            answer.resize(3);

          return answer;
        }catch(...)
        {
          throw std::runtime_error("No nuclide, element, or reaction found for name '" + nuclide + "'");
        }//try / catch
      }//if( !el && !contained_xray )

      // For uranium and lead we'll just hard-code the result
      if( el->atomicNumber == 92 ) //Uranium
        return {98.4f};
      if( el->atomicNumber == 82 ) //Lead
        return {75.0f};

      // We'll return the highest "importance" x-ray for the element
      vector<tuple<float,float>> gamma_energies_and_yields;
      for( const auto &xray : el->xrays )
        gamma_energies_and_yields.push_back({xray.energy, xray.intensity});

      vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(gamma_energies_and_yields);
      std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
        return std::get<2>(a) > std::get<2>(b);
      });

      if( importances.empty() )
        return {};

      return {std::get<0>(importances[0])};
    }//if( !nuc )
    
    // TODO: should we use `photopeak_lis()`?
    
    const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &nuc_characteristics = IsotopeId::characteristicGammas();
    
    vector<float> answer;
    for( const auto &nuc_characteristic : nuc_characteristics )
    {
      if( get<0>(nuc_characteristic) == nuclide )
        answer.push_back( get<2>(nuc_characteristic) );
    }

    if( answer.empty() )
    {
      // We'll return up to a few of the highest "importance" gammas for the element.
      const double age = PeakDef::defaultDecayTime( nuc );
      SandiaDecay::NuclideMixture mixture;
      mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuc,1.0E4) );
      const double parent_activity = mixture.activity( age, nuc );
      if( !(parent_activity > 0.0) )
        return {};

      // The photons actually emitted.  No absolute yield floor: a line can be characteristic of a
      //  source however weak it is, if it is that source's strongest (Po210's 803 keV, at 1E-5 per
      //  decay); lines are ranked against the source's own strongest, below.
      vector<tuple<float,float>> gamma_energies_and_yields;
      const auto add_photons = [&]( const vector<SandiaDecay::EnergyRatePair> &photons ){
        for( const SandiaDecay::EnergyRatePair &photon : photons )
        {
          if( photon.numPerSecond > 0.0 )
            gamma_energies_and_yields.push_back( {photon.energy, photon.numPerSecond} );
        }
      };
      add_photons( mixture.gammas( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, false ) );
      add_photons( mixture.xrays( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy ) );

      // Annihilation radiation is not specific to any source (any positron emitter, or pair
      //  production, gives it), so it is a characteristic line only for a substantial positron
      //  emitter - the same per-decay floor `positron_emitters_for_annihilation(...)` uses.  (Cl36's
      //  ~3E-4 annihilation photons per decay are its only photons above 10 keV, so a floor relative
      //  to its own emission would not exclude it.)
      double positrons = 0.0;
      for( const SandiaDecay::EnergyRatePair &erp : mixture.betaPlusses( age ) )
        positrons += erp.numPerSecond;
      if( (2.0*positrons / parent_activity) >= sm_min_annihilation_photons_per_decay )
        gamma_energies_and_yields.push_back( {510.9989f, static_cast<float>(2.0*positrons)} );

      vector<tuple<float,float,float>> importances = cacl_estimated_gamma_importance(gamma_energies_and_yields);
      std::sort( begin(importances), end(importances), [](const tuple<float,float,float> &a, const tuple<float,float,float> &b) {
        return std::get<2>(a) > std::get<2>(b);
      });

      if( importances.empty() )
        return {};

      const float max_importance = std::get<2>(importances[0]);
      for( const auto &importance : importances )
      {
        if( std::get<2>(importance) < max_importance * 0.25 )
          break;
        answer.push_back( std::get<0>(importance) );
      }

      if( answer.size() > 3 )
        answer.resize(3);
    }//if( answer.empty() )

    return answer;
  }//std::vector<float> get_characteristic_gammas( const std::string &nuclide )


  std::vector<AnnihilationSource> positron_emitters_for_annihilation( InterSpec *interspec,
                                                                      const double min_photons_per_decay )
  {
    const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
    if( !db )
      throw std::runtime_error( "No decay database available" );

    // The nuclides an analyst would usually consider, plus the common positron emitters, as only
    //  Na22 of these is in the search categories.
    set<const SandiaDecay::Nuclide *> nuclides;

    std::vector<IsotopeSearchByEnergy::NucSearchCategory> nuc_categories;
    if( interspec && interspec->nuclideSearch() )
      nuc_categories = interspec->nuclideSearch()->search_categories();
    else
      IsotopeSearchByEnergy::init_category_info( nuc_categories );

    for( const IsotopeSearchByEnergy::NucSearchCategory &category : nuc_categories )
      nuclides.insert( begin(category.m_specific_nuclides), end(category.m_specific_nuclides) );

    for( const auto &nuc_characteristic : IsotopeId::characteristicGammas() )
    {
      if( get<1>(nuc_characteristic) )
        nuclides.insert( get<1>(nuc_characteristic) );
    }

    for( const char * const name : { "Na22", "F18", "Ga68", "Ge68", "Cu64", "Rb82", "N13", "C11",
                                     "O15", "Zr89", "I124", "Co58", "Co56" } )
    {
      if( const SandiaDecay::Nuclide * const nuc = db->nuclide( name ) )
        nuclides.insert( nuc );
    }

    std::vector<AnnihilationSource> answer;
    for( const SandiaDecay::Nuclide * const nuc : nuclides )
    {
      if( !nuc || nuc->isStable() )
        continue;

      // Positrons of the nuclide and its progeny, at its usual age, per decay of the nuclide.
      const double activity = 1.0*SandiaDecay::Bq;
      SandiaDecay::NuclideMixture mixture;
      mixture.addAgedNuclideByActivity( nuc, activity, PeakDef::defaultDecayTime( nuc ) );

      double positrons = 0.0;
      for( const SandiaDecay::EnergyRatePair &erp : mixture.betaPlusses( 0.0 ) )
        positrons += erp.numPerSecond;

      const double photons_per_decay = 2.0 * positrons / activity;
      if( photons_per_decay >= min_photons_per_decay )
        answer.push_back( AnnihilationSource{ nuc, photons_per_decay } );
    }//for( loop over nuclides )

    std::sort( begin(answer), end(answer), []( const AnnihilationSource &a, const AnnihilationSource &b ){
      return (a.photons511PerDecay > b.photons511PerDecay)
              || ((a.photons511PerDecay == b.photons511PerDecay) && (a.nuclide->symbol < b.nuclide->symbol));
    } );

    return answer;
  }//positron_emitters_for_annihilation(...)
  
  
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_nuclides_with_characteristics_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec  )
  {
    std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>> answer;
    
    if( lower_energy > upper_energy )
      std::swap( lower_energy, upper_energy );
    
    // First, we'll identify all the common nuclides, elements, and reactions that are 
    //  typically encountered.
    set<const SandiaDecay::Nuclide *> nuclides;
    set<const SandiaDecay::Element *> elements;
    set<const ReactionGamma::Reaction *> reactions;
    
    std::vector<IsotopeSearchByEnergy::NucSearchCategory> nuc_categories;
    if( interspec && interspec->nuclideSearch() )
      nuc_categories = interspec->nuclideSearch()->search_categories();
    else
      IsotopeSearchByEnergy::init_category_info( nuc_categories );

    for( const IsotopeSearchByEnergy::NucSearchCategory &category : nuc_categories )
    {
      nuclides.insert( begin(category.m_specific_nuclides), end(category.m_specific_nuclides) );
      elements.insert( begin(category.m_specific_elements), end(category.m_specific_elements) );
      reactions.insert( begin(category.m_specific_reactions), end(category.m_specific_reactions) );
    }

    const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &nuc_characteristics = IsotopeId::characteristicGammas();
    for( const auto &nuc_characteristic : nuc_characteristics )
    {
      if( get<1>(nuc_characteristic) )
        nuclides.insert( get<1>(nuc_characteristic) );
    }

    // Now add all the nuclides that have characteristics in the energy range
    for( const auto &nuclide : nuclides )
    {
      const vector<float> gammas = get_characteristic_gammas( nuclide->symbol );
      for( const float gamma : gammas )
      {
        if( gamma >= lower_energy && gamma <= upper_energy )
          answer.push_back( nuclide );
      }
    }

    for( const auto &element : elements )
    {
      const vector<float> xrays = get_characteristic_gammas( element->symbol );
      for( const float xray : xrays )
      {
        if( xray >= lower_energy && xray <= upper_energy )
          answer.push_back( element );
      }
    }

    for( const auto &reaction : reactions )
    {
      const vector<float> gammas = get_characteristic_gammas( reaction->name() );
      for( const float gamma : gammas )
      {
        if( gamma >= lower_energy && gamma <= upper_energy )
          answer.push_back( reaction );
      }
    }

    return answer;
  }//get_nuclides_with_characteristics_in_energy_range(...)
  
  
  FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec )
  {
    // NOTE: this function will eventually be replaced or powered by FitPeaksForNuclideDev::fit_peaks_for_nuclides(...)
    // TODO: This function tries a few different RelActAuto configuration, and chooses the one with the best Chi2. It should actually be Chi2/DOF, and then maybe also take into account what peaks are chosen
    // TODO: The order of Ln(x) equation used should be based on how many peaks, and the energy range - should do something more inteligently
    // TODO: Should split range up - e.g. above and below 120 keV, if nuclide has peaks on both sides
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    FitPeaksForNuclideStatus result;
    
    try
    {
      const SpecUtils::SpectrumType specType = options.specType;

      // Get the target spectrum (the one we are fitting peaks to)
      std::shared_ptr<SpecMeas> meas = interspec->measurment( specType );
      if( !meas )
        throw std::runtime_error( std::string("No ") + SpecUtils::descriptionText(specType) + " measurement available" );

      std::shared_ptr<const SpecUtils::Measurement> target_spectrum = interspec->displayedHistogram( specType );
      if( !target_spectrum )
        throw std::runtime_error( std::string("No ") + SpecUtils::descriptionText(specType) + " spectrum available" );

      // Background subtraction only makes sense when fitting to the foreground
      std::shared_ptr<const SpecUtils::Measurement> background;
      if( specType == SpecUtils::SpectrumType::Foreground )
        background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
      
      // Get detector response function; fall back to foreground SpecMeas's DRF if needed
      std::shared_ptr<const DetectorPeakResponse> drf = meas->detector();
      std::shared_ptr<const PeakFitDetPrefs> peak_fit_prefs = meas ? meas->peakFitDetPrefs() : nullptr;
      if( !drf && (specType != SpecUtils::SpectrumType::Foreground) )
      {
        std::shared_ptr<SpecMeas> fg_meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
        if( fg_meas )
        {
          drf = fg_meas->detector();
          if( !peak_fit_prefs )
            peak_fit_prefs = fg_meas->peakFitDetPrefs();
        }
      }

      DetectedPeaksOptions det_peaks_options;
      det_peaks_options.specType = specType;
      det_peaks_options.nonBackgroundPeaksOnly = false; // TODO: as we improve `FitPeaksForNuclides::fit_peaks_for_nuclides(...)` this will need to be change
      const DetectedPeakStatus detected_peaks = AnalystChecks::detected_peaks( det_peaks_options, interspec );
      const std::vector<std::shared_ptr<const PeakDef>> &auto_search_peaks = detected_peaks.peaks;

      vector<RelActCalcAuto::SrcVariant> sources;
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      if( !db )
        throw runtime_error( "No decay database available" );

      for( const string &src : options.sources )
      {
        if( SpecUtils::iequals_ascii( src, "null" ) )
          continue;

        RelActCalcAuto::SrcVariant source = RelActCalcAuto::source_from_string(src);

        if( RelActCalcAuto::is_null(source) )
          throw runtime_error( "Unrecognized source '" + src + "'" );

        // Skip duplicate sources (e.g., LLM may pass the same nuclide twice)
        const bool is_duplicate = std::any_of( sources.begin(), sources.end(),
          [&source]( const RelActCalcAuto::SrcVariant &s ){ return s == source; } );
        if( is_duplicate )
          continue;

        sources.push_back( source );
      }//for( const string &src : options.sources )

      const bool fit_norm_peaks = options.fitSrcPeaksOptions.test( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks )
                                  || options.fitSrcPeaksOptions.test( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse );
      if( sources.empty() && !fit_norm_peaks )
        throw runtime_error( "No sources specified" );

      GetUserPeakOptions user_peak_options;
      user_peak_options.specType = specType;
      //user_peak_options.lowerEnergy = ...;
      //user_peak_options.upperEnergy = ...;
      
      const GetUserPeakStatus user_peaks_status = AnalystChecks::get_user_peaks( user_peak_options, interspec );
      const std::vector<std::shared_ptr<const PeakDef>> &user_peaks = user_peaks_status.peaks;
      
      const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> fit_options = options.fitSrcPeaksOptions;
      
      const PeakFitUtils::CoarseResolutionType det_type = peak_fit_prefs
        ? peak_fit_prefs->m_det_type
        : PeakFitUtils::coarse_det_type( target_spectrum, meas );
      FitPeaksForNuclides::PeakFitForNuclideConfig fit_config = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( det_type );

      // The skew type comes from the peak-fit prefs (e.g., a GADRAS detector's peak shape), as in the GUI
      if( peak_fit_prefs )
        fit_config.skew_type = peak_fit_prefs->m_peak_skew_type;

      // TODO: we will need to update `config` from default in the future

      const FitPeaksForNuclides::PeakFitResult fit_results = FitPeaksForNuclides::fit_peaks_for_nuclides(
        auto_search_peaks, target_spectrum, sources, user_peaks, background, drf, fit_options, fit_config, peak_fit_prefs
      );

      switch( fit_results.status )
      {
        case RelActCalcAuto::RelActAutoSolution::Status::Success:
        case RelActCalcAuto::RelActAutoSolution::Status::UsableWithWarnings:
          break;

        case RelActCalcAuto::RelActAutoSolution::Status::NotInitiated:
          throw runtime_error( "Failed to initialize setting up the peak fit"
                              + (fit_results.error_message.empty() ? string("") : (string("; error: ") + fit_results.error_message) ) );

        case RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem:
          throw runtime_error( "Failed to setup up the peak fit"
                              + (fit_results.error_message.empty() ? string("") : (string("; error: ") + fit_results.error_message) ) );

        case RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem:
          throw runtime_error( "Failed to solve the peak fit"
                              + (fit_results.error_message.empty() ? string("") : (string("; error: ") + fit_results.error_message) ) );

        case RelActCalcAuto::RelActAutoSolution::Status::UserCanceled:
          throw runtime_error( "Peak fit failed"
                              + (fit_results.error_message.empty() ? string("") : (string("; error: ") + fit_results.error_message) ) );
      }//switch( fit_results.status )

      result.warnings = fit_results.warnings;
      map<string,Wt::WColor> pending_colors;  // So each new source gets its own color
      for( const PeakDef &p : fit_results.observable_peaks )
      {
        auto sp = make_shared<PeakDef>(p);

        // Keep colors the fitter gave (e.g., gray for NORM lines it added on its own)
        if( sp->lineColor().isDefault() )
          sp->setLineColor( color_for_new_peak( *sp, interspec, pending_colors ) );

        result.fitPeaks.push_back( sp );
      }//for( const PeakDef &p : fit_results.observable_peaks )

      if( result.fitPeaks.empty() )
        result.noObservablePeaksFit = true;

      result.peaksToRemove = fit_results.original_peaks_to_remove;

      // If doNotAddPeaksToUserSession is false, update the peak model using setPeaks
      //  to atomically replace the affected peaks while preserving unrelated ones.
      if( !options.doNotAddPeaksToUserSession && !result.fitPeaks.empty() )
      {
        const std::shared_ptr<PeakModel> pmodel = interspec->peakModel();
        if( pmodel )
        {
          std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> existing_peaks
            = pmodel->peaks( specType );

          std::deque<std::shared_ptr<const PeakDef>> combined_peaks;

          if( existing_peaks )
          {
            // Keep existing peaks that are NOT in peaksToRemove
            for( const std::shared_ptr<const PeakDef> &p : *existing_peaks )
            {
              const bool should_remove = std::any_of( result.peaksToRemove.begin(),
                result.peaksToRemove.end(),
                [&p]( const std::shared_ptr<const PeakDef> &rm ) -> bool { return rm == p; } );

              if( !should_remove )
                combined_peaks.push_back( p );
            }
          }//if( existing_peaks )

          // Add the newly fit peaks
          for( const std::shared_ptr<const PeakDef> &p : result.fitPeaks )
            combined_peaks.push_back( p );

#if( PERFORM_DEVELOPER_CHECKS )
          // Verify every removed sibling has a corresponding replacement peak in combined_peaks
          // whose ROI covers the sibling's energy.
          for( const std::shared_ptr<const PeakDef> &rm : result.peaksToRemove )
          {
            const double rm_energy = rm->mean();
            bool has_replacement = false;
            for( const std::shared_ptr<const PeakDef> &cp : combined_peaks )
            {
              if( (rm_energy >= cp->lowerX()) && (rm_energy <= cp->upperX()) )
              {
                has_replacement = true;
                break;
              }
            }

            if( !has_replacement )
            {
              std::cerr << "WARNING: removed peak at " << rm_energy
                   << " keV (source: " << rm->sourceName()
                   << ") has no replacement peak in combined_peaks." << std::endl;
              assert( has_replacement );
            }
          }
#endif

          pmodel->setPeaks( combined_peaks, specType );
          result.peaksWereRemovedFromSession = true;
        }
      }//if( !options.doNotAddPeaksToUserSession )
    }catch( const std::exception &e )
    {
      throw std::runtime_error("Error in fit_peaks_for_nuclides: " + string(e.what()));
    }//try / catch

    return result;
  }//FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec )


  void fit_peaks_for_nuclides_async( const FitPeaksForNuclideOptions &options, InterSpec *interspec,
                                      FitPeaksCallback callback )
  {
    assert( callback );
    assert( interspec );
    assert( wApp );

    if( !interspec )
    {
      callback( string("No InterSpec session available") );
      return;
    }

#if( BUILD_AS_UNIT_TEST_SUITE )
    // See fit_user_peak_async: run synchronously in the unit-test build (no GUI event loop) so the
    //  async tool can be driven through the normal executeTool() path.
    try
    {
      const FitPeaksForNuclideStatus result = fit_peaks_for_nuclides( options, interspec );
      callback( result );
    }catch( const std::exception &e )
    {
      callback( string("Error in fit_peaks_for_nuclides: ") + e.what() );
    }
    return;
#endif

    Wt::WServer *server = Wt::WServer::instance();
    if( !server )
    {
      callback( string("No WServer instance available") );
      return;
    }

    // Stage A (GUI thread): capture all data needed for background computation
    try
    {
      const SpecUtils::SpectrumType specType = options.specType;

      shared_ptr<SpecMeas> meas = interspec->measurment( specType );
      if( !meas )
        throw runtime_error( string("No ") + SpecUtils::descriptionText(specType) + " measurement available" );

      shared_ptr<const SpecUtils::Measurement> target_spectrum = interspec->displayedHistogram( specType );
      if( !target_spectrum )
        throw runtime_error( string("No ") + SpecUtils::descriptionText(specType) + " spectrum available" );

      shared_ptr<const SpecUtils::Measurement> background;
      if( specType == SpecUtils::SpectrumType::Foreground )
        background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );

      shared_ptr<const DetectorPeakResponse> drf = meas->detector();
      shared_ptr<const PeakFitDetPrefs> peak_fit_prefs = meas->peakFitDetPrefs();
      if( !drf && (specType != SpecUtils::SpectrumType::Foreground) )
      {
        shared_ptr<SpecMeas> fg_meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
        if( fg_meas )
        {
          drf = fg_meas->detector();
          if( !peak_fit_prefs )
            peak_fit_prefs = fg_meas->peakFitDetPrefs();
        }
      }

      DetectedPeaksOptions det_peaks_options;
      det_peaks_options.specType = specType;
      det_peaks_options.nonBackgroundPeaksOnly = false;
      const DetectedPeakStatus detected_peaks = AnalystChecks::detected_peaks( det_peaks_options, interspec );
      const vector<shared_ptr<const PeakDef>> auto_search_peaks = detected_peaks.peaks;

      vector<RelActCalcAuto::SrcVariant> sources;
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      if( !db )
        throw runtime_error( "No decay database available" );

      for( const string &src : options.sources )
      {
        if( SpecUtils::iequals_ascii( src, "null" ) )
          continue;

        RelActCalcAuto::SrcVariant source = RelActCalcAuto::source_from_string( src );
        if( RelActCalcAuto::is_null( source ) )
          throw runtime_error( "Unrecognized source '" + src + "'" );

        const bool is_duplicate = any_of( sources.begin(), sources.end(),
          [&source]( const RelActCalcAuto::SrcVariant &s ){ return s == source; } );
        if( is_duplicate )
          continue;

        sources.push_back( source );
      }

      const bool fit_norm_peaks = options.fitSrcPeaksOptions.test( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks )
                                  || options.fitSrcPeaksOptions.test( FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse );
      if( sources.empty() && !fit_norm_peaks )
        throw runtime_error( "No sources specified" );

      GetUserPeakOptions user_peak_options;
      user_peak_options.specType = specType;
      const GetUserPeakStatus user_peaks_status = AnalystChecks::get_user_peaks( user_peak_options, interspec );
      const vector<shared_ptr<const PeakDef>> user_peaks = user_peaks_status.peaks;

      const Wt::WFlags<FitPeaksForNuclides::FitSrcPeaksOptions> fit_options = options.fitSrcPeaksOptions;
      const PeakFitUtils::CoarseResolutionType det_type = peak_fit_prefs
        ? peak_fit_prefs->m_det_type
        : PeakFitUtils::coarse_det_type( target_spectrum, meas );
      FitPeaksForNuclides::PeakFitForNuclideConfig fit_config = FitPeaksForNuclides::PeakFitForNuclideConfig::default_config( det_type );

      // The skew type comes from the peak-fit prefs (e.g., a GADRAS detector's peak shape), as in the GUI
      if( peak_fit_prefs )
        fit_config.skew_type = peak_fit_prefs->m_peak_skew_type;

      // Make copies of spectra for the background thread (same pattern as FitPeaksForNuclidesGui.cpp)
      shared_ptr<SpecUtils::Measurement> target_copy = make_shared<SpecUtils::Measurement>( *target_spectrum );
      shared_ptr<SpecUtils::Measurement> bg_copy = background
                                    ? make_shared<SpecUtils::Measurement>( *background ) : nullptr;

      weak_ptr<const SpecUtils::Measurement> weak_target = target_spectrum;
      const string sessionid = wApp->sessionId();
      const bool doNotAddPeaksToUserSession = options.doNotAddPeaksToUserSession;

      // Stage B (background thread): run the expensive peak fitting
      shared_ptr<FitPeaksForNuclides::PeakFitResult> fit_result
        = make_shared<FitPeaksForNuclides::PeakFitResult>();

      server->ioService().boost::asio::io_service::post( std::bind( [=](){
        fit_result->status = RelActCalcAuto::RelActAutoSolution::Status::Success;
        fit_result->error_message.clear();

        try
        {
          *fit_result = FitPeaksForNuclides::fit_peaks_for_nuclides(
            auto_search_peaks, target_copy, sources, user_peaks, bg_copy, drf, fit_options, fit_config, peak_fit_prefs
          );
        }catch( const exception &e )
        {
          fit_result->status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
          fit_result->error_message = string("Fit failed: ") + e.what();
        }catch( ... )
        {
          fit_result->status = RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem;
          fit_result->error_message = "Fit failed with unknown error";
        }

        // Stage C (GUI thread): post-process results and call callback
        Wt::WServer::instance()->post( sessionid, std::bind( [=](){
          InterSpec *viewer = InterSpec::instance();
          if( !wApp || !viewer )
          {
            // Session is gone; silently discard results
            return;
          }

          // Stale check: if the spectrum has changed, discard results
          shared_ptr<const SpecUtils::Measurement> current_target = weak_target.lock();
          if( !current_target || viewer->displayedHistogram( specType ) != current_target )
          {
            callback( string("Spectrum changed while fitting peaks - results discarded") );
            return;
          }

          // Check fit status
          switch( fit_result->status )
          {
            case RelActCalcAuto::RelActAutoSolution::Status::Success:
            case RelActCalcAuto::RelActAutoSolution::Status::UsableWithWarnings:
              break;
            case RelActCalcAuto::RelActAutoSolution::Status::NotInitiated:
              callback( string("Failed to initialize peak fit")
                + (fit_result->error_message.empty() ? "" : ("; error: " + fit_result->error_message)) );
              return;
            case RelActCalcAuto::RelActAutoSolution::Status::FailedToSetupProblem:
              callback( string("Failed to setup peak fit")
                + (fit_result->error_message.empty() ? "" : ("; error: " + fit_result->error_message)) );
              return;
            case RelActCalcAuto::RelActAutoSolution::Status::FailToSolveProblem:
              callback( string("Failed to solve peak fit")
                + (fit_result->error_message.empty() ? "" : ("; error: " + fit_result->error_message)) );
              return;
            case RelActCalcAuto::RelActAutoSolution::Status::UserCanceled:
              callback( string("Peak fit failed")
                + (fit_result->error_message.empty() ? "" : ("; error: " + fit_result->error_message)) );
              return;
          }//switch( fit_result->status )

          // Build FitPeaksForNuclideStatus from results
          FitPeaksForNuclideStatus result;
          result.warnings = fit_result->warnings;

          map<string,Wt::WColor> pending_colors;  // So each new source gets its own color

          for( const PeakDef &p : fit_result->observable_peaks )
          {
            shared_ptr<PeakDef> sp = make_shared<PeakDef>( p );

            // Keep colors the fitter gave (e.g., gray for NORM lines it added on its own)
            if( sp->lineColor().isDefault() )
              sp->setLineColor( color_for_new_peak( *sp, viewer, pending_colors ) );

            result.fitPeaks.push_back( sp );
          }

          if( result.fitPeaks.empty() )
          {
            result.noObservablePeaksFit = true;
            const string lost_rois_warning = "Lost all ROIs while iterationing to refine solution - stopped early.";
            result.warnings.erase(
              remove( result.warnings.begin(), result.warnings.end(), lost_rois_warning ),
              result.warnings.end()
            );
          }

          result.peaksToRemove = fit_result->original_peaks_to_remove;

          if( !doNotAddPeaksToUserSession && !result.fitPeaks.empty() )
          {
            const std::shared_ptr<PeakModel> pmodel = viewer->peakModel();
            if( pmodel )
            {
              shared_ptr<const deque<shared_ptr<const PeakDef>>> existing_peaks
                = pmodel->peaks( specType );

              deque<shared_ptr<const PeakDef>> combined_peaks;

              if( existing_peaks )
              {
                for( const shared_ptr<const PeakDef> &p : *existing_peaks )
                {
                  const bool should_remove = any_of( result.peaksToRemove.begin(),
                    result.peaksToRemove.end(),
                    [&p]( const shared_ptr<const PeakDef> &rm ) -> bool { return rm == p; } );

                  if( !should_remove )
                    combined_peaks.push_back( p );
                }
              }

              for( const shared_ptr<const PeakDef> &p : result.fitPeaks )
                combined_peaks.push_back( p );

              pmodel->setPeaks( combined_peaks, specType );
              result.peaksWereRemovedFromSession = true;
            }
          }//if( !doNotAddPeaksToUserSession && !result.fitPeaks.empty() )

          callback( result );

          wApp->triggerUpdate();
        } ) );//Stage C post
      } ) );//Stage B post

    }catch( const exception &e )
    {
      // Stage A failed; call callback with error synchronously
      callback( string("Error in fit_peaks_for_nuclides: ") + e.what() );
    }
  }//void fit_peaks_for_nuclides_async(...)


  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_characteristics_near_energy( const double energy, InterSpec *interspec )
  {
    const float fwhm = get_expected_fwhm( energy, interspec );

    return get_nuclides_with_characteristics_in_energy_range( energy - fwhm, energy + fwhm, interspec );
  }//get_characteristics_near_energy(...)
  
  
  std::function<float(double)> make_expected_fwhm_fcn( const shared_ptr<const DetectorPeakResponse> &drf,
                                                       const shared_ptr<const deque<shared_ptr<const PeakDef>>> &peaks,
                                                       const shared_ptr<const SpecUtils::Measurement> &spectrum,
                                                       const PeakFitUtils::CoarseResolutionType det_type,
                                                       std::string *source )
  {
    using DrfResForm = DetectorPeakResponse::ResolutionFnctForm;

    const bool isHPGe = (det_type == PeakFitUtils::CoarseResolutionType::High);
    const bool drf_has_res = (drf && drf->hasResolutionInfo());

    vector<shared_ptr<const PeakDef>> gauss_peaks;
    if( peaks )
    {
      for( const shared_ptr<const PeakDef> &p : *peaks )
      {
        if( p && p->gausPeak() && (p->amplitude() > 0.0) && (p->fwhm() > 0.0) )
          gauss_peaks.push_back( p );
      }
    }//if( peaks )

    // Peaks that measure the resolution well enough to overrule a DRF: significant, and not the
    //  Doppler-broadened 511 keV annihilation peak.  Their energy span must be decent too.
    size_t num_measuring = 0;
    double min_measuring_energy = std::numeric_limits<double>::max(), max_measuring_energy = 0.0;
    for( const shared_ptr<const PeakDef> &p : gauss_peaks )
    {
      const double nsigma = (p->amplitudeUncert() > 0.0) ? (p->amplitude() / p->amplitudeUncert()) : 0.0;
      if( (nsigma >= 6.0) && (fabs( p->mean() - 510.999 ) > 3.0*p->fwhm()) )
      {
        num_measuring += 1;
        min_measuring_energy = std::min( min_measuring_energy, p->mean() );
        max_measuring_energy = std::max( max_measuring_energy, p->mean() );
      }
    }//for( const shared_ptr<const PeakDef> &p : gauss_peaks )
    const bool peaks_measure_resolution = (num_measuring >= 3) && (max_measuring_energy >= 2.0*min_measuring_energy);

    // A resolution curve fit to the peaks, anchored to the shape of the DRF's curve (or, without
    //  one, the detector class's), so a few peaks give a sensible curve everywhere; the scale of
    //  the DRF that best matches the peaks (`drf_scale`) says if the DRF is for this detector.
    const DrfResForm fit_form = DrfResForm::kSqrtPolynomial;
    vector<float> fit_coefs;
    double fit_lower = 0.0, fit_upper = 0.0, drf_scale = 0.0;
    if( !gauss_peaks.empty() && spectrum )
    {
      try
      {
        std::function<double(double)> shape_prior;
        if( drf_has_res )
          shape_prior = [drf]( const double e ) -> double { return drf->peakResolutionFWHM( static_cast<float>(e) ); };

        // The curve is fit above this floor (it must be positive), and held flat below it.
        const double analysis_floor_kev = isHPGe ? 10.0 : 20.0;

        vector<float> fit_uncerts;
        string note;
        FitPeaksForNuclides::detail::fit_fwhm_function_robust( gauss_peaks, spectrum, det_type, analysis_floor_kev, shape_prior,
                                     fit_form, fit_coefs, fit_uncerts, fit_lower, fit_upper, note, nullptr,
                                     drf_has_res ? &drf_scale : nullptr );
      }catch( std::exception & )
      {
        fit_coefs.clear();
      }
    }//if( we have peaks to fit to )

    // A DRF whose FWHM is off by more than a third from what the peaks show is not for this
    //  detector (e.g., a generic DRF), so the peaks win; otherwise the DRF is used unchanged.
    const double max_drf_scale_off = 1.33;
    const bool drf_overruled = drf_has_res && !fit_coefs.empty() && peaks_measure_resolution
                               && (drf_scale > 0.0) && (fabs( std::log( drf_scale ) ) > std::log( max_drf_scale_off ));

    shared_ptr<const DetectorPeakResponse> res_drf = (drf_has_res && !drf_overruled) ? drf : nullptr;
    if( res_drf || fit_coefs.empty() )
      fit_coefs.clear();

    // If the robust fit could not be made, a plain resolution fit to the peaks, as before.
    shared_ptr<DetectorPeakResponse> fit_drf;
    if( !res_drf && fit_coefs.empty() && !gauss_peaks.empty() && spectrum )
    {
      try
      {
        auto peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(gauss_peaks), end(gauss_peaks) );
        fit_drf = make_shared<DetectorPeakResponse>();
        fit_drf->setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                               0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic);
        fit_drf->fitResolution( peak_deque, spectrum, DrfResForm::kSqrtPolynomial );
      }catch( std::exception & )
      {
        fit_drf.reset();
      }
    }//if( need a plain resolution fit )

    // GADRAS-parameter wag by detection type, as the final fallback.
    const vector<float> wag_pars = isHPGe ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };

    if( source )
    {
      char buffer[256] = { '\0' };
      if( res_drf )
      {
        *source = "detector response function";
      }else if( !fit_coefs.empty() || fit_drf )
      {
        snprintf( buffer, sizeof(buffer), "fit to the widths of %i peaks in the spectrum", static_cast<int>(gauss_peaks.size()) );
        *source = buffer;
        if( drf_overruled )
        {
          snprintf( buffer, sizeof(buffer), " (the detector response function's FWHM is %.2g times what the peaks show,"
                    " so it is likely not for this detector)", 1.0/drf_scale );
          *source += buffer;
        }
      }else if( !gauss_peaks.empty() )
      {
        *source = "width of the nearest peak, if within 20% of the energy, else a generic "
                  + string(isHPGe ? "HPGe" : "NaI") + " estimate";
      }else
      {
        *source = "generic " + string(isHPGe ? "HPGe" : "NaI") + " resolution estimate";
      }
    }//if( source )

    return [res_drf, fit_coefs, fit_form, fit_lower, fit_upper, fit_drf, gauss_peaks, wag_pars]( const double energy ) -> float
    {
      float fwhm = -1.0f;

      if( res_drf )
        fwhm = res_drf->peakResolutionFWHM( static_cast<float>(energy) );

      if( (fwhm <= 0.0f) && !fit_coefs.empty() )
        fwhm = static_cast<float>( FitPeaksForNuclides::detail::model_fwhm_at( energy, fit_form, fit_coefs,
                                                                               fit_lower, fit_upper ) );

      if( (fwhm <= 0.0f) && fit_drf )
        fwhm = fit_drf->peakResolutionFWHM( static_cast<float>(energy) );

      // Find nearest peak, and if within 20% of the energy, use it.
      if( (fwhm <= 0.0f) && !gauss_peaks.empty() )
      {
        double smallest_energy_diff = std::numeric_limits<double>::max();
        double nearest_fwhm = 0.0;
        for( const auto &peak : gauss_peaks )
        {
          const double energy_diff = std::abs(peak->mean() - energy);
          if( energy_diff < smallest_energy_diff )
          {
            smallest_energy_diff = energy_diff;
            nearest_fwhm = peak->fwhm();
          }
        }

        if( smallest_energy_diff < 0.2*energy )
          fwhm = static_cast<float>( nearest_fwhm );
      }//if( (fwhm <= 0.0) && !gauss_peaks.empty() )

      if( !(fwhm > 0.0f) || IsNan(fwhm) || IsInf(fwhm) )
        fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, DrfResForm::kGadrasResolutionFcn, wag_pars );

      return fwhm;
    };
  }//make_expected_fwhm_fcn(...)


  std::function<float(double)> make_expected_fwhm_fcn( InterSpec *interspec, std::string *source )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");

    shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !meas || !spectrum )
      throw std::runtime_error("No foreground loaded");

    // The detected peaks, to fit a resolution to, or to check the DRF's resolution against (which
    //  may run the automated peak search, once for the spectrum).
    shared_ptr<const deque<shared_ptr<const PeakDef>>> peak_deque;
    try
    {
      DetectedPeaksOptions opts;
      opts.specType = SpecUtils::SpectrumType::Foreground;
      opts.nonBackgroundPeaksOnly = false;
      const DetectedPeakStatus peaks = detected_peaks( opts, interspec );
      peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks.peaks), end(peaks.peaks) );
    }catch( std::exception & )
    {
    }

    return make_expected_fwhm_fcn( meas->detector(), peak_deque, spectrum,
                                   resolve_det_type( interspec, meas, spectrum ), source );
  }//make_expected_fwhm_fcn( InterSpec * )


  float get_expected_fwhm( const double energy, InterSpec *interspec, std::string *source )
  {
    if( !(energy > 0.0) || IsInf(energy) )
      throw std::runtime_error( "Energy must be positive to get an expected FWHM." );

    const std::function<float(double)> fwhm_fcn = make_expected_fwhm_fcn( interspec, source );

    const float fwhm = fwhm_fcn( energy );

    if( fwhm <= 0.0 )
      throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy) );

    return fwhm;
  }//get_expected_fwhm(...)
  
  
  /** Helper to combine channel counts by summing every `factor` channels.
   * Input `counts` must have size divisible by `factor` (or the last bin covers fewer channels).
   */
  static vector<float> combine_channels( const vector<float> &counts, const size_t factor )
  {
    assert( factor >= 1 );
    if( factor <= 1 )
      return counts;

    const size_t ncombined = (counts.size() + factor - 1) / factor;
    vector<float> result( ncombined, 0.0f );
    for( size_t i = 0; i < counts.size(); ++i )
      result[i / factor] += counts[i];

    return result;
  }//combine_channels(...)


  /** Helper to extract every `factor`-th lower-edge energy from the sub-range energies.
   * `sub_energies` has (raw_num_channels + 1) entries (lower edges + final upper edge).
   * Returns just the lower-edge of each combined bin (ncombined entries).
   */
  static vector<float> combine_energies( const vector<float> &sub_energies, const size_t factor )
  {
    assert( factor >= 1 );
    assert( sub_energies.size() >= 2 );

    const size_t raw_count = sub_energies.size() - 1; // number of raw channels
    const size_t ncombined = (raw_count + factor - 1) / factor;

    vector<float> result( ncombined );
    for( size_t i = 0; i < ncombined; ++i )
      result[i] = sub_energies[i * factor];

    return result;
  }//combine_energies(...)


  /** Helper to fill in BackgroundSecondaryInfo for a comparison spectrum.
   *
   * @param comp The comparison spectrum (background or secondary)
   * @param foreground The foreground spectrum
   * @param fore_counts Total foreground counts in the energy range
   * @param first_channel First foreground channel index of the range
   * @param last_channel Last foreground channel index of the range (inclusive)
   * @param factor The 2^N combining factor
   * @param include_channels Whether to include per-channel data
   */
  static SpectrumCountsInEnergyRange::BackgroundSecondaryInfo
  fill_comparison_info( const shared_ptr<const SpecUtils::Measurement> &comp,
                        const shared_ptr<const SpecUtils::Measurement> &foreground,
                        const double fore_counts,
                        const size_t first_channel,
                        const size_t last_channel,
                        const size_t factor,
                        const bool include_channels )
  {
    assert( comp && foreground );

    SpectrumCountsInEnergyRange::BackgroundSecondaryInfo info;

    // Determine live-time scale factor; default to 1.0 if either is invalid
    const float fore_lt = foreground->live_time();
    const float comp_lt = comp->live_time();
    const double lt_sf = ( (fore_lt > 0.0f) && (comp_lt > 0.0f) )
                          ? (static_cast<double>(fore_lt) / static_cast<double>(comp_lt))
                          : 1.0;
    info.live_time_scale_factor = lt_sf;

    // Determine if we need to rebin the comparison spectrum to foreground calibration
    const shared_ptr<const SpecUtils::EnergyCalibration> fore_cal = foreground->energy_calibration();
    const shared_ptr<const SpecUtils::EnergyCalibration> comp_cal = comp->energy_calibration();

    const bool same_cal = ( fore_cal && comp_cal && (*fore_cal == *comp_cal) );

    // Get the comparison channel counts aligned to foreground channels
    vector<float> comp_sub_counts;
    if( same_cal )
    {
      // Same calibration: extract directly from same channel indices
      const shared_ptr<const vector<float>> &comp_counts = comp->gamma_counts();
      if( comp_counts && (comp_counts->size() > last_channel) )
      {
        comp_sub_counts.assign( comp_counts->begin() + first_channel,
                                comp_counts->begin() + last_channel + 1 );
      }
    }
    else
    {
      // Different calibration: rebin entire comparison spectrum into foreground channels,
      // then extract sub-range.  We rebin the full spectrum because rebin_by_lower_edge
      // requires at least 4 channels in both inputs.
      const shared_ptr<const vector<float>> &comp_counts = comp->gamma_counts();
      const shared_ptr<const vector<float>> &comp_energies = comp->gamma_channel_energies();
      const shared_ptr<const vector<float>> &fore_energies = foreground->gamma_channel_energies();

      if( comp_counts && comp_energies && fore_energies
          && (comp_energies->size() >= 4) && (fore_energies->size() >= 4) )
      {
        vector<float> rebinned;
        SpecUtils::rebin_by_lower_edge( *comp_energies, *comp_counts, *fore_energies, rebinned );

        if( rebinned.size() > last_channel )
        {
          comp_sub_counts.assign( rebinned.begin() + first_channel,
                                  rebinned.begin() + last_channel + 1 );
        }
      }//if( valid calibration data for rebinning )
    }//if( same_cal ) / else

    // Compute total comparison counts in range
    double comp_total = 0.0;
    for( const float c : comp_sub_counts )
      comp_total += std::max( c, 0.0f );

    info.counts = comp_total;
    info.cps = (comp_lt > 0.0f) ? (comp_total / static_cast<double>(comp_lt)) : std::numeric_limits<double>::quiet_NaN();

    // Statistical significance (Poisson)
    const double scaled_comp = comp_total * lt_sf;
    const double comp_sigma = sqrt( comp_total );
    const double fore_sigma = sqrt( std::max(fore_counts, 0.0) );
    const double scaled_comp_sigma = lt_sf * comp_sigma;
    const double total_sigma = sqrt( scaled_comp_sigma * scaled_comp_sigma + fore_sigma * fore_sigma );
    info.num_sigma_rel_foreground = (total_sigma > 0.0) ? ((fore_counts - scaled_comp) / total_sigma) : 0.0;

    // Per-channel data
    if( include_channels && !comp_sub_counts.empty() )
    {
      // Combine channels if needed
      vector<float> combined = combine_channels( comp_sub_counts, factor );

      // Scale by live-time ratio
      for( float &val : combined )
        val = static_cast<float>( static_cast<double>(val) * lt_sf );

      info.channel_counts = std::move( combined );
    }

    return info;
  }//fill_comparison_info(...)


  SpectrumCountsInEnergyRange get_counts_in_energy_range( double lower_energy, double upper_energy,
                                                          int max_channels, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error( "No InterSpec session available" );

    const shared_ptr<const SpecUtils::Measurement> foreground
      = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
    if( !foreground )
      throw runtime_error( "get_counts_in_energy_range: no foreground loaded." );

    if( lower_energy > upper_energy )
      std::swap( lower_energy, upper_energy );

    const shared_ptr<const vector<float>> &channel_energies = foreground->gamma_channel_energies();
    const shared_ptr<const vector<float>> &gamma_counts = foreground->gamma_counts();
    const size_t num_channels = foreground->num_gamma_channels();

    if( !channel_energies || !gamma_counts || (num_channels < 2) || (channel_energies->size() < 3) )
      throw runtime_error( "get_counts_in_energy_range: foreground has no valid channel data." );

    // channel_energies has (num_channels + 1) entries: lower edge of each channel + final upper edge
    assert( channel_energies->size() == (num_channels + 1) );

    // Snap to channel boundaries
    size_t first_channel = foreground->find_gamma_channel( static_cast<float>(lower_energy) );
    size_t last_channel = foreground->find_gamma_channel( static_cast<float>(upper_energy) );

    // Make sure we include at least one channel
    if( last_channel < first_channel )
      last_channel = first_channel;

    SpectrumCountsInEnergyRange answer;
    answer.lower_energy = (*channel_energies)[first_channel];
    answer.upper_energy = (*channel_energies)[last_channel + 1];

    // Compute foreground total for the range
    answer.foreground_counts = foreground->gamma_channels_sum( first_channel, last_channel );
    if( foreground->live_time() > 0.0f )
      answer.foreground_cps = answer.foreground_counts / foreground->live_time();
    else
      answer.foreground_cps = std::numeric_limits<double>::quiet_NaN();

    // Treat max_channels <= 1 as a request for just gross counts (single bin)
    if( max_channels < 1 )
      max_channels = 1;

    const bool include_channels = (max_channels > 1);

    if( include_channels )
    {
      size_t raw_count = last_channel - first_channel + 1;

      // Find smallest 2^N combining factor so output has <= max_channels bins
      size_t factor = 1;
      while( (raw_count / factor) > static_cast<size_t>(max_channels) )
        factor *= 2;

      // Expand range so raw_count is divisible by factor
      if( (factor > 1) && ((raw_count % factor) != 0) )
      {
        const size_t deficit = factor - (raw_count % factor);
        size_t expand_low = deficit / 2;
        size_t expand_high = deficit - expand_low;

        // Apply low expansion with clamping
        const size_t actual_expand_low = std::min( expand_low, first_channel );
        first_channel -= actual_expand_low;
        const size_t leftover_low = expand_low - actual_expand_low;

        // Apply high expansion (absorb leftover from low side)
        const size_t total_expand_high = expand_high + leftover_low;
        const size_t max_expand_high = (num_channels - 1) - last_channel;
        const size_t actual_expand_high = std::min( total_expand_high, max_expand_high );
        last_channel += actual_expand_high;
        const size_t leftover_high = total_expand_high - actual_expand_high;

        // If there's still leftover (pinned at upper boundary), try expanding low more
        if( leftover_high > 0 )
        {
          const size_t extra_low = std::min( leftover_high, first_channel );
          first_channel -= extra_low;
        }

        // Update raw_count and actual energies
        raw_count = last_channel - first_channel + 1;

#if( PERFORM_DEVELOPER_CHECKS )
        if( (first_channel > 0) && ((last_channel + 1) < num_channels) )
        {
          // If we're not pinned at both boundaries, divisibility should hold
          if( (raw_count % factor) != 0 )
          {
            const string msg = "get_counts_in_energy_range: raw_count=" + std::to_string(raw_count)
                               + " not divisible by factor=" + std::to_string(factor)
                               + " (first_channel=" + std::to_string(first_channel)
                               + ", last_channel=" + std::to_string(last_channel)
                               + ", num_channels=" + std::to_string(num_channels) + ")";
            log_developer_error( __func__, msg.c_str() );
          }
        }
#endif
      }//if( need to expand for divisibility )

      answer.lower_energy = (*channel_energies)[first_channel];
      answer.upper_energy = (*channel_energies)[last_channel + 1];
      answer.channel_combine_factor = factor;

      // Recompute foreground total for the (possibly expanded) range
      answer.foreground_counts = foreground->gamma_channels_sum( first_channel, last_channel );
      if( foreground->live_time() > 0.0f )
        answer.foreground_cps = answer.foreground_counts / foreground->live_time();

      // Extract sub-vectors for the range
      const vector<float> sub_energies( channel_energies->begin() + first_channel,
                                        channel_energies->begin() + last_channel + 2 );
      const vector<float> sub_counts( gamma_counts->begin() + first_channel,
                                      gamma_counts->begin() + last_channel + 1 );

      answer.channel_lower_energies = combine_energies( sub_energies, factor );
      answer.foreground_channel_counts = combine_channels( sub_counts, factor );
    }//if( include_channels )

    // Handle background
    const shared_ptr<const SpecUtils::Measurement> background
      = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
    if( background )
    {
      answer.background_info = fill_comparison_info( background, foreground,
        answer.foreground_counts, first_channel, last_channel,
        include_channels ? answer.channel_combine_factor : 1, include_channels );
    }

    // Handle secondary
    const shared_ptr<const SpecUtils::Measurement> secondary
      = interspec->displayedHistogram( SpecUtils::SpectrumType::SecondForeground );
    if( secondary )
    {
      answer.secondary_info = fill_comparison_info( secondary, foreground,
        answer.foreground_counts, first_channel, last_channel,
        include_channels ? answer.channel_combine_factor : 1, include_channels );
    }

    return answer;
  }//get_counts_in_energy_range(...)


  const char* to_string( EscapePeakType type )
  {
    switch( type )
    {
      case EscapePeakType::SingleEscape: return "SingleEscape";
      case EscapePeakType::DoubleEscape: return "DoubleEscape";
    }//switch( type )

    throw runtime_error( "Invalid EscapePeakType value" );
  }//const char* to_string( EscapePeakType type )


  const char* to_string( SumPeakType type )
  {
    switch( type )
    {
      case SumPeakType::NotASumPeak: return "NotASumPeak";
      case SumPeakType::RandomSum:   return "RandomSum";
      case SumPeakType::CascadeSum:  return "CascadeSum";
      case SumPeakType::CascadeAndRandomSum: return "CascadeAndRandomSum";
      case SumPeakType::Unknown:     return "Unknown";
    }//switch( type )

    throw runtime_error( "Invalid SumPeakType value" );
  }//const char* to_string( SumPeakType type )


  const char *to_string( EditPeakStructuralAction action )
  {
    switch( action )
    {
      case EditPeakStructuralAction::DeletePeak:    return "DeletePeak";
      case EditPeakStructuralAction::SplitFromRoi:  return "SplitFromRoi";
      case EditPeakStructuralAction::MergeWithLeft:  return "MergeWithLeft";
      case EditPeakStructuralAction::MergeWithRight: return "MergeWithRight";
    }
    throw runtime_error( "Invalid EditPeakStructuralAction value" );
  }//to_string( EditPeakStructuralAction )


  EditPeakStructuralAction structural_action_from_string( const std::string &str )
  {
    if( str == "DeletePeak" )    return EditPeakStructuralAction::DeletePeak;
    if( str == "SplitFromRoi" )  return EditPeakStructuralAction::SplitFromRoi;
    if( str == "MergeWithLeft" ) return EditPeakStructuralAction::MergeWithLeft;
    if( str == "MergeWithRight" )return EditPeakStructuralAction::MergeWithRight;
    throw runtime_error( "Invalid structural action: " + str );
  }//structural_action_from_string(...)


  /** Helper to parse a continuum type string, with "None" as an alias for "NoOffset". */
  PeakContinuum::OffsetType parse_continuum_type_str( const std::string &str )
  {
    if( str == "None" || str == "NoOffset" )
      return PeakContinuum::OffsetType::NoOffset;

    // str_to_offset_type_str expects a char* and length; it handles the standard names
    return PeakContinuum::str_to_offset_type_str( str.c_str(), str.size() );
  }//parse_continuum_type_str(...)


  EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "edit_analysis_peak: No InterSpec session available" );

    // Get the PeakModel for the specified spectrum type
    const std::shared_ptr<PeakModel> peakModel = interspec->peakModel();
    if( !peakModel )
      throw runtime_error( "edit_analysis_peak: No peak model available" );

    // Find nearest peak
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( options.energy );
    if( !peak )
    {
      return EditAnalysisPeakStatus{
        false,
        "No peak found near " + std::to_string( options.energy ) + " keV",
        false,
        std::nullopt,
        {}
      };
    }

    // Validate peak is within 1 FWHM (or 1 keV, for energies given to the nearest keV) of the
    //  requested energy
    const double max_energy_diff = std::max( 1.0, static_cast<double>( get_expected_fwhm( options.energy, interspec ) ) );
    const double energy_diff = std::fabs( peak->mean() - options.energy );
    if( energy_diff > max_energy_diff )
    {
      return EditAnalysisPeakStatus{
        false,
        "Nearest peak at " + std::to_string( peak->mean() ) + " keV is more than 1 FWHM ("
          + std::to_string( max_energy_diff ) + " keV) away from requested "
          + std::to_string( options.energy ) + " keV",
        false,
        std::nullopt,
        {}
      };
    }

    // ---- Handle structural actions ----
    if( options.structuralAction.has_value() )
    {
      const EditPeakStructuralAction action = *options.structuralAction;

      if( action == EditPeakStructuralAction::DeletePeak )
      {
        peakModel->removePeak( peak );
        return EditAnalysisPeakStatus{
          true,
          "Deleted peak at " + std::to_string( peak->mean() ) + " keV",
          false,
          std::nullopt,
          {}
        };
      }//if( DeletePeak )

      const std::shared_ptr<const PeakContinuum> originalContinuum = peak->continuum();

      if( action == EditPeakStructuralAction::SplitFromRoi )
      {
        // Create independent copy with its own continuum; siblings keep original continuum.
        PeakDef modifiedPeak = *peak;
        modifiedPeak.makeUniqueNewContinuum();
        peakModel->updatePeak( peak, modifiedPeak );

        return EditAnalysisPeakStatus{
          true,
          "Split peak at " + std::to_string( peak->mean() ) + " keV into its own ROI",
          false,
          make_shared<const PeakDef>( modifiedPeak ),
          { make_shared<const PeakDef>( modifiedPeak ) }
        };
      }//if( SplitFromRoi )

      // MergeWithLeft / MergeWithRight
      {
        shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
        if( !all_peaks_ptr )
          throw runtime_error( "No peaks available in model" );

        const deque<shared_ptr<const PeakDef>> &all_peaks = *all_peaks_ptr;

        assert( std::is_sorted( all_peaks.begin(), all_peaks.end(),
          []( const shared_ptr<const PeakDef> &a, const shared_ptr<const PeakDef> &b ){
            return a->mean() < b->mean();
          } ) );

        const deque<shared_ptr<const PeakDef>>::const_iterator peak_iter
          = find( all_peaks.begin(), all_peaks.end(), peak );
        if( peak_iter == all_peaks.end() )
          throw runtime_error( "Could not find peak in model" );

        shared_ptr<const PeakDef> adjacent_peak;
        if( action == EditPeakStructuralAction::MergeWithLeft )
        {
          if( peak_iter == all_peaks.begin() )
            throw runtime_error( "No peak to the left to merge with" );
          adjacent_peak = *(peak_iter - 1);
        }
        else
        {
          if( (peak_iter + 1) == all_peaks.end() )
            throw runtime_error( "No peak to the right to merge with" );
          adjacent_peak = *(peak_iter + 1);
        }

        const shared_ptr<const PeakContinuum> adjContinuum = adjacent_peak->continuum();
        if( adjContinuum == originalContinuum )
          throw runtime_error( "Adjacent peak already shares the same ROI" );

        // Collect peaks from both ROIs
        vector<shared_ptr<const PeakDef>> oldPeaksToRemove;
        for( const shared_ptr<const PeakDef> &p : all_peaks )
        {
          if( (p->continuum() == originalContinuum) || (p->continuum() == adjContinuum) )
            oldPeaksToRemove.push_back( p );
        }

        // Create merged continuum with union of ranges
        std::shared_ptr<PeakContinuum> mergedContinuum
          = make_shared<PeakContinuum>( *originalContinuum );
        const double newLower = (std::min)( originalContinuum->lowerEnergy(),
                                            adjContinuum->lowerEnergy() );
        const double newUpper = (std::max)( originalContinuum->upperEnergy(),
                                            adjContinuum->upperEnergy() );
        mergedContinuum->setRange( newLower, newUpper );

        vector<PeakDef> mergedPeaks;
        for( const shared_ptr<const PeakDef> &p : oldPeaksToRemove )
        {
          PeakDef peakCopy = *p;
          peakCopy.setContinuum( mergedContinuum );
          mergedPeaks.push_back( peakCopy );
        }

        // Refit the merged ROI
        const shared_ptr<const SpecUtils::Measurement> data
          = interspec->displayedHistogram( options.specType );
        shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
        const shared_ptr<const DetectorPeakResponse> det = meas ? meas->detector() : nullptr;

        bool used_refit = false;
        if( data )
        {
          vector<shared_ptr<const PeakDef>> peaksToFit;
          for( const PeakDef &p : mergedPeaks )
            peaksToFit.push_back( make_shared<const PeakDef>( p ) );

          shared_ptr<const PeakFitDetPrefs> refitPrefs = meas ? meas->peakFitDetPrefs() : nullptr;
          if( !refitPrefs && det )
            refitPrefs = det->peakFitDetPrefs();
          const PeakFitUtils::CoarseResolutionType refitDetType
            = refitPrefs ? refitPrefs->m_det_type : PeakFitUtils::coarse_det_type( data, nullptr );

          const vector<shared_ptr<const PeakDef>> refitResult
            = refitPeaksThatShareROI( data, det, peaksToFit, refitDetType,
                                      Wt::WFlags<PeakFitLM::PeakFitLMOptions>{} );

          if( refitResult.size() == mergedPeaks.size() )
          {
            mergedPeaks.clear();
            for( const shared_ptr<const PeakDef> &p : refitResult )
              mergedPeaks.push_back( *p );
            used_refit = true;
          }
        }//if( data )

        peakModel->updatePeaks( oldPeaksToRemove, mergedPeaks );

        // Collect results
        const std::shared_ptr<const PeakContinuum> resultContinuum
          = mergedPeaks.empty() ? mergedContinuum : mergedPeaks.front().continuum();

        vector<shared_ptr<const PeakDef>> peaks_in_roi;
        shared_ptr<const PeakDef> updated_peak;

        all_peaks_ptr = peakModel->peaks();
        if( all_peaks_ptr )
        {
          for( const shared_ptr<const PeakDef> &p : *all_peaks_ptr )
          {
            if( p->continuum() == resultContinuum )
            {
              peaks_in_roi.push_back( p );
              if( !updated_peak && fabs( p->mean() - peak->mean() ) < 0.5 )
                updated_peak = p;
            }
          }
        }

        if( !updated_peak && !peaks_in_roi.empty() )
          updated_peak = peaks_in_roi.front();

        string msg = "Merged ROIs into combined region [" + std::to_string( newLower )
          + ", " + std::to_string( newUpper ) + "] keV with "
          + std::to_string( peaks_in_roi.size() ) + " peaks";

        if( used_refit )
          msg += " (refit to data)";
        else
          msg += ". Warning: refit failed or peaks became insignificant; un-fit parameters used.";

        const double avg_fwhm = 0.5 * (peak->fwhm() + adjacent_peak->fwhm());
        const double peak_separation = fabs( peak->mean() - adjacent_peak->mean() );
        if( peak_separation > 6.0 * avg_fwhm )
        {
          msg += " Warning: peaks are " + std::to_string( peak_separation / avg_fwhm )
            + " FWHM apart; merging may not be appropriate.";
        }

        return EditAnalysisPeakStatus{
          true,
          msg,
          used_refit,
          updated_peak,
          peaks_in_roi
        };
      }//MergeWithLeft / MergeWithRight
    }//if( structuralAction )


    // ---- Property-based editing ----
    try
    {
      // Create modified copy and collect ROI siblings
      PeakDef modifiedPeak = *peak;
      const std::shared_ptr<const PeakContinuum> originalContinuum = peak->continuum();

      vector<shared_ptr<const PeakDef>> roiSiblings;
      {
        shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
        if( all_peaks_ptr )
        {
          for( const shared_ptr<const PeakDef> &p : *all_peaks_ptr )
          {
            if( (p != peak) && (p->continuum() == originalContinuum) )
              roiSiblings.push_back( p );
          }
        }
      }

      modifiedPeak.makeUniqueNewContinuum();

      // Track whether shape-affecting values changed (for auto-refit determination)
      bool shapeValueChanged = false;
      bool typeChanged = false; // skew or continuum type changed (triggers full refit)

      // ---- Apply skew type ----
      if( options.skewType.has_value() )
      {
        const PeakDef::SkewType skew_type = PeakDef::skew_from_string( *options.skewType );
        const bool skew_type_changed = (skew_type != modifiedPeak.skewType());
        modifiedPeak.setSkewType( skew_type );

        // `setSkewType` doesnt touch the coefficients, so the peak keeps the previous skew types
        //  values - which mean something else for the new type (e.g. a DoubleSidedCrystalBall `n`
        //  would become a GADRAS power).  So when the type changes, seed the new types values from
        //  the peak-fit prefs (e.g., a GADRAS detector's tail powers/extents) or else its defaults,
        //  with its default fit-for flags; this is what #PeakEdit does when the user changes skew
        //  type.  An explicitly requested `skewParN`, below, overrides this.
        if( skew_type_changed )
        {
          const shared_ptr<const SpecMeas> skew_meas = interspec->measurment( options.specType );
          const shared_ptr<const PeakFitDetPrefs> meas_prefs = skew_meas ? skew_meas->peakFitDetPrefs() : nullptr;
          const shared_ptr<const DetectorPeakResponse> skew_drf = skew_meas ? skew_meas->detector() : nullptr;
          const shared_ptr<const PeakFitDetPrefs> drf_prefs = skew_drf ? skew_drf->peakFitDetPrefs() : nullptr;

          const size_t num_skew_par = PeakDef::num_skew_parameters( skew_type );
          for( size_t i = 0; i < num_skew_par; ++i )
          {
            const PeakDef::CoefficientType ct
                            = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );
            modifiedPeak.set_coefficient( skew_starting_value( skew_type, ct, meas_prefs.get(), drf_prefs.get() ), ct );
            modifiedPeak.setFitFor( ct, PeakDef::skew_parameter_fit_by_default( skew_type, ct ) );
          }//for( each skew parameter of the new skew type )
        }//if( skew_type_changed )

        typeChanged = true;
      }

      // ---- Apply skew parameter values ----
      const std::optional<double> *skew_pars[] = {
        &options.skewPar0, &options.skewPar1, &options.skewPar2,
        &options.skewPar3, &options.skewPar4, &options.skewPar5
      };
      const PeakDef::CoefficientType skew_coef_types[] = {
        PeakDef::CoefficientType::SkewPar0, PeakDef::CoefficientType::SkewPar1,
        PeakDef::CoefficientType::SkewPar2, PeakDef::CoefficientType::SkewPar3,
        PeakDef::CoefficientType::SkewPar4, PeakDef::CoefficientType::SkewPar5
      };
      static_assert( std::size(skew_pars) == std::size(skew_coef_types) );
      const size_t num_skew = PeakDef::num_skew_parameters( modifiedPeak.skewType() );

      for( size_t i = 0; i < std::size(skew_pars); ++i )
      {
        if( skew_pars[i]->has_value() )
        {
          if( i >= num_skew )
          {
            throw runtime_error( "skewPar" + std::to_string( i ) + " is not valid for skew type "
              + string( PeakDef::to_string( modifiedPeak.skewType() ) )
              + " (has " + std::to_string( num_skew ) + " parameters)" );
          }
          modifiedPeak.set_coefficient( *(*skew_pars[i]), skew_coef_types[i] );
          shapeValueChanged = true;
        }
      }

      // ---- Apply skew fit-for flags ----
      const std::optional<bool> *skew_fit_flags[] = {
        &options.fitForSkewPar0, &options.fitForSkewPar1, &options.fitForSkewPar2,
        &options.fitForSkewPar3, &options.fitForSkewPar4, &options.fitForSkewPar5
      };
      static_assert( std::size(skew_fit_flags) == std::size(skew_coef_types) );
      for( size_t i = 0; i < std::size(skew_fit_flags); ++i )
      {
        if( skew_fit_flags[i]->has_value() )
        {
          if( i >= num_skew )
          {
            throw runtime_error( "fitForSkewPar" + std::to_string( i )
              + " is not valid for skew type "
              + string( PeakDef::to_string( modifiedPeak.skewType() ) ) );
          }
          modifiedPeak.setFitFor( skew_coef_types[i], *(*skew_fit_flags[i]) );
        }
      }

      // ---- Apply continuum type ----
      shared_ptr<PeakContinuum> continuum = modifiedPeak.getContinuum();

      if( options.continuumType.has_value() )
      {
        const PeakContinuum::OffsetType cont_type = parse_continuum_type_str( *options.continuumType );
        continuum->setType( cont_type );
        typeChanged = true;
      }

      // ---- Apply continuum coefficient values ----
      const std::optional<double> *cont_coefs[] = {
        &options.continuumCoef0, &options.continuumCoef1,
        &options.continuumCoef2, &options.continuumCoef3
      };
      const size_t num_cont_pars = PeakContinuum::num_parameters( continuum->type() );

      for( size_t i = 0; i < 4; ++i )
      {
        if( cont_coefs[i]->has_value() )
        {
          if( i >= num_cont_pars )
          {
            throw runtime_error( "continuumCoef" + std::to_string( i )
              + " is not valid for continuum type "
              + string( PeakContinuum::offset_type_name( continuum->type() ) )
              + " (has " + std::to_string( num_cont_pars ) + " parameters)" );
          }
          continuum->setPolynomialCoef( i, *(*cont_coefs[i]) );
          shapeValueChanged = true;
        }
      }

      // ---- Apply continuum fit-for flags ----
      const std::optional<bool> *cont_fit_flags[] = {
        &options.fitForContinuumCoef0, &options.fitForContinuumCoef1,
        &options.fitForContinuumCoef2, &options.fitForContinuumCoef3
      };
      for( size_t i = 0; i < 4; ++i )
      {
        if( cont_fit_flags[i]->has_value() )
        {
          if( i >= num_cont_pars )
          {
            throw runtime_error( "fitForContinuumCoef" + std::to_string( i )
              + " is not valid for continuum type "
              + string( PeakContinuum::offset_type_name( continuum->type() ) ) );
          }
          continuum->setPolynomialCoefFitFor( i, *(*cont_fit_flags[i]) );
        }
      }

      // ---- Apply ROI bounds ----
      if( options.roiLowerEnergy.has_value() || options.roiUpperEnergy.has_value() )
      {
        const double new_lower = options.roiLowerEnergy.value_or( continuum->lowerEnergy() );
        const double new_upper = options.roiUpperEnergy.value_or( continuum->upperEnergy() );

        if( new_lower >= new_upper )
          throw runtime_error( "ROI lower energy must be less than upper energy" );

        continuum->setRange( new_lower, new_upper );
        shapeValueChanged = true;
      }

      // ---- Apply Gaussian properties ----
      if( options.meanEnergy.has_value() )
      {
        modifiedPeak.setMean( *options.meanEnergy );
        shapeValueChanged = true;
      }

      if( options.fwhm.has_value() )
      {
        modifiedPeak.setSigma( (*options.fwhm) / 2.35482 );
        shapeValueChanged = true;
      }

      if( options.amplitude.has_value() )
      {
        modifiedPeak.setAmplitude( *options.amplitude );
        shapeValueChanged = true;
      }

      // ---- Apply Gaussian fit-for flags ----
      if( options.fitForMean.has_value() )
        modifiedPeak.setFitFor( PeakDef::CoefficientType::Mean, *options.fitForMean );
      if( options.fitForFwhm.has_value() )
        modifiedPeak.setFitFor( PeakDef::CoefficientType::Sigma, *options.fitForFwhm );
      if( options.fitForAmplitude.has_value() )
        modifiedPeak.setFitFor( PeakDef::CoefficientType::GaussAmplitude, *options.fitForAmplitude );

      // ---- Apply uncertainties ----
      if( options.meanUncertainty.has_value() )
        modifiedPeak.setMeanUncert( *options.meanUncertainty );
      if( options.fwhmUncertainty.has_value() )
        modifiedPeak.setSigmaUncert( (*options.fwhmUncertainty) / 2.35482 );
      if( options.amplitudeUncertainty.has_value() )
        modifiedPeak.setAmplitudeUncert( *options.amplitudeUncertainty );

      // ---- Apply metadata ----
      if( options.source.has_value() )
      {
        PeakModel::SetGammaSource result
          = PeakModel::setNuclideXrayReaction( modifiedPeak, *options.source, 4.0 );
        if( result == PeakModel::SetGammaSource::FailedSourceChange )
          throw runtime_error( "Failed to set source: " + *options.source );

        // Recolor for the new source, like the peak table does: take the color of the new source's
        //  other peaks, else a new color for the source (or none, if the source was cleared).  A peak
        //  keeps its own color - one without a source, or one differing from the other peaks of its old
        //  source; a lone peak's color is taken to be its source's, so it gets the new source's color.
        const string old_src = peak->sourceName();
        const string new_src = modifiedPeak.sourceName();
        if( !options.color.has_value() && (new_src != old_src) )
        {
          // The distinct colors of the other peaks of a source
          const shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks = peakModel->peaks();
          const auto other_peak_colors = [&all_peaks, &peak]( const string &src ) -> vector<Wt::WColor> {
            vector<Wt::WColor> colors;
            for( size_t i = 0; !src.empty() && all_peaks && (i < all_peaks->size()); ++i )
            {
              const shared_ptr<const PeakDef> &p = (*all_peaks)[i];
              const Wt::WColor &c = p->lineColor();
              if( (p != peak) && (p->sourceName() == src) && !c.isDefault()
                 && std::none_of( begin(colors), end(colors), [&c]( const Wt::WColor &o ){ return same_rgb( c, o ); } ) )
              {
                colors.push_back( c );
              }
            }
            return colors;
          };//other_peak_colors lambda

          const Wt::WColor &old_color = peak->lineColor();
          const vector<Wt::WColor> old_src_colors = other_peak_colors( old_src );
          const bool matches_old_src = std::any_of( begin(old_src_colors), end(old_src_colors),
                                         [&old_color]( const Wt::WColor &c ){ return same_rgb( c, old_color ); } );
          const bool own_color = !old_color.isDefault()
                                 && (old_src.empty() || (!old_src_colors.empty() && !matches_old_src));

          if( !own_color )
          {
            const vector<Wt::WColor> new_src_colors = other_peak_colors( new_src );
            map<string,Wt::WColor> pending_colors;
            const Wt::WColor color = (new_src_colors.size() == 1)
                                       ? new_src_colors[0]
                                       : color_for_new_peak( modifiedPeak, interspec, pending_colors );

            // If there is no color for the new source (e.g., the color theme doesnt color peaks by
            //  source), a lone peak keeps its color, as in the peak table.
            if( !color.isDefault() || new_src.empty() || !old_src_colors.empty() )
              modifiedPeak.setLineColor( color );
          }//if( !own_color )
        }//if( source changed, and no color given )
      }

      if( options.color.has_value() )
        modifiedPeak.setLineColor( Wt::WColor( *options.color ) );

      if( options.userLabel.has_value() )
        modifiedPeak.setUserLabel( *options.userLabel );

      // ---- Apply usage flags ----
      if( options.useForEnergyCalibration.has_value() )
        modifiedPeak.useForEnergyCalibration( *options.useForEnergyCalibration );
      if( options.useForShieldingSourceFit.has_value() )
        modifiedPeak.useForShieldingSourceFit( *options.useForShieldingSourceFit );
      if( options.useForManualRelEff.has_value() )
        modifiedPeak.useForManualRelEff( *options.useForManualRelEff );


      // ---- Determine whether to refit ----
      bool doRefit = false;
      if( options.refit.has_value() )
        doRefit = *options.refit;
      else
        doRefit = shapeValueChanged || typeChanged;

      bool refitPerformed = false;
      string refitWarning;

      if( doRefit )
      {
        const shared_ptr<const SpecUtils::Measurement> data
          = interspec->displayedHistogram( options.specType );
        shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
        const shared_ptr<const DetectorPeakResponse> det = meas ? meas->detector() : nullptr;

        if( data )
        {
          // Build vector of all peaks in the ROI for refitting
          vector<shared_ptr<const PeakDef>> peaksToFit;
          peaksToFit.push_back( make_shared<const PeakDef>( modifiedPeak ) );

          const std::shared_ptr<PeakContinuum> newContinuum = modifiedPeak.getContinuum();
          for( const shared_ptr<const PeakDef> &sibling : roiSiblings )
          {
            PeakDef siblingCopy = *sibling;
            siblingCopy.setContinuum( newContinuum );
            peaksToFit.push_back( make_shared<const PeakDef>( siblingCopy ) );
          }

          // Use full refit for type changes or when explicitly forced with no shape changes.
          // Use medium refinement when parameter values were explicitly set (to stay close to them).
          Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options;
          if( !typeChanged && shapeValueChanged )
            fit_options |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;

          shared_ptr<const PeakFitDetPrefs> refitPrefs = meas ? meas->peakFitDetPrefs() : nullptr;
          if( !refitPrefs && det )
            refitPrefs = det->peakFitDetPrefs();
          const PeakFitUtils::CoarseResolutionType refitDetType
            = refitPrefs ? refitPrefs->m_det_type : PeakFitUtils::coarse_det_type( data, nullptr );

          const vector<shared_ptr<const PeakDef>> refitResult
            = refitPeaksThatShareROI( data, det, peaksToFit, refitDetType, fit_options );

          if( refitResult.size() == peaksToFit.size() )
          {
            // Refit succeeded — find the target peak and siblings in the result
            for( const shared_ptr<const PeakDef> &p : refitResult )
            {
              if( fabs( p->mean() - modifiedPeak.mean() ) < 0.5 )
                modifiedPeak = *p;
            }

            // Update roiSiblings with refitted versions
            roiSiblings.clear();
            for( const shared_ptr<const PeakDef> &p : refitResult )
            {
              if( fabs( p->mean() - modifiedPeak.mean() ) >= 0.5 )
              {
                // These are siblings — we need shared_ptr<const PeakDef> for updatePeaks
                // We'll rebuild them below when updating the model
              }
            }

            // All refit peaks share the same continuum — use them directly
            // Rebuild the modified peak and sibling list from refitResult
            vector<shared_ptr<const PeakDef>> oldPeaks;
            vector<PeakDef> newPeaks;

            oldPeaks.push_back( peak );
            // Find our target peak in refitResult
            bool found_target = false;
            for( const shared_ptr<const PeakDef> &p : refitResult )
            {
              if( !found_target && fabs( p->mean() - peak->mean() ) < 0.5 )
              {
                modifiedPeak = *p;
                newPeaks.push_back( *p );
                found_target = true;
              }
            }
            if( !found_target )
            {
              // Shouldn't happen, but fall back to our un-refit peak
              newPeaks.push_back( modifiedPeak );
            }

            // Now handle siblings
            // Build a map of original sibling peaks to their refit versions
            shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
            if( all_peaks_ptr )
            {
              for( const shared_ptr<const PeakDef> &origSibling : *all_peaks_ptr )
              {
                if( (origSibling != peak) && (origSibling->continuum() == originalContinuum) )
                {
                  oldPeaks.push_back( origSibling );

                  // Find matching refit peak for this sibling
                  bool found_sibling = false;
                  for( const shared_ptr<const PeakDef> &refitPeak : refitResult )
                  {
                    if( fabs( refitPeak->mean() - origSibling->mean() ) < 0.5
                        && fabs( refitPeak->mean() - modifiedPeak.mean() ) >= 0.01 )
                    {
                      newPeaks.push_back( *refitPeak );
                      found_sibling = true;
                      break;
                    }
                  }

                  if( !found_sibling )
                  {
                    // Fall back: copy sibling with new continuum
                    PeakDef siblingCopy = *origSibling;
                    siblingCopy.setContinuum( modifiedPeak.getContinuum() );
                    newPeaks.push_back( siblingCopy );
                  }
                }
              }
            }

            peakModel->updatePeaks( oldPeaks, newPeaks );
            refitPerformed = true;

            // Build return value
            vector<shared_ptr<const PeakDef>> peaks_in_roi;
            shared_ptr<const PeakDef> updated_peak;
            const std::shared_ptr<const PeakContinuum> resultCont = modifiedPeak.continuum();

            all_peaks_ptr = peakModel->peaks();
            if( all_peaks_ptr )
            {
              for( const shared_ptr<const PeakDef> &p : *all_peaks_ptr )
              {
                if( p->continuum() == resultCont )
                {
                  peaks_in_roi.push_back( p );
                  if( !updated_peak && fabs( p->mean() - modifiedPeak.mean() ) < 0.01 )
                    updated_peak = p;
                }
              }
            }

            if( !updated_peak )
              updated_peak = make_shared<const PeakDef>( modifiedPeak );

            return EditAnalysisPeakStatus{
              true,
              "Edited and refit peak at " + std::to_string( modifiedPeak.mean() ) + " keV",
              true,
              updated_peak,
              peaks_in_roi
            };
          }
          else
          {
            // Refit returned different number of peaks — fall back to un-refit
            refitWarning = " Warning: refit was attempted but one or more peaks became"
              " insignificant; un-refit parameters were used instead.";
          }
        }
        else
        {
          refitWarning = " Warning: no spectrum data available for refitting.";
        }
      }//if( doRefit )


      // ---- Update peaks in model (no refit, or refit failed) ----
      if( roiSiblings.empty() )
      {
        peakModel->updatePeak( peak, modifiedPeak );
      }
      else
      {
        vector<shared_ptr<const PeakDef>> oldPeaks;
        vector<PeakDef> newPeaks;

        oldPeaks.push_back( peak );
        newPeaks.push_back( modifiedPeak );

        const std::shared_ptr<PeakContinuum> newContinuum = modifiedPeak.getContinuum();
        for( const shared_ptr<const PeakDef> &sibling : roiSiblings )
        {
          oldPeaks.push_back( sibling );

          PeakDef siblingCopy = *sibling;
          siblingCopy.setContinuum( newContinuum );
          newPeaks.push_back( siblingCopy );
        }

        peakModel->updatePeaks( oldPeaks, newPeaks );
      }

      // ---- Build return value ----
      vector<shared_ptr<const PeakDef>> peaks_in_roi;
      shared_ptr<const PeakDef> updated_peak;
      {
        shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
        if( all_peaks_ptr )
        {
          const std::shared_ptr<const PeakContinuum> resultCont = modifiedPeak.continuum();
          for( const shared_ptr<const PeakDef> &p : *all_peaks_ptr )
          {
            if( p->continuum() == resultCont )
            {
              peaks_in_roi.push_back( p );
              if( !updated_peak && fabs( p->mean() - modifiedPeak.mean() ) < 0.01 )
                updated_peak = p;
            }
          }
        }
      }

      if( !updated_peak )
        updated_peak = make_shared<const PeakDef>( modifiedPeak );

      string msg = "Edited peak at " + std::to_string( modifiedPeak.mean() ) + " keV";
      if( !refitWarning.empty() )
        msg += refitWarning;

      return EditAnalysisPeakStatus{
        true,
        msg,
        refitPerformed,
        updated_peak,
        peaks_in_roi
      };
    }
    catch( const exception &e )
    {
      return EditAnalysisPeakStatus{
        false,
        string( "Error editing peak: " ) + e.what(),
        false,
        std::nullopt,
        {}
      };
    }
  }//EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec )


  EscapePeakCheckStatus escape_peak_check( const EscapePeakCheckOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "escape_peak_check: No InterSpec session available" );

    const double input_energy = options.energy;

    // Constants for escape peak calculations
    const double electron_rest_mass = 510.9989; // keV
    const double single_escape_offset = electron_rest_mass;
    const double double_escape_offset = 2.0 * electron_rest_mass;

    // Get the spectrum
    shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
    if( !meas )
    {
      throw runtime_error( "No measurement loaded for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );
    }

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram( options.specType );
    if( !spectrum )
    {
      throw runtime_error( "No spectrum displayed for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );
    }

    // Determine detector type for pair production threshold
    const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );
    const double pair_prod_thresh = isHPGe ? 1255.0 : 2585.0;  // escape peaks are rarely visible below these

    // Initialize result
    EscapePeakCheckStatus result;

    const double spectrum_min_energy = std::max( 0.0, static_cast<double>( spectrum->gamma_energy_min() ) );
    const double spectrum_max_energy = spectrum->gamma_energy_max();

    // Calculate potential energies.  This peak can only have escape peaks of its own if it is above
    //  the pair-production threshold, and the escape energies are in the spectrum.
    result.belowPairProductionThreshold = (input_energy < double_escape_offset);
    if( !result.belowPairProductionThreshold )
    {
      if( (input_energy - single_escape_offset) > spectrum_min_energy )
        result.potentialSingleEscapePeakEnergy = input_energy - single_escape_offset;
      if( (input_energy - double_escape_offset) > spectrum_min_energy )
        result.potentialDoubleEscapePeakEnergy = input_energy - double_escape_offset;
    }
    result.potentialParentPeakSingleEscape = input_energy + single_escape_offset;
    result.potentialParentPeakDoubleEscape = input_energy + double_escape_offset;

    // Calculate search windows for each potential energy
    // Window = max(1.0, min(0.5*fwhm, 15.0)) - fairly arbitrary and untested
    const std::function<float(double)> fwhm_fcn = make_expected_fwhm_fcn( interspec );
    const auto calc_window = [&fwhm_fcn]( const double energy ) -> double {
      if( energy <= 0.0 )
        return 1.0; // If energy is not positive, just use minimum window

      const double fwhm = fwhm_fcn( energy );
      return std::max( 1.0, std::min( 0.5 * fwhm, 15.0 ) );
    };

    // Calculate windows for potential escape peaks of the input energy
    result.singleEscapeSearchWindow = calc_window( input_energy - single_escape_offset );
    result.doubleEscapeSearchWindow = calc_window( input_energy - double_escape_offset );

    // Calculate windows for potential parent peaks (if input is an escape peak)
    result.singleEscapeParentSearchWindow = calc_window( result.potentialParentPeakSingleEscape );
    result.doubleEscapeParentSearchWindow = calc_window( result.potentialParentPeakDoubleEscape );

    // Get all peaks in the spectrum (user + auto-search)
    DetectedPeaksOptions peak_options;
    peak_options.specType = options.specType;
    peak_options.nonBackgroundPeaksOnly = false;
    const DetectedPeakStatus peak_status = detected_peaks( peak_options, interspec );
    const vector<shared_ptr<const PeakDef>> &all_peaks = peak_status.peaks;

    // Helper lambda to find the peak nearest a given energy, within the window
    const auto find_peak_near_energy = [&all_peaks]( const double target_energy, const double window )
        -> shared_ptr<const PeakDef> {
      shared_ptr<const PeakDef> nearest;
      for( const shared_ptr<const PeakDef> &peak : all_peaks )
      {
        if( !peak )
          continue;

        const double diff = fabs( peak->mean() - target_energy );
        if( (diff <= window) && (!nearest || (diff < fabs( nearest->mean() - target_energy ))) )
          nearest = peak;
      }
      return nearest;
    };

    const auto source_of = []( const shared_ptr<const PeakDef> &peak ) -> std::optional<std::string> {
      if( peak->parentNuclide() )
        return peak->parentNuclide()->symbol;
      if( peak->xrayElement() )
        return peak->xrayElement()->symbol + " x-ray";
      if( peak->reaction() )
        return peak->reaction()->name();
      return std::nullopt;
    };

    const auto found_peak = [&find_peak_near_energy, &source_of]( const double energy, const double window )
        -> std::optional<EscapeCheckPeak> {
      const shared_ptr<const PeakDef> peak = find_peak_near_energy( energy, window );
      if( !peak )
        return std::nullopt;
      return EscapeCheckPeak{ peak, source_of( peak ) };
    };

    char note_buffer[512] = { '\0' };

    // This peak as a parent: are its escape peaks there?
    if( result.belowPairProductionThreshold )
    {
      snprintf( note_buffer, sizeof(note_buffer), "%.2f keV is below the 1022 keV pair-production threshold, so"
                " a gamma of this energy produces no single or double escape peaks.", input_energy );
      result.notes.push_back( note_buffer );
    }else
    {
      if( result.potentialSingleEscapePeakEnergy )
        result.singleEscapePeakFound = found_peak( *result.potentialSingleEscapePeakEnergy, result.singleEscapeSearchWindow );
      if( result.potentialDoubleEscapePeakEnergy )
        result.doubleEscapePeakFound = found_peak( *result.potentialDoubleEscapePeakEnergy, result.doubleEscapeSearchWindow );
    }//if( result.belowPairProductionThreshold ) / else

    // This peak as an escape peak: is there a peak where its parent would be?
    result.peakAtPotentialSingleEscapeParent = found_peak( result.potentialParentPeakSingleEscape,
                                                           result.singleEscapeParentSearchWindow );
    result.peakAtPotentialDoubleEscapeParent = found_peak( result.potentialParentPeakDoubleEscape,
                                                           result.doubleEscapeParentSearchWindow );

    if( result.potentialParentPeakSingleEscape < double_escape_offset )
    {
      const std::optional<EscapeCheckPeak> &p = result.peakAtPotentialSingleEscapeParent;
      if( p )
        snprintf( note_buffer, sizeof(note_buffer), "There is a peak at %.2f keV, where a single-escape parent would"
                  " be, but a %.2f keV gamma is below the 1022 keV pair-production threshold, so %.2f keV can not be"
                  " its single escape peak.", p->peak->mean(), p->peak->mean(), input_energy );
      else
        snprintf( note_buffer, sizeof(note_buffer), "A single-escape parent would be at %.2f keV, below the 1022 keV"
                  " pair-production threshold, so %.2f keV can not be a single escape peak.",
                  result.potentialParentPeakSingleEscape, input_energy );
      result.notes.push_back( note_buffer );
    }

    for( const EscapePeakType type : { EscapePeakType::SingleEscape, EscapePeakType::DoubleEscape } )
    {
      const bool single = (type == EscapePeakType::SingleEscape);
      const std::optional<EscapeCheckPeak> &p = single ? result.peakAtPotentialSingleEscapeParent
                                                       : result.peakAtPotentialDoubleEscapeParent;
      if( p && (p->peak->mean() >= double_escape_offset) && (p->peak->mean() <= pair_prod_thresh) )
      {
        snprintf( note_buffer, sizeof(note_buffer), "There is a peak at %.2f keV, where a %s-escape parent would be,"
                  " but escape peaks are rarely visible for gammas below about %.0f keV on this type of detector.",
                  p->peak->mean(), (single ? "single" : "double"), pair_prod_thresh );
        result.notes.push_back( note_buffer );
      }
    }//for( single and double escape )

    // Helper lambda to create label for an escape peak
    const auto create_escape_label = []( const shared_ptr<const PeakDef> &parent_peak,
                                          const EscapePeakType escape_type,
                                          ParentPeakInfo &parent_info ) {
      const char *escape_prefix = (escape_type == EscapePeakType::SingleEscape) ? "S.E." : "D.E.";
      char buffer[256];

      // If the parent peak has a source (nuclide, reaction, or x-ray), use it
      if( parent_peak->parentNuclide() )
      {
        snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                 parent_peak->parentNuclide()->symbol.c_str(),
                 escape_prefix,
                 parent_peak->mean() );
        parent_info.sourceLabel = buffer;
      }
      else if( parent_peak->reaction() )
      {
        snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                 parent_peak->reaction()->name().c_str(),
                 escape_prefix,
                 parent_peak->mean() );
        parent_info.sourceLabel = buffer;
      }
      else if( parent_peak->xrayElement() )
      {
        snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                 parent_peak->xrayElement()->symbol.c_str(),
                 escape_prefix,
                 parent_peak->mean() );
        parent_info.sourceLabel = buffer;
      }
      else
      {
        // No source - create user label
        snprintf( buffer, sizeof(buffer), "%.2f %s", parent_peak->mean(), escape_prefix );
        parent_info.userLabel = buffer;
      }
    };

    // Check if the input energy is an escape peak of another peak
    // IMPORTANT: Check double escape (1022 keV) FIRST to avoid misidentifying
    // a double escape peak as a single escape peak

    // Check for double escape parent (input + 1022 keV)
    const double potential_de_parent = result.potentialParentPeakDoubleEscape;
    if( potential_de_parent > pair_prod_thresh )
    {
      shared_ptr<const PeakDef> parent_peak = find_peak_near_energy(
        potential_de_parent,
        result.doubleEscapeParentSearchWindow
      );

      if( parent_peak )
      {
        ParentPeakInfo parent_info;
        parent_info.parentPeakEnergy = parent_peak->mean();
        parent_info.escapeType = EscapePeakType::DoubleEscape;
        create_escape_label( parent_peak, EscapePeakType::DoubleEscape, parent_info );
        result.parentPeak = parent_info;
      }
    }//if( potential_de_parent > pair_prod_thresh )

    // Check for single escape parent (input + 511 keV) - only if we haven't found a D.E. parent
    if( !result.parentPeak.has_value() )
    {
      const double potential_se_parent = result.potentialParentPeakSingleEscape;
      if( potential_se_parent > pair_prod_thresh )
      {
        shared_ptr<const PeakDef> parent_peak = find_peak_near_energy(
          potential_se_parent,
          result.singleEscapeParentSearchWindow
        );

        if( parent_peak )
        {
          ParentPeakInfo parent_info;
          parent_info.parentPeakEnergy = parent_peak->mean();
          parent_info.escapeType = EscapePeakType::SingleEscape;
          create_escape_label( parent_peak, EscapePeakType::SingleEscape, parent_info );
          result.parentPeak = parent_info;
        }
      }//if( potential_se_parent > pair_prod_thresh )
    }//if( !result.parentPeak.has_value() )

    // Edge case: Check if parent peak would be above detector range
    // This happens when the parent gamma is too high energy to be detected, but its
    // escape peaks are within the detector range
    if( !result.parentPeak.has_value() )
    {
      const double detector_max_energy = spectrum->gamma_energy_max();

      // Collect all unique nuclides and reactions from existing peaks
      set<const SandiaDecay::Nuclide *> nuclides_in_spectrum;
      set<const ReactionGamma::Reaction *> reactions_in_spectrum;

      for( const shared_ptr<const PeakDef> &peak : all_peaks )
      {
        if( !peak )
          continue;

        if( peak->parentNuclide() )
          nuclides_in_spectrum.insert( peak->parentNuclide() );

        if( peak->reaction() )
          reactions_in_spectrum.insert( peak->reaction() );
      }

      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      if( !db )
        throw runtime_error( "No decay database available" );

      // Helper lambda to check if an above-range gamma matches our input as an escape peak
      const auto check_above_range_gamma = [&]( const double gamma_energy,
                                                 const EscapePeakType escape_type,
                                                 const string &source_name ) -> bool {
        if( gamma_energy <= detector_max_energy )
          return false; // Gamma is in range, not an edge case

        if( gamma_energy < pair_prod_thresh )
          return false; // Can't produce escape peaks

        const double expected_escape_energy = (escape_type == EscapePeakType::SingleEscape)
                                               ? (gamma_energy - single_escape_offset)
                                               : (gamma_energy - double_escape_offset);

        const double window = (escape_type == EscapePeakType::SingleEscape)
                              ? result.singleEscapeParentSearchWindow
                              : result.doubleEscapeParentSearchWindow;

        if( fabs(expected_escape_energy - input_energy) <= window )
        {
          // Found a match!
          ParentPeakInfo parent_info;
          parent_info.parentPeakEnergy = gamma_energy;
          parent_info.escapeType = escape_type;

          const char *escape_prefix = (escape_type == EscapePeakType::SingleEscape) ? "S.E." : "D.E.";
          char buffer[256];
          snprintf( buffer, sizeof(buffer), "%s %s %.2f keV",
                   source_name.c_str(), escape_prefix, gamma_energy );
          parent_info.sourceLabel = buffer;

          result.parentPeak = parent_info;
          return true;
        }

        return false;
      };

      // Check nuclides for gammas above detector range
      for( const SandiaDecay::Nuclide *nuc : nuclides_in_spectrum )
      {
        if( !nuc )
          continue;

        const double age = PeakDef::defaultDecayTime( nuc );
        SandiaDecay::NuclideMixture mixture;
        mixture.addNuclide( SandiaDecay::NuclideActivityPair(nuc, 1.0) );
        const vector<SandiaDecay::EnergyRatePair> photons = mixture.photons(
          age,
          SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy
        );

        // Check double escape first, then single escape
        for( const auto &photon : photons )
        {
          const float gamma_energy = photon.energy;

          if( check_above_range_gamma(gamma_energy, EscapePeakType::DoubleEscape, nuc->symbol) )
            break; // Found a match

          if( check_above_range_gamma(gamma_energy, EscapePeakType::SingleEscape, nuc->symbol) )
            break; // Found a match
        }

        if( result.parentPeak.has_value() )
          break; // Found a match, stop searching
      }//for( const SandiaDecay::Nuclide *nuc : nuclides_in_spectrum )

      // Check reactions for gammas above detector range
      if( !result.parentPeak.has_value() )
      {
        for( const ReactionGamma::Reaction *rxn : reactions_in_spectrum )
        {
          if( !rxn )
            continue;

          for( const auto &gamma : rxn->gammas )
          {
            const float gamma_energy = gamma.energy;

            if( check_above_range_gamma(gamma_energy, EscapePeakType::DoubleEscape, rxn->name()) )
              break; // Found a match

            if( check_above_range_gamma(gamma_energy, EscapePeakType::SingleEscape, rxn->name()) )
              break; // Found a match
          }

          if( result.parentPeak.has_value() )
            break; // Found a match, stop searching
        }//for( const ReactionGamma::Reaction *rxn : reactions_in_spectrum )
      }//if( !result.parentPeak.has_value() )
    }//if( !result.parentPeak.has_value() ) - edge case check

    // Detector K x-ray escape: a photon absorbed near the crystal surface whose fluorescence K x-ray
    //  escapes deposits (energy - x-ray energy).  Which x-rays depends on the crystal material.
    struct DetXray
    {
      const char *element;
      const char *line;
      double energy;  // intensity-weighted line energy, keV
      double k_edge;  // keV
    };

    const DetXray ge_lines[] = { {"Ge", "Ka", 9.876, 11.10}, {"Ge", "Kb", 10.98, 11.10} };
    const DetXray i_lines[] = { {"I", "Ka", PeakFitUtils::sm_iodine_kalpha_escape_kev, 33.17},
                                {"I", "Kb", PeakFitUtils::sm_iodine_kbeta_escape_kev, 33.17} };
    const DetXray czt_lines[] = { {"Cd", "Ka", 23.11, 26.71}, {"Cd", "Kb", 26.1, 26.71},
                                  {"Te", "Ka", 27.38, 31.81}, {"Te", "Kb", 31.0, 31.81} };
    const DetXray labr_lines[] = { {"La", "Ka", 33.30, 38.92}, {"La", "Kb", 37.8, 38.92},
                                   {"Br", "Ka", 11.91, 13.47}, {"Br", "Kb", 13.29, 13.47} };

    shared_ptr<const PeakFitDetPrefs> fit_prefs = meas->peakFitDetPrefs();
    if( !fit_prefs && meas->detector() )
      fit_prefs = meas->detector()->peakFitDetPrefs();
    const PeakFitUtils::CoarseResolutionType det_type = PeakFitUtils::effective_det_type( fit_prefs, spectrum, meas );

    vector<DetXray> xray_lines;
    switch( det_type )
    {
      case PeakFitUtils::CoarseResolutionType::High:
        result.detectorMaterialAssumed = "Ge (HPGe)";
        xray_lines.assign( std::begin(ge_lines), std::end(ge_lines) );
        break;

      case PeakFitUtils::CoarseResolutionType::Low:
        result.detectorMaterialAssumed = "NaI (iodine; a CsI detector would also have Cs K x-ray escapes)";
        xray_lines.assign( std::begin(i_lines), std::end(i_lines) );
        break;

      case PeakFitUtils::CoarseResolutionType::CZT:
        result.detectorMaterialAssumed = "CdTe/CdZnTe";
        xray_lines.assign( std::begin(czt_lines), std::end(czt_lines) );
        break;

      case PeakFitUtils::CoarseResolutionType::LaBr:
        result.detectorMaterialAssumed = "LaBr3";
        xray_lines.assign( std::begin(labr_lines), std::end(labr_lines) );
        break;

      case PeakFitUtils::CoarseResolutionType::MedRes:
      case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
      case PeakFitUtils::CoarseResolutionType::Unknown:
        result.detectorMaterialAssumed = isHPGe ? "Ge (HPGe)" : "unknown, so the common materials (Ge, NaI, CdTe/CZT, LaBr3) are all checked";
        xray_lines.assign( std::begin(ge_lines), std::end(ge_lines) );
        if( !isHPGe )
        {
          xray_lines.insert( end(xray_lines), std::begin(i_lines), std::end(i_lines) );
          xray_lines.insert( end(xray_lines), std::begin(czt_lines), std::end(czt_lines) );
          xray_lines.insert( end(xray_lines), std::begin(labr_lines), std::end(labr_lines) );
        }
        break;
    }//switch( det_type )

    // A generous upper bound on the escape-peak to parent-peak area ratio.  The escape is a surface
    //  effect: largest just above the K edge, then falling roughly as E^-2 as photons penetrate
    //  deeper.  NaI uses the GADRAS-derived iodine fractions (x2); other materials a 15% plateau to
    //  3x the K edge (Ge: ~10% near 20 keV, ~1% near 60 keV, measured), K-beta at 0.3 of K-alpha.
    const auto max_escape_fraction = []( const DetXray &x, const double parent_energy ) -> double {
      if( parent_energy <= x.k_edge )
        return 0.0;

      const bool kbeta = (string(x.line) == "Kb");
      if( string(x.element) == "I" )
      {
        double ka = 0.0, kb = 0.0;
        PeakFitUtils::nai_iodine_escape_fractions( std::min( parent_energy, 250.0 ), ka, kb );
        double frac = 2.0 * (kbeta ? kb : ka);
        if( parent_energy > 250.0 )
          frac *= (250.0*250.0) / (parent_energy*parent_energy);
        return frac;
      }//if( iodine )

      const double plateau_end = 3.0 * x.k_edge;
      double frac = 0.15;
      if( parent_energy > plateau_end )
        frac *= (plateau_end*plateau_end) / (parent_energy*parent_energy);
      return kbeta ? 0.3*frac : frac;
    };

    const shared_ptr<const PeakDef> queried_peak = find_peak_near_energy( input_energy, calc_window( input_energy ) );

    for( const DetXray &x : xray_lines )
    {
      XrayEscapeInfo info;
      info.element = x.element;
      info.line = x.line;
      info.xrayEnergy = x.energy;
      info.kEdgeEnergy = x.k_edge;
      info.parentEnergy = input_energy + x.energy;
      if( info.parentEnergy > spectrum_max_energy )
        continue;

      info.parentPeak = found_peak( info.parentEnergy, calc_window( info.parentEnergy ) );
      info.maxEscapeFraction = max_escape_fraction( x, info.parentEnergy );

      if( info.parentPeak && queried_peak && (info.parentPeak->peak->peakArea() > 0.0) )
        info.areaRatio = queried_peak->peakArea() / info.parentPeak->peak->peakArea();

      info.plausible = info.parentPeak.has_value() && info.areaRatio.has_value()
                       && (info.maxEscapeFraction > 0.0) && (*info.areaRatio <= info.maxEscapeFraction);

      result.xrayEscapes.push_back( std::move(info) );
    }//for( const DetXray &x : xray_lines )

    if( !queried_peak && !result.xrayEscapes.empty() )
      result.notes.push_back( "No peak was found at the queried energy, so whether it is an x-ray escape peak can not"
                              " be judged from areas; the x-ray escape entries only list where parents would be." );

    return result;
  }//EscapePeakCheckStatus escape_peak_check( const EscapePeakCheckOptions &options, InterSpec *interspec )


  SumPeakCheckStatus sum_peak_check( const SumPeakCheckOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "sum_peak_check: No InterSpec session available" );

    const double target_energy = options.energy;

    // Get the spectrum
    shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
    if( !meas )
    {
      throw runtime_error( "No measurement loaded for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );
    }

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram( options.specType );
    if( !spectrum )
    {
      throw runtime_error( "No spectrum displayed for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );
    }

    // Determine detector type for dead time threshold
    const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );
    const double dead_time_threshold = isHPGe ? 0.20 : 0.10;

    // Calculate dead time fraction
    double dead_time_fraction = 0.0;
    const bool has_valid_timing = (spectrum->live_time() > 0.0) && (spectrum->real_time() > 0.0);
    if( has_valid_timing )
      dead_time_fraction = (spectrum->real_time() - spectrum->live_time()) / spectrum->real_time();
    const double live_time = has_valid_timing ? spectrum->live_time() : -1.0;

    // Check if cascade-sum detection should be disabled based on distance
    // Cascade sums typically only occur at distances < 15 cm
    const double cascade_distance_threshold = 15.0 * PhysicalUnits::cm;
    const bool check_cascade_sums = !options.distance.has_value() || (*options.distance <= cascade_distance_threshold);

    // Initialize result
    SumPeakCheckStatus result;

    // Calculate search window for finding contributing peaks
    const std::function<float(double)> fwhm_fcn = make_expected_fwhm_fcn( interspec );
    const auto window_at = [&fwhm_fcn]( const double energy ) -> double {
      return std::max( 1.0, std::min( 0.5 * fwhm_fcn( std::max( energy, 1.0 ) ), 15.0 ) );
    };
    result.searchWindow = window_at( target_energy );

    // Get all peaks in the spectrum
    DetectedPeaksOptions peak_options;
    peak_options.specType = options.specType;
    peak_options.nonBackgroundPeaksOnly = false;
    const DetectedPeakStatus peak_status = detected_peaks( peak_options, interspec );
    const vector<shared_ptr<const PeakDef>> &all_peaks = peak_status.peaks;

    const auto nearest_peak = [&all_peaks]( const double energy, const double window ) -> shared_ptr<const PeakDef> {
      shared_ptr<const PeakDef> nearest;
      for( const shared_ptr<const PeakDef> &peak : all_peaks )
      {
        if( !peak || !peak->gausPeak() )
          continue;
        const double diff = fabs( peak->mean() - energy );
        if( (diff <= window) && (!nearest || (diff < fabs( nearest->mean() - energy ))) )
          nearest = peak;
      }
      return nearest;
    };

    const auto significance = []( const shared_ptr<const PeakDef> &p ) -> double {
      return (p->peakAreaUncert() > 0.0) ? (p->peakArea() / p->peakAreaUncert()) : 0.0;
    };

    // The peak at the queried energy, if there is one.
    const shared_ptr<const PeakDef> target_peak = nearest_peak( target_energy, result.searchWindow );

    // Decay cascades (with x-rays) of the nuclides assigned to peaks, built once each.
    struct NucCascades
    {
      std::vector<ceelo::DecayCascade> cascades;
      double per_parent_decay = 0.0;
    };
    std::map<const SandiaDecay::Nuclide *, NucCascades> cascade_cache;
    const auto cascades_of = [&cascade_cache]( const SandiaDecay::Nuclide *nuc ) -> const NucCascades & {
      auto pos = cascade_cache.find( nuc );
      if( pos == end(cascade_cache) )
      {
        NucCascades c;
        try
        {
          c.cascades = nuclide_cascades( nuc, c.per_parent_decay );
        }catch( std::exception & )
        {
          c.cascades.clear();
        }
        pos = cascade_cache.emplace( nuc, std::move(c) ).first;
      }
      return pos->second;
    };

    // True if the line of `peak` can be emitted twice in one decay (e.g., two transitions of the
    //  same energy in a cascade), so its self-sum is not necessarily random.
    const auto line_self_coincident = [&]( const shared_ptr<const PeakDef> &peak ) -> bool {
      const SandiaDecay::Nuclide * const nuc = peak->parentNuclide();
      if( !nuc )
        return false;
      const NucCascades &c = cascades_of( nuc );
      const double tol = std::max( 1.0, 0.5*peak->fwhm() );
      string type;
      return (cascade_pair_probability( c.cascades, c.per_parent_decay, peak->mean(), tol, peak->mean(), tol, type ) > 0.0);
    };

    // True if `nuc` (with its progeny, at its default age) emits a photon within `tol` of `energy`
    //  in at least `min_yield` of its decays.
    const auto nuclide_has_line_near = []( const SandiaDecay::Nuclide *nuc, const double energy,
                                           const double tol, const double min_yield = 1.0E-4 ) -> bool {
      if( !nuc )
        return false;
      const double activity = 1.0*SandiaDecay::Bq;
      SandiaDecay::NuclideMixture mix;
      mix.addAgedNuclideByActivity( nuc, activity, PeakDef::defaultDecayTime( nuc, nullptr ) / PhysicalUnits::second );
      for( const SandiaDecay::EnergyRatePair &erp : mix.photons( 0.0 ) )
      {
        if( ((erp.numPerSecond / activity) >= min_yield) && (fabs( erp.energy - energy ) <= tol) )
          return true;
      }
      return false;
    };

    // True if a peak has a source assigned with a real line at its energy (so it is not a sum peak).
    const auto is_assigned_real_line = [&]( const shared_ptr<const PeakDef> &peak ) -> bool {
      const SandiaDecay::Nuclide * const nuc = peak->parentNuclide();
      if( !nuc )
        return (peak->xrayElement() || peak->reaction());
      return nuclide_has_line_near( nuc, peak->mean(), std::max( 1.0, 0.5*peak->fwhm() ) );
    };

    // Pile-up time constant, tau: a line of rate R sums with itself at tau*R^2.  Measured from
    //  self-sum peaks of strong lines (other than the queried energy, which would be circular),
    //  else bounded from the dead time.
    bool tau_measured = false;
    if( has_valid_timing )
    {
      // The dead time per count bounds the pile-up resolving time from above
      const double total_counts = spectrum->gamma_count_sum();
      const double dead_time = spectrum->real_time() - spectrum->live_time();
      const double tau_upper_bound = ((total_counts > 0.0) && (dead_time > 0.0)) ? (dead_time / total_counts) : -1.0;

      vector<double> taus;
      string tau_sources;
      for( const shared_ptr<const PeakDef> &parent : all_peaks )
      {
        if( !parent || !parent->gausPeak() || (significance( parent ) < 10.0) )
          continue;

        const double sum_energy = 2.0 * parent->mean();
        const double window = window_at( sum_energy );
        if( fabs( sum_energy - target_energy ) <= std::max( window, result.searchWindow ) )
          continue;

        // A peak at 2E may instead be a real line - of the peak's assigned source, or of the line's
        //  own nuclide (e.g., Eu152's 688.7 keV next to 2 x 344.3 keV).
        const shared_ptr<const PeakDef> sum_peak = nearest_peak( sum_energy, window );
        if( !sum_peak || (significance( sum_peak ) < 3.0) || (sum_peak->peakArea() >= parent->peakArea())
            || is_assigned_real_line( sum_peak ) || line_self_coincident( parent )
            || nuclide_has_line_near( parent->parentNuclide(), sum_peak->mean(), window ) )
          continue;

        // ... and a "pile-up" time longer than the dead time per count means it is something else.
        const double tau = sum_peak->peakArea() * live_time / (parent->peakArea() * parent->peakArea());
        if( (tau_upper_bound > 0.0) && (tau > tau_upper_bound) )
          continue;

        taus.push_back( tau );

        char buffer[96];
        snprintf( buffer, sizeof(buffer), "%s%.2f keV = 2 x %.2f keV", (tau_sources.empty() ? "" : ", "),
                  sum_peak->mean(), parent->mean() );
        tau_sources += buffer;
      }//for( loop over potential self-summing peaks )

      if( !taus.empty() )
      {
        std::sort( begin(taus), end(taus) );
        const size_t n = taus.size();
        result.pileupTau = (n % 2) ? taus[n/2] : 0.5*(taus[n/2 - 1] + taus[n/2]);
        result.pileupTauIsUpperBound = false;
        result.pileupTauSource = "measured from the self-sum peak" + string(n > 1 ? "s " : " ") + tau_sources;
        tau_measured = true;
      }else
      {
        if( tau_upper_bound > 0.0 )
        {
          result.pileupTau = tau_upper_bound;
          result.pileupTauIsUpperBound = true;
          result.pileupTauSource = "an upper bound from the dead time per count (no self-sum peak of a strong"
                                   " line was found to measure it; pile-up is usually much rarer than this)";
        }
      }//if( measured tau ) / else
    }//if( has_valid_timing )

    // Structure to hold candidate sum peak pairs
    struct CandidatePair
    {
      shared_ptr<const PeakDef> peak1;
      shared_ptr<const PeakDef> peak2;
      double sum_energy = 0.0;
      double energy_error = 0.0;  // How far from target energy
      SumPeakType sum_type = SumPeakType::Unknown;
      double coincidence_fraction = 1.0;  // For cascade sums
      double selection_metric = 0.0;      // For ranking candidates
      std::string cascade_nuclide, coincidence_type, reason;
      bool cascade_also_possible = false;
      std::optional<double> expected_random_cps;
    };

    vector<CandidatePair> candidates;

    const std::optional<double> observed_cps = (target_peak && has_valid_timing)
                                ? std::optional<double>( target_peak->peakArea() / live_time ) : std::nullopt;
    const double observed_cps_uncert = (target_peak && has_valid_timing && (target_peak->peakAreaUncert() > 0.0))
                                ? (target_peak->peakAreaUncert() / live_time) : 0.0;

    char buffer[512];

    // Search for all pairs of peaks that sum to the target energy
    for( size_t i = 0; i < all_peaks.size(); ++i )
    {
      const shared_ptr<const PeakDef> &peak1 = all_peaks[i];
      if( !peak1 )
        continue;

      const double energy1 = peak1->mean();

      // Peak can sum with itself (for random summing)
      for( size_t j = i; j < all_peaks.size(); ++j )
      {
        const shared_ptr<const PeakDef> &peak2 = all_peaks[j];
        if( !peak2 )
          continue;

        const double energy2 = peak2->mean();
        const double sum_energy = energy1 + energy2;
        const double energy_error = fabs( sum_energy - target_energy );

        // Check if this pair sums to the target energy within the search window
        if( energy_error > result.searchWindow )
          continue;

        // A sum peak should always be smaller than the peaks creating it
        if( target_peak && ((target_peak->peakArea() >= peak1->peakArea())
                             || (target_peak->peakArea() >= peak2->peakArea())) )
          continue;

        CandidatePair candidate;
        candidate.peak1 = peak1;
        candidate.peak2 = peak2;
        candidate.sum_energy = sum_energy;
        candidate.energy_error = energy_error;

        const bool same_peak = (i == j);

        // Is this pair a cascade of a nuclide assigned to either peak?  (An x-ray peak may be
        //  assigned to the element, or to nothing.)
        if( check_cascade_sums )
        {
          for( const SandiaDecay::Nuclide *nuc : { peak1->parentNuclide(), peak2->parentNuclide() } )
          {
            if( !nuc || (nuc->symbol == candidate.cascade_nuclide) )
              continue;
            const NucCascades &c = cascades_of( nuc );
            const double tol1 = std::max( 1.0, 0.5*peak1->fwhm() ), tol2 = std::max( 1.0, 0.5*peak2->fwhm() );
            string type;
            const double prob = cascade_pair_probability( c.cascades, c.per_parent_decay, energy1, tol1, energy2, tol2, type );
            if( prob > 0.0 && ((candidate.cascade_nuclide.empty()) || (prob > candidate.coincidence_fraction)) )
            {
              candidate.cascade_nuclide = nuc->symbol;
              candidate.coincidence_fraction = prob;
              candidate.coincidence_type = type;
            }
          }//for( nuclides of the two peaks )
        }//if( check_cascade_sums )

        const bool is_cascade = !candidate.cascade_nuclide.empty();

        if( result.pileupTau && has_valid_timing )
        {
          const double r1 = peak1->peakArea() / live_time, r2 = peak2->peakArea() / live_time;
          candidate.expected_random_cps = same_peak ? ((*result.pileupTau) * r1 * r1)
                                                    : (2.0 * (*result.pileupTau) * r1 * r2);
        }

        // Why a pair is not a cascade, for the reason text
        string not_cascade_why;
        if( !check_cascade_sums )
          not_cascade_why = "cascade sums were not considered, as the source is more than 15 cm away";
        else if( !peak1->parentNuclide() && !peak2->parentNuclide() )
          not_cascade_why = "neither peak has a nuclide assigned, so cascades could not be checked (assign the"
                            " nuclide to a peak, e.g., with add_analysis_peaks_for_source, and check again)";
        else
          not_cascade_why = "the lines are not a known cascade of the assigned source";

        const bool can_compare = tau_measured && candidate.expected_random_cps && observed_cps;
        const double random_frac = can_compare ? (*candidate.expected_random_cps / *observed_cps) : 0.0;
        const bool significant_excess = can_compare
                                        && ((*observed_cps - *candidate.expected_random_cps) > 3.0*observed_cps_uncert);
        const bool legacy_random = has_valid_timing && (dead_time_fraction > dead_time_threshold);

        // When the observed and random-sum rates can not be compared, say why.
        string no_compare_why;
        if( !can_compare )
        {
          char why[256] = { '\0' };
          if( !observed_cps && candidate.expected_random_cps )
            snprintf( why, sizeof(why), "there is no fitted peak at this energy, so no observed rate to compare with"
                      " the %s%.2g cps that random summing would give", (tau_measured ? "" : "at most "),
                      *candidate.expected_random_cps );
          else if( observed_cps && candidate.expected_random_cps )
            snprintf( why, sizeof(why), "random summing would give at most %.2g cps of the %.2g cps observed (the"
                      " pile-up rate is only bounded, from the dead time)", *candidate.expected_random_cps, *observed_cps );
          else if( !observed_cps )
            snprintf( why, sizeof(why), "there is no fitted peak at this energy, and the pile-up rate could not be found" );
          else
            snprintf( why, sizeof(why), "the pile-up rate could not be found" );
          no_compare_why = why;
        }//if( !can_compare )

        // Is the queried energy a real line - of the peak's assigned source, or of a nuclide assigned
        //  to one of the pair (e.g., Eu152's 688.67 keV at 2 x 344.28, or a crossover transition)?
        //  Then it need not be a sum peak.  (At least 1E-3 per decay, so a dense chain's trace lines
        //  do not count.)
        string real_line_of;
        const double line_energy = target_peak ? target_peak->mean() : target_energy;
        const double line_tol = target_peak ? std::max( 1.0, 0.5*target_peak->fwhm() ) : result.searchWindow;
        const double min_line_yield = 1.0E-3;
        if( target_peak && is_assigned_real_line( target_peak ) )
          real_line_of = target_peak->sourceName();
        else if( nuclide_has_line_near( peak1->parentNuclide(), line_energy, line_tol, min_line_yield ) )
          real_line_of = peak1->parentNuclide()->symbol;
        else if( nuclide_has_line_near( peak2->parentNuclide(), line_energy, line_tol, min_line_yield ) )
          real_line_of = peak2->parentNuclide()->symbol;

        // With only an upper bound on the pile-up rate, random summing is ruled out if even that
        //  bound is far below the observed rate.
        const bool random_ruled_out = !tau_measured && observed_cps && candidate.expected_random_cps
                                      && (*observed_cps > 3.0*(*candidate.expected_random_cps));

        if( is_cascade )
        {
          if( can_compare && (random_frac >= 0.7) )
          {
            candidate.sum_type = SumPeakType::RandomSum;
            candidate.cascade_also_possible = true;
            snprintf( buffer, sizeof(buffer), "the two lines are emitted together in %.2g of %s decays, but random"
                      " summing alone accounts for %.0f%% of the observed rate", candidate.coincidence_fraction,
                      candidate.cascade_nuclide.c_str(), 100.0*random_frac );
          }else if( can_compare && (random_frac >= 0.15) && significant_excess )
          {
            candidate.sum_type = SumPeakType::CascadeAndRandomSum;
            snprintf( buffer, sizeof(buffer), "the two lines are emitted together (%s, in %.2g of %s decays), and"
                      " random summing accounts for about %.0f%% of the observed rate", candidate.coincidence_type.c_str(),
                      candidate.coincidence_fraction, candidate.cascade_nuclide.c_str(), 100.0*random_frac );
          }else
          {
            candidate.sum_type = SumPeakType::CascadeSum;
            char random_part[128] = { '\0' };
            if( can_compare )
              snprintf( random_part, sizeof(random_part), ", and random summing accounts for about %.0f%% of the"
                        " observed rate%s", 100.0*random_frac,
                        (random_frac >= 0.15) ? " (not a significant part, given the uncertainty)" : "" );
            snprintf( buffer, sizeof(buffer), "the two lines are emitted together (%s, in %.2g of %s decays)%s%s%s",
                      candidate.coincidence_type.c_str(), candidate.coincidence_fraction, candidate.cascade_nuclide.c_str(),
                      random_part, (can_compare ? "" : "; "), no_compare_why.c_str() );
          }
        }else if( !real_line_of.empty() )
        {
          candidate.sum_type = SumPeakType::Unknown;
          snprintf( buffer, sizeof(buffer), "the energies add up, but %s has a line at this energy, so a peak here"
                    " need not be a sum peak (summing may add to it); %s", real_line_of.c_str(),
                    (can_compare ? not_cascade_why.c_str() : no_compare_why.c_str()) );
        }else if( same_peak && !line_self_coincident( peak1 ) )
        {
          if( random_ruled_out )
          {
            candidate.sum_type = SumPeakType::Unknown;
            snprintf( buffer, sizeof(buffer), "a line can only sum with itself randomly, but even the upper bound on"
                      " random summing (%.2g cps) is far below the %.2g cps observed, so this is likely a real line",
                      *candidate.expected_random_cps, *observed_cps );
          }else if( can_compare && (random_frac < 0.2) && significant_excess )
          {
            candidate.sum_type = SumPeakType::Unknown;
            snprintf( buffer, sizeof(buffer), "a line can only sum with itself randomly, but the observed rate is %.2g"
                      " times what random summing gives, so this may be a real line", 1.0/std::max(random_frac, 1.0E-9) );
          }else
          {
            candidate.sum_type = SumPeakType::RandomSum;
            snprintf( buffer, sizeof(buffer), "a decay emits this line at most once, so it can only sum with itself"
                      " by random coincidence (pile-up)%s%s",
                      can_compare ? ", and the observed rate is consistent with that" : "; ", no_compare_why.c_str() );
          }
        }else if( can_compare )
        {
          if( random_frac >= (1.0/3.0) )
          {
            candidate.sum_type = SumPeakType::RandomSum;
            snprintf( buffer, sizeof(buffer), "random summing of these lines accounts for %.0f%% of the observed rate;"
                      " %s", 100.0*random_frac, not_cascade_why.c_str() );
          }else
          {
            candidate.sum_type = SumPeakType::Unknown;
            snprintf( buffer, sizeof(buffer), "random summing accounts for only %.0f%% of the observed rate, and %s",
                      100.0*random_frac, not_cascade_why.c_str() );
          }
        }else if( legacy_random )
        {
          candidate.sum_type = SumPeakType::RandomSum;
          snprintf( buffer, sizeof(buffer), "the dead time (%.0f%%) is high enough for random summing; %s",
                    100.0*dead_time_fraction, no_compare_why.c_str() );
        }else
        {
          candidate.sum_type = SumPeakType::Unknown;
          snprintf( buffer, sizeof(buffer), "the energies add up, but %s; %s", not_cascade_why.c_str(),
                    no_compare_why.c_str() );
        }//if( is_cascade ) / else ...
        candidate.reason = buffer;
        if( is_cascade && !real_line_of.empty() )
          candidate.reason += "; " + real_line_of + " also has a line at this energy (e.g., a crossover"
                              " transition), which may be most of a peak here";

        // For cascade sums: area1 * area2 * coincidence_fraction; for others: area1 * area2
        candidate.selection_metric = peak1->peakArea() * peak2->peakArea() * (is_cascade ? candidate.coincidence_fraction : 1.0);

        candidates.push_back( candidate );
      }//for( j : all_peaks )
    }//for( i : all_peaks )

    // If we found candidates, select the best one
    if( !candidates.empty() )
    {
      // First, sort by energy error (closest match first)
      std::stable_sort( candidates.begin(), candidates.end(),
                []( const CandidatePair &lhs, const CandidatePair &rhs ) {
                  return lhs.energy_error < rhs.energy_error;
                } );

      // If we have multiple candidates with similar energies (within 0.5 keV),
      // select the one with the largest selection metric
      const double best_energy_error = candidates[0].energy_error;
      const double energy_similarity_threshold = 0.5; // keV

      const CandidatePair *best_candidate = &candidates[0];
      for( const CandidatePair &candidate : candidates )
      {
        if( fabs(candidate.energy_error - best_energy_error) > energy_similarity_threshold )
          break; // Candidates are sorted by energy error
        if( candidate.selection_metric > best_candidate->selection_metric )
          best_candidate = &candidate;
      }

      CandidatePair best = *best_candidate;

      // Without a measured pile-up rate, a dead-time-based random-sum call is checked against the
      //  largest peak: if it does not visibly sum with itself, random summing is unlikely.
      if( (best.sum_type == SumPeakType::RandomSum) && !tau_measured && (best.peak1 != best.peak2) )
      {
        double max_area = 0.0;
        double max_area_energy = 0.0;
        for( const shared_ptr<const PeakDef> &peak : all_peaks )
        {
          if( peak && (peak->peakArea() > max_area) )
          {
            max_area = peak->peakArea();
            max_area_energy = peak->mean();
          }
        }

        const double expected_self_sum_energy = 2.0 * max_area_energy;
        const bool self_sum_in_range = (expected_self_sum_energy <= spectrum->gamma_energy_max());
        const bool checking_self_sum = (fabs(target_energy - expected_self_sum_energy) <= result.searchWindow);
        if( self_sum_in_range && !checking_self_sum && !nearest_peak( expected_self_sum_energy, result.searchWindow ) )
        {
          best.sum_type = SumPeakType::Unknown;
          snprintf( buffer, sizeof(buffer), "although the dead time is high, the largest peak (%.2f keV) does not"
                    " visibly sum with itself, so random summing is unlikely", max_area_energy );
          best.reason = buffer;
        }
      }//if( dead-time based random-sum call )

      // Create the result
      SumPeakInfo sum_info;
      sum_info.sumType = best.sum_type;
      sum_info.firstPeak = best.peak1;
      sum_info.secondPeak = best.peak2;
      sum_info.reason = best.reason;
      sum_info.cascadeAlsoPossible = best.cascade_also_possible;
      sum_info.observedSumCps = observed_cps;
      sum_info.expectedRandomSumCps = best.expected_random_cps;

      if( !best.cascade_nuclide.empty() )
      {
        sum_info.coincidenceFraction = best.coincidence_fraction;
        sum_info.cascadeNuclide = best.cascade_nuclide;
        sum_info.coincidenceType = best.coincidence_type;
      }

      // Generate labels
      const double e1 = best.peak1->mean();
      const double e2 = best.peak2->mean();

      const SandiaDecay::Nuclide * const nuc1 = best.peak1->parentNuclide();
      const SandiaDecay::Nuclide * const nuc2 = best.peak2->parentNuclide();
      const ReactionGamma::Reaction * const rxn1 = best.peak1->reaction();
      const ReactionGamma::Reaction * const rxn2 = best.peak2->reaction();
      const SandiaDecay::Element * const xray1 = best.peak1->xrayElement();
      const SandiaDecay::Element * const xray2 = best.peak2->xrayElement();

      if( (best.sum_type == SumPeakType::CascadeSum) || (best.sum_type == SumPeakType::CascadeAndRandomSum) )
      {
        snprintf( buffer, sizeof(buffer), "%s cascade-sum %.2f + %.2f keV", best.cascade_nuclide.c_str(), e1, e2 );
        sum_info.userLabel = buffer;
      }
      else if( (nuc1 || rxn1 || xray1) && (nuc2 || rxn2 || xray2) )
      {
        // Both peaks have sources - include in label
        string src1, src2;
        if( nuc1 ) src1 = nuc1->symbol;
        else if( rxn1 ) src1 = rxn1->name();
        else if( xray1 ) src1 = xray1->symbol;

        if( nuc2 ) src2 = nuc2->symbol;
        else if( rxn2 ) src2 = rxn2->name();
        else if( xray2 ) src2 = xray2->symbol;

        const char *sum_type_str = (best.sum_type == SumPeakType::RandomSum) ? "random-sum" : "sum";
        snprintf( buffer, sizeof(buffer), "%s + %s %s %.2f + %.2f keV",
                 src1.c_str(), src2.c_str(), sum_type_str, e1, e2 );
        sum_info.userLabel = buffer;
      }
      else
      {
        // Create user label
        snprintf( buffer, sizeof(buffer), "Sum %.2f + %.2f keV", e1, e2 );
        sum_info.userLabel = buffer;
      }

      result.sumPeakInfo = sum_info;
    }//if( we have candidates )

    // Other sources with a line here, which could explain the peak instead
    try
    {
      std::set<string> names;
      for( const auto &item : get_nuclides_with_characteristics_in_energy_range( target_energy - result.searchWindow,
                                                                                target_energy + result.searchWindow,
                                                                                interspec ) )
      {
        if( std::holds_alternative<const SandiaDecay::Nuclide *>(item) )
          names.insert( std::get<const SandiaDecay::Nuclide *>(item)->symbol );
        else if( std::holds_alternative<const SandiaDecay::Element *>(item) )
          names.insert( std::get<const SandiaDecay::Element *>(item)->symbol + " x-ray" );
        else if( std::holds_alternative<const ReactionGamma::Reaction *>(item) )
          names.insert( std::get<const ReactionGamma::Reaction *>(item)->name() );
      }
      result.otherSourcesAtThisEnergy.assign( begin(names), end(names) );
    }catch( std::exception & )
    {
    }

    return result;
  }//SumPeakCheckStatus sum_peak_check( const SumPeakCheckOptions &options, InterSpec *interspec )


  const char* to_string( ComptonScatterConfidence c )
  {
    switch( c )
    {
      case ComptonScatterConfidence::No:         return "No";
      case ComptonScatterConfidence::Unlikely:   return "Unlikely";
      case ComptonScatterConfidence::Likely:     return "Likely";
      case ComptonScatterConfidence::VeryLikely: return "VeryLikely";
    }
    return "Unknown";
  }//const char* to_string( ComptonScatterConfidence c )


  const char* to_string( ComptonScatterType t )
  {
    switch( t )
    {
      case ComptonScatterType::None:               return "None";
      case ComptonScatterType::ComptonPeak:        return "ComptonPeak";
      case ComptonScatterType::ComptonEdge:        return "ComptonEdge";
      case ComptonScatterType::ComptonPeakOrEdge: return "ComptonPeakOrEdge";
    }
    return "Unknown";
  }//const char* to_string( ComptonScatterType t )


  ComptonPeakCheckStatus compton_peak_check(
      const shared_ptr<const PeakDef> &candidate,
      const shared_ptr<const deque<shared_ptr<const PeakDef>>> &user_peaks,
      const shared_ptr<const SpecUtils::Measurement> &spectrum,
      const shared_ptr<const DetectorPeakResponse> &detector,
      const std::function<float(double)> &expected_fwhm )
  {
    ComptonPeakCheckStatus result;

    if( !candidate || !candidate->gausPeak() || !user_peaks || user_peaks->empty() )
    {
      result.reason = "no Gaussian candidate peak, or no other peaks";
      return result;
    }

    const double electron_rest_mass = 510.9989; // keV; matches escape_peak_check at AnalystChecks.cpp:2025

    result.candidateMean = candidate->mean();
    result.candidateFwhm = candidate->fwhm();
    const double cand_amp = candidate->amplitude();

    // Cheap-path FWHM only: use DRF resolution info if available; otherwise
    // -1 sentinel triggers the sqrt(E)-vs-parent fallback below.
    double expected_det_fwhm = -1.0;
    if( expected_fwhm )
    {
      expected_det_fwhm = expected_fwhm( candidate->mean() );
      if( !(expected_det_fwhm > 0.0) )
        expected_det_fwhm = -1.0;
    }else if( detector && detector->hasResolutionInfo() )
    {
      try
      {
        expected_det_fwhm = detector->peakResolutionFWHM( static_cast<float>(candidate->mean()) );
      }catch( std::exception & )
      {
        expected_det_fwhm = -1.0;
      }
    }
    result.expectedDetectorFwhm = expected_det_fwhm;

    // Search higher-energy user peaks for kinematic match.
    struct ParentMatch
    {
      shared_ptr<const PeakDef> parent;
      ComptonScatterType type = ComptonScatterType::None;
      double scatter_angle_deg = 0.0;
      double energy_match_score = std::numeric_limits<double>::max(); // |dE| / candidateFwhm
    };

    // Hoist candidate-side constants out of the parent loop and guard against
    // a degenerate zero-mean candidate (would divide by zero in the kinematic
    // inversion below).
    const double Eprime = candidate->mean();
    const double cand_fwhm = std::max( candidate->fwhm(), 1e-6 );
    if( Eprime <= 0.0 )
      return result;

    // Every higher-energy peak whose scatter features match the candidate's energy
    vector<ParentMatch> matches;

    for( const shared_ptr<const PeakDef> &p : *user_peaks )
    {
      if( !p || !p->gausPeak() || (p == candidate) )
        continue;
      if( p->mean() <= Eprime + 0.5*cand_fwhm )
        continue;

      const double Ep = p->mean();
      const double alpha = Ep / electron_rest_mass;

      // Energy-distance window for accepting a candidate as a kinematic match.
      // Detector resolution at the feature energy is roughly the parent's FWHM
      // (HPGe ~constant, NaI ~sqrt(E)) and the candidate's own FWHM bounds the
      // centroid uncertainty for broad backscatter / edge fits, so take the
      // looser of the two scaled appropriately.
      const double match_window = std::max( 2.5 * p->fwhm(), 1.25 * cand_fwhm );

      // Kinematic ComptonPeak inversion:
      //   E' = Ep / ( 1 + alpha (1 - cos theta) )
      //     => cos theta = 1 - (1/alpha)*( Ep/E' - 1 )
      const double cos_theta_raw = 1.0 - (1.0/alpha) * ( Ep / Eprime - 1.0 );

      double dE_peak_keV = std::numeric_limits<double>::max();
      double theta_deg = 0.0;
      if( cos_theta_raw >= -1.0 && cos_theta_raw <= 1.0 )
      {
        // Geometrically possible scatter angle.
        const double Eprime_at_theta = Ep / ( 1.0 + alpha*(1.0 - cos_theta_raw) );
        dE_peak_keV = std::abs( Eprime - Eprime_at_theta );
        theta_deg = std::acos( cos_theta_raw ) * 180.0 / 3.14159265358979323846;
      }else if( cos_theta_raw < -1.0 && Ep > 100.0 )
      {
        // Candidate is below the kinematic minimum (E_bs at theta=180 deg).
        // Multiple-scatter buildup and broad-peak fit effects can pull the
        // centroid below E_bs; treat as a backscatter (theta=180) match.
        // Gate on Ep > 100 keV: below ~100 keV the kinematic separation
        // E_p - E_bs collapses (e.g., 18 keV at Ep=80, 12 keV at Ep=60), so
        // backscatter and parent are no longer distinguishable features and
        // any candidate near E_bs is more plausibly a different real peak.
        const double E_bs = Ep / (1.0 + 2.0*alpha);
        dE_peak_keV = std::abs( Eprime - E_bs );
        theta_deg = 180.0;
      }// (cos_theta_raw > 1 is impossible given the parent-above-candidate filter.)

      const double dpeak = dE_peak_keV / cand_fwhm;  // bucket_match scoring stays in candidate-FWHM units

      // ComptonEdge: E_edge = Ep * 2 alpha / ( 1 + 2 alpha )
      const double E_edge = Ep * 2.0 * alpha / ( 1.0 + 2.0 * alpha );
      const double dE_edge_keV = std::abs( Eprime - E_edge );
      const double dedge = dE_edge_keV / cand_fwhm;

      ComptonScatterType this_type = ComptonScatterType::None;
      double this_score = std::numeric_limits<double>::max();
      double this_angle = 0.0;
      const bool peak_in_window = (dE_peak_keV <= match_window);
      const bool edge_in_window = (dE_edge_keV <= match_window);
      if( peak_in_window && edge_in_window
          && (Ep >= 200.0) && (Ep <= 320.0) )
      {
        this_type = ComptonScatterType::ComptonPeakOrEdge;
        this_score = std::min( dpeak, dedge );
        this_angle = theta_deg;  // 180 deg if below-kinematic branch matched
      }else if( peak_in_window && (!edge_in_window || dpeak <= dedge) )
      {
        this_type = ComptonScatterType::ComptonPeak;
        this_score = dpeak;
        this_angle = theta_deg;
      }else if( edge_in_window )
      {
        this_type = ComptonScatterType::ComptonEdge;
        this_score = dedge;
        this_angle = 180.0;
      }else
      {
        continue;
      }

      ParentMatch match;
      match.parent = p;
      match.type = this_type;
      match.scatter_angle_deg = this_angle;
      match.energy_match_score = this_score;
      matches.push_back( match );
    }//for( each potential parent )

    if( matches.empty() )
    {
      result.reason = "no higher-energy peak has a backscatter peak or Compton edge at the candidate's energy";
      return result;
    }

    // The best parent that could have produced the candidate - a parent must be larger than its
    //  scatter feature, which takes at most half its area.  The energy match alone does not rank
    //  them: inside a parent's backscatter range the kinematics are solved for the candidate's own
    //  energy, so every such parent matches exactly.  Prefer parents whose match forms a visible
    //  feature - a backscatter peak (scatters of ~120-180 degrees bunch up in energy) or an edge -
    //  then larger parents, then the closer energy match.
    const auto forms_feature = []( const ParentMatch &m ) -> bool {
      return (m.type == ComptonScatterType::ComptonEdge) || (m.type == ComptonScatterType::ComptonPeakOrEdge)
             || (m.scatter_angle_deg >= 120.0);
    };
    std::stable_sort( begin(matches), end(matches), [&forms_feature]( const ParentMatch &a, const ParentMatch &b ){
      if( forms_feature(a) != forms_feature(b) )
        return forms_feature(a);
      if( a.parent->amplitude() != b.parent->amplitude() )
        return a.parent->amplitude() > b.parent->amplitude();
      return a.energy_match_score < b.energy_match_score;
    } );

    char reason_buffer[256] = { '\0' };
    const ParentMatch *best_ptr = nullptr;
    for( const ParentMatch &match : matches )
    {
      const double amp = match.parent->amplitude();
      if( (amp > cand_amp) && ((cand_amp / amp) <= 0.5) )
      {
        best_ptr = &match;
        break;
      }
    }//for( loop over matches, closest first )

    if( !best_ptr )
    {
      const shared_ptr<const PeakDef> &closest = matches.front().parent;
      const double closest_amp = closest->amplitude();
      if( closest_amp <= cand_amp )
        snprintf( reason_buffer, sizeof(reason_buffer), "the candidate peak is larger than its kinematic parent"
                  " at %.2f keV%s, so it is not a scatter feature", closest->mean(),
                  (matches.size() > 1) ? " (and is too large for every other matching parent)" : "" );
      else
        snprintf( reason_buffer, sizeof(reason_buffer), "the candidate peak's area is %.2g of its kinematic parent's"
                  " (at %.2f keV)%s - too large for a scatter feature", cand_amp / closest_amp, closest->mean(),
                  (matches.size() > 1) ? ", and too large for every other matching parent" : "" );
      result.reason = reason_buffer;
      return result;
    }//if( no parent could have produced the candidate )

    const ParentMatch best = *best_ptr;
    const shared_ptr<const PeakDef> parent = best.parent;
    const double parent_amp = parent->amplitude();
    const double area_fraction = (parent_amp > 0.0) ? (cand_amp / parent_amp) : 0.0;

    // FWHM check.  Two parallel metrics:
    //   - ratio (cand/expected): captures resolution-relative broadening.
    //     Informative for HPGe where physics broadening (a few keV from
    //     angular spread, multiple-scatter buildup) dominates the small
    //     detector resolution and yields large ratios.
    //   - absolute excess (cand - expected) in keV: captures physically
    //     meaningful absolute broadening.  At low-resolution detectors (NaI)
    //     the detector FWHM dominates the total in quadrature so the same
    //     physics-driven broadening yields only a modest ratio (~1.1-1.4),
    //     but the absolute keV excess is still meaningful.
    // We bucket both, then take the more confident verdict.
    //
    // Detector FWHM at the candidate energy: prefer DRF cheap-path; otherwise
    // sqrt(E)-scale parent FWHM (HPGe ~ sqrt(a+bE), NaI ~ sqrt(E)).
    double expected_at_cand = expected_det_fwhm;
    if( expected_at_cand <= 0.0 )
    {
      const double pmean = parent->mean();
      const double pfwhm = parent->fwhm();
      expected_at_cand = (pmean > 0.0)
                          ? pfwhm * std::sqrt( candidate->mean() / pmean )
                          : pfwhm;
    }

    const double fwhm_ratio = (expected_at_cand > 0.0)
                              ? (candidate->fwhm() / expected_at_cand)
                              : 0.0;
    const double fwhm_excess_keV = std::max( 0.0, candidate->fwhm() - expected_at_cand );

    // Hard guard: a real photopeak (ratio close to 1, near-zero absolute excess)
    // is rejected.  Require failure on BOTH axes to bail out.
    if( fwhm_ratio < 1.15 && fwhm_excess_keV < 0.5 )
    {
      snprintf( reason_buffer, sizeof(reason_buffer), "the candidate peak is as narrow as a photopeak (FWHM %.3g keV,"
                " expected %.3g keV), while scatter features are broader", candidate->fwhm(), expected_at_cand );
      result.reason = reason_buffer;
      return result;
    }

    // Parent S/C ratio over a FWHM-scaled window centred on the parent.
    // Using the full continuum ROI would include continuum far from the peak
    // (especially for wide multi-peak ROIs), giving a misleadingly low ratio;
    // a +/- FWHM window captures the peak's main body and the continuum
    // immediately under it.  We also collect all peaks sharing the continuum
    // so CDF-step continuum integrations get the correct sibling contributions.
    double parent_sc = 0.0;
    assert( parent->continuum() );
    if( parent->continuum() )
    {
      const shared_ptr<const PeakContinuum> pcont = parent->continuum();

      const double parent_mean = parent->mean();
      const double parent_fwhm = parent->fwhm();
      const double sc_lower = parent_mean - parent_fwhm;
      const double sc_upper = parent_mean + parent_fwhm;
      
      // The peak list only affects offset_integral for step continua.  For
      // polynomial / linear / etc. types the peak list is ignored by the
      // integration routine, so skip the sibling-peak scan in that case.
      std::vector<const PeakDef *> roi_peaks;
      if( PeakContinuum::is_step_continuum( pcont->type() ) )
      {
        roi_peaks.reserve( 4 );
        for( const shared_ptr<const PeakDef> &q : *user_peaks )
        {
          if( q && (q->continuum() == pcont) )
            roi_peaks.push_back( q.get() );
        }
      }
      if( roi_peaks.empty() )
        roi_peaks.push_back( parent.get() );
      
      try
      {
        const double cont_integral = pcont->offset_integral( sc_lower, sc_upper, spectrum, roi_peaks.data(), roi_peaks.size() );
        if( cont_integral > 0.0 )
          parent_sc = parent_amp / cont_integral;
      }catch( std::exception & )
      {
        parent_sc = 0.0;
      }
    }
    result.parentSignalToContinuum = parent_sc;

    // Bucket each criterion conservatively; final confidence = min over criteria.
    // FWHM bucket: take the more confident of two parallel verdicts (ratio
    // and absolute keV excess).  See FWHM-check comment above for rationale.
    auto bucket_fwhm = []( double r, double excess_keV ) {
      ComptonScatterConfidence ratio_b = ComptonScatterConfidence::No;
      if(      r >= 2.00 ) ratio_b = ComptonScatterConfidence::VeryLikely;
      else if( r >= 1.40 ) ratio_b = ComptonScatterConfidence::Likely;
      else if( r >= 1.15 ) ratio_b = ComptonScatterConfidence::Unlikely;

      ComptonScatterConfidence excess_b = ComptonScatterConfidence::No;
      if(      excess_keV >= 5.0 ) excess_b = ComptonScatterConfidence::VeryLikely;
      else if( excess_keV >= 2.0 ) excess_b = ComptonScatterConfidence::Likely;
      else if( excess_keV >= 0.5 ) excess_b = ComptonScatterConfidence::Unlikely;

      return (static_cast<int>(ratio_b) >= static_cast<int>(excess_b)) ? ratio_b : excess_b;
    };
    auto bucket_area = []( double f ) {
      if( f >= 0.02 && f <= 0.20 ) return ComptonScatterConfidence::VeryLikely;
      if( f >= 0.005 && f <= 0.30 ) return ComptonScatterConfidence::Likely;
      if( f < 0.005 || (f > 0.30 && f <= 0.50) ) return ComptonScatterConfidence::Unlikely;
      return ComptonScatterConfidence::No;
    };
    auto bucket_sc = []( double s ) {
      if( s >= 3.0 ) return ComptonScatterConfidence::VeryLikely;
      if( s >= 1.0 ) return ComptonScatterConfidence::Likely;
      if( s >= 0.3 ) return ComptonScatterConfidence::Unlikely;
      return ComptonScatterConfidence::No;
    };
    auto bucket_angle = [&]( double theta_deg ) {
      if( best.type == ComptonScatterType::ComptonEdge )
        return ComptonScatterConfidence::VeryLikely; // angle implicit in edge match
      // Thresholds loosened to accommodate angle uncertainty from ~0.5 FWHM
      // drift in the candidate mean (continuum-shape effects on broad scatter
      // peaks).  In the typical scatter kinematic regime dtheta/dE' is several
      // deg/keV, so a half-FWHM drift in fitted mean translates to tens of
      // deg of angle uncertainty -- tight angle bins would be too fragile.
      if( theta_deg >= 120.0 ) return ComptonScatterConfidence::VeryLikely;
      if( theta_deg >=  75.0 ) return ComptonScatterConfidence::Likely;
      if( theta_deg >=  30.0 ) return ComptonScatterConfidence::Unlikely;
      return ComptonScatterConfidence::No;
    };
    auto bucket_match = []( double s ) {
      if( s <= 0.5 ) return ComptonScatterConfidence::VeryLikely;
      if( s <= 1.0 ) return ComptonScatterConfidence::Likely;
      if( s <= 2.0 ) return ComptonScatterConfidence::Unlikely;
      return ComptonScatterConfidence::No;
    };

    const ComptonScatterConfidence fwhm_b = bucket_fwhm( fwhm_ratio, fwhm_excess_keV );

    const ComptonScatterConfidence buckets[] = {
      fwhm_b,
      bucket_area( area_fraction ),
      bucket_sc( parent_sc ),
      bucket_angle( best.scatter_angle_deg ),
      bucket_match( best.energy_match_score ),
    };

    // What each criterion measured, in `buckets` order, to say which ones held the verdict down.
    char criteria_buffer[5][96] = {};
    snprintf( criteria_buffer[0], sizeof(criteria_buffer[0]), "width %.2gx the detector FWHM (%.2g keV wider)",
              fwhm_ratio, fwhm_excess_keV );
    snprintf( criteria_buffer[1], sizeof(criteria_buffer[1]), "area %.2g of the parent's (0.02-0.2 is typical)",
              area_fraction );
    snprintf( criteria_buffer[2], sizeof(criteria_buffer[2]), "parent signal-to-continuum %.2g", parent_sc );
    if( best.type == ComptonScatterType::ComptonEdge )
      snprintf( criteria_buffer[3], sizeof(criteria_buffer[3]), "Compton edge match" );
    else
      snprintf( criteria_buffer[3], sizeof(criteria_buffer[3]), "scatter angle %.0f degrees", best.scatter_angle_deg );
    snprintf( criteria_buffer[4], sizeof(criteria_buffer[4]), "energy offset %.2g FWHM from the expected feature",
              best.energy_match_score );

    ComptonScatterConfidence overall = ComptonScatterConfidence::VeryLikely;
    for( const ComptonScatterConfidence b : buckets )
    {
      if( static_cast<int>(b) < static_cast<int>(overall) )
        overall = b;
    }

    // FWHM floor-raiser: a candidate that is much wider than detector
    // resolution is hard to explain as a real photopeak, so strong FWHM
    // evidence lifts the overall verdict by one tier.  We still respect a
    // hard No on any other axis -- a No anywhere keeps the overall as No.
    if( overall != ComptonScatterConfidence::No )
    {
      ComptonScatterConfidence floor = ComptonScatterConfidence::No;
      if( fwhm_b == ComptonScatterConfidence::VeryLikely )
        floor = ComptonScatterConfidence::Likely;
      else if( fwhm_b == ComptonScatterConfidence::Likely )
        floor = ComptonScatterConfidence::Unlikely;
      if( static_cast<int>(floor) > static_cast<int>(overall) )
        overall = floor;
    }

    result.confidence = overall;
    result.type = best.type;

    if( overall == ComptonScatterConfidence::No )
    {
      snprintf( reason_buffer, sizeof(reason_buffer), "the candidate peak matches the kinematics of the %.2f keV peak,"
                " but its width, area, the parent's prominence, or the scatter angle make that unlikely", parent->mean() );
      result.reason = reason_buffer;
    }else
    {
      snprintf( reason_buffer, sizeof(reason_buffer), "the candidate peak's energy, width, and area (%.2g of the"
                " parent's) fit a %s of the %.2f keV peak", area_fraction,
                (best.type == ComptonScatterType::ComptonEdge) ? "Compton edge"
                  : ((best.type == ComptonScatterType::ComptonPeakOrEdge) ? "backscatter peak or Compton edge"
                                                                         : "backscatter peak"),
                parent->mean() );
      result.reason = reason_buffer;

      // The verdict is the least confident criterion; name those below VeryLikely.
      string limits;
      for( size_t i = 0; i < 5; ++i )
      {
        if( buckets[i] != ComptonScatterConfidence::VeryLikely )
          limits += (limits.empty() ? "" : ", ") + string( criteria_buffer[i] );
      }
      if( !limits.empty() )
        result.reason += "; " + string( to_string( overall ) ) + " rather than VeryLikely because of: " + limits;
    }

    char label_buf[180] = { '\0' };
    if( (overall != ComptonScatterConfidence::No) && (best.type != ComptonScatterType::None) )
    {
      string type_str;
      switch( best.type )
      {
        case ComptonScatterType::None:
          break;

        case ComptonScatterType::ComptonPeak:
          type_str = "Compton backscatter peak";
          break;
        case ComptonScatterType::ComptonEdge:
          type_str = "Compton edge";
          break;
        case ComptonScatterType::ComptonPeakOrEdge:
          type_str = "Compton backscatter or edge";
          break;
      }//switch( best.type )

      if( parent->parentNuclide() )
      {
        snprintf( label_buf, sizeof(label_buf), "%s from %s %.2f keV peak",
                 type_str.c_str(), parent->parentNuclide()->symbol.c_str(), parent->mean() );
      }else if( parent->reaction() )
      {
        snprintf( label_buf, sizeof(label_buf), "%s from %s %.2f keV peak",
                 type_str.c_str(), parent->reaction()->name().c_str(), parent->mean() );
      }else
      {
        snprintf( label_buf, sizeof(label_buf), "%s from %.2f keV peak", type_str.c_str(), parent->mean() );
      }
    }//if( (overall != ComptonScatterConfidence::No) && (best.type != ComptonScatterType::None) )

    ComptonScatterParentInfo pinfo;
    pinfo.parentPeak = parent;
    pinfo.parentEnergy = parent->mean();
    pinfo.parentFwhm = parent->fwhm();
    pinfo.scatterAngle = best.scatter_angle_deg;
    pinfo.areaFraction = area_fraction;
    pinfo.userLabel = string( label_buf );
    result.parentInfo = pinfo;

    return result;
  }//ComptonPeakCheckStatus compton_peak_check( candidate, user_peaks, spectrum, detector )


  std::vector<ComptonKinematicParent> compton_kinematic_parents( const double energy,
                                    const double window,
                                    const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                                    const std::shared_ptr<const SpecUtils::Measurement> &spectrum )
  {
    const double electron_rest_mass = 510.9989;
    const double min_parent_energy = 50.0, min_parent_nsigma = 5.0;

    vector<ComptonKinematicParent> answer;
    if( !(energy > 0.0) || !(window > 0.0) )
      return answer;

    const double live_time = (spectrum && (spectrum->live_time() > 0.0f)) ? spectrum->live_time() : -1.0;

    // Counts in [lower, upper], or -1 if the range is not all within the spectrum.
    const auto counts_in = [&spectrum]( const double lower, const double upper ) -> double {
      if( !spectrum || !spectrum->channel_energies() || (upper <= lower)
          || (lower < spectrum->gamma_energy_min()) || (upper > spectrum->gamma_energy_max()) )
        return -1.0;
      return spectrum->gamma_integral( static_cast<float>(lower), static_cast<float>(upper) );
    };

    double max_area_above = 0.0;
    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      if( p && p->gausPeak() && (p->mean() > energy) )
        max_area_above = std::max( max_area_above, p->peakArea() );
    }

    for( const shared_ptr<const PeakDef> &p : peaks )
    {
      if( !p || !p->gausPeak() )
        continue;

      const double Ep = p->mean();
      const double area = p->peakArea();
      const double nsigma = (p->peakAreaUncert() > 0.0) ? (area / p->peakAreaUncert()) : 0.0;
      if( (Ep < min_parent_energy) || (Ep <= (energy + window)) || !(area > 0.0) || (nsigma < min_parent_nsigma) )
        continue;

      const double alpha = Ep / electron_rest_mass;
      const double E_bs = Ep / (1.0 + 2.0*alpha);         // 180 deg scatter
      const double E_bs120 = Ep / (1.0 + 1.5*alpha);      // 120 deg: the backscatter bump spans [E_bs, E_bs120]
      const double E_edge = Ep - E_bs;

      const bool near_backscatter = (energy >= (E_bs - window)) && (energy <= (E_bs120 + window));
      const bool near_edge = (std::fabs( energy - E_edge ) <= window);
      if( !near_backscatter && !near_edge )
        continue;

      ComptonKinematicParent parent;
      parent.parentPeak = p;
      parent.parentEnergy = Ep;
      parent.parentArea = area;
      parent.parentCps = (live_time > 0.0) ? (area / live_time) : -1.0;
      if( p->parentNuclide() )
        parent.parentSource = p->parentNuclide()->symbol;
      else if( p->xrayElement() )
        parent.parentSource = p->xrayElement()->symbol;
      else if( p->reaction() )
        parent.parentSource = p->reaction()->name();
      parent.backscatterEnergy = E_bs;
      parent.comptonEdgeEnergy = E_edge;
      parent.window = window;

      if( near_backscatter )
      {
        const double cos_theta = (energy >= E_bs) ? std::clamp( 1.0 - (1.0/alpha)*(Ep/energy - 1.0), -1.0, 1.0 ) : -1.0;
        parent.scatterAngle = std::acos( cos_theta ) * 180.0 / PhysicalUnits::pi;
      }

      if( near_backscatter && near_edge )
        parent.feature = ComptonScatterType::ComptonPeakOrEdge;
      else if( near_edge )
        parent.feature = ComptonScatterType::ComptonEdge;
      else
        parent.feature = ComptonScatterType::ComptonPeak;

      parent.prominence = (max_area_above > 0.0) ? (area / max_area_above) : 0.0;

      // Are the counts at the energy above a straight line through the regions just outside the
      //  parent's scatter features?  (The edge is a step down, so the feature region sits above
      //  such a line.)
      const double feature_lower = std::min( E_bs, E_edge ) - window;
      const double feature_upper = std::max( E_bs120, E_edge ) + window;
      double band_width = 2.0*window;
      if( spectrum && (spectrum->num_gamma_channels() > 1) )
        band_width = std::max( band_width, 4.0*(spectrum->gamma_energy_max() - spectrum->gamma_energy_min())
                                             / spectrum->num_gamma_channels() );
      const double lower_counts = counts_in( feature_lower - band_width, feature_lower );
      const double upper_counts = counts_in( feature_upper, feature_upper + band_width );
      const double observed = counts_in( energy - window, energy + window );
      if( (lower_counts >= 0.0) && (upper_counts >= 0.0) && (observed >= 0.0) )
      {
        const double lower_mid = feature_lower - 0.5*band_width, upper_mid = feature_upper + 0.5*band_width;
        const double frac = std::clamp( (energy - lower_mid) / (upper_mid - lower_mid), 0.0, 1.0 );
        const double density = ((1.0 - frac)*lower_counts + frac*upper_counts) / band_width;
        const double density_var = ((1.0 - frac)*(1.0 - frac)*lower_counts + frac*frac*upper_counts)
                                   / (band_width*band_width);
        const double expected = density * 2.0 * window;
        const double sigma = std::sqrt( std::max( expected, 1.0 ) + 4.0*window*window*density_var );
        parent.excessSigma = (sigma > 0.0) ? ((observed - expected) / sigma) : 0.0;
      }//if( the interpolation can be made )

      if( (parent.prominence >= 0.25) && (parent.excessSigma >= 3.0) )
        parent.confidence = ComptonScatterConfidence::VeryLikely;
      else if( (parent.prominence >= 0.25) || ((parent.prominence >= 0.05) && (parent.excessSigma >= 2.0)) )
        parent.confidence = ComptonScatterConfidence::Likely;
      else
        parent.confidence = ComptonScatterConfidence::Unlikely;

      answer.push_back( std::move(parent) );
    }//for( const shared_ptr<const PeakDef> &p : peaks )

    std::sort( begin(answer), end(answer), []( const ComptonKinematicParent &a, const ComptonKinematicParent &b ) {
      if( a.confidence != b.confidence )
        return static_cast<int>(a.confidence) > static_cast<int>(b.confidence);
      return a.prominence > b.prominence;
    } );

    return answer;
  }//compton_kinematic_parents(...)


  ComptonPeakCheckStatus compton_peak_check( const ComptonPeakCheckOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "compton_peak_check: No InterSpec session available" );

    shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
    if( !meas )
      throw runtime_error( "compton_peak_check: No measurement loaded for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram( options.specType );
    if( !spectrum )
      throw runtime_error( "compton_peak_check: No spectrum displayed for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );

    // The combination of user and auto-search peaks, for both the candidate and its parents.
    DetectedPeaksOptions detectable_peaks_options;
    detectable_peaks_options.specType = options.specType;
    detectable_peaks_options.nonBackgroundPeaksOnly = false;
    const DetectedPeakStatus detectable_peaks = detected_peaks( detectable_peaks_options, interspec );

    auto all_peaks = make_shared<const deque<shared_ptr<const PeakDef>>>( std::begin(detectable_peaks.peaks),
                                                                          std::end(detectable_peaks.peaks) );

    // The detector resolution at the energy, checked against this spectrum's peaks.
    const std::function<float(double)> fwhm_fcn = make_expected_fwhm_fcn( meas->detector(), all_peaks, spectrum,
                                                                          resolve_det_type( interspec, meas, spectrum ) );
    const double expected_fwhm = fwhm_fcn( options.energy );

    // The candidate: the Gaussian peak nearest the energy.
    shared_ptr<const PeakDef> nearest;
    double smallest_diff = std::numeric_limits<double>::max();
    for( const shared_ptr<const PeakDef> &p : *all_peaks )
    {
      if( !p || !p->gausPeak() )
        continue;
      const double diff = std::abs( p->mean() - options.energy );
      if( diff < smallest_diff )
      {
        smallest_diff = diff;
        nearest = p;
      }
    }

    const double tol = std::max( 0.5, nearest ? nearest->fwhm() : expected_fwhm );
    const bool have_candidate = (nearest && (smallest_diff <= tol));

    ComptonPeakCheckStatus result;
    if( have_candidate )
      result = compton_peak_check( nearest, all_peaks, spectrum, meas->detector(), fwhm_fcn );

    result.searchWindow = tol;
    result.noCandidatePeak = !have_candidate;
    if( nearest )
    {
      result.nearestPeakEnergy = nearest->mean();
      result.nearestPeakDistance = smallest_diff;
    }

    // From kinematics alone - which is all there is to go on for a broad feature never fit as a peak.
    const double kinematic_energy = have_candidate ? nearest->mean() : options.energy;
    const double kinematic_window = std::max( (expected_fwhm > 0.0) ? expected_fwhm : 1.0,
                                              have_candidate ? nearest->fwhm() : 0.0 );
    result.kinematicParents = compton_kinematic_parents( kinematic_energy, kinematic_window,
                                                         detectable_peaks.peaks, spectrum );

    // The verdict on a fitted candidate is about that peak (a narrow photopeak can sit in a
    //  backscatter region); say so when the region alone looks more like a scatter feature.
    if( have_candidate && !result.kinematicParents.empty() )
    {
      const ComptonKinematicParent &region = result.kinematicParents.front();
      if( (static_cast<int>(region.confidence) >= static_cast<int>(ComptonScatterConfidence::Likely))
          && (static_cast<int>(result.confidence) <= static_cast<int>(ComptonScatterConfidence::Unlikely)) )
      {
        char buffer[384] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%sthe energy is in the backscatter/Compton-edge region of the %.2f keV peak"
                  " (backscatter %.1f keV, edge %.1f keV; see kinematicParents), but the fitted peak here is judged on"
                  " its own width and area, which look %s a scatter feature", (result.reason.empty() ? "" : "; "),
                  region.parentEnergy, region.backscatterEnergy, region.comptonEdgeEnergy,
                  (result.confidence == ComptonScatterConfidence::No) ? "unlike" : "only somewhat like" );
        result.reason += buffer;
      }
    }//if( a candidate was judged, and the region suggests a scatter feature )

    if( !have_candidate )
    {
      char buffer[256] = { '\0' };
      if( nearest )
        snprintf( buffer, sizeof(buffer), "No fitted peak within %.2f keV of %.2f keV (the nearest is %.2f keV away),"
                  " so only the Compton kinematics of the higher-energy peaks were checked", tol, options.energy, smallest_diff );
      else
        snprintf( buffer, sizeof(buffer), "No fitted peaks in the spectrum, so only Compton kinematics were checked" );
      result.reason = buffer;

      if( result.kinematicParents.empty() )
      {
        result.reason += "; no higher-energy peak has a backscatter peak or Compton edge near this energy.";
      }else
      {
        const ComptonKinematicParent &best = result.kinematicParents.front();
        result.type = best.feature;
        result.confidence = best.confidence;
        snprintf( buffer, sizeof(buffer), "; the energy is in the backscatter/Compton-edge region of the %.2f keV"
                  " peak (backscatter %.1f keV, edge %.1f keV).", best.parentEnergy, best.backscatterEnergy,
                  best.comptonEdgeEnergy );
        result.reason += buffer;
      }
    }//if( !have_candidate )

    return result;
  }//ComptonPeakCheckStatus compton_peak_check( const ComptonPeakCheckOptions &options, InterSpec *interspec )


  // ---------------------------------------------------------------------------------------------
  // beta_continuum_check: deterministic pure-beta-emitter (bremsstrahlung) continuum analysis.
  //
  // Thresholds below were tuned against the ~800-spectrum injected benchmark (35 detector configs,
  // NaI through HPGe; P-32/Sr-90/Tl-204/Y-90 vs shielded-gamma negatives and background nulls);
  // kept as named file-locals so future tuning is a one-line change.
  // ---------------------------------------------------------------------------------------------
  static constexpr double sk_beta_coarse_bin_kev      = 24.0;   // coarse analysis bin width
  static constexpr double sk_beta_xray_max_kev        = 120.0;  // peaks below this are reported as x-rays, and dont veto
  static constexpr double sk_beta_analysis_max_kev    = 3000.0; // upper end of analysis range
  static constexpr double sk_beta_mode_max_kev        = 1500.0; // continuum hump mode must be below this
  static constexpr double sk_beta_min_elevation_sigma = 8.0;    // below this: NoContinuumElevation
  static constexpr double sk_beta_peak_veto_sigma     = 5.0;    // gamma-peak elevation bar for the veto
  static constexpr double sk_beta_term_ratio_max      = 0.05;   // cliff detector: intensity at termination / max
  static constexpr double sk_beta_echar_max_kev       = 350.0;  // max brem-like characteristic (1/e) energy
  static constexpr double sk_beta_term_max_kev        = 2450.0; // max plausible beta endpoint (Y-90 is 2280 keV)
  static constexpr double sk_beta_seg_resid_max_sigma = 8.0;    // smoothness: max |segment residual|
  static constexpr double sk_beta_outlier_frac_max    = 0.25;   // smoothness: max fit-outlier fraction
  static constexpr double sk_beta_sys_frac            = 0.03;   // fractional systematic floor on the model
  static constexpr double sk_beta_511_dominant_sigma  = 8.0;    // 511 significance for the positron demotion
  static constexpr double sk_beta_max_fwhm_frac       = 0.12;   // cap on expectedFwhm(E)/E - protects against a
                                                                //  garbage resolution fit (e.g. from two spurious
                                                                //  peaks) blowing up the exclusion regions
  static constexpr double sk_beta_fit_count_quantile  = 0.99;   // shape fitting/gating is bounded at this net-count
                                                                //  quantile; a trace tail beyond it (commonly random
                                                                //  summing / pile-up) is reported but must not drive
                                                                //  the verdict or the endpoint lower bound
  static constexpr double sk_beta_cliff_meaningful_sigma = 9.0; // a cliff-like termination ratio is categorical
                                                                //  evidence only when the smoothed intensity at the
                                                                //  termination is this many sigma above zero (3x the
                                                                //  3-sigma detection floor); otherwise termination is
                                                                //  statistics-limited and the ratio is inconclusive
  static constexpr double sk_beta_endpoint_slop_frac  = 0.05;   // candidate-nuclide matching slop, as a fraction of
                                                                //  the termination energy - resolution smearing and
                                                                //  residual summing can push the observed termination
                                                                //  slightly past the true endpoint (e.g. Y-90 phantom
                                                                //  terminating at ~2390 vs the 2280 keV endpoint)

  /** (Nearly) pure beta-minus emitters considered as consistency candidates; endpoints are pulled
   from SandiaDecay at runtime (including in-equilibrium daughters, e.g. Sr90 -> Y90). */
  static const char * const sk_pure_beta_nuclide_names[] = {
    "H3", "C14", "P32", "S35", "Cl36", "Ca45", "Ni63", "Sr89", "Sr90", "Y90",
    "Tc99", "Ru106", "Pm147", "Tl204", "Bi210"
  };

  /** One coarse bin of the live-time-normalized, background-subtracted continuum. */
  struct BetaCoarseBin
  {
    double center = 0.0;    // keV
    double net = 0.0;       // foreground - scaled background, counts
    double var = 0.0;       // variance of net
    bool excluded = false;  // below the x-ray floor, or inside a peak-exclusion region
  };//struct BetaCoarseBin

  /** Max beta endpoint (keV) over `nuc` and any descendants with half-life shorter than the
   parent (i.e., nuclides that reach equilibrium with it); `note` names the descendant when
   the endpoint comes from one.  Returns 0 if no beta branch found. */
  static double max_equilibrium_beta_endpoint( const SandiaDecay::Nuclide * const nuc, std::string &note )
  {
    note.clear();
    if( !nuc )
      return 0.0;

    double best_kev = 0.0;
    const SandiaDecay::Nuclide *best_nuc = nullptr;

    const vector<const SandiaDecay::Nuclide *> descendants = nuc->descendants();
    for( const SandiaDecay::Nuclide * const d : descendants )
    {
      if( !d || ((d != nuc) && (d->halfLife > nuc->halfLife)) )
        continue;

      for( const SandiaDecay::Transition * const trans : d->decaysToChildren )
      {
        if( !trans )
          continue;

        for( const SandiaDecay::RadParticle &particle : trans->products )
        {
          // For beta particles, RadParticle::energy is the endpoint energy (in SandiaDecay units of keV)
          if( (particle.type == SandiaDecay::BetaParticle) && (particle.energy > best_kev) )
          {
            best_kev = particle.energy;
            best_nuc = d;
          }
        }//for( particles )
      }//for( transitions )
    }//for( descendants )

    if( best_nuc && (best_nuc != nuc) )
      note = "endpoint is from in-equilibrium daughter " + best_nuc->symbol;

    return best_kev;
  }//max_equilibrium_beta_endpoint(...)


  /** Rebin `meas` onto the uniform grid given by `edges` (channel lower energies in keV, plus the
   final upper edge); returns per-bin counts (clamped non-negative). */
  static vector<double> beta_rebin_counts( const shared_ptr<const SpecUtils::Measurement> &meas,
                                           const vector<float> &edges )
  {
    assert( meas && (edges.size() > 1) );

    const shared_ptr<const vector<float>> &counts = meas->gamma_counts();
    const shared_ptr<const vector<float>> &energies = meas->gamma_channel_energies();

    // rebin_by_lower_edge requires at least 4 channels in both inputs.
    if( !counts || !energies || (counts->size() < 4) || (energies->size() < 4) )
      throw runtime_error( "beta_continuum_check: spectrum has too few channels" );

    vector<float> rebinned;
    SpecUtils::rebin_by_lower_edge( *energies, *counts, edges, rebinned );

    vector<double> result( edges.size() - 1, 0.0 );
    const size_t nfill = std::min( result.size(), rebinned.size() );
    for( size_t i = 0; i < nfill; ++i )
      result[i] = std::max( 0.0, static_cast<double>(rebinned[i]) );

    return result;
  }//beta_rebin_counts(...)


  /** A [lower,upper] keV interval excluded from continuum-shape analysis (peak regions, 511). */
  struct BetaExclusionInterval
  {
    double lower = 0.0;
    double upper = 0.0;
  };//struct BetaExclusionInterval

  /** Builds the merged, sorted peak-exclusion intervals: +-2 FWHM around every foreground peak,
   plus the 511 keV annihilation region.

   Background-only peaks are deliberately NOT excluded: the analysis runs on the live-time
   normalized net spectrum, where background lines cancel (their statistical residual is already
   in the per-bin variance), and any background line that does not cancel will have been found as
   a foreground peak anyway.  Excluding background peaks would tile most of a low-resolution
   spectrum with exclusions (wide NORM peaks every few hundred keV), blinding the tail analysis.
   */
  static vector<BetaExclusionInterval> beta_exclusion_intervals(
                                  const vector<shared_ptr<const PeakDef>> &fore_peaks,
                                  const std::function<float(double)> &fwhm_fcn )
  {
    vector<BetaExclusionInterval> intervals;

    const auto add_peak = [&intervals, &fwhm_fcn]( const shared_ptr<const PeakDef> &p ){
      if( !p )
        return;
      const double mean = p->mean();
      double fwhm = p->gausPeak() ? p->fwhm() : 0.0;
      if( fwhm <= 0.0 )
        fwhm = std::max( 1.0, std::min( static_cast<double>(fwhm_fcn( mean )),
                                        sk_beta_max_fwhm_frac*mean ) );
      intervals.push_back( { mean - 2.0*fwhm, mean + 2.0*fwhm } );
    };

    for( const shared_ptr<const PeakDef> &p : fore_peaks )
      add_peak( p );

    // Always exclude the annihilation region - 511 is allowed for a beta source (pair production
    // of high-energy brems, or a small beta+ branch), so it must not distort the shape fit.
    const double fwhm511 = std::max( 1.0, std::min( static_cast<double>(fwhm_fcn( 511.0 )),
                                                    sk_beta_max_fwhm_frac*511.0 ) );
    const double half511 = std::max( 2.0*fwhm511, 10.0 );
    intervals.push_back( { 511.0 - half511, 511.0 + half511 } );

    std::sort( begin(intervals), end(intervals),
               []( const BetaExclusionInterval &a, const BetaExclusionInterval &b ){
                 return a.lower < b.lower;
               } );

    vector<BetaExclusionInterval> merged;
    for( const BetaExclusionInterval &iv : intervals )
    {
      if( merged.empty() || (iv.lower > merged.back().upper) )
        merged.push_back( iv );
      else
        merged.back().upper = std::max( merged.back().upper, iv.upper );
    }

    return merged;
  }//beta_exclusion_intervals(...)


  /** Result of the iterative outlier-robust weighted polynomial fit of ln(net) vs energy. */
  struct BetaLogFitResult
  {
    bool valid = false;
    double chi2Dof = 0.0;
    double outlierFraction = 0.0;
    double maxSegmentResidual = 0.0;      // max |summed residual| over 5 segments, in sigma
    double slopePerKev = 0.0;             // d ln(net) / dE at the requested midpoint energy
  };//struct BetaLogFitResult

  /** Iterative (3-sigma rejection) weighted cubic (quadratic if few bins) fit of ln(net) vs
   energy over `fit_indices` of `bins`; residuals are normalized with a `sk_beta_sys_frac`
   fractional systematic floor added in quadrature, so high-statistics spectra are not
   penalized for percent-level wiggles. */
  static BetaLogFitResult beta_robust_log_fit( const vector<BetaCoarseBin> &bins,
                                               const vector<size_t> &fit_indices,
                                               const double emid_kev )
  {
    BetaLogFitResult result;

    const size_t nfit = fit_indices.size();
    if( nfit < 6 )
      return result;

    const int order = (nfit >= 10) ? 3 : 2;

    // Energies are scaled to MeV in the design matrix for conditioning.
    const auto poly_ln = [order]( const Eigen::VectorXd &coefs, const double e_kev ) -> double {
      const double x = e_kev / 1000.0;
      double ln_val = 0.0, xk = 1.0;
      for( int k = 0; k <= order; ++k, xk *= x )
        ln_val += coefs[k] * xk;
      return std::max( -50.0, std::min( 50.0, ln_val ) );
    };

    vector<char> live( bins.size(), 0 );
    for( const size_t idx : fit_indices )
      live[idx] = 1;

    Eigen::VectorXd coefs;

    for( int iter = 0; iter < 5; ++iter )
    {
      vector<size_t> pts;
      for( const size_t idx : fit_indices )
      {
        if( live[idx] )
          pts.push_back( idx );
      }

      if( pts.size() < 6 )
        break;  //keep coefficients from the previous iteration

      Eigen::MatrixXd design( pts.size(), order + 1 );
      Eigen::VectorXd rhs( pts.size() );
      for( size_t row = 0; row < pts.size(); ++row )
      {
        const BetaCoarseBin &bin = bins[pts[row]];
        // Weight of ln(net): sigma_ln = sqrt(var)/net, so sqrt(weight) = net/sqrt(var)
        const double sqrt_w = bin.net / std::sqrt( std::max( bin.var, 1.0e-9 ) );
        const double x = bin.center / 1000.0;
        double xk = 1.0;
        for( int k = 0; k <= order; ++k, xk *= x )
          design(row, k) = sqrt_w * xk;
        rhs[row] = sqrt_w * std::log( bin.net );
      }

      coefs = design.colPivHouseholderQr().solve( rhs );

      bool changed = false;
      for( const size_t idx : fit_indices )
      {
        if( !live[idx] )
          continue;  //once rejected, stays rejected
        const BetaCoarseBin &bin = bins[idx];
        const double model = std::exp( poly_ln( coefs, bin.center ) );
        const double sig_eff = std::sqrt( bin.var + std::pow( sk_beta_sys_frac * model, 2.0 ) );
        if( std::fabs( (bin.net - model) / sig_eff ) > 3.0 )
        {
          live[idx] = 0;
          changed = true;
        }
      }//for( fit_indices )

      if( !changed )
        break;
    }//for( iter )

    if( coefs.size() != (order + 1) )
      return result;

    result.valid = true;

    // chi2/dof over surviving bins, outlier fraction, and segment systematics over all fit bins.
    size_t nlive = 0;
    double chi2 = 0.0;
    for( const size_t idx : fit_indices )
    {
      if( !live[idx] )
        continue;
      const BetaCoarseBin &bin = bins[idx];
      const double model = std::exp( poly_ln( coefs, bin.center ) );
      const double sig_eff = std::sqrt( bin.var + std::pow( sk_beta_sys_frac * model, 2.0 ) );
      chi2 += std::pow( (bin.net - model) / sig_eff, 2.0 );
      nlive += 1;
    }
    result.chi2Dof = chi2 / std::max( 1.0, static_cast<double>(nlive) - order - 1.0 );
    result.outlierFraction = 1.0 - (static_cast<double>(nlive) / static_cast<double>(nfit));

    // Max |summed residual| over 5 contiguous segments - catches broad systematic structure
    // (Compton bumps, backscatter humps) that per-bin outlier rejection can miss.
    const size_t nseg = 5;
    for( size_t seg = 0; seg < nseg; ++seg )
    {
      const size_t seg_begin = (seg * nfit) / nseg;
      const size_t seg_end = ((seg + 1) * nfit) / nseg;
      double sum_resid = 0.0, sum_var = 0.0;
      for( size_t i = seg_begin; i < seg_end; ++i )
      {
        const BetaCoarseBin &bin = bins[fit_indices[i]];
        const double model = std::exp( poly_ln( coefs, bin.center ) );
        sum_resid += (bin.net - model);
        sum_var += bin.var + std::pow( sk_beta_sys_frac * model, 2.0 );
      }
      if( sum_var > 0.0 )
        result.maxSegmentResidual = std::max( result.maxSegmentResidual,
                                              std::fabs( sum_resid ) / std::sqrt( sum_var ) );
    }//for( segments )

    // Local slope of ln(net) at the fit-range midpoint (per keV).
    const double xmid = emid_kev / 1000.0;
    double slope_per_mev = 0.0, xk = 1.0;
    for( int k = 1; k <= order; ++k, xk *= xmid )
      slope_per_mev += k * coefs[k] * xk;
    result.slopePerKev = slope_per_mev / 1000.0;

    return result;
  }//beta_robust_log_fit(...)


  const char *to_string( const BetaContinuumVerdict verdict )
  {
    switch( verdict )
    {
      case BetaContinuumVerdict::BackgroundNotLoaded:    return "BackgroundNotLoaded";
      case BetaContinuumVerdict::NoContinuumElevation:   return "NoContinuumElevation";
      case BetaContinuumVerdict::BremLike:               return "BremLike";
      case BetaContinuumVerdict::Ambiguous:              return "Ambiguous";
      case BetaContinuumVerdict::NotBremLike:            return "NotBremLike";
      case BetaContinuumVerdict::AnnihilationDominated:  return "AnnihilationDominated";
    }
    assert( 0 );
    return "InvalidBetaContinuumVerdict";
  }//to_string( BetaContinuumVerdict )


  BetaContinuumCheckStatus beta_continuum_check( const BetaContinuumCheckInput &input )
  {
    BetaContinuumCheckStatus status;

    const shared_ptr<const SpecUtils::Measurement> &fore = input.foreground;
    if( !fore || !fore->gamma_counts() || (fore->num_gamma_channels() < 4) )
      throw runtime_error( "beta_continuum_check: no usable foreground spectrum" );

    if( !input.expectedFwhm )
      throw runtime_error( "beta_continuum_check: no expected-FWHM function provided" );

    const shared_ptr<const SpecUtils::Measurement> &back = input.background;
    if( !back || !back->gamma_counts() || (back->num_gamma_channels() < 4) )
    {
      status.verdict = BetaContinuumVerdict::BackgroundNotLoaded;
      status.warnings.push_back( "A background spectrum is required to test for a beta-emitter"
                                 " (bremsstrahlung) continuum; please load a background measurement"
                                 " taken with the same detector." );
      return status;
    }

    const std::function<float(double)> &fwhm_fcn = input.expectedFwhm;

    const double fore_lt = (fore->live_time() > 0.0f) ? fore->live_time() : 1.0;
    const double back_lt = (back->live_time() > 0.0f) ? back->live_time() : 1.0;
    const double lt_sf = fore_lt / back_lt;

    // ---- Classify foreground peaks: x-ray region / annihilation / gamma-source evidence ----
    // FWHM capped at a physical fraction of energy, in case the resolution estimate is garbage
    // (e.g. fit through only a couple spurious wide peaks).
    const double fwhm511 = std::max( 1.0, std::min( static_cast<double>(fwhm_fcn( 511.0 )),
                                                    sk_beta_max_fwhm_frac*511.0 ) );
    double sig511 = 0.0;

    for( const shared_ptr<const PeakDef> &peak : input.foregroundPeaks )
    {
      if( !peak || !peak->gausPeak() )
        continue;

      const double mean = peak->mean();
      const double fg_cps = peak->amplitude() / fore_lt;

      // Elevation significance over the closest matching background peak (if any); same matching
      // and arithmetic as detected_peaks, but with the fractional systematic floor included so
      // huge-statistics spectra dont produce meaninglessly-large sigmas.
      shared_ptr<const PeakDef> closest_bg;
      double smallest_de = std::numeric_limits<double>::infinity();
      for( const shared_ptr<const PeakDef> &bg_peak : input.backgroundPeaks )
      {
        if( !bg_peak || !bg_peak->gausPeak() )
          continue;
        const double de = fabs( mean - bg_peak->mean() );
        const double avg_fwhm = 0.75 * (peak->fwhm() + bg_peak->fwhm()) / 2.0;
        if( (de < avg_fwhm) && (de < smallest_de) )
        {
          closest_bg = bg_peak;
          smallest_de = de;
        }
      }//for( background peaks )

      const double bg_cps = closest_bg ? (closest_bg->amplitude() / back_lt) : 0.0;
      const double fg_amp_uncert = elevation_test_uncert( *peak, input.foregroundPeaks, fore );
      const double bg_amp_uncert = closest_bg ? elevation_test_uncert( *closest_bg, input.backgroundPeaks, back ) : 0.0;
      const double fg_cps_uncert = (fg_amp_uncert > 0.0)
                                     ? (fg_amp_uncert / fore_lt)
                                     : (std::sqrt( std::max( peak->amplitude(), 1.0 ) ) / fore_lt);
      const double bg_cps_uncert = (closest_bg && (bg_amp_uncert > 0.0))
                                     ? (bg_amp_uncert / back_lt)
                                     : (closest_bg ? (std::sqrt( std::max( closest_bg->amplitude(), 1.0 ) ) / back_lt) : 0.0);
      const double combined_uncert = std::sqrt( fg_cps_uncert*fg_cps_uncert
                                                + bg_cps_uncert*bg_cps_uncert
                                                + std::pow( sk_beta_sys_frac * fg_cps, 2.0 ) );
      const double sigma_elev = (combined_uncert > 0.0) ? ((fg_cps - bg_cps) / combined_uncert) : 0.0;

      BetaContinuumPeakInfo info;
      info.energy = mean;
      info.amplitude = peak->amplitude();
      info.numSigmaElevated = sigma_elev;
      info.isAnnihilation = (fabs( mean - 511.0 ) < 1.5*fwhm511);

      if( info.isAnnihilation )
        sig511 = std::max( sig511, sigma_elev );
      else if( mean <= sk_beta_xray_max_kev )
        status.xrayRegionPeaks.push_back( info );
      else if( sigma_elev > sk_beta_peak_veto_sigma )
        status.elevatedGammaPeaks.push_back( info );
    }//for( foreground peaks )

    // ---- Coarse net-continuum grid, with peak regions excluded ----
    const double upper_kev = std::min( sk_beta_analysis_max_kev,
                                       static_cast<double>(fore->gamma_energy_max()) );
    if( upper_kev < (sk_beta_xray_max_kev + 10.0*sk_beta_coarse_bin_kev) )
      throw runtime_error( "beta_continuum_check: spectrum energy range (up to "
                           + std::to_string( static_cast<int>(upper_kev) )
                           + " keV) is too small to analyze the continuum shape" );

    const size_t nbins = static_cast<size_t>( std::floor( upper_kev / sk_beta_coarse_bin_kev ) );
    vector<float> edges( nbins + 1 );
    for( size_t i = 0; i <= nbins; ++i )
      edges[i] = static_cast<float>( i * sk_beta_coarse_bin_kev );

    const vector<double> fore_binned = beta_rebin_counts( fore, edges );
    const vector<double> back_binned = beta_rebin_counts( back, edges );

    const vector<BetaExclusionInterval> exclusions
        = beta_exclusion_intervals( input.foregroundPeaks, fwhm_fcn );

    vector<BetaCoarseBin> bins( nbins );
    for( size_t i = 0; i < nbins; ++i )
    {
      BetaCoarseBin &bin = bins[i];
      bin.center = (i + 0.5) * sk_beta_coarse_bin_kev;
      bin.net = fore_binned[i] - lt_sf*back_binned[i];
      bin.var = std::max( fore_binned[i], 1.0 ) + lt_sf*lt_sf*std::max( back_binned[i], 1.0 );
      bin.excluded = (bin.center <= sk_beta_xray_max_kev);
      for( size_t j = 0; !bin.excluded && (j < exclusions.size()); ++j )
        bin.excluded = ((bin.center >= exclusions[j].lower) && (bin.center <= exclusions[j].upper));
    }//for( coarse bins )

    // ---- Total continuum-elevation significance ----
    double sum_net = 0.0, sum_var = 0.0;
    for( const BetaCoarseBin &bin : bins )
    {
      if( bin.excluded )
        continue;
      sum_net += bin.net;
      sum_var += bin.var;
    }
    status.continuumElevationSigma = sum_net / std::sqrt( std::max( sum_var, 1.0e-9 ) );

    if( status.continuumElevationSigma < sk_beta_min_elevation_sigma )
    {
      status.verdict = BetaContinuumVerdict::NoContinuumElevation;
      return status;
    }

    // ---- Smoothed profile: hump mode, termination energy, termination (cliff) ratio ----
    vector<double> smooth( nbins, 0.0 ), smooth_var( nbins, 0.0 );
    for( size_t i = 0; i < nbins; ++i )
    {
      const size_t lo = (i >= 2) ? (i - 2) : size_t(0);
      const size_t hi = std::min( i + 2, nbins - 1 );
      const double cnt = static_cast<double>( hi - lo + 1 );
      for( size_t j = lo; j <= hi; ++j )
      {
        smooth[i] += bins[j].net / cnt;
        smooth_var[i] += bins[j].var / (cnt*cnt);
      }
    }//for( coarse bins )

    size_t imode = nbins;
    for( size_t i = 0; i < nbins; ++i )
    {
      if( bins[i].excluded || (bins[i].center >= sk_beta_mode_max_kev) )
        continue;
      if( (imode == nbins) || (smooth[i] > smooth[imode]) )
        imode = i;
    }

    if( imode == nbins )
    {
      status.verdict = BetaContinuumVerdict::Ambiguous;
      status.warnings.push_back( "Continuum is elevated, but no usable (non-peak) region was found"
                                 " to locate the continuum maximum - too many peak regions?" );
      return status;
    }

    status.modeEnergy = bins[imode].center;

    size_t iterm_raw = imode;
    for( size_t i = imode + 1; i < nbins; ++i )
    {
      if( !bins[i].excluded && (smooth[i] > 3.0*std::sqrt( std::max( smooth_var[i], 1.0e-9 ) )) )
        iterm_raw = i;
    }

    // Bound the shape analysis at the sk_beta_fit_count_quantile net-count quantile: a trace
    // high-energy tail (random-summing/pile-up, or tiny residuals) can be statistically
    // significant while holding <1% of the counts, and letting it stretch the fit range makes
    // the shape gates grade the tail instead of the continuum the counts are actually in.
    double total_pos_net = 0.0;
    for( const BetaCoarseBin &bin : bins )
    {
      if( !bin.excluded )
        total_pos_net += std::max( bin.net, 0.0 );
    }
    size_t iquant = iterm_raw;
    double cum_net = 0.0;
    for( size_t i = 0; i < nbins; ++i )
    {
      if( bins[i].excluded )
        continue;
      cum_net += std::max( bins[i].net, 0.0 );
      if( cum_net >= sk_beta_fit_count_quantile*total_pos_net )
      {
        iquant = i;
        break;
      }
    }//for( coarse bins )

    // ---- Two-pass shape evaluation ----
    // The raw (last-significant-bin) termination is tried first, preserving full-range grading
    // for clean spectra where the high tail is simply more continuum.  When a weak tail extends
    // past the count-quantile bound, a second, clamped pass is tried: a trace tail (commonly
    // random-summing/pile-up) can be statistically significant while holding <1% of the counts,
    // and letting it stretch the fit range makes the gates grade the tail instead of the
    // continuum the counts are actually in.  The better verdict wins.
    struct BetaRangeEval
    {
      bool valid = false;
      int rank = -1;   // 2 = BremLike, 1 = Ambiguous, 0 = NotBremLike
      BetaContinuumVerdict verdict = BetaContinuumVerdict::NotBremLike;
      size_t iterm = 0;
      double termRatio = 0.0;
      std::optional<double> eChar;
      BetaLogFitResult fit;
      bool echarOk = false;
    };//struct BetaRangeEval

    const bool gate_peaks = status.elevatedGammaPeaks.empty();

    const auto evaluate_range = [&]( size_t iterm_eval ) -> BetaRangeEval {
      BetaRangeEval eval;
      while( (iterm_eval > imode) && bins[iterm_eval].excluded )
        iterm_eval -= 1;
      if( iterm_eval <= (imode + 3) )
        return eval;  //too little range beyond the maximum to analyze

      eval.iterm = iterm_eval;
      eval.termRatio = std::max( smooth[iterm_eval], 0.0 ) / std::max( smooth[imode], 1.0e-9 );

      vector<size_t> fit_indices;
      for( size_t i = imode + 1; i <= iterm_eval; ++i )
      {
        if( !bins[i].excluded && (bins[i].net > 0.0) )
          fit_indices.push_back( i );
      }

      const double emid_kev = 0.5*(bins[imode].center + bins[iterm_eval].center);
      eval.fit = beta_robust_log_fit( bins, fit_indices, emid_kev );
      if( !eval.fit.valid )
        return eval;

      eval.valid = true;
      if( eval.fit.slopePerKev < 0.0 )
        eval.eChar = -1.0 / eval.fit.slopePerKev;

      // A large termination ratio only counts as a hard (photopeak/Compton-edge style) cliff
      // when the intensity at termination is well above the detection floor; for weak spectra
      // the 3-sigma termination is statistics-limited and the ratio is meaningless.
      const bool cliff_bad = (eval.termRatio >= sk_beta_term_ratio_max);
      const bool cliff_meaningful = (smooth[iterm_eval] > (sk_beta_cliff_meaningful_sigma
                          * std::sqrt( std::max( smooth_var[iterm_eval], 1.0e-9 ) )));

      const bool gate_term = (bins[iterm_eval].center < sk_beta_term_max_kev);
      const bool gate_falling = eval.eChar.has_value();
      eval.echarOk = gate_falling && (*eval.eChar < sk_beta_echar_max_kev);
      const bool gate_smooth = (eval.fit.maxSegmentResidual < sk_beta_seg_resid_max_sigma)
                               && (eval.fit.outlierFraction < sk_beta_outlier_frac_max);

      // Hard gates are categorical evidence AGAINST a beta source (gamma peaks, a real cliff,
      // termination beyond any plausible beta endpoint, a non-falling continuum); any failure
      // => NotBremLike.  Soft gates are shape-quality evidence (how exponential-clean the fall
      // is); a failure only demotes to Ambiguous - heavily scattered/hardened but genuine beta
      // continua (e.g. Y-90 in a body phantom) fail these without being gamma-like.
      const bool hard_ok = gate_peaks && gate_term && gate_falling
                           && !(cliff_bad && cliff_meaningful);
      const bool soft_ok = eval.echarOk && gate_smooth && !cliff_bad;

      if( !hard_ok )
      {
        eval.verdict = BetaContinuumVerdict::NotBremLike;
        eval.rank = 0;
      }else if( soft_ok )
      {
        eval.verdict = BetaContinuumVerdict::BremLike;
        eval.rank = 2;
      }else
      {
        eval.verdict = BetaContinuumVerdict::Ambiguous;
        eval.rank = 1;
      }

      return eval;
    };//evaluate_range lambda

    const BetaRangeEval full_eval = evaluate_range( iterm_raw );

    const size_t iterm_clamped = std::min( iterm_raw, std::max( iquant, imode + 4 ) );
    BetaRangeEval clamped_eval;
    if( (iterm_clamped + 2) < iterm_raw )
      clamped_eval = evaluate_range( iterm_clamped );

    const bool use_clamped = clamped_eval.valid
                             && (!full_eval.valid || (clamped_eval.rank > full_eval.rank));
    const BetaRangeEval &chosen = use_clamped ? clamped_eval : full_eval;

    if( !chosen.valid )
    {
      status.verdict = BetaContinuumVerdict::Ambiguous;
      status.warnings.push_back( "Continuum is elevated, but extends too little beyond its"
                                 " maximum (or has too few usable bins) to analyze its shape." );
      return status;
    }

    status.terminationEnergy = bins[chosen.iterm].center;
    status.terminationRatio = chosen.termRatio;
    status.fitChi2Dof = chosen.fit.chi2Dof;
    status.fitOutlierFraction = chosen.fit.outlierFraction;
    status.fitMaxSegmentResidual = chosen.fit.maxSegmentResidual;
    status.characteristicEnergy = chosen.eChar;
    status.verdict = chosen.verdict;

    // Report (but dont grade on) the weak tail beyond the fitted range, when clamping mattered.
    if( use_clamped && (iterm_raw > (chosen.iterm + 2)) )
    {
      double tail_net = 0.0, tail_var = 0.0;
      for( size_t i = chosen.iterm + 1; i <= iterm_raw; ++i )
      {
        if( bins[i].excluded )
          continue;
        tail_net += bins[i].net;
        tail_var += bins[i].var;
      }
      const double tail_sigma = tail_net / std::sqrt( std::max( tail_var, 1.0e-9 ) );
      if( tail_sigma > 2.0 )
      {
        status.highEnergyTailUpperEnergy = bins[iterm_raw].center;
        status.highEnergyTailSigma = tail_sigma;

        const double dead_time_frac = ((fore->real_time() > 0.0f) && (fore->live_time() > 0.0f)
                                       && (fore->real_time() > fore->live_time()))
                    ? (1.0 - fore->live_time()/fore->real_time()) : 0.0;
        char buffer[512];
        snprintf( buffer, sizeof(buffer),
                  "A weak high-energy net tail (%.1f sigma, <%.0f%% of net counts) extends beyond"
                  " the fitted continuum, up to ~%.0f keV.  It was not used for shape grading or"
                  " the beta-endpoint lower bound; such a tail is commonly random-summing"
                  " (pile-up)%s.",
                  tail_sigma, 100.0*(1.0 - sk_beta_fit_count_quantile),
                  bins[iterm_raw].center,
                  (dead_time_frac > 0.005) ? " - consistent with the foreground dead time" : "" );
        status.warnings.push_back( buffer );
      }//if( tail is significant )
    }//if( clamped pass chosen and a tail extends beyond it )

    if( (status.verdict == BetaContinuumVerdict::Ambiguous) && !chosen.echarOk )
      status.warnings.push_back( "The continuum falls more slowly than typical unshielded-beta"
                                 " bremsstrahlung; this can be a heavily scattered/shielded hard"
                                 " beta source (e.g. Y-90 through a person or thick shielding),"
                                 " but also a peak-less scattered gamma continuum." );

    // 511-dominance: a strong annihilation peak with the continuum dying just above it suggests
    // a positron source (e.g. F-18) whose "continuum" is 511-scatter, not beta-minus brems.
    status.annihilationDominant = (sig511 > sk_beta_511_dominant_sigma)
                                  && (*status.terminationEnergy < (511.0 + 5.0*fwhm511));
    if( status.annihilationDominant && (status.verdict == BetaContinuumVerdict::BremLike) )
      status.verdict = BetaContinuumVerdict::AnnihilationDominated;

    // ---- Candidate pure-beta nuclides consistent with the observed termination energy ----
    if( (status.verdict == BetaContinuumVerdict::BremLike)
        || (status.verdict == BetaContinuumVerdict::Ambiguous)
        || (status.verdict == BetaContinuumVerdict::AnnihilationDominated) )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      if( db )
      {
        for( const char * const name : sk_pure_beta_nuclide_names )
        {
          const SandiaDecay::Nuclide * const nuc = db->nuclide( name );
          std::string note;
          const double endpoint = max_equilibrium_beta_endpoint( nuc, note );
          // terminationEnergy is a lower bound on the endpoint, but resolution smearing and
          // residual summing can push it slightly past the true endpoint - allow proportional slop.
          const double slop = std::max( sk_beta_coarse_bin_kev,
                                        sk_beta_endpoint_slop_frac * (*status.terminationEnergy) );
          if( nuc && (endpoint >= (*status.terminationEnergy - slop)) )
          {
            BetaContinuumCandidateNuclide candidate;
            candidate.nuclide = nuc->symbol;
            candidate.betaEndpoint = endpoint;
            candidate.note = note;
            status.candidateNuclides.push_back( candidate );
          }
        }//for( candidate nuclide names )

        std::sort( begin(status.candidateNuclides), end(status.candidateNuclides),
                   []( const BetaContinuumCandidateNuclide &a, const BetaContinuumCandidateNuclide &b ){
                     return a.betaEndpoint < b.betaEndpoint;
                   } );
      }//if( db )
    }//if( verdict allows candidates )

#if( PERFORM_DEVELOPER_CHECKS )
    // Consistency of the emitted verdict with the stored metrics.
    {
      assert( status.modeEnergy.has_value() && status.terminationEnergy.has_value() );
      assert( *status.modeEnergy <= *status.terminationEnergy );
      assert( status.continuumElevationSigma >= sk_beta_min_elevation_sigma );
      // A BremLike verdict must be fully consistent with the reported metrics; the Ambiguous /
      // NotBremLike distinction additionally depends on range-internal quantities (cliff
      // meaningfulness, which of the two evaluation ranges won), so it isnt re-derivable from
      // the status alone.
      if( (status.verdict == BetaContinuumVerdict::BremLike)
          || (status.verdict == BetaContinuumVerdict::AnnihilationDominated) )
      {
        assert( status.elevatedGammaPeaks.empty() );
        assert( *status.terminationRatio < sk_beta_term_ratio_max );
        assert( *status.terminationEnergy < sk_beta_term_max_kev );
        assert( status.characteristicEnergy.has_value()
                && (*status.characteristicEnergy < sk_beta_echar_max_kev) );
      }
      assert( (status.verdict != BetaContinuumVerdict::AnnihilationDominated)
              || status.annihilationDominant );
      assert( (status.verdict != BetaContinuumVerdict::BackgroundNotLoaded)
              && (status.verdict != BetaContinuumVerdict::NoContinuumElevation) );
    }
#endif //PERFORM_DEVELOPER_CHECKS

    return status;
  }//beta_continuum_check( const BetaContinuumCheckInput &input )


  BetaContinuumCheckStatus beta_continuum_check( const BetaContinuumCheckOptions &options,
                                                 InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "beta_continuum_check: No InterSpec session available" );

    if( options.specType == SpecUtils::SpectrumType::Background )
      throw runtime_error( "beta_continuum_check: the Background spectrum cannot be the candidate"
                           " spectrum - it is the subtraction baseline" );

    shared_ptr<SpecMeas> meas = interspec->measurment( options.specType );
    if( !meas )
      throw runtime_error( "beta_continuum_check: No measurement loaded for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );

    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram( options.specType );
    if( !spectrum )
      throw runtime_error( "beta_continuum_check: No spectrum displayed for "
                          + string(SpecUtils::descriptionText(options.specType))
                          + " spectrum" );

    const set<int> sample_nums = interspec->displayedSamples( options.specType );

    BetaContinuumCheckInput input;
    input.foreground = spectrum;

    // Without a displayed background we cannot analyze - the natural background continuum is
    // itself a smooth falling curve, so a no-subtraction shape analysis has an unacceptable
    // false-positive rate (empirically verified) - return the clear status instead of throwing.
    shared_ptr<SpecMeas> back_meas = interspec->measurment( SpecUtils::SpectrumType::Background );
    shared_ptr<const SpecUtils::Measurement> back_spectrum
        = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
    if( !back_meas || !back_spectrum )
    {
      BetaContinuumCheckStatus status;
      status.verdict = BetaContinuumVerdict::BackgroundNotLoaded;
      status.warnings.push_back( "A background spectrum is required to test for a beta-emitter"
                                 " (bremsstrahlung) continuum; please load a background measurement"
                                 " taken with the same detector." );
      return status;
    }

    // Foreground user + auto-search peaks, de-duplicated - same set detected_peaks reports.
    DetectedPeaksOptions fg_opts;
    fg_opts.specType = options.specType;
    fg_opts.nonBackgroundPeaksOnly = false;
    fg_opts.lowerEnergy = std::nullopt;
    fg_opts.upperEnergy = std::nullopt;
    const DetectedPeakStatus fg_peaks = detected_peaks( fg_opts, interspec );
    input.foregroundPeaks = fg_peaks.peaks;

    // Background peaks: recover NORM peaks hidden under foreground peaks first, then fetch the
    // (augmented) auto-search peaks - the same sequence detected_peaks uses; we are on the GUI
    // thread, where blocking on the searches is acceptable (they are cached/coalesced).
    const set<int> back_samples = interspec->displayedSamples( SpecUtils::SpectrumType::Background );
    const shared_ptr<const PeakFitDetPrefs> fitPrefs = meas ? meas->peakFitDetPrefs() : nullptr;
    PeakSearchGuiUtils::ensure_background_peaks_recovered( meas, sample_nums, spectrum,
                                                          back_meas, back_samples,
                                                          back_spectrum, fitPrefs );
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> back_auto_peaks
        = PeakSearchGuiUtils::get_or_launch_automated_search_peaks( back_meas, back_samples,
                                back_spectrum, back_meas->detector(), fitPrefs ).get();
    if( back_auto_peaks )
      input.backgroundPeaks.insert( end(input.backgroundPeaks),
                                    begin(*back_auto_peaks), end(*back_auto_peaks) );
    const shared_ptr<const deque<shared_ptr<const PeakDef>>> back_user_peaks = back_meas->peaks( back_samples );
    if( back_user_peaks )
      input.backgroundPeaks.insert( end(input.backgroundPeaks),
                                    begin(*back_user_peaks), end(*back_user_peaks) );

    input.background = back_spectrum;
    input.expectedFwhm = make_expected_fwhm_fcn( interspec );

    return beta_continuum_check( input );
  }//beta_continuum_check( const BetaContinuumCheckOptions &options, InterSpec *interspec )

} // namespace AnalystChecks
