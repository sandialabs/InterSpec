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
#include <string>
#include <limits>
#include <algorithm>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AnalystChecks.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/IsotopeId.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/FitPeaksForNuclides.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/RelActCalcAuto.h"

#include "SpecUtils/SpecUtilsAsync.h"

#include "SandiaDecay/SandiaDecay.h"

#include <Wt/WApplication>

using namespace std;

namespace AnalystChecks
{
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
    
    shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = meas->automatedSearchPeaks(sample_nums);
    shared_ptr<const deque<shared_ptr<const PeakDef>>> user_peaks = meas->peaks(sample_nums);
    
    if (!auto_peaks) {
      // Search for peaks
      const bool singleThreaded = false;
      shared_ptr<const PeakFitDetPrefs> fitPrefsAuto = meas->peakFitDetPrefs();
      const auto det = meas->detector();
      const vector<shared_ptr<const PeakDef>> found_auto_peaks
                = ExperimentalAutomatedPeakSearch::search_for_peaks(spectrum, det, user_peaks, singleThreaded, fitPrefsAuto);
        
      auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>(begin(found_auto_peaks), end(found_auto_peaks));
      auto_peaks = autopeaksdeque;
      meas->setAutomatedSearchPeaks(sample_nums, autopeaksdeque);
    }

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

      // Get or search for background auto peaks
      shared_ptr<const deque<shared_ptr<const PeakDef>>> background_auto_peaks = background_meas->automatedSearchPeaks(background_sample_nums);

      if( !background_auto_peaks && background_spectrum )
      {
        const bool singleThreaded = false;
        shared_ptr<const PeakFitDetPrefs> bgFitPrefs = background_meas->peakFitDetPrefs();
        const auto det = background_meas->detector();

        shared_ptr<const deque<shared_ptr<const PeakDef>>> background_user_peaks;
        if( background_meas->sampleNumsWithPeaks().count(background_sample_nums) )
          background_user_peaks = background_meas->peaks(background_sample_nums);

        const vector<shared_ptr<const PeakDef>> found_background_auto_peaks
                    = ExperimentalAutomatedPeakSearch::search_for_peaks(background_spectrum, det, background_user_peaks, singleThreaded, bgFitPrefs);

        auto background_autopeaksdeque = make_shared<deque<shared_ptr<const PeakDef>>>(begin(found_background_auto_peaks), end(found_background_auto_peaks));
        background_auto_peaks = background_autopeaksdeque;
        background_meas->setAutomatedSearchPeaks(background_sample_nums, background_autopeaksdeque);
      }

      // Filter foreground peaks based on background
      vector<shared_ptr<const PeakDef>> filtered_peaks;
      vector<shared_ptr<const PeakDef>> filtered_analysis_peaks;
      const double foreground_live_time = (spectrum->live_time() > 0.0f) ? spectrum->live_time() : 1.0f;
      const double background_live_time = (background_spectrum->live_time() > 0.0f ? background_spectrum->live_time() : 1.0f);

      for( const shared_ptr<const PeakDef> &fg_peak : all_peaks )
      {
        bool include_peak = true;

        // Check if there is a corresponding background peak
        if( background_auto_peaks )
        {
          // Find the closest matching background peak
          shared_ptr<const PeakDef> closest_bg_peak;
          double smallest_energy_diff = std::numeric_limits<double>::infinity();

          for( const shared_ptr<const PeakDef> &bg_peak : *background_auto_peaks )
          {
            assert( fg_peak != bg_peak );
            
            // Check if peaks are at roughly the same energy (within their widths)
            const double energy_diff = fabs(fg_peak->mean() - bg_peak->mean());
            const double avg_fwhm = 0.75 * (fg_peak->fwhm() + bg_peak->fwhm()) / 2.0;

            if( (energy_diff < avg_fwhm) && (energy_diff < smallest_energy_diff) )
            {
              closest_bg_peak = bg_peak;
              smallest_energy_diff = energy_diff;
            }
          }

          // If we found a matching background peak, compare CPS
          if( closest_bg_peak )
          {
            const double fg_cps = (fg_peak->amplitude() / foreground_live_time);
            const double bg_cps = (closest_bg_peak->amplitude() / background_live_time);

            // Check if foreground peak is elevated by at least 20% relative to background
            const double elevation_threshold = 1.20;
            if( fg_cps > (bg_cps * elevation_threshold) )
            {
              include_peak = true;

              // Foreground is elevated by >20%, now check statistical significance if uncertainties are available
              const double fg_amp_uncert = fg_peak->amplitudeUncert();
              const double bg_amp_uncert = closest_bg_peak->amplitudeUncert();

              if( (fg_amp_uncert > 0.0) && (bg_amp_uncert > 0.0) )
              {
                // Calculate CPS uncertainties
                const double fg_cps_uncert = fg_amp_uncert / foreground_live_time;
                const double bg_cps_uncert = bg_amp_uncert / background_live_time;

                // Calculate combined uncertainty
                const double combined_uncert = sqrt(fg_cps_uncert*fg_cps_uncert + bg_cps_uncert*bg_cps_uncert);

                // Check if foreground is elevated by at least 2.25 sigma
                const double min_sigma = 2.25;
                const double sigma_elevation = ((fg_cps - bg_cps) / combined_uncert);
                include_peak = (sigma_elevation > min_sigma);
              }
            }else
            {
              // Foreground peak is not elevated by 20% - exclude it
              include_peak = false;
            }//if( fg_cps > (bg_cps * elevation_threshold) ) / else
          }else
          {
            // If we didnt detect the peak in the background, but this is a pretty insignificant peak that
            //  could be a background peak, then dont include it
            if( (fg_peak->amplitude() < 100) && (fg_peak->amplitudeUncert() > 0.0) )
            {
              // TODO: check if peak potentually lines up with a background peak, and if so, estimate the amplitude we would expect for it in the background spectrum, and if we wouldnt totally expect to detect it in the background spectrum, and the foreground peak significance isnt that high or its amplitude is consistent with what is expected from other confirmed background peaks of the same series, then dont include it.
            }
          }//if( closest_bg_peak )
        }//if( background_auto_peaks )

        if( include_peak )
        {
          filtered_peaks.push_back(fg_peak);

          // Check if this peak is also in the user_analysis_peaks
          const bool is_user_peak = std::find(user_analysis_peaks.begin(), user_analysis_peaks.end(), fg_peak) != user_analysis_peaks.end();
          if( is_user_peak )
            filtered_analysis_peaks.push_back(fg_peak);
        }
      }//for( const shared_ptr<const PeakDef> &fg_peak : all_peaks )

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
  
  
  FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec )
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
    foundPeaks = searchForPeakFromUser( options.energy, pixelPerKev, data, origPeaks, det, auto_search_peaks, fitPrefs );

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

      string error_msg = "Could not fit peak at " + SpecUtils::printCompact(options.energy, 4) + " keV";
      if( within_existing_roi )
        error_msg += " - the energy is within the ROI of an existing peak, so maybe the peak you wanted already exists in the analysis peak list, or another peak just cant be fit in the same ROI";
      else
        error_msg += " - please try a different energy";
      error_msg += ".";

      throw runtime_error( error_msg );
    }
    
    // The new peak may share a ROI with previously existing peaks, so we have to figure out
    //  what peak we just fit.
    //  TODO: For the moment, we'll just match up to nearest peak, but in a normal user double click for a peak, we match up peaks in ROI to previous peaks, and then take whats left-over, which is maybe a better way.
    shared_ptr<const PeakDef> fit_peak;
    for( shared_ptr<const PeakDef> peak : peaks_to_add_in )
    {
      if( peak->parentNuclide() || peak->xrayElement() || peak->reaction() )
        continue;  //We know new peak doesnt have a source assigned yet
      
      const double diff = fabs(peak->mean() - options.energy);
      if( !fit_peak || (diff < fabs(fit_peak->mean() - options.energy)) )
        fit_peak = peak;
    }//for( shared_ptr<const PeakDef> peak : peaks_to_add_in )
    
    
    assert( fit_peak );
    if( !fit_peak )
      throw runtime_error( "Could not identify newly fit peak." );
    
    // get fit peak
    const string source = options.source.has_value() ? options.source.value() : ""s;

    if( source.empty() )
    {
      // Nothing to do here
    }else if( SpecUtils::icontains(source, "unknown") || SpecUtils::istarts_with(source, "unk") )
    {
      // `ToolRegistry::executeGetUnidentifiedDetectedPeaks(...)` will treat the peak as unknown if user label is set to "UNKNOWN"
      shared_ptr<PeakDef> new_fit_peak = make_shared<PeakDef>( *fit_peak );
      new_fit_peak->setUserLabel( "UNKNOWN" );
      fit_peak = new_fit_peak;
    }else
    {
      shared_ptr<PeakDef> new_fit_peak = make_shared<PeakDef>( *fit_peak );
      
      PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( *new_fit_peak, source, 4.0 );
      
      if( res == PeakModel::SetGammaSource::FailedSourceChange )
        throw runtime_error( "Peak was fit, but source string '" + source + "' was not valid to set source." );
      
      bool replaced_peak = false;
      for( size_t index = 0; index < peaks_to_add_in.size(); ++index )
      {
        if( peaks_to_add_in[index] == fit_peak )
        {
          replaced_peak = true;
          fit_peak = new_fit_peak;
          peaks_to_add_in[index] = new_fit_peak;
          break;
        }
      }
      assert( replaced_peak );
      
      vector<Wt::WColor> src_colors;
      vector<shared_ptr<const PeakDef>> other_src_peaks;
      for( const shared_ptr<const PeakDef> &srcPeak : origPeaks )
      {
        if( (srcPeak->parentNuclide() || srcPeak->xrayElement() || srcPeak->reaction())
            && srcPeak->parentNuclide()==new_fit_peak->parentNuclide()
            && srcPeak->xrayElement()==new_fit_peak->xrayElement()
            && srcPeak->reaction()==new_fit_peak->reaction() )
        {
          other_src_peaks.push_back( srcPeak );
          
          if( !srcPeak->lineColor().isDefault()
             && (std::find( begin(src_colors), end(src_colors),srcPeak->lineColor()) == end(src_colors)) )
          {
            src_colors.push_back( srcPeak->lineColor() );
          }
        }
      }//for( const auto &p : *m_peaks )
          
      if( src_colors.size() == 1 )
      {
        new_fit_peak->setLineColor( src_colors[0] );
      }else if( other_src_peaks.empty() )
      {
        ReferencePhotopeakDisplay *ref_lines = interspec->referenceLinesWidget();
        if( ref_lines )
        {
          string src_name;
          if( new_fit_peak->parentNuclide() )
            src_name = new_fit_peak->parentNuclide()->symbol;
          else if( new_fit_peak->xrayElement() )
            src_name = new_fit_peak->xrayElement()->symbol;
          else if( new_fit_peak->reaction() )
            src_name = new_fit_peak->reaction()->name();
          
          Wt::WColor c = ref_lines->suggestColorForSource( src_name );
          
          if( c.isDefault() && !src_name.empty() )
          {
            c = ref_lines->nextGenericSourceColor();
            if( !c.isDefault() )
              ref_lines->updateColorCacheForSource( src_name, c );
          }
          
          new_fit_peak->setLineColor( c );
        }
      }
    }//if( options.source.has_value() )
    
    PeakModel *pmodel = interspec->peakModel();
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
    float sum_yields_sqrt_energy = 0.0f;
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      sum_yields_sqrt_energy += std::get<1>(gamma) * sqrt(std::get<0>(gamma));
    }
    
    vector<tuple<float,float,float>> answer;
    answer.reserve(gamma_energies_and_yields.size());
    for( const tuple<float,float> &gamma : gamma_energies_and_yields )
    {
      float importance = std::get<1>(gamma) * sqrt(std::get<0>(gamma)) / sum_yields_sqrt_energy;
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
      const vector<SandiaDecay::EnergyRatePair> photons
                                        = mixture.photons(age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy);

      vector<tuple<float,float>> gamma_energies_and_yields;
      for( const auto &photon : photons )
        gamma_energies_and_yields.push_back({photon.energy, photon.numPerSecond});

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

      if( sources.empty() )
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

      // TODO: we will need to update `config` from default in the future

      const FitPeaksForNuclides::PeakFitResult fit_results = FitPeaksForNuclides::fit_peaks_for_nuclides(
        auto_search_peaks, target_spectrum, sources, user_peaks, background, drf, fit_options, fit_config, peak_fit_prefs
      );

      switch( fit_results.status )
      {
        case RelActCalcAuto::RelActAutoSolution::Status::Success:
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

      ReferencePhotopeakDisplay * const ref_lines = interspec ? interspec->referenceLinesWidget() : nullptr;

      result.warnings = fit_results.warnings;
      for( const PeakDef &p : fit_results.observable_peaks )
      {
        auto sp = make_shared<PeakDef>(p);

        string src_str;
        if( sp->parentNuclide() )
          src_str = sp->parentNuclide()->symbol;
        else if( sp->xrayElement() )
          src_str = sp->xrayElement()->symbol;
        else if( sp->reaction() )
          src_str = sp->reaction()->name();

        if( ref_lines && !src_str.empty() )
        {
          Wt::WColor color = ref_lines->suggestColorForSource( src_str );
          if( color.isDefault() )
          {
            color = ref_lines->nextGenericSourceColor();
            ref_lines->updateColorCacheForSource( src_str, color );
          }
          sp->setLineColor( color );
        }//if( ref_lines && !src_str.empty() )

        result.fitPeaks.push_back( sp );
      }//for( const PeakDef &p : fit_results.observable_peaks )

      result.peaksToRemove = fit_results.original_peaks_to_remove;

      // If doNotAddPeaksToUserSession is false, update the peak model using setPeaks
      //  to atomically replace the affected peaks while preserving unrelated ones.
      if( !options.doNotAddPeaksToUserSession && !result.fitPeaks.empty() )
      {
        PeakModel * const pmodel = interspec->peakModel();
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

  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_characteristics_near_energy( const double energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");
    
    shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !meas || !spectrum )
      throw std::runtime_error("No foreground loaded");
    
    
    float fwhm = -1.0;
    
    // Try the currently loaded peak detector response
    if( meas->detector() && meas->detector()->hasResolutionInfo() )
      fwhm = meas->detector()->peakResolutionFWHM( static_cast<float>(energy) );
    
    // Estimate the FWHM response...
    shared_ptr<deque<shared_ptr<const PeakDef>>> peak_deque;
    if( fwhm <= 0.0 )
    {
      try
      {
        DetectedPeaksOptions opts;
        opts.specType = SpecUtils::SpectrumType::Foreground;
        const DetectedPeakStatus peaks = detected_peaks( opts, interspec );
        peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks.peaks), end(peaks.peaks) );
        
        DetectorPeakResponse drf;
        drf.setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                               0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic);
        
        drf.fitResolution( peak_deque, spectrum, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
        
        fwhm = drf.peakResolutionFWHM( static_cast<float>(energy) );
      }catch( std::exception & )
      {
      }
    }//if( fwhm <= 0.0 )
    
    // Find nearest peak, and if within 20% of the energy, use it.
    if( (fwhm <= 0.0) && peak_deque )
    {
      double smallest_energy_diff = std::numeric_limits<double>::max();
      double nearest_fwhm = 0.0;
      for( const auto &peak : *peak_deque )
      {
        double energy_diff = std::abs(peak->mean() - energy);
        if( energy_diff < smallest_energy_diff )
        {
          smallest_energy_diff = energy_diff;
          nearest_fwhm = peak->fwhm();
        }
      }

      if( smallest_energy_diff < 0.2*energy )
        fwhm = nearest_fwhm;
    }//if( (fwhm <= 0.0) && peak_deque )
    
    // Finally, a wag based on detection type.
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
    const vector<float> pars = isHPGe ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
    fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars );

    if( fwhm <= 0.0 )
      throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy) );
    
    return get_nuclides_with_characteristics_in_energy_range( energy - fwhm, energy + fwhm, interspec );
  }//get_characteristics_near_energy(...)
  
  
  float get_expected_fwhm( const double energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");
    
    shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
    shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !meas || !spectrum )
      throw std::runtime_error("No foreground loaded");
    
    float fwhm = -1.0;
    
    // Try the currently loaded peak detector response
    if( meas->detector() && meas->detector()->hasResolutionInfo() )
      fwhm = meas->detector()->peakResolutionFWHM( static_cast<float>(energy) );
    
    // Estimate the FWHM response using detected peaks...
    shared_ptr<deque<shared_ptr<const PeakDef>>> peak_deque;
    if( fwhm <= 0.0 )
    {
      try
      {
        DetectedPeaksOptions opts;
        opts.specType = SpecUtils::SpectrumType::Foreground;
        const DetectedPeakStatus peaks = detected_peaks( opts, interspec );
        peak_deque = make_shared<deque<shared_ptr<const PeakDef>>>( begin(peaks.peaks), end(peaks.peaks) );
        
        DetectorPeakResponse drf;
        drf.setIntrinsicEfficiencyFormula( "1.0", 2.54*PhysicalUnits::cm, PhysicalUnits::keV,
                                               0.0f, 0.0f, DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic);
        
        drf.fitResolution( peak_deque, spectrum, DetectorPeakResponse::ResolutionFnctForm::kSqrtPolynomial );
        
        fwhm = drf.peakResolutionFWHM( static_cast<float>(energy) );
      }catch( std::exception & )
      {
      }
    }//if( fwhm <= 0.0 )
    
    // Find nearest peak, and if within 20% of the energy, use it.
    if( (fwhm <= 0.0) && peak_deque )
    {
      double smallest_energy_diff = std::numeric_limits<double>::max();
      double nearest_fwhm = 0.0;
      for( const auto &peak : *peak_deque )
      {
        double energy_diff = std::abs(peak->mean() - energy);
        if( energy_diff < smallest_energy_diff )
        {
          smallest_energy_diff = energy_diff;
          nearest_fwhm = peak->fwhm();
        }
      }

      if( smallest_energy_diff < 0.2*energy )
        fwhm = nearest_fwhm;
    }//if( (fwhm <= 0.0) && peak_deque )
    
    // Finally, a wag based on detection type.
    if( fwhm <= 0.0 )
    {
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
      const vector<float> pars = isHPGe ? vector<float>{ 1.54f, 0.264f, 0.33f } : vector<float>{ -6.5f, 7.5f, 0.55f };
      fwhm = DetectorPeakResponse::peakResolutionFWHM( energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars );
    }

    if( fwhm <= 0.0 )
      throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy) );
    
    return fwhm;
  }//get_expected_fwhm(...)
  
  
  SpectrumCountsInEnergyRange get_counts_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec )
  {
    if( !interspec )
      throw std::runtime_error("No InterSpec session available");

    shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
    if( !foreground )
      throw runtime_error( "get_counts_in_energy_range: no foreground loaded." );

    if( lower_energy > upper_energy )
      std::swap( lower_energy, upper_energy );

    SpectrumCountsInEnergyRange answer;
    answer.lower_energy = lower_energy;
    answer.upper_energy = upper_energy;

    answer.foreground_counts = foreground->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
    if( foreground->live_time() > 0.0 )
      answer.foreground_cps = answer.foreground_counts / foreground->live_time();
    else 
      answer.foreground_cps = std::numeric_limits<double>::quiet_NaN();
    
    shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);
    if( background )
    {
      SpectrumCountsInEnergyRange::CountsWithComparisonToForeground back;

      back.counts = background->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
      if( background->live_time() > 0.0 )
        back.cps = back.counts / background->live_time();
      else
        back.cps = std::numeric_limits<double>::quiet_NaN();

      const double backSF = foreground->live_time() / background->live_time();
      const double nfore = answer.foreground_counts;
      const double nback = back.counts;
      const double scaleback = nback * backSF;
      const double backsigma = sqrt(nback);
      const double forsigma = sqrt(nfore);
      const double backscalesigma = backSF * backsigma;
      const double total_back_fore_sigma = sqrt(backscalesigma*backscalesigma + forsigma*forsigma);
      const double nsigma = (nfore - scaleback) / total_back_fore_sigma;

      back.num_sigma_rel_foreground = nsigma;

      answer.background_info = std::move( back );
    }//if( background )

    shared_ptr<const SpecUtils::Measurement> secondary = interspec->displayedHistogram(SpecUtils::SpectrumType::SecondForeground);
    if( secondary )
    {
      SpectrumCountsInEnergyRange::CountsWithComparisonToForeground sec;

      sec.counts = secondary->gamma_integral( static_cast<float>(lower_energy), static_cast<float>(upper_energy) );
      if( secondary->live_time() > 0.0 )
        sec.cps = sec.counts / secondary->live_time();
      else
        sec.cps = std::numeric_limits<double>::quiet_NaN();

      const double secSF = foreground->live_time() / secondary->live_time();
      const double nfore = answer.foreground_counts;
      const double nsec = sec.counts;
      const double scalesec = nsec * secSF;
      const double secsigma = sqrt(nsec);
      const double forsigma = sqrt(nfore);
      const double secscalesigma = secSF * secsigma;
      const double total_sec_fore_sigma = sqrt(secscalesigma*secscalesigma + forsigma*forsigma);
      const double nsigma = (nfore - scalesec) / total_sec_fore_sigma;

      sec.num_sigma_rel_foreground = nsigma;

      answer.secondary_info = std::move( sec );
    }//if( secondary )

    return answer;
  }//SpectrumCountsInEnergyRange get_counts_in_energy_range( const double lower_energy, const double upper_energy, InterSpec *interspec )


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
      case SumPeakType::Unknown:     return "Unknown";
    }//switch( type )

    throw runtime_error( "Invalid SumPeakType value" );
  }//const char* to_string( SumPeakType type )


  const char* to_string( EditPeakAction action )
  {
    switch( action )
    {
      case EditPeakAction::SetEnergy:                    return "SetEnergy";
      case EditPeakAction::SetFwhm:                      return "SetFwhm";
      case EditPeakAction::SetAmplitude:                 return "SetAmplitude";
      case EditPeakAction::SetEnergyUncertainty:         return "SetEnergyUncertainty";
      case EditPeakAction::SetFwhmUncertainty:           return "SetFwhmUncertainty";
      case EditPeakAction::SetAmplitudeUncertainty:      return "SetAmplitudeUncertainty";
      case EditPeakAction::SetRoiLower:                  return "SetRoiLower";
      case EditPeakAction::SetRoiUpper:                  return "SetRoiUpper";
      case EditPeakAction::SetSkewType:                  return "SetSkewType";
      case EditPeakAction::SetContinuumType:             return "SetContinuumType";
      case EditPeakAction::SetSource:                    return "SetSource";
      case EditPeakAction::SetColor:                     return "SetColor";
      case EditPeakAction::SetUserLabel:                 return "SetUserLabel";
      case EditPeakAction::SetUseForEnergyCalibration:   return "SetUseForEnergyCalibration";
      case EditPeakAction::SetUseForShieldingSourceFit:  return "SetUseForShieldingSourceFit";
      case EditPeakAction::SetUseForManualRelEff:        return "SetUseForManualRelEff";
      case EditPeakAction::DeletePeak:                   return "DeletePeak";
      case EditPeakAction::SplitFromRoi:                 return "SplitFromRoi";
      case EditPeakAction::MergeWithLeft:                return "MergeWithLeft";
      case EditPeakAction::MergeWithRight:               return "MergeWithRight";
    }//switch( action )

    throw runtime_error( "Invalid EditPeakAction value" );
  }//const char* to_string( EditPeakAction action )


  EditPeakAction edit_peak_action_from_string( const std::string &str )
  {
    if( str == "SetEnergy" )                   return EditPeakAction::SetEnergy;
    if( str == "SetFwhm" )                     return EditPeakAction::SetFwhm;
    if( str == "SetAmplitude" )                return EditPeakAction::SetAmplitude;
    if( str == "SetEnergyUncertainty" )        return EditPeakAction::SetEnergyUncertainty;
    if( str == "SetFwhmUncertainty" )          return EditPeakAction::SetFwhmUncertainty;
    if( str == "SetAmplitudeUncertainty" )     return EditPeakAction::SetAmplitudeUncertainty;
    if( str == "SetRoiLower" )                 return EditPeakAction::SetRoiLower;
    if( str == "SetRoiUpper" )                 return EditPeakAction::SetRoiUpper;
    if( str == "SetSkewType" )                 return EditPeakAction::SetSkewType;
    if( str == "SetContinuumType" )            return EditPeakAction::SetContinuumType;
    if( str == "SetSource" )                   return EditPeakAction::SetSource;
    if( str == "SetColor" )                    return EditPeakAction::SetColor;
    if( str == "SetUserLabel" )                return EditPeakAction::SetUserLabel;
    if( str == "SetUseForEnergyCalibration" )  return EditPeakAction::SetUseForEnergyCalibration;
    if( str == "SetUseForShieldingSourceFit" ) return EditPeakAction::SetUseForShieldingSourceFit;
    if( str == "SetUseForManualRelEff" )       return EditPeakAction::SetUseForManualRelEff;
    if( str == "DeletePeak" )                  return EditPeakAction::DeletePeak;
    if( str == "SplitFromRoi" )                return EditPeakAction::SplitFromRoi;
    if( str == "MergeWithLeft" )               return EditPeakAction::MergeWithLeft;
    if( str == "MergeWithRight" )              return EditPeakAction::MergeWithRight;

    throw runtime_error( "Invalid edit action: " + str );
  }//EditPeakAction edit_peak_action_from_string( const std::string &str )


  EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec )
  {
    if( !interspec )
      throw runtime_error( "edit_analysis_peak: No InterSpec session available" );

    // Get the PeakModel for the specified spectrum type
    PeakModel *peakModel = interspec->peakModel();
    if( !peakModel )
      throw runtime_error( "edit_analysis_peak: No peak model available" );

    // Find nearest peak
    PeakModel::PeakShrdPtr peak = peakModel->nearestPeak( options.energy );
    if( !peak )
    {
      return EditAnalysisPeakStatus{
        false,
        "No peak found near " + std::to_string(options.energy) + " keV",
        std::nullopt,
        {}
      };
    }//if( !peak )

    // Validate peak is within 1 FWHM of requested energy
    const double fwhm = get_expected_fwhm( options.energy, interspec );
    const double energy_diff = std::fabs( peak->mean() - options.energy );
    if( energy_diff > fwhm )
    {
      return EditAnalysisPeakStatus{
        false,
        "Nearest peak at " + std::to_string(peak->mean()) + " keV is more than 1 FWHM ("
          + std::to_string(fwhm) + " keV) away from requested " + std::to_string(options.energy) + " keV",
        std::nullopt,
        {}
      };
    }//if( energy_diff > fwhm )

    // Handle delete action separately since it doesn't need a modified peak
    if( options.editAction == EditPeakAction::DeletePeak )
    {
      peakModel->removePeak( peak );
      return EditAnalysisPeakStatus{
        true,
        "Deleted peak at " + std::to_string(peak->mean()) + " keV",
        std::nullopt,
        {}
      };
    }//if( DeletePeak )

    // Create a modified copy of the peak.
    //  Note: PeakDef copy shares the same PeakContinuum pointer, so we must create a new
    //  independent continuum.  If other peaks share the original continuum (i.e., are in the
    //  same ROI), we need to copy those peaks too and point them to the new continuum, so the
    //  ROI remains intact but fully independent from the originals.
    PeakDef modifiedPeak = *peak;

    // Collect all peaks sharing the same ROI continuum as the original peak
    const std::shared_ptr<const PeakContinuum> originalContinuum = peak->continuum();
    vector<shared_ptr<const PeakDef>> roiSiblings; // other peaks in same ROI (not including the target peak)
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
    }//end collect ROI siblings

    // Create a new independent continuum for the modified peak
    modifiedPeak.makeUniqueNewContinuum();

    // Perform the edit based on action
    try
    {
      switch( options.editAction )
      {
        case EditPeakAction::SetEnergy:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetEnergy requires doubleValue parameter" );
          modifiedPeak.setMean( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetFwhm:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetFwhm requires doubleValue parameter" );
          // Convert FWHM to sigma (sigma = FWHM / 2.35482)
          modifiedPeak.setSigma( (*options.doubleValue) / 2.35482 );
          break;
        }

        case EditPeakAction::SetAmplitude:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetAmplitude requires doubleValue parameter" );
          modifiedPeak.setAmplitude( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetEnergyUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetEnergyUncertainty requires doubleValue parameter" );
          modifiedPeak.setMeanUncert( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetFwhmUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetFwhmUncertainty requires doubleValue parameter" );
          // Convert FWHM uncertainty to sigma uncertainty
          modifiedPeak.setSigmaUncert( (*options.doubleValue) / 2.35482 );
          break;
        }

        case EditPeakAction::SetAmplitudeUncertainty:
        {
          if( !options.doubleValue )
            throw runtime_error( "SetAmplitudeUncertainty requires doubleValue parameter" );
          modifiedPeak.setAmplitudeUncert( *options.doubleValue );
          break;
        }

        case EditPeakAction::SetRoiLower:
        case EditPeakAction::SetRoiUpper:
        {
          if( !options.doubleValue )
            throw runtime_error( string(to_string(options.editAction)) + " requires doubleValue parameter" );

          shared_ptr<PeakContinuum> continuum = modifiedPeak.getContinuum();
          const double current_lower = continuum->lowerEnergy();
          const double current_upper = continuum->upperEnergy();

          const double new_lower = (options.editAction == EditPeakAction::SetRoiLower)
                                    ? *options.doubleValue : current_lower;
          const double new_upper = (options.editAction == EditPeakAction::SetRoiUpper)
                                    ? *options.doubleValue : current_upper;

          if( new_lower >= new_upper )
            throw runtime_error( "ROI lower energy must be less than upper energy" );

          continuum->setRange( new_lower, new_upper );
          break;
        }

        case EditPeakAction::SetSkewType:
        {
          if( !options.stringValue )
            throw runtime_error( "SetSkewType requires stringValue parameter" );

          // Parse skew type string
          PeakDef::SkewType skew_type;
          const string &skew_str = *options.stringValue;
          if( skew_str == "NoSkew" )                      skew_type = PeakDef::SkewType::NoSkew;
          else if( skew_str == "Bortel" )                 skew_type = PeakDef::SkewType::Bortel;
          else if( skew_str == "GaussExp" )               skew_type = PeakDef::SkewType::GaussExp;
          else if( skew_str == "CrystalBall" )            skew_type = PeakDef::SkewType::CrystalBall;
          else if( skew_str == "ExpGaussExp" )            skew_type = PeakDef::SkewType::ExpGaussExp;
          else if( skew_str == "DoubleSidedCrystalBall" ) skew_type = PeakDef::SkewType::DoubleSidedCrystalBall;
          else throw runtime_error( "Invalid skew type: " + skew_str );

          modifiedPeak.setSkewType( skew_type );
          break;
        }

        case EditPeakAction::SetContinuumType:
        {
          if( !options.stringValue )
            throw runtime_error( "SetContinuumType requires stringValue parameter" );

          // Parse continuum type string
          PeakContinuum::OffsetType continuum_type;
          const string &cont_str = *options.stringValue;
          if( cont_str == "None" )           continuum_type = PeakContinuum::OffsetType::NoOffset;
          else if( cont_str == "Constant" )  continuum_type = PeakContinuum::OffsetType::Constant;
          else if( cont_str == "Linear" )    continuum_type = PeakContinuum::OffsetType::Linear;
          else if( cont_str == "Quadratic" ) continuum_type = PeakContinuum::OffsetType::Quadratic;
          else if( cont_str == "Cubic" )     continuum_type = PeakContinuum::OffsetType::Cubic;
          else if( cont_str == "FlatStep" )  continuum_type = PeakContinuum::OffsetType::FlatStep;
          else if( cont_str == "LinearStep" )continuum_type = PeakContinuum::OffsetType::LinearStep;
          else if( cont_str == "BiLinearStep" )continuum_type = PeakContinuum::OffsetType::BiLinearStep;
          else if( cont_str == "FlatStepCDF" )continuum_type = PeakContinuum::OffsetType::FlatStepCDF;
          else if( cont_str == "LinearStepCDF" )continuum_type = PeakContinuum::OffsetType::LinearStepCDF;
          else if( cont_str == "BiLinearStepCDF" )continuum_type = PeakContinuum::OffsetType::BiLinearStepCDF;
          else if( cont_str == "External" )  continuum_type = PeakContinuum::OffsetType::External;
          else throw runtime_error( "Invalid continuum type: " + cont_str );

          shared_ptr<PeakContinuum> continuum = modifiedPeak.getContinuum();
          continuum->setType( continuum_type );
          break;
        }

        case EditPeakAction::SetSource:
        {
          if( !options.stringValue )
            throw runtime_error( "SetSource requires stringValue parameter" );

          // Use PeakModel::setNuclideXrayReaction to parse and set the source
          PeakModel::SetGammaSource result = PeakModel::setNuclideXrayReaction( modifiedPeak, *options.stringValue, 4.0 );
          if( result == PeakModel::SetGammaSource::FailedSourceChange )
            throw runtime_error( "Failed to set source: " + *options.stringValue );
          break;
        }

        case EditPeakAction::SetColor:
        {
          if( !options.stringValue )
            throw runtime_error( "SetColor requires stringValue parameter" );

          // Parse CSS color string to Wt::WColor
          Wt::WColor color( *options.stringValue );
          modifiedPeak.setLineColor( color );
          break;
        }

        case EditPeakAction::SetUserLabel:
        {
          if( !options.stringValue )
            throw runtime_error( "SetUserLabel requires stringValue parameter" );

          modifiedPeak.setUserLabel( *options.stringValue );
          break;
        }

        case EditPeakAction::SetUseForEnergyCalibration:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForEnergyCalibration requires boolValue parameter" );

          modifiedPeak.useForEnergyCalibration( *options.boolValue );
          break;
        }

        case EditPeakAction::SetUseForShieldingSourceFit:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForShieldingSourceFit requires boolValue parameter" );

          modifiedPeak.useForShieldingSourceFit( *options.boolValue );
          break;
        }

        case EditPeakAction::SetUseForManualRelEff:
        {
          if( !options.boolValue )
            throw runtime_error( "SetUseForManualRelEff requires boolValue parameter" );

          modifiedPeak.useForManualRelEff( *options.boolValue );
          break;
        }

        case EditPeakAction::SplitFromRoi:
        {
          // Make the peak have its own unique continuum, separate from the ROI.
          //  We already called makeUniqueNewContinuum() above, so we just need to make sure the
          //  sibling peaks are NOT updated to share the new continuum - they should keep
          //  their original shared continuum.
          roiSiblings.clear();
          break;
        }

        case EditPeakAction::MergeWithLeft:
        case EditPeakAction::MergeWithRight:
        {
          // Merging two ROIs: create a new combined ROI whose energy range is the union of
          //  the two original ROIs, copy all peaks from both ROIs into the new continuum,
          //  then refit the combined ROI to the data.

          shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
          if( !all_peaks_ptr )
            throw runtime_error( "No peaks available in model" );

          const deque<shared_ptr<const PeakDef>> &all_peaks = *all_peaks_ptr;

          // PeakModel::peaks() returns peaks sorted by mean energy
          assert( std::is_sorted( all_peaks.begin(), all_peaks.end(),
            []( const shared_ptr<const PeakDef> &a, const shared_ptr<const PeakDef> &b ){
              return a->mean() < b->mean();
            } ) );

          // Find the adjacent peak to merge with
          const deque<shared_ptr<const PeakDef>>::const_iterator peak_iter
            = find( all_peaks.begin(), all_peaks.end(), peak );
          if( peak_iter == all_peaks.end() )
            throw runtime_error( "Could not find peak in model" );

          shared_ptr<const PeakDef> adjacent_peak;
          if( options.editAction == EditPeakAction::MergeWithLeft )
          {
            if( peak_iter == all_peaks.begin() )
              throw runtime_error( "No peak to the left to merge with" );
            adjacent_peak = *(peak_iter - 1);
          }
          else // MergeWithRight
          {
            if( (peak_iter + 1) == all_peaks.end() )
              throw runtime_error( "No peak to the right to merge with" );
            adjacent_peak = *(peak_iter + 1);
          }

          const shared_ptr<const PeakContinuum> adjContinuum = adjacent_peak->continuum();

          if( adjContinuum == originalContinuum )
            throw runtime_error( "Adjacent peak already shares the same ROI" );

          // Collect all peaks from both ROIs (original and adjacent)
          vector<shared_ptr<const PeakDef>> oldPeaksToRemove;
          for( const shared_ptr<const PeakDef> &p : all_peaks )
          {
            if( (p->continuum() == originalContinuum) || (p->continuum() == adjContinuum) )
              oldPeaksToRemove.push_back( p );
          }

          // Create new continuum whose range is the union of both ROIs
          std::shared_ptr<PeakContinuum> mergedContinuum
            = make_shared<PeakContinuum>( *originalContinuum );
          const double newLower = (std::min)( originalContinuum->lowerEnergy(),
                                            adjContinuum->lowerEnergy() );
          const double newUpper = (std::max)( originalContinuum->upperEnergy(),
                                            adjContinuum->upperEnergy() );
          mergedContinuum->setRange( newLower, newUpper );

          // Copy all peaks from both ROIs, pointing them to the new merged continuum
          vector<PeakDef> mergedPeaks;
          for( const shared_ptr<const PeakDef> &p : oldPeaksToRemove )
          {
            PeakDef peakCopy = *p;
            peakCopy.setContinuum( mergedContinuum );
            mergedPeaks.push_back( peakCopy );
          }

          // Try to refit the combined ROI to the data
          const shared_ptr<const SpecUtils::Measurement> data
            = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
          shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
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
                                        Wt::WFlags<PeakFitLM::PeakFitLMOptions>(0) );

            // Only use refit result if no peaks were lost (e.g., became insignificant)
            if( refitResult.size() == mergedPeaks.size() )
            {
              mergedPeaks.clear();
              for( const shared_ptr<const PeakDef> &p : refitResult )
                mergedPeaks.push_back( *p );
              used_refit = true;
            }
          }//if( data )

          // Replace old peaks from both ROIs with the new merged set
          peakModel->updatePeaks( oldPeaksToRemove, mergedPeaks );

          // Collect the resulting peaks from the model for the return value.
          //  The mergedPeaks all share a common continuum (either mergedContinuum, or the
          //  refit's continuum), so we can identify them by matching that pointer.
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

          // Build the result message with any relevant warnings
          string msg = "Merged ROIs into combined region [" + std::to_string( newLower )
            + ", " + std::to_string( newUpper ) + "] keV with "
            + std::to_string( peaks_in_roi.size() ) + " peaks";

          if( used_refit )
          {
            msg += " (refit to data)";
          }
          else
          {
            msg += ". Warning: the combined ROI could not be refit to the data because"
                   " one or more peaks became insignificant during the fit; the un-fit"
                   " peak parameters were used instead.";
          }

          // Warn if the peak being merged is far from the adjacent peak
          const double avg_fwhm = 0.5 * (peak->fwhm() + adjacent_peak->fwhm());
          const double peak_separation = fabs( peak->mean() - adjacent_peak->mean() );
          if( peak_separation > 6.0 * avg_fwhm )
          {
            msg += " Warning: the peak at " + std::to_string( peak->mean() )
              + " keV is " + std::to_string( peak_separation / avg_fwhm )
              + " FWHM away from the adjacent peak at "
              + std::to_string( adjacent_peak->mean() )
              + " keV; merging peaks this far apart may not be appropriate.";
          }

          return EditAnalysisPeakStatus{
            true,
            msg,
            updated_peak,
            peaks_in_roi
          };
        }

        case EditPeakAction::DeletePeak:
          // Already handled above
          assert( false );
          break;
      }//switch( options.editAction )

      // Update the peak(s) in the model.
      //  If there are sibling peaks that shared the original continuum, we need to copy them
      //  to share the new continuum, so the ROI remains intact but fully independent from the
      //  originals.
      if( roiSiblings.empty() )
      {
        peakModel->updatePeak( peak, modifiedPeak );
      }else
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

      // Get all peaks in the same ROI for the return value
      vector<shared_ptr<const PeakDef>> peaks_in_roi;
      shared_ptr<const PeakDef> updated_peak;
      shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks_ptr = peakModel->peaks();
      if( all_peaks_ptr )
      {
        const deque<shared_ptr<const PeakDef>> &all_peaks = *all_peaks_ptr;
        const std::shared_ptr<const PeakContinuum> continuum = modifiedPeak.continuum();

        for( const shared_ptr<const PeakDef> &p : all_peaks )
        {
          if( p->continuum() == continuum )
          {
            peaks_in_roi.push_back( p );
            // Find the updated peak by comparing energies
            if( !updated_peak && fabs( p->mean() - modifiedPeak.mean() ) < 0.01 )
              updated_peak = p;
          }
        }
      }

      // If we couldn't find the updated peak in the model, use our local copy
      if( !updated_peak )
        updated_peak = make_shared<const PeakDef>( modifiedPeak );

      return EditAnalysisPeakStatus{
        true,
        "Successfully performed " + string(to_string(options.editAction)) + " on peak at " + std::to_string(modifiedPeak.mean()) + " keV",
        updated_peak,
        peaks_in_roi
      };
    }
    catch( const exception &e )
    {
      return EditAnalysisPeakStatus{
        false,
        string("Error performing ") + to_string(options.editAction) + ": " + e.what(),
        std::nullopt,
        {}
      };
    }//try / catch
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
    const double pair_prod_thresh = isHPGe ? 1255.0 : 2585.0;

    // Initialize result
    EscapePeakCheckStatus result;

    // Calculate potential energies (independent of actual peaks)
    result.potentialSingleEscapePeakEnergy = input_energy - single_escape_offset;
    result.potentialDoubleEscapePeakEnergy = input_energy - double_escape_offset;
    result.potentialParentPeakSingleEscape = input_energy + single_escape_offset;
    result.potentialParentPeakDoubleEscape = input_energy + double_escape_offset;

    // Calculate search windows for each potential energy
    // Window = max(1.0, min(0.5*fwhm, 15.0)) - fairly arbitrary and untested
    const auto calc_window = [&]( const double energy ) -> double {
      if( energy < 0.0 )
        return 1.0; // If energy is negative, just use minimum window

      const double fwhm = get_expected_fwhm( energy, interspec );
      return std::max( 1.0, std::min( 0.5 * fwhm, 15.0 ) );
    };

    // Calculate windows for potential escape peaks of the input energy
    if( result.potentialSingleEscapePeakEnergy > 0.0 )
      result.singleEscapeSearchWindow = calc_window( result.potentialSingleEscapePeakEnergy );
    else
      result.singleEscapeSearchWindow = 1.0;

    if( result.potentialDoubleEscapePeakEnergy > 0.0 )
      result.doubleEscapeSearchWindow = calc_window( result.potentialDoubleEscapePeakEnergy );
    else
      result.doubleEscapeSearchWindow = 1.0;

    // Calculate windows for potential parent peaks (if input is an escape peak)
    result.singleEscapeParentSearchWindow = calc_window( result.potentialParentPeakSingleEscape );
    result.doubleEscapeParentSearchWindow = calc_window( result.potentialParentPeakDoubleEscape );

    // Get all peaks in the spectrum (user + auto-search)
    DetectedPeaksOptions peak_options;
    peak_options.specType = options.specType;
    peak_options.nonBackgroundPeaksOnly = false;
    const DetectedPeakStatus peak_status = detected_peaks( peak_options, interspec );
    const vector<shared_ptr<const PeakDef>> &all_peaks = peak_status.peaks;

    // Helper lambda to find a peak near a given energy
    const auto find_peak_near_energy = [&all_peaks]( const double target_energy, const double window )
        -> shared_ptr<const PeakDef> {
      for( const shared_ptr<const PeakDef> &peak : all_peaks )
      {
        if( !peak )
          continue;

        const double peak_energy = peak->mean();
        if( fabs(peak_energy - target_energy) <= window )
          return peak;
      }
      return nullptr;
    };

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

    // Check if cascade-sum detection should be disabled based on distance
    // Cascade sums typically only occur at distances < 15 cm
    const double cascade_distance_threshold = 15.0 * PhysicalUnits::cm;
    const bool check_cascade_sums = !options.distance.has_value() || (*options.distance <= cascade_distance_threshold);

    // Initialize result
    SumPeakCheckStatus result;

    // Calculate search window for finding contributing peaks
    const double fwhm = get_expected_fwhm( target_energy, interspec );
    result.searchWindow = std::max( 1.0, std::min( 0.5 * fwhm, 15.0 ) );

    // Get all peaks in the spectrum
    DetectedPeaksOptions peak_options;
    peak_options.specType = options.specType;
    peak_options.nonBackgroundPeaksOnly = false;
    const DetectedPeakStatus peak_status = detected_peaks( peak_options, interspec );
    const vector<shared_ptr<const PeakDef>> &all_peaks = peak_status.peaks;

    // Structure to hold candidate sum peak pairs
    struct CandidatePair
    {
      shared_ptr<const PeakDef> peak1;
      shared_ptr<const PeakDef> peak2;
      double sum_energy;
      double energy_error;  // How far from target energy
      SumPeakType sum_type;
      double coincidence_fraction;  // For cascade sums
      double selection_metric;      // For ranking candidates
    };

    vector<CandidatePair> candidates;

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
        if( energy_error <= result.searchWindow )
        {
          // Check if there's a peak at the target energy, and if so, validate that
          // the sum peak is smaller in amplitude than both contributing peaks
          // (A sum peak should always be smaller than the peaks creating it)
          shared_ptr<const PeakDef> target_peak;
          for( const shared_ptr<const PeakDef> &peak : all_peaks )
          {
            if( !peak )
              continue;
            if( fabs(peak->mean() - target_energy) <= result.searchWindow )
            {
              target_peak = peak;
              break;
            }
          }

          // If we found a peak at the target energy, check amplitudes
          if( target_peak )
          {
            const double target_area = target_peak->peakArea();
            const double area1 = peak1->peakArea();
            const double area2 = peak2->peakArea();

            // Sum peak should be smaller than both contributing peaks
            // Skip this candidate if the target peak is larger than either contributor
            if( target_area >= area1 || target_area >= area2 )
              continue;
          }

          CandidatePair candidate;
          candidate.peak1 = peak1;
          candidate.peak2 = peak2;
          candidate.sum_energy = sum_energy;
          candidate.energy_error = energy_error;
          candidate.sum_type = SumPeakType::Unknown;
          candidate.coincidence_fraction = 1.0;
          candidate.selection_metric = 0.0;

          // Check if this could be a cascade-sum peak
          // Only check if cascade detection is enabled (distance <= 15 cm or not specified)
          const SandiaDecay::Nuclide * const nuc1 = peak1->parentNuclide();
          const SandiaDecay::Nuclide * const nuc2 = peak2->parentNuclide();

          if( check_cascade_sums && nuc1 && nuc2 && (nuc1 == nuc2) )
          {
            // Both peaks have the same nuclide source - check for cascade coincidence
            // Get cascade information directly from SandiaDecay
            try
            {
              const double age_in_seconds = PeakDef::defaultDecayTime( nuc1, nullptr ) / PhysicalUnits::second;

              SandiaDecay::NuclideMixture mix;
              const double parent_activity = 1.0 * SandiaDecay::Bq;
              mix.addAgedNuclideByActivity( nuc1, parent_activity, age_in_seconds * SandiaDecay::second );

              const vector<SandiaDecay::NuclideActivityPair> activities = mix.activity( 0.0 );

              // Search for cascade coincidences
              for( const SandiaDecay::NuclideActivityPair &nap : activities )
              {
                const SandiaDecay::Nuclide * const nuclide = nap.nuclide;
                const double activity = nap.activity;

                for( const SandiaDecay::Transition * const transition : nuclide->decaysToChildren )
                {
                  for( const SandiaDecay::RadParticle &particle : transition->products )
                  {
                    if( particle.type != SandiaDecay::GammaParticle )
                      continue;

                    const double br = activity * particle.intensity * transition->branchRatio / parent_activity;

                    for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
                    {
                      const unsigned short int part_ind = particle.coincidences[coinc_index].first;
                      const float fraction = particle.coincidences[coinc_index].second;

                      if( part_ind < transition->products.size() )
                      {
                        const SandiaDecay::RadParticle &coinc_part = transition->products[part_ind];

                        if( coinc_part.type == SandiaDecay::ProductType::GammaParticle )
                        {
                          const float e1 = particle.energy;
                          const float e2 = coinc_part.energy;

                          // Check if this cascade matches our peak pair (allowing for 1 keV tolerance)
                          const bool match1 = (fabs(e1 - energy1) < 1.0 && fabs(e2 - energy2) < 1.0);
                          const bool match2 = (fabs(e1 - energy2) < 1.0 && fabs(e2 - energy1) < 1.0);

                          if( match1 || match2 )
                          {
                            candidate.sum_type = SumPeakType::CascadeSum;
                            candidate.coincidence_fraction = br * fraction;
                            break;
                          }
                        }//if( coinc_part is gamma )
                      }//if( part_ind valid )

                      if( candidate.sum_type == SumPeakType::CascadeSum )
                        break;
                    }//for( coincidences )

                    if( candidate.sum_type == SumPeakType::CascadeSum )
                      break;
                  }//for( particles )

                  if( candidate.sum_type == SumPeakType::CascadeSum )
                    break;
                }//for( transitions )

                if( candidate.sum_type == SumPeakType::CascadeSum )
                  break;
              }//for( activities )
            }
            catch( const std::exception &e )
            {
              // If we can't get cascade information, just continue
              // The peak pair will remain as Unknown type
            }
          }//if( both peaks have same nuclide source )

          // Check if this could be a random-sum peak
          if( (candidate.sum_type != SumPeakType::CascadeSum) && has_valid_timing )
          {
            if( dead_time_fraction > dead_time_threshold )
            {
              // High dead time - could be random summing
              // We'll mark it as such, but further validation would check if the
              // highest amplitude peaks produce expected sum peaks
              candidate.sum_type = SumPeakType::RandomSum;
            }
          }

          // Calculate selection metric
          // For cascade sums: area1 * area2 * coincidence_fraction
          // For others: area1 * area2
          const double area1 = peak1->peakArea();
          const double area2 = peak2->peakArea();
          candidate.selection_metric = area1 * area2 * candidate.coincidence_fraction;

          candidates.push_back( candidate );
        }//if( energy matches within window )
      }//for( j : all_peaks )
    }//for( i : all_peaks )

    // If we found candidates, select the best one
    if( !candidates.empty() )
    {
      // First, sort by energy error (closest match first)
      std::sort( candidates.begin(), candidates.end(),
                []( const CandidatePair &lhs, const CandidatePair &rhs ) {
                  return lhs.energy_error < rhs.energy_error;
                } );

      // If we have multiple candidates with similar energies (within 0.5 keV),
      // select the one with the largest selection metric
      const double best_energy_error = candidates[0].energy_error;
      const double energy_similarity_threshold = 0.5; // keV

      // Find all candidates within the similarity threshold
      vector<CandidatePair> similar_candidates;
      for( const CandidatePair &candidate : candidates )
      {
        if( fabs(candidate.energy_error - best_energy_error) <= energy_similarity_threshold )
          similar_candidates.push_back( candidate );
        else
          break; // Candidates are sorted by energy error
      }

      // Among similar candidates, select the one with the largest metric
      const CandidatePair *best_candidate = &similar_candidates[0];
      for( const CandidatePair &candidate : similar_candidates )
      {
        if( candidate.selection_metric > best_candidate->selection_metric )
          best_candidate = &candidate;
      }

      // Additional validation for random-sum peaks
      if( best_candidate->sum_type == SumPeakType::RandomSum )
      {
        // Find the highest amplitude peak
        double max_area = 0.0;
        double max_area_energy = 0.0;
        for( const shared_ptr<const PeakDef> &peak : all_peaks )
        {
          if( !peak )
            continue;
          if( peak->peakArea() > max_area )
          {
            max_area = peak->peakArea();
            max_area_energy = peak->mean();
          }
        }

        // Check if we should expect to see the largest peak summing with itself
        const double expected_self_sum_energy = 2.0 * max_area_energy;
        const bool self_sum_in_range = (expected_self_sum_energy <= spectrum->gamma_energy_max());
        const bool checking_self_sum = (fabs(target_energy - expected_self_sum_energy) <= result.searchWindow);

        // If the largest peak should sum with itself but we're not seeing that,
        // and we're not currently checking for the self-sum peak, then this is probably
        // not a random-sum peak
        if( self_sum_in_range && !checking_self_sum )
        {
          // Check if a peak exists at twice the highest energy
          bool found_self_sum = false;
          for( const shared_ptr<const PeakDef> &peak : all_peaks )
          {
            if( !peak )
              continue;
            if( fabs(peak->mean() - expected_self_sum_energy) <= result.searchWindow )
            {
              found_self_sum = true;
              break;
            }
          }

          if( !found_self_sum )
          {
            // The highest peak doesn't sum with itself, so random summing is unlikely
            // Downgrade to Unknown
            const_cast<CandidatePair*>(best_candidate)->sum_type = SumPeakType::Unknown;
          }
        }
      }//if( random sum validation )

      // Create the result
      SumPeakInfo sum_info;
      sum_info.sumType = best_candidate->sum_type;
      sum_info.firstPeak = best_candidate->peak1;
      sum_info.secondPeak = best_candidate->peak2;

      if( best_candidate->sum_type == SumPeakType::CascadeSum )
        sum_info.coincidenceFraction = best_candidate->coincidence_fraction;

      // Generate labels
      char buffer[512];
      const double e1 = best_candidate->peak1->mean();
      const double e2 = best_candidate->peak2->mean();

      // Try to create a source label if both peaks have sources
      const SandiaDecay::Nuclide * const nuc1 = best_candidate->peak1->parentNuclide();
      const SandiaDecay::Nuclide * const nuc2 = best_candidate->peak2->parentNuclide();
      const ReactionGamma::Reaction * const rxn1 = best_candidate->peak1->reaction();
      const ReactionGamma::Reaction * const rxn2 = best_candidate->peak2->reaction();
      const SandiaDecay::Element * const xray1 = best_candidate->peak1->xrayElement();
      const SandiaDecay::Element * const xray2 = best_candidate->peak2->xrayElement();

      if( nuc1 && nuc2 && (nuc1 == nuc2) && (best_candidate->sum_type == SumPeakType::CascadeSum) )
      {
        snprintf( buffer, sizeof(buffer), "%s cascade-sum %.2f + %.2f keV",
                 nuc1->symbol.c_str(), e1, e2 );
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

        const char *sum_type_str = (best_candidate->sum_type == SumPeakType::RandomSum) ? "random-sum" : "sum";
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
      const shared_ptr<const DetectorPeakResponse> &detector )
  {
    ComptonPeakCheckStatus result;

    if( !candidate || !candidate->gausPeak() || !user_peaks || user_peaks->empty() )
      return result;

    const double electron_rest_mass = 510.9989; // keV; matches escape_peak_check at AnalystChecks.cpp:2025

    result.candidateMean = candidate->mean();
    result.candidateFwhm = candidate->fwhm();
    const double cand_amp = candidate->amplitude();

    // Cheap-path FWHM only: use DRF resolution info if available; otherwise
    // -1 sentinel triggers the sqrt(E)-vs-parent fallback below.
    double expected_det_fwhm = -1.0;
    if( detector && detector->hasResolutionInfo() )
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

    ParentMatch best;

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

      if( this_score < best.energy_match_score )
      {
        best.parent = p;
        best.type = this_type;
        best.scatter_angle_deg = this_angle;
        best.energy_match_score = this_score;
      }
    }//for( each potential parent )

    if( !best.parent )
      return result;

    const shared_ptr<const PeakDef> parent = best.parent;
    const double parent_amp = parent->amplitude();

    // Hard guard: parent must be at least as large as candidate.
    if( parent_amp <= cand_amp )
      return result;

    const double area_fraction = (parent_amp > 0.0) ? (cand_amp / parent_amp) : 0.0;
    if( area_fraction > 0.5 )
      return result;

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
      return result;

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

    char label_buf[160];
    if( parent->parentNuclide() )
    {
      snprintf( label_buf, sizeof(label_buf), "Compton scatter from %s %.2f keV peak",
                parent->parentNuclide()->symbol.c_str(), parent->mean() );
    }else if( parent->reaction() )
    {
      snprintf( label_buf, sizeof(label_buf), "Compton scatter from %s %.2f keV peak",
                parent->reaction()->name().c_str(), parent->mean() );
    }else
    {
      snprintf( label_buf, sizeof(label_buf), "Compton scatter from %.2f keV peak", parent->mean() );
    }

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

    PeakModel *peak_model = interspec->peakModel();
    if( !peak_model )
      return ComptonPeakCheckStatus{};

    const shared_ptr<const deque<shared_ptr<const PeakDef>>> user_peaks
                                                  = peak_model->peaks( options.specType );
    if( !user_peaks || user_peaks->empty() )
      return ComptonPeakCheckStatus{};

    // Find the candidate peak: nearest user peak to options.energy.
    shared_ptr<const PeakDef> candidate;
    double smallest_diff = std::numeric_limits<double>::max();
    for( const shared_ptr<const PeakDef> &p : *user_peaks )
    {
      if( !p || !p->gausPeak() )
        continue;
      const double diff = std::abs( p->mean() - options.energy );
      if( diff < smallest_diff )
      {
        smallest_diff = diff;
        candidate = p;
      }
    }

    if( !candidate )
      return ComptonPeakCheckStatus{};

    const double tol = std::max( 0.5, candidate->fwhm() );

    if( smallest_diff > tol )
    {
      ComptonPeakCheckStatus empty;
      empty.searchWindow = tol;
      return empty;
    }

    ComptonPeakCheckStatus result = compton_peak_check( candidate, user_peaks, spectrum,
                                                        meas->detector() );
    result.searchWindow = tol;
    return result;
  }//ComptonPeakCheckStatus compton_peak_check( const ComptonPeakCheckOptions &options, InterSpec *interspec )

} // namespace AnalystChecks
