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
#include <deque>
#include <set>
#include <algorithm>
#include <limits>

#include "InterSpec/AnalystChecks.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"

#include <Wt/WApplication>

using namespace std;

namespace AnalystChecks
{
  std::string hello_world()
  {
    return "Hello World from AnalystChecks namespace!";
  }
  
  DetectedPeakStatus detected_peaks(const DetectedPeaksOptions& options, InterSpec* interspec)
  {
    if (!interspec) {
      throw std::runtime_error("No InterSpec session available");
    }
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
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
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
      const auto det = meas->detector();
      const vector<shared_ptr<const PeakDef>> found_auto_peaks
                = ExperimentalAutomatedPeakSearch::search_for_peaks(spectrum, det, user_peaks, singleThreaded, isHPGe);
        
      auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>(begin(found_auto_peaks), end(found_auto_peaks));
      auto_peaks = autopeaksdeque;
      meas->setAutomatedSearchPeaks(sample_nums, autopeaksdeque);
    }

    vector<shared_ptr<const PeakDef>> all_peaks;

    if (user_peaks)
      all_peaks.insert(all_peaks.end(), user_peaks->begin(), user_peaks->end());

    // We will add auto-search peaks only if there isnt already a user peak at essentially the same energy
    if (auto_peaks) {
      //Lambda to find the nearest peak, so far, to a given energy.
      auto nearest_peak = [&all_peaks](const float energy) -> std::shared_ptr<const PeakDef> {
        std::shared_ptr<const PeakDef> nearest;
        double minDE = std::numeric_limits<double>::infinity();

        for (const auto &peak : all_peaks) {
          const double dE = fabs(peak->mean() - energy);
          if ((dE < minDE)
              && ((energy > peak->lowerX()) && (energy < peak->upperX()))) {
            minDE = dE;
            nearest = peak;
          }
        }

        return nearest;
      };

      for (const shared_ptr<const PeakDef> &peak : *auto_peaks) {
        auto nearpeak = nearest_peak(peak->mean());
        const double peak_sigma = peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth();
        if (!nearpeak || (fabs(nearpeak->mean() - peak->mean()) > peak_sigma))
          all_peaks.push_back(peak);
      }
    }
    
    std::sort(begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr);
    
    DetectedPeakStatus result;
    result.userSession = user_session;
    result.peaks = std::move(all_peaks);
    
    return result;
  }
  
  
  FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
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
    
    
    const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );
    
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
    foundPeaks = searchForPeakFromUser( options.energy, pixelPerKev, data, origPeaks, det, isHPGe );
    
    vector<shared_ptr<const PeakDef>> &peaks_to_add_in = foundPeaks.first;
    const vector<shared_ptr<const PeakDef>> &peaks_to_remove = foundPeaks.second;
    
    if( peaks_to_add_in.empty() )
      throw runtime_error( "Could not fit peak." );
    
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
    
    if( options.source.has_value() )
    {
      const string source = options.source.value();
      
      shared_ptr<PeakDef> new_fit_peak = make_shared<PeakDef>( *fit_peak );
      
      PeakModel::SetGammaSource res = PeakModel::setNuclideXrayReaction( *new_fit_peak, source, 4.0 );
      
      if( res == PeakModel::SetGammaSource::NoSourceChange )
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
          new_fit_peak->setLineColor( c );
        }
      }
    }//if( options.source.has_value() )
    
    
    if( options.addToUsersPeaks )
    {
      if( options.specType == SpecUtils::SpectrumType::Foreground )
      {
        PeakModel *pmodel = interspec->peakModel();
        assert( pmodel );
        if( pmodel )
        {
          pmodel->removePeaks( peaks_to_remove );
          pmodel->addPeaks( peaks_to_add_in );
        }//if( pmodel )
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
        interspec->setPeaks( options.specType, new_deque );
      }
    }//if( options.addToUsersPeaks )
    
    
    
    FitPeakStatus result;
    
    result.userSession = user_session;
    result.fitPeak = fit_peak;
    result.peaksInRoi = peaks_to_add_in;
    
    return result;
  }//FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );
  
  
  GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec )
  {
    if (!interspec)
      throw std::runtime_error("No InterSpec session available");
    
    // Get session ID from InterSpec (using Wt application)
    string user_session = "direct_call";
    if (Wt::WApplication* app = Wt::WApplication::instance()) {
      user_session = app->sessionId();
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
    
    shared_ptr<deque<shared_ptr<const PeakDef>>> meas_peaks = meas->peaks(sample_nums);
    assert( meas_peaks );
    if( !meas_peaks )
      throw runtime_error( "Unexpected error getting existing peak list" );
    
    
    GetUserPeakStatus status;
    status.userSession = user_session;
    if( meas_peaks )
      status.peaks.insert( end(status.peaks), begin(*meas_peaks), end(*meas_peaks) );
    
    return status;
  }
  
  
  std::vector<float> get_characteristic_gammas( const std::string &nuclide )
  {
    //
    //const vector<tuple<string, const SandiaDecay::Nuclide *, float>> &nuc_characteristics = IsotopeId::characteristicGammas();
    //blah blah blah
  }
} // namespace AnalystChecks
