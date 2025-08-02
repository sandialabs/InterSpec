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
  
} // namespace AnalystChecks 