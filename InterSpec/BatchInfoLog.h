#ifndef BatchInfoLog_h
#define BatchInfoLog_h
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

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

// Forward declarations
class SpecMeas;
class DetectorPeakResponse;


namespace ShieldingSourceFitCalc
{
  struct ModelFitResults;
  struct ShieldingSourceFitOptions;
}

namespace GammaInteractionCalc
{
  struct PeakDetail;
  struct PeakDetailSrc;
  struct SourceDetails;
  struct ShieldingDetails;
  struct ShieldingSourceFitOptions;
  enum class GeometryType : int;
}

namespace BatchPeak
{
  struct BatchPeakFitResult;
}


namespace SpecUtils
{
  class Measurement;
}

namespace BatchActivity
{
  struct BatchActivityFitOptions;
}


namespace BatchInfoLog
{
  void add_basic_src_details( const GammaInteractionCalc::SourceDetails &src,
                            const std::shared_ptr<const DetectorPeakResponse> &drf,
                            const bool useBq,
                            const std::vector<GammaInteractionCalc::ShieldingDetails> *shield_details,
                             nlohmann::basic_json<> &src_json );
  
  void add_act_shield_fit_options_to_json( const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options,
                               const double distance,
                               const GammaInteractionCalc::GeometryType geometry,
                               const std::shared_ptr<const DetectorPeakResponse> &drf,
                                          nlohmann::basic_json<> &data );
  
  /** Adds basic information about a peak (energy, fwhm, counts, etc), but not any information
     about gammas that contribute to it, etc
   */
  void add_basic_peak_info( const GammaInteractionCalc::PeakDetail &peak, nlohmann::basic_json<> &peak_json );
  
  void shield_src_fit_results_to_json( const ShieldingSourceFitCalc::ModelFitResults &results,
                                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                                      const bool useBq,
                                      nlohmann::basic_json<> &data );
  
  void add_gamma_info_for_peak( const GammaInteractionCalc::PeakDetailSrc &ps,
                    const GammaInteractionCalc::SourceDetails * const src,
                    const std::shared_ptr<const DetectorPeakResponse> &drf,
                    const bool useBq,
                    const std::vector<GammaInteractionCalc::ShieldingDetails> * const shield_details,
                               nlohmann::basic_json<> &gamma_json );
  
  
  void add_hist_to_json( nlohmann::basic_json<> &data,
                       const bool is_background,
                       const std::shared_ptr<const SpecUtils::Measurement> &spec_ptr,
                       const std::shared_ptr<const SpecMeas> &spec_file,
                       const std::set<int> &sample_numbers,
                       const std::string &filename,
                        const std::shared_ptr<const BatchPeak::BatchPeakFitResult> &peak_fit );
  
  
  void add_activity_fit_options_to_json( nlohmann::basic_json<> &data,
                                        const BatchActivity::BatchActivityFitOptions &options );
  
  void add_exe_info_to_json( nlohmann::basic_json<> &data );
  
  void add_peak_fit_options_to_json( nlohmann::basic_json<> &data, const BatchPeak::BatchPeakFitOptions &options );
}//namespace BatchInfoLog

#endif //BatchInfoLog_h
