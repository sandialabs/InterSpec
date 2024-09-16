#ifndef BatchActivity_h
#define BatchActivity_h
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
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "InterSpec/BatchPeak.h"

// Forward declarations
class PeakDef;
class SpecMeas;
class DetectorPeakResponse;
namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils

namespace ShieldingSourceFitCalc
{
  struct ModelFitResults;
}

/** The functions necessary to batch-fit activity and shielding.  */
namespace BatchActivity
{
  /** Returns the DRF specified by the user on the command line.
   
   Throws exception if DRF file is invalid, or specified detector could not be loaded.
   If both input strings are empty, returns nullptr.
   */
  InterSpec_API std::shared_ptr<DetectorPeakResponse> init_drf_from_name( std::string drf_file, std::string drf_name );
  
  struct InterSpec_API BatchActivityFitOptions
    : public BatchPeak::BatchPeakFitOptions
  {
    bool use_bq = false;
    std::shared_ptr<DetectorPeakResponse> drf_override;
    boost::optional<double> distance_override;
    bool hard_background_sub;
  };//struct BatchActivityFitOptions
  
  
  InterSpec_API void fit_activities_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          const BatchActivityFitOptions &options );
  
  struct InterSpec_API BatchActivityFitResult
  {
    enum class InterSpec_API ResultCode
    {
      CouldntInitializeStaticResources,
      NoExemplar,
      CouldntOpenExemplar,
      ErrorPickingSpectrumFromExemplar,
      CouldntOpenInputFile,
      CouldntOpenBackgroundFile,
      NoInputSrcShieldModel,
      ForegroundSampleNumberUnderSpecified,
      BackgroundSampleNumberUnderSpecified,
      InvalidLiveTimeForHardBackSub,
      SpecifiedDistanceWithFixedGeomDet,
      ErrorWithHardBackgroundSubtract,
      ErrorApplyingExemplarEneCalToFore,
      
      ForegroundPeakFitFailed,
      BackgroundPeakFitFailed,
      NoExistingBackgroundPeaks,
      NoFitForegroundPeaks,
      NoDetEffFnct,
      InvalidDistance,
      InvalidGeometry,
      InvalidFitOptions,
      ExemplarUsedBackSubButNoBackground,
      NoShieldingsNode,
      ErrorParsingShielding,
      MissingNuclidesNode,
      InvalidNuclideNode,
      NoSourceNuclides,
      
      FitNotSuccessful,
      DidNotFitAllSources,
      
      UnknownStatus,
      Success
    };//enum class ResultCode
    
    static const char *to_str( const ResultCode code );
    
    ResultCode m_result_code;
    std::string m_error_msg;
    std::vector<std::string> m_warnings;
    
    std::string m_filename;
    std::shared_ptr<const SpecMeas> m_foreground_file;
    std::set<int> m_foreground_sample_numbers;  ///!< right now limited to single sample - may change in future
    std::shared_ptr<const SpecUtils::Measurement> m_foreground;
    
    std::shared_ptr<const SpecMeas> m_exemplar_file;
    std::set<int> m_exemplar_sample_numbers;    ///!< right now limited to single sample - may change in future
    std::shared_ptr<const SpecUtils::Measurement> m_exemplar;
    
    std::shared_ptr<const SpecMeas> m_background_file;  ///!< May be same as foreground or exemplar
    std::set<int> m_background_sample_numbers;  ///!< right now limited to single sample - may change in future
    std::shared_ptr<const SpecUtils::Measurement> m_background;
    
    BatchActivityFitOptions m_options;
    
    std::shared_ptr<const BatchPeak::BatchPeakFitResult> m_peak_fit_results;
    std::shared_ptr<const BatchPeak::BatchPeakFitResult> m_background_peak_fit_results;
    
    std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> m_fit_results;
  };//struct BatchActivityFitResult
  
  
  
  /** Fits the peaks, activities, and shieldins
   
   exemplar peaks for a given file.
   
   @param exemplar_filename The file-path of the N42-2012 file with the example peaks, or the file-path of the CSV with peak info
   @param exemplar_sample_nums If a N42-2012 file is used for exemplar, and which peaks to use is ambiguous, these
          sample numbers specify which peaks to use.  Must be blank if exemplar is CSV file, or if N42-2012 file, this combination
          of sample numbers must specify peaks to use.
   @param cached_exemplar_n42 If non-null, then `exemplar_filename` will be ignored, and this file will be used; to avoid re-parsing
          of the exemplar file over-and-over again.
   @param filename The name of the spectrum file to fit peaks to.
   @param options The options to use for fitting peaks; note, not all options are used, as some of them are only applicable to
          #fit_peaks_in_files
   */
  InterSpec_API BatchActivityFitResult fit_activities_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          const BatchActivityFitOptions &options );
  
}//namespace BatchPeak

#endif //BatchPeak_h
