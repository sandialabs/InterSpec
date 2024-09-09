#ifndef BatchPeak_h
#define BatchPeak_h
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

// Forward declarations
class PeakDef;
class SpecMeas;

namespace SpecUtils
{
  class Measurement;
}//namespace SpecUtils


namespace BatchPeak
{
  /** Fits the peaks in an 'exemplar' file, in a number of other similar file.
   A work in progress.
   
   Currently the results aren't too great, could use:
   - Limits, that are maybe adjustable, for how significant a peak has to be, before it is kept
   - The FWHM needs to be enforced to be reasonable - like after fitting the peaks, fit the FWHM, and then constrain the FWHM to be
      something reasonable around that (hopefully the higher-statistics peaks will dominate the FWHM - or instead could do auto
      peak search, and use those to get the FWHM - this is probably the better idea - or maybe some combination)
   - The FWHM within a single ROI should be constrained
   
   TODO:
   - specify a struct that contains all the options for the fit - we will likely be adding more options, like statistical significance
   - have a .ini file able to back the command line options, so users can specify their own default options
   - When peaks are not fit for, print out their Currie detection limit
   */
  
  struct InterSpec_API BatchPeakFitOptions
  {
    bool to_stdout;
    bool refit_energy_cal;
    bool use_exemplar_energy_cal;
    bool write_n42_with_results;
    bool show_nonfit_peaks;
    bool overwrite_output_files;
    bool create_csv_output;
    std::string output_dir;
    std::string background_subtract_file;
    std::set<int> background_subtract_samples;
    
    /** The directory to allow report template to look in to include other templates.
     If specified, then the standard report directory cant be used.
     */
    std::string template_include_dir;
    
    /** File paths to report templates, that will be saved for each input files. */
    std::vector<std::string> report_templates;
    
    /** File path to report templates that summarizes all input files. */
    std::vector<std::string> summary_report_templates;
  };//struct BatchPeakFitOptions
  
  
  InterSpec_API void fit_peaks_in_files( const std::string &exemplar_filename,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          const BatchPeakFitOptions &options );
  
  struct InterSpec_API BatchPeakFitResult
  {
    std::string file_path;
    BatchPeakFitOptions options;
    
    std::shared_ptr<const SpecMeas> exemplar;
    std::set<int> exemplar_sample_nums;
    std::deque<std::shared_ptr<const PeakDef>> exemplar_peaks;
    std::shared_ptr<const SpecUtils::Measurement> exemplar_spectrum;
    std::vector<std::shared_ptr<const PeakDef>> unfit_exemplar_peaks;  //Exemplar peaks not found in the spectrum
    
    std::shared_ptr<SpecMeas> measurement;
    std::shared_ptr<SpecUtils::Measurement> spectrum;
    std::set<int> sample_numbers;
    std::deque<std::shared_ptr<const PeakDef>> fit_peaks;
    
    /** Background spectrum that was subtracted from the foreground, to make `spectrum`, if any. */
    std::shared_ptr<const SpecUtils::Measurement> background;
    
    bool success;
    std::vector<std::string> warnings;
  };//struct BatchPeakFitResult
  
  
  /** Fits the exemplar peaks for a given file.
   
   @param exemplar_filename The file-path of the N42-2012 file with the example peaks, or the file-path of the CSV with peak info
   @param exemplar_sample_nums If a N42-2012 file is used for exemplar, and which peaks to use is ambiguous, these
          sample numbers specify which peaks to use.  Must be blank if exemplar is CSV file, or if N42-2012 file, this combination
          of sample numbers must specify peaks to use.
   @param cached_exemplar_n42 If non-null, then `exemplar_filename` will be ignored, and this file will be used; to avoid re-parsing
          of the exemplar file over-and-over again.
   @param filename The name of the spectrum file to fit peaks to.
   @param cached_spectrum If you have already parsed/opened the `filename` spectrum file, you can provide it here to
          avoid overhead of re-parsing it.
   @param foreground_sample_numbers The sample numbers to fit the peaks to.  If left empty, will try to automatically determine.
   @param options The options to use for fitting peaks; note, not all options are used, as some of them are only applicable to
          #fit_peaks_in_files
   */
  InterSpec_API BatchPeakFitResult fit_peaks_in_file( const std::string &exemplar_filename,
                          std::set<int> exemplar_sample_nums,
                          std::shared_ptr<const SpecMeas> cached_exemplar_n42,
                          const std::string &filename,
                          std::shared_ptr<SpecMeas> cached_spectrum,
                          std::set<int> foreground_sample_numbers,
                          const BatchPeakFitOptions &options );
  
}//namespace BatchPeak

#endif //BatchPeak_h
