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
  class EnergyCalibration;
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
    /** If specified to be true, then instead of looking for just the peaks in the exemplar file, a search for all peaks in the specrtum will be performed. */
    bool fit_all_peaks;
    bool to_stdout;
    bool refit_energy_cal;
    bool use_exemplar_energy_cal;
    bool write_n42_with_results;
    bool show_nonfit_peaks;
    bool overwrite_output_files;
    bool create_csv_output;
    bool create_json_output;
    std::string output_dir;
    std::string background_subtract_file;
    std::set<int> background_subtract_samples;
    std::shared_ptr<SpecMeas> cached_background_subtract_spec;
    bool use_existing_background_peaks;
    bool use_exemplar_energy_cal_for_background;
    
    /** The improvement to the Chi2 of a peak fit required, over just fitting the continuum, to the ROI.
     
     A negative or zero value indicates no requirement (and default, since we are asserting peak
     is likely in the spectrum for batch analysis), and for general peak searching, reasonable
     values are between ~1 (a weak peak) and ~5 (a significant peak).
     */
    double peak_stat_threshold;
    
    /** Specifies how well the peak must match in shape to a gaussian in order to keep the peak.
     
     The higher this number, the more like a gaussian the fit peak is.
     It is the ratio of the null hypothesis chi2 (continuum only, no Gaussian),
     to the test hypothesis (continuum + Gaussian) chi2.
     A reasonable value for this seems to be ~4.
     A zero or negative value will mean no requirement, and also no
     `peak_stat_threshold` requirement.
     */
    double peak_hypothesis_threshold;
    
    /** The directory to allow report template to look in to include other templates.
     If specified, then the standard report directory cant be used.
     */
    std::string template_include_dir;
    
    
    /** File paths to report templates, that will be saved for each input files. 
     
      If the string contains the string ':--DisplayName--:', then everything before this string will
      be the path to the template, and everything after will be the display name of the template.
      \sa sm_report_display_name_marker
    */
    std::vector<std::string> report_templates;
    
    /** File path to report templates that summarizes all input files. 
     
     Similar to `report_templates`, the string may be delimited by ':--DisplayName--:', to specify 
     the display name of the template.
     \sa sm_report_display_name_marker
    */
    std::vector<std::string> summary_report_templates;

    /** Delimeter used within report template filesystem paths to seperate the filesystem path, and the display name of the template.
     If this delimeter is not present, then the full string is used for both these properties.
     
     This is used for report templates uploaded via HTML, and spooled to disk, so we can still generate reasonably named reports.
     */
    static const char * const sm_report_display_name_marker; // = ":--DisplayName--:"
  };//struct BatchPeakFitOptions

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

    /** Background spectrum that was subtracted from the foreground, to make `spectrum`, if any.

     The background subtraction can either be on a peak-by-peak basis, or a hard
     background subtraction, see `BatchPeakFitOptions::use_exemplar_energy_cal_for_background`.
     */
    std::shared_ptr<SpecUtils::Measurement> background;

    bool success;
    std::vector<std::string> warnings;

    /** The original energy calibration of the spectrum, before re-fitting it (if done). */
    std::shared_ptr<const SpecUtils::EnergyCalibration> original_energy_cal;

    /** The energy calibration after fitting for it - will only be non-null if energy calibration was performed */
    std::shared_ptr<const SpecUtils::EnergyCalibration> refit_energy_cal;
  };//struct BatchPeakFitResult


  struct BatchPeakFitSummaryResults
  {
    BatchPeakFitOptions options;
    std::string exemplar_filename;
    std::shared_ptr<const SpecMeas> exemplar;
    std::set<int> exemplar_sample_nums;

    /** Each of these next four `file_*` variables will have the same number of entries, as the input number of files. */
    std::vector<BatchPeakFitResult> file_results;
    std::vector<std::string> file_json;
    std::vector<std::string> file_peak_csvs;
    std::vector<std::vector<std::string>> file_reports;

    std::string summary_json;
    std::vector<std::string> summary_reports;
    std::vector<std::string> warnings;
  };//struct BatchPeakFitSummaryResults


  /** Fits the peaks in a number of spectrum files, prodicing individual reports (if wanted) as well as a summary report.
   @param exemplar_filename The name of the spectrum file to use as the exemplar file.  This may be an N42

   */
  InterSpec_API void fit_peaks_in_files( const std::string &exemplar_filename,
                          std::shared_ptr<const SpecMeas> optional_parsed_exemplar_n42,
                          const std::set<int> &exemplar_sample_nums,
                          const std::vector<std::string> &files,
                          std::vector<std::shared_ptr<SpecMeas>> cached_files,
                          const BatchPeakFitOptions &options,
                          BatchPeakFitSummaryResults *results = nullptr );


  
  
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
  
  /** Function that applies the energy calibration from the exemplar spectrum, to a spectrum from a different file.
   
   @param energy_cal The energy calibration to apply to `to_spectrum`, and optionally `to_specfile`
   @param to_spectrum The spectrum, which may or may not be in `to_specfile`, to apply the energy calibration from `from_spectrum`
   @param to_specfile The (optional) spectrum file to apply the energy calibration to; this will also take care of shifting peak energies
   @param used_sample_nums The sample numbers used to create the `to_spectrum` from the `to_specfile` - used to keep peaks from
          being moved twice.
   */
  void propagate_energy_cal( const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
                                          std::shared_ptr<SpecUtils::Measurement> &to_spectrum,
                                          std::shared_ptr<SpecMeas> &to_specfile,
                                          const std::set<int> &used_sample_nums );
  
  /** Finds the spectrum, peaks, and sample numbers to use from the exemplar file.
   
   @param [out] exemplar_spectrum Will be set to the spectrum to use as the exemplar spectrum (may be a sum of
          multiple Measurements)
   @param [out] exemplar_peaks Will be set to the peaks to use from the exemplar file.
   @param [in|out] exemplar_sample_nums Sample numbers to use in the exemplar file.  If non-empty, these sample
          number will be used to retrieve the spectrum to use.  Contents will be set to the used sample numbers.
   @param [in] exemplar_n42 The spectrum file to retrieve the spectrum/peaks for
   */
  void get_exemplar_spectrum_and_peaks(
                                       std::shared_ptr<const SpecUtils::Measurement> &exemplar_spectrum,
                                       std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &exemplar_peaks,
                                       std::set<int> &exemplar_sample_nums,
                                       const std::shared_ptr<const SpecMeas> &exemplar_n42 );
}//namespace BatchPeak

#endif //BatchPeak_h
