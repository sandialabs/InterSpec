#ifndef FarmAnalysis_h
#define FarmAnalysis_h

#include "InterSpec_config.h"

#include <string>
#include <memory>
#include <vector>

// Forward declarations
class PeakDef;
class DetectorPeakResponse;
struct SpecFileInfoToQuery;

namespace SpecUtils { class Measurement; class SpecFile; enum class DetectorType : int; }
namespace Farm { struct FarmOptions; struct EnrichmentResults; }

namespace Farm
{

/** Performs automated peak search and stores results as JSON in info.farm_peaks_json.

 @param foreground The foreground spectrum to search for peaks in
 @param drf Optional detector response function (can be nullptr)
 @param is_hpge Whether this is an HPGe detector (affects peak filtering)
 @param info The SpecFileInfoToQuery to store results in (farm_peaks_json field)
 */
void perform_peak_search(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::shared_ptr<const DetectorPeakResponse> &drf,
    const bool is_hpge,
    SpecFileInfoToQuery &info );


/** Runs GADRAS Full Spectrum Isotope ID analysis.

 @param exe_path Path to the GADRAS executable
 @param spec_file SpecFile containing foreground (and optionally background)
 @param detector_type The detector type for DRF selection
 @param synthesize_background If true, synthesize background if not present
 @returns JSON string of ExternalRidResults, or empty on failure
 */
std::string run_gadras_full_spectrum_id_analysis(
    const std::string &exe_path,
    std::shared_ptr<SpecUtils::SpecFile> spec_file,
    const SpecUtils::DetectorType detector_type,
    const bool synthesize_background );


/** Check if peaks indicate uranium isotopics should be performed.

 Checks for: 185 keV peak (within 0.75 FWHM) AND (205.3 OR 143.8 keV peak)

 @param peaks_json JSON array of peak objects with "mean" and "fwhm" fields
 @returns true if U isotopics should be performed
 */
bool should_do_uranium_isotopics( const std::string &peaks_json );


/** Check if peaks indicate plutonium isotopics should be performed.

 Checks for: 375 keV peak AND one of (129.3, 413.7, 662 keV)

 @param peaks_json JSON array of peak objects with "mean" and "fwhm" fields
 @returns true if Pu isotopics should be performed
 */
bool should_do_plutonium_isotopics( const std::string &peaks_json );


/** Run RelActCalcAuto isotopics using the given options XML file.

 @param foreground Foreground spectrum
 @param background Background spectrum (can be nullptr)
 @param peaks Peaks to use (can be empty for auto-detection)
 @param options_xml_path Full path to the RelActCalcAuto options XML file
 @param drf Detector response function (can be nullptr for auto-detection)
 @returns EnrichmentResults struct (check warnings for errors)
 */
EnrichmentResults run_relact_isotopics(
    const std::shared_ptr<const SpecUtils::Measurement> &foreground,
    const std::shared_ptr<const SpecUtils::Measurement> &background,
    const std::vector<std::shared_ptr<const PeakDef>> &peaks,
    const std::string &options_xml_path,
    const std::shared_ptr<const DetectorPeakResponse> &drf );


/** Run FRAM executable for isotopics analysis.

 TODO: This function signature will need to be refined once FRAM calling
 conventions are specified by the user.

 @param fram_exe_path Path to FRAM executable
 @param fram_output_path Directory where FRAM writes output
 @param input_n42_path Path to N42 file to analyze
 @param is_uranium true for uranium, false for plutonium
 @returns EnrichmentResults struct (check warnings for errors)
 */
EnrichmentResults run_fram_isotopics(
    const std::string &fram_exe_path,
    const std::string &fram_output_path,
    const std::string &input_n42_path,
    const bool is_uranium,
    const bool is_plutonium );


/** Maps SpecUtils::DetectorType to GADRAS DRF name.

 @param type The detector type
 @returns DRF name string, or empty if detector type is not known/supported
 */
std::string detector_type_to_gadras_drf( const SpecUtils::DetectorType type );



/** Write a .farm.fertilized.n42 file with FARM results as JSON remark.

 @param original_file_path Path to original spectrum file
 @param meas The SpecFile to write (should be copy with fore/back only)
 @param info The SpecFileInfoToQuery containing FARM results
 */
void write_fertilized_n42(
    const std::string &original_file_path,
    const SpecUtils::SpecFile &meas,
    const SpecFileInfoToQuery &info );

} // namespace Farm

#endif // FarmAnalysis_h
