#ifndef AnalystChecks_h
#define AnalystChecks_h
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
#include <vector>
#include <optional>
#include <memory>
#include <variant>

#include "SpecUtils/SpecFile.h"

#include "InterSpec/ReactionGamma.h" //for ReactionGamma::Reaction

// Forward declarations
class PeakDef;
class InterSpec;

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}

namespace AnalystChecks
{
  /** Options for peak detection analysis. */
  struct DetectedPeaksOptions {
    SpecUtils::SpectrumType specType;
    bool nonBackgroundPeaksOnly;
    std::optional<double> lowerEnergy;
    std::optional<double> upperEnergy;
  };
  
  /** Results of peak detection analysis. */
  struct DetectedPeakStatus {
    std::vector<std::shared_ptr<const PeakDef>> peaks;

    /** Subset of `peaks` that are user-defined analysis peaks (not auto-detected peaks). */
    std::vector<std::shared_ptr<const PeakDef>> analysis_peaks;
  };
  
  /** Perform automated peak detection on the specified spectrum.
   
   This function combines user-defined peaks with automatically detected peaks,
   avoiding duplicates by checking peak overlap.
   
   @param options Specifies which spectrum to analyze and optional user session
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return DetectedPeakStatus containing the found peaks and session info
   @throws std::runtime_error if spectrum is not available or other errors occur
   */
  InterSpec_API DetectedPeakStatus detected_peaks(const DetectedPeaksOptions& options, InterSpec* interspec);
  
  
  struct FitPeakOptions {
    bool doNotAddToAnalysisPeaks;
    double energy;
    SpecUtils::SpectrumType specType;
    std::optional<std::string> source;
  };
  
  struct FitPeakStatus {
    std::shared_ptr<const PeakDef> fitPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;
  };
  
  InterSpec_API FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );
  
  struct GetUserPeakOptions {
    SpecUtils::SpectrumType specType;
    std::optional<double> lowerEnergy;
    std::optional<double> upperEnergy;
  };
  
  struct GetUserPeakStatus {
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  InterSpec_API GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec );

  /** Get a list of unique sources assigned to peaks in the specified spectrum.

   Returns a vector of source strings (nuclide symbols, element symbols, or reaction names)
   that are assigned to user peaks in the specified spectrum. Each source appears only once
   in the result, even if multiple peaks are assigned to that source.

   @param specType Which spectrum to get sources from (Foreground, Background, or SecondForeground)
   @param interspec Pointer to the InterSpec session
   @return Vector of unique source strings assigned to peaks
   @throws std::runtime_error if InterSpec is null or spectrum not loaded
   */
  InterSpec_API std::vector<std::string> get_identified_sources( const SpecUtils::SpectrumType specType, InterSpec *interspec );

  struct FitPeaksForNuclideOptions {
    std::vector<std::string> sources;
    bool doNotAddPeaksToUserSession;
  };
  
  struct FitPeaksForNuclideStatus
  {
    std::vector<std::shared_ptr<const PeakDef>> fitPeaks;
    std::vector<std::string> warnings;
  };
  
  InterSpec_API FitPeaksForNuclideStatus fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options, InterSpec *interspec );
  
  /** Calculates an approximate importance that a peak in a spectrum will have for a given nuclide.
   * 
   * Importance is defined by `yield(i)*sqrt(energy(i))/sum(yield*sqrt(energy))`.
   * A very rough approximation, but can be useful for matching things.
   * This is what PhotoPeak.lis uses for its fourth column.
   * 
   * @param gamma_energies_and_yields A vector of tuples containing the energy and yield of each gamma.
   * @return A tuple containing the energy, yield, and estimated importance.
   */
  std::vector<std::tuple<float,float,float>> cacl_estimated_gamma_importance( const std::vector<std::tuple<float,float>> &gamma_energies_and_yields );

  InterSpec_API std::vector<float> get_characteristic_gammas( const std::string &nuclide );
  
  /** Returns the Nuclides, x-rays, and reaction gammas from CharacteristicGammas.txt, in the energy range. */
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_nuclides_with_characteristics_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec );
  
  /** Gets the approximate FWHM of the foreground spectrum for the provided energy, and returns
   `get_nuclides_with_characteristics_in_energy_range(energy-fwhm, energy+fwhm)`
   
   Throws exception on error (no foreground, invalid InterSpec pointer, etc).
   */
  std::vector<std::variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>>
  get_characteristics_near_energy( const double energy, InterSpec *interspec );

  /** Get the expected Full Width at Half Maximum (FWHM) for a peak at the specified energy.
   
   Uses the detector response function if available, otherwise fits detected peaks to estimate
   resolution, or falls back to detector type-based estimation.
   
   @param energy The energy (in keV) to calculate the expected FWHM for
   @param interspec Pointer to the InterSpec session
   @return The expected FWHM in keV
   @throws std::runtime_error if no InterSpec session, no foreground loaded, or FWHM cannot be determined
   */
  float get_expected_fwhm( const double energy, InterSpec *interspec );
  
  struct SpectrumCountsInEnergyRange
  {
    double lower_energy;
    double upper_energy;

    double foreground_counts;

    /** Counts divided by live time.  If live time of spectrum is not valid, this will be set to NaN. */
    double foreground_cps;

    struct CountsWithComparisonToForeground
    {
      double counts;

      /** Counts divided by live time.  If live time of spectrum is not valid, this will be set to NaN. */
      double cps;

      /** The number of statistical sigma the foreground is is elevated, 
       * relative to this spectrum (the background or secondary), when live-time normalized. 
       * Positive numbers indicate foreground is elevated, negative numbers indicate 
       * deficit in foreground.
       * */
      double num_sigma_rel_foreground;
    };//struct CountsWithComparisonToForeground

    std::optional<CountsWithComparisonToForeground> background_info;
    std::optional<CountsWithComparisonToForeground> secondary_info;
  };//struct SpectrumCountsInEnergyRange

  SpectrumCountsInEnergyRange get_counts_in_energy_range( double lower_energy, double upper_energy, InterSpec *interspec );

  /** Enum specifying the action to perform when editing a peak. */
  enum class EditPeakAction
  {
    SetEnergy,
    SetFwhm,
    SetAmplitude,
    SetEnergyUncertainty,
    SetFwhmUncertainty,
    SetAmplitudeUncertainty,
    SetRoiLower,
    SetRoiUpper,
    SetSkewType,
    SetContinuumType,
    SetSource,
    SetColor,
    SetUserLabel,
    SetUseForEnergyCalibration,
    SetUseForShieldingSourceFit,
    SetUseForManualRelEff,
    DeletePeak,
    SplitFromRoi,
    MergeWithLeft,
    MergeWithRight
  };//enum class EditPeakAction

  /** Converts EditPeakAction enum to string representation.

   @param action The edit action to convert
   @return String representation of the action
   */
  const char* to_string( EditPeakAction action );

  /** Converts string to EditPeakAction enum.

   @param str String representation of the action
   @return The corresponding EditPeakAction enum value
   @throws std::runtime_error if string does not match a valid action
   */
  EditPeakAction edit_peak_action_from_string( const std::string &str );

  /** Options for editing a peak at a specified energy. */
  struct EditAnalysisPeakOptions
  {
    double energy;
    EditPeakAction editAction;
    SpecUtils::SpectrumType specType;
    std::optional<double> doubleValue;
    std::optional<std::string> stringValue;
    std::optional<bool> boolValue;
    std::optional<double> uncertainty;
  };//struct EditAnalysisPeakOptions

  /** Results of peak edit operation. */
  struct EditAnalysisPeakStatus
  {
    bool success;
    std::string message;
    std::optional<std::shared_ptr<const PeakDef>> modifiedPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;
  };//struct EditAnalysisPeakStatus

  /** Edit or delete a peak at the specified energy.

   This function allows modifying peak properties (energy, FWHM, amplitude),
   ROI bounds, continuum/skew types, assigning sources, setting colors/labels,
   deleting peaks, or splitting/merging ROIs.

   @param options Specifies the energy, action, and values for the edit
   @param interspec Pointer to the InterSpec session containing the peak
   @return EditAnalysisPeakStatus containing success/failure and modified peak info
   @throws std::runtime_error if InterSpec is null or other errors occur
   */
  InterSpec_API EditAnalysisPeakStatus edit_analysis_peak( const EditAnalysisPeakOptions &options, InterSpec *interspec );

  /** Enum specifying the type of escape peak. */
  enum class EscapePeakType
  {
    SingleEscape,  // Parent energy minus 511 keV (one positron escapes)
    DoubleEscape   // Parent energy minus 1022 keV (both positrons escape)
  };//enum class EscapePeakType

  /** Converts EscapePeakType enum to string representation.

   @param type The escape peak type to convert
   @return String representation of the type
   */
  const char* to_string( EscapePeakType type );

  /** Information about the parent peak if the queried energy is an escape peak. */
  struct ParentPeakInfo
  {
    double parentPeakEnergy;
    EscapePeakType escapeType;
    std::optional<std::string> sourceLabel;  // e.g., "Th232 S.E. 2614.52 keV" or "Co56 D.E. 3451.15 keV"
    std::optional<std::string> userLabel;    // e.g., "2614.53 S.E." if no source assigned
  };//struct ParentPeakInfo

  /** Options for checking if a peak is an escape peak. */
  struct EscapePeakCheckOptions
  {
    double energy;
    SpecUtils::SpectrumType specType;  // defaults to Foreground if not specified
  };//struct EscapePeakCheckOptions

  /** Results of escape peak analysis. */
  struct EscapePeakCheckStatus
  {
    // If the specified peak IS an escape peak of another peak in the spectrum
    std::optional<ParentPeakInfo> parentPeak;

    // Energies where escape peaks COULD be found (calculated independent of actual peaks)
    double potentialSingleEscapePeakEnergy;     // input - 511 keV
    double potentialDoubleEscapePeakEnergy;     // input - 1022 keV
    double potentialParentPeakSingleEscape;     // input + 511 keV
    double potentialParentPeakDoubleEscape;     // input + 1022 keV

    // Search windows (±keV) for finding peaks at those energies
    // Calculated as max(1.0, min(0.5*fwhm, 15.0)) - fairly arbitrary and untested
    double singleEscapeSearchWindow;            // window for finding S.E. peak of input energy
    double doubleEscapeSearchWindow;            // window for finding D.E. peak of input energy
    double singleEscapeParentSearchWindow;      // window for finding parent if input is S.E.
    double doubleEscapeParentSearchWindow;      // window for finding parent if input is D.E.
  };//struct EscapePeakCheckStatus

  /** Check if a peak at a specified energy is a single or double escape peak.

   This function analyzes whether the peak at the given energy is an escape peak of another
   peak in the spectrum. It checks for parent peaks at energy+511 keV (single escape) and
   energy+1022 keV (double escape), verifying that parent peaks are above the pair production
   threshold (~1255 keV for HPGe). The function also handles the edge case where the parent
   peak is above the detector's energy range by checking nuclide photon energies.

   Note: The function checks for double escape parents (1022 keV above) BEFORE single escape
   parents (511 keV above) to avoid misidentifying a double escape peak as a single escape.

   @param options Specifies the energy to check and which spectrum to use
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return EscapePeakCheckStatus containing parent peak info (if found) and calculated potential energies
   @throws std::runtime_error if InterSpec is null, no spectrum loaded, or other errors occur
   */
  InterSpec_API EscapePeakCheckStatus escape_peak_check( const EscapePeakCheckOptions &options, InterSpec *interspec );

  /** Enum specifying the type of sum peak. */
  enum class SumPeakType
  {
    NotASumPeak,   // Not identified as a sum peak
    RandomSum,     // Sum from random coincidence (pile-up) due to high dead time
    CascadeSum,    // Sum from cascade decay photons emitted nearly simultaneously
    Unknown        // Could be a sum peak but type cannot be determined
  };//enum class SumPeakType

  /** Converts SumPeakType enum to string representation.

   @param type The sum peak type to convert
   @return String representation of the type
   */
  const char* to_string( SumPeakType type );

  /** Information about a sum peak if the queried energy is identified as one. */
  struct SumPeakInfo
  {
    SumPeakType sumType;
    std::shared_ptr<const PeakDef> firstPeak;   // First contributing peak
    std::shared_ptr<const PeakDef> secondPeak;  // Second contributing peak
    std::optional<std::string> userLabel;       // e.g., "Sum 661.66 + 661.66 keV" or "Co60 cascade-sum 1173.23+1332.49 keV"
    std::optional<double> coincidenceFraction;  // For cascade sums: fraction of time photons are emitted together
  };//struct SumPeakInfo

  /** Options for checking if a peak is a sum peak. */
  struct SumPeakCheckOptions
  {
    double energy;
    SpecUtils::SpectrumType specType;  // defaults to Foreground if not specified
    std::optional<double> distance;     // Optional: source-to-detector distance (in PhysicalUnits, e.g., PhysicalUnits::cm). If > 15 cm, cascade-sums are not considered
  };//struct SumPeakCheckOptions

  /** Results of sum peak analysis. */
  struct SumPeakCheckStatus
  {
    // If the specified peak IS a sum peak
    std::optional<SumPeakInfo> sumPeakInfo;

    // Energy windows used for searching for contributing peaks
    double searchWindow;  // calculated as max(1.0, min(0.5*fwhm, 15.0))
  };//struct SumPeakCheckStatus

  /** Check if a peak at a specified energy is a sum peak (random-sum or cascade-sum).

   This function analyzes whether the peak at the given energy is the result of summing
   two lower-energy photons. There are two types of sum peaks:

   1. Random-sum peaks: Occur when unrelated photons arrive at the detector nearly
      simultaneously, typically when dead time is high (>20% for HPGe, >10% for low-res).
      The function validates by checking if the highest-amplitude peaks produce expected
      sum peaks.

   2. Cascade-sum peaks: Occur when a nuclide emits two photons in quick succession and
      the detector is close enough (<~15 cm) to detect both. Both contributing peaks must
      have the same nuclide source assigned, and the nuclide's decay scheme must list
      those gammas as cascade coincidences. If a distance parameter is specified and is
      greater than 15 cm, cascade-sum detection is disabled.

   When multiple candidate peak pairs are found, the function selects the best match by:
   - First prioritizing the closest energy match to the target
   - If multiple pairs have similar energies (within ~0.5 keV), selecting the pair with
     the largest combined area (for cascade sums: area1 × area2 × coincidenceFraction)

   @param options Specifies the energy to check, which spectrum to use, and optionally
                  the source-to-detector distance (in PhysicalUnits)
   @param interspec Pointer to the InterSpec session containing the spectrum
   @return SumPeakCheckStatus containing sum peak info (if identified) and search parameters
   @throws std::runtime_error if InterSpec is null, no spectrum loaded, or other errors occur
   */
  InterSpec_API SumPeakCheckStatus sum_peak_check( const SumPeakCheckOptions &options, InterSpec *interspec );

} // namespace AnalystChecks

#endif // AnalystChecks_h 
