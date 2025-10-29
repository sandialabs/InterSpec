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
  /** A simple "hello world" type function for testing purposes.
   
   Returns a greeting message indicating the AnalystChecks module is working.
   
   @return A string containing a hello world message.
   */
  InterSpec_API std::string hello_world();
  
  /** Options for peak detection analysis. */
  struct DetectedPeaksOptions {
    SpecUtils::SpectrumType specType;
    std::optional<std::string> userSession;
  };
  
  /** Results of peak detection analysis. */
  struct DetectedPeakStatus {
    std::string userSession;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
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
    bool addToUsersPeaks;
    double energy;
    SpecUtils::SpectrumType specType;
    std::optional<std::string> source;
    std::optional<std::string> userSession;
  };
  
  struct FitPeakStatus {
    std::string userSession;
    std::shared_ptr<const PeakDef> fitPeak;
    std::vector<std::shared_ptr<const PeakDef>> peaksInRoi;
  };
  
  InterSpec_API FitPeakStatus fit_user_peak( const FitPeakOptions &options, InterSpec *interspec );
  
  struct GetUserPeakOptions {
    SpecUtils::SpectrumType specType;
    std::optional<std::string> userSession;
    std::optional<double> lowerEnergy;
    std::optional<double> upperEnergy;
  };
  
  struct GetUserPeakStatus {
    std::string userSession;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  InterSpec_API GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec );

  struct FitPeaksForNuclideOptions {
    std::vector<std::string> nuclides;
    bool doNotAddPeaksToUserSession;
  };
  
  struct FitPeaksForNuclideStatus {
    std::vector<std::shared_ptr<const PeakDef>> fitPeaks;
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
    std::optional<std::string> userSession;
  };//struct EditAnalysisPeakOptions

  /** Results of peak edit operation. */
  struct EditAnalysisPeakStatus
  {
    std::string userSession;
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

} // namespace AnalystChecks

#endif // AnalystChecks_h 
