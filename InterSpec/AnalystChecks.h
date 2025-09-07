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

#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <optional>
#include <functional>

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

namespace RelActCalcAuto
{
  struct RelActAutoSolution;
}//namespace RelActCalcAuto

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
  };
  
  struct GetUserPeakStatus {
    std::string userSession;
    std::vector<std::shared_ptr<const PeakDef>> peaks;
  };
  
  InterSpec_API GetUserPeakStatus get_user_peaks( const GetUserPeakOptions &options, InterSpec *interspec );

  struct FitPeaksForNuclideOptions {
    std::vector<std::string> nuclides;
    bool doNotAddPeaksToUserSession = true;
    std::optional<std::string> userSession;
    bool computeAsync = false;
    double numSigmaSignifigance = 1.5;
  };
  
  enum class FitPeaksForNuclideStatus
  {
    Success,
    FailureSettingUp,
    FailedComputation
  };
  
  struct FitPeaksForNuclideResult {
    FitPeaksForNuclideStatus status;
    std::string error_message;

    std::shared_ptr<RelActCalcAuto::RelActAutoSolution> solution;
    std::vector<std::shared_ptr<const PeakDef>> fitPeaks;
  };
  
  InterSpec_API void fit_peaks_for_nuclides( const FitPeaksForNuclideOptions &options,
                                             InterSpec *interspec,
                                             std::function<void(const FitPeaksForNuclideResult &result)> callback );
  
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
} // namespace AnalystChecks

#endif // AnalystChecks_h 
