#ifndef PeakFitDetPrefs_h
#define PeakFitDetPrefs_h
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

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <optional>

#include <Wt/WFlags>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"

// Forward declarations
namespace rapidxml
{
  template<class Ch> class xml_node;
  template<class Ch> class xml_document;
}//namespace rapidxml

namespace SpecUtils
{
  class Measurement;
}

class DetectorPeakResponse;


/** Holds coarse detector resolution type, default peak skew type, and optional
 fixed skew parameters for peak fitting.

 The skew parameters (SkewPar0..SkewPar3) may be energy-dependent or fixed:
 - If a parameter is energy-dependent (see `PeakDef::is_energy_dependent`), both
   lower and upper energy values must be specified. The parameter value will be
   linearly interpolated between these for any given peak energy.
 - If a parameter is NOT energy-dependent, only the lower value is used.
 - If a parameter is not specified (nullopt), it will be fit per-ROI.
 */
struct PeakFitDetPrefs
{
  /** The coarse detector type to assume when fitting peaks - affects the FWHM
   limits of peaks that could be fit.
   */
  PeakFitUtils::CoarseResolutionType m_det_type;

  /** Skew type to default to fitting new peaks with. */
  PeakDef::SkewType m_peak_skew_type;

  /** Lower-energy skew parameter values (index 0..3 => SkewPar0..SkewPar3).
   See class-level comment for semantics of energy-dependent vs fixed params.
   */
  std::optional<double> m_lower_energy_skew[4];

  /** Upper-energy skew parameter values (index 0..3 => SkewPar0..SkewPar3).
   Only meaningful for energy-dependent parameters.
   */
  std::optional<double> m_upper_energy_skew[4];

  /** If true, each ROI will have its own independently fit skew parameters.
   If false, skew parameters are shared/related across ROIs.
   */
  bool m_roi_independent_skew;

  /** How peak FWHM (sigma) is determined during fitting. */
  enum class FwhmMethod : int
  {
    /** FWHM is freely fit from the data (default, current behavior). */
    Normal,

    /** FWHM is fixed to the DRF prediction for the peak energy.
     Requires the DetectorPeakResponse to have resolution info.
     Falls back to Normal if DRF lacks FWHM info.
     */
    DetFwhm,

    /** FWHM starts at the DRF prediction but is allowed to refine within
     a small range during fitting.
     Requires the DetectorPeakResponse to have resolution info.
     Falls back to Normal if DRF lacks FWHM info.
     */
    DetPlusRefine
  };//enum class FwhmMethod

  /** How peak FWHM is determined during fitting; default Normal. */
  FwhmMethod m_fwhm_method;

  enum class LoadingSource : int
  {
    /** No specific source; using application default. */
    Default,

    /** User explicitly set these values in the GUI. */
    UserInputInGui,

    /** Loaded from the DetectorPeakResponse. */
    FromDetectorPeakResponse,

    /** Derived from a hard-coded mapping of SpecUtils::DetectorType. */
    DefaultForDetectorType,

    /** Guessed from spectral data features (future). */
    FromSpectralData
  };//enum class LoadingSource

  LoadingSource m_source;


  PeakFitDetPrefs();

  bool operator==( const PeakFitDetPrefs &rhs ) const;
  bool operator!=( const PeakFitDetPrefs &rhs ) const;


  // -- String conversion for CoarseResolutionType --

  /** Returns a short string identifier: "Low", "LaBr", "CZT", "High", "Unknown". */
  static const char *to_str( PeakFitUtils::CoarseResolutionType type );

  /** Parses a CoarseResolutionType from string (case-insensitive).
   Throws std::exception on unrecognized input.
   */
  static PeakFitUtils::CoarseResolutionType coarse_res_from_str( const std::string &str );


  // -- String conversion for FwhmMethod --

  static const char *to_str( FwhmMethod method );
  static FwhmMethod fwhm_method_from_str( const std::string &str );

  // -- String conversion for LoadingSource --

  static const char *to_str( LoadingSource src );
  static LoadingSource loading_source_from_str( const std::string &str );


  // -- XML serialization --

  static const int sm_xmlSerializationMajorVersion; // = 0
  static const int sm_xmlSerializationMinorVersion; // = 0

  /** Appends a `<PeakFitDetPrefs>` child node to `parent_node`. */
  void toXml( ::rapidxml::xml_node<char> *parent_node,
              ::rapidxml::xml_document<char> *doc ) const;

  /** Reads from a `<PeakFitDetPrefs>` XML node. Throws on error. */
  void fromXml( const ::rapidxml::xml_node<char> *node );


  // -- URL serialization (for embedding in DetectorPeakResponse URL) --

  /** Returns "&"-separated key=value pairs for embedding in a URL query string.
   Returns empty string if all defaults.
   Keys: DT (det type), SK (skew type), LS0..LS3 / US0..US3 (skew params).
   */
  std::string toUrlQueryParts() const;

  /** Parses from key=value pairs previously extracted from a URL query.
   Throws on error.
   */
  void fromUrlQueryParts( const std::map<std::string,std::string> &parts );


  // -- Factory methods --

  /** Returns a default PeakFitDetPrefs for a given SpecUtils::DetectorType,
   optionally checking keywords in manufacturer/model strings.
   Returns nullptr if no default mapping exists.
   */
  static std::shared_ptr<const PeakFitDetPrefs>
    defaultForDetectorType( const int specutils_det_type,
                            const std::string &manufacturer,
                            const std::string &model );

  /** Stub for future: guess PeakFitDetPrefs from spectral data features.
   Currently returns nullptr.
   */
  static std::shared_ptr<const PeakFitDetPrefs>
    guessFromSpectralData( const std::shared_ptr<const SpecUtils::Measurement> &meas );
};//struct PeakFitDetPrefs


/** Applies the FwhmMethod from prefs to peaks before fitting.

 For FwhmMethod::DetFwhm: sets each peak's sigma to the DRF prediction and
   sets fitFor(Sigma) = false so sigma is held constant during the fit.
 For FwhmMethod::DetPlusRefine: sets each peak's sigma to the DRF prediction
   and adds SmallFwhmRefinementOnly to fit_options so sigma can refine within +-15%.
 For FwhmMethod::Normal: no-op.

 Silently falls back to Normal if drf is null or lacks resolution info.
 */
void apply_fwhm_method_to_peaks(
  std::vector<std::shared_ptr<PeakDef>> &peaks,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const PeakFitDetPrefs &prefs,
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> &fit_options );


#endif //PeakFitDetPrefs_h
