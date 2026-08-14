#ifndef ExportSpecFileCAM_h
#define ExportSpecFileCAM_h
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
#include <set>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <deque>
#include <optional>

#include <Wt/WSignal.h>
#include <Wt/WModelIndex.h>
#include <Wt/WContainerWidget.h>
#include <Wt/WAbstractItemModel.h>

#include "SpecUtils/CAMIO.h"

class SpecMeas;
class PeakDef;
class InterSpec;
class RowStretchTreeView;
class NativeFloatSpinBox;
class DetectorPeakResponse;

namespace Wt
{
  class WText;
  class WLabel;
  class WLineEdit;
  class WCheckBox;
  class WComboBox;
  class WPushButton;
  class WContainerWidget;
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
}

namespace SpecUtils
{
  struct EnergyCalibration;
}

/** Logic (independent of the GUI) for building a GENIE nuclide library, FWHM shape
 calibration, and efficiency curve to write into an exported CNF file.

 See `ExportSpecFile.cpp` for the "Spectrum File Export" dialog options that drive this,
 and `CAMInputOutput::CnfGenieExtras` for what actually gets written into the file.
 */
namespace ExportSpecFileCAM
{
  /** Whether the nuclide library should only contain the exact lines the foreground's peaks
   are assigned to, or all lines (above some yield threshold) of every identified source.
   */
  enum class GenieLibraryLineMode
  {
    PeakLinesOnly,
    AllLinesAboveThreshold
  };//enum class GenieLibraryLineMode


  /** One candidate gamma or x-ray line for a `GenieLibrarySource`. */
  struct GenieLibraryLine
  {
    /** Energy, in keV. */
    float energy = 0.0f;
    /** Negative if not known (Genie will estimate an uncertainty from the value itself). */
    float energy_uncert = -1.0f;

    /** Yield, in photons/xrays per decay of the parent nuclide - a fraction, not a percent (it is
     converted to the percent GENIE wants in `to_library_lines(...)`).

     Usually 0 to 1, but can exceed 1: annihilation lines get a contribution from every positron,
     so e.g. Na22's 511 keV yield is about 1.8.
     */
    float yield = 0.0f;
    /** Negative if not known (Genie will estimate an uncertainty from the value itself). */
    float yield_uncert = -1.0f;

    bool is_xray = false;

    /** Set once, after all sources/lines have been built and merged; matches (as closely as
     practical, for display purposes) what `CAMInputOutput::CAMIO::AssignKeyLines()` will
     independently compute when the file is actually written.
     */
    bool is_key_line = false;

    /** True when no peak was fit to this line, and it is only in the library because it falls
     within `sm_interference_num_fwhm` FWHM of a line another source *does* have a peak at - so
     Genie can correct that peak for this source's contribution.  See `build_genie_library(...)`.
     */
    bool is_interference = false;

    /** Whether this line should be written to the library; defaults to true (checked). */
    bool included = true;
  };//struct GenieLibraryLine


  /** One source (nuclide) and its candidate library lines. */
  struct GenieLibrarySource
  {
    /** The nuclide this source represents; kept (in addition to `name`) so callers can key a
     `nuclide_ages` map (see `build_genie_library(...)`) when the user changes this source's age.
     */
    const SandiaDecay::Nuclide *nuclide = nullptr;

    /** E.g. "Co-60". */
    std::string name;

    float half_life_seconds = 0.0f;
    /** Negative if not known (Genie will estimate an uncertainty from the value itself). */
    float half_life_uncert_seconds = -1.0f;

    /** Whether the user can usefully vary this nuclide's age (i.e.,
     `!PeakDef::ageFitNotAllowed(nuclide)`); if false, the GUI shouldn't offer an age input for
     this source - but `age_seconds` may still be non-zero (see `age_seconds`).
     */
    bool is_ageable = false;

    /** The age (in seconds) used to compute `lines`' yields; always defaults to
     `PeakDef::defaultDecayTime(nuclide)`, even when `!is_ageable` - e.g. Cs137's spectrum
     doesn't meaningfully change with further aging once its short-lived Ba137m daughter
     reaches secular equilibrium, but an age of exactly 0 would (incorrectly) give its 661.7 keV
     line zero yield, since that daughter population hasn't grown in yet at t=0.  Changing this
     requires re-calling `build_genie_library(...)` (passing the desired age in `nuclide_ages`,
     only honored for `is_ageable` sources) to recompute `lines`' yields.
     */
    double age_seconds = 0.0;

    std::vector<GenieLibraryLine> lines;

    /** Whether this source (and by default, all its lines) should be written to the
     library; defaults to true (checked). Un-checking this un-checks all of `lines`.
     */
    bool included = true;
  };//struct GenieLibrarySource


  /** Builds the candidate GENIE library from the given peaks.

   Peaks assigned to a `ReactionGamma::Reaction`, or to an x-ray with no associated parent
   nuclide (i.e., no half-life is available to report), cannot be represented in a Genie
   nuclide library and are skipped (a warning is added to `warnings`, if provided).

   @param peaks The fitted peaks (typically the foreground's peaks) to build the library from.
   @param mode Whether to only include the exact lines the peaks are assigned to, or all lines
          of each identified source above `yield_threshold_percent`.
   @param yield_threshold_percent Only used when `mode == AllLinesAboveThreshold`; only lines
          with yield at least this percent of the source's most intense line are included.
   @param combine_unresolvable_lines If true, lines a detector with `fwhm_coeffs` resolution
          could not tell apart (within a source) are merged into a single yield-weighted line,
          mirroring how the GENIE Library Editor (and the "Peak-Map" tool) combine unresolvable
          lines.  See `sm_cluster_num_sigma` for the clustering criterion.
   @param include_interference_lines Only used when `mode == PeakLinesOnly`.  If true, each source
          additionally gets the lines it emits within `sm_interference_num_fwhm` FWHM of *another*
          source's peak-assigned line, flagged `GenieLibraryLine::is_interference`.  Without them a
          library says only one nuclide emits at that energy, and Genie attributes the whole peak to
          it rather than splitting it between the two.  Only nuclides that already have a peak of
          their own become sources, so this never adds a nuclide the user did not identify, and the
          added lines are subject to `yield_threshold_percent` (relative to their own source's most
          intense line) so a chain's hundreds of trace lines do not flood the library.
   @param fwhm_coeffs The `{A0, A1}` of `FWHM = A0 + A1*sqrt(energy)` used to decide which lines
          are unresolvable, and how far an interference line may be from a peak's line; this should
          be the same shape calibration being written into the CNF file, so the library matches what
          Genie itself will be able to resolve.  Only used when `combine_unresolvable_lines` or
          `include_interference_lines` is true.
   @param energy_range The `{lower, upper}` energy (keV) the spectrum covers; lines outside it are
          left out, since Genie could never match a peak to them.  Pass `{0,0}` to disable.  Only
          used when `mode == AllLinesAboveThreshold`.
   @param nuclide_ages Optional per-nuclide age overrides (in seconds) to use when computing
          line yields, honored only for nuclides where `!PeakDef::ageFitNotAllowed(nuclide)`; a
          nuclide not present in this map (or where `ageFitNotAllowed(nuclide)`) uses
          `PeakDef::defaultDecayTime(nuclide)`.  Re-call this function (with an updated age here)
          to recompute a source's line yields after the user changes its age in the GUI.
   @param warnings If non-null, human-readable warnings are appended here (e.g., for peaks
          that could not be represented in the library).
   @returns One entry per distinct source (nuclide) found in `peaks`, each with candidate lines,
            all defaulting to `included = true`.
   */
  std::vector<GenieLibrarySource> build_genie_library(
                                  const std::deque<std::shared_ptr<const PeakDef>> &peaks,
                                  const GenieLibraryLineMode mode,
                                  const double yield_threshold_percent,
                                  const bool combine_unresolvable_lines,
                                  const bool include_interference_lines,
                                  const std::pair<float,float> &fwhm_coeffs,
                                  const std::pair<float,float> &energy_range,
                                  const std::map<const SandiaDecay::Nuclide *, double> &nuclide_ages = {},
                                  std::vector<std::string> *warnings = nullptr );

  /** How many peak sigma apart two gamma lines have to be before they are treated as separately
   resolvable; matches the `photopeak_cluster_sigma` default that
   `ShieldingSourceFitCalc::ShieldingSourceFitOptions` (and hence the Activity/Shielding fit and
   the manual rel-eff calculations) use to decide which lines contribute to the same peak.
   */
  inline constexpr double sm_cluster_num_sigma = 1.25;

  /** How many FWHM from a peak-assigned line another source's line has to be within before it is
   written into the library as an interference line; see `build_genie_library(...)`.

   Deliberately wider than `sm_cluster_num_sigma` - lines that far apart *are* separately
   resolvable, but they still overlap enough to bias each other's areas, which is exactly the case
   Genie's interference correction exists for.
   */
  inline constexpr double sm_interference_num_fwhm = 1.5;

  /** Flattens the checked (`included`) sources/lines of `sources` into the line-list that
   `SpecUtils::SpecFile::write_cnf(...)`'s `CnfGenieExtras` expects.
   */
  std::vector<CAMInputOutput::CnfGenieExtras::LibraryLine>
      to_library_lines( const std::vector<GenieLibrarySource> &sources );


  /** Fits `FWHM = a + b*sqrt(energy)` (Genie's FWHM equation form) to the given DRF's FWHM
   function, sampled log-spaced across its valid energy range (or 59-2614 keV, if the DRF doesn't
   specify one), minimizing *relative* error.

   Throws if the DRF's FWHM cannot be fit to this form.

   @returns {a, b}
   */
  std::pair<float,float> fit_genie_fwhm_from_drf( const DetectorPeakResponse &drf );

  /** A default `{FWHMOFF, FWHMSLOPE}` for `FWHM = FWHMOFF + FWHMSLOPE*sqrt(energy)`, for when
   no better information is available.

   For HPGe this is the same Ge default `CAMIO::AddDetectorType(...)` writes.  For everything
   else it is a least-squares fit of the Genie equation form to `PeakFitUtils::nai_fwhm_fcn(...)`
   (InterSpec's nominal "NaI 3x3" resolution), which tracks a real NaI detector much better than
   the Genie manual's own NaI defaults do at low energy.
   */
  std::pair<float,float> default_genie_fwhm( const bool is_hpge );


  /** The result of converting a `DetectorPeakResponse`'s efficiency function into the closest
   available Genie efficiency model.

   Genie derives its own efficiency fit from the calibration points, weighting each by
   `1/sigma_rel^2`, so `points` always carries a non-zero uncertainty (see
   `genie_default_eff_uncertainty(...)`) - with zero uncertainties that fit is undefined and Genie
   can only offer its "Interpolated" model.  `fit_coeffs` is that same fit, written into the file
   so the curve is there whether or not Genie recomputes it; see `CAMIO::AddEfficiencyFit(...)`.
   */
  struct GenieEfficiencyResult
  {
    CAMInputOutput::CAMIO::EfficiencyModel model = CAMInputOutput::CAMIO::EfficiencyModel::Unknown;
    std::vector<CAMInputOutput::EfficiencyPoint> points;

    /** Coefficients of `ln(eff) = sum_i{ fit_coeffs[i] * ln(energy/keV)^i }`; empty when the
     points could not be fit well, in which case `model` is `INTERPOL`.
     */
    std::vector<float> fit_coeffs;
    /** `E0` the fit's Empirical display form is expressed about; see
     `CAMIO::sm_geom_default_fit_ref_energy`.
     */
    float fit_reference_energy = 0.0f;
    /** Reduced chi-square of that fit. */
    float fit_chi_square = 1.0f;
    /** Name of the detector response the curve came from; GENIE keeps a description alongside its
     own efficiency calibrations.
     */
    std::string detector_name;
  };//struct GenieEfficiencyResult


  /** The relative efficiency uncertainty to report for a calibration point at `energy` (keV),
   when the detector response function does not carry one of its own (none currently do).

   Genie weights its efficiency fit by `1/sigma_rel^2`, so these do more than document precision -
   they set the shape of the curve Genie fits, and a zero would make that fit undefined.  The
   steps track what Genie itself carries in a real efficiency calibration: 15% below 55 keV, then
   10%, 8%, 6% and 4% as energy rises.
   */
  float genie_default_eff_uncertainty( const float energy );


  /** Fits `ln(eff) = sum_i{ coeffs[i]*ln(energy)^i }` to efficiency calibration points, weighted
   by `1/rel_uncerts^2` - the curve GENIE derives from the points of its own GEOM block, and what
   `CAMIO::AddEfficiencyFit(...)` writes.

   @param energies Point energies, in keV.
   @param efficiencies Absolute efficiencies at those energies; all must be > 0.
   @param rel_uncerts Relative uncertainty of each efficiency; pass empty for an unweighted fit.
   @param order Polynomial order; there must be at least `order+1` points.
   */
  std::vector<double> fit_genie_efficiency_curve( const std::vector<double> &energies,
                                                  const std::vector<double> &efficiencies,
                                                  const std::vector<double> &rel_uncerts,
                                                  const size_t order );

  /** The GENIE low-tail parameter `T` (in keV) equivalent to a peak's skew, or a negative value
   if the peak has no skew this maps onto.

   GENIE's low tail is the same function as InterSpec's `GaussExp` skew, with `T = skew*sigma`;
   `ExpGaussExp`'s lower-side parameter uses the same convention.  Other skew forms are different
   functions and are reported as "no tail" rather than approximated.
   */
  double genie_low_tail( const PeakDef &peak );

  /** Fits GENIE's low-tail calibration `T(E) = B2 + B3*E` to whichever of `peaks` have a tail
   (see `genie_low_tail(...)`); returns nothing when none do, in which case no low-tail
   calibration should be written.
   */
  std::optional<std::pair<float,float>> fit_genie_low_tail_cal(
                              const std::deque<std::shared_ptr<const PeakDef>> &peaks );


  /** Converts InterSpec's fitted peaks into the `CAMInputOutput::Peak` records
   `SpecUtils::SpecFile::write_cnf(...)` writes into a CNF file's PEAK block.

   Only Gaussian peaks are converted - a "data defined" peak has no centroid or width to report.
   Channel numbers and the centroid come from `energy_cal`; if it is null or invalid, the
   channel-based fields are left zero.

   Note: see `CAMIO::AddPeak(...)` - real GENIE reads the regions of interest out of files written
   this way, but what the per-peak values look like in its own peak report has not been checked.

   @param peaks The fitted peaks to convert.
   @param energy_cal The energy calibration of the spectrum being written, used to convert
          energies to channel numbers.
   @param live_time_s Live time in seconds, used for the count-rate fields; count rates are left
          zero if this is not positive.
   @param efficiency The efficiency curve being written into the same file, if any; each peak
          records the curve's value at its energy, exactly as GENIE does.  Pass null when no
          efficiency is being written, and those fields are left zero - which is also what a
          GENIE file without an efficiency calibration has.
   @param data The spectrum the peaks were fit to.  A stepped continuum's area cannot be computed
          without it (`PeakContinuum::offset_integral(...)` throws), so ROIs with a stepped
          continuum are written with a continuum area of zero when this is null.
   @param warnings If non-null, human-readable warnings are appended here - notably when a peak
          falls outside the efficiency curve's energy range, where its efficiency is written as
          zero and GENIE would divide the peak area by it.
   */
  std::vector<CAMInputOutput::Peak> to_cam_peaks(
              const std::deque<std::shared_ptr<const PeakDef>> &peaks,
              const std::shared_ptr<const SpecUtils::EnergyCalibration> &energy_cal,
              const float live_time_s,
              const GenieEfficiencyResult * const efficiency = nullptr,
              const std::shared_ptr<const SpecUtils::Measurement> &data = nullptr,
              std::vector<std::string> *warnings = nullptr );


  /** Evaluates the efficiency curve `convert_efficiency_to_genie(...)` produced, at `energy` (keV):
   from its fitted coefficients when it has them, otherwise by interpolating its points in log-log
   space (which is what GENIE's "Interpolated" model does).  Returns 0 outside the points' range,
   or if there is nothing to evaluate.
   */
  double genie_efficiency_at( const GenieEfficiencyResult &efficiency, const double energy );


  /** Converts a DetectorPeakResponse's efficiency into the closest available Genie efficiency
   model, writing out (energy, efficiency, uncertainty) points in all cases, plus the fitted
   curve when the points support one:
   - `kEnergyEfficiencyPairs` -> the DRF's own tabulated point energies.
   - `kExpOfLogPowerSeries` and `kFunctialEfficienyForm` -> points sampled log-spaced across the
     DRF's valid range.
   In every case a weighted ln(eff)-vs-ln(energy) polynomial is fit to the points; if it describes
   them well the model is tagged `EMPIRICAL` and its coefficients are written, otherwise the model
   is `INTERPOL` and Genie interpolates between the points.

   Genie's efficiency curve is an *absolute* full-energy-peak efficiency - counts recorded per
   gamma emitted by the source, at the geometry the calibration was taken at - not InterSpec's
   intrinsic efficiency (per gamma striking the detector face).  So for a far-field DRF the
   points written are `drf.efficiency(energy, distance)`, and the caller must supply the source
   distance the resulting CNF file is meant to describe.  For a fixed-geometry DRF the intrinsic
   efficiency already *is* the absolute efficiency, and `distance` is ignored.

   @param distance Source-to-detector distance (in InterSpec's internal length units) at which to
          evaluate the absolute efficiency; ignored when `drf.isFixedGeometry()`.  Must be > 0 for
          a far-field DRF, or an exception is thrown.
   @param warnings If non-null, human-readable warnings are appended here.
   */
  GenieEfficiencyResult convert_efficiency_to_genie( const DetectorPeakResponse &drf,
                                                     const double distance,
                                                     std::vector<std::string> *warnings = nullptr );


  /** Two-level (source -> line) table model backing the GENIE library preview/selection table
   in the CNF export options: top-level rows are `GenieLibrarySource`s, and their children are
   that source's `GenieLibraryLine`s.  Each row (source or line) has a checkbox in the
   `Column::Include` column; checking/unchecking a source checks/unchecks all its lines.

   Line rows additionally have a `Column::Key` checkbox.  Exactly one included line of each source
   is the key line, so that column behaves like a radio group: checking a line moves the key line
   off whatever line held it, and unchecking the key line is refused.
   */
  class GenieLibraryModel : public Wt::WAbstractItemModel
  {
  public:
    enum class Column : int
    {
      /** Nuclide name (source rows) or energy, in keV, (line rows) - the header carries the unit. */
      Name,
      /** Half-life (source rows) or yield percent (line rows). */
      Info,
      /** Editable age (source rows of ageable nuclides only; blank/unused otherwise). */
      Age,
      /** Key-line checkbox; line rows of included lines only (blank on source rows). */
      Key,
      Include,
      NumColumns
    };//enum class Column

    explicit GenieLibraryModel( Wt::WObject *parent = nullptr );

    /** Replaces the underlying data, e.g. after `build_genie_library(...)` is re-run because the
     user changed the line-selection mode, yield threshold, combine-lines option, or a source's
     age.  Resets the model (all previous Wt::WModelIndex's for this model become invalid).
     */
    void setSources( std::vector<GenieLibrarySource> sources );

    const std::vector<GenieLibrarySource> &sources() const { return m_sources; }

    /** Whether any source's age can actually be edited; if not, the Age column is hidden. */
    bool anySourceIsAgeable() const;

    virtual Wt::WModelIndex index( int row, int column, const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const override;
    virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual Wt::cpp17::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const override;
    virtual bool setData( const Wt::WModelIndex &index, const Wt::cpp17::any &value, int role = Wt::EditRole ) override;
    virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const override;
    virtual Wt::cpp17::any headerData( int section, Wt::Orientation orientation = Wt::Horizontal, int role = Wt::DisplayRole ) const override;

    /** Emitted after a source's age is edited by the user (with the new age, in seconds); the
     owning widget should re-run `build_genie_library(...)` with an updated `nuclide_ages` map
     and call `setSources(...)` again to refresh the line yields.
     */
    Wt::Signal<const SandiaDecay::Nuclide *, double> &ageEdited();

  private:
    /** Makes exactly one *included* line of the given source the key line, keeping the current
     one if it is still included, and otherwise falling back to the highest-yield included line.
     Sources with no included lines are left with no key line.  Returns whether anything changed.
     */
    bool ensureSingleKeyLine( const int source_index );

    std::vector<GenieLibrarySource> m_sources;
    Wt::Signal<const SandiaDecay::Nuclide *, double> m_ageEdited;
  };//class GenieLibraryModel


  /** Which Genie shape-calibration ("FWHM") this dialog will write. */
  enum class GenieFwhmSource
  {
    None,
    DefaultHPGe,
    DefaultNaI,
    FromDrf,
    /** Manually fit from the foreground's peaks, via a `MakeFwhmForDrfWindow` opened by this
     widget; see `fwhmCoefficients()`.
     */
    FromPeaks
  };//enum class GenieFwhmSource


  /** The CNF-format-only export options panel: whether to write a Genie nuclide library, FWHM
   shape calibration, efficiency curve, and/or energy calibration, plus the controls (including
   the `GenieLibraryModel`-backed preview/selection table) needed to configure each.

   Shown/hidden by `ExportSpecFileTool` only when the CNF format is selected; see
   `ExportSpecFile.cpp`.
   */
  class GenieCnfOptionsWidget : public Wt::WContainerWidget
  {
  public:
    explicit GenieCnfOptionsWidget( InterSpec *viewer, Wt::WContainerWidget *parent = nullptr );

    /** Call whenever the selected file/samples change, so the library table (and FWHM/efficiency
     defaults) can be rebuilt/updated from the (possibly new) peaks and detector response.
     */
    /** @param detectors The detector names the export will actually write.  The live time, count
     rates and continuum areas are all computed over these, so passing the wrong set makes those
     numbers disagree with the spectrum in the same file.  Pass empty for all detectors.
     */
    void updateForFile( const std::shared_ptr<const SpecMeas> &spec, const std::set<int> &samples,
                        const std::vector<std::string> &detectors = {} );

    bool writeLibrary() const;
    /** Whether the fitted peaks should be written into the file's PEAK block. */
    bool writePeaks() const;
    bool writeFwhm() const;
    bool writeEfficiency() const;

    /** Whether the energy calibration should be written.

     Always true while a spectrum is being written (the option is hidden then - channel counts
     are not interpretable without it); only then does the user's checkbox decide.
     */
    bool writeEnergyCal() const;

    /** The `{FWHMOFF, FWHMSLOPE}` the current options would write; also what the library
     line-clustering uses, so the library only distinguishes lines Genie could resolve.
     Falls back to `default_genie_fwhm(...)` when the user has turned FWHM writing off.
     */
    std::pair<float,float> currentFwhmCoefficients() const;

    /** The distance entered in the efficiency distance field, in InterSpec's internal length
     units; returns a negative value if the field is empty, hidden, or unparseable.
     */
    double currentEfficiencyDistance() const;

    /** Whether the spectrum itself should be written.  When false the result is a
     nuclide-library/calibration-only CAM file, and the "write energy calibration" option becomes
     meaningful (and is shown).
     */
    bool writeSpectrum() const;

    /** Whether the current combination of options can actually produce a file.

     False when the options as set would write nothing at all, or when an input needed by a
     checked option is missing/invalid (e.g. an unparseable efficiency distance).  `reason`, if
     non-null, is set to a localized explanation.
     */
    bool canExport( Wt::WString *reason = nullptr ) const;

    /** Emitted when `canExport()` may have changed. */
    Wt::Signal<> &exportableChanged();

    /** The `GenieFwhmSource` currently selected in the combo (which only offers usable ones). */
    GenieFwhmSource currentFwhmSource() const;

    /** Builds the `CnfGenieExtras` reflecting the current widget state (checked library
     sources/lines, chosen FWHM source, whether to include efficiency/energy-cal), for the
     spectrum/samples last passed to `updateForFile(...)`.
     */
    CAMInputOutput::CnfGenieExtras currentExtras( std::vector<std::string> *warnings = nullptr ) const;

    /** The warnings writing the efficiency curve would produce for the current options - points
     dropped as negligible, a fit rejected in favour of plain interpolation, or no curve at all.

     Recomputed whenever the efficiency options change so the dialog can show them *before* the
     user exports; the conversion is a couple of dozen samples and one small fit.
     */
    std::vector<std::string> efficiencyWarnings() const;

  protected:
    void rebuildLibraryTable();
    void handleLibraryOptionsChanged();
    void handleSourceAgeEdited( const SandiaDecay::Nuclide *nuc, double age_seconds );
    void handleFwhmSourceChanged();
    void handleFitFwhmFromPeaksClicked();
    void handleFwhmFitFromToolUpdated( std::shared_ptr<DetectorPeakResponse> new_drf );

    /** Shows/hides the efficiency distance field, and validates what is in it.  Also serves to
     force a round-trip, so the checkbox state actually reaches the server before the
     (link-based, and hence non-round-tripping) "Export" button is followed.
     */
    void handleEfficiencyOrEnergyCalChanged();

    void handleDetectorChanged( std::shared_ptr<DetectorPeakResponse> drf );
    void handleSourceExpanded( const Wt::WModelIndex &index );
    void handleSourceCollapsed( const Wt::WModelIndex &index );

    /** Shows/hides the library options that only apply to some of the other options: the
     interference checkbox (peak-lines mode, and only worth offering with more than one source),
     and the yield threshold (which gates both the "all lines" mode and the interference lines).

     Has to run after `rebuildLibraryTable()`, since the source count is only known then.
     */
    void updateLibraryOptionVisibility();

    /** Whether the library being built depends on the FWHM shape calibration - i.e. whether
     changing the FWHM has to trigger a `rebuildLibraryTable()`.
     */
    bool libraryUsesFwhm() const;

    /** Shows/hides the options that only make sense for one of spectrum / library-only output. */
    void handleWriteSpectrumChanged();

    /** Re-populates the FWHM-source combo with only the sources that can actually produce a
     shape calibration right now, and updates the warning text.
     */
    void updateAvailableFwhmSources();

    /** Updates `m_warningsTxt` from the library warnings plus any "this option cannot be
     honored" message, and emits `exportableChanged()`.
     */
    void updateWarningsAndExportable();

    InterSpec *m_interspec;

    /** The spectrum last passed to `updateForFile(...)`.
     Note: a `shared_ptr`, not a `weak_ptr` - `ExportSpecFileTool::currentlySelectedFile()` can
     return a freshly-created `SpecMeas` (e.g. the "foreground + background" combination) that
     nothing else holds a reference to.
     */
    std::shared_ptr<const SpecMeas> m_spec;
    std::set<int> m_samples;

    /** The detectors the export will write; empty means all of them.  See `updateForFile(...)`. */
    std::vector<std::string> m_detectors;
    std::map<const SandiaDecay::Nuclide *, double> m_nuclide_ages;

    /** Whether the current spectrum has an energy calibration that can be written at all. */
    bool m_energy_cal_ok = false;

    /** Names of the sources the user has expanded in the table, so expansion survives a rebuild. */
    std::set<std::string> m_expanded_sources;

    /** Warnings from the last `build_genie_library(...)`, kept so they can be re-combined with
     option-level warnings without rebuilding.
     */
    std::vector<std::string> m_library_warnings;

    /** Which `GenieFwhmSource` each entry of `m_fwhmSourceCb` corresponds to.  The combo only
     offers sources that can actually produce a shape calibration right now, so its index is not
     the enum value.
     */
    std::vector<GenieFwhmSource> m_available_fwhm_sources;

    /** Emitted when whether the current options can actually be exported changes; the owning
     `ExportSpecFileTool` uses this to enable/disable its "Export" button.
     */
    Wt::Signal<> m_exportableChanged;

    Wt::WCheckBox *m_writeSpectrumCb;
    Wt::WCheckBox *m_writePeaksCb;
    Wt::WCheckBox *m_writeLibraryCb;
    Wt::WComboBox *m_libraryModeCb;
    /** Only shown in peak-lines mode, and only when there is more than one source to interfere. */
    Wt::WCheckBox *m_interferenceLinesCb;
    Wt::WLabel *m_thresholdLabel;
    NativeFloatSpinBox *m_thresholdEdit;
    Wt::WCheckBox *m_combineLinesCb;
    RowStretchTreeView *m_libraryTable;
    GenieLibraryModel *m_libraryModel;

    Wt::WCheckBox *m_writeFwhmCb;
    Wt::WComboBox *m_fwhmSourceCb;
    Wt::WPushButton *m_fitFwhmFromPeaksBtn;
    /** Shows the coefficients a completed "From peaks..." fit produced. */
    Wt::WText *m_fwhmFitTxt;
    bool m_have_manual_fwhm;
    std::pair<float,float> m_manual_fwhm_coeffs;

    Wt::WCheckBox *m_writeEfficiencyCb;
    /** Only shown for a non-fixed-geometry DRF, where an absolute efficiency needs a distance. */
    Wt::WLabel *m_effDistanceLabel;
    Wt::WLineEdit *m_effDistanceEdit;

    Wt::WCheckBox *m_writeEnergyCalCb;

    Wt::WText *m_warningsTxt;
  };//class GenieCnfOptionsWidget

}//namespace ExportSpecFileCAM

#endif //ExportSpecFileCAM_h
