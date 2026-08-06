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

#include <Wt/WSignal>
#include <Wt/WModelIndex>
#include <Wt/WContainerWidget>
#include <Wt/WAbstractItemModel>

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
  class WCheckBox;
  class WComboBox;
  class WPushButton;
  class WContainerWidget;
}//namespace Wt

namespace SandiaDecay
{
  struct Nuclide;
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

    /** Yield, in photons/xrays per decay (i.e., a value from 0 to 1, not a percent). */
    float yield = 0.0f;
    /** Negative if not known (Genie will estimate an uncertainty from the value itself). */
    float yield_uncert = -1.0f;

    bool is_xray = false;

    /** Set once, after all sources/lines have been built and merged; matches (as closely as
     practical, for display purposes) what `CAMInputOutput::CAMIO::AssignKeyLines()` will
     independently compute when the file is actually written.
     */
    bool is_key_line = false;

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
   @param combine_unresolvable_lines If true, lines within about 0.5 keV of each other (within
          a source) are merged into a single yield-weighted line, mirroring how the GENIE
          Library Editor (and the "Peak-Map" tool) combine unresolvable lines.
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
                                  const std::map<const SandiaDecay::Nuclide *, double> &nuclide_ages = {},
                                  std::vector<std::string> *warnings = nullptr );

  /** Flattens the checked (`included`) sources/lines of `sources` into the line-list that
   `SpecUtils::SpecFile::write_cnf(...)`'s `CnfGenieExtras` expects.
   */
  std::vector<CAMInputOutput::CnfGenieExtras::LibraryLine>
      to_library_lines( const std::vector<GenieLibrarySource> &sources );


  /** Fits `FWHM = a + b*sqrt(energy)` (Genie's FWHM equation form) to the given DRF's FWHM
   function, sampled at points across its valid energy range (or 59-2614 keV, if the DRF
   doesn't specify a valid range) using an unweighted linear least-squares fit.

   @returns {a, b}
   */
  std::pair<float,float> fit_genie_fwhm_from_drf( const DetectorPeakResponse &drf );

  /** The default `{FWHMOFF, FWHMSLOPE}` Genie itself uses for HPGe or NaI detectors (see the
   Genie 2000 Customization Tools Manual's detector-type-dependent initial parameters table).
   */
  std::pair<float,float> default_genie_fwhm( const bool is_hpge );


  /** The result of converting a `DetectorPeakResponse`'s intrinsic efficiency function into
   the closest available Genie efficiency model.

   Note: `CAMIO`'s GEOM-block write support (see `AddEfficiencyPoint(s)`) only writes tabulated
   (energy, efficiency) points - not persisted polynomial coefficients (the Genie 2000
   Customization Tools Manual describes coefficient parameters like `CAM_F_DHCALFAC1`, but their
   on-disk layout isn't known to this code, and `CAMIO`'s read-side implementation of
   `GetEfficiencyPoints()` only ever extracts points too) - so regardless of `model`, Genie is
   expected to (re)derive its internal fit from `points` using the algorithm named by `model`.
   */
  struct GenieEfficiencyResult
  {
    CAMInputOutput::CAMIO::EfficiencyModel model = CAMInputOutput::CAMIO::EfficiencyModel::Unknown;
    std::vector<CAMInputOutput::EfficiencyPoint> points;
  };//struct GenieEfficiencyResult

  /** Converts a DetectorPeakResponse's intrinsic efficiency into the closest available Genie
   efficiency model, writing out sampled (energy, efficiency) points in all cases (see note on
   `GenieEfficiencyResult`), tagged with the best-matching Genie model name:
   - `kEnergyEfficiencyPairs` -> `SPLINE`, using the DRF's own tabulated points directly.
   - `kExpOfLogPowerSeries` -> `DUAL`, Genie's own default model form for Ge/NaI detectors, and
     algebraically the same functional family (points are sampled from the exact formula).
   - `kFunctialEfficienyForm` -> points are sampled across the DRF's valid range, and tagged
     `DUAL` if a ln(eff)-vs-ln(energy) polynomial fits those points well, else `SPLINE`.

   @param warnings If non-null, human-readable warnings are appended here.
   */
  GenieEfficiencyResult convert_efficiency_to_genie( const DetectorPeakResponse &drf,
                                                     std::vector<std::string> *warnings = nullptr );


  /** Two-level (source -> line) table model backing the GENIE library preview/selection table
   in the CNF export options: top-level rows are `GenieLibrarySource`s, and their children are
   that source's `GenieLibraryLine`s.  Each row (source or line) has a checkbox in the
   `Column::Include` column; checking/unchecking a source checks/unchecks all its lines.
   */
  class GenieLibraryModel : public Wt::WAbstractItemModel
  {
  public:
    enum class Column : int
    {
      /** Nuclide name (source rows) or energy in keV (line rows). */
      Name,
      /** Half-life (source rows) or yield/key-line/x-ray info (line rows). */
      Info,
      /** Editable age (source rows of ageable nuclides only; blank/unused otherwise). */
      Age,
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

    virtual Wt::WModelIndex index( int row, int column, const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual Wt::WModelIndex parent( const Wt::WModelIndex &index ) const override;
    virtual int rowCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual int columnCount( const Wt::WModelIndex &parent = Wt::WModelIndex() ) const override;
    virtual boost::any data( const Wt::WModelIndex &index, int role = Wt::DisplayRole ) const override;
    virtual bool setData( const Wt::WModelIndex &index, const boost::any &value, int role = Wt::EditRole ) override;
    virtual Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const override;
    virtual boost::any headerData( int section, Wt::Orientation orientation = Wt::Horizontal, int role = Wt::DisplayRole ) const override;

    /** Emitted after a source's age is edited by the user (with the new age, in seconds); the
     owning widget should re-run `build_genie_library(...)` with an updated `nuclide_ages` map
     and call `setSources(...)` again to refresh the line yields.
     */
    Wt::Signal<const SandiaDecay::Nuclide *, double> &ageEdited();

  private:
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
    void updateForFile( const std::shared_ptr<const SpecMeas> &spec, const std::set<int> &samples );

    bool writeLibrary() const;
    bool writeFwhm() const;
    bool writeEfficiency() const;
    bool writeEnergyCal() const;

    /** Builds the `CnfGenieExtras` reflecting the current widget state (checked library
     sources/lines, chosen FWHM source, whether to include efficiency/energy-cal), for the
     spectrum/samples last passed to `updateForFile(...)`.
     */
    CAMInputOutput::CnfGenieExtras currentExtras() const;

  protected:
    void rebuildLibraryTable();
    void handleLibraryOptionsChanged();
    void handleSourceAgeEdited( const SandiaDecay::Nuclide *nuc, double age_seconds );
    void handleFwhmSourceChanged();
    void handleFitFwhmFromPeaksClicked();
    void handleFwhmFitFromToolUpdated( std::shared_ptr<DetectorPeakResponse> new_drf );

    InterSpec *m_interspec;
    std::weak_ptr<const SpecMeas> m_spec;
    std::set<int> m_samples;
    std::map<const SandiaDecay::Nuclide *, double> m_nuclide_ages;

    Wt::WCheckBox *m_writeLibraryCb;
    Wt::WComboBox *m_libraryModeCb;
    NativeFloatSpinBox *m_thresholdEdit;
    Wt::WCheckBox *m_combineLinesCb;
    RowStretchTreeView *m_libraryTable;
    GenieLibraryModel *m_libraryModel;

    Wt::WCheckBox *m_writeFwhmCb;
    Wt::WComboBox *m_fwhmSourceCb;
    Wt::WPushButton *m_fitFwhmFromPeaksBtn;
    bool m_have_manual_fwhm;
    std::pair<float,float> m_manual_fwhm_coeffs;

    Wt::WCheckBox *m_writeEfficiencyCb;
    Wt::WCheckBox *m_writeEnergyCalCb;

    Wt::WText *m_warningsTxt;
  };//class GenieCnfOptionsWidget

}//namespace ExportSpecFileCAM

#endif //ExportSpecFileCAM_h
