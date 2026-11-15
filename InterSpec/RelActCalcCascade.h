#ifndef RelActCalcCascade_h
#define RelActCalcCascade_h
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
#include <utility>
#include <optional>

struct Material;
class DetectorPeakResponse;
namespace SandiaDecay{ struct Nuclide; }
namespace ceelo{ class SummingOutProgram; }
namespace RelActCalc{ template<typename T> struct PhysModelShield; }
namespace GammaInteractionCalc{ class ShieldScatterAugment; }

/** True-coincidence (cascade) summing for a source described the way the "Physical" relative-efficiency model
 describes it: a self-attenuating slab of one material, behind external attenuators, in front of a detector
 with a known response - plus one geometric quantity the relative-efficiency model does not have, the
 fractional solid angle the detector subtends from the front of the source (i.e. the distance).

 The summing efficiencies are physical, not the fitted relative-efficiency curve.  For a photon of energy E
 emitted at a node `n` of the source (a depth into the slab), the full-energy-peak or total efficiency is

     eps_fep(E; n) = g * eps_int,fep(E) * T_n(E)
     eps_tot(E; n) = g * eps_int,tot(E) * ( T_n(E) + A_n(E) )
     T_n(E) = exp( -mu_self(E) * depth_n ) * prod_k exp( -mu_k(E) * AD_k )

 with `g` the fractional solid angle from the front of the source, `eps_int` the detector's intrinsic
 full-energy or total efficiency, and `A_n` the photons that scatter in the material in front of the node and
 still reach the detector (`GammaInteractionCalc::ShieldScatterAugment`, as the Activity/Shielding fit's
 summing correction uses): they count for the total (any energy deposited) efficiency, so attenuation
 lowers summing-out by less than it lowers the peaks.  Both photons of a coincident pair leave the same
 node, so a partner's efficiency is conditional on the primary having escaped from there; the factors are
 full-energy-efficiency weighted averages over the nodes.

 Not modeled: sums of three or more photons (at most a few percent of a two-photon peak for uranium, but
 up to a fifth of one for a cascade-rich nuclide like Eu-152 touching the detector); summing-out of a line
 that is itself an x-ray (only gammas get a factor, though x-rays are counted as partners); sum peaks of a
 decay whose level scheme CeeLo does not accept (its lines still get summing-out factors, from pairwise
 coincidence data; see `SummingModel::Source::line_complete`); and photons below the energy the detector
 response is characterized from (see `DetectorEfficiencies::lower_energy`), which are taken as not detected.

 The attenuation enters smoothly except through `A_n`, which is interpolated linearly in a table, so has
 kinks (it is continuous) at the table's atomic numbers and areal densities.

 Summing-out is the exact level-scheme expectation CeeLo's analytic engine computes
 (`ceelo::SummingOutProgram`); summing-in and coincidence-sum peaks come from the engine's photon-pair table
 (`ceelo::coincidence_sum_pairs`), with the pairs' coincident survival the engine's exact one for pairs of
 level-scheme transitions, and estimated for the rest (see `SummingModel::evaluate`).

 Everything is a continuous function of the inputs, and is templated on their type, so it can be evaluated
 inside a Ceres cost function with `ceres::Jet`s.  Units are PhysicalUnits, except energies, which are keV.
 */
namespace RelActCalcCascade
{
  /** Where the detector's total (any energy deposited) efficiency came from. */
  enum class TotalEffSource
  {
    /** The detector response function carries one. */
    DetectorResponse,

    /** The response has none (or an unphysical one), so a reference total-to-peak ratio was applied to the
     response's own full-energy efficiency: a blend of a small (50 x 30 mm) and a large (InterSpec's example
     response, data/cascade_summing_example_response.xml) HPGe crystal's, chosen by how fast the response's
     full-energy efficiency falls with energy - which is what the crystal's size shows up as.
     */
    ExampleDetectorRatio
  };//enum class TotalEffSource


  /** The detector's far-field intrinsic efficiencies at a list of energies, and its face geometry. */
  struct DetectorEfficiencies
  {
    /** Intrinsic full-energy-peak and total efficiency at each energy; zero below #lower_energy.  Above the
     response's stated range its own (extrapolated) values are used, limited to [0,1].
     */
    std::vector<double> fep, total;

    /** The energy (keV) below which photons are taken as not detected: the lower end of the range the
     response says it is characterized over, or 40 keV if it does not say.  A response's extrapolation below
     that is commonly an exponential of a polynomial in log-energy, which can be anything.
     */
    double lower_energy = 0.0;

    TotalEffSource total_source = TotalEffSource::DetectorResponse;

    /** Why `total_source` is not the detector response; empty if it is. */
    std::string total_source_reason;

    /** The diameter and setback the response's far-field intrinsic efficiencies are quoted against. */
    double detector_diameter = 0.0;
    double detector_setback = 0.0;

    /** The crystal disk the response models (see `DetectorPeakResponse::modeledCrystalDiameter`); the same
     as the above unless the response is geometry-aware.
     */
    double crystal_diameter = 0.0;
    double crystal_setback = 0.0;
  };//struct DetectorEfficiencies


  /** Tabulates the detector's intrinsic efficiencies at `energies_kev`.

   The total efficiency is the response's own when it has one that is physical (not above one, not below the
   full-energy efficiency) over the response's energy range; otherwise a reference total-to-peak ratio is
   used (see `TotalEffSource::ExampleDetectorRatio`).  Which, is a property of the response - it does not
   depend on `energies_kev` - so every use of one response agrees.

   Throws if the response is null, invalid, fixed-geometry, or has no diameter; or if the fall-back is needed
   and the example response can not be loaded (from InterSpec's static data directory).
   */
  DetectorEfficiencies detector_efficiencies( const std::shared_ptr<const DetectorPeakResponse> &drf,
                                              const std::vector<double> &energies_kev );


  /** Fractional solid angle of a disk of radius `detector_radius`, from a point on its axis a distance
   `distance` away.  Templated so the distance may carry derivatives.
   */
  template<typename T>
  T on_axis_solid_angle( const T &distance, const double detector_radius );

  /** Inverse of #on_axis_solid_angle. `solid_angle` must be in (0, 0.5]. */
  template<typename T>
  T on_axis_distance( const T &solid_angle, const double detector_radius );


  /** The depth quadrature of the self-attenuating slab: fractions of its areal density, and weights that sum
   to one.  Fixed (it does not depend on the areal density), so results are smooth in the fit parameters;
   graded toward the front face, where low-energy lines come from.
   */
  struct DepthNodes
  {
    std::vector<double> fraction, weight;
  };//struct DepthNodes

  const DepthNodes &depth_nodes();


  /** An energy range to predict coincidence-sum peaks in (e.g. a fit ROI). */
  struct SumRange
  {
    double lower_energy = 0.0, upper_energy = 0.0;

    /** Sum peaks whose rate, per decay of the source, could not reach this even with the source touching the
     detector and nothing attenuated, are left out.
     */
    double min_rate_per_decay = 0.0;
  };//struct SumRange


  /** One source nuclide: the gamma lines to get factors for, and where to look for its sum peaks. */
  struct SourceInput
  {
    const SandiaDecay::Nuclide *nuclide = nullptr;

    /** Age of the nuclide (its progeny are grown in from a pure parent, as RelActCalcAuto ages its lines). */
    double age = 0.0;

    /** Energies (keV) of the gamma lines of this nuclide the caller models. */
    std::vector<double> line_energies;

    std::vector<SumRange> sum_ranges;
  };//struct SourceInput


  /** A photon pair of one decay that adds into a sum peak. */
  struct SumPair
  {
    /** Indices into `SummingModel::energies()` of the two photons (a is the lower energy one). */
    size_t energy_a = 0, energy_b = 0;

    /** Pairs emitted per (current) decay of the source nuclide; includes the angular-correlation factor. */
    double weight = 0.0;

    /** P(b emitted | a emitted), and the reverse - without the angular factor. */
    double prob_b_given_a = 0.0, prob_a_given_b = 0.0;

    /** The pair's index among the summing-out program's pairs (see `SummingModel::Source::program`), which
     gives its exact coincident survival; -1 for a pair with a member that is not on a level-scheme transition
     (an electron-capture x-ray, annihilation photon, or unplaced gamma), whose survival is estimated.
     */
    int program_pair = -1;

    /** For a member that is a gamma: its line in the source's summing-out program, used to estimate the
     coincident survival when `program_pair` is -1; else -1.
     */
    int program_line_a = -1, program_line_b = -1;
  };//struct SumPair


  /** Coincident photon pairs of one source nuclide whose energies add to (nearly) the same value. */
  struct SumPeak
  {
    /** Mean sum energy (keV) of the pairs, weighted by their rate at the detector's intrinsic efficiencies. */
    double energy = 0.0;

    /** If the sum lands on one of the source's modeled lines (summing-in to it), that line's index in
     `SourceInput::line_energies`; else -1, and this is a peak of its own.
     */
    int line_index = -1;

    /** Upper bound of the rate, per decay of the source nuclide: source touching the detector, no attenuation. */
    double max_rate_per_decay = 0.0;

    std::vector<SumPair> pairs;
  };//struct SumPeak


  /** The geometry and attenuation the summing is evaluated for. */
  template<typename T>
  struct EvalInput
  {
    /** Fractional solid angle (of 4 pi) the detector subtends from the front of the source, on its axis. */
    T front_solid_angle = T(0.0);

    /** Areal density of the self-attenuating slab; zero (and `self_mu` empty) if there is none. */
    T self_areal_density = T(0.0);

    /** Mass-attenuation coefficient of the self-attenuating material at each of `SummingModel::energies()`;
     empty if there is no self-attenuation.
     */
    std::vector<T> self_mu;

    /** Transmission through all external attenuators at each of `SummingModel::energies()`. */
    std::vector<T> external_transmission;

    /** For the photons that scatter in the attenuators and still reach the detector: the (mass-weighted)
     atomic number of the self-attenuating material, and the external attenuators' summed areal density and
     summed (atomic number times areal density).
     */
    T self_atomic_number = T(0.0);
    T external_areal_density = T(0.0);
    T external_an_areal_density = T(0.0);
  };//struct EvalInput


  /** What the full-energy-peak of each member of a sum pair contributes, at this geometry. */
  template<typename T>
  struct SumPairValue
  {
    /** Absolute full-energy efficiency of photon b, averaged over where photon a's detected counts come from,
     times the probability nothing else of the decay deposits (the pair's coincident survival); and the same
     with the roles swapped.  A sum peak's counts are then

         (counts of line a per unit yield) * weight * b_given_a  ==  (same for b) * weight * a_given_b

     in a perfectly consistent model; a caller with a fitted relative-efficiency curve picks or combines them.
     */
    T b_given_a = T(0.0), a_given_b = T(0.0);
  };//struct SumPairValue


  template<typename T>
  struct EvalResult
  {
    /** `[source][line]`: factor the full-energy counts of the line are multiplied by for summing-out. */
    std::vector<std::vector<T>> summing_out;

    /** `[source][sum peak][pair]`. */
    std::vector<std::vector<std::vector<SumPairValue<T>>>> sum_pairs;
  };//struct EvalResult


  /** The parameter-independent part of the summing calculation for a set of source nuclides and a detector:
   built once, then evaluated for many geometries.  Immutable once built, so it may be shared between threads.
   */
  class SummingModel
  {
  public:
    /** A sum of two photons within this (keV) of one of the source's own lines counts as adding to that line
     (`SumPeak::line_index`), rather than as a peak of its own.
     */
    static constexpr double sm_sum_on_line_kev = 0.25;

    struct Source
    {
      SourceInput input;

      /** Whether each line was found among the nuclide's cascade gammas (its factor is 1 if not), and whether
       its decay branches all have an accepted level scheme (if not, part of its factor uses CeeLo's simpler
       pairwise form).
       */
      std::vector<char> line_found, line_complete;

      /** Emission probability of each line per (current) decay of the nuclide, from the cascade model. */
      std::vector<double> line_emission;

      std::vector<SumPeak> sum_peaks;

      /** Summing-out program; its first `input.line_energies.size()` lines are those lines, the rest are
       gamma members of sum pairs.
       */
      std::shared_ptr<const ceelo::SummingOutProgram> program;

      /** For each of the program's energies, its index in `SummingModel::energies()`, or -1 if a photon of
       that energy can not matter (its efficiency is taken as zero).
       */
      std::vector<int> program_energy;

      /** Index in `SummingModel::energies()` of each of the program's lines. */
      std::vector<size_t> program_line_energy;
    };//struct Source

    /** @param sources The source nuclides (of one relative-efficiency curve).
     @param drf The detector response.
     @param shield_scatter Whether to count the photons that scatter in the attenuators and still reach the
            detector toward the total efficiency (see #models_shield_scatter); without, only unscattered
            photons count, which has closed forms to test against.

     Throws on an unusable detector response (see #detector_efficiencies), or a null nuclide.
     */
    SummingModel( const std::vector<SourceInput> &sources,
                  const std::shared_ptr<const DetectorPeakResponse> &drf,
                  const bool shield_scatter = true );

    const std::vector<Source> &sources() const { return m_sources; }

    /** Every energy (keV, sorted) the evaluation needs attenuation values at. */
    const std::vector<double> &energies() const { return m_energies; }

    const DetectorEfficiencies &detector() const { return m_detector; }

    /** The largest front solid angle a source can have: a point touching the detector's face. */
    double max_front_solid_angle() const { return m_max_solid_angle; }

    /** Whether the geometry changes anything: false if no line has a coincident partner the detector can
     see, and there is no sum peak - for example, for sources with a single gamma.
     */
    bool has_effect() const { return m_has_effect; }

    /** The front solid angle of a source whose near surface is `distance` from the detector's face: the
     fraction that, times the detector's far-field intrinsic efficiencies, is its absolute efficiency there.

     For a geometry-aware response this is read off the response's own (near-field) total efficiency, at a
     reference energy; otherwise it is the solid angle of the detector's disk.  A source face that is a disk of
     `source_radius` (zero: a point on the detector's axis) sees less of the detector toward its edge; what a
     summing correction then needs is not the mean solid angle `<g>` over the face but `<g^2>/<g>` - a
     coincident partner comes from where the detected photon did - and that is what this returns.
     */
    double front_solid_angle( const double distance, const double source_radius = 0.0 ) const;

    /** Inverse of #front_solid_angle: the distance from the detector's face to the front of the source, between
     zero and 10 m.
     */
    double front_distance( const double solid_angle, const double source_radius = 0.0 ) const;

    /** The absolute full-energy and total efficiencies, at each of #energies, for an unattenuated point on the
     detector's axis whose front solid angle is `solid_angle` - what #evaluate uses.

     For a geometry-aware response these are its own near-field efficiencies (see `m_near_solid_angle`); for
     any other, the solid angle times the far-field intrinsic efficiencies.
     */
    template<typename T>
    void point_efficiencies( const T &solid_angle, std::vector<T> &fep, std::vector<T> &total ) const;

    /** Whether the photons that scatter in the attenuators and still reach the detector are modeled; they
     are not if GADRAS's shield-scatter data (sandia.shieldscatter.db, in InterSpec's static data directory)
     could not be loaded, and attenuation then lowers the summing corrections by too much.
     */
    bool models_shield_scatter() const;

    /** Mass-attenuation coefficients of `material` at #energies; computed on first use for a material
     (thread-safe).
     */
    std::shared_ptr<const std::vector<double>> material_mu( const std::shared_ptr<const Material> &material ) const;

    /** Fills the attenuation members of `input` from the Physical model's shields. */
    template<typename T>
    void set_attenuation( const std::optional<RelActCalc::PhysModelShield<T>> &self_atten,
                          const std::vector<RelActCalc::PhysModelShield<T>> &external_attens,
                          EvalInput<T> &input ) const;

    /** Evaluates the summing-out factors and the sum-pair values.

     The coincident survival of a sum pair (nothing else emitted with the pair deposits energy) is the
     engine's exact one for pairs of level-scheme transitions - nearly all of them.  For the rest it is
     estimated from the summing-out factors of its gamma members: for member a with partner b,
     `c_out(a) / (1 - P(b|a) eps_tot(b))` - the summing-out of a, less the part b itself is responsible for.
     Sums of three photons are not modeled.
     */
    template<typename T>
    void evaluate( const EvalInput<T> &input, EvalResult<T> &result ) const;

  private:
    std::vector<Source> m_sources;
    std::vector<double> m_energies;
    DetectorEfficiencies m_detector;
    double m_max_solid_angle = 0.0;
    bool m_has_effect = false;

    std::shared_ptr<const DetectorPeakResponse> m_drf;

    /** For a geometry-aware response: the energy its near-field total efficiency is read at, and its far-field
     intrinsic total efficiency there (zero if distances are not taken from the response).
     */
    double m_ref_energy = 0.0, m_ref_total_intrinsic = 0.0;

    /** For a geometry-aware response: how its near-field efficiency, at each of #energies, departs from
     "front solid angle times far-field intrinsic efficiency".  Low-energy photons stop near the crystal's
     face and high-energy ones deeper, where the solid angle is smaller, so close to the detector the
     efficiency does not scale with one solid angle at every energy; the departure grows with nearness.

     `m_near_solid_angle` is a grid of front solid angles ascending from zero (far away, where the ratios
     are one); the ratio tables are [energy][grid point], with the slopes of their monotone cubic
     interpolants.  Empty for responses that are not geometry-aware.

     The ratios are interpolated in energy from the response's `NearFieldGrid`, which is made once per
     response: each near-field efficiency is a ray-trace, and a fit makes several models.
     */
    std::vector<double> m_near_solid_angle;
    std::vector<std::vector<double>> m_near_fep_ratio, m_near_fep_slope, m_near_tot_ratio, m_near_tot_slope;

    struct NearFieldGrid;
    std::shared_ptr<const NearFieldGrid> near_field_grid() const;
    void build_near_field_tables();
    void build_scatter_table();

    /** The near-field ratio at `solid_angle`, from one energy's row of the tables. */
    template<typename T>
    T near_field_ratio( const std::vector<double> &ratios, const std::vector<double> &slopes,
                        const T &solid_angle ) const;

    /** The on-axis front solid angle of a point source at `distance`. */
    double axis_solid_angle( const double distance ) const;

    struct MaterialMuCache;
    std::shared_ptr<MaterialMuCache> m_material_mu;

    /** The scattered-photon part of the total efficiency behind attenuators, at #energies; null if the
     scatter data could not be loaded.
     */
    std::shared_ptr<const GammaInteractionCalc::ShieldScatterAugment> m_scatter;
  };//class SummingModel
}//namespace RelActCalcCascade

#endif //RelActCalcCascade_h
