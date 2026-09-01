#ifndef GammaInteractionCalc_h
#define GammaInteractionCalc_h
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
#include <array>
#include <deque>
#include <tuple>
#include <memory>
#include <atomic>
#include <vector>
#include <utility>

#include <boost/asio/deadline_timer.hpp>

#include <Wt/WColor.h>

#include "Minuit2/FCNBase.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

class PeakDef;
struct Material;
class DetectorPeakResponse;

namespace rapidxml
{
  template<typename Ch> class xml_node;
}

namespace SpecUtils
{
  class Measurement;
}

namespace GammaInteractionCalc
{

class CascadeSummingCalc;

/** Maximum areal density allowed for computations, in units of g/cm2 -
 
 Not in units of PhysicalUnits; e.g., you need to multiple by PhysicalUnits::cm2 / PhysicalUnits::g before using for computation..
 
 Currently value of 500 is fairly arbitrary.
 */
static const double sm_max_areal_density_g_cm2 = 500.0;

/** We will let trace source activities be specified in a few different ways, since it makes sense from a users-perspective, and because
 we may want total activity to vary independently from physical volume and then this just makes things a little easier to track in the Minuit
 variables.
 */
enum class TraceActivityType : int
{
  TotalActivity,
  ActivityPerCm3,
  ExponentialDistribution,  //Activity of Bq/m2 of the entire column of soil
  ActivityPerGram, //Needs to come last as wont be available if a "void" material
  //ActivityPPM,
  NumTraceActivityType
};//enum class TraceActivityType

/** Gives string representation to a TraceActivityType value. */
const char *to_str( const TraceActivityType type );


/** Enum that lists the geometry types we can compute volumetric sources for. */
enum class GeometryType : int
{
  Spherical,
  CylinderEndOn,
  CylinderSideOn,
  Rectangular,
  NumGeometryType
};//enum class GeometryType

/** Gives string representation to a GeometryType value. */
const char *to_str( const GeometryType type );

  

//Returned in units of 1.0/[Length], so that
//  exp( -transmition_length_coefficient(...) * thickness)
//  gives you the probability a gamma of given energy will go through the
//  material of given thickness.
//Energy units should be in SandiaDecay/PhysicalUnits units.
double transmition_length_coefficient( const Material *material, float energy );


/** The FULL-ENERGY-PEAK survival-removal coefficient (plan 3.4):

      mu_rem = transmition_length_coefficient(material, energy) - f_win * mu_Compton

 Rayleigh scattering is elastic, so it cannot move a photon out of the peak window, and forward
 Compton scatters whose energy loss stays inside the window are still counted in the peak.  Using
 plain mu_total therefore OVER-attenuates the peak; CeeLo measures that at -8..-16% at 60 keV for
 low-Z contents at depth.

 IMPORTANT: InterSpec's mu already excludes Rayleigh - `massAttenuationCoefficientElement` returns
 compton+photoelectric+pair and the SNL path returns 0.0 for RayleighScatter - so CeeLo's
 `mu_total - mu_Rayleigh - f_win*mu_Compton` has no Rayleigh term to subtract here.  Subtracting one
 would double-count, and would do so silently because that call currently returns zero.

 `window_keV` is a HALF-width (CeeLo's +-win convention), i.e. FWHM/2 of the detector the response
 is anchored on.  Pass <= 0 to get no credit (returns mu_total unchanged).

 The in-window fraction is a quadrature over the Klein-Nishina distribution weighted by the
 material's incoherent scattering function, so this is NOT cheap: it is memoised per
 (material, energy, window) and must never be called inside a per-ray or per-element loop.

 Headroom, so nobody expects more of this than it can give: the whole correction is bounded by
 f_win * (mu_Compton/mu_total).  Measured at a 0.35 keV half-window and 60 keV that is 0.30% of mu
 for steel (photoelectric dominates, so there is almost nothing to credit) versus 3.1% for water;
 above 344 keV it is under 0.3% for both.
 */
double fep_survival_removal_coefficient( const Material *material, float energy,
                                         double window_keV );
double transmition_coefficient_material( const Material *material, float energy,
                                float length );


/** A convenience call to #transmition_coefficient_material that uses a static (compile-time defined) definition of air.
 
 Example use of this function:
 \code {.cpp}
 double energy = 661.0*PhysicalUnits::keV;
 double distance = 3.2*PhysicalUnits::cm;
 double mu = transmission_coefficient_air( energy, distance );
 double transmission_fraction = exp( -mu );
 \endcode
 */
double transmission_coefficient_air( float energy, float length );

/** Similar to #transmission_coefficient_air, but not including length.
 
 Example use of this function:
 \code {.cpp}
 double energy = 661.0*PhysicalUnits::keV;
 double distance = 3.2*PhysicalUnits::cm;
 double mu = transmission_length_coefficient_air( energy );
 double transmission_fraction = exp( - mu * distance );
 \endcode
 */
double transmission_length_coefficient_air( float energy );


//Returned in units of [Length]^2/[mass], so that
//  exp( -mass_attenuation_coef * areal_density )
//  gives you the probability a gamma of given energy will go through the
//  material with given atomic number and areal_density.
//  The quantity retuned by this function is commonly labeled μ
//Energy units should be in SandiaDecay/PhysicalUnits units.
double mass_attenuation_coef( float atomic_number, float energy );
double transmition_coefficient_generic( float atomic_number, float areal_density,
                                float energy );

/** Returns the average gammas per second during a measurement, when the decay of the input nuclides is accounted for.
 @param mixture The nuclide mixture - currently may only have a single parent nuclide.
 @param intial_age The mixture nuclides age when the measurement is started.
 @param measurement_duration The duration of the measurement.
 
 @returns The average rates for each of the energies, ordered by energy.
 */
std::vector<SandiaDecay::EnergyRatePair> decay_during_meas_corrected_gammas(
                                                    const SandiaDecay::NuclideMixture &mixture,
                                                    const double initial_age,
                                                    const double measurement_duration );
    
void example_integration();

// When debugging we will grab a static mutex so we dont get jumbled stdout
#define DEBUG_RAYTRACE_CALCS 0



//exit_point_of_sphere_z(...): Makes exit_point[3]=<x1,y1,z1> be the point of
//  intersection of a sphere or radius 'sphere_rad' centered at the origin, with
//  the line pointing from source_point[3]=<x0,y0,z0> towards
//  <0,0,observation_dist>.  Returns the distance between source_point[3]
//  and exit_point[3].
//  If postiveSolution==false, then the solution at a negative z will be given
//Note: CAN be used with a 2-dimensional integral which holds phi constant
//Note: it is safe to pass the same array as source_point[3] and exit_point[3]
double exit_point_of_sphere_z( const double source_point[3],
                               double exit_point[3],
                               double sphere_rad,
                               double observation_dist,
                               bool postiveSolution = true );

/** An enum to to tell #cylinder_line_intersection which exit point from sphere you want.

 A better name would be CylinderIntersectionDirection, but thats too long.

 Defined in terms of the ray parameter: writing the ray as P(t) = source + t*(detector - source), t
 increases monotonically toward the detector, so #TowardDetector is the larger-t crossing and
 #AwayFromDetector the smaller-t one.  For a source inside the volume this means #AwayFromDetector is
 behind the source (t < 0); for a source outside it is the point the ray enters the volume.
 */
enum class CylExitDir
{
  TowardDetector,
  AwayFromDetector
};

/** Starting from a 'source_point' within the volume of the cylinder, and heading towards the 'detector_point' (think center
 of the detector face), returns the total attenuation coefficient along the path, including recursing into any sub-tubes, as well as
 sets the point where the ray leaves the cylinder.
 
 Note: the cylinder is always oriented along the z-axis, and centered at {0,0,0}.
 
 @param[in] radius The outer radius of the cylinder
 @param[in] half_length The half-length of the cylinder
 @param[in] source The {x, y, z} source location; must be within cylinder volume, and transformed so that {0,0,0} is center of cylinder.
 @param[in] detector The {x, y, z} point on the detector face we care about (so center of detector, unless you are integrating over
 the detector face), in the coordinate system where cylinder is centered at {0,0,0}.
 @param[out] exit_point The final exit point from the cylinder, where the path will no longer go through the volume.
 @returns The distance from source location to exit point.  Returns 0.0 if line does not intersect cylinder.  Note that this is not the
 distance in the cylinder, but the total distance from source to exit point, so if source is outside volume, may be larger than cylinder
 dimensions.  Note that if source is on radius of cylinder, the value returned will be 0.0., for all values of CylExitDir.
 \code{.cpp}
 double exit_point[3];
 const double distance_in_m = cylinder_line_intersection( 0.5*m, 100*cm, {0,0.1*m,20*cm}, {0,0,10*m}, exit_point );
 const double trans_fraction = exp( -trans_coef ); //trans_fraction will be between 0 and 1.
 \endcode
 
 TODO: this function should be broken into two separate functions.  One to handle finding the exit
 point when you know the source is inside the volume.  And one to find both intersection points (if any) of external points.  This would
 both increase the efficiency of the function, and also make the use cleaner/easier.
 */
double cylinder_line_intersection( const double radius, const double half_length,
                              const double source[3],
                              const double detector[3],
                              const CylExitDir direction,
                              double exit_point[3] ) noexcept;

/** Provides the distance and exit point of a ray, originating inside a rectangle, when it goes from \p source to \p detector.
 
 @param[in] half_width The half-width of the rectangle; e.g., the x-extent of the rectangle.
 @param[in] half_height The half-height of the rectangle; e.g., the y-extent of the rectangle.
 @param[in] half_depth The half-depth of the rectangle; e.g., the z-extent of the rectangle.
 @param[in] source The {x,y,z} location the ray originates from (e.g., center of source voxel when integrating over).  This point must
            be inside, or on surface of rectangle, or else results are not defined.
 @param[in] detector The {x,y,z} location the ray terminates (e.g., the center face of the detector). This point must be outside, or on
            surface of rectangle - and not equal in values to \p source, or else results are not defined.
 @param[out] exit_point The {x,y,z} location where the ray exits the rectangle.  Note that this can be same array as \p source.
 @returns The distance inside the rectangle the ray traverses.
 */
double rectangle_exit_location( const double half_width, const double half_height,
                               const double half_depth,
                               const double source[3],
                               const double detector[3],
                               double exit_point[3] ) noexcept;

/** Provides the intersection points on the rectangle when the ray originates from \p source and terminates at \p detector.
 
 @param[in] half_width The half-width of the rectangle; e.g., the x-extent of the rectangle.
 @param[in] half_height The half-height of the rectangle; e.g., the y-extent of the rectangle.
 @param[in] half_depth The half-depth of the rectangle; e.g., the z-extent of the rectangle.
 @param[in] source The {x,y,z} location the ray originates from (e.g., center of source voxel when integrating over).  This point must
 be outside, or on surface of rectangle, or else results are not defined.
 @param[in] detector The {x,y,z} location the ray terminates (e.g., the center face of the detector)
 @param[out] near_source_intersection The {x,y,z} location where the ray enters the rectangle.  If ray does not intersect the
             rectangle, the values in the array are not specified.
 @param[out] near_detector_intersection The {x,y,z} location where the ray exits the rectangle.  If ray does not intersect the
 rectangle, the values in the array are not specified.
 @returns If the ray intersects the rectangle.
 */
bool rectangle_intersections( const double half_width, const double half_height,
                               const double half_depth,
                               const double source[3],
                               const double detector[3],
                               double near_source_intersection[3],
                               double near_detector_intersection[3] ) noexcept;


//distance(...): returns distance between two points specified in terms of x, y,
//  and z.
double distance( const double point_a[3], const double point_b[3] );


/*
 *The following functions are depreciated, but left commented out since I will
 *likely need to refer to them later - wcjohns 20121127
//exit_point_of_sphere_x(...): similar to exit_point_of_sphere_z(...), but
//  for the detector at <observation_dist,0,0>
//Makes x0,y0,z0 be the point of intersection of a sphere or radius
//  'sphere_rad' centered at the origin, with the line pointing from <x0,y0,z0>
//  towards <observation_dist,0,0>
//  Note: can not be used with a 2-dimensional integral which holds phi or theta
//        constant
void exit_point_of_sphere_x( double &x0, double &y0, double &z0,
                            double sphere_rad, double observation_dist );

double len_in_sphere_x( double r, double theta,
                        double phi, double sphere_rad, double observation_dist );
*/

/** The volumetric (trace / self-attenuating) source calculator, over `double` for the display and
 `ceres::Jet<>` for the auto-differentiated fit - defined in GammaInteractionCalc_imp.hpp. */
template<typename T>
struct DistributedSrcCalcT;

/** Which kind of shell a #DistributedSrcCalcT layer is: an attenuating material of finite extent,
 or a zero-extent "generic" shielding given directly as an areal density.

 Namespace-scope (it used to live inside the now-deleted double-only `DistributedSrcCalc`) so the
 templated calculator, which is the single volumetric model, owns no vestigial dependency.
 */
enum class ShellType
{
  Material, Generic
};//enum class ShellType
  

/** A struct to hold enough information for each point on the pulls chart that shows comparison of
 detected vs model, for each peak
 */
struct PeakResultPlotInfo
{
  /** Energy of the gamma associated with the peak. */
  double energy;
  
  /** `(observed_counts - expected_counts) / observed_uncertainty` */
  double numSigmaOff;
  
  /** `observed_counts / expected_counts` */
  double observedOverExpected;
  double observedOverExpectedUncert;

  /** The foreground peak area (before background subtraction), and its uncertainty. */
  double foregroundCounts = 0.0;
  double foregroundUncert = 0.0;

  /** The observed counts in the peak, after subtracting background. */
  double observedCounts = 0.0;
  /** Combined uncertainty from foreground and background. */
  double observedUncert = 0.0;

  /** The expected number of counts from the model. */
  double expectedCounts = 0.0;

  /** The background counts subtracted from the observed peak area, and its uncertainty. */
  double backgroundCounts = 0.0;
  double backgroundUncert = 0.0;

  /** `peak.lineColor()` */
  Wt::WColor peakColor;
};//struct PeakResultPlotInfo
  
  
  /**
   A struct to capture the details of each source that contributed to a peak peak.
   
   This is primarily to later turn to JSON, and allow customizing log files, through inja templating.
   */
struct PeakDetailSrc
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  
  /** Energy of the nuclides gamma, in keV */
  double energy = 0.0;
  /** The number of this energy gamma, per second, for each Bq of parent nuclide. */
  double br = 0.0;
  double cpsAtSource = 0.0;
  double age = 0.0;
  
  /** For point sources, the activity of the point source.
   For self-atten sources, the activity per volume.
   For trace sources, the activity used for calculations - activity per volume, except for exponential distributions, then its activity per m2.
   */
  double calcActivity = 0.0; //Not used - can be removed
  
  /** If decay during measurement is being accounted for, then this is the rate of this gamma
   after decay correction, divided by the rate before decay correction.
   */
  double decayCorrection = 0.0;
  
  bool isTraceSource = false;//Not used - can be removed
  TraceActivityType traceSourceType = TraceActivityType::NumTraceActivityType;//Not used
  
  bool isSelfAttenSource = false;//Not used - can be removed
  
  double countsAtSource = 0.0;
  double ageUncert = 0.0;//Not used - can be removed
  //bool ageIsFit = false;
  //bool canFitAge = false;
  
  double activity = 0.0;//Not used - can be removed
  double activityUncert = 0.0;//Not used - can be removed
  double displayActivity = 0.0;//Not used - can be removed
  double displayActivityUncert = 0.0;//Not used - can be removed
  
  double massFraction = 0.0;//Not used - can be removed
  double massFractionUncert = 0.0;//Not used - can be removed
  bool isFittingMassFraction = false;//Not used - can be removed
  
  /** The expected number of counts contributed by this source, to the peak. */
  double modelContribToPeak = 0.0;
};//struct PeakDetailSrc
  
  
/** A struct to capture the details of each peak detected vs what parts of the model contributed.
 
 This is primarily to later turn to JSON, and allow customizing log files, through inja templating.
 */
struct PeakDetail
{
  double energy, decayParticleEnergy, fwhm;
  /** counts=peak.peakArea(); countsUncert=peak.peakAreaUncert();  cps=peak.peakArea()/LiveTime*/
  double counts, countsUncert, cps, cpsUncert;
  double expectedCounts, observedCounts, observedUncert, numSigmaOff;
  double observedOverExpected, observedOverExpectedUncert;
  //float modelInto4Pi, modelInto4PiCps;
  double detSolidAngle, detIntrinsicEff, detEff;

  /** Fractional (1-sigma) detector-efficiency uncertainty at this peaks energy,
   from the DRFs uncertainty info; 0 when the DRF has none or the
   `account_for_drf_uncert` option is off.  When non-zero, #observedUncert and
   #numSigmaOff include this component (added in quadrature as
   expectedCounts*drfEffFracUncert), matching the GLS-whitened fit.
   */
  double drfEffFracUncert = 0.0;

  /** Detector-efficiency validity flag at this peaks energy and the fit
   geometry (DetectorPeakResponse::EffFlag; Ok when inside the responses
   validated regime).  Non-Ok values mean the efficiency is a transfer /
   extrapolation with honestly inflated uncertainty - surfaced in the fit
   warnings and the calculation log.
   */
  DetectorPeakResponse::EffFlag drfEffFlag = DetectorPeakResponse::EffFlag::Ok;

  /** True-coincidence (cascade) summing correction applied to this peaks
   expected counts (`correct_for_cascade_summing` option).  cascadeNetMult is
   the emission-weighted net multiplier over the peaks source nuclides
   (= C_out + C_in); the out/in decomposition is per the analytic engine.
   All 1.0/0.0 when the correction was not applied.
   */
  bool cascadeCorrApplied = false;
  double cascadeNetMult = 1.0;
  double cascadeSummingOut = 1.0;
  double cascadeSummingIn = 0.0;
  
  float backgroundCounts, backgroundCountsUncert;
  
  std::string assignedNuclide;
  
  /** Further information about the sources that contribute to this peak.
   Please note that the `PeakDetailSrc` are not entirely filled out at creation, so do not rely on the member variables
   being accurate until the end of the calculations.
   */
  std::vector<PeakDetailSrc> m_sources;
  
  /** The fractional attenuation by this material (e.g., no attenuation is 1.0. Not valid for volumetric sources.
   The index of this vector is the same index as shielding in the model
   */
  std::vector<double> m_attenuations;
  
  /** The total attenuation factor (i.e., fraction no-interactions; 1.0 is no attenuation, 0.0 is total attentuation), that
   all the shieldings cause.
   Only valid for point sources.
   */
  double m_totalShieldAttenFactor;
  /** The attenuation, due to air, between the last shielding and the detector.
   For volumetric sources, this is only approximate
   */
  double m_airAttenFactor;
  
  /** The total attenuation of the source from all shielding and air.  For volumetric sources, this will only be approximate.
   */
  double m_totalAttenFactor;
  
  // TODO: for self-attenuating sources, could repeat the computation, put with only the source shell, then shell+1-other, then shell+different-1-other, and so on, to get the effect of each shell.
  
  struct VolumeSrc
  {
    bool trace; // Trace or intrinsic
    
    /** The integral of the efficiency to make it to the detector, over the source area.
     i.e., the shielding area times average efficiency to make it to the detector face.
     Does not include detector intrinsic efficiency.
     */
    double integral;
    
    /** The volume of the shielding. */
    double volume;
    
    /** `integral / volume` */
    double averageEfficiencyPerSourceGamma;
    
    /** The gammas per second, per unit-volume, for this source energy.
     This value times `VolumeSrc::integral` gives the expected number of gammas, at this energy, to strike the detector face.
     */
    double srcVolumetricActivity;
    
    bool inSituExponential;        //Not used - can be removed
    double inSituRelaxationLength; //Not used - can be removed
    
    double detIntrinsicEff;
    
    std::string sourceName;
  };//struct VolumeSrc
  
  std::vector<VolumeSrc> m_volumetric_srcs;
  
  PeakDetail()
  : energy( 0.0 ), decayParticleEnergy( 0.0 ), fwhm( 0.0 ), counts( 0.0 ),
  countsUncert( 0.0 ), cps( 0.0 ), cpsUncert( 0.0 ),
  expectedCounts( 0.0 ), observedCounts( 0.0 ), observedUncert( 0.0 ),
  numSigmaOff( 0.0 ), observedOverExpected( 0.0 ),
  //modelInto4Pi( 0.0f ), modelInto4PiCps( 0.0f ),
  detSolidAngle( 0.0 ), detIntrinsicEff( 0.0 ), detEff( 0.0 ),
  backgroundCounts( 0.0f ), backgroundCountsUncert( 0.0f ),
  assignedNuclide{}, m_sources{}, m_attenuations{},
  m_totalShieldAttenFactor( 1.0 ), m_airAttenFactor( 1.0 ),
  m_totalAttenFactor( 1.0 )
  {
  }
};//struct PeakDetail
  
  
/** A struct to capture the details of each shielding.
   
   This is primarily to later turn to JSON, and allow customizing log files, through inja templating.
*/
struct ShieldingDetails
{
  struct SelfAttenComponent
  {
    const SandiaDecay::Nuclide *m_nuclide;
    bool m_is_fit = false;
    double m_mass_frac = 0.0, m_mass_frac_uncert = 0.0;
  };
  
  struct TraceSrcDetail
  {
    const SandiaDecay::Nuclide *m_nuclide;
    GammaInteractionCalc::TraceActivityType m_trace_type = GammaInteractionCalc::TraceActivityType::NumTraceActivityType;
    bool m_is_exp_dist = false;
    double m_relaxation_length = -1.0;
  };//struct TraceSrcDetail
  
  /*
  struct PeakAttenuation
  {
    double energy;
    double attenuation;
    double incomingModelCps;
    double outgoingModelCps;
  };//struct PeakAttenuation
  */
  
  std::string m_name;
  std::string m_chemical_formula;
  
  bool m_is_generic = false;
  double m_an = 0.0, m_ad = 0.0;
  double m_density = 0.0;
  
  unsigned int m_num_dimensions = 0;
  GammaInteractionCalc::GeometryType m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
  double m_thickness = 0.0;
  double m_volume = 0.0;
  double m_volume_uncert = 0.0;
  double m_inner_rad = 0.0;
  double m_inner_dimensions[3] = {0.0, 0.0, 0.0};
  double m_outer_dimensions[3] = {0.0, 0.0, 0.0};
  bool m_fit_dimension[3] = { false, false, false };
  double m_dimension_uncert[3] = { 0.0, 0.0, 0.0 };
  
  std::vector<SelfAttenComponent> m_mass_fractions;
  std::vector<TraceSrcDetail> m_trace_sources;
  //std::vector<PeakAttenuation> m_peak_attens;
};//struct ShieldingDetails
  
  
struct SourceDetails
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  
  /** The total activity of the source. */
  double activity = 0.0;
  double activityUncertainty = 0.0;
  /** If the activity is fit for.
   
   Note: if source type is `ShieldingSourceFitCalc::ModelSourceType::Intrinsic`, then this value will be false,
   even if a shielding dimension is being fit for.
   */
  bool activityIsFit = false;
  double nuclideMass = 0.0;
  double age = 0.0;
  double ageUncertainty = 0.0;
  bool ageIsFittable = false;
  bool ageIsFit = false;
  const SandiaDecay::Nuclide *ageDefiningNuc = nullptr;
  
  bool isTraceSource;
  TraceActivityType traceActivityType = TraceActivityType::NumTraceActivityType;
  double traceSrcDisplayAct;
  double traceSrcDisplayActUncertainty;
  double traceRelaxationLength;
  
  bool isSelfAttenSource;
  size_t selfAttenShieldIndex;
  std::string selfAttenShieldName;
  bool isSelfAttenVariableMassFrac;
  /** This is the fraction of the element in the shielding, that this nuclide is for. */
  double selfAttenMassFrac;
  double selfAttenMassFracUncertainty;
  
  // We wont put peaks into this struct, but instead when we make the JSON, we'll
  //  insert peaks from `PeakDetail` as `PeakDetailSrc::nuclide` match this nuclide.
};//struct SourceDetails
  
  
/** The effective shielding seen by the gammas of one source at one energy that reach
 the detector - the inputs `GadrasShieldScatter::getContinuum` needs to estimate the
 scattered continuum (for future cascade-summing corrections).

 The averages are attenuation-weighted: each ray from the source volume is weighted by
 its `exp(-sum(mu*x))*solid_angle` contribution, so they reflect the shielding seen by
 the photons that actually reach the detector.
 */
struct EffectiveShieldingInfo
{
  const SandiaDecay::Nuclide *nuclide = nullptr;
  double energy = 0.0;

  /** Gammas of this energy emitted into 4pi over the measurement live time;
   -1 if not computed (in-situ exponentially-distributed trace sources).
   */
  double gammas_into_4pi = -1.0;

  /** Attenuation-weighted traversed areal density, in PhysicalUnits. */
  double effective_ad = 0.0;

  /** Effective atomic number, weighted by traversed mass (sum AN_l*AD_l / sum AD_l). */
  double effective_an_mass = 0.0;

  /** Effective atomic number, weighted by each materials share of the attenuation. */
  double effective_an_xs = 0.0;

  /** Hydrogen mass fraction of the traversed areal density. */
  double hydrogen_frac_of_ad = 0.0;

  bool is_point_source = false;
};//struct EffectiveShieldingInfo


/** This is the setup of the problem - distances, shielding, how to calculate. */
struct ShieldSourceConfig
{
  /** Distance to the CENTER of the item being measured.
   
   Note: If detector is fixed geometry, this value will not be used.
   */
  double distance;

  /** Lateral offsets of the source/shielding assembly from the detector axis -
   user-set, never fit.

   source_offsets[0] is the horizontal offset perpendicular to the detector axis,
   source_offsets[1] the vertical one; for Spherical geometry only the magnitude
   matters (the second component is ignored by the UI).  In assembly-centered
   coordinates the detector sits at {-dx, -dy, distance} (CylinderSideOn:
   {distance, -dx, -dy}) - see #detector_geom_from_config.

   Not used for fixed-geometry detectors.  Serialized as an optional
   <SourceOffsets dx=".." dy=".."> node (absent means zero).
   */
  double source_offsets[2] = { 0.0, 0.0 };
  
  GammaInteractionCalc::GeometryType geometry;
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> shieldings;
  std::vector<ShieldingSourceFitCalc::SourceFitDef> sources;
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions options;
  
  rapidxml::xml_node<char> *serialize( rapidxml::xml_node<char> *base_node ) const;
  void deSerialize( const rapidxml::xml_node<char> *base_node );
};//struct ShieldSourceCalcInput

  
class ShieldingSourceChi2Fcn
    : public ROOT::Minuit2::FCNBase
{
//This class evaluated the chi2 of a given hypothesis, where it is assumed the
//  radioactive source is a point source located at the center of concentric
//  spherical shells consisting of non-radioactive materials, that may either
//  be defined via the 'Material' class, or be a generic material defined by
//  atomic number and areal density which is indicated by a NULL pointer.
//
//parameters:
//-Activity (in MBq) nuclide 0, (nuclides are sorted alphabetically by name)
//-Age of nuclide 0
//-Activity (in MBq) nuclide 1
//-Age of nuclide 1 (if negative, must be negative value of one plus index of
//                   defining nuclide; e.g., a negative int.  Hacky, but whatever.
//                   To get age for this case, do `age = 1 + 2*((-1*parameter) - 1)` )
// ...
//if material 0 normal material (if Material* is non-NULL pointer)
//    -Material { spherical thickness | cylindrical radius thickness | rectangular width thickness }
//    -Material { ignored | cylindrical length thickness | rectangular height thickness }
//    -Material { ignored | ignored | rectangular depth thickness }
//else if generic material (if Material* is NULL)
//    -atomic number
//    -areal density  (in units of PhysicalUnits,
//                     e.g. to print out to user divide by g/cm2)
//    - ignored
//if material 1 normal material (if Material* is non-NULL pointer)
//    -Material { spherical thickness | cylindrical radius thickness | rectangular width thickness }
// ...
//
//could add another member variable that holds pointer to source isotopes to
//  fit for mass fractions of.  Then there would be M-1 additional parameters,
//  where M is the number of isotopes to fit for in a material.
//  Mass fraction is then ...
  
public:
  /** In order to keep numbers roughly where Minuit2 can handle them, we have to work in units of
   1.0E6 becquerel.  E.g., for like plutonium problems, some of the branching rations of 1E-11 are
   still significant to the problem, so to avoid losing accuracy in decay calculations, we will
   multiply things  a bit.
   
   \sa ns_decay_act_mult in RelActCalcAuto.cpp
   */
  static const double sm_activityUnits;  //SandiaDecay::MBq


  
  /** Options controlling the supplemental per-peak information (Currie detection-limit checks,
   and nominal implied activities for peaks not used in the fit) that
   `ShieldingSourceFitCalc::fit_model(...)` computes after a successful fit.

   These are runtime-only options - they are not serialized with the model.
   */
  struct SupplementalInfoOptions
  {
    /** Whether to compute the supplemental info at all. */
    bool compute = true;

    /** Confidence level for the Currie-style limits; e.g., 0.95 for 95%. */
    double confidence_level = 0.95;

    /** Number of channels on either side of the peak region used to estimate the continuum. */
    size_t num_side_channels = 4;

    /** Width of the peak region, in multiples of the peak FWHM.
     The default of 2.5 (i.e., +-1.25 FWHM) is what ISO 11929:2010 recommends.
     */
    double roi_num_fwhm = 2.5;

    /** Whether activities in the descriptions should be written in curie, rather than becquerel. */
    bool use_curie = true;
  };//struct SupplementalInfoOptions


  struct ShieldSourceInput
  {
    ShieldSourceConfig config;

    std::shared_ptr<const DetectorPeakResponse> detector;
    std::deque<std::shared_ptr<const PeakDef>> foreground_peaks;
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> background_peaks;

    std::shared_ptr<const SpecUtils::Measurement> foreground;
    std::shared_ptr<const SpecUtils::Measurement> background;

    /** The background scale factor for the realtive peak areas.
     Normally this is the live time dormalization factor (e.g., background real time divided by foreground real time), but the user may have
     scaled the background for some bespoke reason, so we need to respect this for subtracting background peak areas.

     If zero, negative, inf, or NaN, will use live-time scaling.
     */
    double background_sf = -1.0;

    /** Options for the supplemental per-peak info computed after the fit; see the struct. */
    SupplementalInfoOptions supplemental_options;

    /** Peaks in `foreground_peaks` that do not correspond to a peak actually observed/fit in the
     spectrum (e.g., exemplar peaks a batch analysis could not fit).  Their fitted amplitude is
     meaningless, but including them lets the supplemental info computation predict the counts the
     fitted model implies for them (e.g., to convert a counts-based detection limit to activity).
     They must never have `useForShieldingSourceFit()` set.
     */
    std::vector<std::shared_ptr<const PeakDef>> synthetic_peaks;

    /** An already-built cascade calculator to take the (expensive) cascade enumeration and
     shield-scatter table from, rather than rebuilding them.

     Only the peak windows are rebuilt from this input's peaks, so it is safe to pass one built
     for a different peak list - but everything else (source nuclides and their starting ages,
     detector, shield stack) must match, since none of it is re-derived.  Used by
     `ShieldingSourceFitCalc::compute_supplemental_peak_info(...)`, which re-creates the function
     object with extra peaks right after a fit.  Ignored unless
     `config.options.correct_for_cascade_summing` is set.

     Note that when `supplemental_options.compute` is set, the cascade calculator is given a
     window for every peak in `foreground_peaks` rather than only the fitted ones - so the peaks
     the post-fit pass adds are already covered either way.
     */
    std::shared_ptr<const CascadeSummingCalc> reuse_cascade_calc;
  };//struct ShieldSourceInput

  static std::pair<std::shared_ptr<ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> create(
                                     const ShieldSourceInput &input );
  
protected:
  ShieldingSourceChi2Fcn(
                      const double distance,
                      const double liveTime,
                      const double realTime,
                      const std::vector<PeakDef> &peaks,
                      std::shared_ptr<const DetectorPeakResponse> detector,
                      const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                      const GeometryType geometry,
                      const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options );

public:
  virtual ~ShieldingSourceChi2Fcn();

  /** Returns the geometry of this ShieldingSourceChi2Fcn */
  const GeometryType geometry() const;
  
  /** Causes exception to be thrown if DoEval() is called afterwards. */
  void cancelFit();
  
  /** Similar to #cancelFit, but status will be set to #CalcStatus::CanceledNoUpdate */
  void cancelFitWithNoUpdate();
  

  
  /** Information tracked during test evaluation of the Chi2.  Set this object
   using #setGuiProgressUpdater so you can properly bind things to the GUI.
   */
  struct GuiProgressUpdateInfo
  {
    GuiProgressUpdateInfo( const size_t updateFreqMs,
                          std::function<void(size_t ncalls, double elapsed_time, double chi2, std::vector<double> pars)> updater );
    
    void fitting_starting();
    void completed_eval( const double chi2, const std::vector<double> &pars );
    
    size_t numFunctionCallsSoFar();
    double bestChi2SoFar();
    std::vector<double> bestParametersSoFar();
    
  private:
    /** A mutex that protects all member variables.
     Lock is taken in #fitting_starting and #completed_eval functions, and while m_gui_updater is being called.
     */
    std::mutex m_mutex;
    
    /** The function called periodically to update the best chi2 found. */
    const std::function<void(size_t, double, double, std::vector<double>)> m_gui_updater;
    
    const size_t m_update_frequency_ms;
    
    size_t m_num_fcn_calls;
    
    double m_fitStartTime;
    double m_currentTime;
    double m_lastGuiUpdateTime;
    
    double m_bestChi2;
    std::vector<double> m_bestParameters;
  };//struct GuiProgressUpdateInfo
  
  /** Call this function to have the Chi2 call the specified callback with the
      best solution found so far, at the specified intervals.  The callback will
      be posted to be executed in the specified Wt session using WServer.
   */
  void setGuiProgressUpdater( std::shared_ptr<GuiProgressUpdateInfo> updateInfo );
  
  /** Call this function *just* before starting to fit; it will set the zombie
      timer for the specified delay.  The zombie timer is how long the fit can
      go on before its considered a failure and is aborted.
      If the guiProgressUpdater() has bee set, then the best chi2 seen so far
      and, the timer, and parameters are reset.
   */
  void fittingIsStarting( const size_t deadlineMs );
  
  /** Call this function as soon as fitting is done - it cancels the zombie
      timer.
   */
  void fittingIsFinished();
  
  /** Sets whether to use a SpecUtilsAsync::ThreadPool when calculating contributions for each peak from self-attenuating and/or
   trace sources.
   If you are doing multiple parallel fits, you may want to disable multithread to better use the cpu.
   
   Default is true.
   
   \sa ShieldingSourceFitCalc::ShieldingSourceFitOptions::multithread_self_atten
   */
  void setSelfAttMultiThread( const bool do_multithread );
  
  /** Used to set initial source definitions, and add the needed fitting parameters.
   
   Returns number of fitting parameters added.
   */
  size_t setInitialSourceDefinitions( const std::vector<ShieldingSourceFitCalc::SourceFitDef> &src_definitions,
                                     const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                     ROOT::Minuit2::MnUserParameters &inputPrams );
  
  const std::vector<ShieldingSourceFitCalc::SourceFitDef> &initialSourceDefinitions() const;
  
  /* Returns the mass-fraction of the specified nuclide, in the specified shielding.
 
   `errors` must either be empty (in which case uncert will be set to zero), or the
   same size as `pars`.
 
   Note: the returned mass fraction is mass fraction for the nuclides element - and not of the
   entire shielding.
   
   Nuclide may, or may not, have its mass-fraction being fitted for.
   
   The element is required when nuc is nullptr, which indicates the mass fraction of the "other"
   non-source nuclides.
   */
  void massFractionOfElement( double &massFrac, double &uncert,
                     const size_t material_index,
                     const SandiaDecay::Nuclide *nuc,
                     const SandiaDecay::Element *el,
                     const std::vector<double> &pars,
                     const std::vector<double> &errors ) const;
  
  double massFractionOfElement( const size_t material_index,
                          const SandiaDecay::Nuclide *nuc,
                          const std::vector<double> &pars ) const;
  double massFractionOfElementUncertainty( const size_t material_index,
                             const SandiaDecay::Nuclide *nuc,
                             const std::vector<double> &pars,
                             const std::vector<double> &error ) const;
  
  /** Return if the passed in nuclide is having its mass-fraction fitted.
   */
  bool isVariableMassFraction( const size_t material_index,
                               const SandiaDecay::Nuclide *nuc ) const;
  /** Returns if the elements "other nuclide" fraction is being fit. */
  bool isVariableOtherMassFraction( const size_t material_index,
                               const SandiaDecay::Element *el ) const;
  
  bool hasVariableMassFraction( const size_t material_index ) const;
  
  std::vector<std::shared_ptr<const Material>> materialsFittingMassFracsFor() const;
  
  /** Returns nuclides that are self-attenuating, wether fitting the mass-fraction for them or not. */
  std::vector<const SandiaDecay::Nuclide *> selfAttenuatingNuclides( const size_t material_index ) const;
  
  /** Returns the nuclides fitting mass fraction for, grouped by element.
   
   nullptr nuclides represent the "other" non-source nuclides, if that component is being fit.
   
   The mass-fraction parameter order is dictated by this result (which is also the same as `m_initial_shieldings`), of looping
   over the Element, and then all the nuclides for that element ordered as in the vector..
   */
  std::map<const SandiaDecay::Element *,std::vector<const SandiaDecay::Nuclide *>> nuclideFittingMassFracFor(
                                                                const size_t material_index ) const;
  /** */
  std::vector<const SandiaDecay::Element *> elementsFittingOtherFracFor( const size_t material_index ) const;

  /** Returns information about all self-attenuating sources, grouped by element, wether mass-fraction for them are being fit or not.
   
   A nullptr nuclide indicates "other" non-source component, and may or may not be present if not fit for (if fit for, it will be present).
   
   Information returned about each nuclide is mass-fraction, uncertainty (valid if >0), and if it was fit for.
   */
  std::map<const SandiaDecay::Element *,std::vector<std::tuple<const SandiaDecay::Nuclide *,double,double,bool>>> selfAttenSrcInfo(
                                                            const size_t material_index,
                                                            const std::vector<double> &pars,
                                                            const std::vector<double> &error ) const;
  
  //setBackgroundPeaks(...): if you wish to correct for background counts, you
  //  can set that here.  The peaks you pass in should be the original
  //  background peaks for the background; similarly for the live time.
  //  The `scale_factor` is normally the "live time" normalization factor
  //  (e.g. foreground live time divided by background live time), but sometimes users use a custom
  //  scale factor (wrong or missing live-times usually), so we will use this instead; if this is
  //  zero, negative, inf, or NaN, we'll use live-time normalization to scale the peaks by.
  void setBackgroundPeaks( const std::vector<PeakDef> &peaks, const double liveTime, const double scale_factor );

  
  /** The calculation status for ShieldingSourceChi2Fcn. */
  enum class CalcStatus : int
  {
    NotCanceled = 0,
    UserCanceled = 1,
    Timeout = 2,
    CanceledNoUpdate = 3
  };//enum class CalcStatus : int
  
  
  /** Exception thrown from #DoEval when fitting is canceled by the user, or times-out.
   Useful for propagating reason
   
   */
  class CancelException : public std::exception
  {
  public:
    CancelException( const CalcStatus cancel_code );
    
    CalcStatus m_code;
  };//class CancelException
  
  /** Performs evaluation of Chi2, for parameters x.

   May through CancelException (if user or time limit cancelled computation), or other std::exception (on other error type).
   */
  virtual double DoEval( const std::vector<double> &x ) const;

  /** The current cancel/timeout status - lets external optimizer drivers (e.g., the
   Ceres cost function) poll for cancellation, instead of relying on #DoEval throwing.
   */
  CalcStatus currentCancelStatus() const;

  /** Reports a completed chi2 evaluation to the GUI progress updater (if a fit is in
   progress and an updater was set) - the same reporting #DoEval does internally, for
   use by optimizer drivers that dont evaluate through #DoEval.
   */
  void reportCompletedEval( const double chi2, const std::vector<double> &params ) const;

  
  /** For interface compatibility; calls directly to #DoEval */
  virtual double operator()( const std::vector<double> &x ) const;
  
  
  /** Gives the chi2 contributions for each peak
   
   @param params The parameters currently defining the model.
   @param error_params The errors on the parameters - only used if `log_info` is non-null
   @param mixturecache Used to speed up multiple calls to this function, and may be empty at first.
   @param log_info If non-null, filled out with details about the computation.
   */
  typedef std::map< const SandiaDecay::Nuclide *, SandiaDecay::NuclideMixture> NucMixtureCache;
  std::vector<PeakResultPlotInfo> energy_chi_contributions(
                                  const std::vector<double> &params,
                                  const std::vector<double> &error_params,
                                  NucMixtureCache &mixturecache,
                                  std::vector<PeakDetail> *log_info = nullptr ) const;

  /** Templated (double or ceres::Jet<>) core of the expected-counts computation,
   used by the auto-differentiated (Ceres) fit; definitions are in
   GammaInteractionCalc_imp.hpp, which callers must include (after ceres
   headers, when using Jets).

   Computes the model-expected counts for each peak that
   #expected_observed_chis would include (i.e., peaks with a decay particle or
   annihilation gamma), in that same order, with derivative information
   propagated through shielding transmission, volumetric-source integration
   (via #DistributedSrcCalcT and #self_shielding_integration_imp - NOT Cuhre),
   and - by numeric differencing re-seeded into the Jet - nuclide age.

   No logging; unlike #DoEval, exceptions (including CancelException) propagate.
   */
  template<typename T>
  std::vector<T> expected_peak_counts_imp( const std::vector<T> &params,
                                           NucMixtureCache &mixturecache ) const;

  /** Builds the volumetric-source (self-attenuating and trace) calculators for the
   given parameters - one per (source nuclide, clustered gamma energy) - with their
   `m_srcVolumetricActivity` set, ready to integrate.  The per-shell composition
   metadata (density, effective AN, hydrogen fraction) is also filled in.

   Used by #expected_peak_counts_imp, and by #computeEffectiveShielding; defined in
   GammaInteractionCalc_imp.hpp.
   */
  template<typename T>
  /** `log_info`, when non-null (only meaningful for the T=double instantiation), routes the
   per-source clustering through the logging #cluster_peak_activities overload so the
   `PeakDetailSrc` entries (cpsAtSource / countsAtSource) get filled.  That is the only thing the
   display path needs from this that the fit path does not - so both can share one builder rather
   than keeping a second, cascade-blind copy of the volumetric model.
   */
  std::vector<std::unique_ptr<DistributedSrcCalcT<T>>> build_volumetric_calculators(
                          const std::vector<T> &params,
                          NucMixtureCache &mixturecache,
                          const std::vector<std::pair<double,double>> &energie_widths,
                          std::vector<PeakDetail> *log_info = nullptr ) const;

  /** The user-set lateral offsets of the source assembly from the detector axis -
   see #ShieldSourceConfig::source_offsets.
   */
  double sourceOffsetX() const;
  double sourceOffsetY() const;

  /** The true line-of-sight distance from the assembly center to the detector,
   sqrt( distance^2 + offsetX^2 + offsetY^2 ).
   */
  double trueSourceToDetectorDistance() const;

  /** Computes the effective (attenuation-weighted) shielding each sources gammas
   traverse, per gamma energy of the fit peaks - the inputs needed for GADRAS
   shield-scatter based cascade-summing corrections.

   A post-fit diagnostic; call with the final fit parameters.  Defined in
   GammaInteractionCalc_imp.hpp.
   */
  std::vector<EffectiveShieldingInfo> computeEffectiveShielding( const std::vector<double> &params,
                                                                 NucMixtureCache &mixturecache ) const;

  /** Templated equivalents of the same-named double-valued accessors below;
   defined in GammaInteractionCalc_imp.hpp.
   */
  template<typename T>
  T activity_imp( const SandiaDecay::Nuclide *nuclide, const std::vector<T> &params ) const;
  template<typename T>
  T age_imp( const SandiaDecay::Nuclide *nuclide, const std::vector<T> &params ) const;
  template<typename T>
  T activityOfSelfAttenSource_imp( const SandiaDecay::Nuclide *nuclide, const std::vector<T> &params ) const;
  template<typename T>
  T massFractionOfElement_imp( const size_t material_index, const SandiaDecay::Nuclide *nuc,
                               const std::vector<T> &params ) const;
  template<typename T>
  T volumeOfMaterial_imp( const size_t matn, const std::vector<T> &params ) const;
  template<typename T>
  T sphericalThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T cylindricalRadiusThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T cylindricalLengthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T rectangularWidthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T rectangularHeightThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T rectangularDepthThickness_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T arealDensity_imp( const size_t materialNum, const std::vector<T> &params ) const;
  template<typename T>
  T atomicNumber_imp( const size_t materialNum, const std::vector<T> &params ) const;

  /** Templated equivalent of #cluster_peak_activities (no logging).

   When T is a ceres::Jet and `age` carries derivative information, the
   derivative of each gammas yield with respect to age is computed by numeric
   differencing (the decay calculations are not templatable), and chained into
   the result - same approach as RelActCalcAuto.
   */
  template<typename T>
  static void cluster_peak_activities_imp( std::map<double,T> &energy_count_map,
                  const std::vector<std::pair<double,double>> &energie_widths,
                  SandiaDecay::NuclideMixture &mixture,
                  const T &act, const T &age,
                  const double photopeakClusterSigma,
                  const double energyToCluster,
                  const bool accountForDecayDuringMeas,
                  const double measDuration );

  void log_shield_info( const std::vector<double> &params,
                        const std::vector<double> &error_params,
                        const std::vector<ShieldingSourceFitCalc::SourceFitDef> &fit_src_info,
                        std::vector<ShieldingDetails> &info ) const;
  
  void log_source_info( const std::vector<double> &params,
                        const std::vector<double> &error_params,
                        const std::vector<ShieldingSourceFitCalc::SourceFitDef> &fit_src_info,
                        std::vector<SourceDetails> &info ) const;
  
  
  ShieldingSourceChi2Fcn&	operator=( const ShieldingSourceChi2Fcn & );
  virtual double Up() const;

  size_t numExpectedFitParameters() const;

  //nuclide(...) throws std::runtime_error if an invalid number is passed in
  const SandiaDecay::Nuclide *nuclide( const size_t nuclide_index ) const;

  /** Returns activity for a given nuclide.
   
   If a self-atten source, will be calculated using `activityOfSelfAttenSource(...)`.
   If a trace source, will return activity in the units the trace source is defined (so total act, or per cm3, or m2, or per gram).
   
   Will throw `runtime_exception` if params is wrong size or invalid nuclide passed in.
   */
  double activity( const SandiaDecay::Nuclide *nuclide,
                   const std::vector<double> &params ) const;
  double activityUncertainty( const SandiaDecay::Nuclide *nuclide,
                  const std::vector<double> &params,
                  const std::vector<double> &errors ) const;
  /**
   
   Will throw `runtime_exception` if params is wrong size or invalid nuclide passed in.
   */
  double age( const SandiaDecay::Nuclide *nuclide,
                   const std::vector<double> &params ) const;
  size_t nuclideIndex( const SandiaDecay::Nuclide *nuclide ) const;

  /** Returns whether or not the nuclide is a self-attenuating source. */
  bool isSelfAttenSource( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns whether or not the nuclide is a trace source. */
  bool isTraceSource( const SandiaDecay::Nuclide *nuc ) const;
  
  std::vector<const SandiaDecay::Nuclide *> traceNuclidesForMaterial( const size_t material_index ) const;
  
  /** Returns the trace source activity type for a nuc.
   
   Throws exception if nuc is not a trace source.
   */
  GammaInteractionCalc::TraceActivityType traceSourceActivityType(
                                                          const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns relaxation length for a nuclide.
   
   Throws exception if nuc is not a trace source, or not TraceActivityType::ExponentialDistribution.
   */
  double relaxationLength( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns whether or not the nuclide is a self-attenuating OR trace source. */
  bool isVolumetricSource( const SandiaDecay::Nuclide *nuc ) const;
  
  /** Returns the volume of a material; note does not include the volume of inner shieldings. */
  double volumeOfMaterial( const size_t materialN, const std::vector<double> &params ) const;
  
  /** Returns the uncertainty on the volume, taking inner dimensions, outer dimensions, and all dimensions to be independent.  */
  double volumeUncertaintyOfMaterial( const int materialN, const std::vector<double> &params,
                                      const std::vector<double> &errors ) const;
  
  /** Returns the activity of the specified nuclide that is a self-attenuating source.
   
   Throws exception if nuclide is nullptr, or not a self-attenuating source.
   */
  double activityOfSelfAttenSource( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const;
  
  /** Returns the total activity of the specified nuclide.
   
   Note: for trace sources, the #activity function returns the display activity (so either total, per cc, or per g) for the trace sources.
   
   Throws exception if nuclide is nullptr.
   */
  double totalActivity( const SandiaDecay::Nuclide *nuc, const std::vector<double> &params ) const;
  
  /** Similar to totalActivity, but returns its uncertainty. */
  double totalActivityUncertainty( const SandiaDecay::Nuclide *nuc,
                                   const std::vector<double> &params,
                                   const std::vector<double> &paramErrors ) const;
  
  
  size_t numNuclides() const;
  size_t numMaterials() const;
  
  bool isSpecificMaterial( const size_t materialNum ) const;
  bool isGenericMaterial( const size_t materialNum ) const;
  
  //material(): will throw exception if invalid materialNum, and will return
  //  nullptr if a generic material.
  std::shared_ptr<const Material> material( const size_t materialNum ) const;
  
  //sphericalThickness(...): will throw std::runtime_exception if material is a generic
  //  material
  double sphericalThickness( const size_t materialNum,
                    const std::vector<double> &params ) const;

  double cylindricalRadiusThickness( const size_t materialNum,
                            const std::vector<double> &params ) const;
  double cylindricalLengthThickness( const size_t materialNum,
                                    const std::vector<double> &params ) const;
  double rectangularWidthThickness( const size_t materialNum,
                                    const std::vector<double> &params ) const;
  double rectangularHeightThickness( const size_t materialNum,
                                   const std::vector<double> &params ) const;
  double rectangularDepthThickness( const size_t materialNum,
                                    const std::vector<double> &params ) const;
  
  
  //arealDensity(...): will throw std::runtime_exception if material is a
  //  specific material
  double arealDensity( const size_t materialNum,
                       const std::vector<double> &params ) const;

  //atomicNumber(...): will throw std::runtime_exception if material is a
  //  specific material
  double atomicNumber( const size_t materialNum,
                       const std::vector<double> &params ) const;

  const std::vector<PeakDef> &peaks() const;

  /** Note: the `setBackgroundPeaks(...)` function has scaled these peaks by the background normalization factor. */
  const std::vector<PeakDef> &backgroundPeaks() const;
  double backgroundNormalizationFactor() const;

  /** Energies (keV) of the peaks included in the fit - the same inclusion
   rule and ordering as the fit residuals, so per-peak quantities line up
   one-to-one with them.
   */
  std::vector<double> includedPeakEnergies() const;

  /** Row-major NxN fractional detector-efficiency covariance among the
   included peaks (see #includedPeakEnergies), evaluated on-axis at the fit
   distance; strongly correlated between nearby energies for grounded
   Monte-Carlo-parameterized responses.  Empty if the detector has no
   uncertainty information.

   Not yet folded into the fit residuals - available for the chi2 / reported
   parameter covariance to consume (a common-mode efficiency error maps ~1:1
   onto activity and must not be averaged down by sqrt(num-peaks)).
   */
  std::vector<double> peakEffFracCovariance() const;

  /** Square root of the diagonal of #peakEffFracCovariance - the per-peak
   1-sigma fractional efficiency uncertainty envelope.  Empty if the detector
   has no uncertainty information.
   */
  std::vector<double> peakEffFracUncerts() const;


  /** Detector-efficiency validity flag for each included peak (same ordering
   as #includedPeakEnergies), evaluated at the fit geometry (distance and any
   source offset) - i.e. the DetectorPeakResponse::EffEval::flag the forward
   model sees.  Non-Ok flags mark queries outside the response's validated
   regime (near-field below the validity floor, refuse-grade off-axis,
   collimator shadowing, energy clamping) and should be surfaced as fit
   warnings.  Empty for fixed-geometry detectors or when no detector is set.
   */
  std::vector<std::pair<double,DetectorPeakResponse::EffFlag>> peakDrfEffFlags() const;

  double distance() const;

  double liveTime() const;
  double realTime() const;

  /** The input this object was created from, kept verbatim.

   Notably `ShieldSourceInput::foreground_peaks` holds _all_ the peaks that were given, including
   the ones not used in the fit; `peaks()`, in contrast, holds only the peaks the fit uses.

   Copying this, changing which peaks are flagged `useForShieldingSourceFit()`, and calling
   `create(...)` again gives a function object with an identical parameter layout, so a fitted
   parameter vector can be applied to it directly.
   */
  const ShieldSourceInput &createInput() const;

  /** The cascade-summing calculator built for this object, or null when the correction is off.
   Hand this to `ShieldSourceInput::reuse_cascade_calc` when re-creating for the same scene.
   */
  const std::shared_ptr<const CascadeSummingCalc> &cascadeCalc() const;

  const std::shared_ptr<const DetectorPeakResponse> &detector() const;

  const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options() const;

  /** The volumetric-efficiency method actually used, after resolving
   `options().volumetric_eff_method` against the DRF - never #VolumetricEffMethod::Auto.
   Reported so a fit that silently fell back (e.g. to flat-disk) is visible in the results.
   */
  ShieldingSourceFitCalc::VolumetricEffMethod resolvedVolumetricEffMethod() const;

  /** Human-readable note on how #resolvedVolumetricEffMethod was arrived at (e.g.
   "Auto -> MC transfer", or a fallback reason); empty when the request was used as-is.
   */
  const std::string &volumetricEffResolveNote() const;

  /** Non-empty when an EXPLICITLY requested near-field method (MCTransfer / EffTran) could not be
   honored and the fit fell back to flat-disk.  Distinct from #volumetricEffResolveNote: a note is
   informational (including every `Auto` resolution), whereas this says the user asked for something
   specific and did not get it, and is surfaced in `ModelFitResults::errormsgs`.
   */
  const std::string &volumetricEffResolveError() const;

  /** Which volumetric-efficiency method `requested` resolves to for `drf`, considering only what the
   DRF can offer - no transfer response is built, nothing is cached, and it is safe to call on every
   UI change.

   Shared by #resolveVolumetricEffMethod and by the GUI's "-> EFFTRAN transfer" status text, so the
   status cannot claim one thing while the fit does another.  The fit may still downgrade further if
   *building* the chosen response throws, which this (deliberately) cannot predict.

   `Auto` policy: the transfer is strictly the more accurate model, so it is preferred by default and
   only stepped down to flat-disk where flat-disk is *known* to be right - at or beyond the distance
   the DRF's efficiency was characterized at (see `CeeLoUtils::pinnedEfficiencyDistanceCm`), because
   there the flat-disk extrapolation is precisely what the measured curve already encodes.

   "Far enough" is judged against everything that sets the angular scale, not distance alone:
     - at or beyond the DRF's characterization distance,
     - on-axis (zero lateral offset),
     - and far compared to the DETECTOR's own half-width plus the SOURCE's transverse half-extent,
       using the response's near-field gate (`SigmaFloors::near_regime_a`, d < a*that is near).
   Off-axis eta(E,theta) is not a near-field effect and does not fall off with distance, so a wide
   source still subtends real angles however far away it is - which is why its extent enters here.

   Bias is deliberately one-way: any input that is unknown or non-positive keeps the correction.  A
   needless transfer only costs time, whereas a wrongly-chosen flat-disk is a silent accuracy loss.

   @param source_distance Source-to-detector distance (PhysicalUnits); <= 0 when unknown.
   @param off_axis_offset Radial offset of the source from the detector axis (PhysicalUnits); 0 is
          on-axis.
   @param source_half_extent Transverse half-extent of the source assembly (PhysicalUnits), from
          #assemblyTransverseHalfExtent; <= 0 when unknown.
   */
  static ShieldingSourceFitCalc::VolumetricEffMethod resolveVolumetricEffMethodForDrf(
                      const std::shared_ptr<const DetectorPeakResponse> &drf,
                      const ShieldingSourceFitCalc::VolumetricEffMethod requested,
                      const double source_distance,
                      const double off_axis_offset,
                      const double source_half_extent );

  /** Largest transverse (perpendicular to the detector axis) half-extent of the source assembly, in
   PhysicalUnits - i.e. how far off the axis a source element can possibly sit.

   Per-layer `m_dimensions` are thicknesses (half-dimensions for the innermost layer), so this sums
   them; generic layers have no physical extent and are skipped.  Transverse means radial for a
   sphere or an end-on cylinder, but spans the length for a side-on cylinder and the width/height for
   a box, since the detector axis differs.

   Shared by the fit and the GUI so both feed #resolveVolumetricEffMethodForDrf the same number.
   */
  static double assemblyTransverseHalfExtent(
                      const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                      const GeometryType geometry );

  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &initialShieldings() const;
  

  //observedPeakEnergyWidths(): get energy sorted pairs of peak means and widths
  static std::vector< std::pair<double,double> > observedPeakEnergyWidths(
                                                              const std::vector<PeakDef> &peaks );

  //cluster_peak_activities(...): clusters the number of decays per second
  //  by energies.  If photopeakClusterSigma>0.0, then all gamma lines nearby
  //  a peaks mean will be considered to contribute to that peak.  'Nearby'
  //  is defined by 'energie_widths' - the energy and sigma of fit peaks.
  //  The same 'energie_widths' must be used with all subsequent calls using
  //  the same 'energy_count_map' (this is unchecked, so dont violate it).
  //  energie_widths should consist of the true photopeak energies, and not
  //  the detected energy.
  //  It is also assumed that the Nuclides have been added to the mixture
  //  with an activity of 1.0*sm_activityUnits (this will be divided out then
  //  multiplied by 'act'), at an initial age of 0.0.
  //  There must be exactly one parent nuclide in 'mixture' or an exception will
  //  be thrown.
  //  If energyToCluster is > 0.0, then only photopeaks within
  //  photopeakClusterSigma in mixture will be clustered and added to
  //  energy_count_map.  If energyToCluster <= 0.0, then all photopeaks in
  //  mixture will be clustered and added to energy_count_map.
  // @param accountForDecayDuringMeas If true, a correction will be made for the decay
  //        of the nuclides during the measurement (i.e., for nuclides without prodginy,
  //        the total counts for each energy that are output will be lower than without
  //        this correction, causing data to fit a higher activity; if prodginy are
  //        involved, the correction could go in either direction), with the specified
  //        activity cooresponding to the begining of the measurement.
  // @param measDuration The duration of the measurement, in seconds; only used if
  //        `accountForDecayDuringMeas` is also used.
  static void cluster_peak_activities( std::map<double,double> &energy_count_map,
                  const std::vector< std::pair<double,double> > &energie_widths,
                  SandiaDecay::NuclideMixture &mixture,
                  const double act,
                  const double thisAge,
                  const double photopeakClusterSigma,
                  const double energyToCluster,
                  const bool accountForDecayDuringMeas,
                  const double measDuration,
                  std::vector<GammaInteractionCalc::PeakDetail> *log_info
              );


  //returns the chi computed from the expected verses observed counts; one
  //  chi2 for each peak energy.  Each returned entry is {energy,chi,scale,PeakColor,ScaleUncert},
  //  where scale is observed/expected
  //
  //  `eff_frac_uncerts`, when non-null, holds the per-peak fractional detector-
  //  efficiency uncertainties (same inclusion rule and order as
  //  #includedPeakEnergies, i.e. one entry per included peak); each peaks
  //  denominator is then inflated to sqrt(stat^2 + (expected*frac)^2) so the
  //  displayed chi agrees with the GLS-whitened (Ceres) fit that accounts for
  //  the DRF uncertainty.  Pass null (or when the DRF has no uncertainty info)
  //  for the historical statistics-only behavior.
  //  `eff_flags`, when non-null, holds the per-peak detector-efficiency
  //  validity flags (from #peakDrfEffFlags, same inclusion rule and order) to
  //  record on the PeakDetail log entries.
  //  `eff_frac_uncerts` governs the displayed observed-uncertainty column, and
  //  hence `numSigmaOff`, which is the MARGINAL pull
  //  (obs - exp)/sqrt(stat^2 + obs^2*C_ii) - a genuine per-peak sigma taken from
  //  the diagonal of the same covariance the correlated fit minimizes.  A GLS
  //  whitened residual would NOT be usable here: it mixes peaks in Cholesky
  //  order, so it is not a property of any one peak and cannot be read as a
  //  pull-vs-energy trend (see ShieldSourcePullTrend).  The full whitening
  //  matrix is still reported, on `ModelFitResults::efficiency_whitening`.
  static std::vector<PeakResultPlotInfo> expected_observed_chis(
                              const std::vector<PeakDef> &peaks,
                              const std::vector<PeakDef> &backgroundPeaks,
                              const std::map<double,double> &energy_count_map,
                              std::vector<GammaInteractionCalc::PeakDetail> *log_info = nullptr,
                              const std::vector<double> *eff_frac_uncerts = nullptr,
                              const std::vector<std::pair<double,DetectorPeakResponse::EffFlag>> *eff_flags = nullptr );
protected:
  
  void zombieCallback( const boost::system::error_code &ec );


protected:
  
  std::atomic<CalcStatus> m_cancel;
  
  //Used to determine if m_guiUpdateInfo related stuff should be checked on
  std::atomic<bool> m_isFitting;
  
  std::shared_ptr<GuiProgressUpdateInfo> m_guiUpdateInfo;
  
  //Sometimes self attenuating fits can run-away and never end - protect against it
  std::mutex m_zombieCheckTimerMutex;
  std::shared_ptr<boost::asio::deadline_timer> m_zombieCheckTimer;
  
  double m_distance;

  /** Lateral offsets of the source assembly from the detector axis - see
   #ShieldSourceConfig::source_offsets.  Set by #create; never fit.
   */
  double m_sourceOffsets[2] = { 0.0, 0.0 };
  double m_liveTime;

  std::vector<PeakDef> m_peaks;

  /** The live time normalization factor (or user overrided value) for the background spectrum. */
  double m_backgroundPeakScale;

  /** The peak area and area uncertainty have been scaled by `m_backgroundPeakScale` (see `setBackgroundPeaks(...)` function).
   TODO: stop scaling peaks by this factor, and just use this factor during the computations.
   */
  std::vector<PeakDef> m_backgroundPeaks;


  std::shared_ptr<const DetectorPeakResponse> m_detector;

  /** The input `create(...)` was given, kept verbatim.

   The chi2 evaluation itself does not use this - it is kept so that after the fit, the per-peak
   supplemental information (detection-limit checks, and the activities implied by peaks not used in
   the fit) can be computed by building a second function object over an augmented peak list.
   Keeping the input exactly as it came in is what guarantees the second object lays its parameters
   out identically, so the fitted parameter values can be applied to it directly.
   */
  ShieldSourceInput m_create_input;

  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_initial_shieldings;
  
  std::vector<const SandiaDecay::Nuclide *> m_nuclides; //sorted alphabetically and unique
  
  const GeometryType m_geometry;
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions m_options;

  /** Cascade (true-coincidence) summing corrections; non-null only when
   `m_options.correct_for_cascade_summing` is set (see
   ShieldingSourceChi2Fcn::create, which validates the DRF has the needed
   total-efficiency info and throws otherwise).  Immutable after creation;
   safe for the fit worker threads.
   */
  std::shared_ptr<const CascadeSummingCalc> m_cascadeCalc;

  /** The volumetric-efficiency method actually used for this problem, after resolving
   `m_options.volumetric_eff_method` (Auto -> a concrete method) against the DRF.  Never
   #VolumetricEffMethod::Auto once #create has run.  Copied into every
   #DistributedSrcCalcT built for the integration. */
  ShieldingSourceFitCalc::VolumetricEffMethod m_resolvedVolEffMethod
                                    = ShieldingSourceFitCalc::VolumetricEffMethod::FlatDisk;

  /** The resolved absolute-FEP-efficiency response for volumetric sources (null for FlatDisk); built
   once by #resolveVolumetricEffMethod and shared (const) with the fit worker threads. */
  std::shared_ptr<const ceelo::DetectorResponse> m_volEffResponse;

  /** Human-readable note describing how the volumetric-efficiency method was resolved (e.g.
   "Auto -> MC transfer", or a fallback reason); empty when the request was used as-is.
   Exposed by #volumetricEffResolveNote, and carried into the reports via
   `ModelFitResults::volumetric_eff_note`. */
  std::string m_volEffResolveNote;

  /** See #volumetricEffResolveError - set only when an explicitly requested near-field method could
   not be honored. */
  std::string m_volEffResolveError;

  /** Resolves #m_options.volumetric_eff_method against #m_detector into #m_resolvedVolEffMethod +
   #m_volEffResponse (+ #m_volEffResolveNote / #m_volEffResolveError).  Builds the EFFTRAN transfer
   response when needed - so must run once, single-threaded, at the end of #create (the EFFTRAN build
   is not safe to construct concurrently, though the const evaluation is).

   NOT gated on there being a volumetric source: an off-axis or near-field point source needs the
   same correction.  No-op (leaving flat-disk) for a missing/invalid/fixed-geometry DRF, and whenever
   #resolveVolumetricEffMethodForDrf resolves to FlatDisk - which is what keeps the cost off the
   common far-field case.  The EFFTRAN anchor comes from `DetectorPeakResponse::geometry()`, which a
   DRF can carry without ever having been Monte-Carlo characterized, so a measured curve plus a known
   detector is enough. */
  void resolveVolumetricEffMethod();

public:
  /** Per-material attenuation-coefficient functions (energy -> mu*chord) along
   the center->detector ray at the given parameters, plus the air path length
   and the stacks along-ray (effective-AN, areal-density) for the GADRAS
   scatter augmentation.  Replicates the point-source propagation setup of
   #expected_peak_counts_imp / #energy_chi_contributions without touching any
   counts map - used to build the cascade-correction efficiency functors.
   Empty (no attenuation, zero air) for fixed-geometry DRFs.
   Defined in GammaInteractionCalc_imp.hpp.
   */
  template <typename T>
  struct PointSrcAttenContext
  {
    std::vector<std::function<T(float)>> att_fcns;
    T air_dist = T(0.0);        //length units
    T total_ad_gcm2 = T(0.0);   //along-ray areal density, g/cm2
    T eff_an = T(0.0);          //AD-weighted effective atomic number
  };

  template <typename T>
  PointSrcAttenContext<T> pointSourceAttenContext( const std::vector<T> &params ) const;

  /** Multiplies each entry of a single nuclides cluster map (keyed exactly by
   the #observedPeakEnergyWidths energies) by that nuclides cascade-summing
   net correction at the fit geometry.  No-op when #m_cascadeCalc is null.
   When `log_info` is provided (the logging double path), the per-peak
   cascade fields are also filled.
   Defined in GammaInteractionCalc_imp.hpp.
   */
  template <typename T>
  void applyCascadeToClusterMap( std::map<double,T> &cluster_map,
                                 const SandiaDecay::Nuclide *nuclide,
                                 const T &age,
                                 const PointSrcAttenContext<T> &atten_ctx,
                                 std::vector<GammaInteractionCalc::PeakDetail> *log_info ) const;

protected:
  
  /** The real-time of the measurement; only used if decay during measurement is being accounted for. */
  double m_realTime;
  
  std::vector<ShieldingSourceFitCalc::SourceFitDef> m_initialSrcDefinitions;
  
  //A cache of nuclide mixtures to
  mutable NucMixtureCache m_mixtureCache;
  static const size_t sm_maxMixtureCacheSize = 10000;
};//class ShieldingSourceChi2Fcn

}//namespace GammaInteractionCalc

#endif  //GammaInteractionCalc_h






