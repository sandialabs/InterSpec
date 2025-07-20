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
#include <tuple>
#include <array>
#include <atomic>
#include <vector>
#include <utility>

#include <boost/asio/deadline_timer.hpp>

#include <Wt/WColor>

#include "Minuit2/FCNBase.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/ShieldingSourceFitCalc.h"

class PeakDef;
struct Material;
class DetectorPeakResponse;

namespace SpecUtils
{
  class Measurement;
}

namespace GammaInteractionCalc
{

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


/** Integrand for Cuba library; just calls #DistributedSrcCalc::eval_spherical .
 @param userdata must be pointer to a DistributedSrcCalc object.
 */
int DistributedSrcCalc_integrand_spherical( const int *ndim, const double xx[],
                         const int *ncomp, double ff[], void *userdata );

/** Integrand for Cuba library; just calls #DistributedSrcCalc::eval_cylinder.
 @param userdata must be pointer to a DistributedSrcCalc object.
 */
int DistributedSrcCalc_integrand_cylindrical( const int *ndim, const double xx[],
                                           const int *ncomp, double ff[], void *userdata );

/** Integrand for Cuba library; just calls #DistributedSrcCalc::eval_single_cyl_end_on.
 @param userdata must be pointer to a DistributedSrcCalc object.
 
 Requires ((DistributedSrcCalc *)userdata)->m_dimensionsTransLenAndType.size() == 1 (checked by assert on debug builds)
 */
int DistributedSrcCalc_integrand_single_cyl_end_on( const int *ndim, const double xx[],
                                             const int *ncomp, double ff[], void *userdata );

/** Integrand for Cuba library; just calls #DistributedSrcCalc::eval_rect.
 @param userdata must be pointer to a DistributedSrcCalc object.
 */
int DistributedSrcCalc_integrand_rectangular( const int *ndim, const double xx[],
                                                   const int *ncomp, double ff[], void *userdata );


//point_to_line_dist(...): calculates the distance from 'point' to the line
//  passing through 'l0' and 'l1'
double point_to_line_dist( const double point[3],
                           const double l0[3], const double l1[3] );

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

struct DistributedSrcCalc
{
  //Right now this struct assumes sources are solid, in terms of the attenuation
  //  calculation
  DistributedSrcCalc();
  
  /** Evaluates the probability of gamma originating at specified coordinates making it to the detector.
   
   @param [in] xx Array of coordinates for the location to evaluate.  The coordinates are specified on interval [0,1] (which this function maps to physical location)
   @param [in] ndimptr A pointer that specifies how many dimensions is being integrated over.  The integer value must either be 2, or three.  If 2 then the corrdinates coorespond to {radius, theta}, and if 3 then {radius, theta, phi}, where the input coordinate values of [0,1] are mapped to go from [0 to sphere_rad], [0, pi] and optionally [0 to 2pi].
   @par [out] ff Pointer to where the result is placed.
   @par ncompptr unused.
   */
  void eval_spherical( const double xx[], const int *ndimptr,
                       double ff[], const int *ncompptr ) const;
  
  void eval_single_cyl_end_on( const double xx[], const int *ndimptr,
                      double ff[], const int *ncompptr ) const noexcept;
  
  void eval_cylinder( const double xx[], const int *ndimptr,
                       double ff[], const int *ncompptr ) const noexcept;
  
  void eval_rect( const double xx[], const int *ndimptr,
                        double ff[], const int *ncompptr ) const noexcept;

  GeometryType m_geometry;
  
  size_t m_materialIndex;
  double m_detectorRadius;
  double m_observationDist;
  
  /** Whether to account for attenuation in air.  If you want this, you must also set m_airTransLenCoef to the appropriate value. */
  bool m_attenuateForAir;
  
  /** The length attenuation factor for air at the energy of this object.
   
   If you want attenuation in air you must set this value when you set m_attenuateForAir to true.
   
   To get the attenuation from air using this factor you would do: exp( -m_airTransLenCoef * air_dist );
   */
  double m_airTransLenCoef;
  
  /** Wether the #m_srcVolumetricActivity should be interpreted as a surface contamination, per unit area, divided by relaxation length (L),
   with an exponential distribution (e.g. exp(-r/L) ) in the depth from surface.  E.g., for InSitu soil contamination measurements.
   If false, #m_srcVolumetricActivity is interpreted as activity per cubic-area (as its name suggests).
   */
  bool m_isInSituExponential;
  
  /** The relaxation length for in-situ exponential distribution.
   
   E.g., where ~63% of the contamination is within this distance of the surface.
   
   Not used when #m_isInSituExponential is false.
   */
  double m_inSituRelaxationLength;
  
  enum class ShellType
  {
    Material, Generic
  };
  
  /** The dimensions of this shielding, and also the transmission length coefficient for it, and if it is a material shielding or a
   generic shielding.
   
   Note that the dimensions for each layer of shielding are the outer dimensions, and not the thicknesses.
   */
  std::vector<std::tuple<std::array<double,3>,double,ShellType> > m_dimensionsTransLenAndType;

  /** The photons (i.e., sum of activities times br) per volume of the shielding.
   Is not used during integration; used to multiple integral by to get number of expected peak counts.
   */
  double m_srcVolumetricActivity;
  
  /** The energy of gamma being integrated over.
   
   Doesnt look to be used in the integration, but used for debug writing out.
   */
  double m_energy;
  
  /** The nuclide responsible for gamma being integrated over.
   
   Not used during integration - only for debug writing out; may be nullptr.
   */
  const SandiaDecay::Nuclide *m_nuclide;
  
  
  /** TODO: Setting the integral value as part of the DistributedSrcCalc is poor form - need to fix */
  double integral;
};//struct DistributedSrcCalc
  

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


  /** A struct to hold information about the material shieldings are made out of, and their respective self-attenuating and trace sources.
   
   Thicknesses and activity levels are specified by the fitting parameters, so not tracked in this struct.
   
   TODO: We could *almost* switch to using ShieldingSourceFitCalc::ShieldingInfo instead of ShieldLayerInfo,
         but this would remove our logic to use `nucs_to_fit_mass_fraction_for`, since ShieldingInfo fits all mass-fractions, or none
         (but since we only always do all or nothing, we arent actually losing much) - so we could upgrade ShieldingInfo to return mass-fractions
         to fit, and just always have it return all self-atten sources, if fitting for them is selected.
   */
  struct ShieldLayerInfo_remove_mea
  {
    /** The material this struct holds info for.
     If nullptr it represents generic shielding (e.g., AN/AD), and neither #self_atten_sources or #trace_sources may have entries.
     */
    std::shared_ptr<const Material> material;
    
    /** The nuclides which act as self-attenuating sources within the material, and their mass-fraction of the material.
     The Material object must contain these nuclides as components, as well as at least one peak being used have this
     nuclide as its assigned nuclide.  Note that some of these nuclides may also be in `nucs_to_fit_mass_fraction_for`
     (via setNuclidesToFitMassFractionFor(...)).
     */
    std::vector<std::pair<const SandiaDecay::Nuclide *,double>> self_atten_sources;
    
    /** Nuclides to fit mass fractions for; the sum of the mass fraction of all these nuclides will remain constant in the fit.
     
     Currently, our logic has it so, no matter what, all self-attenuating nuclides will be
     fit for their fractions, or none.  In principal we could allow something different,
     but this hasn't been tested.
     
     Nuclides will stored sorted alphabetically
     */
    std::vector<const SandiaDecay::Nuclide *> nucs_to_fit_mass_fraction_for;
    
    /** The trace sources info within this shielding.
     
     The 'double' argument of the tuple is the relaxation parameter, iff TraceActivityType==TraceActivityType::ExponentialDistribution.
     
     At least one peak being used have this nuclide as its assigned nuclide.
     */
    std::vector<std::tuple<const SandiaDecay::Nuclide *,TraceActivityType,double>> trace_sources;
  };//struct ShieldLayerInfo

  static std::pair<std::shared_ptr<ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> create(
                                     const double distance,
                                     const GammaInteractionCalc::GeometryType geometry,
                                     const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                     const std::vector<ShieldingSourceFitCalc::SourceFitDef> &src_definitions,
                                     std::shared_ptr<const DetectorPeakResponse> detector,
                                     std::shared_ptr<const SpecUtils::Measurement> foreground,
                                     std::shared_ptr<const SpecUtils::Measurement> background,
                                     std::deque<std::shared_ptr<const PeakDef>> foreground_peaks,
                                     std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> background_peaks,
                                     const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options );
  
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
  
  std::vector<const Material *> materialsFittingMassFracsFor() const;
  
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
  //  background peaks for the background; similarly for the live time.  This
  //  function will scale peak areas/uncertainties for the live time
  void setBackgroundPeaks( const std::vector<PeakDef> &peaks, double liveTime );
  
  
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

  
  /** For interface compatibility; calls directly to #DoEval */
  virtual double operator()( const std::vector<double> &x ) const;
  
  
  /** Gives the chi2 contributions for each peak
   
   @param params The parameters currently defining the model.
   @param error_params The errors on the parameters - only used if `log_info` is non-null
   @param mixturecache Used to speed up multiple calls to this function, and may be empty at first.
   @param info If non-null, will be filled with computation notes.
   @param log_info If non-null, filled out with details about the computation.
   */
  typedef std::map< const SandiaDecay::Nuclide *, SandiaDecay::NuclideMixture> NucMixtureCache;
  std::vector<PeakResultPlotInfo> energy_chi_contributions(
                                  const std::vector<double> &params,
                                  const std::vector<double> &error_params,
                                  NucMixtureCache &mixturecache,
                                  std::vector<std::string> *info = nullptr,
                                  std::vector<PeakDetail> *log_info = nullptr ) const;

  void log_shield_info( const std::vector<double> &params,
                        const std::vector<double> &error_params,
                        const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &fit_src_info,
                        std::vector<ShieldingDetails> &info ) const;
  
  void log_source_info( const std::vector<double> &params,
                        const std::vector<double> &error_params,
                        const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &fit_src_info,
                        std::vector<SourceDetails> &info ) const;
  
  
  ShieldingSourceChi2Fcn&	operator=( const ShieldingSourceChi2Fcn & );
  virtual double Up() const;

  size_t numExpectedFitParameters() const;

  //nuclide(...) throws std::runtime_error if an invalid number is passed in
  const SandiaDecay::Nuclide *nuclide( const size_t nuclide_index ) const;

  //activity(...) and age(...) will throw runtime_exception if params is wrong
  //  size or invalid nuclide passed in
  double activity( const SandiaDecay::Nuclide *nuclide,
                   const std::vector<double> &params ) const;
  double activityUncertainty( const SandiaDecay::Nuclide *nuclide,
                  const std::vector<double> &params,
                  const std::vector<double> &errors ) const;
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
  //  NULL if a generic material.
  const Material *material( const size_t materialNum ) const;
  
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
  const std::vector<PeakDef> &backgroundPeaks() const;

  double distance() const;
  
  const std::shared_ptr<const DetectorPeakResponse> &detector() const;
  
  const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options() const;
  
  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &initialShieldings() const;
  
  static void selfShieldingIntegration( DistributedSrcCalc &calculator );

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
                  std::vector<std::string> *info,
                  std::vector<GammaInteractionCalc::PeakDetail> *log_info
              );


  //returns the chi computed from the expected verses observed counts; one
  //  chi2 for each peak energy.  Each returned entry is {energy,chi,scale,PeakColor,ScaleUncert},
  //  where scale is observed/expected
  static std::vector<PeakResultPlotInfo> expected_observed_chis(
                              const std::vector<PeakDef> &peaks,
                              const std::vector<PeakDef> &backgroundPeaks,
                              const std::map<double,double> &energy_count_map,
                              std::vector<std::string> *info = 0,
                              std::vector<GammaInteractionCalc::PeakDetail> *log_info = nullptr );
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
  double m_liveTime;

  std::vector<PeakDef> m_peaks;
  std::vector<PeakDef> m_backgroundPeaks;
  std::shared_ptr<const DetectorPeakResponse> m_detector;
  
  const std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_initial_shieldings;
  
  std::vector<const SandiaDecay::Nuclide *> m_nuclides; //sorted alphabetically and unique
  
  const GeometryType m_geometry;
  
  ShieldingSourceFitCalc::ShieldingSourceFitOptions m_options;
  
  /** The real-time of the measurement; only used if decay during measurement is being accounted for. */
  double m_realTime;
  
  std::vector<ShieldingSourceFitCalc::SourceFitDef> m_initialSrcDefinitions;
  
  //A cache of nuclide mixtures to
  mutable NucMixtureCache m_mixtureCache;
  static const size_t sm_maxMixtureCacheSize = 10000;
};//class ShieldingSourceChi2Fcn

}//namespace GammaInteractionCalc

#endif  //GammaInteractionCalc_h






