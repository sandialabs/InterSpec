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

#include <set>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <cctype>
#include <istream>
#include <fstream>
#include <sstream>
#include <utility>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/asio/deadline_timer.hpp>

#include <Wt/WServer>
#include <Wt/WIOService>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/DateTime.h"
#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/DrfSelect.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/Integrate.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ShieldingSourceFitCalc.h"

using namespace std;

#if( DEBUG_RAYTRACE_CALCS )
#include <mutex>
static std::recursive_mutex s_stdout_raytrace_mutex;

#  ifdef NDEBUG
static_assert( 0, "Disable DEBUG_RAYTRACE_CALCS for release builds" );
#  else
#    ifndef _MSC_VER
#      warning "DEBUG_RAYTRACE_CALCS is enabled"
#    else
#      pragma message("Warning: DEBUG_RAYTRACE_CALCS is enabled")
#    endif //_MSC_VER / else
#  endif //NDEBUG / else
#endif //DEBUG_RAYTRACE_CALCS

namespace GammaInteractionCalc
{

const char *to_str( const TraceActivityType type )
{
  switch( type )
  {
    case TraceActivityType::TotalActivity:
      return "TotalActivity";
      
    case TraceActivityType::ActivityPerCm3:
      return "ActivityPerCm3";
      
    case TraceActivityType::ActivityPerGram:
      return "ActivityPerGram";
      
    case TraceActivityType::ExponentialDistribution:
      return "ExponentialDistribution";
      
    case TraceActivityType::NumTraceActivityType:
      return "NumTraceActivityType";
  }//switch( type )
  
  return "InvalidTraceActivityType";
}//to_str( TraceActivityType )


const char *to_str( const GeometryType type )
{
  switch( type )
  {
    case GeometryType::Spherical:       return "Spherical";
    case GeometryType::CylinderEndOn:   return "CylinderEndOn";
    case GeometryType::CylinderSideOn:  return "CylinderSideOn";
    case GeometryType::Rectangular:     return "Rectangular";
    case GeometryType::NumGeometryType: return "NumGeometryType";
  }//switch( type )
  
  return "InvalidGeometryType";
}//to_str( GeometryType )



const double ShieldingSourceChi2Fcn::sm_activityUnits = SandiaDecay::MBq;

  
//Returned in units of 1.0/[Length], so that
//  exp( -transmition_length_coefficient(...) * thickness)
//  gives you the probability a gamma of given energy will go through the
//  material of given thickness.
double transmition_length_coefficient( const Material *material, float energy )
{
  double mu = 0.0;
  
  for( const Material::ElementFractionPair &p : material->elements )
  {
    const SandiaDecay::Element *el = p.first;
    const int atomicNumber = el->atomicNumber;
    if( atomicNumber > MassAttenuation::sm_max_xs_atomic_number )
      throw std::runtime_error( "transmition_length_coefficient(...): invalid "
                                "atomic number" );
    const double density = p.second * material->density;
    const double xs_per_mass = MassAttenuation::massAttenuationCoefficientElement( atomicNumber, energy );
    
    mu += (density*xs_per_mass);
  }//for( const Material::ElementFractionPair &p : material->elements )

  for( const Material::NuclideFractionPair &p : material->nuclides )
  {
    const int atomicNumber = p.first->atomicNumber;

    if( atomicNumber > MassAttenuation::sm_max_xs_atomic_number )
      throw std::runtime_error( "transmition_length_coefficient(...): invalid "
                                "atomic number" );
    double density = p.second * material->density;
    const double xs_per_mass = MassAttenuation::massAttenuationCoefficientElement( atomicNumber, energy );
    
    mu += (density*xs_per_mass);
  }//for( const Material::ElementFractionPair &p : material->elements )

  return mu;
}//double transmition_length_coefficient(...)


double transmition_coefficient_material( const Material *material, float energy,
                                float length )
{
  return length * transmition_length_coefficient( material, energy );
}


double transmission_length_coefficient_air( float energy )
{
  double mu = 0.0;
  
  // Instead of using the components, we could use the mass-weighted atomic number
  //const float air_an = 7.3737f;  //Gadras uses 7.2
  //const float air_density = static_cast<float>( 0.00129 * PhysicalUnits::g / PhysicalUnits::cm3 );
  //const double mu = MassAttenuation::massAttenuationCoefficientFracAN( air_an, energy );
  //return mu * air_density;
  
  //Air (taken from definition used in InterSpec v1.0.8 20211017):
  //  AN=7, Density=6.08133411407 (0.000974337052716 g/cm3)
  //  AN=8, density=1.86634886265 (0.000299022026427 g/cm3)
  //  AN=18, density=0.103864975274 (1.66410021206e-05 g/cm3)
  const pair<int,double> air_def[] = {{7,6.08133411407}, {8,1.86634886265}, {18,0.103864975274}};
  
  for( const auto &i : air_def )
  {
    const int atomicNumber = i.first;
    const double density = i.second;
    const double xs_per_mass = MassAttenuation::massAttenuationCoefficientElement( atomicNumber, energy );
    
    mu += (density*xs_per_mass);
  }//for( const Material::ElementFractionPair &p : material->elements )
  
  return mu;
}//transmission_length_coefficient_air(...)


double transmission_coefficient_air( float energy, float length )
{
  return length * transmission_length_coefficient_air( energy );
}


//Returned in units of [Length]^2/[mass], so that
//  exp( -mass_attenuation_coef * areal_density )
//  gives you the probability a gamma of given energy will go through the
//  material with given atomic number and areal_density.
//  The quantity retuned by this function is commonly labeled μ
double mass_attenuation_coef( float atomic_number, float energy )
{
  const double xs_per_mass = MassAttenuation::massAttenuationCoefficientFracAN( atomic_number, energy );
  
  return xs_per_mass;
}

double transmition_coefficient_generic( float atomic_number, float areal_density,
                                float energy )
{
  return areal_density * mass_attenuation_coef( atomic_number, energy );
}
  
  
vector<SandiaDecay::EnergyRatePair> decay_during_meas_corrected_gammas(
                                                      const SandiaDecay::NuclideMixture &mixture,
                                                      const double age,
                                                      const double measDuration )
{
  // TODO: SandiaDecay had `SandiaDecay::NuclideMixture::decayPhotonsInInterval(...)` added 20241124
  //       so we could use that code, however, not doing that now because it is single threaded.
  // This section of code takes up a good amount of CPU time, and is really begging for optimizations
  //
  // TODO: if nuclide decays to stable children, then just use standard formula to correct for decay (see below developer check for this)
  //
  // We want to get the activity at the beginning of the measurement, but we want
  //  to account for decay (or build up!) throughout the measurement; rather than
  //  being intelligent about it, we'll just average throughout the dwell time
  //
  // TODO: figure out how many time slices are needed for reasonable accuracy.
  //       A first quick go at this, for In110 (hl=4.9h), over a 2.8 hour measurement,
  //       gave corrections of
  //            For 250: 0.826922  (which matches analytical answer of 0.826922)
  //            For 50:  0.82692
  //            For 25:  0.826913
  //            For 10:  0.826869
  //
  //  A spot check using Mn56 (hl=9283.8s), and a meas duration of 86423.86s, with 50 timeslices
  //    gives a correction factor of 0.15473519892011101, vs analytical answer of 0.1547364350135815
  
  const size_t characteristic_time_slices = 50; // Maybe at least ~5 sig figs
  
  if( mixture.numInitialNuclides() != 1 )
    throw runtime_error( "ShieldingSourceChi2Fcn::decayCorrectedGammas():"
                        " passed in mixture must have exactly one parent nuclide" );
  const SandiaDecay::Nuclide *nuclide = mixture.initialNuclide(0);
  assert( nuclide );
  if( !nuclide )
    throw std::logic_error( "decayCorrectedGammas: nullptr nuc" );
  
  
  // TODO: the parent half-life may not be the relevant one; should check down the chain
  const double halflife = std::max( nuclide->halfLife, 0.01*PhysicalUnits::second );
  const double characteristicTime = std::min( measDuration, halflife );
  const double dt = characteristicTime / characteristic_time_slices;
  
  const int num_time_slices = std::max( 50, static_cast<int>( std::ceil( measDuration / dt ) ) );
  const int nthread = SpecUtilsAsync::num_logical_cpu_cores();

  
  vector<SandiaDecay::EnergyRatePair> gammas
                             = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );
  
  assert( std::is_sorted(begin(gammas), end(gammas),
         []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
    return lhs.energy < rhs.energy;
  }) );
  
  vector<vector<SandiaDecay::EnergyRatePair>> energy_rate_pairs( nthread,
                                                                vector<SandiaDecay::EnergyRatePair>(gammas.size(), {0.0,0.0}) );
  
  
  SpecUtilsAsync::ThreadPool pool;
  
  for( int threadnum = 0; threadnum < nthread; ++threadnum )
  {
    pool.post( [threadnum, nthread, &energy_rate_pairs, num_time_slices, measDuration, age, &mixture](){
      
      vector<SandiaDecay::EnergyRatePair> &result = energy_rate_pairs[threadnum];
      
      // TODO: right now using simple midpoint integration; should boost::math::quadrature::trapezoidal, or boost::math::quadrature::gauss::integrate, or something, but then evaluating things becomes a little tricky (I think we need to decay to the set number of points beforehand, and then make a function to retireve these inside the boost integration routine
      for( int timeslice = threadnum; timeslice < num_time_slices; timeslice += nthread )
      {
        const double this_age = age + measDuration*(2.0*timeslice + 1.0)/(2.0*num_time_slices);
        vector<SandiaDecay::EnergyRatePair> these_gammas
        = mixture.photons( this_age, SandiaDecay::NuclideMixture::OrderByEnergy );
        
        assert( these_gammas.size() == result.size() );
        if( these_gammas.size() != result.size() )
          throw std::logic_error( "gamma result size doesnt match expected." );
        
        for( size_t i = 0; i < these_gammas.size(); ++i )
          result[i].numPerSecond += these_gammas[i].numPerSecond;
      }//for( loop over timeslices )
      
      for( size_t i = 0; i < result.size(); ++i )
        result[i].numPerSecond /= num_time_slices;
    } );
  }//for( int threadnum = 0; threadnum < nthread; ++threadnum )
  
  pool.join();
  
  vector<SandiaDecay::EnergyRatePair> corrected_gammas = gammas;
  for( size_t i = 0; i < corrected_gammas.size(); ++i )
  {
    corrected_gammas[i].numPerSecond = 0.0;
    for( size_t threadnum = 0; threadnum < nthread; ++threadnum )
    {
      assert( corrected_gammas.size() == energy_rate_pairs[threadnum].size() );
      corrected_gammas[i].numPerSecond += energy_rate_pairs[threadnum][i].numPerSecond;
    }
  }//for( size_t i = 0; i < corrected_gammas.size(); ++i )
  
  //cout << "For " << nuclide->symbol << " the time corrections are:" << endl;
  //for( size_t i = 0; i < corrected_gammas.size(); ++i )
  //  cout << std::setw(15) << corrected_gammas[i].energy << ": "
  //       << corrected_gammas[i].numPerSecond/gammas[i].numPerSecond
  //       << " (" << corrected_gammas[i].numPerSecond << " vs orig "
  //       << gammas[i].numPerSecond << ")" << endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
  {// begin check against SandiaDecay calc
    const vector<SandiaDecay::EnergyCountPair> sdanswer
            = mixture.decayPhotonsInInterval( age, measDuration,
                                             SandiaDecay::NuclideMixture::OrderByEnergy,
                                            characteristic_time_slices );
    assert( corrected_gammas.size() == sdanswer.size() );
    for( size_t i = 0; i < corrected_gammas.size(); ++i )
    {
      const double sdRate = sdanswer[i].count / measDuration;
      const double localRate = corrected_gammas[i].numPerSecond;
      assert( sdanswer[i].energy == corrected_gammas[i].energy );
      assert( (fabs(sdRate - localRate) < 1.0E-10*std::max(sdRate, localRate))
             || (fabs(sdRate - localRate) < 1.0E-14) );
    }
  }// end check against SandiaDecay calc
  
  if( nuclide->decaysToStableChildren() )
  {
    const double lambda = nuclide->decayConstant();
    const double corr_factor = (1.0 - exp(-1.0*lambda*measDuration)) / (lambda * measDuration);
    
    //cout << "corr_factor=" << corr_factor << endl;
    assert( gammas.size() == corrected_gammas.size() );
    for( size_t i = 0; i < corrected_gammas.size(); ++i )
    {
      const double numerical_answer = corrected_gammas[i].numPerSecond;
      const double uncorrected_answer = gammas[i].numPerSecond;
      if( (uncorrected_answer > DBL_EPSILON) && (uncorrected_answer > DBL_EPSILON) )
      {
        const double numerical_corr_factor = numerical_answer / uncorrected_answer;
        if( fabs(numerical_corr_factor - corr_factor) > 0.001 )
        {
          string msg = "Found decay correction value of "
          + std::to_string(numerical_corr_factor) + ", when a true value of "
          + std::to_string(corr_factor) + " was expected for " + nuclide->symbol
          + " with half life "
          + PhysicalUnits::printToBestTimeUnits(nuclide->halfLife,6)
          + " and a measurement time "
          + PhysicalUnits::printToBestTimeUnits(measDuration,6)
          + " (at energy " + std::to_string(corrected_gammas[i].energy) + " keV)";
          log_developer_error( __func__, msg.c_str() );
        }//if( error is larger than expected )
      }//if( not a zero BR gamma )
    }//for( loop over gammas )
  }//if( we can use standard formula to correct )
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  //cout << endl << endl << endl;
  
  assert( std::is_sorted(begin(corrected_gammas), end(corrected_gammas),
         []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
    return lhs.energy < rhs.energy;
  }) );
  
  return corrected_gammas;
}//decay_during_meas_corrected_gammas(...)
  


double exit_point_of_sphere_z( const double source_point[3],
                               double exit_point[3],
                               double sphere_rad,
                               double observation_dist,
                               bool postiveSolution )
{
  /*
   *Makes x0,y0,z0 be the point of intersection of a sphere or radius
   *'sphere_rad' centered at the origin, with the line pointing from <x0,y0,z0>
   *towards <0,0,observation_dist>
   *Note: CAN be used with a 2-dimensional integral which holds phi constant
   */

  using namespace std;

  const double a = source_point[0];
  const double b = source_point[1];
  const double c = source_point[2];
  const double &S = sphere_rad;
  const double &R = observation_dist;

  const double r = sqrt( a*a + b*b + c*c );
  if( observation_dist < sphere_rad )
  {
    cerr << "observation_dist=" << observation_dist
         << ", sphere_rad=" << sphere_rad << endl;
    throw runtime_error( "exit_point_of_sphere_z(...): obs_dist < sphere_rad" );
  }//if( observation_dist < sphere_rad )

  if( r > sphere_rad )
  {
    if( ((r-sphere_rad)/sphere_rad) < 0.0001 )
    {
      exit_point[0] = source_point[0];
      exit_point[1] = source_point[1];
      exit_point[2] = source_point[2];
      return 0.0;
    }//if( this is just a rounding error )
    
    cerr << "sphere_rad=" << sphere_rad << ", r=" << r << endl;
    throw runtime_error( "exit_point_of_sphere_z(...): r > sphere_rad" );
  }//if( r > sphere_rad )

  //TODO: factor the math below to save CPU time
  if( postiveSolution )
  {
    exit_point[0] =-(a*sqrt((R*R-2*c*R+c*c+b*b+a*a)*S*S+(-b*b-a*a)*R*R)-a*R*R+a*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[1] =(-b*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*R*R-b*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[2] =(R*(sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*b+a*a)-c*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R))/(R*R-2*c*R+c*c+b*b+a*a);
    
    //  const double factor_1 = R*R-2*c*R+c*c+b*b+a*a;
    //  const double factor_2 = R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R;
    //  const double x_pos_new =-(a*sqrt((factor_1)*S*S+(-b*b-a*a)*R*R)-a*R*R+a*c*R)/factor_1;
    //  const double y_pos_new =(-b*sqrt(factor_2)+b*R*R-b*c*R)/factor_1;
    //  const double z_pos_new =(R*(sqrt(factor_2)+b*b+a*a)-c*sqrt(factor_2))/(R*R-2*c*R+c*c+b*b+a*a);
    //  assert( x_pos_new == exit_point[0] );
    //  assert( y_pos_new == exit_point[1] );
    //  assert( z_pos_new == exit_point[2] );
  }else
  {
    exit_point[0] = (a*sqrt((R*R-2*c*R+c*c+b*b+a*a)*S*S+(-b*b-a*a)*R*R)+a*R*R-a*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[1] = (b*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*R*R-b*c*R)/(R*R-2*c*R+c*c+b*b+a*a);
    exit_point[2] = (c*sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+R*(-sqrt(R*R*S*S-2*c*R*S*S+c*c*S*S+b*b*S*S+a*a*S*S-b*b*R*R-a*a*R*R)+b*b+a*a))/(R*R-2*c*R+c*c+b*b+a*a);
  }//if( postiveSolution ) / else

  const double dx = a - exit_point[0];
  const double dy = b - exit_point[1];
  const double dz = c - exit_point[2];

  return sqrt( dx*dx + dy*dy + dz*dz );
}//exit_point_of_sphere_z


double distance( const double a[3], const double b[3] )
{
  const double dx = a[0]-b[0];
  const double dy = a[1]-b[1];
  const double dz = a[2]-b[2];

  return sqrt( dx*dx + dy*dy + dz*dz );
}//double distance( const double a[3], b[3] )
  

int DistributedSrcCalc_integrand_spherical( const int *ndim, const double xx[],
                                           const int *ncomp, double ff[], void *userdata )
{
  const DistributedSrcCalc * const objToIntegrate = (DistributedSrcCalc *)userdata;
  
  assert( objToIntegrate );
  assert( objToIntegrate->m_geometry == GeometryType::Spherical );
  
  objToIntegrate->eval_spherical( xx, ndim, ff, ncomp );
  
  return 0;
}//DistributedSrcCalc_integrand_spherical(...)


int DistributedSrcCalc_integrand_cylindrical( const int *ndim, const double xx[],
                                             const int *ncomp, double ff[], void *userdata )
{
  const DistributedSrcCalc * const objToIntegrate = (DistributedSrcCalc *)userdata;
  
  assert( objToIntegrate );
  assert( (objToIntegrate->m_geometry == GeometryType::CylinderEndOn)
          || (objToIntegrate->m_geometry == GeometryType::CylinderSideOn) );
  
  objToIntegrate->eval_cylinder( xx, ndim, ff, ncomp );
  
  return 0;
}//DistributedSrcCalc_integrand_cylindrical(...)


int DistributedSrcCalc_integrand_single_cyl_end_on( const int *ndim, const double xx[],
                                                   const int *ncomp, double ff[], void *userdata )
{
  const DistributedSrcCalc * const objToIntegrate = (DistributedSrcCalc *)userdata;
  
  assert( objToIntegrate );
  assert( objToIntegrate->m_geometry == GeometryType::CylinderEndOn );
  assert( objToIntegrate->m_dimensionsTransLenAndType.size() == 1 );
  
  objToIntegrate->eval_single_cyl_end_on( xx, ndim, ff, ncomp );
  
  return 0;
}//DistributedSrcCalc_integrand_single_cyl_end_on(...)


int DistributedSrcCalc_integrand_rectangular( const int *ndim, const double xx[],
                                          const int *ncomp, double ff[], void *userdata )
{
  const DistributedSrcCalc * const objToIntegrate = (DistributedSrcCalc *)userdata;
  
  assert( objToIntegrate );
  assert( objToIntegrate->m_geometry == GeometryType::Rectangular );
  
  objToIntegrate->eval_rect( xx, ndim, ff, ncomp );
  
  return 0;
}//DistributedSrcCalc_integrand_rectangular(...)




DistributedSrcCalc::DistributedSrcCalc()
{
  m_geometry = GeometryType::NumGeometryType;
  m_materialIndex = 0;
  m_detectorRadius = -1.0;
  m_observationDist = -1.0;
  m_attenuateForAir = false;
  m_airTransLenCoef = 0.0;
  m_isInSituExponential = false;
  m_inSituRelaxationLength = 0.0;
  m_nuclide = NULL;
}//DistributedSrcCalc()

  
  
double point_to_line_dist( const double point[3],
                           const double p0[3], const double p1[3] )
{
  double distsqrd = 0.0;
  double v[3], w[3];
  for( int i = 0; i < 3; ++i )
  {
    v[i] = p1[i] - p0[i];
    w[i] = point[i] - p0[i];
  }
  
  const double d1 = w[0]*v[0] + w[1]*v[1] + w[2]*v[2];
  const double d2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  const double frac = d1 / d2;
  
  for( int i = 0; i < 3; ++i )
  {
    const double t = point[i] - (p0[i] + frac*v[i]);
    distsqrd += t*t;
  }
  
  return sqrt( distsqrd );
}//double point_to_line_dist(...)
  
  

void DistributedSrcCalc::eval_spherical( const double xx[], const int *ndimptr,
                        double ff[], const int *ncompptr ) const
{
  assert( m_geometry == GeometryType::Spherical );
  assert( m_materialIndex < m_dimensionsTransLenAndType.size() );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_materialIndex]) == ShellType::Material );
  
  const double pi = PhysicalUnits::pi;
  
  const int ndim = (ndimptr ? (*ndimptr) : 3);
  assert( (ndim == 2) || (ndim == 3) );

  const double source_inner_rad = ((m_materialIndex>0)
                            ? std::get<0>(m_dimensionsTransLenAndType[m_materialIndex-1])[0]
                            : 0.0);
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex]);
  const double source_outer_rad = dimensions[0];

  //cuba goes from 0 to one, so we have to scale the variables
  //r:     goes from zero to sphere_rad
  //theta: goes from 0 to pi
  //phi:   goes from 0 to 2pi

  double max_theta = pi;
  
  const double x_r = xx[0];
  const double x_theta = xx[1];
  const double x_phi = (ndim<3 ? 0.5 : xx[2]);
  
  const double r = source_inner_rad + x_r * (source_outer_rad - source_inner_rad);
  const double theta =  x_theta * max_theta;
  
  //phi should go from zero to two pi, but we could reduce this to just 0 to
  //  pi, due to symmetry, but not doing right now (cause I'm to lazy to test it)
  const double max_phi = 2.0 * pi;
  const double phi = x_phi * max_phi;
  const double j = (source_outer_rad - source_inner_rad) * max_theta * max_phi;
  const double dV = j * r * r * sin(theta);

  double source_point[3];
  source_point[0] = r * sin( theta ) * cos( phi );
  source_point[1] = r * sin( theta ) * sin( phi );
  source_point[2] = r * cos( theta );

  // We'll the source location since 'source_point' will get modified.
  const double x = source_point[0], y = source_point[1], z = source_point[2];

  double trans = 0.0;
  
  {//begin code-block to compute distance through source
    // - this could probably be cleaned up and made more efficient
    double exit_point[3];
    const double srcRad = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex])[0];
    const double srcTransCoef = std::get<1>(m_dimensionsTransLenAndType[m_materialIndex]);
    double dist_in_src = exit_point_of_sphere_z( source_point, exit_point,
                                                 srcRad, m_observationDist );
    
    double min_rad = 0.0;
    bool needShellCompute = false;
    double inner_shell_point[3] = { 0.0 };
    
    if( m_materialIndex > 0 )
    {
      exit_point_of_sphere_z( source_point, inner_shell_point, srcRad,
                                           m_observationDist, false );
      
      //Take advantage of symmetry, move the position to half way between the
      //  positive and negative exit points, and then check if the radius of
      //  this point is both in a sub shell, and that its closer to the detector
      //  then the source point.
      //We will then use this same point to calc distances through the shells
      //  to the detectors, and then double that distance for all but the source
      //  shell (of which we will acount for dist from originating pos to first
      //  sub-shell)
      inner_shell_point[0] = 0.5*(inner_shell_point[0] + exit_point[0]);
      inner_shell_point[1] = 0.5*(inner_shell_point[1] + exit_point[1]);
      inner_shell_point[2] = 0.5*(inner_shell_point[2] + exit_point[2]);
      min_rad = sqrt( inner_shell_point[0]*inner_shell_point[0]
                      + inner_shell_point[1]*inner_shell_point[1]
                      + inner_shell_point[2]*inner_shell_point[2] );
      
      const array<double,3> &sub_dims = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex-1]);
      const double subrad = sub_dims[0];
      needShellCompute = ( (min_rad < subrad) && (inner_shell_point[2] > source_point[2]) );
      // Note that if this shell is a generic shell, and the ray *just* touches it, then we wont
      //  include the generic attenuation.
    }//if( m_materialIndex > 0 )
    
    if( needShellCompute )
    {
      double original_point[3];
      memcpy( original_point, source_point, 3*sizeof(double) );
      memcpy( source_point, inner_shell_point, 3*sizeof(double) );
      
      //Calc how far from the gammas original position it was, to the first
      //  inner shell it hit
      const array<double,3> &sub_dims = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex-1]);
      const double innerRad = sub_dims[0];
      exit_point_of_sphere_z( source_point, exit_point,
                                        innerRad, m_observationDist, false );
      double srcDist2 = 0.0;
      for( int i = 0; i < 3; ++i )
      {
        const double diff = exit_point[i] - original_point[i];
        srcDist2 += diff*diff;
      }//for( int i = 0; i < 3; ++i )
      trans += (srcTransCoef * sqrt(srcDist2));
     
      //find inner most sphere the ray passes through
      size_t start_index = 0;
      while( std::get<0>(m_dimensionsTransLenAndType[start_index])[0] < min_rad
            && start_index < m_dimensionsTransLenAndType.size() )
        ++start_index;
      
      //Some hopefully un-needed logic checks
      assert( m_materialIndex != 0 );
      assert( start_index < m_materialIndex );
      assert( start_index < m_dimensionsTransLenAndType.size() );
        
      if( start_index == m_dimensionsTransLenAndType.size() )
        throw runtime_error( "Logic error 1 in DistributedSrcCalc::eval(...)" );
      if( start_index >= m_materialIndex )
        throw runtime_error( "Logic error 2 in DistributedSrcCalc::eval(...)" );
      if( m_materialIndex == 0 )
        throw runtime_error( "Logic error 3 in DistributedSrcCalc::eval(...)" );
      
      //calc distance it travels through the inner spheres
      for( size_t index = start_index; index <= m_materialIndex; ++index )
      {
        const std::array<double,3> &dims = std::get<0>(m_dimensionsTransLenAndType[index]);
        const double transCoef = std::get<1>(m_dimensionsTransLenAndType[index]);
        const ShellType type = std::get<2>(m_dimensionsTransLenAndType[index]);
        assert( (type == ShellType::Material) || (type == ShellType::Generic) );
        
        switch( type )
        {
          case ShellType::Material:
          {
            const double shellRad = dims[0];
            double dist = exit_point_of_sphere_z( source_point, source_point,
                                                 shellRad, m_observationDist );
            if( index != m_materialIndex )
              dist = 2.0*dist;
            
            trans += (transCoef * dist);
            break;
          }//case ShellType::Material:
            
          case ShellType::Generic:
          {
            trans += transCoef;
            
            // Make sure this is a zero dimension shell
            assert( !index || (dims == std::get<0>(m_dimensionsTransLenAndType[index-1])) );
            break;
          }//case ShellType::Generic:
        }//switch( type )
        
        
      }//for( ++sphere_index; sphere_index < m_materialIndex; ++sphere_index )
    }else
    {
      memcpy( source_point, exit_point, 3*sizeof(double) );
      trans += (srcTransCoef * dist_in_src);
    }//if( line actually goes into child sphere ) / else
  }//end codeblock to compute distance through source
  
  for( size_t i = m_materialIndex+1; i < m_dimensionsTransLenAndType.size(); ++i )
  {
    const std::array<double,3> &dims = std::get<0>(m_dimensionsTransLenAndType[i]);
    const double transLenCoef = std::get<1>(m_dimensionsTransLenAndType[i]);
    const ShellType type = std::get<2>(m_dimensionsTransLenAndType[i]);
    
    const double sphereRad = dims[0];
    
    switch( type )
    {
      case ShellType::Material:
      {
        const double dist_in_sphere = exit_point_of_sphere_z( source_point,
                                                             source_point, sphereRad, m_observationDist );
        trans += (transLenCoef * dist_in_sphere);
        
        break;
      }//case ShellType::Material:
        
      case ShellType::Generic:
      {
        trans += transLenCoef;
        
        // Make sure this is a zero dimension shell
        assert( !i || (dims == std::get<0>(m_dimensionsTransLenAndType[i-1])) );
        break;
      }//case ShellType::Generic:
    }//switch( type )
  }//for( size_t i = 0; i < m_dimensionsTransLenAndType.size(); ++i )

  trans = exp( -trans );
  
  
  if( m_attenuateForAir )
  {
    // source_point is the exit point on the last of the shielding, and the detector is at
    //  [0,0,m_observationDist]
    const double dx = -source_point[0];
    const double dy = -source_point[1];
    const double dz = source_point[2] - m_observationDist;
    
    const double air_dist = sqrt( dx*dx + dy*dy + dz*dz );
    
    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )
  
  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(source_outer_rad - r) / m_inSituRelaxationLength );
  }
  
  const double z_dist = (z - m_observationDist);
  const double dist_to_det = sqrt( x*x + y*y + z_dist*z_dist );
  
  // Note: previous to 20211104 m_observationDist was used instead of dist_to_det; I believe this
  //       change is an appropriate correction, but still needs to be validated/ensured.
  trans *= DetectorPeakResponse::fractionalSolidAngle( 2.0*m_detectorRadius, dist_to_det );

  ff[0] = trans * dV;
}//eval_spherical(...)


void DistributedSrcCalc::eval_single_cyl_end_on( const double xx[], const int *ndimptr,
                                          double ff[], const int *ncompptr ) const noexcept
{
#if( DEBUG_RAYTRACE_CALCS )
  std::lock_guard<std::recursive_mutex> scoped_lock( s_stdout_raytrace_mutex );
#endif

  assert( m_materialIndex == 0 );
  assert( m_geometry == GeometryType::CylinderEndOn );
  assert( m_dimensionsTransLenAndType.size() == 1 );
  assert( std::get<2>(m_dimensionsTransLenAndType[0]) == ShellType::Material );
  
  //if( m_dimensionsTransLenAndType.size() != 1 )
  //  throw runtime_error( "eval_single_cyl_end_on only supports a single shielding" );
  
  
  // This just integrates a right circular cylinder
  const int ndim = (ndimptr ? (*ndimptr) : 2);
  assert( ndim == 2 ); // ndim==3 is also valid and this function will work for that as well
  
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex]);
  const double source_outer_rad = dimensions[0];
  const double source_half_z = dimensions[1];
  const double total_height = 2.0 * source_half_z;
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_materialIndex]);
  
  
  // Right now we are only dealing with a single cylinder; however, when we move to nesting
  //  cylinders, I think we'll integrate over the entire source cylinder volume, but just set the
  //  source term for the inner cylinders to zero, meaning I think we will just remove this
  //  'source_inner_rad' variable, but just leaving it in for them moment as a reminder to think
  //  about exactly how to do things a bit more.
  const double source_inner_rad = 0.0;
  
  
  //cuba goes from 0 to one for each dimension, so we have to scale the variables
  //  r:     goes from cylindrical inner radius, to outer radius
  //  theta: goes from 0 to 2pi
  //  z:     goes from the negative half-height to positive half-height
  
  double max_theta = 2.0 * PhysicalUnits::pi;
  
  const double x_r = xx[0];
  const double x_theta = ((ndim==3) ? xx[1] : 0.0);
  const double x_z = xx[ ((ndim==3) ? 2 : 1) ];
  
  const double r = source_inner_rad + x_r * (source_outer_rad - source_inner_rad);
  const double theta = x_theta * max_theta;
  const double z = total_height * (x_z - 0.5);
  
  const double j = (source_outer_rad - source_inner_rad) * max_theta * total_height;
  const double dV = j * r;
  
  const double eval_z_dist_to_det = m_observationDist - z;
  
  //const double eval_point[3] = { r * cos(theta), r * sin(theta), z };

  double exit_radius = 0.0;
  double trans = 0.0;
  
  {//begin code-block to compute distance through source
    // Take advantage of theta symmetry here
    const double z_dist_in_src = source_half_z - z;
    const double r_dist_in_src = r * z_dist_in_src / eval_z_dist_to_det;
    exit_radius = r - r_dist_in_src;
    
    const double dist_in_src = sqrt(z_dist_in_src*z_dist_in_src + r_dist_in_src*r_dist_in_src);
    
    trans += (trans_len_coef * dist_in_src);
  }//end codeblock to compute distance through source

  trans = exp( -trans );
  
  if( m_attenuateForAir )
  {
    const double dz = m_observationDist - source_half_z;
    
    const double air_dist = sqrt( exit_radius*exit_radius + dz*dz );
    
    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )
  
  
  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(source_half_z - z) / m_inSituRelaxationLength );
  }
  
  // Finally toss in the geometric factor (e.g., 1/r2 from where we are evaluating to to detector).
  const double eval_dist_to_det = sqrt( r*r + eval_z_dist_to_det*eval_z_dist_to_det );
  trans *= DetectorPeakResponse::fractionalSolidAngle( 2.0*m_detectorRadius, eval_dist_to_det );
  
  
  // For debug builds, also check against the more general transport
#ifndef NDEBUG
  double test_ff[1];
  eval_cylinder( xx, ndimptr, test_ff, ncompptr );
  const double this_answer = trans * dV;
  assert( fabs(test_ff[0] - this_answer) < 1.0E-6*std::max(0.001,std::max( fabs(test_ff[0]), fabs(this_answer))) );
#endif

  
  ff[0] = trans * dV;
}//void eval_single_cyl_end_on(...)


// TODO: define/use a proper simple 3-point vector struct
//struct ThreeVec : public std::array<double,3>
//{
//  ThreeVec() : std::array<double,3>( { 0.0, 0.0, 0.0 } ){}
//  double x() const { return this->operator[](0); }
//  double y() const { return this->operator[](1); }
//  double z() const { return this->operator[](2); }
//};//struct ThreeVec



double cylinder_line_intersection( const double radius, const double half_length,
                              const double source[3],
                              const double detector[3],
                              const CylExitDir direction,
                              double exit_point[3] ) noexcept
{
  // TODO: this function should be broken into two separate functions.  One to handle finding the exit
  //  point when you know the source is inside the volume.  And one to find both intersection points
  //  (if any) of external points.  This would both increase the efficiency of the function, and also
  //  make the use cleaner/easier.
  
  // TODO: need to clearly define what happens for points exactly on boundary, or just over or
  //       whatever (I think things are set up so equality is counted as being inside volume,
  //       but there may be edge cases (pun realized) that dont obey this.
  
// When debugging we will grab a static mutex so we dont get jumbled stdout  
#if( DEBUG_RAYTRACE_CALCS )
  std::lock_guard<std::recursive_mutex> scoped_lock( s_stdout_raytrace_mutex );
#endif
  
  
  assert( radius >= 0.0 );
  assert( half_length >= 0.0 );
  assert( (direction == CylExitDir::TowardDetector)
          || (direction == CylExitDir::AwayFromDetector) );
  
  // A convenience function for handling case where the line never enters our volume.
  auto handle_line_outside_volume = [&exit_point,&source]() -> double {
    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = source[2];
    
#if( DEBUG_RAYTRACE_CALCS )
    cout << "\tDoes not intersect volume" << endl << endl;
#endif
    
    return 0.0;
  };//handle_line_outside_volume lamda
  
  
  if( (radius <= 0.0) || (half_length <= 0.0) )
    return handle_line_outside_volume();
  
#if( DEBUG_RAYTRACE_CALCS )
  cout << "Rad=" << radius << ", half-z=" << half_length
       << ", src={" << source[0] << "," << source[1] << "," << source[2] << "}"
       << ", det={" << detector[0] << "," << detector[1] << "," << detector[2] << "}" << endl;
#endif
  

  // Get unit direction vector from source to final position
  double unit[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };
  
  {// begin scope on norm
    const double norm = sqrt( unit[0]*unit[0] + unit[1]*unit[1] + unit[2]*unit[2] );
    unit[0] /= norm;
    unit[1] /= norm;
    unit[2] /= norm;
  }// end scope on norm
  
  // Check if parallel to z-axis
  //  We should probably compare to DBL_MIN, but realistically anything less than DBL_EPSILON
  //  is close enough to zero for our purposes (the DRF will fail far before this assumption fails)
  if( (fabs(unit[0]) < DBL_EPSILON) && (fabs(unit[1]) < DBL_EPSILON) )
  {
    assert( fabs(unit[2]) > DBL_EPSILON );
    
    const double r = sqrt(source[0]*source[0] + source[1]*source[1]);
    
    if( r > radius )
      return handle_line_outside_volume();
    
    double exit_z = (unit[2] > 0.0) ? half_length : -half_length;
    switch( direction )
    {
      case CylExitDir::TowardDetector:
        break;
        
      case CylExitDir::AwayFromDetector:
        exit_z *= -1.0;
        break;
    }//switch( direction )
    
   
    const double distance = fabs( exit_z - source[2] );
    
    exit_point[0] = source[0];
    exit_point[1] = source[1];
    exit_point[2] = exit_z;
    
#if( DEBUG_RAYTRACE_CALCS )
    cout << "\tParrallel on z\n"
         << "\tExit={" << exit_point[0] << "," << exit_point[1] << "," << exit_point[2] << "}\n"
         << "\tDistance=" << distance << endl << endl;
#endif
    
    return distance;
  }//if( (fabs(unit[0]) < DBL_EPSILON) && (fabs(unit[1]) < DBL_EPSILON) )
  
  
  // Make sure both source and detector z-coordinates are not on the same end of the cylinder,
  //  and larger than the half length.
  if( (signbit(source[2]) == signbit(detector[2]))
     && (fabs(source[2]) > half_length)
     && (fabs(detector[2]) > half_length) )
  {
    return handle_line_outside_volume();
  }
  
  
  // Check that the source point isnt on the same side of the circle as the detector, but outside
  //  of the circles radius.  This will prevent the case where the infinite line would intersect
  //  the cylinder, but not between source and detector (which we should return zero for)
  // TODO: get rid of these sqrts, and also probably other ones throughout this function
  const double src_rad = sqrt( source[0]*source[0] + source[1]*source[1] );
  if( src_rad >= radius )
  {
    //There is probably a better way to check if source and detector are within 90 degrees of each
    //  other
    const double src_unit_x = source[0] / src_rad;
    const double src_unit_y = source[1] / src_rad;
    
    const double det_rad = sqrt( detector[0]*detector[0] + detector[1]*detector[1] );
    const double det_unit_x = detector[0] / det_rad;
    const double det_unit_y = detector[1] / det_rad;
    
    const double unit_dx = src_unit_x - det_unit_x;
    const double unit_dy = src_unit_y - det_unit_y;
    
    const double unit_dist_2 = unit_dx*unit_dx + unit_dy*unit_dy;

    // TODO: for radius=1, half_length=1, source={1,-1,0}, detector{1,1,0}, because of numerical
    //       rounding (unit_dist_2 = (2 - 4E-16) - i.e. 2*DBL_EPSILON) we will return in this next
    //       statement, which we probably shouldnt since it is an exact intersection.  So we should
    //       probably implement a check that accounts for these rounding errors
    if( unit_dist_2 <= 2.0 )
      return handle_line_outside_volume();
  }//if( src_rad > radius )
  
  
  // Using the x-y component of unit, find where it intersects with the radius, which could be
  //  well after the half-length
  // x2 + y2 = r2
  // y = mx + c
  
  double x_exit, y_exit, z_exit;  //intersection in the direction of detector
  
  // Lets find x, y, and z
  if( detector[0] != source[0] )
  {
    //The original implementation of finding the x,y intersection is, I believe, less numerically stable than
    //  the newer one (20250719).  And actually, the original implementation is a little upgraded for one
    //  egrogiously unstable aspect where previous to 20250709, the `y_exit` computation used what is
    //  now used for `y_exit_direct` - and was quite bad for some test cases where lines just glanced the circle.
#define USE_ORIG_CYL_INTERSECT 0
    
#if( USE_ORIG_CYL_INTERSECT )
    //(1 + m2)x2 + 2cmx + c2 – r2 = 0
    const double m = unit[1] / unit[0];
    assert( fabs( m - ((detector[1] - source[1]) / (detector[0] - source[0])) ) < std::max(1.0E-6,1.0E-6 * m) );
    const double c = source[1] - m*source[0];
    double other_x_exit;
    
    {// begin scope to solve for x of intersection
      const double a_q = 1 + m*m;
      const double b_q = 2.0 * c * m;
      const double c_q = c*c - radius*radius;
      
      if( 4.0*a_q*c_q > b_q*b_q )
        return handle_line_outside_volume();
      
      const double sqrt_quantity = sqrt(b_q*b_q -4.0*a_q*c_q);
      
      const double x_sol_1 = 0.5*(-b_q + sqrt_quantity) / a_q;
      const double x_sol_2 = 0.5*(-b_q - sqrt_quantity) / a_q;
      
      // Pick the solution towards/away-from the detector
      const bool sol_1_toward = (fabs(x_sol_1 - detector[0]) < fabs(x_sol_2 - detector[0]));
      
      switch( direction )
      {
        case CylExitDir::TowardDetector:
          x_exit       = sol_1_toward ? x_sol_1 : x_sol_2;
          other_x_exit = sol_1_toward ? x_sol_2 : x_sol_1;
          break;
            
        case CylExitDir::AwayFromDetector:
          x_exit       = sol_1_toward ? x_sol_2 : x_sol_1;
          other_x_exit = sol_1_toward ? x_sol_1 : x_sol_2;
          break;
      }//switch( direction )
    }// end scope to solve for x of intersection
    
    
    y_exit = sqrt( radius*radius - x_exit*x_exit );
    // y_exit could need multiplying by -1
    const double y_exit_direct = m*x_exit + c; //This can be numerically quite unstable for lines just glancing the cylinder
    // Choose the sign that matches the expected direction
    if( y_exit_direct < 0.0 )
      y_exit = -y_exit;
#else
    double other_x_exit;
    
    {// begin scope to solve for x,y of intersection
      // Line parametrically: P(t) = source + t*(detector-source)
      // x = source.x + t*(detector.x - source.x)
      // y = source.y + t*(detector.y - source.y)
      
      const double dx = detector[0] - source[0];
      const double dy = detector[1] - source[1];
      
      // Substitute into circle equation x² + y² = r²:
      // (source.x + t*dx)² + (source.y + t*dy)² = r²
      //
      // Expanding and rearranging into quadratic form: at² + bt + c = 0
      const double a = dx*dx + dy*dy;
      const double b = 2.0*(source[0]*dx + source[1]*dy);
      const double c = source[0]*source[0] + source[1]*source[1] - radius*radius;
      
      // Calculate discriminant
      const double discriminant = b*b - 4.0*a*c;
      
      if( discriminant < 0 )
        return handle_line_outside_volume();
      
      // Calculate the parameter values for intersection points
      const double sqrt_discriminant = sqrt(discriminant);
      const double t1 = (-b - sqrt_discriminant) / (2.0*a);
      const double t2 = (-b + sqrt_discriminant) / (2.0*a);
      
      // Calculate intersection points
      const double point1[2] = {source[0] + t1*dx, source[1] + t1*dy};
      const double point2[2] = {source[0] + t2*dx, source[1] + t2*dy};
      
      // Calculate distances from source to order the points
      const double dist1_sq = (point1[0] - detector[0])*(point1[0] - detector[0]) + (point1[1] - detector[1])*(point1[1] - detector[1]);
      const double dist2_sq = (point2[0] - detector[0])*(point2[0] - detector[0]) + (point2[1] - detector[1])*(point2[1] - detector[1]);
      
      
      // Pick the solution towards/away-from the detector
      const bool sol_1_toward = (dist1_sq < dist2_sq);
      
      switch( direction )
      {
        case CylExitDir::TowardDetector:
          x_exit       = sol_1_toward ? point1[0] : point2[0];
          y_exit       = sol_1_toward ? point1[1] : point2[1];
          other_x_exit = sol_1_toward ? point2[0] : point1[0];
          break;
          
        case CylExitDir::AwayFromDetector:
          x_exit       = sol_1_toward ? point2[0] : point1[0];
          y_exit       = sol_1_toward ? point2[1] : point1[1];
          other_x_exit = sol_1_toward ? point1[0] : point2[0];
          break;
      }//switch( direction )
    }// end scope to solve for x,y of intersection
#endif //if USE_ORIG_CYL_INTERSECT / else
    
    
    // Now find the z corresponding to {x_exit, y_exit}.
    // We'll make the z-equation as a function of x, and then plug in.
    const double m_zx = unit[2] / unit[0];
    assert( fabs( m_zx - ((detector[2] - source[2]) / (detector[0] - source[0])) ) < std::max(1.0E-6,1.0E-6 * m_zx) );
    
    const double c_zx = source[2] - m_zx*source[0];
    z_exit = m_zx*x_exit + c_zx;
    
    // If the z_exit is past the half_length, lets check that the other intersection to the
    //  infinite cylinder happens either in the volume, or on the other side of our volume
    if( fabs(z_exit) > half_length )
    {
      const double other_z_exit = m_zx*other_x_exit + c_zx;
      if( (fabs(other_z_exit) > half_length) && (signbit(z_exit) == signbit(other_z_exit)) )
        return handle_line_outside_volume();
    }//if( fabs(z_exit) > half_length )
  
  }else
  {
    // If we are here, the x coordinate of the source and detector location are the same.
    // x2 + y2 = r2
    
    assert( detector[0] == source[0] );
    
    //if( x is outside of cylinder - there is no way we would intersect cylinder )
    if( detector[0] > radius )
      return handle_line_outside_volume();
    
    x_exit = detector[0];
    
    y_exit = sqrt(radius*radius - x_exit*x_exit);
    
    switch( direction )
    {
      case CylExitDir::TowardDetector:   break;
      case CylExitDir::AwayFromDetector: y_exit *= -1.0; break;
    }//switch( direction )
    
    if( unit[1] < 0.0 )
      y_exit *= -1.0;
    
    // If detector and source x and y are the same, that means is parallel to z-axis, which we've
    //  already dealt with that case
    assert( detector[1] != source[1] );
    
    const double m_zy = unit[2] / unit[1];
    assert( fabs( m_zy - ((detector[2] - source[2]) / (detector[1] - source[1])) ) < std::max(1.0E-6,1.0E-6 * m_zy) );
    
    const double c_zy = source[2] - m_zy*source[1];
    z_exit = m_zy*y_exit + c_zy;
    
    // If the z_exit is past the half_length, lets check that the other intersection to the
    //  infinite cylinder happens either in the volume, or on the other side of our volume
    if( fabs(z_exit) > half_length )
    {
      const double other_y_exit = -y_exit;
      const double other_z_exit = m_zy*other_y_exit + c_zy;
      
      if( (fabs(other_z_exit) > half_length) && (signbit(z_exit) == signbit(other_z_exit)) )
        return handle_line_outside_volume();
    }//if( fabs(z_exit) > half_length )
  }//if( detector[0] != source[0] ) / else
  
  // If we are here, we are guaranteed the line does go through our volume
  
  // If we are exiting the cylinder on the ends, we need to figure out where on those disks we exit
  if( fabs(z_exit) > half_length )
  {
    assert( detector[2] != source[2] );
    
    const double z = ((z_exit < 0.0) ? -half_length : half_length);
    
    // Just solve for x and y as a function of z, and fill those in
    const double m_xz = unit[0] / unit[2];
    assert( fabs( m_xz - ((detector[0] - source[0]) / (detector[2] - source[2])) ) < std::max(1.0E-6,1.0E-6 * m_xz) );
    const double c_xz = source[0] - m_xz*source[2];
    
    const double m_yz = unit[1] / unit[2];
    assert( fabs( m_yz - ((detector[1] - source[1]) / (detector[2] - source[2])) ) < std::max(1.0E-6,1.0E-6 * m_yz) );
    const double c_yz = source[1] - m_yz*source[2];
    
    // A sanity check that we get source x and y locations correct, when evaluating a source z
    //  location
    assert( fabs((m_xz*source[2] + c_xz) - source[0]) < std::max(1.0,fabs(source[0])) );
    assert( fabs((m_yz*source[2] + c_yz) - source[1]) < 1.0E-6*std::max(1.0,fabs(source[1])) );
    
    x_exit = m_xz*z + c_xz;
    y_exit = m_yz*z + c_yz;
    z_exit = z;
  }else
  {
    // nothing to do here, we know we are exiting along the cylinder
  }//if( fabs(z_exit) > half_length ) / else
  
  const double dx = (source[0] - x_exit);
  const double dy = (source[1] - y_exit);
  const double dz = (source[2] - z_exit);
  const double dist_scaler = sqrt( dx*dx + dy*dy + dz*dz );
  
  exit_point[0] = x_exit;
  exit_point[1] = y_exit;
  exit_point[2] = z_exit;
  
#if( DEBUG_RAYTRACE_CALCS )
  cout << "\tCylInter={" << x_exit << "," << y_exit << "," << z_exit << "}" << endl
       << "\tExit={" << exit_point[0] << "," << exit_point[1] << "," << exit_point[2] << "}\n"
       << "\tDistance=" << dist_scaler << endl << endl;
#endif
  
  return dist_scaler;
}//double cylinder_line_intersection(...)


void DistributedSrcCalc::eval_cylinder( const double xx[], const int *ndimptr,
                                           double ff[], const int *ncompptr ) const noexcept
{
  assert( m_materialIndex < m_dimensionsTransLenAndType.size() );
  assert( (m_geometry == GeometryType::CylinderSideOn)
          || (m_geometry == GeometryType::CylinderEndOn) );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_materialIndex]) == ShellType::Material );
  
  //assert( m_materialIndex == 0 );
  //if( m_materialIndex != 0 )
  //  throw runtime_error( "eval_cylinder currently only supports trace source in inner most shielding shielding" );
  
  // This just integrates a right circular cylinder
  const int ndim = (ndimptr ? (*ndimptr) : 3);
  assert( ((m_geometry == GeometryType::CylinderEndOn) && ((ndim == 2) || (ndim == 3)))
         || ((m_geometry == GeometryType::CylinderSideOn) && (ndim == 3)) );
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex]);
  const double source_outer_rad = dimensions[0];
  const double source_half_z = dimensions[1];
  const double total_height = 2.0 * source_half_z;
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_materialIndex]);
  
  
  //cuba goes from 0 to one for each dimension, so we have to scale the variables
  //  r:     goes from cylindrical inner radius, to outer radius
  //  theta: goes from 0 to 2pi
  //  z:     goes from the negative half-height to positive half-height
  
  double max_theta = 2.0 * PhysicalUnits::pi;
  
  const double r = xx[0] * source_outer_rad;
  const double theta = ((ndim==3) ? xx[1] : 0.0) * max_theta;
  const double z = total_height * (xx[((ndim==3) ? 2 : 1)] - 0.5);
  
  // If point to evaluate is within inner-cylinder, set the source term to zero and return.
  if( m_materialIndex > 0 )
  {
    const array<double,3> &inner_dims = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex-1]);
    const double &inner_rad = inner_dims[0];
    const double &inner_half_height = inner_dims[1];
    
    if( (r < inner_rad) && (fabs(z) < inner_half_height) )
    {
      ff[0] = 0.0;
      return;
    }
  }//if( m_materialIndex > 0 )
  
  
  const double j = source_outer_rad * max_theta * total_height;
  const double dV = j * r;
  
  const double eval_point[3] = { r * cos(theta), r * sin(theta), z };  //cos is hitting 5% of total function time for debug build at least
  const bool is_side_on = (m_geometry == GeometryType::CylinderSideOn);
  const double detector_pos[3] = {
    (is_side_on ? m_observationDist : 0.0), 0.0, (is_side_on ? 0.0 : m_observationDist)
  };
  
  
  double exit_point[3];
  double dist_in_cyl = cylinder_line_intersection( source_outer_rad, source_half_z, eval_point, detector_pos, CylExitDir::TowardDetector, exit_point );
  
  //- make so cylinder_line_intersection take either a near or away from detector argument, and have it return zero if it doesnt intersect at all
  //    - Look through cylinder_line_intersection and make sure it can be mostly independent of the two points actual values
  //- Accoutn for inner cylinders by getting their near and far intersections, and make sure they are closer to detector than source location
  //- Also, update the recommended dimensions txt
  
  double trans = 0.0;
  
  // Do transport through inner cylinders, and also subtract off that distance through source
  //  cylinder
  double inner_distance = 0.0;
  for( size_t i = 0; i < m_materialIndex; ++i )
  {
    const std::array<double,3> &local_dims = std::get<0>(m_dimensionsTransLenAndType[i]);
    const double local_rad = local_dims[0];
    const double local_half_z = local_dims[1];
    const double local_trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[i]);
    const ShellType type = std::get<2>(m_dimensionsTransLenAndType[i]);
    
    double local_exit_point[3];
    const double dist_to_exit = cylinder_line_intersection( local_rad, local_half_z, eval_point,
                                      detector_pos, CylExitDir::TowardDetector, local_exit_point );
    
    if( dist_to_exit <= 0.0 ) //Line doesnt intersect cylinder
    {
      // TODO: consider eval_point is exactly at local_rad - which would return zero; do we want to
      //       attribute this point to the inner or outer cylinder?
#ifndef NDEBUG
      const double dist_src_to_exit = cylinder_line_intersection( local_rad, local_half_z, eval_point,
                                    detector_pos, CylExitDir::AwayFromDetector, local_exit_point);
      const double eval_rad = sqrt(eval_point[0]*eval_point[0] + eval_point[1]*eval_point[1]);
      
      assert( (0.0 == dist_src_to_exit)
             || (fabs(eval_rad - local_rad) < (0.000001*local_rad))
             || (fabs(local_half_z - fabs(eval_point[2])) < (0.000001*local_half_z) ) );
#endif
      
      continue;
    }//if( we dont intersect, or source is exactly on radius )
    
    // Check we either exit on the end, or on radius of the cylinder
    assert( (fabs(fabs(local_exit_point[2]) - local_half_z) < 1.0E-6*std::max(1.0,local_half_z))
      || (fabs( local_rad - sqrt(local_exit_point[0]*local_exit_point[0]
                                + local_exit_point[1]*local_exit_point[1]))
            < 1.0E-9*std::max(1.0,local_rad))
    );
    
    assert( fabs(local_exit_point[2]) < (1.0 + 1.0E-9)*std::max(local_half_z,1.0) );
    
    double local_enter_point[3];
    const double dist_to_enter = cylinder_line_intersection( local_rad, local_half_z, eval_point,
                                   detector_pos, CylExitDir::AwayFromDetector, local_enter_point );
    
    assert( dist_to_enter <= dist_to_exit );
    
    // Check we either enter on the end, or on radius of the cylinder
    assert( (fabs(fabs(local_enter_point[2]) - local_half_z) < 1.0E-6*std::max(1.0,local_half_z))
      || (fabs( local_rad - sqrt(local_enter_point[0]*local_enter_point[0]
                                + local_enter_point[1]*local_enter_point[1]))
              < 1.0E-9*std::max(1.0,local_rad))
    );
    assert( fabs(local_enter_point[2]) < (1.0 + 1.0E-9)*std::max(local_half_z,1.0) );
    
    const double local_distance = dist_to_exit - dist_to_enter;
    assert( fabs(local_distance - distance(local_enter_point,local_exit_point)) < 1.0E-9*std::max(1.0,local_distance) );
    
    switch( type )
    {
      case ShellType::Material:
        trans += (local_trans_len_coef * (local_distance - inner_distance));
        break;
        
      case ShellType::Generic:
        trans += local_trans_len_coef;
        
        // Make sure zero thickness shell
        assert( !i || (local_dims == std::get<0>(m_dimensionsTransLenAndType[i-1])) );
        break;
    }//switch( type )

    inner_distance = local_distance;
  }//for( size_t i = 0; i < m_materialIndex; ++i )
  
  trans += (trans_len_coef * (dist_in_cyl - inner_distance));
  
  // Do transport through outer cylinders
  for( size_t i = m_materialIndex + 1; i < m_dimensionsTransLenAndType.size(); ++i )
  {
    const std::array<double,3> &shield_dims = std::get<0>(m_dimensionsTransLenAndType[i]);
    const double shield_outer_rad = shield_dims[0];
    const double shield_half_z = shield_dims[1];
    const double shield_trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[i]);
    const ShellType type = std::get<2>(m_dimensionsTransLenAndType[i]);
    
    switch( type )
    {
      case ShellType::Generic:
      {
        trans += shield_trans_len_coef;
        
        // Make sure zero thickness shell
        assert( shield_dims == std::get<0>(m_dimensionsTransLenAndType[i-1]) );
        
        break;
      }//case ShellType::Generic:
        
      case ShellType::Material:
      {
        // TODO: #cylinder_line_intersection will be safe to re-use exit_point for source location and exit_point in near future - make sure to update this here by getting rid of outer_exit_point
        double outer_exit_point[3];
        const double dist_in_shield = cylinder_line_intersection( shield_outer_rad, shield_half_z,
                                                                 exit_point, detector_pos,
                                                                 CylExitDir::TowardDetector, outer_exit_point );
        
        trans += (shield_trans_len_coef * dist_in_shield);
        
        exit_point[0] = outer_exit_point[0];
        exit_point[1] = outer_exit_point[1];
        exit_point[2] = outer_exit_point[2];
        
        break;
      }//case ShellType::Material:
    }//switch( type )
  }//for( loop over outer shielding )
  
  
  trans = exp( -trans );
  
  if( m_attenuateForAir )
  {
    const double dz = m_observationDist - source_half_z;
    const double air_dist = distance( exit_point, detector_pos );
    
    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )
  
  
  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    if( is_side_on )
      trans *= exp( -(source_outer_rad - r) / m_inSituRelaxationLength );
    else
      trans *= exp( -(source_half_z - z) / m_inSituRelaxationLength );
  }//if( m_isInSituExponential )
  
  
  // Finally toss in the geometric factor (e.g., 1/r2 from where we are evaluating to to detector).
  const double eval_dist_to_det = distance( eval_point, detector_pos );
  trans *= DetectorPeakResponse::fractionalSolidAngle( 2.0*m_detectorRadius, eval_dist_to_det );
  
  ff[0] = trans * dV;
}//void eval_cylinder(...)



double rectangle_exit_location( const double half_width, const double half_height,
                               const double half_depth,
                               const double source[3],
                               const double detector[3],
                               double exit_point[3] ) noexcept
{
  assert( half_width > 0.0 );
  assert( half_height > 0.0 );
  assert( half_depth > 0.0 );
  
  assert( fabs(source[0]) <= (half_width + 1.0E-12) );
  assert( fabs(source[1]) <= (half_height + 1.0E-12) );
  assert( fabs(source[2]) <= (half_depth + 1.0E-12) );
  
  assert(fabs(detector[0]) > (half_width - 1.0E-12)
         || fabs(detector[1]) > (half_height - 1.0E-12)
         || fabs(detector[2]) > (half_depth - 1.0E-12) );
  
  assert( !((source[0] == detector[0])
             && (source[1] == detector[1])
             && (source[2] == detector[2])) );
  
#ifndef NDEBUG
  // We only get here for debug builds - source and exit_point may be same array, so we'll make a
  //  copy for development checks
  const double src_copy[3] = { source[0], source[1], source[2] };
#endif

  
  double norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };
  
  // Currently (20211126) the detector will be [0.0,0.0,m_observationDist], so we know which face
  //  of the detector the ray will exit, so we could first check for this case with this commented
  //  out code, and then return an answer efficiently; need to benchmark before bothering to
  //  uncomment this special case.
  //if( (fabs(detector[0]) <= half_width) && (fabs(detector[1]) <= half_height) )
  //{
  //  const double z_frac_inside = fabs( (half_depth - source[2]) / (detector[2] - source[2]) );
  //  exit_point[0] = source[0] + z_frac_inside * norm[0];
  //  exit_point[1] = source[1] + z_frac_inside * norm[1];
  //  exit_point[2] = ((detector[2] > 0.0) ? half_depth : -half_depth); //std::copysign(half_depth,detector[2]);
  //
  //  return distance( source, exit_point );
  //}//if( we know ray is exiting the face on positive/negative z )
  
  
  const double total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;
  
  // We'll find the intersection of all the possible planes, and then choose the one that
  //  is the closest to source
  // Recall equation of a line in three-space is:
  //   x = x_0 + t*a --> source[0] + (t * norm[0])
  //   y = z_0 + t*b --> source[1] + (t * norm[1])
  //   z = z_0 + t*c --> source[2] + (t * norm[2])
  
  const double inv_slope_x = (norm[0] == 0.0) ? DBL_MAX : (1.0 / norm[0]);
  const double x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const double t_intersect_x = (x_intersect - source[0])*inv_slope_x;
  
  const double inv_slope_y = (norm[1] == 0.0) ? DBL_MAX : (1.0 / norm[1]);
  const double y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const double t_intersect_y = (y_intersect - source[1])*inv_slope_y;
  
  const double inv_slope_z = (norm[2] == 0.0) ? DBL_MAX : (1.0 / norm[2]);
  const double z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const double t_intersect_z = (z_intersect - source[2])*inv_slope_z;
  
  if( (t_intersect_x <= t_intersect_y) && (t_intersect_x <= t_intersect_z) )
  {
    // We are exiting through the plane perpendicular to x-axis
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x * norm[2]);
    
    assert( fabs( fabs((src_copy[0] + (t_intersect_x * norm[0]))) - half_width ) < half_width*1.0E-9 );
    assert( fabs(t_intersect_x - distance(src_copy, exit_point)) < 1.0E-9*std::max(1.0,t_intersect_x) );
    
    return t_intersect_x;
  }else if( t_intersect_y <= t_intersect_z )
  {
    // We are exiting through the plane perpendicular to y-axis
    exit_point[0] = source[0] + (t_intersect_y * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y * norm[2]);
    
    assert( fabs( fabs((src_copy[1] + (t_intersect_y * norm[1]))) - half_height ) < half_height*1.0E-9 );
    assert( fabs(t_intersect_y - distance(src_copy, exit_point)) < 1.0E-9*std::max(1.0,t_intersect_y) );
    
    return t_intersect_y;
  }else
  {
    // We are exiting through the plane perpendicular to z-axis
    
    exit_point[0] = source[0] + (t_intersect_z * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);
    
    assert( fabs( fabs((src_copy[2] + (t_intersect_z * norm[2]))) - half_depth ) < half_depth*1.0E-9 );
    assert( (fabs(t_intersect_z - distance(src_copy, exit_point)) < 1.0E-9*std::max(1.0,t_intersect_z)) );
    
    return t_intersect_z;
  }// if / else figure out where we are exiting.
  
  assert( 0 );
  
  return 0.0;
}//rectangle_exit_location(...)



bool rectangle_intersections( const double half_width, const double half_height,
                             const double half_depth,
                             const double source[3],
                             const double detector[3],
                             double enter_point[3],
                             double exit_point[3] ) noexcept
{
  // Only checking inputs sanity on debug builds since this is a hot-path function, and is only
  //  called from one spot, so it should be good to only check inputs on development builds.
  
  // Make sure we arent passing garbage dimensions in ever
  assert( (half_width > 0.0) && (half_height > 0.0) && (half_depth > 0.0) );
  
  // Dev test that both the source and detector points are outside of the box; only checking on
  //  debug builds because this should really be the case
  assert( (fabs(source[0]) >= (half_width - 1.0E-12))
          || (fabs(source[1]) >= (half_height - 1.0E-12))
          || (fabs(source[2]) >= (half_depth - 1.0E-12)) );
  assert( (fabs(detector[0]) >= (half_width - 1.0E-12))
         || (fabs(detector[1]) >= (half_height - 1.0E-12))
         || (fabs(detector[2]) >= (half_depth - 1.0E-12)) );
  
  // Make sure detector and source arent in same position.
  assert( (detector[0] != source[0])
         || (detector[1] != source[1])
         || (detector[2] != source[2]) );
  
  // See notes in #rectangle_exit_location about
  
  double norm[3] = { detector[0] - source[0], detector[1] - source[1], detector[2] - source[2] };
  const double total_dist = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
  norm[0] /= total_dist;
  norm[1] /= total_dist;
  norm[2] /= total_dist;
  
  const double inv_slope_x = (norm[0] == 0.0) ? DBL_MAX : (1.0 / norm[0]);
  const double x_intersect = (inv_slope_x >= 0.0) ? half_width : -half_width;
  const double t_intersect_x_det = (x_intersect - source[0])*inv_slope_x;
  const double t_intersect_x_src = (-x_intersect - source[0])*inv_slope_x;
  
  
  const double inv_slope_y = (norm[1] == 0.0) ? DBL_MAX : (1.0 / norm[1]);
  const double y_intersect = (inv_slope_y >= 0.0) ? half_height : -half_height;
  const double t_intersect_y_det = (y_intersect - source[1])*inv_slope_y;
  const double t_intersect_y_src = (-y_intersect - source[1])*inv_slope_y;
  
  const double inv_slope_z = (norm[2] == 0.0) ? DBL_MAX : (1.0 / norm[2]);
  const double z_intersect = (inv_slope_z >= 0.0) ? half_depth : -half_depth;
  const double t_intersect_z_det = (z_intersect - source[2])*inv_slope_z;
  const double t_intersect_z_src = (-z_intersect - source[2])*inv_slope_z;
  
  const bool intersects_x_src = (t_intersect_x_src >= 0.0);
  const bool intersects_y_src = (t_intersect_y_src >= 0.0);
  const bool intersects_z_src = (t_intersect_z_src >= 0.0);
  
  const bool x_before_y_src = (!intersects_y_src || (t_intersect_x_src <= t_intersect_y_src));
  const bool x_before_z_src = (!intersects_y_src || (t_intersect_x_src <= t_intersect_z_src));
  const bool y_before_z_src = (!intersects_z_src || (t_intersect_y_src <= t_intersect_z_src));
  
  if( intersects_x_src && x_before_y_src && x_before_z_src )
  {
    // We are exiting through the plane perpendicular to x-axis
    const double src_intersect_x = ((norm[0] >= 0.0) ? -half_width : half_width);
    const double src_intersect_y = source[1] + (t_intersect_x_src * norm[1]);
    const double src_intersect_z = source[2] + (t_intersect_x_src * norm[2]);
      
    if( (fabs(src_intersect_y) > half_height)
       || (fabs(src_intersect_z) > half_depth) )
    {
      assert( std::min(std::min(t_intersect_x_src,t_intersect_y_src),t_intersect_z_src)
             <= (std::min(std::min(t_intersect_x_det,t_intersect_y_det),t_intersect_z_det)+1.0E-9) );
      return false;
    }
    
    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
    
    assert( fabs( fabs((source[0] + (t_intersect_x_src * norm[0]))) - half_width ) < half_width*1.0E-9 );
  }else if( intersects_y_src && y_before_z_src )
  {
    // We are exiting through the plane perpendicular to y-axis
    const double src_intersect_x = source[0] + (t_intersect_y_src * norm[0]);
    const double src_intersect_y = ((norm[1] >= 0.0) ? -half_height : half_height);
    const double src_intersect_z = source[2] + (t_intersect_y_src * norm[2]);
    
    if( (fabs(src_intersect_x) > half_width)
       || (fabs(src_intersect_z) > half_depth) )
    {
      assert( std::min(std::min(t_intersect_x_src,t_intersect_y_src),t_intersect_z_src)
             <= (std::min(std::min(t_intersect_x_det,t_intersect_y_det),t_intersect_z_det)+1.0E-9) );
      return false;
    }
    
    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
    
    assert( fabs( fabs((source[1] + (t_intersect_y_src * norm[1]))) - half_height ) < half_height*1.0E-9 );
  }else if( intersects_z_src )
  {
    // We are exiting through the plane perpendicular to z-axis
    const double src_intersect_x = source[0] + (t_intersect_z_src * norm[0]);
    const double src_intersect_y = source[1] + (t_intersect_z_src * norm[1]);
    const double src_intersect_z = ((norm[2] >= 0.0) ? -half_depth : half_depth);
    
    if( (fabs(src_intersect_x) > half_width)
       || (fabs(src_intersect_y) > half_height) )
    {
      assert( std::min(std::min(t_intersect_x_src,t_intersect_y_src),t_intersect_z_src)
             <= (std::min(std::min(t_intersect_x_det,t_intersect_y_det),t_intersect_z_det)+1.0E-9) );
      return false;
    }
    
    enter_point[0] = src_intersect_x;
    enter_point[1] = src_intersect_y;
    enter_point[2] = src_intersect_z;
    
    assert( fabs( fabs((source[2] + (t_intersect_z_src * norm[2]))) - half_depth ) < half_depth*1.0E-9 );
  }else
  {
    return false;
  }// if / else figure out where we are exiting.
  
  const bool intersects_x_det = (t_intersect_x_det >= 0.0);
  const bool intersects_y_det = (t_intersect_y_det >= 0.0);
  const bool intersects_z_det = (t_intersect_z_det >= 0.0);
  
  const bool x_before_y_det = (!intersects_y_det || (t_intersect_x_det <= t_intersect_y_det));
  const bool x_before_z_det = (!intersects_z_det || (t_intersect_x_det <= t_intersect_z_det));
  const bool y_before_z_det = (!intersects_z_det || (t_intersect_y_det <= t_intersect_z_det));
  
  
  
  if( intersects_x_det && x_before_y_det && x_before_z_det )
  {
    // We are exiting through the plane perpendicular to x-axis
    exit_point[0] = ((norm[0] >= 0.0) ? half_width : -half_width);
    exit_point[1] = source[1] + (t_intersect_x_det * norm[1]);
    exit_point[2] = source[2] + (t_intersect_x_det * norm[2]);
    
    assert( fabs( fabs((source[0] + (t_intersect_x_det * norm[0]))) - half_width ) < half_width*1.0E-9 );
  }else if( intersects_y_det && y_before_z_det )
  {
    // We are exiting through the plane perpendicular to y-axis
    exit_point[0] = source[0] + (t_intersect_y_det * norm[0]);
    exit_point[1] = ((norm[1] >= 0.0) ? half_height : -half_height);
    exit_point[2] = source[2] + (t_intersect_y_det * norm[2]);
    
    assert( fabs( fabs((source[1] + (t_intersect_y_det * norm[1]))) - half_height ) < half_height*1.0E-9 );
  }else if( intersects_z_det )
  {
    // We are exiting through the plane perpendicular to z-axis
    exit_point[0] = source[0] + (t_intersect_z_det * norm[0]);
    exit_point[1] = source[1] + (t_intersect_z_det * norm[1]);
    exit_point[2] = ((norm[2] >= 0.0) ? half_depth : -half_depth);
    
    assert( fabs( fabs((source[2] + (t_intersect_z_det * norm[2]))) - half_depth ) < half_depth*1.0E-9 );
  }else
  {
    assert( 0 );
  }// if / else figure out where we are exiting.
  

  assert( fabs(enter_point[0]) <= (half_width + 1.0E-9) );
  assert( fabs(enter_point[1]) <= (half_height + 1.0E-9) );
  assert( fabs(enter_point[2]) <= (half_depth + 1.0E-9) );
  
  assert( fabs(exit_point[0]) <= (half_width + 1.0E-9) );
  assert( fabs(exit_point[1]) <= (half_height + 1.0E-9) );
  assert( fabs(exit_point[2]) <= (half_depth + 1.0E-9) );
  
  return true;
}//rectangle_intersections( ... )



void DistributedSrcCalc::eval_rect( const double xx[], const int *ndimptr,
               double ff[], const int *ncompptr ) const noexcept
{
  assert( m_geometry == GeometryType::Rectangular );
  assert( m_materialIndex < m_dimensionsTransLenAndType.size() );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_materialIndex]) == ShellType::Material );
  
  const int ndim = (ndimptr ? (*ndimptr) : 2);
  assert( ndim == 3 );
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex]);
  const double half_width  = dimensions[0];
  const double half_height = dimensions[1];
  const double half_depth  = dimensions[2];
  
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_materialIndex]);
  
  // Translate the [0,1.0] coordinates from Cuba, to the world coordinates we are integrating over.
  //  (note: we would get the same answer if we integrated over just half the width/height, but it
  //   also ends up taking same amount of evaluations, so we will leave fully integrating over
  //   each dimensions, to not cause problems in the future if we add offsets or whatever)
  const double eval_x = (xx[0] - 0.5) * 2.0 * half_width;
  const double eval_y = (xx[1] - 0.5) * 2.0 * half_height;
  const double eval_z = (xx[2] - 0.5) * 2.0 * half_depth;

  
  // Check to see if [eval_x,eval_y,eval_z] is inside an inner volume, and if so set value to zero
  //  and return
  if( m_materialIndex > 0 )
  {
    const array<double,3> &dims = std::get<0>(m_dimensionsTransLenAndType[m_materialIndex-1]);
    if( (fabs(eval_x) < dims[0]) && (fabs(eval_y) < dims[1]) && (fabs(eval_z) < dims[2]) )
    {
      ff[0] = 0.0;
      return;
    }
  }//if( m_materialIndex > 0 )
  
  
  const double dV = 8.0 * half_width * half_height * half_depth;
  
  const double eval_loc[3] = { eval_x, eval_y, eval_z };
  const double detector_loc[3] = { 0.0, 0.0, m_observationDist };

  double exit_point[3];
  double dist_in_src = rectangle_exit_location( half_width, half_height, half_depth,
                                               eval_loc, detector_loc, exit_point );
  
  double trans = 0.0;
  
  // Do transport through inner shielding's, and also subtract off that distance through source
  //  Note: integrating over a 4cmx4cmx4cm void takes 381 evaluations; making the inner 2x2x2cm
  //        portion an inner void, and integrating over the outer cube (with outer dimensions still
  //        4x4x4cm, but just the inner 2x2x2cm removed) takes 11049 evaluations - 30 times longer!
  double inner_rect_dist = 0.0;
  for( size_t i = 0; i < m_materialIndex; ++i )
  {
    const std::array<double,3> &dims = std::get<0>(m_dimensionsTransLenAndType[i]);
    const double trans_len_coef_shield = std::get<1>(m_dimensionsTransLenAndType[i]);
    const ShellType type = std::get<2>(m_dimensionsTransLenAndType[i]);
    
    double enter_loc[3], exit_loc[3];
    const bool intersects = rectangle_intersections( dims[0], dims[1], dims[2],
                                                    eval_loc, detector_loc, enter_loc, exit_loc );
    
    if( intersects )
    {
      const double dist = distance( enter_loc, exit_loc );
      
      switch( type )
      {
        case ShellType::Generic:
        {
          trans += trans_len_coef_shield;
          assert( !i || (dims == std::get<0>(m_dimensionsTransLenAndType[i-1])) );
          break;
        }//case ShellType::Generic:
          
        case ShellType::Material:
        {
          trans += (trans_len_coef_shield * (dist - inner_rect_dist));
          break;
        }//case ShellType::Material:
      }//switch( type )
      
      inner_rect_dist = dist;
    }//if( intersects )
  }//for( size_t i = 0; i < m_materialIndex; ++i )
  
  
  trans += (trans_len_coef * (dist_in_src - inner_rect_dist));
  
  
  // Account for additional external shielding's
  for( size_t i = m_materialIndex + 1; i < m_dimensionsTransLenAndType.size(); ++i )
  {
    const std::array<double,3> &outer_dims = std::get<0>(m_dimensionsTransLenAndType[i]);
    const double trans_len_coef_shield = std::get<1>(m_dimensionsTransLenAndType[i]);
    const ShellType type = std::get<2>(m_dimensionsTransLenAndType[i]);
  
    switch( type )
    {
      case ShellType::Generic:
      {
        trans += trans_len_coef_shield;
        assert( outer_dims == std::get<0>(m_dimensionsTransLenAndType[i-1]) );
        break;
      }
        
      case ShellType::Material:
      {
        const double dist_in_shield = rectangle_exit_location( outer_dims[0], outer_dims[1],
                                                              outer_dims[2], exit_point, detector_loc, exit_point );
        
        trans += (trans_len_coef_shield * dist_in_shield);
        break;
      }//case ShellType::Material:
    }//switch( type )
    
    
  }//for( loop over outer shieldings )
  
  
  trans = exp( -trans );
  
  if( m_attenuateForAir )
  {
    const double air_dist = distance( exit_point, detector_loc ); 
    trans *= exp( -m_airTransLenCoef * air_dist );
  }//if( m_attenuateForAir )
  
  
  if( m_isInSituExponential )
  {
    assert( m_inSituRelaxationLength > 0.0 );
    trans *= exp( -(half_depth - eval_z) / m_inSituRelaxationLength );
  }
  
  const double eval_dist_to_det = distance(eval_loc, detector_loc);
  trans *= DetectorPeakResponse::fractionalSolidAngle( 2.0*m_detectorRadius, eval_dist_to_det );
  
  ff[0] = trans * dV;
}//void eval_rect(...)


std::pair<std::shared_ptr<ShieldingSourceChi2Fcn>, ROOT::Minuit2::MnUserParameters> ShieldingSourceChi2Fcn::create(
                                const double distance,
                                const GeometryType geom,
                                const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                const std::vector<ShieldingSourceFitCalc::SourceFitDef> &src_definitions,
                                std::shared_ptr<const DetectorPeakResponse> detector,
                                std::shared_ptr<const SpecUtils::Measurement> foreground,
                                std::shared_ptr<const SpecUtils::Measurement> background,
                                std::deque<std::shared_ptr<const PeakDef>> foreground_peaks,
                                std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> background_peaks,
                                const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options )
{
  using GammaInteractionCalc::ShieldingSourceChi2Fcn;
    
  ROOT::Minuit2::MnUserParameters inputPrams;
  
  //Get the peaks we'll be using in the fit
  vector<PeakDef> peaks;
  std::set<const SandiaDecay::Nuclide *> nuclides;
  for( const shared_ptr<const PeakDef> &peak : foreground_peaks )
  {
    if( peak && peak->useForShieldingSourceFit() )
    {
      peaks.push_back( *peak );
      nuclides.insert( peak->parentNuclide() );
    }
  }//for(...)
  
  if( peaks.empty() )
    throw runtime_error( "There are no peaks selected for the fit" );
  
  double liveTime = foreground ? foreground->live_time() : 1.0f;
  double realTime = foreground ? foreground->real_time() : 0.0f;
  
  if( liveTime <= 0.0 )
  {
    passMessage( "There was no defined detector live time, so assuming 300 seconds",
                WarningWidget::WarningMsgHigh );
    realTime = liveTime = 300.0 * PhysicalUnits::second;
  }//if( liveTime <= 0.0 )
    
  
  shared_ptr<GammaInteractionCalc::ShieldingSourceChi2Fcn> answer;
  answer.reset( new GammaInteractionCalc::ShieldingSourceChi2Fcn( distance, liveTime, realTime,
                                                    peaks, detector, shieldings, geom, options ) );
  
  //I think num_fit_params will end up same as inputPrams.VariableParameters()
  size_t num_fit_params = 0;
  
  // TODO: need to check that if we are fitting a shielding thickness we have an
  //  isotope that has more than one peak in the fit, and if it is also fitting
  //  age should have more than two peaks.  There are probably a few more things
  //  to check here...
  
  //Setup the parameters from the sources
  num_fit_params += answer->setInitialSourceDefinitions( src_definitions, inputPrams );
  
  
  //setup the parameters for the shieldings
  for( size_t i = 0; i < shieldings.size(); ++i )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &info = shieldings[i];
    
    if( info.m_isGenericMaterial )
    {
      const double an = info.m_dimensions[0]; //select->atomicNumber();
      const double ad = info.m_dimensions[1]; //select->arealDensity();
      const bool fitAn = info.m_fitDimensions[0];
      const bool fitAD = info.m_fitDimensions[1];
      
      const string name = "Generic_" + std::to_string(i);
      const double adUnits = PhysicalUnits::g/PhysicalUnits::cm2;
      
      num_fit_params += fitAn + fitAD;
      
      if( fitAn )
        inputPrams.Add( name + "_AN", an, std::max(0.1*an,2.5),
                       1.0*MassAttenuation::sm_min_xs_atomic_number,
                       1.0*MassAttenuation::sm_max_xs_atomic_number );
      else
        inputPrams.Add( name + "_AN_FIXED", an );
      
      if( fitAD )
        inputPrams.Add( name + "_AD", ad, std::max(5.0*adUnits, 0.1*ad), 0.0, 400.0*adUnits );  //400g/cm2 is about 35cm Pb
      else
        inputPrams.Add( name + "_AD", ad );
      
      inputPrams.Add( name + "_dummyshielding2", 0.0 );
    }else
    {
      std::shared_ptr<const Material> mat = info.m_material;
      string name;
      if( mat )
        name = mat->name + std::to_string(i);
      else
        name = "unspecifiedmat_" + std::to_string(i);
      
      switch( geom )
      {
        case GeometryType::Spherical:
        {
          const double thickness = info.m_dimensions[0];
          const bool fitThickness = mat ? info.m_fitDimensions[0] : false;
          
          num_fit_params += fitThickness;
          
          if( fitThickness )
            inputPrams.Add( name + "_thickness", thickness, std::max(10.0*PhysicalUnits::mm,0.25*thickness), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_thickness", thickness );
          
          inputPrams.Add( name + "_dummyshielding1", 0.0 );
          inputPrams.Add( name + "_dummyshielding2", 0.0 );
          
          break;
        }//case GeometryType::Spherical:
          
        case GeometryType::CylinderEndOn:
        case GeometryType::CylinderSideOn:
        {
          const double rad = info.m_dimensions[0];
          const double len = info.m_dimensions[1];
          
          const bool fitRad = mat ? info.m_fitDimensions[0] : false;
          const bool fitLen = mat ? info.m_fitDimensions[1] : false;
          
          num_fit_params += fitRad;
          num_fit_params += fitLen;
          
          if( fitRad )
            inputPrams.Add( name + "_dr", rad, std::max(2.5*PhysicalUnits::mm,0.25*rad), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_dr", rad );
          
          if( fitLen )
            inputPrams.Add( name + "_dz", len, std::max(2.5*PhysicalUnits::mm,0.25*len), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_dz", len );
          
          inputPrams.Add( name + "_dummyshielding2", 0.0 );
          
          break;
        }//case GeometryType::CylinderEndOn and CylinderSideOn:
          
          
        case GeometryType::Rectangular:
        {
          const double width = info.m_dimensions[0];
          const double height = info.m_dimensions[1];
          const double depth = info.m_dimensions[2];
          
          const bool fitWidth = mat ? info.m_fitDimensions[0] : false;
          const bool fitHeight = mat ? info.m_fitDimensions[1] : false;
          const bool fitDepth = mat ? info.m_fitDimensions[2] : false;
          
          num_fit_params += fitWidth;
          num_fit_params += fitHeight;
          num_fit_params += fitDepth;
          
          if( fitWidth )
            inputPrams.Add( name + "_dx", width, std::max(2.5*PhysicalUnits::mm,0.25*width), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_dx", width );
          
          if( fitHeight )
            inputPrams.Add( name + "_dy", height, std::max(2.5*PhysicalUnits::mm,0.25*height), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_dy", height );
          
          if( fitDepth )
            inputPrams.Add( name + "_dz", depth, std::max(2.5*PhysicalUnits::mm,0.25*depth), 0, 1000.0*PhysicalUnits::m );
          else
            inputPrams.Add( name + "_dz", depth );
          
          break;
        }//case GeometryType::Rectangular:
          
        case GeometryType::NumGeometryType:
        {
          assert( 0 );
          break;
        }//case GeometryType::NumGeometryType:
      }//switch( geometry() )
    }//if( generic material ) / else
  }//for( size_t i = 0; i < shieldings.size(); ++i )
  
  //  if( num_fit_params < 1 )
  //    throw runtime_error( "There is nothing being fit for" );
  
  if( background && background_peaks && !background_peaks->empty() )
  {
    if( options.background_peak_subtract )
    {
      vector<PeakDef> backgroundpeaks;
      for( const auto &p : *background_peaks )
        backgroundpeaks.push_back( *p );
      answer->setBackgroundPeaks( backgroundpeaks, background->live_time() );
    }else
    {
      cerr << __FUNCTION__ << ": background peaks were passed in, but options said background peak subtraction was not wanted!" << endl;
      // In principle, this is fine - I just want to make sure things were being treated consistently before
      //  the check for `options.background_peak_subtract` was added 20240917
      assert( 0 );
    }
  }//if( background && background_peaks && !background_peaks->empty() )
  
  
  for( size_t shielding_index = 0; shielding_index < shieldings.size(); ++shielding_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = shieldings[shielding_index];
    
    for( const auto &el_nucs : shield.m_nuclideFractions_ )
    {
      const SandiaDecay::Element * const el = el_nucs.first;
      const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &nuc_frac_fits = el_nucs.second;
      
      assert( el );
      
      // If we arent fitting any mass fractions, we wont need add any parameters for this element
      bool is_fitting_any_el = false;
      
      //Get the isotopes to fit mass fractions of, for this element, and their initial values
      double fracmaterial_fit = 0.0;
      vector<pair<const SandiaDecay::Nuclide *,double>> massfracs;
      for( const auto &i : nuc_frac_fits )
      {
        if( std::get<2>(i) )
        {
          is_fitting_any_el = true;
          fracmaterial_fit += std::get<1>(i);
          massfracs.emplace_back( std::get<0>(i), std::get<1>(i) );
        }
      }//for( const auto &i : nuc_frac_fits )
      
      if( !is_fitting_any_el )
        continue; //just keep in mind may still be a self-atten source, that isnt being fit
      
      shared_ptr<const Material> mat = shield.m_material;
      if( !mat )
        throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...)"
                            " serious logic error when fitting for mass frac - invalid material");
      
      if( massfracs.size() == 1 )
      {
        string msg = "Only a single nuclide of " + (el ? el->symbol : "nullptr") + " is selected"
        " to vary (" + (massfracs[0].first ? massfracs[0].first->symbol : string("other nuclides"))
        + ") - you need to select either multiple nuclides for each element, or one nuclide"
        " and vary the non-source nuclides.";
        
        passMessage( msg, WarningWidget::WarningMsgHigh );
        throw runtime_error( "Error fitting mass fraction - " + msg );
      }//if( massfracs.size() == 1 )
      
      if( fracmaterial_fit <= 0.0 )
      {
        passMessage( "When fitting for mass fractions of source nuclides, the "
                    "sum of the fit for mass fractions equal the sum of the "
                    "initial values, therefore the initial sum of mass "
                    "fractions must be greater than 0.0",
                    WarningWidget::WarningMsgHigh );
        throw runtime_error( "Error fitting mass fraction - sum fit fractions zero or less." );
      }//if( fracmaterial <= 0.0 )
      
      double usedmassfrac = 0.0;
      const size_t num_frac_fit_pars = (massfracs.empty() ? 0 : massfracs.size()-1);
      for( size_t frac_index = 0; frac_index < num_frac_fit_pars; ++frac_index )
      {
        const pair<const SandiaDecay::Nuclide *,double> &nmf = massfracs[frac_index];
        
        string name = mat->name + "_" + (nmf.first ? nmf.first->symbol : (el->symbol + "_other"))
                      + "_" + std::to_string(shielding_index);
        double val = 0.0;
        const double remaining_frac = (fracmaterial_fit - usedmassfrac);
        if( remaining_frac > nmf.second )
          val = nmf.second / remaining_frac;
        
        usedmassfrac += nmf.second;
        inputPrams.Add( name, val, max(0.1*val,0.01), 0, 1.0 );
        ++num_fit_params;
      }//for( size_t j = 0; i < nmassfrac; ++j )
    }//for( const auto &el_nucs : shield.m_nuclideFractions_ )
  }//for( size_t shielding_index = 0; shielding_index < shieldings.size(); ++shielding_index )
  
  
  if( num_fit_params != inputPrams.VariableParameters() )
    throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...): "
                        "there is a serious logic error in this function, "
                        "please let wcjohns@sandia.gov know about this." );
  
  const size_t num_expected = answer->numExpectedFitParameters();
  const size_t num_input_pars = inputPrams.Parameters().size();
  assert( num_expected == num_input_pars );
  if( num_expected != num_input_pars )
    throw logic_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...): "
                      "mismatch between num expected fit parameters, and number actual fit pars." );
  return {answer, inputPrams};
}//pair<shared_ptr<ShieldingSourceChi2Fcn>,ROOT::Minuit2::MnUserParameters> create(...)
  

//This class evaluated the chi2 of a given hypothesis, where it is assumed the
//  radioactive source is a point source located at the center of concentric
//  spherical shells consisting of non-radioactive materials, that may either
//  be defined via the 'Material' class, or be a generic material defined by
//  atomic number and areal density which is indicated by a NULL pointer.
//
//paramaters:
//-Activity nuclide 0, where nuclides are sorted alphebaetically by name
//-Age of nuclide 0
//-Activity nuclide 1
//-Age of nuclide 1
// ...
//if material 0 normal material (if Material* is non-NULL pointer)
//    -Material Thinckness
//else if generic material (if Material* is NULL)
//    -atomic number
//    -areal density
//if material 0 normal material (if Material* is non-NULL pointer)
//    -Material Thinckness
//else if generic material (if Material* is NULL)
//    -atomic number
//    -areal density
// ...

ShieldingSourceChi2Fcn::ShieldingSourceChi2Fcn(
                                 const double distance, const double liveTime, const double realTime,
                                 const std::vector<PeakDef> &peaks,
                                 std::shared_ptr<const DetectorPeakResponse> detector,
                                 const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                 const GeometryType geometry,
                                 const ShieldingSourceFitCalc::ShieldingSourceFitOptions &options )
  : ROOT::Minuit2::FCNBase(),
    m_cancel( CalcStatus::NotCanceled ),
    m_isFitting( false ),
    m_distance( distance ),
    m_liveTime( liveTime ),
    m_peaks( peaks ),
    m_detector( detector ),
    m_initial_shieldings( shieldings ),
    m_nuclides( 0, (const SandiaDecay::Nuclide *)NULL ),
    m_geometry( geometry ),
    m_options( options ),
    m_realTime( realTime )
{
  set<const SandiaDecay::Nuclide *> nucs;
  for( const PeakDef &p : m_peaks )
  {
    if( p.parentNuclide() )
      nucs.insert( p.parentNuclide() );
  }
      
  for( const SandiaDecay::Nuclide *n : nucs )
    m_nuclides.push_back( n );
  
  std::sort( m_nuclides.begin(), m_nuclides.end(),
            []( const SandiaDecay::Nuclide *lhs, const SandiaDecay::Nuclide *rhs ) -> bool {
              if( !lhs ) return false;
              if( !rhs ) return true;
              return (lhs->symbol < rhs->symbol);
  } );
    
  // Go through and sanity-check self attenuating and trace sources
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    const shared_ptr<const Material> &material = info.m_material;
    
    const map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,bool>>> &nuclideFractions
                                                                        = info.m_nuclideFractions_;
    size_t num_self_atten = 0;
    for( const auto &el_nucs : nuclideFractions )
      num_self_atten += el_nucs.second.size();
    
    if( !material && (num_self_atten || !info.m_traceSources.empty()) )
    {
      throw runtime_error( "ShieldingSourceChi2Fcn: Self attenuating and trace sources must"
                          " be empty for a generic shielding" );
    }
    
    // Check trace sources.
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      const SandiaDecay::Nuclide *nuc = trace.m_nuclide;
      if( !nuc )
        throw runtime_error( "ShieldingSourceChi2Fcn: null trace source" );
      
      if( !nucs.count(nuc) )
        throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                            + nuc->symbol + " doesnt have a peak with that as an assigned nuclide" );
    }//for( loop over trace sources )
    
    // Check self-attenuating sources
    std::map<short,double> self_atten_mass_fracs;
    for( const auto &el_nucs : nuclideFractions )
    {
      const SandiaDecay::Element * const el = el_nucs.first;
      const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &els_nucs = el_nucs.second;
      assert( el );
      assert( !els_nucs.empty() );
      if( nucs.empty() )
        continue;
      
      if( !el )
        throw std::logic_error( "ShieldingSourceChi2Fcn: got a nullptr element for self-atten nuclide list" );
      
      size_t num_fit = 0;
      set<const SandiaDecay::Nuclide *> nucs_seen;
      for( const auto &nuc_frac_fit : els_nucs )
      {
        assert( material );
        
        const SandiaDecay::Nuclide * const nuc = std::get<0>(nuc_frac_fit);
        assert( !nuc || (nuc->atomicNumber == el->atomicNumber) );
        
        const double mass_frac = std::get<1>(nuc_frac_fit);
        const bool fit = std::get<2>(nuc_frac_fit);
        
        if( nucs_seen.count(nuc) )
          throw runtime_error( "ShieldingSourceChi2Fcn: duplicate self-atten nuclide, "
                              + (nuc ? nuc->symbol : string("other fraction")) + ", seen." );
        
        nucs_seen.insert(nuc);
        num_fit += fit;
        
        if( nuc && !nucs.count(nuc) )
          throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                              + nuc->symbol + " doesnt have a peak with that as an assigned nuclide" );
        
        if( !self_atten_mass_fracs.count(el->atomicNumber) )
          self_atten_mass_fracs[el->atomicNumber] = 0.0;
        self_atten_mass_fracs[el->atomicNumber] += mass_frac;
        
        int hasElement = 0;
        for( const auto &p : material->elements )
          hasElement += (p.first->atomicNumber == el->atomicNumber);
        
        for( const auto &p : material->nuclides )
          hasElement += nuc ? (p.first == nuc) : 0;
        
        if( !hasElement )
          throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                              + (nuc ? nuc->symbol : "other nuclides") + " not in shielding " + material->name );
      }//for( const auto &nuc_frac_fit : nucs )
      
      if( num_fit == 1 )
      {
        throw runtime_error( "ShieldingSourceChi2Fcn: element " + el->name + " is fitting for a"
                            " single nuclide fraction. You must fit for either no nuclides in an"
                            " element, or more than one." );
      }
    }//for( const auto &el_nucs : nuclideFractions )
  }//for( const ShieldingSourceFitCalc::ShieldingInfo &info : shieldings )
}//ShieldingSourceChi2Fcn


ShieldingSourceChi2Fcn::~ShieldingSourceChi2Fcn()
{
  {//begin lock on m_zombieCheckTimerMutex
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
    
    if( m_zombieCheckTimer )
    {
      boost::system::error_code ec;
      m_zombieCheckTimer->cancel( ec );
      if( ec )
        cerr << "ShieldingSourceChi2Fcn::fittingIsFinished(): Got error cancelling m_zombieCheckTimer: "
        << ec.message() << endl;
      m_zombieCheckTimer.reset();
    }//if( m_zombieCheckTimer )
  }//end lock on m_zombieCheckTimerMutex
}


const GeometryType ShieldingSourceChi2Fcn::geometry() const
{
  return m_geometry;
}


void ShieldingSourceChi2Fcn::cancelFit()
{
  m_cancel = CalcStatus::UserCanceled;
}
  
  
void ShieldingSourceChi2Fcn::cancelFitWithNoUpdate()
{
  m_cancel = CalcStatus::CanceledNoUpdate;
}
  


void ShieldingSourceChi2Fcn::setGuiProgressUpdater( std::shared_ptr<GuiProgressUpdateInfo> updateInfo )
{
  m_guiUpdateInfo = updateInfo;
}//void setGuiUpdater(...)
  
  
void ShieldingSourceChi2Fcn::zombieCallback( const boost::system::error_code &ec )
{
  if( ec )  //timer was cancelled
    return;
  
  m_cancel = CalcStatus::Timeout;

  cerr << "In Zombie callback!" << endl;
}//void zombieCallback( const boost::system::error_code &ec )
  
  
void ShieldingSourceChi2Fcn::fittingIsStarting( const size_t deadlineMs )
{
  m_isFitting = true;
  
  Wt::WServer *server = Wt::WServer::instance();
  if( deadlineMs && server )
  {
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
    m_zombieCheckTimer = make_shared<boost::asio::deadline_timer>( server->ioService() );
    m_zombieCheckTimer->expires_from_now( boost::posix_time::milliseconds(deadlineMs) );
    m_zombieCheckTimer->async_wait( [this](const boost::system::error_code &ec){ zombieCallback(ec); } );
  }//if( deadlineMs )

  if( m_guiUpdateInfo )
    m_guiUpdateInfo->fitting_starting();
}//void fittingIsStarting()

  
void ShieldingSourceChi2Fcn::fittingIsFinished()
{
  {//begin lock on m_zombieCheckTimerMutex
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
  
    if( m_zombieCheckTimer )
    {
      boost::system::error_code ec;
      m_zombieCheckTimer->cancel( ec );
      if( ec )
        cerr << "ShieldingSourceChi2Fcn::fittingIsFinished(): Got error cancelling m_zombieCheckTimer: "
             << ec.message() << endl;
      m_zombieCheckTimer.reset();
    }//if( m_zombieCheckTimer )
  }//end lock on m_zombieCheckTimerMutex
  
  m_isFitting = false;
}//void fittingIsFinished()
  

void ShieldingSourceChi2Fcn::setSelfAttMultiThread( const bool do_multithread )
{
  m_options.multithread_self_atten = do_multithread;
}
  
  
size_t ShieldingSourceChi2Fcn::setInitialSourceDefinitions(
                        const std::vector<ShieldingSourceFitCalc::SourceFitDef> &src_definitions,
                        ROOT::Minuit2::MnUserParameters &inputPrams )
{
  assert( m_initialSrcDefinitions.empty() );
  if( !m_initialSrcDefinitions.empty() )
    throw runtime_error( "ShieldingSourceChi2Fcn::setInitialSourceDefinitions: already called." );
 
  const size_t niso = src_definitions.size();
  if( niso < 1 )
    throw runtime_error( "There are no isotopes being fit for" );
  
  
  set<const SandiaDecay::Nuclide *> unique_nucs;
  for( const auto &i : src_definitions )
  {
    if( !i.nuclide || i.nuclide->isStable() )
      throw runtime_error( "Invalid input nuclide set as source" );
    unique_nucs.insert( i.nuclide );
  }
  
  // We _could_ support having a nuclide as a point source + volumetric source, or
  //  even as a volumetric source in multiple shieldings, but we wont for now.
  if( unique_nucs.size() != src_definitions.size() )
    throw runtime_error( "Source nuclide specified more than once as source" );
  
  
  if( src_definitions.size() != numNuclides() )
    throw runtime_error( "In setting initial source definitions, there are not the same number"
                        " of source definitions (" + std::to_string(src_definitions.size())
                        + ") as nuclides (" + std::to_string(numNuclides())
                        + ") in peaks being used." );
  
  
  size_t num_fit_params = 0;
  
  for( size_t i = 0; i < src_definitions.size(); ++i )
  {
    const SandiaDecay::Nuclide * const nuclide = this->nuclide(i);
    assert( nuclide );
    if( !nuclide )
      throw std::logic_error( "Invalid nuclide in ShieldingSourceChi2Fcn" );
    
    const ShieldingSourceFitCalc::SourceFitDef *srcdef = nullptr;
    for( const auto &src : src_definitions )
    {
      if( src.nuclide == nuclide )
      {
        srcdef = &src;
        break;
      }
    }//for( const auto &src : src_definitions )
    
    if( !srcdef )
      throw runtime_error( "There was a peak with nuclide " + nuclide->symbol
                          + " but no SourceFitDef for this nuclide" );
    
    double activity = srcdef->activity;
    
    //We could do a lot better here by estimating the activity of the sources
    //  the first time they are fit for
    //XXX - should make better guesses for source activity the first time a fit
    //      is performed
    
    //    cerr << "Initial activity is " << m_sourceModel->activity( ison )/PhysicalUnits::curie
    //         << " which is a minuit value of " << activity << endl;
    //    activity = 10.0*PhysicalUnits::curie*(1.0E-6) / ShieldingSourceChi2Fcn::sm_activityUnits;
    
    assert( !srcdef->fitActivity
           || (srcdef->sourceType != ShieldingSourceFitCalc::ModelSourceType::Intrinsic) );
    
    const bool fitAct = srcdef->fitActivity && (srcdef->sourceType != ShieldingSourceFitCalc::ModelSourceType::Intrinsic);
    const bool fitAge = srcdef->fitAge;
    
    const SandiaDecay::Nuclide *ageDefiningNuc = srcdef->ageDefiningNuc;
    const bool hasOwnAge = (!ageDefiningNuc || (ageDefiningNuc == nuclide));
    
    num_fit_params += fitAct + (fitAge && hasOwnAge);
    
    // Put activity into units of ShieldingSourceChi2Fcn
    activity /= ShieldingSourceChi2Fcn::sm_activityUnits;
    
    if( fitAct )
    {
      //We cant have both a specified lower and upper limit on activity, or
      //  such a large range will make Minuit2 choke and give a completely
      //  in-accurate answer (returns not even the best chi2 it found), if only
      //  one fit parameter.
      //      inputPrams.Add( nuclide->symbol + "Strength", activity, activityStep, 0.0,
      //                     10000.0*PhysicalUnits::curie/ShieldingSourceChi2Fcn::sm_activityUnits );
      const string name = nuclide->symbol + "Strength";
      const double activityStep = (activity < 0.0001 ? 0.0001 : 0.1*activity);
      inputPrams.Add( name, activity, activityStep );
      inputPrams.SetLowerLimit( name, 0.0 );
    }else
    {
      inputPrams.Add( nuclide->symbol + "Strength", activity );
    }
    
    
    if( fitAge && hasOwnAge )
    {
      //We could do a lot better on creating the age range - there must be some way to easily
      //  determine how old an isotope has to get before it essentially doesn't change any more
      //  (prompt HL, etc.).  I guess we could look at the peaks being used to fit for and age them
      //  until their ratios don't change within the available statistical precision.
      //  But for the moment, we'll do something much simpler and use the maximum of either the
      //  longest progenies half-life, or the sum half life of all progeny
      //  \TODO: improve this max decay time estimate; some possibilities are:
      //    - Could probably ignore the parents half-life, or only partially take into account
      //    - For common nuclides could define reasonable fixed values
      //    - Could look at gamma spectrum produced over time, and pick the time when the selected
      //      photopeak ratios change little enough as to not be statistically significant to the
      //      observed data (or even just hard-coded limits).
      auto maxNuclideDecayHL = []( const SandiaDecay::Nuclide * const nuc ) -> double {
        double maxhl = 0.0, sumhlfs = 0.0;
        
        for( auto n : nuc->descendants() )
        {
          if( !n->isStable() )
          {
            sumhlfs += n->halfLife;
            maxhl = std::max( maxhl, n->halfLife );
          }
        }//for( auto n : nuc->descendants() )
        
        //cout << "For nuc=" << nuc->symbol << " maxhl=" << PhysicalUnits::printToBestTimeUnits(maxhl)
        //     << ", sumhlfs=" << PhysicalUnits::printToBestTimeUnits(sumhlfs)
        //     << " - will set max age to " << PhysicalUnits::printToBestTimeUnits(std::max( 7*maxhl, 4*sumhlfs ))
        //     << endl;
        
        //return 100.0*nuc->halfLife;
        return std::max( 7*maxhl, 4*sumhlfs );
      };//maxNuclideDecayHL
      
      double maxAge = -1.0;
      if( ageDefiningNuc == nuclide )
      {
        //We are
        
        for( size_t trialInd = 0; trialInd < this->numNuclides(); ++trialInd )
        {
          const SandiaDecay::Nuclide *trialNuc = this->nuclide( int(trialInd) );
          if( trialNuc->atomicNumber == nuclide->atomicNumber )
          {
            const double thisMaxAge = maxNuclideDecayHL( nuclide );
            maxAge = std::max( maxAge, thisMaxAge );
          }
        }//for( loop over all nuclides being fit for )
      }else
      {
        maxAge = maxNuclideDecayHL( nuclide );
      }
      assert( maxAge > 0.0 );
      
      double age = srcdef->age;
      double ageStep = 0.25 * nuclide->halfLife;
      
      // Limit the maximum age to be the larger of ten times the current age, or 200 years.  This
      //  is both to be reasonable in terms of answers we get, and because for really long-lived
      //  isotopes, we could have a max age at this point so large it will cause Minuit to give
      //  NaNs for age, even on first iteration.
      maxAge = std::min( maxAge, std::max(10.0*age, 200.0*PhysicalUnits::year) );
      
      // If the age is currently over 10000 years, it is just really getting unreasonable, so
      //  larger than Minuit can handle, so will impose a tougher 100k year limit, but let grow past
      //  this, but require the user to hit "fit" over and over again.
      if( maxAge > 10000.0*PhysicalUnits::year )
        maxAge = std::max(2.0*age, 10000.0*PhysicalUnits::year);
      
      // But no matter what we'll limit to the age of earth, which at least for a few select
      //  examples tried, Minuit was okay with (it wasnt okay with like 1.2E20 years that some of
      //  the uraniums would give)
      age = std::min( age, 4.543e+9 * PhysicalUnits::year );
      maxAge = std::min( maxAge, 4.543e+9 * PhysicalUnits::year );
      
      ageStep = std::min( ageStep, 0.1*maxAge );
      if( age > 0 )
        ageStep = std::min( 0.1*age, ageStep );
      
      //cout << "For nuclide " << nuclide->symbol << " adding age=" << age << ", with step " << ageStep << " and max age " << maxAge << endl;
      
      inputPrams.Add( nuclide->symbol + "Age", age, ageStep, 0, maxAge  );
    }else if( hasOwnAge )
    {
      const double age = srcdef->age;
      //cout << nuclide->symbol << " has own non-fitting age going in as param " << inputPrams.Parameters().size() << endl;
      inputPrams.Add( nuclide->symbol + "Age", age );
    }else  //see if defining nuclide age is fixed, if so use it, else put in negative integer of index of age...
    {
      assert( ageDefiningNuc );
      
      // Previous to 20210825, we were doing the next commented-out line, which is not correct; we
      //  want the relative nuclide index to be for the order nuclides are added into 'inputPrams',
      //  which may be different m_sourceModel->m_nuclides.
      //const int age_defining_index = m_sourceModel->nuclideIndex( ageDefiningNuc );
      
      int age_defining_index = -1;
      for( size_t ansnucn = 0; ansnucn < this->numNuclides(); ++ansnucn )
      {
        const SandiaDecay::Nuclide *answnuclide = this->nuclide( static_cast<int>(ansnucn) );
        if( answnuclide == ageDefiningNuc )
        {
          age_defining_index = static_cast<int>(ansnucn);
          break;
        }
      }//for( size_t ansnucn = 0; ansnucn < answer->numNuclides(); ++ansnucn )
      
      //assert( age_defining_index >= 0 );
      if( age_defining_index < 0 )  //shouldnt ever happen, but JIC
        throw runtime_error( "Error finding age defining nuclide for " + nuclide->symbol
                            + " (should have been " + ageDefiningNuc->symbol + ")" );
      
      const double ageIndexVal = -1.0*(age_defining_index + 1);
      //cout << nuclide->symbol << ": ageIndexVal=" << ageIndexVal
      //<< " going in as param " << inputPrams.Parameters().size() << endl;
      inputPrams.Add( nuclide->symbol + "Age", ageIndexVal );
    }
  }//for( size_t src_index = 0; src_index < src_definitions.size(); ++src_index )
  
  m_initialSrcDefinitions = src_definitions;
  
  return num_fit_params;
}//size_t setInitialSourceDefinitions(...)
  
  
const std::vector<ShieldingSourceFitCalc::SourceFitDef> &ShieldingSourceChi2Fcn::initialSourceDefinitions() const
{
  return m_initialSrcDefinitions;
}
  
  
const SandiaDecay::Nuclide *ShieldingSourceChi2Fcn::nuclide( const size_t nuc_index ) const
{
  if( nuc_index >= m_nuclides.size() )
    throw std::runtime_error( "ShieldingSourceChi2Fcn::nuclide(int): invalid index" );
  return m_nuclides[nuc_index];
}//const Nuclide *nuclide( const int nucN ) const


double ShieldingSourceChi2Fcn::operator()( const std::vector<double> &x ) const
{
  return DoEval( x );
}//double operator()(...)


ShieldingSourceChi2Fcn &ShieldingSourceChi2Fcn::operator=( const ShieldingSourceChi2Fcn &rhs )
{
  m_cancel = rhs.m_cancel.load();
  m_distance = rhs.m_distance;
  m_liveTime = rhs.m_liveTime;
  m_peaks = rhs.m_peaks;
  m_detector = rhs.m_detector;
  m_nuclides = rhs.m_nuclides;
  m_options = rhs.m_options;
  
  //m_isFitting
  //m_guiUpdateInfo
  //m_zombieCheckTimer
  
  return *this;
}//operator=


double ShieldingSourceChi2Fcn::Up() const
{
  return 1.0;
}//double Up();

  
bool ShieldingSourceChi2Fcn::isVariableMassFraction( const size_t material_index,
                                        const SandiaDecay::Nuclide *nuc ) const 
{
  if( material_index >= m_initial_shieldings.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::isVariableMassFraction: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  if( !nuc )
    return false;
  
  const SandiaDecay::Element * const el = db->element( nuc->atomicNumber );
  const auto pos = shield.m_nuclideFractions_.find( el );
  if( pos == end(shield.m_nuclideFractions_) )
    return false;
  
  const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &fracs = pos->second;
  for( const auto &n : fracs )
  {
    if( std::get<0>(n) == nuc )
      return std::get<2>(n);
  }
  
  return false;
}//isVariableMassFraction(...)
  
  
bool ShieldingSourceChi2Fcn::isVariableOtherMassFraction( const size_t material_index,
                               const SandiaDecay::Element *el ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::isVariableOtherMassFraction: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
  
  const auto pos = shield.m_nuclideFractions_.find( el );
  if( pos == end(shield.m_nuclideFractions_) )
    return false;
  
  const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> &fracs = pos->second;
  for( const auto &n : fracs )
  {
    if( !std::get<0>(n) ) //Look for nullptr nuclide, and return its value
      return std::get<2>(n);
  }
  
  return false;
}//isVariableOtherMassFraction(...)
  
  
bool ShieldingSourceChi2Fcn::hasVariableMassFraction( const size_t material_index ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::hasVariableMassFraction: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
  
  for( const auto &el_nucs : shield.m_nuclideFractions_ )
  {
    for( const auto &nucs : el_nucs.second )
    {
      assert( el_nucs.first );
      assert( !get<0>(nucs) || (get<0>(nucs)->atomicNumber == el_nucs.first->atomicNumber) );
      if( get<2>(nucs) )
        return true;
    }//for( const auto &nucs : el_nucs.second )
  }//for( const auto &el_nucs : shield.m_nuclideFractions_ )

  return false;
}//bool hasVariableMassFraction( const Material *material ) const

  
vector<const Material *> ShieldingSourceChi2Fcn::materialsFittingMassFracsFor() const
{
  vector<const Material *> answer;
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &ms : m_initial_shieldings )
  {
    bool isFittingAny = false;
    for( const auto &el_nucs : ms.m_nuclideFractions_ )
    {
      for( const auto &nuc_info : el_nucs.second )
        isFittingAny |= get<2>(nuc_info);
    }
    
    if( isFittingAny )
      answer.push_back( ms.m_material.get() );
  }//for( ShieldLayerInfo &ms : m_initial_shieldings )
  
  return answer;
}//vector<const Material *material> materialsFittingMassFracsFor() const
  

vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::selfAttenuatingNuclides( const size_t material_index ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );

  const ShieldingSourceFitCalc::ShieldingInfo &ms = m_initial_shieldings[material_index];
  
  if( !ms.m_material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::selfAttenuatingNuclides():"
                        " material index points to generic shielding" );
  
  vector<const SandiaDecay::Nuclide *> answer;
  for( const auto &el_nucs : ms.m_nuclideFractions_ )
  {
    for( const auto &nuc_info : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
      if( nuc ) //"other" fraction of element is nullptr
        answer.push_back( nuc );
    }//for( const auto &nuc_info : el_nucs.second )
  }//for( const auto &el_nucs : ms.m_nuclideFractions_ )
  
  return answer;
}//vector<const SandiaDecay::Nuclide *> selfAttenuatingNuclides( const size_t material_index ) const
  
  
vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::traceNuclidesForMaterial( const size_t material_index ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );

  const ShieldingSourceFitCalc::ShieldingInfo &ms = m_initial_shieldings[material_index];
  
  if( !ms.m_material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::traceNuclidesForMaterial():"
                        " material index points to generic shielding" );
  
  vector<const SandiaDecay::Nuclide *> answer;
  for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : ms.m_traceSources )
  {
    assert( trace.m_nuclide );
    if( trace.m_nuclide )
      answer.push_back( trace.m_nuclide );
  }
  
  return answer;
}//std::vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::traceNuclidesForMaterial( const size_t material_index ) const


std::map<const SandiaDecay::Element *,std::vector<const SandiaDecay::Nuclide *>>
        ShieldingSourceChi2Fcn::nuclideFittingMassFracFor( const size_t material_index ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &ms = m_initial_shieldings[material_index];
  
  if( !ms.m_material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor():"
                        " material index points to generic shielding" );

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  map<const SandiaDecay::Element *,vector<const SandiaDecay::Nuclide *>> answer;
  for( const auto &el_nucs : ms.m_nuclideFractions_ )
  {
    for( const auto &nuc_info : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
      const bool fit = get<2>(nuc_info);
      
      if( nuc && fit ) //"other" fraction of element is nullptr
        answer[el_nucs.first].push_back( nuc );
    }//for( const auto &nuc_info : el_nucs.second )
  }//for( const auto &el_nucs : ms.m_nuclideFractions_ )
  
  return answer;
}//nuclideFittingMassFracFor(...)

  
vector<const SandiaDecay::Element *> ShieldingSourceChi2Fcn::elementsFittingOtherFracFor(
                                                              const size_t material_index ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &ms = m_initial_shieldings[material_index];
  
  if( !ms.m_material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor():"
                        " material index points to generic shielding" );

  vector<const SandiaDecay::Element *> answer;
  for( const auto &el_nucs : ms.m_nuclideFractions_ )
  {
    const SandiaDecay::Element * const el = el_nucs.first;
    assert( el );
    
    for( const auto &nuc_info : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
      const bool fit = get<2>(nuc_info);
      
      if( el && !nuc && fit ) //"other" fraction of element is nullptr
      {
        assert( std::count(begin(answer), end(answer), el) == 0 );
        answer.push_back( el );
#ifdef NDEBUG
        break; //on debug builds we'll keep looping to make sure only one entry with nullptr for nuclide
#endif
      }
    }//for( const auto &nuc_info : el_nucs.second )
  }//for( const auto &el_nucs : ms.m_nuclideFractions_ )
  
  return answer;
}//elementsFittingOtherFracFor(...)
  
  
map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,double,bool>>> 
                  ShieldingSourceChi2Fcn::selfAttenSrcInfo( const size_t material_index,
                                                          const vector<double> &pars,
                                                          const vector<double> &error ) const
{
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::selfAttenSrcInfo: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &ms = m_initial_shieldings[material_index];
  
  map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,double,bool>>> answer;
  
  assert( ms.m_material );
  if( !ms.m_material )
    return answer;
  
  for( const auto &el_nucs : ms.m_nuclideFractions_ )
  {
    const SandiaDecay::Element * const el = el_nucs.first;
    assert( el );
    
    for( const auto &nuc_inf : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const nuc = get<0>(nuc_inf);
      double frac = get<1>(nuc_inf), uncert = 0.0;
      const bool fit_frac = get<2>(nuc_inf);
      if( fit_frac )
        massFractionOfElement( frac, uncert, material_index, nuc, el, pars, error );
      answer[el].emplace_back( nuc, frac, uncert, fit_frac );
    }//for( const auto &nuc_inf : el_nucs.second )
  }//for( const auto &el_nucs : ms.m_nuclideFractions_ )
  
  return answer;
}//selfAttenSrcInfo(...)
  
  
double ShieldingSourceChi2Fcn::massFractionOfElement( const size_t material_index,
                      const SandiaDecay::Nuclide *nuc,
                      const std::vector<double> &pars ) const
{
  assert( nuc );
  if( !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): nullptr nuc" );
  double massfrac = 0.0, uncert = 0.0;
  massFractionOfElement( massfrac, uncert, material_index, nuc, nullptr, pars, vector<double>() );
  return massfrac;
}
  
  
void ShieldingSourceChi2Fcn::massFractionOfElement( double &massFrac, double &uncert,
                                          const size_t material_index, 
                                          const SandiaDecay::Nuclide *nuc,
                                          const SandiaDecay::Element *el,
                                          const vector<double> &pars,
                                          const vector<double> &errors ) const
{
  massFrac = uncert = 0.0;
  
  if( material_index >= m_initial_shieldings.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::massFraction: invalid material index" );
  
  const ShieldingSourceFitCalc::ShieldingInfo &shielding = m_initial_shieldings[material_index];
  
  assert( nuc || el );
  assert( !nuc || !el || (nuc->atomicNumber == el->atomicNumber) );
  
  if( !shielding.m_material || (!el && !nuc) )
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): invalid input" );
  
  if( nuc && el && (nuc->atomicNumber != el->atomicNumber) )
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): invalid element for nuclide" );
  
  //We may not be fitting for this mass fraction; lets check this, and if we arent, return its
  //  initial mass fraction
  const short int atomic_num = nuc ? nuc->atomicNumber : el->atomicNumber;
  for( const auto &el_nucs : shielding.m_nuclideFractions_ )
  {
    const auto * const this_el = el_nucs.first;
    assert( this_el );
    if( !this_el || (this_el->atomicNumber != atomic_num) )
      continue;
    
    for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const this_nuc = get<0>(nuc_info);
      if( this_nuc == nuc )
      {
        const bool fit_frac = get<2>(nuc_info);
        if( !fit_frac )
        {
          uncert = 0.0;
          massFrac = get<1>(nuc_info);
          return;
        }//
        
        break; // We found the nuclide, no need to continue
      }//if( we found the nuclide of interest )
    }//for( loop over nuclide information for this Element )
  }//for( loop over to check if we are actually fitting this mass fraction )
  
  
  // We need to know how many mass fraction parameters there are for shieldings before the
  //  shielding we are interested in and how many parameters in the shielding of interest
  //  come before the element of the nuclide we care about, in the shielding we care about.
  size_t num_pre_el_coefs = 0;
  
  //Get how many mass fraction coefficients are for materials that come before this one
  for( size_t pre_mat_index = 0; pre_mat_index < material_index; ++pre_mat_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &pre_shield = m_initial_shieldings[pre_mat_index];
    assert( !!pre_shield.m_material || pre_shield.m_nuclideFractions_.empty() );
    if( !pre_shield.m_material ) //skip if generic material
      continue;
    
    for( const auto &el_nucs : pre_shield.m_nuclideFractions_ )
    {
      size_t num_fit_this_el = 0;
      for( const auto &nuc : el_nucs.second )
        num_fit_this_el += get<2>(nuc);
      
      assert( num_fit_this_el != 1 );
      if( num_fit_this_el > 1 )
        num_pre_el_coefs += (num_fit_this_el - 1);
    }//for( loop over nuclide fractions )
  }//for( loop over materials before current one to count number mass fraction parameters )
  
  
  const map<const SandiaDecay::Element *,vector<tuple<const SandiaDecay::Nuclide *,double,bool>>>
  &self_atten_nucs = shielding.m_nuclideFractions_;
  
  const tuple<const SandiaDecay::Nuclide *,double,bool> *initial_nuc_info = nullptr;
  const vector<tuple<const SandiaDecay::Nuclide *,double,bool>> *initial_el_info = nullptr;
  
  // We also need to know what index the nuclide of interest is in, in its element.
  size_t nuc_index_in_el = 0, num_variable_nucs_in_el = 0;
  
  for( const auto &el_nucs : self_atten_nucs )
  {
    const SandiaDecay::Element * const this_el = el_nucs.first;
    assert( this_el );
    if( !this_el )
      throw std::logic_error( "massFraction: nullptr element" );
    
    if( (nuc && (this_el->atomicNumber == nuc->atomicNumber)) 
       || (el && (el->atomicNumber == this_el->atomicNumber)) )
    {
      assert( !el || (el == this_el) );
      if( !el )
        el = this_el;
      initial_el_info = &el_nucs.second;
    }else
    {
      size_t num_fit_fracs = 0;
      for( const auto &nuc_info : el_nucs.second )
        num_fit_fracs += get<2>(nuc_info);
      
      assert( num_fit_fracs != 1 );
      
      if( num_fit_fracs > 0 )
        num_pre_el_coefs += (num_fit_fracs - 1);
      
      continue;
    }//if( this is the correct
    
    for( const auto &nuc_info : el_nucs.second )
    {
      num_variable_nucs_in_el += get<2>(nuc_info);
      
      if( !initial_nuc_info )
      {
        if( get<0>(nuc_info) == nuc )
          initial_nuc_info = &nuc_info;
        else
          nuc_index_in_el += get<2>(nuc_info);
      }//if( !initial_nuc_info )
    }//for( const auto &nuc_info : el_nucs.second )
    
    assert( el );
    assert( nuc_index_in_el < el_nucs.second.size() );
    assert( !nuc || (this_el->atomicNumber == nuc->atomicNumber) );
      
    if( (nuc && (this_el->atomicNumber == nuc->atomicNumber))
       || (el && (this_el->atomicNumber == el->atomicNumber)) )
    {
      assert( initial_el_info );
      assert( !nuc || initial_nuc_info );
      
      break;
    }
  }//for( const auto &el_nucs : self_atten_nucs )
  
  assert( initial_el_info );
  assert( initial_nuc_info || !nuc );
  
  if( !initial_el_info )
    throw logic_error( "Failed to find element in shielding???" );
  
  assert( num_variable_nucs_in_el > 1 );
  if( num_variable_nucs_in_el <= 1 )
    throw logic_error( "Found invalid number of variable nucs for an element" );
  
  // If we want the other non-source fraction for this element, it _may_ not be explicitly
  //  included in the material, in that case we know we arent fitting for it, so we'll just
  //  return 1.0 minus all the other source components.
  if( !nuc && !initial_nuc_info )
  {
    uncert = 0.0;
    double other_frac = 0.0;
    for( const auto &nuc_info : *initial_el_info )
    {
      assert( get<0>(nuc_info) );
      if( get<0>(nuc_info) )
        other_frac += get<1>(nuc_info);
    }
      
    massFrac = std::max( 1.0 - other_frac, 0.0 );
    return;
  }//if( we wanted other nuc info, but that wasnt explicitly included )
  
  if( !initial_nuc_info || !initial_el_info )
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): "
                         + nuc->symbol + " was not a self-attenuating nuclide in "
                         + shielding.m_material->name );
  
  // If we arent fitting mass-fraction, we can return here.
  if( !get<2>(*initial_nuc_info) )
  {
    massFrac = get<1>(*initial_nuc_info);
    uncert = 0.0;
    return;
  }//if( we are not fitting this mass fraction )
  
  double totalfrac = 0.0;//The total mass fraction of elements being fit, for this nuclide
  for( const tuple<const SandiaDecay::Nuclide *,double,bool> &i : *initial_el_info )
  {
    const double frac = get<1>(i);
    const bool fit = get<2>(i);
    if( fit )
      totalfrac += frac;
  }//for( loop over nuclides of this element )

  
  size_t matmassfracstart = 0;
  matmassfracstart += 2 * m_nuclides.size();
  matmassfracstart += 3 * m_initial_shieldings.size();
  matmassfracstart += num_pre_el_coefs;
  
  
  const size_t numfraccoefs = num_variable_nucs_in_el - 1;
  const size_t thismatmassfrac = matmassfracstart + nuc_index_in_el;
  
  double frac = 1.0;
  for( size_t index = matmassfracstart; index < thismatmassfrac; ++index )
    frac *= (1.0 - pars.at(index));
  double prefrac = frac;
  
  if( nuc_index_in_el != numfraccoefs )
    frac *= pars.at(thismatmassfrac);
  
  massFrac = totalfrac * frac;
  
  if( IsNan(massFrac) || IsInf(massFrac) )
  {
    cerr << "Got invalid mass frac:" << endl;
    cerr << "\ttotalfrac=" << totalfrac << endl;
    cerr << "\tprefrac=" << prefrac << endl;
    if( nuc_index_in_el != numfraccoefs )
      cerr << "\tpars.at(thismatmassfrac)=" << pars.at(thismatmassfrac) << endl;
//    cerr << "\tPars: ";
//    for( size_t index = matmassfracstart; index < thismatmassfrac; ++index )
//      cerr << "{" << index << "," << pars.at(index) << "},";
//    cerr << endl;
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(...): invalid parameters" );
  }
  
  if( !errors.empty() )
  {
    if( errors.size() != pars.size() )
      throw runtime_error( "ShieldingSourceChi2Fcn::massFraction():"
                           " invalid error parameter vector size" );
    double fracuncert = 0.0;
    if( nuc_index_in_el != numfraccoefs )
      fracuncert = errors.at(thismatmassfrac) / pars.at(thismatmassfrac);
    else
      fracuncert = errors.at(thismatmassfrac-1) / pars.at(thismatmassfrac-1);
    uncert = fracuncert * massFrac;
  }//if( !errors.empty() )
}//double massFraction(...) const

  
double ShieldingSourceChi2Fcn::massFractionOfElementUncertainty( const size_t material_index,
                            const SandiaDecay::Nuclide *nuc,
                            const std::vector<double> &pars,
                            const std::vector<double> &error ) const
{
  const SandiaDecay::Element * const el = nullptr; //`massFraction(...)` should pick this up
  double massfrac = 0.0, uncert = 0.0;
  massFractionOfElement( massfrac, uncert, material_index, nuc, el, pars, error );
  return uncert;
}//massFractionUncert(...)
  

size_t ShieldingSourceChi2Fcn::numExpectedFitParameters() const
{
  size_t npar = 2 * m_nuclides.size();
  
  npar += 3 * m_initial_shieldings.size();
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    for( const auto &el_nucs : info.m_nuclideFractions_ )
    {
      size_t num_nucs_this_el = 0;
      for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucs.second )
        num_nucs_this_el += get<2>(nuc_info);
      
      assert( num_nucs_this_el != 1 );
      
      if( num_nucs_this_el > 1 )
        npar += (num_nucs_this_el - 1);
    }//for( const auto &el_nucs : info.m_nuclideFractions_ )
  }//for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  
  return npar;
}//int numExpectedFitParameters() const


bool ShieldingSourceChi2Fcn::isVolumetricSource( const SandiaDecay::Nuclide *nuc ) const
{
  return isSelfAttenSource(nuc) || isTraceSource(nuc);
}


bool ShieldingSourceChi2Fcn::isSelfAttenSource( const SandiaDecay::Nuclide *nuclide ) const
{
  if( !nuclide )
    return false;

  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    for( const auto &el_nucs : info.m_nuclideFractions_ )
    {
      const SandiaDecay::Element * const el = el_nucs.first;
      assert( el );
      if( !el || (el->atomicNumber != nuclide->atomicNumber) )
        continue;
      
      for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucs.second )
      {
        const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
        if( nuc == nuclide )
          return true;
      }//for( loop over self-atten nucs in this element )
    }//for( loop over self-atten elements )
  }//for( loop over shieldings )
  
  return false;
}//bool isSelfAttenSource(...) const;


bool ShieldingSourceChi2Fcn::isTraceSource( const SandiaDecay::Nuclide *nuclide ) const
{
  if( !nuclide )
    return false;
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      if( trace.m_nuclide == nuclide )
        return true;
    }//for( loop over trace sources )
  }//for( loop over shieldings )
   
  return false;
}//bool isTraceSource(...) const;


TraceActivityType ShieldingSourceChi2Fcn::traceSourceActivityType(
                                                          const SandiaDecay::Nuclide *nuc ) const
{
  if( !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: null nuclide" );
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      if( trace.m_nuclide == nuc )
        return trace.m_type;
    }//for( loop over trace sources )
  }//for( loop over shieldings )
  
  throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: " + nuc->symbol
                       + " not a trace source" );
  
  return TraceActivityType::NumTraceActivityType;
}//traceSourceActivityType(...)


double ShieldingSourceChi2Fcn::relaxationLength( const SandiaDecay::Nuclide *nuc ) const
{
  if( !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: null nuclide" );
  
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : m_initial_shieldings )
  {
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      if( trace.m_nuclide == nuc )
      {
        if( trace.m_type != TraceActivityType::ExponentialDistribution )
          throw runtime_error( "ShieldingSourceChi2Fcn::relaxationLength: " + nuc->symbol
                               + " is not an exponential distribution trace source." );
        return trace.m_relaxationDistance;
      }
    }//for( loop over trace sources )
  }//for( loop over shieldings )
  
  throw runtime_error( "ShieldingSourceChi2Fcn::relaxationLength: " + nuc->symbol
                      + " not a trace source" );
  
  return -1.0;
}//double relaxationLength( const SandiaDecay::Nuclide *nuc ) const;


double ShieldingSourceChi2Fcn::volumeOfMaterial( const size_t matn, const vector<double> &params ) const
{
  if( matn >= m_initial_shieldings.size() )
    throw runtime_error( "volumeOfMaterial: invalid material index" );
  
  if( !m_initial_shieldings[matn].m_material )
    throw runtime_error( "volumeOfMaterial: cant be called for generic material" );
  
  double inner_dim_1 = 0.0, inner_dim_2 = 0.0, inner_dim_3 = 0.0;
  double outer_dim_1 = 0.0, outer_dim_2 = 0.0, outer_dim_3 = 0.0;
  
  for( int index = 0; index <= matn; ++index )
  {
    if( !m_initial_shieldings[index].m_material )  //if a generic shielding
      continue;
    
    inner_dim_1 = outer_dim_1;
    inner_dim_2 = outer_dim_2;
    inner_dim_3 = outer_dim_3;
    
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        outer_dim_1 += sphericalThickness( index, params );
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        outer_dim_1 += cylindricalRadiusThickness( index, params );
        outer_dim_2 += cylindricalLengthThickness( index, params );
        break;
        
      case GeometryType::Rectangular:
        outer_dim_1 += rectangularWidthThickness( index, params );
        outer_dim_2 += rectangularHeightThickness( index, params );
        outer_dim_3 += rectangularDepthThickness( index, params );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
  }//for( int index = 0; index < matn; ++index )
  
  switch( m_geometry )
  {
    case GeometryType::Spherical:
    {
      return (4.0/3.0)*PhysicalUnits::pi*( pow(outer_dim_1,3.0) - pow(inner_dim_1,3.0));
    }
      
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      const double inner_vol = PhysicalUnits::pi * inner_dim_1 * inner_dim_1 * 2.0 * inner_dim_2;
      const double outer_vol = PhysicalUnits::pi * outer_dim_1 * outer_dim_1 * 2.0 * outer_dim_2;
      
      return outer_vol - inner_vol;
    }
      
    case GeometryType::Rectangular:
    {
      const double inner_vol = 8.0 * inner_dim_1 * inner_dim_2 * inner_dim_3;
      const double outer_vol = 8.0 * outer_dim_1 * outer_dim_2 * outer_dim_3;
      
      return outer_vol - inner_vol;
    }
      
    case GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch( m_geometry )
  
  assert( 0 );
  return 0.0;
}//double volumeOfMaterial(...)


double ShieldingSourceChi2Fcn::volumeUncertaintyOfMaterial( const int matn,
                                  const vector<double> &params, const vector<double> &errors ) const
{
  // We will take into account the uncertainties of the inner layers to our current shell, the
  //  uncertainty of the thickness, and the uncertainty of the mass fractions.
  if( (matn < 0) || (matn >= m_initial_shieldings.size()) )
    throw runtime_error( "volumeUncertaintyOfMaterial: invalid material index" );
  
  if( !m_initial_shieldings[matn].m_material )
    throw runtime_error( "volumeUncertaintyOfMaterial: cant be called for generic material" );
  
  // TODO: clean this up to be a little neater/shorter, after testing
  double dim_1 = 0.0, dim_2 = 0.0, dim_3 = 0.0;
  double dim_1_uncert_sq = 0.0, dim_2_uncert_sq = 0.0, dim_3_uncert_sq = 0.0;
  
  for( int mat_index = 0; mat_index <= matn; ++mat_index )
  {
    if( !m_initial_shieldings[mat_index].m_material )  //if a generic shielding
      continue;
    
    switch( m_geometry )
    {
      case GeometryType::Spherical:
      {
        const double thick = sphericalThickness( mat_index, params );
        const double thickUncert = sphericalThickness( mat_index, errors );
        
        if( mat_index == matn )
        {
          const double outerRad = dim_1 + thick;
          
          const double volumeUncertDueToInnerRad = 4.0 * PhysicalUnits::pi * pow(dim_1,2.0) * sqrt(dim_1_uncert_sq);
          // We take the uncertainty in volume due to the thickness uncertainty to be independent of
          //  the inner-radius uncertainty - which is approximately, probably, about right to kinda
          //  fairly cover all the uncertainties - maybe
          const double volumeUncertDueToThickness = 4.0 * PhysicalUnits::pi * pow(outerRad,2.0) * thickUncert;
          
          // Take the volume uncert to be the sqrt of sum of inner and outer volumes
          const double volUncertSquared = pow(volumeUncertDueToInnerRad,2.0) + pow(volumeUncertDueToThickness,2.0);
          
          return sqrt( volUncertSquared );
        }//if( mat_index == matn )
        
        // For the layer of shielding we are concerned with, we will assume the uncertainty on the
        //  inner radius is the squared sum of all inner layer uncertainties - this is actually a bit
        //  over the top - these thickness uncertainties are likely very correlated - so much so we
        //  could probably use just the thickness uncertainty on just the previous layer, but we'll be
        //  conservative, or something.
        dim_1_uncert_sq += thickUncert*thickUncert;
        dim_1 += thick;
        
        break;
      }//case GeometryType::Spherical:
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
      {
        const double dr = cylindricalRadiusThickness( mat_index, params );
        const double dr_uncert = cylindricalRadiusThickness( mat_index, errors );
        
        const double dz = cylindricalLengthThickness( mat_index, params );
        const double dz_uncert = cylindricalLengthThickness( mat_index, errors );
        
        if( mat_index == matn )
        {
          const double outerRad = dim_1 + dr;
          const double outerZ = dim_2 + dz;
          
          const double volUncertInnerRad = 2.0 * PhysicalUnits::pi * dim_1 * 2.0 * dim_2 * sqrt(dim_1_uncert_sq);
          const double volUncertOuterRad = 2.0 * PhysicalUnits::pi * outerRad *  2.0 * outerZ * dr_uncert;
          
          const double volUncertInnerZ = PhysicalUnits::pi * dim_1 * dim_1 * 2.0 * sqrt(dim_2_uncert_sq);
          const double volUncertOuterZ = PhysicalUnits::pi * dim_1 * dim_1 * 2.0 * dz_uncert;
        
          const double volUncertSquared = pow(volUncertInnerRad,2.0) + pow(volUncertOuterRad,2.0)
                                          + pow(volUncertInnerZ,2.0) + pow(volUncertOuterZ,2.0);
          
          return sqrt( volUncertSquared );
        }//if( mat_index == matn )
        
        dim_1_uncert_sq += dr_uncert * dr_uncert;
        dim_2_uncert_sq += dz_uncert * dz_uncert;
        
        dim_1 += dr;
        dim_2 += dz;
        
        break;
      }//case GeometryType::CylinderEndOn or CylinderSideOn
        
      case GeometryType::Rectangular:
      {
        const double dw = rectangularWidthThickness( mat_index, params );
        const double dw_uncert = 2.0 * rectangularWidthThickness( mat_index, errors );
        
        const double dh = rectangularHeightThickness( mat_index, params );
        const double dh_uncert = 2.0 * rectangularHeightThickness( mat_index, errors );
        
        const double dd = rectangularDepthThickness( mat_index, params );
        const double dd_uncert = 2.0 * rectangularDepthThickness( mat_index, errors );
        
        if( mat_index == matn )
        {
          const double inner_w = 2.0 * dim_1;
          const double inner_h = 2.0 * dim_2;
          const double inner_d = 2.0 * dim_3;
          
          const double outer_w = 2.0 * (dim_1 + dw);
          const double outer_h = 2.0 * (dim_2 + dh);
          const double outer_d = 2.0 * (dim_3 + dd);
          
          
          const double vol_uncert_inner_1 = inner_h * inner_d * sqrt(dim_1_uncert_sq);
          const double vol_uncert_inner_2 = inner_w * inner_d * sqrt(dim_2_uncert_sq);
          const double vol_uncert_inner_3 = inner_w * inner_h * sqrt(dim_3_uncert_sq);
          
          const double vol_uncert_outer_1 = outer_h * outer_d * dw_uncert;
          const double vol_uncert_outer_2 = outer_w * outer_d * dh_uncert;
          const double vol_uncert_outer_3 = outer_w * outer_h * dd_uncert;
          
          const double inner_vol_uncert_sq = vol_uncert_inner_1*vol_uncert_inner_1 + vol_uncert_inner_2*vol_uncert_inner_2 + vol_uncert_inner_3*vol_uncert_inner_3;
          const double outer_vol_uncert_sq = vol_uncert_outer_1*vol_uncert_outer_1 + vol_uncert_outer_2*vol_uncert_outer_2 + vol_uncert_outer_3*vol_uncert_outer_3;
          
          return sqrt( inner_vol_uncert_sq + outer_vol_uncert_sq);
        }//if( mat_index == matn )
        
        dim_1_uncert_sq += dw_uncert * dw_uncert;
        dim_2_uncert_sq += dh_uncert * dh_uncert;
        dim_3_uncert_sq += dd_uncert * dd_uncert;
        
        dim_1 += dw;
        dim_2 += dh;
        dim_3 += dd;
        
        break;
      }
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
  }//for( int index = 0; index < matn; ++index )
  
  assert( 0 );
  return 0.0;
}//volumeUncertaintyOfMaterial(...)


double ShieldingSourceChi2Fcn::activityOfSelfAttenSource(
                                       const SandiaDecay::Nuclide *nuclide,
                                       const std::vector<double> &params ) const
{
  assert( nuclide );
  if( !nuclide )
    throw runtime_error( "ShieldingSourceChi2Fcn::activityOfSelfAttenSource: called with null nuclide" );
  
  assert( isSelfAttenSource(nuclide) );
  
  bool foundSrc = false;
  double activity = 0.0;

  const size_t num_mataterials = m_initial_shieldings.size();
  for( size_t material_index = 0; material_index < num_mataterials; ++material_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
    
    const shared_ptr<const Material> &mat = shield.m_material;
    if( !mat )  //if a generic shielding
      continue;
    
    // Determine if nuclide is a self-attenuating source of this shielding
    const tuple<const SandiaDecay::Nuclide *,double,bool> *self_att_pos = nullptr;
    for( const auto &el_nucs : shield.m_nuclideFractions_ )
    {
      assert( el_nucs.first );
      if( el_nucs.first->atomicNumber != nuclide->atomicNumber )
        continue;
      
      for( const auto &nuc_info : el_nucs.second )
      {
        if( get<0>(nuc_info) == nuclide )
        {
          self_att_pos = &(nuc_info);
          break;
        }//if( get<0>(nuc_info) == nuclide )
      }//for( const auto &nuc_info : el_nucs.second )
    }//for( const auto &el_nucs : shield.m_nuclideFractions_ )
    
    // Check that nuclide is a self-atten source of this shielding, if not lets keep looking.
    if( !self_att_pos )
      continue;
    
    const SandiaDecay::Nuclide *src = get<0>(*self_att_pos);
    const double initial_frac = get<1>(*self_att_pos);
    const bool fit = get<2>(*self_att_pos);
    
    foundSrc = true;
    
    double massFrac = 0.0;
    
    if( fit )
      massFrac = ShieldingSourceChi2Fcn::massFractionOfElement( material_index, nuclide, params );
    else
      massFrac = initial_frac;
    
    const double vol = volumeOfMaterial( material_index, params );
    const double mass_grams = massFrac * mat->density * vol / PhysicalUnits::gram;
    const double activity_per_gram = nuclide->activityPerGram();
    activity += mass_grams * activity_per_gram;
    
    //        cerr << "activityOfSelfAttenSource: " << nuclide->symbol << " (material_index=" << material_index
    //         << "), mat=" << mat->name
    //        << ", thick=" << thick/PhysicalUnits::cm
    //            << " cm, inner_rad=" << (radius-thick)/PhysicalUnits::cm
    //        << " cm, vol=" << vol/(PhysicalUnits::cm*PhysicalUnits::cm*PhysicalUnits::cm)
    //        << " cm3, density=" << mat->density*(PhysicalUnits::cm*PhysicalUnits::cm*PhysicalUnits::cm)/PhysicalUnits::g
    //        << " g/cm3, mass=" << mass_grams << " g, activity_per_gram=" << activity_per_gram/PhysicalUnits::ci
    //        << " ci, mass_grams * activity_per_gram=" << mass_grams * activity_per_gram/PhysicalUnits::ci << " ci" << endl;
    
    // We could probably break the loop here, but I guess JIC there are multiple objects using this
    //  nuclide as a self-attenuating source we will keep going.
  }//for( int material_index = 0; material_index < num_mataterials; ++material_index )
  
  if( !foundSrc )
    throw runtime_error( "ShieldingSourceChi2Fcn::activityOfSelfAttenSource: "
                        + nuclide->symbol + " is not a self-attenuating source" );
  
//  cerr << "activityOfSelfAttenSource: " << nuclide->symbol << " activity=" << activity/PhysicalUnits::ci << endl;
  
  return activity;
}//double activityOfSelfAttenSource(...) const;


double ShieldingSourceChi2Fcn::totalActivity( const SandiaDecay::Nuclide *nuclide,
                                  const std::vector<double> &params ) const
{
  assert( nuclide );
  if( !nuclide )
    throw runtime_error( "ShieldingSourceChi2Fcn::totalActivity: called with null nuclide" );
  
  const bool isTrace = isTraceSource(nuclide);
  const bool isSelfAtten = !isTrace && isSelfAttenSource(nuclide);
  
  if( isSelfAtten )
    return activityOfSelfAttenSource( nuclide, params );
  
  const size_t ind = nuclideIndex( nuclide );
  const double activity = params[2*ind] * sm_activityUnits;

  if( !isTrace )
    return activity;
  
  bool foundSrc = false;
  
  for( size_t matn = 0; matn < m_initial_shieldings.size(); ++matn )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[matn];
    
    const shared_ptr<const Material> &mat = shield.m_material;
    if( !mat )  //if a generic shielding
      continue;
    
    const double vol = volumeOfMaterial( matn, params );
    
    // Determine if nuclide is a trace source of this shielding
    TraceActivityType type = TraceActivityType::NumTraceActivityType;
    const std::vector<ShieldingSourceFitCalc::TraceSourceInfo> &trace_srcs = shield.m_traceSources;
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : trace_srcs )
    {
      if( trace.m_nuclide == nuclide )
      {
        type = trace.m_type;
        break;
      }
    }//for( const auto &trace : trace_srcs )
    
    if( type == TraceActivityType::NumTraceActivityType )
      continue;  //we didnt find it
    
    foundSrc = true;
    
    switch( type )
    {
      case TraceActivityType::TotalActivity:
        return activity;
        
      case TraceActivityType::ActivityPerCm3:
        return activity * vol / PhysicalUnits::cm3;
        break;
        
      case TraceActivityType::ActivityPerGram:
        return activity * mat->density * vol / PhysicalUnits::gram;
        break;
        
      case TraceActivityType::ExponentialDistribution:
      {
        switch( m_geometry )
        {
          case GeometryType::Spherical:
          {
            const double r = sphericalThickness(matn, params);
            return activity * 4.0 * PhysicalUnits::pi * r * r / PhysicalUnits::m2;
          }
            
          case GeometryType::CylinderEndOn:
          {
            const double r = cylindricalRadiusThickness(matn, params);
            return activity * PhysicalUnits::pi * r * r / PhysicalUnits::m2;
          }
            
          case GeometryType::CylinderSideOn:
          {
            const double r = cylindricalRadiusThickness(matn, params);
            const double z = cylindricalLengthThickness(matn, params);
            return activity * PhysicalUnits::pi * 2.0 * r * z / PhysicalUnits::m2;
          }
            
          case GeometryType::Rectangular:
          {
            const double w = 2.0 * rectangularWidthThickness(matn, params);
            const double h = 2.0 * rectangularHeightThickness(matn, params);
            return activity * w * h / PhysicalUnits::m2;
          }
            
          case GeometryType::NumGeometryType:
            assert( 0 );
            break;
        }//switch( m_geometry )
        
        break;
      }//case TraceActivityType::ExponentialDistribution:
        
      case TraceActivityType::NumTraceActivityType:
        assert( 0 );
        break;
    }//switch( type )
  }//for( int matn = 0; matn < nmat; ++matn )
  
  throw runtime_error( "ShieldingSourceChi2Fcn::totalActivity: "
                        + nuclide->symbol + " is not a self-attenuating source" );
  
  return activity;
}//double totalActivity(...)



double ShieldingSourceChi2Fcn::totalActivityUncertainty( const SandiaDecay::Nuclide *nuclide,
                                const std::vector<double> &params,
                                const std::vector<double> &errors ) const
{
  assert( nuclide );
  if( !nuclide )
    throw runtime_error( "ShieldingSourceChi2Fcn::totalActivityUncertainty: called with null nuclide" );
  
  assert( params.size() == errors.size() );
  
  const bool isTrace = isTraceSource(nuclide);
  
  double uncert = activityUncertainty( nuclide, params, errors );
  
  if( !isTrace )
    return uncert;
  
  const double display_activity = activity( nuclide, params );
  const double total_activity = totalActivity( nuclide, params );
  
  if( (display_activity > FLT_EPSILON) && (total_activity > FLT_EPSILON) )
    uncert *= (total_activity / display_activity);
  
  return uncert;
}//double totalActivityUncertainty(...)



double ShieldingSourceChi2Fcn::activity( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
  {
    cerr << "params.size()=" << params.size()
         << ", numExpectedFitParameters()=" << numExpectedFitParameters()
         << endl;
    throw runtime_error( "ShieldingSourceChi2Fcn::activity(...) invalid params size" );
  }//if( params.size() != numExpectedFitParameters() )

  if( isSelfAttenSource(nuclide) )
    return activityOfSelfAttenSource( nuclide, params );
  
  const size_t ind = nuclideIndex( nuclide );  
  return params[2*ind] * sm_activityUnits;
}//double activity(...)



double ShieldingSourceChi2Fcn::activityUncertainty( const SandiaDecay::Nuclide *nuclide,
                                                       const std::vector<double> &params,
                                                        const std::vector<double> &errors ) const
{
  if( params.size() != numExpectedFitParameters() )
  {
    cerr << "params.size()=" << params.size()
    << ", numExpectedFitParameters()=" << numExpectedFitParameters()
    << endl;
    throw runtime_error( "ShieldingSourceChi2Fcn::activityUncertainty(...) invalid params size" );
  }//if( params.size() != numExpectedFitParameters() )
  
  if( params.size() != errors.size() )
    throw runtime_error( "ShieldingSourceChi2Fcn::activityUncertainty: size(params) != size(errors)" );
  
  // If it is a simple point source, we can get the uncertainty by just passing in the
  //  parameter uncertainties and calculating activity from those.  However, if its a
  //  self-attenuating source, or trace source we have to do a little more propagation and such.
  const bool isTrace = isTraceSource(nuclide);
  const bool isSelfAtten = !isTrace && isSelfAttenSource(nuclide);
  
  const double niaveActUncert = activity( nuclide, errors );
  
  if( !isTrace && !isSelfAtten )
  {
    return niaveActUncert;
  }//
    
  
  // We'll track activity just as a runtime check we are computing things consistent with
  //  #ShieldingSourceChi2Fcn::totalActivity.
  double activity = 0.0;
  
  // We will take into account the uncertainties of the inner layers to our current shell, the
  //  uncertainty of the thickness, and the uncertainty of the mass fractions.
  double activityUncertSquared = 0.0;
  
  const size_t num_materials = m_initial_shieldings.size();
  for( size_t material_index = 0; material_index < num_materials; ++material_index )
  {
    if( isGenericMaterial(material_index) )
      continue;
    
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
    const shared_ptr<const Material> &mat = shield.m_material;
    assert( mat );
    
    const double vol = volumeOfMaterial(material_index, params);
    const double volUncert = volumeUncertaintyOfMaterial(static_cast<int>(material_index), params, errors);
  
    // Add in uncertainty contributions if this is a self attenuating source
    for( const auto el_nucs : shield.m_nuclideFractions_ )
    {
      const SandiaDecay::Element * const el = el_nucs.first;
      assert( el && el_nucs.second.size() );
      
      for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucs.second )
      {
        const SandiaDecay::Nuclide *src = get<0>(nuc_info);
        const double initial_frac = get<1>(nuc_info);
        const bool fit = get<2>(nuc_info);
        
        if( src != nuclide )
          continue;
        
        double massFrac = 0.0, massFracUncert = 0.0;
        if( fit )
          massFractionOfElement( massFrac, massFracUncert, material_index, src, el, params, errors );
        else
          massFrac = initial_frac;
        
        
        const double iso_density = massFrac * mat->density / PhysicalUnits::gram;
        const double activity_per_gram = nuclide->activityPerGram();
        const double mass_grams = iso_density * vol;
        const double massUncertaintySquared_grams = iso_density * volUncert * volUncert;
        const double thisActivity = mass_grams * activity_per_gram;
        const double thisActivityUncertaintySquared = massUncertaintySquared_grams * activity_per_gram;
        
        //cout << src->symbol << ": sqrt(thisActivityUncertaintySquared)=" << sqrt(thisActivityUncertaintySquared) << endl
        //<< "\tvolumeUncertDueToInnerRad=" << volumeUncertDueToInnerRad
        //<< " (radUncert=" << sqrt(radiusUncertSquared)/PhysicalUnits::cm << " cm)" << endl
        //<< "\tvolumeUncertDueToThickness=" << volumeUncertDueToThickness
        //<< " (thickUncert=" << sqrt(thickUncert)/PhysicalUnits::cm << " cm)" << endl
        //<< "\tmassFrac=" << massFrac << ", massFracUncert=" << massFracUncert
        //<< endl;
        
        activityUncertSquared += thisActivityUncertaintySquared;
        if( !IsNan(massFrac) && !IsInf(massFrac) && (massFrac > FLT_EPSILON) )
          activityUncertSquared += std::pow( thisActivity * massFracUncert / massFrac, 2.0 );
        
        activity += thisActivity;
      }//for( const SandiaDecay::Nuclide *src : srcs )
    }//for( const auto el_nucs : shield.m_nuclideFractions_ )
    
  
    // Add in uncertainty contributions if this is a trace source
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : shield.m_traceSources )
    {
      if( trace.m_nuclide != nuclide )
        continue;
      
      const TraceActivityType trace_type = trace.m_type;
      
      // Trace activity may be total, per cc, or per gram - but what the user really cares about
      //  is the uncertainty on the total (e.g., they dont care if it comes from trace activity
      //  or from dimension uncertainties), so we'll be a bit conservative here and just throw in
      //  all uncertainties, even if we are returning an uncertainty per cc or gram, this way when
      //  the the per cc or per gram gets multiplied by the volume, the expected uncertainty will
      //  be given.
      const double thisActivity = ShieldingSourceChi2Fcn::activity( nuclide, params );
      
      double thisActUncertSquared = niaveActUncert * niaveActUncert;
      switch( trace_type )
      {
        case TraceActivityType::TotalActivity:
        {
          // I *guess* we can ignore effects of volume uncertainty; this isnt strictly true due to
          //  1/r2 effects, but actually the activity uncertainty parameter should already take this
          //  into account - so I think we are good.
          activity += thisActivity;
          break;
        }//case TraceActivityType::TotalActivity:
          
          
        case TraceActivityType::ActivityPerCm3:
        case TraceActivityType::ActivityPerGram:
        {
          // Again, we are probably over-estimating the statistical error here, but since
          //  statistical errors probably always pale in comparison to systematic, I dont think this
          //  is an issue (or usually would people would just say this is being conservative, which
          //  no one ever questions).
          const double volFracUncert = volUncert / vol;
          
          if( trace_type == TraceActivityType::ActivityPerCm3 )
            activity += thisActivity * (vol/PhysicalUnits::cm3);
          else
            activity += thisActivity * vol * mat->density / PhysicalUnits::gram;
          
          thisActUncertSquared += volFracUncert*volFracUncert*thisActivity*thisActivity;
          break;
        }//case TraceActivityType::ActivityPerCm3 or ActivityPerGram
          
          
        case TraceActivityType::ExponentialDistribution:
        {
          switch( m_geometry )
          {
            case GeometryType::Spherical:
            {
              const double r = sphericalThickness(material_index, params);
              activity += thisActivity * 4.0 * PhysicalUnits::pi * r * r / PhysicalUnits::m2;
              break;
            }
              
            case GeometryType::CylinderEndOn:
            {
              const double r = cylindricalRadiusThickness(material_index, params);
              activity += thisActivity * PhysicalUnits::pi * r * r / PhysicalUnits::m2;
              break;
            }
              
            case GeometryType::CylinderSideOn:
            {
              const double r = cylindricalRadiusThickness(material_index, params);
              const double z = cylindricalLengthThickness(material_index, params);
              activity += thisActivity * PhysicalUnits::pi * 2.0 * r * z / PhysicalUnits::m2;
              break;
            }
              
            case GeometryType::Rectangular:
            {
              const double w = 2.0 * rectangularWidthThickness(material_index, params);
              const double h = 2.0 * rectangularHeightThickness(material_index, params);
              activity += thisActivity * w * h / PhysicalUnits::m2;
              break;
            }
              
            case GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( m_geometry )
          
          break;
        }//case TraceActivityType::ExponentialDistribution:
          
        case TraceActivityType::NumTraceActivityType:
          assert( 0 );
          throw runtime_error( "ShieldingSourceChi2Fcn::activityUncertainty: unexpected TraceActivityType" );
          break;
      }//switch( trace.second )
      
      
      activityUncertSquared += thisActUncertSquared;
    }//for( loop over trace sources )
  }//for( size_t material_index = 0; material_index < num_materials; ++material_index )
  
  if( activityUncertSquared < 0.0 || IsNan(activityUncertSquared) || IsInf(activityUncertSquared) )
    throw runtime_error( "error calculating activity uncertainty for "
                        + string(isTrace ? "trace" : (isSelfAtten ? "self-attenuating" : "uncertain!?!"))
                        + " source; squared value calculated is " + std::to_string(activityUncertSquared) );
  
  const double normalCalcAct = ShieldingSourceChi2Fcn::totalActivity( nuclide, params );
  assert( fabs(normalCalcAct - activity) < std::max(1.0E-6*PhysicalUnits::bq, 0.001*std::max(normalCalcAct,activity)) );
  
  return sqrt( activityUncertSquared );
}//double activityUncertainty(...)



double ShieldingSourceChi2Fcn::age( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
    throw runtime_error( "ShieldingSourceChi2Fcn::age(...) invalid params size" );

  const size_t ind = nuclideIndex( nuclide );

  //If parameter is zero or positive, we will use this age.
  if( params[2*ind+1] >= -0.00001 )
    return std::max( params[2*ind+1], 0.0 );
  
  //Else the parameter value indicates the "defining" nuclide whose age should be
  // specified.  To get defining index add one to value, and take negative
  const double findex = -1.0*params[2*ind+1];
  const int nearFIndex = static_cast<int>( std::round(findex) );
  
  if( fabs(findex - nearFIndex) > 0.01 || nearFIndex < 1 )
    throw runtime_error( "Got a negative age value that is not indicating a"
                         " defining nuclide age: " + std::to_string(params[2*ind+1]) );
  
  const int defining_nuclide_index = nearFIndex - 1;
  const int defining_nuc_age_index = 2*defining_nuclide_index + 1;
  
  if( static_cast<size_t>(defining_nuc_age_index) >= params.size() )
    throw runtime_error( "Got a negative age value that is larger than could be"
                         " for indicating a defining nuclide age: " + std::to_string(params[2*ind+1]) );
  
  const double defining_age = params[defining_nuc_age_index];
  if( defining_age < -0.00001 )
    throw runtime_error( "Defining age is also less than zero (shouldnt happen)" );
  
  return std::max( defining_age, 0.0 );
}//double age(...)


size_t ShieldingSourceChi2Fcn::nuclideIndex( const SandiaDecay::Nuclide *nuclide ) const
{
  vector<const SandiaDecay::Nuclide *>::const_iterator pos;
  pos = find( m_nuclides.begin(), m_nuclides.end(), nuclide );
  if( pos == m_nuclides.end() )
    throw std::runtime_error( "ShieldingSourceChi2Fcn::nuclideIndex(...): invalid nuclide" );

  return (pos - m_nuclides.begin());
}//size_t index( const SandiaDecay::Nuclide *nuclide ) const

  
size_t ShieldingSourceChi2Fcn::numNuclides() const
{
  return m_nuclides.size();
}
  
size_t ShieldingSourceChi2Fcn::numMaterials() const
{
  return m_initial_shieldings.size();
}//int numMaterials() const


const Material *ShieldingSourceChi2Fcn::material( const size_t materialNum ) const
{
  return m_initial_shieldings.at(materialNum).m_material.get();
}

bool ShieldingSourceChi2Fcn::isSpecificMaterial( const size_t materialNum ) const
{
  return (m_initial_shieldings.at(materialNum).m_material != nullptr );
}//bool isSpecificMaterial( int materialNum ) const


bool ShieldingSourceChi2Fcn::isGenericMaterial( const size_t materialNum ) const
{
  return (m_initial_shieldings.at(materialNum).m_material == nullptr );
}//bool isGenericMaterial( int materialNum ) const
  

double ShieldingSourceChi2Fcn::sphericalThickness( const size_t materialNum, const vector<double> &params ) const
{
  if( isGenericMaterial( materialNum ) )
    throw std::runtime_error( "ShieldingSourceChi2Fcn::sphericalThickness(...): "
                              "You should not call this function for generic "
                              "materials" );

  return params.at( 2*m_nuclides.size() + 3*materialNum );
}//double sphericalThickness(...)


double ShieldingSourceChi2Fcn::cylindricalRadiusThickness( const size_t materialNum, const std::vector<double> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum );
}//double cylindricalRadiusThickness(...)


double ShieldingSourceChi2Fcn::cylindricalLengthThickness( const size_t materialNum, const std::vector<double> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}//double cylindricalLengthThickness(...)


double ShieldingSourceChi2Fcn::rectangularWidthThickness( const size_t materialNum, const std::vector<double> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum );
}//double rectangularWidthThickness(...)


double ShieldingSourceChi2Fcn::rectangularHeightThickness( const size_t materialNum, const std::vector<double> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}//double rectangularHeightThickness(...)


double ShieldingSourceChi2Fcn::rectangularDepthThickness( const size_t materialNum, const std::vector<double> &params ) const
{
  assert( !isGenericMaterial( materialNum ) );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum + 2 );
}//double rectangularDepthThickness(...)


//arealDensity(...): will throw std::runtime_exception if material is a
//  specific material
double ShieldingSourceChi2Fcn::arealDensity( const size_t materialNum,
                       const std::vector<double> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "ShieldingSourceChi2Fcn::arealDensity(...): "
                              "You can only call this function for generic materials" );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum + 1 );
}//double arealDensity(...)


//atomicNumber(...): will throw std::runtime_exception if material is a
//  specific material
double ShieldingSourceChi2Fcn::atomicNumber( const size_t materialNum,
                       const std::vector<double> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "ShieldingSourceChi2Fcn::atomicNumber(...): "
                              "You can only call this function for generic materials" );
  
  return params.at( 2*m_nuclides.size() + 3*materialNum );
}//double atomicNumber(...)


const std::vector<PeakDef> &ShieldingSourceChi2Fcn::peaks() const
{
  return m_peaks;
}


const std::vector<PeakDef> &ShieldingSourceChi2Fcn::backgroundPeaks() const
{
  return m_backgroundPeaks;
}
  

double ShieldingSourceChi2Fcn::distance() const
{
  return m_distance;
}
  

const std::shared_ptr<const DetectorPeakResponse> &ShieldingSourceChi2Fcn::detector() const
{
  return m_detector;
}
  
const ShieldingSourceFitCalc::ShieldingSourceFitOptions &ShieldingSourceChi2Fcn::options() const
{
  return m_options;
}
  
const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &ShieldingSourceChi2Fcn::initialShieldings() const
{
  return m_initial_shieldings;
}//const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &initialShieldings() const
  

ShieldingSourceChi2Fcn::CancelException::CancelException( const ShieldingSourceChi2Fcn::CalcStatus cancel_code )
 : std::exception(),
  m_code( cancel_code )
{
}


double ShieldingSourceChi2Fcn::DoEval( const std::vector<double> &x ) const
{
  const CalcStatus cancelCode = m_cancel.load();
  switch( cancelCode )
  {
    case ShieldingSourceChi2Fcn::CalcStatus::NotCanceled:
      // Keep going
      break;
      
    case ShieldingSourceChi2Fcn::CalcStatus::UserCanceled:
    case ShieldingSourceChi2Fcn::CalcStatus::Timeout:
    case ShieldingSourceChi2Fcn::CalcStatus::CanceledNoUpdate:
      throw CancelException( cancelCode );
      break;
  }//switch( m_cancel.load() )

  try
  {
    // Check we're not being passed total non-sense here
    for( size_t i = 0; i < x.size(); ++i )
    {
      if( IsNan(x[i]) || IsInf(x[i]) )
      {
        throw runtime_error( "Invalid parameter (" + std::to_string(i) + ", "
                            + std::to_string(x[i]) + ") passed to chi2 fit fcn" );
      }
    }//for( size_t i = 0; i < x.size(); ++i )
    
    if( m_mixtureCache.size() > sm_maxMixtureCacheSize )
      m_mixtureCache.clear();
    
    const vector<PeakResultPlotInfo> chi2s
                          = energy_chi_contributions( x, {}, m_mixtureCache, nullptr, nullptr );
    double chi2 = 0.0;
    
    const size_t npoints = chi2s.size();
    for( size_t i = 0; i < npoints; ++i )
      chi2 += pow( chi2s[i].numSigmaOff, 2.0 );
    
    if( m_isFitting && m_guiUpdateInfo )
      m_guiUpdateInfo->completed_eval( chi2, x );
    
    //cout << "Returning chi2=" << chi2 << " for {" ;
    //for( size_t i = 0; i < x.size(); ++i )
    //  cout << (i ? "," : "") << x[i];
    //cout << "}" << endl;
    
    return chi2;
  }catch(...)
  {
  }

  return std::numeric_limits<double>::max();
}//double DoEval( const std::vector<double> &x ) const


namespace
{
//Would c++11 lambdas be awesome?
bool first_lessthan( const pair<double,double> &lhs, const pair<double,double> &rhs )
{
  return lhs.first < rhs.first;
}
}

//ToDo: add ability to give summary about 
void ShieldingSourceChi2Fcn::cluster_peak_activities( std::map<double,double> &energy_count_map,
                                                           const std::vector< pair<double,double> > &energie_widths,
                                                           SandiaDecay::NuclideMixture &mixture,
                                                           const double act,
                                                           const double age,
                                                           const double photopeakClusterSigma,
                                                           const double energyToCluster,
                                                           const bool accountForDecayDuringMeas,
                                                           const double measDuration,
                                                           vector<string> *info,
                                                           vector<PeakDetail> *log_info )
{
  typedef pair<double,double> DoublePair;

  if( energie_widths.empty() )
    return;
  
  if( energy_count_map.empty() )
  {
    for( const DoublePair &dp : energie_widths )
      energy_count_map[dp.first] = 0.0;
  }//if( energy_count_map.empty() )

  if( info )
  {
    stringstream msg;
    msg << "For ";
    for( int n = 0; n < mixture.numInitialNuclides(); ++n )
      msg << (n ? ", " : "") << mixture.initialNuclide(n)->symbol;
    msg << " at age " << PhysicalUnitsLocalized::printToBestTimeUnits(age) << ":";
    info->push_back( msg.str() );
  }//if( info )

  /*
  if( log_info )
  {
    for( int n = 0; n < mixture.numInitialNuclides(); ++n )
    {
      GammaInteractionCalc::ActShieldCalcLogInfo::SrcDef &src_def
                                          = log_info->m_sources[mixture.initialNuclide(n)->symbol];
      src_def.act = act;
      src_def.age = age;
    }
  }//if( log_info )
   */
  
  if( mixture.numInitialNuclides() != 1 )
    throw runtime_error( "ShieldingSourceChi2Fcn::cluster_peak_activities():"
                         " passed in mixture must have exactly one parent nuclide" );
  const SandiaDecay::Nuclide *nuclide = mixture.initialNuclide(0);
  assert( nuclide );
  if( !nuclide )
    throw std::logic_error( "nullptr nuc" );
  
  // We'll calculate our gammas differently if we are correcting for decays during the measurement 
  //  or not.
  //  We will only fill the `non_decay_cor_gammas` if we are correcting for decays during
  //  measurement, and we are creating log info (we only use non-corrected gammas to put
  //  the correction factor into the log file).
  vector<SandiaDecay::EnergyRatePair> gammas, non_decay_cor_gammas;
  
  assert( !accountForDecayDuringMeas || (measDuration > 0.0) );  //This could happen, I guess if foreground real-time is zero, but it shouldnt, so I'll leave this assert in to check for things.
  
  if( accountForDecayDuringMeas && (measDuration > 0.0) && !IsNan(measDuration) && !IsInf(measDuration) )
  {
    gammas = decay_during_meas_corrected_gammas( mixture, age, measDuration );
    
    // We will only use non-decay-corrected gammas if we are logging information
    if( info || log_info )
      non_decay_cor_gammas = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );
  }else
  {
    gammas = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );
  }//if( accountForDecayDuringMeas ) / else
  
  
  //The problem we have is that 'gammas' have the activity of the original
  //  parent ('nuclide') decreased by agining by 'age', however we want the
  //  parent to have 'sm_activityUnits' activity at 'age', so we we'll add a
  //  correction factor.
  const double age_sf = sm_activityUnits / mixture.activity(age, nuclide);
  
#if( PERFORM_DEVELOPER_CHECKS )
  double old_way_age_sf = 1.0;
  const vector<SandiaDecay::NuclideActivityPair> aged_activities
                                                      = mixture.activity( age );

  for( const SandiaDecay::NuclideActivityPair &red : aged_activities )
  {
    if( red.nuclide == nuclide )
    {
      assert( fabs(mixture.activity(age, nuclide) - red.activity) < 0.001*red.activity );
      
      old_way_age_sf = 1.0*sm_activityUnits / red.activity;
      break;
    }
  }//for( size_t i = 0; i < aged_activities.size(); ++i )
  
  assert( fabs(old_way_age_sf - age_sf) < 0.0001*old_way_age_sf );
#endif //#if( PERFORM_DEVELOPER_CHECKS )
  
  for( const SandiaDecay::EnergyRatePair &aep : gammas )
  {
    const pair<double,double> epair(aep.energy,0.0);
    vector< pair<double,double> >::const_iterator epos;
    epos = lower_bound( energie_widths.begin(), energie_widths.end(),
                        epair, &first_lessthan );

    if( photopeakClusterSigma<=0.0
        && (epos==energie_widths.end() || epos->first!=aep.energy) )  //queezy double compare
      continue;

    if( epos == energie_widths.end() )
    {
      //aep.energy is larger than the largest fit peak energy
      if( (aep.energy - energie_widths.back().first)
          > (energie_widths.back().second*photopeakClusterSigma) )
        continue;
      epos--;
    }else if( epos == energie_widths.begin() && epos->first!=aep.energy ) //queezy double compare
    {
      //aep.energy is less than the largest fit peak energy
      if( (epos->first - aep.energy) > (epos->second*photopeakClusterSigma) )
        continue;
    }else if( epos->first != aep.energy )
    {
      //see if the nearest peaks are close enough; if so, assign to closest
      //  peak
      const double next_width = epos->second;
      const double next_energy = epos->first;
      const double prev_width = (epos-1)->second;
      const double prev_energy = (epos-1)->first;

      const double prev_sigma = (aep.energy - prev_energy) / prev_width;
      const double next_sigma = (next_energy - aep.energy) / next_width;
/*
      cerr << "aep.energy=" << aep.energy << ", next_energy=" << next_energy
           << ", next_width=" << next_width << ", prev_energy="
           << prev_energy << ", prev_width=" << prev_width << endl
           << "prev_sigma=" << prev_sigma << ", next_sigma=" << next_sigma
           << " m_options.photopeak_cluster_sigma=" << photopeakClusterSigma << endl;
*/
      
      if( prev_sigma < next_sigma && prev_sigma < photopeakClusterSigma )
      {
        epos--;
//          cout << "Assigning (prev) " << aep.energy << " to fit peak " << epos->first << endl << endl;
      }else if( next_sigma < photopeakClusterSigma )
      {
        //epos is already the correct value
//          cout << "Assigning (next) " << aep.energy << " to fit peak " << epos->first << endl << endl;
      }else
      {
//          cout << "not assigning " << aep.energy << " to a fit photopeak" << endl << endl;
        continue;
      }
    }//if / else

    if( energyToCluster > 0.0 )
    {
      const double dist = fabs((energyToCluster - epos->first) / epos->second);
      if( dist > photopeakClusterSigma )
        continue;
    }//if( energyToCluster > 0.0 )

    
    const double energy = epos->first;

    if( energy_count_map.count( energy ) == 0 )
    {
      stringstream msg;
      msg << "There is a programming logic error in "
             "ShieldingSourceChi2Fcn::cluster_peak_activities(...)"
             " which is keeping this activity/shielding fit from happening "
             "(place a) - please complain to Will Johnson about this";
      passMessage( msg.str(), WarningWidget::WarningMsgHigh );
      throw std::runtime_error( msg.str() );
    }//if( energy_count_map.count( energy ) == 0 )

//    cerr << "Letting " << mixture.initialNuclide(0)->symbol << " contribute "
//         << (aep.numPerSecond * act * age_sf / sm_activityUnits) << " counts to "
//         << energy << " keV" << endl;
    
    const double contribution = aep.numPerSecond * act * age_sf / sm_activityUnits;
    energy_count_map[energy] += contribution;
    
    if( info )
    {
      stringstream msg;
      msg << "\tPeak attributed to " << energy << " keV received "
          << contribution*PhysicalUnits::second
          << " cps from " << aep.energy << " keV line, which has I="
          << age_sf * aep.numPerSecond/sm_activityUnits;
      
      if( !non_decay_cor_gammas.empty() )
      {
        //Find same-energy gamma, and get correction factor
        const auto non_corr_pos = std::find_if( begin(non_decay_cor_gammas), end(non_decay_cor_gammas),
          [&aep]( const SandiaDecay::EnergyRatePair &v ) {
            return fabs(v.energy - aep.energy) < 0.00001;
        });
        
        assert( non_corr_pos != end(non_decay_cor_gammas) );
        if( non_corr_pos != end(non_decay_cor_gammas) )
        {
          msg << " (decay correction " << aep.numPerSecond/non_corr_pos->numPerSecond << ")";
        }
      }//if( we are correcting for decays during measurement )
      
      info->push_back( msg.str() );
    }//if( info )
    
    if( log_info )
    {
      assert( mixture.numInitialNuclides() == 1 );
      const SandiaDecay::Nuclide * const nuc = mixture.initialNuclide(0);
      assert( nuc );
      
      auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const PeakDetail &p ) {
        return (p.decayParticleEnergy == energy);
      });
      
      assert( pos != end(*log_info) );
      if( pos != end(*log_info) )
      {
        PeakDetailSrc src;
        src.nuclide = nuc;
        src.energy = aep.energy;
        src.br = age_sf * aep.numPerSecond / sm_activityUnits;
        src.cpsAtSource = contribution*PhysicalUnits::second;
        src.age = age;
        src.calcActivity = act;
        src.decayCorrection = 0.0;
          
        if( !non_decay_cor_gammas.empty() )
        {
          //Find same-energy gamma, and get correction factor
          const auto non_corr_pos = std::find_if( begin(non_decay_cor_gammas), end(non_decay_cor_gammas),
                                                 [&aep]( const SandiaDecay::EnergyRatePair &v ) {
            return fabs(v.energy - aep.energy) < 0.00001;
          });
          
          assert( non_corr_pos != end(non_decay_cor_gammas) );
          if( non_corr_pos != end(non_decay_cor_gammas) )
            src.decayCorrection = aep.numPerSecond / non_corr_pos->numPerSecond;
        }//if( we are correcting for decays during measurement )
        
        // Note: we are only partually filling out the `PeakDetailSrc` object here; we will pick up the rest
        //       later, but this means we have to be careful to not rely on any of the member variables, notably:
        /*
        src.isTraceSource = false;
        src.traceSourceType = TraceActivityType::NumTraceActivityType;
        src.isSelfAttenSource = false;
        src.countsAtSource = 0.0;
        src.ageUncert = 0.0;
        src.activity = 0.0;
        src.activityUncert = 0.0;
        src.displayActivity = 0.0;
        src.displayActivityUncert = 0.0;
        src.massFraction = 0.0;
        src.massFractionUncert = 0.0;
        src.isFittingMassFraction = false;
        */
        
        pos->m_sources.push_back( src );
      }//if( pos != end(*log_info) )
    }//if( log_info )
    
  }//for( const SandiaDecay::AbundanceEnergyPair &aep : gammas )
/*
  cout << "For " << nuclide->symbol << " unshielded " << endl;
  for( EnergyCountMap::value_type &energy_count : energy_count_map )
    cout << "\t" << energy_count.first << "kev br="
         << energy_count.second/act
         << " rate=" << energy_count.second*PhysicalUnits::second << "/s"
         << endl;
*/
}//cluster_peak_activities(...)


vector<PeakResultPlotInfo> ShieldingSourceChi2Fcn::expected_observed_chis(
                                           const std::vector<PeakDef> &peaks,
                                           const std::vector<PeakDef> &backPeaks,
                                           const std::map<double,double> &energy_count_map,
                                           vector<string> *info,
                                           vector<GammaInteractionCalc::PeakDetail> *log_info )
{
  typedef map<double,double> EnergyCountMap;

  
  if( info )
    info->push_back( "Chi2 Contributions Of Peaks" );
  
  //Go through and match the predicted number of counts to the observed number
  //  of counts and get the chi2.
  //Note that matching between expected and observed peaks is done via energy
  //  which make me a bit queasy for some reason
  vector<PeakResultPlotInfo> answer;

  for( const PeakDef &peak : peaks )
  {
    if( !peak.decayParticle() && (peak.sourceGammaType()!=PeakDef::AnnihilationGamma) )
      continue;
    
    const double energy = peak.gammaParticleEnergy();
    
    double observed_counts = peak.peakArea();
    double observed_uncertainty = peak.peakAreaUncert();
    double backCounts = 0.0, backUncert2 = 0.0;
    const double nsigmaNear = 1.0;
    
    //XXX - this linear search is bad, and should be fixed
    for( const PeakDef &backPeak : backPeaks )
    {
      const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
      if( fabs(backPeak.mean()-peak.mean()) < (nsigmaNear*sigma) )
      {
        backCounts += backPeak.peakArea();
        backUncert2 += backPeak.peakAreaUncert()*backPeak.peakAreaUncert();
      }//if( fabs(backPeak.mean()-peak.mean()) < sigma )
    }//for( const PeakDef &peak : backPeaks )
    
    
    if( energy_count_map.count(energy) == 0 )
    {
      stringstream msg;
      msg << "There is a programming logic error in "
             "ShieldingSourceChi2Fcn::expected_observed_chis(...)"
             " which is keeping this activity/shielding fit from happening "
             "(place b) - please complain to Will Johnson about this";
      passMessage( msg.str(), WarningWidget::WarningMsgHigh );
      throw std::runtime_error( msg.str() );
    }//if( energy_count_map.count( energy ) == 0 )

    double expected_counts = 0.0;
    for( const EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
      if( fabs(energy_count.first-energy) < (0.1*sigma) )  //XXX - in principle we have already clustered phootpopeaks, and could just quicly access the energies
        expected_counts += energy_count.second;
    }

    if( backCounts > 0.0 )
    {
      // Background peak area and uncertainty have already been live-time normalized to foreground
      //  when set in `setBackgroundPeaks(...)`.
      observed_counts -= backCounts;
      observed_uncertainty = sqrt( observed_uncertainty*observed_uncertainty + backUncert2 );
    }
    
    const double chi = (observed_counts - expected_counts) / observed_uncertainty;
    const double scale = observed_counts / expected_counts;
    const double scale_uncert = observed_uncertainty / expected_counts;
    
    PeakResultPlotInfo peak_info;
    peak_info.energy = energy;
    peak_info.numSigmaOff = chi;
    peak_info.observedOverExpected = scale;
    peak_info.peakColor = peak.lineColor();
    peak_info.observedOverExpectedUncert = scale_uncert;
    answer.push_back( peak_info );
    
    if( info )
    {
      stringstream msg;
      msg << "\tAt " << energy/SandiaDecay::keV << " expected "
          << expected_counts << " counts, received " << observed_counts
          << " +- " << observed_uncertainty << " counts.";
      if( backCounts > 0.0 )
        msg << " (after correcting for " << backCounts << " +- "
            << sqrt(backUncert2) << " counts in background)";
      msg << " giving (observed-expected)/uncert=" << chi;
      info->push_back( msg.str() );
    }//if( info )
    
    if( log_info )
    {
      try
      {
        const double energy = peak.gammaParticleEnergy();
        
        auto pos = std::find_if( begin(*log_info), end(*log_info),
                                [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        
        assert( pos != end(*log_info) );
        
        if( pos != end(*log_info) )
        {
          GammaInteractionCalc::PeakDetail &log_peak = *pos;
          
          assert( log_peak.energy == peak.mean() );
          assert( log_peak.decayParticleEnergy == peak.gammaParticleEnergy() );
          assert( (peak.type() != PeakDef::GaussianDefined) || (log_peak.fwhm == peak.fwhm()) );
          assert( log_peak.counts == peak.peakArea() );
          assert( log_peak.countsUncert == peak.peakAreaUncert() );
          
          log_peak.expectedCounts = expected_counts;
          log_peak.observedCounts = observed_counts;
          log_peak.observedUncert = observed_uncertainty;
          
          log_peak.numSigmaOff = chi;
          log_peak.observedOverExpected = scale;
          log_peak.observedOverExpectedUncert = scale_uncert;
          
          //log_peak.modelInto4Pi = ;
          //log_peak.modelInto4PiCps = ;
          
          if( backCounts > 0 )
          {
            log_peak.backgroundCounts = backCounts;
            log_peak.backgroundCountsUncert = sqrt( backUncert2 );
          }
          
          double totalAtSourceCounts = 0.0;
          for( const PeakDetailSrc &psrc : log_peak.m_sources )
            totalAtSourceCounts += psrc.countsAtSource;
          
          for( PeakDetailSrc &psrc : log_peak.m_sources )
          {
            const double src_counts = expected_counts * psrc.countsAtSource / totalAtSourceCounts;
            psrc.modelContribToPeak = src_counts;
          }
        }//if( pos != end(*log_info) )
      }catch( std::exception & )
      {
        assert( 0 );
      }
    }//if( log_info )
  }//for( const PeakDef &peak : m_peaks )

  return answer;
}//expected_observed_chis(...)


vector< pair<double,double> > ShieldingSourceChi2Fcn::observedPeakEnergyWidths( const std::vector<PeakDef> &peaks )
{
  //Populate energy_count_map keys for the energies of interest (so we dont have
  //  to compute quantites for all photopeak energies, which would be really
  //  inefficint)
  vector< pair<double,double> > energie_widths;
  for( const PeakDef &peak : peaks )
  {
    try
    {
      const double energy = peak.gammaParticleEnergy();
      const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
      energie_widths.push_back( make_pair( energy, sigma ) );
    }catch( std::exception & )
    {
      cerr << "ShieldingSourceChi2Fcn::observedPeakEnergyWidths(...)\n\tIssue here" << endl;
    }
  }//for( const PeakDef &peak : peaks )

  sort( energie_widths.begin(), energie_widths.end(), &first_lessthan );

  return energie_widths;
}//observedPeakEnergyWidths()


void ShieldingSourceChi2Fcn::selfShieldingIntegration( DistributedSrcCalc &calculator )
{
  // TODO: if an exponential distribution with zero relaxation length, then just integrate over
  //       the surface
  //if( calculator.m_isInSituExponential && (calculator.m_inSituRelaxationLength < FLT_EPSILON) )
  //{
  //  // Integrate over surface instead of volume
  //}
  
  int ndim = -1;  //the number of dimensions of the integral.
  
  switch ( calculator.m_geometry )
  {
    case GeometryType::Spherical:
    case GeometryType::CylinderEndOn:
      ndim = 2;
      break;
      
    case GeometryType::CylinderSideOn:
    case GeometryType::Rectangular:
      ndim = 3;
      break;
      
    case GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch ( calculator.m_geometry )
  
  assert( (ndim == 2) || (ndim == 3) );
  
  void *userdata = (void *)&calculator;
  
  /** For an example problem, found setting the requested relative accuracy to larger values did decrease number of evaluations, but
   only ~30% less for a relative accuracy of 100 times less than our nominal
         epsrel {1e-6, 1e-5, 1e-4,1e-3, 1e-2},
         neval { 2413, 2413, 2159, 1905, 1397}
  
   TODO: maybe it makes sense to lower relative accuracy when the current Chi2 is pretty bad.
   */
  const double epsrel = 1e-4;  //the requested relative accuracy
  const double epsabs = -1.0;//1e-12; //the requested absolute accuracy
  //const int verbose = 0;
  //const int last = 4;  //use the importance function without smoothing (*I think*)
  const int mineval = 0; //the minimum number of integrand evaluations required.
  const int maxeval = 5000000; //the (approximate) maximum number of integrand evaluations allowed.

  int nregions, neval, fail;
  double error, prob;

  calculator.integral = 0.0;

  try
  {
    // For the moment, we know cylinders and rectange wont throw exception.
    // TODO: need to make it so DistributedSrcCalc::eval_spherical doesnt ever throw exception
    
    switch( calculator.m_geometry )
    {
      case GeometryType::Spherical:
        Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                                  Integrate::LastImportanceFcnt,
                                  mineval, maxeval, nregions, neval,
                                  fail, calculator.integral, error, prob );
        break;
        
      case GeometryType::CylinderEndOn:
        if( calculator.m_dimensionsTransLenAndType.size() == 1 )
        {
          // For a single end-on cylinder we can use the ever-so-slightly faster function
          //  to evaluate this (debug builds will validate gives same answer as
          //  DistributedSrcCalc::eval_cylinder)
          
          Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_single_cyl_end_on, userdata, epsrel, epsabs,
                                    Integrate::LastImportanceFcnt,
                                    mineval, maxeval, nregions, neval,
                                    fail, calculator.integral, error, prob );
        }else
        {
          Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_cylindrical, userdata, epsrel, epsabs,
                                    Integrate::LastImportanceFcnt,
                                    mineval, maxeval, nregions, neval,
                                    fail, calculator.integral, error, prob );
        }
        break;
        
      case GeometryType::CylinderSideOn:
        Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_cylindrical, userdata, epsrel, epsabs,
                                  Integrate::LastImportanceFcnt,
                                  mineval, maxeval, nregions, neval,
                                  fail, calculator.integral, error, prob );
        break;
        
      case GeometryType::Rectangular:
        Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_rectangular, userdata, epsrel, epsabs,
                                  Integrate::LastImportanceFcnt,
                                  mineval, maxeval, nregions, neval,
                                  fail, calculator.integral, error, prob );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( objToIntegrate->m_geometry )
  }catch( std::exception &e )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    stringstream msg;
    msg << "Integration failed after " << neval
    << " evaluations got integral " << calculator.integral
    << "\n\t calculator.energy: " << calculator.m_energy
    << "\n\t m_srcVolumetricActivity: " << calculator.m_srcVolumetricActivity
    << "\n\t m_geometry: " << to_str(calculator.m_geometry)
    << "\n\t m_materialIndex: " << calculator.m_materialIndex
    << "\n\t m_detectorRadius: " << calculator.m_detectorRadius
    << "\n\t m_observationDist: " << calculator.m_observationDist
    << "\n\t m_attenuateForAir: " << calculator.m_attenuateForAir
    << "\n\t m_airTransLenCoef: " << calculator.m_airTransLenCoef
    << "\n\t m_isInSituExponential: " << calculator.m_isInSituExponential
    << "\n\t m_inSituRelaxationLength: " << calculator.m_inSituRelaxationLength
    << "\n\t m_dimensionsTransLenAndType: {";
    for( const auto &i : calculator.m_dimensionsTransLenAndType )
    {
      const auto &dims = std::get<0>(i);
      msg << "{[" << dims[0] << "," << dims[1] << "," << dims[2] << "],"
          << std::get<1>(i) << "," << static_cast<int>(std::get<2>(i)) << "}, ";
    }
    
    msg << "}"
    << "\n\t nuclide: " << (calculator.m_nuclide ? calculator.m_nuclide->symbol : string())
    << "\n\t exception: '" << e.what() << "'"
    << endl;
      
    log_developer_error( __func__, msg.str().c_str() );
    cerr << __func__ << msg.str() << endl;
#endif
    
    cerr << "Integration failed: " << e.what() << endl;
    calculator.integral = 0.0;
  }//try / catch
  

  
  //should check for fails and stuff all over here
  
/*
  static std::mutex m;
  std::lock_guard<std::mutex> scoped_lock( m );
  cerr << "After " << neval << " evaluations got integral " << calculator.integral
       << endl << "\t calculator.energy=" << calculator.energy
       << " calculator.m_srcVolumetricActivity=" << calculator.m_srcVolumetricActivity
       << ", rad=" << get<0>(calculator.m_dimensionsTransLenAndType[0])[0]
       << ", transCoef=" << calculator.m_sphereRadAndTransLenCoef[0].second
       << endl;
*/
}//void selfShieldingIntegration(...)

  
void ShieldingSourceChi2Fcn::setBackgroundPeaks(
                                              const std::vector<PeakDef> &peaks,
                                              double liveTime )
{
  m_backgroundPeaks.clear();
  const double scale = m_liveTime / liveTime;
  
  for( PeakDef p : peaks )
  {
    if( p.gausPeak() )
    {
      p.setPeakArea( scale * p.peakArea() );
      p.setPeakAreaUncert( scale * p.peakAreaUncert() );
      m_backgroundPeaks.push_back( p );
    }else
    {
      stringstream msg;
      msg << "The non-gaussian background peak at " << p.mean() << " keV "
          << " will not be used for background peak area subtraction;"
          << " non-gaussian peaks may be supported for this in the future";
      passMessage( msg.str(), WarningWidget::WarningMsgHigh );
    }//if( p.gausPeak() ) / else
  }//for( PeakDef &p : m_backgroundPeaks )
}//void setBackgroundPeaks(...)
  
  
vector<PeakResultPlotInfo>
       ShieldingSourceChi2Fcn::energy_chi_contributions( const std::vector<double> &x,
                                                        const std::vector<double> &error_params,
                                         ShieldingSourceChi2Fcn::NucMixtureCache &mixturecache,
                                         std::vector<std::string> *info,
                                         std::vector<GammaInteractionCalc::PeakDetail> *log_info ) const
{
  //XXX - this function compares a lot of doubles, and this always makes me
  //      queezy - this should be checked on!
  typedef map<double,double> EnergyCountMap;

  // Make sure attenuate_for_air isnt set for fixed geometry DRFs
  assert( !m_options.attenuate_for_air || !m_detector || !m_detector->isFixedGeometry() );

  assert( !log_info || (x.size() == error_params.size()) );
  
//  cerr << "energy_chi_contributions: vals={ ";
//  for( size_t i = 0; i < x.size(); ++i )
//    cerr << x[i] << ", ";
//  cerr << "}; size=" << x.size() << endl;
  
  if( info )
  {
    info->push_back( "LiveTime="
                    + std::to_string(m_liveTime/PhysicalUnits::second)
                    + " s (" + PhysicalUnitsLocalized::printToBestTimeUnits(m_liveTime)
                    + ")" );
    info->push_back( "Distance to source center from detector: "
                      + PhysicalUnits::printToBestLengthUnits(m_distance) );
    if( m_detector && m_detector->isValid() )
    {
      if( m_detector->isFixedGeometry() )
      {
        info->push_back( "Detector: " + m_detector->name() + ", fixed geometry" );
      }else
      {
        info->push_back( "Detector: " + m_detector->name() + " radius "
                        + std::to_string( 0.5*m_detector->detectorDiameter()/PhysicalUnits::cm )
                        + " cm" );
      }
    }//if( m_detector && m_detector->isValid() )
    
    if( m_options.multiple_nucs_contribute_to_peaks )
      info->push_back( "Allowing multiple nuclides being fit for to potentially contribute to the same photopeak" );
    else
      info->push_back( "Not allowing multiple nuclides being fit for to contribute to the same photopeak" );
    
    if( m_options.account_for_decay_during_meas )
      info->push_back( "Branching ratios are being corrected for nuclide decay during measurement" );
    
    //Should put in information about the shielding here
  }//if( info )
  
  if( log_info )
  {
    log_info->resize( m_peaks.size() );
    
    for( size_t i = 0; i < m_peaks.size(); ++i )
    {
      const PeakDef &peak = m_peaks[i];
      
      try
      {
        GammaInteractionCalc::PeakDetail &log_peak = (*log_info)[i];
        
        log_peak.energy = peak.mean();
        log_peak.decayParticleEnergy = peak.gammaParticleEnergy();
        if( peak.type() == PeakDef::GaussianDefined )
        {
          log_peak.fwhm = peak.fwhm();
        }
        
        log_peak.counts = peak.peakArea();
        log_peak.countsUncert = peak.peakAreaUncert();
        
        log_peak.cps = log_peak.counts / m_liveTime;
        log_peak.cpsUncert = log_peak.countsUncert / m_liveTime;
        
        if( peak.parentNuclide() )
          log_peak.assignedNuclide = peak.m_parentNuclide->symbol;
        else if( peak.xrayElement() )
          log_peak.assignedNuclide = peak.xrayElement()->name;
        else if( peak.reaction() )
          log_peak.assignedNuclide = peak.reaction()->name();
        else
        {
          assert( 0 );
          log_peak.assignedNuclide = "null";
        }
      }catch( std::exception & )
      {
        assert( 0 );
      }
    }//for( const PeakDef &peak : m_peaks )
  }//if( log_info )
  
  
  EnergyCountMap energy_count_map;
  const vector<pair<double,double> > energie_widths = observedPeakEnergyWidths( m_peaks );
  
  if( m_options.multiple_nucs_contribute_to_peaks )
  {
    //Get the number of source gammas from each nuclide
    for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
    {
      if( isVolumetricSource(nuclide) )
        continue;
      
      const double act = activity( nuclide, x );
      const double thisage = age( nuclide, x );
      
      //{
      //  static std::mutex s_debug_mutex;
      //  std::lock_guard<std::mutex> lock( s_debug_mutex );
      //  cout << "\tActivity of " << nuclide->symbol << " is "
      //       << PhysicalUnits::printToBestActivityUnits(act)
      //       << " with age " << PhysicalUnits::printToBestTimeUnits(thisage) << endl;
      //}
      
      
      if( mixturecache.find(nuclide) == mixturecache.end() )
        mixturecache[nuclide].addNuclideByActivity( nuclide, sm_activityUnits );
      
      cluster_peak_activities( energy_count_map, energie_widths,
                               mixturecache[nuclide], act, thisage,
                               m_options.photopeak_cluster_sigma, -1.0,
                               m_options.account_for_decay_during_meas, m_realTime,
                               info, log_info );
    }//for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
  }else
  {
    for( const PeakDef &peak : m_peaks )
    {
      const SandiaDecay::Nuclide *nuclide = peak.parentNuclide();
      const SandiaDecay::RadParticle *particle = peak.decayParticle();
      
      if( !nuclide
          || (!particle && (peak.sourceGammaType() != PeakDef::AnnihilationGamma))
          || isVolumetricSource(nuclide) )
        continue;
      
      const double act = activity( nuclide, x );
      const double thisage = age( nuclide, x );
      
      if( mixturecache.find(nuclide) == mixturecache.end() )
        mixturecache[nuclide].addNuclideByActivity( nuclide, sm_activityUnits );
      
      const float energy = peak.gammaParticleEnergy();
      cluster_peak_activities( energy_count_map, energie_widths,
                               mixturecache[nuclide], act, thisage,
                               m_options.photopeak_cluster_sigma, energy,
                               m_options.account_for_decay_during_meas, m_realTime,
                               info, log_info );
    }//for( const PeakDef &peak : m_peaks )
  }//if( m_options.multiple_nucs_contribute_to_peaks )

  //Propagate the gammas through each material - note we are using the fit peak
  //  mean here, and not the (pre-cluster) photopeak energy
  double shield_outer_rad = 0.0;
  const size_t nMaterials = m_initial_shieldings.size();
  
  for( size_t materialN = 0; materialN < nMaterials; ++materialN )
  {
    boost::function<double(float)> att_coef_fcn;
    const ShieldingSourceFitCalc::ShieldingInfo &shielding = m_initial_shieldings[materialN];
    const shared_ptr<const Material> &material = shielding.m_material;

    if( !material )
    {
      // Generic material here.
      
      float atomic_number = static_cast<float>(atomicNumber( materialN, x ));
      float areal_density = static_cast<float>(arealDensity( materialN, x ));
      
      // Even though we always set AD and AN limits on the parameters when creating them in
      //  ShieldingSourceDisplay::shieldingFitnessFcn(...), sometimes Minuit will go (way!) outside
      //  those bounds, and then everything becomes wack - so we'll clamp the values of AN/AD.
      if( atomic_number < MassAttenuation::sm_min_xs_atomic_number )
        atomic_number = MassAttenuation::sm_min_xs_atomic_number;
      if( atomic_number > MassAttenuation::sm_max_xs_atomic_number )
        atomic_number = MassAttenuation::sm_max_xs_atomic_number;
      
      const double ad_in_gcm2 = areal_density * PhysicalUnits::cm2 / PhysicalUnits::g;
      if( ad_in_gcm2 < 0.0 )
        areal_density = 0.0;
      if( ad_in_gcm2 > sm_max_areal_density_g_cm2 )
        areal_density = static_cast<float>(sm_max_areal_density_g_cm2*PhysicalUnits::g/PhysicalUnits::cm2);
  
      att_coef_fcn = boost::bind( &transmition_coefficient_generic, atomic_number, areal_density,
                                 boost::placeholders::_1 );
    }else
    {
      double thickness = 0.0;
      switch( m_geometry )
      {
        case GeometryType::Spherical:
          thickness = sphericalThickness(materialN,x);
          break;
          
        case GeometryType::CylinderEndOn:
          thickness = cylindricalLengthThickness(materialN,x);
          break;
          
        case GeometryType::CylinderSideOn:
          thickness = cylindricalRadiusThickness(materialN,x);
          break;
          
        case GeometryType::Rectangular:
          thickness = rectangularDepthThickness(materialN,x);
          break;
          
        case GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
      
      shield_outer_rad += thickness;
      
      att_coef_fcn = boost::bind( &transmition_coefficient_material, material.get(),
                                 boost::placeholders::_1, static_cast<float>(thickness) );
    }//if( generic material ) / else

/*
    //Uncomment out this section if you want to print out the transmission
    //  coefficients
    cout << "On material " << materialN << endl;
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      //adxs is in
      double adxs = transmition_length_coefficient( material, energy_count.first );
      adxs /= material->density;
      adxs *= PhysicalUnits::gram / PhysicalUnits::cm2;
      cout << energy_count.first << ", " << adxs << endl;
    }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
*/

    assert( att_coef_fcn );
    
    if( info )
    {
      if( isGenericMaterial(materialN) )
      {
        const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
        const double ad = arealDensity(materialN, x) / adUnits;
        const double an = atomicNumber(materialN, x);
        
        char buffer[256];
        snprintf( buffer, sizeof(buffer),
                 "Shielding %i: AN=%.2f, AD=%.3f g/cm2", static_cast<int>(materialN), an, ad );
        info->push_back( buffer );
      }else
      {
        stringstream title;
        title << "Shielding: " << material->name << ", "
              << material->chemicalFormula() << ", density="
              << material->density *PhysicalUnits::cm3/PhysicalUnits::gram
              << " g/cm3, ";
        
        title << to_str(m_geometry) << ": {";
        switch( m_geometry )
        {
          case GeometryType::Spherical:
          {
            const double thickness = sphericalThickness(materialN,x);
            title << "rad_thickness=" << PhysicalUnits::printToBestLengthUnits(thickness);
            break;
          }
          
          case GeometryType::CylinderSideOn:
          case GeometryType::CylinderEndOn:
          {
            const double r = cylindricalRadiusThickness(materialN,x);
            const double z = cylindricalLengthThickness(materialN,x);
            
            title << "rad_thickness=" << PhysicalUnits::printToBestLengthUnits(r)
                  << ", len_thickness=" << PhysicalUnits::printToBestLengthUnits(z);
            break;
          }
            
          case GeometryType::Rectangular:
          {
            const double w = rectangularWidthThickness(materialN,x);
            const double h = rectangularHeightThickness(materialN,x);
            const double d = rectangularDepthThickness(materialN,x);
            
            title << "width_thickness=" << PhysicalUnits::printToBestLengthUnits(w)
                  << ", height_thickness=" << PhysicalUnits::printToBestLengthUnits(h)
                  << ", depth_thickness=" << PhysicalUnits::printToBestLengthUnits(d);
            break;
          }
            
          case GeometryType::NumGeometryType:
            assert( 0 );
            break;
        }//switch( m_geometry )
        
        title << "}";
        
        info->push_back( title.str() );
      }//if( isGenericMaterial(materialN) ) / else
      
      for( EnergyCountMap::value_type &energy_count : energy_count_map )
      {
        if( energy_count.second <= 0.0 )
          continue;
        
        const double f = exp( -1.0 * att_coef_fcn( energy_count.first ) );
        stringstream msg;
        msg << "\tReduced counts of " << energy_count.first
            << " keV photopeak by " << f << " from "
            << energy_count.second * PhysicalUnits::second << " cps, to "
            << f * energy_count.second * PhysicalUnits::second << " cps";
        info->push_back( msg.str() );
      }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
    }//if( info )
    
    
    if( log_info )
    {
      for( EnergyCountMap::value_type &energy_count : energy_count_map )
      {
        if( energy_count.second <= 0.0 )
          continue;
       
        const double energy = energy_count.first;
        auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        
        assert( pos != end(*log_info) );
        
        if( pos != end(*log_info) )
        {
          const double f = exp( -1.0 * att_coef_fcn( energy_count.first ) );
          pos->m_attenuations.resize( nMaterials, 0.0 );
          pos->m_attenuations[materialN] = f;
          pos->m_totalAttenFactor *= f;
          pos->m_totalShieldAttenFactor *= f;
        }
      }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
    }//if( log_info )
    
    // Note: only non-volumetric sources are included in `energy_count_map`
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      const double energy = energy_count.first;
      const double frac = exp( -1.0 * att_coef_fcn( energy ) );
      
      energy_count.second *= frac;
    }
  }//for( int materialN = 0; materialN < nMaterials; ++materialN )

  const double air_dist = std::max( 0.0, m_distance - shield_outer_rad );
  if( m_options.attenuate_for_air && (!m_detector || !m_detector->isFixedGeometry()) )
  {
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      const double coef = transmission_length_coefficient_air( energy_count.first );
      const double atten = exp( -1.0 * coef * air_dist );
      energy_count.second *= atten;
      
      if( log_info )
      {
        const double energy = energy_count.first;
        auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        assert( pos != end(*log_info) );
        if( pos != end(*log_info) )
        {
          pos->m_airAttenFactor = atten;
          pos->m_totalAttenFactor *= atten;
        }
      }//if( log_info )
    }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
  }//if( m_options.attenuate_for_air )
  
  
  //Fold in the detector response
  if( m_detector && m_detector->isValid() )
  {
    const bool fixed_geom = m_detector->isFixedGeometry();
    if( info )
      info->push_back( "Detector Efficiency Effects"
                      + string(fixed_geom ? " (Fixed Geometry)" : "") );
    
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
//      cerr << "Absolute efficiency at " << energy_count.first << " keV is "
//           << m_detector->intrinsicEfficiency( energy_count.first ) << " and the "
//           << " total efficiency is " << m_detector->efficiency( energy_count.first, m_distance ) << endl;
      
      const double eff = fixed_geom ? m_detector->intrinsicEfficiency(energy_count.first)
                                    : m_detector->efficiency( energy_count.first, m_distance );
      
      if( info )
      {
        if( energy_count.second > 0.0 )
        {
          double deteff = m_detector->intrinsicEfficiency( energy_count.first );
          stringstream msg;
          msg << "\t" << energy_count.first << " keV photopeak reduced by "
              << eff/deteff << " * " << deteff
              << " (solid angle)*(det intrinsic eff) from "
              << energy_count.second*PhysicalUnits::second << " cps "
              << "to " << energy_count.second*PhysicalUnits::second*eff << " cps";
          info->push_back( msg.str() );
        }else
        {
          double deteff = m_detector->intrinsicEfficiency( energy_count.first );
          stringstream msg;
          msg << "\t" << energy_count.first << " keV photopeak reduced by "
              << deteff << " by the detectors absolute efficiency.";
          info->push_back( msg.str() );
        }
      }//if( info )
      
      
      if( log_info )
      {
        const double energy = energy_count.first;
        auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        
        assert( pos != end(*log_info) );
        
        if( pos != end(*log_info) )
        {
          const double deteff = m_detector->intrinsicEfficiency( energy_count.first );
          pos->detEff = eff;
          pos->detIntrinsicEff = deteff;
          pos->detSolidAngle = eff / deteff;
        }
      }//if( log_info )
      
      energy_count.second *= eff;
    }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
  }else 
  {
    const double detDiam = 1.0 * PhysicalUnits::cm;
    const double r = 0.5 * detDiam;
    const double D = m_distance;
    const double fracAngle = 0.5*(1.0 - (D/sqrt(D*D+r*r)));
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
      energy_count.second *= fracAngle;
    
    if( info )
      info->push_back( "Solid angle reduces counts by a factor of " + std::to_string(fracAngle) );
    
    if( log_info )
    {
      for( EnergyCountMap::value_type &energy_count : energy_count_map )
      {
        const double energy = energy_count.first;
        auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        
        assert( pos != end(*log_info) );
        
        if( pos != end(*log_info) )
        {
          pos->detEff = 1.0;
          pos->detIntrinsicEff = 1.0;
          pos->detSolidAngle = 1.0;
        }
      }//for( EnergyCountMap::value_type &energy_count : energy_count_map )
    }//if( log_info )
  }//if( m_detector && m_detector->isValid() ) / else


  //This is where contributions from self-attenuating and traces source are calculated
  if( info )
  {
    for( const auto &shield : m_initial_shieldings )
    {
      if( !shield.m_material )
        continue;
      
      size_t num_volume_src = 0;
      for( const auto &el_nucs : shield.m_nuclideFractions_ )
      {
        for( const auto &nucs : el_nucs.second )
          num_volume_src += (get<0>(nucs) != nullptr);
      }
      num_volume_src += shield.m_traceSources.size();
      
      if( num_volume_src )
      {
        info->push_back( "Self Attenuating Source Info (shielding, detector,"
                        " distance, live time, and amount of material not-accounted for):" );
        break;
      }
    }//for( int materialN = 0; materialN < nMaterials; ++materialN )
  }//if( info )

  using GammaInteractionCalc::transmition_length_coefficient;
  
  // For mass-varied materials, particularly involving multiple elements, we need to
  //  create a version of the material that will give the correct cross-section, for
  //  the current mass-fraction variation
  vector<std::shared_ptr<const Material>> local_materials;
  
  vector<std::unique_ptr<DistributedSrcCalc>> calculators;
  std::map<const DistributedSrcCalc *, bool> is_trace_src; //only used for logging
  std::map<const DistributedSrcCalc *, size_t> vol_src_material_index; //only used for logging
  
  
  bool has_trace = false, has_self_atten = false;
  
  for( size_t material_index = 0; material_index < nMaterials; ++material_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = m_initial_shieldings[material_index];
    shared_ptr<const Material> material = shield.m_material;
    
    if( !material )
      continue;
    
    vector<const SandiaDecay::Nuclide *> combined_srcs;
    const vector<ShieldingSourceFitCalc::TraceSourceInfo> &trace_srcs = shield.m_traceSources;
    for( const auto &p : trace_srcs )
    {
      assert( p.m_nuclide );
      combined_srcs.push_back( p.m_nuclide );
    }
    
    for( const auto &el_nucs : shield.m_nuclideFractions_ )
    {
      for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc : el_nucs.second )
      {
        const SandiaDecay::Nuclide * const nuclide = get<0>(nuc);
        if( nuclide )
          combined_srcs.push_back( nuclide );
      }//for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc : el_nucs.second )
    }//for( const auto &el_nucs : shield.m_nuclideFractions_ )
    
    const auto &srcs = combined_srcs;
    
    if( srcs.empty() )
      continue;
     
    assert( !m_detector || !m_detector->isFixedGeometry() );
    if( m_detector && m_detector->isFixedGeometry() )
      throw logic_error( "Self-attenuating and trace sources are not allowed for fixed geometry detector response functions." );
    
    DistributedSrcCalc baseCalculator;
    baseCalculator.m_geometry = m_geometry;
    
    if( m_detector )
      baseCalculator.m_detectorRadius = 0.5 * m_detector->detectorDiameter();
    else
      baseCalculator.m_detectorRadius = 0.5 * PhysicalUnits::cm;

    baseCalculator.m_observationDist = m_distance;
    baseCalculator.m_attenuateForAir = m_options.attenuate_for_air;
    baseCalculator.m_materialIndex = material_index;
    
    baseCalculator.m_isInSituExponential = false;
    baseCalculator.m_inSituRelaxationLength = -1.0;
    
    for( size_t src_index = 0; src_index < srcs.size(); ++src_index )
    {
      const SandiaDecay::Nuclide *src = srcs[src_index];
      const bool is_trace = (src_index < trace_srcs.size());
      
#if( PERFORM_DEVELOPER_CHECKS )
      {// begin quick sanity check
        bool trace_check = false;
        for( const auto &p : trace_srcs )
          trace_check = (trace_check || (p.m_nuclide == src));
        assert( trace_check == is_trace );
      }// end quick sanity check
#endif
      
      double actPerVol = 0.0;
      EnergyCountMap local_energy_count_map;
      
      if( is_trace )
      {
        has_trace = true;
        assert( src_index < trace_srcs.size() );
        
        const double act = activity(src, x);
        
        switch( trace_srcs[src_index].m_type )
        {
          case TraceActivityType::TotalActivity:
          {
            const double vol = volumeOfMaterial(material_index, x);
            actPerVol = act / vol;
            break;
          }
            
          case TraceActivityType::ActivityPerCm3:
            actPerVol = act / PhysicalUnits::cm3;
            break;
            
          case TraceActivityType::ActivityPerGram:
            actPerVol = act * material->density / PhysicalUnits::g;
            break;
            
          case TraceActivityType::ExponentialDistribution:
          {
            const double relaxation_len = trace_srcs[src_index].m_relaxationDistance;
            assert( relaxation_len > 0.0 );
            
            //actually activity of the entire soil-column, all the way down, so not an actPerVol,
            //  strictly, but per surface area
            actPerVol = act / PhysicalUnits::m2;
            
            // We need to normalize the activity by integrating over the entire volume
            switch( m_geometry )
            {
              case GeometryType::Spherical:
              {
                // Integral from 0.0 to shielding radius R, of "r*r*exp(-(R-r)/L)" where L is
                //  relaxation length, is L*(L*L*(2-2*exp(-R/L)) - 2*L*R + R*R)
                //  So we will normalize by this factor (times 4pi) so the surface contamination
                //  level (in activity/m2), multiplied by surface are will give total activity
                
                //4 π integral_0^R e^(-(R - ρ)/L) ρ^2 dρ = 4 π L (L^2 (2 - 2 e^(-R/L)) - 2 L R + R^2)
                
                const double R = sphericalThickness(material_index, x);
                const double L = relaxation_len;
                const double norm = 4*PhysicalUnits::pi * L * (L*L*(2 - 2*exp(-R/L)) - 2*L*R + R*R);
                actPerVol /= norm;
                
                break;
              }//case GeometryType::Spherical:
                
              case GeometryType::Rectangular:
              case GeometryType::CylinderEndOn:
              {
                // Integral from 0.0 to shielding depth R, of relaxation length L, is L - L*exp(-R/L)
                //  So we will normalize by this factor so the activity represents the entire activity
                //  of the column
                
                double R = -1.0;
                if( m_geometry == GeometryType::Rectangular )
                  R = 2.0 * rectangularDepthThickness(material_index, x);
                else
                  R = 2.0 * cylindricalLengthThickness(material_index, x);

                assert( R >= 0.0 );
                
                const double L = relaxation_len;
                const double norm = L * (1.0 - exp(-R / L) );
                actPerVol /= norm;
                break;
              }// case GeometryType::Rectangular or CylinderEndOn
                
              case GeometryType::CylinderSideOn:
              {
                //integrate 2*pi*r*exp(-(R-r)/L) wrt r, from 0 to R
                const double R = cylindricalRadiusThickness(material_index, x);
                const double L = relaxation_len;
                
                const double norm = 2 * L * PhysicalUnits::pi * (L*(exp(-R/L) - 1) + R);
                actPerVol /= norm;
                break;
              }//case GeometryType::CylinderSideOn:
              
              case GeometryType::NumGeometryType:
                assert( 0 );
                break;
            }//switch( m_geometry )
            
            break;
          }//case TraceActivityType::ExponentialDistribution:
            
          case TraceActivityType::NumTraceActivityType:
            assert(0);
            throw runtime_error("");
            break;
        }//switch ( trace_srcs[src_index].second )
      }else //if( is_trace )
      {
        has_self_atten = true;
        
        const double actPerMass = src->activityPerGram() / PhysicalUnits::gram;
        const double massFract = massFractionOfElement( material_index, src, x );
        
        actPerVol = actPerMass * massFract * material->density;
      }//if( is_trace ) / else
      
      
      const double thisage = age( src, x );
      
      //      cerr << "For " << src->symbol << " actPerVol=" << actPerVol
      //           << ", material->density=" << material->density << "("
      //           << material->density *PhysicalUnits::cm3/PhysicalUnits::g
      //           << "), actPerMass=" << actPerMass << ", massFraction="
      //           << massFraction << endl;
      
      if( mixturecache.find(src) == mixturecache.end() )
        mixturecache[src].addNuclideByActivity( src, sm_activityUnits );
      
      if( m_options.multiple_nucs_contribute_to_peaks )
      {
        cluster_peak_activities( local_energy_count_map, energie_widths,
                                 mixturecache[src], actPerVol, thisage,
                                 m_options.photopeak_cluster_sigma, -1.0,
                                 m_options.account_for_decay_during_meas, m_realTime,
                                 info, log_info );
      }else
      {
        for( const PeakDef &peak : m_peaks )
        {
          if( (peak.parentNuclide() == src)
             && (peak.decayParticle() || (peak.sourceGammaType() == PeakDef::AnnihilationGamma)) )
          {
            cluster_peak_activities( local_energy_count_map, energie_widths,
                                    mixturecache[src], actPerVol, thisage,
                                    m_options.photopeak_cluster_sigma,
                                    peak.gammaParticleEnergy(),
                                    m_options.account_for_decay_during_meas, m_realTime,
                                    info, log_info );
          }
        }//for( const PeakDef &peak : m_peaks )
      }//if( m_options.multiple_nucs_contribute_to_peaks ) / else

      for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
      {
        std::unique_ptr<DistributedSrcCalc> calculator( new DistributedSrcCalc(baseCalculator) );

        calculator->m_nuclide = src;
        calculator->m_energy = energy_count.first;
        calculator->m_srcVolumetricActivity = energy_count.second;
        
        if( m_options.attenuate_for_air )
          calculator->m_airTransLenCoef = transmission_length_coefficient_air( energy_count.first );
        else
          calculator->m_airTransLenCoef = 0.0;
        
        if( is_trace )
        {
          assert( src_index < trace_srcs.size() );
          
          if( trace_srcs[src_index].m_type == TraceActivityType::ExponentialDistribution )
          {
            calculator->m_isInSituExponential = true;
            calculator->m_inSituRelaxationLength = trace_srcs[src_index].m_relaxationDistance;
            assert( calculator->m_inSituRelaxationLength > 0.0 );
          }
        }//if( is_trace )
        
        std::array<double,3> outer_dims = { 0.0, 0.0, 0.0 };
        
        for( int subMat = 0; subMat < nMaterials; ++subMat )
        {
          if( isGenericMaterial( subMat ) )
          {
            // A generic material at the very center has zero volume, so skip it
            if( subMat == 0 )
              continue;
            
            const double an = atomicNumber(subMat, x);
            const double ad = arealDensity(subMat, x);
            const double transLenCoef = transmition_coefficient_generic( an, ad, calculator->m_energy );
            
            //cout << "Adding generic material (index=" << subMat << ") with AN=" << an
            //     << " and AD=" << ad / (PhysicalUnits::g/PhysicalUnits::cm2) << " g/cm2"
            //     << " atten(" << calculator->m_energy << " keV-->" << transLenCoef << ") = " << exp(-1.0*transLenCoef)
            //     << endl;
            
#if( defined(__GNUC__) && __GNUC__ < 5 )
            calculator->m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Generic} );
#else
            calculator->m_dimensionsTransLenAndType.push_back( {outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Generic} );
#endif
            
            continue;
          }//if( isGenericMaterial( subMat ) )

          switch( m_geometry )
          {
            case GeometryType::Spherical:
              outer_dims[0] += sphericalThickness( subMat, x );
              break;
              
            case GeometryType::CylinderEndOn:
            case GeometryType::CylinderSideOn:
              outer_dims[0] += cylindricalRadiusThickness( subMat, x );
              outer_dims[1] += cylindricalLengthThickness( subMat, x );
              break;
              
            case GeometryType::Rectangular:
              outer_dims[0] += rectangularWidthThickness( subMat, x );
              outer_dims[1] += rectangularHeightThickness( subMat, x );
              outer_dims[2] += rectangularDepthThickness( subMat, x );
              break;
              
            case GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( m_geometry )
          
          const shared_ptr<const Material> &material = m_initial_shieldings[subMat].m_material;
          
          
          bool pastDetector = false;
          switch( m_geometry )
          {
            case GeometryType::Spherical:
            case GeometryType::CylinderSideOn:
              pastDetector = (outer_dims[0] > m_distance);
              break;
              
            case GeometryType::CylinderEndOn:
              pastDetector = (outer_dims[1] > m_distance);
              break;
              
            case GeometryType::Rectangular:
              pastDetector = (outer_dims[2] > m_distance);
              break;
              
            case GeometryType::NumGeometryType:
              assert( 0 );
              break;
          }//switch( m_geometry )
          
          if( pastDetector )
            throw runtime_error( "energy_chi_contributions: radius > distance" );
          
          const double transLenCoef = transmition_length_coefficient( material.get(), calculator->m_energy );

#if( defined(__GNUC__) && __GNUC__ < 5 )
          calculator->m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Material} );
#else
          calculator->m_dimensionsTransLenAndType.push_back( {outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Material} );
#endif
        }//for( int subMat = 0; subMat < nMaterials; ++subMat )

        if( calculator->m_dimensionsTransLenAndType.empty() )
          throw std::logic_error( "No source/shielding sphere for calculator" );
        
        const DistributedSrcCalc * const raw_calc = calculator.get();
        
        if( log_info )
        {
          is_trace_src[raw_calc] = is_trace;
          vol_src_material_index[raw_calc] = material_index;
        }
        
        calculators.push_back( std::move(calculator) );
      }//for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
    }//for( const SandiaDecay::Nuclide *src : self_atten_srcs )
  }//for( int material_index = 0; material_index < nMaterials; ++material_index )

  if( calculators.size() )
  {
    if( m_options.multithread_self_atten )
    {
      SpecUtilsAsync::ThreadPool pool;
      for( const unique_ptr<DistributedSrcCalc> &calculator : calculators )
        pool.post( boost::bind( &ShieldingSourceChi2Fcn::selfShieldingIntegration, boost::ref(*calculator) ) );
      pool.join();
    }else
    {
      for( const unique_ptr<DistributedSrcCalc> &calculator : calculators )
        selfShieldingIntegration(*calculator);
    }
    
//    vector<boost::function<void()> > workers;
//    for( DistributedSrcCalc &calculator : calculators )
//    {
//      boost::function<void()> worker = boost::bind( &ShieldingSourceChi2Fcn::selfShieldingIntegration, boost::ref(calculator) );
//      workers.push_back( worker );
//    }//for( size_t i = 0; i < calculators.size(); ++i )
//    SpecUtils::do_asyncronous_work( workers, false );

    if( info )
    {
      string msg;
      if( !has_self_atten )
        msg += "Trace";
      else if( !has_trace )
        msg += "Self Attenuating";
      else
        msg += "Trace and Self Attenuating";
      
      msg += " Source Info (after accounting for shielding, distance, detector and live time):";
      
      info->push_back( std::move(msg) );
    }//if( info )
    
    for( const unique_ptr<DistributedSrcCalc> &calculator : calculators )
    {
      double contrib = calculator->integral * calculator->m_srcVolumetricActivity;

      
      if( m_detector && m_detector->isValid() )
        contrib *= m_detector->intrinsicEfficiency( calculator->m_energy );

      if( energy_count_map.find( calculator->m_energy ) != energy_count_map.end() )
      {
//        cerr << "Adding " << contrib << " to energy " << calculator->m_energy
//             << ", calculator->integral=" << calculator->integral
//             << ", calculator->m_srcVolumetricActivity=" << calculator->m_srcVolumetricActivity
//             << endl;
        energy_count_map[calculator->m_energy] += contrib;
      }else
      {
//        cerr << "Setting " << contrib*m_liveTime << " counts to energy " << calculator->energy
//             << " for thickness=" << calculator->m_dimensionsTransLenAndType[calculator->m_materialIndex].first[0] / PhysicalUnits::cm
//             << " cm" << endl;
        energy_count_map[calculator->m_energy] = contrib;
      }
      
      if( info )
      {
        assert( calculator->m_materialIndex < m_initial_shieldings.size() );
        const shared_ptr<const Material> &material = m_initial_shieldings[calculator->m_materialIndex].m_material;
        const int index = static_cast<int>( calculator->m_materialIndex );
        
        stringstream msg;
        msg << "\tAttributing " << contrib*PhysicalUnits::second << " cps to "
            << calculator->m_energy/PhysicalUnits::keV << " keV photopeak ";
        if( calculator->m_nuclide )
          msg << "(from " << calculator->m_nuclide->symbol << ") ";
        
        msg << "for thicknesses {";
        const array<double,3> &dims = std::get<0>(calculator->m_dimensionsTransLenAndType[index]);
        double dx = dims[0];
        double dy = dims[1];
        double dz = dims[2];
        if( index > 0 )
        {
          const array<double,3> &inner_dims = std::get<0>(calculator->m_dimensionsTransLenAndType[index-1]);
          dx -= inner_dims[0];
          dy -= inner_dims[1];
          dz -= inner_dims[2];
        }
        
        switch( m_geometry )
        {
          case GeometryType::Spherical:
            msg << "rad=" << PhysicalUnits::printToBestLengthUnits(dx);
            break;
            
          case GeometryType::CylinderEndOn:
          case GeometryType::CylinderSideOn:
            msg << "rad=" << PhysicalUnits::printToBestLengthUnits(dx)
                << ", len=" << PhysicalUnits::printToBestLengthUnits(dy);
            break;
          
          case GeometryType::Rectangular:
            msg << "width=" << PhysicalUnits::printToBestLengthUnits(dx)
                << ", height=" << PhysicalUnits::printToBestLengthUnits(dy)
                << ", depth=" << PhysicalUnits::printToBestLengthUnits(dy);
            break;
          case GeometryType::NumGeometryType:
            break;
        }//switch( m_geometry )
        
        msg << "} " << material->name;
        
        info->push_back( msg.str() );
      }//if( info )
      
      if( log_info )
      {
        const double energy = calculator->m_energy;
        auto pos = std::find_if( begin(*log_info), end(*log_info), [energy]( const GammaInteractionCalc::PeakDetail &val ) {
          return energy == val.decayParticleEnergy;
        });
        
        assert( pos != end(*log_info) );
        
        assert( is_trace_src.count(calculator.get()) );
        assert( vol_src_material_index.count(calculator.get()) );
        
        if( pos != end(*log_info) )
        {
          GammaInteractionCalc::PeakDetail &peak = *pos;
         
          // If a volumetric source, we need to update the CPS of the source, to be from per volume, to total
          //  The trick is here we only want to correct each quantity once... which should be the case.
          const double volume = volumeOfMaterial(calculator->m_materialIndex, x);
          for( PeakDetailSrc &src : peak.m_sources )
          {
            const bool isTrace = isTraceSource(src.nuclide);         //note: `src.isTraceSource` is not filled out yet
            const bool isSelfAtten = isSelfAttenSource(src.nuclide);// note: `src.isSelfAttenSource` is not filled out yet
            
            if( (isTrace || isSelfAtten) && (src.nuclide == calculator->m_nuclide) )
            {
              if( !calculator->m_isInSituExponential )
              {
                src.cpsAtSource *= volume;
                src.countsAtSource *= volume;
              }else
              {
                // For in-situ exponential, we are tracking per surface area
                const size_t mat_index = calculator->m_materialIndex;
                const double L = relaxationLength(src.nuclide);
                
                switch( m_geometry )
                {
                  case GeometryType::Spherical:
                  {
                    // TODO: The same approach as for GeometryType::CylinderEndOn, doesnt seem to work here; ran out
                    //  of time to figure out proper answer, so just setting to zero for the moment
                    src.cpsAtSource = 0;
                    src.countsAtSource = 0;
                    //const double R = sphericalThickness(mat_index, x);
                    //const double norm = 4*PhysicalUnits::pi * L * (L*L*(2 - 2*exp(-R/L)) - 2*L*R + R*R);
                    //const double sa = 4.0*PhysicalUnits::pi*R*R;
                    //src.cpsAtSource *= (sa * norm);
                    //src.countsAtSource *= (sa * norm);
                    break;
                  }//case GeometryType::Spherical:
                    
                  case GeometryType::CylinderSideOn:
                  {
                    // TODO: The same approach as for GeometryType::CylinderEndOn, doesnt seem to work here; ran out
                    //  of time to figure out proper answer, so just setting to zero for the moment
                    src.cpsAtSource = 0;
                    src.countsAtSource = 0;
                    //const double R = cylindricalRadiusThickness(mat_index, x);
                    //const double norm = 2 * L * PhysicalUnits::pi * (L*(exp(-R/L) - 1) + R);
                    //const double h = 2.0*cylindricalLengthThickness(mat_index, x);
                    //const double sa = h * 2.0*R*PhysicalUnits::pi;
                    //src.cpsAtSource *= (sa * norm);
                    //src.countsAtSource *= (sa * norm);
                    break;
                  }//case GeometryType::CylinderSideOn:
                  
                  case GeometryType::CylinderEndOn:
                  {
                    const double R = 2.0 * cylindricalLengthThickness(mat_index, x);
                    const double norm = L * (1.0 - exp(-R / L) );
                    const double r = cylindricalRadiusThickness(mat_index, x);
                    const double sa = PhysicalUnits::pi*r*r;
                    src.cpsAtSource *= (sa * norm);
                    src.countsAtSource *= (sa * norm);
                    break;
                  }//case GeometryType::CylinderEndOn:
                  
                  case GeometryType::Rectangular:
                  {
                    // TODO: The same approach as for GeometryType::CylinderEndOn, doesnt seem to work here; ran out
                    //  of time to figure out proper answer, so just setting to zero for the moment
                    src.cpsAtSource = 0;
                    src.countsAtSource = 0;
                    //const double R = 2.0 * rectangularDepthThickness(mat_index, x);
                    //const double norm = L * (1.0 - exp(-R / L) );
                    //const double w = rectangularWidthThickness(mat_index, x);
                    //const double h = rectangularHeightThickness(mat_index, x);
                    //const double sa = w*h;
                    //src.cpsAtSource *= (sa * norm);
                    //src.countsAtSource *= (sa * norm);
                    break;
                  }
                    
                  case GeometryType::NumGeometryType:
                    assert( 0 );
                    break;
                }//switch( m_geometry )
              }//if( !calculator->m_isInSituExponential ) / else
            }//if( (src.isTraceSource || src.isSelfAttenSource) && (src.nuclide == calculator->m_nuclide) )
          }//for( PeakDetailSrc &src : peak.m_sources )
          
          PeakDetail::VolumeSrc src;
          src.trace = is_trace_src[calculator.get()];
          src.integral = calculator->integral;
          src.volume = volume;
          src.averageEfficiencyPerSourceGamma = src.integral / src.volume;
          src.srcVolumetricActivity = calculator->m_srcVolumetricActivity;
          src.inSituExponential = calculator->m_isInSituExponential;
          src.inSituRelaxationLength = calculator->m_inSituRelaxationLength;
          src.detIntrinsicEff = 1.0;
          if( m_detector && m_detector->isValid() )
            src.detIntrinsicEff = m_detector->intrinsicEfficiency( calculator->m_energy );
          src.sourceName = calculator->m_nuclide ? calculator->m_nuclide->symbol : string("null");
          
          
          peak.m_airAttenFactor = 1.0;
          if( m_options.attenuate_for_air )
          {
            const double coef = transmission_length_coefficient_air( peak.energy );
            peak.m_airAttenFactor = exp( -1.0 * coef * air_dist );
          }
              
          double det_intrinsic = 1.0, det_total_eff = 1.0;
          if( m_detector && m_detector->isValid() )
          {
            det_intrinsic = m_detector->intrinsicEfficiency( peak.energy );
            det_total_eff = m_detector->isFixedGeometry() ? det_intrinsic : m_detector->efficiency(peak.energy, m_distance);
          }
          const double geom_factor = det_total_eff / det_intrinsic;
          peak.detIntrinsicEff = det_intrinsic;
          peak.detEff = det_total_eff;
          peak.detSolidAngle = det_total_eff / det_intrinsic;
          peak.m_totalAttenFactor = src.averageEfficiencyPerSourceGamma / geom_factor;
          peak.m_totalShieldAttenFactor = peak.m_totalAttenFactor / peak.m_airAttenFactor;
          
          peak.m_volumetric_srcs.push_back( std::move(src) );
        }//if( pos != end(*log_info) )
      }//if( log_info )
    }//for( DistributedSrcCalc &calculator : calculators )
  }//if( calculators.size() )

  
  //Account for live time
  for( EnergyCountMap::value_type &energy_count : energy_count_map )
    energy_count.second *= m_liveTime;

  if( log_info )
  {
    const size_t nnucs = numNuclides();
    for( size_t i = 0; i < nnucs; ++i )
    {
      const SandiaDecay::Nuclide * const nuc = nuclide(i);
      assert( nuc );
      if( !nuc )
        continue;
      
      const double activityVal = activity(nuc, x);
      double activityUncertVal = 0.0, massFractionVal = 0.0, massFractionUncertVal = 0.0;
      double ageUncert = 0.0;
      
      if( !error_params.empty() )
      {
        ageUncert = age(nuc, error_params);
        activityUncertVal = activityUncertainty(nuc, x, error_params);
      }
      
      if( isSelfAttenSource(nuc) )
      {
        for( size_t shielding_index = 0; shielding_index < numMaterials(); ++shielding_index )
        {
          const Material * const mat = material(shielding_index);
          if( !mat )
            continue;
          
          const vector<const SandiaDecay::Nuclide *> nucs = selfAttenuatingNuclides(shielding_index);
          const auto pos = std::find( begin(nucs), end(nucs), nuc );
          if( pos != end(nucs) )
          {
            massFractionVal = massFractionOfElement(shielding_index, nuc, x);
            if( !error_params.empty() )
              massFractionUncertVal = massFractionOfElementUncertainty(shielding_index, nuc, x, error_params);
          }
        }//for( loop over shielding_index )
      }//if( isSelfAtten )
        
      
      bool foundNuc = false;
      for( GammaInteractionCalc::PeakDetail &p : *log_info )
      {
        for( PeakDetailSrc &psrc : p.m_sources )
        {
          if( psrc.nuclide == nuc )
          {
            foundNuc = true;
            
            psrc.isTraceSource = isTraceSource(nuc);
            if( psrc.isTraceSource )
              psrc.traceSourceType = traceSourceActivityType(nuc);
            psrc.isSelfAttenSource = isSelfAttenSource(nuc);
            //double ShieldingSourceChi2Fcn::relaxationLength(nuc);
            for( size_t i = 0; i < numMaterials(); ++i )
              psrc.isFittingMassFraction |= isVariableMassFraction( i, nuc );
            
            psrc.countsAtSource = psrc.cpsAtSource * m_liveTime;
            //psrc.countsUncert = countsUncert
            psrc.ageUncert = ageUncert;
            //psrc.ageIsFit = ageIsFit
            //psrc.canFitAge = canFitAge
            psrc.activity = activityVal;
            psrc.activityUncert = activityUncertVal;
            psrc.displayActivity = psrc.isTraceSource ? totalActivity(nuc,x) : activityVal;
            psrc.displayActivityUncert = (psrc.isTraceSource && !error_params.empty())
                                          ? totalActivityUncertainty(nuc,x,error_params)
                                          : activityUncertVal;
            psrc.massFraction = massFractionVal;
            psrc.massFractionUncert = massFractionUncertVal;
          }//if( psrc.nuclide == nuc )
        }//for( PeakDetailSrc &psrc : p.m_sources )
      }//for( GammaInteractionCalc::PeakDetail &p : *log_info )
      
      assert( foundNuc );
    }//for( size_t i = 0; i < nnucs; ++i )
  }//if( log_info )
  
  
  return expected_observed_chis( m_peaks, m_backgroundPeaks, energy_count_map, info, log_info );
}//vector<PeakResultPlotInfo> energy_chi_contributions(...) const

  
void ShieldingSourceChi2Fcn::log_shield_info( const vector<double> &params,
                                              const vector<double> &errors,
                              const vector<ShieldingSourceFitCalc::IsoFitStruct> &fit_src_info,
                              vector<ShieldingDetails> &shielding_details ) const
{
  const ShieldingSourceChi2Fcn * const chi2Fcn = this;
  
  const shared_ptr<const DetectorPeakResponse> &det = chi2Fcn->detector();
  
  const DetectorPeakResponse::EffGeometryType detType = (det && det->isValid())
                                                  ? det->geometryType()
                                                  : DetectorPeakResponse::EffGeometryType::FarField;
  try
  {
    double shield_outer_rad = 0.0;
    double outer_dims[3] = { 0.0, 0.0, 0.0 };  //Only filled out `if( log_info != nullptr )`
    
    const GammaInteractionCalc::GeometryType geometry = chi2Fcn->geometry();
    const size_t nMaterials = chi2Fcn->numMaterials();
    
    for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
    {
      ShieldingDetails shieldInfo;
      shieldInfo.m_geometry = geometry;
      shieldInfo.m_inner_rad = shield_outer_rad;
      shieldInfo.m_is_generic = chi2Fcn->isGenericMaterial(shielding_index);
      
      if( shieldInfo.m_is_generic )
      {
        const float atomic_number = static_cast<float>(chi2Fcn->atomicNumber( shielding_index, params ));
        const float areal_density = static_cast<float>(chi2Fcn->arealDensity( shielding_index, params ));
        const double ad_in_gcm2 = areal_density * PhysicalUnits::cm2 / PhysicalUnits::g;
        //auto att_coef_fcn = boost::bind( &transmition_coefficient_generic, atomic_number, areal_density, boost::placeholders::_1 );
        
        shieldInfo.m_name = "Generic (an=" + SpecUtils::printCompact(atomic_number,3)
                      + ", ad=" + SpecUtils::printCompact(ad_in_gcm2,4) + ")";
        shieldInfo.m_ad = areal_density;
        shieldInfo.m_an = atomic_number;
        
        shieldInfo.m_thickness = 0.0;
        shieldInfo.m_volume = 0.0;
        shieldInfo.m_volume_uncert = 0.0;
        shieldInfo.m_num_dimensions = 0;
        for( size_t i = 0; i < 3; ++i )
        {
          shieldInfo.m_inner_dimensions[i] = outer_dims[i];
          shieldInfo.m_outer_dimensions[i] = outer_dims[i];
        }
      }else
      {
        const Material * const material = chi2Fcn->material( shielding_index );
        assert( material );
        if( !material )
          throw std::logic_error( "Unexpected null material" );
        
        double thickness = 0.0;
        double dim_thickness[3] = { 0.0, 0.0, 0.0 }; //Only filled out `if( log_info != nullptr )`
        switch( geometry )
        {
          case GammaInteractionCalc::GeometryType::Spherical:
            thickness = chi2Fcn->sphericalThickness(shielding_index, params);
            dim_thickness[0] = thickness;
            shieldInfo.m_num_dimensions = 1;
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderEndOn:
            thickness = chi2Fcn->cylindricalLengthThickness(shielding_index, params);
            dim_thickness[0] = chi2Fcn->cylindricalRadiusThickness(shielding_index, params);
            dim_thickness[1] = thickness;
            shieldInfo.m_num_dimensions = 2;
            break;
            
          case GammaInteractionCalc::GeometryType::CylinderSideOn:
            thickness = chi2Fcn->cylindricalRadiusThickness(shielding_index, params);
            dim_thickness[0] = thickness;
            dim_thickness[1] = chi2Fcn->cylindricalLengthThickness(shielding_index, params);
            shieldInfo.m_num_dimensions = 2;
            break;
            
          case GammaInteractionCalc::GeometryType::Rectangular:
            thickness = chi2Fcn->rectangularDepthThickness(shielding_index, params);
            dim_thickness[0] = chi2Fcn->rectangularWidthThickness(shielding_index, params);
            dim_thickness[1] = chi2Fcn->rectangularHeightThickness(shielding_index, params);
            dim_thickness[2] = thickness;
            shieldInfo.m_num_dimensions = 3;
            break;
            
          case GammaInteractionCalc::GeometryType::NumGeometryType:
            assert( 0 );
            break;
        }//switch( m_geometry )
        
        
        shieldInfo.m_name = material ? material->name : string("null");
        shieldInfo.m_chemical_formula = material ? material->chemicalFormula() : string("null");
        shieldInfo.m_density = material ? material->density : 0.0f;
        
        shieldInfo.m_ad = (shieldInfo.m_density * thickness);
        shieldInfo.m_an = material ? material->massWeightedAtomicNumber() : 0.0f;
        shieldInfo.m_thickness = thickness;
        
        for( size_t i = 0; i < 3; ++i )
        {
          shieldInfo.m_inner_dimensions[i] = outer_dims[i];
          outer_dims[i] += dim_thickness[i];
          shieldInfo.m_outer_dimensions[i] = outer_dims[i];
        }
        
        shield_outer_rad += thickness;
        shieldInfo.m_volume = chi2Fcn->volumeOfMaterial( shielding_index, params );
        shieldInfo.m_volume_uncert = chi2Fcn->volumeUncertaintyOfMaterial( static_cast<int>(shielding_index), params, errors );
        
        const vector<const SandiaDecay::Nuclide *> self_atten_nucs = chi2Fcn->selfAttenuatingNuclides( shielding_index );
        
        map<const SandiaDecay::Element *,vector<const SandiaDecay::Nuclide *>> fit_frac_el_to_nucs
            = chi2Fcn->nuclideFittingMassFracFor( shielding_index );
        vector<const SandiaDecay::Nuclide *> fit_frac_nucs;
        for( const auto &el_nucs : fit_frac_el_to_nucs )
        {
          for( const SandiaDecay::Nuclide *nuc : el_nucs.second )
          {
            if( !nuc )
              continue;
            
            fit_frac_nucs.push_back( nuc );
            assert( std::find(begin(self_atten_nucs), end(self_atten_nucs), nuc) != end(self_atten_nucs) );
          }//for( const const SandiaDecay::Nuclide *nuc : el_nucs.second )
        }//for( const auto &el_nucs : fit_frac_el_to_nucs )
        
        
        for( const SandiaDecay::Nuclide *nuc : self_atten_nucs )
        {
          assert( chi2Fcn->isVolumetricSource(nuc) );
          assert( chi2Fcn->isSelfAttenSource(nuc) );
          assert( !chi2Fcn->isTraceSource(nuc) );
          
          ShieldingDetails::SelfAttenComponent comp;
          comp.m_nuclide = nuc;

#ifndef NDEBUG
          const auto src_pos = find_if( begin(fit_src_info), end(fit_src_info),
                  [nuc]( const ShieldingSourceFitCalc::IsoFitStruct &info ){
            return info.nuclide == nuc;
          });
          assert( src_pos != end(fit_src_info) );
#endif
          
          const auto fit_pos = std::find(begin(fit_frac_nucs), end(fit_frac_nucs), nuc);
          comp.m_is_fit = (fit_pos != end(fit_frac_nucs));
          
          comp.m_mass_frac = chi2Fcn->massFractionOfElement(shielding_index, nuc, params);
          comp.m_mass_frac_uncert = chi2Fcn->massFractionOfElementUncertainty(shielding_index, nuc, params, errors );
          
          //chi2Fcn->activityOfSelfAttenSource( nuc, params );
          
          shieldInfo.m_mass_fractions.push_back( comp );
        }//for( const SandiaDecay::Nuclide *nuc : self_atten_nucs )
        
        const vector<const SandiaDecay::Nuclide *> trace_srcs = chi2Fcn->traceNuclidesForMaterial( shielding_index );
        
        for( const SandiaDecay::Nuclide *nuc : trace_srcs )
        {
          assert( chi2Fcn->isVolumetricSource(nuc) );
          assert( chi2Fcn->isTraceSource(nuc) );
          assert( !chi2Fcn->isSelfAttenSource(nuc) );
        
          ShieldingDetails::TraceSrcDetail src;
          src.m_nuclide = nuc;
          //src.m_age = chi2Fcn->age( nuc, params );
          src.m_trace_type = chi2Fcn->traceSourceActivityType( nuc );
          src.m_is_exp_dist = (src.m_trace_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution);
          
#ifndef NDEBUG
          const auto src_pos = find_if( begin(fit_src_info), end(fit_src_info),
                  [nuc]( const ShieldingSourceFitCalc::IsoFitStruct &info ){
            return info.nuclide == nuc;
          });
          assert( src_pos != end(fit_src_info) );
#endif
          
          if( src.m_is_exp_dist )
            src.m_relaxation_length = chi2Fcn->relaxationLength( nuc );
          
          shieldInfo.m_trace_sources.push_back( src );
        }//for( const SandiaDecay::Nuclide *nuc : trace_srcs )
      }//if( shieldInfo.m_is_generic ) / else
      
      shielding_details.push_back( std::move(shieldInfo) );
    }//for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
  
  }catch( std::exception &e )
  {
    throw std::runtime_error( "There was an error and log info about shielding is not be complete: "
                         + string(e.what()) );
  }
}//log_shield_info(...)
  
  
void ShieldingSourceChi2Fcn::log_source_info( const std::vector<double> &params,
                        const std::vector<double> &errors,
                        const vector<ShieldingSourceFitCalc::IsoFitStruct> &fit_src_info,
                        std::vector<SourceDetails> &info ) const
{
  const ShieldingSourceChi2Fcn * const chi2Fcn = this;
  
  try
  {
    info.clear();
    
    const size_t nnuc = chi2Fcn->numNuclides();
    for( size_t nucn = 0; nucn < nnuc; ++nucn )
    {
      const SandiaDecay::Nuclide *nuc = chi2Fcn->nuclide( nucn );
      assert( nuc );
      if( !nuc )
        continue;
      
      const auto pos = std::find_if( begin(fit_src_info), end(fit_src_info), [nuc]( const ShieldingSourceFitCalc::IsoFitStruct &iso ) {
        return iso.nuclide == nuc;
      });
      
      assert( pos != end(fit_src_info) );
      if( pos == end(fit_src_info) )
        throw runtime_error( "Missing ShieldingSourceFitCalc::IsoFitStruct for " + nuc->symbol );
      
      const ShieldingSourceFitCalc::IsoFitStruct &fit_info = *pos;
      
      SourceDetails src;
      src.nuclide = nuc;
      src.activity = chi2Fcn->activity( nuc, params );
      src.activityUncertainty = chi2Fcn->activityUncertainty( nuc, params, errors );
      src.activityIsFit = fit_info.fitActivity;
      src.nuclideMass = (src.activity / nuc->activityPerGram()) * PhysicalUnits::gram;
      src.age = chi2Fcn->age( nuc, params );
      src.ageUncertainty = chi2Fcn->age( nuc, errors );;
      src.ageIsFittable = fit_info.ageIsFittable;
      src.ageIsFit = fit_info.fitAge;
      src.ageDefiningNuc = fit_info.ageDefiningNuc;
      src.isTraceSource = chi2Fcn->isTraceSource(nuc);
      if( src.isTraceSource )
      {
        src.traceActivityType = chi2Fcn->traceSourceActivityType(nuc);
        src.traceSrcDisplayAct = chi2Fcn->totalActivity(nuc,params);
        src.traceSrcDisplayActUncertainty = chi2Fcn->totalActivityUncertainty(nuc, params, errors );
        src.traceRelaxationLength = 0.0;
        
        TraceActivityType traceType = chi2Fcn->traceSourceActivityType(nuc);
        if( traceType == TraceActivityType::ExponentialDistribution )
          src.traceRelaxationLength = chi2Fcn->relaxationLength(nuc);
      }//
      
      src.isSelfAttenSource = chi2Fcn->isSelfAttenSource( nuc );
      if( src.isSelfAttenSource )
      {
        bool found_shield = false;
        for( size_t shielding_index = 0; shielding_index < chi2Fcn->numMaterials(); ++shielding_index )
        {
          const Material * const mat = chi2Fcn->material(shielding_index);
          if( !mat )
            continue;
          
          const vector<const SandiaDecay::Nuclide *> nucs = chi2Fcn->selfAttenuatingNuclides(shielding_index);
          const auto pos = std::find( begin(nucs), end(nucs), nuc );
          if( pos != end(nucs) )
          {
            found_shield = true;
            
            src.selfAttenShieldIndex = shielding_index;
            src.selfAttenShieldName = mat->name;
            
            src.isSelfAttenVariableMassFrac = chi2Fcn->isVariableMassFraction(shielding_index, nuc);
            src.selfAttenMassFrac = chi2Fcn->massFractionOfElement(shielding_index, nuc, params);
            src.selfAttenMassFracUncertainty = chi2Fcn->massFractionOfElementUncertainty(shielding_index, nuc, params, errors);
          }
        }//for( loop over shielding_index )
        
        assert( found_shield );
      }//if( src.isSelfAttenSource )
      
      info.push_back( std::move(src) );
    }//for( size_t nucn = 0; nucn < nnuc; ++nucn )
  }catch( std::exception &e )
  {
    throw runtime_error( "There was an error and source information log is not complete: "
                                  + string(e.what()) );
  }//try / catch
}//void log_source_info(...)
  
  
ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::GuiProgressUpdateInfo( const size_t updateFreqMs,
                        std::function<void(size_t, double, double, std::vector<double>)> updater )
  : m_gui_updater( updater ),
  m_update_frequency_ms( updateFreqMs ),
  m_num_fcn_calls( 0 ),
  m_fitStartTime( 0.0 ),
  m_currentTime( 0.0 ),
  m_lastGuiUpdateTime( 0.0 ),
  m_bestChi2( std::numeric_limits<double>::infinity() )
{
}


void ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::fitting_starting()
{
  std::lock_guard<std::mutex> scoped_lock( m_mutex );

  m_bestChi2 = std::numeric_limits<double>::max();
  m_num_fcn_calls = 0;
  m_bestParameters.clear();
  m_lastGuiUpdateTime = SpecUtils::get_wall_time();
  m_fitStartTime = m_lastGuiUpdateTime;
  m_currentTime = m_lastGuiUpdateTime;
}// void GuiProgressUpdateInfo::fitting_starting()


void ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::completed_eval( const double chi2, const std::vector<double> &pars )
{
  std::lock_guard<std::mutex> scoped_lock( m_mutex );
  
  m_num_fcn_calls += 1;
  const double prevBestChi2 = m_bestChi2;
  if( chi2 < prevBestChi2 )
  {
    m_bestChi2 = chi2;
    m_bestParameters = pars;
  }
  
  m_currentTime = SpecUtils::get_wall_time();
  const double elapsed_time = m_currentTime - m_fitStartTime;
  if( elapsed_time > 0.001*m_update_frequency_ms )
  {
    m_lastGuiUpdateTime = m_currentTime;
    m_gui_updater( m_num_fcn_calls, elapsed_time, m_bestChi2, m_bestParameters );
  }
}//completed_eval(...)


size_t ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::numFunctionCallsSoFar()
{
  std::lock_guard<std::mutex> scoped_lock( m_mutex );
  return m_num_fcn_calls;
}

double ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::bestChi2SoFar()
{
  std::lock_guard<std::mutex> scoped_lock( m_mutex );
  return m_bestChi2;
}

std::vector<double> ShieldingSourceChi2Fcn::GuiProgressUpdateInfo::bestParametersSoFar()
{
  std::lock_guard<std::mutex> scoped_lock( m_mutex );
  return m_bestParameters;
}

}//namespace GammaInteractionCalc
