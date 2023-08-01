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
    const double xs_per_mass = MassAttenuation::massAttenuationCoeficient( atomicNumber, energy );
    
    mu += (density*xs_per_mass);
  }//for( const Material::ElementFractionPair &p : material->elements )

  for( const Material::NuclideFractionPair &p : material->nuclides )
  {
    const int atomicNumber = p.first->atomicNumber;

    if( atomicNumber > MassAttenuation::sm_max_xs_atomic_number )
      throw std::runtime_error( "transmition_length_coefficient(...): invalid "
                                "atomic number" );
    double density = p.second * material->density;
    const double xs_per_mass = MassAttenuation::massAttenuationCoeficient( atomicNumber, energy );
    
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
  //const double mu = MassAttenuation::massAttenuationCoeficient( air_an, energy );
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
    const double xs_per_mass = MassAttenuation::massAttenuationCoeficient( atomicNumber, energy );
    
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
  const double xs_per_mass = MassAttenuation::massAttenuationCoeficient( atomic_number, energy );
  
  return xs_per_mass;
}

double transmition_coefficient_generic( float atomic_number, float areal_density,
                                float energy )
{
  return areal_density * mass_attenuation_coef( atomic_number, energy );
}


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


void example_integration()
{
  
  //cout << "Done in example_integration()" << endl;
  //return;
  
  
  /*
   //  double point[3], x1[3], x2[3];
   //  point[0] = point[1] = point[2] = 0.0;
   //  x1[0] = x1[1] = x1[2] = 0.0;
   //  x2[0] = x2[1] = x2[2] = 0.0;
   
   //  x1[0] = -1;
   //  x1[1] = 5;
   
   //  x2[0] = 1;
   //  x2[1] = 1;
   //  cerr << GammaInteractionCalc::point_to_line_dist( point, x1, x2 ) << endl;
   //  return 1;

   
    double x0=0.0, y0=-0.5, z0 =-0.5;
    double x=x0, y=y0, z=z0;
    double sphere_rad = 1.0, observation_dist = 4.0;

    const double pi = 3.14159265358979;
    double r = 0.5;
    double theta = pi/2.0;
    double phi = 0.0;
    x0 = r * sin( theta ) * cos( phi );
    y0 = r * sin( theta ) * sin( phi );
    z0 = r * cos( theta );
    x=x0; y=y0; z=z0;
    exit_point_of_sphere_z( x0, y0, z0, sphere_rad, observation_dist );
    cerr << "Exit at x0=" << x0 << ", y0=" << y0 << ", z0=" << z0 << endl << endl;

    phi = pi/2.0;
    x0 = r * sin( theta ) * cos( phi );
    y0 = r * sin( theta ) * sin( phi );
    z0 = r * cos( theta );
    x=x0; y=y0; z=z0;
    exit_point_of_sphere_z( x0, y0, z0, sphere_rad, observation_dist );
    cerr << "Exit at x0=" << x0 << ", y0=" << y0 << ", z0=" << z0 << endl << endl;



    phi = pi/2.0 + 0.1;
    theta = pi/2.0 + 0.1;
    x0 = r * sin( theta ) * cos( phi );
    y0 = r * sin( theta ) * sin( phi );
    z0 = r * cos( theta );
    x=x0; y=y0; z=z0;
    exit_point_of_sphere_z( x0, y0, z0, sphere_rad, observation_dist );
    cerr << "Exit at x0=" << x0 << ", y0=" << y0 << ", z0=" << z0 << endl << endl;

    exit(1);
  */

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    MaterialDB materialdb;
    materialdb.parseGadrasMaterialFile( "../data/MaterialDataBase.txt", db, false );
    //const Material *material = materialdb.material( "Air" );
    DistributedSrcCalc ObjectToIntegrate;
  
    const double energy = 185.0*PhysicalUnits::keV;
    ObjectToIntegrate.m_geometry = GeometryType::Spherical;
    ObjectToIntegrate.m_sourceIndex = 0;
    ObjectToIntegrate.m_attenuateForAir = false;
    ObjectToIntegrate.m_airTransLenCoef = 0.0;
    ObjectToIntegrate.m_isInSituExponential = false;
    ObjectToIntegrate.m_inSituRelaxationLength = 0.0;
    ObjectToIntegrate.m_detectorRadius  = 2.0 * PhysicalUnits::cm;
    ObjectToIntegrate.m_observationDist = 400.0 * PhysicalUnits::cm;

    double sphereRad = 0.0, transLenCoef = 0.0;

    const Material *material = materialdb.material( "void" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 99.5* PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
    ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif
  
    material = materialdb.material( "U" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 0.5 * PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
    ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif
    ObjectToIntegrate.m_sourceIndex = ObjectToIntegrate.m_dimensionsTransLenAndType.size() - 1;

    material = materialdb.material( "Fe" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 0.5 * PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
  ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#else
    ObjectToIntegrate.m_dimensionsTransLenAndType.push_back( {{sphereRad,0.0,0.0},transLenCoef,DistributedSrcCalc::ShellType::Material} );
#endif

//Original implementation (pre 20140816) with call directly to Cuhre library
//  gave:
//  ndim=2 CUHRE RESULT:	nregions 18	neval 2275	fail 0
//  CUHRE RESULT:	1.83917814 +- 0.00001307	p = 0.000
 
//  ndim=3 CUHRE RESULT:	nregions 23	neval 5715	fail 0
//  CUHRE RESULT:	1.83917736 +- 0.00001600	p = 0.000
  
  
    int ndim = 2;  //the number of dimensions of the integral.
    void *userdata = (void *)&ObjectToIntegrate;
    const double epsrel = 1e-5;  //the requested relative accuracy
    const double epsabs = -1.0;//1e-12; //the requested absolute accuracy
    const int mineval = 0; //the minimum number of integrand evaluations required.
    const int maxeval = 5000000; //the (approximate) maximum number of integrand evaluations allowed.

    int nregions, neval, fail;
    double integral, error, prob;

    ndim = 2;
  Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                            Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                            fail, integral, error, prob );

  printf("ndim=%d CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
      ndim, nregions, neval, fail);
  printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob);
  printf("\n\n" );
  
  ndim = 3;
  Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_spherical, userdata, epsrel, epsabs,
                             Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                            fail, integral, error, prob );
  
    printf("ndim=%d CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
       ndim, nregions, neval, fail);
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob );
    cout << endl << endl;

  /*
    ndim = 2;
  #define SEED 0
  #define NNEW 1000
  #define FLATNESS 25.
    Suave(ndim, ncomp, DistributedSrcCalc_integrand_spherical, userdata,
      epsrel, epsabs, verbose | last, SEED,
      mineval, maxeval, NNEW, FLATNESS,
      &nregions, &neval, &fail, integral, error, prob);
    printf("ndim=%d SAUVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
       ndim, nregions, neval, fail);
    for( comp = 0; comp < ncomp; ++comp )
      printf("SAUVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
        integral[comp], error[comp], prob[comp]);
  */

}//void example_integration()
  

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
  m_sourceIndex = 0;
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
  assert( m_sourceIndex < m_dimensionsTransLenAndType.size() );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_sourceIndex]) == ShellType::Material );
  
  const double pi = PhysicalUnits::pi;
  
  const int ndim = (ndimptr ? (*ndimptr) : 3);

  const double source_inner_rad = ((m_sourceIndex>0)
                            ? std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex-1])[0]
                            : 0.0);
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex]);
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
    const double srcRad = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex])[0];
    const double srcTransCoef = std::get<1>(m_dimensionsTransLenAndType[m_sourceIndex]);
    double dist_in_src = exit_point_of_sphere_z( source_point, exit_point,
                                                 srcRad, m_observationDist );
    
    double min_rad = 0.0;
    bool needShellCompute = false;
    double inner_shell_point[3] = { 0.0 };
    
    if( m_sourceIndex > 0 )
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
      
      const array<double,3> &sub_dims = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex-1]);
      const double subrad = sub_dims[0];
      needShellCompute = ( (min_rad < subrad) && (inner_shell_point[2] > source_point[2]) );
      // Note that if this shell is a generic shell, and the ray *just* touches it, then we wont
      //  include the generic attenuation.
    }//if( m_sourceIndex > 0 )
    
    if( needShellCompute )
    {
      double original_point[3];
      memcpy( original_point, source_point, 3*sizeof(double) );
      memcpy( source_point, inner_shell_point, 3*sizeof(double) );
      
      //Calc how far from the gammas original position it was, to the first
      //  inner shell it hit
      const array<double,3> &sub_dims = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex-1]);
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
      assert( m_sourceIndex != 0 );
      assert( start_index < m_sourceIndex );
      assert( start_index < m_dimensionsTransLenAndType.size() );
        
      if( start_index == m_dimensionsTransLenAndType.size() )
        throw runtime_error( "Logic error 1 in DistributedSrcCalc::eval(...)" );
      if( start_index >= m_sourceIndex )
        throw runtime_error( "Logic error 2 in DistributedSrcCalc::eval(...)" );
      if( m_sourceIndex == 0 )
        throw runtime_error( "Logic error 3 in DistributedSrcCalc::eval(...)" );
      
      //calc distance it travels through the inner spheres
      for( size_t index = start_index; index <= m_sourceIndex; ++index )
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
            if( index != m_sourceIndex )
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
        
        
      }//for( ++sphere_index; sphere_index < m_sourceIndex; ++sphere_index )
    }else
    {
      memcpy( source_point, exit_point, 3*sizeof(double) );
      trans += (srcTransCoef * dist_in_src);
    }//if( line actually goes into child sphere ) / else
  }//end codeblock to compute distance through source
  
  for( size_t i = m_sourceIndex+1; i < m_dimensionsTransLenAndType.size(); ++i )
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

  assert( m_sourceIndex == 0 );
  assert( m_geometry == GeometryType::CylinderEndOn );
  assert( m_dimensionsTransLenAndType.size() == 1 );
  assert( std::get<2>(m_dimensionsTransLenAndType[0]) == ShellType::Material );
  
  //if( m_dimensionsTransLenAndType.size() != 1 )
  //  throw runtime_error( "eval_single_cyl_end_on only supports a single shielding" );
  
  
  // This just integrates a right circular cylinder
  const int ndim = (ndimptr ? (*ndimptr) : 2);
  assert( ndim == 2 ); // ndim==3 is also valid and this function will work for that as well
  
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex]);
  const double source_outer_rad = dimensions[0];
  const double source_half_z = dimensions[1];
  const double total_height = 2.0 * source_half_z;
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_sourceIndex]);
  
  
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
  assert( fabs(test_ff[0] - this_answer) < 1.0E-9*std::max(0.001,std::max( fabs(test_ff[0]), fabs(this_answer))) );
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


#if( DEBUG_RAYTRACE_CALCS )
void test_cylinder_line_intersection()
{
  // This function doesnt exhaustively test #cylinder_line_intersection, but I think it hits all the
  //  lines of its code, so is decent, although there is probably some permutation or possible test
  //  case left out.
  
  const double oneOverSqrt2 = 1.0 / sqrt(2.0);
  
  double dist;
  double radius = 1.0;
  double half_length = 0.5;
  double source[3], detector[3], exit_point[3], start_point[3];
  
  
  // Start tests of a point in the cylinder, and going to the detector
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 2.0; detector[2] = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.0) < 1.0E-14 );
  assert( fabs(exit_point[1] - 1.0) < 1.0E-14 );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-14 );
  assert( fabs(dist - 1.0) < 1.0E-14 );
  
  
  source[0] = 0.5; source[1] = 0.5; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 2.0; detector[2] = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - oneOverSqrt2) < 1.0E-14 );
  assert( fabs(exit_point[1] - oneOverSqrt2) < 1.0E-14 );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-14 );
  assert( fabs(dist - 0.292893) < 1.0E-5 );
  
  
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = -2.0; detector[1] = 0.0; detector[2] = -2.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - -0.5) < 1.0E-14 );
  assert( fabs(exit_point[1] - 0) < 1.0E-14 );
  assert( fabs(exit_point[2] - -0.5) < 1.0E-14 );
  assert( fabs(dist - oneOverSqrt2) < 1.0E-14 );
  
  
  source[0] = 0.5; source[1] = 0.5; source[2] = 0.0;
  detector[0] = 0.5; detector[1] = 0.5; detector[2] = 1.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.5) < 1.0E-14 );
  assert( fabs(exit_point[1] - 0.5) < 1.0E-14 );
  assert( fabs(exit_point[2] - 0.5) < 1.0E-14 );
  assert( fabs(dist - 0.5) < 1.0E-14 );
  
  
  source[0] = 0.5; source[1] = 0.5; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 1.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.25) < 1.0E-14 );
  assert( fabs(exit_point[1] - 0.25) < 1.0E-14 );
  assert( fabs(exit_point[2] - 0.5) < 1.0E-14 );
  assert( fabs(dist - sqrt(2*0.25*0.25 + 0.5*0.5)) < 1.0E-14 );
  
  
  radius = 1.0;
  half_length = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0; detector[1] = 5; detector[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.0) < 1.0E-14 );
  assert( fabs(exit_point[1] - 1.0) < 1.0E-14 );
  assert( fabs(exit_point[2] - 1.0) < 1.0E-14 );
  assert( fabs(dist - sqrt(2.0)) < 1.0E-14 );
  
  
  radius = 1.0;
  half_length = 1.0;
  source[0] = 0.25; source[1] = 0.25; source[2] = 0.25;
  detector[0] = 0.25; detector[1] = 0.25; detector[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.25) < 1.0E-14 );
  assert( fabs(exit_point[1] - 0.25) < 1.0E-14 );
  assert( fabs(exit_point[2] - 1.0) < 1.0E-14 );
  assert( fabs(dist - 0.75) < 1.0E-14 );
  

  radius = 225000;
  half_length = 1000;
  source[0] = 112500.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2000;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 0.5*source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - half_length) < 1.0E-12 );
  assert( fabs(dist - sqrt(half_length*half_length + 0.25*source[0]*source[0])) < 1.0E-12 );
  
  
  // Start tests of a point outside the cylinder, and going through to the detector on other side
  radius = 1;
  half_length = 1;
  source[0] = 0.0; source[1] = 0.0; source[2] = -2.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - half_length) < 1.0E-12 );
  assert( fabs(dist - 3.0) < 1.0E-12 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] + half_length) < 1.0E-12 );
  assert( fabs(dist - 1.0) < 1.0E-12 );
  
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - half_length) < 1.0E-12 );
  assert( fabs(dist - 2.0) < 1.0E-12 );
  
  
  radius = 1;
  half_length = 1;
  source[0] = 10.0; source[1] = 10.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  radius = 1;
  half_length = 1;
  source[0] = -10.0; source[1] = 0.0; source[2] = -10.0;
  detector[0] = 10.0; detector[1] = 0.0; detector[2] = 10;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 1.0) < 1.0E-12 );
  assert( fabs(dist - sqrt(11.0*11.0 + 11.0*11.0)) < 1.0E-12 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( fabs(exit_point[0] - -1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-12 );
  assert( fabs(exit_point[2] - -1.0) < 1.0E-12 );
  assert( fabs(dist - sqrt(2*81.0)) < 1.0E-12 );
  
  // {-1,0,-1} to {1,0,1}
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 1.0) < 1.0E-12 );
  assert( fabs(dist - sqrt(8.0)) < 1.0E-12 );
  
  
  radius = 1;
  half_length = 1;
  source[0] = -5; source[1] = 5; source[2] = 0.0;
  detector[0] = 5; detector[1] = 4; detector[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( dist == 0.0 );
  
  radius = 1;
  half_length = 1;
  source[0] = 5; source[1] = 5; source[2] = -10.0;
  detector[0] = 5; detector[1] = 5; detector[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( dist == 0.0 );
  
  
  radius = 0.0;
  half_length = 0.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source[0] = 2; source[1] = 2; source[2] = -10.0;
  detector[0] = 2; detector[1] = 3; detector[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( dist == 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source[0] = 0.5; source[1] = 0; source[2] = 2;
  detector[0] = 0.5; detector[1] = 2; detector[2] = 10.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( dist == 0.0 );
  
  
  radius = 1;
  half_length = 1;
  source[0] = 0.0; source[1] = 2; source[2] = -2.0;
  detector[0] = 0.0; detector[1] = 1; detector[2] = 2.0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - source[0]) < 1.0E-12 );
  assert( fabs(exit_point[1] - source[1]) < 1.0E-12 );
  assert( fabs(exit_point[2] - source[2]) < 1.0E-12 );
  assert( dist == 0.0 );
  
  
  // Now do some cases to simulate what will happen when integrating over a volume and there is a
  //  sub-volume
  
  radius = 1;
  half_length = 1;
  source[0] = 1.5; source[1] = 0; source[2] = 0;
  detector[0] = -1.5; detector[1] = 0; detector[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, start_point );
  assert( fabs(start_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(start_point[1] - 0) < 1.0E-12 );
  assert( fabs(start_point[2] - 0) < 1.0E-12 );
  assert( dist == 0.5 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] + 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 0) < 1.0E-12 );
  assert( dist == 2.5 );
  
  dist = cylinder_line_intersection( radius, half_length, start_point, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] + 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 0) < 1.0E-12 );
  assert( fabs(dist - 2.0) < 1.0E-12 );
  
  
  
  
  radius = 1;
  half_length = 1;
  source[0] = -1.5; source[1] = 0; source[2] = 0;
  detector[0] = 1.5; detector[1] = 0; detector[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, start_point );
  assert( fabs(start_point[0] + 1.0) < 1.0E-12 );
  assert( fabs(start_point[1] - 0) < 1.0E-12 );
  assert( fabs(start_point[2] - 0) < 1.0E-12 );
  assert( dist == 0.5 );
  
  dist = cylinder_line_intersection( radius, half_length, start_point, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 0) < 1.0E-12 );
  assert( fabs(dist - 2.0) < 1.0E-12 );
  
  
  // Test on boundaries
  
  // Here we'll make the line pass *just* inside the circle to get around how exact intersections
  //  arent totally worked out yet.
  radius = 1.0;
  half_length = 1;
  source[0] = 1.0; source[1] = -1.0; source[2] = 0;
  detector[0] = (1.0 - 2.0E-6); detector[1] = 1; detector[2] = 0;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 3.0E-6 );
  assert( fabs(exit_point[1] - 0) < 3.0E-6 );
  assert( fabs(exit_point[2] - 0) < 3.0E-6 );
  assert( fabs(dist - 1.0) < 3.0E-5 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 3.0E-6 );
  assert( fabs(exit_point[1] - 0) < 3.0E-6 );
  assert( fabs(exit_point[2] - 0) < 3.0E-6 );
  assert( fabs(dist - 1.0) < 3.0E-5);
  
  dist = cylinder_line_intersection( radius, half_length, exit_point, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 3.0E-6 );
  assert( fabs(exit_point[1] - 0) < 3.0E-6 );
  assert( fabs(exit_point[2] - 0) < 3.0E-6 );
  assert( fabs(dist) < 3.0E-5 );

  
  radius = 1;
  half_length = 1;
  source[0] = 1.0; source[1] = 0; source[2] = -5;
  detector[0] = 1.0; detector[1] = 0; detector[2] = 5;
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0) < 1.0E-12 );
  assert( fabs(exit_point[2] - -1.0) < 1.0E-12 );
  assert( dist == 4.0 );
  
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-12 );
  assert( fabs(exit_point[1] - 0) < 1.0E-12 );
  assert( fabs(exit_point[2] - 1.0) < 1.0E-12 );
  assert( dist == 6.0 );
  
   
  // TODO: need more tests for exactly on boundary, or whatever
  
  // Check case where source is between cylinder and detector, so line segment doesnt actually
  //  intersect the cylinder, although it would if line was infinitely long.
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source[0] = 12.691863768989704; source[1] = 0.050459131221286244; source[2] = 0.0;
  detector[0] = 269.23999999999995; detector[1] = 0; detector[2] = 5;
  double to_enter_dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  double to_exit_dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( to_enter_dist == 0.0 );
  assert( to_exit_dist == 0.0 );
  
  
  
  //radius = 6.3499999999999996;
  //half_length = 106.67999999999999;
  //// source is at radius 6.35002389242
  //source[0] = -2.430014265260013; source[1] = 5.8666714672787954; source[2] = 106.68019720704258;
  //detector[0] = 269.23999999999995; detector[1] = 0; detector[2] = 0;
  //dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  //dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );

  
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source[0] = 4.4901280605345768; source[1] = 4.4901280605345759; source[2] = 0;
  detector[0] = 269.23999999999995; detector[1] = 0; detector[2] = 0;
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( dist == 0.0 );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( dist == 0.0 );
  
  
  radius = 6.3499999999999996;
  half_length = 106.67999999999999;
  source[0] = -2.430014265260013; source[1] = 5.8666714672787954; source[2] = 106.68019720704258;
  detector[0] = 269.23999999999995; detector[1] = 0 ;detector[2] = 0;
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::TowardDetector, exit_point );
  assert( fabs( radius - sqrt(exit_point[0]*exit_point[0] + exit_point[1]*exit_point[1]) ) < 1.0E-9*std::max(1.0,radius) );
  
  dist = cylinder_line_intersection( radius, half_length, source, detector, CylExitDir::AwayFromDetector, exit_point );
  assert( fabs( exit_point[2] - half_length ) < 1.0E-9*std::max(1.0,half_length) ); //exit on end
  
}//void test_cylinder_line_intersection()

#endif //DEBUG_RAYTRACE_CALCS


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
    
    
    y_exit = m*x_exit + c;
    
    
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
  assert( m_sourceIndex < m_dimensionsTransLenAndType.size() );
  assert( (m_geometry == GeometryType::CylinderSideOn)
          || (m_geometry == GeometryType::CylinderEndOn) );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_sourceIndex]) == ShellType::Material );
  
  //assert( m_sourceIndex == 0 );
  //if( m_sourceIndex != 0 )
  //  throw runtime_error( "eval_cylinder currently only supports trace source in inner most shielding shielding" );
  
  // This just integrates a right circular cylinder
  const int ndim = (ndimptr ? (*ndimptr) : 3);
  assert( ((m_geometry == GeometryType::CylinderEndOn) && ((ndim == 2) || (ndim == 3)))
         || ((m_geometry == GeometryType::CylinderSideOn) && (ndim == 3)) );
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex]);
  const double source_outer_rad = dimensions[0];
  const double source_half_z = dimensions[1];
  const double total_height = 2.0 * source_half_z;
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_sourceIndex]);
  
  
  //cuba goes from 0 to one for each dimension, so we have to scale the variables
  //  r:     goes from cylindrical inner radius, to outer radius
  //  theta: goes from 0 to 2pi
  //  z:     goes from the negative half-height to positive half-height
  
  double max_theta = 2.0 * PhysicalUnits::pi;
  
  const double r = xx[0] * source_outer_rad;
  const double theta = ((ndim==3) ? xx[1] : 0.0) * max_theta;
  const double z = total_height * (xx[((ndim==3) ? 2 : 1)] - 0.5);
  
  // If point to evaluate is within inner-cylinder, set the source term to zero and return.
  if( m_sourceIndex > 0 )
  {
    const array<double,3> &inner_dims = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex-1]);
    const double &inner_rad = inner_dims[0];
    const double &inner_half_height = inner_dims[1];
    
    if( (r < inner_rad) && (fabs(z) < inner_half_height) )
    {
      ff[0] = 0.0;
      return;
    }
  }//if( m_sourceIndex > 0 )
  
  
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
  for( size_t i = 0; i < m_sourceIndex; ++i )
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
  }//for( size_t i = 0; i < m_sourceIndex; ++i )
  
  trans += (trans_len_coef * (dist_in_cyl - inner_distance));
  
  // Do transport through outer cylinders
  for( size_t i = m_sourceIndex + 1; i < m_dimensionsTransLenAndType.size(); ++i )
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

#if( DEBUG_RAYTRACE_CALCS )
void test_rectangular_intersections()
{
  bool intersected;
  double half_width, half_height, half_depth, dist_in_shape;
  double source[3], detector[3], exit_point[3], enter_point[3];
  
  /*
  // First we'll test the simple case where we know the ray exits the volume on the plane at
  //  +-half_depth on the z-axis
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 10.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 1.0) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - half_depth) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - sqrt(1.0*1.0 + 0.5*0.5)) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.5) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - half_depth) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - -0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - half_depth) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.5; source[1] = -0.5; source[2] = -1.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = 3.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 2) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.5) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - -0.5) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - half_depth) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.5; detector[1] = -0.5; detector[2] = -2.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 1.0*1.0)) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - -0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - -half_depth) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  // Now test the other cases of the ray exiting the rectangle on an arbitrary face.
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 1.0) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = -2.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 1.0) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - -1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 1.0) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = -2.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - 1.0) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - -1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - 0.0) < 1.0E-9*std::max(half_depth,fabs(exit_point[2])) );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -1.0; source[1] = 0.5; source[2] = -0.5;
  detector[0] = 3.0; detector[1] = 0.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - sqrt(0.25*0.25 + 0.25*0.25 + 2.0*2.0)) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - half_width) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - -0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 1.0; source[1] = -1.0; source[2] = -1.0;
  detector[0] = 0.0; detector[1] = 3.0; detector[2] = 0.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(dist_in_shape - sqrt(0.5*0.5 + 0.5*0.5 + 2.0*2.0)) < 1.0E-9*std::max(1.0,dist_in_shape) );
  assert( fabs(exit_point[0] - 0.5) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - -0.5) < 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -0.5; source[1] = -0.5; source[2] = 0.5;
  detector[0] = 2.5; detector[1] = 1.0; detector[2] = 1.0;
  dist_in_shape = rectangle_exit_location( half_width, half_height, half_depth, source, detector, exit_point );
  assert( fabs(exit_point[0] - 1.0) < 1.0E-9*std::max(1.0,fabs(exit_point[0])) );
  assert( fabs(exit_point[1] - 0.25) < 1.0E-9*std::max(1.0,fabs(exit_point[1])) );
  assert( fabs(exit_point[2] - 0.75) < 1.0E-9*std::max(1.0,fabs(exit_point[2])) );
  
  // TODO: add more general direction test cases for rectangle_exit_location
  
  
  //
   */
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -2.0; source[1] = 2.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 2.0; detector[2] = 0.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  assert( !intersected );
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = -2.0; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 2.0; detector[1] = 0.0; detector[2] = 0.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  assert( intersected );
  assert( enter_point[0] == -1.0 );
  assert( enter_point[1] == 0.0 );
  assert( enter_point[2] == 0.0 );
  assert( exit_point[0] == 1.0 );
  assert( exit_point[1] == 0.0 );
  assert( exit_point[2] == 0.0 );
  
  
  
  half_width = 1.0; half_height = 1.0; half_depth = 1.0;
  source[0] = 0.0; source[1] = 0.0; source[2] = -2.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 2.0;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  assert( intersected );
  assert( enter_point[0] == 0.0 );
  assert( enter_point[1] == 0.0 );
  assert( enter_point[2] == -1.0 );
  assert( exit_point[0] == 0.0 );
  assert( exit_point[1] == 0.0 );
  assert( exit_point[2] == 1.0 );
  
  half_width = 10.0; half_height = 10.0; half_depth = 10.0;
  source[0] = -10; source[1] = 0.0; source[2] = 0.0;
  detector[0] = 0.0; detector[1] = 0.0; detector[2] = 250;
  intersected = rectangle_intersections( half_width, half_height, half_depth,
                                        source, detector, enter_point, exit_point );
  assert( intersected );
  //...
  
  // TODO Add more test cases here
  
  
  {// Begin test integration over box
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    MaterialDB materialdb;
    materialdb.parseGadrasMaterialFile( "../data/MaterialDataBase.txt", db, false );
    const double energy = 185*PhysicalUnits::keV;
    
    //const Material *material = materialdb.material( "Uranium" );
    const Material *material = materialdb.material( "void" );
    assert( material );
    
    
    DistributedSrcCalc calc;
    calc.m_geometry = GeometryType::Rectangular;
    calc.m_sourceIndex = 1;
    calc.m_detectorRadius = 1.0*PhysicalUnits::cm;
    calc.m_observationDist = 25*PhysicalUnits::cm;;
    calc.m_attenuateForAir = false;
    calc.m_airTransLenCoef = transmission_length_coefficient_air( energy );
    calc.m_isInSituExponential = false;
    calc.m_inSituRelaxationLength = 3*PhysicalUnits::cm;
    calc.m_srcVolumetricActivity = 100*PhysicalUnits::bq/PhysicalUnits::cm3;
    calc.m_energy = energy;
    calc.m_nuclide = nullptr;
    
    double half_width = 1*PhysicalUnits::cm;
    double half_height = 1*PhysicalUnits::cm;
    double half_depth = 1*PhysicalUnits::cm;
    double trans_coef = transmition_length_coefficient( material, energy );
#if( defined(__GNUC__) && __GNUC__ < 5 )
    calc.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{half_width,half_height,half_depth},trans_coef,ShellType::Material} );
#else
    calc.m_dimensionsTransLenAndType.push_back( {{half_width,half_height,half_depth},trans_coef,ShellType::Material} );
#endif
    
    half_width = 2*PhysicalUnits::cm;
    half_height = 2*PhysicalUnits::cm;
    half_depth = 2*PhysicalUnits::cm;
#if( defined(__GNUC__) && __GNUC__ < 5 )
    calc.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{{half_width,half_height,half_depth},trans_coef,ShellType::Material} );
#else
    calc.m_dimensionsTransLenAndType.push_back( {{half_width,half_height,half_depth},trans_coef,ShellType::Material} );
#endif
    
    calc.integral = 0.0;
    
    void *userdata = (void *)&calc;
    const double epsrel = 1e-4;  //the requested relative accuracy
    const double epsabs = -1.0;//1e-12; //the requested absolute accuracy
    const int mineval = 0; //the minimum number of integrand evaluations required.
    const int maxeval = 5000000; //the (approximate) maximum number of integrand evaluations allowed.
    
    int ndim = 3;
    int nregions, neval, fail;
    double error, prob;
    
    Integrate::CuhreIntegrate( ndim, DistributedSrcCalc_integrand_rectangular, userdata, epsrel, epsabs,
                              Integrate::LastImportanceFcnt,
                              mineval, maxeval, nregions, neval,
                              fail, calc.integral, error, prob );
    
    printf("Rectangle CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", calc.integral, error, prob);
    printf("\n\n" );
  }// End test integration over box
  
  
  
  cout << "Done in test_rectangular_intersections() - no issues found." << endl;
}//void test_rectangular_intersections()
#endif //#if( DEBUG_RAYTRACE_CALCS )



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
  assert( m_sourceIndex < m_dimensionsTransLenAndType.size() );
  assert( std::get<2>(m_dimensionsTransLenAndType[m_sourceIndex]) == ShellType::Material );
  
  const int ndim = (ndimptr ? (*ndimptr) : 2);
  assert( ndim == 3 );
  
  const std::array<double,3> &dimensions = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex]);
  const double half_width  = dimensions[0];
  const double half_height = dimensions[1];
  const double half_depth  = dimensions[2];
  
  const double trans_len_coef = std::get<1>(m_dimensionsTransLenAndType[m_sourceIndex]);
  
  // Translate the [0,1.0] coordinates from Cuba, to the world coordinates we are integrating over.
  //  (note: we would get the same answer if we integrated over just half the width/height, but it
  //   also ends up taking same amount of evaluations, so we will leave fully integrating over
  //   each dimensions, to not cause problems in the future if we add offsets or whatever)
  const double eval_x = (xx[0] - 0.5) * 2.0 * half_width;
  const double eval_y = (xx[1] - 0.5) * 2.0 * half_height;
  const double eval_z = (xx[2] - 0.5) * 2.0 * half_depth;

  
  // Check to see if [eval_x,eval_y,eval_z] is inside an inner volume, and if so set value to zero
  //  and return
  if( m_sourceIndex > 0 )
  {
    const array<double,3> &dims = std::get<0>(m_dimensionsTransLenAndType[m_sourceIndex-1]);
    if( (fabs(eval_x) < dims[0]) && (fabs(eval_y) < dims[1]) && (fabs(eval_z) < dims[2]) )
    {
      ff[0] = 0.0;
      return;
    }
  }//if( m_sourceIndex > 0 )
  
  
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
  for( size_t i = 0; i < m_sourceIndex; ++i )
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
  }//for( size_t i = 0; i < m_sourceIndex; ++i )
  
  
  trans += (trans_len_coef * (dist_in_src - inner_rect_dist));
  
  
  // Account for additional external shielding's
  for( size_t i = m_sourceIndex + 1; i < m_dimensionsTransLenAndType.size(); ++i )
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
    throw runtime_error( "There are not peaks selected for the fit" );
  
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
    vector<PeakDef> backgroundpeaks;
    for( const auto &p : *background_peaks )
      backgroundpeaks.push_back( *p );
    answer->setBackgroundPeaks( backgroundpeaks, background->live_time() );
  }
  
  
  for( size_t shielding_index = 0; shielding_index < shieldings.size(); ++shielding_index )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = shieldings[shielding_index];
    
    if( shield.m_fitMassFrac )
    {
      //Get the isotopes to fit mass fractions of
      vector<const SandiaDecay::Nuclide *> nucstofit;
      
      double fracmaterial = 0.0;
      vector<pair<const SandiaDecay::Nuclide *,double>> massfracs;
      for( const auto &i : shield.m_nuclideFractions )
      {
        fracmaterial += i.second;
        nucstofit.push_back( i.first );
        massfracs.emplace_back( i.first, i.second );
      }
      
      shared_ptr<const Material> mat = shield.m_material;
      if( !mat )
        throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...)"
                            " serious logic error when fitting for mass frac");
      
      if( fracmaterial <= 0.0 )
      {
        passMessage( "When fitting for mass fractions of source nuclides, the "
                    "sum of the fit for mass fractions equal the sum of the "
                    "initial values, therefore the initial sum of mass "
                    "fractions must be greater than 0.0",
                    WarningWidget::WarningMsgHigh );
        throw runtime_error( "Error fitting mass fraction" );
      }//if( fracmaterial <= 0.0 )
      
      // TODO: we could/should replace the below call to `setNuclidesToFitMassFractionFor`, with setting this directorly to `ShieldingSourceChi2Fcn::ShieldLayerInfo::nucs_to_fit_mass_fraction_for` above - and then adding the paramters here...
      answer->setNuclidesToFitMassFractionFor( shielding_index, nucstofit );
      
      nucstofit = answer->nuclideFittingMassFracFor( shielding_index );
      
      
      vector<std::pair<const SandiaDecay::Nuclide *,double>> orderedmassfracs;
      const size_t nfitnucs = (nucstofit.empty() ? 0 : nucstofit.size()-1);
      for( size_t j = 0; j < nfitnucs; ++j )
      {
        const SandiaDecay::Nuclide *nuc = nucstofit[j];
        for( size_t k = 0; k < nfitnucs; ++k )
          if( massfracs[k].first == nuc )
            orderedmassfracs.push_back( massfracs[k] );
      }//for( size_t j = 0; i < nfitnucs; ++j )
      
      if( nfitnucs != orderedmassfracs.size() )
        throw runtime_error( "nfitnucs != orderedmassfracs.size()" );
      
      double usedmassfrac = 0.0;
      for( size_t j = 0; j < nfitnucs; ++j )
      {
        const std::pair<const SandiaDecay::Nuclide *,double> &nmf = orderedmassfracs[j];
        string name = mat->name + "_" + nmf.first->symbol
        + "_" + std::to_string(shielding_index);
        double val = 0.0;
        const double remaining_frac = (fracmaterial - usedmassfrac);
        if( remaining_frac > nmf.second )
          val = nmf.second / remaining_frac;
        
        usedmassfrac += nmf.second;
        inputPrams.Add( name, val, max(0.1*val,0.01), 0, 1.0 );
        ++num_fit_params;
      }//for( size_t j = 0; i < nmassfrac; ++j )
    }//if( fit for mass fractions to )
  }//for( size_t shielding_index = 0; shielding_index < shieldings.size(); ++shielding_index )
  
  
  if( num_fit_params != inputPrams.VariableParameters() )
    throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...): "
                        "there is a serious logic error in this function, "
                        "please let wcjohns@sandia.gov know about this." );
  
  
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
    m_materials{},
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

      
      
  // Convert from input shieldings to what we use internally
  for( const ShieldingSourceFitCalc::ShieldingInfo &info : shieldings )
  {
    ShieldingSourceChi2Fcn::ShieldLayerInfo materialAndSrc;
    materialAndSrc.material = info.m_material;
    
    for( const pair<const SandiaDecay::Nuclide *,double> &nf : info.m_nuclideFractions )
      materialAndSrc.self_atten_sources.emplace_back( nf.first, nf.second );
    
    //if( info.m_fitMassFrac )
    //{
    // TODO: / NOTE: In `create(...)`, we will call `setNuclidesToFitMassFractionFor(...)`,
    //               but we could do that here instead, which would make more sense, but would
    //               also need to make sure the logic worked out.
    //}
    
    for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
    {
      const SandiaDecay::Nuclide * const nuc = trace.m_nuclide;
      GammaInteractionCalc::TraceActivityType type = trace.m_type;
      
      const double relax_len = (type == TraceActivityType::ExponentialDistribution)
                                ? trace.m_relaxationDistance : -1.0;
      
#if( defined(__GNUC__) && __GNUC__ < 5 )
      materialAndSrc.trace_sources.push_back( tuple<const SandiaDecay::Nuclide *,TraceActivityType,double>{nuc, type, relax_len} );
#else
      materialAndSrc.trace_sources.push_back( {nuc, type, relax_len} );
#endif
    }//for( loop over trace source nuclides )
    
    m_materials.push_back( materialAndSrc );
  }//for( const ShieldingSourceFitCalc::ShieldingInfo &info : shieldings )
      
  assert( m_materials.size() == shieldings.size() );

    
  // Go through and sanity-check traceSources
  for( const ShieldingSourceChi2Fcn::ShieldLayerInfo &info : m_materials )
  {
    const shared_ptr<const Material> &material = info.material;
    
    if( !material )
    {
      if( !info.self_atten_sources.empty() || !info.trace_sources.empty() )
        throw runtime_error( "ShieldingSourceChi2Fcn: Self attenuating and trace sources must"
                             " be empty for a generic shielding" );
      continue;
    }//if( !material )
    
    std::map<short,double> self_atten_mass_fracs;
    for( const pair<const SandiaDecay::Nuclide *,double> &nuc_frac : info.self_atten_sources )
    {
      const SandiaDecay::Nuclide * const nuc = nuc_frac.first;
      const double mass_frac = nuc_frac.second;
      
      if( !nucs.count(nuc) )
        throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                            + nuc->symbol + " doesnt have a peak with that as an assigned nuclide" );
      
      if( !self_atten_mass_fracs.count(nuc->atomicNumber) )
        self_atten_mass_fracs[nuc->atomicNumber] = 0.0;
      self_atten_mass_fracs[nuc->atomicNumber] += mass_frac;
      
      int hasElement = 0;
      for( const auto &p : material->elements )
        hasElement += (p.first->atomicNumber == nuc->atomicNumber);
      
      for( const auto &p : material->nuclides )
        hasElement += (p.first == nuc);
      
      if( !hasElement )
        throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                            + nuc->symbol + " not in shielding " + material->name );
    }//for( const SandiaDecay::Nuclide *nuc : info.self_atten_sources )
  
/*
    // Check source nuclides dont add up to a larger fraction of the material, then
    //  the individual element of the nuclide should.
    if( !self_atten_mass_fracs.empty() )
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      for( const auto &i : self_atten_mass_fracs )
      {
        const SandiaDecay::Element * const el = db->element( i.first );
        assert( el );
        if( el )
        {
          const double elFrac = material->massFractionOfElementInMaterial(el);
          assert( i.second <= (elFrac + 1.0E-6) );
          if( i.second > (elFrac + 1.0E-6) )
            throw runtime_error( el->symbol + " has source-nuclide fraction (" + to_string(elFrac)
                                + ") greater than elements fraction (" + to_string(i.second)
                                + ") in material " + material->name  );
        }//if( el )
      }//for( const auto &i : self_atten_mass_fracs )
    }//if( !self_atten_mass_fracs.empty() )
*/
    
    for( const tuple<const SandiaDecay::Nuclide *,TraceActivityType,double> &trace : info.trace_sources )
    {
      const SandiaDecay::Nuclide *nuc = std::get<0>(trace);
      if( !nuc )
        throw runtime_error( "ShieldingSourceChi2Fcn: null trace source" );
      
      if( !nucs.count(nuc) )
        throw runtime_error( "ShieldingSourceChi2Fcn: self attenuating nuclide "
                  + nuc->symbol + " doesnt have a peak with that as an assigned nuclide" );
    }//for( loop over trace sources )
    
    // Not checking if trace and self-atten sources share the same nuclide, because although the
    //  gui doesnt currently allow this, I *think* the computation will work out fine, maybe.
    
  }//for( const TraceSourceInfo &info : traceSources )
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
    throw runtime_error( "There are not the same number of source definitions as nuclides in peaks being used." );
  
  
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
  m_materials = rhs.m_materials;
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
  if( material_index >= m_materials.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::isVariableMassFraction: invalid material index" );
  
  const ShieldLayerInfo &shield = m_materials[material_index];
  const vector<const SandiaDecay::Nuclide *> &nucs = shield.nucs_to_fit_mass_fraction_for;
  return (std::find( begin(nucs), end(nucs), nuc) != end(nucs));
}//isVariableMassFraction(...)
  
  
bool ShieldingSourceChi2Fcn::hasVariableMassFraction( const size_t material_index ) const
{
  if( material_index >= m_materials.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::hasVariableMassFraction: invalid material index" );
  
  const vector<const SandiaDecay::Nuclide *> &nucs = m_materials[material_index].nucs_to_fit_mass_fraction_for;
  return !nucs.empty();
}//bool hasVariableMassFraction( const Material *material ) const

  
std::shared_ptr<Material> ShieldingSourceChi2Fcn::variedMassFracMaterial(
                                                              const size_t material_index,
                                                              const std::vector<double> &x ) const
{
  if( material_index >= m_materials.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::variedMassFracMaterial: invalid material index" );
  
  const ShieldLayerInfo &shield = m_materials[material_index];
  
  if( !shield.material || shield.self_atten_sources.empty() )
      throw runtime_error( "variedMassFracMaterial(): "
                           + (shield.material ? shield.material->name : string("generic"))
                           + " does not have self-attenuating nuclides" );
    
  std::shared_ptr<Material> answer = std::make_shared<Material>( *shield.material );
  
  //if( !SpecUtils::icontains( answer->name, " - varied") )
  //  answer->name += " - varied";
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  vector<pair<const SandiaDecay::Nuclide *,float>> &nuclides = answer->nuclides;
  vector<pair<const SandiaDecay::Element *,float>> &elements = answer->elements;
  
  // First get mass-fractions of initial material
  map<const SandiaDecay::Element *,double> new_elements;
  for( const auto &el_frac : elements )
  {
    assert( el_frac.first );
    if( !el_frac.first )
      continue;
    
    if( !new_elements.count(el_frac.first) )
      new_elements[el_frac.first] = 0.0;
    new_elements[el_frac.first] += el_frac.second;
  }//for( const auto &el_frac : elements )
  
  for( const auto &nuc_frac : nuclides )
  {
    assert( nuc_frac.first );
    if( !nuc_frac.first )
      continue;
    
    const SandiaDecay::Element * const el = db->element(nuc_frac.first->atomicNumber);
    assert( el );
    if( !el )
      continue;
    
    if( !new_elements.count(el) )
      new_elements[el] = 0.0;
    new_elements[el] += nuc_frac.second;
  }//for( const auto &nuc_frac : nuclides )
  
  // Get mass fractions of elements being used as self-attenuating sources.
  double self_atten_sum = 0.0;
  map<const SandiaDecay::Element *,double> self_atten_elements;
  for( const auto &nuc_frac : shield.self_atten_sources )
  {
    assert( nuc_frac.first );
    if( !nuc_frac.first )
      continue;
    
    const SandiaDecay::Element * const el = db->element(nuc_frac.first->atomicNumber);
    assert( el );
    if( !el )
      continue;
    
    if( !new_elements.count(el) )
      throw runtime_error( "variedMassFracMaterial(): Nuclide " + nuc_frac.first->symbol
                          + " does not have a Element in the material" );
    
    if( !self_atten_elements.count(el) )
      self_atten_elements[el] = 0.0;
    self_atten_elements[el] += nuc_frac.second;
    self_atten_sum += nuc_frac.second;
  }//for( const auto &nuc_frac : shield.self_atten_sources )
  
  assert( self_atten_sum < (1.0 + 1.0E-6) );
  if( self_atten_sum >= 1.0 )
  {
    for( auto &i : self_atten_elements )
      i.second /= self_atten_sum;
    self_atten_sum = 1.0;
    
    // Remove all elements that aren't self-attenuating
    vector<const SandiaDecay::Element *> to_remove;
    for( const auto &i : new_elements )
    {
      if( !self_atten_elements.count(i.first) )
        to_remove.push_back( i.first );
    }//for( const auto &i : new_elements )
    
    for( const auto el : to_remove )
      new_elements.erase( el );
  }//if( self_atten_sum > 1.0 )
  
  double initial_non_self_atten_el_sum = 0.0;
  for( const auto &i : new_elements )
  {
    if( !self_atten_elements.count(i.first) )
      initial_non_self_atten_el_sum += i.second;
  }//for( const auto &i : new_elements )
  
  // Scale non-self-attenuating elements, so the sum of everything is 1.0
  const double final_non_self_atten_frac = std::max( 1.0 - self_atten_sum, 0.0 );
  const double non_self_atten_mult = final_non_self_atten_frac / initial_non_self_atten_el_sum;
  for( auto &i : new_elements )
  {
    if( !self_atten_elements.count(i.first) )
      i.second /= (self_atten_sum);
  }//for( const auto &i : new_elements )
  
  // Now alter mass-fractions for the varied nuclides
  double prefrac = 0.0, postfrac = 0.0;
  for( const SandiaDecay::Nuclide * const nuc : shield.nucs_to_fit_mass_fraction_for )
  {
    assert( nuc );
    if( !nuc )
      continue;
    
    double initial_frac = -1.0;
    for( const auto &i : shield.self_atten_sources )
    {
      if( i.first == nuc )
      {
        initial_frac = i.second;
        break;
      }
    }//for( const auto &i : shield.self_atten_sources )
    assert( initial_frac >= 0.0 );
    if( initial_frac < 0.0 )
      throw runtime_error( "variedMassFracMaterial(): Nuc fitting mass-fraction for doesnt appear in self-atten. sources?" );
    
    const double fit_frac = massFraction( material_index, nuc, x );
    
    prefrac += initial_frac;
    postfrac += fit_frac;
    
    const SandiaDecay::Element * const el = db->element(nuc->atomicNumber);
    assert( el );
    if( !el )
      continue;
    
    new_elements[el] += (fit_frac - initial_frac);
    assert( new_elements[el] >= 0.0 );
  }//for( const SandiaDecay::Nuclide * const nuc : shield.nucs_to_fit_mass_fraction_for )
  
  
  if( IsNan(postfrac) || IsInf(postfrac)
      || ((fabs(prefrac-postfrac)/max(prefrac,postfrac)) > 0.001) )
  {
    cerr << "ShieldingSourceChi2Fcn::variedMassFracMaterial: prefrac="
          << prefrac << ", postfrac=" << postfrac << endl;

    throw runtime_error( "ShieldingSourceChi2Fcn::variedMassFracMaterial(...)"
                         " prefrac did not match postfrac" );
  }//if( invalid results )
  
  
  nuclides.clear();
  elements.clear();
  
  for( auto &i : new_elements )
    elements.emplace_back( i.first, static_cast<float>(i.second) );
      
  return answer;
}//variedMassFracMaterial(...)

  
void ShieldingSourceChi2Fcn::setNuclidesToFitMassFractionFor( const size_t material_index,
                          const vector<const SandiaDecay::Nuclide *> &nuclides )
{
  if( material_index >= m_materials.size() )
    throw std::logic_error( "ShieldingSourceChi2Fcn::setNuclidesToFitMassFractionFor: invalid material index" );
  
  ShieldLayerInfo &shielding = m_materials[material_index];
  const shared_ptr<const Material> &material = shielding.material;
  
  if( !material )
    throw std::logic_error( "ShieldingSourceChi2Fcn::setNuclidesToFitMassFractionFor: material index points to generic shielding" );
  
  // TODO: 20230708 for the moment, our logic is that we fit the fractions for all self-attenuating sources... we could change this - I think things are all in place, but it is untested.
  assert( nuclides.size() == shielding.self_atten_sources.size() );
  
  auto &curr_var_srcs = shielding.nucs_to_fit_mass_fraction_for;
  
  for( const SandiaDecay::Nuclide *nuc : nuclides )
  {
    if( std::count( nuclides.begin(), nuclides.end(), nuc ) != 1 )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): input must"
                           " only be unique nuclides" );
    int innuc = 0;
    for( const Material::NuclideFractionPair &nfp : material->nuclides )
      innuc += int(nfp.first == nuc);
    
    for( const auto &efp : material->elements )
      innuc += int(efp.first->atomicNumber == nuc->atomicNumber);
    
    if( !innuc )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): passed in"
                           " nuclide " + nuc->symbol + " not in material "
                           + material->name );
    
    if( std::count( m_nuclides.begin(), m_nuclides.end(), nuc ) != 1 )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): input must "
                           "be in source m_nuclides" );
    
    const auto &self_atten_srcs = shielding.self_atten_sources;
    
    bool srcIsSelfAtten = false;
    for( size_t i = 0; !srcIsSelfAtten && (i < self_atten_srcs.size()); ++i )
      srcIsSelfAtten = (self_atten_srcs[i].first == nuc);
    
    if( !srcIsSelfAtten )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): "
                           + nuc->symbol
                           + " is not in a self-attenuating shielding" );
    
    if( std::find(begin(curr_var_srcs), end(curr_var_srcs), nuc) != end(curr_var_srcs) )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): "
                           + nuc->symbol
                           + " is already a variable mass-fraction self-attenuating nuclide" );
  }//for( const SandiaDecay::Nuclide *nuc : nuclides )
  
  curr_var_srcs.insert( end(curr_var_srcs), begin(nuclides), end(nuclides) );
  
  std::sort( begin(curr_var_srcs), end(curr_var_srcs),
    []( const SandiaDecay::Nuclide *lhs, const SandiaDecay::Nuclide *rhs ) -> bool {
      if( !lhs ) return false;
      if( !rhs ) return true;
      return (lhs->symbol < rhs->symbol);
    } );
}//setNuclidesToFitMassFractionFor(...)
  
vector<const Material *> ShieldingSourceChi2Fcn::materialsFittingMassFracsFor() const
{
  vector<const Material *> answer;
  
  for( const ShieldLayerInfo &ms : m_materials )
  {
    if( !ms.nucs_to_fit_mass_fraction_for.empty() )
      answer.push_back( ms.material.get() );
  }//for( ShieldLayerInfo &ms : m_materials )
  
  return answer;
}//vector<const Material *material> materialsFittingMassFracsFor() const
  

vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::selfAttenuatingNuclides( const size_t material_index ) const
{
  if( material_index >= m_materials.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );

  const ShieldingSourceChi2Fcn::ShieldLayerInfo &ms = m_materials[material_index];
  
  if( !ms.material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::selfAttenuatingNuclides():"
                        " material index points to generic shielding" );
  
  vector<const SandiaDecay::Nuclide *> answer;
  for( const auto &i : ms.self_atten_sources )
    answer.push_back( i.first );
  
  return answer;
}//vector<const SandiaDecay::Nuclide *> selfAttenuatingNuclides( const size_t material_index ) const
  
  
vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::traceNuclidesForMaterial( const size_t material_index ) const
{
  if( material_index >= m_materials.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );

  const ShieldingSourceChi2Fcn::ShieldLayerInfo &ms = m_materials[material_index];
  
  if( !ms.material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::traceNuclidesForMaterial():"
                        " material index points to generic shielding" );
  
  vector<const SandiaDecay::Nuclide *> answer;
  for( const auto &i : ms.trace_sources )
    answer.push_back( get<0>(i) );
  
  return answer;
}//std::vector<const SandiaDecay::Nuclide *> ShieldingSourceChi2Fcn::traceNuclidesForMaterial( const size_t material_index ) const
  
const std::vector<const SandiaDecay::Nuclide *> &
        ShieldingSourceChi2Fcn::nuclideFittingMassFracFor( const size_t material_index ) const
{
  if( material_index >= m_materials.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor: invalid material index" );
  
  const ShieldLayerInfo &ms = m_materials[material_index];
  
  if( !ms.material ) // JIC
    throw logic_error( "ShieldingSourceChi2Fcn::nuclideFittingMassFracFor():"
                        " material index points to generic shielding" );

  return ms.nucs_to_fit_mass_fraction_for;
}//nuclideFittingMassFracFor(...)

  
double ShieldingSourceChi2Fcn::massFraction( const size_t material_index,
                      const SandiaDecay::Nuclide *nuc,
                      const std::vector<double> &pars ) const
{
  double massfrac = 0.0, uncert = 0.0;
  massFraction( massfrac, uncert, material_index, nuc, pars, vector<double>() );
  return massfrac;
}
  
  
void ShieldingSourceChi2Fcn::massFraction( double &massFrac, double &uncert,
                                  const size_t material_index, const SandiaDecay::Nuclide *nuc,
                                  const vector<double> &pars, const vector<double> &errors ) const
{
  massFrac = uncert = 0.0;
  
  if( material_index >= m_materials.size() )
    throw logic_error( "ShieldingSourceChi2Fcn::massFraction: invalid material index" );
  
  const ShieldLayerInfo &shielding = m_materials[material_index];
  
  if( !shielding.material || !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): invalid input" );
  
  //if( shielding.nucs_to_fit_mass_fraction_for.empty() )
  //  throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): "
  //                       + shielding.material->name + " is not a material with a variable"
  //                       " mass fraction" );
  
  const vector<const SandiaDecay::Nuclide *> &nucs = shielding.nucs_to_fit_mass_fraction_for;
  
  //nucs is actually sorted by symbol name, could do better than linear search
  const auto pos = std::find( begin(nucs), end(nucs), nuc );
  if( pos == end(nucs) )
  {
    // Not a variable mass-fraction
    for( const auto &i : shielding.self_atten_sources )
    {
      if( i.first == nuc )
      {
        massFrac = i.second;
        uncert = 0.0;
        return;
      }//if( i.first == nuc )
    }//for( const auto &i : shielding.self_atten_sources )
    
    throw runtime_error( "ShieldingSourceChi2Fcn::massFraction(): "
                         + nuc->symbol + " was not a self-attenuating nuclide in "
                         + shielding.material->name );
  }//if( pos == end(nucs) )
    
  
  double totalfrac = 0.0;
  for( const pair<const SandiaDecay::Nuclide *,double> &i : shielding.self_atten_sources )
  {
    const SandiaDecay::Nuclide * const this_nuc = i.first;
    const auto this_pos = std::find( begin(nucs), end(nucs), this_nuc );
    if( this_pos != end(nucs) )
      totalfrac += i.second;
  }

  size_t matmassfracstart = 0;
  for( size_t index = 0; index < material_index; ++index )
  {
    const ShieldLayerInfo &info = m_materials[index];
    if( !info.nucs_to_fit_mass_fraction_for.empty() )
      matmassfracstart += (info.nucs_to_fit_mass_fraction_for.size() - 1);
  }
  
  matmassfracstart += 2 * m_nuclides.size();
  matmassfracstart += 3 * m_materials.size();
  
  const size_t numfraccoefs = size_t(nucs.size() - 1);
  const size_t massfracnum = (pos - nucs.begin());
  const size_t thismatmassfrac = matmassfracstart + massfracnum;
  
  double frac = 1.0;
  for( size_t index = matmassfracstart; index < thismatmassfrac; ++index )
    frac *= (1.0 - pars.at(index));
  double prefrac = frac;
  
  if( massfracnum != numfraccoefs )
    frac *= pars.at(thismatmassfrac);
  
  massFrac = totalfrac * frac;
  
  if( IsNan(massFrac) || IsInf(massFrac) )
  {
    cerr << "Got invalid mass frac:" << endl;
    cerr << "\ttotalfrac=" << totalfrac << endl;
    cerr << "\tprefrac=" << prefrac << endl;
    if( massfracnum != numfraccoefs )
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
    if( massfracnum != numfraccoefs )
      fracuncert = errors.at(thismatmassfrac) / pars.at(thismatmassfrac);
    else
      fracuncert = errors.at(thismatmassfrac-1) / pars.at(thismatmassfrac-1);
    uncert = fracuncert * massFrac;
  }//if( !errors.empty() )
}//double massFraction(...) const

  
double ShieldingSourceChi2Fcn::massFractionUncert( const size_t material_index,
                            const SandiaDecay::Nuclide *nuc,
                            const std::vector<double> &pars,
                            const std::vector<double> &error ) const
{
  double massfrac = 0.0, uncert = 0.0;
  massFraction( massfrac, uncert, material_index, nuc, pars, error );
  return uncert;
}//massFractionUncert(...)
  

size_t ShieldingSourceChi2Fcn::numExpectedFitParameters() const
{
  size_t npar = 2 * m_nuclides.size();
  
  npar += 3 * m_materials.size();
  
  for( const ShieldLayerInfo &info : m_materials )
    npar += (info.nucs_to_fit_mass_fraction_for.empty() ? size_t(0) : size_t(info.nucs_to_fit_mass_fraction_for.size() - 1) );
  
  
  return npar;
}//int numExpectedFitParameters() const


bool ShieldingSourceChi2Fcn::isVolumetricSource( const SandiaDecay::Nuclide *nuc ) const
{
  return isSelfAttenSource(nuc) || isTraceSource(nuc);
}


bool ShieldingSourceChi2Fcn::isSelfAttenSource( const SandiaDecay::Nuclide *nuclide ) const
{
  //TODO: this could probably be made a little more efficient
  if( !nuclide )
    return false;

  for( const ShieldLayerInfo &shield : m_materials )
  {
    for( const auto &nuc_frac : shield.self_atten_sources )
    {
      if( nuc_frac.first == nuclide )
        return true;
    }
  }//for( const ShieldLayerInfo &material : m_materials )

  return false;
}//bool isSelfAttenSource(...) const;


bool ShieldingSourceChi2Fcn::isTraceSource( const SandiaDecay::Nuclide *nuclide ) const
{
  //TODO: this could probably be made a little more efficient
  if( !nuclide )
    return false;
  
  for( const ShieldLayerInfo &shield : m_materials )
  {
    for( const tuple<const SandiaDecay::Nuclide *,TraceActivityType,double> &trace : shield.trace_sources )
    {
      if( get<0>(trace) == nuclide )
        return true;
    }
  }//for( const ShieldLayerInfo &material : m_materials )
  
  return false;
}//bool isTraceSource(...) const;


TraceActivityType ShieldingSourceChi2Fcn::traceSourceActivityType(
                                                          const SandiaDecay::Nuclide *nuc ) const
{
  //TODO: this could probably be made a little more efficient
  if( !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: null nuclide" );
  
  for( const ShieldLayerInfo &shield : m_materials )
  {
    for( const tuple<const SandiaDecay::Nuclide *,TraceActivityType,double> &trace : shield.trace_sources )
    {
      if( get<0>(trace) == nuc )
        return get<1>(trace);
    }
  }//for( const ShieldLayerInfo &material : m_materials )
  
  throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: " + nuc->symbol
                       + " not a trace source" );
  
  return TraceActivityType::NumTraceActivityType;
}//traceSourceActivityType(...)


double ShieldingSourceChi2Fcn::relaxationLength( const SandiaDecay::Nuclide *nuc ) const
{
  //TODO: this could probably be made a little more efficient
  if( !nuc )
    throw runtime_error( "ShieldingSourceChi2Fcn::traceSourceActivityType: null nuclide" );
  
  for( const ShieldLayerInfo &shield : m_materials )
  {
    for( const tuple<const SandiaDecay::Nuclide *,TraceActivityType,double> &trace : shield.trace_sources )
    {
      if( get<0>(trace) == nuc )
      {
        if( get<1>(trace) != TraceActivityType::ExponentialDistribution )
          throw runtime_error( "ShieldingSourceChi2Fcn::relaxationLength: " + nuc->symbol
                               + " is not an exponential distribution trace source." );
        return get<2>(trace);
      }
    }
  }//for( const ShieldLayerInfo &material : m_materials )
  
  throw runtime_error( "ShieldingSourceChi2Fcn::relaxationLength: " + nuc->symbol
                      + " not a trace source" );
  
  return -1.0;
}//double relaxationLength( const SandiaDecay::Nuclide *nuc ) const;


double ShieldingSourceChi2Fcn::volumeOfMaterial( const size_t matn, const vector<double> &params ) const
{
  if( matn >= m_materials.size() )
    throw runtime_error( "volumeOfMaterial: invalid material index" );
  
  if( !m_materials[matn].material )
    throw runtime_error( "volumeOfMaterial: cant be called for generic material" );
  
  double inner_dim_1 = 0.0, inner_dim_2 = 0.0, inner_dim_3 = 0.0;
  double outer_dim_1 = 0.0, outer_dim_2 = 0.0, outer_dim_3 = 0.0;
  
  for( int index = 0; index <= matn; ++index )
  {
    if( !m_materials[index].material )  //if a generic shielding
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
  if( (matn < 0) || (matn >= m_materials.size()) )
    throw runtime_error( "volumeUncertaintyOfMaterial: invalid material index" );
  
  if( !m_materials[matn].material )
    throw runtime_error( "volumeUncertaintyOfMaterial: cant be called for generic material" );
  
  // TODO: clean this up to be a little neater/shorter, after testing
  double dim_1 = 0.0, dim_2 = 0.0, dim_3 = 0.0;
  double dim_1_uncert_sq = 0.0, dim_2_uncert_sq = 0.0, dim_3_uncert_sq = 0.0;
  
  for( int mat_index = 0; mat_index <= matn; ++mat_index )
  {
    if( !m_materials[mat_index].material )  //if a generic shielding
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

  const size_t num_mataterials = m_materials.size();
  for( size_t material_index = 0; material_index < num_mataterials; ++material_index )
  {
    const ShieldLayerInfo &shield = m_materials[material_index];
    
    const shared_ptr<const Material> &mat = shield.material;
    if( !mat )  //if a generic shielding
      continue;
    
    // Determine if nuclide is a self-attenuating source of this shielding
    const auto &self_atten_srcs = shield.self_atten_sources;
    vector<pair<const SandiaDecay::Nuclide *,double>>::const_iterator self_att_pos;
    for( self_att_pos = begin(self_atten_srcs); self_att_pos != end(self_atten_srcs); ++self_att_pos )
    {
      if( self_att_pos->first == nuclide )
        break;
    }
    
    // Check that nuclide is a self-atten source of this shielding, if not lets keep looking.
    if(self_att_pos == end(self_atten_srcs) )
      continue;
    
    foundSrc = true;
    
    double massFrac = 0.0;
    
    if( hasVariableMassFraction(material_index) )
    {
      massFrac = ShieldingSourceChi2Fcn::massFraction( material_index, nuclide, params );
    }else
    {
      massFrac = self_att_pos->second;
    }//if( hasVariableMassFraction(mat) ) / else
    
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
  
  const int nmat = static_cast<int>( numMaterials() );
  for( int matn = 0; matn < nmat; ++matn )
  {
    const ShieldLayerInfo &shield = m_materials[matn];
    
    const shared_ptr<const Material> &mat = shield.material;
    if( !mat )  //if a generic shielding
      continue;
    
    const double vol = volumeOfMaterial( matn, params );
    
    // Determine if nuclide is a trace source of this shielding
    TraceActivityType type = TraceActivityType::NumTraceActivityType;
    const auto &trace_srcs = shield.trace_sources;
    for( const auto &trace : trace_srcs )
    {
      if( get<0>(trace) == nuclide )
      {
        type = get<1>(trace);
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
  
  const int num_materials = static_cast<int>( numMaterials() );
  for( int material_index = 0; material_index < num_materials; ++material_index )
  {
    if( isGenericMaterial(material_index) )
      continue;
    
    const ShieldLayerInfo &shield = m_materials[material_index];
    const shared_ptr<const Material> &mat = shield.material;
    
    if( shield.self_atten_sources.empty() && shield.trace_sources.empty() )
      continue;
    
    //const auto &self_attens = shield.self_atten_sources;
    //const auto &traces = shield.trace_sources;
    //const bool is_self_atten = (std::find(begin(self_attens), end(self_attens), nuclide) != end(self_attens));
    
    const double vol = volumeOfMaterial(material_index, params);
    const double volUncert = volumeUncertaintyOfMaterial(material_index, params, errors);
  
    // Add in uncertainty contributions if this is a self attenuating source
    for( const auto &src_frac : shield.self_atten_sources )
    {
      const SandiaDecay::Nuclide *src = src_frac.first;
      
      if( src != nuclide )
        continue;
      
      double massFrac = 0.0, massFracUncert = 0.0;
      
      if( hasVariableMassFraction(material_index) )
      {
        massFraction( massFrac, massFracUncert, material_index, src, params, errors );
      }else
      {
        massFrac = src_frac.second;
      }//if( hasVariableMassFraction(material_index) ) / else
      
      
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
  
    // Add in uncertainty contributions if this is a trace source
    for( const tuple<const SandiaDecay::Nuclide *,TraceActivityType,double> &trace : shield.trace_sources )
    {
      if( get<0>(trace) != nuclide )
        continue;
      
      const TraceActivityType trace_type = get<1>(trace);
      
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
              const double r =  sphericalThickness(material_index, params);
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
    throw runtime_error( "error calculating activity uncertainty for self-attenuating source;"
                         " squared value calculated is " + std::to_string(activityUncertSquared) );
  
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
  return m_materials.size();
}//int numMaterials() const


const Material *ShieldingSourceChi2Fcn::material( const size_t materialNum ) const
{
  return m_materials.at(materialNum).material.get();
}

bool ShieldingSourceChi2Fcn::isSpecificMaterial( const size_t materialNum ) const
{
  return (m_materials.at(materialNum).material != nullptr );
}//bool isSpecificMaterial( int materialNum ) const


bool ShieldingSourceChi2Fcn::isGenericMaterial( const size_t materialNum ) const
{
  return (m_materials.at(materialNum).material == nullptr );
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
    
    const vector< tuple<double,double,double,Wt::WColor,double> > chi2s
                                           = energy_chi_contributions( x, m_mixtureCache, nullptr );
    double chi2 = 0.0;
    
    const size_t npoints = chi2s.size();
    for( size_t i = 0; i < npoints; ++i )
      chi2 += pow( std::get<1>(chi2s[i]), 2.0 );
    
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
                                                           vector<string> *info )
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
    msg << " at age " << PhysicalUnits::printToBestTimeUnits(age) << ":";
    info->push_back( msg.str() );
  }//if( info )

  if( mixture.numInitialNuclides() != 1 )
    throw runtime_error( "ShieldingSourceChi2Fcn::cluster_peak_activities():"
                         " passed in mixture must have exactly one parent nuclide" );
  const SandiaDecay::Nuclide *nuclide = mixture.initialNuclide(0);
  assert( nuclide );
  if( !nuclide )
    throw std::logic_error( "nullptr nuc" );
  
  vector<SandiaDecay::EnergyRatePair> gammas
                = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );
  
  assert( !accountForDecayDuringMeas || (measDuration > 0.0) );  //This could happen, I guess if foreground real-time is zero, but it shouldnt, so I'll leave this assert in to check for things.
  
  if( accountForDecayDuringMeas && (measDuration > 0.0)
     && !IsNan(measDuration) && !IsInf(measDuration) )
  {
    // This section of code takes up a good amount of CPU time, and is really begging for optimizations
    //
    // TODO: if nuclide decays to stable children, then just use standard formula to correct for decay (see below developer check for this)
    //
    // We want to get the activity at the beggining of the measurement, but we want
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
    const size_t characteristic_time_slices = 50; // Maybe at least ~5 sig figs
    
    // TODO: the parent half-life may not be the relevant one; should check down the chain
    const double halflife = std::max( nuclide->halfLife, 0.01*PhysicalUnits::second );
    const double characteristicTime = std::min( measDuration, halflife );
    const double dt = characteristicTime / characteristic_time_slices;
    
    const int num_time_slices = std::max( 50, static_cast<int>( std::ceil( measDuration / dt ) ) );
    const int nthread = SpecUtilsAsync::num_logical_cpu_cores();
    
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
        corrected_gammas[i].numPerSecond += energy_rate_pairs[threadnum][i].numPerSecond;
    }//for( size_t i = 0; i < corrected_gammas.size(); ++i )
    
    //cout << "For " << nuclide->symbol << " the time corrections are:" << endl;
    //for( size_t i = 0; i < corrected_gammas.size(); ++i )
    //  cout << std::setw(15) << corrected_gammas[i].energy << ": "
    //       << corrected_gammas[i].numPerSecond/gammas[i].numPerSecond
    //       << " (" << corrected_gammas[i].numPerSecond << " vs orig "
    //       << gammas[i].numPerSecond << ")" << endl;
    
#if( PERFORM_DEVELOPER_CHECKS )
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
    
    if( info )
    {
      // TODO: we could add something about effect of accounting for nuclide aging
    }//if( info )
    
    gammas = std::move( corrected_gammas );
  }//if( accountForDecayDuringMeas )
  
  
  //The problem we have is that 'gammas' have the activity of the original
  //  parent ('nuclide') decreased by agining by 'age', however we want the
  //  parent to have 'sm_activityUnits' activity at 'age', so we we'll add a
  //  correction factor.
  const vector<SandiaDecay::NuclideActivityPair> aged_activities
                                                      = mixture.activity( age );

  
  double age_sf = 1.0;
  for( const SandiaDecay::NuclideActivityPair &red : aged_activities )
  {
    if( red.nuclide == nuclide )
    {
      age_sf = 1.0*sm_activityUnits / red.activity;
      break;
    }
  }//for( size_t i = 0; i < aged_activities.size(); ++i )
  
  
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
          << age_sf * aep.numPerSecond/sm_activityUnits << "";
      info->push_back( msg.str() );
    }//if( info )
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


vector< tuple<double,double,double,Wt::WColor,double> > ShieldingSourceChi2Fcn::expected_observed_chis(
                                           const std::vector<PeakDef> &peaks,
                                           const std::vector<PeakDef> &backPeaks,
                                           const std::map<double,double> &energy_count_map,
                                           vector<string> *info )
{
  typedef map<double,double> EnergyCountMap;

  
  if( info )
    info->push_back( "Chi2 Contributions Of Peaks" );
  
  //Go through and match the predicted number of counts to the observed number
  //  of counts and get the chi2.
  //Note that matching betoween expected and observed peaks is done via energy
  //  which make me a bit queezy for some reason
  vector< tuple<double,double,double,Wt::WColor,double> > answer;

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
      observed_counts -= backCounts;
      observed_uncertainty = sqrt( observed_uncertainty*observed_uncertainty + backUncert2 );
    }
    
    const double chi = (observed_counts - expected_counts) / observed_uncertainty;
    const double scale = observed_counts / expected_counts;
    const double scale_uncert = observed_uncertainty / expected_counts;
    answer.emplace_back( make_tuple(energy, chi, scale, peak.lineColor(), scale_uncert) );
    
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
    << "\n\t m_sourceIndex: " << calculator.m_sourceIndex
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
  
  
vector< tuple<double,double,double,Wt::WColor,double> >
       ShieldingSourceChi2Fcn::energy_chi_contributions( const std::vector<double> &x,
                                         ShieldingSourceChi2Fcn::NucMixtureCache &mixturecache,
                                         std::vector<std::string> *info ) const
{
  //XXX - this function compares a lot of doubles, and this always makes me
  //      queezy - this should be checked on!
  typedef map<double,double> EnergyCountMap;

//  cerr << "energy_chi_contributions: vals={ ";
//  for( size_t i = 0; i < x.size(); ++i )
//    cerr << x[i] << ", ";
//  cerr << "}; size=" << x.size() << endl;
  
  if( info )
  {
    info->push_back( "LiveTime="
                    + std::to_string(m_liveTime/PhysicalUnits::second)
                    + " s (" + PhysicalUnits::printToBestTimeUnits(m_liveTime)
                    + ")" );
    info->push_back( "Distance to source center from detector: "
                      + PhysicalUnits::printToBestLengthUnits(m_distance) );
    if( m_detector && m_detector->isValid() )
      info->push_back( "Detector: "
                        + m_detector->name() + " radius "
                        + std::to_string( 0.5*m_detector->detectorDiameter()/PhysicalUnits::cm )
                      + " cm" );
    if( m_options.multiple_nucs_contribute_to_peaks )
      info->push_back( "Allowing multiple nuclides being fit for to potentially contribute to the same photopeak" );
    else
      info->push_back( "Not allowing multiple nuclides being fit for to contribute to the same photopeak" );
    
    if( m_options.account_for_decay_during_meas )
      info->push_back( "Branching ratios are being corrected for nuclide decay during measurment" );
    
    //Should put in information about the shielding here
  }//if( info )
  
  
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
                               info );
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
                               info );
    }//for( const PeakDef &peak : m_peaks )
  }//if( m_options.multiple_nucs_contribute_to_peaks )

  //Propagate the gammas through each material - note we are using the fit peak
  //  mean here, and not the (pre-cluster) photopeak energy
  double shield_outer_rad = 0.0;
  const size_t nMaterials = m_materials.size();
  for( size_t materialN = 0; materialN < nMaterials; ++materialN )
  {
    boost::function<double(float)> att_coef_fcn;
    const ShieldLayerInfo &shielding = m_materials[materialN];
    const shared_ptr<const Material> &material = shielding.material;

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
    
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      const double energy = energy_count.first;
      const double frac = exp( -1.0 * att_coef_fcn( energy ) );
      
      energy_count.second *= frac;
    }
  }//for( int materialN = 0; materialN < nMaterials; ++materialN )

  
  if( m_options.attenuate_for_air )
  {
    const double air_dist = std::max( 0.0, m_distance - shield_outer_rad );
    
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
      const double coef = transmission_length_coefficient_air( energy_count.first );
      energy_count.second *= exp( -1.0 * coef * air_dist );
    }
  }//if( m_options.attenuate_for_air )
  
  
  //Fold in the detector response
  if( m_detector && m_detector->isValid() )
  {
    if( info )
      info->push_back( "Detector Efficiency Effects" );
    
    for( EnergyCountMap::value_type &energy_count : energy_count_map )
    {
//      cerr << "Absolute efficiency at " << energy_count.first << " keV is "
//           << m_detector->intrinsicEfficiency( energy_count.first ) << " and the "
//           << " total efficiency is " << m_detector->efficiency( energy_count.first, m_distance ) << endl;
      
      const double eff = m_detector->efficiency( energy_count.first, m_distance );
      
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
  }//if( m_detector && m_detector->isValid() ) / else


  //This is where contributions from self-attenuating and traces source are calculated
  
  // We'll make a copy of materials since we may mass-fraction vary the isotopics
  vector<ShieldLayerInfo> materials = m_materials;
  
  if( info )
  {
    for( const auto &shield : materials )
    {
      if( shield.material && (!shield.self_atten_sources.empty() || !shield.trace_sources.empty()) )
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
  
  vector<DistributedSrcCalc> calculators;
  bool has_trace = false, has_self_atten = false;
  
  for( size_t material_index = 0; material_index < nMaterials; ++material_index )
  {
    const ShieldLayerInfo &shield = materials[material_index];
    shared_ptr<const Material> material = shield.material;
    
    if( !material )
      continue;
    
    const vector<pair<const SandiaDecay::Nuclide *,double>> &self_atten_srcs = shield.self_atten_sources;
    const vector<tuple<const SandiaDecay::Nuclide *,TraceActivityType,double>> &trace_srcs = shield.trace_sources;
    
    if( !self_atten_srcs.empty() )
    {
      // TODO: We could avoid creating a new material if we are only fitting for mass-fractions of
      //       nuclides within one element, and their sum is less than the elements mass-fraction.
      material = variedMassFracMaterial( material_index, x );
      local_materials.push_back( material );
    }//if( !self_atten_srcs.empty() )
    
    
    vector<const SandiaDecay::Nuclide *> combined_srcs;
    for( const auto &p : self_atten_srcs )
      combined_srcs.push_back( p.first );
    for( const auto &p : trace_srcs )
      combined_srcs.push_back( std::get<0>(p) );
    
    const auto &srcs = combined_srcs;
    
    if( srcs.empty() )
      continue;
    
    DistributedSrcCalc baseCalculator;
    baseCalculator.m_geometry = m_geometry;
    
    if( m_detector )
      baseCalculator.m_detectorRadius = 0.5 * m_detector->detectorDiameter();
    else
      baseCalculator.m_detectorRadius = 0.5 * PhysicalUnits::cm;

    baseCalculator.m_observationDist = m_distance;
    baseCalculator.m_attenuateForAir = m_options.attenuate_for_air;
    baseCalculator.m_sourceIndex = material_index;
    
    baseCalculator.m_isInSituExponential = false;
    baseCalculator.m_inSituRelaxationLength = -1.0;
    
    for( size_t src_index = 0; src_index < srcs.size(); ++src_index )
    {
      const SandiaDecay::Nuclide *src = srcs[src_index];
      const bool is_trace = (src_index >= self_atten_srcs.size());
      
#if( PERFORM_DEVELOPER_CHECKS )
      {// begin quick sanity check
        bool trace_check = false;
        for( const auto &p : trace_srcs )
          trace_check = (trace_check || (std::get<0>(p) == src));
        assert( trace_check == is_trace );
      }// end quick sanity check
#endif
      
      double actPerVol = 0.0;
      EnergyCountMap local_energy_count_map;
      
      if( is_trace )
      {
        has_trace = true;
        const size_t trace_index = src_index - self_atten_srcs.size();
        assert( trace_index < trace_srcs.size() );
        
        const double act = activity(src, x);
        
        switch( std::get<1>(trace_srcs[trace_index]) )
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
            const double relaxation_len = std::get<2>(trace_srcs[trace_index]);
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
        }//switch ( trace_srcs[trace_index].second )
      }else //if( is_trace )
      {
        has_self_atten = true;
        
        const double actPerMass = src->activityPerGram() / PhysicalUnits::gram;
        double massFract = 0.0;
        
        
        const ShieldLayerInfo &shielding = m_materials[material_index];
        const vector<const SandiaDecay::Nuclide *> &mass_frac_nucs = shielding.nucs_to_fit_mass_fraction_for;
        
        if( std::count(begin(mass_frac_nucs), end(mass_frac_nucs), src) )
        {
          assert( shielding.material );
          massFract = massFraction( material_index, src, x );
        }else
        {
          for( const auto &src_frac : self_atten_srcs )
          {
            if( src_frac.first == src )
              massFract += src_frac.second;
          }
        }
        
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
                                 info );
      }else
      {
        for( const PeakDef &peak : m_peaks )
        {
          if( peak.parentNuclide()==src
             && (peak.decayParticle() || (peak.sourceGammaType()==PeakDef::AnnihilationGamma)) )
            cluster_peak_activities( local_energy_count_map, energie_widths,
                                     mixturecache[src], actPerVol, thisage,
                                     m_options.photopeak_cluster_sigma,
                                     peak.gammaParticleEnergy(),
                                     m_options.account_for_decay_during_meas, m_realTime,
                                     info );
        }//for( const PeakDef &peak : m_peaks )
      }//if( m_options.multiple_nucs_contribute_to_peaks ) / else

      for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
      {
        DistributedSrcCalc calculator = baseCalculator;

        calculator.m_nuclide = src;
        calculator.m_energy = energy_count.first;
        calculator.m_srcVolumetricActivity = energy_count.second;
        
        if( m_options.attenuate_for_air )
          calculator.m_airTransLenCoef = transmission_length_coefficient_air( energy_count.first );
        else
          calculator.m_airTransLenCoef = 0.0;
        
        if( is_trace )
        {
          const size_t trace_index = src_index - self_atten_srcs.size();
          assert( trace_index < trace_srcs.size() );
          
          if( std::get<1>(trace_srcs[trace_index]) == TraceActivityType::ExponentialDistribution )
          {
            calculator.m_isInSituExponential = true;
            calculator.m_inSituRelaxationLength = std::get<2>(trace_srcs[trace_index]);
            assert( calculator.m_inSituRelaxationLength > 0.0 );
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
            const double transLenCoef = transmition_coefficient_generic( an, ad, calculator.m_energy );
            
            //cout << "Adding generic material (index=" << subMat << ") with AN=" << an
            //     << " and AD=" << ad / (PhysicalUnits::g/PhysicalUnits::cm2) << " g/cm2"
            //     << " atten(" << calculator.m_energy << " keV-->" << transLenCoef << ") = " << exp(-1.0*transLenCoef)
            //     << endl;
            
#if( defined(__GNUC__) && __GNUC__ < 5 )
            calculator.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Generic} );
#else
            calculator.m_dimensionsTransLenAndType.push_back( {outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Generic} );
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
          
          const shared_ptr<const Material> &material = materials[subMat].material;
          
          
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
          
          const double transLenCoef = transmition_length_coefficient( material.get(), calculator.m_energy );

#if( defined(__GNUC__) && __GNUC__ < 5 )
          calculator.m_dimensionsTransLenAndType.push_back( tuple<array<double,3>,double,DistributedSrcCalc::ShellType>{outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Material} );
#else
          calculator.m_dimensionsTransLenAndType.push_back( {outer_dims, transLenCoef, DistributedSrcCalc::ShellType::Material} );
#endif
        }//for( int subMat = 0; subMat < nMaterials; ++subMat )

        if( calculator.m_dimensionsTransLenAndType.empty() )
          throw std::logic_error( "No source/shielding sphere for calculator" );
        
        calculators.push_back( calculator );
      }//for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
    }//for( const SandiaDecay::Nuclide *src : self_atten_srcs )
  }//for( int material_index = 0; material_index < nMaterials; ++material_index )

  if( calculators.size() )
  {
    if( m_options.multithread_self_atten )
    {
      SpecUtilsAsync::ThreadPool pool;
      for( DistributedSrcCalc &calculator : calculators )
        pool.post( boost::bind( &ShieldingSourceChi2Fcn::selfShieldingIntegration, boost::ref(calculator) ) );
      pool.join();
    }else
    {
      for( DistributedSrcCalc &calculator : calculators )
        selfShieldingIntegration(calculator);
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
    
    
    for( const DistributedSrcCalc &calculator : calculators )
    {
      double contrib = calculator.integral * calculator.m_srcVolumetricActivity;

      
      if( m_detector && m_detector->isValid() )
        contrib *= m_detector->intrinsicEfficiency( calculator.m_energy );

      if( energy_count_map.find( calculator.m_energy ) != energy_count_map.end() )
      {
//        cerr << "Adding " << contrib << " to energy " << calculator.m_energy
//             << ", calculator.integral=" << calculator.integral
//             << ", calculator.m_srcVolumetricActivity=" << calculator.m_srcVolumetricActivity
//             << endl;
        energy_count_map[calculator.m_energy] += contrib;
      }else
      {
//        cerr << "Setting " << contrib*m_liveTime << " counts to energy " << calculator.energy
//             << " for thickness=" << calculator.m_dimensionsTransLenAndType[calculator.m_sourceIndex].first[0] / PhysicalUnits::cm
//             << " cm" << endl;
        energy_count_map[calculator.m_energy] = contrib;
      }
      
      if( info )
      {
        const shared_ptr<const Material> &material = materials[calculator.m_sourceIndex].material;
        const int index = static_cast<int>( calculator.m_sourceIndex );
        
        stringstream msg;
        msg << "\tAttributing " << contrib*PhysicalUnits::second << " cps to "
            << calculator.m_energy/PhysicalUnits::keV << " keV photopeak ";
        if( calculator.m_nuclide )
          msg << "(from " << calculator.m_nuclide->symbol << ") ";
        
        msg << "for thicknesses {";
        const array<double,3> &dims = std::get<0>(calculator.m_dimensionsTransLenAndType[index]);
        double dx = dims[0];
        double dy = dims[1];
        double dz = dims[2];
        if( index > 0 )
        {
          const array<double,3> &inner_dims = std::get<0>(calculator.m_dimensionsTransLenAndType[index-1]);
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
    }//for( DistributedSrcCalc &calculator : calculators )
  }//if( calculators.size() )

  
  //Account for live time
  for( EnergyCountMap::value_type &energy_count : energy_count_map )
    energy_count.second *= m_liveTime;

  return expected_observed_chis( m_peaks, m_backgroundPeaks, energy_count_map, info );
}//vector<tuple<double,double,double> > energy_chi_contributions(...) const


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
