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

#include "SpecUtils/DateTime.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/Integrate.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace GammaInteractionCalc
{

const double PointSourceShieldingChi2Fcn::sm_activityUnits = SandiaDecay::MBq;

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
    using namespace std;
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
    SelfAttCalc ObjectToIntegrate;
  
    const double energy = 185.0*PhysicalUnits::keV;
    ObjectToIntegrate.m_sourceIndex = 0;
    ObjectToIntegrate.m_detectorRadius  = 2.0 * PhysicalUnits::cm;
    ObjectToIntegrate.m_observationDist = 400.0 * PhysicalUnits::cm;

    double sphereRad = 0.0, transLenCoef = 0.0;

    const Material *material = materialdb.material( "void" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 99.5* PhysicalUnits::cm;
    ObjectToIntegrate.m_sphereRadAndTransLenCoef.push_back( std::pair<double,double>(sphereRad,transLenCoef) );
  
    material = materialdb.material( "U" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 0.5 * PhysicalUnits::cm;
    ObjectToIntegrate.m_sphereRadAndTransLenCoef.push_back( std::pair<double,double>(sphereRad,transLenCoef) );
    ObjectToIntegrate.m_sourceIndex = ObjectToIntegrate.m_sphereRadAndTransLenCoef.size() - 1;

    material = materialdb.material( "Fe" );
    transLenCoef = GammaInteractionCalc::transmition_length_coefficient( material, energy );
    sphereRad += 0.5 * PhysicalUnits::cm;
    ObjectToIntegrate.m_sphereRadAndTransLenCoef.push_back( std::pair<double,double>(sphereRad,transLenCoef) );


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
  Integrate::CuhreIntegrate( ndim, Integrand, userdata, epsrel, epsabs,
                            Integrate::LastImportanceFcnt, mineval, maxeval, nregions, neval,
                            fail, integral, error, prob );

  printf("ndim=%d CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
      ndim, nregions, neval, fail);
  printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", integral, error, prob);
  printf("\n\n" );
  
  ndim = 3;
  Integrate::CuhreIntegrate( ndim, Integrand, userdata, epsrel, epsabs,
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
    Suave(ndim, ncomp, Integrand, userdata,
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
  

int Integrand( const int *ndim, const double xx[],
                      const int *ncomp, double ff[], void *userdata )
{
  const SelfAttCalc *objToIntegrate = (SelfAttCalc *)userdata;

  try
  {
    objToIntegrate->eval( xx, ndim, ff, ncomp );
  }catch( const std::runtime_error e )
  {
    cerr << "Caught runtime_error: " << e.what() << " setting evaluation to zero!" << endl;
    ff[0] = 0.0;
  }catch( std::exception e )
  {
    ff[0] = 0.0;
    cerr << "Caught exception: " << e.what() << " setting evaluation to zero!" << endl;
  }

  return 0;
}//int Integrand(...)


SelfAttCalc::SelfAttCalc()
{
  m_sourceIndex = 0;
  m_detectorRadius = -1.0;
  m_observationDist = -1.0;
  nuclide = NULL;
}//SelfAttCalc()

  
  
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
  
  

void SelfAttCalc::eval( const double xx[], const int *ndimptr,
                        double ff[], const int *ncompptr ) const
{
  using namespace std;
  const double pi = 3.14159265358979;

  const int ndim = (ndimptr ? (*ndimptr) : 3);

//  assert( ndim == 2 || ndim == 3 );
//  const int ncomp = (ncompptr ? *ncompptr : 1);
//  assert( ncomp==1 );
//  assert( m_sourceIndex == 0 );  //Right now the source must be the inner most material!
//  assert( m_sourceIndex >= 0 && m_sourceIndex < m_sphereRadAndTransLenCoef.size() );

  double source_inner_rad = ((m_sourceIndex>0)
                            ? m_sphereRadAndTransLenCoef[m_sourceIndex-1].first
                            : 0.0);
  const double surce_outer_rad = m_sphereRadAndTransLenCoef[m_sourceIndex].first;

  //cuba goes from 0 to one, so we have to scale the variables
  //r:     goes from zero to sphere_rad
  //theta: goes from 0 to pi
  //phi:   goes from 0 to 2pi

  double max_theta = pi;

  /*
   // XXX - 20180416 - It was found the efficiency saver in this comment block
   //       was actually making the answer for a 1 kg DU ball significantly
   //       wrong!  Revisit at some point or just delete it.
  //Lets not waste time integrating over area in the source which wont
  //  contribute to radiation making it outside the ball
  const double MAX_ATT_LEN = 7.5;
  const double srcTransLenCoef = m_sphereRadAndTransLenCoef[m_sourceIndex].second;
  if( (MAX_ATT_LEN*srcTransLenCoef) < surce_outer_rad )
  {
    source_inner_rad = std::max( source_inner_rad,
                                 surce_outer_rad-MAX_ATT_LEN*srcTransLenCoef );

    //Well use a kinda course approximation to limit theta so we dont bother
    //  integrating over the "back side" of the sphere when radiation from this
    //  area wont actually make it through the sphere
    //For a 10cm sphere of Uranium and a 2D integral, reduces neval from
    //  6435 to 6045; for 3D integral, reduce neval from 16383 to 15621.
    const double opposite = 0.5*MAX_ATT_LEN*srcTransLenCoef;
    const double hypotenuse = surce_outer_rad;
    const double angle = asin( opposite/hypotenuse );
    max_theta = 0.5*pi + angle;
  }//if( not the whole source voluyme will contribute )
  */
  
  const double x_r = xx[0];
  const double x_theta = xx[1];
  const double x_phi = (ndim<3 ? 0.5 : xx[2]);
  
  const double r = source_inner_rad + x_r * (surce_outer_rad-source_inner_rad);
  const double theta =  x_theta * max_theta;
  
  //phi should go from zero to two pi, but we could reduce this to just 0 to
  //  pi, due to symetrty, but not doing right now (cause I'm to lazy to test it)
  const double max_phi = 2.0 * pi;
  const double phi = x_phi * max_phi;
  const double j = (surce_outer_rad - source_inner_rad) * max_theta * max_phi;
  const double dV = j * r * r * sin(theta);

  double source_point[3];
  source_point[0] = r * sin( theta ) * cos( phi );
  source_point[1] = r * sin( theta ) * sin( phi );
  source_point[2] = r * cos( theta );

  double trans = 0.0;
  
  {//begin codeblock to compute distance through source
    // - this could probably be cleaned up and made more effieicent
    double exit_point[3];
    const double srcRad = m_sphereRadAndTransLenCoef[m_sourceIndex].first;
    const double srcTransCoef = m_sphereRadAndTransLenCoef[m_sourceIndex].second;
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
      //  than the source point.
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
      const double subrad = m_sphereRadAndTransLenCoef[m_sourceIndex-1].first;
      needShellCompute = ( (min_rad < subrad)
                           && (inner_shell_point[2] > source_point[2]) );
    }//if( m_sourceIndex > 0 )
    
    if( needShellCompute )
    {
      double original_point[3];
      memcpy( original_point, source_point, 3*sizeof(double) );
      memcpy( source_point, inner_shell_point, 3*sizeof(double) );
      
      //Calc how far from the gammas original position it was, to the first
      //  inner shell it hit
      const double innerRad = m_sphereRadAndTransLenCoef[m_sourceIndex-1].first;
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
      while( m_sphereRadAndTransLenCoef[start_index].first < min_rad
            && start_index < m_sphereRadAndTransLenCoef.size() )
        ++start_index;
      
      //Some hopefully un-needed logic checks
      if( start_index == m_sphereRadAndTransLenCoef.size() )
        throw runtime_error( "Logic error 1 in SelfAttCalc::eval(...)" );
      if( start_index >= m_sourceIndex )
        throw runtime_error( "Logic error 2 in SelfAttCalc::eval(...)" );
      if( m_sourceIndex == 0 )
        throw runtime_error( "Logic error 3 in SelfAttCalc::eval(...)" );
      
      //calc distance it travels through the inner spheres
      for( size_t index = start_index; index <= m_sourceIndex; ++index )
      {
        const double shellRad = m_sphereRadAndTransLenCoef[index].first;
        const double transCoef = m_sphereRadAndTransLenCoef[index].second;
        double dist = exit_point_of_sphere_z( source_point, source_point,
                                                  shellRad, m_observationDist );
        if( index != m_sourceIndex )
          dist = 2.0*dist;
        
        trans += (transCoef * dist);
      }//for( ++sphere_index; sphere_index < m_sourceIndex; ++sphere_index )
    }else
    {
      memcpy( source_point, exit_point, 3*sizeof(double) );
      trans += (srcTransCoef * dist_in_src);
    }//if( line actually goes into daughter sphere ) / else
  }//end codeblock to compute distance through source
  
  for( size_t i = m_sourceIndex+1; i < m_sphereRadAndTransLenCoef.size(); ++i )
  {
    const double sphereRad = m_sphereRadAndTransLenCoef[i].first;
    const double transLenCoef = m_sphereRadAndTransLenCoef[i].second;
    const double dist_in_sphere = exit_point_of_sphere_z( source_point,
                                   source_point, sphereRad, m_observationDist );
    trans += (transLenCoef * dist_in_sphere);
  }//for( size_t i = 0; i < m_sphereRadAndTransLenCoef.size(); ++i )

  trans = exp( -trans );
  trans *= DetectorPeakResponse::fractionalSolidAngle( 2.0*m_detectorRadius,
                                                       m_observationDist );

  ff[0] = trans * dV;
}//eval(...)


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

PointSourceShieldingChi2Fcn::PointSourceShieldingChi2Fcn(
                                 double distance, double liveTime,
                                 const std::vector<PeakDef> &peaks,
                                 std::shared_ptr<const DetectorPeakResponse> detector,
                                 const std::vector<PointSourceShieldingChi2Fcn::MaterialAndSources> &materials,
                                 bool allowMultipleNucsContribToPeaks )
  : ROOT::Minuit2::FCNBase(),
    m_cancel( 0 ),
    m_isFitting( false ),
    m_guiUpdateFrequencyMs( 0 ),
    m_distance( distance ),
    m_liveTime( liveTime ),
    m_photopeakClusterSigma( 1.25 ),
    m_peaks( peaks ),
    m_detector( detector ),
    m_materials( materials ),
    m_nuclides( 0, (const SandiaDecay::Nuclide *)NULL ),
    m_allowMultipleNucsContribToPeaks( allowMultipleNucsContribToPeaks )
{
  set<const SandiaDecay::Nuclide *> nucs;
  for( const PeakDef &p : m_peaks )
    if( p.parentNuclide() )
    nucs.insert( p.parentNuclide() );
  for( const SandiaDecay::Nuclide *n : nucs )
    m_nuclides.push_back( n );

  std::sort( m_nuclides.begin(), m_nuclides.end(),
            []( const SandiaDecay::Nuclide *lhs, const SandiaDecay::Nuclide *rhs ) -> bool {
              if( !lhs ) return false;
              if( !rhs ) return true;
              return (lhs->symbol < rhs->symbol);
            }  );
}//PointSourceShieldingChi2Fcn


PointSourceShieldingChi2Fcn::~PointSourceShieldingChi2Fcn()
{
  {//begin lock on m_zombieCheckTimerMutex
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
    
    if( m_zombieCheckTimer )
    {
      boost::system::error_code ec;
      m_zombieCheckTimer->cancel( ec );
      if( ec )
        cerr << "PointSourceShieldingChi2Fcn::fittindIsFinished(): Got error cancelling m_zombieCheckTimer: "
        << ec.message() << endl;
      m_zombieCheckTimer.reset();
    }//if( m_zombieCheckTimer )
  }//end lock on m_zombieCheckTimerMutex
}

void PointSourceShieldingChi2Fcn::cancelFit()
{
  m_cancel = 1;
}


void PointSourceShieldingChi2Fcn::setGuiProgressUpdater( const size_t updateFrequencyMs,
                                                  std::shared_ptr<GuiProgressUpdateInfo> updateInfo )

{
  m_guiUpdateFrequencyMs = updateFrequencyMs;
  m_guiUpdateInfo = updateInfo;
}//void setGuiUpdater(...)
  
  
void PointSourceShieldingChi2Fcn::zombieCallback( const boost::system::error_code &ec )
{
  if( ec )  //timer was cancelled
    return;
  
  m_cancel = 2;

  cerr << "In Zombie callback!" << endl;
}//void zombieCallback( const boost::system::error_code &ec )
  
  
void PointSourceShieldingChi2Fcn::fittingIsStarting( const size_t deadlineMs )
{
  m_isFitting = true;
  
  if( deadlineMs )
  {
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
    m_zombieCheckTimer = make_shared<boost::asio::deadline_timer>( Wt::WServer::instance()->ioService() );
    m_zombieCheckTimer->expires_from_now( boost::posix_time::milliseconds(deadlineMs) );
    m_zombieCheckTimer->async_wait( [this](const boost::system::error_code &ec){ zombieCallback(ec); } );
  }//if( deadlineMs )

  if( m_guiUpdateInfo )
  {
    std::lock_guard<std::mutex> scoped_lock( m_guiUpdateInfo->m_bestParsMutex );
    m_guiUpdateInfo->m_bestChi2 = std::numeric_limits<double>::max();
    m_guiUpdateInfo->m_num_fcn_calls = 0;
    m_guiUpdateInfo->m_bestParameters.clear();
    m_guiUpdateInfo->m_lastGuiUpdateTime = SpecUtils::get_wall_time();
    m_guiUpdateInfo->m_fitStartTime = m_guiUpdateInfo->m_lastGuiUpdateTime.load();
    m_guiUpdateInfo->m_currentTime = m_guiUpdateInfo->m_lastGuiUpdateTime.load();
  }
}//void fittingIsStarting()

  
void PointSourceShieldingChi2Fcn::fittindIsFinished()
{
  {//begin lock on m_zombieCheckTimerMutex
    std::lock_guard<std::mutex> scoped_lock( m_zombieCheckTimerMutex );
  
    if( m_zombieCheckTimer )
    {
      boost::system::error_code ec;
      m_zombieCheckTimer->cancel( ec );
      if( ec )
        cerr << "PointSourceShieldingChi2Fcn::fittindIsFinished(): Got error cancelling m_zombieCheckTimer: "
             << ec.message() << endl;
      m_zombieCheckTimer.reset();
    }//if( m_zombieCheckTimer )
  }//end lock on m_zombieCheckTimerMutex
  
  m_isFitting = false;
}//void fittindIsFinished()
  
  
const SandiaDecay::Nuclide *PointSourceShieldingChi2Fcn::nuclide( const int nuc ) const
{
  if( nuc < 0 || nuc > static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "PointSourceShieldingChi2Fcn::nuclide(int): invalid index" );
  return m_nuclides[nuc];
}//const Nuclide *nuclide( const int nucN ) const


double PointSourceShieldingChi2Fcn::operator()( const std::vector<double> &x ) const
{
  return DoEval( x );
}//double operator()(...)


PointSourceShieldingChi2Fcn &PointSourceShieldingChi2Fcn::operator=( const PointSourceShieldingChi2Fcn &rhs )
{
  m_cancel = rhs.m_cancel.load();
  m_distance = rhs.m_distance;
  m_liveTime = rhs.m_liveTime;
  m_photopeakClusterSigma = rhs.m_photopeakClusterSigma;
  m_peaks = rhs.m_peaks;
  m_detector = rhs.m_detector;
  m_materials = rhs.m_materials;
  m_nuclides = rhs.m_nuclides;
  m_allowMultipleNucsContribToPeaks = rhs.m_allowMultipleNucsContribToPeaks;
  m_nuclidesToFitMassFractionFor    = rhs.m_nuclidesToFitMassFractionFor;
  
  return *this;
}//operator=


double PointSourceShieldingChi2Fcn::Up() const
{
  return 1.0;
}//double Up();

#if( ALLOW_MASS_FRACTION_FIT )
  
bool PointSourceShieldingChi2Fcn::isVariableMassFraction( const Material *mat,
                                        const SandiaDecay::Nuclide *nuc ) const 
{
  MaterialToNucsMap::const_iterator pos
                                   = m_nuclidesToFitMassFractionFor.find( mat );
  if( pos == m_nuclidesToFitMassFractionFor.end() )
    return false;
  const vector<const SandiaDecay::Nuclide *> &nucs = pos->second;
  return ( std::find( nucs.begin(), nucs.end(), nuc) != nucs.end() );
}//isVariableMassFraction(...)
  
  
bool PointSourceShieldingChi2Fcn::hasVariableMassFraction(
                                                const Material *mat ) const
{
  MaterialToNucsMap::const_iterator pos
                                   = m_nuclidesToFitMassFractionFor.find( mat );
  return (pos != m_nuclidesToFitMassFractionFor.end());
}//bool hasVariableMassFraction( const Material *material ) const



std::shared_ptr<Material> PointSourceShieldingChi2Fcn::variedMassFracMaterial(
                                          const Material *material,
                                          const std::vector<double> &x ) const
{
  if( !material )
    throw runtime_error( "variedMassFracMaterial(): invalid input material" );
  
  if( !hasVariableMassFraction(material) )
    throw runtime_error( "variedMassFracMaterial(): " + material->name
                         + " does not have a variable mass fraction nuclide" );
  
  std::shared_ptr<Material> answer = std::make_shared<Material>( *material );
  MaterialToNucsMap::const_iterator iter = m_nuclidesToFitMassFractionFor.find( material );
  
  if( iter == end(m_nuclidesToFitMassFractionFor) )
    throw runtime_error( "variedMassFracMaterial(): material not found in m_nuclidesToFitMassFractionFor" );
  
  const vector<const SandiaDecay::Nuclide *> &nucs = iter->second;
  
  double prefrac = 0.0, postfrac = 0.0;
  for( const SandiaDecay::Nuclide *nuc : nucs )
  {
    const double frac = massFraction( material, nuc, x );
    
    for( Material::NuclideFractionPair &nfp : answer->nuclides )
    {
      if( nfp.first == nuc )
      {
        prefrac += nfp.second;
        postfrac += frac;
        nfp.second = frac;
      }//if( nfp.first == nuc )
    }//for( Material::NuclideFractionPair &nfp : answer->nuclides )
  }//for( const SandiaDecay::Nuclide *nuc : nucs )
  
  if( IsNan(postfrac) || IsInf(postfrac)
      || ((fabs(prefrac-postfrac)/max(prefrac,postfrac)) > 0.01) )
  {
    cerr << "PointSourceShieldingChi2Fcn::variedMassFracMaterial: prefrac="
          << prefrac << ", postfrac=" << postfrac << endl;

    throw runtime_error( "PointSourceShieldingChi2Fcn::variedMassFracMaterial(...)"
                         " prefrac did not match postfrac" );
  }//if( invalid results )
  
  return answer;
}//variedMassFracMaterial(...)


void PointSourceShieldingChi2Fcn::setNuclidesToFitMassFractionFor(
                          const Material *material,
                          const vector<const SandiaDecay::Nuclide *> &nuclides )
{
  if( !material )
    throw runtime_error( "setNuclidesToFitMassFractionFor(...): "
                         "invalid material" );
  bool validMaterial = false;
  for( const MaterialAndSources &ms : m_materials )
    validMaterial |= (ms.first == material);
  
  if( !validMaterial )
    throw runtime_error( "setNuclidesToFitMassFractionFor(...): "
                        "material passed in is not a shielding material" );
  
  for( const SandiaDecay::Nuclide *nuc : nuclides )
  {
    if( std::count( nuclides.begin(), nuclides.end(), nuc ) != 1 )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): input must"
                           " only be unique nuclides" );
    int innuc = 0;
    for( const Material::NuclideFractionPair &nfp : material->nuclides )
      innuc += int(nfp.first == nuc);
    if( !innuc )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): passed in"
                           " nuclide " + nuc->symbol + " not in material "
                           + material->name );
    if( std::count( nuclides.begin(), nuclides.end(), nuc ) != 1 )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): input must"
                          " only be unique nuclides" );
    if( std::count( m_nuclides.begin(), m_nuclides.end(), nuc ) != 1 )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): input must "
                           "be in source m_nuclides" );
    if( !isActivityDefinedFromShielding( nuc ) )
      throw runtime_error( "setNuclidesToFitMassFractionFor(...): "
                           + nuc->symbol
                           + " is not in a self-attenuating shielding" );
    
  }//for( const SandiaDecay::Nuclide *nuc : nuclides )
  
  m_nuclidesToFitMassFractionFor[material] = nuclides;
  
  std::sort( m_nuclidesToFitMassFractionFor[material].begin(),
             m_nuclidesToFitMassFractionFor[material].end(),
    []( const SandiaDecay::Nuclide *lhs, const SandiaDecay::Nuclide *rhs ) -> bool {
      if( !lhs ) return false;
      if( !rhs ) return true;
      return (lhs->symbol < rhs->symbol);
    } );
}//setNuclidesToFitMassFractionFor(...)
  
vector<const Material *>
               PointSourceShieldingChi2Fcn::materialsFittingMassFracsFor() const
{
  vector<const Material *> answer;
  for( const MaterialToNucsMap::value_type &n : m_nuclidesToFitMassFractionFor )
  {
    if( !n.second.empty() )
      answer.push_back( n.first );
  }
  return answer;
}//vector<const Material *material> materialsFittingMassFracsFor() const
  
  
const std::vector<const SandiaDecay::Nuclide *> &
PointSourceShieldingChi2Fcn::nuclideFittingMassFracFor( const Material *material ) const
{
  if( !material )
    throw runtime_error( "PointSourceShieldingChi2Fcn::nuclideFittingMassFracFor():"
                        " invalid input" );

  typedef MaterialToNucsMap::const_iterator MaterialToNucsMapIter;
  MaterialToNucsMapIter matpos = m_nuclidesToFitMassFractionFor.find(material);
    
  if( matpos == m_nuclidesToFitMassFractionFor.end() )
    throw runtime_error( "PointSourceShieldingChi2Fcn::nuclideFittingMassFracFor(): "
                        + material->name + " is not a material with a variable"
                        " mass fraction" );
  return matpos->second;
}//nuclideFittingMassFracFor(...)

double PointSourceShieldingChi2Fcn::massFraction( const Material *material,
                      const SandiaDecay::Nuclide *nuc,
                      const std::vector<double> &pars ) const
{
  double massfrac = 0.0, uncert = 0.0;
  massFraction( massfrac, uncert, material, nuc, pars, vector<double>() );
  return massfrac;
}
  
void PointSourceShieldingChi2Fcn::massFraction( double &massFrac,
                                                   double &uncert,
                                                  const Material *material,
                                        const SandiaDecay::Nuclide *nuc,
                                        const std::vector<double> &pars,
                                        const std::vector<double> &errors ) const
{
  massFrac = uncert = 0.0;
  typedef MaterialToNucsMap::const_iterator MaterialToNucsMapIter;
  if( !material || !nuc )
    throw runtime_error( "PointSourceShieldingChi2Fcn::massFraction(): invalid "
                         "input" );
  MaterialToNucsMapIter matpos = m_nuclidesToFitMassFractionFor.find(material);
  
  if( matpos == m_nuclidesToFitMassFractionFor.end() )
    throw runtime_error( "PointSourceShieldingChi2Fcn::massFraction(): "
                         + material->name + " is not a material with a variable"
                         " mass fraction" );
  const std::vector<const SandiaDecay::Nuclide *> &nucs
                        = m_nuclidesToFitMassFractionFor.find(material)->second;
  vector<const SandiaDecay::Nuclide *>::const_iterator pos;
  //nucs is actually sorted by symbol name, could do better than linear search
  pos = std::find( nucs.begin(), nucs.end(), nuc );
  if( pos == nucs.end() )
    throw runtime_error( "PointSourceShieldingChi2Fcn::massFraction(): "
                         +nuc->symbol + " was not a nuc fit for mass fraction"
                         " in " + material->name );
  
  double totalfrac = 0.0;
  for( const SandiaDecay::Nuclide *n : nucs )
  {
    for( const Material::NuclideFractionPair &nfp : material->nuclides )
    {
      if( nfp.first == n )
      {
        totalfrac += nfp.second;
      }//if( nfp.first == src )
    }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )
  }//for( const SandiaDecay::Nuclide *n : nucs )
  
  size_t matmassfracstart = 0;
  for( MaterialToNucsMapIter iter = m_nuclidesToFitMassFractionFor.begin();
      iter != matpos; ++iter )
  {
    if( iter->second.size() )
      matmassfracstart += (iter->second.size()-1);
  }
  
  matmassfracstart += 2 * m_nuclides.size();
#if(USE_CONSISTEN_NUM_SHIELDING_PARS)
  matmassfracstart += 2*m_materials.size();
#else
  for( const MaterialAndSources &material : m_materials )
    matmassfracstart += (material.first ? size_t(1) : size_t(2));
#endif
  
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
    throw runtime_error( "PointSourceShieldingChi2Fcn::massFraction(...): invalid parameters" );
  }
  
  if( !errors.empty() )
  {
    if( errors.size() != pars.size() )
      throw runtime_error( "PointSourceShieldingChi2Fcn::massFraction():"
                           " invalid error parameter vector size" );
    double fracuncert = 0.0;
    if( massfracnum != numfraccoefs )
      fracuncert = errors.at(thismatmassfrac) / pars.at(thismatmassfrac);
    else
      fracuncert = errors.at(thismatmassfrac-1) / pars.at(thismatmassfrac-1);
    uncert = fracuncert * massFrac;
  }//if( !errors.empty() )
}//double massFraction(...) const

  
double PointSourceShieldingChi2Fcn::massFractionUncert( const Material *material,
                            const SandiaDecay::Nuclide *nuc,
                            const std::vector<double> &pars,
                            const std::vector<double> &error ) const
{
  double massfrac = 0.0, uncert = 0.0;
  massFraction( massfrac, uncert, material, nuc, pars, error );
  return uncert;
}//massFractionUncert(...)(
#endif
  
  

size_t PointSourceShieldingChi2Fcn::numExpectedFitParameters() const
{
  size_t npar = 2 * m_nuclides.size();

#if(USE_CONSISTEN_NUM_SHIELDING_PARS)
  npar += 2*m_materials.size();
#else
  for( const MaterialAndSources &material : m_materials )
    npar += (material.first ? size_t(1) : size_t(2));
#endif
  
#if( ALLOW_MASS_FRACTION_FIT )
  for( const MaterialToNucsMap::value_type &vt : m_nuclidesToFitMassFractionFor )
    npar += ( vt.second.empty() ? size_t(0) : size_t(vt.second.size()-1) );
#endif
  
  return npar;
}//int numExpectedFitParameters() const


bool PointSourceShieldingChi2Fcn::isActivityDefinedFromShielding( const SandiaDecay::Nuclide *nuclide ) const
{
  //XXX - this is fairly comutationally ineficienct
  if( !nuclide )
    return false;

  for( const MaterialAndSources &material : m_materials )
  {
    for( const SandiaDecay::Nuclide *nuc : material.second )
    {
      if( nuc == nuclide )
        return true;
    }
  }//for( const MaterialAndSources &material : m_materials )

  return false;
}//bool isActivityDefinedFromShielding(...) const;


double PointSourceShieldingChi2Fcn::activityDefinedFromShielding(
                                       const SandiaDecay::Nuclide *nuclide,
                                       const std::vector<double> &params ) const
{
  double radius = 0.0, activity = 0.0;

  const int nmat = static_cast<int>(numMaterials());
  for( int matn = 0; matn < nmat; ++matn )
  {
    if( isGenericMaterial( matn ) )
      continue;

    const double thick = thickness( matn, params );
    const Material *mat = m_materials[matn].first;
    
//    cerr << mat->name << "(matn=" << matn << "), thick=" << thick/PhysicalUnits::cm << endl;
//    cerr << "\tshould be" << params.at( 2*m_nuclides.size() + 2*matn ) << endl;
    
    for( const SandiaDecay::Nuclide *src : m_materials[matn].second )
    {
      if( src == nuclide )
      {
        double massFrac = 0.0;
        
#if( ALLOW_MASS_FRACTION_FIT )
        if( hasVariableMassFraction(mat) )
        {
          massFrac = PointSourceShieldingChi2Fcn::massFraction( mat, src, params );
        }else
#endif
        {
          for( const Material::NuclideFractionPair &nfp : mat->nuclides )
          {
            if( nfp.first == src )
              massFrac += nfp.second;
          }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )
        }//if( hasVariableMassFraction(mat) ) / else

/*
        if( massFrac <= 0.0 )
        {
          //This condition can happen during fitting if the fit goes down to
          //  a mass fraction of zero, so its not an absolute error, just not
          //  ever practically expected.
          stringstream msg;
          msg << "Potential serious error in activityDefinedFromShielding for "
              << src->symbol << ", massFrac=" << massFrac;
          cerr << msg.str() << endl;
#if( PERFORM_DEVELOPER_CHECKS )
          log_developer_error( __func__, msg.str().c_str() );
#endif
        }//if( massFrac <= 0.0 )
*/
        
        const double vol = (4.0/3.0)*PhysicalUnits::pi*( pow(radius+thick,3.0) - pow(radius,3.0));
        const double mass_grams = massFrac * mat->density * vol / PhysicalUnits::gram;
        const double activity_per_gram = nuclide->activityPerGram();
        activity += mass_grams * activity_per_gram;
        
//        cerr << "activityDefinedFromShielding: " << nuclide->symbol << " (matn=" << matn
//         << "), mat=" << mat->name
//        << ", thick=" << thick/PhysicalUnits::cm
//            << " cm, inner_rad=" << radius/PhysicalUnits::cm
//        << " cm, vol=" << vol/(PhysicalUnits::cm*PhysicalUnits::cm*PhysicalUnits::cm)
//        << " cm3, density=" << mat->density*(PhysicalUnits::cm*PhysicalUnits::cm*PhysicalUnits::cm)/PhysicalUnits::g
//        << " g/cm3, mass=" << mass_grams << " g, activity_per_gram=" << activity_per_gram/PhysicalUnits::ci
//        << " ci, mass_grams * activity_per_gram=" << mass_grams * activity_per_gram/PhysicalUnits::ci << " ci" << endl;
      }//if( src == nuclide )
    }//for( const SandiaDecay::Nuclide *src : srcs )

    radius += thick;
  }//for( size_t matn = 0; matn < nmat; ++matn )

//  cerr << "activityDefinedFromShielding: " << nuclide->symbol << " activity=" << activity/PhysicalUnits::ci << endl;
  
  return activity;
}//bool activityDefinedFromShielding(...) const;




double PointSourceShieldingChi2Fcn::activity( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
  {
    cerr << "params.size()=" << params.size()
         << ", numExpectedFitParameters()=" << numExpectedFitParameters()
         << endl;
    throw runtime_error( "PointSourceShieldingChi2Fcn::activity(...) invalid params size" );
  }//if( params.size() != numExpectedFitParameters() )


  if( isActivityDefinedFromShielding( nuclide ) )
    return activityDefinedFromShielding( nuclide, params );
  
  const size_t ind = nuclideIndex( nuclide );  
  return params[2*ind] * sm_activityUnits;
}//double activity(...)



double PointSourceShieldingChi2Fcn::activityUncertainty( const SandiaDecay::Nuclide *nuclide,
                                                       const std::vector<double> &params,
                                                        const std::vector<double> &errors ) const
{
  if( params.size() != numExpectedFitParameters() )
  {
    cerr << "params.size()=" << params.size()
    << ", numExpectedFitParameters()=" << numExpectedFitParameters()
    << endl;
    throw runtime_error( "PointSourceShieldingChi2Fcn::activityUncertainty(...) invalid params size" );
  }//if( params.size() != numExpectedFitParameters() )
  
  if( params.size() != errors.size() )
    throw runtime_error( "PointSourceShieldingChi2Fcn::activityUncertainty: size(params) != size(errors)" );
  
  // If it is a simple point source, we can get the uncertainty by just passing in the
  //  parameter uncertainties and calculating activity from those.  However, if its a
  //  self-attenuating source, we have to do a little more propagation and such.
  
  if( !isActivityDefinedFromShielding( nuclide ) )
    return activity( nuclide, errors );
    
  
  // We will take into account the uncertainties of the inner layers to our current shell, the
  //  uncertainty of the thickness, and the uncertainty of the mass fractions.
  double radius = 0.0, radiusUncertSquared = 0.0, activity = 0.0, activityUncertSquared = 0.0;
  const int nmat = static_cast<int>(numMaterials());
  for( int matn = 0; matn < nmat; ++matn )
  {
    if( isGenericMaterial( matn ) )
      continue;
    
    const double thick = thickness( matn, params );
    const double thickUncert = thickness( matn, errors );
    const Material * const mat = m_materials[matn].first;
    
    for( const SandiaDecay::Nuclide *src : m_materials[matn].second )
    {
      if( src == nuclide )
      {
        double massFrac = 0.0, massFracUncert = 0.0;
        
#if( ALLOW_MASS_FRACTION_FIT )
        if( hasVariableMassFraction(mat) )
        {
          massFraction( massFrac, massFracUncert, mat, src, params, errors );
        }else
#endif
        {
          for( const Material::NuclideFractionPair &nfp : mat->nuclides )
          {
            if( nfp.first == src )
              massFrac += nfp.second;
          }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )
        }//if( hasVariableMassFraction(mat) ) / else

        const double outerRad = radius + thick;
        
        const double iso_density = massFrac * mat->density / PhysicalUnits::gram;
        const double activity_per_gram = nuclide->activityPerGram();
        
        const double volumeUncertDueToInnerRad = 4.0 * PhysicalUnits::pi * pow(radius,2.0) * sqrt(radiusUncertSquared);
        
        // We take the uncertainty in volume due to the thickness uncertainty to be independent of
        //  the inner-radius uncertainty - which is approximately, probably, about right to kinda
        //  fairly cover all the uncertainties - maybe
        const double volumeUncertDueToThickness = 4.0 * PhysicalUnits::pi * pow(outerRad,2.0) * thickUncert;
        
        const double vol = (4.0/3.0)*PhysicalUnits::pi*( pow(outerRad,3.0) - pow(radius,3.0));
        const double volUncertSquared = pow(volumeUncertDueToInnerRad,2.0) + pow(volumeUncertDueToThickness,2.0);
        const double mass_grams = iso_density * vol;
        const double massUncertaintySquared_grams = iso_density * volUncertSquared;
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
      }//if( src == nuclide )
    }//for( const SandiaDecay::Nuclide *src : srcs )
    
    radius += thick;
    
    // For the layer of shielding we are concerned with, we will assume the uncertainty on the
    //  inner radius is the squared sum of all inner layer uncertainties - this is actually a bit
    //  over the top - these thickness uncertainties are likely very correlated - so much so we
    //  could probably use just the thickness uncertainty on just the previous layer, but we'll be
    //  conservative, or something.
    const double radiusUncert = thickness( matn, errors );
    radiusUncertSquared += radiusUncert*radiusUncert;
  }//for( size_t matn = 0; matn < nmat; ++matn )
  
  if( activityUncertSquared < 0.0 || IsNan(activityUncertSquared) || IsInf(activityUncertSquared) )
    throw runtime_error( "error calculating activity uncertainty for self-attenuating source;"
                         " squared value calculated is " + std::to_string(activityUncertSquared) );
  
  const double normalCalcAct = PointSourceShieldingChi2Fcn::activity( nuclide, params );
  assert( fabs(normalCalcAct - activity) < 0.001*std::max(normalCalcAct,activity) );
  
  return sqrt( activityUncertSquared );
}//double activityUncertainty(...)



double PointSourceShieldingChi2Fcn::age( const SandiaDecay::Nuclide *nuclide,
                                      const std::vector<double> &params ) const
{
  if( params.size() != numExpectedFitParameters() )
    throw runtime_error( "PointSourceShieldingChi2Fcn::age(...) invalid params size" );

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


size_t PointSourceShieldingChi2Fcn::nuclideIndex( const SandiaDecay::Nuclide *nuclide ) const
{
  vector<const SandiaDecay::Nuclide *>::const_iterator pos;
  pos = find( m_nuclides.begin(), m_nuclides.end(), nuclide );
  if( pos == m_nuclides.end() )
    throw std::runtime_error( "PointSourceShieldingChi2Fcn::nuclideIndex(...): invalid nuclide" );

  return (pos - m_nuclides.begin());
}//size_t index( const SandiaDecay::Nuclide *nuclide ) const

  
size_t PointSourceShieldingChi2Fcn::numNuclides() const
{
  return m_nuclides.size();
}
  
size_t PointSourceShieldingChi2Fcn::numMaterials() const
{
  return m_materials.size();
}//int numMaterials() const


const Material *PointSourceShieldingChi2Fcn::material( int materialNum ) const
{
  return m_materials.at(materialNum).first;
}

bool PointSourceShieldingChi2Fcn::isSpecificMaterial( int materialNum ) const
{
  return (NULL != m_materials.at( materialNum ).first );
}//bool isSpecificMaterial( int materialNum ) const


bool PointSourceShieldingChi2Fcn::isGenericMaterial( int materialNum ) const
{
  return (NULL == m_materials.at( materialNum ).first );
}//bool isGenericMaterial( int materialNum ) const

  
bool PointSourceShieldingChi2Fcn::hasSelfAttenuatingSource() const
{
  const size_t nmat = numMaterials();
  for( size_t i = 0; i < nmat; ++i )
    if( materialIsSelfAttenuatingSource( static_cast<int>(i) ) )
      return true;
  return false;
}//bool hasSelfAttenuatingSource() const


bool PointSourceShieldingChi2Fcn::materialIsSelfAttenuatingSource( int materialNum ) const
{
  return !m_materials.at(materialNum).second.empty();
}//bool materialIsSelfAttenuatingSource( int materialNum ) const;
  

vector<const SandiaDecay::Nuclide *> PointSourceShieldingChi2Fcn::sourceNulcidesForShielding( int materialNum ) const
{
  return m_materials.at(materialNum).second;
}//vector<SandiaDecay::Nuclide *> sourceNulcidesForShielding( int materialNum ) const


//thickness(...): will throw std::runtime_exception if material is a generic
//  material
double PointSourceShieldingChi2Fcn::thickness( int materialNum,
                    const std::vector<double> &params ) const
{
  if( isGenericMaterial( materialNum ) )
    throw std::runtime_error( "PointSourceShieldingChi2Fcn::thickness(...): "
                              "You should not call this function for generic "
                              "materials" );

#if( USE_CONSISTEN_NUM_SHIELDING_PARS )
  return params.at( 2*m_nuclides.size() + 2*materialNum );
#else
  size_t index = 2*m_nuclides.size();
  for( int n = 0; n < materialNum; ++n )
    index += (m_materials[n].first ? 1 : 2);
  return params.at( index );
#endif
}//double thickness(...)


//arealDensity(...): will throw std::runtime_exception if material is a
//  specific material
double PointSourceShieldingChi2Fcn::arealDensity( int materialNum,
                       const std::vector<double> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "PointSourceShieldingChi2Fcn::arealDensity(...): "
                              "You can only call this function for generic "
                              "materials" );
#if( USE_CONSISTEN_NUM_SHIELDING_PARS )
  return params.at( 2*m_nuclides.size() + 2*materialNum + 1 );
#else
  size_t index = 2*m_nuclides.size();
  for( int n = 0; n < materialNum; ++n )
    index += (m_materials[n].first ? 1 : 2);
  return params.at( index + 1 );
#endif
}//double arealDensity(...)


//atomicNumber(...): will throw std::runtime_exception if material is a
//  specific material
double PointSourceShieldingChi2Fcn::atomicNumber( int materialNum,
                       const std::vector<double> &params ) const
{
  if( !isGenericMaterial( materialNum ) )
    throw std::runtime_error( "PointSourceShieldingChi2Fcn::atomicNumber(...): "
                              "You can only call this function for generic "
                              "materials" );
#if( USE_CONSISTEN_NUM_SHIELDING_PARS )
  return params.at( 2*m_nuclides.size() + 2*materialNum );
#else
  size_t index = 2*m_nuclides.size();
  for( int n = 0; n < materialNum; ++n )
    index += (m_materials[n].first ? 1 : 2);
  return params.at( index );
#endif
}//double atomicNumber(...)


const std::vector<PeakDef> &PointSourceShieldingChi2Fcn::peaks() const
{
  return m_peaks;
}


double PointSourceShieldingChi2Fcn::DoEval( const std::vector<double> &x ) const
{
  const int cancelCode = m_cancel.load();
  switch( cancelCode )
  {
    case 0:  break;
    case 1:  throw runtime_error( "User cancelled" );            break;
    case 2:  throw runtime_error( "Calculation took to long." ); break;
    default: throw runtime_error( "Other evaluation error, code " + std::to_string(cancelCode) );
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
      = energy_chi_contributions( x, m_mixtureCache,
                               m_allowMultipleNucsContribToPeaks );
    double chi2 = 0.0;
    
    const size_t npoints = chi2s.size();
    for( size_t i = 0; i < npoints; ++i )
    chi2 += pow( std::get<1>(chi2s[i]), 2.0 );
    
    if( m_isFitting && m_guiUpdateFrequencyMs.load() && m_guiUpdateInfo )
    {
      //Need best chi2, and its cooresponding errors and
      m_guiUpdateInfo->m_num_fcn_calls++;
      const double prevBestChi2 = m_guiUpdateInfo->m_bestChi2.load();
      
      if( chi2 < prevBestChi2 )
      {
        std::lock_guard<std::mutex> scoped_lock( m_guiUpdateInfo->m_bestParsMutex );
        m_guiUpdateInfo->m_bestChi2 = chi2;
        m_guiUpdateInfo->m_bestParameters = x;
      }
      
      std::lock_guard<std::mutex> scoped_lock( m_guiUpdateInfo->m_guiUpdaterMutex );
      const double lastUpdateTime = m_guiUpdateInfo->m_lastGuiUpdateTime;
      const double currentTime = SpecUtils::get_wall_time();
      if( (currentTime - lastUpdateTime) > 0.001*m_guiUpdateFrequencyMs.load() )
      {
        m_guiUpdateInfo->m_lastGuiUpdateTime = currentTime;
        m_guiUpdateInfo->m_currentTime = currentTime;
        
        m_guiUpdateInfo->m_guiUpdater();
      }
    }//if( m_guiUpdateFrequencyMs.load() && m_guiUpdateInfo )
    
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

//ToDo: add ability to give summarry about 
void PointSourceShieldingChi2Fcn::cluster_peak_activities( std::map<double,double> &energy_count_map,
                                                           const std::vector< pair<double,double> > &energie_widths,
                                                           SandiaDecay::NuclideMixture &mixture,
                                                           const double act,
                                                           const double age,
                                                           const double photopeakClusterSigma,
                                                           const double energyToCluster,
                                                           vector<string> *info )
{
  typedef pair<double,double> DoublePair;

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

  const vector<SandiaDecay::EnergyRatePair> gammas
                = mixture.photons( age, SandiaDecay::NuclideMixture::OrderByEnergy );
  
  const vector<SandiaDecay::NuclideActivityPair> aged_activities
                                                      = mixture.activity( age );

  //The problem we have is that 'gammas' have the activity of the original
  //  parent ('nuclide') decreased by agining by 'age', however we want the
  //  parent to have 'sm_activityUnits' activity at 'age', so we we'll add a
  //  correction factor.
  if( mixture.numInitialNuclides() != 1 )
    throw runtime_error( "PointSourceShieldingChi2Fcn::cluster_peak_activities():"
                         " passed in mixture must have exactly one parent nuclide" );
  const SandiaDecay::Nuclide *nuclide = mixture.initialNuclide(0);
  
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
           << " m_photopeakClusterSigma=" << photopeakClusterSigma << endl;
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
             "PointSourceShieldingChi2Fcn::cluster_peak_activities(...)"
             " which is keeping this activity/shielding fit from happening "
             "(place a) - please complain to Will Johnson about this";
      passMessage( msg.str(),
                "PointSourceShieldingChi2Fcn", WarningWidget::WarningMsgHigh );
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


vector< tuple<double,double,double,Wt::WColor,double> > PointSourceShieldingChi2Fcn::expected_observed_chis(
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
             "PointSourceShieldingChi2Fcn::expected_observed_chis(...)"
             " which is keeping this activity/shielding fit from happening "
             "(place b) - please complain to Will Johnson about this";
      passMessage( msg.str(),
                   "PointSourceShieldingChi2Fcn", WarningWidget::WarningMsgHigh );
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


vector< pair<double,double> > PointSourceShieldingChi2Fcn::observedPeakEnergyWidths( const std::vector<PeakDef> &peaks )
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
      cerr << "PointSourceShieldingChi2Fcn::observedPeakEnergyWidths(...)\n\tIssue here" << endl;
    }
  }//for( const PeakDef &peak : peaks )

  sort( energie_widths.begin(), energie_widths.end(), &first_lessthan );

  return energie_widths;
}//observedPeakEnergyWidths()


void PointSourceShieldingChi2Fcn::selfShieldingIntegration( SelfAttCalc &calculator )
{
  const int ndim = 2;  //the number of dimensions of the integral.
  void *userdata = (void *)&calculator;
  const double epsrel = 1e-4;  //the requested relative accuracy
  const double epsabs = -1.0;//1e-12; //the requested absolute accuracy
  //const int verbose = 0;
  //const int last = 4;  //use the importance function without smoothing (*I think*)
  const int mineval = 0; //the minimum number of integrand evaluations required.
  const int maxeval = 5000000; //the (approximate) maximum number of integrand evaluations allowed.

  int nregions, neval, fail;
  double error, prob;

  calculator.integral = 0.0;

  Integrate::CuhreIntegrate( ndim, Integrand, userdata, epsrel, epsabs,
                             Integrate::LastImportanceFcnt,
                             mineval, maxeval, nregions, neval,
                             fail, calculator.integral, error, prob );
  
  //should check for fails and stuff all over here
  
  
/*
  static std::mutex m;
  std::lock_guard<std::mutex> scoped_lock( m );
  cerr << "After " << neval << " evaluations got integral " << calculator.integral
       << endl << "\t calculator.energy=" << calculator.energy
       << " calculator.srcVolumumetricActivity=" << calculator.srcVolumumetricActivity
       << ", rad=" << calculator.m_sphereRadAndTransLenCoef[0].first
       << ", transCoef=" << calculator.m_sphereRadAndTransLenCoef[0].second
       << endl;
*/
}//void selfShieldingIntegration(...)

  
void PointSourceShieldingChi2Fcn::setBackgroundPeaks(
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
          << " will not be used for backgorund peak area subtraction;"
          << " non-gaussian peaks may be supported for this in the future";
      passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    }//if( p.gausPeak() ) / else
  }//for( PeakDef &p : m_backgroundPeaks )
}//void setBackgroundPeaks(...)
  
  
vector< tuple<double,double,double,Wt::WColor,double> >
       PointSourceShieldingChi2Fcn::energy_chi_contributions(
           const std::vector<double> &x,
           PointSourceShieldingChi2Fcn::NucMixtureCache &mixturecache,
                                        const bool allow_multiple_iso_contri,
                                        std::vector<std::string> *info                     
                                                             ) const
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
    if( allow_multiple_iso_contri )
      info->push_back( "Allowing multiple nuclides being fit for to potentially contribute to the same photopeak" );
    else
      info->push_back( "Not allowing multiple nuclides being fit for to contribute to the same photopeak" );
    
    //SHould put in information about the shielding here
  }//if( info )
  
  
  EnergyCountMap energy_count_map;
  const vector<pair<double,double> > energie_widths = observedPeakEnergyWidths( m_peaks );
  
  if( allow_multiple_iso_contri )
  {
    //Get the number of source gammas from each nuclide
    for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
    {
      if( isActivityDefinedFromShielding( nuclide ) )
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
                               m_photopeakClusterSigma, -1.0, info );
    }//for( const SandiaDecay::Nuclide *nuclide : m_nuclides )
  }else
  {
    for( const PeakDef &peak : m_peaks )
    {
      const SandiaDecay::Nuclide *nuclide = peak.parentNuclide();
      const SandiaDecay::RadParticle *particle = peak.decayParticle();
      
      if( !nuclide
          || (!particle && (peak.sourceGammaType()!=PeakDef::AnnihilationGamma))
          || isActivityDefinedFromShielding(nuclide) )
        continue;
      
      const double act = activity( nuclide, x );
      const double thisage = age( nuclide, x );
      
      if( mixturecache.find(nuclide) == mixturecache.end() )
        mixturecache[nuclide].addNuclideByActivity( nuclide, sm_activityUnits );
      
      const float energy = peak.gammaParticleEnergy();
      cluster_peak_activities( energy_count_map, energie_widths,
                               mixturecache[nuclide], act, thisage,
                               m_photopeakClusterSigma, energy, info);
    }//for( const PeakDef &peak : m_peaks )
  }//if( allow_multiple_iso_contri )

  //Propogate the gammas through each material - note we are using the fit peak
  //  mean here, and not the (pre-cluster) photopeak energy
  const int nMaterials = static_cast<int>( m_materials.size() );
  for( int materialN = 0; materialN < nMaterials; ++materialN )
  {
    boost::function<double(float)> att_coef_fcn;
    const Material * const material = m_materials[materialN].first;
//    const vector<const SandiaDecay::Nuclide *> &srcs
//                                                = m_materials[materialN].second;

    if( isGenericMaterial( materialN ) )
    {
      const float atomic_number = static_cast<float>(atomicNumber( materialN, x ));
      const float areal_density = static_cast<float>(arealDensity( materialN, x ));
      
      att_coef_fcn = boost::bind( &transmition_coefficient_generic,
                                  atomic_number, areal_density, _1 );
    }else// if( srcs.empty() )  //XXX - I cant figure out why this if statment is here
    {
      float thick = static_cast<float>( thickness( materialN, x ) );

      //I think this next step isnt even needed, but just to make sure...
      if( SpecUtils::iequals_ascii( material->name, "void") )
        thick = 0.0;

      att_coef_fcn = boost::bind( &transmition_coefficient_material,
                                  material, _1, thick );
    }//if( isGenericMaterial( materialN ) ) / else

/*
    //Uncomment out this section if you want to print out the transmition
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

    if( att_coef_fcn )
    {
      if( info )
      {
        if( isGenericMaterial(materialN) )
        {
          const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
          const double ad = arealDensity(materialN, x) / adUnits;
          const double an = atomicNumber(materialN, x);
          
          char buffer[256];
          snprintf( buffer, sizeof(buffer),
                    "Shielding %i: AN=%.2f, AD=%.2f g/cm2", materialN, an, ad );
          info->push_back( buffer );
        }else
        {
          const double thick = thickness( materialN, x );
          stringstream title;
          title << "Shielding: " << material->name << ", "
                << material->chemicalFormula() << ", density="
                << material->density *PhysicalUnits::cm3/PhysicalUnits::gram
                << " g/cm3, thickness=" << (thick/PhysicalUnits::cm) << " cm";
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
        energy_count.second *= exp( -1.0 * att_coef_fcn( energy_count.first ) );
    }//if( att_coef_fcn )
  }//for( int materialN = 0; materialN < nMaterials; ++materialN )


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
      
      if( info  )
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
      info->push_back( "Solid angle reduces counts by a factor of "
                      + std::to_string(fracAngle) );
  }//if( m_detector && m_detector->isValid() ) / else


  //This is where contributions from sources defined from the shielding should
  //  be calculated
#if( ALLOW_MASS_FRACTION_FIT )
  vector<MaterialAndSources> materials = m_materials;
  vector<std::shared_ptr<Material> > customMaterials;
#else
  const std::vector<MaterialAndSources> &materials = m_materials;
#endif
  
  if( info )
  {
    for( int materialN = 0; materialN < nMaterials; ++materialN )
    {
      const vector<const SandiaDecay::Nuclide *> &srcs
                                         = materials[materialN].second;
      if( !isGenericMaterial(materialN) && !srcs.empty() )
      {
        info->push_back( "Self Attenuating Source Info (shielding, detector,"
                        " distance, live time, and amount of material not-accounted for):" );
        break;
      }
    }//for( int materialN = 0; materialN < nMaterials; ++materialN )
  }//if( info )

  using GammaInteractionCalc::transmition_length_coefficient;
  
  vector<SelfAttCalc> calculators;
  for( int materialN = 0; materialN < nMaterials; ++materialN )
  {
    const Material *material = materials[materialN].first;
    const vector<const SandiaDecay::Nuclide *> &srcs
                                                = materials[materialN].second;

    if( isGenericMaterial( materialN ) || srcs.empty() )
      continue;

#if( ALLOW_MASS_FRACTION_FIT )
    if( hasVariableMassFraction( material ) )
    {
      std::shared_ptr<Material> mat = variedMassFracMaterial( material, x );
      customMaterials.push_back( mat );
      materials[materialN].first = material = mat.get();
    }//if( hasVariableMassFraction( material ) )
#endif
    
    SelfAttCalc baseCalculator;
    if( m_detector )
      baseCalculator.m_detectorRadius = 0.5 * m_detector->detectorDiameter();
    else
      baseCalculator.m_detectorRadius = 0.5 * PhysicalUnits::cm;

    baseCalculator.m_observationDist = m_distance;
    baseCalculator.m_sourceIndex = materialN;

    
    for( const SandiaDecay::Nuclide *src : srcs )
    {
      EnergyCountMap local_energy_count_map;
      const double actPerMass = src->activityPerGram() / PhysicalUnits::gram;
      double massFract = 0.0;
      
      for( const Material::NuclideFractionPair &nfp : material->nuclides )
      {
        if( nfp.first == src )
          massFract += nfp.second;
      }//for( nuclide in this material )
      
      const double actPerVol = actPerMass * massFract * material->density;
      const double thisage = age( src, x );
      
      //      cerr << "For " << src->symbol << " actPerVol=" << actPerVol
      //           << ", material->density=" << material->density << "("
      //           << material->density *PhysicalUnits::cm3/PhysicalUnits::g
      //           << "), actPerMass=" << actPerMass << ", massFraction="
      //           << massFraction << endl;
      
      if( mixturecache.find(src) == mixturecache.end() )
        mixturecache[src].addNuclideByActivity( src, sm_activityUnits );
      
      if( allow_multiple_iso_contri )
      {
        cluster_peak_activities( local_energy_count_map, energie_widths,
                                 mixturecache[src], actPerVol, thisage,
                                 m_photopeakClusterSigma, -1.0, info );
      }else
      {
        for( const PeakDef &peak : m_peaks )
        {
          if( peak.parentNuclide()==src
             && (peak.decayParticle() || (peak.sourceGammaType()==PeakDef::AnnihilationGamma)) )
            cluster_peak_activities( local_energy_count_map, energie_widths,
                                     mixturecache[src], actPerVol, thisage,
                                      m_photopeakClusterSigma,
                                     peak.gammaParticleEnergy(), info );
        }//for( const PeakDef &peak : m_peaks )
      }//if( allow_multiple_iso_contri ) / else

      for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
      {
        SelfAttCalc calculator = baseCalculator;

        double radius = 0.0;
        calculator.nuclide = src;
        calculator.energy = energy_count.first;
        calculator.srcVolumumetricActivity = energy_count.second;

        for( int subMat = 0; subMat < nMaterials; ++subMat )
        {
          if( isGenericMaterial( subMat ) )
            continue;

          const Material *const material = materials[subMat].first;
          radius += thickness( subMat, x );
          if( radius > m_distance )
            throw runtime_error( "PointSourceShieldingChi2Fcn::"
                                "energy_chi_contributions: radius > distance" );

          
          const double transLenCoef = transmition_length_coefficient( material, calculator.energy );

          calculator.m_sphereRadAndTransLenCoef.push_back( make_pair(radius,transLenCoef) );
        }//for( int subMat = 0; subMat < nMaterials; ++subMat )

        if( calculator.m_sphereRadAndTransLenCoef.empty() )
          throw std::logic_error( "No source/shielding sphere for calculator" );
        
        calculators.push_back( calculator );
      }//for( const EnergyCountMap::value_type &energy_count : local_energy_count_map )
    }//for( const SandiaDecay::Nuclide *src : srcs )
  }//for( int materialN = 0; materialN < nMaterials; ++materialN )

  if( calculators.size() )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( SelfAttCalc &calculator : calculators )
      pool.post( boost::bind( &PointSourceShieldingChi2Fcn::selfShieldingIntegration, boost::ref(calculator) ) );
    pool.join();
    
//    vector<boost::function<void()> > workers;
//    for( SelfAttCalc &calculator : calculators )
//    {
//      boost::function<void()> worker = boost::bind( &PointSourceShieldingChi2Fcn::selfShieldingIntegration, boost::ref(calculator) );
//      workers.push_back( worker );
//    }//for( size_t i = 0; i < calculators.size(); ++i )
//    SpecUtils::do_asyncronous_work( workers, false );

    if( info )
      info->push_back( "Self Attenuating Source Info (after accounting for "
                       "shielding, distanec, detector and live time):" );
    
    for( const SelfAttCalc &calculator : calculators )
    {
      double contrib = calculator.integral * calculator.srcVolumumetricActivity;

      if( m_detector && m_detector->isValid() )
        contrib *= m_detector->intrinsicEfficiency( calculator.energy );

      if( energy_count_map.find( calculator.energy ) != energy_count_map.end() )
      {
//        cerr << "Adding " << contrib << " to energy " << calculator.energy
//             << ", calculator.integral=" << calculator.integral
//             << ", calculator.srcVolumumetricActivity=" << calculator.srcVolumumetricActivity
//             << endl;
        energy_count_map[calculator.energy] += contrib;
      }else
      {
//        cerr << "Setting " << contrib*m_liveTime << " counts to energy " << calculator.energy
//             << " for thickness=" << calculator.m_sphereRadAndTransLenCoef[calculator.m_sourceIndex].first / PhysicalUnits::cm
//             << " cm" << endl;
        energy_count_map[calculator.energy] = contrib;
      }
      
      if( info )
      {
        const Material *const material = materials[calculator.m_sourceIndex].first;
        const int index = static_cast<int>( calculator.m_sourceIndex );
        double thick = calculator.m_sphereRadAndTransLenCoef[index].first;
        if( index > 0 )
          thick -= calculator.m_sphereRadAndTransLenCoef[index-1].first;
        
        stringstream msg;
        msg << "\tAttributing " << contrib*PhysicalUnits::second << " cps to "
            << calculator.energy/PhysicalUnits::keV << " keV photopeak ";
        if( calculator.nuclide )
          msg << "(from " << calculator.nuclide->symbol << ") "; 
        msg << "for thickness " << PhysicalUnits::printToBestLengthUnits(thick)
            << " " << material->name;
        info->push_back( msg.str() );
      }//if( info )
    }//for( SelfAttCalc &calculator : calculators )
  }//if( calculators.size() )

  
  //Account for live time
  for( EnergyCountMap::value_type &energy_count : energy_count_map )
    energy_count.second *= m_liveTime;

  return expected_observed_chis( m_peaks, m_backgroundPeaks, energy_count_map, info );
}//vector<tuple<double,double,double> > energy_chi_contributions(...) const

}//namespace GammaInteractionCalc
