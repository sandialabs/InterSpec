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

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <algorithm>

#include <boost/math/tools/roots.hpp>

#include <Wt/WString.h>

#include "InterSpec/DoseCalc.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/GadrasShieldScatter.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"


using namespace std;


namespace DoseCalc
{
  /** Applies coefficients in ANSI/ANS-6.1.1-1991 to calculate fluence-to-dose 
   * factor (Sv-cm2).
   *
   * \param energy Energy in MeV
   *
   */
  double dose_equiv( const float energy, const float *C )
  {
    const double X = std::log(energy);
    return 1.0E-12*exp(C[0]+X*(C[1]+X*(C[2]+X*(C[3]+X*C[4]))));  //equiv to 1.0E-12 * exp(C[0] + C[1]*X + C[2]*X*X + C[3]*X*X*X + C[4]*X*X*X*X);
  }
  
  float gamma_dose( const float energy_kev )
  {
    double dose = 0.0;
    
    const float below_cut[5] = { 2.30402f, 0.75167f, -1.04725f, -0.50091f, -0.06302f };
    const float above_cut[5] = { 1.52070f, 0.79329f, -0.07261f,  0.01228f,  0.00347f };
    
    const float energy = energy_kev / 1000.0f;

    if( energy <= 3.0 )
    {
      const float *coefs = (energy <= 0.15f) ? static_cast<const float *>(below_cut) : static_cast<const float *>(above_cut);
      dose = dose_equiv( energy, coefs );
    }else if( energy <= 6.9 )
    {
      dose = 8.760E-13 + energy*(3.988E-12 - energy*2.895E-13);
    }else
    {
      dose = 1.461E-11;
    }
    
    return static_cast<float>( dose * PhysicalUnits::cm2 * PhysicalUnits::sievert );
  }//float gamma_dose( const float energy )
  
  
  float neutron_dose( const float energy_kev )
  {
    float dose = 0.0f;
    const float below_coefs[5] = { 3.430895f, 0.7725710f,  9.834081E-2f, 4.903466E-3f, 8.149667E-5f };
    const float above_coefs[5] = { 4.952167f, 0.6644235f, -0.1017445f,  -1.496004E-3f, 3.636748E-3f };
    
    const float energy = energy_kev / 1.0E3f; //convert from keV to MeV
    
    if( energy <= 17.7)          //dose only defined to 17.7 MeV
    {
      const float *coefs = ((energy < 0.01f) ? below_coefs : above_coefs);
      dose = dose_equiv(energy,coefs);
    }else                        //extrapolate beyond 17.7 MeV
    {
      dose = 5.0966831E-10f;  //Dose per neutron
    }
    
    return dose * PhysicalUnits::cm2 * PhysicalUnits::sievert;
  }//float neutron_dose( const float energy )

  bool less_than_first( const std::pair<float,float> &lhs, const float rhs )
  {
    return lhs.first < rhs;
  }

   
  
  double gamma_dose_with_shielding( const vector<float> &energies,
                                    const vector<float> &intensity,
                                    const float areal_density,
                                    const float atomic_number,
                                    const float distance,
                                    const GadrasShieldScatter &scatter )
  {
    if( energies.empty() || (intensity.size() != energies.size()) )
      throw runtime_error( "gamma_dose_with_shielding(): No gamma lines to calculate dose from." );
    
    const float maxenergy = *std::max_element( energies.begin(), energies.end() );
    
    vector<float> energy_groups( 1024 );
    vector<double>  continuum( 1024, 0.0 );  //double to avoid numerical issues with repeated adding
    for( size_t i = 0; i < energy_groups.size(); ++i )
      energy_groups[i] = (maxenergy*i) / (energy_groups.size() - 1);
    
    if( distance <= 0.0 )
      throw runtime_error( "gamma_dose_with_shielding(): Distance must be greater than zero." );
    
    // The upper areal-density bound comes from the scatter database itself
    // (sandia.shieldscatter.db). Values above this would be extrapolated in
    // log-scaled scatter space, so we reject them rather than returning a
    // bogus dose.
    const double max_ad = scatter.maxArealDensity();
    if( atomic_number < 1 || atomic_number > 98
        || areal_density < 0.0
        || (areal_density*PhysicalUnits::cm2/PhysicalUnits::g) > max_ad )
      throw runtime_error( "gamma_dose_with_shielding(): Invalid shielding;"
                           " atomic number must be between 1 and 98, and areal"
                           " density between zero and " + std::to_string(max_ad)
                           + " g/cm2." );
    
    //Calculate atomic number and areal density, taking into account how the air
    //  effects the shielding.
    float an = atomic_number;
    float ad = areal_density;
    
    //TODO: not taking into account
    // At 3m, 59 keV gammas are attenuated by 6.5%
    if( distance > 300.0f*PhysicalUnits::cm )
    {
      const float air_an = 7.3737f;  //Gadras uses 7.2
      const float air_density = static_cast<float>( 0.00129 * PhysicalUnits::g / PhysicalUnits::cm3 );
      const float air_ad = air_density * distance; //0.00129*A(61)*DISTANCE
      an = (atomic_number*areal_density + air_an*air_ad) / (areal_density + air_ad);
      ad = areal_density + air_ad;
    }//if( distance > 300.0f*PhysicalUnits::cm )
    
    const double areal_density_g_cm2 = ad*PhysicalUnits::cm2/PhysicalUnits::g;
    
    double dose_nonscatter = 0.0f;
    vector< pair<float,float> > uncollided_lines;
    
    for( size_t energyInd = 0; energyInd < energies.size(); ++energyInd )
    {
      const float energy = energies[energyInd];
      const float I_0 = intensity[energyInd];
//      const double mu = MassAttenuation::massAttenuationCoefficientFracAN( an, energy );
//      const double transmittion = exp( -mu * areal_density );
      
      const float hydrogen_frac_ad = 0.0f;
      vector<float> thiscontinuum;
      const float uncollided
             = scatter.getContinuum( thiscontinuum, energy, I_0,
                                     an, areal_density_g_cm2, hydrogen_frac_ad,
                                     energy_groups );
      
      uncollided_lines.push_back( make_pair(energy, uncollided) ) ;
      
      assert( thiscontinuum.size() == energy_groups.size() );
      assert( continuum.size() == thiscontinuum.size() );
      
      //Add in the dose from the uncoloided gamma, pre geom factor
      dose_nonscatter += uncollided * gamma_dose( energy );
      
      if( ad > 0 )
      {
        for( size_t i = 0; i < thiscontinuum.size(); ++i )
          continuum[i] += thiscontinuum[i];
      }
    }//for( size_t energyGroup = 0; energyGroup < energies.size(); ++energyGroup )
    
    if( false )
    {//begin write a debug output...
      std::cerr << "\nWrtingin flux_out.csv file!\n";
      std::ofstream output( "flux_out.csv", ios::binary | ios::out );
      output << "energy,counts" << endl;
      vector<double> specoutput = continuum;
      
//      for( size_t i = 0; i < energies.size(); ++i )
      for( size_t i = 0; i < uncollided_lines.size(); ++i )
      {
//        const double mu = MassAttenuation::massAttenuationCoeficient( (int)atomic_number, energies[i] );
        const double uncollided_intensity = uncollided_lines[i].second; //intensity[i] * exp( -mu * areal_density );
//        const size_t pos = std::lower_bound( energy_groups.begin(), energy_groups.end(), energies[i] ) - energy_groups.begin();
        const size_t pos = std::lower_bound( energy_groups.begin(), energy_groups.end(), uncollided_lines[i].first ) - energy_groups.begin();
        specoutput[pos] += uncollided_intensity;
      }//for( size_t i = 0; i < energies.size(); ++i )
      
      for( size_t i = 0; i < energy_groups.size(); ++i )
        output << energy_groups[i] << "," << specoutput[i] << endl;
    }//end write a debug output...
    
    //Compute the dose from the continuum.
    double dose_continuum = 0.0;
    const size_t nbin = energy_groups.size();
    for( size_t i = 0; i < (nbin-1); ++i )
      dose_continuum += continuum[i]*gamma_dose( 0.5*(energy_groups[i] + energy_groups[i+1]) );
    
    //Get the last energy froup in the continuum
    dose_continuum += continuum[nbin-1]*gamma_dose( energy_groups[nbin-1] + 0.5*(energy_groups[nbin-1] - energy_groups[nbin-2]) );
    
    const double surfacearea = (4.0*PhysicalUnits::pi*distance*distance);
    
    return (dose_nonscatter + dose_continuum) / surfacearea;
  }//double gamma_dose_with_shielding(...)


  double fit_areal_density( const vector<float> &energies,
                            const vector<float> &intensities,
                            const float atomic_number,
                            const double user_entered_dose,
                            const double distance,
                            const GadrasShieldScatter &scatter )
  {
    // Dose is monotone-decreasing in areal density across all tested
    // (source, shield-Z) pairs — see DoseMonotonicityScan in
    // target/testing/test_DoseCalc.cpp. The root is bracketed by f(0) > 0
    // and f(max_ad) < 0; if either guard fails the inverse problem has no
    // physical solution and we throw.
    const double dose_no_shielding = gamma_dose_with_shielding(
        energies, intensities, 0.0f, atomic_number, distance, scatter );

    if( dose_no_shielding < user_entered_dose )
      throw runtime_error( "Dose from source specified is less than"
                          " entered dose, even with out shielding." );

    // A few g/cm^2 below the table top: scatter values right at the edge
    // are unreliable, and an answer pinned there is unphysical for the
    // user's purpose.
    const double max_ad_g_cm2 = scatter.maxArealDensity() - 5.0;
    const float max_ad_internal
      = static_cast<float>( max_ad_g_cm2 * PhysicalUnits::g / PhysicalUnits::cm2 );

    const double dose_at_max_ad = gamma_dose_with_shielding(
        energies, intensities, max_ad_internal, atomic_number, distance, scatter );

    if( dose_at_max_ad > user_entered_dose )
      throw runtime_error( "Over " + std::to_string(max_ad_g_cm2)
                          + " g/cm2 shielding is required - can not compute." );

    auto residual = [&]( const double ad_g_cm2 ) -> double {
      const float ad_internal
        = static_cast<float>( ad_g_cm2 * PhysicalUnits::g / PhysicalUnits::cm2 );
      const double dose = gamma_dose_with_shielding(
          energies, intensities, ad_internal, atomic_number, distance, scatter );
      return dose - user_entered_dose;
    };

    // Bracket the root to within 1e-3 g/cm^2; on a 251 g/cm^2 range that
    // is ~18 bisect iterations. Cap well above for pathological cases.
    const double ad_tol_g_cm2 = 1.0e-3;
    auto term_cond = [ad_tol_g_cm2]( double lo, double hi ) -> bool {
      return std::fabs( lo - hi ) < ad_tol_g_cm2;
    };

    boost::uintmax_t num_iters = 100;
    const boost::uintmax_t max_iters = num_iters;
    const std::pair<double, double> bracket
      = boost::math::tools::bisect( residual, 0.0, max_ad_g_cm2, term_cond, num_iters );

    if( num_iters >= max_iters )
      throw runtime_error( Wt::WString::tr("dcw-err-failed-fit-AD").toUTF8() );

    const double ad_g_cm2 = 0.5 * ( bracket.first + bracket.second );
    return ad_g_cm2 * PhysicalUnits::g / PhysicalUnits::cm2;
  }//double fit_areal_density(...)
}//namespace DoseCalc



