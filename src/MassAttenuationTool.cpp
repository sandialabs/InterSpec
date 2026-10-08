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

#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>

#include "cross_sections/CrossSectionData.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/MassAttenuationTool_imp.hpp"

using namespace std;

static_assert( MassAttenuation::sm_max_xs_atomic_number == ceelo::kMaxZ,
               "MassAttenuation must cover the same elements as CeeLo" );
static_assert( (MassAttenuation::sm_min_xs_energy_keV == ceelo::kPhotonDataMinEnergy_keV)
               && (MassAttenuation::sm_max_xs_energy_keV == ceelo::kPhotonDataMaxEnergy_keV),
               "MassAttenuation must advertise CeeLo's photon-data energy range" );

namespace
{
  /** Checks the inputs, and returns the energy, in MeV, clamped into CeeLo's photon-data range. */
  double checked_energy_mev( const int atomic_number, const double energy )
  {
    if( (atomic_number < MassAttenuation::sm_min_xs_atomic_number)
       || (atomic_number > MassAttenuation::sm_max_xs_atomic_number) )
      throw runtime_error( "Invalid atomic number (" + std::to_string(atomic_number) + ") for mass attenuation" );

    if( !std::isfinite(energy) )
      throw runtime_error( "Invalid energy for mass attenuation" );

    const double energy_kev = std::clamp( energy / PhysicalUnits::keV,
                                          MassAttenuation::sm_min_xs_energy_keV, MassAttenuation::sm_max_xs_energy_keV );
    return energy_kev / 1000.0;
  }//checked_energy_mev(...)


  /** Multiplier from barns/atom to a mass attenuation coefficient in PhysicalUnits, for the given (already validated)
   atomic number.  Uses CeeLo's atomic weights, so the result matches CeeLo's own transport.
   */
  double barns_to_mass_atten( const int atomic_number )
  {
    static const std::array<double,MassAttenuation::sm_max_xs_atomic_number + 1> s_factors = [](){
      const double avogadro = 6.02214076E23;
      const double barn = 1.0E-24 * PhysicalUnits::cm2;
      const ceelo::CrossSectionData &xs = ceelo::CrossSectionData::instance();
      std::array<double,MassAttenuation::sm_max_xs_atomic_number + 1> factors{};
      for( int z = MassAttenuation::sm_min_xs_atomic_number; z <= MassAttenuation::sm_max_xs_atomic_number; ++z )
        factors[z] = barn * avogadro / (xs.atomic_weight(z) * PhysicalUnits::gram);
      return factors;
    }();

    return s_factors[atomic_number];
  }//barns_to_mass_atten(...)
}//namespace


namespace MassAttenuation
{
  double massAttenuationCoefficientElement( const int atomic_number, const double energy )
  {
    const double energy_mev = checked_energy_mev( atomic_number, energy );
    const ceelo::CrossSectionData &xs = ceelo::CrossSectionData::instance();

    // sigma_pair_production(...) is zero below threshold
    const double barns = xs.sigma_compton( atomic_number, energy_mev )
                         + xs.sigma_photoelectric( atomic_number, energy_mev )
                         + xs.sigma_pair_production( atomic_number, energy_mev );

    return barns * barns_to_mass_atten( atomic_number );
  }//double massAttenuationCoefficientElement( const int atomic_number, const double energy )


  double massAttenuationCoefficientElement( const int atomic_number, const double energy, GammaEmProcces process )
  {
    const double energy_mev = checked_energy_mev( atomic_number, energy );
    const ceelo::CrossSectionData &xs = ceelo::CrossSectionData::instance();

    double barns = 0.0;
    switch( process )
    {
      case GammaEmProcces::ComptonScatter:  barns = xs.sigma_compton( atomic_number, energy_mev );         break;
      case GammaEmProcces::RayleighScatter: barns = xs.sigma_rayleigh( atomic_number, energy_mev );        break;
      case GammaEmProcces::PhotoElectric:   barns = xs.sigma_photoelectric( atomic_number, energy_mev );   break;
      case GammaEmProcces::PairProduction:  barns = xs.sigma_pair_production( atomic_number, energy_mev ); break;
      case GammaEmProcces::NumGammaEmProcces:
      default:
        throw runtime_error( "Invalid EM process" );
    }//switch( process )

    return barns * barns_to_mass_atten( atomic_number );
  }//double massAttenuationCoefficientElement( atomic_number, energy, process )


  double massAttenuationCoefficientFracAN( const double atomic_number, const double energy )
  {
    return mass_atten_coef_frac_an<double>( atomic_number, energy );
  }
}//namespace MassAttenuation
