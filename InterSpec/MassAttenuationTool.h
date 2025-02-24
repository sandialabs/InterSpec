#ifndef MassAttenuationTool_h
#define MassAttenuationTool_h
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

#include <string>
#include <vector>
#include <atomic>
#include <stdexcept>

#define USE_SNL_GAMMA_ATTENUATION_VALUES 0

namespace MassAttenuation
{
  /** Sets the base directory that the 'em_xs_data' dir and GadrasContinuum.lib
      should be in.  Defaults to 'data' but incase you cant change current
      working directory to somewhere you can give a relative path to, use this
      function.  You can only call this function before getting any
      cross-sections, or else exception will be thrown.
   */
#ifdef _WIN32
  void set_data_directory( const std::wstring &dir );
#else
  void set_data_directory( const std::string &dir );
#endif
  
  /** Tests if the directory pointed to by 'dir' has a sub-directory
      'em_xs_data' with valid data files.
   
       Throws exception if directory is invalid.
   */
  /*
#ifdef _WIN32
  void test_data_directory_validity( const std::wstring &dir );
#else
  void test_data_directory_validity( const std::string &dir );
#endif
  */
  
  
  static const int sm_min_xs_atomic_number = 1;
  
#if( USE_SNL_GAMMA_ATTENUATION_VALUES )
  static const int sm_max_xs_atomic_number = 94;
#else
  static const int sm_max_xs_atomic_number = 98;
#endif

/** As of 20201124, the valid energy ranges of the various process, for each element is at least the following:
    - Compton has range 1 to 100000
    - Rayleigh has range 1.00331 to 35436.5
    - PhotoElectric has range 1 to 73282.4
    - PairProd has range 1025.12 to 100000
 */

  enum class GammaEmProcces : int
  {
    ComptonScatter,
    RayleighScatter,
    PhotoElectric,
    PairProduction,
    NumGammaEmProcces
  };//enum GammaEmProcces
  
  
  /** Gives the mass attenuation coefficient in units of PhysicalUnits, that
   * is to print out to familiar units you would divide by
   * (PhysicalUnits::cm2 / PhysicalUnits::g) to put the result into units
   * of cm2/g.
   * To calculate the
   *
   *  Will throw ErrorLoadingDataException if the XS data cannot be loaded.
   *  Will throw std::runtime_error if atomic_number or energy is less than 1.01 keV, or larger than 100 MeV.
   *
   * \param atomic_number Atomic number ranging from 1 to 98, inclusive
   * \param energy Energy (in keV), between 1.01 and 100,000.
   * \returns The mass attenuation coefficient for:
   *          compton + pair production + photo electric.
   */
  float massAttenuationCoefficientElement( const int atomic_number, const float energy );
  
  /** Similar to the other #massAttenuationCoefficientElement function, but instead
   * only for a specific sub-proccess.
   */
  float massAttenuationCoefficientElement( const int atomic_number, const float energy, MassAttenuation::GammaEmProcces process );
  
  /** Similar to #massAttenuationCoeficient, but interpolated the result between
   * floor(atomic_number) and ceil(atomic_number) linearly.
   */
  float massAttenuationCoefficientFracAN( const float atomic_number, const float energy );
  
  
  /** Compute the total attenuation coefficient using GADRASs CrossSection.lib.
   * Assumes "data/CrossSection.lib" (from GADRAS) exists, and upon first calling
   * of this function will read it in; if reading fails, will throw
   * std::runtime_exception.  Function is thread safe, but if called
   * simultaneously from multiple threads at first, each may read in
   * CrossSection.lib, reulsting in a slight innefficncy (this is so we can not
   * have the data in memory unless will use it, while still being thread safe).
   *
   * \param energy Gamma energy in keV.
   * \param atomic_number Atomic number of material, inclusively between 1 and 94.
   * \param scatter_mu The non-coherent scatter cross section in cm2/g
   * \param photoelectric_mu The photoelectric cross section in cm2/g
   * \param pair_prod_mu The pair production cross section in cm2/g
   *
   * \returns total cross-section in units cm2/g - not in PhysicalUnits.
   
   from:   Frank Biggs and Ruth Lighthill, "Analytical Approximations for X-ray
   Cross Sections II," SNL Rept.,
   SC-RR-71 0507 (1971).
   (note there is an update to this in 1988 that includes more fitting
   intervals, see
   http://www2.cose.isu.edu/~tforest/Classes/NucSim/Day9/sand87_0070.pdf)
   
   and:    Frank Biggs and Ruth Lighthill, "Analytical Approximations for Total
   Pair-Production Cross Sections,"
   SNL Rept. SC-RR-68-619 (1968).
   */
  float AttCoef( const float energy, const int atomic_number,
                float &scatter_mu, float &photoelectric_mu, float &pair_prod_mu );
  
  
  
  /** Exception thrown when there is an error loading the cross-section data. */
  class ErrorLoadingDataException
  : public std::runtime_error
  {
  public:
    ErrorLoadingDataException( const std::string &msg )
    : std::runtime_error( msg )
    {}
  };//class ErrorLoadingDataException
  

/** A debug function to printout cross-sections to a CSV for plotting. */
  //void print_xs_csv( std::ostream &output )
}//namespace MassAttenuation


#endif  //MassAttenuationTool_h
