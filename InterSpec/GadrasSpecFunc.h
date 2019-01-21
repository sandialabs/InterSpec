#ifndef GadrasSpecFunc_h
#define GadrasSpecFunc_h
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

#include <string>
#include <vector>


/** \brief Computes the scattered continuum for a given energy gamma, through
 * shielding of a given areal density and atomic number.
 *
 * Adapted from SpecFunc.f (sent to wcjohns 20150922) code.
 * Uses Continuum.lib that is included with GadrasDRF v18.6.1.
 */
class GadrasScatterTable
{
public:
  /** Constructor reads in Continuum.lib from file.
   *
   * \param datafile Path and filename to Continuum.lib file that comes with 
   * GADRAS.
   *
   * Throws std::runtime_error() if file cant be opened, or is incorrectly 
   * sized.
   */
  GadrasScatterTable( const std::string &datafile );
  
  /** Destructor */
  ~GadrasScatterTable();
  
  /** Calculates scattered continuum for given energy gamma, through specified
   * shielding.
   *
   * \param answer Where results are placed. Each element cooresponds in energy
   *               to #binning.  Will be resized and zeroed out before computing
   *               answer.  Results will be normalized to coorespond to the 
   *               specified source intensity and shielding into 4 pi.
   * \param sourceEnergy Energy of incomming gamma, in keV.
   * \param sourceIntensity Number of gammas emmited from the source (so 
   *        pre-shielding) you want continuum to coorispond to.
   * \param atomicNumber Atomic number of shielding. H=1, Al=16, etc.
   * \param arealDensity Areal density of shielding, in g/cm2.  
   * \param fractionHydrogen Fractional areal density of hydrogen in shield
   * \param xstool attenuation tool
   * \param binning The binning structure you would like #answer placed in.
   *
   * \returns Number of un-collided gammas at source energy through the
   *          shielding specified.
   *
   * May throw exception if #binning is not monotonically increasing in energy.
   * Will throw exception if energy, areal density, or atomic number are out
   * of range.
   */
  float getContinuum( std::vector<float> &answer,
                     const float sourceEnergy,
                     const float sourceIntensity,
                     const float atomicNumber,
                     const float arealDensity,
                     const float fractionAdHydrogen,
                     const std::vector<float> &binning ) const;
  
protected:
  static const int sm_num_areal_density = 9;
  static const int sm_num_atomic_number = 6;
  static const int sm_num_input_energy = 16;
  static const int sm_num_out_energy = 32;
  
  float m_data[9][6][16][32];  //[areal density][atomic number][incomming energy][fractional output continuum energy]
};//class GadrasScatterTable


#endif  //GadrasSpecFunc_h
