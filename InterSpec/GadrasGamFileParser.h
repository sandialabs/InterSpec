#ifndef GadrasGamFileParser_h
#define GadrasGamFileParser_h
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
#include <iostream>

/** Opens GADRAS GAM files.
 *
 * Currently only opens .gam files that start with two 0's, called version 1, and GamFileVersion = 2.0 and  2.1.
 *
 * Does not currently parse all information, such as model geometry, largest dimension, multiplicity, etc.
 */
class GadrasGamFile
{
public:
  GadrasGamFile();
  
  
  /** Parses the corresponding stream.
   * Throws runtime_error exception when the file is not formatted properly.
   */
  void parse_data( std::istream &strm );
  
  
  enum class GamFileVersion
  {
    /** The .gam files that start with two 0's.
     The gamma lines do not include shielding information.
     */
    Version_1,
    
    /** Indicated by a "GamFileVersion = 2.0" header in the file. */
    Version_2_0,
    
    /** Indicated by a "GamFileVersion = 2.1" header in the file. */
    Version_2_1,
    
    /**  */
    NumGamFileVersions
  };//enum class GamFileVersion
  
  
  /** Writes the .gam file to the provided stream.
   *
   * Currently only supports writing GamFileVersion::Version_2_0, and will throw exception for other versions.
   *
   * Throws runtime_error exception if output fails or errors encountered.
   */
  void write_gam( std::ostream &output, const GamFileVersion version );
  
  
  /*
   A relevant note is:
   If the width, W, or a photon group with lower energy EL, meets any of the following conditions,
   GADRAS transport will treat the group as a discrete line:
   - W < 0.25 keV
   - W < 0.4 keV and EL > 200 keV
   - W < 1.0 keV and EL > 10 MeV
   If any of these conditions are met, GADRAS will look at the (assumed) continuum groups above and
   below the potential discrete group. The counts-per-keV in the potential discrete group must be
   greater than the adjacent continuum groups. If it is, it uses the higher-energy continuum group
   and subtracts those counts from the discrete group to get the net discrete counts, with energy
   equal to the midpoint of the group.
   */
  
  /** The gamma file version parsed.
   */
  GamFileVersion m_version;
  
  /** Energy, in keV of the un-collided gamma rays.
   * Will have same number of entries as, and coorespond index-wise to
   * #m_photon_lines_flux, #m_photon_lines_an, and #m_photon_lines_ad.
   */
  std::vector<float> m_photon_lines_energy;
  
  /** The flux into 4pi for the corresponding m_photon_lines_energy energy.*/
  std::vector<float> m_photon_lines_flux;
  
  /** The yield-effective weighted atomic number of intervening material. */
  std::vector<float> m_photon_lines_an;
  
  /** The yield-effective weighted areal density of intervening material. */
  std::vector<float> m_photon_lines_ad;
  
  
  /** Lower energy of energy group for the corresponding (index-wise)
   * #m_photon_group_flux.  Energy is in keV.
   */
  std::vector<float> m_photon_group_boundries;
  
  /** Flux of the continuum for the given energy group.
   * Note that the last last entry should have a flux of 0, in order to allow
   * defining the upper energy of the last group.
   */
  std::vector<float> m_photon_group_flux;
  
  
  /** Lower energy of the neutron energy group.  Will have same number of
   * entries as #m_neutron_group_flux.
   * Energy is in keV (note that .gam files use eV).
   */
  std::vector<float> m_neutron_group_boundries;
  
  /** Flux in corresponding energy group.
   */
  std::vector<float> m_neutron_group_flux;
};//class GadrasGamFile

#endif //GadrasGamFileParser_h




