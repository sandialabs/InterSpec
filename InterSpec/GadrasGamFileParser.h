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
 * Opens only the newer format of GAM files (that start with two 0's).
 */
class GadrasGamFile
{
public:
  GadrasGamFile();
  
  
  /** Parses the cooresponding stream.
   * Throws runtime_error exception when the file is not formatted properly.
   */
  void parse_data( std::istream &strm );
  
  
  
  /* A relevant note from the GADRAS users manual:
   * If the width of the group exceeds 1% of the average energy of the group,
   * the radiation within the group is assumed to be distributed uniformly 
   * within the group.  For groups that have widths under 1% of the average 
   * energy, any excess flux within the group relative to a linear interpolation 
   * of surrounding groups is assumed to be emitted at a discrete energy, 
   * corresponding to the midpoint of the group.
   */
  
  /** Energy, in keV of the un-collided gamma rays. 
   * Will have same number of entries as, and coorespond index-wise to
   * #m_photon_lines_flux, #m_photon_lines_an, and #m_photon_lines_ad.
   */
  std::vector<float> m_photon_lines_energy;
  
  /** The flux into 4pi for the cooresponding m_photon_lines_energy energy.*/
  std::vector<float> m_photon_lines_flux;
  
  /** The yeild-effective weighted atomic number of interveining material. */
  std::vector<float> m_photon_lines_an;
  
  /** The yield-effective weighted areal density of interveining material. */
  std::vector<float> m_photon_lines_ad;
  
  
  /** Lower energy of energy group for the cooresponding (index-wise) 
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
  
  /** Flux in cooresponding energy group.
   */
  std::vector<float> m_neutron_group_flux;
};//class GadrasGamFile

#endif //GadrasGamFileParser_h




