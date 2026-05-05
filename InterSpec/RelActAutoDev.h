#ifndef RelActAutoDev_h
#define RelActAutoDev_h
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
#include <cstddef>

#include "InterSpec_config.h"



namespace RelActAutoDev
{
  int dev_code();

  enum class IdbMaterialType { Uranium, Plutonium, MOX };

  /** Runs RelActCalcAuto::solve on each HPGe spectrum in idb_dir, comparing results
   against the known mass fractions in metadata_xml_path. Writes a single HTML
   summary to output_html_path. Processes at most max_spectra HPGe spectra
   (pass std::numeric_limits<size_t>::max() to process all). */
  void run_idb_comparison( const std::string &idb_dir,
                           const std::string &metadata_xml_path,
                           IdbMaterialType material_type,
                           const std::string &output_html_path,
                           size_t max_spectra = 15 );
}//namespace RelActAutoDev


#endif //RelActAutoDev_h
