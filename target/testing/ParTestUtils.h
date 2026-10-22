#ifndef ParTestUtils_h
#define ParTestUtils_h
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
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>

#include "io/DetectorResponse.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/AngleOutxImport.h"
#include "InterSpec/DetectorEffG2kPar.h"

/** Test-only helpers shared by the tests that need `.par` bytes: the layout
 `DetEffG2kPar::parseParFile` decodes, written back out.  InterSpec itself never
 writes this format.
 */
namespace ParTestUtils
{
  inline void append_u16( std::string &out, const uint16_t v )
  {
    out.push_back( static_cast<char>( v & 0xFF ) );
    out.push_back( static_cast<char>( (v >> 8) & 0xFF ) );
  }

  template<class T>
  void append_le( std::string &out, const T value )
  {
    static_assert( (sizeof(T) == 4) || (sizeof(T) == 8), "float or double only" );
    uint8_t bytes[sizeof(T)];
    memcpy( bytes, &value, sizeof(T) );  // all supported platforms are little-endian
    out.append( reinterpret_cast<const char *>( bytes ), sizeof(T) );
  }

  /** The bytes of `par`: an energy header (emin, emax, count, energies), then one record per energy
   - a 28-byte record header (with the marker double at offset 16) followed by the uint16 cells.
   The file stores the angular and radial steps as floats, so a `ParFile` that is to compare equal
   to its own parse should hold float-representable steps.
   */
  inline std::string par_file_bytes( const DetEffG2kPar::ParFile &par )
  {
    if( par.energies_keV.size() != par.grids.size() )
      throw std::runtime_error( "par_file_bytes: energies and grids differ in count." );

    std::string out;
    append_le<double>( out, par.emin_keV );
    append_le<double>( out, par.emax_keV );
    append_u16( out, static_cast<uint16_t>( par.energies_keV.size() ) );
    for( const double energy : par.energies_keV )
      append_le<double>( out, energy );

    for( const DetEffG2kPar::ParGrid &g : par.grids )
    {
      if( g.V.size() != (static_cast<size_t>(g.nrows) * g.ncols) )
        throw std::runtime_error( "par_file_bytes: grid cell count does not match its dimensions." );

      append_u16( out, g.ncols );
      append_u16( out, g.nrows );
      out.append( 8, '\0' );
      append_le<float>( out, static_cast<float>( g.theta_step_rad ) );
      append_le<double>( out, 4707532.0 );
      append_le<float>( out, static_cast<float>( g.r_step ) );
      for( const uint16_t v : g.V )
        append_u16( out, v );
    }//for( each energy's grid )

    return out;
  }//par_file_bytes(...)


  /** The geometry test_data/det_eff/mock.par was generated from (see WriteMockParFiles in
   test_DetectorEffG2kPar.cpp): what the importer makes of mock_DETECTOR.txt (`def`), with the core
   and bulletizing of the ANGLE detector at `detx_path` (Angle-detector-only.detx) in place of the
   importer's guesses - the two things a DETECTOR.txt cannot state.
   */
  inline ceelo::GeometryDescriptor mock_truth_geometry( const DetEffG2kPar::DetectorDef &def,
                                                        const std::string &detx_path )
  {
    std::ifstream detx( detx_path.c_str(), std::ios::in | std::ios::binary );
    if( !detx.is_open() )
      throw std::runtime_error( "mock_truth_geometry: could not open " + detx_path );

    std::vector<std::string> warnings;
    const ceelo::GeometryDescriptor angle = CeeLoUtils::buildAngleGeometry( AngleOutx::parse( detx ), warnings );

    ceelo::GeometryDescriptor gd = DetEffG2kPar::geometryFromDetectorDef( def, warnings );
    gd.bore = angle.bore;
    gd.bullet_radius_cm = angle.bullet_radius_cm;
    if( !gd.bore )
      throw std::runtime_error( "mock_truth_geometry: the ANGLE detector has no core." );
    return gd;
  }//mock_truth_geometry(...)
}//namespace ParTestUtils

#endif //ParTestUtils_h
