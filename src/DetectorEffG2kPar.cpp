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
#include <array>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <Eigen/Core>

#include "io/Pchip.h"
#include "io/ResponseKernel.h"
#include "io/DetectorResponse.h"
#include "materials/Material.h"
#include "geometry/Geometry.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "SandiaDecay.h"

#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DetectorEffG2kPar.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;

namespace
{
/* ---------------------------------------------------------------------------
   Portable little-endian readers.  SpecUtils' read_binary_data does not
   byte-swap portably, so - as in SpecFile_aspect.cpp:105/111 - we assemble the
   value byte-by-byte and memcpy the integer bit-pattern into the float/double.
   ------------------------------------------------------------------------- */
uint16_t read_u16le( const uint8_t *p )
{
  return static_cast<uint16_t>( p[0] | (static_cast<uint16_t>(p[1]) << 8) );
}

int32_t read_i32le( const uint8_t *p )
{
  const uint32_t u = static_cast<uint32_t>(p[0])
                   | (static_cast<uint32_t>(p[1]) << 8)
                   | (static_cast<uint32_t>(p[2]) << 16)
                   | (static_cast<uint32_t>(p[3]) << 24);
  int32_t v;
  memcpy( &v, &u, sizeof(v) );
  return v;
}

float read_f32le( const uint8_t *p )
{
  const uint32_t u = static_cast<uint32_t>(p[0])
                   | (static_cast<uint32_t>(p[1]) << 8)
                   | (static_cast<uint32_t>(p[2]) << 16)
                   | (static_cast<uint32_t>(p[3]) << 24);
  float v;
  memcpy( &v, &u, sizeof(v) );
  return v;
}

double read_f64le( const uint8_t *p )
{
  uint64_t u = 0;
  for( int i = 0; i < 8; ++i )
    u |= (static_cast<uint64_t>(p[i]) << (8 * i));
  double v;
  memcpy( &v, &u, sizeof(v) );
  return v;
}

// The per-record marker double, and its offset within each record.  Both are
//  derived facts from the files (see the header comment), not vendor constants.
const double sm_record_marker = 4707532.0;
const size_t sm_marker_in_record = 16;
const size_t sm_dir_entry_stride = 20;
}//anonymous namespace


namespace DetEffG2kPar
{

//=============================================================================
//  DETECTOR.txt  (ASCII geometry record)   - a port of detector_txt_decode.py
//=============================================================================

double DetectorDef::d( int i ) const
{
  switch( i )
  {
    case 0: return d1_crystal_diam_mm;
    case 1: return d2_crystal_len_mm;
    case 2: return d3_window_diam_mm;
    case 3: return d4_endcap_od_mm;
    case 4: return d5_endcap_len_mm;
    case 5: return d6_front_gap_mm;
    case 6: return d7_side_gap_mm;
  }
  return std::numeric_limits<double>::quiet_NaN();
}//DetectorDef::d(int)


namespace
{
/** Parses one comma-separated token to a number, or returns NaN (mirrors the
 Python `_num`, which returns None). */
double parse_num( const std::string &tok )
{
  std::string t = tok;
  SpecUtils::trim( t );
  if( t.empty() )
    return std::numeric_limits<double>::quiet_NaN();
  try
  {
    size_t pos = 0;
    const double v = std::stod( t, &pos );
    // Reject trailing garbage so "12abc" is not silently accepted as 12.
    while( pos < t.size() && std::isspace( static_cast<unsigned char>(t[pos]) ) )
      ++pos;
    return (pos == t.size()) ? v : std::numeric_limits<double>::quiet_NaN();
  }catch( const std::exception & )
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
}//parse_num(...)


/** Best-effort serial number from a "# ..." comment line: the token following
 an "S/N", "SN" or "s/n" marker.  Returns "" when none is found. */
std::string serial_from_comment( const std::string &comment )
{
  const std::string lc = SpecUtils::to_lower_ascii_copy( comment );
  size_t pos = std::string::npos;
  size_t skip = 0;
  const size_t sn_slash = lc.find( "s/n" );
  const size_t sn_plain = lc.find( "sn" );
  if( sn_slash != std::string::npos )
  {
    pos = sn_slash;
    skip = 3;
  }else if( sn_plain != std::string::npos )
  {
    pos = sn_plain;
    skip = 2;
  }

  if( pos == std::string::npos )
    return "";

  size_t i = pos + skip;
  while( i < comment.size()
         && !std::isalnum( static_cast<unsigned char>(comment[i]) ) )
    ++i;
  std::string out;
  while( i < comment.size()
         && (std::isalnum( static_cast<unsigned char>(comment[i]) )
             || comment[i] == '-') )
  {
    out.push_back( comment[i] );
    ++i;
  }
  return out;
}//serial_from_comment(...)


/** Classifies a record from its D3 window and layer stack (port of the Python
 `classify_detector`). */
DetectorDef::Kind classify_detector( const DetectorDef &rec )
{
  const double d3 = rec.d3_window_diam_mm;

  const std::string slot0_mat = SpecUtils::to_lower_ascii_copy( rec.layers[0].material );
  const double slot0_th = rec.layers[0].thickness_mm;
  const std::string slot1_mat = SpecUtils::to_lower_ascii_copy( rec.layers[1].material );

  if( (d3 == 0.0) || std::isnan( d3 ) )
  {
    if( slot1_mat.empty() && slot0_mat == "ge" && slot0_th > 0.1 )
      return DetectorDef::Kind::NCoax;
    if( slot1_mat.empty() && slot0_mat == "ge" )
      return DetectorDef::Kind::NCoax;
    return slot1_mat.empty() ? DetectorDef::Kind::NCoax : DetectorDef::Kind::Generic;
  }

  // D3 > 0: thin-window detector.
  if( slot0_mat == "ge" && slot0_th > 0.0 && !slot1_mat.empty() )
  {
    const double back_al = rec.layers[8].thickness_mm;
    if( back_al >= 4.0 )
      return DetectorDef::Kind::Falcon;
    return DetectorDef::Kind::PType;
  }

  if( slot0_mat.empty() && !slot1_mat.empty() )
    return DetectorDef::Kind::PType;

  if( slot0_mat == "ge" && (slot1_mat == "c" || slot1_mat == "be" || slot1_mat == "al") )
    return DetectorDef::Kind::PType;

  return DetectorDef::Kind::Generic;
}//classify_detector(...)
}//anonymous namespace


std::vector<DetectorDef> parseDetectorTxt( std::istream &input )
{
  // Read all lines, stripping trailing CR/LF (files use CRLF).
  std::vector<std::string> lines;
  {
    std::string ln;
    while( std::getline( input, ln ) )
    {
      while( !ln.empty() && (ln.back() == '\r' || ln.back() == '\n') )
        ln.pop_back();
      lines.push_back( ln );
    }
  }

  std::vector<DetectorDef> recs;
  std::string pending_comment;
  size_t i = 0;

  auto ends_with_par = []( const std::string &tok ) -> bool
  {
    return SpecUtils::iends_with( tok, ".par" );
  };

  while( i < lines.size() )
  {
    std::string s = lines[i];
    SpecUtils::trim( s );
    if( s.empty() )
    {
      ++i;
      continue;
    }
    if( s[0] == '#' )
    {
      pending_comment = s;
      ++i;
      continue;
    }

    // A definition line: split on ',' and find the token ending in ".par".
    std::vector<std::string> parts;
    SpecUtils::split_no_delim_compress( parts, lines[i], "," );
    for( std::string &p : parts )
      SpecUtils::trim( p );

    size_t par_idx = parts.size();
    for( size_t k = 0; k < parts.size(); ++k )
    {
      if( ends_with_par( parts[k] ) )
      {
        par_idx = k;
        break;
      }
    }
    if( par_idx >= parts.size() )
    {
      ++i;
      continue;   // not a definition line
    }

    DetectorDef rec;
    rec.comment = pending_comment;
    rec.serial = serial_from_comment( pending_comment );
    pending_comment.clear();
    rec.name = parts.empty() ? std::string() : parts[0];
    rec.parFile = parts[par_idx];

    // Numeric fields between the name and the .par token: D1..D7 then the
    //  "26" type code (8 values in every observed file; keep whatever is there).
    std::vector<double> scalars;
    for( size_t k = 1; k < par_idx; ++k )
    {
      const double v = parse_num( parts[k] );
      if( !std::isnan( v ) )
        scalars.push_back( v );
    }
    auto scalar = [&]( size_t idx ) -> double
    {
      return (idx < scalars.size()) ? scalars[idx]
                                    : std::numeric_limits<double>::quiet_NaN();
    };
    rec.d1_crystal_diam_mm = scalar( 0 );
    rec.d2_crystal_len_mm  = scalar( 1 );
    rec.d3_window_diam_mm  = scalar( 2 );
    rec.d4_endcap_od_mm    = scalar( 3 );
    rec.d5_endcap_len_mm   = scalar( 4 );
    rec.d6_front_gap_mm    = scalar( 5 );
    rec.d7_side_gap_mm     = scalar( 6 );
    {
      const double tc = scalar( 7 );
      rec.typeCode = std::isnan( tc ) ? 0 : static_cast<int>( tc );
    }
    if( par_idx + 1 < parts.size() )
    {
      const double kv = parse_num( parts[par_idx + 1] );
      rec.kCode = std::isnan( kv ) ? 0 : static_cast<int>( kv );
    }

    ++i;

    // Layer lines: until a blank line, the next record's def line, or a
    //  trailing comment/banner.  A bare ",,," line is a kept empty slot.
    std::vector<LayerEntry> parsed_layers;
    while( i < lines.size() )
    {
      std::string trimmed = lines[i];
      SpecUtils::trim( trimmed );
      if( trimmed.empty() )
        break;
      if( trimmed[0] == '#' )
        break;

      std::vector<std::string> lp;
      SpecUtils::split_no_delim_compress( lp, lines[i], "," );
      for( std::string &p : lp )
        SpecUtils::trim( p );

      bool is_next_def = false;
      for( const std::string &p : lp )
        is_next_def |= ends_with_par( p );
      if( is_next_def )
        break;

      LayerEntry entry;
      entry.material = lp.empty() ? std::string() : lp[0];
      if( lp.size() > 1 )
      {
        const double th = parse_num( lp[1] );
        entry.thickness_mm = std::isnan( th ) ? 0.0 : th;
      }
      if( lp.size() > 2 )
      {
        const double de = parse_num( lp[2] );
        entry.density = std::isnan( de ) ? 0.0 : de;
      }
      parsed_layers.push_back( entry );
      ++i;
    }//while( layer lines )

    for( size_t k = 0; k < rec.layers.size() && k < parsed_layers.size(); ++k )
      rec.layers[k] = parsed_layers[k];

    rec.kind = classify_detector( rec );
    recs.push_back( std::move(rec) );
  }//while( i < lines.size() )

  return recs;
}//parseDetectorTxt(...)


DetectorDef selectDetectorDef( const std::vector<DetectorDef> &defs,
                               const std::string &parFileName )
{
  if( defs.empty() )
    throw std::runtime_error( "selectDetectorDef: no detector records were parsed." );

  const std::string want = SpecUtils::to_lower_ascii_copy(
                              SpecUtils::filename( parFileName ) );

  for( const DetectorDef &d : defs )
  {
    const std::string have = SpecUtils::to_lower_ascii_copy(
                                SpecUtils::filename( d.parFile ) );
    if( !want.empty() && have == want )
      return d;
  }

  if( defs.size() == 1 )
    return defs.front();

  throw std::runtime_error( "selectDetectorDef: none of the " + std::to_string(defs.size())
                            + " records match '" + parFileName + "'." );
}//selectDetectorDef(...)


//=============================================================================
//  .PAR  (binary spatial-efficiency grid)   - a port of par_decode.py
//=============================================================================

ParFile parseParFile( const std::vector<uint8_t> &bytes )
{
  const size_t len = bytes.size();
  if( len < 32 )
    throw std::runtime_error( "parseParFile: file too small to be a parameter grid." );

  // Bounds-checked accessors into the byte buffer.
  auto need = [&]( size_t off, size_t n )
  {
    if( off + n > len )
      throw std::runtime_error( "parseParFile: read past end of file (corrupt or truncated)." );
  };
  auto u16 = [&]( size_t o ){ need(o,2); return read_u16le( &bytes[o] ); };
  auto i32 = [&]( size_t o ){ need(o,4); return read_i32le( &bytes[o] ); };
  auto f32 = [&]( size_t o ){ need(o,4); return read_f32le( &bytes[o] ); };
  auto f64 = [&]( size_t o ){ need(o,8); return read_f64le( &bytes[o] ); };

  // Every record marker (a fixed double) sits `sm_marker_in_record` into its
  //  record; the marker offsets give the framing without hardcoding sizes.
  std::vector<size_t> markers;
  for( size_t o = 0; o + 8 <= len; ++o )
  {
    if( std::fabs( read_f64le( &bytes[o] ) - sm_record_marker ) < 1e-3 )
      markers.push_back( o );
  }
  if( markers.empty() )
    throw std::runtime_error( "parseParFile: no record markers found (not a recognized grid file)." );

  const size_t header_bytes = markers[0] - sm_marker_in_record;
  const size_t record_bytes = (markers.size() > 1)
                                ? (markers[1] - markers[0])
                                : (len - header_bytes);

  ParFile out;
  out.emin_keV = f64( 0x00 );
  out.emax_keV = f64( 0x08 );
  const uint16_t n = u16( 0x10 );
  if( n == 0 )
    throw std::runtime_error( "parseParFile: zero energies in header." );

  out.energies_keV.resize( n );
  for( uint16_t e = 0; e < n; ++e )
    out.energies_keV[e] = f64( 0x12 + 8 * e );

  // Size identity from the derived framing.
  if( len != header_bytes + static_cast<size_t>(n) * record_bytes )
    throw std::runtime_error( "parseParFile: size identity failed (header + N*record != file size)." );

  if( markers.size() != n )
    throw std::runtime_error( "parseParFile: marker count does not match the energy count." );

  // Energies must be strictly ascending (interpolation and PCHIP require it).
  for( uint16_t e = 0; e + 1 < n; ++e )
  {
    if( !(out.energies_keV[e] < out.energies_keV[e + 1]) )
      throw std::runtime_error( "parseParFile: energy grid is not strictly ascending." );
  }

  // Decode each record: a 28-byte header then nrows*ncols uint16, row-major.
  out.grids.resize( n );
  for( uint16_t e = 0; e < n; ++e )
  {
    const size_t rec_off = header_bytes + static_cast<size_t>(e) * record_bytes;

    ParGrid &g = out.grids[e];
    g.ncols = u16( rec_off + 0 );
    g.nrows = u16( rec_off + 2 );
    g.theta_step_rad = f32( rec_off + 12 );
    g.r_step = f32( rec_off + 24 );

    const double marker = f64( rec_off + sm_marker_in_record );
    if( std::fabs( marker - sm_record_marker ) > 1e-3 )
      throw std::runtime_error( "parseParFile: record marker missing (framing derivation failed)." );

    if( g.ncols == 0 || g.nrows == 0 )
      throw std::runtime_error( "parseParFile: degenerate grid dimensions." );
    if( g.theta_step_rad <= 0.0 || g.r_step <= 0.0 )
      throw std::runtime_error( "parseParFile: non-positive grid step." );

    const size_t ncell = static_cast<size_t>(g.nrows) * g.ncols;
    if( 28 + 2 * ncell != record_bytes )
      throw std::runtime_error( "parseParFile: record grid size does not fill the record." );

    const size_t grid_off = rec_off + 28;
    need( grid_off, 2 * ncell );
    g.V.resize( ncell );
    for( size_t c = 0; c < ncell; ++c )
      g.V[c] = read_u16le( &bytes[grid_off + 2 * c] );
  }//for( each record )

  return out;
}//parseParFile( bytes )


ParFile parseParFile( const std::string &path )
{
  // Read the exact bytes (SpecUtils::load_file_data appends a trailing NUL,
  //  which would break the strict file-size identity check).
  std::ifstream input( path.c_str(), std::ios::binary | std::ios::ate );
  if( !input.is_open() )
    throw std::runtime_error( "parseParFile: could not open '" + path + "'." );

  const std::streamsize size = input.tellg();
  if( size <= 0 )
    throw std::runtime_error( "parseParFile: '" + path + "' is empty." );
  input.seekg( 0, std::ios::beg );

  std::vector<uint8_t> bytes( static_cast<size_t>(size) );
  if( !input.read( reinterpret_cast<char *>(bytes.data()), size ) )
    throw std::runtime_error( "parseParFile: could not read '" + path + "'." );

  return parseParFile( bytes );
}//parseParFile( path )


//=============================================================================
//  ParEfficiency  - the ground-truth grid evaluator
//=============================================================================

namespace
{
/** Bilinear interpolation of one grid in (ln_d, theta), in the raw-V domain -
 exactly par_decode.py's `interp_grid`, except that "no data here" cells are
 excluded from the blend rather than being averaged in as if they were data.

 A cell of V == 0 is a SENTINEL, not a measurement: eff = 10^(-V/1000) would be
 exactly 1.0, i.e. 100% absolute efficiency, which no real geometry reaches.
 Those cells mark source positions INSIDE THE ENDCAP - physically unreachable,
 hence a sentinel rather than a number.  Re-expressing the zone boundary in
 cylindrical coords recovers the can to within the grid quantization (42.5 mm
 lateral vs a 43.8 mm endcap radius, 81.0 mm deep vs an 83.8 mm endcap length).
 The zone is energy-INDEPENDENT, is contiguous from the innermost row outward,
 and exists only behind the face plane: zero such cells at theta <= 90 deg, with
 onset at the first column past it (92.5 deg).

 Sentinels are excluded and the remaining weight renormalized, so approaching
 the body from outside degrades toward the nearest real data instead of jumping
 to 100%.  `all_sentinel` reports the case where there is no data to fall back
 on, letting callers refuse to serve a value rather than invent one. */
double interp_grid_V( const ParGrid &g, const double ln_d, const double theta_rad,
                      bool *all_sentinel = nullptr )
{
  const int nrows = static_cast<int>( g.nrows );
  const int ncols = static_cast<int>( g.ncols );

  const double row_f = ln_d / g.r_step;
  const double col_f = std::fabs( theta_rad ) / g.theta_step_rad;

  int r0 = static_cast<int>( row_f );
  int c0 = static_cast<int>( col_f );
  r0 = std::max( 0, std::min( r0, nrows - 2 ) );
  c0 = std::max( 0, std::min( c0, ncols - 2 ) );
  const int r1 = r0 + 1;
  const int c1 = c0 + 1;

  double fr = row_f - r0;
  double fc = col_f - c0;
  fr = std::max( 0.0, std::min( 1.0, fr ) );
  fc = std::max( 0.0, std::min( 1.0, fc ) );

  const double vs[4] = { static_cast<double>( g.V[r0 * ncols + c0] ),
                         static_cast<double>( g.V[r0 * ncols + c1] ),
                         static_cast<double>( g.V[r1 * ncols + c0] ),
                         static_cast<double>( g.V[r1 * ncols + c1] ) };
  const double ws[4] = { (1.0 - fr) * (1.0 - fc), (1.0 - fr) * fc,
                         fr * (1.0 - fc), fr * fc };

  double sum = 0.0, wsum = 0.0;
  for( size_t i = 0; i < 4; ++i )
  {
    if( vs[i] == 0.0 )     // sentinel: inside the detector body, no data
      continue;
    sum += ws[i] * vs[i];
    wsum += ws[i];
  }

  if( all_sentinel )
    *all_sentinel = (wsum <= 0.0);

  if( wsum <= 0.0 )
    return 0.0;            // caller decides; 0 keeps the old (eff = 1) shape

  return sum / wsum;
}//interp_grid_V(...)
}//anonymous namespace


ParEfficiency::ParEfficiency( ParFile par )
  : m_par( std::move(par) )
{
  if( m_par.energies_keV.size() < 2 || m_par.grids.size() != m_par.energies_keV.size() )
    throw std::runtime_error( "ParEfficiency: need at least two parallel energy/grid entries." );
}//ParEfficiency(...)


double ParEfficiency::efficiency( double energy_keV, double dist_from_face_mm,
                                  double theta_rad, bool *no_data ) const
{
  const std::vector<double> &energies = m_par.energies_keV;
  const size_t ne = energies.size();

  const double ln_d = std::log( std::max( dist_from_face_mm, 1.0 ) );

  // Interpolated V at every energy node for this (distance, angle), then a
  //  shape-preserving monotone cubic (PCHIP) across ln(E).  Log-log is the same
  //  domain the reference tool uses; PCHIP (vs a plain cubic spline, which
  //  overshoots the low-E efficiency knee) is our deliberate, overshoot-safe
  //  choice.  We therefore differ from the reference ".ecc" between nodes, worst
  //  ~1.4% at 200 keV - the reference behaves like a local forward-quadratic in
  //  (V, ln E), and an independent Monte-Carlo arbiter says PCHIP is the more
  //  physical of the two.  Matching the reference more closely would mean being
  //  less physically accurate; see the DetectorEffG2kPar.h header for the
  //  measured comparison and the ~0.23% quantization floor on any such match.
  // The V == 0 sentinel zone is energy-independent (verified across the corpus:
  //  identical row extents at every energy node), so "no data" is a property of
  //  the position alone and any node reporting it means all of them do.
  std::vector<double> lnE( ne ), Vnode( ne );
  bool any_no_data = false;
  for( size_t e = 0; e < ne; ++e )
  {
    bool cell_no_data = false;
    lnE[e] = std::log( energies[e] );
    Vnode[e] = interp_grid_V( m_par.grids[e], ln_d, theta_rad, &cell_no_data );
    any_no_data |= cell_no_data;
  }

  if( no_data )
    *no_data = any_no_data;

  double V;
  if( ne == 1 )
  {
    V = Vnode[0];
  }else
  {
    const ceelo::Pchip curve( lnE, Vnode );      // clamps outside [Emin,Emax]
    V = curve( std::log( std::max( energy_keV, 1e-9 ) ) );
  }

  return std::pow( 10.0, -V / 1000.0 );
}//ParEfficiency::efficiency(...)


double ParEfficiency::airTransmission( double energy_keV, double dist_mm,
                                       double pressure_atm ) const
{
  if( pressure_atm <= 0.0 )
    return 1.0;

  // Rayleigh (coherent) scatter is EXCLUDED from the removal coefficient: it is
  //  elastic, so a deflected photon still arrives at full energy and still
  //  lands in the full-energy peak.  Over an air path it is a redirection, not
  //  a removal, and for a detector subtending a small solid angle the in- and
  //  out-scatter very nearly cancel.  Three independent reasons this is right:
  //   * Measured.  The corpus has PAIRED air/vacuum runs at the same position,
  //     so `ecc_air/ecc_vac` isolates the reference tool's air factor exactly
  //     (grid, spatial and energy interpolation all cancel).  Fitting
  //     -ln(T) = d*(mu_noRS + w*mu_rs) over 45-2000 keV at 15/20/100/250/999.9
  //     mm gives w ~ 0, not w = 1: rms residual 0.0007-0.008% at w=0 versus
  //     0.0007-0.066% at w=1, the gap growing with path length exactly as a
  //     mu error must.
  //   * Physical.  This is CeeLo's own position - see AirAttenuation::
  //     AnalyticNoScatter ("excludes Rayleigh (isotropic cancellation
  //     argument)") and fep_survival_removal_mu() = mu_total - mu_rs - f*mu_cs.
  //   * Consistent.  InterSpec's own transmission_length_coefficient_air()
  //     already omits it (MassAttenuationTool returns 0 for RayleighScatter).
  //  The second-order effect - elastic deflection lengthening the remaining
  //  path - is negligible in air: CeeLo's rayleigh_deflection_loss_fraction()
  //  gives h ~ 0.03 at 45 keV over 1 m, i.e. a 0.009% loss (0.16% even at
  //  10 m).  It only matters for thick high-Z layers, not for air.
  //  Including mu_rs would UNDER-report efficiency (-0.28% at 45 keV/1 m,
  //  -1.4% at 5 m), biasing any derived activity high.
  //
  // The coefficients are 1/cm at the reference dry-air density; linear
  //  attenuation scales with density, i.e. with pressure.
  static const ceelo::Material air = ceelo::make_Air();
  const ceelo::MacroscopicXS xs = air.macroscopic_xs( energy_keV / 1000.0 ); // MeV
  const double mu_per_cm = xs.mu_total() - xs.mu_rs;
  const double d_cm = dist_mm / 10.0;
  return std::exp( -mu_per_cm * pressure_atm * d_cm );
}//ParEfficiency::airTransmission(...)


//=============================================================================
//  CeeLo geometry + response builder  (no Monte Carlo)
//=============================================================================

namespace
{
/** Builds a MaterialSpec for a layer-stack material token.

 The tokens are elemental symbols, so rather than keep a hand-written table of
 the ones our sample files happen to use, we resolve the symbol through
 SandiaDecay (which knows every element, and rejects non-elements - "brass",
 "air", "mylar" and friends all return null).  That matters: a hand table
 silently dropped a whole magnesium endcap on a detector in the sample set,
 which is the sort of thing that only shows up when someone imports a detector
 nobody tested with.

 Density comes from the file, which supplies one on every layer line seen; only
 if it is absent do we need a fallback, and then we use CeeLo's built-in for the
 handful of shielding elements it ships (also picking up their exact reference
 densities) and otherwise fail rather than invent a number.  Returns false for
 anything unresolvable - see the caller, which surfaces the skipped layer. */
bool material_for_token( const std::string &token_in, double parsed_density,
                         ceelo::MaterialSpec &out )
{
  const std::string t = SpecUtils::to_lower_ascii_copy( token_in );
  if( t.empty() )
    return false;

  // The file's own density wins; these are the fallbacks when it has none.
  // `make_HPGe` is also how the crystal itself is specified, so a "ge" layer
  //  must map to the same MaterialSpec name or the material table de-dupe
  //  (by name) would add a second, near-identical germanium.
  const ceelo::Material *fallback = nullptr;
  static const ceelo::Material ge = ceelo::make_HPGe();
  static const ceelo::Material al = ceelo::make_Aluminum();
  static const ceelo::Material cu = ceelo::make_Copper();
  static const ceelo::Material fe = ceelo::make_Iron();
  static const ceelo::Material sn = ceelo::make_Tin();
  static const ceelo::Material pb = ceelo::make_Lead();
  static const ceelo::Material w  = ceelo::make_Tungsten();
  if( t == "ge" ) fallback = &ge;
  else if( t == "al" ) fallback = &al;
  else if( t == "cu" ) fallback = &cu;
  else if( t == "fe" ) fallback = &fe;
  else if( t == "sn" ) fallback = &sn;
  else if( t == "pb" ) fallback = &pb;
  else if( t == "w" ) fallback = &w;

  if( fallback && !(parsed_density > 0.0) )
  {
    out = ceelo::MaterialSpec::from( *fallback );
    return true;
  }

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element * const el = db ? db->element( t ) : nullptr;
  if( !el || (el->atomicNumber <= 0) || (el->atomicNumber > 255) )
    return false;

  if( !(parsed_density > 0.0) )
    return false;      //an element we have no reference density for

  out = ceelo::MaterialSpec();
  out.name = fallback ? fallback->name() : el->symbol;
  out.density_g_per_cm3 = parsed_density;
  out.composition.push_back(
      ceelo::MaterialComponent{ static_cast<uint8_t>(el->atomicNumber), 1.0 } );
  return true;
}//material_for_token(...)


/** Drops geometry features until the descriptor is buildable.  Extends
 CeeLoUtils::relaxGeometryFeatures (fillet / round-tip) by also dropping a bore
 hole that does not fit - legitimate here because the bore is an estimate not
 present in DETECTOR.txt, and the node identity cancels the kernel anyway. */
void sanitize_geometry( ceelo::GeometryDescriptor &gd )
{
  std::vector<std::string> warnings;
  CeeLoUtils::relaxGeometryFeatures( gd, warnings );

  if( !gd.problems().empty() && gd.bore )
  {
    gd.bore.reset();
    CeeLoUtils::relaxGeometryFeatures( gd, warnings );
  }

  if( !gd.problems().empty() && gd.bullet_radius_cm != 0.0 )
    gd.bullet_radius_cm = 0.0;

  const std::vector<ceelo::GeometryProblem> remaining = gd.problems();
  if( !remaining.empty() )
  {
    std::string msg = "makeDrf: the crystal geometry is unusable:";
    for( const ceelo::GeometryProblem p : remaining )
      msg += std::string(" ") + ceelo::to_string( p );
    throw std::runtime_error( msg );
  }
}//sanitize_geometry(...)


/** Builds a CeeLo GeometryDescriptor from a parsed DetectorDef.  Reuses only
 the field->CeeLo structure of the research prototype; every number comes from
 `def`.  Bore + fillet are not in DETECTOR.txt, so conservative defaults are
 used (and dropped by sanitize_geometry if they don't fit).

 Appends to `warnings` for any housing layer that could not be modeled, so the
 omission reaches the user instead of vanishing (see makeDrf, which folds these
 into the DRF description). */
ceelo::GeometryDescriptor build_geometry( const DetectorDef &def,
                                         std::vector<std::string> &warnings )
{
  using namespace ceelo;

  const double diam_mm = def.d1_crystal_diam_mm;
  const double len_mm = def.d2_crystal_len_mm;
  if( !(diam_mm > 0.0) || !(len_mm > 0.0) )
    throw std::runtime_error( "makeDrf: DETECTOR.txt is missing the crystal diameter/length." );

  const double radius_cm = 0.05 * diam_mm;    // (diam/2)/10
  const double length_cm = 0.1 * len_mm;

  GeometryDescriptor gd;
  gd.set_dimensions( CylinderDims{ radius_cm, length_cm } );
  gd.reference_point = ReferencePoint::EndcapFront;

  // Material table, deduped by name; HPGe (crystal) is index 0.
  std::vector<std::string> mat_names;
  auto add_material = [&]( const MaterialSpec &spec ) -> int
  {
    for( size_t i = 0; i < gd.materials.size(); ++i )
    {
      if( gd.materials[i].name == spec.name )
        return static_cast<int>( i );
    }
    gd.materials.push_back( spec );
    return static_cast<int>( gd.materials.size() - 1 );
  };

  gd.crystal_material_index = add_material( MaterialSpec::from( make_HPGe() ) );

  // Dead layers are the Ge entries in slots 0 / 3 / 6 (front / side / back);
  //  they are already folded into the grid, but shape the crystal solid.
  const double dead_front_cm = 0.1 * def.layers[0].thickness_mm;
  const double dead_side_cm  = 0.1 * def.layers[3].thickness_mm;
  const double dead_back_cm  = 0.1 * def.layers[6].thickness_mm;
  if( dead_front_cm > 0.0 || dead_side_cm > 0.0 || dead_back_cm > 0.0 )
    gd.dead_layer = DeadLayerConfig{ dead_front_cm, dead_side_cm, dead_back_cm };

  // Front vacuum gap (D6), widening to the side gap (D7) - matches the
  //  prototype's spacer and sets the endcap-front offset.
  const double front_gap_cm = 0.1 * def.d6_front_gap_mm;
  const double side_gap_cm = 0.1 * def.d7_side_gap_mm;
  if( front_gap_cm > 0.0 )
  {
    MaterialSpec vac;
    vac.name = "Vacuum";
    vac.density_g_per_cm3 = 1.0e-25;
    vac.composition.push_back( MaterialComponent{ 1, 1.0 } );
    const int vac_idx = add_material( vac );

    LayerSpec gap;
    gap.material_index = vac_idx;
    gap.front_thickness_cm = front_gap_cm;
    gap.side_thickness_cm = (side_gap_cm > front_gap_cm) ? (side_gap_cm - front_gap_cm) : 0.0;
    gap.z_start_cm = 0.0;
    gap.z_end_cm = length_cm;
    gd.layers.push_back( gap );
  }

  // Housing layers: front path uses slots 1,2 (front thickness); side path
  //  uses slots 4,5 (side thickness).  Back path (slots 7,8) is not modeled by
  //  CeeLo (theta>90 is handled from the grid directly).
  auto add_housing = [&]( int slot, bool is_front )
  {
    const LayerEntry &le = def.layers[slot];
    if( le.material.empty() || !(le.thickness_mm > 0.0) )
      return;
    MaterialSpec spec;
    if( !material_for_token( le.material, le.density, spec ) )
    {
      // Not fatal - the grid holds the true efficiency at its own nodes either
      //  way - but it does degrade the kernel that shapes the between-node and
      //  near-field model, so say so rather than drop it silently.
      warnings.push_back( "layer material '" + le.material + "' ("
                          + SpecUtils::printCompact( le.thickness_mm, 4 )
                          + " mm) is not a recognized element, so that layer is"
                            " not modeled" );
      return;
    }
    const int idx = add_material( spec );
    LayerSpec ls;
    ls.material_index = idx;
    ls.front_thickness_cm = is_front ? (0.1 * le.thickness_mm) : 0.0;
    ls.side_thickness_cm  = is_front ? 0.0 : (0.1 * le.thickness_mm);
    ls.z_start_cm = 0.0;
    ls.z_end_cm = length_cm;
    gd.layers.push_back( ls );
  };
  add_housing( 1, true );    // front window
  add_housing( 2, true );    // front endcap
  add_housing( 4, false );   // side absorber / liner
  add_housing( 5, false );   // side casing

  // Bore hole + fillet: not in DETECTOR.txt.  Add a conservative coaxial bore
  //  for the coaxial families only; sanitize_geometry drops it if it doesn't
  //  fit (safe - the node identity cancels the kernel).
  if( (def.kind == DetectorDef::Kind::NCoax || def.kind == DetectorDef::Kind::PType)
      && length_cm > 1.5 )
  {
    const double bore_depth_cm = std::max( 0.5, length_cm - 0.8 );
    gd.bore = BoreHoleConfig{ 0.4, bore_depth_cm, true };
  }
  gd.bullet_radius_cm = 0.8;

  sanitize_geometry( gd );
  return gd;
}//build_geometry(...)
}//anonymous namespace


std::shared_ptr<DetectorPeakResponse> makeDrf( const ParFile &par, const DetectorDef &def )
{
  using namespace ceelo;

  if( par.energies_keV.size() < 2 || par.grids.size() != par.energies_keV.size() )
    throw std::runtime_error( "makeDrf: the parameter grid has fewer than two energies." );

  const std::vector<double> &energies = par.energies_keV;
  const size_t ne = energies.size();

  // ---- 1) Geometry from the parsed record ---------------------------------
  std::vector<std::string> geom_warnings;
  GeometryDescriptor gd = build_geometry( def, geom_warnings );

  // ---- 2) Response shell + stored mu tables (self-contained) --------------
  auto resp = std::make_shared<DetectorResponse>();
  resp->descriptor = gd;
  for( size_t i = 0; i < gd.materials.size(); ++i )
    resp->mu_tables.push_back( MuTable::sample( gd.materials[i].to_material(),
                                                static_cast<int>(i) ) );

  resp->provenance.method = ProductionMethod::CurveTransfer;   // no Monte Carlo
  resp->provenance.profile = ResponseProfile::General;         // we carry a near model
  resp->provenance.kernel_n_rays = 2048;
  resp->provenance.detector_name = def.name;
  resp->provenance.valid_e_min_keV = energies.front();
  resp->provenance.valid_e_max_keV = energies.back();
  // min_distance_cm is set in step 7, once the near-field ladder's floor is
  //  known: it is that floor, so a closer query is flagged rather than served
  //  from an extrapolation of the bottom node.

  resp->finalize();   // builds geometry + material->mu lookups (empty tables tolerated)

  const ParEfficiency parEff( par );

  // Angular nodes: cosines of the grid's polar columns on [0, 90 deg], stored
  //  ascending in cos(theta) (theta=90 -> 0 first ... theta=0 -> 1 last).  The
  //  same node set is shared by eta_fep and near_field so the far term cancels.
  const double theta_step = par.grids.front().theta_step_rad;
  const double half_pi = 0.5 * 3.14159265358979323846;
  size_t nhalf = static_cast<size_t>( std::floor( half_pi / theta_step ) ) + 1;
  nhalf = std::min<size_t>( nhalf, par.grids.front().ncols );
  if( nhalf < 2 )
    throw std::runtime_error( "makeDrf: too few angular columns to build a response." );

  std::vector<double> cos_thetas( nhalf );
  for( size_t j = 0; j < nhalf; ++j )
  {
    const double angle = static_cast<double>( nhalf - 1 - j ) * theta_step;
    cos_thetas[j] = std::cos( angle );
  }
  const size_t nc = cos_thetas.size();

  // Far reference (crystal frame) beyond the largest expected query distance,
  //  so every corpus query lands inside the filled near-field region.
  const double d_ref_cm = 300.0;
  const Eigen::Vector3d face = CeeLoUtils::detectorFacePosition( gd );

  // Convert a crystal-frame source position to the grid's face-frame index.
  //  The grid is SPHERICAL: the row is the RADIAL range |v| from the endcap-face
  //  centre to the source, and the column is the polar angle from the axis.  It
  //  is NOT the axial height (the perpendicular projection -v.z()); see the
  //  header's "distance-axis convention" note for the three independent lines of
  //  evidence.  On the axis the two coincide, which is why on-axis points match
  //  under either reading.
  auto face_sample = [&]( const Eigen::Vector3d &src, double &d_face_mm, double &theta_face )
  {
    const Eigen::Vector3d v = src - face;
    const double n = v.norm();
    d_face_mm = n * 10.0;   // radial range from the face centre (mm)
    theta_face = (n > 0.0) ? std::acos( std::max( -1.0, std::min( 1.0, -v.z() / n ) ) ) : 0.0;
  };

  const std::vector<double> edges = gd.crystal_k_edges( energies.front(), energies.back() );

  // ---- 6) eta_fep at the far reference; cache ln_eta for the near term ------
  EtaTable &eta = resp->eta_fep;
  eta.energies_keV = energies;
  eta.cos_thetas = cos_thetas;
  eta.phis_deg.clear();
  eta.edges_keV = edges;
  eta.ln_eta.assign( ne * nc, 0.0 );
  eta.frac_sigma.assign( ne * nc, 0.0 );

  for( size_t c = 0; c < nc; ++c )
  {
    const Eigen::Vector3d src = source_position( d_ref_cm, cos_thetas[c], 0.0 );
    const ApertureQuadrature q = resp->make_quadrature( src );   // energy-independent

    double d_face_mm = 0.0, theta_face = 0.0;
    face_sample( src, d_face_mm, theta_face );

    for( size_t e = 0; e < ne; ++e )
    {
      const double E = energies[e];
      const double K = resp->kernel_K( E, q, MuChoice::Total );
      const double eps = parEff.efficiency( E, d_face_mm, theta_face );
      eta.ln_eta[eta.index( e, c, 0 )] = std::log( std::max( eps, 1e-300 ) / std::max( K, 1e-300 ) );
    }
  }//for( each cos-theta )
  eta.finalize();

  // Read the far term back through eval_ln (exact at nodes) so the near-field
  //  subtraction below uses the SAME value eval-time will add - making the far
  //  term cancel algebraically no matter how finalize()/eval_ln evolve.
  std::vector<double> ln_eta_far( ne * nc, 0.0 );   // [e*nc + c]
  for( size_t c = 0; c < nc; ++c )
  {
    for( size_t e = 0; e < ne; ++e )
    {
      bool clamped = false;
      ln_eta_far[e * nc + c] = eta.eval_ln( energies[e], cos_thetas[c], 0.0, clamped );
    }
  }

  // ---- 7) Near field: reproduce close-in points via the node identity ------
  // Distance ladder (crystal frame), geometric from near-contact up to the far
  //  reference; the top node coincides with the far reference so ln_N -> 0 there.
  //
  // The floor must clear the ENDCAP FACE PLANE, not merely be small.  The ladder
  //  is in the crystal frame while the grid is indexed in the face frame, and
  //  those origins differ by `endcap_front_offset_cm`.  Writing the face-frame
  //  offset vector as v = src - face, its length is minimized on-axis:
  //      |v|^2 = d^2 - 2*d*cos(theta)*off + off^2   ->   |v| = |d - off| at theta=0,
  //  so a ladder node with d < off sits BEHIND the face plane: the radial range
  //  collapses toward zero (clamping onto the grid's 1 mm bottom row) and
  //  `theta_face` flips past 90 deg - on-axis at d=0.5 cm with off=0.57 cm it
  //  reaches 180 deg.  Such nodes are not merely inaccurate, they are anchors of
  //  the PCHIP ladder, so they corrupt legitimate neighbouring queries too:
  //  measured before this floor existed, dispatch disagreed with the raw grid by
  //  158% at 2 cm and 27% at 5 cm, and returned FEP efficiencies near 0.6 at
  //  7 MeV - physically impossible for this crystal, and unflagged.
  //
  // Flooring d at off + margin makes |v| >= margin everywhere on the ladder.
  //  The margin is in the FACE frame and is kept well above the grid's 1 mm
  //  bottom row so the closest node still lands on real data rather than the
  //  row clamp.
  //
  // What the floor does NOT fix: the ladder is keyed on crystal-frame
  //  (d, cos_theta) - those are the axes `common_eval` reduces a query to - while
  //  the value at each node is read from the grid in face-frame coords.  The map
  //  between the frames is a translation by `off`, which is near-identity in the
  //  far field but strongly nonlinear once d ~ off: grazing ladder columns reach
  //  theta_face of ~126 deg (still real grid data - `interp_grid_V` interpolates
  //  the full column span, it does not clamp at 90 deg).  So ln_n, expressed on
  //  the crystal-frame axes it is interpolated over, has curvature near contact
  //  that this lattice cannot resolve.  Measured against the raw grid, dispatch
  //  is exact AT the nodes but errs by tens to thousands of percent between them
  //  for face-frame queries inside ~5 cm at grazing angles.  See the header note;
  //  callers needing contact geometry off-axis must not use this response yet.
  const size_t nd = 18;
  const double d_face_margin_cm = 0.2;   // >= 2 mm of face-frame standoff
  const double d_min_cm = -face.z() + d_face_margin_cm;   // face is (0,0,-offset)
  std::vector<double> dists_cm( nd );
  for( size_t i = 0; i < nd; ++i )
  {
    const double f = static_cast<double>(i) / (nd - 1);
    dists_cm[i] = std::exp( std::log(d_min_cm) + f * (std::log(d_ref_cm) - std::log(d_min_cm)) );
  }
  dists_cm.back() = d_ref_cm;

  NearFieldModel &nf = resp->near_field;
  nf.energies_keV = energies;
  nf.cos_thetas = cos_thetas;
  nf.dists_cm = dists_cm;
  nf.ln_n.assign( ne * nc * nd, 0.0 );
  nf.frac_sigma.assign( ne * nc * nd, 0.0 );

  // Each angular column's ladder is TRUNCATED at its first valid rung: rungs
  //  whose face-frame position falls in the grid's no-data zone are not filled
  //  from the file, because there is nothing there to reproduce (the grid stores
  //  V = 0, which decodes to eff = 1.0 - 100%, impossible).  Filling them with
  //  that literal value anchors the PCHIP ladder on an impossible number and
  //  corrupts legitimate neighbouring queries; that is how off-axis dispatch
  //  came to disagree with the raw grid by thousands of percent.
  //
  // `NearFieldModel` has ONE shared distance axis (index(e,c,d) is rectangular,
  //  and is serialized that way), so a literally ragged ladder is not
  //  representable.  It does not need to be: `Pchip::operator()` clamps flat
  //  below its first node (`xq <= x_.front() -> y_.front()`), so truncating a
  //  column is EXACTLY equivalent to filling its inner rungs with the first
  //  valid value.  We therefore write that hold, which yields a truncated
  //  interpolant without touching CeeLo's structure or file format.
  //
  // The no-data zone is the ENDCAP VOLUME itself - source positions inside the
  //  housing, which are physically unreachable, so the reference tool stores a
  //  sentinel instead of a number.  Converting the zone's boundary from the
  //  grid's (radial, theta) into cylindrical coords recovers the can: lateral
  //  extent 42.5 mm vs the 43.8 mm endcap radius, depth 81.0 mm vs the 83.8 mm
  //  endcap length (LAB06: 31.9 vs 38.1 mm, and 111.3 vs 133.4 mm).  It
  //  therefore only exists behind the face plane - zero such cells at
  //  theta <= 90 deg, onset at 92.5 deg.
  //
  // So this costs near-field reach only in the grazing columns, where the
  //  crystal-frame ladder position translates to a face-frame point inside that
  //  can: on the corpus 18/37 columns lose nothing at all, the on-axis ladder
  //  still reaches 2 mm from the endcap face, and only the last four columns
  //  (82.5-90 deg) are held back to ~4.5 cm.  48 of 666 rungs on one detector,
  //  39 of 666 on the other.
  size_t n_truncated = 0;
  for( size_t c = 0; c < nc; ++c )
  {
    bool have_valid = false;
    for( size_t k = 0; k < nd; ++k )
    {
      const size_t di = nd - 1 - k;   // far -> near, so the hold source exists
      const Eigen::Vector3d src = source_position( dists_cm[di], cos_thetas[c], 0.0 );

      double d_face_mm = 0.0, theta_face = 0.0;
      face_sample( src, d_face_mm, theta_face );

      // Every rung must clear the endcap face by the margin; the d_min_cm floor
      //  above guarantees it.  theta_face is deliberately NOT bounded by 90 deg:
      //  translating a grazing near-contact point into the face frame
      //  legitimately swings it behind the face plane, and the grid does carry
      //  real columns out there (its span is the full 180 deg).
      assert( d_face_mm >= 10.0 * 0.5 * d_face_margin_cm );

      bool no_data = false;
      parEff.efficiency( energies[0], d_face_mm, theta_face, &no_data );

      if( no_data )
      {
        // Below this column's first valid rung: hold it (== truncation, see
        //  above).  Keep the held rung's sigma inflated - the position is
        //  uncharacterized by the file, not modelled by us, so a query landing
        //  here should be honestly uncertain rather than confidently wrong.
        ++n_truncated;
        for( size_t e = 0; e < ne; ++e )
        {
          const size_t idx = nf.index( e, c, di );
          nf.ln_n[idx] = have_valid ? nf.ln_n[nf.index( e, c, di + 1 )] : 0.0;
          nf.frac_sigma[idx] = 1.0;   // 100%: no data underlies this rung
        }
        continue;
      }

      const ApertureQuadrature q = resp->make_quadrature( src );
      for( size_t e = 0; e < ne; ++e )
      {
        const double E = energies[e];
        const double K = resp->kernel_K( E, q, MuChoice::Total );
        const double eps = parEff.efficiency( E, d_face_mm, theta_face );
        const double ln_ratio = std::log( std::max( eps, 1e-300 ) / std::max( K, 1e-300 ) );
        nf.ln_n[nf.index( e, c, di )] = ln_ratio - ln_eta_far[e * nc + c];
      }
      have_valid = true;
    }//for( each distance, far -> near )
  }//for( each cos-theta )

  if( n_truncated )
    geom_warnings.push_back( "near-field ladder truncated at "
                             + std::to_string(n_truncated) + " of "
                             + std::to_string(nc * nd) + " rungs, where the"
                             " face-frame position lies behind the endcap in the"
                             " region the file does not characterize (grazing"
                             " angles at close range only)" );

  // Near gate covers the whole filled region (break at the far reference), so
  //  every query inside d_ref_cm picks up the near-field term.
  nf.break_cos_thetas = { 0.0, 1.0 };
  nf.break_d_cm.assign( ne * 2, d_ref_cm );
  nf.finalize();

  // Below the ladder's first node there is no data to interpolate, only the
  //  bottom-node extrapolation `ln_boost` clamps to - and a source that close is
  //  at or inside the endcap.  Record the floor so the validity limit travels
  //  with the response (it is serialized, and callers can read it).
  //
  // NOTE this does NOT by itself make such a query come back flagged: the only
  //  `NearFieldUnmodeled` path (DetectorResponse.cpp, the d < d_gate branch)
  //  fires when `near_field` is EMPTY, and ours never is.  min_distance_cm only
  //  widens d_gate, i.e. which queries take the near-field branch at all.  An
  //  honest flag for close/grazing queries still needs to be added - see the
  //  no-data (endcap-volume) note in the ladder fill above.
  //
  // Also note d_min_cm is a floor on the ON-AXIS reach only, not a validity
  //  radius: it is measured along the axis from the endcap FACE, whereas the
  //  no-data region is the endcap CAN (~42 mm lateral, ~81 mm deep on 18211381).
  //  A point further than d_min_cm from the face centre can therefore still be
  //  inside the housing once it is past 90 deg - validity off-axis is bounded by
  //  the can's envelope, not by a sphere of radius d_min_cm.  Measured first
  //  valid radial range versus face-frame polar angle (18211381): 1 mm at
  //  0-90 deg, then 42 mm at 95-100 deg, 50 mm at 120 deg, 84 mm at 150 deg.
  resp->provenance.min_distance_cm = d_min_cm;

  // ---- 8) Total efficiency: NOT characterized by these files ---------------
  // The grid is full-energy-peak only, so there is nothing to build a total
  //  efficiency from.  Say so explicitly: the default `KernelExact` tier would
  //  serve a BARE-CRYSTAL kernel - no housing, no peak-to-total - and callers
  //  cannot distinguish that from a verified total.  It is not a small error
  //  either; at low energy the bare kernel falls below this response's own FEP
  //  efficiency, which is impossible.  `NotCharacterized` makes eps_total
  //  refuse, which in turn makes `hasAnyTotalEfficiencyInfo()` false, so
  //  cascade summing correctly declines to run on this DRF.
  resp->tot_eff.tier = TotEffTier::NotCharacterized;
  resp->tot_eff.finalize();

  // ---- Attach to a DetectorPeakResponse ------------------------------------
  auto drf = std::make_shared<DetectorPeakResponse>();

  std::string name = def.name;
  if( name.empty() )
    name = SpecUtils::filename( def.parFile );
  std::string descrip = "Imported from a detector-characterization parameter file";
  if( !def.serial.empty() )
    descrip += " (S/N " + def.serial + ")";
  descrip += ".  Full-energy-peak efficiency only; total efficiency is not"
             " characterized, so cascade-summing corrections are unavailable.";
  for( const std::string &w : geom_warnings )
    descrip += "  Note: " + w + ".";

  drf->setName( name );
  drf->setDescription( descrip );
  drf->setDrfSource( DetectorPeakResponse::DrfSource::CharacterizationParFile );

  // Sample an on-axis far-field legacy curve so the DRF is valid/serializable
  //  even without the CeeLo response (setEfficiencyPoints sets the energy range
  //  from the sampled points; the DRF source we just set is preserved).
  const std::shared_ptr<const DetectorResponse> const_resp = resp;
  CeeLoUtils::setLegacyEfficiencyFromResponse( *drf, const_resp );

  // With a CeeLo response attached its descriptor is authoritative, so the
  //  explicit setGeometry is redundant; kept for the response-less code path.
  drf->setGeometry( std::make_shared<const GeometryDescriptor>( gd ) );
  drf->setCeeloResponse( const_resp );

  return drf;
}//makeDrf(...)


std::shared_ptr<DetectorPeakResponse> makeDrfFromFiles( const std::string &parPath,
                                                        const std::string &detectorTxtPath )
{
  const ParFile par = parseParFile( parPath );

  std::ifstream txt( detectorTxtPath.c_str(), std::ios::binary );
  if( !txt.is_open() )
    throw std::runtime_error( "makeDrfFromFiles: could not open '" + detectorTxtPath + "'." );
  const std::vector<DetectorDef> defs = parseDetectorTxt( txt );

  const DetectorDef def = selectDetectorDef( defs, SpecUtils::filename( parPath ) );
  return makeDrf( par, def );
}//makeDrfFromFiles(...)

}//namespace DetEffG2kPar
