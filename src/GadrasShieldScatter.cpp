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

 The interpolation algorithm and binary file layout are ported from
 gadras_19.5.3/CGADFunc/ShieldScatter/src/ (Computer.cpp, Database.cpp,
 DatabaseParser.cpp). The energy-resonance handling, log-space
 interpolation, geometry blending, and below-minimum-AD extrapolation
 are reproduced as in the GADRAS source.
 */

#include "InterSpec_config.h"

#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/GadrasShieldScatter.h"

using namespace std;


// All implementation lives in this anonymous namespace; the public
// interface (GadrasShieldScatter and its private Database forward declaration)
// is in the header.
namespace
{
  /** Mirrors the GADRAS Database class: holds the parsed sandia.shieldscatter.db
   contents in memory. Sizes:
     energies[E], atomicNumbers[Z], arealDensities[AD], hydrogenMassFractions[H]
     k/l/mEdgeEnergies[Z], k/l/mShellRegionUpperBounds[Z]
     logScatter{Point,Basis,Slab}[H][AD][Z][E][G]   (always stored as ln)
     {point,basis,slab}ShapeFactor[Z][AD]
   */
  struct ShieldScatterDb
  {
    int groupCount = 0;
    int energyCount = 0;
    int atomicNumberCount = 0;
    int arealDensityCount = 0;
    int hydrogenMassFractionCount = 0;

    vector<double> energies;
    vector<double> atomicNumbers;
    vector<double> arealDensities;
    vector<double> hydrogenMassFractions;

    vector<double> kEdgeEnergies;
    vector<double> kShellRegionUpperBounds;
    vector<double> lEdgeEnergies;
    vector<double> lShellRegionUpperBounds;
    vector<double> mEdgeEnergies;
    vector<double> mShellRegionUpperBounds;

    // Indexed [hydrogenIdx][arealDensityIdx][atomicNumberIdx][energyIdx][groupIdx]
    // (matches GADRAS access order, even though the file stores in a different
    // index order, see ConvertColumnMajorArrayToLogMatrix).
    vector< vector< vector< vector< vector<double> > > > > logScatterPoint;
    vector< vector< vector< vector< vector<double> > > > > logScatterBasis;
    vector< vector< vector< vector< vector<double> > > > > logScatterSlab;

    // Indexed [atomicNumberIdx][arealDensityIdx]
    vector< vector<double> > pointShapeFactor;
    vector< vector<double> > basisShapeFactor;
    vector< vector<double> > slabShapeFactor;
  };//struct ShieldScatterDb


  /** Sequential little-endian read of a single value from `in`.
   Throws on failure.
   */
  template<typename T>
  T read_value( ifstream &in, const string &filename )
  {
    T v;
    if( !in.read( reinterpret_cast<char *>( &v ), sizeof(T) ) )
      throw runtime_error( "GadrasShieldScatter: short read in " + filename );
    return v;
  }


  void read_floats( ifstream &in, const int count, vector<float> &out, const string &filename )
  {
    out.assign( count, 0.0f );
    if( count <= 0 )
      return;
    const std::streamsize bytes = static_cast<std::streamsize>( sizeof(float) ) * count;
    if( !in.read( reinterpret_cast<char *>( out.data() ), bytes ) )
      throw runtime_error( "GadrasShieldScatter: short read in " + filename );
  }


  void floats_to_doubles( const vector<float> &src, vector<double> &dst )
  {
    dst.assign( src.begin(), src.end() );
  }


  /** Convert a flat column-major array (as written by Fortran) of dimensions
   I*J*K*L*M into a 5-D nested vector indexed [i][j][k][l][m]. If
   `dataIsLogScaled` is false, take ln() of each element so the returned
   matrix is always log-scaled.
   */
  void column_major_to_log_5d( const bool dataIsLogScaled,
                               const vector<float> &flat,
                               const int I, const int J,
                               const int K, const int L, const int M,
                               vector< vector< vector< vector< vector<double> > > > > &out )
  {
    out.assign( I, vector< vector< vector< vector<double> > > >(
                       J, vector< vector< vector<double> > >(
                            K, vector< vector<double> >(
                                 L, vector<double>( M, 0.0 ) ) ) ) );

    for( int i = 0; i < I; ++i )
    {
      for( int j = 0; j < J; ++j )
      {
        for( int k = 0; k < K; ++k )
        {
          for( int l = 0; l < L; ++l )
          {
            for( int m = 0; m < M; ++m )
            {
              const int idx = m + M * ( l + L * ( k + K * ( j + J * i ) ) );
              const float v = flat[idx];
              out[i][j][k][l][m] = dataIsLogScaled ? static_cast<double>( v )
                                                  : std::log( static_cast<double>( v ) );
            }
          }
        }
      }
    }
  }


  /** Convert a flat column-major array of size I*J into a 2-D nested vector
   indexed [j][i] (i.e. transpose to row-major access order, matching
   GADRAS GetPointShapeFactor( atomicNumberIdx, arealDensityIdx )).
   */
  void column_major_to_2d( const vector<float> &flat,
                           const int I, const int J,
                           vector< vector<double> > &out )
  {
    out.assign( J, vector<double>( I, 0.0 ) );
    for( int i = 0; i < I; ++i )
    {
      for( int j = 0; j < J; ++j )
      {
        const int idx = ( j + J * i );
        out[j][i] = flat[idx];
      }
    }
  }


  /** Parse `filePath` into `db`. File layout matches gadras_19.5.3
   sandia.shieldscatter.db (see DatabaseParser.cpp in the GADRAS C++ source).
   */
  void parse_db_file( const string &filePath, ShieldScatterDb &db )
  {
#ifdef _WIN32
    const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16( filePath );
    ifstream in( wpath.c_str(), ios::in | ios::binary );
#else
    ifstream in( filePath.c_str(), ios::in | ios::binary );
#endif
    if( !in.is_open() )
      throw runtime_error( "GadrasShieldScatter: failed to open " + filePath );

    db.groupCount = read_value<int>( in, filePath );
    db.energyCount = read_value<int>( in, filePath );

    // Sanity bounds; the file is ~1.4 MB and a corrupt header could otherwise
    // request many gigabytes of allocation below.
    if( db.groupCount <= 0 || db.groupCount > 4096
       || db.energyCount <= 0 || db.energyCount > 1024 )
    {
      throw runtime_error( "GadrasShieldScatter: implausible counts in " + filePath );
    }

    vector<float> energies;
    read_floats( in, db.energyCount, energies, filePath );

    db.atomicNumberCount = read_value<int>( in, filePath );
    if( db.atomicNumberCount <= 0 || db.atomicNumberCount > 256 )
      throw runtime_error( "GadrasShieldScatter: implausible AN count in " + filePath );

    vector<float> atomicNumbers;
    read_floats( in, db.atomicNumberCount, atomicNumbers, filePath );

    vector<float> kEdge, kUpper, lEdge, lUpper, mEdge, mUpper;
    read_floats( in, db.atomicNumberCount, kEdge, filePath );
    read_floats( in, db.atomicNumberCount, kUpper, filePath );
    read_floats( in, db.atomicNumberCount, lEdge, filePath );
    read_floats( in, db.atomicNumberCount, lUpper, filePath );
    read_floats( in, db.atomicNumberCount, mEdge, filePath );
    read_floats( in, db.atomicNumberCount, mUpper, filePath );

    db.arealDensityCount = read_value<int>( in, filePath );
    if( db.arealDensityCount <= 0 || db.arealDensityCount > 1024 )
      throw runtime_error( "GadrasShieldScatter: implausible AD count in " + filePath );

    vector<float> arealDensities;
    read_floats( in, db.arealDensityCount, arealDensities, filePath );

    db.hydrogenMassFractionCount = read_value<int>( in, filePath );
    if( db.hydrogenMassFractionCount <= 0 || db.hydrogenMassFractionCount > 64 )
      throw runtime_error( "GadrasShieldScatter: implausible H count in " + filePath );

    vector<float> hydrogenMassFractions;
    read_floats( in, db.hydrogenMassFractionCount, hydrogenMassFractions, filePath );

    const int dataIsLogScaledFlag = read_value<int>( in, filePath );
    const bool dataIsLogScaled = ( dataIsLogScaledFlag == 1 );

    const int totalScatter = db.hydrogenMassFractionCount * db.arealDensityCount
                              * db.atomicNumberCount * db.energyCount * db.groupCount;
    const int totalShape = db.atomicNumberCount * db.arealDensityCount;

    vector<float> primaryScatterPoint, primaryScatterBasis, primaryScatterSlab;
    read_floats( in, totalScatter, primaryScatterPoint, filePath );
    read_floats( in, totalScatter, primaryScatterBasis, filePath );
    read_floats( in, totalScatter, primaryScatterSlab, filePath );

    vector<float> pointShape, basisShape, slabShape;
    read_floats( in, totalShape, pointShape, filePath );
    read_floats( in, totalShape, basisShape, filePath );
    read_floats( in, totalShape, slabShape, filePath );

    in.close();

    // Promote scalars to doubles.
    floats_to_doubles( energies, db.energies );
    floats_to_doubles( atomicNumbers, db.atomicNumbers );
    floats_to_doubles( arealDensities, db.arealDensities );
    floats_to_doubles( hydrogenMassFractions, db.hydrogenMassFractions );
    floats_to_doubles( kEdge, db.kEdgeEnergies );
    floats_to_doubles( kUpper, db.kShellRegionUpperBounds );
    floats_to_doubles( lEdge, db.lEdgeEnergies );
    floats_to_doubles( lUpper, db.lShellRegionUpperBounds );
    floats_to_doubles( mEdge, db.mEdgeEnergies );
    floats_to_doubles( mUpper, db.mShellRegionUpperBounds );

    column_major_to_log_5d( dataIsLogScaled, primaryScatterPoint,
                            db.hydrogenMassFractionCount, db.arealDensityCount,
                            db.atomicNumberCount, db.energyCount, db.groupCount,
                            db.logScatterPoint );
    column_major_to_log_5d( dataIsLogScaled, primaryScatterBasis,
                            db.hydrogenMassFractionCount, db.arealDensityCount,
                            db.atomicNumberCount, db.energyCount, db.groupCount,
                            db.logScatterBasis );
    column_major_to_log_5d( dataIsLogScaled, primaryScatterSlab,
                            db.hydrogenMassFractionCount, db.arealDensityCount,
                            db.atomicNumberCount, db.energyCount, db.groupCount,
                            db.logScatterSlab );

    column_major_to_2d( pointShape, db.arealDensityCount, db.atomicNumberCount,
                        db.pointShapeFactor );
    column_major_to_2d( basisShape, db.arealDensityCount, db.atomicNumberCount,
                        db.basisShapeFactor );
    column_major_to_2d( slabShape, db.arealDensityCount, db.atomicNumberCount,
                        db.slabShapeFactor );
  }//parse_db_file


  // ---------- Interpolation, ported from GADRAS Computer.cpp ---------------

  struct Interpolants
  {
    int Idx1 = 0;
    int Idx2 = 1;
    double Weight = 0.0;
  };


  struct ShapeWeights
  {
    double Point = 0.0;
    double Basis = 0.0;
    double Slab = 0.0;
  };


  /** Linear-search 1-D bracketing. The mesh arrays in the database all have
   fewer than 50 entries, so binary search is not worth the branch overhead.
   The output Weight is *not* clamped here; energy-domain callers clamp
   themselves to disable extrapolation.
   */
  Interpolants get_interpolants( const int count, const vector<double> &arr,
                                 const double value )
  {
    Interpolants r;
    for( r.Idx2 = 1; r.Idx2 < count; ++r.Idx2 )
    {
      if( arr[r.Idx2] > value )
        break;
    }
    r.Idx2 = std::min( count - 1, r.Idx2 );
    r.Idx1 = r.Idx2 - 1;
    r.Weight = ( value - arr[r.Idx1] ) / ( arr[r.Idx2] - arr[r.Idx1] );
    return r;
  }


  /** If the source energy lies inside one element's K/L/M-shell resonance
   region, replace that element's energy interpolants with ones that bracket
   *outside* the region, so we don't interpolate across a sharp edge.
   */
  Interpolants move_out_of_resonance( const int sourceEnergyCount,
                                      const vector<double> &sourceEnergies,
                                      const double edge, const double upperBound,
                                      const Interpolants &ergInt,
                                      const double sourceEnergy )
  {
    Interpolants result = ergInt;

    const bool within1 = ( sourceEnergies[ergInt.Idx1] > edge
                          && sourceEnergies[ergInt.Idx1] < upperBound );
    const bool within2 = ( sourceEnergies[ergInt.Idx2] > edge
                          && sourceEnergies[ergInt.Idx2] < upperBound );
    if( !within1 && !within2 )
      return result;

    if( within1 )
    {
      int ii;
      for( ii = 0; ii < sourceEnergyCount - 1; ++ii )
      {
        if( sourceEnergies[ii + 1] > edge )
          break;
      }
      result.Idx1 = ii;
    }

    if( within2 )
    {
      int ii;
      for( ii = 0; ii < sourceEnergyCount; ++ii )
      {
        if( sourceEnergies[ii] > upperBound )
          break;
      }
      result.Idx2 = ii;
    }

    result.Weight = ( sourceEnergy - sourceEnergies[result.Idx1] )
                    / ( sourceEnergies[result.Idx2] - sourceEnergies[result.Idx1] );
    result.Weight = std::max( 0.0, std::min( 1.0, result.Weight ) );
    return result;
  }


  /** Compute the geometry blending weights for a given (Z, AD) cell.
   The shape-factor scalars stored in the database are nominal "inner radius
   over outer radius" values for each geometry. shapeFactor=0 picks point,
   shapeFactor=1 picks slab, intermediate values blend.
   */
  ShapeWeights compute_shape_weights( const ShieldScatterDb &db,
                                      const double shapeFactor,
                                      const Interpolants &anInt,
                                      const Interpolants &adInt )
  {
    double pointSf = 0.0, basisSf = 0.0, slabSf = 0.0;

    if( adInt.Idx1 == -1 )
    {
      // Below-minimum AD: use AD index 0 only.
      const double w = anInt.Weight;
      pointSf = ( 1.0 - w ) * db.pointShapeFactor[anInt.Idx1][adInt.Idx2]
                + w * db.pointShapeFactor[anInt.Idx2][adInt.Idx2];
      basisSf = ( 1.0 - w ) * db.basisShapeFactor[anInt.Idx1][adInt.Idx2]
                + w * db.basisShapeFactor[anInt.Idx2][adInt.Idx2];
      slabSf = ( 1.0 - w ) * db.slabShapeFactor[anInt.Idx1][adInt.Idx2]
               + w * db.slabShapeFactor[anInt.Idx2][adInt.Idx2];
    }
    else
    {
      const double aw = anInt.Weight;
      const double dw = adInt.Weight;
      pointSf = ( 1.0 - aw ) * ( ( 1.0 - dw ) * db.pointShapeFactor[anInt.Idx1][adInt.Idx1]
                                  + dw * db.pointShapeFactor[anInt.Idx1][adInt.Idx2] )
                + aw * ( ( 1.0 - dw ) * db.pointShapeFactor[anInt.Idx2][adInt.Idx1]
                          + dw * db.pointShapeFactor[anInt.Idx2][adInt.Idx2] );
      basisSf = ( 1.0 - aw ) * ( ( 1.0 - dw ) * db.basisShapeFactor[anInt.Idx1][adInt.Idx1]
                                  + dw * db.basisShapeFactor[anInt.Idx1][adInt.Idx2] )
                + aw * ( ( 1.0 - dw ) * db.basisShapeFactor[anInt.Idx2][adInt.Idx1]
                          + dw * db.basisShapeFactor[anInt.Idx2][adInt.Idx2] );
      slabSf = ( 1.0 - aw ) * ( ( 1.0 - dw ) * db.slabShapeFactor[anInt.Idx1][adInt.Idx1]
                                 + dw * db.slabShapeFactor[anInt.Idx1][adInt.Idx2] )
               + aw * ( ( 1.0 - dw ) * db.slabShapeFactor[anInt.Idx2][adInt.Idx1]
                         + dw * db.slabShapeFactor[anInt.Idx2][adInt.Idx2] );
    }

    ShapeWeights w;
    if( shapeFactor < pointSf )
    {
      w.Point = 1.0;
    }
    else if( shapeFactor > slabSf )
    {
      w.Slab = 1.0;
    }
    else if( shapeFactor > basisSf )
    {
      w.Slab = ( shapeFactor - basisSf ) / ( slabSf - basisSf );
      w.Basis = 1.0 - w.Slab;
    }
    else
    {
      w.Basis = ( shapeFactor - pointSf ) / ( basisSf - pointSf );
      w.Point = 1.0 - w.Basis;
    }
    return w;
  }


  /** Per-call scratch buffers for the nested interpolation. The buffers are
   instance-local so multiple GadrasShieldScatter calls on different threads
   can each use their own.
   */
  struct ScratchBuffers
  {
    vector<double> shape;
    vector<double> hyd1, hyd2;
    vector<double> erg1, erg2;
    vector<double> ad1, ad2;
    vector<double> an1, an2;

    void resize( const int groupCount )
    {
      shape.resize( groupCount );
      hyd1.resize( groupCount );
      hyd2.resize( groupCount );
      erg1.resize( groupCount );
      erg2.resize( groupCount );
      ad1.resize( groupCount );
      ad2.resize( groupCount );
      an1.resize( groupCount );
      an2.resize( groupCount );
    }
  };


  void linear_interp_into( const int groupCount,
                           const vector<double> &a, const vector<double> &b,
                           const double weight, vector<double> &out )
  {
    const double rw = 1.0 - weight;
    for( int g = 0; g < groupCount; ++g )
      out[g] = rw * a[g] + weight * b[g];
  }


  void scaled_copy( const vector<double> &src, const double scalar, vector<double> &out )
  {
    for( size_t i = 0; i < src.size(); ++i )
      out[i] = scalar * src[i];
  }


  void scaled_add( const vector<double> &src, const double scalar, vector<double> &out )
  {
    for( size_t i = 0; i < src.size(); ++i )
      out[i] += scalar * src[i];
  }


  /** Pick log-scatter for one (E, Z, AD, H) cell from the requested geometry,
   scale by `weight`, and either initialize or accumulate into `out`.
   */
  void contribute_geometry( const vector<double> &cell, const double weight,
                            bool &outputInitialized, vector<double> &out )
  {
    if( outputInitialized )
    {
      scaled_add( cell, weight, out );
    }
    else
    {
      scaled_copy( cell, weight, out );
      outputInitialized = true;
    }
  }


  void interpolate_shape( const ShieldScatterDb &db,
                          const int eIdx, const int zIdx,
                          const int adIdx, const int hIdx,
                          const ShapeWeights &sw, vector<double> &out )
  {
    bool initialized = false;
    if( sw.Point > 1e-30 )
      contribute_geometry( db.logScatterPoint[hIdx][adIdx][zIdx][eIdx],
                           sw.Point, initialized, out );
    if( sw.Basis > 1e-30 )
      contribute_geometry( db.logScatterBasis[hIdx][adIdx][zIdx][eIdx],
                           sw.Basis, initialized, out );
    if( sw.Slab > 1e-30 )
      contribute_geometry( db.logScatterSlab[hIdx][adIdx][zIdx][eIdx],
                           sw.Slab, initialized, out );
  }


  void interpolate_hydrogen_shape( const ShieldScatterDb &db,
                                   const int eIdx, const int zIdx, const int adIdx,
                                   const Interpolants &hInt,
                                   const ShapeWeights &sw,
                                   ScratchBuffers &buf, vector<double> &out )
  {
    if( std::abs( hInt.Weight ) < 1e-5 )
    {
      interpolate_shape( db, eIdx, zIdx, adIdx, hInt.Idx1, sw, out );
    }
    else
    {
      interpolate_shape( db, eIdx, zIdx, adIdx, hInt.Idx1, sw, buf.hyd1 );
      interpolate_shape( db, eIdx, zIdx, adIdx, hInt.Idx2, sw, buf.hyd2 );
      linear_interp_into( db.groupCount, buf.hyd1, buf.hyd2, hInt.Weight, out );
    }
  }


  void interpolate_energy_hydrogen_shape( const ShieldScatterDb &db,
                                          const Interpolants &ergInt,
                                          const int zIdx, const int adIdx,
                                          const Interpolants &hInt,
                                          const ShapeWeights &sw,
                                          ScratchBuffers &buf, vector<double> &out )
  {
    interpolate_hydrogen_shape( db, ergInt.Idx1, zIdx, adIdx, hInt, sw, buf, buf.erg1 );
    interpolate_hydrogen_shape( db, ergInt.Idx2, zIdx, adIdx, hInt, sw, buf, buf.erg2 );
    linear_interp_into( db.groupCount, buf.erg1, buf.erg2, ergInt.Weight, out );
  }


  void interpolate_ad_energy_hydrogen_shape( const ShieldScatterDb &db,
                                             const Interpolants &ergInt,
                                             const int zIdx,
                                             const Interpolants &adInt,
                                             const Interpolants &hInt,
                                             const ShapeWeights &sw,
                                             ScratchBuffers &buf,
                                             vector<double> &out )
  {
    if( adInt.Idx1 == -1 )
    {
      // Linear extrapolation to zero scatter at AD=0; intermediate buffers
      // need to be in linear space for this single hop, then re-logged for
      // the outer Z interpolation.
      std::fill( buf.ad1.begin(), buf.ad1.end(), 0.0 );
      interpolate_energy_hydrogen_shape( db, ergInt, zIdx, adInt.Idx2,
                                         hInt, sw, buf, buf.ad2 );
      for( double &v : buf.ad2 )
      {
        if( std::abs( v ) < 1e-20 )
          continue;
        v = std::exp( v );
      }
      linear_interp_into( db.groupCount, buf.ad1, buf.ad2, adInt.Weight, out );
      for( double &v : out )
      {
        if( v < 1e-20 )
          v = -200.0;
        else
          v = std::log( v );
      }
    }
    else
    {
      interpolate_energy_hydrogen_shape( db, ergInt, zIdx, adInt.Idx1,
                                         hInt, sw, buf, buf.ad1 );
      interpolate_energy_hydrogen_shape( db, ergInt, zIdx, adInt.Idx2,
                                         hInt, sw, buf, buf.ad2 );
      linear_interp_into( db.groupCount, buf.ad1, buf.ad2, adInt.Weight, out );
    }
  }


  void interpolate_all( const ShieldScatterDb &db,
                        const Interpolants &anInt,
                        const Interpolants &an1ErgInt,
                        const Interpolants &an2ErgInt,
                        const Interpolants &adInt,
                        const Interpolants &hInt,
                        const ShapeWeights &sw,
                        ScratchBuffers &buf,
                        vector<double> &scatter )
  {
    interpolate_ad_energy_hydrogen_shape( db, an1ErgInt, anInt.Idx1, adInt, hInt, sw,
                                          buf, buf.an1 );
    interpolate_ad_energy_hydrogen_shape( db, an2ErgInt, anInt.Idx2, adInt, hInt, sw,
                                          buf, buf.an2 );
    linear_interp_into( db.groupCount, buf.an1, buf.an2, anInt.Weight, scatter );

    for( double &v : scatter )
    {
      if( std::abs( v ) < 1e-20 )
        continue;
      v = std::exp( v );
    }
  }


  void compute_scatter( const ShieldScatterDb &db,
                        double sourceEnergy,
                        double atomicNumber,
                        const double arealDensity,
                        double hydrogenADFrac,
                        const double shapeFactor,
                        vector<double> &scatter )
  {
    const int groupCount = db.groupCount;
    scatter.assign( groupCount, 0.0 );

    if( atomicNumber < 1e-20 || arealDensity < 1e-20 )
      return;

    if( hydrogenADFrac > 1e-10 )
    {
      if( std::abs( hydrogenADFrac - 1.0 ) < 1e-10 )
        atomicNumber = 1.0;
      else
        atomicNumber = ( atomicNumber - hydrogenADFrac ) / ( 1.0 - hydrogenADFrac );
    }

    const double maxAn = db.atomicNumbers[db.atomicNumberCount - 1];
    atomicNumber = std::max( 1.0, std::min( maxAn, atomicNumber ) );

    Interpolants ergInt = get_interpolants( db.energyCount, db.energies, sourceEnergy );
    ergInt.Weight = std::max( 0.0, std::min( 1.0, ergInt.Weight ) );

    const Interpolants anInt = get_interpolants( db.atomicNumberCount, db.atomicNumbers,
                                                 atomicNumber );

    Interpolants an1ErgInt = ergInt;
    Interpolants an2ErgInt = ergInt;

    // K/L/M shell resonance handling. Each shell is checked independently;
    // if the source energy lies within one element's resonance region, that
    // element's energy interpolants are pushed outside the region.
    if( sourceEnergy < db.kShellRegionUpperBounds[db.atomicNumberCount - 1] )
    {
      an1ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.kEdgeEnergies[anInt.Idx1],
                                         db.kShellRegionUpperBounds[anInt.Idx1],
                                         ergInt, sourceEnergy );
      an2ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.kEdgeEnergies[anInt.Idx2],
                                         db.kShellRegionUpperBounds[anInt.Idx2],
                                         ergInt, sourceEnergy );
    }
    if( sourceEnergy < db.lShellRegionUpperBounds[db.atomicNumberCount - 1] )
    {
      an1ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.lEdgeEnergies[anInt.Idx1],
                                         db.lShellRegionUpperBounds[anInt.Idx1],
                                         ergInt, sourceEnergy );
      an2ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.lEdgeEnergies[anInt.Idx2],
                                         db.lShellRegionUpperBounds[anInt.Idx2],
                                         ergInt, sourceEnergy );
    }
    if( sourceEnergy < db.mShellRegionUpperBounds[db.atomicNumberCount - 1] )
    {
      an1ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.mEdgeEnergies[anInt.Idx1],
                                         db.mShellRegionUpperBounds[anInt.Idx1],
                                         ergInt, sourceEnergy );
      an2ErgInt = move_out_of_resonance( db.energyCount, db.energies,
                                         db.mEdgeEnergies[anInt.Idx2],
                                         db.mShellRegionUpperBounds[anInt.Idx2],
                                         ergInt, sourceEnergy );
    }

    Interpolants adInt;
    if( arealDensity < db.arealDensities[0] )
    {
      // Special flag: linearly extrapolate to zero scatter at AD=0.
      adInt.Idx1 = -1;
      adInt.Idx2 = 0;
      adInt.Weight = arealDensity / db.arealDensities[0];
    }
    else
    {
      adInt = get_interpolants( db.arealDensityCount, db.arealDensities, arealDensity );
    }

    const Interpolants hInt = get_interpolants( db.hydrogenMassFractionCount,
                                                db.hydrogenMassFractions, hydrogenADFrac );

    const ShapeWeights sw = compute_shape_weights( db, shapeFactor, anInt, adInt );

    ScratchBuffers buf;
    buf.resize( groupCount );

    interpolate_all( db, anInt, an1ErgInt, an2ErgInt, adInt, hInt, sw, buf, scatter );

    for( double &v : scatter )
    {
      if( std::isnan( v ) )
        v = 0.0;
      v = std::max( 0.0, v );
    }
  }//compute_scatter


  // ---------- Uncollided transmission ---------------------------------------

  /** Exponential attenuation through the shield, with optional special handling
   for high hydrogen content.
   */
  float transmission_with_h( const float E, const float AN,
                              const float AD, const float AHIN )
  {
    if( AN <= 0.5f || AD <= 0.0f )
      return 1.0f;
    if( E < 1.0f )
      return 0.0f;

    const float FR = std::max( 0.0f, std::min( 1.0f, AHIN * AD - 1.0f ) );
    float AH = FR * AHIN;
    const float gcm2 = static_cast<float>( PhysicalUnits::g / PhysicalUnits::cm2 );

    if( AH > 0.0f )
    {
      if( E > 400.0f )
        AH *= std::exp( -( ( E - 400.0f ) / 1000.0f ) );

      const float ANother = ( AN * AD - AH * AD ) / ( ( 1.0f - AH ) * AD );
      const float mu = MassAttenuation::massAttenuationCoefficientFracAN( ANother, E );
      const float mu_h = MassAttenuation::massAttenuationCoefficientElement( 1, E );
      return std::exp( -( 1.0f - AH ) * AD * gcm2 * mu )
             * std::exp( -AH * AD * gcm2 * mu_h );
    }

    const float mu = MassAttenuation::massAttenuationCoefficientFracAN( AN, E );
    return std::exp( -std::max( 0.0f, AD * gcm2 * mu ) );
  }
}//anonymous namespace


// The PIMPL Database struct just wraps ShieldScatterDb so the header doesn't
// have to expose the implementation details.
struct GadrasShieldScatter::Database
{
  ShieldScatterDb db;
};


GadrasShieldScatter::GadrasShieldScatter( const std::string &datafile )
: m_db( new Database() )
{
  parse_db_file( datafile, m_db->db );
}


GadrasShieldScatter::~GadrasShieldScatter()
{
  // unique_ptr<Database> handles cleanup; needs out-of-line dtor since
  // Database is incomplete in the header.
}


int GadrasShieldScatter::groupCount() const
{
  return m_db->db.groupCount;
}


void GadrasShieldScatter::groupBounds( const double sourceEnergy,
                                       std::vector<double> &bounds ) const
{
  const int gc = m_db->db.groupCount;
  bounds.assign( gc + 1, 0.0 );
  const double scalar = sourceEnergy / static_cast<double>( gc );
  for( int i = 0; i < gc + 1; ++i )
    bounds[i] = scalar * i;
}


void GadrasShieldScatter::computeShieldScatter( const double sourceEnergy,
                                                const double atomicNumber,
                                                const double arealDensity,
                                                const double hydrogenMassFraction,
                                                const double shapeFactor,
                                                std::vector<double> &scatter ) const
{
  compute_scatter( m_db->db, sourceEnergy, atomicNumber, arealDensity,
                   hydrogenMassFraction, shapeFactor, scatter );
}


float GadrasShieldScatter::getContinuum( std::vector<float> &answer,
                                         const float sourceEnergy,
                                         const float sourceIntensity,
                                         const float atomicNumber,
                                         const float arealDensity,
                                         const float fractionAdHydrogen,
                                         const std::vector<float> &output_binning,
                                         const float shapeFactor ) const
{
  answer.assign( output_binning.size(), 0.0f );

  vector<double> scatter;
  compute_scatter( m_db->db, sourceEnergy, atomicNumber, arealDensity,
                   fractionAdHydrogen, shapeFactor, scatter );

  const int gc = m_db->db.groupCount;
  const float trans = transmission_with_h( sourceEnergy, atomicNumber,
                                           arealDensity, fractionAdHydrogen );
  const float scale = sourceIntensity * trans;

  // Build the source-energy-scaled lower-edge bin axis the database returned
  // scatter on, and rebin onto the caller's binning. Using the same
  // (one-lower-edge-per-bin) convention the legacy GadrasScatterTable uses.
  vector<float> orig_binning( gc );
  vector<float> scaled_scatter( gc );
  const double sc = sourceEnergy / static_cast<double>( gc );
  for( int g = 0; g < gc; ++g )
  {
    orig_binning[g] = static_cast<float>( sc * g );
    scaled_scatter[g] = static_cast<float>( scatter[g] ) * scale;
  }

  if( output_binning.size() >= 4 && gc >= 4 )
  {
    SpecUtils::rebin_by_lower_edge( orig_binning, scaled_scatter,
                                    output_binning, answer );
  }
  else
  {
    // rebin_by_lower_edge requires at least 4 channels on each side; for
    // very short binnings just spread the integral uniformly. This case
    // shouldn't happen in practice but guards against the throw.
    float total = 0.0f;
    for( float v : scaled_scatter )
      total += v;
    const size_t n = answer.size();
    if( n > 0 )
    {
      const float per = total / static_cast<float>( n );
      for( size_t i = 0; i < n; ++i )
        answer[i] = per;
    }
  }

  return trans * sourceIntensity;
}//GadrasShieldScatter::getContinuum
