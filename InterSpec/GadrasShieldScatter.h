#ifndef GadrasShieldScatter_h
#define GadrasShieldScatter_h
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

#include <memory>
#include <string>
#include <vector>


/** Computes the scattered continuum reaching a detector from a shielded
 point gamma source, using GADRAS's modern shielding-scatter table
 (`sandia.shieldscatter.db`).

 Parallel to GadrasScatterTable, which uses the legacy `GadrasContinuum.lib`
 file. The two classes are intentionally interchangeable for InterSpec's
 dose-calc / activity-calc usage; the new model interpolates a richer
 5-D log-scaled table (energy, Z, areal density, hydrogen mass fraction,
 output group) and supports geometry blending between point/basis/slab
 source-in-shield configurations.
 */
class GadrasShieldScatter
{
public:
  /** Loads `sandia.shieldscatter.db` from `datafile`.
   Throws std::runtime_error if the file cannot be opened or is malformed.
   */
  explicit GadrasShieldScatter( const std::string &datafile );

  ~GadrasShieldScatter();

  /** Drop-in replacement for GadrasScatterTable::getContinuum.

   Computes the per-uncollided-leakage-photon scatter spectrum on the
   database's native group structure, scales it by `sourceIntensity`
   times the uncollided transmission fraction, and rebins it to
   `binning`. Writes results into `answer` (resized to binning.size()).

   \param answer  Output bin counts; resized and zeroed.
   \param sourceEnergy  Incident gamma energy, in keV.
   \param sourceIntensity  Number of source gammas before shielding.
   \param atomicNumber  Effective atomic number of the shielding (1..max in db).
   \param arealDensity  Areal density of shielding, g/cm^2.
   \param fractionAdHydrogen  Mass fraction of hydrogen in the shield (0..1).
   \param binning  Lower-edge energies (keV) for the output bins.
   \param shapeFactor  0 = point source in solid sphere, 1 = slab geometry.
                       Defaults to 0 to match GadrasScatterTable's geometry.
   \returns Uncollided source intensity reaching the detector
            (i.e. sourceIntensity * uncollided_transmission).
   */
  float getContinuum( std::vector<float> &answer,
                      const float sourceEnergy,
                      const float sourceIntensity,
                      const float atomicNumber,
                      const float arealDensity,
                      const float fractionAdHydrogen,
                      const std::vector<float> &binning,
                      const float shapeFactor = 0.0f ) const;

  /** Direct interface mirroring the GADRAS C function. Returns scatter
   intensities PER UNCOLLIDED LEAKAGE PHOTON on the database's native group
   structure. No transmission scaling, no rebinning. `scatter` is resized
   to groupCount().
   */
  void computeShieldScatter( const double sourceEnergy,
                             const double atomicNumber,
                             const double arealDensity,
                             const double hydrogenMassFraction,
                             const double shapeFactor,
                             std::vector<double> &scatter ) const;

  /** Number of scatter output groups in the database. */
  int groupCount() const;

  /** Returns the source-energy-dependent group bin edges, sized
   groupCount()+1, equally spaced from 0 to `sourceEnergy`.
   */
  void groupBounds( const double sourceEnergy,
                    std::vector<double> &bounds ) const;

private:
  // Internal database; defined in the .cpp file.
  struct Database;
  std::unique_ptr<Database> m_db;
};//class GadrasShieldScatter

#endif  //GadrasShieldScatter_h
