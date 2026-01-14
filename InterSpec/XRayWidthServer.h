#ifndef XRayWidthServer_h
#define XRayWidthServer_h
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

#include <map>
#include <mutex>
#include <memory>
#include <string>
#include <vector>


/** Namespace for x-ray natural linewidth and Doppler broadening data.

 This server provides access to natural linewidths (Lorentzian HWHM) and
 thermal Doppler broadening widths (Gaussian HWHM) for x-ray transitions
 from elements Z=1-98 (through Californium).

 Data is loaded from the external XML file `data/xray_widths.xml` following
 InterSpec's established singleton server pattern (similar to MoreNuclideInfo
 and DecayDataBaseServer).

 Data sources:
 - Campbell & Papp (2001) - Widths of atomic K-N7 levels
 - Krause & Oliver (1979) - Natural widths of K and L levels
 - LBNL X-ray Data Booklet
 */
namespace XRayWidths
{

enum class LoadStatus
{
  NotLoaded,      // Database not yet initialized
  FailedToLoad,   // Failed to load xray_widths.xml
  Loaded          // Successfully loaded
};//enum class LoadStatus


/** Structure to hold x-ray width data for a single transition.
 */
struct XRayWidthEntry
{
  int atomic_number;         // Atomic number (Z)
  double energy_kev;         // X-ray energy in keV
  double hwhm_natural_ev;    // Natural linewidth (Lorentzian HWHM) in eV
  double hwhm_doppler_ev;    // Doppler width at 295K (Gaussian HWHM) in eV
  std::string line_label;    // X-ray transition label (e.g., "Kα1", "Lβ2")

  XRayWidthEntry();
  XRayWidthEntry( const int z, const double energy, const double natural,
                  const double doppler, const std::string &label );
};//struct XRayWidthEntry


/** X-ray width database singleton class.

 Loads and provides access to natural linewidths and Doppler broadening data
 for x-ray transitions. Thread-safe singleton with lazy initialization.

 Example usage:
 @code
   // Get natural width for U Kα1 at 98.4 keV
   const auto db = XRayWidthDatabase::instance();
   if( db )
   {
     const double width_kev = db->get_natural_width_hwhm_kev( 92, 98.4 );
     // Use for VoigtWithExpTail peak fitting...
   }
 @endcode
 */
class XRayWidthDatabase
{
public:
  /** Returns singleton instance of the database.

   First call will load the database from xray_widths.xml. Subsequent calls
   return the cached instance.

   @returns Shared pointer to database, or nullptr if loading failed
   */
  static std::shared_ptr<const XRayWidthDatabase> instance();


  /** Removes the global instance (for testing or reloading).
   */
  static void remove_global_instance();


  /** Returns current load status of the database.

   @returns LoadStatus indicating whether database is loaded, failed, or not yet initialized
   */
  static LoadStatus status();


  /** Returns natural linewidth (Lorentzian HWHM) for x-ray transition.

   Use this function for fluorescence x-rays where Doppler broadening is not
   relevant (e.g., XRF analysis, electron/synchrotron-induced x-rays).

   @param z Atomic number (1-98)
   @param energy_kev X-ray energy in keV
   @param tolerance_kev Energy matching tolerance (default 0.5 keV)

   @returns Natural width HWHM in keV, or -1.0 if no match found
   */
  double get_natural_width_hwhm_kev( const int z, const double energy_kev,
                                      const double tolerance_kev = 0.5 ) const;


  /** Returns Doppler broadening width (Gaussian HWHM) for x-ray transition.

   Doppler broadening arises from thermal motion of the nucleus emitting the
   x-ray. Width scales with temperature as HWHM ∝ √T.

   For decay x-rays at room temperature (295K), Doppler width is typically
   0.3-0.7 eV for heavy elements and high-energy K-shell x-rays.

   @param z Atomic number (1-98)
   @param energy_kev X-ray energy in keV
   @param tolerance_kev Energy matching tolerance (default 0.5 keV)
   @param temperature_k Sample temperature in Kelvin (default 295K = room temp)

   @returns Doppler width HWHM in keV at specified temperature, or -1.0 if no match
   */
  double get_doppler_width_hwhm_kev( const int z, const double energy_kev,
                                      const double tolerance_kev = 0.5,
                                      const double temperature_k = 295.0 ) const;


  /** Returns total effective width for decay x-rays (natural + Doppler).

   For decay x-rays from radioactive sources, the total width combines natural
   and Doppler components in quadrature:
     Total HWHM = sqrt(natural² + doppler²)

   Use this function for characteristic x-rays from radioactive decay (U-238,
   Pu-239, etc.) measured at room temperature.

   @param z Atomic number (1-98)
   @param energy_kev X-ray energy in keV
   @param tolerance_kev Energy matching tolerance (default 0.5 keV)
   @param temperature_k Sample temperature in Kelvin (default 295K = room temp)

   @returns Total width HWHM in keV, or -1.0 if no match found
   */
  double get_total_width_hwhm_kev( const int z, const double energy_kev,
                                    const double tolerance_kev = 0.5,
                                    const double temperature_k = 295.0 ) const;


  /** Returns all x-ray width entries for a given element.

   Useful for debugging, validation, or displaying available data.

   @param z Atomic number (1-98)
   @returns Vector of all width entries for this element (may be empty)
   */
  std::vector<XRayWidthEntry> get_all_widths_for_element( const int z ) const;


  /** Returns total number of x-ray width entries in database.

   @returns Total entry count (should be ~700-1000 for complete database)
   */
  size_t num_entries() const;


private:
  /** Private constructor - use instance() to access singleton.
   */
  XRayWidthDatabase();


  /** Initializes database by loading xray_widths.xml.

   Searches for file in:
   1. User writable directory
   2. Static data directory

   If file not found or parsing fails, marks status as FailedToLoad.
   */
  void init();


  /** Parses xray_widths.xml file and populates database.

   @param filepath Full path to xray_widths.xml
   @throws std::runtime_error on parsing errors
   */
  void parse_xml( const std::string &filepath );


  /** Finds best matching x-ray entry for given Z and energy.

   @param z Atomic number
   @param energy_kev X-ray energy in keV
   @param tolerance_kev Energy matching tolerance

   @returns Pointer to best matching entry, or nullptr if no match within tolerance
   */
  const XRayWidthEntry * find_best_match( const int z, const double energy_kev,
                                           const double tolerance_kev ) const;


  // Storage: map from atomic number to vector of width entries
  std::map<int, std::vector<XRayWidthEntry>> m_widths_by_element;

  // Global singleton management (static members)
  static std::mutex sm_mutex;
  static LoadStatus sm_status;
  static std::shared_ptr<const XRayWidthDatabase> sm_database;

};//class XRayWidthDatabase


/** Computes nuclear recoil Doppler broadening HWHM for x-rays following alpha decay.

 For alpha-emitting nuclides, characteristic x-rays can be emitted by the recoiling
 daughter nucleus. The recoil velocity causes Doppler broadening of the x-ray line.

 Recoil energy: E_recoil = E_alpha × (m_alpha / m_daughter)
 Recoil velocity: v_recoil = sqrt(2 × E_recoil / m_daughter)
 Doppler HWHM: HWHM = E_xray × (v_recoil / c) / sqrt(2×ln2)

 For U-238 → Th-234 alpha decay (E_alpha ≈ 4.2 MeV):
   - Recoil energy ≈ 72 keV
   - Recoil velocity ≈ 350 km/s
   - Doppler HWHM for U Kα (98.4 keV) ≈ 70-100 eV

 @param nuclide_symbol Parent nuclide symbol (e.g., "U238", "Pu239")
 @param xray_energy_kev X-ray energy in keV
 @param alpha_energy_kev Alpha particle energy in keV (optional - will use highest energy alpha if not specified)

 @returns Recoil Doppler HWHM in keV, or -1.0 if computation fails (nuclide not found, no alpha decay, etc.)

 Note: This computes recoil from alpha decay only. For x-rays from internal conversion
       or thermal equilibrium, this contribution should not be included.
 */
double compute_alpha_recoil_doppler_hwhm( const std::string &nuclide_symbol,
                                           const double xray_energy_kev,
                                           const double alpha_energy_kev = -1.0 );

}//namespace XRayWidths

#endif //XRayWidthServer_h
