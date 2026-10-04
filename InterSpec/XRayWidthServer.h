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

namespace SandiaDecay
{
  struct Element;
  struct Transition;
}


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

  /** Vacancy shell (destination shell of the transition), from XML `vs` attribute.
   This is the correct grouping key for shell-amplitude fitting (i.e., to group all
   radiative decays that fill the same vacancy shell).
   Empty means unknown / not provided.
   */
  std::string vacancy_shell;

  /** Natural linewidth (Lorentzian HWHM) in eV.
   If negative, then the width was not available for this line in the XML file.
   */
  double hwhm_natural_ev;

  /** Thermal Doppler broadening (Gaussian HWHM) in eV at \( T_{ref} = 295 \) K.
   Computed in code from `energy_kev` and element atomic mass.
   If negative, then the width was not available for this line in the XML file.
   */
  double hwhm_doppler_ev;

  /** Width provenance source id (from XML `src` attribute).
   Only meaningful when `hwhm_natural_ev >= 0`.
   Negative means "not present / unknown".
   */
  int source_id;

  /** xraylib "RadRate" radiative branching ratio within the vacancy shell (dimensionless),
   from XML `rr` attribute.
   Missing `rr` means unknown / not provided.
   Negative means unknown / not present.
   */
  double rad_rate;

  std::string line_label;    // X-ray transition label (e.g., "Kα1", "Lβ2")

  XRayWidthEntry();
  XRayWidthEntry( const int z, const double energy, const double natural,
                  const double doppler, const int source_id,
                  const std::string &vacancy_shell, const double rad_rate,
                  const std::string &label );

  bool has_width_data() const { return ( hwhm_natural_ev >= 0.0 ); }
  bool has_rad_rate() const { return ( rad_rate >= 0.0 ); }
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
     // Use for VoigtPlusBortel peak fitting...
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



/** Helper function to look up natural line width (Lorentzian HWHM in keV) for fluorescent x-ray transitions.

 This function provides natural line widths for fluorescent x-ray transitions from elements Z=1-98
 (hydrogen through californium) for x-rays with energies ≥10 keV. The natural line width
 represents the Lorentzian broadening due to the finite lifetime of the atomic state
 (Heisenberg uncertainty principle).

 Used for fluorescent x-rays (e.g., XRF, or x-rays a source induces in itself or its shielding).  For decay
 x-rays, `get_xray_total_width_for_decay()` finds the emitting element from the transition.

 Data is loaded from the external XML file `data/xray_widths.xml`, which contains
 comprehensive coverage of K-shell and L-shell x-ray transitions. Users can override
 the database by placing a custom xray_widths.xml in the writable data directory.

 **Usage for VoigtPlusBortel peaks:**
 Peak creators (e.g., RelActAuto) should call this function when creating peaks for fluorescent
 x-rays with the VoigtPlusBortel distribution. If a width is found (return value >= 0), set
 SkewPar0 to that value and mark it as fixed. If not found (return value < 0), leave SkewPar0
 as a fittable parameter with a reasonable starting value (~0.005 keV).

 Data sources:
 - Campbell & Papp (2001) - Widths of atomic K-N7 levels
 - Krause & Oliver (1979) - Natural widths of atomic K and L levels
 - X-ray Data Booklet (Lawrence Berkeley National Laboratory)

 @param element Pointer to the SandiaDecay::Element (must not be nullptr)
 @param energy_kev The x-ray energy in keV to find the matching transition
 @param tolerance_kev The energy tolerance for matching transitions (default 0.5 keV)

 @returns The natural linewidth Lorentzian HWHM (half-width at half-maximum) in keV if
          found, or -1.0 if the element/transition is not in the database, element is nullptr,
          or the database failed to load.

 Note: The returned value is HWHM (half-width at half-maximum), not FWHM.
 To convert to FWHM: FWHM = 2 * HWHM
 */
double get_xray_lorentzian_width( const SandiaDecay::Element *element, const double energy_kev,
                                  const double tolerance_kev = 0.5 );


/** Returns the Lorentzian linewidth (HWHM, keV) for an x-ray emitted in a decay.

 This is the natural linewidth of the x-ray of the emitting atom (the daughter, or the parent for an isomeric
 transition), as `get_xray_lorentzian_width()` gives it.  No alpha-recoil Doppler term is added: the K vacancies
 of alpha decays come almost entirely from internal conversion in the daughter, whose recoil (v/c ~ 1E-3) has
 stopped (~0.1 ps) before most levels de-excite.  LANL micro-calorimeter spectra of U-235 standards give the
 Th K-alpha1 and K-beta1 HWHM as 48-51 eV (natural 48 eV); adding the full recoil in quadrature gave 82-90 eV.

 @param transition Pointer to the SandiaDecay::Transition that produces the x-ray
                  (must not be nullptr and must contain an x-ray particle)
 @param xray_energy_kev The x-ray energy in keV (should match one of the x-rays in the transition)

 @returns Lorentzian HWHM in keV, or -1.0 if the transition is nullptr, has no x-ray at the specified energy,
          or the natural width lookup fails.
 */
double get_xray_total_width_for_decay( const SandiaDecay::Transition *transition,
                                       const double xray_energy_kev );

}//namespace XRayWidths

#endif //XRayWidthServer_h
