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
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/XmlUtils.hpp"
#include "InterSpec/XRayWidthServer.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;


namespace
{
  // Static singleton management
  std::mutex sm_mutex;
  XRayWidths::LoadStatus sm_status = XRayWidths::LoadStatus::NotLoaded;
  std::shared_ptr<const XRayWidths::XRayWidthDatabase> sm_database;

}//namespace


namespace XRayWidths
{

// ============================================================================
// XRayWidthEntry Implementation
// ============================================================================

XRayWidthEntry::XRayWidthEntry()
  : atomic_number( 0 ),
    energy_kev( 0.0 ),
    hwhm_natural_ev( 0.0 ),
    hwhm_doppler_ev( 0.0 ),
    line_label()
{
}


XRayWidthEntry::XRayWidthEntry( const int z, const double energy,
                                const double natural, const double doppler,
                                const std::string &label )
  : atomic_number( z ),
    energy_kev( energy ),
    hwhm_natural_ev( natural ),
    hwhm_doppler_ev( doppler ),
    line_label( label )
{
}


// ============================================================================
// XRayWidthDatabase Implementation
// ============================================================================

// Static member initialization
std::mutex XRayWidthDatabase::sm_mutex;
LoadStatus XRayWidthDatabase::sm_status = LoadStatus::NotLoaded;
std::shared_ptr<const XRayWidthDatabase> XRayWidthDatabase::sm_database;


std::shared_ptr<const XRayWidthDatabase> XRayWidthDatabase::instance()
{
  std::lock_guard<std::mutex> lock( sm_mutex );

  if( sm_status == LoadStatus::Loaded )
  {
    assert( sm_database );
    return sm_database;
  }

  if( sm_status == LoadStatus::FailedToLoad )
    return nullptr;

  // First time initialization
  assert( sm_status == LoadStatus::NotLoaded );
  assert( !sm_database );

  try
  {
    std::shared_ptr<XRayWidthDatabase> db( new XRayWidthDatabase() );
    db->init();

    sm_database = db;
    sm_status = LoadStatus::Loaded;

    return sm_database;
  }
  catch( std::exception &e )
  {
    cerr << "XRayWidthDatabase::instance() - Failed to load database: "
         << e.what() << endl;
    sm_status = LoadStatus::FailedToLoad;
    return nullptr;
  }
}//instance()


void XRayWidthDatabase::remove_global_instance()
{
  std::lock_guard<std::mutex> lock( sm_mutex );

  if( sm_database )
  {
    assert( sm_status == LoadStatus::Loaded );
    sm_database.reset();
    sm_status = LoadStatus::NotLoaded;
  }
  else
  {
    assert( sm_status != LoadStatus::Loaded );
  }
}//remove_global_instance()


LoadStatus XRayWidthDatabase::status()
{
  std::lock_guard<std::mutex> lock( sm_mutex );
  return sm_status;
}//status()


XRayWidthDatabase::XRayWidthDatabase()
  : m_widths_by_element()
{
}


void XRayWidthDatabase::init()
{
  m_widths_by_element.clear();

  const string filename = "xray_widths.xml";

  ifstream infile;

  // First try user writable data directory
#if( BUILD_AS_ELECTRON_APP || IOS || ANDROID || BUILD_AS_OSX_APP || BUILD_AS_LOCAL_SERVER || BUILD_AS_WX_WIDGETS_APP )
  if( !infile.is_open() )
  {
    try
    {
      const string data_dir = InterSpec::writableDataDirectory();
      const string file_path = SpecUtils::append_path( data_dir, filename );
#ifdef _WIN32
      infile.open( SpecUtils::convert_from_utf8_to_utf16(file_path).c_str(), ios::in | ios::binary );
#else
      infile.open( file_path.c_str(), ios::in | ios::binary );
#endif
    }
    catch( std::exception & )
    {
      // InterSpec::writableDataDirectory throws exception if it hasn't been set
    }
  }
#endif

  // Try static data directory
  if( !infile.is_open() )
  {
    const string data_dir = InterSpec::staticDataDirectory();
    const string file_path = SpecUtils::append_path( data_dir, filename );
#ifdef _WIN32
    infile.open( SpecUtils::convert_from_utf8_to_utf16( file_path ).c_str(), ios::in | ios::binary );
#else
    infile.open( file_path.c_str(), ios::in | ios::binary );
#endif
  }

  // Finally try current working directory
  if( !infile.is_open() )
  {
#ifdef _WIN32
    infile.open( SpecUtils::convert_from_utf8_to_utf16( filename ).c_str(), ios::in | ios::binary );
#else
    infile.open( filename.c_str(), ios::in | ios::binary );
#endif
  }

  if( !infile.is_open() )
    throw runtime_error( "XRayWidthDatabase::init: could not open '" + filename + "'" );

  rapidxml::file<char> input_data( infile );
  rapidxml::xml_document<char> doc;

  try
  {
    // Parse non-destructively
    doc.parse<rapidxml::parse_trim_whitespace | rapidxml::parse_non_destructive>( input_data.data() );
  }
  catch( rapidxml::parse_error &e )
  {
    const string msg = "XRayWidthDatabase::init: Error parsing '" + filename + "': " + string(e.what());
    throw runtime_error( msg );
  }

  const rapidxml::xml_node<char> *base_node = doc.first_node( "XRayWidths" );
  if( !base_node )
    throw runtime_error( "XRayWidthDatabase::init: Could not find root <XRayWidths> node in '" + filename + "'" );

  // Parse all <Element> nodes
  XML_FOREACH_CHILD( element_node, base_node, "Element" )
  {
    // Parse Z attribute
    const rapidxml::xml_attribute<char> *z_attr = element_node->first_attribute( "Z", 1 );
    if( !z_attr )
    {
      cerr << "Warning: XRayWidthDatabase - <Element> node missing Z attribute, skipping" << endl;
      continue;
    }

    int z = 0;
    if( !SpecUtils::parse_int( z_attr->value(), z_attr->value_size(), z ) || z < 1 || z > 98 )
    {
      cerr << "Warning: XRayWidthDatabase - invalid Z=" << z_attr->value() << ", skipping" << endl;
      continue;
    }

    // Parse all <XRay> child nodes
    XML_FOREACH_CHILD( xray_node, element_node, "XRay" )
    {
      try
      {
        // Parse energy attribute
        const rapidxml::xml_attribute<char> *energy_attr = xray_node->first_attribute( "energy", 6 );
        if( !energy_attr )
        {
          cerr << "Warning: XRayWidthDatabase - <XRay> missing energy attribute for Z=" << z << ", skipping" << endl;
          continue;
        }

        double energy_kev = 0.0;
        if( !SpecUtils::parse_double( energy_attr->value(), energy_attr->value_size(), energy_kev ) )
        {
          cerr << "Warning: XRayWidthDatabase - invalid energy for Z=" << z << ", skipping" << endl;
          continue;
        }

        // Validate energy threshold (≥10 keV)
        if( energy_kev < 10.0 )
        {
          cerr << "Warning: XRayWidthDatabase - energy " << energy_kev
               << " keV below 10 keV threshold for Z=" << z << ", skipping" << endl;
          continue;
        }

        // Parse hwhm_natural attribute
        const rapidxml::xml_attribute<char> *natural_attr = xray_node->first_attribute( "hwhm_natural", 12 );
        if( !natural_attr )
        {
          cerr << "Warning: XRayWidthDatabase - <XRay> missing hwhm_natural for Z=" << z
               << " at " << energy_kev << " keV, skipping" << endl;
          continue;
        }

        double hwhm_natural_ev = 0.0;
        if( !SpecUtils::parse_double( natural_attr->value(), natural_attr->value_size(), hwhm_natural_ev ) )
        {
          cerr << "Warning: XRayWidthDatabase - invalid hwhm_natural for Z=" << z << ", skipping" << endl;
          continue;
        }

        // Validate natural width range (0.1-200 eV)
        if( hwhm_natural_ev < 0.1 || hwhm_natural_ev > 200.0 )
        {
          cerr << "Warning: XRayWidthDatabase - unreasonable natural width " << hwhm_natural_ev
               << " eV for Z=" << z << " at " << energy_kev << " keV, skipping" << endl;
          continue;
        }

        // Parse hwhm_doppler attribute
        const rapidxml::xml_attribute<char> *doppler_attr = xray_node->first_attribute( "hwhm_doppler", 12 );
        if( !doppler_attr )
        {
          cerr << "Warning: XRayWidthDatabase - <XRay> missing hwhm_doppler for Z=" << z
               << " at " << energy_kev << " keV, skipping" << endl;
          continue;
        }

        double hwhm_doppler_ev = 0.0;
        if( !SpecUtils::parse_double( doppler_attr->value(), doppler_attr->value_size(), hwhm_doppler_ev ) )
        {
          cerr << "Warning: XRayWidthDatabase - invalid hwhm_doppler for Z=" << z << ", skipping" << endl;
          continue;
        }

        // Validate Doppler width range (0.01-5 eV at 295K)
        if( hwhm_doppler_ev < 0.01 || hwhm_doppler_ev > 5.0 )
        {
          cerr << "Warning: XRayWidthDatabase - unreasonable Doppler width " << hwhm_doppler_ev
               << " eV for Z=" << z << " at " << energy_kev << " keV, skipping" << endl;
          continue;
        }

        // Parse optional line label
        string line_label;
        const rapidxml::xml_attribute<char> *line_attr = xray_node->first_attribute( "line", 4 );
        if( line_attr )
          line_label = SpecUtils::xml_value_str( line_attr );

        // Create entry and add to database
        XRayWidthEntry entry( z, energy_kev, hwhm_natural_ev, hwhm_doppler_ev, line_label );
        m_widths_by_element[z].push_back( entry );

      }
      catch( std::exception &e )
      {
        cerr << "Warning: XRayWidthDatabase - exception parsing <XRay> for Z=" << z
             << ": " << e.what() << ", skipping entry" << endl;
      }
    }//foreach <XRay> node
  }//foreach <Element> node

  // Report loading statistics
  size_t total_entries = 0;
  for( const auto &elem_pair : m_widths_by_element )
    total_entries += elem_pair.second.size();

  cout << "XRayWidthDatabase loaded " << total_entries << " x-ray width entries for "
       << m_widths_by_element.size() << " elements" << endl;

}//init()


const XRayWidthEntry * XRayWidthDatabase::find_best_match( const int z,
                                                            const double energy_kev,
                                                            const double tolerance_kev ) const
{
  const auto elem_iter = m_widths_by_element.find( z );
  if( elem_iter == m_widths_by_element.end() )
    return nullptr;

  const vector<XRayWidthEntry> &entries = elem_iter->second;
  if( entries.empty() )
    return nullptr;

  const XRayWidthEntry *best_match = nullptr;
  double best_diff = tolerance_kev;

  for( const XRayWidthEntry &entry : entries )
  {
    const double energy_diff = fabs( entry.energy_kev - energy_kev );
    if( energy_diff < best_diff )
    {
      best_diff = energy_diff;
      best_match = &entry;
    }
  }

  return best_match;
}//find_best_match()


double XRayWidthDatabase::get_natural_width_hwhm_kev( const int z,
                                                       const double energy_kev,
                                                       const double tolerance_kev ) const
{
  const XRayWidthEntry *entry = find_best_match( z, energy_kev, tolerance_kev );
  if( !entry )
    return -1.0;

  // Convert eV to keV
  return entry->hwhm_natural_ev / 1000.0;
}//get_natural_width_hwhm_kev()


double XRayWidthDatabase::get_doppler_width_hwhm_kev( const int z,
                                                       const double energy_kev,
                                                       const double tolerance_kev,
                                                       const double temperature_k ) const
{
  const XRayWidthEntry *entry = find_best_match( z, energy_kev, tolerance_kev );
  if( !entry )
    return -1.0;

  // Scale Doppler width with temperature: HWHM ∝ sqrt(T)
  // Database stores width at 295K
  const double temperature_scale = sqrt( temperature_k / 295.0 );
  const double hwhm_doppler_ev = entry->hwhm_doppler_ev * temperature_scale;

  // Convert eV to keV
  return hwhm_doppler_ev / 1000.0;
}//get_doppler_width_hwhm_kev()


double XRayWidthDatabase::get_total_width_hwhm_kev( const int z,
                                                     const double energy_kev,
                                                     const double tolerance_kev,
                                                     const double temperature_k ) const
{
  const XRayWidthEntry *entry = find_best_match( z, energy_kev, tolerance_kev );
  if( !entry )
    return -1.0;

  // Scale Doppler width with temperature
  const double temperature_scale = sqrt( temperature_k / 295.0 );
  const double hwhm_doppler_ev = entry->hwhm_doppler_ev * temperature_scale;

  // Combine in quadrature: total² = natural² + doppler²
  const double hwhm_total_ev = sqrt( entry->hwhm_natural_ev * entry->hwhm_natural_ev +
                                      hwhm_doppler_ev * hwhm_doppler_ev );

  // Convert eV to keV
  return hwhm_total_ev / 1000.0;
}//get_total_width_hwhm_kev()


std::vector<XRayWidthEntry> XRayWidthDatabase::get_all_widths_for_element( const int z ) const
{
  const auto elem_iter = m_widths_by_element.find( z );
  if( elem_iter == m_widths_by_element.end() )
    return std::vector<XRayWidthEntry>();

  return elem_iter->second;
}//get_all_widths_for_element()


size_t XRayWidthDatabase::num_entries() const
{
  size_t total = 0;
  for( const auto &elem_pair : m_widths_by_element )
    total += elem_pair.second.size();

  return total;
}//num_entries()


double compute_alpha_recoil_doppler_hwhm( const std::string &nuclide_symbol,
                                           const double xray_energy_kev,
                                           const double alpha_energy_kev )
{
  // Get SandiaDecay database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    cerr << "compute_alpha_recoil_doppler_hwhm: SandiaDecay database not available" << endl;
    return -1.0;
  }

  // Look up parent nuclide
  const SandiaDecay::Nuclide * const parent = db->nuclide( nuclide_symbol );
  if( !parent )
  {
    cerr << "compute_alpha_recoil_doppler_hwhm: Nuclide '" << nuclide_symbol << "' not found" << endl;
    return -1.0;
  }

  // Find all alpha decay transitions and collect alpha particles with their intensities
  struct AlphaInfo
  {
    double energy_kev;
    double intensity;  // Relative intensity (branch ratio × alpha intensity)
    const SandiaDecay::Nuclide *daughter;
  };

  std::vector<AlphaInfo> alpha_list;
  double total_alpha_intensity = 0.0;

  for( const SandiaDecay::Transition *transition : parent->decaysToChildren )
  {
    if( !transition )
      continue;

    // Check if this is an alpha decay
    if( transition->mode != SandiaDecay::AlphaDecay &&
        transition->mode != SandiaDecay::BetaAndAlphaDecay &&
        transition->mode != SandiaDecay::ElectronCaptureAndAlphaDecay &&
        transition->mode != SandiaDecay::BetaPlusAndAlphaDecay )
      continue;

    const SandiaDecay::Nuclide * const daughter = transition->child;
    if( !daughter )
      continue;

    const double branch_ratio = transition->branchRatio;

    // Find alpha particles in this transition
    for( const SandiaDecay::RadParticle &particle : transition->products )
    {
      if( particle.type != SandiaDecay::AlphaParticle )
        continue;

      const double this_alpha_energy = particle.energy;
      const double this_alpha_intensity = particle.intensity;

      // If user specified alpha energy, only include matching alphas (within 10 keV tolerance)
      if( alpha_energy_kev > 0.0 )
      {
        if( fabs( this_alpha_energy - alpha_energy_kev ) > 10.0 )
          continue;
      }

      // Add this alpha to the list
      AlphaInfo info;
      info.energy_kev = this_alpha_energy;
      info.intensity = branch_ratio * this_alpha_intensity;
      info.daughter = daughter;

      alpha_list.push_back( info );
      total_alpha_intensity += info.intensity;
    }
  }

  if( alpha_list.empty() || total_alpha_intensity <= 0.0 )
  {
    cerr << "compute_alpha_recoil_doppler_hwhm: No alpha decay found for " << nuclide_symbol << endl;
    return -1.0;
  }

  // Compute intensity-weighted Doppler broadening
  // Since different alpha energies give different recoil velocities, we compute the
  // Doppler width for each alpha and then combine them weighted by intensity

  const double alpha_mass_amu = 4.002603;  // Alpha particle mass in amu
  const double amu_to_kev = 931494.0;  // keV per amu (1 amu·c² = 931.494 MeV)
  const double sqrt_2ln2 = 1.177410022515;  // sqrt(2×ln(2))

  double weighted_hwhm_sum = 0.0;

  for( const AlphaInfo &alpha : alpha_list )
  {
    const double daughter_mass_amu = alpha.daughter->massNumber;

    // Compute recoil energy from momentum conservation
    // E_recoil = E_alpha × (m_alpha / m_daughter)
    const double recoil_energy_kev = alpha.energy_kev * (alpha_mass_amu / daughter_mass_amu);

    // Compute recoil velocity: E = 0.5 × m × v²  →  v = sqrt(2 × E / m)
    const double recoil_velocity_over_c = sqrt( 2.0 * recoil_energy_kev / (daughter_mass_amu * amu_to_kev) );

    // Compute Doppler broadening HWHM for this alpha
    // For random recoil directions (isotropic): HWHM ≈ E × (v/c) / sqrt(2×ln2)
    const double this_hwhm_kev = xray_energy_kev * recoil_velocity_over_c / sqrt_2ln2;

    // Weight by intensity
    weighted_hwhm_sum += this_hwhm_kev * alpha.intensity;
  }

  const double hwhm_doppler_kev = weighted_hwhm_sum / total_alpha_intensity;

  #if( PERFORM_DEVELOPER_CHECKS )
    // Sanity check: for U-238 Kα (98.4 keV), expect ~70-100 eV HWHM
    if( hwhm_doppler_kev < 0.001 || hwhm_doppler_kev > 0.300 )
    {
      cout << "compute_alpha_recoil_doppler_hwhm: Warning - unusual recoil Doppler width "
           << (hwhm_doppler_kev * 1000.0) << " eV for " << nuclide_symbol
           << " x-ray at " << xray_energy_kev << " keV" << endl;
      cout << "  Found " << alpha_list.size() << " alpha particle(s):" << endl;
      for( const AlphaInfo &alpha : alpha_list )
      {
        cout << "    E_alpha = " << alpha.energy_kev << " keV, intensity = "
             << (alpha.intensity * 100.0 / total_alpha_intensity) << "%" << endl;
      }
    }
  #endif

  return hwhm_doppler_kev;
}//compute_alpha_recoil_doppler_hwhm()


}//namespace XRayWidths
