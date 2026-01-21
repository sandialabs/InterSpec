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

  constexpr double sm_doppler_reference_temperature_k = 295.0;

  /** Doppler HWHM (Gaussian) in eV at 295 K.
   Formula: hwhm_eV = 0.0855 * ( E_keV / sqrt( M_amu ) )
   */
  double doppler_hwhm_ev_at_reference_temp( const double energy_kev, const double atomic_mass_amu )
  {
    if( energy_kev <= 0.0 || atomic_mass_amu <= 0.0 )
      return -1.0;

    return 0.0855 * ( energy_kev / sqrt( atomic_mass_amu ) );
  }

}//namespace


namespace XRayWidths
{

// ============================================================================
// XRayWidthEntry Implementation
// ============================================================================

XRayWidthEntry::XRayWidthEntry()
  : atomic_number( 0 ),
    energy_kev( 0.0 ),
    vacancy_shell(),
    hwhm_natural_ev( -1.0 ),
    hwhm_doppler_ev( -1.0 ),
    source_id( -1 ),
    rad_rate( -1.0 ),
    line_label()
{
}


XRayWidthEntry::XRayWidthEntry( const int z, const double energy,
                                const double natural, const double doppler,
                                const int src_id,
                                const std::string &vacancy_shell_,
                                const double rad_rate_,
                                const std::string &label )
  : atomic_number( z ),
    energy_kev( energy ),
    vacancy_shell( vacancy_shell_ ),
    hwhm_natural_ev( natural ),
    hwhm_doppler_ev( doppler ),
    source_id( src_id ),
    rad_rate( rad_rate_ ),
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

  // Resolve energy threshold (default 10 keV if missing)
  double energy_threshold_kev = 10.0;
  const rapidxml::xml_node<char> *meta_node = base_node->first_node( "Meta", 4 );
  if( meta_node )
  {
    const rapidxml::xml_node<char> *threshold_node = meta_node->first_node( "EnergyThreshold", 15 );
    if( threshold_node && threshold_node->value() && threshold_node->value_size() )
    {
      double parsed_threshold = 0.0;
      if( SpecUtils::parse_double( threshold_node->value(), threshold_node->value_size(), parsed_threshold )
          && parsed_threshold > 0.0 )
      {
        energy_threshold_kev = parsed_threshold;
      }
    }
  }

  // We compute Doppler widths in code using SandiaDecay element atomic masses.
  // If SandiaDecay is not available, we still load natural widths, but Doppler
  // and total widths will be unavailable.
  const SandiaDecay::SandiaDecayDataBase *decay_db = nullptr;
  try
  {
    decay_db = DecayDataBaseServer::database();
  }
  catch( std::exception &e )
  {
    cerr << "Warning: XRayWidthDatabase - SandiaDecay database not available ("
         << e.what() << "); Doppler widths will be unavailable" << endl;
  }

  // Parse all <Element> nodes.
  //
  // XML fields of interest:
  // - `vs`: vacancy shell (destination shell of the transition); correct grouping key for shell-amplitude fitting.
  // - `rr`: xraylib "RadRate", a radiative branching ratio among radiative decays within that vacancy shell (dimensionless).
  // - `hwhm`: Lorentzian HWHM (eV). Missing `hwhm` means width unavailable.
  // - Missing `rr` means unknown / not provided.
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

    // Element atomic mass (amu) for Doppler computation
    const SandiaDecay::Element *sd_element = decay_db ? decay_db->element( z ) : nullptr;
    if( !sd_element )
    {
      cerr << "Warning: XRayWidthDatabase - SandiaDecay element missing for Z=" << z
           << ", Doppler widths will be unavailable" << endl;
    }

    const double atomic_mass_amu = sd_element ? sd_element->atomicMass() : -1.0;

    // Parse all <XRay> child nodes
    XML_FOREACH_CHILD( xray_node, element_node, "XRay" )
    {
      try
      {
        // Parse energy attribute (new: "e", old: "energy")
        const rapidxml::xml_attribute<char> *energy_attr = xray_node->first_attribute( "e", 1 );
        if( !energy_attr )
          energy_attr = xray_node->first_attribute( "energy", 6 );
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

        // Validate energy threshold (≥ EnergyThreshold keV)
        if( energy_kev < energy_threshold_kev )
        {
          cerr << "Warning: XRayWidthDatabase - energy " << energy_kev
               << " keV below " << energy_threshold_kev << " keV threshold for Z=" << z << ", skipping" << endl;
          continue;
        }

        // Parse optional width attribute (new: "hwhm", old: "hwhm_natural")
        const rapidxml::xml_attribute<char> *hwhm_attr = xray_node->first_attribute( "hwhm", 4 );
        if( !hwhm_attr )
          hwhm_attr = xray_node->first_attribute( "hwhm_natural", 12 );

        double hwhm_natural_ev = -1.0;
        double hwhm_doppler_ev = -1.0;
        int src_id = -1;
        double rr = -1.0;

        // Parse optional vacancy shell (destination shell), `vs`
        std::string vacancy_shell;
        const rapidxml::xml_attribute<char> *vs_attr = xray_node->first_attribute( "vs", 2 );
        if( vs_attr )
          vacancy_shell = SpecUtils::xml_value_str( vs_attr );

        // Parse optional radiative branching ratio within that vacancy shell, `rr`
        const rapidxml::xml_attribute<char> *rr_attr = xray_node->first_attribute( "rr", 2 );
        if( rr_attr )
        {
          double parsed_rr = -1.0;
          if( SpecUtils::parse_double( rr_attr->value(), rr_attr->value_size(), parsed_rr ) )
          {
            if( parsed_rr >= 0.0 )
              rr = parsed_rr;
          }
        }

        if( hwhm_attr )
        {
          double parsed_hwhm = 0.0;
          if( !SpecUtils::parse_double( hwhm_attr->value(), hwhm_attr->value_size(), parsed_hwhm ) )
          {
            cerr << "Warning: XRayWidthDatabase - invalid hwhm for Z=" << z << " at "
                 << energy_kev << " keV, treating width as unavailable" << endl;
          }
          else if( parsed_hwhm <= 0.0 || parsed_hwhm > 2000.0 )
          {
            cerr << "Warning: XRayWidthDatabase - unreasonable hwhm " << parsed_hwhm
                 << " eV for Z=" << z << " at " << energy_kev
                 << " keV, treating width as unavailable" << endl;
          }
          else
          {
            hwhm_natural_ev = parsed_hwhm;

            if( atomic_mass_amu > 0.0 )
              hwhm_doppler_ev = doppler_hwhm_ev_at_reference_temp( energy_kev, atomic_mass_amu );

            // Parse optional provenance source id (only meaningful when width exists)
            const rapidxml::xml_attribute<char> *src_attr = xray_node->first_attribute( "src", 3 );
            if( src_attr )
            {
              int parsed_src_id = -1;
              if( SpecUtils::parse_int( src_attr->value(), src_attr->value_size(), parsed_src_id ) )
                src_id = parsed_src_id;
            }
          }
        }

        // Parse optional line label (new: "l", old: "line")
        string line_label;
        const rapidxml::xml_attribute<char> *line_attr = xray_node->first_attribute( "l", 1 );
        if( !line_attr )
          line_attr = xray_node->first_attribute( "line", 4 );
        if( line_attr )
          line_label = SpecUtils::xml_value_str( line_attr );

        // Create entry and add to database
        XRayWidthEntry entry( z, energy_kev, hwhm_natural_ev, hwhm_doppler_ev, src_id,
                              vacancy_shell, rr, line_label );
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

  // If the XML explicitly contains a line at the requested energy but without
  // width data, we must treat that as "no width" and not substitute a nearby line.
  constexpr double exact_match_eps_kev = 0.001; // 1 eV

  const XRayWidthEntry *best_match = nullptr;
  double best_diff = tolerance_kev;
  const XRayWidthEntry *missing_exact_match = nullptr;

  for( const XRayWidthEntry &entry : entries )
  {
    const double energy_diff = fabs( entry.energy_kev - energy_kev );

    if( energy_diff <= exact_match_eps_kev )
    {
      if( entry.has_width_data() )
        return &entry;

      missing_exact_match = &entry;
      continue;
    }

    if( !entry.has_width_data() )
      continue;

    if( energy_diff < best_diff )
    {
      best_diff = energy_diff;
      best_match = &entry;
    }
  }

  if( missing_exact_match )
    return missing_exact_match;

  return best_match;
}//find_best_match()


double XRayWidthDatabase::get_natural_width_hwhm_kev( const int z,
                                                       const double energy_kev,
                                                       const double tolerance_kev ) const
{
  const XRayWidthEntry *entry = find_best_match( z, energy_kev, tolerance_kev );
  if( !entry || !entry->has_width_data() )
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
  if( !entry || !entry->has_width_data() || entry->hwhm_doppler_ev < 0.0 )
    return -1.0;

  // Scale Doppler width with temperature: HWHM ∝ sqrt(T)
  // Database stores width at reference temperature
  const double temperature_scale = sqrt( temperature_k / sm_doppler_reference_temperature_k );
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
  if( !entry || !entry->has_width_data() || entry->hwhm_doppler_ev < 0.0 )
    return -1.0;

  // Scale Doppler width with temperature
  const double temperature_scale = sqrt( temperature_k / sm_doppler_reference_temperature_k );
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


double get_xray_lorentzian_width( const SandiaDecay::Element *element, const double energy_kev,
                                   const double tolerance_kev )
{
  // Look up natural linewidth (Lorentzian HWHM) from external XML database
  // Data loaded from xray_widths.xml (Campbell & Papp 2001, Krause & Oliver 1979, LBNL)
  // This function is for fluorescent x-rays only (XRF, synchrotron, electron beam)

  if( !element )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_lorentzian_width called with nullptr element" );
    #endif
    return -1.0;
  }

  const std::shared_ptr<const XRayWidthDatabase> db =
    XRayWidthDatabase::instance();

  if( !db )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "XRayWidthDatabase failed to load - xray_widths.xml missing or invalid" );
    #endif
    return -1.0;
  }

  return db->get_natural_width_hwhm_kev( element->atomicNumber, energy_kev, tolerance_kev );
}//double get_xray_lorentzian_width(...)


double get_xray_total_width_for_decay( const SandiaDecay::Transition *transition,
                                       const double xray_energy_kev )
{
  // Returns total x-ray linewidth (natural + alpha recoil Doppler) for decay x-rays
  // Accounts for both natural linewidth and nuclear recoil from alpha decay

  if( !transition )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_total_width_for_decay called with nullptr transition" );
    #endif
    return -1.0;
  }

  // Find the x-ray particle in the transition
  const SandiaDecay::RadParticle *xray_particle = nullptr;
  const double energy_tolerance = 0.5; // keV tolerance for energy matching

  for( const SandiaDecay::RadParticle &particle : transition->products )
  {
    if( particle.type == SandiaDecay::XrayParticle
        && fabs( particle.energy - xray_energy_kev ) < energy_tolerance )
    {
      xray_particle = &particle;
      break;
    }
  }

  if( !xray_particle )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_total_width_for_decay: no x-ray found at specified energy in transition" );
    #endif
    return -1.0;
  }

  // Get the SandiaDecay database to look up the element
  const SandiaDecay::SandiaDecayDataBase * const decay_db = DecayDataBaseServer::database();
  if( !decay_db )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_total_width_for_decay: SandiaDecay database not available" );
    #endif
    return -1.0;
  }

  // Check if this is an alpha decay and compute recoil Doppler broadening
  const bool is_alpha_decay = (transition->mode == SandiaDecay::AlphaDecay
                                || transition->mode == SandiaDecay::BetaAndAlphaDecay
                                || transition->mode == SandiaDecay::ElectronCaptureAndAlphaDecay
                                || transition->mode == SandiaDecay::BetaPlusAndAlphaDecay);

  // Determine which element emits the x-ray
  // For most decays (alpha, beta, electron capture, etc.), the x-ray is emitted by the
  // daughter atom after the decay process creates an atomic vacancy. The only exception
  // is isomeric transitions where the atomic number doesn't change (parent = daughter).
  const SandiaDecay::Nuclide *xray_emitting_nuclide = transition->child;
  
  // If no daughter (e.g., spontaneous fission), fall back to parent
  // For isomeric transitions, parent and child have same atomic number, so either works
  if( !xray_emitting_nuclide )
  {
    xray_emitting_nuclide = transition->parent;
  }

  if( !xray_emitting_nuclide )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_total_width_for_decay: could not determine x-ray emitting nuclide" );
    #endif
    return -1.0;
  }

  const SandiaDecay::Element *element = decay_db->element( xray_emitting_nuclide->atomicNumber );
  if( !element )
  {
    #if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "get_xray_total_width_for_decay: could not find element for atomic number" );
    #endif
    return -1.0;
  }

  // Get natural linewidth (same for parent and daughter if same element)
  const double natural_hwhm_kev = get_xray_lorentzian_width( element, xray_energy_kev, energy_tolerance );
  if( natural_hwhm_kev < 0.0 )
    return -1.0;

  if( !is_alpha_decay )
  {
    // For non-alpha decays, return only natural width
    return natural_hwhm_kev;
  }

  // For alpha decays, compute recoil Doppler broadening
  // The x-ray is emitted by the daughter nucleus after recoil
  const SandiaDecay::Nuclide *parent = transition->parent;
  const SandiaDecay::Nuclide *daughter = transition->child;
  
  if( !parent || !daughter )
  {
    // If no parent or daughter specified, use natural width only
    return natural_hwhm_kev;
  }

  // Compute intensity-weighted recoil Doppler broadening for all alpha energies
  // Different alpha energies give different recoil velocities, so we compute the
  // Doppler width for each alpha and combine them weighted by intensity
  const double alpha_mass_amu = 4.002603;  // Alpha particle mass in amu
  const double amu_to_kev = 931494.0;  // keV per amu (1 amu·c² = 931.494 MeV)
  const double sqrt_2ln2 = 1.177410022515;  // sqrt(2×ln(2))
  const double daughter_mass_amu = daughter->massNumber;

  double weighted_hwhm_sum = 0.0;
  double total_alpha_intensity = 0.0;
  size_t num_alpha_particles = 0;

  // Single loop: find alpha particles and compute their recoil Doppler contributions
  for( const SandiaDecay::RadParticle &particle : transition->products )
  {
    if( particle.type == SandiaDecay::AlphaParticle )
    {
      ++num_alpha_particles;
      const double alpha_intensity = particle.intensity;
      total_alpha_intensity += alpha_intensity;

      // Compute recoil energy from momentum conservation
      // E_recoil = E_alpha × (m_alpha / m_daughter)
      const double recoil_energy_kev = particle.energy * (alpha_mass_amu / daughter_mass_amu);

      // Compute recoil velocity: E = 0.5 × m × v²  →  v = sqrt(2 × E / m)
      const double recoil_velocity_over_c = sqrt( 2.0 * recoil_energy_kev / (daughter_mass_amu * amu_to_kev) );

      // Compute Doppler broadening HWHM for this alpha
      // For random recoil directions (isotropic): HWHM ≈ E × (v/c) / sqrt(2×ln2)
      const double this_hwhm_kev = xray_energy_kev * recoil_velocity_over_c / sqrt_2ln2;

      // Weight by intensity
      weighted_hwhm_sum += this_hwhm_kev * alpha_intensity;
    }
  }

  if( num_alpha_particles == 0 || total_alpha_intensity <= 0.0 )
  {
    // No alpha particles found, return natural width only
    return natural_hwhm_kev;
  }

  const double recoil_hwhm_kev = weighted_hwhm_sum / total_alpha_intensity;

  // Combine natural and recoil widths in quadrature
  const double total_hwhm_kev = sqrt( natural_hwhm_kev * natural_hwhm_kev + recoil_hwhm_kev * recoil_hwhm_kev );

  return total_hwhm_kev;
}//double get_xray_total_width_for_decay(...)


}//namespace XRayWidths
