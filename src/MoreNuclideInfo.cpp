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

#include <chrono> // just for timing things
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

namespace
{
  std::mutex sm_mutex;
  MoreNuclideInfo::InfoStatus sm_status = MoreNuclideInfo::InfoStatus::NotInited;
  std::shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> sm_db;
}


namespace MoreNuclideInfo
{
  InfoStatus more_nuc_info_db_status()
  {
    std::lock_guard<std::mutex> lock(sm_mutex);
    return sm_status;
  }


  RefInfo::RefInfo( string &&key, string &&url, string &&desc )
    : m_key( std::move( key ) ), m_url( std::move( url ) ), m_desc( std::move( desc ) )
  {
  }

  
  size_t NucInfo::memsize() const
  {
    size_t len = sizeof( this )
      + m_nuclide.capacity()
      + m_notes.capacity()
      + (m_associated.capacity() * sizeof( string ));

    for( const string &val : m_associated )
      len += val.capacity();

    return len;
  }//memsize()


  size_t RefInfo::memsize() const
  {
    return sizeof( this )
      + m_key.capacity()
      + m_url.capacity()
      + m_desc.capacity();
  }//memsize()
  

  size_t MoreNucInfoDb::memsize() const
  {
    size_t len = sizeof( this )
      + (m_nuc_infos.size() * sizeof( const SandiaDecay::Nuclide * ))
      + (m_other_infos.size() * sizeof( string ))
      + (m_references.size() * sizeof( string ));
    
    for( const auto &v : m_nuc_infos )
      len += v.second.memsize();

    for( const auto &v : m_other_infos )
      len += v.first.capacity() + v.second.memsize();

    for( const auto &v : m_references )
      len += v.first.capacity() + v.second.memsize();

    return len;
  }//size_t MoreNucInfoDb::memsize()


  std::shared_ptr<const MoreNucInfoDb> MoreNucInfoDb::instance()
  {
    std::lock_guard<std::mutex> lock( sm_mutex );
    if( sm_db )
    {
      assert( sm_status == MoreNuclideInfo::InfoStatus::Inited );
      return sm_db;
    }

    assert( sm_status != MoreNuclideInfo::InfoStatus::Inited );

    if( sm_status == MoreNuclideInfo::InfoStatus::FailedToInit )
      return nullptr;

    assert( sm_status == MoreNuclideInfo::InfoStatus::NotInited );

    try
    {
      sm_status = MoreNuclideInfo::InfoStatus::FailedToInit;

      auto db = std::make_shared<MoreNucInfoDb>();
      db->init();
      
      sm_db = db;
      sm_status = MoreNuclideInfo::InfoStatus::Inited;
    }catch( std::exception &e )
    {
      cerr << "MoreNucInfoDb - failed to init: " << e.what() << endl;
    }// try / catch

    return sm_db;
  }//std::shared_ptr<const MoreNucInfoDb> MoreNucInfoDb::instance()

  
  void MoreNucInfoDb::remove_global_instance()
  {
    std::lock_guard<std::mutex> lock( sm_mutex );
    if( sm_db )
    {
      assert( sm_status == MoreNuclideInfo::InfoStatus::Inited );
      sm_db.reset();
      sm_status = MoreNuclideInfo::InfoStatus::NotInited;
    }else
    {
      assert( sm_status != MoreNuclideInfo::InfoStatus::Inited );
    }
  }//void MoreNucInfoDb::remove_global_instance()
  
  
  MoreNucInfoDb::MoreNucInfoDb()
  {
  }


  void MoreNucInfoDb::init()
  {
    // This function seemingly only takes:
    //  "Took 2112 micro-seconds to initualize MoreNucInfoDb, and it takes up 43.6357 kb memory."
    // To go through, but will leave debug metric printing out till I try it on a few more
    //  computers.  Debug mode takes 24 ms to initualize.
    const auto start_time = std::chrono::steady_clock::now();
    
    m_nuc_infos.clear();
    m_other_infos.clear();
    m_references.clear();

    const string filename = "more_nuclide_info.xml";

    const SandiaDecay::SandiaDecayDataBase *const db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw runtime_error( "Error getting DecayDataBaseServer" );


    ifstream infile;

    // First we'll look in the users writable data directory
    //  e.g., 'C:\Users\<username>\AppData\Roaming\InterSpec' for the more info
    //  file, and if not there, we'll look for the file distributed with 
    //  InterSpec
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
      }catch(std::exception &)
      {
        //InterSpec::writableDataDirectory throws exception if it hasnt been set
      }
    }//if( try to open users more_nuclide_info.xml )
#endif

    // If we didnt open the file, try the one that comes with InterSpec
    if( !infile.is_open() )
    {
      const string data_dir = InterSpec::staticDataDirectory();
      const string file_path = SpecUtils::append_path( data_dir, filename );
#ifdef _WIN32
      infile.open( SpecUtils::convert_from_utf8_to_utf16( file_path ).c_str(), ios::in | ios::binary );
#else
      infile.open( file_path.c_str(), ios::in | ios::binary );
#endif
    }//if( try to open default more_nuclide_info.xml )
    
    // Finally, try the CWD
    if( !infile.is_open() )
    {
#ifdef _WIN32
      infile.open( SpecUtils::convert_from_utf8_to_utf16( filename ).c_str(), ios::in | ios::binary );
#else
      infile.open( filename.c_str(), ios::in | ios::binary );
#endif
    }//if( try to open default more_nuclide_info.xml )

    
    if( !infile.is_open() )
      throw runtime_error( "MoreNucInfoDb::init: couldnt open '" + filename + "'" );

    rapidxml::file<char> input_data( infile );
    rapidxml::xml_document<char> doc;

    try
    {
      // We'll parse non-destructively, meaning we will have to use
      // xml_base::name_size() and xml_base::value_size() and also entities
      // will not be translated
      doc.parse<rapidxml::parse_trim_whitespace | rapidxml::parse_non_destructive>( input_data.data() );
    }catch( rapidxml::parse_error &e )
    {
      string msg = "Error parsing MoreNucInfoDb XML: " + string( e.what() );
      const char *const position = e.where<char>();
      if( position && *position )
      {
        const char *end_pos = (const char *)input_data.data() + input_data.size();
        end_pos = std::min( position + 80, end_pos );

        msg += "\n\tAt: " + std::string( position, end_pos );
      }//if( position )

      throw runtime_error( msg );
    }//try / catch
 
    const rapidxml::xml_node<char> * const base_node = XML_FIRST_NODE( &doc, "AddNucInfo" );
    if( !base_node )
      throw runtime_error( filename + " doesnt have base-node 'AddNucInfo'" );

    const rapidxml::xml_node<char> * const references_node = XML_FIRST_NODE( base_node, "References" );
    if( !references_node )
      throw runtime_error( filename + " doesnt have node 'References'" );

    const rapidxml::xml_node<char> * const nuclides_node = XML_FIRST_NODE( base_node, "Nuclides" );
    if( !nuclides_node )
      throw runtime_error( filename + " doesnt have node 'Nuclides'" );


    XML_FOREACH_CHILD(ref_node, references_node, "Ref")
    {
      const rapidxml::xml_attribute<char> *const key_att = XML_FIRST_ATTRIB( ref_node, "key");
      // The following assert is for developing the XML, not coding logic
      assert( key_att && key_att->value_size() );

      if( !key_att || !key_att->value_size() )
        continue;
      
      const rapidxml::xml_node<char> *const desc_node = XML_FIRST_NODE( ref_node, "desc" );
      const rapidxml::xml_node<char> *const url_node = XML_FIRST_NODE( ref_node, "url" );

      std::string key = SpecUtils::xml_value_str( key_att );
      std::string key_cpy = key;
      std::string url = SpecUtils::xml_value_str( url_node );
      std::string desc = SpecUtils::xml_value_str( desc_node );

      // The following assert is for developing the XML, not coding logic
      assert( !(url.empty() && desc.empty()) );

      m_references.emplace( std::move(key), RefInfo{std::move(key_cpy), std::move(url), std::move(desc)} );
    }//for( loop over <Ref> nodes )


    XML_FOREACH_CHILD( nuc_node, nuclides_node, "Nuc" )
    {
      const rapidxml::xml_attribute<char> *const name_att = XML_FIRST_ATTRIB( nuc_node, "name" );
      
      // The following assert is for developing the XML, not coding logic
      assert( name_att && name_att->value_size() );

      if( !name_att || !name_att->value_size() )
        continue;

      const rapidxml::xml_node<char> *const assoc_node = XML_FIRST_NODE( nuc_node, "Associated" );
      const rapidxml::xml_node<char> *const notes_node = XML_FIRST_NODE( nuc_node, "Notes" );

      NucInfo info;
      info.m_nuclide = SpecUtils::xml_value_str( name_att );
      info.m_notes = SpecUtils::xml_value_str( notes_node );

      const string associated = SpecUtils::xml_value_str( assoc_node );
      SpecUtils::split( info.m_associated, associated, ";" );
      for( std::string &val : info.m_associated )
        SpecUtils::trim( val );

      // The following assert is for developing the XML, not coding logic
      assert( !info.m_notes.empty() || !info.m_associated.empty() );

      const SandiaDecay::Nuclide * nuc = db->nuclide( info.m_nuclide );
      if( nuc )
      {
        m_nuc_infos.emplace( nuc, std::move( info ) );
      } else
      {
        // I think the next line would be fine, but havent checked it
        //m_other_infos[info.m_nuclide] = std::move( info );
        const string nucstr = info.m_nuclide;
        m_other_infos.emplace( std::move(nucstr), std::move( info ) );
      }
    }//XML_FOREACH_CHILD( nuc_node, nuclides_node, "Nuc" )

    if( m_nuc_infos.empty() )
      throw runtime_error( "MoreNucInfoDb::init: Failed to read in any nuclides" );

    const auto finish_time = std::chrono::steady_clock::now();
    const auto total_time = chrono::duration_cast<chrono::microseconds>(finish_time - start_time);

    //cout << "Took " << total_time.count() << " micro-seconds to initualize MoreNucInfoDb, and it takes up "
    //  << (memsize()/1024.0) << " kb memory." << endl;
  }//void init()

  
  const NucInfo *MoreNucInfoDb::info( const std::string &nucstr ) const
  {
    const SandiaDecay::SandiaDecayDataBase *const db = DecayDataBaseServer::database();
    assert( db );
    if( db )
    {
      const auto nuc = db->nuclide( nucstr );
      if( nuc )
        return info( nuc );
    }

    const auto pos = m_other_infos.find( nucstr );
    if( pos == end( m_other_infos ) )
      return nullptr;

    return &(pos->second);
  }//info(...)


  const NucInfo *MoreNucInfoDb::info( const SandiaDecay::Nuclide *nuc ) const
  {
    if( !nuc )
      return nullptr;

    const auto pos = m_nuc_infos.find( nuc );
    if( pos == end( m_nuc_infos ) )
      return nullptr;

    return &(pos->second);
  }//info(...)

}//namespace MoreNuclideInfo
