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

#include <Wt/WString>
#include <Wt/WTemplate>
#include <Wt/WApplication>

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
  std::string more_info_html( const SandiaDecay::Nuclide *const nuc, const bool useBq )
  {
    using namespace Wt;

    if( !nuc || nuc->isStable() )
      return "Invalid Nuclide";

    const SandiaDecay::SandiaDecayDataBase *const db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw runtime_error( "Error getting DecayDataBaseServer" );


    try
    {
      WTemplate tmplt;

      {// begin reading in more_nuc_info_tmplt.xml
        assert( wApp );
        const string filename = "InterSpec_resources/static_text/more_nuc_info_tmplt.xml";
        const string docroot = wApp ? wApp->docRoot() : "";
        const string file_path = SpecUtils::append_path( docroot, filename );

#ifdef _WIN32
        ifstream infile( SpecUtils::convert_from_utf8_to_utf16( file_path ).c_str(), ios::in | ios::binary );
#else
        ifstream infile( file_path.c_str(), ios::in | ios::binary );
#endif
        assert( infile.is_open() );
        if( !infile.is_open() )
          throw runtime_error( "Failed to open more_nuc_info_tmplt.xml" );

        std::stringstream buffer;
        buffer << infile.rdbuf();
        WString tmplt_txt = WString::fromUTF8( buffer.str() );
        tmplt.setTemplateText( tmplt_txt, Wt::TextFormat::XHTMLText );
      }// end reading in more_nuc_info_tmplt.xml

      const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
      assert( el );
      string elname = el ? el->name : string();
      if( !elname.empty() )
        elname[0] = static_cast<char>(std::toupper( static_cast<unsigned char>(elname[0]) ));

      tmplt.bindString( "symbol", nuc->symbol, Wt::TextFormat::XHTMLText );
      tmplt.bindString( "element", elname, Wt::TextFormat::XHTMLText );
      tmplt.bindInt( "mass-number", nuc->massNumber );
      
      const string halflife = PhysicalUnits::printToBestTimeUnits( nuc->halfLife, 4, SandiaDecay::second );
      tmplt.bindString( "half-life", halflife, Wt::TextFormat::XHTMLText );

      const string amu = PhysicalUnits::printCompact( nuc->atomicMass, 6 ) + " amu";
      tmplt.bindString( "atomic-mass", amu, Wt::TextFormat::PlainText );

      const double specificActivity = nuc->activityPerGram() / PhysicalUnits::gram;
      const string sa = PhysicalUnits::printToBestSpecificActivityUnits( specificActivity, 3, !useBq );
      tmplt.bindString( "specific-activity", sa, Wt::TextFormat::XHTMLText );

      double natural_abundance = -1;
      if( el )
      {
        for( const auto &n : el->isotopes )
        {
          if( n.nuclide == nuc )
          {
            natural_abundance = n.abundance;
            break;
          }
        }
      }//if( el )

      tmplt.setCondition( "if-is-natural", (natural_abundance > 0.0) );
      const string natabun = PhysicalUnits::printCompact( natural_abundance, 4 );
      tmplt.bindString( "natural-abundance", natabun, Wt::TextFormat::XHTMLText );

      // I think we need to bind in the JS as XHTMLUnsafeText so it doesnt get stripped out
      tmplt.bindString( "decay-chain-btn",
        "<button type=\"button\" "
          "class=\"Wt-btn with-label\" "
          "onclick=\"Wt.emit($('.specviewer')[0].id,{name:'miscSignal'}, 'decayChain-" + nuc->symbol + "');\""
          ">Show Decay Chain</button>",
        Wt::TextFormat::XHTMLUnsafeText );

      

      tmplt.setCondition( "if-has-parents", !nuc->decaysFromParents.empty() );
      if( !nuc->decaysFromParents.empty() )
      {
        tmplt.bindString( "decay-through-btn",
          "<button type=\"button\" "
          "class=\"Wt-btn with-label\" "
          "onclick=\"Wt.emit($('.specviewer')[0].id,{name:'miscSignal'}, 'decaysThrough-" + nuc->symbol + "');\""
          ">Show Decays Through</button>",
          Wt::TextFormat::XHTMLUnsafeText );
      }

      string decaysToHtml;
      for( const SandiaDecay::Transition *transition : nuc->decaysToChildren )
      {
        if( transition->child )
        {
          decaysToHtml += 
            "<div>Decays to " + transition->child->symbol
            + " by "
            + SandiaDecay::to_str( transition->mode )
            + " decay, BR=" + PhysicalUnits::printCompact( transition->branchRatio, 5 )
            + "</div>";
        }//if( transition->child )
      }//for( const SandiaDecay::Transition * transition : nuc->decaysToChildren)

      tmplt.bindString( "decays-to-html", decaysToHtml, Wt::TextFormat::XHTMLText );
      

      string decaysFromHtml;
      for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )
      {
        const SandiaDecay::Nuclide *parentNuclide = parentTrans->parent;

        if( parentNuclide && (nuc != parentNuclide) )
        {
          const float br_from = parentNuclide->branchRatioFromForebear( nuc );
          const float br_to = parentNuclide->branchRatioToDecendant( nuc );

          if( br_from <= 0.0 || br_to > 0.0 )
            decaysFromHtml += "<div>Branch Ratio from " + parentNuclide->symbol + ": " + PhysicalUnits::printCompact( br_to, 5 ) + "</div>";
          else
            decaysFromHtml += "<div>Branch Ratio through " + parentNuclide->symbol + ": " + PhysicalUnits::printCompact( br_from, 5 ) + "</div>";
        }
      }//for( const SandiaDecay::Transition *parentTrans : nuc->decaysFromParents )

      tmplt.bindString( "decays-from-html", decaysFromHtml, Wt::TextFormat::XHTMLText );

      const auto more_info_db = MoreNucInfoDb::instance();
      if( more_info_db )
      {
        const NucInfo *const more_info = more_info_db->info( nuc );

        if( more_info )
        { 
          string related;
          for( size_t i = 0; i < more_info->m_associated.size(); ++i )
          {
            related += (i ? ", " : "") + more_info->m_associated[i];
            const SandiaDecay::Nuclide *associated = db->nuclide( more_info->m_associated[i] );
            if( associated )
              related += " (hl=" + PhysicalUnits::printToBestTimeUnits( associated->halfLife, 2 ) + ")";
          }
          tmplt.setCondition( "if-has-related", !related.empty() );
          tmplt.bindString( "related-nucs", related, Wt::TextFormat::PlainText );

          string notes = more_info->m_notes;
          SpecUtils::trim( notes ); //just to make sure (although I think the XML parser already did this)

          tmplt.setCondition( "if-has-more-info", !notes.empty() );
          SpecUtils::ireplace_all( notes, "\r", "" );
          vector<string> more_info_lines;
      
          boost::algorithm::split( more_info_lines, notes, 
            boost::algorithm::is_any_of( "\n" ),
            boost::algorithm::token_compress_off );

          notes = "<p>";
          for( size_t i = 0; i < more_info_lines.size(); ++i )
          {
            string line = more_info_lines[i];
            SpecUtils::trim( line );
            
            notes += line;

            if( line.empty() )
            {
              notes += "</p><p>";

              // Skip all the next empty lines
              for( ; (i+1) < more_info_lines.size(); ++i )
              {
                string nextline = more_info_lines[i+1];
                SpecUtils::trim( nextline );
                if( !nextline.empty() )
                  break;
              }//for( skip consequitive empty lines )
            }//if( line.empty() )
          }//for( size_t i = 0; i < more_info_lines.size(); ++i )
          
          notes += "</p>";
          tmplt.bindString( "more-info", notes, Wt::TextFormat::XHTMLText );
        }//if( more_info )
      }//if( more_info_db )

      //if( nuc->canObtainSecularEquilibrium() )
      //  information.push_back( "Can reach secular equilibrium" );
      //else
      //  information.push_back( "Cannot reach secular equilibrium" );
      

      std::stringstream output;
      tmplt.renderTemplate( output );
      
      return output.str();
    }catch( std::exception &e )
    {
      return "Error creating nuclide information: " + string( e.what() );
    }//try / catch

    return "Logic Error";
  }//std::string more_info_html( const SandiaDecay::Nuclide *const nuc );


  InfoStatus more_nuc_info_db_status()
  {
    std::lock_guard<std::mutex> lock(sm_mutex);
    return sm_status;
  }


  RefInfo::RefInfo( string &&key, string &&url, string &&desc )
    : m_key( move( key ) ), m_url( move( url ) ), m_desc( move( desc ) )
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

    const string filename = "OUO_more_nuclide_info.xml";

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

      m_references.emplace( move(key), RefInfo{move(key_cpy), move(url), move(desc)} );
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

    cout << "Took " << total_time.count() << " micro-seconds to initualize MoreNucInfoDb, and it takes up "
      << (memsize()/1024.0) << " kb memory." << endl;
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
