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

#include "InterSpec/ReferenceLinePredef.h"

#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace std;

namespace ReferenceLinePredef
{

namespace
{
  void sanitize_label_str( string &label )
  {
    SpecUtils::trim( label );
    SpecUtils::ireplace_all( label, " ", "" );
    SpecUtils::to_lower_ascii( label );
  }//void sanitize_label_str( string &label )
}//namespace


void load_ref_line_file( const string& filepath, 
                        map<string,NucMix>& nuc_mixes, 
                        map<string,CustomSrcLines>& custom_lines,
                        std::vector<IndividualSource> *individual_sources )
{
  try
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      throw std::logic_error( "invalid SandiaDecayDataBase" );
    
    std::vector<char> data;
    SpecUtils::load_file_data( filepath.c_str(), data );
    
    rapidxml::xml_document<char> doc;
    const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
    doc.parse<flags>( &data.front() );
    
    typedef rapidxml::xml_node<char>      XmlNode;
    typedef rapidxml::xml_attribute<char> XmlAttribute;
    
    const XmlNode *ref_lines = doc.first_node( "RefLineDefinitions" );
    if( !ref_lines )
      throw runtime_error( "No RefLineDefinitions node." );
    
    XML_FOREACH_CHILD( nuc_mix, ref_lines, "NucMixture" )
    {
      const XmlAttribute *mix_name = XML_FIRST_ATTRIB( nuc_mix, "name" );
      const XmlAttribute *def_age = XML_FIRST_ATTRIB( nuc_mix, "default-age" );
      
      string mix_name_str = SpecUtils::xml_value_str( mix_name );
      if( mix_name_str.empty() )
        throw runtime_error( "No mixture name" );
      
      NucMix mix;
      mix.m_name = mix_name_str;
      mix.m_default_age = 0.0;
      mix.m_default_age_str = "";
      mix.m_fixed_act_fractions = true;
      if( def_age )
      {
        mix.m_default_age_str = SpecUtils::xml_value_str( def_age );
        mix.m_default_age = PhysicalUnits::stringToTimeDuration( mix.m_default_age_str );
      }
      
      const XmlAttribute *ref_age_attrib = XML_FIRST_ATTRIB( nuc_mix, "reference-age" );
      const string ref_age_str = SpecUtils::xml_value_str( ref_age_attrib );
      double ref_age = -1.0;
      if( !ref_age_str.empty() )
      {
        mix.m_fixed_act_fractions = false;
        ref_age = PhysicalUnits::stringToTimeDuration( ref_age_str );
      }//if( !ref_age_str.empty() )
      
      
      const XmlAttribute *weight_attrib = XML_FIRST_ATTRIB( nuc_mix, "weight" );
      if( weight_attrib && weight_attrib->value_size() )
      {
        if( !SpecUtils::parse_float(weight_attrib->value(), weight_attrib->value_size(), mix.m_weight ) )
          cerr << "Failed to parse mixture '" << mix_name_str << "' weight attribute." << endl;
      }
      
      XML_FOREACH_CHILD( nuc, nuc_mix, "Nuc" )
      {
        const XmlAttribute *nuc_name = XML_FIRST_ATTRIB( nuc, "name" );
        const XmlAttribute *nuc_act_frac = XML_FIRST_ATTRIB( nuc, "act-frac" );
        const XmlAttribute *age_offset = XML_FIRST_ATTRIB( nuc, "age-offset" );
        const XmlAttribute *color = XML_FIRST_ATTRIB( nuc, "color" );
        
        if( !nuc_name || !nuc_name->value_size() )
          throw runtime_error( "No nuclide name" );
        
        const string nuc_name_str = SpecUtils::xml_value_str( nuc_name );
        
        if( !nuc_act_frac || !nuc_act_frac->value_size() )
          throw runtime_error( "No activity fraction for " + nuc_name_str
                              + " in " + mix_name_str );
        
        NucMixComp comp;
        comp.m_age_offset = 0.0;
        comp.m_nuclide = db->nuclide( nuc_name_str );
        if( !comp.m_nuclide )
          throw runtime_error( "Invalid nuclide: " + nuc_name_str );
        
        if( !SpecUtils::parse_double( nuc_act_frac->value(), nuc_act_frac->value_size(),
                                     comp.m_rel_act )
           || (comp.m_rel_act < 0.0) )
          throw runtime_error( "Invalid activity fraction: " + nuc_name_str );
        
        if( age_offset && age_offset->value_size() )
        {
          const string age_offset_str = SpecUtils::xml_value_str( age_offset );
          comp.m_age_offset = PhysicalUnits::stringToTimeDuration( age_offset_str );
        }//if( age_offset && age_offset->value_size() )
        
        if( color && color->value_size() )
        {
          const string color_str = SpecUtils::xml_value_str(color);
          comp.m_color = Wt::WColor( color_str );
          if( comp.m_color.isDefault() )
          {
            const string msg = "NucMixture named '" + mix_name_str + "'"
            " has invalid color value for '" + nuc_name_str + "': "
            "'" + color_str + "' - not CSS color string.";
            cerr << msg << endl;
            //throw runtime_error( msg ); // we wont disregard the whole file over a single color
          }//if( parsing of color failed )
        }//if( user specified color )
        
        mix.m_components.push_back( std::move(comp) );
      }//XML_FOREACH_CHILD( nuc, nuc_mix, "Nuc" )
      
      assert( mix.m_fixed_act_fractions == (ref_age < 0.0) );
      if( !mix.m_fixed_act_fractions && (ref_age > 0.0) )
      {
        // Figure out activities at t=0:
        //  specified_act = A_0 * exp( -ref_age * nuc->decayConstant() );
        //  A_0 = specified_act / exp( -ref_age * nuc->decayConstant() );
        
        for( NucMixComp &m : mix.m_components )
        {
          const double given_act = m.m_rel_act;
          const SandiaDecay::Nuclide *nuc = m.m_nuclide;
          m.m_rel_act = given_act / exp( -ref_age * nuc->decayConstant() );
          
#if( PERFORM_DEVELOPER_CHECKS )
          SandiaDecay::NuclideMixture decay_mix;
          decay_mix.addNuclideByActivity( nuc, m.m_rel_act );
          const double ref_act = decay_mix.activity(ref_age, nuc);
          if( (given_act > 0.0)
             && (fabs(ref_act - given_act) > 1.0E-5*std::max(ref_act,given_act)) ) //1.0E-5 arbitrary
          {
            log_developer_error( __func__, "Decay correction calculation has failed." );
            assert( fabs(ref_act - given_act) < 1.0E-5*std::max(ref_act,given_act) );
          }
#endif
        }//for( NucMixComp &m : mix.m_components )
      }//if( mix.m_fixed_act_fractions )
      
      double act_fraction_sum = 0.0;
      for( const NucMixComp &m : mix.m_components )
        act_fraction_sum += m.m_rel_act;
      
      if( (act_fraction_sum <= 0.0) || IsNan(act_fraction_sum) || IsInf(act_fraction_sum) )
        throw runtime_error( "Invalid activity fraction sum" );
      
      for( NucMixComp &m : mix.m_components )
        m.m_rel_act /= act_fraction_sum;
      
      sanitize_label_str( mix_name_str );
      nuc_mixes[mix_name_str] = std::move(mix);
    }//XML_FOREACH_CHILD( nuc_mix, ref_lines, "NucMixture" )
    
    
    XML_FOREACH_CHILD( source, ref_lines, "SourceLines" )
    {
      const XmlAttribute *src_name = XML_FIRST_ATTRIB( source, "name" );
      string src_name_str = SpecUtils::xml_value_str( src_name );
      SpecUtils::trim( src_name_str );
      
      if( src_name_str.empty() )
        throw runtime_error( "No name specified for a SourceLines element" );
      
      CustomSrcLines src_lines;
      src_lines.m_name = src_name_str;
      src_lines.m_max_branch_ratio = 0.0;
      
      const XmlAttribute *weight_attrib = XML_FIRST_ATTRIB( source, "weight" );
      if( weight_attrib && weight_attrib->value_size() )
      {
        if( !SpecUtils::parse_float(weight_attrib->value(), weight_attrib->value_size(), src_lines.m_weight ) )
          cerr << "Failed to parse source '" << src_name_str << "' weight attribute." << endl;
      }
      
      XML_FOREACH_CHILD( line, source, "Line" )
      {
        const XmlAttribute *info = XML_FIRST_ATTRIB( line, "info" );
        string info_str = SpecUtils::xml_value_str( info );
        SpecUtils::trim( info_str );
        
        const string values_str = SpecUtils::xml_value_str( line );
        
        vector<float> values;
        SpecUtils::split_to_floats( values_str, values );
        if( values.size() != 2 )
          throw runtime_error( "SourceLines named '" + src_name_str + "' provided "
                              + std::to_string(values.size()) + " values (expected two numbers)" );
        
        const float energy = values[0];
        const float br = values[1];
        
        if( (energy <= 0.f) || (br < 0.0f) )
          throw runtime_error( "SourceLines named '" + src_name_str + "' has a negative value." );
        
        src_lines.m_max_branch_ratio = std::max( src_lines.m_max_branch_ratio, br );
        
        const SandiaDecay::Nuclide *nuc = nullptr;
        const SandiaDecay::Transition *trans = nullptr;
        const XmlAttribute *nuclide = XML_FIRST_ATTRIB( line, "nuc" );
        if( nuclide && nuclide->value_size() )
        {
          const string nuc_name = SpecUtils::xml_value_str(nuclide);
          nuc = db->nuclide( nuc_name );
          if( !nuc )
            throw runtime_error( "SourceLines named '" + src_name_str + "' has an invalid nuclide ('" + nuc_name + "')." );
          
          size_t transition_index = 0;
          PeakDef::SourceGammaType nearestGammaType;
          PeakDef::findNearestPhotopeak( nuc, energy, 0.25, false, -1.0,
                                        trans, transition_index, nearestGammaType );
          if( !trans )
            throw runtime_error( "SourceLines named '" + src_name_str + "' with nuclide '"
                                + nuc_name + "', couldnt be matched to source data for energy "
                                + std::to_string(energy) + " keV." );
        }//if( nuclide && nuclide->value_size() )
        
        bool atten_applies = true;
        const XmlAttribute *atten = XML_FIRST_ATTRIB( line, "atten" );
        if( atten && atten->value_size() )
        {
          const string atten_str = SpecUtils::xml_value_str(atten);
          if( SpecUtils::iequals_ascii(atten_str, "0")
             || SpecUtils::iequals_ascii(atten_str, "false")
             || SpecUtils::iequals_ascii(atten_str, "no") )
          {
            atten_applies = false;
          }else if( !SpecUtils::iequals_ascii(atten_str, "1")
                   && !SpecUtils::iequals_ascii(atten_str, "true")
                   && !SpecUtils::iequals_ascii(atten_str, "yes") )
          {
            throw runtime_error( "SourceLines named '" + src_name_str + "' has invalid value "
                                " for the 'atten' attribute ('" + atten_str + "')" );
          }
        }//if( atten && atten->value_size() )
        
        Wt::WColor color;
        const XmlAttribute *color_attrib = XML_FIRST_ATTRIB( line, "color" );
        if( color_attrib && color_attrib->value_size() )
        {
          const string color_str = SpecUtils::xml_value_str(color_attrib);
          color = Wt::WColor( color_str );
          if( color.isDefault() )
          {
            // Parsing of color failed
            const string msg =  "SourceLines named '" + src_name_str + "' has invalid value "
            " for color attribute ('" + color_str
            + "') - could not parse as CSS color.";
            cerr << msg << endl;
            //throw runtime_error( msg ); //We wont disregard the whole file over an invalid file
          }
        }//if( color_attrib && color_attrib->value_size() )
        
        src_lines.m_lines.emplace_back( energy, br, std::move(info_str), nuc,
                                       trans, atten_applies, std::move(color) );
      }//XML_FOREACH_CHILD( line, source, "Line" )
      
      if( src_lines.m_lines.empty() )
        throw runtime_error( "No lines specified for SourceLines named '" + src_name_str + "'" );
      
      if( src_lines.m_max_branch_ratio <= 0.0f )
        throw runtime_error( "Lines specified for SourceLines named '" + src_name_str + "' were all zero amplitude." );
      
      sanitize_label_str( src_name_str );
      custom_lines[std::move(src_name_str)] = std::move(src_lines);
    }//XML_FOREACH_CHILD( source, ref_lines, "SourceLines" )
    
    if( individual_sources )
    {
      XML_FOREACH_CHILD( ind_src, ref_lines, "IndividualSource" )
      {
        const XmlAttribute *src_name = XML_FIRST_ATTRIB( ind_src, "name" );
        const XmlAttribute *weight_attrib = XML_FIRST_ATTRIB( ind_src, "weight" );
        const XmlAttribute *def_age = XML_FIRST_ATTRIB( ind_src, "default-age" );
        const XmlAttribute *color_attrib = XML_FIRST_ATTRIB( ind_src, "color" );
        const XmlAttribute *is_background_attrib = XML_FIRST_ATTRIB( ind_src, "is-background" );
        const XmlAttribute *shield_material_attrib = XML_FIRST_ATTRIB( ind_src, "shield-material" );
        const XmlAttribute *shield_thickness_attrib = XML_FIRST_ATTRIB( ind_src, "shield-thickness" );
        
        string src_name_str = SpecUtils::xml_value_str( src_name );
        SpecUtils::trim( src_name_str );
        if( src_name_str.empty() )
          continue;
        
        IndividualSource src;
        src.m_name = src_name_str;
        
        if( weight_attrib && weight_attrib->value_size() )
        {
          if( !SpecUtils::parse_float(weight_attrib->value(), weight_attrib->value_size(), src.m_weight ) )
            cerr << "Failed to parse source '" << src_name_str << "' weight attribute." << endl;
        }
        
        src.m_nuclide = db->nuclide( src_name_str );
        if( !src.m_nuclide )
          src.m_element = db->element( src_name_str );
        
        if( !src.m_nuclide && !src.m_element )
        {
          const ReactionGamma * const reaction_db = ReactionGammaServer::database();
          assert( reaction_db );
          if( !reaction_db )
            throw runtime_error( "Couldnt load reaction DB" );
          
          vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
          reaction_db->gammas( src_name_str, possible_rctns );
          
          // TODO: we are currently taking the first reaction; however, in principle there could be multiple - however, `ReactionGamma` doesnt have an interface to just return a reaction by name, I guess because
          for( size_t i = 0; !src.m_reaction && (i < possible_rctns.size()); ++i )
            src.m_reaction = possible_rctns[i].reaction;
          
          if( !src.m_reaction )
            throw runtime_error( "Invalid source not cooresponding to a nuclide, element, or reaction: '" + src_name_str + "'" );
        }//if( !src.m_nuclide && !src.m_element )
        
        if( src.m_nuclide && def_age && def_age->value_size() )
        {
          const string default_age_str = SpecUtils::xml_value_str( def_age );
          src.m_age = PhysicalUnits::stringToTimeDuration( default_age_str );
        }//if( def_age && def_age->value_size() )
        
        if( color_attrib && color_attrib->value_size() )
        {
          const string color_str = SpecUtils::xml_value_str(color_attrib);
          Wt::WColor color = Wt::WColor( color_str );
          if( !color.isDefault() )
          {
            src.m_color = color;
          }else
          {
            // Parsing of color failed
            const string msg =  "SourceLines named '" + src_name_str + "' has invalid value "
            " for color attribute ('" + color_str
            + "') - could not parse as CSS color.";
            cerr << msg << endl;
            //throw runtime_error( msg ); //We wont disregard the whole file over an invalid color
          }
        }//if( color_attrib && color_attrib->value_size() )
        
        
        if( is_background_attrib )
        {
          const string atten_str = SpecUtils::xml_value_str(is_background_attrib);
          if( !SpecUtils::iequals_ascii(atten_str, "1")
             && !SpecUtils::iequals_ascii(atten_str, "true")
             && !SpecUtils::iequals_ascii(atten_str, "yes") )
          {
            src.is_background = true;
          }
        }//if( is_background_attrib )
        
        
        const string material_str = SpecUtils::xml_value_str(shield_material_attrib);
        if( !material_str.empty() )
        {
          src.shielding_material = material_str;
          const string thickness_str = SpecUtils::xml_value_str(shield_thickness_attrib);
          if( thickness_str.empty() )
            throw runtime_error( "Shielding thickness must be specified if shielding material is" );
          src.shielding_thickness = PhysicalUnits::stringToDistance( thickness_str );
        }//if( shield_material_attrib )
        
        individual_sources->push_back( std::move(src) );
      }//XML_FOREACH_CHILD( ind_src, ref_lines, "IndividualSource" )
    }//if( individual_sources )
  }catch( std::exception &e )
  {
    cerr << "Failed to load '" << filepath << "' as custom ref lines: "
    << e.what() << endl;
  }//try / catch to load XML data
}//load_ref_line_file(...)

} //namespace ReferenceLinePredef
