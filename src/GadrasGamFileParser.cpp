/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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

#include <string>
#include <vector>
#include <cstdio>

#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/GadrasGamFileParser.h"

using namespace std;

//std::istream& UtilityFunctions::safe_get_line(std::istream& is, std::string& t)

GadrasGamFile::GadrasGamFile()
{
  //nothing to do here
}//GadrasGamFile constructor



void GadrasGamFile::parse_data( std::istream &strm )
{
  string line;
  if( !UtilityFunctions::safe_get_line( strm, line) )
    throw runtime_error( "GadrasGamFile::parse_data(): Invalid input." );
  
  vector<string> fields;
  UtilityFunctions::split( fields, line, " \t" );
  if( fields.size() < 2 )
    throw runtime_error( "Invalid first line" );
  
  if( (fields[0] == "0") && (fields[0] == "0") )
  {
    //Pre GADRAS version 18.6.0.  Nothing to do here
  }else if( fields[0] == "GamFileVersion" )
  {
    //Post GADRAS version 18.6.0
    /* Header of file will look like:
     GamFileVersion = 2.0
     DataType = SnTransport 18.6.1.0
     Geometry = Spherical
     LargestDimension = 2.200E+00
     Data =
     */
    if( !UtilityFunctions::safe_get_line(strm, line) || !UtilityFunctions::starts_with(line, "DataType ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid second line." );
    if( !UtilityFunctions::safe_get_line( strm, line) || !UtilityFunctions::starts_with(line, "Geometry ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid third line." );
    if( !UtilityFunctions::safe_get_line( strm, line) || !UtilityFunctions::starts_with(line, "LargestDimension ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid fourth line." );
    if( !UtilityFunctions::safe_get_line( strm, line) || !UtilityFunctions::starts_with(line, "Data ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid fifth line." );
  }else
  {
    //A really old gam file?
    throw runtime_error( "Invalid first line for known gam files" );
  }
  
  int photon_lines, photon_groups, neutron_groups;
  if( !UtilityFunctions::safe_get_line( strm, line ) )
    throw runtime_error( "GadrasGamFile::parse_data(): Must have at least two lines." );
  
  if( line.find("! photon lines, photon groups, neutron groups") == string::npos )
    cerr << "Warning, second line of .gam is suspect: '" << line << "'" << endl;
  
  if( sscanf( line.c_str(), "%i %i %i", &photon_lines, &photon_groups, &neutron_groups ) != 3 )
    throw runtime_error( "GadrasGamFile::parse_data(): Second line must have three ints." );

  if( photon_lines < 0 || photon_groups < 0 || neutron_groups < 0
      || photon_lines > 32000 || photon_groups > 32000 || neutron_groups > 32000 )
  {
    throw runtime_error( "GadrasGamFile::parse_data(): Invalid number of lines." );
  }
  
  m_photon_lines_energy.resize( photon_lines );
  m_photon_lines_flux.resize( photon_lines );
  m_photon_lines_an.resize( photon_lines );
  m_photon_lines_ad.resize( photon_lines );
  
  for( int i = 0; i < photon_lines; ++i )
  {
    if( !UtilityFunctions::safe_get_line( strm, line) )
      throw runtime_error( "GadrasGamFile::parse_data(): Claimed wrong number of photon lines." );
    
    if( !i )
    {
      const size_t pos = line.find( "! photon lines" );
      if( pos == string::npos )
        cerr << "Warning, didnt find \"! photon lines\" like expected\n";
      else
        line = line.substr( 0, pos );
    }//if( !i )
    
    vector<float> linevaleus;
    const bool splitres = UtilityFunctions::split_to_floats( line.c_str(), line.size(), linevaleus );
    if( !splitres || (linevaleus.size() != 4) )
      throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 4 gamma line values." );
    
    m_photon_lines_energy[i] = linevaleus[0];
    m_photon_lines_flux[i] = linevaleus[1];
    m_photon_lines_an[i] = linevaleus[2];
    m_photon_lines_ad[i] = linevaleus[3];
    
/*
    const int nitem = sscanf( line.c_str(), "%f %f %f %f",
                              &(m_photon_lines_energy[i]),
                              &(m_photon_lines_flux[i]),
                              &(m_photon_lines_an[i]),
                              &(m_photon_lines_ad[i]) );
    if( nitem != 4 )
      throw runtime_error( "GadrasGamFile::parse_data(): Didnt get 4 gamma lines." );
*/
  }//for( int i = 0; i < photon_lines; ++i )
  
  //The file lists one more line than claimed number of photon or neutron groups
  //  in order to give the upper energy of the last group.  So we will incrememt
  //  the number of groups here if there are any groups
  photon_groups += (photon_groups ? 1 : 0);
  neutron_groups += (neutron_groups ? 1 : 0);
  
  m_photon_group_boundries.resize( photon_groups );
  m_photon_group_flux.resize( photon_groups );
  for( int i = 0; i < photon_groups; ++i )
  {
    if( !UtilityFunctions::safe_get_line( strm, line) )
      throw runtime_error( "GadrasGamFile::parse_data(): Claimed wrong number of photon groups." );
    
    if( !i )
    {
      const size_t pos = line.find( "! photon groups" );
      if( pos == string::npos )
        cerr << "Warning, didnt find \"! photon groups\" like expected\n";
      else
        line = line.substr( 0, pos );
    }//if( !i )
    
    vector<float> linevaleus;
    const bool splitres = UtilityFunctions::split_to_floats( line.c_str(), line.size(), linevaleus );
    if( !splitres || (linevaleus.size() != 2) )
      throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 2 gamma group values." );
    
    m_photon_group_boundries[i] = linevaleus[0];
    m_photon_group_flux[i] = linevaleus[1];
  }//for( int i = 0; i < photon_groups; ++i )
  
  
  
  m_neutron_group_boundries.resize( neutron_groups );
  m_neutron_group_flux.resize( neutron_groups );
  for( int i = 0; i < neutron_groups; ++i )
  {
    if( !UtilityFunctions::safe_get_line( strm, line) )
      throw runtime_error( "GadrasGamFile::parse_data(): Claimed wrong number of neutron groups." );
    
    if( !i )
    {
      const size_t pos = line.find( "! neutron groups" );
      if( pos == string::npos )
        cerr << "Warning, didnt find \"! neutron groups\" like expected\n";
      else
        line = line.substr( 0, pos );
    }//if( !i )
    
    vector<float> linevaleus;
    const bool splitres = UtilityFunctions::split_to_floats( line.c_str(), line.size(), linevaleus );
    if( !splitres || (linevaleus.size() != 2) )
      throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 2 neutron group values." );
    
    m_neutron_group_boundries[i] = 0.001*linevaleus[0];
    m_neutron_group_flux[i] = linevaleus[1];
  }//for( int i = 0; i < neutron_groups; ++i )
  
  /*  File will then possible have:
   !
   0.767265975475311 ! k-eff
   4.29675034427086 ! multiplication
   */
  
  cerr << "photon_lines=" << photon_lines << ", photon_groups=" << photon_groups << ", neutron_groups=" << neutron_groups << endl;
}//void parse_data( std::istream &strm )
