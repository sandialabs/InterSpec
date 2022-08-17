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

#include <string>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <stdexcept>

#include "SpecUtils/ParseUtils.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/GadrasGamFileParser.h"

using namespace std;

//std::istream& SpecUtils::safe_get_line(std::istream& is, std::string& t)

GadrasGamFile::GadrasGamFile()
{
  m_version = GadrasGamFile::GamFileVersion::NumGamFileVersions;
}//GadrasGamFile constructor



void GadrasGamFile::parse_data( std::istream &strm )
{
  GamFileVersion version = GamFileVersion::NumGamFileVersions;
  
  string line;
  if( !SpecUtils::safe_get_line( strm, line) )
    throw runtime_error( "GadrasGamFile::parse_data(): Invalid input." );
  
  vector<string> fields;
  SpecUtils::split( fields, line, " \t" );
  if( fields.size() < 2 )
    throw runtime_error( "Invalid first line" );
  
  if( (fields[0] == "0") && (fields[1] == "0") )
  {
    //Pre GADRAS version 18.6.0.  Nothing to do here
    version = GadrasGamFile::GamFileVersion::Version_1;
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
    
    if( fields.size() < 3 || ((fields[2] != "2.0") && (fields[2] != "2.1")) )
      throw runtime_error( "GadrasGamFile::parse_data(): Unknown gam file version: " + line );
    
    if( fields[2] == "2.0" )
      version = GadrasGamFile::GamFileVersion::Version_2_0;
    else
      version = GadrasGamFile::GamFileVersion::Version_2_1;
    
    if( !SpecUtils::safe_get_line(strm, line) || !SpecUtils::starts_with(line, "DataType ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid \"DataType\" line." );
    
    if( !SpecUtils::safe_get_line( strm, line) || !SpecUtils::starts_with(line, "Geometry ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid \"Geometry\" line." );
    
    if( version == GadrasGamFile::GamFileVersion::Version_2_1 )
    {
      if( !SpecUtils::safe_get_line( strm, line) || !SpecUtils::starts_with(line, "Description ") )
        throw runtime_error( "GadrasGamFile::parse_data(): Invalid \"Description\" line." );
    }
    
    if( !SpecUtils::safe_get_line( strm, line) || !SpecUtils::starts_with(line, "LargestDimension ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid \"LargestDimension\" line." );
    
    if( !SpecUtils::safe_get_line( strm, line) || !SpecUtils::starts_with(line, "Data ") )
      throw runtime_error( "GadrasGamFile::parse_data(): Invalid \"Data\" line." );
    
  }else
  {
    //A really old gam file?
    throw runtime_error( "Invalid first line for known gam files" );
  }
  
  int photon_lines, photon_groups, neutron_groups;
  if( !SpecUtils::safe_get_line( strm, line ) )
    throw runtime_error( "GadrasGamFile::parse_data(): Must have at least two lines." );
  
  if( line.find("! photon lines, photon groups, neutron groups") == string::npos )
    cerr << "Warning, second line of .gam is suspect: '" << line << "'" << endl;
  
  if( sscanf( line.c_str(), "%i %i %i", &photon_lines, &photon_groups, &neutron_groups ) != 3 )
    throw runtime_error( "GadrasGamFile::parse_data(): Second line must have three ints." );
  
  if( photon_lines < 0 || photon_groups < 0 || neutron_groups < 0
     || photon_lines > 32000 /* || photon_groups > 32000 */ || neutron_groups > 32000 )
  {
    throw runtime_error( "GadrasGamFile::parse_data(): Invalid number of lines." );
  }
  
  m_photon_lines_energy.resize( photon_lines, 0.0f );
  m_photon_lines_flux.resize( photon_lines, 0.0f );
  m_photon_lines_an.resize( photon_lines, 0.0f );
  m_photon_lines_ad.resize( photon_lines, 0.0f );
  
  for( int i = 0; i < photon_lines; ++i )
  {
    if( !SpecUtils::safe_get_line( strm, line) )
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
    const bool splitres = SpecUtils::split_to_floats( line.c_str(), line.size(), linevaleus );
    
    switch( version )
    {
      case GamFileVersion::Version_1:
      case GamFileVersion::Version_2_0:
      case GamFileVersion::NumGamFileVersions:
        if( !splitres || (linevaleus.size() != 4) )
          throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 4 gamma line values." );
        break;
        
      case GamFileVersion::Version_2_1:
        if( !splitres || (linevaleus.size() != 6) )
          throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 6 gamma line values." );
        break;
    }//switch( version )
    
    
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
    if( !SpecUtils::safe_get_line( strm, line) )
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
    const bool splitres = SpecUtils::split_to_floats( line.c_str(), line.size(), linevaleus );
    
    if( !splitres )
      throw runtime_error( "GadrasGamFile::parse_data(): Not all characters on a gamma group line were numbers." );
    
    switch( version )
    {
      case GamFileVersion::Version_1:
      case GamFileVersion::Version_2_0:
      case GamFileVersion::NumGamFileVersions:
        if( linevaleus.size() != 2 )
          throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 2 gamma group values." );
        break;
        
      case GamFileVersion::Version_2_1:
        if( linevaleus.size() != 6 )
          throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 6 gamma group values." );
        break;
    }//switch( version )
    
    
    m_photon_group_boundries[i] = linevaleus[0];
    m_photon_group_flux[i] = linevaleus[1];
  }//for( int i = 0; i < photon_groups; ++i )
  
  
  
  m_neutron_group_boundries.resize( neutron_groups );
  m_neutron_group_flux.resize( neutron_groups );
  for( int i = 0; i < neutron_groups; ++i )
  {
    if( !SpecUtils::safe_get_line( strm, line) )
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
    const bool splitres = SpecUtils::split_to_floats( line.c_str(), line.size(), linevaleus );
    if( !splitres || (linevaleus.size() != 2) )
      throw runtime_error( "GadrasGamFile::parse_data(): Failed to parse 2 neutron group values." );
    
    m_neutron_group_boundries[i] = 0.001*linevaleus[0];
    m_neutron_group_flux[i] = linevaleus[1];
  }//for( int i = 0; i < neutron_groups; ++i )
  
  m_version = version;
  
  /*  File will then possible have:
   !
   0.767265975475311 ! k-eff
   4.29675034427086 ! multiplication
   */
  
  cerr << "photon_lines=" << photon_lines << ", photon_groups=" << photon_groups << ", neutron_groups=" << neutron_groups << endl;
}//void parse_data( std::istream &strm )


void GadrasGamFile::write_gam( std::ostream &output, const GadrasGamFile::GamFileVersion version )
{
  if( version != GadrasGamFile::GamFileVersion::Version_2_0 )
    throw runtime_error( "GadrasGamFile::write_gam currently only supports writing Version1 gam files" );
  
  if( m_photon_lines_energy.size() != m_photon_lines_flux.size() )
    throw runtime_error( "Photon lines energie has different number of entries than fluxes" );
  
  if( m_photon_group_boundries.size() != m_photon_group_flux.size() )
    throw runtime_error( "Photon group energie has different number of entries than fluxes" );
  
  if( m_neutron_group_boundries.size() != m_neutron_group_flux.size() )
    throw runtime_error( "Neutron group energie has different number of entries than fluxes" );
  
  if( m_photon_lines_energy.empty()
     && m_photon_group_boundries.empty()
     && m_neutron_group_boundries.empty() )
    throw std::runtime_error( "No gamma or neutron information to write." );
  
  
  const char * const endline = "\r\n";
  
  /*
   output << "0 0 ! NewFormat, ModelGeometry" << endline;
   output << m_photon_lines_energy.size()
   << " " << m_photon_group_boundries.size()
   << " " << m_neutron_group_boundries.size()
   << " ! photon lines, photon groups, neutron groups "
   << endline;
   */
  
  output << "GamFileVersion = 2.0" << endline
  << "DataType = SnTransport 18.6.6.0" << endline
  << "Geometry = Spherical" << endline
  << "LargestDimension = 3.5E+01" << endline
  << "Data = " << endline;
  
  const size_t numlines = m_photon_lines_energy.size();
  const size_t num_photon_groups = m_photon_group_boundries.empty() ? 0u : (m_photon_group_boundries.size() - 1);
  const size_t num_neutron_groups = m_neutron_group_boundries.empty() ? 0u : (m_neutron_group_boundries.size() - 1);
  
  output << std::right << std::setw(11)
  << numlines
  << " "
  << std::right << std::setw(11)
  << num_photon_groups
  << " "
  << std::right << std::setw(11)
  << num_neutron_groups
  << " ! photon lines, photon groups, neutron groups "
  << endline;
  
  for( size_t i = 0; i < m_photon_lines_energy.size(); ++i )
  {
    output << std::right << std::setw(11) << std::fixed << std::setprecision(3)
    << m_photon_lines_energy[i]
    << " " << std::scientific << std::setprecision(3) << std::setw(11)
    << m_photon_lines_flux[i];
    
    const float an = (i < m_photon_lines_an.size()) ? m_photon_lines_an[i] : 0.0;
    const float ad = (i < m_photon_lines_ad.size()) ? m_photon_lines_ad[i] : 0.0;
    
    output << " " << std::right << std::fixed << std::setprecision(3)
    << an
    << " " << std::right << std::fixed << std::setprecision(3)
    << ad;
    
    
    if( !i )
      output << " ! photon lines";
    //  output << " ! photon groups (lower bound (keV), intensity(g/s)) ";
    
    output << endline;
  }//for( loop over photon lines )
  
  
  for( size_t i = 0; i < m_photon_group_flux.size(); ++i )
  {
    output << std::right << std::setw(11) << std::fixed << std::setprecision(3)
    << m_photon_group_boundries[i]
    << " " << std::scientific << std::setprecision(3)
    << m_photon_group_flux[i];
    
    if( !i )
      output << " ! photon groups";
    //  output << " ! photon groups (lower bound (keV), intensity(g/s)) ";
    output << endline;
  }
  
  for( size_t i = 0; i < m_neutron_group_boundries.size(); ++i )
  {
    output << std::right << std::setw(11) << std::scientific << std::setprecision(4)
    << 1000.0*m_neutron_group_boundries[i]
    << " " << std::scientific << std::setprecision(3)
    << m_neutron_group_flux[i];
    if( !i )
      output << " ! neutron groups (upper bound (eV), intensity (n/s))";
    output << endline;
  }
  
  // I have not checked if these next lines are actually needed.
  output << "!" << endline;
  output << "0.5 ! k-eff" << endline;
  output << "2.2 ! multiplication" << endline;
  output << endline;
}//void write_gam( std::ofstream &output, const GamFileVersion version )
