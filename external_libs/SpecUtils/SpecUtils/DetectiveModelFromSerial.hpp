/* SpecUtils: a library to parse, save, and manipulate gamma spectrum data files.
 
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


#if( PRINT_DETECTIVE_DB_TO_PDF )
#define source(str) ,str
#else
//If we dont want the source information, we will leave it out to reduce bloat
#define source(str)
#endif



namespace
{
//This file implements the detective_model_from_serial() function that is
//  relied on in MeasurementInfo::load_from_binary_spc(...).  
//  This publically released version has had all known serial numbers scrubbed,
//  but still tries to guess model.
  
  
  struct det_info
  {
    int32_t serial;
    DetectorType type;
#if( PRINT_DETECTIVE_DB_TO_PDF )
    const char *description;
    static bool less_than(const det_info &lhs, const det_info &rhs){ return lhs.serial < rhs.serial; }
#endif
  };//struct det_info


  const det_info ns_detectivetypes[] =
  {
    {12345,kDetectiveExDetector source("Example serial number of 12345 mapping to a Detective-EX") },
    {6789,kDetectiveEx100Detector source("Example serial number of 6789 mapping to a Detective-EX100") }
  };
  
  static const size_t num_detectivetypes = sizeof(ns_detectivetypes) / sizeof(ns_detectivetypes[0]);
  
#if( PRINT_DETECTIVE_DB_TO_PDF )
  void create_latex_detective_info( std::ostream &strm )
  {
    static_assert( 0, "PRINT_DETECTIVE_DB_TO_PDF code removed for public release" );
    
  }//void create_latex_info( std::ostream &strm )
#endif
  
//Returns kDetectiveExDetector, kDetectiveEx100Detector,
//  or kMicroDetectiveDetector if the serial number cooresponds to one we know
//  the model of, or if it can be guessed from the serial number
// Returns kUnknownDetector if it cant be guessed.
DetectorType detective_model_from_serial( const std::string &instrument_id )
{
  //Now find runs of digits, so we can search out the numerical serial
  //  number
  std::vector<std::string> serialdigits;
  for( size_t i = 0; i < instrument_id.size(); ++i )
  {
    if( !isdigit(instrument_id[i]) )
      continue;
    size_t j = 1;
    while( isdigit(instrument_id[i+j]) && (i+j)< instrument_id.size() )
      ++j;
    const std::string strval = instrument_id.substr(i,j);
    if( j > 2 && strval != "100" )
      serialdigits.push_back( instrument_id.substr(i,j) );
    i += (j-1);
  }//for( size_t i = 0; i < instrument_id.size(); ++i )
  
  
  for( const string &serialstr : serialdigits )
  {
    int32_t val = 0;
    try
    {
      val = static_cast<int32_t>( std::stoi(serialstr) );
    }catch(...){}
    
    if( serialdigits.size() >= 2 && val < 100 )
      continue;
    
    for( size_t i = 0; i < num_detectivetypes; ++i )
    {
      if( ns_detectivetypes[i].serial == val )
        return DetectorType( ns_detectivetypes[i].type );
    }
  }
  
#if( PERFORM_DEVELOPER_CHECKS || INCLUDE_ANALYSIS_TEST_SUITE )
  char msg[256];
  snprintf( msg, sizeof(msg), "Detective model for detector with serial"
           " number '%s' is not known", instrument_id.c_str() );
  passMessage( msg, "", 3 );
#endif
  
  //The case independance does not appear to be necassarry, but JIC
  if( UtilityFunctions::icontains( instrument_id, "EX100" )
     || UtilityFunctions::icontains( instrument_id, "EX 100" ) )
  {
    return kDetectiveEx100Detector;
  }
  
  if( UtilityFunctions::icontains( instrument_id, "Micro" )
     || UtilityFunctions::icontains( instrument_id, "uDet" )
     || UtilityFunctions::icontains( instrument_id, "HX" )
     || UtilityFunctions::icontains( instrument_id, "uDX") )
  {
    return kMicroDetectiveDetector;
  }
  
  
  // Make a not very educated guess based off of limited examples
  for( const string &serialstr : serialdigits )
  {
    int32_t val = 0;
    try
    {
      val = static_cast<int32_t>( std::stoi(serialstr) );
    }catch(...){}
    
    if( val >= 500 && val < 4000 )
      return kDetectiveExDetector;
    
    if( val >= 4000 && val < 5000 )
      return kDetectiveEx100Detector;
  }//for( const string &serialstr : serialdigits )
  
  return kUnknownDetector;
}//DetectorType detective_model_from_serial( const string &instrument_id )
}//namespace Detectiv
