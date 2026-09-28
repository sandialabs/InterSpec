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

// Native driver for the light app: reads JSON requests (one per line) from a file or stdin, and
//  writes each response on its own line - the same dispatcher the browser uses, for debugging.
//  Lines starting with '#' are ignored.
//  Usage: light_cli [requests.jsonl] [ref_lines.json]
//  If a reference-line library is given, it is loaded (like the page does) before the requests.

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

extern "C" const char *light_call( const char *request );

int main( int argc, char **argv )
{
  std::ifstream file;
  if( argc > 1 )
  {
    file.open( argv[1] );
    if( !file )
    {
      std::cerr << "Could not open '" << argv[1] << "'" << std::endl;
      return 1;
    }
  }
  std::istream &input = (argc > 1) ? static_cast<std::istream &>(file) : std::cin;

  int nerrors = 0;
  if( argc > 2 )
  {
    std::ifstream lib( argv[2] );
    std::stringstream strm;
    strm << "{\"method\":\"setRefLibrary\",\"params\":{\"lib\":" << lib.rdbuf() << "}}";
    const std::string response = light_call( strm.str().c_str() );
    nerrors += (response.compare( 0, 9, "{\"error\":" ) == 0);
    std::cout << response << std::endl;
  }

  std::string line;
  while( std::getline( input, line ) )
  {
    if( line.empty() || (line[0] == '#') )
      continue;
    const std::string response = light_call( line.c_str() );
    // Test scripts mark requests that should fail with "expectError" (see tests/check_expectations.py)
    nerrors += ((response.compare( 0, 9, "{\"error\":" ) == 0) && (line.find( "\"expectError\"" ) == std::string::npos));
    std::cout << response << std::endl;
  }

  return nerrors ? 2 : 0;
}//main(...)
