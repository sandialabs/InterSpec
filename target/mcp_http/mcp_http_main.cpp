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

/** The InterSpec multi-session MCP analysis server: MCP clients create sessions, exchange files,
 and use InterSpec (and its LLM assistant) through `/mcp-api`; no session can be opened from a
 browser.  See README.md in this directory.
 */

#include "InterSpec_config.h"

#include <string>
#include <cstdlib>
#include <csignal>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecServer.h"

static_assert( BUILD_AS_MCP_HTTP_SERVER, "This file is only for the MCP HTTP server build" );


int main( int argc, char **argv )
{
#ifdef _WIN32
  AppUtils::getUtf8Args( argc, argv );
#else
  // A write to the pipe of a browser that just died (e.g., killed for memory) must fail with
  //  EPIPE, not kill the server, and every session with it.
  std::signal( SIGPIPE, SIG_IGN );
#endif

  std::cout << "InterSpec MCP analysis server, compiled " << __DATE__ << std::endl;

  int server_port_num = 8078;
  std::string http_address = "127.0.0.1", docroot = ".", wt_config, user_data_dir, llm_config_dir,
              static_data_dir;

  namespace po = boost::program_options;

  const unsigned term_width = AppUtils::terminal_width();
  const unsigned min_description_length = ((term_width < 80u) ? term_width/2u : 40u);

  po::options_description cl_desc( "Allowed options", term_width, min_description_length );
  cl_desc.add_options()
  ("help,h",  "produce this help message")
  ("http-port", po::value<int>(&server_port_num)->default_value(8078), "The HTTP port to serve on.")
  ("http-address", po::value<std::string>(&http_address)->default_value("127.0.0.1"),
   "The network address to bind to; '0.0.0.0' serves on all interfaces (e.g., in a container)."
   "  Clients must present the MCP bearer token from llm_config.xml."
   )
  ("llm-config-dir", po::value<std::string>(&llm_config_dir),
   "Required: the directory with llm_config.xml (and optionally llm_tools_config.xml,"
   " llm_agent_*.xml, mcp_agent_instructions.md, and an llm_knowledge directory).  The"
   " INTERSPEC_LLM_CONFIG_DIR environment variable is used if not specified."
   )
  ("userdatadir", po::value<std::string>(&user_data_dir),
   "Required: a writable directory for InterSpec's per-user data (the sessions' own files go in"
   " the exchange directory set in llm_config.xml)."
   )
  ("docroot", po::value<std::string>(&docroot)->default_value("."),
   "The directory containing the 'InterSpec_resources', 'resources', and 'data' directories."
   )
  ("config,c", po::value<std::string>(&wt_config),
   "The Wt config XML file; defaults to <docroot>/data/config/wt_config_mcp_http.xml."
   )
  ("static-data-dir", po::value<std::string>(&static_data_dir),
   "The static data directory (nuclear data, cross-sections, ...); defaults to <docroot>/data."
   )
  ;

  po::variables_map cl_vm;
  try
  {
    po::store( po::command_line_parser(argc,argv).options(cl_desc).run(), cl_vm );
    po::notify( cl_vm );
  }catch( std::exception &e )
  {
    std::cerr << "Command line argument error: " << e.what() << std::endl << std::endl;
    std::cout << cl_desc << std::endl;
    return 1;
  }//try catch

  if( cl_vm.count("help") )
  {
    std::cout << "Available command-line options for the InterSpec MCP analysis server are:\n"
              << cl_desc << std::endl;
    return 0;
  }

  if( llm_config_dir.empty() )
  {
    const char * const env_dir = std::getenv( "INTERSPEC_LLM_CONFIG_DIR" );
    if( env_dir )
      llm_config_dir = env_dir;
  }

  // The server must run with the configuration it was deployed with, never the inert default in
  //  the static data directory.
  if( llm_config_dir.empty() || !SpecUtils::is_directory( llm_config_dir ) )
  {
    std::cerr << "An LLM config directory must be given with --llm-config-dir (or"
              << " INTERSPEC_LLM_CONFIG_DIR)." << std::endl;
    return -27;
  }
  LlmConfig::setConfigDirectory( llm_config_dir );

  if( user_data_dir.empty() || !SpecUtils::is_directory( user_data_dir ) )
  {
    std::cerr << "A writable user data directory must be given with --userdatadir." << std::endl;
    return -25;
  }

  if( (server_port_num <= 0) || (server_port_num > 65535) )
  {
    std::cerr << "Invalid server port number: " << server_port_num << std::endl;
    return -23;
  }

  if( static_data_dir.empty() )
    static_data_dir = SpecUtils::append_path( docroot, "data" );
  if( !SpecUtils::is_directory( static_data_dir ) )
  {
    std::cerr << "The static data directory ('" << static_data_dir << "') does not exist;"
              << " specify --docroot or --static-data-dir." << std::endl;
    return -26;
  }
  InterSpec::setStaticDataDirectory( static_data_dir );

  if( wt_config.empty() )
    wt_config = SpecUtils::append_path( static_data_dir, "config/wt_config_mcp_http.xml" );

  const int rval = InterSpecServer::start_server( argv[0], user_data_dir.c_str(), docroot.c_str(),
                                                  wt_config.c_str(),
                                                  static_cast<unsigned short int>(server_port_num),
                                                  http_address.c_str() );
  if( rval < 0 )
  {
    std::cerr << "Failed to start server, val=" << rval << std::endl;
    InterSpecServer::killServer();  // it may have started listening before failing
    return rval;
  }

  const int shutdown_status = InterSpecServer::wait_for_shutdown();

  // Ends the sessions (removing their directories, unless configured to keep them), and stops the
  //  browser that hosts them.
  InterSpecServer::killServer();

  return shutdown_status;
}//int main( int argc, char **argv )
