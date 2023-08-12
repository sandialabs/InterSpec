//  Created by wcjohns on 20110322.
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

/**
 This file is for creating a InterSpec web-server; either localhost (primarily for development
 purposes), or as a web-server (has only been superficially tested inside docker containers).
 
 To build stand-alone "native" apps, see the sub-directories in the "target" directory.
 */

#include "InterSpec_config.h"

#include <string>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecServer.h"

#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
#include "testing/developcode.h"
#endif

#if( USE_BATCH_TOOLS )
#include "InterSpec/BatchCommandLine.h"
#endif

//#include "InterSpec/QRSpectrum.h"

int main( int argc, char **argv )
{
#ifdef _WIN32
  AppUtils::getUtf8Args( argc, argv );
#endif
  
  //return QRSpectrum::dev_code();
  
#if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  return developcode::run_development_code();
#endif
  
  std::cout << std::showbase << std::hex << "Running with Wt version "
  << WT_VERSION << std::dec << ", from executable compiled on "
  << __DATE__ << std::endl;
  
#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "Developer tests are being performed" << std::endl;
#endif
  
  std::cout << std::endl;
  
  int server_port_num;
  std::string docroot, wt_config, user_data_dir;
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
  std::string http_address = "127.0.0.1";
  static_assert( !BUILD_AS_LOCAL_SERVER, "Web and local server should not both be enabled");
#endif
  
  namespace po = boost::program_options;
  
  unsigned term_width = AppUtils::terminal_width();
  unsigned min_description_length = ((term_width < 80u) ? term_width/2u : 40u);
  
  po::options_description cl_desc("Allowed options", term_width, min_description_length);
  cl_desc.add_options()
  ("help,h",  "produce this help message")
  ("http-port", po::value<int>(&server_port_num)->default_value(8080),
   "The HTTP port to bind the web-server too.  Ports below 1024 are not recommended, and require elevated privileges.")
#if( BUILD_FOR_WEB_DEPLOYMENT )
  ( "http-address", po::value<std::string>(&http_address),
   "The network HTTP address to bind the web-server too; '127.0.0.1' is localhost, while '0.0.0.0' will serve the web-app to the external network."
   )
#endif
  ("userdatadir", po::value<std::string>(&user_data_dir),
   "The directory to store user data to, or to look in for custom user data (serial_to_model.csv, etc)."
   )
  ("config", po::value<std::string>(&wt_config),
   "The Wt config XML file to use."
   )
  ("docroot", po::value<std::string>(&docroot),
   "The directory that contains the 'InterSpec_resources' and 'data' directories.\n"
   "All files in the docroot directory, and its subdirectories are available via HTTP.\n"
#if( !BUILD_FOR_WEB_DEPLOYMENT )
   "Defaults to current working directory."
#endif
   )
  ("static-data-dir", "The static data directory (e.g., 'data' dir that holds cross-sections, "
   "nuclear-data, etc) to use.  If not specified, uses 'data' in the `docroot` directory."
   )
#if( USE_BATCH_TOOLS )
  ("batch-peak-fit", "Batch-fit peaks.")
  ("batch-act-fit", "Batch shielding/source fit.")
#endif
  ;
  
  po::variables_map cl_vm;
  try
  {
    po::parsed_options parsed_opts
      = po::command_line_parser(argc,argv)
#if( USE_BATCH_TOOLS )
       .allow_unregistered()
#endif
       .options(cl_desc)
       .run();
    
    po::store( parsed_opts, cl_vm );
    po::notify( cl_vm );
  }catch( std::exception &e )
  {
    std::cerr << "Command line argument error: " << e.what() << std::endl << std::endl;
    std::cout << cl_desc << std::endl;
    return 1;
  }//try catch
  
  
#if( USE_BATCH_TOOLS )
  const bool is_batch = (cl_vm.count("batch-peak-fit") || cl_vm.count("batch-act-fit"));
#else
  const bool is_batch = false;
#endif
  
  if( cl_vm.count("help") && !is_batch )
  {
    std::cout << "Available command-line options for starting the InterSpec web-server are:\n";
    std::cout << cl_desc << std::endl;
    return 0;
  }//if( cl_vm.count("help") )
  
  
#if( BUILD_FOR_WEB_DEPLOYMENT )
  if( cl_vm.count("config") )
  {
    std::cerr << "You must specify the Wt config file to use (the 'config' option)" << std::endl;
    return -20;
  }
  
  if( cl_vm.count("http-address") )
  {
    std::cerr << "You must specify the network adapter address to bind to"
    << " (the 'http-address' option)." << std::endl;
    return -21;
  }
  
  if( cl_vm.count("docroot") )
  {
    std::cerr << "You must specify the HTTP document root directory to use (the 'docroot' option)" << std::endl;
    return -22;
  }
#endif
  
  if( (server_port_num <= 0) || (server_port_num > 65535) )
  {
    std::cerr << "Invalid server port number: " << server_port_num
    << ", must be between 1 and 65535" << std::endl;
    return -23;
  }
  
  
  if( server_port_num <= 1024 )
  {
#if( !BUILD_FOR_WEB_DEPLOYMENT )
    std::cerr << "Server port number below 1024 not allowed." << std::endl;
    return -23;
#else
    std::cerr << "Warning: using a privileged port is not recommended." << std::endl;
#endif
  }//if( server_port_num <= 1024 )
  
  
#if( _WIN32 && !BUILD_FOR_WEB_DEPLOYMENT && BUILD_AS_LOCAL_SERVER )
  if( docroot.empty() )
  {
    // I cant get MSVC to set CWD to anywhere besides InterSpec/out/build/x64-Debug/,
    //  so we'll look for our resources up to three levels up.
    //  However, we def dont want to do this for anything other than localhost development
    const std::string targetfile = "InterSpec_resources/InterSpec.css";
    std::string pardir = "";
    while( pardir.size() < 9 )
    {
      const std::string testfile = SpecUtils::append_path(pardir, targetfile);
      if( SpecUtils::is_file(testfile) )
      {
        docroot = pardir;
        break;
      }//if( SpecUtils::is_file(testfile) )
      pardir = SpecUtils::append_path(pardir, "..");
    }//while (pardir.size() < 6)
    
    if( !SpecUtils::is_file( SpecUtils::append_path(docroot, targetfile) ) )
    {
      std::cerr << "Unable to find base directory that contains 'InterSpec_resources' directory."
      << std::endl;
      return -24;
    }
  }//if( docroot.empty() )
#endif //_WIN32 && !BUILD_FOR_WEB_DEPLOYMENT
  
  if( user_data_dir.empty() )
  {
#if( BUILD_AS_LOCAL_SERVER )
    // If there is a "user_data" directory in the CWD, we'll set this as the writeable data
    //  directory to simulate desktop app behavior of saving DRFs and similar
    const std::string cwd = SpecUtils::get_working_path();
    const std::string dev_user_data = SpecUtils::append_path( cwd, "user_data" );
    if( SpecUtils::is_directory( dev_user_data ) )
    {
      user_data_dir = dev_user_data;
    }else
    {
      std::cerr << "No '" << dev_user_data << "' - you must specify writeable data directory,"
                << " or there must be a 'user_data' directory in the current working directory."
                << std::endl;
      return -25;
    }
#else
    std::cerr << "You must specify the directory to store user data to (the 'userdatadir' option)."
              << std::endl;
    return -25;
#endif
  }//if( user_data_dir.empty() )

  
#if( !BUILD_FOR_WEB_DEPLOYMENT )
  if( docroot.empty() )
    docroot = ".";
  
  if( wt_config.empty() )
    wt_config = SpecUtils::append_path( docroot, "data/config/wt_config_localweb.xml" );
#endif
  
  if( cl_vm.count("static-data-dir") )
  {
    std::string datadir = cl_vm["static-data-dir"].as<std::string>();
    if( !SpecUtils::is_directory(datadir) )
      datadir = SpecUtils::append_path( docroot, datadir );
    
    if( !SpecUtils::is_directory(datadir) )
    {
      std::cerr << "Specified 'static-data-dir' ('" << cl_vm["static-data-dir"].as<std::string>()
                << "') is not a directory." << std::endl;
      return -26;
    }//if( !SpecUtils::is_directory(datadir) )
    
    // We wont make datadir path absolute, to avoid possible long name hassles on Windows, and I
    //  dont think the code changes current working directory anywhere.
    //datadir = SpecUtils::make_canonical_path(datadir);
    
    InterSpec::setStaticDataDirectory( datadir );
  }//if( cl_vm.count("static-data-dir") )
  
  
#if( USE_BATCH_TOOLS )
  if( is_batch )
    return BatchCommandLine::run_batch_command( argc, argv );
#endif
  
  
  // Start the InterSpec server
  const int rval = InterSpecServer::start_server( argv[0], user_data_dir.c_str(),
                                                 docroot.c_str(),
                                                 wt_config.c_str(),
                                                 static_cast<short int>(server_port_num) );
  if( rval < 0 )
  {
    std::cerr << "Failed to start server, val=" << rval << std::endl;
    return rval;
  }
  
  std::cout << "\nYou can now point your browser to: " << InterSpecServer::urlBeingServedOn()
            << std::endl;
  
  return InterSpecServer::wait_for_shutdown();
}//int main( int argc, const char * argv[] )


