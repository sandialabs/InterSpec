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
#include <cstdlib>
#include <iostream>

#include <boost/program_options.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/InterSpecServer.h"
#include "InterSpec/BatchCommandLine.h"
#include "InterSpec/DecayDataBaseServer.h"

#if( BUILD_AS_OSX_APP || defined(__APPLE__) )
#include "target/osx/macOsUtils.h"
#endif

static_assert( USE_BATCH_CLI_TOOLS, "You must have USE_BATCH_CLI_TOOLS enabled to build this code" );


#ifdef _WIN32
int wmain( int argc, wchar_t *wargv[] )
{
  char **argv;
  AppUtils::getUtf8Args( argc, wargv, argv );
#else
int main( int argc, char *argv[] )
{
#endif
  std::string user_data_dir, docroot;
  bool batch_peak_fit = false, batch_act_fit = false;
#if( USE_REL_ACT_TOOL )
  bool batch_iso_from_nucs = false;
#endif
  
  namespace po = boost::program_options;
  
  unsigned term_width = AppUtils::terminal_width();
  unsigned min_description_length = ((term_width < 80u) ? term_width/2u : 40u);
  
  po::options_description cl_desc("Allowed options", term_width, min_description_length);
  cl_desc.add_options()
  ("help,h",  "produce this help message")
  ("userdatadir", po::value<std::string>(&user_data_dir),
   "The directory to store user data to, or to look in for custom user data (detector efficiency functions, etc)."
   )
  ("docroot", po::value<std::string>(&docroot),
   "The directory that contains the 'InterSpec_resources' and 'data' directories.\n"
   )
  ("batch-peak-fit", po::value<bool>(&batch_peak_fit)->implicit_value(true)->default_value(false),
   "Batch-fit peaks.\n"
   "\tUse '--batch-peak-fit --help' to see available options."
   )
  ("batch-act-fit", po::value<bool>(&batch_act_fit)->implicit_value(true)->default_value(false),
   "Batch shielding/source fit.\n"
   "\tUse '--batch-act-fit --help' to see available options."
   )
#if( USE_REL_ACT_TOOL )
  ("batch-iso-from-nucs", po::value<bool>(&batch_iso_from_nucs)->implicit_value(true)->default_value(false),
   "Batch isotopics-by-nuclides (RelActCalcAuto) analysis.\n"
   "\tUse '--batch-iso-from-nucs --help' to see available options."
   )
#endif
  ;
  
  po::variables_map cl_vm;
  try
  {
    po::parsed_options parsed_opts
      = po::command_line_parser(argc,argv)
       .allow_unregistered()
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

#if( USE_REL_ACT_TOOL )
  const bool is_batch = (batch_peak_fit || batch_act_fit || batch_iso_from_nucs);
#else
  const bool is_batch = (batch_peak_fit || batch_act_fit);
#endif
  
  if( cl_vm.count("help") )
  {
    if( is_batch )
    {
      return BatchCommandLine::run_batch_command( argc, argv );
    }else
    {
      std::cout << "Executable compiled on " << __DATE__ << std::endl;
      std::cout << "Available command-line options for starting the InterSpec web-server are:\n";
      std::cout << cl_desc << std::endl;
      return EXIT_SUCCESS;
    }
  }//if( cl_vm.count("help") )
  
  if( !is_batch )
  {
#if( USE_REL_ACT_TOOL )
    std::cerr << "You must specify to use 'batch-peak-fit', 'batch-act-fit',"
                 " or 'batch-iso-from-nucs'.\n"
              << "With additional available options of:\n";
#else
    std::cerr << "You must specify to use either 'batch-peak-fit', or 'batch-act-fit'.\n"
              << "With additional available options of:\n";
#endif
    std::cout << cl_desc << std::endl;
    
    return EXIT_FAILURE;
  }

  if( docroot.empty() )
  {
#if( BUILD_AS_OSX_APP )
    docroot = macOsUtils::static_data_base_dir();
#else
    // I cant get MSVC to set CWD to anywhere besides InterSpec/out/build/x64-Debug/,
    //  so we'll look for our resources up to three levels up.
    //  However, we def dont want to do this for anything other than localhost development
    std::string targetfile = "InterSpec_resources/InterSpec.css";
    if( !AppUtils::locate_file( targetfile, false, 3, false ) )
    {
      std::cerr << "Unable to find base directory that contains 'InterSpec_resources' directory."
      << std::endl;
      return -24;
    }
    
    docroot = SpecUtils::parent_path( SpecUtils::parent_path( targetfile ) );
    if( docroot.empty() && !targetfile.empty() ) //running from same dir as exe.
      docroot = ".";
#endif
  }//if( docroot.empty() )

  if( !SpecUtils::is_directory(docroot) )
  {
    std::cerr << "Sorry, the directory '" << docroot << "' is not a valid directory.\n"
    "Please specify the 'docroot' argument with the path pointing to where InterSpecs static"
    " data is located (i.e., the 'InterSpec_resources' and 'data' directories)." << std::endl;
    
    return -5;
  }
  
  const std::string datadir = SpecUtils::append_path( docroot, "data" );
  if( !SpecUtils::is_directory(datadir) )
  {
    std::cerr << "There is not a 'data' directory in docroot ('" << docroot << "') as required"
    << std::endl;
    return -6;
  }
  
  try
  {
    InterSpec::setStaticDataDirectory( datadir );
  }catch( std::exception &e )
  {
    std::cerr << "Failed to set static data directory: " << e.what() << std::endl;
    return -7;
  }
  
  
  if( user_data_dir.empty() )
  {
#if( BUILD_AS_OSX_APP || defined(__APPLE__) )
    user_data_dir = macOsUtils::user_data_dir();
    if( !SpecUtils::is_directory(user_data_dir) )
    {
      user_data_dir = "";
      std::cerr << "Warning: not setting user data directory, to what the GUI version of InterSpec"
      " uses - use the '--userdatadir' option to set this if you would like to share detectors,"
      " or other files." << std::endl;
    }
#elif defined(_WIN32)
    user_data_dir = AppUtils::user_data_dir();
#else
    // Linux: try directories in priority order, falling through to the next
    // candidate if the previous one does not exist.  The --userdatadir flag
    // overrides this when sharing with a specific install.
    //
    // Priority:
    //   1. ~/.local/share/InterSpec        - wxWidgets GUI default
    //                                        (wxStandardPaths::GetUserDataDir,
    //                                        XDG-style, app name "InterSpec").
    //   2. ~/.config/InterSpecAddOn        - target/electron GUI default
    //                                        (Electron's app.getPath('userData')).
    //   3. ~/.local/share/sandia.InterSpec - plain XDG with "sandia.InterSpec".
    {
      const char * const xdg_data   = std::getenv( "XDG_DATA_HOME" );
      const char * const xdg_config = std::getenv( "XDG_CONFIG_HOME" );
      const char * const home       = std::getenv( "HOME" );

      const std::string data_base = (xdg_data && xdg_data[0])
        ? std::string( xdg_data )
        : ((home && home[0]) ? SpecUtils::append_path( home, ".local/share" )
                             : std::string());
      const std::string config_base = (xdg_config && xdg_config[0])
        ? std::string( xdg_config )
        : ((home && home[0]) ? SpecUtils::append_path( home, ".config" )
                             : std::string());

      const std::string candidates[] = {
        data_base.empty()   ? std::string() : SpecUtils::append_path( data_base,   "InterSpec" ),
        config_base.empty() ? std::string() : SpecUtils::append_path( config_base, "InterSpecAddOn" ),
        data_base.empty()   ? std::string() : SpecUtils::append_path( data_base,   "sandia.InterSpec" ),
      };

      for( const std::string &c : candidates )
      {
        if( !c.empty() && SpecUtils::is_directory(c) )
        {
          user_data_dir = c;
          break;
        }
      }

      if( user_data_dir.empty() )
      {
        std::cerr << "Warning: not setting user data directory, to what the GUI version of InterSpec"
        " uses - use the '--userdatadir' option to set this if you would like to share detectors,"
        " or other files." << std::endl;
      }
    }
#endif
  }//if( user_data_dir.empty() )
  
  // We need to somehow set the static data directory
  
  if( !user_data_dir.empty() )
  {
    try
    {
      InterSpec::setWritableDataDirectory( user_data_dir );
    }catch( std::exception &e )
    {
      std::cerr << "Warning: Failed to set user data directory: " << e.what() << std::endl
      << "Use the '--userdatadir' option to set this to the same one the InterSpec GUI app"
      << " uses, if you would like to share detectors, or other files." << std::endl;
      //return -8;
    }
    
    try
    {
      const std::string user_decay = SpecUtils::append_path(user_data_dir, "sandia.decay.xml");
      if( SpecUtils::is_file(user_decay) )
        DecayDataBaseServer::setDecayXmlFile(user_decay);  //overwrites what `InterSpec::setStaticDataDirectory` set
    }catch (std::exception& e)
    {
      std::cerr << "Failed in setting user sandia.decay.xml file: " << e.what() << std::endl;
      return -9;
    }
    
    try
    {
      const std::string user_reaction = SpecUtils::append_path(user_data_dir, "sandia.reactiongamma.xml");
      if( SpecUtils::is_file(user_reaction) )
        ReactionGammaServer::set_xml_file_location(user_reaction); //overwrites what `InterSpec::setStaticDataDirectory` set
    }catch (std::exception& e)
    {
      std::cerr << "Failed in setting user sandia.reactiongamma.xml file: " << e.what() << std::endl;
      return -10;
    }
  }//if( !user_data_dir.empty() )

  
#if( PERFORM_DEVELOPER_CHECKS )
  std::cout << "Developer tests are being performed" << std::endl;
#endif

  
  return BatchCommandLine::run_batch_command( argc, argv );
}//int main( int argc, const char * argv[] )


