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

#include <mutex>
#include <string>

#include <Wt/Dbo/Session>
#include <Wt/Dbo/Exception>
#include <Wt/Dbo/Transaction>
#include <Wt/Dbo/SqlConnection>

#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/ResourceUpdate.h"
#include "SpecUtils/UtilityFunctions.h"

using namespace std;

namespace
{
  std::string ns_user_data_dir = "";
  std::mutex ns_user_data_dir_mutex;
}//namespace


namespace ResourceUpdate
{
  void setUserDataDirectory( const std::string &dir )
  {
    std::unique_lock<mutex> lock( ns_user_data_dir_mutex );
    
    if( !UtilityFunctions::is_directory(dir) )
    {
      const auto parent = UtilityFunctions::parent_path(dir);
      if( !UtilityFunctions::is_directory(parent) )
        throw runtime_error( "setUserDataDirectory(): Not a valid directory: '" + dir + "'" );
      if( UtilityFunctions::create_directory(dir) == 0 )
        throw runtime_error( "setUserDataDirectory(): Unable to create directory: '" + dir + "'" );
    }
    
    ns_user_data_dir = dir;
  }//setUserDataDirectory(...)
  
  
  std::string userDataDirectory()
  {
    std::unique_lock<mutex> lock( ns_user_data_dir_mutex );
    if( ns_user_data_dir.empty() )
      throw runtime_error( "userDataDirectory(): user directory hasnt been set yet." );
    return ns_user_data_dir;
  }
  
  
  void setupGlobalPrefsFromDb()
  {
    /*
     std::unique_ptr<DataBaseUtils::DbSession> session;
     std::unique_ptr<DataBaseUtils::DbTransaction> transaction;
     
     try
     {
     session.reset( new DataBaseUtils::DbSession() );
     transaction.reset( new DataBaseUtils::DbTransaction(*session) );
     }catch(...)
     {
     throw runtime_error( "setupGlobalPrefsFromDb(): unable to open database session" );
     }//try / catch
     
     
     //SandiaDecayXml
     try
     {
     auto result = session->session()->find<InterSpecGlobalSetting>()
     .where( "settingCatagory = 'SandiaDecayXml'" )
     .orderBy("id desc").limit( 1 ).resultList();
     if( result.size() > 0 )
     {
     
     }
     }catch(...)
     {
     throw runtime_error( "Error retrieving SandiaDecayXml gloabal session" );
     }//try / catch load SAndiaDecayXml
     */
    
    //DetectorSerialToModelInfo,
    
    
    //GadrasDetectorResponseFunction,
    //RelativeEfficiencyDetectorResponseFunctionsCsv,
    //SandiaReactionXml,
    //GammaCrossSectionData, FileQueryFields,  //MassAttenuation::test_data_directory_validity( const std::string &dir )
    //InvalidResourceUpdateType
    
    
    
  }//void setupGlobalPrefsFromDb()
  
  
  ResourceUpdateCollection::ResourceUpdateCollection()
  {
    //TODO
  }
  
  
  ResourceUpdateCollection::~ResourceUpdateCollection()
  {
    //TODO
  }
  
}
