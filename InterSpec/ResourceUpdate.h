#ifndef ResourceUpdate_h
#define ResourceUpdate_h
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

/**
 This is only a first hack-in attempt at getting things working - I expect the
 file format and XML schema defined to be somewhat stable, but imrpovements to
 the code may take a little while.
 
 Run after DataBaseVersionUpgrade::checkAndUpgradeVersion();
 
 Allows updating
 SerialToDetectorModel::set_detector_model_input_csv( "data/OUO_detective_serial_to_model.csv" );
 DecayDataBaseServer::setDecayXmlFile( const std::string &path_and_file )
 MassAttenuation::set_data_directory( const std::string &dir );
 
 
 DataBaseUtils::preferenceDatabaseFile();
 */

namespace ResourceUpdate
{
  /** Directory to expand uplaoded resources to for saving across sessions.
      Must be set before anyother calls into the ResourceUpdate.
   
      Throws exception if directory doesnt exist (and cant be created) or user
      cant access it, or an empty directory is specified.
   */
  void setUserDataDirectory( const std::string &dir );
  
  /** Returns the user data directory.
   
      Throws exception if directory hasnt been set.
   */
  std::string userDataDirectory();
  
  
  /**
   
   DecayDataBaseServer::setDecayXmlFile(...)
   MassAttenuation::set_data_directory(...)  (make sure there is a 'em_xs_data' directory in the directory specified, which in turn has all th info)
   
   Will throw exception on error.
   */
  void setupGlobalPrefsFromDb();
  
  enum class ResourceUpdateType
  {
    GadrasDetectorResponseFunction,
    RelativeEfficiencyDetectorResponseFunctionsCsv,
    SandiaDecayXml,
    SandiaReactionXml,
    DetectorSerialToModelInfo,
    ColorThemeJson,
    GammaCrossSectionData,  //MassAttenuation::test_data_directory_validity( const std::string &dir )
    FileQueryFields,
    InvalidResourceUpdateType
  };//enum class ResourceUpdateType
  
  /** Represents one item inside the manifest. */
  struct ResourceUpdateInstance
  {
    bool valid;
    ResourceUpdateType type;
    std::string display_name;
    std::string message;
    std::string path_in_iruz;  //iruz = InterSpec Resource Update Zip
  };//struct ResourceUpdateInstance
  
  
  struct ResourceUpdateCollection
  {
    ResourceUpdateCollection();
    ~ResourceUpdateCollection();
    
    void open_iruz_file( const std::string &filename );  //iruz = InterSpec Resource Update Zip
    
    bool valid;
    std::string errors_opening;
    
    std::string author;
    std::string author_message;
    //TODO: std::string sha256; //to present to user for verification, could use somethign like: https://github.com/okdshin/PicoSHA2/blob/master/picosha2.h
    
    std::string original_iruz_file_location;
    std::string temp_dir_location;
    
    std::vector<ResourceUpdateInstance> resource_instances;
  };//struct ResourceUpdateCollection
  
}//namespace ResourceUpdate


#endif  //ResourceUpdate_h
