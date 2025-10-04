#ifndef LLM_CONFIG_H
#define LLM_CONFIG_H
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
#include <memory>

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

/** Configuration settings for LLM interface and MCP server.
 
 This class handles loading and saving LLM configuration from XML files.
 It follows the InterSpec pattern of looking first in writableDataDirectory()
 then falling back to staticDataDirectory().
 */
class LlmConfig
{
  static const int sm_xmlSerializationVersion = 0;
  
public:
  
  /** LLM API settings - this is for the "LLM Assistant" tool. */
  struct LlmApi
  {
    bool enabled = false;
    std::string apiEndpoint;    // Example: "https://api.openai.com/v1/chat/completions"
    std::string bearerToken;    //
    std::string model;          // Example: "gpt-4"
    int maxTokens = 0;          // Example: 4000
    int contextLengthLimit = 0; // Example: 128000;
    std::string systemPrompt;   // Example: "You are an expert assistant for InterSpec...";
  };//struct LlmApi
  
  
  /** MCP server settings */
  struct McpServer
  {
    static const std::string sm_invalid_bearer_token;
    
    bool enabled = false;
    
#if( MCP_ENABLE_AUTH )
    /** The bearer token that requestors must supply.
     
     `sm_invalid_bearer_token` is a canary, in principle it should never be set to this value if `enabled` is true,
     but just to be sure,  then you should refuse to create a server if set to this value.
     */
    std::string bearerToken = McpServer::sm_invalid_bearer_token;
#endif
  };//struct McpServer
  
public:
  LlmApi llmApi;
  McpServer mcpServer;
  
  /** Load configuration from XML files with fallback logic.
   
   Tries to load from:
   1. InterSpec::writableDataDirectory() + "/llm_config.xml" 
   2. InterSpec::staticDataDirectory() + "/llm_config.xml"
   3. Uses hardcoded defaults if neither exists
   
   If either config file exists, but is invalid, throws exception.
   
   @return LlmConfig instance with loaded settings, or nullptr if no config file is present.
   */
  static std::shared_ptr<LlmConfig> load();
  
  /** Load configuration from a specific XML file.
   
   @param filename Path to the XML configuration file
   @return LlmConfig instance with loaded settings
   @throws std::runtime_error if file cannot be parsed
   */
  static std::shared_ptr<LlmConfig> loadFromFile( const std::string &filename );
  
  /** Save configuration to a specific XML file.
   
   @param config The configuration to save
   @param filename Path where to save the XML configuration file
   @return true if save was successful, false otherwise
   
   Suggest saving to InterSpec::writableDataDirectory() + "/llm_config.xml" - e.g., what `LlmConfig::getUserConfigPath()` will return.
   */
  static bool saveToFile( const LlmConfig &config, const std::string &filename );


  /** Get the user-writable config file path - e.g., `InterSpec::writableDataDirectory() + "/llm_config.xml"` */
  static std::string getUserConfigPath();
  
  /** Get the default config file path - e.g., `InterSpec::staticDataDirectory() + "/llm_config.xml"` */
  static std::string getDefaultConfigPath();
};

#endif // LLM_CONFIG_H
