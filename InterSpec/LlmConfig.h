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

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

/** Configuration settings for LLM interface and MCP server.
 
 This class handles loading and saving LLM configuration from XML files.
 It follows the InterSpec pattern of looking first in writableDataDirectory()
 then falling back to staticDataDirectory().
 */
class LlmConfig {
public:
  /** LLM API settings */
  struct LlmApi {
    std::string apiEndpoint = "https://api.openai.com/v1/chat/completions";
    std::string bearerToken = "";
    std::string model = "gpt-4";
    int maxTokens = 4000;
    int contextLengthLimit = 128000;
    std::string systemPrompt = "You are an expert assistant for InterSpec, a gamma-ray spectrum analysis application. You can help users analyze spectra, identify peaks, fit nuclides, and understand results. Use the available tools to access spectrum data and perform analysis operations.";
  };
  
  /** MCP server settings */
  struct McpServer {
    int port = 8081;
    std::string bearerToken = "";
    bool enabled = true;
  };
  
  /** UI settings */
  struct Interface {
    bool defaultVisible = false;
    int panelWidth = 400;
  };
  
public:
  LlmApi llmApi;
  McpServer mcpServer;
  Interface interface;
  
  /** Load configuration from XML files with fallback logic.
   
   Tries to load from:
   1. InterSpec::writableDataDirectory() + "/llm_config.xml" 
   2. InterSpec::staticDataDirectory() + "/llm_config.xml"
   3. Uses hardcoded defaults if neither exists
   
   @return LlmConfig instance with loaded settings
   */
  static LlmConfig load();
  
  /** Save configuration to writable data directory.
   
   Saves to InterSpec::writableDataDirectory() + "/llm_config.xml"
   
   @param config The configuration to save
   @return true if save was successful, false otherwise
   */
  static bool save(const LlmConfig& config);
  
  /** Load configuration from a specific XML file.
   
   @param filename Path to the XML configuration file
   @return LlmConfig instance with loaded settings
   @throws std::runtime_error if file cannot be parsed
   */
  static LlmConfig loadFromFile(const std::string& filename);
  
  /** Save configuration to a specific XML file.
   
   @param config The configuration to save
   @param filename Path where to save the XML configuration file
   @return true if save was successful, false otherwise
   */
  static bool saveToFile(const LlmConfig& config, const std::string& filename);

private:
  /** Get the user-writable config file path */
  static std::string getUserConfigPath();
  
  /** Get the default config file path */
  static std::string getDefaultConfigPath();
};

#endif // LLM_CONFIG_H