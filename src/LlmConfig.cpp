#include "InterSpec_config.h"
#include "InterSpec/LlmConfig.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/XmlUtils.hpp"



using namespace std;

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

const std::string LlmConfig::McpServer::sm_invalid_bearer_token = "INVALID-BEARER-TOKEN";


std::shared_ptr<LlmConfig> LlmConfig::load()
{
  // Try user config first
  string userPath = getUserConfigPath();
  if( SpecUtils::is_file(userPath) )
    return loadFromFile(userPath);
  
  // Try default config
  string defaultPath = getDefaultConfigPath();
  if( SpecUtils::is_file(defaultPath) )
    return loadFromFile(defaultPath);
  
  cout << "Using hardcoded LLM config defaults" << endl;
  LlmConfig config;
  assert( !config.llmApi.enabled );
  assert( !config.mcpServer.enabled );
  
  config.llmApi.enabled = false;   //just to make sure
  config.mcpServer.enabled = false;
  
  return nullptr;
}//LlmConfig LlmConfig::load()


std::shared_ptr<LlmConfig> LlmConfig::loadFromFile( const std::string &filename )
{
  try
  {
    std::vector<char> xmlContent;
    SpecUtils::load_file_data( filename.c_str(), xmlContent );
    
    rapidxml::xml_document<char> doc;
    const int flags = (rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace);
    
    doc.parse<flags>( &xmlContent[0] );
    
    const rapidxml::xml_node<char> * const root = XML_FIRST_NODE(&doc, "LlmConfig");
    if( !root )
      throw std::runtime_error("Missing LlmConfig root node");
    
    
    // A reminder double check these logics when changing RoiRange::sm_xmlSerializationVersion
    static_assert( LlmConfig::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    
    XmlUtils::check_xml_version( root, LlmConfig::sm_xmlSerializationVersion );
    
    std::shared_ptr<LlmConfig> config = make_shared<LlmConfig>();
    
    {// Begin load LLM API settings
      const rapidxml::xml_node<char> * const llmApi = XML_FIRST_NODE(root, "LlmApi");
      if( !llmApi )
        throw runtime_error( "Missing 'LlmApi' node." );
      
      config->llmApi.enabled = XmlUtils::get_bool_node_value(llmApi, "Enabled");
      config->llmApi.apiEndpoint = XmlUtils::get_string_node_value(llmApi, "ApiEndpoint");
      config->llmApi.bearerToken = XmlUtils::get_string_node_value(llmApi, "BearerToken");
      config->llmApi.model = XmlUtils::get_string_node_value(llmApi, "Model");
      config->llmApi.maxTokens = XmlUtils::get_int_node_value(llmApi, "MaxTokens");
      config->llmApi.contextLengthLimit = XmlUtils::get_int_node_value(llmApi, "ContextLengthLimit");
      config->llmApi.systemPrompt = XmlUtils::get_string_node_value(llmApi, "SystemPrompt");
    }// End load LLM API settings
    
    {// Begin load MCP server settings
      const rapidxml::xml_node<char> * const mcpServer = XML_FIRST_NODE(root, "McpServer");
      if( !mcpServer )
        throw runtime_error( "Missing 'McpServer' node." );
      config->mcpServer.enabled = XmlUtils::get_bool_node_value(mcpServer, "Enabled");
#if( MCP_ENABLE_AUTH )
      config->mcpServer.bearerToken = XmlUtils::get_string_node_value( mcpServer, "BearerToken" );
#endif
    }// End load MCP server settings
    
    return config;
  }catch( const std::exception &e )
  {
    throw std::runtime_error("Failed to parse config XML: " + string(e.what()));
  }//try / catch
  
  assert( 0 );
  return nullptr;
}//LlmConfig LlmConfig::loadFromFile( const std::string &filename )


bool LlmConfig::saveToFile( const LlmConfig &config, const std::string &filename )
{
  try
  {
    rapidxml::xml_document<char> doc;
    
    // Create root node
    rapidxml::xml_node<char>* root = doc.allocate_node(rapidxml::node_element, "LlmConfig");
    doc.append_node(root);
    
    XmlUtils::append_version_attrib( root, LlmConfig::sm_xmlSerializationVersion );
    
    // LLM API settings
    rapidxml::xml_node<char>* llmApi = doc.allocate_node(rapidxml::node_element, "LlmApi");
    root->append_node(llmApi);
    XmlUtils::append_bool_node(   llmApi, "Enabled",            config.llmApi.enabled );
    XmlUtils::append_string_node( llmApi, "ApiEndpoint",        config.llmApi.apiEndpoint );
    XmlUtils::append_string_node( llmApi, "BearerToken",        config.llmApi.bearerToken );
    XmlUtils::append_string_node( llmApi, "Model",              config.llmApi.model );
    XmlUtils::append_int_node(    llmApi, "MaxTokens",          config.llmApi.maxTokens );
    XmlUtils::append_int_node(    llmApi, "ContextLengthLimit", config.llmApi.contextLengthLimit );
    XmlUtils::append_string_node( llmApi, "SystemPrompt",       config.llmApi.systemPrompt );
    
    
    // MCP server settings
    rapidxml::xml_node<char> *mcpServer = doc.allocate_node(rapidxml::node_element, "McpServer");
    root->append_node(mcpServer);
    XmlUtils::append_bool_node(   mcpServer, "Enabled",     config.mcpServer.enabled );
#if( MCP_ENABLE_AUTH )
    XmlUtils::append_string_node( mcpServer, "BearerToken", config.mcpServer.bearerToken );
#endif
    
    // Write to file
#ifdef _WIN32
    const std::wstring woutcsv = SpecUtils::convert_from_utf8_to_utf16(filename);
    std::ofstream file( woutcsv.c_str(), ios::binary | ios::out );
#else
    std::ofstream file( filename.c_str(), ios::binary | ios::out);
#endif
    
    if( !file.is_open() )
      return false;
    
    //file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << doc;
    
    return file.good();
  }catch( const std::exception &e )
  {
    cout << "Failed to save config: " << e.what() << endl;
    return false;
  }//try / catch
  
  assert( 0 );
  return false; //cant actually get here
}//bool LlmConfig::saveToFile( const LlmConfig &config, const std::string &filename )


std::string LlmConfig::getUserConfigPath()
{
  return SpecUtils::append_path(InterSpec::writableDataDirectory(), "llm_config.xml");
}


std::string LlmConfig::getDefaultConfigPath()
{
  return SpecUtils::append_path(InterSpec::staticDataDirectory(), "llm_config.xml");
}
