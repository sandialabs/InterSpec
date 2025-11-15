#include "InterSpec_config.h"
#include "InterSpec/LlmConfig.h"

#include <iostream>
#include <fstream>
#include <stdexcept>

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

#include <nlohmann/json.hpp>

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/XmlUtils.hpp"



using namespace std;

static_assert( USE_LLM_INTERFACE, "You should not include this library unless USE_LLM_INTERFACE is enabled" );

std::string agentTypeToString( AgentType type )
{
  switch( type )
  {
    case AgentType::MainAgent:    return "MainAgent";
    case AgentType::NuclideId:    return "NuclideId";
    case AgentType::ActivityFit:  return "ActivityFit";
  }

  throw std::invalid_argument( "Unknown AgentType" );
}//agentTypeToString(...)


AgentType stringToAgentType( const std::string &name )
{
  if( name == "MainAgent" )    return AgentType::MainAgent;
  if( name == "NuclideId" )    return AgentType::NuclideId;
  if( name == "ActivityFit" )  return AgentType::ActivityFit;

  throw std::invalid_argument( "Unknown agent name: " + name );
}//stringToAgentType(...)

const std::string LlmConfig::McpServer::sm_invalid_bearer_token = "INVALID-BEARER-TOKEN";


std::shared_ptr<LlmConfig> LlmConfig::load()
{
  // Determine which path to use for each of the three config files
  // User files override default files
  
  auto get_config_file_path = []( const string &filename ){
    const string user_data_dir = ([](){ try{ return InterSpec::writableDataDirectory(); }catch(...){ return ""s; } })();
    if( !user_data_dir.empty() && SpecUtils::is_file(SpecUtils::append_path(user_data_dir,filename)) )
      return SpecUtils::append_path(user_data_dir,filename);
    const string static_data_dir = InterSpec::staticDataDirectory();
    return SpecUtils::append_path(static_data_dir,filename);
  };
  
  const string configPath = get_config_file_path( "llm_config.xml" );
  const string agentsPath = get_config_file_path( "llm_agents.xml" );
  const string toolsPath = get_config_file_path( "llm_tools_config.xml" );
  
  // If we have config file, try to load it
  std::shared_ptr<LlmConfig> config = make_shared<LlmConfig>();
  
  try
  {
    std::pair<LlmApi, McpServer> apiAndMcp = loadApiAndMcpConfigs(configPath);
    config->llmApi = std::move(apiAndMcp.first);
    config->mcpServer = std::move(apiAndMcp.second);
  }catch( const std::exception &e )
  {
    // If user provided any override file and it failed to parse, return nullptr
    cerr << "Failed to load LLM config from '" << configPath << "': " << e.what() << endl;
    throw; // rethrow
  }
  
  if( !config->llmApi.enabled && !config->mcpServer.enabled )
  {
    cout << "Not enabling LLM assistant or MCP server." << endl;
    return config;
  }
  
  try
  {
    config->agents = loadAgentsFromFile(agentsPath);

    // Verify that all required agents have been loaded
    auto requireAgent = []( const AgentType type, const LlmConfig &config ) {
      for( const AgentConfig &agent : config.agents )
      {
        if( agent.type == type )
          return;
      }
      throw std::runtime_error( "Required agent '" + agentTypeToString(type) + "' not found in agent configuration" );
    };

    requireAgent( AgentType::MainAgent, *config );
    requireAgent( AgentType::NuclideId, *config );
    requireAgent( AgentType::ActivityFit, *config );
  }catch( const std::exception &e )
  {
    // If user provided any override file and it failed to parse, return nullptr
    cerr << "Failed to load LLM agent config from '" << agentsPath << "': " << e.what() << endl;
    throw; // rethrow
  }
   
  try
  {
    config->tools = loadToolConfigsFromFile(toolsPath);
  }catch( const std::exception &e )
  {
    // If user provided any override file and it failed to parse, return nullptr
    cerr << "Failed to load LLM tool config from '" << toolsPath << "': " << e.what() << endl;
    throw; // rethrow
  }
  
#if( BUILD_AS_UNIT_TEST_SUITE )
  // We dont have any unit-tests that call out to a LLM, but we do have some tests that require a valid LLM config,
  //  so we will stub a dummy one in here.
  config->llmApi.enabled = true;
  config->llmApi.apiEndpoint = "http://localhost:1234556787";
  config->llmApi.bearerToken = "AnInvalidBearerToken...";
  config->llmApi.model = "AnInvalidMdoel";
  config->llmApi.maxTokens = 4000;
  config->llmApi.contextLengthLimit = 128000;
#endif
  
  return config;
}//LlmConfig LlmConfig::load()


std::pair<LlmConfig::LlmApi, LlmConfig::McpServer> LlmConfig::loadApiAndMcpConfigs( const std::string &llmConfigPath )
{
  try
  {
    std::vector<char> xmlContent;
    SpecUtils::load_file_data( llmConfigPath.c_str(), xmlContent );

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

    LlmApi llmApi;
    McpServer mcpServer;

    {// Begin load LLM API settings
      const rapidxml::xml_node<char> * const llmApiNode = XML_FIRST_NODE(root, "LlmApi");
      if( !llmApiNode )
        throw runtime_error( "Missing 'LlmApi' node." );

      static_assert( LlmApi::sm_xmlSerializationVersion == 0,
                    "needs to be updated for new serialization version." );
      XmlUtils::check_xml_version( llmApiNode, LlmApi::sm_xmlSerializationVersion );

      llmApi.enabled = XmlUtils::get_bool_node_value(llmApiNode, "Enabled");
      llmApi.apiEndpoint = XmlUtils::get_string_node_value(llmApiNode, "ApiEndpoint");
      llmApi.bearerToken = XmlUtils::get_string_node_value(llmApiNode, "BearerToken");
      llmApi.model = XmlUtils::get_string_node_value(llmApiNode, "Model");
      llmApi.maxTokens = XmlUtils::get_int_node_value(llmApiNode, "MaxTokens");
      llmApi.contextLengthLimit = XmlUtils::get_int_node_value(llmApiNode, "ContextLengthLimit");
    }// End load LLM API settings

    {// Begin load MCP server settings
      const rapidxml::xml_node<char> * const mcpServerNode = XML_FIRST_NODE(root, "McpServer");
      if( !mcpServerNode )
        throw runtime_error( "Missing 'McpServer' node." );

      static_assert( McpServer::sm_xmlSerializationVersion == 0,
                    "needs to be updated for new serialization version." );
      XmlUtils::check_xml_version( mcpServerNode, McpServer::sm_xmlSerializationVersion );

      mcpServer.enabled = XmlUtils::get_bool_node_value(mcpServerNode, "Enabled");
#if( MCP_ENABLE_AUTH )
      mcpServer.bearerToken = XmlUtils::get_string_node_value( mcpServerNode, "BearerToken" );
#endif
    }// End load MCP server settings

    return std::make_pair( llmApi, mcpServer );
  }catch( const std::exception &e )
  {
    throw std::runtime_error("Failed to parse config XML: " + string(e.what()));
  }//try / catch
}//LlmConfig::loadApiAndMcpConfigs( const std::string &llmConfigPath )


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
    XmlUtils::append_version_attrib( llmApi, LlmApi::sm_xmlSerializationVersion );
    XmlUtils::append_bool_node(   llmApi, "Enabled",            config.llmApi.enabled );
    XmlUtils::append_string_node( llmApi, "ApiEndpoint",        config.llmApi.apiEndpoint );
    XmlUtils::append_string_node( llmApi, "BearerToken",        config.llmApi.bearerToken );
    XmlUtils::append_string_node( llmApi, "Model",              config.llmApi.model );
    XmlUtils::append_int_node(    llmApi, "MaxTokens",          config.llmApi.maxTokens );
    XmlUtils::append_int_node(    llmApi, "ContextLengthLimit", config.llmApi.contextLengthLimit );

    // NOTE: Agents and tools are now saved separately in llm_agents.xml and llm_tools_config.xml,
    // not in llm_config.xml


    // MCP server settings
    rapidxml::xml_node<char> *mcpServer = doc.allocate_node(rapidxml::node_element, "McpServer");
    root->append_node(mcpServer);
    XmlUtils::append_version_attrib( mcpServer, McpServer::sm_xmlSerializationVersion );
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


std::vector<LlmConfig::ToolConfig> LlmConfig::loadToolConfigsFromFile( const std::string &toolsConfigPath )
{
  if( !SpecUtils::is_file(toolsConfigPath) )
  {
    throw std::runtime_error("Tool configuration file not found: " + toolsConfigPath);
  }

  try
  {
    std::vector<char> xmlContent;
    SpecUtils::load_file_data( toolsConfigPath.c_str(), xmlContent );

    rapidxml::xml_document<char> doc;
    const int flags = (rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace);

    doc.parse<flags>( &xmlContent[0] );

    // First check for LlmConfig wrapper (for future compatibility)
    const rapidxml::xml_node<char> * llmConfigNode = XML_FIRST_NODE(&doc, "LlmConfig");
    if( llmConfigNode )
    {
      static_assert( LlmConfig::sm_xmlSerializationVersion == 0,
                    "needs to be updated for new serialization version." );
      XmlUtils::check_xml_version( llmConfigNode, LlmConfig::sm_xmlSerializationVersion );
    }

    // Find ToolDefinitions node (either at root or inside LlmConfig)
    const rapidxml::xml_node<char> * const root = llmConfigNode
      ? XML_FIRST_NODE(llmConfigNode, "ToolDefinitions")
      : XML_FIRST_NODE(&doc, "ToolDefinitions");

    if( !root )
      throw std::runtime_error("Missing ToolDefinitions node in " + toolsConfigPath);

    static_assert( ToolConfig::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    XmlUtils::check_xml_version( root, ToolConfig::sm_xmlSerializationVersion );

    std::vector<LlmConfig::ToolConfig> tools;

    // Parse each Tool node
    XML_FOREACH_CHILD(toolNode, root, "Tool")
    {
      LlmConfig::ToolConfig tool;

      const rapidxml::xml_attribute<char> * const nameAttr = toolNode->first_attribute("name");
      if( !nameAttr || !nameAttr->value() || (nameAttr->value_size() == 0) )
        throw runtime_error( "Tool node missing 'name' attribute in " + toolsConfigPath );

      tool.name = SpecUtils::xml_value_str( nameAttr );

      // Load descriptions
      XML_FOREACH_CHILD(descNode, toolNode, "Description")
      {
        const rapidxml::xml_attribute<char> * const roleAttr = descNode->first_attribute("role");
        const string roleName = SpecUtils::xml_value_str( roleAttr );

        const string descText = SpecUtils::xml_value_str( descNode );

        if( roleName.empty() )
        {
          // No role attribute - this is the default description
          tool.defaultDescription = descText;
        }else
        {
          try
          {
            const AgentType agentType = stringToAgentType( roleName );
            tool.roleDescriptions[agentType] = descText;
          }catch( const std::exception &e )
          {
            cerr << "Warning: Unknown agent type '" << roleName << "' in tool description for '" << tool.name << "', skipping" << endl;
          }
        }
      }//for( loop over Description nodes )

      // Load AvailableFor
      const rapidxml::xml_node<char> * const availableForNode = XML_FIRST_NODE(toolNode, "AvailableFor");
      if( availableForNode )
      {
        for( const rapidxml::xml_node<char> *agentNode = availableForNode->first_node("Agent");
            agentNode; agentNode = agentNode->next_sibling("Agent") )
        {
          if( agentNode->value() && (agentNode->value_size() > 0) )
          {
            const string agentName = agentNode->value();
            try
            {
              const AgentType agentType = stringToAgentType( agentName );
              tool.availableForAgents.push_back( agentType );
            }catch( const std::exception &e )
            {
              cerr << "Warning: Unknown agent type '" << agentName << "' in tool availableFor for '" << tool.name << "', skipping" << endl;
            }
          }
        }//for( loop over Agent nodes in AvailableFor )
      }//if( availableForNode )

      // Load ParametersSchema
      const rapidxml::xml_node<char> * const schemaNode = XML_FIRST_NODE(toolNode, "ParametersSchema");
      if( schemaNode && schemaNode->value() && (schemaNode->value_size() > 0) )
      {
        const string schemaStr = SpecUtils::xml_value_str( schemaNode );
        try
        {
          tool.parametersSchema = nlohmann::json::parse( schemaStr );
        }catch( const std::exception &e )
        {
          throw std::runtime_error( "Failed to parse JSON schema for tool '" + tool.name + "': " + string(e.what()) );
        }
      }else
      {
        // No schema provided - set to empty object
        tool.parametersSchema = nlohmann::json::object();
      }

      tools.push_back(tool);
    }//for( loop over Tool nodes )

    cout << "Loaded " << tools.size() << " tool configurations from " + toolsConfigPath << endl;
    return tools;
  }catch( const std::exception &e )
  {
    throw std::runtime_error("Failed to parse tool config XML from " + toolsConfigPath + ": " + string(e.what()));
  }//try / catch
}//loadToolConfigsFromFile(...)


std::vector<LlmConfig::AgentConfig> LlmConfig::loadAgentsFromFile( const std::string &agentsConfigPath )
{
  if( !SpecUtils::is_file(agentsConfigPath) )
    throw std::runtime_error("Agent configuration file not found: " + agentsConfigPath);

  try
  {
    std::vector<char> xmlContent;
    SpecUtils::load_file_data( agentsConfigPath.c_str(), xmlContent );

    rapidxml::xml_document<char> doc;
    const int flags = (rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace);

    doc.parse<flags>( &xmlContent[0] );

    // First check for LlmConfig wrapper (for future compatibility)
    const rapidxml::xml_node<char> * llmConfigNode = XML_FIRST_NODE(&doc, "LlmConfig");
    if( llmConfigNode )
    {
      static_assert( LlmConfig::sm_xmlSerializationVersion == 0,
                    "needs to be updated for new serialization version." );
      XmlUtils::check_xml_version( llmConfigNode, LlmConfig::sm_xmlSerializationVersion );
    }

    // Find AgentDefinitions node (either at root or inside LlmConfig)
    const rapidxml::xml_node<char> * const root = llmConfigNode
      ? XML_FIRST_NODE(llmConfigNode, "AgentDefinitions")
      : XML_FIRST_NODE(&doc, "AgentDefinitions");

    if( !root )
      throw std::runtime_error("Missing AgentDefinitions node in " + agentsConfigPath);

    static_assert( AgentConfig::sm_xmlSerializationVersion == 0,
                  "needs to be updated for new serialization version." );
    XmlUtils::check_xml_version( root, AgentConfig::sm_xmlSerializationVersion );

    std::vector<LlmConfig::AgentConfig> agents;

    // Parse each Agent node
    XML_FOREACH_CHILD( agentNode, root, "Agent" )
    {
      LlmConfig::AgentConfig agent;

      const rapidxml::xml_attribute<char> * const nameAttr = agentNode->first_attribute("name");
      if( !nameAttr || !nameAttr->value() || (nameAttr->value_size() == 0) )
        throw runtime_error( "Agent node missing 'name' attribute in " + agentsConfigPath );

      agent.name = SpecUtils::xml_value_str( nameAttr );

      try
      {
        agent.type = stringToAgentType( agent.name );
      }catch( const std::exception &e )
      {
        throw std::runtime_error( "Unknown agent type '" + agent.name + "' in agent config" );
      }

      // Load Description
      const rapidxml::xml_node<char> * const descNode = XmlUtils::get_required_node( agentNode, "Description" );
      agent.description = SpecUtils::xml_value_str( descNode );

      // Load SystemPrompt
      const rapidxml::xml_node<char> * const promptNode = XmlUtils::get_required_node( agentNode, "SystemPrompt" );
      agent.systemPrompt = SpecUtils::xml_value_str( promptNode );

      agents.push_back(agent);
    }//for( loop over Agent nodes )

    cout << "Loaded " << agents.size() << " agent configurations from " + agentsConfigPath << endl;
    return agents;
  }catch( const std::exception &e )
  {
    throw std::runtime_error("Failed to parse agent config XML from " + agentsConfigPath + ": " + string(e.what()));
  }//try / catch
}//loadAgentsFromFile(...)
