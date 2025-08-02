#include "InterSpec_config.h"
#include "InterSpec/LlmConfig.h"

#if( USE_LLM_INTERFACE )

#include <iostream>
#include <fstream>
#include <stdexcept>

#include "InterSpec/InterSpec.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include <rapidxml/rapidxml.hpp>
#include <rapidxml/rapidxml_print.hpp>

using namespace std;

LlmConfig LlmConfig::load() {
  LlmConfig config;
  
  // Try user config first
  string userPath = getUserConfigPath();
  if (SpecUtils::is_file(userPath)) {
    try {
      cout << "Loading LLM config from: " << userPath << endl;
      return loadFromFile(userPath);
    } catch (const std::exception& e) {
      cout << "Failed to load user config: " << e.what() << endl;
    }
  }
  
  // Try default config
  string defaultPath = getDefaultConfigPath();
  if (SpecUtils::is_file(defaultPath)) {
    try {
      cout << "Loading LLM config from default: " << defaultPath << endl;
      return loadFromFile(defaultPath);
    } catch (const std::exception& e) {
      cout << "Failed to load default config: " << e.what() << endl;
    }
  }
  
  cout << "Using hardcoded LLM config defaults" << endl;
  return config; // Return defaults
}

bool LlmConfig::save(const LlmConfig& config) {
  string userPath = getUserConfigPath();
  return saveToFile(config, userPath);
}

LlmConfig LlmConfig::loadFromFile(const std::string& filename) {
  LlmConfig config;
  
  ifstream file(filename.c_str());
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open config file: " + filename);
  }
  
  // Read entire file
  string xmlContent((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
  
  try {
    rapidxml::xml_document<char> doc;
    doc.parse<0>(&xmlContent[0]);
    
    rapidxml::xml_node<char>* root = doc.first_node("LlmConfig");
    if (!root) {
      throw std::runtime_error("Missing LlmConfig root node");
    }
    
    // Load LLM API settings
    if (rapidxml::xml_node<char>* llmApi = root->first_node("LlmApi")) {
      if (rapidxml::xml_node<char>* node = llmApi->first_node("ApiEndpoint")) {
        config.llmApi.apiEndpoint = node->value();
      }
      if (rapidxml::xml_node<char>* node = llmApi->first_node("BearerToken")) {
        config.llmApi.bearerToken = node->value();
      }
      if (rapidxml::xml_node<char>* node = llmApi->first_node("Model")) {
        config.llmApi.model = node->value();
      }
      if (rapidxml::xml_node<char>* node = llmApi->first_node("MaxTokens")) {
        string val = node->value();
        if (!val.empty()) config.llmApi.maxTokens = std::stoi(val);
      }
      if (rapidxml::xml_node<char>* node = llmApi->first_node("ContextLengthLimit")) {
        string val = node->value();
        if (!val.empty()) config.llmApi.contextLengthLimit = std::stoi(val);
      }
      if (rapidxml::xml_node<char>* node = llmApi->first_node("SystemPrompt")) {
        config.llmApi.systemPrompt = node->value();
      }
    }
    
    // Load MCP server settings
    if (rapidxml::xml_node<char>* mcpServer = root->first_node("McpServer")) {
      if (rapidxml::xml_node<char>* node = mcpServer->first_node("Port")) {
        string val = node->value();
        if (!val.empty()) config.mcpServer.port = std::stoi(val);
      }
      if (rapidxml::xml_node<char>* node = mcpServer->first_node("BearerToken")) {
        config.mcpServer.bearerToken = node->value();
      }
      if (rapidxml::xml_node<char>* node = mcpServer->first_node("Enabled")) {
        string val = node->value();
        config.mcpServer.enabled = (val == "true");
      }
    }
    
    // Load interface settings
    if (rapidxml::xml_node<char>* interface = root->first_node("Interface")) {
      if (rapidxml::xml_node<char>* node = interface->first_node("DefaultVisible")) {
        string val = node->value();
        config.interface.defaultVisible = (val == "true");
      }
      if (rapidxml::xml_node<char>* node = interface->first_node("PanelWidth")) {
        string val = node->value();
        if (!val.empty()) config.interface.panelWidth = std::stoi(val);
      }
    }
    
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to parse config XML: " + string(e.what()));
  }
  
  return config;
}

bool LlmConfig::saveToFile(const LlmConfig& config, const std::string& filename) {
  try {
    rapidxml::xml_document<char> doc;
    
    // Create root node
    rapidxml::xml_node<char>* root = doc.allocate_node(rapidxml::node_element, "LlmConfig");
    doc.append_node(root);
    
    // LLM API settings
    rapidxml::xml_node<char>* llmApi = doc.allocate_node(rapidxml::node_element, "LlmApi");
    root->append_node(llmApi);
    
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "ApiEndpoint", 
                        doc.allocate_string(config.llmApi.apiEndpoint.c_str())));
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "BearerToken", 
                        doc.allocate_string(config.llmApi.bearerToken.c_str())));
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "Model", 
                        doc.allocate_string(config.llmApi.model.c_str())));
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "MaxTokens", 
                        doc.allocate_string(to_string(config.llmApi.maxTokens).c_str())));
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "ContextLengthLimit", 
                        doc.allocate_string(to_string(config.llmApi.contextLengthLimit).c_str())));
    llmApi->append_node(doc.allocate_node(rapidxml::node_element, "SystemPrompt", 
                        doc.allocate_string(config.llmApi.systemPrompt.c_str())));
    
    // MCP server settings
    rapidxml::xml_node<char>* mcpServer = doc.allocate_node(rapidxml::node_element, "McpServer");
    root->append_node(mcpServer);
    
    mcpServer->append_node(doc.allocate_node(rapidxml::node_element, "Port", 
                          doc.allocate_string(to_string(config.mcpServer.port).c_str())));
    mcpServer->append_node(doc.allocate_node(rapidxml::node_element, "BearerToken", 
                          doc.allocate_string(config.mcpServer.bearerToken.c_str())));
    mcpServer->append_node(doc.allocate_node(rapidxml::node_element, "Enabled", 
                          doc.allocate_string(config.mcpServer.enabled ? "true" : "false")));
    
    // Interface settings
    rapidxml::xml_node<char>* interface = doc.allocate_node(rapidxml::node_element, "Interface");
    root->append_node(interface);
    
    interface->append_node(doc.allocate_node(rapidxml::node_element, "DefaultVisible", 
                          doc.allocate_string(config.interface.defaultVisible ? "true" : "false")));
    interface->append_node(doc.allocate_node(rapidxml::node_element, "PanelWidth", 
                          doc.allocate_string(to_string(config.interface.panelWidth).c_str())));
    
    // Write to file
    ofstream file(filename.c_str());
    if (!file.is_open()) {
      return false;
    }
    
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << doc;
    
    return file.good();
    
  } catch (const std::exception& e) {
    cout << "Failed to save config: " << e.what() << endl;
    return false;
  }
}

std::string LlmConfig::getUserConfigPath() {
  return SpecUtils::append_path(InterSpec::writableDataDirectory(), "llm_config.xml");
}

std::string LlmConfig::getDefaultConfigPath() {
  return SpecUtils::append_path(InterSpec::staticDataDirectory(), "llm_config.xml");
}

#endif // USE_LLM_INTERFACE