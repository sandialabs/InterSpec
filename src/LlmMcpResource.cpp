#include "InterSpec_config.h"
#include "InterSpec/LlmMcpResource.h"


#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <optional>
#include <functional>

#include <Wt/WServer>
#include <Wt/WResource>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/Http/Request>
#include <Wt/Http/Response>

#include "external_libs/SpecUtils/3rdparty/nlohmann/json.hpp"

#include "SpecUtils/SpecFile.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/LlmToolRegistry.h"


static_assert( USE_LLM_INTERFACE, "This file should not be being compiled unless USE_LLM_INTERFACE is enabled" );
#if( BUILD_FOR_WEB_DEPLOYMENT || defined(IOS) || defined(ANDROID) )
static_assert( 0, "MCP server can not be enabled for web or mobile builds" );
#endif


// for convenience
using namespace std;
using namespace Wt;
using json = nlohmann::json;

/*
nlohmann::json detected_peaks( const DetectedPeaksOptions options )
{
  set<InterSpecApp *> instances = InterSpecApp::runningInstances();
  InterSpecApp *app = nullptr;
  string user_session;
  if( options.userSession.has_value() && !options.userSession->empty() )
  {
    for( InterSpecApp *app_instance : instances )
    {
      Wt::WApplication::UpdateLock lock(app_instance);
      if( lock &&
          ((app_instance->sessionId() == options.userSession.value())
           || (app_instance->externalToken() == options.userSession.value())) )
      {
        app = app_instance;
        user_session = options.userSession.value();
        break;
      }
    }   
  }else
  {
    std::chrono::steady_clock::time_point earliest_start_time{};
    for( InterSpecApp *app_instance : instances )
    {
      Wt::WApplication::UpdateLock lock(app_instance);
      if( lock &&
          (!app || (earliest_start_time > app_instance->startTime())) )
      {
        app = app_instance;
        user_session = app->sessionId();
        earliest_start_time = app_instance->startTime();
      }
    }
  }//if( options.userSession.has_value() && !options.userSession->empty() ) / else
  
  if( !app )
    throw std::runtime_error("No InterSpecApp instance found for user session: " + options.userSession.value());

  WApplication::UpdateLock lock(app);
  if( !lock )
    throw std::runtime_error("Failed to get app lock for user session: " + options.userSession.value());
  
  
  
  InterSpec *viewer = app->viewer();
  if( !viewer )
    throw std::runtime_error("Failed to get viewer for user session: " + options.userSession.value());

  std::shared_ptr<SpecMeas> meas = viewer->measurment( options.specType );
  if( !meas )
    throw std::runtime_error("No measurement loaded for " 
        + std::string(SpecUtils::descriptionText(options.specType)) 
        + " spectrum");

  const set<int> sample_nums = viewer->displayedSamples( options.specType );
  if( sample_nums.empty() )
    throw std::runtime_error("No samples displayed for " 
        + std::string(SpecUtils::descriptionText(options.specType)) 
        + " spectrum");

  shared_ptr<const SpecUtils::Measurement> spectrum = viewer->displayedHistogram( options.specType );
  if( !spectrum )
    throw std::runtime_error("No spectrum displayed for " 
        + std::string(SpecUtils::descriptionText(options.specType)) 
        + " spectrum");
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> auto_peaks = meas->automatedSearchPeaks( sample_nums );
  shared_ptr<const deque<shared_ptr<const PeakDef>>> user_peaks = meas->peaks( sample_nums );
  
  if( !auto_peaks )
  {
    // Search for peaks
    const bool singleThreaded = false;
    const bool isHPGe = PeakFitUtils::is_likely_high_res( viewer );
    const auto det = meas->detector();
    const vector<shared_ptr<const PeakDef>> found_auto_peaks
              = ExperimentalAutomatedPeakSearch::search_for_peaks( spectrum, det, user_peaks, singleThreaded, isHPGe );
      
    auto autopeaksdeque = make_shared<std::deque<std::shared_ptr<const PeakDef>>>( begin(found_auto_peaks), end(found_auto_peaks) );
    auto_peaks = autopeaksdeque;
    meas->setAutomatedSearchPeaks( sample_nums, autopeaksdeque );
  }//if( !auto_peaks )
  
  vector<shared_ptr<const PeakDef>> all_peaks;

  if( user_peaks )
    all_peaks.insert( all_peaks.end(), user_peaks->begin(), user_peaks->end() );

  // We will add auto-search peaks only if there isnt already a user peak at essentually the same energy
  if( auto_peaks )
  {
    //Lambda to find the nearest peak, so far, to a given energy.
    auto nearest_peak = [&all_peaks]( const float energy) -> std::shared_ptr<const PeakDef> {
      std::shared_ptr<const PeakDef> nearest;
      double minDE = std::numeric_limits<double>::infinity();

      for( const auto &peak : all_peaks )
      {
        const double dE = fabs( peak->mean() - energy );
        if( (dE < minDE)
            && ((energy > peak->lowerX()) && (energy < peak->upperX())) )
        {
        minDE = dE;
        nearest = peak;
        }//if( dE < minDE )
      }//for( const auto &peak : all_peaks )

      return nearest;
    };//nearest_peak(...)

    for( const shared_ptr<const PeakDef> &peak : *auto_peaks )
    {
      auto nearpeak = nearest_peak( peak->mean() );
      const double peak_sigma = peak->gausPeak() ? peak->sigma() : 0.25*peak->roiWidth();
      if( !nearpeak || (fabs(nearpeak->mean() - peak->mean()) > peak_sigma) )
        all_peaks.push_back( peak );
    }
  }//if( auto_peaks )
  
  std::sort( begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr );
  
  DetectedPeakStatus result;
  result.userSession = user_session;
  result.peaks = std::move(all_peaks);
  
  json result_json;
  to_json(result_json, result);

  return result_json;
}//nlohmann::json detected_peaks( SpecUtils::SpectrumType specType )
*/


std::vector<std::string> nuclide_id_candidates() {
    std::cout << "-> Called nuclide_id_candidates" << std::endl;
    return {"U-235", "U-238", "Pu-239"};
}

struct EnrichmentStatus {
    bool success;
    std::string nuclide;
    double mass_frac;
    double mass_frac_uncert;
};

void to_json(json& j, const EnrichmentStatus& p) {
    j = json{{"success", p.success},
             {"nuclide", p.nuclide},
             {"mass_frac", p.mass_frac},
             {"mass_frac_uncert", p.mass_frac_uncert}};
}

struct EnrichOptions {
    bool background_subtract;
    std::optional<std::string> parameters_set;
    double peak_threshold;
};

void from_json(const json& j, EnrichOptions& p) {
    j.at("background_subtract").get_to(p.background_subtract);
    p.parameters_set = j.value("parameters_set", std::optional<std::string>{});
    j.at("peak_threshold").get_to(p.peak_threshold);
}

EnrichmentStatus get_enrichment(const EnrichOptions& options) {
    std::cout << "-> Called get_enrichment with background_subtract: " << options.background_subtract
              << ", threshold: " << options.peak_threshold << std::endl;
    if (options.parameters_set) {
        std::cout << "   with parameter_set: " << *options.parameters_set << std::endl;
    }
    return {true, "U-235", 0.045, 0.001};
}

// --- LlmMcpResource Implementation ---


LlmMcpResource::LlmMcpResource( const std::shared_ptr<const LlmConfig> &config )
: Wt::WResource(),
  llm_config_( config ),
  tool_registry_( config ? make_unique<LlmTools::ToolRegistry>(*config) : unique_ptr<LlmTools::ToolRegistry>() )  //This can throw
{
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource must be initialized with a valid config" );
  
  if( !config->mcpServer.enabled )
    throw std::logic_error( "Cannot initialize LlmMcpResource if MCP isnt enabled in the config" );
  
  register_default_tools();
}


LlmMcpResource::~LlmMcpResource() = default;

void LlmMcpResource::register_tool(Tool tool) {
    registered_tools_[tool.name] = std::move(tool);
}

void LlmMcpResource::handleRequest(const Wt::Http::Request& request, Wt::Http::Response& response)
{
  //cout << "McpResource: handleRequest, path='" << request.pathInfo() << "', method='" << request.method() << "'" << endl;
  //for( auto& header : request.headers() ) {
  //  cout << "  Header: " << header.name() << " = " << header.value() << endl;
  //}
  
  // Set headers required for some MCP clients
  response.addHeader("Access-Control-Allow-Origin", "*");
  response.addHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
  response.addHeader("Access-Control-Allow-Headers", "Content-Type, Authorization");
  
  // Handle CORS preflight requests
  if( request.method() == "OPTIONS" )
  {
    response.setStatus(204); // No Content
    return;
  }//if( request.method() == "OPTIONS" )
  
  try
  {
    const std::string& path = request.pathInfo();
    
    // Handle root path (empty string) - this is what LM Studio sends
    if( (path == "") || (path == "/") )
    {
      if( request.method() == "GET" )
      {
        // Check if client wants SSE (Server-Sent Events)
        std::string accept = request.headerValue("Accept");
        if( accept.find("text/event-stream") != std::string::npos )
        {
          // Handle SSE connection
          cout << "Handling SSE connection" << endl;
          response.setMimeType("text/event-stream");
          response.addHeader("Cache-Control", "no-cache");
          response.addHeader("Connection", "keep-alive");
          response.addHeader("Access-Control-Allow-Origin", "*");
          
          // Send initial SSE event
          response.out() << "data: {\"jsonrpc\":\"2.0\",\"method\":\"notifications/initialized\",\"params\":{}}\n\n";
          response.out().flush();
          return;
        } else {
          // Regular GET request - return auth info (MCP discovery)
          handle_get_auth(request, response);
          return;
        }
      } else if (request.method() == "POST") {
        // POST request - handle MCP operations
#if( MCP_ENABLE_AUTH )
        if( !validate_request( request,response ) )
          return;
#endif
        handle_call_tool(request, response);
        return;
      }
    }
    
    // Handle specific endpoint paths (for backward compatibility)
    if (path == "/get_auth" || path == "/mcp/get_auth") {
      handle_get_auth(request, response);
      return;
    }
    
#if( MCP_ENABLE_AUTH )
    // For all other endpoints, validate the request if auth is enabled.
    if( !validate_request( request,response ) )
      return;
#endif
    
    if (path == "/list_tools" || path == "/mcp/list_tools") {
      handle_list_tools(request, response);
    } else if (path == "/call_tool" || path == "/mcp/call_tool") {
      handle_call_tool(request, response);
    } else {
      response.setStatus(404);
      response.out() << json{{"error", "Not Found"}, {"path", path}}.dump();
    }
  } catch (const std::exception& e) {
    response.setStatus(500);
    response.setMimeType("application/json");
    response.out() << json{{"error", e.what()}}.dump();
  }
}//void handleRequest(...)


#if( MCP_ENABLE_AUTH )
bool LlmMcpResource::validate_request(const Wt::Http::Request& request, Wt::Http::Response &response)
{
  if( secret_token_.empty() )
    return true;
  
    std::string auth_header = request.headerValue("Authorization");
    if( auth_header.empty() )
    {
      response.setStatus(401);
      response.out() << "Missing Authorization header";
      return false;
    }

    std::string auth_str = auth_header;
    std::string prefix = "Bearer ";
    if (auth_str.rfind(prefix, 0) != 0) {
      response.setStatus(401);
      response.out() << "Invalid Authorization header format";
      return false;
    }

    std::string token = auth_str.substr(prefix.length());
    if (token != secret_token_) {
      response.setStatus(401);
      response.out() << "Invalid token";
      return false;
    }
  
  return true;
}
#endif // MCP_ENABLE_AUTH

void LlmMcpResource::handle_get_auth(const Wt::Http::Request& request, Wt::Http::Response& response)
{
  response.setMimeType("application/json");
#ifdef MCP_ENABLE_AUTH
  // If auth is compiled in, tell the client we require a bearer token, if `secret_token_` not empty
  if( secret_token_.empty() )
    response.out() << json{{"type", "none"}}.dump();
  else
    response.out() << json{{"type", "bearer"}}.dump();
#else
  // Otherwise, tell the client no authentication is needed.
  response.out() << json{{"type", "none"}}.dump();
#endif // MCP_ENABLE_AUTH
}

void LlmMcpResource::handle_list_tools(const Wt::Http::Request& request, Wt::Http::Response& response)
{
    json tool_definitions = json::array();
    for (const auto& pair : registered_tools_) {
        const auto& tool = pair.second;
        tool_definitions.push_back({
            {"name", tool.name},
            {"description", tool.description},
            {"parameters", tool.parameters_schema}
        });
    }
    response.setMimeType("application/json");
    response.out() << tool_definitions.dump();
}

void LlmMcpResource::handle_call_tool(const Wt::Http::Request& request, Wt::Http::Response& response)
{
    if (request.method() != "POST") {
         response.setStatus(405); // Method Not Allowed
         response.setMimeType("application/json");
         response.out() << json{{"error", "Method must be POST for tool operations"}}.dump();
         return;
    }

    try {
        std::string request_body;
        request.in().seekg(0, std::ios::end);
        size_t length = request.in().tellg();
        request.in().seekg(0, std::ios::beg);
        request_body.resize(length);
        request.in().read(&request_body[0], length);
        
        cout << "Request body: " << request_body << endl;
        
        json mcp_request = json::parse(request_body);
        cout << "Parsed JSON: " << mcp_request.dump(2) << endl;
        
        // Handle JSON-RPC MCP requests
        if (mcp_request.contains("method")) {
            std::string method = mcp_request["method"];
            json id = mcp_request.value("id", json(nullptr));
            
            if (method == "initialize") {
                // Handle MCP initialization
                cout << "Handling initialize method" << endl;
                response.setMimeType("application/json");
                response.out() << json{
                    {"jsonrpc", "2.0"},
                    {"id", id},
                    {"result", {
                        {"protocolVersion", "2024-11-05"},
                        {"capabilities", {
                            {"tools", json::object()}
                        }},
                        {"serverInfo", {
                            {"name", "InterSpec MCP Server"},
                            {"version", "1.0.0"}
                        }}
                    }}
                                 }.dump();
                 return;
                 
                         } else if (method == "notifications/initialized") {
                // Handle initialization notification - no response needed
                cout << "Handling notifications/initialized - sending 204" << endl;
                response.setStatus(204); // No Content
                return;
                 
             } else if (method == "tools/list") {
                // Return list of available tools
                json tool_definitions = json::array();
                for (const auto& pair : registered_tools_) {
                    const auto& tool = pair.second;
                    tool_definitions.push_back({
                        {"name", tool.name},
                        {"description", tool.description},
                        {"inputSchema", tool.parameters_schema}
                    });
                }
                
                response.setMimeType("application/json");
                response.out() << json{
                    {"jsonrpc", "2.0"},
                    {"id", id},
                    {"result", {{"tools", tool_definitions}}}
                }.dump();
                return;
                
            } else if (method == "tools/call") {
                // Handle tool call
                if (!mcp_request.contains("params")) {
                    throw std::runtime_error("Missing params for tools/call");
                }
                
                json params = mcp_request["params"];
                std::string tool_name = params.at("name");
                json arguments = params.value("arguments", json::object());
                
                auto it = registered_tools_.find(tool_name);
                if (it == registered_tools_.end()) {
                    response.setStatus(404);
                    response.setMimeType("application/json");
                    response.out() << json{
                        {"jsonrpc", "2.0"},
                        {"id", id},
                        {"error", {
                            {"code", -32601},
                            {"message", "Unknown tool: " + tool_name}
                        }}
                    }.dump();
                    return;
                }

                // Call the registered executor function for the tool.
                json result_content = it->second.executor(arguments);

                response.setMimeType("application/json");
                response.out() << json{
                    {"jsonrpc", "2.0"},
                    {"id", id},
                    {"result", {
                        {"content", {{
                            {"type", "text"},
                            {"text", result_content.dump()}
                        }}}
                    }}
                }.dump();
                return;
            }
        }
        
        // Fallback: Handle legacy format
        if (mcp_request.contains("tool_name")) {
            std::string tool_name = mcp_request.at("tool_name");
            json parameters = mcp_request.value("parameters", json::object());

            auto it = registered_tools_.find(tool_name);
            if (it == registered_tools_.end()) {
                throw std::runtime_error("Unknown tool name: " + tool_name);
            }

            // Call the registered executor function for the tool.
            json result_content = it->second.executor(parameters);

            response.setMimeType("application/jsonl");
            response.out() << json{
                {"type", "tool_response"},
                {"tool_name", tool_name},
                {"content", result_content}
            }.dump() << std::endl;
            response.out() << json{
                {"type", "tool_response_final"}
            }.dump() << std::endl;
            return;
        }
        
        // Unknown request format
        response.setStatus(400);
        response.setMimeType("application/json");
        response.out() << json{{"error", "Unrecognized request format"}}.dump();
        
    } catch (const std::exception& e) {
        response.setStatus(500);
        response.setMimeType("application/json");
        response.out() << json{{"error", e.what()}}.dump();
    }
}


void LlmMcpResource::register_default_tools()
{
  assert( tool_registry_ && llm_config_ );
  
  if( !tool_registry_ )
    throw std::logic_error( "LlmMcpResource: should not have null tool registry" );
  
  if( !llm_config_ )
    throw std::logic_error( "LlmMcpResource: should not have null llm config" );
  
  // Create MCP-compatible tools by adapting the shared tools
  for( const auto &toolPair : tool_registry_->getTools() )
  {
    const std::string& name = toolPair.first;
    const LlmTools::SharedTool& sharedTool = toolPair.second;
    
    Tool mcpTool;
    mcpTool.name = sharedTool.name;
    mcpTool.description = sharedTool.description;
    mcpTool.parameters_schema = sharedTool.parameters_schema;
    
    // Add userSession parameter to schema if not already present
    if( mcpTool.parameters_schema.is_object() )
    {
      // Ensure "properties" object exists
      if( !mcpTool.parameters_schema.contains("properties") )
        mcpTool.parameters_schema["properties"] = nlohmann::json::object();

      // Add userSession to properties if not already there
      if( !mcpTool.parameters_schema["properties"].contains("userSession") )
      {
        mcpTool.parameters_schema["properties"]["userSession"] = {
          {"type", "string"},
          {"description", "Optional: the user session identifier.  If not specified, will use most recent session."}
        };
      }
    }
    
    
    // Create an adapter that finds the InterSpec session and calls the shared tool
    mcpTool.executor = [this, name]( json params ) -> json {
      // Extract session ID from params if available
      string sessionId;
      if (params.contains("userSession") && params["userSession"].is_string()) {
        sessionId = params["userSession"];
      }
      
      if( params.contains("userSession") )
        params.erase( "userSession" ); //Tool registry doesnt need to know that this parameter exists
      
      // Find the appropriate InterSpec instance
      InterSpec* interspec = findInterSpecBySessionId(sessionId);
      if (!interspec) {
        throw std::runtime_error("No active InterSpec session found for session: " + 
                                (sessionId.empty() ? "(no session ID provided)" : sessionId));
      }
      
      // Use the shared tool registry to execute the tool
      return tool_registry_->executeTool(name, params, interspec);
    };
    
    register_tool(mcpTool);
  }
  
  cout << "MCP Resource registered " << registered_tools_.size() << " tools from shared registry" << endl;
}//void LlmMcpResource::register_default_tools()

InterSpec* LlmMcpResource::findInterSpecBySessionId(const std::string& sessionId) {
  set<InterSpecApp *> instances = InterSpecApp::runningInstances();
  
  if (sessionId.empty()) {
    // If no session ID provided, find the most recently started instance
    InterSpecApp *mostRecentApp = nullptr;
    std::chrono::steady_clock::time_point earliestStartTime{};
    
    for (InterSpecApp *app_instance : instances) {
      Wt::WApplication::UpdateLock lock(app_instance);
      if (lock && (!mostRecentApp || (earliestStartTime > app_instance->startTime()))) {
        mostRecentApp = app_instance;
        earliestStartTime = app_instance->startTime();
      }
    }
    
    if (mostRecentApp) {
      Wt::WApplication::UpdateLock lock(mostRecentApp);
      if (lock) {
        return mostRecentApp->viewer();
      }
    }
  } else {
    // Look for specific session ID
    for (InterSpecApp *app_instance : instances) {
      Wt::WApplication::UpdateLock lock(app_instance);
      if (lock && 
          ((app_instance->sessionId() == sessionId) || 
           (app_instance->externalToken() == sessionId))) {
        return app_instance->viewer();
      }
    }
  }
  
  return nullptr;
}
