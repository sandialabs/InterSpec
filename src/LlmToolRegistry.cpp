#include "InterSpec_config.h"
#include "InterSpec/LlmToolRegistry.h"

#if( USE_LLM_INTERFACE )

#include <iostream>
#include <stdexcept>

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/AnalystChecks.h"

#include <Wt/WApplication>

#include "SpecUtils/SpecFile.h"

using namespace std;
using json = nlohmann::json;

namespace {
  // JSON conversion for SpecUtils::SpectrumType enum  
  NLOHMANN_JSON_SERIALIZE_ENUM(SpecUtils::SpectrumType, {
      {SpecUtils::SpectrumType::Foreground, "Foreground"},
      {SpecUtils::SpectrumType::Background, "Background"},
      {SpecUtils::SpectrumType::SecondForeground, "Secondary"},
  })

  void from_json(const json& j, AnalystChecks::DetectedPeaksOptions& p) {
    std::string specTypeStr = j.at("specType").get<std::string>();
    if (specTypeStr == "Foreground") {
      p.specType = SpecUtils::SpectrumType::Foreground;
    } else if (specTypeStr == "Background") {
      p.specType = SpecUtils::SpectrumType::Background;
    } else if (specTypeStr == "Secondary") {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    } else {
      throw std::runtime_error("Invalid spectrum type: " + specTypeStr);
    }
    p.userSession = j.value("userSession", std::optional<std::string>{});
  }

  void to_json(json& j, const PeakDef& p) {
    j = json{{"lowerEnergy", p.lowerX()}, {"upperEnergy", p.upperX()} };

    if( p.type() == PeakDef::DefintionType::GaussianDefined )
    {
      j["fwhm"] = p.fwhm();
      j["energy"] = p.mean();
      j["amplitude"] = p.peakArea();
      if( p.amplitudeUncert() > 0.0 )
      {
        j["num_sigma"] = p.amplitudeUncert() / p.peakArea();
      }
    }else if( p.type() == PeakDef::DefintionType::DataDefined )
    {
      j["type"] = "DataDefined";
    }
    
    const uintptr_t ptr_val = reinterpret_cast<uintptr_t>(p.continuum().get());
    j["roiID"] = static_cast<uint64_t>( ptr_val );
  }

  void to_json(json& j, const AnalystChecks::DetectedPeakStatus& p) {
    json peaks = json::array();
    for( const auto& peak : p.peaks ) {
      json peak_json;
      to_json(peak_json, *peak);
      peaks.push_back(peak_json);
    }
    j = json{{"userSession", p.userSession},
      {"peaks", peaks}};
  }
}

namespace LlmTools {

ToolRegistry& ToolRegistry::instance() {
  static ToolRegistry registry;
  return registry;
}

void ToolRegistry::registerTool(const SharedTool& tool) {
  m_tools[tool.name] = tool;
}

void ToolRegistry::registerDefaultTools() {
  if (m_defaultToolsRegistered) {
    return;
  }
  
  cout << "Registering default LLM tools..." << endl;
  
  // Register detected_peaks tool
  registerTool({
    "detected_peaks",
    "Returns all peaks detected by automated peak search. Does not add peaks to the user working session.",
    json::parse(R"({
      "type": "object",
      "properties": {
          "specType": { 
            "type": "string", 
            "description": "Which displayed spectrum to search for peaks in; the user is almost always interested in the Foreground, except to check if a peak is in both the foreground and background.", 
            "enum": ["Foreground", "Background", "Secondary"] 
          },
          "userSession": { 
            "type": "string", 
            "description": "Optional: the user session identifier.  If not specified, will use most recent session, if still alive, otherwise the most recently started session." 
          }
      },
      "required": ["specType"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executePeakDetection(params, interspec);
    }
  });
  
  // Register spectrum_info tool
  registerTool({
    "get_spectrum_info", 
    "Get basic information about the currently loaded spectrum",
    json::parse(R"({
      "type": "object",
      "properties": {
          "specType": { 
            "type": "string", 
            "description": "Which spectrum to get info for", 
            "enum": ["Foreground", "Background", "Secondary"] 
          }
      },
      "required": ["specType"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return getSpectrumInfo(params, interspec);
    }
  });
  
  // Register test tool
  registerTool({
    "test_tool",
    "A simple test tool that returns session information",
    json::parse(R"({
      "type": "object",
      "properties": {},
      "required": []
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      json result;
      result["message"] = "Test tool executed successfully";
      result["sessionId"] = interspec ? "valid_session" : "null_session";
      result["timestamp"] = chrono::duration_cast<chrono::seconds>(
        chrono::system_clock::now().time_since_epoch()).count();
      return result;
    }
  });
  
  m_defaultToolsRegistered = true;
  cout << "Registered " << m_tools.size() << " default tools" << endl;
}

const std::map<std::string, SharedTool>& ToolRegistry::getTools() const {
  return m_tools;
}

const SharedTool* ToolRegistry::getTool(const std::string& name) const {
  auto it = m_tools.find(name);
  return (it != m_tools.end()) ? &it->second : nullptr;
}

nlohmann::json ToolRegistry::executeTool(const std::string& toolName, 
                                       const nlohmann::json& parameters, 
                                       InterSpec* interspec) {
  const SharedTool* tool = getTool(toolName);
  if (!tool) {
    throw std::runtime_error("Tool not found: " + toolName);
  }
  
  try {
    cout << "Executing tool: " << toolName << " with params: " << parameters.dump() << endl;
    json result = tool->executor(parameters, interspec);
    cout << "Tool result: " << result.dump() << endl;
    return result;
  } catch (const std::exception& e) {
    throw std::runtime_error("Tool execution failed for " + toolName + ": " + e.what());
  }
}

void ToolRegistry::clearTools() {
  m_tools.clear();
  m_defaultToolsRegistered = false;
}

// Implementation of specific tool functions
json ToolRegistry::executePeakDetection(const json& params, InterSpec* interspec) {
  if (!interspec) {
    throw std::runtime_error("No InterSpec session available");
  }
  
  // Parse parameters into DetectedPeaksOptions
  AnalystChecks::DetectedPeaksOptions options;
  from_json(params, options);
  
  // Call the AnalystChecks function to perform the actual peak detection
  AnalystChecks::DetectedPeakStatus result = AnalystChecks::detected_peaks(options, interspec);
  
  // Convert the result to JSON and return
  json result_json;
  to_json(result_json, result);

  return result_json;
}

json ToolRegistry::getSpectrumInfo(const json& params, InterSpec* interspec) {
  if (!interspec) {
    throw std::runtime_error("No InterSpec session available");
  }
  
  string specTypeStr = params.at("specType").get<string>();
  SpecUtils::SpectrumType specType;
  
  if (specTypeStr == "Foreground") specType = SpecUtils::SpectrumType::Foreground;
  else if (specTypeStr == "Background") specType = SpecUtils::SpectrumType::Background;
  else if (specTypeStr == "Secondary") specType = SpecUtils::SpectrumType::SecondForeground;
  else throw std::runtime_error("Invalid spectrum type: " + specTypeStr);
  
  std::shared_ptr<SpecMeas> meas = interspec->measurment(specType);
  if (!meas) {
    throw std::runtime_error("No measurement loaded for " + specTypeStr + " spectrum");
  }
  
  shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(specType);
  if (!spectrum) {
    throw std::runtime_error("No spectrum displayed for " + specTypeStr + " spectrum");
  }
  
  json result;
  result["specType"] = specTypeStr;
  result["detectorName"] = spectrum->detector_name();
  result["liveTime"] = spectrum->live_time();
  result["realTime"] = spectrum->real_time();
  result["numChannels"] = spectrum->num_gamma_channels();
  result["energyCalibration"] = spectrum->calibration_coeffs();
  
  if (spectrum->gamma_counts() && !spectrum->gamma_counts()->empty()) {
    result["totalCounts"] = spectrum->gamma_count_sum();
  }
  
  return result;
}

} // namespace LlmTools

#endif // USE_LLM_INTERFACE
