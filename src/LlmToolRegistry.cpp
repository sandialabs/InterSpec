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
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/PhysicalUnits.h"

#include <Wt/WApplication>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SandiaDecay/SandiaDecay.h"

using namespace std;
using json = nlohmann::json;

namespace {
  // JSON conversion for SpecUtils::SpectrumType enum  
  NLOHMANN_JSON_SERIALIZE_ENUM(SpecUtils::SpectrumType, {
      {SpecUtils::SpectrumType::Foreground, "Foreground"},
      {SpecUtils::SpectrumType::Background, "Background"},
      {SpecUtils::SpectrumType::SecondForeground, "Secondary"},
  })
  
  double rount_to_hundredth(double val){ return 0.01*std::round(100.0*val); }

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

  
  void from_json(const json& j, AnalystChecks::FitPeakOptions& p) {
    p.energy = j.at("energy").get<double>();
    
    p.addToUsersPeaks = true;
    if( j.contains("addToUsersPeaks") )
      p.addToUsersPeaks = j.at("addToUsersPeaks").get<bool>();
    
    std::string specTypeStr = j.value("specType", std::string());
    if (specTypeStr.empty() || specTypeStr == "Foreground") {
      p.specType = SpecUtils::SpectrumType::Foreground;
    } else if (specTypeStr == "Background") {
      p.specType = SpecUtils::SpectrumType::Background;
    } else if (specTypeStr == "Secondary") {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    } else {
      throw std::runtime_error("Invalid spectrum type: " + specTypeStr);
    }
    
    p.userSession = j.value("userSession", std::optional<std::string>{});
    p.source = j.value("source", std::optional<std::string>{});
  }
  
  /*
  void to_json(json& j, const PeakDef& p) {
    j = json{{"lowerEnergy", rount_to_hundredth(p.lowerX())}, {"upperEnergy", rount_to_hundredth(p.upperX())} };

    if( p.type() == PeakDef::DefintionType::GaussianDefined )
    {
      j["fwhm"] = p.fwhm();
      j["energy"] = rount_to_hundredth(p.mean());
      j["amplitude"] = rount_to_hundredth(p.peakArea());
      if( p.amplitudeUncert() > 0.0 )
        j["numSigma"] = rount_to_hundredth( p.peakArea() / p.amplitudeUncert() );
    }else if( p.type() == PeakDef::DefintionType::DataDefined )
    {
      j["type"] = "DataDefined";
    }
    
    const uintptr_t ptr_val = reinterpret_cast<uintptr_t>(p.continuum().get());
    j["roiID"] = static_cast<uint64_t>( ptr_val );
  }
   */
  
  void to_json( json &peak_json, const shared_ptr<const PeakDef> &peak, const shared_ptr<const SpecUtils::Measurement> &meas ){
    if( peak->type() == PeakDef::DefintionType::GaussianDefined )
    {
      peak_json["fwhm"] = rount_to_hundredth(peak->fwhm());
      peak_json["energy"] = rount_to_hundredth(peak->mean());
      peak_json["amplitude"] = rount_to_hundredth(peak->peakArea());
      if( peak->amplitudeUncert() > 0.0 )
      {
        peak_json["numSigma"] = rount_to_hundredth( peak->peakArea() / peak->amplitudeUncert() );
        peak_json["amplitudeUncert"] = rount_to_hundredth(peak->amplitudeUncert());
      }
      if( meas && (meas->live_time() > 0.0) )
      {
        const double cps = peak->peakArea() / meas->live_time();
        peak_json["cps"] = rount_to_hundredth( cps );
        if( peak->amplitudeUncert() > 0.0 )
        {
          const double cpsUncert = cps * peak->amplitudeUncert() / peak->peakArea();
          peak_json["cpsUncert"] = rount_to_hundredth( cpsUncert );
        }
      }
    }else if( peak->type() == PeakDef::DefintionType::DataDefined )
    {
      peak_json["type"] = "DataDefined";
      // TODO: add in more infor using `meas` here
    }
    
    
    if( const SandiaDecay::Nuclide * const nuc = peak->parentNuclide() )
    {
      auto &src = peak_json["source"];
      src["nuclide"] = nuc->symbol;
      const SandiaDecay::Transition * const trans = peak->nuclearTransition();
      if( trans )
      {
        src["transition"] = (trans->parent ? trans->parent->symbol : string("null"))
        + "->" + (trans->child ? trans->child->symbol : string("null"));
      }
      const SandiaDecay::RadParticle * const particle = peak->decayParticle(); //may be null
      if( particle )
      {
        src["photonType"] = SandiaDecay::to_str( particle->type );
        src["photonEnergy"] = particle->energy;
      }
      
      const char *gamma_type = nullptr;
      switch( peak->sourceGammaType() )
      {
        case PeakDef::NormalGamma:
          gamma_type = "gamma";
          break;
        case PeakDef::AnnihilationGamma:
          gamma_type = "Annih.";
          break;
        case PeakDef::SingleEscapeGamma:
          gamma_type = "S.E.";
          break;
        case PeakDef::DoubleEscapeGamma:
          gamma_type = "D.E.";
          break;
        case PeakDef::XrayGamma:
          gamma_type = "x-ray";
          break;
      }//
      
      if( gamma_type )
        src["photonType"] = gamma_type;
      
      src["energy"] = peak->gammaParticleEnergy();
    }else if( const SandiaDecay::Element * const el = peak->xrayElement() )
    {
      peak_json["element"] = el->symbol;
      peak_json["source"]["photonType"] = "x-ray";
      peak_json["source"]["energy"] = peak->xrayEnergy();
    }else if( const ReactionGamma::Reaction * const rctn = peak->reaction() )
    {
      peak_json["source"]["reaction"] = rctn->name();
      peak_json["source"]["energy"] = peak->reactionEnergy();
    }
  }
  
  void to_json( json &roi_json,
               const shared_ptr<const PeakContinuum> &cont,
               const vector<shared_ptr<const PeakDef>> &peaks,
               const shared_ptr<const SpecUtils::Measurement> &meas ){
    roi_json = json{
      {"lowerEnergy", rount_to_hundredth(cont->lowerEnergy())},
      {"upperEnergy", rount_to_hundredth(cont->upperEnergy())},
      {"continuumType", PeakContinuum::offset_type_str(cont->type()) }
      //("continuumCounts", cont->offset_integral( cont->lowerEnergy(), cont->upperEnergy(), dataH ) }
    };
    
    for( const shared_ptr<const PeakDef> &peak : peaks ) {
      json peak_json;
      
      to_json( peak_json, peak, meas );
      
      roi_json["peaks"].push_back( peak_json );
    }
  }//void to_json( json &roi_json, const shared_ptr<const PeakContinuum> &cont, const vector<shared_ptr<const PeakDef>> &peaks )

  
  void to_json( json &peak_rois,
               const std::vector<std::shared_ptr<const PeakDef>> &peaks,
               const shared_ptr<const SpecUtils::Measurement> &meas ){
    peak_rois = json::array();
    
    vector<pair<shared_ptr<const PeakContinuum>,vector<shared_ptr<const PeakDef>>>> rois;
    for( const shared_ptr<const PeakDef> &peak : peaks ) {
      auto pos = std::find_if( begin(rois), end(rois), [&peak]( const auto &val){ return val.first == peak->continuum(); } );
      if( pos == end(rois) )
        rois.push_back( make_pair(peak->continuum(), vector<shared_ptr<const PeakDef>>(1,peak) ) );
      else
        pos->second.push_back( peak );
    }
    
    
    for( size_t roi_index = 0; roi_index < rois.size(); ++roi_index )
    {
      const shared_ptr<const PeakContinuum> &cont = rois[roi_index].first;
      const vector<shared_ptr<const PeakDef>> &peaks = rois[roi_index].second;
      
      json roi_json;
      
      to_json( roi_json, cont, peaks, meas );
      
      roi_json["roiID"] = static_cast<int>(roi_index);
      
      peak_rois.push_back(roi_json);
    }
  }
  
  
  void to_json(json& j,
               const AnalystChecks::DetectedPeakStatus& p,
               const shared_ptr<const SpecUtils::Measurement> &meas ) {
    json peak_rois;
    to_json( peak_rois, p.peaks, meas );
    
    j = json{{"userSession", p.userSession},
      {"rois", peak_rois}};
  }//void to_json(json& j, const AnalystChecks::DetectedPeakStatus& p) {
  
  void to_json(json& j, const AnalystChecks::FitPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {
    
    const std::shared_ptr<const PeakDef> &fitPeak = p.fitPeak;
    const std::vector<std::shared_ptr<const PeakDef>> &peaksInRoi = p.peaksInRoi;
    
    j = json{{"userSession", p.userSession}};
    
    if( fitPeak )
    {
      json roi_json;
      
      to_json( roi_json, fitPeak->continuum(), peaksInRoi, meas );
      
      j["roi"] = roi_json;
      
      //json peak_json;
      //to_json( peak_json, fitPeak );
      j["fitPeakEnergy"] = rount_to_hundredth(fitPeak->mean());
    }else
    {
      j["error"] = "No peak fit.";
    }
  }//void to_json(json& j, const AnalystChecks::DetectedPeakStatus& p)
  
  
  void from_json(const json& j, AnalystChecks::GetUserPeakOptions& p) {
    p.userSession = j.value("userSession", std::optional<std::string>{});
    
    std::string specTypeStr = j.value("specType", std::string());
    if (specTypeStr.empty() || specTypeStr == "Foreground") {
      p.specType = SpecUtils::SpectrumType::Foreground;
    } else if (specTypeStr == "Background") {
      p.specType = SpecUtils::SpectrumType::Background;
    } else if (specTypeStr == "Secondary") {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    } else {
      throw std::runtime_error("Invalid spectrum type: " + specTypeStr);
    }
  }
   
  
  
  // Call the AnalystChecks function to actually get the peaks
  void to_json(json& j, const AnalystChecks::GetUserPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {
    json peak_rois;
    to_json( peak_rois, p.peaks, meas );
    
    j = json{{"userSession", p.userSession},
      {"rois", peak_rois}};
  }
}//namespace

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
    "Returns all Regions Of Interest (ROI) with peaks detected by automated peak search. For ROI gives lower and upper energies, and for each peak it gives energy (in keV), FWHM, amplitude (area), and statistical significance (numSigma); if the peak is associated with a source, will also give information on that. Does not add peaks to the user peaks.",
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
            "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
          }
      },
      "required": ["specType"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executePeakDetection(params, interspec);
    }
  });
  
  // Register detected_peaks tool
  registerTool({
    "fit_peak",
    "Fit and add a peak to the users peaks, at approximately the specified energy, optionally associating a source with it.  Returns Region Of Interest that was either created, or the peak was added to.  If fit failed, reason will be described in 'error' field.",
    json::parse(R"({
      "type": "object",
      "properties": {
          "energy": { 
            "type": "number", 
            "description": "Approximate energy (in keV) to look for a peak to fit." 
          },
          "source": {
            "type": "string",
            "description": "Optional: The parent nuclide (ex U235, I131, Ba133) or x-ray flourescense element (ex Pb, U, W) or nuclear reaction (ex H(n,g)) assigned to the peak."
          },
          "specType": { 
            "type": "string", 
            "description": "Optional: Which displayed spectrum to search for peaks in; if not specified will use foreground (which is what user usually wants).", 
            "enum": ["Foreground", "Background", "Secondary"] 
          },
          "addToUsersPeaks": {
            "type": "boolean",
            "description": "Optional: if fit peak should be added to users peaks; defaults to true."
          },
          "userSession": { 
            "type": "string", 
            "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
          }
      },
      "required": ["energy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executePeakFit(params, interspec);
    }
  });
  
  
  // Register detected_peaks tool
  registerTool({
    "get_user_peaks",
    "Gets the user peaks that have either been manually fit by the user, or by the 'fit_peak' tool call, or similar.",
    json::parse(R"({
      "type": "object",
      "properties": {
          "specType": { 
            "type": "string", 
            "description": "Optional: Which displayed spectrum to search for peaks in; if not specified will use foreground (which is what user usually wants).", 
            "enum": ["Foreground", "Background", "Secondary"] 
          },
          "userSession": { 
            "type": "string", 
            "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
          }
      }
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetUserPeaks(params, interspec);
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
          },
          "userSession": { 
            "type": "string", 
            "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
          }
      },
      "required": ["specType"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetSpectrumInfo(params, interspec);
    }
  });

  registerTool({
    "primary_gammas_for_nuclide",
    "Get the most likely energies a peak will be detected at for a nuclide, x-ray element, or nuclear reaction. Returns one to a few energies (in keV).  This is only a rough guess at the most prominent or unique peaks for the source; not detecting a peak at these energies does not rule out a source, but it makes the source a candidate.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "type": "string",
          "description": "The nuclide, x-ray element, or nuclear reaction to get the characteristic gammas for."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetCharacteristicGammasForNuclide(params);
    }
  });


  registerTool({
    "nuclides_with_primary_gammas_in_energy_range",
    "Get the commonly encountered field nuclides, x-ray elements, or nuclear reactions with primary gammas in the specified energy range.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "lowerEnergy": {
          "type": "number",
          "description": "The lower energy (in keV) to search for primary gammas."
        },
        "upperEnergy": {
          "type": "number",
          "description": "The upper energy (in keV) to search for primary gammas."
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        }
      },
      "required": ["lowerEnergy", "upperEnergy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
    }
  });

  registerTool({
    "nuclides_with_primary_gammas_near_energy",
    "Get the commonly encountered field nuclides, x-ray elements, or nuclear reactions with primary gammas near the specified energy.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "energy": {
          "type": "number",
          "description": "The energy (in keV) to search for primary gammas."
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        }
      },
      "required": ["energy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
    }
  });

  registerTool({
    "associated_nuclides",
    "Gets other nuclides that are commonly detected along with the specified nuclide, or nuclides that might be mis-identified as the specified nuclide.  You you should check for these associated nuclides being present or not when the specified nuclide is observed.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "type": "string",
          "description": "The nuclide to get associated nuclides for."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetAssociatedNuclides(params);
    }
  });

  registerTool({
    "analyst_notes_for_nuclide",
    "Gets analyst notes for the specified nuclide, such as when the nuclide is used, typical activity ranges, etc.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "type": "string",
          "description": "The nuclide to get analyst notes for."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetNuclideAnalystNotes(params);
    }
  });
  
  registerTool({
    "nuclide_info",
    "Gets nuclear data information about the specified nuclide, such as half-life, decay modes, etc.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "type": "string",
          "description": "The nuclide to get data for."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetNuclideInfo(params);
    }
  });

  registerTool({
    "nuclide_decay_chain",
    "Returns a list of all progeny nuclides, of the specified nuclide, and thier decay modes and branching ratios.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "type": "string",
          "description": "The nuclide to get data for."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetNuclideDecayChain(params);
    }
  });


  registerTool({
    "loaded_spectra",
    "Returns array of the currently loaded spectra from [\"Foreground\", \"Background\", \"Secondary\"].  If no spectra are loaded, returns empty array.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "userSession": { 
            "type": "string", 
            "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
          }
      }
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetLoadedSpectra(params, interspec);
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
  if( !interspec )
    throw std::runtime_error("No InterSpec session available");
  
  // Parse parameters into DetectedPeaksOptions
  AnalystChecks::DetectedPeaksOptions options;
  from_json(params, options);
  
  shared_ptr<const SpecUtils::Measurement> meas;
  if( interspec )
    meas = interspec->displayedHistogram( options.specType );
  
  // Call the AnalystChecks function to perform the actual peak detection
  AnalystChecks::DetectedPeakStatus result = AnalystChecks::detected_peaks(options, interspec);
  
  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result, meas );

  return result_json;
}
  
  
nlohmann::json ToolRegistry::executePeakFit(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  // Parse parameters into DetectedPeaksOptions
  AnalystChecks::FitPeakOptions options;
  from_json(params, options);
  
  // Call the AnalystChecks function to perform the actual peak detection
  const AnalystChecks::FitPeakStatus result = AnalystChecks::fit_user_peak( options, interspec );
  
  shared_ptr<const SpecUtils::Measurement> meas;
  if( interspec )
    meas = interspec->displayedHistogram( options.specType );
  
  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result, meas );
  
  return result_json;
}//nlohmann::json ToolRegistry::executePeakFit(const nlohmann::json& params, InterSpec* interspec)

  
nlohmann::json ToolRegistry::executeGetUserPeaks(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  // Parse parameters into DetectedPeaksOptions
  AnalystChecks::GetUserPeakOptions options;
  from_json(params, options);
  
  // Call the AnalystChecks function to actually get the peaks
  const AnalystChecks::GetUserPeakStatus result = AnalystChecks::get_user_peaks( options, interspec);
  
  shared_ptr<const SpecUtils::Measurement> meas;
  if( interspec )
    meas = interspec->displayedHistogram( options.specType );
  
  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result, meas );
  
  return result_json;
}//nlohmann::json executeGetUserPeaks(const nlohmann::json& params, InterSpec* interspec)
  
json ToolRegistry::executeGetSpectrumInfo(const json& params, InterSpec* interspec) {
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

  const set<int> &displayedSamples = interspec->displayedSamples(specType);
  
  shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(specType);
  if (!spectrum)
    throw std::runtime_error("No spectrum displayed for " + specTypeStr + " spectrum");
  
  json result;
  result["specType"] = specTypeStr;
  result["detectorName"] = spectrum->detector_name();
  result["fileName"] = meas->filename();
  result["liveTime"] = spectrum->live_time();
  result["realTime"] = spectrum->real_time();
  result["startTime"] = SpecUtils::to_iso_string( spectrum->start_time() );
  result["numChannels"] = spectrum->num_gamma_channels();
  //result["energyCalibration"] = spectrum->calibration_coeffs();
  result["displayedSamples"] = displayedSamples;
  result["includesNeutron"] = spectrum->contained_neutron();
  if( spectrum->contained_neutron() )
  {
    result["neutronCounts"] = spectrum->neutron_counts_sum();
    result["neutronLiveTime"] = spectrum->neutron_live_time();
    result["neutronCPS"] = spectrum->neutron_counts_sum() / spectrum->neutron_live_time();
  }
  if( !spectrum->title().empty() )
    result["title"] = spectrum->title();
  if( meas->detector_type() != SpecUtils::DetectorType::Unknown )
    result["detectorType"] = SpecUtils::detectorTypeToString(meas->detector_type());
  if( !meas->manufacturer().empty() )
    result["detectorManufacturer"] = meas->manufacturer();
  if( !meas->instrument_model().empty() )
    result["detectorModel"] = meas->instrument_model();
  if( !meas->instrument_id().empty() )
    result["detectorSerialNumber"] = meas->instrument_id();
  if( spectrum->has_gps_info() )
  {
    result["gpsLatitude"] = spectrum->latitude();
    result["gpsLongitude"] = spectrum->longitude();
  }

  result["minimumGammaEnergy"] = spectrum->gamma_energy_min();
  result["maximumGammaEnergy"] = spectrum->gamma_energy_max();

  if (spectrum->gamma_counts() && !spectrum->gamma_counts()->empty()) {
    result["totalCounts"] = spectrum->gamma_count_sum();
  }
  
  return result;
}//json ToolRegistry::executeGetSpectrumInfo(const json& params, InterSpec* interspec)
  
nlohmann::json ToolRegistry::executeGetCharacteristicGammasForNuclide( const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();
  json result;
  result["nuclide"] = nuclide;
  result["characteristicGammas"] = AnalystChecks::get_characteristic_gammas( nuclide );
  return result;
}

nlohmann::json ToolRegistry::executeGetLoadedSpectra( const nlohmann::json& params, InterSpec* interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  vector<SpecUtils::SpectrumType> specTypes = { SpecUtils::SpectrumType::Foreground, SpecUtils::SpectrumType::Background, SpecUtils::SpectrumType::SecondForeground };
  vector<string> loadedSpectra;
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground ) )
    loadedSpectra.push_back("Foreground");
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::Background ) )
    loadedSpectra.push_back("Background");
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::SecondForeground ) )
    loadedSpectra.push_back("Secondary");
  
  return json{ loadedSpectra };
}
  
nlohmann::json ToolRegistry::executeGetNuclidesWithCharacteristicsInEnergyRange( const nlohmann::json& params, InterSpec* interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  vector<variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>> result;
  if( params.contains("lowerEnergy") && params.contains("upperEnergy") )
  {
    const double lower_energy = params.at("lowerEnergy").get<double>();
    const double upper_energy = params.at("upperEnergy").get<double>();
    result = AnalystChecks::get_nuclides_with_characteristics_in_energy_range( lower_energy, upper_energy, interspec );
  }else if( params.contains("energy") )
  {
    const double energy = params.at("energy").get<double>();
    result = AnalystChecks::get_characteristics_near_energy( energy, interspec );
  }else
  {
    throw std::runtime_error("Missing lowerEnergy, upperEnergy, or energy parameter");
  }
  
  json result_json = json::array();
  for( const auto &item : result )
  {
    if( std::holds_alternative<const SandiaDecay::Nuclide *>(item) )
    {
      result_json.push_back( std::get<const SandiaDecay::Nuclide *>(item)->symbol );
    }else if( std::holds_alternative<const SandiaDecay::Element *>(item) )
    {
      result_json.push_back( std::get<const SandiaDecay::Element *>(item)->symbol );
    }else if( std::holds_alternative<const ReactionGamma::Reaction *>(item) )
    {
      result_json.push_back( std::get<const ReactionGamma::Reaction *>(item)->name() );
    }
  }
  return result_json;
}
  

nlohmann::json ToolRegistry::executeGetAssociatedNuclides( const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();

  try
  {
    const shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> info_db = MoreNuclideInfo::MoreNucInfoDb::instance();
    if( !info_db )
      throw std::runtime_error("No MoreNucInfoDb instance available");
    
    const MoreNuclideInfo::NucInfo *nuc_info = info_db->info(nuclide);
    if( !nuc_info )
      throw runtime_error( "No info for " + nuclide );
    
    return json{{"status", "success"}, {"nuclide", nuclide}, {"associatedNuclides", nuc_info->m_associated} };
  }catch( const std::exception &e )
  {
    return json{{"status", "failed"}, {"nuclide", nuclide}, {"error", e.what()}};
  }
  assert( 0 );
  return {};
}//nlohmann::json executeGetAssociatedNuclides(const nlohmann::json& params, InterSpec* interspec)


nlohmann::json ToolRegistry::executeGetNuclideAnalystNotes( const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();
  
  try
  {
    const shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> info_db = MoreNuclideInfo::MoreNucInfoDb::instance();
    if( !info_db )
      throw std::runtime_error("No MoreNucInfoDb instance available");
    
    const MoreNuclideInfo::NucInfo *nuc_info = info_db->info(nuclide);
    if( !nuc_info || nuc_info->m_notes.empty() )
      throw runtime_error( "No info for " + nuclide );
    
    return json{{"status", "success"}, {"nuclide", nuclide}, {"notes", nuc_info->m_notes} };
  }catch( const std::exception &e )
  {
    return json{{"status", "failed"}, {"nuclide", nuclide}, {"error", e.what()}};
  }
  assert( 0 );
  return {};
}//nlohmann::json executeGetNuclideAnalystNotes(const nlohmann::json& params, InterSpec* interspec)

  
nlohmann::json ToolRegistry::executeGetNuclideInfo(const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();
  
  nlohmann::json result;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    result["error"] = "Could not initialize nuclide DecayDataBase.";
    return result;
  }
  
  const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
  if( !nuc )
  {
    result["error"] = "Nuclide '" + nuclide + "' is not a valid nuclide.";
    return result;
  }
  
  try
  {
    const nlohmann::json associated = ToolRegistry::executeGetAssociatedNuclides( params );
    if( associated.contains("associatedNuclides") )
      result["associatedNuclides"] = associated["associatedNuclides"];
  }catch( std::exception & )
  {
  }
  
  result["symbol"] = nuc->symbol;
  result["atomicNumber"] = static_cast<int>(nuc->atomicNumber);
  result["massNumber"] = static_cast<int>(nuc->massNumber);
  result["isomerNumber"] = static_cast<int>(nuc->isomerNumber);
  result["atomicMass"] = nuc->atomicMass;
  result["halfLife"] = PhysicalUnits::printToBestTimeUnits(nuc->halfLife, 6);
  if( nuc->canObtainPromptEquilibrium() )
    result["promptEquilibriumHalfLife"] = PhysicalUnits::printToBestTimeUnits(nuc->promptEquilibriumHalfLife(), 6);
  if( nuc->canObtainSecularEquilibrium() )
    result["secularEquilibriumHalfLife"] = PhysicalUnits::printToBestTimeUnits(nuc->secularEquilibriumHalfLife(), 6);
  result["atomsPerGram"] = nuc->atomsPerGram();
  result["activityPerGram"] = nuc->activityPerGram();
  result["isStable"] = nuc->isStable();
  result["decaysToStableChildren"] = nuc->decaysToStableChildren();
  result["defaultAge"] = PhysicalUnits::printToBestTimeUnits( PeakDef::defaultDecayTime(nuc), 6 );
  
  for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
  {
    result["decays"].push_back( json{{"child", trans->child ? trans->child->symbol : ""},
      {"branchingRatio",trans->branchRatio},
      {"decayType", SandiaDecay::to_str(trans->mode)}}
    );
  }
  
  return result;
}//nlohmann::json ToolRegistry::executeGetNuclideInfo(const nlohmann::json& params )

nlohmann::json ToolRegistry::executeGetNuclideDecayChain(const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();

  nlohmann::json result;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
  {
    result["error"] = "Could not initialize nuclide DecayDataBase.";
    return result;
  }
  
  const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
  if( !nuc )
  {
    result["error"] = "Nuclide '" + nuclide + "' is not a valid nuclide.";
    return result;
  }
  
  const vector<const SandiaDecay::Nuclide *> descendants = nuc->descendants();
  for( const SandiaDecay::Nuclide * const kid : descendants )
  {
    nlohmann::json kid_info;
    kid_info["nuclide"] = kid->symbol;
    kid_info["halfLife"] = PhysicalUnits::printToBestTimeUnits(kid->halfLife, 6);
    
    for( const SandiaDecay::Transition *trans : nuc->decaysToChildren )
    {
      kid_info["decays"].push_back( json{{"child", trans->child ? trans->child->symbol : ""},
        {"branchingRatio",trans->branchRatio},
        {"decayType", SandiaDecay::to_str(trans->mode)}}
      );
    }
    
    result.push_back( std::move(kid_info) );
  }//for( const SandiaDecay::Nuclide * const kid : descendants )
  
  return result;
}//nlohmann::json ToolRegistry::executeGetNuclideDecayChain(const nlohmann::json& params )
  
} // namespace LlmTools

#endif // USE_LLM_INTERFACE
