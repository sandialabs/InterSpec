#include "InterSpec_config.h"
#include "InterSpec/LlmToolRegistry.h"

#if( USE_LLM_INTERFACE )

#include <sstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/MoreNuclideInfo.h"
#include "InterSpec/ExternalRidResult.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"

#include <Wt/WApplication>
#include <Wt/Dbo/Dbo>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
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

  /** Some LLMs will give number values as strings, so this function will check types and return the correct answer.

   @param parent The parent JSON object.
   @param name The field name of the number to parse.

   An exception will be thrown if cant be converted to int.
   */
  double get_number( const json& parent, const string &name )
  {
    if( !parent.contains(name) )
      throw runtime_error( "'" + name + "' parameter must be specified." );

    if( parent[name].is_number() )
      return parent.at(name).get<double>();

    if( parent[name].is_string() )
    {
      string strval = parent.at(name).get<string>();
      double val;
      if( !(stringstream(strval) >> val) )
        throw runtime_error( "'" + name + "' parameter must be a number." );
      return val;
    }//if( parent[name].is_string() )

    throw runtime_error( "'" + name + "' parameter must be a number." );
  }//double get_number( const json& parent, const string &name )

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

    p.energy = get_number( j, "energy" );
    
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

  void from_json(const json& j, AnalystChecks::FitPeaksForNuclideOptions& p) {
    const json& nuclideParam = j.at("nuclide");
    if (nuclideParam.is_string()) {
      p.nuclides = {nuclideParam.get<std::string>()};
    } else if (nuclideParam.is_array()) {
      p.nuclides = nuclideParam.get<std::vector<std::string>>();
    } else {
      throw std::runtime_error("Invalid nuclide parameter: must be string or array of strings");
    }
    
    p.doNotAddPeaksToUserSession = j.value("doNotAddPeaksToUserSession", false);
  }

  void to_json(json& j, const AnalystChecks::SpectrumCountsInEnergyRange::CountsWithComparisonToForeground& c) {
    j = json{
      {"counts", c.counts},
      {"cps", c.cps},
      {"numSigmaCpsRelForeground", c.num_sigma_rel_foreground}
    };
  }

  void to_json(json& j, const AnalystChecks::SpectrumCountsInEnergyRange& c) {
    j = json{
      {"lowerEnergy", c.lower_energy},
      {"upperEnergy", c.upper_energy},
      {"foregroundCounts", c.foreground_counts},
      {"foregroundCps", c.foreground_cps}
    };

    if (c.background_info.has_value()) {
      json background_json;
      to_json(background_json, c.background_info.value());
      j["backgroundInfo"] = background_json;
    }

    if (c.secondary_info.has_value()) {
      json secondary_json;
      to_json(secondary_json, c.secondary_info.value());
      j["secondaryInfo"] = secondary_json;
    }
  }


  void to_json(json& j, const DetectionLimitCalc::CurrieMdaResult& result) {
    j = json{
      {"gammaEnergy", result.input.gamma_energy},
      {"roiLowerEnergy", result.input.roi_lower_energy},
      {"roiUpperEnergy", result.input.roi_upper_energy},
      {"numLowerSideChannels", result.input.num_lower_side_channels},
      {"numUpperSideChannels", result.input.num_upper_side_channels},
      {"detectionProbability", result.input.detection_probability},
      {"additionalUncertainty", result.input.additional_uncertainty},
      {"firstPeakRegionChannel", static_cast<int>(result.first_peak_region_channel)},
      {"lastPeakRegionChannel", static_cast<int>(result.last_peak_region_channel)},
      {"peakRegionCountsSum", result.peak_region_counts_sum},
      {"estimatedPeakContinuumCounts", result.estimated_peak_continuum_counts},
      {"estimatedPeakContinuumUncert", result.estimated_peak_continuum_uncert},
      {"decisionThreshold", result.decision_threshold},
      {"detectionLimit", result.detection_limit},
      {"sourceCounts", result.source_counts},
      {"lowerLimit", result.lower_limit},
      {"upperLimit", result.upper_limit},
      {"peakPresentInData", (result.source_counts > result.decision_threshold)}
    };
  }
   
  
  
  // Call the AnalystChecks function to actually get the peaks
  void to_json(json& j, const AnalystChecks::GetUserPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {
    json peak_rois;
    to_json( peak_rois, p.peaks, meas );
    
    j = json{{"userSession", p.userSession},
      {"rois", peak_rois}};
  }

  void to_json(json& j, const AnalystChecks::FitPeaksForNuclideStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {
    json peak_rois;
    to_json( peak_rois, p.fitPeaks, meas );
    
    j = json{{"rois", peak_rois}};
  }
  
  /** Returns a `nlohmann::json::array` containing the source catagories (Medical, Industrial, NORM, etc) that a source (Nuclide, Element, or Reaction),
   belongs to.
   
   The `src` must be a `const SandiaDecay::Element *`, `const ReactionGamma::Reaction *` , or `const SandiaDecay::Nuclide *`
   */
  template<class T>
  typename std::enable_if<
    std::is_same<T, const SandiaDecay::Nuclide *>::value ||
    std::is_same<T, const SandiaDecay::Element *>::value ||
    std::is_same<T, const ReactionGamma::Reaction *>::value,
    nlohmann::json
  >::type source_categories( T src, InterSpec *interspec )
  {
    nlohmann::json sourceCatagories = nlohmann::json::array();
    
    IsotopeSearchByEnergy * const search = interspec ? interspec->nuclideSearch() : nullptr;
    if( !search || !src )
      return sourceCatagories;
    
    const std::vector<IsotopeSearchByEnergy::NucSearchCategory> &categories = search->search_categories();
    
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_medical_category_key, categories) )
      sourceCatagories.push_back( "Medical" );
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_industrial_category_key, categories) )
      sourceCatagories.push_back( "Industrial" );
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_norm_category_key, categories) )
      sourceCatagories.push_back( "NORM" );
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_snm_category_key, categories) )
      sourceCatagories.push_back( "SNM" );
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_common_category_key, categories) )
      sourceCatagories.push_back( "Common" );
    if( IsotopeSearchByEnergy::is_in_category(src, IsotopeSearchByEnergy::sm_fission_category_key, categories) )
      sourceCatagories.push_back( "Fission" );
    
    return sourceCatagories;
  }//source_categories(...)
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
    "get_analysis_peaks",
    "Gets the peaks to use for further analysis that have either been manually fit by the user, or by the 'fit_peak' tool call, or similar.  These peaks tracked by an internal state.",
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
    "primary_gammas_for_source",
    "Get the most likely energies a peak will be detected at for a source (nuclide, x-ray element, or nuclear reaction). Returns one to a few energies (in keV).  This is only a rough guess at the most prominent or unique peaks for the source; not detecting a peak at these energies does not rule out a source, but it makes the source a candidate.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "source": {
          "type": "string",
          "description": "The source (nuclide, x-ray element, or nuclear reaction) to get the characteristic gammas for."
        }
      },
      "required": ["source"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetCharacteristicGammasForSource(params);
    }
  });

// TODO: combine sources_with_primary_gammas_in_energy_range and sources_with_primary_gammas_near_energy, and make them return the nuclide search tool results, not just primary nuclides
  registerTool({
    "sources_with_primary_gammas_in_energy_range",
    "Get the commonly encountered field nuclides, x-ray elements, or reactions with primary gammas in the specified energy range.",
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
    "sources_with_primary_gammas_near_energy",
    "Get the commonly encountered field nuclides, x-ray elements, or reactions with primary gammas near the specified energy.",
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

#if( !INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
  registerTool({
    "sources_associated_with_source",
    "Gets other sources (nuclides, reactions, x-rays) that are commonly detected along with the specified source, or sources that might be mis-identified as the specified source.  You you should check for these associated sources being present or not when the specified source is observed.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "source": {
          "type": "string",
          "description": "The source to get associated sources for."
        }
      },
      "required": ["source"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetAssociatedSources(params);
    }
  });

  registerTool({
    "analyst_notes_for_source",
    "Gets analyst notes for the specified source (nuclide, x-ray, reaction), such as when the source is used/seen, typical activity ranges of nuclides, etc.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "source": {
          "type": "string",
          "description": "The source to get analyst notes for."
        }
      },
      "required": ["source"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetSourceAnalystNotes(params);
    }
  });
#endif //!INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO
  
  registerTool({
    "source_info",
#if( !INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
    "Gets nuclear data information about the specified source (nuclide, x-ray element, reaction), such as half-life, decay modes, etc.",
#else
    "Gets relevant information about the specified source (nuclide, x-ray element, reaction), such as half-life, decay modes, analyst notes (such as when the source is used/seen, typical activity ranges of nuclides, etc.), associated sources (i.e. sources that are commonly detected along with the specified source, or sources that might be mis-identified as the specified source - you should check for these associated sources being present or not when the specified source is observed), source catagories (medical, NORM, industrial, etc), and other nuclear data.",
#endif
    json::parse(R"({
      "type": "object",
      "properties": {
        "source": {
          "type": "string",
          "description": "The source to get data for."
        }
      },
      "required": ["source"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetSourceInfo(params, interspec);
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
    "automated_isotope_id_results",
    "Get the isotope ID from the avaiable automated ID algorithms. Will return the detection systems on-board nuclide ID results, if present, as well as the GADRAS Full Spectrum Isotope ID results if the app is configured to get.",
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
      return executeGetAutomatedRiidId(params, interspec);
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

  registerTool({
    "fit_peaks_for_nuclide",
    "Fits peaks for the one or more specified nuclide(s). Returns the peaks that were fit.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "nuclide": {
          "anyOf": [
            {"type": "string"},
            {"type": "array", "items": {"type": "string"}}
          ],
          "description": "The nuclide or array of nuclides to fit peaks for (ex U235, I131, Ba133)."
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        },
        "doNotAddPeaksToUserSession": {
          "type": "boolean",
          "description": "Optional: if true, fitted peaks will not be added to the user's peak collection; defaults to false."
        }
      },
      "required": ["nuclide"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return ToolRegistry::executeFitPeaksForNuclide(params, interspec);
    }
  });

  registerTool({
    "get_counts_in_energy_range",
    "Get the counts in the spectra for the specified energy range as well as compares statistical significance of differences between foreground and background and/or secondary count-rates, if those spectra are loaded.  This function can be used to check if the count rate of an energy range is elevated in the foreground relative to the background, especially above the 2614 keV peak.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "lowerEnergy": {
          "type": "number",
          "description": "The lower energy bound in keV."
        },
        "upperEnergy": {
          "type": "number",
          "description": "The upper energy bound in keV."
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        }
      },
      "required": ["lowerEnergy", "upperEnergy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return ToolRegistry::executeGetCountsInEnergyRange(params, interspec);
    }
  });

  registerTool({
    "get_expected_fwhm",
    "Calculate the expected Full Width at Half Maximum (FWHM) for a peak at the specified energy, based on the detector response function, detected peaks, or detector type. This is useful for determining the expected width of a peak for ROI calculations.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "energy": {
          "type": "number",
          "description": "The energy (in keV) to calculate the expected FWHM for."
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        }
      },
      "required": ["energy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return ToolRegistry::executeGetExpectedFwhm(params, interspec);
    }
  });

  registerTool({
    "currie_mda_calc",
    "Calculate Minimum Detectable Activity (MDA) and detection confidence intervals using Currie-style (ISO 11929) methodology. Determines if a peak is detectable at a specified energy. Returns decision threshold, detection limit, and confidence intervals - as well as if it looks like a peak is at the energy (see `peakPresentInData` in output).",
    json::parse(R"({
      "type": "object",
      "properties": {
        "energy": {
          "type": "number",
          "description": "The energy (in keV) of the photopeak for which the detection limit is being calculated."
        },
        "detectionProbability": {
          "type": "number",
          "description": "Optional: Detection probability (confidence level), typically 0.95 for 95% confidence. Defaults to 0.95.",
          "minimum": 0.5,
          "maximum": 0.999,
          "default": 0.95
        },
        "additionalUncertainty": {
          "type": "number",
          "description": "Optional: Additional relative uncertainty to include (e.g., from detector response function uncertainties). Defaults to 0.0.",
          "minimum": 0.0,
          "default": 0.0
        },
        "userSession": { 
          "type": "string", 
          "description": "Optional: the user session identifier.  If not specified, will use most recent session." 
        }
      },
      "required": ["energy"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return ToolRegistry::executeCurrieMdaCalc(params, interspec);
    }
  });


  registerTool({
    "source_photons",
    "Returns a list of energy/intensity pairs associated with a nuclide, elemental fluorescence x-rays, or reaction.  The energies are expressed in keV, and intensities are given as photons per source Becquerel for nuclides, or normalized to 1 for flourescence or reactions.  Results will be sorted by energy.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "Source": {
          "type": "string",
          "description": "The nuclide, element, or reaction for which to retrieve decay product data. Examples: 'U238' (nuclide), 'Pb' (element), 'H(n,g)' (reaction)."
        },
        "Age": {
          "type": "string",
          "description":  "The age of the nuclide at which to return decay products. Applicable only for nuclides; specifying this for elements or reactions will result in an error. If not specified for nuclides, a default age will be used. Examples: '23.25 y', '23 years 60 days', '0s', '5 half-lives'. An age of '0 seconds' will return products from the specified nuclide only, while a positive value will include products from the nuclide's full decay chain, aged to the specified amount."
        },
        "MaxResults": {
          "type": "integer",
          "minimum": 1,
          "description": "The maximum number of energy/intensity pairs to return; if the source has more than this number of photons, then information for this number of the largest intensity photons will be returned.  If not specified, a value of 125 will be assumed."
        }
      },
      "required": ["Source"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetSourcePhotons(params);
    }
  });

  /*
  registerTool({
    "attenuation_factor_of_shielding",
    "Given an array of photon energies, returns the fraction of photons that pass through the shielding without being attenuated. For each energy, the returned value will be between 0 and 1, where 1 indicates no attenuation (all photons pass through) and 0 indicates complete attenuation (all photons interact).",
    json::parse(R"({
  "type": "object",
  "properties": {
    "Shielding": {
      "type": "object",
      "description": "Defines the shielding properties used to attenuate gammas and x-rays. Shielding must be specified as an object with one of the following formats:\n\n1. **Areal density and effective atomic number**: Provide an object with the keys `AD` (areal density in g/cm²) and `AN` (effective atomic number). Example: `{ \"AD\": 20.25, \"AN\": 26 }`.\n\n2. **Material and thickness**: Provide an object with the keys `Material` (element symbol or name) and `Thickness` (thickness in cm). Example: `{ \"Material\": \"Fe\", \"Thickness\": \"1.25cm\" }`.\n\nAvailable materials include element symbols or names and materials returned by the `get_materials` tool.",
      "properties": {
        "AD": {
          "type": "number",
          "description": "Areal density of the shielding material in g/cm²."
        },
        "AN": {
          "type": "number",
          "description": "Effective atomic number of the shielding material."
        },
        "Material": {
          "type": "string",
          "description": "The shielding material, specified as an element symbol or name. Example: 'Fe' or 'Iron'."
        },
        "Thickness": {
          "type": "string",
          "description": "The thickness of the shielding material, specified as a string with units. Example: '1.25cm'."
        }
      },
      "additionalProperties": false
    },
    "Energies": {
      "type": "array",
      "items": {
        "type": "number"
      },
      "description": "An array of photon energy values in keV to apply the attenuation calculation to. Each value must be a positive float representing the energy of a photon. Example: [511.0, 1460.8, 2614.5]. The returned attenuation fractions will correspond 1:1 to these input energies."
    }
  },
  "required": ["Shielding", "Energies"]
})"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetAttenuationOfShielding(params, interspec);
    }
  });
*/
  
  const string shielding_desc = R"("Optional: Defines the shielding used to attenuate gammas and x-rays. Shielding must be specified as an object with one of the following formats:\n\n1. **Areal density and effective atomic number**: Provide an object with the keys `AD` (areal density in g/cm²) and `AN` (effective atomic number). Example: `{ \"AD\": 20.25, \"AN\": 26 }`.\n\n2. **Material and thickness**: Provide an object with the keys `Material` (element symbol or name) and `Thickness` (thickness in cm). Example: `{ \"Material\": \"Fe\", \"Thickness\": \"1.25cm\" }`.\n\nAvailable materials include element symbols or names and materials returned by the `get_materials` tool.")";
  
  // Note: If using OpenAI's `strict: true` mode for structured outputs, `oneOf` is not supported - instead use `anyOf`.
  
  const string shielding_obj_properties = R"({
      "type": "object",
      "description": "Optional: Defines the shielding used to attenuate gammas and x-rays. Shielding must be specified in one of two formats.",
      "oneOf": [
        {
          "type": "object",
          "properties": {
            "AD": {
              "type": "number",
              "description": "Areal density of the shielding material in g/cm²."
            },
            "AN": {
              "type": "number",
              "description": "Effective atomic number of the shielding material."
            }
          },
          "required": ["AD", "AN"],
          "additionalProperties": false
        },
        {
          "type": "object",
          "properties": {
            "Material": {
              "type": "string",
              "description": "The shielding material, specified as an element symbol or material name. Available material names are returned by the `get_materials` tool. Example: 'Fe', which is equivalent to 'Iron', or 'glass_plate', which is returned from `get_materials`."
            },
            "Thickness": {
              "type": "string",
              "description": "The thickness of the shielding material, specified as a string with units. Example: '1.25cm'."
            }
          },
          "required": ["Material", "Thickness"],
          "additionalProperties": false
        }
      ]
    })";
  
  const string distance_obj_properties = R"({
    "type": "string",
    "description": "Optional: The distance the detector is from the center of the radioactive source.  If not specified, the gemoetric detection factor will not be accounted for. If a detector efficiency is not currently loaded for the foreground, or the detector efficiency is for a 'fixed geometry', then distance may not be specified. Example distances: '1.25 cm', '3 ft, 2 inches', etc." 
  })";
  
  const string par_schema = R"({
"type": "object",
"properties": {
  "Shielding": {
    "type": "array",
    "items": )" + shielding_obj_properties +
R"(,
    "description": "An array of shielding objects, applied in order from source to detector. Each shielding layer attenuates the photons that pass through it."
  },
  "Distance": )" + distance_obj_properties +
R"(,
  "IncludeAirAttenuation": {
    "type": "boolean",
    "description": "Optional: If true, includes attenuation from air between the source and detector. Only applies when Distance is specified. Defaults to false."
  },
  "Energies": {
    "type": "array",
    "items": {
      "type": "number"
    },
    "description": "An array of photon energy values in keV to apply the calculation to. Each value must be a positive float representing the energy of a photon. Example: [511.0, 1460.8, 2614.5]. The returned attenuation fractions will correspond 1:1 to these input energies."
  }
},
"required": ["Energies"]
})";
  
  registerTool({
    "photopeak_detection_efficiency",
    "Given an array of source photon energies, returns the fraction of photons that will contribute to a peak in data. The returned answer will include the attenuation factor of shielding (if specified), the fraction of photons making it to the detector (if distance is specified), the detection probability of photons that are incident upon the detector (i.e., the intrinsic detection efficiency) if a detector efficiency function is loaded (use 'detector_efficiency_function_info' tool call with no arguments to see if a detector efficiency function is loaded, and 'avaiable_detector_efficiency_functions' and 'load_detector_efficiency_function' to load an efficiency function), and gives total detection probability.",
    json::parse(par_schema),
    [](const json& params, InterSpec* interspec) -> json {
      return executePhotopeakDetectionCalc(params, interspec);
    }
  });

  registerTool({
    "get_materials",
    "Returns a list of avaiable materials for shielding. The materials name is what you will specify to other functions to use it as a shielding.",
    json{},
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetMaterials(interspec);
    }
  });

  registerTool({
    "get_material_info",
    "Returns shielding material information (density, effective atomic number, elemental mass fractions, etc).",
  json::parse(R"({
      "type": "object",
      "properties": {
        "material": { 
          "type": "string", 
          "description": "The shielding material name to get information about." 
          }
      },
      "required": ["material"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetMaterialInfo( params, interspec);
    }
  });


  registerTool({
    "avaiable_detector_efficiency_functions",
    "Returns a list of detector efficiency function names that can be loaded to the foreground spectrum.",
    json::parse(R"({
      "type": "object",
      "properties": {},
      "required": []
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeAvailableDetectors( params, interspec);
    }
  });
  
  registerTool({
    "load_detector_efficiency_function",
    "Loads a detector efficiency function to use for calculations. Can load from default available detectors, user previously used, filesystem path, GADRAS directory, or URI.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "identifier": {
          "type": "string",
          "description": "The detector efficiency function identifier. Can be a name from avaiable_detector_efficiency_functions, a filesystem path to a detector file or GADRAS directory, or a URI. Required."
        },
        "detectorName": {
          "type": "string",
          "description": "Optional: The specific detector name to use when the identifier points to a file containing multiple detector efficiencies (e.g., RelEff CSV files). If not specified and the file contains only one detector, that detector will be used. If not specified and the file contains multiple detectors, an error will be returned."
        },
        "source": {
          "type": "string",
          "enum": ["DefaultAvailable", "UserPreviouslyUsed", "FilePath", "GadrasDirectory", "URI", "AnySource"],
          "description": "Optional: Source type hint for loading. DefaultAvailable: built-in detectors. UserPreviouslyUsed: user's database. FilePath: filesystem path to detector file. GadrasDirectory: GADRAS detector directory. URI: web URI. AnySource: try all options in order. Defaults to AnySource if not specified."
        }
      },
      "required": ["identifier"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeLoadDetectorEfficiency( params, interspec);
    }
  });
  
  registerTool({
    "detector_efficiency_function_info",
    "Returns information (name, description, if has FWHM info, or if is fixed geometry, etc) about either the currently loaded detector efficiency function, or if a name is specified, that detectors efficiency function.",
    json::parse(R"({
      "type": "object",
      "properties": {
        "name": {
          "type": "string",
          "description": "Optional: The name of the detector efficiency function to return information about. If not specified, will return information about the currently loaded detector efficiency function."
          }
      },
      "required": []
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeGetDetectorInfo( params, interspec);
    }
  });

  registerTool({
    "search_sources_by_energy",
    "Searches for radiation sources (nuclides, x-rays, reactions) that have gammas, x-rays, alphas, and/or beta-end-point at specified energies. Returns candidate sources ranked by how well they match the input energies and current foreground spectrum.\n\nIt is best to only search on a single energy at a time, and only for foreground peaks whose source is not identified; do not search on multiple energies unless you are investigating or suspect they are from the same source, or are really having a hard time identifying a source.\n\nThe 'profile_score' metric (larger = better match) is the best indicator of fit quality for the loaded spectrum, considering peak heights and detector response.\n\nFurther investigation of top returned candidates is required to make sure they match the spectrum (e.g., the spectrum is not missing a peak it should have from that source, or the peaks in the spectrum have expected relative amplitudes, etc.)",
    json::parse(R"({
      "type": "object",
      "properties": {
        "energies": {
          "type": "array",
          "description": "The energies the source MUST match.  Since the source must match ALL entries, it is best to only search on one to a few unidentified peaks at a time, ususally starting with only a single peak, and then moving to two only if no satisfactory answer can be found. You usually will not want to search on a foreground peak that roughly aligns in energy with a background peak (either auto-search peak, or user/analysis fit peak), unless the foreground peak CPS is substantially elevated relative to the background peak.",
          "items": {
            "type": "object",
            "properties": {
              "energy": {
                "type": "number",
                "description": "Energy in keV to search for"
              },
              "window": {
                "type": "number",
                "description": "Search window (±keV) around this energy. If not specified, defaults to the recomended window of 1.27 times the expected FWHM at that energy, or if expected FWHM cant be determined, 1.5 keV for HPGe detectors or 10 keV for other detector types.  If the expected FWHM is not known, then it is recomended to not specify this field."
              }
            },
            "required": ["energy"]
          },
          "minItems": 1,
          "description": "Array of energies to search for. Each source must have emissions matching ALL specified energies."
        },
        "source_category": {
          "type": "string",
          "description": "Category filter for source types to search. Can use either the display name (e.g., 'Nuclides + X-rays') or the i18n key (e.g., 'isbe-category-nuc-xray'). Default is 'Nuclides + X-rays'. Common values: 'Nuclides + X-rays', 'Nuclides, X-rays, Reactions', 'Fluorescence x-rays', 'Reactions', 'Medical Nuclides', 'Industrial Nuclides', 'SNM', 'NORM', 'Field Encountered Sources', 'Field Encountered Nuclides', 'Fission Products', 'Alphas', 'Beta Endpoint', 'Un-aged Nuclides'."
        },
        "min_half_life": {
          "type": "string",
          "description": "Minimum nuclide half-life as a string (e.g., '10 s', '5 min', '1 h', '2 y'). Default is '100 m'."
        },
        "min_branching_ratio": {
          "type": "number",
          "description": "Minimum relative branching ratio (0.0 to 1.0). This is the aged gamma line BR divided by max intensity line. Default is 0.0."
        },
        "max_results": {
          "type": "integer",
          "description": "Maximum number of results to return. Default is 10."
        },
        "sort_by": {
          "type": "string",
          "enum": ["ProfileScore", "SumEnergyDifference"],
          "description": "Sort criterion: 'ProfileScore' (larger=better, considers spectrum fit and peak heights) or 'SumEnergyDifference' (smaller=better, simple energy distance). Default is 'ProfileScore'."
        }
      },
      "required": ["energies"]
    })"),
    [](const json& params, InterSpec* interspec) -> json {
      return executeSearchSourcesByEnergy( params, interspec);
    }
  });

  /*
   TODO: tool calls to add
     - executeGetNuclidesWithCharacteristicsInEnergyRange - this should be _replaced_ with an energyrange search - like the GUI tool, and not just the characteristic gamma search
     - `primary_gammas_for_source` has been renamed from `primary_gammas_for_nuclide` to better reflect it handles sources (nuclides, x-rays, reactions)
     - Activity/Shielding fit
     - RelActAuto
     - Loading of spectra from filesystem to foreground/background
     - Energy calibration...
     - Dose calc
   
   */
  
  

  /*
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
  */
  
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
  
  try
  {
#if( !defined(NDEBUG) && !BUILD_AS_UNIT_TEST_SUITE )
    cout << "Executing tool: " << toolName << " with params: " << parameters.dump() << endl;
#endif
    
    json result = tool->executor(parameters, interspec);
    
#if( !defined(NDEBUG) && !BUILD_AS_UNIT_TEST_SUITE )
    cout << "Tool result: " << result.dump() << endl;
#endif
    
    return result;
  }catch( const std::exception &e )
  {
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
  
nlohmann::json ToolRegistry::executeGetCharacteristicGammasForSource( const nlohmann::json& params )
{
  const string source = params.at("source").get<string>();
  json result;
  result["source"] = source;
  result["characteristicGammas"] = AnalystChecks::get_characteristic_gammas( source );
  return result;
}

nlohmann::json ToolRegistry::executeGetLoadedSpectra( const nlohmann::json& params, InterSpec* interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  vector<string> loadedSpectra;
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground ) )
    loadedSpectra.push_back("Foreground");
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::Background ) )
    loadedSpectra.push_back("Background");
  if( interspec->displayedHistogram( SpecUtils::SpectrumType::SecondForeground ) )
    loadedSpectra.push_back("Secondary");
  
  return json(loadedSpectra);
}
  
nlohmann::json ToolRegistry::executeGetNuclidesWithCharacteristicsInEnergyRange( const nlohmann::json& params, InterSpec* interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  vector<variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>> result;
  if( params.contains("lowerEnergy") && params.contains("upperEnergy") )
  {
    const double lower_energy = get_number( params, "lowerEnergy" );
    const double upper_energy = get_number( params, "upperEnergy" );
    result = AnalystChecks::get_nuclides_with_characteristics_in_energy_range( lower_energy, upper_energy, interspec );
  }else if( params.contains("energy") )
  {
    const double energy = get_number( params, "energy" );
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
  

nlohmann::json ToolRegistry::executeGetAssociatedSources( const nlohmann::json& params )
{
  const string nuclide = params.at("source").get<string>();

  const shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> info_db = MoreNuclideInfo::MoreNucInfoDb::instance();
  if( !info_db )
    throw std::runtime_error("No MoreNucInfoDb instance available");
    
  const MoreNuclideInfo::NucInfo *nuc_info = info_db->info(nuclide);
  if( !nuc_info )
    throw runtime_error( "No info for " + nuclide );
    
  vector<string> associated = nuc_info->m_associated;
    
  // We'll normalize assicated nuclides (e.g. return "Co60" instead of "Co-60"
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Could not open sandia.decay.xml" );
    
  for( string &val : associated )
  {
    const SandiaDecay::Nuclide *nuc = db->nuclide( val );
    if( nuc )
      val = nuc->symbol;
  }
    
  return json{{"source", nuclide}, {"associatedSources", associated} };
}//nlohmann::json executeGetAssociatedSources(const nlohmann::json& params, InterSpec* interspec)


nlohmann::json ToolRegistry::executeGetSourceAnalystNotes( const nlohmann::json& params )
{
  const string nuclide = params.at("source").get<string>();
  
  const shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> info_db = MoreNuclideInfo::MoreNucInfoDb::instance();
  if( !info_db )
    throw std::runtime_error("No MoreNucInfoDb instance available");
    
  const MoreNuclideInfo::NucInfo *nuc_info = info_db->info(nuclide);
  if( !nuc_info || nuc_info->m_notes.empty() )
    throw runtime_error( "No info for " + nuclide );
    
  return json{{"source", nuclide}, {"analystNotes", nuc_info->m_notes}};
}//nlohmann::json executeGetSourceAnalystNotes(const nlohmann::json& params, InterSpec* interspec)

  
nlohmann::json ToolRegistry::executeGetSourceInfo(const nlohmann::json& params, InterSpec* interspec )
{
  const string nuclide = params.at("source").get<string>();
  
  nlohmann::json result;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Could not initialize nuclide DecayDataBase." );
  
  const SandiaDecay::Element * el = nullptr;
  const ReactionGamma::Reaction *rctn = nullptr;
  const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
  if( nuc )
  {
    result["type"] = "nuclide";
    result["symbol"] = nuc->symbol;
    result["source"] = nuc->symbol;
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
  }//if( nuc )
  
  
  if( !nuc )
  {
    el = db->element( nuclide );
    if( el )
    {
      result["type"] = "x-ray";
      result["symbol"] = el->symbol;
      result["source"] = el->symbol;
      result["name"] = el->name;
      result["atomicNumber"] = static_cast<int>(el->atomicNumber);
      result["atomicMass"] = el->atomicMass();

      // We wont add x-rays so this way it keeps things consistent for the LLM to have to call "source_photons" for this info
      //nlohmann::json xraysArray = nlohmann::json::array();
      //for( const SandiaDecay::EnergyIntensityPair &xray : el->xrays )
      //{
      //  nlohmann::json xrayObj;
      //  xrayObj["energy"] = xray.energy;
      //  xrayObj["intensity"] = xray.intensity;
      //  xraysArray.push_back( xrayObj );
      //}
      //result["xrays"] = xraysArray;

      // Add naturally occurring isotopes array
      nlohmann::json isotopesArray = nlohmann::json::array();
      for( const SandiaDecay::NuclideAbundancePair &iso : el->isotopes )
      {
        nlohmann::json isoObj;
        isoObj["nuclide"] = iso.nuclide ? iso.nuclide->symbol : "";
        isoObj["abundance"] = iso.abundance;
        isotopesArray.push_back( isoObj );
      }
      result["naturallyOccuringIsotopes"] = isotopesArray;
    }
  }
  
  if( !nuc && !el )
  {
    const ReactionGamma * const rctn_db = ReactionGammaServer::database();
    
    std::vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
    if( rctn_db )
      rctn_db->gammas( nuclide, possible_rctns );
    
    if( !possible_rctns.empty() && possible_rctns[0].reaction ) // Take the first reaction found
      rctn = possible_rctns[0].reaction;
    
    if( rctn )
    {
      result["type"] = "reaction";
      result["name"] = rctn->name();
      result["source"] = rctn->name();
      if( !rctn->remark.empty() )
        result["remark"] = rctn->remark;
      
      if( rctn->targetNuclide )
        result["targetNuclide"] = rctn->targetNuclide->symbol;
      if( rctn->targetElement )
        result["targetElement"] = rctn->targetElement->symbol;
      if( rctn->productNuclide )
        result["productNuclide"] = rctn->productNuclide->symbol;
      
      switch( rctn->type )
      {
        case ReactionType::AlphaNeutron:    //gammas produced by alpha,n sources; X(a,n).
          result["process"] = "alpha-neutron";
          result["processDescription"] = "nucleus + alpha -> nucleus_{+2p+1n} + n + gammas";
          break;
        case ReactionType::NeutronAlpha:    //gammas produced by n,alpha sources; X(n,a).  ex. B10 + neutron -> Li7 + alpha + gammas
          result["process"] = "neutron-alpha";
          result["processDescription"] = "nucleus + neutron -> nucleus_{-2p-1n} + alpha + gammas";
          break;
        case ReactionType::AlphaProton:     //gammas produced by alpha,p sources; X(a,p).  ex. N14 + alpha -> O17 + proton + gammas
          result["process"] = "alpha-proton";
          result["processDescription"] = "nucleus + alpha -> nucleus_{+1p+2n} + proton + gammas";
          break;
        case ReactionType::NeutronCapture:
          result["process"] = "neutron capture";
          result["processDescription"] = "nucleus + neutron -> nucleus_{+1n} + gammas";
          break;
        case ReactionType::NeutronInelasticScatter:
          result["process"] = "neutron inelastic scatter";
          result["processDescription"] = "nucleus + neutron -> nucleus + neutron + gammas";
          break;
        case ReactionType::AlphaInelasticScatter:
          result["process"] = "alpha inelastic scatter";
          result["processDescription"] = "nucleus + alpha -> nucleus + alpha + gammas";
          break;
        case ReactionType::AnnihilationReaction:
          result["process"] = "Annihilation";
          result["processDescription"] = "gamma -> two 511 keV gammas";
          break;
        case ReactionType::NumReactionType:
          break;
      };//switch( rctn->type )
    }//if( rctn )
  }//if( !nuc && !el )
  
  
  if( !nuc && !el && !rctn )
    throw runtime_error( "Source '" + nuclide + "' is not a valid nuclide, x-ray element, or reaction." );
  
  try
  {
    const nlohmann::json associated = ToolRegistry::executeGetAssociatedSources( params );
    if( associated.contains("associatedSources") )
      result["associatedSources"] = associated["associatedSources"];
  }catch( std::exception & )
  {
  }
  
  try
  {
    const nlohmann::json analystNotes = ToolRegistry::executeGetSourceAnalystNotes( params );
    if( analystNotes.contains("analystNotes") )
      result["analystNotes"] = analystNotes["analystNotes"];
  }catch( std::exception & )
  {
  }

  try
  {
    nlohmann::json sourceCatagories;
    
    if( nuc )
      sourceCatagories = source_categories( nuc, interspec );
    else if( el )
      sourceCatagories = source_categories( el, interspec );
    else if( rctn )
      sourceCatagories = source_categories( rctn, interspec );
    
    if( !sourceCatagories.empty() )
    {
      assert( sourceCatagories.is_array() );
      result["sourceCatagories"] = sourceCatagories;
    }
  }catch( std::exception & )
  {
  }


  return result;
}//nlohmann::json executeGetSourceInfo(const nlohmann::json& params, InterSpec* interspec )

nlohmann::json ToolRegistry::executeGetNuclideDecayChain(const nlohmann::json& params )
{
  const string nuclide = params.at("nuclide").get<string>();

  nlohmann::json result;
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw std::runtime_error( "Could not initialize nuclide DecayDataBase." );
  
  const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide );
  if( !nuc )
    throw std::runtime_error( "Nuclide '" + nuclide + "' is not a valid nuclide." );
  
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
  

nlohmann::json ToolRegistry::executeGetAutomatedRiidId(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw std::runtime_error("Could not initialize nuclide DecayDataBase.");
  

  std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !meas )
    throw std::runtime_error("No foreground spectrum loaded.");
  
  nlohmann::json result;

  shared_ptr<const SpecUtils::DetectorAnalysis> riid_ana = meas->detectors_analysis();
  
  if( riid_ana )
  {
    vector<pair<string, string>> riid_nucs;  //<description, nuclide>
  
    for( const SpecUtils::DetectorAnalysisResult &res : riid_ana->results_ )
    {
      if( res.nuclide_.empty() || res.isEmpty() )
        continue;
  
      auto pos = std::find_if(riid_nucs.begin(), riid_nucs.end(),
          [&res](const auto& v) { return v.first == res.nuclide_; });
      if( pos != riid_nucs.end() )
        continue;
  
      string nuc_name = res.nuclide_;
  
      const SandiaDecay::Nuclide *nuc = db->nuclide(nuc_name);
      if( !nuc )
      {
        vector<string> fields;
        SpecUtils::split(fields, nuc_name, " \t,");
        for( const auto &v : fields )
        {
          nuc = db->nuclide(v);
          if( nuc )
          {
            nuc_name = v;
            break;
          }
        }//for( const auto &v : fields )
      }//if( !nuc )
  
      riid_nucs.push_back( {res.nuclide_, nuc ? nuc->symbol : ""} );
    }//for( loop over RIID results )

    if( riid_nucs.empty() )
    {
      result["detectorSystemId"]["description"] = "No nuclides found.";
    }else
    {
      result["detectorSystemId"]["description"] = riid_ana? riid_ana->algorithm_name_ : "On-board RIID algorithm";
      for( const auto &v : riid_nucs )
        result["detectorSystemId"]["nuclides"].push_back( v.second.empty() ? v.first : v.second );
    }
  }else
  {
    result["detectorSystemId"]["description"] = "No ID algorithm present.";
  }//if( riid_ana )


  const ReferencePhotopeakDisplay * const refWidget = interspec->referenceLinesWidget();
  if( refWidget )
  {
    shared_ptr<const ExternalRidResults> riid = refWidget->currentExternalRidResults();
    
    if( riid && !riid->isotopes.empty() )
    {
      result["externalRiidTool"]["description"] = riid->algorithmName;
      for( const ExternalRidIsotope &v : riid->isotopes )
        result["externalRiidTool"]["nuclides"].push_back( v.name );
    }
  }//if( refWidget )

  return result;
}//nlohmann::json executeGetAutomatedRiidId(const nlohmann::json& params, InterSpec* interspec)

nlohmann::json ToolRegistry::executeFitPeaksForNuclide(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  // Parse parameters into FitPeaksForNuclideOptions
  AnalystChecks::FitPeaksForNuclideOptions options;
  from_json(params, options);
  
  // Call the AnalystChecks function to perform the actual peak fitting
  const AnalystChecks::FitPeaksForNuclideStatus result = AnalystChecks::fit_peaks_for_nuclides( options, interspec );
  
  shared_ptr<const SpecUtils::Measurement> meas;
  if( interspec )
    meas = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  
  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result, meas );
  
  return result_json;
}//nlohmann::json executeFitPeaksForNuclide(const nlohmann::json& params, InterSpec* interspec)

nlohmann::json ToolRegistry::executeGetCountsInEnergyRange(const nlohmann::json& params, InterSpec* interspec)
{
  if (!interspec) {
    throw std::runtime_error("No InterSpec session available.");
  }


  const double lowerEnergy = get_number( params, "lowerEnergy" );
  const double upperEnergy = get_number( params, "upperEnergy" );
  if( upperEnergy < lowerEnergy )
    throw runtime_error( "lowerEnergy is larger than upperEnergy" );

  // Call the AnalystChecks function to get the counts in energy range
  const AnalystChecks::SpectrumCountsInEnergyRange result = AnalystChecks::get_counts_in_energy_range(lowerEnergy, upperEnergy, interspec);

  // Convert the result to JSON and return
  json result_json;
  to_json(result_json, result);

  return result_json;
}//nlohmann::json executeGetCountsInEnergyRange(const nlohmann::json& params, InterSpec* interspec)

nlohmann::json ToolRegistry::executeGetExpectedFwhm(const nlohmann::json& params, InterSpec* interspec)
{
  if (!interspec)
    throw std::runtime_error("No InterSpec session available.");

  const double energy = get_number( params, "energy" );

  const float fwhm = AnalystChecks::get_expected_fwhm( energy, interspec );
  
  json result;
  result["energy"] = energy;
  result["fwhm"] = fwhm;
  return result;
}//nlohmann::json executeGetExpectedFwhm(const nlohmann::json& params, InterSpec* interspec)

nlohmann::json ToolRegistry::executeCurrieMdaCalc(const nlohmann::json& params, InterSpec* interspec)
{
  if (!interspec)
    throw std::runtime_error("No InterSpec session available.");

  shared_ptr<const SpecUtils::Measurement> spectrum = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  if (!spectrum)
    throw std::runtime_error("No foreground spectrum loaded");

  // Parse only the selected parameters
  const double energy = get_number( params, "energy" );
  const double detection_probability = params.value("detectionProbability", 0.95);
  const float additional_uncertainty = params.value("additionalUncertainty", 0.0f);

  // Get expected FWHM for the energy to determine ROI width
  float fwhm = -1.0;
  shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
  if (meas && meas->detector() && meas->detector()->hasResolutionInfo())
    fwhm = meas->detector()->peakResolutionFWHM(static_cast<float>(energy));
  
  // Fallback FWHM estimation if needed
  if (fwhm <= 0.0) {
    const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
    const vector<float> pars = isHPGe ? vector<float>{1.54f, 0.264f, 0.33f} : vector<float>{-6.5f, 7.5f, 0.55f};
    fwhm = DetectorPeakResponse::peakResolutionFWHM(energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars);
  }

  if (fwhm <= 0.0)
    throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy));

  // Set up CurrieMdaInput with fixed values for ROI and side channels
  DetectionLimitCalc::CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = static_cast<float>(energy);
  
  // Set ROI to be ±1.5 FWHM around the energy
  const float roi_half_width = 1.25f * fwhm; // recommended by ISO 11929:2010, could instead use 1.19
  input.roi_lower_energy = static_cast<float>(energy) - roi_half_width;
  input.roi_upper_energy = static_cast<float>(energy) + roi_half_width;
  
  // Use fixed values for side channels (typical values)
  input.num_lower_side_channels = 4;
  input.num_upper_side_channels = 4;
  
  // Use the parsed parameters
  input.detection_probability = detection_probability;
  input.additional_uncertainty = additional_uncertainty;

  // Call the DetectionLimitCalc function to perform the calculation
  const DetectionLimitCalc::CurrieMdaResult result = DetectionLimitCalc::currie_mda_calc(input);

  // Convert the result to JSON and return
  json result_json;
  to_json(result_json, result);

  return result_json;
}//nlohmann::json executeCurrieMdaCalc(const nlohmann::json& params, InterSpec* interspec)


/** Returns a list of energy/intensity pairs associated with a nuclide, elemental fluorescence x-ray, or reaction.
 The energies are expressed in keV, and intensities are given as photons per source Becquerel.
 */
nlohmann::json ToolRegistry::executeGetSourcePhotons(const nlohmann::json& params){

  if( !params.contains("Source") )
    throw runtime_error( "'Source' parameter must be specified." );

  // The nuclide, element, or reaction for which to retrieve decay product data.
  //   Examples: 'U238' (nuclide), 'Pb' (element), 'H(n,g)' (reaction).
  const string &src = params["Source"];

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Could not initialize nuclide DecayDataBase - can not continue." );

  const SandiaDecay::Nuclide *nuc = nullptr;
  const SandiaDecay::Element *el = nullptr;
  const ReactionGamma::Reaction *rctn = nullptr;

  nuc = db->nuclide(src);

  if( !nuc )
    el = db->element(src);

  if( !nuc && !el )
  {
    const ReactionGamma * const rctn_db = ReactionGammaServer::database();

    std::vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
    if( rctn_db )
      rctn_db->gammas( src, possible_rctns );

    if( !possible_rctns.empty() && possible_rctns[0].reaction ) // Take the first reaction found
      rctn = possible_rctns[0].reaction;
  }//if( !nuc && !el )

  if( !nuc && !el && !rctn )
    throw runtime_error( "Could not interpret '" + src + "' as a nuclide, x-ray, or reaction." );

  double age_in_seconds = 0.0;
  const bool has_age = params.contains("Age");
  if( !nuc && has_age )
    throw runtime_error( "You can only specify 'Age' for a nuclide source" );

  if( params.contains("Age") )
  {
    const string &age_str = params["Age"];

    try
    {
      age_in_seconds = PhysicalUnits::stringToTimeDuration(age_str) / PhysicalUnits::second;
      if( age_in_seconds < 0.0 )
        throw runtime_error( "Nuclide age ('" + age_str + "')must be larger than zero." );
    }catch( std::exception &e )
    {
      throw runtime_error( "Could not interpret '" + age_str + "' as a time duration for nuclide age." );
    }
  }else if( nuc )
  {
    age_in_seconds = PeakDef::defaultDecayTime( nuc, nullptr ) / PhysicalUnits::second;
  }

  const int max_results = params.value("MaxResults", 125);
  if( max_results < 1 )
    throw runtime_error( "'MaxResults' must be 1 or larger." );

  vector<pair<double,double>> result;
  if( nuc )
  {
    SandiaDecay::NuclideMixture mix;

    // We will add the nuclide as aged, so this way the gammas will be for 1 bq parent activity
    mix.addAgedNuclideByActivity( nuc, 1.0*SandiaDecay::Bq, age_in_seconds * SandiaDecay::second );
    const vector<SandiaDecay::EnergyRatePair> energy_rate = mix.photons(0.0);
    result.reserve( energy_rate.size() );
    for( const SandiaDecay::EnergyRatePair &erp : energy_rate )
      result.emplace_back( erp.energy, erp.numPerSecond );
  }else if( el )
  {
    result.reserve( el->xrays.size() );
    for( const SandiaDecay::EnergyIntensityPair &eip : el->xrays )
      result.emplace_back( eip.energy, eip.intensity );
  }else
  {
    assert( rctn );
    result.reserve( rctn->gammas.size() );
    double max_abundance = 0.0;
    for( const ReactionGamma::Reaction::EnergyYield &g : rctn->gammas )
    {
      const double abundance = static_cast<double>(g.abundance);
      max_abundance = std::max( max_abundance, abundance );
      result.emplace_back( static_cast<double>(g.energy), abundance );
    }
    if( (max_abundance != 1.0) && (max_abundance > 0.0) )
    {
      for( pair<double, double> &energy_rate : result )
        energy_rate.second /= max_abundance;
    }
  }

  if( static_cast<int>(result.size()) > max_results )
  {
    std::sort( begin(result), end(result), []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
      return lhs.second > rhs.second;
    } );
    result.resize( static_cast<size_t>(max_results) );
  }

  std::sort( begin(result), end(result), []( const pair<double,double> &lhs, const pair<double,double> &rhs ){
    return lhs.first < rhs.first;
  } );


  nlohmann::json json_array;

  for( const pair<double,double> &val : result )
    json_array.push_back({val.first, val.second});

  return json_array;
}//nlohmann::json executeGetSourcePhotons(const nlohmann::json& params, InterSpec* interspec)



nlohmann::json ToolRegistry::executeGetAttenuationOfShielding( const nlohmann::json& params, InterSpec* interspec )
{
  using namespace PhysicalUnits;

  // Parse the energies array
  if( !params.contains("Energies") || !params["Energies"].is_array() )
    throw runtime_error( "Energies parameter must be specified as an array." );

  const json& energies_json = params["Energies"];
  vector<double> energies;
  for( const auto& energy_val : energies_json )
  {
    if( energy_val.is_number() )
      energies.push_back( energy_val.get<double>() * keV );
    else
      throw runtime_error( "Each energy must be a number in keV." );
  }

  if( energies.empty() )
    throw runtime_error( "At least one energy must be specified." );

  // Parse the shielding specification
  if( !params.contains("Shielding") || !params["Shielding"].is_object() )
    throw runtime_error( "Shielding parameter must be specified as an object." );

  const json& shielding = params["Shielding"];

  // Determine which shielding format is being used
  const bool has_ad = shielding.contains("AD");
  const bool has_an = shielding.contains("AN");
  const bool has_material = shielding.contains("Material");
  const bool has_thickness = shielding.contains("Thickness");

  // Validate that only one format is specified
  const bool is_generic = has_ad && has_an;
  const bool is_material = has_material && has_thickness;

  if( is_generic && is_material )
    throw runtime_error( "Cannot specify both {AD, AN} and {Material, Thickness}. Use only one format." );

  if( !is_generic && !is_material )
    throw runtime_error( "Shielding must be specified either as {AD, AN} or {Material, Thickness}." );

  // Check for extra fields
  if( is_generic && (has_material || has_thickness) )
    throw runtime_error( "When using {AD, AN} format, do not specify Material or Thickness." );

  if( is_material && (has_ad || has_an) )
    throw runtime_error( "When using {Material, Thickness} format, do not specify AD or AN." );

  json result = json::array();

  // Format 1: Areal density and effective atomic number
  if( is_generic )
  {
    const double areal_density = get_number( shielding, "AD" ) * (g / cm2);
    const double atomic_number = get_number( shielding, "AN" );

    if( atomic_number < 1.0 || atomic_number > 100.0 )
      throw runtime_error( "Atomic number must be between 1 and 100." );

    if( areal_density < 0.0 )
      throw runtime_error( "Areal density must be positive." );

    for( const double energy : energies )
    {
      const double mu = GammaInteractionCalc::transmition_coefficient_generic(
        static_cast<float>(atomic_number),
        static_cast<float>(areal_density),
        static_cast<float>(energy)
      );
      const double transmission_fraction = std::exp( -mu );
      result.push_back( transmission_fraction );
    }
  }
  // Format 2: Material and thickness
  else if( is_material )
  {
    if( !interspec )
      throw runtime_error( "InterSpec instance is required for material-based shielding." );

    MaterialDB *db = interspec->materialDataBase();
    if( !db )
      throw runtime_error( "Material database not available." );

    const string material_name = shielding["Material"].get<string>();
    const Material *material = db->material( material_name );

    if( !material )
      throw runtime_error( "Material '" + material_name + "' not found." );

    const string thickness_str = shielding["Thickness"].get<string>();
    const double thickness = PhysicalUnits::stringToDistance( thickness_str );

    if( thickness < 0.0 )
      throw runtime_error( "Thickness must be positive." );

    for( const double energy : energies )
    {
      const double mu = GammaInteractionCalc::transmition_coefficient_material(
        material,
        static_cast<float>(energy),
        static_cast<float>(thickness)
      );
      const double transmission_fraction = std::exp( -mu );
      result.push_back( transmission_fraction );
    }
  }

  return result;
}


nlohmann::json ToolRegistry::executeGetMaterials( InterSpec* interspec )
{
  using namespace PhysicalUnits;

  if( !interspec )
    throw std::runtime_error( "InterSpec instance is required for retrieving materials." );

  // To keep from overwhelming the context, we will just return a subset of possbile shieldings
  return nlohmann::json{ "acetone", "adipose tissue", "air", "aluminum oxide",
    "amber", "ammonia", "bakelite", "baratol high explosive",
    "benzene", "beryllium oxide", "bgo", "blood",
    "bone", "boracitol high explosive", "borax",  "brass",
    "bronze", "butane", "calcium_carbonate", "carbon_dioxide",
    "cellulose_cellophane", "cellulose_nitrate", "celotex", "cesium_iodide",
    "COMPB high explosive", "concrete", "dacron", "dry soil",
    "electronic soup", "ethane", "ethyl_alcohol", "explosive/simulant",
    "fiberglass", "glass_lead", "glass_plate", "glucose",
    "glycerol", "granite", "graphite", "gypsum",
    "Jet fuel", "kapton", "kevlar",  "lipoly",
    "mock he", "Monocalcium phosphate", "muscle", "mylar",
    "nylon 6,6", "nylon-8062", "nylon-6-10", "paraffin", "PBX-9404 high explosive", "PBX-9501 high explosive",
    "pine wood", "plastic_sc_vinyltoluene", "plexiglass", "Plutonium dioxide",
    "plywood", "polyacrylonitrile", "polycarbonate",
    "polychlorostyrene", "polyethylene", "polypropylene", "polystyrene", 
    "Polyurethane foam", "propane", "pvc - polyvinyl chloride", "pvt - polyvinyl toluene",
    "pyrex_glass", "rubber_butyl", "rubber_natural", "rubber_neoprene",
    "Salt water", "silver_iodide", "skin", "sodium_iodide", 
    "soft tissue", "soil", "stainless steel ss-304", "stainless-steel nist",
    "Steel - AISI 1040 Med. Carbon", "teflon",
    "Thorium oxide", "titanium_dioxide", "Uranium hexafluoride",
    "Uranium metal", "void", "water","wet soil"
  };

  /*  
  MaterialDB * const db = interspec->materialDataBase();
  if( !db )
    throw std::runtime_error( "Material database not available." );

  const std::vector<const Material *> materials = db->materials();

  json result = json::array();
  for( const Material * const material : materials )
  {
    json material_json;
    material_json["name"] = material->name;
    material_json["density"] = material->density / (g / cm3);
    material_json["effAN"] = material->massWeightedAtomicNumber();

    //if( !material->description.empty() )
    //  material_json["description"] = material->description;

    result.push_back( material_json );
  }

  return result;
  */
}


nlohmann::json ToolRegistry::executeGetMaterialInfo( const nlohmann::json& params, InterSpec* interspec )
{
  using namespace PhysicalUnits;

  if( !interspec )
    throw std::runtime_error( "InterSpec instance is required for retrieving material information." );

  MaterialDB * const db = interspec->materialDataBase();
  if( !db )
    throw std::runtime_error( "Material database not available." );

  const std::string material_name = params.at("material").get<std::string>();
  const Material * const material = db->material( material_name );

  if( !material )
    throw std::runtime_error( "Material '" + material_name + "' not found." );

  json result;
  result["name"] = material->name;
  result["density"] = material->density / (g / cm3);
  result["effectiveAtomicNumber"] = material->massWeightedAtomicNumber();

  if( !material->description.empty() )
    result["description"] = material->description;

  // Add chemical formula
  result["massFractionChemicalFormula"] = material->chemicalFormula();

  // Add element composition
  if( !material->elements.empty() )
  {
    json elements_json = json::array();
    for( const Material::ElementFractionPair &element_pair : material->elements )
    {
      json element_json;
      element_json["symbol"] = element_pair.first->symbol;
      element_json["name"] = element_pair.first->name;
      element_json["atomicNumber"] = static_cast<int>(element_pair.first->atomicNumber);
      element_json["massFraction"] = element_pair.second;
      elements_json.push_back( element_json );
    }
    result["elements"] = elements_json;
  }

  // Add nuclide composition if present
  if( !material->nuclides.empty() )
  {
    json nuclides_json = json::array();
    for( const Material::NuclideFractionPair &nuclide_pair : material->nuclides )
    {
      json nuclide_json;
      nuclide_json["symbol"] = nuclide_pair.first->symbol;
      nuclide_json["massFraction"] = nuclide_pair.second;
      nuclides_json.push_back( nuclide_json );
    }
    result["nuclides"] = nuclides_json;
  }

  return result;
}//nlohmann::json executeGetMaterialInfo( const nlohmann::json& params, InterSpec* interspec )
  
  
  /* Cooresponds to the "photopeak_detection_efficiency" callback, whose description is:
    "Given an array of source photon energies, returns the fraction of photons that will contribute to a peak in data. The returned answer will include the attenuation factor of shielding (if specified), the fraction of photons making it to the detector (if distance is specified), the detection probability of photons that are incident upon the detector (i.e., the intrinsic detection efficiency) if a detector efficiency function is loaded (use 'detector_efficiency_function_info' tool call with no arguments to see if a detector efficiency function is loaded, and 'avaiable_detector_efficiency_functions' and 'load_detector_efficiency_function' to load an efficiency function), and gives total detection probability.",

   And paramaters accepted is:
   ```
   {
      "type": "object",
      "properties": {
      "Shielding": {
      "type": "array",
      "description": "Optional: An array of shielding objects, applied in order from source to detector. Each shielding layer attenuates the photons that pass through it. Each shielding object must be specified with one of the following formats:\n\n1. **Areal density and effective atomic number**: Provide an object with the keys `AD` (areal density in g/cm²) and `AN` (effective atomic number). Example: `{ \"AD\": 20.25, \"AN\": 26 }`.\n\n2. **Material and thickness**: Provide an object with the keys `Material` (element symbol or name) and `Thickness` (thickness in cm). Example: `{ \"Material\": \"Fe\", \"Thickness\": \"1.25cm\" }`.\n\nAvailable materials include element symbols or names and materials returned by the `get_materials` tool.",
   "items": {
     "type": "object",
     "properties": {
       "AD": {
         "type": "number",
         "description": "Areal density of the shielding material in g/cm²."
       },
       "AN": {
         "type": "number",
         "description": "Effective atomic number of the shielding material."
       },
       "Material": {
         "type": "string",
         "description": "The shielding material, specified as an element symbol or name. Example: 'Fe' or 'Iron'."
       },
       "Thickness": {
         "type": "string",
         "description": "The thickness of the shielding material, specified as a string with units. Example: '1.25cm'."
       }
     },
     "additionalProperties": false
   }
    },
    "Distance": {
      "type": "string",
      "description": "Optional: The distance the detector is from the center of the radioactive source.  If not specified, the gemoetric detection factor will not be accounted for. If a detector efficiency is not currently loaded for the foreground, or the detector efficiency is for a 'fixed geometry', then distance may not be specified. Example distances: '1.25 cm', '3 ft, 2 inches', etc."
    },
    "IncludeAirAttenuation": {
      "type": "boolean",
      "description": "Optional: If true, includes attenuation from air between the source and detector. Only applies when Distance is specified. Defaults to false."
    },
    "Energies": {
      "type": "array",
      "items": {
        "type": "number"
      },
      "description": "An array of photon energy values in keV to apply the calculation to. Each value must be a positive float representing the energy of a photon. Example: [511.0, 1460.8, 2614.5]. The returned attenuation fractions will correspond 1:1 to these input energies."
    }
  },
  "required": ["Energies"]
  ```
*/
nlohmann::json ToolRegistry::executePhotopeakDetectionCalc(const nlohmann::json& params, InterSpec* interspec)
{
  using namespace std;
  using namespace nlohmann;

  // Step 1: Parse and normalize inputs

  // Parse Energies - allow single number or array
  vector<double> energies;
  if( params.contains("Energies") )
  {
    if( params["Energies"].is_number() )
    {
      energies.push_back( params["Energies"].get<double>() );
    }
    else if( params["Energies"].is_array() )
    {
      energies = params["Energies"].get<vector<double>>();
    }
    else
    {
      throw runtime_error( "Energies must be a number or array of numbers" );
    }
  }
  else
  {
    throw runtime_error( "Energies parameter is required" );
  }

  if( energies.empty() )
    throw runtime_error( "At least one energy must be specified" );

  for( const double energy : energies )
  {
    if( energy <= 0.0 )
      throw runtime_error( "All energies must be positive, got " + to_string(energy) + " keV" );
  }

  // Parse Shielding - allow single object or array
  vector<json> shieldings;
  if( params.contains("Shielding") )
  {
    if( params["Shielding"].is_object() )
    {
      shieldings.push_back( params["Shielding"] );
    }
    else if( params["Shielding"].is_array() )
    {
      shieldings = params["Shielding"].get<vector<json>>();
    }
    else
    {
      throw runtime_error( "Shielding must be an object or array of objects" );
    }
  }

  // Parse Distance (optional)
  double distance = 0.0;
  bool has_distance = false;
  if( params.contains("Distance") && params["Distance"].is_string() )
  {
    const string distance_str = params["Distance"].get<string>();
    distance = PhysicalUnits::stringToDistance( distance_str );
    has_distance = true;
  }

  // Parse IncludeAirAttenuation (optional, default false)
  const bool include_air = params.value("IncludeAirAttenuation", false);

  // Step 2: Validate parameters

  if( include_air && !has_distance )
    throw runtime_error( "IncludeAirAttenuation requires Distance to be specified" );

  // Get detector and material database
  shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  shared_ptr<DetectorPeakResponse> detector;
  if( meas )
    detector = meas->detector();

  if( has_distance )
  {
    if( !detector || !detector->isValid() )
      throw runtime_error( "Distance specified but no detector efficiency function is currently loaded" );

    if( detector->isFixedGeometry() )
      throw runtime_error( "Distance cannot be specified when detector efficiency function is for fixed geometry" );
  }

  MaterialDB *materialDB = interspec->materialDataBase();
  if( !materialDB )
    throw runtime_error( "Material database not available" );

  // Step 3: Parse shielding objects and calculate thicknesses

  struct ShieldingInfo
  {
    bool is_generic;  // true for AD/AN, false for Material/Thickness
    double atomic_number;
    double areal_density;
    const Material *material;
    double thickness;
  };

  vector<ShieldingInfo> parsed_shieldings;
  double total_shielding_thickness = 0.0;

  for( const json &shield_json : shieldings )
  {
    ShieldingInfo info;
    info.is_generic = false;
    info.atomic_number = 0.0;
    info.areal_density = 0.0;
    info.material = nullptr;
    info.thickness = 0.0;

    const bool has_ad = shield_json.contains("AD");
    const bool has_an = shield_json.contains("AN");
    const bool has_material = shield_json.contains("Material");
    const bool has_thickness = shield_json.contains("Thickness");

    if( has_ad && has_an && !has_material && !has_thickness )
    {
      // Generic shielding (AD/AN format)
      info.is_generic = true;
      info.areal_density = shield_json["AD"].get<double>();
      info.atomic_number = shield_json["AN"].get<double>();
      info.thickness = 0.0;  // Generic shielding has no physical thickness

      if( info.areal_density < 0.0 )
        throw runtime_error( "Areal density (AD) must be non-negative, got " + to_string(info.areal_density) );

      if( info.atomic_number <= 0.0 )
        throw runtime_error( "Atomic number (AN) must be positive, got " + to_string(info.atomic_number) );

      // Convert AD from g/cm2 to PhysicalUnits
      info.areal_density *= (PhysicalUnits::g / PhysicalUnits::cm2);
    }
    else if( has_material && has_thickness && !has_ad && !has_an )
    {
      // Material shielding
      info.is_generic = false;
      const string material_name = shield_json["Material"].get<string>();
      const string thickness_str = shield_json["Thickness"].get<string>();

      info.material = materialDB->material( material_name );
      if( !info.material )
        throw runtime_error( "Material '" + material_name + "' not found in material database" );

      info.thickness = PhysicalUnits::stringToDistance( thickness_str );

      if( info.thickness < 0.0 )
        throw runtime_error( "Thickness must be non-negative for material '" + material_name + "', got " + thickness_str );

      total_shielding_thickness += info.thickness;
    }
    else
    {
      throw runtime_error( "Each shielding must have either (AD and AN) or (Material and Thickness), not a mix" );
    }

    parsed_shieldings.push_back( info );
  }

  // Step 4: Calculate air distance
  double air_distance = 0.0;
  if( include_air )
  {
    air_distance = distance - total_shielding_thickness;
    if( air_distance < 0.0 )
      throw runtime_error( "Total shielding thickness (" + to_string(total_shielding_thickness/PhysicalUnits::cm)
                          + " cm) exceeds specified distance (" + to_string(distance/PhysicalUnits::cm) + " cm)" );
  }

  // Step 5: Calculate results for each energy
  json results = json::array();

  for( const double energy_kev : energies )
  {
    const float energy = static_cast<float>( energy_kev * PhysicalUnits::keV );

    json result;
    result["energy"] = energy_kev;

    double final_efficiency = 1.0;

    // Calculate shielding attenuations
    json shielding_attenuations = json::array();
    for( const ShieldingInfo &shield : parsed_shieldings )
    {
      double attenuation;

      if( shield.is_generic )
      {
        const double mu = GammaInteractionCalc::transmition_coefficient_generic(
          static_cast<float>(shield.atomic_number),
          static_cast<float>(shield.areal_density),
          energy
        );
        attenuation = exp( -mu );
      }
      else
      {
        const double mu = GammaInteractionCalc::transmition_coefficient_material(
          shield.material,
          energy,
          static_cast<float>(shield.thickness)
        );
        attenuation = exp( -mu );
      }

      shielding_attenuations.push_back( attenuation );
      final_efficiency *= attenuation;
    }

    if( !shielding_attenuations.empty() )
      result["shieldingAttenuations"] = shielding_attenuations;

    // Calculate air attenuation
    if( include_air )
    {
      const double mu = GammaInteractionCalc::transmission_coefficient_air( energy, static_cast<float>(air_distance) );
      const double air_attenuation = exp( -mu );
      result["airAttenuation"] = air_attenuation;
      final_efficiency *= air_attenuation;
    }

    // Calculate distance geometry factor
    if( has_distance && detector )
    {
      const double detector_diameter = detector->detectorDiameter();
      const double solid_angle = DetectorPeakResponse::fractionalSolidAngle( detector_diameter, distance );
      result["distanceGeometryFactor"] = solid_angle;
      final_efficiency *= solid_angle;
    }

    // Calculate detector intrinsic efficiency
    if( detector && detector->isValid() )
    {
      const float intrinsic_eff = detector->intrinsicEfficiency( energy );
      result["detectorIntrinsicEfficiency"] = intrinsic_eff;
      final_efficiency *= intrinsic_eff;
    }

    result["finalEfficiency"] = final_efficiency;
    results.push_back( result );
  }

  // Step 6: Build and return result
  json response;
  response["energies"] = energies;
  response["results"] = results;

  return response;
}//nlohmann::json executePhotopeakDetectionCalc(const nlohmann::json& params, InterSpec* interspec)
  

  /** cooresponds to the `avaiable_detector_efficiency_functions` tool call, which has the description:
   "Returns a list of detector efficiency function names, and the detector efficiency source, that can be loaded to the foreground spectrum."
 */
nlohmann::json ToolRegistry::executeAvailableDetectors(const nlohmann::json& params, InterSpec* interspec)
{
  json result = json::array();

  // Track detector names we've already added to avoid duplicates
  std::set<std::string> addedDetectors;

  // Get currently loaded detector to add at the end if not already in list
  std::shared_ptr<DetectorPeakResponse> currentDet;
  if( interspec )
  {
    std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
    if( meas )
      currentDet = meas->detector();
  }

  // Get available RelEff detector files
  try
  {
    const std::vector<std::string> relEffFiles = DrfSelect::potential_rel_eff_det_files();

    for( const std::string &filepath : relEffFiles )
    {
      try
      {
        std::vector<std::string> credits;
        std::vector<std::shared_ptr<DetectorPeakResponse>> drfs;

#ifdef _WIN32
        const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16(filepath);
        std::ifstream input( wfilepath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
        std::ifstream input( filepath.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

        if( !input.is_open() )
          continue;

        DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );

        for( const std::shared_ptr<DetectorPeakResponse> &drf : drfs )
        {
          if( !drf || !drf->isValid() || drf->name().empty() )
            continue;

          // Skip if we've already added this detector
          if( addedDetectors.count(drf->name()) )
            continue;

          json detectorInfo;
          detectorInfo["name"] = drf->name();
          detectorInfo["source"] = "DefaultAvailable";
          if( !drf->description().empty() )
            detectorInfo["description"] = drf->description();

          result.push_back( detectorInfo );
          addedDetectors.insert( drf->name() );
        }
      }catch( std::exception &e )
      {
        // Skip files that can't be parsed
        std::cerr << "Error parsing RelEff file " << filepath << ": " << e.what() << std::endl;
      }
    }
  }catch( std::exception &e )
  {
    std::cerr << "Error getting RelEff detector files: " << e.what() << std::endl;
  }

  // Get available GADRAS detector directories
  try
  {
    const std::vector<std::string> gadrasDirs = DrfSelect::potential_gadras_det_dirs( interspec );

    for( const std::string &dir : gadrasDirs )
    {
      try
      {
        const std::vector<std::string> drfDirs = DrfSelect::recursive_list_gadras_drfs( dir );

        for( const std::string &drfDir : drfDirs )
        {
          try
          {
            std::shared_ptr<DetectorPeakResponse> drf = DrfSelect::initAGadrasDetectorFromDirectory( drfDir );

            if( !drf || !drf->isValid() || drf->name().empty() )
              continue;

            // Skip if we've already added this detector
            if( addedDetectors.count(drf->name()) )
              continue;

            json detectorInfo;
            detectorInfo["name"] = drf->name();
            detectorInfo["source"] = "DefaultAvailable";
            if( !drf->description().empty() )
              detectorInfo["description"] = drf->description();

            result.push_back( detectorInfo );
            addedDetectors.insert( drf->name() );
          }catch( std::exception &e )
          {
            // Skip directories that can't be parsed
            std::cerr << "Error parsing GADRAS DRF from " << drfDir << ": " << e.what() << std::endl;
          }
        }
      }catch( std::exception &e )
      {
        std::cerr << "Error listing GADRAS DRFs in " << dir << ": " << e.what() << std::endl;
      }
    }
  }catch( std::exception &e )
  {
    std::cerr << "Error getting GADRAS detector directories: " << e.what() << std::endl;
  }

  // Add currently loaded detector at the end if it exists and hasn't been added yet
  if( currentDet && currentDet->isValid() && !currentDet->name().empty()
      && !addedDetectors.count(currentDet->name()) )
  {
    json detectorInfo;
    detectorInfo["name"] = currentDet->name();
    detectorInfo["source"] = "CurrentlyLoaded";
    if( !currentDet->description().empty() )
      detectorInfo["description"] = currentDet->description();

    result.push_back( detectorInfo );
    addedDetectors.insert( currentDet->name() );
  }

  // Add previously used detectors from database
  if( interspec )
  {
    try
    {
      std::shared_ptr<DataBaseUtils::DbSession> sql = interspec->sql();
      if( sql )
      {
        DataBaseUtils::DbTransaction transaction( *sql );

        Wt::Dbo::collection<Wt::Dbo::ptr<DetectorPeakResponse>> det_effs
          = sql->session()->find<DetectorPeakResponse>()
                           .where( "InterSpecUser_id = ?" ).bind( interspec->user().id() )
                           .orderBy( "-1*m_lastUsedUtc" );

        for( auto iter = det_effs.begin(); iter != det_effs.end(); ++iter )
        {
          Wt::Dbo::ptr<DetectorPeakResponse> det_ptr = *iter;
          if( !det_ptr )
            continue;

          const std::string detName = det_ptr->name();
          if( detName.empty() || addedDetectors.count(detName) )
            continue;

          json detectorInfo;
          detectorInfo["name"] = detName;
          detectorInfo["source"] = "UserPreviouslyUsed";

          const std::string detDesc = det_ptr->description();
          if( !detDesc.empty() )
            detectorInfo["description"] = detDesc;

          result.push_back( detectorInfo );
          addedDetectors.insert( detName );
        }

        transaction.commit();
      }
    }catch( std::exception &e )
    {
      std::cerr << "Error getting previously used detectors from database: " << e.what() << std::endl;
    }
  }

  return result;
}//nlohmann::json executeAvailableDetectors(const nlohmann::json& params, InterSpec* interspec)


std::shared_ptr<DetectorPeakResponse> ToolRegistry::findDetectorByIdentifier(
  const std::string& identifier,
  const std::string& detectorName,
  const std::string& sourceHint,
  InterSpec* interspec,
  std::string& loadedFrom
)
{
  std::shared_ptr<DetectorPeakResponse> drf;
  loadedFrom.clear();

  // Helper lambda to try loading from a file path
  auto tryLoadFromFilePath = [&](const std::string& path) -> std::shared_ptr<DetectorPeakResponse> {
    if( !SpecUtils::is_file(path) && !SpecUtils::is_directory(path) )
      return nullptr;

    // Try as GADRAS directory
    if( SpecUtils::is_directory(path) )
    {
      try
      {
        std::shared_ptr<DetectorPeakResponse> det = DrfSelect::initAGadrasDetectorFromDirectory( path );
        if( det && det->isValid() )
        {
          loadedFrom = "GadrasDirectory";
          return det;
        }
      }catch( std::exception & )
      {
      }
    }

    // Read file contents to check format
    try
    {
#ifdef _WIN32
      const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16(path);
      std::ifstream input( wpath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
      std::ifstream input( path.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

      if( input.is_open() )
      {
        std::string contents;
        input.seekg( 0, std::ios::end );
        contents.resize( input.tellg() );
        input.seekg( 0, std::ios::beg );
        input.read( &contents[0], contents.size() );
        input.close();

        // Check if contents look like a URI (starts with VER= or contains DRF parameters)
        if( contents.find("VER=") != std::string::npos
            || (contents.find("DIAM=") != std::string::npos && contents.find("EFFT=") != std::string::npos) )
        {
          try
          {
            std::shared_ptr<DetectorPeakResponse> det = std::make_shared<DetectorPeakResponse>();
            det->fromAppUrl( contents );
            if( det->isValid() )
            {
              loadedFrom = "FilePath";
              return det;
            }
          }catch( std::exception & )
          {
            // Not a valid URI format, continue to other formats
          }
        }

        // Check if it looks like XML
        if( contents.find("<?xml") != std::string::npos || contents.find("<DetectorPeakResponse") != std::string::npos )
        {
          std::shared_ptr<DetectorPeakResponse> det = std::make_shared<DetectorPeakResponse>();
          rapidxml::xml_document<char> doc;
          doc.parse<rapidxml::parse_trim_whitespace | rapidxml::parse_normalize_whitespace>( &contents[0] );
          rapidxml::xml_node<char> *node = doc.first_node();
          if( node )
          {
            det->fromXml( node );
            if( det->isValid() )
            {
              loadedFrom = "FilePath";
              return det;
            }
          }
        }
      }
    }catch( std::exception & )
    {
    }

    // Try as RelEff CSV file
    try
    {
      std::vector<std::string> credits;
      std::vector<std::shared_ptr<DetectorPeakResponse>> drfs;

#ifdef _WIN32
      const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16(path);
      std::ifstream input( wpath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
      std::ifstream input( path.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

      if( input.is_open() )
      {
        DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );

        if( !detectorName.empty() )
        {
          // Find specific detector by name
          for( const auto& det : drfs )
          {
            if( det && det->isValid() && det->name() == detectorName )
            {
              loadedFrom = "FilePath";
              return det;
            }
          }
        }
        else if( drfs.size() == 1 && drfs[0] && drfs[0]->isValid() )
        {
          // Only one detector in file
          loadedFrom = "FilePath";
          return drfs[0];
        }
        else if( drfs.size() > 1 )
        {
          throw std::runtime_error( "File contains multiple detectors. Please specify detectorName parameter." );
        }
      }
    }catch( std::exception & )
    {
    }

    // Try as ECC file
    try
    {
#ifdef _WIN32
      const std::wstring wpath = SpecUtils::convert_from_utf8_to_utf16(path);
      std::ifstream input( wpath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
      std::ifstream input( path.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

      if( input.is_open() )
      {
        auto result = DetectorPeakResponse::parseEccFile( input );
        std::shared_ptr<DetectorPeakResponse> det = std::get<0>( result );
        if( det && det->isValid() )
        {
          loadedFrom = "FilePath";
          return det;
        }
      }
    }catch( std::exception & )
    {
    }

    return nullptr;
  };

  // Try loading based on source hint
  if( sourceHint == "DefaultAvailable" || sourceHint == "AnySource" )
  {
    // Try from RelEff files
    try
    {
      const std::vector<std::string> relEffFiles = DrfSelect::potential_rel_eff_det_files();
      for( const std::string& filepath : relEffFiles )
      {
        try
        {
          std::vector<std::string> credits;
          std::vector<std::shared_ptr<DetectorPeakResponse>> drfs;

#ifdef _WIN32
          const std::wstring wfilepath = SpecUtils::convert_from_utf8_to_utf16(filepath);
          std::ifstream input( wfilepath.c_str(), std::ios_base::binary | std::ios_base::in );
#else
          std::ifstream input( filepath.c_str(), std::ios_base::binary | std::ios_base::in );
#endif

          if( input.is_open() )
          {
            DetectorPeakResponse::parseMultipleRelEffDrfCsv( input, credits, drfs );
            for( const auto& det : drfs )
            {
              if( det && det->isValid() && det->name() == identifier )
              {
                drf = det;
                loadedFrom = "DefaultAvailable";
                break;
              }
            }
          }
          if( drf )
            break;
        }catch( std::exception & )
        {
        }
      }
    }catch( std::exception & )
    {
    }

    // Try from GADRAS directories
    if( !drf && interspec )
    {
      try
      {
        const std::vector<std::string> gadrasDirs = DrfSelect::potential_gadras_det_dirs( interspec );
        for( const std::string& dir : gadrasDirs )
        {
          try
          {
            const std::vector<std::string> drfDirs = DrfSelect::recursive_list_gadras_drfs( dir );
            for( const std::string& drfDir : drfDirs )
            {
              try
              {
                std::shared_ptr<DetectorPeakResponse> det = DrfSelect::initAGadrasDetectorFromDirectory( drfDir );
                if( det && det->isValid() && det->name() == identifier )
                {
                  drf = det;
                  loadedFrom = "DefaultAvailable";
                  break;
                }
              }catch( std::exception & )
              {
              }
            }
            if( drf )
              break;
          }catch( std::exception & )
          {
          }
        }
      }catch( std::exception & )
      {
      }
    }
  }

  // Try from user previously used
  if( !drf && interspec && (sourceHint == "UserPreviouslyUsed" || sourceHint == "AnySource") )
  {
    try
    {
      std::shared_ptr<DataBaseUtils::DbSession> sql = interspec->sql();
      if( sql )
      {
        DataBaseUtils::DbTransaction transaction( *sql );
        Wt::Dbo::collection<Wt::Dbo::ptr<DetectorPeakResponse>> det_effs
          = sql->session()->find<DetectorPeakResponse>()
                           .where( "InterSpecUser_id = ?" ).bind( interspec->user().id() )
                           .where( "m_name = ?" ).bind( identifier );

        if( det_effs.size() > 0 )
        {
          Wt::Dbo::ptr<DetectorPeakResponse> det_ptr = *det_effs.begin();
          if( det_ptr )
          {
            drf = std::make_shared<DetectorPeakResponse>( *det_ptr );
            loadedFrom = "UserPreviouslyUsed";
          }
        }
        transaction.commit();
      }
    }catch( std::exception & )
    {
    }
  }

  // Try as file path
  if( !drf && (sourceHint == "FilePath" || sourceHint == "GadrasDirectory" || sourceHint == "AnySource") )
  {
    drf = tryLoadFromFilePath( identifier );
  }

  // Try as URI
  if( !drf && (sourceHint == "URI" || sourceHint == "AnySource") )
  {
    try
    {
      std::shared_ptr<DetectorPeakResponse> det = std::make_shared<DetectorPeakResponse>();
      det->fromAppUrl( identifier );
      if( det->isValid() )
      {
        drf = det;
        loadedFrom = "URI";
      }
    }catch( std::exception & )
    {
      // Not a valid URI or failed to parse
    }
  }

  return drf;
}//findDetectorByIdentifier()


  /** Corresponds to the `load_detector_efficiency_function` callback.
   Description "Loads a detector efficiency function to use for calculations. Can load from default available detectors, user previously used, filesystem path, GADRAS directory, or URI."
   Parameter description:
```
   {
      "type": "object",
      "properties": {
        "identifier": {
          "type": "string",
          "description": "The detector efficiency function identifier. Can be a name from avaiable_detector_efficiency_functions, a filesystem path to a detector file or GADRAS directory, or a URI. Required."
        },
        "detectorName": {
          "type": "string",
          "description": "Optional: The specific detector name to use when the identifier points to a file containing multiple detector efficiencies (e.g., RelEff CSV files). If not specified and the file contains only one detector, that detector will be used. If not specified and the file contains multiple detectors, an error will be returned."
        },
        "source": {
          "type": "string",
          "enum": ["DefaultAvailable", "UserPreviouslyUsed", "FilePath", "GadrasDirectory", "URI", "AnySource"],
          "description": "Optional: Source type hint for loading. DefaultAvailable: built-in detectors. UserPreviouslyUsed: user's database. FilePath: filesystem path to detector file. GadrasDirectory: GADRAS detector directory. URI: web URI. AnySource: try all options in order. Defaults to AnySource if not specified."
        }
      },
      "required": ["identifier"]
    }
   ```
 */
nlohmann::json ToolRegistry::executeLoadDetectorEfficiency(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Check if a foreground spectrum is loaded
  std::shared_ptr<SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    throw std::runtime_error( "Cannot load detector efficiency function: no foreground spectrum is currently loaded. Please load a spectrum file first." );

  if( !params.contains("identifier") || !params["identifier"].is_string() )
    throw std::runtime_error( "Missing required 'identifier' parameter." );

  const std::string identifier = params["identifier"].get<std::string>();
  const std::string detectorName = params.value("detectorName", std::string());
  const std::string sourceHint = params.value("source", std::string("AnySource"));

  std::string loadedFrom;
  std::shared_ptr<DetectorPeakResponse> drf = ToolRegistry::findDetectorByIdentifier( identifier, detectorName, sourceHint, interspec, loadedFrom );

  if( !drf || !drf->isValid() )
  {
    std::string errorMsg = "Failed to load detector efficiency function '" + identifier + "'";
    if( !detectorName.empty() )
      errorMsg += " with detector name '" + detectorName + "'";
    errorMsg += " from source hint: " + sourceHint;
    throw std::runtime_error( errorMsg );
  }

  // Load the detector
  interspec->detectorChanged().emit( drf );

  // Return success response
  json result;
  result["success"] = true;
  result["detectorName"] = drf->name();
  result["source"] = loadedFrom;
  if( !drf->description().empty() )
    result["description"] = drf->description();

  return result;
}//nlohmann::json executeLoadDetectorEfficiency(const nlohmann::json& params, InterSpec* interspec)
  

/** Cooresponds to the `detector_efficiency_function_info` tool call.
 Description: "Returns information (name, description, if has FWHM info, or if is fixed geometry, etc) about either the currently loaded detector efficiency function, or if a name is specified, that detectors efficiency function.",
 Parameter description:
 ```
 {
   "type": "object",
   "properties": {
     "name": {
       "type": "string",
       "description": "Optional: The name of the detector efficiency function to return information about. If not specified, will return information about the currently loaded detector efficiency function."
     }
   },
   "required": []
 }
 ```
*/
nlohmann::json ToolRegistry::executeGetDetectorInfo(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  std::shared_ptr<DetectorPeakResponse> drf;

  // Check if a specific detector name was requested
  if( params.contains("name") && params["name"].is_string() )
  {
    const std::string name = params["name"].get<std::string>();
    std::string loadedFrom;

    // Use helper function to find the detector
    drf = ToolRegistry::findDetectorByIdentifier( name, "", "AnySource", interspec, loadedFrom );

    if( !drf )
    {
      throw std::runtime_error( "Detector efficiency function '" + name + "' not found." );
    }
  }
  else
  {
    // Get currently loaded detector
    std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
    if( meas )
      drf = meas->detector();

    if( !drf )
    {
      throw std::runtime_error( "No detector efficiency function currently loaded." );
    }
  }

  // Build the response JSON
  json result;
  result["name"] = drf->name();
  result["description"] = drf->description();
  result["isValid"] = drf->isValid();
  result["hasResolutionInfo"] = drf->hasResolutionInfo();

  // Geometry type
  const DetectorPeakResponse::EffGeometryType geomType = drf->geometryType();
  std::string geomTypeStr;
  std::string geomTypeDesc;

  switch( geomType )
  {
    case DetectorPeakResponse::EffGeometryType::FarField:
      geomTypeStr = "FarField";
      geomTypeDesc = "Detection efficiency varies with ~1/r²";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
      geomTypeStr = "FixedGeomTotalAct";
      geomTypeDesc = "Fixed geometry, full-energy efficiency per source decay (total activity)";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
      geomTypeStr = "FixedGeomActPerCm2";
      geomTypeDesc = "Fixed geometry, efficiency per cm² surface area";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
      geomTypeStr = "FixedGeomActPerM2";
      geomTypeDesc = "Fixed geometry, efficiency per m² surface area";
      break;
    case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
      geomTypeStr = "FixedGeomActPerGram";
      geomTypeDesc = "Fixed geometry, efficiency per gram of source";
      break;
  }

  result["geometryType"] = geomTypeStr;
  result["geometryTypeDescription"] = geomTypeDesc;
  result["isFixedGeometry"] = (geomType != DetectorPeakResponse::EffGeometryType::FarField);

  // Energy range
  result["lowerEnergy"] = drf->lowerEnergy();
  result["upperEnergy"] = drf->upperEnergy();

  // Detector diameter (only for far field)
  if( geomType == DetectorPeakResponse::EffGeometryType::FarField )
  {
    result["detectorDiameter"] = drf->detectorDiameter();
  }

  // Hash and source
  result["hash"] = drf->hashValue();

  const DetectorPeakResponse::DrfSource drfSource = drf->drfSource();
  std::string drfSourceStr;
  switch( drfSource )
  {
    case DetectorPeakResponse::DrfSource::DefaultGadrasDrf:
      drfSourceStr = "DefaultGadrasDrf";
      break;
    case DetectorPeakResponse::DrfSource::UserAddedGadrasDrf:
      drfSourceStr = "UserAddedGadrasDrf";
      break;
    case DetectorPeakResponse::DrfSource::UserAddedRelativeEfficiencyDrf:
      drfSourceStr = "UserAddedRelativeEfficiencyDrf";
      break;
    case DetectorPeakResponse::DrfSource::UserSpecifiedFormulaDrf:
      drfSourceStr = "UserSpecifiedFormulaDrf";
      break;
    case DetectorPeakResponse::DrfSource::FromSpectrumFileDrf:
      drfSourceStr = "FromSpectrumFileDrf";
      break;
    case DetectorPeakResponse::DrfSource::UserCreatedDrf:
      drfSourceStr = "UserCreatedDrf";
      break;
    case DetectorPeakResponse::DrfSource::IsocsEcc:
      drfSourceStr = "IsocsEcc";
      break;
    case DetectorPeakResponse::DrfSource::UnknownDrfSource:
    default:
      drfSourceStr = "UnknownDrfSource";
      break;
  }
  result["drfSource"] = drfSourceStr;

  // Calculate efficiencies and FWHM at standard energies if valid
  if( drf->isValid() )
  {
    const std::vector<float> energies = { 59.54f, 122.04f, 185.71f, 344.28f, 511.0f,
                                          661.66f, 1000.99f, 1173.23f, 1332.49f, 1836.06f, 2614.53f };

    json efficiencies = json::array();
    for( float energy : energies )
    {
      json effData;
      effData["energy"] = energy;

      try
      {
        const float eff = drf->intrinsicEfficiency( energy );
        effData["intrinsicEfficiency"] = eff;
      }catch( std::exception &e )
      {
        effData["intrinsicEfficiency"] = nullptr;
        effData["efficiencyError"] = e.what();
      }

      if( drf->hasResolutionInfo() )
      {
        try
        {
          const float fwhm = drf->peakResolutionFWHM( energy );
          effData["fwhm"] = fwhm;
        }catch( std::exception &e )
        {
          effData["fwhm"] = nullptr;
          effData["fwhmError"] = e.what();
        }
      }

      efficiencies.push_back( effData );
    }

    result["efficienciesAtStandardEnergies"] = efficiencies;
  }

  return result;
}//nlohmann::json executeGetDetectorInfo(InterSpec* interspec)


nlohmann::json ToolRegistry::executeSearchSourcesByEnergy(const nlohmann::json& params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Parse energies array (required)
  if( !params.contains("energies") || !params["energies"].is_array() || params["energies"].empty() )
    throw std::runtime_error( "The 'energies' parameter must be a non-empty array." );

  const json& energies_array = params["energies"];
  std::vector<double> energies;
  std::vector<double> windows;

  // Get detector type for fallback window calculation
  const bool isHPGe = PeakFitUtils::is_likely_high_res( interspec );

  for( const auto& energy_obj : energies_array )
  {
    if( !energy_obj.contains("energy") )
      throw std::runtime_error( "Each energy object must contain an 'energy' field." );

    const double energy = get_number( energy_obj, "energy" );
    if( energy <= 0.0 )
      throw std::runtime_error( "Energy values must be positive." );

    energies.push_back( energy );

    double window = 10.0; // fallback default
    if( energy_obj.contains("window") )
    {
      window = get_number( energy_obj, "window" );
    }
    else
    {
      // Calculate default window as (3.0/2.35482) * FWHM
      // This corresponds to about ±1.27 sigma, which is kinda within an "okay" energy calbration.
      try
      {
        const double fwhm = AnalystChecks::get_expected_fwhm( energy, interspec );
        window = (3.0 / 2.35482) * fwhm;
      }catch( std::exception & )
      {
        // Failed to get FWHM, will use detector-type-based fallback
        window = isHPGe ? 1.5 : 10.0;
      }
    }//if( window specified ) / else

    if( window < 0.0 )
      throw std::runtime_error( "Window values must be non-negative." );

    windows.push_back( window );
  }//for( loop over energies_array )

  // Parse optional parameters
  // Use i18n key as default since translations may not be available in all environments
  std::string category_str = "isbe-category-nuc-xray";  // "Nuclides + X-rays"
  if( params.contains("source_category") && params["source_category"].is_string() )
    category_str = params["source_category"].get<std::string>();

  std::string min_half_life_str = "100 m";
  if( params.contains("min_half_life") && params["min_half_life"].is_string() )
    min_half_life_str = params["min_half_life"].get<std::string>();

  double min_br = 0.0;
  if( params.contains("min_branching_ratio") )
    min_br = get_number( params, "min_branching_ratio" );

  int max_results = 10;
  if( params.contains("max_results") )
    max_results = static_cast<int>( get_number( params, "max_results" ) );

  if( max_results <= 0 )
    throw std::runtime_error( "max_results must be positive." );

  std::string sort_by = "ProfileScore";
  if( params.contains("sort_by") && params["sort_by"].is_string() )
    sort_by = params["sort_by"].get<std::string>();

  if( (sort_by != "ProfileScore") && (sort_by != "SumEnergyDifference") )
    throw std::runtime_error( "sort_by must be either 'ProfileScore' or 'SumEnergyDifference'." );

  // Parse minimum half-life
  double min_half_life = 0.0;
  try
  {
    min_half_life = PhysicalUnits::stringToTimeDuration( min_half_life_str );
  }catch( std::exception &e )
  {
    throw std::runtime_error( std::string("Invalid min_half_life format: ") + e.what()
                             + ". Examples: '10 s', '5 min', '1 h', '2 y'" );
  }

  // Load category information
  std::vector<IsotopeSearchByEnergy::NucSearchCategory> categories;
  try
  {
    IsotopeSearchByEnergy::init_category_info( categories );
  }catch( std::exception &e )
  {
    throw std::runtime_error( std::string("Failed to initialize search categories: ") + e.what() );
  }

  // Find the matching category
  const IsotopeSearchByEnergy::NucSearchCategory *selected_category = nullptr;
  for( const auto& cat : categories )
  {
    std::string cat_name = cat.m_name.toUTF8();
    std::string cat_key = cat.m_name.key();

    if( (cat_name == category_str) || (cat_key == category_str) )
    {
      selected_category = &cat;
      break;
    }
  }//for( loop over categories )

  if( !selected_category )
  {
    std::string error_msg = "Invalid source_category '" + category_str + "'. Valid categories: ";
    for( size_t i = 0; i < categories.size(); ++i )
    {
      if( i > 0 )
        error_msg += ", ";
      error_msg += "'" + categories[i].m_name.toUTF8() + "'";
    }
    throw std::runtime_error( error_msg );
  }

  // Build radiation source flags from category
  Wt::WFlags<IsotopeSearchByEnergyModel::RadSource> srcs;
  if( selected_category->m_nuclides )
    srcs |= IsotopeSearchByEnergyModel::RadSource::NuclideGammaOrXray;
  if( selected_category->m_fluorescence_xrays )
    srcs |= IsotopeSearchByEnergyModel::RadSource::kFluorescenceXRay;
  if( selected_category->m_reactions )
    srcs |= IsotopeSearchByEnergyModel::RadSource::kReaction;
  if( selected_category->m_alphas )
    srcs |= IsotopeSearchByEnergyModel::RadSource::kAlpha;
  if( selected_category->m_beta_endpoint )
    srcs |= IsotopeSearchByEnergyModel::RadSource::kBetaEndpoint;
  if( selected_category->m_no_progeny )
    srcs |= IsotopeSearchByEnergyModel::RadSource::kNoNuclideProgeny;

  // Create working space for search
  auto workingspace = std::make_shared<IsotopeSearchByEnergyModel::SearchWorkingSpace>();
  workingspace->energies = energies;
  workingspace->windows = windows;
  workingspace->isHPGe = isHPGe;

  // Get foreground spectrum data if available
  std::shared_ptr<SpecMeas> foreground = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( foreground )
  {
    const std::set<int> samplenums = interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );
    auto userpeaks = foreground->peaks( samplenums );
    auto autopeaks = foreground->automatedSearchPeaks( samplenums );

    workingspace->foreground = foreground;
    workingspace->foreground_samplenums = samplenums;

    if( userpeaks )
      workingspace->user_peaks.insert( end(workingspace->user_peaks),
                                       begin(*userpeaks), end(*userpeaks) );
    if( autopeaks )
      workingspace->automated_search_peaks.insert( end(workingspace->automated_search_peaks),
                                                   begin(*autopeaks), end(*autopeaks) );

    workingspace->detector_response_function = foreground->detector();
    workingspace->displayed_measurement = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  }//if( foreground )

  // Perform search synchronously (no threading for LLM tool)
  const std::string appid = "";  // Empty since we're not posting to event loop
  boost::function< void(void) > updatefcn;  // No callback needed

  try
  {
    IsotopeSearchByEnergyModel::setSearchEnergies( workingspace, min_br, min_half_life, srcs,
                                                   selected_category->m_specific_elements,
                                                   selected_category->m_specific_nuclides,
                                                   selected_category->m_specific_reactions,
                                                   appid, updatefcn );
  }catch( std::exception &e )
  {
    throw std::runtime_error( std::string("Search failed: ") + e.what() );
  }

  // Check for errors
  if( !workingspace->error_msg.empty() )
    throw std::runtime_error( "Search error: " + workingspace->error_msg );

  // Sort results
  const IsotopeSearchByEnergyModel::Column sort_column =
    (sort_by == "ProfileScore")
      ? IsotopeSearchByEnergyModel::Column::ProfileDistance
      : IsotopeSearchByEnergyModel::Column::Distance;

  const Wt::SortOrder sort_order = (sort_by == "ProfileScore")
                                     ? Wt::DescendingOrder  // Larger ProfileScore is better
                                     : Wt::AscendingOrder;  // Smaller distance is better

  IsotopeSearchByEnergyModel::sortData( workingspace->matches, energies, sort_column, sort_order );

  // Build JSON output
  json result;
  json sources_array = json::array();

  const size_t num_results = std::min( static_cast<size_t>(max_results), workingspace->matches.size() );

  for( size_t i = 0; i < num_results; ++i )
  {
    const std::vector<IsotopeSearchByEnergyModel::IsotopeMatch>& match_set = workingspace->matches[i];

    if( match_set.empty() )
      continue;

    const IsotopeSearchByEnergyModel::IsotopeMatch& first_match = match_set[0];

    json source_obj;

    // Determine source type and name
    nlohmann::json catagories;
    if( first_match.m_nuclide )
    {
      source_obj["source"] = first_match.m_nuclide->symbol;
      source_obj["source_type"] = "nuclide";
      source_obj["half_life"] = first_match.m_displayData[IsotopeSearchByEnergyModel::ParentHalfLife].toUTF8();
      source_obj["assumed_age"] = first_match.m_displayData[IsotopeSearchByEnergyModel::AssumedAge].toUTF8();
      catagories = source_categories( first_match.m_nuclide, interspec );
    }
    else if( first_match.m_element )
    {
      source_obj["source"] = first_match.m_element->symbol + " x-ray";
      source_obj["source_type"] = "x-ray";
      catagories = source_categories( first_match.m_element, interspec );
    }
    else if( first_match.m_reaction )
    {
      source_obj["source"] = first_match.m_reaction->name();
      source_obj["source_type"] = "reaction";
      catagories = source_categories( first_match.m_reaction, interspec );
    }
    else
    {
      continue; // Unknown source type
    }
    
    source_obj["sum_energy_difference"] = rount_to_hundredth( first_match.m_distance );
    source_obj["profile_score"] = rount_to_hundredth( first_match.m_profileDistance );
    
    if( catagories.is_array() && !catagories.empty() )
      source_obj["sourceCatagories"] = catagories;

    // Build energy matches array
    json energy_matches_array = json::array();

    for( size_t j = 0; j < match_set.size(); ++j )
    {
      const IsotopeSearchByEnergyModel::IsotopeMatch& match = match_set[j];

      json energy_match;
      energy_match["search_energy"] = energies[j];

      const std::string energy_str = match.m_displayData[IsotopeSearchByEnergyModel::Energy].toUTF8();
      energy_match["matched_energy"] = std::stod( energy_str );

      const std::string br_str = match.m_displayData[IsotopeSearchByEnergyModel::BranchRatio].toUTF8();
      energy_match["relative_br"] = std::stod( br_str );

      if( match.m_nuclide && match.m_transition )
      {
        const std::string transition_str =
          match.m_displayData[IsotopeSearchByEnergyModel::SpecificIsotope].toUTF8();
        energy_match["transition"] = transition_str;
      }

      energy_matches_array.push_back( energy_match );
    }//for( loop over match_set )

    source_obj["energy_matches"] = energy_matches_array;
    sources_array.push_back( source_obj );
  }//for( loop over workingspace->matches )

  result["sources"] = sources_array;

  // Add search parameters for reference
  json search_params;
  search_params["energies"] = json::array();
  for( size_t i = 0; i < energies.size(); ++i )
  {
    json e;
    e["energy"] = energies[i];
    e["window"] = windows[i];
    search_params["energies"].push_back( e );
  }
  search_params["category"] = category_str;
  search_params["min_half_life"] = min_half_life_str;
  search_params["min_branching_ratio"] = min_br;
  search_params["max_results"] = max_results;
  search_params["sort_by"] = sort_by;

  result["search_parameters"] = search_params;

  return result;
}//nlohmann::json executeSearchSourcesByEnergy(...)


} // namespace LlmTools

#endif // USE_LLM_INTERFACE
