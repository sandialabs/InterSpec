#include "InterSpec_config.h"
#include "InterSpec/LlmToolRegistry.h"

#if( USE_LLM_INTERFACE )

#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmDeepResearchAgent.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/LlmConversationHistory.h"
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
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/ExternalRidResult.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSearchByEnergy.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/IsotopeSearchByEnergyModel.h"

#include "InterSpec/DoseCalc.h"
#include "InterSpec/GadrasSpecFunc.h"

#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/LlmIsotopicsTool.h"
#include "InterSpec/LlmActivityFitTool.h"
#include "InterSpec/LlmRelActManualTool.h"

#include "Minuit2/MnUserParameters.h"

#include <Wt/WDialog>
#include <Wt/WTextArea>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/Dbo/Dbo>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/EnergyCalibration.h"

#include "InterSpec/EnergyCal.h"
#include "InterSpec/PeakModel.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/ReactionGamma.h"

using namespace std;
using json = nlohmann::json;

namespace {
  // JSON conversion for SpecUtils::SpectrumType enum
  NLOHMANN_JSON_SERIALIZE_ENUM(SpecUtils::SpectrumType, {
      {SpecUtils::SpectrumType::Foreground, "Foreground"},
      {SpecUtils::SpectrumType::Background, "Background"},
      {SpecUtils::SpectrumType::SecondForeground, "Secondary"},
  })

  // JSON conversion for AnalystChecks::EditPeakAction enum
  NLOHMANN_JSON_SERIALIZE_ENUM(AnalystChecks::EditPeakAction, {
      {AnalystChecks::EditPeakAction::SetEnergy, "SetEnergy"},
      {AnalystChecks::EditPeakAction::SetFwhm, "SetFwhm"},
      {AnalystChecks::EditPeakAction::SetAmplitude, "SetAmplitude"},
      {AnalystChecks::EditPeakAction::SetEnergyUncertainty, "SetEnergyUncertainty"},
      {AnalystChecks::EditPeakAction::SetFwhmUncertainty, "SetFwhmUncertainty"},
      {AnalystChecks::EditPeakAction::SetAmplitudeUncertainty, "SetAmplitudeUncertainty"},
      {AnalystChecks::EditPeakAction::SetRoiLower, "SetRoiLower"},
      {AnalystChecks::EditPeakAction::SetRoiUpper, "SetRoiUpper"},
      {AnalystChecks::EditPeakAction::SetSkewType, "SetSkewType"},
      {AnalystChecks::EditPeakAction::SetContinuumType, "SetContinuumType"},
      {AnalystChecks::EditPeakAction::SetSource, "SetSource"},
      {AnalystChecks::EditPeakAction::SetColor, "SetColor"},
      {AnalystChecks::EditPeakAction::SetUserLabel, "SetUserLabel"},
      {AnalystChecks::EditPeakAction::SetUseForEnergyCalibration, "SetUseForEnergyCalibration"},
      {AnalystChecks::EditPeakAction::SetUseForShieldingSourceFit, "SetUseForShieldingSourceFit"},
      {AnalystChecks::EditPeakAction::SetUseForManualRelEff, "SetUseForManualRelEff"},
      {AnalystChecks::EditPeakAction::DeletePeak, "DeletePeak"},
      {AnalystChecks::EditPeakAction::SplitFromRoi, "SplitFromRoi"},
      {AnalystChecks::EditPeakAction::MergeWithLeft, "MergeWithLeft"},
      {AnalystChecks::EditPeakAction::MergeWithRight, "MergeWithRight"},
  })

  // JSON conversion for AnalystChecks::EscapePeakType enum
  NLOHMANN_JSON_SERIALIZE_ENUM(AnalystChecks::EscapePeakType, {
      {AnalystChecks::EscapePeakType::SingleEscape, "SingleEscape"},
      {AnalystChecks::EscapePeakType::DoubleEscape, "DoubleEscape"},
  })

  double rount_to_hundredth(double val){ return 0.01*std::round(100.0*val); }

  /** Some LLMs will give number values as strings, so this function will check types and return the correct answer.

   @param parent The parent JSON object.
   @param name The field name of the number to parse.

   An exception will be thrown if cant be converted to int.
   */
  double get_number( const json& parent, const string &name )
  {
    if( !parent.contains(name) || parent[name].is_null() )
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


  /** Some LLMs will give boolean values as strings, so this function will check types and return the correct answer.

   @param parent The parent JSON object.
   @param name The field name of the boolean to parse.

   An exception will be thrown if cant be converted to boolean.
   */
  bool get_boolean( const json& parent, const string &name )
  {
    if( !parent.contains(name) || parent[name].is_null() )
      throw runtime_error( "'" + name + "' parameter must be specified." );

    if( parent[name].is_boolean() )
      return parent.at(name).get<bool>();

    if( parent[name].is_string() )
    {
      string strval = parent.at(name).get<string>();

      // Convert to lowercase for case-insensitive comparison
      for( char &c : strval )
        c = std::tolower( static_cast<unsigned char>(c) );

      if( strval == "true" || strval == "1" || strval == "yes" )
        return true;

      if( strval == "false" || strval == "0" || strval == "no" )
        return false;

      throw runtime_error( "'" + name + "' parameter must be a boolean (received '" + strval + "')." );
    }//if( parent[name].is_string() )

    // Handle numeric values as booleans (0 = false, non-zero = true)
    if( parent[name].is_number() )
    {
      const double val = parent.at(name).get<double>();
      return (val != 0.0);
    }

    throw runtime_error( "'" + name + "' parameter must be a boolean." );
  }//bool get_boolean( const json& parent, const string &name )


  /** Some LLMs will give number values as strings, so this function will check types and return the correct answer.

   This is similar to get_number, but works on a JSON value directly rather than a named field.

   @param val The JSON value to parse as a double.

   An exception will be thrown if cant be converted to double.
   */
  double get_double( const json& val )
  {
    if( val.is_number() )
      return val.get<double>();

    if( val.is_string() )
    {
      const string strval = val.get<string>();
      double result;
      if( !(stringstream(strval) >> result) )
        throw runtime_error( "Value must be a number (received '" + strval + "')." );
      return result;
    }//if( val.is_string() )

    throw runtime_error( "Value must be a number." );
  }//double get_double( const json& val )


  /** Overload of get_number that returns a default value if the field is not present.

   @param parent The parent JSON object.
   @param name The field name of the number to parse.
   @param default_value The value to return if the field is not present.

   An exception will be thrown if the field exists but cant be converted to a number.
   */
  double get_number( const json& parent, const string &name, const double default_value )
  {
    if( !parent.contains(name) || parent[name].is_null() )
      return default_value;

    return get_number( parent, name );
  }//double get_number( const json& parent, const string &name, const double default_value )


  /** Overload of get_boolean that returns a default value if the field is not present.

   @param parent The parent JSON object.
   @param name The field name of the boolean to parse.
   @param default_value The value to return if the field is not present.

   An exception will be thrown if the field exists but cant be converted to a boolean.
   */
  bool get_boolean( const json& parent, const string &name, const bool default_value )
  {
    if( !parent.contains(name) || parent[name].is_null() )
      return default_value;

    return get_boolean( parent, name );
  }//bool get_boolean( const json& parent, const string &name, const bool default_value )


  /** Overload of get_double that returns a default value if the value is null or missing.

   @param val The JSON value to parse as a double.
   @param default_value The value to return if val is null.

   An exception will be thrown if the value exists but cant be converted to a double.
   */
  double get_double( const json& val, const double default_value )
  {
    if( val.is_null() )
      return default_value;

    return get_double( val );
  }//double get_double( const json& val, const double default_value )

  // Returns the actual key in `params` that matches `key` case-insensitively,
  // or `key` itself if not found (so callers can let contains/operator[] produce a natural error).
  string find_case_insensitive_key( const string &key, const nlohmann::json &params )
  {
    assert( !key.empty() );
    if( params.contains(key) || key.empty() ) //Go for exact match first
      return key;
    
    const bool key_is_padded = (!key.empty() && (std::isspace( static_cast<unsigned char>(key.front()) )
                                                 || std::isspace( static_cast<unsigned char>(key.back()) )));
    
    for( const auto &item : params.items() )
    {
      if( SpecUtils::iequals_ascii( item.key(), key ) )
        return item.key();
      
      if( !key_is_padded
         && !item.key().empty()
         && (std::isspace( static_cast<unsigned char>(item.key().front()) )
             || std::isspace( static_cast<unsigned char>(item.key().back()) )) )
      {
        if( SpecUtils::iequals_ascii( key, SpecUtils::trim_copy(item.key()) ) )
          return item.key();
      }
    }//for( const auto &item : params.items() )
    
    return key;
  };//find_case_insensitive_key

  /** Some LLMs will give complex values (arrays, objects) as JSON strings instead of native types.
   This function normalizes a field by parsing stringified JSON if needed.

   @param params The parent JSON object (non-const reference, will be modified if field is a string).
   @param field_name The name of the field to normalize.

   If the field is a string, it will attempt to parse it as JSON and replace the field value.
   If parsing fails, the original string value remains.
   */
  void normalize_json_field( json& params, const string& field_name )
  {
    if( !params.contains(field_name) )
      return;

    if( !params[field_name].is_string() )
      return; // Already a native type, no normalization needed

    string str_value = params[field_name].get<string>();
    SpecUtils::trim( str_value );

    // Try to parse as JSON
    try
    {
      json parsed = json::parse( str_value );
      params[field_name] = std::move(parsed);
    }catch( const json::parse_error &e )
    {
      // If parsing fails, leave as string (might be intended as a string value)
      // Let the downstream code handle validation
      cerr << "Failed to convert string into a JSON object: trimmed='" << str_value << "'." << endl;
      cerr << endl;
    }
  }//void normalize_json_field( json& params, const string& field_name )

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
    p.nonBackgroundPeaksOnly = get_boolean( j, "NonBackgroundPeaksOnly", false );

    // Optional energy bounds
    p.lowerEnergy.reset();
    p.upperEnergy.reset();

    if( j.contains("lowerEnergy") && !j["lowerEnergy"].is_null() )
      p.lowerEnergy = get_number( j, "lowerEnergy" );

    if( j.contains("upperEnergy") && !j["upperEnergy"].is_null() )
      p.upperEnergy = get_number( j, "upperEnergy" );

    if( p.lowerEnergy && p.upperEnergy && (*p.upperEnergy < *p.lowerEnergy) )
    {
      throw std::runtime_error( "upperEnergy (" + std::to_string(*p.upperEnergy)
                                + " keV) must be greater than or equal to lowerEnergy ("
                                + std::to_string(*p.lowerEnergy) + " keV)" );
    }
  }

  void from_json(const json& j, AnalystChecks::FitPeakOptions& p) {

    p.energy = get_number( j, "energy" );

    p.doNotAddToAnalysisPeaks = get_boolean( j, "DoNotAddToAnalysisPeaks", false );

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

    p.source = std::nullopt;
    if( j.contains("source") )
      p.source = j["source"].get<std::string>();
    else if( j.contains("nuclide") )
      p.source = j["nuclide"].get<std::string>();
    else if( j.contains("element") )
      p.source = j["element"].get<std::string>();
    else if( j.contains("xray") )
      p.source = j["xray"].get<std::string>();
    else if( j.contains("x-ray") )
      p.source = j["x-ray"].get<std::string>();
    else if( j.contains("reaction") )
      p.source = j["reaction"].get<std::string>();
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
    
    auto source_type_str = []( const PeakDef::SourceGammaType type ) -> const char * {
      switch( type )
      {
        case PeakDef::NormalGamma:
          return "gamma";
        case PeakDef::AnnihilationGamma:
          return "Annih.";
        case PeakDef::SingleEscapeGamma:
          return "S.E.";
        case PeakDef::DoubleEscapeGamma:
          return "D.E.";
        case PeakDef::XrayGamma:
          return "x-ray";
      }//
      assert( 0 );
      return nullptr;
    };
    
    
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
      
      const char * const gamma_type = source_type_str( peak->sourceGammaType() );
      if( gamma_type )
        src["photonType"] = gamma_type;
      
      src["energy"] = peak->gammaParticleEnergy();
    }else if( const SandiaDecay::Element * const el = peak->xrayElement() )
    {
      peak_json["source"]["element"] = el->symbol;
      peak_json["source"]["photonType"] = "x-ray";
      peak_json["source"]["energy"] = peak->xrayEnergy();
    }else if( const ReactionGamma::Reaction * const rctn = peak->reaction() )
    {
      peak_json["source"]["reaction"] = rctn->name();
      peak_json["source"]["energy"] = peak->reactionEnergy();
      const char * const gamma_type = source_type_str( peak->sourceGammaType() );
      if( gamma_type )
        peak_json["source"]["photonType"] = gamma_type;
    }

    // Add boolean flags for peak usage - use raw user preferences, not computed values
    peak_json["useForEnergyCalibration"] = peak->useForEnergyCalibrationUserPreference();
    peak_json["useForShieldingSourceFit"] = peak->useForShieldingSourceFitUserPreference();
    peak_json["useForManualRelEff"] = peak->useForManualRelEffUserPreference();
  }
  
  void to_json( json &roi_json,
               const shared_ptr<const PeakContinuum> &cont,
               const vector<shared_ptr<const PeakDef>> &peaks,
               const shared_ptr<const SpecUtils::Measurement> &meas,
               const std::vector<std::shared_ptr<const PeakDef>> *analysis_peaks = nullptr ){
    roi_json = json{
      {"lowerEnergy", rount_to_hundredth(cont->lowerEnergy())},
      {"upperEnergy", rount_to_hundredth(cont->upperEnergy())},
      {"continuumType", PeakContinuum::offset_type_str(cont->type()) }
      //("continuumCounts", cont->offset_integral( cont->lowerEnergy(), cont->upperEnergy(), dataH ) }
    };

    for( const shared_ptr<const PeakDef> &peak : peaks ) {
      json peak_json;

      to_json( peak_json, peak, meas );

      // Check if this peak is in the analysis_peaks list
      if( analysis_peaks )
      {
        const bool is_analysis_peak = std::find(analysis_peaks->begin(), analysis_peaks->end(), peak) != analysis_peaks->end();
        peak_json["isAnalysisPeak"] = is_analysis_peak;
      }

      roi_json["peaks"].push_back( peak_json );
    }
  }//void to_json( json &roi_json, const shared_ptr<const PeakContinuum> &cont, const vector<shared_ptr<const PeakDef>> &peaks )

  
  void to_json( json &peak_rois,
               const std::vector<std::shared_ptr<const PeakDef>> &peaks,
               const shared_ptr<const SpecUtils::Measurement> &meas,
               const std::vector<std::shared_ptr<const PeakDef>> *analysis_peaks = nullptr ){
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

      to_json( roi_json, cont, peaks, meas, analysis_peaks );

      roi_json["roiID"] = static_cast<int>(roi_index);

      peak_rois.push_back(roi_json);
    }
  }
  
  
  void to_json(json& j,
               const AnalystChecks::DetectedPeakStatus& p,
               const shared_ptr<const SpecUtils::Measurement> &meas ) {
    json peak_rois;
    to_json( peak_rois, p.peaks, meas, &p.analysis_peaks );

    j = json{{"rois", peak_rois}};
  }//void to_json(json& j, const AnalystChecks::DetectedPeakStatus& p) {
  
  void to_json(json& j, const AnalystChecks::FitPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {

    const std::shared_ptr<const PeakDef> &fitPeak = p.fitPeak;
    const std::vector<std::shared_ptr<const PeakDef>> &peaksInRoi = p.peaksInRoi;

    j = json::object();

    if( fitPeak )
    {
      json roi_json;

      to_json( roi_json, fitPeak->continuum(), peaksInRoi, meas );

      j["roi"] = roi_json;

      //json peak_json;
      //to_json( peak_json, fitPeak );
      j["fitPeakEnergy"] = rount_to_hundredth(fitPeak->mean());

      j["status"] = "success";

      //"response":"OK"
      //"status":200
    }else
    {
      j["error"] = "No peak fit.";
      j["status"] = "failed";
    }
  }//void to_json(json& j, const AnalystChecks::DetectedPeakStatus& p)
  
  
  SpecUtils::SpectrumType parse_spectrum_type( const json &j, const std::string &field_name = "specType" )
  {
    std::string specTypeStr = j.value(field_name, std::string());
    if( specTypeStr.empty() || specTypeStr == "Foreground" )
    {
      return SpecUtils::SpectrumType::Foreground;
    }else if( specTypeStr == "Background" )
    {
      return SpecUtils::SpectrumType::Background;
    }else if( specTypeStr == "Secondary" )
    {
      return SpecUtils::SpectrumType::SecondForeground;
    }else
    {
      throw std::runtime_error( "Invalid spectrum type: " + specTypeStr );
    }
  }

  void from_json(const json& j, AnalystChecks::GetUserPeakOptions& p) {
    p.specType = parse_spectrum_type( j );

    // Parse optional energy range
    if( j.contains("lowerEnergy") && !j["lowerEnergy"].is_null() )
      p.lowerEnergy = get_number( j, "lowerEnergy" );

    if( j.contains("upperEnergy") && !j["upperEnergy"].is_null() )
      p.upperEnergy = get_number( j, "upperEnergy" );
  }

  void from_json(const json& j, AnalystChecks::FitPeaksForNuclideOptions& p) {
    // source may be absent, null, empty string, empty array, or an array containing nulls
    // when FitNormBkgrndPeaks is set (fit only NORM background peaks).
    const bool has_source = j.contains("source") && !j["source"].is_null();
    if( has_source )
    {
      const json& sourcesParam = j["source"];
      if( sourcesParam.is_string() && !sourcesParam.get<std::string>().empty() )
      {
        p.sources = { sourcesParam.get<std::string>() };
      }else if( sourcesParam.is_array() )
      {
        for( const json &entry : sourcesParam )
        {
          if( entry.is_string() && !entry.get<std::string>().empty() )
            p.sources.push_back( entry.get<std::string>() );
          else if( !entry.is_null() )
            throw std::runtime_error( "source array entries must be strings or null" );
        }
      }else if( !sourcesParam.is_string() ) // non-empty string already handled above
      {
        throw std::runtime_error( "Invalid sources parameter: must be string or array of strings" );
      }
    }//if( has_source )

    // Parse options array; also accept legacy field names for backwards compatibility
    const std::string options_field = j.contains("options") ? "options"
                                    : j.contains("fitSrcPeaksOptions") ? "fitSrcPeaksOptions"
                                    : std::string();

    if( !options_field.empty() && j[options_field].is_array() )
    {
      static const std::map<std::string, FitPeaksForNuclides::FitSrcPeaksOptions> flag_map = {
        { "DoNotUseExistingRois",          FitPeaksForNuclides::FitSrcPeaksOptions::DoNotUseExistingRois },
        { "RefitInterferingAnalysisPeaks", FitPeaksForNuclides::FitSrcPeaksOptions::ExistingPeaksAsFreePeak },
        { "DoNotVaryEnergyCal",            FitPeaksForNuclides::FitSrcPeaksOptions::DoNotVaryEnergyCal },
        { "DoNotRefineEnergyCal",          FitPeaksForNuclides::FitSrcPeaksOptions::DoNotRefineEnergyCal },
        { "FitNormPeaks",                  FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks },
        { "FitNormBkgrndPeaksDontUse",     FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaksDontUse },
      };

      for( const json &flag_val : j[options_field] )
      {
        if( !flag_val.is_string() )
          throw std::runtime_error( "options entries must be strings" );
        const std::string flag_str = flag_val.get<std::string>();

        if( flag_str == "DoNotAddToUserSession" )
        {
          p.doNotAddPeaksToUserSession = true;
          continue;
        }

        const auto it = flag_map.find( flag_str );
        if( it == flag_map.end() )
          throw std::runtime_error( "Unknown options flag: '" + flag_str + "'" );
        p.fitSrcPeaksOptions |= it->second;
      }
    }//if options present

    // Legacy standalone field; options array takes precedence if both present
    if( options_field != "options" )
      p.doNotAddPeaksToUserSession = get_boolean( j, "doNotAddPeaksToUserSession", false );

    // Validate: sources may only be empty when FitNormBkgrndPeaks is set
    const bool fit_norm = (p.fitSrcPeaksOptions & FitPeaksForNuclides::FitSrcPeaksOptions::FitNormBkgrndPeaks);
    if( p.sources.empty() && !fit_norm )
      throw std::runtime_error( "source must be specified unless FitNormBkgrndPeaks is included in fitSrcPeaksOptions" );
  }

  void from_json( const json &j, AnalystChecks::EditAnalysisPeakOptions &p )
  {
    p.energy = get_number( j, "energy" );

    // Parse editAction manually from string
    const std::string action_str = j.at("editAction").get<std::string>();
    p.editAction = AnalystChecks::edit_peak_action_from_string( action_str );

    // Parse spectrum type (default to Foreground)
    std::string specTypeStr = j.value("specType", std::string());
    if( specTypeStr.empty() || specTypeStr == "Foreground" )
    {
      p.specType = SpecUtils::SpectrumType::Foreground;
    }else if( specTypeStr == "Background" )
    {
      p.specType = SpecUtils::SpectrumType::Background;
    }else if( specTypeStr == "Secondary" )
    {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    }else
    {
      throw std::runtime_error( "Invalid spectrum type: " + specTypeStr );
    }

    // Get optional values
    if( j.contains("doubleValue") && !j["doubleValue"].is_null() )
      p.doubleValue = get_double( j["doubleValue"] );

    if( j.contains("stringValue") && !j["stringValue"].is_null() )
      p.stringValue = j["stringValue"].get<std::string>();

    if( j.contains("boolValue") && !j["boolValue"].is_null() )
      p.boolValue = get_boolean( j, "boolValue" );

    if( j.contains("uncertainty") && !j["uncertainty"].is_null() )
      p.uncertainty = get_number( j, "uncertainty" );
  }//void from_json( const json &j, AnalystChecks::EditAnalysisPeakOptions &p )

  void to_json( json &j, const AnalystChecks::EditAnalysisPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas )
  {
    j = json{
      {"success", p.success},
      {"message", p.message}
    };

    if( p.modifiedPeak )
    {
      json peak_info;
      to_json( peak_info, *p.modifiedPeak, meas );
      j["modifiedPeak"] = peak_info;
    }

    if( !p.peaksInRoi.empty() )
    {
      json roi_peaks = json::array();
      for( const auto &peak : p.peaksInRoi )
      {
        json peak_json;
        to_json( peak_json, peak, meas );
        roi_peaks.push_back( peak_json );
      }
      j["peaksInRoi"] = roi_peaks;
    }
  }//void to_json( json &j, const AnalystChecks::EditAnalysisPeakStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas )

  void from_json( const json &j, AnalystChecks::EscapePeakCheckOptions &p )
  {
    p.energy = get_number( j, "energy" );

    // Parse spectrum type (default to Foreground)
    std::string specTypeStr = j.value("specType", std::string());
    if( specTypeStr.empty() || specTypeStr == "Foreground" )
    {
      p.specType = SpecUtils::SpectrumType::Foreground;
    }else if( specTypeStr == "Background" )
    {
      p.specType = SpecUtils::SpectrumType::Background;
    }else if( specTypeStr == "Secondary" )
    {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    }else
    {
      throw std::runtime_error( "Invalid spectrum type: " + specTypeStr );
    }
  }//void from_json( const json &j, AnalystChecks::EscapePeakCheckOptions &p )

  void to_json( json &j, const AnalystChecks::EscapePeakCheckStatus &p )
  {
    j = json{
      {"potentialSingleEscapePeakEnergy", rount_to_hundredth(p.potentialSingleEscapePeakEnergy)},
      {"potentialDoubleEscapePeakEnergy", rount_to_hundredth(p.potentialDoubleEscapePeakEnergy)},
      {"potentialParentPeakSingleEscape", rount_to_hundredth(p.potentialParentPeakSingleEscape)},
      {"potentialParentPeakDoubleEscape", rount_to_hundredth(p.potentialParentPeakDoubleEscape)},
      //{"singleEscapeSearchWindow", rount_to_hundredth(p.singleEscapeSearchWindow)},
      //{"doubleEscapeSearchWindow", rount_to_hundredth(p.doubleEscapeSearchWindow)},
      //{"singleEscapeParentSearchWindow", rount_to_hundredth(p.singleEscapeParentSearchWindow)},
      //{"doubleEscapeParentSearchWindow", rount_to_hundredth(p.doubleEscapeParentSearchWindow)}
    };

    if( p.parentPeak.has_value() )
    {
      json parent_info = json{
        {"parentPeakEnergy", rount_to_hundredth(p.parentPeak->parentPeakEnergy)},
        {"escapeType", AnalystChecks::to_string(p.parentPeak->escapeType) }
      };

      if( p.parentPeak->sourceLabel.has_value() )
        parent_info["sourceLabel"] = *p.parentPeak->sourceLabel;

      if( p.parentPeak->userLabel.has_value() )
        parent_info["userLabel"] = *p.parentPeak->userLabel;

      j["parentPeak"] = parent_info;
    }
  }//void to_json( json &j, const AnalystChecks::EscapePeakCheckStatus &p )

  void from_json( const json &j, AnalystChecks::SumPeakCheckOptions &p )
  {
    p.energy = get_number( j, "energy" );

    // Parse spectrum type (default to Foreground)
    std::string specTypeStr = j.value("specType", std::string());
    if( specTypeStr.empty() || specTypeStr == "Foreground" )
    {
      p.specType = SpecUtils::SpectrumType::Foreground;
    }else if( specTypeStr == "Background" )
    {
      p.specType = SpecUtils::SpectrumType::Background;
    }else if( specTypeStr == "Secondary" )
    {
      p.specType = SpecUtils::SpectrumType::SecondForeground;
    }else
    {
      throw std::runtime_error( "Invalid spectrum type: " + specTypeStr );
    }

    // Parse optional distance parameter
    if( j.contains("distance") && !j["distance"].is_null() )
    {
      const std::string distanceStr = j["distance"].get<std::string>();
      try
      {
        p.distance = PhysicalUnits::stringToDistance( distanceStr );
      }catch( const std::exception &e )
      {
        throw std::runtime_error( "Invalid distance value '" + distanceStr + "': " + e.what() );
      }
    }
  }//void from_json( const json &j, AnalystChecks::SumPeakCheckOptions &p )

  void to_json( json &j, const AnalystChecks::SumPeakCheckStatus &p )
  {
    j = json{
      {"searchWindow", rount_to_hundredth(p.searchWindow)}
    };

    if( p.sumPeakInfo.has_value() )
    {
      json sum_info = json{
        {"sumType", AnalystChecks::to_string(p.sumPeakInfo->sumType) }
      };

      // Serialize first peak information
      if( p.sumPeakInfo->firstPeak )
      {
        json first_peak = json{
          {"energy", rount_to_hundredth(p.sumPeakInfo->firstPeak->mean())},
          {"area", rount_to_hundredth(p.sumPeakInfo->firstPeak->peakArea())}
        };

        // Add source information if available
        if( p.sumPeakInfo->firstPeak->parentNuclide() )
          first_peak["source"] = p.sumPeakInfo->firstPeak->parentNuclide()->symbol;
        else if( p.sumPeakInfo->firstPeak->reaction() )
          first_peak["source"] = p.sumPeakInfo->firstPeak->reaction()->name();
        else if( p.sumPeakInfo->firstPeak->xrayElement() )
          first_peak["source"] = p.sumPeakInfo->firstPeak->xrayElement()->symbol;

        sum_info["firstPeak"] = first_peak;
      }

      // Serialize second peak information
      if( p.sumPeakInfo->secondPeak )
      {
        json second_peak = json{
          {"energy", rount_to_hundredth(p.sumPeakInfo->secondPeak->mean())},
          {"area", rount_to_hundredth(p.sumPeakInfo->secondPeak->peakArea())}
        };

        // Add source information if available
        if( p.sumPeakInfo->secondPeak->parentNuclide() )
          second_peak["source"] = p.sumPeakInfo->secondPeak->parentNuclide()->symbol;
        else if( p.sumPeakInfo->secondPeak->reaction() )
          second_peak["source"] = p.sumPeakInfo->secondPeak->reaction()->name();
        else if( p.sumPeakInfo->secondPeak->xrayElement() )
          second_peak["source"] = p.sumPeakInfo->secondPeak->xrayElement()->symbol;

        sum_info["secondPeak"] = second_peak;
      }

      // Add labels if available
      if( p.sumPeakInfo->userLabel.has_value() )
        sum_info["userLabel"] = *p.sumPeakInfo->userLabel;

      // Add coincidence fraction for cascade sums
      if( p.sumPeakInfo->coincidenceFraction.has_value() )
        sum_info["coincidenceFraction"] = *p.sumPeakInfo->coincidenceFraction;

      j["sumPeakInfo"] = sum_info;
    }
  }//void to_json( json &j, const AnalystChecks::SumPeakCheckStatus &p )

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

    j = json{{"rois", peak_rois}};
  }

  void to_json( json &j, const std::vector<std::string> &sources )
  {
    if( sources.empty() )
      j = json{{"sources", "none"}};
    else
      j = json{{"sources", sources}};
  }

  /** Helper function to create simplified peak JSON with essential fields.
   
   Creates a JSON object for a peak containing: energy, fwhm, amplitude, amplitudeUncert, numSigma, and cps.
   
   @param peak The peak to convert to simplified JSON
   @param meas Optional measurement data (used to calculate cps from live_time)
   @param include_source If true, includes source information (nuclide, element, or reaction) if assigned to the peak
   @returns JSON object with simplified peak data, or empty JSON object if peak is not GaussianDefined
   */
  json peak_to_simplified_json( const shared_ptr<const PeakDef> &peak, 
                                 const shared_ptr<const SpecUtils::Measurement> &meas,
                                 const bool include_source )
  {
    json peak_json;
    
    if( !peak || (peak->type() != PeakDef::DefintionType::GaussianDefined) )
      return peak_json;
    
    peak_json["energy"] = rount_to_hundredth( peak->mean() );
    peak_json["fwhm"] = rount_to_hundredth( peak->fwhm() );
    peak_json["amplitude"] = rount_to_hundredth( peak->amplitude() );
    peak_json["amplitudeUncert"] = rount_to_hundredth( peak->amplitudeUncert() );
    
    if( peak->amplitudeUncert() > 0.0 )
      peak_json["numSigma"] = rount_to_hundredth( peak->amplitude() / peak->amplitudeUncert() );
    
    if( meas && (meas->live_time() > 0.0) )
    {
      const double cps = peak->amplitude() / meas->live_time();
      peak_json["cps"] = rount_to_hundredth( cps );
    }
    
    if( include_source )
    {
      // Include source information if assigned
      if( const SandiaDecay::Nuclide * const nuc = peak->parentNuclide() )
        peak_json["source"] = nuc->symbol;
      else if( const SandiaDecay::Element * const el = peak->xrayElement() )
        peak_json["source"] = el->symbol;
      else if( const ReactionGamma::Reaction * const rctn = peak->reaction() )
        peak_json["source"] = rctn->name();
    }
    
    return peak_json;
  }

  void to_json(json& j, const AnalystChecks::FitPeaksForNuclideStatus &p, const shared_ptr<const SpecUtils::Measurement> &meas ) {
    // Return simplified peak information instead of full verbose format
    json peaks_array = json::array();
    for( const shared_ptr<const PeakDef> &peak : p.fitPeaks )
    {
      json peak_json = peak_to_simplified_json( peak, meas, true );
      if( !peak_json.empty() )
        peaks_array.push_back( peak_json );
    }
    
    j = json{{"peaks", peaks_array}};
    if( !p.warnings.empty() )
      j["warnings"] = p.warnings;
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


  /** Helper function to get the GadrasScatterTable for dose calculations.
   *  Uses a static shared pointer with mutex protection for thread-safe lazy initialization.
   *  The scatter table is loaded from GadrasContinuum.lib in the static data directory.
   */
  mutex sm_dose_scatter_mutex;
  bool sm_tried_init_dose_scatter = false;
  shared_ptr<const GadrasScatterTable> sm_dose_scatter;

  shared_ptr<const GadrasScatterTable> getDoseCalcScatterTable()
  {
    lock_guard<mutex> lock( sm_dose_scatter_mutex );
    if( sm_tried_init_dose_scatter )
      return sm_dose_scatter;

    sm_tried_init_dose_scatter = true;

    try
    {
      const string data = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GadrasContinuum.lib" );
      sm_dose_scatter = make_shared<GadrasScatterTable>( data );
    }catch( exception &e )
    {
      cerr << "LlmToolRegistry: Failed to init GADRAS scatter for dose calc: " << e.what() << endl;
    }

    return sm_dose_scatter;
  }//getDoseCalcScatterTable()



}//namespace

namespace LlmTools
{

ToolRegistry::ToolRegistry( const LlmConfig &config )
{
  registerDefaultTools( config );
}
  

void ToolRegistry::registerTool(const SharedTool& tool) {
  m_tools[tool.name] = tool;
}

SharedTool ToolRegistry::createToolWithExecutor( const std::string &toolName )
{
  SharedTool tool;
  tool.name = toolName;

  // Assign executor based on tool name
  if( toolName == "get_detected_peaks" )
  {
    tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
      return executePeakDetection(params, interspec);
      };
    }else if( toolName == "add_analysis_peak" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executePeakFit(params, interspec);
      };
    }else if( toolName == "edit_analysis_peak" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeEditAnalysisPeak(params, interspec);
      };
    }else if( toolName == "escape_peak_check" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeEscapePeakCheck(params, interspec);
      };
    }else if( toolName == "sum_peak_check" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeSumPeakCheck(params, interspec);
      };
    }else if( toolName == "get_analysis_peaks" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetUserPeaks(params, interspec);
      };
    }else if( toolName == "get_identified_sources" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetIdentifiedSources(params, interspec);
      };
    }else if( toolName == "get_unidentified_peaks" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetUnidentifiedDetectedPeaks(params, interspec);
      };
    }else if( toolName == "get_spectrum_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetSpectrumInfo(params, interspec);
      };
    }else if( toolName == "primary_gammas_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetCharacteristicGammasForSource(params);
      };
    }else if( toolName == "sources_with_primary_gammas_in_energy_range" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
      };
    }else if( toolName == "sources_with_primary_gammas_near_energy" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
      };
    }
#if( !INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
    else if( toolName == "sources_associated_with_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetAssociatedSources(params);
      };
    }else if( toolName == "analyst_notes_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetSourceAnalystNotes(params);
      };
    }
#endif
    else if( toolName == "source_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetSourceInfo(params, interspec);
      };
    }else if( toolName == "nuclide_decay_chain" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetNuclideDecayChain(params);
      };
    }else if( toolName == "automated_source_id_results" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetAutomatedRiidId(params, interspec);
      };
    }else if( toolName == "loaded_spectra" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetLoadedSpectra(params, interspec);
      };
    }else if( toolName == "add_analysis_peaks_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeFitPeaksForNuclide(params, interspec);
      };
    }else if( toolName == "get_counts_in_energy_range" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetCountsInEnergyRange(params, interspec);
      };
    }else if( toolName == "get_expected_fwhm" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetExpectedFwhm(params, interspec);
      };
    }else if( toolName == "currie_mda_calc" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeCurrieMdaCalc(params, interspec);
      };
    }else if( toolName == "calculate_dose" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeCalculateDose(params, interspec);
      };
    }else if( toolName == "source_photons" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetSourcePhotons(params);
      };
    }else if( toolName == "photopeak_detection_efficiency" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executePhotopeakDetectionCalc(params, interspec);
      };
    }else if( toolName == "get_materials" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetMaterials(interspec);
      };
    }else if( toolName == "get_material_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetMaterialInfo(params, interspec);
      };
    }else if( toolName == "available_detector_efficiency_functions" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeAvailableDetectors(params, interspec);
      };
    }else if( toolName == "load_detector_efficiency_function" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeLoadDetectorEfficiency(params, interspec);
      };
    }else if( toolName == "detector_efficiency_function_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeGetDetectorInfo(params, interspec);
      };
    }else if( toolName == "search_sources_by_energy" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeSearchSourcesByEnergy(params, interspec);
      };
    }else if( toolName == "activity_fit" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeActivityFit(params, interspec);
      };
    }else if( toolName == "activity_fit_one_off" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeActivityFitOneOff(params, interspec);
      };
    }else if( toolName == "get_shielding_source_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeGetShieldingSourceConfig(params, interspec);
      };
    }else if( toolName == "modify_shielding_source_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeModifyShieldingSourceConfig(params, interspec);
      };
    }else if( toolName == "mark_peaks_for_activity_fit" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeMarkPeaksForActivityFit(params, interspec);
      };
    }else if( toolName == "close_activity_shielding_display" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeCloseActivityShieldingDisplay(params, interspec);
      };
    }else if( toolName == "ask_user_question" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return ActivityFitTool::executeAskUserQuestion(params, interspec);
      };
    }else if( toolName == "set_workflow_state" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return executeSetWorkflowState(params, interspec);
      };
    }else if( toolName == "list_isotopics_presets" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeListIsotopicsPresets(params, interspec);
      };
    }else if( toolName == "get_isotopics_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeGetIsotopicsConfig(params, interspec);
      };
    }else if( toolName == "reset_isotopics_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeResetIsotopicsConfig(params, interspec);
      };
    }else if( toolName == "load_isotopics_preset" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeLoadIsotopicsPreset(params, interspec);
      };
    }else if( toolName == "perform_isotopics_calculation" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executePerformIsotopics(params, interspec);
      };
    }else if( toolName == "modify_isotopics_nuclides" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeModifyIsotopicsNuclides(params, interspec);
      };
    }else if( toolName == "modify_isotopics_rois" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeModifyIsotopicsRois(params, interspec);
      };
    }else if( toolName == "modify_isotopics_curve_settings" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeModifyIsotopicsCurveSettings(params, interspec);
      };
    }else if( toolName == "modify_isotopics_options" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeModifyIsotopicsOptions(params, interspec);
      };
    }else if( toolName == "modify_isotopics_constraints" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeModifyIsotopicsConstraints(params, interspec);
      };
    }else if( toolName == "get_isotopics_config_schema" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return IsotopicsTool::executeGetIsotopicsConfigSchema(params, interspec);
      };
    }else if( toolName == "peak_based_relative_efficiency" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return RelActManualTool::executePeakBasedRelativeEfficiency(params, interspec);
      };
    }else if( toolName == "get_rel_act_manual_state" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
        return RelActManualTool::executeGetRelActManualState(params, interspec);
      };
    }else if( toolName == "create_peak_checkpoint" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory* history) -> json {
        return executeCreatePeakCheckpoint(params, interspec, history);
      };
    }else if( toolName == "restore_peaks_to_checkpoint" )
    {
      tool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory* history) -> json {
        return executeRestorePeaksToCheckpoint(params, interspec, history);
      };
    }else
    {
      throw std::runtime_error( "Unknown tool name: " + toolName );
    }

    return tool;
}//ToolRegistry::createToolWithExecutor(...)


void ToolRegistry::registerDefaultTools( const LlmConfig &config )
{
  if( !m_tools.empty() )
    return;

  // Get tool configurations from config (already loaded from llm_tools_config.xml)
  if( config.tools.empty() )
    throw std::runtime_error("No tools configuration provided - cannot initialize LLM interface");
  
  //cout << "Registering default LLM tools..." << endl;

  // Register agent-specific invoke tools dynamically from config
  for( const LlmConfig::AgentConfig &agent : config.agents )
  {
    // Skip MainAgent - it doesn't get its own invoke tool
    if( agent.type == AgentType::MainAgent )
      continue;
    
    // Create invoke_<AgentName> tool
    const string toolName = "invoke_" + agent.name;

    SharedTool invokeTool;
    invokeTool.name = toolName;
    invokeTool.description = agent.description;
    
    // Schema for invoke tool
    invokeTool.parameters_schema = json::parse(R"({
        "type": "object",
        "properties": {
          "context": {
            "type": "string",
            "description": "Context about the current analysis state."
          },
          "task": {
            "type": "string",
            "description": "The specific task for the sub-agent to perform."
          }
        },
        "required": ["context", "task"]
      })");
    
    // Executor is a placeholder - actual invocation handled in executeToolCalls
    invokeTool.executor = [](const json& params, InterSpec* interspec, shared_ptr<LlmInteraction>, LlmConversationHistory*) -> json {
      assert( 0 );
      throw std::runtime_error( "invoke_* tools should be handled by executeToolCalls, not called directly" );
    };

    // Set which agents can invoke this sub-agent
    // If agent.availableForAgents is empty, default to MainAgent only for backward compatibility
    if( agent.availableForAgents.empty() )
      invokeTool.availableForAgents = {AgentType::MainAgent};
    else
      invokeTool.availableForAgents = agent.availableForAgents;

    // Prevent recursive invoke loops by removing self-invocation.
    invokeTool.availableForAgents.erase(
      std::remove( begin(invokeTool.availableForAgents), end(invokeTool.availableForAgents), agent.type ),
      end(invokeTool.availableForAgents)
    );

    registerTool(invokeTool);
  }//for( loop over agents in config )
  
  // Register tools from configs
  for( const LlmConfig::ToolConfig &toolConfig : config.tools )
  {
    try
    {
      // Create tool with executor
      SharedTool tool = createToolWithExecutor( toolConfig.name );

      // Apply description and schema from config
      if( !toolConfig.defaultDescription.empty() )
        tool.description = toolConfig.defaultDescription;

      // Apply parameter schema (already parsed and validated; {} is a valid empty schema)
      if( !toolConfig.parametersSchema.is_null() )
        tool.parameters_schema = toolConfig.parametersSchema;

      // Apply role-specific descriptions
      tool.roleDescriptions = toolConfig.roleDescriptions;

      // Apply agent restrictions
      tool.availableForAgents = toolConfig.availableForAgents;

      registerTool(tool);
    }catch( const std::exception &e )
    {
      cerr << "Warning: Failed to create tool '" << toolConfig.name << "': " << e.what() << endl;
    }
  }//for( loop over tool configs )

  bool loaded_deep_research_skills = false;
  LlmDeepResearch::registerDeepResearchTools( config.llmApi.deep_research_url,
                                              [this]( const SharedTool &tool ){ registerTool(tool); },
                                              loaded_deep_research_skills );

  if( m_tools.find("invoke_DeepResearch") == end(m_tools) )
    cerr << "Warning: DeepResearch agent invoke tool was not registered. Check llm_agents.xml configuration." << endl;

  if( !loaded_deep_research_skills )
    cerr << "Info: No DeepResearch SKILL.md files were found under llm_knowledge directories." << endl;

  // NOTE: Hard-coded fallback tool definitions have been removed.
  // If no tools are loaded from XML, we throw an exception above.
  // This ensures that:
  //   1. The XML configuration is the single source of truth for tool definitions
  //   2. Tool descriptions and schemas can be updated without recompilation
  //   3. Missing configuration is caught early rather than silently using outdated hardcoded values

  // Runtime validation: Check that all expected tools are registered
  std::vector<std::string> expectedTools = {
    "get_detected_peaks",
    "add_analysis_peak",
    "edit_analysis_peak",
    "escape_peak_check",
    "sum_peak_check",
    "get_analysis_peaks",
    "get_identified_sources",
    "get_unidentified_peaks",
    "get_spectrum_info",
    "primary_gammas_for_source",
    "sources_with_primary_gammas_in_energy_range",
    "sources_with_primary_gammas_near_energy",
#if( !INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
    "sources_associated_with_source",
    "analyst_notes_for_source",
#endif
    "source_info",
    "nuclide_decay_chain",
    "automated_source_id_results",
    "loaded_spectra",
    "add_analysis_peaks_for_source",
    "get_counts_in_energy_range",
    "get_expected_fwhm",
    "currie_mda_calc",
    "calculate_dose",
    "source_photons",
    "photopeak_detection_efficiency",
    "get_materials",
    "get_material_info",
    "available_detector_efficiency_functions",
    "load_detector_efficiency_function",
    "detector_efficiency_function_info",
    "search_sources_by_energy",
    "activity_fit",
    "activity_fit_one_off",
    "ask_user_question",
    "close_activity_shielding_display",
    "get_shielding_source_config",
    "mark_peaks_for_activity_fit",
    "modify_shielding_source_config",
    "get_isotopics_config",
    "get_isotopics_config_schema",
    "list_isotopics_presets",
    "load_isotopics_preset",
    "modify_isotopics_constraints",
    "modify_isotopics_curve_settings",
    "modify_isotopics_nuclides",
    "modify_isotopics_options",
    "modify_isotopics_rois",
    "perform_isotopics_calculation",
    "reset_isotopics_config",
    "set_workflow_state",
    "get_rel_act_manual_state",
    "peak_based_relative_efficiency"
  };

  if( !config.llmApi.deep_research_url.empty() )
    expectedTools.push_back( "query_deep_research_endpoint" );

  std::vector<std::string> missingTools;
  for( const std::string &toolName : expectedTools )
  {
    if( m_tools.find(toolName) == m_tools.end() )
    {
      missingTools.push_back(toolName);
    }
  }

  if( !missingTools.empty() )
  {
    cerr << "WARNING: The following expected tools were not registered:" << endl;
    for( const std::string &toolName : missingTools )
    {
      cerr << "  - " << toolName << endl;
    }
  }

  // Check for unexpected tools (tools in registry but not in expected list)
  std::vector<std::string> unexpectedTools;
  for( const auto &[toolName, tool] : m_tools )
  {
    // Skip invoke_ tools as they are dynamically created
    if( toolName.find("invoke_") == 0 )
      continue;

    if( toolName.find("deepresearch_skill_") == 0 )
      continue;

    bool found = false;
    for( const std::string &expected : expectedTools )
    {
      if( toolName == expected )
      {
        found = true;
        break;
      }
    }
    if( !found )
    {
      unexpectedTools.push_back(toolName);
    }
  }

  if( !unexpectedTools.empty() )
  {
    cerr << "WARNING: The following unexpected tools were registered:" << endl;
    for( const std::string &toolName : unexpectedTools )
    {
      cerr << "  - " << toolName << endl;
    }
  }

  //cout << "Registered " << m_tools.size() << " default tools" << endl;
}

const std::map<std::string, SharedTool>& ToolRegistry::getTools() const {
  return m_tools;
}

const SharedTool* ToolRegistry::getTool(const std::string& name) const {
  auto it = m_tools.find(name);
  return (it != m_tools.end()) ? &it->second : nullptr;
}

std::map<std::string, SharedTool> ToolRegistry::getToolsForAgent( AgentType agentType ) const
{
  std::map<std::string, SharedTool> filteredTools;

  for( const auto &[toolName, tool] : m_tools )
  {
    if( (agentType == AgentType::DeepResearch) && (toolName.find("invoke_") == 0) )
      continue;

    if( (agentType == AgentType::DeepResearch) && tool.availableForAgents.empty() )
      continue;

    // If availableForAgents is empty, tool is available to all agents
    if( tool.availableForAgents.empty() )
    {
      filteredTools[toolName] = tool;
      continue;
    }

    // Check if this agent is in the availableForAgents list
    const bool agentAllowed = std::find( tool.availableForAgents.begin(),
                                         tool.availableForAgents.end(),
                                         agentType ) != tool.availableForAgents.end();

    if( agentAllowed )
      filteredTools[toolName] = tool;
  }//for( loop over all tools )

  return filteredTools;
}//getToolsForAgent(...)


std::string ToolRegistry::getDescriptionForAgent( const std::string &toolName, AgentType agentType ) const
{
  const SharedTool * const tool = getTool(toolName);
  if( !tool )
    return "";

  // Check if there's a role-specific description for this agent
  const auto iter = tool->roleDescriptions.find( agentType );
  if( iter != tool->roleDescriptions.end() )
    return iter->second;

  // Fall back to default description
  return tool->description;
}//getDescriptionForAgent(...)


nlohmann::json ToolRegistry::executeTool(const std::string& toolName,
                                       const nlohmann::json& parameters,
                                       InterSpec* interspec,
                                       shared_ptr<LlmInteraction> conversation,
                                       LlmConversationHistory* history) const
{
  const SharedTool* tool = getTool(toolName);
  if (!tool) {
    throw std::runtime_error("Tool not found: " + toolName + ". Perhaps something got garbled somewhere and retrying would help." );
  }

  try
  {
#if( !defined(NDEBUG) && !BUILD_AS_UNIT_TEST_SUITE )
    cout << "Executing tool: " << toolName << " with params: " << parameters.dump() << endl;
#endif

    json result = tool->executor(parameters, interspec, conversation, history);
    
#if( !defined(NDEBUG) && !BUILD_AS_UNIT_TEST_SUITE )
    std::string resultStr = result.dump();
    if( resultStr.length() > 100 )
      resultStr = resultStr.substr(0, 100) + "...";
    cout << "Tool result: " << resultStr << endl;
#endif
    
    return result;
  }catch( const std::exception &e )
  {
    const string err_msg = e.what();
    throw std::runtime_error("Tool execution failed for " + toolName + ": " + err_msg);
  }
}

void ToolRegistry::clearTools() {
  m_tools.clear();
}

// Implementation of specific tool functions
json ToolRegistry::executePeakDetection(const json& params, InterSpec* interspec) {
  if( !interspec )
    throw std::runtime_error("No InterSpec session available");
  
  // Parse parameters into DetectedPeaksOptions
  AnalystChecks::DetectedPeaksOptions options;
  try
  {
    from_json(params, options);
  }catch( std::exception &e )
  {
    cerr << "executePeakDetection: Failed to parse params: " << params.dump() << endl;
    throw;
  }

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
  
  UndoRedoManager::PeakModelChange undo_sentry;
  
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


nlohmann::json ToolRegistry::executeGetIdentifiedSources( const nlohmann::json &params, InterSpec *interspec )
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Parse spectrum type parameter
  const SpecUtils::SpectrumType specType = parse_spectrum_type( params );

  // Call the AnalystChecks function to get the sources
  const std::vector<std::string> sources = AnalystChecks::get_identified_sources( specType, interspec );

  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, sources );

  return result_json;
}//nlohmann::json ToolRegistry::executeGetIdentifiedSources(const nlohmann::json& params, InterSpec* interspec)


nlohmann::json ToolRegistry::executeGetUnidentifiedDetectedPeaks( const nlohmann::json &params, InterSpec *interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");

  // Parse max_results parameter (default: 5)
  size_t max_results = 5;
  if( params.contains("max_results") && !params["max_results"].is_null() )
  {
    const int max_val = params["max_results"].get<int>();
    if( max_val < 1 )
      throw std::runtime_error("max_results must be at least 1");
    max_results = static_cast<size_t>( max_val );
  }

  // Parse verbose parameter (default: false)
  const bool verbose = get_boolean( params, "verbose", false );

  shared_ptr<const SpecUtils::Measurement> meas = interspec->displayedHistogram(SpecUtils::SpectrumType::Foreground);
  shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram(SpecUtils::SpectrumType::Background);

  AnalystChecks::DetectedPeaksOptions options;
  options.specType = SpecUtils::SpectrumType::Foreground;
  options.nonBackgroundPeaksOnly = !!background;
  AnalystChecks::DetectedPeakStatus result = AnalystChecks::detected_peaks(options, interspec);

  vector<shared_ptr<const PeakDef>> unidentified_peaks;
  for( const shared_ptr<const PeakDef> &peak : result.peaks )
  {
    if( !peak->hasSourceGammaAssigned()
       && (peak->userLabel().empty() || SpecUtils::icontains(peak->userLabel(), "unknown") ) )
    {
      unidentified_peaks.push_back( peak );
    }
  }//for( const shared_ptr<const PeakDef> &peak : result.peaks )

  // Sort peaks by amplitude (largest first)
  std::sort( unidentified_peaks.begin(), unidentified_peaks.end(),
    []( const shared_ptr<const PeakDef> &lhs, const shared_ptr<const PeakDef> &rhs ) -> bool {
      return lhs->amplitude() > rhs->amplitude();
    }
  );

  // Apply max_results limit
  if( unidentified_peaks.size() > max_results )
    unidentified_peaks.resize( max_results );

  result.peaks = unidentified_peaks;

  json result_json;

  if( verbose )
  {
    // Verbose mode: return full peak description with all aspects
    to_json( result_json, result, meas );
  }
  else
  {
    // Non-verbose mode: return simplified JSON with only essential fields
    json peaks_array = json::array();
    for( const shared_ptr<const PeakDef> &peak : unidentified_peaks )
    {
      json peak_json = peak_to_simplified_json( peak, meas, false );
      if( !peak_json.empty() )
        peaks_array.push_back( peak_json );
    }
    result_json["peaks"] = peaks_array;
  }

  // Add elevation metrics if background is available and peaks passed the elevation filter
  if( background && options.nonBackgroundPeaksOnly )
  {
    // Get background peaks
    shared_ptr<SpecMeas> background_meas = interspec->measurment(SpecUtils::SpectrumType::Background);
    if( background_meas )
    {
      const set<int> bg_sample_nums = interspec->displayedSamples(SpecUtils::SpectrumType::Background);
      shared_ptr<const deque<shared_ptr<const PeakDef>>> bg_peaks = background_meas->automatedSearchPeaks(bg_sample_nums);
      
      if( bg_peaks && !bg_peaks->empty() )
      {
        const double fg_live_time = (meas && (meas->live_time() > 0.0f)) ? meas->live_time() : 1.0f;
        const double bg_live_time = (background->live_time() > 0.0f) ? background->live_time() : 1.0f;
        
        // Lambda to find matching background peak and calculate elevation metrics
        auto calculate_elevation = [&]( const shared_ptr<const PeakDef> &fg_peak ) 
          -> std::optional<std::pair<double, double>>
        {
          if( !fg_peak )
            return std::nullopt;
          
          // Find the closest matching background peak (same logic as detected_peaks(...) in AnalystChecks.cpp)
          shared_ptr<const PeakDef> closest_bg_peak;
          double smallest_energy_diff = std::numeric_limits<double>::infinity();
          
          for( const shared_ptr<const PeakDef> &bg_peak : *bg_peaks )
          {
            if( !bg_peak )
              continue;
            
            const double energy_diff = fabs(fg_peak->mean() - bg_peak->mean());
            const double avg_fwhm = 0.75 * (fg_peak->fwhm() + bg_peak->fwhm()) / 2.0;
            
            if( (energy_diff < avg_fwhm) && (energy_diff < smallest_energy_diff) )
            {
              closest_bg_peak = bg_peak;
              smallest_energy_diff = energy_diff;
            }
          }
          
          if( !closest_bg_peak )
            return std::nullopt;
          
          // Calculate CPS for both peaks
          const double fg_cps = fg_peak->amplitude() / fg_live_time;
          const double bg_cps = closest_bg_peak->amplitude() / bg_live_time;
          
          if( bg_cps <= 0.0 )
            return std::nullopt;
          
          // Calculate percent elevation: ((fg_cps / bg_cps) - 1.0) * 100.0
          const double percent_elevation = ((fg_cps / bg_cps) - 1.0) * 100.0;
          
          // Calculate sigma elevation if uncertainties are available
          double sigma_elevation = 0.0;
          const double fg_amp_uncert = fg_peak->amplitudeUncert();
          const double bg_amp_uncert = closest_bg_peak->amplitudeUncert();
          
          if( (fg_amp_uncert > 0.0) && (bg_amp_uncert > 0.0) )
          {
            const double fg_cps_uncert = fg_amp_uncert / fg_live_time;
            const double bg_cps_uncert = bg_amp_uncert / bg_live_time;
            const double combined_uncert = sqrt(fg_cps_uncert*fg_cps_uncert + bg_cps_uncert*bg_cps_uncert);
            
            if( combined_uncert > 0.0 )
              sigma_elevation = (fg_cps - bg_cps) / combined_uncert;
          }
          
          return std::make_pair(percent_elevation, sigma_elevation);
        };//calculate_elevation lambda
        
        if( verbose )
        {
          // Verbose mode: navigate ROI structure
          assert( result_json.contains("rois") && result_json["rois"].is_array() );
          
          if( result_json.contains("rois") && result_json["rois"].is_array() )
          {
            json &rois_array = result_json["rois"];
            
            for( json &roi : rois_array )
            {
              assert( roi.is_object() && roi.contains("peaks") && roi["peaks"].is_array() );
              
              if( !roi.is_object() || !roi.contains("peaks") || !roi["peaks"].is_array() )
                continue;
              
              json &peaks_in_roi = roi["peaks"];
              
              for( json &peak_json : peaks_in_roi )
              {
                assert( peak_json.is_object() && peak_json.contains("energy") );
                if( !peak_json.is_object() || !peak_json.contains("energy") )
                  continue;
                
                const double peak_energy = peak_json["energy"].get<double>();
                
                // Find matching peak in unidentified_peaks by energy
                for( const shared_ptr<const PeakDef> &fg_peak : unidentified_peaks )
                {
                  if( !fg_peak || (fabs(fg_peak->mean() - peak_energy) > 0.01) )
                    continue;
                  
                  const std::optional<std::pair<double, double>> elevation = calculate_elevation(fg_peak);
                  
                  if( elevation )
                  {
                    peak_json["elevatedOverBackgroundPeakPercent"] = rount_to_hundredth(elevation->first);
                    
                    if( elevation->second > 0.0 )
                      peak_json["elevatedOverBackgroundPeakNumSigma"] = rount_to_hundredth(elevation->second);
                  }
                  
                  break;
                }//for( const shared_ptr<const PeakDef> &fg_peak : unidentified_peaks )
              }//for( json &peak_json : peaks_in_roi )
            }//for( json &roi : rois_array )
          }//if( result_json.contains("rois") && result_json["rois"].is_array() )
        }else
        {
          // Non-verbose mode: direct peaks array
          assert( result_json.contains("peaks") && result_json["peaks"].is_array() );
          
          if( result_json.contains("peaks") && result_json["peaks"].is_array() )
          {
            json &peaks_array = result_json["peaks"];
            
            // The peaks in the JSON correspond 1-to-1 with unidentified_peaks
            assert( peaks_array.size() <= unidentified_peaks.size() );
            
            for( size_t i = 0; i < peaks_array.size() && i < unidentified_peaks.size(); ++i )
            {
              json &peak_json = peaks_array[i];
              assert( peak_json.is_object() );
              if( !peak_json.is_object() )
                continue;
              
              const std::optional<std::pair<double, double>> elevation = calculate_elevation(unidentified_peaks[i]);
              
              if( elevation )
              {
                peak_json["elevatedOverBackgroundPeakPercent"] = rount_to_hundredth(elevation->first);
                
                if( elevation->second > 0.0 )
                  peak_json["elevatedOverBackgroundPeakNumSigma"] = rount_to_hundredth(elevation->second);
              }
            }//for( size_t i = 0; i < peaks_array.size() && i < unidentified_peaks.size(); ++i )
          }//if( result_json.contains("peaks") && result_json["peaks"].is_array() )
        }//if( verbose ) / else
      }//if( bg_peaks && !bg_peaks->empty() )
    }//if( background_meas )
  }//if( background && options.nonBackgroundPeaksOnly )

  return result_json;
}//nlohmann::json ToolRegistry::executeGetUnidentifiedDetectedPeaks( const nlohmann::json &params, InterSpec *interspec )
  

nlohmann::json ToolRegistry::executeEditAnalysisPeak( const nlohmann::json &params, InterSpec *interspec )
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Parse parameters
  AnalystChecks::EditAnalysisPeakOptions options;
  from_json( params, options );

  UndoRedoManager::PeakModelChange undo_sentry;

  // Execute the edit operation
  const AnalystChecks::EditAnalysisPeakStatus result = AnalystChecks::edit_analysis_peak( options, interspec );

  // Get the spectrum for JSON conversion
  shared_ptr<const SpecUtils::Measurement> meas;
  if( interspec )
    meas = interspec->displayedHistogram( options.specType );

  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result, meas );


  return result_json;
}//nlohmann::json ToolRegistry::executeEditAnalysisPeak( const nlohmann::json &params, InterSpec *interspec )


nlohmann::json ToolRegistry::executeEscapePeakCheck( const nlohmann::json &params, InterSpec *interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");

  // Parse parameters
  AnalystChecks::EscapePeakCheckOptions options;
  from_json( params, options );

  // Execute the escape peak check
  const AnalystChecks::EscapePeakCheckStatus result = AnalystChecks::escape_peak_check( options, interspec );

  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result );

  return result_json;
}//nlohmann::json ToolRegistry::executeEscapePeakCheck( const nlohmann::json &params, InterSpec *interspec )


nlohmann::json ToolRegistry::executeSumPeakCheck( const nlohmann::json &params, InterSpec *interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");

  // Parse parameters
  AnalystChecks::SumPeakCheckOptions options;
  from_json( params, options );

  // Execute the sum peak check
  const AnalystChecks::SumPeakCheckStatus result = AnalystChecks::sum_peak_check( options, interspec );

  // Convert the result to JSON and return
  json result_json;
  to_json( result_json, result );

  return result_json;
}//nlohmann::json ToolRegistry::executeSumPeakCheck( const nlohmann::json &params, InterSpec *interspec )


json ToolRegistry::executeGetSpectrumInfo(const json& params, InterSpec* interspec) {
  if (!interspec) {
    throw std::runtime_error("No InterSpec session available");
  }
  
  string specTypeStr = params.at( find_case_insensitive_key("specType", params) ).get<string>();
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
  if( spectrum->live_time() > 0 )
    result["liveTime"] = spectrum->live_time();
  else
    result["liveTime"] = "N/A";
  
  if( spectrum->real_time() > 0.0 )
    result["realTime"] = spectrum->real_time();
  else
    result["realTime"] = "N/A";
  
  if( (spectrum->live_time() > 0.0) && (spectrum->real_time() > 0.0) )
    result["deadTimeFraction"] = (spectrum->real_time() - spectrum->live_time()) / spectrum->real_time();
  else
    result["deadTimeFraction"] = "N/A";
  
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
    result["neutronCPS_uncertainty"] = sqrt( std::max(2.0, spectrum->neutron_counts_sum()) ) / spectrum->neutron_live_time();
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
  const string source = params.at( find_case_insensitive_key("source", params) ).get<string>();
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
  
  nlohmann::json result;
  result["loaded_spectra"] = json(loadedSpectra);
  
  return result;
}
  
nlohmann::json ToolRegistry::executeGetNuclidesWithCharacteristicsInEnergyRange( const nlohmann::json& params, InterSpec* interspec )
{
  if( !interspec )
    throw std::runtime_error("No InterSpec session available.");
  
  vector<variant<const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *>> result;
  const string lower_energy_key = find_case_insensitive_key( "lowerEnergy", params );
  const string upper_energy_key = find_case_insensitive_key( "upperEnergy", params );
  const string energy_key = find_case_insensitive_key( "energy", params );
  if( params.contains(lower_energy_key) && params.contains(upper_energy_key) )
  {
    const double lower_energy = get_number( params, lower_energy_key );
    const double upper_energy = get_number( params, upper_energy_key );
    result = AnalystChecks::get_nuclides_with_characteristics_in_energy_range( lower_energy, upper_energy, interspec );
  }else if( params.contains(energy_key) )
  {
    const double energy = get_number( params, energy_key );
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
  const string nuclide = params.at( find_case_insensitive_key("source", params) ).get<string>();

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
  const string nuclide = params.at( find_case_insensitive_key("source", params) ).get<string>();
  
  const shared_ptr<const MoreNuclideInfo::MoreNucInfoDb> info_db = MoreNuclideInfo::MoreNucInfoDb::instance();
  if( !info_db )
    throw std::runtime_error("No MoreNucInfoDb instance available");
    
  const MoreNuclideInfo::NucInfo *nuc_info = info_db->info(nuclide);
  if( !nuc_info || nuc_info->m_notes.empty() )
    throw runtime_error( "No info for " + nuclide );
    
  string notes = nuc_info->m_notes;
  SpecUtils::ireplace_all( notes, "\r\n", "\n" );
  
  return json{{"source", nuclide}, {"analystNotes", notes}};
}//nlohmann::json executeGetSourceAnalystNotes(const nlohmann::json& params, InterSpec* interspec)

  
nlohmann::json ToolRegistry::executeGetSourceInfo(const nlohmann::json& params, InterSpec* interspec )
{
  const string nuclide = params.at( find_case_insensitive_key("source", params) ).get<string>();
  
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
    string xray_str = nuclide;
    SpecUtils::ireplace_all(xray_str, "x-ray", "");
    SpecUtils::ireplace_all(xray_str, "xray", "");
    SpecUtils::trim(xray_str);
    
    el = db->element( xray_str );
    if( el )
    {
      result["type"] = "x-ray";
      result["symbol"] = el->symbol;
      result["source"] = el->symbol;
      result["name"] = el->name;
      result["atomicNumber"] = static_cast<int>(el->atomicNumber);
      result["atomicMass"] = el->atomicMass();

      if( xray_str.size() != nuclide.size() )
        result["note"] = "When specifying a x-ray source, only use the atomic symbol (ex Fe, Pb, Cd, etc), and dont include 'x-ray', 'xray', etc.";
      
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
    
    try
    {
      std::vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      if( rctn_db )
        rctn_db->gammas( nuclide, possible_rctns );
      
      if( !possible_rctns.empty() && possible_rctns[0].reaction ) // Take the first reaction found
        rctn = possible_rctns[0].reaction;
    }catch( std::exception & )
    {
      // Happens for example when there isnt an opening and closing paranthesis in the string
    }
    
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
  const string nuclide = params.at( find_case_insensitive_key("nuclide", params) ).get<string>();

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

  // Parse basic parameters
  const double energy = get_number( params, find_case_insensitive_key("energy", params) );
  const double detection_probability = get_number( params, find_case_insensitive_key("detectionProbability", params), 0.95 );
  const float additional_uncertainty = static_cast<float>( get_number( params, find_case_insensitive_key("additionalUncertainty", params), 0.0 ) );
  const bool assert_background_spectrum = get_boolean( params, find_case_insensitive_key("assertBackgroundSpectrum", params), false );

  // Parse optional nuclide, distance, age, and shielding parameters
  const SandiaDecay::Nuclide *nuclide = nullptr;
  double distance = -1.0;
  double age = -1.0;
  double branch_ratio = 0.0;
  double shield_transmission = 1.0;
  bool has_nuclide = false;
  bool has_distance = false;
  bool has_shielding = false;

  const string nuclide_key = find_case_insensitive_key( "nuclide", params );
  if( params.contains(nuclide_key) && params[nuclide_key].is_string() )
  {
    string nuclide_str = params[nuclide_key].get<string>();
    SpecUtils::trim( nuclide_str );
    if( !nuclide_str.empty() )
    {
      const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
      if( !db )
        throw runtime_error( "Could not initialize nuclide DecayDataBase." );
      
      nuclide = db->nuclide( nuclide_str );
      if( !nuclide )
        throw runtime_error( "Could not find nuclide '" + nuclide_str + "' in decay database." );
      
      has_nuclide = true;
    }
  }

  const string distance_key = find_case_insensitive_key( "distance", params );
  if( params.contains(distance_key) && params[distance_key].is_string() )
  {
    string distance_str = params[distance_key].get<string>();
    SpecUtils::trim( distance_str );
    if( !distance_str.empty() )
    {
      distance = PhysicalUnits::stringToDistance( distance_str );
      if( distance < 0.0 )
        throw runtime_error( "Distance must be non-negative, got '" + distance_str + "'." );
      has_distance = true;
    }
  }

  // Get expected FWHM for the energy (needed for ROI width and branch ratio calculation)
  float fwhm = -1.0;
  shared_ptr<SpecMeas> meas = interspec->measurment(SpecUtils::SpectrumType::Foreground);
  
  // If nuclide or distance is provided, require a valid detector
  if( has_nuclide || has_distance )
  {
    if( !meas || !meas->detector() || !meas->detector()->isValid() )
      throw runtime_error( "A valid detector efficiency function is required when nuclide or distance is specified." );
  }
  
  // Try to get FWHM using AnalystChecks::get_expected_fwhm first
  try
  {
    fwhm = AnalystChecks::get_expected_fwhm( energy, interspec );
  }
  catch( std::exception & )
  {
    // Fallback FWHM estimation if AnalystChecks::get_expected_fwhm fails
    if (meas && meas->detector() && meas->detector()->hasResolutionInfo())
      fwhm = meas->detector()->peakResolutionFWHM(static_cast<float>(energy));
    
    if (fwhm <= 0.0) {
      const bool isHPGe = PeakFitUtils::is_likely_high_res(interspec);
      const vector<float> pars = isHPGe ? vector<float>{1.54f, 0.264f, 0.33f} : vector<float>{-6.5f, 7.5f, 0.55f};
      fwhm = DetectorPeakResponse::peakResolutionFWHM(energy, DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn, pars);
    }
  }

  if (fwhm <= 0.0)
    throw std::runtime_error("Could not determine FWHM for energy " + std::to_string(energy));

  if( has_nuclide )
  {
    // Parse age if provided
    const string age_key = find_case_insensitive_key( "age", params );
    if( params.contains(age_key) && params[age_key].is_string() )
    {
      string age_str = params[age_key].get<string>();
      SpecUtils::trim( age_str );
      if( !age_str.empty() )
      {
        try
        {
          age = PhysicalUnits::stringToTimeDuration( age_str );
          if( age < 0.0 )
            throw runtime_error( "Age must be non-negative, got '" + age_str + "'." );
        }catch( std::exception &e )
        {
          throw runtime_error( "Could not interpret '" + age_str + "' as a time duration for nuclide age: " + e.what() );
        }
      }
    }
    
    // If age not provided, use default age for the nuclide
    if( age < 0.0 )
      age = PeakDef::defaultDecayTime( nuclide );

    // Calculate branch ratio for the specified energy
    const double dummy_activity = 0.001 * SandiaDecay::curie;
    SandiaDecay::NuclideMixture mixture;
    mixture.addAgedNuclideByActivity( nuclide, dummy_activity, age );
    
    const vector<SandiaDecay::EnergyRatePair> photons = mixture.photons( 0.0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
    
    // Find the largest amplitude line within 1.0 keV of the specified energy
    const double energy_tolerance = 1.0; // keV
    const SandiaDecay::EnergyRatePair *largest_amplitude_erp = nullptr;
    double largest_amplitude = 0.0;
    
    for( const SandiaDecay::EnergyRatePair &erp : photons )
    {
      if( fabs(erp.energy - energy) < energy_tolerance )
      {
        const double amplitude = erp.numPerSecond / dummy_activity;
        if( amplitude > largest_amplitude )
        {
          largest_amplitude = amplitude;
          largest_amplitude_erp = &erp;
        }
      }
    }
    
    if( !largest_amplitude_erp )
      throw runtime_error( "Could not find a gamma or x-ray within " + to_string(energy_tolerance) + " keV of energy " + to_string(energy) + " keV for nuclide " + nuclide->symbol + "." );
    
    // Calculate peak_sigma from FWHM (FWHM = 2.355 * sigma)
    const double peak_sigma = static_cast<double>(fwhm) / 2.355;
    const double peak_window = 1.25 * peak_sigma;
    const double largest_energy = largest_amplitude_erp->energy;
    
    // Sum all gammas within 1.25*peak_sigma of the largest amplitude gamma line
    for( const SandiaDecay::EnergyRatePair &erp : photons )
    {
      if( fabs(erp.energy - largest_energy) < peak_window )
      {
        branch_ratio += erp.numPerSecond / dummy_activity;
      }
    }
    
    if( branch_ratio <= 0.0 )
      throw runtime_error( "Could not calculate branch ratio for energy " + to_string(energy) + " keV for nuclide " + nuclide->symbol + "." );
  }

  // Parse shielding if provided
  const string shielding_key = find_case_insensitive_key( "shielding", params );
  if( params.contains(shielding_key) && params[shielding_key].is_object() )
  {
    has_shielding = true;
    const json &shielding_json = params[shielding_key];
    
    const bool has_ad = shielding_json.contains("AD");
    const bool has_an = shielding_json.contains("AN");
    const bool has_material = shielding_json.contains("Material");
    const bool has_thickness = shielding_json.contains("Thickness");
    
    if( has_ad && has_an && !has_material && !has_thickness )
    {
      // Generic shielding (AD/AN format)
      const double areal_density = get_double( shielding_json["AD"] ) * (PhysicalUnits::g / PhysicalUnits::cm2);
      const double atomic_number = get_double( shielding_json["AN"] );
      
      if( areal_density < 0.0 )
        throw runtime_error( "Areal density (AD) must be non-negative." );
      
      if( atomic_number <= 0.0 || atomic_number > 100.0 )
        throw runtime_error( "Atomic number (AN) must be between 1 and 100." );
      
      const float energy_float = static_cast<float>( energy * PhysicalUnits::keV );
      const double mu = GammaInteractionCalc::transmition_coefficient_generic(
        static_cast<float>(atomic_number),
        static_cast<float>(areal_density),
        energy_float
      );
      shield_transmission = exp( -mu );
    }
    else if( has_material && has_thickness && !has_ad && !has_an )
    {
      // Material shielding
      const string material_name = shielding_json["Material"].get<string>();
      const string thickness_str = shielding_json["Thickness"].get<string>();
      
      MaterialDB *materialDB = interspec->materialDataBase();
      if( !materialDB )
        throw runtime_error( "Material database not available." );
      
      const Material *material = materialDB->material( material_name );
      if( !material )
        throw runtime_error( "Material '" + material_name + "' not found in material database." );
      
      const double thickness = PhysicalUnits::stringToDistance( thickness_str );
      if( thickness < 0.0 )
        throw runtime_error( "Thickness must be non-negative for material '" + material_name + "', got '" + thickness_str + "'." );
      
      const float energy_float = static_cast<float>( energy * PhysicalUnits::keV );
      const double mu = GammaInteractionCalc::transmition_coefficient_material(
        material,
        energy_float,
        static_cast<float>(thickness)
      );
      shield_transmission = exp( -mu );
    }
    else
    {
      throw runtime_error( "Shielding must have either (AD and AN) or (Material and Thickness), not a mix." );
    }
  }

  // Set up CurrieMdaInput
  DetectionLimitCalc::CurrieMdaInput input;
  input.spectrum = spectrum;
  input.gamma_energy = static_cast<float>(energy);
  
  // Set ROI to be ±1.25 FWHM around the energy (recommended by ISO 11929:2010)
  const float roi_half_width = 1.25f * fwhm;
  input.roi_lower_energy = static_cast<float>(energy) - roi_half_width;
  input.roi_upper_energy = static_cast<float>(energy) + roi_half_width;
  
  // Handle side channels based on assertBackgroundSpectrum
  if( assert_background_spectrum )
  {
    input.num_lower_side_channels = 0;
    input.num_upper_side_channels = 0;
  }
  else
  {
    // Use fixed values for side channels (typical values)
    input.num_lower_side_channels = 4;
    input.num_upper_side_channels = 4;
  }
  
  // Use the parsed parameters
  input.detection_probability = detection_probability;
  input.additional_uncertainty = additional_uncertainty;

  // Call the DetectionLimitCalc function to perform the calculation
  const DetectionLimitCalc::CurrieMdaResult result = DetectionLimitCalc::currie_mda_calc(input);

  // Convert the result to JSON
  json result_json;
  to_json(result_json, result);

  // If nuclide and distance are provided, calculate and add activity information
  if( has_nuclide && has_distance )
  {
    shared_ptr<const DetectorPeakResponse> drf = meas ? meas->detector() : nullptr;
    if( drf && !drf->isValid() )
      drf.reset();
    
    if( !drf )
      throw runtime_error( "Distance specified but no detector efficiency function is currently loaded." );
    
    const bool fixed_geom = drf->isFixedGeometry();
    if( fixed_geom && has_distance )
      throw runtime_error( "Distance cannot be specified when detector efficiency function is for fixed geometry." );
    
    const float energy_float = static_cast<float>(energy);
    const double det_eff = fixed_geom ? drf->intrinsicEfficiency(energy_float)
                                      : drf->efficiency(energy_float, distance);
    
    const float live_time = spectrum->live_time();
    const double air_transmission = (distance > 0.0) 
      ? exp( -GammaInteractionCalc::transmission_coefficient_air( energy_float, static_cast<float>(distance) ) )
      : 1.0;
    
    const double counts_per_bq_into_4pi = branch_ratio * live_time * shield_transmission;
    const double counts_per_bq_into_4pi_with_air = air_transmission * counts_per_bq_into_4pi;
    const double counts_4pi = fixed_geom ? counts_per_bq_into_4pi : counts_per_bq_into_4pi_with_air;
    const double gammas_per_bq = counts_4pi * det_eff;
    
    if( gammas_per_bq > 0.0 )
    {
      // Add activity-related fields
      result_json["gammasPerBq"] = gammas_per_bq;
      result_json["branchRatio"] = branch_ratio;
      
      if( has_shielding )
        result_json["shieldingTransmission"] = shield_transmission;
      
      if( distance > 0.0 && !fixed_geom )
        result_json["airTransmission"] = air_transmission;
      
      if( drf )
      {
        result_json["detectorIntrinsicEfficiency"] = drf->intrinsicEfficiency(energy_float);
        
        if( distance >= 0.0 && !fixed_geom )
        {
          const double geom_eff = DetectorPeakResponse::fractionalSolidAngle( drf->detectorDiameter(), distance );
          result_json["solidAngleFraction"] = geom_eff;
        }
      }
      
      // Convert counts to activity
      if( result.source_counts > result.decision_threshold )
      {
        // Signal detected - provide observed activity and range
        const double nominal_act = result.source_counts / gammas_per_bq;
        result_json["observedActivity"] = nominal_act;
        
        const double lower_act = result.lower_limit / gammas_per_bq;
        const double upper_act = result.upper_limit / gammas_per_bq;
        result_json["activityRange"] = json{{"lower", lower_act}, {"upper", upper_act}};
      }
      else if( result.upper_limit >= 0.0 )
      {
        // No signal detected - provide upper bound
        const double simple_mda = result.upper_limit / gammas_per_bq;
        result_json["activityUpperBound"] = simple_mda;
      }
      
      // Always provide decision threshold and detection limit in activity
      const double decision_threshold_act = result.decision_threshold / gammas_per_bq;
      result_json["decisionThresholdActivity"] = decision_threshold_act;
      
      const double detection_limit_act = result.detection_limit / gammas_per_bq;
      result_json["detectionLimitActivity"] = detection_limit_act;
    }
  }

  return result_json;
}//nlohmann::json executeCurrieMdaCalc(const nlohmann::json& params, InterSpec* interspec)


nlohmann::json ToolRegistry::executeCalculateDose(const nlohmann::json& params, InterSpec* interspec)
{
  using namespace PhysicalUnits;

  // 1. Validate InterSpec instance
  if( !interspec )
    throw runtime_error( "No InterSpec session available." );

  // 2. Parse required parameters
  const string nuclide_str = params.at( find_case_insensitive_key("nuclide", params) ).get<string>();
  const string distance_str = params.at( find_case_insensitive_key("distance", params) ).get<string>();
  const string activity_str = params.at( find_case_insensitive_key("activity", params) ).get<string>();

  // 3. Parse distance and activity
  const double distance = stringToDistance( distance_str );
  const double activity = stringToActivity( activity_str );

  // 4. Get nuclide from database
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Nuclide database not available." );

  const SandiaDecay::Nuclide * const nuc = db->nuclide( nuclide_str );
  if( !nuc )
    throw runtime_error( "Unknown nuclide: " + nuclide_str );

  // 5. Parse optional age parameter (use default if not provided)
  const string age_key = find_case_insensitive_key( "age", params );
  double age;
  if( !params.contains(age_key) )
  {
    age = PeakDef::defaultDecayTime( nuc );
  }else if( params.contains( age_key ) )
  { 
    if( !params[age_key].is_string() )
      throw runtime_error( "the `age` parameter must be a string." );
    
    const string age_str = params[age_key].get<string>();
    try
    {
      age = stringToTimeDuration( age_str );
    }catch( std::runtime_error & )
    {
      throw runtime_error( "Can not interpret age string ('" + age_str + "') as a time duration." );
    }
  }

  // 6. Create nuclide mixture at specified age
  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, activity, age );

  // 7. Extract photon energies and intensities
  vector<float> energies, intensities;
  for( const auto &i : mix.photons(0) )
  {
    energies.push_back( i.energy );
    intensities.push_back( i.numPerSecond );
  }

  // 8. Parse shielding (if provided)
  float areal_density = 0.0f;
  float atomic_number = 26.0f; // iron (doesn't matter when AD=0)

  // Check if material-based shielding is provided
  const string material_key = find_case_insensitive_key( "material", params );
  const string areal_density_key = find_case_insensitive_key( "arealDensity", params );
  const string thickness_key = find_case_insensitive_key( "thickness", params );
  const string atomic_number_key = find_case_insensitive_key( "atomicNumber", params );
  if( params.contains( material_key ) && !params[material_key].is_null() )
  {
    const string material_name = params[material_key].get<string>();

    // Get material from database
    MaterialDB * const materialDB = interspec->materialDataBase();
    const Material * const mat = materialDB->material( material_name );
    if( !mat )
      throw runtime_error( "Unknown material: " + material_name );

    atomic_number = static_cast<float>( mat->massWeightedAtomicNumber() );

    // Check if arealDensity is specified (alternative to thickness)
    if( params.contains( areal_density_key ) && !params[areal_density_key].is_null() )
    {
      areal_density = static_cast<float>( get_number( params, areal_density_key ) );
      areal_density *= static_cast<float>( gram / cm2 );
    }
    else if( params.contains( thickness_key ) && !params[thickness_key].is_null() )
    {
      // Parse thickness and compute areal density
      const string thickness_str = params[thickness_key].get<string>();
      const double thickness = stringToDistance( thickness_str );
      const double density = mat->density; // in PhysicalUnits (g/cm³)
      areal_density = static_cast<float>( thickness * density );
    }
    else
    {
      throw runtime_error( "When material is specified, either thickness or arealDensity must be provided" );
    }
  }
  // Check if direct areal density/atomic number is provided
  else if( params.contains( areal_density_key ) && !params[areal_density_key].is_null() )
  {
    areal_density = static_cast<float>( get_number( params, areal_density_key ) );
    areal_density *= static_cast<float>( gram / cm2 );
    atomic_number = static_cast<float>( get_number( params, atomic_number_key ) );
  }

  // 9. Get scatter table (using helper function)
  shared_ptr<const GadrasScatterTable> scatter = getDoseCalcScatterTable();
  if( !scatter )
    throw runtime_error( "Failed to load scatter table for dose calculations." );

  // 10. Calculate dose
  const double dose = DoseCalc::gamma_dose_with_shielding(
    energies, intensities, areal_density, atomic_number, static_cast<float>(distance), *scatter );

  // 11. Format results
  const string sv_hr = printToBestEquivalentDoseRateUnits( dose, 5, true );
  const string rem_hr = printToBestEquivalentDoseRateUnits( dose, 5, false );

  // 12. Return JSON
  json result;
  result["success"] = true;
  result["dose_rate_Sv_per_hr"] = (dose / (sievert/hour));
  result["dose_rate_si_str"] = sv_hr;
  result["dose_rate_REM_per_hr"] = (dose / (rem/hour));
  result["dose_rate_REM_str"] = rem_hr;
  result["nuclide"] = nuclide_str;
  result["distance"] = distance_str;
  result["activity"] = activity_str;
  result["age"] = printToBestTimeUnits( age, 4 );

  if( areal_density > 0.0f )
  {
    const double ad_g_cm2 = areal_density * cm2 / gram;
    result["shielding"] = json::object();
    result["shielding"]["arealDensity_g_cm2"] = ad_g_cm2;
    result["shielding"]["atomicNumber"] = atomic_number;
  }

  return result;
}//nlohmann::json executeCalculateDose(const nlohmann::json& params, InterSpec* interspec)


/** Returns a list of energy/intensity pairs associated with a nuclide, elemental fluorescence x-ray, or reaction.
 The energies are expressed in keV, and intensities are given as photons per source Becquerel.
 */
nlohmann::json ToolRegistry::executeGetSourcePhotons(const nlohmann::json& params){

  const string source_key = find_case_insensitive_key( "source", params );
  if( !params.contains(source_key) )
    throw runtime_error( "'source' parameter must be specified." );

  // The nuclide, element, or reaction for which to retrieve decay product data.
  //   Examples: 'U238' (nuclide), 'Pb' (element), 'H(n,g)' (reaction).
  string src = params[source_key];

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Could not initialize nuclide DecayDataBase - can not continue." );

  const SandiaDecay::Nuclide *nuc = nullptr;
  const SandiaDecay::Element *el = nullptr;
  const ReactionGamma::Reaction *rctn = nullptr;

  nuc = db->nuclide(src);

  if( !nuc )
  {
    const bool contained_xray = (SpecUtils::icontains(src, "x-ray")
                                 || SpecUtils::icontains(src, "x ray")
                                 || SpecUtils::icontains(src, "xray")
                                 || SpecUtils::icontains(src, "element")
                                 || SpecUtils::icontains(src, "fluorescence") );
    if( contained_xray )
    {
      SpecUtils::ireplace_all( src, "x-ray", "");
      SpecUtils::ireplace_all( src, "x ray", "");
      SpecUtils::ireplace_all( src, "xray", "");
      SpecUtils::ireplace_all( src, "element", "");
      SpecUtils::ireplace_all( src, "fluorescence", "");
    }
    SpecUtils::trim( src );
    
    el = db->element(src);
  }

  if( !nuc && !el )
  {
    const ReactionGamma * const rctn_db = ReactionGammaServer::database();

    try
    {
      std::vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
      if( rctn_db )
        rctn_db->gammas( src, possible_rctns );
      
      if( !possible_rctns.empty() && possible_rctns[0].reaction ) // Take the first reaction found
        rctn = possible_rctns[0].reaction;
    }catch( std::exception & )
    {
      // Happens for example when there isnt an opening and closing parenthesis
    }
  }//if( !nuc && !el )

  if( !nuc && !el && !rctn )
    throw runtime_error( "Could not interpret '" + src + "' as a nuclide, x-ray, or reaction." );

  const string age_key = find_case_insensitive_key( "age", params );
  
  double age_in_seconds = 0.0;
  const bool has_age = params.contains(age_key);
  if( !nuc && has_age )
    throw runtime_error( "You can only specify 'age' for a nuclide source" );
  
  if( has_age )
  {
    const string &age_str = params[age_key];

    try
    {
      age_in_seconds = PhysicalUnits::stringToTimeDuration(age_str) / PhysicalUnits::second;
      if( age_in_seconds < 0.0 )
        throw runtime_error( "Nuclide age ('" + age_str + "') must be larger than zero." );
    }catch( std::exception &e )
    {
      throw runtime_error( "Could not interpret '" + age_str + "' as a time duration for nuclide age." );
    }
  }else if( nuc )
  {
    age_in_seconds = PeakDef::defaultDecayTime( nuc, nullptr ) / PhysicalUnits::second;
  }
  
  const string max_results_key = find_case_insensitive_key( "maxResults", params );
  
  const int max_results = static_cast<int>( get_number( params, max_results_key, 125.0 ) );
  if( max_results < 1 )
    throw runtime_error( "'maxResults' must be 1 or larger." );

  const string cascade_key = find_case_insensitive_key( "includeCascadeSumEnergies", params );

  const bool cascade = get_boolean( params, cascade_key, false );
  if( !nuc && cascade )
    throw runtime_error( "You can only request cascade decay information for nuclides" );
  
  
  vector<pair<double,double>> result;
  
  //transition, first gamma BR, first gamma energy, second gamma energy, coincidence fraction
  typedef tuple<const SandiaDecay::Transition *, double, float, float, float> coincidence_info_t;
  vector<coincidence_info_t> nuc_gamma_coincidences;
  
  
  if( nuc )
  {
    SandiaDecay::NuclideMixture mix;

    // We will add the nuclide as aged, so this way the gammas will be for 1 bq parent activity
    const double parent_activity = 1.0*SandiaDecay::Bq;
    mix.addAgedNuclideByActivity( nuc, parent_activity, age_in_seconds * SandiaDecay::second );
    const vector<SandiaDecay::EnergyRatePair> energy_rate = mix.photons(0.0);
    result.reserve( energy_rate.size() );
    for( const SandiaDecay::EnergyRatePair &erp : energy_rate )
      result.emplace_back( erp.energy, erp.numPerSecond );
    
    if( cascade )
    {
      const vector<SandiaDecay::NuclideActivityPair> activities = mix.activity( 0.0 );
      
      for( const SandiaDecay::NuclideActivityPair &nap : activities )
      {
        const SandiaDecay::Nuclide *nuclide = nap.nuclide;
        const double activity = nap.activity;
        
        for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
        {
          for( const SandiaDecay::RadParticle &particle : transition->products )
          {
            if( (particle.type != SandiaDecay::GammaParticle) )
              continue;
            
            const double br = activity * particle.intensity * transition->branchRatio / parent_activity;
            
            for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
            {
              const unsigned short int part_ind = particle.coincidences[coinc_index].first;
              const float fraction = particle.coincidences[coinc_index].second;
              assert( part_ind < transition->products.size() );
              if( part_ind < transition->products.size() )
              {
                const SandiaDecay::RadParticle &coinc_part = transition->products[part_ind];
                
                // The BR of second gamma is just for debugging
                const double second_br = activity * coinc_part.intensity * transition->branchRatio / parent_activity;
                
                if( coinc_part.type == SandiaDecay::ProductType::GammaParticle )
                  nuc_gamma_coincidences.emplace_back( transition, br, particle.energy, coinc_part.energy, fraction );
              }//if (part_ind < transition->products.size())
            }//for( loop over coincidences )
          }//for( const SandiaDecay::RadParticle &particle : transition->products )
        }//for( const SandiaDecay::Transition *transition : nuclide->decaysToChildren )
      }//for( const SandiaDecay::NuclideActivityPair &nap : activities )
      
      // Sort by max intensity
      std::sort( begin(nuc_gamma_coincidences), end(nuc_gamma_coincidences), []( const coincidence_info_t &lhs, const coincidence_info_t &rhs ) {
        return ( (std::get<1>(lhs)*std::get<4>(lhs)) > (std::get<1>(rhs)*std::get<4>(rhs)) );
      });
      
      if( static_cast<int>(nuc_gamma_coincidences.size()) > max_results )
        nuc_gamma_coincidences.resize( static_cast<size_t>(max_results) );
    }//if( cascade )
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
  
  nlohmann::json answer = nlohmann::json::object();
  answer["photons"] = std::move(json_array);
  if( nuc )
    answer["photonsDescription"] = "Photons emmitted per Bq of parent activity.";
  else if( el )
    answer["photonsDescription"] = "Fluorescent x-rays, and relative intensities, emmitted by " + el->name + ".";
  else if( rctn )
    answer["photonsDescription"] = "Gammas emmitted by the " + rctn->name() + " reaction.";
    
  
  if( nuc && cascade )
  {
    nlohmann::json cascades_array = nlohmann::json::array();
    
    //double max_coinc_amp = 1.0;
    //if( !nuc_gamma_coincidences.empty() )
    //  max_coinc_amp = std::get<1>(nuc_gamma_coincidences.front()) * std::get<4>(nuc_gamma_coincidences.front());
    //assert( !IsNan(max_coinc_amp) && !IsInf(max_coinc_amp) && (max_coinc_amp > 0.0) );
    
    for( const coincidence_info_t &coinc : nuc_gamma_coincidences )
    {
      const double coinc_amp = std::get<1>(coinc) * std::get<4>(coinc);// / max_coinc_amp;
      const SandiaDecay::Transition * const trans = std::get<0>(coinc);
      string transition;
      if( trans && trans->parent )
        transition = trans->parent->symbol + " -> " + (trans->child ? trans->child->symbol : "various"s);
      
      cascades_array.emplace_back( nlohmann::json{
        {"Intensity", coinc_amp},
        {"Energies", {std::get<2>(coinc), std::get<3>(coinc)}},
        {"Transition", std::move(transition) }
      } );
    }//for( const coincidence_info_t &coinc : nuc_gamma_coincidences )
    
    
    answer["cascadeSums"] = std::move(cascades_array);
    answer["cascadeSumsDescription"] = "Photons emitted at the same time during nuclear decay and may be detected together.  The Intensity gives how often the photons are emitted together - not an actual rate, detection probability, or anything comparible to the photons rate. If you applying attenuation or detection probability to the results, first apply to each photon individually, then multiple those results by intensity to get overall relative probabbility.";
    
    return answer;
  }//if( nuc && cascade )
  

  return answer;
}//nlohmann::json executeGetSourcePhotons(const nlohmann::json& params, InterSpec* interspec)



nlohmann::json ToolRegistry::executeGetAttenuationOfShielding( nlohmann::json params, InterSpec* interspec )
{
  using namespace PhysicalUnits;

  // Normalize fields in case they're stringified JSON (resolve keys case-insensitively first)
  const string energies_key = find_case_insensitive_key( "Energies", params );
  const string shielding_key = find_case_insensitive_key( "Shielding", params );
  normalize_json_field( params, energies_key );
  normalize_json_field( params, shielding_key );

  // Parse the energies array
  if( !params.contains(energies_key) || !params[energies_key].is_array() )
    throw runtime_error( "Energies parameter must be specified as an array." );

  const json& energies_json = params[energies_key];
  vector<double> energies;
  for( const auto& energy_val : energies_json )
  {
    energies.push_back( get_double( energy_val ) * keV );
  }

  if( energies.empty() )
    throw runtime_error( "At least one energy must be specified." );

  // Parse the shielding specification
  if( !params.contains(shielding_key) || !params[shielding_key].is_object() )
    throw runtime_error( "Shielding parameter must be specified as an object." );

  const json& shielding = params[shielding_key];

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

  const std::string material_name = params.at( find_case_insensitive_key("material", params) ).get<std::string>();
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
    "Given an array of source photon energies, returns the fraction of photons that will contribute to a peak in data. The returned answer will include the attenuation factor of shielding (if specified), the fraction of photons making it to the detector (if distance is specified), the detection probability of photons that are incident upon the detector (i.e., the intrinsic detection efficiency) if a detector efficiency function is loaded (use 'detector_efficiency_function_info' tool call with no arguments to see if a detector efficiency function is loaded, and 'available_detector_efficiency_functions' and 'load_detector_efficiency_function' to load an efficiency function), and gives total detection probability.",

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
nlohmann::json ToolRegistry::executePhotopeakDetectionCalc(nlohmann::json params, InterSpec* interspec)
{
  using namespace std;
  using namespace nlohmann;

  // Step 1: Parse and normalize inputs

  // Normalize fields in case they're stringified JSON (resolve keys case-insensitively first)
  const string energies_key = find_case_insensitive_key( "Energies", params );
  const string shielding_key = find_case_insensitive_key( "Shielding", params );
  normalize_json_field( params, energies_key );
  normalize_json_field( params, shielding_key );

  // Parse Energies - allow single number or array
  vector<double> energies;
  if( params.contains(energies_key) )
  {
    if( params[energies_key].is_number() || params[energies_key].is_string() )
    {
      energies.push_back( get_double( params[energies_key] ) );
    }
    else if( params[energies_key].is_array() )
    {
      for( const auto& energy_val : params[energies_key] )
        energies.push_back( get_double( energy_val ) );
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
  if( params.contains(shielding_key) )
  {
    if( params[shielding_key].is_object() )
    {
      shieldings.push_back( params[shielding_key] );
    }
    else if( params[shielding_key].is_array() )
    {
      shieldings = params[shielding_key].get<vector<json>>();
    }
    else
    {
      throw runtime_error( "Shielding must be an object or array of objects" );
    }
  }

  // Parse Distance (optional)
  const string distance_key = find_case_insensitive_key( "Distance", params );
  double distance = 0.0;
  bool has_distance = false;
  if( params.contains(distance_key) && params[distance_key].is_string() )
  {
    string distance_str = params[distance_key].get<string>();
    SpecUtils::trim( distance_str );
    if( !distance_str.empty() )
    {
      distance = PhysicalUnits::stringToDistance( distance_str );
      has_distance = true;
    }//if( !distance_str.empty() )
  }//if( params.contains("Distance") && params["Distance"].is_string() )

  // Parse IncludeAirAttenuation (optional, default false)
  const bool include_air = get_boolean( params, find_case_insensitive_key("IncludeAirAttenuation", params), false );

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
      info.areal_density = get_double( shield_json["AD"] );
      info.atomic_number = get_double( shield_json["AN"] );
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
  

  /** cooresponds to the `available_detector_efficiency_functions` tool call, which has the description:
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
          "description": "The detector efficiency function identifier. Can be a name from available_detector_efficiency_functions, a filesystem path to a detector file or GADRAS directory, or a URI. Required."
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

  const string identifier_key = find_case_insensitive_key( "identifier", params );
  if( !params.contains(identifier_key) || !params[identifier_key].is_string() )
    throw std::runtime_error( "Missing required 'identifier' parameter." );

  const std::string identifier = params[identifier_key].get<std::string>();
  const string detector_name_key = find_case_insensitive_key( "detectorName", params );
  const std::string detectorName = params.contains(detector_name_key) ? params[detector_name_key].get<std::string>() : std::string();
  const string source_hint_key = find_case_insensitive_key( "source", params );
  const std::string sourceHint = params.contains(source_hint_key) ? params[source_hint_key].get<std::string>() : std::string("AnySource");

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
  const string name_key = find_case_insensitive_key( "name", params );
  if( params.contains(name_key) && params[name_key].is_string() && !params[name_key].get<std::string>().empty() )
  {
    const std::string name = params[name_key].get<std::string>();
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
    case DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic:
    case DetectorPeakResponse::EffGeometryType::FarFieldAbsolute:
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
  result["isFixedGeometry"] = (geomType != DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic) && (geomType != DetectorPeakResponse::EffGeometryType::FarFieldAbsolute);

  // Energy range
  result["lowerEnergy"] = drf->lowerEnergy();
  result["upperEnergy"] = drf->upperEnergy();

  // Detector diameter (only for far field)
  if( (geomType == DetectorPeakResponse::EffGeometryType::FarFieldIntrinsic) || (geomType == DetectorPeakResponse::EffGeometryType::FarFieldAbsolute) )
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


nlohmann::json ToolRegistry::executeSearchSourcesByEnergy(nlohmann::json params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Normalize energies field in case it's a stringified JSON array (resolve key case-insensitively first)
  const string energies_key = find_case_insensitive_key( "energies", params );
  normalize_json_field( params, energies_key );

  // Parse energies array (required)
  if( !params.contains(energies_key) || !params[energies_key].is_array() || params[energies_key].empty() )
    throw std::runtime_error( "The 'energies' parameter must be a non-empty array." );

  const json& energies_array = params[energies_key];
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
    if( energy_obj.contains("window") && !energy_obj["window"].is_null() )
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
  const string category_key = find_case_insensitive_key( "source_category", params );
  std::string category_str = "isbe-category-nuc-xray";  // "Nuclides + X-rays"
  if( params.contains(category_key) && params[category_key].is_string() )
    category_str = params[category_key].get<std::string>();

  const string min_half_life_key = find_case_insensitive_key( "min_half_life", params );
  std::string min_half_life_str = "100 m";
  if( params.contains(min_half_life_key) && params[min_half_life_key].is_string() )
    min_half_life_str = params[min_half_life_key].get<std::string>();

  const string min_br_key = find_case_insensitive_key( "min_branching_ratio", params );
  double min_br = 0.0;
  if( params.contains(min_br_key) && !params[min_br_key].is_null() )
    min_br = get_number( params, min_br_key );

  const string max_results_key = find_case_insensitive_key( "max_results", params );
  int max_results = 10;
  if( params.contains(max_results_key) && !params[max_results_key].is_null() )
    max_results = static_cast<int>( get_number( params, max_results_key ) );

  if( max_results <= 0 )
    throw std::runtime_error( "max_results must be positive." );

  const string sort_by_key = find_case_insensitive_key( "sort_by", params );
  std::string sort_by = "ProfileScore";
  if( params.contains(sort_by_key) && params[sort_by_key].is_string() )
    sort_by = params[sort_by_key].get<std::string>();

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




nlohmann::json ToolRegistry::executeSetWorkflowState(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;

  if( !interspec )
    throw runtime_error( "InterSpec instance required for set_workflow_state" );

  // Get parameters
  const string state_key = find_case_insensitive_key( "state", params );
  const string notes_key = find_case_insensitive_key( "notes", params );
  const string new_state = params.contains(state_key) ? params[state_key].get<string>() : string();
  const string notes = params.contains(notes_key) ? params[notes_key].get<string>() : string();

  if( new_state.empty() )
    throw runtime_error( "state parameter is required" );

  json result;
  result["success"] = true;
  result["new_state"] = new_state;

  if( !notes.empty() )
    result["notes"] = notes;

  // Get current LLM interface and conversation
  LlmToolGui *llm_gui = interspec->currentLlmTool();
  if( !llm_gui )
  {
    result["warning"] = "No active LLM interface - state tracking not available";
    return result;
  }

  LlmInterface *interface = llm_gui->llmInterface();
  if( !interface )
  {
    result["warning"] = "No LLM interface found - state tracking not available";
    return result;
  }

  shared_ptr<LlmInteraction> conversation = interface->getCurrentConversation();
  if( !conversation )
  {
    result["warning"] = "No active conversation - state tracking not available";
    return result;
  }

  // Check if conversation has a state machine
  if( !conversation->state_machine )
  {
    result["note"] = "Agent does not have a state machine defined";
    return result;
  }

  AgentStateMachine *sm = conversation->state_machine.get();
  const string previous_state = sm->getCurrentState();
  result["previous_state"] = previous_state;

  // Check if state exists
  if( !sm->hasState(new_state) )
  {
    result["success"] = false;
    result["error"] = "Unknown state: " + new_state;
    return result;
  }

  // Check if transition is valid (soft enforcement)
  if( !sm->canTransitionTo(new_state) )
  {
    result["warning"] = "Unexpected state transition";
    result["expected_transitions"] = sm->getAllowedTransitions();
    result["reason"] = "Transition from " + previous_state +
                       " to " + new_state + " not in allowed transitions";
    result["transition_forced"] = true;
    // Continue anyway - soft enforcement
  }

  // Perform transition
  sm->transitionTo(new_state);

  // Build comprehensive response with narrative context
  const AgentStateMachine::StateDefinition& state_def = sm->getStateDefinition(new_state);

  // Add narrative transition summary
  string transition_summary = "You are now in state **" + new_state + "**";
  if( !previous_state.empty() && previous_state != new_state )
  {
    transition_summary += ", having transitioned from **" + previous_state + "**";
  }
  transition_summary += ".";
  result["transition_summary"] = transition_summary;

  // Include state-specific information
  result["description"] = state_def.description;

  // Enhanced guidance format with clearer structure
  string guidance_text = "\n**Current State Guidance:**\n" + state_def.prompt_guidance;
  result["guidance"] = guidance_text;

  result["is_final_state"] = state_def.is_final;

  // Allowed transitions
  if( !state_def.allowed_transitions.empty() )
  {
    result["allowed_next_states"] = state_def.allowed_transitions;

    string transitions_text = "\n**Allowed Next States:** ";
    for( size_t i = 0; i < state_def.allowed_transitions.size(); ++i )
    {
      if( i > 0 )
        transitions_text += ", ";
      transitions_text += state_def.allowed_transitions[i];
    }
    result["allowed_transitions_text"] = transitions_text;
  }
  else
  {
    result["allowed_next_states"] = json::array();
  }

  // Include suggested tools if any
  if( !state_def.suggested_tools.empty() )
    result["suggested_tools"] = state_def.suggested_tools;

  return result;
}

json ToolRegistry::executeCreatePeakCheckpoint( const json& params,
                                                InterSpec* interspec,
                                                LlmConversationHistory* history )
{
  if( !interspec )
    throw runtime_error( "No InterSpec session available." );

  if( !history )
    throw runtime_error( "Peak checkpoints require a conversation history context." );

  // Get the optional name hint
  string name_hint;
  if( params.contains( "name_hint" ) && params["name_hint"].is_string() )
    name_hint = params["name_hint"].get<string>();

  // Generate a unique checkpoint name
  string checkpoint_name;
  if( name_hint.empty() )
  {
    const auto ms = chrono::duration_cast<chrono::milliseconds>(
      chrono::system_clock::now().time_since_epoch() ).count();
    checkpoint_name = "checkpoint_" + to_string( ms );
  }else
  {
    checkpoint_name = name_hint;
  }

  // Ensure uniqueness by appending a suffix if needed
  {
    const vector<PeakCheckpoint> &existing = history->m_peak_checkpoints;
    bool name_taken = false;
    for( const PeakCheckpoint &cp : existing )
    {
      if( cp.m_checkpoint_name == checkpoint_name )
      {
        name_taken = true;
        break;
      }
    }

    if( name_taken )
    {
      for( int suffix = 2; suffix < 1000; ++suffix )
      {
        const string candidate = checkpoint_name + "_" + to_string( suffix );
        bool found = false;
        for( const PeakCheckpoint &cp : existing )
        {
          if( cp.m_checkpoint_name == candidate )
          {
            found = true;
            break;
          }
        }

        if( !found )
        {
          checkpoint_name = candidate;
          break;
        }
      }//for( int suffix = 2; suffix < 1000; ++suffix )
    }//if( name_taken )
  }

  PeakCheckpoint checkpoint;
  checkpoint.m_checkpoint_name = checkpoint_name;
  checkpoint.m_creation_time = chrono::system_clock::now();

  PeakModel * const peak_model = interspec->peakModel();
  assert( peak_model );
  if( !peak_model )
    throw runtime_error( "No PeakModel available." );

  // Snapshot peaks and energy calibration for each spectrum type
  const SpecUtils::SpectrumType spec_types[] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };

  for( const SpecUtils::SpectrumType type : spec_types )
  {
    const shared_ptr<const SpecUtils::Measurement> meas = interspec->displayedHistogram( type );
    if( !meas )
      continue;

    const shared_ptr<const deque<shared_ptr<const PeakDef>>> current_peaks = peak_model->peaks( type );

    // Deep-copy the peaks deque (the shared_ptr from PeakModel may be mutated later)
    shared_ptr<deque<shared_ptr<const PeakDef>>> peaks_copy;
    if( current_peaks && !current_peaks->empty() )
      peaks_copy = make_shared<deque<shared_ptr<const PeakDef>>>( *current_peaks );

    const shared_ptr<const SpecUtils::EnergyCalibration> cal = meas->energy_calibration();

    switch( type )
    {
      case SpecUtils::SpectrumType::Foreground:
        checkpoint.m_foreground_peaks = peaks_copy;
        checkpoint.m_foreground_cal = cal;
        break;
      case SpecUtils::SpectrumType::Background:
        checkpoint.m_background_peaks = peaks_copy;
        checkpoint.m_background_cal = cal;
        break;
      case SpecUtils::SpectrumType::SecondForeground:
        checkpoint.m_secondary_peaks = peaks_copy;
        checkpoint.m_secondary_cal = cal;
        break;
    }//switch( type )
  }//for( const SpecUtils::SpectrumType type : spec_types )

  history->m_peak_checkpoints.push_back( std::move( checkpoint ) );

  json result = json::object();
  result["checkpoint_name"] = checkpoint_name;
  result["status"] = "created";

  return result;
}//executeCreatePeakCheckpoint(...)


json ToolRegistry::executeRestorePeaksToCheckpoint( const json& params,
                                                     InterSpec* interspec,
                                                     LlmConversationHistory* history )
{
  if( !interspec )
    throw runtime_error( "No InterSpec session available." );

  if( !history )
    throw runtime_error( "Peak checkpoints require a conversation history context." );

  if( !params.contains( "checkpoint_name" ) || !params["checkpoint_name"].is_string() )
    throw runtime_error( "Missing required parameter 'checkpoint_name'." );

  const string checkpoint_name = params["checkpoint_name"].get<string>();

  // Find the checkpoint (search from newest to oldest)
  const PeakCheckpoint *found_cp = nullptr;
  for( auto it = history->m_peak_checkpoints.rbegin();
       it != history->m_peak_checkpoints.rend(); ++it )
  {
    if( it->m_checkpoint_name == checkpoint_name )
    {
      found_cp = &(*it);
      break;
    }
  }

  if( !found_cp )
  {
    string available_names;
    for( const PeakCheckpoint &cp : history->m_peak_checkpoints )
    {
      if( !available_names.empty() )
        available_names += ", ";
      available_names += "'" + cp.m_checkpoint_name + "'";
    }

    throw runtime_error( "Checkpoint '" + checkpoint_name + "' not found."
      + (available_names.empty() ? " No checkpoints exist." : " Available: " + available_names) );
  }

  const SpecUtils::SpectrumType spec_types[] = {
    SpecUtils::SpectrumType::Foreground,
    SpecUtils::SpectrumType::Background,
    SpecUtils::SpectrumType::SecondForeground
  };

  vector<string> restored_types;

  for( const SpecUtils::SpectrumType type : spec_types )
  {
    const shared_ptr<const SpecUtils::Measurement> meas = interspec->displayedHistogram( type );
    if( !meas )
      continue;

    // Get the checkpoint's peaks and calibration for this spectrum type
    shared_ptr<deque<shared_ptr<const PeakDef>>> cp_peaks;
    shared_ptr<const SpecUtils::EnergyCalibration> cp_cal;

    switch( type )
    {
      case SpecUtils::SpectrumType::Foreground:
        cp_peaks = found_cp->m_foreground_peaks;
        cp_cal = found_cp->m_foreground_cal;
        break;
      case SpecUtils::SpectrumType::Background:
        cp_peaks = found_cp->m_background_peaks;
        cp_cal = found_cp->m_background_cal;
        break;
      case SpecUtils::SpectrumType::SecondForeground:
        cp_peaks = found_cp->m_secondary_peaks;
        cp_cal = found_cp->m_secondary_cal;
        break;
    }//switch( type )

    // If checkpoint had no data for this type (spectrum wasn't loaded at checkpoint time), skip
    if( !cp_cal )
      continue;

    shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks_to_set;

    if( !cp_peaks || cp_peaks->empty() )
    {
      // Checkpoint had no peaks for this type - set empty
      peaks_to_set = make_shared<const deque<shared_ptr<const PeakDef>>>();
    }else
    {
      // Check if energy calibration has changed
      const shared_ptr<const SpecUtils::EnergyCalibration> current_cal = meas->energy_calibration();
      assert( current_cal );

      if( !current_cal || (cp_cal == current_cal) || (*cp_cal == *current_cal) )
      {
        // Calibration unchanged - use peaks directly
        peaks_to_set = cp_peaks;
      }else
      {
        // Calibration changed - translate peaks
        const deque<shared_ptr<const PeakDef>> translated =
          EnergyCal::translatePeaksForCalibrationChange( *cp_peaks, cp_cal, current_cal );
        peaks_to_set = make_shared<const deque<shared_ptr<const PeakDef>>>( translated );
      }
    }//if( !cp_peaks || cp_peaks->empty() ) / else

    interspec->setPeaks( type, peaks_to_set );

    switch( type )
    {
      case SpecUtils::SpectrumType::Foreground:      restored_types.push_back( "Foreground" );      break;
      case SpecUtils::SpectrumType::Background:       restored_types.push_back( "Background" );       break;
      case SpecUtils::SpectrumType::SecondForeground: restored_types.push_back( "Secondary" ); break;
    }
  }//for( const SpecUtils::SpectrumType type : spec_types )

  json result = json::object();
  result["checkpoint_name"] = checkpoint_name;
  result["status"] = "restored";
  result["restored_spectrum_types"] = restored_types;

  return result;
}//executeRestorePeaksToCheckpoint(...)


} // namespace LlmTools

#endif // USE_LLM_INTERFACE
