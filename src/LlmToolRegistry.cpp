#include "InterSpec_config.h"
#include "InterSpec/LlmToolRegistry.h"

#if( USE_LLM_INTERFACE )

#include <algorithm>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/LlmInterface.h"
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

#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"

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


  /** Some LLMs will give boolean values as strings, so this function will check types and return the correct answer.

   @param parent The parent JSON object.
   @param name The field name of the boolean to parse.

   An exception will be thrown if cant be converted to boolean.
   */
  bool get_boolean( const json& parent, const string &name )
  {
    if( !parent.contains(name) )
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
    if( !parent.contains(name) )
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
    if( !parent.contains(name) )
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
    const json& nuclideParam = j.at("nuclide");
    if (nuclideParam.is_string()) {
      p.nuclides = {nuclideParam.get<std::string>()};
    } else if (nuclideParam.is_array()) {
      p.nuclides = nuclideParam.get<std::vector<std::string>>();
    } else {
      throw std::runtime_error("Invalid nuclide parameter: must be string or array of strings");
    }

    p.doNotAddPeaksToUserSession = get_boolean( j, "doNotAddPeaksToUserSession", false );
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
    j = json{{"sources", sources}};
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
    tool.executor = [](const json& params, InterSpec* interspec) -> json {
      return executePeakDetection(params, interspec);
      };
    }else if( toolName == "add_analysis_peak" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executePeakFit(params, interspec);
      };
    }else if( toolName == "edit_analysis_peak" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeEditAnalysisPeak(params, interspec);
      };
    }else if( toolName == "escape_peak_check" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeEscapePeakCheck(params, interspec);
      };
    }else if( toolName == "sum_peak_check" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeSumPeakCheck(params, interspec);
      };
    }else if( toolName == "get_analysis_peaks" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetUserPeaks(params, interspec);
      };
    }else if( toolName == "get_identified_sources" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetIdentifiedSources(params, interspec);
      };
    }else if( toolName == "get_unidentified_peaks" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetUnidentifiedDetectedPeaks(params, interspec);
      };
    }else if( toolName == "get_spectrum_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetSpectrumInfo(params, interspec);
      };
    }else if( toolName == "primary_gammas_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetCharacteristicGammasForSource(params);
      };
    }else if( toolName == "sources_with_primary_gammas_in_energy_range" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
      };
    }else if( toolName == "sources_with_primary_gammas_near_energy" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetNuclidesWithCharacteristicsInEnergyRange(params, interspec);
      };
    }
#if( !INCLUDE_NOTES_AND_ASSOCIATED_SRCS_WITH_SRC_INFO )
    else if( toolName == "sources_associated_with_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetAssociatedSources(params);
      };
    }else if( toolName == "analyst_notes_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetSourceAnalystNotes(params);
      };
    }
#endif
    else if( toolName == "source_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetSourceInfo(params, interspec);
      };
    }else if( toolName == "nuclide_decay_chain" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetNuclideDecayChain(params);
      };
    }else if( toolName == "automated_source_id_results" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetAutomatedRiidId(params, interspec);
      };
    }else if( toolName == "loaded_spectra" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetLoadedSpectra(params, interspec);
      };
    }
  /*
   // Not using `add_analysis_peaks_for_source` until we get it working properly...
    else if( toolName == "add_analysis_peaks_for_source" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeFitPeaksForNuclide(params, interspec);
      };
    }
    */
    else if( toolName == "get_counts_in_energy_range" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetCountsInEnergyRange(params, interspec);
      };
    }else if( toolName == "get_expected_fwhm" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetExpectedFwhm(params, interspec);
      };
    }else if( toolName == "currie_mda_calc" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeCurrieMdaCalc(params, interspec);
      };
    }else if( toolName == "source_photons" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetSourcePhotons(params);
      };
    }else if( toolName == "photopeak_detection_efficiency" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executePhotopeakDetectionCalc(params, interspec);
      };
    }else if( toolName == "get_materials" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetMaterials(interspec);
      };
    }else if( toolName == "get_material_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetMaterialInfo(params, interspec);
      };
    }else if( toolName == "available_detector_efficiency_functions" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeAvailableDetectors(params, interspec);
      };
    }else if( toolName == "load_detector_efficiency_function" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeLoadDetectorEfficiency(params, interspec);
      };
    }else if( toolName == "detector_efficiency_function_info" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetDetectorInfo(params, interspec);
      };
    }else if( toolName == "search_sources_by_energy" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeSearchSourcesByEnergy(params, interspec);
      };
    }else if( toolName == "activity_fit" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeActivityFit(params, interspec);
      };
    }else if( toolName == "get_shielding_source_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeGetShieldingSourceConfig(params, interspec);
      };
    }else if( toolName == "modify_shielding_source_config" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeModifyShieldingSourceConfig(params, interspec);
      };
    }else if( toolName == "mark_peaks_for_activity_fit" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeMarkPeaksForActivityFit(params, interspec);
      };
    }else if( toolName == "close_activity_shielding_display" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeCloseActivityShieldingDisplay(params, interspec);
      };
    }else if( toolName == "ask_user_question" )
    {
      tool.executor = [](const json& params, InterSpec* interspec) -> json {
        return executeAskUserQuestion(params, interspec);
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

  // Helper lambda to apply config overrides to a tool
  auto applyToolConfig = [&config](SharedTool &tool) {
    // Look for this tool in the config
    for( const LlmConfig::ToolConfig &toolConfig : config.tools )
    {
      if( toolConfig.name == tool.name )
      {
        // Apply agent restrictions
        if( !toolConfig.availableForAgents.empty() )
          tool.availableForAgents = toolConfig.availableForAgents;

        // Apply default description override
        if( !toolConfig.defaultDescription.empty() )
          tool.description = toolConfig.defaultDescription;

        // Apply role-specific descriptions
        tool.roleDescriptions = toolConfig.roleDescriptions;

        // Apply parameter schema override (already parsed and validated; {} is a valid empty schema)
        if( !toolConfig.parametersSchema.is_null() )
        {
          tool.parameters_schema = toolConfig.parametersSchema;
        }

        cout << "  Applied config overrides for tool: " << tool.name << endl;
        break;
      }
    }
  };

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
            "description": "Background context about the current analysis state"
          },
          "task": {
            "type": "string",
            "description": "The specific task for the sub-agent to perform"
          }
        },
        "required": ["context", "task"]
      })");
    
    // Executor is a placeholder - actual invocation handled in executeToolCalls
    invokeTool.executor = [](const json& params, InterSpec* interspec) -> json {
      assert( 0 );
      throw std::runtime_error( "invoke_* tools should be handled by executeToolCalls, not called directly" );
    };

    // Set which agents can invoke this sub-agent
    // If agent.availableForAgents is empty, default to MainAgent only for backward compatibility
    if( agent.availableForAgents.empty() )
      invokeTool.availableForAgents = {AgentType::MainAgent};
    else
      invokeTool.availableForAgents = agent.availableForAgents;

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
      {
        tool.parameters_schema = toolConfig.parametersSchema;
      }

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

  // NOTE: Hard-coded fallback tool definitions have been removed.
  // If no tools are loaded from XML, we throw an exception above.
  // This ensures that:
  //   1. The XML configuration is the single source of truth for tool definitions
  //   2. Tool descriptions and schemas can be updated without recompilation
  //   3. Missing configuration is caught early rather than silently using outdated hardcoded values

  // Runtime validation: Check that all expected tools are registered
  const std::vector<std::string> expectedTools = {
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
    //"add_analysis_peaks_for_source",
    "get_counts_in_energy_range",
    "get_expected_fwhm",
    "currie_mda_calc",
    "source_photons",
    "photopeak_detection_efficiency",
    "get_materials",
    "get_material_info",
    "available_detector_efficiency_functions",
    "load_detector_efficiency_function",
    "detector_efficiency_function_info",
    "search_sources_by_energy",
    "activity_fit",
    "ask_user_question",
    "close_activity_shielding_display",
    "get_shielding_source_config",
    "mark_peaks_for_activity_fit",
    "modify_shielding_source_config"
  };

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
                                       InterSpec* interspec)  const
{
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
  to_json( result_json, result, meas );

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
    
    try
    {
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
  const double detection_probability = get_number( params, "detectionProbability", 0.95 );
  const float additional_uncertainty = static_cast<float>( get_number( params, "additionalUncertainty", 0.0 ) );

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
  string src = params["Source"];

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

  const int max_results = static_cast<int>( get_number( params, "MaxResults", 125.0 ) );
  if( max_results < 1 )
    throw runtime_error( "'MaxResults' must be 1 or larger." );


  const bool cascade = get_boolean( params, "IncludeCascadeSumEnergies", false );
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
  
  
  if( nuc && cascade )
  {
    nlohmann::json answer = nlohmann::json::object();
    answer["photons"] = std::move(json_array);
    answer["photonsDescription"] = "Photons emmitted per Bq of parent activity.";
    
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
  

  return json_array;
}//nlohmann::json executeGetSourcePhotons(const nlohmann::json& params, InterSpec* interspec)



nlohmann::json ToolRegistry::executeGetAttenuationOfShielding( nlohmann::json params, InterSpec* interspec )
{
  using namespace PhysicalUnits;

  // Normalize fields in case they're stringified JSON
  normalize_json_field( params, "Energies" );
  normalize_json_field( params, "Shielding" );

  // Parse the energies array
  if( !params.contains("Energies") || !params["Energies"].is_array() )
    throw runtime_error( "Energies parameter must be specified as an array." );

  const json& energies_json = params["Energies"];
  vector<double> energies;
  for( const auto& energy_val : energies_json )
  {
    energies.push_back( get_double( energy_val ) * keV );
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

  // Normalize fields in case they're stringified JSON
  normalize_json_field( params, "Energies" );
  normalize_json_field( params, "Shielding" );

  // Parse Energies - allow single number or array
  vector<double> energies;
  if( params.contains("Energies") )
  {
    if( params["Energies"].is_number() || params["Energies"].is_string() )
    {
      energies.push_back( get_double( params["Energies"] ) );
    }
    else if( params["Energies"].is_array() )
    {
      for( const auto& energy_val : params["Energies"] )
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
    string distance_str = params["Distance"].get<string>();
    SpecUtils::trim( distance_str );
    if( !distance_str.empty() )
    {
      distance = PhysicalUnits::stringToDistance( distance_str );
      has_distance = true;
    }//if( !distance_str.empty() )
  }//if( params.contains("Distance") && params["Distance"].is_string() )

  // Parse IncludeAirAttenuation (optional, default false)
  const bool include_air = get_boolean( params, "IncludeAirAttenuation", false );

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


nlohmann::json ToolRegistry::executeSearchSourcesByEnergy(nlohmann::json params, InterSpec* interspec)
{
  if( !interspec )
    throw std::runtime_error( "No InterSpec session available." );

  // Normalize energies field in case it's a stringified JSON array
  normalize_json_field( params, "energies" );

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


// ============================================================================
// Activity/Shielding Fit Tools Implementation
// ============================================================================

namespace {
  // Helper functions for activity/shielding fitting

  /** Parse a distance/length string with units (e.g., "100 cm", "3 ft", "1.5 m")
   * Returns value in PhysicalUnits (mm base unit).
   * If input is numeric without units, assumes cm.
   */
  double parse_distance_string( const string &distance_str )
  {
    string trimmed = distance_str;
    SpecUtils::trim( trimmed );

    // Try to parse with units first
    try
    {
      return PhysicalUnits::stringToDistance( trimmed );
    }catch( std::exception & )
    {
      // If that failed, check if it's a plain number
    }

    // Check if string contains only valid numeric characters
    bool has_invalid_chars = false;
    for( const char c : trimmed )
    {
      if( !std::isdigit(c) && c != '.' && c != 'e' && c != 'E' &&
          c != '+' && c != '-' && c != ' ' )
      {
        has_invalid_chars = true;
        break;
      }
    }

    if( has_invalid_chars )
      throw runtime_error( "Could not parse distance '" + distance_str + "': invalid format" );

    // Parse as plain number and assume cm
    double val = 0.0;
    if( !SpecUtils::parse_double( trimmed.c_str(), trimmed.size(), val ) )
      throw runtime_error( "Could not parse distance '" + distance_str + "' as a number" );

    return val * PhysicalUnits::cm;
  }//parse_distance_string(...)


  /** Parse an activity string with units (e.g., "10 uCi", "1 MBq", "100 Bq")
   * Returns value in PhysicalUnits (becquerels).
   */
  double parse_activity_string( const string &activity_str )
  {
    string trimmed = activity_str;
    SpecUtils::trim( trimmed );

    try
    {
      return PhysicalUnits::stringToActivity( trimmed );
    }catch( std::exception &e )
    {
      throw runtime_error( "Could not parse activity '" + activity_str + "': " + e.what() );
    }
  }//parse_activity_string(...)


  /** Parse an age/time string with units (e.g., "20 years", "6 months", "100 days")
   * Returns value in PhysicalUnits (seconds).
   * If input is numeric, assumes seconds.
   */
  double parse_age_string( const string &age_str )
  {
    string trimmed = age_str;
    SpecUtils::trim( trimmed );

    // Try to parse as plain number first (assume seconds)
    try
    {
      const double val = std::stod( trimmed );
      return val * PhysicalUnits::second;
    }catch(...)
    {
    }

    // Parse with units
    try
    {
      return PhysicalUnits::stringToTimeDuration( trimmed );
    }catch( std::exception &e )
    {
      throw runtime_error( "Could not parse age/time '" + age_str + "': " + e.what() );
    }
  }//parse_age_string(...)


  /** Find a peak by energy with tolerance of max(FWHM, 1 keV).
   * Returns pointer to peak, or nullptr if not found.
   * Throws exception if multiple peaks match (ambiguous).
   */
  std::shared_ptr<const PeakDef> find_peak_by_energy(
    const double energy,
    const std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &peaks,
    const string &context_msg = ""
  )
  {
    if( !peaks || peaks->empty() )
      return nullptr;

    std::shared_ptr<const PeakDef> found_peak;

    for( const auto &peak : *peaks )
    {
      if( !peak )
        continue;

      const double peak_energy = peak->mean();
      const double tolerance = std::max( peak->fwhm(), 1.0 );

      if( fabs(energy - peak_energy) <= tolerance )
      {
        if( found_peak )
        {
          // Multiple peaks match - ambiguous
          string msg = "Multiple peaks found near " + std::to_string(energy) + " keV";
          if( !context_msg.empty() )
            msg += " (" + context_msg + ")";
          msg += ". Found peaks at " + std::to_string(found_peak->mean())
                + " keV and " + std::to_string(peak_energy) + " keV.";
          throw runtime_error( msg );
        }
        found_peak = peak;
      }
    }

    return found_peak;
  }//find_peak_by_energy(...)


  /** Check if age fitting is allowed for a nuclide (equivalent to !PeakDef::ageFitNotAllowed()).
   * Age fitting is not allowed if:
   * - Nuclide is null or decays to stable children
   * - Nuclide reaches equilibrium too quickly (e.g., Cs-137)
   * - No gamma-emitting progeny exist
   */
  bool is_age_fit_allowed( const SandiaDecay::Nuclide *nuc )
  {
    if( !nuc )
      return false;

    // Use the PeakDef function if available
    return !PeakDef::ageFitNotAllowed( nuc );
  }//is_age_fit_allowed(...)


  /** Count unique progeny nuclides with peaks assigned.
   * Only counts peaks with useForShieldingSourceFit() == true.
   */
  size_t count_progeny_peaks(
    const SandiaDecay::Nuclide *parent_nuc,
    const std::deque<std::shared_ptr<const PeakDef>> &peaks
  )
  {
    if( !parent_nuc )
      return 0;

    std::set<const SandiaDecay::Nuclide*> progeny_nuclides;

    for( const auto &peak : peaks )
    {
      if( !peak || !peak->useForShieldingSourceFit() )
        continue;

      if( peak->parentNuclide() != parent_nuc )
        continue;

      const SandiaDecay::Transition *trans = peak->nuclearTransition();
      if( trans && trans->parent )
      {
        // This is a progeny peak
        if( trans->parent != parent_nuc )
          progeny_nuclides.insert( trans->parent );
      }
    }

    return progeny_nuclides.size();
  }//count_progeny_peaks(...)

}//namespace


nlohmann::json ToolRegistry::executeCloseActivityShieldingDisplay(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for close_activity_shielding_display" );

  json result;
  result["success"] = false;

  // Get the ShieldingSourceDisplay from InterSpec
  ShieldingSourceDisplay *display = interspec->shieldingSourceFit();

  if( !display )
  {
    result["message"] = "Activity/Shielding fit GUI is not currently open";
    result["success"] = true; // Not an error - just wasn't open
    return result;
  }

  // Close the display (this will delete it)
  interspec->closeShieldingSourceFit();

  result["success"] = true;
  result["message"] = "Activity/Shielding fit GUI closed successfully";

  return result;
}//executeCloseActivityShieldingDisplay(...)


nlohmann::json ToolRegistry::executeAskUserQuestion(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for ask_user_question" );

  // Get parameters
  const string question = params.value( "question", string() );
  if( question.empty() )
    throw runtime_error( "question parameter is required" );

  const string default_response = params.value( "default_response", string() );

  json result;

  // Create a modal dialog
  Wt::WDialog *dialog = new Wt::WDialog( "Question from Analysis Agent" );
  dialog->setModal( true );
  dialog->rejectWhenEscapePressed( false );

  // Add question text
  new Wt::WText( question, dialog->contents() );
  dialog->contents()->addWidget( new Wt::WBreak() );
  dialog->contents()->addWidget( new Wt::WBreak() );

  // Add text area for response
  Wt::WTextArea *response_area = new Wt::WTextArea( dialog->contents() );
  response_area->setColumns( 60 );
  response_area->setRows( 5 );
  if( !default_response.empty() )
    response_area->setText( default_response );

  dialog->contents()->addWidget( new Wt::WBreak() );

  // Add submit button
  Wt::WPushButton *submit_btn = new Wt::WPushButton( "Submit", dialog->contents() );
  submit_btn->setDefault( true );

  // Result will be stored here
  string user_response;
  bool dialog_finished = false;

  // Connect submit button (note: clicked() provides a WMouseEvent parameter)
  submit_btn->clicked().connect( [dialog, response_area, &user_response, &dialog_finished](const Wt::WMouseEvent &) {
    user_response = response_area->text().toUTF8();
    dialog_finished = true;
    dialog->accept();
  });

  // Execute the dialog (blocking)
  dialog->exec();

  // Return the user's response
  result["response"] = user_response;
  result["success"] = true;

  return result;
}//executeAskUserQuestion(...)


nlohmann::json ToolRegistry::executeMarkPeaksForActivityFit(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for mark_peaks_for_activity_fit" );

  // Get parameters
  if( !params.contains("peak_energies") || !params["peak_energies"].is_array() )
    throw runtime_error( "peak_energies array parameter is required" );

  const bool use_for_fit = params.value( "use_for_fit", true );
  const string spec_type_str = params.value( "specType", string("Foreground") );

  SpecUtils::SpectrumType spec_type = SpecUtils::SpectrumType::Foreground;
  if( spec_type_str == "Background" )
    spec_type = SpecUtils::SpectrumType::Background;
  else if( spec_type_str == "Secondary" )
    spec_type = SpecUtils::SpectrumType::SecondForeground;

  // Get the PeakModel - this is the correct way to modify peaks
  PeakModel *peakModel = interspec->peakModel();
  if( !peakModel )
    throw runtime_error( "PeakModel not available" );

  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
  if( !peaks || peaks->empty() )
    throw runtime_error( "No peaks available for " + spec_type_str );

  json result;
  result["marked_peaks"] = json::array();
  result["errors"] = json::array();

  int num_marked = 0;

  // Process each energy
  for( const auto &energy_json : params["peak_energies"] )
  {
    const double energy = energy_json.get<double>();

    try
    {
      // Find the peak by energy
      bool found = false;
      for( size_t i = 0; i < peaks->size(); ++i )
      {
        const std::shared_ptr<const PeakDef> &peak = (*peaks)[i];
        if( !peak )
          continue;

        const double peak_energy = peak->mean();
        const double tolerance = std::max( peak->fwhm(), 1.0 );

        if( std::fabs( peak_energy - energy ) <= tolerance )
        {
          // Check if we need to actually change the flag
          if( peak->useForShieldingSourceFit() != use_for_fit )
          {
            // Use PeakModel::setData to properly modify the peak flag
            Wt::WModelIndex model_index = peakModel->index( i, PeakModel::kUseForShieldingSourceFit );
            peakModel->setData( model_index, boost::any(use_for_fit) );

            json marked_peak;
            marked_peak["energy"] = peak_energy;
            marked_peak["use_for_fit"] = use_for_fit;
            if( peak->parentNuclide() )
              marked_peak["nuclide"] = peak->parentNuclide()->symbol;
            result["marked_peaks"].push_back( marked_peak );

            ++num_marked;
          }

          found = true;
          break;
        }
      }//for( loop over peaks to find match )

      if( !found )
      {
        json error;
        error["energy"] = energy;
        error["error"] = "No peak found near " + std::to_string(energy) + " keV";
        result["errors"].push_back( error );
      }

    }catch( std::exception &e )
    {
      json error;
      error["energy"] = energy;
      error["error"] = e.what();
      result["errors"].push_back( error );
    }
  }//for( loop over peak energies )

  result["num_marked"] = num_marked;
  result["success"] = (num_marked > 0);

  return result;
}//executeMarkPeaksForActivityFit(...)


/** Helper function to convert ShieldingSourceDisplayState to JSON.
 * This function converts the configuration and peaks from a ShieldingSourceDisplayState
 * into a JSON format consistent with the modify_shielding_source_config tool.
 *
 * @param state The ShieldingSourceDisplayState to convert
 * @param peaks The peaks to use for counting progeny peaks; if nullptr, progeny peak counts will not be included
 * @return JSON object with configuration and peaks
 */
nlohmann::json format_shielding_source_state_to_json(
  const ShieldingSourceDisplay::ShieldingSourceDisplayState &state,
  const std::deque<std::shared_ptr<const PeakDef>> *peaks
)
{
  json result;

  if( !state.config )
    return result;

  const GammaInteractionCalc::ShieldSourceConfig &config = *state.config;

  // Format distance
  const double distance_cm = config.distance / PhysicalUnits::cm;
  result["distance_cm"] = distance_cm;

  // Format geometry
  switch( config.geometry )
  {
    case GammaInteractionCalc::GeometryType::Spherical:
      result["geometry"] = "Spherical";
      break;
    case GammaInteractionCalc::GeometryType::CylinderSideOn:
      result["geometry"] = "CylinderSideOn";
      break;
    case GammaInteractionCalc::GeometryType::CylinderEndOn:
      result["geometry"] = "CylinderEndOn";
      break;
    case GammaInteractionCalc::GeometryType::Rectangular:
      result["geometry"] = "Rectangular";
      break;
    case GammaInteractionCalc::GeometryType::NumGeometryType:
      result["geometry"] = "Unknown";
      break;
  }

  // Format shielding layers
  result["shielding"] = json::array();
  for( size_t i = 0; i < config.shieldings.size(); ++i )
  {
    const ShieldingSourceFitCalc::ShieldingInfo &shield = config.shieldings[i];
    json shield_json;

    shield_json["index"] = i;

    if( shield.m_material )
    {
      shield_json["material"] = shield.m_material->name;
      shield_json["is_generic"] = false;
    }
    else if( shield.m_isGenericMaterial )
    {
      shield_json["material"] = "Generic (AN=" + std::to_string(shield.m_dimensions[0]) + ")";
      shield_json["is_generic"] = true;
    }

    // Format dimensions based on geometry
    switch( shield.m_geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        shield_json["radial_thickness_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;
        shield_json["fit_thickness"] = shield.m_fitDimensions[0];
        break;

      case GammaInteractionCalc::GeometryType::CylinderSideOn:
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
        shield_json["radius_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;;
        shield_json["length_cm"] = shield.m_dimensions[1] / PhysicalUnits::cm;
        shield_json["fit_radius"] = shield.m_fitDimensions[0];
        shield_json["fit_length"] = shield.m_fitDimensions[1];
        break;

      case GammaInteractionCalc::GeometryType::Rectangular:
        shield_json["width_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;;
        shield_json["height_cm"] = shield.m_dimensions[1] / PhysicalUnits::cm;;
        shield_json["depth_cm"] = shield.m_dimensions[2] / PhysicalUnits::cm;;
        shield_json["fit_width"] = shield.m_fitDimensions[0];
        shield_json["fit_height"] = shield.m_fitDimensions[1];
        shield_json["fit_depth"] = shield.m_fitDimensions[2];
        break;

      case GammaInteractionCalc::GeometryType::NumGeometryType:
        break;
    }

    // Add trace sources if any
    if( !shield.m_traceSources.empty() )
    {
      shield_json["trace_sources"] = json::array();
      for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : shield.m_traceSources )
      {
        json trace_json;

        if( trace.m_nuclide )
          trace_json["nuclide"] = trace.m_nuclide->symbol;

        trace_json["activity"] = trace.m_activity;
        trace_json["fit_activity"] = trace.m_fitActivity;

        // Format trace activity type
        switch( trace.m_type )
        {
          case GammaInteractionCalc::TraceActivityType::TotalActivity:
            trace_json["activity_type"] = "TotalActivity";
            trace_json["units"] = "Bq";
            break;
          case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
            trace_json["activity_type"] = "ActivityPerCm3";
            trace_json["units"] = "Bq/cm3";
            break;
          case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
            trace_json["activity_type"] = "ExponentialDistribution";
            trace_json["units"] = "Bq/m2";
            trace_json["relaxation_distance_cm"] = trace.m_relaxationDistance;
            break;
          case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
            trace_json["activity_type"] = "ActivityPerGram";
            trace_json["units"] = "Bq/g";
            break;
          case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
            break;
        }

        shield_json["trace_sources"].push_back( trace_json );
      }
    }

    // Add self-attenuating (intrinsic) source nuclide fractions if any
    if( !shield.m_nuclideFractions_.empty() )
    {
      shield_json["self_atten_sources"] = json::array();
      for( const auto &elem_fractions : shield.m_nuclideFractions_ )
      {
        const SandiaDecay::Element *element = elem_fractions.first;
        const auto &nuclide_vec = elem_fractions.second;

        for( const auto &nuc_tuple : nuclide_vec )
        {
          const SandiaDecay::Nuclide *nuclide = std::get<0>( nuc_tuple );
          const double mass_fraction = std::get<1>( nuc_tuple );
          const bool fit_fraction = std::get<2>( nuc_tuple );

          json self_atten_json;
          self_atten_json["element"] = element->symbol;

          if( nuclide )
            self_atten_json["nuclide"] = nuclide->symbol;
          else
            self_atten_json["nuclide"] = "Other"; // Represents stable/other isotopes

          self_atten_json["mass_fraction"] = mass_fraction;
          self_atten_json["fit_mass_fraction"] = fit_fraction;

          shield_json["self_atten_sources"].push_back( self_atten_json );
        }
      }
    }

    result["shielding"].push_back( shield_json );
  }

  // Format sources
  result["sources"] = json::array();
  for( const ShieldingSourceFitCalc::SourceFitDef &src : config.sources )
  {
    json src_json;

    if( src.nuclide )
      src_json["nuclide"] = src.nuclide->symbol;

    src_json["activity_bq"] = src.activity;
    src_json["fit_activity"] = src.fitActivity;

    if( src.age > 0.0 )
    {
      src_json["age_years"] = src.age / PhysicalUnits::year;
      src_json["fit_age"] = src.fitAge;
    }

    // Format source type
    switch( src.sourceType )
    {
      case ShieldingSourceFitCalc::ModelSourceType::Point:
        src_json["source_type"] = "Point";
        break;
      case ShieldingSourceFitCalc::ModelSourceType::Trace:
        src_json["source_type"] = "Trace";
        break;
      case ShieldingSourceFitCalc::ModelSourceType::Intrinsic:
        src_json["source_type"] = "Intrinsic";
        break;
    }

    // Add age fit allowability if we have a nuclide
    if( src.nuclide )
    {
      src_json["age_potentially_fittable"] = is_age_fit_allowed( src.nuclide );

      // Add progeny peak count if peaks are provided
      if( peaks )
      {
        const size_t num_progeny = count_progeny_peaks( src.nuclide, *peaks );
        src_json["num_progeny_peaks"] = num_progeny;
      }
    }

    result["sources"].push_back( src_json );
  }

  // Format peaks
  result["peaks"] = json::array();
  for( const ShieldingSourceDisplay::ShieldingSourceDisplayState::Peak &peak : state.peaks )
  {
    json peak_json;
    peak_json["energy_kev"] = peak.energy;
    peak_json["use_for_fit"] = peak.use;
    if( !peak.nuclideSymbol.empty() )
      peak_json["nuclide"] = peak.nuclideSymbol;
    result["peaks"].push_back( peak_json );
  }

  // Format fit options
  result["options"] = json::object();
  result["options"]["attenuate_for_air"] = config.options.attenuate_for_air;
  result["options"]["multiple_nucs_contribute_to_peaks"] = config.options.multiple_nucs_contribute_to_peaks;
  result["options"]["background_peak_subtract"] = config.options.background_peak_subtract;
  result["options"]["same_age_isotopes"] = config.options.same_age_isotopes;
  result["options"]["account_for_decay_during_meas"] = config.options.account_for_decay_during_meas;

  return result;
}//format_shielding_source_state_to_json(...)


nlohmann::json ToolRegistry::executeGetShieldingSourceConfig(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for get_shielding_source_config" );

  json result;
  result["has_config"] = false;
  
  // We will always grab the configuration from the GUI - if we grab it from the foreground SpecMeas, it
  //  wont have picked up any source peaks we may have added since the last time the GUI was open - we could
  //  probably work our way around this and manually update things, but just to keep it all consistent...

  // Check if GUI is open (without creating it)
  ShieldingSourceDisplay *display = interspec->shieldingSourceFit( false );

  const bool was_open = !!display;
  if( !display )
    display = interspec->shieldingSourceFit( true );
  
  if( !display )
    throw runtime_error( "Unable to crate Shield/Source fit GUI to alter its configuration" );
  
  // Extract configuration from GUI by serializing it
  ShieldingSourceDisplay::ShieldingSourceDisplayState state = display->serialize();

  if( !state.config )
  {
    result["message"] = "GUI is open but has no configuration";
    result["source"] = "GUI (ShieldingSourceDisplay)";
    return result;
  }

  // Get peaks from InterSpec for progeny counting
  PeakModel *peakModel = interspec->peakModel();
  const std::deque<std::shared_ptr<const PeakDef>> *peaks = nullptr;
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks_ptr;
  if( peakModel )
  {
    peaks_ptr = peakModel->peaks();
    if( peaks_ptr )
      peaks = peaks_ptr.get();
  }

  // Convert state to JSON
  const json config_json = format_shielding_source_state_to_json( state, peaks );
  
  result["has_config"] = true;
  //result["source"] = "GUI (ShieldingSourceDisplay)";
  result["config"] = config_json;
  
  if( !was_open )
    interspec->closeShieldingSourceFit();
  
  return result;
  
  /*
  // Check saved model in SpecMeas
  std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !meas )
  {
    result["message"] = "No foreground spectrum loaded";
    return result;
  }

  const rapidxml::xml_document<char> *model_xml = meas->shieldingSourceModel();
  if( !model_xml )
  {
    result["message"] = "No saved shielding/source model found";
    return result;
  }

  // Parse the XML to get ShieldingSourceDisplayState
  try
  {
    const rapidxml::xml_node<char> *base_node = model_xml->first_node( "ShieldingSourceFit" );
    if( !base_node )
      throw runtime_error( "XML does not contain ShieldingSourceFit node" );

    MaterialDB *materialDb = interspec->materialDataBase();
    if( !materialDb )
      throw runtime_error( "Material database not available" );

    // Deserialize the XML into a state object
    ShieldingSourceDisplay::ShieldingSourceDisplayState state;
    state.deSerialize( base_node, materialDb );

    if( !state.config )
      throw runtime_error( "Failed to parse configuration from XML" );

    // Get peaks from InterSpec for progeny counting
    PeakModel *peakModel = interspec->peakModel();
    const std::deque<std::shared_ptr<const PeakDef>> *peaks = nullptr;
    std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks_ptr;
    if( peakModel )
    {
      peaks_ptr = peakModel->peaks();
      if( peaks_ptr )
        peaks = peaks_ptr.get();
    }

    // Convert state to JSON
    const json config_json = format_shielding_source_state_to_json( state, peaks );

    result["has_config"] = true;
    result["source"] = "Saved model (SpecMeas::m_shieldingSourceModel)";
    result["config"] = config_json;
  }
  catch( const std::exception &e )
  {
    result["has_config"] = false;
    result["source"] = "Saved model (SpecMeas::m_shieldingSourceModel)";
    result["message"] = string("Failed to parse configuration: ") + e.what();
  }
   */
  
  return result;
}//executeGetShieldingSourceConfig(...)


nlohmann::json ToolRegistry::executeModifyShieldingSourceConfig(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for modify_shielding_source_config" );

  // Get operation
  const string operation = params.value( "operation", string() );
  if( operation.empty() )
    throw runtime_error( "operation parameter is required" );

  json result;
  result["operation"] = operation;

  // Check if GUI is already open (without creating it)
  ShieldingSourceDisplay *display = interspec->shieldingSourceFit( false );
  const bool gui_was_open = (display != nullptr);

  // If GUI wasn't open, open it temporarily to get/modify state
  if( !display )
  {
    display = interspec->shieldingSourceFit( true );
    if( !display )
      throw runtime_error( "Failed to open Activity/Shielding Fit tool" );
  }

  // Get current state
  ShieldingSourceDisplay::ShieldingSourceDisplayState state = display->serialize();

  // Make a copy of the config to modify
  if( !state.config )
    throw runtime_error( "No configuration available to modify" );

  GammaInteractionCalc::ShieldSourceConfig modified_config = *state.config;

  // Apply the requested operation
  if( operation == "set_distance" )
  {
    const string distance_str = params.value( "distance", string() );
    if( distance_str.empty() )
      throw runtime_error( "distance parameter is required for set_distance operation" );

    modified_config.distance = parse_distance_string( distance_str );
    const double distance_cm = modified_config.distance / PhysicalUnits::cm;
    result["new_distance_cm"] = distance_cm;
  }
  else if( operation == "set_geometry" )
  {
    const string geometry_str = params.value( "geometry", string() );
    if( geometry_str.empty() )
      throw runtime_error( "geometry parameter is required for set_geometry operation" );

    if( geometry_str == "Spherical" )
      modified_config.geometry = GammaInteractionCalc::GeometryType::Spherical;
    else if( geometry_str == "CylinderSideOn" )
      modified_config.geometry = GammaInteractionCalc::GeometryType::CylinderSideOn;
    else if( geometry_str == "CylinderEndOn" )
      modified_config.geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
    else if( geometry_str == "Rectangular" )
      modified_config.geometry = GammaInteractionCalc::GeometryType::Rectangular;
    else
      throw runtime_error( "Invalid geometry: " + geometry_str );

    result["new_geometry"] = geometry_str;
  }
  else if( operation == "add_shielding" )
  {
    string material = params.value( "material", string() );
    if( material.empty() )
      throw runtime_error( "material parameter is required for add_shielding operation" );

    MaterialDB *materialDb = interspec->materialDataBase();
    if( !materialDb )
      throw runtime_error( "Material database not available" );

    ShieldingSourceFitCalc::ShieldingInfo shielding;
    shielding.m_geometry = modified_config.geometry;

    // Try to find the material and normalize the name to match what's in the database
    const Material *mat = materialDb->material( material );
    if( !mat )
      throw runtime_error( "Material '" + material + "' not found in database" );

    // Use the canonical material name from the database (e.g., "Fe (iron)" instead of "Fe")
    material = mat->name;

    shielding.m_material = std::make_shared<const Material>( *mat );
    shielding.m_isGenericMaterial = false;
    shielding.m_forFitting = true;

    // Handle dimensions based on geometry
    const string radial_thickness_str = params.value( "radial_thickness", string() );
    const string radius_str = params.value( "radius", string() );
    const string length_str = params.value( "length", string() );

    // Get fit flags
    const bool fit_thickness = params.value( "fit_thickness", radial_thickness_str.empty() );
    const bool fit_radius = params.value( "fit_radius", radius_str.empty() );
    const bool fit_length = params.value( "fit_length", length_str.empty() );

    switch( modified_config.geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
      {
        if( !radial_thickness_str.empty() )
          shielding.m_dimensions[0] = parse_distance_string( radial_thickness_str );
        else
          shielding.m_dimensions[0] = 0.5 * PhysicalUnits::cm; // Default to 0.5 cm

        shielding.m_fitDimensions[0] = fit_thickness;
        const double thickness_cm = shielding.m_dimensions[0] / PhysicalUnits::cm;
        result["thickness_cm"] = thickness_cm;
        result["fit_thickness"] = fit_thickness;
        break;
      }

      case GammaInteractionCalc::GeometryType::CylinderSideOn:
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      {
        if( !radius_str.empty() )
          shielding.m_dimensions[0] = parse_distance_string( radius_str );
        else
          shielding.m_dimensions[0] = 5.0 * PhysicalUnits::cm; // Default radius 5 cm

        if( !length_str.empty() )
          shielding.m_dimensions[1] = parse_distance_string( length_str );
        else
          shielding.m_dimensions[1] = 10.0 * PhysicalUnits::cm; // Default length 10 cm

        shielding.m_fitDimensions[0] = fit_radius;
        shielding.m_fitDimensions[1] = fit_length;
        const double radius_cm = shielding.m_dimensions[0] / PhysicalUnits::cm;
        const double length_cm = shielding.m_dimensions[1] / PhysicalUnits::cm;
        result["radius_cm"] = radius_cm;
        result["length_cm"] = length_cm;
        result["fit_radius"] = fit_radius;
        result["fit_length"] = fit_length;
        break;
      }

      case GammaInteractionCalc::GeometryType::Rectangular:
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        throw runtime_error( "Rectangular geometry not yet supported for add_shielding" );
    }

    modified_config.shieldings.push_back( shielding );
    result["material_added"] = material;
    result["num_shieldings"] = modified_config.shieldings.size();
  }
  else if( operation == "remove_shielding" )
  {
    string material = params.value( "material", string() );
    const int index = params.value( "index", -1 );

    if( material.empty() && index < 0 )
      throw runtime_error( "Either material or index parameter is required for remove_shielding operation" );

    if( !material.empty() && index >= 0 )
      throw runtime_error( "Cannot specify both material and index for remove_shielding - use one or the other" );

    if( index >= 0 )
    {
      // Remove by index
      if( index >= static_cast<int>(modified_config.shieldings.size()) )
        throw runtime_error( "Shielding index " + std::to_string(index) + " out of range (have " +
                           std::to_string(modified_config.shieldings.size()) + " shielding layers)" );

      const string removed_material = modified_config.shieldings[index].m_material
                                       ? modified_config.shieldings[index].m_material->name
                                       : "Generic";
      modified_config.shieldings.erase( modified_config.shieldings.begin() + index );
      result["index_removed"] = index;
      result["material_removed"] = removed_material;
    }
    else
    {
      // Try to normalize material name using database
      MaterialDB *materialDb = interspec->materialDataBase();
      if( materialDb )
      {
        const Material *mat = materialDb->material( material );
        if( mat )
          material = mat->name;  // Use canonical name from database
      }

      // Remove by material name - check name, description, or substring match
      auto it = std::find_if( modified_config.shieldings.begin(), modified_config.shieldings.end(),
        [&material]( const ShieldingSourceFitCalc::ShieldingInfo &s ) {
          if( !s.m_material )
            return false;
          // Check exact match of name or description
          if( s.m_material->name == material || s.m_material->description == material )
            return true;
          // Check if material is a substring of the name (e.g., "Fe" matches "Fe (iron)")
          if( s.m_material->name.find(material) != string::npos )
            return true;
          return false;
        });

      if( it == modified_config.shieldings.end() )
        throw runtime_error( "Shielding material '" + material + "' not found in configuration" );

      const string removed_material = it->m_material ? it->m_material->name : "Generic";
      const int removed_index = std::distance( modified_config.shieldings.begin(), it );
      modified_config.shieldings.erase( it );
      result["material_removed"] = removed_material;
      result["index_removed"] = removed_index;
    }

    result["num_shieldings"] = modified_config.shieldings.size();
  }
  else if( operation == "add_source" )
  {
    const string nuclide_str = params.value( "nuclide", string() );
    if( nuclide_str.empty() )
      throw runtime_error( "nuclide parameter is required for add_source operation" );

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuclide = db->nuclide( nuclide_str );
    if( !nuclide )
      throw runtime_error( "Nuclide '" + nuclide_str + "' not found in database" );

    // Check if source already exists
    auto existing = std::find_if( modified_config.sources.begin(), modified_config.sources.end(),
      [nuclide]( const ShieldingSourceFitCalc::SourceFitDef &s ) {
        return s.nuclide == nuclide;
      });

    if( existing != modified_config.sources.end() )
      throw runtime_error( "Source '" + nuclide_str + "' already exists in configuration" );

    ShieldingSourceFitCalc::SourceFitDef source;
    source.nuclide = nuclide;

    // Parse activity
    const string activity_str = params.value( "activity", string() );
    if( !activity_str.empty() )
      source.activity = parse_activity_string( activity_str );
    else
      source.activity = 1.0 * PhysicalUnits::microCi; // Default

    source.fitActivity = true; // Always fit activity by default

    // Parse age if provided
    const string age_str = params.value( "age", string() );
    if( !age_str.empty() )
    {
      source.age = parse_age_string( age_str );
      source.fitAge = params.value( "fit_age", false );
    }
    else
    {
      source.age = 0.0;
      source.fitAge = false;
    }

    // Parse source type
    const string source_type_str = params.value( "source_type", "Point" );
    if( source_type_str == "Point" )
      source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
    else if( source_type_str == "Trace" )
      source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Trace;
    else if( source_type_str == "Intrinsic" )
      source.sourceType = ShieldingSourceFitCalc::ModelSourceType::Intrinsic;
    else
      throw runtime_error( "Invalid source_type: " + source_type_str );

    modified_config.sources.push_back( source );

    // Mark all peaks with this nuclide to be used for fitting
    PeakModel *peakModel = interspec->peakModel();
    if( peakModel )
    {
      std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
      if( peaks )
      {
        int num_peaks_marked = 0;
        int total_peaks_for_nuclide = 0;
        for( size_t i = 0; i < peaks->size(); ++i )
        {
          const std::shared_ptr<const PeakDef> &peak = (*peaks)[i];
          if( peak && peak->parentNuclide() == nuclide )
          {
            ++total_peaks_for_nuclide;
            if( !peak->useForShieldingSourceFit() )
            {
              // Create model index for this peak's kUseForShieldingSourceFit column
              Wt::WModelIndex model_index = peakModel->index( i, PeakModel::kUseForShieldingSourceFit );
              // Set the value to true (checked)
              peakModel->setData( model_index, boost::any(true) );
              ++num_peaks_marked;
            }
          }
        }
        result["num_peaks_marked"] = num_peaks_marked;
      }
    }

    result["nuclide_added"] = nuclide_str;
    const double activity_bq = source.activity / PhysicalUnits::bq;
    result["activity_bq"] = activity_bq;
    result["num_sources"] = modified_config.sources.size();
  }
  else if( operation == "remove_source" )
  {
    const string nuclide_str = params.value( "nuclide", string() );
    if( nuclide_str.empty() )
      throw runtime_error( "nuclide parameter is required for remove_source operation" );

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuclide = db->nuclide( nuclide_str );
    if( !nuclide )
      throw runtime_error( "Nuclide '" + nuclide_str + "' not found in database" );

    // Find and remove the source
    auto it = std::find_if( modified_config.sources.begin(), modified_config.sources.end(),
      [nuclide]( const ShieldingSourceFitCalc::SourceFitDef &s ) {
        return s.nuclide == nuclide;
      });

    if( it == modified_config.sources.end() )
      throw runtime_error( "Source '" + nuclide_str + "' not found in configuration" );

    modified_config.sources.erase( it );

    // Unmark all peaks with this nuclide from being used for fitting
    PeakModel *peakModel = interspec->peakModel();
    if( peakModel )
    {
      std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
      if( peaks )
      {
        int num_peaks_unmarked = 0;
        for( size_t i = 0; i < peaks->size(); ++i )
        {
          const std::shared_ptr<const PeakDef> &peak = (*peaks)[i];
          if( peak && peak->parentNuclide() == nuclide && peak->useForShieldingSourceFit() )
          {
            // Create model index for this peak's kUseForShieldingSourceFit column
            Wt::WModelIndex model_index = peakModel->index( i, PeakModel::kUseForShieldingSourceFit );
            // Set the value to false (unchecked)
            peakModel->setData( model_index, boost::any(false) );
            ++num_peaks_unmarked;
          }
        }
        result["num_peaks_unmarked"] = num_peaks_unmarked;
      }
    }

    result["nuclide_removed"] = nuclide_str;
    result["num_sources"] = modified_config.sources.size();
  }
  else if( operation == "set_source_age" )
  {
    const string nuclide_str = params.value( "nuclide", string() );
    if( nuclide_str.empty() )
      throw runtime_error( "nuclide parameter is required for set_source_age operation" );

    const string age_str = params.value( "age", string() );
    if( age_str.empty() )
      throw runtime_error( "age parameter is required for set_source_age operation" );

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const SandiaDecay::Nuclide *nuclide = db->nuclide( nuclide_str );
    if( !nuclide )
      throw runtime_error( "Nuclide '" + nuclide_str + "' not found in database" );

    // Find the source
    auto it = std::find_if( modified_config.sources.begin(), modified_config.sources.end(),
      [nuclide]( const ShieldingSourceFitCalc::SourceFitDef &s ) {
        return s.nuclide == nuclide;
      });

    if( it == modified_config.sources.end() )
      throw runtime_error( "Source '" + nuclide_str + "' not found in configuration" );

    it->age = parse_age_string( age_str );
    it->fitAge = params.value( "fit_age", false );
    result["nuclide"] = nuclide_str;
    result["age_years"] = it->age / PhysicalUnits::year;
    result["fit_age"] = it->fitAge;
  }
  else if( operation == "set_fit_options" )
  {
    // Set fitting options
    bool options_changed = false;

    if( params.contains("attenuate_for_air") )
    {
      const bool new_value = get_boolean( params, "attenuate_for_air" );
      modified_config.options.attenuate_for_air = new_value;
      result["attenuate_for_air"] = new_value;
      options_changed = true;
    }

    if( params.contains("multiple_nucs_contribute_to_peaks") )
    {
      const bool new_value = get_boolean( params, "multiple_nucs_contribute_to_peaks" );
      modified_config.options.multiple_nucs_contribute_to_peaks = new_value;
      result["multiple_nucs_contribute_to_peaks"] = new_value;
      options_changed = true;
    }

    if( params.contains("background_peak_subtract") )
    {
      const bool new_value = get_boolean( params, "background_peak_subtract" );
      modified_config.options.background_peak_subtract = new_value;
      result["background_peak_subtract"] = new_value;
      options_changed = true;
    }

    if( params.contains("same_age_isotopes") )
    {
      const bool new_value = get_boolean( params, "same_age_isotopes" );
      modified_config.options.same_age_isotopes = new_value;
      result["same_age_isotopes"] = new_value;
      options_changed = true;
    }

    if( !options_changed )
      throw runtime_error( "set_fit_options requires at least one option parameter (attenuate_for_air, multiple_nucs_contribute_to_peaks, background_peak_subtract, or same_age_isotopes)" );
  }
  else
  {
    throw runtime_error( "Unknown operation: " + operation );
  }

  // Update the state with the modified config
  state.config = std::make_shared<const GammaInteractionCalc::ShieldSourceConfig>( modified_config );

  // Deserialize the modified state back to the GUI
  display->deSerialize( state, 0 );

  // Save to SpecMeas
  interspec->saveShieldingSourceModelToForegroundSpecMeas();

  // If GUI wasn't originally open, close it
  if( !gui_was_open )
    interspec->closeShieldingSourceFit();

  result["success"] = true;
  result["gui_was_open"] = gui_was_open;
  result["message"] = "Configuration modified successfully";

  return result;
}//executeModifyShieldingSourceConfig(...)


/** Helper function to format ModelFitResults into a JSON object.
 * This function converts fit results (activities, ages, shielding, fit quality)
 * into a structured JSON format for returning to the LLM.
 *
 * @param fit_results Pointer to ModelFitResults (can be null)
 * @return JSON object with formatted results, or empty object if fit_results is null
 */
nlohmann::json format_fit_results_to_json(
  const std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> &fit_results
)
{
  json result;

  if( !fit_results )
    return result;

  // Format fit quality
  result["chi2"] = fit_results->chi2;
  result["dof"] = fit_results->numDOF;
  result["chi2_per_dof"] = (fit_results->numDOF > 0) ? (fit_results->chi2 / fit_results->numDOF) : fit_results->chi2;

  // Format source results
  result["sources"] = json::array();
  for( const auto &src : fit_results->fit_src_info )
  {
    json src_json;
    src_json["nuclide"] = src.nuclide->symbol;
    src_json["activity_bq"] = src.activity;

    if( src.activityUncertainty.has_value() )
      src_json["activity_bq_uncertainty"] = *src.activityUncertainty;

    
    // Format as human-readable string
    InterSpec * const interspec = InterSpec::instance();
    const bool use_curries = interspec ? true : !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );
    if( src.activityUncertainty.has_value() && (src.activityUncertainty.value() > 0.0) )
      src_json["activity_str"] = PhysicalUnits::printToBestActivityUnitsWithUncert( src.activity, src.activityUncertainty.value(), 5, use_curries );
    else
      src_json["activity_str"] = PhysicalUnits::printToBestActivityUnits(src.activity, 5, use_curries);

    // Add age if fitted
    if( src.ageUncertainty.has_value() && *src.ageUncertainty > 0.0 )
    {
      const double age_years = src.age / PhysicalUnits::year;
      const double age_unc_years = *src.ageUncertainty / PhysicalUnits::year;
      src_json["age_years"] = age_years;
      src_json["age_years_uncertainty"] = age_unc_years;
      src_json["age_str"] = std::to_string(age_years) + " ± " + std::to_string(age_unc_years) + " years";
    }

    result["sources"].push_back( src_json );
  }

  // Format shielding results if any
  if( !fit_results->final_shieldings.empty() )
  {
    result["shielding"] = json::array();
    for( const auto &shield : fit_results->final_shieldings )
    {
      json shield_json;

      if( shield.m_material )
        shield_json["material"] = shield.m_material->name;
      else if( shield.m_isGenericMaterial )
        shield_json["material"] = "Generic (AN=" + std::to_string(shield.m_dimensions[0]) + ")";

      // Format dimensions based on geometry
      switch( shield.m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
        {
          shield_json["thickness_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;
          shield_json["thickness_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[0], 4 );
          if( shield.m_dimensionUncerts[0] > 0.0 )
          {
            shield_json["thickness_cm_uncertainty"] = shield.m_dimensionUncerts[0] / PhysicalUnits::cm;
            shield_json["thickness_str_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[0], 4 );
          }
          break;
        }//case GammaInteractionCalc::GeometryType::Spherical:

        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
        {
          shield_json["radius_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;
          shield_json["radius_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[0], 4 );
          shield_json["length_cm"] = shield.m_dimensions[1] / PhysicalUnits::cm;;
          shield_json["length_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[1], 4 );
          
          if( shield.m_dimensionUncerts[0] > 0.0 )
          {
            shield_json["radius_cm_uncertainty"] = shield.m_dimensionUncerts[0] / PhysicalUnits::cm;
            shield_json["radius_str_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[0], 4 );
          }
          
          if( shield.m_dimensionUncerts[1] > 0.0 )
          {
            shield_json["length_cm_uncertainty"] = shield.m_dimensionUncerts[1] / PhysicalUnits::cm;
            shield_json["length_str_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[1], 4 );
          }
         
          break;
        }//case GammaInteractionCalc::GeometryType::CylinderEndOn: / CylinderSideOn:

        case GammaInteractionCalc::GeometryType::Rectangular:
        {
          shield_json["width_cm"] = shield.m_dimensions[0] / PhysicalUnits::cm;
          shield_json["height_cm"] = shield.m_dimensions[1] / PhysicalUnits::cm;
          shield_json["depth_cm"] = shield.m_dimensions[2] / PhysicalUnits::cm;
          
          shield_json["width_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[0], 4 );
          shield_json["height_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[1], 4 );
          shield_json["depth_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensions[2], 4 );
          
          
          if( shield.m_dimensionUncerts[0] > 0.0 )
          {
            shield_json["width_cm_uncertainty"] = shield.m_dimensionUncerts[0] / PhysicalUnits::cm;
            shield_json["width_cm_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[0], 4 );
          }
          
          if( shield.m_dimensionUncerts[1] > 0.0 )
          {
            shield_json["height_cm_uncertainty"] = shield.m_dimensionUncerts[1] / PhysicalUnits::cm;
            shield_json["height_cm_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[1], 4 );
          }
          
          if( shield.m_dimensionUncerts[2] > 0.0 )
          {
            shield_json["depth_cm_uncertainty"] = shield.m_dimensionUncerts[2] / PhysicalUnits::cm;
            shield_json["depth_cm_uncertainty"] = PhysicalUnits::printToBestLengthUnits( shield.m_dimensionUncerts[2], 4 );
          }
          
          break;
        }//case GammaInteractionCalc::GeometryType::Rectangular:

        case GammaInteractionCalc::GeometryType::NumGeometryType:
          break;
      }

      result["shielding"].push_back( shield_json );
    }
  }

  result["num_peaks_used"] = fit_results->foreground_peaks.size();

  return result;
}//format_fit_results_to_json(...)


nlohmann::json ToolRegistry::executeActivityFit(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for activity_fit" );

  // Get mode
  const string mode = params.value( "mode", string() );
  if( mode.empty() )
    throw runtime_error( "mode parameter is required" );

  if( mode != "from_app_state" && mode != "single_peak" && mode != "custom" )
    throw runtime_error( "mode must be 'from_app_state', 'single_peak', or 'custom'" );

  json result;
  result["status"] = "failed";
  result["mode"] = mode;

  // Handle from_app_state mode - use/create GUI
  if( mode == "from_app_state" )
  {
    // Check if GUI is already open
    ShieldingSourceDisplay *display = interspec->shieldingSourceFit();

    if( !display )
    {
      // For now, creating the GUI requires more complex setup with viewer, etc.
      // Return an error asking user to open the GUI manually first
      throw runtime_error( "Activity/Shielding fit GUI is not currently open. "
                         "Please open it manually from the Tools menu, or use 'custom' or 'single_peak' modes for direct fitting." );
    }

    // TODO: Apply any parameter overrides from params (distance, add shielding, etc.)
    // For now, just trigger the fit with existing configuration

    // Trigger the fit - doModelFit returns results when fit completes
    std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> fit_results = display->doModelFit( false );

    result["status"] = "success";
    result["gui_displayed"] = true;

    // If we got results back, format and return them using the helper function
    if( fit_results )
    {
      result["message"] = "Activity/Shielding fit completed successfully in GUI.";

      // Format the results using the helper function
      const json formatted_results = format_fit_results_to_json( fit_results );

      // Merge the formatted results into our result object
      result["chi2"] = formatted_results["chi2"];
      result["dof"] = formatted_results["dof"];
      result["chi2_per_dof"] = formatted_results["chi2_per_dof"];
      result["sources"] = formatted_results["sources"];
      result["num_peaks_used"] = formatted_results["num_peaks_used"];

      if( formatted_results.contains("shielding") )
        result["shielding"] = formatted_results["shielding"];
    }
    else
    {
      result["message"] = "Activity/Shielding fit started in GUI. Results will be displayed when fit completes.";
    }

    return result;
  }

  // For single_peak and custom modes: build ShieldSourceInput and fit directly
  // Get foreground measurement
  std::shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    throw runtime_error( "No foreground spectrum loaded" );

  std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !meas )
    throw runtime_error( "No foreground SpecMeas available" );

  const std::set<int> &sample_nums = interspec->displayedSamples( SpecUtils::SpectrumType::Foreground );

  // Get background if available
  std::shared_ptr<const SpecUtils::Measurement> background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> background_peaks;
  if( background )
  {
    std::shared_ptr<SpecMeas> back_meas = interspec->measurment( SpecUtils::SpectrumType::Background );
    if( back_meas )
    {
      const std::set<int> &back_samples = interspec->displayedSamples( SpecUtils::SpectrumType::Background );
      background_peaks = back_meas->peaks( back_samples );
    }
  }

  // Get detector response function
  std::shared_ptr<const DetectorPeakResponse> detector = meas->detector();
  if( !detector || !detector->isValid() )
    throw runtime_error( "No valid detector efficiency function loaded. Please load one before fitting." );

  // Get peaks
  std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> all_peaks = meas->peaks( sample_nums );
  if( !all_peaks || all_peaks->empty() )
    throw runtime_error( "No peaks available for fitting" );

  // Filter peaks to use for fitting
  std::deque<std::shared_ptr<const PeakDef>> fitting_peaks;

  if( params.contains("peak_energies") && params["peak_energies"].is_array() )
  {
    // Use specified peaks
    for( const auto &energy_json : params["peak_energies"] )
    {
      const double energy = energy_json.get<double>();
      auto peak = find_peak_by_energy( energy, all_peaks, "activity_fit" );
      if( !peak )
        throw runtime_error( "No peak found near " + std::to_string(energy) + " keV" );

      if( !peak->parentNuclide() )
        throw runtime_error( "Peak at " + std::to_string(peak->mean())
                           + " keV does not have a nuclide assigned" );

      fitting_peaks.push_back( peak );
    }
  }else
  {
    // Use peaks marked for fitting
    for( const auto &peak : *all_peaks )
    {
      if( peak && peak->useForShieldingSourceFit() && peak->parentNuclide() )
        fitting_peaks.push_back( peak );
    }
  }

  if( fitting_peaks.empty() )
    throw runtime_error( "No peaks available for fitting. Either specify peak_energies or mark peaks with useForShieldingSourceFit" );

  // Build ShieldSourceInput
  GammaInteractionCalc::ShieldingSourceChi2Fcn::ShieldSourceInput chi_input;

  // Set measurements
  chi_input.detector = detector;
  chi_input.foreground = foreground;
  chi_input.background = background;
  chi_input.foreground_peaks.assign( fitting_peaks.begin(), fitting_peaks.end() );
  chi_input.background_peaks = background_peaks;

  // Parse distance
  if( detector->isFixedGeometry() )
  {
    chi_input.config.distance = 0.0;
  }else if( !params.contains("distance") )
  {
    throw runtime_error( "distance parameter is required for non-fixed-geometry detectors" );
  }else
  {
    const string distance_str = params["distance"].is_string()
                                 ? params["distance"].get<string>()
                                 : std::to_string( params["distance"].get<double>() );
    chi_input.config.distance = parse_distance_string( distance_str );
  }

  // Parse geometry
  const string geometry_str = params.value( "geometry", string("Spherical") );
  if( geometry_str == "Spherical" )
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::Spherical;
  else if( geometry_str == "CylinderSideOn" )
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::CylinderSideOn;
  else if( geometry_str == "CylinderEndOn" )
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::CylinderEndOn;
  else if( geometry_str == "Rectangular" )
    chi_input.config.geometry = GammaInteractionCalc::GeometryType::Rectangular;
  else
    throw runtime_error( "Invalid geometry: " + geometry_str );

  // Parse options
  if( params.contains("options") && params["options"].is_object() )
  {
    const json &opts = params["options"];
    chi_input.config.options.multiple_nucs_contribute_to_peaks = opts.value( "multiple_nucs_contribute_to_peaks", true );
    chi_input.config.options.attenuate_for_air = opts.value( "attenuate_for_air", true );
    chi_input.config.options.account_for_decay_during_meas = opts.value( "account_for_decay_during_meas", false );
    chi_input.config.options.photopeak_cluster_sigma = opts.value( "photopeak_cluster_sigma", 1.25 );
    chi_input.config.options.background_peak_subtract = opts.value( "background_peak_subtract", true );
    chi_input.config.options.same_age_isotopes = opts.value( "same_age_isotopes", true );
  }

  // Build source definitions from peaks (simple case - one source per nuclide)
  std::map<const SandiaDecay::Nuclide*, ShieldingSourceFitCalc::SourceFitDef> sources_map;
  for( const auto &peak : fitting_peaks )
  {
    const SandiaDecay::Nuclide *nuc = peak->parentNuclide();
    if( !nuc )
      continue;

    if( sources_map.find(nuc) == sources_map.end() )
    {
      ShieldingSourceFitCalc::SourceFitDef srcdef;
      srcdef.nuclide = nuc;
      srcdef.activity = 1.0E-6 * PhysicalUnits::curie; // Initial guess: 1 uCi
      srcdef.fitActivity = true;
      srcdef.age = PeakDef::defaultDecayTime( nuc );
      srcdef.fitAge = false; // Default: don't fit age
      srcdef.ageDefiningNuc = nullptr;
      srcdef.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;

      sources_map[nuc] = srcdef;
    }
  }

  // Convert map to vector and set sources
  chi_input.config.sources.clear();
  for( auto &pair : sources_map )
    chi_input.config.sources.push_back( pair.second );

  // TODO: Parse shielding from parameters (for now, no shielding)
  chi_input.config.shieldings.clear();

  // Create the chi2 function and parameters
  auto [chi2Fcn, inputParams] = GammaInteractionCalc::ShieldingSourceChi2Fcn::create( chi_input );

  if( !chi2Fcn )
    throw runtime_error( "Failed to create chi2 function for fitting" );

  // Perform the fit (synchronous)
  auto progress = std::make_shared<ShieldingSourceFitCalc::ModelFitProgress>();
  auto fit_results = std::make_shared<ShieldingSourceFitCalc::ModelFitResults>();
  fit_results->initial_shieldings = chi_input.config.shieldings;

  ShieldingSourceFitCalc::fit_model(
    "",  // Empty wtsession = synchronous
    chi2Fcn,
    std::make_shared<ROOT::Minuit2::MnUserParameters>(inputParams),
    progress,
    [](){},  // No progress callback
    fit_results,
    [](){}   // No finished callback
  );

  // Format results
  result["status"] = "success";
  result["chi2"] = fit_results->chi2;
  result["dof"] = fit_results->numDOF;
  result["chi2_per_dof"] = (fit_results->numDOF > 0) ? (fit_results->chi2 / fit_results->numDOF) : 0.0;

  // Add source results
  result["sources"] = json::array();
  for( const auto &src : fit_results->fit_src_info )
  {
    json src_json;
    src_json["nuclide"] = src.nuclide->symbol;
    src_json["activity_bq"] = src.activity;
    if( src.activityUncertainty.has_value() )
      src_json["activity_bq_uncertainty"] = *src.activityUncertainty;

    // Format as human-readable string
    InterSpec * const interspec = InterSpec::instance();
    const bool use_curries = interspec ? true : !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );
    if( src.activityUncertainty.has_value() && (src.activityUncertainty.value() > 0.0) )
      src_json["activity_str"]= PhysicalUnits::printToBestActivityUnitsWithUncert( src.activity, src.activityUncertainty.value(), 5, use_curries );
    else
      src_json["activity_str"] = PhysicalUnits::printToBestActivityUnits( src.activity, 5, use_curries );

    result["sources"].push_back( src_json );
  }

  result["message"] = "Fit completed successfully";
  result["num_peaks_used"] = fitting_peaks.size();

  return result;
}//executeActivityFit(...)


} // namespace LlmTools

#endif // USE_LLM_INTERFACE
