#include "InterSpec_config.h"
#include "InterSpec/LlmActivityFitTool.h"

#if( USE_LLM_INTERFACE )

#include <algorithm>
#include <map>
#include <set>
#include <deque>
#include <memory>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/LlmInterface.h"
#include "InterSpec/LlmToolGui.h"
#include "InterSpec/LlmConversationHistory.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"

#include "Minuit2/MnUserParameters.h"

#include <Wt/WDialog>
#include <Wt/WTextArea>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

using namespace std;
using json = nlohmann::json;

namespace LlmTools {
namespace ActivityFitTool {

namespace {
  // Helper functions for activity/shielding fitting

  /** Parse a boolean value from a JSON object field.

   Handles boolean, string ("true"/"false"), and numeric (0/non-zero) values.

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


  /** Parse a trace source activity string that may include per-volume, per-mass, or per-area units.
   *
   * Examples:
   *   "10 uCi" -> TotalActivity, 10 uCi
   *   "1 Bq/cm3" -> ActivityPerCm3, 1 Bq
   *   "5 kBq per cm3" -> ActivityPerCm3, 5000 Bq
   *   "100 Bq/g" -> ActivityPerGram, 100 Bq
   *   "50 Bq per gram" -> ActivityPerGram, 50 Bq
   *   "20 Bq/m2" -> ActivityPerM2, 20 Bq
   *
   * Returns a pair of (activity_type, activity_value_in_bq)
   */
  std::pair<GammaInteractionCalc::TraceActivityType, double> parse_trace_activity_string( const string &activity_str )
  {
    string trimmed = activity_str;
    SpecUtils::trim( trimmed );

    // Convert to lowercase for case-insensitive matching
    string lower = trimmed;
    std::transform( lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c){ return std::tolower(c); } );

    // Check for per-volume units: /cm3, /cm^3, per cm3, per cubic centimeter, etc.
    const bool has_per_cm3 = (lower.find("/cm3") != string::npos) ||
                             (lower.find("/cm^3") != string::npos) ||
                             (lower.find("/cc") != string::npos) ||
                             (lower.find("per cm3") != string::npos) ||
                             (lower.find("per cm^3") != string::npos) ||
                             (lower.find("per cc") != string::npos) ||
                             (lower.find("per cubic centimeter") != string::npos);

    // Check for per-mass units: /g, /gram, per gram, per g, /kg, per kilogram, etc.
    const bool has_per_gram = (lower.find("/g") != string::npos) ||
                              (lower.find("/gram") != string::npos) ||
                              (lower.find("per g") != string::npos) ||
                              (lower.find("per gram") != string::npos) ||
                              (lower.find("/kg") != string::npos) ||
                              (lower.find("per kg") != string::npos) ||
                              (lower.find("per kilogram") != string::npos);

    // Check for per-area units: /m2, /m^2, per m2, per square meter, etc.
    const bool has_per_m2 = (lower.find("/m2") != string::npos) ||
                            (lower.find("/m^2") != string::npos) ||
                            (lower.find("per m2") != string::npos) ||
                            (lower.find("per m^2") != string::npos) ||
                            (lower.find("per square meter") != string::npos);

    // Determine activity type based on units found
    GammaInteractionCalc::TraceActivityType activity_type = GammaInteractionCalc::TraceActivityType::TotalActivity;

    if( has_per_cm3 )
      activity_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
    else if( has_per_gram )
      activity_type = GammaInteractionCalc::TraceActivityType::ActivityPerGram;
    else if( has_per_m2 )
    {
      throw runtime_error( "ActivityPerM2 is not currently supported for trace sources. "
                          "Use TotalActivity, ActivityPerCm3, ActivityPerGram, or ExponentialDistribution." );
    }

    // Extract just the activity portion (before the "/" or "per")
    string activity_only = trimmed;

    // Find the position where per-unit portion starts
    size_t per_pos = string::npos;

    // Try to find "/" first
    per_pos = activity_only.find('/');

    // If no "/" found, try "per" (case-insensitive)
    if( per_pos == string::npos )
    {
      size_t per_word_pos = lower.find(" per ");
      if( per_word_pos != string::npos )
        per_pos = per_word_pos;
    }

    // Extract activity portion
    if( per_pos != string::npos )
    {
      activity_only = activity_only.substr( 0, per_pos );
      SpecUtils::trim( activity_only );
    }

    // Parse the activity value
    double activity_value = 0.0;
    try
    {
      activity_value = PhysicalUnits::stringToActivity( activity_only );
    }
    catch( std::exception &e )
    {
      throw runtime_error( "Could not parse trace activity '" + activity_str + "': " + e.what() );
    }

    return std::make_pair( activity_type, activity_value );
  }//parse_trace_activity_string(...)



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


nlohmann::json executeCloseActivityShieldingDisplay(
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


nlohmann::json executeAskUserQuestion(
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


nlohmann::json executeMarkPeaksForActivityFit(
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
            Wt::WModelIndex peak_index = peakModel->indexOfPeak(peak);
            assert( peak_index.isValid() );
            Wt::WModelIndex model_index = peakModel->index( peak_index.row(), PeakModel::kUseForShieldingSourceFit );
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


nlohmann::json executeGetShieldingSourceConfig(
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


nlohmann::json executeModifyShieldingSourceConfig(
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
      const long int removed_index = std::distance( modified_config.shieldings.begin(), it );
      modified_config.shieldings.erase( it );
      result["material_removed"] = removed_material;
      result["index_removed"] = removed_index;
    }

    result["num_shieldings"] = modified_config.shieldings.size();
  }
  /*
   // 20251208: adding a source through this command maybe introduces more arbitraryness than is useful, so we'll skip htis for now
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
    if( !peakModel )
      throw runtime_error( "No currently active peak model - this is really not expected" );
    
    shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
    if( !peaks )
      throw std::runtime_error( "Could not add " + nuclide->symbol + " as a source - there are currently no analysis peaks" );
    
    int num_peaks_marked = 0, total_peaks_for_nuclide = 0;
    for( const shared_ptr<const PeakDef> &peak : *peaks )
    {
      assert( peak );
      total_peaks_for_nuclide += (peak && (peak->parentNuclide() == nuclide));
      assert( !peak || (peak->parentNuclide() != nuclide) || !peak->useForShieldingSourceFit() );
      
      if( !peak || (peak->parentNuclide() != nuclide) || peak->useForShieldingSourceFit() )
        continue;
      
      if( !peak->useForShieldingSourceFit() )
      {
        Wt::WModelIndex peak_index = peakModel->indexOfPeak(peak);
        Wt::WModelIndex model_index = peakModel->index( peak_index.row(), PeakModel::kUseForShieldingSourceFit );
        peakModel->setData( model_index, boost::any(true) );
        ++num_peaks_marked;
      }//if( peak && (peak->parentNuclide() == nuclide) )
    }//for( size_t i = 0; i < peaks->size(); ++i )
    
    assert( num_peaks_marked == total_peaks_for_nuclide );
    
    if( !total_peaks_for_nuclide )
      throw std::runtime_error( "Could not add " + nuclide->symbol + " as a source - no peaks with this nuclide associated with it are present." );

    result["num_peaks_marked"] = num_peaks_marked;
    result["nuclide_added"] = nuclide_str;
    const double activity_bq = source.activity / PhysicalUnits::bq;
    result["activity_bq"] = activity_bq;
    result["num_sources"] = modified_config.sources.size();
  }
   */
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

/** Add activity information to JSON for LLM consumption.
 * Handles Point, Trace, and Self-Attenuating sources differently.
 *
 * @param src Source details from fit results
 * @param use_bq If true, use Becquerel for display; if false, use Curies
 * @param shield_details Shielding details (needed for self-atten fit status logic), can be nullptr
 * @param activity_json JSON object to populate with activity fields
 */
void add_activity_to_json(
  const GammaInteractionCalc::SourceDetails &src,
  const bool use_bq,
  const std::vector<GammaInteractionCalc::ShieldingDetails> * const shield_details,
  nlohmann::json &activity_json
)
{
  // Handle Trace sources (contamination in shielding)
  if( src.isTraceSource )
  {
    // Get postfix based on trace activity type
    std::string trace_postfix;
    switch( src.traceActivityType )
    {
      case GammaInteractionCalc::TraceActivityType::TotalActivity:
        trace_postfix = "";
        break;
      case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
        trace_postfix = "/cm³";
        break;
      case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
        trace_postfix = "/m² exp";
        break;
      case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
        trace_postfix = "/g";
        break;
      case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
        break;
    }//switch( src.traceActivityType )

    activity_json["type"] = GammaInteractionCalc::to_str( src.traceActivityType );
    activity_json["postfix"] = trace_postfix;

    // Display activity (per-unit: Bq/cm³, Bq/g, etc.)
    activity_json["display_bq"] = src.traceSrcDisplayAct / PhysicalUnits::bq;
    activity_json["display_ci"] = src.traceSrcDisplayAct / PhysicalUnits::ci;
    activity_json["display_str"] = PhysicalUnits::printToBestActivityUnits( src.traceSrcDisplayAct, 4, !use_bq ) + trace_postfix;

    // Total activity in entire shielding volume
    activity_json["total_bq"] = src.activity / PhysicalUnits::bq;
    activity_json["total_ci"] = src.activity / PhysicalUnits::ci;
    activity_json["total_str"] = PhysicalUnits::printToBestActivityUnits( src.activity, 4, !use_bq );

    activity_json["was_fit"] = src.activityIsFit;

    if( src.activityIsFit )
    {
      activity_json["display_uncertainty_bq"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::bq;
      activity_json["display_uncertainty_ci"] = src.traceSrcDisplayActUncertainty / PhysicalUnits::ci;
      activity_json["display_uncertainty_str"] = PhysicalUnits::printToBestActivityUnits( src.traceSrcDisplayActUncertainty, 4, !use_bq ) + trace_postfix;

      activity_json["total_uncertainty_bq"] = src.activityUncertainty / PhysicalUnits::bq;
      activity_json["total_uncertainty_ci"] = src.activityUncertainty / PhysicalUnits::ci;
      activity_json["total_uncertainty_str"] = PhysicalUnits::printToBestActivityUnits( src.activityUncertainty, 4, !use_bq );

      const double display_uncert_percent = 100.0 * src.traceSrcDisplayActUncertainty / src.traceSrcDisplayAct;
      activity_json["uncertainty_percent"] = display_uncert_percent;
    }//if( src.activityIsFit )

    // Add nuclide mass
    activity_json["nuclide_mass_grams"] = src.nuclideMass / PhysicalUnits::gram;
    activity_json["nuclide_mass_str"] = PhysicalUnits::printToBestMassUnits( src.nuclideMass, 4 );

    return;
  }//if( src.isTraceSource )


  // Handle Self-Attenuating/Intrinsic sources and Point sources
  // Both use the same activity fields, but fit status logic differs

  activity_json["bq"] = src.activity / PhysicalUnits::bq;
  activity_json["ci"] = src.activity / PhysicalUnits::ci;
  activity_json["str"] = PhysicalUnits::printToBestActivityUnits( src.activity, 4, !use_bq );

  // Determine if activity was fit (complex logic for self-atten sources)
  bool activityIsFit = src.activityIsFit;

  if( src.isSelfAttenSource )
  {
    // For self-atten sources, activity is considered fit if:
    // 1. Activity itself was fit, OR
    // 2. Mass fraction was fit (changes activity), OR
    // 3. Any shielding dimension was fit (volume changes, affects total activity)
    activityIsFit |= src.isSelfAttenVariableMassFrac;

    if( shield_details && (src.selfAttenShieldIndex < shield_details->size()) )
    {
      const GammaInteractionCalc::ShieldingDetails &shield = (*shield_details)[src.selfAttenShieldIndex];
      for( unsigned int i = 0; i < shield.m_num_dimensions; ++i )
        activityIsFit |= shield.m_fit_dimension[i];
    }//if( shield_details )
  }//if( src.isSelfAttenSource )

  activity_json["was_fit"] = activityIsFit;

  if( activityIsFit )
  {
    activity_json["uncertainty_bq"] = src.activityUncertainty / PhysicalUnits::bq;
    activity_json["uncertainty_ci"] = src.activityUncertainty / PhysicalUnits::ci;
    activity_json["uncertainty_str"] = PhysicalUnits::printToBestActivityUnits( src.activityUncertainty, 4, !use_bq );

    const double act_uncert_percent = 100.0 * src.activityUncertainty / src.activity;
    activity_json["uncertainty_percent"] = act_uncert_percent;
  }//if( activityIsFit )

  // Add nuclide mass
  activity_json["nuclide_mass_grams"] = src.nuclideMass / PhysicalUnits::gram;
  activity_json["nuclide_mass_str"] = PhysicalUnits::printToBestMassUnits( src.nuclideMass, 4 );
}//add_activity_to_json(...)


/** Add age information to JSON for LLM consumption.
 *
 * @param src Source details from fit results
 * @param age_json JSON object to populate with age fields
 */
void add_age_to_json(
  const GammaInteractionCalc::SourceDetails &src,
  nlohmann::json &age_json
)
{
  age_json["seconds"] = src.age / PhysicalUnits::second;
  age_json["str"] = PhysicalUnits::printToBestTimeUnits( src.age, 4 );
  age_json["was_fit"] = src.ageIsFit;
  age_json["is_fittable"] = src.ageIsFittable;

  if( src.ageIsFit )
  {
    age_json["uncertainty_seconds"] = src.ageUncertainty / PhysicalUnits::second;
    age_json["uncertainty_str"] = PhysicalUnits::printToBestTimeUnits( src.ageUncertainty, 4 );
  }//if( src.ageIsFit )

  if( src.ageDefiningNuc )
    age_json["defining_nuclide"] = src.ageDefiningNuc->symbol;
}//add_age_to_json(...)


/** Add peak information to JSON for LLM consumption.
 *
 * @param peak Peak details from fit results
 * @param peak_json JSON object to populate with peak fields
 */
void add_peak_info_to_json(
  const GammaInteractionCalc::PeakDetail &peak,
  nlohmann::json &peak_json
)
{
  peak_json["peak_energy_kev"] = peak.energy;
  peak_json["peak_fwhm_kev"] = peak.fwhm;

  peak_json["observed_counts"] = peak.observedCounts;
  peak_json["observed_counts_uncert"] = peak.observedUncert;
  peak_json["observed_cps"] = peak.cps;
  peak_json["predicted_counts"] = peak.expectedCounts;
  peak_json["num_sigma_off"] = peak.numSigmaOff;
  peak_json["observed_over_predicted"] = peak.observedOverExpected;
  peak_json["observed_over_predicted_uncert"] = peak.observedOverExpectedUncert;

  if( peak.backgroundCounts > 0.0f )
  {
    peak_json["background_counts"] = peak.backgroundCounts;
    if( peak.backgroundCountsUncert > 0.0f )
      peak_json["background_counts_uncert"] = peak.backgroundCountsUncert;
  }//if( peak.backgroundCounts > 0.0f )

  // Detector efficiency information
  peak_json["detector_solid_angle"] = peak.detSolidAngle;
  peak_json["detector_intrinsic_eff"] = peak.detIntrinsicEff;
  peak_json["detector_efficiency"] = peak.detEff;

  // Attenuation information
  auto &atten = peak_json["attenuation"];
  atten["total_factor"] = peak.m_totalAttenFactor;
  atten["air_factor"] = peak.m_airAttenFactor;
  atten["shielding_factor"] = peak.m_totalShieldAttenFactor;
  atten["by_shield"] = peak.m_attenuations;  // Array of attenuation by each shield
}//add_peak_info_to_json(...)


/** Add gamma contribution information to JSON for LLM consumption.
 *
 * @param gamma_src Gamma source contribution to peak
 * @param gamma_json JSON object to populate with gamma fields
 */
void add_gamma_contrib_to_json(
  const GammaInteractionCalc::PeakDetailSrc &gamma_src,
  nlohmann::json &gamma_json
)
{
  gamma_json["energy_kev"] = gamma_src.energy;
  gamma_json["branching_ratio"] = gamma_src.br;
  gamma_json["branching_ratio_str"] = SpecUtils::printCompact( gamma_src.br, 5 );

  gamma_json["source_photons_per_sec"] = gamma_src.cpsAtSource;
  gamma_json["source_photons_in_spectrum"] = gamma_src.countsAtSource;
  gamma_json["predicted_counts_contributed"] = gamma_src.modelContribToPeak;

  if( gamma_src.decayCorrection > 0.0 )
  {
    gamma_json["has_decay_correction"] = true;
    gamma_json["decay_correction"] = gamma_src.decayCorrection;
  }else
  {
    gamma_json["has_decay_correction"] = false;
  }
}//add_gamma_contrib_to_json(...)


/** Add shielding information to JSON for LLM consumption.
 *
 * @param shield Shielding details from fit results
 * @param shield_index Index of this shield in the shielding array
 * @param source_details Source details (needed to get activity for trace/self-atten sources), can be nullptr
 * @param use_bq If true, use Becquerel for display; if false, use Curies
 * @param shield_json JSON object to populate with shielding fields
 */
void add_shielding_to_json(
  const GammaInteractionCalc::ShieldingDetails &shield,
  const size_t shield_index,
  const std::vector<GammaInteractionCalc::SourceDetails> * const source_details,
  const bool use_bq,
  nlohmann::json &shield_json
)
{
  shield_json["name"] = shield.m_name;
  shield_json["shielding_number"] = static_cast<int>( shield_index );
  shield_json["is_generic"] = shield.m_is_generic;

  if( !shield.m_is_generic )
  {
    shield_json["chemical_formula"] = shield.m_chemical_formula;
    const double density = shield.m_density * PhysicalUnits::cm3 / PhysicalUnits::g;
    shield_json["density_g_per_cm3"] = density;
    shield_json["density_str"] = SpecUtils::printCompact( density, 5 ) + " g/cm³";
  }//if( !shield.m_is_generic )

  // Thickness/dimension information
  shield_json["thickness_cm"] = shield.m_thickness / PhysicalUnits::cm;
  shield_json["thickness_str"] = PhysicalUnits::printToBestLengthUnits( shield.m_thickness, 3 );

  shield_json["volume_cm3"] = shield.m_volume / PhysicalUnits::cm3;

  // Dimension-specific information based on geometry
  shield_json["num_dimensions"] = static_cast<int>( shield.m_num_dimensions );
  shield_json["geometry"] = GammaInteractionCalc::to_str( shield.m_geometry );

  // Check if any dimension was fit
  bool any_dim_fit = false;
  for( unsigned int i = 0; i < shield.m_num_dimensions; ++i )
    any_dim_fit |= shield.m_fit_dimension[i];

  shield_json["any_dimension_fit"] = any_dim_fit;

  if( any_dim_fit )
  {
    shield_json["volume_uncert_cm3"] = shield.m_volume_uncert / PhysicalUnits::cm3;
  }//if( any_dim_fit )

  // Trace sources in this shielding
  if( !shield.m_trace_sources.empty() && source_details )
  {
    auto &trace_srcs = shield_json["trace_sources"];
    trace_srcs = json::array();
    for( const auto &trace : shield.m_trace_sources )
    {
      json trace_json;
      trace_json["nuclide"] = trace.m_nuclide ? trace.m_nuclide->symbol : std::string("null");
      trace_json["trace_activity_type"] = GammaInteractionCalc::to_str( trace.m_trace_type );

      // Find corresponding source in source_details to get activity
      for( const auto &src : *source_details )
      {
        if( src.nuclide == trace.m_nuclide && src.isTraceSource && src.selfAttenShieldIndex == shield_index )
        {
          trace_json["total_activity_str"] = PhysicalUnits::printToBestActivityUnits( src.activity, 4, !use_bq );
          trace_json["display_activity_str"] = PhysicalUnits::printToBestActivityUnits( src.traceSrcDisplayAct, 4, !use_bq );
          trace_json["nuclide_mass_str"] = PhysicalUnits::printToBestMassUnits( src.nuclideMass, 4 );
          break;
        }//if( matching trace source )
      }//for( loop over source_details )

      if( trace.m_is_exp_dist )
      {
        trace_json["is_exponential_distribution"] = true;
        trace_json["relaxation_length_cm"] = trace.m_relaxation_length / PhysicalUnits::cm;
      }//if( trace.m_is_exp_dist )
      trace_srcs.push_back( trace_json );
    }//for( trace sources )
  }//if( trace sources )

  // Self-attenuating sources in this shielding
  if( !shield.m_mass_fractions.empty() && source_details )
  {
    auto &self_atten_srcs = shield_json["self_atten_sources"];
    self_atten_srcs = json::array();
    for( const auto &comp : shield.m_mass_fractions )
    {
      json comp_json;
      comp_json["nuclide"] = comp.m_nuclide ? comp.m_nuclide->symbol : std::string("null");
      comp_json["mass_fraction"] = comp.m_mass_frac;
      comp_json["mass_fraction_percent"] = comp.m_mass_frac * 100.0;
      comp_json["is_fit"] = comp.m_is_fit;
      if( comp.m_is_fit )
      {
        comp_json["mass_fraction_uncert"] = comp.m_mass_frac_uncert;
        comp_json["mass_fraction_uncert_percent"] = comp.m_mass_frac_uncert * 100.0;
      }//if( is_fit )

      // Find corresponding source in source_details to get activity and mass
      for( const auto &src : *source_details )
      {
        if( src.nuclide == comp.m_nuclide && src.isSelfAttenSource && src.selfAttenShieldIndex == shield_index )
        {
          comp_json["total_activity_str"] = PhysicalUnits::printToBestActivityUnits( src.activity, 4, !use_bq );
          comp_json["nuclide_mass_str"] = PhysicalUnits::printToBestMassUnits( src.nuclideMass, 4 );
          break;
        }//if( matching self-atten source )
      }//for( loop over source_details )

      self_atten_srcs.push_back( comp_json );
    }//for( mass fractions )
  }//if( mass fractions )
}//add_shielding_to_json(...)


/** Add fit configuration information to JSON for LLM consumption.
 *
 * @param results Fit results containing configuration
 * @param drf Detector response function
 * @param config_json JSON object to populate with configuration fields
 */
void add_fit_config_to_json(
  const ShieldingSourceFitCalc::ModelFitResults &results,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  nlohmann::json &config_json
)
{
  const bool fixedGeom = (drf && drf->isFixedGeometry());
  config_json["fixed_geometry_detector"] = fixedGeom;

  if( !fixedGeom )
  {
    config_json["distance_cm"] = results.distance / PhysicalUnits::cm;
    config_json["distance_str"] = PhysicalUnits::printToBestLengthUnits( results.distance, 3 );
    config_json["geometry"] = GammaInteractionCalc::to_str( results.geometry );
  }else
  {
    std::string geom_desc;
    switch( drf->geometryType() )
    {
      case DetectorPeakResponse::EffGeometryType::FarField:
        break;
      case DetectorPeakResponse::EffGeometryType::FixedGeomTotalAct:
        geom_desc = "total activity";
        break;
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerCm2:
        geom_desc = "activity per square centimeter";
        break;
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerM2:
        geom_desc = "activity per square meter";
        break;
      case DetectorPeakResponse::EffGeometryType::FixedGeomActPerGram:
        geom_desc = "activity per gram";
        break;
    }//switch( drf->geometryType() )

    config_json["fixed_geometry_type"] = geom_desc;
  }//if( !fixedGeom ) / else

  // Fit options
  auto &fit_options = config_json["fit_options"];
  fit_options["interference_correction"] = results.options.multiple_nucs_contribute_to_peaks;
  fit_options["attenuate_for_air"] = (results.options.attenuate_for_air && !fixedGeom);
  fit_options["decay_during_measurement"] = results.options.account_for_decay_during_meas;
  fit_options["photopeak_cluster_sigma"] = results.options.photopeak_cluster_sigma;
  fit_options["background_peak_subtract"] = results.options.background_peak_subtract;
  fit_options["element_nuclides_same_age"] = results.options.same_age_isotopes;
}//add_fit_config_to_json(...)


/** Comprehensive helper function to format ModelFitResults into JSON for LLM consumption.
 * This function uses the rich SourceDetails, PeakDetail, ShieldingDetails, and PeakResultPlotInfo
 * structures instead of the simpler SourceFitDef data.
 *
 * @param fit_results Pointer to ModelFitResults (can be null)
 * @param drf Detector response function (for geometry type)
 * @param use_bq If true, use Becquerel for display; if false, use Curies
 * @return JSON object with comprehensive fit results, or empty object if fit_results is null
 */
nlohmann::json fit_results_to_comprehensive_json(
  const std::shared_ptr<ShieldingSourceFitCalc::ModelFitResults> &fit_results,
  const std::shared_ptr<const DetectorPeakResponse> &drf,
  const bool use_bq
)
{
  json result;

  if( !fit_results )
    return result;

  // 1. Fit quality metrics
  auto &fit_quality = result["fit_quality"];
  fit_quality["chi2"] = fit_results->chi2;
  fit_quality["dof"] = fit_results->numDOF;
  fit_quality["chi2_per_dof"] = (fit_results->numDOF > 0)
                                  ? (fit_results->chi2 / fit_results->numDOF)
                                  : 0.0;
  fit_quality["num_peaks_used"] = fit_results->foreground_peaks.size();
  fit_quality["edm"] = fit_results->edm;
  fit_quality["num_fcn_calls"] = fit_results->num_fcn_calls;

  // 2. Fit configuration
  add_fit_config_to_json( *fit_results, drf, result["fit_configuration"] );

  // 3. Sources (using SourceDetails, not SourceFitDef)
  if( fit_results->source_calc_details )
  {
    result["sources"] = json::array();

    for( const auto &src : *fit_results->source_calc_details )
    {
      json src_json;

      src_json["nuclide"] = src.nuclide->symbol;
      src_json["half_life_str"] = PhysicalUnits::printToBestTimeUnits( src.nuclide->halfLife, 5 );

      // Determine source type
      if( src.isTraceSource )
        src_json["source_type"] = "Trace";
      else if( src.isSelfAttenSource )
        src_json["source_type"] = "SelfAttenuating";
      else
        src_json["source_type"] = "Point";

      // Add activity
      add_activity_to_json( src, use_bq, fit_results->shield_calc_details.get(), src_json["activity"] );

      // Add age
      add_age_to_json( src, src_json["age"] );

      // Add type-specific details
      if( src.isTraceSource || src.isSelfAttenSource )
      {
        src_json["shielding_name"] = src.selfAttenShieldName;
        src_json["shielding_index"] = static_cast<int>( src.selfAttenShieldIndex );

        if( src.isTraceSource )
        {
          auto &trace_details = src_json["trace_details"];
          trace_details["trace_activity_type"] = GammaInteractionCalc::to_str( src.traceActivityType );

          if( src.traceActivityType == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
          {
            trace_details["relaxation_length_cm"] = src.traceRelaxationLength / PhysicalUnits::cm;
            trace_details["relaxation_length_str"] = PhysicalUnits::printToBestLengthUnits( src.traceRelaxationLength, 4 );
          }//if( ExponentialDistribution )
        }//if( isTraceSource )

        if( src.isSelfAttenSource )
        {
          auto &self_atten = src_json["self_atten_details"];
          self_atten["mass_fraction"] = src.selfAttenMassFrac;
          self_atten["mass_fraction_percent"] = src.selfAttenMassFrac * 100.0;
          self_atten["mass_fraction_percent_str"] = SpecUtils::printCompact( 100.0 * src.selfAttenMassFrac, 5 ) + "%";
          self_atten["mass_fraction_was_fit"] = src.isSelfAttenVariableMassFrac;

          if( src.isSelfAttenVariableMassFrac )
          {
            self_atten["mass_fraction_uncert"] = src.selfAttenMassFracUncertainty;
            self_atten["mass_fraction_uncert_percent"] = src.selfAttenMassFracUncertainty * 100.0;
          }//if( mass fraction was fit )

        }//if( isSelfAttenSource )
      }//if( Trace or SelfAtten )

      // Add peaks this source contributes to
      if( fit_results->peak_calc_details )
      {
        src_json["peaks_this_source_contributes_to"] = json::array();

        for( const auto &peak : *fit_results->peak_calc_details )
        {
          // Check if this source contributes to this peak
          bool contributes = false;
          for( const auto &peak_src : peak.m_sources )
          {
            if( peak_src.nuclide == src.nuclide )
            {
              contributes = true;
              break;
            }//if( this source contributes )
          }//for( peak sources )

          if( !contributes )
            continue;

          json peak_json;
          add_peak_info_to_json( peak, peak_json );

          // Add gamma contributions from this source
          peak_json["gammas_from_this_source"] = json::array();
          for( const auto &peak_src : peak.m_sources )
          {
            if( peak_src.nuclide == src.nuclide )
            {
              json gamma_json;
              add_gamma_contrib_to_json( peak_src, gamma_json );
              peak_json["gammas_from_this_source"].push_back( gamma_json );
            }//if( from this source )
          }//for( gamma contributions )

          src_json["peaks_this_source_contributes_to"].push_back( peak_json );
        }//for( peaks )
      }//if( peak_calc_details )

      result["sources"].push_back( src_json );
    }//for( sources )
  }//if( source_calc_details )

  // 4. Shielding
  if( fit_results->shield_calc_details )
  {
    auto &shielding = result["shielding"];
    shielding["geometry"] = GammaInteractionCalc::to_str( fit_results->geometry );
    shielding["num_shieldings"] = static_cast<int>( fit_results->shield_calc_details->size() );

    // Dimension meanings based on geometry
    switch( fit_results->geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        shielding["dimension_meanings"] = json::array({"Radius"});
        shielding["num_dimensions"] = 1;
        break;
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        shielding["dimension_meanings"] = json::array({"Radius", "Length"});
        shielding["num_dimensions"] = 2;
        break;
      case GammaInteractionCalc::GeometryType::Rectangular:
        shielding["dimension_meanings"] = json::array({"Width", "Height", "Depth"});
        shielding["num_dimensions"] = 3;
        break;
      default:
        break;
    }//switch( geometry )

    shielding["shields"] = json::array();
    bool any_fit = false;

    for( size_t i = 0; i < fit_results->shield_calc_details->size(); ++i )
    {
      const auto &shield = (*fit_results->shield_calc_details)[i];
      json shield_json;

      add_shielding_to_json( shield, i, fit_results->source_calc_details.get(), use_bq, shield_json );

      // Check if any dimension was fit
      for( unsigned int d = 0; d < shield.m_num_dimensions; ++d )
        any_fit |= shield.m_fit_dimension[d];

      shielding["shields"].push_back( shield_json );
    }//for( shields )

    shielding["any_shielding_fit"] = any_fit;
  }//if( shield_calc_details )

  // 5. Peak comparison (solution_to_peak_comparison as requested)
  if( fit_results->peak_comparisons )
  {
    result["solution_to_peak_comparison"] = json::array();

    for( const auto &comp : *fit_results->peak_comparisons )
    {
      json comp_json;
      comp_json["energy_kev"] = comp.energy;
      comp_json["num_sigma_off"] = comp.numSigmaOff;
      comp_json["observed_over_predicted"] = comp.observedOverExpected;
      comp_json["observed_over_predicted_uncert"] = comp.observedOverExpectedUncert;

      result["solution_to_peak_comparison"].push_back( comp_json );
    }//for( peak comparisons )
  }//if( peak_comparisons )

  return result;
}//fit_results_to_comprehensive_json(...)


nlohmann::json executeActivityFit(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for activity_fit" );

  json result;
  result["status"] = "failed";

  // Check if GUI is already open
  ShieldingSourceDisplay *display = interspec->shieldingSourceFit();

  if( !display )
  {
    throw runtime_error( "Activity/Shielding fit GUI is not currently open. "
                       "Please open it manually from the Tools menu, or use 'activity_fit_one_off' for direct fitting without the GUI." );
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

    // Get detector from foreground spectrum (may be nullptr)
    std::shared_ptr<SpecMeas> meas = interspec->measurment( SpecUtils::SpectrumType::Foreground );
    std::shared_ptr<const DetectorPeakResponse> detector = meas ? meas->detector() : nullptr;

    // Get user preference for Bq vs Ci
    const bool use_bq = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec );

    // Format the results using the comprehensive JSON function
    const json formatted_results = fit_results_to_comprehensive_json( fit_results, detector, use_bq );

    // Merge the formatted results into our result object
    result["fit_quality"] = formatted_results["fit_quality"];
    result["fit_configuration"] = formatted_results["fit_configuration"];
    result["sources"] = formatted_results["sources"];
    result["shielding"] = formatted_results["shielding"];
    result["solution_to_peak_comparison"] = formatted_results["solution_to_peak_comparison"];

    // Keep legacy fields for backward compatibility
    result["chi2"] = formatted_results["fit_quality"]["chi2"];
    result["dof"] = formatted_results["fit_quality"]["dof"];
    result["chi2_per_dof"] = formatted_results["fit_quality"]["chi2_per_dof"];
    result["num_peaks_used"] = formatted_results["fit_quality"]["num_peaks_used"];
  }
  else
  {
    result["message"] = "Activity/Shielding fit started in GUI. Results will be displayed when fit completes.";
  }

  return result;
}//executeActivityFit(...)


nlohmann::json executeActivityFitOneOff(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required for activity_fit_one_off" );

  json result;
  result["status"] = "failed";

  // Build ShieldSourceInput and fit directly
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
  
  

  // Build source definitions from peaks or parse from parameters
  chi_input.config.sources.clear();

  if( params.contains("sources") && params["sources"].is_array() && !params["sources"].empty() )
  {
    // Parse sources from parameters and validate against peaks
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    if( !db )
      throw runtime_error( "Nuclide database not available" );

    for( const auto &source_json : params["sources"] )
    {
      // Required field: nuclide
      if( !source_json.contains("nuclide") || !source_json["nuclide"].is_string() )
        throw runtime_error( "Each source must have a 'nuclide' field (string)" );

      const string nuclide_name = source_json["nuclide"].get<string>();
      const SandiaDecay::Nuclide *nuc = db->nuclide( nuclide_name );

      if( !nuc )
        throw runtime_error( "Nuclide '" + nuclide_name + "' not found in database" );

      // Validate nuclide appears in fitting_peaks
      bool found_in_peaks = false;
      for( const auto &peak : fitting_peaks )
      {
        if( peak->parentNuclide() == nuc )
        {
          found_in_peaks = true;
          break;
        }
      }

      if( !found_in_peaks )
        throw runtime_error( "Source nuclide '" + nuclide_name + "' not found in fitting peaks. Sources must correspond to peaks." );

      // Build source definition
      ShieldingSourceFitCalc::SourceFitDef srcdef;
      srcdef.nuclide = nuc;
      srcdef.sourceType = ShieldingSourceFitCalc::ModelSourceType::Point;
      srcdef.ageDefiningNuc = nullptr; // Can be enhanced later for same_age_isotopes

      // Optional field: activity (string with units or number)
      if( source_json.contains("activity") )
      {
        if( source_json["activity"].is_string() )
        {
          const string activity_str = source_json["activity"].get<string>();
          srcdef.activity = parse_activity_string( activity_str );
        }
        else if( source_json["activity"].is_number() )
        {
          const double activity_val = source_json["activity"].get<double>();
          srcdef.activity = activity_val * PhysicalUnits::becquerel;
        }
        else
        {
          throw runtime_error( "Source 'activity' must be a string with units or a number" );
        }
      }
      else
      {
        srcdef.activity = 1.0E-6 * PhysicalUnits::curie; // Default: 1 uCi
      }

      // Optional field: fit_activity (boolean)
      if( source_json.contains("fit_activity") )
      {
        if( !source_json["fit_activity"].is_boolean() )
          throw runtime_error( "Source 'fit_activity' must be a boolean" );
        srcdef.fitActivity = source_json["fit_activity"].get<bool>();
      }
      else
      {
        srcdef.fitActivity = true; // Default: fit activity
      }

      // Optional field: age (string with units or number)
      if( source_json.contains("age") )
      {
        if( source_json["age"].is_string() )
        {
          const string age_str = source_json["age"].get<string>();
          try
          {
            srcdef.age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, nuc->halfLife );
          }
          catch( std::exception &e )
          {
            throw runtime_error( "Could not parse source age '" + age_str + "': " + e.what() );
          }
        }
        else if( source_json["age"].is_number() )
        {
          const double age_val = source_json["age"].get<double>();
          srcdef.age = age_val * PhysicalUnits::second;
        }
        else
        {
          throw runtime_error( "Source 'age' must be a string with units or a number" );
        }
      }
      else
      {
        srcdef.age = PeakDef::defaultDecayTime( nuc ); // Default age
      }

      // Optional field: fit_age (boolean)
      if( source_json.contains("fit_age") )
      {
        if( !source_json["fit_age"].is_boolean() )
          throw runtime_error( "Source 'fit_age' must be a boolean" );
        srcdef.fitAge = source_json["fit_age"].get<bool>();
      }
      else
      {
        srcdef.fitAge = false; // Default: don't fit age
      }

      chi_input.config.sources.push_back( srcdef );
    }
  }
  else
  {
    // Default behavior: auto-generate sources from peaks (backward compatibility)
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
    for( auto &pair : sources_map )
      chi_input.config.sources.push_back( pair.second );
  }

  // Parse shielding from parameters or clear for no shielding
  chi_input.config.shieldings.clear();

  if( params.contains("shielding") && params["shielding"].is_array() )
  {
    MaterialDB *materialDb = interspec->materialDataBase();
    if( !materialDb )
      throw runtime_error( "Material database not available" );

    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    if( !db )
      throw runtime_error( "Nuclide database not available" );

    for( const auto &shield_json : params["shielding"] )
    {
      // Required field: material
      if( !shield_json.contains("material") || !shield_json["material"].is_string() )
        throw runtime_error( "Each shielding must have a 'material' field (string)" );

      const string material_name = shield_json["material"].get<string>();
      const Material *mat = materialDb->material( material_name );

      if( !mat )
        throw runtime_error( "Material '" + material_name + "' not found in database" );

      // Create ShieldingInfo object
      ShieldingSourceFitCalc::ShieldingInfo shielding;
      shielding.m_geometry = chi_input.config.geometry;
      shielding.m_material = std::make_shared<const Material>( *mat );
      shielding.m_isGenericMaterial = false;
      shielding.m_forFitting = true;

      // Parse dimensions based on geometry
      switch( chi_input.config.geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
        {
          // Field: radial_thickness or thickness (string with units)
          string thickness_str;
          if( shield_json.contains("radial_thickness") && shield_json["radial_thickness"].is_string() )
            thickness_str = shield_json["radial_thickness"].get<string>();
          else if( shield_json.contains("thickness") && shield_json["thickness"].is_string() )
            thickness_str = shield_json["thickness"].get<string>();

          if( !thickness_str.empty() )
            shielding.m_dimensions[0] = parse_distance_string( thickness_str );
          else
            shielding.m_dimensions[0] = 1.0 * PhysicalUnits::cm; // Default: 1 cm

          // Fit flag: default true if dimension not specified
          if( shield_json.contains("fit_thickness") )
          {
            if( !shield_json["fit_thickness"].is_boolean() )
              throw runtime_error( "Shielding 'fit_thickness' must be a boolean" );
            shielding.m_fitDimensions[0] = shield_json["fit_thickness"].get<bool>();
          }
          else
          {
            shielding.m_fitDimensions[0] = thickness_str.empty();
          }
          break;
        }

        case GammaInteractionCalc::GeometryType::CylinderSideOn:
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        {
          // Field: radius (string with units) - this is the radius thickness
          string radius_str;
          if( shield_json.contains("radius") && shield_json["radius"].is_string() )
            radius_str = shield_json["radius"].get<string>();

          if( !radius_str.empty() )
            shielding.m_dimensions[0] = parse_distance_string( radius_str );
          else
            shielding.m_dimensions[0] = 5.0 * PhysicalUnits::cm; // Default: 5 cm

          // Field: length (string with units) - this is the length thickness
          // For first (inner) cylinder, this is half-length
          // For subsequent layers, this is how much longer on each side (half of total increase)
          string length_str;
          if( shield_json.contains("length") && shield_json["length"].is_string() )
            length_str = shield_json["length"].get<string>();

          if( !length_str.empty() )
            shielding.m_dimensions[1] = parse_distance_string( length_str );
          else
            shielding.m_dimensions[1] = 10.0 * PhysicalUnits::cm; // Default: 10 cm

          // Fit flags: default true if dimensions not specified
          if( shield_json.contains("fit_radius") )
          {
            if( !shield_json["fit_radius"].is_boolean() )
              throw runtime_error( "Shielding 'fit_radius' must be a boolean" );
            shielding.m_fitDimensions[0] = shield_json["fit_radius"].get<bool>();
          }
          else
          {
            shielding.m_fitDimensions[0] = radius_str.empty();
          }

          if( shield_json.contains("fit_length") )
          {
            if( !shield_json["fit_length"].is_boolean() )
              throw runtime_error( "Shielding 'fit_length' must be a boolean" );
            shielding.m_fitDimensions[1] = shield_json["fit_length"].get<bool>();
          }
          else
          {
            shielding.m_fitDimensions[1] = length_str.empty();
          }
          break;
        }

        case GammaInteractionCalc::GeometryType::Rectangular:
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          throw runtime_error( "Rectangular geometry not yet supported for shielding" );
      }

      // Parse trace sources (nested array)
      if( shield_json.contains("trace_sources") && shield_json["trace_sources"].is_array() )
      {
        for( const auto &trace_json : shield_json["trace_sources"] )
        {
          // Required field: nuclide
          if( !trace_json.contains("nuclide") || !trace_json["nuclide"].is_string() )
            throw runtime_error( "Each trace source must have a 'nuclide' field (string)" );

          const string trace_nuclide_name = trace_json["nuclide"].get<string>();
          const SandiaDecay::Nuclide *trace_nuc = db->nuclide( trace_nuclide_name );

          if( !trace_nuc )
            throw runtime_error( "Trace source nuclide '" + trace_nuclide_name + "' not found in database" );

          // Validate nuclide appears in fitting_peaks
          bool found_in_peaks = false;
          for( const auto &peak : fitting_peaks )
          {
            if( peak->parentNuclide() == trace_nuc )
            {
              found_in_peaks = true;
              break;
            }
          }

          if( !found_in_peaks )
            throw runtime_error( "Trace source nuclide '" + trace_nuclide_name + "' not found in fitting peaks" );

          // Create TraceSourceInfo
          ShieldingSourceFitCalc::TraceSourceInfo trace_info;
          trace_info.m_nuclide = trace_nuc;

          // Optional field: fit_activity (boolean)
          if( trace_json.contains("fit_activity") )
          {
            if( !trace_json["fit_activity"].is_boolean() )
              throw runtime_error( "Trace source 'fit_activity' must be a boolean" );
            trace_info.m_fitActivity = trace_json["fit_activity"].get<bool>();
          }
          else
          {
            trace_info.m_fitActivity = true; // Default: fit activity
          }

          // Optional field: activity_type (string)
          // If not specified explicitly, will be inferred from activity string units
          GammaInteractionCalc::TraceActivityType activity_type = GammaInteractionCalc::TraceActivityType::TotalActivity;
          bool activity_type_explicit = false;

          if( trace_json.contains("activity_type") && trace_json["activity_type"].is_string() )
          {
            const string type_str = trace_json["activity_type"].get<string>();
            activity_type_explicit = true;

            if( type_str == "TotalActivity" )
              activity_type = GammaInteractionCalc::TraceActivityType::TotalActivity;
            else if( type_str == "ActivityPerCm3" || type_str == "ActivityPerVolume" )
              activity_type = GammaInteractionCalc::TraceActivityType::ActivityPerCm3;
            else if( type_str == "ActivityPerGram" || type_str == "ActivityPerMass" )
              activity_type = GammaInteractionCalc::TraceActivityType::ActivityPerGram;
            else if( type_str == "ExponentialDistribution" )
              activity_type = GammaInteractionCalc::TraceActivityType::ExponentialDistribution;
            else
              throw runtime_error( "Invalid trace activity_type: '" + type_str + "'. Must be TotalActivity, ActivityPerCm3, ActivityPerGram, or ExponentialDistribution" );
          }

          // Optional field: activity (string with units or number)
          // If activity contains "/cm3", "/g", etc., it will override activity_type
          if( trace_json.contains("activity") )
          {
            if( trace_json["activity"].is_string() )
            {
              const string activity_str = trace_json["activity"].get<string>();

              // Parse activity string which may contain per-volume or per-mass units
              std::pair<GammaInteractionCalc::TraceActivityType, double> parsed = parse_trace_activity_string( activity_str );

              // If activity string contains per-unit specifier, use it (unless ExponentialDistribution was explicitly set)
              if( !activity_type_explicit || activity_type != GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
              {
                // Only override if parsed type is not TotalActivity or if no explicit type was given
                if( parsed.first != GammaInteractionCalc::TraceActivityType::TotalActivity || !activity_type_explicit )
                {
                  activity_type = parsed.first;
                }
              }

              trace_info.m_activity = parsed.second;
            }
            else if( trace_json["activity"].is_number() )
            {
              const double activity_val = trace_json["activity"].get<double>();
              trace_info.m_activity = activity_val * PhysicalUnits::becquerel;
            }
            else
            {
              throw runtime_error( "Trace source 'activity' must be a string with units or a number" );
            }
          }
          else
          {
            // Defaults by type
            switch( activity_type )
            {
              case GammaInteractionCalc::TraceActivityType::TotalActivity:
                trace_info.m_activity = 1.0E-6 * PhysicalUnits::curie; // 1 uCi
                break;
              case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
              case GammaInteractionCalc::TraceActivityType::ExponentialDistribution:
                trace_info.m_activity = 1.0 * PhysicalUnits::bq / PhysicalUnits::cm3;
                break;
              case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
                trace_info.m_activity = 1.0 * PhysicalUnits::bq / PhysicalUnits::gram;
                break;
              case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
                break;
            }
          }

          trace_info.m_type = activity_type;

          // Conditional field: relaxation_length (only for ExponentialDistribution)
          if( activity_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
          {
            if( trace_json.contains("relaxation_length") && trace_json["relaxation_length"].is_string() )
            {
              const string relax_str = trace_json["relaxation_length"].get<string>();
              trace_info.m_relaxationDistance = parse_distance_string( relax_str );
            }
            else
            {
              trace_info.m_relaxationDistance = 1.0 * PhysicalUnits::cm; // Default: 1 cm
            }
          }

          shielding.m_traceSources.push_back( trace_info );
        }
      }

      chi_input.config.shieldings.push_back( shielding );
    }
  }

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

  // Format results using the comprehensive JSON function
  result["status"] = "success";

  // Get user preference for Bq vs Ci
  InterSpec *interspec_inst = InterSpec::instance();
  const bool use_bq = !interspec_inst ? true
                      : !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", interspec_inst );

  // Create comprehensive JSON
  const json fit_json = fit_results_to_comprehensive_json( fit_results, detector, use_bq );

  // Merge into result
  result["fit_quality"] = fit_json["fit_quality"];
  result["fit_configuration"] = fit_json["fit_configuration"];
  result["sources"] = fit_json["sources"];
  result["shielding"] = fit_json["shielding"];
  result["solution_to_peak_comparison"] = fit_json["solution_to_peak_comparison"];

  // Keep legacy fields for backward compatibility
  result["chi2"] = fit_json["fit_quality"]["chi2"];
  result["dof"] = fit_json["fit_quality"]["dof"];
  result["chi2_per_dof"] = fit_json["fit_quality"]["chi2_per_dof"];
  result["num_peaks_used"] = fit_json["fit_quality"]["num_peaks_used"];

  result["message"] = "Fit completed successfully";

  return result;
}//executeActivityFitOneOff(...)

}//namespace ActivityFitTool
}//namespace LlmTools

#endif //USE_LLM_INTERFACE
