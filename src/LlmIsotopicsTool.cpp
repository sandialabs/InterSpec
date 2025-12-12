#include "InterSpec_config.h"
#include "InterSpec/LlmIsotopicsTool.h"

#if( USE_LLM_INTERFACE )

#include <algorithm>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "InterSpec/InterSpec.h"
#include "InterSpec/LlmConfig.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/RelActAutoGui.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcAuto.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

using namespace std;
using json = nlohmann::json;

namespace LlmTools {
namespace IsotopicsTool {

// ============================================================================
// Isotopics Tools Helper Functions
// ============================================================================

namespace {

/** Helper function to convert RelActAutoGuiState to a summary JSON representation.

 This provides a brief overview suitable for listing presets or showing what was loaded.

 @param state The state to convert
 @returns JSON object with summary information (note, description, config counts, basic options)
 */
nlohmann::json stateToJsonSummary( const RelActCalcAuto::RelActAutoGuiState &state )
{
  using namespace std;
  using json = nlohmann::json;

  json result;

  // Basic metadata
  if( !state.note.empty() )
    result["note"] = state.note;

  if( !state.description.empty() )
    result["description"] = state.description;

  // Configuration summary
  const RelActCalcAuto::Options &options = state.options;

  result["roi_count"] = options.rois.size();

  // Count nuclides across all rel eff curves
  size_t total_nuclides = 0;
  vector<string> nuclide_names;
  for( const RelActCalcAuto::RelEffCurveInput &curve : options.rel_eff_curves )
  {
    for( const RelActCalcAuto::NucInputInfo &nuc : curve.nuclides )
    {
      total_nuclides++;
      nuclide_names.push_back( nuc.name() );
    }
  }
  result["total_nuclides"] = total_nuclides;
  result["nuclides"] = nuclide_names;

  result["rel_eff_curve_count"] = options.rel_eff_curves.size();

  // Basic options
  result["fit_energy_cal"] = options.fit_energy_cal;
  result["fwhm_form"] = RelActCalcAuto::to_str(options.fwhm_form);
  result["background_subtract"] = state.background_subtract;

  return result;
}//stateToJsonSummary()


/** Helper function to convert RelActAutoGuiState to a detailed JSON representation.

 This provides complete configuration details suitable for reviewing the full configuration.

 @param state The state to convert
 @returns JSON object with detailed configuration information
 */
nlohmann::json stateToJsonDetailed( const RelActCalcAuto::RelActAutoGuiState &state )
{
  using namespace std;
  using json = nlohmann::json;

  json result;

  // Basic metadata
  if( !state.note.empty() )
    result["note"] = state.note;

  if( !state.description.empty() )
    result["description"] = state.description;

  const RelActCalcAuto::Options &opts = state.options;

  // ROIs with full details
  result["rois"] = json::array();
  for( const auto &roi : opts.rois )
  {
    json roi_json;
    roi_json["lower_energy"] = roi.lower_energy;
    roi_json["upper_energy"] = roi.upper_energy;
    roi_json["continuum_type"] = PeakContinuum::offset_type_label_tr(roi.continuum_type);

    const char *range_type_str = RelActCalcAuto::RoiRange::to_str( roi.range_limits_type );
    roi_json["range_type"] = range_type_str;

    result["rois"].push_back( roi_json );
  }

  // Lambda to create JSON for a rel eff curve
  auto curveToJson = [](const RelActCalcAuto::RelEffCurveInput &curve, bool include_index, size_t index) -> json {
    json curve_json;
    
    curve_json["rel_eff_eqn_type"] = RelActCalc::to_str( curve.rel_eff_eqn_type );
    curve_json["nucs_of_el_same_age"] = curve.nucs_of_el_same_age;

    const bool is_physical = (curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel);

    if( is_physical )
    {
      // Physical model - include shielding info and Hoerl option
      curve_json["use_hoerl"] = curve.phys_model_use_hoerl;

      if( curve.phys_model_self_atten )
      {
        json self_atten;
        if( curve.phys_model_self_atten->material )
          self_atten["material"] = curve.phys_model_self_atten->material->name;
        else
          self_atten["atomic_number"] = curve.phys_model_self_atten->atomic_number;
        self_atten["areal_density"] = curve.phys_model_self_atten->areal_density;
        self_atten["fit_atomic_number"] = curve.phys_model_self_atten->fit_atomic_number;
        self_atten["fit_areal_density"] = curve.phys_model_self_atten->fit_areal_density;
        curve_json["self_attenuation"] = self_atten;
      }

      if( !curve.phys_model_external_atten.empty() )
      {
        json ext_atten_array = json::array();
        for( const auto &ext : curve.phys_model_external_atten )
        {
          if( ext )
          {
            json ext_json;
            if( ext->material )
              ext_json["material"] = ext->material->name;
            else
              ext_json["atomic_number"] = ext->atomic_number;
            ext_json["areal_density"] = ext->areal_density;
            ext_json["fit_atomic_number"] = ext->fit_atomic_number;
            ext_json["fit_areal_density"] = ext->fit_areal_density;
            ext_atten_array.push_back( ext_json );
          }
        }
        if( !ext_atten_array.empty() )
          curve_json["external_attenuation"] = ext_atten_array;
      }
    }
    else
    {
      // Non-physical model - include equation order
      curve_json["rel_eff_eqn_order"] = curve.rel_eff_eqn_order;
    }

    if( include_index )
    {
      curve_json["index"] = index;
      curve_json["name"] = curve.name;
    }

    return curve_json;
  };

  const bool multiple_curves = opts.rel_eff_curves.size() > 1;

  // Rel eff curve(s) details
  if( multiple_curves )
  {
    // Multiple curves - use array with indices
    result["rel_eff_curves"] = json::array();
    for( size_t i = 0; i < opts.rel_eff_curves.size(); ++i )
      result["rel_eff_curves"].push_back( curveToJson(opts.rel_eff_curves[i], true, i) );
  }
  else if( !opts.rel_eff_curves.empty() )
  {
    // Single curve - use object (not array), no index needed
    result["rel_eff_curve"] = curveToJson( opts.rel_eff_curves[0], false, 0 );
  }

  // Nuclides with full details (from all rel eff curves)
  result["nuclides"] = json::array();
  for( size_t curve_idx = 0; curve_idx < opts.rel_eff_curves.size(); ++curve_idx )
  {
    const RelActCalcAuto::RelEffCurveInput &curve = opts.rel_eff_curves[curve_idx];
    for( const RelActCalcAuto::NucInputInfo &nuc_input : curve.nuclides )
    {
      json nuc_json;
      nuc_json["name"] = nuc_input.name();

      // Only include curve info if there are multiple curves
      if( multiple_curves )
      {
        nuc_json["rel_eff_curve"] = curve.name;
        nuc_json["rel_eff_index"] = curve_idx;
      }

      if( nuc_input.age >= 0.0 )
      {
        nuc_json["age_days"] = nuc_input.age / PhysicalUnits::day;
        nuc_json["age"] = PhysicalUnits::printToBestTimeUnits( nuc_input.age, 6 );
        nuc_json["fit_age"] = nuc_input.fit_age;
      }

      result["nuclides"].push_back( nuc_json );
    }
  }

  // All options
  result["fit_energy_cal"] = opts.fit_energy_cal;
  result["fwhm_form"] = RelActCalcAuto::to_str(opts.fwhm_form);
  result["fwhm_estimation_method"] = RelActCalcAuto::to_str(opts.fwhm_estimation_method);
  result["skew_type"] = PeakDef::to_string(opts.skew_type);
  result["additional_br_uncert"] = opts.additional_br_uncert;
  result["background_subtract"] = state.background_subtract;
  result["show_ref_lines"] = state.show_ref_lines;

  if( state.lower_display_energy < state.upper_display_energy )
  {
    result["lower_display_energy"] = state.lower_display_energy;
    result["upper_display_energy"] = state.upper_display_energy;
  }

  if( !opts.spectrum_title.empty() )
    result["spectrum_title"] = opts.spectrum_title;

  return result;
}//stateToJsonDetailed()


/** Helper function to extract metadata from RelActAutoGuiState for listing purposes.

 @param state The state to extract info from
 @returns JSON object with metadata (description, sources, energy_range, rois, rel_eff_form)
 */
nlohmann::json extractStateMetadata( const RelActCalcAuto::RelActAutoGuiState &state )
{
  using namespace std;
  using json = nlohmann::json;

  json metadata;

  // Extract description
  if( !state.description.empty() )
    metadata["description"] = state.description;

  // Extract sources (nuclides) from all rel eff curves - remove duplicates
  set<string> unique_sources;
  for( const auto &curve : state.options.rel_eff_curves )
  {
    for( const auto &nuc_info : curve.nuclides )
    {
      // Extract source name from variant
      if( const SandiaDecay::Nuclide * const* nuc = std::get_if<const SandiaDecay::Nuclide*>(&nuc_info.source) )
      {
        if( *nuc )
          unique_sources.insert( (*nuc)->symbol );
      }
      else if( const SandiaDecay::Element * const* el = std::get_if<const SandiaDecay::Element*>(&nuc_info.source) )
      {
        if( *el )
          unique_sources.insert( (*el)->symbol );
      }
      else if( const ReactionGamma::Reaction * const* rxn = std::get_if<const ReactionGamma::Reaction*>(&nuc_info.source) )
      {
        if( *rxn )
          unique_sources.insert( (*rxn)->name() );
      }
    }
  }

  if( !unique_sources.empty() )
  {
    json sources_array = json::array();
    for( const string &source : unique_sources )
      sources_array.push_back( source );
    metadata["sources"] = sources_array;
  }

  // Extract energy range from ROIs
  double min_energy = std::numeric_limits<double>::max();
  double max_energy = std::numeric_limits<double>::lowest();

  for( const auto &roi : state.options.rois )
  {
    min_energy = std::min( min_energy, roi.lower_energy );
    max_energy = std::max( max_energy, roi.upper_energy );
  }

  if( min_energy < max_energy && !state.options.rois.empty() )
  {
    metadata["energy_range"] = SpecUtils::printCompact( min_energy, 5 ) + " - "
                              + SpecUtils::printCompact( max_energy, 5 ) + " keV";
  }

  // Extract ROIs list
  json rois_array = json::array();
  for( const auto &roi : state.options.rois )
  {
    json roi_obj;
    roi_obj["lower_energy"] = roi.lower_energy;
    roi_obj["upper_energy"] = roi.upper_energy;
    rois_array.push_back( roi_obj );
  }
  if( !rois_array.empty() )
    metadata["rois"] = rois_array;

  // Extract relative efficiency form(s)
  if( state.options.rel_eff_curves.size() == 1 )
  {
    // Single curve - return as string
    metadata["rel_eff_form"] = RelActCalc::to_str( state.options.rel_eff_curves[0].rel_eff_eqn_type );
  }
  else if( state.options.rel_eff_curves.size() > 1 )
  {
    // Multiple curves - return as array
    json rel_eff_forms_array = json::array();
    for( const auto &curve : state.options.rel_eff_curves )
      rel_eff_forms_array.push_back( RelActCalc::to_str( curve.rel_eff_eqn_type ) );
    metadata["rel_eff_forms"] = rel_eff_forms_array;
  }

  return metadata;
}//extractStateMetadata()


/** Helper function to get the current RelActAutoGuiState from SpecMeas XML.

 Returns a copy of the state, so modifications to the returned object will not affect
 the global state unless you explicitly save it back using saveRelActState().

 @param interspec InterSpec instance
 @param create_if_missing If true, creates a new empty state if none exists.
                          If false and no state exists, throws std::runtime_error.
 @returns RelActAutoGuiState by value (a copy of the current state)
 @throws std::runtime_error if interspec is null, no foreground spectrum loaded,
                            or if create_if_missing is false and no state exists
 */
RelActCalcAuto::RelActAutoGuiState getOrCreateRelActState(
  InterSpec* interspec,
  bool create_if_missing = true
)
{
  if( !interspec )
    throw std::runtime_error( "InterSpec instance required" );

  RelActAutoGui *gui = interspec->relActAutoWindow(false);
  if( gui )
  {
    try
    {
      RelActCalcAuto::RelActAutoGuiState state;
      gui->serialize( state );
      return state;
    }catch( std::exception &e )
    {
      cerr << "Failed to serialize RelActAutoGui to a state!" << endl;
    }
  }

  std::shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw std::runtime_error( "No foreground spectrum loaded" );

  MaterialDB *materialDb = interspec->materialDataBase();
  std::unique_ptr<RelActCalcAuto::RelActAutoGuiState> state_ptr = spec->getRelActAutoGuiState( materialDb );

  if( !state_ptr )
  {
    if( create_if_missing )
      return RelActCalcAuto::RelActAutoGuiState();
    else
      throw std::runtime_error( "No isotopics configuration exists" );
  }

  // Return by value (move from unique_ptr)
  return std::move( *state_ptr );
}//getOrCreateRelActState()


/** Helper function to save RelActAutoGuiState to SpecMeas as XML.

 @param interspec InterSpec instance
 @param state The state to save
 @throws std::runtime_error if interspec is null or no foreground spectrum loaded
 */
void saveRelActState(
  InterSpec* interspec,
  const RelActCalcAuto::RelActAutoGuiState &state
)
{
  if( !interspec )
    throw std::runtime_error( "InterSpec instance required" );

  std::shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw std::runtime_error( "No foreground spectrum loaded" );

  // Use new convenience setter
  spec->setRelActAutoGuiState( &state );

  RelActAutoGui *gui = interspec->relActAutoWindow(false);
  if( gui )
    gui->deSerialize( state );
}//saveRelActState()

}//namespace (anonymous)


// ============================================================================
// Isotopics Discovery and State Management Tools
// ============================================================================

nlohmann::json executeListIsotopicsPresets(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  using namespace std;

  json result;
  result["presets"] = json::array();

  // Check if there's a current configuration
  try
  {
    RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

    json current_preset;
    current_preset["name"] = "Current Configuration";
    current_preset["is_current"] = true;

    // Extract metadata from current state
    json metadata = extractStateMetadata( state );

    // Append description to base string if present
    string base_description = "Current isotopics configuration in memory";
    if( metadata.contains("description") && !metadata["description"].get<string>().empty() )
    {
      base_description += ": " + metadata["description"].get<string>();
      metadata.erase("description"); // Remove from metadata since we've merged it
    }
    current_preset["description"] = base_description;

    // Merge remaining metadata fields
    current_preset.update(metadata);

    result["presets"].push_back( current_preset );
  }catch( ... )
  {
    // No current state - that's fine (getOrCreateRelActState throws if create_if_missing=false and no state exists)
  }

  // Scan both writable and static data/rel_act directories for preset files
  const string writable_data_dir = InterSpec::writableDataDirectory();
  const string static_data_dir = InterSpec::staticDataDirectory();

  const string writable_rel_act_dir = SpecUtils::append_path( writable_data_dir, "rel_act" );
  const string static_rel_act_dir = SpecUtils::append_path( static_data_dir, "rel_act" );

  // Collect preset files from both directories
  vector<string> preset_files;

  // Scan writable directory first (user-created presets)
  if( SpecUtils::is_directory(writable_rel_act_dir) )
  {
    vector<string> writable_files = SpecUtils::recursive_ls( writable_rel_act_dir, ".xml" );
    preset_files.insert( preset_files.end(), writable_files.begin(), writable_files.end() );
  }

  // Scan static directory (built-in presets)
  if( SpecUtils::is_directory(static_rel_act_dir) )
  {
    vector<string> static_files = SpecUtils::recursive_ls( static_rel_act_dir, ".xml" );
    preset_files.insert( preset_files.end(), static_files.begin(), static_files.end() );
  }

  if( preset_files.empty() )
  {
    result["warning"] = "No isotopics preset files found in data directories";
    return result;
  }

  // Filter out unwanted presets
  const vector<string> exclude_keywords = { "U inside U", "multiple U" };

  for( const string &filepath : preset_files )
  {
    const string filename = SpecUtils::filename(filepath);

    // Check if filename contains any exclude keywords (case-insensitive)
    bool should_exclude = false;
    for( const string &keyword : exclude_keywords )
    {
      if( SpecUtils::icontains(filename, keyword) )
      {
        should_exclude = true;
        break;
      }
    }

    if( should_exclude )
      continue;

    json preset;
    // Remove .xml extension from name for cleaner display
    string preset_name = filename;
    if( SpecUtils::iends_with(preset_name, ".xml") )
      preset_name = preset_name.substr(0, preset_name.size() - 4);
    preset["name"] = preset_name;
    preset["is_current"] = false;
    // We wont include the full path to the XML file, as the LLM doesnt need to know that
    //preset["path"] = filepath;

    // Indicate whether this is a user preset or built-in preset
    if( SpecUtils::istarts_with(filepath, writable_rel_act_dir) )
      preset["location"] = "user";
    else
      preset["location"] = "built-in";

    // Try to extract basic info from filename
    if( SpecUtils::icontains(filename, "Pu") )
      preset["material_type"] = "Plutonium";
    else if( SpecUtils::icontains(filename, "U") )
      preset["material_type"] = "Uranium";

    // Try to parse the preset file to extract metadata
    try
    {
      // Parse XML file
      vector<char> xml_data;
      SpecUtils::load_file_data( filepath.c_str(), xml_data );

      rapidxml::xml_document<char> doc;
      const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
      doc.parse<flags>( &xml_data[0] );

      const rapidxml::xml_node<char> *base_node = doc.first_node( "RelActCalcAuto" );
      if( base_node )
      {
        // Deserialize as RelActAutoGuiState
        RelActCalcAuto::RelActAutoGuiState file_state;
        MaterialDB *materialDb = interspec ? interspec->materialDataBase() : nullptr;
        file_state.deSerialize( base_node, materialDb );

        // Extract metadata
        json metadata = extractStateMetadata( file_state );

        // Merge metadata into preset
        preset.update( metadata );
      }
    }
    catch( std::exception &e )
    {
      // Failed to parse - that's okay, just skip the metadata extraction
      preset["parse_warning"] = string("Could not extract metadata: ") + e.what();
    }

    result["presets"].push_back( preset );
  }

  result["count"] = result["presets"].size();
  return result;
}//executeListIsotopicsPresets()


nlohmann::json executeGetIsotopicsConfig(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  using namespace std;

  try
  {
    RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

    json result;
    result["has_config"] = true;

    // Use the detailed conversion helper
    json details = stateToJsonDetailed( state );
    result.update( details );

    return result;
  }catch( std::exception &e )
  {
    // No configuration exists
    json result;
    result["has_config"] = false;
    result["message"] = "No isotopics configuration exists";
    return result;
  }
}//executeGetIsotopicsConfig()


nlohmann::json executeResetIsotopicsConfig(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw std::runtime_error( "InterSpec instance required" );

  std::shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw std::runtime_error( "No foreground spectrum loaded" );

  // Clear the state
  spec->setRelActAutoGuiState( nullptr );
  
  RelActAutoGui *gui = interspec->relActAutoWindow(false);
  if( gui )
  {
    try
    {
      RelActCalcAuto::RelActAutoGuiState empty_state;
      gui->deSerialize( empty_state );
    }catch( std::exception &e )
    {
      cerr << "Failed to deSerialize an empty state to RelActAutoGui!  " << e.what() << endl;
    }
  }

  json result;
  result["success"] = true;
  result["message"] = "Isotopics configuration cleared";

  return result;
}//executeResetIsotopicsConfig()


nlohmann::json executeLoadIsotopicsPreset(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  using namespace std;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  string preset_name = params.value( "preset", string() );

  if( preset_name.empty() )
    throw runtime_error( "preset parameter is required" );

  json result;

  // Extract just the filename (security: prevent directory traversal)
  preset_name = SpecUtils::filename( preset_name );

  if( preset_name.empty() )
    throw runtime_error( "Invalid preset name" );

  // Add .xml extension if not present
  if( !SpecUtils::iends_with(preset_name, ".xml") )
    preset_name += ".xml";

  // Search for the preset file in writable directory first, then static directory
  const string writable_data_dir = InterSpec::writableDataDirectory();
  const string static_data_dir = InterSpec::staticDataDirectory();

  const string writable_rel_act_dir = SpecUtils::append_path( writable_data_dir, "rel_act" );
  const string static_rel_act_dir = SpecUtils::append_path( static_data_dir, "rel_act" );

  string preset_path;

  // Check writable directory first
  string candidate = SpecUtils::append_path( writable_rel_act_dir, preset_name );
  if( SpecUtils::is_file(candidate) )
  {
    preset_path = candidate;
  }
  else
  {
    // Check static directory
    candidate = SpecUtils::append_path( static_rel_act_dir, preset_name );
    if( SpecUtils::is_file(candidate) )
    {
      preset_path = candidate;
    }
  }

  if( preset_path.empty() )
    throw runtime_error( "Could not find preset: " + preset_name );

  // Load and parse the XML file
  vector<char> xml_data;
  SpecUtils::load_file_data( preset_path.c_str(), xml_data );

  if( xml_data.size() < 10 )
    throw runtime_error( "Failed to read preset file: " + preset_path );

  // Parse XML into a new document
  auto xml_doc = make_unique<rapidxml::xml_document<char>>();

  try
  {
    const int flags = rapidxml::parse_normalize_whitespace | rapidxml::parse_trim_whitespace;
    xml_doc->parse<flags>( &(xml_data[0]) );
  }catch( std::exception &e )
  {
    throw runtime_error( "Failed to parse preset XML: " + string(e.what()) );
  }

  const rapidxml::xml_node<char> *base_node = xml_doc->first_node( "RelActCalcAuto" );
  if( !base_node )
    throw runtime_error( "Invalid preset file - missing RelActCalcAuto root node" );

  // Parse the state to provide summary (after we've saved the XML)
  RelActCalcAuto::RelActAutoGuiState state;
  MaterialDB *materialDb = interspec->materialDataBase();

  try
  {
    state.deSerialize( base_node, materialDb );
  }catch( std::exception &e )
  {
    throw runtime_error( "Failed to load preset configuration: " + string(e.what()) );
  }

  // Save XML document to SpecMeas
  shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw runtime_error( "No foreground spectrum loaded" );

  spec->setRelActAutoGuiState( &state );

  RelActAutoGui *gui = interspec->relActAutoWindow(false);
  if( gui )
  {
    try
    {
      gui->deSerialize( state );
    }catch( std::exception &e )
    {
      cerr << "Failed to de-serialize preset '" << preset_name << "' to the GUI!  error: " << e.what() << endl;
    }
  }//if( gui )

  // Build result with configuration summary
  result["success"] = true;
  result["preset_name"] = SpecUtils::filename(preset_path);
  result["preset_path"] = preset_path;

  // Use summary helper to get configuration info
  //json summary = stateToJsonSummary( state );
  json summary = stateToJsonDetailed( state );
  result["configuration"] = summary;

  return result;
}//executeLoadIsotopicsPreset()


nlohmann::json executePerformIsotopics(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  using namespace std;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  // Get current configuration - REQUIRED (throws if not present)
  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

  json result;

  // Extract options from state
  const RelActCalcAuto::Options &options = state.options;

  // Get foreground spectrum
  shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw runtime_error( "No foreground spectrum loaded" );

  shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    throw runtime_error( "No foreground histogram available" );

  // Get background if requested
  shared_ptr<const SpecUtils::Measurement> background;
  if( state.background_subtract )
  {
    background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
    // background can be null - solve() will handle it
  }

  // Get detector response function
  shared_ptr<const DetectorPeakResponse> drf = spec->detector();

  // Get peaks
  PeakModel *peakModel = interspec->peakModel();
  vector<shared_ptr<const PeakDef>> all_peaks;

  if( peakModel )
  {
    shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = peakModel->peaks();
    if( peaks )
    {
      for( const auto &peak : *peaks )
      {
        if( peak )
          all_peaks.push_back( peak );
      }
    }
  }

  // Execute the solve
  try
  {
    RelActCalcAuto::RelActAutoSolution solution = RelActCalcAuto::solve(
      options,
      foreground,
      background,
      drf,
      all_peaks,
      nullptr  // no cancel callback
    );

    // Check if solve was successful
    if( solution.m_status != RelActCalcAuto::RelActAutoSolution::Status::Success )
    {
      result["success"] = false;
      result["status"] = "Failed";
      result["error"] = solution.m_error_message;

      if( !solution.m_warnings.empty() )
        result["warnings"] = solution.m_warnings;

      return result;
    }

    // Build success result
    result["success"] = true;
    result["status"] = "Success";

    // Quality metrics
    json quality;
    quality["chi2"] = solution.m_chi2;
    quality["dof"] = solution.m_dof;
    quality["chi2_per_dof"] = (solution.m_dof > 0) ? (solution.m_chi2 / solution.m_dof) : -1.0;
    quality["num_function_eval"] = solution.m_num_function_eval_solution;
    result["quality"] = quality;

    // Extract isotopics results (first rel eff curve)
    if( solution.m_rel_activities.empty() || solution.m_rel_activities[0].empty() )
      throw runtime_error( "No relative activities in solution" );

    const vector<RelActCalcAuto::NuclideRelAct> &rel_acts = solution.m_rel_activities[0];

    // Calculate mass fractions for all isotopes
    json isotopics = json::array();

    for( const RelActCalcAuto::NuclideRelAct &nuc_act : rel_acts )
    {
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( nuc_act.source );
      if( !nuc )
        continue;

      json nuc_result;
      nuc_result["nuclide"] = nuc->symbol;

      // Get mass fraction
      try
      {
        const pair<double,optional<double>> mass_frac = solution.mass_enrichment_fraction( nuc, 0 );
        nuc_result["mass_fraction"] = mass_frac.first;
        if( mass_frac.second )
          nuc_result["mass_fraction_uncert"] = *mass_frac.second;
      }catch( exception &e )
      {
        nuc_result["mass_fraction_error"] = e.what();
      }

      // Relative activity
      nuc_result["rel_activity"] = nuc_act.rel_activity;
      nuc_result["rel_activity_uncert"] = nuc_act.rel_activity_uncertainty;

      // Age
      if( nuc_act.age >= 0.0 )
      {
        nuc_result["age_days"] = nuc_act.age / PhysicalUnits::day;
        nuc_result["age_uncert_days"] = nuc_act.age_uncertainty / PhysicalUnits::day;
        nuc_result["age_was_fit"] = nuc_act.age_was_fit;
      }

      isotopics.push_back( nuc_result );
    }

    result["isotopics"] = isotopics;

    // Include Pu242 correction if available
    if( !solution.m_corrected_pu.empty() && solution.m_corrected_pu[0] )
    {
      json pu242_corr;
      const auto &corr_pu = solution.m_corrected_pu[0];

      pu242_corr["pu238"] = corr_pu->pu238_mass_frac;
      pu242_corr["pu239"] = corr_pu->pu239_mass_frac;
      pu242_corr["pu240"] = corr_pu->pu240_mass_frac;
      pu242_corr["pu241"] = corr_pu->pu241_mass_frac;
      pu242_corr["pu242"] = corr_pu->pu242_mass_frac;

      result["pu242_corrected"] = pu242_corr;
    }

    // Warnings
    if( !solution.m_warnings.empty() )
      result["warnings"] = solution.m_warnings;

    // Save updated state back to SpecMeas
    saveRelActState( interspec, state );

  }catch( exception &e )
  {
    result["success"] = false;
    result["error"] = e.what();
    return result;
  }

  return result;
}//executePerformIsotopics()


namespace
{
  // Helper to get the single relative efficiency curve from the state.
  // Throws if no curves exist or if multiple curves exist.
  RelActCalcAuto::RelEffCurveInput &getSingleCurve( RelActCalcAuto::RelActAutoGuiState &state )
  {
    if( state.options.rel_eff_curves.empty() )
      throw std::runtime_error( "No relative efficiency curves in configuration" );

    if( state.options.rel_eff_curves.size() > 1 )
      throw std::runtime_error( "Editing multi-curve configurations is not supported. "
                               "Current configuration has " + std::to_string( state.options.rel_eff_curves.size() ) + " curves." );

    return state.options.rel_eff_curves[0];
  }//getSingleCurve()


  // Helper to parse a nuclide/element/reaction from JSON and create a NucInputInfo entry.
  // Handles the "name" field which can be a nuclide symbol (e.g., "Pu239"),
  // an element symbol for x-rays (e.g., "Pu"), or a reaction name (e.g., "Fe(n,n')").
  RelActCalcAuto::NucInputInfo parseNuclideFromJson( const nlohmann::json &nuc_json )
  {
    using namespace std;

    const string source_name = nuc_json.at( "name" ).get<string>();

    // Try to parse as nuclide, element, or reaction
    RelActCalcAuto::SrcVariant source = RelActCalcAuto::source_from_string( source_name );

    if( RelActCalcAuto::is_null( source ) )
      throw runtime_error( "Unknown source '" + source_name + "'. Must be a valid nuclide (e.g., 'Pu239'), "
                          "element for x-rays (e.g., 'Pu'), or reaction (e.g., 'Fe(n,n')')." );

    RelActCalcAuto::NucInputInfo nuc_info;
    nuc_info.source = source;

    // Get the nuclide pointer to check if this is a nuclide vs element/reaction
    const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( source );
    const SandiaDecay::Element *el = RelActCalcAuto::element( source );
    const ReactionGamma::Reaction *rxn = RelActCalcAuto::reaction( source );

    // Handle age - only applicable for nuclides
    if( nuc )
    {
      if( nuc_json.contains( "age" ) )
      {
        const string age_str = nuc_json.at( "age" ).get<string>();
        nuc_info.age = PhysicalUnitsLocalized::stringToTimeDuration( age_str );
        if( nuc_info.age < 0 )
          throw runtime_error( "Age must be >= 0 for nuclide '" + source_name + "'" );
      }
      else
      {
        // Use default age for the nuclide
        nuc_info.age = PeakDef::defaultDecayTime( nuc );
      }

      // Handle fit_age
      if( nuc_json.contains( "fit_age" ) )
        nuc_info.fit_age = nuc_json.at( "fit_age" ).get<bool>();

      // Handle age bounds - only if fit_age is true
      if( nuc_info.fit_age )
      {
        if( nuc_json.contains( "age_min" ) )
        {
          const string age_min_str = nuc_json.at( "age_min" ).get<string>();
          nuc_info.fit_age_min = PhysicalUnitsLocalized::stringToTimeDuration( age_min_str );
          if( *nuc_info.fit_age_min < 0 )
            throw runtime_error( "age_min must be >= 0" );
        }

        if( nuc_json.contains( "age_max" ) )
        {
          const string age_max_str = nuc_json.at( "age_max" ).get<string>();
          nuc_info.fit_age_max = PhysicalUnitsLocalized::stringToTimeDuration( age_max_str );
          if( nuc_info.fit_age_min && (*nuc_info.fit_age_max <= *nuc_info.fit_age_min) )
            throw runtime_error( "age_max must be > age_min" );
        }
      }//if( fit_age )
    }
    else if( el || rxn )
    {
      // For elements (x-rays) and reactions, age must be 0 or negative (indicating N/A)
      nuc_info.age = -1.0;
      nuc_info.fit_age = false;
    }

    // Handle activity bounds (applicable to all source types)
    if( nuc_json.contains( "min_rel_act" ) )
      nuc_info.min_rel_act = nuc_json.at( "min_rel_act" ).get<double>();

    if( nuc_json.contains( "max_rel_act" ) )
    {
      nuc_info.max_rel_act = nuc_json.at( "max_rel_act" ).get<double>();
      if( nuc_info.min_rel_act && (*nuc_info.max_rel_act < *nuc_info.min_rel_act) )
        throw runtime_error( "max_rel_act must be >= min_rel_act" );
    }

    return nuc_info;
  }//parseNuclideFromJson()


  // Helper to parse an ROI from JSON
  RelActCalcAuto::RoiRange parseRoiFromJson( const nlohmann::json &roi_json )
  {
    RelActCalcAuto::RoiRange roi;

    roi.lower_energy = roi_json.at( "lower_energy" ).get<double>();
    roi.upper_energy = roi_json.at( "upper_energy" ).get<double>();

    if( roi.lower_energy >= roi.upper_energy )
      throw std::runtime_error( "lower_energy (" + std::to_string( roi.lower_energy )
                                + ") must be less than upper_energy (" + std::to_string( roi.upper_energy ) + ")" );

    if( roi_json.contains( "continuum_type" ) )
    {
      const std::string continuum_str = roi_json.at( "continuum_type" ).get<std::string>();
      roi.continuum_type = PeakContinuum::str_to_offset_type_str( continuum_str.c_str(), continuum_str.size() );
    }

    if( roi_json.contains( "range_type" ) )
    {
      const std::string range_str = roi_json.at( "range_type" ).get<std::string>();
      roi.range_limits_type = RelActCalcAuto::RoiRange::range_limits_type_from_str( range_str.c_str() );
    }else
    {
      roi.range_limits_type = RelActCalcAuto::RoiRange::RangeLimitsType::Fixed;
    }

    return roi;
  }//parseRoiFromJson()


  // Helper to validate that ROIs do not overlap
  void validateRoisNonOverlapping( const std::vector<RelActCalcAuto::RoiRange> &rois )
  {
    if( rois.size() < 2 )
      return;

    // Create a sorted copy by lower_energy
    std::vector<const RelActCalcAuto::RoiRange *> sorted_rois;
    for( const auto &roi : rois )
      sorted_rois.push_back( &roi );

    std::sort( sorted_rois.begin(), sorted_rois.end(),
      []( const RelActCalcAuto::RoiRange *a, const RelActCalcAuto::RoiRange *b ) {
        return a->lower_energy < b->lower_energy;
      });

    for( size_t i = 1; i < sorted_rois.size(); ++i )
    {
      if( sorted_rois[i]->lower_energy < sorted_rois[i-1]->upper_energy )
      {
        throw std::runtime_error( "ROI [" + std::to_string( sorted_rois[i-1]->lower_energy )
                                 + ", " + std::to_string( sorted_rois[i-1]->upper_energy )
                                 + "] overlaps with ROI [" + std::to_string( sorted_rois[i]->lower_energy )
                                 + ", " + std::to_string( sorted_rois[i]->upper_energy ) + "]" );
      }
    }
  }//validateRoisNonOverlapping()


  // Helper to parse shielding from JSON
  std::shared_ptr<RelActCalc::PhysicalModelShieldInput> parseShieldingFromJson(
    const nlohmann::json &shield_json,
    MaterialDB *materialDb )
  {
    auto shield = std::make_shared<RelActCalc::PhysicalModelShieldInput>();

    if( shield_json.contains( "material" ) )
    {
      const std::string mat_name = shield_json.at( "material" ).get<std::string>();
      const Material *mat = materialDb ? materialDb->material( mat_name ) : nullptr;
      if( !mat )
        throw std::runtime_error( "Unknown material '" + mat_name + "'" );
      shield->material = std::shared_ptr<const Material>( mat, [](const Material*){} ); // Non-owning
      shield->atomic_number = 0.0; // Must be 0 when material is specified

      // Handle thickness (only valid with material)
      if( shield_json.contains( "thickness" ) )
      {
        const std::string thickness_str = shield_json.at( "thickness" ).get<std::string>();
        const double thickness = PhysicalUnits::stringToDistance( thickness_str );
        shield->areal_density = thickness * mat->density;
      }
      else if( shield_json.contains( "areal_density" ) )
      {
        shield->areal_density = shield_json.at( "areal_density" ).get<double>() * PhysicalUnits::g / PhysicalUnits::cm2;
      }
    }
    else if( shield_json.contains( "atomic_number" ) )
    {
      shield->atomic_number = shield_json.at( "atomic_number" ).get<double>();
      if( shield->atomic_number < 1.0 || shield->atomic_number > 98.0 )
        throw std::runtime_error( "atomic_number must be in range [1, 98]" );

      if( shield_json.contains( "areal_density" ) )
        shield->areal_density = shield_json.at( "areal_density" ).get<double>() * PhysicalUnits::g / PhysicalUnits::cm2;
    }
    else
    {
      throw std::runtime_error( "Shielding must specify either 'material' or 'atomic_number'" );
    }

    if( shield_json.contains( "fit_areal_density" ) )
      shield->fit_areal_density = shield_json.at( "fit_areal_density" ).get<bool>();

    if( shield_json.contains( "fit_atomic_number" ) )
    {
      shield->fit_atomic_number = shield_json.at( "fit_atomic_number" ).get<bool>();
      if( shield->fit_atomic_number && shield->material )
        throw std::runtime_error( "fit_atomic_number cannot be true when a material is specified" );
    }

    return shield;
  }//parseShieldingFromJson()

}//namespace


nlohmann::json executeModifyIsotopicsNuclides(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  // Get or create current state
  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, true );

  // If we just created the state, add a default rel eff curve
  if( state.options.rel_eff_curves.empty() )
  {
    RelActCalcAuto::RelEffCurveInput default_curve;
    default_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
    default_curve.rel_eff_eqn_order = 3;
    state.options.rel_eff_curves.push_back( default_curve );
  }

  // Get the single curve (throws if multiple)
  RelActCalcAuto::RelEffCurveInput &curve = getSingleCurve( state );

  const string action = params.at( "action" ).get<string>();

  if( !params.contains( "nuclides" ) || !params["nuclides"].is_array() )
    throw runtime_error( "nuclides array parameter is required" );

  const json &nuclides_array = params["nuclides"];

  json result;
  result["modified"] = json::array();

  if( action == "add" )
  {
    for( const auto &nuc_json : nuclides_array )
    {
      RelActCalcAuto::NucInputInfo nuc_info = parseNuclideFromJson( nuc_json );

      // Check if nuclide already exists
      const string nuc_name = nuc_info.name();
      bool found = false;
      for( const auto &existing : curve.nuclides )
      {
        if( existing.name() == nuc_name )
        {
          found = true;
          break;
        }
      }

      if( found )
        throw runtime_error( "Source '" + nuc_name + "' already exists in configuration. Use 'update' action to modify." );

      curve.nuclides.push_back( nuc_info );
      result["modified"].push_back( nuc_name );
    }
  }
  else if( action == "remove" )
  {
    for( const auto &nuc_json : nuclides_array )
    {
      const string nuc_name = nuc_json.at( "name" ).get<string>();

      auto it = std::find_if( curve.nuclides.begin(), curve.nuclides.end(),
        [&nuc_name]( const RelActCalcAuto::NucInputInfo &n ) {
          return n.name() == nuc_name;
        });

      if( it == curve.nuclides.end() )
        throw runtime_error( "Source '" + nuc_name + "' not found in configuration" );

      curve.nuclides.erase( it );
      result["modified"].push_back( nuc_name );
    }
  }
  else if( action == "update" )
  {
    for( const auto &nuc_json : nuclides_array )
    {
      const string nuc_name = nuc_json.at( "name" ).get<string>();

      auto it = std::find_if( curve.nuclides.begin(), curve.nuclides.end(),
        [&nuc_name]( const RelActCalcAuto::NucInputInfo &n ) {
          return n.name() == nuc_name;
        });

      if( it == curve.nuclides.end() )
        throw runtime_error( "Source '" + nuc_name + "' not found in configuration. Use 'add' action to add new sources." );

      // Update only the fields that are specified
      const SandiaDecay::Nuclide *nuc = RelActCalcAuto::nuclide( it->source );

      if( nuc_json.contains( "age" ) )
      {
        if( !nuc )
          throw runtime_error( "Cannot set age for x-ray or reaction source '" + nuc_name + "'" );

        const string age_str = nuc_json.at( "age" ).get<string>();
        it->age = PhysicalUnitsLocalized::stringToTimeDuration( age_str );
        if( it->age < 0 )
          throw runtime_error( "Age must be >= 0" );
      }

      if( nuc_json.contains( "fit_age" ) )
      {
        if( !nuc )
          throw runtime_error( "Cannot set fit_age for x-ray or reaction source '" + nuc_name + "'" );
        it->fit_age = nuc_json.at( "fit_age" ).get<bool>();
      }

      if( nuc_json.contains( "age_min" ) )
      {
        const string age_min_str = nuc_json.at( "age_min" ).get<string>();
        it->fit_age_min = PhysicalUnitsLocalized::stringToTimeDuration( age_min_str );
      }

      if( nuc_json.contains( "age_max" ) )
      {
        const string age_max_str = nuc_json.at( "age_max" ).get<string>();
        it->fit_age_max = PhysicalUnitsLocalized::stringToTimeDuration( age_max_str );
      }

      if( nuc_json.contains( "min_rel_act" ) )
        it->min_rel_act = nuc_json.at( "min_rel_act" ).get<double>();

      if( nuc_json.contains( "max_rel_act" ) )
        it->max_rel_act = nuc_json.at( "max_rel_act" ).get<double>();

      result["modified"].push_back( nuc_name );
    }
  }
  else
  {
    throw runtime_error( "Unknown action '" + action + "'. Valid actions: add, remove, update" );
  }

  // Validate constraints after modification
  curve.check_nuclide_constraints();

  // Save state back
  saveRelActState( interspec, state );

  result["success"] = true;
  result["message"] = "Nuclides " + action + " successful";
  result["nuclide_count"] = curve.nuclides.size();

  return result;
}//executeModifyIsotopicsNuclides()


nlohmann::json executeModifyIsotopicsRois(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, true );

  // If we just created the state, add a default rel eff curve
  if( state.options.rel_eff_curves.empty() )
  {
    RelActCalcAuto::RelEffCurveInput default_curve;
    default_curve.rel_eff_eqn_type = RelActCalc::RelEffEqnForm::LnX;
    default_curve.rel_eff_eqn_order = 3;
    state.options.rel_eff_curves.push_back( default_curve );
  }

  const string action = params.at( "action" ).get<string>();

  json result;
  result["modified"] = json::array();

  if( action == "add" )
  {
    if( !params.contains( "rois" ) || !params["rois"].is_array() )
      throw runtime_error( "rois array parameter is required for add action" );

    for( const auto &roi_json : params["rois"] )
    {
      RelActCalcAuto::RoiRange roi = parseRoiFromJson( roi_json );
      state.options.rois.push_back( roi );

      json roi_result;
      roi_result["lower_energy"] = roi.lower_energy;
      roi_result["upper_energy"] = roi.upper_energy;
      result["modified"].push_back( roi_result );
    }

    // Validate no overlaps
    validateRoisNonOverlapping( state.options.rois );
  }
  else if( action == "remove" )
  {
    if( !params.contains( "rois" ) || !params["rois"].is_array() )
      throw runtime_error( "rois array parameter is required for remove action" );

    for( const auto &roi_json : params["rois"] )
    {
      const double lower = roi_json.at( "lower_energy" ).get<double>();
      const double upper = roi_json.at( "upper_energy" ).get<double>();

      // Find ROI by energy bounds (with small tolerance)
      auto it = std::find_if( state.options.rois.begin(), state.options.rois.end(),
        [lower, upper]( const RelActCalcAuto::RoiRange &roi ) {
          return (fabs( roi.lower_energy - lower ) < 0.1) && (fabs( roi.upper_energy - upper ) < 0.1);
        });

      if( it == state.options.rois.end() )
        throw runtime_error( "ROI [" + to_string( lower ) + ", " + to_string( upper ) + "] not found" );

      json roi_result;
      roi_result["lower_energy"] = it->lower_energy;
      roi_result["upper_energy"] = it->upper_energy;
      result["modified"].push_back( roi_result );

      state.options.rois.erase( it );
    }
  }
  else if( action == "update" )
  {
    if( !params.contains( "rois" ) || !params["rois"].is_array() )
      throw runtime_error( "rois array parameter is required for update action" );

    for( const auto &roi_json : params["rois"] )
    {
      const double lower = roi_json.at( "lower_energy" ).get<double>();
      const double upper = roi_json.at( "upper_energy" ).get<double>();

      auto it = std::find_if( state.options.rois.begin(), state.options.rois.end(),
        [lower, upper]( const RelActCalcAuto::RoiRange &roi ) {
          return (fabs( roi.lower_energy - lower ) < 0.1) && (fabs( roi.upper_energy - upper ) < 0.1);
        });

      if( it == state.options.rois.end() )
        throw runtime_error( "ROI [" + to_string( lower ) + ", " + to_string( upper ) + "] not found" );

      // Update bounds if new values specified
      if( roi_json.contains( "new_lower_energy" ) )
        it->lower_energy = roi_json.at( "new_lower_energy" ).get<double>();

      if( roi_json.contains( "new_upper_energy" ) )
        it->upper_energy = roi_json.at( "new_upper_energy" ).get<double>();

      if( it->lower_energy >= it->upper_energy )
        throw runtime_error( "lower_energy must be less than upper_energy" );

      if( roi_json.contains( "continuum_type" ) )
      {
        const string continuum_str = roi_json.at( "continuum_type" ).get<string>();
        it->continuum_type = PeakContinuum::str_to_offset_type_str( continuum_str.c_str(), continuum_str.size() );
      }

      if( roi_json.contains( "range_type" ) )
      {
        const string range_str = roi_json.at( "range_type" ).get<string>();
        it->range_limits_type = RelActCalcAuto::RoiRange::range_limits_type_from_str( range_str.c_str() );
      }

      json roi_result;
      roi_result["lower_energy"] = it->lower_energy;
      roi_result["upper_energy"] = it->upper_energy;
      result["modified"].push_back( roi_result );
    }

    // Validate no overlaps after updates
    validateRoisNonOverlapping( state.options.rois );
  }
  else if( action == "replace_all" )
  {
    if( !params.contains( "rois" ) || !params["rois"].is_array() )
      throw runtime_error( "rois array parameter is required for replace_all action" );

    state.options.rois.clear();

    for( const auto &roi_json : params["rois"] )
    {
      RelActCalcAuto::RoiRange roi = parseRoiFromJson( roi_json );
      state.options.rois.push_back( roi );

      json roi_result;
      roi_result["lower_energy"] = roi.lower_energy;
      roi_result["upper_energy"] = roi.upper_energy;
      result["modified"].push_back( roi_result );
    }

    // Validate no overlaps
    validateRoisNonOverlapping( state.options.rois );
  }
  else
  {
    throw runtime_error( "Unknown action '" + action + "'. Valid actions: add, remove, update, replace_all" );
  }

  // Save state back
  saveRelActState( interspec, state );

  result["success"] = true;
  result["message"] = "ROIs " + action + " successful";
  result["roi_count"] = state.options.rois.size();

  return result;
}//executeModifyIsotopicsRois()


nlohmann::json executeModifyIsotopicsCurveSettings(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  // Get current configuration - REQUIRED (throws if not present)
  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

  // Get the single curve (throws if multiple)
  RelActCalcAuto::RelEffCurveInput &curve = getSingleCurve( state );

  MaterialDB *materialDb = interspec->materialDataBase();

  json result;
  json changes = json::array();

  // Handle rel_eff_eqn_type
  if( params.contains( "rel_eff_eqn_type" ) )
  {
    const string eqn_type_str = params.at( "rel_eff_eqn_type" ).get<string>();
    curve.rel_eff_eqn_type = RelActCalc::rel_eff_eqn_form_from_str( eqn_type_str.c_str() );
    changes.push_back( "rel_eff_eqn_type=" + eqn_type_str );
  }

  // Handle rel_eff_eqn_order (only for non-physical models)
  if( params.contains( "rel_eff_eqn_order" ) )
  {
    if( curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      throw runtime_error( "rel_eff_eqn_order is not applicable for FramPhysicalModel. Use shielding settings instead." );

    curve.rel_eff_eqn_order = params.at( "rel_eff_eqn_order" ).get<size_t>();
    changes.push_back( "rel_eff_eqn_order=" + to_string( curve.rel_eff_eqn_order ) );
  }

  // Handle nucs_of_el_same_age
  if( params.contains( "nucs_of_el_same_age" ) )
  {
    curve.nucs_of_el_same_age = params.at( "nucs_of_el_same_age" ).get<bool>();
    changes.push_back( string( "nucs_of_el_same_age=" ) + (curve.nucs_of_el_same_age ? "true" : "false") );
  }

  // Physical model specific settings
  const bool is_physical = (curve.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel);

  if( params.contains( "use_hoerl" ) )
  {
    if( !is_physical )
      throw runtime_error( "use_hoerl is only applicable for FramPhysicalModel equation type" );

    curve.phys_model_use_hoerl = params.at( "use_hoerl" ).get<bool>();
    changes.push_back( string( "use_hoerl=" ) + (curve.phys_model_use_hoerl ? "true" : "false") );
  }

  // Handle self_attenuation
  if( params.contains( "self_attenuation" ) )
  {
    if( !is_physical )
      throw runtime_error( "self_attenuation is only applicable for FramPhysicalModel equation type" );

    const json &self_atten_json = params["self_attenuation"];
    if( self_atten_json.is_null() )
    {
      curve.phys_model_self_atten = nullptr;
      changes.push_back( "self_attenuation=null" );
    }
    else
    {
      curve.phys_model_self_atten = parseShieldingFromJson( self_atten_json, materialDb );
      changes.push_back( "self_attenuation updated" );
    }
  }

  // Handle external_attenuation
  if( params.contains( "external_attenuation" ) )
  {
    if( !is_physical )
      throw runtime_error( "external_attenuation is only applicable for FramPhysicalModel equation type" );

    const json &ext_atten_json = params["external_attenuation"];
    curve.phys_model_external_atten.clear();

    if( ext_atten_json.is_array() )
    {
      for( const auto &ext_json : ext_atten_json )
      {
        auto ext = parseShieldingFromJson( ext_json, materialDb );
        curve.phys_model_external_atten.push_back( ext );
      }
      changes.push_back( "external_attenuation=" + to_string( curve.phys_model_external_atten.size() ) + " layers" );
    }
    else if( ext_atten_json.is_null() )
    {
      changes.push_back( "external_attenuation=null" );
    }
  }

  // Validate state after modifications
  state.options.check_same_hoerl_and_external_shielding_specifications();

  // Save state back
  saveRelActState( interspec, state );

  result["success"] = true;
  result["changes"] = changes;
  result["rel_eff_eqn_type"] = RelActCalc::to_str( curve.rel_eff_eqn_type );

  return result;
}//executeModifyIsotopicsCurveSettings()


nlohmann::json executeModifyIsotopicsOptions(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  // Get current configuration - REQUIRED (throws if not present)
  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

  json result;
  json changes = json::array();

  // Handle fit_energy_cal
  if( params.contains( "fit_energy_cal" ) )
  {
    state.options.fit_energy_cal = params.at( "fit_energy_cal" ).get<bool>();
    changes.push_back( string( "fit_energy_cal=" ) + (state.options.fit_energy_cal ? "true" : "false") );
  }

  // Handle fwhm_form - map Polynomial_X to Berstein_X internally
  if( params.contains( "fwhm_form" ) )
  {
    string fwhm_str = params.at( "fwhm_form" ).get<string>();

    // Map Polynomial_X to Berstein_X (the internal form uses Berstein)
    if( fwhm_str.find( "Polynomial_" ) == 0 )
    {
      fwhm_str = "Berstein_" + fwhm_str.substr( 11 ); // Replace "Polynomial_" with "Berstein_"
    }

    state.options.fwhm_form = RelActCalcAuto::fwhm_form_from_str( fwhm_str.c_str() );
    changes.push_back( "fwhm_form=" + fwhm_str );
  }

  // Handle fwhm_estimation_method
  if( params.contains( "fwhm_estimation_method" ) )
  {
    const string method_str = params.at( "fwhm_estimation_method" ).get<string>();
    state.options.fwhm_estimation_method = RelActCalcAuto::fwhm_estimation_method_from_str( method_str.c_str() );
    changes.push_back( "fwhm_estimation_method=" + method_str );
  }

  // Handle skew_type
  if( params.contains( "skew_type" ) )
  {
    const string skew_str = params.at( "skew_type" ).get<string>();
    state.options.skew_type = PeakDef::skew_from_string( skew_str.c_str() );
    changes.push_back( "skew_type=" + skew_str );
  }

  // Handle background_subtract
  if( params.contains( "background_subtract" ) )
  {
    const bool bg_sub = params.at( "background_subtract" ).get<bool>();

    // Validate that a background spectrum is loaded if enabling
    if( bg_sub )
    {
      shared_ptr<const SpecMeas> bgMeas = interspec->measurment( SpecUtils::SpectrumType::Background );
      shared_ptr<const SpecUtils::Measurement> bgSpec = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );
      if( !bgMeas || !bgSpec || bgSpec->num_gamma_channels() < 16 )
        throw runtime_error( "background_subtract=true requires a background spectrum to be loaded" );
    }

    state.background_subtract = bg_sub;
    changes.push_back( string( "background_subtract=" ) + (bg_sub ? "true" : "false") );
  }

  // Handle note
  if( params.contains( "note" ) )
  {
    state.note = params.at( "note" ).get<string>();
    changes.push_back( "note updated" );
  }

  // Handle description
  if( params.contains( "description" ) )
  {
    state.description = params.at( "description" ).get<string>();
    changes.push_back( "description updated" );
  }

  // Handle show_ref_lines
  if( params.contains( "show_ref_lines" ) )
  {
    state.show_ref_lines = params.at( "show_ref_lines" ).get<bool>();
    changes.push_back( string( "show_ref_lines=" ) + (state.show_ref_lines ? "true" : "false") );
  }

  // Handle additional_br_uncert
  if( params.contains( "additional_br_uncert" ) )
  {
    state.options.additional_br_uncert = params.at( "additional_br_uncert" ).get<double>();
    changes.push_back( "additional_br_uncert=" + to_string( state.options.additional_br_uncert ) );
  }

  // Save state back
  saveRelActState( interspec, state );

  result["success"] = true;
  result["changes"] = changes;

  return result;
}//executeModifyIsotopicsOptions()


nlohmann::json executeModifyIsotopicsConstraints(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  if( !interspec )
    throw runtime_error( "InterSpec instance required" );

  // Get current configuration - REQUIRED (throws if not present)
  RelActCalcAuto::RelActAutoGuiState state = getOrCreateRelActState( interspec, false );

  // Get the single curve (throws if multiple)
  RelActCalcAuto::RelEffCurveInput &curve = getSingleCurve( state );

  const string action = params.at( "action" ).get<string>();
  const string constraint_type = params.at( "constraint_type" ).get<string>();

  json result;
  result["modified"] = json::array();

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Decay database not available" );

  if( constraint_type == "activity_ratio" )
  {
    if( action == "add" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string constrained_name = con_json.at( "nuclide" ).get<string>();
        const string controlling_name = con_json.at( "relative_to" ).get<string>();
        const double ratio = con_json.at( "ratio" ).get<double>();

        RelActCalcAuto::SrcVariant constrained_src = RelActCalcAuto::source_from_string( constrained_name );
        RelActCalcAuto::SrcVariant controlling_src = RelActCalcAuto::source_from_string( controlling_name );

        if( RelActCalcAuto::is_null( constrained_src ) )
          throw runtime_error( "Unknown source: " + constrained_name );
        if( RelActCalcAuto::is_null( controlling_src ) )
          throw runtime_error( "Unknown source: " + controlling_name );

        // Check that nuclide doesn't already have a constraint
        for( const auto &existing : curve.act_ratio_constraints )
        {
          if( RelActCalcAuto::to_name( existing.constrained_source ) == constrained_name )
            throw runtime_error( "Source '" + constrained_name + "' already has an activity ratio constraint. Use 'update' to modify." );
        }

        RelActCalcAuto::RelEffCurveInput::ActRatioConstraint constraint;
        constraint.constrained_source = constrained_src;
        constraint.controlling_source = controlling_src;
        constraint.constrained_to_controlled_activity_ratio = ratio;

        curve.act_ratio_constraints.push_back( constraint );
        result["modified"].push_back( constrained_name );
      }
    }
    else if( action == "remove" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string constrained_name = con_json.at( "nuclide" ).get<string>();

        auto it = std::find_if( curve.act_ratio_constraints.begin(), curve.act_ratio_constraints.end(),
          [&constrained_name]( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &c ) {
            return RelActCalcAuto::to_name( c.constrained_source ) == constrained_name;
          });

        if( it == curve.act_ratio_constraints.end() )
          throw runtime_error( "No activity ratio constraint found for '" + constrained_name + "'" );

        curve.act_ratio_constraints.erase( it );
        result["modified"].push_back( constrained_name );
      }
    }
    else if( action == "update" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string constrained_name = con_json.at( "nuclide" ).get<string>();

        auto it = std::find_if( curve.act_ratio_constraints.begin(), curve.act_ratio_constraints.end(),
          [&constrained_name]( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &c ) {
            return RelActCalcAuto::to_name( c.constrained_source ) == constrained_name;
          });

        if( it == curve.act_ratio_constraints.end() )
          throw runtime_error( "No activity ratio constraint found for '" + constrained_name + "'" );

        if( con_json.contains( "relative_to" ) )
        {
          const string controlling_name = con_json.at( "relative_to" ).get<string>();
          RelActCalcAuto::SrcVariant controlling_src = RelActCalcAuto::source_from_string( controlling_name );
          if( RelActCalcAuto::is_null( controlling_src ) )
            throw runtime_error( "Unknown source: " + controlling_name );
          it->controlling_source = controlling_src;
        }

        if( con_json.contains( "ratio" ) )
          it->constrained_to_controlled_activity_ratio = con_json.at( "ratio" ).get<double>();

        result["modified"].push_back( constrained_name );
      }
    }
    else if( action == "clear_all" )
    {
      const size_t count = curve.act_ratio_constraints.size();
      curve.act_ratio_constraints.clear();
      result["cleared"] = count;
    }
    else
    {
      throw runtime_error( "Unknown action '" + action + "'. Valid actions: add, remove, update, clear_all" );
    }
  }
  else if( constraint_type == "mass_fraction" )
  {
    if( action == "add" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string nuc_name = con_json.at( "nuclide" ).get<string>();
        const double lower = con_json.value( "lower_fraction", 0.0 );
        const double upper = con_json.at( "upper_fraction" ).get<double>();

        const SandiaDecay::Nuclide *nuc = db->nuclide( nuc_name );
        if( !nuc )
          throw runtime_error( "Mass fraction constraints only apply to nuclides, not '" + nuc_name + "'" );

        // Check for existing constraint
        for( const auto &existing : curve.mass_fraction_constraints )
        {
          if( existing.nuclide && existing.nuclide->symbol == nuc_name )
            throw runtime_error( "Nuclide '" + nuc_name + "' already has a mass fraction constraint. Use 'update' to modify." );
        }

        RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
        constraint.nuclide = nuc;
        constraint.lower_mass_fraction = lower;
        constraint.upper_mass_fraction = upper;

        curve.mass_fraction_constraints.push_back( constraint );
        result["modified"].push_back( nuc_name );
      }
    }
    else if( action == "remove" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string nuc_name = con_json.at( "nuclide" ).get<string>();

        auto it = std::find_if( curve.mass_fraction_constraints.begin(), curve.mass_fraction_constraints.end(),
          [&nuc_name]( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &c ) {
            return c.nuclide && c.nuclide->symbol == nuc_name;
          });

        if( it == curve.mass_fraction_constraints.end() )
          throw runtime_error( "No mass fraction constraint found for '" + nuc_name + "'" );

        curve.mass_fraction_constraints.erase( it );
        result["modified"].push_back( nuc_name );
      }
    }
    else if( action == "update" )
    {
      if( !params.contains( "constraints" ) || !params["constraints"].is_array() )
        throw runtime_error( "constraints array parameter is required" );

      for( const auto &con_json : params["constraints"] )
      {
        const string nuc_name = con_json.at( "nuclide" ).get<string>();

        auto it = std::find_if( curve.mass_fraction_constraints.begin(), curve.mass_fraction_constraints.end(),
          [&nuc_name]( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &c ) {
            return c.nuclide && c.nuclide->symbol == nuc_name;
          });

        if( it == curve.mass_fraction_constraints.end() )
          throw runtime_error( "No mass fraction constraint found for '" + nuc_name + "'" );

        if( con_json.contains( "lower_fraction" ) )
          it->lower_mass_fraction = con_json.at( "lower_fraction" ).get<double>();

        if( con_json.contains( "upper_fraction" ) )
          it->upper_mass_fraction = con_json.at( "upper_fraction" ).get<double>();

        result["modified"].push_back( nuc_name );
      }
    }
    else if( action == "clear_all" )
    {
      const size_t count = curve.mass_fraction_constraints.size();
      curve.mass_fraction_constraints.clear();
      result["cleared"] = count;
    }
    else
    {
      throw runtime_error( "Unknown action '" + action + "'. Valid actions: add, remove, update, clear_all" );
    }
  }
  else
  {
    throw runtime_error( "Unknown constraint_type '" + constraint_type + "'. Valid types: activity_ratio, mass_fraction" );
  }

  // Validate constraints after modification
  curve.check_nuclide_constraints();

  // Save state back
  saveRelActState( interspec, state );

  result["success"] = true;
  result["message"] = constraint_type + " constraints " + action + " successful";

  return result;
}//executeModifyIsotopicsConstraints()


nlohmann::json executeGetIsotopicsConfigSchema(
  const nlohmann::json& params,
  InterSpec* interspec )
{
  using namespace std;
  using json = nlohmann::json;

  json result;

  // ============================================================================
  // Configuration Structure Overview
  // ============================================================================
  result["structure"] = json::object();
  result["structure"]["description"] = "Isotopics configuration has these main sections: "
    "1) 'rel_eff_curve' - relative efficiency curve settings (equation type, shielding for physical model), "
    "2) 'nuclides' - list of nuclides/sources to analyze with ages and constraints, "
    "3) 'rois' - energy regions of interest to fit, "
    "4) Global options (fwhm_form, skew_type, fit_energy_cal, etc.). "
    "Use get_isotopics_config to see the current configuration structure.";

  // ============================================================================
  // Relative Efficiency Curve Settings
  // ============================================================================
  result["rel_eff_curve"] = json::object();
  result["rel_eff_curve"]["description"] = "Settings for the relative efficiency curve. "
    "For non-physical models (LnX, LnY, LnXLnY), use rel_eff_eqn_order to set polynomial order. "
    "For FramPhysicalModel, use self_attenuation/external_attenuation for shielding and use_hoerl option.";

  result["rel_eff_curve"]["fields"] = json::object();
  result["rel_eff_curve"]["fields"]["rel_eff_eqn_type"] = "The relative efficiency equation form. "
    "Non-physical models (LnX, LnY, LnXLnY) fit a polynomial; FramPhysicalModel uses physics-based attenuation.";
  result["rel_eff_curve"]["fields"]["rel_eff_eqn_order"] = "Polynomial order for non-physical models (typically 2-6). "
    "Higher orders allow more flexibility but may overfit. Ignored for FramPhysicalModel.";
  result["rel_eff_curve"]["fields"]["nucs_of_el_same_age"] = "If true, nuclides of the same element are constrained to have the same age.";
  result["rel_eff_curve"]["fields"]["use_hoerl"] = "(FramPhysicalModel only) Whether to use Hoerl form for the relative efficiency equation.";
  result["rel_eff_curve"]["fields"]["self_attenuation"] = "(FramPhysicalModel only) Self-attenuation of the source material. "
    "Specify material+thickness or atomic_number+areal_density.";
  result["rel_eff_curve"]["fields"]["external_attenuation"] = "(FramPhysicalModel only) Array of external shielding layers between source and detector.";

  // Rel eff equation types
  result["rel_eff_eqn_types"] = json::array();
  result["rel_eff_eqn_types"].push_back( {{"value", "LnX"}, {"description", "Log of energy polynomial - efficiency = exp(sum of a_i * ln(E)^i)"}} );
  result["rel_eff_eqn_types"].push_back( {{"value", "LnY"}, {"description", "Log of efficiency polynomial - ln(efficiency) = sum of a_i * E^i"}} );
  result["rel_eff_eqn_types"].push_back( {{"value", "LnXLnY"}, {"description", "Log-log polynomial - ln(efficiency) = sum of a_i * ln(E)^i. Common for HPGe."}} );
  result["rel_eff_eqn_types"].push_back( {{"value", "FramPhysicalModel"}, {"description", "Physical model using self-attenuation and external shielding. "
    "Requires shielding parameters. Use when source material attenuation is significant."}} );

  // ============================================================================
  // Nuclide/Source Fields
  // ============================================================================
  result["nuclides"] = json::object();
  result["nuclides"]["description"] = "Array of sources (nuclides, x-ray elements, or reactions) to include in the analysis. "
    "Each source contributes gamma lines that are fit within the ROIs.";

  result["nuclides"]["fields"] = json::object();
  result["nuclides"]["fields"]["name"] = "Source identifier. Can be: "
    "nuclide symbol (e.g., 'Pu239', 'U235', 'Am241'), "
    "element for x-rays (e.g., 'Pu', 'U'), or "
    "reaction (e.g., 'Fe(n,n')', 'H(n,n')').";
  result["nuclides"]["fields"]["age"] = "(Nuclides only) Time since chemical purification, with units (e.g., '5 years', '20 days', '1826.25 days'). "
    "Affects relative intensities of gamma lines from decay chain. Not applicable for x-rays or reactions.";
  result["nuclides"]["fields"]["age_days"] = "(Nuclides only) Alternative to 'age' - numeric age in days.";
  result["nuclides"]["fields"]["fit_age"] = "(Nuclides only) If true, age will be fitted during analysis within age_min/age_max bounds.";
  result["nuclides"]["fields"]["age_min"] = "(Nuclides only) Minimum age bound for fitting (with units). Required if fit_age=true.";
  result["nuclides"]["fields"]["age_max"] = "(Nuclides only) Maximum age bound for fitting (with units). Required if fit_age=true.";
  result["nuclides"]["fields"]["min_rel_act"] = "Minimum relative activity constraint for fitting.";
  result["nuclides"]["fields"]["max_rel_act"] = "Maximum relative activity constraint for fitting.";

  // ============================================================================
  // ROI (Region of Interest) Fields
  // ============================================================================
  result["rois"] = json::object();
  result["rois"]["description"] = "Array of energy regions to fit. Each ROI defines an energy range containing peaks to be analyzed. "
    "ROIs must not overlap. Peaks within ROIs are fit simultaneously.";

  result["rois"]["fields"] = json::object();
  result["rois"]["fields"]["lower_energy"] = "Lower energy bound of the ROI in keV.";
  result["rois"]["fields"]["upper_energy"] = "Upper energy bound of the ROI in keV. Must be > lower_energy.";
  result["rois"]["fields"]["continuum_type"] = "Type of continuum (background) model for this ROI.";
  result["rois"]["fields"]["range_type"] = "How the ROI bounds behave during fitting.";

  // Continuum types
  result["continuum_types"] = json::array();
  result["continuum_types"].push_back( {{"value", "Linear"}, {"description", "Linear continuum - 2 parameters"}} );
  result["continuum_types"].push_back( {{"value", "Quadratic"}, {"description", "Quadratic continuum - 3 parameters"}} );
  result["continuum_types"].push_back( {{"value", "Cubic"}, {"description", "Cubic continuum - 4 parameters"}} );
  result["continuum_types"].push_back( {{"value", "FlatStep"}, {"description", "Flat step continuum - accounts for Compton edge"}} );
  result["continuum_types"].push_back( {{"value", "LinearStep"}, {"description", "Linear step continuum"}} );
  result["continuum_types"].push_back( {{"value", "BiLinearStep"}, {"description", "Bi-linear step continuum"}} );

  // ROI range types
  result["range_types"] = json::array();
  result["range_types"].push_back( {{"value", "Fixed"}, {"description", "ROI bounds are fixed at specified energies"}} );
  result["range_types"].push_back( {{"value", "CanExpandForFwhm"}, {"description", "ROI can expand to accommodate peak FWHM"}} );
  result["range_types"].push_back( {{"value", "CanBeBrokenUp"}, {"description", "ROI can be split into multiple regions, based on gamma lines of source and FWHM of detector, if needed"}} );

  // ============================================================================
  // Global Options
  // ============================================================================
  result["global_options"] = json::object();
  result["global_options"]["description"] = "Settings that apply to the entire analysis.";

  result["global_options"]["fields"] = json::object();
  result["global_options"]["fields"]["fit_energy_cal"] = "If true, energy calibration is refined during fitting.";
  result["global_options"]["fields"]["fwhm_form"] = "Functional form for peak FWHM vs energy.";
  result["global_options"]["fields"]["fwhm_estimation_method"] = "How initial FWHM parameters are determined.";
  result["global_options"]["fields"]["skew_type"] = "Peak shape skew model (asymmetry).";
  result["global_options"]["fields"]["additional_br_uncert"] = "Additional fractional uncertainty added to branching ratios (e.g., 0.05 for 5%).";
  result["global_options"]["fields"]["background_subtract"] = "If true, subtract background spectrum before analysis. Requires background to be loaded.";
  result["global_options"]["fields"]["show_ref_lines"] = "If true, show reference lines in the display.";
  result["global_options"]["fields"]["note"] = "User notes for this configuration.";
  result["global_options"]["fields"]["description"] = "Description of this configuration.";

  // FWHM forms
  result["fwhm_forms"] = json::array();
  result["fwhm_forms"].push_back( {{"value", "Gadras"}, {"description", "GADRAS functional form: sqrt(a + b*E + c/E)"}} );
  result["fwhm_forms"].push_back( {{"value", "SqrtEnergyPlusInverse"}, {"description", "sqrt(a + b*E + c/E) - similar to Gadras"}} );
  result["fwhm_forms"].push_back( {{"value", "ConstantPlusSqrtEnergy"}, {"description", "a + b*sqrt(E)"}} );
  result["fwhm_forms"].push_back( {{"value", "Polynomial_2"}, {"description", "2nd order polynomial in energy"}} );
  result["fwhm_forms"].push_back( {{"value", "Polynomial_3"}, {"description", "3rd order polynomial in energy"}} );
  result["fwhm_forms"].push_back( {{"value", "Polynomial_4"}, {"description", "4th order polynomial in energy"}} );
  result["fwhm_forms"].push_back( {{"value", "Polynomial_5"}, {"description", "5th order polynomial in energy"}} );
  result["fwhm_forms"].push_back( {{"value", "Polynomial_6"}, {"description", "6th order polynomial in energy"}} );

  // FWHM estimation methods
  result["fwhm_estimation_methods"] = json::array();
  result["fwhm_estimation_methods"].push_back( {{"value", "StartFromDetEffOrPeaksInSpectrum"},
    {"description", "Use detector efficiency FWHM if available, otherwise estimate from peaks in spectrum, then refine during fit"}} );
  result["fwhm_estimation_methods"].push_back( {{"value", "StartFromPeaksInSpectrum"},
    {"description", "Estimate FWHM from peaks in spectrum, then refine during fit"}} );
  result["fwhm_estimation_methods"].push_back( {{"value", "StartFromDetectorEfficiency"},
    {"description", "Start from detector efficiency FWHM, then refine during fit"}} );
  result["fwhm_estimation_methods"].push_back( {{"value", "FixedToDetectorEfficiency"},
    {"description", "Use detector efficiency FWHM and keep fixed (do not fit)"}} );

  // Skew types
  result["skew_types"] = json::array();
  result["skew_types"].push_back( {{"value", "NoSkew"}, {"description", "Pure Gaussian peaks (symmetric)"}} );
  result["skew_types"].push_back( {{"value", "Bortel"}, {"description", "Bortel skew model - low-energy tail"}} );
  result["skew_types"].push_back( {{"value", "GaussExp"}, {"description", "Gaussian-Exponential hybrid - exponential tail on low side"}} );
  result["skew_types"].push_back( {{"value", "CrystalBall"}, {"description", "Crystal Ball function - power-law tail"}} );
  result["skew_types"].push_back( {{"value", "ExpGaussExp"}, {"description", "Exponential-Gaussian-Exponential - tails on both sides"}} );
  result["skew_types"].push_back( {{"value", "DoubleSidedCrystalBall"}, {"description", "Double-sided Crystal Ball - power-law tails on both sides"}} );

  // ============================================================================
  // Constraint Fields
  // ============================================================================
  result["constraints"] = json::object();
  result["constraints"]["description"] = "Activity ratio and mass fraction constraints between nuclides.";

  result["constraints"]["activity_ratio"] = json::object();
  result["constraints"]["activity_ratio"]["description"] = "Constrain activity ratio between two nuclides (e.g., Pu240/Pu239 = 0.06).";
  result["constraints"]["activity_ratio"]["fields"] = json::object();
  result["constraints"]["activity_ratio"]["fields"]["nuclide"] = "The constrained nuclide (e.g., 'Pu240').";
  result["constraints"]["activity_ratio"]["fields"]["relative_to"] = "The reference nuclide (e.g., 'Pu239').";
  result["constraints"]["activity_ratio"]["fields"]["ratio"] = "Activity ratio: constrained_activity / reference_activity.";

  result["constraints"]["mass_fraction"] = json::object();
  result["constraints"]["mass_fraction"]["description"] = "Constrain mass fraction of a nuclide within its element (e.g., U235 is 3-20% of total uranium).";
  result["constraints"]["mass_fraction"]["fields"] = json::object();
  result["constraints"]["mass_fraction"]["fields"]["nuclide"] = "The nuclide to constrain (e.g., 'U235').";
  result["constraints"]["mass_fraction"]["fields"]["lower_fraction"] = "Minimum mass fraction (0-1).";
  result["constraints"]["mass_fraction"]["fields"]["upper_fraction"] = "Maximum mass fraction (0-1).";

  // ============================================================================
  // Shielding Fields (for FramPhysicalModel)
  // ============================================================================
  result["shielding"] = json::object();
  result["shielding"]["description"] = "Shielding parameters for FramPhysicalModel. Can specify either material+thickness or atomic_number+areal_density.";

  result["shielding"]["fields"] = json::object();
  result["shielding"]["fields"]["material"] = "Material name from the materials database (e.g., 'Iron', 'Lead', 'plutonium_dioxide'). "
    "Use get_materials tool to list available materials.";
  result["shielding"]["fields"]["atomic_number"] = "Effective atomic number (1-98). Use instead of material for generic shielding. Cannot be used with 'material'.";
  result["shielding"]["fields"]["areal_density"] = "Areal density in g/cm2. This is mass per unit area = density * thickness.";
  result["shielding"]["fields"]["thickness"] = "Thickness with units (e.g., '1 cm', '2.54 mm'). Only valid when 'material' is specified. "
    "Internally converted to areal_density using material density.";
  result["shielding"]["fields"]["fit_areal_density"] = "If true, areal density will be fitted during analysis.";
  result["shielding"]["fields"]["fit_atomic_number"] = "If true, atomic number will be fitted. Only valid with 'atomic_number', not 'material'. "
    "Usually not recommended - prefer fitting areal_density instead.";

  result["success"] = true;

  return result;
}//executeGetIsotopicsConfigSchema()


}//namespace IsotopicsTool
}//namespace LlmTools

#endif //USE_LLM_INTERFACE
