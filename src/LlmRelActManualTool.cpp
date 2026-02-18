/* InterSpec: an application to analyze spectral gamma radiation data.

 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov.

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "InterSpec_config.h"
#include "InterSpec/LlmRelActManualTool.h"

#if( USE_LLM_INTERFACE )

#include <map>
#include <set>
#include <deque>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "rapidxml/rapidxml.hpp"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/AnalystChecks.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActManualGui.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using json = nlohmann::json;

namespace LlmTools
{
namespace RelActManualTool
{

namespace
{

/** Parse AddUncert enum from string.
 * @param str The string value (e.g., "StatOnly", "FiftyPercent")
 * @returns The corresponding AddUncert enum value
 * @throws runtime_error if invalid string
 */
RelActManualGui::AddUncert parseAddUncert( const string &str )
{
  if( str == "Unweighted" )
    return RelActManualGui::AddUncert::Unweighted;
  if( str == "StatOnly" )
    return RelActManualGui::AddUncert::StatOnly;
  if( str == "OnePercent" )
    return RelActManualGui::AddUncert::OnePercent;
  if( str == "FivePercent" )
    return RelActManualGui::AddUncert::FivePercent;
  if( str == "TenPercent" )
    return RelActManualGui::AddUncert::TenPercent;
  if( str == "TwentyFivePercent" )
    return RelActManualGui::AddUncert::TwentyFivePercent;
  if( str == "FiftyPercent" )
    return RelActManualGui::AddUncert::FiftyPercent;
  if( str == "SeventyFivePercent" )
    return RelActManualGui::AddUncert::SeventyFivePercent;
  if( str == "OneHundredPercent" )
    return RelActManualGui::AddUncert::OneHundredPercent;
  
  throw runtime_error( "Invalid additional_uncertainty value: '" + str + "'" );
}//parseAddUncert


/** Convert AddUncert enum to string.
 */
string addUncertToString( const RelActManualGui::AddUncert val )
{
  switch( val )
  {
    case RelActManualGui::AddUncert::Unweighted:         return "Unweighted";
    case RelActManualGui::AddUncert::StatOnly:           return "StatOnly";
    case RelActManualGui::AddUncert::OnePercent:         return "OnePercent";
    case RelActManualGui::AddUncert::FivePercent:        return "FivePercent";
    case RelActManualGui::AddUncert::TenPercent:         return "TenPercent";
    case RelActManualGui::AddUncert::TwentyFivePercent:  return "TwentyFivePercent";
    case RelActManualGui::AddUncert::FiftyPercent:       return "FiftyPercent";
    case RelActManualGui::AddUncert::SeventyFivePercent: return "SeventyFivePercent";
    case RelActManualGui::AddUncert::OneHundredPercent:  return "OneHundredPercent";
    case RelActManualGui::AddUncert::NumAddUncert:       break;
  }
  return "Unknown";
}//addUncertToString


/** Convert AddUncert enum to the numeric uncertainty value.
 * @returns The fractional uncertainty value, or -1.0 for Unweighted
 */
double addUncertToValue( const RelActManualGui::AddUncert val )
{
  switch( val )
  {
    case RelActManualGui::AddUncert::Unweighted:         return -1.0;
    case RelActManualGui::AddUncert::StatOnly:           return 0.0;
    case RelActManualGui::AddUncert::OnePercent:         return 0.01;
    case RelActManualGui::AddUncert::FivePercent:        return 0.05;
    case RelActManualGui::AddUncert::TenPercent:         return 0.1;
    case RelActManualGui::AddUncert::TwentyFivePercent:  return 0.25;
    case RelActManualGui::AddUncert::FiftyPercent:       return 0.5;
    case RelActManualGui::AddUncert::SeventyFivePercent: return 0.75;
    case RelActManualGui::AddUncert::OneHundredPercent:  return 1.0;
    case RelActManualGui::AddUncert::NumAddUncert:       break;
  }
  return 0.0;
}//addUncertToValue


/** Parse RelEffEqnForm from string.
 */
RelActCalc::RelEffEqnForm parseEqnForm( const string &str )
{
  if( str == "LnX" )
    return RelActCalc::RelEffEqnForm::LnX;
  if( str == "LnY" )
    return RelActCalc::RelEffEqnForm::LnY;
  if( str == "LnXLnY" )
    return RelActCalc::RelEffEqnForm::LnXLnY;
  if( str == "FramEmpirical" )
    return RelActCalc::RelEffEqnForm::FramEmpirical;
  if( str == "FramPhysicalModel" )
    return RelActCalc::RelEffEqnForm::FramPhysicalModel;
  
  throw runtime_error( "Invalid eqn_form value: '" + str + "'" );
}//parseEqnForm


/** Convert RelEffEqnForm to string.
 */
string eqnFormToString( const RelActCalc::RelEffEqnForm form )
{
  switch( form )
  {
    case RelActCalc::RelEffEqnForm::LnX:               return "LnX";
    case RelActCalc::RelEffEqnForm::LnY:               return "LnY";
    case RelActCalc::RelEffEqnForm::LnXLnY:            return "LnXLnY";
    case RelActCalc::RelEffEqnForm::FramEmpirical:     return "FramEmpirical";
    case RelActCalc::RelEffEqnForm::FramPhysicalModel: return "FramPhysicalModel";
  }
  return "Unknown";
}//eqnFormToString


/** Get default age for a nuclide based on its half-life.
 */
double getDefaultNuclideAge( const SandiaDecay::Nuclide *nuc )
{
  if( !nuc )
    return 0.0;
  
  // Use 5 half-lives as the default age for long-lived nuclides,
  // but cap at a reasonable value
  const double halfLife = nuc->halfLife;
  
  // For very long-lived nuclides (like U-238), use a reasonable age
  if( halfLife > 1.0E10 * PhysicalUnits::year )
    return 20.0 * PhysicalUnits::year;
  
  // For medium-lived nuclides, use 5 half-lives
  if( halfLife > PhysicalUnits::year )
    return std::min( 5.0 * halfLife, 100.0 * PhysicalUnits::year );
  
  // For short-lived nuclides, use 5 half-lives
  return 5.0 * halfLife;
}//getDefaultNuclideAge

}//anonymous namespace


nlohmann::json executePeakBasedRelativeEfficiency(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  using namespace RelActCalcManual;
  
  if( !interspec )
    throw runtime_error( "InterSpec instance required" );
  
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  if( !db )
    throw runtime_error( "Nuclear decay database not available" );
  
  // Get foreground spectrum
  shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( !spec )
    throw runtime_error( "No foreground spectrum loaded" );
  
  shared_ptr<const SpecUtils::Measurement> foreground = interspec->displayedHistogram( SpecUtils::SpectrumType::Foreground );
  if( !foreground )
    throw runtime_error( "No foreground histogram available" );
  
  // Parse input parameters
  const bool has_peak_energies = params.contains("peak_energies") && params["peak_energies"].is_array();
  const bool has_sources = params.contains("sources") && params["sources"].is_array();
  
  if( !has_peak_energies && !has_sources )
    throw runtime_error( "Either 'peak_energies' or 'sources' must be specified" );
  
  if( has_peak_energies && has_sources )
    throw runtime_error( "'peak_energies' and 'sources' are mutually exclusive" );

  const bool has_energy_ranges = params.contains( "energy_ranges" ) && params["energy_ranges"].is_array();
  if( has_energy_ranges && has_peak_energies )
    throw runtime_error( "'energy_ranges' can only be used with 'sources', not with 'peak_energies'" );

  const bool has_exclude_peaks = params.contains( "exclude_peak_energies" ) && params["exclude_peak_energies"].is_array();
  if( has_exclude_peaks && has_peak_energies )
    throw runtime_error( "'exclude_peak_energies' can only be used with 'sources', not with 'peak_energies'" );

  // Parse optional parameters
  RelActCalc::RelEffEqnForm eqn_form = RelActCalc::RelEffEqnForm::LnY;
  if( params.contains("eqn_form") && params["eqn_form"].is_string() )
    eqn_form = parseEqnForm( params["eqn_form"].get<string>() );
  
  size_t eqn_order = 3;
  if( params.contains("eqn_order") && params["eqn_order"].is_number_integer() )
    eqn_order = params["eqn_order"].get<size_t>();
  
  double match_tolerance = 1.5;
  if( params.contains("match_tolerance") && params["match_tolerance"].is_number() )
    match_tolerance = params["match_tolerance"].get<double>();
  
  bool background_subtract = false;
  if( params.contains("background_subtract") && params["background_subtract"].is_boolean() )
    background_subtract = params["background_subtract"].get<bool>();
  
  RelActManualGui::AddUncert add_uncert = RelActManualGui::AddUncert::StatOnly;
  if( params.contains("additional_uncertainty") && params["additional_uncertainty"].is_string() )
    add_uncert = parseAddUncert( params["additional_uncertainty"].get<string>() );
  
  bool save_to_state = false;
  if( params.contains("save_to_state") && params["save_to_state"].is_boolean() )
    save_to_state = params["save_to_state"].get<bool>();

  // Parse energy_ranges if provided (only valid with sources)
  vector<pair<double, double>> energy_ranges;
  if( has_energy_ranges )
  {
    const double spec_lower = static_cast<double>( foreground->gamma_energy_min() );
    const double spec_upper = static_cast<double>( foreground->gamma_energy_max() );

    for( const auto &range : params["energy_ranges"] )
    {
      if( !range.is_object() )
        continue;

      const bool has_lower = range.contains( "lower_energy" ) && range["lower_energy"].is_number();
      const bool has_upper = range.contains( "upper_energy" ) && range["upper_energy"].is_number();

      if( !has_lower && !has_upper )
        throw runtime_error( "Each energy range must specify at least one of 'lower_energy' or 'upper_energy'" );

      const double lower = has_lower ? range["lower_energy"].get<double>() : spec_lower;
      const double upper = has_upper ? range["upper_energy"].get<double>() : spec_upper;

      if( lower > upper )
        throw runtime_error( "Energy range lower bound (" + to_string( lower )
                            + " keV) must not exceed upper bound (" + to_string( upper ) + " keV)" );

      energy_ranges.emplace_back( lower, upper );
    }//for( each range in energy_ranges )

    if( energy_ranges.empty() )
      throw runtime_error( "'energy_ranges' was specified but contained no valid range objects" );
  }//if( has_energy_ranges )

  // Parse exclude_peak_energies if provided (only valid with sources)
  vector<double> exclude_peak_energies;
  if( has_exclude_peaks )
  {
    for( const auto &e : params["exclude_peak_energies"] )
    {
      if( e.is_number() )
        exclude_peak_energies.push_back( e.get<double>() );
    }
  }//if( has_exclude_peaks )

  // Parse nuclide_ages if provided
  map<string, double> nuclide_ages_input;
  if( params.contains("nuclide_ages") && params["nuclide_ages"].is_object() )
  {
    for( const auto &[nuc_name, age_val] : params["nuclide_ages"].items() )
    {
      if( age_val.is_string() )
      {
        const string age_str = age_val.get<string>();
        const SandiaDecay::Nuclide *nuc = db->nuclide( nuc_name );
        
        try
        {
          const double age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( age_str, nuc ? nuc->halfLife : -1.0 );
          nuclide_ages_input[nuc_name] = age;
        }
        catch( exception &e )
        {
          throw runtime_error( "Failed to parse age for " + nuc_name + ": " + e.what() );
        }
      }
      else if( age_val.is_number() )
      {
        // Assume seconds if numeric
        nuclide_ages_input[nuc_name] = age_val.get<double>();
      }
    }
  }//if( nuclide_ages provided )
  
  // Get peaks from PeakModel
  PeakModel *peakModel = interspec->peakModel();
  if( !peakModel )
    throw runtime_error( "Peak model not available" );
  
  shared_ptr<const deque<shared_ptr<const PeakDef>>> all_peaks = peakModel->peaks();
  if( !all_peaks || all_peaks->empty() )
    throw runtime_error( "No peaks are fit in the spectrum"
                        " - please add analysis peaks (i.e., using `add_analysis_peak` or `add_analysis_peaks_for_source`) "
                        "for the sources you want to analyze, before calling this tool." );

  // Get background peaks if needed
  deque<shared_ptr<const PeakDef>> background_peaks;
  shared_ptr<const SpecUtils::Measurement> background;
  if( background_subtract )
  {
    background = interspec->displayedHistogram( SpecUtils::SpectrumType::Background );

    AnalystChecks::DetectedPeaksOptions back_opts;
    back_opts.specType = SpecUtils::SpectrumType::Background;
    back_opts.nonBackgroundPeaksOnly = false;

    try
    {
      const AnalystChecks::DetectedPeakStatus back_status
        = AnalystChecks::detected_peaks( back_opts, interspec );

      for( const shared_ptr<const PeakDef> &p : back_status.peaks )
        background_peaks.push_back( p );
    }catch( ... )
    {
      // If background peak detection fails, just proceed without background peaks
    }
  }//if( background_subtract )
  
  // Collect peaks to use based on input specification
  vector<shared_ptr<const PeakDef>> selected_peaks;
  set<string> target_sources;
  
  if( has_peak_energies )
  {
    vector<double> target_energies;
    for( const auto &e : params["peak_energies"] )
    {
      if( e.is_number() )
        target_energies.push_back( e.get<double>() );
    }
    
    for( const double energy : target_energies )
    {
      // Find peak closest to this energy
      shared_ptr<const PeakDef> best_peak = nullptr;
      double best_dist = numeric_limits<double>::max();
      
      for( const shared_ptr<const PeakDef> &p : *all_peaks )
      {
        if( !p )
          continue;
        
        const double dist = fabs( p->mean() - energy );
        const double tolerance = p->gausPeak() ? 1.5*p->fwhm() : 0.25*(p->upperX() - p->lowerX());
        
        if( dist < tolerance && dist < best_dist )
        {
          best_peak = p;
          best_dist = dist;
        }
      }
      
      if( best_peak )
      {
        selected_peaks.push_back( best_peak );
        
        if( best_peak->parentNuclide() )
          target_sources.insert( best_peak->parentNuclide()->symbol );
        else if( best_peak->reaction() )
          target_sources.insert( best_peak->reaction()->name() );
      }
    }//for( target energies )
    
    if( selected_peaks.empty() )
      throw runtime_error( "No peaks found matching the specified energies" );
  }
  else // has_sources
  {
    for( const auto &s : params["sources"] )
    {
      if( s.is_string() )
        target_sources.insert( s.get<string>() );
    }
    
    for( const shared_ptr<const PeakDef> &p : *all_peaks )
    {
      if( !p )
        continue;
      
      string peak_source;
      if( p->parentNuclide() )
        peak_source = p->parentNuclide()->symbol;
      else if( p->reaction() )
        peak_source = p->reaction()->name();
      else
        continue;
      
      if( target_sources.count( peak_source ) )
      {
        if( !energy_ranges.empty() )
        {
          const double peak_energy = p->mean();

          bool in_range = false;
          for( const pair<double, double> &range : energy_ranges )
          {
            if( (peak_energy >= range.first) && (peak_energy <= range.second) )
            {
              in_range = true;
              break;
            }
          }//for( each energy range )

          if( !in_range )
            continue;
        }//if( energy_ranges specified )

        // Check if this peak should be excluded by energy
        if( !exclude_peak_energies.empty() )
        {
          const double peak_mean = p->mean();
          const double fwhm = p->gausPeak() ? p->fwhm() : 0.25 * (p->upperX() - p->lowerX());
          const double tol = 1.25 * fwhm;

          bool excluded = false;
          for( const double excl_energy : exclude_peak_energies )
          {
            if( fabs( peak_mean - excl_energy ) < tol )
            {
              excluded = true;
              break;
            }
          }//for( each exclude energy )

          if( excluded )
            continue;
        }//if( exclude_peak_energies specified )

        selected_peaks.push_back( p );
      }
    }//for( peaks )

    if( selected_peaks.empty() )
    {
      string sources_csv;
      for( const string &s : target_sources )
        sources_csv += (sources_csv.empty() ? "" : ", ") + s;

      string msg = "No peaks found with the specified sources (" + sources_csv + ")";
      if( !energy_ranges.empty() )
        msg += " within the specified energy range(s)";
      msg += " assigned to them. Please add analysis peaks for these sources and try again.";

      throw runtime_error( msg );
    }
  }//if( has_peak_energies ) / else
  
  // Verify all peaks have sources
  for( const auto &p : selected_peaks )
  {
    if( !p->parentNuclide() && !p->reaction() )
      throw runtime_error( "Peak at " + to_string(p->mean()) + " keV has no source assigned" );
  }
  
  // Build nuclide age map
  map<const SandiaDecay::Nuclide *, double> nuclide_ages;
  for( const auto &p : selected_peaks )
  {
    if( p->parentNuclide() )
    {
      const SandiaDecay::Nuclide *nuc = p->parentNuclide();
      if( nuclide_ages.find(nuc) == nuclide_ages.end() )
      {
        // Check if age was specified in input
        auto it = nuclide_ages_input.find( nuc->symbol );
        if( it != nuclide_ages_input.end() )
          nuclide_ages[nuc] = it->second;
        else
          nuclide_ages[nuc] = getDefaultNuclideAge( nuc );
      }
    }
  }
  
  // Convert selected peaks to GenericPeakInfo
  vector<GenericPeakInfo> peak_infos;
  const double add_uncert_val = addUncertToValue( add_uncert );
  const double match_tol_sigma = 2.35482 * match_tolerance;
  const float foreground_live_time = foreground ? foreground->live_time() : 1.0f;
  const float background_live_time = background ? background->live_time() : 1.0f;
  
  for( const shared_ptr<const PeakDef> &p : selected_peaks )
  {
    GenericPeakInfo peak;
    peak.m_mean = peak.m_energy = p->mean();
    peak.m_fwhm = p->gausPeak() ? p->fwhm() : (2.35482 * 0.25 * p->roiWidth());
    peak.m_counts = p->amplitude();
    peak.m_counts_uncert = p->amplitudeUncert();
    peak.m_base_rel_eff_uncert = add_uncert_val;
    
    // Use gamma particle energy instead of fit mean if available
    if( p->gammaParticleEnergy() > 0.0 )
      peak.m_energy = p->gammaParticleEnergy();
    
    // Background subtraction
    if( background_subtract )
    {
      const double sigma = p->gausPeak() ? p->sigma() : 0.25 * p->roiWidth();
      const double scale = foreground_live_time / background_live_time;
      
      double back_counts = 0.0, back_uncert_2 = 0.0;
      for( const shared_ptr<const PeakDef> &back_peak : background_peaks )
      {
        if( fabs(back_peak->mean() - p->mean()) < sigma )
        {
          back_counts += scale * back_peak->peakArea();
          back_uncert_2 += scale * scale * back_peak->peakAreaUncert() * back_peak->peakAreaUncert();
        }
      }
      
      if( back_counts > 0.0 )
      {
        peak.m_counts -= back_counts;
        peak.m_counts_uncert = sqrt( peak.m_counts_uncert * peak.m_counts_uncert + back_uncert_2 );
      }
      
      if( peak.m_counts <= 0.0 )
        continue;  // Skip peaks with zero or negative counts after background subtraction
    }//if( background_subtract )
    
    peak_infos.push_back( peak );
  }//for( selected peaks )
  
  if( peak_infos.empty() )
    throw runtime_error( "No valid peaks after processing" );
  
  // Build nuclides_to_match_to
  vector<SandiaDecayNuc> nuclides_to_match_to;
  for( const auto &p : selected_peaks )
  {
    SandiaDecayNuc nuc;
    bool already_have = false;
    
    if( p->parentNuclide() )
    {
      for( const auto &existing : nuclides_to_match_to )
      {
        if( existing.nuclide == p->parentNuclide() )
        {
          already_have = true;
          break;
        }
      }
      
      if( !already_have )
      {
        nuc.nuclide = p->parentNuclide();
        auto it = nuclide_ages.find( nuc.nuclide );
        nuc.age = (it != nuclide_ages.end()) ? it->second : getDefaultNuclideAge(nuc.nuclide);
        nuc.correct_for_decay_during_meas = false;
        nuclides_to_match_to.push_back( nuc );
      }
    }
    else if( p->reaction() )
    {
      for( const auto &existing : nuclides_to_match_to )
      {
        if( existing.reaction == p->reaction() )
        {
          already_have = true;
          break;
        }
      }
      
      if( !already_have )
      {
        nuc.reaction = p->reaction();
        nuclides_to_match_to.push_back( nuc );
      }
    }
  }//for( selected peaks )
  
  // Fill in nuclide info using PeakCsvInput::fill_in_nuclide_info
  {
    vector<PeakCsvInput::NucAndAge> isotopes;
    for( const auto &n : nuclides_to_match_to )
    {
      if( n.nuclide )
        isotopes.emplace_back( n.nuclide->symbol, n.age, n.correct_for_decay_during_meas );
      else if( n.reaction )
        isotopes.emplace_back( n.reaction->name(), -1.0, false );
    }
    
    const float meas_time = foreground ? foreground->real_time() : -1.0f;
    const double tol = std::max( match_tol_sigma, 0.0001 );
    
    const PeakCsvInput::NucMatchResults matched_res = PeakCsvInput::fill_in_nuclide_info(
      peak_infos,
      PeakCsvInput::NucDataSrc::SandiaDecay,
      {},  // energy ranges
      isotopes,
      tol,
      {},  // excluded energies
      meas_time
    );
    
    if( !matched_res.unused_isotopes.empty() )
    {
      string unused_nucs;
      for( size_t i = 0; i < matched_res.unused_isotopes.size(); ++i )
        unused_nucs += string(i ? ", " : "") + matched_res.unused_isotopes[i];
      throw runtime_error( "Failed to match nuclide(s) to peaks: " + unused_nucs );
    }
    
    peak_infos = matched_res.peaks_matched;
  }
  
  // Build RelEffInput
  RelEffInput input;
  input.peaks = peak_infos;
  input.eqn_form = eqn_form;
  input.eqn_order = eqn_order;
  input.use_ceres_to_fit_eqn = (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel);
  
  // Setup physical model if needed
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    input.phys_model_use_hoerl = true;
    input.phys_model_detector = spec->detector();
    
    // If no detector loaded, we'll skip the physical model - can't proceed without a valid DRF
    if( !input.phys_model_detector || !input.phys_model_detector->isValid() )
    {
      throw runtime_error( "FramPhysicalModel requires a valid detector efficiency function (DRF) - none is loaded" );
    }
    
    MaterialDB *materialDb = interspec->materialDataBase();
    
    // Helper lambda to parse a shielding JSON object
    auto parseShielding = [materialDb]( const json &shield, const string &context ) -> shared_ptr<RelActCalc::PhysicalModelShieldInput>
    {
      auto shield_input = make_shared<RelActCalc::PhysicalModelShieldInput>();
      
      const bool has_material = shield.contains("material") && shield["material"].is_string();
      const bool has_atomic_number = shield.contains("atomic_number") && shield["atomic_number"].is_number();
      const bool has_areal_density = shield.contains("areal_density") && shield["areal_density"].is_number();
      const bool has_thickness = shield.contains("thickness") && shield["thickness"].is_string();
      
      if( has_material && has_atomic_number )
        throw runtime_error( context + ": 'material' and 'atomic_number' are mutually exclusive" );
      
      if( !has_material && !has_atomic_number )
        throw runtime_error( context + ": requires either 'material' or 'atomic_number'" );
      
      if( has_thickness && !has_material )
        throw runtime_error( context + ": 'thickness' can only be used with 'material'" );
      
      if( has_thickness && has_areal_density )
        throw runtime_error( context + ": 'thickness' and 'areal_density' are mutually exclusive" );
      
      // Parse material or atomic number
      if( has_material )
      {
        const string mat_name = shield["material"].get<string>();
        if( !materialDb )
          throw runtime_error( "Material database not available" );
        
        const Material *mat = materialDb->material( mat_name );
        if( !mat )
          throw runtime_error( context + ": Unknown material '" + mat_name + "'" );
        
        shield_input->material = make_shared<Material>( *mat );
        shield_input->atomic_number = 0.0;
        
        // Get areal density from either areal_density or thickness
        if( has_areal_density )
        {
          shield_input->areal_density = shield["areal_density"].get<double>() * PhysicalUnits::g / PhysicalUnits::cm2;
        }
        else if( has_thickness )
        {
          const string thickness_str = shield["thickness"].get<string>();
          try
          {
            const double thickness = PhysicalUnits::stringToDistance( thickness_str );
            // areal_density = thickness * density
            shield_input->areal_density = thickness * mat->density;
          }
          catch( exception &e )
          {
            throw runtime_error( context + ": Failed to parse thickness '" + thickness_str + "': " + e.what() );
          }
        }
        else
        {
          throw runtime_error( context + ": requires either 'areal_density' or 'thickness' when using 'material'" );
        }
      }
      else // has_atomic_number
      {
        shield_input->atomic_number = shield["atomic_number"].get<double>();
        
        if( !has_areal_density )
          throw runtime_error( context + ": 'areal_density' is required when using 'atomic_number'" );
        
        shield_input->areal_density = shield["areal_density"].get<double>() * PhysicalUnits::g / PhysicalUnits::cm2;
      }
      
      // Fit options
      if( shield.contains("fit_areal_density") && shield["fit_areal_density"].is_boolean() )
        shield_input->fit_areal_density = shield["fit_areal_density"].get<bool>();
      
      if( shield.contains("fit_atomic_number") && shield["fit_atomic_number"].is_boolean() )
      {
        if( shield_input->material )
          throw runtime_error( context + ": Cannot fit atomic_number when using a material" );
        shield_input->fit_atomic_number = shield["fit_atomic_number"].get<bool>();
      }
      
      return shield_input;
    };//parseShielding lambda
    
    // Parse self-attenuation shielding if provided
    if( params.contains("self_atten_shielding") && params["self_atten_shielding"].is_object() )
    {
      input.phys_model_self_atten = parseShielding( params["self_atten_shielding"], "self_atten_shielding" );
    }//if( self_atten_shielding )
    
    // Parse external shieldings if provided
    if( params.contains("external_shieldings") && params["external_shieldings"].is_array() )
    {
      size_t idx = 0;
      for( const auto &shield : params["external_shieldings"] )
      {
        if( !shield.is_object() )
          continue;
        
        input.phys_model_external_attens.push_back( 
          parseShielding( shield, "external_shieldings[" + to_string(idx) + "]" ) 
        );
        idx++;
      }//for( external shieldings )
    }//if( external_shieldings )
  }//if( FramPhysicalModel )
  
  // Solve
  RelEffSolution solution = solve_relative_efficiency( input );
  
  // Build result JSON
  json result;
  
  if( solution.m_status != ManualSolutionStatus::Success )
  {
    result["success"] = false;
    result["error"] = solution.m_error_message;
    if( !solution.m_warnings.empty() )
      result["warnings"] = solution.m_warnings;
    return result;
  }
  
  result["success"] = true;

  // If sources were specified, return the peak energies actually used
  if( has_sources )
  {
    json used_energies = json::array();
    for( const shared_ptr<const PeakDef> &p : selected_peaks )
      used_energies.push_back( p->mean() );
    result["peak_energies_used"] = used_energies;
  }

  // Sources with relative activities and mass fractions
  json sources_arr = json::array();
  for( const auto &rel_act : solution.m_rel_activities )
  {
    json src;
    src["name"] = rel_act.m_isotope;
    src["rel_activity"] = rel_act.m_rel_activity;
    src["rel_activity_uncert"] = rel_act.m_rel_activity_uncert;
    
    // Try to get mass fraction
    try
    {
      const double mass_frac = solution.mass_fraction( rel_act.m_isotope );
      src["mass_fraction"] = mass_frac;
    }
    catch( ... )
    {
      // Mass fraction not available (e.g., for reactions)
    }
    
    sources_arr.push_back( src );
  }
  result["sources"] = sources_arr;
  
  // Activity ratios
  json activity_ratios = json::array();
  json mass_ratios = json::array();
  
  for( size_t i = 0; i < solution.m_rel_activities.size(); ++i )
  {
    for( size_t j = i + 1; j < solution.m_rel_activities.size(); ++j )
    {
      try
      {
        json ratio_obj;
        ratio_obj["numerator"] = solution.m_rel_activities[i].m_isotope;
        ratio_obj["denominator"] = solution.m_rel_activities[j].m_isotope;
        ratio_obj["ratio"] = solution.activity_ratio( i, j );
        ratio_obj["ratio_uncert"] = solution.activity_ratio_uncert( i, j );
        activity_ratios.push_back( ratio_obj );
      }
      catch( ... )
      {
        // Skip if ratio can't be computed
      }
      
      try
      {
        const double mass_frac_i = solution.mass_fraction( solution.m_rel_activities[i].m_isotope );
        const double mass_frac_j = solution.mass_fraction( solution.m_rel_activities[j].m_isotope );
        
        if( mass_frac_j > 0 )
        {
          json mass_ratio_obj;
          mass_ratio_obj["numerator"] = solution.m_rel_activities[i].m_isotope;
          mass_ratio_obj["denominator"] = solution.m_rel_activities[j].m_isotope;
          mass_ratio_obj["ratio"] = mass_frac_i / mass_frac_j;
          // Note: mass ratio uncertainty would require propagation from activity uncertainties
          mass_ratios.push_back( mass_ratio_obj );
        }
      }
      catch( ... )
      {
        // Skip if mass fraction not available
      }
    }//for( j )
  }//for( i )
  
  result["activity_ratios"] = activity_ratios;
  result["mass_ratios"] = mass_ratios;
  
  // Peak information with observed vs fit efficiency
  json peaks_arr = json::array();
  
  for( const auto &peak : solution.m_input.peaks )
  {
    json peak_obj;
    peak_obj["energy"] = peak.m_energy;
    peak_obj["counts"] = peak.m_counts;
    peak_obj["counts_uncert"] = peak.m_counts_uncert;
    
    // Get source name
    if( !peak.m_source_gammas.empty() )
      peak_obj["source"] = peak.m_source_gammas[0].m_isotope;
    
    // Calculate observed efficiency: eff = counts / (rel_activity * yield)
    double total_expected = 0.0;
    for( const auto &line : peak.m_source_gammas )
    {
      // Find the relative activity for this source
      for( const auto &rel_act : solution.m_rel_activities )
      {
        if( rel_act.m_isotope == line.m_isotope )
        {
          total_expected += rel_act.m_rel_activity * line.m_yield;
          break;
        }
      }
    }
    
    if( total_expected > 0.0 )
    {
      const double obs_eff = peak.m_counts / total_expected;
      const double obs_eff_uncert = peak.m_counts_uncert / total_expected;
      const double fit_eff = solution.relative_efficiency( peak.m_energy );
      
      peak_obj["observed_efficiency"] = obs_eff;
      peak_obj["observed_efficiency_uncert"] = obs_eff_uncert;
      peak_obj["fit_efficiency"] = fit_eff;
      
      // Calculate residual sigma
      if( obs_eff_uncert > 0.0 )
      {
        const double residual_sigma = fabs(obs_eff - fit_eff) / obs_eff_uncert;
        peak_obj["residual_sigma"] = residual_sigma;
      }
    }
    
    peaks_arr.push_back( peak_obj );
  }//for( peaks )
  
  result["peaks"] = peaks_arr;
  
  // Fit quality
  json quality;
  quality["chi2"] = solution.m_chi2;
  quality["dof"] = solution.m_dof;
  quality["chi2_per_dof"] = (solution.m_dof > 0) ? (solution.m_chi2 / solution.m_dof) : -1.0;
  result["fit_quality"] = quality;
  
  // Relative efficiency equation text
  try
  {
    result["rel_eff_equation"] = solution.rel_eff_eqn_txt( false );
  }
  catch( ... )
  {
    result["rel_eff_equation"] = nullptr;
  }
  
  // Warnings
  result["warnings"] = solution.m_warnings;
  
  // Save to state if requested
  result["saved_to_state"] = false;
  if( save_to_state )
  {
    // Build GUI state and save
    RelActManualGui::GuiState gui_state;
    gui_state.m_relEffEqnFormIndex = eqn_form;
    gui_state.m_relEffEqnOrderIndex = static_cast<int>(eqn_order);
    gui_state.m_nucDataSrcIndex = PeakCsvInput::NucDataSrc::SandiaDecay;
    gui_state.m_matchToleranceValue = static_cast<float>(match_tolerance);
    gui_state.m_addUncertIndex = add_uncert;
    gui_state.m_backgroundSubtract = background_subtract;
    gui_state.m_physModelUseHoerl = true;
    
    // Copy nuclide ages
    for( const auto &[nuc, age] : nuclide_ages )
    {
      if( nuc )
        gui_state.nucAge[nuc->symbol] = age;
    }
    
    // Serialize and save to SpecMeas
    try
    {
      unique_ptr<rapidxml::xml_document<char>> doc( new rapidxml::xml_document<char>() );
      gui_state.serialize( doc.get() );
      
      // Update GUI if open (do this before moving doc)
      RelActManualGui *gui = interspec->relActManualWidget( false );
      if( gui && doc->first_node() )
        gui->deSerialize( doc->first_node() );
      
      spec->setRelActManualGuiState( std::move(doc) );
      
      result["saved_to_state"] = true;
    }
    catch( exception &e )
    {
      result["warnings"].push_back( string("Failed to save state: ") + e.what() );
    }
  }//if( save_to_state )
  
  return result;
}//executePeakBasedRelativeEfficiency


/** Helper function to convert GuiState to JSON. */
json guiStateToJson( const RelActManualGui::GuiState &gui_state )
{
  json state;
  state["eqn_form"] = eqnFormToString( gui_state.m_relEffEqnFormIndex );
  state["eqn_order"] = gui_state.m_relEffEqnOrderIndex;
  state["match_tolerance"] = gui_state.m_matchToleranceValue;
  state["additional_uncertainty"] = addUncertToString( gui_state.m_addUncertIndex );
  state["background_subtract"] = gui_state.m_backgroundSubtract;
  
  // Nuclide ages
  json nuc_ages;
  for( const auto &[name, age] : gui_state.nucAge )
  {
    nuc_ages[name] = to_string(age) + " s";
  }
  state["nuclide_ages"] = nuc_ages;
  
  // Decay correction flags
  json decay_correct;
  for( const auto &[name, correct] : gui_state.nucDecayCorrect )
  {
    decay_correct[name] = correct;
  }
  state["nuclide_decay_correct"] = decay_correct;
  
  return state;
}//guiStateToJson


nlohmann::json executeGetRelActManualState(
  const nlohmann::json& params,
  InterSpec* interspec
)
{
  if( !interspec )
    throw runtime_error( "InterSpec instance required" );
  
  json result;
  result["success"] = true;
  result["has_state"] = false;
  result["source_from_gui"] = false;
  
  // Try to get state from GUI first
  RelActManualGui *gui = interspec->relActManualWidget( false );
  if( gui )
  {
    try
    {
      // Serialize the GUI state to XML, then deserialize to GuiState
      rapidxml::xml_document<char> doc;
      rapidxml::xml_node<char> *node = gui->serialize( &doc );
      
      if( node )
      {
        RelActManualGui::GuiState gui_state;
        gui_state.deSerialize( node );
        
        result["has_state"] = true;
        result["source_from_gui"] = true;
        result["state"] = guiStateToJson( gui_state );
        return result;
      }
    }
    catch( ... )
    {
      // Fall through to try SpecMeas
    }
  }//if( gui )
  
  // Try to get state from SpecMeas
  shared_ptr<SpecMeas> spec = interspec->measurment( SpecUtils::SpectrumType::Foreground );
  if( spec )
  {
    rapidxml::xml_document<char> *xml_state = spec->relActManualGuiState();
    if( xml_state && xml_state->first_node() )
    {
      try
      {
        RelActManualGui::GuiState gui_state;
        gui_state.deSerialize( xml_state->first_node() );
        
        result["has_state"] = true;
        result["source_from_gui"] = false;
        result["state"] = guiStateToJson( gui_state );
      }
      catch( exception &e )
      {
        result["error"] = string("Failed to parse stored state: ") + e.what();
      }
    }
  }//if( spec )
  
  return result;
}//executeGetRelActManualState


}//namespace RelActManualTool
}//namespace LlmTools

#endif // USE_LLM_INTERFACE

