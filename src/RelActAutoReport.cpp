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

#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PeakDef.h"
#include "InterSpec/BatchInfoLog.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/DateTime.h"
#include "SpecUtils/EnergyCalibration.h"

#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/ReactionGamma.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <set>

using namespace std;

namespace RelActAutoReport
{

nlohmann::json solution_to_json(const RelActCalcAuto::RelActAutoSolution& solution, 
                                const ReportConfig& config)
{
  nlohmann::json json_data;
  
  // Basic solution status and info
  json_data["status"] = {
    {"success", solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success},
    {"status_code", static_cast<int>(solution.m_status)},
    {"error_message", solution.m_error_message}
  };
  
  json_data["spectrum_title"] = solution.m_options.spectrum_title;
  
  // Analysis results (only if successful)
  if (solution.m_status == RelActCalcAuto::RelActAutoSolution::Status::Success) {
    json_data["chi2"] = solution.m_chi2;
    json_data["dof"] = solution.m_dof;
    json_data["chi2_per_dof"] = solution.m_dof > 0 ? solution.m_chi2 / solution.m_dof : 0.0;
    
    // Relative efficiency equations
    json_data["rel_eff_curves"] = nlohmann::json::array();
    for (size_t i = 0; i < solution.m_rel_eff_forms.size(); ++i) {
      nlohmann::json curve_info;
      curve_info["index"] = i;
      curve_info["name"] = solution.m_options.rel_eff_curves[i].name;
      
      try {
        curve_info["equation_text"] = solution.rel_eff_txt(false, i);
        curve_info["equation_html"] = solution.rel_eff_txt(true, i);
        curve_info["js_function"] = solution.rel_eff_eqn_js_function(i);
      } catch (const std::exception& e) {
        curve_info["equation_error"] = e.what();
      }
      
      json_data["rel_eff_curves"].push_back(curve_info);
    }
    
    // Relative activities by curve
    json_data["relative_activities"] = nlohmann::json::array();
    for (size_t curve_idx = 0; curve_idx < solution.m_rel_activities.size(); ++curve_idx) {
      nlohmann::json curve_activities;
      curve_activities["curve_index"] = curve_idx;
      curve_activities["nuclides"] = nlohmann::json::array();
      
      double sum_rel_mass = 0.0;
      for (const auto& act : solution.m_rel_activities[curve_idx]) {
        const SandiaDecay::Nuclide* nuc = RelActCalcAuto::nuclide(act.source);
        if (nuc) {
          sum_rel_mass += act.rel_activity / nuc->activityPerGram();
        }
      }
      
      for (const auto& act : solution.m_rel_activities[curve_idx]) {
        nlohmann::json nuc_info;
        nuc_info["name"] = act.name();
        nuc_info["rel_activity"] = act.rel_activity;
        nuc_info["rel_activity_uncertainty"] = act.rel_activity_uncertainty;
        nuc_info["age"] = act.age;
        nuc_info["age_uncertainty"] = act.age_uncertainty;
        nuc_info["age_was_fit"] = act.age_was_fit;
        
        const SandiaDecay::Nuclide* nuc = RelActCalcAuto::nuclide(act.source);
        if (nuc) {
          const double rel_mass = act.rel_activity / nuc->activityPerGram();
          nuc_info["total_mass_fraction"] = sum_rel_mass > 0 ? (rel_mass / sum_rel_mass) : 0.0;
          
          // Enrichment calculations
          try {
            pair<double,optional<double>> enrich = solution.mass_enrichment_fraction(nuc, curve_idx );
            nuc_info["enrichment"] = enrich.first;
            if( enrich.second.has_value() )
            {
              nuc_info["has_enrichment_uncert"] = true;
              nuc_info["enrichment_uncert"] = *enrich.second;
              nuc_info["enrichment_minus_2sigma"] = enrich.first - 2.0*enrich.second.value();
              nuc_info["enrichment_plus_2sigma"] = enrich.first + 2.0*enrich.second.value();
            }else
            {
              nuc_info["has_enrichment_uncert"] = false;
            }
          } catch (const std::exception& e) {
            nuc_info["enrichment_error"] = e.what();
          }
          
          // Detector counts
          try {
            nuc_info["detector_counts"] = solution.nuclide_counts(act.source, curve_idx);
          } catch (const std::exception& e) {
            nuc_info["counts_error"] = e.what();
          }
        }
        
        curve_activities["nuclides"].push_back(nuc_info);
      }
      
      json_data["relative_activities"].push_back(curve_activities);
    }
    
    // Mass and activity ratios
    json_data["ratios"] = nlohmann::json::array();
    for (size_t curve_idx = 0; curve_idx < solution.m_rel_activities.size(); ++curve_idx) {
      const auto& rel_activities = solution.m_rel_activities[curve_idx];
      
      nlohmann::json curve_ratios;
      curve_ratios["curve_index"] = curve_idx;
      curve_ratios["pairs"] = nlohmann::json::array();
      
      for (size_t i = 1; i < rel_activities.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
          const auto& nuc_i = rel_activities[i];
          const auto& nuc_j = rel_activities[j];
          
          nlohmann::json ratio_pair;
          ratio_pair["numerator"] = nuc_i.name();
          ratio_pair["denominator"] = nuc_j.name();
          ratio_pair["activity_ratio"] = nuc_i.rel_activity / nuc_j.rel_activity;
          
          const SandiaDecay::Nuclide* nuc_i_ptr = RelActCalcAuto::nuclide(nuc_i.source);
          const SandiaDecay::Nuclide* nuc_j_ptr = RelActCalcAuto::nuclide(nuc_j.source);
          
          if (nuc_i_ptr && nuc_j_ptr) {
            const double mass_i = nuc_i.rel_activity / nuc_i_ptr->activityPerGram();
            const double mass_j = nuc_j.rel_activity / nuc_j_ptr->activityPerGram();
            ratio_pair["mass_ratio"] = mass_i / mass_j;
          }
          
          // Activity ratio uncertainty
          try {
            const double uncert = solution.activity_ratio_uncertainty(nuc_i.source, curve_idx, 
                                                                     nuc_j.source, curve_idx);
            ratio_pair["activity_ratio_uncertainty"] = uncert;
            ratio_pair["activity_ratio_uncertainty_percent"] = 
              100.0 * uncert / ratio_pair["activity_ratio"].get<double>();
          } catch (...) {
            // Uncertainty calculation failed
          }
          
          curve_ratios["pairs"].push_back(ratio_pair);
        }
      }
      
      json_data["ratios"].push_back(curve_ratios);
    }
    
    // Plutonium correction data
    if (!solution.m_corrected_pu.empty()) {
      assert( solution.m_corrected_pu.size() == solution.m_uncorrected_pu.size() );
      
      json_data["plutonium_corrections"] = nlohmann::json::array();
      for (size_t i = 0; i < solution.m_corrected_pu.size(); ++i) {
        if (solution.m_corrected_pu[i]) {
          const auto& pu_corr = *solution.m_corrected_pu[i];
          nlohmann::json pu_data;
          pu_data["curve_index"] = i;
          pu_data["mass_fractions"] = {
            {"pu238", pu_corr.pu238_mass_frac},
            {"pu239", pu_corr.pu239_mass_frac},
            {"pu240", pu_corr.pu240_mass_frac},
            {"pu241", pu_corr.pu241_mass_frac},
            {"pu242", pu_corr.pu242_mass_frac},
            {"IsWithinRange", pu_corr.is_within_range},
            {"pu242Uncert", pu_corr.pu242_uncert}
          };
          if( (i < solution.m_uncorrected_pu.size()) && solution.m_uncorrected_pu[i] ) {
            pu_data["uncorrected_mass_fractions"] = {
              {"pu238", solution.m_uncorrected_pu[i]->pu238_rel_mass},
              {"pu239", solution.m_uncorrected_pu[i]->pu239_rel_mass},
              {"pu240", solution.m_uncorrected_pu[i]->pu240_rel_mass},
              {"pu241", solution.m_uncorrected_pu[i]->pu241_rel_mass},
              {"puOther", solution.m_uncorrected_pu[i]->other_pu_mass},
              {"Age", PhysicalUnits::printToBestTimeUnits(solution.m_uncorrected_pu[i]->pu_age,6) },
              {"AgeSeconds", solution.m_uncorrected_pu[i]->pu_age }
            };
          }
          json_data["plutonium_corrections"].push_back(pu_data);
        }
      }
    }
  }
  
  // Peak information
  if( config.include_peak_details && !solution.m_fit_peaks.empty() )
  {
    json_data["peaks"] = nlohmann::json::array();

    vector<pair<PeakDef,size_t>> peak_rel_eff;
    for( size_t rel_eff_index = 0; rel_eff_index < solution.m_fit_peaks_for_each_curve.size(); ++rel_eff_index )
    {
      for( const PeakDef &p : solution.m_fit_peaks_for_each_curve[rel_eff_index] )
        peak_rel_eff.emplace_back( p, rel_eff_index );
    }

    // Now grab the free-floating peaks; we will just assign them releff index 0 for the moment
    for( const PeakDef &peak : solution.m_fit_peaks )
    {
      if( !peak.parentNuclide() && !peak.xrayElement() && !peak.reaction() )
        peak_rel_eff.emplace_back( peak, size_t(0) );
    }

    std::sort( begin(peak_rel_eff), end(peak_rel_eff), []( const pair<PeakDef,size_t> &lhs, const pair<PeakDef,size_t> &rhs ) -> bool {
      return lhs.first.mean() < rhs.first.mean();
    } );


    for (const pair<PeakDef,size_t> &peak_releff_index : peak_rel_eff )
    {
      const PeakDef &peak = peak_releff_index.first;
      const size_t rel_eff_index = peak_releff_index.second;

      nlohmann::json peak_info;
      peak_info["energy"] = peak.mean();
      peak_info["amplitude"] = peak.amplitude();
      peak_info["amplitude_uncertainty"] = peak.amplitudeUncert();
      peak_info["amplitude_uncertainty_percent"] = 
        peak.amplitude() != 0 ? 100.0 * peak.amplitudeUncert() / peak.amplitude() : 0.0;
      peak_info["rel_eff_index"] = static_cast<int>(rel_eff_index);

      // Nuclide information
      const SandiaDecay::Nuclide* nuc = peak.parentNuclide();
      const SandiaDecay::Element* el = peak.xrayElement();
      const ReactionGamma::Reaction* reaction = peak.reaction();
      
      if (nuc) {
        peak_info["nuclide"] = nuc->symbol;
        peak_info["source_type"] = "nuclide";
      } else if (el) {
        peak_info["nuclide"] = el->symbol;
        peak_info["source_type"] = "element";
      } else if (reaction) {
        peak_info["nuclide"] = reaction->name();
        peak_info["source_type"] = "reaction";
      }
      
      // Continuum information
      if (peak.continuum()) {
        peak_info["continuum"] = {
          {"type", PeakContinuum::offset_type_label_tr(peak.continuum()->type())},
          {"lower_energy", peak.continuum()->lowerEnergy()},
          {"upper_energy", peak.continuum()->upperEnergy()}
        };
      }
      
      json_data["peaks"].push_back(peak_info);
    }
  }
  
  // Energy calibration adjustments
  json_data["energy_calibration"] = {
    {"was_fit", solution.m_options.fit_energy_cal},
    {"adjustments", nlohmann::json::array()}
  };
  
  for (size_t i = 0; i < solution.m_energy_cal_adjustments.size(); ++i) {
    if (solution.m_fit_energy_cal[i]) {
      nlohmann::json adj;
      adj["parameter_index"] = i;
      adj["raw_value"] = solution.m_energy_cal_adjustments[i];
      adj["was_fit"] = solution.m_fit_energy_cal[i];
      
      // Convert to physical units
      double physical_value = (solution.m_energy_cal_adjustments[i] / 
                              RelActCalcAuto::RelActAutoSolution::sm_energy_par_offset - 1.0);
      
      const size_t num_chan = solution.m_spectrum && solution.m_spectrum->energy_calibration() 
                             ? solution.m_spectrum->energy_calibration()->num_channels() : 1;
      
      if (i == 0) {
        physical_value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple;
        adj["type"] = "offset";
        adj["units"] = "keV";
      } else if (i == 1) {
        physical_value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / num_chan;
        adj["type"] = "gain";
        adj["units"] = "keV/channel";
      } else if (i == 2) {
        physical_value *= RelActCalcAuto::RelActAutoSolution::sm_energy_cal_multiple / (num_chan * num_chan);
        adj["type"] = "quadratic";
        adj["units"] = "keV/channel²";
      }
      
      adj["physical_value"] = physical_value;
      json_data["energy_calibration"]["adjustments"].push_back(adj);
    }
  }
  
  // Analysis options
  json_data["options"] = {
    {"fit_energy_cal", solution.m_options.fit_energy_cal},
    {"fwhm_form", RelActCalcAuto::to_str(solution.m_options.fwhm_form)},
    {"fwhm_estimation_method", RelActCalcAuto::to_str(solution.m_options.fwhm_estimation_method)},
    {"skew_type", static_cast<int>(solution.m_options.skew_type)},
    {"additional_br_uncert", solution.m_options.additional_br_uncert}
  };
  
  // Add curve-specific options
  json_data["curve_options"] = nlohmann::json::array();
  for (const auto& curve : solution.m_options.rel_eff_curves) {
    nlohmann::json curve_opts;
    curve_opts["name"] = curve.name;
    curve_opts["nucs_of_el_same_age"] = curve.nucs_of_el_same_age;
    curve_opts["rel_eff_eqn_type"] = RelActCalc::to_str(curve.rel_eff_eqn_type);
    if (curve.rel_eff_eqn_type != RelActCalc::RelEffEqnForm::FramPhysicalModel) {
      curve_opts["rel_eff_eqn_order"] = curve.rel_eff_eqn_order;
    }
    json_data["curve_options"].push_back(curve_opts);
  }
  
  // Warnings
  if (config.include_warnings && !solution.m_warnings.empty()) {
    json_data["warnings"] = solution.m_warnings;
  }
  
  // Timing information
  if (config.include_timing_info) {
    json_data["timing"] = {
      {"function_eval_solution", solution.m_num_function_eval_solution},
      {"function_eval_total", solution.m_num_function_eval_total},
      {"microseconds_eval", solution.m_num_microseconds_eval},
      {"microseconds_in_eval", solution.m_num_microseconds_in_eval},
      {"duration_formatted", format_duration(solution.m_num_microseconds_eval)}
    };
  }
  
  // Timestamps
  json_data["timestamps"] = {
    {"utc", SpecUtils::to_common_string(chrono::time_point_cast<chrono::microseconds>(chrono::system_clock::now()), true)},
    {"local", ""} // TODO: Implement local time formatting
  };
  
  // Multiple rel eff curves flag
  json_data["have_multiple_rel_eff"] = solution.m_rel_eff_forms.size() > 1;
  
  // Add spectrum histograms using BatchInfoLog function
  if (config.include_spectrum_chart) {
    // Add foreground spectrum data
    if (solution.m_foreground) {
      // Use empty set for sample numbers since we don't have that info from RelActAutoSolution
      set<int> sample_numbers;
      const string filename = solution.m_foreground->title();
      
      auto &spec_obj = json_data["foreground"];
      
      const deque<std::shared_ptr<const PeakDef>> * const fit_peaks = nullptr;
      BatchInfoLog::add_hist_to_json( json_data, false, solution.m_foreground,  nullptr, sample_numbers, filename,  fit_peaks );
    }
    
    // Add background spectrum data
    if (solution.m_background) {
      set<int> sample_numbers;
      const string filename = solution.m_background->title();
      
      auto &spec_obj = json_data["background"];
      
      const deque<std::shared_ptr<const PeakDef>> * const fit_peaks = nullptr;
      
      BatchInfoLog::add_hist_to_json(spec_obj, true, solution.m_background, nullptr, sample_numbers, filename, fit_peaks );
    }
    
    
    if (solution.m_spectrum ){
      set<int> sample_numbers;
      const string filename = solution.m_background->title();
      
      auto &spec_obj = json_data["fit_spectrum"];
      deque<std::shared_ptr<const PeakDef>> fit_peaks;
      for( const auto &peak : solution.m_fit_peaks )
        fit_peaks.push_back( std::make_shared<const PeakDef>(peak) );
      
      BatchInfoLog::add_hist_to_json(spec_obj, false, solution.m_spectrum, nullptr, sample_numbers, filename, &fit_peaks );


      auto peaks_to_json = [&]( const std::vector<PeakDef> &peaks ) -> nlohmann::json {
        vector<std::shared_ptr<const PeakDef>> peaks_deque;
        for( const auto &peak : peaks )
          peaks_deque.push_back( std::make_shared<const PeakDef>(peak) );
        
        const string peaks_json = PeakDef::peak_json( peaks_deque, solution.m_spectrum, Wt::WColor(0,51,255), 255 ); //will [{},{},...]
        return nlohmann::json::parse( peaks_json );
      };

      json_data["fit_peaks"] = peaks_to_json( solution.m_fit_peaks );
      json_data["fit_peaks_in_spectrums_cal"] = peaks_to_json( solution.m_fit_peaks_in_spectrums_cal );

      auto fit_peaks_for_each_curve = nlohmann::json::array();
      for( const auto &peaks : solution.m_fit_peaks_for_each_curve )
        fit_peaks_for_each_curve.push_back( peaks_to_json( peaks ) );
      spec_obj["fit_peaks_for_each_curve"] = fit_peaks_for_each_curve;

      auto fit_peaks_in_spectrums_cal_for_each_curve = nlohmann::json::array();
      for( const auto &peaks : solution.m_fit_peaks_in_spectrums_cal_for_each_curve )
        fit_peaks_in_spectrums_cal_for_each_curve.push_back( peaks_to_json( peaks ) );
      json_data["fit_peaks_in_spectrums_cal_for_each_curve"] = fit_peaks_in_spectrums_cal_for_each_curve;
    }//if( solution.m_spectrum )
  }//if( config.include_spectrum_chart )
  
  return json_data;
}

std::string generate_report(const RelActCalcAuto::RelActAutoSolution& solution,
                           const ReportConfig& config)
{
  auto env = InjaCallbacks::create_environment();
  
  // Convert solution to JSON
  auto json_data = solution_to_json(solution, config);
  
  // Get template
  string template_str;
  if (!config.custom_template.empty()) {
    template_str = config.custom_template;
  } else {
    template_str = get_builtin_template(config.format);
  }
  
  // Render template
  return env.render(template_str, json_data);
}

void generate_report(std::ostream& out,
                    const RelActCalcAuto::RelActAutoSolution& solution,
                    const ReportConfig& config)
{
  out << generate_report(solution, config);
}

// Forward declarations
std::string get_html_template();
std::string get_text_template();

std::string get_builtin_template(ReportFormat format)
{
  switch (format) {
    case ReportFormat::Html:
      return get_html_template();
    case ReportFormat::Text:
      return get_text_template();
    case ReportFormat::Json:
      return "{{ . }}";  // Just output the raw JSON
    default:
      throw std::invalid_argument("Unknown report format");
  }
}

std::string get_builtin_css()
{
  return R"(
    body { font-family: Arial, sans-serif; margin: 20px; }
    .results { margin: 20px 0; }
    .releffeqn { font-family: monospace; background: #f5f5f5; padding: 10px; margin: 10px 0; }
    table { border-collapse: collapse; width: 100%; margin: 15px 0; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #f2f2f2; }
    .resulttable { margin: 20px 0; }
    .warning { color: #d9534f; background: #f2dede; border: 1px solid #ebccd1; padding: 10px; margin: 10px 0; }
    .error { color: #a94442; background: #f2dede; border: 1px solid #ebccd1; padding: 10px; margin: 10px 0; }
    .success { color: #3c763d; background: #dff0d8; border: 1px solid #d6e9c6; padding: 10px; margin: 10px 0; }
    .anatime { font-size: 0.9em; color: #666; margin: 20px 0; }
    .anacomputetime { font-size: 0.9em; color: #666; margin: 10px 0; }
  )";
}

std::string get_builtin_js()
{
  return R"(
    // Basic JavaScript for interactive features in HTML reports
    document.addEventListener('DOMContentLoaded', function() {
      // Add click handlers for tables to make them sortable (optional enhancement)
      var tables = document.querySelectorAll('table.resulttable');
      tables.forEach(function(table) {
        // Future: Add sorting functionality
      });
    });
  )";
}



std::string get_html_template()
{
  return R"(<!DOCTYPE html>
<html>
<head>
    <title>{{ spectrum_title | default("Relative Activity Analysis Report") }}</title>
    <meta charset="utf-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .results { margin: 20px 0; }
        .releffeqn { font-family: monospace; background: #f5f5f5; padding: 10px; margin: 10px 0; }
        table { border-collapse: collapse; width: 100%; margin: 15px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .resulttable { margin: 20px 0; }
        .warning { color: #d9534f; background: #f2dede; border: 1px solid #ebccd1; padding: 10px; margin: 10px 0; }
        .error { color: #a94442; background: #f2dede; border: 1px solid #ebccd1; padding: 10px; margin: 10px 0; }
        .success { color: #3c763d; background: #dff0d8; border: 1px solid #d6e9c6; padding: 10px; margin: 10px 0; }
        .anatime { font-size: 0.9em; color: #666; margin: 20px 0; }
        .anacomputetime { font-size: 0.9em; color: #666; margin: 10px 0; }
    </style>
</head>
<body>
    <h1>{{ spectrum_title | default("Relative Activity Analysis Report") }}</h1>
    
    {% if status.success %}
        <div class="success">
            <strong>Analysis completed successfully</strong><br>
            χ² = {{ "%.3f" | format(chi2) }}, DOF = {{ dof }}, χ²/DOF = {{ "%.3f" | format(chi2_per_dof) }}
        </div>
    {% else %}
        <div class="error">
            <strong>Analysis failed:</strong> {{ status.error_message }}
        </div>
    {% endif %}
    
    {% if rel_eff_curves %}
        <h2>Relative Efficiency Equations</h2>
        {% for curve in rel_eff_curves %}
            <div class="releffeqn">
                Rel. Eff. Eqn{% if have_multiple_rel_eff %} {{ curve.index }}{% endif %}: 
                y = {{ curve.equation_html | safe }}
            </div>
        {% endfor %}
    {% endif %}
    
    {% if plutonium_corrections %}
        <h2>Plutonium Mass Fractions</h2>
        {% for pu_data in plutonium_corrections %}
            <table class="resulttable">
                <caption>Plutonium mass fractions{% if have_multiple_rel_eff %} Rel. Eff. {{ pu_data.curve_index }}{% endif %}</caption>
                <thead>
                    <tr>
                        <th>Nuclide</th>
                        <th>% Pu Mass</th>
                    </tr>
                </thead>
                <tbody>
                    <tr><td>Pu238</td><td>{{ "%.4f" | format(pu_data.mass_fractions.pu238 * 100) }}</td></tr>
                    <tr><td>Pu239</td><td>{{ "%.4f" | format(pu_data.mass_fractions.pu239 * 100) }}</td></tr>
                    <tr><td>Pu240</td><td>{{ "%.4f" | format(pu_data.mass_fractions.pu240 * 100) }}</td></tr>
                    <tr><td>Pu241</td><td>{{ "%.4f" | format(pu_data.mass_fractions.pu241 * 100) }}</td></tr>
                    <tr><td>Pu242 (by corr)</td><td>{{ "%.4f" | format(pu_data.mass_fractions.pu242 * 100) }}</td></tr>
                </tbody>
            </table>
        {% endfor %}
    {% endif %}
    
    {% if relative_activities %}
        <h2>Relative Activities and Mass Fractions</h2>
        {% for curve_data in relative_activities %}
            <table class="resulttable">
                <caption>Relative activities and mass fractions{% if have_multiple_rel_eff %} Rel. Eff. {{ curve_data.curve_index }}{% endif %}</caption>
                <thead>
                    <tr>
                        <th>Nuclide</th>
                        <th>Rel. Act.</th>
                        <th>Total Mass Frac.</th>
                        <th>Enrichment</th>
                        <th>Enrich 2σ</th>
                        <th>Det. Counts</th>
                    </tr>
                </thead>
                <tbody>
                    {% for nuclide in curve_data.nuclides %}
                        <tr>
                            <td>{{ nuclide.name }}</td>
                            <td>{{ nuclide.rel_activity | scientific(6) }} ± {{ nuclide.rel_activity_uncertainty | scientific(6) }}</td>
                            <td>{% if nuclide.total_mass_fraction %}{{ "%.2f" | format(nuclide.total_mass_fraction * 100) }}%{% endif %}</td>
                            <td>{% if nuclide.enrichment %}{{ "%.2f" | format(nuclide.enrichment * 100) }}%{% endif %}</td>
                            <td>
                                {% if nuclide.enrichment_minus_2sigma and nuclide.enrichment_plus_2sigma %}
                                    {{ "%.2f" | format(nuclide.enrichment_minus_2sigma * 100) }}%, {{ "%.2f" | format(nuclide.enrichment_plus_2sigma * 100) }}%
                                {% else %}
                                    --
                                {% endif %}
                            </td>
                            <td>{% if nuclide.detector_counts %}{{ nuclide.detector_counts | scientific(6) }}{% endif %}</td>
                        </tr>
                    {% endfor %}
                </tbody>
            </table>
        {% endfor %}
    {% endif %}
    
    {% if ratios %}
        <h2>Mass and Activity Ratios</h2>
        {% for curve_data in ratios %}
            <table class="resulttable">
                <caption>Mass and Activity Ratios{% if have_multiple_rel_eff %} Rel. Eff. {{ curve_data.curve_index }}{% endif %}</caption>
                <thead>
                    <tr>
                        <th>Nuclides</th>
                        <th>Mass Ratio</th>
                        <th>Activity Ratio</th>
                        <th>Uncertainty</th>
                    </tr>
                </thead>
                <tbody>
                    {% for pair in curve_data.pairs %}
                        <tr>
                            <td>{{ pair.numerator }}/{{ pair.denominator }}</td>
                            <td>{% if pair.mass_ratio %}{{ pair.mass_ratio | scientific(6) }}{% else %}--{% endif %}</td>
                            <td>{{ pair.activity_ratio | scientific(6) }}</td>
                            <td>{% if pair.activity_ratio_uncertainty_percent %}{{ "%.3f" | format(pair.activity_ratio_uncertainty_percent) }}%{% else %}--{% endif %}</td>
                        </tr>
                        <tr>
                            <td>{{ pair.denominator }}/{{ pair.numerator }}</td>
                            <td>{% if pair.mass_ratio %}{{ (1.0 / pair.mass_ratio) | scientific(6) }}{% else %}--{% endif %}</td>
                            <td>{{ (1.0 / pair.activity_ratio) | scientific(6) }}</td>
                            <td>{% if pair.activity_ratio_uncertainty_percent %}{{ "%.3f" | format(pair.activity_ratio_uncertainty_percent) }}%{% else %}--{% endif %}</td>
                        </tr>
                    {% endfor %}
                </tbody>
            </table>
        {% endfor %}
    {% endif %}
    
    {% if energy_calibration.was_fit and energy_calibration.adjustments %}
        <h2>Energy Calibration Adjustments</h2>
        <div>
            Fit 
            {% for adj in energy_calibration.adjustments %}
                {% if not loop.first %} and {% endif %}
                {{ adj.type }} adjustment of {{ "%.3f" | format(adj.physical_value) }} {{ adj.units }}
            {% endfor %}
        </div>
    {% endif %}
    
    {% if warnings %}
        <h2>Warnings</h2>
        {% for warning in warnings %}
            <div class="warning">{{ warning }}</div>
        {% endfor %}
    {% endif %}
    
    {% if timing %}
        <div class="anacomputetime">
            Computation took {{ timing.duration_formatted }} with {{ timing.function_eval_solution }} function calls to solve, 
            and {{ timing.function_eval_total - timing.function_eval_solution }} more for covariance
        </div>
    {% endif %}
    
    <div class="anatime">
        Analysis performed {{ timestamps.local }}{% if timestamps.utc %} ({{ timestamps.utc }} UTC){% endif %}
    </div>
</body>
</html>)";
}

std::string get_text_template()
{
  return R"(
{{ spectrum_title | default("Relative Activity Analysis Report") }}
================================================================

{% if status.success -%}
ANALYSIS RESULTS:
χ² = {{ "%.3f" | format(chi2) }}, DOF = {{ dof }}, χ²/DOF = {{ "%.3f" | format(chi2_per_dof) }}

{% else -%}
ANALYSIS FAILED: {{ status.error_message }}

{% endif -%}

{% if rel_eff_curves -%}
RELATIVE EFFICIENCY EQUATIONS:
{% for curve in rel_eff_curves -%}
Rel. Eff. Eqn{% if have_multiple_rel_eff %} {{ curve.index }}{% endif %}: y = {{ curve.equation_text }}
{% endfor -%}

{% endif -%}

{% if relative_activities -%}
{% for curve_data in relative_activities -%}
RELATIVE ACTIVITIES AND MASS FRACTIONS{% if have_multiple_rel_eff %} (Rel. Eff. {{ curve_data.curve_index }}){% endif %}:
{% for nuclide in curve_data.nuclides -%}
{{ nuclide.name }}: 
  Rel. Activity: {{ nuclide.rel_activity | scientific(6) }} ± {{ nuclide.rel_activity_uncertainty | scientific(6) }}
  {% if nuclide.enrichment -%}Enrichment: {{ "%.2f" | format(nuclide.enrichment * 100) }}%{% endif %}
  {% if nuclide.detector_counts -%}Detector Counts: {{ nuclide.detector_counts | scientific(6) }}{% endif %}
{% endfor -%}

{% endfor -%}
{% endif -%}

{% if energy_calibration.was_fit and energy_calibration.adjustments -%}
ENERGY CALIBRATION ADJUSTMENTS:
{% for adj in energy_calibration.adjustments -%}
{{ adj.type | title }} adjustment: {{ "%.3f" | format(adj.physical_value) }} {{ adj.units }}
{% endfor -%}

{% endif -%}

{% if warnings -%}
WARNINGS:
{% for warning in warnings -%}
- {{ warning }}
{% endfor -%}

{% endif -%}

{% if timing -%}
TIMING:
Computation took {{ timing.duration_formatted }} with {{ timing.function_eval_solution }} function calls to solve
{% endif -%}

Analysis performed: {{ timestamps.local }}{% if timestamps.utc %} ({{ timestamps.utc }} UTC){% endif %}
)";
}

std::string load_template_file(const std::string& template_path)
{
  std::ifstream file(template_path);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open template file: " + template_path);
  }
  
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

std::string format_number(double value, int precision)
{
  std::stringstream ss;
  ss << std::scientific << std::setprecision(precision) << value;
  return ss.str();
}

std::string format_percentage(double value, int precision)
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << (value * 100.0) << "%";
  return ss.str();
}

std::string format_duration(int microseconds)
{
  const double seconds = microseconds * 1.0e-6;
  return PhysicalUnits::printToBestTimeUnits(seconds);
}

namespace InjaCallbacks
{
  inja::Environment create_environment()
  {
    inja::Environment env;
    
    // Add custom functions for scientific notation
    env.add_callback("scientific", 2, [](inja::Arguments& args) {
      double value = args.at(0)->get<double>();
      int precision = args.at(1)->get<int>();
      return format_number(value, precision);
    });
    
    // Add function for formatting with printf-style
    env.add_callback("format", 2, [](inja::Arguments& args) {
      std::string format_str = args.at(0)->get<std::string>();
      double value = args.at(1)->get<double>();
      
      char buffer[256];
      std::snprintf(buffer, sizeof(buffer), format_str.c_str(), value);
      return std::string(buffer);
    });
    
    // Add default filter
    env.add_callback("default", 2, [](inja::Arguments& args) {
      if (args.at(0)->is_null() || args.at(0)->empty()) {
        return args.at(1)->get<std::string>();
      }
      return args.at(0)->get<std::string>();
    });
    
    return env;
  }
}

} // namespace RelActAutoReport 
