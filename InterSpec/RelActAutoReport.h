#ifndef RelActAutoReport_h
#define RelActAutoReport_h
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

#include <string>
#include <ostream>

// Include inja for templating
#include "external_libs/SpecUtils/3rdparty/inja/inja.hpp"

// Forward declarations
namespace RelActCalcAuto
{
  struct RelActAutoSolution;
}

namespace RelActAutoReport
{
  /** Enumeration for different report formats */
  enum class ReportFormat : int
  {
    Html,
    Text,
    Json  // Raw JSON data
  };

  /** Configuration for report generation */
  struct ReportConfig
  {
    /** Format of the output report */
    ReportFormat format = ReportFormat::Html;
    
    /** Include spectrum chart in HTML reports */
    bool include_spectrum_chart = true;
    
    /** Include relative efficiency chart in HTML reports */
    bool include_rel_eff_chart = true;
    
    /** Include detailed peak information */
    bool include_peak_details = true;
    
    /** Include computation timing information */
    bool include_timing_info = true;
    
    /** Include warnings in the report */
    bool include_warnings = true;
    
    /** Custom template string (if empty, uses built-in templates) */
    std::string custom_template;
    
    /** Custom CSS for HTML reports */
    std::string custom_css;
    
    /** Custom JavaScript for HTML reports */
    std::string custom_js;
  };

  /** Convert RelActAutoSolution to JSON data suitable for templating
   *
   * @param solution The solution to convert
   * @param config Configuration for what data to include
   * @returns JSON object containing all the templating data
   */
  nlohmann::json solution_to_json(const RelActCalcAuto::RelActAutoSolution& solution, 
                                  const ReportConfig& config = ReportConfig{});

  /** Generate a report using Inja templating
   *
   * @param solution The solution to report on
   * @param config Configuration for the report generation
   * @returns The formatted report as a string
   */
  std::string generate_report(const RelActCalcAuto::RelActAutoSolution& solution,
                             const ReportConfig& config = ReportConfig{});

  /** Generate a report and write it to a stream
   *
   * @param out The output stream to write to
   * @param solution The solution to report on
   * @param config Configuration for the report generation
   */
  void generate_report(std::ostream& out,
                      const RelActCalcAuto::RelActAutoSolution& solution,
                      const ReportConfig& config = ReportConfig{});

  /** Get the built-in template for a specific format
   *
   * @param format The format to get the template for
   * @returns The template string
   */
  std::string get_builtin_template(ReportFormat format);

  /** Get the built-in CSS for HTML reports
   *
   * @returns The CSS string
   */
  std::string get_builtin_css();

  /** Get the built-in JavaScript for HTML reports
   *
   * @returns The JavaScript string
   */
  std::string get_builtin_js();

  /** Load external template file
   *
   * @param template_path Path to the template file
   * @returns The template content
   * @throws std::runtime_error if file cannot be loaded
   */
  std::string load_template_file(const std::string& template_path);

  /** Helper function to format numbers with appropriate precision
   *
   * @param value The number to format
   * @param precision Number of significant digits
   * @returns Formatted string
   */
  std::string format_number(double value, int precision = 6);

  /** Helper function to format percentages
   *
   * @param value The value to format as percentage
   * @param precision Number of decimal places
   * @returns Formatted percentage string
   */
  std::string format_percentage(double value, int precision = 2);

  /** Helper function to format time durations
   *
   * @param microseconds Time in microseconds
   * @returns Human-readable time string
   */
  std::string format_duration(int microseconds);

  /** Custom Inja callbacks for specialized formatting */
  namespace InjaCallbacks
  {
    /** Initialize custom Inja environment with InterSpec-specific functions */
    inja::Environment create_environment();
  }

} // namespace RelActAutoReport

#endif // RelActAutoReport_h 