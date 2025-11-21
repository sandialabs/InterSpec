/* Example usage of RelActAutoReport templating system
 * 
 * This example shows how to use the new Inja-based templating system
 * to generate reports from RelActAutoSolution objects.
 */

#include "InterSpec/RelActAutoReport.h"
#include "InterSpec/RelActCalcAuto.h"

using namespace RelActAutoReport;

void example_basic_usage()
{
  // Assume you have a RelActAutoSolution from your analysis
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Basic HTML report generation
  {
    ReportConfig config;
    config.format = ReportFormat::Html;
    config.include_spectrum_chart = true;
    config.include_rel_eff_chart = true;
    
    std::string html_report = generate_report(solution, config);
    
    // Write to file
    std::ofstream html_file("report.html");
    html_file << html_report;
  }
  
  // Basic text report generation
  {
    ReportConfig config;
    config.format = ReportFormat::Text;
    
    std::string text_report = generate_report(solution, config);
    
    // Write to file or cout
    std::cout << text_report << std::endl;
  }
  
  // JSON data export
  {
    ReportConfig config;
    config.format = ReportFormat::Json;
    
    std::string json_report = generate_report(solution, config);
    
    // You can now use this JSON with external tools
    std::ofstream json_file("data.json");
    json_file << json_report;
  }
}

void example_custom_template()
{
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Using a custom template
  ReportConfig config;
  config.format = ReportFormat::Html;
  
  // Load custom template from file
  try {
    config.custom_template = load_template_file("my_custom_template.inja");
    
    std::string custom_report = generate_report(solution, config);
    
    std::ofstream output("custom_report.html");
    output << custom_report;
    
  } catch (const std::exception& e) {
    std::cerr << "Error loading custom template: " << e.what() << std::endl;
  }
}

void example_external_template_files()
{
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Using external template files
  ReportConfig config;
  
  // For HTML reports with external template
  try {
    config.custom_template = load_template_file("InterSpec_resources/static_text/rel_act_auto_html_report.inja");
    config.format = ReportFormat::Html;
    
    std::string report = generate_report(solution, config);
    
    std::ofstream html_output("external_template_report.html");
    html_output << report;
    
  } catch (const std::exception& e) {
    std::cerr << "Error using external template: " << e.what() << std::endl;
    
    // Fallback to built-in template
    config.custom_template.clear();
    std::string fallback_report = generate_report(solution, config);
    
    std::ofstream fallback_output("fallback_report.html");
    fallback_output << fallback_report;
  }
}

void example_json_data_access()
{
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Get the JSON data structure for custom processing
  ReportConfig config;
  nlohmann::json data = solution_to_json(solution, config);
  
  // Access specific data
  if (data["status"]["success"].get<bool>()) {
    std::cout << "Analysis was successful!" << std::endl;
    std::cout << "Chi-squared: " << data["chi2"].get<double>() << std::endl;
    std::cout << "DOF: " << data["dof"].get<size_t>() << std::endl;
    
    // Access relative activities
    if (data.contains("relative_activities")) {
      for (const auto& curve : data["relative_activities"]) {
        std::cout << "Curve " << curve["curve_index"] << ":" << std::endl;
        for (const auto& nuclide : curve["nuclides"]) {
          std::cout << "  " << nuclide["name"].get<std::string>() 
                   << ": " << nuclide["rel_activity"].get<double>() << std::endl;
        }
      }
    }
  } else {
    std::cout << "Analysis failed: " << data["status"]["error_message"].get<std::string>() << std::endl;
  }
}

void example_stream_output()
{
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Direct output to stream
  ReportConfig config;
  config.format = ReportFormat::Text;
  
  // Output to cout
  std::cout << "=== ANALYSIS REPORT ===" << std::endl;
  generate_report(std::cout, solution, config);
  
  // Output to file stream
  std::ofstream report_file("stream_report.txt");
  generate_report(report_file, solution, config);
}

// Example of integrating with existing RelActAutoSolution methods
void example_replacement_usage()
{
  RelActCalcAuto::RelActAutoSolution solution; // = your_analysis_result;
  
  // Instead of calling solution.print_html_report(stream)
  {
    std::ofstream html_file("new_html_report.html");
    ReportConfig config;
    config.format = ReportFormat::Html;
    generate_report(html_file, solution, config);
  }
  
  // Instead of calling solution.print_summary(stream)
  {
    std::ofstream text_file("new_summary.txt");
    ReportConfig config;
    config.format = ReportFormat::Text;
    config.include_spectrum_chart = false;  // Simplified summary
    config.include_rel_eff_chart = false;
    config.include_peak_details = false;
    generate_report(text_file, solution, config);
  }
} 