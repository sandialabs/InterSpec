#ifndef FRAMRESULTS_H_
#define FRAMRESULTS_H_

//C++ includes
#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>
//External includes
#ifdef isinf
#undef isinf
#endif
#ifdef isnan
#undef isnan
#endif
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct FRAMResults
{
  struct FRAMHeader
  {
    std::string m_analysis_program;  // "FRAM" or "RelActCalcAuto"
    std::string m_programVersion;
    std::string m_spectrumFile;
    std::time_t m_spectrumDate;
    double m_liveTime;
    double m_realTime;
    double m_numChannels;
    double m_ecalOffset;
    double m_ecalSlope;
    std::time_t m_analysisTime = static_cast<std::time_t>(-1);
    std::string m_parameterSet;
    int m_efficiencyModel;
    std::string m_comment;
  };
  FRAMHeader m_FRAMHeader;

  struct NuclideResult
  {
    std::string m_nuclide;
    double m_massFraction           = 0.0;
    double m_massFractionSigmaStat  = 0.0;
    double m_massFractionSigmaSys   = 0.0;
    double m_massFractionRSDStat    = 0.0;
    double m_massFractionRSDSys     = 0.0;
    double m_totalPower             = 0.0;
  };
  std::vector<NuclideResult> m_nuclideResults;

  struct MassRatios
  {
    std::string m_numeratorNuclide;
    std::string m_denominatorNuclide;
    double m_massRatio        =  0.0;
    double m_massRatioRSDStat = -1.0; // –1 = not applicable
    double m_massRatioRSDSys  = -1.0;
  };
  std::vector<MassRatios> m_massRatios;

  bool   m_diagnosticsPassed;
  double m_effective240Frac;
  double m_effective240SigmaStat;
  double m_effective240RSDStat;
  double m_specificPower           =  0.0;
  double m_specificPowerSigmaStat  =  0.0;
  double m_specificPowerSigmaSys   =  0.0;
  double m_specificPowerRSDStat    =  0.0;
  double m_specificPowerRSDSys     =  0.0;

  double m_timeSinceSep           = -1.0; // –1 = not applicable
  double m_timeSinceSepSigmaStat  = -1.0;
  double m_timeSinceSepRSDStat    = -1.0;

  std::vector<std::string> m_warnings;
}; // struct FRAMResults


/**
 * Helper: parse the date‑time string
 **/
inline std::time_t parse_fram_datetime(const std::string& s)
{
  std::istringstream iss(s);
  std::tm tm = {};
  iss >> std::get_time(&tm, "%d-%b-%Y %H:%M:%S");
  if (iss.fail())
  {
    return static_cast<std::time_t>(-1);
  }
  return std::mktime(&tm);
}


/**
 *   ADL conversion – JSON -> FRAMResults
 **/
inline void from_json(const json& a_json, FRAMResults& a_framResult)
{
  // ---------- Header ----------
  if (a_json.contains("Header")) 
  {
    const json& jHeader = a_json.at("Header");
    FRAMResults::FRAMHeader framHeader;
    if (jHeader.contains("spectrum_file"))
    {
      jHeader.at("spectrum_file").get_to(framHeader.m_spectrumFile);
    }

    if (jHeader.contains("spectrum_date"))
    {
      framHeader.m_spectrumDate = parse_fram_datetime(jHeader.at("spectrum_date").get<std::string>());
    }

    if (jHeader.contains("version"))
    {
      jHeader.at("version").get_to(framHeader.m_programVersion);
    }

    if (jHeader.contains("live_time"))
    {
      jHeader.at("live_time").get_to(framHeader.m_liveTime);
    }

    if (jHeader.contains("real_time"))
    {
      jHeader.at("real_time").get_to(framHeader.m_realTime);
    }

    if (jHeader.contains("num_channels"))
    {
      jHeader.at("num_channels").get_to(framHeader.m_numChannels);
    }

    if (jHeader.contains("ecal_offset"))
    {
      jHeader.at("ecal_offset").get_to(framHeader.m_ecalOffset);
    }

    if (jHeader.contains("ecal_slope"))
    {
      jHeader.at("ecal_slope").get_to(framHeader.m_ecalSlope);
    }

    if (jHeader.contains("analysis_date"))
    {
      framHeader.m_analysisTime = parse_fram_datetime(jHeader.at("analysis_date").get<std::string>());
    }

    if (jHeader.contains("parameter_set"))
    {
      jHeader.at("parameter_set").get_to(framHeader.m_parameterSet);
    }

    if (jHeader.contains("efficiency_model"))
    {
      jHeader.at("efficiency_model").get_to(framHeader.m_efficiencyModel);
    }

    if (jHeader.contains("comment"))
    {
      jHeader.at("comment").get_to(framHeader.m_comment);
    }

    a_framResult.m_FRAMHeader = framHeader;
  }

  // ---------- Diagnostics  ----------
  if (a_json.contains("Diagnostics")) 
  {
    const json& jDiagnostics = a_json.at("Diagnostics");

    if (jDiagnostics.contains("diagnostics_passed"))
    {
      jDiagnostics.at("diagnostics_passed").get_to(a_framResult.m_diagnosticsPassed);
    }
  }

  // ---------- Isotopics ----------
  if (a_json.contains("Isotopics") && a_json.at("Isotopics").is_array()) 
  {
    for (const auto& jIso : a_json.at("Isotopics")) 
    {
      FRAMResults::NuclideResult nuclideResult;
      nuclideResult.m_nuclide = jIso.value("nuclide", "");
      nuclideResult.m_massFraction = jIso.value("mass_percent", 0.0);
      nuclideResult.m_massFractionSigmaStat = jIso.value("mass_sigma_stat", 0.0);
      nuclideResult.m_massFractionSigmaSys = jIso.value("mass_sigma_sys", 0.0);
      nuclideResult.m_massFractionRSDStat = jIso.value("mass_rsd_stat", 0.0);
      nuclideResult.m_massFractionRSDSys = jIso.value("mass_rsd_sys", 0.0);
      nuclideResult.m_totalPower = jIso.value("total_power", 0.0);
      a_framResult.m_nuclideResults.emplace_back(std::move(nuclideResult));
    }
  }

  // ---------- Relative_masses ----------
  if (a_json.contains("Relative_masses") && a_json.at("Relative_masses").is_array()) 
  {
    for (const auto& rel : a_json.at("Relative_masses")) 
    {
      FRAMResults::MassRatios massRatio;
      massRatio.m_numeratorNuclide   = rel.value("numerator",   "");
      massRatio.m_denominatorNuclide = rel.value("denominator", "");
      massRatio.m_massRatio          = rel.value("relative_mass", 0.0);
      massRatio.m_massRatioRSDStat   = rel.value("relative_mass_rsd_stat", 0.0);
      a_framResult.m_massRatios.emplace_back(std::move(massRatio));
    }
  }

  //---------- OtherShortResults ----------
  if (a_json.contains("OtherShortResults")) 
  {
    const json& jOther = a_json.at("OtherShortResults");

    if (jOther.contains("effective_240_frac"))
    {
      jOther.at("effective_240_frac").get_to(a_framResult.m_effective240Frac);
    }

    if (jOther.contains("effective_240_sigma_stat"))
    {
      jOther.at("effective_240_sigma_stat").get_to(a_framResult.m_effective240SigmaStat);
    }

      if (jOther.contains("effective_240_rsd_stat"))
    {
      jOther.at("effective_240_rsd_stat").get_to(a_framResult.m_effective240RSDStat);
    }

    if (jOther.contains("specific_power"))
    {
      jOther.at("specific_power").get_to(a_framResult.m_specificPower);
    }

    if (jOther.contains("specific_power_sigma_stat"))
    {
      jOther.at("specific_power_sigma_stat").get_to(a_framResult.m_specificPowerSigmaStat);
    }

      if (jOther.contains("specific_power_sigma_sys"))
    {
      jOther.at("specific_power_sigma_sys").get_to(a_framResult.m_specificPowerSigmaSys);
    }

    if (jOther.contains("specific_power_rsd_stat"))
    {
      jOther.at("specific_power_rsd_stat").get_to(a_framResult.m_specificPowerRSDStat);
    }

    if (jOther.contains("specific_power_rsd_sys"))
    {
      jOther.at("specific_power_rsd_sys").get_to(a_framResult.m_specificPowerRSDSys);
    }

    if (jOther.contains("time_since_sep"))
    {
      jOther.at("time_since_sep").get_to(a_framResult.m_timeSinceSep);
    }

    if (jOther.contains("time_since_sep_sigma_stat"))
    {
      jOther.at("time_since_sep_sigma_stat").get_to(a_framResult.m_timeSinceSepSigmaStat);
    }

      if (jOther.contains("time_since_sep_rsd_stat"))
    {
      jOther.at("time_since_sep_rsd_stat").get_to(a_framResult.m_timeSinceSepRSDStat);
    }
  }

  // ---------- Warnings ----------
  if (a_json.contains("Warnings") && a_json.at("Warnings").is_array()) 
  {
    for (const auto& jWarn : a_json.at("Warnings")) 
    {
      std::string warning;
      warning = jWarn.value("warning", "");
      a_framResult.m_warnings.emplace_back(warning);
    }
  }

}

// -----------------------------------------------------------------------------
//   ADL conversion – FRAMResults -> JSON
// -----------------------------------------------------------------------------

inline std::string fram_time_to_string(std::time_t t)
{
  if (t == static_cast<std::time_t>(-1))
  {
    return "";
  }
  std::tm* tm = std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(tm, "%d-%b-%Y %H:%M:%S");
  return oss.str();
}


inline void to_json(json& j, const FRAMResults::FRAMHeader& h)
{
  j = json
  {
    {"program",          h.m_analysis_program},
    {"version",          h.m_programVersion},
    {"spectrum_file",    h.m_spectrumFile},
    {"spectrum_date",    fram_time_to_string(h.m_spectrumDate)},
    {"live_time",        h.m_liveTime},
    {"real_time",        h.m_realTime},
    {"num_channels",     h.m_numChannels},
    {"ecal_offset",      h.m_ecalOffset},
    {"ecal_slope",       h.m_ecalSlope},
    {"analysis_date",    fram_time_to_string(h.m_analysisTime)},
    {"parameter_set",    h.m_parameterSet},
    {"efficiency_model", h.m_efficiencyModel},
    {"comment",          h.m_comment}
  };
}


inline void to_json(json& j, const FRAMResults::NuclideResult& n)
{
  j = json
  {
    {"nuclide",          n.m_nuclide},
    {"mass_percent",     n.m_massFraction},
    {"mass_sigma_stat",  n.m_massFractionSigmaStat},
    {"mass_sigma_sys",   n.m_massFractionSigmaSys},
    {"mass_rsd_stat",    n.m_massFractionRSDStat},
    {"mass_rsd_sys",     n.m_massFractionRSDSys},
    {"total_power",      n.m_totalPower}
  };
}


inline void to_json(json& j, const FRAMResults::MassRatios& m)
{
  j = json
  {
    {"numerator",                m.m_numeratorNuclide},
    {"denominator",              m.m_denominatorNuclide},
    {"relative_mass",            m.m_massRatio},
    {"relative_mass_rsd_stat",   m.m_massRatioRSDStat},
    {"relative_mass_rsd_sys",    m.m_massRatioRSDSys}
  };
}


inline void to_json(json& j, const FRAMResults& r)
{
  // Warnings
  json warnings = json::array();

  for (const auto& w : r.m_warnings)
  {
     warnings.push_back( {{"warning", w}} );
  }
  j = json
  {
    {"Header", r.m_FRAMHeader},

    {"Diagnostics", {
        {"diagnostics_passed", r.m_diagnosticsPassed}
    }},

    {"Isotopics", r.m_nuclideResults},

    {"Relative_masses", r.m_massRatios},

    {"OtherShortResults", 
    {
      {"effective_240_frac",             r.m_effective240Frac},
      {"effective_240_sigma_stat",       r.m_effective240SigmaStat},
      {"effective_240_rsd_stat",         r.m_effective240RSDStat},
      {"specific_power",                 r.m_specificPower},
      {"specific_power_sigma_stat",      r.m_specificPowerSigmaStat},
      {"specific_power_sigma_sys",       r.m_specificPowerSigmaSys},
      {"specific_power_rsd_stat",        r.m_specificPowerRSDStat},
      {"specific_power_rsd_sys",         r.m_specificPowerRSDSys},
      {"time_since_sep",                 r.m_timeSinceSep},
      {"time_since_sep_sigma_stat",      r.m_timeSinceSepSigmaStat},
      {"time_since_sep_rsd_stat",        r.m_timeSinceSepRSDStat}
    }},

      {"Warnings", warnings}
  };
}

#endif // FRAMRESULTS_H_
