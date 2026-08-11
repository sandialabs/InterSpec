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


#include <set>
#include <map>
#include <memory>
#include <vector>


#include <Wt/Utils>
#include <Wt/WColor>
#include <Wt/WPoint>
#include <Wt/WLength>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WCssStyleSheet>
#include <Wt/WContainerWidget>


#include "SandiaDecay/SandiaDecay.h"


#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/RelEffChart.h"
#include "InterSpec/RelActCalcManual.h"


using namespace Wt;
using namespace std;

namespace
{
/** The source assigned to a peak, as a `SrcVariant`; a monostate variant if the peak has no source. */
RelActCalcAuto::SrcVariant peak_src_variant( const PeakDef &p )
{
  if( const SandiaDecay::Nuclide *nuc = p.parentNuclide() )
    return nuc;
  if( const SandiaDecay::Element *el = p.xrayElement() )
    return el;
  if( const ReactionGamma::Reaction *rctn = p.reaction() )
    return rctn;
  return RelActCalcAuto::SrcVariant{};
}//peak_src_variant(...)


/** User-facing text for why a rel. eff. data point was left off the chart.

 Kept in sync with `RelActCalcAuto::RelActAutoSolution::ObsEff::ExclusionReason` - see the comment on that enum.
 Falls back to English when there is no Wt session (e.g., a batch-mode HTML report).
 */
WString exclusion_reason_text( const RelActCalcAuto::RelActAutoSolution::ObsEff::ExclusionReason reason )
{
  using ExclusionReason = RelActCalcAuto::RelActAutoSolution::ObsEff::ExclusionReason;

  const char *key = "", *fallback = "";

  switch( reason )
  {
    case ExclusionReason::NotExcluded:
      return WString();

    case ExclusionReason::NoCountsForCurve:
      key = "rec-excl-no-counts";
      fallback = "This relative efficiency curve has no gammas here.";
      break;

    case ExclusionReason::NonPositiveEff:
      key = "rec-excl-non-pos-eff";
      fallback = "The freely-fit peak area came out negative, so there is no efficiency to show.";
      break;

    case ExclusionReason::Insignificant:
      key = "rec-excl-insignificant";
      fallback = "Not statistically significant compared to its uncertainty.";
      break;

    case ExclusionReason::TailLeakage:
      key = "rec-excl-tail-leakage";
      fallback = "Most of the counts here come from neighboring peaks' tails, so a small peak-shape"
                 " mismatch would dominate the value.";
      break;

    case ExclusionReason::OutsideRoi:
      key = "rec-excl-outside-roi";
      fallback = "The peak is not fully contained within its ROI, so its area is not fully measured.";
      break;
  }//switch( reason )

  return wApp ? WString::tr(key) : WString::fromUTF8(fallback);
}//exclusion_reason_text(...)
}//namespace


RelEffChart::RelEffChart( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_jsgraph( jsRef() + ".chart" ),
  m_xAxisTitle{},
  m_yAxisTitle{},
  m_topMargin( 5 ),
  m_rightMargin( 5 ),
  m_bottomMargin( 0 ),
  m_leftMargin( 5 ),
  m_titlePadding( -3 )
{
  addStyleClass( "RelEffChart" );
  setOverflow(Overflow::OverflowHidden);
  
  // We will actually totally define the CSS through WApplication, rather than using the same CSS
  //  file we use for the self-contained HTML reports.
  //wApp->useStyleSheet( "InterSpec_resources/RelEffPlot.css" );
  
  setCssRules();
  InterSpec *interspec = InterSpec::instance();
  if( interspec )
  {
    interspec->useMessageResourceBundle( "RelEffChart" );
    interspec->colorThemeChanged().connect( this, &RelEffChart::setCssRules );
  }

  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/RelEffPlot.js" );
}//RelEffChart constructor


void RelEffChart::defineJavaScript()
{
  string options = "{ "
  "margins: {"
  " top: " + std::to_string(m_topMargin) + ","
  " right: " + std::to_string(m_rightMargin) + ","
  " bottom: " + std::to_string(m_bottomMargin) + ","
  " left: " + std::to_string(m_leftMargin) +
  " }, "
  " titleToAxisPadding: " + std::to_string(m_titlePadding);
  
  if( !m_xAxisTitle.empty() )
    options += ", xAxisTitle: \"" + m_xAxisTitle + "\"";
  if( !m_yAxisTitle.empty() )
    options += ", yAxisTitle: \"" + m_yAxisTitle + "\"";
  options += " }";
  
  setJavaScriptMember( "chart", "new RelEffPlot(" + jsRef() + ", " + options + ");");
  
  // ResizeObserver is Safari 13.1+; on a 13.0 WKWebView an unguarded `new ResizeObserver` throws,
  //  killing the rest of this script block - including the pending setRelEffData flush below - so
  //  the chart would render empty with dead interactions.  Degrading to no resize handling there
  //  is the acceptable fallback.
  setJavaScriptMember( "resizeObserver",
    "(typeof ResizeObserver === 'undefined') ? null : new ResizeObserver(entries => {"
      "for (let entry of entries) {"
        "if( entry.target && (entry.target.id === '" + id() + "') ){"
          // When we "Clear Session", jsRef() will give a null result temporarily, so we'll protect against that
          "const c=" + jsRef() + ";"
          "if(c && c.chart)"
            "c.chart.handleResize();"
        "}"
      "}"
    "});"
  );

  doJavaScript( "{const c=" + jsRef() + "; if(c && c.resizeObserver) c.resizeObserver.observe(c);}" );
  
  // Set dataset colors after chart is initialized
  setRelEffCurveColors();
  
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void RelEffChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}


void RelEffChart::setData(const std::vector<RelEffChart::ReCurveInfo> &infoSets)
{
  std::string jsonData = jsonForData(infoSets);
  sendDataToJavaScript(jsonData);
}

// For backward compatibility
void RelEffChart::setData(const RelEffChart::ReCurveInfo &info)
{
  std::vector<RelEffChart::ReCurveInfo> infoSets = {info};
  setData(infoSets);
}

std::string RelEffChart::jsonForData(const std::vector<ReCurveInfo> &infoSets)
{
  std::vector<RelEffChartDataset> datasets;
  for(const ReCurveInfo &info : infoSets)
  {
    const double live_time = info.live_time;
    const vector<RelActCalcAuto::RelActAutoSolution::ObsEff> &obs_eff_data = info.obs_eff_data;
    const std::vector<RelActCalcAuto::NuclideRelAct> &rel_acts = info.rel_acts;
    const std::string &relEffEqn = info.js_rel_eff_eqn;
    const std::string &relEffUncertEqn = info.js_rel_eff_uncert_eqn;
    const WString &re_name = info.re_curve_name;
    const WString &chi2_title_str = info.re_curve_eqn_txt;
    
    RelEffChartDataset dataset;
    dataset.relEffEqn = relEffEqn;
    dataset.curveName = re_name;
    dataset.chi2_title_str = chi2_title_str;
    dataset.relEffEqnUncert = relEffUncertEqn.empty() ? "null" : relEffUncertEqn;
    
    std::vector<RelActCalcManual::GenericPeakInfo> peaks;
    map<string,pair<double,string>> relActsColors;
    
    for(RelActCalcAuto::RelActAutoSolution::ObsEff obsEff : obs_eff_data)
    {
      // Decide if this point is trustworthy enough to plot; if not, record it (and why) so the chart can
      //  show the user what was left off, rather than silently dropping it.
      if( !RelActCalcAuto::RelActAutoSolution::show_obs_eff_point( obsEff ) )
      {
        RelEffChart::OmittedPoint omitted;
        omitted.energy = obsEff.energy;
        omitted.reason = exclusion_reason_text( obsEff.exclusion_reason );
        omitted.fit_counts = obsEff.curve_fit_amplitude;
        omitted.expected_counts = obsEff.curve_model_fraction * obsEff.initial_clustered_peak_amplitude;

        string srcs;
        set<string> seen_srcs;
        for( const PeakDef &p : obsEff.fit_peaks ) //already sorted largest peak first
        {
          const string name = RelActCalcAuto::to_name( peak_src_variant(p) );
          if( !name.empty() && seen_srcs.insert(name).second )
            srcs += (srcs.empty() ? "" : ", ") + name;
        }
        omitted.sources = WString::fromUTF8( srcs );

        dataset.omittedPoints.push_back( std::move(omitted) );
        continue;
      }//if( point is not shown )

      // Create a GenericPeakInfo from the ObsEff data
      RelActCalcManual::GenericPeakInfo peak;
      peak.m_energy = obsEff.energy;
      peak.m_mean = obsEff.energy;
      peak.m_fwhm = 2.35482 * obsEff.effective_sigma; // Convert sigma to FWHM
      peak.m_counts = 0.0;

      // The fractional uncertainty of this curves share of the cluster; includes both the free-fit area
      //  uncertainty, and (when more than one rel. eff. curve has gammas here) the unmeasurable split
      //  between curves - see `ObsEff::curve_model_fraction`.
      const double frac_uncert = (obsEff.observed_efficiency > 0.0)
                                   ? (obsEff.observed_efficiency_uncert / obsEff.observed_efficiency)
                                   : 0.0;

      // Process all peaks in this ObsEff to create GenericLineInfo entries
      for(const PeakDef &p : obsEff.fit_peaks)
      {
        const RelActCalcAuto::SrcVariant peak_src = peak_src_variant( p );
        if( RelActCalcAuto::to_name(peak_src).empty() )
          continue; // Skip free-floating peaks, although I dont think they should be here

        // Find corresponding NuclideRelAct
        const RelActCalcAuto::NuclideRelAct *nuc_info = nullptr;
        for( size_t i = 0; !nuc_info && (i < rel_acts.size()); ++i )
        {
          if( rel_acts[i].source == peak_src )
            nuc_info = &(rel_acts[i]);
        }
        
        assert( nuc_info );
        if( !nuc_info )
          continue;
        
        RelActCalcManual::GenericLineInfo line;
        line.m_isotope = RelActCalcAuto::to_name(peak_src);
        
        // Calculate m_yield = amplitude / (observed_efficiency * nuclide_activity)
        if( (obsEff.observed_efficiency > 0.0) && (nuc_info->rel_activity > 0.0) )
        {
          line.m_yield = p.amplitude() / (obsEff.observed_efficiency * nuc_info->rel_activity);
        }else
        {
          line.m_yield = 0.0;
          cerr << "Warning: Invalid observed_efficiency (" << obsEff.observed_efficiency
               << ") or rel_activity (" << nuc_info->rel_activity
               << ") for " << line.m_isotope << " at " << obsEff.energy << " keV" << endl;
        }
        peak.m_source_gammas.push_back(line);

        // Set up the source -> (activity, color) mapping; all peaks of a source get the color of the primary
        //  peak.  Every source must get an entry, even without a color: `jsonForDataset(...)` sums the source
        //  activities from this map to get the denominator of the plotted efficiency, so leaving a source out
        //  would drop its counts from the denominator but not the numerator, and the point would plot too high.
        //  `p.lineColor()` is set by RelActCalcAuto from `NucInputInfo::peak_color_css`; when that was not
        //  populated the color is left empty here, and the chart falls back to the curves own color.
        const string src_name = RelActCalcAuto::to_name(peak_src);
        const auto pos = relActsColors.find(src_name);
        if( pos == std::end(relActsColors) )
        {
          const string css_color = p.lineColor().isDefault() ? string() : p.lineColor().cssText();
          relActsColors[src_name] = std::make_pair(nuc_info->rel_activity, css_color);
        }

        peak.m_counts += p.amplitude();
      }//for(const PeakDef &p : obsEff.fit_peaks)

      peak.m_counts_uncert = frac_uncert * peak.m_counts;

      peaks.push_back( std::move(peak) );
      dataset.peakCurveFractions.push_back( obsEff.curve_model_fraction );
    }//for(const RelActAutoSolution::ObsEff &obsEff : obs_eff_data)

    assert( peaks.size() == dataset.peakCurveFractions.size() );

    dataset.peaks = peaks;
    dataset.relActsColors = relActsColors;
    
    datasets.push_back(dataset);
  }
  
  return jsonForData(datasets);
}

// Implementation of the setData method that handles multiple datasets
void RelEffChart::setData(const std::vector<RelEffChartDataset> &datasets)
{
  std::string jsonData = jsonForData(datasets);
  sendDataToJavaScript(jsonData);
  
  // Previous to 20250520, we used to use CSS to color data markers.
  //  Then we switched to just doing it in the JSON/JS directly, to make it
  //  easier to keep things in sync acros doing HTML reports and within interactive
  //  InterSpec - leaving the code in, but commented out for for the moment incase we
  //  want to go back
  /*
  assert(wApp);
  WCssStyleSheet &style = wApp->styleSheet();
  
  set<string> nucs_with_colors;
  for( size_t i = 0; i < datasets.size(); ++i )
  {
    const auto &dataset = datasets[i];
    
    for(const auto &nuc_act_color : dataset.relActsColors)
    {
      string nuc = nuc_act_color.first;
      
      
      const string &color = nuc_act_color.second.second;
      
      if(nuc.empty() || color.empty() || nucs_with_colors.count(nuc))
        continue;
      
      const string rulename = "#" + id() + " .RelEffPlot circle." + nuc;
      if(m_cssRules.count(rulename))
        style.removeRule(m_cssRules[rulename]);
      
      m_cssRules[rulename] = style.addRule(rulename, "fill: " + color + ";");
      nucs_with_colors.insert(nuc);
    }
  }
  */
}

std::string RelEffChart::jsonForData(const std::vector<RelEffChartDataset> &datasets)
{
  // Build the JSON array of dataset objects for JavaScript
  stringstream datasetsJson;
  datasetsJson << "[";
  
  for(size_t datasetIndex = 0; datasetIndex < datasets.size(); ++datasetIndex) {
    const RelEffChartDataset &dataset = datasets[datasetIndex];
    
    if(datasetIndex > 0) {
      datasetsJson << ", ";
    }
    
    datasetsJson << jsonForDataset(dataset, datasetIndex == 0);
  }
  
  datasetsJson << "]";
  
  return datasetsJson.str();
}

std::string RelEffChart::jsonForDataset(const RelEffChartDataset &dataset, bool isFirstDataset)
{
  char buffer[512] = { '\0' };

  // printf's %g/%f render NaN/Inf as bare `nan`/`inf` tokens - not valid JavaScript - and a single
  //  such token makes the whole `...setRelEffData(...)` statement a SyntaxError, silently dropping
  //  the entire chart update (stale/empty chart, and on a first render a dead drag-zoom because
  //  the mouse-capture rect never gets sized).  So every numeric field goes through this.
  const auto finite_or = []( const double val, const double fallback ) -> double {
    return (IsNan(val) || IsInf(val)) ? fallback : val;
  };

  // Start dataset object
  stringstream datasetJson;
  datasetJson << "{";
  
  // Add the data values JSON
  stringstream rel_eff_plot_values;
  size_t njson_entries = 0;
  rel_eff_plot_values << "\"data_vals\": [";
  for(size_t index = 0; index < dataset.peaks.size(); ++index)
  {
    const RelActCalcManual::GenericPeakInfo &peak = dataset.peaks[index];
    
    string isotopes_json;
    double src_counts = 0.0;
    for(const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas)
    {
      const auto pos = dataset.relActsColors.find(line.m_isotope);
      if(pos == end(dataset.relActsColors))
      {
        snprintf(buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G}",
                (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(),
                finite_or(line.m_yield, 0.0));
      }else
      {
        const double rel_act = finite_or( pos->second.first, 0.0 );
        const string &color = pos->second.second;
        src_counts += rel_act * finite_or( line.m_yield, 0.0 );
        if( color.empty() ) //no color assigned to this source; the chart will use the curves own color
          snprintf(buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G, \"rel_act\": %1.6G}",
                   (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(),
                   finite_or(line.m_yield, 0.0), rel_act );
        else
          snprintf(buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G, \"rel_act\": %1.6G, \"color\": \"%s\"}",
                   (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(),
                   finite_or(line.m_yield, 0.0), rel_act, color.c_str() );
      }

      isotopes_json += buffer;
    }//for(const RelEff::GammaLineInfo &line : peak.m_source_gammas)

    const double eff = peak.m_counts / src_counts;
    double eff_uncert = peak.m_counts_uncert / src_counts;
    
    if(IsNan(eff) || IsInf(eff))
    {
      cerr << "RelEffChart::jsonForDataset: Got invalid eff at " << peak.m_energy << ", peak.m_counts="
          << peak.m_counts << ", src_counts=" << src_counts << endl;
      
      continue;
    }
    
    if(IsNan(eff_uncert) || IsInf(eff_uncert))
    {
      cerr << "RelEffChart::jsonForDataset: Got invalid eff_uncert at " << peak.m_energy
          << ", peak.m_counts=" << peak.m_counts_uncert << ", src_counts=" << src_counts << endl;
      
      eff_uncert = 0.0;
    }
    
    // How much of this energy's counts this curve is assigned; < 1 means another rel. eff. curve also has
    //  gammas here, and the split between them is fitted through the attenuation physics / curve shapes
    //  (constrained by the rest of the spectrum) rather than resolved by this peak itself.
    const double blend_frac = (index < dataset.peakCurveFractions.size())
                                ? finite_or( dataset.peakCurveFractions[index], 1.0 ) : 1.0;

    snprintf(buffer, sizeof(buffer),
            "%s{\"energy\": %.2f, \"mean\": %.2f, \"counts\": %1.7g, \"counts_uncert\": %1.7g,"
            " \"eff\": %1.6g, \"eff_uncert\": %1.6g, \"blend_frac\": %1.4g, \"nuc_info\": ",
            (njson_entries ? ", " : ""), finite_or(peak.m_energy, 0.0), finite_or(peak.m_mean, 0.0),
            finite_or(peak.m_counts, 0.0), finite_or(peak.m_counts_uncert, 0.0),
            eff, eff_uncert, blend_frac);

    rel_eff_plot_values << buffer;
    rel_eff_plot_values << "[" << isotopes_json.c_str() << "]}";

    njson_entries += 1;
  }//for(size_t index = 0; index < input_peaks.size(); ++index)
  
  rel_eff_plot_values << "]";
  
  // Add data_vals to the dataset
  if(njson_entries < 1)
    datasetJson << "\"data_vals\": null";
  else
    datasetJson << rel_eff_plot_values.str();
  
  // Add fit_eqn to the dataset
  datasetJson << ", \"fit_eqn\": " << (dataset.relEffEqn.empty() ? "null" : dataset.relEffEqn);
  
  // Add chi2_txt to the dataset
  datasetJson << ", \"chi2_txt\": " << (dataset.chi2_title_str.empty() ? "null" : dataset.chi2_title_str.jsStringLiteral());

  // The curves name; labels its section of the omitted-points panel when there is more than one curve
  datasetJson << ", \"curve_name\": " << (dataset.curveName.empty() ? "null" : dataset.curveName.jsStringLiteral());
  
  // Add fit_uncert_fcn to the dataset
  datasetJson << ", \"fit_uncert_fcn\": " << (dataset.relEffEqnUncert.empty() ? "null" : dataset.relEffEqnUncert);

  // The points intentionally left off the chart, and why; the chart shows these behind an info icon, so the
  //  user can check that what was omitted is believable, rather than the points just silently disappearing.
  if( dataset.omittedPoints.empty() )
  {
    datasetJson << ", \"omitted_pts\": null";
  }else
  {
    datasetJson << ", \"omitted_pts\": [";
    for( size_t i = 0; i < dataset.omittedPoints.size(); ++i )
    {
      const OmittedPoint &pt = dataset.omittedPoints[i];
      snprintf( buffer, sizeof(buffer),
               "%s{\"energy\": %.2f, \"fit_counts\": %1.7g, \"expected_counts\": %1.7g, \"srcs\": ",
               (i ? ", " : ""), finite_or(pt.energy, 0.0), finite_or(pt.fit_counts, 0.0),
               finite_or(pt.expected_counts, 0.0) );
      datasetJson << buffer << pt.sources.jsStringLiteral()
                  << ", \"reason\": " << pt.reason.jsStringLiteral() << "}";
    }
    datasetJson << "]";
  }//if( no omitted points ) / else

  // Close dataset object
  datasetJson << "}";
  
  return datasetJson.str();
}

void RelEffChart::sendDataToJavaScript(const std::string &jsonData)
{
  const string js = m_jsgraph + ".setRelEffData(" + jsonData + ");";
  
  if(isRendered())
    doJavaScript(js);
  else
    m_pendingJs.push_back(js);
}

// Backwards compatibility method
void RelEffChart::setData(const std::vector<RelActCalcManual::GenericPeakInfo> &peaks,
                          const map<string,pair<double,string>> &relActsColors,
                          string relEffEqn,
                          const Wt::WString &chi2_title_str,
                          const string &relEffEqnUncert)
{
  RelEffChartDataset dataset;
  dataset.peaks = peaks;
  dataset.relActsColors = relActsColors;
  dataset.relEffEqn = relEffEqn;
  dataset.chi2_title_str = chi2_title_str;
  dataset.relEffEqnUncert = relEffEqnUncert;
  
  std::vector<RelEffChartDataset> datasets = {dataset};
  setData(datasets);
}

void RelEffChart::setCssRules()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  std::shared_ptr<const ColorTheme> theme = interspec ? interspec->getColorTheme() : nullptr;

  assert( theme );
  // NOTE: a null theme must NOT skip the structural rules below (the app does not load
  //  RelEffPlot.css) - without e.g. the tooltip's `position: fixed; pointer-events: none;` the
  //  tooltip lays out in-flow inside an overflow:hidden div and hover appears broken until a
  //  colorThemeChanged() re-runs this.  Only the theme-color calls at the end need the theme.

  WCssStyleSheet &style = wApp->styleSheet();
  
  //m_cssRules[".RelEffPlot"] = style.addRule( ".RelEffPlot", "" );
  string rulename = ".RelEffPlot circle - cursor";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlot circle", "cursor: pointer;" );
  
  rulename = "div.RelEffPlotTooltip";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " div.RelEffPlotTooltip",
                                         "position: fixed;"
                                         " padding: 6px;"
                                         " font: 12px sans-serif;"
                                         " background: #ffffcc;"
                                         " border: 0px;"
                                         " border-radius: 8px;"
                                         " pointer-events: none;"
                                         " color: #444422;" );
  
  // Points whose counts are shared with another rel. eff. curve are drawn hollow - the split between curves
  //  is fitted through the attenuation physics / curve shapes rather than resolved by the peak itself, so
  //  such a point is not independent evidence for this curve.  See `ObsEff::curve_model_fraction`.
  rulename = ".RelEffPlot circle.blended";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlot circle.blended",
                                         "fill-opacity: 0.15; stroke-width: 1.5;" );

  // The "some points were left off this chart" affordance; see `RelEffPlot.prototype.updateOmittedPointsInfo`.
  //  Note the icon/panel are positioned against this widgets div (which has class "RelEffChart"); the
  //  "RelEffPlot" class is on the <svg> inside it, so it must not be used as the positioning ancestor.
  rulename = "RelEffChart - position";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id(), "position: relative;" );

  rulename = ".RelEffPlotOmittedIcon";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedIcon",
                                         "position: absolute; left: 2px; top: 2px;"
                                         " width: 15px; height: 15px; line-height: 15px;"
                                         " text-align: center; font: italic bold 11px serif;"
                                         " border-radius: 50%; border: 1px solid currentColor;"
                                         " opacity: 0.55; cursor: pointer; user-select: none;" );

  rulename = ".RelEffPlotOmittedIcon:hover";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedIcon:hover", "opacity: 1;" );

  rulename = ".RelEffPlotOmittedPanel - layout";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedPanel",
                                         "position: absolute; left: 17px; top: 2px;" //adjoins the icon, so
                                         //  there is no gap for the pointer to fall through on its way over
                                         " max-width: 95%; max-height: 90%; overflow: auto; z-index: 10;"
                                         " padding: 6px 8px; font: 11px sans-serif; border-radius: 6px;"
                                         " border-style: solid; border-width: 1px;"
                                         " box-shadow: 0 1px 4px rgba(0,0,0,0.25);" );

  rulename = ".RelEffPlotOmittedPanel td";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedPanel td",
                                         "padding: 1px 6px 1px 0; vertical-align: top;" );

  rulename = ".RelEffPlotOmittedHeader";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedHeader td",
                                         "font-weight: bold; border-bottom: 1px solid currentColor;" );

  rulename = ".RelEffPlotOmittedTitle";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedTitle",
                                         "font-weight: bold; margin-bottom: 2px;" );

  rulename = ".RelEffPlotOmittedExplain";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedExplain",
                                         "margin-bottom: 4px; opacity: 0.8;" );

  rulename = ".RelEffPlotOmittedCurve";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedCurve",
                                         "font-weight: bold; padding-top: 4px;" );

  rulename = ".RelEffPlotBlendNote";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotBlendNote",
                                         "margin-top: 3px; max-width: 260px; opacity: 0.85;" );

  // The panel and icon sit over the chart, so they need colors of their own.  Use the apps global color
  //  theme variables (set on :root by `InterSpec::applyColorTheme`), so they follow the theme like the rest
  //  of the non-chart UI; the fallbacks are for the self-contained HTML reports, which have no variables.
  rulename = ".RelEffPlotOmittedPanel - colors";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlotOmittedPanel, "
                                          "#" + id() + " .RelEffPlotOmittedIcon",
                                         "background: var(--interspec-background-color, #ffffee);"
                                         " color: var(--interspec-text-color, #333322);"
                                         " border-color: var(--interspec-border-color, #cccc99);" );

  if( theme )
  {
    setLineColor( theme->foregroundLine );
    setDefaultMarkerColor( theme->backgroundLine );
    setChartBackgroundColor( theme->spectrumChartBackground );
    setAxisLineColor( theme->spectrumAxisLines );
    setTextColor( theme->spectrumChartText );

    // Update dataset colors when theme changes
    setRelEffCurveColors();
  }//if( theme )
}//void setCssRules()


void RelEffChart::setLineColor( const Wt::WColor &color )
{
  WColor lineColor = color.isDefault() ? WColor("steelblue") : color;
  
  WCssStyleSheet &style = wApp->styleSheet();
  string rulename = ".RelEffPlot path";
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlot path",
                                       "stroke: " + lineColor.cssText() + ";"
                                       " stroke-width: 2; fill: none;" );
  
  
  rulename = ".RelEffPlot .errorbar";
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( "#" + id() + " .RelEffPlot .errorbar",
                                       "fill: none;"
                                       " stroke-width: 1;"
                                       " stroke: " + lineColor.cssText() + ";" );
}//setLineColor(...)


void RelEffChart::setTextColor( const Wt::WColor &color )
{
  WColor textColor = color.isDefault() ? WColor(0,0,0) : color;
  const string c = textColor.cssText();
  
  const string rulename = "TextColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  const string div_sel = "#" + id() + " ";
  m_cssRules[rulename] = style.addRule( div_sel + ".RelEffPlot .xaxistitle, "
                                        + div_sel + ".RelEffPlot .yaxistitle, "
                                        + div_sel + ".RelEffPlot .yaxis, .RelEffPlot .yaxislabel, "
                                        + div_sel + ".RelEffPlot .xaxis, "
                                        + div_sel + ".RelEffPlot .tick text", "fill: " + c );
}


void RelEffChart::setAxisLineColor( const Wt::WColor &color )
{
  WColor axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  string rulename = "AxisColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  const string div_sel = "#" + id() + " ";
  const string selector = div_sel + ".RelEffPlot .xAxis path, "
                        + div_sel + ".RelEffPlot .xAxis line, "
                        + div_sel + ".RelEffPlot .yAxis path, "
                        + div_sel + ".RelEffPlot .yAxis line";
  const string css = "fill: none;"
                     " stroke-width: 1;"
                     " shape-rendering: crispEdges;"
                     " stroke: " + axisColor.cssText() + ";";
  m_cssRules[rulename] = style.addRule( selector, css );
}


void RelEffChart::setChartBackgroundColor( const Wt::WColor &color )
{
  const string c = color.isDefault() ? "rgba(0,0,0,0)" : color.cssText();

  const string rulename = "BackgroundColor";

  WCssStyleSheet &style = wApp->styleSheet();

  if( color.isDefault() )
  {
    if( m_cssRules.count(rulename) )
      style.removeRule( m_cssRules[rulename] );
    
    return;
  }//if( color.isDefault() )
  
  m_cssRules[rulename] = style.addRule( "#" + id() + " > svg", "background: " + c + ";" );
}


void RelEffChart::setDefaultMarkerColor( const Wt::WColor &color )
{
  WColor markerColor = color.isDefault() ? WColor(0,51,255,155) : color;
  
  WCssStyleSheet &style = wApp->styleSheet();
  const string rulename = ".RelEffPlot circle - def-color";
  
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  const string div_sel = "#" + id() + " ";
  const string selector = div_sel + ".RelEffPlot circle, "
                        + div_sel + ".RelEffPlot circle.noiso, "
                        + div_sel + ".RelEffPlot circle.multiiso";
  m_cssRules[rulename] = style.addRule( selector, "fill: " + markerColor.cssText() + ";" );
}

void RelEffChart::setXAxisTitle( const std::string &title )
{
  m_xAxisTitle = Wt::Utils::htmlEncode( title );
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle(\"" + m_xAxisTitle + "\",true);" );
}//


void RelEffChart::setYAxisTitle( const std::string &title )
{
  m_yAxisTitle = Wt::Utils::htmlEncode( title );
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setYAxisTitle(\"" + m_yAxisTitle + "\",true);" );
}


void RelEffChart::setContentMargins( int top, int right, int bottom, int left )
{
  m_topMargin = top;
  m_rightMargin = right;
  m_bottomMargin = bottom;
  m_leftMargin = left;
  
  if( isRendered() )
  {
    const string margins = "{"
    " top: " + std::to_string(m_topMargin) + ","
    " right: " + std::to_string(m_rightMargin) + ","
    " bottom: " + std::to_string(m_bottomMargin) + ","
    " left: " + std::to_string(m_leftMargin) +
    " }";
    doJavaScript( m_jsgraph + ".setMargins(" + margins + ");" );
  }//if( isRendered() )
}//void setContentMargins( int top, int right, int bottom, int left );


RelEffChart::~RelEffChart()
{
  Wt::WApplication *app = Wt::WApplication::instance();
  if( app )
  {
    WCssStyleSheet &style = app->styleSheet();
    for( const auto &rule : m_cssRules )
      style.removeRule( rule.second );
    m_cssRules.clear();
  }//if( app )
}//~RelEffChart()

void RelEffChart::setRelEffCurveColors()
{
  InterSpec *interspec = InterSpec::instance();
  if(!interspec)
    return;
    
  std::shared_ptr<const ColorTheme> theme = interspec->getColorTheme();
  if(!theme)
    return;
    
  // Get colors from theme
  WColor foregroundColor = theme->foregroundLine;
  WColor backgroundColor = theme->backgroundLine;
  WColor secondaryColor = theme->secondaryLine;
  
  // Convert to javascript array of colors
  std::vector<std::string> colors;
  
  // Add colors to array if not default
  if(!foregroundColor.isDefault())
    colors.push_back(foregroundColor.cssText());
      
  if(!secondaryColor.isDefault())
    colors.push_back(secondaryColor.cssText());
    
  if(!backgroundColor.isDefault())
    colors.push_back(backgroundColor.cssText());
  
  if(colors.empty())
    return;
    
  // Build JSON array of colors
  std::stringstream colorsJson;
  colorsJson << "[";
  for(size_t i = 0; i < colors.size(); ++i) {
    if(i > 0)
      colorsJson << ", ";
    colorsJson << "\"" << colors[i] << "\"";
  }
  colorsJson << "]";
  
  // Send to JavaScript
  const string js = m_jsgraph + ".setRelEffCurveColors(" + colorsJson.str() + ");";
  
  if(isRendered())
    doJavaScript(js);
  else
    m_pendingJs.push_back(js);
}//void RelEffChart::setRelEffCurveColors()


