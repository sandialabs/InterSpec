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


#include <map>
#include <memory>
#include <vector>


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


RelEffChart::RelEffChart( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_jsgraph( jsRef() + ".chart" )
{
  addStyleClass( "RelEffChart" );
  setOverflow(Overflow::OverflowHidden);
  
  // We will actually totally define the CSS through WApplication, rather than using the same CSS
  //  file we use for the self-contained HTML reports.
  //wApp->useStyleSheet( "InterSpec_resources/RelEffPlot.css" );
  
  setCssRules();
  InterSpec *interspec = InterSpec::instance();
  if( interspec )
    interspec->colorThemeChanged().connect( this, &RelEffChart::setCssRules );
  
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/RelEffPlot.js" );
}//RelEffChart constructor


void RelEffChart::defineJavaScript()
{
  string options = "{ margins: { top: 5, right: 20, bottom: 20, left: 40 } }";
  
  setJavaScriptMember( "chart", "new RelEffPlot(" + jsRef() + ", " + options + ");");
  
  setJavaScriptMember( "resizeObserver",
    "new ResizeObserver(entries => {"
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
  
  callJavaScriptMember( "resizeObserver.observe", jsRef() );
  
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



void RelEffChart::setData( const double live_time,
                          const vector<PeakDef> &fit_peaks,
                          const std::vector<RelActCalcAuto::NuclideRelAct> &rel_acts,
                          const std::string &relEffEqn )
{
  std::vector<RelActCalcManual::GenericPeakInfo> peaks;
  map<string,pair<double,string>> relActsColors;
  
  for( const PeakDef &p : fit_peaks )
  {
    const SandiaDecay::Nuclide *nuc = p.parentNuclide();
  
    // A free-floating peak wont have a nuclide associated
    if( !nuc )
      continue;
    
    RelActCalcManual::GenericPeakInfo peak;
    peak.m_energy = p.mean();
    peak.m_fwhm = p.fwhm();
    peak.m_counts = p.amplitude();
    
    // Amplitude uncertainties arent accurate/relevant, so set to zero.
    peak.m_counts_uncert = 0.0; //p.amplitudeUncert();
    //peak.m_base_rel_eff_uncert = ...;
    
    RelActCalcManual::GenericLineInfo line;
    line.m_isotope = nuc->symbol;
    line.m_yield = 0.0;
    
    const RelActCalcAuto::NuclideRelAct *nuc_info = nullptr;
    for( const auto &rel_act : rel_acts )
      nuc_info = (rel_act.nuclide == nuc) ? &rel_act : nuc_info;
    
    assert( nuc_info );
    if( !nuc_info )
      continue;
    
    for( const pair<double,double> &energy_br : nuc_info->gamma_energy_br )
    {
      if( fabs(energy_br.first - p.gammaParticleEnergy()) < 1.0E-6 )
        line.m_yield += live_time * energy_br.second; 
    }
    
    peak.m_source_gammas.push_back( line );
    
    const auto pos = relActsColors.find(nuc->symbol);
    if( pos == std::end(relActsColors) )
    {
      const string css_color = p.lineColor().isDefault() ? string() : p.lineColor().cssText();
      relActsColors[nuc->symbol] = std::make_pair(nuc_info->rel_activity, css_color);
    }
      
    peaks.push_back( peak );
  }//for( const PeakDef &p : fit_peaks )
  
  setData( peaks, relActsColors, relEffEqn );
}//void setData( std::vector<PeakDef> m_fit_peaks, std::string relEffEqn )


void RelEffChart::setData( const std::vector<RelActCalcManual::GenericPeakInfo> &peaks,
                          const map<string,pair<double,string>> &relActsColors,
                          string relEffEqn )
{
  char buffer[512] = { '\0' };
  
  assert( wApp );
  WCssStyleSheet &style = wApp->styleSheet();
  
  set<string> nucs_with_colors;
  for( const auto nuc_act_color : relActsColors )
  {
    const string &nuc = nuc_act_color.first;
    const string &color = nuc_act_color.second.second;
    
    if( nuc.empty() || color.empty() || nucs_with_colors.count(nuc) )
      continue;
    
    const string rulename = ".RelEffPlot circle." + nuc;
    if( m_cssRules.count(rulename) )
      style.removeRule( m_cssRules[rulename] );
    
    m_cssRules[rulename] = style.addRule( rulename, "fill: " + color + ";" );
    nucs_with_colors.insert( nuc );
  }//for( const auto nuc_act_color : relActsColors )
  
  //Write out the data JSON
  stringstream rel_eff_plot_values, add_rel_eff_plot_css;
  size_t njson_entries = 0;
  rel_eff_plot_values << "[";
  for( size_t index = 0; index < peaks.size(); ++index )
  {
    const RelActCalcManual::GenericPeakInfo &peak = peaks[index];
    
    string isotopes_json;
    double src_counts = 0.0;
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
    {
      const auto pos = relActsColors.find( line.m_isotope );
      if( pos == end(relActsColors) )
      {
        //const double meas_rel_eff = info.m_counts / (info.m_source_gammas[i].m_yield * rel_act);
        
        snprintf( buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G}",
                 (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(), line.m_yield );
      }else
      {
        const double rel_act = pos->second.first;
        src_counts += rel_act * line.m_yield;
        snprintf( buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G, \"rel_act\": %1.6G}",
                 (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(), line.m_yield, rel_act );
      }
      
      isotopes_json += buffer;
    }//for( const RelEff::GammaLineInfo &line : peak.m_source_gammas )
    
    const double eff = peak.m_counts / src_counts;
    double eff_uncert = peak.m_counts_uncert / src_counts;
    
    if( IsNan(eff) || IsInf(eff) )
    {
      cerr << "RelEffChart::setData: Got invalid eff at " << peak.m_energy << ", peak.m_counts="
           << peak.m_counts << ", src_counts=" << src_counts << endl;
      
      continue;
    }
    
    if( IsNan(eff_uncert) || IsInf(eff_uncert) )
    {
      cerr << "RelEffChart::setData: Got invalid eff_uncert at " << peak.m_energy
           << ", peak.m_counts=" << peak.m_counts_uncert << ", src_counts=" << src_counts << endl;
      
      eff_uncert = 0.0;
    }
    
    snprintf( buffer, sizeof(buffer),
             "%s{\"energy\": %.2f, \"counts\": %1.7g, \"counts_uncert\": %1.7g,"
             " \"eff\": %1.6g, \"eff_uncert\": %1.6g, \"nuc_info\": ",
             (njson_entries ? ", " : ""), peak.m_energy, peak.m_counts, peak.m_counts_uncert,
             eff, eff_uncert );
    
    rel_eff_plot_values << buffer;
    rel_eff_plot_values << "[" << isotopes_json.c_str() << "]}";

    njson_entries += 1;
  }//for( size_t index = 0; index < input_peaks.size(); ++index )
  
  rel_eff_plot_values << "]";
  
  if( relEffEqn.empty() )
    relEffEqn = "null";
  
  const string js = m_jsgraph + ".setRelEffData(" + rel_eff_plot_values.str() + "," + relEffEqn + ");";
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void setData( std::shared_ptr<Measurement> data_hist )


void RelEffChart::setCssRules()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  std::shared_ptr<const ColorTheme> theme = interspec ? interspec->getColorTheme() : nullptr;
  
  assert( theme );
  if( !theme )
    return;
  
  WCssStyleSheet &style = wApp->styleSheet();
  
  //m_cssRules[".RelEffPlot"] = style.addRule( ".RelEffPlot", "" );
  string rulename = ".RelEffPlot circle - cursor";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( ".RelEffPlot circle", "cursor: pointer;" );
  
  rulename = "div.RelEffPlotTooltip";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "div.RelEffPlotTooltip", "position: fixed; padding: 6px; font: 12px sans-serif; background: #ffffcc; border: 0px; border-radius: 8px; pointer-events: none; color: #444422;" );
  
  setLineColor( theme->foregroundLine );
  setDefaultMarkerColor( theme->backgroundLine );
  setChartBackgroundColor( theme->spectrumChartBackground );
  setAxisLineColor( theme->spectrumAxisLines );
  setTextColor( theme->spectrumChartText );
}//void setCssRules()


void RelEffChart::setLineColor( const Wt::WColor &color )
{
  WColor lineColor = color.isDefault() ? WColor("steelblue") : color;
  
  WCssStyleSheet &style = wApp->styleSheet();
  string rulename = ".RelEffPlot path";
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( ".RelEffPlot path", "stroke: " + lineColor.cssText() + "; stroke-width: 2; fill: none;" );
  
  
  rulename = ".RelEffPlot .errorbar";
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  m_cssRules[rulename] = style.addRule( ".RelEffPlot .errorbar", "fill: none; stroke-width: 1; stroke: " + lineColor.cssText() + ";" );
}//setLineColor(...)


void RelEffChart::setTextColor( const Wt::WColor &color )
{
  WColor textColor = color.isDefault() ? WColor(0,0,0) : color;
  const string c = textColor.cssText();
  
  const string rulename = "TextColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  m_cssRules[rulename] = style.addRule( ".RelEffPlot .xaxistitle, .RelEffPlot .yaxistitle, .RelEffPlot .yaxis, .RelEffPlot .yaxislabel, .RelEffPlot .xaxis, .RelEffPlot .tick text", "fill: " + c );
}


void RelEffChart::setAxisLineColor( const Wt::WColor &color )
{
  WColor axisColor = color.isDefault() ? WColor(0,0,0) : color;
  
  string rulename = "AxisColor";
  
  WCssStyleSheet &style = wApp->styleSheet();
  if( m_cssRules.count(rulename) )
    style.removeRule( m_cssRules[rulename] );
  
  const char * const selector = ".RelEffPlot .xAxis path, .RelEffPlot .xAxis line, .RelEffPlot .yAxis path, .RelEffPlot .yAxis line";
  const string css = "fill: none; stroke-width: 1; shape-rendering: crispEdges; stroke: " + axisColor.cssText() + ";";
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
  
  m_cssRules[rulename] = style.addRule( ".RelEffPlot circle, .RelEffPlot circle.noiso, .RelEffPlot circle.multiiso",
                                        "fill: " + markerColor.cssText() + ";" );
}


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


