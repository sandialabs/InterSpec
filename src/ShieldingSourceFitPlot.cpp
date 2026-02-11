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
#include <vector>
#include <sstream>

#include "rapidxml/rapidxml.hpp"

#include "3rdparty/date/include/date/date.h"

#include <Wt/Utils>
#include <Wt/WColor>
#include <Wt/WString>
#include <Wt/WApplication>
#include <Wt/WJavaScript>
#include <Wt/WStringStream>
#include <Wt/WCssStyleSheet>
#include <Wt/WContainerWidget>

#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ShieldingSourceFitPlot.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "3rdparty/inja/inja.hpp"

using namespace Wt;
using namespace std;


ShieldingSourceFitPlot::ShieldingSourceFitPlot( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_jsgraph( jsRef() + ".chart" ),
  m_xAxisTitle{},
  m_yAxisTitleChi{},
  m_yAxisTitleMult{},
  m_topMargin( 5 ),
  m_rightMargin( 5 ),
  m_bottomMargin( 0 ),
  m_leftMargin( 5 ),
  m_titlePadding( -3 ),
  m_showChi( true )
{
  addStyleClass( "ShieldingSourceFitPlot" );
  setOverflow( Overflow::OverflowHidden );

  wApp->useStyleSheet( "InterSpec_resources/ShieldingSourceFitPlot.css" );
    
  // Set default localized titles
  m_xAxisTitle = WString::tr("Energy (keV)");
  m_yAxisTitleChi = WString::tr("x2g-yaxis-title-chi");
  m_yAxisTitleMult = WString::tr("x2g-yaxis-title-mult");

  setCssRules();
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  if( interspec )
  {
    interspec->useMessageResourceBundle( "ShieldingSourceDisplay" );
    interspec->colorThemeChanged().connect( this, &ShieldingSourceFitPlot::setCssRules );
  }
    
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/ShieldingSourceFitPlot.js" );
}//ShieldingSourceFitPlot constructor


void ShieldingSourceFitPlot::defineJavaScript()
{
  const WString yAxisTitle = m_showChi ? m_yAxisTitleChi : m_yAxisTitleMult;

  string options = "{ "
  "margins: {"
  " top: " + std::to_string(m_topMargin) + ","
  " right: " + std::to_string(m_rightMargin) + ","
  " bottom: " + std::to_string(m_bottomMargin) + ","
  " left: " + std::to_string(m_leftMargin) +
  " }, "
  " titleToAxisPadding: " + std::to_string(m_titlePadding) + ","
  " showChi: " + (m_showChi ? "true" : "false");

  if( !m_xAxisTitle.empty() )
    options += ", xAxisTitle: " + m_xAxisTitle.jsStringLiteral();
  if( !yAxisTitle.empty() )
    options += ", yAxisTitle: " + yAxisTitle.jsStringLiteral();
  options += " }";

  setJavaScriptMember( "chart", "new ShieldingSourceFitPlot(" + jsRef() + ", " + options + ");" );

  // Send all localized strings to JavaScript on first render
  doJavaScript( m_jsgraph + ".setLocalizations( " + localizedStringsJson() + " );" );

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

  // Initialize JSignals for communication from JavaScript to C++
  m_displayModeChangedJS.reset( new JSignal<bool>( this, "displayModeChanged", true ) );
  m_displayModeChangedJS->connect( boost::bind( &ShieldingSourceFitPlot::handleDisplayModeChanged, this, boost::placeholders::_1 ) );

  m_dataPointClickedJS.reset( new JSignal<double>( this, "dataPointClicked", true ) );
  m_dataPointClickedJS->connect( boost::bind( &ShieldingSourceFitPlot::handleDataPointClicked, this, boost::placeholders::_1 ) );

  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void defineJavaScript()


void ShieldingSourceFitPlot::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);

  WContainerWidget::render( flags );

  if( renderFull )
    defineJavaScript();
}


void ShieldingSourceFitPlot::refresh()
{
  WContainerWidget::refresh();

  // Update localized titles
  m_xAxisTitle = WString::tr("Energy (keV)");
  m_yAxisTitleChi = WString::tr("x2g-yaxis-title-chi");
  m_yAxisTitleMult = WString::tr("x2g-yaxis-title-mult");

  if( isRendered() )
    doJavaScript( m_jsgraph + ".setLocalizations( " + localizedStringsJson() + " );" );
}


std::string ShieldingSourceFitPlot::localizedStringsJson()
{
  stringstream json;
  json << "{"
       << "\"xAxisTitle\": " << m_xAxisTitle.jsStringLiteral() << ","
       << "\"yAxisTitleChi\": " << m_yAxisTitleChi.jsStringLiteral() << ","
       << "\"yAxisTitleMult\": " << m_yAxisTitleMult.jsStringLiteral() << ","
       << "\"tooltipChi\": " << WString::tr("x2g-tt-chi-axis").jsStringLiteral() << ","
       << "\"tooltipMult\": " << WString::tr("x2g-tt-scale-axis").jsStringLiteral() << ","
       << "\"ttNumObserved\": " << WString::tr("x2g-tt-num-observed").jsStringLiteral() << ","
       << "\"ttNumExpected\": " << WString::tr("x2g-tt-num-expected").jsStringLiteral() << ","
       << "\"ttNumSigmaOff\": " << WString::tr("x2g-tt-num-sigma-off").jsStringLiteral() << ","
       << "\"ttObsOverExp\": " << WString::tr("x2g-tt-obs-over-exp").jsStringLiteral() << ","
       << "\"ttNumForeground\": " << WString::tr("x2g-tt-num-foreground").jsStringLiteral() << ","
       << "\"ttNumBackground\": " << WString::tr("x2g-tt-num-background").jsStringLiteral() << ","
       << "\"ttSourcesContributing\": " << WString::tr("x2g-tt-sources-contributing").jsStringLiteral() << ","
       << "\"ttCountsAbbrev\": " << WString::tr("x2g-tt-counts-abbrev").jsStringLiteral() << ","
       << "\"ttBrAbbrev\": " << WString::tr("x2g-tt-br-abbrev").jsStringLiteral() << ","
       << "\"devLabel\": " << WString::tr("x2g-dev-label").jsStringLiteral()
       << "}";
  return json.str();
}


void ShieldingSourceFitPlot::setData( const ShieldingSourceFitCalc::ModelFitResults &results )
{
  std::string jsonData = jsonForData( results );
  sendDataToJavaScript( jsonData );
}


std::string ShieldingSourceFitPlot::jsonForData( const ShieldingSourceFitCalc::ModelFitResults &results )
{
  nlohmann::json json_obj;
  json_obj["data_points"] = nlohmann::json::array();

  const vector<GammaInteractionCalc::PeakResultPlotInfo> *peak_comparisons = results.peak_comparisons.get();
  const vector<GammaInteractionCalc::PeakDetail> *peak_details = results.peak_calc_details.get();

  if( !peak_comparisons || !peak_details )
  {
    return json_obj.dump();
  }

  // Ensure we have matching data
  if( peak_comparisons->size() != peak_details->size() )
  {
    cerr << "ShieldingSourceFitPlot::jsonForData: peak_comparisons and peak_details size mismatch: "
         << peak_comparisons->size() << " vs " << peak_details->size() << endl;
  }

  const size_t npoints = std::min( peak_comparisons->size(), peak_details->size() );

  for( size_t i = 0; i < npoints; ++i )
  {
    const GammaInteractionCalc::PeakResultPlotInfo &comparison = (*peak_comparisons)[i];
    const GammaInteractionCalc::PeakDetail &detail = (*peak_details)[i];

    nlohmann::json point;
    point["energy"] = comparison.energy;
    point["chi"] = comparison.numSigmaOff;
    point["mult"] = comparison.observedOverExpected;
    point["mult_uncert"] = comparison.observedOverExpectedUncert;
    point["color"] = comparison.peakColor.cssText(false);
    point["nuclide"] = detail.assignedNuclide;
    point["numForeground"] = comparison.foregroundCounts;
    point["numForegroundUncert"] = comparison.foregroundUncert;
    point["numObserved"] = comparison.observedCounts;
    point["numObservedUncert"] = comparison.observedUncert;
    point["numExpected"] = comparison.expectedCounts;
    point["numBackground"] = comparison.backgroundCounts;
    point["numBackgroundUncert"] = comparison.backgroundUncert;

    // Add source contributions
    point["sources"] = nlohmann::json::array();
    for( size_t j = 0; j < detail.m_sources.size(); ++j )
    {
      const GammaInteractionCalc::PeakDetailSrc &src = detail.m_sources[j];

      // Calculate fraction that this source contributes to the peak
      double totalCountsAtSource = 0.0;
      for( const GammaInteractionCalc::PeakDetailSrc &s : detail.m_sources )
        totalCountsAtSource += s.countsAtSource;

      const double fraction = (totalCountsAtSource > 0.0) ? (src.countsAtSource / totalCountsAtSource) : 0.0;

      const string nuclideSymbol = src.nuclide ? src.nuclide->symbol : string("unknown");

      nlohmann::json source;
      source["nuclide"] = nuclideSymbol;
      source["fraction"] = fraction;
      source["counts"] = src.modelContribToPeak;
      source["br"] = src.br;
      source["energy"] = src.energy;
      point["sources"].push_back( source );
    }

    json_obj["data_points"].push_back( point );
  }

  return json_obj.dump();
}


void ShieldingSourceFitPlot::sendDataToJavaScript( const std::string &jsonData )
{
  const string js = m_jsgraph + ".setData(" + jsonData + ");";

  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}


void ShieldingSourceFitPlot::setShowChi( bool show_chi )
{
  m_showChi = show_chi;

  const string js = m_jsgraph + ".setShowChi(" + (show_chi ? "true" : "false") + ");";

  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}


void ShieldingSourceFitPlot::setXAxisTitle( const Wt::WString &title )
{
  m_xAxisTitle = title;
  if( isRendered() )
    doJavaScript( m_jsgraph + ".setXAxisTitle(" + title.jsStringLiteral() + ", true);" );
}


void ShieldingSourceFitPlot::setYAxisTitleChi( const Wt::WString &title )
{
  m_yAxisTitleChi = title;
  if( isRendered() && m_showChi )
    doJavaScript( m_jsgraph + ".setYAxisTitle(" + title.jsStringLiteral() + ", true);" );
}


void ShieldingSourceFitPlot::setYAxisTitleMult( const Wt::WString &title )
{
  m_yAxisTitleMult = title;
  if( isRendered() && !m_showChi )
    doJavaScript( m_jsgraph + ".setYAxisTitle(" + title.jsStringLiteral() + ", true);" );
}


void ShieldingSourceFitPlot::setContentMargins( int top, int right, int bottom, int left )
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
}//void setContentMargins( int top, int right, int bottom, int left )


Wt::Signal<bool> &ShieldingSourceFitPlot::displayModeChanged()
{
  return m_displayModeChanged;
}


Wt::Signal<double> &ShieldingSourceFitPlot::dataPointClicked()
{
  return m_dataPointClicked;
}


void ShieldingSourceFitPlot::setCssRules()
{
  InterSpec *interspec = InterSpec::instance();
  assert( interspec );
  std::shared_ptr<const ColorTheme> theme = interspec ? interspec->getColorTheme() : nullptr;

  assert( theme );
  if( !theme )
    return;

  WCssStyleSheet &style = wApp->styleSheet();

  // Set tooltip styling
  string rulename = "div.ShieldingSourceFitPlotTooltip";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " div.ShieldingSourceFitPlotTooltip",
                                         "position: fixed;"
                                         " padding: 6px;"
                                         " font: 12px sans-serif;"
                                         " background: #ffffcc;"
                                         " border: 0px;"
                                         " border-radius: 8px;"
                                         " pointer-events: none;"
                                         " color: #444422;" );

  // Set circle cursor
  rulename = ".ShieldingSourceFitPlot circle - cursor";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .ShieldingSourceFitPlot circle", "cursor: pointer;" );

  // Set y-axis area cursor
  rulename = ".ShieldingSourceFitPlot .yaxisarea - cursor";
  if( !m_cssRules.count(rulename) )
    m_cssRules[rulename] = style.addRule( "#" + id() + " .ShieldingSourceFitPlot .yaxisarea", "cursor: pointer;" );
}//void setCssRules()


void ShieldingSourceFitPlot::handleDisplayModeChanged( bool showChi )
{
  m_displayModeChanged.emit( showChi );
}//void handleDisplayModeChanged( bool showChi )


void ShieldingSourceFitPlot::handleDataPointClicked( double energy )
{
  m_dataPointClicked.emit( energy );
}//void handleDataPointClicked( double energy )


ShieldingSourceFitPlot::~ShieldingSourceFitPlot()
{
  Wt::WApplication *app = Wt::WApplication::instance();
  if( app )
  {
    WCssStyleSheet &style = app->styleSheet();
    for( const auto &rule : m_cssRules )
      style.removeRule( rule.second );
    m_cssRules.clear();
  }//if( app )
}//~ShieldingSourceFitPlot()
