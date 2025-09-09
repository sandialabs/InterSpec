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


#include <memory>
#include <sstream>
#include <string>

#include <Wt/Utils>
#include <Wt/WColor>
#include <Wt/WLength>
#include <Wt/WJavaScript>
#include <Wt/WApplication>
#include <Wt/WStringStream>
#include <Wt/WCssStyleSheet>
#include <Wt/WContainerWidget>

#include "InterSpec/DrfChart.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/DetectorPeakResponse.h"

using namespace std;
using namespace Wt;


DrfChart::DrfChart( WContainerWidget *parent )
: WContainerWidget( parent ),
  m_detector( nullptr ),
  m_jsgraph( jsRef() + ".chart" ),
  m_minEnergy( 0.0 ),
  m_maxEnergy( 3000.0 )
{
  addStyleClass( "DrfChart" );
  setOverflow( Overflow::OverflowHidden );
  
  // Setup CSS rules
  setCssRules();
  
  // Connect to color theme changes
  InterSpec *interspec = InterSpec::instance();
  if( interspec )
    interspec->colorThemeChanged().connect( this, &DrfChart::setCssRules );
  
  // Require JavaScript resources
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DrfChart.js" );
  
  // Apply CSS
  wApp->useStyleSheet( "InterSpec_resources/DrfChart.css" );
}//DrfChart constructor


void DrfChart::handleColorThemeChange( std::shared_ptr<const ColorTheme> theme )
{
  setCssRules();
}//handleColorThemeChange(...)


void DrfChart::updateChart( std::shared_ptr<const DetectorPeakResponse> det )
{
  m_detector = det;
  
  if( !det || !det->isValid() )
  {
    // Send empty data to JavaScript
    const string js = m_jsgraph + ".setData([]);"; 
    if( isRendered() )
      doJavaScript( js );
    else
      m_pendingJs.push_back( js );
    return;
  }
  
  // Generate JSON data and send to JavaScript
  const string jsonData = generateJsonData();
  sendDataToJavaScript( jsonData );
} //DrfChart::updateChart()


void DrfChart::defineJavaScript()
{
  string options = "{ "
    "margins: {"
    " top: 5,"
    " right: 50,"
    " bottom: 40,"
    " left: 60"
    " } }";
  
  setJavaScriptMember( "chart", "new DrfChart(" + jsRef() + ", " + options + ");");
  
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
  
  // Execute any pending JS calls
  for( const string &js : m_pendingJs )
    doJavaScript( js );
  m_pendingJs.clear();
  m_pendingJs.shrink_to_fit();
}//void DrfChart::defineJavaScript()


void DrfChart::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( renderFull )
    defineJavaScript();
}//void DrfChart::render(...)


void DrfChart::setCssRules()
{
  InterSpec *interspec = InterSpec::instance();
  std::shared_ptr<const ColorTheme> theme = interspec ? interspec->getColorTheme() : nullptr;
  
  if( !theme )
    return;
  
  // Set CSS variables for the theme colors
  WCssStyleSheet &style = wApp->styleSheet();
  
  // Remove any existing CSS rules we may have added
  for( auto &p : m_cssRules )
    style.removeRule( p.second );
  m_cssRules.clear();
  
  // Set CSS variables that the DrfChart.css will use
  const string selector = "#" + id();
  string cssProps;
  cssProps += "--d3spec-fore-line-color: " + (theme->foregroundLine.isDefault() ? "#B02B2C" : theme->foregroundLine.cssText()) + "; ";
  cssProps += "--d3spec-back-line-color: " + (theme->backgroundLine.isDefault() ? "#3F4C6B" : theme->backgroundLine.cssText()) + "; ";
  cssProps += "--d3spec-axis-color: " + (theme->spectrumAxisLines.isDefault() ? "black" : theme->spectrumAxisLines.cssText()) + "; ";
  cssProps += "--d3spec-text-color: " + (theme->spectrumChartText.isDefault() ? "black" : theme->spectrumChartText.cssText()) + "; ";
  cssProps += "--d3spec-background-color: " + (theme->spectrumChartBackground.isDefault() ? "transparent" : theme->spectrumChartBackground.cssText()) + "; ";
  cssProps += "--d3spec-chart-area-color: " + (theme->spectrumChartBackground.isDefault() ? "rgba(0,0,0,0)" : theme->spectrumChartBackground.cssText()) + ";";
  
  m_cssRules["css-variables"] = style.addRule( selector, cssProps );
}//void DrfChart::setCssRules()


void DrfChart::setLineColor( const Wt::WColor &color )
{
  // This method is no longer needed since we use CSS variables
  // Keep it for API compatibility but don't do anything
}//void DrfChart::setLineColor(...)


void DrfChart::setTextColor( const Wt::WColor &color )
{
  // This method is no longer needed since we use CSS variables
  // Keep it for API compatibility but don't do anything
}//void DrfChart::setTextColor(...)


void DrfChart::setAxisLineColor( const Wt::WColor &color )
{
  // This method is no longer needed since we use CSS variables
  // Keep it for API compatibility but don't do anything
}//void DrfChart::setAxisLineColor(...)


void DrfChart::setChartBackgroundColor( const Wt::WColor &color )
{
  // This method is no longer needed since we use CSS variables
  // Keep it for API compatibility but don't do anything
}//void DrfChart::setChartBackgroundColor(...)


void DrfChart::setLineColors()
{
  // Line colors are now handled entirely via CSS variables
  // No JavaScript calls needed since CSS will handle the styling automatically
}//void DrfChart::setLineColors()


void DrfChart::sendDataToJavaScript( const std::string &jsonData )
{
  const string js = m_jsgraph + ".setData(" + jsonData + ");";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void DrfChart::sendDataToJavaScript(...)


std::string DrfChart::generateJsonData() const
{
  if( !m_detector || !m_detector->isValid() )
    return "[]";
    
  stringstream json;
  json << "[";
  
  const bool hasResolution = m_detector->hasResolutionInfo();
  
  float minEnergy = 0.0f, maxEnergy = 3000.0f;
  
  // If DRF has a defined energy range, of at least 100 keV, use it.
  if( m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0) )
  {
    minEnergy = m_detector->lowerEnergy();
    maxEnergy = m_detector->upperEnergy();
  }
  
  // Set energy range for x-axis
  const_cast<DrfChart*>(this)->m_minEnergy = minEnergy;
  const_cast<DrfChart*>(this)->m_maxEnergy = maxEnergy;
  
  // TODO: pick this better using the chart width or whatever
  int numEnergyPoints = static_cast<int>(floor(maxEnergy - minEnergy) / 2.5);
  if( numEnergyPoints > 4500 ) //4500 chosen arbitrarily
    numEnergyPoints = 4500;
    
  bool first = true;
  for( int i = 0; i < numEnergyPoints; ++i )
  {
    const float energy = minEnergy + (float(i)/float(numEnergyPoints)) * (maxEnergy-minEnergy);
    const float efficiency = static_cast<float>( m_detector->intrinsicEfficiency( energy ) );
    
    // Skip any points outside where we would expect.
    if( IsNan(efficiency) || IsInf(efficiency) || efficiency < 0.0f )
      continue;
      
    if( !first )
      json << ",";
    first = false;
    
    json << "{\"energy\":" << energy << ",\"efficiency\":" << efficiency;
    
    if( hasResolution )
    {
      const float fwhm = m_detector->peakResolutionFWHM(energy);
      if( !IsNan(fwhm) && !IsInf(fwhm) && fwhm >= 0.0f && fwhm < 9999.9f )
      {
        json << ",\"fwhm\":" << fwhm;
      }
    }
    
    json << "}";
  }
  
  json << "]";
  return json.str();
}//std::string DrfChart::generateJsonData()


void DrfChart::setXAxisRange( double minEnergy, double maxEnergy )
{
  m_minEnergy = minEnergy;
  m_maxEnergy = maxEnergy;
  
  const string js = m_jsgraph + ".setXRange(" + std::to_string(minEnergy) + ", " + std::to_string(maxEnergy) + ");";
  
  if( isRendered() )
    doJavaScript( js );
  else
    m_pendingJs.push_back( js );
}//void DrfChart::setXAxisRange(...)


std::pair<double, double> DrfChart::getXAxisRange() const
{
  return std::make_pair( m_minEnergy, m_maxEnergy );
}//std::pair<double, double> DrfChart::getXAxisRange()


DrfChart::~DrfChart()
{
  Wt::WApplication *app = Wt::WApplication::instance();
  if( app )
  {
    WCssStyleSheet &style = app->styleSheet();
    for( const auto &rule : m_cssRules )
      style.removeRule( rule.second );
    m_cssRules.clear();
  }//if( app )
}//~DrfChart()

