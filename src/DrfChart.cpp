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
  
  // Require JavaScript resources
  wApp->require( "InterSpec_resources/d3.v3.min.js", "d3.v3.js" );
  wApp->require( "InterSpec_resources/DrfChart.js" );
  
  wApp->useStyleSheet( "InterSpec_resources/DrfChart.css" );
}//DrfChart constructor


void DrfChart::updateChart( std::shared_ptr<const DetectorPeakResponse> det )
{
  m_detector = det;
  
  if( !det || !det->isValid() )
  {
    // Send empty data to JavaScript
    const string efficiencyJs = m_jsgraph + ".setEfficiencyData(null);";
    const string fwhmJs = m_jsgraph + ".setFwhmData(null);"; 
    if( isRendered() )
    {
      doJavaScript( efficiencyJs );
      doJavaScript( fwhmJs );
    }else
    {
      m_pendingJs.push_back( efficiencyJs );
      m_pendingJs.push_back( fwhmJs );
    }
    return;
  }
  
  // Generate and send efficiency data for JavaScript calculation
  const string efficiencyData = generateEfficiencyData();
  const string efficiencyJs = m_jsgraph + ".setEfficiencyData(" + efficiencyData + ");";
  if( isRendered() )
    doJavaScript( efficiencyJs );
  else
    m_pendingJs.push_back( efficiencyJs );
  
  // Generate and send FWHM data if available
  if( det->hasResolutionInfo() )
  {
    const string fwhmData = generateFwhmData();
    const string fwhmJs = m_jsgraph + ".setFwhmData(" + fwhmData + ");";
    if( isRendered() )
      doJavaScript( fwhmJs );
    else
      m_pendingJs.push_back( fwhmJs );
  }
  else
  {
    // Clear FWHM data
    const string fwhmJs = m_jsgraph + ".setFwhmData(null);";
    if( isRendered() )
      doJavaScript( fwhmJs );
    else
      m_pendingJs.push_back( fwhmJs );
  }
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


std::string DrfChart::generateEfficiencyData() const
{
  if( !m_detector || !m_detector->isValid() )
    return "null";
    
  stringstream json;
  json << "{";
  
  // Get efficiency form
  const auto efficiencyForm = m_detector->efficiencyFcnType();
  
  // Get energy units
  float energyUnits = m_detector->efficiencyEnergyUnits();
  
  if( efficiencyForm == DetectorPeakResponse::kFunctialEfficienyForm )
  {
    // For functional form, generate 100 points and send as kEnergyEfficiencyPairs with energyUnits=1
    json << "\"form\":\"kEnergyEfficiencyPairs\"";
    json << ",\"energyUnits\":1";
    
    float minEnergy = 50.0f, maxEnergy = 3000.0f;
    if( m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0) )
    {
      minEnergy = m_detector->lowerEnergy();
      maxEnergy = m_detector->upperEnergy();
    }
    
    // Add energy extent to JSON
    json << ",\"energyExtent\":[" << minEnergy << "," << maxEnergy << "]";
    
    json << ",\"pairs\":[";
    bool first = true;
    for( int i = 0; i < 100; ++i )
    {
      const float energy = minEnergy + (float(i)/99.0f) * (maxEnergy-minEnergy);
      const float efficiency = static_cast<float>( m_detector->intrinsicEfficiency( energy ) );
      
      // Skip invalid efficiency values
      if( IsNan(efficiency) || IsInf(efficiency) || efficiency < 0.0f )
        continue;
        
      if( !first )
        json << ",";
      first = false;
      json << "{\"energy\":" << energy << ",\"efficiency\":" << efficiency << "}";
    }
    json << "]";
  }
  else if( efficiencyForm == DetectorPeakResponse::kEnergyEfficiencyPairs )
  {
    json << "\"form\":\"kEnergyEfficiencyPairs\"";
    json << ",\"energyUnits\":" << energyUnits;
    
    const auto& pairs = m_detector->getEnergyEfficiencyPair();
    
    // Calculate energy extent from pairs data
    float minEnergy = 0.0f, maxEnergy = 3000.0f;
    if( !pairs.empty() )
    {
      minEnergy = pairs[0].energy;
      maxEnergy = pairs[0].energy;
      for( const auto& pair : pairs )
      {
        if( pair.energy < minEnergy )
          minEnergy = pair.energy;
        if( pair.energy > maxEnergy )
          maxEnergy = pair.energy;
      }
    }
    
    // Use detector range if available and reasonable, but constrain by actual data range
    if( m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0) )
    {
      minEnergy = std::max(minEnergy, static_cast<float>(m_detector->lowerEnergy()));
      maxEnergy = std::min(maxEnergy, static_cast<float>(m_detector->upperEnergy()));
    }
    
    json << ",\"energyExtent\":[" << minEnergy << "," << maxEnergy << "]";
    
    json << ",\"pairs\":[";
    for( size_t i = 0; i < pairs.size(); ++i )
    {
      if( i > 0 )
        json << ",";
      json << "{\"energy\":" << pairs[i].energy << ",\"efficiency\":" << pairs[i].efficiency << "}";
    }
    json << "]";
  }
  else if( efficiencyForm == DetectorPeakResponse::kExpOfLogPowerSeries )
  {
    json << "\"form\":\"kExpOfLogPowerSeries\"";
    json << ",\"energyUnits\":" << energyUnits;
    
    // Add energy extent to JSON
    float minEnergy = 0.0f, maxEnergy = 3000.0f;
    if( m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0) )
    {
      minEnergy = m_detector->lowerEnergy();
      maxEnergy = m_detector->upperEnergy();
    }
    json << ",\"energyExtent\":[" << minEnergy << "," << maxEnergy << "]";
    
    const auto& coeffs = m_detector->efficiencyExpOfLogsCoeffs();
    json << ",\"coefficients\":[";
    for( size_t i = 0; i < coeffs.size(); ++i )
    {
      if( i > 0 )
        json << ",";
      json << coeffs[i];
    }
    json << "]";
  }
  else
  {
    json << "\"form\":\"unknown\"";
    json << ",\"energyUnits\":1";
  }
  
  json << "}";
  return json.str();
}//std::string DrfChart::generateEfficiencyData()

std::pair<float, float> DrfChart::getEnergyRange() const
{
  if( !m_detector || !m_detector->isValid() )
    return std::make_pair(0.0f, 3000.0f);
    
  float minEnergy = 50.0f, maxEnergy = 3000.0f;
  
  // First, try to get range from EnergyEfficiencyPair data if available
  const auto efficiencyForm = m_detector->efficiencyFcnType();
  if( efficiencyForm == DetectorPeakResponse::kEnergyEfficiencyPairs )
  {
    const auto& pairs = m_detector->getEnergyEfficiencyPair();
    if( !pairs.empty() )
    {
      minEnergy = pairs[0].energy;
      maxEnergy = pairs[0].energy;
      for( const auto& pair : pairs )
      {
        if( pair.energy < minEnergy )
          minEnergy = pair.energy;
        if( pair.energy > maxEnergy )
          maxEnergy = pair.energy;
      }
    }
  }
  
  // Use detector range if available and reasonable, but constrain by actual data range
  if( m_detector->upperEnergy() > (m_detector->lowerEnergy() + 100.0) )
  {
    minEnergy = std::max(minEnergy, static_cast<float>(m_detector->lowerEnergy()));
    maxEnergy = std::min(maxEnergy, static_cast<float>(m_detector->upperEnergy()));
  }
  
  return std::make_pair(minEnergy, maxEnergy);
}

std::string DrfChart::generateFwhmData() const
{
  if( !m_detector || !m_detector->isValid() || !m_detector->hasResolutionInfo() )
    return "null";
    
  stringstream json;
  json << "[";
  
  // Use the same energy range as efficiency data
  const auto energyRange = getEnergyRange();
  const float minEnergy = energyRange.first;
  const float maxEnergy = energyRange.second;
  
  // Generate FWHM data points using the same range as efficiency data
  int numEnergyPoints = static_cast<int>(floor(maxEnergy - minEnergy) / 2.5);
  if( numEnergyPoints > 4500 ) //4500 chosen arbitrarily
    numEnergyPoints = 4500;
    
  bool first = true;
  for( int i = 0; i < numEnergyPoints; ++i )
  {
    const float energy = minEnergy + (float(i)/float(numEnergyPoints)) * (maxEnergy-minEnergy);
    const float fwhm = m_detector->peakResolutionFWHM(energy);
    
    // Skip any points outside where we would expect.
    if( IsNan(fwhm) || IsInf(fwhm) || fwhm < 0.0f || fwhm >= 9999.9f )
      continue;
      
    if( !first )
      json << ",";
    first = false;
    
    json << "{\"energy\":" << energy << ",\"fwhm\":" << fwhm << "}";
  }
  
  json << "]";
  return json.str();
}//std::string DrfChart::generateFwhmData()


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

