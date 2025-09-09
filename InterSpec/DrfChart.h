#ifndef DrfChart_h
#define DrfChart_h
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

#include <map>
#include <vector>

#include <Wt/WContainerWidget>

struct ColorTheme;
class DetectorPeakResponse;

namespace Wt
{
  class WColor;
  class WCssTextRule;
}

class DrfChart : public Wt::WContainerWidget
{
protected:
  void defineJavaScript();
  void setCssRules();
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );
  
  std::shared_ptr<const DetectorPeakResponse> m_detector;
  
  /** The javascript variable name used to refer to the DrfChart object.
   Currently is `jsRef() + ".chart"`.
   */
  const std::string m_jsgraph;
  
  /** JS calls requested before the widget has been rendered, so wouldnt have
     ended up doing anything are saved here, and then executed once the widget
     is rendered.
   */
  std::vector<std::string> m_pendingJs;
  
  std::map<std::string,Wt::WCssTextRule *> m_cssRules;

  
public:
  DrfChart( Wt::WContainerWidget *parent = 0 );
  virtual ~DrfChart();
  
  void handleColorThemeChange( std::shared_ptr<const ColorTheme> theme );
  
  void updateChart( std::shared_ptr<const DetectorPeakResponse> det );
  
  /** Set the x-axis range from C++ */
  void setXAxisRange( double minEnergy, double maxEnergy );
  
  /** Get the current x-axis range */
  std::pair<double, double> getXAxisRange() const;
  
  /** Set custom colors based on ColorTheme */
  void setLineColors();
  
  void setLineColor( const Wt::WColor &color );
  void setTextColor( const Wt::WColor &color );
  void setAxisLineColor( const Wt::WColor &color );
  void setChartBackgroundColor( const Wt::WColor &color );
  
protected:
  /** Send data to JavaScript chart */
  void sendDataToJavaScript( const std::string &jsonData );
  
  /** Convert detector data to JSON format */
  std::string generateJsonData() const;
  
private:
  double m_minEnergy, m_maxEnergy;
};//class DrfChart

#endif //DrfChart_h

