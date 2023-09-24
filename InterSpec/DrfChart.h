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

#include <Wt/WColor>
#include <Wt/WContainerWidget>
#include <Wt/Chart/WCartesianChart>

class ColorTheme;
class DetectorPeakResponse;

namespace Wt
{
  class WStandardItemModel;
}

class DrfChart : public Wt::Chart::WCartesianChart
{
protected:
  Wt::WColor m_chartEnergyLineColor;
  Wt::WColor m_chartFwhmLineColor;
  Wt::WStandardItemModel* m_efficiencyModel;
  
  std::shared_ptr<const DetectorPeakResponse> m_detector;

  
public:
  DrfChart( Wt::WContainerWidget *parent = 0 );
  
  void handleColorThemeChange( std::shared_ptr<const ColorTheme> theme );
  
  void updateChart( std::shared_ptr<const DetectorPeakResponse> det );
};//class DrfChart

#endif //DrfChart_h

