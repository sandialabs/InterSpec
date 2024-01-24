#ifndef SearchMode3DChart_H
#define SearchMode3DChart_H
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
#include <memory>
#include <string>
#include <vector>

#include <Wt/WContainerWidget>

//Some forward declarations
class SpecMeas;
class InterSpec;
class NativeFloatSpinBox;
class SearchMode3DDataModel;

namespace SpecUtils{ enum class SpectrumType : int; }

namespace Wt
{
  class WText;
  namespace Chart
  {
    class WSpinBox;
    class WGridData;
    class WGridLayout;
    class WDoubleSpinBox;
    class WCartesian3DChart;
  }//namespace Chart
}//namespace Wt


class SearchMode3DChart : public Wt::WContainerWidget
{
public:
  /** Constructor.
     \param viewer Must be a valid pointer to thw 'owning' InterSpec instance.
   */
  SearchMode3DChart( InterSpec *viewer, Wt::WContainerWidget *parent = 0 );

  /** Destructor.  Nothing fancy (no-op) */
  virtual ~SearchMode3DChart();
  
  /** Redraws the chart, and refreshes the various widgets any time the data
      being viewed in InterSpec is changed.
   */
  void updateDisplay();
  
  /** Sets the z-axis (observed counts) as log or linear. */
  void setLogZ( const bool log );
  
  /** Updates the chart to match the energy and time ranges specified (after
      validating the user input).
   */
  void updateRange();

  /** Custom loader to start the loading of the chart, after the widget actually
      starts being rendered.  See initChart() for further explanation.
   */
  virtual void load();
  
protected:
  /** Performs the construction of the widget hierarchy, and initializes
    all necassary components.
   */
  void init();
  
  /** For a weird reason (at least with Wt 3.4.4) if we load the chart to begin
      with, there is a JS exception having to do with WSpinBoxs.  So we will
      load the chart after the other widgets.
   */
  void initChart();
  
  
  /** Sets m_width and m_height so we know the (client side) rendered dimensions
   of this widget, incase this would influence anything.
   */
  virtual void layoutSizeChanged( int width, int height );
  
  /** Sets appropriate limits for the time range input widgets, and makes sure
      they have valid values.
   */
  void setTimeLimits();
  
  /** Sets appropriate limits for the energy range input widgets, and makes sure
      the have valid values.
   */
  void setEnergyLimits();
  
  /** Function that gets called when the user loads a new spectrum. */
  void newSpectralDataSet( const SpecUtils::SpectrumType type,
                          const std::shared_ptr<SpecMeas> &meas,
                          const std::set<int> &sample_numbers,
                          const std::vector<std::string> &detectors );
  
  void setBinningLimits();
  void handleNumTimeDivsChanged();
  void handleNumEnergyDivsChanged();
  
protected:
  InterSpec *m_viewer;
  
  Wt::WGridLayout *m_layout;
  
  int m_width;   //the rendered width of this widget in in pixels
  int m_height;  //the rendered height of this widget in in pixels
  
  bool m_loaded;  //set to true once load() function is called.
  SearchMode3DDataModel *m_model;
  Wt::Chart::WGridData  *m_data;
  Wt::Chart::WCartesian3DChart *m_chart;
  
  NativeFloatSpinBox *m_inputMinEnergy;
  NativeFloatSpinBox *m_inputMaxEnergy;
  NativeFloatSpinBox *m_inputMinTime;
  NativeFloatSpinBox *m_inputMaxTime;
  
  Wt::WSpinBox *m_timeDivisions;
  Wt::WSpinBox *m_energyDivisions;
  
  Wt::WText *m_loadingTxt;
  
  /** Some absolute limits on min and max number of divisions to display in the 3D chart. */
  static const int sm_maxTimeDivs;
  static const int sm_maxEnergyDivs;
  static const int sm_minTimeOrEnergyDivs;
};//class SearchMode3DChart


#endif
