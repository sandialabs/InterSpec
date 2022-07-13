#ifndef FeatureMarkerWidget_h
#define FeatureMarkerWidget_h
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
#include <memory>

#include <Wt/WContainerWidget>

#include "InterSpec/AuxWindow.h"

class InterSpec;
class FeatureMarkerWidget;


namespace Wt
{
  class WSpinBox;
  class WCheckBox;
}//namespace Wt


class FeatureMarkerWindow : public AuxWindow
{
public:
  FeatureMarkerWindow( InterSpec* viewer );
  
  virtual ~FeatureMarkerWindow();
  
protected:
  FeatureMarkerWidget *m_feature;
};//class FeatureMarkerWindow



class FeatureMarkerWidget : public Wt::WContainerWidget
{
public:
  FeatureMarkerWidget( InterSpec* viewer, Wt::WContainerWidget *parent = 0 );
  
  virtual ~FeatureMarkerWidget();
  
protected:
  void init();
  
  void handleComptonAngleChanged();
  
  InterSpec *m_viewer;
  
  Wt::WCheckBox *m_escapePeaks;
  Wt::WCheckBox *m_comptonPeak;
  Wt::WSpinBox  *m_comptonAngle;
  Wt::WCheckBox *m_comptonEdge;
  Wt::WCheckBox *m_sumPeaks;
};//class FeatureMarkerWidget

#endif //FeatureMarkerWidget_h
