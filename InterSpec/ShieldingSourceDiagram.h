#ifndef Shielding2DView_h
#define Shielding2DView_h
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

#include <Wt/WContainerWidget>

#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/SimpleDialog.h"

namespace Wt
{
  class WComboBox;
  class WGridLayout;
}

// Forward declarations
class Shielding2DView;
class Shielding3DView;

// Dialog class for displaying shielding diagrams with 2D/3D view switching
class ShieldingDiagramDialog : public SimpleDialog
{
public:
  // Static factory method to create a dialog with 2D/3D view switcher
  static ShieldingDiagramDialog *createShieldingDiagram(
                                     const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                     const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                     GammaInteractionCalc::GeometryType geometry,
                                     double detectorDistance,
                                     double detectorDiameter
                                     );
  
  // Switch between 2D and 3D views
  void switchView( bool show3D );
  
private:
  // Private constructor - use createShieldingDiagram() instead
  ShieldingDiagramDialog(
                         const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                         const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                         GammaInteractionCalc::GeometryType geometry,
                         double detectorDistance,
                         double detectorDiameter
                         );
  
  void handleViewTypeToggle();
  
  Shielding2DView *m_2DView;
  Shielding3DView *m_3DView;
  Wt::WComboBox *m_select;
  Wt::WGridLayout *m_layout;
  
  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::IsoFitStruct> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
};

// Create JSON representation of shielding data
std::string createShieldingDiagramJson(
                                       const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                       const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                                       GammaInteractionCalc::GeometryType geometry,
                                       double detectorDistance,
                                       double detectorDiameter
                                       );

class Shielding2DView : public Wt::WContainerWidget
{
public:
  Shielding2DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   Wt::WContainerWidget *parent = 0 );
  
private:
  void defineJavaScript();
  std::string createJsonData() const;

  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::IsoFitStruct> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
};

// Shielding3DView class (merged from Shielding3DView.h)
class Shielding3DView : public Wt::WContainerWidget
{
public:
  Shielding3DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::IsoFitStruct> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   Wt::WContainerWidget *parent = 0 );
  
private:
  void defineJavaScript();
  std::string createJsonData() const;

  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::IsoFitStruct> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
};

#endif // Shielding2DView_h

