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

#include <memory>
#include <string>
#include <vector>

#include <Wt/WSignal.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceFitCalc.h"
#include "InterSpec/SimpleDialog.h"

namespace Wt
{
  class WText;
  class WCheckBox;
  class WComboBox;
  class WVBoxLayout;
}

// Forward declarations
class Shielding2DView;
class Shielding3DView;
class DetectorPeakResponse;

// Dialog class for displaying shielding diagrams with 2D/3D view switching
class ShieldingDiagramDialog : public SimpleDialog
{
  friend class SimpleDialog;
public:
  // Static factory method to create a dialog with 2D/3D view switcher.
  // sourceOffset0/sourceOffset1 are the user-set off-axis offsets
  // (GammaInteractionCalc::ShieldSourceConfig::source_offsets[0] and [1]).
  // `drf` lets the 3D view draw the detector's actual geometry, when the DRF records one.
  static ShieldingDiagramDialog *createShieldingDiagram(
                                     const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                     const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                                     GammaInteractionCalc::GeometryType geometry,
                                     double detectorDistance,
                                     double detectorDiameter,
                                     double sourceOffset0 = 0.0,
                                     double sourceOffset1 = 0.0,
                                     std::shared_ptr<const DetectorPeakResponse> drf = nullptr
                                     );

  // Switch between 2D and 3D views
  void switchView( bool show3D );

  // Update the data for both views
  void updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   double sourceOffset0 = 0.0,
                   double sourceOffset1 = 0.0,
                   std::shared_ptr<const DetectorPeakResponse> drf = nullptr );

  /** Emitted with a gamma energy (keV; <= 0 for the default) when the 3D view wants the fit's
   integration lines - the owner answers with #setVolumetricLines or #setVolumetricLinesError. */
  Wt::Signal<double> &volumetricLinesRequested();

  /** Shows these lines (see GammaInteractionCalc::ShieldingSourceChi2Fcn::sampleVolumetricLines). */
  void setVolumetricLines( const GammaInteractionCalc::VolumetricLineSample &sample );

  /** The lines could not be computed: says why, and turns the option off (back to the default
   energy, should that energy be the problem). */
  void setVolumetricLinesError( const Wt::WString &message );

protected:
  // Constructor is protected; use SimpleDialog::make<ShieldingDiagramDialog>() to create.
  ShieldingDiagramDialog(
                         const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                         const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                         GammaInteractionCalc::GeometryType geometry,
                         double detectorDistance,
                         double detectorDiameter,
                         double sourceOffset0,
                         double sourceOffset1,
                         std::shared_ptr<const DetectorPeakResponse> drf
                         );

private:
  void handleViewTypeToggle();
  void handleShowLinesToggled();
  void handleLineEnergyChanged();
  void requestLines();
  bool linesPossible() const;

  Shielding2DView *m_2DView;
  Shielding3DView *m_3DView;
  Wt::WComboBox *m_select;
  Wt::WVBoxLayout *m_layout;
  Wt::WContainerWidget *m_viewHolder;   //the stretching cell the current view fills

  // The 3D view's integration-line controls
  Wt::WContainerWidget *m_linesControls;
  Wt::WCheckBox *m_showLines;
  Wt::WComboBox *m_lineEnergy;
  Wt::WText *m_linesMsg;
  std::vector<double> m_lineEnergies;   //keV, the entries of m_lineEnergy
  std::string m_linesJson;              //what the 3D view was last given; resent when it is rebuilt
  Wt::Signal<double> m_linesRequested;

  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::SourceFitDef> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
  double m_sourceOffsets[2];
  std::shared_ptr<const DetectorPeakResponse> m_drf;
};

// Create JSON representation of shielding data
std::string createShieldingDiagramJson(
                                       const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                                       const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                                       GammaInteractionCalc::GeometryType geometry,
                                       double detectorDistance,
                                       double detectorDiameter,
                                       double sourceOffset0 = 0.0,
                                       double sourceOffset1 = 0.0
                                       );

class Shielding2DView : public Wt::WContainerWidget
{
public:
  Shielding2DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   double sourceOffset0 = 0.0,
                   double sourceOffset1 = 0.0 );

  // Update the data and refresh the display
  void updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   double sourceOffset0 = 0.0,
                   double sourceOffset1 = 0.0 );

private:
  void defineJavaScript();
  std::string createJsonData() const;

  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::SourceFitDef> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
  double m_sourceOffsets[2];
};

// Shielding3DView class (merged from Shielding3DView.h)
class Shielding3DView : public Wt::WContainerWidget
{
public:
  Shielding3DView( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   double sourceOffset0 = 0.0,
                   double sourceOffset1 = 0.0,
                   std::shared_ptr<const DetectorPeakResponse> drf = nullptr );

  // Update the data and refresh the display
  void updateData( const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                   const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources,
                   GammaInteractionCalc::GeometryType geometry,
                   double detectorDistance,
                   double detectorDiameter,
                   double sourceOffset0 = 0.0,
                   double sourceOffset1 = 0.0,
                   std::shared_ptr<const DetectorPeakResponse> drf = nullptr );

  /** The detector's geometry for the 3D view: `null` when the DRF records none (the JS then draws
   a placeholder), else each region of `DetectorGeometryDiagram::buildModel` as a polycone profile
   in the crystal frame (cm), plus the endcap-front offset that places the crystal face behind the
   detector face the source distance is measured to.
   */
  static std::string createDetectorJson( const std::shared_ptr<const DetectorPeakResponse> &drf );

  /** Shows integration lines (JSON from the dialog), or none for "null". */
  void setVolumetricLines( const std::string &json );

private:
  void defineJavaScript();
  std::string createJsonData() const;

  std::vector<ShieldingSourceFitCalc::ShieldingInfo> m_shieldings;
  std::vector<ShieldingSourceFitCalc::SourceFitDef> m_sources;
  GammaInteractionCalc::GeometryType m_geometry;
  double m_detectorDistance;
  double m_detectorDiameter;
  double m_sourceOffsets[2];
  std::shared_ptr<const DetectorPeakResponse> m_drf;
};

#endif // Shielding2DView_h

