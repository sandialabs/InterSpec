#ifndef DetectorGeometryInput_h
#define DetectorGeometryInput_h
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
#include <vector>

#include <Wt/WContainerWidget.h>

class InterSpec;
class DetectorPeakResponse;
class ShieldMaterialSuggestion;

namespace Wt
{
  class WText;
  class WTable;
  class WTableRow;
  class WLabel;
  class WCheckBox;
  class WComboBox;
  class WLineEdit;
  class WPushButton;
}//namespace Wt

namespace ceelo{ struct GeometryDescriptor; }

/** A form for entering/editing the physical geometry of a detector - shape,
 crystal dimensions and material, bore hole (coax HPGe), dead layers, one or
 more concentric endcap/housing layers, an optional collimator, and the
 distance reference point - i.e., everything a Monte-Carlo characterization
 of the detector response needs.

 Converts to/from ceelo::GeometryDescriptor.  Used by the
 "Make MC Response" tool (see MakeMcResponseForDrf.h), but reusable anywhere
 detector geometry entry is needed.
 */
class DetectorGeometryInput : public Wt::WContainerWidget
{
public:
  DetectorGeometryInput( InterSpec *viewer );
  virtual ~DetectorGeometryInput() override;

  /** Builds the geometry descriptor from the current inputs.
   Throws std::runtime_error with a user-displayable message when input is
   missing or invalid.
   */
  ceelo::GeometryDescriptor toDescriptor() const;

  /** Sets the GUI state from a descriptor (e.g., an existing MC responses
   geometry).  Layer materials are matched back to their stored
   name/density/composition.
   */
  void setFromDescriptor( const ceelo::GeometryDescriptor &descriptor );

  /** Seeds the form from a legacy DRF: a cylinder with diameter =
   `drf->detectorDiameter()` and (estimated) length = diameter; a note is
   shown that the length is a guess.  Does nothing for null/invalid DRFs.
   */
  void seedFromDrf( std::shared_ptr<const DetectorPeakResponse> drf );

  /** Whether toDescriptor() would currently succeed. */
  bool isValid() const;

  /** Emitted on any user edit (after validity re-evaluation). */
  Wt::Signal<> &changed();

protected:
  void init();
  void handleShapeChange();
  void handleUserInput();
  void addLayerRow( const Wt::WString &material, const Wt::WString &frontThick,
                    const Wt::WString &sideThick );
  void removeLayerRow();

  /** One concentric endcap/housing layer input row. */
  struct LayerRow
  {
    Wt::WLineEdit *material;
    Wt::WLineEdit *frontThickness;
    Wt::WLineEdit *sideThickness;
  };//struct LayerRow

  InterSpec *m_interspec;

  Wt::WComboBox *m_shape;            //Cylinder | Coaxial HPGe | Rectangular
  Wt::WComboBox *m_crystalMaterial;  //HPGe | NaI | CZT | LaBr3

  Wt::WLabel *m_dim1Label, *m_dim2Label, *m_dim3Label;
  Wt::WLineEdit *m_dim1, *m_dim2, *m_dim3;
  Wt::WTableRow *m_dim3Row;   //shown only for boxes

  Wt::WLabel *m_bulletLabel;
  Wt::WLineEdit *m_bulletRadius;
  Wt::WTableRow *m_bulletRow;   //hidden for boxes; a fillet is a cylinder feature

  Wt::WLabel *m_boreLabel;
  Wt::WLineEdit *m_boreDiam, *m_boreDepth;
  Wt::WCheckBox *m_boreRounded;  //in the bore row, so it follows its visibility
  Wt::WTableRow *m_boreRow;   //shown only for coaxial
  Wt::WLabel *m_deadLabel;
  Wt::WLineEdit *m_deadFront, *m_deadSide;
  Wt::WTableRow *m_deadRow;   //hidden for boxes

  Wt::WTable *m_layersTable;
  std::vector<LayerRow> m_layers;
  Wt::WPushButton *m_addLayer, *m_removeLayer;

  Wt::WCheckBox *m_hasCollimator;
  Wt::WContainerWidget *m_collimatorRow;
  Wt::WLineEdit *m_collimatorMaterial, *m_collimatorThickness, *m_collimatorExtension;

  Wt::WComboBox *m_referencePoint;   //Endcap front | Crystal face

  Wt::WText *m_note;

  ShieldMaterialSuggestion *m_materialSuggestion;

  Wt::Signal<> m_changed;
};//class DetectorGeometryInput

#endif //DetectorGeometryInput_h
