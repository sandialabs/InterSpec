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

#include <map>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

#include <Wt/WText.h>
#include <Wt/WTable.h>
#include <Wt/WTableRow.h>
#include <Wt/WTableCell.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/WComboBox.h>
#include <Wt/WLineEdit.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WSuggestionPopup.h>
#include <Wt/WContainerWidget.h>

#include "SandiaDecay.h"

// CeeLo (external_libs/CeeLo/src)
#include "io/DetectorResponse.h"
#include "materials/Material.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/CeeLoUtils.h"
#include "InterSpec/DetectorGeometryInput.h"
#include "InterSpec/ShieldMaterialSuggestion.h"

using namespace Wt;
using namespace std;

namespace
{
  /** Parses a user distance string to cm; throws std::runtime_error with the
   field name on invalid input.
   */
  double distance_cm( const Wt::WLineEdit *edit, const char *field,
                      const bool allow_zero_or_empty )
  {
    const string txt = edit->text().toUTF8();
    if( txt.empty() )
    {
      if( allow_zero_or_empty )
        return 0.0;
      throw runtime_error( WString::tr("dgi-err-missing").arg(field).toUTF8() );
    }

    try
    {
      const double dist = PhysicalUnits::stringToDistance( txt );
      if( (dist < 0.0) || (!allow_zero_or_empty && (dist <= 0.0)) )
        throw runtime_error( "" );
      return dist / PhysicalUnits::cm;
    }catch( std::exception & )
    {
      throw runtime_error( WString::tr("dgi-err-invalid").arg(field).toUTF8() );
    }
  }//distance_cm(...)





  /** The crystal-material combo entries, mapped to CeeLo built-ins. */
  ceelo::MaterialSpec crystal_material_spec( const int index )
  {
    switch( index )
    {
      case 0: return ceelo::MaterialSpec::from( ceelo::make_HPGe() );
      case 1: return ceelo::MaterialSpec::from( ceelo::make_NaI() );
      case 2: return ceelo::MaterialSpec::from( ceelo::make_CZT() );
      case 3: return ceelo::MaterialSpec::from( ceelo::make_LaBr3() );
    }
    throw runtime_error( "Please select a crystal material." );
  }//crystal_material_spec(...)


  /** Best-effort match of a stored MaterialSpec back to the crystal combo. */
  int crystal_index_for_name( const std::string &name )
  {
    for( int i = 0; i < 4; ++i )
    {
      try
      {
        if( crystal_material_spec(i).name == name )
          return i;
      }catch( std::exception & )
      {
      }
    }
    return -1;
  }//crystal_index_for_name(...)


  string cm_to_str( const double cm )
  {
    return PhysicalUnits::printToBestLengthUnits( cm * PhysicalUnits::cm, 4 );
  }
}//namespace


DetectorGeometryInput::DetectorGeometryInput( InterSpec *viewer )
  : WContainerWidget(),
    m_interspec( viewer ),
    m_shape( nullptr ),
    m_crystalMaterial( nullptr ),
    m_dim1Label( nullptr ), m_dim2Label( nullptr ), m_dim3Label( nullptr ),
    m_dim1( nullptr ), m_dim2( nullptr ), m_dim3( nullptr ),
    m_dim3Row( nullptr ),
    m_boreLabel( nullptr ),
    m_boreDiam( nullptr ), m_boreDepth( nullptr ),
    m_boreRow( nullptr ),
    m_deadLabel( nullptr ),
    m_deadFront( nullptr ), m_deadSide( nullptr ),
    m_deadRow( nullptr ),
    m_layersTable( nullptr ),
    m_addLayer( nullptr ), m_removeLayer( nullptr ),
    m_hasCollimator( nullptr ),
    m_collimatorRow( nullptr ),
    m_collimatorMaterial( nullptr ), m_collimatorThickness( nullptr ),
    m_collimatorExtension( nullptr ),
    m_referencePoint( nullptr ),
    m_note( nullptr ),
    m_materialSuggestion( nullptr ),
    m_changed()
{
  assert( m_interspec );
  if( m_interspec )
    m_interspec->useMessageResourceBundle( "MakeMcResponseForDrf" );

  addStyleClass( "DetectorGeometryInput" );

  init();
}//DetectorGeometryInput constructor


DetectorGeometryInput::~DetectorGeometryInput()
{
}


void DetectorGeometryInput::init()
{
  m_materialSuggestion = addChild( std::make_unique<ShieldMaterialSuggestion>() );

  WTable *table = addNew<WTable>();
  table->addStyleClass( "DgiTable" );

  int row = 0;
  auto add_row = [table,&row]( const WString &label ) -> WTableCell * {
    table->elementAt( row, 0 )->addNew<WLabel>( label );
    WTableCell *cell = table->elementAt( row, 1 );
    ++row;
    return cell;
  };

  auto make_dist_edit = []( WTableCell *cell ) -> WLineEdit * {
    WLineEdit *edit = cell->addNew<WLineEdit>();
    edit->setAutoComplete( false );
    edit->setAttributeValue( "spellcheck", "off" );
    edit->setValidator( std::make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceUnitOptionalRegex ) );
    edit->setTextSize( 8 );
    return edit;
  };

  //Shape
  {
    WTableCell *cell = add_row( WString::tr("dgi-shape") );
    m_shape = cell->addNew<WComboBox>();
    m_shape->addItem( WString::tr("dgi-shape-cyl") );
    m_shape->addItem( WString::tr("dgi-shape-coax") );
    m_shape->addItem( WString::tr("dgi-shape-box") );
    m_shape->setCurrentIndex( 0 );
    m_shape->activated().connect( this, &DetectorGeometryInput::handleShapeChange );
  }

  //Crystal material
  {
    WTableCell *cell = add_row( WString::tr("dgi-crystal-material") );
    m_crystalMaterial = cell->addNew<WComboBox>();
    m_crystalMaterial->addItem( "HPGe" );
    m_crystalMaterial->addItem( "NaI(Tl)" );
    m_crystalMaterial->addItem( "CZT" );
    m_crystalMaterial->addItem( "LaBr3" );
    m_crystalMaterial->setCurrentIndex( 1 );
    m_crystalMaterial->activated().connect( this, &DetectorGeometryInput::handleUserInput );
  }

  //Dimensions (labels swap with shape)
  {
    WTableCell *cell = add_row( "" );
    m_dim1Label = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_dim1 = make_dist_edit( cell );

    cell = add_row( "" );
    m_dim2Label = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_dim2 = make_dist_edit( cell );

    cell = add_row( "" );
    m_dim3Label = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_dim3Row = table->rowAt( row-1 );
    m_dim3 = make_dist_edit( cell );

    for( WLineEdit *edit : { m_dim1, m_dim2, m_dim3 } )
    {
      edit->changed().connect( this, &DetectorGeometryInput::handleUserInput );
      edit->enterPressed().connect( this, &DetectorGeometryInput::handleUserInput );
    }
  }

  //Bore (coax only)
  {
    WTableCell *cell = add_row( WString::tr("dgi-bore") );
    m_boreLabel = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_boreRow = table->rowAt( row-1 );
    m_boreDiam = make_dist_edit( cell );
    m_boreDiam->setPlaceholderText( WString::tr("dgi-bore-diam-ph") );
    m_boreDepth = make_dist_edit( cell );
    m_boreDepth->setPlaceholderText( WString::tr("dgi-bore-depth-ph") );
    m_boreDiam->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_boreDepth->changed().connect( this, &DetectorGeometryInput::handleUserInput );
  }

  //Dead layer (coax/planar HPGe-style; optional)
  {
    WTableCell *cell = add_row( WString::tr("dgi-dead-layer") );
    m_deadLabel = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_deadRow = table->rowAt( row-1 );
    m_deadFront = make_dist_edit( cell );
    m_deadFront->setPlaceholderText( WString::tr("dgi-dead-front-ph") );
    m_deadSide = make_dist_edit( cell );
    m_deadSide->setPlaceholderText( WString::tr("dgi-dead-side-ph") );
    m_deadFront->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_deadSide->changed().connect( this, &DetectorGeometryInput::handleUserInput );
  }

  //Concentric layers (endcap/window/housing); multiple allowed.
  {
    WContainerWidget *layersHolder = addNew<WContainerWidget>();
    layersHolder->addStyleClass( "DgiLayers" );
    layersHolder->addNew<WText>( WString::tr("dgi-layers-label") );

    m_layersTable = layersHolder->addNew<WTable>();
    m_layersTable->addStyleClass( "DgiLayersTable" );
    m_layersTable->elementAt( 0, 0 )->addNew<WText>( WString::tr("dgi-layer-material") );
    m_layersTable->elementAt( 0, 1 )->addNew<WText>( WString::tr("dgi-layer-front") );
    m_layersTable->elementAt( 0, 2 )->addNew<WText>( WString::tr("dgi-layer-side") );

    WContainerWidget *btns = layersHolder->addNew<WContainerWidget>();
    m_addLayer = btns->addNew<WPushButton>( WString::tr("dgi-add-layer") );
    m_addLayer->addStyleClass( "LinkBtn" );
    m_addLayer->clicked().connect( this, [this](){ addLayerRow( "", "", "" ); } );
    m_removeLayer = btns->addNew<WPushButton>( WString::tr("dgi-remove-layer") );
    m_removeLayer->addStyleClass( "LinkBtn" );
    m_removeLayer->clicked().connect( this, &DetectorGeometryInput::removeLayerRow );

    addLayerRow( "Al", "1 mm", "1 mm" );
  }

  //Collimator (optional)
  {
    m_hasCollimator = addNew<WCheckBox>( WString::tr("dgi-collimator-cb") );
    m_hasCollimator->changed().connect( this, &DetectorGeometryInput::handleUserInput );

    m_collimatorRow = addNew<WContainerWidget>();
    m_collimatorRow->addStyleClass( "DgiCollimator" );
    m_collimatorRow->addNew<WLabel>( WString::tr("dgi-layer-material") );
    m_collimatorMaterial = m_collimatorRow->addNew<WLineEdit>();
    m_collimatorMaterial->setTextSize( 8 );
    m_materialSuggestion->forEdit( m_collimatorMaterial,
                          PopupTrigger::Editing | PopupTrigger::DropDownIcon );
    m_collimatorRow->addNew<WLabel>( WString::tr("dgi-collimator-thickness") );
    m_collimatorThickness = m_collimatorRow->addNew<WLineEdit>();
    m_collimatorThickness->setTextSize( 6 );
    m_collimatorRow->addNew<WLabel>( WString::tr("dgi-collimator-extension") );
    m_collimatorExtension = m_collimatorRow->addNew<WLineEdit>();
    m_collimatorExtension->setTextSize( 6 );
    for( WLineEdit *edit : { m_collimatorMaterial, m_collimatorThickness, m_collimatorExtension } )
      edit->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_collimatorRow->hide();
  }

  //Reference point
  {
    WContainerWidget *refRow = addNew<WContainerWidget>();
    refRow->addNew<WLabel>( WString::tr("dgi-ref-point") );
    m_referencePoint = refRow->addNew<WComboBox>();
    m_referencePoint->addItem( WString::tr("dgi-ref-endcap") );
    m_referencePoint->addItem( WString::tr("dgi-ref-crystal") );
    m_referencePoint->setCurrentIndex( 0 );  //InterSpec convention: detector face
    m_referencePoint->activated().connect( this, &DetectorGeometryInput::handleUserInput );
  }

  m_note = addNew<WText>( "" );
  m_note->addStyleClass( "DgiNote" );
  m_note->setInline( false );

  handleShapeChange();
}//init()


Wt::Signal<> &DetectorGeometryInput::changed()
{
  return m_changed;
}


void DetectorGeometryInput::handleShapeChange()
{
  const int shape = m_shape->currentIndex();
  const bool coax = (shape == 1);
  const bool box = (shape == 2);

  m_dim1Label->setText( WString::tr( box ? "dgi-dim-width" : "dgi-dim-diameter" ) );
  m_dim2Label->setText( WString::tr( box ? "dgi-dim-height" : "dgi-dim-length" ) );
  m_dim3Label->setText( WString::tr("dgi-dim-length") );

  //Show/hide whole table rows (WTableRow::setHidden) - NOT by walking widget
  //  parents, which for a table cell is not a reliable row handle.
  m_dim3Row->setHidden( !box );          //third dimension only for boxes
  m_boreRow->setHidden( !coax );         //bore only for coaxial HPGe
  m_deadRow->setHidden( box );           //dead layer for cylinders/HPGe, not boxes

  handleUserInput();
}//handleShapeChange()


void DetectorGeometryInput::handleUserInput()
{
  //init() wires this handler into widgets it creates along the way (e.g. the initial
  //  layer row), so it can fire before later widgets exist; m_note is created last.
  if( !m_note )
    return;

  m_collimatorRow->setHidden( !m_hasCollimator->isChecked() );

  try
  {
    toDescriptor();
    m_note->setText( "" );
  }catch( std::exception &e )
  {
    m_note->setText( WString::fromUTF8( e.what() ) );
  }

  m_changed.emit();
}//handleUserInput()


void DetectorGeometryInput::addLayerRow( const Wt::WString &material,
                                         const Wt::WString &frontThick,
                                         const Wt::WString &sideThick )
{
  const int row = m_layersTable->rowCount();

  LayerRow layer;
  layer.material = m_layersTable->elementAt( row, 0 )->addNew<WLineEdit>( material );
  layer.material->setTextSize( 10 );
  m_materialSuggestion->forEdit( layer.material,
                        PopupTrigger::Editing | PopupTrigger::DropDownIcon );
  layer.frontThickness = m_layersTable->elementAt( row, 1 )->addNew<WLineEdit>( frontThick );
  layer.frontThickness->setTextSize( 6 );
  layer.sideThickness = m_layersTable->elementAt( row, 2 )->addNew<WLineEdit>( sideThick );
  layer.sideThickness->setTextSize( 6 );

  for( WLineEdit *edit : { layer.material, layer.frontThickness, layer.sideThickness } )
    edit->changed().connect( this, &DetectorGeometryInput::handleUserInput );

  m_layers.push_back( layer );
  m_removeLayer->setEnabled( m_layers.size() > 1 );

  handleUserInput();
}//addLayerRow(...)


void DetectorGeometryInput::removeLayerRow()
{
  if( m_layers.size() <= 1 )
    return;

  m_layersTable->removeRow( m_layersTable->rowCount() - 1 );
  m_layers.pop_back();
  m_removeLayer->setEnabled( m_layers.size() > 1 );

  handleUserInput();
}//removeLayerRow()


bool DetectorGeometryInput::isValid() const
{
  try
  {
    toDescriptor();
    return true;
  }catch( std::exception & )
  {
    return false;
  }
}//isValid()


ceelo::GeometryDescriptor DetectorGeometryInput::toDescriptor() const
{
  ceelo::GeometryDescriptor gd;

  const int shape = m_shape->currentIndex();
  const bool coax = (shape == 1);
  const bool box = (shape == 2);

  gd.materials.push_back( crystal_material_spec( m_crystalMaterial->currentIndex() ) );
  gd.crystal_material_index = 0;

  double crystal_len = 0.0;
  if( box )
  {
    const double width = distance_cm( m_dim1, "width", false );
    const double height = distance_cm( m_dim2, "height", false );
    crystal_len = distance_cm( m_dim3, "length", false );
    gd.shape = ceelo::DetectorShape::Box;
    gd.dimensions_cm = { 0.5*width, 0.5*height, crystal_len };
    gd.symmetry = ceelo::ResponseSymmetry::Quadrant;
  }else
  {
    const double diam = distance_cm( m_dim1, "diameter", false );
    crystal_len = distance_cm( m_dim2, "length", false );
    gd.shape = ceelo::DetectorShape::Cylinder;
    gd.dimensions_cm = { 0.5*diam, crystal_len };
    gd.symmetry = ceelo::ResponseSymmetry::Axial;

    if( coax )
    {
      const double bore_diam = distance_cm( m_boreDiam, "bore diameter", true );
      const double bore_depth = distance_cm( m_boreDepth, "bore depth", true );
      if( (bore_diam > 0.0) != (bore_depth > 0.0) )
        throw runtime_error( WString::tr("dgi-err-bore").toUTF8() );
      if( bore_diam > 0.0 )
      {
        if( (0.5*bore_diam >= 0.5*diam) || (bore_depth >= crystal_len) )
          throw runtime_error( WString::tr("dgi-err-bore-size").toUTF8() );
        gd.bore = ceelo::BoreHoleConfig{ 0.5*bore_diam, bore_depth };
      }
    }//if( coax )
  }//if( box ) / else

  if( !box )
  {
    const double dead_front = distance_cm( m_deadFront, "front dead layer", true );
    const double dead_side = distance_cm( m_deadSide, "side dead layer", true );
    if( (dead_front > 0.0) || (dead_side > 0.0) )
      gd.dead_layer = ceelo::DeadLayerConfig{ dead_front, dead_side, 0.0 };
  }

  const shared_ptr<const MaterialDB> matDb = MaterialDB::initialized()
                                              ? MaterialDB::instance() : nullptr;

  for( const LayerRow &layer : m_layers )
  {
    const string mat_name = layer.material->text().toUTF8();
    const double front = distance_cm( layer.frontThickness, "layer front thickness", true );
    const double side = distance_cm( layer.sideThickness, "layer side thickness", true );

    if( mat_name.empty() && (front <= 0.0) && (side <= 0.0) )
      continue;  //blank row

    if( mat_name.empty() || ((front <= 0.0) && (side <= 0.0)) )
      throw runtime_error( WString::tr("dgi-err-layer").toUTF8() );

    shared_ptr<const Material> mat = matDb ? matDb->material( mat_name ) : nullptr;
    if( !mat && matDb )
    {
      try
      {
        // Allows chemical formulas like "C0.5H0.5 d=1.2"
        const SandiaDecay::SandiaDecayDataBase * const nucDb = DecayDataBaseServer::database();
        mat = matDb->materialFromChemicalFormula( mat_name, nucDb );
      }catch( std::exception & )
      {
      }
    }
    if( !mat )
      throw runtime_error( WString::tr("dgi-err-material").arg(mat_name).toUTF8() );

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( CeeLoUtils::to_ceelo_material(*mat) );
    spec.front_thickness_cm = front;
    spec.side_thickness_cm = side;
    spec.z_start_cm = 0.0;
    spec.z_end_cm = crystal_len;
    gd.layers.push_back( spec );
  }//for( const LayerRow &layer : m_layers )

  if( m_hasCollimator->isChecked() )
  {
    const string mat_name = m_collimatorMaterial->text().toUTF8();
    const double thickness = distance_cm( m_collimatorThickness, "collimator thickness", false );
    const double extension = distance_cm( m_collimatorExtension, "collimator extension", true );

    shared_ptr<const Material> mat = matDb ? matDb->material( mat_name ) : nullptr;
    if( !mat )
      throw runtime_error( WString::tr("dgi-err-material").arg(mat_name).toUTF8() );

    ceelo::CollimatorSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( CeeLoUtils::to_ceelo_material(*mat) );
    spec.side_thickness_cm = thickness;
    spec.z_start_cm = -extension;
    spec.z_end_cm = crystal_len;
    gd.collimator = spec;
  }//if( m_hasCollimator->isChecked() )

  gd.reference_point = (m_referencePoint->currentIndex() == 1)
                          ? ceelo::ReferencePoint::CrystalFace
                          : ceelo::ReferencePoint::EndcapFront;

  return gd;
}//toDescriptor()


void DetectorGeometryInput::setFromDescriptor( const ceelo::GeometryDescriptor &gd )
{
  const bool box = (gd.shape == ceelo::DetectorShape::Box);
  const bool coax = !box && gd.bore.has_value();
  m_shape->setCurrentIndex( box ? 2 : (coax ? 1 : 0) );

  if( gd.crystal_material_index >= 0
      && (gd.crystal_material_index < static_cast<int>(gd.materials.size())) )
  {
    const int idx = crystal_index_for_name(
                gd.materials[static_cast<size_t>(gd.crystal_material_index)].name );
    if( idx >= 0 )
      m_crystalMaterial->setCurrentIndex( idx );
  }

  if( box && (gd.dimensions_cm.size() >= 3) )
  {
    m_dim1->setText( cm_to_str( 2.0*gd.dimensions_cm[0] ) );
    m_dim2->setText( cm_to_str( 2.0*gd.dimensions_cm[1] ) );
    m_dim3->setText( cm_to_str( gd.dimensions_cm[2] ) );
  }else if( gd.dimensions_cm.size() >= 2 )
  {
    m_dim1->setText( cm_to_str( 2.0*gd.dimensions_cm[0] ) );
    m_dim2->setText( cm_to_str( gd.dimensions_cm[1] ) );
    m_dim3->setText( "" );
  }

  if( gd.bore )
  {
    m_boreDiam->setText( cm_to_str( 2.0*gd.bore->radius ) );
    m_boreDepth->setText( cm_to_str( gd.bore->depth ) );
  }else
  {
    m_boreDiam->setText( "" );
    m_boreDepth->setText( "" );
  }

  if( gd.dead_layer )
  {
    m_deadFront->setText( cm_to_str( gd.dead_layer->front ) );
    m_deadSide->setText( cm_to_str( gd.dead_layer->side ) );
  }else
  {
    m_deadFront->setText( "" );
    m_deadSide->setText( "" );
  }

  while( m_layers.size() > 1 )
    removeLayerRow();

  bool first = true;
  for( const ceelo::LayerSpec &layer : gd.layers )
  {
    if( (layer.material_index < 0)
        || (layer.material_index >= static_cast<int>(gd.materials.size())) )
      continue;
    const string name = gd.materials[static_cast<size_t>(layer.material_index)].name;
    const string front = cm_to_str( layer.front_thickness_cm );
    const string side = cm_to_str( layer.side_thickness_cm );
    if( first )
    {
      m_layers[0].material->setText( WString::fromUTF8(name) );
      m_layers[0].frontThickness->setText( WString::fromUTF8(front) );
      m_layers[0].sideThickness->setText( WString::fromUTF8(side) );
      first = false;
    }else
    {
      addLayerRow( WString::fromUTF8(name), WString::fromUTF8(front),
                   WString::fromUTF8(side) );
    }
  }//for( const ceelo::LayerSpec &layer : gd.layers )

  if( first )  //no layers
  {
    m_layers[0].material->setText( "" );
    m_layers[0].frontThickness->setText( "" );
    m_layers[0].sideThickness->setText( "" );
  }

  m_hasCollimator->setChecked( gd.collimator.has_value() );
  if( gd.collimator
      && (gd.collimator->material_index >= 0)
      && (gd.collimator->material_index < static_cast<int>(gd.materials.size())) )
  {
    m_collimatorMaterial->setText( WString::fromUTF8(
              gd.materials[static_cast<size_t>(gd.collimator->material_index)].name ) );
    m_collimatorThickness->setText( cm_to_str( gd.collimator->side_thickness_cm ) );
    m_collimatorExtension->setText( cm_to_str( std::max(0.0, -gd.collimator->z_start_cm) ) );
  }

  m_referencePoint->setCurrentIndex(
      (gd.reference_point == ceelo::ReferencePoint::CrystalFace) ? 1 : 0 );

  handleShapeChange();
}//setFromDescriptor(...)


void DetectorGeometryInput::seedFromDrf( std::shared_ptr<const DetectorPeakResponse> drf )
{
  if( !drf || !drf->isValid() || (drf->detectorDiameter() <= 0.0f) )
    return;

  if( drf->ceeloResponse() )
  {
    setFromDescriptor( drf->ceeloResponse()->descriptor );
    return;
  }

  const double diam_cm = drf->detectorDiameter() / PhysicalUnits::cm;
  m_shape->setCurrentIndex( 0 );
  m_dim1->setText( cm_to_str( diam_cm ) );
  m_dim2->setText( cm_to_str( diam_cm ) );  //length unknown: guess = diameter
  m_note->setText( WString::tr("dgi-seeded-note") );
  handleShapeChange();
}//seedFromDrf(...)
