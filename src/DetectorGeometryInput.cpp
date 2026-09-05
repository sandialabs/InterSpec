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
#include <cctype>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <Wt/WText.h>
#include <Wt/WTable.h>
#include <Wt/WTableRow.h>
#include <Wt/WTableCell.h>
#include <Wt/WLabel.h>
#include <Wt/WCheckBox.h>
#include <Wt/Utils.h>
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
#include "InterSpec/GadrasDetectorDat.h"
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





  /** The crystal-material combo entries, taken from the GADRAS material table.

   Driven by the table rather than hand-listed so the two cannot drift: any
   crystal a `Detector.dat` can name is a crystal this form offers, so importing
   one never has to substitute something else for it.
   `CeeLoUtils::gadrasCrystalMaterial` resolves each name - preferring a CeeLo
   built-in where one exists, else building it from the table's formula and
   density - which is the same path the import itself takes, so an imported
   descriptor matches an entry here by name.
   */
  const std::vector<std::string> &crystal_names()
  {
    static const std::vector<std::string> names = []() -> std::vector<std::string> {
      std::vector<std::string> n;
      for( const GadrasDetectorDat::MaterialInfo &info : GadrasDetectorDat::materialTable() )
      {
        // Only what can actually be built; a table entry whose formula does not
        //  resolve would throw the moment it was selected.
        try
        {
          CeeLoUtils::gadrasCrystalMaterial( info.name );
          n.push_back( info.name );
        }catch( std::exception & )
        {
        }
      }
      return n;
    }();

    return names;
  }//crystal_names()


  ceelo::MaterialSpec crystal_material_spec( const int index )
  {
    const std::vector<std::string> &names = crystal_names();
    if( (index >= 0) && (index < static_cast<int>(names.size())) )
      return CeeLoUtils::gadrasCrystalMaterial( names[static_cast<size_t>(index)] );

    throw runtime_error( "Please select a crystal material." );
  }//crystal_material_spec(...)


  /** The crystal a DRFs free text names, as an index into the crystal combo, or -1.

   A DRF seeded from a GADRAS `Detector.dat` records the crystal in its description ("GADRAS
   geometry, LaBr3; ..."), and shipped detectors carry it in their name ("LaBr 10%").  That text is
   the only place the crystal survives for a DRF that reaches this form without a
   `ceelo::GeometryDescriptor` - which is every DRF that was handed to the app as an efficiency
   curve rather than characterized here.  Guessing from it beats silently rebuilding a LaBr detector
   out of NaI.

   Matched on the GADRAS material-table names, longest first, and only on whole tokens - so "Si"
   does not match inside "Since", and "LaBr3" wins over any shorter prefix.
   */
  int crystal_index_from_text( const std::string &text )
  {
    if( text.empty() )
      return -1;

    std::string haystack = text;
    SpecUtils::to_lower_ascii( haystack );

    // Longest name first, so a name that contains a shorter one cannot shadow it.
    std::vector<std::string> names = crystal_names();
    std::sort( begin(names), end(names),
      []( const std::string &lhs, const std::string &rhs ){ return lhs.size() > rhs.size(); } );

    for( const std::string &name : names )
    {
      std::string needle = name;
      SpecUtils::to_lower_ascii( needle );

      for( size_t pos = haystack.find(needle);
           pos != std::string::npos;
           pos = haystack.find(needle, pos + 1) )
      {
        const bool start_ok = (pos == 0) || !std::isalnum( static_cast<unsigned char>(haystack[pos-1]) );
        const size_t after = pos + needle.size();
        const bool end_ok = (after >= haystack.size())
                            || !std::isalnum( static_cast<unsigned char>(haystack[after]) );
        if( !start_ok || !end_ok )
          continue;

        // The combo is populated from `crystal_names()` in order, so its index is the entrys.
        const std::vector<std::string> &all = crystal_names();
        for( int i = 0; i < static_cast<int>(all.size()); ++i )
        {
          if( all[i] == name )
            return i;
        }
      }//for( each occurrence )
    }//for( const std::string &name : names )

    return -1;
  }//crystal_index_from_text(...)


  /** Exact match of a stored MaterialSpec back to a crystal combo entry, or -1. */
  int crystal_index_for_name( const std::string &name )
  {
    for( int i = 0; i < static_cast<int>(crystal_names().size()); ++i )
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
    m_bulletLabel( nullptr ),
    m_bulletRadius( nullptr ),
    m_bulletRow( nullptr ),
    m_boreLabel( nullptr ),
    m_boreDiam( nullptr ), m_boreDepth( nullptr ),
    m_boreRounded( nullptr ),
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
    m_note( nullptr ),
    m_importNotes( nullptr ),
    m_materialSuggestion( nullptr ),
    m_restoringState( false ),
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
    for( const std::string &name : crystal_names() )
      m_crystalMaterial->addItem( WString::fromUTF8(name) );

    // Default to NaI, the most common crystal, when the table offers it.
    const int nai = crystal_index_for_name( "NaI" );
    m_crystalMaterial->setCurrentIndex( (nai >= 0) ? nai : 0 );
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

  //Front-edge fillet ("bulletization"); any cylindrical crystal, coax or not
  {
    WTableCell *cell = add_row( WString::tr("dgi-bullet") );
    m_bulletLabel = dynamic_cast<WLabel *>( table->elementAt(row-1,0)->widget(0) );
    m_bulletRow = table->rowAt( row-1 );
    m_bulletRadius = make_dist_edit( cell );
    m_bulletRadius->setPlaceholderText( WString::tr("dgi-bullet-ph") );
    m_bulletRadius->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_bulletRadius->enterPressed().connect( this, &DetectorGeometryInput::handleUserInput );
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
    m_boreRounded = cell->addNew<WCheckBox>( WString::tr("dgi-bore-rounded") );
    m_boreRounded->addStyleClass( "DgiBoreRounded" );
    m_boreDiam->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_boreDepth->changed().connect( this, &DetectorGeometryInput::handleUserInput );
    m_boreRounded->changed().connect( this, &DetectorGeometryInput::handleUserInput );
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

  m_note = addNew<WText>( "" );
  m_note->addStyleClass( "DgiNote" );
  m_note->setInline( false );

  // Import warnings - what the source file could not express.  Sits with the
  //  form the user is being asked to correct, rather than only in a toast that
  //  is gone by the time they look at the geometry.
  m_importNotes = addNew<WText>( "" );
  m_importNotes->addStyleClass( "DgiImportNotes" );
  m_importNotes->setInline( false );
  m_importNotes->hide();

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
  m_bulletRow->setHidden( box );         //a fillet is a cylinder-crystal feature
  m_boreRow->setHidden( !coax );         //bore (and its rounded tip) only for coaxial HPGe
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

  if( !m_restoringState )
    m_changed.emit();
}//handleUserInput()


void DetectorGeometryInput::addLayerRow( const Wt::WString &material,
                                         const Wt::WString &frontThick,
                                         const Wt::WString &sideThick,
                                         const shared_ptr<const ceelo::MaterialSpec> &seeded )
{
  const int row = m_layersTable->rowCount();

  LayerRow layer;
  layer.seeded = seeded;
  layer.seededName = material.toUTF8();
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


bool DetectorGeometryInput::State::Layer::operator==( const Layer &rhs ) const
{
  return (material == rhs.material)
         && (frontThickness == rhs.frontThickness)
         && (sideThickness == rhs.sideThickness)
         && (seeded == rhs.seeded)
         && (seededName == rhs.seededName);
}//State::Layer::operator==


bool DetectorGeometryInput::State::operator==( const State &rhs ) const
{
  return (shape == rhs.shape)
         && (crystalMaterial == rhs.crystalMaterial)
         && (dim1 == rhs.dim1) && (dim2 == rhs.dim2) && (dim3 == rhs.dim3)
         && (bulletRadius == rhs.bulletRadius)
         && (boreDiam == rhs.boreDiam) && (boreDepth == rhs.boreDepth)
         && (deadFront == rhs.deadFront) && (deadSide == rhs.deadSide)
         && (boreRounded == rhs.boreRounded)
         && (hasCollimator == rhs.hasCollimator)
         && (collimatorMaterial == rhs.collimatorMaterial)
         && (collimatorThickness == rhs.collimatorThickness)
         && (collimatorExtension == rhs.collimatorExtension)
         && (layers == rhs.layers);
}//State::operator==


DetectorGeometryInput::State DetectorGeometryInput::currentState() const
{
  State state;

  state.shape = m_shape->currentIndex();
  state.crystalMaterial = m_crystalMaterial->currentIndex();

  state.dim1 = m_dim1->text().toUTF8();
  state.dim2 = m_dim2->text().toUTF8();
  state.dim3 = m_dim3->text().toUTF8();
  state.bulletRadius = m_bulletRadius->text().toUTF8();
  state.boreDiam = m_boreDiam->text().toUTF8();
  state.boreDepth = m_boreDepth->text().toUTF8();
  state.deadFront = m_deadFront->text().toUTF8();
  state.deadSide = m_deadSide->text().toUTF8();

  state.boreRounded = m_boreRounded->isChecked();
  state.hasCollimator = m_hasCollimator->isChecked();
  state.collimatorMaterial = m_collimatorMaterial->text().toUTF8();
  state.collimatorThickness = m_collimatorThickness->text().toUTF8();
  state.collimatorExtension = m_collimatorExtension->text().toUTF8();

  for( const LayerRow &row : m_layers )
  {
    State::Layer layer;
    layer.material = row.material->text().toUTF8();
    layer.frontThickness = row.frontThickness->text().toUTF8();
    layer.sideThickness = row.sideThickness->text().toUTF8();
    layer.seeded = row.seeded;
    layer.seededName = row.seededName;
    state.layers.push_back( std::move(layer) );
  }//for( const LayerRow &row : m_layers )

  return state;
}//DetectorGeometryInput::currentState()


void DetectorGeometryInput::setState( const State &state )
{
  m_restoringState = true;

  m_shape->setCurrentIndex( state.shape );
  m_crystalMaterial->setCurrentIndex( state.crystalMaterial );

  m_dim1->setText( WString::fromUTF8(state.dim1) );
  m_dim2->setText( WString::fromUTF8(state.dim2) );
  m_dim3->setText( WString::fromUTF8(state.dim3) );
  m_bulletRadius->setText( WString::fromUTF8(state.bulletRadius) );
  m_boreDiam->setText( WString::fromUTF8(state.boreDiam) );
  m_boreDepth->setText( WString::fromUTF8(state.boreDepth) );
  m_deadFront->setText( WString::fromUTF8(state.deadFront) );
  m_deadSide->setText( WString::fromUTF8(state.deadSide) );

  m_boreRounded->setChecked( state.boreRounded );
  m_hasCollimator->setChecked( state.hasCollimator );
  m_collimatorMaterial->setText( WString::fromUTF8(state.collimatorMaterial) );
  m_collimatorThickness->setText( WString::fromUTF8(state.collimatorThickness) );
  m_collimatorExtension->setText( WString::fromUTF8(state.collimatorExtension) );

  // Rebuild the layer rows; `addLayerRow` is the only place that wires their signals up.
  while( m_layersTable->rowCount() > 1 )
    m_layersTable->removeRow( m_layersTable->rowCount() - 1 );
  m_layers.clear();
  for( const State::Layer &layer : state.layers )
    addLayerRow( WString::fromUTF8(layer.material), WString::fromUTF8(layer.frontThickness),
                 WString::fromUTF8(layer.sideThickness), layer.seeded );
  for( size_t i = 0; (i < m_layers.size()) && (i < state.layers.size()); ++i )
    m_layers[i].seededName = state.layers[i].seededName;
  m_removeLayer->setEnabled( m_layers.size() > 1 );

  // Sync the shape-dependent row visibility (and, through `handleUserInput`, the validity note and
  //  collimator row) to the restored values.
  handleShapeChange();

  m_restoringState = false;
}//DetectorGeometryInput::setState(...)


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

    //Read before the bore, so the bore checks below see the final crystal profile
    gd.bullet_radius_cm = distance_cm( m_bulletRadius, "front edge fillet radius", true );

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
        gd.bore = ceelo::BoreHoleConfig{ 0.5*bore_diam, bore_depth,
                                         m_boreRounded->isChecked() };
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

  //Everything CeeLo asserts on (fillet vs radius/length/dead layer, bore vs
  //  fillet and vs the active radius), checked once from the library's own
  //  list, so this form and the ANGLE import can't disagree about what's legal.
  //  Those preconditions are asserts, so without this a release build would
  //  silently trace garbage instead of telling the user what's wrong.
  {
    const vector<ceelo::GeometryProblem> problems = gd.problems();
    if( !problems.empty() )
      throw runtime_error(
          WString::tr( CeeLoUtils::geometryProblemMsgId( problems.front() ) ).toUTF8() );
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

    // A layer seeded from an imported descriptor may name a material this
    //  session's MaterialDB has never heard of (an ANGLE file defines its own
    //  materials inline, for instance).  As long as the user has not changed the
    //  name, the composition it came in with is the right answer.
    const bool use_seeded = layer.seeded && (mat_name == layer.seededName);

    // A generic attenuator - "AN=13.2, AD=1.35 g/cm2" - is not a material name
    //  and never will resolve as one; it is how a GADRAS import expresses a
    //  layer the file gives only an effective atomic number and an areal density
    //  for.  Rebuilding it from the edited text is what lets the user change
    //  either number and have it mean something.
    double gen_an = 0.0, gen_ad = 0.0;
    if( !use_seeded
        && CeeLoUtils::parseGenericAttenuatorName( mat_name, gen_an, gen_ad ) )
    {
      const double thickness = (front > 0.0) ? front : side;
      if( thickness <= 0.0 )
        throw runtime_error( WString::tr("dgi-err-layer").toUTF8() );

      ceelo::LayerSpec gspec;
      gspec.material_index = static_cast<int>( gd.materials.size() );
      gd.materials.push_back(
            CeeLoUtils::genericAttenuatorMaterial( gen_an, gen_ad, thickness ) );
      gspec.front_thickness_cm = front;
      gspec.side_thickness_cm = side;
      gspec.z_start_cm = 0.0;
      gspec.z_end_cm = crystal_len;
      gd.layers.push_back( gspec );
      continue;
    }//if( a generic AN/AD attenuator )

    shared_ptr<const Material> mat;
    if( !use_seeded )
    {
      mat = matDb ? matDb->material( mat_name ) : nullptr;
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
    }//if( !use_seeded )

    ceelo::LayerSpec spec;
    spec.material_index = static_cast<int>( gd.materials.size() );
    gd.materials.push_back( use_seeded ? *layer.seeded
                                       : CeeLoUtils::to_ceelo_material(*mat) );
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

  // InterSpec measures every distance from the detector face, and forms CeeLo query positions
  //  itself (CeeLoUtils::sourcePositionFromFace), so the descriptor's reference point is never
  //  consulted by InterSpec; it is set to the face for anything else that reads the descriptor.
  gd.reference_point = ceelo::ReferencePoint::EndcapFront;

  return gd;
}//toDescriptor()


void DetectorGeometryInput::setFromDescriptor( const ceelo::GeometryDescriptor &gd,
                                               const std::vector<std::string> &notes )
{
  const bool box = (gd.shape == ceelo::DetectorShape::Box);
  const bool coax = !box && gd.bore.has_value();
  m_shape->setCurrentIndex( box ? 2 : (coax ? 1 : 0) );

  m_substitutedCrystal.clear();
  if( gd.crystal_material_index >= 0
      && (gd.crystal_material_index < static_cast<int>(gd.materials.size())) )
  {
    const ceelo::MaterialSpec &crystal
                  = gd.materials[static_cast<size_t>(gd.crystal_material_index)];
    const int idx = crystal_index_for_name( crystal.name );
    if( idx >= 0 )
    {
      m_crystalMaterial->setCurrentIndex( idx );
    }else
    {
      // Nothing in the list is this material.  The list covers every crystal
      //  GADRAS knows, so in practice this is only reachable from an ANGLE file,
      //  which defines its materials inline and so can name anything at all.
      //  Leaving the combo where it was would quietly rebuild the crystal as
      //  whatever it happened to show, so fall back to NaI and say so.
      const int nai = crystal_index_for_name( "NaI" );
      m_crystalMaterial->setCurrentIndex( (nai >= 0) ? nai : 0 );
      m_substitutedCrystal = crystal.name;
    }
  }//if( the descriptor names a crystal )

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

  m_bulletRadius->setText( (gd.bullet_radius_cm > 0.0)
                           ? cm_to_str( gd.bullet_radius_cm ) : string() );

  if( gd.bore )
  {
    m_boreDiam->setText( cm_to_str( 2.0*gd.bore->radius ) );
    m_boreDepth->setText( cm_to_str( gd.bore->depth ) );
    m_boreRounded->setChecked( gd.bore->rounded_tip );
  }else
  {
    m_boreDiam->setText( "" );
    m_boreDepth->setText( "" );
    m_boreRounded->setChecked( false );
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
    const ceelo::MaterialSpec &mspec = gd.materials[static_cast<size_t>(layer.material_index)];
    const string name = mspec.name;
    const string front = cm_to_str( layer.front_thickness_cm );
    const string side = cm_to_str( layer.side_thickness_cm );
    const shared_ptr<const ceelo::MaterialSpec> seeded
                            = make_shared<const ceelo::MaterialSpec>( mspec );
    if( first )
    {
      m_layers[0].seeded = seeded;
      m_layers[0].seededName = name;
      m_layers[0].material->setText( WString::fromUTF8(name) );
      m_layers[0].frontThickness->setText( WString::fromUTF8(front) );
      m_layers[0].sideThickness->setText( WString::fromUTF8(side) );
      first = false;
    }else
    {
      addLayerRow( WString::fromUTF8(name), WString::fromUTF8(front),
                   WString::fromUTF8(side), seeded );
    }
  }//for( const ceelo::LayerSpec &layer : gd.layers )

  if( first )  //no layers
  {
    m_layers[0].seeded.reset();
    m_layers[0].seededName.clear();
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

  // Import warnings, plus the crystal substitution if there was one - the
  //  substitution changes what gets simulated, so it belongs with them.
  vector<string> all_notes = notes;
  if( !m_substitutedCrystal.empty() )
  {
    const int idx = m_crystalMaterial->currentIndex();
    const string used = ((idx >= 0) && (idx < static_cast<int>(crystal_names().size())))
                        ? crystal_names()[static_cast<size_t>(idx)] : string();
    all_notes.push_back( "This detector's crystal is " + m_substitutedCrystal
                         + ", which is not a material this form offers, so " + used
                         + " was selected instead.  Its efficiency will differ -"
                         " please set the crystal material yourself." );
  }//if( a crystal was substituted )

  if( all_notes.empty() )
  {
    m_importNotes->setText( "" );
    m_importNotes->hide();
  }else
  {
    string html;
    for( const string &n : all_notes )
      html += "<div>" + Wt::Utils::htmlEncode(n) + "</div>";
    m_importNotes->setText( WString::fromUTF8(html) );
    m_importNotes->show();
  }

  handleShapeChange();
}//setFromDescriptor(...)


void DetectorGeometryInput::seedFromDrf( std::shared_ptr<const DetectorPeakResponse> drf )
{
  if( !drf || !drf->isValid() || (drf->detectorDiameter() <= 0.0f) )
    return;

  // The DRFs own geometry, when it has one - from a generated response, or set by an importer
  //  (a GADRAS Detector.dat, an ANGLE model) that knew the detector's shape.
  if( const std::shared_ptr<const ceelo::GeometryDescriptor> gd = drf->geometry() )
  {
    setFromDescriptor( *gd );
    return;
  }

  const double diam_cm = drf->detectorDiameter() / PhysicalUnits::cm;
  m_shape->setCurrentIndex( 0 );
  m_dim1->setText( cm_to_str( diam_cm ) );
  m_dim2->setText( cm_to_str( diam_cm ) );  //length unknown: guess = diameter

  // No geometry at all: the description (and often the name) still says what the crystal is, and
  //  taking the default NaI for, say, a LaBr detector would simulate the wrong material entirely.
  const int crystal = crystal_index_from_text( drf->description() + " " + drf->name() );
  if( crystal >= 0 )
    m_crystalMaterial->setCurrentIndex( crystal );

  m_note->setText( WString::tr("dgi-seeded-note") );
  handleShapeChange();
}//seedFromDrf(...)
