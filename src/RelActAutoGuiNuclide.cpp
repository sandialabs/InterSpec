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
#include <limits>
#include <random>

#include <Wt/WLabel>
#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/RelActAutoGui.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/RelActAutoGuiNuclide.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"


using namespace std;
using namespace Wt;

namespace
{
  const WColor ns_default_color = WColor("#FF6633");
}//namespace

enum class NucConstraintType : int
{
  None, RelActRange, MassFraction, ActRatio, NumTypes
};//enum class NucConstraintType

class RelActAutoGuiNuclideConstraint : public WContainerWidget
{
  RelActAutoGui *m_gui;
  RelActAutoGuiNuclide *m_nuc;
  Wt::WComboBox *m_constraint_type;
  Wt::WStackedWidget *m_stacked_widget;

  NativeFloatSpinBox *m_min_rel_act_edit;
  NativeFloatSpinBox *m_max_rel_act_edit;

  NativeFloatSpinBox *m_min_mass_frac_edit;
  NativeFloatSpinBox *m_max_mass_frac_edit;

  Wt::WComboBox *m_act_control_nuc_combo;
  NativeFloatSpinBox *m_activity_ratio_edit;

  Wt::Signal<> m_remove;
  Wt::Signal<> m_changed;
  

public:
  RelActAutoGuiNuclideConstraint( RelActAutoGui *gui, RelActAutoGuiNuclide *nuc, WContainerWidget *parent )
    : WContainerWidget( parent ),
    m_gui( gui ),
    m_nuc( nuc ),
    m_constraint_type( nullptr ),
    m_stacked_widget( nullptr ),
    m_min_rel_act_edit( nullptr ),
    m_max_rel_act_edit( nullptr ),
    m_min_mass_frac_edit( nullptr ),
    m_max_mass_frac_edit( nullptr ),
    m_act_control_nuc_combo( nullptr ),
    m_activity_ratio_edit( nullptr ),
    m_remove( this ),
    m_changed( this )
  {
    addStyleClass( "RelActAutoGuiNuclideConstraint" );

    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );

    m_constraint_type = new WComboBox( this );
    m_constraint_type->setNoSelectionEnabled( false );
    m_constraint_type->changed().connect( this, &RelActAutoGuiNuclideConstraint::handleConstraintTypeChange );
    m_stacked_widget = new WStackedWidget( this );

   
    for( NucConstraintType type = NucConstraintType::None; 
        type < NucConstraintType::NumTypes; 
        type = NucConstraintType(static_cast<int>(type) + 1) )
    {
      WContainerWidget *container = new WContainerWidget( m_stacked_widget );
      container->addStyleClass( "ConstraintInputDiv" );
      m_stacked_widget->addWidget( container );
      const int model_row = m_constraint_type->count();

      switch( type )
      {
        case NucConstraintType::None:
        {
          //m_constraint_type->addItem( "None" );
          WLabel *label = new WLabel( "Select a constraint type.", container );
          label->setInline( false );
          break;
        }
        
        case NucConstraintType::RelActRange:
        {
          //m_constraint_type->addItem( "Rel. Act" );
          WLabel *label = new WLabel( "min:", container );
          m_min_rel_act_edit = new NativeFloatSpinBox( container );
          label->setBuddy( m_min_rel_act_edit );
          m_min_rel_act_edit->setSpinnerHidden( true );
          m_min_rel_act_edit->setMinimum( std::numeric_limits<float>::min() );
          m_min_rel_act_edit->setWidth( WLength(35.0, WLength::Pixel) );
          m_min_rel_act_edit->valueChanged().connect( this, &RelActAutoGuiNuclideConstraint::handleRelActRangeChange );

          label = new WLabel( "max:", container );
          m_max_rel_act_edit = new NativeFloatSpinBox( container );
          label->setBuddy( m_max_rel_act_edit );
          m_max_rel_act_edit->setSpinnerHidden( true );
          m_max_rel_act_edit->setWidth( WLength(35.0, WLength::Pixel) );
          m_max_rel_act_edit->setMinimum( std::numeric_limits<float>::min() );
          m_max_rel_act_edit->valueChanged().connect( this, &RelActAutoGuiNuclideConstraint::handleRelActRangeChange );
          break;
        }//case NucConstraintType::RelActRange
        
        case NucConstraintType::MassFraction:
        {
          //m_constraint_type->addItem( "Mass Frac." );
          WLabel *label = new WLabel( "min:", container );
          m_min_mass_frac_edit = new NativeFloatSpinBox( container );
          label->setBuddy( m_min_mass_frac_edit );
          m_min_mass_frac_edit->setSpinnerHidden( true );
          m_min_mass_frac_edit->setWidth( WLength(35.0, WLength::Pixel) );
          m_min_mass_frac_edit->setRange( 0.0, 1.0 );
          m_min_mass_frac_edit->valueChanged().connect( this, &RelActAutoGuiNuclideConstraint::handleMassFractionChange );

          label = new WLabel( "max:", container );
          m_max_mass_frac_edit = new NativeFloatSpinBox( container );
          label->setBuddy( m_max_mass_frac_edit );
          m_max_mass_frac_edit->setSpinnerHidden( true );
          m_max_mass_frac_edit->setWidth( WLength(35.0, WLength::Pixel) );
          m_max_mass_frac_edit->setRange( 0.0, 1.0 );
          m_max_mass_frac_edit->valueChanged().connect( this, &RelActAutoGuiNuclideConstraint::handleMassFractionChange );
          break;
        }//case NucConstraintType::MassFraction
        
        case NucConstraintType::ActRatio:
        {
          //m_constraint_type->addItem( "Act Ratio" );
          WLabel *label = new WLabel( "src:", container );
          m_act_control_nuc_combo = new WComboBox( container );
          label->setBuddy( m_act_control_nuc_combo );
          m_act_control_nuc_combo->setMaximumSize( WLength(60, WLength::Pixel), WLength::Auto );
          m_act_control_nuc_combo->changed().connect( this, &RelActAutoGuiNuclideConstraint::handleActivityRatioChange );
          
          label = new WLabel( "val:", container );
          m_activity_ratio_edit = new NativeFloatSpinBox( container );
          label->setBuddy( m_activity_ratio_edit );
          m_activity_ratio_edit->setSpinnerHidden( true );
          m_activity_ratio_edit->setWidth( WLength(35.0, WLength::Pixel) );
          m_activity_ratio_edit->setMinimum( std::numeric_limits<float>::min() );
          m_activity_ratio_edit->valueChanged().connect( this, &RelActAutoGuiNuclideConstraint::handleActivityRatioChange );
          break;
        }//case NucConstraintType::ActRatio
        
        case NucConstraintType::NumTypes:
          assert( 0 );
          break;
      }//switch( type )
    }//for( loop over NucConstraintType )


    m_constraint_type->setCurrentIndex( static_cast<int>(NucConstraintType::None) );

    m_constraint_type->changed().connect( this, &RelActAutoGuiNuclideConstraint::handleConstraintTypeChange );

    WContainerWidget *spacer = new WContainerWidget( this );
    spacer->setStyleClass( "RelActAutoSpacer" );

    WPushButton *removeBtn = new WPushButton( this );
    removeBtn->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
    removeBtn->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
    removeBtn->clicked().connect( this, &RelActAutoGuiNuclideConstraint::handleRemove );

    HelpSystem::attachToolTipOn( removeBtn, "Remove this constraint", showToolTips );

    updateAllowedConstraints();
  }//RelActAutoGuiNuclideConstraint


  void setAvailableConstraintTypes( set<NucConstraintType> types )
  {
    if( !types.count(NucConstraintType::None) )
      types.insert( NucConstraintType::None );

    const NucConstraintType current_type = constraintType();

    m_constraint_type->clear();
    WAbstractItemModel *constraint_type_model = m_constraint_type->model();
    assert( constraint_type_model );
    if( !constraint_type_model )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::RelActAutoGuiNuclideConstraint() called when no constraint type model is present" );


    for( const NucConstraintType &type : types )
    {
      const int model_row = m_constraint_type->count();

      const char *title = nullptr;
      switch( type )
      {
        case NucConstraintType::None:
        {
          title = "None";
          break;
        }
        
        case NucConstraintType::RelActRange:
        {
          title = "Rel. Act";
          break;
        }//case NucConstraintType::MassFraction
        
        case NucConstraintType::ActRatio:
        {
          title = "Act Ratio";
          break;
        }//case NucConstraintType::ActRatio
        
        case NucConstraintType::MassFraction:
        {
          title = "Mass Frac.";
          break;
        }//case NucConstraintType::MassFraction
        
        case NucConstraintType::NumTypes:
          assert( 0 );
          break;
      }//switch( type )

      assert( title );
      if( !title )
        throw runtime_error( "RelActAutoGuiNuclideConstraint::setAvailableConstraintTypes() called when no title is available" );

      m_constraint_type->addItem( WString::fromUTF8(title) );
      
      if( type == current_type )
        m_constraint_type->setCurrentIndex( model_row );
      
      const WModelIndex index = m_constraint_type->model()->index( model_row, 0 );

      assert( index.isValid() );
      if( !index.isValid() )
        throw runtime_error( "RelActAutoGuiNuclideConstraint::RelActAutoGuiNuclideConstraint() called when no constraint type model index is valid" );

      constraint_type_model->setData( index, boost::any(type), Wt::ItemDataRole::UserRole );
    }//for( const NucConstraintType &type : types )
  }//void setAvaiableConstraintTypes( const std::set<NucConstraintType> &types )


  void handleConstraintTypeChange()
  {
    m_stacked_widget->setCurrentIndex( static_cast<int>(constraintType()) );

    m_changed.emit();
  }

  void handleRelActRangeChange()
  {
    m_changed.emit();
  }

  void handleMassFractionChange()
  {
    m_changed.emit();
  }

  void handleActivityRatioChange()
  {
    m_changed.emit();
  }

  NucConstraintType constraintType() const
  {
    if( m_constraint_type->count() == 0 )
      return NucConstraintType::None;

    const std::string current_text = m_constraint_type->currentText().toUTF8();

    const int row = m_constraint_type->currentIndex();
    assert( (row >= 0) && (row < static_cast<int>(NucConstraintType::NumTypes)) );
    if( row < 0 )
      return NucConstraintType::None;

    WAbstractItemModel *model = m_constraint_type->model();
    assert( model );
    if( !model )
      return NucConstraintType::None;

    const WModelIndex index = model->index( row, 0 );
    assert( index.isValid() );
    if( !index.isValid() )
      return NucConstraintType::None;

    const boost::any any_data = model->data( index, Wt::ItemDataRole::UserRole);
    assert( !any_data.empty() );
    if( any_data.empty() )
      return NucConstraintType::None;

    const NucConstraintType type = boost::any_cast<NucConstraintType>(any_data);
    return type;
  }//NucConstraintType constraintType() const

  
  Wt::Signal<> &remove()
  {
    return m_remove;
  }

  Wt::Signal<> &changed()
  {
    return m_changed;
  }

  void handleRemove()
  {
    m_remove.emit();
  }


  void setRelActRange( std::optional<double> min_rel_act, std::optional<double> max_rel_act )
  {
    if( min_rel_act.has_value() )
      m_min_rel_act_edit->setText( WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(min_rel_act.value(), 6)) );
    else
      m_min_rel_act_edit->setText( "" );

    if( max_rel_act.has_value() )
      m_max_rel_act_edit->setText( WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(max_rel_act.value(), 6)) );
    else
      m_max_rel_act_edit->setText( "" );
  }

  std::optional<double> minRelAct() const
  {
    if( constraintType() != NucConstraintType::RelActRange )
      return std::nullopt;

    if( m_min_rel_act_edit->text().empty() )
      return std::nullopt;
    return static_cast<double>(m_min_rel_act_edit->value());
  }

  std::optional<double> maxRelAct() const
  {
    if( constraintType() != NucConstraintType::RelActRange )
      return std::nullopt;

    if( m_max_rel_act_edit->text().empty() )
      return std::nullopt;
    return static_cast<double>(m_max_rel_act_edit->value());
  }

  void setActRatio( const RelActCalcAuto::SrcVariant &controlling_src, double constrained_to_controlled_activity_ratio )
  {
    if( RelActCalcAuto::is_null(controlling_src) )
      throw runtime_error( "Invalid ActRatioConstraint - controlling source must be specified" );

    if( constrained_to_controlled_activity_ratio <= 0.0 )
      throw runtime_error( "Invalid ActRatioConstraint - constrained to controlled activity ratio must be > 0.0" );

    m_constraint_type->setCurrentIndex( static_cast<int>(NucConstraintType::ActRatio) );
    m_stacked_widget->setCurrentIndex( static_cast<int>(NucConstraintType::ActRatio) );

    const int rel_eff_curve_index = m_nuc->relEffCurveIndex(); //May throw exception

    //m_act_control_nuc_combo should already be up to date, but just in case, we'll update it
    const std::vector<RelActCalcAuto::NucInputInfo> nucs = m_gui->getNucInputInfo( rel_eff_curve_index );
    m_act_control_nuc_combo->clear();
    int src_index = -1;
    for( int index = 0; index < static_cast<int>(nucs.size()); ++index )
    {
      const RelActCalcAuto::NucInputInfo &nuc = nucs[index];

      assert( !RelActCalcAuto::is_null(nuc.source) );
      if( RelActCalcAuto::is_null(nuc.source) )
        continue;

      if( nuc.source == controlling_src )
        src_index = index;

      // TODO: there are some restrictions on which sources can be used for an act ratio constraint - so we should check that here
      m_act_control_nuc_combo->addItem( WString::fromUTF8(RelActCalcAuto::to_name(nuc.source)) );
    }//for( const RelActCalcAuto::NucInputInfo &nuc : nucs )

    if( src_index == -1 )
      throw runtime_error( "Invalid ActRatioConstraint - controlling nuclide not found" );

    m_act_control_nuc_combo->setCurrentIndex( src_index ); 
    m_activity_ratio_edit->setValue( constrained_to_controlled_activity_ratio );
  }//setActRatio


  RelActCalcAuto::RelEffCurveInput::ActRatioConstraint actRatioConstraint() const 
  {
    if( constraintType() != NucConstraintType::ActRatio )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::actRatioConstraint() called when no constraint is present" );

    if( m_act_control_nuc_combo->currentIndex() < 0 )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::actRatioConstraint() called when no controlling nuclide is selected" );

    const float ratio = static_cast<float>(m_activity_ratio_edit->value());
    if( ratio <= 0.0 )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::actRatioConstraint() called when ratio is <= 0.0" );

    RelActCalcAuto::RelEffCurveInput::ActRatioConstraint constraint;
    constraint.controlling_source = std::monostate();
    constraint.constrained_source = m_nuc->source();
    constraint.constrained_to_controlled_activity_ratio = ratio;

    const string src_name = m_act_control_nuc_combo->currentText().toUTF8();

    const int rel_eff_curve_index = m_nuc->relEffCurveIndex();  //May throw exception
    const std::vector<RelActCalcAuto::NucInputInfo> nucs = m_gui->getNucInputInfo( rel_eff_curve_index );
    for( const RelActCalcAuto::NucInputInfo &src : nucs )
    {
      const RelActCalcAuto::SrcVariant &src_variant = src.source;
      assert( !RelActCalcAuto::is_null(src_variant) );
      if( RelActCalcAuto::is_null(src_variant) )
        continue;

      const string nuc_name = RelActCalcAuto::to_name(src_variant);
      if( nuc_name == src_name )
      {
        constraint.controlling_source = src_variant;
        break;
      }
    }//for( const RelActCalcAuto::NucInputInfo &src : nucs )
    
    if( RelActCalcAuto::is_null(constraint.controlling_source) )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::actRatioConstraint() called when no controlling nuclide is found" );

    return constraint;
  }//RelActCalcAuto::RelEffCurveInput::ActRatioConstraint actRatioConstraint() const


  void setMassFraction( const double lower_mass_fraction, const double upper_mass_fraction )
  {
    m_constraint_type->setCurrentIndex( static_cast<int>(NucConstraintType::MassFraction) );
    m_stacked_widget->setCurrentIndex( static_cast<int>(NucConstraintType::MassFraction) );
    
    const SandiaDecay::Nuclide *nuc = m_nuc->nuclide();
    assert( nuc );
    if( !nuc )
      throw runtime_error( "Invalid MassFractionConstraint - nuclide must be specified" );

    m_min_mass_frac_edit->setValue( static_cast<float>(lower_mass_fraction) );
    m_max_mass_frac_edit->setValue( static_cast<float>(upper_mass_fraction) );
  }//void setMassFraction( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint )


  RelActCalcAuto::RelEffCurveInput::MassFractionConstraint massFractionConstraint() const
  {
    if( constraintType() != NucConstraintType::MassFraction )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::massFractionConstraint() called when no constraint is present" );


    RelActCalcAuto::RelEffCurveInput::MassFractionConstraint constraint;
    constraint.nuclide = m_nuc->nuclide();
    constraint.lower_mass_fraction = static_cast<float>(m_min_mass_frac_edit->value());
    constraint.upper_mass_fraction = static_cast<float>(m_max_mass_frac_edit->value());

    if( !constraint.nuclide )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::massFractionConstraint() called when no nuclide is specified" );

    if( constraint.lower_mass_fraction > constraint.upper_mass_fraction )
      std::swap( constraint.lower_mass_fraction, constraint.upper_mass_fraction );

    if( constraint.lower_mass_fraction < 0.0 )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::massFractionConstraint() called when lower mass fraction is < 0.0" );

    if( constraint.upper_mass_fraction > 1.0 )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::massFractionConstraint() called when upper mass fraction is <= 0.0" );

    if( (constraint.lower_mass_fraction == constraint.upper_mass_fraction) 
        && ((constraint.lower_mass_fraction == 1.0) || (constraint.lower_mass_fraction == 0.0)) )
      throw runtime_error( "RelActAutoGuiNuclideConstraint::massFractionConstraint() called when lower and upper mass fractions are the same and not 0.0 or 1.0" );

    return constraint;
  }//RelActCalcAuto::RelEffCurveInput::MassFractionConstraint massFractionConstraint() const


  void updateAllowedConstraints()
  {
    int rel_eff_curve_index = -1;

    try
    {
      rel_eff_curve_index = m_nuc->relEffCurveIndex();
    }catch( std::exception &e )
    {
      cerr << "Note from RelActAutoGuiNuclideConstraint: not updating allowed constraints"
      << " due to not being able to determine Rel Eff curve index." << endl;
      assert( 0 );
      return;
    }//try / catch

    const vector<RelActCalcAuto::NucInputInfo> nucs = m_gui->getNucInputInfo( rel_eff_curve_index );
    const RelActCalcAuto::SrcVariant src_variant = m_nuc->source();
    const SandiaDecay::Nuclide * const nuc = m_nuc->nuclide();
    const SandiaDecay::Element * const el = m_nuc->element();
    const ReactionGamma::Reaction * const rctn = m_nuc->reaction();

    const string current_text = m_act_control_nuc_combo->currentText().toUTF8();

    //None, RelActRange, MassFraction, ActRatio, 
    set<NucConstraintType> types;

    if( !RelActCalcAuto::is_null(src_variant) )
    {
      size_t same_an_nucs = 0, num_valid_sources = 0;

      for( const RelActCalcAuto::NucInputInfo &other_src : nucs )
      {
        if( RelActCalcAuto::is_null(other_src.source) )
          continue;
        const SandiaDecay::Nuclide *other_nuc = RelActCalcAuto::nuclide(other_src.source);
        if( nuc && other_nuc && (nuc->atomicNumber == other_nuc->atomicNumber) )
          ++same_an_nucs;

        ++num_valid_sources;
      }//

      if( num_valid_sources > 1 )
      {
        types.insert( NucConstraintType::ActRatio );
        types.insert( NucConstraintType::RelActRange );

        if( m_nuc->nuclide() && (same_an_nucs > 1) )
          types.insert( NucConstraintType::MassFraction );
      }
    }//if( m_nuc->nuclide() || m_nuc->element() || m_nuc->reaction() )
    
    setAvailableConstraintTypes( types );
    
    m_act_control_nuc_combo->clear();
    int src_index = -1;
    for( int index = 0; index < static_cast<int>(nucs.size()); ++index )
    {
      const RelActCalcAuto::NucInputInfo &src = nucs[index];

      if( src.source == src_variant )
        continue;

      if( RelActCalcAuto::is_null(src.source) )
        continue;

      const string src_txt = RelActCalcAuto::to_name(src.source);

      if( src_txt == current_text )
        src_index = index;

      // TODO: there are some restrictions on which sources can be used for an act ratio constraint - so we should check that here
      m_act_control_nuc_combo->addItem( WString::fromUTF8(src_txt) );
    }//for( const RelActCalcAuto::NucInputInfo &nuc : nucs )

    if( src_index != -1 )
      m_act_control_nuc_combo->setCurrentIndex( src_index );
  }//void updateAllowedConstraints()
};//class RelActAutoGuiNuclideConstraint



RelActAutoGuiNuclide::RelActAutoGuiNuclide( RelActAutoGui *gui, WContainerWidget *parent )
  : WContainerWidget( parent ),
  m_gui( gui ),
  m_nuclide_edit( nullptr ),
  m_age_container( nullptr ),
  m_age_edit( nullptr ),
  m_fit_age( nullptr ),
  m_color_select( nullptr ),
  m_lower_container( nullptr ),
  m_add_constraint_btn( nullptr ),
  m_constraint( nullptr ),
  m_age_range_container( nullptr ),
  m_fit_age_min_edit( nullptr ),
  m_fit_age_max_edit( nullptr ),
  m_summary_text( nullptr ),
  m_updated( this ),
  m_remove( this ),
  m_fit_age_changed( this ),
  m_age_changed( this ),
  m_src_info( std::monostate() )
{
  addStyleClass( "RelActAutoGuiNuclide" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  
  WContainerWidget *upper_container = new WContainerWidget( this );
  upper_container->setStyleClass( "UpperRow" );


  WLabel *label = new WLabel( "Nuclide:", upper_container );
  m_nuclide_edit = new WLineEdit( "", upper_container );
  
  m_nuclide_edit->setAutoComplete( false );
  m_nuclide_edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclide_edit->setAttributeValue( "autocorrect", "off" );
  m_nuclide_edit->setAttributeValue( "spellcheck", "off" );
#endif
  label->setBuddy( m_nuclide_edit );
  m_nuclide_edit->setWidth( WLength(75.0, WLength::Pixel) );
  
  m_nuclide_edit->changed().connect( this, &RelActAutoGuiNuclide::handleIsotopeChange );
  
  // If we are typing in this box, we want to let app-hotkeys propagate up, but not arrow keys and
  //  stuff
  const string jsAppKeyDownFcn = wApp->javaScriptClass() + ".appKeyDown";
  const string keyDownJs = "function(s1,e1){"
  "if(e1 && e1.ctrlKey && e1.key && " + jsAppKeyDownFcn + ")"
  + jsAppKeyDownFcn + "(e1);"
  "}";
  m_nuclide_edit->keyWentDown().connect( keyDownJs );
  
  const char *tooltip = "ex. <b>U235</b>, <b>235 Uranium</b>"
  ", <b>U-235m</b> (meta stable state)"
  ", <b>Cs137</b>, Pb, Fe(n,n), etc.";
  HelpSystem::attachToolTipOn( m_nuclide_edit, tooltip, showToolTips );
  
  string replacerJs, matcherJs;
  IsotopeNameFilterModel::replacerJs( replacerJs );
  IsotopeNameFilterModel::nuclideNameMatcherJs( matcherJs );
  IsotopeNameFilterModel *isoSuggestModel = new IsotopeNameFilterModel( this );
  isoSuggestModel->excludeXrays( false );
  isoSuggestModel->excludeEscapes( false );
  isoSuggestModel->excludeReactions( false );
  
  WSuggestionPopup *nuclideSuggest = new WSuggestionPopup( matcherJs, replacerJs, this );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  nuclideSuggest->setJavaScriptMember("wtNoReparent", "true");
#endif
  nuclideSuggest->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
  nuclideSuggest->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
  
  IsotopeNameFilterModel::setQuickTypeFixHackjs( nuclideSuggest );
  
  isoSuggestModel->filter( "" );
  nuclideSuggest->setFilterLength( -1 );
  nuclideSuggest->setModel( isoSuggestModel );
  nuclideSuggest->filterModel().connect( isoSuggestModel, &IsotopeNameFilterModel::filter );
  nuclideSuggest->forEdit( m_nuclide_edit, WSuggestionPopup::Editing );  // | WSuggestionPopup::DropDownIcon
  
  m_age_container = new WContainerWidget( upper_container );
  m_age_container->setStyleClass( "RelActAutoGuiNuclideAgeContainer" );
  
  label = new WLabel( "Age:", m_age_container );
  m_age_edit = new WLineEdit( "", m_age_container );
  m_age_edit->setWidth( WLength(70.0, WLength::Pixel) );
  label->setBuddy( m_age_edit );
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), upper_container );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_age_edit->setValidator(validator);
  m_age_edit->setAutoComplete( false );
  m_age_edit->setAttributeValue( "ondragstart", "return false" );
  m_age_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeChange );
  m_age_edit->enterPressed().connect( this, &RelActAutoGuiNuclide::handleAgeChange );
  
  m_fit_age = new WCheckBox( "Fit Age", m_age_container );
  m_fit_age->addStyleClass( "CbNoLineBreak" );
  m_fit_age->setWordWrap( false );
  m_fit_age->checked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  m_fit_age->unChecked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  
  m_summary_text = new WText( upper_container );
  m_summary_text->setStyleClass( "RelActAutoGuiNuclideSummaryText" );

  //WContainerWidget *spacer = new WContainerWidget( upper_container );
  //spacer->addStyleClass( "RelActAutoSpacer" );
  
  
  m_color_select = new ColorSelect( ColorSelect::PrefferNative, upper_container );
  m_color_select->setColor( ns_default_color );
  
  m_color_select->cssColorChanged().connect( this, &RelActAutoGuiNuclide::handleColorChange );
  
  
  WPushButton *removeEnergyRange = new WPushButton( upper_container );
  removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
  removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeEnergyRange->clicked().connect( this, &RelActAutoGuiNuclide::handleRemoveSelf );


  m_lower_container = new WContainerWidget( this );
  m_lower_container->setStyleClass( "LowerRow" );

  m_add_constraint_btn = new WPushButton( "Add Constraint", m_lower_container );
  m_add_constraint_btn->addStyleClass( "LinkBtn AddConstraintBtn" );
  m_add_constraint_btn->clicked().connect( this, &RelActAutoGuiNuclide::handleAddConstraint );
  HelpSystem::attachToolTipOn( m_add_constraint_btn, "Add a constraint to this nuclide", showToolTips );
  
  
  m_constraint = new RelActAutoGuiNuclideConstraint( m_gui, this, m_lower_container );
  m_constraint->remove().connect( this, &RelActAutoGuiNuclide::handleRemoveConstraint );
  m_constraint->changed().connect( this, &RelActAutoGuiNuclide::handleConstraintChanged );
  m_constraint->hide();


  m_age_range_container = new WContainerWidget( m_lower_container );
  m_age_range_container->setStyleClass( "NucAgeRangeContainer" );

  WRegExpValidator *min_max_validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), upper_container );
  min_max_validator->setFlags(Wt::MatchCaseInsensitive);

  label = new WLabel( "Min Age:", m_age_range_container );
  m_fit_age_min_edit = new WLineEdit( m_age_range_container );
  label->setBuddy( m_fit_age_min_edit );
  m_fit_age_min_edit->setWidth( WLength(40.0, WLength::Pixel) );
  m_fit_age_min_edit->setValidator(min_max_validator);
  m_fit_age_min_edit->setAutoComplete( false );
  m_fit_age_min_edit->setAttributeValue( "ondragstart", "return false" );
  m_fit_age_min_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );
  m_fit_age_min_edit->enterPressed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );


  label = new WLabel( "Max Age:", m_age_range_container );
  m_fit_age_max_edit = new WLineEdit( m_age_range_container );
  label->setBuddy( m_fit_age_max_edit );
  m_fit_age_max_edit->setWidth( WLength(40.0, WLength::Pixel) );
  m_fit_age_max_edit->setValidator(min_max_validator);
  m_fit_age_max_edit->setAutoComplete( false );
  m_fit_age_max_edit->setAttributeValue( "ondragstart", "return false" );
  m_fit_age_max_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );
  m_fit_age_max_edit->enterPressed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );


  WContainerWidget *spacer = new WContainerWidget( m_lower_container );
  spacer->addStyleClass( "RelActAutoSpacer" );
  

  // Hide the age stuff by default, until a nuclide is selected
  m_age_container->hide();
  m_age_range_container->hide();
}//RelActAutoGuiNuclide
  

int RelActAutoGuiNuclide::relEffCurveIndex() const
{
  return m_gui->relEffCurveIndex( this );
}//size_t relEffCurveIndex() const


bool RelActAutoGuiNuclide::hasActRatioConstraint() const
{
  try
  {
    actRatioConstraint();
    return true;
  }catch( std::exception & )
  {
  }

  return false;
}

RelActCalcAuto::RelEffCurveInput::ActRatioConstraint RelActAutoGuiNuclide::actRatioConstraint() const
{
  if( !m_constraint || m_constraint->isHidden()  )
    throw runtime_error( "RelActAutoGuiNuclide::actRatioConstraint() called when no constraint is present" );

  return m_constraint->actRatioConstraint();
}

bool RelActAutoGuiNuclide::hasMassFractionConstraint() const
{
  try
  {
    massFractionConstraint();
    return true;
  }catch( std::exception & )
  {
  }

  return false;
}//bool RelActAutoGuiNuclide::hasMassFractionConstraint() const


RelActCalcAuto::RelEffCurveInput::MassFractionConstraint RelActAutoGuiNuclide::massFractionConstraint() const
{
  if( !m_constraint || m_constraint->isHidden() )
    throw runtime_error( "RelActAutoGuiNuclide::massFractionConstraint() called when no constraint is present" );

  return m_constraint->massFractionConstraint();
}//RelActCalcAuto::RelEffCurveInput::MassFractionConstraint massFractionConstraint() const


Wt::WColor RelActAutoGuiNuclide::getColorForSource( const RelActCalcAuto::SrcVariant &source ) const
{
  const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide( source );
  const SandiaDecay::Element * const el = RelActCalcAuto::element( source );
  const ReactionGamma::Reaction * const rctn = RelActCalcAuto::reaction( source );
  const string src_name = RelActCalcAuto::is_null(source) ? string() : RelActCalcAuto::to_name(source);

  const ReferencePhotopeakDisplay * const refdisp = InterSpec::instance()->referenceLinesWidget();
  if( refdisp )
  {
    const Wt::WColor c = refdisp->suggestColorForSource( src_name );
    if( !c.isDefault() )
      return c;
  }//if( refdisp )


  // Check peaks for the same source, and if so, use that color
  InterSpec * const interspec = InterSpec::instance();
  const PeakModel * const pmodel = interspec->peakModel();
  assert( pmodel );

  const shared_ptr<const deque<shared_ptr<const PeakDef>>> peaks = pmodel->peaks();
  assert( peaks );
  if( peaks )
  {
    for( const PeakModel::PeakShrdPtr &peak : *peaks )
    {
      assert( peak );
      if( !peak || peak->lineColor().isDefault() )
        continue;

      const SandiaDecay::Nuclide * const peak_nuc = peak->parentNuclide();
      const SandiaDecay::Element * const peak_el = peak->xrayElement();
      const ReactionGamma::Reaction * const peak_rctn = peak->reaction();
      if( (peak_nuc == nuc) && (peak_el == el) && (peak_rctn == rctn) )
        return peak->lineColor();
    }//for( const PeakModel::PeakShrdPtr &peak : *peaks )
  }//if( peaks )


  shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();
  if( theme )
  {
    const map<string,WColor> &defcolors = theme->referenceLineColorForSources;
    const auto pos = defcolors.find( src_name );
    if( pos != end(defcolors) )
      return pos->second;
  }//if( theme )


  // If we still havent found a color for this source, and the color picker is the default color,
  //  We'll use the first avaiable color of a predifined set of colors (the same one reference lines widget uses).

  //Use the same pre-defined colors as the reference lines widget, and use the first one that is not already used
  vector<Wt::WColor> def_line_colors;

  if( theme )
    def_line_colors = theme->referenceLineColor;

  if( def_line_colors.empty() )
    def_line_colors = ReferencePhotopeakDisplay::sm_def_line_colors;

  try
  {
    const int rel_eff_curve_index = relEffCurveIndex();
    const std::vector<RelActCalcAuto::NucInputInfo> nucs = m_gui->getNucInputInfo( rel_eff_curve_index );

    for( const RelActCalcAuto::NucInputInfo &other_nuc : nucs )
    {
      if( source == other_nuc.source )
        continue;

      const Wt::WColor c = WColor( other_nuc.peak_color_css );
      const auto pos = std::find( begin(def_line_colors), end(def_line_colors), c );
      if( pos != end(def_line_colors) )
        def_line_colors.erase( pos );
    }

    if( !def_line_colors.empty() )
      return def_line_colors[0];
  }catch( std::exception &e )
  {
    cerr << "RelActAutoGuiNuclide::handleIsotopeChange(): failed to get Rel Eff index to get nuc color." << endl;
  }


  // We've failed at everything else, so generate a random color...
  //  (we could use some better scheme to get a better color than random, but whatever at this point)
  const char letters[] = "0123456789ABCDEF";

  string hex = "#";
  for( size_t i = 0; i < 6; ++i )
    hex += letters[std::rand() % 16];
  return WColor( hex );
}//Wt::WColor getColorForSource() const


void RelActAutoGuiNuclide::handleIsotopeChange()
{
  //This function is called when the user changes the nuclide via the gui; not meant to be called if loading a
  //  serialized state.
  const auto prev_src_info = m_src_info;

  const auto nuc_input = source();
  m_src_info = nuc_input;
  
  const auto hide_age_stuff = [this,&nuc_input](){
    m_age_container->hide();
    m_age_range_container->hide();
    m_fit_age->setUnChecked();
    
    const string nucstr = m_nuclide_edit->text().toUTF8();
    if( !nucstr.empty() && std::holds_alternative<std::monostate>(nuc_input) )
      passMessage( nucstr + " is not a valid nuclide, x-ray, or reaction.", WarningWidget::WarningMsgHigh );
  };//hide_age_stuff
  
  if( nuc_input != prev_src_info )
  {
    updateAllowedConstraints();
    handleRemoveConstraint();
  }
  
  const string src_name = source_name();
  const SandiaDecay::Nuclide * const nuc = nuclide();
  const SandiaDecay::Element * const el = element();
  const ReactionGamma::Reaction * const rctn = reaction();
  
  if( !nuc && !el && !rctn )
  {
    hide_age_stuff();
    handleRemoveConstraint();
    m_updated.emit();
    return;
  }

  // We should only have at most, one of the source pointers be non-null
  assert( static_cast<int>(!!nuc) + static_cast<int>(!!el) + static_cast<int>(!!rctn) <= 1 );

  if( nuc )
  {
    if( IsInf(nuc->halfLife) )
    {
      const string nucstr = m_nuclide_edit->text().toUTF8();
      passMessage( nucstr + " is a stable nuclide.", WarningWidget::WarningMsgHigh );
      
      m_nuclide_edit->setText( "" );
      m_nuclide_edit->validate();
      m_fit_age->setUnChecked();
      m_fit_age->hide();
      
      m_updated.emit();
      return;
    }//if( IsInf(nuc->halfLife) )

    bool fit_age = false;
    const bool age_is_fittable = !PeakDef::ageFitNotAllowed( nuc );

    if( nuc->decaysToStableChildren() )
    {
      m_age_edit->setText( "0y" );
    }else
    {
      string agestr, fit_lower_age, fit_upper_age;

      try
      {
        const int rel_eff_index = relEffCurveIndex();

        if( m_gui->suggestInitialNuclideAge( rel_eff_index, nuc, agestr, fit_age, fit_lower_age, fit_upper_age ) )
        {
          m_fit_age_min_edit->setText( WString::fromUTF8(fit_lower_age) );
          m_fit_age_max_edit->setText( WString::fromUTF8(fit_upper_age) );
        }else
        {
          PeakDef::defaultDecayTime( nuc, &agestr );
          m_fit_age_min_edit->setText( "" );
          m_fit_age_max_edit->setText( "" );
        }

        m_age_edit->setText( WString::fromUTF8(agestr) );
      }catch( std::exception &e )
      {
        cerr << "RelActAutoGuiNuclide::handleIsotopeChange(): unable to get Rel Eff index - skipping getting" << endl;

        PeakDef::defaultDecayTime( nuc, &agestr );
        m_age_edit->setText( WString::fromUTF8(agestr) );
        m_fit_age_min_edit->setText( "" );
        m_fit_age_max_edit->setText( "" );
      }//try / catch get ages users want
    }//if( nuc->decaysToStableChildren() ) / else

    m_fit_age->setChecked( fit_age && age_is_fittable );
    m_age_container->setHidden( !age_is_fittable );
    m_age_range_container->setHidden( !age_is_fittable || !fit_age );
  }//if( nuc )
  
  if( el || rctn || !nuc )
    hide_age_stuff();

  if( RelActCalcAuto::is_null(nuc_input) )
  {
    m_color_select->setColor( ns_default_color );
  }else if( nuc_input != prev_src_info )
  {
    const WColor color = getColorForSource( nuc_input );
    m_color_select->setColor( color );

    ReferencePhotopeakDisplay * const refdisp = InterSpec::instance()->referenceLinesWidget();
    if( refdisp )
      refdisp->updateColorCacheForSource( src_name, color );
  }//if( nuc_input != prev_src_info )

  updateAllowedConstraints();
  
  m_updated.emit();
}//void handleIsotopeChange()
  
  
void RelActAutoGuiNuclide::setFitAge( const bool do_fit )
{
  m_fit_age->setChecked( do_fit );
  m_age_range_container->setHidden( !do_fit );
}//setFitAge(...)


bool RelActAutoGuiNuclide::fitAge() const
{
  return m_fit_age->isChecked();
}


void RelActAutoGuiNuclide::setAge( const string &age )
{
  m_age_edit->setValueText( WString::fromUTF8(age) );
}


void RelActAutoGuiNuclide::setAge( const Wt::WString &age )
{
  m_age_edit->setValueText( age );
}


void RelActAutoGuiNuclide::setAgeRange( Wt::WString min_age, Wt::WString max_age )
{
  m_fit_age_min_edit->setValueText( min_age );
  m_fit_age_max_edit->setValueText( max_age );
  
  validateAndCorrectAgeRange();
}//void setAgeRange( Wt::WString min_age, Wt::WString max_age )


std::pair<Wt::WString,Wt::WString> RelActAutoGuiNuclide::ageRangeStr() const
{
  return std::make_pair( m_fit_age_min_edit->text(), m_fit_age_max_edit->text() );
}


std::pair<std::optional<double>,std::optional<double>> RelActAutoGuiNuclide::ageRange() const
{
  if( !m_fit_age->isChecked() )
    return std::make_pair( std::nullopt, std::nullopt );

  const SandiaDecay::Nuclide * const nuc = nuclide();
  if( !nuc )
    return std::make_pair( std::nullopt, std::nullopt );
      
  const auto str_to_age = [nuc]( const string &str ) -> double {
    if( nuc )
      return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( str, nuc->halfLife );
    return PhysicalUnitsLocalized::stringToTimeDuration( str );
  };

  std::optional<double> min_age, max_age;

  try
  {
    if( !m_fit_age_min_edit->text().empty() )
      min_age = str_to_age( m_fit_age_min_edit->text().toUTF8() );
  }catch( const std::exception & )
  {
  }

  try
  {
    if( !m_fit_age_max_edit->text().empty() )
      max_age = str_to_age( m_fit_age_max_edit->text().toUTF8() );
  }catch( const std::exception & )
  {
  }
  
  return std::make_pair( min_age, max_age );
}//std::pair<std::optional<double>,std::optional<double>> ageRange() const


void RelActAutoGuiNuclide::setNuclideEditFocus()
{
  m_nuclide_edit->setFocus();
}//void setNuclideEditFocus()


void RelActAutoGuiNuclide::handleAgeChange()
{
  if( m_fit_age->isChecked() )
    validateAndCorrectAgeRange();
  

  m_age_changed.emit( this );
  m_updated.emit();
}//void handleAgeChange()


void RelActAutoGuiNuclide::handleAgeRangeChange()
{
  validateAndCorrectAgeRange();

  m_age_changed.emit( this );
  m_updated.emit();
}//void handleAgeRangeChange()


void RelActAutoGuiNuclide::validateAndCorrectAgeRange()
{
  const SandiaDecay::Nuclide * const nuc = nuclide();

  if( !nuc )
  {
    // We'll just make sure things are consistent - we probably don't need to do this
    assert( m_age_container->isHidden() );
    assert( m_age_range_container->isHidden() );
    assert( !m_fit_age->isChecked() );

    m_age_edit->setText( "" );
    m_fit_age_min_edit->setText( "" );
    m_fit_age_max_edit->setText( "" );

    m_age_container->hide();
    m_age_range_container->hide();
    m_fit_age->setUnChecked();
    return;
  }//if( !nuc )
  
      
  const auto str_to_age = [nuc]( const string &str ) -> double {
    if( nuc )
      return PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( str, nuc->halfLife );
    return PhysicalUnitsLocalized::stringToTimeDuration( str );
  };


  WString age_str = m_age_edit->text();
  WString min_age = m_fit_age_min_edit->text();
  WString max_age = m_fit_age_max_edit->text();
  double age = 0.0, min_age_val = 0.0, max_age_val = std::numeric_limits<double>::infinity();

  try
  {
    if( !age_str.empty() )
      age = str_to_age( age_str.toUTF8() );
  }catch( const std::exception & )
  {
  }   
  
  try
  {
    if( !min_age.empty() )
      min_age_val = str_to_age( min_age.toUTF8() );
  }catch( const std::exception & )
  {
  }
      
  try
  {
    if( !max_age.empty() )
      max_age_val = str_to_age( max_age.toUTF8() );
  }catch( const std::exception & )
  {
  }

  if( min_age_val > max_age_val )
  {
    std::swap( min_age, max_age );
    std::swap( min_age_val, max_age_val );
    m_fit_age_min_edit->setValueText( min_age );
    m_fit_age_max_edit->setValueText( max_age );
  }


  if( m_fit_age->isChecked() )
  {
    if( min_age_val > age )
      m_age_edit->setText( min_age );

    if( max_age_val < age )
      m_age_edit->setText( max_age );
  }//if( m_fit_age->isChecked() )
}//void handleAgeRangeChange()


void RelActAutoGuiNuclide::handleFitAgeChange()
{
  m_age_range_container->setHidden( !m_fit_age->isChecked() );
  if( m_fit_age->isChecked() )
    validateAndCorrectAgeRange();

  m_fit_age_changed.emit( this, m_fit_age->isChecked() );
  m_updated.emit();
}//void handleFitAgeChange()


void RelActAutoGuiNuclide::handleColorChange()
{
  const Wt::WColor c = m_color_select->color();
  const string src_name = source_name();
  if( !c.isDefault() && !src_name.empty() )
  {
    ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
    if( refdisp )
      refdisp->updateColorCacheForSource( src_name, c );
  }//if( !c.isDefault() )

  m_updated.emit();
}//void handleColorChange()


// Will return std::monostate if invalid text entered.
std::variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> RelActAutoGuiNuclide::source() const
{
  const string nucstr = m_nuclide_edit->text().toUTF8();
  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    throw runtime_error( "Couldnt load decay DB" );
  
  const SandiaDecay::Nuclide *nuc = db->nuclide( nucstr );
  if( nuc )
    return nuc;
  
  const SandiaDecay::Element *el = db->element( nucstr );
  if( el )
    return el;
  
  const ReactionGamma * const reaction_db = ReactionGammaServer::database();
  assert( reaction_db );
  if( !reaction_db )
    throw runtime_error( "Couldnt load reaction DB" );
  
  try
  {
    const ReactionGamma::Reaction *reaction = nullptr;
    vector<ReactionGamma::ReactionPhotopeak> possible_rctns;
    reaction_db->gammas( nucstr, possible_rctns );
    
    // TODO: we are currently taking the first reaction; however, in principle there could be multiple - however, `ReactionGamma` doesnt have an interface to just return a reaction by name, I guess because
    for( size_t i = 0; !reaction && (i < possible_rctns.size()); ++i )
      reaction = possible_rctns[i].reaction;
    
    if( reaction )
      return reaction;
  }catch( std::exception & )
  {
    //ReactionGamma::gammas(...) throws if not a valid reaction
  }//try / catch
  
  return std::monostate{};
}//variant<...> source() const


const SandiaDecay::Nuclide *RelActAutoGuiNuclide::nuclide() const
{
  const auto src = source();
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(src) )
    return std::get<const SandiaDecay::Nuclide *>(src);
  return nullptr;
}


const SandiaDecay::Element *RelActAutoGuiNuclide::element() const
{
  const auto src = source();
  if( std::holds_alternative<const SandiaDecay::Element *>(src) )
    return std::get<const SandiaDecay::Element *>(src);
  return nullptr;
}


const ReactionGamma::Reaction *RelActAutoGuiNuclide::reaction() const
{
  const auto src = source();
  if( std::holds_alternative<const ReactionGamma::Reaction *>(src) )
    return std::get<const ReactionGamma::Reaction *>(src);
  return nullptr;
}


string RelActAutoGuiNuclide::source_name() const
{
  const RelActCalcAuto::SrcVariant src = source();
  return RelActCalcAuto::is_null(src) ? string() : RelActCalcAuto::to_name(src);
}//string source_name() const


WColor RelActAutoGuiNuclide::color() const
{
  return m_color_select->color();
}


void RelActAutoGuiNuclide::setColor( const WColor &color )
{
  m_color_select->setColor( color );

  const string src_name = source_name();
  if( color.isDefault() || src_name.empty() )
    return;

  ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
  if( refdisp )
    refdisp->updateColorCacheForSource( src_name, color );
}//void setColor( const WColor &color )


double RelActAutoGuiNuclide::age() const
{ 
  const SandiaDecay::Nuclide * const nuc = nuclide();
  if( !nuc )
    return 0.0;
  
  if( m_age_container->isHidden() || m_age_edit->text().empty() )
    return PeakDef::defaultDecayTime( nuc );
  
  double age = 0.0;
  try
  {
    const string agestr = m_age_edit->text().toUTF8();
    age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, nuc->halfLife );
  }catch( std::exception & )
  {
    age = PeakDef::defaultDecayTime( nuc );
  }//try / catch
  
  return age;
}//double age() const


WString RelActAutoGuiNuclide::ageStr() const
{
  return m_age_edit->text();
}


RelActCalcAuto::NucInputInfo RelActAutoGuiNuclide::toNucInputInfo() const
{
  const auto nuc_input = source();
  
  if( RelActCalcAuto::is_null(nuc_input) )
    throw runtime_error( "No valid nuclide" );

  RelActCalcAuto::NucInputInfo nuc_info;
  nuc_info.source = nuc_input;

  if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
  {
    nuc_info.age = age(); // Must not be negative.
    nuc_info.fit_age = m_fit_age->isChecked();

    if( nuc_info.fit_age )
    {
      const auto age_range = ageRange();
      nuc_info.fit_age_min = age_range.first;
      nuc_info.fit_age_max = age_range.second;
    }
  }//if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
  

  switch(  m_constraint->constraintType() )
  {
    case NucConstraintType::None:
    case NucConstraintType::NumTypes:
    case NucConstraintType::MassFraction:
    case NucConstraintType::ActRatio:
      //These are given as constraints, not as input
      break;
    
    case NucConstraintType::RelActRange:
      nuc_info.min_rel_act = m_constraint->minRelAct();
      nuc_info.max_rel_act = m_constraint->maxRelAct();
      break;
  }//switch(  m_constraint->constraintType() )

  // TODO: Not implemented: vector<double> gammas_to_exclude;
  //nuc_info.gammas_to_exclude = ;
  
  nuc_info.peak_color_css = m_color_select->color().cssText();
  
  return nuc_info;
}//RelActCalcAuto::RoiRange toRoiRange() const


void RelActAutoGuiNuclide::fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info )
{
  handleRemoveConstraint();

  if( RelActCalcAuto::is_null(info.source) )
  {
    m_nuclide_edit->setText( "" );
    if( !info.peak_color_css.empty() )
      m_color_select->setColor( WColor(info.peak_color_css) );
    
    m_age_container->hide();
    m_age_edit->setText( "0s" );
    m_fit_age->setUnChecked();
    m_age_range_container->hide();

    updateAllowedConstraints();

    m_updated.emit();
    
    return;
  }//if( !info.nuclide && !info.element && !info.reaction )

  const std::string src_name = RelActCalcAuto::to_name(info.source);
  m_nuclide_edit->setText( WString::fromUTF8(src_name) );

  bool set_color = false;
  if( !info.peak_color_css.empty() )
  {
    const WColor c(info.peak_color_css);
    m_color_select->setColor( c );

    set_color = !c.isDefault();
    if( !c.isDefault() && !src_name.empty() )
    {
      ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
      if( refdisp )
        refdisp->updateColorCacheForSource( src_name, c );
    }//if( !c.isDefault() )
  }//if( !info.peak_color_css.empty() )

  const SandiaDecay::Nuclide * const nuc = RelActCalcAuto::nuclide(info.source);
  if( nuc )
  {
    const bool age_is_fittable = !PeakDef::ageFitNotAllowed(nuc);
    m_age_container->setHidden( !age_is_fittable );
    m_fit_age->setChecked( age_is_fittable && info.fit_age );
    m_age_range_container->setHidden( !age_is_fittable || !m_fit_age->isChecked() );
    
    if( !age_is_fittable || (info.age < 0.0) )
    {
      string agestr = "0s";
      if( nuc )
        PeakDef::defaultDecayTime( nuc, &agestr );
      m_age_edit->setText( WString::fromUTF8(agestr) );
    }else
    {
      const string agestr = PhysicalUnitsLocalized::printToBestTimeUnits(info.age);
      m_age_edit->setText( WString::fromUTF8(agestr) );
    }
    
    // TODO: it would be nice to print the times compact.  e.g., "20 y" instead of "20.000 y"
    WString min_str, max_str;
    if( info.fit_age_min.has_value() )
      min_str = WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(info.fit_age_min.value()));
    if( info.fit_age_max.has_value() )
      max_str = WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(info.fit_age_max.value()));

    m_fit_age_min_edit->setText( min_str );
    m_fit_age_max_edit->setText( max_str );
  }else
  {
    m_age_container->hide();
    m_fit_age->setUnChecked();
    m_age_range_container->hide();
    m_fit_age_min_edit->setText( "" );
    m_fit_age_max_edit->setText( "" );
  }//if( info.nuclide )
  
  m_constraint->setRelActRange( info.min_rel_act, info.max_rel_act );
  
  WString min_age_str, max_age_str;
  if( info.fit_age_min.has_value() )
    min_age_str = WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(info.fit_age_min.value(), 6) );
  if( info.fit_age_max.has_value() )
    max_age_str = WString::fromUTF8( PhysicalUnitsLocalized::printToBestTimeUnits(info.fit_age_max.value(), 6) );
  
  setAgeRange( min_age_str, max_age_str );
  
  //std::optional<double> info.starting_rel_act;
  //std::vector<double> info.gammas_to_exclude;
  // Not currently supported: info.gammas_to_exclude -> vector<double>;

  updateAllowedConstraints();

  const SandiaDecay::Element * const el = RelActCalcAuto::element(info.source);
  const ReactionGamma::Reaction * const rctn = RelActCalcAuto::reaction(info.source);

  if( !nuc )
  {
    m_age_container->hide();
    m_age_edit->setText( "0s" );
    m_fit_age->setUnChecked();
    m_age_range_container->hide();
  }

  if( !set_color && (el || rctn || nuc) )
  {
    assert( !src_name.empty() );
    const Wt::WColor c = getColorForSource( info.source );
    m_color_select->setColor( c );
    ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
    if( refdisp )
      refdisp->updateColorCacheForSource( src_name, c );
  }

  m_updated.emit();
}//void fromNucInputInfo( const RelActCalcAuto::NucInputInfo &info )


void RelActAutoGuiNuclide::handleRemoveSelf()
{
  m_remove.emit();
}


Wt::Signal<> &RelActAutoGuiNuclide::updated()
{
  return m_updated;
}


Wt::Signal<> &RelActAutoGuiNuclide::remove()
{
  return m_remove;
}


Wt::Signal<RelActAutoGuiNuclide *,bool> &RelActAutoGuiNuclide::fit_age_changed()
{
  return m_fit_age_changed;
}


Wt::Signal<RelActAutoGuiNuclide *> &RelActAutoGuiNuclide::age_changed()
{
  return m_age_changed;
}


void RelActAutoGuiNuclide::handleAddConstraint()
{
  if( m_constraint->isVisible() )
    return;
  
  m_add_constraint_btn->hide();
  m_constraint->show();

  m_updated.emit();
}//void handleAddConstraint()


void RelActAutoGuiNuclide::handleRemoveConstraint()
{
  if( !m_constraint->isVisible() )
    return;
  
  m_constraint->hide();
  m_add_constraint_btn->show();
  m_updated.emit();
}


void RelActAutoGuiNuclide::handleConstraintChanged()
{
  m_updated.emit();
}


void RelActAutoGuiNuclide::updateAllowedConstraints()
{
  m_constraint->updateAllowedConstraints();
}


void RelActAutoGuiNuclide::setSummaryText( const Wt::WString &text )
{
  m_summary_text->setText( text );
}


void RelActAutoGuiNuclide::addActRatioConstraint( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint )
{
  m_constraint->show();
  m_add_constraint_btn->show();

  assert( !RelActCalcAuto::is_null(constraint.constrained_source) );
  assert( !RelActCalcAuto::is_null(constraint.controlling_source) );

  const RelActCalcAuto::SrcVariant src = source();

  assert( !RelActCalcAuto::is_null(src) );
  if( RelActCalcAuto::is_null(src) )
    throw runtime_error( "Invalid ActRatioConstraint - not allowed when a source isnt defined." );

  assert( constraint.constrained_source == src );  
  if( constraint.constrained_source != src )
    throw runtime_error( "Invalid ActRatioConstraint - not not intended for this source." );
  
  m_constraint->setActRatio( constraint.controlling_source, constraint.constrained_to_controlled_activity_ratio );
}//void addActRatioConstraint( const ActRatioConstraint & )


void RelActAutoGuiNuclide::addMassFractionConstraint( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint )
{
  m_constraint->show();
  m_add_constraint_btn->show();

  const SandiaDecay::Nuclide *nuc = nuclide();  
  assert( constraint.nuclide == nuc );
  
  assert( nuc );
  if( !nuc )
    throw runtime_error( "Invalid MassFractionConstraint - not allowed for not a nuclide." );
  
  assert( constraint.nuclide == nuc );
  if( constraint.nuclide != nuc )
    throw runtime_error( "Invalid MassFractionConstraint - meant for another nuclide." );

  m_constraint->setMassFraction( constraint.lower_mass_fraction, constraint.upper_mass_fraction );
}//void addMassFractionConstraint( const MassFractionConstraint & )


