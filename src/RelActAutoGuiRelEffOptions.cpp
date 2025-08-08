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
#include <string>
#include <iostream>

#include <Wt/WLabel>
#include <Wt/WCheckBox>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WInPlaceEdit>
#include <Wt/WContainerWidget>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActAutoGui.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/SpectraFileModel.h"
#include "InterSpec/RelEffShieldWidget.h"
#include "InterSpec/RelActAutoGuiRelEffOptions.h"


using namespace Wt;
using namespace std;

RelActAutoGuiRelEffOptions::RelActAutoGuiRelEffOptions(RelActAutoGui *gui, Wt::WString name, Wt::WContainerWidget *parent)
    : Wt::WContainerWidget(parent),
      m_gui(gui),
      m_rel_eff_eqn_form(nullptr),
      m_eqn_order_div(nullptr),
      m_rel_eff_eqn_order_label(nullptr),
      m_rel_eff_eqn_order(nullptr),
      m_rel_eff_curve_name(nullptr),
      m_pu_corr_enabled( false ),
      m_pu_corr_div(nullptr),
      m_pu_corr_method(nullptr),
      m_phys_model_opts(nullptr),
      m_phys_model_shields(nullptr),
      m_phys_model_self_atten(nullptr),
      m_phys_ext_attens(nullptr),
      m_phys_model_use_hoerl(nullptr),
      m_phys_model_same_hoerl_on_all_curves(nullptr),
      m_phys_model_same_ext_shield_all_curves(nullptr),
      m_phys_model_shielded_by_other_curves(nullptr),
      m_eqn_txt(nullptr),
      m_add_del_rel_eff_div(nullptr),
      m_add_rel_eff_btn(nullptr),
      m_del_rel_eff_btn(nullptr),
      m_add_rel_eff_curve_signal(this),
      m_del_rel_eff_curve_signal(this),
      m_name_changed_signal(this),
      m_eqn_form_changed(this),
      m_same_hoerl_on_all_curves(this),
      m_same_ext_shield_on_all_curves(this),
      m_shielded_by_other_curves(this),
      m_options_changed_signal(this),
      m_has_multiple_phys_models(false)
{
  addStyleClass( "RelEffCurveOptions" );

  InterSpec * const viewer = InterSpec::instance();
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", viewer );

  m_rel_eff_curve_name = new WInPlaceEdit( this );
  m_rel_eff_curve_name->setStyleClass( "RelEffCurveName" );
  m_rel_eff_curve_name->setPlaceholderText( "Curve Name" );
  m_rel_eff_curve_name->setEmptyText( "Curve Name" );
  m_rel_eff_curve_name->setButtonsEnabled( false );
  m_rel_eff_curve_name->setText( name );
  m_rel_eff_curve_name->valueChanged().connect( this, &RelActAutoGuiRelEffOptions::emitNameChanged );

  WContainerWidget *eqnTypeDiv = new WContainerWidget( this );
  eqnTypeDiv->addStyleClass( "EqnTypeDiv" );
  
  WLabel *label = new WLabel( "Eqn Type", eqnTypeDiv );
  
  m_rel_eff_eqn_form = new WComboBox( eqnTypeDiv );
  m_rel_eff_eqn_form->addStyleClass( "GridSecondCol GridFirstRow" );
  label->setBuddy( m_rel_eff_eqn_form );
  m_rel_eff_eqn_form->activated().connect( this, &RelActAutoGuiRelEffOptions::handleRelEffEqnTypeChanged );
  
  
  const char *tooltip = "The functional form to use for the relative efficiciency curve.<br />"
  "Options are:"
  "<table style=\"margin-left: 10px;\">"
  "<tr><th>Log(energy):</th>               <th>y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...</th></tr>"
  "<tr><th>Log(rel. eff.):</th>            <th>y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )</th></tr>"
  "<tr><th>Log(energy)Log(rel. eff.):</th> <th>y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )</th></tr>"
  "<tr><th>Empirical:</th>            <th>y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )</th></tr>"
  "<tr><th>Physical:</th>             <th>"
  "  y = <span style=\"display: inline-block; vertical-align: middle;\">"
  "        <span style=\"display: block; text-align: center;\">"
  "         (1 - exp(-AD<sub>0</sub>*μ<sub>0</sub>))"
  "        </span>"
  "        <span style=\"display: block; border-top: 1px solid black; text-align: center;\">"
  "          (AD<sub>0</sub>*μ<sub>0</sub>)"
  "        </span>"
  "      </span>"
  "       * [exp(-AD<sub>1</sub>*μ<sub>1</sub>) * exp(-AD<sub>2</sub>*μ<sub>2</sub>) + ...]"
  "       * [Det. Eff.]"
  "       * [E<sup>b</sup> * C<sup>1/E</sup>]"
  "  </th></tr>"
  "</table>";
  HelpSystem::attachToolTipOn( {eqnTypeDiv}, tooltip, showToolTips );
  
  // Will assume FramEmpirical is the highest
  static_assert( static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel)
                > static_cast<int>(RelActCalc::RelEffEqnForm::LnXLnY),
                "RelEffEqnForm was changed!"
                );
  
  
  for( int i = 0; i <= static_cast<int>(RelActCalc::RelEffEqnForm::FramPhysicalModel); ++i )
  {
    const auto eqn_form = RelActCalc::RelEffEqnForm( i );
    
    const char *txt = "";
    switch( eqn_form )
    {
      case RelActCalc::RelEffEqnForm::LnX:
        //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
        txt = "Log(x)";
        break;
        
      case RelActCalc::RelEffEqnForm::LnY:
        //y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
        txt = "Log(y)";
        break;
        
      case RelActCalc::RelEffEqnForm::LnXLnY:
        //y = exp( a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
        txt = "Log(x)Log(y)";
        break;
        
      case RelActCalc::RelEffEqnForm::FramEmpirical:
        //y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
        txt = "Empirical";
        break;
        
      case RelActCalc::RelEffEqnForm::FramPhysicalModel:
        txt = "Physical";
        break;
    }
    
    m_rel_eff_eqn_form->addItem( txt );
  }//for( loop over RelEffEqnForm )
  
  m_rel_eff_eqn_form->setCurrentIndex( static_cast<int>(RelActCalc::RelEffEqnForm::LnX) );
  
  m_eqn_order_div = new WContainerWidget( this );
  m_eqn_order_div->addStyleClass( "EqnOrderDiv" );
  
  label = new WLabel( "Eqn Order", m_eqn_order_div );
  m_rel_eff_eqn_order = new WComboBox( m_eqn_order_div );
  label->setBuddy( m_rel_eff_eqn_order );
  m_rel_eff_eqn_order->activated().connect( m_gui, &RelActAutoGui::handleRelEffEqnOrderChanged );
  
  m_rel_eff_eqn_order->addItem( "0" );
  m_rel_eff_eqn_order->addItem( "1" );
  m_rel_eff_eqn_order->addItem( "2" );
  m_rel_eff_eqn_order->addItem( "3" );
  m_rel_eff_eqn_order->addItem( "4" );
  m_rel_eff_eqn_order->addItem( "5" );
  m_rel_eff_eqn_order->addItem( "6" );
  m_rel_eff_eqn_order->setCurrentIndex( 3 );
  
  tooltip = "The order (how many energy-dependent terms) relative efficiency equation to use.";
  HelpSystem::attachToolTipOn( {m_eqn_order_div}, tooltip, showToolTips );

  m_pu_corr_div = new WContainerWidget( this );
  m_pu_corr_div->addStyleClass( "PuCorrDiv" );
  
  label = new WLabel( "Pu242 corr", m_pu_corr_div );
  m_pu_corr_method = new WComboBox( m_pu_corr_div );
  label->setBuddy( m_pu_corr_method );
  m_pu_corr_method->activated().connect( m_gui, &RelActAutoGui::handlePuByCorrelationChanged );
  tooltip = "Pu-242 is often not directly observable in gamma spectra.  However, to"
  " correct for this isotope when calculating enrichment, the expected contributions of this"
  " isotope can be inferred from the other Pu isotopes."
  "  This form allows you to select the correction method.";
  HelpSystem::attachToolTipOn( {m_pu_corr_div}, tooltip, showToolTips );
  m_pu_corr_div->hide();
  m_pu_corr_enabled = false;

  
    
  m_phys_model_opts = new WContainerWidget( this );
  m_phys_model_opts->addStyleClass( "PhysicalModelOpts" );

  // We will pu the Use 
  WContainerWidget *phys_opt_row = new WContainerWidget( m_phys_model_opts );
  phys_opt_row->setStyleClass( "PhysicalModelOptRow" );

  SpecMeasManager *spec_manager = viewer->fileManager();
  SpectraFileModel *spec_model = spec_manager ? spec_manager->model() : nullptr;
  DetectorDisplay *det_disp = new DetectorDisplay( viewer, spec_model, phys_opt_row );

  WContainerWidget *phys_opt_cb_row = new WContainerWidget( phys_opt_row );
  phys_opt_cb_row->addStyleClass( "PhysOptionCbDiv" );
  m_phys_model_use_hoerl = new WCheckBox( "Use corr. fcn.", phys_opt_cb_row );
  m_phys_model_use_hoerl->addStyleClass( "UseCorrFcnCb CbNoLineBreak PhysOptionCb GridFirstRow GridFirstCol" );
  m_phys_model_use_hoerl->setChecked( true );
  m_phys_model_use_hoerl->checked().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
  m_phys_model_use_hoerl->unChecked().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );

  // Add the new checkboxes for shared settings
  m_phys_model_same_hoerl_on_all_curves = new WCheckBox( "Share corr. fcn.", phys_opt_cb_row );
  m_phys_model_same_hoerl_on_all_curves->addStyleClass( "SameHoerlAllCurvesCb CbNoLineBreak PhysOptionCb GridFirstRow GridSecondCol" );
  m_phys_model_same_hoerl_on_all_curves->setChecked( false );
  m_phys_model_same_hoerl_on_all_curves->checked().connect( this, &RelActAutoGuiRelEffOptions::handleSameHoerlOnAllCurvesChanged );
  m_phys_model_same_hoerl_on_all_curves->unChecked().connect( this, &RelActAutoGuiRelEffOptions::handleSameHoerlOnAllCurvesChanged );
  m_phys_model_same_hoerl_on_all_curves->hide();
  
  tooltip = "When checked, the Correlation Function setting will be synchronized across all Physical model curves.";
  HelpSystem::attachToolTipOn( {m_phys_model_same_hoerl_on_all_curves}, tooltip, showToolTips );


  m_phys_model_same_ext_shield_all_curves = new WCheckBox( "Share ext. atten.", phys_opt_cb_row );
  m_phys_model_same_ext_shield_all_curves->addStyleClass( "SameExtShieldAllCurvesCb PhysOptionCb GridSecondRow GridFirstCol" );
  m_phys_model_same_ext_shield_all_curves->setChecked( false );
  m_phys_model_same_ext_shield_all_curves->checked().connect( this, &RelActAutoGuiRelEffOptions::handleSameExternalShieldingChanged );
  m_phys_model_same_ext_shield_all_curves->unChecked().connect( this, &RelActAutoGuiRelEffOptions::handleSameExternalShieldingChanged );
  m_phys_model_same_ext_shield_all_curves->hide();
  
  tooltip = "When checked, the external shielding configuration will be synchronized across all Physical model curves.";
  HelpSystem::attachToolTipOn( {m_phys_model_same_ext_shield_all_curves}, tooltip, showToolTips );

  m_phys_model_shielded_by_other_curves = new WCheckBox( "Shielded by other curves", phys_opt_cb_row );
  m_phys_model_shielded_by_other_curves->addStyleClass( "ShieldedByOtherCurvesCb PhysOptionCb GridSecondRow GridSecondCol" );
  m_phys_model_shielded_by_other_curves->setChecked( false );
  m_phys_model_shielded_by_other_curves->checked().connect( this, &RelActAutoGuiRelEffOptions::handleShieldedByOtherCurvesChanged );
  m_phys_model_shielded_by_other_curves->unChecked().connect( this, &RelActAutoGuiRelEffOptions::handleShieldedByOtherCurvesChanged );
  m_phys_model_shielded_by_other_curves->hide();

  tooltip = "When checked, the sources of this curve will be shielded by all other Physical model curves.";
  HelpSystem::attachToolTipOn( {m_phys_model_shielded_by_other_curves}, tooltip, showToolTips );

  m_phys_model_shields = new WContainerWidget( m_phys_model_opts );
  m_phys_model_shields->addStyleClass( "PhysicalModelShields" );
    
  m_phys_ext_attens = new WContainerWidget( m_phys_model_shields );
  m_phys_model_opts->hide();

  //WContainerWidget *spacer = new WContainerWidget( this );
  //spacer->addStyleClass( "RelActAutoSpacer" );

  m_eqn_txt = new WText( this );
  m_eqn_txt->addStyleClass( "RelEffEqnTxt" );

  m_add_del_rel_eff_div = new WContainerWidget( this );
  m_add_del_rel_eff_div->addStyleClass( "AddDelRelEffDiv" );

  m_del_rel_eff_btn = new WPushButton( m_add_del_rel_eff_div );
  m_del_rel_eff_btn->setIcon( Wt::WLink( "InterSpec_resources/images/minus_min_black.svg" ) );
  m_del_rel_eff_btn->setStyleClass( "Wt-icon DelRelEffCurve" );
  m_del_rel_eff_btn->clicked().connect( this, &RelActAutoGuiRelEffOptions::emitDelRelEffCurve );
  m_del_rel_eff_btn->hide();
  
  tooltip = "Removes this relative efficiency curve.";
  HelpSystem::attachToolTipOn( {m_del_rel_eff_btn}, tooltip, showToolTips );

  m_add_rel_eff_btn = new WPushButton( m_add_del_rel_eff_div );
  m_add_rel_eff_btn->setIcon( Wt::WLink( "InterSpec_resources/images/plus_min_black.svg" ) );
  m_add_rel_eff_btn->setStyleClass( "Wt-icon AddRelEffCurve" );
  m_add_rel_eff_btn->clicked().connect( this, &RelActAutoGuiRelEffOptions::emitAddRelEffCurve );
  tooltip = "Add another relative efficiency curve.";
  HelpSystem::attachToolTipOn( {m_add_rel_eff_btn}, tooltip, showToolTips );
}

RelActAutoGuiRelEffOptions::~RelActAutoGuiRelEffOptions() {
    std::cout << "\n\n\n~RelActAutoGuiRelEffOptions\n\n\n";
}

void RelActAutoGuiRelEffOptions::showAndHideOptionsForEqnType()
{
  const RelActCalc::RelEffEqnForm eqn_type
                  = RelActCalc::RelEffEqnForm( std::max( 0, m_rel_eff_eqn_form->currentIndex() ) );
  
  const bool is_physical = (eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel);
  
  m_eqn_order_div->setHidden( is_physical );
  m_phys_model_opts->setHidden( !is_physical );
  
  // Hide multi-curve options if not a physical model
  if( !is_physical )
  {
    m_phys_model_same_hoerl_on_all_curves->hide();
    m_phys_model_same_ext_shield_all_curves->hide();
    m_phys_model_shielded_by_other_curves->hide();
  }else
  {
    // Show multi-curve options only if there are multiple physical models
    const bool show_shared_options = (m_has_multiple_phys_models && is_physical);
    m_phys_model_same_hoerl_on_all_curves->setHidden(!show_shared_options);
    m_phys_model_same_ext_shield_all_curves->setHidden(!show_shared_options);
    m_phys_model_shielded_by_other_curves->setHidden(!show_shared_options);
  }
  
  if( is_physical && !m_phys_model_self_atten )
    initPhysModelShields();
}

void RelActAutoGuiRelEffOptions::setHasMultiplePhysicalModels(const bool has_multiple_phys_models)
{
  if( m_has_multiple_phys_models != has_multiple_phys_models )
  {
    m_has_multiple_phys_models = has_multiple_phys_models;
    
    // Only show shared options if this is a physical model and there are multiple physical models
    const bool is_physical = (rel_eff_eqn_form() == RelActCalc::RelEffEqnForm::FramPhysicalModel);
    const bool show_shared_options = (m_has_multiple_phys_models && is_physical);

    m_phys_model_same_hoerl_on_all_curves->setHidden(!show_shared_options);
    m_phys_model_same_ext_shield_all_curves->setHidden(!show_shared_options);
    m_phys_model_shielded_by_other_curves->setHidden(!show_shared_options);
  }
}

bool RelActAutoGuiRelEffOptions::physModelSameHoerlOnAllCurves() const
{
  return m_phys_model_same_hoerl_on_all_curves->isChecked();
}

void RelActAutoGuiRelEffOptions::setPhysModelSameHoerlOnAllCurves(const bool same_hoerl_all_curves)
{
  m_phys_model_same_hoerl_on_all_curves->setChecked(same_hoerl_all_curves);
}

bool RelActAutoGuiRelEffOptions::physModelSameExtShieldAllCurves() const
{
  return m_phys_model_same_ext_shield_all_curves->isChecked();
}

void RelActAutoGuiRelEffOptions::setPhysModelSameExtShieldAllCurves(const bool same_ext_shield_all_curves)
{
  m_phys_model_same_ext_shield_all_curves->setChecked(same_ext_shield_all_curves);

  if( same_ext_shield_all_curves )
    m_phys_model_shielded_by_other_curves->setChecked( false );
}


bool RelActAutoGuiRelEffOptions::physModelShieldedByOtherCurves() const
{
  return m_phys_model_shielded_by_other_curves->isChecked();
}


void RelActAutoGuiRelEffOptions::setPhysModelShieldedByOtherCurves( const bool shielded_by_others )
{
  m_phys_model_shielded_by_other_curves->setChecked( shielded_by_others );
  if( shielded_by_others )
    m_phys_model_same_ext_shield_all_curves->setChecked( false );
}

void RelActAutoGuiRelEffOptions::handleRelEffEqnTypeChanged()
{
  showAndHideOptionsForEqnType();
  m_eqn_form_changed.emit( this );
}//void handleRelEffEqnTypeChanged()

void RelActAutoGuiRelEffOptions::handleSameHoerlOnAllCurvesChanged()
{
  m_same_hoerl_on_all_curves.emit(this);
}


void RelActAutoGuiRelEffOptions::handleSameExternalShieldingChanged()
{
  m_same_ext_shield_on_all_curves.emit(this);
}

void RelActAutoGuiRelEffOptions::handleShieldedByOtherCurvesChanged()
{
  m_shielded_by_other_curves.emit(this);
}

void RelActAutoGuiRelEffOptions::initPhysModelShields()
{
  if( m_phys_model_self_atten )
  {
    m_phys_model_self_atten->resetState();
  }else
  {
    m_phys_model_self_atten = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::SelfAtten );
    m_phys_model_shields->insertWidget( 0, m_phys_model_self_atten );
    m_phys_model_self_atten->changed().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
  }
  
  vector<RelEffShieldWidget *> starting_ext_shields;
  const vector<WWidget *> &ext_atten_widgets = m_phys_ext_attens->children();
  for( WWidget *w : ext_atten_widgets )
  {
    RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
    assert( sw );
    if( sw )
      starting_ext_shields.push_back( sw );
  }
  
  if( starting_ext_shields.empty() )
  {
    RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten,
                                                     m_phys_ext_attens );
    sw->changed().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
  }else
  {
    starting_ext_shields[0]->resetState();
    for( size_t i = 1; i < starting_ext_shields.size(); ++i )
      delete starting_ext_shields[i];
  }//if( starting_ext_shields.empty() )
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::addRelEffCurve()
{
  return m_add_rel_eff_curve_signal;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::delRelEffCurve()
{
  return m_del_rel_eff_curve_signal;
}

Wt::Signal<RelActAutoGuiRelEffOptions *,Wt::WString> &RelActAutoGuiRelEffOptions::nameChanged()
{
  return m_name_changed_signal;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::equationTypeChanged()
{
  return m_eqn_form_changed;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::sameHoerlOnAllCurvesChanged()
{
  return m_same_hoerl_on_all_curves;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::sameExternalShieldingChanged()
{
  return m_same_ext_shield_on_all_curves;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::shieldedByOtherCurvesChanged()
{
  return m_shielded_by_other_curves;
}

Wt::Signal<RelActAutoGuiRelEffOptions *> &RelActAutoGuiRelEffOptions::optionsChanged()
{
  return m_options_changed_signal;
}

void RelActAutoGuiRelEffOptions::emitAddRelEffCurve()
{
  m_add_rel_eff_curve_signal.emit( this );
}

void RelActAutoGuiRelEffOptions::emitDelRelEffCurve()
{
  m_del_rel_eff_curve_signal.emit( this );
}

void RelActAutoGuiRelEffOptions::emitNameChanged()
{
  m_name_changed_signal.emit( this, m_rel_eff_curve_name->text() );
}

void RelActAutoGuiRelEffOptions::emitOptionsChanged()
{
  m_options_changed_signal.emit( this );
}


void RelActAutoGuiRelEffOptions::setIsOnlyOneRelEffCurve( const bool is_only_rel_eff_curve )
{
  m_del_rel_eff_btn->setHidden( is_only_rel_eff_curve );
  m_rel_eff_curve_name->setHidden( is_only_rel_eff_curve );
  m_rel_eff_curve_name->setHidden( is_only_rel_eff_curve );
}

Wt::WString RelActAutoGuiRelEffOptions::name() const
{
  return m_rel_eff_curve_name->text();
}

void RelActAutoGuiRelEffOptions::setName( const Wt::WString &name )
{
  m_rel_eff_curve_name->setText( name );
  m_name_changed_signal.emit( this, name );
}

void RelActAutoGuiRelEffOptions::setRelEffCurveInput( const RelActCalcAuto::RelEffCurveInput &rel_eff )
{
  if( !rel_eff.name.empty() && (rel_eff.name != m_rel_eff_curve_name->text().toUTF8()) )
  {
    m_rel_eff_curve_name->setText( WString::fromUTF8(rel_eff.name) );
    m_name_changed_signal.emit( this, m_rel_eff_curve_name->text() );
  }

  m_rel_eff_eqn_form->setCurrentIndex( static_cast<int>(rel_eff.rel_eff_eqn_type) );
  showAndHideOptionsForEqnType();

  if( rel_eff.rel_eff_eqn_type == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    m_phys_model_use_hoerl->setChecked( rel_eff.phys_model_use_hoerl );
    initPhysModelShields();
    assert( m_phys_model_self_atten );
    
    if( rel_eff.phys_model_self_atten )
    {
      RelEffShieldState state;
      state.setStateFromFitInput( *rel_eff.phys_model_self_atten );
      m_phys_model_self_atten->setState( state );
    }else
    {
      m_phys_model_self_atten->resetState();
    }
    
    const size_t num_ext_atten = std::max( size_t(1), rel_eff.phys_model_external_atten.size() );

    // Remove any extra external attenuation widgets
    const vector<WWidget *> starting_ext_atten_widgets = m_phys_ext_attens->children();
    for( size_t i = num_ext_atten; i < starting_ext_atten_widgets.size(); ++i )
      delete starting_ext_atten_widgets[i];

    // Reset the state of any remaining external attenuation widgets (JIC - probably dont really need to do this)
    for( int i = 0; i < m_phys_ext_attens->children().size(); ++i )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[i] );
      assert( sw );
      if( sw )
        sw->resetState();
    }

    // Add any new external attenuation widgets
    while( m_phys_ext_attens->children().size() < num_ext_atten )
    {
      RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten,
                                                     m_phys_ext_attens );
      sw->changed().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
    }

    assert( m_phys_ext_attens->children().size() == num_ext_atten );

    for( size_t i = 0; i < std::min(num_ext_atten, m_phys_ext_attens->children().size()); ++i )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[i] );
      assert( sw );
      if( !sw )
        continue;
      
      if( (i < rel_eff.phys_model_external_atten.size()) && rel_eff.phys_model_external_atten[i] )
      {
        RelEffShieldState state;
        state.setStateFromFitInput( *rel_eff.phys_model_external_atten[i] );
        sw->setState( state );
      }else
      {
        sw->resetState();
      }
    }//for( size_t i = 0; i < std::min(num_ext_atten, m_phys_ext_attens->children().size()); ++i )

    m_phys_model_shielded_by_other_curves->setChecked( !rel_eff.shielded_by_other_phys_model_curve_shieldings.empty() );
  }else
  {
    const int eqn_order = std::max( 0, std::min(m_rel_eff_eqn_order->count() - 1, static_cast<int>(rel_eff.rel_eff_eqn_order) ) );
    m_rel_eff_eqn_order->setCurrentIndex( eqn_order );
  }

  
  if( rel_eff.pu242_correlation_method != RelActCalc::PuCorrMethod::NotApplicable )
  {
    // We will be a bit cheap here, and set one entry in pu correlation selector, and then
    //  rely on `updateDuringRenderForNuclideChange()` to input all the actual options
    m_pu_corr_method->clear();
    m_pu_corr_method->addItem( RelActCalc::to_description(rel_eff.pu242_correlation_method) );
    m_pu_corr_method->setCurrentIndex( m_pu_corr_method->count() - 1 );
  }

  updatePuCorrelationOptions( rel_eff.nuclides );
    

  // We are not dealing with below constraints, as we will display them in the nuclide areas
  //std::vector<ActRatioConstraint> rel_eff.act_ratio_constraints;
  //std::vector<MassFractionConstraint> rel_eff.mass_fraction_constraints;
}

void RelActAutoGuiRelEffOptions::updatePuCorrelationOptions( const vector<RelActCalcAuto::NucInputInfo> &nuclides )
{
  set<string> nuc_names;
  for( const RelActCalcAuto::NucInputInfo &src : nuclides )
  {
    const SandiaDecay::Nuclide *src_nuc = RelActCalcAuto::nuclide(src.source);
    if( src_nuc )
      nuc_names.insert( src_nuc->symbol );
  }//for( const RelActCalcAuto::NucInputInfo &src : nuclides )


  // We need Pu238, Pu239, and Pu240 for Bignan95_PWR and Bignan95_BWR.
  // We need Pu239 and at least one other isotope for ByPu239Only.
  const bool can_bignan = (nuc_names.count("Pu238") && nuc_names.count("Pu239")
                           && nuc_names.count("Pu240"));
  const bool can_239only = (nuc_names.count("Pu239")
                             && (nuc_names.count("Pu238") || nuc_names.count("Pu240")
                                  || nuc_names.count("Pu241") || nuc_names.count("Am241")) );
  
  
  const bool show_pu_corr = (can_bignan || can_239only);
  m_pu_corr_div->setHidden( !show_pu_corr );
  m_pu_corr_enabled = show_pu_corr;

  if( !show_pu_corr )
    return;

  const int nentries = m_pu_corr_method->count();
  const int nentries_needed = 1 + (can_239only ? 1 : 0) + (can_bignan ? 2 : 0);
  if( nentries != nentries_needed )
  {
    string current_txt;
    if( m_pu_corr_method->count() )
      current_txt = m_pu_corr_method->currentText().toUTF8();
      
    m_pu_corr_method->clear();
      
    const string &none = RelActCalc::to_description( RelActCalc::PuCorrMethod::NotApplicable );
    m_pu_corr_method->addItem( WString::fromUTF8(none) );
    const string &bignanBwr = RelActCalc::to_description( RelActCalc::PuCorrMethod::Bignan95_BWR );
    const string &bignanPwr = RelActCalc::to_description( RelActCalc::PuCorrMethod::Bignan95_PWR );
    const string &byPu239Only = RelActCalc::to_description( RelActCalc::PuCorrMethod::ByPu239Only );

    int next_index = 0;
    if( can_239only )
    {
      if( current_txt == byPu239Only )
        next_index = m_pu_corr_method->count();
      m_pu_corr_method->addItem( WString::fromUTF8(byPu239Only) );
    }//if( can_239only )
      
    if( can_bignan )
    {
      if( current_txt == bignanPwr )
        next_index = m_pu_corr_method->count();
      m_pu_corr_method->addItem( WString::fromUTF8(bignanPwr) );
        
      if( current_txt == bignanBwr )
        next_index = m_pu_corr_method->count();
      m_pu_corr_method->addItem( WString::fromUTF8(bignanBwr) );
    }//if( can_bignan )
      
    m_pu_corr_method->setCurrentIndex( next_index );
  }//if( nentries != nentries_needed )
}


void RelActAutoGuiRelEffOptions::setEqnTxt( const Wt::WString &txt )
{
  m_eqn_txt->setText( txt );
}//setEqnTxt


void RelActAutoGuiRelEffOptions::update_shield_widget( 
          const RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo &shield,
                                      RelEffShieldWidget *w ) 
{
  if( !w )
    return;
      
  const Material * const widget_mat = w->material();
  assert( (!!shield.material) == (!!widget_mat) );
  assert( !widget_mat || !shield.material || (widget_mat->name == shield.material->name) );
      
  if( (!!shield.material) != (!!widget_mat) )
    w->setMaterialSelected( !!shield.material );
      
  if( shield.material && (!widget_mat || (widget_mat->name != shield.material->name)) )
    w->setMaterial( shield.material->name );
      
  if( shield.atomic_number )
  {
    assert( shield.atomic_number_was_fit == w->fitAtomicNumber() );
    assert( shield.atomic_number_was_fit || ((*shield.atomic_number) == w->atomicNumber()));
    if( (*shield.atomic_number) != w->atomicNumber() )
      w->setAtomicNumber( *shield.atomic_number );
        
    if( shield.atomic_number_was_fit != w->fitAtomicNumber() )
      w->setFitAtomicNumber( shield.atomic_number_was_fit );
  }//if( shield.atomic_number )
      
  
  if( shield.material )
  {
    assert( shield.areal_density_was_fit == w->fitThickness() );
    const double thick = shield.areal_density / shield.material->density;
    
    if( !shield.areal_density_was_fit && (thick != w->thickness()) )
    {
      cerr << "Found fixed shielding thickness has changed for '" << shield.material->name << "'; was "
      << PhysicalUnits::printToBestLengthUnits(w->thickness()) << " but returned answer is "
      << PhysicalUnits::printToBestLengthUnits(thick) << endl;
      cerr << endl;
    }
    assert( shield.areal_density_was_fit || (thick == w->thickness()));
    
    w->setThickness( thick );
    if( shield.areal_density_was_fit != w->fitThickness() )
      w->setFitThickness( shield.areal_density_was_fit );
  }else
  {
    const double ad_g_cm2 = shield.areal_density / PhysicalUnits::g_per_cm2;
    
    assert( shield.areal_density_was_fit == w->fitArealDensity() );
    assert( shield.areal_density_was_fit || (ad_g_cm2 == w->arealDensity()));
    
    if( ad_g_cm2 != w->arealDensity() )
      w->setArealDensity( ad_g_cm2 );
    
    if( shield.areal_density_was_fit != w->fitArealDensity() )
      w->setFitArealDensity( shield.areal_density_was_fit );
  }//if( shield.material ) / else
}

void RelActAutoGuiRelEffOptions::update_self_atten_shield_widget(
                    const optional<RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo> &self_atten )
{
  if( !m_phys_model_self_atten )
    initPhysModelShields();
  
  if( self_atten )
    update_shield_widget( *self_atten, m_phys_model_self_atten );
  else
    m_phys_model_self_atten->resetState();
}


void RelActAutoGuiRelEffOptions::update_external_atten_shield_widget(
                     const vector<RelActCalcAuto::RelActAutoSolution::PhysicalModelFitInfo::ShieldInfo> &ext_shields )
{
  if( !m_phys_model_self_atten )
    initPhysModelShields();

  const size_t num_ext_atten = std::max( size_t(1), ext_shields.size() );
  while( m_phys_ext_attens->children().size() > num_ext_atten )
  {
    delete m_phys_ext_attens->children().back();
  }

  while( m_phys_ext_attens->children().size() < num_ext_atten )
  {
    RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten,
                                                     m_phys_ext_attens );
    sw->changed().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
  }

  if( ext_shields.empty() )
  {
    assert( m_phys_ext_attens->children().size() == 1 );
    RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[0] );
    assert( sw );
    if( sw )
      sw->resetState();
  }else
  {
    const vector<WWidget *> kids = m_phys_ext_attens->children();
    assert( kids.size() == num_ext_atten );
    for( size_t i = 0; i < kids.size(); ++i )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( kids[i] );
      assert( sw );
      if( sw && (i < ext_shields.size()) )
        update_shield_widget( ext_shields[i], sw );
      else if( sw )
        sw->resetState();
    }//for( size_t i = 0; i < num_ext_atten; ++i )
  }//if( ext_shields.empty() ) / else
}//update_external_atten_shield_widget(...)


void RelActAutoGuiRelEffOptions::update_self_atten_shield_widget( const RelEffShieldWidget *shield )
{
  if( !m_phys_model_self_atten )
    initPhysModelShields();
  
  assert( m_phys_model_self_atten );
  if( !m_phys_model_self_atten )
    return;
  
  if( shield )
  {
    std::unique_ptr<RelEffShieldState> state = shield->state();
    if( state )
      m_phys_model_self_atten->setState( *state );
  }else
  {
    m_phys_model_self_atten->resetState();
  }
}//void update_self_atten_shield_widget(...)


void RelActAutoGuiRelEffOptions::update_external_atten_shield_widget( const std::vector<const RelEffShieldWidget *> &ext_shields )
{
  if( !m_phys_model_self_atten )
    initPhysModelShields();
  
  const size_t num_ext_atten = std::max( size_t(1), ext_shields.size() );
  while( m_phys_ext_attens->children().size() > num_ext_atten )
  {
    delete m_phys_ext_attens->children().back();
  }
  
  while( m_phys_ext_attens->children().size() < num_ext_atten )
  {
    RelEffShieldWidget *sw = new RelEffShieldWidget( RelEffShieldWidget::ShieldType::ExternalAtten,
                                                    m_phys_ext_attens );
    sw->changed().connect( this, &RelActAutoGuiRelEffOptions::emitOptionsChanged );
  }
  
  if( ext_shields.empty() )
  {
    assert( m_phys_ext_attens->children().size() == 1 );
    RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[0] );
    assert( sw );
    if( sw )
      sw->resetState();
  }else
  {
    assert( m_phys_ext_attens->children().size() == num_ext_atten );
    const vector<WWidget *> childs = m_phys_ext_attens->children();
    
    for( size_t i = 0; i < childs.size(); ++i )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[i] );
      assert( sw );

      if( sw && (i < ext_shields.size()) && ext_shields[i] )
      {
        unique_ptr<RelEffShieldState> state = ext_shields[i]->state();
        if( state )
          sw->setState( *state );
      }else if( sw )
      {
        sw->resetState();
      }
    }//for( size_t i = 0; i < num_ext_atten; ++i )
  }//if( ext_shields.empty() ) / else
}//void update_external_atten_shield_widget(...)


const RelEffShieldWidget *RelActAutoGuiRelEffOptions::selfAttenWidget() const
{
  return m_phys_model_self_atten;
}


vector<const RelEffShieldWidget *> RelActAutoGuiRelEffOptions::externalAttenWidgets() const
{
  vector<const RelEffShieldWidget *> answer;
  const vector<WWidget *> childs = m_phys_ext_attens->children();
  
  for( size_t i = 0; i < childs.size(); ++i )
  {
    RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( m_phys_ext_attens->children()[i] );
    assert( sw );
    if( sw )
      answer.push_back( sw );
  }
  
  return answer;
}//externalAttenWidgets()


RelActCalc::RelEffEqnForm RelActAutoGuiRelEffOptions::rel_eff_eqn_form() const
{
  return RelActCalc::RelEffEqnForm( std::max( 0, m_rel_eff_eqn_form->currentIndex() ) );
}

size_t RelActAutoGuiRelEffOptions::rel_eff_eqn_order() const
{
  const RelActCalc::RelEffEqnForm form = rel_eff_eqn_form();
  if( form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return 0;

  return std::max( 0, m_rel_eff_eqn_order->currentIndex() );
}

std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> RelActAutoGuiRelEffOptions::phys_model_self_atten() const
{
  const RelActCalc::RelEffEqnForm form = rel_eff_eqn_form();
  if( form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return nullptr;

  if( !m_phys_model_self_atten || !m_phys_model_self_atten->nonEmpty() )
    return nullptr;

  return m_phys_model_self_atten->fitInput();
}

vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>> RelActAutoGuiRelEffOptions::phys_model_external_atten() const
{
  vector<shared_ptr<const RelActCalc::PhysicalModelShieldInput>> answer;
  const RelActCalc::RelEffEqnForm form = rel_eff_eqn_form();
  if( form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return answer;

  for( WWidget *w : m_phys_ext_attens->children() )
    {
      RelEffShieldWidget *sw = dynamic_cast<RelEffShieldWidget *>( w );
      if( sw && sw->nonEmpty() )
      {
        auto input = sw->fitInput();
        if( input )
          answer.push_back( input );
      }
    }//for( WWidget *w : m_phys_ext_attens->children() )

  return answer;
}

bool RelActAutoGuiRelEffOptions::phys_model_use_hoerl() const
{
  const RelActCalc::RelEffEqnForm form = rel_eff_eqn_form();
  if( form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return false;
  return m_phys_model_use_hoerl->isChecked();
}

void RelActAutoGuiRelEffOptions::setPhysModelUseHoerl( const bool use_hoerl )
{
  m_phys_model_use_hoerl->setChecked( use_hoerl );
}

RelActCalc::PuCorrMethod RelActAutoGuiRelEffOptions::pu242_correlation_method() const
{
  if( !m_pu_corr_enabled )
    return RelActCalc::PuCorrMethod::NotApplicable;

  
  const string currtxt = m_pu_corr_method->currentText().toUTF8();
  for( int i = 0; i <= static_cast<int>(RelActCalc::PuCorrMethod::NotApplicable); ++i )
  {
    const auto method = RelActCalc::PuCorrMethod(i);
    const string &desc = RelActCalc::to_description( method );
    if( desc == currtxt )
      return method;
  }//for( loop over RelActCalc::PuCorrMethod )
  
  assert( 0 );
  return RelActCalc::PuCorrMethod::NotApplicable;
} 
