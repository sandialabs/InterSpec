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


#include <Wt/WLabel>
#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WLineEdit>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/ColorSelect.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/RelActCalcAuto.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/RelActAutoGuiNuclide.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"
#include "InterSpec/ReferencePhotopeakDisplay.h"


using namespace std;
using namespace Wt;


enum class NucConstraintType
{
  None, RelActRange, MassFraction, ActRatio
};//enum class NucConstraintType

class RelActAutoGuiNuclideConstraint : public WContainerWidget
{
  RelActAutoGuiNuclide *m_nuc;
  Wt::Signal<> m_remove;
  Wt::Signal<> m_changed;
  

public:
  RelActAutoGuiNuclideConstraint( RelActAutoGuiNuclide *nuc,WContainerWidget *parent )
    : WContainerWidget( parent ),
    m_remove( this ),
    m_nuc( nuc )
  {
    addStyleClass( "RelActAutoGuiNuclideConstraint" );

    const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );

    WContainerWidget *spacer = new WContainerWidget( this );
    spacer->setStyleClass( "RelActAutoSpacer" );

    WPushButton *removeBtn = new WPushButton( this );
    removeBtn->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
    removeBtn->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
    removeBtn->clicked().connect( this, &RelActAutoGuiNuclideConstraint::handleRemove );

    HelpSystem::attachToolTipOn( removeBtn, "Remove this constraint", showToolTips );



  }//RelActAutoGuiNuclideConstraint

  NucConstraintType constraintType() const
  {
    return NucConstraintType::None;
  }
  
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
    //m_min_rel_act_edit->setText( WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(min_rel_act)) );
    //m_max_rel_act_edit->setText( WString::fromUTF8(PhysicalUnitsLocalized::printToBestTimeUnits(max_rel_act)) );
  }

  std::optional<double> minRelAct() const
  {
    //return m_min_rel_act;
    return std::nullopt;
  }

  std::optional<double> maxRelAct() const
  {
    //return m_max_rel_act;
    return std::nullopt;
  }
  
};//class RelActAutoGuiNuclideConstraint



RelActAutoGuiNuclide::RelActAutoGuiNuclide( WContainerWidget *parent )
  : WContainerWidget( parent ),
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
  m_updated( this ),
  m_remove( this ),
  m_fit_age_changed( this ),
  m_age_changed( this )
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
  m_fit_age->setWordWrap( false );
  m_fit_age->checked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  m_fit_age->unChecked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  
  WContainerWidget *spacer = new WContainerWidget( upper_container );
  spacer->addStyleClass( "RelActAutoSpacer" );
  
  
  m_color_select = new ColorSelect( ColorSelect::PrefferNative, upper_container );
  m_color_select->setColor( WColor("#FF6633") );
  
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
  
  
  m_constraint = new RelActAutoGuiNuclideConstraint( this, m_lower_container );
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
  m_fit_age_min_edit->setWidth( WLength(50.0, WLength::Pixel) );
  m_fit_age_min_edit->setValidator(min_max_validator);
  m_fit_age_min_edit->setAutoComplete( false );
  m_fit_age_min_edit->setAttributeValue( "ondragstart", "return false" );
  m_fit_age_min_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );
  m_fit_age_min_edit->enterPressed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );


  label = new WLabel( "Max Age:", m_age_range_container );
  m_fit_age_max_edit = new WLineEdit( m_age_range_container );
  label->setBuddy( m_fit_age_max_edit );
  m_fit_age_max_edit->setWidth( WLength(50.0, WLength::Pixel) );
  m_fit_age_max_edit->setValidator(min_max_validator);
  m_fit_age_max_edit->setAutoComplete( false );
  m_fit_age_max_edit->setAttributeValue( "ondragstart", "return false" );
  m_fit_age_max_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );
  m_fit_age_max_edit->enterPressed().connect( this, &RelActAutoGuiNuclide::handleAgeRangeChange );


  spacer = new WContainerWidget( m_lower_container );
  spacer->addStyleClass( "RelActAutoSpacer" );
  

  // Hide the age stuff by default, until a nuclide is selected
  m_age_container->hide();
  m_age_range_container->hide();
}//RelActAutoGuiNuclide
  
  
void RelActAutoGuiNuclide::handleIsotopeChange()
{
  const auto nuc_input = nuclide();
  
  const auto hide_age_stuff = [this,&nuc_input](){
    m_age_container->hide();
    m_age_range_container->hide();
    m_fit_age->setUnChecked();
    
    const string nucstr = m_nuclide_edit->text().toUTF8();
    if( !nucstr.empty() && std::holds_alternative<std::monostate>(nuc_input) )
      passMessage( nucstr + " is not a valid nuclide, x-ray, or reaction.", WarningWidget::WarningMsgHigh );
  };//hide_age_stuff
  
  updateAllowedConstraints();

  if( std::holds_alternative<std::monostate>(nuc_input) )
  {
    hide_age_stuff();
    m_updated.emit();
    return;
  }
  
  std::string src_name;
  const SandiaDecay::Nuclide *nuc = nullptr;
  const SandiaDecay::Element *el = nullptr;
  const ReactionGamma::Reaction * reaction = nullptr;
  
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
  {
    nuc = std::get<const SandiaDecay::Nuclide *>(nuc_input);
    
    assert( nuc );
    if( !nuc )
      throw runtime_error( "No valid nuclide" );
    
    src_name = nuc->symbol;
    
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
    
    const bool age_is_fittable = !PeakDef::ageFitNotAllowed( nuc );
    
    if( nuc->decaysToStableChildren() )
    {
      m_age_edit->setText( "0y" );
    }else
    {
      string agestr;
      PeakDef::defaultDecayTime( nuc, &agestr );
      m_age_edit->setText( agestr );
    }
    
    m_age_container->setHidden( !age_is_fittable );
    m_fit_age->setChecked( false );
    m_age_range_container->setHidden( !age_is_fittable || !m_fit_age->isChecked() );
  }
  
  if( std::holds_alternative<const SandiaDecay::Element *>(nuc_input) )
  {
    el = std::get<const SandiaDecay::Element *>(nuc_input);
    assert( el );
    src_name = el->symbol;
    hide_age_stuff();
  }
  
  if( std::holds_alternative<const ReactionGamma::Reaction *>(nuc_input) )
  {
    reaction = std::get<const ReactionGamma::Reaction *>(nuc_input);
    assert( reaction );
    src_name = reaction->name();
    hide_age_stuff();
  }
  
  bool haveFoundColor = false;
  
  // Check user is showing reference lines, displayed peaks, and previous user-selected colors
  if( !haveFoundColor )
  {
    const ReferencePhotopeakDisplay *refdisp = InterSpec::instance()->referenceLinesWidget();
    if( refdisp )
    {
      const Wt::WColor c = refdisp->suggestColorForSource( src_name );
      if( !c.isDefault() )
      {
        haveFoundColor = true;
        m_color_select->setColor( c );
      }
    }//if( refdisp )
  }//if( !haveFoundColor )
  
  // Check if the theme explicitly specifies a color for this nuclide
  if( !haveFoundColor )
  {
    shared_ptr<const ColorTheme> theme = InterSpec::instance()->getColorTheme();
    if( theme )
    {
      const map<string,WColor> &defcolors = theme->referenceLineColorForSources;
      const auto pos = defcolors.find( src_name );
      if( pos != end(defcolors) )
        m_color_select->setColor( pos->second );
    }//if( theme )
  }//if( !haveFoundColor )
  
  if( !haveFoundColor )
  {
    // TODO: create a cache of previous user-selected colors
  }//if( !haveFoundColor )
  
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

  const variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> 
    src_info = nuclide();

  const SandiaDecay::Nuclide * const nuc = std::holds_alternative<const SandiaDecay::Nuclide *>(src_info) 
                                             ? std::get<const SandiaDecay::Nuclide *>(src_info) 
                                             : nullptr;
      
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
}


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
}


void RelActAutoGuiNuclide::validateAndCorrectAgeRange()
{
  const variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> 
    src_info = nuclide();

  const SandiaDecay::Nuclide * const nuc = std::holds_alternative<const SandiaDecay::Nuclide *>(src_info) 
                                             ? std::get<const SandiaDecay::Nuclide *>(src_info) 
                                             : nullptr;
      
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
}


void RelActAutoGuiNuclide::handleColorChange()
{
  m_updated.emit();
}


// Will return std::monostate if invalid text entered.
std::variant<std::monostate, const SandiaDecay::Nuclide *, const SandiaDecay::Element *, const ReactionGamma::Reaction *> RelActAutoGuiNuclide::nuclide() const
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
}


WColor RelActAutoGuiNuclide::color() const
{
  return m_color_select->color();
}


void RelActAutoGuiNuclide::setColor( const WColor &color )
{
  m_color_select->setColor( color );
}


double RelActAutoGuiNuclide::age() const
{
  const auto nuc_input = nuclide();
  if( !std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
    return 0.0;
  
  const SandiaDecay::Nuclide * const nuc = std::get<const SandiaDecay::Nuclide *>(nuc_input);
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
  const auto nuc_input = nuclide();
  
  
  RelActCalcAuto::NucInputInfo nuc_info;
  
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
  {
    nuc_info.nuclide = std::get<const SandiaDecay::Nuclide *>(nuc_input);
    nuc_info.age = age(); // Must not be negative.
    nuc_info.fit_age = m_fit_age->isChecked();

    if( nuc_info.fit_age )
    {
      const auto age_range = ageRange();
      nuc_info.fit_age_min = age_range.first;
      nuc_info.fit_age_max = age_range.second;
    }
  }else if( std::holds_alternative<const SandiaDecay::Element *>(nuc_input) )
  {
    nuc_info.element = std::get<const SandiaDecay::Element *>(nuc_input);
  }else if( std::holds_alternative<const ReactionGamma::Reaction *>(nuc_input) )
  {
    nuc_info.reaction = std::get<const ReactionGamma::Reaction *>(nuc_input);
  }else
  {
    throw runtime_error( "No valid nuclide" );
  }
  

  switch(  m_constraint->constraintType() )
  {
    case NucConstraintType::None:
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
  if( !info.nuclide && !info.element && !info.reaction )
  {
    m_nuclide_edit->setText( "" );
    if( !info.peak_color_css.empty() )
      m_color_select->setColor( WColor(info.peak_color_css) );
    
    m_age_container->hide();
    m_age_edit->setText( "0s" );
    m_fit_age->setUnChecked();
    m_age_range_container->hide();
    m_updated.emit();
    
    return;
  }//if( !info.nuclide && !info.element && !info.reaction )
  
  if( info.nuclide )
    m_nuclide_edit->setText( WString::fromUTF8(info.nuclide->symbol) );
  else if( info.element )
    m_nuclide_edit->setText( WString::fromUTF8(info.element->symbol) );
  else if( info.reaction )
    m_nuclide_edit->setText( WString::fromUTF8(info.reaction->name()) );
  
  if( !info.peak_color_css.empty() )
    m_color_select->setColor( WColor(info.peak_color_css) );
  
  if( info.nuclide )
  {
    const SandiaDecay::Nuclide * const nuc = info.nuclide;
    
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
    
    // TODO: blah blah blah - implement the below
    //std::optional<double> info.fit_age_min;
    //std::optional<double> info.fit_age_max;

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
  
  handleIsotopeChange();
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
  // TODO: blah blah blah - implement this
}


void RelActAutoGuiNuclide::addActRatioConstraint( const RelActCalcAuto::RelEffCurveInput::ActRatioConstraint &constraint )
{
  // TODO: blah blah blah - implement this
}


void RelActAutoGuiNuclide::addMassFractionConstraint( const RelActCalcAuto::RelEffCurveInput::MassFractionConstraint &constraint )
{
  // TODO: blah blah blah - implement this
}


void RelActAutoGuiNuclide::setIsInCurves( const std::set<size_t> &curves_with_nuc, size_t num_rel_eff_curves )
{
  // TODO: blah blah blah - implement this
}

