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


RelActAutoGuiNuclide::RelActAutoGuiNuclide( WContainerWidget *parent )
  : WContainerWidget( parent ),
  m_nuclide_edit( nullptr ),
  m_age_label( nullptr ),
  m_age_edit( nullptr ),
  m_fit_age( nullptr ),
  m_color_select( nullptr ),
  m_updated( this ),
  m_remove( this )
{
  addStyleClass( "RelActAutoGuiNuclide" );
  
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", InterSpec::instance() );
  
  WLabel *label = new WLabel( "Nuclide:", this );
  m_nuclide_edit = new WLineEdit( "", this );
  
  m_nuclide_edit->setAutoComplete( false );
  m_nuclide_edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclide_edit->setAttributeValue( "autocorrect", "off" );
  m_nuclide_edit->setAttributeValue( "spellcheck", "off" );
#endif
  label->setBuddy( m_nuclide_edit );
  
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
  
  m_age_label = new WLabel( "Age:", this );
  m_age_edit = new WLineEdit( "", this );
  m_age_label->setBuddy( m_age_edit );
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_age_edit->setValidator(validator);
  m_age_edit->setAutoComplete( false );
  m_age_edit->setAttributeValue( "ondragstart", "return false" );
  m_age_edit->changed().connect( this, &RelActAutoGuiNuclide::handleAgeChange );
  
  m_fit_age = new WCheckBox( "Fit Age", this );
  m_fit_age->setWordWrap( false );
  m_fit_age->checked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  m_fit_age->unChecked().connect( this, &RelActAutoGuiNuclide::handleFitAgeChange );
  
  WContainerWidget *spacer = new WContainerWidget( this );
  spacer->addStyleClass( "RelActAutoGuiNuclideSpacer" );
  
  
  m_color_select = new ColorSelect( ColorSelect::PrefferNative, this );
  m_color_select->setColor( WColor("#FF6633") );
  
  m_color_select->cssColorChanged().connect( this, &RelActAutoGuiNuclide::handleColorChange );
  
  
  WPushButton *removeEnergyRange = new WPushButton( this );
  removeEnergyRange->setStyleClass( "DeleteEnergyRangeOrNuc Wt-icon" );
  removeEnergyRange->setIcon( "InterSpec_resources/images/minus_min_black.svg" );
  removeEnergyRange->clicked().connect( this, &RelActAutoGuiNuclide::handleRemoveSelf );
}//RelActAutoGuiNuclide
  
  
void RelActAutoGuiNuclide::handleIsotopeChange()
{
  const auto nuc_input = nuclide();
  
  const auto hide_age_stuff = [this,&nuc_input](){
    m_age_label->hide();
    m_age_edit->hide();
    m_fit_age->setUnChecked();
    m_fit_age->hide();
    
    const string nucstr = m_nuclide_edit->text().toUTF8();
    if( !nucstr.empty() && std::holds_alternative<std::monostate>(nuc_input) )
      passMessage( nucstr + " is not a valid nuclide, x-ray, or reaction.", WarningWidget::WarningMsgHigh );
  };//hide_age_stuff
  
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
    
    m_age_label->setHidden( !age_is_fittable );
    m_age_edit->setHidden( !age_is_fittable );
    m_fit_age->setHidden( !age_is_fittable );
    m_fit_age->setChecked( false );
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
  
  
void RelActAutoGuiNuclide::setFitAgeVisible( bool visible, bool do_fit )
{
  if( visible )
  {
    m_fit_age->setHidden( false );
    m_fit_age->setChecked( do_fit );
  }else
  {
    m_fit_age->setHidden( true );
    m_fit_age->setChecked( false );
  }
}//setFitAgeVisible(...)
  

void RelActAutoGuiNuclide::setAgeDisabled( bool disabled )
{
  m_fit_age->setDisabled( disabled );
  m_age_edit->setDisabled( disabled );
}//void setAgeDisabled( bool disabled )


void RelActAutoGuiNuclide::setAge( const string &age )
{
  m_age_edit->setValueText( WString::fromUTF8(age) );
}


void RelActAutoGuiNuclide::setNuclideEditFocus()
{
  m_nuclide_edit->setFocus();
}


void RelActAutoGuiNuclide::handleAgeChange()
{
  m_updated.emit();
}//void handleAgeChange()


void RelActAutoGuiNuclide::handleFitAgeChange()
{
  
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
  
  if( m_age_edit->isHidden() || m_age_edit->text().empty() )
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


RelActCalcAuto::NucInputInfo RelActAutoGuiNuclide::toNucInputInfo() const
{
  const auto nuc_input = nuclide();
  
  
  RelActCalcAuto::NucInputInfo nuc_info;
  
  if( std::holds_alternative<const SandiaDecay::Nuclide *>(nuc_input) )
  {
    nuc_info.nuclide = std::get<const SandiaDecay::Nuclide *>(nuc_input);
    nuc_info.age = age(); // Must not be negative.
    nuc_info.fit_age = m_fit_age->isChecked();
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
    
    m_age_label->hide();
    m_age_edit->hide();
    m_age_edit->setText( "0s" );
    m_fit_age->setUnChecked();
    m_fit_age->hide();
    
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
    m_age_label->setHidden( !age_is_fittable );
    m_age_edit->setHidden( !age_is_fittable );
    m_fit_age->setHidden( !age_is_fittable );
    m_fit_age->setChecked( age_is_fittable && info.fit_age );
    
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
  }else
  {
    m_age_label->hide();
    m_age_edit->hide();
    m_fit_age->setUnChecked();
    m_fit_age->hide();
  }//if( info.nuclide )
  
  // Not currently supported: info.gammas_to_exclude -> vector<double>;
  
  // TODO: blah blah blah - implement the below
  //std::optional<double> info.min_rel_act;
  //std::optional<double> info.max_rel_act;
  //std::optional<double> info.starting_rel_act;
  //std::vector<double> info.gammas_to_exclude;
  
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

