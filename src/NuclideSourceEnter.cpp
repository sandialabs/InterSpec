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

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WContainerWidget>

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/NuclideSourceEnter.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"


using namespace std;
using namespace Wt;



NuclideSourceEnterController::NuclideSourceEnterController( Wt::WLineEdit *nuclideEdit,
                               Wt::WLineEdit *nuclideAgeEdit,
                               Wt::WText *halfLifeTxt,
                               Wt::WObject *parent )
  : WObject( parent ),
  m_nuclideEdit( nuclideEdit ),
  m_nuclideAgeEdit( nuclideAgeEdit ),
  m_halfLifeTxt( halfLifeTxt ),
  m_currentNuc( nullptr ),
  m_prevAgeTxt{},
  m_changed( this )
{
  assert( nuclideEdit );
  assert( nuclideAgeEdit );
    
  m_nuclideEdit->setAutoComplete( false );
  m_nuclideEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclideEdit->setAttributeValue( "autocorrect", "off" );
  m_nuclideEdit->setAttributeValue( "spellcheck", "off" );
#endif
  m_nuclideEdit->setAttributeValue( "spellcheck", "false" );
    
  m_nuclideAgeEdit->setAutoComplete( false );
  m_nuclideAgeEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_nuclideAgeEdit->setAttributeValue( "autocorrect", "off" );
  m_nuclideAgeEdit->setAttributeValue( "spellcheck", "off" );
#endif
    
    
  string replacerJs, matcherJs;
  PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
  PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
    
  WSuggestionPopup *suggestions = new WSuggestionPopup( matcherJs, replacerJs, this );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
  suggestions->setJavaScriptMember("wtNoReparent", "true");
#endif
    
  suggestions->addStyleClass( "nuclide-suggest" );
  suggestions->forEdit( m_nuclideEdit, WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

    
  IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel( this );
    
  filterModel->excludeNuclides( false );
  filterModel->excludeXrays( true );
  filterModel->excludeEscapes( true );
  filterModel->excludeReactions( true );
  
  filterModel->filter( "" );
  suggestions->setFilterLength( -1 );
  suggestions->setModel( filterModel );
  suggestions->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
    
  m_nuclideEdit->changed().connect( this, &NuclideSourceEnterController::handleNuclideUserInput );
  m_nuclideEdit->enterPressed().connect( this, &NuclideSourceEnterController::handleNuclideUserInput );
//    m_nuclideEdit->blurred().connect( this, &NuclideSourceEnterController::handleNuclideUserInput );
    
  IsotopeNameFilterModel::setQuickTypeFixHackjs( suggestions );
  IsotopeNameFilterModel::setEnterKeyMatchFixJs( suggestions, m_nuclideEdit );
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
  validator->setFlags(Wt::MatchCaseInsensitive);
  m_nuclideAgeEdit->setValidator(validator);
  
  m_nuclideAgeEdit->changed().connect( this, &NuclideSourceEnterController::handleAgeUserChange );
  m_nuclideAgeEdit->blurred().connect( this, &NuclideSourceEnterController::handleAgeUserChange );
  m_nuclideAgeEdit->enterPressed().connect( this, &NuclideSourceEnterController::handleAgeUserChange );
    
}//NuclideSourceEnterController constructor
  
  
  
void NuclideSourceEnterController::handleNuclideUserInput()
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const string isotopeLabel = m_nuclideEdit->text().toUTF8();
  const SandiaDecay::Nuclide *nuc = db->nuclide( isotopeLabel );
  
  if( nuc == m_currentNuc )
    return;

  m_currentNuc = nuc;
  
  if( m_halfLifeTxt )
  {
    WString hlstr;
    if( nuc )
      hlstr = PhysicalUnitsLocalized::printToBestTimeUnits( nuc->halfLife, 2, SandiaDecay::second );
    m_halfLifeTxt->setText( WString("{1}={2}").arg(WString::tr("T1/2")).arg(hlstr) );
  }//if( m_halfLifeTxt )

  bool useCurrentAge = false;
  const bool showPromptOnly = false; //m_promptLinesOnly
  
  string agestr = useCurrentAge ? m_nuclideAgeEdit->text().toUTF8() : string();
  if( m_currentNuc && m_prevAgeTxt.count(m_currentNuc) )
  {
    useCurrentAge = true;
    agestr = m_prevAgeTxt[m_currentNuc];
    try
    {
      const double hl = (m_currentNuc ? m_currentNuc->halfLife : -1.0);
      double x = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
      if( x < 0.0 )
        throw runtime_error("");
    }catch(...)
    {
      useCurrentAge = false;
      agestr = "";
      m_prevAgeTxt.erase( m_currentNuc );
    }
  }else
  {
    agestr = "";
    useCurrentAge = false;
  }
  
  m_nuclideAgeEdit->enable();
  
  if( nuc )
  {
    const bool prompt = (nuc->canObtainPromptEquilibrium() && showPromptOnly);
    const bool notMuchEvolution = PeakDef::ageFitNotAllowed(nuc);
    
    if( (prompt || notMuchEvolution) && m_nuclideAgeEdit->isEnabled() )
    {
      m_nuclideAgeEdit->disable();
      m_nuclideAgeEdit->setText( "" );
    }
  }//if( nuc )
  
  try
  {
    if( nuc && !IsInf(nuc->halfLife) && !nuc->decaysToChildren.empty() && !useCurrentAge )
    {
      if( agestr.size() )
        m_nuclideAgeEdit->setText( agestr );
    
      if( PeakDef::ageFitNotAllowed(nuc) )
      {
        m_nuclideAgeEdit->setText( "" );
        m_nuclideAgeEdit->disable();
      }else if( nuc->canObtainPromptEquilibrium() && showPromptOnly )
      {
        WString hlstr = PhysicalUnitsLocalized::printToBestTimeUnits(
                                                            5.0*nuc->promptEquilibriumHalfLife(),
                                                            2, SandiaDecay::second );
        m_nuclideAgeEdit->setText( hlstr );
      }else if( agestr.empty() )
      {
        string agstr;
        PeakDef::defaultDecayTime( nuc, &agstr );
        m_nuclideAgeEdit->setText( agstr );
      }else
      {
        m_nuclideAgeEdit->setText( agestr );
      }
    }else if( nuc && useCurrentAge )
    {
      const double hl = (nuc ? nuc->halfLife : -1.0);
      double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
      if( age > 100.0*nuc->halfLife || age < 0.0 )
        throw std::runtime_error( "" );
      m_nuclideAgeEdit->setText( agestr );
    }else if( nuc )
    {
      if( IsInf(nuc->halfLife) )
      {
        passMessage( WString::tr("dcw-nuc-is-stable").arg(isotopeLabel), WarningWidget::WarningMsgHigh );
      }else
      {
        passMessage( isotopeLabel + " is missing decay data", WarningWidget::WarningMsgHigh );
      }
      
      m_nuclideEdit->setText( "" );
      m_nuclideAgeEdit->setText( "" );
      m_nuclideAgeEdit->disable();
    }else if( !nuc )
    {
      m_nuclideAgeEdit->setText( "" );
      m_nuclideAgeEdit->disable();
    }//if( nuc && !IsInf(nuc->halfLife) ) / else
  }catch(...)
  {
    if( nuc )
    {
      string defagstr;
      PeakDef::defaultDecayTime( nuc, &defagstr );
      WString msg = WString::tr("dcw-made-age-reasonable")
                      .arg(nuc->symbol)
                      .arg(agestr)
                      .arg(defagstr);
      passMessage( msg, WarningWidget::WarningMsgLow );
      m_nuclideAgeEdit->setText( defagstr );
    }else
    {
      m_nuclideAgeEdit->setText( "0y" );
    }//if( nuc ) / else
  }//try / catch
  
  m_changed.emit();
}//void NuclideSourceEnterController::handleNuclideUserInput()

  
void NuclideSourceEnterController::handleAgeUserChange()
{
  if( !m_currentNuc )
  {
    m_nuclideAgeEdit->setText( "" );
    m_nuclideAgeEdit->disable();
    return;
  }//if( !m_currentNuc )
  
  
  //Validate the age is something reasonable - if not, change it.
  if( m_nuclideAgeEdit->text().empty() )
    m_nuclideAgeEdit->setText( "0y" );
  
  
  double age = std::numeric_limits<double>::max();
  string agestr = m_nuclideAgeEdit->text().toUTF8();
  
  try
  {
    const double hl = (m_currentNuc ? m_currentNuc->halfLife : -1.0);
    age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
  }catch(...) {}
  
  if( age > 50.0*m_currentNuc->halfLife || age < 0.0 )
  {
    passMessage( WString::tr("dcw-age-to-big").arg(agestr).arg(m_currentNuc->symbol),
                 WarningWidget::WarningMsgHigh );
    if( m_prevAgeTxt.count(m_currentNuc) )
      m_nuclideAgeEdit->setText( m_prevAgeTxt[m_currentNuc] );
  }//if( age > 50.0*m_currentNuc->halfLife )
  
  m_prevAgeTxt[m_currentNuc] = m_nuclideAgeEdit->text().toUTF8();
  
  m_changed.emit();
}//handleAgeUserChange()
  

Wt::WString NuclideSourceEnterController::nuclideAgeStr() const
{
  if( m_nuclideAgeEdit->isDisabled() )
    return WString();
  return m_nuclideAgeEdit->text();
}


double NuclideSourceEnterController::nuclideAge() const
{
  if( !m_currentNuc )
    throw runtime_error( "No current nuclide" );
  
  string agestr = m_nuclideAgeEdit->text().toUTF8();
  
  //m_nuclideAgeEdit should only be disabled if the isotope returns true for
  //  PeakDef::ageFitNotAllowed(nuc).
  if( m_nuclideAgeEdit->isDisabled() )
  {
    double ageval = 5.0*m_currentNuc->halfLife;
    if( m_currentNuc->canObtainPromptEquilibrium() )
      ageval = log(10000.0)/log(2.0) * m_currentNuc->promptEquilibriumHalfLife();
    agestr = PhysicalUnitsLocalized::printToBestTimeUnits( ageval );
    
//      cerr << "Niave Agestr=" << agestr << endl;
//      ageval = PeakDef::defaultDecayTime( m_currentNuc, &agestr );
//      cerr << "defaultDecayTime Agestr=" << agestr << endl;
    
//      double gregsage = 0.0;
//      EquilibriumType actual_type = getEquilibrium( m_currentNuc, gregsage, SecularEquilibrium );
//      cerr << "gregsage=" << PhysicalUnits::printToBestTimeUnits(gregsage) << endl;
    
//      cerr << "Agestr=" << agestr << endl;
  }//if( m_nuclideAgeEdit->isDisabled() )
  
  const double hl = (m_currentNuc ? m_currentNuc->halfLife : -1.0);
  const double age = PhysicalUnitsLocalized::stringToTimeDurationPossibleHalfLife( agestr, hl );
  if( age < 0 )
    throw runtime_error( "Invalid negative age" );
  
  return age;
}//double nuclideAge() const


const SandiaDecay::Nuclide *NuclideSourceEnterController::nuclide() const
{
  return m_currentNuc;
}


void NuclideSourceEnterController::setNuclideText( const std::string &txt )
{
  m_nuclideEdit->setText( WString::fromUTF8(txt) );
  handleNuclideUserInput();
}

void NuclideSourceEnterController::setNuclideAgeTxt( const std::string &txt )
{
  m_nuclideAgeEdit->setText( txt );
  handleAgeUserChange();
}


Wt::Signal<> &NuclideSourceEnterController::changed()
{
  return m_changed;
}




NuclideSourceEnter::NuclideSourceEnter( const bool showHalfLife, const bool showToolTips, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_nuclideEdit( nullptr ),
    m_nuclideAgeEdit( nullptr ),
    m_halfLifeTxt( nullptr ),
    m_controller( nullptr )
{
  InterSpec::instance()->useMessageResourceBundle( "NuclideSourceEnter" );
  
  wApp->useStyleSheet( "InterSpec_resources/NuclideSourceEnter.css" );
  addStyleClass( "NuclideSourceEnter" );
  
  WLabel *nucLabel = new WLabel( WString::tr("nuclide-label") );
  m_nuclideEdit = new WLineEdit();
  
  m_nuclideEdit->setMinimumSize( 30, WLength::Auto );
  nucLabel->setBuddy( m_nuclideEdit );
  
  WLabel *ageLabel = new WLabel( WString::tr("age-label") );
  m_nuclideAgeEdit = new WLineEdit();
  m_nuclideAgeEdit->setMinimumSize( 30, WLength::Auto );
  m_nuclideAgeEdit->setPlaceholderText( WString::tr("N/A") );
  ageLabel->setBuddy( m_nuclideAgeEdit );
  
  if( showHalfLife )
  {
    m_halfLifeTxt = new WText();
    m_halfLifeTxt->addStyleClass( "HLTxt" );
  }//if( showHalfLife )
  
#define USE_WLAYOUT_FOR_THIS 0
#if( USE_WLAYOUT_FOR_THIS )
  WGridLayout *layout = new WGridLayout( this );
  layout->setContentsMargins( 0, 0, 0, 0 );
  layout->addWidget( nucLabel, 0, 0, AlignMiddle );
  layout->addWidget( m_nuclideEdit, 0, 1 );
  layout->addWidget( ageLabel, 1, 0, AlignMiddle );
  layout->addWidget( m_nuclideAgeEdit, 1, 1 );
  if( m_halfLifeTxt )
    layout->addWidget( m_halfLifeTxt, 0, 2, AlignMiddle );
  layout->setColumnStretch( 1, 1 );
#else
  wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
  nucLabel->addStyleClass( "GridFirstRow GridFirstCol GridVertCenter" );
  addWidget( nucLabel );
  
  m_nuclideEdit->addStyleClass( "GridFirstRow GridSecondCol GridStretchCol" );
  addWidget( m_nuclideEdit );
  
  if( m_halfLifeTxt )
  {
    m_halfLifeTxt->addStyleClass( "GridFirstRow GridThirdCol GridVertCenter" );
    addWidget( m_halfLifeTxt );
  }//if( m_halfLifeTxt )
  
  ageLabel->addStyleClass( "GridSecondRow GridFirstCol GridVertCenter" );
  addWidget( ageLabel );
  
  m_nuclideAgeEdit->addStyleClass( "GridSecondRow GridSecondCol GridStretchCol" );
  addWidget( m_nuclideAgeEdit );
#endif
  
  HelpSystem::attachToolTipOn( {nucLabel, m_nuclideEdit},
                              WString::tr("dcw-tt-nuc-edit"), showToolTips );
  
  HelpSystem::attachToolTipOn( {ageLabel, m_nuclideAgeEdit},
                              WString::tr("dcw-tt-age-edit"), showToolTips );
  
  m_controller = new NuclideSourceEnterController( m_nuclideEdit, m_nuclideAgeEdit, m_halfLifeTxt, this );
}//NuclideSourceEnter constructor


NuclideSourceEnter::~NuclideSourceEnter()
{
}


void NuclideSourceEnter::setNuclideText( const string &txt )
{
  m_controller->setNuclideText( txt );
}


void NuclideSourceEnter::setNuclideAgeTxt( const string &txt )
{
  m_controller->setNuclideAgeTxt( txt );
}


double NuclideSourceEnter::getMaximumHalfLife( const SandiaDecay::Nuclide* nuclide )
{
  double halfLife = nuclide->halfLife;
  
  if( nuclide->halfLife < std::numeric_limits<double>::infinity() )
  {
    const vector<const SandiaDecay::Nuclide *> children = nuclide->descendants();
    for( size_t i = 1; i < children.size(); ++i )
    {
      const double childHalfLife = children[i]->halfLife;
      if( childHalfLife < std::numeric_limits<double>::infinity() )
      {
        if( childHalfLife > halfLife )
          halfLife = childHalfLife;
      }
    }
  }
  
  return halfLife;
}//getMaximumHalfLife(...)


void NuclideSourceEnter::getMaximumHalfLifeInDescension( const SandiaDecay::Nuclide * const nuclide, double parentHalfLife, double &halfLife )
{
  if( !nuclide )
  {
    halfLife = 0;
    return;
  }
  
  if (nuclide->halfLife < parentHalfLife)
  {
    if( nuclide->halfLife > halfLife )
      halfLife = nuclide->halfLife;
    
    for( size_t decN = 0; decN < nuclide->decaysToChildren.size(); ++decN )
    {
      const SandiaDecay::Transition * const transition = nuclide->decaysToChildren[decN];
      if( !transition )
        continue; // spontaneous fission
      
      const SandiaDecay::Nuclide * const child = transition->child;
      if( !child )
        continue; // not sure...
      
      const double childHalfLife = child->halfLife;
      if( childHalfLife < nuclide->halfLife && child->decayConstant() > 1e-50 )
      {
        getMaximumHalfLifeInDescension(child, parentHalfLife, halfLife);
      }
    }
  }
}//getMaximumHalfLifeInDescension(...)


bool NuclideSourceEnter::trySecularEquilibrium( const SandiaDecay::Nuclide * const nuclide, double &age )
{
  if( !nuclide )
    return false;
  
  double maxHalfLife = getMaximumHalfLife(nuclide);
  
  if( maxHalfLife < nuclide->halfLife )
  {
    age = log(10000.0)/log(2.0) * maxHalfLife; // 99.99% of longest-lived descendant activity reached
    return true;
  }
  return false;
}//trySecularEquilibrium(...)


bool NuclideSourceEnter::tryPromptEquilibrium( const SandiaDecay::Nuclide * const nuclide, double &age )
{
  if( !nuclide )
    return false;
  
  // find maximum half life in monotonically decreasing half lives
  double maxHalfLife = 0.0;
  for( size_t decN = 0; decN < nuclide->decaysToChildren.size(); ++decN )
  {
    const SandiaDecay::Transition * const transition = nuclide->decaysToChildren[decN];
    if( !transition )
      continue; // spontaneous fission
    
    const SandiaDecay::Nuclide * const child = transition->child;
    if( !child )
      continue; // not sure...
    
    if( child->decayConstant() > 1e-50 )
      getMaximumHalfLifeInDescension(child, nuclide->halfLife, maxHalfLife);
  }//for( size_t decN = 0; decN < nuclide->decaysToChildren.size(); ++decN )
  
  if( maxHalfLife < nuclide->halfLife )
  {
    age = log(10000.0) * maxHalfLife; // 99.99% of longest-lived descendants activity reached
    return true;
  }
  
  return false;
}//tryPromptEquilibrium(...)


NuclideSourceEnter::EquilibriumType NuclideSourceEnter::getEquilibrium( const SandiaDecay::Nuclide* nuclide, double &age, EquilibriumType eqCode )
{
  switch( eqCode )
  {
    case SecularEquilibrium:
      if( trySecularEquilibrium(nuclide, age) )
        return SecularEquilibrium;
      
      if( tryPromptEquilibrium(nuclide, age) )
        return PromptEquilibrium;
      
      age = 0.0;
      return NoEquilibrium;
      break;
      
    case PromptEquilibrium:
      if( tryPromptEquilibrium(nuclide, age) )
        return PromptEquilibrium;
      
      if( trySecularEquilibrium(nuclide, age) )
        return SecularEquilibrium;
      
      age = 0.0;
      return NoEquilibrium;
      break;
      
    default:
      return NoEquilibrium;
  }
}//getEquilibrium(...)


Wt::WString NuclideSourceEnter::nuclideAgeStr() const
{
  return m_controller->nuclideAgeStr();
}


double NuclideSourceEnter::nuclideAge() const
{
  return m_controller->nuclideAge();
}//double nuclideAge()


std::vector< std::pair<float,float> > NuclideSourceEnter::photonEnergyAndIntensity( const double activity ) const
{
  std::vector< std::pair<float,float> > answer;
  
  const SandiaDecay::Nuclide *nuc = m_controller->nuclide();
  
  if( !nuc || activity <= 0.0 )
    return answer;
  
  double age = 0.0;
  try
  {
    age = nuclideAge();
  }catch(...)
  {
    return answer;
  }
  
  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, activity, age );
  
  const auto results = mix.photons(0.0, SandiaDecay::NuclideMixture::OrderByEnergy);
  
  for( const SandiaDecay::EnergyRatePair &aep : results )
    answer.push_back( make_pair( static_cast<float>(aep.energy), static_cast<float>(aep.numPerSecond) ) );
  
  return answer;
}//std::vector< std::pair<float,float> > photonEnergyAndIntensity() const


Wt::Signal<> &NuclideSourceEnter::changed()
{
  return m_controller->changed();
}


const SandiaDecay::Nuclide *NuclideSourceEnter::nuclide() const
{
  return m_controller->nuclide();
}

