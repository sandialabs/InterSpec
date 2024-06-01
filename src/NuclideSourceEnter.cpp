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


NuclideSourceEnter::NuclideSourceEnter( const bool showToolTips, Wt::WContainerWidget *parent )
   : WContainerWidget( parent ),
     m_nuclideEdit( 0 ),
     m_nuclideAgeEdit( 0 ),
     m_halfLifeTxt( 0 ),
     m_currentNuc( 0 ),
     m_changed( this )
  {
    InterSpec::instance()->useMessageResourceBundle( "NuclideSourceEnter" );

    wApp->useStyleSheet( "InterSpec_resources/NuclideSourceEnter.css" );
    addStyleClass( "NuclideSourceEnter" );
    
    WGridLayout *layout = new WGridLayout( this );
//    layout->setHorizontalSpacing( 0 );
//    layout->setVerticalSpacing( 0 );
    layout->setContentsMargins( 0, 0, 0, 0 );
    
    WLabel *label = new WLabel( WString("{1}:").arg(WString::tr("Nuclide")) );
    layout->addWidget( label, 0, 0, AlignMiddle );
    label->addStyleClass( "DoseFieldLabel" );
    m_nuclideEdit = new WLineEdit();
    m_nuclideEdit->setAutoComplete( false );
    m_nuclideEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_nuclideEdit->setAttributeValue( "autocorrect", "off" );
    m_nuclideEdit->setAttributeValue( "spellcheck", "off" );
#endif
    layout->addWidget( m_nuclideEdit, 0, 1 );
    m_nuclideEdit->addStyleClass( "DoseEnterTxt" );
    m_nuclideEdit->setAttributeValue( "spellcheck", "false" );
    m_nuclideEdit->setMinimumSize( 30, WLength::Auto );
    label->setBuddy( m_nuclideEdit );
    
    label = new WLabel( WString("{1}:").arg(WString::tr("Age")) );
    WLabel *ageLabel = label;
    layout->addWidget( label, 1, 0, AlignMiddle );
    label->addStyleClass( "DoseFieldLabel" );
    m_nuclideAgeEdit = new WLineEdit();
    m_nuclideAgeEdit->setMinimumSize( 30, WLength::Auto );
    
    m_nuclideAgeEdit->setAutoComplete( false );
    m_nuclideAgeEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_nuclideAgeEdit->setAttributeValue( "autocorrect", "off" );
    m_nuclideAgeEdit->setAttributeValue( "spellcheck", "off" );
#endif
    m_nuclideAgeEdit->setPlaceholderText( WString::tr("N/A") );
    layout->addWidget( m_nuclideAgeEdit, 1, 1 );
    label->setBuddy( m_nuclideAgeEdit );
    
    m_halfLifeTxt = new WText();
    m_halfLifeTxt->addStyleClass( "DoseHLTxt" );
    layout->addWidget( m_halfLifeTxt, 0, 2, AlignMiddle );
    
//    layout->setColumnStretch( 0, 1 );
    layout->setColumnStretch( 1, 1 );
    
    string replacerJs, matcherJs;
    PhotopeakDelegate::EditWidget::replacerJs( replacerJs );
    PhotopeakDelegate::EditWidget::nuclideNameMatcherJs( matcherJs );
    
    WSuggestionPopup *suggestions = new WSuggestionPopup( matcherJs, replacerJs );
#if( WT_VERSION < 0x3070000 ) //I'm not sure what version of Wt "wtNoReparent" went away.
    suggestions->setJavaScriptMember("wtNoReparent", "true");
#endif
    
    suggestions->setMaximumSize( WLength::Auto, WLength(15, WLength::FontEm) );
    suggestions->setWidth( WLength(70, Wt::WLength::Unit::Pixel) );
    suggestions->forEdit( m_nuclideEdit,
                         WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );
    
    
    IsotopeNameFilterModel *filterModel = new IsotopeNameFilterModel( this );
    
    filterModel->excludeNuclides( false );
    filterModel->excludeXrays( true );
    filterModel->excludeEscapes( true );
    filterModel->excludeReactions( true );
    
    filterModel->filter( "" );
    suggestions->setFilterLength( -1 );
    suggestions->setModel( filterModel );
    suggestions->filterModel().connect( filterModel, &IsotopeNameFilterModel::filter );
    
    m_nuclideEdit->changed().connect( this, &NuclideSourceEnter::handleNuclideUserInput );
//    m_nuclideEdit->enterPressed().connect( this, &NuclideSourceEnter::handleNuclideUserInput );
//    m_nuclideEdit->blurred().connect( this, &NuclideSourceEnter::handleNuclideUserInput );

    HelpSystem::attachToolTipOn( m_nuclideEdit, WString::tr("dcw-tt-nuc-edit"), showToolTips );
    
    
    WRegExpValidator *validator = new WRegExpValidator( PhysicalUnitsLocalized::timeDurationHalfLiveOptionalRegex(), this );
    validator->setFlags(Wt::MatchCaseInsensitive);
    m_nuclideAgeEdit->setValidator(validator);
    
    m_nuclideAgeEdit->changed().connect( this, &NuclideSourceEnter::handleAgeUserChange );
    m_nuclideAgeEdit->blurred().connect( this, &NuclideSourceEnter::handleAgeUserChange );
    m_nuclideAgeEdit->enterPressed().connect( this, &NuclideSourceEnter::handleAgeUserChange );
    
    HelpSystem::attachToolTipOn( {ageLabel, m_nuclideAgeEdit}, 
                                WString::tr("dcw-tt-age-edit"), showToolTips );
  }//NuclideSourceEnter constructor
  
  
  NuclideSourceEnter::~NuclideSourceEnter()
  {
  }
  
  void NuclideSourceEnter::setNuclideText( const string &txt )
  {
    m_nuclideEdit->setText( WString::fromUTF8(txt) );
    handleNuclideUserInput();
  }
  
  void NuclideSourceEnter::setNuclideAgeTxt( const string &txt )
  {
    m_nuclideAgeEdit->setText( txt );
    handleAgeUserChange();
  }
  
  void NuclideSourceEnter::handleNuclideUserInput()
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const string isotopeLabel = m_nuclideEdit->text().toUTF8();
    const SandiaDecay::Nuclide *nuc = db->nuclide( isotopeLabel );
    
    if( nuc == m_currentNuc )
      return;

    m_currentNuc = nuc;
    
    WString hlstr;
    if( nuc )
      hlstr = PhysicalUnitsLocalized::printToBestTimeUnits( nuc->halfLife, 2, SandiaDecay::second );
    m_halfLifeTxt->setText( WString("{1}={2}").arg(WString::tr("T1/2")).arg(hlstr) );

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
  }//void handleNuclideUserInput()
  
  
  void NuclideSourceEnter::handleAgeUserChange()
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
  }//void handleAgeUserChange()
  
  
  
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
  }
  
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
  }
  
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
  }
  
  
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
  }
  
  Wt::WString NuclideSourceEnter::nuclideAgeStr() const
  {
    if( m_nuclideAgeEdit->isDisabled() )
      return WString();
    return m_nuclideAgeEdit->text();
  }
  
  double NuclideSourceEnter::nuclideAge() const
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
  }//double nuclideAge()
  
  /** 
   * \returns empty results if no valid isotope, an invalid age, or negative or 
   *          zero activity. Other wise returns <energy,gamma/sec> pairs.
   */
  std::vector< std::pair<float,float> > NuclideSourceEnter::photonEnergyAndIntensity( const double activity ) const
  {
    std::vector< std::pair<float,float> > answer;
    
    if( !m_currentNuc || activity <= 0.0 )
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
    mix.addAgedNuclideByActivity( m_currentNuc, activity, age );
    
    const auto results = mix.photons(0.0, SandiaDecay::NuclideMixture::OrderByEnergy);
    
    for( const SandiaDecay::EnergyRatePair &aep : results )
      answer.push_back( make_pair( static_cast<float>(aep.energy), static_cast<float>(aep.numPerSecond) ) );
    
    return answer;
  }//std::vector< std::pair<float,float> > photonEnergyAndIntensity() const
  
  Wt::Signal<> &NuclideSourceEnter::changed()
  {
    return m_changed;
  }
  
  const SandiaDecay::Nuclide *NuclideSourceEnter::nuclide() const
  { 
    return m_currentNuc; 
  }
  
