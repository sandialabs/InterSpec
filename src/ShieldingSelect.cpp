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
#include <cstdlib>
#include <sstream>
#include <iostream>

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WImage>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WGroupBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WDoubleValidator>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>


#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"


#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "SpecUtils/RapidXmlUtils.hpp"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace Wt;
using namespace std;

using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::TraceActivityType;

const int ShieldingSelect::sm_xmlSerializationMajorVersion = 0;
const int ShieldingSelect::sm_xmlSerializationMinorVersion = 1;


namespace
{
double distance_of_input_text( const WLineEdit *edit )
{
  const string text = edit->text().toUTF8();
  
  try
  {
    return PhysicalUnits::stringToDistance( text );
  }catch( std::exception & )
  {
  }
  
  throw runtime_error( "Error converting '" + text + "' to a distance" );
  
  return 0.0;
}//distance_of_input_text(...)

}//namespace


/**
 
 */
class TraceSrcDisplay : public WGroupBox
{
  const SandiaDecay::Nuclide *m_currentNuclide;
  
  /** The activity the user has entered.  E.g., maybe total, but maybe activity per cm3 or per gram. */
  double m_currentDisplayActivity;
  
  /** The current total activity of the shielding. */
  double m_currentTotalActivity;
  
  ShieldingSelect *m_parent;
  Wt::WComboBox *m_isoSelect;
  Wt::WLineEdit *m_activityInput;
  Wt::WComboBox *m_activityType;
  Wt::WCheckBox *m_allowFitting;
  Wt::Signal<const SandiaDecay::Nuclide *,double> m_activityUpdated;
  Wt::Signal<TraceSrcDisplay *, const SandiaDecay::Nuclide * /* old nuclide */> m_nucChangedSignal;
  
  
public:
  TraceSrcDisplay( ShieldingSelect *parent )
  : WGroupBox( "Trace Source", parent->m_traceSources ),
    m_currentNuclide( nullptr ),
    m_currentDisplayActivity( 0.0 ),
    m_currentTotalActivity( 0.0 ),
    m_parent( parent ),
    m_isoSelect( nullptr ),
    m_activityInput( nullptr ),
    m_activityType( nullptr ),
    m_allowFitting( nullptr ),
    m_activityUpdated( this ),
    m_nucChangedSignal( this )
  {
    const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    assert( !parent->isGenericMaterial() );
    addStyleClass( "TraceSrcDisplay" );
    
    WPushButton *closeIcon = new WPushButton( this );
    closeIcon->addStyleClass( "closeicon-wtdefault ThirdCol FirstRow" );
    closeIcon->clicked().connect( boost::bind( &ShieldingSelect::removeTraceSourceWidget, m_parent, this) );
    closeIcon->setToolTip( "Remove this trace source." );
    
    
    WLabel *label = new WLabel( "Nuclide", this );
    label->addStyleClass( "FirstCol SecondRow" );
    
    m_isoSelect = new WComboBox( this );
    m_isoSelect->addItem( "Select" );
    m_isoSelect->activated().connect( this, &TraceSrcDisplay::handleUserNuclideChange );
    m_isoSelect->addStyleClass( "SecondCol SecondRow SpanTwoCol" );
    label->setBuddy( m_isoSelect );
    
    
    label = new WLabel( "Activity", this );
    label->addStyleClass( "FirstCol ThirdRow" );
    
    
    m_activityInput = new WLineEdit( this );
    m_activityInput->addStyleClass( "SecondCol ThirdRow" );
    m_activityInput->setAutoComplete( false );
    
    WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityRegex, m_activityInput );
    val->setFlags( Wt::MatchCaseInsensitive );
    m_activityInput->setValidator( val );
    m_activityInput->changed().connect( this, &TraceSrcDisplay::handleUserActivityChange );
    m_activityInput->enterPressed().connect( this, &TraceSrcDisplay::handleUserActivityChange );
    m_activityInput->setText( (useBq ? "37 MBq" : "1 mCi") );
    m_currentTotalActivity = m_currentDisplayActivity = 0.001*PhysicalUnits::ci;
    label->setBuddy( m_activityInput );
    
    m_activityType = new WComboBox( this );
    m_activityType->activated().connect( this, &TraceSrcDisplay::handleUserChangeActivityType );
    m_activityType->addStyleClass( "ThirdCol ThirdRow" );
    
    m_allowFitting = new WCheckBox( "Fit activity value", this );
    m_allowFitting->addStyleClass( "SecondCol FourthRow SpanTwoCol" );
    m_allowFitting->changed().connect( this, &TraceSrcDisplay::handleUserChangeAllowFit );
      
    updateAvailableActivityTypes();
    updateAvailableIsotopes();
  }//TraceSrcDisplay constructor
  
  void serialize( rapidxml::xml_node<char> * const parent_node ) const
  {
    rapidxml::xml_document<char> *doc = parent_node->document();
    
    rapidxml::xml_node<char> * const base_node
                                      = doc->allocate_node( rapidxml::node_element, "TraceSource" );
    parent_node->append_node( base_node );
    
    const char *value = m_currentNuclide ? m_currentNuclide->symbol.c_str() : "None";
    rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, "Nuclide", value );
    base_node->append_node( node );
    
    
    value = doc->allocate_string( m_activityInput->text().toUTF8().c_str() );
    node = doc->allocate_node( rapidxml::node_element, "DisplayActivity", value );
    base_node->append_node( node );
    
    // Put the total activity as an attribute, just to make a little more human readable/checkable
    const string total_act = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity );
    value = doc->allocate_string( total_act.c_str() );
    rapidxml::xml_attribute<char> *attr = doc->allocate_attribute( "TotalActivity", value );
    node->append_attribute( attr );
    
    const TraceActivityType type = TraceActivityType( m_activityType->currentIndex() );
    value = GammaInteractionCalc::to_str( type );
    
    node = doc->allocate_node( rapidxml::node_element, "TraceActivityType", value );
    base_node->append_node( node );
      
    value = m_allowFitting->isChecked() ? "1" : "0";
    node = doc->allocate_node( rapidxml::node_element, "AllowFitting", value );
    base_node->append_node( node );
  }//void serialize( rapidxml::xml_node<> *parent )
  
  
  void deSerialize( const rapidxml::xml_node<char> * const trace_node )
  {
    //Should check that the model has the nuclide, otherwise dont select a nuclide
    //For development builds should check nuclide is already marked as a trace source in the source fitting model.activityUpdated()
    
    if( !trace_node || !trace_node->value()
       || rapidxml::internal::compare(trace_node->value(), trace_node->value_size(),"TraceSource", 11, true) )
      throw runtime_error( "TraceSrcDisplay::deSerialize: called with invalid node" );
      
    const rapidxml::xml_node<char> *nuc_node = trace_node->first_node( "Nuclide" );
    const rapidxml::xml_node<char> *disp_act_node = trace_node->first_node( "DisplayActivity" );
    const rapidxml::xml_node<char> *act_type_node = trace_node->first_node( "TraceActivityType" );
    const rapidxml::xml_node<char> *allow_fit_node = trace_node->first_node( "AllowFitting" );
    
    if( !nuc_node || !disp_act_node || !act_type_node || !allow_fit_node )
      throw runtime_error( "TraceSrcDisplay::deSerialize: missing node" );

    updateAvailableIsotopes();
    updateAvailableActivityTypes();
    
    const string nuclide = SpecUtils::xml_value_str(nuc_node);
    
    int nucSelectIndex = -1;
    for( int i = 0; i < m_isoSelect->count(); ++i )
    {
      const string nuctxt = m_isoSelect->itemText(i).toUTF8();
      if( nuctxt == nuclide )
      {
        nucSelectIndex = i;
        break;
      }
    }//for( int i = 0; i < m_isoSelect->count(); ++i )
    
    if( nucSelectIndex < 0 )
    {
      m_isoSelect->setCurrentIndex(0); //"Select"
      cerr << "TraceSrcDisplay::deSerialize: Failed to match the expected nuclide ('"
           << nuclide << "') to avaiable ones."
           << endl;
      // Should we throw an error here or something???
    }else
    {
      m_isoSelect->setCurrentIndex(nucSelectIndex);
    }
    
    handleUserNuclideChange();
    
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    string acttxt = SpecUtils::xml_value_str(disp_act_node);
    
    // Check activity text is actually valid, and convert to users current prefered units
    try
    {
      const double activity = PhysicalUnits::stringToActivity(acttxt);
      const bool isInBq = SpecUtils::icontains(acttxt, "bq");
      if( isInBq == useCi )
        acttxt = PhysicalUnits::printToBestActivityUnits( activity, 6, useCi );
    }catch(...)
    {
      cerr << "TraceSrcDisplay::deSerialize: Activity from XML ('" << acttxt << "') is invalid.\n";
      acttxt = (useCi ? "0 uCi" : "0 bq");
      // Should we throw an error here or something???
    }
    
    m_activityInput->setText( WString::fromUTF8(acttxt) );
    
    const string act_type = SpecUtils::xml_value_str(act_type_node);
    
    TraceActivityType type = TraceActivityType::NumTraceActivityType;
    for( TraceActivityType t = TraceActivityType(0);
        t != TraceActivityType::NumTraceActivityType;
        t = TraceActivityType(static_cast<int>(t) + 1) )
    {
      if( act_type == GammaInteractionCalc::to_str(t) )
      {
        type = t;
        break;
      }
    }//for( loop over TraceActivityTypes )
    
    if( type == TraceActivityType::NumTraceActivityType )
    {
      cerr << "TraceSrcDisplay::deSerialize: Trace activity type in XML ('" << act_type << "'),"
           << " is invalid, setting to toal activity." << endl;
      type = TraceActivityType::TotalActivity;
    }
    
    m_activityType->setCurrentIndex( static_cast<int>(type) );
    
    const string do_fit_xml = SpecUtils::xml_value_str(allow_fit_node);
    const bool allow_fit = ((do_fit_xml == "1") || SpecUtils::iequals_ascii(do_fit_xml,"true"));
    m_allowFitting->setChecked( allow_fit );
    
    updateTotalActivityFromDisplayActivity();
    handleUserActivityChange();
    handleUserChangeAllowFit();
  }//void deSerialize( rapidxml::xml_node<> *trace_src_node )
  
  
  bool allowFittingActivity() const
  {
    const bool allow = m_allowFitting->isChecked();;
    return allow;
  }
  
  void modelSourceAdded( const SandiaDecay::Nuclide * const nuc )
  {
    if( !nuc )
      return;
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      return;
    
    bool hasNuc = false;
    for( int index = 1; !hasNuc && (index < m_isoSelect->count()); ++index )
    {
      const string label = m_isoSelect->itemText(index).toUTF8();
      const SandiaDecay::Nuclide * const thisNuc = db->nuclide( label );
      
      hasNuc = (hasNuc ||(thisNuc == nuc));
    }//for( loop over current nuclides )
    
    if( hasNuc )
      return;
    
    m_isoSelect->addItem( nuc->symbol );
  }//void modelSourceAdded( const SandiaDecay::Nuclide * )
  
  
  void modelSourceRemoved( const SandiaDecay::Nuclide * const nuc )
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      return;
    
    if( nuc == m_currentNuclide )
    {
      const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      m_isoSelect->removeItem( m_isoSelect->currentIndex() );
      m_isoSelect->setCurrentIndex( 0 );
      m_currentNuclide = nullptr;
      m_currentTotalActivity = m_currentDisplayActivity = 0.0;
      m_activityInput->setText( (useCi ? "0 uCi" : "0 bq") );
      
      // TODO: do we need to emit that we are removing this trace source? If we're here, the model
      //       already knows about this
    }else
    {
      const int nucIndex = m_isoSelect->findText( nuc->symbol );
      if( nucIndex >= 1 )
      {
        // zeroth entry should always be "Select"
        const WString currentText = m_isoSelect->currentText();
        m_isoSelect->removeItem( nucIndex );
        m_isoSelect->setCurrentIndex( m_isoSelect->findText( currentText ) );
      }//if( nucIndex >= 1 )
    }// if( nuc == m_currentNuclide ) / else
  }//void modelSourceRemoved( const SandiaDecay::Nuclide * )
  
  
  void handleUserActivityChange()
  {
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    if( m_parent->isGenericMaterial() )
    {
      assert( 0 );
      cerr << "handleUserActivityChange() called for a generic material - shouldnt happen" << endl;
      return;
    }
    
    assert( m_parent->m_shieldSrcDisp );
    
    const double shieldVolume = m_parent->shieldingVolume();
    const double shieldMass = m_parent->shieldingMass();
    
    const double shieldVolumeCm3 = shieldVolume / PhysicalUnits::cm3;
    const double shieldMassGram = shieldMass / PhysicalUnits::gram;
    
    const string userActTxt = m_activityInput->text().toUTF8();
    const TraceActivityType type = TraceActivityType( m_activityType->currentIndex() );
    
    try
    {
      m_currentDisplayActivity = PhysicalUnits::stringToActivity( userActTxt );
      
      switch( type )
      {
        case TraceActivityType::TotalActivity:
          m_currentTotalActivity = m_currentDisplayActivity;
          break;
          
        case TraceActivityType::ActivityPerCm3:
          m_currentTotalActivity = m_currentDisplayActivity * shieldVolumeCm3;
          break;
          
        case TraceActivityType::ActivityPerGram:
          m_currentTotalActivity = m_currentDisplayActivity * shieldMassGram;
          break;
        
        case TraceActivityType::NumTraceActivityType:
          assert( 0 );
          m_activityInput->setText( useCi ? "0 uCi" : "0 bq");
          m_currentTotalActivity = m_currentDisplayActivity = 0.0;
          break;
      }//switch( type )
    }catch( std::exception & )
    {
      cerr << "User entered activity '" << userActTxt << "' is invalid - reverting" << endl;
      
      string txt;
      switch( type )
      {
        case TraceActivityType::TotalActivity:
          txt = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity, 3, useCi );
          break;
          
        case TraceActivityType::ActivityPerCm3:
          if( shieldVolumeCm3 <= FLT_EPSILON )
          {
            m_currentTotalActivity = m_currentDisplayActivity = 0.0;
            txt = (useCi ? "0 uCi" : "0 bq");
          }else
          {
            txt = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity/shieldVolumeCm3, 3, useCi );
          }
          break;
          
        case TraceActivityType::ActivityPerGram:
          if( shieldMassGram <= FLT_EPSILON )
          {
            m_currentTotalActivity = m_currentDisplayActivity = 0.0;
            txt = (useCi ? "0 uCi" : "0 bq");
          }else
          {
            txt = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity/shieldMassGram, 3, useCi );
          }
          break;
        
        case TraceActivityType::NumTraceActivityType:
          assert( 0 );
          txt = (useCi ? "0 uCi" : "0 bq");
          break;
      }//switch( type )
      
      m_activityInput->setText( txt );
    }// try / catch
    
    if( m_currentNuclide )
      m_activityUpdated.emit( m_currentNuclide, m_currentTotalActivity );
  }//void handleUserActivityChange()
  
  /** Returns total activity, or activity per gram, or activity per cm3, based on user input*/
  double displayActivity()
  {
    return m_currentDisplayActivity;
  }//double displayActivity()
  
  
  TraceActivityType activityType() const
  {
    const int currentIndex = m_activityType->currentIndex();
    if( (currentIndex < 0)
       || (currentIndex >= static_cast<int>(TraceActivityType::NumTraceActivityType)) )
      return TraceActivityType::NumTraceActivityType;
  
    return static_cast<TraceActivityType>( currentIndex );
  }
  
  
  void updateAvailableActivityTypes()
  {
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    const int previous = m_activityType->currentIndex();
    assert( previous < static_cast<int>(TraceActivityType::NumTraceActivityType) );
    
    m_activityType->clear();
    
    if( m_parent->isGenericMaterial() )
      return; //Probably wont ever get here
    
    std::shared_ptr<Material> material = m_parent->material();
    if( !material )
      return;
    
    // If less dense than 1% of nominal air, then dont allow activity per gram
    const double minDensity = 0.000013*PhysicalUnits::g/PhysicalUnits::cm3;
    const bool allowActPerGram = (material->density > minDensity);
    
    for( TraceActivityType type = TraceActivityType(0);
        type < TraceActivityType::NumTraceActivityType;
        type = TraceActivityType(static_cast<int>(type) + 1) )
    {
      switch( type )
      {
        case TraceActivityType::TotalActivity:
          m_activityType->addItem( "Total" );
          break;
          
        case TraceActivityType::ActivityPerCm3:
          m_activityType->addItem( "per cm^3" );
          break;
          
        case TraceActivityType::ActivityPerGram:
          if( allowActPerGram )
            m_activityType->addItem( "per gram" );
          break;
          
        case TraceActivityType::NumTraceActivityType:
          break;
      }//switch( type )
    }//for( loop over TraceActivityTypes )
    
    bool changedIndex = false;
    int currentIndex = 0;
    if( (previous >= 0)
        && (allowActPerGram
            || (previous != static_cast<int>(TraceActivityType::ActivityPerGram))) )
    {
      changedIndex = true;
      currentIndex = previous;
    }
    
    m_activityType->setCurrentIndex( currentIndex );
    

    if( changedIndex )
    {
      // If we are fitting for shielding thickness, then we probably dont also want to fit
      //  for activity per cm3 or per gram
      //  TODO: investigate how well fitting activity per gram or cm3 works when fitting thickness
      //if( m_parent->fitThickness()
      //   && ((currentIndex == static_cast<int>(TraceActivityType::ActivityPerCm3))
      //     || (currentIndex == static_cast<int>(TraceActivityType::ActivityPerGram))) )
      //{
      //  m_allowFitting->setChecked( false );
      //}
      
      // We will keep total activity the same, but update the display
      updateDispActivityFromTotalActivity();
    }//if( changedIndex )
    
  }//void updateAvailableActivityTypes()
  
  
  /** We will keep total activity the same, but update the display activity based on current value of m_activityType. */
  void updateDispActivityFromTotalActivity()
  {
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    const int currentIndex = m_activityType->currentIndex();
    
    const double shieldVolume = m_parent->shieldingVolume();
    const double shieldMass = m_parent->shieldingMass();
    
    const double shieldVolumeCm3 = shieldVolume / PhysicalUnits::cm3;
    const double shieldMassGram = shieldMass / PhysicalUnits::gram;
  
    m_currentDisplayActivity = 0.0;
    
    switch( TraceActivityType(currentIndex) )
    {
      case TraceActivityType::TotalActivity:
        m_currentDisplayActivity = m_currentTotalActivity;
        break;
        
      case TraceActivityType::ActivityPerCm3:
        if( shieldVolumeCm3 <= FLT_EPSILON )
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        else
          m_currentDisplayActivity = m_currentTotalActivity / shieldVolumeCm3;
        break;
        
      case TraceActivityType::ActivityPerGram:
        if( shieldMassGram <= FLT_EPSILON )
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        else
          m_currentDisplayActivity = m_currentTotalActivity / shieldMassGram;
        break;
      
      case TraceActivityType::NumTraceActivityType:
        assert( 0 );
        m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        break;
    }//switch( type )
    
    const string actTxt = PhysicalUnits::printToBestActivityUnits( m_currentDisplayActivity, 4, useCi );
    m_activityInput->setText( actTxt );
  }//void updateDispActivityFromTotalActivity()
  
  
  void updateTotalActivityFromDisplayActivity()
  {
    const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    
    if( !m_currentNuclide )
    {
      m_currentTotalActivity = m_currentDisplayActivity = 0.0;
      m_activityInput->setText( (useCi ? "0 uCi" : "0 bq") );
      return;
    }//if( !m_currentNuclide )
    
    
    const int currentIndex = m_activityType->currentIndex();
    
    const double shieldVolume = m_parent->shieldingVolume();
    const double shieldMass = m_parent->shieldingMass();
    
    const double shieldVolumeCm3 = shieldVolume / PhysicalUnits::cm3;
    const double shieldMassGram = shieldMass / PhysicalUnits::gram;
  
    
    m_currentTotalActivity = 0.0;
    
    switch( TraceActivityType(currentIndex) )
    {
      case TraceActivityType::TotalActivity:
        m_currentTotalActivity = m_currentDisplayActivity;
        break;
        
      case TraceActivityType::ActivityPerCm3:
        if( shieldVolumeCm3 <= FLT_EPSILON )
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        else
          m_currentTotalActivity = m_currentDisplayActivity * shieldVolumeCm3;
        break;
        
      case TraceActivityType::ActivityPerGram:
        if( shieldMassGram <= FLT_EPSILON )
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        else
          m_currentTotalActivity = m_currentDisplayActivity * shieldMassGram;
        break;
      
      case TraceActivityType::NumTraceActivityType:
        assert( 0 );
        m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        break;
    }//switch( type )
  }//void updateTotalActivityFromDisplayActivity()
  
  
  void updateForMaterialChange()
  {
    shared_ptr<Material> mat = m_parent->material();
    if( !mat )
    {
      const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      m_currentDisplayActivity = m_currentTotalActivity = 0.0;
      m_activityInput->setText( (useBq ? "0 Bq" : "0 uCi") );
      m_isoSelect->setCurrentIndex( 0 );
      m_allowFitting->setChecked( false );
      
      if( m_currentNuclide )
      {
        const SandiaDecay::Nuclide * const oldNuc = m_currentNuclide;
        m_currentNuclide = nullptr;
        m_nucChangedSignal.emit( this, oldNuc );
      }
    }//if( !mat )
    
    updateAvailableActivityTypes();
    updateAvailableIsotopes();
    updateTotalActivityFromDisplayActivity();
  }//void updateForMaterialChange()
  
  
  Wt::Signal<const SandiaDecay::Nuclide *,double> &activityUpdated()
  {
    return m_activityUpdated;
  }
  
  Wt::Signal<TraceSrcDisplay *, const SandiaDecay::Nuclide *> &nucChangedSignal()
  {
    return m_nucChangedSignal;
  }
  
  const SandiaDecay::Nuclide *nuclide() const
  {
    return m_currentNuclide;
  }//nuclide()
  
  double totalActivity() const
  {
    return m_currentTotalActivity;
  }
  
  double displayActivity() const
  {
    return m_currentDisplayActivity;
  }
  
  void setTraceSourceTotalActivity( const double total_activity )
  {
    assert( m_currentNuclide );
    if( !m_currentNuclide )
      throw runtime_error( "setTraceSourceTotalActivity: no current nuclide" );
    
    m_currentTotalActivity = total_activity;
    updateDispActivityFromTotalActivity();
    
    m_activityUpdated.emit( m_currentNuclide, m_currentTotalActivity );
  }//void setTraceSourceTotalActivity(...)
  
  
  void handleUserNuclideChange()
  {
    const SandiaDecay::Nuclide * const oldNuc = m_currentNuclide;
    
    const string nuctxt = m_isoSelect->currentText().toUTF8();
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    if( !db )
      m_currentNuclide = nullptr;
    else
      m_currentNuclide = db->nuclide( nuctxt );
    
    // Start with the models activity.
    SourceFitModel *srcmodel = m_parent->m_sourceModel;
    if( m_currentNuclide && srcmodel )
    {
      // This next line will throw exception if invalid nuclide ...
      //  but that should never happen, right?
      const int nucnum = srcmodel->nuclideIndex(m_currentNuclide);
      double act = srcmodel->activity(nucnum);
      const TraceActivityType type = TraceActivityType( m_activityType->currentIndex() );
      
      switch( type )
      {
        case GammaInteractionCalc::TraceActivityType::TotalActivity:
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
          act /= (m_parent->shieldingVolume() / PhysicalUnits::cm3);
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
          act /= (m_parent->shieldingMass() / PhysicalUnits::gram);
          break;
          
        case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
          assert( 0 );
          break;
      }//switch( type )
      
      m_currentDisplayActivity = act;
      
      const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      const string actstr = PhysicalUnits::printToBestActivityUnits(act,4,useCi);
      m_activityInput->setValueText( WString::fromUTF8(actstr) );
      
      updateTotalActivityFromDisplayActivity();
    }//if( srcmodel )
    
    handleUserChangeAllowFit();
    
    if( m_currentNuclide != oldNuc )
      m_nucChangedSignal.emit( this, oldNuc );
  }//void handleUserNuclideChange()
  
  
  void handleUserChangeActivityType()
  {
    updateDispActivityFromTotalActivity();
  }//handleUserChangeActivityType()
  
  
  void handleUserChangeAllowFit()
  {
    //We dont actually display the fit checkbox in the model table for trace sources, and we only
    //  consult this widget if we are fitting for it, but lets still keep the model in sync, for
    //  the moment (although this "in sync" only goes one direction; if you change the model value
    //  it wont change this widgets value, atm).
    if( !m_currentNuclide )
      return;
    
    SourceFitModel *srcmodel = m_parent->m_sourceModel;
    assert( srcmodel );
    if( !srcmodel )
      return;
    
    try
    {
      WModelIndex index = srcmodel->index( m_currentNuclide, SourceFitModel::Columns::kFitActivity );
      const bool doFit = m_allowFitting->isChecked();
      srcmodel->setData( index, doFit, Wt::CheckStateRole );
    }catch( std::exception &e )
    {
      assert( 0 ); // shouldnt ever happen, right???
      cerr << "handleUserChangeAllowFit(): Unexpected exception getting nuclide index from"
           << " SourceFitModel" << endl;
      return;
    }//try / catch
  }//void handleUserChangeAllowFit()
  
  /** Doesn not cause change signals to be emitted */
  void unselectNuclide( const SandiaDecay::Nuclide *nuc )
  {
    if( !nuc )
      return;
    
    if( nuc == m_currentNuclide )
    {
      m_isoSelect->setCurrentIndex( 0 );
      m_currentNuclide = nullptr;
    }
  }//void unselectNuclide( const SandiaDecay::Nuclide *nuc )
  
  /** Updates which nuclides the user can choose.
   Still hashing out the behaviour I want, but currently, if calling this
   
   */
  void updateAvailableIsotopes()
  {
    assert( m_parent );
    
    SourceFitModel *model = m_parent->m_sourceModel;
    assert( model );
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    assert( db );
    
    const WString currentNucTxt = m_isoSelect->currentText();
    
    const SandiaDecay::Nuclide *currentNuc = nullptr;
    vector<const SandiaDecay::Nuclide *> currentNuclides;
    for( int index = 1; index < m_isoSelect->count(); ++index )
    {
      const string nucName = m_isoSelect->itemText(index).toUTF8();
      const SandiaDecay::Nuclide *nuc = db->nuclide( nucName );
      if( !nuc )
      {
#if( PERFORM_DEVELOPER_CHECKS )
        const string msg = "Unexpected non-nuclide entry in m_isoSelect: " + nucName;
        log_developer_error( __func__, msg.c_str() );
#endif
        continue;
      }//if( !nuc )
      
      if( index == m_isoSelect->currentIndex() )
        currentNuc = nuc;
      
      currentNuclides.push_back( nuc );
    }//for( int index = 1; index < m_isoSelect->count(); ++index )
    
    
    const int numNucs = model->numNuclides();

    m_isoSelect->clear();
    m_isoSelect->addItem( "Select" );
    
    int indexToSelect = -1, numNucAdded = 0;
    for( int index = 0; index < numNucs; ++index )
    {
      const SandiaDecay::Nuclide *nuc = model->nuclide( index );
      if( !nuc )
      {
        cerr << "Unexpected null nuclide!!!" << endl;
        continue;
      }
      
      const double activity = model->activity( index );
      //const bool fitActivity = model->fitActivity( index );
      
      if( nuc == currentNuc )
      {
        indexToSelect = index + 1;
      }else
      {
        // Next line checks for self-attenuating source - should we extend to trace sources? If we
        //  do we need to modify our serialization logic a bit.
        if( model->sourceType(index) == ModelSourceType::Intrinsic )
          continue;
      }
      
      numNucAdded += 1;
      m_isoSelect->addItem( nuc->symbol );
    }//for( int index = 0; index < numNucs; ++index )
    
    if( indexToSelect > 0 )
      m_isoSelect->setCurrentIndex( indexToSelect );
    // else emit signal saying nuclide changed.
  }//updateAvailableIsotopes()
};//class TraceSrcDisplay



SourceCheckbox::SourceCheckbox( const SandiaDecay::Nuclide *nuclide,
                               double massFrac, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_useAsSourceCb( NULL ),
    m_massFraction( NULL ),
    m_nuclide( nuclide )
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSelect.css" );
  
  new WText( "&nbsp;&nbsp;&nbsp;", Wt::XHTMLText, this );

  if( nuclide )
    m_useAsSourceCb = new WCheckBox( nuclide->symbol, this );
  else
    m_useAsSourceCb = new WCheckBox( "Non Source Frac", this );

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  if(!db)
    return;

  string labeltxt = " - ";
  if( nuclide )
    labeltxt += db->element(nuclide->atomicNumber)->symbol;
  labeltxt += " Mass Frac:";

  WLabel *label = new WLabel( labeltxt, this );
  
  m_massFraction = new NativeFloatSpinBox( this );
  label->setBuddy( m_massFraction );
  m_massFraction->setAutoComplete( false );
  //m_massFraction->setDecimals( 3 );
  //m_massFraction->setSingleStep( 0.01 );
  m_massFraction->setRange( 0.0, 1.0 );
  //m_massFraction->setTextSize( 5 );
  m_massFraction->setWidth( 80 );
  m_massFraction->setMargin( 3, Wt::Left );
  m_massFraction->setSpinnerHidden( true );
  m_massFraction->setValue( massFrac );

//  m_massFraction->disable();

  if( !nuclide )
  {
    m_useAsSourceCb->setUnChecked();
    m_useAsSourceCb->hide();
    m_useAsSourceCb->disable();
  }//if( !nuclide )
}//SourceCheckbox constructor

SourceCheckbox::~SourceCheckbox()
{
}

double SourceCheckbox::massFraction() const
{
  return m_massFraction->value();
}

void SourceCheckbox::setMassFraction( double frac )
{
  m_massFraction->setValue( frac );
}



bool SourceCheckbox::useAsSource() const
{
  return m_useAsSourceCb->isChecked();
}

void SourceCheckbox::setUseAsSource( bool use )
{
  m_useAsSourceCb->setChecked( use );
}

const SandiaDecay::Nuclide *SourceCheckbox::isotope() const
{
  return m_nuclide;
}

Wt::EventSignal<> &SourceCheckbox::checked()
{
  return m_useAsSourceCb->checked();
}

Wt::EventSignal<> &SourceCheckbox::changed()
{
  return m_useAsSourceCb->changed();
}

Wt::EventSignal<> &SourceCheckbox::unChecked()
{
  return m_useAsSourceCb->unChecked();
}

Wt::Signal<float> &SourceCheckbox::massFractionChanged()
{
  return m_massFraction->valueChanged();
}


ShieldingSelect::ShieldingSelect( MaterialDB *materialDB,
                 Wt::WSuggestionPopup *materialSuggest,
                 Wt::WContainerWidget *parent )
: WContainerWidget( parent ),
  m_toggleImage( nullptr ),
  m_shieldSrcDisp( nullptr ),
  m_forFitting( false ),
  m_materialDB( materialDB ),
  m_sourceModel( nullptr ),
  m_geometry( GeometryType::Spherical ),
  m_materialSuggest( materialSuggest ),
  m_materialEdit( nullptr ),
  m_isGenericMaterial( false ),
  m_materialSummary( nullptr ),
  m_closeIcon( nullptr ),
  m_addIcon( nullptr ),
  m_addTraceSourceItem( nullptr ),
  m_dimensionsStack( nullptr ),
  m_genericDiv( nullptr ),
  m_arealDensityEdit( nullptr ),
  m_fitArealDensityCB( nullptr ),
  m_atomicNumberEdit( nullptr ),
  m_fitAtomicNumberCB( nullptr ),
  m_sphericalDiv( nullptr ),
  m_thicknessEdit( nullptr ),
  m_fitThicknessCB( nullptr ),
  m_cylindricalDiv( nullptr ),
  m_cylRadiusEdit( nullptr ),
  m_fitCylRadiusCB( nullptr ),
  m_cylLengthEdit( nullptr ),
  m_fitCylLengthCB( nullptr ),
  m_rectangularDiv( nullptr ),
  m_rectWidthEdit( nullptr ),
  m_fitRectWidthCB( nullptr ),
  m_rectHeightEdit( nullptr ),
  m_fitRectHeightCB( nullptr ),
  m_rectDepthEdit( nullptr ),
  m_fitRectDepthCB( nullptr ),
  m_fitMassFrac( nullptr ),
  m_asSourceCBs( nullptr ),
  m_traceSources( nullptr )
{
  init();
}


ShieldingSelect::ShieldingSelect( MaterialDB *materialDB,
                                  SourceFitModel *sourceModel,
                                  Wt::WSuggestionPopup *materialSuggest,
                                  ShieldingSourceDisplay *shieldSource,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_toggleImage( nullptr ),
    m_shieldSrcDisp( shieldSource ),
    m_forFitting( (shieldSource != nullptr) ),
    m_materialDB( materialDB ),
    m_sourceModel( sourceModel ),
    m_geometry( GeometryType::Spherical ),
    m_materialSuggest( materialSuggest ),
    m_materialEdit( NULL ),
    m_isGenericMaterial( false ),
    m_materialSummary( NULL ),
    m_closeIcon( NULL ),
    m_addIcon( NULL ),
    m_addTraceSourceItem( nullptr ),
    m_dimensionsStack( nullptr ),
    m_genericDiv( nullptr ),
    m_arealDensityEdit( nullptr ),
    m_fitArealDensityCB( nullptr ),
    m_atomicNumberEdit( nullptr ),
    m_fitAtomicNumberCB( nullptr ),
    m_sphericalDiv( nullptr ),
    m_cylindricalDiv( nullptr ),
    m_cylRadiusEdit( nullptr ),
    m_fitCylRadiusCB( nullptr ),
    m_cylLengthEdit( nullptr ),
    m_fitCylLengthCB( nullptr ),
    m_rectangularDiv( nullptr ),
    m_thicknessEdit( NULL ),
    m_fitThicknessCB( NULL ),
    m_fitMassFrac( NULL ),
    m_asSourceCBs( NULL ),
    m_traceSources( nullptr )
{
  init();
}


void ShieldingSelect::setClosableAndAddable( bool closeable , WGridLayout* layout )
{
  if( closeable )
  {
    if( m_closeIcon )
      return;
    
    m_closeIcon = new WPushButton(); //WText();
    m_addIcon = new WPushButton();
    m_closeIcon->setStyleClass( "ShieldingDelete Wt-icon" );
    m_closeIcon->setIcon("InterSpec_resources/images/minus_min_black.svg");
    m_addIcon->setStyleClass( "ShieldingAdd Wt-icon" );
    m_addIcon->setIcon("InterSpec_resources/images/plus_min_black.svg");
    

//    m_closeIcon->setToolTip( "Remove this shielding" );
//    m_addIcon->setToolTip( "Add a shielding" );
//    WPopupMenuItem *item = popup->addItem( "Before this shielding" );
//    item->triggered().connect( this, &ShieldingSelect::emitAddBeforeSignal )
//    item = popup->addItem( "After this shielding" );
//    item->triggered().connect( this, &ShieldingSelect::emitAddBeforeSignal )
//    m_addIcon->setMenu( popup );
    
    PopupDivMenu *popup = new PopupDivMenu( NULL, PopupDivMenu::TransientMenu );
    PopupDivMenuItem *item = popup->addMenuItem( "Add shielding before" );
    item->triggered().connect( this, &ShieldingSelect::emitAddBeforeSignal );
    item = popup->addMenuItem( "Add Shielding after" );
    item->triggered().connect( this, &ShieldingSelect::emitAddAfterSignal );
    m_addTraceSourceItem = popup->addMenuItem( "Add Trace Source" );
    m_addTraceSourceItem->triggered().connect( this, &ShieldingSelect::addTraceSource );
    
    m_addIcon->setMenu( popup );
    
    layout->addWidget( m_closeIcon, 0, 2, AlignMiddle | AlignRight );
    layout->addWidget( m_addIcon, 1, 2, AlignTop | AlignRight );
    
    m_closeIcon->clicked().connect( this, &ShieldingSelect::emitRemoveSignal );
  }else
  {
    if( m_closeIcon )
    {
      delete m_closeIcon;
      m_closeIcon = NULL;
    }
    
    if( m_addIcon )
    {
      delete m_addIcon;
      m_addIcon = 0;
    }
  }//if( closeable ) / else
}//void ShieldingSelect::setClosableAndAddable( bool closeable )


bool ShieldingSelect::fitForMassFractions() const
{
  if( !m_fitMassFrac || !m_fitMassFrac->isChecked() )
    return false;
  
  if( isGenericMaterial() )
    return false;
  
  int nchecked = 0;
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    for( WWidget *widget : etnm.second->children() )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      nchecked += (src && (src->useAsSource()));
    }//for( WWidget *widget : isotopeDiv->children() )
  }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )

  return (nchecked > 1);
}//bool fitForMassFractions() const


void ShieldingSelect::setMassFraction( const SandiaDecay::Nuclide *nuc,
                                       double fraction )
{
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    for( WWidget *widget : etnm.second->children() )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      if( src && (nuc == src->isotope()) )
      {
        src->setMassFraction( fraction );
        std::shared_ptr<Material> mat = material();
        
        for( Material::NuclideFractionPair &nfp : mat->nuclides )
        {
          if( nfp.first == nuc )
            nfp.second = fraction;
        }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )

        return;
      }//if( nuc == src )
    }//for( WWidget *widget : isotopeDiv->children() )
  }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
 
  
  throw runtime_error( "ShieldingSelect::setMassFraction(...): could not"
                       " match a source isotope with nuclide passed in." );
}//void setMassFraction( const SandiaDecay::Nuclide *nuc, double fraction );


void ShieldingSelect::setMaterialNameAndThickness( const string &name,
                                                   const string &thickness )
{
  if( m_isGenericMaterial )
    handleToggleGeneric();
  
  const Material *mat = material( name );
  
  if( !mat && name.size() )
    throw runtime_error( "'" + name + "' was not a recognized material." );
  
  if( thickness.size() )
    PhysicalUnits::stringToDistance( thickness );
  
  m_materialEdit->setText( name );
  m_thicknessEdit->setText( thickness );
  
  handleMaterialChange();
}//void setMaterialNameAndThickness(...)


void ShieldingSelect::setAtomicNumberAndArealDensity( const double an, const double ad )
{
  if( an < 1.0 || an > 100.0 )
    throw runtime_error( "setAtomicNumberAndArealDensity: Atomic number must be between 1 and 100." );
  
  const double ad_gcm2 = ad * PhysicalUnits::cm2/PhysicalUnits::g;
  if( ad_gcm2 < 0.0 || ad_gcm2 > 1000.0 )
    throw runtime_error( "setAtomicNumberAndArealDensity: Areal density must be between 0 and 1000 g/cm2." );
  
  if( !m_isGenericMaterial )
    handleToggleGeneric();
  
  m_atomicNumberEdit->setText( std::to_string(an) );
  m_arealDensityEdit->setText( std::to_string(ad_gcm2) );
  
  handleMaterialChange();
}//void setAtomicNumberAndArealDensity( const double an, const double ad )



WLineEdit *ShieldingSelect::materialEdit()
{
  return m_materialEdit;
}



void ShieldingSelect::setSphericalThickness( const double thickness )
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  if( thickness < 0.0 )
    m_thicknessEdit->setText( "" );
  else
    m_thicknessEdit->setText( PhysicalUnits::printToBestLengthUnits(thickness,3) );

  handleMaterialChange();
}//void setSphericalThickness( const double thickness )


void ShieldingSelect::setCylindricalRadiusThickness( const double radius )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  if( radius < 0.0 )
    m_cylRadiusEdit->setText( "" );
  else
    m_cylRadiusEdit->setText( PhysicalUnits::printToBestLengthUnits(radius,3) );
  
  handleMaterialChange();
}//void setCylindricalRadiusThickness( const double radius )


void ShieldingSelect::setCylindricalLengthThickness( const double length )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  if( length < 0.0 )
    m_cylLengthEdit->setText( "" );
  else
    m_cylLengthEdit->setText( PhysicalUnits::printToBestLengthUnits(length,3) );
  
  handleMaterialChange();
}//void setCylindricalLengthThickness( const double length )


void ShieldingSelect::setRectangularWidthThickness( const double width )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( width < 0.0 )
    m_rectWidthEdit->setText( "" );
  else
    m_rectWidthEdit->setText( PhysicalUnits::printToBestLengthUnits(width,3) );
  
  handleMaterialChange();
}//void setRectangularWidth( const double width )


void ShieldingSelect::setRectangularHeightThickness( const double height )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( height < 0.0 )
    m_rectHeightEdit->setText( "" );
  else
    m_rectHeightEdit->setText( PhysicalUnits::printToBestLengthUnits(height,3) );
  
  handleMaterialChange();
}//void setRectangularHeight( const double height )


void ShieldingSelect::setRectangularDepthThickness( const double depth )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( depth < 0.0 )
    m_rectDepthEdit->setText( "" );
  else
    m_rectDepthEdit->setText( PhysicalUnits::printToBestLengthUnits(depth,3) );
  
  handleMaterialChange();
}//void setRectangularDepth( const double depth )



void ShieldingSelect::setSphericalThicknessEditEnabled( const bool enabled )
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  m_thicknessEdit->setEnabled(enabled);
}


Wt::WLineEdit *ShieldingSelect::arealDensityEdit()
{
  return m_arealDensityEdit;
}


Wt::WLineEdit *ShieldingSelect::atomicNumberEdit()
{
  return m_atomicNumberEdit;
}


vector<const SandiaDecay::Nuclide *> ShieldingSelect::traceSourceNuclides() const
{
  vector<const SandiaDecay::Nuclide *> answer;
  
  if( m_traceSources )
  {
    for( WWidget *w : m_traceSources->children() )
    {
      const TraceSrcDisplay *src = dynamic_cast<const TraceSrcDisplay *>( w );
      assert( src );
      
      if( src && src->nuclide() )
        answer.push_back( src->nuclide() );
    }//for( WWidget *w : m_traceSources->children() )
  }//if( m_traceSources )
  return answer;
}//std::vector<const SandiaDecay::Nuclide *> traceSourceNuclides() const;


bool ShieldingSelect::isTraceSourceForNuclide( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  
  return (w != nullptr);
}


double ShieldingSelect::traceSourceTotalActivity( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "traceSourceTotalActivity: called with invalid nuclide" );
  
  return w->totalActivity();
}

void ShieldingSelect::setTraceSourceTotalActivity( const SandiaDecay::Nuclide *nuc,
                                                  const double activity )
{
  TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "setTraceSourceTotalActivity: called with invalid nuclide" );
  
  w->setTraceSourceTotalActivity( activity );
}//void setTraceSourceTotalActivity( nuclide, activity );


double ShieldingSelect::traceSourceDisplayActivity( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "traceSourceDisplayActivity: called with invalid nuclide" );
  
  return w->displayActivity();
}


TraceActivityType ShieldingSelect::traceSourceType( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "traceSourceType: called with invalid nuclide" );
  
  return w->activityType();
}//TraceActivityType traceSourceType( const SandiaDecay::Nuclide *nuc ) const


bool ShieldingSelect::fitTraceSourceActivity( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "traceSourceType: called with invalid nuclide" );
  
  const bool allow = w->allowFittingActivity();
  
  return allow;
}//TraceActivityType traceSourceType( const SandiaDecay::Nuclide *nuc ) const



double ShieldingSelect::updateTotalTraceSourceActivityForGeometryChange( const SandiaDecay::Nuclide *nuc )
{
  TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "updateTotalTraceSourceActivityForGeometryChange: called with invalid nuclide" );
  
  w->updateTotalActivityFromDisplayActivity();
  
  return w->totalActivity();
}//double updateTotalTraceSourceActivityForGeometryChange( const SandiaDecay::Nuclide *nuc )


const TraceSrcDisplay *ShieldingSelect::traceSourceWidgetForNuclide( const SandiaDecay::Nuclide *nuc ) const
{
  if( !m_traceSources || !nuc )
    return nullptr;
  
  for( WWidget *w : m_traceSources->children() )
  {
    const TraceSrcDisplay *src = dynamic_cast<const TraceSrcDisplay *>( w );
    assert( src );

    if( src && (src->nuclide() == nuc) )
      return src;
  }//for( WWidget *w : m_traceSources->children() )
  
  return nullptr;
}//traceSourceWidgetForNuclide(...) const


TraceSrcDisplay *ShieldingSelect::traceSourceWidgetForNuclide( const SandiaDecay::Nuclide *nuc )
{
  if( !m_traceSources || !nuc )
    return nullptr;
  
  for( WWidget *w : m_traceSources->children() )
  {
    TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
    assert( src );
    
    if( src && (src->nuclide() == nuc) )
      return src;
  }//for( WWidget *w : m_traceSources->children() )
  
  return nullptr;
}//traceSourceWidgetForNuclide(...) non-const


void ShieldingSelect::init()
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSelect.css" );
  
  //TODO/NOTE: had to hard code this as false because there is no way
  //to easily get the preference via InterSpec because
  //is still initializing when calling at this moment.
  bool showToolTips = false;
  if( auto interspec = InterSpec::instance() )
    showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", interspec );
  
  addStyleClass( "ShieldingSelect" );
  
  if( m_materialSuggest )
  {
    if( m_materialSuggest->objectName().empty() )
      m_materialSuggest->setObjectName( "ShieldingSuggest" + id() );
    m_materialSuggestName = m_materialSuggest->objectName();
  }//if( m_materialSuggest )
  
//  int voidIndex = -1;
  WContainerWidget *materialDiv = new WContainerWidget( this );
  WGridLayout* materialDivLayout = new WGridLayout();
  materialDivLayout->setContentsMargins(2,2,2,2);
  materialDiv->setLayout(materialDivLayout);
  
  
  
  m_toggleImage = new Wt::WImage(Wt::WLink("InterSpec_resources/images/shield.png"));
  m_toggleImage->clicked().connect( this,&ShieldingSelect::handleToggleGeneric );
  m_toggleImage->decorationStyle().setCursor(PointingHandCursor);
  m_toggleImage->addStyleClass( "Wt-icon" );
 
  HelpSystem::attachToolTipOn( m_toggleImage,
    "Toggle between material and generic shielding",
                              showToolTips, HelpSystem::ToolTipPosition::Top );
  
  materialDivLayout->addWidget( m_toggleImage, 0, 0, AlignLeft );
  
  m_materialEdit = new WLineEdit( "" );
  m_materialEdit->setAutoComplete( false );

  m_materialEdit->changed().connect( this, &ShieldingSelect::handleMaterialChange );
  m_materialEdit->enterPressed().connect( this, &ShieldingSelect::handleMaterialChange );
  //m_materialEdit->blurred().connect( this, &ShieldingSelect::handleMaterialChange );

  
  if( m_forFitting )
  {
    materialDivLayout->addWidget( m_materialEdit, 0, 1, AlignmentFlag::AlignMiddle );
  }else
  {
    materialDivLayout->addWidget( m_materialEdit, 0, 1 );
  }
  
  HelpSystem::attachToolTipOn( m_materialEdit,
    "You can either enter the name of a pre-defined material or element here"
    " (clear form text and click arrow on right of form to see all predefined"
    " options), or you can specify the atomic make up of the material, similar"
    " to C0.5H0.2Ni0.6, where the numbers are the density in g/cm3 of the"
    " preceding element in the material, so the example would have a total"
    " density of 0.5+0.2+0.6=1.3 g/cm3."
    " To enter materials with isotopic components, you should single or double"
    " quote the nuclide, ex: 'U238'0.2'U235'0.8",
                              showToolTips, HelpSystem::ToolTipPosition::Top );

  
  if( m_materialSuggest )
    m_materialSuggest->forEdit( m_materialEdit,
                   WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  
  m_materialSummary = new WText( "", XHTMLText );
  if( m_forFitting )
  {
    materialDivLayout->addWidget( m_materialSummary, 1, 1, AlignMiddle );
  }else
  {
    m_materialSummary->addStyleClass( "MaterialSummary" );
  }
  
  materialDivLayout->setColumnStretch(1,1);
  if( m_forFitting )
    setClosableAndAddable( true,  materialDivLayout );

  m_dimensionsStack = new WStackedWidget( this );
  
  
  // Begin setting up generic material widgets
  m_genericDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_genericDiv );
  
  WGridLayout *genericMatLayout = new WGridLayout( m_genericDiv );
  genericMatLayout->setContentsMargins(3,3,3,3);
  WLabel *label = new WLabel( "AD" );
  label->setAttributeValue( "style", "padding-left: 1em;" );
  label->setToolTip( "Areal Density of the shielding in g/cm2" );
  genericMatLayout->addWidget( label, 0, 2+m_forFitting, AlignMiddle );
  
  m_arealDensityEdit = new WLineEdit();
  m_arealDensityEdit->setToolTip( "Areal Density of the shielding in g/cm2" );
  if( m_forFitting )
    m_arealDensityEdit->setText( "15.0" );
  else
    m_arealDensityEdit->setEmptyText( "Areal Density" );
  
  genericMatLayout->addWidget( m_arealDensityEdit, 0, 3+m_forFitting, AlignMiddle );
  genericMatLayout->setColumnStretch( 3+m_forFitting, 1 );
  
  m_arealDensityEdit->setTextSize( 5 );
  label->setBuddy( m_arealDensityEdit );
  
  WDoubleValidator *adValidator = new WDoubleValidator( this );
  adValidator->setRange( 0.0, 500.0 );
  adValidator->setInvalidTooSmallText( "0.0" );
  adValidator->setInvalidTooLargeText( "500.0" );
  adValidator->setInvalidNotANumberText( "0.0" );
  adValidator->setInvalidBlankText( "0.0" );
  
  m_arealDensityEdit->setValidator( adValidator );
  label = new WLabel( "g/cm<sup>2</sup>");
  label->setAttributeValue( "style", "font-size: 75%;" );
  genericMatLayout->addWidget( label, 0, 4+m_forFitting, AlignMiddle );
  
  m_arealDensityEdit->addStyleClass( "numberValidator" ); //used to detect mobile keyboard
  
  if( m_forFitting )
  {
    m_fitArealDensityCB = new WCheckBox( "Fit" );
    m_fitArealDensityCB->setChecked( true );
    genericMatLayout->addWidget( m_fitArealDensityCB, 0, 6, AlignMiddle );
  }
  
  label = new WLabel( "AN" );
  label->setToolTip( "Atomic Number of the shielding" );
  genericMatLayout->addWidget( label, 0, 0, AlignMiddle );
  
  m_atomicNumberEdit = new WLineEdit();
  m_atomicNumberEdit->setToolTip( "Atomic Number of the shielding" );
  if( m_forFitting )
    m_atomicNumberEdit->setText( "15.0" );
  else
    m_atomicNumberEdit->setEmptyText( "Atomic Number" );
  
  genericMatLayout->addWidget( m_atomicNumberEdit, 0, 1, AlignMiddle );
  genericMatLayout->setColumnStretch( 1, 1 );
  
  m_atomicNumberEdit->setTextSize( 5 );
  label->setBuddy( m_atomicNumberEdit );
  WDoubleValidator *dblValidator = new WDoubleValidator( MassAttenuation::sm_min_xs_atomic_number, MassAttenuation::sm_max_xs_atomic_number, m_atomicNumberEdit );
  m_atomicNumberEdit->setValidator( dblValidator );
  m_atomicNumberEdit->addStyleClass( "numberValidator"); //used to detect mobile keyboard
  dblValidator->setInvalidTooSmallText( std::to_string(MassAttenuation::sm_min_xs_atomic_number) );
  dblValidator->setInvalidTooLargeText( std::to_string(MassAttenuation::sm_max_xs_atomic_number) );
  dblValidator->setInvalidNotANumberText( "1.0" );
  dblValidator->setInvalidBlankText( "1.0" );
  
  if( m_forFitting )
  {
    m_fitAtomicNumberCB = new WCheckBox( "Fit" );
    m_fitAtomicNumberCB->setChecked( false );
    genericMatLayout->addWidget( m_fitAtomicNumberCB, 0, 2, AlignMiddle );
    
    m_asSourceCBs = new WContainerWidget( this );
    m_asSourceCBs->addStyleClass( "ShieldingAsSourceCBDiv" );
    label = new WLabel( "Source for:", m_asSourceCBs );
    label->setInline( false );
    m_asSourceCBs->hide();
    const char *tooltip = "When these nuclides are used as sources they are"
    " treated as uniformly distributed in the material"
    " (which is assumed spherical), so self attenuation"
    " and other factors are accounted for.";
    
    HelpSystem::attachToolTipOn( m_asSourceCBs,tooltip, showToolTips );
    m_fitMassFrac = new WCheckBox( "Fit Mass Fractions", m_asSourceCBs );
    m_fitMassFrac->hide();
    m_fitMassFrac->setInline( false );
  }//if( m_forFitting )
  
  
  // Note: the changed() signal should catch when you change the value and then the input is blurred
  m_atomicNumberEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_atomicNumberEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  
  m_arealDensityEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_arealDensityEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  
  
  // A validator for the various geometry distances
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  
  
  // A lamda for setting up dimension edits for the various geometry dimensions
  auto setupDimEdit = [this,distValidator]( const char *labelTxt, WLineEdit *&edit, WCheckBox *&fitCb, WGridLayout *grid ){
    const int row = grid->rowCount();
    
    WLabel *label = new WLabel( labelTxt );
    grid->addWidget( label, row, 0, AlignMiddle );
    
    edit = new WLineEdit( "1.0 cm" );
    edit->setValidator( distValidator );
    edit->setAutoComplete( false );
    label->setBuddy( edit );
    
    //From a very brief experiment, it looks like the below JS would remove the uncertainty text,
    //  however it appears when the server pushes the value of the edit to the client, the 'value'
    //  tag of the element isnt updated... I have no clue.
    //  const char *focusjs = "function(s,e){try{"
    //  "s.value = s.value.replace('\\(.*\\)','').value.replace('  ',' ');"
    //  "}catch(e){}}";
    //  edit->changed().connect( focusjs );
    
    edit->changed().connect( boost::bind( &ShieldingSelect::removeUncertFromDistanceEdit, this, edit) );
    edit->enterPressed().connect( boost::bind( &ShieldingSelect::removeUncertFromDistanceEdit, this, edit) );
    //edit->blurred().connect( boost::bind( &ShieldingSelect::removeUncertFromDistanceEdit, this, edit) );
    
    edit->changed().connect( this, &ShieldingSelect::handleMaterialChange );
    edit->enterPressed().connect( this, &ShieldingSelect::handleMaterialChange );
    //edit->blurred().connect( this, &ShieldingSelect::handleMaterialChange );
    
    if( m_forFitting )
    {
      edit->setWidth( 150 );
      grid->addWidget( edit, row, 1, AlignMiddle );
      
      fitCb = new WCheckBox( "Fit" );
      fitCb->setChecked( false );
      grid->addWidget( fitCb, row, 2, AlignMiddle | AlignRight );
      grid->setColumnStretch( 3, 1 );
    }else
    {
      grid->addWidget( edit, row, 1 );
    }//if( m_forFitting ) / else
  };//setupDimEdit(...)
  
  
  
  // Begin setting up spherical widgets
  m_sphericalDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_sphericalDiv );
  
  WGridLayout *sphericalLayout = new WGridLayout();
  m_sphericalDiv->setLayout( sphericalLayout );
  sphericalLayout->setContentsMargins( 3, 3, 3, 3 );
  
  setupDimEdit( "Thickness", m_thicknessEdit, m_fitThicknessCB, sphericalLayout );
  
  if( !m_forFitting )
  {
    // TODO: right now the only place a ShieldingSelect will be non-generic or non-spherical is
    //       in the Shielding/Activity fit tool, which we always have (m_forFitting == true), but
    //       if geometry is cylindrical or rectangular, and (m_forFitting == false), then the
    //       material summary wont be visible.
    sphericalLayout->addWidget( m_materialSummary, 0, 2, AlignMiddle );
    sphericalLayout->setColumnStretch( 1, 1 );
  }//if( !m_forFitting )
  
  // Begin setting up cylindrical widgets
  m_cylindricalDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_cylindricalDiv );
  
  WGridLayout *cylindricalLayout = new WGridLayout();
  m_cylindricalDiv->setLayout( cylindricalLayout );
  cylindricalLayout->setContentsMargins( 3, 3, 3, 3 );
  
  setupDimEdit( "Radius", m_cylRadiusEdit, m_fitCylRadiusCB, cylindricalLayout );
  setupDimEdit( "Length", m_cylLengthEdit, m_fitCylLengthCB, cylindricalLayout );

  
  // Begin setting up rectangular widgets
  m_rectangularDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_rectangularDiv );
  
  WGridLayout *rectangularLayout = new WGridLayout();
  m_rectangularDiv->setLayout( rectangularLayout );
  rectangularLayout->setContentsMargins( 3, 3, 3, 3 );
  
  setupDimEdit( "Width", m_rectWidthEdit, m_fitRectWidthCB, rectangularLayout );
  setupDimEdit( "Height", m_rectHeightEdit, m_fitRectHeightCB, rectangularLayout );
  setupDimEdit( "Depth", m_rectDepthEdit, m_fitRectDepthCB, rectangularLayout );
  
  // We're all done creating the widgets
  
  handleMaterialChange();
  
  if( m_forFitting )
  {
    InterSpec *interspec = InterSpec::instance();
    const bool isMobile = (interspec && interspec->isMobile());
    if( !isMobile )
      m_materialEdit->setFocus();
  }
}//void ShieldingSelect::init()


ShieldingSelect::~ShieldingSelect()
{
  if( m_materialSuggest && m_materialEdit )
  {
    WApplication *app = wApp;
    
    WContainerWidget *root = (app ? app->domRoot() : (WContainerWidget *)0);
    
    if( root )
    {
      //Testing for root apears to be enough, and m_materialSuggestName does not
      //  need to be done for our current specific use case, but leaving in the
      //  checking of the name for the futire, JIC.
      WWidget *w = app->findWidget( m_materialSuggestName );
      if( w || m_materialSuggestName.empty() )
        m_materialSuggest->removeEdit( m_materialEdit );
      else
       cerr << "~ShieldingSelect(): Suggest not in DOM, not removing form from suggestion" << endl;
    }else
    {
      cerr << "~ShieldingSelect(): no DOM root, not removing form from suggestion" << endl;
    }
  }//if( m_materialSuggest && m_materialEdit )
  
}//~ShieldingSelect()


void ShieldingSelect::emitRemoveSignal()
{
  if( m_asSourceCBs )
  {
    for( const ElementToNuclideMap::value_type &etnp : m_sourceIsotopes )
    {
//      const SandiaDecay::Element *el = etnp.first;
      WContainerWidget *cont = etnp.second;

      const vector<WWidget *> &children = cont->children();
      for( WWidget *child : children )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( !cb )
          continue;
        if( cb->useAsSource() && cb->isotope() )
        {
          cb->setUseAsSource( false );
          removingIsotopeAsSource().emit( cb->isotope(), ModelSourceType::Intrinsic );
        }//if( cb->useAsSource() )
      }//for( WWidget *child : children )
    }//for(...)
  }//if( m_asSourceCBs )

  m_removeSignal.emit( this );
}//void emitRemoveSignal()


void ShieldingSelect::emitAddBeforeSignal()
{
  m_addShieldingBefore.emit( this );
}


void ShieldingSelect::emitAddAfterSignal()
{
  m_addShieldingAfter.emit( this );
}



void ShieldingSelect::checkIsCorrectCurrentGeometry( const GeometryType wanted, const char *fcn ) const
{
  bool okayGeom = false;
  
  const bool is_generic = isGenericMaterial();
  
  switch( wanted )
  {
    case GeometryType::Spherical:
      okayGeom = (!is_generic && (m_geometry == GeometryType::Spherical));
      break;
      
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
      okayGeom = (!is_generic
                  && ((m_geometry == GeometryType::CylinderEndOn)
                      || (m_geometry == GeometryType::CylinderSideOn)));
      break;
      
    case GeometryType::Rectangular:
      okayGeom = (!is_generic && (m_geometry == GeometryType::Rectangular));
      break;
      
    case GeometryType::NumGeometryType:
      okayGeom = is_generic;
      break;
  }//switch( wanted )
  
  
  if( okayGeom )
    return;
  
  const string current_geom = GammaInteractionCalc::to_str(m_geometry);
  const string wanted_geom = GammaInteractionCalc::to_str(wanted);
  //if( (wanted == GeometryType::CylinderEndOn) || (wanted == GeometryType::CylinderSideOn) )
  //  remove "EndOn" or "SideOn"
  
  if( wanted == GeometryType::NumGeometryType )
    throw runtime_error( "ShieldingSelect::" + string(fcn) + ": can not be called when"
                        " a not a generic material (currently " + current_geom + ")" );
  
  throw runtime_error( "ShieldingSelect::" + string(fcn) + ": can not be called when"
                        " a generic material or non-" + wanted_geom
                        + " (currently " + current_geom + ")" );
}//void checkIsCorrectCurrentGeometry(...)




void ShieldingSelect::addTraceSource()
{
  if( !m_traceSources )
  {
    m_traceSources = new WContainerWidget( this );
    m_traceSources->addStyleClass( "TraceSrcContainer" );
  }
  
  TraceSrcDisplay *src = new TraceSrcDisplay( this );
  src->nucChangedSignal().connect( this, &ShieldingSelect::handleTraceSourceNuclideChange );
  src->activityUpdated().connect( this, &ShieldingSelect::handleTraceSourceActivityChange );
}//void addTraceSource()


void ShieldingSelect::removeTraceSourceWidget( TraceSrcDisplay *toRemove )
{
  if( !toRemove )
    return;
  
  for( WWidget *w : m_traceSources->children() )
  {
    TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
    assert( src );
    if( !src ) // JIC
      continue;
    
    if( src == toRemove )
    {
      handleTraceSourceWidgetAboutToBeRemoved( src );
      
      delete src;
      
      setTraceSourceMenuItemStatus();
      
      return;
    }//if( src == toRemove )
  }//for( WWidget *w : m_traceSources->children() )
  
  
  throw runtime_error( "ShieldingSelect::removeTraceSourceWidget(TraceSrcDisplay *): called with invalid display." );
}//void removeTraceSourceWidget( TraceSrcDisplay *src )



void ShieldingSelect::setTraceSourceMenuItemStatus()
{
  if( !m_addTraceSourceItem )
    return;
  
  if( m_isGenericMaterial || !m_sourceModel )
  {
    m_addTraceSourceItem->setDisabled( true );
    return;
  }//if( m_isGenericMaterial )
  
  int numAvailableNuclides = 0;
  const int numNuclides = m_sourceModel->numNuclides();
  for( int nucIndex = 0; nucIndex < numNuclides; ++nucIndex )
  {
    const bool selfAttenSrc = m_sourceModel->isVolumetricSource( nucIndex );
    if( !selfAttenSrc )
      numAvailableNuclides += 1;
  }
  
  m_addTraceSourceItem->setDisabled( !numAvailableNuclides );
}//void setTraceSourceMenuItemStatus()


void ShieldingSelect::handleTraceSourceNuclideChange( TraceSrcDisplay *changedSrc, const SandiaDecay::Nuclide *oldNuc )
{
  assert( changedSrc );
  const SandiaDecay::Nuclide *nuc = changedSrc->nuclide();
  
  if( oldNuc )
    removingIsotopeAsSource().emit( oldNuc, ModelSourceType::Trace );
  
  if( changedSrc && changedSrc->nuclide() )
    addingIsotopeAsSource().emit( changedSrc->nuclide(), ModelSourceType::Trace );
  
  vector<pair<const SandiaDecay::Nuclide *,double>> traceNucs;
  
  const vector<WWidget *> &traceSources = m_traceSources->children();
  for( WWidget *w : traceSources )
  {
    TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
    assert( src );
    if( !src ) // JIC
      continue;
      
    if( src == changedSrc )
      continue;
    
    src->unselectNuclide( nuc );
    src->updateAvailableIsotopes();
    
    const SandiaDecay::Nuclide *srcNuc = src->nuclide();
    if( srcNuc )
    {
      traceNucs.push_back( { srcNuc, src->totalActivity() } );
    }
  }//for( WWidget *w : traceSources )
  
  setTraceSourceMenuItemStatus();
}//void handleTraceSourceNuclideChange( TraceSrcDisplay *src );


void ShieldingSelect::handleTraceSourceActivityChange( const SandiaDecay::Nuclide *nuc, const double activity )
{
  m_activityFromVolumeNeedUpdating.emit( this, nuc );
}

void ShieldingSelect::handleTraceSourceWidgetAboutToBeRemoved( TraceSrcDisplay *src )
{
  if( !src || !src->nuclide() )
    return;
    
  removingIsotopeAsSource().emit( src->nuclide(), ModelSourceType::Trace );
}//handleTraceSourceWidgetAboutToBeRemoved(...)


double ShieldingSelect::shieldingVolume() const
{
  if( isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::shieldingVolume(): You should not call this function"
                         " for generic materials" );
  
  double volume = -1.0;
  
  const double pi = PhysicalUnits::pi;
  
  switch( m_geometry )
  {
    case GeometryType::Spherical:
    {
      // Spheres are defined by their thickness relative to their inner sphere, so we have to loop
      //  a bit to get our current radius.
      double inner_rad = 0.0;
      if( m_shieldSrcDisp )
      {
        for( const ShieldingSelect *inner = m_shieldSrcDisp->innerShielding(this);
            inner; inner = m_shieldSrcDisp->innerShielding(inner) )
        {
          assert( inner->geometry() == m_geometry );
          inner_rad += inner->thickness();
        }
      }//if( m_shieldSrcDisp )
      
      const double outer_rad = inner_rad + thickness();
      
      const double innerVolume = (4.0/3.0) * pi * inner_rad * inner_rad * inner_rad;
      const double outerVolume = (4.0/3.0) * pi * outer_rad * outer_rad * outer_rad;
      assert( outerVolume > innerVolume );
      
      volume = outerVolume - innerVolume;
      
      break;
    }//case GeometryType::Spherical:
      
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
    {
      double inner_rad = 0.0, inner_half_length = 0.0;
      if( m_shieldSrcDisp )
      {
        for( const ShieldingSelect *inner = m_shieldSrcDisp->innerShielding(this);
            inner; inner = m_shieldSrcDisp->innerShielding(inner) )
        {
          assert( inner->geometry() == m_geometry );
          inner_rad += inner->cylindricalRadiusThickness();
          inner_half_length += inner->cylindricalLengthThickness();
        }
      }//if( m_shieldSrcDisp )
      
      const double innerVolume = pi * inner_rad * inner_rad * 2.0 * inner_half_length;
      
      
      const double rad = inner_rad + cylindricalRadiusThickness();
      const double len = inner_half_length + cylindricalLengthThickness();
      const double outerVolume = pi * rad * rad * 2.0 * len;
      
      assert( outerVolume >= innerVolume );
      
      volume = outerVolume - innerVolume;
      
      break;
    }//case CylinderEndOn and CylinderSideOn
      
    case GeometryType::Rectangular:
    {
      double inner_half_width = 0.0, inner_half_height = 0.0, inner_half_depth = 0.0;
      if( m_shieldSrcDisp )
      {
        for( const ShieldingSelect *inner = m_shieldSrcDisp->innerShielding(this);
            inner; inner = m_shieldSrcDisp->innerShielding(inner) )
        {
          assert( inner->geometry() == m_geometry );
          inner_half_width  += inner->rectangularWidthThickness();
          inner_half_height += inner->rectangularHeightThickness();
          inner_half_depth  += inner->rectangularDepthThickness();
        }
      }//if( m_shieldSrcDisp )
      
      
      const double innerVolume = 8.0 * inner_half_width * inner_half_height * inner_half_depth;
      
      const double half_width  = inner_half_width  + rectangularWidthThickness();
      const double half_height = inner_half_height + rectangularHeightThickness();
      const double half_depth  = inner_half_depth  + rectangularDepthThickness();
      
      const double outerVolume = 8.0 * half_width * half_height * half_depth;
      
      assert( outerVolume >= innerVolume );
      
      volume = outerVolume - innerVolume;
      
      break;
    }//case Rectangular:
      
    case GeometryType::NumGeometryType:
      assert(0);
      throw runtime_error("shieldingVolume(): invalid geometry");
      break;
  }//switch( m_geometry )

  assert( volume >= 0.0 );
  
  return volume;
}//double ShieldingSelect::shieldingVolume()


double ShieldingSelect::shieldingMass() const
{
  if( isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::shieldingMass(): You should not call this function"
                         " for generic materials" );

  if( !m_currentMaterial )
    return 0.0;
  
  const double volume = shieldingVolume();
  
  return m_currentMaterial->density * volume;
}//double ShieldingSelect::shieldingMass()







Wt::Signal<ShieldingSelect *> &ShieldingSelect::remove()
{
  return m_removeSignal;
}

Wt::Signal<ShieldingSelect *> &ShieldingSelect::addShieldingBefore()
{
  return m_addShieldingBefore;
}


Wt::Signal<ShieldingSelect *> &ShieldingSelect::addShieldingAfter()
{
  return m_addShieldingAfter;
}


Wt::Signal<ShieldingSelect *> &ShieldingSelect::materialModified()
{
  return m_materialModifiedSignal;
}


Wt::Signal<ShieldingSelect *> &ShieldingSelect::materialChanged()
{
  return m_materialChangedSignal;
}


Wt::Signal<ShieldingSelect *,const SandiaDecay::Nuclide *> &ShieldingSelect::activityFromVolumeNeedUpdating()
{
  return m_activityFromVolumeNeedUpdating;
}


Wt::Signal<const SandiaDecay::Nuclide *,ModelSourceType> &ShieldingSelect::addingIsotopeAsSource()
{
  return m_addingIsotopeAsSource;
}


Wt::Signal<const SandiaDecay::Nuclide *,ModelSourceType> &ShieldingSelect::removingIsotopeAsSource()
{
  return m_removingIsotopeAsSource;
}

void ShieldingSelect::setGeometry( GammaInteractionCalc::GeometryType type )
{
  assert( type != GeometryType::NumGeometryType );
  if( type == GeometryType::NumGeometryType )
    throw runtime_error( "setGeometry: invalid geometry" );
  
  m_geometry = type;
  
  displayInputsForCurrentGeometry();
}//void setGeometry( GammaInteractionCalc::GeometryType type )


GammaInteractionCalc::GeometryType ShieldingSelect::geometry() const
{
  return m_geometry;
}

bool ShieldingSelect::isGenericMaterial() const
{
  return m_isGenericMaterial;
}//bool isGenericMaterial() const




double ShieldingSelect::thickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  return distance_of_input_text( m_thicknessEdit );
}//double thickness() const


double ShieldingSelect::cylindricalRadiusThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderSideOn, __func__ );
  
  return distance_of_input_text( m_cylRadiusEdit );
}//double cylindricalRadiusThickness() const


double ShieldingSelect::cylindricalLengthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderSideOn, __func__ );
  
  return distance_of_input_text( m_cylLengthEdit );
}//double cylindricalLength() const


double ShieldingSelect::rectangularWidthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderSideOn, __func__ );
  
  return distance_of_input_text( m_rectWidthEdit );
}//double rectangularWidth() const


double ShieldingSelect::rectangularHeightThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  return distance_of_input_text( m_rectHeightEdit );
}//double rectangularHeight() const


double ShieldingSelect::rectangularDepthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  return distance_of_input_text( m_rectDepthEdit );
}//double rectangularDepth() const


double ShieldingSelect::atomicNumber() const
{
  if( !isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::atomicNumber() can not be called when a material" );
  
  double answer = 0;
  const string text = m_atomicNumberEdit->text().toUTF8();
  
  if( text.empty() )
    return 13.0;
  
  if( !(stringstream(text) >> answer) )
    throw std::runtime_error( "Error converting '" + text + "' to a atomic number");
  return answer;
}//double atomicNumber() const


double ShieldingSelect::arealDensity() const
{
  if( !isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::arealDensity() can not be called when a material" );
  
  const WString &text = m_arealDensityEdit->text();
  string txtstr = text.toUTF8();
  SpecUtils::trim( txtstr );
  
  if( txtstr.empty() )
    return 0.0;

  double answer = 0;
  if( !(stringstream(txtstr) >> answer) )
    throw std::runtime_error( "Error converting '" + txtstr + "' to an areal density");

  answer *= (PhysicalUnits::gram / PhysicalUnits::cm2);

  return answer;
}//double arealDensity() const;


bool ShieldingSelect::fitThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( string(__func__) + ": can't be called when not for fitting." );
  
  return m_fitThicknessCB->isChecked();
}//bool fitThickness() const


bool ShieldingSelect::fitCylindricalRadiusThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitCylRadiusCB->isChecked();
}//bool fitCylindricalRadiusThickness() const


bool ShieldingSelect::fitCylindricalLengthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitCylLengthCB->isChecked();
}//bool fitCylindricalLengthThickness() const


bool ShieldingSelect::fitRectangularWidthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitRectWidthCB->isChecked();
}//bool fitRectangularWidthThickness() const


bool ShieldingSelect::fitRectangularHeightThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitRectHeightCB->isChecked();
}//bool fitRectangularHeightThickness() const


bool ShieldingSelect::fitRectangularDepthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitRectDepthCB->isChecked();
}//bool fitRectangularDepthThickness() const


bool ShieldingSelect::fitAtomicNumber() const
{
  if( !isGenericMaterial() || !m_forFitting )
    throw std::runtime_error( "ShieldingSelect::fitAtomicNumber() can only be "
                              "called for a generic material" );
  return m_fitAtomicNumberCB->isChecked();
}//bool fitAtomicNumber() const


bool ShieldingSelect::fitArealDensity() const
{
  if( !isGenericMaterial() || !m_forFitting )
    throw std::runtime_error( "ShieldingSelect::fitArealDensity() can only be "
                              "called for a generic material" );
  return m_fitArealDensityCB->isChecked();
}//bool fitArealDensity() const


const Material *ShieldingSelect::material( const std::string &text )
{
  //See if 'text' is the name of a material in the database
  try
  {
    const Material *answer = m_materialDB->material( text );
    
    cerr << "ShieldingSelect::material(...)\n\tPotential err here, should account for"
         << " possibly variaed mass fractions!" << endl;
    
    return answer;
  }catch(...)
  {
    //material wasnt in the database
  }
  
  //See if 'text' is a chemical formula, if so add it to possible suggestions
  try
  {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    const Material *mat = m_materialDB->parseChemicalFormula( text, db );
    m_materialSuggest->addSuggestion( mat->name, mat->name );
    return mat;
  }catch(...)
  {
    // No luck
  }

  return nullptr;
}//std::shared_ptr<Material> material( const std::string &text )


std::shared_ptr<Material> ShieldingSelect::material()
{
  if( m_isGenericMaterial )
    return nullptr;

  const string text = SpecUtils::trim_copy( m_materialEdit->text().toUTF8() );
  if( m_currentMaterial && (text == m_currentMaterialDescrip) )
    return m_currentMaterial;

  if( text.empty() )
  {
    m_currentMaterial.reset();
    m_currentMaterialDescrip = "";
    return m_currentMaterial;
  }
  
  const Material *mat = material( text );
  
  if( mat )
  {
    m_currentMaterialDescrip = text;
    m_currentMaterial = std::make_shared<Material>( *mat );
  }else
  {
    m_currentMaterial.reset();
    m_currentMaterialDescrip = "";
  }

  return m_currentMaterial;
}//const Material *material() const;


void ShieldingSelect::updateIfMassFractionCanFit()
{
  if( !m_fitMassFrac )
    return;
  
  int nchecked = 0;
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    for( WWidget *widget : etnm.second->children() )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      nchecked += (src && (src->useAsSource()));
    }//for( WWidget *widget : isotopeDiv->children() )
  }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  
  const bool shouldHide = (nchecked < 2);
  
  if( shouldHide == m_fitMassFrac->isHidden() )
    return;
  
  if( shouldHide )
  {
    m_fitMassFrac->setUnChecked();
    m_fitMassFrac->hide();
  }else
  {
    m_fitMassFrac->show();
  }
}//void updateIfMassFractionCanFit()


void ShieldingSelect::isotopeCheckedCallback( const SandiaDecay::Nuclide *nuc )
{
  updateIfMassFractionCanFit();
  m_addingIsotopeAsSource.emit( nuc, ModelSourceType::Intrinsic );
}//void isotopeCheckedCallback( const std::string symbol )



void ShieldingSelect::isotopeUnCheckedCallback( const SandiaDecay::Nuclide *iso )
{
  updateIfMassFractionCanFit();
  m_removingIsotopeAsSource.emit( iso, ModelSourceType::Intrinsic );
}//void isotopeUnCheckedCallback( const std::string symbol )



void ShieldingSelect::uncheckSourceIsotopeCheckBox( const SandiaDecay::Nuclide *iso )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  
  if( m_traceSources )
  {
    for( WWidget *w : m_traceSources->children() )
    {
      TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
      assert( src );
      if( src ) // JIC
        src->unselectNuclide( iso );
    }//for( WWidget *w : traceSources )
  }//if( m_traceSources )
  
  if( !m_asSourceCBs || !iso || !db )
    return;

  const SandiaDecay::Element *el = db->element( iso->atomicNumber );

  if( m_sourceIsotopes.find(el) == m_sourceIsotopes.end() )
    return;

  WContainerWidget *cont = m_sourceIsotopes[el];

  const vector<WWidget *> &children = cont->children();
  for( WWidget *child : children )
  {
    SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );

    if( !cb || (cb->isotope() != iso) )
      continue;

    if( cb->useAsSource() )
    {
//      removingIsotopeAsSource().emit( symbol, ModelSourceType::Intrinsic );
      cb->setUseAsSource( false );
    }//if( cb->isChecked() )
  }//for( WWidget *child : children )
  
  updateIfMassFractionCanFit();
}//void uncheckSourceIsotopeCheckBox( const std::string &symol )


void ShieldingSelect::sourceRemovedFromModel( const SandiaDecay::Nuclide *nuc )
{
  
  // Make sure no trace sources are using nuc, and if they are, remove that trace source
  if( m_traceSources )
  {
    vector<TraceSrcDisplay *> todel;
    for( WWidget *w : m_traceSources->children() )
    {
      TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
      assert( src );
      if( !src ) // JIC
        continue;
      
      const SandiaDecay::Nuclide *srcNuc = src->nuclide();
      if( srcNuc && (nuc == srcNuc) )
        todel.push_back( src );
      else
        src->modelSourceRemoved( nuc );
    }//for( WWidget *w : traceSources )
    
    for( TraceSrcDisplay *w : todel )
    {
      handleTraceSourceWidgetAboutToBeRemoved( w );
      delete w;
    }
  }//if( m_traceSources )
  
  // Now deal with self-attenuating sources.
  if( !m_asSourceCBs )
    return;

  set<const SandiaDecay::Element *> elsToRemove;

  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const vector<WWidget *> &children = elDiv.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb/* || !cb->useAsSource() */)
        continue;

      const SandiaDecay::Nuclide *iso = cb->isotope();
      if( !iso )
        continue;

      if( iso == nuc )
      {
        elsToRemove.insert( elDiv.first );
        delete cb;
      }//if( iso == nuc )
    }//for( WWidget *child : children )
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )


  for( const SandiaDecay::Element *el : elsToRemove )
  {
    WContainerWidget *w = m_sourceIsotopes[el];
    if( w->children().empty() )
    {
      delete w;
      m_sourceIsotopes.erase( el );
    }//if( w->children().empty() )
  }//for( const SandiaDecay::Element *el : elsToRemove )


  if( m_sourceIsotopes.empty() )
    m_asSourceCBs->hide();
  
  updateIfMassFractionCanFit();
  
  //call updateMassFractionDisplays() to update the "Assuming XX% other U isos"
  updateMassFractionDisplays( m_currentMaterial );
}//void sourceRemovedFromModel( const std::string &symbol )


vector< ShieldingSelect::NucMasFrac > ShieldingSelect::sourceNuclideMassFractions()
{
  vector< ShieldingSelect::NucMasFrac > answer;
  
  const vector<const SandiaDecay::Nuclide *> nucs = selfAttenNuclides();
  std::shared_ptr<const Material> mat = material();
  
  if( !mat || nucs.empty() )
    return answer;
  
  for( const Material::NuclideFractionPair &nfp : mat->nuclides )
  {
    if( std::find( nucs.begin(), nucs.end(), nfp.first ) != nucs.end() )
      answer.push_back( nfp );
  }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )
  
  return answer;
}//std::vector< NucMasFrac > sourceNuclideMassFractions()


vector<const SandiaDecay::Nuclide *> ShieldingSelect::selfAttenNuclides()
{
  set<const SandiaDecay::Nuclide *> answer;

  std::shared_ptr<const Material> mat = material();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  if( !mat || !db )
    return vector<const SandiaDecay::Nuclide *>();

  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const vector<WWidget *> &children = elDiv.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb || !cb->useAsSource() )
        continue;

      const SandiaDecay::Nuclide *iso = cb->isotope();
      if( !iso )
        continue;

      for( const Material::NuclideFractionPair &nfp : mat->nuclides )
      {
        if( nfp.first && (nfp.first==iso) )
          answer.insert( nfp.first );
      }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )

      const SandiaDecay::Element *el = db->element( iso->atomicNumber );

      for( const Material::ElementFractionPair &nfp : mat->elements )
      {
        if( nfp.first && (nfp.first==el) )
          answer.insert( iso );
      }//for( const Material::NuclideFractionPair &nfp : mat->nuclides )
    }//for( WWidget *child : children )
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )

  vector<const SandiaDecay::Nuclide *> nucs;
  for( const SandiaDecay::Nuclide *n : answer )
    nucs.push_back( n );

  return nucs;
}//std::vector<const SandiaDecay::Nuclide *> selfAttenNuclides() const


double ShieldingSelect::massFractionOfElement( const SandiaDecay::Nuclide *iso,
                                               std::shared_ptr<const Material> mat )
{
  if( !mat )
    return 0.0;

  //Make sure the material has the isotope reequested to add
  const vector< Material::NuclideFractionPair >  &nuclides = mat->nuclides;
  const vector< Material::ElementFractionPair > &elements = mat->elements;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  if( !db || !iso )
    return 0.0;

  const SandiaDecay::Element *element = db->element( iso->atomicNumber );

  double massFracFromNuclide = 0.0, massFracOfElement = 0.0, sumMassFrac = 0.0;
  bool hasNuclide = false, hasElement = false;
  for( const Material::ElementFractionPair &efp : elements )
  {
    if( efp.first == element )
    {
      hasElement = true;
      massFracOfElement += efp.second;
    }//if( efp.first == element )
  }//for( const Material::ElementFractionPair &efp : elements )

  for( const Material::NuclideFractionPair &efp : nuclides )
  {
    if( efp.first->atomicNumber == iso->atomicNumber )
      sumMassFrac += efp.second;

    if( efp.first == iso )
    {
      hasNuclide = true;
      massFracFromNuclide += efp.second;
    }//if( efp.first == iso )
  }//for( const Material::NuclideFractionPair &efp : nuclides )

  if( !hasNuclide && !hasElement )
    throw runtime_error( "Material Doesnt Contain Isotope" );

  if( hasElement && !hasNuclide )
  {
    sumMassFrac += 1.0;
    const vector<SandiaDecay::NuclideAbundancePair> &isos = element->isotopes;

    bool hasNaturalAbundance = false;
    for( const SandiaDecay::NuclideAbundancePair &i : isos )
    {
      hasNaturalAbundance |= (i.abundance!=0.0);
      if( i.nuclide == iso )
        massFracFromNuclide += i.abundance;
    }

    if( !hasNaturalAbundance )
    {
      for( const SandiaDecay::NuclideAbundancePair &i : isos )
      {
        if( i.nuclide == iso )
          massFracFromNuclide += 1.0/isos.size();
      }
    }//if( !hasNaturalAbundance )


  }//if( hasElement )


  if( sumMassFrac == 0.0 )
    return 0.0;

  return massFracFromNuclide/sumMassFrac;
}//double massFractionOfElement( const SandiaDecay::Nuclide *iso )



void ShieldingSelect::modelNuclideAdded( const SandiaDecay::Nuclide *iso )
{
  // Deal with trace sources
  // Make sure no trace sources are using nuc, and if they are, remove that trace source
  if( m_traceSources )
  {
    for( WWidget *w : m_traceSources->children() )
    {
      TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
      assert( src );
      if( src ) // JIC
        src->modelSourceAdded( iso );
    }//for( WWidget *w : traceSources )
  }//if( m_traceSources )
  
  
  if( !m_asSourceCBs )
    return;

  std::shared_ptr<const Material> mat = material();

  if( !mat )
    return;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  if( !db || !iso )
    return;

  const SandiaDecay::Element *element = db->element( iso->atomicNumber );

  if( m_sourceIsotopes.find(element) != m_sourceIsotopes.end() )
  {
    const vector<WWidget *> children = m_sourceIsotopes[element]->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( cb && (cb->isotope()==iso) )
        return;
    }//for( WWidget *child : children )
  }//if( m_sourceIsotopes.find(element) != m_sourceIsotopes.end() )

  double massFrac = 0.0;
  try
  {
    massFrac = massFractionOfElement( iso, mat );
  }catch(...)
  {
    return;
  }

  if( m_sourceIsotopes.find(element) == m_sourceIsotopes.end() )
    m_sourceIsotopes[element] = new WContainerWidget( this );
  WContainerWidget *isotopeDiv = m_sourceIsotopes[element];
  for( WWidget *widget : isotopeDiv->children() )
  {
    SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
    if( src && (src->isotope()==iso) )
      return;
  }//for( WWidget *widget : isotopeDiv->children() )


  SourceCheckbox *cb = new SourceCheckbox( iso, massFrac );
  isotopeDiv->insertWidget( isotopeDiv->count(), cb );
  cb->checked().connect( boost::bind( &ShieldingSelect::isotopeCheckedCallback, this, iso ) );
  cb->unChecked().connect( boost::bind( &ShieldingSelect::isotopeUnCheckedCallback, this, iso ) );
  cb->massFractionChanged().connect( boost::bind( &ShieldingSelect::handleIsotopicChange, this, _1, iso ) );

  handleIsotopicChange( static_cast<float>(massFrac), iso );
  

/*  //Commenting out since I'm guessing most elements have at least 2 isotopes
  //Make sure the material has the isotope reequested to add
  const vector< Material::NuclideFractionPair > &nuclides = mat->nuclides;
  const vector< Material::ElementFractionPair > &elements = mat->elements;

  //Now make sure that there is less than 2 isotopes, then dont let mass
  //  fraction be editable.
  int isoCount = 0;
  for( const Material::NuclideFractionPair &nfp : nuclides )
    isoCount += (nfp.first && (nfp.first->atomicNumber==iso->atomicNumber));
  if( isoCount < 2 )
  {
    const SandiaDecay::Element *el = db->element( iso->atomicNumber );
    isoCount = static_cast<int>( db->nuclides( el ).size() );
  }//if( isoCount < 2 )
  cb->m_massFraction->disable();
  cb->m_massFraction->hide();
*/


  if( m_asSourceCBs->isHidden() )
    m_asSourceCBs->show();
}//void modelNuclideAdded( const std::string &symol )


void ShieldingSelect::handleIsotopicChange( float fraction, const SandiaDecay::Nuclide *nuc )
{
/*
 *handleIsotopicChange(...): This functions is a bit long-winded, and should
 *  probably be refactored or re-gone-through at some point.  Its fairly
 *  computationally inefiecient (_lots_ of loops), but I'm _guessing_ that
 *  this function still isnt the bottleneck in changing mass fractions, but
 *  instead its probably the bandwidth and user rendering, but I still havent
 *  benchmarked this funtion.
*/

  
  
  if( fraction < 0.0f )
    fraction = 0.0f;
  else if( fraction > 1.0f )
    fraction = 1.0f;

  std::shared_ptr<Material> mat = material();
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  if( !mat || !nuc || !db )
    return;

  const SandiaDecay::Element *element = db->element( nuc->atomicNumber );

  ElementToNuclideMap::iterator isos = m_sourceIsotopes.find(element);
  if( isos == m_sourceIsotopes.end() )
  {
    cerr << "\nShieldingSelect::handleIsotopicChange(...)\n\tSerious programming error I will ignore"  << endl;
    return;
  }//if( isos == m_sourceIsotopes.end() )

  vector< Material::NuclideFractionPair > &nuclides = mat->nuclides;
  vector< Material::ElementFractionPair > &elements = mat->elements;

  //First, go through and make sure that the material specifiecs the nuclide in
  //  its formula.  If not, make sure its element is specified in the formula,
  //  and then change from specifieng the the element to each of its individual
  //  isotopes, according to its natural composition - we will insert all
  //  isotopes into the 'nuclides' vector, even if its mass fraction will be 0.0
  bool hasIsotope = false, hasElement = false;
  for( const Material::NuclideFractionPair &nfp : nuclides )
    hasIsotope |= (nfp.first==nuc);
  for( const Material::ElementFractionPair &efp : elements )
    hasElement |= (efp.first==element);

  if( !hasIsotope && !hasElement )
    throw runtime_error( "ShieldingSelect::handleIsotopicChange(...): ran into"
                         " unexpected error" );

  //If the material has the element coorespoding to 'nuc', we'll transfer all
  //  that material to isotopes, so we can vary the mass fractions of the
  //  isotopes - and stuff
  if( hasElement )
  {
    double elMassFrac = 0.0;
    vector< Material::ElementFractionPair > newElements;
    for( const Material::ElementFractionPair &efp : elements )
    {
      if(efp.first==element)
        elMassFrac += efp.second;
      else
        newElements.push_back( efp );
    }//for(...)

    newElements.swap( elements );

    //Now get the natural abundance of the isotopes, and for all other isotopes
    //  set there mass fraction as 0.0
    typedef map<const SandiaDecay::Nuclide *, double> NucToAbundanceMap;
    NucToAbundanceMap nucAbunMap;
    const vector<const SandiaDecay::Nuclide *> nucs = db->nuclides( element );
    for( const SandiaDecay::Nuclide *n : nucs )
      nucAbunMap[n] = 0.0;

    bool hasNaturalAbundance = false;
    for( const SandiaDecay::NuclideAbundancePair &nap : element->isotopes )
    {
      hasNaturalAbundance |= (nap.abundance!=0.0);
      nucAbunMap[nap.nuclide] = nap.abundance*elMassFrac;
    }

    //If we didnt have a natural abundance, then just assume all isotopes are
    //  equally probable
    if( !hasNaturalAbundance )
    {
      const vector<const SandiaDecay::Nuclide *> nucs = db->nuclides( element );

      for( const SandiaDecay::Nuclide *nuc : nucs )
        nucAbunMap[nuc] = elMassFrac / nucs.size();
    }//if( !hasNaturalAbundance )

    //Now add to 'nuclides' all the isotopes for element
    for( const NucToAbundanceMap::value_type &vt : nucAbunMap )
      nuclides.push_back( make_pair(vt.first, static_cast<float>(vt.second)) );

    //Now we'll go through and consolidate all same isotopes.
    //  We could probably skip this, but I'll leave in - for the edge case
    //  were a user defines a material with both the element, and one or more
    //  of its isotopes.
    typedef map<const SandiaDecay::Nuclide *,float> NucToCoefMap;
    NucToCoefMap nucCoefs;
    for( const Material::NuclideFractionPair &nfp : nuclides )
    {
      NucToCoefMap::iterator pos = nucCoefs.find(nfp.first);
      if( pos == nucCoefs.end() )
        nucCoefs[nfp.first] = nfp.second;
      else
        pos->second += nfp.second;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )

    if( nucCoefs.size() != nuclides.size() )
    {
      nuclides.clear();
      for( const NucToCoefMap::value_type &vt : nucCoefs )
      {
        nuclides.push_back( make_pair(vt.first,vt.second) );
      }
    }//if( nucCoefs.size() != nuclides.size() )
  }//if( hasElement )

  //Now that were here, we are garunteed that 'nuc' is somewhere in 'nuclides'
  bool nucIsSrc = false;
  set<const SandiaDecay::Nuclide *> srcNucs, visibleNucs, allNucs;
  for( WWidget *child : isos->second->children() )
  {
    SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
    if( !cb )
      continue;

    if( cb->isotope() )
      visibleNucs.insert( cb->isotope() );

    if( cb->useAsSource() && cb->isotope() )
      srcNucs.insert( cb->isotope() );
    if( cb->isotope() == nuc )
      nucIsSrc |= cb->useAsSource();
  }//for( WWidget *child : children )


  double elMassFrac = 0.0, origNucFrac = 0.0,
      srcNucsFrac = 0.0, visibleNucFrac = 0.0;
  for( const Material::NuclideFractionPair &nfp : nuclides )
  {
    if( nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
      allNucs.insert( nfp.first );

    if( nfp.first==nuc )
      origNucFrac += nfp.second;
    if( nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
      elMassFrac += nfp.second;
    if( srcNucs.count(nfp.first) )
      srcNucsFrac += nfp.second;
    if( visibleNucs.count(nfp.first) )
      visibleNucFrac += nfp.second;
  }//for( const Material::NuclideFractionPair &nfp : nuclides )

  origNucFrac /= elMassFrac;
  srcNucsFrac /= elMassFrac;
  visibleNucFrac /= elMassFrac;

  const double nonVisibleFrac = 1.0 - visibleNucFrac;
  const double nonSrcFrac = 1.0 - srcNucsFrac - (nucIsSrc ? 0.0 : origNucFrac);
  const double fracDiff = fraction - origNucFrac;

  if( fabs(origNucFrac-fraction) < 0.000001 )
  {
    updateMassFractionDisplays( mat );
    return;
  }//if( fabs(origNucFrac-fraction) < 0.00001 )

  //Alright, were gonna go through a bunch of logic so that when the user
  //  changes mass fractions, the other mass fractions adjust roughly as I would
  //  like them to, and remain consistent (e.g. total mass fraction is 1.0)
  if( (fracDiff < 0.0) && fabs(nonVisibleFrac)>0.0001 )
  {
    //User lowered mass fraction, and there are some isotopes for this element
    //  which are not visible to the user (e.g. they are not fitting for these
    //  isotopes, at all).
    const double nonVisibleSF = 1.0 - fracDiff/nonVisibleFrac;

    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( !visibleNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second *= nonVisibleSF;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }else if( fracDiff < 0.0 && allNucs.size()!=visibleNucs.size() )
  {
    //User lowered mass fraction, but there is some isotopes for this element
    //  that are not visible, but have (at least near) zero contibution to the
    //  element.
    if( visibleNucs.size() > allNucs.size() )
      throw runtime_error( "ShieldingSelect::handleIsotopicChange():"
                           " Invalid material nuclide state." );

    const size_t nNonVisible = allNucs.size() - visibleNucs.size();
    const double nonVisibleFrac = -fracDiff/nNonVisible;

    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( !visibleNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second += elMassFrac*nonVisibleFrac;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }else if( fracDiff < 0.0 )
  {
    //User lowered mass fraction, and all isotopes for this element are visible
    double origFrac = 1.0 - srcNucsFrac;
    if( nucIsSrc )
      origFrac -= origNucFrac;
    const double wantedFrac = 1.0 - fraction;

    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( srcNucs.count(nfp.first) )
        nfp.second *= wantedFrac/origFrac;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }else if( fracDiff < nonVisibleFrac )
  {
    //User increased mass fraction, but small enough so we can take this from
    //  the isotopes not visible to the user
    const double nonVisibleSF = (nonVisibleFrac - fracDiff) / nonVisibleFrac;

    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( !visibleNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second *= nonVisibleSF;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }else if( fracDiff < (nonVisibleFrac+nonSrcFrac) )
  {
    //User increased mass fraction, but enough so we have to both take this from
    //  the isotopes not visible, as well as the non-source visible isotopes
//    const double diffNeeded = fracDiff - nonVisibleFrac;
//    const double srcFracMult = (nonSrcFrac - diffNeeded) / nonSrcFrac;

    double nonNucSrcFrac = 0.0;
    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first!=nuc && !srcNucs.count(nfp.first)
          && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nonNucSrcFrac += nfp.second;
      else if( nfp.first==nuc )
        nonNucSrcFrac -= (elMassFrac*fraction-nfp.second);
    }//for( Material::NuclideFractionPair &nfp : nuclides )
    
    //nonNucSrcFrac may be slightly below zero here - it might just be due to
    //  float roundoff, but I'm not entirely sure
    
    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( !srcNucs.count(nfp.first) && nfp.second!=0.0
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
      {
        if( nonNucSrcFrac > 0.0 )
        {
          cerr << "ShieldingSelect::handleIsotopicChange(...): I dont thing the mass fraction statment is corect!" << endl;
          nfp.second = elMassFrac*(1.0-fraction)*nfp.second/nonNucSrcFrac;
        }else
          nfp.second = 0.0;
      }else if( !visibleNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second = 0.0;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }else
  {
    //User increased mass fraction, but enough so we have to both take this from
    //  the isotopes not visible, as well as the non-source visible isotopes,
    //  as well as the other source isotopes
    const double diffNeeded = fracDiff - nonVisibleFrac - nonSrcFrac;

    double otherSrcFrac = srcNucsFrac;
    if( nucIsSrc )
      otherSrcFrac -= origNucFrac;

    const double srcFracMult = (otherSrcFrac - diffNeeded) / otherSrcFrac;

    for( Material::NuclideFractionPair &nfp : nuclides )
    {
      if( nfp.first==nuc )
        nfp.second = elMassFrac*fraction;
      else if( srcNucs.count(nfp.first) )
        nfp.second *= srcFracMult;
      else if( !srcNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second = 0.0;
      else if( !visibleNucs.count(nfp.first)
               && nfp.first && (nfp.first->atomicNumber==nuc->atomicNumber) )
        nfp.second = 0.0;
    }//for( const Material::NuclideFractionPair &nfp : nuclides )
  }//if() / else to figure out how to deal with mass fraction change.


/*
  //Code to test that the total mass fraction still adds up to 1.0.
  //  I have noticed a _small_ amount of leaking (ex 1.00366 instead of 1.0).
  double totalMassFrac = 0.0;
  for( Material::NuclideFractionPair &nfp : nuclides )
    totalMassFrac += nfp.second;
  for( const Material::ElementFractionPair &efp : elements )
    totalMassFrac += efp.second;
  cerr << "totalMassFrac=" << totalMassFrac << endl;
*/

  updateMassFractionDisplays( mat );
  
  //Need to make sure mass fraction for the element passed in is at most 1.0

  materialModified().emit( this );
}//void handleIsotopicChange( double fraction, const SandiaDecay::Nuclide *nuc )



void ShieldingSelect::updateMassFractionDisplays( std::shared_ptr<const Material> mat )
{
  if( !mat )
    return;

  //Lets go through and update the values displayed to the user of the
  //  isotopics
  for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  {
    int nisos = 0;
    double frac_accounted_for = 0.0;
    
    WText *otherfractxt = 0;
    const vector<WWidget *> children = vt.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb )
      {
        if( child->hasStyleClass( "MassFracNoFit" ) )
          otherfractxt = dynamic_cast<WText *>( child );
        continue;
      }//if( !cb )

      double massFrac = 0.0;
      try
      {
        massFrac = massFractionOfElement( cb->isotope(), mat );
        ++nisos;
      }catch(...)  //hopefully this never happens
      {
        passMessage( "There has been an unexpected internal error in"
                     " dealing with the material change, you may want to"
                     " try re-selecting the material, as well as checking the"
                     " calculation log before trusting the results.",
                     "",
                     WarningWidget::WarningMsgHigh );
        cerr << endl << "ShieldingSelect::updateMassFractionDisplays(...)\n\tSerious programming error here"
             << endl << endl;
      }//try/catch

      frac_accounted_for += massFrac;
      
      if( fabs(cb->massFraction()-massFrac) > 0.000001 )
        cb->setMassFraction( massFrac );
    }//for( WWidget *child : children )
    
    if( !nisos || (otherfractxt != children.back()) )
    {
      if( otherfractxt )
        delete otherfractxt;
      otherfractxt = 0;
    }
    
    if( nisos && !otherfractxt )
    {
      otherfractxt = new WText( vt.second );
      otherfractxt->addStyleClass( "MassFracNoFit" );
      otherfractxt->setToolTip( "If mass fractions are fit for, their fractional"
                                " sum will remain constant in the fit." );
    }
    
    if( otherfractxt )
    {
      if( fabs(1.0-frac_accounted_for) < 1.0E-6 )
      {
        //Pretty much all the mass of this element is accounted for, so dont
        //  display this text
        delete otherfractxt;
        otherfractxt = 0;
      }else
      {
        char buffer[128];
        const bool fit = (m_fitMassFrac && m_fitMassFrac->isVisible() && m_fitMassFrac->isChecked());
        const double percent_other = 100.0*(1.0-frac_accounted_for);
        const char *elsym = vt.first->symbol.c_str();
      
        snprintf( buffer, sizeof(buffer), "Assuming %s%.2f%% other %s isos",
                  (fit ? "fixed " : ""), percent_other, elsym );
        otherfractxt->setText( buffer );
      }
    }//if( otherfractxt )
  }//for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  
  updateIfMassFractionCanFit();
}//void updateMassFractionDisplays()


void ShieldingSelect::removeUncertFromDistanceEdit( Wt::WLineEdit *edit )
{
  if( !edit )
    return;
  
  //Get rid of the uncertainty text if the value is edited; this currently
  //  will cause the uncertainty to be removed
  string thickstr = edit->text().toUTF8();
  
  SpecUtils::trim( thickstr );
  
  if( thickstr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
  {
    thickstr += " cm";
    edit->setText( thickstr );
  }
  
  const size_t open_pos = thickstr.find( '(' );
  if( open_pos == string::npos )
    return;
  
  const size_t close_pos = thickstr.find( ')', open_pos );
  if( close_pos == string::npos )
    return;
  
  thickstr.erase( thickstr.begin()+open_pos, thickstr.begin()+close_pos+1 );
  size_t pos;
  while( (pos = thickstr.find( "  " )) != string::npos )
    thickstr.erase( thickstr.begin()+pos, thickstr.begin()+pos+1 );
  
  edit->setText( thickstr );
}//void ShieldingSelect::removeUncertFromThickness()


void ShieldingSelect::handleToggleGeneric()
{
  m_isGenericMaterial = (!m_isGenericMaterial);
  
  
  if( m_isGenericMaterial )
  {
    const string oldmaterial = m_materialEdit->text().toUTF8();
    
    //See if we can convert the current material into AN, AD
    m_dimensionsStack->setCurrentWidget( m_genericDiv );
    m_materialSummary->setText( "" );
    m_materialEdit->setText( "Generic" );
    m_materialEdit->disable();
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/atom_black.png") );
    
    const Material *mat = nullptr;
    try
    {
      mat = m_materialDB->material( oldmaterial );
    }catch(std::exception &)
    {}
    
    double ad = -1.0;
    if( mat )
    {
      const string thick = m_thicknessEdit->text().toUTF8();
      try
      {
        ad = mat->density * PhysicalUnits::stringToDistance(thick);
      }catch( std::exception & )
      {
      }
    }//if( mat )
    
    if( ad >= 0.0 && mat )
    {
      char an_buffer[32], ad_buffer[32];
      snprintf( an_buffer, sizeof(an_buffer), "%.2f", mat->massWeightedAtomicNumber() );
      snprintf( ad_buffer, sizeof(ad_buffer), "%.2f", ad*PhysicalUnits::cm2/PhysicalUnits::g );
      
      m_atomicNumberEdit->setText( an_buffer );
      m_arealDensityEdit->setText( ad_buffer );
    }else
    {
      m_atomicNumberEdit->setText( "26" );
      m_arealDensityEdit->setText( "0.0" );
    }
  }else
  {
    m_toggleImage->setImageLink(Wt::WLink("InterSpec_resources/images/shield.png"));
    m_materialEdit->enable();
    displayInputsForCurrentGeometry();
    
    string aNstr = SpecUtils::trim_copy( m_atomicNumberEdit->text().toUTF8() );
    string aDstr = SpecUtils::trim_copy( m_arealDensityEdit->text().toUTF8() );
    
    WLineEdit *dist_edit = nullptr;
    switch( m_geometry )
    {
      case GeometryType::Spherical:      dist_edit = m_thicknessEdit; break;
      case GeometryType::CylinderEndOn:  dist_edit = m_cylLengthEdit; break;
      case GeometryType::CylinderSideOn: dist_edit = m_cylRadiusEdit; break;
      case GeometryType::Rectangular:    dist_edit = m_rectDepthEdit; break;
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    assert( dist_edit );
    
    if( aNstr.empty() || aDstr.empty()
        || std::atof(aDstr.c_str())<=0.0
        || std::atof(aNstr.c_str())<=0.0 )
    {
      m_materialEdit->setText( "" );
      dist_edit->setText( "1.0 cm" );
    }else
    {
      //Get the atomic number closest to what the generic material was
      int atomicNum = (m_forFitting ? 16 : -1);
      double ad = -1;
      try
      {
        atomicNum = static_cast<int>( 0.5+std::stod(aNstr) );
        ad = std::stod(aDstr) * PhysicalUnits::g / PhysicalUnits::cm2;
      }catch(...){};
      
      if( !m_forFitting )
      {
        if( atomicNum > 97 )
          atomicNum = -1;
      }else
      {
        atomicNum = max( atomicNum, 1 );
        atomicNum = min( atomicNum, 97 );
      }
      
      const SandiaDecay::Element *el = 0;
      
      if( atomicNum >= 1 )
      {
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        el = db->element( atomicNum );
      }
      
      if( el )
      {
        m_materialEdit->setText( el->symbol );
        const Material *mat = m_materialDB->material( el->symbol );
        if( ad >= 0.0 )
        {
          const double dist = ad / mat->density;
          const string diststr = PhysicalUnits::printToBestLengthUnits( dist );
          dist_edit->setText( diststr );
        }else
        {
          dist_edit->setText( "0 cm" );
        }
      }else
      {
        m_materialEdit->setText( "" );
        dist_edit->setText( "1 cm" );
      }
    }//if( input is empty ) / else
  }//if( m_isGenericMaterial ) / else
  
  
  handleMaterialChange();
}//void ShieldingSelect::handleToggleGeneric()


void ShieldingSelect::displayInputsForCurrentGeometry()
{
  WContainerWidget *displayWidget = nullptr;
  if( m_isGenericMaterial )
  {
    displayWidget = m_genericDiv;
  }else
  {
    switch( m_geometry )
    {
      case GeometryType::Spherical:       displayWidget = m_sphericalDiv;   break;
      case GeometryType::CylinderEndOn:   displayWidget = m_cylindricalDiv; break;
      case GeometryType::CylinderSideOn:  displayWidget = m_cylindricalDiv; break;
      case GeometryType::Rectangular:     displayWidget = m_rectangularDiv; break;
      case GeometryType::NumGeometryType: break;
    }//switch( m_geometry )
  }//if( is generic ) / else
  
  assert( displayWidget );
  if( !displayWidget )
    throw runtime_error( "ShieldingSelect: invalid m_geometry." );
  
  m_dimensionsStack->setCurrentWidget( displayWidget );
}//void displayInputsForCurrentGeometry()


void ShieldingSelect::handleMaterialChange()
{
  typedef pair<const SandiaDecay::Element *,float> ElementFrac;
  typedef pair<const SandiaDecay::Nuclide *,float> NuclideFrac;

  std::shared_ptr<Material> newMaterial;
  std::shared_ptr<Material> previousMaterial = m_currentMaterial;
  
  displayInputsForCurrentGeometry();
  setTraceSourceMenuItemStatus();
  
  if( m_isGenericMaterial )
  {
    m_materialSummary->setText( "" );
    m_materialEdit->setText( "Generic" );
    m_materialEdit->disable();
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/atom_black.png") );
    
    
    if( m_traceSources )
    {
      for( WWidget *w : m_traceSources->children() )
      {
        TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
        assert( src );
        handleTraceSourceWidgetAboutToBeRemoved( src );
      }//for( WWidget *w : traceSources )
      
      m_traceSources->clear();
    }//if( m_traceSources )
  }else
  {
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/shield.png") );
    m_materialEdit->enable();
    
    
    string tooltip = "nothing";
    char summary[128];
    summary[0] = '\0';
    
    newMaterial = material();
    
    if( !!newMaterial )
    {
//      if( SpecUtils::iequals_ascii(newMaterial->name, "void") )
//      {
//        if( m_forFitting )
//        {
//          m_fitThicknessCB->setChecked( false );
//          m_fitThicknessCB->hide();
//        }//if( m_forFitting )
//        tooltip = "";
//      }else
      {
        if( m_forFitting )
          m_fitThicknessCB->show();

        const double density = newMaterial->density
                               * PhysicalUnits::cm3 / PhysicalUnits::gram;
        
        if( m_forFitting )
        {
          const float effAtomicNumber = newMaterial->massWeightedAtomicNumber();
          snprintf( summary, sizeof(summary),
                    "&rho;=%.2g g/cm<sup>3</sup>, <span style=\"text-decoration:overline\">AN</span>&#126;%.1f",
                    density, effAtomicNumber );
        }else
        {
          snprintf( summary, sizeof(summary),
                    "&rho;=%.2g g/cm<sup>3</sup>", density );
        }//if( m_forFitting ) / else
        
        
        
        tooltip += newMaterial->name + " consist of (mass fraction, element):\n";

        for( const ElementFrac &ef : newMaterial->elements )
        {
          if( ef.first )
          {
            char buffer[256];
            snprintf( buffer, sizeof(buffer), "%.4f %s\n",
                      ef.second, ef.first->name.c_str() );
            tooltip += buffer;
          }
        }//for( const ElementFrac &ef : newMaterial->elements )

        for( const NuclideFrac &ef : newMaterial->nuclides )
        {
          if( ef.first )
          {
            char buffer[256];
            snprintf( buffer, sizeof(buffer), "\t%.4f %s\n",
                      ef.second, ef.first->symbol.c_str() );
            tooltip += buffer;
          }
        }//for( const NuclideFrac &ef : newMaterial->nuclides )

        //Could consider putting attenuiation coefficients here...
      }//if( newMaterial->name == "void" )
    }else if( m_materialEdit->text().narrow().length() )
    {
      snprintf( summary, sizeof(summary), "invalid mat." );
      WStringStream msg;
      msg << "'" << m_materialEdit->text().toUTF8()
          << "' is not a valid material";
      passMessage( msg.str(), "", WarningWidget::WarningMsgInfo );
    }//if( material ) / else
    

//NOTE: can't add tooltip to this, causes WT error when toggling.  Can't fix.
//    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
//    const bool showToolTips = true;//InterSpecUser::preferenceValue<bool>( "ShowTooltips", app->viewer() );
//    HelpSystem::attachToolTipOn( this,tooltip, showToolTips );
    
    m_materialSummary->setText( summary );
    
    if( m_traceSources )
    {
      const vector<WWidget *> &traceSources = m_traceSources->children();
      for( WWidget *w : traceSources )
      {
        TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
        assert( src );
        if( !src || !src->nuclide() )
          continue;
        
#if( PERFORM_DEVELOPER_CHECKS )
        // Make sure we havent messed the logic up and let the same nuclide be both a
        //  self-attenuating source, and a trace source.
        for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
        {
          const vector<WWidget *> children = vt.second->children();
          for( WWidget *child : children )
          {
            SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
            if( cb && cb->useAsSource() )
            {
              assert( cb->isotope() != src->nuclide() );
            }
          }//for( WWidget *child : children )
        }//for(...)
#endif
        
        src->updateForMaterialChange();
        m_activityFromVolumeNeedUpdating.emit( this, src->nuclide() );
      }//for( WWidget *w : traceSources )
    }//if( m_traceSources )
  }//if( generic material ) / else
  
  
  //Now we need to update the activities for any isotopes that are
  if( !!newMaterial && m_asSourceCBs && (previousMaterial == newMaterial) )
  {
    updateMassFractionDisplays( newMaterial );

    for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
    {
      const vector<WWidget *> children = vt.second->children();
      for( WWidget *child : children )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->useAsSource() )
          m_activityFromVolumeNeedUpdating.emit( this, cb->isotope() );
      }//for( WWidget *child : children )
    }//for(...)
  }//if( (previousMaterial == newMaterial) && m_asSourceCBs )


  if( (previousMaterial != newMaterial) && m_asSourceCBs )
  {
    for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
    {
      const vector<WWidget *> children = vt.second->children();
      for( WWidget *child : children )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->useAsSource() )
          removingIsotopeAsSource().emit( cb->isotope(), ModelSourceType::Intrinsic );
      }//for( WWidget *child : children )

      delete vt.second;
    }//for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )

    m_sourceIsotopes.clear();

    if( newMaterial && m_sourceModel )
    {
      const int nrow = m_sourceModel->rowCount();
      for( int row = 0; row < nrow; ++row )
      {
        const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
        modelNuclideAdded( nuc );
      }//for( int row = 0; row < nrow; ++row )
    }//if( newMaterial )

    if( m_sourceIsotopes.size() )
      m_asSourceCBs->show();
    else
      m_asSourceCBs->hide();
  }//if( previousMaterial != newMaterial )

  
#if( PERFORM_DEVELOPER_CHECKS )
  for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  {
    const vector<WWidget *> children = vt.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb || !newMaterial || !cb->isotope() )
        continue;
      
      double massFrac = 0.0;
      try
      {
        massFrac = massFractionOfElement( cb->isotope(), newMaterial );
      }catch(...)  //hopefully this never happens
      {
        stringstream msg;
        msg << "Failed to get massFractionOfElement " << cb->isotope()->symbol
            << " from " << newMaterial->name;
        log_developer_error( __func__, msg.str().c_str() );
      }//try/catch
      
      if( fabs(cb->massFraction()-massFrac) > 0.00001 )
      {
        stringstream msg;
        msg << "Mass cb fraction for " << cb->isotope()->symbol << " in "
            << newMaterial->name << " was " << cb->massFraction()
            << " but was expected to be " << massFrac << " from material.";
        log_developer_error( __func__, msg.str().c_str() );
      }
    }//for( WWidget *child : children )
  }//for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
#endif
  
  cerr << "\nShieldingSelect::handleMaterialChange()\n\tShould remove this call to "
       << "updateMassFractionDisplays(...) its verified the developer checks always pass\n" << endl;
  updateMassFractionDisplays( newMaterial );

  
  if( previousMaterial != newMaterial )
    m_materialModifiedSignal.emit( this );
  else
    m_materialChangedSignal.emit( this );
}//void handleMaterialChange()



void ShieldingSelect::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  const char *name, *value;
  rapidxml::xml_node<char> *base_node, *node;
  rapidxml::xml_attribute<char> *attr;
  
  name = "Shielding";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  string versionstr = std::to_string(ShieldingSelect::sm_xmlSerializationMajorVersion)
                      + "." + std::to_string(ShieldingSelect::sm_xmlSerializationMinorVersion);
  value = doc->allocate_string( versionstr.c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );

  name = "Geometry";
  value = GammaInteractionCalc::to_str(m_geometry);
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );
  
  name = "ForFitting";
  value = m_forFitting ? "1" : "0";
  node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( node );

  if( m_isGenericMaterial )
  {
    rapidxml::xml_node<> *generic_node;
    
    name = "Generic";
    generic_node = doc->allocate_node( rapidxml::node_element, name );
    base_node->append_node( generic_node );
    
    name = "ArealDensity";
    value = doc->allocate_string( m_arealDensityEdit->valueText().toUTF8().c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    generic_node->append_node( node );
    if( m_forFitting )
    {
      value = m_fitArealDensityCB->isChecked() ? "1" : "0";
      attr = doc->allocate_attribute( "Fit", value );
      node->append_attribute( attr );
    }//if( m_fitArealDensityCB )
    
    name = "AtomicNumber";
    value = doc->allocate_string( m_atomicNumberEdit->valueText().toUTF8().c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    generic_node->append_node( node );
    if( m_forFitting )
    {
      value = m_fitAtomicNumberCB->isChecked() ? "1" : "0";
      attr = doc->allocate_attribute( "Fit", value );
      node->append_attribute( attr );
    }//if( m_fitAtomicNumberCB )
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    auto addTruth = [doc,generic_node]( const char *truthName, const boost::optional<double> &value ){
      if( value )
      {
        const string strval = std::to_string(*value);
        const char *value = doc->allocate_string( strval.c_str() );
        rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, value );
        generic_node->append_node( node );
      }
    };//addTruth(...)
    addTruth( "TruthAD", truthAD );
    addTruth( "TruthADTolerance", truthADTolerance );
    addTruth( "TruthAN", truthAN );
    addTruth( "TruthANTolerance", truthANTolerance );
#endif
  }else
  {
    rapidxml::xml_node<> *material_node, *mass_frac_node, *iso_node;
    
    name = "Material";
    material_node = doc->allocate_node( rapidxml::node_element, name );
    base_node->append_node( material_node );
    
    name = "Name";
    value = doc->allocate_string( m_materialEdit->valueText().toUTF8().c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    material_node->append_node( node );
    
    //Lambda to add in dimension elements
    auto addDimensionNode = [this,doc,material_node]( const char *name, WLineEdit *edit, WCheckBox *cb )
                                                     -> rapidxml::xml_node<char> * {
      const char *value = doc->allocate_string( edit->valueText().toUTF8().c_str() );
      auto node = doc->allocate_node( rapidxml::node_element, name, value );
      material_node->append_node( node );
      if( m_forFitting )
      {
        value = cb->isChecked() ? "1" : "0";
        auto attr = doc->allocate_attribute( "Fit", value );
        node->append_attribute( attr );
      }//if( m_forFitting )
      
      return node;
    };//addDimensionNode lambda
    
    
    // For backward compatibility with XML serialization version 0.0, we will add in a <Thickness>
    //  element, corresponding to dimension along detector axis
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        node = addDimensionNode( "Thickness", m_thicknessEdit, m_fitThicknessCB );
        node->append_attribute( doc->allocate_attribute( "Remark", "CylinderRadiusThickness" ) );
        break;
        
      case GeometryType::CylinderEndOn:
        node = addDimensionNode( "Thickness", m_cylLengthEdit, m_fitCylLengthCB );
        node->append_attribute( doc->allocate_attribute( "Remark", "CylinderLengthThickness" ) );
        break;
        
      case GeometryType::CylinderSideOn:
        node = addDimensionNode( "Thickness", m_cylRadiusEdit, m_fitCylRadiusCB );
        node->append_attribute( doc->allocate_attribute( "Remark", "CylinderRadiusThickness" ) );
        break;
        
      case GeometryType::Rectangular:
        node = addDimensionNode( "Thickness", m_rectDepthEdit, m_fitRectDepthCB );
        node->append_attribute( doc->allocate_attribute( "Remark", "RectangularDepthThickness" ) );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    // We could write dimensions from all the WLineEdits into the XML, which would kinda keep state
    //  across changing geometries, but maybe for the moment we'll just write the relevant
    //  dimensions.  TODO: think about this a little more
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        addDimensionNode( "SphericalThickness", m_thicknessEdit, m_fitThicknessCB );
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        addDimensionNode( "CylinderRadiusThickness", m_cylRadiusEdit, m_fitCylRadiusCB );
        addDimensionNode( "CylinderLengthThickness", m_cylLengthEdit, m_fitCylLengthCB );
        break;
        
      case GeometryType::Rectangular:
        addDimensionNode( "RectangularWidthThickness",  m_rectWidthEdit,  m_fitRectWidthCB );
        addDimensionNode( "RectangularHeightThickness", m_rectHeightEdit, m_fitRectHeightCB );
        addDimensionNode( "RectangularDepthThickness",  m_rectDepthEdit,  m_fitRectDepthCB );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    
    if( m_forFitting )
    {
      name = "FitMassFraction";
      value = (m_fitMassFrac && m_fitMassFrac->isChecked()) ? "1" : "0";
      mass_frac_node = doc->allocate_node( rapidxml::node_element, name, value );
      material_node->append_node( mass_frac_node );
    }//if( m_forFitting )
    
    for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
    {
      for( WWidget *widget : etnm.second->children() )
      {
        SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
        
        if( src && src->useAsSource() && src->isotope() )
        {
          iso_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
          material_node->append_node( iso_node );
          
          value = doc->allocate_string( src->isotope()->symbol.c_str() );
          node = doc->allocate_node( rapidxml::node_element, "Name", value );
          iso_node->append_node( node );
          
          value = doc->allocate_string( std::to_string(src->massFraction()).c_str() );
          node = doc->allocate_node( rapidxml::node_element, "MassFrac", value );
          iso_node->append_node( node );
        }//if( src && src->useAsSource() )
      }//for( WWidget *widget : isotopeDiv->children() )
    }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
    
    
    if( m_traceSources )
    {
      for( WWidget *w : m_traceSources->children() )
      {
        const TraceSrcDisplay *src = dynamic_cast<const TraceSrcDisplay *>( w );
        assert( src );
        if( src && src->nuclide() )
          src->serialize( material_node );
      }//for( WWidget *w : m_traceSources->children() )
    }//if( m_traceSources )
    
    
    #if( INCLUDE_ANALYSIS_TEST_SUITE )
        auto addTruth = [doc,material_node]( const char *truthName, const boost::optional<double> &value ){
          if( value )
          {
            const string strval = PhysicalUnits::printToBestLengthUnits(*value,6);
            const char *value = doc->allocate_string( strval.c_str() );
            rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, value );
            material_node->append_node( node );
          }
        };//addTruth(...)
    
        addTruth( "TruthThickness", truthThickness );
        addTruth( "TruthThicknessTolerance", truthThicknessTolerance );
    #endif
  }//if( m_isGenericMaterial ) / else
}//void serialize( rapidxml::xml_document<> &doc ) const;


void ShieldingSelect::deSerialize( const rapidxml::xml_node<char> *shield_node )
{
  rapidxml::xml_attribute<char> *attr;
  rapidxml::xml_node<char> *node, *geom_node, *generic_node, *material_node;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
 
  if( shield_node->name() != string("Shielding") )
    throw runtime_error( "ShieldingSelects XML node should be 'Shielding'" );
  
  attr = shield_node->first_attribute( "version", 7 );
  int version;
  if( !attr || !attr->value() || !(stringstream(attr->value())>>version) )
    throw runtime_error( "ShieldingSelects should be versioned" );
  
  if( version != sm_xmlSerializationMajorVersion )
    throw runtime_error( "Invalid XML version for ShieldingSelect" );
  
  // Note that version is either "0" (implied minor version 0), or "0.1" or something; we could read
  //  in the minor version to sm_xmlSerializationMinorVersion and use it, but theres really no need.
  
  GeometryType geometry = GeometryType::Spherical;
  geom_node = shield_node->first_node( "Geometry", 8 );
  if( geom_node && geom_node->value() )
  {
    bool matched = false;
    const string geomstr( geom_node->value(), geom_node->value() + geom_node->value_size() );
    for( GeometryType type = GeometryType(0);
        !matched && (type != GeometryType::NumGeometryType);
        type = GeometryType( static_cast<int>(type) + 1 ) )
    {
      matched = (geomstr == GammaInteractionCalc::to_str(type));
      if( matched )
        geometry = type;
    }//for( loop over GeometryTypes )
    
    if( !matched )
      throw runtime_error( "Invalid geometry type: '" + geomstr + "'" );
  }//if( geom_node )
  
  
  bool forFitting;
  node = shield_node->first_node( "ForFitting", 10 );
  if( !node || !node->value() || !(stringstream(node->value())>>forFitting) )
    throw runtime_error( "Missing/invalid for fitting node" );

  if( m_forFitting != forFitting )
    throw runtime_error( "ShieldingSelect m_forFitting must be same as "
                         "XML being deserialized" );
  
  generic_node = shield_node->first_node( "Generic", 7 );
  material_node = shield_node->first_node( "Material", 8 );
  
  if( generic_node )
  {
    const rapidxml::xml_node<char> *ad_node, *an_node;
    ad_node = generic_node->first_node( "ArealDensity", 12 );
    an_node = generic_node->first_node( "AtomicNumber", 12 );
    
    if( !ad_node || !ad_node->value() || !ad_node || !ad_node->value() )
      throw runtime_error( "Generic material must have ArealDensity and"
                           " AtomicNumber nodes" );
    
    m_geometry = geometry;
    
    if( !m_isGenericMaterial )
      handleToggleGeneric();
    
    m_arealDensityEdit->setValueText( WString::fromUTF8(ad_node->value()) );
    m_atomicNumberEdit->setValueText( WString::fromUTF8(an_node->value()) );
    
    if( m_forFitting && m_fitArealDensityCB && m_fitAtomicNumberCB )
    {
      bool fit;
      attr = ad_node->first_attribute( "Fit", 3 );
      if( attr && attr->value() && (stringstream(attr->value())>>fit) )
        m_fitArealDensityCB->setChecked( fit );
      
      attr = an_node->first_attribute( "Fit", 3 );
      if( attr && attr->value() && (stringstream(attr->value())>>fit) )
        m_fitAtomicNumberCB->setChecked( fit );
    }//if( m_fitArealDensityCB )
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    auto getTruth = [generic_node]( const char *truthName, boost::optional<double> &value ){
      value.reset();
      auto node = generic_node->first_node( truthName );
      if( !node || !node->value() )
        return;
      
      double dblvalue;
      if( (stringstream(node->value()) >> dblvalue) )
        value = dblvalue;
      else
        cerr << "\n\nFailed to deserialize shielding " << truthName << " from " << node->value()
        << "\n\n" << endl;
    };//getTruth(...)
    
    getTruth( "TruthAD", truthAD );
    getTruth( "TruthADTolerance", truthADTolerance );
    getTruth( "TruthAN", truthAN );
    getTruth( "TruthANTolerance", truthANTolerance );
#endif
  }//if( generic_node )
  
  if( material_node )
  {
    bool fitMassFrac = false;
    vector<const SandiaDecay::Nuclide *> srcnuclides;
    const rapidxml::xml_node<> *frac_node, *iso_node, *name_node;
    const rapidxml::xml_node<> *dim_nodes[3] = { nullptr, nullptr, nullptr };
    WLineEdit *dim_edits[3] = { nullptr, nullptr, nullptr };
    WCheckBox *dim_cb[3] = { nullptr, nullptr, nullptr };
    
    
    name_node = material_node->first_node( "Name", 4 );
    if( !name_node || !name_node->value() )
      throw runtime_error( "Material node didnt have name node" );
    
    int required_dim = 0;
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        required_dim = 1;
        
        dim_nodes[0] = XML_FIRST_NODE(material_node, "SphericalThickness");
        if( !dim_nodes[0] )
          dim_nodes[0] = XML_FIRST_NODE(material_node, "Thickness");
        
        dim_edits[0] = m_thicknessEdit;
        dim_cb[0] = m_fitThicknessCB;
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        required_dim = 2;
        
        dim_nodes[0] = XML_FIRST_NODE(material_node, "CylinderRadiusThickness");
        dim_nodes[1] = XML_FIRST_NODE(material_node, "CylinderLengthThickness");
        
        dim_edits[0] = m_cylRadiusEdit;
        dim_edits[1] = m_cylLengthEdit;
        dim_cb[0]    = m_fitCylRadiusCB;
        dim_cb[1]    = m_fitCylLengthCB;
        break;
        
      case GeometryType::Rectangular:
        required_dim = 3;
        
        dim_nodes[0] = XML_FIRST_NODE(material_node, "RectangularWidthThickness");
        dim_nodes[1] = XML_FIRST_NODE(material_node, "RectangularHeightThickness");
        dim_nodes[2] = XML_FIRST_NODE(material_node, "RectangularDepthThickness");
        
        dim_edits[0] = m_rectWidthEdit;
        dim_edits[1] = m_rectHeightEdit;
        dim_edits[2] = m_rectDepthEdit;
        dim_cb[0]    = m_fitRectWidthCB;
        dim_cb[1]    = m_fitRectHeightCB;
        dim_cb[2]    = m_fitRectDepthCB;
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    // Check and make sure all the dimension node values are valid
    bool fit_dim[3] = { false, false, false };
    for( int i = 0; i < required_dim; ++i )
    {
      if( !dim_nodes[i] || !dim_nodes[i]->value() )
        throw runtime_error( "Missing required dimension node" );
      
      const string val( dim_nodes[i]->value(), dim_nodes[i]->value() + dim_nodes[i]->value_size() );
      try
      {
        const double dist = PhysicalUnits::stringToDistance(val);
        if( dist < 0.0 )
          throw runtime_error( "" );
      }catch( std::exception & )
      {
        throw runtime_error( "Invalid dimension given '" + val + "'" );
      }
      
      if( m_forFitting )
      {
        const auto attr = dim_nodes[i]->first_attribute( "Fit", 3 );
        
        if( !attr || !(stringstream(attr->value()) >> fit_dim[i]) )
          throw runtime_error( "Material node expected thickness Fit attribute" );
      }
    }//for( int i = 0; i < required_dim; ++i )
    
    
    const bool geom_changed = (m_geometry != geometry);
    
    // Now that we've mostly validated the XML, start actually changing the widgets values.
    m_geometry = geometry;
    
    if( m_isGenericMaterial )
      handleToggleGeneric();
    
    if( geom_changed )
      displayInputsForCurrentGeometry();
    
    const string material_name( name_node->value(), name_node->value() + name_node->value_size() );
    m_materialEdit->setValueText( WString::fromUTF8(material_name) );
    
    for( int i = 0; i < required_dim; ++i )
    {
      assert( dim_nodes[i] );
      assert( dim_edits[i] );
      
      const string val( dim_nodes[i]->value(), dim_nodes[i]->value() + dim_nodes[i]->value_size() );
      dim_edits[i]->setValueText( WString::fromUTF8(val) );
      if( m_forFitting )
      {
        assert( dim_cb[i] );
        dim_cb[i]->setChecked( fit_dim[i] );
      }
    }//for( int i = 0; i < required_dim; ++i )
    
    handleMaterialChange();
    
    if( m_forFitting )
    {
      const rapidxml::xml_node<> *fitmassfrac_node = XML_FIRST_NODE(material_node, "FitMassFraction");
      if( fitmassfrac_node && fitmassfrac_node->value() )
      {
        stringstream(fitmassfrac_node->value()) >> fitMassFrac;
      }//if( m_forFitting )
    }//if( m_forFitting )

    double last_frac = 0.0;
    const SandiaDecay::Nuclide *last_nuc = NULL;
    
    for( iso_node = material_node->first_node( "Nuclide", 7 );
        iso_node; iso_node = iso_node->next_sibling( "Nuclide", 7 ) )
    {
      name_node = iso_node->first_node( "Name", 4 );
      frac_node = iso_node->first_node( "MassFrac", 8 );
      
      if( !name_node || !name_node->value()
          || !frac_node || !frac_node->value() )
        throw runtime_error( "Missing invalid name/mass frac node form iso" );
      
      const SandiaDecay::Nuclide *nuc = db->nuclide( name_node->value() );
      if( !nuc )
        throw runtime_error( string(name_node->value()) + " is not a valid isotope" );
      
      srcnuclides.push_back( nuc );
      
      double fraction;
      if( !(stringstream(frac_node->value()) >> fraction) )
        throw runtime_error( "Invalid mass fraction: " + string(frac_node->value()) );
      
//      modelNuclideAdded( nuc );
//      const SandiaDecay::Element *el = db->element( nuc->atomicNumber );
//      if( m_sourceIsotopes.find( el ) == m_sourceIsotopes.end() )
//      {
//      }//if(...)
      try
      {
        setMassFraction( nuc, fraction );
        
        last_nuc = nuc;
        last_frac = fraction;
      }catch( std::exception &e )
      {
        cerr << "ShieldingSelect::deSerialize(...)\n\tCaught: " << e.what()
             << " but continuuing anyway" << endl;
      }//try / catch
    }//for( loop over isotope nodes )
    
    //Now set the check boxes to make all the source nuclides called out in the
    //  XML as actual src nuclides, since we only saved nuclides we actually
    //  wanted to use as source nuclides.
    if( m_currentMaterial )
    {
      for( const SandiaDecay::Nuclide *nuc : srcnuclides )
      {
        for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
        {
          for( WWidget *widget : etnm.second->children() )
          {
            SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
            if( src && (nuc == src->isotope()) )
              src->setUseAsSource( true );
          }
        }
      }//for( const SandiaDecay::Nuclide *nuc : srcnuclides )
      
      
      // Get rid of all the trace sources - I dont *think* we need to emit signals and all that since
      if( m_traceSources )
        m_traceSources->clear();
      
      for( auto trace_node = XML_FIRST_NODE(material_node, "TraceSource");
          trace_node; trace_node = XML_NEXT_TWIN(trace_node) )
      {
        addTraceSource();
        assert( m_traceSources );
        
        TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( m_traceSources->children().back() );
        assert( src );
        
        // When we de-serialize we assume the source fitting model already has the nuclides
        //  available to become trace sources
        if( src )  //JIC
          src->deSerialize( trace_node );
      }//for( loop over trace nodes )
    }//if( m_currentMaterial )
    
    if( m_fitMassFrac )
      m_fitMassFrac->setChecked( fitMassFrac );
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    auto getTruth = [material_node]( const char *truthName, boost::optional<double> &value ){
      value.reset();
      auto node = material_node->first_node( truthName );
      if( !node || !node->value() )
        return;
      
      try
      {
        value = PhysicalUnits::stringToDistance( node->value() );
      }catch( std::exception &e )
      {
        cerr << "\n\nFailed to deserialize a shielding " << truthName << ": " << e.what() << "\n\n";
      }
    };//getTruth(...)
    
    getTruth( "TruthThickness", truthThickness );
    getTruth( "TruthThicknessTolerance", truthThicknessTolerance );
#endif
    
    
    //Calling handleIsotopicChange(...) will perform some normailizations
    //  and stuff (I think), that setMassFraction(...) doesnt do
    if( last_nuc )
      handleIsotopicChange( static_cast<float>(last_frac), last_nuc );
  }//if( material_node )
}//void deSerialize( const rapidxml::xml_node<char> *shielding_node ) const


