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
#include <Wt/WServer>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WGroupBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
#include <Wt/WStackedWidget>
#include <Wt/WRegExpValidator>
#include <Wt/WSuggestionPopup>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/RapidXmlUtils.hpp"

#include "InterSpec/AppUtils.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace Wt;
using namespace std;

using GammaInteractionCalc::GeometryType;
using GammaInteractionCalc::TraceActivityType;



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
  
  /** This uncertainty is only displayed as the tool-tip of `m_activityInput`, and is only set from
   `setTraceSourceTotalActivity(...)`.  Everywhere else `m_currentTotalActivity`
   gets set, the uncertainty will be set to a negative value.
   A negative or zero value indicates no uncertainty.
   */
  double m_currentTotalActivityUncert;
  
  ShieldingSelect *m_parent;
  Wt::WComboBox *m_isoSelect;
  /** The tool tip on `m_activityInput` will convey the uncertatainty information. */
  Wt::WLineEdit *m_activityInput;
  Wt::WComboBox *m_activityType;
  Wt::WCheckBox *m_allowFitting;
  
  Wt::WContainerWidget *m_relaxationDiv;
  Wt::WLineEdit *m_relaxationDistance;
  Wt::WText *m_relaxationDescription;
  
  Wt::Signal<> m_userChanged;
  Wt::Signal<const SandiaDecay::Nuclide *,double> m_activityUpdated;
  Wt::Signal<TraceSrcDisplay *, const SandiaDecay::Nuclide * /* old nuclide */> m_nucChangedSignal;
  
  
public:
  TraceSrcDisplay( ShieldingSelect *parent )
  : WGroupBox( "Trace Source", parent->m_traceSources ),
    m_currentNuclide( nullptr ),
    m_currentDisplayActivity( 0.0 ),
    m_currentTotalActivity( 0.0 ),
    m_currentTotalActivityUncert( -1.0 ),
    m_parent( parent ),
    m_isoSelect( nullptr ),
    m_activityInput( nullptr ),
    m_activityType( nullptr ),
    m_allowFitting( nullptr ),
    m_relaxationDiv( nullptr ),
    m_relaxationDistance( nullptr ),
    m_relaxationDescription( nullptr ),
    m_userChanged( this ),
    m_activityUpdated( this ),
    m_nucChangedSignal( this )
  {
    wApp->useStyleSheet( "InterSpec_resources/GridLayoutHelpers.css" );
    
    const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    assert( !parent->isGenericMaterial() );
    addStyleClass( "TraceSrcDisplay" );
    
    WPushButton *closeIcon = new WPushButton( this );
    closeIcon->addStyleClass( "closeicon-wtdefault GridThirdCol GridFirstRow" );
    closeIcon->clicked().connect( boost::bind( &ShieldingSelect::removeTraceSourceWidget, m_parent, this) );
    closeIcon->setToolTip( "Remove this trace source." );
    
    
    WLabel *label = new WLabel( "Nuclide", this );
    label->addStyleClass( "GridFirstCol GridSecondRow" );
    
    m_isoSelect = new WComboBox( this );
    m_isoSelect->addItem( "Select" );
    m_isoSelect->activated().connect( this, &TraceSrcDisplay::handleUserNuclideChange );
    m_isoSelect->addStyleClass( "GridSecondCol GridSecondRow GridSpanTwoCol" );
    label->setBuddy( m_isoSelect );
    
    
    label = new WLabel( "Activity", this );
    label->addStyleClass( "GridFirstCol GridThirdRow" );
    
    
    m_activityInput = new WLineEdit( this );
    m_activityInput->addStyleClass( "GridSecondCol GridStretchCol GridThirdRow" );
    
    m_activityInput->setAutoComplete( false );
    m_activityInput->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_activityInput->setAttributeValue( "autocorrect", "off" );
    m_activityInput->setAttributeValue( "spellcheck", "off" );
#endif
    
    WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityRegex, m_activityInput );
    val->setFlags( Wt::MatchCaseInsensitive );
    m_activityInput->setValidator( val );
    m_activityInput->changed().connect( this, &TraceSrcDisplay::handleUserActivityChange );
    m_activityInput->enterPressed().connect( this, &TraceSrcDisplay::handleUserActivityChange );
    m_activityInput->setText( (useBq ? "37 MBq" : "1 mCi") );
    m_activityInput->setToolTip( WString() );
    m_currentTotalActivity = m_currentDisplayActivity = 0.001*PhysicalUnits::ci;
    m_currentTotalActivityUncert = -1.0;
    label->setBuddy( m_activityInput );
    
    m_activityType = new WComboBox( this );
    m_activityType->activated().connect( this, &TraceSrcDisplay::handleUserChangeActivityType );
    m_activityType->addStyleClass( "GridThirdCol GridThirdRow" );
    
    m_allowFitting = new WCheckBox( "Fit activity value", this );
    m_allowFitting->addStyleClass( "GridSecondCol GridStretchCol GridFourthRow GridSpanTwoCol CbNoLineBreak" );
    
    
    m_relaxationDiv = new WContainerWidget( this );
    m_relaxationDiv->addStyleClass( "GridFirstCol GridFifthRow GridSpanThreeCol RelaxDistDisplay" );
    
    WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, m_relaxationDiv );
    distValidator->setFlags( Wt::MatchCaseInsensitive );
    
    label = new WLabel( "Relaxation&nbsp;Distance", m_relaxationDiv );
    label->addStyleClass( "GridFirstCol GridFirstRow" );
    
    m_relaxationDistance = new WLineEdit( m_relaxationDiv );
    m_relaxationDistance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_relaxationDistance->setAttributeValue( "autocorrect", "off" );
    m_relaxationDistance->setAttributeValue( "spellcheck", "off" );
#endif
    m_relaxationDistance->setText( "1.0 cm" );
    m_relaxationDistance->setValidator( distValidator );
    m_relaxationDistance->addStyleClass( "GridSecondCol GridStretchCol GridFirstRow" );
    m_relaxationDistance->changed().connect( this, &TraceSrcDisplay::handleUserRelaxDistChange );
    m_relaxationDistance->enterPressed().connect( this, &TraceSrcDisplay::handleUserRelaxDistChange );
    
    m_relaxationDescription = new WText("", m_relaxationDiv);
    m_relaxationDescription->addStyleClass( "GridFirstCol GridSecondRow GridSpanTwoCol RelaxDescription" );
    m_relaxationDiv->hide();
    
    // Using WCheckBox::changed() instead of WCheckBox::checked()/unChecked() causes some odd
    //  issues in Wt 3.3.4 (at least) where m_allowFitting->isChecked() isnt always up to date in
    //  the immediate render cycle (it is in the immediate call to to the connected signal, but
    //  seemingly not in calls during render()... havent looked into this much yet.
    //m_allowFitting->changed().connect( this, &TraceSrcDisplay::handleUserChangeAllowFit );
    m_allowFitting->checked().connect( this, &TraceSrcDisplay::handleUserChangeAllowFit );
    m_allowFitting->unChecked().connect( this, &TraceSrcDisplay::handleUserChangeAllowFit );
    
    
    m_allowFitting->checked().connect( this, &TraceSrcDisplay::emitUserChaged );
    m_allowFitting->unChecked().connect( this, &TraceSrcDisplay::emitUserChaged );
    m_activityType->activated().connect( this, &TraceSrcDisplay::emitUserChaged );
    m_relaxationDistance->changed().connect( this, &TraceSrcDisplay::emitUserChaged );
    m_activityInput->changed().connect( this, &TraceSrcDisplay::emitUserChaged );
    m_isoSelect->activated().connect( this, &TraceSrcDisplay::emitUserChaged );
    
      
    updateAvailableActivityTypes();
    updateAvailableIsotopes();
  }//TraceSrcDisplay constructor
  
  
  void serialize( rapidxml::xml_node<char> * const parent_node ) const
  {
    const ShieldingSourceFitCalc::TraceSourceInfo trace = toTraceSourceInfo();
    trace.serialize( parent_node );
  }//void serialize( rapidxml::xml_node<> *parent )
  
  
  void deSerialize( const rapidxml::xml_node<char> *shielding_node )
  {
    ShieldingSourceFitCalc::TraceSourceInfo trace;
    trace.deSerialize( shielding_node );

#if( PERFORM_DEVELOPER_CHECKS )
    const ShieldingSourceFitCalc::TraceSourceInfo checktrace = toTraceSourceInfo();
    
    try
    {
      ShieldingSourceFitCalc::TraceSourceInfo::equalEnough( trace, checktrace );
    }catch( std::exception &e )
    {
      cerr << "Failed to roundtrip TraceSrcDisplay to TraceSourceInfo: " << e.what() << endl;
      assert( 0 );
    }
#endif
  }//void deSerialize( const rapidxml::xml_node<char> *shielding_node )
  
  
  ShieldingSourceFitCalc::TraceSourceInfo toTraceSourceInfo() const
  {
    ShieldingSourceFitCalc::TraceSourceInfo trace;
    trace.m_type = activityType();
    trace.m_fitActivity = allowFittingActivity();
    trace.m_nuclide = nuclide();
    trace.m_activity = displayActivity();
    if( trace.m_type == GammaInteractionCalc::TraceActivityType::ExponentialDistribution )
      trace.m_relaxationDistance = relaxationLength();
    
    return trace;
  }//ShieldingSourceFitCalc::TraceSourceInfo toTraceSourceInfo() const
  
  
  void fromTraceSourceInfo( const ShieldingSourceFitCalc::TraceSourceInfo &trace )
  {
    //Should check that the model has the nuclide, otherwise dont select a nuclide
    //For development builds should check nuclide is already marked as a trace source in the source fitting model.activityUpdated()
    TraceActivityType type = trace.m_type;
    if( type == TraceActivityType::NumTraceActivityType )
    {
      cerr << "TraceSrcDisplay::deSerialize: Trace activity type in XML ('" << to_str(type) << "'),"
      << " is invalid, setting to total activity." << endl;
      type = TraceActivityType::TotalActivity;
    }
    
    updateAvailableIsotopes();
    updateAvailableActivityTypes();
    
    m_isoSelect->setCurrentIndex(0); //"Select"
    if( trace.m_nuclide )
    {
      bool found = false;
      for( int i = 0; !found && (i < m_isoSelect->count()); ++i )
      {
        const string nuctxt = m_isoSelect->itemText(i).toUTF8();
        if( nuctxt == trace.m_nuclide->symbol )
        {
          found = true;
          m_isoSelect->setCurrentIndex(i);
        }
      }//for( int i = 0; i < m_isoSelect->count(); ++i )
      
      if( !found )
      {
        cerr << "TraceSrcDisplay::deSerialize: Failed to match the expected nuclide ('"
             << trace.m_nuclide->symbol << "') to avaiable ones."
             << endl;
        // Should we throw an error here or something???
      }
    }//if( trace.m_nuclide )
    
    handleUserNuclideChange();
    
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    const string acttxt = PhysicalUnits::printToBestActivityUnits( trace.m_activity, 6, useCi );
    m_activityInput->setText( WString::fromUTF8(acttxt) );
    m_activityInput->setToolTip( WString() );
    
    m_activityType->setCurrentIndex( static_cast<int>(type) );
    
    if( type == TraceActivityType::ExponentialDistribution )
      m_relaxationDistance->setValueText( PhysicalUnits::printToBestLengthUnits(trace.m_relaxationDistance) );
    
    m_allowFitting->setChecked( trace.m_fitActivity );
    
    updateRelaxationDisplay();
    updateTotalActivityFromDisplayActivity();
    handleUserActivityChange();
    handleUserChangeAllowFit();
    
    
#if( PERFORM_DEVELOPER_CHECKS )
    const ShieldingSourceFitCalc::TraceSourceInfo roundtrip = toTraceSourceInfo();
    try
    {
      ShieldingSourceFitCalc::TraceSourceInfo::equalEnough( trace, roundtrip );
    }catch( std::exception &e )
    {
      cerr << "Failed to rount-trip TraceSrcDisplay: " << e.what() << endl;
      assert( 0 );
    }
#endif
  }// void fromTraceSourceInfo( const ShieldingSourceFitCalc::TraceSourceInfo &trace )
  
  
  bool allowFittingActivity() const
  {
    const bool allow = m_allowFitting->isChecked();
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
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      m_isoSelect->removeItem( m_isoSelect->currentIndex() );
      m_isoSelect->setCurrentIndex( 0 );
      m_currentNuclide = nullptr;
      m_currentTotalActivity = m_currentDisplayActivity = 0.0;
      m_currentTotalActivityUncert = -1.0;
      m_activityInput->setText( (useCi ? "0 uCi" : "0 bq") );
      m_activityInput->setToolTip( WString() );
      
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
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
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
    
    m_currentTotalActivityUncert = -1.0; // TODO: we could parse out the +-...
    
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
        
        case TraceActivityType::ExponentialDistribution:
          m_currentTotalActivity = m_currentDisplayActivity * m_parent->inSituSurfaceArea();
          break;
          
        case TraceActivityType::NumTraceActivityType:
          assert( 0 );
          m_activityInput->setText( useCi ? "0 uCi" : "0 bq" );
          m_activityInput->setToolTip( WString() );
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
        
        case TraceActivityType::ExponentialDistribution:
        {
          const double surface_area_m2 = m_parent->inSituSurfaceArea() / PhysicalUnits::m2;
          
          if( surface_area_m2 <= FLT_EPSILON )
          {
            m_currentTotalActivity = m_currentDisplayActivity = 0.0;
            txt = (useCi ? "0 uCi" : "0 bq");
          }else
          {
            txt = PhysicalUnits::printToBestActivityUnits( m_currentTotalActivity/surface_area_m2, 3, useCi );
          }
          break;
        }//case TraceActivityType::ExponentialDistribution:
          
        case TraceActivityType::NumTraceActivityType:
          assert( 0 );
          txt = (useCi ? "0 uCi" : "0 bq");
          break;
      }//switch( type )
      
      m_activityInput->setText( txt );
      m_activityInput->setToolTip( WString() );
    }// try / catch
    
    if( m_currentNuclide )
      m_activityUpdated.emit( m_currentNuclide, m_currentTotalActivity );
  }//void handleUserActivityChange()
  
  
  void handleUserRelaxDistChange()
  {
    double distance = 0.0;
    try
    {
      distance = PhysicalUnits::stringToDistance( m_relaxationDistance->text().toUTF8() );
      if( distance < 10.0*PhysicalUnits::um )
        throw runtime_error( "Relaxation distance must be at least 10 um." );
    }catch( std::exception &e )
    {
      distance = 1.0 * PhysicalUnits::cm;
      m_relaxationDistance->setText( "1.0 cm" );
      passMessage( e.what(), WarningWidget::WarningMsgHigh );
    }// try / catch
    
    // TODO: put in actual recommended dimensions for FOV and depth, based on current detector distance and energy peaks.
    string geom_msg;
    switch( m_parent->geometry() )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        geom_msg = "The in-situ exponential distribution is not recommended for spherical geometry, but may be used.";
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
        geom_msg = "Please make sure cylinder radius is large enough to cover field of view, and depth is at least a number of relaxation lengths large.";
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        geom_msg = "The in-situ exponential distribution is not recommended for cylindrical side-on geometry, but may be used.";
        break;
        
      case GammaInteractionCalc::GeometryType::Rectangular:
        geom_msg = "Please make sure width and height values are large enough to cover field of view, and depth is at least a number of relaxation lengths.";
        break;
        
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_parent->geometry() )
    
    string msg = "This exponentially distributed in-situ surface contamination is usually used "
                 "for soil contamination."
                 " The relaxation distance is the depth from the surface at which ~63% of the"
                 " contamination is above, with the contamination decaying exponentially with"
                 " depth."
                 "<p>" + geom_msg + "</p>";
    m_relaxationDescription->setText( msg );
    
    // We need to update chi2 chart and stuff, so lets cause that to happen.
    m_parent->materialModified().emit(m_parent);
  }//void handleUserRelaxDistChange()
  
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
  
  double relaxationLength() const
  {
    const int currentIndex = m_activityType->currentIndex();
    if( (currentIndex < 0)
       || (currentIndex >= static_cast<int>(TraceActivityType::NumTraceActivityType)) )
      return -1.0;
    
    const TraceActivityType type = static_cast<TraceActivityType>(currentIndex);
    if( type != TraceActivityType::ExponentialDistribution )
      throw runtime_error( "TraceSrcDisplay::relaxationLength() not a ExponentialDistribution" );
    
    return PhysicalUnits::stringToDistance( m_relaxationDistance->text().toUTF8() );
  }//double relaxationLength()
  
  
  void updateAvailableActivityTypes()
  {
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    const int previous = m_activityType->currentIndex();
    assert( previous < static_cast<int>(TraceActivityType::NumTraceActivityType) );
    
    m_activityType->clear();
    
    if( m_parent->isGenericMaterial() )
      return; //Probably wont ever get here
    
    std::shared_ptr<const Material> material = m_parent->material();
    
    // If less dense than 1% of nominal air, then dont allow activity per gram
    const double minDensity = 0.000013*PhysicalUnits::g/PhysicalUnits::cm3;
    const bool allowActPerGram = (material ? (material->density > minDensity) : false);
    
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
        
        case TraceActivityType::ExponentialDistribution:
          m_activityType->addItem( "per m^2 exp" );
          break;
          
        case TraceActivityType::ActivityPerGram:
          if( allowActPerGram )
            m_activityType->addItem( "per gram" );
          break;
          
        case TraceActivityType::NumTraceActivityType:
          break;
      }//switch( type )
    }//for( loop over TraceActivityTypes )
    
    int currentIndex = 0;
    if( (previous >= 0)
        && (previous < static_cast<int>(TraceActivityType::NumTraceActivityType) )
        && (allowActPerGram
            || (previous != static_cast<int>(TraceActivityType::ActivityPerGram))) )
    {
      currentIndex = previous;
    }
    
    m_activityType->setCurrentIndex( currentIndex );
  
    if( previous != currentIndex )
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
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    const int currentIndex = m_activityType->currentIndex();
    
    const double shieldVolume = m_parent->shieldingVolume();
    const double shieldMass = m_parent->shieldingMass();
    
    const double shieldVolumeCm3 = shieldVolume / PhysicalUnits::cm3;
    const double shieldMassGram = shieldMass / PhysicalUnits::gram;
  
    double displayUncert = 0.0;
    m_currentDisplayActivity = 0.0;
    
    switch( TraceActivityType(currentIndex) )
    {
      case TraceActivityType::TotalActivity:
        m_currentDisplayActivity = m_currentTotalActivity;
        if( m_currentTotalActivityUncert > 0.0 )
          displayUncert = m_currentTotalActivityUncert;
        break;
        
      case TraceActivityType::ActivityPerCm3:
        if( shieldVolumeCm3 <= FLT_EPSILON )
        {
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        }else
        {
          m_currentDisplayActivity = m_currentTotalActivity / shieldVolumeCm3;
          if( m_currentTotalActivityUncert > 0.0 )
            displayUncert = m_currentTotalActivityUncert / shieldVolumeCm3;
        }
        break;
        
      case TraceActivityType::ExponentialDistribution:
      {
        const double surface_area_m2 = m_parent->inSituSurfaceArea() / PhysicalUnits::m2;
        if( surface_area_m2 <= FLT_EPSILON )
        {
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        }else
        {
          m_currentDisplayActivity = m_currentTotalActivity / surface_area_m2;
          if( m_currentTotalActivityUncert > 0.0 )
            displayUncert = m_currentTotalActivityUncert / surface_area_m2;
        }
        break;
      }//case TraceActivityType::ActivityPerGram:
        
      case TraceActivityType::ActivityPerGram:
        if( shieldMassGram <= FLT_EPSILON )
        {
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        }else
        {
          m_currentDisplayActivity = m_currentTotalActivity / shieldMassGram;
          if( m_currentTotalActivityUncert > 0.0 )
            displayUncert = m_currentTotalActivityUncert / shieldMassGram;
        }
        break;
      
      case TraceActivityType::NumTraceActivityType:
        assert( 0 );
        m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        break;
    }//switch( type )
    
    const string actTxt = PhysicalUnits::printToBestActivityUnits( m_currentDisplayActivity, 4, useCi );
    m_activityInput->setText( actTxt );
    
    WString tt;
    if( displayUncert > 0.0 )
    {
      const PhysicalUnits::UnitNameValuePair &units
                          = PhysicalUnits::bestActivityUnitHtml( m_currentDisplayActivity, useCi );
      const double value = m_currentDisplayActivity / units.second;
      const double uncert = displayUncert / units.second;
      string txt = PhysicalUnits::printValueWithUncertainty(value, uncert, 6) + " " + units.first;
      tt = WString::fromUTF8(txt);
    }//if( displayUncert > 0.0 )
    
    m_activityInput->setToolTip( tt, TextFormat::XHTMLText );
  }//void updateDispActivityFromTotalActivity()
  
  
  void updateTotalActivityFromDisplayActivity()
  {
    const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    
    
    if( !m_currentNuclide )
    {
      m_currentTotalActivity = m_currentDisplayActivity = 0.0;
      m_currentTotalActivityUncert = -1.0;
      m_activityInput->setText( (useCi ? "0 uCi" : "0 bq") );
      m_activityInput->setToolTip( WString() );
      return;
    }//if( !m_currentNuclide )
    
    
    const int currentIndex = m_activityType->currentIndex();
    
    const double shieldVolume = m_parent->shieldingVolume();
    const double shieldMass = m_parent->shieldingMass();
    
    const double shieldVolumeCm3 = shieldVolume / PhysicalUnits::cm3;
    const double shieldMassGram = shieldMass / PhysicalUnits::gram;
  
    const double origTotalAct = m_currentTotalActivity;
    const double origActUncert = m_currentTotalActivityUncert;
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
        
      case TraceActivityType::ExponentialDistribution:
      {
        const double surface_area_m2 = m_parent->inSituSurfaceArea() / PhysicalUnits::m2;
        if( surface_area_m2 <= FLT_EPSILON )
          m_currentDisplayActivity = m_currentTotalActivity = 0.0;
        else
          m_currentTotalActivity = m_currentDisplayActivity * surface_area_m2;
        
        break;
      }//case TraceActivityType::ExponentialDistribution:
        
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
    
    if( origActUncert )
      m_currentTotalActivityUncert = origActUncert * (m_currentTotalActivity / origTotalAct);
  }//void updateTotalActivityFromDisplayActivity()
  
  
  void deSelectNuclideNoEmit()
  {
    // Note: this function does not emit that its nuclide changed
    assert( (!m_currentNuclide) == (!m_isoSelect->currentIndex()) );
    
    m_isoSelect->setCurrentIndex( 0 );
    if( m_currentNuclide )
    {
      m_currentTotalActivity = m_currentDisplayActivity = 0.0;
      m_currentTotalActivityUncert = -1.0;
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      m_activityInput->setText( (useCi ? "0 uCi" : "0 bq") );
      m_activityInput->setToolTip( WString() );
      
      //const SandiaDecay::Nuclide * const oldNuc = m_currentNuclide;
      m_currentNuclide = nullptr;
      //m_nucChangedSignal.emit( this, oldNuc );
    }
  }//deSelectNuclideNoEmit()
  
  
  void updateForMaterialChange()
  {
    shared_ptr<const Material> mat = m_parent->material();
    if( !mat )
    {
      const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      
      m_currentDisplayActivity = m_currentTotalActivity = 0.0;
      m_currentTotalActivityUncert = -1.0;
      m_activityInput->setText( (useBq ? "0 Bq" : "0 uCi") );
      m_activityInput->setToolTip( WString() );
      m_isoSelect->setCurrentIndex( 0 );
      m_allowFitting->setChecked( false );
      
      if( m_currentNuclide )
      {
        const SandiaDecay::Nuclide * const oldNuc = m_currentNuclide;
        m_currentNuclide = nullptr;
        m_nucChangedSignal.emit( this, oldNuc );
      }
    }//if( !mat )
    
    if( m_isoSelect->isDisabled() != (!mat) )
    {
      // This doesnt actually seem to work, but assuming bug in Wt 3.3.4, so leaving in
      setDisabled( !mat );
      m_isoSelect->setDisabled( !mat );
      m_activityInput->setDisabled( !mat );
      m_allowFitting->setDisabled( !mat );
    }//if( m_isoSelect->isDisabled() != (!mat) )
    
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
  
  void emitUserChaged()
  {
    m_userChanged.emit();
  }
  
  Wt::Signal<> &userChanged()
  {
    return m_userChanged;
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
  
  void setTraceSourceTotalActivity( const double total_activity,
                                   const double total_activity_uncert,
                                   const bool emit_change )
  {
    assert( m_currentNuclide );
    if( !m_currentNuclide )
      throw runtime_error( "setTraceSourceTotalActivity: no current nuclide" );
    
    m_currentTotalActivity = total_activity;
    m_currentTotalActivityUncert = total_activity_uncert;
    
    updateDispActivityFromTotalActivity();
    
    if( emit_change )
      m_activityUpdated.emit( m_currentNuclide, m_currentTotalActivity );
  }//void setTraceSourceTotalActivity(...)
  
  
  void handleUserNuclideChange()
  {
    UndoRedoManager::BlockUndoRedoInserts block;
    
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
      const bool fitAct = srcmodel->fitActivity(nucnum);
      const TraceActivityType type = TraceActivityType( m_activityType->currentIndex() );
      
      switch( type )
      {
        case GammaInteractionCalc::TraceActivityType::TotalActivity:
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerCm3:
          act /= (m_parent->shieldingVolume() / PhysicalUnits::cm3);
          break;
          
        case TraceActivityType::ExponentialDistribution:
          act /= (m_parent->inSituSurfaceArea() / PhysicalUnits::m2);
          break;
          
        case GammaInteractionCalc::TraceActivityType::ActivityPerGram:
          act /= (m_parent->shieldingMass() / PhysicalUnits::gram);
          break;
          
        case GammaInteractionCalc::TraceActivityType::NumTraceActivityType:
          assert( 0 );
          break;
      }//switch( type )
      
      m_currentDisplayActivity = act;
      
      const bool useCi = !UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
      const string actstr = PhysicalUnits::printToBestActivityUnits(act,4,useCi);
      m_activityInput->setValueText( WString::fromUTF8(actstr) );
      
      m_allowFitting->setChecked( fitAct );
      
      updateTotalActivityFromDisplayActivity();
    }//if( srcmodel )
    
    handleUserChangeAllowFit();
    
    if( m_currentNuclide != oldNuc )
      m_nucChangedSignal.emit( this, oldNuc );
  }//void handleUserNuclideChange()
  
  
  void updateRelaxationDisplay()
  {
    const bool show = (m_activityType->currentIndex()
                       == static_cast<int>(TraceActivityType::ExponentialDistribution));
    const bool isChanging = (show == m_relaxationDiv->isHidden());
    
    m_relaxationDiv->setHidden( !show );
    if( show )
      handleUserRelaxDistChange(); // This emits the material modified signal, which will wipe out uncertainties, which is appropriate.
    else if( isChanging )
      m_parent->materialModified().emit(m_parent); //Also emit that material changed when going from exp dist to other one
  }//void showOrHideRelaxation()
  
  
  void handleUserChangeActivityType()
  {
    updateRelaxationDisplay();
    
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
  
  
  void selectNuclideForUser( const SandiaDecay::Nuclide *nuc )
  {
    // Emits signals and such, as if user made the change
    if( !nuc )
    {
      m_isoSelect->setCurrentIndex(0);
      handleUserNuclideChange();
      return;
    }
    
    int nucSelectIndex = -1;
    for( int i = 0; i < m_isoSelect->count(); ++i )
    {
      const string nuctxt = m_isoSelect->itemText(i).toUTF8();
      if( nuctxt == nuc->symbol )
      {
        nucSelectIndex = i;
        break;
      }
    }//for( int i = 0; i < m_isoSelect->count(); ++i )
    
    assert( nucSelectIndex >= 0 );
    
    if( nucSelectIndex >= 0 )
    {
      m_isoSelect->setCurrentIndex(nucSelectIndex);
      handleUserNuclideChange();
    }else
    {
      // Shouldnt ever get here.
      cerr << "selectNuclideForUser: Failed to find '" << nuc->symbol
           << "' as a choice in m_isoSelect." << endl;
    }
  }//void selectNuclideForUser( const SandiaDecay::Nuclide *nuc )


  
  
  
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
        if( model->sourceType(index) == ShieldingSourceFitCalc::ModelSourceType::Intrinsic )
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
    m_label( nullptr ),
    m_massFraction( NULL ),
    m_fitFraction( nullptr ),
    m_nuclide( nuclide )
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSelect.css" );
  
  addStyleClass( "SourceCheckbox" );

  WString txt = nuclide ? WString::fromUTF8(nuclide->symbol) : WString();
  m_useAsSourceCb = new WCheckBox( txt, this );
  m_useAsSourceCb->addStyleClass( "UseAsSrcCb CbNoLineBreak" );
  m_label = new WLabel( "--", this );
  m_label->addStyleClass( "SelfAttenSrcLabel" );
  
  m_massFraction = new NativeFloatSpinBox( this );
  m_label->setBuddy( m_massFraction );
  m_massFraction->setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  m_massFraction->setAttributeValue( "autocorrect", "off" );
  m_massFraction->setAttributeValue( "spellcheck", "off" );
#endif
  m_massFraction->setRange( 0.0, 1.0 );
  m_massFraction->setWidth( 80 );
  m_massFraction->setSpinnerHidden( true );
  m_massFraction->setValue( massFrac );

  m_fitFraction = new WCheckBox( WString::tr("Fit"), this );
  m_fitFraction->addStyleClass( "FitFractionCb CbNoLineBreak" );
  if( !nuclide )
  {
    m_useAsSourceCb->setUnChecked();
    m_useAsSourceCb->disable();
    m_useAsSourceCb->hide();
    m_useAsSourceCb->setHiddenKeepsGeometry( true );
    m_massFraction->disable();
    
    m_label->setText( WString::tr("ss-non-src-frac") );
  }//if( !nuclide )
  
  m_useAsSourceCb->checked().connect( this, &SourceCheckbox::handleUseCbChange );
  m_useAsSourceCb->unChecked().connect( this, &SourceCheckbox::handleUseCbChange );
  
  m_fitFraction->checked().connect( this, &SourceCheckbox::handleFitMassFractionChanged );
  m_fitFraction->unChecked().connect( this, &SourceCheckbox::handleFitMassFractionChanged );
  
  handleUseCbChange();
}//SourceCheckbox constructor

SourceCheckbox::~SourceCheckbox()
{
}


void SourceCheckbox::setLabelText( const Wt::WString &label )
{
  m_label->setText( label );
}


void SourceCheckbox::handleUseCbChange()
{
  if( !m_nuclide )
    return;
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  if( m_useAsSourceCb->isChecked() )
  {
    m_massFraction->show();
    m_fitFraction->show();
    
    WString labeltxt = WString::tr("ss-src-mass-frac-label");
    if( m_nuclide )
      labeltxt.arg(db->element(m_nuclide->atomicNumber)->symbol);
    m_label->setText( labeltxt );
  }else
  {
    m_massFraction->hide();
    m_label->setText( "" );
    m_fitFraction->hide();
  }
}//void handleUseCbChange()

double SourceCheckbox::massFraction() const
{
  return m_massFraction->value();
}

void SourceCheckbox::setMassFraction( double frac, double uncert )
{
  m_massFraction->setValue( std::max( 0.0, std::min( 1.0, frac ) ) );
  
  if( uncert > 0.0 )
  {
    const string txt = PhysicalUnits::printValueWithUncertainty( frac, uncert, 6 );
    setToolTip( WString::tr("ss-tt-with-uncert").arg(txt) );
  }else
  {
    setToolTip( WString::tr("ss-tt-without") );
  }
}

bool SourceCheckbox::useAsSource() const
{
  return m_useAsSourceCb->isChecked();
}

void SourceCheckbox::setUseAsSource( bool use )
{
  m_useAsSourceCb->setChecked( use );
  handleUseCbChange();
}

const SandiaDecay::Nuclide *SourceCheckbox::isotope() const
{
  return m_nuclide;
}

Wt::EventSignal<> &SourceCheckbox::checked()
{
  return m_useAsSourceCb->checked();
}

//Wt::EventSignal<> &SourceCheckbox::changed()
//{
//  return m_useAsSourceCb->changed();
//}

Wt::EventSignal<> &SourceCheckbox::unChecked()
{
  return m_useAsSourceCb->unChecked();
}

Wt::Signal<float> &SourceCheckbox::massFractionChanged()
{
  return m_massFraction->valueChanged();
}

bool SourceCheckbox::fitMassFraction() const
{
  return m_fitFraction->isChecked();
};

void SourceCheckbox::setFitMassFraction( const bool fit )
{
  m_fitFraction->setChecked( fit );
  handleFitMassFractionChanged();
}


Wt::EventSignal<> &SourceCheckbox::fitMassFractionChecked()
{
  return m_fitFraction->checked();
}

Wt::EventSignal<> &SourceCheckbox::fitMassFractionUnChecked()
{
  return m_fitFraction->unChecked();
}

void SourceCheckbox::handleFitMassFractionChanged()
{
  // Nothing to do here?
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
  m_addTraceSrcBtn( nullptr ),
  m_fixedGeometry( false ),
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
  m_asSourceCBs( nullptr ),
  m_traceSources( nullptr )
{
  init();
}


ShieldingSelect::ShieldingSelect( MaterialDB *materialDB,
                                  SourceFitModel *sourceModel,
                                  Wt::WSuggestionPopup *materialSuggest,
                                  const ShieldingSourceDisplay *shieldSource,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_toggleImage( nullptr ),
    m_shieldSrcDisp( shieldSource ),
    m_forFitting( (shieldSource != nullptr) ),
    m_materialDB( materialDB ),
    m_sourceModel( sourceModel ),
    m_geometry( GeometryType::Spherical ),
    m_materialSuggest( materialSuggest ),
    m_materialEdit( nullptr ),
    m_isGenericMaterial( false ),
    m_materialSummary( nullptr ),
    m_closeIcon( nullptr ),
    m_addIcon( nullptr ),
    m_addTraceSrcBtn( nullptr ),
    m_fixedGeometry( false ),
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
    m_asSourceCBs( nullptr ),
    m_traceSources( nullptr )
{
  init();
}


void ShieldingSelect::setClosableAndAddable( bool closeable, WGridLayout *layout )
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
        
    PopupDivMenu *popup = new PopupDivMenu( m_addIcon, PopupDivMenu::TransientMenu );
    PopupDivMenuItem *item = popup->addMenuItem( WString::tr("ss-add-shield-before") );
    item->triggered().connect( this, &ShieldingSelect::emitAddBeforeSignal );
    item = popup->addMenuItem( WString::tr("ss-add-shield-after") );
    item->triggered().connect( this, &ShieldingSelect::emitAddAfterSignal );
    
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


bool ShieldingSelect::fitForAnyMassFractions() const
{
  if( isGenericMaterial() )
    return false;
  
  int nchecked = 0;
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    for( WWidget *widget : etnm.second->children() )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      nchecked += (src && (src->fitMassFraction()));
    }//for( WWidget *widget : isotopeDiv->children() )
  }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )

  return (nchecked > 1);
}//bool fitForAnyMassFractions() const


void ShieldingSelect::setMassFractions( map<const SandiaDecay::Element *,vector<ShieldingSelect::MassFracInfo>> fractions )
{
  if( !m_currentMaterial )
  {
    if( !fractions.empty() )
      throw runtime_error( "setMassFractions: no current material" );
    return;
  }//if( !m_currentMaterial )
    
  // Make sure each input element has at least some nuclides, the nuclides are grouped with correct
  //  element, and nuclides arent repeated.
  for( const auto &el_nucs : fractions )
  {
    const SandiaDecay::Element * const el = el_nucs.first;
    if( !el )
      throw runtime_error( "ShieldingSelect::setMassFractions(): nullptr Element" );
  
    if( el_nucs.second.empty() )
      throw runtime_error( "ShieldingSelect::setMassFractions(): empty Element" );
  
    bool material_has_el = false;
    for( const auto &el_frac : m_currentMaterial->elements )
      material_has_el |= (el_frac.first == el);
    for( const auto &nuc_frac : m_currentMaterial->nuclides )
      material_has_el |= (nuc_frac.first->atomicNumber == el->atomicNumber);
    
    if( !material_has_el )
      throw runtime_error( "ShieldingSelect::setMassFractions(): input has Element not in material" );
    
    set<const SandiaDecay::Nuclide *> seen_nucs;
    for( const auto &nuc_info : el_nucs.second )
    {
      const SandiaDecay::Nuclide * const nuc = nuc_info.m_nuclide;
      if( seen_nucs.count(nuc) )
        throw runtime_error( "ShieldingSelect::setMassFractions(): repeated nuclide." );
      seen_nucs.insert( nuc );
      if( nuc && (nuc->atomicNumber != el->atomicNumber) )
        throw runtime_error( "ShieldingSelect::setMassFractions(): nuclide grouped into wrong element." );
    }//for( const auto &nuc_info : el_nucs.second )
  }//for( const auto &el_nucs : fractions )
  
  //We'll do a bit of a sanity check, and make sure the input sets all the nuclides that the GUI
  //  currently has listed; we'll ignore "other" non-source component.
  set<const SandiaDecay::Nuclide *> used_input_nuclides;
  
  // If a sum of nuclides for an element is more than one, normalize it down to 1.
  for( auto &el_nucs : fractions )
  {
    double frac_sum = 0.0;
    for( const MassFracInfo &nuc : el_nucs.second )
    {
      frac_sum += nuc.m_fraction;
      if( nuc.m_nuclide )
        used_input_nuclides.insert( nuc.m_nuclide );
    }
    
    if( frac_sum > 1.0 )
    {
      for( auto &nuc : el_nucs.second )
      {
        nuc.m_fraction /= frac_sum;
        if( nuc.m_frac_uncert > 0.0 )
          nuc.m_frac_uncert /= frac_sum;
      }
    }//if( s.second > 1.0 )
  }//for( const auto &el_nucs : fractions )
  
  
  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const SandiaDecay::Element * const element = elDiv.first;
    assert( element );
    if( !element )
      continue;
    
    const vector<ShieldingSelect::MassFracInfo> *nuc_infos = nullptr;
    
    const auto frac_pos = fractions.find(element);
    //assert( frac_pos != end(fractions) );
    if( frac_pos != end(fractions) )
       nuc_infos = &(frac_pos->second);

    const vector<WWidget *> &children = elDiv.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb )
        continue;
      
      const SandiaDecay::Nuclide * const nuc = cb->isotope();
      used_input_nuclides.erase( nuc );

      double src_frac_sum = 0.0; 
      const ShieldingSelect::MassFracInfo *nuc_info = nullptr;
      if( nuc_infos )
      {
        for( const ShieldingSelect::MassFracInfo &i : *nuc_infos )
        {
          if( i.m_nuclide == nuc )
            nuc_info = &i;
          
          if( i.m_nuclide && i.m_use_as_source )
            src_frac_sum += i.m_fraction;
        }//for( const auto &i : *nuc_infos )
      }//if( nuc_infos )
      
      
      const bool source_nuc = (nuc_info && nuc_info->m_use_as_source);
      const bool fit_frac = (nuc_info && nuc_info->m_fit_mass_frac);
      
      cb->setUseAsSource( source_nuc );
      if( source_nuc )
      {
        cb->setMassFraction( nuc_info->m_fraction, nuc_info->m_frac_uncert );
      }else if( !nuc_info )
      {
        if( !nuc )
          cb->setMassFraction( std::max(1.0 - src_frac_sum, 0.0), 0.0 );
        cb->setFitMassFraction( false );
        cb->setUseAsSource( false );
      }
    }//for( WWidget *child : children )
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  
  checkAndUpdateMassFractionCanFit();
  updateSelfAttenOtherNucFraction();
  
  assert( used_input_nuclides.empty() );
  if( !used_input_nuclides.empty() )
    throw runtime_error( "ShieldingSelect::setMassFractions: didnt set all the input mass-fractions." );
}//void setMassFractions( map<Element,vector<tuple<Nuclide,double,double,bool>>> )


void ShieldingSelect::setMassFraction( const SandiaDecay::Nuclide * const nuc,
                                       const double fraction,
                                      const double uncert )
{
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    const vector<WWidget *> kids = etnm.second->children();
    
    for( WWidget *widget : kids )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      if( src && (nuc == src->isotope()) )
      {
        src->setMassFraction( fraction, uncert );
        updateSelfAttenOtherNucFraction();
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
  if( ad_gcm2 < 0.0 || ad_gcm2 > GammaInteractionCalc::sm_max_areal_density_g_cm2 )
    throw runtime_error( "setAtomicNumberAndArealDensity: Areal density must be between 0 and "
                         + std::to_string(GammaInteractionCalc::sm_max_areal_density_g_cm2)+ " g/cm2." );
  
  if( !m_isGenericMaterial )
    handleToggleGeneric();
  
  m_atomicNumberEdit->setValue( an );
  m_arealDensityEdit->setValue( ad_gcm2 );
  
  handleMaterialChange();
}//void setAtomicNumberAndArealDensity( const double an, const double ad )


void ShieldingSelect::setAtomicNumberAndArealDensity( const std::string &an, const std::string &ad )
{
  double an_value = 1, ad_value = 0;
  if( !an.empty() && !(stringstream(an) >> an_value) )
    throw runtime_error( "setAtomicNumberAndArealDensity: Atomic number not valid number." );
  
  if( !ad.empty() && !(stringstream(ad) >> ad_value) )
    throw runtime_error( "setAtomicNumberAndArealDensity: Areal density not valid number." );
  
  if( (an_value < 1.0) || (an_value > 100.0) )
    throw runtime_error( "setAtomicNumberAndArealDensity: Atomic number value must be between 1 and 100." );
  
  if( ad_value < 0.0 || ad_value > GammaInteractionCalc::sm_max_areal_density_g_cm2 )
    throw runtime_error( "setAtomicNumberAndArealDensity: Areal density value must be between 0 and "
                        + std::to_string(GammaInteractionCalc::sm_max_areal_density_g_cm2)+ " g/cm2." );
  
  if( !m_isGenericMaterial )
    handleToggleGeneric();
  
  m_atomicNumberEdit->setText( an );
  m_arealDensityEdit->setText( ad );
  
  handleMaterialChange();
}//setAtomicNumberAndArealDensity(...)


void ShieldingSelect::setToNoShielding()
{
  bool wasEmpty = true;
  if( m_isGenericMaterial )
  {
    wasEmpty = (m_atomicNumberEdit->text().empty() || m_arealDensityEdit->text().empty() || (m_arealDensityEdit->value() == 0.0));
    
    m_atomicNumberEdit->setText( "" );
    m_arealDensityEdit->setText( "" );
  }else
  {
    wasEmpty = m_materialEdit->text().empty();
    m_materialEdit->setText( "" );
  }
  
  if( !wasEmpty )
    handleMaterialChange();
}//void setToNoShielding()


bool ShieldingSelect::isNoShielding()
{
  if( m_isGenericMaterial )
  {
    if( m_arealDensityEdit->text().empty() || m_atomicNumberEdit->text().empty() )
      return true;
    if( arealDensity() <= 0.0 )
      return true;
  
    return false;
  }//if( m_isGenericMaterial )
  
  if( m_materialEdit->text().empty() )
    return true;
  
  if( !m_currentMaterial )
    return true;
  
  switch( m_geometry )
  {
    case GammaInteractionCalc::GeometryType::Spherical:
      if( m_thicknessEdit->text().empty() || (thickness() <= 0.0) )
        return true;
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderEndOn:
      if( m_cylLengthEdit->text().empty() || (cylindricalLengthThickness() <= 0.0) )
        return true;
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderSideOn:
      if( m_cylRadiusEdit->text().empty() || (cylindricalRadiusThickness() <= 0.0) )
        return true;
      break;
      
    case GammaInteractionCalc::GeometryType::Rectangular:
      if( m_rectDepthEdit->text().empty() || (rectangularDepthThickness() <= 0.0) )
        return true;
      break;
      
    case GammaInteractionCalc::GeometryType::NumGeometryType:
      break;
  }//switch( m_geometry )
  
  return false;
}//bool isNoShielding();


WLineEdit *ShieldingSelect::materialEdit()
{
  return m_materialEdit;
}


const WLineEdit *ShieldingSelect::materialEdit() const
{
  return m_materialEdit;
}


Wt::WLineEdit *ShieldingSelect::thicknessEdit()
{
  return m_thicknessEdit;
}


const Wt::WLineEdit *ShieldingSelect::thicknessEdit() const
{
  return m_thicknessEdit;
}


vector<WLineEdit *> ShieldingSelect::distanceEdits()
{
  vector<WLineEdit *> answer;
  
  switch( m_geometry )
  {
    case GeometryType::Spherical:
      answer.push_back( m_thicknessEdit );
      break;
      
    case GeometryType::CylinderEndOn:
    case GeometryType::CylinderSideOn:
      answer.push_back( m_cylRadiusEdit );
      answer.push_back( m_cylLengthEdit );
      break;
      
    case GeometryType::Rectangular:
      answer.push_back( m_rectWidthEdit );
      answer.push_back( m_rectHeightEdit );
      answer.push_back( m_rectDepthEdit );
      break;
      
    case GeometryType::NumGeometryType:
      assert(0);
      throw runtime_error("shieldingVolume(): invalid geometry");
      break;
  }//switch( m_geometry )
  
  return answer;
}//vector<WLineEdit *> distanceEdits()


void ShieldingSelect::setSphericalThickness( const double thickness )
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  // The DoseCalc tool ends up calling this function while reacting to a materialModified() or
  //  materialChanged() signal - however, handleMaterialChange() will emit this signal as well,
  //  causing a stack overflow and crash.  So if we arent actually changing the value, we'll break
  //  the change and just return without making a change.
  const string newval = (thickness < 0.0) ?  string("")
                                          : PhysicalUnits::printToBestLengthUnits(thickness,3);
  
  const string oldval = m_thicknessEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_thicknessEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setSphericalThickness( const double thickness )


void ShieldingSelect::setCylindricalRadiusThickness( const double radius )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  // See comments from #setSphericalThickness
  const string newval = (radius < 0.0) ? string("")
                        : PhysicalUnits::printToBestLengthUnits(radius,3);
  
  const string oldval = m_cylRadiusEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_cylRadiusEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setCylindricalRadiusThickness( const double radius )


void ShieldingSelect::setCylindricalLengthThickness( const double length )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  // See comments from #setSphericalThickness
  const string newval = (length < 0.0) ? string("")
                                       : PhysicalUnits::printToBestLengthUnits(length,3);
  
  const string oldval = m_cylLengthEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_cylLengthEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setCylindricalLengthThickness( const double length )


void ShieldingSelect::setRectangularWidthThickness( const double width )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  // See comments from #setSphericalThickness
  const string newval = (width < 0.0) ? string("")
                                      : PhysicalUnits::printToBestLengthUnits(width,3);
  
  const string oldval = m_rectWidthEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_rectWidthEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setRectangularWidth( const double width )


void ShieldingSelect::setRectangularHeightThickness( const double height )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  // See comments from #setSphericalThickness
  const string newval = (height < 0.0) ? string("")
                                       : PhysicalUnits::printToBestLengthUnits(height,3);
  
  const string oldval = m_rectHeightEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_rectHeightEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setRectangularHeight( const double height )


void ShieldingSelect::setRectangularDepthThickness( const double depth )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  // See comments from #setSphericalThickness
  const string newval = (depth < 0.0) ? string("")
                                      : PhysicalUnits::printToBestLengthUnits(depth,3);
  
  const string oldval = m_rectDepthEdit->text().toUTF8();
  if( newval != oldval )
  {
    m_rectDepthEdit->setText( newval );
    //handleMaterialChange();
  }
}//void setRectangularDepth( const double depth )



void ShieldingSelect::setSphericalThicknessEditEnabled( const bool enabled )
{
  checkIsCorrectCurrentGeometry( GeometryType::Spherical, __func__ );
  
  m_thicknessEdit->setEnabled(enabled);
}


NativeFloatSpinBox *ShieldingSelect::arealDensityEdit()
{
  return m_arealDensityEdit;
}


NativeFloatSpinBox *ShieldingSelect::atomicNumberEdit()
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
                                                  const double activity,
                                                  const double uncertainty,
                                                  const bool emit_change )
{
  TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "setTraceSourceTotalActivity: called with invalid nuclide" );
  
  w->setTraceSourceTotalActivity( activity, uncertainty, emit_change );
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


double ShieldingSelect::relaxationLength( const SandiaDecay::Nuclide *nuc ) const
{
  const TraceSrcDisplay *w = traceSourceWidgetForNuclide( nuc );
  if( !w )
    throw runtime_error( "relaxationLength: called with invalid nuclide" );
  
  if( w->activityType() != TraceActivityType::ExponentialDistribution )
    throw runtime_error( "relaxationLength: called for non-exponential distrinbution" );
  
  return w->relaxationLength();
}//double relaxationLength( const SandiaDecay::Nuclide *nuc ) const;


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
  assert( !m_fixedGeometry );
  
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


shared_ptr<const string> ShieldingSelect::getStateAsXml() const
{
  try
  {
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    auto xml_data = make_shared<string>();
    rapidxml::print( std::back_inserter(*xml_data), doc, rapidxml::print_no_indenting );
    
    return xml_data;
  }catch( std::exception &e )
  {
    Wt::log("error") << "Failed to serialize ShieldingSelect during undo/redo step: " << e.what();
  }
  
  return nullptr;
}//ShieldingSelect::getStateAsXml()


void ShieldingSelect::init()
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSelect.css" );
 
  InterSpec *interspec = InterSpec::instance();
  
  //TODO/NOTE: had to hard code this as false because there is no way
  //to easily get the preference via InterSpec because
  //is still initializing when calling at this moment.
  const bool showToolTips = interspec ? UserPreferences::preferenceValue<bool>( "ShowTooltips", interspec ) : false;
  if( interspec )
    interspec->useMessageResourceBundle( "ShieldingSelect" );
  
  addStyleClass( "ShieldingSelect" );
  
//  int voidIndex = -1;
  WContainerWidget *materialDiv = new WContainerWidget( this );
  WGridLayout* materialDivLayout = new WGridLayout();
  materialDivLayout->setContentsMargins(2,2,2,2);
  materialDiv->setLayout(materialDivLayout);
  
  
  
  m_toggleImage = new Wt::WImage(Wt::WLink("InterSpec_resources/images/shield.png"));
  m_toggleImage->clicked().connect( this,&ShieldingSelect::handleToggleGeneric );
  m_toggleImage->decorationStyle().setCursor(PointingHandCursor);
  m_toggleImage->addStyleClass( "Wt-icon" );
 
  m_toggleImage->clicked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  
  HelpSystem::attachToolTipOn( m_toggleImage, WString::tr("ss-tt-shield-type-toggle"),
                              showToolTips, HelpSystem::ToolTipPosition::Top );
  
  materialDivLayout->addWidget( m_toggleImage, 0, 0, AlignLeft );
  
  m_materialEdit = new WLineEdit( "" );
  
  m_materialEdit->setAutoComplete( false );
  m_materialEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_materialEdit->setAttributeValue( "autocorrect", "off" );
  m_materialEdit->setAttributeValue( "spellcheck", "off" );
#endif

  m_materialEdit->changed().connect( this, &ShieldingSelect::handleUserChangedMaterialName );
  m_materialEdit->enterPressed().connect( this, &ShieldingSelect::handleUserChangedMaterialName );

  // Instead of hooking undo/redo to changed(), we'll have handleUserChangedMaterialName() handle it
  // We will only insert an undo/redo step when the field losses focus.
  //m_materialEdit->changed().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  
  if( m_forFitting )
  {
    materialDivLayout->addWidget( m_materialEdit, 0, 1, AlignmentFlag::AlignMiddle );
  }else
  {
    materialDivLayout->addWidget( m_materialEdit, 0, 1 );
  }
  
  WString material_name_tt;
  if( m_forFitting )
  {
    material_name_tt = WString("<p>{1}</p><p>{2}</p>")
      .arg( WString::tr("ss-tt-material-name") )
      .arg( WString::tr("ss-tt-material-name-fit") );
  }else
  {
    material_name_tt = WString::tr("ss-tt-material-name");
  }
  
  HelpSystem::attachToolTipOn( m_materialEdit, material_name_tt,
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
  //label->setToolTip( WString::tr("ss-tt-areal-density") );
  genericMatLayout->addWidget( label, 0, 2+m_forFitting, AlignMiddle );
  
  m_arealDensityEdit = new NativeFloatSpinBox();
  m_arealDensityEdit->setSpinnerHidden( true );
  m_arealDensityEdit->setFormatString( "%.4g" );
  m_arealDensityEdit->setTextSize( 5 );
  m_arealDensityEdit->setRange( 0.0f,
                            static_cast<float>(GammaInteractionCalc::sm_max_areal_density_g_cm2) );
  label->setBuddy( m_arealDensityEdit );
  //m_arealDensityEdit->setToolTip( WString::tr("ss-tt-areal-density") );
  
  if( m_forFitting )
    m_arealDensityEdit->setValue( 15.0f );
  else
    m_arealDensityEdit->setPlaceholderText( WString::tr("ss-areal-density-empty-txt") );
  
  genericMatLayout->addWidget( m_arealDensityEdit, 0, 3+m_forFitting, AlignMiddle );
  genericMatLayout->setColumnStretch( 3+m_forFitting, 1 );
  
  HelpSystem::attachToolTipOn( {label, m_arealDensityEdit}, WString::tr("ss-tt-areal-density"),
                              showToolTips, HelpSystem::ToolTipPosition::Top );
  
  label = new WLabel( "g/cm<sup>2</sup>");
  label->setAttributeValue( "style", "font-size: 75%;" );
  genericMatLayout->addWidget( label, 0, 4+m_forFitting, AlignMiddle );
  
  if( m_forFitting )
  {
    m_fitArealDensityCB = new WCheckBox( WString::tr("Fit") );
    m_fitArealDensityCB->setChecked( true );
    m_fitArealDensityCB->addStyleClass( "CbNoLineBreak" );
    genericMatLayout->addWidget( m_fitArealDensityCB, 0, 6, AlignMiddle );
  }
  
  label = new WLabel( "AN" );
  //label->setToolTip( WString::tr("ss-tt-atomic-number") );
  genericMatLayout->addWidget( label, 0, 0, AlignMiddle );
  
  m_atomicNumberEdit = new NativeFloatSpinBox();
  m_atomicNumberEdit->setSpinnerHidden( true );
  m_atomicNumberEdit->setFormatString( "%.3g" );
  m_atomicNumberEdit->setTextSize( 5 );
  m_atomicNumberEdit->setRange( MassAttenuation::sm_min_xs_atomic_number,
                                 MassAttenuation::sm_max_xs_atomic_number );
  label->setBuddy( m_atomicNumberEdit );
  //m_atomicNumberEdit->setToolTip( WString::tr("ss-tt-atomic-number") );
  
  if( m_forFitting )
    m_atomicNumberEdit->setValue( 15.0f );
  else
    m_atomicNumberEdit->setPlaceholderText( WString::tr("ss-atomic-number-empty-txt") );
  
  genericMatLayout->addWidget( m_atomicNumberEdit, 0, 1, AlignMiddle );
  genericMatLayout->setColumnStretch( 1, 1 );
  
  HelpSystem::attachToolTipOn( {label, m_atomicNumberEdit}, WString::tr("ss-tt-atomic-number"),
                              showToolTips, HelpSystem::ToolTipPosition::Top );
  
  if( m_forFitting )
  {
    m_fitAtomicNumberCB = new WCheckBox( "Fit" );
    m_fitAtomicNumberCB->setChecked( false );
    m_fitAtomicNumberCB->addStyleClass( "CbNoLineBreak" );
    genericMatLayout->addWidget( m_fitAtomicNumberCB, 0, 2, AlignMiddle );
    
    m_asSourceCBs = new WContainerWidget( this );
    m_asSourceCBs->addStyleClass( "ShieldingAsSourceCBDiv" );
    label = new WLabel( "Source for:", m_asSourceCBs );
    label->setInline( false );
    m_asSourceCBs->hide();
    HelpSystem::attachToolTipOn( m_asSourceCBs, WString::tr("ss-source-for-cb"), showToolTips );
  }//if( m_forFitting )
  
  
  m_atomicNumberEdit->valueChanged().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_arealDensityEdit->valueChanged().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  
  m_atomicNumberEdit->valueChanged().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  m_arealDensityEdit->valueChanged().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  
  // A validator for the various geometry distances
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  
  
  // A lamda for setting up dimension edits for the various geometry dimensions
  auto setupDimEdit = [this,distValidator]( const WString &labelTxt, WLineEdit *&edit, WCheckBox *&fitCb, WGridLayout *grid ){
    const int row = grid->rowCount();
    
    WLabel *label = new WLabel( labelTxt );
    grid->addWidget( label, row, 0, AlignMiddle );
    
    edit = new WLineEdit( "1.0 cm" );
    edit->setValidator( distValidator );
    
    edit->setAutoComplete( false );
    edit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    edit->setAttributeValue( "autocorrect", "off" );
    edit->setAttributeValue( "spellcheck", "off" );
#endif
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
    
    edit->changed().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
    edit->enterPressed().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
    
    if( m_forFitting )
    {
      //edit->setWidth( 150 );
      grid->addWidget( edit, row, 1, AlignMiddle );
      
      fitCb = new WCheckBox( WString::tr("Fit") );
      fitCb->setChecked( false );
      fitCb->addStyleClass( "CbNoLineBreak" );
      grid->addWidget( fitCb, row, 2, AlignMiddle | AlignRight );
      
      fitCb->checked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
      fitCb->unChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
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
  sphericalLayout->setColumnStretch( 1, 1 );
  
  const bool fistShield = (!m_shieldSrcDisp || !m_shieldSrcDisp->numberShieldings());
  
  const char *lbltxt_key = "ss-thickness"; //fistShield ? "Radius" : "Thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_thicknessEdit, m_fitThicknessCB, sphericalLayout );
  
  if( !m_forFitting )
  {
    // TODO: right now the only place a ShieldingSelect will be non-generic or non-spherical is
    //       in the Shielding/Activity fit tool, which we always have (m_forFitting == true), but
    //       if geometry is cylindrical or rectangular, and (m_forFitting == false), then the
    //       material summary wont be visible.
    sphericalLayout->addWidget( m_materialSummary, 0, 2, AlignMiddle );
    sphericalLayout->setColumnStretch( 2, 1 );
  }//if( !m_forFitting )
  
  // Begin setting up cylindrical widgets
  m_cylindricalDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_cylindricalDiv );
  
  WGridLayout *cylindricalLayout = new WGridLayout();
  m_cylindricalDiv->setLayout( cylindricalLayout );
  cylindricalLayout->setContentsMargins( 3, 3, 3, 3 );
  cylindricalLayout->setColumnStretch( 1, 1 );
  
  lbltxt_key = fistShield ? "ss-radius" : "ss-radial-thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_cylRadiusEdit, m_fitCylRadiusCB, cylindricalLayout );
  lbltxt_key = fistShield ? "ss-half-length" : "ss-length-thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_cylLengthEdit, m_fitCylLengthCB, cylindricalLayout );

  
  // Begin setting up rectangular widgets
  m_rectangularDiv = new WContainerWidget();
  m_dimensionsStack->addWidget( m_rectangularDiv );
  
  WGridLayout *rectangularLayout = new WGridLayout();
  m_rectangularDiv->setLayout( rectangularLayout );
  rectangularLayout->setContentsMargins( 3, 3, 3, 3 );
  rectangularLayout->setColumnStretch( 1, 1 );
  
  lbltxt_key = fistShield ? "ss-half-width"  : "ss-width-thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_rectWidthEdit, m_fitRectWidthCB, rectangularLayout );
  lbltxt_key = fistShield ? "ss-half-height" : "ss-height-thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_rectHeightEdit, m_fitRectHeightCB, rectangularLayout );
  lbltxt_key = fistShield ? "ss-half-depth"  : "ss-depth-thickness";
  setupDimEdit( WString::tr(lbltxt_key), m_rectDepthEdit, m_fitRectDepthCB, rectangularLayout );
  
  // Finally we have to
  if( m_forFitting )
  {
    m_addTraceSrcBtn = new WPushButton( WString::tr("ss-add-trace-src"), this );
    m_addTraceSrcBtn->addStyleClass( "LinkBtn AddTrcSrcBtn" );
    m_addTraceSrcBtn->clicked().connect( this, &ShieldingSelect::addTraceSource );
    m_addTraceSrcBtn->clicked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  }//if( m_forFitting )
  
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
    // We will double check that m_materialSuggest is still owned by the InterSpec instance.
    //  If this ShieldingSelect is in a WDialog, it may be destructing after the InterSpec instance
    //  destructs!  In that case, we shouldnt access `m_materialSuggest`, as its already been
    //  deleted.
    InterSpec *interspec = InterSpec::instance();
    if( interspec )
    {
      const vector<Wt::WObject *> &kids = interspec->Wt::WObject::children();
      // kids.size() is usually just 2 or 3, so this isnt that heavy of an operation.
      auto pos = std::find( begin(kids), end(kids), static_cast<WObject *>(m_materialSuggest) );
      if( pos != end(kids) )
      {
        // We get here, for example, when you manually close the Activity/Shielding fit window (or
        //  remove a shielding from within it).
        m_materialSuggest->removeEdit( m_materialEdit );
      }else
      {
        // I dont think we get here - but I'll leave the assert in, just to see
        //  -- we get here sometimes if DoseCalc tool is showing and we clear the session
        cerr << "~ShieldingSelect(): Suggest not in DOM, not removing edit from suggestion" << endl;
//        assert( 0 );
      }
    }else
    {
      // We get here when you close a tab, or the window or whatever (e.g., WApplication is being
      //  destructed); I think because WApplication::instance() is nullptr, so so is `interspec`
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
          removingIsotopeAsSource().emit( cb->isotope(), ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
        }//if( cb->useAsSource() )
      }//for( WWidget *child : children )
    }//for(...)
  }//if( m_asSourceCBs )

  // Remove trace sources so model will get updated
  if( m_traceSources )
  {
    const vector<WWidget *> kids = m_traceSources->children();
    for( WWidget *w : kids )
    {
      TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
      assert( src );
      if( src ) // JIC
        removeTraceSourceWidget( src );
    }
  }//if( m_traceSources )
  
  
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
  UndoRedoManager::BlockUndoRedoInserts block;
  
  // TODO: maybe we shouldnt add a trace source if material() is nullptr; also when material()
  //       becomes nullptr maybe we should remove trace sources
  if( !m_traceSources )
  {
    m_traceSources = new WContainerWidget( this );
    m_traceSources->addStyleClass( "TraceSrcContainer" );
  }
  
  TraceSrcDisplay *src = new TraceSrcDisplay( this );
  src->nucChangedSignal().connect( this, &ShieldingSelect::handleTraceSourceNuclideChange );
  src->activityUpdated().connect( this, &ShieldingSelect::handleTraceSourceActivityChange );
  
  src->userChanged().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  
  // Check how many nuclides arent already volumetric sources, and if exactly one, just select it
  const SandiaDecay::Nuclide *nuc = nullptr;
  const int nnuc = m_sourceModel->numNuclides();
  for( int i = 0; i < nnuc; ++i )
  {
    if( !m_sourceModel->isVolumetricSource(i) )
    {
      if( nuc )
      {
        // We have more than one possible nuclide the user could select; so dont select any of them
        nuc = nullptr;
        break;
      }
      
      nuc = m_sourceModel->nuclide(i);
    }//if( not a volumetric source )
  }//for( int i = 0; i < nnuc; ++i )
  
  if( nuc && material() )
    src->selectNuclideForUser( nuc );
}//void addTraceSource()


void ShieldingSelect::removeTraceSourceWidget( TraceSrcDisplay *toRemove )
{
  if( !toRemove )
    return;
  
  UndoRedoManager::BlockUndoRedoInserts block;
  
  const vector<WWidget *> kids = m_traceSources->children();
  
  for( WWidget *w : kids )
  {
    TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
    assert( src );
    if( !src ) // JIC
      continue;
    
    if( src == toRemove )
    {
      handleTraceSourceWidgetAboutToBeRemoved( src );
      
      delete src;
      
      setTraceSourceBtnStatus();
      
      handleUserChangeForUndoRedo();
      
      return;
    }//if( src == toRemove )
  }//for( WWidget *w : m_traceSources->children() )
  
  
  throw runtime_error( "ShieldingSelect::removeTraceSourceWidget(TraceSrcDisplay *): called with invalid display." );
}//void removeTraceSourceWidget( TraceSrcDisplay *src )



void ShieldingSelect::setTraceSourceBtnStatus()
{
  if( !m_addTraceSrcBtn )
    return;
  
  if( m_isGenericMaterial || !m_sourceModel || m_fixedGeometry )
  {
    m_addTraceSrcBtn->setDisabled( true );
    m_addTraceSrcBtn->setHidden( true );
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
  
  m_addTraceSrcBtn->setDisabled( !numAvailableNuclides );
  m_addTraceSrcBtn->setHidden( !numAvailableNuclides );
}//void setTraceSourceBtnStatus()


void ShieldingSelect::handleTraceSourceNuclideChange( TraceSrcDisplay *changedSrc, const SandiaDecay::Nuclide *oldNuc )
{
  assert( changedSrc );
  
  if( !material() )
    changedSrc->deSelectNuclideNoEmit();
  
  const SandiaDecay::Nuclide *nuc = changedSrc->nuclide();
  
  if( nuc == oldNuc )
    return;
  
  if( oldNuc )
    removingIsotopeAsSource().emit( oldNuc, ShieldingSourceFitCalc::ModelSourceType::Trace );
  
  if( changedSrc && changedSrc->nuclide() )
    addingIsotopeAsSource().emit( changedSrc->nuclide(), ShieldingSourceFitCalc::ModelSourceType::Trace );
  
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
  
  setTraceSourceBtnStatus();
}//void handleTraceSourceNuclideChange( TraceSrcDisplay *src );


void ShieldingSelect::handleTraceSourceActivityChange( const SandiaDecay::Nuclide *nuc, const double activity )
{
  m_activityFromVolumeNeedUpdating.emit( this, nuc );
}

void ShieldingSelect::handleTraceSourceWidgetAboutToBeRemoved( TraceSrcDisplay *src )
{
  if( !src || !src->nuclide() )
    return;
    
  removingIsotopeAsSource().emit( src->nuclide(), ShieldingSourceFitCalc::ModelSourceType::Trace );
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


double ShieldingSelect::inSituSurfaceArea() const
{
  if( isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::shieldingVolume(): You should not call this function"
                        " for generic materials" );
  
  double surface_area = -1.0;
  
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
      
      surface_area = 4 * PhysicalUnits::pi * outer_rad * outer_rad;
      
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
      
      const double rad = inner_rad + cylindricalRadiusThickness();
      const double half_len = inner_half_length + cylindricalLengthThickness();
      
      if( m_geometry == GeometryType::CylinderEndOn )
        surface_area = PhysicalUnits::pi * rad * rad;
      else
        surface_area = PhysicalUnits::pi * 2.0 * rad * 2.0 * half_len;
    
      break;
    }//case CylinderEndOn and CylinderSideOn
      
    case GeometryType::Rectangular:
    {
      double inner_half_width = 0.0, inner_half_height = 0.0;
      if( m_shieldSrcDisp )
      {
        for( const ShieldingSelect *inner = m_shieldSrcDisp->innerShielding(this);
            inner; inner = m_shieldSrcDisp->innerShielding(inner) )
        {
          assert( inner->geometry() == m_geometry );
          inner_half_width  += inner->rectangularWidthThickness();
          inner_half_height += inner->rectangularHeightThickness();
        }
      }//if( m_shieldSrcDisp )
      
      
      const double half_width  = inner_half_width  + rectangularWidthThickness();
      const double half_height = inner_half_height + rectangularHeightThickness();
      
      surface_area = 4.0 * half_width * half_height;
      
      break;
    }//case Rectangular:
      
    case GeometryType::NumGeometryType:
      assert(0);
      throw runtime_error("shieldingVolume(): invalid geometry");
      break;
  }//switch( m_geometry )
  
  assert( surface_area >= 0.0 );
  
  return surface_area;
}//double ShieldingSelect::inSituSurfaceArea() const


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


Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> &ShieldingSelect::addingIsotopeAsSource()
{
  return m_addingIsotopeAsSource;
}


Wt::Signal<const SandiaDecay::Nuclide *,ShieldingSourceFitCalc::ModelSourceType> &ShieldingSelect::removingIsotopeAsSource()
{
  return m_removingIsotopeAsSource;
}


Wt::Signal<ShieldingSelect *, shared_ptr<const string>, shared_ptr<const string>>
           &ShieldingSelect::userChangedStateSignal()
{
  if( !m_prevState )
    m_prevState = getStateAsXml();
    
  return m_userChangedStateSignal;
}


void ShieldingSelect::setGeometry( GammaInteractionCalc::GeometryType type )
{
  assert( !m_fixedGeometry || (type == GammaInteractionCalc::GeometryType::Spherical) );
  if( m_fixedGeometry && (type != GammaInteractionCalc::GeometryType::Spherical) )
    throw std::logic_error( "Geometry type must be spherical when fixed geometry" );
  
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


void ShieldingSelect::setFixedGeometry( const bool fixed_geom )
{
  if( !m_forFitting )
    throw std::logic_error( "ShieldingSelect::setFixedGeometry: should not be called for non-fitting instances." );
  
  //if( fixed_geom == m_fixedGeometry )
  //  return;
  
  m_fixedGeometry = fixed_geom;
  if( m_fixedGeometry )
  {
    if( m_geometry != GammaInteractionCalc::GeometryType::Spherical )
      setGeometry( GammaInteractionCalc::GeometryType::Spherical );
    
    vector<const SandiaDecay::Nuclide *> self_atten_nucs = selfAttenNuclides();
    vector<const SandiaDecay::Nuclide *> trace_srcs = traceSourceNuclides();
    
    if( m_traceSources )
    {
      const vector<WWidget *> kids = m_traceSources->children();
      for( WWidget *w : kids )
      {
        TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
        if( src )
        {
          //src->deSelectNuclideNoEmit();
          removeTraceSourceWidget( src );
        }
      }
    }//if( m_traceSources )
    
    if( m_asSourceCBs )
    {
      for( const ElementToNuclideMap::value_type &etnp : m_sourceIsotopes )
      {
        const vector<WWidget *> kids = etnp.second->children();
        for( WWidget *child : kids )
        {
          SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
          if( cb && cb->useAsSource() && cb->isotope() )
          {
            cb->setUseAsSource( false );
            removingIsotopeAsSource().emit( cb->isotope(), ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
          }//if( cb->useAsSource() )
        }//for( WWidget *child : children )
        
        etnp.second->hide();
      }//for(...)
      
      m_asSourceCBs->hide();
    }//if( m_asSourceCBs )
    
    assert( selfAttenNuclides().empty() );
  }else
  {
    if( m_asSourceCBs )
    {
      if( !m_sourceIsotopes.empty() )
        m_asSourceCBs->show();
      
      for( const ElementToNuclideMap::value_type &etnp : m_sourceIsotopes )
      {
        etnp.second->show();
      }//for(...)
    }//if( m_asSourceCBs )
  }//if( m_fixedGeometry ) / else
  
  setTraceSourceBtnStatus();
}//void setFixedGeometry( const bool fixed_geom );


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
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
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
  
  switch( m_atomicNumberEdit->validate() )
  {
    case Wt::WValidator::Invalid:
      throw runtime_error( "ShieldingSelect::atomicNumber() invalid value entered" );
      
    case Wt::WValidator::InvalidEmpty:
      return 13.0;
      
    case Wt::WValidator::Valid:
      break;
  }//switch( m_arealDensityEdit->validate() )
  
  return m_atomicNumberEdit->value();
}//double atomicNumber() const


double ShieldingSelect::arealDensity() const
{
  if( !isGenericMaterial() )
    throw runtime_error( "ShieldingSelect::arealDensity() can not be called when a material" );
  
  switch( m_arealDensityEdit->validate() )
  {
    case Wt::WValidator::Invalid:
      throw runtime_error( "ShieldingSelect::arealDensity() invalid value entered" );
      
    case Wt::WValidator::InvalidEmpty:
      return 0.0;
      
    case Wt::WValidator::Valid:
      break;
  }//switch( m_arealDensityEdit->validate() )
  
  double answer = m_arealDensityEdit->value();
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
  
  return (m_fitCylRadiusCB->isVisible() && m_fitCylRadiusCB->isChecked());
}//bool fitCylindricalRadiusThickness() const


bool ShieldingSelect::fitCylindricalLengthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderEndOn, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return (m_fitCylLengthCB->isVisible() && m_fitCylLengthCB->isChecked());
}//bool fitCylindricalLengthThickness() const


bool ShieldingSelect::fitRectangularWidthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return (m_fitRectWidthCB->isVisible() && m_fitRectWidthCB->isChecked());
}//bool fitRectangularWidthThickness() const


bool ShieldingSelect::fitRectangularHeightThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return (m_fitRectHeightCB->isVisible() && m_fitRectHeightCB->isChecked());
}//bool fitRectangularHeightThickness() const


bool ShieldingSelect::fitRectangularDepthThickness() const
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  return m_fitRectDepthCB->isChecked();
}//bool fitRectangularDepthThickness() const


void ShieldingSelect::setFitCylindricalRadiusEnabled( const bool allow )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderSideOn, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  WCheckBox *cb = m_fitCylRadiusCB;
  const bool previous = cb->isVisible();
  if( previous == allow )
    return;
  
  cb->setChecked( false );
  cb->setHidden( !allow );
  
  assert( allow || !fitCylindricalRadiusThickness() );
}//void setFitCylindricalRadiusEnabled( const bool allow )


void ShieldingSelect::setFitCylindricalLengthEnabled( const bool allow )
{
  checkIsCorrectCurrentGeometry( GeometryType::CylinderSideOn, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  WCheckBox *cb = m_fitCylLengthCB;
  const bool previous = cb->isVisible();
  if( previous == allow )
    return;
  
  cb->setChecked( false );
  cb->setHidden( !allow );
  
  assert( allow || !fitCylindricalLengthThickness() );
}//void setFitCylindricalLengthEnabled( const bool allow )


void ShieldingSelect::setFitRectangularHeightEnabled( const bool allow )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  WCheckBox *cb = m_fitRectHeightCB;
  const bool previous = cb->isVisible();
  if( previous == allow )
    return;
  
  cb->setChecked( false );
  cb->setHidden( !allow );
  
  assert( allow || !fitRectangularHeightThickness() );
}//void setFitRectangularHeightEnabled( const bool allow )


void ShieldingSelect::setFitRectangularWidthEnabled( const bool allow )
{
  checkIsCorrectCurrentGeometry( GeometryType::Rectangular, __func__ );
  
  if( !m_forFitting )
    throw runtime_error( __func__ + string(": can't be called when not for fitting.") );
  
  WCheckBox *cb = m_fitRectWidthCB;
  const bool previous = cb->isVisible();
  if( previous == allow )
    return;
  
  cb->setChecked( false );
  cb->setHidden( !allow );
  
  assert( allow || !fitRectangularWidthThickness() );
}//void setFitRectangularWidthEnabled( const bool allow )


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
    if( m_materialSuggest )
      m_materialSuggest->addSuggestion( mat->name, mat->name );
    return mat;
  }catch(...)
  {
    // No luck
  }

  return nullptr;
}//std::shared_ptr<Material> material( const std::string &text )


std::shared_ptr<const Material> ShieldingSelect::material()
{
  if( m_isGenericMaterial )
    return nullptr;

  updateMaterialFromUserInputTxt();

  return m_currentMaterial;
}//const Material *material() const;


std::shared_ptr<const Material> ShieldingSelect::currentMaterial() const
{
  return m_currentMaterial;
}

void ShieldingSelect::checkAndUpdateMassFractionCanFit()
{
  for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
  {
    int nchecked = 0, num_fit_frac = 0;
    SourceCheckbox *non_src_cb = nullptr;
    
    for( WWidget *widget : etnm.second->children() )
    {
      SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
      nchecked += (src && src->isotope() && src->useAsSource());
      num_fit_frac += (src && src->isotope() && src->fitMassFraction());
      if( src && !src->isotope() )
        non_src_cb = src;
    }//for( WWidget *widget : isotopeDiv->children() )
    
    assert( non_src_cb );
    if( non_src_cb )
      non_src_cb->setHidden( !nchecked );
    
    if( (num_fit_frac == 0) && non_src_cb )
      non_src_cb->setFitMassFraction( false );
    else if( num_fit_frac == 1 )
      non_src_cb->setFitMassFraction( true );
  }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
}//void checkAndUpdateMassFractionCanFit()


void ShieldingSelect::isotopeCheckedCallback( const SandiaDecay::Nuclide *nuc )
{
  checkAndUpdateMassFractionCanFit();
  
  // Make sure sum of source isotopes are not above 1.0;
  if( nuc )
  {
    SourceCheckbox *this_src_cb = nullptr;
    double other_nuc_sum = 0.0;
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    
    for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
    {
      for( WWidget *child : elDiv.second->children() )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->isotope() && cb->useAsSource() )
        {
          if( cb->isotope() == nuc )
            this_src_cb = cb;
          else
            other_nuc_sum += cb->massFraction();
        }
      }//for( loop over mass fraction widgets for this element )
    }//for( loop over mapping of element to their mass fractions )
    
    assert( this_src_cb );
    if( this_src_cb )
    {
      const double initialMassFrac = this_src_cb->massFraction();
      if( IsInf(initialMassFrac) || IsNan(initialMassFrac) || (initialMassFrac < 0.0) )
        this_src_cb->setMassFraction( 0.0, -1.0 );
      else if( initialMassFrac > 1.0 )
        this_src_cb->setMassFraction( 1.0, -1.0 );
        
      if( other_nuc_sum > 1.0 )
      {
        this_src_cb->setMassFraction( 0.0, -1.0 );
        
        for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
        {
          for( WWidget *child : elDiv.second->children() )
          {
            SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
            if( cb && cb->useAsSource() && (cb->isotope() != nuc) )
              cb->setMassFraction( cb->massFraction() / other_nuc_sum, -1.0 );
          }//for( loop over mass fraction widgets for this element )
        }//for( loop over mapping of element to their mass fractions )
      }else if( (other_nuc_sum + this_src_cb->massFraction()) > 1.0 )
      {
        this_src_cb->setMassFraction( 1.0 - other_nuc_sum, -1.0 );
      }
    }//if( this_src_cb )
    
    updateSelfAttenOtherNucFraction();
  }//if( nuc )
  
  m_addingIsotopeAsSource.emit( nuc, ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
}//void isotopeCheckedCallback( const std::string symbol )



void ShieldingSelect::isotopeUnCheckedCallback( const SandiaDecay::Nuclide *iso )
{
  checkAndUpdateMassFractionCanFit();
  setTraceSourceBtnStatus();
  m_removingIsotopeAsSource.emit( iso, ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
}//void isotopeUnCheckedCallback( const std::string symbol )



void ShieldingSelect::uncheckSourceIsotopeCheckBox( const SandiaDecay::Nuclide *iso )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  
  if( m_traceSources )
  {
    const vector<WWidget *> kids = m_traceSources->children();
    for( WWidget *w : kids )
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
//      removingIsotopeAsSource().emit( symbol, ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
      cb->setUseAsSource( false );
    }//if( cb->isChecked() )
  }//for( WWidget *child : children )
  
  checkAndUpdateMassFractionCanFit();
  setTraceSourceBtnStatus();
}//void uncheckSourceIsotopeCheckBox( const std::string &symol )


void ShieldingSelect::sourceRemovedFromModel( const SandiaDecay::Nuclide *nuc )
{
  
  // Make sure no trace sources are using nuc, and if they are, remove that trace source
  if( m_traceSources )
  {
    vector<TraceSrcDisplay *> todel;
    const vector<WWidget *> kids = m_traceSources->children();
    for( WWidget *w : kids )
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

  checkAndUpdateMassFractionCanFit();
  setTraceSourceBtnStatus();
  updateSelfAttenOtherNucFraction();
}//void sourceRemovedFromModel( const std::string &symbol )


map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> ShieldingSelect::sourceNuclideMassFractions() const
{
  map<const SandiaDecay::Element *,vector<ShieldingSelect::NucMasFrac>> answer;
  
  shared_ptr<const Material> mat = m_currentMaterial;
  
  if( !mat )
    return answer;
  
  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const SandiaDecay::Element * const element = elDiv.first;
    assert( element );
    if( !element )
      continue;
    
    vector< ShieldingSelect::NucMasFrac > el_answer;
    const vector<WWidget *> &children = elDiv.second->children();

    // We dont want to include the "other" non-source component into the results, if we 
    //  have no actial sources selected - so lets check for that first.
    int num_use_as_src = 0;
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      num_use_as_src += (cb && cb->isotope() && cb->useAsSource());
    }

    if( !num_use_as_src )
      continue;

    double frac_sum = 0.0;
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );

      const SandiaDecay::Nuclide *iso = cb ? cb->isotope() : nullptr;
      assert( !iso || (iso->atomicNumber == element->atomicNumber));
      
      if( !cb || (cb->isotope() && !cb->useAsSource()) || (!cb->isotope() && !cb->fitMassFraction()) )  
        continue;

      const double mass_frac = cb->massFraction();
      const bool fit_mass_frac = cb->fitMassFraction();
      frac_sum += mass_frac;

      el_answer.emplace_back( iso, mass_frac, fit_mass_frac );
    }//for( WWidget *child : children )
    
    assert( frac_sum < (1.0 + 1.0E-6) );

    if( frac_sum > 1.0 )
    {
      const double frac_over = frac_sum - 1.0;
      for( auto &v : el_answer )
      {
        if( !get<0>(v) )
        {
          double other_frac = get<1>(v);
          const double to_sub = std::min( other_frac, frac_over );  
          get<1>(v) -= to_sub;
          frac_sum -= to_sub;
          break;
        }
      }

      for( auto &v : el_answer )
      {
        if( get<0>(v) )
          get<1>(v) /= frac_sum;
      }
    }//if( frac_sum > 1.0 )
    
    if( !el_answer.empty() )
      answer[element] = el_answer;
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  
  return answer;
}//std::vector< NucMasFrac > sourceNuclideMassFractions()


vector<const SandiaDecay::Nuclide *> ShieldingSelect::selfAttenNuclides() const
{
  // We'll use a set at first to make sure we dont return duplicates (not necessary, just
  //  a sanity check).
  set<const SandiaDecay::Nuclide *> answer;
  const shared_ptr<const Material> &mat = m_currentMaterial;

  if( !mat )
    return vector<const SandiaDecay::Nuclide *>();

  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const SandiaDecay::Element * const el = elDiv.first;
    const vector<WWidget *> &children = elDiv.second->children();
    assert( el );
    
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb || !cb->useAsSource() )
        continue;

      const SandiaDecay::Nuclide *iso = cb->isotope();
      if( !iso )
        continue;

      assert( !el || (el->atomicNumber == iso->atomicNumber) );
      
      //We'll check to make sure this nuclide might be in the material before adding, as a
      //  sanity check, that should never fail...
      bool material_has_nuc = false;
      for( const Material::NuclideFractionPair &nfp : mat->nuclides )
        material_has_nuc |= (nfp.first && (nfp.first==iso));
      
      for( const Material::ElementFractionPair &nfp : mat->elements )
        material_has_nuc |= (nfp.first && (nfp.first==el));
      
      assert( material_has_nuc );
      
      if( material_has_nuc )
      {
        //We shouldnt be encountering the same nuclide multiple times I dont think, so we'll check
        assert( !answer.count(iso) );
        answer.insert( iso );
      }
    }//for( WWidget *child : children )
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )

  return vector<const SandiaDecay::Nuclide *>( begin(answer), end(answer) );
}//std::vector<const SandiaDecay::Nuclide *> selfAttenNuclides() const


vector<const SandiaDecay::Element *> ShieldingSelect::elementsFittingNonSourceComponent() const
{
  vector<const SandiaDecay::Element *> answer;
  
  for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  {
    const SandiaDecay::Element * const el = elDiv.first;
    const vector<WWidget *> &children = elDiv.second->children();
    assert( el );
    
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( cb && cb->isotope() && cb->useAsSource() )
      {
        answer.push_back( el );
        break;
      }
    }//for( WWidget *child : children )
  }//for( const ElementToNuclideMap::value_type &elDiv : m_sourceIsotopes )
  
  return vector<const SandiaDecay::Element *>( begin(answer), end(answer) );
}//vector<const SandiaDecay::Element *> ShieldingSelect::elementsFittingNonSourceComponent() const


double ShieldingSelect::nuclidesMassFractionInElementOfMaterial( const SandiaDecay::Nuclide * const iso,
                                                        const std::shared_ptr<const Material> &mat )
{
  assert( mat );
  assert( iso );
  if( !mat || !iso )
    throw std::runtime_error( "nuclidesMassFractionInElementOfMaterial: Invalid material or isotope." );

  const SandiaDecay::SandiaDecayDataBase * const db = DecayDataBaseServer::database();
  const SandiaDecay::Element * const element = db->element( iso->atomicNumber );
  
  //Make sure the material has the isotope requested to add
  const vector< Material::NuclideFractionPair > &nuclides = mat->nuclides;
  const vector< Material::ElementFractionPair > &elements = mat->elements;

  double nuclidesMassFraction = 0.0;
  double elementsMassFraction = 0.0;
  
  bool hasNuclide = false, hasElement = false;
  for( const Material::ElementFractionPair &efp : elements )
  {
    if( efp.first == element )
    {
      hasElement = true;
      elementsMassFraction += efp.second;
    }//if( efp.first == element )
  }//for( const Material::ElementFractionPair &efp : elements )

  for( const Material::NuclideFractionPair &efp : nuclides )
  {
    if( efp.first->atomicNumber == iso->atomicNumber )
      elementsMassFraction += efp.second;

    if( efp.first == iso )
    {
      hasNuclide = true;
      nuclidesMassFraction += efp.second;
    }//if( efp.first == iso )
  }//for( const Material::NuclideFractionPair &efp : nuclides )

  if( !hasNuclide && !hasElement )
    throw runtime_error( "Material doesnt contain nuclide" );

  if( hasElement && !hasNuclide )
  {
    bool hasNaturalAbundance = false;
    const vector<SandiaDecay::NuclideAbundancePair> &natural_isos = element->isotopes;

    for( const SandiaDecay::NuclideAbundancePair &i : natural_isos )
    {
      hasNaturalAbundance |= (i.abundance != 0.0);
      if( i.nuclide == iso )
        nuclidesMassFraction += i.abundance * elementsMassFraction;
    }

    // If no natural abundance, then assign to 
    if( !hasNaturalAbundance )
    {
      vector<const SandiaDecay::Nuclide *> nucs_of_el = db->nuclides(element);
      nuclidesMassFraction = elementsMassFraction / nucs_of_el.size();
    }
  }//if( hasElement )

  if( elementsMassFraction == 0.0 )
    return 0.0;

  return nuclidesMassFraction;
}//double nuclidesMassFractionInElementOfMaterial( const SandiaDecay::Nuclide *iso )


double ShieldingSelect::elementsMassFractionInMaterial( const SandiaDecay::Element * const element,
                                             const std::shared_ptr<const Material> &mat )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  assert( mat );
  assert( element );
  
  if( !mat || !element || !db )
    throw std::runtime_error( "elementsMassFractionInMaterial: Invalid material or isotope." );
  
  //Make sure the material has the isotope requested to add
  const vector<Material::NuclideFractionPair> &nuclides = mat->nuclides;
  const vector<Material::ElementFractionPair> &elements = mat->elements;

  double elementsMassFraction = 0.0;
  bool hasNuclide = false, hasElement = false;
  
  for( const Material::ElementFractionPair &efp : elements )
  {
    if( efp.first == element )
    {
      hasElement = true;
      elementsMassFraction += efp.second;
    }
  }//for( const Material::ElementFractionPair &efp : elements )

  for( const Material::NuclideFractionPair &efp : nuclides )
  {
    if( efp.first->atomicNumber == element->atomicNumber )
    {
      hasNuclide = true;
      elementsMassFraction += efp.second;
    }
  }//for( const Material::NuclideFractionPair &efp : nuclides )

  if( !hasNuclide && !hasElement )
    throw runtime_error( "Material doesnt contain element" );

  return elementsMassFraction;
}//double elementsMassFractionInMaterial( el, mat )


void ShieldingSelect::modelNuclideAdded( const SandiaDecay::Nuclide *iso )
{
  // Deal with trace sources
  // Make sure no trace sources are using nuc, and if they are, remove that trace source
  if( m_traceSources )
  {
    const vector<WWidget *> kids = m_traceSources->children();
    for( WWidget *w : kids )
    {
      TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( w );
      assert( src );
      if( src ) // JIC
        src->modelSourceAdded( iso );
    }//for( WWidget *w : traceSources )
  }//if( m_traceSources )
  
  setTraceSourceBtnStatus();
  
  if( !m_asSourceCBs )
    return;

  std::shared_ptr<const Material> mat = m_currentMaterial;

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
    massFrac = nuclidesMassFractionInElementOfMaterial( iso, mat );
  }catch(...)
  {
    return;
  }

  // We may need to add "other" non-src checkbox, and we also want to add all other nuclides
  //  before this one (e.g., keep it on the bottom)
  SourceCheckbox *other_src_cb = nullptr;
  
  if( m_sourceIsotopes.find(element) == m_sourceIsotopes.end() )
    m_sourceIsotopes[element] = new WContainerWidget( m_asSourceCBs );
  
  WContainerWidget *isotopeDiv = m_sourceIsotopes[element];
  
  double accounted_for_frac = 0.0;
  for( WWidget *widget : isotopeDiv->children() )
  {
    SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
    if( src && (src->isotope()==iso) )
      return;
    
    if( src && !src->isotope() )
      other_src_cb = src;
    
    if( src && src->isotope() && src->useAsSource() )
      accounted_for_frac += src->massFraction();
  }//for( WWidget *widget : isotopeDiv->children() )

  if( !other_src_cb )
  {
    const double other_frac = std::max( 1.0 - accounted_for_frac, 0.0 );
    other_src_cb = new SourceCheckbox( nullptr, other_frac, isotopeDiv );
    other_src_cb->setFitMassFraction( false );
    other_src_cb->setLabelText( Wt::WString::tr("ss-non-src-frac-el").arg(element->symbol) );

    // The user input into mass fraction of this "other" non-src nuclide is disabled, so we
    //  dont need to hook it up to other_src_cb->handleIsotopicChange()....
    other_src_cb->fitMassFractionChecked().connect( 
                                boost::bind( &ShieldingSelect::handleFitMassFractionChanged, this,
                                            true, iso, element ) );
    other_src_cb->fitMassFractionUnChecked().connect(
                                boost::bind( &ShieldingSelect::handleFitMassFractionChanged, this,
                                            false, iso, element ) );
    
    other_src_cb->fitMassFractionChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
    other_src_cb->fitMassFractionUnChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
    other_src_cb->massFractionChanged().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  }//if( !other_src_cb )
  
  SourceCheckbox *cb = new SourceCheckbox( iso, massFrac );
  isotopeDiv->insertWidget( isotopeDiv->count() - 1, cb );
  //isotopeDiv->insertBefore( cb, other_src_cb );
  cb->checked().connect( boost::bind( &ShieldingSelect::isotopeCheckedCallback, this, iso ) );
  cb->unChecked().connect( boost::bind( &ShieldingSelect::isotopeUnCheckedCallback, this, iso ) );
  cb->massFractionChanged().connect( boost::bind( &ShieldingSelect::handleIsotopicChange, this,
                                                 boost::placeholders::_1, iso ) );

  cb->fitMassFractionChecked().connect( boost::bind( &ShieldingSelect::handleFitMassFractionChanged, 
                                                    this, true, iso, element ) );
  cb->fitMassFractionUnChecked().connect( boost::bind( &ShieldingSelect::handleFitMassFractionChanged, 
                                                      this, false, iso, element ) );
  
  handleIsotopicChange( static_cast<float>(massFrac), iso );
  handleFitMassFractionChanged( cb->fitMassFraction(), iso, element );
  
  cb->checked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  cb->unChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  cb->massFractionChanged().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  cb->fitMassFractionChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );
  cb->fitMassFractionUnChecked().connect( this, &ShieldingSelect::handleUserChangeForUndoRedo );

  if( m_asSourceCBs->isHidden() )
    m_asSourceCBs->show();

  updateSelfAttenOtherNucFraction();
}//void modelNuclideAdded( const std::string &symol )


void ShieldingSelect::handleIsotopicChange( const float input_fraction,
                                           const SandiaDecay::Nuclide * const nuc )
{
  assert( nuc );
  if( !nuc )
    throw runtime_error( "ShieldingSelect::handleIsotopicChange: invalid nuclide." );
  
  // Clamp the fraction to be between 0, and 1
  const float fraction = std::min( std::max( input_fraction, 0.0f ), 1.0f );
  
  const shared_ptr<const Material> &mat = m_currentMaterial;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Element * const element = db->element( nuc->atomicNumber );
  assert( element && (element->atomicNumber == nuc->atomicNumber) );
  
  if( !mat || !nuc || !element || !db )
    return;

  ElementToNuclideMap::iterator isos = m_sourceIsotopes.find(element);
  assert( isos != end(m_sourceIsotopes) );
  if( isos == end(m_sourceIsotopes) )
    throw runtime_error( "Nuclide '" + nuc->symbol + "' is not in material '" + mat->name + "'" );

  
  {//Begin adjust things
    SourceCheckbox *non_src_cb = nullptr;
    double total_other_src_in_el = 0.0;
    bool setNucFraction = false, useNuclideAsSrc = false;
    for( WWidget *child : isos->second->children() )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb )
        continue;
      
      if( cb->isotope() == nuc )
      {
        setNucFraction = true;
        useNuclideAsSrc = cb->useAsSource();
      }else if( !cb->isotope() )
      {
        non_src_cb = cb;
      }else if( cb->useAsSource() )
      {
        if( cb->isotope() != nuc )
          total_other_src_in_el += cb->massFraction();
      }//if( cb && !cb->isotope() ) / else
    }//for( loop over all self-atten source inputs for this element )
    
    assert( setNucFraction );
    if( !setNucFraction )
      throw std::logic_error( "Failed to set mass fraction for " + nuc->symbol + " - shouldnt have happened" );
    
    assert( non_src_cb );
    if( !non_src_cb )
      throw runtime_error( "ShieldingSelect::handleIsotopicChange(): failed to find SourceCheckbox for 'other' non-src nuclides" );
    
    if( !useNuclideAsSrc )
      return;
    
    if( (total_other_src_in_el + fraction) > 1.0 )
    {
      non_src_cb->setMassFraction( 0.0, 0.0 );
      const double other_src_amount = std::max(0.0, 1.0 - fraction);
      
      for( WWidget *child : isos->second->children() )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->useAsSource() && cb->isotope() && (cb->isotope() != nuc) )
        {
          const double prev_frac = cb->massFraction();
          const double multiple = ((fraction == 1.0) || (total_other_src_in_el <= 0.0))
                                    ? 0.0 : ((1.0 - fraction) / total_other_src_in_el);
          const double new_frac = multiple * prev_frac;
          cb->setMassFraction( new_frac, 0.0 );
        }
      }//for( const auto &el_widgets = m_sourceIsotopes )
      
      total_other_src_in_el = other_src_amount;
    }else
    {
      // We dont really have to do this, `updateSelfAttenOtherNucFraction()` will do it as well
      const double non_src_amount = std::max(0.0, 1.0 - fraction - total_other_src_in_el);
      non_src_cb->setMassFraction( non_src_amount, 0.0 );
    }//if( total_src_in_mat > 1.0 ) / else
  }//for( loop over all elements that have a self-atten source )
  
  updateSelfAttenOtherNucFraction();

  materialModified().emit( this );
}//void handleIsotopicChange( double fraction, const SandiaDecay::Nuclide *nuc )


void ShieldingSelect::handleFitMassFractionChanged( const bool fit,
                                  const SandiaDecay::Nuclide * const nuclide,
                                  const SandiaDecay::Element *el )
{
  //Here, we need to make sure:
  //  - At least two sources are selected to fit mass fraction, or zero are selected
  //  - The "other" non-source component is either showing or not showing
  //  - Adjust "other" amount based on what is currently checked
  //  - Hide "other" non-src cb, if no sources for this element are being used.
  
  updateSelfAttenOtherNucFraction();
}//void handleFitMassFractionChanged(...)


void ShieldingSelect::updateSelfAttenOtherNucFraction()
{
  for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  {
    const SandiaDecay::Element *el = vt.first;
    assert( el );
    if( !el )
      continue;

    Wt::WText *warning = nullptr;
    SourceCheckbox *non_src_frac_cb = nullptr;

    int num_as_src = 0, num_fit_frac = 0;
    double nuc_frac_in_el = 0.0;
    
    const vector<WWidget *> children = vt.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb )
      {
        WText *w = dynamic_cast<Wt::WText *>( child );
        if( w && w->hasStyleClass( "warning" ) )
          warning = w;
        continue;
      }

      num_fit_frac += cb->fitMassFraction();

      if( cb->isotope() && cb->useAsSource() )
      {
        ++num_as_src;
        nuc_frac_in_el += cb->massFraction();
      }
      
      if( !cb->isotope() )
      {
        assert( !non_src_frac_cb ); // Make sure we only have one non-src cb
        non_src_frac_cb = cb;
      }
    }//for( WWidget *child : children )
    
    assert( !num_as_src || non_src_frac_cb );

    // Check if we have more source fraction than 1.0, and if so, normalize all the source fractions
    //  and set the "other" non-src fraction to be zero
    if( num_as_src && (nuc_frac_in_el > 1.0) )
    {
      nuc_frac_in_el = 1.0;
      for( WWidget *child : children )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->useAsSource() && cb->isotope() )
        {
          const double new_frac = cb->massFraction() / nuc_frac_in_el;
          cb->setMassFraction( new_frac, 0.0 );
        }
      }//for( WWidget *child : children )
    }//if( num_as_src && (nuc_frac_in_el > 1.0) )

    if( !num_as_src && non_src_frac_cb )
      non_src_frac_cb->setFitMassFraction(false);

    assert( non_src_frac_cb || (num_as_src == 0) );
    if( non_src_frac_cb )
    {
      const double other_frac = 1.0 - nuc_frac_in_el;
      non_src_frac_cb->setMassFraction( other_frac, 0.0 );
    }

    if( (num_fit_frac == 1) && !warning )
    {
      warning = new Wt::WText( WString::tr("ss-warning-only-one-src-frac").arg(el->name), vt.second );
      warning->setStyleClass( "warning" );
      // Alternatively, if there is only one source fraction that is being fit, then we need to fit the "other" non-src fraction
      //  to make sure the total source fraction is 1.0
      //if( non_src_frac_cb )
      //  non_src_frac_cb->setFitMassFraction( true );
    }
    
    if( warning )
      warning->setHidden( (num_fit_frac != 1) );

    if( non_src_frac_cb )
      non_src_frac_cb->setHidden( (num_as_src == 0) );
  }//for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
}//void updateSelfAttenOtherNucFraction();


void ShieldingSelect::setMassFractionDisplaysToMaterial( std::shared_ptr<const Material> mat )
{
  if( !mat )
    return;

  assert( !m_fixedGeometry );
  if( m_fixedGeometry )
    throw logic_error( "Cant set mass fraction when fixed geometry" );
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  //Lets go through and update the values displayed to the user of the
  //  isotopics
  for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  {
    double frac_accounted_for = 0.0;
    
    const vector<WWidget *> children = vt.second->children();
    for( WWidget *child : children )
    {
      SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
      if( !cb )
        continue;

      const SandiaDecay::Nuclide * const nuc = cb->isotope();
      assert( nuc );
      if( !nuc )
        continue;
      
      double massFrac = 0.0;
      try
      {
        massFrac = nuclidesMassFractionInElementOfMaterial( nuc, mat ); //
      }catch(...)  //hopefully this never happens
      {
        passMessage( "There has been an unexpected internal error in"
                     " dealing with the material change, you may want to"
                     " try re-selecting the material, as well as checking the"
                     " calculation log before trusting the results.",
                     WarningWidget::WarningMsgHigh );
        cerr << "\nShieldingSelect::setMassFractionDisplaysToMaterial(...)\n"
        "\tSerious programming error here.\n" << endl;
      }//try/catch

      frac_accounted_for += massFrac;
      
      cb->setMassFraction( massFrac, 0.0 );
      
      // Lets set a tool tip with a little more information
      string tt = nuc->symbol + " is " + SpecUtils::printCompact(100.0*massFrac, 4) + "%";
      if( massFrac < 0.01 )
      {
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), " (%.0f PPM)", 1.0E6*massFrac );
        tt += buffer;
      }
      tt += " of the shielding material by mass.";
      
      double el_mass_frac = 0.0;
      const SandiaDecay::Element * const el = db->element( nuc->atomicNumber );
      if( el )
      {
        el_mass_frac = elementsMassFractionInMaterial( el, mat );
        const double nuc_frac_of_el = massFrac / el_mass_frac;
        const string percent_of_el_by_mass_str = SpecUtils::printCompact(100.0*nuc_frac_of_el, 4);
        tt += "\t" + nuc->symbol + " is " + percent_of_el_by_mass_str + "%";
        if( nuc_frac_of_el < 0.01 )
        {
          char buffer[64] = { '\0' };
          snprintf( buffer, sizeof(buffer), " (%.0f PPM)", 1.0E6*nuc_frac_of_el );
          tt += buffer;
        }
        tt += " of the " + el->name + " by mass.";
      }//if( el )
      
      // TODO: add PPM by number of atoms, if less than
      
      cb->setToolTip( tt );
    }//for( WWidget *child : children )
  }//for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
  
  checkAndUpdateMassFractionCanFit();
  updateSelfAttenOtherNucFraction();
}//void setMassFractionDisplaysToMaterial()


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
    m_materialEdit->setText( WString::tr("ss-generic") );
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
      m_atomicNumberEdit->setValue( mat->massWeightedAtomicNumber() );
      m_arealDensityEdit->setValue( ad*PhysicalUnits::cm2/PhysicalUnits::g );
    }else
    {
      m_atomicNumberEdit->setValue( 26.0f );
      m_arealDensityEdit->setValue( 0.0f );
    }
    
    if( m_forFitting )
    {
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          assert( m_fitArealDensityCB && m_fitThicknessCB );
          m_fitArealDensityCB->setChecked( m_fitThicknessCB->isChecked() );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          assert( m_fitArealDensityCB && m_fitCylRadiusCB && m_fitCylLengthCB );
          m_fitArealDensityCB->setChecked( m_fitCylRadiusCB->isChecked()
                                          || m_fitCylLengthCB->isChecked() );
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          assert( m_fitArealDensityCB && m_fitRectWidthCB && m_fitRectHeightCB && m_fitRectDepthCB );
          m_fitArealDensityCB->setChecked( m_fitRectWidthCB->isChecked()
                                          || m_fitRectHeightCB->isChecked()
                                          || m_fitRectDepthCB->isChecked() );
          break;
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
    }//if( m_forFitting )
  }else
  {
    m_toggleImage->setImageLink(Wt::WLink("InterSpec_resources/images/shield.png"));
    m_materialEdit->enable();
    displayInputsForCurrentGeometry();
    
    double ad = -1, an = -1;
    
    try
    {
      an = m_atomicNumberEdit->value();
      ad = m_arealDensityEdit->value() * PhysicalUnits::g / PhysicalUnits::cm2;
    }catch(...)
    {
      //
    }// try / catch
    
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
    
    if( (an <= 0.0) || (ad <=0.0) )
    {
      m_materialEdit->setText( "" );
      dist_edit->setText( "1.0 cm" );
    }else
    {
      //Get the atomic number closest to what the generic material was
      int atomicNum = (an > 0.5) ? static_cast<int>(0.5 + an) : (m_forFitting ? 16 : -1);
      
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


void ShieldingSelect::updateMaterialFromUserInputTxt()
{
  if( m_isGenericMaterial )
  {
    if( m_currentMaterial || !m_currentMaterialDescrip.empty() )
    {
      m_currentMaterial.reset();
      m_currentMaterialDescrip = "";
    }
    
    return;
  }//if( m_isGenericMaterial )
  
  const string text = SpecUtils::trim_copy( m_materialEdit->text().toUTF8() );
  if( m_currentMaterial && (text == m_currentMaterialDescrip) )
    return;

  if( text.empty() )
  {
    if( m_currentMaterial || !m_currentMaterialDescrip.empty() )
    {
      m_currentMaterial.reset();
      m_currentMaterialDescrip = "";
    }
    
    return;
  }//if( text.empty() )
  
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
}//void updateMaterialFromUserInputTxt()


void ShieldingSelect::handleMaterialChange()
{
  typedef pair<const SandiaDecay::Element *,float> ElementFrac;
  typedef pair<const SandiaDecay::Nuclide *,float> NuclideFrac;

  std::shared_ptr<const Material> newMaterial;
  std::shared_ptr<const Material> previousMaterial = m_currentMaterial;
  
  displayInputsForCurrentGeometry();
  setTraceSourceBtnStatus();
  
  if( m_isGenericMaterial )
  {
    m_materialSummary->setText( "" );
    m_materialEdit->setText( WString::tr("ss-generic") );
    m_materialEdit->disable();
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/atom_black.png") );
    
    if( m_traceSources )
    {
      const vector<WWidget *> kids = m_traceSources->children();
      for( WWidget *w : kids )
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
   
    updateMaterialFromUserInputTxt();
    
    newMaterial = m_currentMaterial;
    
    if( newMaterial )
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
                    "&rho;=%.3g g/cm<sup>3</sup>, <span style=\"text-decoration:overline\">AN</span>&#126;%.1f",
                    density, effAtomicNumber );
        }else
        {
          snprintf( summary, sizeof(summary),
                    "&rho;=%.3g g/cm<sup>3</sup>", density );
        }//if( m_forFitting ) / else
        
        
        
        tooltip += newMaterial->name + " " + WString::tr("ss-consists-of-mass-frac").toUTF8() + ":\n";

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
      passMessage( WString::tr("ss-err-invalid-material").arg(m_materialEdit->text()),
                  WarningWidget::WarningMsgInfo );
    }//if( material ) / else
    

//NOTE: can't add tooltip to this, causes WT error when toggling.  Can't fix.
//    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
//    const bool showToolTips = true;//UserPreferences::preferenceValue<bool>( "ShowTooltips", app->viewer() );
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
        
        if( !newMaterial )
        {
          const SandiaDecay::Nuclide * const oldNuc = src->nuclide();
          if( oldNuc )
          {
            src->deSelectNuclideNoEmit();
            src->nucChangedSignal().emit( src, oldNuc );
          }
        }else
        {
          src->updateForMaterialChange();
          m_activityFromVolumeNeedUpdating.emit( this, src->nuclide() );
        }
      }//for( WWidget *w : traceSources )
    }//if( m_traceSources )
  }//if( generic material ) / else
  
  
  //Now we need to update the activities for any isotopes that are
  if( !!newMaterial && m_asSourceCBs && (previousMaterial == newMaterial) )
  {
    //setMassFractionDisplaysToMaterial( newMaterial );

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
          removingIsotopeAsSource().emit( cb->isotope(), ShieldingSourceFitCalc::ModelSourceType::Intrinsic );
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

    m_asSourceCBs->setHidden( m_sourceIsotopes.empty() || m_fixedGeometry );
  }//if( previousMaterial != newMaterial )

  
#if( PERFORM_DEVELOPER_CHECKS )
  /*
   // This check is outdated - we no longer use the Material to track mass-fractions.
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
        massFrac = nuclidesMassFractionInElementOfMaterial( cb->isotope(), newMaterial );
      }catch(...)  //hopefully this never happens
      {
        stringstream msg;
        msg << "Failed to get nuclidesMassFractionInElementOfMaterial " << cb->isotope()->symbol
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
   */
#endif
  
  if( previousMaterial != newMaterial )
    m_materialModifiedSignal.emit( this );
  else
    m_materialChangedSignal.emit( this );
}//void handleMaterialChange()


void ShieldingSelect::handleUserChangedMaterialName()
{
  // When the user clicks on a  popup suggestions, first `m_materialEdit->changed()` is emitted
  // with the current text being an invalid material (just whatever they have typed so far),
  // then the `m_materialSuggest->activated()` signal is emitted (which updates the edits text),
  // then `m_materialEdit->changed()` changed again, with the now correct text in the edit.
  //
  // This function tries to avoid showing the user an error from the first `changed()` signal
  // where the text is just partial (by just putting back the previous valid text).
  
  const string text = SpecUtils::trim_copy( m_materialEdit->text().toUTF8() );
  if( text == m_currentMaterialDescrip )
    return;
  
  if( !text.empty()  )
  {
    const Material *mat = material( text );
    
    if( !mat )
    {
      //Check if `text` is the beginning of a shielding name, and if so, just set the text
      //  to the previous value.  If we dont do this, and its an invalid material, then the
      //  user will see an error message.
      const vector<string> &material_names = m_materialDB->names();
      for( const string &mat_name : material_names )
      {
        if( SpecUtils::istarts_with(mat_name, text) )
        {
          m_materialEdit->setText( WString::fromUTF8(m_currentMaterialDescrip) );
          return;
        }
      }//for( const string &mat_name : material_names )
    }//if( !mat )
  }//if( !text.empty() && (text != m_currentMaterialDescrip) )
  
  handleMaterialChange();
  
  handleUserChangeForUndoRedo();
}//void ShieldingSelect::handleUserChangedMaterialName()


void ShieldingSelect::handleUserChangeForUndoRedoWorker( const bool emit_change )
{
  auto xml_data = make_shared<string>();
  try
  {
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    rapidxml::print( std::back_inserter(*xml_data), doc, rapidxml::print_no_indenting );
    
    if( !m_prevState || ((*m_prevState) != (*xml_data)) )
    {
      if( emit_change )
        m_userChangedStateSignal.emit( this, m_prevState, xml_data );
      m_prevState = xml_data;
      //cout << "Storing ShieldingSelect taking " << xml_data->size() << " bytes mem" << endl;
    }else
    {
      cerr << "ShieldingSelect::handleUserChangeForUndoRedo(): identical state" << endl;
    }
  }catch( std::exception &e )
  {
    Wt::log("error") << "Failed to serialize in ShieldingSelect::handleUserChangeForUndoRedo().";
    // m_prevState = nullptr;
  }//try / catch
}//void handleUserChangeForUndoRedoWorker();


void ShieldingSelect::handleUserChangeForUndoRedo()
{
  if( !m_userChangedStateSignal.isConnected() )
    return;
  
  // We'll wait until after the current event loop to grab this widgets state;
  //  we could be here before all the changes have taken effect (e.g., I havent
  //  taken the care to make sure this function is called last in signal slots).

  Wt::WServer *server = Wt::WServer::instance();
  assert( server && wApp );
  if( server )
  {
    auto worker = wApp->bind( boost::bind(&ShieldingSelect::handleUserChangeForUndoRedoWorker, this, true) );
    server->post( wApp->sessionId(), worker );
  }
}//void handleUserChangeForUndoRedo()


ShieldingSourceFitCalc::ShieldingInfo ShieldingSelect::toShieldingInfo() const
{
  ShieldingSourceFitCalc::ShieldingInfo answer;
  
  answer.m_isGenericMaterial = m_isGenericMaterial;
  answer.m_forFitting = m_forFitting;
  
  if( m_isGenericMaterial )
  {
    answer.m_geometry = GammaInteractionCalc::GeometryType::NumGeometryType;
    
    answer.m_dimensions[0] = atomicNumber();
    answer.m_dimensions[1] = arealDensity();

    if( m_forFitting )
    {
      answer.m_fitDimensions[0] = fitAtomicNumber();
      answer.m_fitDimensions[1] = fitArealDensity();
    }
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    answer.m_truthDimensions[0] = truthAN;
    answer.m_truthDimensionsTolerances[0] = truthANTolerance;
    answer.m_truthDimensions[1] = truthAD;
    answer.m_truthDimensionsTolerances[1] = truthADTolerance;
#endif
  }else
  {
    answer.m_geometry = m_geometry;
    answer.m_material = m_currentMaterial;
    
    // Note: not calling `bool ShieldingSelect::fitThickness() const`, and similar, below, because
    //       some of these check if the checkbox is visible as well, which wont be the case if
    //       we are doing the sanity check of round tripping without the widget in the DOM
    switch( m_geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        answer.m_dimensions[0] = thickness();
        if( m_forFitting )
          answer.m_fitDimensions[0] = (m_fitThicknessCB && m_fitThicknessCB->isChecked());  //fitThickness();
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        answer.m_dimensions[0] = cylindricalRadiusThickness();
        answer.m_dimensions[1] = cylindricalLengthThickness();
        if( m_forFitting )
        {
          answer.m_fitDimensions[0] = (m_fitCylRadiusCB && m_fitCylRadiusCB->isChecked()); //fitCylindricalRadiusThickness();
          answer.m_fitDimensions[1] = (m_fitCylLengthCB && m_fitCylLengthCB->isChecked()); //fitCylindricalLengthThickness();
        }
        break;
        
      case GammaInteractionCalc::GeometryType::Rectangular:
        answer.m_dimensions[0] = rectangularWidthThickness();
        answer.m_dimensions[1] = rectangularHeightThickness();
        answer.m_dimensions[2] = rectangularDepthThickness();
        if( m_forFitting )
        {
          answer.m_fitDimensions[0] = (m_fitRectWidthCB && m_fitRectWidthCB->isChecked());   //fitRectangularWidthThickness();
          answer.m_fitDimensions[1] = (m_fitRectHeightCB && m_fitRectHeightCB->isChecked()); //fitRectangularHeightThickness();
          answer.m_fitDimensions[2] = (m_fitRectDepthCB && m_fitRectDepthCB->isChecked());   //fitRectangularDepthThickness();
        }
        break;
        
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )

#if( INCLUDE_ANALYSIS_TEST_SUITE )
    answer.m_truthDimensions[0] = truthThickness;
    answer.m_truthDimensionsTolerances[0] = truthThicknessTolerance;
    answer.m_truthDimensions[1] = truthThicknessD2;
    answer.m_truthDimensionsTolerances[1] = truthThicknessD2Tolerance;
    answer.m_truthDimensions[2] = truthThicknessD3;
    answer.m_truthDimensionsTolerances[2] = truthThicknessD3Tolerance;
    answer.m_truthFitMassFractions = truthFitMassFractions;
#endif
  }//if( m_isGenericMaterial ) / else
  
  
  // Self-atten source stuff
  const map<const SandiaDecay::Element *,vector<NucMasFrac>> selfAttens = sourceNuclideMassFractions();
  answer.m_nuclideFractions_ = selfAttens;
  
  // Trace source stuff
  if( m_traceSources )
  {
    for( WWidget *w : m_traceSources->children() )
    {
      const TraceSrcDisplay *src = dynamic_cast<const TraceSrcDisplay *>( w );
      assert( src );
      
      if( src && src->nuclide() )
        answer.m_traceSources.push_back( src->toTraceSourceInfo() );
    }//for( WWidget *w : m_traceSources->children() )
  }//if( m_traceSources )
  
  return answer;
}//ShieldingSourceFitCalc::ShieldingInfo toShieldingInfo() const


void ShieldingSelect::fromShieldingInfo( const ShieldingSourceFitCalc::ShieldingInfo &info )
{
  
  if( m_forFitting != info.m_forFitting )
    throw runtime_error( "ShieldingSelect m_forFitting must be same as XML being deserialized" );
  

  if( info.m_isGenericMaterial )
  {
    m_geometry = GeometryType::Spherical;
    
    if( !m_isGenericMaterial )
      handleToggleGeneric();
    
    m_atomicNumberEdit->setValue( info.m_dimensions[0] );
    m_arealDensityEdit->setValue( info.m_dimensions[1] * PhysicalUnits::cm2 / PhysicalUnits::gram );
    
    if( m_forFitting && m_fitArealDensityCB && m_fitAtomicNumberCB )
    {
      m_fitAtomicNumberCB->setChecked( info.m_fitDimensions[0] );
      m_fitArealDensityCB->setChecked( info.m_fitDimensions[1] );
    }//if( m_fitArealDensityCB )
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    truthAN  = info.m_truthDimensions[0];
    truthANTolerance = info.m_truthDimensionsTolerances[0];
    truthAD  = info.m_truthDimensions[1];
    truthADTolerance = info.m_truthDimensionsTolerances[1];
#endif
  }else
  {
    const bool wasGeneric = m_isGenericMaterial;
    const bool geom_changed = (m_geometry != info.m_geometry);
    if( info.m_geometry == GeometryType::NumGeometryType )
      throw runtime_error( "Invalid geometry type." );
    
    m_geometry = info.m_geometry;

    if( m_isGenericMaterial )
      handleToggleGeneric();

    if( geom_changed || wasGeneric )
      displayInputsForCurrentGeometry();
    
    // TODO: we could maually go in and get thickness text into the XML, so it is exact same as the user entered
    WLineEdit *dim_edits[3] = { nullptr, nullptr, nullptr };
    WCheckBox *dim_cb[3] = { nullptr, nullptr, nullptr };
    
    
    int required_dim = 0;
    switch( info.m_geometry )
    {
      case GeometryType::Spherical:
        required_dim = 1;
        dim_edits[0] = m_thicknessEdit;
        dim_cb[0] = m_fitThicknessCB;
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        required_dim = 2;
        dim_edits[0] = m_cylRadiusEdit;
        dim_edits[1] = m_cylLengthEdit;
        dim_cb[0]    = m_fitCylRadiusCB;
        dim_cb[1]    = m_fitCylLengthCB;
        break;
        
      case GeometryType::Rectangular:
        required_dim = 3;
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
    for( int i = 0; i < required_dim; ++i )
    {
      assert( dim_edits[i] );
      assert( !m_forFitting || dim_cb[i] );
      
      const string valstr = PhysicalUnits::printToBestLengthUnits( info.m_dimensions[i], 4 );
      dim_edits[i]->setText( WString::fromUTF8(valstr) );
      if( m_forFitting )
        dim_cb[i]->setChecked( info.m_fitDimensions[i] );
    }//for( int i = 0; i < required_dim; ++i )
    
      
    const string material_name = info.m_material ? info.m_material->name : string();
    m_materialEdit->setValueText( WString::fromUTF8(material_name) );
    
    handleMaterialChange();
    
    // We dont have uncertainties of mass fraction, so to use `setMassFractions(...)` we need to
    //  slightly convert formats from `info.m_nuclideFractions` to the format with uncertainty
    map<const SandiaDecay::Element *,vector<MassFracInfo>> fractions;
    for( const auto &i : info.m_nuclideFractions_ )
    {
      for( const auto &vals : i.second )
      {
        MassFracInfo info;
        info.m_nuclide = get<0>(vals);
        info.m_fraction = get<1>(vals);
        info.m_frac_uncert = 0.0;
        info.m_fit_mass_frac = get<2>(vals);
        info.m_use_as_source = true;  //We are only tracking nuclides being used as sources, so this will be true
        
        fractions[i.first].push_back( info );
      }
    }
    setMassFractions( fractions );
    
    //Now set the check boxes to make all the source nuclides called out in the
    //  XML as actual src nuclides, since we only saved nuclides we actually
    //  wanted to use as source nuclides.
    if( m_currentMaterial )
    {
      const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
      set<SourceCheckbox *> setSources, allSources;
      for( const auto &el_nucfrac : info.m_nuclideFractions_ )
      {
        const SandiaDecay::Element * const el = el_nucfrac.first;
        assert( el );
        assert( !el_nucfrac.second.empty() );
        
        for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucfrac.second )
        {
          const SandiaDecay::Nuclide * const nuc = get<0>(nuc_info);
          const double massFrac = get<1>(nuc_info);
          const bool fit_mass_frac = get<2>(nuc_info);
          
          for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
          {
            // The "other" non-src fraction is nullptr for all elements, so lets
            //  make sure to match element, as well as nuclide.
            if( etnm.first != el )
              continue;

            for( WWidget *widget : etnm.second->children() )
            {
              SourceCheckbox *src = dynamic_cast<SourceCheckbox *>( widget );
              if( src )
                allSources.insert( src );
              
              if( src && (nuc == src->isotope()) )
              {
                setSources.insert( src );
                if( nuc )
                  src->setUseAsSource( true );
                src->setFitMassFraction( fit_mass_frac );
                src->setMassFraction( massFrac, 0.0 );
              }
            }
          }//for( const ElementToNuclideMap::value_type &etnm : m_sourceIsotopes )
        }//for( const tuple<const SandiaDecay::Nuclide *,double,bool> &nuc_info : el_nucfrac.second )
      }//for( const auto &el_nucfrac : info.m_nuclideFractions_ )
      
      // Make sure SourceCheckbox present, that havent been set, are all set to not-sources.
      for( SourceCheckbox *src : allSources )
      {
        if( !setSources.count(src) )
        {
          src->setUseAsSource( false );
          src->setFitMassFraction( false );
          if( !src->isotope() )
            src->setMassFraction( 0.0, 0.0 );
        }//if( !setSources.count(src) )
      }//for( SourceCheckbox *src : allSources )
      
      // Get rid of all the trace sources - I dont *think* we need to emit signals and all that
      if( m_traceSources )
        m_traceSources->clear();
      
      for( const ShieldingSourceFitCalc::TraceSourceInfo &trace : info.m_traceSources )
      {
        addTraceSource();
        assert( m_traceSources );
        
        TraceSrcDisplay *src = dynamic_cast<TraceSrcDisplay *>( m_traceSources->children().back() );
        assert( src );
        
        // When we de-serialize we assume the source fitting model already has the nuclides
        //  available to become trace sources
        if( src )  //JIC
          src->fromTraceSourceInfo( trace );
      }//for( const TraceSourceInfo &trace : info.m_traceSources )
    }//if( m_currentMaterial )
    
    /*
     blah blah blah blah blah
    if( m_forFitting && m_fitMassFrac )
      m_fitMassFrac->setChecked( info.m_fitMassFrac );
    */
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    truthThickness = info.m_truthDimensions[0];
    truthThicknessTolerance = info.m_truthDimensionsTolerances[0];
    truthThicknessD2 = info.m_truthDimensions[1];
    truthThicknessD2Tolerance = info.m_truthDimensionsTolerances[1];
    truthThicknessD3 = info.m_truthDimensions[2];
    truthThicknessD3Tolerance = info.m_truthDimensionsTolerances[2];
    truthFitMassFractions = info.m_truthFitMassFractions;
#endif
  }//if( info.m_isGenericMaterial ) / else
}//void fromShieldingInfo( const ShieldingSourceFitCalc::ShieldingInfo &info )


void ShieldingSelect::serialize( rapidxml::xml_node<char> *parent_node ) const
{
  const ShieldingSourceFitCalc::ShieldingInfo info = toShieldingInfo();
  rapidxml::xml_node<char> * const added_node = info.serialize( parent_node );
  
  // Lets update the XML to have the user inputted text, to help exactly round-trip things
  if( !m_isGenericMaterial )
  {
    rapidxml::xml_node<char> * const material_node = XML_FIRST_NODE(added_node, "Material");
    
    auto setval = [material_node]( const char *name, Wt::WLineEdit *edit ){
      assert( material_node && name && edit );
      rapidxml::xml_node<char> *node = material_node ? material_node->first_node(name) : nullptr;
      assert( node );
      if( node )
      {
        const string valstr = edit->text().toUTF8();
        char *val = node->document()->allocate_string( valstr.c_str(), valstr.size() + 1 );
        node->value( val, valstr.size() );
      }
    };
    
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        setval( "Thickness", m_thicknessEdit );
        setval( "SphericalThickness", m_thicknessEdit );
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        setval( "CylinderRadiusThickness", m_cylRadiusEdit );
        setval( "CylinderLengthThickness", m_cylLengthEdit );
        break;
        
      case GeometryType::Rectangular:
        setval( "RectangularWidthThickness", m_rectWidthEdit );
        setval( "RectangularHeightThickness", m_rectHeightEdit );
        setval( "RectangularDepthThickness", m_rectDepthEdit );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
  }//if( !m_isGenericMaterial )
  
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  try
  {
    ShieldingSourceFitCalc::ShieldingInfo decoded_info;
    decoded_info.deSerialize( added_node, m_materialDB );
    ShieldingSourceFitCalc::ShieldingInfo::equalEnough( info, decoded_info );
  }catch( std::exception &e )
  {
    const string msg = "Error roundtripping ShieldingSelect to XML and back: " + string(e.what());
    log_developer_error( __func__, msg.c_str() );
    assert( 0 );
  }//try / catch
  
  /** Function to test that a ShieldingSelect is same as *this* - at least the components that get encoded to a URL. */
  auto testShieldingSelectPartiallySameAsOrig = [this]( const ShieldingSelect &test ){
    // We dont currently encode trace or self-attenuating source info
    //  Wt::WContainerWidget *m_asSourceCBs;
    //  ElementToNuclideMap m_sourceIsotopes;
    //  Wt::WContainerWidget *m_traceSources;
    // So we wont test for, ATM
    
    if( test.m_isGenericMaterial != m_isGenericMaterial )
      throw runtime_error( "Generic vs Material doesnt match" );
    
    if( test.m_forFitting != m_forFitting )
      throw runtime_error( "For fitting indicator doesnt match" );
    
    if( m_isGenericMaterial )
    {
      if( test.m_arealDensityEdit->text() != m_arealDensityEdit->text() )
        throw runtime_error( "Areal Density text doesnt match" );
      
      if( test.m_atomicNumberEdit->text() != m_atomicNumberEdit->text() )
        throw runtime_error( "Atomic number text doesnt match" );
      
      if( m_forFitting && test.m_fitArealDensityCB->isChecked() != m_fitArealDensityCB->isChecked() )
        throw runtime_error( "Fitting for areal density doesnt match" );
      
      if( m_forFitting && test.m_fitAtomicNumberCB->isChecked() != m_fitAtomicNumberCB->isChecked() )
        throw runtime_error( "Fitting for atomic number doesnt match" );
    }else
    {
      const string test_material = test.m_materialEdit->text().toUTF8();
      const string orig_material = m_materialEdit->text().toUTF8();
      
      if( test_material != orig_material )
        throw runtime_error( "Material name text doesnt match" );
      
      if( test.m_geometry != m_geometry )
        throw runtime_error( "geometries dont match" );
      
      switch( m_geometry )
      {
        case GammaInteractionCalc::GeometryType::Spherical:
          if( test.m_thicknessEdit->text() != m_thicknessEdit->text() )
            throw runtime_error( "Thickness text doesnt match ('"
                                + test.m_thicknessEdit->text().toUTF8() + "' vs '"
                                + m_thicknessEdit->text().toUTF8() + "'" );
          if( m_forFitting && (test.m_fitThicknessCB->isChecked() != m_fitThicknessCB->isChecked()) )
            throw runtime_error( "Fitting for thickness doesnt match" );
          break;
          
        case GammaInteractionCalc::GeometryType::CylinderEndOn:
        case GammaInteractionCalc::GeometryType::CylinderSideOn:
          if( test.m_cylRadiusEdit->text() != m_cylRadiusEdit->text() )
            throw runtime_error( "Cyl Rad text doesnt match" );
          if( test.m_cylLengthEdit->text() != m_cylLengthEdit->text() )
            throw runtime_error( "Cyl Len text doesnt match" );
          if( m_forFitting && test.m_fitCylRadiusCB->isChecked() != m_fitCylRadiusCB->isChecked() )
            throw runtime_error( "Fitting for Cyl Rad doesnt match" );
          if( m_forFitting && test.m_fitCylLengthCB->isChecked() != m_fitCylLengthCB->isChecked() )
            throw runtime_error( "Fitting for Cyl Len doesnt match" );
          break;
          
        case GammaInteractionCalc::GeometryType::Rectangular:
          if( test.m_rectWidthEdit->text() != m_rectWidthEdit->text() )
            throw runtime_error( "Rect width text doesnt match" );
          if( test.m_rectHeightEdit->text() != m_rectHeightEdit->text() )
            throw runtime_error( "Rect height text doesnt match" );
          if( test.m_rectDepthEdit->text() != m_rectDepthEdit->text() )
            throw runtime_error( "Rect depth text doesnt match" );
          if( m_forFitting && test.m_fitRectWidthCB->isChecked() != m_fitRectWidthCB->isChecked() )
            throw runtime_error( "Fitting for Rect width doesnt match" );
          if( m_forFitting && test.m_fitRectHeightCB->isChecked() != m_fitRectHeightCB->isChecked() )
            throw runtime_error( "Fitting for Rect Height doesnt match" );
          if( m_forFitting && test.m_fitRectDepthCB->isChecked() != m_fitRectDepthCB->isChecked() )
            throw runtime_error( "Fitting for Rect Depth doesnt match" );
          break;
          
        case GammaInteractionCalc::GeometryType::NumGeometryType:
          assert( 0 );
          break;
      }//switch( m_geometry )
    }//}//if( m_isGenericMaterial ) / else
  };//testShieldingSelectPartiallySameAsOrig
    
  try
  {
    const string uri = encodeStateToUrl();
    
    ShieldingSelect test( m_materialDB, m_sourceModel, m_materialSuggest, m_shieldSrcDisp  );
    test.handleAppUrl( uri );
    testShieldingSelectPartiallySameAsOrig( test );
  }catch( std::exception &e )
  {
    const string msg = "Error roundtripping ShieldingSelect to URI and back: " + string(e.what());
    log_developer_error( __func__, msg.c_str() );
//    assert( 0 );
  }//try
  
  /*
   //This test will fail for trace-sources, as getting their total activity depends
   //  on them being inserted into m_shieldSrcDisp
  try
  {
    ShieldingSelect test( m_materialDB, m_sourceModel, m_materialSuggest, m_shieldSrcDisp );
    test.deSerialize(added_node);
    testShieldingSelectPartiallySameAsOrig( test );
  }catch( std::exception &e )
  {
    const string msg = "Error roundtripping ShieldingSelect to XML and back: " + string(e.what());
    log_developer_error( __func__, msg.c_str() );
    assert( 0 );
  }//try
   */
#endif  //PERFORM_DEVELOPER_CHECKS
  
}//void serialize( rapidxml::xml_document<> &doc ) const;


void ShieldingSelect::deSerialize( const rapidxml::xml_node<char> *shield_node,
                                  const bool is_fixed_geom_det )
{
  ShieldingSourceFitCalc::ShieldingInfo info;
  info.deSerialize( shield_node, m_materialDB );
  
  // If not for fitting, we shouldnt have intrinsic or trace sources
  assert( info.m_forFitting || info.m_traceSources.empty() );
  assert( info.m_forFitting || info.m_nuclideFractions_.empty() );
  
  if( is_fixed_geom_det )
  {
    info.m_nuclideFractions_.clear();
    info.m_traceSources.clear();
    info.m_geometry = GammaInteractionCalc::GeometryType::Spherical;
  }//if( is_fixed_geom_det )
  
  fromShieldingInfo( info );
  
  // Reach into the XML and set distances exactly equal to user input text
  if( !m_isGenericMaterial )
  {
    rapidxml::xml_node<char> * const material_node = XML_FIRST_NODE(shield_node, "Material");
    assert( material_node );
    
    auto getval = [material_node]( const char *name, Wt::WLineEdit *edit ){
      assert( material_node && name && edit );
      rapidxml::xml_node<char> *node = material_node ? material_node->first_node(name) : nullptr;
      assert( node );
      if( node && edit )
      {
        const string valstr = SpecUtils::xml_value_str( node );
        edit->setText( WString::fromUTF8(valstr) );
      }
    };
    
    switch( m_geometry )
    {
      case GeometryType::Spherical:
        getval( "SphericalThickness", m_thicknessEdit );
        break;
        
      case GeometryType::CylinderEndOn:
      case GeometryType::CylinderSideOn:
        getval( "CylinderRadiusThickness", m_cylRadiusEdit );
        getval( "CylinderLengthThickness", m_cylLengthEdit );
        break;
        
      case GeometryType::Rectangular:
        getval( "RectangularWidthThickness", m_rectWidthEdit );
        getval( "RectangularHeightThickness", m_rectHeightEdit );
        getval( "RectangularDepthThickness", m_rectDepthEdit );
        break;
        
      case GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    handleMaterialChange();
  }//if( !m_isGenericMaterial )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  const ShieldingSourceFitCalc::ShieldingInfo roundtripped = toShieldingInfo();
  try
  {
    ShieldingSourceFitCalc::ShieldingInfo::equalEnough( info, roundtripped );
  }catch( std::exception &e )
  {
    cerr << "ShieldingSelect::deSerialize - Failed to round-trip setting/getting: " << e.what() << endl;
    assert( 0 );
  }
#endif
  
  // Make sure the `m_prevState` will get updated
  Wt::WServer *server = Wt::WServer::instance();
  assert( server && wApp );
  if( m_userChangedStateSignal.isConnected() && server )
  {
    auto worker = wApp->bind( boost::bind(&ShieldingSelect::handleUserChangeForUndoRedoWorker, this, false) );
    server->post( wApp->sessionId(), worker );
  }//if( m_userChangedStateSignal.isConnected() && server )
}//void deSerialize( const rapidxml::xml_node<char> *shielding_node ) const


std::string ShieldingSelect::encodeStateToUrl() const
{
  const ShieldingSourceFitCalc::ShieldingInfo info = toShieldingInfo();
  string uri = info.encodeStateToUrl();
  
  // Update distances to exactly match user input for round-tripping
  auto updateval = [&uri]( const string name, WLineEdit *edit ){
    assert( edit );
    const auto pos = uri.find( "&" + name + "=" );
    assert( pos != string::npos );
    if( pos == string::npos )
      return;
    const string newval = edit ? edit->text().toUTF8() : string();
    
    const string pre = uri.substr(0, pos + 2 + name.size() );
    const auto nextstart = uri.find( "&", pos + 2 + name.size() );
    const string post = (nextstart == string::npos) ? string() : uri.substr(nextstart);
    
    uri = pre + newval + post;
  };//updateval...
  
  switch( m_geometry )
  {
    case GammaInteractionCalc::GeometryType::Spherical:
      if( m_isGenericMaterial )
      {
        updateval( "AD", m_arealDensityEdit );
        updateval( "AN", m_atomicNumberEdit );
      }else
      {
        updateval( "D1", m_thicknessEdit );
      }
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderEndOn:
    case GammaInteractionCalc::GeometryType::CylinderSideOn:
      updateval( "D1", m_cylRadiusEdit );
      updateval( "D2", m_cylLengthEdit );
      break;
      
    case GammaInteractionCalc::GeometryType::Rectangular:
      updateval( "D1", m_rectWidthEdit );
      updateval( "D2", m_rectHeightEdit );
      updateval( "D3", m_rectDepthEdit );
      break;
      
    case GammaInteractionCalc::GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch( m_geometry )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  // TODO: 20230705: remove this section of code after initial testing
  
  // "V=1&G=S&D1=1.2cm&N=Fe&NTRACE=3&TRACEN=1&V=1&N=U238&A=1.2uCi&T=total&F=1"
  // "V": version
  // "G": geometry
  // "F": for fitting; if not specified than false
  // "D1": "D2": Thickness, depth, etc
  // "FD1": "FD2": fit the corresponding dimensions
  // "N": material name
  // "AN": atomic number
  // "FAN": fit atomic number - if not specified than false
  // "AD": areal density
  // "FAD": fit areal density - if not specified than false
  // ...
  
  string answer = "V=1";
  
  if( m_forFitting )
    answer += "&F=1";
  
  if( m_isGenericMaterial )
  {
    answer += "&AD=" + m_arealDensityEdit->text().toUTF8();
    answer += "&AN=" + m_atomicNumberEdit->text().toUTF8();
    
    if( m_forFitting && m_fitAtomicNumberCB->isChecked() )
      answer += "&FAN=1";
    if( m_forFitting && m_fitArealDensityCB->isChecked() )
      answer += "&FAD=1";
  }else
  {
    std::string material_name = m_materialEdit->text().toUTF8();
    SpecUtils::ireplace_all(material_name, "#", "%23" );
    SpecUtils::ireplace_all(material_name, "&", "%26" );
    
    if( !m_currentMaterial )
      material_name = "";
    
    answer += "&N=" + material_name;
      
    switch( m_geometry )
    {
      case GammaInteractionCalc::GeometryType::Spherical:
        answer += "&G=S&D1=" + m_thicknessEdit->text().toUTF8();
        if( m_forFitting && m_fitThicknessCB->isChecked() )
          answer += "&FD1=1";
        break;
        
      case GammaInteractionCalc::GeometryType::CylinderEndOn:
      case GammaInteractionCalc::GeometryType::CylinderSideOn:
        if( m_geometry == GammaInteractionCalc::GeometryType::CylinderEndOn )
          answer += "&G=CE";
        else
          answer += "&G=CS";
        answer += "&D1=" + m_cylRadiusEdit->text().toUTF8();
        answer += "&D2=" + m_cylLengthEdit->text().toUTF8();
        if( m_forFitting && m_fitCylRadiusCB->isChecked() )
          answer += "&FD1=1";
        if( m_forFitting && m_fitCylLengthCB->isChecked() )
          answer += "&FD2=1";
        break;
        
      case GammaInteractionCalc::GeometryType::Rectangular:
        answer += "&G=R";
        answer += "&D1=" + m_rectWidthEdit->text().toUTF8();
        answer += "&D2=" + m_rectHeightEdit->text().toUTF8();
        answer += "&D3=" + m_rectDepthEdit->text().toUTF8();
        if( m_forFitting && m_fitRectWidthCB->isChecked() )
          answer += "&FD1=1";
        if( m_forFitting && m_fitRectHeightCB->isChecked() )
          answer += "&FD2=1";
        if( m_forFitting && m_fitRectDepthCB->isChecked() )
          answer += "&FD3=1";
        break;
        
      case GammaInteractionCalc::GeometryType::NumGeometryType:
        assert( 0 );
        break;
    }//switch( m_geometry )
    
    // TODO: encode self-attenuating and trace sources, and maybe "truth" value; with a string like "NTRACE=3&TRACEN=1&V=1&N=U238&A=1.2uCi&T=total&F=1"
  }//if( m_isGenericMaterial ) / else
  
  assert( uri == answer );
#endif
  
  return uri;
}//std::string encodeStateToUrl() const


void ShieldingSelect::handleAppUrl( std::string query_str )
{
  ShieldingSourceFitCalc::ShieldingInfo info;
  info.handleAppUrl( query_str, m_materialDB );
  fromShieldingInfo( info );
  
  // Update distances to exactly match user input for round-tripping
  auto updateval = [&query_str]( const string name, WLineEdit *edit ){
    assert( edit );
    const auto pos = query_str.find( "&" + name + "=" );
    assert( pos != string::npos );
    if( pos == string::npos )
      return;
    
    string val = query_str.substr( pos + 2 + name.size() );
    val = val.substr( 0, val.find( "&" ) );
    if( edit )
      edit->setText( WString::fromUTF8(val) );
  };//updateval...
  
  switch( m_geometry )
  {
    case GammaInteractionCalc::GeometryType::Spherical:
      if( m_isGenericMaterial )
      {
        updateval( "AD", m_arealDensityEdit );
        updateval( "AN", m_atomicNumberEdit );
      }else
      {
        updateval( "D1", m_thicknessEdit );
      }
      break;
      
    case GammaInteractionCalc::GeometryType::CylinderEndOn:
    case GammaInteractionCalc::GeometryType::CylinderSideOn:
      updateval( "D1", m_cylRadiusEdit );
      updateval( "D2", m_cylLengthEdit );
      break;
      
    case GammaInteractionCalc::GeometryType::Rectangular:
      updateval( "D1", m_rectWidthEdit );
      updateval( "D2", m_rectHeightEdit );
      updateval( "D3", m_rectDepthEdit );
      break;
      
    case GammaInteractionCalc::GeometryType::NumGeometryType:
      assert( 0 );
      break;
  }//switch( m_geometry )
  
  
#if( PERFORM_DEVELOPER_CHECKS )
  ShieldingSourceFitCalc::ShieldingInfo roundtrip = toShieldingInfo();
  try
  {
    ShieldingSourceFitCalc::ShieldingInfo::equalEnough( info, roundtrip );
  }catch( std::exception &e )
  {
    cerr << "ShieldingSelect::handleAppUrl: failed to roundtrip: " << e.what() << endl;
    //assert( 0 );
  }
#endif //#if( PERFORM_DEVELOPER_CHECKS )
}//void handleAppUrl( std::string query_str )
