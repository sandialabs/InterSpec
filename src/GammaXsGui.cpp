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
#include <vector>
#include <sstream>

#include <Wt/WText.h>
#include <Wt/WLabel.h>
#include <Wt/WLineEdit.h>
#include <Wt/WGridLayout.h>
#include <Wt/WPushButton.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WDoubleValidator.h>
#include <Wt/WRegExpValidator.h>
#include <Wt/WSuggestionPopup.h>
#include <Wt/WAbstractItemModel.h>

#include "InterSpec/AppUtils.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"


#if( USE_QR_CODES )
#include <Wt/Utils.h>

#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;


GammaXsGui::GammaXsGui( Wt::WSuggestionPopup *materialSuggestion,
                        InterSpec* viewer )
  : WContainerWidget(),
    m_energyEdit( nullptr ),
    m_energyValidator( nullptr ),
    m_materialEdit( nullptr ),
    m_materialSuggestion( materialSuggestion ),
    m_effectiveZ( nullptr ),
    m_totalAttenuation( nullptr ),
    m_compton( nullptr ),
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
    m_rayleigh( nullptr ),
#endif
    m_photoElectric( nullptr ),
    m_conversion( nullptr ),
    m_density( nullptr ),
    m_thickness( nullptr ),
    m_transmissionFraction( nullptr ),
    m_transmissionFractionVal( -1.0f ),
    m_layout( nullptr ),
    m_specViewer( viewer ),
    m_detectorDisplay( nullptr ),
    m_detectorDistanceLabel( nullptr ),
    m_detectorDistance( nullptr ),
    m_efficiencyLabel( nullptr ),
    m_efficiency( nullptr ),
    m_totalEfficiencyLabel( nullptr ),
    m_totalEfficiency( nullptr ),
    m_intrinsicEfficiencyLabel( nullptr ),
    m_intrinsicEfficiency( nullptr ),
    m_fractionalAngleLabel( nullptr ),
    m_fractionalAngle( nullptr )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;

  addStyleClass( "GammaXsGui" );

  wApp->useStyleSheet( "InterSpec_resources/GammaXsGui.css" );

  m_layout = setLayout( std::make_unique<WGridLayout>() );

  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_specViewer );

  m_specViewer->useMessageResourceBundle( "GammaXsGui" );

  {
    auto edit = std::make_unique<WLineEdit>( "100" );
    m_energyEdit = edit.get();
    m_energyEdit->setAutoComplete( false );
    m_energyEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_energyEdit->setAttributeValue( "autocorrect", "off" );
    m_energyEdit->setAttributeValue( "spellcheck", "off" );
#endif
    m_energyEdit->addStyleClass( "numberValidator"); //used to detect mobile keyboard
    m_energyEdit->changed().connect( this, &GammaXsGui::calculateCrossSections );
    m_energyEdit->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
    m_energyEdit->focussed().connect( this, &GammaXsGui::calculateCrossSections );
    m_energyEdit->blurred().connect( this, &GammaXsGui::calculateCrossSections );
    // Validator is a child of the edit widget, so create after parenting
    m_layout->addWidget( std::move(edit), 0, 1, 1, 1 );
  }
  m_energyValidator = std::make_shared<WDoubleValidator>( 1.0, 10000.0 );
  m_energyEdit->setValidator( m_energyValidator );

  int row = 0;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("energy-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  // m_energyEdit already added above at (0, 1)
  m_layout->addWidget( std::make_unique<WLabel>( "keV" ), row, 2, 1, 1, AlignmentFlag::Left );

  ++row;

  WLabel *matLabel = m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-material-formula") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto edit = std::make_unique<WLineEdit>( "C0.5H0.2Ni0.3" );
    m_materialEdit = edit.get();
    m_materialEdit->setAutoComplete( false );
    m_materialEdit->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_materialEdit->setAttributeValue( "autocorrect", "off" );
    m_materialEdit->setAttributeValue( "spellcheck", "off" );
#endif
    m_materialEdit->changed().connect( this, &GammaXsGui::handleMaterialChange );
    m_materialEdit->enterPressed().connect( this, &GammaXsGui::handleMaterialChange );
//  m_materialEdit->focussed().connect( this, &GammaXsGui::handleMaterialChange );
    m_materialEdit->blurred().connect( this, &GammaXsGui::handleMaterialChange );
    m_layout->addWidget( std::move(edit), row, 1, 1, 2 );
  }
  HelpSystem::attachToolTipOn( {matLabel, m_materialEdit}, WString::tr("gxsg-tt-material"), showToolTips );
  m_materialSuggestion->forEdit( m_materialEdit,
                  PopupTrigger::Editing | PopupTrigger::DropDownIcon );

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-total-atten-xs") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_totalAttenuation = w.get();
    m_totalAttenuation->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"cm\x00B2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "cm2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-compton-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_compton = w.get();
    m_compton->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"cm\x00B2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "cm2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif

#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-rayleigh-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_rayleigh = w.get();
    m_rayleigh->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"cm\x00B2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "cm2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif
#endif //#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-photoelec-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_photoElectric = w.get();
    m_photoElectric->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"cm\x00B2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "cm2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-pp-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_conversion = w.get();
    m_conversion->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"cm\x00B2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "cm2/g" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-mass-avrg-an") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>( "" );
    m_effectiveZ = w.get();
    m_effectiveZ->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 1 );
  }

  ++row;
  m_layout->addWidget( std::make_unique<WText>( WString::tr("gxsg-opt-atten") ), row, 0, 1, 3, AlignmentFlag::Bottom );

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-density-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto edit = std::make_unique<WLineEdit>();
    m_density = edit.get();
    m_density->setAutoComplete( false );
    m_density->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_density->setAttributeValue( "autocorrect", "off" );
    m_density->setAttributeValue( "spellcheck", "off" );
#endif
    m_density->addStyleClass( "numberValidator" ); //used to detect mobile keyboard
    m_density->changed().connect( this, &GammaXsGui::calculateCrossSections );
    m_density->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
    m_density->focussed().connect( this, &GammaXsGui::calculateCrossSections );
    m_density->blurred().connect( this, &GammaXsGui::calculateCrossSections );
    m_layout->addWidget( std::move(edit), row, 1, 1, 1 );
  }
  {
    auto doubValidator = std::make_shared<WDoubleValidator>();
    m_density->setValidator( doubValidator );
  }
#ifndef WT_NO_STD_WSTRING
  m_layout->addWidget( std::make_unique<WLabel>( L"g/cm\x00B3" ), row, 2, 1, 1, AlignmentFlag::Left );
#else
  m_layout->addWidget( std::make_unique<WLabel>( "g/cm3" ), row, 2, 1, 1, AlignmentFlag::Left );
#endif

  // Create shared distance validator for thickness and detector distance fields
  auto distValidator = std::make_shared<WRegExpValidator>( PhysicalUnits::sm_distanceRegex );
  distValidator->setFlags( Wt::RegExpFlag::MatchCaseInsensitive );

  ++row;
  m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-thickness-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto edit = std::make_unique<WLineEdit>( "1 cm" );
    m_thickness = edit.get();
    m_thickness->setAutoComplete( false );
    m_thickness->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_thickness->setAttributeValue( "autocorrect", "off" );
    m_thickness->setAttributeValue( "spellcheck", "off" );
#endif
    m_thickness->setValidator( distValidator );
    m_thickness->changed().connect( this, &GammaXsGui::calculateCrossSections );
    m_thickness->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
    m_thickness->focussed().connect( this, &GammaXsGui::calculateCrossSections );
    m_thickness->blurred().connect( this, &GammaXsGui::calculateCrossSections );
    m_layout->addWidget( std::move(edit), row, 1, 1, 1 );
  }
  HelpSystem::attachToolTipOn( m_thickness, WString::tr("gxsg-tt-thickness"), showToolTips );

  ++row;
  WLabel *transFracLabel = m_layout->addWidget( std::make_unique<WLabel>( WString::tr("gxsg-trans-frac-label") ), row, 0, 1, 1, AlignmentFlag::Left );
  {
    auto w = std::make_unique<WText>();
    m_transmissionFraction = w.get();
    m_transmissionFraction->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 2 );
  }
  HelpSystem::attachToolTipOn( {transFracLabel, m_transmissionFraction}, WString::tr("gxsg-tt-trans-frac"),
                              showToolTips );
  m_layout->setColumnStretch( 1, row );

  ++row;
  {
    SpecMeasManager *fileManager = m_specViewer->fileManager();
    auto dd = std::make_unique<DetectorDisplay>( m_specViewer, fileManager->model() );
    m_detectorDisplay = dd.get();
    m_detectorDisplay->setMinimumSize( WLength::Auto, WLength(1.5,WLength::Unit::FontEm) );
    m_layout->addWidget( std::move(dd), row, 0, 1, 3 );
  }

  ++row;
  {
    auto lbl = std::make_unique<WLabel>( WString::tr("Distance") );
    m_detectorDistanceLabel = lbl.get();
    m_layout->addWidget( std::move(lbl), row, 0, 1, 1, AlignmentFlag::Left );
  }
  {
    auto edit = std::make_unique<WLineEdit>("2 cm");
    m_detectorDistance = edit.get();
    m_detectorDistance->setAutoComplete( false );
    m_detectorDistance->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
    m_detectorDistance->setAttributeValue( "autocorrect", "off" );
    m_detectorDistance->setAttributeValue( "spellcheck", "off" );
#endif
    m_detectorDistance->setValidator( distValidator );
    m_detectorDistance->changed().connect( this, &GammaXsGui::updateDetectorCalc );
    m_detectorDistance->enterPressed().connect( this, &GammaXsGui::updateDetectorCalc );
    m_detectorDistance->focussed().connect( this, &GammaXsGui::updateDetectorCalc );
    m_detectorDistance->blurred().connect( this, &GammaXsGui::updateDetectorCalc );
    m_layout->addWidget( std::move(edit), row, 1, 1, 2 );
  }
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"This is the distance of the selected detector.", showToolTips );
  HelpSystem::attachToolTipOn( m_detectorDistance, WString::tr("gxsg-tt-distance"), showToolTips );

  ++row;
  {
    auto lbl = std::make_unique<WLabel>( WString::tr("gxsg-intrinsic-eff-label") );
    m_intrinsicEfficiencyLabel = lbl.get();
    m_layout->addWidget( std::move(lbl), row, 0, 1, 1, AlignmentFlag::Left );
  }
  {
    auto w = std::make_unique<WText>();
    m_intrinsicEfficiency = w.get();
    m_intrinsicEfficiency->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 2 );
  }
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"Intrinsic efficiency (in-peak detection efficiency of gammas striking detector face).", showToolTips );
  HelpSystem::attachToolTipOn( m_intrinsicEfficiency, WString::tr("gxsg-tt-intrinsic-eff"),
                              showToolTips );

  ++row;
  {
    auto lbl = std::make_unique<WLabel>( WString::tr("gxsg-solid-angle-frac-label") );
    m_fractionalAngleLabel = lbl.get();
    m_layout->addWidget( std::move(lbl), row, 0, 1, 1, AlignmentFlag::Left );
  }
  {
    auto w = std::make_unique<WText>();
    m_fractionalAngle = w.get();
    m_fractionalAngle->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 2 );
  }
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"Fractional solid angle of the selected detector at specified distance." , showToolTips );
  HelpSystem::attachToolTipOn( m_fractionalAngle, WString::tr("gxsg-tt-solid-angle-frac"),
                              showToolTips );

  ++row;
  {
    auto lbl = std::make_unique<WLabel>( WString::tr("gxsg-det-eff-label") );
    m_efficiencyLabel = lbl.get();
    m_layout->addWidget( std::move(lbl), row, 0, 1, 1, AlignmentFlag::Left );
  }
  {
    auto w = std::make_unique<WText>();
    m_efficiency = w.get();
    m_efficiency->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 2 );
  }
  HelpSystem::attachToolTipOn( {m_efficiencyLabel, m_efficiency}, WString::tr("gxsg-tt-det-eff"),
                              showToolTips );

  ++row;
  {
    auto lbl = std::make_unique<WLabel>( WString::tr("gxsg-total-eff-label") );
    m_totalEfficiencyLabel = lbl.get();
    m_layout->addWidget( std::move(lbl), row, 0, 1, 1, AlignmentFlag::Left );
  }
  {
    auto w = std::make_unique<WText>();
    m_totalEfficiency = w.get();
    m_totalEfficiency->addStyleClass( "GxsResult" );
    m_layout->addWidget( std::move(w), row, 1, 1, 2 );
  }
  HelpSystem::attachToolTipOn( m_totalEfficiencyLabel,
    "Transmission fraction times detection efficiency." , showToolTips );
  HelpSystem::attachToolTipOn( m_totalEfficiency, WString::tr("gxsg-tt-total-eff"),
                              showToolTips );

  m_specViewer->detectorChanged().connect( this, &GammaXsGui::handleDetectorChange );
  m_specViewer->detectorModified().connect( this, &GammaXsGui::handleDetectorChange );
  handleDetectorChange( m_detectorDisplay->detector() );

  calculateCrossSections();
}//GammaXsGui


GammaXsGui::~GammaXsGui()
{
  m_materialSuggestion->removeEdit( m_materialEdit );
}//~GammaXsGui()


void GammaXsGui::handleDetectorChange( std::shared_ptr<DetectorPeakResponse> det )
{
  m_detector = det;
  updateDetectorCalc();
}


void GammaXsGui::updateDetectorCalc()
{
  const bool hasDetector = !!m_detector;
  const bool fixed_geom = (m_detector && m_detector->isFixedGeometry());
  
  m_detectorDistanceLabel->setHidden( !hasDetector || fixed_geom );
  m_detectorDistance->setHidden( !hasDetector || fixed_geom );
  
  m_intrinsicEfficiencyLabel->setHidden( !hasDetector || fixed_geom );
  m_intrinsicEfficiency->setHidden( !hasDetector || fixed_geom );
  
  m_fractionalAngleLabel->setHidden( !hasDetector || fixed_geom );
  m_fractionalAngle->setHidden( !hasDetector || fixed_geom );
  
  m_efficiencyLabel->setHidden( !hasDetector );
  m_efficiency->setHidden( !hasDetector );
  
  m_totalEfficiencyLabel->setHidden( !hasDetector );
  m_totalEfficiency->setHidden( !hasDetector );
  
  if( !hasDetector )
    return;
  
  try
  {
    // efficiency
    char buffer[32];
    const string energystr = m_energyEdit->text().toUTF8();
    const float energy = static_cast<float>( std::stod( energystr ) );
    const double intrinsic_eff =  m_detector->intrinsicEfficiency( energy );
    
    if( m_detector->isFixedGeometry() )
    {
      snprintf( buffer, sizeof(buffer), "%.4g", intrinsic_eff );
      m_efficiency->setText( buffer );
      
      if( m_transmissionFractionVal >= 0.0f )
      {
        snprintf( buffer, sizeof(buffer), "%.4g", intrinsic_eff*m_transmissionFractionVal );
        m_totalEfficiency->setText( buffer );
      }else
      {
        m_totalEfficiency->setText( "---" );
      }
    }else
    {
      const string diststr = m_detectorDistance->text().toUTF8();
      const double dist = PhysicalUnits::stringToDistance( diststr );
      const double det_eff = m_detector->efficiency( energy, dist );
      snprintf( buffer, sizeof(buffer), "%.4g", det_eff );
      m_efficiency->setText( buffer );
      
      if( m_transmissionFractionVal >= 0.0f )
      {
        snprintf( buffer, sizeof(buffer), "%.4g", det_eff*m_transmissionFractionVal );
        m_totalEfficiency->setText( buffer );
      }else
      {
        m_totalEfficiency->setText( "---" );
      }
      
      // fractional solid angle
      const float diameter = m_detector->detectorDiameter();
      const double gfactor =  m_detector->fractionalSolidAngle( diameter, dist + m_detector->detectorSetback() );
      snprintf( buffer, sizeof(buffer), "%.4g", gfactor );
      m_fractionalAngle->setText( buffer );
      
      // absolute efficiency
      snprintf( buffer, sizeof(buffer), "%.4g", intrinsic_eff );
      m_intrinsicEfficiency->setText( buffer );
    }//if( fixed geometry ) / else
    
  }catch( std::runtime_error & )
  {
    //could not get all the measurements for the detector
    m_efficiency->setText( WString::tr("n/a") );
    m_intrinsicEfficiency->setText( WString::tr("n/a") );
    m_fractionalAngle->setText( WString::tr("n/a") );
    m_totalEfficiency->setText( WString::tr("n/a") );
  }// try / catch(std::runtime_error& e)
}//updateDetectorCalc()


vector<pair<const SandiaDecay::Element *, float> > GammaXsGui::parseMaterial()
{
  vector<pair<const SandiaDecay::Element *, float> > answer;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  //first look to see if its in the database
  string text = m_materialEdit->text().toUTF8();
  std::shared_ptr<const Material> material;

  const std::shared_ptr<const MaterialDB> matdb = MaterialDB::instance();

  try
  {
    material = matdb->material( text );
  }catch(...){}

  if( !material )
  {
    try
    {
      material = MaterialDB::materialFromChemicalFormula( text, db );

      // Check if this suggestion already exists before adding
      std::shared_ptr<Wt::WAbstractItemModel> mdl = m_materialSuggestion->model();
      const Wt::WString suggName = Wt::WString::fromUTF8( material->name );
      bool alreadyHave = false;
      for( int row = 0; !alreadyHave && (row < mdl->rowCount()); ++row )
        alreadyHave = (Wt::asString( mdl->data( row, 0 ) ) == suggName);
      if( !alreadyHave )
        m_materialSuggestion->addSuggestion( material->name, material->name );
    }catch(...){}
  }//if( !material )

  if( material )
  {
    answer.insert( answer.end(),
                   material->elements.begin(), material->elements.end() );

    for( const Material::NuclideFractionPair &nf : material->nuclides )
    {
      pair<const SandiaDecay::Element *, float> thisone( (const SandiaDecay::Element *)0, 0.0f );
      if( nf.first )
      {
        //XXX - below isnt strictly okay I dont think, since it messes up the
        //      mass fractions a little converting the isotope into the element,
        //      but whatever
        thisone.first = db->element( nf.first->atomicNumber );
        thisone.second = nf.second;
      }//
      if( thisone.first )
        answer.push_back( thisone );
      else
        cerr << "Warning in GammaXsGui::parseMaterial() - unexpected error "
             << "finding the element for a nuclide - wrong answer presented!"
             << endl;
    }//for( const Material::NuclideFractionPair &nf : material->nuclides )

    return answer;
  }//if( material )

  return answer;
}//parseMaterial()



void GammaXsGui::handleAppUrl( std::string query_str )
{
  // Do we want to add an undo/redo step here?
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  if( query_str.size() && query_str[0] == '?' )
    query_str = query_str.substr(1);
  
  const map<string,string> values = AppUtils::query_str_key_values( query_str );
  
  // Check version is appropriate
  const auto ver_iter = values.find( "V" );
  if( (ver_iter == end(values))
     || ((ver_iter->second != "1") && !SpecUtils::istarts_with(ver_iter->second, "1.")) )
    throw runtime_error( "Gamma XS Calc: URI not compatible version" );
  
  // A boilerplate lambda for grapping and validating the fields
  auto getField = [&]( const string &name, const int type, Wt::WLineEdit *edit ){
    const auto i = values.find( name );
    if( i == end(values) )
      throw runtime_error( "Gamma XS Calc: missing field '" + name + "' in URI." );
    string value = i->second;
    SpecUtils::ireplace_all( value, "%23", "#" );
    SpecUtils::ireplace_all( value, "%26", "&" );
    
    // Let the value be empty,
    if( !value.empty() )
    {
      try
      {
        switch( type )
        {
          case 0: break;
          case 1:
            std::stod(value); 
            break;
          case 2: 
            PhysicalUnits::stringToDistance(value); 
            break;
        }
      }catch( std::exception & )
      {
        throw runtime_error( "Gamma XS Calc: invalid field '" + name + "' value '" + value + "' in URI." );
      }// try / catch
    }//if( !value.empty() )
    
    edit->setValueText( WString::fromUTF8(value) );
  };//getField(...)
  
  getField( "E", 1, m_energyEdit );
  getField( "M", 0, m_materialEdit );
  getField( "D", 1, m_density );
  getField( "T", 2, m_thickness );
  getField( "R", 2, m_detectorDistance );
  
  handleMaterialChange();
}//void handleAppUrl( std::string query_str )


string GammaXsGui::encodeStateToUrl() const
{
  // "interspec://gammaxs?e=1001&m=Fe&d=1.2&t=0.1cm&r=1.2m"
  string answer = "V=1";
  
  auto addField = [&answer]( const string &name, Wt::WLineEdit *edit ){
    string val = edit->text().toUTF8();
    SpecUtils::ireplace_all(val, "#", "%23" );
    SpecUtils::ireplace_all(val, "&", "%26" );
    answer += "&" + name + "=" + val;
  };
  
  addField( "E", m_energyEdit );
  addField( "M", m_materialEdit );
  addField( "D", m_density );
  addField( "T", m_thickness );
  addField( "R", m_detectorDistance );
  
  return answer;
}//string GammaXsGui::encodeStateToUrl() const


void GammaXsGui::resetAnserFields()
{
  m_effectiveZ->setText( "---" );
  m_totalAttenuation->setText( "---" );
  m_compton->setText( "---" );
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  m_rayleigh->setText( "---" );
#endif
  m_photoElectric->setText( "---" );
  m_conversion->setText( "---" );
  m_transmissionFraction->setText( "---" );
  m_efficiency->setText( "---" );
  m_intrinsicEfficiency->setText( "---" );
  m_fractionalAngle->setText( "---" );
  m_totalEfficiency->setText( "---" );
  m_transmissionFractionVal = -1.0f;
}//void resetAnserFields()


void GammaXsGui::handleMaterialChange()
{
  const string text = m_materialEdit->text().toUTF8();
  try
  {
    const std::shared_ptr<const Material> material = MaterialDB::instance()->material( text );
    const double density = material->density*PhysicalUnits::cm3/PhysicalUnits::g;
    
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.4g", density );
    m_density->setText( buffer );
  }catch(...)
  {
//    m_density->setText( "" );
  }

  calculateCrossSections();
}//void GammaXsGui::handleMaterialChange()


void GammaXsGui::calculateCrossSections()
{
  checkAndAddUndoRedo();
  
  float energy = -999.0;
  vector<pair<const SandiaDecay::Element *, float> > chemFormula;

  if( m_energyEdit->validate() == ValidationState::Valid )
  {
    try
    {
      energy = static_cast<float>( std::stod( m_energyEdit->text().narrow() ) );
    }catch(...){}
  }//if( m_energyEdit->validate() == ValidationState::Valid )

  chemFormula = parseMaterial();

  if( (energy <= 0.0f) || chemFormula.empty() )
  {
    resetAnserFields();
    return;
  }//if( energy <= 0.0 )

  double atomicMass = 0.0, totalMass = 0.0;
  double comptonMu = 0.0, photoMu = 0.0, pairMu = 0.0;
  double totalMu = 0.0;

#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  double rayleighMu = 0.0;
#endif
  
  energy *= static_cast<float>(PhysicalUnits::keV);

  for( Material::ElementFractionPair &nf : chemFormula )
  {
    const SandiaDecay::Element *el = nf.first;
    const double xsmult = nf.second;
    const int AN = el->atomicNumber;
    
    atomicMass +=  xsmult * AN;
    totalMass += xsmult;
    
    try
    {
      comptonMu  += xsmult * MassAttenuation::massAttenuationCoefficientElement( AN, energy, MassAttenuation::GammaEmProcces::ComptonScatter ) ;
    }catch(exception &e)
    {
      passMessage( WString("gxsg-warn-suspect").arg(e.what()) , 3 );
    }

#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
    try
    {
      rayleighMu += xsmult * MassAttenuation::massAttenuationCoefficientElement( AN, energy, MassAttenuation::GammaEmProcces::RayleighScatter );
    }catch(exception &e)
    {
      passMessage( WString("gxsg-warn-suspect").arg(e.what()) , 3 );
    }
#endif

    try
    {
      photoMu    += xsmult * MassAttenuation::massAttenuationCoefficientElement( AN, energy, MassAttenuation::GammaEmProcces::PhotoElectric );
    }catch(exception &e)
    {
      passMessage( WString("gxsg-warn-suspect").arg(e.what()) , 3 );
    }

    try
    {
      pairMu     += xsmult * MassAttenuation::massAttenuationCoefficientElement( AN, energy, MassAttenuation::GammaEmProcces::PairProduction );
    }catch(exception &e)
    {
      passMessage( WString("gxsg-warn-suspect").arg(e.what()) , 3 );
    }
    
    try
    {
      totalMu += xsmult * MassAttenuation::massAttenuationCoefficientElement( AN, energy );
    }catch(exception &e)
    {
      passMessage( WString("gxsg-warn-suspect").arg(e.what()) , 3 );
    }
  }//for( Material::NuclideFractionPair &nf : chemFormula )

  comptonMu  *= PhysicalUnits::g / PhysicalUnits::cm2;
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  rayleighMu *= PhysicalUnits::g / PhysicalUnits::cm2;
#endif
  photoMu    *= PhysicalUnits::g / PhysicalUnits::cm2;
  pairMu     *= PhysicalUnits::g / PhysicalUnits::cm2;
  totalMu    *= PhysicalUnits::g / PhysicalUnits::cm2;

//  double totalMu = comptonMu /*+ rayleighMu*/ + photoMu + pairMu;
//  double totalMu = MassAttenuation::massAttenuationCoeficient( AN, energy );

  char buffer[32];
  snprintf( buffer, sizeof(buffer), "%.4g", totalMu );
  m_totalAttenuation->setText( buffer );
  
  snprintf( buffer, sizeof(buffer), "%.4g", comptonMu );
  m_compton->setText( buffer );
  
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  snprintf( buffer, sizeof(buffer), "%.4g", rayleighMu );
  m_rayleigh->setText( buffer );
#endif
  
  snprintf( buffer, sizeof(buffer), "%.4g", photoMu );
  m_photoElectric->setText( buffer );
  
  snprintf( buffer, sizeof(buffer), "%.4g", pairMu );
  m_conversion->setText( buffer );

  const double averageAN = atomicMass / totalMass;
  snprintf( buffer, sizeof(buffer), "%.2f", averageAN );
  m_effectiveZ->setText( buffer );
  
  try
  {
    const string densitystr = m_density->text().narrow();
    const string thicknessstr = m_thickness->text().narrow();

    const double density = std::stod(densitystr)
                           * PhysicalUnits::g / PhysicalUnits::cm3;
    const double distance = PhysicalUnits::stringToDistance( thicknessstr );
    totalMu *= PhysicalUnits::cm2 / PhysicalUnits::g;
    const double transmittion = exp( -totalMu * density * distance );
    m_transmissionFractionVal = static_cast<float>( transmittion );
    snprintf( buffer, sizeof(buffer), "%.4g", transmittion );
    m_transmissionFraction->setText( buffer );
  }catch(...)
  {
    m_transmissionFraction->setText( "---" );
    m_transmissionFractionVal = -1.0f;
  }//try / catch
  
  updateDetectorCalc(); //update the Detector values too!
}//void GammaXsGui::calculateCrossSections()


std::array<Wt::WString,4> GammaXsGui::inputs() const
{
  return array<Wt::WString,4>{
    m_energyEdit->text(), m_materialEdit->text(), m_density->text(), m_thickness->text(),
  };
}//std::array<Wt::WString,4> GammaXsWindow::inputs() const


void GammaXsGui::fromInputs( const std::array<Wt::WString,4> &inputs )
{
  m_energyEdit->setText( inputs[0] );
  m_materialEdit->setText( inputs[1] );
  m_density->setText( inputs[2] );
  m_thickness->setText( inputs[3] );
  
  m_prevInputs = inputs;
  
  calculateCrossSections();
}//void fromInputs( const std::array<Wt::WString,4> &inputs );


void GammaXsGui::checkAndAddUndoRedo()
{
  array<Wt::WString,4> current = inputs();
  
  bool match = true;
  for( size_t i = 0; match && (i < m_prevInputs.size()); ++i )
    match = (match && (current[i] == m_prevInputs[i]));
  
  if( match )
    return;
  
  array<Wt::WString,4> prev = std::move(m_prevInputs);
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() )
  {
    auto undo = [prev](){
      InterSpec *viewer = InterSpec::instance();
      GammaXsWindow *xswin = viewer ? viewer->showGammaXsTool() : nullptr;
      GammaXsGui *tool = xswin ? xswin->xstool() : nullptr;
      if( tool )
        tool->fromInputs( prev );
    };
    
    auto redo = [current](){
      InterSpec *viewer = InterSpec::instance();
      GammaXsWindow *xswin = viewer ? viewer->showGammaXsTool() : nullptr;
      GammaXsGui *tool = xswin ? xswin->xstool() : nullptr;
      if( tool )
        tool->fromInputs( current );
    };
    
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Change XS tool value." );
  }//if( can add undo )
  
  m_prevInputs = std::move( current );
}//void checkAndAddUndoRedo();


GammaXsWindow::GammaXsWindow( Wt::WSuggestionPopup *materialSuggestion ,
                              InterSpec* viewer)
  : AuxWindow( WString::tr("window-title-xs-calc"),
              (AuxWindowProperties::PhoneNotFullScreen
               | AuxWindowProperties::SetCloseable
               | AuxWindowProperties::DisableCollapse) ),
  m_tool( nullptr )
{
  rejectWhenEscapePressed( true );

  contents()->setOverflow( Wt::Overflow::Auto, Wt::Orientation::Vertical );

  m_tool = contents()->addNew<GammaXsGui>( materialSuggestion, viewer );
  
  AuxWindow::addHelpInFooter( footer(), "gamma-xs-dialog" );
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = footer()->addNew<WPushButton>();
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( [this](){
    try
    {
      const string url = "interspec://gammaxs/?" + Wt::Utils::urlEncode(m_tool->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("gxsg-tool-state-title"), WString::tr("gxsg-tool-state-text") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  } );
#endif //USE_QR_CODES
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
  int w = viewer->renderedWidth();
  int h = viewer->renderedHeight();
  if( viewer->isPhone() && (w < 100) )
  {
    w = wApp->environment().screenWidth();
    h = wApp->environment().screenHeight();
  }
  
  if( viewer->isPhone() && (w > h) && (w > 100) )
    titleBar()->hide();
  
  show();
  
  if( viewer && (viewer->renderedHeight() > 100) )
  {
    float safeAreas[4] = { 0.0f };
    
#if( IOS )
    InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
    if( app )
    {
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
    }//if( app )
#endif
    //repositionWindow( -32768, static_cast<int>(std::max(3.0f,0.5f*safeAreas[0])) );
    setMaximumSize( WLength::Auto, viewer->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f) );
  }//if( we know the screens size )
  
  resizeToFitOnScreen();
  centerWindowHeavyHanded();
  
  //If mobile take focus away from text field so the keyboard doesnt
  //  automatically show - doesnt always work
  if( viewer->isMobile() )
    closeButton->setFocus();
  
  centerWindow();
}//GammaXsWindow(...) constrctor


GammaXsWindow::~GammaXsWindow()
{
}

GammaXsGui *GammaXsWindow::xstool()
{
  return m_tool;
}
