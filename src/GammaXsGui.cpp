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

#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WDoubleValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>

#include "InterSpec/SpecMeas.h"
#include "InterSpec/DrfSelect.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"

using namespace Wt;
using namespace std;


GammaXsGui::GammaXsGui( MaterialDB *materialDB,
                        Wt::WSuggestionPopup *materialSuggestion,
                        InterSpec* viewer,
                        WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_energyEdit( NULL ),
    m_energyValidator( NULL ),
    m_materialEdit( NULL ),
    m_materialSuggestion( materialSuggestion ),
    m_effectiveZ( NULL ),
    m_totalAttenuation( NULL ),
    m_compton( NULL ),
#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
    m_rayleigh( NULL ),
#endif
    m_photoElectric( NULL ),
    m_conversion( NULL ),
    m_density( NULL ),
    m_distance( NULL ),
    m_transmissionFraction( NULL ),
    m_transmissionFractionVal( -1.0f ),
    m_layout( NULL ),
    m_materialDB( materialDB ),
    m_specViewer( viewer ),
    m_detectorDisplay( NULL ),
    m_efficiency( NULL ),
    m_totalEfficiency( NULL ),
    m_fractionalAngle( NULL )
{
  m_layout = new WGridLayout();
  setLayout( m_layout );

  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_specViewer );

  m_energyEdit = new WLineEdit( "100" );
  m_energyEdit->setAutoComplete( false );
  m_energyValidator = new WDoubleValidator( 1.0, 10000.0, m_energyEdit );
  m_energyEdit->setValidator( m_energyValidator );
  m_energyEdit->addStyleClass( "numberValidator"); //used to detect mobile keyboard
  m_energyEdit->changed().connect( this, &GammaXsGui::calculateCrossSections );
  m_energyEdit->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
  m_energyEdit->focussed().connect( this, &GammaXsGui::calculateCrossSections );
  m_energyEdit->blurred().connect( this, &GammaXsGui::calculateCrossSections );

  int row = 0;
  WLabel *label = new WLabel( "Energy:" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_layout->addWidget( m_energyEdit, row, 1, 1, 1 );
  label = new WLabel( "keV" );
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  ++row;
  
  const char *tooltip =  "You can either enter the name of a pre-defined"
  " material or element (clear form text"
  " and click arrow on right of form to see all predefined options), or you can"
  " enter a chemical formula such as 'H0.112 O0.88' where the numbers are the"
  " density in g/cm3, (or use mass fraction and Density field), of the"
  " preceding element.";
  
  label = new WLabel( "Material/mass-formula" );
  HelpSystem::attachToolTipOn( label, tooltip, showToolTips );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_materialEdit = new WLineEdit( "C0.5H0.2Ni0.3" );
  m_materialEdit->setAutoComplete( false );
  m_materialEdit->changed().connect( this, &GammaXsGui::handleMaterialChange );
  m_materialEdit->enterPressed().connect( this, &GammaXsGui::handleMaterialChange );
//  m_materialEdit->focussed().connect( this, &GammaXsGui::handleMaterialChange );
  m_materialEdit->blurred().connect( this, &GammaXsGui::handleMaterialChange );

  HelpSystem::attachToolTipOn( m_materialEdit, tooltip, showToolTips );
  
  
  m_layout->addWidget( m_materialEdit, row, 1, 1, 2 );
  m_materialSuggestion->forEdit( m_materialEdit,
                  WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  ++row;
  label = new WLabel( "Total att. cross section" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_totalAttenuation = new WText();
  m_layout->addWidget( m_totalAttenuation, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"cm\x00B2/g" );
#else
  label = new WLabel( "cm2/g" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  ++row;
  label = new WLabel( "Compton" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_compton = new WText();
  m_layout->addWidget( m_compton, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"cm\x00B2/g" );
#else
  label = new WLabel( "cm2/g" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  ++row;
  label = new WLabel( "Rayleigh" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_rayleigh = new WText();
  m_layout->addWidget( m_rayleigh, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"cm\x00B2/g" );
#else
  label = new WLabel( "cm2/g" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );
#endif //#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
  
  ++row;
  label = new WLabel( "Photoelectric" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_photoElectric = new WText();
  m_layout->addWidget( m_photoElectric, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"cm\x00B2/g" );
#else
  label = new WLabel( "cm2/g" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  ++row;
  label = new WLabel( "Pair production" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_conversion = new WText();
  m_layout->addWidget( m_conversion, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"cm\x00B2/g" );
#else
  label = new WLabel( "cm2/g" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

  ++row;
  label = new WLabel( "Mass avrg atomic num" );
  m_effectiveZ = new WText( "" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_layout->addWidget( m_effectiveZ, row, 1, 1, 1 );
  
  ++row;
  WText *attText = new WText( "Attenuation (optional):" );
//  attText->setMargin( WLength(1.0,WLength::FontEm), Top );
//  attText->setPadding( WLength(1.0,WLength::FontEm), Top );
  m_layout->addWidget( attText, row, 0, 1, 3, AlignBottom );

  ++row;
  label = new WLabel( "Density:" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_density = new WLineEdit();
  m_density->setAutoComplete( false );
  WDoubleValidator *doubValidator = new WDoubleValidator( m_density );
  m_density->setValidator( doubValidator );
  m_density->addStyleClass( "numberValidator" ); //used to detect mobile keyboard
    
  m_layout->addWidget( m_density, row, 1, 1, 1 );
#ifndef WT_NO_STD_WSTRING
  label = new WLabel( L"g/cm\x00B3" );
#else
  label = new WLabel( "g/cm3" );
#endif
  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );

    
  m_density->changed().connect( this, &GammaXsGui::calculateCrossSections );
  m_density->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
  m_density->focussed().connect( this, &GammaXsGui::calculateCrossSections );
  m_density->blurred().connect( this, &GammaXsGui::calculateCrossSections );

  ++row;
  label = new WLabel( "Thickness:" );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_distance = new WLineEdit( "1 cm" );
  m_distance->setAutoComplete( false );
    
  WRegExpValidator *distValidator
                = new WRegExpValidator( PhysicalUnits::sm_distanceRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  m_distance->setValidator( distValidator );
    
  m_layout->addWidget( m_distance, row, 1, 1, 1 );
//  label = new WLabel( "cm" );
//  m_layout->addWidget( label, row, 2, 1, 1, AlignLeft );
  m_distance->changed().connect( this, &GammaXsGui::calculateCrossSections );
  m_distance->enterPressed().connect( this, &GammaXsGui::calculateCrossSections );
  m_distance->focussed().connect( this, &GammaXsGui::calculateCrossSections );
  m_distance->blurred().connect( this, &GammaXsGui::calculateCrossSections );
  HelpSystem::attachToolTipOn( m_distance, "Thickness of the attenuator.", showToolTips );
  
  ++row;
  label = new WLabel( "Trans. Frac." );
  m_layout->addWidget( label, row, 0, 1, 1, AlignLeft );
  m_transmissionFraction = new WText();
  m_layout->addWidget( m_transmissionFraction, row, 1, 1, 2 );
  
//  HelpSystem::attachToolTipOn( label,"This is the fraction of gammas that will make it through"
//                                     " the specified shielding.", showToolTips );
  
  HelpSystem::attachToolTipOn( m_transmissionFraction,
    "This is the fraction of gammas that will make it through the specified"
    " shielding.", showToolTips );
  m_layout->setColumnStretch( 1, row );

  ++row;
  SpecMeasManager *fileManager = m_specViewer->fileManager();
  m_detectorDisplay = new DetectorDisplay( m_specViewer, fileManager->model() );
  m_detectorDisplay->setMinimumSize( WLength::Auto, WLength(1.5,WLength::FontEm) );
  m_layout->addWidget( m_detectorDisplay, row, 0, 1, 3 );
  
  ++row;
  int detectorCount=0;
  m_detectorLabel[detectorCount] = new WLabel( "Distance" );
  m_layout->addWidget( m_detectorLabel[detectorCount] , row, 0, 1, 1, AlignLeft );
  m_detectorDistance = new WLineEdit("2 cm");
  m_detectorDistance->setAutoComplete( false );
  m_detectorDistance->setValidator( distValidator );
  m_layout->addWidget( m_detectorDistance, row, 1, 1, 2 );
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"This is the distance of the selected detector.", showToolTips );
  HelpSystem::attachToolTipOn( m_detectorDistance,
    "This is the distance from the source center to the detector.",
                              showToolTips );
  
  m_detectorDistance->changed().connect( this, &GammaXsGui::updateDetectorCalc );
  m_detectorDistance->enterPressed().connect( this, &GammaXsGui::updateDetectorCalc );
  m_detectorDistance->focussed().connect( this, &GammaXsGui::updateDetectorCalc );
  m_detectorDistance->blurred().connect( this, &GammaXsGui::updateDetectorCalc );
  
  ++row;
  ++detectorCount;
  m_detectorLabel[detectorCount] = new WLabel( "Intrinsic Efficiency" );
  m_layout->addWidget( m_detectorLabel[detectorCount] , row, 0, 1, 1, AlignLeft );
  m_intrinsicEfficiency = new WText();
  m_layout->addWidget( m_intrinsicEfficiency, row, 1, 1, 2 );
  
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"Intrinsic efficiency (in-peak detection efficiency of gammas striking detector face).", showToolTips );
  HelpSystem::attachToolTipOn( m_intrinsicEfficiency,
    "Intrinsic efficiency (in-peak detection efficiency of gammas striking"
    " detector face).", showToolTips );
  
  ++row;
  ++detectorCount;
  m_detectorLabel[detectorCount] = new WLabel( "Solid Angle Fraction" );
  m_layout->addWidget( m_detectorLabel[detectorCount], row, 0, 1, 1, AlignLeft );
  m_fractionalAngle = new WText();
  m_layout->addWidget( m_fractionalAngle, row, 1, 1, 2 );
  
//  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],"Fractional solid angle of the selected detector at specified distance." , showToolTips );
  HelpSystem::attachToolTipOn( m_fractionalAngle,
    "Fractional solid angle of the selected detector at specified distance.",
                              showToolTips );
  
  
  ++row;
  ++detectorCount;
  m_detectorLabel[detectorCount] = new WLabel( "Detection Efficiency" );
  m_layout->addWidget( m_detectorLabel[detectorCount], row, 0, 1, 1, AlignLeft );
  m_efficiency = new WText();
  m_layout->addWidget( m_efficiency, row, 1, 1, 2 );
  
  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],
    "Intrinsic efficiency times solid angle fraction. E.g. the fraction of"
    " gammas making it out of the shielding that will be detected.",
                              showToolTips );
  HelpSystem::attachToolTipOn( m_efficiency,
    "Intrinsic efficiency times the solid angle fraction. E.g. the fraction of"
    " gammas making it out of the shielding that will be detected. Does not"
    " include attenuation." , showToolTips );
  
  
  ++row;
  ++detectorCount;
  m_detectorLabel[detectorCount] = new WLabel( "Total Efficiency" );
  m_layout->addWidget( m_detectorLabel[detectorCount], row, 0, 1, 1, AlignLeft );
  m_totalEfficiency = new WText();
  m_layout->addWidget( m_totalEfficiency, row, 1, 1, 2 );
  
  HelpSystem::attachToolTipOn( m_detectorLabel[detectorCount],
    "Transmition fraction times detection efficiency." , showToolTips );
  HelpSystem::attachToolTipOn( m_totalEfficiency,
    "Transmition fraction times detection efficiency. E.g. the fraction of"
    " gammas emmitted from the source that will be detected.",
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
  if( m_efficiency->isHidden() == hasDetector )
  {
    const size_t nlabel = (sizeof(m_detectorLabel)/sizeof(*m_detectorLabel));
    for( size_t i = 0; i < nlabel; ++i )
      m_detectorLabel[i]->setHidden( !hasDetector );
    m_efficiency->setHidden( !hasDetector );
    m_intrinsicEfficiency->setHidden( !hasDetector );
    m_fractionalAngle->setHidden( !hasDetector );
    m_detectorDistance->setHidden( !hasDetector );
    m_totalEfficiency->setHidden( !hasDetector );
  }//if( m_efficiency->isHidden() == hasDetector )
  
  if( !hasDetector )
    return;
  
  try
  {
    // efficiency
    const string energystr = m_energyEdit->text().narrow();
    const string diststr = m_detectorDistance->text().narrow();
    const float energy = static_cast<float>( std::stod( energystr ) );
    const double dist = PhysicalUnits::stringToDistance( diststr );
    double val = m_detector->efficiency( energy, dist );
    char buffer[32];
    snprintf( buffer, sizeof(buffer), "%.4g", val );
    m_efficiency->setText( buffer );
    
    if( m_transmissionFractionVal >= 0.0f )
    {
      snprintf( buffer, sizeof(buffer), "%.4g", val*m_transmissionFractionVal );
      m_totalEfficiency->setText( buffer );
    }else
    {
      m_totalEfficiency->setText( "---" );
    }
    
    // absolute efficiency
    val =  m_detector->intrinsicEfficiency( energy );
    snprintf( buffer, sizeof(buffer), "%.4g", val );
    m_intrinsicEfficiency->setText( buffer );
    
    // fractional solid angle
    const float diameter = m_detector->detectorDiameter();
    val =  m_detector->fractionalSolidAngle( diameter, dist );
    snprintf( buffer, sizeof(buffer), "%.4g", val );
    m_fractionalAngle->setText( buffer );
  }catch( std::runtime_error & )
  {
    //could not get all the measurements for the detector
    m_efficiency->setText( "n/a" );
    m_intrinsicEfficiency->setText( "n/a" );
    m_fractionalAngle->setText( "n/a" );
    m_totalEfficiency->setText( "n/a" );
  }// try / catch(std::runtime_error& e)
}//updateDetectorCalc()


vector<pair<const SandiaDecay::Element *, float> > GammaXsGui::parseMaterial()
{
  vector<pair<const SandiaDecay::Element *, float> > answer;

  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

  //first look to see if its in the database
  string text = m_materialEdit->text().toUTF8();
  const Material *material = NULL;

  try
  {
    material = m_materialDB->material( text );
  }catch(...){}

  if( !material )
  {
    try
    {
      material = m_materialDB->parseChemicalFormula( text, db );
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
        cerr << "Waring in GammaXsGui::parseMaterial() - unexpected error "
             << "finding the element for a nuclide - wrong anser presented!"
             << endl;
    }//for( const Material::NuclideFractionPair &nf : material->nuclides )

    return answer;
  }//if( material )

  return answer;
}//parseMaterial()



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
    const Material *material = m_materialDB->material( text );
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
  float energy = -999.0;
  vector<pair<const SandiaDecay::Element *, float> > chemFormula;

  if( m_energyEdit->validate() == WValidator::Valid )
  {
    try
    {
      energy = static_cast<float>( std::stod( m_energyEdit->text().narrow() ) );
    }catch(...){}
  }//if( m_energyEdit->validate() == WValidator::Valid )

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
      comptonMu  += xsmult * MassAttenuation::massAttenuationCoeficient( AN, energy, MassAttenuation::GammaEmProcces::ComptonScatter ) ;
    }catch(exception &e)
    {
      passMessage( string("XS Calculation may be suspect: ")+e.what(), "", 3 );
    }

#if( !USE_SNL_GAMMA_ATTENUATION_VALUES )
    try
    {
      rayleighMu += xsmult * MassAttenuation::massAttenuationCoeficient( AN, energy, MassAttenuation::GammaEmProcces::RayleighScatter );
    }catch(exception &e)
    {
      passMessage( string("XS Calculation may be suspect: ")+e.what(), "", 3 );
    }
#endif

    try
    {
      photoMu    += xsmult * MassAttenuation::massAttenuationCoeficient( AN, energy, MassAttenuation::GammaEmProcces::PhotoElectric );
    }catch(exception &e)
    {
      passMessage( string("XS Calculation may be suspect: ")+e.what(), "", 3 );
    }

    try
    {
      pairMu     += xsmult * MassAttenuation::massAttenuationCoeficient( AN, energy, MassAttenuation::GammaEmProcces::PairProduction );
    }catch(exception &e)
    {
      passMessage( string("XS Calculation may be suspect: ")+e.what(), "", 3 );
    }
    
    try
    {
      totalMu += xsmult * MassAttenuation::massAttenuationCoeficient( AN, energy );
    }catch(exception &e)
    {
      passMessage( string("XS Calculation may be suspect: ")+e.what(), "", 3 );
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
    const string distancestr = m_distance->text().narrow();

    const double density = std::stod(densitystr)
                           * PhysicalUnits::g / PhysicalUnits::cm3;
    const double distance = PhysicalUnits::stringToDistance( distancestr );
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



GammaXsWindow::GammaXsWindow( MaterialDB *materialDB,
                              Wt::WSuggestionPopup *materialSuggestion ,
                              InterSpec* viewer)
  : AuxWindow( "Gamma XS Calc",
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
               | AuxWindowProperties::SetCloseable
               | AuxWindowProperties::DisableCollapse) )
{
  rejectWhenEscapePressed( true );

  new GammaXsGui( materialDB, materialSuggestion, viewer, contents() );
  //gui->setHeight( WLength(100,WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "gamma-xs-dialog" );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( this, &AuxWindow::hide );
  finished().connect( boost::bind( &GammaXsWindow::deleteWindow, this ) );
  
  show();
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(wApp);
  
  if( app && app->isPhone() )
  {
    titleBar()->hide();
  
    if( viewer )
    {
      float safeAreas[4] = { 0.0f };
#if( IOS )
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
#endif
      repositionWindow( -32768, static_cast<int>(std::max(3.0f,0.5f*safeAreas[0])) );
      setMaximumSize( WLength::Auto, viewer->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f) );
    }
  }else
  {
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }//if( isPhone ) / else

  
  //If mobile take focus away from text field so the keyboard doesnt
  //  automatically show - doesnt always work
  if( viewer->isMobile() )
    closeButton->setFocus();
  
  centerWindow();
}//GammaXsWindow(...) constrctor


GammaXsWindow::~GammaXsWindow()
{
}


void GammaXsWindow::deleteWindow( GammaXsWindow *window )
{
  if( window )
    delete window;
}
