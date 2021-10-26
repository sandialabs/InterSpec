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
#include <Wt/WLineEdit>
#include <Wt/WCheckBox>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WApplication>
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
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/ShieldingSourceDisplay.h"


using namespace Wt;
using namespace std;


SourceCheckbox::SourceCheckbox( const SandiaDecay::Nuclide *nuclide,
                               double massFrac, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_useAsSourceCb( NULL ),
    m_massFraction( NULL ),
    m_nuclide( nuclide )
{
//  WLabel *label = NULL;
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

Wt::EventSignal<> &SourceCheckbox::unChecked()
{
  return m_useAsSourceCb->unChecked();
}

Wt::Signal<float> &SourceCheckbox::massFractionChanged()
{
  return m_massFraction->valueChanged();
}



ShieldingSelect::ShieldingSelect( MaterialDB *materialDB,
                                  SourceFitModel *sourceModel,
                                  Wt::WSuggestionPopup *materialSuggest,
                                  bool forFitting,
                                  Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_forFitting( forFitting ),
    m_materialDB( materialDB ),
    m_sourceModel( sourceModel ),
    m_materialSuggest( materialSuggest ),
    m_materialEdit( NULL ),
    m_isGenericMaterial( false ),
    m_materialSummarry( NULL ),
    m_closeIcon( NULL ),
    m_addIcon( NULL ),
    m_thicknessEdit( NULL ),
    m_fitThicknessCB( NULL ),
    m_thicknessDiv( NULL ),
    m_arealDensityEdit( NULL ),
    m_fitArealDensityCB( NULL ),
    m_atomicNumberEdit( NULL ),
    m_fitAtomicNumberCB( NULL ),
    m_fitMassFrac( NULL ),
    m_genericMaterialDiv( NULL ),
    m_asSourceCBs( NULL )
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
    PopupDivMenuItem *item = popup->addMenuItem( "Before this shielding" );
    item->triggered().connect( this, &ShieldingSelect::emitAddBeforeSignal );
    item = popup->addMenuItem( "After this shielding" );
    item->triggered().connect( this, &ShieldingSelect::emitAddAfterSignal );
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


Wt::WLineEdit *ShieldingSelect::thicknessEdit()
{
  return m_thicknessEdit;
}


Wt::WLineEdit *ShieldingSelect::arealDensityEdit()
{
  return m_arealDensityEdit;
}


Wt::WLineEdit *ShieldingSelect::atomicNumberEdit()
{
  return m_atomicNumberEdit;
}



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

  // m_materialEdit->setTextSize( 22 );
  // m_materialEdit->setWidth( 155 );
  
  if( m_materialSuggest )
    m_materialSuggest->forEdit( m_materialEdit,
                   WSuggestionPopup::Editing | WSuggestionPopup::DropDownIcon );

  
  m_materialSummarry = new WText( "", XHTMLText );
  if( m_forFitting )
  {
    materialDivLayout->addWidget( m_materialSummarry, 1, 1, AlignMiddle );
  }else
  {
    //m_materialSummarry->setWidth( WLength(100,WLength::Unit::Pixel) );
    m_materialSummarry->addStyleClass( "MaterialSummary" );
    //materialDivLayout->addWidget( m_materialSummarry, 1, 2 );
  }
  
  materialDivLayout->setColumnStretch(1,1);
  if( m_forFitting )
    setClosableAndAddable( true,  materialDivLayout );

  m_thicknessDiv = new WContainerWidget( this );
  
  WGridLayout * thicknessDivLayout = new WGridLayout();
  m_thicknessDiv->setLayout(thicknessDivLayout);
  thicknessDivLayout->setContentsMargins(3,3,3,3);
  WLabel *label = new WLabel( "Thickness" );
  thicknessDivLayout->addWidget(label,0,0,AlignMiddle);
  m_thicknessEdit = new WLineEdit( "1.0 cm" );
  m_thicknessEdit->setAutoComplete( false );
  
  label->setBuddy( m_thicknessEdit );
  
  WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_thicknessEdit->setValidator( validator );
  
  if( m_forFitting )
  {
    m_thicknessEdit->setWidth( 150 );
    thicknessDivLayout->addWidget(m_thicknessEdit,0,1,AlignMiddle);
    
    m_fitThicknessCB = new WCheckBox( "Fit" );
    m_fitThicknessCB->setChecked( true );
    thicknessDivLayout->addWidget(m_fitThicknessCB,0,2,AlignMiddle | AlignRight);
    thicknessDivLayout->setColumnStretch(3,1);
  }else
  {
    thicknessDivLayout->addWidget( m_thicknessEdit, 0, 1 );
    thicknessDivLayout->addWidget( m_materialSummarry, 0, 2, AlignMiddle );
    thicknessDivLayout->setColumnStretch( 1, 1 );
  }//if( m_forFitting ) / else
  
  
  const WBorder anAdBorder( WBorder::Solid, WBorder::None, black );

  m_genericMaterialDiv = new WContainerWidget( this );
  WGridLayout *genericMatLayout = new WGridLayout( m_genericMaterialDiv );
  genericMatLayout->setContentsMargins(3,3,3,3);
//  m_genericMaterialDiv->decorationStyle().setBorder( anAdBorder, Wt::All );

//  WContainerWidget *anDiv = new WContainerWidget( m_genericMaterialDiv );
//  WContainerWidget *adDiv = new WContainerWidget( m_genericMaterialDiv );
  
//  anDiv->decorationStyle().setBorder( anAdBorder, Wt::All );
//  adDiv->decorationStyle().setBorder( anAdBorder, Wt::All );
  
  
//  WContainerWidget *anAdDiv = new WContainerWidget( m_genericMaterialDiv );
//  WGridLayout *anAdDivLayout = new WGridLayout( anAdDiv );
//  anAdDivLayout->setContentsMargins(3,3,3,3);
//  anAdDivLayout->setColumnStretch( 1, 1 );
  
  label = new WLabel( "AD" );
  label->setAttributeValue( "style", "padding-left: 1em;" );
//  label = new WLabel( "Areal Density" );
  label->setToolTip( "Areal Density of the shielding in g/cm2" );
//  anAdDivLayout->addWidget(label,1,0,AlignMiddle);
  genericMatLayout->addWidget( label, 0, 2+m_forFitting, AlignMiddle );
  
  m_arealDensityEdit = new WLineEdit();
  m_arealDensityEdit->setToolTip( "Areal Density of the shielding in g/cm2" );
  if( m_forFitting )
    m_arealDensityEdit->setText( "15.0" );
  else
    m_arealDensityEdit->setEmptyText( "Areal Density" );
  
//  anAdDivLayout->addWidget(m_arealDensityEdit,1,1,AlignMiddle);
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
//  anAdDivLayout->addWidget(label,1,2,AlignMiddle);
  genericMatLayout->addWidget( label, 0, 4+m_forFitting, AlignMiddle );
  
  m_arealDensityEdit->addStyleClass( "numberValidator" ); //used to detect mobile keyboard
  
  if( m_forFitting )
  {
    m_fitArealDensityCB = new WCheckBox( "Fit" );
    m_fitArealDensityCB->setChecked( true );
//    anAdDivLayout->addWidget(m_fitArealDensityCB,1,3,AlignRight | AlignMiddle);
    genericMatLayout->addWidget( m_fitArealDensityCB, 0, 6, AlignMiddle );
  }

//  label = new WLabel( "Atomic Num." );
//  anAdDivLayout->addWidget(label,0,0,AlignMiddle);
  label = new WLabel( "AN" );
  label->setToolTip( "Atomic Number of the shielding" );
  genericMatLayout->addWidget( label, 0, 0, AlignMiddle );
  
  m_atomicNumberEdit = new WLineEdit();
  m_atomicNumberEdit->setToolTip( "Atomic Number of the shielding" );
  if( m_forFitting )
    m_atomicNumberEdit->setText( "15.0" );
  else
    m_atomicNumberEdit->setEmptyText( "Atomic Number" );

//  anAdDivLayout->addWidget(m_atomicNumberEdit,0,1,AlignMiddle);
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
//    anAdDivLayout->addWidget(m_fitAtomicNumberCB,0,3,AlignRight |AlignMiddle);
    genericMatLayout->addWidget( m_fitAtomicNumberCB, 0, 2, AlignMiddle );
    
    m_asSourceCBs = new WContainerWidget( this );
    m_asSourceCBs->addStyleClass( "ShieldingAsSourceCBDiv" );
    WLabel *label = new WLabel( "Source for:", m_asSourceCBs );
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

  
  handleMaterialChange();
  
  m_atomicNumberEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_atomicNumberEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_atomicNumberEdit->blurred().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  
  m_arealDensityEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_arealDensityEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_arealDensityEdit->blurred().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  
  m_materialEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_materialEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_materialEdit->blurred().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );

  m_thicknessEdit->changed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_thicknessEdit->enterPressed().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );
  m_thicknessEdit->blurred().connect( boost::bind( &ShieldingSelect::handleMaterialChange, this ) );

  //From a very brief experiment, it looks like the below JS would remove the
  //  uncertainty text, however it appears when the server pushes the value of
  //  the edit to the client, the 'value' tag of the element isnt updated... I
  //  have no clue.
  //  const char *focusjs = "function(s,e){try{"
  //  "s.value = s.value.replace('\\(.*\\)','').value.replace('  ',' ');"
  //  "}catch(e){}}";
  //  m_thicknessEdit->changed().connect( focusjs );
  m_thicknessEdit->changed().connect( this, &ShieldingSelect::removeUncertFromThickness );
  m_thicknessEdit->enterPressed().connect( this, &ShieldingSelect::removeUncertFromThickness );
  m_thicknessEdit->blurred().connect( this, &ShieldingSelect::removeUncertFromThickness );
  
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
          removingIsotopeAsSource().emit( cb->isotope() );
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


Wt::Signal<ShieldingSelect *,const SandiaDecay::Nuclide *> &ShieldingSelect::activityFromThicknessNeedUpdating()
{
  return m_activityFromThicknessNeedUpdating;
}


Wt::Signal<const SandiaDecay::Nuclide *> &ShieldingSelect::addingIsotopeAsSource()
{
  return m_addingIsotopeAsSource;
}


Wt::Signal<const SandiaDecay::Nuclide *> &ShieldingSelect::removingIsotopeAsSource()
{
  return m_removingIsotopeAsSource;
}


bool ShieldingSelect::isGenericMaterial() const
{
  return m_isGenericMaterial;
}//bool isGenericMaterial() const


double ShieldingSelect::thickness() const
{
  if( isGenericMaterial() )
    throw std::runtime_error( "ShieldingSelect::thickness() can not be called "
                              "for a generic material" );
  
  const string text = m_thicknessEdit->text().toUTF8();

  try
  {
    return PhysicalUnits::stringToDistance( text );
  }catch( std::exception & )
  {
  }
  
  throw runtime_error( "Error converting '" + text + "' to a distance" );
  return 0.0;
}//double thickness() const


double ShieldingSelect::atomicNumber() const
{
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
  if( isGenericMaterial() || !m_forFitting )
    throw std::runtime_error( "ShieldingSelect::fitThickness() can not be "
                              "called for a generic material" );
  return m_fitThicknessCB->isChecked();
}//bool fitThickness() const


//fitAtomicNumber():
//  throws std::runtime_error if not a GenericMaterial
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
  //See if 'text' is the name of a meterial in the database
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
  }catch(...){}

  return 0;
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
#if( !ALLOW_MASS_FRACTION_FIT )
  assert(0);
#endif
  
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
#if( ALLOW_MASS_FRACTION_FIT )
  updateIfMassFractionCanFit();
#endif
  m_addingIsotopeAsSource.emit( nuc );
}//void isotopeCheckedCallback( const std::string symbol )



void ShieldingSelect::isotopeUnCheckedCallback( const SandiaDecay::Nuclide *iso )
{
#if( ALLOW_MASS_FRACTION_FIT )
  updateIfMassFractionCanFit();
#endif
  m_removingIsotopeAsSource.emit( iso );
}//void isotopeUnCheckedCallback( const std::string symbol )



void ShieldingSelect::uncheckSourceIsotopeCheckBox( const SandiaDecay::Nuclide *iso )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();

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
//      removingIsotopeAsSource().emit( symbol );
      cb->setUseAsSource( false );
    }//if( cb->isChecked() )
  }//for( WWidget *child : children )
  
#if( !ALLOW_MASS_FRACTION_FIT )
  updateIfMassFractionCanFit();
#endif
}//void uncheckSourceIsotopeCheckBox( const std::string &symol )


void ShieldingSelect::removeSourceIsotopeCheckBox( const SandiaDecay::Nuclide *nuc )
{
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
  
#if( !ALLOW_MASS_FRACTION_FIT )
  updateIfMassFractionCanFit();
#endif
  
  //call updateMassFractionDisplays() to update the "Assuming XX% other U isos"
  updateMassFractionDisplays( m_currentMaterial );
}//void removeSourceIsotopeCheckBox( const std::string &symbol )


vector< ShieldingSelect::NucMasFrac > ShieldingSelect::sourceNuclideMassFractions()
{
  vector< ShieldingSelect::NucMasFrac > answer;
  
  const vector<const SandiaDecay::Nuclide *> nucs = nuclidesToUseAsSources();
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


vector<const SandiaDecay::Nuclide *> ShieldingSelect::nuclidesToUseAsSources()
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
}//std::vector<const SandiaDecay::Nuclide *> nuclidesToUseAsSources() const


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



void ShieldingSelect::addSourceIsotopeCheckBox( const SandiaDecay::Nuclide *iso )
{
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
}//void addSourceIsotopeCheckBox( const std::string &symol )


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
  
#if( ALLOW_MASS_FRACTION_FIT )
  updateIfMassFractionCanFit();
#endif
}//void updateMassFractionDisplays()


void ShieldingSelect::removeUncertFromThickness()
{
  cerr << "removeUncertFromThickness()" << endl;
  //Get rid of the uncertainty text if the value is edited; this currently
  //  will cause the uncertainty to be removed
  string thickstr = m_thicknessEdit->text().toUTF8();
  
  SpecUtils::trim( thickstr );
  
  if( thickstr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
  {
    thickstr += " cm";
    m_thicknessEdit->setText( thickstr );
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
  
  m_thicknessEdit->setText( thickstr );
}//void ShieldingSelect::removeUncertFromThickness()


void ShieldingSelect::handleToggleGeneric()
{
  m_isGenericMaterial = (!m_isGenericMaterial);
  
  
  if( m_isGenericMaterial )
  {
    const string oldmaterial = m_materialEdit->text().toUTF8();
    
    //See if we can convert the current material into AN, AD
    m_thicknessDiv->hide();
    m_genericMaterialDiv->show();
    m_materialSummarry->setText( "" );
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
    m_thicknessDiv->show();
    m_genericMaterialDiv->hide();
    
    string aNstr = SpecUtils::trim_copy( m_atomicNumberEdit->text().toUTF8() );
    string aDstr = SpecUtils::trim_copy( m_arealDensityEdit->text().toUTF8() );
    
    if( aNstr.empty() || aDstr.empty()
        || std::atof(aDstr.c_str())<=0.0
        || std::atof(aNstr.c_str())<=0.0 )
    {
      m_materialEdit->setText( "" );
      m_thicknessEdit->setText( "1.0 cm" );
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
          m_thicknessEdit->setText( diststr );
        }else
        {
          m_thicknessEdit->setText( "0 cm" );
        }
      }else
      {
        m_materialEdit->setText( "" );
        m_thicknessEdit->setText( "1 cm" );
      }
    }//if( input is empty ) / else
  }//if( m_isGenericMaterial ) / else
  
  
  handleMaterialChange();
}//void ShieldingSelect::handleToggleGeneric()


void ShieldingSelect::handleMaterialChange()
{
  typedef pair<const SandiaDecay::Element *,float> ElementFrac;
  typedef pair<const SandiaDecay::Nuclide *,float> NuclideFrac;

  std::shared_ptr<Material> newMaterial;
  std::shared_ptr<Material> previousMaterial = m_currentMaterial;
  
  if( m_isGenericMaterial )
  {
    m_thicknessDiv->hide();
    m_genericMaterialDiv->show();
    m_materialSummarry->setText( "" );
    m_materialEdit->setText( "Generic" );
    m_materialEdit->disable();
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/atom_black.png") );
  }else
  {
    m_toggleImage->setImageLink( Wt::WLink("InterSpec_resources/images/shield.png") );
    m_materialEdit->enable();
    m_thicknessDiv->show();
    m_genericMaterialDiv->hide();
    
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
    
    m_materialSummarry->setText( summary );
  }//if( generic material ) / else
  
  
  //Now we need to update the activites for any isotopes that are
  if( !!newMaterial && m_asSourceCBs && (previousMaterial==newMaterial) )
  {
    updateMassFractionDisplays( newMaterial );

    for( ElementToNuclideMap::value_type &vt : m_sourceIsotopes )
    {
      const vector<WWidget *> children = vt.second->children();
      for( WWidget *child : children )
      {
        SourceCheckbox *cb = dynamic_cast<SourceCheckbox *>( child );
        if( cb && cb->useAsSource() )
          m_activityFromThicknessNeedUpdating.emit( this, cb->isotope() );
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
          removingIsotopeAsSource().emit( cb->isotope() );
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
        addSourceIsotopeCheckBox( nuc );
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



