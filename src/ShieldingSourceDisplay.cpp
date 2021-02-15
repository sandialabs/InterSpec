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

#include <boost/scope_exit.hpp>

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnUserParameters.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameterState.h"


#include <Wt/WText>
#include <Wt/WTime>
#include <Wt/WImage>
#include <Wt/WTable>
#include <Wt/WLabel>
#include <Wt/WAnchor>
#include <Wt/WBorder>
#include <Wt/WServer>
#include <Wt/WPainter>
#include <Wt/WTextArea>
#include <Wt/WRectArea>
#include <Wt/WResource>
#include <Wt/WSvgImage>
#include <Wt/WTableCell>
#include <Wt/WTabWidget>
#include <Wt/WIOService>
#include <Wt/WGridLayout>
#include <Wt/WPushButton>
#include <Wt/WJavaScript>
#include <Wt/WFileUpload>
#include <Wt/WSplitButton>
#include <Wt/Http/Response>
#include <Wt/WSelectionBox>
#include <Wt/WStandardItem>
#include <Wt/WItemDelegate>
#include <Wt/WDoubleSpinBox>
#include <Wt/Dbo/QueryModel>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>
#include <Wt/WDoubleValidator>
#include <Wt/WStandardItemModel>
#include <Wt/WCssDecorationStyle>
#include <Wt/Chart/WChartPalette>
#include <Wt/Chart/WCartesianChart>

#include "InterSpec/PeakDef.h"
#include "InterSpec/SpecMeas.h"
#include "InterSpec/PopupDiv.h"
#include "InterSpec/PeakModel.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/MaterialDB.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/DetectorEdit.h"
#include "InterSpec/InterSpecUser.h"
#include "InterSpec/DataBaseUtils.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/PhysicalUnits.h"
#include "SandiaDecay/SandiaDecay.h"
#include "InterSpec/SpecMeasManager.h"
#include "InterSpec/RowStretchTreeView.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/GammaInteractionCalc.h"
#include "InterSpec/ShieldingSourceDisplay.h"

using namespace Wt;
using namespace std;

typedef std::shared_ptr<const PeakDef> PeakShrdPtr;

const int ShieldingSelect::sm_xmlSerializationVersion = 0;
const int ShieldingSourceDisplay::sm_xmlSerializationVersion = 0;

namespace
{
  const std::string ns_no_uncert_info_txt = "Perform model fit to update and get uncertainties.";
}


WT_DECLARE_WT_MEMBER
(ShowChi2Info, Wt::JavaScriptFunction, "ShowChi2Info",
function(id,info)
{
  $("<div id=\"" + id + "inf\" class=\"PointInfo\"></div>").html(info).appendTo($('#'+id));
}
);


namespace
{
  class StringDownloadResource : public Wt::WResource
  {
    ShieldingSourceDisplay *m_display;
  public:
    StringDownloadResource( ShieldingSourceDisplay *parent )
      : WResource( parent ), m_display( parent )
    {}
    
    virtual ~StringDownloadResource()
    {
      beingDeleted();
    }
    
    virtual void handleRequest( const Wt::Http::Request &request,
                               Wt::Http::Response &response )
    {
      suggestFileName( "shielding_source_fit_model.xml", WResource::Attachment );
      response.setMimeType( "application/xml" );
      if( !m_display )
        return;
      rapidxml::xml_document<char> doc;
      m_display->serialize( &doc );
      rapidxml::print( response.out(), doc, 0 );
    }
  };//class PeakCsvResource : public Wt::WResource
  
}//namespace

SourceCheckbox::SourceCheckbox( const SandiaDecay::Nuclide *nuclide,
                               double massFrac, Wt::WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_useAsSourceCb( NULL ),
    m_massFraction( NULL ),
    m_nuclide( nuclide )
{
//  WLabel *label = NULL;
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
  
  m_materialEdit->setText( WString::fromUTF8(name) );
  m_thicknessEdit->setText( WString::fromUTF8(thickness) );
  
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


const Wt::WLineEdit *ShieldingSelect::materialEdit() const
{
  return m_materialEdit;
}

const Wt::WLineEdit *ShieldingSelect::thicknessEdit() const
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
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSourceDisplay.css" );
  
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
    InterSpecApp *app = dynamic_cast<InterSpecApp *>( wApp );
    const bool isMobile = (app && app->isMobile());
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
    return std::shared_ptr<Material>();

  const string text = SpecUtils::trim_copy( m_materialEdit->text().toUTF8() );
  if( m_currentMaterial && text==m_currentMaterialDescrip )
    return m_currentMaterial;

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
  
  if( !m_fitMassFrac
      || (m_fitMassFrac->isVisible() && !m_fitMassFrac->isChecked()) )
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
  
  if( nchecked < 2 )
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
      
      if( fabs(cb->massFraction()-massFrac) > 0.00001 )
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


SourceFitModel::SourceFitModel( PeakModel *peakModel,
                                const bool sameAgeIsotopes,
                                Wt::WObject *parent )
  : WAbstractItemModel( parent ),
    m_sortOrder( Wt::AscendingOrder ),
    m_sortColumn( kIsotope ),
    m_peakModel( peakModel ),
    m_sameAgeForIsotopes( sameAgeIsotopes )
{
  auto interspec = InterSpec::instance();
  if( !interspec )
  {
    m_displayCurries = true;
  }else
  {
    m_displayCurries = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", interspec );
    //auto callback = wApp->bind( boost::bind( &SourceFitModel::displayUnitsChanged, this) );
    //InterSpecUser::addCallbackWhenChanged( interspec->m_user, "DisplayBecquerel", callback, interspec );
    InterSpecUser::addCallbackWhenChanged( "DisplayBecquerel", this, &SourceFitModel::displayUnitsChanged );
  }//if( !interspec ) / else
  
  peakModel->rowsAboutToBeRemoved().connect( this, &SourceFitModel::peakModelRowsRemovedCallback );
  peakModel->rowsInserted().connect( this, &SourceFitModel::peakModelRowsInsertedCallback );
  peakModel->dataChanged().connect( this, &SourceFitModel::peakModelDataChangedCallback );
  peakModel->modelReset().connect( this, &SourceFitModel::peakModelResetCallback );
  repopulateIsotopes();
}//SourceFitModel(...)


SourceFitModel::~SourceFitModel()
{
  //nothing to do here
}//~SourceFitModel()


void SourceFitModel::displayUnitsChanged( boost::any value )
{
  //cout << "in SourceFitModel::displayUnitsChanged" << endl;
  try
  {
    //const bool useBq = InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
    const bool useBq = boost::any_cast<bool>(value);
    
    if( useBq == m_displayCurries )
    {
      m_displayCurries = !useBq;
      const int nrow = rowCount();
      if( nrow > 0 )
        dataChanged().emit( createIndex(0,0,nullptr), createIndex(nrow-1,kNumColumns-1,nullptr) );
    }
  }catch( std::exception &e )
  {
    //Shouldnt ever happen, but print a little something out JIC
    cerr << "SourceFitModel::displayUnitsChanged: Failed to convert boost any: " << e.what() << endl;
  }
  
  //cout << "m_displayCurries is now: " << m_displayCurries << endl;
}//void SourceFitModel::displayUnitsChanged( boost::any value )


int SourceFitModel::numNuclides() const
{
  return static_cast<int>( m_nuclides.size() );
}//int numNuclides() const


int SourceFitModel::row( const SandiaDecay::Nuclide *nuclide ) const
{
  if( !nuclide )
    return -1;

  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int i = 0; i < nrow; ++i )
    if( m_nuclides[i].nuclide == nuclide )
      return i;
  return -1;
}//int row( const SandiaDecay::Nuclide *nuclide ) const;


int SourceFitModel::nuclideIndex( const SandiaDecay::Nuclide *nuc ) const
{
  for( size_t i = 0; i < m_nuclides.size(); ++i )
    if( m_nuclides[i].nuclide == nuc )
      return static_cast<int>(i);
  throw std::runtime_error( "SourceFitModel::nuclideIndex with invalid nuclide (" + (nuc ? nuc->symbol : string("null")) + ")" );
  return -1;
}//int nuclideNumber( const SandiaDecay::Nuclide *nuc ) const


const SandiaDecay::Nuclide *SourceFitModel::nuclide( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].nuclide;
}//const SandiaDecay::Nuclide *nuclide( int nuc ) const


double SourceFitModel::activity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].activity * GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
}//double activity( int nuc ) const

double SourceFitModel::activityUncert( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].activityUncertainty * GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
}//double activityUncert( int nuc ) const


bool SourceFitModel::fitActivity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].fitActivity;
}//bool fitActivity( int nuc ) const


double SourceFitModel::age( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].age;
}//double age( int nuc ) const


double SourceFitModel::ageUncert( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].ageUncertainty;
}//double ageUncert( int nuc ) const


bool SourceFitModel::fitAge( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  const auto &nucibj = m_nuclides[nuc];
  
  return ((nucibj.numProgenyPeaksSelected > 1 || nucibj.ageDefiningNuc) && nucibj.fitAge);
}//bool fitAge( int nuc ) const

#if( INCLUDE_ANALYSIS_TEST_SUITE )
boost::optional<double> SourceFitModel::truthActivity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthActivity;
}

boost::optional<double> SourceFitModel::truthActivityTolerance( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthActivityTolerance;
}

boost::optional<double> SourceFitModel::truthAge( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthAge;
}

boost::optional<double> SourceFitModel::truthAgeTolerance( int nuc ) const
{
  if( (nuc < 0) || (nuc >= static_cast<int>(m_nuclides.size())) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].truthAgeTolerance;
}
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )


const SandiaDecay::Nuclide *SourceFitModel::ageDefiningNuclide(
                                  const SandiaDecay::Nuclide *dependantNuc ) const
{
  const int nuc = nuclideIndex( dependantNuc );
  if( (nuc < 0) || (nuc >= static_cast<int>(m_nuclides.size())) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  
  if( m_nuclides[nuc].ageDefiningNuc )
    return m_nuclides[nuc].ageDefiningNuc;
  return m_nuclides[nuc].nuclide;
}//bool ageDefiningNuclide( int nuc ) const


bool SourceFitModel::shieldingDeterminedActivity( int nuc ) const
{
  if( nuc<0 || nuc>=static_cast<int>(m_nuclides.size()) )
    throw std::runtime_error( "SourceFitModel: called with invalid index" );
  return m_nuclides[nuc].shieldingIsSource;
}//bool shieldingDeterminedActivity( int nuc ) const


void SourceFitModel::setSharredAgeNuclide( const SandiaDecay::Nuclide *dependantNuc,
                                           const SandiaDecay::Nuclide *definingNuc )
{
  const int row = nuclideIndex( dependantNuc );
  if( row < 0 )
    return;
  
  if( definingNuc )
  {
    const int definingRow = nuclideIndex( definingNuc );
    if( definingRow < 0 )
      throw runtime_error( "SourceFitModel::setSharredAgeNuclide: defining"
                           " nuclide must also be in the model" );
  }//if( definingNuc )
 
  if( definingNuc == dependantNuc )
    definingNuc = nullptr;
  
  IsoFitStruct &iso = m_nuclides[row];
  if( iso.ageDefiningNuc == definingNuc )
    return;
  
  if( definingNuc && dependantNuc && (definingNuc->atomicNumber != dependantNuc->atomicNumber) )
    throw runtime_error( "SourceFitModel::setSharredAgeNuclide: dependant and"
                         " defining nuclides must have same atomic number" );
  
  iso.ageDefiningNuc = definingNuc;
  dataChanged().emit( createIndex(row,kAge,(void *)0),
                      createIndex(row,kAgeUncertainty,(void *)0) );
}//void makeAgeFitable( const SandiaDecay::Nuclide *nuc, bool fit )


void SourceFitModel::makeActivityEditable( const SandiaDecay::Nuclide *nuc )
{
  const int nrows = static_cast<int>( m_nuclides.size() );

  for( int row = 0; row < nrows; ++row )
  {
    IsoFitStruct &iso = m_nuclides[row];

    if( iso.nuclide && (iso.nuclide==nuc) )
    {
      if( !iso.shieldingIsSource )
        continue;
      iso.shieldingIsSource = false;

      dataChanged().emit( createIndex(row,kActivity,(void *)0),
                          createIndex(row,kFitActivity,(void *)0) );
    }//if( this is the nuclide we want )
  }//for( const IsoFitStruct &iso : m_nuclides )
}//void makeActivityEditable( const std::string &symbol )


void SourceFitModel::makeActivityNonEditable( const SandiaDecay::Nuclide *nuc )
{
  const int nrows = static_cast<int>( m_nuclides.size() );

  for( int row = 0; row < nrows; ++row )
  {
    IsoFitStruct &iso = m_nuclides[row];

    if( iso.nuclide && (iso.nuclide==nuc) )
    {
      if( iso.shieldingIsSource )
        continue;
      iso.shieldingIsSource = true;

      dataChanged().emit( createIndex(row,kActivity,(void *)0),
                          createIndex(row,kFitActivity,(void *)0) );
    }//if( this is the nuclide we want )
  }//for( const IsoFitStruct &iso : m_nuclides )
}//void makeActivityNonEditable( const std::string &symbol )


void SourceFitModel::setUseSameAgeForIsotopes( bool useSame )
{
  if( m_sameAgeForIsotopes == useSame )
    return;

  m_sameAgeForIsotopes = useSame;
  
  if( !m_sameAgeForIsotopes )
  {
    for( size_t i = 0; i < m_nuclides.size(); ++i )
      setSharredAgeNuclide( m_nuclides[i].nuclide, NULL );
    return;
  }//if( !m_sameAgeForIsotopes )
  
  //construct map from element number to isotopes being fit for
  //if there is more than one isotope, choose the youngest age to assign to
  //all of the isotopes.  Do the assignemnet and update the chart
  typedef map< int, vector<const SandiaDecay::Nuclide *> > ElToNucMap_t;
  ElToNucMap_t elToNucMap;
  
  for( const IsoFitStruct &n : m_nuclides )
  {
    if( n.nuclide )
      elToNucMap[n.nuclide->atomicNumber].push_back( n.nuclide );
  }
  
  for( const ElToNucMap_t::value_type &vt : elToNucMap )
  {
    const vector<const SandiaDecay::Nuclide *> &nucs = vt.second;
    if( nucs.size() < 2 )
      continue;
    
    bool fitAgeWanted = false;
    size_t minIndex = 0;
    double minAge = std::numeric_limits<double>::infinity();
    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const SandiaDecay::Nuclide *nuc = nucs[i];
      const int ind = nuclideIndex( nuc );
      const double thisage = age( ind );
      if( thisage < minAge )
      {
        minIndex = i;
        minAge = thisage;
      }//if( age < minAge )
      fitAgeWanted = (fitAgeWanted || fitAge(ind));
    }//for( const SandiaDecay::Nuclide *nuc : nucs )
    
    for( size_t i = 0; i < nucs.size(); ++i )
    {
      const SandiaDecay::Nuclide *nuc = nucs[i];
      
      if( i == minIndex )
      {
        WModelIndex ind = index( nuc, SourceFitModel::kAge );
        setData( ind, fitAgeWanted );
        ind = index( nuc, SourceFitModel::kAgeUncertainty );
        setData( ind, boost::any() );
        setSharredAgeNuclide( nuc, NULL );
      }else
      {
        setSharredAgeNuclide( nuc, nucs[minIndex] );
      }//if( i == 0 )
    }//for( const SandiaDecay::Nuclide *nuc : nucs )
  }//for( const ElToNucMap_t::value_type &vt : elToNucMap )
}//void setUseSameAgeForIsotopes( bool useSame )


void SourceFitModel::insertPeak( const PeakShrdPtr peak )
{
  if( !peak )
  {  //probably wont ever happen
    cerr << "SourceFitModel::insertPeak(...)\n\tNot a valid row being added" << endl;
    return;
  }//if( !peak )

  const SandiaDecay::Nuclide *nuclide = peak->parentNuclide();
  if( !nuclide )
    return;
  
  if( !peak->useForShieldingSourceFit() )
    return;

  for( const IsoFitStruct &iso : m_nuclides )
    if( iso.nuclide == peak->parentNuclide() )
      return;

  IsoFitStruct newIso;
  newIso.activity = 0.001*PhysicalUnits::curie / GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
  newIso.fitActivity = true;
  newIso.nuclide = peak->parentNuclide();
  newIso.age = PeakDef::defaultDecayTime( newIso.nuclide );
  newIso.fitAge = false;
  newIso.ageIsFittable = !PeakDef::ageFitNotAllowed( newIso.nuclide );
  
  std::map<const SandiaDecay::Nuclide *, IsoFitStruct>::iterator oldval
                                    = m_previousResults.find( newIso.nuclide );
  if( oldval != m_previousResults.end() )
  {
    const IsoFitStruct &oldIso = oldval->second;
    newIso.activity = oldIso.activity;
    newIso.fitActivity = oldIso.fitActivity;
    newIso.age = oldIso.age;
    newIso.fitAge = oldIso.fitAge;
    newIso.ageIsFittable = oldIso.ageIsFittable;
    
    m_previousResults.erase( oldval );
  }//if( m_previousResults.count(newIso.nuclide) )
  
  if( m_previousResults.size() > 100 )
    m_previousResults.clear();
  
  if( m_sameAgeForIsotopes )
  {
    vector<size_t> thisElementIndexs;
    for( size_t i = 0; i < m_nuclides.size(); ++i )
    {
      const IsoFitStruct &iso = m_nuclides[i];
      if( iso.nuclide->atomicNumber == newIso.nuclide->atomicNumber )
      {
        thisElementIndexs.push_back( i );
        if( !iso.ageDefiningNuc )
          newIso.ageDefiningNuc = iso.nuclide;
      }
    }//for( const IsoFitStruct &iso : m_nuclides )
    
    if( !newIso.ageDefiningNuc && !thisElementIndexs.empty() )
    {
      const IsoFitStruct &previso = m_nuclides[thisElementIndexs[0]];
      newIso.ageDefiningNuc = previso.nuclide;
      //There is a slight hickup with emitting a datachanged and then inserting
      //  a row; if there wasnt this hickup, it would be nice to use the below
//      if( newIso.age < previso.age )
//        setSharredAgeNuclide( previso.nuclide, newIso.nuclide );
//      else
//        newIso.ageDefiningNuc = previso.nuclide;
    }//if( !newIso.ageDefiningNuc && !thisElementIndexs.empty() )
  }//if( m_sameAgeForIsotopes )
  
  
  
  const int npeaks = m_peakModel->rowCount();
  set<const SandiaDecay::Nuclide *> progeny;
  for( int peakn = 0; peakn < npeaks; ++peakn )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );
    if( peak && (peak->parentNuclide() == newIso.nuclide) && peak->useForShieldingSourceFit() )
    {
      const SandiaDecay::Transition * const trans = peak->nuclearTransition();
      if( trans && trans->parent )
        progeny.insert( trans->parent );
    }
  }//for( loop over peaks)
  newIso.numProgenyPeaksSelected = progeny.size();
  
  std::vector<IsoFitStruct>::iterator pos;
  pos = std::lower_bound( m_nuclides.begin(), m_nuclides.end(), newIso,
                          boost::bind( &SourceFitModel::compare, _1, _2, m_sortColumn, m_sortOrder ) );

  const int row = static_cast<int>( pos - m_nuclides.begin() );
  beginInsertRows( WModelIndex(), row, row );
  m_nuclides.insert( pos, newIso );
  endInsertRows();
}//void SourceFitModel::insertPeak( const PeakShrdPtr peak )


void SourceFitModel::peakModelRowsInsertedCallback( Wt::WModelIndex /*parent*/,
                                                     int firstRow, int lastRow )
{
//  if( firstRow != lastRow )
//  {
//    WStringStream msg;
//    msg << "SourceFitModel::peakModelRowsInsertedCallback(...)"
//            " first row should equal last row: "
//         << firstRow << " != " << lastRow;
//    cerr << msg.str() << endl;
//    throw std::runtime_error( msg.str() );
//  }//if( firstRow != lastRow )

  for( int row = firstRow; row <= lastRow; ++row )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );
    if( peak && peak->parentNuclide() )
      insertPeak( peak );
   }//for( int row = firstRow; row <= lastRow; ++row )
}//void peakModelRowsInsertedCallback(...)


void SourceFitModel::peakModelRowsRemovedCallback( Wt::WModelIndex /*index*/,
                                                    int firstRow, int lastRow )
{
  if( m_nuclides.empty() )
    return;

  vector<PeakShrdPtr> livingPeaks;
  const int nrow = m_peakModel->rowCount();
  
//  for( int row = firstRow; row <= lastRow; ++row )
  for( int row = 0; row < nrow; ++row )
  {
    if( row >= firstRow && row <= lastRow )
      continue;

    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );

    if( !peak )
    {  //probably wont ever happen
      cerr << "SourceFitModel::peakModelRowsRemovedCallback(...)\n\tNot a valid row being removed" << endl;
      continue;
    }//if( !peak )

    livingPeaks.push_back( peak );
  }//for( int row = firstRow; row <= lastRow; ++row )
  
  vector<int> removedRows;
  vector<const SandiaDecay::Nuclide *> removed;
  
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(row,0) );

    if( !peak )
      continue;

    const SandiaDecay::Nuclide *nuc = peak->parentNuclide();

    if( !nuc )
      continue;

    bool shouldRemove = true;
    for( PeakShrdPtr p : livingPeaks )
    {
      if( (p != peak) && (p->parentNuclide() == nuc) )
      {
        shouldRemove = false;
        break;
      }//if( (p != peak) && (p->parentNuclide() == nuc) )
    }//for( PeakShrdPtr p : livingPeaks )

    if( !shouldRemove )
      continue;

    if( std::find(removed.begin(), removed.end(), nuc) != removed.end() )
      continue;

    size_t index = m_nuclides.size();
    for( size_t i = 0; i < index; ++i )
      if( m_nuclides[i].nuclide == peak->parentNuclide() )
        index = i;
    
    if( index == m_nuclides.size() )
    {
      cerr << "SourceFitModel::peakModelRowsRemovedCallback(...)\n\tSerious logic error, fix this ish" << endl;
      cerr << "m_nuclides.size()=" << m_nuclides.size()
           << ", m_peakModel->npeaks()=" << m_peakModel->npeaks()
           << " nuc->symbol=" << nuc->symbol << endl;
      continue;
    }//if( index == m_nuclides.size() )

    removed.push_back( nuc );
    
    m_previousResults[nuc] = m_nuclides[index];
    
    const int iindex = static_cast<int>(index);
    removedRows.push_back( iindex );
    beginRemoveRows( WModelIndex(), iindex, iindex );
    m_nuclides.erase( m_nuclides.begin() + iindex );
    endRemoveRows();
  }//for( int row = firstRow; row <= lastRow; ++row )
  
  //Now we have to go back through and make sure the removed nuclides werent
  //  actually the age defining for an existing nuclide
  for( const SandiaDecay::Nuclide *nuc : removed )
  {
    const SandiaDecay::Nuclide *newDefining = nullptr;
    double minage = std::numeric_limits<double>::infinity();
    vector<const SandiaDecay::Nuclide *> nucstochange;
    
    for( IsoFitStruct &ifs : m_nuclides )
    {
      if( ifs.ageDefiningNuc == nuc )
      {
        nucstochange.push_back( ifs.nuclide );
        if( ifs.age < minage )
        {
          minage = ifs.age;
          newDefining = ifs.nuclide;
        }
      }//if( ifs.ageDefiningNuc == nuc )
    }//for( IsoFitStruct &ifs : m_nuclides )
    
    for( const SandiaDecay::Nuclide *nuc : nucstochange )
      setSharredAgeNuclide( nuc, newDefining );
  }//for( const SandiaDecay::Nuclide *nuc : removed )
  
//  for( const int iindex : removedRows )
//  {
//    beginRemoveRows( WModelIndex(), iindex, iindex );
//    m_nuclides.erase( m_nuclides.begin() + iindex );
//    endRemoveRows();
//  }//for( const int iindex : removeRows )
}//void peakModelRowsRemovedCallback( Wt::WModelIndex index, int firstRow, int lastRow )


bool is_within( int a, int start, int end )
{
  return (a>=start && a<=end);
}


void SourceFitModel::peakModelDataChangedCallback( Wt::WModelIndex topLeft,
                                                    Wt::WModelIndex bottomRight )
{
  const PeakModel::Columns colstart = PeakModel::Columns( topLeft.column() );
  const PeakModel::Columns colend = PeakModel::Columns( bottomRight.column() );
  const bool changedIso = is_within( PeakModel::kIsotope, colstart, colend);
  const bool changedFit = is_within( PeakModel::kUseForShieldingSourceFit, colstart, colend);
  
  const bool dirty = (changedIso || changedFit);
  
  if( !dirty )
  {
    const bool changedPhotoPeak
                    = is_within( PeakModel::kPhotoPeakEnergy, colstart, colend);
    const bool changedAmp = is_within( PeakModel::kAmplitude, colstart, colend);
    
    //Not dirty as far as this SourceFitModel is concerned, but as a convience
    //  to things hooked to this model, lets see if we chould redraw the chi2
    //  chart, and if do, emit a dataChanged() signal to indicate this.
    if( changedAmp || changedPhotoPeak )
    {
      try
      {
        cerr << "topLeft.isValid()=" << topLeft.isValid() << endl;
        const PeakModel::PeakShrdPtr peak = m_peakModel->peak(topLeft);
        cerr << "!!peak=" << !!peak << endl;
        if( !!peak )
          cerr << "peak->useForShieldingSourceFit()=" << peak->useForShieldingSourceFit() << endl;
        
        const bool beingUsed = (!!peak && peak->useForShieldingSourceFit());
        if( beingUsed )
          dataChanged().emit( WModelIndex(), WModelIndex() );
      }catch(...)
      {
        cerr << "Unexpected exception" << endl;
      }
    }//if( peak amplitude changed )
    
    return;
  }//if( !dirty )

  vector<const SandiaDecay::Nuclide *> preisotopes, postisotopes;
  for( const IsoFitStruct &ifs : m_nuclides )
    preisotopes.push_back( ifs.nuclide );
  
  const size_t npreisotopes = m_nuclides.size();
  
  repopulateIsotopes();
  
  for( const IsoFitStruct &ifs : m_nuclides )
    postisotopes.push_back( ifs.nuclide );
  
  if( m_sameAgeForIsotopes && (preisotopes.size()!=postisotopes.size()) )
  {
    vector<const SandiaDecay::Nuclide *> removednucs, addednucs;
    for( const SandiaDecay::Nuclide *nuc : preisotopes )
    {
      if( !std::count(postisotopes.begin(),postisotopes.end(),nuc) )
        removednucs.push_back( nuc );
    }
    
    for( const SandiaDecay::Nuclide *nuc : postisotopes )
    {
      if( !std::count(preisotopes.begin(),preisotopes.end(),nuc) )
        addednucs.push_back( nuc );
    }
    
    for( const SandiaDecay::Nuclide *nuc : removednucs )
    {
      bool removedADefining = false;
      const SandiaDecay::Nuclide *defining = NULL;
      double minage = std::numeric_limits<double>::infinity();
      
      for( IsoFitStruct &ifs : m_nuclides )
      {
        if( ifs.nuclide->atomicNumber == nuc->atomicNumber )
        {
          if( ifs.age < minage )
          {
            minage = ifs.age;
            defining = ifs.nuclide;
          }
          removedADefining = (removedADefining || (ifs.ageDefiningNuc==nuc));
        }//if( ifs.nuclide->atomicNumber == nuc->atomicNumber )
      }//for( IsoFitStruct &ifs : m_nuclides )
      
      if( removedADefining )
      {
        for( IsoFitStruct &ifs : m_nuclides )
        {
          if( ifs.nuclide->atomicNumber == nuc->atomicNumber )
            setSharredAgeNuclide( ifs.nuclide, defining );
        }//for( IsoFitStruct &ifs : m_nuclides )
      }//if( removedADefining )
    }//for( const SandiaDecay::Nuclide *nuc : removednucs )
    
    for( const SandiaDecay::Nuclide *nuc : addednucs )
    {
      const SandiaDecay::Nuclide *defining = NULL;
      double minage = std::numeric_limits<double>::infinity();
      
      for( IsoFitStruct &ifs : m_nuclides )
      {
        if( (ifs.age<minage) && (ifs.nuclide->atomicNumber==nuc->atomicNumber) )
        {
          minage = ifs.age;
          defining = ifs.nuclide;
        }
      }//for( IsoFitStruct &ifs : m_nuclides )
      
      if( !defining )
        defining = nuc;
      
      for( IsoFitStruct &ifs : m_nuclides )
        if( ifs.nuclide->atomicNumber == nuc->atomicNumber )
          setSharredAgeNuclide( ifs.nuclide, defining );
    }//for( const SandiaDecay::Nuclide *nuc : addednucs )
  }//if( m_sameAgeForIsotopes && (preisotopes!=postisotopes) )
  
  
  if( npreisotopes == m_nuclides.size() )
  { //No signals got emmited to update the chi2Chart, so lets check if we should
    try
    {
      const PeakModel::PeakShrdPtr &peak = m_peakModel->peak(topLeft);
      if( changedFit || (peak && peak->useForShieldingSourceFit()) )
        dataChanged().emit( index(0,0), index(rowCount()-1,columnCount()-1) );
      
    }catch(...){}
  }//if( npreisotopes == m_nuclides.size() )
  
}//void peakModelDataChangedCallback(...);


void SourceFitModel::peakModelResetCallback()
{
  if( m_nuclides.empty() )
    return;

  for( const IsoFitStruct &iss : m_nuclides )
    m_previousResults[iss.nuclide] = iss;
  
  beginRemoveRows( WModelIndex(), 0, static_cast<int>(m_nuclides.size())-1 );
  m_nuclides.clear();
  endRemoveRows();
}//void peakModelResetCallback();


void SourceFitModel::repopulateIsotopes()
{
  //go through and see if we have to delete any peaks
  vector<int> indexs_to_remove;
  const int norigoiso = static_cast<int>( m_nuclides.size() );

  for( int ison = 0; ison < norigoiso; ++ison )
  {
    IsoFitStruct &isof = m_nuclides[ison];

    size_t numSourcePeaks = 0; // I think we can remove this and instead use only progeny.size(), but havent tested
    set<const SandiaDecay::Nuclide *> progeny;
    const int npeaks = m_peakModel->rowCount();
    for( int peakn = 0; peakn < npeaks; ++peakn )
    {
      const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );

      if( (peak->parentNuclide() == isof.nuclide) && peak->useForShieldingSourceFit() )
      {
        const SandiaDecay::Transition * const trans = peak->nuclearTransition();
        if( trans && trans->parent )
          progeny.insert( trans->parent );
        numSourcePeaks += 1;
      }
    }//for( loop over peaks)

    isof.numProgenyPeaksSelected = progeny.size();
    
    if( numSourcePeaks == 0 )
    {
      m_previousResults[isof.nuclide] = isof;
      indexs_to_remove.push_back( ison );
    }
  }//for( const IsoFitStruct &isof : m_nuclides )

//  cerr << "indexs_to_remove.size()=" << indexs_to_remove.size() << endl;

  for( vector<int>::reverse_iterator iter = indexs_to_remove.rbegin();
       iter != indexs_to_remove.rend(); ++iter )
  {
    const int iindex = *iter;
    beginRemoveRows( WModelIndex(), iindex, iindex );
    m_nuclides.erase( m_nuclides.begin() + iindex );
    endRemoveRows();
  }//if( found != m_peakModel->rowCount() )


  //go through and see if we have to add any peaks
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    const PeakShrdPtr peak = m_peakModel->peak( m_peakModel->index(peakn,0) );
    if( !peak || !peak->useForShieldingSourceFit() )
       continue;

    bool found = false;
    for( const IsoFitStruct &isof : m_nuclides )
    {
      found = (peak->parentNuclide() == isof.nuclide);
      if( found )
        break;
    }//for( const IsoFitStruct &isof : m_nuclides )

    if( !found )
      insertPeak( peak );
  }//for( loop over peaks)
    
}//void SourceFitModel::repopulateIsotopes()



int SourceFitModel::columnCount( const Wt::WModelIndex & ) const
{
  return static_cast<int>( kNumColumns );
}//int columnCount( const Wt::WModelIndex & ) const


int SourceFitModel::rowCount( const Wt::WModelIndex & parent) const
{
  if (parent.isValid())
  {
    return 0;
  } //parent.isValid()
  return static_cast<int>( m_nuclides.size() );
}//int rowCount( const Wt::WModelIndex & ) const


Wt::WFlags<Wt::ItemFlag> SourceFitModel::flags( const Wt::WModelIndex &index ) const
{
  const int row = index.row();
  const int nrow = static_cast<int>( m_nuclides.size() );

  if( row>=nrow || row<0 )
    return WFlags<ItemFlag>();

  const IsoFitStruct &iso = m_nuclides[row];

  switch( index.column() )
  {
    case kIsotope:
      break;
    case kActivity:
      if( iso.shieldingIsSource )
        return WFlags<ItemFlag>();
      return WFlags<ItemFlag>(Wt::ItemIsEditable);
    case kFitActivity:
      if( iso.shieldingIsSource )
        return WFlags<ItemFlag>();
      return WFlags<ItemFlag>(ItemIsUserCheckable);
    case kAge:
//      if( iso.shieldingIsSource )
//        return WFlags<ItemFlag>();
      if( iso.ageDefiningNuc && iso.ageDefiningNuc!=iso.nuclide )
        return WFlags<ItemFlag>();
      return WFlags<ItemFlag>(Wt::ItemIsEditable);
    case kFitAge:
//      if( iso.shieldingIsSource )
//        return WFlags<ItemFlag>();
      return WFlags<ItemFlag>(ItemIsUserCheckable);
    case kIsotopeMass:
      return WFlags<ItemFlag>();
    case kActivityUncertainty:
    case kAgeUncertainty:
      return WFlags<ItemFlag>();
      
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:
    case kTruthActivityTolerance:
    case kTruthAge:
    case kTruthAgeTolerance:
      return WFlags<ItemFlag>(Wt::ItemIsEditable);
#endif
      
    case kNumColumns:
      break;
  }//switch( section )

  return WFlags<ItemFlag>();
}//Wt::WFlags<Wt::ItemFlag> flags( const Wt::WModelIndex &index ) const




boost::any SourceFitModel::headerData( int section, Orientation orientation, int role ) const
{
  //When orientation is Horizontal, section is a column number,
  //  when orientation is Vertical, section is a row (peak) number.

  if( role == LevelRole )
    return 0;

  if( (orientation != Horizontal)
      || ((role != DisplayRole)
           &&  (role != Wt::ToolTipRole))  )
    return WAbstractItemModel::headerData( section, orientation, role );

  if( role == Wt::ToolTipRole )
  {
    const char *tooltip = nullptr;
    switch( section )
    {
      case kIsotope:
        tooltip = "Nuclides are added or removed according to the photopeaks "
                  "used for this fit, thus may not be edited in this table";
                  //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kActivity:
        tooltip = "Activity of the nuclide.  Entering a reasonable starting"
            " value may help the fit, or the activity does not have to be"
            " fit for."
            "  Values may be entered in formats similar to '12.4 kBq',"
            " '1.0mCi', '55.3uCi', '15.82Mbq'";
            //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kAge:
        tooltip = "Age of the nuclide.  Values may be entered in format"
            " similar to: '23.3s' '2 hl' - (hf stands for half lives), '5.23y',"
            " '23:14:21.343', '19min', etc.  Note that nuclides that decay to"
            " stable children can not be aged.";
            //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kFitActivity:
        tooltip = "Should the fit try to find the nuclide activity as well, or"
                  " assume the entered value is correct?";
                  //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kFitAge:
        tooltip = "Should the fit try to find the nuclide age as well, or just"
                  " assume the entered value is correct?  This option is not "
                  "available for nuclides which decay to only stable chlidren.";
                  //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kIsotopeMass:
        tooltip = "This is the mass of the bare nuclide assuming the activity"
                  " listed.";
                  //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kActivityUncertainty:
        tooltip = "This is the uncetainty in activity from the fit.";
        //NAZ add to end: "Select arrow to sort table by parameter." 
        break;
      case kAgeUncertainty:
        tooltip = "This is the uncetainty in age from the fit.";
        //NAZ add to end: "Select arrow to sort table by parameter." 
        break;

#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case kTruthActivity:
      case kTruthActivityTolerance:
      case kTruthAge:
      case kTruthAgeTolerance:
        break;
#endif
        
      case kNumColumns:
        break;
    }//switch( section )
    if( !tooltip )
      return boost::any();
    
    return boost::any( WString::fromUTF8(tooltip) );
  }//if( role == Wt::ToolTipRole )

  //If we're here, role==DisplayRole
  switch( section )
  {
    case kIsotope:     return boost::any( WString("Nuclide") );
    case kActivity:    return boost::any( WString("Activity") );
    case kAge:         return boost::any( WString("Age") );
    case kFitActivity: return boost::any( WString("Fit Act.") );
    case kFitAge:      return boost::any( WString("Fit Age") );
    case kIsotopeMass: return boost::any( WString("Mass") );
    case kActivityUncertainty: return boost::any( WString("Act. Uncert.") );
    case kAgeUncertainty:      return boost::any( WString("Age Uncert.") );
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:          return boost::any( WString("Truth Act.") );
    case kTruthActivityTolerance: return boost::any( WString("Truth Act. Tol.") );
    case kTruthAge:               return boost::any( WString("Truth Age") );
    case kTruthAgeTolerance:      return boost::any( WString("Truth Age Tol.") );
#endif
    case kNumColumns:  return boost::any();
  }//switch( section )

  return boost::any();
}//headerData(...)


Wt::WModelIndex SourceFitModel::parent( const Wt::WModelIndex &index ) const
{
  return WModelIndex();
}//Wt::WModelIndex parent(...) const


boost::any SourceFitModel::data( const Wt::WModelIndex &index, int role ) const
{
  //should consider impementing ToolTipRole
  if( role != Wt::DisplayRole && role != Wt::EditRole && role != Wt::ToolTipRole
      && !((role==Wt::CheckStateRole) && ((index.column()==kFitActivity) || (index.column()==kFitAge)))
      )
    return boost::any();


  const int row = index.row();
  const int column = index.column();
  const int nrows = static_cast<int>( m_nuclides.size() );
  
  if( row<0 || column<0 || column>=kNumColumns || row>=nrows )
    return boost::any();

  const IsoFitStruct &isof = m_nuclides[row];

  if( role == Wt::ToolTipRole )
  {
    if( column == kIsotope )
    {
      string msg = "HalfLife="
                 + PhysicalUnits::printToBestTimeUnits( isof.nuclide->halfLife, 2, PhysicalUnits::second );
      //XXX - should put in dominant photopeaks or other information here
      return boost::any( WString(msg) );
    }else if( column == kAge )
    {
      if( !isof.ageIsFittable )
        return boost::any( WString("Nuclides children are stable so aging is not allowed because it doesnt effect gamma spectrum.") );
    }//if / else to detrmine column

    return boost::any();
  }//if( role == Wt::ToolTipRole )

  switch( column )
  {
    case kIsotope:
      return boost::any( WString(isof.nuclide->symbol) );
    case kActivity:
    {
      double act = isof.activity * GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
      const string ans = PhysicalUnits::printToBestActivityUnits( act, 2, m_displayCurries );
      return boost::any( WString(ans) );
    }//case kActivity:

    case kFitActivity:
    {
      if( isof.shieldingIsSource )
        return boost::any();
      return boost::any( isof.fitActivity );
    }//case kFitActivity:

    case kAge:
    {
//      if( isof.shieldingIsSource )
//        return boost::any();
      double age;
      const SandiaDecay::Nuclide *nuc = nullptr;
      if( isof.ageDefiningNuc && (isof.ageDefiningNuc != isof.nuclide) )
      {
        const int ind = nuclideIndex( isof.ageDefiningNuc );
        if( ind >= 0 )
        {
          nuc = isof.ageDefiningNuc;
          age = m_nuclides[ind].age;
        }else
        {
          nuc = isof.nuclide;
          age = isof.age;
          cerr << "SourceFitModel::data: ran into error when retriving"
               << " age for a nuclide with a defining age isotope that isnt in"
               << " the model; charging on!" << endl;
        }//if( ind >= 0 )
      }else
      {
        nuc = isof.nuclide;
        age = isof.age;
      }//if( isof.ageDefiningNuc && (isof.ageDefiningNuc!=isof.nuclide) )
      
      if( !isof.ageIsFittable )
        return boost::any( WString("NA") );
      const string ans = PhysicalUnits::printToBestTimeUnits( age, 2 );
      return boost::any( WString(ans) );
    }//case kAge:

    case kFitAge:
    {
//      if( isof.shieldingIsSource )
//        return boost::any();
      if( isof.ageDefiningNuc && (isof.ageDefiningNuc != isof.nuclide) )
        return boost::any();
      
      // Make sure there is more than two peaks being fitted for this nuclide to enable fitting age
      //  (I guess you could fix activity, and shielding, and select a progeny peak, and fit for
      //   age based on that peak growing in, but this probably isnt realistically ever done, but if
      //   you did want to do it, you could round-about calculate it)
      if( (isof.numProgenyPeaksSelected <= 1) && !isof.ageDefiningNuc )
        return boost::any();
      
      if( !isof.ageIsFittable )
        return boost::any();
      
      return boost::any( isof.fitAge );
    }//case kFitAge:

    case kIsotopeMass:
    {
      const double act = isof.activity * GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
      const double mass_grams = act / isof.nuclide->activityPerGram();

      if( IsInf(mass_grams) || IsNan(mass_grams) )
        return boost::any();

      return boost::any( WString(PhysicalUnits::printToBestMassUnits(mass_grams,4,1.0)) );
    }//case kIsotopeMass:

    case kActivityUncertainty:
    {
      if( isof.activityUncertainty < 0.0 )
        return boost::any();
      
      double act = isof.activityUncertainty * GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
      const string ans = PhysicalUnits::printToBestActivityUnits( act, 2, m_displayCurries );
      return boost::any( WString(ans) );
    }//case kActivityUncertainty:

    case kAgeUncertainty:
    {
      if( (!isof.ageIsFittable) || isof.ageUncertainty < 0.0 )
        return boost::any();
      const string ans = PhysicalUnits::printToBestTimeUnits( isof.ageUncertainty, 2 );
      return boost::any( WString(ans) );
    }//case kAgeUncertainty:

#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:
    {
      if( !isof.truthActivity )
        return boost::any();
      const string ans = PhysicalUnits::printToBestActivityUnits( *isof.truthActivity, 4, m_displayCurries );
      return boost::any( WString(ans) );
    }
      
    case kTruthActivityTolerance:
    {
      if( !isof.truthActivityTolerance )
        return boost::any();
      const string ans = PhysicalUnits::printToBestActivityUnits( *isof.truthActivityTolerance, 4, m_displayCurries );
      return boost::any( WString(ans) );
    }
      
    case kTruthAge:
    {
      if( !isof.truthAge )
        return boost::any();
      const string ans = PhysicalUnits::printToBestTimeUnits( *isof.truthAge, 4 );
      return boost::any( WString(ans) );
    }
      
    case kTruthAgeTolerance:
    {
      if( !isof.truthAgeTolerance )
        return boost::any();
      const string ans = PhysicalUnits::printToBestTimeUnits( *isof.truthAgeTolerance, 4 );
      return boost::any( WString(ans) );
    }
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )
      
    case kNumColumns:
      return boost::any();
  }//switch( column )

  return boost::any();
}//boost::any data(...) const


Wt::WModelIndex SourceFitModel::index( int row, int column,
                                        const Wt::WModelIndex &/*parent*/ ) const
{
  return WAbstractItemModel::createIndex( row, column, (void *)0 );
}//Wt::WModelIndex index(...) const



WModelIndex SourceFitModel::index( const SandiaDecay::Nuclide *nuc,
                                       SourceFitModel::Columns column ) const
{
  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int row = 0; row < nrow; ++row )
    if( m_nuclides[row].nuclide == nuc )
      return index( row, column );
  return WModelIndex();
}//WModelIndex index(...) const


WModelIndex SourceFitModel::index( const string &symbol,
                                   SourceFitModel::Columns column ) const
{
  const int nrow = static_cast<int>( m_nuclides.size() );
  for( int row = 0; row < nrow; ++row )
    if( m_nuclides[row].nuclide && (m_nuclides[row].nuclide->symbol == symbol) )
      return index( row, column );
  return WModelIndex();
}//WModelIndex index(...) const


bool SourceFitModel::setData( const Wt::WModelIndex &index, const boost::any &value, int role )
{
  try
  {
    if( !index.isValid() )
      return false;

    if( role!=Wt::EditRole && role!=Wt::CheckStateRole )
      return false;

    const int row = index.row();
    const int column = index.column();
    const int nrows = static_cast<int>( m_nuclides.size() );

#if( INCLUDE_ANALYSIS_TEST_SUITE )
    if( value.empty() )
    {
      switch( SourceFitModel::Columns(column) )
      {
        case kTruthActivity:
        case kTruthActivityTolerance:
        case kTruthAge:
        case kTruthAgeTolerance:
        case kActivityUncertainty:
        case kAgeUncertainty:
          break;
        default:
          return false;
      }//switch( column )
    }//if( value.empty() )
#else
    if( value.empty() && column!=kAgeUncertainty && column!=kActivityUncertainty )
      return false;
#endif  //#if( INCLUDE_ANALYSIS_TEST_SUITE )
    
    
    if( row<0 || column<0 || column>=kNumColumns || row>=nrows )
      return false;
    
    if( role==Wt::CheckStateRole && column!=kFitActivity && column!=kFitAge )
      return false;
    
    if( role==Wt::EditRole && (column==kFitActivity || column==kFitAge) )
      return false;

    const WString txt_val = asString( value );

    if( column!=kFitActivity && column!=kFitAge && txt_val.empty() )
      return false;

    IsoFitStruct &iso = m_nuclides[row];

    //If were here, all is fine
    switch( column )
    {
      case kIsotope:
        return false;
      case kActivity:
        iso.activity = PhysicalUnits::stringToActivity( txt_val.toUTF8() );
        iso.activity /= GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
        
        if( iso.activityUncertainty >= 0.0 )
          iso.activityUncertainty = -1.0;
      break;

      case kFitActivity:
        iso.fitActivity = boost::any_cast<bool>( value );
      break;

      case kAge:
      {
        if( !iso.ageIsFittable )
        {
          break;
        }else
        {
          const string str = txt_val.toUTF8();
          const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
          iso.age = PhysicalUnits::stringToTimeDurationPossibleHalfLife( str, hl );

          if( iso.ageUncertainty >= 0.0 )
          {
            iso.ageUncertainty = -1.0;
            WModelIndex uncertIndex = createIndex(row, kAgeUncertainty,(void *)0);
            dataChanged().emit( uncertIndex, uncertIndex );
          }//if( iso.ageUncertainty >= 0.0 )
          
          if( m_sameAgeForIsotopes )
          {
            for( size_t i = 0; i < m_nuclides.size(); ++i )
            {
              const int thisrow = static_cast<int>(i);
              if( thisrow == row )
                continue;
              
              if( m_nuclides[i].ageDefiningNuc == iso.nuclide )
              {
                WModelIndex ind = createIndex(thisrow, kAge,(void *)0);
                dataChanged().emit( ind, ind );
                if( iso.ageUncertainty >= 0.0 )
                {
                  ind = createIndex(thisrow, kAgeUncertainty,(void *)0);
                  dataChanged().emit( ind, ind );
                }//if( iso.ageUncertainty >= 0.0 )
              }//if( nuc.ageDefiningNuc == iso.nuclide )
            }//for( IsoFitStruct *nuc : m_nuclides )
          }//if( m_sameAgeForIsotopes )
        }//if( decays to stable children / else )
        break;
      }//case kAge:

      case kFitAge:
        if( !iso.ageIsFittable )
        {
          if( iso.fitAge )
            iso.fitAge = false;
          else
            return false;
        }else
        {
          const bool bvalue = boost::any_cast<bool>( value );
          if( bvalue == iso.fitAge )  //dont change this if we dont have to
            return false;
          iso.fitAge = bvalue;
        }//if( !ageIsFittable ) / else
      break;

      case kIsotopeMass:
        cerr << "SourceFitModel::setData(...)\n\tYou shouldnt be trying to set kIsotopeMass"
             << endl;
      break;

      case kActivityUncertainty:
        iso.activityUncertainty = -1.0;
        if( !value.empty() )
        {
          iso.activityUncertainty = PhysicalUnits::stringToActivity( txt_val.toUTF8() );
          iso.activityUncertainty /= GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
        }//if( !value.empty() )
      break;

      case kAgeUncertainty:
      {
        iso.ageUncertainty = -1.0;
        if( iso.ageIsFittable && !value.empty() )
        {
          const string str = txt_val.toUTF8();
          const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
          iso.ageUncertainty = PhysicalUnits::stringToTimeDurationPossibleHalfLife( str, hl );
        }//if( decays to stable children / else )
        
        if( m_sameAgeForIsotopes )
        {
          for( size_t i = 0; i < m_nuclides.size(); ++i )
          {
            if( m_nuclides[i].ageDefiningNuc == iso.nuclide )
            {
              WModelIndex ind = createIndex( static_cast<int>(i), kAgeUncertainty, (void *)0);
              dataChanged().emit( ind, ind );
            }//if( nuc.ageDefiningNuc == iso.nuclide )
          }//for( IsoFitStruct *nuc : m_nuclides )
        }//if( m_sameAgeForIsotopes )
        
        break;
      }//case kAgeUncertainty:

#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case kTruthActivity:
        if( value.empty() )
          iso.truthActivity.reset();
        else
          iso.truthActivity = PhysicalUnits::stringToActivity( txt_val.toUTF8() );
        break;
        
      case kTruthActivityTolerance:
        if( value.empty() )
          iso.truthActivityTolerance.reset();
        else
          iso.truthActivityTolerance = PhysicalUnits::stringToActivity( txt_val.toUTF8() );
        break;
        
      case kTruthAge:
      {
        const string str = txt_val.toUTF8();
        const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
        if( value.empty() )
          iso.truthAge.reset();
        else
          iso.truthAge = PhysicalUnits::stringToTimeDurationPossibleHalfLife( str, hl );
        break;
      }
         
      case kTruthAgeTolerance:
      {
        const string str = txt_val.toUTF8();
        const double hl = (iso.nuclide ? iso.nuclide->halfLife : -1.0);
        if( value.empty() )
          iso.truthAgeTolerance.reset();
        else
          iso.truthAgeTolerance = PhysicalUnits::stringToTimeDurationPossibleHalfLife( str, hl );
        break;
      }
#endif
        
      case kNumColumns:
        return false;
    }//switch( column )

    if( column==kActivity )
      dataChanged().emit( createIndex(row,0,(void *)0),
                          createIndex(row,kNumColumns-1,(void *)0) );
    else
      dataChanged().emit( index, index );
  }catch( exception &e )
  {
    cerr << "SourceFitModel::setData(...)\n\tWarning: exception caught; what="
         << e.what() << endl;
    return false;
  }//try/catch

  return true;
}//bool setData(...)


bool SourceFitModel::compare( const IsoFitStruct &lhs,
                                const IsoFitStruct &rhs,
                                Columns sortColumn, Wt::SortOrder order )
{
  bool isLess = false;
  
#if( INCLUDE_ANALYSIS_TEST_SUITE )
  auto optionalLess = []( const boost::optional<double> &olhs, const boost::optional<double> &orhs) -> bool{
    if( (!olhs) != (!orhs) )
      return !olhs;
    
    if( !olhs )
      return false;
      
    return ((*olhs) < (*orhs));
  };
#endif
  
  
  switch( sortColumn )
  {
    case kIsotope:    isLess = (lhs.nuclide->symbol<rhs.nuclide->symbol); break;
    case kActivity:   isLess = (lhs.activity < rhs.activity);             break;
    case kFitActivity:isLess = (lhs.fitActivity < rhs.fitActivity);       break;
    case kAge:        isLess = (lhs.age < rhs.age);                       break;
    case kFitAge:     isLess = (lhs.fitAge < rhs.fitAge);                 break;
    case kIsotopeMass:
      isLess = ((lhs.activity/lhs.nuclide->activityPerGram())
              < (rhs.activity/rhs.nuclide->activityPerGram()) );
    break;
    case kActivityUncertainty:
      isLess = (lhs.activityUncertainty<rhs.activityUncertainty);
    break;
    case kAgeUncertainty:
      isLess = (lhs.ageUncertainty<rhs.ageUncertainty);
    break;
      
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    case kTruthActivity:
      isLess = optionalLess( lhs.truthActivity, rhs.truthActivity );
      break;
      
    case kTruthActivityTolerance:
      isLess = optionalLess( lhs.truthActivityTolerance, rhs.truthActivityTolerance );
      break;
      
    case kTruthAge:
      isLess = optionalLess( lhs.truthAge, rhs.truthAge );
      break;
      
    case kTruthAgeTolerance:
      isLess = optionalLess( lhs.truthAgeTolerance, rhs.truthAgeTolerance );
      break;
#endif
      
    case kNumColumns:
      isLess = false;
      break;
  }//switch( sortColumn )

  if(order == Wt::AscendingOrder)
      return !isLess;
  return isLess;
}//bool compare(...);

void SourceFitModel::sort( int column, Wt::SortOrder order )
{
  layoutAboutToBeChanged().emit();
  m_sortOrder = order;
  m_sortColumn = Columns( column );
  std::sort( m_nuclides.begin(), m_nuclides.end(),
             boost::bind( &SourceFitModel::compare, _1, _2, m_sortColumn, m_sortOrder ) );
  layoutChanged().emit();
}//void sort(...)


ShieldingSourceDisplay::Chi2Graphic::Chi2Graphic( Wt::WContainerWidget *parent )
  : Wt::Chart::WCartesianChart( parent ),
    m_nFitForPar( 0 ),
    m_showChi( true ),
    m_textPenColor( Wt::black )
{
  setPreferredMethod( WPaintedWidget::HtmlCanvas );
  LOAD_JAVASCRIPT( wApp, "shieldingSourceDisplay.cpp", "Chi2Graphic", wtjsShowChi2Info);
}

ShieldingSourceDisplay::Chi2Graphic::~Chi2Graphic()
{
}


void ShieldingSourceDisplay::Chi2Graphic::setColorsFromTheme( std::shared_ptr<const ColorTheme> theme )
{
  if( !theme )
    return;
  
  //if( !theme->foregroundLine.isDefault() )
  //  m_chartEnergyLineColor = theme->foregroundLine;
  
  if( !theme->spectrumChartText.isDefault() )
  {
    WPen txtpen(theme->spectrumChartText);
    this->setTextPen( txtpen );
    this->axis(Chart::XAxis).setTextPen( txtpen );
    this->axis(Chart::YAxis).setTextPen( txtpen );
    this->setTextPenColor( theme->spectrumChartText );
  }
  
  if( theme->spectrumChartBackground.isDefault() )
    this->setBackground( Wt::NoBrush );
  else
    this->setBackground( WBrush(theme->spectrumChartBackground) );
  
  if( (theme->spectrumChartMargins.isDefault() && !theme->spectrumChartBackground.isDefault()) )
  {
    //theme->spectrumChartBackground
  }else if( !theme->spectrumChartMargins.isDefault() )
  {
    //theme->spectrumChartMargins
  }
  
  if( !theme->spectrumAxisLines.isDefault() )
  {
    WPen defpen = this->axis(Chart::XAxis).pen();
    defpen.setColor( theme->spectrumAxisLines );
    this->axis(Chart::XAxis).setPen( defpen );
    this->axis(Chart::YAxis).setPen( defpen );
  }
}//ShieldingSourceDisplay::Chi2Graphic::setColorsFromTheme( theme )


void ShieldingSourceDisplay::Chi2Graphic::setNumFitForParams( unsigned int npar )
{
  m_nFitForPar = static_cast<int>( npar );
}

void ShieldingSourceDisplay::Chi2Graphic::setTextPenColor( const Wt::WColor &color )
{
  m_textPenColor = color;
}

void ShieldingSourceDisplay::Chi2Graphic::setShowChiOnChart( const bool show_chi )
{
  m_showChi = show_chi;
  removeSeries(1);
  removeSeries(2);
  setXSeriesColumn( 0 );
  
  Chart::WDataSeries series( (show_chi ? 1 : 2), Chart::PointSeries );
  series.setMarkerSize( 10.0 );
  
  addSeries( series );
}

void ShieldingSourceDisplay::Chi2Graphic::calcAndSetAxisPadding( double yHeightPx )
{
  //Wt newer than ~3.3.2 has WAbstractChart::setAutoLayoutEnabled(true) that
  //  may make this function no longer necassarry.
  
  WAbstractItemModel *theModel = model();

  if( !theModel )
    return;

  double ymin = DBL_MAX, ymax = -DBL_MAX;

//The below doesnt actually return the minimum/maximum of the axis, presumambly
//  since we are using auto range.
//  initLayout();
//  double ymin = axis(Chart::OrdinateAxis).minimum();
//  double ymax = axis(Chart::OrdinateAxis).maximum();
  
  const int nrow = theModel->rowCount();
  const int ycol = m_showChi ? 1 : 2;
  for( int row = 0; row < nrow; ++row )
  {
    try
    {
#if( WT_VERSION >= 0x3030800 )
      const double thischi = theModel->data(row,ycol);
#else
      WModelIndex index = theModel->index(row,ycol);
      const double thischi = boost::any_cast<double>( theModel->data(index) );
#endif
      ymin = std::min( ymin, thischi );
      ymax = std::max( ymax, thischi );
    }catch(...)
    {
    }
  }//for( int row = 0; row < nrow; ++row )
  
  if( !nrow || ymin==DBL_MAX || ymax==-DBL_MAX )
  {
    ymin = 0.0;
    ymax = 100.0;
  }//
  

//Calculate number of pixels we need to pad, for x axis to be at 40 pixels from
//  the bottom of the chart; if less than 10px, pad at least 10px, or at most
//  40px.
  const int topPadding = plotAreaPadding(Top);
  const double fracY = -ymin / (ymax - ymin);
  double pxToXAxis = (fracY >= 1.0) ? 0.0 : (40.0 - fracY*(yHeightPx - topPadding) ) / (1.0-fracY);
  pxToXAxis = std::floor( pxToXAxis + 0.5 );
  pxToXAxis = std::max( pxToXAxis, 10.0 );
  pxToXAxis = std::min( pxToXAxis, 40.0 );
  const int bottomPadding = static_cast<int>(pxToXAxis);
  
  if( bottomPadding != plotAreaPadding(Bottom) )
    setPlotAreaPadding( bottomPadding, Bottom );
  
  yHeightPx -= (topPadding + bottomPadding);
  
  {//Begin codeblock to determin min/max label values, following logic
   //  determined from round*125() and other places in WAxis.C
    const int numLabels = yHeightPx / 25.0;
    const double range = ymax - ymin;
    const double rangePerLabel = range / numLabels;

    double renderinterval = 1.0;
    double n = std::pow(10, std::floor(std::log10(rangePerLabel)));
    double msd = rangePerLabel / n;
    
    if (msd < 1.5)
      renderinterval = n;
    else if (msd < 3.3)
      renderinterval = 2*n;
    else if (msd < 7)
      renderinterval = 5*n;
    else
      renderinterval = 10*n;
    
    ymin += std::numeric_limits<double>::epsilon();
    ymax -=  std::numeric_limits<double>::epsilon();
    ymin = renderinterval * std::floor( ymin / renderinterval);
    ymax = renderinterval * std::ceil( ymax / renderinterval);
  }//End codeblock to determin min/max label values

  
  const int oldleft = plotAreaPadding(Wt::Left);
  const WString minlabel = axis(Chart::Y1Axis).label(ymin);
  const WString maxlabel = axis(Chart::Y1Axis).label(ymax);
  const size_t maxnchars = minlabel.narrow().length();
  const size_t minnchars = maxlabel.narrow().length();
  const size_t nchars = std::max( minnchars, maxnchars );
  const int left = 32 + 6*std::max(static_cast<int>(nchars),2);
  
  if( left != oldleft )
    setPlotAreaPadding( left, Wt::Left );
}//void calcAndSetAxisPadding()


void ShieldingSourceDisplay::Chi2Graphic::paintEvent( WPaintDevice *device )
{
  calcAndSetAxisPadding( device->height().toPixels() );
  
  Wt::Chart::WCartesianChart::paintEvent( device );
}//void paintEvent( Wt::WPaintDevice *paintDevice )


void ShieldingSourceDisplay::Chi2Graphic::paint( Wt::WPainter &painter,
                      const Wt::WRectF &rectangle ) const
{
  WCartesianChart::paint( painter, rectangle );

  //I think removing of the areas() is already done by
  //  WCartesuanChart::paintEvent(...), but jic
  while( !areas().empty() )
    delete areas().front();
  
  WStandardItemModel *chi2Model = dynamic_cast<WStandardItemModel *>( model() );
  if( !chi2Model )  //prob never happen, bust JIC
    return;
  
  const WPointF br = painter.window().bottomRight();
  const int width = br.x();

  const int nrow = chi2Model->rowCount();
  double chi2 = 0.0;
  
  for( int row = 0; row < nrow; ++row )
  {
    try
    {
      WModelIndex index = chi2Model->index(row,1);
      const double thischi = boost::any_cast<double>( chi2Model->data(index) );
      chi2 += thischi*thischi;
      

      WColor color;
      try
      {
        boost::any color_any = chi2Model->data(index, Wt::MarkerPenColorRole );
        color = boost::any_cast<WColor>(color_any);
        if( color.isDefault() )
          throw runtime_error("");
      }catch(...)
      {
        //I dont think we will ever get here, but JIC I guess.
        Wt::Chart::WChartPalette *pal = palette();
        if( pal )
          color = pal->brush(0).color();
        else
          color = WColor( Wt::darkRed );
      }//try / catch, get the color
      
      index = chi2Model->index(row,2);
      const double thisscale = boost::any_cast<double>( chi2Model->data(index) );
      
      index = chi2Model->index(row,0);
      double energy = boost::any_cast<double>( chi2Model->data(index) );
      
      
      const double yval = m_showChi ? thischi : thisscale;
      const WPointF pos = mapToDevice( energy, yval );
      
      index = chi2Model->index(row,3);
      const string nucname = Wt::asString( chi2Model->data(index) ).toUTF8();
      energy = ((100.0*energy+0.5)/100.0);
      
      char mouseoverjs[256];
      if( nucname.empty() || nucname.size() > 8 )
      {
        snprintf( mouseoverjs, sizeof(mouseoverjs),
                  "function(){Wt.WT.ShowChi2Info('%s','%.2f keV, &sigma;<sub>%i</sub>=%.2f, %.2fx model');}",
                  id().c_str(), energy, row, thischi, thisscale );
      }else
      {
        snprintf( mouseoverjs, sizeof(mouseoverjs),
                 "function(){Wt.WT.ShowChi2Info('%s','%s, %.2f keV, &sigma;<sub>%i</sub>=%.2f, %.2fx model');}",
                 id().c_str(), nucname.c_str(), energy, row, thischi, thisscale );
      }
      
      //We could potentially include more info here, but it will be kinda hard,
      //  so lets not worry about it now.
      WRectArea *area = new WRectArea( pos.x()-10.0, pos.y()-10.0, 20.0, 20.0 );
      const string mouseoutjs = "function(){$('#" + id() + "inf').remove();}";
      area->mouseWentOver().connect( string(mouseoverjs) );
      area->mouseWentOut().connect( mouseoutjs );
      //area->clicked()
      
      //I hate doing const_casts, but the Wt source code itself does this
      //  for adding image areas in Wt 3.3.2
      const_cast<Chi2Graphic*>(this)->addArea( area );
      
      if( !m_showChi )
      {
        index = chi2Model->index(row,4);
        const double scale_uncert = boost::any_cast<double>( chi2Model->data(index) );
        
        const WPointF upper_uncert = mapToDevice( energy, yval + scale_uncert );
        const WPointF lower_uncert = mapToDevice( energy, yval - scale_uncert );
        const WPen oldPen = painter.pen();
        painter.setPen( WPen(color) );
        painter.drawLine( upper_uncert, lower_uncert );
        painter.setPen( oldPen );
      }//if( !m_showChi )
    }catch( std::exception &e )
    {
      cerr << "Caught exception drawing Chi2 chart: " << e.what() << endl;
    }
  }//for( int row = 0; row < nrow; ++row )
  
  if( nrow > 0 && !IsNan(sqrt(chi2)) )
  {
    char buffer[64];
    //snprintf( buffer, sizeof(buffer), "&lt;dev&gt;=%.2g&sigma;", (sqrt(chi2)/nrow) );  //displays as literal tesxt, and not the symbols on android
    //snprintf( buffer, sizeof(buffer), "\x3c\xCF\x87\x3E=%.2g\xCF\x83", (sqrt(chi2)/nrow) ); //<χ>=13.2σ
    snprintf( buffer, sizeof(buffer), "\x3c\x64\x65\x76\x3E=%.2g\xCF\x83", (sqrt(chi2)/nrow) ); //<dev>=13.2σ
  
    WString text = WString::fromUTF8(buffer);
    const size_t msglen = SpecUtils::utf8_str_len( buffer, strlen(buffer) );
  
    const double rightPadding = static_cast<double>( plotAreaPadding(Right) );
    const double charwidth = 16;
    const double charheight = 15;
    const double x = width - charwidth*msglen - rightPadding - 5;
    const double y = plotAreaPadding(Top) + 5;
    const double twidth = charwidth*msglen;
    const double theight = charheight + 2;
    
    WPen oldPen = painter.pen();
    painter.setPen( WPen(m_textPenColor) );
    painter.drawText( x, y, twidth, theight, AlignRight, TextSingleLine, text );
    painter.setPen( oldPen );
  }//if( nrow > 0 && !IsNan(sqrt(chi2)) )
  
  if( !m_showChi && axis(Chart::YAxis).minimum() < 1.0 && axis(Chart::YAxis).maximum() > 1.0 )
  {
    WPointF left = mapToDevice( axis(Chart::XAxis).minimum(), 1.0 );
    WPointF right = mapToDevice( axis(Chart::XAxis).maximum(), 1.0 );
    left.setY( left.y() + 0.5 );
    right.setY( right.y() + 0.5 );
    
    WPen pen( GlobalColor::lightGray );
    pen.setWidth( 1 );
    pen.setStyle( Wt::PenStyle::DashLine );
    const WPen oldPen = painter.pen();
    painter.setPen( pen );
    painter.drawLine( left, right );
    painter.setPen( oldPen );
  }//if( m_showChi )
  
  //Draw the y-axis label
  painter.rotate( -90 );
  painter.setPen( WPen() );
  painter.setBrush( WBrush() );
  
  const char *yaxistitle = m_showChi ? "(observed-model)/uncert" : "Mult. of Model";
  
  const char *chi2tooltip  = "The Y-axis is the peak area counts, minus the expected number of counts from the currently displayed"
                             " source and shielding values, all divided by the statistical uncertainty.  The square of this"
                             " is what is used to optimize the shielding and sources.";
  const char *scaletooltip = "The Y-axis is the multiple of observed peak counts, relative to what would be expected from the"
                             " current source and shielding values.";
  
  const char * const yaxistooltip = m_showChi ? chi2tooltip : scaletooltip;
  
  WPen oldPen = painter.pen();
  painter.setPen( WPen(m_textPenColor) );
  painter.drawText( -0.45*painter.viewPort().height(), 0.0, 0.2, 0.1,
                     AlignCenter, yaxistitle );
  painter.setPen( oldPen );
  
  const double yAxisWiddth = plotAreaPadding(Wt::Left);
  const double chartHeight = painter.viewPort().height();
  
  WRectArea *yAxisArea = new WRectArea( 0.0, 0.0, yAxisWiddth, chartHeight );
  yAxisArea->setToolTip( yaxistooltip );
  const_cast<Chi2Graphic*>(this)->addArea( yAxisArea );
  
  painter.restore();
}//Chi2Graphic::paint(


pair<ShieldingSourceDisplay *,AuxWindow *> ShieldingSourceDisplay::createWindow( InterSpec *viewer )
{
  assert( viewer );
  
  AuxWindow *window = nullptr;
  ShieldingSourceDisplay *disp = nullptr;
  
  try
  {
    MaterialDB *matdb = viewer->materialDataBase();
    PeakModel *peakModel = viewer->peakModel();
    WSuggestionPopup *shieldSuggest = viewer->shieldingSuggester();
    
    disp = new ShieldingSourceDisplay( peakModel, viewer, shieldSuggest, matdb );
    window = new AuxWindow( "Activity/Shielding Fit" );

    window->setResizable( true );
    window->contents()->setOffsets(WLength(0,WLength::Pixel));
    window->stretcher()->addWidget( disp, 0, 0 );
    window->stretcher()->setContentsMargins(0,0,0,0);
  //    window->footer()->resize(WLength::Auto, WLength(50.0));
      
    WPushButton *closeButton = window->addCloseButtonToFooter();
    closeButton->clicked().connect(window, &AuxWindow::hide);
      
    AuxWindow::addHelpInFooter( window->footer(), "activity-shielding-dialog" );
    
    window->rejectWhenEscapePressed();
      
    //Should take lock on m_dataMeasurement->mutex_
    
    shared_ptr<SpecMeas> meas = viewer->measurment(SpecUtils::SpectrumType::Foreground);
    rapidxml::xml_document<char> *shield_source = nullptr;
    if( meas )
      shield_source = meas->shieldingSourceModel();
      
    if( shield_source && shield_source->first_node() )
    {
      //string msg = "Will try to deserailize: \n";
      //rapidxml::print( std::back_inserter(msg), *shield_source->first_node(), 0 );
      //cout << msg << endl << endl;
      try
      {
        disp->deSerialize( shield_source->first_node() );
      }catch( std::exception &e )
      {
        string xmlstring;
        rapidxml::print(std::back_inserter(xmlstring), *shield_source, 0);
        stringstream debugmsg;
        debugmsg << "Error loading Shielding/Source model: "
                    "\n\tError Message: " << e.what()
                 << "\n\tModel XML: " << xmlstring;
  #if( PERFORM_DEVELOPER_CHECKS )
        log_developer_error( __func__, debugmsg.str().c_str() );
  #else
        cerr << debugmsg.str() << endl;
  #endif
        passMessage( "There was an error loading the previous shielding/source model for this file - model state is suspect!",
                      "", WarningWidget::WarningMsgHigh );
      }
    }//if( shield_source )
    
  //    m_shieldingSourceFitWindow->resizeScaledWindow( 0.75, 0.75 );
      
    const double windowWidth = 0.95 * viewer->renderedWidth();
    const double windowHeight = 0.95 * viewer->renderedHeight();
      
  //    double footerheight = m_shieldingSourceFitWindow->footer()->height().value();
  //    m_shieldingSourceFitWindow->setMinimumSize( WLength(200), WLength(windowHeight) );
      
    if( (windowHeight > 100) && (windowWidth > 100) )
    {
      if( !viewer->isPhone() )
        window->resizeWindow( windowWidth, windowHeight );

      //Give the m_shieldingSourceFitWindow a hint about what size it will be
      //  rendered at so it can decide what widgets should be rendered - acounting
      //  for borders and stuff (roughly)
      disp->initialSizeHint( windowWidth - 12, windowHeight - 28 - 50 );
    }else if( !viewer->isPhone() )
    {
      //When loading an application state that is showing this window, we may
      //  not know the window size (e.g., windowWidth==windowHeight==0), so
      //  instead skip giving the initial size hint, and instead size things
      //  client side (maybe we should just do this always?)
      window->resizeScaledWindow( 0.95, 0.95 );
    }
        
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    WPushButton *setTruth = new WPushButton( "Set Truth Values" );
    window->footer()->insertWidget( 0, setTruth );
    setTruth->clicked().connect( disp, &ShieldingSourceDisplay::showInputTruthValuesWindow );
#endif
    
    
  //   m_shieldingSourceFitWindow->contents()->  setHeight(WLength(windowHeight));

    window->centerWindow();

    window->finished().connect( viewer, &InterSpec::closeShieldingSourceFitWindow );
      
    window->WDialog::setHidden(false);
    window->show();
    window->centerWindow();
  }catch( std::exception &e )
  {
    passMessage( "Error creating Activity/Shielding fit display: " + string(e.what()),
                 "", WarningWidget::WarningMsgHigh );
    
    if( disp )
      delete disp;
    disp = nullptr;
    
    if( window )
      AuxWindow::deleteAuxWindow( window );
    window = nullptr;
  }//try / catch
  
  return make_pair( disp, window );
}//createWindow( InterSpec *viewer  )




ShieldingSourceDisplay::ShieldingSourceDisplay( PeakModel *peakModel,
                                                InterSpec *specViewer,
                                                WSuggestionPopup *matSuggest,
                                                MaterialDB *materialDB,
                                                WContainerWidget *parent )
  : WContainerWidget( parent ),
    m_chi2ChartNeedsUpdating( true ),
    m_width( 0 ),
    m_height( 0 ),
    m_nResizeSinceHint( 0 ),
    m_modifiedThisForeground( false ),
    m_peakModel( peakModel ),
    m_specViewer( specViewer ),
    m_sourceModel( nullptr ),
    m_peakView( nullptr ),
    m_sourceView( nullptr ),
    m_detectorDisplay( nullptr ),
    m_distanceEdit( nullptr ),
    m_addMaterialShielding( nullptr ),
    m_addGenericShielding( nullptr ),
    m_layout( nullptr ),
    m_addItemMenu( nullptr ),
#if( USE_DB_TO_STORE_SPECTRA )
    m_saveAsNewModelInDb( nullptr ),
#endif
    m_materialSuggest( matSuggest ),
    m_shieldingSelects( nullptr ),
    m_showChi2Text( nullptr ),
    m_chi2Model( nullptr ),
    m_chi2Graphic( nullptr ),
    m_multiIsoPerPeak( nullptr ),
    m_backgroundPeakSub( nullptr ),
    m_sameIsotopesAge( nullptr ),
    m_showChiOnChart( nullptr ),
    m_optionsDiv( nullptr ),
    m_showLog( nullptr ),
    m_logDiv( nullptr ),
    m_materialDB( materialDB ),
    m_fitModelButton( nullptr ),
    m_fitProgressTxt( nullptr ),
    m_cancelfitModelButton( nullptr )
{
  wApp->useStyleSheet( "InterSpec_resources/ShieldingSourceDisplay.css" );
  
  const bool showToolTips = InterSpecUser::preferenceValue<bool>( "ShowTooltips", m_specViewer );
  
  setLayoutSizeAware( true );
  const bool isotopesHaveSameAge = true;
  m_sourceModel = new SourceFitModel( m_peakModel, isotopesHaveSameAge, this );
  m_peakView = new RowStretchTreeView();
  m_peakView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
  
  m_peakView->setAlternatingRowColors( true );
  m_peakView->setEditTriggers( WAbstractItemView::SingleClicked | WAbstractItemView::DoubleClicked );
  m_peakView->setModel( m_peakModel );
//  m_peakView->addStyleClass( "ShieldingSourceDisplay_m_peakView" );

  for( PeakModel::Columns col = PeakModel::Columns(0);
       col < PeakModel::kNumColumns;
       col = PeakModel::Columns(col+1) )
  {
    switch( col )
    {
      case PeakModel::kMean: case PeakModel::kUseForShieldingSourceFit:
      case PeakModel::kPhotoPeakEnergy: case PeakModel::kIsotope:
        m_peakView->setColumnHidden( col, false );
      break;
      default:
        m_peakView->setColumnHidden( col, true );
      break;
    }//switch( col )
  }//for( loop over peak columns )

  WItemDelegate *dblDelagate = new WItemDelegate( this );
  dblDelagate->setTextFormat( "%.2f" );
  m_peakView->setItemDelegateForColumn( PeakModel::kMean, dblDelagate );

  m_peakView->setColumnWidth( PeakModel::kMean, WLength(7.5,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kPhotoPeakEnergy, WLength(10,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kIsotope, WLength(9,WLength::FontEx) );
  m_peakView->setColumnWidth( PeakModel::kUseForShieldingSourceFit, WLength(7,WLength::FontEx) );


  PhotopeakDelegate *nuclideDelegate = new PhotopeakDelegate( PhotopeakDelegate::NuclideDelegate, true, m_peakView );
  m_peakView->setItemDelegateForColumn( PeakModel::kIsotope, nuclideDelegate );


  m_sourceView = new RowStretchTreeView();
  m_sourceView->setRootIsDecorated	(	false); //makes the tree look like a table! :)
  

  m_sourceView->setModel( m_sourceModel );
  m_sourceView->setSortingEnabled( true );
  m_sourceView->setAlternatingRowColors( true );
  m_sourceView->addStyleClass( "ShieldingSourceDisplay_m_sourceView" );

  for( SourceFitModel::Columns col = SourceFitModel::Columns(0);
       col < SourceFitModel::kNumColumns;
       col = SourceFitModel::Columns(col+1) )
  {
    switch( col )
    {
      case SourceFitModel::kActivity:
        //need make custom delegate
      break;
      case SourceFitModel::kAge:
        //need to make custom delegate
      break;

      case SourceFitModel::kFitActivity: case SourceFitModel::kFitAge:
      case SourceFitModel::kIsotope: case SourceFitModel::kNumColumns:
      case SourceFitModel::kIsotopeMass:
      case SourceFitModel::kActivityUncertainty:
      case SourceFitModel::kAgeUncertainty:
#if( INCLUDE_ANALYSIS_TEST_SUITE )
      case SourceFitModel::kTruthActivity: case SourceFitModel::kTruthActivityTolerance:
      case SourceFitModel::kTruthAge: case SourceFitModel::kTruthAgeTolerance:
#endif
      break;
    }//case( col )
  }//for( loop over SourceFitModel columns )


  m_sourceView->setColumnWidth( SourceFitModel::kActivity, WLength(9,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kAge, WLength(9,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kFitAge, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kFitActivity, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kIsotope, WLength(9,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kIsotopeMass, WLength(9,WLength::FontEx) );

  m_sourceView->setColumnWidth( SourceFitModel::kActivityUncertainty, WLength(10,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kAgeUncertainty, WLength(10,WLength::FontEx) );

#if( INCLUDE_ANALYSIS_TEST_SUITE )
  m_sourceView->setColumnWidth( SourceFitModel::kTruthActivity, WLength(13,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthActivityTolerance, WLength(14,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthAge, WLength(13,WLength::FontEx) );
  m_sourceView->setColumnWidth( SourceFitModel::kTruthAgeTolerance, WLength(14,WLength::FontEx) );
#endif
  
  
  m_detectorDisplay = new DetectorDisplay( m_specViewer, m_specViewer->fileManager()->model() );
  m_detectorDisplay->setInline( true );

  Wt::WPushButton *addItemMenubutton = new WPushButton( " " );
  m_addItemMenu = new PopupDivMenu( addItemMenubutton, PopupDivMenu::TransientMenu );
  if( m_specViewer->isMobile() )
  {
     m_addItemMenu->addPhoneBackItem( NULL );
     addItemMenubutton->clicked().connect( m_addItemMenu,
                                            &PopupDivMenu::showFromClick );
  } //mobile
  else
  {
      m_addItemMenu->setButton( addItemMenubutton );
  } //not mobile
  addItemMenubutton->setStyleClass( "ShiledingSourceOptionsBtn" );
  addItemMenubutton->addStyleClass( "GearIcon" );

  //this validates floating point numbers followed by a distance unit
  WRegExpValidator *distValidator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
  distValidator->setFlags( Wt::MatchCaseInsensitive );
  WLabel *distanceLabel = new WLabel( "Distance:" );
  
  m_distanceEdit = new WLineEdit( "100 cm" );
  m_distanceEdit->setTextSize( 5 );

  distanceLabel->setBuddy( m_distanceEdit );
  m_distanceEdit->setValidator( distValidator );
  string tooltip = "Distance from center of source to face of detector. Number must be"
            " followed by units; valid units are: meters, m, cm, mm, km, feet,"
            " ft, ', in, inches, or \".  You may also add multiple distances,"
            " such as '3ft 4in', or '3.6E-2 m 12 cm' which are equivalent to "
            " 40inches and 15.6cm respectively.";

  HelpSystem::attachToolTipOn( m_distanceEdit,tooltip, showToolTips );


//  WContainerWidget *shieldingDiv = new WContainerWidget();
  m_shieldingSelects = new WContainerWidget();
  m_shieldingSelects->setStyleClass( "ShieldingSelectContainer" );

  tooltip = "Shieldings are treated as spherical shells with the innermost"
            " shell being the first one you add (the top one in this box), and"
            " the subsequent shieldings layered around that one.  Any sources"
            " not attributed as components of shieldings are treated as point"
            " sources in the center.";
  HelpSystem::attachToolTipOn( m_shieldingSelects,tooltip, showToolTips );

#define  SPLIT_BUTTON_SHIELDING_ADD 0

#if( SPLIT_BUTTON_SHIELDING_ADD )
  WSplitButton *addShielding = new WSplitButton();
  tooltip = "Adds a shielding that allows you to select from a database, or"
            " enter a chemical formula.  See drop down menu to add a generic"
            " shielding defined by atomic number and areal density.";
  HelpSystem::attachToolTipOn( addShielding->actionButton(), tooltip, showToolTips );
  addShielding->actionButton()->setText( "Add Shielding" );
  addShielding->actionButton()->clicked()
                .connect( this, &ShieldingSourceDisplay::doAddShielding );
  PopupDivMenu *shieldMenu = new PopupDivMenu( addShielding->dropDownButton(), PopupDivMenu::TransientMenu );
  WMenuItem *materialitem = shieldMenu->addItem( "Material" );
  materialitem->triggered().connect( this, &ShieldingSourceDisplay::doAddShielding );
  materialitem->setIcon( "InterSpec_resources/images/shield_add.png" );
  static_assert( 0, "Disabling SPLIT_BUTTON_SHIELDING_ADD during fits not done yet" );
  WMenuItem *genericitem = shieldMenu->addItem( "Generic AN, AD" );
  genericitem->triggered().connect( this, &ShieldingSourceDisplay::addGenericShielding );
  genericitem->setIcon( "InterSpec_resources/images/shape_square_add.png" );
  
  //has style class "dropdown-toggle"
#else
  
  WLabel *addShieldingLabel = new WLabel( "Add Shielding:" );
  m_addMaterialShielding = new WPushButton( "Material" );
  HelpSystem::attachToolTipOn( m_addMaterialShielding,
              "Choose from a library of predefined common shielding materials.",
                              showToolTips, HelpSystem::ToolTipPosition::Top  );
  //m_addMaterialShielding->setStyleClass("ShieldAddIcon");
  m_addMaterialShielding->setIcon( "InterSpec_resources/images/shield_white.png" );
  m_addMaterialShielding->clicked().connect( this,
                                      &ShieldingSourceDisplay::doAddShielding );
  
  m_addGenericShielding = new WPushButton( "Generic" );
  HelpSystem::attachToolTipOn( m_addGenericShielding,
              "Allows you to define and fit for atomic number and areal density.",
                              showToolTips , HelpSystem::ToolTipPosition::Top );
  //m_addGenericShielding->setStyleClass("ShapeSquareAddIcon");
  m_addGenericShielding->setIcon( "InterSpec_resources/images/atom_white.png" );
  m_addGenericShielding->clicked().connect( this,
                                     &ShieldingSourceDisplay::addGenericShielding );
#endif
  
  
  m_fitModelButton = new WPushButton( "Perform Model Fit" );
  m_fitModelButton->clicked().connect( boost::bind(&ShieldingSourceDisplay::doModelFit, this, true) );

  m_fitProgressTxt = new WText();
  m_fitProgressTxt->hide();
  
  m_cancelfitModelButton = new WPushButton( "Cancel Model Fit" );
  m_cancelfitModelButton->clicked().connect( boost::bind( &ShieldingSourceDisplay::cancelModelFit, this ) );
  m_cancelfitModelButton->hide();
  
  m_showLog = m_addItemMenu->addMenuItem( "Calculation Log" );
  m_showLog->triggered().connect( this, &ShieldingSourceDisplay::showCalcLog );
  m_showLog->disable();

  PopupDivMenuItem *item = NULL;
//  PopupDivMenuItem *item = m_addItemMenu->addMenuItem( "Test Serialization" );
//  item->triggered().connect( this, &ShieldingSourceDisplay::testSerialization );

if (m_specViewer->isSupportFile())
{
  item = m_addItemMenu->addMenuItem( "Import Model..." );
  item->triggered().connect( this, &ShieldingSourceDisplay::startModelUpload );
  
  StringDownloadResource *xmlResource = new StringDownloadResource( this );
  item = m_addItemMenu->addMenuItem( "Export Model" );
  item->setLink( WLink( xmlResource ) );
  item->setLinkTarget(Wt::TargetNewWindow);
} //isSupportFile()
  
#if( USE_DB_TO_STORE_SPECTRA )
  item = m_addItemMenu->addMenuItem( "Open From Database..." );
  item->triggered().connect( this,
                          &ShieldingSourceDisplay::startBrowseDatabaseModels );
  
  item = m_addItemMenu->addMenuItem( "Save To Database..." );
  item->triggered().connect(
                boost::bind( &ShieldingSourceDisplay::startSaveModelToDatabase,
                             this, false) );
  
  m_saveAsNewModelInDb = m_addItemMenu->addMenuItem( "Clone In Database..." );
  m_saveAsNewModelInDb->triggered().connect( this,
                            &ShieldingSourceDisplay::saveCloneModelToDatabase );
  m_saveAsNewModelInDb->disable();
#endif //#if( USE_DB_TO_STORE_SPECTRA )
  
  m_showChi2Text = new WText( "(make window taller to show &chi;<sup>2</sup>)",
                              XHTMLUnsafeText );
  m_showChi2Text->setInline( false );
  m_showChi2Text->hide();

  m_chi2Model = new WStandardItemModel( 0, 6, parent );
  m_chi2Graphic = new Chi2Graphic();
  m_chi2Graphic->setModel( m_chi2Model );

  m_chi2Graphic->setPlotAreaPadding( 40, Bottom );
  m_chi2Graphic->setPlotAreaPadding( 12, Right );
  m_chi2Graphic->setPlotAreaPadding(  2, Top );
  m_chi2Graphic->setPlotAreaPadding( 50, Left );
  //m_chi2Graphic->setAutoLayoutEnabled();
  
  m_chi2Graphic->addSeries( Chart::WDataSeries(1, Chart::PointSeries) );
  m_chi2Graphic->setXSeriesColumn( 0 );
  m_chi2Graphic->axis(Chart::XAxis).setTitle( "Energy (keV)" );
  
  WFont font( WFont::Default );
  font.setSize( WFont::Small );
  m_chi2Graphic->axis(Chart::YAxis).setTitleFont(font);
  m_chi2Graphic->axis(Chart::XAxis).setTitleFont(font);
  
  WFont labelFont( WFont::Default );
  labelFont.setSize( WFont::Medium );
  m_chi2Graphic->axis(Chart::YAxis).setLabelFont(font);
  m_chi2Graphic->axis(Chart::XAxis).setLabelFont(font);
  
  
  m_chi2Graphic->axis(Chart::XAxis).setLocation( Chart::ZeroValue );

  m_chi2Graphic->axis(Chart::XAxis).setScale( Chart::LinearScale );
  m_chi2Graphic->axis(Chart::YAxis).setScale( Chart::LinearScale );

  m_chi2Graphic->setType( Chart::ScatterPlot );
//  m_chi2Graphic->setMinimumSize( WLength(200), WLength(175) );
  
  
  //We should check the color theme for colors
  m_specViewer->colorThemeChanged().connect( m_chi2Graphic, &Chi2Graphic::setColorsFromTheme );
  m_chi2Graphic->setColorsFromTheme( m_specViewer->getColorTheme() );
  
  
  //The next line is kinda inefficient because if all that changed was
  //  fit activity or fit age, then we dont really need to update the chi2 chart
  m_sourceModel->dataChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  //XXX -
  //  When the below is connected, the order the isotpes are added to
  //  m_sourceModel is somewhat off, so there is an issue updating the
  //  m_chi2Graphic because of a mismatch, in what isotopes are actually in
  //  m_sourceModel
  m_sourceModel->rowsInserted().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  
  //XXX - this next line causes a runtime_exception to be raised under some
  //      circumstances, causing the Chi2 chart to not be updated - leaving for
  //      now, but I'm not really sure as to the sequence of events that causes
  //      this to happen (model gets out of sync with the chi2 functions)
  m_sourceModel->rowsRemoved().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  m_sourceModel->rowsInserted().connect( this, &ShieldingSourceDisplay::addSourceIsotopesToShieldings );
  m_sourceModel->rowsAboutToBeRemoved().connect( this, &ShieldingSourceDisplay::removeSourceIsotopesFromShieldings );

  m_distanceEdit->changed().connect( this, &ShieldingSourceDisplay::handleUserDistanceChange );
  m_distanceEdit->enterPressed().connect( this, &ShieldingSourceDisplay::handleUserDistanceChange );
  
  m_specViewer->detectorChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  m_specViewer->detectorModified().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  
  m_optionsDiv = new WContainerWidget();
  WGridLayout* optionsLayout = new WGridLayout();
  m_optionsDiv->setLayout(optionsLayout);
  optionsLayout->setContentsMargins(0, 0, 0, 0);
  
  WContainerWidget *allpeaksDiv = new WContainerWidget();
  WCheckBox *allpeaks = new WCheckBox( "All Peaks", allpeaksDiv );
  allpeaks->setAttributeValue( "style", "white-space:nowrap;margin-right:5px;float:right;" + allpeaks->attributeValue("style") );
  optionsLayout->addWidget( allpeaksDiv, 0, 0 );
  allpeaks->setTristate( true );
  allpeaks->changed().connect( boost::bind( &ShieldingSourceDisplay::toggleUseAll, this, allpeaks ) );
  m_peakModel->dataChanged().connect(
                                     boost::bind( &ShieldingSourceDisplay::updateAllPeaksCheckBox, this, allpeaks ) );
  updateAllPeaksCheckBox( allpeaks ); //initialize
  
  
  m_optionsDiv->setOverflow( WContainerWidget::OverflowHidden );
  //The ToolTip of WCheckBoxes is a bit finicky, and only works over the
  //  checkbox itself, so lets make it work over the label to, via lineDiv
  WContainerWidget *lineDiv = new WContainerWidget();
  optionsLayout->addWidget( lineDiv, 1, 0 );
  m_multiIsoPerPeak = new WCheckBox( "Multiple nuclides contribute to peak", lineDiv );
  tooltip = "Checking this option will allow other nuclides being fit for,"
  " other that the one assigned to a peak, to contribute to the"
  " expected counts in a observed peak.  This could for instance be"
  " useful when fitting for both Ho166m and U235 that both have a"
  " peak near 185 keV." ;
  lineDiv->setToolTip( tooltip );
  m_multiIsoPerPeak->setChecked();
  m_multiIsoPerPeak->checked().connect( this, &ShieldingSourceDisplay::multiNucsPerPeakChanged );
  
  lineDiv = new WContainerWidget(  );
  optionsLayout->addWidget( lineDiv, 2, 0 );
  m_backgroundPeakSub = new WCheckBox( "Subtract Background Peaks", lineDiv );
  tooltip = "This forces the isotopes of a element to all be the same age in"
            " the fit.";
  lineDiv->setToolTip( tooltip );
  m_backgroundPeakSub->checked().connect( this, &ShieldingSourceDisplay::backgroundPeakSubChanged );
  
  
  lineDiv = new WContainerWidget(  );
  optionsLayout->addWidget( lineDiv, 3, 0 );
  m_sameIsotopesAge = new WCheckBox( "Isotopes of same element same age", lineDiv );
  tooltip = "Enforce isotopes for the same element should all have the same "
            "age.";
  lineDiv->setToolTip( tooltip );
  m_sameIsotopesAge->checked().connect( this, &ShieldingSourceDisplay::sameIsotopesAgeChanged );
  m_sameIsotopesAge->unChecked().connect( this, &ShieldingSourceDisplay::sameIsotopesAgeChanged );
  m_sameIsotopesAge->setChecked( isotopesHaveSameAge );

  
  lineDiv = new WContainerWidget(  );
  optionsLayout->addWidget( lineDiv, 4, 0 );
  m_showChiOnChart = new WCheckBox( "Show Chi Graphic", lineDiv );
  m_showChiOnChart->setChecked();
  tooltip = "Show the Chi of each peak for the current model on the chart, or"
            " the relative peak area multiple between current model and observed peak.";
  lineDiv->setToolTip( tooltip );
  m_showChiOnChart->checked().connect( this, &ShieldingSourceDisplay::showGraphicTypeChanged );
  m_showChiOnChart->unChecked().connect( this, &ShieldingSourceDisplay::showGraphicTypeChanged );
  
  
  WContainerWidget *detectorDiv = new WContainerWidget();
  detectorDiv->setOverflow(WContainerWidget::OverflowHidden);
  WGridLayout *detectorLayout = new WGridLayout();
  detectorDiv->setLayout( detectorLayout );
  
  
  WContainerWidget* smallerContainer = new WContainerWidget();
  WGridLayout *smallLayout = new WGridLayout();
  smallerContainer->setLayout(smallLayout);
  
  
#if( SPLIT_BUTTON_SHIELDING_ADD )
  distanceLabel->setMargin( 40, Wt::Left );
  smallLayout->addWidget( distanceLabel,           0, 0, AlignRight | AlignMiddle);
  smallLayout->addWidget( m_distanceEdit,          0, 1, 1, 1 );
  smallLayout->addWidget( addShielding,            1, 1, AlignCenter );
  smallLayout->setColumnStretch( 1, 10 );
#else
  smallLayout->addWidget( distanceLabel,           0, 0, AlignRight | AlignMiddle);
  smallLayout->addWidget( m_distanceEdit,          0, 1, 1, 2);
  smallLayout->addWidget( addShieldingLabel,       1, 0, AlignRight | AlignMiddle);
  smallLayout->addWidget( m_addMaterialShielding,  1, 1);
  smallLayout->addWidget( m_addGenericShielding,   1, 2);
  smallLayout->setColumnStretch( 0, 0 );
#endif

  smallLayout->setContentsMargins(0,5,0,5);
  smallerContainer->setPadding(0);
  
  
  
  //---------------
  
  WContainerWidget *peakDiv = new WContainerWidget();
  
  WGridLayout *peakGrid = new Wt::WGridLayout();
  peakGrid->setRowStretch(0, 1);
  peakGrid->setColumnStretch(0, 1);
  peakDiv->setLayout(peakGrid);
  peakGrid->addWidget(m_peakView,0,0);
  
  m_layout = new WGridLayout();
  
  WContainerWidget *bottomLeftDiv = new WContainerWidget();
  WGridLayout *bottomLeftLayout = new WGridLayout();
  bottomLeftDiv->setLayout(bottomLeftLayout);
  bottomLeftLayout->addWidget( peakDiv, 0, 0 );
  bottomLeftLayout->addWidget( m_optionsDiv, 1, 0 );
  bottomLeftLayout->setRowStretch( 0, 1 );
  bottomLeftLayout->setVerticalSpacing( 0 );
  bottomLeftLayout->setHorizontalSpacing( 0 );
  bottomLeftLayout->setContentsMargins( 0, 0, 0, 0 );
  
  WContainerWidget *sourceDiv = new WContainerWidget();
  WGridLayout *sourceGrid = new Wt::WGridLayout();
  sourceGrid->setRowStretch(0, 1);
  sourceGrid->setColumnStretch(0, 1);
  sourceDiv->setLayout(sourceGrid);
  sourceGrid->addWidget(m_sourceView,0,0);
  
  if( m_specViewer->isPhone() )
  {
    //phone layout
    detectorLayout->addWidget( smallerContainer,           0, 0);
    detectorLayout->addWidget( m_shieldingSelects,         1, 0);
    
    detectorLayout->setRowStretch( 1, 1 );
    detectorLayout->setColumnStretch( 0, 1 );
    detectorLayout->setHorizontalSpacing( 0 );
    detectorLayout->setVerticalSpacing( 0 );
    detectorLayout->setContentsMargins( 1, 1, 1, 1 );
    
    WTabWidget* tab = new WTabWidget();
    tab->setOffsets(0);
    tab->setMargin(0);
    m_layout->addWidget(tab, 0,0);
    
    tab->addTab(bottomLeftDiv,"Source Peaks", Wt::WTabWidget::PreLoading);
    tab->addTab(sourceDiv,"Source Isotopes", Wt::WTabWidget::PreLoading);
    tab->addTab(detectorDiv,"Shielding", Wt::WTabWidget::PreLoading);
    
    WContainerWidget * chartDiv = new WContainerWidget();
    chartDiv->setOffsets(0);
    chartDiv->setMargin(0);
    chartDiv->setPadding(5);
    WGridLayout* chartLayout= new WGridLayout();
    chartDiv->setLayout(chartLayout);
    chartLayout->setContentsMargins(0, 0, 0, 0);
    
    chartLayout->addWidget( m_detectorDisplay, 0, 0);
    chartLayout->addWidget( addItemMenubutton, 0, 1);
    chartLayout->addWidget( m_chi2Graphic, 1, 0,1,2);
    chartLayout->addWidget( m_fitModelButton, 2, 0, 1 , 2 );
    chartLayout->addWidget( m_fitProgressTxt, 3, 0, 1, 2 );
    chartLayout->addWidget( m_cancelfitModelButton, 4, 0, 1, 2 );
    
    chartLayout->setRowStretch(1, 1);
    chartLayout->setColumnStretch(0, 1);
    tab->addTab(chartDiv,"Fit", Wt::WTabWidget::PreLoading);
  }else
  {
    //regular layout
    detectorLayout->addWidget( m_detectorDisplay,          0, 0 );
    detectorLayout->addWidget( addItemMenubutton,          0, 1 );
    detectorLayout->addWidget( m_fitModelButton,           3, 0, 1, 2, AlignCenter );
    detectorLayout->addWidget( m_fitProgressTxt,           4, 0, 1, 2 );
    detectorLayout->addWidget( m_cancelfitModelButton,     5, 0, 1, 2, AlignCenter );
    
    detectorLayout->addWidget( smallerContainer,           1, 0, 1, 2);
    detectorLayout->addWidget( m_shieldingSelects,         2, 0, 1, 2 );
    detectorLayout->addWidget( m_showChi2Text,             6, 0, 1, 2 );
    
    detectorLayout->setRowStretch( 2, 1 );
    detectorLayout->setColumnStretch( 0, 1 );
    detectorLayout->setHorizontalSpacing( 0 );
    detectorLayout->setVerticalSpacing( 0 );
    detectorLayout->setContentsMargins( 1, 1, 1, 1 );
    
    WGridLayout *bottomMiddleLayout = new WGridLayout();
    
    bottomMiddleLayout->addWidget( bottomLeftDiv,   0, 0 );
    bottomMiddleLayout->addWidget( sourceDiv, 0, 1 );
    bottomMiddleLayout->setColumnResizable( 0, true, WLength(38.25,WLength::FontEx) );
    bottomMiddleLayout->setHorizontalSpacing( 5 );
    bottomMiddleLayout->setVerticalSpacing( 5 );
    bottomMiddleLayout->setContentsMargins( 0, 0, 0, 0 );
    
    WGridLayout *leftLayout = new WGridLayout();
    leftLayout->addWidget( m_chi2Graphic,      0, 0 );
    leftLayout->addLayout( bottomMiddleLayout,   1, 0 );
    leftLayout->setRowResizable( 0, true, WLength(40.0,WLength::Percentage) );
    leftLayout->setHorizontalSpacing( 5 );
    leftLayout->setVerticalSpacing( 5 );
    leftLayout->setContentsMargins( 0, 0, 0, 0 );
    
    WContainerWidget *leftDiv = new WContainerWidget();
    leftDiv->setLayout(leftLayout);
    leftDiv->setOverflow(WContainerWidget::OverflowHidden);
    
    detectorDiv->setWidth( 290 );
    
    m_layout->addWidget( leftDiv, 0, 0);
    m_layout->addWidget( detectorDiv, 0, 1 );
    m_layout->setColumnStretch( 0, 1 );
    m_layout->setHorizontalSpacing( 0 );
    m_layout->setVerticalSpacing( 0 );
    m_layout->setContentsMargins( 0, 0, 0, 0 );
  } //regular layout
  
  m_layout->setContentsMargins(0, 0, 0, 0);
  setLayout( m_layout );
  setOverflow( WContainerWidget::OverflowVisible);
  setOffsets(WLength(0,WLength::Pixel));
  updateChi2Chart();
}//ShieldingSourceDisplay constructor


//When the button is triggered, update model
void ShieldingSourceDisplay::toggleUseAll( Wt::WCheckBox *button )
{
  const bool useForFit = button->isChecked();
  const size_t npeaks = m_peakModel->npeaks();
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    try
    {
      const PeakModel::PeakShrdPtr peak = m_peakModel->peakPtr( i );
      WModelIndex index = m_peakModel->indexOfPeak( peak );
      
      if( !!peak && index.isValid() && peak->parentNuclide() )
      {
        index = m_peakModel->index( index.row(), PeakModel::kUseForShieldingSourceFit );
        m_peakModel->setData( index, useForFit );
      }
    }catch( std::exception & )
    {
      //shouldnt ever happen, wont worry about
    }
  }//for( size_t i = 0; i < npeaks; ++i )
}//void ShieldingSourceDisplay::toggleUseAll(Wt::WCheckBox* button)


//When model is changed, update checkbox tristate
void ShieldingSourceDisplay::updateAllPeaksCheckBox( WCheckBox *but)
{
    bool allon=true;
    bool alloff=true;
    
    for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    {
        WModelIndex index = m_peakModel->index( peakn,
                                               PeakModel::kUseForShieldingSourceFit );
        const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index);
        if( peak->useForShieldingSourceFit() && alloff)
        {
            alloff=false;
        } //if( peak->useForShieldingSourceFit() && alloff)
        else if( !peak->useForShieldingSourceFit() && allon)
        {
            allon=false;
        } //if( !peak->useForShieldingSourceFit() && allon)
    }//for
    
    if (alloff && !allon)
        but->setCheckState(Wt::Unchecked);
    else if (allon && !alloff)
        but->setCheckState(Wt::Checked);
    else
        but->setCheckState(Wt::PartiallyChecked);
} //updateAllPeaksCheckBox( WCheckBox *but)


void ShieldingSourceDisplay::render( Wt::WFlags<Wt::RenderFlag> flags )
{
  const bool renderFull = (flags & Wt::RenderFlag::RenderFull);
  
  WContainerWidget::render( flags );
  
  if( m_chi2ChartNeedsUpdating )
  {
    updateChi2ChartActual();
    m_chi2ChartNeedsUpdating = false;
  }
}//void render( Wt::WFlags<Wt::RenderFlag> flags )


ShieldingSourceDisplay::~ShieldingSourceDisplay() noexcept(true)
{
  {//begin make sure calculation is cancelled
    std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
    if( m_currentFitFcn )
    {
      try
      {
        m_currentFitFcn->cancelFit();
      }catch( ... )
      {
        cerr << "Caught exception call m_currentFitFcn->cancelFit(), which probably shouldnt happen" << endl;
      }
    }
  }//end make sure calculation is cancelled
  
  if( m_addItemMenu )
  {
    delete m_addItemMenu;
    m_addItemMenu = NULL;
  }//if( m_addItemMenu )
}//ShieldingSourceDisplay destructor constructor


#if( INCLUDE_ANALYSIS_TEST_SUITE )
void ShieldingSourceDisplay::showInputTruthValuesWindow()
{
  //The error handling and display quality is minimal since this is a developer-only tool, and will
  //  not be used by end-users.
  //Also, if you change the model any while this window is open - bad things will happen.
  
  AuxWindow *window = new AuxWindow( "Input Truth Values",
                                     (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                                     | AuxWindowProperties::TabletModal) );
  
  WContainerWidget *contents = window->contents();
  
  try
  {
    WTable *table = new WTable( contents );
    table->setHeaderCount( 1 );
    new WLabel( "Quantitiy", table->elementAt(0, 0) );
    new WLabel( "Value", table->elementAt(0, 1) );
    new WLabel( "Tolerance", table->elementAt(0, 2) );
    
    const int nnuc = m_sourceModel->numNuclides();
    for( int i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
      assert( nuc );
      
      const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
      const bool selfAttNuc = m_sourceModel->shieldingDeterminedActivity( i );
      
//      if( selfAttNuc )
//        throw runtime_error( "Model is not candidate for truth-level info<br />"
//                             "Self-attuating sources not implemented yet" );
//      if( ageNuc && (ageNuc != nuc) )
//        throw runtime_error( "Model is not candidate for truth-level info<br />"
//                             "Shared-age nuclides not allowed" );
      
      // For self-attenuating shieldings, we'll just test the shielding thickness
      // For nuclides whose age is controlled by another nuclide, we dont need to test age.
      if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
        continue;
      
      const bool fitAct = m_sourceModel->fitActivity(i);
      const bool fitAge = m_sourceModel->fitAge(i);
      
      auto setFieldValue = [nuc,this]( WLineEdit *valuefld, const SourceFitModel::Columns type ){
        const int srcrow = m_sourceModel->nuclideIndex(nuc);
        auto index = m_sourceModel->index( srcrow, type );
        valuefld->setText( asString( m_sourceModel->data(index) ) );
      };//setFieldValue(...)
      
      auto fieldUpdate = [this,nuc]( WLineEdit *valuefld, const SourceFitModel::Columns type ){
        const string valtxt = valuefld->text().toUTF8();
        
        const int srcrow = m_sourceModel->nuclideIndex(nuc);
        WModelIndex index = m_sourceModel->index( srcrow, type );
        
        try
        {
          m_sourceModel->setData(index, valtxt);
        }catch(...)
        {
          passMessage( "'" + valtxt + "' is not a valid entry", "", WarningWidget::WarningMsgHigh );
          valuefld->setText( asString( m_sourceModel->data(index) ) );
        }
      };//fieldUpdate
      
      if( fitAct )
      {
        const int row = table->rowCount();
        new WLabel( nuc->symbol + " Activity", table->elementAt(row, 0) );
        WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
        value->setAutoComplete( false );
        setFieldValue( value, SourceFitModel::Columns::kTruthActivity );
        
        auto valueUpdate = [fieldUpdate,value](){
          fieldUpdate( value, SourceFitModel::Columns::kTruthActivity );
        };
        
        value->changed().connect( std::bind(valueUpdate) );
        value->enterPressed().connect( std::bind(valueUpdate) );
        
        
        WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
        tolerance->setAutoComplete( false );
        setFieldValue( tolerance, SourceFitModel::Columns::kTruthActivityTolerance );
        auto toleranceUpdate = [fieldUpdate,tolerance](){
          fieldUpdate( tolerance, SourceFitModel::Columns::kTruthActivityTolerance );
        };
        
        tolerance->changed().connect( std::bind(toleranceUpdate) );
        tolerance->enterPressed().connect( std::bind(toleranceUpdate) );
      }//if( fitAct )
      
      
      if( fitAge )
      {
        const int row = table->rowCount();
        new WLabel( nuc->symbol + " Age", table->elementAt(row, 0) );
        WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
        value->setAutoComplete( false );
        setFieldValue( value, SourceFitModel::Columns::kTruthAge );
        
        auto valueUpdate = [fieldUpdate,value](){
          fieldUpdate( value, SourceFitModel::Columns::kTruthAge );
        };
        
        value->changed().connect( std::bind(valueUpdate) );
        value->enterPressed().connect( std::bind(valueUpdate) );
        
        
        WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
        tolerance->setAutoComplete( false );
        setFieldValue( tolerance, SourceFitModel::Columns::kTruthAgeTolerance );
        auto toleranceUpdate = [fieldUpdate,tolerance](){
          fieldUpdate( tolerance, SourceFitModel::Columns::kTruthAgeTolerance );
        };
        
        tolerance->changed().connect( std::bind(toleranceUpdate) );
        tolerance->enterPressed().connect( std::bind(toleranceUpdate) );
      }//if( fitAge )
    }//for( int i = 0; i < nnuc; ++i )
    
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
      if( !select )
        continue;
      
      if( select->isGenericMaterial() )
      {
        if( select->fitArealDensity() )
        {
          const int row = table->rowCount();
          WLabel *label = new WLabel( "Areal Density", table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          value->setAutoComplete( false );
          WDoubleValidator *dblValidator = new WDoubleValidator( 0, 500, value );
          value->setValidator( dblValidator );
          value->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          label->setBuddy( value );
          
          if( select->truthAD )
            value->setText( std::to_string(*select->truthAD) );
          
          auto updateVal = [select,value](){
            double answer = 0;
            if( (stringstream(value->text().toUTF8()) >> answer) )
              select->truthAD = answer;
            else if( select->truthAD )
              value->setText( std::to_string(*select->truthAD) );
            else
              value->setText( "" );
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          tolerance->setAutoComplete( false );
          dblValidator = new WDoubleValidator( 0, 100, tolerance );
          tolerance->setValidator( dblValidator );
          tolerance->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          
          auto updateTolerance = [select,tolerance](){
            double answer = 0;
            if( (stringstream(tolerance->text().toUTF8()) >> answer) )
              select->truthADTolerance = answer;
            else if( select->truthADTolerance )
              tolerance->setText( std::to_string(*select->truthADTolerance) );
            else
              tolerance->setText( "" );
          };//updateVal(...)
          
          if( select->truthADTolerance )
            tolerance->setText( std::to_string(*select->truthADTolerance) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        }//if( fit AD )
        
        if( select->fitAtomicNumber() )
        {
          const int row = table->rowCount();
          WLabel *label = new WLabel( "Atomic Number", table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          value->setAutoComplete( false );
          WDoubleValidator *dblValidator = new WDoubleValidator( MassAttenuation::sm_min_xs_atomic_number,
                                                                 MassAttenuation::sm_max_xs_atomic_number, value );
          value->setValidator( dblValidator );
          value->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          label->setBuddy( value );
          
          if( select->truthAN )
            value->setText( std::to_string(*select->truthAN) );
          
          auto updateVal = [select,value](){
            double answer = 0;
            if( (stringstream(value->text().toUTF8()) >> answer) )
              select->truthAN = answer;
            else if( select->truthAN )
              value->setText( std::to_string(*select->truthAN) );
            else
              value->setText( "" );
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          tolerance->setAutoComplete( false );
          dblValidator = new WDoubleValidator( 0, 100, tolerance );
          tolerance->setValidator( dblValidator );
          tolerance->addStyleClass( "numberValidator"); //used to detect mobile keyboard
          
          auto updateTolerance = [select,tolerance](){
            double answer = 0;
            const string txt = tolerance->text().toUTF8();
            if( txt.empty() )
            {
              select->truthANTolerance.reset();
            }else if( (stringstream(txt) >> answer) )
            {
              select->truthANTolerance = answer;
            }else if( select->truthANTolerance )
            {
              tolerance->setText( std::to_string(*select->truthANTolerance) );
            }else
            {
              select->truthANTolerance.reset();
              tolerance->setText( "" );
            }
          };//updateVal(...)
          
          if( select->truthANTolerance )
            tolerance->setText( std::to_string(*select->truthANTolerance) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        }//if( fit AN )
      }else
      {
        shared_ptr<Material> mat = select->material();
        if( !mat )
          throw runtime_error( "There is a non-generic material that is blank" );
  
        if( select->fitForMassFractions() )
          throw runtime_error( "A shielding fits for mass-fractions is not implemented yet" );
        
//        vector<const SandiaDecay::Nuclide *> srcnucs = select->nuclidesToUseAsSources();
//        if( !srcnucs.empty() )
//          throw runtime_error( "A shieldings used as sources is not yet implemented" );
        
        if( select->fitThickness() )
        {
          const int row = table->rowCount();
          WLabel *label = new WLabel( "Thickness", table->elementAt(row, 0) );
          WLineEdit *value = new WLineEdit( table->elementAt(row, 1) );
          value->setAutoComplete( false );
          WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, value );
          validator->setFlags( Wt::MatchCaseInsensitive );
          value->setValidator( validator );
          label->setBuddy( value );
          
          if( select->truthThickness )
            value->setText( PhysicalUnits::printToBestLengthUnits( *select->truthThickness, 2 ) );
          
          auto updateVal = [select,value](){
            const string txt = value->text().toUTF8();
            if( txt.empty() )
            {
              select->truthThickness.reset();
              return;
            }
            
            try
            {
              select->truthThickness = PhysicalUnits::stringToDistance( txt );
            }catch( ... )
            {
              if( select->truthThickness )
                value->setText( PhysicalUnits::printToBestLengthUnits(*select->truthThickness,2) );
              else
                value->setText( "" );
            }//try / catch
          };//updateVal(...)
          
          value->changed().connect( std::bind(updateVal) );
          value->enterPressed().connect( std::bind(updateVal) );
          
          WLineEdit *tolerance = new WLineEdit( table->elementAt(row, 2) );
          validator = new WRegExpValidator( PhysicalUnits::sm_distanceUncertaintyUnitsOptionalRegex, tolerance );
          validator->setFlags( Wt::MatchCaseInsensitive );
          tolerance->setValidator( validator );
          
          auto updateTolerance = [select,tolerance](){
            const string txt = tolerance->text().toUTF8();
            if( txt.empty() )
            {
              select->truthThicknessTolerance.reset();
              return;
            }
            
            try
            {
              select->truthThicknessTolerance = PhysicalUnits::stringToDistance( txt );
            }catch( ... )
            {
              if( select->truthThicknessTolerance )
                tolerance->setText( PhysicalUnits::printToBestLengthUnits(*select->truthThicknessTolerance,2) );
              else
                tolerance->setText( "" );
            }//try / catch
          };//updateVal(...)
          
          if( select->truthThicknessTolerance )
            tolerance->setText( PhysicalUnits::printToBestLengthUnits(*select->truthThicknessTolerance) );
          
          tolerance->changed().connect( std::bind(updateTolerance) );
          tolerance->enterPressed().connect( std::bind(updateTolerance) );
        }//if( fit thickness )
      }//if( generic material ) / else
    }//for( WWidget *widget : m_shieldingSelects->children() )
  }catch( std::exception &e )
  {
    contents->clear();
    WText *txt = new WText( e.what() , contents );
    txt->setInline( false );
  }//try / catch

  
  WPushButton *button = window->addCloseButtonToFooter("Okay");
  button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    
  window->centerWindow();
  window->disableCollapse();
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->resizeToFitOnScreen();
  window->show();
}//showInputTruthValuesWindow()


void ShieldingSourceDisplay::setFitQuantitiesToDefaultValues()
{
  const int nnuc = m_sourceModel->numNuclides();
  for( int i = 0; i < nnuc; ++i )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
    assert( nuc );
    
    const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
    const bool selfAttNuc = m_sourceModel->shieldingDeterminedActivity( i );
    
    // For self-attenuating shieldings, we'll just test the shielding thickness
    // For nuclides whose age is controlled by another nuclide, we dont need to test age.
    if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
      continue;
    
    if( m_sourceModel->fitActivity(i) )
    {
      WModelIndex index = m_sourceModel->index( i, SourceFitModel::kActivity );
      m_sourceModel->setData( index, "1 mCi" );
    }//if( fit activity )
    
    if( m_sourceModel->fitAge(i) )
    {
      string agestr;
      PeakDef::defaultDecayTime( nuc, &agestr );
      WModelIndex index = m_sourceModel->index( i, SourceFitModel::kAge );
      m_sourceModel->setData( index, agestr );
    }//if( fit age )
  }//for( int i = 0; i < nnuc; ++i )
  
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( !select )
      continue;
    
    if( select->isGenericMaterial() )
    {
      if( select->fitArealDensity() )
        select->arealDensityEdit()->setText( "0" );
      if( select->fitAtomicNumber() )
        select->atomicNumberEdit()->setText( "26" );
    }else
    {
      shared_ptr<Material> mat = select->material();
      if( !mat || select->fitForMassFractions() )
        continue;
      if( select->fitThickness() )
        select->thicknessEdit()->setText( "1 cm" );
    }//if( generic material ) / else
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void setFitQuantitiesToDefaultValues()



std::tuple<int,int,bool> ShieldingSourceDisplay::numTruthValuesForFitValues()
{
  bool isValid = true;
  int nFitQuantities = 0, nQuantitiesCan = 0;
  
  const int nnuc = m_sourceModel->numNuclides();
  for( int i = 0; i < nnuc; ++i )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
    assert( nuc );
    
    const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
    const bool selfAttNuc = m_sourceModel->shieldingDeterminedActivity( i );
    
    // For self-attenuating shieldings, we'll just test the shielding thickness
    // For nuclides whose age is controlled by another nuclide, we dont need to test age.
    if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
      continue;
    
    if( m_sourceModel->fitActivity(i) )
    {
      const boost::optional<double> activity = m_sourceModel->truthActivity(i);
      const boost::optional<double> tolerance = m_sourceModel->truthActivityTolerance(i);
      nFitQuantities += 1;
      nQuantitiesCan += (activity && tolerance);
      if( !(activity && tolerance) )
      {
        auto actindex = m_sourceModel->index( i, SourceFitModel::kTruthActivity );
        auto tolindex = m_sourceModel->index( i, SourceFitModel::kTruthActivityTolerance );
        
        cerr << "Dont have: (activity && tolerance): (" << !!activity << " && " << !!tolerance << ") -> via data -> ("
        << asString(actindex.data()) << ", " << asString(tolindex.data()) << ")" << endl;
      }
    }//if( fit activity )
    
    if( m_sourceModel->fitAge(i) )
    {
      const boost::optional<double> age = m_sourceModel->truthAge(i);
      const boost::optional<double> tolerance = m_sourceModel->truthAgeTolerance(i);
      nFitQuantities += 1;
      nQuantitiesCan += (age && tolerance);
      if( !(age && tolerance) )
        cerr << "Dont have: (age && tolerance)" << endl;
    }//if( fit age )
  }//for( int i = 0; i < nnuc; ++i )
  
  
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
    if( !select )
      continue;
    
    if( select->isGenericMaterial() )
    {
      if( select->fitArealDensity() )
      {
        if( select->truthAD && select->truthADTolerance )
          nQuantitiesCan += 1;
        nFitQuantities += 1;
      }//if( fit AD )
      
      if( select->fitAtomicNumber() )
      {
        if( select->truthAN && select->truthANTolerance )
          nQuantitiesCan += 1;
        nFitQuantities += 1;
      }//if( fit AN )
    }else
    {
      shared_ptr<Material> mat = select->material();
      if( !mat )
      {
        cerr << "Dont have: Coultn get material" << endl;
        isValid = false;
        continue;
      }
      
      if( select->fitForMassFractions() )
      {
        cerr << "Dont have: Fitting for mass fraction" << endl;
        isValid = false;
        continue;
      }
      
      if( select->fitThickness() )
      {
        if( select->truthThickness && select->truthThicknessTolerance )
          nQuantitiesCan += 1;
        nFitQuantities += 1;
      }//if( fit thickness )
    }//if( generic material ) / else
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  if( nQuantitiesCan != nFitQuantities )
  {
    cerr << "Dont have: nQuantitiesCan != nFitQuantities" << nQuantitiesCan << " != " << nFitQuantities << endl;
    isValid = false;
  }
  
  if( !nQuantitiesCan )
  {
    cerr << "Dont have: !nQuantitiesCan: " << nQuantitiesCan << endl;
    isValid = false;
  }
  
  return std::tuple<int,int,bool>( nQuantitiesCan, nFitQuantities, isValid );
}//bool haveTruthValuesForAllFitValues()


void ShieldingSourceDisplay::renderChi2Chart( Wt::WSvgImage &image )
{
  WPainter p( &image );
  m_chi2Graphic->paint( p );
  p.end();
}//void renderChi2Chart( Wt::WSvgImage &image );


tuple<bool,int,int,vector<string>> ShieldingSourceDisplay::testCurrentFitAgainstTruth()
{
  bool successful = true;
  int numCorrect = 0, numTested = 0;
  vector<string> textInfoLines;

  try
  {
    
    const int nnuc = m_sourceModel->numNuclides();
    for( int i = 0; i < nnuc; ++i )
    {
      const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( i );
      assert( nuc );
      
      const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuc );
      const bool selfAttNuc = m_sourceModel->shieldingDeterminedActivity( i );
      
      // For self-attenuating shieldings, we'll just test the shielding thickness
      // For nuclides whose age is controlled by another nuclide, we dont need to test age.
      if( selfAttNuc || (ageNuc && (ageNuc != nuc)) )
        continue;
      
      if( m_sourceModel->fitActivity(i) )
      {
        boost::optional<double> truthAct = m_sourceModel->truthActivity(i);
        boost::optional<double> tolerance = m_sourceModel->truthActivityTolerance(i);
        
        if( !truthAct || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " activity." );
          continue;
        }
        
        const double fitAct = m_sourceModel->activity(i);
        const bool closeEnough = (fabs(*truthAct - fitAct) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit activity "
                                + PhysicalUnits::printToBestActivityUnits(fitAct) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestActivityUnits(*truthAct)
                                + " and tolerance "
                                + PhysicalUnits::printToBestActivityUnits(*tolerance)
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
      }//if( fit activity )
      
      if( m_sourceModel->fitAge(i) )
      {
        const boost::optional<double> truthAge = m_sourceModel->truthAge(i);
        const boost::optional<double> tolerance = m_sourceModel->truthAgeTolerance(i);
        
        if( !truthAge || !tolerance )
        {
          successful = false;
          textInfoLines.push_back( "Did not have truth value for " + nuc->symbol + " age." );
          continue;
        }
        
        const double fitAge = m_sourceModel->age(i);
        const bool closeEnough = (fabs(*truthAge - fitAge) < *tolerance);
        
        numTested += 1;
        numCorrect += closeEnough;
        
        textInfoLines.push_back( "For " + nuc->symbol + " fit age "
                                + PhysicalUnits::printToBestTimeUnits(fitAge) + " with the"
                                " truth value of "
                                + PhysicalUnits::printToBestTimeUnits(*truthAge)
                                + " and tolerance "
                                + PhysicalUnits::printToBestTimeUnits(*tolerance)
                                + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                );
      }//if( fit age )
    }//for( int i = 0; i < nnuc; ++i )
    
    
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *select = dynamic_cast<ShieldingSelect *>(widget);
      if( !select )
        continue;
      
      if( select->isGenericMaterial() )
      {
        if( select->fitArealDensity() )
        {
          if( !select->truthAD || !select->truthADTolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AD for generic shielding" );
            continue;
          }
          
          const double fitAD = select->arealDensity();
          const bool closeEnough = (fabs(*select->truthAD - fitAD) < *select->truthADTolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAD)
                                  + " with the truth value of " + std::to_string(*select->truthAD)
                                  + " and tolerance " + std::to_string(*select->truthADTolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        }//if( fit AD )
        
        if( select->fitAtomicNumber() )
        {
          if( !select->truthAN || !select->truthANTolerance )
          {
            successful = false;
            textInfoLines.push_back( "Did not have truth AN for generic shielding" );
            continue;
          }
          
          const double fitAN = select->atomicNumber();
          const bool closeEnough = (fabs(*select->truthAN - fitAN) < *select->truthANTolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For Generic Shielding fit AN " + std::to_string(fitAN)
                                  + " with the truth value of " + std::to_string(*select->truthAN)
                                  + " and tolerance " + std::to_string(*select->truthANTolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        }//if( fit AN )
      }else
      {
        shared_ptr<Material> mat = select->material();
        if( !mat )
        {
          successful = false;
          textInfoLines.push_back( "There was an invalid material." );
          continue;
        }
        
        if( select->fitForMassFractions() )
        {
          successful = false;
          textInfoLines.push_back( "Mass fraction is being fit for, which isnt allowed." );
          continue;
        }
        
        if( select->fitThickness() )
        {
          if( !select->truthThickness || !select->truthThicknessTolerance )
          {
            successful = false;
            textInfoLines.push_back( "Missing truth thickness for shielding '" + mat->name + "'" );
            continue;
          }
          
          const double fitThickness = select->thickness();
          const bool closeEnough = (fabs(*select->truthThickness - fitThickness) < *select->truthThicknessTolerance);
          
          numTested += 1;
          numCorrect += closeEnough;
          
          textInfoLines.push_back( "For shielding '" + mat->name + "' fit thickness "
                                  + PhysicalUnits::printToBestLengthUnits(fitThickness,4)
                                  + " with the truth value of "
                                  + PhysicalUnits::printToBestLengthUnits(*select->truthThickness,4)
                                  + " and tolerance "
                                  + PhysicalUnits::printToBestLengthUnits(*select->truthThicknessTolerance)
                                  + (closeEnough ? " - within tolerance." : " - out of tolerance." )
                                  );
        }//if( fit thickness )
      }//if( generic material ) / else
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    successful = (successful && numTested);
  }catch( std::exception &e )
  {
    successful = false;
    textInfoLines.push_back( "Caught exception during testing: " + string(e.what()) );
  }//try / catch
  
  return tuple<bool,int,int,vector<string>>( successful, numCorrect, numTested, textInfoLines );
}//std::tuple<bool,int,int,std::vector<std::string>> testCurrentFitAgainstTruth();



#endif //INCLUDE_ANALYSIS_TEST_SUITE



void ShieldingSourceDisplay::multiNucsPerPeakChanged()
{
  updateChi2Chart();
}//void multiNucsPerPeakChanged()


void ShieldingSourceDisplay::backgroundPeakSubChanged()
{
  if( m_backgroundPeakSub->isChecked() )
  {
    std::shared_ptr<const SpecMeas> back = m_specViewer->measurment(SpecUtils::SpectrumType::Background);
    
    if( !back )
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( "There is no background spectrum loaded, so can not"
                   " subtract background peak areas from foreground",
                   "", WarningWidget::WarningMsgHigh );
      return;
    }//if( !back )
    
    std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks;
    
    const set<int> &displayed = m_specViewer->displayedSamples(SpecUtils::SpectrumType::Background);
    peaks = back->peaks( displayed );
    
    if( !peaks || peaks->empty() )
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( "The background spectrum does not have any peaks"
                   " identified. You should load the background spectrum as"
                   " the foreground, identify the peaks to potentially be"
                   " subtracted, and then switch back to your current"
                   " foreground/background combination.",
                   "", WarningWidget::WarningMsgHigh );
      return;
    }//if( !peaks || peaks->empty() )
  }//if( m_backgroundPeakSub->isChecked() )
  
  updateChi2Chart();
}//void backgroundPeakSubChanged()


void ShieldingSourceDisplay::sameIsotopesAgeChanged()
{
  m_sourceModel->setUseSameAgeForIsotopes( m_sameIsotopesAge->isChecked() );
  updateChi2Chart();
}//void sameIsotopesAgeChanged()

void ShieldingSourceDisplay::showGraphicTypeChanged()
{
  const bool chi = m_showChiOnChart->isChecked();
  m_chi2Graphic->setShowChiOnChart( chi );
  updateChi2Chart();
}

void ShieldingSourceDisplay::handleUserDistanceChange()
{
  string distanceStr = m_distanceEdit->text().toUTF8();
  
  //Default to cm if no distance is given
  SpecUtils::trim( distanceStr );
  if( distanceStr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
  {
    distanceStr += " cm";
    m_distanceEdit->setText( distanceStr );
  }
  
  try
  {
    const double distance = PhysicalUnits::stringToDistance( distanceStr );
    if( distance <= 0.0 )
      throw runtime_error( "Must have a non-zero distance" );
    m_prevDistStr = distanceStr;
  }catch(...)
  {
    m_distanceEdit->setText( m_prevDistStr );
  }
  
  
  updateChi2Chart();
}//void ShieldingSourceDisplay::handleUserDistanceChange()


void ShieldingSourceDisplay::checkDistanceAndThicknessConsistent()
{
  double shieldrad = 0.0, distance = 0.0;
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( thisSelect && !thisSelect->isGenericMaterial() )
      shieldrad += thisSelect->thickness();
  }//for( WWidget *widget : m_shieldingSelects->children() )
  
  const string distanceStr = m_distanceEdit->text().toUTF8();
  distance = PhysicalUnits::stringToDistance( distanceStr );
  if( distance <= 0.0 )
    throw runtime_error( "Distance must be greater than zero" );
  
  if( shieldrad > distance )
  {
    const double scale = 0.95*distance/shieldrad;
    for( WWidget *widget : m_shieldingSelects->children() )
    {
      ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
      if( thisSelect && !thisSelect->isGenericMaterial() )
        thisSelect->thicknessEdit()->setText(
         PhysicalUnits::printToBestLengthUnits(scale*thisSelect->thickness()) );
    }//for( WWidget *widget : m_shieldingSelects->children() )
    
    throw runtime_error( "Shielding thicknesses have been scaled to be less "
                         " than detector distance.");
  }//if( shieldrad < distance )
}//void checkDistanceAndThicknessConsistent()


void ShieldingSourceDisplay::checkAndWarnZeroMassFraction()
{
  //ToDo: will needs to implement, the point of this function will be to allow
  //  making sure the user hasnt selected a shielding to be a self attenuating
  //  source, but only allowed the source isotopes to be a mass fraction of zero.
  //  Should through std::exception with a use-friendly descriptive message on
  //  error.
  
  cerr << "ShieldingSourceDisplay::checkAndWarnZeroMassFraction() unimplemented!" << endl;
  
  
  
}//void checkAndWarnZeroMassFraction()


void ShieldingSourceDisplay::updateChi2Chart()
{
  m_chi2ChartNeedsUpdating = true;
  scheduleRender(); //trigger re-render
}//void updateChi2Chart()


void ShieldingSourceDisplay::updateChi2ChartActual()
{
  try
  {
    checkDistanceAndThicknessConsistent();
  }catch( exception &e )
  {
    passMessage( e.what() + string("<br />Chi2 chart not updated."),
                 "", WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
  try
  {
    vector<ShieldingSelect *> shieldings;
    ROOT::Minuit2::MnUserParameters inputPrams;

    Chi2FcnShrdPtr chi2Fcn = shieldingFitnessFcn( shieldings, inputPrams );
    const unsigned int ndof = inputPrams.VariableParameters();
    
    const vector<double> params = inputPrams.Params();
    const vector<double> errors = inputPrams.Errors();
    GammaInteractionCalc::PointSourceShieldingChi2Fcn::NucMixtureCache mixcache;
    
    
    m_calcLog.clear();
    if( m_logDiv )
    {
      m_logDiv->contents()->clear();
      m_logDiv->hide();
    }//if( m_logDiv )
    
    const vector< tuple<double,double,double,Wt::WColor,double> > chis
                      = chi2Fcn->energy_chi_contributions( params, mixcache,
                                                m_multiIsoPerPeak->isChecked(),
                                                &m_calcLog );
    
    m_showLog->setDisabled( m_calcLog.empty() );

    typedef tuple<double,double,double,Wt::WColor,double> DDPair;
    vector< DDPair > keeper_points;

    for( size_t row = 0; row < chis.size(); ++row )
    {
      const double energy = std::get<0>(chis[row]);
      const double chi = std::get<1>(chis[row]);
      const double scale = std::get<2>(chis[row]);
      const WColor &color = std::get<3>(chis[row]);
      const double scale_uncert = std::get<4>(chis[row]);

      if( fabs(chi) < 1.0E5 && !IsInf(chi) && !IsNan(chi)
          && !IsInf(energy) && !IsNan(energy) )
        keeper_points.push_back( DDPair(energy,chi,scale,color,scale_uncert) );
    }//for( size_t row = 0; row < chis.size(); ++row )

    //If we only have one point, then Wt cant find the y-range, and an assert
    //  statment gets triggered in WAxis.C, so we'll fix this up kinda.
    //  I should probably submit a bug report to Wt...
    if( keeper_points.size() == 1 )
      keeper_points.insert( keeper_points.begin(), DDPair(0.0,0.0,0.0,WColor(),0.0) );

    m_chi2Graphic->setNumFitForParams( ndof );
    if( !m_calcLog.empty() )
    {
      char buffer[64];
      snprintf( buffer, sizeof(buffer), "There %s %i parameter%s fit for",
                (ndof>1 ? "were" : "was"), int(ndof), (ndof>1 ? "s" : "") );
      m_calcLog.push_back( "&nbsp;" );
      m_calcLog.push_back( buffer );
    }//if( !m_calcLog.empty() )
    
    const int nrow = static_cast<int>( keeper_points.size() );
    const int nStartRows = m_chi2Model->rowCount();
    if( nStartRows < nrow )
      m_chi2Model->insertRows( nStartRows, nrow - nStartRows );

    if( nStartRows > nrow )
      m_chi2Model->removeRows( nrow, nStartRows - nrow );

    std::shared_ptr<const deque< PeakModel::PeakShrdPtr > > peaks = m_peakModel->peaks();
    
    for( int row = 0; row < nrow; ++row  )
    {
      const DDPair &p = keeper_points[row];
      
      const double &energy = std::get<0>(p);
      const double &chi = std::get<1>(p);
      const double &scale = std::get<2>(p);
      WColor color = std::get<3>(p);
      const double &scale_uncert = std::get<4>(p);
      
      if( IsNan(energy) || IsInf(chi) )
      {
        stringstream msg;
        msg << "An invalid chi2 was calculated for " << energy
            << " keV, other results may be invalid";
        passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
        continue;
      }//if( IsNan(p.second) || IsInf(p.second) )
      
      m_chi2Model->setData( row, 0, boost::any(energy) );
      m_chi2Model->setData( row, 1, boost::any(chi) );
      m_chi2Model->setData( row, 2, boost::any(scale) );
      
      if( color.isDefault() )
        color = m_specViewer->getColorTheme()->defaultPeakLine;
      color.setRgb( color.red(), color.green(), color.blue(), 255 );
      
      m_chi2Model->setData( row, 1, boost::any(color), Wt::MarkerPenColorRole );
      m_chi2Model->setData( row, 1, boost::any(color), Wt::MarkerBrushColorRole );
      m_chi2Model->setData( row, 2, boost::any(color), Wt::MarkerPenColorRole );
      m_chi2Model->setData( row, 2, boost::any(color), Wt::MarkerBrushColorRole );
      
      //If we wanted to include the nuclide in the model, we would have to loop
      //  over photopeaks in m_peakModel to try and match things up
      WString nuclidename;
      if( !!peaks )
      {
        for( const PeakModel::PeakShrdPtr &peak : *peaks )
        {
          if( peak->useForShieldingSourceFit()
              && peak->parentNuclide() && peak->decayParticle() )
          {
            if( fabs(energy - peak->decayParticle()->energy) < 0.001 )
            {
              nuclidename = peak->parentNuclide()->symbol;
              break;
            }
          }
        }//for( const PeakModel::PeakShrdPtr &p : *peaks )
        
        m_chi2Model->setData( row, 3, boost::any(nuclidename) );
        m_chi2Model->setData( row, 4, boost::any(scale_uncert) );
      }//if( !!peaks )
    }//for( int row = 0; row < nrow; ++row  )
  }catch( std::exception &e )
  {
    //One reason we may have made it here is if there are no peaks selected
    //Another is if we have removed using a peak from the peak, and it was an
    //  age defining, and we havent updated the dependant ages yet (but will in this
    //  event loop).
    if( m_chi2Model->rowCount() )
      m_chi2Model->removeRows( 0, m_chi2Model->rowCount() );
    cerr << "ShieldingSourceDisplay::updateChi2ChartActual()\n\tCaught:" << e.what() << endl;
  }
  
  
  const std::vector<WAbstractArea *> oldareas = m_chi2Graphic->areas();
  for( WAbstractArea *a : oldareas )
  {
    m_chi2Graphic->removeArea( a );
    delete a;
  }

//  cerr << "m_chi2Model->rowCount()=" << m_chi2Model->rowCount()
//       << "m_chi2Model->columnCount()=" << m_chi2Model->columnCount() << endl;

  m_calcLog.push_back( ns_no_uncert_info_txt );
}//void ShieldingSourceDisplay::updateChi2ChartActual()


void ShieldingSourceDisplay::showCalcLog()
{
  if( !m_logDiv )
  {
    m_logDiv = new AuxWindow( "Calculation Log" );
    m_logDiv->contents()->addStyleClass( "CalculationLog" );
    m_logDiv->disableCollapse();
    m_logDiv->rejectWhenEscapePressed();
    m_logDiv->setResizable( true );
  }//if( !m_logDiv )
  
  m_logDiv->contents()->clear();
  for( const string &str : m_calcLog )
    (new WText( str, m_logDiv->contents() ))->setInline( false );
  
  m_logDiv->show();
  
  const int wwidth = m_specViewer->renderedWidth();
  const int wheight = m_specViewer->renderedHeight();
  m_logDiv->setMaximumSize( 0.9*wwidth, wheight );
  m_logDiv->resizeToFitOnScreen();
  m_logDiv->centerWindow();
}//void showCalcLog()


void ShieldingSourceDisplay::finishModelUpload( AuxWindow *window,
                                                WFileUpload *upload )
{
  rapidxml::xml_document<char> original_doc;
  
  try
  {
    serialize( &original_doc );
    
    const std::string filename = upload->spoolFileName();
    
    std::vector<char> data;
    SpecUtils::load_file_data( filename.c_str(), data );
    
    rapidxml::xml_document<char> new_doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;
    new_doc.parse<flags>( &data.front() );
    deSerialize( new_doc.first_node() );
    
    m_modifiedThisForeground = true;
  }catch( std::exception &e )
  {
    stringstream msg;
    msg << "Error opening uploaded Source Shielding Fit Model: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    
    try
    {
      deSerialize( &original_doc );
    }catch( std::exception & )
    {
      passMessage( "Even worse, there was an error trying to recover",
                   "", WarningWidget::WarningMsgHigh );
    }//try / catch
  }//try / catch
  
  delete window;
}//void startModelUpload()


void ShieldingSourceDisplay::modelUploadError( const ::int64_t size_tried,
                                               AuxWindow *window )
{
  stringstream msg;
  const int max_size = static_cast<int>( wApp->maximumRequestSize() );
  msg << "Error uploading Source Shielding Fit Model.  Tried to upload "
      << size_tried << " (max size " << max_size << ")";
  passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
  
  AuxWindow::deleteAuxWindow( window );
}//void modelUploadError( const ::int64_t size_tried );


void ShieldingSourceDisplay::startModelUpload()
{
  AuxWindow *window = new AuxWindow( "Import Source Shielding XML Model",
                      (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                        | AuxWindowProperties::TabletModal) );
  
  WContainerWidget *contents = window->contents();
  WFileUpload *upload = new WFileUpload( contents );
  upload->setInline( false );
  
  upload->uploaded().connect( boost::bind( &ShieldingSourceDisplay::finishModelUpload, this, window, upload ) );
  upload->fileTooLarge().connect( boost::bind( &ShieldingSourceDisplay::modelUploadError, this, _1, window ) );
  upload->changed().connect( upload, &WFileUpload::upload );
  
    
    
  WPushButton *button = window->addCloseButtonToFooter("Cancel");
  button->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  window->centerWindow();
  window->disableCollapse();
//  window->hidden().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );

  window->show();
}//void startModelUpload()

#if( USE_DB_TO_STORE_SPECTRA )
typedef Dbo::QueryModel< Dbo::ptr<ShieldingSourceModel> > QueryModel_t;

void updateDescription( WSelectionBox *select,
                        WSelectionBox *other_select,
                        WTextArea *summary,
                        WPushButton *button,
                        InterSpec *specViewer )
{
  button->disable();
  summary->setValueText( "" );
  if( other_select )
    other_select->clearSelection();
  
  const int row = select->currentIndex();
  if( row < 0 )
    return;
  
  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( select->model() );
  if( !querymodel )
    throw runtime_error( "updateDescription(...): invalid input" );
  
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel )
  {
    cerr << "updateDescription(): Unable to load selected model - sorry :(\n";
    return;
  }//if( !shieldmodel )
  
  WString descrip = shieldmodel->description;
  
  //XXX - I dont understand why the date isnt always valid, unless there is a
  //  bug in Wt with serializing dates to the database somehow.
  if( shieldmodel->creationTime.isValid() )
  {
    descrip += (WString("\nCreated ")
               + shieldmodel->creationTime.toString(DATE_TIME_FORMAT_STR));
    if( abs(shieldmodel->creationTime.secsTo( shieldmodel->serializeTime )) > 30 )
      descrip += (WString(", Saved ")
                 + shieldmodel->serializeTime.toString(DATE_TIME_FORMAT_STR));
  }//if( shieldmodel->creationTime.isValid() )
  
  
  {//begin interact with database
    std::shared_ptr<DataBaseUtils::DbSession> sql = specViewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    const Dbo::collection<Dbo::ptr<UserFileInDb> > &files
                                                  = shieldmodel->filesUsedWith;
    if( files.size() )
    {
      descrip += "\nSaved while working with spectra:";
      for( Dbo::collection<Dbo::ptr<UserFileInDb> >::const_iterator iter = files.begin();
           iter != files.end(); ++iter )
      {
        descrip += "\n    " + (*iter)->filename;
//                + " (" + (*iter)->serializeTime.toString(DATE_TIME_FORMAT_STR)
//                + ")";
      }//for( loop over files )
    }//if( files.size() )
    transaction.commit();
  }//end interact with database
  
  summary->setValueText( descrip );
  
  button->enable();
}//updateDescription(...)


void ShieldingSourceDisplay::removeModelFromDb( WSelectionBox *selec1,
                                                WSelectionBox *selec2 )
{
  WSelectionBox *selec = selec1;
  int row = selec ? selec->currentIndex() : -1;
  if( row < 0 )
  {
    selec = selec2;
    row = selec ? selec->currentIndex() : -1;
  }

  if( !selec || row < 0 )
    return;

  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( selec->model() );  
  if( !querymodel )
    throw runtime_error( "removeModelFromDb(...): invalid input" );
  
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel || !shieldmodel.session() )
  {
    cerr << "Error removeModelFromDb: shieldmodel=" << shieldmodel << " and "
         << " session " << shieldmodel.session() << endl;
    return;
  }

  {
    std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
    DataBaseUtils::DbTransaction transaction( *sql );
    shieldmodel.remove();
    transaction.commit();
    querymodel->reload();
  }
  
  
  if( selec1 )
    querymodel = dynamic_cast<QueryModel_t *>( selec1->model() );
  if( querymodel )
    querymodel->reload();
  if( selec2 )
    querymodel = dynamic_cast<QueryModel_t *>( selec2->model() );
  if( querymodel )
    querymodel->reload();
}//void ShieldingSourceDisplay::removeModelFromDb( WSelectionBox *selec )

bool ShieldingSourceDisplay::loadModelFromDb( Dbo::ptr<ShieldingSourceModel> shieldmodel )
{
  if( !shieldmodel )
    return false;
  
  rapidxml::xml_document<char> original_doc;
  try
  {
    serialize( &original_doc );
    rapidxml::xml_document<char> new_doc;
    const int flags = rapidxml::parse_normalize_whitespace
                      | rapidxml::parse_trim_whitespace;
    string data = shieldmodel->xmlData;
    if( data.size() )
      new_doc.parse<flags>( &(data[0]) );
    deSerialize( new_doc.first_node() );
  }catch( std::exception &e )
  {
    stringstream msg;
    msg << "Error opening Source Shielding Fit Model from database: " << e.what();
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    
    try
    {
      deSerialize( &original_doc );
    }catch( std::exception & )
    {
      passMessage( "Even worse, there was an error trying to recover",
                  "", WarningWidget::WarningMsgHigh );
    }//try / catch
    
    return false;
  }//try / catch
  
  return true;
}//void loadModelFromDb( Wt::Dbo::ptr<ShieldingSourceModel> entry )

void ShieldingSourceDisplay::finishLoadModelFromDatabase( AuxWindow *window,
                                                  WSelectionBox *first_selct,
                                                  WSelectionBox *other_select )
{
  WSelectionBox *selec = first_selct;
  int row = selec ? selec->currentIndex() : -1;
  if( row < 0 )
  {
    selec = other_select;
    row = selec ? selec->currentIndex() : -1;
  }//if( row < 0 )
  
  if( !selec )
    throw runtime_error( "finishLoadModelFromDatabase(...): invalid slection box" );
  
  QueryModel_t *querymodel = dynamic_cast<QueryModel_t *>( selec->model() );
  
  if( row < 0 )
    return;
  if( !querymodel )
    throw runtime_error( "finishLoadModelFromDatabase(...): invalid input" );
    
  Dbo::ptr<ShieldingSourceModel> shieldmodel = querymodel->resultRow( row );
  if( !shieldmodel )
  {
    passMessage( "Unable to load selected model - sorry :(",
                 "", WarningWidget::WarningMsgHigh );
    return;
  }//if( !shieldmodel )
  
  loadModelFromDb( shieldmodel );
  
  delete window;
}//void finishLoadModelFromDatabase()


void ShieldingSourceDisplay::startBrowseDatabaseModels()
{
  if( !m_specViewer || !m_specViewer->m_user )
    throw runtime_error( "startBrowseDatabaseModels(): invalid user" );
  
  WTextArea *summary = NULL;
  WPushButton *accept = NULL, *cancel = NULL, *del = NULL;
  AuxWindow *window = new AuxWindow( "Previously Saved Models",
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal) | AuxWindowProperties::TabletModal) );
  
  try
  {
    WContainerWidget *contents = window->contents();
    summary = new WTextArea();
//    summary->setColumns( 30 );
    summary->setWidth( 316 );
    summary->setHeight( 50 );
    summary->disable();
    accept = new WPushButton( "Load" );
    accept->setIcon( "InterSpec_resources/images/database_go.png" );
    accept->disable();
  
    cancel = new WPushButton( "Cancel" );
    cancel->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  
    del = new WPushButton( "Delete" );
    del->setIcon( "InterSpec_resources/images/minus_min_white.png" );
    del->disable();

    Dbo::ptr<UserFileInDb> dbmeas;
    dbmeas = m_specViewer->measurmentFromDb( SpecUtils::SpectrumType::Foreground, false );
    
    size_t nfileprev[2];
    WSelectionBox *selections[2] = { (WSelectionBox *)0, (WSelectionBox *)0 };
    
    {//begin codeblock for database interaction
      std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
      DataBaseUtils::DbTransaction transaction( *sql );
      nfileprev[0] = dbmeas ? dbmeas->modelsUsedWith.size() : 0;
      nfileprev[1] = m_specViewer->m_user->shieldSrcModels().size();
      transaction.commit();
    }//end codeblock for database interaction
    
    
    for( size_t i = 0; i < 2; ++i )
    {
      if( !nfileprev[i] )
        continue;
      
     	WSelectionBox *selection = new WSelectionBox();
      selections[i] = selection;
      QueryModel_t *model = new QueryModel_t( selection );
      if( i == 0 )
        model->setQuery( dbmeas->modelsUsedWith.find() );
      else
        model->setQuery( m_specViewer->m_user->shieldSrcModels().find() );
      model->addColumn( "Name" );
      selection->setModel( model );
      selection->setModelColumn( 0 );
      selection->setHeight( 100 );
      selection->setWidth( 320 );
      selection->setMaximumSize( 320, 100 );
      selection->setSelectionMode( Wt::SingleSelection );
      
      const char *msg = ((i==0) ? "Models Previously used with this foreground:"
                                : "Models you've previously saved:");
      WText *title = new WText( msg, contents );
      title->setAttributeValue( "style", "font-weight:bold;" );
      title->setMargin( (i?8:12), Wt::Top );
      title->setInline( false );
      selection->setInline( false );
      selection->setMargin( 3, Wt::Top );
      selection->setMargin( 9, Wt::Left );
      contents->addWidget( selection );
    }//for( size_t i = 0; i < 2; ++i )
  
    if( selections[0] )
    {
      selections[0]->activated().connect( boost::bind( &updateDescription,
                selections[0], selections[1], summary, accept, m_specViewer ) );
      selections[0]->doubleClicked().connect(
            boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                         this, window, selections[0], selections[1] ) );
      accept->clicked().connect(
              boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                           this, window, selections[0], selections[1] ) );
    }//if( selections[0] )
    
    if( selections[1] )
    {
      selections[1]->activated().connect( boost::bind( &updateDescription,
                selections[1], selections[0], summary, accept, m_specViewer ) );
      selections[1]->doubleClicked().connect(
             boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                          this, window, selections[1], selections[0] ) );
      accept->clicked().connect(
             boost::bind( &ShieldingSourceDisplay::finishLoadModelFromDatabase,
                          this, window, selections[1], selections[0] ) );
    }//if( selections[1] )
    
    del->clicked().connect(
                        boost::bind( &ShieldingSourceDisplay::removeModelFromDb,
                                       this, selections[0], selections[1] ) );
    
    if( nfileprev[0] || nfileprev[1] )
    {
      WText *title = new WText( "Model Description:" );
      title->setAttributeValue( "style", "font-weight:bold;margin-top:12px;" );

      WCheckBox *cb = new WCheckBox( "Allow delete" );
      
      window->footer()->addWidget( cb );
      
      if( !m_specViewer->isMobile() )
      {
        cb->setFloatSide(Left);
        del->setFloatSide(Left);
      }
        
      window->footer()->addWidget( del );
        
      cb->checked().connect( del, &WPushButton::enable );
      cb->unChecked().connect( del, &WPushButton::disable );
      window->footer()->addWidget( cancel );

      window->footer()->addWidget( accept );

      title->setInline( false );
      summary->setInline( false );
      summary->setMargin( 3, Wt::Top );
      summary->setMargin( 9, Wt::Left );
      
      contents->addWidget( title );
      contents->addWidget( summary );
    }else
    {
      WText *info = new WText( "There are no models in the database for you" );
      info->setInline( false );
      contents->addWidget( info );
      contents->addWidget( cancel );
      
      if( accept )
      {
        delete accept;
        accept = NULL;
      }//if( accept )
      if( summary )
      {
        delete summary;
        summary = NULL;
      }
    }//if( nfileprev[0] || nfileprev[1] )
    
    window->rejectWhenEscapePressed();
    window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
    window->setWidth( 350 );
    window->disableCollapse();
    window->centerWindow();
    window->show();
  }catch( std::exception &e )
  {
    if( accept )
      delete accept;
    if( cancel )
      delete cancel;
    if( summary )
     delete summary;
    if( del )
      delete del;
    delete window;
    passMessage( "Error creating database model browser",
                 "", WarningWidget::WarningMsgHigh );
    cerr << "\n\nShieldingSourceDisplay::startBrowseDatabaseModels() caught: "
         << e.what() << endl << endl;
  }//try / catch
}//void ShieldingSourceDisplay::startBrowseDatabaseModels()


Wt::Dbo::ptr<ShieldingSourceModel> ShieldingSourceDisplay::modelInDb()
{
  saveModelIfAlreadyInDatabase();
  if( m_modelInDb )
    return m_modelInDb;
  
  finishSaveModelToDatabase( defaultModelName(), defaultModelDescription() );
  
  return m_modelInDb;
}//Wt::Dbo::ptr<ShieldingSourceModel> modelInDb()
#endif //#if( USE_DB_TO_STORE_SPECTRA )

std::string ShieldingSourceDisplay::defaultModelName() const
{
  std::shared_ptr<const SpecMeas> meas = m_specViewer->measurment( SpecUtils::SpectrumType::Foreground );
  string name;
  if( meas && !meas->filename().empty() )
  {
    name = "1D model for " + meas->filename() + " ";
    name +=  WDateTime::currentDateTime().toString( "yyyyMMdd" ).toUTF8();
    if( name.length() > 255 )
      name = name.substr( 0, 255 );
  }//if( meas && !meas->filename().empty() )
  return name;
}//defaultModelName()


std::string ShieldingSourceDisplay::defaultModelDescription() const
{
  string descrip;
  for( int row = 0; row < m_sourceModel->rowCount(); ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
    {
      if( row )
        descrip += ", ";
      descrip += nuc->symbol;
    }//if( nuc )
  }//for( int row = 0; row < m_sourceModel->rowCount(); ++row )
  
  descrip += "; shielded by ";
  
  int nshield = 0;
  for( WWidget *w : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( w );
    if( select && select->material() )
    {
      if( nshield )
        descrip += " and ";
      descrip += "'" + select->material()->name + "'";
      ++nshield;
    }//if( select )
  }//for( WWidget *w : m_shieldingSelects->children() )
  
  if( descrip.size() > 511 )
    descrip = descrip.substr( 0, 511 );
  
  return descrip;
}//std::string ShieldingSourceDisplay::defaultModelDescription() 

#if( USE_DB_TO_STORE_SPECTRA )
void ShieldingSourceDisplay::startSaveModelToDatabase( bool prompt )
{
  if( m_modelInDb && !prompt )
  {
    finishSaveModelToDatabase( NULL, NULL, NULL );
    return;
  }//if( m_modelInDb && !prompt )
  
  AuxWindow *window = new AuxWindow( "Import Source Shielding XML Model",
                  (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::IsAlwaysModal)
                   | AuxWindowProperties::TabletModal
                   | AuxWindowProperties::DisableCollapse) );
  WContainerWidget *contents = window->contents();
  window->centerWindow();
  
  window->rejectWhenEscapePressed();
  window->finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, window ) );
  window->show();
 
  WLabel *label = new WLabel( "Enter model name:", contents );
  label->setInline( false );

  WLineEdit *nameEdit = new WLineEdit( contents );
  nameEdit->setInline( false );
  if( m_modelInDb )
    nameEdit->setValueText( m_modelInDb->name );

  if( nameEdit->valueText().empty() )
    nameEdit->setValueText( defaultModelName() );
  
  label = new WLabel( "Enter model description (optional):", contents );
  label->setInline( false );
  
  WLineEdit *descEdit = new WLineEdit( contents );
  descEdit->setInline( false );
  if( m_modelInDb )
    descEdit->setValueText( m_modelInDb->description );
  
  if( descEdit->valueText().empty() )
    descEdit->setValueText( defaultModelDescription() );
  
  nameEdit->enterPressed().connect( boost::bind( &WFormWidget::setFocus, descEdit, true ) );
  
  descEdit->setTextSize( 32 );
  nameEdit->setTextSize( 32 );

 

  WPushButton *button = new WPushButton( "Save", window->footer() );
  button->setIcon( "InterSpec_resources/images/disk2.png" );
  
  button->clicked().connect(
              boost::bind( &ShieldingSourceDisplay::finishSaveModelToDatabase,
                           this, window, nameEdit, descEdit ) );
  descEdit->enterPressed().connect( boost::bind( &WFormWidget::setFocus, button, true ) );
}//void startSaveModelToDatabase()


void ShieldingSourceDisplay::finishSaveModelToDatabase( AuxWindow *window,
                                                        WLineEdit *name_edit,
                                                        WLineEdit *desc_edit )
{
  if( name_edit && name_edit->valueText().empty() )
  {
    WText *txt = new WText( "You must enter a name", window->contents() );
    txt->setInline( false );
    txt->setAttributeValue( "style", "color:red;" );
    return;
  }//if( name_edit && name_edit->valueText().empty() )
  
  //Check that the name is unique
  //  I'm not sure if I want to actually enforce this, or at a minimum should
  //  give a way to overide this
//  if( !m_modelInDb )
//  {
//    DbTransaction transaction( viewer );
//    const size_t nexisting = m_specViewer->m_user->m_shieldSrcModels.find()
//                            .where( "Name = ?" ).bind( name_edit->valueText() )
//                            .resultList().size();
//    transaction.commit();
//    
//    if( nexisting )
//    {
//      vector<WWidget *> children = window->contents()->children();
//      WText *txt = children.size() ? dynamic_cast<WText *>(children.back()) : 0;
//      if( !txt )
//      {
//        txt = new WText( window->contents() );
//        txt->setInline( false );
//        txt->setAttributeValue( "style", "color:red;" );
//      }//if( !txt )
//      txt->setText( "You already have a model with this name" );
//      return;
//    }//if( nexisting )
//  }//if( check if a unique name )
  
  WString name, description;
  if( name_edit )
    name = name_edit->valueText();
  if( desc_edit )
    description = desc_edit->valueText();
  
  finishSaveModelToDatabase( name, description );
  
  if( window )
    delete window;
}//finishSaveModelToDatabase(...)


bool ShieldingSourceDisplay::finishSaveModelToDatabase( const Wt::WString &name,
                                                       const Wt::WString &desc )
{
  ShieldingSourceModel *model = NULL;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );

  try
  {
    if( !m_modelInDb )
    {
      model = new ShieldingSourceModel();
      model->user = m_specViewer->m_user;
      model->serializeTime = WDateTime::currentDateTime();
      m_modelInDb.reset( new ShieldingSourceModel() );
      m_modelInDb = sql->session()->add( model );
    }//if( m_modelInDb ) / else
    
    model = m_modelInDb.modify();
    model->serializeTime = WDateTime::currentDateTime();
    
    if( !name.empty() )
      model->name = name;
    if( !desc.empty() )
      model->description = desc;
    
    const string utfname = model->name.toUTF8();
    const string utfdesc = model->description.toUTF8();
    if( utfname.size() > 255 )
      model->name.fromUTF8( utfname.substr( 0, 255 ) );
    if( utfdesc.size() > 511 )
      model->description.fromUTF8( utfdesc.substr( 0, 511 ) );
    
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    model->xmlData.clear();
    rapidxml::print(std::back_inserter(model->xmlData), doc, 0);

    Dbo::ptr<UserFileInDb> dbmeas;
    dbmeas = m_specViewer->measurmentFromDb( SpecUtils::SpectrumType::Foreground, true );
    
    if( dbmeas )
      model->filesUsedWith.insert( dbmeas );
    
    transaction.commit();
    m_saveAsNewModelInDb->enable();
  }catch( std::exception & )
  {
//    if( m_modelInDb.id() < 0 )
//    {
      m_modelInDb.reset();
      m_saveAsNewModelInDb->disable();
//    }//if( m_modelInDb.id() < 0 )
    transaction.rollback();
    return false;
  }//try / catch
  
  return true;
}//void finishSaveModelToDatabase( AxuWindow *window, Wt::WLineEdit *edit )


void ShieldingSourceDisplay::saveCloneModelToDatabase()
{
  if( !m_modelInDb )
    return;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );

  try
  {
    m_modelInDb.reread();
    ShieldingSourceModel *model = new ShieldingSourceModel();
    model->shallowEquals( *m_modelInDb );
    
    model->user = m_specViewer->m_user;
    model->serializeTime = WDateTime::currentDateTime();
    model->name = model->name + " Clone";
    
    const string utfname = model->name.toUTF8();
    if( utfname.size() > 255 )
      model->name.fromUTF8( utfname.substr(0,255) );
    
    m_modelInDb = sql->session()->add( model );
    
    Dbo::ptr<UserFileInDb> dbmeas;
    dbmeas = m_specViewer->measurmentFromDb( SpecUtils::SpectrumType::Foreground, false );
    if( dbmeas )
      m_modelInDb.modify()->filesUsedWith.insert( dbmeas );
    transaction.commit();
  }catch( std::exception &e )
  {
    m_modelInDb.reset();
    m_saveAsNewModelInDb->disable();
    cerr << "\n\nException caught in ShieldingSourceDisplay::saveCloneModelToDatabase(): "
         << e.what() << endl;
    transaction.rollback();
    passMessage( "Error saving to the database", "", WarningWidget::WarningMsgHigh );
    return;
  }//try / catch
  
  cerr << "\n\nFinishing in ShieldingSourceDisplay::saveCloneModelToDatabase()" << endl;
  startSaveModelToDatabase( true );
}//void saveCloneModelToDatabase()


void ShieldingSourceDisplay::saveModelIfAlreadyInDatabase()
{
  if( !m_modelInDb )
    return;
  
  std::shared_ptr<DataBaseUtils::DbSession> sql = m_specViewer->sql();
  DataBaseUtils::DbTransaction transaction( *sql );
  
  try
  {
//    m_modelInDb.reread();
    m_modelInDb.modify()->serializeTime = WDateTime::currentDateTime();
    rapidxml::xml_document<char> doc;
    serialize( &doc );
    m_modelInDb.modify()->xmlData.clear();
    rapidxml::print(std::back_inserter(m_modelInDb.modify()->xmlData), doc, 0);
    transaction.commit();
    cerr << "\n\nFinishing in ShieldingSourceDisplay::saveModelIfAlreadyInDatabase()" << endl;
  }catch( std::exception &e )
  {
    cerr << "\n\nException caught in ShieldingSourceDisplay::saveModelIfAlreadyInDatabase(): "
    << e.what() << endl;
    transaction.rollback();
  }//try / catch
}//void saveModelIfAlreadyInDatabase()
#endif //#if( USE_DB_TO_STORE_SPECTRA )


void ShieldingSourceDisplay::testSerialization()
{
  rapidxml::xml_document<char> doc;
  
  cerr << "Staring serializtion" << endl;
  serialize( &doc );
  std::string s;
  rapidxml::print(std::back_inserter(s), doc, 0);
  cerr << "Serialized:\n" << doc << endl;
  cerr << "\n\nStarting de-serializtion" << endl;
  deSerialize( &doc );
  cerr << "Done deserializing\n\n" << endl;
}//testSerialization()


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
  value = doc->allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );

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
    
    name = "Thickness";
    value = doc->allocate_string( m_thicknessEdit->valueText().toUTF8().c_str() );
    node = doc->allocate_node( rapidxml::node_element, name, value );
    material_node->append_node( node );
    if( m_forFitting )
    {
      value = m_fitThicknessCB->isChecked() ? "1" : "0";
      attr = doc->allocate_attribute( "Fit", value );
      node->append_attribute( attr );
    }//if( m_forFitting )
    
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
  rapidxml::xml_node<char> *node, *generic_node, *material_node;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
 
  if( shield_node->name() != string("Shielding") )
    throw runtime_error( "ShieldingSelects XML node should be 'Shielding'" );
  
  attr = shield_node->first_attribute( "version", 7 );
  int version;
  if( !attr || !attr->value() || !(stringstream(attr->value())>>version) )
    throw runtime_error( "ShieldingSelects should be versioned" );
  
  if( version != sm_xmlSerializationVersion )
    throw runtime_error( "Invalid XML version for ShieldingSelect" );
  
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
    
    
    
    if( !m_isGenericMaterial )
      handleToggleGeneric();
    
    cout << "AD=" << ad_node->value() << endl;
    cout << "AN=" << an_node->value() << endl;
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
    const rapidxml::xml_node<> *frac_node, *iso_node, *name_node, *thick_node;

    name_node = material_node->first_node( "Name", 4 );
    thick_node = material_node->first_node( "Thickness", 9 );
    if( !name_node || !name_node->value()
        || !thick_node || !thick_node->value() )
      throw runtime_error( "Material node didnt have name or thickness child" );
    
    if( m_isGenericMaterial )
      handleToggleGeneric();
    
    m_materialEdit->setValueText( WString::fromUTF8(name_node->value()) );
    m_thicknessEdit->setValueText( WString::fromUTF8(thick_node->value()) );
    
    handleMaterialChange();
    
    if( m_forFitting )
    {
      bool fit;
      attr = thick_node->first_attribute( "Fit", 3 );
      if( !attr || !attr->value() || !(stringstream(attr->value())>>fit) )
        throw runtime_error( "Material node expected thickness Fit attribute" );
      m_fitThicknessCB->setChecked( fit );
      
      
      const rapidxml::xml_node<> *fitmassfrac_node = material_node->first_node( "FitMassFraction", 15 );
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
        throw runtime_error( string(name_node->value()) + " is not a "
                             "valid isotope" );
      
      srcnuclides.push_back( nuc );
      
      double fraction;
      if( !(stringstream(frac_node->value()) >> fraction) )
        throw runtime_error( "Invalid mass fraction: "
                             + string(frac_node->value()) );
      
//      addSourceIsotopeCheckBox( nuc );
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



void ShieldingSourceDisplay::deSerializePeaksToUse(
                                    const rapidxml::xml_node<char> *peaks_node )
{
  const rapidxml::xml_attribute<char> *attr;
  const rapidxml::xml_node<char> *symbol_node, *energy_node, *peak_node;
  
  for( peak_node = peaks_node->first_node( "Peak", 4 );
       peak_node; peak_node = peak_node->next_sibling( "Peak", 4 ) )
  {
    attr = peak_node->first_attribute( "Use", 3 );
    symbol_node = peak_node->first_node( "Nuclide", 7 );
    energy_node = peak_node->first_node( "Energy", 6 );
    
    bool use = true;
    double energy;
    if( attr && attr->value() && !(stringstream(attr->value())>>use) )
      throw runtime_error( "Invalid 'Use' attribute in Peak XML element" );
    if( !use )
      continue;
    if( !symbol_node || !energy_node
       || !symbol_node->value() || !energy_node->value()
       || !(stringstream(energy_node->value()) >> energy) )
      throw runtime_error( "Invalid or missing node for Peak XML element" );
    
    //    PeakModel::PeakShrdPtr nearestPeak = m_peakModel->nearestPeak( energy );
    
    double nearestE = DBL_MAX;
    WModelIndex nearest_index;
    PeakModel::PeakShrdPtr nearest_peak;
    for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    {
      WModelIndex index = m_peakModel->index( peakn, PeakModel::kUseForShieldingSourceFit );
      PeakModel::PeakShrdPtr peak = m_peakModel->peak( index );
      const SandiaDecay::Nuclide *nuclide = peak->parentNuclide();
      
      try
      {
        const float gamenergy = peak->gammaParticleEnergy();
        const double dE = fabs(gamenergy - energy);
        if( nuclide && (dE < nearestE) )
        {
          nearestE = dE;
          nearest_peak = peak;
          nearest_index = index;
        }//if( this is nearest candidate peak )
      }catch( std::exception & )
      {
        
      }//try / catch
    }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
    
    if( nearest_peak && (nearestE < 0.1) )
      m_peakModel->setData( nearest_index, true );
  }//for( node = peaks_node->first_node(); node; node = node->next_sibling() )
}//void deSerializePeaksToUse( rapidxml::xml_node<char> *peaks_node )


void ShieldingSourceDisplay::deSerializeSourcesToFitFor( const rapidxml::xml_node<char> *sources )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  vector<SourceFitModel::IsoFitStruct> &model_isos = m_sourceModel->m_nuclides;
  vector<bool> model_row_modded( model_isos.size(), false );
  
  for( const rapidxml::xml_node<> *src_node = sources->first_node( "Nuclide", 7 );
       src_node; src_node = src_node->next_sibling( "Nuclide", 7 ) )
  {
    const rapidxml::xml_node<> *name_node = src_node->first_node( "Name", 4 );
    const rapidxml::xml_node<> *activity_node = src_node->first_node( "Activity", 8 );
    const rapidxml::xml_node<> *determined_note = src_node->first_node( "ShieldingDeterminedActivity", 27 );
    
    if( !name_node || !name_node->value() || !activity_node
        || !determined_note || !determined_note->value() )
      throw runtime_error( "Missing necessary element for sources XML" );
    
    const rapidxml::xml_node<> *activity_value_node = activity_node->first_node( "Value", 5 );
    const rapidxml::xml_node<> *activity_uncert_node = activity_node->first_node( "Uncertainty", 11 );
    
    if( !activity_value_node )
      throw runtime_error( "No activity value node" );
    const rapidxml::xml_attribute<char> *fit_activity_attr = activity_value_node->first_attribute( "Fit", 3 );
    
    const rapidxml::xml_node<> *age_node = src_node->first_node( "Age", 3 );
    if( !age_node )
      throw runtime_error( "Missing necessary age element for sources XML" );
    
    const rapidxml::xml_node<> *age_value_node = age_node->first_node( "Value", 5 );
    const rapidxml::xml_attribute<char> *fit_age_attr = age_value_node->first_attribute( "Fit", 3 );
    const rapidxml::xml_attribute<char> *age_defining_attr = age_value_node->first_attribute( "AgeDefiningNuclide", 18 );
    if( !age_defining_attr )
      age_defining_attr = age_value_node->first_attribute( "AgeMaster", 9 ); //sm_xmlSerializationVersion
    
    const rapidxml::xml_node<> *age_uncert_node = age_node->first_node( "Uncertainty", 11 );
    
    if( !activity_value_node || !activity_value_node->value()
        || !activity_uncert_node || !activity_uncert_node->value()
        || !age_value_node || !age_value_node->value()
        || !fit_activity_attr || !fit_activity_attr->value()
        || !age_value_node || !age_value_node->value()
        || !fit_age_attr || !fit_age_attr->value()
        || !age_uncert_node || !age_uncert_node->value() )
      throw runtime_error( "Missing/invalid node for sources XML" );
    
    SourceFitModel::IsoFitStruct row;
    row.nuclide = db->nuclide( name_node->value() );
    if( !row.nuclide )
      throw runtime_error( "Invalid nuclide for sources XML" );
    if( !(stringstream(activity_value_node->value()) >> row.activity) )
      throw runtime_error( "Failed to read activity" );
    if( !(stringstream(activity_uncert_node->value()) >> row.activityUncertainty) )
      throw runtime_error( "Failed to read activity_uncer" );
    if( !(stringstream(age_value_node->value()) >> row.age) )
      throw runtime_error( "Failed to read age" );
    if( !(stringstream(age_uncert_node->value()) >> row.ageUncertainty) )
      throw runtime_error( "Failed to read age_uncert" );
    if( !(stringstream(fit_activity_attr->value()) >> row.fitActivity) )
      throw runtime_error( "Failed to read fit_act" );
    if( !(stringstream(fit_age_attr->value()) >> row.fitAge) )
      throw runtime_error( "Failed to read fit_age" );
    if( !(stringstream(determined_note->value()) >> row.shieldingIsSource) )
      throw runtime_error( "Failed to read shieldingIsSource" );

    row.ageIsFittable = !PeakDef::ageFitNotAllowed( row.nuclide );
    
    if( !age_defining_attr || !age_defining_attr->value() )
      row.ageDefiningNuc = nullptr;
    else
      row.ageDefiningNuc = db->nuclide( age_defining_attr->value() );
    
    row.activity /= GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
    row.activityUncertainty /= GammaInteractionCalc::PointSourceShieldingChi2Fcn::sm_activityUnits;
    
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    auto getTruth = [src_node]( const char *truthName, const bool isActivity,
                                boost::optional<double> &value ){
      value.reset();
      
      auto node = src_node->first_node( truthName );
      if( !node )
        return;
      
      try
      {
        if( isActivity )
          value = PhysicalUnits::stringToActivity( node->value() );
        else
          value = PhysicalUnits::stringToTimeDuration( node->value() );
        
        cout << "Set '" << truthName << "' to value " << *value << " from '" << node->value() <<  "' while deserializing" << endl;
      }catch(...)
      {
        cerr << "Failed to read back in " << truthName << " from " << node->value() << endl;
      }
    };//getTruth(...)
    
    getTruth( "TruthActivity", true, row.truthActivity );
    getTruth( "TruthActivityTolerance", true, row.truthActivityTolerance );
    getTruth( "TruthAge", false, row.truthAge );
    getTruth( "TruthAgeTolerance", false, row.truthAgeTolerance );
#endif
    
    for( size_t i = 0; i < model_isos.size(); ++i )
    {
      SourceFitModel::IsoFitStruct &cand = model_isos[i];
      if( model_row_modded[i] || (cand.nuclide != row.nuclide) )
        continue;
      row.numProgenyPeaksSelected = cand.numProgenyPeaksSelected;
      cand = row;
      model_row_modded[i] = true;
    }//for( size_t i = 0; i < m_sourceModel->m_nuclides.size(); ++i )
  }//for( loop over source isotopes )
  
  
  //We should have modded all rows that are currently in m_sourceModel, but
  //  we're not garunteed all rows that where in the XML are also in m_sourceModel
  //Note 20200929: this check below is a little over-optimistic, so I commented it out; if you
  //  stumble on this code in like a year, and its still commentd out, just delete it, because:
  //  If m_sourceModel had "Tl201" because that is the peaks in our current spectrum (and that is
  //  the peak that is marked to be used in a fit), but the XML model we're loading had "I131" in
  //  it (because thats what we fit for in the previous record of the spectrum file), then
  //  model_row_modded for this model will still be false.
  //for( const bool modded : model_row_modded )
  //{
  //  if( !modded )
  //    throw runtime_error( "Inconsistent state of source model and serialized from XML" );
  //}//for( const bool modded : model_row_modded )
  
  if( m_sourceModel->rowCount() )
  {
    WModelIndex topLeft = m_sourceModel->index( 0, 0 );
    WModelIndex bottomRight = m_sourceModel->index( m_sourceModel->rowCount()-1,
                                              m_sourceModel->columnCount()- 1 );
    m_sourceModel->dataChanged().emit( topLeft, bottomRight );
  }//if( m_sourceModel->rowCount() )
}//void deSerializeSourcesToFitFor( rapidxml::xml_node<char> *sources );


void ShieldingSourceDisplay::deSerializeShieldings( const rapidxml::xml_node<char> *shiledings )
{
  const rapidxml::xml_node<> *shield_node;
  
  for( shield_node = shiledings->first_node( "Shielding", 9 );
      shield_node; shield_node = shield_node->next_sibling( "Shielding", 9 ) )
  {
    ShieldingSelect *select = addShielding( nullptr, false );
    select->deSerialize( shield_node );
  }//for( loop over shieldings )
  
}//void deSerializeShieldings( const rapidxml::xml_node<char> *shiledings )


void ShieldingSourceDisplay::reset()
{
  m_modifiedThisForeground = false;
  m_prevDistStr = "100 cm";
  m_distanceEdit->setValueText( m_prevDistStr );
  const vector<WWidget *> shieldings = m_shieldingSelects->children();
  for( WWidget *child : shieldings )
  {
    const ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( child );
    if( select )
      delete select;
  }//for( WWebWidget *child : shieldings )
  
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    WModelIndex index = m_peakModel->index( peakn,
                                         PeakModel::kUseForShieldingSourceFit );
    const PeakModel::PeakShrdPtr &peak = m_peakModel->peak( index );
    if( peak->useForShieldingSourceFit() )
      m_peakModel->setData( index, false );
  }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  
  //We shouldnt actually need the next line
//  m_sourceModel->repopulateIsotopes();
}//void reset()


void ShieldingSourceDisplay::newForegroundSet()
{
  m_modifiedThisForeground = false;
}//void newForegroundSet()


bool ShieldingSourceDisplay::userChangedDuringCurrentForeground() const
{
  return m_modifiedThisForeground;
}//bool userChangedDuringCurrentForeground() const;


void ShieldingSourceDisplay::deSerialize( const rapidxml::xml_node<char> *base_node )
{
  const rapidxml::xml_node<char> *muti_iso_node, *back_sub_node,
                                 *peaks_node, *isotope_nodes, *shieldings_node,
                                 *dist_node, *same_age_node, *chart_disp_node;
  const rapidxml::xml_attribute<char> *attr;
  
  if( !base_node )
    throw runtime_error( "No ShieldingSourceFit node" );
  
  if( !rapidxml::internal::compare( base_node->name(), base_node->name_size(), "ShieldingSourceFit", 18, true) )
    throw runtime_error( "ShieldingSourceDisplay::deSerialize: invalid node name passed in: '"
                        + std::string(base_node->name(),base_node->name()+base_node->name_size()) + "'" );
  
  muti_iso_node   = base_node->first_node( "MultipleIsotopesPerPeak", 23 );
  back_sub_node   = base_node->first_node( "BackgroundPeakSubtraction", 25 );
  same_age_node   = base_node->first_node( "SameAgeIsotopes", 15 );
  chart_disp_node = base_node->first_node( "ShowChiOnChart", 14 );
  peaks_node      = base_node->first_node( "Peaks", 5 );
  isotope_nodes   = base_node->first_node( "Nuclides", 8 );
  shieldings_node = base_node->first_node( "Shieldings", 10 );
  dist_node       = base_node->first_node( "Distance", 8 );
  
  if( !peaks_node || !isotope_nodes || !shieldings_node )
    throw runtime_error( "Missing necessary XML node" );
  
  int version;
  bool muti_iso, back_sub, same_age = false, show_chi_on_chart = true;
  
  attr = base_node->first_attribute( "version", 7 );
  if( !attr || !attr->value() || !(stringstream(attr->value())>>version) )
    throw runtime_error( "Deserializing requires a version" );
  
  if( version != sm_xmlSerializationVersion )
    throw runtime_error( "Invalid version of ShieldingSourceDisplay XML" );
  
  if( !muti_iso_node || !muti_iso_node->value()
      || !(stringstream(muti_iso_node->value()) >> muti_iso) )
    throw runtime_error( "Invalid or missing MultipleIsotopesPerPeak node" );
  
  
  if( same_age_node && same_age_node->value() ) //not a mandatory element
  {
    if( !(stringstream(muti_iso_node->value()) >> same_age) )
      throw runtime_error( "Invalid SameAgeIsotopes node" );
  }//if( same_age_node && same_age_node->value() )

  if( chart_disp_node && chart_disp_node->value() )
  {
    if( !(stringstream(muti_iso_node->value()) >> show_chi_on_chart) )
      throw runtime_error( "Invalid ShowChiOnChart node" );
  }
  
  if( !back_sub_node || !back_sub_node->value()
     || !(stringstream(back_sub_node->value()) >> back_sub) )
    throw runtime_error( "Invalid or missing BackgroundPeakSubtraction node" );
  back_sub = (back_sub && m_specViewer->measurment(SpecUtils::SpectrumType::Background));
  
  if( !dist_node || !dist_node->value() )
    throw runtime_error( "Invalid or missing Distance node" );
  
  //clear out the GUI
  reset();
  
  m_multiIsoPerPeak->setChecked( muti_iso );
  m_backgroundPeakSub->setChecked( back_sub );
  m_sameIsotopesAge->setChecked( same_age );
  m_showChiOnChart->setChecked( show_chi_on_chart );
  m_chi2Graphic->setShowChiOnChart( show_chi_on_chart );
  m_distanceEdit->setValueText( WString::fromUTF8(dist_node->value()) );
  m_prevDistStr = dist_node->value();
  
  deSerializePeaksToUse( peaks_node );
  m_sourceModel->repopulateIsotopes();
  deSerializeSourcesToFitFor( isotope_nodes );
  deSerializeShieldings( shieldings_node );
  
  m_modifiedThisForeground = true;
  
  updateChi2Chart();
}//void deSerialize( rapidxml::xml_document<char> &doc )


::rapidxml::xml_node<char> *ShieldingSourceDisplay::serialize( rapidxml::xml_node<char> *parent_node )
{
  rapidxml::xml_document<char> *doc = parent_node->document();
  
  const char *name, *value;
  rapidxml::xml_node<> *base_node, *dist_node, *peaks_node, *isotope_nodes,
                       *shieldings_node;
  rapidxml::xml_attribute<> *attr;
  
  name = "ShieldingSourceFit";
  base_node = doc->allocate_node( rapidxml::node_element, name );
  parent_node->append_node( base_node );
  
  //If you change the available options or formatting or whatever, increment the
  //  version field of the XML!
  value = doc->allocate_string( std::to_string(sm_xmlSerializationVersion).c_str() );
  attr = doc->allocate_attribute( "version", value );
  base_node->append_attribute( attr );
  
  name = "MultipleIsotopesPerPeak";
  value = m_multiIsoPerPeak->isChecked() ? "1" : "0";
  dist_node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( dist_node );

  name = "BackgroundPeakSubtraction";
  value = m_backgroundPeakSub->isChecked() ? "1" : "0";
  dist_node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( dist_node );
  
  name = "SameAgeIsotopes";
  value = m_sameIsotopesAge->isChecked() ? "1" : "0";
  dist_node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( dist_node );

  name = "ShowChiOnChart";
  value = m_showChiOnChart->isChecked() ? "1" : "0";
  dist_node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( dist_node );
  
  name = "Distance";
  value = doc->allocate_string( m_distanceEdit->valueText().toUTF8().c_str() );
  dist_node = doc->allocate_node( rapidxml::node_element, name, value );
  base_node->append_node( dist_node );

  name = "Shieldings";
  shieldings_node = doc->allocate_node( rapidxml::node_element, name );
  base_node->append_node( shieldings_node );

  
  const vector<WWidget *> shieldings = m_shieldingSelects->children();
  for( WWidget *child : shieldings )
  {
    const ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( child );
    if( select )
      select->serialize( shieldings_node );
  }//for( WWebWidget *child : shieldings )
  
  peaks_node = doc->allocate_node( rapidxml::node_element, "Peaks" );
  base_node->append_node( peaks_node );
  
  for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  {
    const PeakDef &peak = m_peakModel->peak( peakn );
    const SandiaDecay::Nuclide *nuclide = peak.parentNuclide();
    const SandiaDecay::RadParticle *particle = peak.decayParticle();
    const bool annhilation = peak.sourceGammaType()==PeakDef::AnnihilationGamma;
    
    if( peak.useForShieldingSourceFit() && nuclide && (particle || annhilation) )
    {
      rapidxml::xml_node<> *peak_node, *nuc_node, *energy_node;
      peak_node = doc->allocate_node( rapidxml::node_element, "Peak" );
      peaks_node->append_node( peak_node );
      
      //meh, lets be explicit about using this peak, although were only writing
      //  peaks were using for the fit
      attr = doc->allocate_attribute( "Use", "1" );
      peak_node->append_attribute( attr );
      
      value = doc->allocate_string( nuclide->symbol.c_str() );
      nuc_node = doc->allocate_node( rapidxml::node_element, "Nuclide", value );
      peak_node->append_node( nuc_node );
      
      const float energy = peak.gammaParticleEnergy();
      value = doc->allocate_string( std::to_string(energy).c_str() );
      energy_node = doc->allocate_node( rapidxml::node_element, "Energy", value );
      peak_node->append_node( energy_node );
    }//if( peak.useForShieldingSourceFit() )
  }//for( int peakn = 0; peakn < m_peakModel->rowCount(); ++peakn )
  
  isotope_nodes = doc->allocate_node( rapidxml::node_element, "Nuclides" );
  base_node->append_node( isotope_nodes );
  
  for( int nuc = 0; nuc < m_sourceModel->rowCount(); ++nuc )
  {
    rapidxml::xml_node<> *nuclide_node, *activity_node, *node;
    rapidxml::xml_node<> *age_node, *determined_note, *name_node;
    
    const SandiaDecay::Nuclide *nuclide = m_sourceModel->nuclide( nuc );
    const double activity = m_sourceModel->activity( nuc );
    const double activityUncert = m_sourceModel->activityUncert( nuc );
    const bool fitActivity = m_sourceModel->fitActivity( nuc );
    const double age = m_sourceModel->age( nuc );
    const double ageUncert = m_sourceModel->ageUncert( nuc );
    const bool fitAge = m_sourceModel->fitAge( nuc );
    const SandiaDecay::Nuclide *ageNuc = m_sourceModel->ageDefiningNuclide( nuclide );
    
    const bool shieldingDeterminedActivity
                           = m_sourceModel->shieldingDeterminedActivity( nuc );
    
    nuclide_node = doc->allocate_node( rapidxml::node_element, "Nuclide" );
    isotope_nodes->append_node( nuclide_node );

    value = doc->allocate_string( nuclide->symbol.c_str() );
    name_node = doc->allocate_node( rapidxml::node_element, "Name", value );
    nuclide_node->append_node( name_node );
    
    value = shieldingDeterminedActivity ? "1" : "0";
    determined_note = doc->allocate_node( rapidxml::node_element, "ShieldingDeterminedActivity", value );
    nuclide_node->append_node( determined_note );
    
    activity_node = doc->allocate_node( rapidxml::node_element, "Activity" );
    nuclide_node->append_node( activity_node );
    
    value = doc->allocate_string( std::to_string(activity).c_str() );
    node = doc->allocate_node( rapidxml::node_element, "Value", value );
    activity_node->append_node( node );
    
    value = fitActivity ? "1" : "0";
    attr = doc->allocate_attribute( "Fit", value );
    node->append_attribute( attr );
    
    value = doc->allocate_string( std::to_string(activityUncert).c_str() );
    node = doc->allocate_node( rapidxml::node_element, "Uncertainty", value );
    activity_node->append_node( node );
    
    age_node = doc->allocate_node( rapidxml::node_element, "Age" );
    nuclide_node->append_node( age_node );
    
    value = doc->allocate_string( std::to_string(age).c_str() );
    node = doc->allocate_node( rapidxml::node_element, "Value", value );
    age_node->append_node( node );
    
    value = fitAge ? "1" : "0";
    attr = doc->allocate_attribute( "Fit", value );
    node->append_attribute( attr );
    
    if( ageNuc && (ageNuc != nuclide) )
    {
      value = ageNuc->symbol.c_str();
      
      if( ShieldingSourceDisplay::sm_xmlSerializationVersion == 0 )
      {
        //Depreciating tag 20201201, will remove when XML serialization is updated
        attr = doc->allocate_attribute( "AgeMaster", value );
        node->append_attribute( attr );
      }
      
      attr = doc->allocate_attribute( "AgeDefiningNuclide", value );
      node->append_attribute( attr );
    }//if( ageNuc && (ageNuc != nuclide) )
    
    value = doc->allocate_string( std::to_string(ageUncert).c_str() );
    node = doc->allocate_node( rapidxml::node_element, "Uncertainty", value );
    age_node->append_node( node );
    
#if( INCLUDE_ANALYSIS_TEST_SUITE )
    auto addTruth = [doc,nuclide_node]( const char *truthName, const bool isActivity,
                                        const boost::optional<double> &value ){
      if( !value )
        return;
      string strval;
      const bool useCurries = true;
      if( isActivity )
        strval = PhysicalUnits::printToBestActivityUnits( *value, 6, useCurries );
      else
        strval = PhysicalUnits::printToBestTimeUnits( *value, 6 );
      
      const char *txtvalue = doc->allocate_string( strval.c_str() );
      rapidxml::xml_node<char> *node = doc->allocate_node( rapidxml::node_element, truthName, txtvalue );
      nuclide_node->append_node( node );
    };//addTruth(...)
    
    addTruth( "TruthActivity", true, m_sourceModel->truthActivity(nuc) );
    addTruth( "TruthActivityTolerance", true, m_sourceModel->truthActivityTolerance(nuc) );
    addTruth( "TruthAge", false, m_sourceModel->truthAge(nuc) );
    addTruth( "TruthAgeTolerance", false, m_sourceModel->truthAgeTolerance(nuc) );
#endif
  }//for( int nuc = 0; nuc < m_sourceModel->rowCount(); ++nuc )

  
  
  try
  {
    vector<ShieldingSelect *> shieldings;
    ROOT::Minuit2::MnUserParameters inputPrams;
    Chi2FcnShrdPtr chi2Fcn = shieldingFitnessFcn( shieldings, inputPrams );
    const unsigned int ndof = inputPrams.VariableParameters();
    const vector<double> params = inputPrams.Params();
    const vector<double> errors = inputPrams.Errors();
    GammaInteractionCalc::PointSourceShieldingChi2Fcn::NucMixtureCache mixcache;
    const bool multy_iso = m_multiIsoPerPeak->isChecked();
    const vector< tuple<double,double,double,WColor,double> > chis
     = chi2Fcn->energy_chi_contributions( params, mixcache, multy_iso, nullptr );
    
    if( chis.size() )
    {
      char buffer[64] = { 0 };
      
      //We will write this information - for the record, JIC, but we will not
      //  read it back in later, as we will re-calculate it from the actual
      //  data present when this model gets deserialized.
      rapidxml::xml_node<> *chi2_node = doc->allocate_node( rapidxml::node_element, "Chi2Elements" );
      base_node->append_node( chi2_node );
      
      double chi2 = 0.0;
      rapidxml::xml_node<> *node = nullptr;
      rapidxml::xml_node<> *eval_node = doc->allocate_node( rapidxml::node_element, "EvaluatedEnergies" );
      
      for( const auto &p : chis )
      {
        rapidxml::xml_node<> *point_node = doc->allocate_node( rapidxml::node_element, "EvalPoint" );
        eval_node->append_node( point_node );
        
        const double energy = get<0>(p);
        const double chi = get<1>(p);
        const double scale = get<2>(p);
        const WColor &color = get<3>(p);
        const double scale_uncert = get<4>(p);
        
        if( !IsInf(chi) && !IsNan(chi) )
          chi2 += chi*chi;
        
        node = doc->allocate_node( rapidxml::node_element, "Energy" );
        snprintf( buffer, sizeof(buffer), "%f", energy );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Chi" );
        snprintf( buffer, sizeof(buffer), "%f", chi );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Scale" );
        snprintf( buffer, sizeof(buffer), "%f", scale );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "Color" );
        snprintf( buffer, sizeof(buffer), "%s", (color.isDefault() ? "" : color.cssText().c_str()) );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
        
        node = doc->allocate_node( rapidxml::node_element, "ScaleUncert" );
        snprintf( buffer, sizeof(buffer), "%f", scale_uncert );
        value = doc->allocate_string( buffer );
        node->value( value );
        point_node->append_node( node );
      }//for( const auto &p : chis )
      
      node = doc->allocate_node( rapidxml::node_element, "Chi2" );
      snprintf( buffer, sizeof(buffer), "%f", chi2 );
      value = doc->allocate_string( buffer );
      node->value( value );
      chi2_node->append_node( node );
      
      node = doc->allocate_node( rapidxml::node_element, "NumParamFit" );
      snprintf( buffer, sizeof(buffer), "%u", ndof );
      value = doc->allocate_string( buffer );
      node->value( value );
      chi2_node->append_node( node );
      
      chi2_node->append_node( eval_node );
    }//if( chis.size() )
  }catch( std::exception &e )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, ("Failed to get chi2 info during serialization - caught exception: " + string(e.what())).c_str() );
#endif
  }
  
  return parent_node;
}//::rapidxml::xml_node<char> * serialize()


/*
void ShieldingSourceDisplay::guessDetectorType( SpecUtils::SpectrumType type,
                                                std::shared_ptr<SpecMeas> measurment,
                                                std::set<int> sample_numbers )
{
  if( type != SpecUtils::SpectrumType::Foreground )
    return;

  if( !measurment )
    return;

  const string &manufacturer = measurment->manufacturer();
  const string &instrument_id = measurment->instrument_id();
  const string &instrument_model = measurment->instrument_model();


  cerr << "manufacturer=" << manufacturer << ", instrument_id=" << instrument_id << ", instrument_model=" << instrument_model << endl;

  //manufacturer=ORTEC, instrument_id=Serial #10000033, instrument_model=OSASP


}//void ShieldingSourceDisplay::guessDetectorType()
*/


void ShieldingSourceDisplay::layoutSizeChanged( int width, int height )
{
    //If phone, we just force to show since we don't care about the screen size
    if (m_specViewer->isPhone())
      return;
    
  m_width = width;
  m_height = height;

  if( m_nResizeSinceHint > 3 )
  {
    if( m_height < 320 && !m_chi2Graphic->isHidden() )
    {
      m_chi2Graphic->hide();
//      m_layout->setRowStretch( 1, 0 );
      m_showChi2Text->show();
    }else if( m_chi2Graphic->isHidden() && m_height >= 320 )
    {
      m_chi2Graphic->show();
      m_showChi2Text->hide();
//      m_layout->setRowStretch( 1, 6 );
//      leftLayout->setRowResizable( 0, true );
//      leftLayout->setRowResizable( 1, true );
    }//if / else
  }//if( m_nResizeSinceHint > 3 )

  ++m_nResizeSinceHint;
}//void layoutSizeChanged( int width, int height )


void ShieldingSourceDisplay::initialSizeHint( int width, int height )
{
  m_nResizeSinceHint = 5;
  layoutSizeChanged( width, height );
  m_nResizeSinceHint = 0;
}//void initialSizeHint( int width, int height )


void ShieldingSourceDisplay::doAddShielding()
{
  addShielding( nullptr, true );
}

void ShieldingSourceDisplay::doAddShieldingBefore( ShieldingSelect *select )
{
  addShielding( select, true );
}


void ShieldingSourceDisplay::doAddShieldingAfter( ShieldingSelect *select )
{
  const vector<WWidget *> selects = m_shieldingSelects->children();
  
  int index = -1;
  for( size_t i = 0; i < selects.size(); ++i )
  {
    ShieldingSelect *s = dynamic_cast<ShieldingSelect *>( selects[i] );
      
    if( s == select )
      index = static_cast<int>( i );
    
    if( s && index>=0 && (i < (selects.size()-1)) )
    {
      ShieldingSelect *next = dynamic_cast<ShieldingSelect *>( selects[i+1] );
      addShielding( next, true );
      return;
    }
  }//for( size_t i = 0; i < selects.size(); ++i )
  
  if( index < 0 )
    cerr << "ShieldingSourceDisplay::doAddShieldingAfter(): Warning, failed to"
    << " find passed in select" << endl;
  
  //if index is >= 0, then user requested to add a shielding after the last
  //  shielding
  addShielding( nullptr, true );
}//void doAddShieldingAfter( ShieldingSelect *select )


void ShieldingSourceDisplay::addGenericShielding()
{
  ShieldingSelect *temp = addShielding( nullptr, true );
  if( !temp->isGenericMaterial() )
    temp->handleToggleGeneric();
}//void addGenericShielding()



ShieldingSelect *ShieldingSourceDisplay::addShielding( ShieldingSelect *before, const bool doUpdateChiChart )
{
  m_modifiedThisForeground = true;
  
  ShieldingSelect *select = NULL;

  select = new ShieldingSelect( m_materialDB, m_sourceModel,
                                m_materialSuggest, true );

  if( before && m_shieldingSelects->indexOf(before) >= 0 )
    m_shieldingSelects->insertBefore( select, before );
  else
    m_shieldingSelects->addWidget( select );
  
  
  select->addShieldingBefore().connect( boost::bind( &ShieldingSourceDisplay::doAddShieldingBefore, this, _1 ) );
  select->addShieldingAfter().connect( boost::bind( &ShieldingSourceDisplay::doAddShieldingAfter, this, _1 ) );
  
  //connect up signals of select and such
  select->m_arealDensityEdit->changed().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  select->m_atomicNumberEdit->changed().connect( this, &ShieldingSourceDisplay::updateChi2Chart );

  select->materialChanged().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  select->materialModified().connect( this, &ShieldingSourceDisplay::updateChi2Chart );
  
//  select->m_thicknessEdit->changed().connect( this, &ShieldingSourceDisplay::checkDistanceAndThicknessConsistent );
//  select->m_thicknessEdit->enterPressed().connect( this, &ShieldingSourceDisplay::checkDistanceAndThicknessConsistent );
//  select->m_thicknessEdit->blurred().connect( this, &ShieldingSourceDisplay::checkDistanceAndThicknessConsistent );
  
  select->remove().connect( this, &ShieldingSourceDisplay::removeShielding );
  select->materialModified().connect( this, &ShieldingSourceDisplay::materialModifiedCallback );
  select->materialChanged().connect( this, &ShieldingSourceDisplay::materialChangedCallback );

  select->addingIsotopeAsSource().connect( boost::bind( &ShieldingSourceDisplay::isotopeSelectedAsShieldingCallback, this, select, _1 ) );
  select->removingIsotopeAsSource().connect( boost::bind( &ShieldingSourceDisplay::isotopeDeSelectedAsShieldingCallback, this, select, _1 ) );
  select->activityFromThicknessNeedUpdating().connect( boost::bind( &ShieldingSourceDisplay::updateActivityOfShieldingIsotope, this, _1, _2 ) );

  if( doUpdateChiChart )
    updateChi2Chart();
  
  return select;
}//void addShielding()



void ShieldingSourceDisplay::addSourceIsotopesToShieldings( Wt::WModelIndex, int firstRow, int lastRow )
{
//  m_modifiedThisForeground = true;
  
  firstRow = max( 0, firstRow );
  lastRow = min( lastRow, m_sourceModel->rowCount() );

  vector<const SandiaDecay::Nuclide *> isotopes;
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      isotopes.push_back( nuc );
  }//for( int row = firstRow; row <= lastRow; ++row )

  if( isotopes.empty() )
    return;

  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;
    
    for( const SandiaDecay::Nuclide *iso : isotopes )
      thisSelect->addSourceIsotopeCheckBox( iso );
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void addSourceIsotopesToShieldings( Wt::WModelIndex, int firstRow, int lastRow )


void ShieldingSourceDisplay::removeSourceIsotopesFromShieldings( Wt::WModelIndex, int firstRow, int lastRow )
{
//  m_modifiedThisForeground = true;
  
  const int nrow = m_sourceModel->rowCount();
  if( firstRow<0 || firstRow>lastRow || lastRow<0 || lastRow>=nrow )
  {
    stringstream msg;
    msg << "ShieldingSourceDisplay::removeSourceIsotopesFromShieldings(...)\n\tfunction called with invalid argument "
        << "- this is a major logic error I think";
    throw std::runtime_error( msg.str() );
  }//if( invalid row )

  vector<const SandiaDecay::Nuclide *> isotopes;
  for( int row = firstRow; row <= lastRow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      isotopes.push_back( nuc );
  }//for( int row = firstRow; row <= lastRow; ++row )


  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;

    for( const SandiaDecay::Nuclide *iso : isotopes )
      thisSelect->removeSourceIsotopeCheckBox( iso );
  }//for( WWidget *widget : m_shieldingSelects->children() )
}//void removeSourceIsotopesFromShieldings( Wt::WModelIndex, int firstRow, int lastRow )


void ShieldingSourceDisplay::materialModifiedCallback( ShieldingSelect *select )
{
  if( !select )
  {
    cerr << "ShieldingSourceDisplay::materialModifiedCallback(...)\n\tShouldnt be here!" << endl;
    return;
  }//if( !select )

  m_modifiedThisForeground = true;
  
  cerr << "In ShieldingSourceDisplay::materialModifiedCallback(...)" << endl;
  //I meant to do some more work here...

}//void materialModifiedCallback( ShieldingSelect *select )



void ShieldingSourceDisplay::materialChangedCallback( ShieldingSelect *select )
{
  if( !select )
  {
    cerr << "ShieldingSourceDisplay::materialChangedCallback(...)\n\tShouldnt be here!" << endl;
    return;
  }//if( !select )

  //The select has alread removed any isotopes as sources which arent in the
  //  current material, however, we have to add in all other candidate isotopes

  m_modifiedThisForeground = true;
  
  const int nrow = m_sourceModel->rowCount();
  for( int row = 0; row < nrow; ++row )
  {
    const SandiaDecay::Nuclide *nuc = m_sourceModel->nuclide( row );
    if( nuc )
      select->addSourceIsotopeCheckBox( nuc );
  }//for( int row = 0; row < nrow; ++row )

}//void materialChangedCallback( ShieldingSelect *select )



void ShieldingSourceDisplay::updateActivityOfShieldingIsotope( ShieldingSelect *select,
                                       const SandiaDecay::Nuclide *nuc )
{
  typedef pair<const SandiaDecay::Nuclide *,float> NuclideFrac;

  if( !select || select->isGenericMaterial() || !select->material() )
  {
    stringstream msg;
    msg << "ShieldingSourceDisplay::updateActivityOfShieldingIsotope(...)\n\tShould not be here!";
    throw std::runtime_error( msg.str() );
  }//if( !select || select->isGenericMaterial() )

  double inner_rad = 0.0;
  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *thisSelect = dynamic_cast<ShieldingSelect *>(widget);
    if( !thisSelect || thisSelect->isGenericMaterial() )
      continue;

    if( thisSelect == select )
      break;

    inner_rad += thisSelect->thickness();
  }//for( WWidget *widget : m_shieldingSelects->children() )

  const double outer_rad = inner_rad + select->thickness();

  const double pi = 3.14159265359;
  const double volume = pi* (4.0/3.0) * ( pow(outer_rad,3) - pow(inner_rad,3) );

  std::shared_ptr<const Material> material = select->material();
  const double density = material->density;


  double weight = density * volume;
  const SandiaDecay::Nuclide *nuclide = NULL;

  for( const NuclideFrac &ef : material->nuclides )
  {
    if( ef.first && (ef.first==nuc) )
    {
      nuclide = ef.first;
      weight *= ef.second;
    }
  }//for( const NuclideFrac &ef : material->nuclides )

  /*
   *XXXHere
   */

  if( !nuclide )
  {
    stringstream msg;
    msg << "ShieldingSourceDisplay::updateActivityOfShieldingIsotope(...)\n\tShould not be here! "
        << "nuc=" << nuc << " and material=" << material->name;
    throw std::runtime_error( msg.str() );
  }//if( !nuclide )

  const double weight_in_grams = weight / PhysicalUnits::gram;
  const double activity_per_gram = nuclide->activityPerGram();
  const double activity = activity_per_gram * weight_in_grams / PhysicalUnits::curie;

//  cerr << "updateActivityOfShieldingIsotope: weight_in_grams=" << weight_in_grams
//  << ", activity_per_gram=" << activity_per_gram
//  << ", activity=" << activity << " ci" << endl;
  
  const int row = m_sourceModel->row( nuclide );

  if( row < 0 )
  {
    stringstream msg;
    msg << "ShieldingSourceDisplay::updateActivityOfShieldingIsotope(...)\n\tShould not be here!";
    throw std::runtime_error( msg.str() );
  }//if( row < 0 )

  
  if( IsNan(activity) || IsInf(activity) )
  {
    stringstream msg;
    msg << "An invalid activity was calculated for " << nuclide->symbol
        << ", other results may be invalid";
    passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
    return;
  }//if( IsNan(p.second) || IsInf(p.second) )

  
  stringstream activitystrm;
  activitystrm << activity << "ci";

  WModelIndex index = m_sourceModel->index( row, SourceFitModel::kActivity );
  m_sourceModel->setData( index, activitystrm.str() );
}//void updateActivityOfShieldingIsotope(..)


void ShieldingSourceDisplay::isotopeSelectedAsShieldingCallback(
                            ShieldingSelect *caller, const SandiaDecay::Nuclide *nuc )
{
  //Make sure no other selects have this isotope selected
  const vector<WWidget *> &children = m_shieldingSelects->children();
  for( WWidget *widget : children )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( !select || (select==caller) )
      continue;

    select->uncheckSourceIsotopeCheckBox( nuc );
  }//for( WWidget *widget : children )

  //Set appropriate flags in the SourceFitModel so activity wont be editiable
  m_sourceModel->makeActivityNonEditable( nuc );

  //Change the age to something more reasonable for uranium, plutonium etc,
  //  XXX - this age change should definetly be improved, on a nuclide by
  //        nuclide basis
  string agestr = "";
  if( nuc )  //nuc should always be non-NULL, but just incase
  {
    if( nuc->atomicMass > 200.0 &&  nuc->halfLife>(20.0*PhysicalUnits::year) )
      agestr = "20y";
  }//if( nuc )

  if( agestr.size() )
  {
    WModelIndex index = m_sourceModel->index( nuc, SourceFitModel::kAge );
    m_sourceModel->setData( index, agestr );
  }//if( agestr.size() )

  WModelIndex fitAgeIndex = m_sourceModel->index( nuc, SourceFitModel::kFitAge );
  m_sourceModel->setData( fitAgeIndex, false );

  //Update activity displayed
  updateActivityOfShieldingIsotope( caller, nuc );


  //Update the Chi2
  updateChi2Chart();
}//isotopeSelectedAsShieldingCallback(...)


void ShieldingSourceDisplay::isotopeDeSelectedAsShieldingCallback(
                            ShieldingSelect *select, const SandiaDecay::Nuclide *nuc )
{
  //Set appropriate flags in the SourceFitModel so activity will be editiable
  m_sourceModel->makeActivityEditable( nuc );

  //Update the Chi2
  updateChi2Chart();
}//isotopeDeSelectedAsShieldingCallback(...)


void ShieldingSourceDisplay::removeShielding( ShieldingSelect *select )
{
  if( !select )
    return;

  m_modifiedThisForeground = true;
  
  const int nwidget = m_shieldingSelects->count();
  for( int i = 0; i < nwidget; ++i )
  {
    WWidget *widget = m_shieldingSelects->widget( i );
    ShieldingSelect *shielding = dynamic_cast<ShieldingSelect *>( widget );
    if( shielding == select )
    {
      delete shielding;
      return;
    }//if( shielding == select )
  }//for( int i = 0; i < nwidget; ++i )

  cerr << "\n\nCouldnt finding select to detelet" << endl;
}//void removeShielding( ShieldingSelect *select )





ShieldingSourceDisplay::Chi2FcnShrdPtr ShieldingSourceDisplay::shieldingFitnessFcn(
                                       vector<ShieldingSelect *> &shieldings,
                                       ROOT::Minuit2::MnUserParameters &inputPrams )
{
  //Get the peaks we'll be using in the fit
  vector<PeakDef> peaks;
  std::shared_ptr< const std::deque<PeakShrdPtr> > all_peaks = m_peakModel->peaks();
  if( all_peaks )
  {
    for( const PeakShrdPtr &peak : *all_peaks )
    {
      if( peak && peak->useForShieldingSourceFit() )
        peaks.push_back( *peak );
    }//for(...)
  }//if( all_peaks )

  if( peaks.empty() )
    throw runtime_error( "There are not peaks selected for the fit" );
  

  using GammaInteractionCalc::PointSourceShieldingChi2Fcn;
  double liveTime = m_specViewer->liveTime(SpecUtils::SpectrumType::Foreground) * PhysicalUnits::second;

  if( liveTime <= 0.0 )
  {
    if( m_sourceModel->numNuclides() )
      passMessage( "There was no defined detector live time, so assuming"
                   " 300 seconds", "ShieldingSourceDisplay::doModelFit",
                   WarningWidget::WarningMsgHigh );
    liveTime = 300.0 * PhysicalUnits::second;
  }//if( liveTime <= 0.0 )

  //Get the distance to the object
  const string distanceStr = m_distanceEdit->text().toUTF8();
  const double distance = PhysicalUnits::stringToDistance( distanceStr );

  //Get the shieldings and materials
  shieldings.clear();
  vector<PointSourceShieldingChi2Fcn::MaterialAndSources> materials;

  for( WWidget *widget : m_shieldingSelects->children() )
  {
    ShieldingSelect *select = dynamic_cast<ShieldingSelect *>( widget );
    if( select )
    {
      shieldings.push_back( select );

      std::shared_ptr<const Material> matSPtr = select->material();
      const Material *mat = matSPtr.get();
      if( !mat && !select->isGenericMaterial() )
        mat = m_materialDB->material( "void" );

      PointSourceShieldingChi2Fcn::MaterialAndSources materialAndSrc;
      materialAndSrc.first = mat;
      materialAndSrc.second = select->nuclidesToUseAsSources();

      materials.push_back( materialAndSrc );
    }//if( select )
  }//for( WWidget *widget : m_shieldingSelects->children() )

  
  std::shared_ptr<const DetectorPeakResponse> detector;
  if( m_specViewer->measurment(SpecUtils::SpectrumType::Foreground) )
    detector = m_specViewer->measurment(SpecUtils::SpectrumType::Foreground)->detector();
  
  std::shared_ptr<GammaInteractionCalc::PointSourceShieldingChi2Fcn> answer
    = std::make_shared<GammaInteractionCalc::PointSourceShieldingChi2Fcn>(
                       distance, liveTime, peaks, detector, materials,
                       m_multiIsoPerPeak->isChecked() );

  //I think num_fit_params will end up same as inputPrams.VariableParameters()
  size_t num_fit_params = 0;

  //Setup the parameters from the sources
  const int niso = m_sourceModel->numNuclides();

  if( niso < 1 )
    throw runtime_error( "There are no isotopes being fit for" );
  
  //XXX - need to check that if we are fitting a shielding thickness we have an
  //  isotope that has more than one peak in the fit, and if it is also fitting
  //  age should have more than two peaks.  There are probably a few more things
  //  to check here...
  
  for( size_t nucn = 0; nucn < answer->numNuclides(); ++nucn )
  {
    const SandiaDecay::Nuclide *nuclide = answer->nuclide( int(nucn) );
    const int ison = m_sourceModel->nuclideIndex( nuclide );
    double activity = m_sourceModel->activity( ison )
                            / PointSourceShieldingChi2Fcn::sm_activityUnits;
    
    //We could do a lot better here by estimating the activity of the sources
    //  the first time they are fit for
    //XXX - should make better guesses for source activity the first time a fit
    //      is performed
    
//    cerr << "Initial activity is " << m_sourceModel->activity( ison )/PhysicalUnits::curie
//         << " which is a minuit value of " << activity << endl;
//    activity = 10.0*PhysicalUnits::curie*(1.0E-6) / PointSourceShieldingChi2Fcn::sm_activityUnits;
    
    const bool fitAct = m_sourceModel->fitActivity( ison )
                        && !m_sourceModel->shieldingDeterminedActivity( ison );
    const bool fitAge = m_sourceModel->fitAge( ison );
    
    const SandiaDecay::Nuclide *ageDefiningNuc = m_sourceModel->ageDefiningNuclide( nuclide );
    const bool hasOwnAge = (!ageDefiningNuc || (ageDefiningNuc == nuclide));
    
    //cout << "For nuc=" << nuclide->symbol << " age=" << PhysicalUnits::printToBestTimeUnits(age)
    //     << ", fitage=" << PhysicalUnits::printToBestTimeUnits(fitAge)
    //     << ", hasOwnAge=" << hasOwnAge
    //     << ", ageDefiningNuc=" << (ageDefiningNuc ? ageDefiningNuc->symbol : string("null"))
    //     << endl;
    
    num_fit_params += fitAct + (fitAge && hasOwnAge);

    if( fitAct )
    {
      //We cant have both a specified lower and upper limit on activity, or
      //  such a large range will make Minuit2 choke and give a completely
      //  in-accurate answer (returns not even the best chi2 it found), if only
      //  one fit parameter.
//      inputPrams.Add( nuclide->symbol + "Strength", activity, activityStep, 0.0,
//                     10000.0*PhysicalUnits::curie/PointSourceShieldingChi2Fcn::sm_activityUnits );
      const string name = nuclide->symbol + "Strength";
      const double activityStep = (activity < 0.0001 ? 0.0001 : 0.1*activity);
      inputPrams.Add( name, activity, activityStep );
      inputPrams.SetLowerLimit( name, 0.0 );
    }else
    {
      inputPrams.Add( nuclide->symbol + "Strength", activity );
    }

    
    if( fitAge && hasOwnAge )
    {
      //We could do a lot better on creating the age range - there must be some way to easily
      //  determine how old an isotope has to get before it essentially doesn't change any more
      //  (prompt HL, etc.).  I guess we could look at the peaks being used to fit for and age them
      //  until their ratios don't change within the available statistical precision.
      //  But for the moment, we'll do something much simpler and use the maximum of either the
      //  longest progenies half-life, or the sum half life of all progeny
      //  \TODO: improve this max decay time estimate; some possibilities are:
      //    - Could probably ignore the parents half-life, or only partially take into account
      //    - For common nuclides could define reasonable fixed values
      //    - Could look at gamma spectrum produced over time, and pick the time when the selected
      //      photopeak ratios change little enough as to not be statistically significant to the
      //      observed data (or even just hard-coded limits).
      auto maxNuclideDecayHL = []( const SandiaDecay::Nuclide * const nuc ) -> double {
        double maxhl = 0.0, sumhlfs = 0.0;
        
        for( auto n : nuc->descendants() )
        {
          if( !n->isStable() )
          {
            sumhlfs += n->halfLife;
            maxhl = std::max( maxhl, n->halfLife );
          }
        }//for( auto n : nuc->descendants() )
        
        //cout << "For nuc=" << nuc->symbol << " maxhl=" << PhysicalUnits::printToBestTimeUnits(maxhl)
        //     << ", sumhlfs=" << PhysicalUnits::printToBestTimeUnits(sumhlfs)
        //     << " - will set max age to " << PhysicalUnits::printToBestTimeUnits(std::max( 7*maxhl, 4*sumhlfs ))
        //     << endl;
        
        //return 100.0*nuc->halfLife;
        return std::max( 7*maxhl, 4*sumhlfs );
      };//maxNuclideDecayHL
      
      double maxAge = -1.0;
      if( ageDefiningNuc == nuclide )
      {
        //We are
        
        for( size_t trialInd = 0; trialInd < answer->numNuclides(); ++trialInd )
        {
          const SandiaDecay::Nuclide *trialNuc = answer->nuclide( int(trialInd) );
          if( trialNuc->atomicNumber == nuclide->atomicNumber )
          {
            const double thisMaxAge = maxNuclideDecayHL( nuclide );
            maxAge = std::max( maxAge, thisMaxAge );
          }
        }//for( loop over all nuclides being fit for )
      }else
      {
        maxAge = maxNuclideDecayHL( nuclide );
      }
      assert( maxAge > 0.0 );
      
      double age = m_sourceModel->age( ison );
      double ageStep = 0.25 * nuclide->halfLife;
      
      // Limit the maximum age to be the larger of ten times the current age, or 200 years.  This
      //  is both to be reasonable in terms of answers we get, and because for really long-lived
      //  isotopes, we could have a max age at this point so large it will cause Minuit to give
      //  NaNs for age, even on first iteration.
      maxAge = std::min( maxAge, std::max(10.0*age, 200.0*PhysicalUnits::year) );
      
      // If the age is currently over 10000 years, it is just really getting unreasonable, so
      //  larger than Minuit can handle, so will impose a tougher 100k year limit, but let grow past
      //  this, but require the user to hit "fit" over and over again.
      if( maxAge > 10000.0*PhysicalUnits::year )
        maxAge = std::max(2.0*age, 10000.0*PhysicalUnits::year);
      
      // But no matter what we'll limit to the age of earth, which at least for a few select
      //  examples tried, Minuit was okay with (it wasnt okay with like 1.2E20 years that some of
      //  the uraniums would give)
      age = std::min( age, 4.543e+9 * PhysicalUnits::year );
      maxAge = std::min( maxAge, 4.543e+9 * PhysicalUnits::year );
      
      ageStep = std::min( ageStep, 0.1*maxAge );
      if( age > 0 )
        ageStep = std::min( 0.1*age, ageStep );
      
      inputPrams.Add( nuclide->symbol + "Age", age, ageStep, 0, maxAge  );
    }else if( hasOwnAge )
    {
      const double age = m_sourceModel->age( ison );
      inputPrams.Add( nuclide->symbol + "Age", age );
    }else  //see if defining nuclide age is fixed, if so use it, else put in negative integer of index of age...
    {
      assert( ageDefiningNuc );
      const int age_defining_index = m_sourceModel->nuclideIndex( ageDefiningNuc );
      inputPrams.Add( nuclide->symbol + "Age", -1.0*(age_defining_index + 1) );
    }
  }//for( int ison = 0; ison < niso; ++ison )
  
  
  //setup the parameters for the shieldings
  for( size_t i = 0; i < shieldings.size(); ++i )
  {
    ShieldingSelect *select = shieldings[i];
    if( select->isGenericMaterial() )
    {
      const double an = select->atomicNumber();
      const double ad = select->arealDensity();
      const bool fitAn = select->fitAtomicNumber();
      const bool fitAD = select->fitArealDensity();

      const string name = "Generic_" + std::to_string(i);
      const double adUnits = PhysicalUnits::g/PhysicalUnits::cm2;

      num_fit_params += fitAn + fitAD;

      if( fitAn )
        inputPrams.Add( name + "_AN", an, std::max(0.1*an,2.5),
                        1.0*MassAttenuation::sm_min_xs_atomic_number,
                        1.0*MassAttenuation::sm_max_xs_atomic_number );
      else
        inputPrams.Add( name + "_AN_FIXED", an );

      if( fitAD )
        inputPrams.Add( name + "_AD", ad, std::max(10.0*adUnits, 0.25*ad), 0.0, 400.0*adUnits );
      else
        inputPrams.Add( name + "_AD", ad );
    }else
    {
      const double thickness = select->thickness();
      bool fitThickness = select->fitThickness();

      std::shared_ptr<Material> mat = select->material();
      string name;
      if( mat )
        name = mat->name + std::to_string(i);
      else
        name = "unspecifiedmat_" + std::to_string(i);

      if( !mat /*|| SpecUtils::iequals_ascii(mat->name, "void")*/ )
        fitThickness = false;

      num_fit_params += fitThickness;

      if( fitThickness )
        inputPrams.Add( name + "_thickness", thickness, std::max(10.0*PhysicalUnits::mm,0.25*thickness), 0, 1000.0*PhysicalUnits::m );
      else
        inputPrams.Add( name + "_thickness", thickness );
#if( USE_CONSISTEN_NUM_SHIELDING_PARS )
      inputPrams.Add( name + "_dummyshielding", 0.0 );
#endif
    }//if( generic material ) / else
  }//for( size_t i = 0; i < shieldings.size(); ++i )

//  if( num_fit_params < 1 )
//    throw runtime_error( "There is nothing being fit for" );
  
  if( m_backgroundPeakSub->isChecked() )
  {
    std::shared_ptr<const SpecMeas> back = m_specViewer->measurment(SpecUtils::SpectrumType::Background);
    typedef std::shared_ptr<const PeakDef> PeakPtr;
    typedef std::deque< std::shared_ptr<const PeakDef> > PeakDeque;
    std::shared_ptr<const PeakDeque > backpeaks;
    
    if( back )
    {
      const auto &displayed = m_specViewer->displayedSamples(SpecUtils::SpectrumType::Background);
      backpeaks = back->peaks( displayed );
    }//if( back )
    
    if( backpeaks && !backpeaks->empty() )
    {
      vector<PeakDef> backgroundpeaks;
      for( const PeakPtr &p : *backpeaks )
        backgroundpeaks.push_back( *p );
      answer->setBackgroundPeaks( backgroundpeaks, back->gamma_live_time() );
    }else
    {
      m_backgroundPeakSub->setUnChecked();
      passMessage( "There are no background peaks defined, not subtracting them",
                   "", WarningWidget::WarningMsgInfo );
    }//if( !peaks || peaks->empty() )
  }//if( m_backgroundPeakSub->isChecked() )
  
  
#if( ALLOW_MASS_FRACTION_FIT )  
  for( size_t i = 0; i < shieldings.size(); ++i )
  {
    ShieldingSelect *select = shieldings[i];
    if( select->fitForMassFractions() )
    {
      //Get the isotopes to fit mass fractions of
      const vector<ShieldingSelect::NucMasFrac> massfracs
                                         = select->sourceNuclideMassFractions();
      vector<const SandiaDecay::Nuclide *> nucstofit
                                         = select->nuclidesToUseAsSources();
      
      std::shared_ptr<Material> mat = select->material();
      if( !mat )
        throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...)"
                             " serious logic error when fitting for mass frac");
      answer->setNuclidesToFitMassFractionFor( mat.get(), nucstofit );
      
      nucstofit = answer->nuclideFittingMassFracFor( mat.get() );
      
      double fracmaterial = 0.0;
      for( const ShieldingSelect::NucMasFrac &nmf : massfracs )
        fracmaterial += nmf.second;
      
      if( fracmaterial <= 0.0 )
      {
        passMessage( "When fitting for mass fractions of source nuclides, the "
                    "sum of the fit for mass fractions equal the sum of the "
                    "inital values, therefore the initial sum of mass "
                    "fractions must be greater than 0.0",
                    "", WarningWidget::WarningMsgHigh );
        throw runtime_error( "Error fitting mass fraction" );
      }//if( fracmaterial <= 0.0 )
      
      
      vector<ShieldingSelect::NucMasFrac> orderedmassfracs;
      const size_t nfitnucs = (nucstofit.empty() ? 0 : nucstofit.size()-1);
      for( size_t j = 0; j < nfitnucs; ++j )
      {
        const SandiaDecay::Nuclide *nuc = nucstofit[j];
        for( size_t k = 0; k < nfitnucs; ++k )
          if( massfracs[k].first == nuc )
            orderedmassfracs.push_back( massfracs[k] );
      }//for( size_t j = 0; i < nfitnucs; ++j )
      
      if( nfitnucs != orderedmassfracs.size() )
        throw runtime_error( "nfitnucs != orderedmassfracs.size()" );
      
      double usedmassfrac = 0.0;
      for( size_t j = 0; j < nfitnucs; ++j )
      {
        const ShieldingSelect::NucMasFrac &nmf = orderedmassfracs[j];
        string name = mat->name + "_" + nmf.first->symbol
                      + "_" + std::to_string(i);
        double val = 0.0;
        const double remaining_frac = (fracmaterial - usedmassfrac);
        if( remaining_frac > nmf.second )
          val = nmf.second / remaining_frac;
        
        usedmassfrac += nmf.second;
        inputPrams.Add( name, val, max(0.1*val,0.01), 0, 1.0 );
        ++num_fit_params;
      }//for( size_t j = 0; i < nmassfrac; ++j )
    }//if( fit for mass fractions to )
  }//for( WWidget *widget : m_shieldingSelects->children() )
#endif
  
  
  if( num_fit_params != inputPrams.VariableParameters() )
    throw runtime_error( "ShieldingSourceDisplay::shieldingFitnessFcn(...): "
                        "there is a serious logic error in this function, "
                        "please let wcjohns@sandia.gov know about this." );
  

  return answer;
}//shared_ptr<PointSourceShieldingChi2Fcn> shieldingFitnessFcn()


void ShieldingSourceDisplay::cancelModelFit()
{
  std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
  if( m_currentFitFcn )
    m_currentFitFcn->cancelFit();
}//void cancelModelFit()


void ShieldingSourceDisplay::updateGuiWithModelFitProgress( std::shared_ptr<ModelFitProgress> progress )
{
  ModelFitProgress status;
  
  {
    std::lock_guard<std::mutex> lock( progress->m );
    status.chi2 = progress->chi2;
    status.numFcnCalls = progress->numFcnCalls;
    status.elapsedTime = progress->elapsedTime;
    status.parameters = progress->parameters;
  }

  status.chi2 = floor(10.0*status.chi2 + 0.5) / 10.0;
  status.elapsedTime = floor(10.0*status.elapsedTime + 0.5) / 10.0;
  
  char buffer[128];
  snprintf( buffer, sizeof(buffer),
            "Best &chi;&sup2;: %.1f, %i fcn calls in %.1fs",
            status.chi2, static_cast<int>(status.numFcnCalls), status.elapsedTime );
  
  m_fitProgressTxt->setText( buffer );
  
  wApp->triggerUpdate();
}//void updateGuiWithModelFitProgress( std::shared_ptr<ModelFitProgress> progress )


void ShieldingSourceDisplay::updateGuiWithModelFitResults( std::shared_ptr<ModelFitResults> results )
{
  WApplication *app = wApp;
  WApplication::UpdateLock applock( app );
  
  if( !applock )
  {
    // Shouldnt ever get here
    cerr << "Failed to get application lock!" << endl;
    return;
  }
  
  // Make sure we trigger a app update
  BOOST_SCOPE_EXIT(app){
    if( app )
      app->triggerUpdate();
  } BOOST_SCOPE_EXIT_END
  
  
  assert( results );
  const ModelFitResults::FitStatus status = results->succesful;
  const vector<ShieldingSelect *> &shieldings = results->shieldings;
  const vector<double> &paramValues = results->paramValues;
  const vector<double> &paramErrors = results->paramErrors;
  const vector<string> &errormsgs = results->errormsgs;
  
  m_fitModelButton->show();
  m_fitProgressTxt->setText("");
  m_fitProgressTxt->hide();
  m_cancelfitModelButton->hide();
  
  m_peakView->enable();
  m_sourceView->enable();
  m_optionsDiv->enable();
  m_addItemMenu->enable();
  m_distanceEdit->enable();
  m_detectorDisplay->enable();
  m_shieldingSelects->enable();
  m_addGenericShielding->enable();
  m_addMaterialShielding->enable();
#if( USE_DB_TO_STORE_SPECTRA )
  m_saveAsNewModelInDb->enable();
#endif
  
  std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
  if( !m_currentFitFcn )
  {
    passMessage( "Programming Logic Error - received model fit results at an invalid time.", "", WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, "Programming Logic Error - received model fit results at an invalid time." );
#endif
    return;
  }//if( !m_currentFitFcn )
  
  for( auto shielding : shieldings )
  {
    if( m_shieldingSelects->indexOf(shielding) < 0 )
    {
      passMessage( "Programming Logic Error - shieldings have changed.", "", WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "Programming Logic Error - received fit results when model was no longer valid." );
#endif
      m_currentFitFcn.reset();
      return;
    }
  }//for( auto shielding : shieldings )
  
  if( status == ModelFitResults::FitStatus::Invalid )
  {
    passMessage( "<b>Fit to model failed</b>.",
                  "", WarningWidget::WarningMsgHigh );
    for( auto &s : errormsgs )
      passMessage( s, "", WarningWidget::WarningMsgHigh );
    m_currentFitFcn.reset();
    return;
  }
  
  if( status == ModelFitResults::FitStatus::InterMediate )
  {
    passMessage( "Intermediate Fit status not handled yet.",
                "", WarningWidget::WarningMsgHigh );
    return;
  }
  
  for( auto &s : errormsgs )
    passMessage( s + "<br />Using fit solution anyway." ,
                   "", WarningWidget::WarningMsgHigh );
  
  try
  {
    const size_t nnucs = m_currentFitFcn->numNuclides();
    
    //Go through and set the ages and activities fit for
    for( size_t nucn = 0; nucn < nnucs; ++nucn )
    {
      const SandiaDecay::Nuclide *nuc = m_currentFitFcn->nuclide( int(nucn) );
      
      const double age = m_currentFitFcn->age( nuc, paramValues );
      //XXX - if multiple shieldings or point sources have this isotope, then
      //      this activity below will be the sum of them
      const double activity = m_currentFitFcn->activity( nuc, paramValues );
      
      char actStr[64], ageStr[64];
      snprintf( actStr, sizeof(actStr), "%f bq", (activity/PhysicalUnits::becquerel) );
      snprintf( ageStr, sizeof(ageStr), "%f s", (age/PhysicalUnits::second) );
      
      //      cerr << "activity=" << activity << "-->" << actStr << "-->" << PhysicalUnits::stringToActivity(actStr) << endl;
      
      const int ison = m_sourceModel->row( nuc );
      WModelIndex ageIndex = m_sourceModel->index( ison, SourceFitModel::kAge );
      WModelIndex ageUncerIndex = m_sourceModel->index( ison, SourceFitModel::kAgeUncertainty );
      WModelIndex actIndex = m_sourceModel->index( ison, SourceFitModel::kActivity );
      WModelIndex actUncertIndex = m_sourceModel->index( ison, SourceFitModel::kActivityUncertainty );
      
      if( m_sourceModel->fitAge(ison) )
      {
        bool success = m_sourceModel->setData( ageIndex, WString(ageStr) );
        if( !success )
        {
          stringstream msg;
          msg << "An invalid age was calculated for " << nuc->symbol
          << ", other results may be invalid to";
          passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
        }//if( IsNan(p.second) || IsInf(p.second) )
        
        const double ageUncert = m_currentFitFcn->age( nuc, paramErrors );
        
        char ageUncertStr[64];
        snprintf( ageUncertStr, sizeof(ageUncertStr), "%f s", (ageUncert/PhysicalUnits::second) );
        m_sourceModel->setData( ageUncerIndex, WString(ageUncertStr) );
      }else if( m_sourceModel->ageUncert(ison) >= 0.0 )
      {
        m_sourceModel->setData( ageUncerIndex, boost::any() );
      }//fit( age ) / else
      
      if( m_sourceModel->fitActivity(ison) )
      {
        bool success = m_sourceModel->setData( actIndex, WString(actStr) );
        
        if( !success )
        {
          stringstream msg;
          msg << "An invalid activity was calculated for " << nuc->symbol
          << ", other results may be invalid to";
          passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
        }//if( IsNan(p.second) || IsInf(p.second) )
        
        const double activityUncert = m_currentFitFcn->activity( nuc, paramErrors );
        char actUncertStr[64];
        snprintf( actUncertStr, sizeof(actUncertStr), "%f bq", (activityUncert/PhysicalUnits::becquerel) );
        
        success = m_sourceModel->setData( actUncertIndex, WString(actUncertStr) );
        if( !success )
        {
          stringstream msg;
          msg << "Calculated activity uncertainty for " << nuc->symbol
          << " is invalid, other results may be invalid to";
          passMessage( msg.str(), "", WarningWidget::WarningMsgHigh );
        }//if( IsNan(p.second) || IsInf(p.second) )
      }else if( m_sourceModel->activityUncert(ison) >= 0.0 )
      {
        m_sourceModel->setData( actUncertIndex, boost::any() );
      }
    }//for( int ison = 0; ison < niso; ++ison )
    
    //    double radius = 0.0;
    const int nshieldings = static_cast<int>(shieldings.size());
    for( int i = 0; i < nshieldings; ++i )
    {
      ShieldingSelect *select = shieldings[i];
      if( select->isGenericMaterial() )
      {
        const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
        const double an = m_currentFitFcn->atomicNumber( i, paramValues );
        const double ad = m_currentFitFcn->arealDensity( i, paramValues ) / adUnits;
        
        char buffer[32];
        snprintf( buffer, sizeof(buffer), "%.2f", an );
        select->m_atomicNumberEdit->setText( buffer );
        snprintf( buffer, sizeof(buffer), "%.2g", ad );
        select->m_arealDensityEdit->setText( buffer );
      }else
      {
        WString txt;
        const double thickness = m_currentFitFcn->thickness( i, paramValues );
#ifndef WT_NO_STD_WSTRING
        const double thicknessErr = m_currentFitFcn->thickness( i, paramErrors );
        if( thicknessErr > 0.0 )
          txt = PhysicalUnits::printToBestLengthUnits( thickness, thicknessErr );
        else
#endif
          txt = PhysicalUnits::printToBestLengthUnits( thickness );
        select->m_thicknessEdit->setText( txt );
        //        const double volume = (3.14159265359*4.0/3.0) * ( pow(radius + thickness,3) - pow(radius,3) );
        //        const double density = select->material()->density;
        //        cerr << "From Geometry, material " << select->material()->name << " has "
        //             << "mass " << (volume*density)/PhysicalUnits::gram << " g" << endl;
        //        radius += thickness;
      }//if( genericMaterial ) / else
    }//for( size_t i = 0; i < shieldings.size(); ++i )
    
    
#if( ALLOW_MASS_FRACTION_FIT )
    const vector<const Material *> massfracFitMaterials
    = m_currentFitFcn->materialsFittingMassFracsFor();
    for( int i = 0; i < nshieldings; ++i )
    {
      ShieldingSelect *select = shieldings[i];
      if( select->fitForMassFractions() )
      {
        std::shared_ptr<Material> usrmaterial = select->material();
        
        //Find this material in massfracFitMaterials
        const Material *origMaterial = NULL;
        for( const Material *m : massfracFitMaterials )
        {
          if( m == usrmaterial.get() )
          {
            origMaterial = m;
            break;
          }
        }//for( const Material *m : massfracFitMaterials )
        
        if( !origMaterial )
        {
          cerr << "ShieldingSourceDisplay::updateGuiWithModelFitResults(...)\n\tNo match for usrmaterial in origMaterial" << endl;
          for( const Material *m : massfracFitMaterials )
          {
            if( m->name == usrmaterial->name )
            {
              origMaterial = m;
              break;
            }
          }//for( const Material *m : massfracFitMaterials )
        }//if( !origMaterial )
        
        if( !origMaterial )
          throw runtime_error( "ShieldingSourceDisplay::doModelFit(): logic "
                              "error completing mass fraction fit :(" );
        
        //        std::shared_ptr<Material> newmaterial
        //                = m_currentFitFcn->variedMassFracMaterial( origMaterial, paramValues );
        //        select->updateMassFractionDisplays( newmaterial );
        
        const vector<const SandiaDecay::Nuclide *> &fitnucs
        = m_currentFitFcn->nuclideFittingMassFracFor( origMaterial );
        vector<double> massfracs;
        
        //XXX
        for( const SandiaDecay::Nuclide *nuc : fitnucs )
        {
          double frac = m_currentFitFcn->massFraction( origMaterial, nuc, paramValues );
          massfracs.push_back( frac );
        }//for( const SandiaDecay::Nuclide *nuc : fitnucs )
        
        for( size_t i = 0; i < massfracs.size(); ++i )
        {
          select->setMassFraction( fitnucs[i], massfracs[i] );
          //          cerr << "Setting " << fitnucs[i]->symbol << " mass fraction to "
          //               << massfracs[i] << endl;
        }
        
        select->updateMassFractionDisplays( usrmaterial );
        
        cerr << "ShieldingSourceDisplay::updateGuiWithModelFitResults(...)\n\tThis whole fitting for mass fractions of "
        << "nuclides is sketch: for one thing, I dont like that the "
        << "Material object of ShieldingSourceDisplay is mutable, and it "
        << "appears this is necessary; for another, I feel like just "
        << "everything to do with this is brittle, and it scares me."
        << " Please consider really cleaning this stuff up!"
        << endl << endl;
      }//if( select->fitForMassFractions() )
    }//for( int i = 0; i < nshieldings; ++i )
#endif
    
    updateChi2ChartActual();
    m_chi2ChartNeedsUpdating = false;
    updateCalcLogWithFitResults( m_currentFitFcn, results );
  }catch( std::exception &e )
  {
    passMessage( "Programming issue - caught exception: " + string(e.what())
                 + "<br />Application state may be suspect!", "",
                 WarningWidget::WarningMsgHigh );
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( __func__, ("Programming Issue - caught exception: " + string(e.what())).c_str() );
#endif
  }//try / catch
  
  m_currentFitFcn.reset();
}//void updateGuiWithModelFitResults( std::vector<double> paramValues, paramErrors )



void ShieldingSourceDisplay::doModelFittingWork( const std::string wtsession,
                                          std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams,
                                          std::shared_ptr<ModelFitProgress> progress,
                                          boost::function<void()> progress_fcn,
                                          std::shared_ptr<ModelFitResults> results,
                                          boost::function<void()> update_fcn )
{
  //The self attenuating probing questions are not tested.
  
  assert( results );
  
  results->succesful = ModelFitResults::FitStatus::Invalid;
  
  Chi2FcnShrdPtr chi2Fcn;
  
  {
    std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
    chi2Fcn = m_currentFitFcn;
  }
  
  if( progress_fcn && progress ) //if we are wanted to post updates to the GUI periodically.
  {
    //Create a object that wil be shared by the Chi2 function, and its callback
    //  which will in turn update the variable shared with the function that
    //  gets posted to the GUI thread.  Its a little convoluted, but I think
    //  this lets us better make a consistent handoff of information (e.g., we
    //  can be sure all the member variables of ModelFitProgress are not-changed
    //  by the time the GUI update function gets executed in the Wt event loop).
    auto progressUpdatInfo = std::make_shared<GammaInteractionCalc::PointSourceShieldingChi2Fcn::GuiProgressUpdateInfo>();
  
    //We need to pass in a refernce to progressUpdatInfo to the lamda function
    //  we will assign to progressUpdatInfo->m_guiUpdater, but we need to avoid
    //  creating a circular reference that keeps progressUpdatInfo from ever
    //  getting deleted, so we will pass in a week pointer to m_guiUpdater.
    std::weak_ptr<GammaInteractionCalc::PointSourceShieldingChi2Fcn::GuiProgressUpdateInfo> progress_weak = progressUpdatInfo;
    
    progressUpdatInfo->m_guiUpdater = [progress_fcn,progress_weak,wtsession,progress](){
      auto infoFromChiFcn = progress_weak.lock();
      if( !infoFromChiFcn )
        return;
    
      vector<double> pars;
      
      {
        //As it stands the thread we are currently in is the only thread to access
        // the object pointed to by infoFromChiFcn, so there is no need to bother
        // with the atomic or mutexs, but because Minuit2 I think has a
        // multithreaded mode we might look into later, I'll keep this stuff
        // thread safe for now.
        std::lock_guard<std::mutex> scoped_lock( infoFromChiFcn->m_bestParsMutex );
        pars = infoFromChiFcn->m_bestParameters;
      }
      
      {
        std::lock_guard<std::mutex> scoped_lock( progress->m );
        progress->chi2 = infoFromChiFcn->m_bestChi2;
        progress->elapsedTime = (infoFromChiFcn->m_currentTime - infoFromChiFcn->m_fitStartTime);
        progress->parameters = pars;
        progress->numFcnCalls = infoFromChiFcn->m_num_fcn_calls;
      }
      
      WApplication *app = WApplication::instance();
      if( app && (app->sessionId() == wtsession) )
        progress_fcn();
      else
        WServer::instance()->post( wtsession, progress_fcn );
    };//define progressUpdatInfo->m_guiUpdater
  
    chi2Fcn->setGuiProgressUpdater( sm_model_update_frequency_ms, progressUpdatInfo );
  }//if( progress_fcn && progress )
  
  
  try
  {
    if( !chi2Fcn )
      throw runtime_error( "Programming logic error - Chi2Function pointer is null." );
    
    //I think inputPrams.VariableParameters() == num_fit_params
    if( inputPrams->VariableParameters() > chi2Fcn->peaks().size() )
    {
      WStringStream msg;
      msg << "You are asking to fit " << int(inputPrams->VariableParameters())
      << " parameters, however there are only " << int(chi2Fcn->peaks().size())
      << " peaks, which leads to this being an under-constrained problem";
      throw runtime_error( msg.str() );
    }//if( num_fit_params > peaks.size() )
    
    if( inputPrams->VariableParameters() < 1 )
      throw runtime_error( "No parameters are selected for fitting." );
    
    chi2Fcn->fittingIsStarting( sm_max_model_fit_time_ms );
    
    ROOT::Minuit2::MnUserParameterState inputParamState( *inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
    ROOT::Minuit2::MnMinimize fitter( *chi2Fcn, inputParamState, strategy );
    
    //cout << "Parameters are: {";
    //for( const auto &par : inputPrams->Parameters() )
    //  cout << par.Name() << ", ";
    //cout << endl;
    
    const double tolerance = 2.0*inputPrams->VariableParameters();
    unsigned int maxFcnCall = 50000;  //default minuit2: 200 + 100 * npar + 5 * npar**2
    
    ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
    
    //For some reason the fit to atomic number of generic shielding is horrible!
    //  I'm probably screwing something up in the setup for Minuit, but as a
    //  work around, lets manually coursely scan AN, and then use the best found
    //  AN to then perform another fit to fine tune it.
    //Not a lot of validation went into this, but it does seem to work better.
    vector<string> fit_generic_an;
    for( auto &par : inputPrams->Parameters() )
    {
      const string &name = par.GetName();
      if( (name.find("Generic_")!=string::npos)
         && (name.find("_AN")!=string::npos)
         && (name.find("_FIXED")==string::npos)
         && !par.IsFixed() )  //Looks to be bug in Minuit, so we need the above "_FIXED"
        fit_generic_an.push_back( name );
    }//for( auto &par : inputPrams->Parameters() )
    
    if( !fit_generic_an.empty() )  //should add in a check that we arnt doing self attenuating sources.
    {
      //Re-perform the fit scanning in atomic number to find best.
      const double orig_chi2 = chi2Fcn->DoEval( fitter.Params() );
      
      for( const auto &parname : fit_generic_an )
      {
        vector<double> best_pars = fitter.Params();
        double best_an = inputPrams->Value( parname );
        double best_chi2 = chi2Fcn->DoEval( fitter.Params() );
        
        for( double an = 1.0; an < 101.0; an += 5.0 )
        {
          auto testpar = *inputPrams;
          testpar.SetValue(parname, an);
          testpar.Fix(parname);
          
          ROOT::Minuit2::MnUserParameterState anInputParam( testpar );
          ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
          ROOT::Minuit2::MnMinimize anfitter( *chi2Fcn, anInputParam, strategy );
          
          ROOT::Minuit2::FunctionMinimum anminimum = anfitter( maxFcnCall, tolerance );
          if( !anminimum.IsValid() )
            anminimum = anfitter( maxFcnCall, tolerance );
          
          const double this_chi2 = chi2Fcn->DoEval( anfitter.Params() );
          if( this_chi2 < best_chi2 )
          {
            best_pars = anfitter.Params();
            best_chi2 = this_chi2;
            best_an = an;
          }
        }//for( double an = 1.0; an < 101.0; an += 2.5 )
        
        fitter.SetValue( parname.c_str(), best_an );
        fitter.Fix( parname.c_str() );
      }//for( const auto &parname : fit_generic_an )
      
      for( const auto &parname : fit_generic_an )
        fitter.Release( parname.c_str() );
      
      ROOT::Minuit2::FunctionMinimum postnminimum = fitter( maxFcnCall, tolerance );
      const double post_best_chi2 = chi2Fcn->DoEval( fitter.Params() );
      
      if( post_best_chi2 < orig_chi2 )
      {
        minimum = postnminimum;
        cout << "Found a new minumum of post_best_chi2=" << post_best_chi2 << ", vs old orig_chi2=" << orig_chi2 << endl;
      }//if( post_best_chi2 < orig_chi2 )
      else
      {
        cout << "Old minumum was better: post_best_chi2=" << post_best_chi2 << ", vs old orig_chi2=" << orig_chi2 << endl;
      }
    }else if( !minimum.IsValid() )
    {
      //Try two more times to get a valid fit... a stupid hack
      for( int i = 0; !minimum.IsValid() && i < 2; ++i )
        minimum = fitter( maxFcnCall, tolerance );
    }//if( !fit_generic_an.empty() ) / else
    
    
    
    const ROOT::Minuit2::MnUserParameters &fitParams = minimum.UserParameters();
    
    std::lock_guard<std::mutex> lock( results->m_mutex );
    
    if( !minimum.IsValid() )
    {
      stringstream msg;
      msg << "Fit status is not valid:";
      if( minimum.HasMadePosDefCovar() )
        msg << "<br />&nbsp;&nbsp;-Covariance matrix forced positive-definit";
      if( !minimum.HasAccurateCovar() )
        msg << "<br />&nbsp;&nbsp;-Does not have accurate covariance matrix";
      if( minimum.HasReachedCallLimit() )
        msg << "<br />&nbsp;&nbsp;-Optimization reached call limit.";
      if( !minimum.HasValidCovariance() )
        msg << "<br />&nbsp;&nbsp;-Did not have valid covariance,";
      if( !minimum.HasValidParameters() )
        msg << "<br />&nbsp;&nbsp;-Invalid fit parameters.";
      if( minimum.IsAboveMaxEdm() )
        msg << "<br />&nbsp;&nbsp;-The estimated distance to minimum too large.";
      
      results->errormsgs.push_back( msg.str() );
    }//if( !minimum.IsValid() )
    
    results->succesful = ModelFitResults::FitStatus::Final;
    results->paramValues = fitParams.Params();
    results->paramErrors = fitParams.Errors();
    results->edm = minimum.Edm();
    results->num_fcn_calls = minimum.NFcn();
    results->chi2 = minimum.Fval();  //chi2Fcn->DoEval( results->paramValues );
  }catch( exception &e )
  {
    std::lock_guard<std::mutex> lock( results->m_mutex );
    results->succesful = ModelFitResults::FitStatus::Invalid;
    results->errormsgs.push_back( "Fit not performed: " + string(e.what()) );
  }// try / catch
  
  chi2Fcn->fittindIsFinished();
  
  if( update_fcn )
  {
    WApplication *app = WApplication::instance();
    if( app && (app->sessionId() == wtsession) )
      update_fcn();
    else
      Wt::WServer::instance()->post( wtsession, update_fcn );
  }//if( update_fcn )
}//void doModelFittingWork( std::shared_ptr<ROOT::Minuit2::MnUserParameters> inputPrams )


std::shared_ptr<ShieldingSourceDisplay::ModelFitResults> ShieldingSourceDisplay::doModelFit( const bool fitInBackground )
{
  try
  {
    checkAndWarnZeroMassFraction();
  }catch( std::exception &e )
  {
    passMessage( e.what() + string("<br />Fit not performed."),
                "", WarningWidget::WarningMsgHigh );
    return nullptr;
  }
  
  //Should make a progress dialog
  
  m_modifiedThisForeground = true;
  
  try
  {
    checkDistanceAndThicknessConsistent();
  }catch( exception &e )
  {
    passMessage( e.what() + string("<br />Fit not performed."),
                 "", WarningWidget::WarningMsgHigh );
    return nullptr;
  }//try / catch
  
  Chi2FcnShrdPtr chi2Fcn;
  vector<string> errormsgs;
  vector<ShieldingSelect *> shieldings;
  
  //make sure fitting for at least one nuclide:
  
  auto inputPrams = make_shared<ROOT::Minuit2::MnUserParameters>();
  
  try
  {
    chi2Fcn = shieldingFitnessFcn( shieldings, *inputPrams );
  }catch( std::exception &e )
  {
    passMessage( e.what() + string("<br />Fit not performed (couldnt make Chi2Fcn)."),
                "", WarningWidget::WarningMsgHigh );
    return nullptr;
  }//try / catch
  
  //To Do:
  //  -Add a timer so that after a max amount of time, the calc is canceled (Wt::WTimer or asio deadline_timer)
  //   -see void #include <boost/asio/deadline_timer.hpp> from WIOService
  //  -Maybe, have a progress bar or update function
  //    -See WTimer
  m_fitModelButton->hide();
  m_fitProgressTxt->show();
  m_fitProgressTxt->setText("");
  m_cancelfitModelButton->show();
  
  m_peakView->disable();
  m_sourceView->disable();
  m_optionsDiv->disable();
  m_addItemMenu->disable();
  m_distanceEdit->disable();
  m_detectorDisplay->disable();
  m_shieldingSelects->disable();
  m_addGenericShielding->disable();
  m_addMaterialShielding->disable();
#if( USE_DB_TO_STORE_SPECTRA )
  m_saveAsNewModelInDb->disable();
#endif

  //Need to disable "All Peaks", Detector, Distance, and "Material", and "Generic"
  
  {
    std::lock_guard<std::mutex> lock( m_currentFitFcnMutex );
    m_currentFitFcn = chi2Fcn;
  }
  
  auto results = make_shared<ModelFitResults>();
  results->succesful = ModelFitResults::FitStatus::Invalid;
  results->shieldings = shieldings;
  
  auto progress = std::make_shared<ModelFitProgress>();
  progress->chi2 = std::numeric_limits<double>::max();
  progress->elapsedTime = 0.0;
  progress->numFcnCalls = 0.0;
  boost::function<void()> progress_updater = wApp->bind( boost::bind( &ShieldingSourceDisplay::updateGuiWithModelFitProgress, this, progress ) );
  
  //Wrap the GUI update with WAPplication::bind in case this
  //  ShieldingSourceDisplay widget gets deleted before the computation is over
  boost::function<void()> gui_updater = wApp->bind( boost::bind( &ShieldingSourceDisplay::updateGuiWithModelFitResults, this, results ) );
  
  const string sessionid = wApp->sessionId();
  if( fitInBackground )
  {
    Wt::WServer *server = Wt::WServer::instance();
    server->ioService().post( boost::bind( &ShieldingSourceDisplay::doModelFittingWork,
                            this, sessionid, inputPrams, progress, progress_updater, results, gui_updater ) );
  }else
  {
    doModelFittingWork( sessionid, inputPrams, progress, progress_updater, results, gui_updater );
  }
  
  return results;
}//void doModelFit()



void ShieldingSourceDisplay::updateCalcLogWithFitResults(
                                 ShieldingSourceDisplay::Chi2FcnShrdPtr chi2Fcn,
                                    std::shared_ptr<ModelFitResults> results )
{
  assert( results );
  const std::vector<double> &params = results->paramValues;
  const std::vector<double> &errors = results->paramErrors;
  
  if( m_calcLog.size() && m_calcLog.back() == ns_no_uncert_info_txt )
    m_calcLog.erase( m_calcLog.end()-1, m_calcLog.end() );
  
  try
  {
    
  
#if(ALLOW_MASS_FRACTION_FIT)
  for( const Material *mat : chi2Fcn->materialsFittingMassFracsFor() )
  {
    stringstream msg;
    msg << "Shielding material " << mat->name << " fit mass fractions for "
    << "isotopes:";
    const vector<const SandiaDecay::Nuclide *> &nucs
    = chi2Fcn->nuclideFittingMassFracFor( mat );
    for( const SandiaDecay::Nuclide *n : nucs )
    {
      const double frac = chi2Fcn->massFraction( mat, n, params );
      const double df = chi2Fcn->massFractionUncert( mat, n, params, errors );
      
      msg << " " << n->symbol << "(massfrac=" << frac << "+-" << df << "),";
    }//for( const SandiaDecay::Nuclide *n : nucs )
    
    m_calcLog.push_back( msg.str() );
  }//for( const Material *mat : chi2Fcn->materialsFittingMassFracsFor() )
#endif
    
  {//begin add chi2 line
    stringstream msg;
    msg << "It took " << results->num_fcn_calls
        << " solution trials to reach chi2=" << results->chi2
        << " with an estimated distance to minumum of " << results->edm;
    m_calcLog.push_back( msg.str() );
  }//end add chi2 line
    
  //Need to list fit parameters and uncertainties here
  const size_t nnuc = chi2Fcn->numNuclides();
  for( size_t nucn = 0; nucn < nnuc; ++nucn )
  {
    const SandiaDecay::Nuclide *nuc = chi2Fcn->nuclide( nucn );
    if( nuc )
    {
      const bool useCi = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", m_specViewer );
      const double act = chi2Fcn->activity( nuc, params );
      const string actStr = PhysicalUnits::printToBestActivityUnits( act, 2, useCi );
      
      const double actUncert = chi2Fcn->activity( nuc, errors );
      const string actUncertStr = PhysicalUnits::printToBestActivityUnits( actUncert, 2, useCi );
      
      const double mass = act / nuc->activityPerGram();
      const std::string massStr = PhysicalUnits::printToBestMassUnits( mass, 2, 1.0 );
      
      const double age = chi2Fcn->age( nuc, params );
      const double ageUncert = chi2Fcn->age( nuc, errors );
      const string ageStr = PhysicalUnits::printToBestTimeUnits( age, 2 );
      const string ageUncertStr = PhysicalUnits::printToBestTimeUnits( ageUncert, 2 );
      
      
      const bool isShielding = chi2Fcn->isActivityDefinedFromShielding(nuc);
      
      stringstream msg;
      msg << nuc->symbol << " fit activity " << actStr << " (mass: " << massStr
      << ") with uncertanty " << actUncertStr << " ("
      << floor(0.5 + 10000*actUncert/act)/100.0 << "%)";
      
      if( ageUncert <= DBL_EPSILON )
        msg << " at assumed age " << ageStr;
      else
        msg << " with age " << ageStr << "+- " << ageUncertStr;
      
      msg << (isShielding ? ", a self attenuating source." : ".");
      
      m_calcLog.push_back( msg.str() );
    }//if( nuc )
  }//for( size_t nucn = 0; nucn < nnuc; ++nucn )
  
  const int nmat = static_cast<int>( chi2Fcn->numMaterials() );
  for( int matn = 0; matn < nmat; ++matn )
  {
    stringstream msg;
    if( chi2Fcn->isSpecificMaterial(matn) )
    {
      const double thickness = chi2Fcn->thickness( matn, params );
      const double thicknessUncert = chi2Fcn->thickness( matn, errors );
      
      const Material *mat = chi2Fcn->material( matn );
      if( mat )
      {
        const double density = mat->density * PhysicalUnits::cm3 / PhysicalUnits::gram;
        msg << mat->name << " has density " << std::setprecision(3) << density
            << "g/cm3 ";
        
        if( thicknessUncert > DBL_EPSILON )
        {
#ifndef WT_NO_STD_WSTRING
           msg << "and fit thickness "
              << WString(PhysicalUnits::printToBestLengthUnits( thickness, thicknessUncert )).toUTF8()
              << "." << endl;
#else
          msg << "and fit thickness "
          << PhysicalUnits::printToBestLengthUnits( thickness ) << " (+-"
          << PhysicalUnits::printToBestLengthUnits( thicknessUncert ) << ")."
          << endl;
#endif
        }else
        {
          msg << "and has fixed thickness "
              << PhysicalUnits::printToBestLengthUnits( thickness )
              << "." << endl;
        }
      }//if( mat )
    }else
    {
      const double adUnits = PhysicalUnits::gram / PhysicalUnits::cm2;
      const double ad = chi2Fcn->arealDensity( matn, params ) / adUnits;
      const double adUncert = chi2Fcn->arealDensity( matn, errors ) / adUnits;
      
      const double an = chi2Fcn->atomicNumber( matn, params );
      const double anUncert = chi2Fcn->atomicNumber( matn, errors );
      
      msg << std::setprecision(3) << "Shielding " << matn+1 << " has "
      << "AtomicNumber=" << an;
      if( anUncert > DBL_EPSILON )
        msg << " (+-" << anUncert << ")";
      msg << " and ArealDensity=" << ad;
      if( adUncert > DBL_EPSILON )
        msg << " (+-" << adUncert << ")";
      msg << " g/cm2";
    }//if( specific ) / else( generic material )
    
    m_calcLog.push_back( msg.str() );
  }//for( size_t matn = 0; matn < nmat; ++matn )
  }catch( std::exception & )
  {
    m_calcLog.push_back( "There was an error and log may not be complete." );
  }
}//updateCalcLogWithFitResults(...)









