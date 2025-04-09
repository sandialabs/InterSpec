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

#include <limits>
#include <string>
#include <vector>
#include <sstream>

#include <boost/math/tools/roots.hpp>

#include <Wt/WMenu>
#include <Wt/WText>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WMenuItem>
#include <Wt/WLineEdit>
#include <Wt/WComboBox>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WApplication>
#include <Wt/WButtonGroup>
#include <Wt/WEnvironment>
#include <Wt/WRadioButton>
#include <Wt/WStackedWidget>
#include <Wt/WDoubleValidator>
#include <Wt/WSuggestionPopup>
#include <Wt/WRegExpValidator>


#include "Math/IFunction.h"
#include "Minuit2/FCNBase.h"
//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/CombinedMinimizer.h"

#include "SpecUtils/Filesystem.h"
#include "SpecUtils/StringAlgo.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/AppUtils.h"
#include "InterSpec/DoseCalc.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/GammaXsGui.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/WarningWidget.h"
#include "InterSpec/DoseCalcWidget.h"
#include "InterSpec/GadrasSpecFunc.h"
#include "InterSpec/ShieldingSelect.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/NativeFloatSpinBox.h"
#include "InterSpec/NuclideSourceEnter.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/IsotopeSelectionAids.h"
#include "InterSpec/IsotopeNameFilterModel.h"
#include "InterSpec/PhysicalUnitsLocalized.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#endif

using namespace Wt;
using namespace std;

#define MU_CHARACTER "\xCE\xBC"

namespace
{
  void right_select_item( WMenu *menu, WMenuItem *item )
  {
    menu->select( item );
    item->triggered().emit( item ); //
  }
  
  class FitShieldingAdChi2
  : public ROOT::Minuit2::FCNBase
  {
    const vector<float> &m_energies;
    const vector<float> &m_intensities;
    const float m_atomic_number;
    const double m_user_entered_dose;
    const double m_distance;
    const GadrasScatterTable &m_scatter;

    
  public:
    
    FitShieldingAdChi2( const vector<float> &energies,
                        const vector<float> &intensities,
                        const float atomic_number,
                        const double user_entered_dose,
                        const double distance,
                        const GadrasScatterTable &scatter )
      : m_energies( energies ),
        m_intensities( intensities ),
        m_atomic_number( atomic_number ),
        m_user_entered_dose( user_entered_dose ),
        m_distance( distance ),
        m_scatter( scatter )
    {
    }
    
    virtual double Up() const { return 1.0; }
    size_t nfitPars() const { return 1; }
    
    virtual double operator()( const std::vector<double> &params ) const
    {
      assert( params.size() == 1 );
      try
      {
        const double ad = params[0] * PhysicalUnits::g / PhysicalUnits::cm2;
        const double dose = DoseCalc::gamma_dose_with_shielding( m_energies, m_intensities, ad, m_atomic_number, m_distance, m_scatter );
        
        //return difference in micro-rem per hour from the target
        return fabs(dose - m_user_entered_dose) * 1000000.0 * PhysicalUnits::hour / PhysicalUnits::rem;
      }catch( std::exception & )
      {
        // Can happen if ad is outside of [0,240]
      }
      
      return std::numeric_limits<double>::max();
    }
  };//class FitShieldingAdChi2

  
  
  double fit_ad( const vector<float> &energies, const vector<float> &intensities,
                const float atomic_number, const double user_entered_dose,
                const double distance,
                const GadrasScatterTable &scatter )
  {
    const double dose_no_shielding = DoseCalc::gamma_dose_with_shielding( energies, intensities, 0.0, atomic_number, distance, scatter );
    
    if( dose_no_shielding < user_entered_dose )
      throw runtime_error( "Dose from source specified is less than"
                          " entered dose, even with out shielding." );
    
    FitShieldingAdChi2 chi2Fcn( energies, intensities,
                                atomic_number, user_entered_dose, distance,
                                scatter );
    
    ROOT::Minuit2::MnUserParameters inputPrams;
    inputPrams.Add( "AD", 20, 5 );
    inputPrams.SetLowerLimit( "AD", 0.0 );
    inputPrams.SetUpperLimit( "AD", 500.0 );

    ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 1 ); //0 low, 1 medium, >=2 high
    
    const unsigned int maxFcnCall = 1000;
    const double tolerance = 1.0;
    ROOT::Minuit2::CombinedMinimizer fitter;
    ROOT::Minuit2::FunctionMinimum minimum
    = fitter.Minimize( chi2Fcn, inputParamState,
                      strategy, maxFcnCall, tolerance );
    
    //Not sure why Minuit2 doesnt like converging on the minumum verry well, but
    //  rather than showing the user an error message, we'll give it anither try
    if( minimum.IsAboveMaxEdm() )
    {
      ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
      minimum = fitter( maxFcnCall, tolerance );
    }//if( minimum.IsAboveMaxEdm() )
    
    if( !minimum.IsValid() )
      throw runtime_error( WString::tr("dcw-err-failed-fit-AD").toUTF8() );
    
    const ROOT::Minuit2::MnUserParameters params = minimum.UserState().Parameters();
    const vector<double> pars = params.Params();
    cerr << "Fit " << pars[0] << " g/cm2 with EDM " << minimum.Edm() << endl;
    
    if( pars[0] > 235.0 )
      throw runtime_error( "Over 235 g/cm2 shielding is required - can not compute." );
    
    return pars[0] * PhysicalUnits::g / PhysicalUnits::cm2;
  }//double fit_ad()
}//namespace



DoseCalcWindow::DoseCalcWindow( MaterialDB *materialDB,
                                Wt::WSuggestionPopup *materialSuggestion,
                                InterSpec *viewer )
: AuxWindow( WString::tr("window-title-dose-calc"),
            (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletNotFullScreen)
             | AuxWindowProperties::SetCloseable
             | AuxWindowProperties::DisableCollapse) )
{
  rejectWhenEscapePressed( true );
  
  m_dose = new DoseCalcWidget( materialDB, materialSuggestion, viewer, contents() );
  m_dose->setHeight( WLength(100,WLength::Percentage) );
  
  AuxWindow::addHelpInFooter( footer(), "dose-dialog" );
  
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton();
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://dose/" + Wt::Utils::urlEncode(m_dose->encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("dcw-qr-tool-state-title"),
                                 WString::tr("dcw-qr-tool-state-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
  if( !viewer->isPhone() )
    footer()->addWidget( qr_btn );
#endif //USE_QR_CODES

  
  WPushButton *closeButton = addCloseButtonToFooter( WString::tr("Close") );
  closeButton->clicked().connect( this, &AuxWindow::hide );
  
#if( USE_QR_CODES )
  if( viewer->isPhone() )
    footer()->addWidget( qr_btn );
#endif
  
  show();
  
  // If we are loading this widget, as we  are creating the InterSpec session,
  //  the screen width and height wont be available, so we'll just assume its
  //  big enough, which it should be.
  const int screenW = viewer->renderedWidth();
  const int screenH = viewer->renderedHeight();
  int width = 625, height = 435;
  if( (screenW > 100) && (screenW < width) )
    width = screenW;
  if( (screenH > 100) && (screenH < height) )
    height = screenH;
  
  resizeWindow( width, height );
  
  // But I think this next call should fix things up, even if we do have a tiny screen
  resizeToFitOnScreen();
  
  centerWindowHeavyHanded();
}//GammaXsWindow(...) constructor


DoseCalcWindow::~DoseCalcWindow()
{
}


DoseCalcWidget *DoseCalcWindow::tool()
{
  return m_dose;
}



DoseCalcWidget::DoseCalcWidget( MaterialDB *materialDB,
                                 Wt::WSuggestionPopup *materialSuggestion,
                                 InterSpec *specViewer,
                                 Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
   m_viewer( specViewer ),
   m_materialSuggest( materialSuggestion ),
   m_materialDB( materialDB ),
   m_enterShieldingSelect( NULL ),
   m_answerShieldingSelect( NULL ),
   m_gammaSource( NULL ),
   m_sourceType( NULL ),
   m_gammaSourceDiv( NULL ),
   m_neutronSourceDiv( NULL ),
   m_neutronSourceCombo( NULL ),
   m_activityAnswer( NULL ),
   m_activityEnter( NULL ),
   m_activityEnterUnits( NULL ),
   m_activityAnswerUnits( NULL ),
   m_doseAnswer( NULL ),
   m_doseEnter( NULL ),
   m_doseEnterUnits( NULL ),
   m_doseAnswerUnits( NULL ),
   m_distanceEnter( NULL ),
   m_distanceAnswer( NULL ),
   m_issueTxt( NULL ),
   m_layout( NULL ),
   m_currentCalcQuantity( DoseCalcWidget::NumQuantity ),
   m_stack( NULL ),
   m_menu( nullptr ),
   m_stateUri()
{
  init();
}//DoseCalcWidget constructor


void DoseCalcWidget::init()
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  wApp->useStyleSheet( "InterSpec_resources/DoseCalcWidget.css" );
  m_viewer->useMessageResourceBundle( "DoseCalcWidget" );
      
  addStyleClass( "DoseCalcWidget" );
  m_layout = new WGridLayout( this );
  m_layout->setContentsMargins( 0, 0, 0, 0 );
  
  try
  {
    string continuumData = SpecUtils::append_path( InterSpec::staticDataDirectory(), "GadrasContinuum.lib" );
    
    m_scatter.reset( new GadrasScatterTable( continuumData ) );
  }catch( std::exception &e )
  {
    WString msg = "<div><b>Error initializing resources:</b></div><div>";
//    msg += Wt::Utils::htmlEncode( e.what(), Wt::Utils::EncodeNewLines );
    msg +=  Wt::WWebWidget::escapeText(	WString(e.what()), true );
    msg += "</div>";
    
    m_layout->addWidget( new WText(msg), 0, 0, AlignCenter | AlignMiddle );
    
    return;
  }//try / catch
  
  const bool isPhone = m_viewer->isPhone();
  int w = m_viewer->renderedWidth();
  int h = m_viewer->renderedHeight();
  if(  m_viewer->isMobile() && (w < 100) )
  {
    w = wApp->environment().screenWidth();
    h = wApp->environment().screenHeight();
  }
  
  const bool narrowLayout = ((w > 100) && (w < 450));
  const bool showToolTips = UserPreferences::preferenceValue<bool>( "ShowTooltips", m_viewer );
  const bool useBq = UserPreferences::preferenceValue<bool>( "DisplayBecquerel", InterSpec::instance() );
  
  WContainerWidget *enterDiv = new WContainerWidget();
  WContainerWidget *answerDiv = new WContainerWidget();
  WGridLayout *enterLayout = new WGridLayout( enterDiv );
  WGridLayout *answerLayout = new WGridLayout( answerDiv );
  
  enterDiv->addStyleClass( "DoseEnterDiv" );
  answerDiv->addStyleClass( "DoseAnswerDiv" );
  
  WContainerWidget *workDiv = new WContainerWidget();
  WGridLayout *workLayout = new WGridLayout( workDiv );

  if( narrowLayout )
  {
    addStyleClass( "NarrowDose" );
    answerLayout->setContentsMargins( 5, 1, 5, 1 );
    enterLayout->setContentsMargins( 5, 1, 5, 1 );
    
    workLayout->addWidget( enterDiv, 0, 0 );
    workLayout->addWidget( answerDiv, 1, 0 );
    workLayout->setColumnStretch( 0, 1 );
    workLayout->setRowStretch( 0, 1 );
    workLayout->setRowStretch( 1, 1 );
  }else
  {
    answerLayout->setContentsMargins( 9, 1, 9, 5 );
    enterLayout->setContentsMargins( 9, 1, 9, 5 );
    
    workLayout->addWidget( enterDiv, 0, 0 );
    workLayout->addWidget( answerDiv, 0, 1 );
    workLayout->setColumnStretch( 0, 1 );
    workLayout->setColumnStretch( 1, 1 );
  }
  
  workLayout->setContentsMargins( 0, 0, 0, 0 );
  workLayout->setVerticalSpacing( 0 );
  workLayout->setHorizontalSpacing( 0 );
  
  WText *coltxt = 0;
  if( !isPhone )
  {
    coltxt = new WText( WString::tr("dcw-inputs"), enterDiv );
    coltxt->addStyleClass( "DoseColLabel" );
    coltxt->setInline( false );
    enterLayout->addWidget( coltxt, 0, 0, AlignMiddle | AlignTop );
  }//if( !isPhone )
  
  coltxt = new WText( WString::tr("dcw-answer") );
  coltxt->addStyleClass( "DoseColLabel" );
  coltxt->setInline( false );
  answerLayout->addWidget( coltxt, 0, 0, AlignMiddle | AlignTop );
  
  
  WContainerWidget *introDiv = new WContainerWidget();
  introDiv->addStyleClass( "DoseIntroDiv" );
  WGridLayout *intoTxtLayout = new WGridLayout( introDiv );
  const char * const instr_key = narrowLayout ? "dcw-intro-instructions-vert-phone" : "dcw-intro-instructions";
  WText *txt = new WText( WString::tr(instr_key) );
  txt->addStyleClass( "DoseIntroTxtMain" );
  
  intoTxtLayout->addWidget( txt, 0, 0, AlignMiddle | AlignCenter );
  
  txt = new WText( WString::tr("dcw-intro-calc-desc") );
  txt->addStyleClass( "DoseIntroTxtDetail" );
  
  intoTxtLayout->addWidget( txt, 1, 0, AlignCenter );
  intoTxtLayout->setRowStretch( 0, 10 );
  
  m_menu = new WMenu();
  m_menu->addStyleClass( "DoseCalcSideMenu" );
  m_stack = new WStackedWidget();
//  m_stack->setTransitionAnimation( WAnimation(WAnimation::Fade, WAnimation::Linear, 500), false );
  m_stack->addWidget( introDiv );
  m_stack->addWidget( workDiv );
  
  
  if( narrowLayout )
  {
    m_menu->addStyleClass( "VerticalNavMenuPhone HeavyNavMenuPhone HorizontalMenu" );
    
    m_layout->addWidget( m_menu, 0, 0 );
    m_layout->addWidget( m_stack, 1, 0 );
    m_layout->setColumnStretch( 0, 1 );
    m_layout->setRowStretch( 1, 1 );
  }else
  {
    // Phone is horizontal
    m_menu->addStyleClass( isPhone ? "VerticalNavMenuPhone HeavyNavMenuPhone SideMenuPhone" 
                                    : "VerticalNavMenu HeavyNavMenu SideMenu" );
    m_layout->addWidget( m_menu, 0, 0 );
    m_layout->addWidget( m_stack, 0, 1 );
    m_layout->setColumnStretch( 1, 1 );
  }//if( narrow ) / else horizontal
  
  //First need to put source capability in
  {
    WContainerWidget *sourcDiv = new WContainerWidget();
    sourcDiv->addStyleClass( "DoseEnterEl" );
    
    if( !isPhone )
    {
      txt = new WText( WString::tr("dcw-source-label"), sourcDiv );
      txt->addStyleClass( "DoseEnterInd" );
      txt->setInline( false );
    }//if( !isPhone )
    
    m_sourceType = new WButtonGroup( sourcDiv );
    m_sourceType->checkedChanged().connect( this, &DoseCalcWidget::handleSourceTypeChange );
    
    WContainerWidget *buttonDiv = new WContainerWidget( sourcDiv );
    WRadioButton *gammaButton = new Wt::WRadioButton( WString::tr("Gamma"), buttonDiv );
    gammaButton->setMargin( 7, Wt::Right );
    m_sourceType->addButton( gammaButton, 0 );
    
    WRadioButton *neutronButton = new Wt::WRadioButton( WString::tr("dcw-neutron"), buttonDiv );
    m_sourceType->addButton( neutronButton, 1 );
    neutronButton->disable();
    neutronButton->setToolTip( WString::tr("dcw-tt-neut-not-imp") );
    neutronButton->setAttributeValue( "style", "color: grey;" );
    
    m_sourceType->setSelectedButtonIndex( 0 );
    
    enterLayout->addWidget( sourcDiv, 1, 0, AlignMiddle );
    
    m_gammaSourceDiv = new WContainerWidget( sourcDiv );
    m_neutronSourceDiv = new WContainerWidget( sourcDiv );
    m_neutronSourceDiv->hide();
    m_neutronSourceCombo = new WComboBox( m_neutronSourceDiv );
    m_neutronSourceCombo->addStyleClass( "DoseAnswerCombo" );
    m_neutronSourceCombo->addItem( "Ca252" );
    m_neutronSourceCombo->addItem( "WGPu" );
    m_neutronSourceCombo->addItem( "U238" );
    m_neutronSourceCombo->activated().connect( this, &DoseCalcWidget::updateResult );
    
    m_gammaSource = new NuclideSourceEnter( true, showToolTips, m_gammaSourceDiv );
    m_gammaSource->changed().connect( this, &DoseCalcWidget::updateResult );
  }
  
  
  for( Quantity i = Quantity(0); i < NumQuantity; i = Quantity(i+1) )
  {
    const char *label_key = "";
    switch( i )
    {
      case Activity:  label_key = "Activity";  break;
      case Distance:  label_key = "Distance";  break;
      case Dose:      label_key = "Dose";      break;
      case Shielding: label_key = "Shielding"; break;
      case NumQuantity:                    break;
    }//switch( i )
    
    WMenuItem *item = m_menu->addItem( WString::tr(label_key) );
    item->clicked().connect( boost::bind(&right_select_item, m_menu, item) );
    item->triggered().connect( boost::bind( &DoseCalcWidget::handleQuantityClick, this, i ) );

    m_enterWidgets[i] = new WContainerWidget();
    m_answerWidgets[i] = new WContainerWidget();
    
    enterLayout->addWidget( m_enterWidgets[i], i + 2, 0, AlignMiddle );
    answerLayout->addWidget( m_answerWidgets[i], i + 1, 0, AlignMiddle );

    
    m_enterWidgets[i]->addStyleClass( "DoseEnterEl" );
    m_answerWidgets[i]->addStyleClass( "DoseAnswerEl" );
    
    switch( i )
    {
      case Dose:
      {
        txt = new WText( WString("{1}:").arg(WString::tr("Dose")), m_enterWidgets[i] );
        txt->addStyleClass( "DoseEnterInd" );
        if( !isPhone )
          txt->setInline( false );
        else
          txt->setAttributeValue( "style", "display: inline-block; width: 45px;" );
        
        WContainerWidget *unitdiv = new WContainerWidget( m_answerWidgets[i] );
        WLabel *label = new WLabel( "units: ", unitdiv );
        m_doseAnswerUnits = new WComboBox( unitdiv );
        label->setBuddy( m_doseAnswerUnits );
        
        m_doseAnswer = new WText( m_answerWidgets[i] );
        m_doseAnswer->setInline( false );
        m_doseAnswer->addStyleClass( "DoseAnswerTxt" );
        
        m_doseAnswerValue = 0.1 * PhysicalUnits::rem / PhysicalUnits::hour;
        m_doseAnswer->setText( "100 mR/hr" );
        m_doseAnswerUnits->addItem( "rem/hr" );
        m_doseAnswerUnits->addItem( "sievert/hr" );

        
        m_doseEnter = new WLineEdit( m_enterWidgets[i] );
        
        m_doseEnter->setAutoComplete( false );
        m_doseEnter->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        m_doseEnter->setAttributeValue( "autocorrect", "off" );
        m_doseEnter->setAttributeValue( "spellcheck", "off" );
#endif
        m_doseEnter->addStyleClass( "DoseEnterTxt" );
        m_doseEnter->setText( "100" );
        if( isPhone )
          m_doseEnter->setWidth( 50 );
        
        WDoubleValidator *val = new WDoubleValidator( m_doseEnter );
        val->setBottom( 0.0 );
        m_doseEnter->setValidator( val );
        m_doseEnter->changed().connect( this, &DoseCalcWidget::updateResult );
//        m_doseEnter->blurred().connect( this, &DoseCalcWidget::updateResult );
//        m_doseEnter->enterPressed().connect( this, &DoseCalcWidget::updateResult );

        m_doseEnterUnits = new WComboBox( m_enterWidgets[i] );
        
        
        for( const auto &u : PhysicalUnits::sm_doseRateUnitHtmlNameValues )
        {
          string val = u.first;
          const unsigned char utf8mu[] = { 0xCE, 0xBC, 0 };
          SpecUtils::ireplace_all( val, "&mu;", (const char *)utf8mu /*"\u03BC"*/ );
          m_doseEnterUnits->addItem( WString::fromUTF8(val) );
        }
        
        m_doseEnterUnits->activated().connect( this, &DoseCalcWidget::updateResult );
        m_doseAnswerUnits->activated().connect( this, &DoseCalcWidget::updateResult );
        break;
      }//case Dose:
        
      case Activity:
      {
        txt = new WText( WString("{1}:").arg(WString::tr("Activity")) );
        txt->addStyleClass( "DoseEnterInd" );
        if( !isPhone )
          txt->setInline( false );
        else
          txt->setAttributeValue( "style", "display: inline-block; width: 45px;" );
        
        m_enterWidgets[i]->addWidget( txt );
        
        WContainerWidget *unitdiv = new WContainerWidget( m_answerWidgets[i] );
        WLabel *label = new WLabel( "units: ", unitdiv );
        m_activityAnswerUnits = new WComboBox( unitdiv );
        label->setBuddy( m_activityAnswerUnits );
        
        m_activityEnter = new WLineEdit( m_enterWidgets[i] );
        m_activityEnter->addStyleClass( "DoseEnterTxt" );
        
        m_activityEnter->setAutoComplete( false );
        m_activityEnter->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        m_activityEnter->setAttributeValue( "autocorrect", "off" );
        m_activityEnter->setAttributeValue( "spellcheck", "off" );
#endif
        if( isPhone )
          m_activityEnter->setWidth( 50 );
        
        WRegExpValidator *val = new WRegExpValidator( PhysicalUnits::sm_activityUnitOptionalRegex, this );
        val->setFlags( Wt::MatchCaseInsensitive );
//        WDoubleValidator *val = new WDoubleValidator();
        m_activityEnter->setValidator( val );
//        val->setBottom( 0.0 );
        m_activityEnter->setText( "100" );  
        m_activityEnter->blurred().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_activityEnter->enterPressed().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );

        m_activityEnterUnits = new WComboBox( m_enterWidgets[i] );

        for( const auto &u : PhysicalUnits::sm_activityUnitHtmlNameValues )
        {
          string val = u.first;
          const unsigned char utf8mu[] = { 0xCE, 0xBC, 0 };
          SpecUtils::ireplace_all( val, "&mu;", (const char *)utf8mu /*"\u03BC"*/ );
          m_activityEnterUnits->addItem( WString::fromUTF8(val) );
        }
        
        if( useBq )
          m_activityEnterUnits->setCurrentIndex( 2 ); //MBq
        else
          m_activityEnterUnits->setCurrentIndex( 6 ); //uci
        
        m_activityAnswerUnits->addItem( "curies" );
        m_activityAnswerUnits->addItem( "becquerels" );

        if( useBq )
          m_activityAnswerUnits->setCurrentIndex( 1 );
        else
          m_activityAnswerUnits->setCurrentIndex( 0 );
        
        m_activityAnswer = new WText( m_answerWidgets[i] );
        m_activityAnswer->setInline( false );

		    string actstr = useBq ? "1 MBq" : "200 &mu;Ci";
		    const unsigned char utf8mu[] = { 0xCE, 0xBC, 0 };
		    SpecUtils::ireplace_all(actstr, "&mu;", (char *)utf8mu /*"\u03BC"*/);

        m_activityAnswer->setText( WString::fromUTF8(actstr) );
        m_activityAnswer->addStyleClass( "DoseAnswerTxt" );
        
        m_activityEnterUnits->activated().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_activityAnswerUnits->activated().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        break;
      }//case Activity:
        
      case Distance:
      {
        txt = new WText( WString("{1}:").arg(WString::tr("Distance")) );
        txt->addStyleClass( "DoseEnterInd" );
        if( !isPhone )
          txt->setInline( false );
        else
          txt->setAttributeValue( "style", "display: inline-block; width: 45px;" );
        
        m_enterWidgets[i]->addWidget( txt );
        
        m_distanceAnswer = new WText();
        m_distanceAnswer->setInline( false );
        m_distanceAnswer->addStyleClass( "DoseAnswerTxt" );
        m_distanceAnswer->addStyleClass( "DoseDistanceAnswer" );
        m_answerWidgets[i]->addWidget( m_distanceAnswer );
    
        
        m_distanceEnter = new WLineEdit( "100 cm" );
        
        m_distanceEnter->setAutoComplete( false );
        m_distanceEnter->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
        m_distanceEnter->setAttributeValue( "autocorrect", "off" );
        m_distanceEnter->setAttributeValue( "spellcheck", "off" );
#endif
        m_distanceEnter->addStyleClass( "DoseEnterTxt" );
        
        if( isPhone )
          m_distanceEnter->setWidth( 50 );
        
//        m_distanceEnter->setTextSize( 10 );
        WRegExpValidator *validator = new WRegExpValidator( PhysicalUnits::sm_distanceUnitOptionalRegex, this );
        validator->setFlags( Wt::MatchCaseInsensitive );
        m_distanceEnter->setValidator( validator );
        HelpSystem::attachToolTipOn( m_distanceEnter, WString::tr("dcw-tt-distance"), showToolTips );
        
        m_distanceEnter->changed().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_distanceEnter->blurred().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_distanceEnter->enterPressed().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        
        m_enterWidgets[i]->addWidget( m_distanceEnter );
        break;
      }//case Distance:
        
      case Shielding:
      {
        if( !isPhone )
        {
          txt = new WText( WString("{1}:").arg(WString::tr("Shielding")) );
          txt->addStyleClass( "DoseEnterInd" );
          txt->setInline( false );
          m_enterWidgets[i]->addWidget( txt );
        }
        
        m_enterShieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest );
        m_answerShieldingSelect = new ShieldingSelect( m_materialDB, m_materialSuggest );
        m_answerShieldingSelect->setSphericalThicknessEditEnabled(false);
        m_answerShieldingSelect->arealDensityEdit()->disable();
        m_enterWidgets[i]->addWidget( m_enterShieldingSelect );
        m_answerWidgets[i]->addWidget( m_answerShieldingSelect );
        m_enterShieldingSelect->materialModified().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_enterShieldingSelect->materialChanged().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_answerShieldingSelect->materialModified().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        m_answerShieldingSelect->materialChanged().connect( boost::bind( &DoseCalcWidget::updateResult, this ) );
        
        if( isPhone )
        {
          m_enterShieldingSelect->materialEdit()->setEmptyText( WString::tr("dcw-shield-empty-text") );
        }
        
        break;
      }//case Shielding:
        
      case NumQuantity:
        break;
    }//switch( i )
  }//for( Quantity i = Quantity(0); i < NumQuantity; ++i )
  
//  buttonLayout->setContentsMargins( 0, 0, 0, 0 );
  //Add another widget below the menu buttons to force them all up
//  const int spacerRow = buttonLayout->rowCount();
//  buttonLayout->addWidget( new WContainerWidget(), spacerRow, 0 );
//  buttonLayout->setRowStretch( spacerRow, 10 );
  
  
  
  m_stayTime = new WContainerWidget();
  const int stayRow = answerLayout->rowCount();
  answerLayout->addWidget( m_stayTime, stayRow, 0 );
  
  m_issueTxt = new WText();
  m_issueTxt->addStyleClass( "DoseIssueTxt" );
  answerLayout->addWidget( m_issueTxt, answerLayout->rowCount(), 0, AlignBottom );
  
  if( !narrowLayout )
  {
    for( int i = 1; i < answerLayout->rowCount(); ++i )
      answerLayout->setRowStretch( i, 1 );
    answerLayout->setRowStretch( stayRow, 2 );
    
    for( int i = 1; i < enterLayout->rowCount(); ++i )
      enterLayout->setRowStretch( i, 1 );
  }//if( !narrowLayout )
  
  try
  {
    runtime_sanity_checks( m_scatter.get() );
  }catch( std::exception &e )
  {
    introDiv->clear();
    introDiv->setStyleClass( "DoseCalcRuntimeCheckFailMsg" );
    
    auto errorintro = new WText( WString::tr("dcw-runtime-checks-failed"), introDiv );
    errorintro->setInline( false );
    errorintro->setStyleClass( "DoseCalcRuntimeCheckFailHdr" );
    
    auto errortxt = new WText( e.what(), introDiv );
    errortxt->setInline( false );
    errortxt->setStyleClass( "DoseCalcRuntimeCheckMsgTxt" );
    
    auto furthermsg = new WText( WString::tr("dcw-runtime-check-msg"), introDiv );
    furthermsg->setStyleClass( "DoseCalcRuntimeCheckInstuct" );
    furthermsg->setInline( false );
    
    auto *errordiv = new WText( WString::tr("dcw-runtime-check-err-div") );
    errordiv->setStyleClass( "DoseCalcRuntimeCheckFailBanner" );
    const int row = m_layout->rowCount();
    m_layout->addWidget( errordiv, row, 0, 1, m_layout->columnCount(), AlignCenter );
  }//try / caltch runtime sanity check
  
}//void DoseCalcWidget::init()


DoseCalcWidget::~DoseCalcWidget()
{
  //nothing to do here
}//~DoseCalcWidget()


void DoseCalcWidget::handleAppUrl( std::string path, std::string query_str )
{
#if( PERFORM_DEVELOPER_CHECKS )
  const string expected_uri = path + "?" + query_str;
#endif
  
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  Quantity calcQuantity;
  if( SpecUtils::iequals_ascii(path, "dose") )
    calcQuantity = Quantity::Dose;
  else if( SpecUtils::iequals_ascii(path, "act") )
    calcQuantity = Quantity::Activity;
  else if( SpecUtils::iequals_ascii(path, "dist") )
    calcQuantity = Quantity::Distance;
  else if( SpecUtils::iequals_ascii(path, "shield") )
    calcQuantity = Quantity::Shielding;
  else if( SpecUtils::iequals_ascii(path, "intro") )
  {
    handleQuantityClick( Quantity::NumQuantity );
    return;
  }else
    throw runtime_error( "Dose Calc tool: invalid URI path." );
  
    
  string inshielduri, outshielduri;
  const size_t shieldpos = query_str.find( "&INSHIELD=&" );
  
  if( shieldpos != string::npos )
  {
    inshielduri = query_str.substr(shieldpos + 11);
    query_str = query_str.substr(0, shieldpos);
  }
  
  size_t outshieldpos = inshielduri.find("&OUTSHIELD=&");
  if( outshieldpos != string::npos )
  {
    outshielduri = inshielduri.substr(outshieldpos + 12);
    inshielduri = inshielduri.substr(0, outshieldpos);
  }else if( (outshieldpos = query_str.find("&OUTSHIELD=&")) != string::npos )
  {
    outshielduri = query_str.substr(outshieldpos + 12);
    query_str = query_str.substr(0, outshieldpos);
  }
  
  if( inshielduri.empty() )
    m_enterShieldingSelect->setToNoShielding();
  else
    m_enterShieldingSelect->handleAppUrl( inshielduri );
  
  if( outshielduri.empty() )
    m_answerShieldingSelect->setToNoShielding();
  else
    m_answerShieldingSelect->handleAppUrl( outshielduri );
  
  SpecUtils::ireplace_all( query_str, "%23", "#" );
  SpecUtils::ireplace_all( query_str, "%26", "&" );
  SpecUtils::ireplace_all( query_str, "curries", "curies" ); //fix up me being a bad speller
  
  const map<string,string> parts = AppUtils::query_str_key_values( query_str );
  const auto ver_iter = parts.find( "VER" );
  if( ver_iter == end(parts) )
    Wt::log("warn") << "No 'VER' field in Dose Calc tool URI.";
  else if( ver_iter->second != "1" && !SpecUtils::starts_with(ver_iter->second, "1.") )
    throw runtime_error( "Can not read Dose Calc tool URI version '" + ver_iter->second + "'" );
  
  auto findUnitIndex = [&parts]( const string &key, Wt::WComboBox *combo ) -> int {
    const auto act_unit_iter = parts.find(key);
    if( act_unit_iter == end(parts) )
      return -1;
    
    for( int i = 0; i < combo->count(); ++i )
    {
      if( SpecUtils::iequals_ascii( combo->itemText(i).toUTF8(), act_unit_iter->second ) )
        return i;
    }
    assert( 0 );
    return -1;
  };//findUnitIndex
  
  const int act_in_unit_index = findUnitIndex( "ACTINUNIT", m_activityEnterUnits );
  const int act_out_unit_index = findUnitIndex( "ACTOUTUNIT", m_activityAnswerUnits );
  const int dose_in_unit_index = findUnitIndex( "DOSEINUNIT", m_doseEnterUnits );
  const int dose_out_unit_index = findUnitIndex( "DOSEOUTUNIT", m_doseAnswerUnits );
  
  const auto act_iter = parts.find("ACT");
  if( act_iter != end(parts) && !act_iter->second.empty() && (act_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain activity unit info." );
  
  const auto dose_iter = parts.find("DOSE");
  if( dose_iter != end(parts) && !dose_iter->second.empty() && (dose_in_unit_index < 0) )
    throw runtime_error( "Dose Calc tool URI does not contain dose unit info." );
  
  m_menu->select( static_cast<int>(calcQuantity) );
  handleQuantityClick( calcQuantity );
  
  if( act_in_unit_index >= 0 )
    m_activityEnterUnits->setCurrentIndex( act_in_unit_index );
  if( act_out_unit_index >= 0 )
    m_activityAnswerUnits->setCurrentIndex( act_out_unit_index );
  if( dose_in_unit_index >= 0 )
    m_doseEnterUnits->setCurrentIndex( dose_in_unit_index );
  if( dose_out_unit_index >= 0 )
    m_doseAnswerUnits->setCurrentIndex( dose_out_unit_index );
  
  const auto nuc_iter = parts.find("NUC");
  if( nuc_iter != end(parts) )
    m_gammaSource->setNuclideText( nuc_iter->second );
  
  const auto age_iter = parts.find("AGE");
  if( age_iter != end(parts) )
    m_gammaSource->setNuclideAgeTxt( age_iter->second );
  
  if( act_iter != end(parts) )
    m_activityEnter->setText( WString::fromUTF8(act_iter->second) );
  else
    m_activityEnter->setText( "" );
  
  if( dose_iter != end(parts) )
    m_doseEnter->setText( WString::fromUTF8(dose_iter->second) );
  else
    m_doseEnter->setText( "" );
  
  const auto dist_iter = parts.find("DIST");
  if( dist_iter != end(parts) )
    m_distanceEnter->setText( WString::fromUTF8(dist_iter->second) );
  else
    m_distanceEnter->setText( "" );
  
  updateResult();
  
#if( PERFORM_DEVELOPER_CHECKS )
  if( m_stateUri != expected_uri )
  {
    Wt::log("warn") << "DoseCalcWidget::handleAppUrl: input URI doesnt match current URI.\n\t input: '"
                    << expected_uri.c_str() << "'\n\tresult: '" << m_stateUri.c_str() << "'";
  }
#endif
}//void handleAppUrl( std::string uri )


std::string DoseCalcWidget::encodeStateToUrl() const
{
  // "interspec://dose/act?nuc=u238&dose=1.1ur/h&dist=100cm&..."
  
  string answer;
  
  switch( m_currentCalcQuantity )
  {
    case Dose:        answer += "dose";   break;
    case Activity:    answer += "act";    break;
    case Distance:    answer += "dist";   break;
    case Shielding:   answer += "shield"; break;
    case NumQuantity: answer += "intro"; break;
  }//switch( m_currentCalcQuantity )
  
  answer += "?VER=1";

  if( m_currentCalcQuantity == NumQuantity )
    return answer;
  
  // We could limit what info we put in the URL, based on current m_currentCalcQuantity,
  //  but we might as well put the full state into the URI.
  
  auto addTxtField = [&answer]( const string &key, Wt::WLineEdit *edit ){
    string txt = edit->text().toUTF8();
    SpecUtils::ireplace_all( txt, "#", "%23" );
    SpecUtils::ireplace_all( txt, "&", "%26" );
    answer += "&" + key + "=" + txt;
  };

  const SandiaDecay::Nuclide *nuc = m_gammaSource->nuclide();
  if( nuc )
    answer += "&NUC=" + nuc->symbol;
  
  if( nuc && !PeakDef::ageFitNotAllowed(nuc) )
    answer += "&AGE=" + m_gammaSource->nuclideAgeStr().toUTF8();
  
  addTxtField( "ACT", m_activityEnter );
  answer += "&ACTINUNIT=" + m_activityEnterUnits->currentText().toUTF8();
  answer += "&ACTOUTUNIT=" + m_activityAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DOSE", m_doseEnter );
  answer += "&DOSEINUNIT=" + m_doseEnterUnits->currentText().toUTF8();
  answer += "&DOSEOUTUNIT=" + m_doseAnswerUnits->currentText().toUTF8();
  
  addTxtField( "DIST", m_distanceEnter );
  
  // We'll mark shielding URL starting with the below, and everything after this is the shielding.
  //  Having an order-independent method would be better, but for the moment...
  switch( m_currentCalcQuantity )
  {
    case Dose:
    case Activity:
    case Distance:
      if( !m_enterShieldingSelect->isNoShielding() )
        answer += "&INSHIELD=&" + m_enterShieldingSelect->encodeStateToUrl();
      break;
      
    case Shielding:
      if( !m_answerShieldingSelect->isNoShielding() )
        answer += "&OUTSHIELD=&" + m_answerShieldingSelect->encodeStateToUrl();
      break;
      
    case NumQuantity:
      break;
  }//switch( m_currentCalcQuantity )
  
  
  
  return answer;
}//std::string encodeStateToUrl() const;


void DoseCalcWidget::runtime_sanity_checks( const GadrasScatterTable * const scatter )
{
  if( !scatter )
    throw runtime_error( "Full spectrum transport source matrix not initiated." );
  

  auto check_nuc = [scatter]( const string nuclabel, const double age, const float distance,
                        const float areal_density, const float atomic_number, const double expected ) {
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    if( !db )
      throw runtime_error( "Nuclide database not initiated" );
    
    auto nuc = db->nuclide( nuclabel );
    if( !nuc )
      throw runtime_error( "Error retrieving '" + nuclabel + "' from nuclide database." );
    
    SandiaDecay::NuclideMixture mix;
    mix.addAgedNuclideByActivity( nuc, 100.0E-6*SandiaDecay::curie, age );
    
    vector<float> energies, intensities;
    for( const auto &i : mix.photons(0) )
    {
      energies.push_back( i.energy );
      intensities.push_back( i.numPerSecond );
    }

    const double computed = DoseCalc::gamma_dose_with_shielding( energies, intensities,
                                            areal_density, atomic_number,
                                                  distance, *scatter );
    
    //Check within 1% of expected (1% is arbitrary)
    if( fabs(computed - expected) > 0.02*std::max(computed,expected) )
      throw runtime_error( "Failed " + nuclabel + " with shielding {AN="
                           + std::to_string(atomic_number) + ", AD="
                           + std::to_string(areal_density * PhysicalUnits::cm2 / PhysicalUnits::gram)
                           + "g/cm2}, age=" + PhysicalUnits::printToBestTimeUnits(age)
                           + " and distance=" + PhysicalUnits::printToBestLengthUnits(distance)
                           + ", the sanity check computed "
                           + PhysicalUnits::printToBestEquivalentDoseRateUnits(computed,4,false)
                           + " but expected "
                           + PhysicalUnits::printToBestEquivalentDoseRateUnits(expected,4,false)
                           + ".  Not performing further tests.");
  };//check_nuc lambda
  
  string nuclide = "Co60";
  double age = 0.5*PhysicalUnits::year;
  float distance = 100.0*PhysicalUnits::cm;
  float areal_density = 0.0f;
  float atomic_number = 26.0f;
  double expected = 115.42E-6 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 115.4
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  nuclide = "Cs137";
  age = 0.5*PhysicalUnits::year;
  distance = 100.0*PhysicalUnits::cm;
  areal_density = 0.0f * PhysicalUnits::gram / PhysicalUnits::cm2;
  atomic_number = 26.0f;
  expected = 29.70E-6 * PhysicalUnits::rem/PhysicalUnits::hour;  //Other program gives 29.7 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  areal_density = 5.0f * static_cast<float>(PhysicalUnits::gram / PhysicalUnits::cm2);
  expected = 25.19E-6 * PhysicalUnits::rem/PhysicalUnits::hour;  //Other program gives 23.9 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  areal_density = 50.0f * static_cast<float>(PhysicalUnits::gram / PhysicalUnits::cm2);
  expected = 2.59E-6 * PhysicalUnits::rem/PhysicalUnits::hour;  //Other program gives 0.842 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  nuclide = "U238";
  age = 20*PhysicalUnits::year;
  distance = 100.0*PhysicalUnits::cm;
  areal_density = 0.0f * PhysicalUnits::gram / PhysicalUnits::cm2;
  atomic_number = 60.0f;
  expected = 1.86E-6 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 1.8 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  areal_density = 2.0f * static_cast<float>(PhysicalUnits::gram / PhysicalUnits::cm2);
  expected = 2.34E-6 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 0.78 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  nuclide = "Na22";
  age = 0.0*PhysicalUnits::year;
  distance = 10.0 * static_cast<float>(PhysicalUnits::cm);
  areal_density = 13.0f * static_cast<float>( PhysicalUnits::gram / PhysicalUnits::cm2 );
  atomic_number = 5.0f;
  expected = 7.66E-3 * PhysicalUnits::rem/PhysicalUnits::hour;  //Other program gives 6.1 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  nuclide = "F18";
  age = 0.0*PhysicalUnits::year;
  distance = 200.0*PhysicalUnits::cm;
  areal_density = 0.0f * PhysicalUnits::gram / PhysicalUnits::cm2;
  atomic_number = 5.0f;
  expected = 13.29E-6 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 13.3 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  
  nuclide = "Ba133";
  age = 7*24*3600*PhysicalUnits::second;
  distance = 10.0*PhysicalUnits::cm;
  areal_density = 0.0f * PhysicalUnits::gram / PhysicalUnits::cm2;
  atomic_number = 82.0f;
  expected = 2.44E-3 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 2.4 mRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  areal_density = 10.0f * static_cast<float>(PhysicalUnits::gram / PhysicalUnits::cm2);
  expected = 134.09E-6 * PhysicalUnits::rem/PhysicalUnits::hour; //Other program gives 111 uRem/h
  check_nuc( nuclide, age, distance, areal_density, atomic_number, expected );
  
  
  auto checkPrintDose = [=]( double dose, const string answer, const bool useSv ){
    dose *= useSv ? PhysicalUnits::sievert : PhysicalUnits::rem;
    dose /= PhysicalUnits::hour;
    
    const string value = PhysicalUnits::printToBestEquivalentDoseRateUnits( dose, 2, useSv );
  
    if( value != answer )
      throw runtime_error( "Error printing to dose string; expected '"
                           + answer + "', but got '" + value + "'" );
  };//checkPrintDose(...)
  
  // The printToBest*(...) functions currently internally use snprintf, with a format flag like
  //  '%.2f', which behaves a bit different on the different platforms for trailing zeros, so will
  //  avoid that case for this test - I'm also not sure if rounding is handled different on Windows
  //  e.g., if “round to nearest and ties to even” or "nearest-even"
  checkPrintDose( 1.234, "1.23 sv/hr", true );
  checkPrintDose( 1.234, "1.23 rem/hr", false );
  checkPrintDose( 8.236, "8.24 sv/hr", true );
  checkPrintDose( 8.23688E-7, "823.69 nsv/hr", true );
  checkPrintDose( 8.23688E-5, "82.37 usv/hr", true );
  checkPrintDose( 8.23688E-2, "82.37 msv/hr", true );
  checkPrintDose( 8.23688E2, "823.69 sv/hr", true );
  checkPrintDose( 8.23688E5, "823.69 ksv/hr", true );
  checkPrintDose( 8.23688E6, "8.24 Msv/hr", true );
}//void runtime_sanity_checks()


void DoseCalcWidget::handleQuantityClick( const DoseCalcWidget::Quantity q )
{
  //const bool isDeselect = (q == m_currentCalcQuantity);
  //if( isDeselect )
  //{
  //  m_menu->select( nullptr );
  //  m_stack->setCurrentIndex( 0 );
  //  m_currentCalcQuantity = NumQuantity;
  //  return;
  //}//if( isDeselect )
  
  m_currentCalcQuantity = q;
  
  m_stack->setCurrentIndex( (q == NumQuantity) ? 0 : 1 );
  
  for( Quantity i = Quantity(0); i < NumQuantity; i = Quantity(i+1) )
  {
    m_enterWidgets[i]->setHidden( i==q );
    m_answerWidgets[i]->setHidden( i!=q );
  }
  
  updateResult();
}//void handleQuantityClick( const Quantity q )


void DoseCalcWidget::handleSourceTypeChange()
{
  const int checkid = m_sourceType->checkedId();
  
  m_gammaSourceDiv->setHidden( checkid==1 );
  m_neutronSourceDiv->setHidden( checkid==0 );
  
  if( checkid == 0 )
  {
   
  }else
  {
  
  }//if( gamma source ) / else ( neutron source )
  
  updateResult();
}//void DoseCalcWidget::handleSourceTypeChange()


double DoseCalcWidget::enteredActivity()
{
  double activity = 0.0;
  string activitystr = m_activityEnter->text().toUTF8();
  
  SpecUtils::trim( activitystr );
  
  const size_t pos = activitystr.find_first_not_of( "0123456789Ee.-+" );  //a regex would be better
  
  if( pos == string::npos )
  {
    if( !(stringstream(activitystr) >> activity) )
      throw runtime_error( "Invalid number" );
    const int unitsInd = m_activityEnterUnits->currentIndex();
    activity *= PhysicalUnits::sm_activityUnitHtmlNameValues.at(unitsInd).second;
  }else
  {
    const bool hasb = SpecUtils::icontains( activitystr, "b" );
    const bool hasc = SpecUtils::icontains( activitystr, "c" );
    
    if( hasb && hasc )
      throw runtime_error( "Invalid activity string, couldnt determine units" );
    
    activity = PhysicalUnits::stringToActivity( activitystr );
    
    const PhysicalUnits::UnitNameValuePair &bestunit
                          = PhysicalUnits::bestActivityUnitHtml( activity, hasc );
    
    PhysicalUnits::UnitNameValuePairV::const_iterator pos =
         std::find( PhysicalUnits::sm_activityUnitHtmlNameValues.begin(), PhysicalUnits::sm_activityUnitHtmlNameValues.end(), bestunit );
    assert( pos != PhysicalUnits::sm_activityUnitHtmlNameValues.end() );
    const int activity_index = static_cast<int>( pos - PhysicalUnits::sm_activityUnitHtmlNameValues.begin() );
    m_activityEnterUnits->setCurrentIndex( activity_index );
    
    const double txtdblval = activity / bestunit.second;
    char txtval[32];
    snprintf( txtval, sizeof(txtval), "%f", floor(1000*txtdblval + 0.5)/1000 );  //fix this to 3 sig figs, not 3 decimals
    
    //Get rid of unwanted decimals (they should be all zeros from rounding above)
    //(this should probably be a function call somewhere)
    string durstr = txtval;
    const size_t decpos = durstr.find_last_of( '.' );
    if( decpos != string::npos )
    {
      const string predec = durstr.substr( 0, decpos );
      string postdec = durstr.substr( decpos + 1 );
      size_t lastpos = postdec.find_last_not_of( "0" );
      if( lastpos == string::npos )
        postdec = "";
      else
        postdec = postdec.substr( 0, lastpos + 1 );
      
      durstr = predec;
      if( postdec.size() )
        durstr += "." + postdec;
    }//if( decpos != string::npos )
    
    m_activityEnter->setText( durstr );
  }//if( user entered only numbers ) / else user entered units to
  
  return activity;
}//double enteredActivity()





void DoseCalcWidget::updateResultForGammaSource()
{
  double activity = 0.0;
  
  Quantity quantity = m_currentCalcQuantity;
  
  if( !m_gammaSource->nuclide() )
  {
    m_issueTxt->setText( WString::tr("dcw-err-no-nuc") );
    quantity = NumQuantity;
  }
  
  m_activityAnswer->setText( "" );
  m_doseAnswerValue = 0.0;
  m_doseAnswer->setText( "" );
  m_distanceAnswer->setText( "" );
  if( m_answerShieldingSelect->isGenericMaterial() )
    m_answerShieldingSelect->arealDensityEdit()->setText( "" );
  else
    m_answerShieldingSelect->setSphericalThickness( -1 );
  m_issueTxt->setText( "" );
  

  try
  {
    bool doflash = false;
    switch( quantity )
    {
      case Activity: case Dose: case Distance:
        if( m_enterShieldingSelect->isGenericMaterial() )
        {
          doflash = (m_enterShieldingSelect->atomicNumber() > 0.0
                      && m_enterShieldingSelect->arealDensity() > 0.0);
        }else
        {
          doflash = (m_enterShieldingSelect->material()
                     && m_enterShieldingSelect->thickness() > 0.0);
        }
        break;
      
      case Shielding: case NumQuantity:
      break;
    }
    
    if( doflash )
      doJavaScript( "$('#" + m_enterShieldingSelect->id() + "').fadeIn(100).fadeOut(100).fadeIn(100).fadeOut(100).fadeIn(100);" );
  }catch(...)
  {
    
  }
  
  switch( quantity )
  {
    case Activity:
      activity = 1.0*SandiaDecay::curie;
      break;
      
    case Dose:
    case Distance:
    case Shielding:
    {
      try
      {
        activity = enteredActivity();
      }catch( std::exception & )
      {
        m_issueTxt->setText( WString::tr("dcw-act-invalid") );
        return;
      }
      
      break;
    }//case Dose: case Distance: case Shielding:
      
    case NumQuantity:
      return;
      break;
  }//switch( m_sourceType->checkedId() )
  
  const vector<pair<float,float> > source_gamma
                           = m_gammaSource->photonEnergyAndIntensity( activity );

  if( source_gamma.empty() )
  {
    m_issueTxt->setText( "Error determining source gammas" );
    return;
  }
  
  //XXX - should change this from {AN,AD} to getting mu for the material specifically.
  float shielding_an = 14.0f;
  float shielding_ad = 0.0f;
  
  try
  {
    if( m_enterShieldingSelect->isGenericMaterial() )
    {
      string anstr;
      if( !m_enterShieldingSelect->atomicNumberEdit()->text().empty()
         && !m_enterShieldingSelect->arealDensityEdit()->text().empty() )
      {
        shielding_an = m_enterShieldingSelect->atomicNumber();
        shielding_ad = m_enterShieldingSelect->arealDensity();
      }
    }else
    {
      std::shared_ptr<const Material> mat = m_enterShieldingSelect->material();
      if( mat )
      {
        shielding_an = mat->massWeightedAtomicNumber();
        shielding_ad = mat->density * m_enterShieldingSelect->thickness();
      }//if( mat )
    }//if( m_enterShieldingSelect->isGenericMaterial() ) / else
  }catch(...)
  {
    m_issueTxt->setText( WString::tr("dcw-shield-invalid") );
    return;
  }
  
  vector<float> energies, intensities;
  for( size_t i = 0; i < source_gamma.size(); ++i )
  {
    energies.push_back( source_gamma[i].first );
    intensities.push_back( source_gamma[i].second );
  }
  
  //Set to zero shielding if we are goign to calculate shielding
  //  Or if the user has explicitly requested zero shielding
  if( (quantity == Shielding) || (shielding_an <= 1.0 && shielding_ad <= 0.0) )
  {
    shielding_an = 14.0;
    shielding_ad = 0.0;
  }else if( shielding_an <= 1.0 )
  {
    m_issueTxt->setText( WString::tr("dcw-generic-shield-invalid") );
    return;
  }
  
  double dose_from_source = 0.0, user_entered_dose = 0.0, distance = 10.0*PhysicalUnits::cm;
  
  if( quantity != Distance )
  {
    try
    {
      string diststr = m_distanceEnter->text().toUTF8();
      SpecUtils::trim( diststr );
      if( diststr.find_first_not_of( " \t0123456789.eE+-\n" ) == string::npos )
      {
        diststr += " cm";
        m_distanceEnter->setText( diststr );
      }
      
      distance = PhysicalUnits::stringToDistance( diststr );
    }catch ( std::exception &e )
    {
      m_issueTxt->setText( WString::tr("dcw-dist-invalid") );
      cerr << "Error: " << e.what() << endl;
      return;
    }
  }//if( quantity != Distance )
  
  
  try
  {
    dose_from_source = DoseCalc::gamma_dose_with_shielding( energies, intensities,
                                        shielding_ad, shielding_an, distance, *m_scatter );
  }catch( std::exception &e )
  {
    m_issueTxt->setText( WString::tr("dcw-err-other") );
    cerr << "Error: " << e.what() << endl;
    return;
  }//try / catch
  
  if( quantity != Dose )
  {
    if( !(stringstream(m_doseEnter->text().toUTF8()) >> user_entered_dose) )
    {
      m_issueTxt->setText( WString::tr("dcw-invalid-src-dose") );
      return;
    }
    
    const int doseIndex = m_doseEnterUnits->currentIndex();
    if( doseIndex >= static_cast<int>(PhysicalUnits::sm_doseRateUnitHtmlNameValues.size()) || doseIndex < 0 )
      throw runtime_error( "doseIndex >= PhysicalUnits::sm_doseRateUnitHtmlNameValues.size() || doseIndex < 0" );
    
    user_entered_dose *= PhysicalUnits::sm_doseRateUnitHtmlNameValues[doseIndex].second;
  }//if( quantity != Dose )
  
  
  switch( quantity )
  {
    case Activity:
    {
      const double answer = activity * user_entered_dose / dose_from_source;
      
      
      const bool useCuries = (m_activityAnswerUnits->currentIndex() == 0);
      const string answerstr = PhysicalUnits::printToBestActivityUnits( answer, 2, useCuries );
      m_activityAnswer->setText( answerstr );
      
      break;
    }//case Activity:
      
    case Dose:
    {
      try
      {
        enteredActivity();
      }catch(...)
      {
        m_issueTxt->setText( WString::tr("dcw-invalid-src-act") );
        return;
      }
      
      const bool useRemHr = (m_doseAnswerUnits->currentIndex() == 0);
      const string answerstr = PhysicalUnits::printToBestEquivalentDoseRateUnits( dose_from_source, 2, !useRemHr );
      m_doseAnswerValue = dose_from_source;
      m_doseAnswer->setText( answerstr );

      break;
    }//case Dose:
      
    case Distance:
    {
      if( user_entered_dose <= 0.0 || IsInf(user_entered_dose) || IsNan(user_entered_dose) )
      {
        m_issueTxt->setText( WString::tr("dcw-invalid-src-dose") );
        return;
      }
      
      // Note: we cant actually use 'dist_guess' as is; it doesnt take into account attenuation in
      //       the air, so we'll first take a guess (based on the default distance of 10 cm) that
      //       doesnt account for air attenuation then use this to find a set of distances that
      //       bracket the actual distance that yields the desired dose.
      const double ratio = dose_from_source / user_entered_dose;
      const double dist_guess = distance * sqrt(ratio);
      
      try
      {
        // Choose the accuracy we want to find the answer to; we could probably loosen this up a bit
        const double dist_delta = ((dist_guess < PhysicalUnits::meter) ? 0.1 : 1.0) * PhysicalUnits::mm;
        
        // Make a lambda to easily compute dose for a trial distance
        auto dose_root_calc = [=]( const double radius ) -> double {
          const double dose = DoseCalc::gamma_dose_with_shielding( energies, intensities,
                                                  shielding_ad, shielding_an, radius, *m_scatter );
          return user_entered_dose - dose;
        };//dose_root_calc lamda
        
        auto term_cond = [=]( double x0, double x1 ) -> bool { return fabs(x0 - x1) < dist_delta; };
        
        // Makeup some termination conditions so we dont get stuck in infinite loops
        boost::uintmax_t num_iters;
        const boost::uintmax_t max_bounds_iters = 25;  //seems to take 2 iterations or fewer
        const boost::uintmax_t max_bisect_iters = 250; //seems to take around 20 iterations
        
        // Find lower bounding distance
        double lower_dist = dist_guess, upper_dist = dist_guess;
        for( num_iters = 0; num_iters < max_bounds_iters; ++num_iters )
        {
          if( dose_root_calc(lower_dist) < 0.0 )
            break;
          lower_dist *= 0.5;
        }//for( num_iters = 0; num_iters < max_bounds_iters; ++num_iters )
        
        if( num_iters >= max_bounds_iters )
          throw runtime_error( WString::tr("dcw-failed-find-lower-dist")
                              .arg(lower_dist/PhysicalUnits::cm).toUTF8() );
        
        // Find upper bounding distance
        for( num_iters = 0; num_iters < max_bounds_iters; ++num_iters )
        {
          if( dose_root_calc(upper_dist) > 0.0 )
            break;
          upper_dist *= 2.0;
        }//for( num_iters = 0; num_iters < max_bounds_iters; ++num_iters )
        
        if( num_iters >= max_bounds_iters )
          throw runtime_error( WString::tr("dcw-failed-find-upper-dist")
                              .arg(upper_dist/PhysicalUnits::cm).toUTF8() );
        
        // Now find the actual distance to within 'dist_delta'
        num_iters = max_bisect_iters;
        const pair<double, double> bracketing_dists
           = boost::math::tools::bisect( dose_root_calc, lower_dist, upper_dist, term_cond, num_iters );
        
        if( num_iters >= max_bisect_iters )
          throw runtime_error( WString::tr("dcw-failed-find-dist").toUTF8() );
        
        assert( term_cond(bracketing_dists.first,bracketing_dists.second) );
        
        distance = 0.5*(bracketing_dists.first + bracketing_dists.second);
        
        if( distance <= 0.0 )
          throw runtime_error( WString::tr("dcw-got-neg-dist").arg(distance).toUTF8() );
        
        //const double final_dose = DoseCalc::gamma_dose_with_shielding( energies, intensities,
        //                                        shielding_ad, shielding_an, distance, *m_scatter );
        //cout << "Took " << num_iters << " to get distance in range ["
        //     << PhysicalUnits::printToBestLengthUnits(bracketing_dists.first)
        //     << "," << PhysicalUnits::printToBestLengthUnits(bracketing_dists.second) << "]"
        //     << " for a dose of " << PhysicalUnits::printToBestEquivalentDoseRateUnits(final_dose)
        //     << endl;
        
        assert( distance > 0.0 );
        m_distanceAnswer->setText( PhysicalUnits::printToBestLengthUnits(distance) );
      }catch( std::exception &e )
      {
        m_issueTxt->setText( WString::tr("dcw-err-calc-dose-for-dist") );
        cerr << "Error calculating dose while searching for distance: " << e.what() << endl;
        return;
      }//try / catch
      
      break;
    }//case Distance:
      
    case Shielding:
    {
      try
      {
        std::shared_ptr<const Material> mat;
        const bool isGeneric = m_answerShieldingSelect->isGenericMaterial();
      
        if( isGeneric )
        {
          shielding_an = m_answerShieldingSelect->atomicNumber();
          m_answerShieldingSelect->arealDensityEdit()->setText( "" );
        }else
        {
          m_answerShieldingSelect->setSphericalThickness( -1 );
        
          mat = m_answerShieldingSelect->material();
          if( !mat )
            throw runtime_error( WString::tr("dcw-material-not-valid").toUTF8() );
          
          shielding_an = mat->massWeightedAtomicNumber();
        }//if( isGeneric ) / else
      
  
        if( shielding_an < 1.0 || shielding_an > 98.0 )
        {
          if( isGeneric )
            m_issueTxt->setText( WString::tr("dcw-no-an-entered") );
          else
            m_issueTxt->setText( WString::tr("dcw-no-shield-name") );
          return;
        }//if( shielding_an < 1.0 )
      
      
        const double adfit = fit_ad( energies, intensities, shielding_an, user_entered_dose, distance, *m_scatter );
        
        if( isGeneric )
        {
          m_answerShieldingSelect->arealDensityEdit()->setValue( static_cast<float>(adfit*PhysicalUnits::cm2/PhysicalUnits::g) );
        }else
        {
          const double thickness = adfit / mat->density;
          m_answerShieldingSelect->setSphericalThickness( thickness );
        }
        
      }catch( std::exception &e )
      {
        m_issueTxt->setText( e.what() );
        return;
      }//try / catch
  
      break;
    }//case Shielding:
      
    case NumQuantity:
      return;
      break;
  }//switch( m_sourceType->checkedId() )
  
  checkAndWarnForBrehmSource();
  
  try
  {
    string uri = encodeStateToUrl();
    const bool sameAsPrev = (uri == m_stateUri);
    UndoRedoManager *undoRedo = UndoRedoManager::instance();
    if( undoRedo && undoRedo->canAddUndoRedoNow() && !sameAsPrev )
    {
      const shared_ptr<const string> prev = make_shared<string>( m_stateUri );
      const shared_ptr<const string> current = make_shared<string>( uri );
      
      auto undo_redo = [prev, current]( bool is_undo ){
        DoseCalcWindow *dosewin = InterSpec::instance()->showDoseTool();
        DoseCalcWidget *tool = dosewin ? dosewin->tool() : nullptr;
        const string &uri = is_undo ? *prev : *current;
        if( tool && !uri.empty() )
        {
          const size_t pos = uri.find('?');
          string path = (pos == string::npos) ? string() : uri.substr(0,pos);
          string query = (pos == string::npos) ? uri : uri.substr(pos+1);
          tool->handleAppUrl(path, query);
        }//
      };//undo_redo
      
      auto undo = [undo_redo](){ undo_redo(true); };
      auto redo = [undo_redo](){ undo_redo(false); };
      
      undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Update Dose Calc values." );
    }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
    
    if( !sameAsPrev )
      m_stateUri = std::move(uri);
  }catch( std::exception &e )
  {
    Wt::log("error") << "DoseCalcWidget::updateResultForGammaSource:"
                        " error trying to make undo/redo step: " << e.what();
  }//try / catch
}//void updateResultForGammaSource()


void DoseCalcWidget::checkAndWarnForBrehmSource()
{
  switch( m_currentCalcQuantity )
  {
    case Activity: case Dose:
    case Distance: case Shielding:
      break;
      
    case NumQuantity:
      return;
      break;
  }//switch( m_currentCalcQuantity )
  
  const SandiaDecay::Nuclide *nuc = m_gammaSource->nuclide();
  if( !nuc )
    return;
  
  double age = 0.0;
  try{ age = m_gammaSource->nuclideAge();}catch(...){ return; }
  
  SandiaDecay::NuclideMixture mix;
  mix.addAgedNuclideByActivity( nuc, 0.1*PhysicalUnits::curie, age );
  const vector<SandiaDecay::EnergyRatePair> betas = mix.betas( 0.0 );
  const vector<SandiaDecay::EnergyRatePair> betaplus = mix.betaPlusses( 0.0 );
  const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas( 0.0, SandiaDecay::NuclideMixture::OrderByEnergy, true );
  
  const double betaMaxFrac = 5.0;                        //selected arbitrarily
  const double minBetaEnergy = 500.0*PhysicalUnits::keV;  //selected arbitrarily
  const double minGammaEnergy = 40.0*PhysicalUnits::keV;  //selected arbitrarily
  
  double betaIntensity = 0.0, gammaIntensity = 0.0;
  for( const SandiaDecay::EnergyRatePair &aep : betas )
  {
    if( aep.energy > minBetaEnergy )
      betaIntensity += aep.numPerSecond;
  }
  for( const SandiaDecay::EnergyRatePair &aep : betaplus )
  {
    if( aep.energy > minBetaEnergy )
      betaIntensity += aep.numPerSecond;
  }
  for( const SandiaDecay::EnergyRatePair &aep : gammas )
  {
    if( aep.energy > minGammaEnergy )
      gammaIntensity += aep.numPerSecond;
  }
  
  //P32, S35, Sr90, and Y90 are the main isotopes we want to catch

  if( betas.empty() )
    return;
  
  if( (betaMaxFrac*gammaIntensity > betaIntensity) && !gammas.empty() )
    return;
  
  WString issuetxt = m_issueTxt->text();
  if( issuetxt.empty() )
  {
    issuetxt = WString("<div>{1}</div>")
      .arg( WString::tr("dcw-brehm-warning") );
  }else
  {
    issuetxt = WString("<div>{1}</div><div>{2}</div>")
      .arg(issuetxt)
      .arg( WString::tr("dcw-brehm-warning") );
  }
  
  m_issueTxt->setText( issuetxt );
}//void checkAndWarnForBrehmSource()


void DoseCalcWidget::updateResultForNeutronSource()
{
  cerr << "Need to implement updateResultForNeutronSource()" << endl;
}//void updateResultForNeutronSource()


void DoseCalcWidget::updateResult()
{
  if( m_sourceType->id( m_sourceType->selectedButton() ) == 0 )
    updateResultForGammaSource();
  else
    updateResultForNeutronSource();
  
  updateStayTime();
}//void updateResult()


double DoseCalcWidget::currentDose()
{
  double dose = 0;
  
  switch( m_currentCalcQuantity )
  {
    case Dose:
      if( m_doseAnswer->text().empty() || (m_doseAnswerValue < FLT_MIN) )
        throw runtime_error( WString::tr("dcw-no-dose-entered").toUTF8() );
      dose = m_doseAnswerValue;
      break;
      
    case Activity: case Distance: case Shielding:
    { 
      if( !(stringstream(m_doseEnter->text().toUTF8()) >> dose) )
        throw runtime_error( WString::tr("dcw-no-dose-entered").toUTF8() );
      
      const int doseIndex = m_doseEnterUnits->currentIndex();
      if( doseIndex >= static_cast<int>(PhysicalUnits::sm_doseRateUnitHtmlNameValues.size()) || doseIndex < 0 )
        throw runtime_error( "doseIndex >= PhysicalUnits::sm_doseRateUnitHtmlNameValues.size() || doseIndex < 0" );
      
      dose *= PhysicalUnits::sm_doseRateUnitHtmlNameValues[doseIndex].second;
      
      break;
    }//case Activity: case Distance: case Shielding:
      
    case NumQuantity:
      throw runtime_error( "No dose quantity" );
  }//switch( m_currentCalcQuantity )
  
  return dose;
}//double currentDose()


void DoseCalcWidget::updateStayTime()
{
  double dose = 0;
  
  try
  {
    dose = currentDose();
    
    m_stayTime->clear();
    if( !m_stayTime->hasStyleClass( "DoseStayTime" ) )
      m_stayTime->addStyleClass( "DoseStayTime" );
    
    WText *header = new WText( WString::tr("dcw-stay-time"), m_stayTime );
    header->setInline( false );
    header->addStyleClass( "DoseStayTimeHeader" );
    
    WContainerWidget *tableDiv = new WContainerWidget( m_stayTime );
    tableDiv->addStyleClass( "DoseStayTimeTable" );
    WGridLayout *table = new WGridLayout( tableDiv );
    table->setHorizontalSpacing( 0 );
    table->setVerticalSpacing( 0 );
    table->setContentsMargins( 0, 0, 0, 0 );
    
    map<double,pair<string,WString> > ref_doses;
    
    using PhysicalUnits::rem;
    
    //https://en.wikipedia.org/wiki/Banana_equivalent_dose
    //ref_doses[0.000010*rem]  = make_pair( "10 " MU_CHARACTER "rem", "Eating one banana" );
    
    //http://www.nrc.gov/about-nrc/radiation/around-us/doses-daily-lives.html
    ref_doses[0.0015*rem] = make_pair( "1.5 mrem", WString::tr("dcw-dental-xray") );
    //Flight from NY to LA 2 to 5 mrem
    //Chest xray 10 mrem
    ref_doses[0.620*rem]  = make_pair( "620 mrem", WString::tr("dcw-annual-background") );
  
    //National Council on Radiation Protection (NCRP) dose guidelines for civilian radiation workers  (wcjohns did not verify)
    ref_doses[5.0*rem]    = make_pair( "5 rem",    WString::tr("dcw-annual-occ-limit") );
    ref_doses[10.0*rem]   = make_pair( "10 rem",   WString::tr("dcw-save-property") );
    ref_doses[25.0*rem]   = make_pair( "25 rem",   WString::tr("dcw-life-saving") );
    ref_doses[100.0*rem]  = make_pair( "100 rem",  WString::tr("dcw-radiation-sickness") );
    
    
    //https://en.wikipedia.org/wiki/Sievert#Dose_examples_2
    
    //1 Sv = 100 rem
    //1 uSv = 0.1 mrem = 100 urem.
    
    //ref_doses[0.098*PhysicalUnits::sievert*1.0E-6] = make_pair( "9.8 urem", "banana equivalent dose" );  //0.098  μSv
    //0.25  μSv (25 urem):  U.S. limit on effective dose from a single airport security screening[34]
    //5 to 10  μSv (0.5 to 1 mrem):  one set of dental radiographs[35]
    //1  mSv (100 mrem):  The U.S. 10 CFR § 20.1301(a)(1) dose limit for individual members of the public, total effective dose equivalent, per annum[38]
    //1.5 to 1.7  mSv (150 to 170 mrem):  annual dose for flight attendants[39]
    //10 to 30  mSv (1 to 3 rem):  single full-body CT scan[41][42]
    //50  mSv (5 rem):  The U.S. 10 C.F.R. § 20.1201(a)(1)(i) occupational dose limit, total effective dose equivalent, per annum[43]
    //80  mSv (8 rem):  6-month stay on the International Space Station
    //250  mSv (25 rem):  6-month trip to Mars - radiation due to cosmic rays, which are very difficult to shield against
    //500  mSv (50 rem):  The U.S. 10 C.F.R. § 20.1201(a)(2)(ii) occupational dose limit, shallow-dose equivalent to skin, per annum[43]
    //1  Sv (100 rem):  Maximum allowed radiation exposure for NASA astronauts over their career[29]
    //4 to 5  Sv (400 to 500 rem):  Dose required to kill a human with a 50% risk within 30 days (LD50/30), if the dose is received over a very short duration.[46][47]
    
    int nprinted = 0, ntriedprint = 0;
    
    for( map<double,pair<string,WString> >::const_iterator iter = ref_doses.begin();
        iter != ref_doses.end(); ++iter )
    {
      ++ntriedprint;
      
      const double refdose = iter->first;
      const string &dosestr = iter->second.first;
      const WString &descestr = iter->second.second;
      
      const double time = refdose / dose;
      
      if( time < 1.0*PhysicalUnits::second && ntriedprint < 2 )
        continue;
      
      if( nprinted > 2 && time > 2.5*PhysicalUnits::year )
        break;
      
      const string timestr = PhysicalUnitsLocalized::printToBestTimeUnits(time);
      
      //Add to the table here
      WContainerWidget *labelDiv = new WContainerWidget();
      labelDiv->addStyleClass( "DoseStayLabelDiv" );
      WText *txt = new WText( dosestr, labelDiv );
      txt->setInline( false );
      txt->addStyleClass( "DoseStayDose" );
      
      txt = new WText( descestr, labelDiv );
      txt->setInline( false );
      txt->addStyleClass( "DoseStayDesc" );
      
      const int row = table->rowCount();
      table->addWidget( labelDiv, row, 0 );
      
      txt = new WText( timestr );
      txt->addStyleClass( "DoseStayTimeAnsw" );
      table->addWidget( txt, row, 1 );
      
      ++nprinted;
    }//for( loop over ref_doses )
    
  }catch( std::exception & )
  {
    //reset everything here.
    m_stayTime->removeStyleClass( "DoseStayTime" );
    m_stayTime->clear();
  }
}//void updateStayTime()

