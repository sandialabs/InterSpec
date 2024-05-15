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

#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WTable>
#include <Wt/WComboBox>
#include <Wt/WTableCell>
#include <Wt/WValidator>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WDoubleSpinBox>
#include <Wt/WContainerWidget>

#include "SpecUtils/StringAlgo.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/OneOverR2Calc.h"
#include "InterSpec/UndoRedoManager.h"

#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#include "InterSpec/WarningWidget.h"
#endif

using namespace Wt;
using namespace std;

OneOverR2Calc::OneOverR2Calc()
  : AuxWindow( WString::tr("window-title-1/r2-calc"),
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
               | AuxWindowProperties::SetCloseable
               | AuxWindowProperties::DisableCollapse) ),
    m_nearMeasurement( NULL ),
    m_farMeasurement( NULL ),
    m_backgroundMeasurment( NULL ),
    m_distance( NULL ),
    m_answer( NULL ),
    m_message( NULL ),
    m_prevValues{ 0.0f, 0.0f, 0.0f, 0.0f, 0.0f }
{
  wApp->useStyleSheet( "InterSpec_resources/OneOverR2Calc.css" );
  InterSpec::instance()->useMessageResourceBundle( "OneOverR2Calc" );
  
  WContainerWidget *contentDiv = contents();
  addStyleClass( "OneOverR2CalcWindow" );
  contentDiv->addStyleClass( "OneOverR2Calc" );
  
  WText *message = new WText( WString::tr("oor2c-instructions"), Wt::XHTMLText, contentDiv );
  message->setInline( false );
  message->addStyleClass( "OneOverR2IntroTxt" );

  WTable *layoutTable = new WTable( contentDiv );
  WTableCell *cell = layoutTable->elementAt( 0, 0 );
  WText *label = new WText( WString::tr("oor2c-near-meas-label"), cell );
  
  //TODO:  HelpSystem::attachToolTipOn( label,Intensity can be specified using any unit of measurement (ex. <b>rem</b>, <b>millirem</b>, <b>sievert/hour, gamma counts per second) as long as it is consistent among the fields. , showToolTips );
  
  cell = layoutTable->elementAt( 0, 1 );
  m_nearMeasurement = new WDoubleSpinBox( cell );
  m_nearMeasurement->setDecimals( 2 );
  m_nearMeasurement->setMinimum( 0.0 );
  m_nearMeasurement->setMaximum( 1.0E20 );
  m_nearMeasurement->changed().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->enterPressed().connect( this, &OneOverR2Calc::doCalc );

  cell = layoutTable->elementAt( 1, 0 );
  label = new WText( WString::tr("oor2c-far-meas-label"), cell );

  cell = layoutTable->elementAt( 1, 1 );
  m_farMeasurement = new WDoubleSpinBox( cell );
  m_farMeasurement->setDecimals( 2 );
  m_farMeasurement->setMinimum( 0.0 );
  m_farMeasurement->setMaximum( 1.0E20 );
  m_farMeasurement->changed().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->enterPressed().connect( this, &OneOverR2Calc::doCalc );

  cell = layoutTable->elementAt( 2, 0 );
  label = new WText( WString::tr("oor2c-background-intensity"), cell );

  cell = layoutTable->elementAt( 2, 1 );
  m_backgroundMeasurment = new WDoubleSpinBox( cell );
  m_backgroundMeasurment->setDecimals( 2 );
  m_backgroundMeasurment->setMinimum( 0.0 );
  m_backgroundMeasurment->setMaximum( 1.0E20 );
  m_backgroundMeasurment->changed().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->enterPressed().connect( this, &OneOverR2Calc::doCalc );
  
  cell = layoutTable->elementAt( 3, 0 );
  label = new WText( WString::tr("oor2c-dist-label"), cell );
  
  cell = layoutTable->elementAt( 3, 1 );
  m_distance = new WDoubleSpinBox( cell );
  m_distance->setDecimals( 2 );
  m_distance->setMinimum( 0.0 );
  m_distance->setMaximum( 9.99E6 );
  m_distance->changed().connect( this, &OneOverR2Calc::doCalc );
  m_distance->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_distance->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_distance->enterPressed().connect( this, &OneOverR2Calc::doCalc );

  
  //Power Law
  cell = layoutTable->elementAt( 4, 0 );
  cell->setColumnSpan( 2 );
  WContainerWidget *powerLawDiv = new WContainerWidget( cell );
  WGridLayout *powerLayout = new WGridLayout();
  powerLayout->setContentsMargins( 0, 0, 0, 0 );
  powerLayout->setVerticalSpacing( 0 );
  powerLayout->setHorizontalSpacing( 0 );
  powerLawDiv->setLayout( powerLayout );
  m_powerLawSelect = new WComboBox();
  label = new WText( WString::tr("oor2c-power-law-label") );
  powerLayout->addWidget( label, 0, 0 );
  powerLayout->addWidget( m_powerLawSelect, 0, 1 );
  m_powerLawSelect->addItem( WString::tr("oor2c-low-scatter") );
  m_powerLawSelect->addItem( WString::tr("oor2c-mid-scatter") );
  m_powerLawSelect->addItem( WString::tr("oor2c-high-scatter") );
  m_powerLawSelect->setCurrentIndex( 0 );
  m_powerLawSelect->activated().connect( this, &OneOverR2Calc::powerLawSelected );
  
  cell = layoutTable->elementAt( 5, 0 );
  label = new WText( WString::tr("oor2c-dist-to-near-label"), cell );

  cell = layoutTable->elementAt( 5, 1 );
  m_answer  = new WLineEdit( cell );
  m_answer->setDisabled( true );
  m_answer->setAttributeValue( "ondragstart", "return false" );

  m_backgroundMeasurment->setText( "" );

  message = new WText( WString::tr("oor2c-bottom-message"), Wt::XHTMLText, contents() );
  message->addStyleClass( "OneOverR2TextRow" );
  message->setInline( false );
  
  
  m_message = new WText( "&nbsp;", XHTMLText, contents() );
  m_message->setInline( false );
  m_message->addStyleClass( "OneOverR2Message" );

  m_message->setHiddenKeepsGeometry( true );
  
  AuxWindow::addHelpInFooter( footer(), "1/r2-calc-dialog" );
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://1overr2/?" + Wt::Utils::urlEncode(encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("oor2c-qr-state-title"),
                                 WString::tr("oor2c-qr-state-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("oor2c-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( boost::bind( &AuxWindow::hide, this ) );
  
  rejectWhenEscapePressed();
  show();
  
  InterSpec *viewer = InterSpec::instance();
  const bool isPhone = (viewer && viewer->isPhone());
  
  if( isPhone )
  {
    titleBar()->hide();
    
    if( viewer )
    {
      /* For some reason CSS seems to fail resizing AuxWindows properly, so we have to do it in c++ */
      
      float safeAreas[4] = { 0.0f };
#if( IOS )
      InterSpecApp *app = dynamic_cast<InterSpecApp *>(WApplication::instance());
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
#endif
      //repositionWindow( -32768, static_cast<int>(std::max(3.0f,0.5f*safeAreas[0])) );
      
      // TODO: right now hardcoding width because otherwise width will go to like full-screen
      const double width = 325;
      const double height = viewer->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f);
      setMaximumSize( width, height );
      
      /* ToDo: get safe offsets in c++ land, and then also convert other AuxWindows that are modal on phone to resize correctly. */
      /* Do same for Gamma XS Calc. And Energy Range Sum*/
      centerWindowHeavyHanded();
    }
  }else
  {
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }//if( isPhone ) / else
  
  //Keep the keyboard form popping up
  if( viewer && viewer->isMobile() )
  {
    closeButton->setFocus();
    closeButton->setFloatSide( Wt::Left ); //The "DialogClose" style class defaults to floating to the right, same as the help icon
  }
}//OneOverR2Calc constructor


OneOverR2Calc::~OneOverR2Calc()
{
  //nothing to do here
}//OneOverR2Calc destructor


void OneOverR2Calc::powerLawSelected()
{
  doCalc();
}//void powerLawSelected()


void OneOverR2Calc::doCalc()
{
  m_answer->setText( "" );

  if( m_nearMeasurement->validate() != WValidator::Valid )
  {
    m_message->setText( WString::tr("oor2c-invalid-near") );
    m_message->show();
    return;
  }//if( invalid near measurment )

  if( m_farMeasurement->validate() != WValidator::Valid )
  {
    m_message->setText( WString::tr("oor2c-invalid-far") );
    m_message->show();
    return;
  }//if( invalid far measurment )


  if( m_distance->validate() != WValidator::Valid )
  {
    m_message->setText( WString::tr("oor2c-invalid-dist") );
    m_message->show();
    return;
  }//if( invalid far measurment )

  double background = 0.0;
  if( m_backgroundMeasurment->validate() == WValidator::Valid )
    background = m_backgroundMeasurment->value();

  const double nearStrength = m_nearMeasurement->value() - background;
  const double farStrength = m_farMeasurement->value() - background;
  const double distance = m_distance->value();

  size_t nsigfigs = 4;
  auto nsigfigsLambda = []( WDoubleSpinBox *w ) -> size_t {
    string valstr = w->text().toUTF8();
    size_t ndig = 0;
    for( size_t i = 0; i < valstr.size(); ++i )
      ndig += (valstr[i]>='0' && valstr[i]<='9');
    return ndig;
  };
  nsigfigs = std::max( nsigfigs, nsigfigsLambda(m_nearMeasurement) );
  nsigfigs = std::max( nsigfigs, nsigfigsLambda(m_farMeasurement) );
  nsigfigs = std::max( nsigfigs, nsigfigsLambda(m_distance) );
  
  if( farStrength >= nearStrength )
  {
    m_message->setText( WString::tr("oor2c-far-larger-than-near") );
    m_message->show();
    return;
  }//if( near is larger than far )

  if( (background > 1E-6) && (farStrength <= 1E-6) )
  {
    m_message->setText( WString::tr("oor2c-back-larger-than-far") );
    m_message->show();
    return;
  }//if( background is bigger than the far measurment )

  if( nearStrength <= 1E-6 )
  {
    m_message->setText( WString::tr("oor2c-near-is-zero") );
    m_message->show();
    return;
  }//if( near strength is zero )

  if( farStrength <= 1E-6 )
  {
    m_message->setText( WString::tr("oor2c-far-is-zero") );
    m_message->show();
    return;
  }//if( far strength is zero )

  if( distance < 1E-6 )
  {
    m_message->setText( WString::tr("oor2c-dist-is-zero") );
    m_message->show();
    return;
  }//if( far strength is zero )

  //If weve made it this far, we have valid input
  m_message->setText( "&nbsp;" );
//  m_message->hide();

  /*
  nearStrength = src / r^2
  farStrength = src / (r+distance)^2
  nearStrength * r^2 = farStrength * (r+distance)^2
  nearStrength * r^2 = farStrength * ( r^2 + 2*distance*r + distance*distance)
  (farStrength - nearStrength) * r^2 + 2*farStrength*distance*r + farStrength*distance*distance = 0
  */

  double power = 2.0;
  switch( m_powerLawSelect->currentIndex() )
  {
    case 1: power = 1.85; break;
    case 2: power = 1.65; break;
    default: break;
  }//switch( m_powerLawSelect->currentIndex() )
  

//For 1/r2, would be
//  const double a = farStrength - nearStrength;
//  const double b = 2.0 * farStrength * distance;
//  const double c = farStrength * distance * distance;
//  const double r = (-b - sqrt(b*b - 4.0*a*c)) / (2*a);
////  const double other_r = (-b + sqrt(b*b - 4.0*a*c)) / (2*a);
  
  const double r = distance /( std::pow( (nearStrength/farStrength), 1.0/power) - 1.0 );
  
  char answerStr[64];
  snprintf( answerStr, sizeof(answerStr), ("%." + std::to_string(nsigfigs) + "g").c_str(), r );
  m_answer->setText( answerStr );
  
  // Now deal with undo/redo
  std::array<float,5> currentValues{
    static_cast<float>( m_nearMeasurement->value() ),
    static_cast<float>( m_farMeasurement->value() ),
    static_cast<float>( m_backgroundMeasurment->value() ),
    static_cast<float>( m_distance->value() ),
    static_cast<float>( m_powerLawSelect->currentIndex() )
  };
  
  UndoRedoManager *undoRedo = UndoRedoManager::instance();
  if( undoRedo && undoRedo->canAddUndoRedoNow() && (currentValues != m_prevValues) )
  {
    std::array<float,5> prevValues = m_prevValues;
    auto undo_redo = [currentValues,prevValues]( const bool is_undo ){
      OneOverR2Calc *tool = InterSpec::instance()->createOneOverR2Calculator();
      if( !tool )
        return;
      const std::array<float,5> &val = is_undo ? prevValues : currentValues;
      tool->m_prevValues = val;
      tool->m_nearMeasurement->setValue( val[0] );
      tool->m_farMeasurement->setValue( val[1] );
      tool->m_backgroundMeasurment->setValue( val[2] );
      tool->m_distance->setValue( val[3] );
      tool->m_powerLawSelect->setCurrentIndex( static_cast<int>(val[4]) );
      tool->doCalc();
    };
    
    auto undo = [undo_redo](){ undo_redo(true); };
    auto redo = [undo_redo](){ undo_redo(false); };
    
    undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "1/r2 value change." );
  }//if( value != m_prevValues )
  
  m_prevValues = std::move(currentValues);
}//void doCalc()


void OneOverR2Calc::handleAppUrl( std::string query_str )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const map<string,string> values = AppUtils::query_str_key_values( query_str );
  
  // Check version is appropriate
  const auto ver_iter = values.find( "V" );
  if( (ver_iter == end(values))
     || ((ver_iter->second != "1") && !SpecUtils::istarts_with(ver_iter->second, "1.")) )
    throw runtime_error( "1/r2 Calc.: URI not compatible version." );
  
  
  
  // A boilerplate lambda for grapping and validating the fields
  auto getField = [&]( const string &name, Wt::WLineEdit *edit ){
    const auto i = values.find( name );
    string value = (i == end(values)) ? string() : i->second;
    SpecUtils::ireplace_all( value, "%23", "#" );  //Shouldnt be an issue, but jic
    SpecUtils::ireplace_all( value, "%26", "&" ); //Shouldnt be an issue, but jic
    
    // Let the value be empty,
    if( !value.empty() )
    {
      try
      {
        const double test_val = std::stod(value);
        if( IsInf(test_val) || IsNan(test_val) )
          throw runtime_error( "NaN or Inf" );
      }catch( std::exception & )
      {
        throw runtime_error( "1/r2 Calc.: invalid '" + name + "' value '" + value + "' in URI." );
      }// try / catch
    }//if( !value.empty() )
    
    edit->setValueText( WString::fromUTF8(value) );
  };//getField(...)
  
  
  getField( "NEAR", m_nearMeasurement );
  getField( "FAR", m_farMeasurement );
  getField( "BACK", m_backgroundMeasurment );
  getField( "DIST", m_distance );
  
  int power_index = 0;
  const auto pow_iter = values.find( "POW" );
  if( pow_iter != end(values) )
  {
    if( pow_iter->second == "1.85" )
      power_index = 1;
    else if( pow_iter->second == "1.65" )
      power_index = 2;
  }//if( pow_iter != end(values) )
  
  m_powerLawSelect->setCurrentIndex( power_index );
  
  doCalc();
}//void handleAppUrl( std::string query_str )


std::string OneOverR2Calc::encodeStateToUrl() const
{
  // "near=3.5&far=1.1&back=2&dit=1.9&power=1.75"
  string answer = "V=1";
  
  auto addField = [&answer]( const string &name, Wt::WLineEdit *edit ){
    string val = edit->text().toUTF8();
    SpecUtils::ireplace_all(val, "#", "%23" ); //Shouldnt be an issue, but JIC
    SpecUtils::ireplace_all(val, "&", "%26" ); //Shouldnt be an issue, but JIC
    answer += "&" + name + "=" + val;
  };
  
  addField( "NEAR", m_nearMeasurement );
  addField( "FAR", m_farMeasurement );
  addField( "BACK", m_backgroundMeasurment );
  addField( "DIST", m_distance );

  string power = "2";
  switch( m_powerLawSelect->currentIndex() )
  {
    case 1: power = "1.85"; break;
    case 2: power = "1.65"; break;
    default: break;
  }//switch( m_powerLawSelect->currentIndex() )
  answer += "&POW=" + power;
  
  return answer;
}//string encodeStateToUrl() const
