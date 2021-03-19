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

#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"
#include "InterSpec/HelpSystem.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/OneOverR2Calc.h"

using namespace Wt;
using namespace std;

OneOverR2Calc::OneOverR2Calc()
  : AuxWindow( "1/r<sup>2</sup> Calculator",
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
               | AuxWindowProperties::SetCloseable
               | AuxWindowProperties::DisableCollapse) ),
    m_nearMeasurement( NULL ),
    m_farMeasurement( NULL ),
    m_backgroundMeasurment( NULL ),
    m_distance( NULL ),
    m_answer( NULL ),
    m_message( NULL )
{
  wApp->useStyleSheet( "InterSpec_resources/OneOverR2Calc.css" );
  
  WContainerWidget *contentDiv = contents();
  addStyleClass( "OneOverR2CalcWindow" );
  contentDiv->addStyleClass( "OneOverR2Calc" );
  
  WText *message = NULL;

  const char *topMessage = "Use two measurement at different locations to find "
                           "distance to an unseen source. E.g. when "
                            "the source is behind a wall.";
  
  message = new WText( topMessage, Wt::XHTMLText, contentDiv );
  message->setInline( false );
  message->addStyleClass( "OneOverR2IntroTxt" );

  WTable *layoutTable = new WTable( contentDiv );
  WTableCell *cell = layoutTable->elementAt( 0, 0 );
  WText *label = new WText( "Near Measurement Intensity:", cell );
  
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
  label = new WText( "Far Measurement Intensity:", cell );

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
  label = new WText( "Background Intensity (optional):", cell );

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
  label = new WText( "Distance between measurements:", cell );
  
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
  label = new WText( "Power Law");
  powerLayout->addWidget( label, 0, 0 );
  powerLayout->addWidget( m_powerLawSelect, 0, 1 );
  m_powerLawSelect->addItem( "Low Scatter or using Peak Area, 1/r^2" );
  m_powerLawSelect->addItem( "Mid Scatter Dose Rate, 1/r^1.85" );
  m_powerLawSelect->addItem( "High Scatter Dose Rate, 1/r^1.65" );
  m_powerLawSelect->setCurrentIndex( 0 );
  m_powerLawSelect->activated().connect( this, &OneOverR2Calc::powerLawSelected );
  
  cell = layoutTable->elementAt( 5, 0 );
  label = new WText( "Dist. near measurement to source:", cell );

  cell = layoutTable->elementAt( 5, 1 );
  m_answer  = new WLineEdit( cell );
  m_answer->setDisabled( true );

  m_backgroundMeasurment->setText( "" );

  const char *bottomMessage = "Use the same units for near, background, and far "
                              "measurements. Results are in same units used for "
                              "distance between measurements.";
  message = new WText( bottomMessage, Wt::XHTMLText, contents() );
  message->addStyleClass( "OneOverR2TextRow" );
  message->setInline( false );
  
  
  m_message = new WText( "&nbsp;", XHTMLText, contents() );
  m_message->setInline( false );
  m_message->addStyleClass( "OneOverR2Message" );

  m_message->setHiddenKeepsGeometry( true );
  
  AuxWindow::addHelpInFooter( footer(), "1/r2-calc-dialog" );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  rejectWhenEscapePressed();
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  
  
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(WApplication::instance());
  const bool isPhone = (app && app->isPhone());
  
  if( isPhone )
  {
    titleBar()->hide();
    
    InterSpec *viewer = app->viewer();
    if( viewer )
    {
      /* For some reason CSS seems to fail resizing AuxWindows properly, so we have to do it in c++ */
      
      float safeAreas[4] = { 0.0f };
#if( IOS )
      InterSpecApp::DeviceOrientation orientation = InterSpecApp::DeviceOrientation::Unknown;
      app->getSafeAreaInsets( orientation, safeAreas[0], safeAreas[1], safeAreas[2], safeAreas[3] );
#endif
      repositionWindow( -32768, static_cast<int>(std::max(3.0f,0.5f*safeAreas[0])) );
      setMaximumSize( WLength::Auto, viewer->renderedHeight() - std::max(0.5f*(safeAreas[0]+safeAreas[2]),6.0f) );
      
      /* ToDo: get safe offsets in c++ land, and then also convert other AuxWindows that are modal on phone to resize correctly. */
      /* Do same for Gamma XS Calc. And Energy Range Sum*/
    }
  }else
  {
    resizeToFitOnScreen();
    centerWindowHeavyHanded();
  }//if( isPhone ) / else
  
  //Keep the keyboard form popping up
  if( app && app->isMobile() )
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
    m_message->setText( "Invalid near measurment" );
    m_message->show();
    return;
  }//if( invalid near measurment )

  if( m_farMeasurement->validate() != WValidator::Valid )
  {
    m_message->setText( "Invalid far measurment" );
    m_message->show();
    return;
  }//if( invalid far measurment )


  if( m_distance->validate() != WValidator::Valid )
  {
    m_message->setText( "Invalid distance" );
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
    m_message->setText( "Far measurement must be less than near one" );
    m_message->show();
    return;
  }//if( near is larger than far )

  if( (background > 1E-6) && (farStrength <= 1E-6) )
  {
    m_message->setText( "The background must be less than the far one" );
    m_message->show();
    return;
  }//if( background is bigger than the far measurment )

  if( nearStrength <= 1E-6 )
  {
    m_message->setText( "The near measurement must be non-zero" );
    m_message->show();
    return;
  }//if( near strength is zero )

  if( farStrength <= 1E-6 )
  {
    m_message->setText( "The far measurement must be non-zero" );
    m_message->show();
    return;
  }//if( far strength is zero )

  if( distance < 1E-6 )
  {
    m_message->setText( "Distance between measurements must be non-zero" );
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
}//void doCalc()
