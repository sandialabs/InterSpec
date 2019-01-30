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
#include <Wt/WComboBox>
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
              (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::TabletModal)
               | AuxWindowProperties::SetCloseable
               | AuxWindowProperties::DisableCollapse) ),
    m_nearMeasurement( NULL ),
    m_farMeasurement( NULL ),
    m_backgroundMeasurment( NULL ),
    m_distance( NULL ),
    m_answer( NULL ),
    m_message( NULL )
{
  wApp->useStyleSheet( "InterSpec_resources/OneOverR2CalcInput.css" );
  
//  addStyleClass( "OneOverR2Calc" );
  
  WText *message = NULL;
  WGridLayout * layout = stretcher();

  const char *topMessage = "Use two measurement at different locations to find "
                           "distance to an unknown source location. E.g. when "
                            "the source is behind a wall.";
  message = new WText( topMessage, Wt::XHTMLUnsafeText);
  layout->addWidget(message,0,0,1,2);

  WText *label = new WText( "Near Measurement Intensity:" );
  layout->addWidget(label,1,0);
  
  //TODO:  HelpSystem::attachToolTipOn( label,Intensity can be specified using any unit of measurement (ex. <b>rem</b>, <b>millirem</b>, <b>sievert/hour, gamma counts per second) as long as it is consistent among the fields. , showToolTipInstantly );
  
  m_nearMeasurement = new WDoubleSpinBox(  );
  layout->addWidget(m_nearMeasurement,1,1);
  m_nearMeasurement->setDecimals( 2 );
  m_nearMeasurement->setMinimum( 0.0 );
  m_nearMeasurement->setMaximum( 1.0E20 );
  m_nearMeasurement->changed().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_nearMeasurement->enterPressed().connect( this, &OneOverR2Calc::doCalc );


  label = new WText( "Far Measurement Intensity:");
  layout->addWidget(label,2,0);

  m_farMeasurement = new WDoubleSpinBox(  );
  layout->addWidget(m_farMeasurement,2,1);
  m_farMeasurement->setDecimals( 2 );
  m_farMeasurement->setMinimum( 0.0 );
  m_farMeasurement->setMaximum( 1.0E20 );
  m_farMeasurement->changed().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_farMeasurement->enterPressed().connect( this, &OneOverR2Calc::doCalc );

  label = new WText( "Background Intensity (optional):" );
  layout->addWidget(label,3,0);

  m_backgroundMeasurment = new WDoubleSpinBox();
  layout->addWidget(m_backgroundMeasurment,3,1);
  m_backgroundMeasurment->setDecimals( 2 );
  m_backgroundMeasurment->setMinimum( 0.0 );
  m_backgroundMeasurment->setMaximum( 1.0E20 );
  m_backgroundMeasurment->changed().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_backgroundMeasurment->enterPressed().connect( this, &OneOverR2Calc::doCalc );
  label = new WText( "Distance between measurements:" );
  layout->addWidget(label,4,0);
  m_distance = new WDoubleSpinBox( );
  layout->addWidget(m_distance,4,1);
  m_distance->setDecimals( 2 );
  m_distance->setMinimum( 0.0 );
  m_distance->setMaximum( 9.99E6 );
  m_distance->changed().connect( this, &OneOverR2Calc::doCalc );
  m_distance->valueChanged().connect( this, &OneOverR2Calc::doCalc );
  m_distance->blurred().connect( this, &OneOverR2Calc::doCalc );
  m_distance->enterPressed().connect( this, &OneOverR2Calc::doCalc );

  
  //Power Law
  WContainerWidget *powerLawDiv = new WContainerWidget();
  layout->addWidget( powerLawDiv, 5, 0, 1, 2 );
  WGridLayout *powerLayout = new WGridLayout();
  powerLawDiv->setLayout( powerLayout );
  m_powerLawSelect = new WComboBox();
  label = new WText( "Power Law");
  powerLayout->addWidget( label, 0, 0 );
  powerLayout->addWidget( m_powerLawSelect, 0, 1 );
  m_powerLawSelect->addItem( WString::fromUTF8("Low Scatter or using Peak Area, 1/r\u00B2") );
  m_powerLawSelect->addItem( "Mid Scatter Dose Rate, 1/r^1.85" );
  m_powerLawSelect->addItem( "High Scatter Dose Rate, 1/r^1.65" );
  m_powerLawSelect->setCurrentIndex( 0 );
  m_powerLawSelect->activated().connect( this, &OneOverR2Calc::powerLawSelected );
  
  label = new WText( "Distance of source in front of near measurement");
  layout->addWidget(label,6,0);
  m_answer  = new WLineEdit( );
  layout->addWidget(m_answer,6,1);
  m_answer->setDisabled( true );

  m_backgroundMeasurment->setText( "" );

  const char *bottomMessage = "Use the same units for near, background, and far "
                              "measurements. Results are in same units used for "
                              "distance between measurements.";
  message = new WText( bottomMessage, Wt::XHTMLUnsafeText );
  layout->addWidget( message, 7, 0, 1, 2 );
  
  m_message = new WText( "&nbsp", XHTMLUnsafeText);
  layout->addWidget( m_message, 8, 0, 1, 2 );
  layout->setRowStretch( 8, 1 );
  layout->setColumnStretch( 0, 1 );
  string styleVal = m_message->attributeValue("style").toUTF8();
  m_message->setAttributeValue( "style", "color:blue;" + styleVal );
  m_message->setHiddenKeepsGeometry( true );
  
  AuxWindow::addHelpInFooter( footer(), "1/r2-calc-dialog", this );
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  rejectWhenEscapePressed();
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  show();
  centerWindow();
  
  //Keep the keyboard form popping up
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(WApplication::instance());
  if( app && app->isMobile() )
    closeButton->setFocus();
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

  if( farStrength >= nearStrength )
  {
    m_message->setText( "Far radiation measurement must be smaller than "
                        "near measurment" );
    m_message->show();
    return;
  }//if( near is larger than far )

  if( (background > 1E-6) && (farStrength <= 1E-6) )
  {
    m_message->setText( "The background must be less than the far measurement" );
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
  m_message->setText( "&nbsp" );
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
  snprintf( answerStr, sizeof(answerStr), "%.2g", r );
  m_answer->setText( answerStr );
}//void doCalc()
