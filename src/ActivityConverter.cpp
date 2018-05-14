/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#include <iostream>

#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WComboBox>
#include <Wt/WLineEdit>
#include <Wt/WValidator>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>

#include "InterSpec/AuxWindow.h"
#include "InterSpec/PhysicalUnits.h"
#include "SpecUtils/UtilityFunctions.h"
#include "InterSpec/ActivityConverter.h"

using namespace Wt;
using namespace std;

ActivityConverter::ActivityConverter()
  : AuxWindow( "Activity Converter" ),
    m_bq( NULL ),
    m_ci( NULL ),
    m_message( NULL )
{
  addStyleClass( "ActivityConverter" );
  
  WText *message = NULL;
  WGridLayout* layout =  stretcher();
  setWidth(400);
  const char *topMessage = "Convert between curie and becquerel (ie. 5 MBq, 2 nCi)";
  message = new WText( topMessage, Wt::XHTMLUnsafeText);

  layout->addWidget(message,0,0,1,3);

  WLabel *label = new WLabel( "From: ");
  layout->addWidget(label,1,0);
  label->addStyleClass( "ActivityConverterLabel" );
  
  m_bq = new WLineEdit(  );
    layout->addWidget(m_bq,1,1);

  const char * const bqexp = "^(\\s*\\+?((\\d+(\\.\\d*)?)|(\\.\\d*))"
    "(?:[Ee][+\\-]?\\d+)?\\s*"
    "(b|bq|Bq|k|kB|kBq|M|MB|MBq|G|GB|GBq|T|TB|TBq|c|ci|m|mC|mCi|m|mi|mic|micr|micro|microC|microCi|n|nC|nCi))+\\s*";
    
  WRegExpValidator *validator = new WRegExpValidator( bqexp, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_bq->setValidator(validator);
  
  m_bq->addStyleClass( "ActivityConverterInput" );
  m_bq->setTextSize( 15 );
  m_bq->setWidth( WLength(15.0,WLength::FontEm) );
  m_bq->setText("5 MBq");
  m_bq->changed().connect(  this, &ActivityConverter::doBq);
  m_bq->blurred().connect(   this, &ActivityConverter::doBq);
  m_bq->enterPressed().connect(    this, &ActivityConverter::doBq);
  
  WPushButton *convertButton = new WPushButton("Convert");
  layout->addWidget(convertButton,1,2);
  convertButton->clicked().connect(this, &ActivityConverter::doBq);
    
  label = new WLabel( "To: " );
  layout->addWidget(label,2,0);
  label->addStyleClass( "ActivityConverterLabel" );
  
  m_ci = new WLineEdit( );
  layout->addWidget(m_ci,2,1,1,2);
  m_ci->addStyleClass( "ActivityConverterInput" );
  m_ci->setTextSize( 15 );
  m_ci->setWidth( WLength(15.0,WLength::FontEm) );
  m_ci->setEnabled(false);
  m_ci->setText("135.14 uCi");
  
  m_message = new WText( "&nbsp", XHTMLUnsafeText );
  m_message->setHeight(WLength(50,WLength::Pixel));
  m_message->setAttributeValue( "style", "color:blue;"  );
//  m_message->setHiddenKeepsGeometry( true );
    m_message->hide();
  layout->addWidget(m_message,3,0,1,3);
  layout->setColumnStretch(1, 1);
  layout->setRowStretch(3, 1);
    
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );

  rejectWhenEscapePressed();
  finished().connect( boost::bind( &AuxWindow::deleteAuxWindow, this ) );
  
  centerWindow();
    resizeToFitOnScreen();
  show();

  setResizable( false );
  
}//ActivityConverter constructor


ActivityConverter::~ActivityConverter()
{
  //nothing to do here
}//ActivityConverter destructor


void ActivityConverter::doBq()
{
  try {
    Wt::WString value = m_bq->text();
    std::string val = value.toUTF8();
    UtilityFunctions::trim( val );
    double dbvalue = PhysicalUnits::stringToActivity(val);
    bool curie = false;
    if (val.find("Ci")!= std::string::npos || val.find("ci")!= std::string::npos) {
        curie=false;
    } else if (val.find("Bq")!= std::string::npos || val.find("bq")!= std::string::npos) {
        curie=true;
    } else {
      //Invalid 
      throw val;
    }
    const string ans = PhysicalUnits::printToBestActivityUnits( dbvalue, 2, curie, PhysicalUnits::becquerel );
    m_ci->setText(Wt::WString::fromUTF8(ans));
    //If weve made it this far, we have valid input
    m_message->removeStyleClass("line-above");
    m_message->setText( "&nbsp" );
    m_message->hide();
  } catch (...) {
    m_message->setText("Invalid value.  Allowed units: bq, kBq, MBq, GBq, TBq, ci, mCi, microCi, nCi" );
    m_message->addStyleClass("line-above");
    m_message->show();
  }

}
