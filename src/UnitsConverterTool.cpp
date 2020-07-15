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
#include <iostream>

#include <boost/regex.hpp>

#include <Wt/WText>
#include <Wt/WBreak>
#include <Wt/WLabel>
#include <Wt/WLineEdit>
#include <Wt/WValidator>
#include <Wt/WPushButton>
#include <Wt/WGridLayout>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>

#include "InterSpec/AuxWindow.h"
#include "SpecUtils/StringAlgo.h"
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UnitsConverterTool.h"

using namespace Wt;
using namespace std;

namespace
{

int num_sig_figs( const string &val )
{
  // \TODO: this isnt necassilry correct, because it is being used to determine number of places
  //        after a decimal... we'll worry about this later.
  boost::smatch mtch;
  boost::regex expr( PhysicalUnits::sm_positiveDecimalRegex + string(".*"));
  
  if( boost::regex_match( val, mtch, expr ) )
  {
    assert( mtch.size() == 5 );
    //cout << "Matching fields: " << endl;
    //for( size_t i = 0; i < mtch.size(); ++i )
    //  cout << "\tField " << i << ": " << string( mtch[i].first, mtch[i].second ) << endl;
    
    string numberfield( mtch[1].first, mtch[1].second );
    auto dec_pos = numberfield.find(".");
    if( dec_pos == string::npos )
    {
      while( numberfield.size() && numberfield.back()=='0' )
        numberfield = numberfield.substr(0,numberfield.size()-1);
    }
    
    int nnum = 0;
    for( auto ch : numberfield )
      nnum += isdigit(ch);
    
    return std::min( nnum, 2 );
  }
  
  return 3;
}//num_sig_figs( const string &val )


string convertActivity( const string &val )
{
  try
  {
    const double dbvalue = PhysicalUnits::stringToActivity(val);
    
    boost::smatch mtch;
    boost::regex expr( string("(") + PhysicalUnits::sm_positiveDecimalRegex + string(")" "\\s*([a-zA-Z \\-]+)") );

    bool input_in_bq = false; // input_in_ci = false;
      
    if( boost::regex_match( val, mtch, expr ) )
    {
      assert( mtch.size() == 7 );
      //string number = string( mtch[1].first, mtch[1].second );
      string letters = string( mtch[6].first, mtch[6].second );
        
      //PhysicalUnits::stringToActivity() should make sure is either bq or ci
      input_in_bq = SpecUtils::icontains(letters, "b" );
      //input_in_ci = SpecUtils::icontains(letters, "c" );
    }//if( boost::regex_match( val, mtch, expr ) )
      
    return PhysicalUnits::printToBestActivityUnits( dbvalue, num_sig_figs(val), input_in_bq );
  }catch(...)
  {
  }
  
  return "";
}//string convertActivity( string val )


string convertAbsorbedDose( const string &val )
{
  try
  {
    const double dbvalue = PhysicalUnits::stringToAbsorbedDose(val);
    const bool to_gray = (SpecUtils::icontains(val, "rad") || SpecUtils::icontains(val, "erg"));
    return PhysicalUnits::printToBestAbsorbedDoseUnits( dbvalue, num_sig_figs(val), to_gray );
  }catch(...)
  {
  }
  
  return "";
}//string convertAbsorbedDose( string val )


string convertEquivalentDose( const string &val )
{
  try
  {
    const double dbvalue = PhysicalUnits::stringToEquivalentDose(val);
    const bool to_sievert = SpecUtils::icontains(val, "rem");
    return PhysicalUnits::printToBestEquivalentDoseUnits( dbvalue, num_sig_figs(val), to_sievert );
  }catch(...)
  {
  }
  
  return "";
}//string convertEquivalentDose( string val )

string convertDistance( string val )
{
  try
  {
    const double dbvalue = PhysicalUnits::stringToDistance(val);
    //const bool to_metric = (SpecUtils::icontains(val, "ft") || SpecUtils::icontains(val, "feet")
    //                        || SpecUtils::icontains(val, "in") || SpecUtils::icontains(val, "'")
    //                        || SpecUtils::icontains(val, "\"") );
    return PhysicalUnits::printToBestLengthUnits( dbvalue, num_sig_figs(val) );
  }catch(...)
  {
  }
  
  return "";
}//string convertDistance( string val )
}//namespace


UnitsConverterTool::UnitsConverterTool()
  : AuxWindow( "Units Converter",
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneModal)
               | AuxWindowProperties::DisableCollapse
               | AuxWindowProperties::SetCloseable) ),
    m_input( NULL ),
    m_output( NULL ),
    m_message( NULL )
{
  //addStyleClass( "UnitsConverterTool" );
  
  WGridLayout *layout = stretcher();
  setWidth(400);
  const char *topMessage = "Convert between radiation related units.<br />"
                           "Ex: 5 MBq, 2 nCi, 1.2rad, 15E-3gy, 0.2mrem, 8feet, 9milli-sievert";
  WText *message = new WText( topMessage, Wt::XHTMLUnsafeText);

  layout->addWidget(message,0,0,1,3);

  WLabel *label = new WLabel( "Input: ");
  layout->addWidget(label,1,0);
  //label->addStyleClass( "UnitsConverterToolLabel" );
  
  m_input = new WLineEdit();
  layout->addWidget(m_input,1,1);

  string regex = string("((") + PhysicalUnits::sm_absorbedDoseRegex
                      + ")|(" + PhysicalUnits::sm_activityRegex
                      + ")|(" + PhysicalUnits::sm_equivalentDoseRegex
                      + ")|(" + PhysicalUnits::sm_distanceRegex
                      + "))";
  WRegExpValidator *validator = new WRegExpValidator( regex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_input->setValidator(validator);
  
  //m_input->addStyleClass( "UnitsConverterToolInput" );
  m_input->setTextSize( 15 );
  m_input->setWidth( WLength(15.0,WLength::FontEm) );
  m_input->setText("5 MBq");
  //m_input->changed().connect( this, &UnitsConverterTool::convert ); //the signal is only emitted when the focus is lost
  m_input->blurred().connect( this, &UnitsConverterTool::convert );
  m_input->keyWentUp().connect( this, &UnitsConverterTool::convert );
  //m_input->enterPressed().connect( this, &UnitsConverterTool::convert );
  
  //WPushButton *convertButton = new WPushButton("Convert");
  //layout->addWidget(convertButton,1,2);
  //convertButton->clicked().connect( this, &UnitsConverterTool::convert );
    
  label = new WLabel( "Output: " );
  layout->addWidget(label,2,0);
  //label->addStyleClass( "UnitsConverterToolLabel" );
  
  m_output = new WLineEdit( );
  //layout->addWidget(m_output,2,1,1,2);
  layout->addWidget( m_output, 2, 1 );
  //m_output->addStyleClass( "UnitsConverterToolInput" );
  m_output->setTextSize( 15 );
  m_output->setWidth( WLength(15.0,WLength::FontEm) );
  m_output->setEnabled(false);
  //m_output->setText("135.14 uCi");
  
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
  
  convert();
  
  show();
  centerWindow();
  resizeToFitOnScreen();

  //Keep the keyboard form popping up
  InterSpecApp *app = dynamic_cast<InterSpecApp *>(WApplication::instance());
  if( app && app->isMobile() )
  {
    closeButton->setFocus();
    titleBar()->hide();
  }
  
}//UnitsConverterTool constructor


UnitsConverterTool::~UnitsConverterTool()
{
  //nothing to do here
}//UnitsConverterTool destructor


void UnitsConverterTool::convert()
{
  try
  {
    switch( m_input->validate() )
    {
      case Wt::WValidator::Invalid:
        throw std::runtime_error( "Invalid input" );
        
      case Wt::WValidator::InvalidEmpty:
        throw runtime_error( "" );
        
      case Wt::WValidator::Valid:
        break;
    }//switch( m_input->validate() )
    
    std::string val = m_input->text().toUTF8();
    SpecUtils::trim( val );
    
    if( val.empty() )
      throw runtime_error( "Empty input" );
    
    string ans = convertActivity( val );
    
    if( ans.empty() )
      ans = convertAbsorbedDose( val );
    
    if( ans.empty() )
      ans = convertEquivalentDose( val );
    
    if( ans.empty() )
      ans = convertDistance( val );
    
    if( ans.empty() )
    {
      boost::smatch mtch;
      boost::regex expr( PhysicalUnits::sm_positiveDecimalRegex );
      if( !boost::regex_match( val, mtch, expr ) )
        throw runtime_error( "Couldnt determine units" );
      throw runtime_error( "Couldnt convert number" );
    }//if( ans.empty() )
    
    m_output->setText( Wt::WString::fromUTF8(ans) );
    
    //If weve made it this far, we have valid input
    m_message->removeStyleClass("line-above");
    m_message->setText( "&nbsp" );
    m_message->hide();
  }catch( std::exception &e )
  {
    string errmsg = e.what();
    if( !errmsg.empty() )
      errmsg = "Error: " + errmsg + ".<br />";
    errmsg += "Allowed activity units: bq, becquerel, ci, curie, gray, Gy, rad, sievert, Sv, rem,"
              " roentgen, m, meter, ft, inch, and in.<br />"
              "Allowed (optional) prefixes: n, nano, u, \xc2\xb5, micro, m, milli,"
              " k, killo, M, mega, G, Giga, T, and Tera";
 
    m_output->setText( "" );
    m_message->setText( WString::fromUTF8(errmsg) );
    //m_message->addStyleClass("line-above");
    m_message->show();
  }//try / catch
}//void UnitsConverterTool::convert()
