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
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WContainerWidget>
#include <Wt/WRegExpValidator>

#include "SpecUtils/StringAlgo.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/AppUtils.h"
#include "InterSpec/AuxWindow.h"
#include "InterSpec/InterSpec.h"      // Only for preferenceValue<bool>("DisplayBecquerel")
#include "InterSpec/InterSpecApp.h"
#include "InterSpec/InterSpecUser.h"  // Only for preferenceValue<bool>("DisplayBecquerel")
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/UndoRedoManager.h"
#include "InterSpec/UnitsConverterTool.h"
#include "InterSpec/DecayDataBaseServer.h"


#if( USE_QR_CODES )
#include <Wt/Utils>

#include "InterSpec/QrCode.h"
#include "InterSpec/WarningWidget.h"
#endif


using namespace Wt;
using namespace std;

namespace
{

int num_sig_figs( const string &val )
{
  // We'll loop over matches, and take the one with the most sig-figs (but at least 2)
  int num_sig_figs = 2;
  bool found_match = false;
  
  const char *dec_regex = "[\\+-]?\\s*((\\d+(\\.\\d*)?)|(\\.\\d*))\\s*(?:[Ee][+\\-]?\\d+)?";
  boost::regex expr( dec_regex );
  boost::sregex_iterator match_iter(begin(val), end(val), expr);
  const boost::sregex_iterator match_end{};
  
  //boost::smatch mtch;
  //if( boost::regex_match( val, mtch, expr ) )
  for( ; match_iter != match_end; ++match_iter )
  {
    found_match = true;
    
    boost::smatch mtch = *match_iter;
    assert( mtch.size() == 5 );
    //cout << "Matching fields: " << endl;
    //for( size_t i = 0; i < mtch.size(); ++i )
    //  cout << "\tField " << i << ": " << string( mtch[i].first, mtch[i].second ) << endl;
    
    // Get the number field
    string numberfield( mtch[1].first, mtch[1].second );
    
    // remove leading zeros
    while( !numberfield.empty() && (numberfield[0] == '0') )
      numberfield = numberfield.substr(1);
    
    int nnum = 0;
    for( auto ch : numberfield )
      nnum += isdigit(ch);
    
    num_sig_figs = std::max( nnum, num_sig_figs );
  }//for( ; match_iter != boost::sregex_iterator{}; ++match_iter )
  
  //cout << "'" << val << "' --> " << (found_match ? num_sig_figs : 3) << endl;
  return found_match ? num_sig_figs : 3;
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
    const bool to_metric = (SpecUtils::icontains(val, "ft") || SpecUtils::icontains(val, "feet")
                            || SpecUtils::icontains(val, "in") || SpecUtils::icontains(val, "'")
                            || SpecUtils::icontains(val, "\"") );
    
    if( to_metric )
      return PhysicalUnits::printToBestLengthUnits( dbvalue, num_sig_figs(val) + 1 );
    
    // Else, to inches, ft, etc.
    //  (we wont sully PhysicalUnits with outputting English units)
    double unitval = dbvalue / (2.54 * PhysicalUnits::cm);
    const char *unit = "inch";
    if( unitval < 36*2.54 )
    {
      // Already taken care of
    }else if( unitval < 12*5280 )
    {
      unitval = unitval / 12;
      unit = "feet";
    }else
    {
      unitval = unitval / (12*5280);
      unit = "miles";
    }
    
    char formatflag[32], buffer[64];
    snprintf(formatflag, sizeof(formatflag), "%%.%if %%s", num_sig_figs(val) );
    snprintf(buffer, sizeof(buffer), formatflag, unitval, unit );
    
    return buffer;
  }catch(...)
  {
  }
  
  return "";
}//string convertDistance( string val )


string convertMass( string val )
{
  try
  {
    const double dbvalue = PhysicalUnits::stringToMass(val);
    
    const bool to_metric = (SpecUtils::icontains(val, "ounce") || SpecUtils::icontains(val, "oz")
                            || SpecUtils::icontains(val, "pound") || SpecUtils::icontains(val, "lb")
                            || SpecUtils::icontains(val, "stone") || SpecUtils::icontains(val, "grain") );
    
    if( to_metric )
      return PhysicalUnits::printToBestMassUnits( dbvalue, num_sig_figs(val) + 1 );
    
    
    const double unitval = 2.20462 * dbvalue / PhysicalUnits::kilogram;
    const char *unit = "lbs";
    
    char formatflag[32], buffer[64];
    snprintf(formatflag, sizeof(formatflag), "%%.%if %%s", num_sig_figs(val) );
    snprintf(buffer, sizeof(buffer), formatflag, unitval, unit );
    
    return buffer;
  }catch(...)
  {
  }
  
  return "";
}//string convertMass( string val )


/** Returns pointer for detected nuclide, as well as input string with nuclide portion of it removed.
 
 On failure to find a nuclide, returns nullptr, and original input string.
 
 Will throw exception if multiple nuclides input.
 */
pair<const SandiaDecay::Nuclide *, std::string> detect_nuclide( std::string val )
{
  const SandiaDecay::Nuclide *nuc = nullptr;
  
  const auto db = DecayDataBaseServer::database();
  vector<string> parts;
  SpecUtils::split( parts, val, " \t\n\r" );
  
  // Try single words, Like "U238", "Co-60", etc
  for( size_t i = 0; i < parts.size(); ++i )
  {
    const string &nucstr = parts[i];
    auto nucptr = db->nuclide( nucstr );
    if( nucptr && !nucptr->isStable() )
    {
      // Do a feeble check to make sure the user didnt input multiple nuclides
      if( nuc )
        throw runtime_error( WString::tr("uct-mult-nucs").arg(nuc->symbol).arg(nucptr->symbol).toUTF8() );
      
      nuc = nucptr;
      
      // Reconstruct `val` with the parts we used removed; we arent simply removing all `nucstr`
      //  substrings, jic they are repeated
      val = "";
      for( size_t j = 0; j < parts.size(); ++j )
      {
        if( j != i )
          val += (val.empty() ? "" : " ") + parts[j];
      }
    }//if( nucptr && !nucptr->isStable() )
  }//for( const string &nucstr : parts )
  
  // Try successive works, like "U 238", "60 Cobalt", etc
  //  TODO: not super robust or extensive (but we mostly expect single word input to specify nuclides)
  if( !nuc )
  {
    for( size_t i = 0; (i + 1) < parts.size(); ++i )
    {
      const string &first_word = parts[i];
      const string &second_word = parts[i+1];
      
      auto nucptr = db->nuclide( first_word + " " + second_word );
      if( nucptr && !nucptr->isStable() )
      {
        //Check for "meta", "meta-2", etc - not currently a very robust check
        string meta_str;
        for( const auto &n : parts )
        {
          if( SpecUtils::istarts_with(n, "meta") )
            meta_str = n;
        }
        
        if( !meta_str.empty() )
        {
          auto metaptr = db->nuclide( first_word + " " + second_word + " " + meta_str );
          if( metaptr )
            nucptr = metaptr;
          else
            meta_str.clear();
        }
        
        nuc = nucptr;
        
        // Reconstruct `val` with the parts we used removed
        val = "";
        for( size_t j = 0; j < parts.size(); ++j )
        {
          if( (j != i) && (j != (i+1)) && (meta_str.empty() || (j != (i+2))) )
            val += (val.empty() ? "" : " ") + parts[j];
        }
        
        // Not currently checking the user only inputed a single nuclide here - more hassle than worth
        break;
      }//
    }//for( size_t i = 0; (i + 1) < parts.size(); ++i )
  }//if( !nuc )
  
  SpecUtils::trim( val );
  
  return {nuc, val};
}//detect_nuclide(..)


/** Validate user input server-side, by just seeing if #UnitsConverterTool::convert works for it.
 
 Note: included further down is a WRegExpValidator implementation for units-only conversion, i.e., without a nuclide.
 */
class UnitsInputValidator : public Wt::WValidator
{
public:
  UnitsInputValidator( WObject *parent = nullptr )
    : WValidator( parent )
  {
    
  }
  
  virtual Wt::WValidator::Result validate( const WString &input ) const
  {
    string utfstr = input.toUTF8();
    if( utfstr.empty() )
      return Wt::WValidator::Result( Wt::WValidator::State::InvalidEmpty, "Empty" );
    
    try
    {
      const string result = UnitsConverterTool::convert( utfstr );
      if( result.empty() )//shouldnt happen, but JIC
        throw exception();
      
      return Wt::WValidator::Result( Wt::WValidator::State::Valid, WString::fromUTF8(result) );
    }catch(std::exception &e )
    {
      const string msg = e.what();
      return Wt::WValidator::Result( Wt::WValidator::State::Invalid, WString::fromUTF8(msg) );
    }
    
    return Wt::WValidator::Result( Wt::WValidator::State::Valid );
  }//validate()
};//UnitsInputValidator


#ifndef NDEBUG
void run_tests()
{
  // TODO: move these to the unit tests, and also probably move num_sig_figs to PhysicalUnits.
  
  assert( num_sig_figs("") == 3 );
  assert( num_sig_figs("asdas") == 3 );
  
  assert( num_sig_figs("0") == 2 );
  assert( num_sig_figs("01") == 2 );
  assert( num_sig_figs("00001") == 2 );
  assert( num_sig_figs("1") == 2 );
  assert( num_sig_figs("1.2") == 2 );
  assert( num_sig_figs("1.20") == 3 );
  assert( num_sig_figs("0.2000") == 4 );
  assert( num_sig_figs("1.2000") == 5 );
  assert( num_sig_figs(".2000") == 4 );
  assert( num_sig_figs(".20001") == 5 );
  assert( num_sig_figs("10.2000") == 6 );
  assert( num_sig_figs("010.2000") == 6 );
  assert( num_sig_figs("1.0E-6") == 2 );
  assert( num_sig_figs("1.00E-6") == 3 );
  assert( num_sig_figs("1.01E-6") == 3 );
  assert( num_sig_figs("1.010E-6") == 4 );
  assert( num_sig_figs("01.010E+007") == 4 );
  assert( num_sig_figs("1E+007") == 2 );
  
  assert( num_sig_figs("Co60 1.034 ug") == 4 );
  assert( num_sig_figs("Co60 1. uCi") == 2 );
  assert( num_sig_figs("Co60 1.01 ug") == 3 );
  
  
  auto test_string_to_mass = []( const string str, const double expected_value ){
    try
    {
      const double converted_value = PhysicalUnits::stringToMass( str );
      const double diff = fabs(converted_value - expected_value);
      if( diff > (1.0E-7*expected_value) ) //1.0E-7 is arbitrary
      {
        cerr << "Unexpected value converting '" << str << "' to grams - got " << converted_value
             << ", but expected " << expected_value << " (diff: " << diff << ")" << endl;
        assert( converted_value == expected_value );
      }
    }catch( std::exception &e )
    {
      cerr << "Unexpected failure converting '" << str
           << "' to mass (expected " << expected_value << ")" << endl;
      assert( 0 );
    }
  };//test_string_to_mass(...)
  
  auto test_string_to_mass_fails = []( const string str ){
    try
    {
      const double converted_value = PhysicalUnits::stringToMass( str );
      cerr << "stringToMass eroneously converted '" << str << "' to " << converted_value
           << " - it should not have." << endl;
      assert( 0 );
    }catch(...)
    {
    }
  };//auto test_string_to_mass_fails = []( const string str ){
  
  test_string_to_mass("1g", 1.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 gram", 1.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 kg", 1000.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 kilog", 1000.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 kilogram", 1000.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 kilo-gram", 1000.0*PhysicalUnits::gram );
  test_string_to_mass("1.0 oz", 28.3495*PhysicalUnits::gram );
  test_string_to_mass("0.2E1 lb", 2.0*453.592*PhysicalUnits::gram );
  test_string_to_mass("1.2 stone", 1.2*6350.29*PhysicalUnits::gram );
  test_string_to_mass("5 grain", 5.0*0.0647989*PhysicalUnits::gram );
  test_string_to_mass_fails( "lb" );
  test_string_to_mass_fails( "kg" );
  test_string_to_mass_fails( "1 mCi" );
  test_string_to_mass_fails( "1 ps" );
  test_string_to_mass_fails( "1 gr" );
  test_string_to_mass_fails( "1 gr" );
  test_string_to_mass_fails( "1" );
  test_string_to_mass_fails( "1E-3" );
  test_string_to_mass_fails( "1mm" );
}//void run_tests()
#endif //#ifndef NDEBUG

}//namespace


UnitsConverterTool::UnitsConverterTool()
  : AuxWindow( WString::tr("window-title-units-convert"),
               (Wt::WFlags<AuxWindowProperties>(AuxWindowProperties::PhoneNotFullScreen)
               | AuxWindowProperties::DisableCollapse
               | AuxWindowProperties::SetCloseable) ),
    m_input( NULL ),
    m_output( NULL ),
    m_message( NULL )
{
#ifndef NDEBUG
  run_tests();
#endif
  
  InterSpec *viewer = InterSpec::instance();
  viewer->useMessageResourceBundle( "UnitsConverterTool" );
  
  //addStyleClass( "UnitsConverterTool" );
  
  WGridLayout *layout = stretcher();
  setWidth(400);
  WText *message = new WText( WString::tr("uct-instructions"), Wt::XHTMLText );

  layout->addWidget(message,0,0,1,3);

  WLabel *label = new WLabel( WString::tr("uct-input-label") );
  layout->addWidget(label,1,0);
  //label->addStyleClass( "UnitsConverterToolLabel" );
  
  m_input = new WLineEdit();
  
  m_input->setAutoComplete( false );
  m_input->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_input->setAttributeValue( "autocorrect", "off" );
  m_input->setAttributeValue( "spellcheck", "off" );
#endif
  
  layout->addWidget(m_input,1,1);


  /*
  // To validate just unit conversions, we can use a regular expression
  //  pretty easily, but to validate the more complex input that inlcudes
  //  a nuclide, we have to do that server-side; but left in this comment
  //  is the client side units-only conversion; there is also a start to a
  //  regex to match nuclides, but its likely not worth considering.
  //
  //A regex to match a nuclide input like:
  //  Co60, cobalt-60, Co 60, Co-60, Co 60m, Co60m, Co60m2, Co60 meta, Co 60 meta
  //But not input like
  //  60-Co meta, 60Co, 60-Co
  //  (because then "5uC" would match a nuclide)
  //is:
  //const char *possible_nuc_regex = "(([a-zA-Z]{1,10})(\\s|-)*([0-9]{1,3})(\\s|-)*(m2|m3|m-2|m-3|meta|meta2|meta-2|m2|m){0,1})";
  
  const string regex = string("(") + ""
    "(" + PhysicalUnits::sm_absorbedDoseRegex + ")"
    "|(" + PhysicalUnits::sm_activityRegex + ")"
    "|(" + PhysicalUnits::sm_equivalentDoseRegex + ")"
    "|(" + PhysicalUnits::sm_distanceRegex + ")"
  ")";
  
  WRegExpValidator *validator = new WRegExpValidator( regex, this );
  validator->setFlags( Wt::MatchCaseInsensitive );
  m_input->setValidator(validator);
  */
  
  m_input->setValidator( new UnitsInputValidator(this) );
  
  
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
    
  label = new WLabel( WString::tr("uct-output-label") );
  layout->addWidget(label,2,0);
  //label->addStyleClass( "UnitsConverterToolLabel" );
  
  m_output = new WLineEdit( );
  
  m_output->setAttributeValue( "ondragstart", "return false" );
#if( BUILD_AS_OSX_APP || IOS )
  m_output->setAttributeValue( "autocorrect", "off" );
  m_output->setAttributeValue( "spellcheck", "off" );
#endif

  //layout->addWidget(m_output,2,1,1,2);
  layout->addWidget( m_output, 2, 1 );
  //m_output->addStyleClass( "UnitsConverterToolInput" );
  m_output->setTextSize( 15 );
  m_output->setWidth( WLength(15.0,WLength::FontEm) );
  m_output->setEnabled(false);
  //m_output->setText("135.14 uCi");
  
  m_message = new WText( "&nbsp", XHTMLUnsafeText );
  m_message->setHeight(WLength(50,WLength::Pixel));
  m_message->setAttributeValue( "style", "color: rgb(18,101,200);"  );
//  m_message->setHiddenKeepsGeometry( true );
  m_message->hide();
  layout->addWidget(m_message,3,0,1,3);
  layout->setColumnStretch(1, 1);
  layout->setRowStretch(3, 1);
    
  
#if( USE_QR_CODES )
  WPushButton *qr_btn = new WPushButton( footer() );
  qr_btn->setText( WString::tr("QR Code") );
  qr_btn->setIcon( "InterSpec_resources/images/qr-code.svg" );
  qr_btn->setStyleClass( "LinkBtn DownloadBtn DialogFooterQrBtn" );
  qr_btn->clicked().preventPropagation();
  qr_btn->clicked().connect( std::bind( [this](){
    try
    {
      const string url = "interspec://unit/?" + Wt::Utils::urlEncode(encodeStateToUrl());
      QrCode::displayTxtAsQrCode( url, WString::tr("uct-qr-window-title"),
                                 WString::tr("uct-qr-window-txt") );
    }catch( std::exception &e )
    {
      passMessage( WString::tr("app-qr-err").arg(e.what()), WarningWidget::WarningMsgHigh );
    }
  }) );
#endif //USE_QR_CODES
  
  WPushButton *closeButton = addCloseButtonToFooter();
  closeButton->clicked().connect( boost::bind( &AuxWindow::hide, this ) );

  rejectWhenEscapePressed();
  
  convert();
  
  show();
  centerWindow();
  resizeToFitOnScreen();

  //Keep the keyboard form popping up
  if( viewer->isMobile() )
  {
    if( viewer->isPhone() )
    {
      int w = viewer->renderedWidth();
      int h = viewer->renderedHeight();
      if( w < 100 )
      {
        w = wApp->environment().screenWidth();
        h = wApp->environment().screenHeight();
      }
      
      if(  (w > 100) && (w > h) )
        titleBar()->hide();
    }//if( viewer->isPhone() )
    
    closeButton->setFocus();
  }//if( viewer->isMobile() )
  
}//UnitsConverterTool constructor


UnitsConverterTool::~UnitsConverterTool()
{
  //nothing to do here
}//UnitsConverterTool destructor


std::string UnitsConverterTool::convert( std::string val )
{
  SpecUtils::trim( val );
  
  if( val.empty() )
    throw runtime_error( WString::tr("uct-err-empty-input").toUTF8() );
  
  pair<const SandiaDecay::Nuclide *, std::string> nuc_res = detect_nuclide( val );
  const SandiaDecay::Nuclide *nuc = nuc_res.first;
  	
  
  if( nuc )
  {
    val = nuc_res.second;
    
    // TODO: could specialize things to detect dose or activity, as well as distance to return activity or dose...
    //       (but for the moment we'll just convert activity <---> masses
    
    // Detect activity, to return mass
    try
    {
      const double dbvalue = PhysicalUnits::stringToActivity(val);
      const double grams = dbvalue / nuc->activityPerGram();
      
      return PhysicalUnits::printToBestMassUnits(grams * PhysicalUnits::gram, num_sig_figs(val) + 1 );
    }catch(...)
    {
    }
    
    // Detect mass, to return activity
    try
    {
      const double grams = PhysicalUnits::stringToMass( val ) / PhysicalUnits::gram;
      const double actPerGram = PhysicalUnits::bq * nuc->activityPerGram() / SandiaDecay::Bq;
      const double activity = actPerGram * grams;
      
      bool useCuries = false;
      InterSpec *viewer = InterSpec::instance();
      if( viewer)
        useCuries = !InterSpecUser::preferenceValue<bool>( "DisplayBecquerel", viewer );
      
      return PhysicalUnits::printToBestActivityUnits( activity, num_sig_figs(val), useCuries );
    }catch(...)
    {
    }
    
    if( val.empty() )
    {
      // Return some basic info, half life, specificActivity, anything else
    }
    
    throw runtime_error( WString::tr("uct-nuc-unknown").arg(nuc->symbol).arg(val).toUTF8() );
  }//if( nuc )
  
  
  string ans = convertActivity( val );
  
  if( ans.empty() )
    ans = convertAbsorbedDose( val );
  
  if( ans.empty() )
    ans = convertEquivalentDose( val );
  
  if( ans.empty() )
    ans = convertDistance( val );
  
  if( ans.empty() )
    ans = convertMass( val );
  
  if( ans.empty() )
  {
    boost::smatch mtch;
    boost::regex expr( PhysicalUnits::sm_positiveDecimalRegex );
    if( !boost::regex_match( val, mtch, expr ) )
      throw runtime_error( WString::tr("uct-couldnt-find-units").toUTF8() );
    throw runtime_error( WString::tr("uct-couldnt-find-units").toUTF8() );
  }//if( ans.empty() )
  
  return ans;
}//static std::string convert( std::string input );


void UnitsConverterTool::convert()
{
  try
  {
    switch( m_input->validate() )
    {
      case Wt::WValidator::Invalid:
        throw std::runtime_error( WString::tr("uct-invalid-input").toUTF8() );
        
      case Wt::WValidator::InvalidEmpty:
        throw runtime_error( "" );
        
      case Wt::WValidator::Valid:
        break;
    }//switch( m_input->validate() )
    
    // TODO: make this function static, so we can call from terminal widget, and then add in converting mass, and detecting if a nulcide was specified, and if so, if alone print out summary, but if with another quantity convert for that
      
    WString current = m_input->text();
    std::string val = current.toUTF8();
    std::string ans = convert( val );
    
    m_output->setText( Wt::WString::fromUTF8(ans) );
    
    //If weve made it this far, we have valid input
    m_message->removeStyleClass("line-above");
    m_message->setText( "&nbsp" );
    m_message->hide();
    
    // We will only add an undo/redo step if the input validated, and both the input
    //  and output have changed (i.e., we convert as the user types, so we dont want
    //  the users typing in "5 meters" to be like three seperate steps).
    if( (m_input->text() != m_prevInput) && (ans != m_prevAnswer) )
    {
      UndoRedoManager *undoRedo = UndoRedoManager::instance();
      if( undoRedo && undoRedo->canAddUndoRedoNow() )
      {
        const WString prev = m_prevInput;
        const string prev_ans = m_prevAnswer;
        auto undo = [prev,prev_ans](){
          UnitsConverterTool *tool = InterSpec::instance()->createUnitsConverterTool();
          if( !tool )
            return;
          tool->m_prevInput = prev;
          tool->m_prevAnswer = prev_ans;
          tool->m_input->setText( prev );
          tool->convert();
        };
        
        auto redo = [current,ans](){
          UnitsConverterTool *tool = InterSpec::instance()->createUnitsConverterTool();
          if( !tool )
            return;
          tool->m_prevInput = current;
          tool->m_prevAnswer = ans;
          tool->m_input->setText( current );
          tool->convert();
        };
        
        undoRedo->addUndoRedoStep( std::move(undo), std::move(redo), "Update units converter value." );
      }//if( undoRedo && undoRedo->canAddUndoRedoNow() )
    }//if( m_input->text() != m_prevInput )
    
    m_prevInput = std::move(current);
    m_prevAnswer = std::move(ans);
  }catch( std::exception &e )
  {
    m_output->setText( "" );
    m_message->setText( WString::tr("uct-err-gen").arg(e.what()) );
    //m_message->addStyleClass("line-above");
    m_message->show();
  }//try / catch
}//void UnitsConverterTool::convert()


void UnitsConverterTool::handleAppUrl( std::string query_str )
{
  UndoRedoManager::BlockUndoRedoInserts undo_blocker;
  
  const map<string,string> values = AppUtils::query_str_key_values( query_str );
  
  // Check version is appropriate
  const auto ver_iter = values.find( "V" );
  if( (ver_iter == end(values))
     || ((ver_iter->second != "1") && !SpecUtils::istarts_with(ver_iter->second, "1.")) )
    throw runtime_error( "Units Converter: URI not compatible version." );
  
  const auto iter = values.find( "INPUT" );
  if( iter == end(values) )
    throw runtime_error( "Units Converter: no input value specified." );
  
  string value = iter->second;
  SpecUtils::ireplace_all( value, "%23", "#" );
  SpecUtils::ireplace_all( value, "%26", "&" );
  m_prevInput = WString::fromUTF8(value);
  m_input->setValueText( m_prevInput );
  
  convert();
}//void handleAppUrl( std::string query_str )


std::string UnitsConverterTool::encodeStateToUrl() const
{
  string val = m_input->text().toUTF8();
  SpecUtils::ireplace_all( val, "#", "%23" );
  SpecUtils::ireplace_all( val, "&", "%26" );
  
  return "V=1&INPUT=" + val;
}//std::string encodeStateToUrl() const
