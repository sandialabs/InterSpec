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
#include <cctype>
#include <exception>

#include <Wt/WColor>
#include <Wt/WFormWidget>
#include <Wt/WApplication>
#include <Wt/WEnvironment>

#include "InterSpec/ColorSelect.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


namespace
{
  bool hasNativePicker()
  {
#if( BUILD_AS_OSX_APP || IOS )
    return false;
#else
    const WEnvironment &env = wApp->environment();
    const bool doesntHave = (   env.agentIsMobileWebKit()
                        || env.agentIsIEMobile()
                        || env.agentIsMobileWebKit()
                        || env.agent()==WEnvironment::MobileWebKitiPhone
                        || env.agent()==WEnvironment::MobileWebKitAndroid
                        || env.agent()==WEnvironment::MobileWebKit
                        || env.userAgent().find("Opera Mobi") != std::string::npos
                        || env.userAgent().find("Android") != std::string::npos
                        || env.userAgent().find("RIM ") != std::string::npos
                        || env.userAgent().find("iPad") != std::string::npos
                        || env.userAgent().find("iPhone") != std::string::npos
                        || env.userAgent().find("iPod") != std::string::npos
                        || env.userAgent().find("Mobile") != std::string::npos
                        || wApp->environment().agentIsIE()
                        );
    return !doesntHave;
#endif
  }//bool hasNativePicker() const
  
  std::string toHexFormat( const WColor &c )
  {
    if( c.isDefault() )
      return "#000000";
    
    //The clamping probably overkill
    const int r = std::min( std::max(c.red(),0), 255 );
    const int g = std::min( std::max(c.green(),0), 255 );
    const int b = std::min( std::max(c.blue(),0), 255 );
    
    //Note: we are (silently) not supporting alpha right now; a limitation of the
    //      HTML5 color picker atm.
    //const int a = std::min( std::max(c.alpha(),0), 255 );
    
    //Put the RGB in hex format like "#F90800"
    char buffer[16];
    snprintf( buffer, sizeof(buffer), "#%02X%02X%02X", r, g, b );
    
    return buffer;
  }
  
#if( defined(_MSC_VER) || ANDROID )
  bool ishexnumber(const char c)
  {
	  static const char hexchars[] = { '0','1','2','3','4','5','6','7','8','9','a','A','b','B','c','C','d','D','e','E','f','F' };
	  for (const char h : hexchars)
		  if (c == h)
			  return true;
	  return false;
  }
#endif //#if( defined(_MSC_VER) )

  /** The value must be in seven-character hexadecimal notation, meaning the "#"
   character followed by two digits each representing red, green, and blue,
   like this: "#rrggbb"
   
   Throws exception if format is not valid.
   */
  void checkColorStringFormat( const std::string &val )
  {
    if( val.empty() )
      return;
    
    if( val.size() != 7 )
      throw runtime_error( "Css Color must be 7 characters long" );
    
    if( val[0] != '#' )
      throw runtime_error( "Css Color must have first character of '#'" );
    
    for( size_t i = 1; i < val.size(); ++i )
      if( !ishexnumber(val[i]) )
        throw runtime_error( "For Css color, all characters after first must be hex digits" );
  }//checkColorStringFormat(...)

  
  const std::string ns_init_js = INLINE_JAVASCRIPT
  (
    function(id,colorstr)
    {
      var el = $('#'+id);
      
      var onChange = function(event){
        try{
          if( el.data('color') === el.get(0).value )
          {
            console.log( "No change in color, returning" );
            return;
          }else
          {
            el.data('color',el.get(0).value);
            Wt.emit( id, {name: 'colorchange', eventObject: event}, el.get(0).value );
          }
        }catch(err){
          console.log( "Caught error in onChange" );
        }
      };
      
      try
      {
        el.data('color',colorstr);
        el.get(0).addEventListener("input", onChange, false);  //is fired on the <input> element every time the color changes.
        //colorWell.addEventListener("change", updateAll, false);  //The change event is fired when the user dismisses the color picker.
      }catch(err){
        console.log( "Got error initing ColorSelect in JS" );
      }
    }
   );
}//namespace


ColorSelect::ColorSelect( Wt::WFlags<ColorSelectOptions> options, Wt::WContainerWidget *parent )
: WFormWidget( parent ),
  m_hasBeenSet( false ),
  m_usingNative( options.testFlag(PrefferNative) && hasNativePicker() ),
  m_userSelectedColor( this, "colorchange", true )
{
  if( m_usingNative )
  {
    setAttributeValue( "type", "color" );
    setAttributeValue( "value", "#000000" );
    
    doJavaScript( "var fcn = " + ns_init_js + "; fcn('" + id() + "','#000000');" );
  }else
  {
    wApp->require( "InterSpec_resources/assets/js/spectrum_1.8.0/spectrum.min.js" );
    wApp->useStyleSheet( "InterSpec_resources/assets/js/spectrum_1.8.0/spectrum.min.css" );
    
    string init_js = "$('#" + id() + "').spectrum({"
      "showSelectionPalette: true, "
      //"showInitial: true, "
      //"maxSelectionSize: 10,"
      "localStorageKey: \"colorpicker.palette\","
      "change: function(color){ Wt.emit('" + id() + "', {name: 'colorchange'}, (color ? color.toHexString() : 'null') );  }";
    if( options.testFlag(AllowNoColor) )
    {
      init_js += ",color: null, allowEmpty: true";
      setAttributeValue( "value", "" );
    }else
    {
      init_js += ",color: \"#000000\"";
      setAttributeValue( "value", "#000000" );
    }
    
    // className: "full-spectrum",
    
    //ToDo: Support alpha - requires patching Wt.
    //init_js += ",showAlpha: true";
    
    if( options.testFlag(DontShowPalete) )
    {
      init_js += ",showPalette: false";
    }else
    {
      //ToDo: come up with a way to pass in the colors in a color theme, and put
      //      those in the palete first, followed by the following.
      init_js += ",showPalette: true";
      init_js += R"tok( ,palette: [
      ["rgb(0, 0, 0)",       "rgb(67, 67, 67)",      "rgb(102, 102, 102)",   "rgb(138, 138, 138)",  "rgb(204, 204, 204)",   "rgb(217, 217, 217)",   "rgb(255, 255, 255)"],
      ["rgb(152, 0, 0)",     "rgb(255, 0, 0)",       "rgb(255, 153, 0)",     "rgb(255, 255, 0)",    "rgb(0, 255, 0)",
       "rgb(0, 255, 255)",   "rgb(74, 134, 232)",    "rgb(0, 0, 255)",       "rgb(153, 0, 255)",    "rgb(255, 0, 255)"],
      ["rgb(230, 184, 175)", "rgb(244, 204, 204)",   "rgb(252, 229, 205)",   "rgb(255, 242, 204)",  "rgb(217, 234, 211)",
       "rgb(208, 224, 227)", "rgb(201, 218, 248)",   "rgb(207, 226, 243)",   "rgb(217, 210, 233)",  "rgb(234, 209, 220)"],
      ["rgb(221, 126, 107)", "rgb(234, 153, 153)",   "rgb(249, 203, 156)",   "rgb(255, 229, 153)",  "rgb(182, 215, 168)",
       "rgb(162, 196, 201)", "rgb(164, 194, 244)",   "rgb(159, 197, 232)",   "rgb(180, 167, 214)",  "rgb(213, 166, 189)"],
      ["rgb(204, 65, 37)",   "rgb(224, 102, 102)",   "rgb(246, 178, 107)",   "rgb(255, 217, 102)",  "rgb(147, 196, 125)",
       "rgb(118, 165, 175)", "rgb(109, 158, 235)",   "rgb(111, 168, 220)",   "rgb(142, 124, 195)",  "rgb(194, 123, 160)"],
      ["rgb(166, 28, 0)",    "rgb(204, 0, 0)",       "rgb(230, 145, 56)",    "rgb(241, 194, 50)",   "rgb(106, 168, 79)",
       "rgb(69, 129, 142)",  "rgb(60, 120, 216)",    "rgb(61, 133, 198)",    "rgb(103, 78, 167)",   "rgb(166, 77, 121)"],
      ["rgb(91, 15, 0)",     "rgb(102, 0, 0)",       "rgb(120, 63, 4)",      "rgb(127, 96, 0)",     "rgb(39, 78, 19)",
       "rgb(12, 52, 61)",    "rgb(28, 69, 135)",     "rgb(7, 55, 99)",       "rgb(32, 18, 77)",     "rgb(76, 17, 48)"]] )tok";
    }
    
    if( options.testFlag(DontShowTextInput) )
      init_js += ",showInput: false";
    else
      init_js += ",showInput: true, preferredFormat: \"rgb\"";  //hex other option
    
    init_js += "});";
    
    doJavaScript( init_js );
  }//if( m_usingNative ) / else
  
  m_userSelectedColor.connect( boost::bind( &ColorSelect::colorSetCallback, this, _1 ) );
}//ColorSelect(...)


ColorSelect::~ColorSelect()
{
}


Wt::WColor ColorSelect::color() const
{
  const std::string c = attributeValue("value").toUTF8();
  if( !m_hasBeenSet || c.empty() )
    return WColor();
  return WColor(c);
}//WColor color() const


void ColorSelect::setColor( const Wt::WColor &c )
{
  m_hasBeenSet = !c.isDefault();
  if( m_hasBeenSet )
    setValueText( toHexFormat(c) );
  else
    setValueText( "" );
}//void setColor( const Wt::WColor &c ) const;


void ColorSelect::setDisabled( bool disabled )
{
  Wt::WFormWidget::setDisabled(disabled);
  if( !m_usingNative )
  {
    if( disabled )
      doJavaScript("$('#" + id() + "').spectrum(\"disable\");" );
    else
      doJavaScript("$('#" + id() + "').spectrum(\"enable\");" );
  }
}

void ColorSelect::colorSetCallback( const std::string &val )
{
  const string old = attributeValue("value").toUTF8();
  
  if( val == "null" )
  {
    m_hasBeenSet = true;
    setAttributeValue( "value", "" );
    m_cssColorChanged.emit( WColor() );
    return;
  }
  
  try
  {
    WColor newcolor( val );
    m_hasBeenSet = true;
    if( old != val )
    {
      setAttributeValue( "value", val );
      m_cssColorChanged.emit( newcolor );
    }
  }catch( std::exception & )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    log_developer_error( BOOST_CURRENT_FUNCTION, "Invalid CSS color called back " );
#endif
    setAttributeValue( "value", "#000000" );  //to force a repaint ...
    setAttributeValue( "value", old );
    if( m_usingNative )
    {
      doJavaScript( "try{" + jsRef() + ".value='" + old + "';}catch(e){}" );  //just to make sure
    }else
    {
      doJavaScript("$('#" + id() + "').spectrum(\"set\",\"" + old + "\");" );
    }
  }//try / catch
}//void colorSetCallback( const std::string &val )


WT_USTRING ColorSelect::valueText() const
{
  if( !m_hasBeenSet )
    return "";
  return attributeValue("value");
}

void ColorSelect::setValueText( const WT_USTRING &value )
{
  const string color = value.toUTF8();
  checkColorStringFormat( color );

  const string old = attributeValue("value").toUTF8();
  if( !m_hasBeenSet || old != value )
  {
    m_hasBeenSet = !value.empty();
    setAttributeValue( "value", value );

    if( m_usingNative )
    {
      //For some reason we need this next line, or else when the color of reference
      //  lines are changed, it wont show up in this widget... not sure, weird!
      doJavaScript( "try{" + jsRef() + ".value='" + value.toUTF8() + "';}catch(e){}" );
    }else
    {
      doJavaScript("$('#" + id() + "').spectrum(\"set\"," + (value.empty() ? string("null") : ("\"" + value.toUTF8() + "\"")) + ");" );
    }
  }//if( old != value )
}//setValueText(...)


Wt::Signal<Wt::WColor> &ColorSelect::cssColorChanged()
{
  return m_cssColorChanged;
}


DomElementType ColorSelect::domElementType() const
{
  return Wt::DomElement_INPUT;
}

bool ColorSelect::willUseNativeColorPicker()
{
  return false;
  //return hasNativePicker();
}

