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
#include <memory>

#include <Wt/WColor.h>
#include <Wt/WFormWidget.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include "InterSpec/InterSpec.h"
#include "InterSpec/ColorTheme.h"
#include "InterSpec/ColorSelect.h"

using namespace std;
using namespace Wt;

#define INLINE_JAVASCRIPT(...) #__VA_ARGS__


namespace
{
  /** All modern browsers (Chrome, Firefox, Safari 15.4+, Edge) support
      the HTML5 color input.  However, the native picker does not support
      custom swatches, "no color" option, or other features that Coloris provides.
      On desktop builds, we currently always use Coloris instead.
   */
  bool hasNativePicker()
  {
#if( BUILD_AS_OSX_APP || IOS )
    return false;
#else
    // All browsers we support have the native color picker, but we return
    // false to always use Coloris for consistent UX and feature support.
    return false;
#endif
  }//bool hasNativePicker()

  /** Returns the CSS text representation of a color (e.g., "rgb(255, 0, 0)").
      Returns "#000000" for default/invalid colors.
   */
  std::string toCssText( const WColor &c )
  {
    if( c.isDefault() )
      return "#000000";
    return c.cssText( false );
  }


  // JS for native HTML5 color picker initialization
  const std::string ns_init_js = INLINE_JAVASCRIPT
  (
    function(id,colorstr)
    {
      var el = document.getElementById(id);
      if( !el ) return;
      el._isData = el._isData || {};

      var onChange = function(event){
        try{
          if( el._isData.color === el.value )
          {
            console.log( "No change in color, returning" );
            return;
          }else
          {
            el._isData.color = el.value;
            Wt.emit( id, {name: 'colorchange', eventObject: event}, el.value );
          }
        }catch(err){
          console.log( "Caught error in onChange" );
        }
      };

      try
      {
        el._isData.color = colorstr;
        el.addEventListener("input", onChange, false);
      }catch(err){
        console.log( "Got error initing ColorSelect in JS" );
      }
    }
   );


  // JS to set up a ColorSelect input before Coloris is loaded.
  // Sets up click handler to request lazy-load, and change handler for color updates.
  const std::string ns_coloris_pre_init_js = INLINE_JAVASCRIPT
  (
    function( id, loadSig, changeSig )
    {
      var el = document.getElementById( id );
      if( !el ) return;

      el.readOnly = true;
      el.style.cursor = 'pointer';
      el.style.backgroundColor = el.value || 'transparent';
      el.style.border = '0';
      el.style.boxShadow = 'inset 0 0 1px rgba(0,0,0,.5)';
      el.style.borderRadius = '5px';
      el.style.width = '30px';
      el.style.height = '22px';
      el.style.padding = '0';
      el.style.caretColor = 'transparent';
      el.style.color = 'transparent';

      el.addEventListener( 'click', function( e )
      {
        if( typeof Coloris === 'undefined' )
        {
          e.preventDefault();
          e.stopPropagation();
          Wt.emit( id, loadSig );
        }else
        {
          var picker = document.getElementById( 'clr-picker' );
          if( picker )
          {
            var maxZ = 0;
            var node = el;
            while( node && node !== document.body )
            {
              var z = parseInt( window.getComputedStyle( node ).zIndex, 10 );
              if( !isNaN(z) && z > maxZ )
                maxZ = z;
              node = node.parentElement;
            }
            picker.style.zIndex = ( maxZ > 0 ? maxZ + 1 : 10000 );
          }
        }
      }, true );

      el.addEventListener( 'change', function()
      {
        el.style.backgroundColor = el.value || 'transparent';
        Wt.emit( id, changeSig, el.value || 'null' );
      });
    }
  );


  // Builds the JS expression that returns the swatches array, including user palette from localStorage
  std::string buildSwatchesJs()
  {
    // The fixed palette converted from the original spectrum.js rgb() values to hex
    return
      "(function(){"
        "var user=[];"
        "try{user=JSON.parse(localStorage.getItem('colorpicker.palette')||'[]');}catch(e){}"
        "var fixed=["
          "'#000000','#434343','#666666','#8A8A8A','#CCCCCC','#D9D9D9','#FFFFFF',"
          "'#980000','#FF0000','#FF9900','#FFFF00','#00FF00',"
          "'#00FFFF','#4A86E8','#0000FF','#9900FF','#FF00FF',"
          "'#E6B8AF','#F4CCCC','#FCE5CD','#FFF2CC','#D9EAD3',"
          "'#D0E0E3','#C9DAF8','#CFE2F3','#D9D2E9','#EAD1DC',"
          "'#DD7E6B','#EA9999','#F9CB9C','#FFE599','#B6D7A8',"
          "'#A2C4C9','#A4C2F4','#9FC5E8','#B4A7D6','#D5A6BD',"
          "'#CC4125','#E06666','#F6B26B','#FFD966','#93C47D',"
          "'#76A5AF','#6D9EEB','#6FA8DC','#8E7CC3','#C27BA0',"
          "'#A61C00','#CC0000','#E69138','#F1C232','#6AA84F',"
          "'#45818E','#3C78D8','#3D85C6','#674EA7','#A64D79',"
          "'#5B0F00','#660000','#783F04','#7F6000','#274E13',"
          "'#0C343D','#1C4587','#073763','#20124D','#4C1130'"
        "];"
        "return user.concat(fixed.filter(function(c){return user.indexOf(c.toUpperCase())<0&&user.indexOf(c.toLowerCase())<0;}));"
      "}())";
  }
}//namespace


ColorSelect::ColorSelect( Wt::WFlags<ColorSelectOptions> options )
: WFormWidget(),
  m_hasBeenSet( false ),
  m_usingNative( options.test(PrefferNative) && hasNativePicker() ),
  m_allowNoColor( options.test(AllowNoColor) ),
  m_userSelectedColor( this, "colorchange", true ),
  m_requestColorisLoad( this, "loadcoloris", true ),
  m_cssColorChanged()
{
  if( m_usingNative )
  {
    setAttributeValue( "type", "color" );
    setAttributeValue( "value", "#000000" );

    doJavaScript( "var fcn = " + ns_init_js + "; fcn('" + id() + "','#000000');" );
  }else
  {
    setAttributeValue( "type", "text" );
    setAttributeValue( "data-coloris", "" );

    if( m_allowNoColor )
    {
      setAttributeValue( "value", "" );
    }else
    {
      setAttributeValue( "value", "#000000" );
    }

    // Set up click handler (for lazy-load or z-index) and change handler
    doJavaScript( "var fcn = " + ns_coloris_pre_init_js + ";"
      " fcn('" + id() + "', 'loadcoloris', 'colorchange');" );

    const bool colorisAlreadyLoaded = wApp->styleSheet().isDefined( "ColorisPickerStyle" );
    if( colorisAlreadyLoaded )
    {
      // Coloris is already loaded for this session - directly initialize this widget
      initColorisForWidget();
    }else
    {
      // Coloris not yet loaded - lazy-load on first click via server round-trip
      m_requestColorisLoad.connect( this, &ColorSelect::loadColoris );
    }
  }//if( m_usingNative ) / else

  m_userSelectedColor.connect( this, [this]( const std::string &val ){ colorSetCallback( val ); } );
}//ColorSelect(...)


ColorSelect::~ColorSelect()
{
  if( !m_usingNative && m_allowNoColor )
  {
    doJavaScript( "try{Coloris.removeInstance('#" + id() + "');}catch(e){}" );
  }
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
    setValueText( toCssText(c) );
  else
    setValueText( "" );
}//void setColor( const Wt::WColor &c ) const;


void ColorSelect::setDisabled( bool disabled )
{
  Wt::WFormWidget::setDisabled( disabled );
  if( !m_usingNative )
  {
    doJavaScript( "try{var el=" + jsRef() + ";"
      "if(el) el.disabled=" + (disabled ? "true" : "false") + ";"
      "}catch(e){}" );
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
    log_developer_error( __func__, "Invalid CSS color called back " );
#endif
    setAttributeValue( "value", "#000000" );  //to force a repaint ...
    setAttributeValue( "value", old );
    if( m_usingNative )
    {
      doJavaScript( "try{" + jsRef() + ".value='" + old + "';}catch(e){}" );
    }else
    {
      // Set the input value and update the Coloris preview wrapper if present
      doJavaScript( "try{var el=" + jsRef() + "; if(el){"
        "el.value='" + old + "';"
        "el.style.backgroundColor='" + old + "';"
        "var w=el.parentNode;"
        "if(w&&w.classList.contains('clr-field')) w.style.color='" + old + "';"
        "}}catch(e){}" );
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

  const string old = attributeValue("value").toUTF8();
  if( !m_hasBeenSet || old != value )
  {
    m_hasBeenSet = !value.empty();
    setAttributeValue( "value", value );

    if( m_usingNative )
    {
      doJavaScript( "try{" + jsRef() + ".value='" + value.toUTF8() + "';}catch(e){}" );
    }else
    {
      const string colorVal = value.empty() ? string("transparent") : value.toUTF8();
      doJavaScript( "try{var el=" + jsRef() + "; if(el){"
        "el.value='" + value.toUTF8() + "';"
        "el.style.backgroundColor='" + colorVal + "';"
        "var w=el.parentNode;"
        "if(w&&w.classList.contains('clr-field')) w.style.color='" + colorVal + "';"
        "}}catch(e){}" );
    }
  }//if( old != value )
}//setValueText(...)


Wt::Signal<Wt::WColor> &ColorSelect::cssColorChanged()
{
  return m_cssColorChanged;
}


DomElementType ColorSelect::domElementType() const
{
  return Wt::DomElementType::INPUT;
}

bool ColorSelect::willUseNativeColorPicker()
{
  return false;
  //return hasNativePicker();
}


void ColorSelect::loadColoris()
{
  // Note: Coloris uses delegated event handlers on `document` that call
  // Element.prototype.matches.call(e.target, selector).  Wt's event handling
  // can cause e.target to not be a proper Element, resulting in a couple
  // "Illegal invocation" console errors during the initial lazy-load click.
  // These are non-fatal and do not affect picker functionality (we bind
  // directly to elements, avoiding the delegation path for our inputs).
  //
  // The proper fix is to prepend a try/catch wrapper around matches() at the
  // top of coloris.min.js, before Coloris's IIFE captures the reference:
  //   !function(){if(!Element.prototype._origMatches){
  //     var o=Element.prototype.matches;
  //     Element.prototype._origMatches=o;
  //     Element.prototype.matches=function(s){
  //       try{return o.call(this,s)}catch(e){return false}
  //     };
  //   }}();
  // This is not done here to avoid modifying the vendored Coloris library and
  // to avoid the global performance impact of wrapping a native DOM method.

  wApp->require( "InterSpec_resources/assets/js/Coloris-0.25.0/coloris.min.js" );
  wApp->useStyleSheet( "InterSpec_resources/assets/js/Coloris-0.25.0/coloris.min.css" );

  // Determine dark/light mode from the current color theme
  string themeMode = "light";
  InterSpec *const interspec = InterSpec::instance();
  if( interspec )
  {
    const std::shared_ptr<const ColorTheme> theme = interspec->getColorTheme();
    if( theme && (theme->nonChartAreaTheme == "dark") )
      themeMode = "dark";
  }

  // Add CSS to reduce picker size/padding
  WCssStyleSheet &style = wApp->styleSheet();
  const string rulename = "ColorisPickerStyle";
  if( !style.isDefined(rulename) )
  {
    style.addRule( ".clr-picker",
      "width: 220px !important;"
      "border-radius: 6px !important;",
      rulename );

    style.addRule( ".clr-picker .clr-gradient",
      "height: 80px !important;"
      "margin: 0 0 6px !important;"
      "border-radius: 6px 6px 0 0 !important;" );

    style.addRule( ".clr-picker .clr-hue",
      "margin: 2px 8px !important;" );

    style.addRule( ".clr-picker .clr-alpha",
      "margin: 2px 8px !important;" );

    style.addRule( ".clr-picker .clr-swatches",
      "margin: 2px 4px !important;" );

    style.addRule( ".clr-picker .clr-swatches div",
      "justify-content: center !important;" );

    style.addRule( ".clr-picker .clr-swatches button",
      "width: 16px !important;"
      "height: 16px !important;" );

    style.addRule( ".clr-picker .clr-color",
      "margin: 6px 8px 6px 4px !important;" );

    style.addRule( ".clr-picker .clr-preview",
      "margin: 6px 0 6px 8px !important;" );

    style.addRule( ".clr-picker .clr-close, .clr-picker .clr-clear",
      "margin: 0 8px 4px 8px !important;" );

    style.addRule( ".clr-picker .clr-close",
      "margin-left: auto !important;" );

    style.addRule( ".clr-picker .clr-format",
      "margin: 0 8px 4px !important;" );
  }

  // Build global Coloris configuration.
  // We pass an array of actual DOM elements instead of a CSS selector because
  // Wt's event handling wraps events in a way that breaks Coloris's delegated
  // event handler (Element.prototype.matches fails with "Illegal invocation").
  // Passing elements directly binds click/input handlers to each element,
  // avoiding the delegation path entirely.
  string init_js =
    "var allColorInputs=document.querySelectorAll('[data-coloris]');"
    "Coloris({"
      "el:Array.from(allColorInputs),"
      "theme:'default',"
      "themeMode:'" + themeMode + "',"
      "wrap:true,"
      "alpha:false,"
      "format:'rgb',"
      "formatToggle:true,"
      "closeButton:true,"
      "swatches:" + buildSwatchesJs() + ","
      "onChange:function(color,el){"
        "if(color){"
          "try{"
            "var p=JSON.parse(localStorage.getItem('colorpicker.palette')||'[]');"
            "p=p.filter(function(c){return c.toUpperCase()!==color.toUpperCase();});"
            "p.unshift(color);"
            "if(p.length>10) p.length=10;"
            "localStorage.setItem('colorpicker.palette',JSON.stringify(p));"
          "}catch(e){}"
        "}"
      "}"
    "});";

  // Add Wt-btn class to Coloris Close/Clear buttons so they match app styling
  init_js +=
    "var cp=document.getElementById('clr-picker');"
    "if(cp){"
      "var btns=cp.querySelectorAll('.clr-close,.clr-clear');"
      "btns.forEach(function(b){b.classList.add('ColorisBtn');});"
    "}";

  // Per-instance config for AllowNoColor widgets
  if( m_allowNoColor )
    init_js += "Coloris.setInstance('#" + id() + "',{clearButton:true});";

  // Re-trigger click to open the picker now that Coloris is loaded
  init_js += "setTimeout(function(){"
    "var el=" + jsRef() + ";"
    "if(el) el.click();"
  "},50);";

  doJavaScript( init_js );
}//void loadColoris()


void ColorSelect::initColorisForWidget()
{
  // Coloris is already loaded globally - just bind and wrap this specific element.
  // Coloris.set({el: [element]}) calls D(element) to bind click/input handlers,
  // and R(element) to wrap the input in a clr-field div.
  string init_js = "Coloris.set({el:[" + jsRef() + "]});";

  if( m_allowNoColor )
    init_js += "Coloris.setInstance('#" + id() + "',{clearButton:true});";

  doJavaScript( init_js );
}//void initColorisForWidget()
