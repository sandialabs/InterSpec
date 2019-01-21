#ifndef ColorSelect_h
#define ColorSelect_h
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

#include <Wt/WColor>
#include <Wt/WFormWidget>

/** Simple class to allow the user to select a color using a native HTML5 color
    choosing input.
 
    ToDo: check out what happens on old browsers, and maybe fallback to
          something like http://bgrins.github.io/spectrum/ or
          https://github.com/wieringen/tinycolorpicker
          Also, currently alpha isnt supported.
 */
class ColorSelect : public Wt::WFormWidget
{
public:
  enum ColorSelectOptions
  {
    /** If selected, a native HTML5 color chooser will be used, if the client
     supports it; if not the JavaScript chooser ('Spectrum',
     https://bgrins.github.io/spectrum/) will be used.
     */
    PrefferNative = 0x01,
    
    /** Allows user to select no color at all.  This option is not supported by
     native color choosers, so will be ignored in those cases.
     */
    AllowNoColor = 0x02,
    
    /** Do not show a palete of colors user can pick from
     Not supported for native color picker (it will always give a palete option)
     */
    DontShowPalete = 0x04,
    
    /** Do not show text input filed user can enter text color.
     Not supported for native color picker (it will always give a palete option)
     */
    DontShowTextInput = 0x08,
  };
  
  /** Constructor: If a native color picker is used, or AllowNoColor is not
   specified, color defaults to #000000 (black).  If non-native and AllowNoColor
   is specified, defaults to no color.
   */
  ColorSelect( Wt::WFlags<ColorSelectOptions> options, Wt::WContainerWidget *parent = nullptr );
  
  virtual ~ColorSelect();
  
  /** */
  Wt::WColor color() const;
  
  /** */
  void setColor( const Wt::WColor &c );
  
  /** */
  virtual void setDisabled( bool disabled );
  
  /** Signal emmitted whenever the user selects a new color.  Argument should
      always be a valid css string that has seven characters, and analogous to
      '#000000', where the zeroes could be any valid hexadecimal digit.
   */
  Wt::Signal<Wt::WColor> &cssColorChanged();
  
  /** The current color. Ex. '#000000' */
  virtual WT_USTRING valueText() const;
  
  /** Sets the color of the select.  Value must be in seven-character
      hexadecimal notation, meaning the "#" character followed by two digits
      each representing red, green, and blue, like this: "#rrggbb", or else
      exception is thrown.
  */
  virtual void setValueText( const WT_USTRING &value );
  
  /** Returns Wt::DomElement_INPUT so this element will be something like:
      <input type="color" value="#000000"></input> (for native HTML5 picker)
   */
  virtual Wt::DomElementType domElementType() const;

  /** Returns whether a native HTML5 color picker will be used, or `spectrum`
    (a js library) color picker will be used.
    This is useful to know when you are inserting the color picker inside of a
    WLayout, because when it is *not* a native color picker, you need to wrap
    this widget in a WContainerWidget so it will be positioned correctly.
   
   ToDo: possibly take care of this wrapping using domElementType() and such
   */
  static bool willUseNativeColorPicker();
  
protected:
  /** Function that gets called by the JSignal to update the 'value' attribute
      in c++ land.
   */
  void colorSetCallback( const std::string &val );
  
protected:
  bool m_hasBeenSet;
  const bool m_usingNative;
  Wt::JSignal<std::string> m_userSelectedColor;
  Wt::Signal<Wt::WColor> m_cssColorChanged; //JSignal<string> could probably be made a WColor, but whatever for now.
};//class ColorSelect

#endif //ColorSelect_h
