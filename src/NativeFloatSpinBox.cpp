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
#include <limits>

#include <Wt/WSignal>
#include <Wt/WString>
#include <Wt/WLineEdit>
#include <Wt/WContainerWidget>

#include "InterSpec/NativeFloatSpinBox.h"

using namespace std;
using namespace Wt;


NativeFloatSpinBox::NativeFloatSpinBox( Wt::WContainerWidget *parent )
    : WLineEdit( parent ),
      m_value( 0.0f ),
      m_min( -std::numeric_limits<float>::max() ),
      m_max( std::numeric_limits<float>::max() ),
      m_format( "%.6G" )
{
  addStyleClass( "FloatInput" );
  setAutoComplete( false );
#if( BUILD_AS_OSX_APP || IOS )
  setAttributeValue( "autocorrect", "off" );
  setAttributeValue( "spellcheck", "off" );
#endif
  setAttributeValue( "type", "number" );
  setAttributeValue( "step", "any" ); //For FF
  setValue( 0.0f );
          
  // @TODO: Since #handleChanged is connected first, it will be called last after the subsequent connections.
  //        Should figure out how to make sure #handleChanged is always called
  //        first so we can go down to only parsing the text to a float when
  //        its actually changed, and so we can get rid of m_txt.
  // Note: boost::signals2 offers call slot groups
  //       (see https://www.boost.org/doc/libs/1_72_0/doc/html/signals2/tutorial.html#id-1.3.37.4.4.3)
  //       but Wt::EventSignal doesnt seem to keep this.
  // Note: if clients of this widget always connect to #valueChanged, instead
  //       of WFormWidget::changed, then this problem also disappears.
  changed().connect( this, &NativeFloatSpinBox::handleChanged );
  enterPressed().connect( this, &NativeFloatSpinBox::handleChanged );
}
 

NativeFloatSpinBox::~NativeFloatSpinBox()
{
}
  

void NativeFloatSpinBox:: setValue( const float value )
{
  m_value = value;
  char buffer[64];
  snprintf( buffer, sizeof(buffer), m_format.c_str(), value );
  m_txt = buffer;
  setValueText( WString::fromUTF8(m_txt) );
}
 

void NativeFloatSpinBox::setSingleStep( const float step )
{
  char buffer[64];
  snprintf( buffer, sizeof(buffer), m_format.c_str(), step );
  setAttributeValue( "step", buffer );
}
  
 
void NativeFloatSpinBox::setMinimum( const float minval )
{
  m_min = minval;
  
  char buffer[64];
  snprintf( buffer, sizeof(buffer), m_format.c_str(), minval );
  setAttributeValue( "min", buffer );
}
  

void NativeFloatSpinBox::setMaximum( const float maxval )
{
  m_max = maxval;
  
  char buffer[64];
  snprintf( buffer, sizeof(buffer), m_format.c_str(), maxval );
  setAttributeValue( "max", buffer );
}

void NativeFloatSpinBox::setRange( float minval, float maxval )
{
  if( minval > maxval )
    std::swap( minval, maxval );
  setMinimum( minval );
  setMaximum( maxval );
}


float NativeFloatSpinBox::value()
{
  updateValueFromText();
  return m_value;
}


void NativeFloatSpinBox::setPlaceholderText(const WString& placeholder)
{
  //The Wt setPlaceholderText doesnt appear to be working, so manually set
  //  the attribute value
  setAttributeValue( "placeholder", placeholder );
  Wt::WLineEdit::setPlaceholderText( placeholder );
}
          

void NativeFloatSpinBox::setFormatString( const std::string &format )
{
  // Do a check that the format flag is actually valid.
  char buffer[64];
  const int nchars = snprintf( buffer, sizeof(buffer), format.c_str(), 1.23f );
  
  if( nchars <= 0 )
    throw runtime_error( "NativeFloatSpinBox::setFormatString(\"" + format
                          + "\"): invalid format flag." );
  
  m_format = format;
}


void NativeFloatSpinBox::setSpinnerHidden( const bool hidden )
{
  const char * const styleClass = "HideSpinners";
  const bool hasClass = hasStyleClass( styleClass );
  
  if( hidden && !hasClass )
    addStyleClass( styleClass );
  
  if( !hidden && hasClass )
    removeStyleClass( styleClass );
}


Signal<float> &NativeFloatSpinBox::valueChanged()
{
  return m_valueChanged;
}
          

void NativeFloatSpinBox::updateValueFromText()
{
  // Updated m_value only if displayed text has changed from what was used
  //  to set m_value.
  //  This is to avoid roundoff errors in float<-->txt.  This whole function
  //  and m_txt can be removed if we can ensure #handleChanged is always
  //  called first when the user changes value. (or client code just never
  //  connects to changed())
          
  const string valstr = valueText().toUTF8();
  if( m_txt == valstr )
    return;
  //const float newval = WLocale::currentLocale().toDouble( valueText() );
          
  float newval = 0.0f;
  const int nread = sscanf( valstr.c_str(), "%f", &newval );
  if( nread != 1 )
  {
    if( !placeholderText().empty() && valstr.empty() )
    {
      m_txt = "";
      m_value = 0.0f;
    }else
    {
      setValue( m_value );
    }
          
    return;
  }//if( failed to read in a value )
          
  if( newval < m_min || newval > m_max )
  {
    setValue( m_value );
    return;
  }

  m_txt = valstr;
  m_value = newval;
}//void updateValueFromText()

          
void NativeFloatSpinBox::handleChanged()
{
  updateValueFromText();
  m_valueChanged.emit( m_value );
}//void handleChanged()
