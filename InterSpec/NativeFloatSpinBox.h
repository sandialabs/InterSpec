#ifndef NativeFloatSpinBox_h
#define NativeFloatSpinBox_h
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

#include <Wt/WSignal>
#include <Wt/WLineEdit>

//Forward declerations
namespace Wt { class WString; }

/**
 Similar to WDoubleSpinBox if you called WDoubleSpinBox::setNativeControl(true),
 but probably not implemented nearly as well.  This is more an experiment to
 control formatting and precision, but should be a drop-in replacement for
 going back to WDoubleSpinBox.
 (As of Wt 3.3.4, calling setNativeControl causes a client-side exception)
 */
class NativeFloatSpinBox : public Wt::WLineEdit
{
public:
  NativeFloatSpinBox( Wt::WContainerWidget *parent = 0 );
          
  virtual ~NativeFloatSpinBox();
  
  void setValue( const float value );
  
  void setSingleStep( const float step );
  
  void setMinimum( const float minval );
  
  void setMaximum( const float maxval );

  void setRange( float minval, float maxval );
  
  float value();
  
  virtual void setPlaceholderText( const Wt::WString &placeholder );
          
  /** Signal emitted when the entered value changes from the client side.
   Not triggered when #setValue is called from the c++.
   
   Note: Do not connect to WFormWidget::changed().
         If you do the the value returned by #value will not be updated until after all connections to
         #WFormWidget::changed have been called.
   */
  Wt::Signal<float> &valueChanged();

protected:
  void updateValueFromText();
  void handleChanged();
          
  float m_value;
  float m_min, m_max;
  Wt::Signal<float> m_valueChanged;
  
  std::string m_txt;
};//NativeFloatSpinBox

#endif //NativeFloatSpinBox_h
