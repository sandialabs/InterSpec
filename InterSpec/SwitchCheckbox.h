#ifndef SwitchCheckbox_h
#define SwitchCheckbox_h
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
#include <Wt/WContainerWidget>

//Forward declerations
namespace Wt
{
class WString;
class WCheckBox;
}//namespace Wt

/**  An html checkbox, but styled to display as a toggle switch.
 */
class SwitchCheckbox : public Wt::WContainerWidget
{
public:
  
  /** Constructor that will cause a label to be displayed on the right side of the widget, and for the toggle to function as a on/off, where
   the background of the switch will change color to indicate on/off.
   */
  SwitchCheckbox( const Wt::WString &rightlabel, Wt::WContainerWidget *parent = nullptr );
  
  /** Constructor that will cause there to be a label on the left side, and the right side.  The background will not change when switched,
   just the position of the round indicator.
   */
  SwitchCheckbox( const Wt::WString &left_off_label, const Wt::WString &right_on_label,
                 Wt::WContainerWidget *parent = nullptr );
  
  // Function that are straight pass-throughs to the WCheckBox functions
  bool isChecked() const;
  void setChecked( bool checked );
  virtual void setUnChecked();
  virtual void setChecked();
  
  // Using WCheckBox::changed() instead of WCheckBox::checked()/unChecked() causes some odd
  //  issues in Wt 3.3.4 (at least) where m_allowFitting->isChecked() isnt always up to date in
  //  the immediate render cycle (it is in the immediate call to to the connected signal, but
  //  seemingly not in calls during render()... havent looked into this much yet, but not enabling
  //  connecting to this signal until I figure it out
  //Wt::EventSignal<> &changed();
  Wt::EventSignal<> &checked();
  Wt::EventSignal<> &unChecked();
  
protected:
  void init( const Wt::WString &leftlabel, const Wt::WString &rightlabel );
  
  Wt::WCheckBox *m_cb;
};//SwitchCheckbox

#endif //SwitchCheckbox_h
