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


#include <Wt/WText>
#include <Wt/WString>
#include <Wt/WCheckBox>
#include <Wt/WApplication>
#include <Wt/WContainerWidget>

#include "InterSpec/SwitchCheckbox.h"

using namespace std;
using namespace Wt;

SwitchCheckbox::SwitchCheckbox( const Wt::WString &rightlabel, Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
   m_cb( nullptr )
{
  init( "", rightlabel );
}

SwitchCheckbox::SwitchCheckbox( const Wt::WString &leftlabel, const Wt::WString &rightlabel,
               Wt::WContainerWidget *parent )
 : WContainerWidget( parent ),
   m_cb( nullptr )
{
  init( leftlabel, rightlabel );
}



void SwitchCheckbox::init( const Wt::WString &leftlabel, const Wt::WString &rightlabel )
{
  wApp->useStyleSheet( "InterSpec_resources/SwitchCheckbox.css" );
  addStyleClass( "SwitchCheckbox" );
  
  // Wt only lets you call WLabel::setBuddy for one label per form widget, so well use some
  //  Javascript to put some labels before/after the checkbox.
  //  Also, using WLabel::setAttributeValue("for",m_cb->id()) doesnt seem to work (WLabel is a
  //  <span> element until you set it as a buddy).
  //  So instead we'll put left/right label html element into a WText - ehh, it works
  WText *leftText = nullptr, *rightText = nullptr;
  
  if( !leftlabel.empty() )
    leftText = new WText( this );
  
  // Note that if we create a WCheckBox with some text for its label, the HTML structure will be
  //  totally different, so we will always make our own labels.
  m_cb = new WCheckBox( this );
  
  if( !rightlabel.empty() )
    rightText = new WText( this );
  
  if( leftText )
    leftText->setText( "<label for=\"" + m_cb->id() + "\">" + leftlabel.toUTF8() + "</label>" );
  
  if( rightText )
    rightText->setText( "<label for=\"" + m_cb->id() + "\">" + rightlabel.toUTF8() + "</label>" );
  
  if( !leftText || !rightText )
    m_cb->addStyleClass( "onoff" );
  else
    m_cb->addStyleClass( "twoopt" );
}//void init()


bool SwitchCheckbox::isChecked() const
{
  return m_cb->isChecked();
}

void SwitchCheckbox::setChecked( bool checked )
{
  m_cb->setChecked( checked );
}

void SwitchCheckbox::setUnChecked()
{
  m_cb->setUnChecked();
}

void SwitchCheckbox::setChecked()
{
  m_cb->setChecked();
}

Wt::EventSignal<> &SwitchCheckbox::changed()
{
  return m_cb->changed();
}

Wt::EventSignal<> &SwitchCheckbox::checked()
{
  return m_cb->checked();
}

Wt::EventSignal<> &SwitchCheckbox::unChecked()
{
  return m_cb->unChecked();
}

