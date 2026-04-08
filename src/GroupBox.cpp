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

#include <Wt/WText.h>
#include <Wt/WConfig.h>

#include "InterSpec/GroupBox.h"


GroupBox::GroupBox()
{
  init();
}


GroupBox::GroupBox( const Wt::WString &title )
{
  init();
  setTitle( title );
}


void GroupBox::init()
{
  // Set the same JS layout members that WGroupBox sets, so the Wt JS layout
  // system can properly compute preferred sizes for this widget.
  // This is important when the widget is inside a flex or layout-managed container.
  setJavaScriptMember( "wtGetPS", WT_CLASS ".LastGetPS" );
  setJavaScriptMember( "wtGetExtraMS", WT_CLASS ".LastGetPS" );
}


Wt::DomElementType GroupBox::domElementType() const
{
  return Wt::DomElementType::FIELDSET;
}


void GroupBox::setTitle( const Wt::WString &title )
{
  // The legend is always the first child; create or update it.
  if( count() == 0 || !dynamic_cast<Wt::WText *>( widget( 0 ) ) )
  {
    Wt::WText *legend = insertNew<Wt::WText>( 0, title );
    legend->setHtmlTagName( "legend" );
  }
  else
  {
    static_cast<Wt::WText *>( widget( 0 ) )->setText( title );
  }
}


const Wt::WString &GroupBox::title() const
{
  static const Wt::WString empty;

  if( count() > 0 )
  {
    const Wt::WText *legend = dynamic_cast<const Wt::WText *>( widget( 0 ) );
    if( legend )
      return legend->text();
  }

  return empty;
}
