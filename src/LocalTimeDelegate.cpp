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

#include <Wt/WAny.h>
#include <Wt/WText.h>
#include <Wt/WDateTime.h>
#include <Wt/WModelIndex.h>
#include <Wt/WApplication.h>
#include <Wt/WEnvironment.h>
#include <Wt/WItemDelegate.h>

#include "InterSpec/LocalTimeDelegate.h"


using namespace std;
using namespace Wt;

LocalTimeDelegate::LocalTimeDelegate()
    : WItemDelegate(),
      m_timeZoneOffset( 0 )
{
  if( wApp )
    m_timeZoneOffset = wApp->environment().timeZoneOffset();
}


LocalTimeDelegate::~LocalTimeDelegate()
{
}


std::unique_ptr<WWidget> LocalTimeDelegate::update( WWidget *widget,
                            const WModelIndex &index,
                            WFlags< ViewItemRenderFlag > flags )
{
  if( flags.test( ViewItemRenderFlag::Editing ) )
    throw runtime_error( "LocalTimeDelegate not for editing" );

  std::unique_ptr<WWidget> newWidget;

  WText *text = dynamic_cast<WText *>( widget );

  if( !text )
  {
    auto textPtr = std::make_unique<WText>();
    text = textPtr.get();
    newWidget = std::move( textPtr );
    text->setObjectName( "t" );
    if( !index.isValid() || (index.isValid()
       && !(index.flags().test( ItemFlag::XHTMLText ))) )
      text->setTextFormat( TextFormat::Plain );
  }//if( !text )

  try
  {
    WDateTime val = Wt::cpp17::any_cast<WDateTime>( index.data() );
    const int offsetSecs = 60 * static_cast<int>( m_timeZoneOffset.count() );
    val = val.addSecs( offsetSecs );
    text->setText( val.toString( DATE_TIME_FORMAT_STR ) );
  }catch( std::exception &e )
  {
    cerr << "LocalTimeDelegate caught: " << e.what() << endl;
    text->setText( "" );
    return newWidget;
  }//try / catch

  WWidget *w = newWidget ? newWidget.get() : widget;

  WString tooltip = asString( index.data( ItemDataRole::ToolTip ) );
  if( !tooltip.empty() )
    w->setToolTip( tooltip );

  w->setStyleClass( asString( index.data( ItemDataRole::StyleClass ) ) );
  if( flags.test( ViewItemRenderFlag::Selected ) )
    w->addStyleClass( "Wt-selected" );

  return newWidget;
}//std::unique_ptr<WWidget> update(...)
