/* InterSpec: an application to analyze spectral gamma radiation data.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative emails of interspec@sandia.gov, or srb@sandia.gov.
 
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
#include <Wt/WBoostAny>
#include <Wt/WDateTime>
#include <Wt/WModelIndex>
#include <Wt/WApplication>
#include <Wt/WEnvironment>
#include <Wt/WItemDelegate>

#include "InterSpec/LocalTimeDelegate.h"


using namespace std;
using namespace Wt;

LocalTimeDelegate::LocalTimeDelegate( WObject *parent )
    : WItemDelegate( parent ),
      m_timeZoneOffset( 0 )
{
  if( wApp )
    m_timeZoneOffset = wApp->environment().timeZoneOffset();
}


LocalTimeDelegate::~LocalTimeDelegate()
{
}


WWidget *LocalTimeDelegate::update( WWidget *widget,
                            const WModelIndex &index,
                            WFlags< ViewItemRenderFlag > flags )
{
  if( flags & RenderEditing )
    throw runtime_error( "LocalTimeDelegate not for editing" );
      
  if( !(flags & RenderEditing) )
  {
    WText *text = dynamic_cast<WText *>( widget );
        
    if( !text )
    {
      widget = text = new WText();
      text->setObjectName( "t" );
      if( !index.isValid() || (index.isValid()
         && !(index.flags() & ItemIsXHTMLText)) )
        text->setTextFormat( PlainText );
    }//if( !text )
        
    try
    {
      WDateTime val = boost::any_cast<WDateTime>( index.data() );
//WLocalDateTime - not yet implemented by Wt
//        WLocalDateTime localval( val.date(), val.time(),
//                                 WLocale::currentLocale() );
//        cerr << "\n\nPutting " << val.toString( DATE_TIME_FORMAT_STR ).toUTF8()
//             << " instead of " << localval.toString( DATE_TIME_FORMAT_STR ).toUTF8()
//             << endl;
//        cerr << "timeZone()=" << WLocale::currentLocale().timeZone() << endl;
//        text->setText( localval.toString( DATE_TIME_FORMAT_STR ) );
      val = val.addSecs( 60*m_timeZoneOffset );
      text->setText( val.toString( DATE_TIME_FORMAT_STR ) );
    }catch( std::exception &e )
    {
      cerr << "LocalTimeDelegate caught: " << e.what() << endl;
      text->setText( "" );
      return widget;
    }//try / catch
  }//if( !(flags & RenderEditing) )
      
  WString tooltip = asString( index.data(ToolTipRole) );
  if( !tooltip.empty() )
    widget->setToolTip( tooltip );
      
  widget->setStyleClass( asString( index.data(StyleClassRole) ) );
  if( flags & RenderSelected )
      widget->addStyleClass(  "Wt-selected" );
      
  return widget;
}//WWidget *update(...)
