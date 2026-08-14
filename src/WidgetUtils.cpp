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

#include <memory>

#include <Wt/WServer.h>
#include <Wt/WWidget.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/WidgetUtils.h"


namespace WidgetUtils
{

void removeWidgetNow( Wt::WWidget *child )
{
  if( !child )
    return;

  // Dropping the returned unique_ptr destroys `child` exactly once.
  const std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();
}//void removeWidgetNow( Wt::WWidget *child )


void removeWidgetLater( Wt::WWidget *child )
{
  if( !child )
    return;

  Wt::WApplication * const app = Wt::WApplication::instance();
  std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();

  // Without a session we cant post a task; destroy now rather than leak.
  if( !app || !Wt::WServer::instance() )
    return;

  std::shared_ptr<Wt::WWidget> keep_alive( doomed.release() );
  Wt::WServer::instance()->post( app->sessionId(), [keep_alive](){} );
}//void removeWidgetLater( Wt::WWidget *child )


void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child )
{
  if( !parent || !child )
    return;

  Wt::WApplication * const app = Wt::WApplication::instance();
  std::unique_ptr<Wt::WWidget> doomed = parent->removeWidget( child );

  if( !app || !Wt::WServer::instance() )
    return;

  std::shared_ptr<Wt::WWidget> keep_alive( doomed.release() );
  Wt::WServer::instance()->post( app->sessionId(), [keep_alive](){} );
}//void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child )


WidgetHandle::WidgetHandle( Wt::WWidget *widget )
 : m_id( widget ? widget->id() : std::string() )
{
}


Wt::WWidget *WidgetHandle::resolve() const
{
  if( m_id.empty() )
    return nullptr;

  Wt::WApplication * const app = Wt::WApplication::instance();
  if( !app || !app->domRoot() )
    return nullptr;

  // domRoot() reaches dialogs too: WDialog derives from WPopupWidget, whose ctor adds it to
  //  domRoot_, and WContainerWidget::iterateChildren walks both children_ and the layouts widgets.
  return app->domRoot()->findById( m_id );
}

}//namespace WidgetUtils
