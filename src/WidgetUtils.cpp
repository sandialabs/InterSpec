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
#include <cassert>

#include <Wt/WServer.h>
#include <Wt/WWidget.h>
#include <Wt/WApplication.h>
#include <Wt/WContainerWidget.h>

#include "InterSpec/InterSpec.h"
#include "InterSpec/WidgetUtils.h"


namespace WidgetUtils
{

namespace
{
/** Destroy `doomed` on the session thread, on the next event-loop iteration.

 The widget must be destroyed *inside* the posted lambda's body, not by letting the lambda's own
 destructor drop the last reference: `WServer::post()` -> `WServer::schedule(0ms)` ->
 `asio::post(strand_, fn)`, and asio destroys the handler on the **io-service thread** after
 invoking it.  Letting `~WWidget` run there would be a disaster - it happens with the session lock
 released and `WApplication::instance()` null, so `WWebWidget::repaint()` (reached from
 `~WContainerWidget::clear()`) dereferences a null application, and `renderOk()` /
 `EventSignalBase::prepareDestruct()` silently skip un-registering the widget from the renderers
 update map and the applications exposed-signal map, leaving dangling entries.

 Running `reset()` in the body destroys the widget on the session thread under the update lock;
 whatever copies of the (now empty) holder asio later destroys are inert.
 */
void destroy_on_session_thread( std::unique_ptr<Wt::WWidget> doomed )
{
  if( !doomed )
    return;

  Wt::WApplication * const app = Wt::WApplication::instance();
  Wt::WServer * const server = Wt::WServer::instance();

  // No session/server to post to - destroy now, on this (session) thread, rather than leak.
  if( !app || !server )
    return;

  // Park ownership on wApp rather than inside the posted task.  If the task is never delivered -
  //  the session dies between the post and its pickup, so WebController::handleApplicationEvent
  //  returns without running it - asio destroys the undelivered task on an io-service thread, and
  //  a widget owned by that task would be destroyed there.  Parked on wApp instead, the worst case
  //  is that it lives until ~WApplication destroys it, on the session thread.
  Wt::WWidget * const raw = doomed.get();
  app->addChild( std::move(doomed) );

  server->post( app->sessionId(), [raw](){
    Wt::WApplication * const a = Wt::WApplication::instance();
    if( a )
    {
      // Dropping the returned unique_ptr destroys `raw` here, on the session thread, under the
      //  update lock.
      const std::unique_ptr<Wt::WObject> owned = a->removeChild( raw );
    }
  } );
}//void destroy_on_session_thread( std::unique_ptr<Wt::WWidget> )
}//namespace


void destroyLater( std::unique_ptr<Wt::WWidget> widget )
{
  destroy_on_session_thread( std::move(widget) );
}//void destroyLater( std::unique_ptr<Wt::WWidget> )


void removeWidgetNow( Wt::WWidget *child )
{
  if( !child )
    return;

  // Dropping the returned unique_ptr destroys `child` exactly once.
  const std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();

  // A null return means nobody owned `child` through its widget parent, so nothing was destroyed
  //  (see the header for which parent shapes do this).  Callers of this function are meant to be
  //  replacing a `delete child;`, so that is virtually always a bug at the call site.
  assert( doomed );
}//void removeWidgetNow( Wt::WWidget *child )


void removeWidgetLater( Wt::WWidget *child )
{
  if( !child )
    return;

  std::unique_ptr<Wt::WWidget> doomed = child->removeFromParent();
  assert( doomed );   //see #removeWidgetNow
  destroy_on_session_thread( std::move(doomed) );
}//void removeWidgetLater( Wt::WWidget *child )


void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child )
{
  if( !parent || !child )
    return;

  std::unique_ptr<Wt::WWidget> doomed = parent->removeWidget( child );
  assert( doomed );   //see #removeWidgetNow
  destroy_on_session_thread( std::move(doomed) );
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
  if( !app )
    return nullptr;

  // WApplication::findById searches domRoot_ and domRoot2_.  It reaches widgets inside dialogs
  //  (WCompositeWidget::findById forwards to its implementation, and WTemplate::iterateChildren
  //  walks the bound widgets) and inside layouts/WStackedWidgets (WContainerWidget::iterateChildren
  //  walks children_ *and* layout_->iterateWidgets).
  return app->findById( m_id );
}


void trackSessionDialog( Wt::WDialog *dialog )
{
  InterSpec * const viewer = InterSpec::instance();
  if( viewer && dialog )
    viewer->trackToolDialog( dialog );
}//void trackSessionDialog( Wt::WDialog *dialog )

}//namespace WidgetUtils
