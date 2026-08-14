#ifndef WidgetUtils_h
#define WidgetUtils_h
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
#include <memory>

namespace Wt
{
  class WWidget;
  class WContainerWidget;
}//namespace Wt


/** Helpers for destroying widgets under Wt 4 ownership rules.

 In Wt 4 a widget's parent owns it through a `std::unique_ptr` (either in
 `WObject::children_`, or in a `WWidgetItem` when a layout manages it), so the Wt 3
 idiom of `delete widget;` is a double-free: `~WWebWidget` calls `removeFromParent()`,
 which hands the still-owning `unique_ptr` back, and that temporary destroys the object
 a second time.  Use these helpers (or the parent's own `removeWidget()`) instead.
 */
namespace WidgetUtils
{
  /** Detach `child` from its parent and destroy it immediately.

   The Wt 4 replacement for a plain `delete child;`.  Safe to call with a null `child`.  Do NOT use
   this while the widget being removed is emitting a signal that the current call stack is inside of
   - use #removeWidgetLater.

   IMPORTANT - these rely on `WWidget::removeFromParent()`, which returns null (i.e. destroys
   NOTHING, silently) for three parent shapes.  Asserts catch them in a developer build:
     - the widget has no *widget* parent, e.g. it is owned only through `WObject::addChild`
       (the "parked" Cat-C tool shape).  Use `owner->removeChild(w)` instead.
     - the widget is a *global* widget - anything deriving from `WPopupWidget` (`WDialog`,
       `AuxWindow`, `SimpleDialog`, `WSuggestionPopup`) or a `WPopupMenu`.  `addGlobalWidget()`
       leaves it parented to `domRoot_` but hands ownership back, so this would unparent it
       without destroying it.  Use `AuxWindow::deleteAuxWindow()` / the owning `removeChild()`.
     - the widget was added with `addWidget`/`addNew` to a container that has a *layout*; the
       layout does not manage it, so `WContainerWidget::removeWidget` finds nothing to return.
   */
  void removeWidgetNow( Wt::WWidget *child );

  /** Detach `child` from its parent now, but destroy it on the NEXT event-loop iteration.

   Required when a widget is removed in response to its OWN signal (e.g. a rows "remove"
   button): destroying it synchronously would free the emitting signal while Wt is still
   emitting it, which is a use-after-free.  Detaching immediately keeps the model/UI in
   sync; the returned unique_ptr is kept alive in a posted task and dropped (destroying the
   widget) once the emit has unwound.
   */
  void removeWidgetLater( Wt::WWidget *child );

  /** Same as #removeWidgetLater, for when the parent is conveniently in scope. */
  void removeWidgetLater( Wt::WContainerWidget *parent, Wt::WWidget *child );

  /** Take ownership of an ALREADY-DETACHED widget and destroy it on the next event-loop iteration.

   For the APIs that hand you a `unique_ptr` directly - `WMenu::removeItem()`,
   `WTabWidget::removeTab()`, `WContainerWidget::removeWidget()` - when the removal happens inside
   the removed widget's own signal emission.  Do NOT simply capture the `unique_ptr`/`shared_ptr`
   into a `WServer::post` lambda and let the lambda's destructor free it: asio destroys the task
   handler on an io-service thread, so `~WWidget` would run there.
   */
  void destroyLater( std::unique_ptr<Wt::WWidget> widget );


  /** A handle to a widget that is safe to capture into a callback which crosses a thread boundary.

   `Wt::Core::observing_ptr` (and `observable::bindSafe`) must NOT be used for this:
   `Wt::Core::observable::observers_` is a bare `std::vector` with no mutex or atomic anywhere in
   `Wt/Core/observable.cpp`, and every copy/move/destroy of an `observing_ptr` mutates it.  Since
   `Wt::WServer::post()` takes its task by const-ref and **copy-constructs** the `std::function`
   (see `WebController.h`'s `ApplicationEvent`), a completion that captures an `observing_ptr` gets
   that observer list mutated from the worker thread - racing the session thread's `~observable`.

   This handle instead stores only the widget's DOM id - an inert `std::string` - and looks the
   widget back up when it is used.

   \sa CLAUDE.md, "Tool window lifecycle"
   */
  class WidgetHandle
  {
  public:
    /** Construct from a widget (may be null).  Must be called on the session thread. */
    explicit WidgetHandle( Wt::WWidget *widget );

    /** Resolve back to the widget, or null if it no longer exists.

     **Must be called on the session thread** (i.e. from inside a `WServer::post(sessionId, …)`
     continuation, which holds the `WApplication::UpdateLock`).
     */
    Wt::WWidget *resolve() const;

    /** #resolve, dynamic_cast to `T`; null if the widget is gone or is not a `T`. */
    template<class T>
    T *resolve_as() const { return dynamic_cast<T *>( resolve() ); }

  private:
    std::string m_id;
  };//class WidgetHandle
}//namespace WidgetUtils

#endif //WidgetUtils_h
