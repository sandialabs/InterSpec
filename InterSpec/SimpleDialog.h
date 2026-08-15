#ifndef SimpleDialog_h
#define SimpleDialog_h
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

#include <Wt/WDialog.h>
#include <Wt/WString.h>
#include <Wt/WApplication.h>

#include "InterSpec/WidgetUtils.h"

namespace Wt
{
  class WPushButton;
  class WContainerWidget;
}//namespace Wt

/** A simple minimal dialog meant to ask the user a question modal question where user should respond before continuing.
   Kinda similar to a iOS dialog asking a question.
 
 Shown centered in the middle of the screen.
 
 TODO:
 - If you create a SimpleDialog by clicking on a button that is in a AuxWindow, the SimpleDialog will be in back of the AuxWindow.  Need
    to implement raising the SimpleDialog to be on top; in Wt 4.4.0 at least there is a bringToFront() call, but not in 3.3.4.  A work around
    for this is to post creating the SimpleDialog to WServer, and then SimpleDialog will be created on top..
 - Maybe add a way of cancelling dialog if user clicks outside of the dialog
 - Add option to not show grey cover over the rest of the window
 - Test out more
 */
class SimpleDialog : public Wt::WDialog
{
public:
  enum SimpleDialogProperties
  {
    // Show modal background
    // Allow clicking on background to dismiss
    // Allow escape key to dismiss
  };//enum SimpleDialogProperties


  /** Create a SimpleDialog (or derived class) with proper Wt4 ownership.

   The created dialog is owned by the WApplication instance via addChild().
   The dialog will auto-delete when the user clicks a button or the dialog
   is otherwise finished.
  */
  template<typename T = SimpleDialog, typename... Args>
  static T *make( Args&&... args )
  {
    std::unique_ptr<T> ptr( new T( std::forward<Args>(args)... ) );
    T * const dialog = Wt::WApplication::instance()->addChild( std::move( ptr ) );

    // Ownership is wApp's, not the InterSpec instance's, so nothing would otherwise stop this
    //  dialog outliving the viewer on "Clear Session..." - register it so ~InterSpec can sweep it.
    WidgetUtils::trackSessionDialog( dialog );

    return dialog;
  }

  ~SimpleDialog();

  /** See notes for \c m_multipleBringToFront, but basically this is an over-ride to avoid
   jank when creating multiple SimpleDialogs at the same time (or close together anyway).
   
   You must call this function before initial render of the dialog for it to have any effect.
   */
  void doNotUseMultpleBringstoFront();
  
  /** Add a button to the footer.
   
   Buttons are added left-to-right, and clicking on them will cause the dialog to hide and become deleted, so you dont need to worry
   about cleaning up the dialog.
   
   Hookup to the returned button to trigger actions after clicking.
   */
  Wt::WPushButton *addButton( const Wt::WString &txt );
  
  /** Enables Escape-to-reject behavior, with a fallback for dialogs whose focused
      child widgets otherwise swallow Escape.
   */
  virtual void rejectWhenEscapePressed( bool enable = true );

  /** Sets the dialog's maximum width, overriding the default responsive ~50vw.

   Wt4's WLength supports viewport units, e.g. `WLength(95, Wt::LengthUnit::ViewportWidth)`.  Pass
   `WLength::Auto` to clear an override.  This is a width-only convenience around setMaximumSize()
   that preserves any max height.  The scrollable body follows the dialog's width on its own (the
   dialog layout gives it `max-width: 100%`), so only the height needs any help - see
   #updateBodySizeForWindow.
   */
  void setMaxWidth( const Wt::WLength &width );

  /** Like WDialog::setMaximumSize(), but also keeps the scrollable `.body` element's height in sync
   (the base call only reaches the outer dialog and its inner layout).  Accepts viewport units.
   */
  virtual void setMaximumSize( const Wt::WLength &width, const Wt::WLength &height ) override;

  /** Sets how much of the dialog is *not* the scrollable body - the title bar, the footer, and any
   margins - so #updateBodySizeForWindow knows how much of the window is left for the body.

   Defaults to 90 px, which is what the stock title bar plus footer take.  A dialog that hides the
   title bar or has no footer buttons (typically the phone layouts) should reduce it.
   */
  void setBodyChromeHeight( int pixels );

  /** Asks for the scrollable body to be this tall, rather than sizing to its content.

   The request is clamped to what the browser window allows, and re-clamped when the window is
   resized.  Pass a value <= 0 (the default) to size to content.
   */
  void setBodyPreferredHeight( double pixels );

  /** Sizes the scrollable `.body` so the dialog cannot outgrow the browser window.

   This arithmetic lives in C++ rather than in `SimpleDialog.css` because Wt 4 lays a dialog's
   title/body/footer out with a flex layout, whose JavaScript rewrites the body's `max-height` to
   `100%` on every reflow after copying the *inline* value onto the flex item that wraps it.  A
   limit coming from a stylesheet is therefore discarded, while one set on the widget is honoured.
   Without it a dialog with tall content grows until it hits the 95vh cap on `.simple-dialog` and
   then has its overflow - including the footer buttons - clipped, with no scrollbar.

   Called on construction and from InterSpec whenever the browser window changes size.  AuxWindow
   needs no equivalent; it re-fits itself in JavaScript (`AuxWindowOnDomResize` in AuxWindow.cpp).
   */
  void updateBodySizeForWindow();

  /** Force a dialog to be destroyed *now*, without emitting `finished()`.

   The SimpleDialog counterpart of `AuxWindow::deleteAuxWindow()`, and the only correct way to tear a
   dialog down from an owner's destructor: `done()`/`accept()`/`reject()` synchronously emit
   `finished()`, and those handlers typically call back into the (half-destroyed) owner.  Safe to
   call with a null `dialog`.

   Normal, user-driven dismissal should keep going through the buttons/`done()`, which self-destruct
   via `startDeleteSelf()`.
   */
  static void deleteSimpleDialog( SimpleDialog *dialog );

protected:
  /** Constructors are protected to enforce use of the SimpleDialog::make() factory,
   which ensures proper Wt4 widget ownership via WApplication::addChild().
  */
  SimpleDialog();
  SimpleDialog( const Wt::WString &title );
  SimpleDialog( const Wt::WString &title, const Wt::WString &content );

  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

  void init( const Wt::WString &title, const Wt::WString &content );

private:
  void startDeleteSelf();
  void deleteSelf();

protected:
  
  /** Holds the title text.
   Will have CSS style class "title".
   
   Will be null if no title text is passed in.
   
   This pointer is never accessed from this class, but may be updated by derived classes.
   */
  Wt::WText *m_title;
  
  /** Holds the message contents.
   Will have CSS style class "content".
   
   Will be null if no contents text is passed in.
   
   This pointer is never accessed from this class, but may be updated by derived classes.
   */
  Wt::WText *m_msgContents;
  
  /** When this dialog is initially rendered, some javascript will run on a delay, and a few times
   to make sure it is the front-most window (specifically for on mobile when you create a QR
   code by tapping on a AuxWindows title button, which then calls that windows bring to front
   on a delay).
   
   However, this can create some jank, specifically when you upload a picture that finds a QR
   code.  So this variable allows us to override the default behavior.
   */
  bool m_multipleBringToFront;
  
  Wt::Signals::connection m_escapeConnection1;

private:
  /** Height taken by the title bar, footer and margins; see #setBodyChromeHeight. */
  int m_bodyChromeHeight;

  /** Requested body height, or <= 0 to size to content; see #setBodyPreferredHeight. */
  double m_bodyPreferredHeight;
};//class SimpleDialog


#endif //SimpleDialog_h

