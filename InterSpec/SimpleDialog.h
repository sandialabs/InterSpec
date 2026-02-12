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

#include <Wt/WDialog>
#include <Wt/WString>

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


  SimpleDialog();
  SimpleDialog( const Wt::WString &title );
  SimpleDialog( const Wt::WString &title, const Wt::WString &content );
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
  
protected:
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
};//class SimpleDialog


#endif //SimpleDialog_h

