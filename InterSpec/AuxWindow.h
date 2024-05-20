#ifndef AuxWindow_h
#define AuxWindow_h
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

#include <Wt/WDialog>
#include <Wt/WString>
#include <Wt/WConfig.h>
#include <Wt/WJavaScript>

namespace Wt
{
  class WText;
  class WImage;
  class WPanel;
  class WLength;
  class WPushButton;
  class WGridLayout;
  class WVBoxLayout;
  class WContainerWidget;
}//namespace Wt


enum AuxWindowProperties
{
  /** Window can be modal on phones and tablets. Without this all AuxWindows
      are full screen on phones.
   */
  PhoneNotFullScreen = 0x01,
  
  /** Window can be modal on tablets.  Without this all AuxWindows
   are full screen on tablets (unless PhoneNotFullScreen is set).
   */
  TabletNotFullScreen = 0x02,
  
  /** When set, this window will be modal on PCs (e.g., no interaction with anything behind this window), if if not full screen on tablets
   and phones, will be modal then.
   */
  IsModal = 0x04,
  
  /** Styles the window to be an InterSpec help window. */
  IsHelpWindow = 0x08,
  
  /** Dont allow the window to be collapsed.  TODO: invert this to default to not allowing collapse. */
  DisableCollapse = 0x10,
  
  /** TODO: SetCloseable not actually implemented - actually the opposite should be implemented - e.g., NotCloseable. */
  SetCloseable = 0x20,
  
  /** Enables resize.
   
   Note: currently there is a bug in Wt 3.3.4 that results in the window not
   being allowed to be resized smaller than the initial size, unless minimum
   width is explicitly specified before setting window resizable.
   Thus, if this flag is used, a minimum size of 200px by 50px will be set on
   the window.
   This can be overridden by you later.
   */
  EnableResize = 0x40
};//enum AuxWindowProperties

/** AuxWindow is a resizeable, moveable, collapsable window, which is a
    specialization of a WDialog window to suit our purposes.
    AuxWindow was created due to the observation that using jQuery to make
    a moveable, resizable window suffered from performance issues on the
    client side when contents where complex.
    By default and AuxWindow is shown upon initialization, and from the Wt
    standpoint it will always be shown to aviod issues with lazy-loading and
    optimizations normally employed by Wt - javascript is used to show/hide the
    window.
 
 Note: For an AuxWindow to play nicely with the tool tabs hiding/showing
       mechanism, the following must be satisfied:
 -# place all contents in the layout returned by stretcher() OR place contents
    into the widget returned by WDialog::contents(), but not both.
 -# if you call both() then DO NOT call WDialog::contents() (or at a minimum
    add or remove widgets from it).
 -# if you add widgets to WDialog::contents(), then calling stretcher() will
    remove the widgets added through WDialog::contents()
 -# See setToolTabsVisible in InterSpec.cpp for the syntax of putting it in the
    tools tabs
     -- Collapsing/expanding memory requires hooking it up properly.
 
 At this point the AuxWindow class is growing long in the tooth, and should
 probably be re-designed, perhaps by a few new classes.  Some of the original
 motivations for this class, such as forcing Wt to always load the contents
 (even when its hidden), allowing contents to be minimized, setting window
 to be closable, or custom title bar, or a footer, are either no longer needed
 (e.g., features have been added to Wt), or there are better ways to handle
 things.
 
 Having not yet looked at the WDialog class in Wt 4, I think the next approach
 to creating this class should probably go the route of starting with the
 WDialog source code, and modifying it to give the functionality we want.  The
 subclassing of WDialog doesnt give ideal results in terms of programming
 interface, or rendering options.
 
 \todo implement the animation for showing/hiding
 \todo how the window is initially centered could be improved
 \todo Fix fact that in Internet Explorer, the close button is as far left
       as it can be - when in fact it should float to the right
 
*/
class AuxWindow : public Wt::WDialog
{
public:
  //By default AuxWindow will be shown, centered in the window, at 50% of
  //  browser size.  Wt will assume that the window is visible as well, for
  //  the purposes of lazy loading of content.
  AuxWindow( const Wt::WString &windowTitle, Wt::WFlags<AuxWindowProperties> properties = Wt::WFlags<AuxWindowProperties>(0) );
  virtual ~AuxWindow();

  
  //Override WDialog's setResizable
  void setResizable(bool resizable);
    
  //Override WDialog's setModal
  void setModal(bool modal);
  //Override footer
  Wt::WContainerWidget* footer();

  //Override footer
//  void resize(Wt::WLength width, Wt::WLength height);
    
  //deleteAuxWindow(...): just a convienience funtion to delete the window; ex.
  // AuxWindow *w = new AuxWindow( "Example Window" );
  // w->rejectWhenEscapePressed();
  // w->setClosable( true );
  // w->finished( boost::bind( &AuxWindow:::deleteAuxWindow, window ) );
  static void deleteAuxWindow( AuxWindow *window );

  /** A convenience function to call #deleteAuxWindow via binding to a signal. */
  void deleteSelf();
  
  //rejectWhenEscapePressed(): sets it up so hitting escape will cause the
  //  finished() signal to be emmitted; also hides the dialog
  virtual void rejectWhenEscapePressed( bool enable = true );
  
  //setWindowTitle
  virtual void setWindowTitle( const Wt::WString& windowTitle );

  /** If a full-screen window on a mobile device, this function does nothing;
      otherwise if just calls WDialog::setMaximumSize(width,height).
   */
  virtual void setMaximumSize( const Wt::WLength &width, const Wt::WLength &height );
  
  //The show(), hide(), and setHidden() functions call the javascript versions
  //  of these functions, and do not modify Wt's view of if this window is shown
  //  or not, which can be important since Wt uses lazy loading of contents.
  //  If you wish to enable lazy loading (and then subsequent showing), use
  //  the WDialog::setHidden() function.
  virtual void show();
  virtual void hide();
  virtual void expand(); //executes jsExpand()
  
  //TODO: currently the WAnimation is not used
  void setHidden( bool hide, const Wt::WAnimation &animation=Wt::WAnimation() );

  //isVisible() and isHidden() return whether or not the window is visible
  //  on the client side (not what Wt thinks)
  virtual bool isVisible() const;
  virtual bool isHidden() const;

  //centerWindow(): centers the window in the browser viewport.  Only works
  //  if it is called after the window has been made visible.
  void centerWindow();
  
  void centerWindowHeavyHanded();
  
  //resizeToFitOnScreen(): executes some javascript that checks if width/height
  //  are larger than the screen, and if so, resizes the window, and sets the
  //  contents() overflow to auto.
  void resizeToFitOnScreen();

  
  //repositionWindow(...): repositions window to coordinates <x,y> of the
  //  browsers viewport.
  //If x==-32768 or y==-32768, than that dimension will bot be repositioned
  void repositionWindow( int x, int y );
  

  /** Returns the JavaScript #repositionWindow will execute. */
  std::string repositionWindowJs( int x, int y );
  
  //resizeWindow(...): resizes this AuxWindow to the specified size in pixels
  void resizeWindow( int width, int height );

  
  //resizeScaledWindow(...): resizes this AuxWindow as a fraction of the browser
  //  window size.
  //  If a (zerro or) negative value is specified for one of the dimensions,
  //  that dimension wont be resized.
  void resizeScaledWindow( double xRatio, double yRatio );

  /** Returns the JavaScript #resizeScaledWindow will execute. */
  std::string resizeScaledWindowJs( double xRatio, double yRatio ) const;
  
  //contents() is actually hiding WDialog::contents(), not overiding it, so
  //  dont expect a call from a WDialog pointer to give you what you want.
  //Also, this function will generate an assert if you have already called
  //  stretcher().
  Wt::WContainerWidget *contents() const;

  //stretcher(): initializes (if necessary) and returns the WGridLayout used to
  //  stretch contents to the body of theis dialog.
  virtual Wt::WGridLayout *stretcher();

  //disableCollapse(): hides the icon, and clears the m_collapseSlot function
  void disableCollapse();

  //The JSlot's dont modify Wt's server-side view of if it's is visible or not.
  //  You can execute them by calling exec( "null", "{quietly: true}" );
  //  where the quietly indicates dont have the client javascript code notify
  //  the server of this call.
  virtual Wt::JSlot &jsHide(); //expands window before hiding
  virtual Wt::JSlot &jsShow(); //doesnt expand window if previously collapsed
  virtual Wt::JSlot &jsExpand();
  virtual Wt::JSlot &jsCollapse();
  virtual Wt::JSlot &jsToggleExpandState();

  //Signals you can connect to to know when the window is
  //  opened/collapsed/expanded
  //  To listen for the window being closed, connect to the finished() signal.
//  Wt::JSignal<> &closed();
  Wt::JSignal<> &opened();
  Wt::JSignal<> &collapsed();
  Wt::JSignal<> &expanded();

  //should we overide WWidget::enable() and WWidget::disable() ?

  virtual void setClosable( bool closable );

  //addCloseButtonToFooter(): adds a style appropriate (phone/desktop)
  //  close button to the footer.  Desktop will display the text "Close", while
  //  mobile is a "Back" button.
  //  The button contains either the styleclass "PhoneDialogClose" or
  //  "DialogClose"
  //The button is not connected to any slots (e.g. it does not actually close
  //  this dialog).
  Wt::WPushButton *addCloseButtonToFooter( Wt::WString override_txt = "Close",
                                           const bool float_right = false,
                                           Wt::WContainerWidget *footerOverride = nullptr );
  
  //Help button add to footer
  static void addHelpInFooter(Wt::WContainerWidget *footer, std::string page );
//  static void openHelpWindow(std::string page,AuxWindow *parent);
  
  void emitReject();
  
  /** Returns if dialog is customized for phones. */
  bool isPhone() const;
  
protected:
  virtual void render( Wt::WFlags<Wt::RenderFlag> flags );

protected:
  bool m_auxIsHidden;
  bool m_modalOrig;  //original modal value
  Wt::WText       *m_titleText;
  Wt::WImage      *m_collapseIcon;
  Wt::WImage      *m_expandIcon;
  Wt::WInteractWidget *m_closeIcon;

  std::unique_ptr<Wt::JSlot>      m_collapseSlot;
  std::unique_ptr<Wt::JSlot>      m_expandSlot;
  std::unique_ptr<Wt::JSlot>      m_showSlot;
  std::unique_ptr<Wt::JSlot>      m_hideSlot;
  std::unique_ptr<Wt::JSlot>      m_toggleExpandedStatusSlot;
  std::unique_ptr<Wt::JSignal<> > m_closedSignal;
  std::unique_ptr<Wt::JSignal<> > m_openedSignal;
  std::unique_ptr<Wt::JSignal<> > m_collapsedSignal;
  std::unique_ptr<Wt::JSignal<> > m_expandedSignal;

#define USE_NEW_AUXWINDOW_ISH 0
#if( USE_NEW_AUXWINDOW_ISH )
  std::unique_ptr<Wt::JSlot> m_repositionSlot;   //->exec(id(),"{top:5,left:10}");
  std::unique_ptr<Wt::JSlot> m_centerSlot;       //->exec(id(),"null");
  std::unique_ptr<Wt::JSlot> m_resizeSlot;       //->exec(id(),"{x:200,y:150}");
  std::unique_ptr<Wt::JSlot> m_resizeScaledSlot; //->exec(id(),"{x:0.6,y:0.75}");
#endif

  
  Wt::WGridLayout *m_contentStretcher;
  
  bool m_destructing;
  
  bool m_escapeIsReject;
  Wt::Signals::connection m_escapeConnection1, m_escapeConnection2;
  
  /** Set to true if on a phone, and contructor wasnt called with PhoneNotFullScreen. */
  bool m_isPhone;
  
  /** Set to true if on a tablet, and contructor wasnt called with PhoneNotFullScreen or
   TabletNotFullScreen.
   */
  bool m_isTablet;
  
  /** Set to true if user agent string contains "Android" */
  bool m_isAndroid;
  
  Wt::WContainerWidget* m_footer;
};//class AuxWindow

#endif // AuxWindow_h

