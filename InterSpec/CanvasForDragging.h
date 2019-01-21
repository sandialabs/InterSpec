#ifndef CanvasForDragging_h
#define CanvasForDragging_h
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

#include <Wt/WEvent>
#include <Wt/WPaintedWidget>


class CanvasForDragging;

namespace Wt
{
  class JSlot;
  class WWidget;
  class WPaintDevice;
  namespace Chart{ class WAbstractChart; }
}//namespace Wt

#define USE_OverlayDragEvent 1

#if( USE_OverlayDragEvent )
struct OverlayDragEvent
{
  int dt, x0, x1, y0, y1, offsetLeft, offsetTop;
  int button, keyCode, charCode;
  Wt::WFlags<Wt::KeyboardModifier> keyModifiers;
  int wheelDelta;
  void clear();
};//struct OverlayDragEvent

bool operator>>( std::istream &, OverlayDragEvent& t);
#endif

//Should consider ingeriting from a WInteractWidget instead of WPaintedWidget
//  to maybe solve issue with Wt 3.3 where canvas elements drawn in JS get
//  erased when client and server communicate.
class CanvasForDragging : public Wt::WPaintedWidget
{
public:
  //This HTML5 Canvas object will be initialized to be same size/location
  //  as 'parent', and will have an id of parent->id()+"Cover".
  //Passing in 'highlight' = true, will mean when the user drags of over the
  //  canvas, the dragged area will be highlighted in yellow; 'false'outline'=true
  //  means some grey lines will be creaded to outline an area while its being
  //  dragged.
  //
  //In part to minimize network traffic, this widget disconnects the standard
  //  Wt mouseMoved(), mouseWentUp(), and mouseWentDown() signals, and instead
  //  replaces them with the userDragged and userSingleClicked signals that
  //  are only fired when a user does a drag (or rather zoom in) operation,
  //  or a single click (in the standard sense, not HTML5 where any mouseup/down
  //  combo is a click {e.g. where the mouse stayed in one position}).
  //  All other signals are propogated through to parent widget passed in.
  //  Additionally, there is a signal for if a key is pressed while the mouse
  //  is over the widget.  When a file is dragged from the filesystem onto this
  //  canvas the file will be uploaded, and a signal will be emitted from the
  //  fileDrop() signal, notifying connected slots of its availablilty.

  CanvasForDragging( Wt::Chart::WAbstractChart *parent,
                     bool outline, bool highlight, bool altShiftHighlight );
  virtual ~CanvasForDragging();
  
  //setScrollingParent(...): sets the scrolling frame which contains the chart
  //  from which the overlay canvas should not extend beyond.  Calling
  //  this function with a NULL argument removes this parent.
  void setScrollingParent( Wt::WContainerWidget *parent );

#if( USE_OverlayDragEvent )
  Wt::JSignal<OverlayDragEvent> &userDragged();
#else
  Wt::JSignal<int,int,Wt::WMouseEvent> &userDragged();
#endif
  Wt::JSignal<int,int,int,int,int> &userSingleClicked();
  Wt::JSignal<int,int,int,int,int> &rightClickSignal();
  Wt::JSignal<Wt::WKeyEvent> &keyPressWhileMousedOver();
  Wt::JSignal<int> &controlMouseDown();
  Wt::JSignal<int,int> &controlMouseMove();  //Only emitted if highlighting is not set
  Wt::JSignal<int,int> *dblTap();
  //jsException() is mostly for debugging and will probably be removed in the
  //  future
  Wt::JSignal<std::string> *jsException();  //[only partially implemented client side] notifies you of javascript exceptions - probably is not needed for non-debug releases

  //When the user control-drags, communication is fired back to the server every
  //  time the mouse moves.  To limit how often this will actually happen, you
  //  can set a debounce period so that events will be fired back to the server
  //  no more often than the rate you set.
  //  By default, no debounce period is set, defaulting to firing signals every
  //  time mouse is moved.  Calling this function with milliseconds<=0 returns
  //  to this default behavior.
  void setControlDragDebounceTimePeriod( int milliseconds );
  int controlDragDebounceTimePeriod() const;

  //setChartPadding(): sets the chart padding in JS form m_parent
  void setChartPadding();
  
  //setInteractionMode(...): sets the chart interaction mode to javascript
  void setInteractionMode( bool outline, bool highlight, bool altShiftHighlight );
  
  //setShowRefLineInfoForMouseOver(): set wether or not the text information
  //  should be shown for the line that the mouse is currently over.  Default is
  //  to show the information.
  void setShowRefLineInfoForMouseOver( const bool show );

  //setNoSpectrumManipulationMode(): makes it so things like identiying a ROI
  //  or deleteing a peak will not be possible if called with true argument.
  //  Defaults to false (e.g. defaults to allowing manipulating interactions).
  //  Not that the signals for these events may still be emitted, but the
  //  drawing on the canvas will not be done.
  void setNoSpectrumManipulationMode( const bool nomanip );
  
  
#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/>  &mousedown();
  Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/,int/*dt in ms*/>  &mouseup();
  Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  &leftDownMouseMove();
  Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  &altLeftDownMouseMove();
  //pinchZoomChange(): {x0 of touch 1, x_curr of touch1, x0 of touch 2, x_curr of touch 2}
  Wt::JSignal<int,int/*x_t0*/,int/*x0_t1*/,int/*x_t1*/,int/*dir*/>  &pinchZoomChange();
#endif

  
protected:
  void loadJs();
  virtual void paintEvent( Wt::WPaintDevice *paintDevice );

protected:
  
  //Lets make sure this thing can only be rendered as a canvas!
  //  May cause trouble for old IE - untested as of 20140403
  virtual Wt::WPaintedWidget::Method getMethod() const;
  
#if( USE_OverlayDragEvent )
  Wt::JSignal<OverlayDragEvent> *m_userDraggedSignal;
#else
  Wt::JSignal<int/*x0*/,int/*y0*/,Wt::WMouseEvent /*mouseup event*/> *m_userDraggedSignal;
#endif
  
  Wt::JSignal<int/*x0*/,int/*y0*/,int/*modifyer*/, int/*pageX*/, int/*pageY*/> *m_userSingleClickedSignal;
  Wt::JSignal<int/*x0*/,int/*y0*/,int/*modifyer*/, int/*pageX*/, int/*pageY*/> *m_rightClickSignal;
  //m_keyPressWhileMousedOverSignal: not created until keyPressWhileMousedOver()
  //  is called
  std::unique_ptr<Wt::JSignal<Wt::WKeyEvent> > m_keyPressWhileMousedOverSignal;
  Wt::JSignal<int>             *m_controlMouseDown;
  Wt::JSignal<int/*x0*/,int/*x1*/>                  *m_controlMouseMove;
  Wt::JSignal<int,int>         *m_dblTap;
  Wt::JSignal<std::string>     *m_jsException;

#if( USE_HIGH_BANDWIDTH_INTERACTIONS )
  Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/>  *m_mousedown;
  Wt::JSignal<int/*x*/,int/*y*/,int/*mods*/,int/*dt in ms*/>  *m_mouseup;
  Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  *m_leftDownMouseMove;
  Wt::JSignal<int/*x*/,int/*y*/,int/*dt in ms*/>  *m_altLeftDownMouseMove;
  Wt::JSignal<int/*x0_t0*/,int/*x_t0*/,int/*x0_t1*/,int/*x_t1*/,int/*dir*/>  *m_pinchZoomChange;
#endif
  
  //JSlots are not controlled under the Wt memmorry managment, but we need them
  //  to stick around for as long as this class, so rather than adding a ton
  //  of member variables, well just store shared_ptrs to them (scoped_ptrs
  //  not compatible w/ stl)
  std::vector< std::shared_ptr<Wt::JSlot> > m_jslots;

  std::shared_ptr<Wt::JSlot> m_setChartPaddingSlot;
  std::shared_ptr<Wt::JSlot> m_setInteractionModeSlot;
  
  int m_controlDragDebounceTimePeriod;
  
private:
  Wt::Chart::WAbstractChart* m_parent;
};//class CanvasForDragging

#endif
