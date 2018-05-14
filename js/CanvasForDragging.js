/* Note: this is at the same time valid JavaScript and C++. */


WT_DECLARE_WT_MEMBER
(OverlayOnMouseUp, Wt::JavaScriptFunction, "OverlayOnMouseUp",
function( sender, e )
{
  var id = sender.id;
  var can = $('#c'+id).eq(0);
  var canEl = this.getElement( 'c' + id );
  
  if( !can )
  return;
  
  var co = can.offset();
  
  //FF will still fire the mouse up event for this canvas event if the mouse
  //  has gone out of the canvas (but mouse button stil held down).
  if( !can.data('hasMouse') )
  return;
  //  this.cancelEvent(e, this.CancelPropagate);
  
  if( canEl )
  {
    canEl.width = canEl.width;
    Wt.WT.DrawGammaLines('c'+id,false);
  }
 
  //make sure the user didnt hit escape
  if( can.data('startDragX')===null && can.data('startDragY')===null )
    return;
 
  var highbandwidth = can.data('HighBandwidth');
  var x0 = Math.round(can.data('startDragX'));
  var y0 = Math.round(can.data('startDragY'));
  var t0 = can.data('startDragT');
  var t1 = (new Date()).getTime();
  var x1 = Math.round(e.pageX - can.offset().left);
  var y1 = Math.round(e.pageY - can.offset().top);
  var dt = t1 - t0;
  
  var result = "" + dt + '&' + x0 + '&' + x1 + '&' + y0 + '&' + y1 +'&'+co.left+'&'+co.top;
  var button = this.button(e);
  if(!button)
  {
    if (this.buttons & 1)
    button = 1;
    else if (this.buttons & 2)
    button = 2;
    else if (this.buttons & 4)
    button = 4;
    else
    button = -1;
  }
  
  // Extra check for OSX, need to make Ctrl-left emulation as right click be
  // understood as a left click
  var mac=!!navigator.platform.match(/(Mac)/i);
  if (mac && e.ctrlKey && button===4)
  {
    button=1;
  }
  
  result += '&' + button;
  if (typeof e.keyCode !== 'undefined')
  result += '&' + e.keyCode;
  else result += '&0';
  if (typeof e.charCode !== 'undefined')
  result += '&' + e.charCode;
  else result += '&0';
  if (e.altKey)   result += '&1';
  else            result += '&0';
  if (e.ctrlKey)  result += '&1';
  else            result += '&0';
  if (e.metaKey)  result += '&1';
  else            result += '&0';
  if (e.shiftKey) result += '&1';
  else            result += '&0';
  var delta = this.wheelDelta(e);
  result += '&' + delta;
  
  
  if( x0!==null && y0!==null && can.data('mouseWasDrugged') )
  {
    var ady = Math.abs(y1-y0);
    var adx = Math.abs(x1-x0);
    var highbzoomout = (highbandwidth && button===1
                         && !e.shiftKey && !e.ctrlKey && x1<x0 && adx>ady && dt>100);
    if( !highbzoomout )
      Wt.emit( id, {name: 'userDragged', eventObject: sender}, result );
  }
  
  var modifiers = 0x0;
  if(e.shiftKey) modifiers |= 0x1;
  if(e.ctrlKey)  modifiers |= 0x2;
  if(e.altKey)   modifiers |= 0x4;
  if(e.metaKey)  modifiers |= 0x8;
  
  if( button===4 && !can.data('mouseWasDrugged') )
    Wt.emit( id, {name: 'rightClick', eventObject: sender}, x1, y1, modifiers, e.pageX, e.pageY );
  
  can.data('ContEst',null);
  can.data('cntrldwn',null);
  can.data('startDragX',null);
  can.data('startDragY',null);
  can.data('startDragT',null);
  
  if( highbandwidth )
  {
    window.clearTimeout(can.data('HBTimer'));
    can.data('RenderWaiting',false);
    
    Wt.emit( id, {name:'UserMouseUp'}, Math.round(x1), Math.round(y1), modifiers, dt );
    if( button===1 && !e.ctrlKey && !e.metaKey && !e.shiftKey )
      return;
  }
  
  Wt.WT.OverlayUpdateMouseCoords( can, canEl, x1, y1 );
}
);


WT_DECLARE_WT_MEMBER
(OverlayOnMouseOut, Wt::JavaScriptFunction, "OverlayOnMouseOut",
function( id )
{
  var can = $('#c'+id);
  var canElement = this.getElement('c'+id);
  if(!can||!canElement)
  return;
  
  canElement.width = canElement.width;
  can.data('startDragX',null);
  can.data('startDragY',null);
  can.data('startDragT',null);
  can.data('hasMouse',false);
  can.data('ContEst',null);
  //Update gamma lines we should have drawn
  Wt.WT.DrawGammaLines('c'+id,false);
}
);


WT_DECLARE_WT_MEMBER
(OverlayTouchBegin, Wt::JavaScriptFunction, "OverlayTouchBegin",
function(s,e)
{
  e.preventDefault();
  e.stopPropagation();
  
  try
  {
    var id = s.id;
    var can = $('#c'+id);
    var cEl = this.getElement('c'+id);
    if( !can || !cEl )
      throw 'OverlayTouchBegin Error';
    
    if( e.changedTouches.length < 1 )
      return;
    
    if( !can.data('Touches') )
      can.data('Touches', {});
    var touches = can.data('Touches');
    
    var numInitialTouch = Object.keys(touches).length;
    
    if( numInitialTouch < 1 )
    {
      if(can.data('HighBandwidth'))
        can.data('OrigZoom',true);
      
      var fcn = function()
      {
        can.data('rightTouchTimer',null);
        var tchs = can.data('Touches');
        if( !tchs || Object.keys(tchs).length!==1 )
        return;
        
        var t = tchs[Object.keys(tchs)[0]];
        var x1 = Math.round( t.lastX );
        var y1 = Math.round( t.lastY );
        if( Math.abs(x1-t.startX)>15 || Math.abs(y1-t.startY)>15 )
        return;
        
        var pagex = x1 + can.offset().left;
        if( pagex > 50 )
        pagex -= 50;
        var pagey = y1 + can.offset().top;
        Wt.emit( id, {name: 'rightClick', eventObject: s}, x1, y1, 0, pagex, pagey );
      };
      
      can.data('rightTouchTimer', setTimeout(fcn,600));
    }else
    {
      window.clearTimeout(can.data('rightTouchTimer'));
      can.data('rightTouchTimer',null);
    }//if( touches.length < 1 )
    
    for( var i = 0; i < e.changedTouches.length; i++ )
    {
      var touch = e.changedTouches[i];
      if( !touches.hasOwnProperty(touch.identifier) )
      {
        var info = new Object();
        info.id = touch.identifier;
        info.startX = touch.pageX - can.offset().left;
        info.startY = touch.pageY - can.offset().top;
        info.lastX = info.startX;
        info.lastY = info.startY;
        info.startTime = (new Date()).getTime();
        touches[info.id] = info;
      }//else, presumably there will be more than one touch and OverlayTouchChange will fix things up
    }//for( var i = 0; i < e.changedTouches.length; i++ )
    
    
    if( can.data('HighBandwidth') )
    {
      var keys = Object.keys(touches);
      var mods = 1000*keys.length;
      var info = touches[keys[keys.length-1]];
      Wt.emit( id, {name:'UserMouseDown'}, Math.round(info.startX), Math.round(info.startY), mods );
    }//if( do highbandwidth stuff )
    
    Wt.WT.OverlayTouchChange(s,e);
  }catch(prob)
  {
    Wt.emit( id, {name: 'jsException', eventObject: s}, 'OverlayTouchEnd: ' + prob );
  }
}
);



WT_DECLARE_WT_MEMBER
(IsDeletePeakSwipe, Wt::JavaScriptFunction, "IsDeletePeakSwipe",
function(touches)
{
  var keys=Object.keys(touches);
  if( keys.length !== 2 ) return false;
  var t1 = touches[keys[0]];
  var t2 = touches[keys[1]];
  
  var dy1 = t1.startY-t1.lastY;
  var dy2 = t2.startY-t2.lastY;
  var dydiff = Math.abs(dy2-dy1);
  if( dydiff > 15 ) return false;
  if( Math.abs(t1.lastX-t2.lastX) < 25 )   return false;
  var dy = Math.min(dy1,dy2);
  var dx1 = t1.startX-t1.lastX;
  var dx2 = t2.startX-t2.lastX;
  var dx = Math.abs(dx1-dx2);
  return (dy>dx&&dy>15);
}
);


WT_DECLARE_WT_MEMBER
(IsControlDragSwipe, Wt::JavaScriptFunction, "IsControlDragSwipe",
function(touches)
{
  var keys=Object.keys(touches);
  if( keys.length !== 2 ) return false;
  
  var t1 = touches[keys[0]];
  var t2 = touches[keys[1]];
  
  if( t1.startX > t1.lastX ) return false;
  if( t2.startX > t2.lastX ) return false;
  
  if( !isFinite(t1.startX) || !isFinite(t1.lastX)
  || !isFinite(t2.startX) || !isFinite(t2.lastX)
  || !isFinite(t1.startY) || !isFinite(t1.lastY)
  || !isFinite(t2.startY) || !isFinite(t2.lastY) )
  return false;
  
  var startdx = t1.startX - t2.startX;
  var nowdx = t1.lastX - t2.lastX;
  var yavrg = 0.5*(t1.startY+t2.startY);
  if( Math.abs(yavrg-t1.lastY) > 20 ) return false;
  if( Math.abs(yavrg-t2.lastY) > 20 )  return false;
  if( Math.abs(startdx-nowdx) > 20 ) return false;
  return (Math.abs(t1.lastX - t1.startX) > 30 );
}
);


WT_DECLARE_WT_MEMBER
(IsAltShiftSwipe, Wt::JavaScriptFunction, "IsAltShiftSwipe",
function(touches)
{
  var keys=Object.keys(touches);
  if( keys.length !== 2 ) return false;
  var t1 = touches[keys[0]];
  var t2 = touches[keys[1]];
  if( Math.abs(t1.startX-t2.startX) > 20 ) return false;
  if( Math.abs(t1.lastX-t2.lastX) > 25 )  return false;
  return ( (t1.lastX - t1.startX) > 30 );
}
);


WT_DECLARE_WT_MEMBER
(OverlayTouchEnd, Wt::JavaScriptFunction, "OverlayTouchEnd",
function(s,e)
{
  try
  {
    var id = s.id;
    var can = $('#c'+id);
    var cEl = this.getElement('c'+id);
    if( !can || !cEl )
      throw 'Couldnt get can element';
    var co = can.offset();
    // e.preventDefault();
    
    var touches = can.data('Touches');
    if( !touches )
    {
      console.log( 'Couldnt get Touches data' );
      return;
    }
    
    for( var i = 0; i < e.changedTouches.length; i++ )
    {
      var removed = e.changedTouches[i];
      var tid = removed.identifier;
      if( tid in touches )
      {
        touches[tid].lastX = removed.pageX - co.left;
        touches[tid].lastY = removed.pageY - co.top;
        touches[tid].lastTime = (new Date()).getTime();
      }
    }
    
    if( e.touches.length > 0 )
    {
      if(can.data('HighBandwidth'))
        can.data('OrigZoom',true);
      return;
    }
    
    cEl.width = cEl.width;
    
    var drawMode = can.data('drawMode');
    
    var ntouches = Object.keys(touches).length;
    var tip = $('#' + id + '_tip');
    var isdel = Wt.WT.IsDeletePeakSwipe(touches);
    var isCtrl = Wt.WT.IsControlDragSwipe(touches);
    var isAltShift = Wt.WT.IsAltShiftSwipe(touches);
    
    if( isdel || isCtrl || isAltShift )
    {
      tip.hide();
      var t1 = touches[Object.keys(touches)[0]];
      var t2 = touches[Object.keys(touches)[1]];
      var x0 = Math.min(t1.lastX,t2.lastX);
      var x1 = Math.max(t1.lastX,t2.lastX);
      
      Wt.WT.OverlayUpdateMouseCoords( can, cEl, x1, y1 );
      Wt.WT.DrawGammaLines('c'+id, false, x1);
      var shiftMod = (isdel || isAltShift) ? '1' : '0';
      var ctrlMod = (isCtrl && !isdel) ? '1' : '0';
      var altMod = (isAltShift && !isCtrl && !isdel) ? '1' : '0';
      
      if( altMod === '1' )
      {
        x0 = Math.min(t1.startX,t2.startX);
      }else if( ctrlMod === '1' )
      {
        x0 = Math.min(t1.startX,t2.startX);
        x1 = Math.max(t1.lastX,t2.lastX);
        if( x0 > x1 )
        x1 = [x0, x0 = x1][0];
      }//if( ctrlMod === '1' )
      
      var dt = Math.max(t1.lastTime,t2.lastTime) - Math.min(t1.startTime,t2.startTime);
      
      var result = "" + dt +'&'+x0+'&'+x1+'&'+t1.startY+'&'+t1.startY+'&'+co.left
      +'&'+co.top+"&1&0&0&" + altMod + "&" + ctrlMod +"&"
      + "0&" + shiftMod + "&0";
      Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
    }else if( ntouches===1 )
    {
      var tid = Object.keys(touches)[0];
      var touch = touches[tid];
      var x0 = Math.round( touch.startX );
      var y0 = Math.round( touch.startY );
      var x1 = Math.round( touch.lastX );
      var y1 = Math.round( touch.lastY );
      var dt = touch.lastTime - touch.startTime;
      var dx = Math.abs(x1-x0);
      var dy = Math.abs(y1-y0);
      var dragged = (dx>15 || dy>15);
      var button = 1;
      if( dt>500 && !dragged )
      button = 4;
      
      Wt.WT.OverlayUpdateMouseCoords( can, cEl, x1, y1 );
      Wt.WT.DrawGammaLines('c'+id, false, x1);
      
      var xeqn;
      if( can.data('sumpeak') && !can.data('sumpeakclick') )
      {
        xeqn = can.data('xeqn');
        if( xeqn )
          can.data('sumpeakclick',xeqn(touches[Object.keys(touches)[0]].lastX));
      }
      
      Wt.WT.DrawExpectedPeakConsequences(can,cEl,x1,y1);
      
      //XXX - this next line isnt doing what I would think it should
      if( xeqn )
      can.data('sumpeakclick',xeqn(x1));
      
      window.clearTimeout(can.data('rightTouchTimer'));
      can.data('rightTouchTimer',null);
      
      //Lets try to detect double tap event...
      var wasDoubleTap = false;
      if( dt < 200 )
      {
        //All the time period and locality limits are all arbitrary and just guessed
        var lastTUp = can.data( 'LastTouchUp' );
        if( lastTUp && ((touch.lastTime-lastTUp.time) < 700)
            && Math.abs(touch.lastX-lastTUp.x) < 35
            && Math.abs(touch.lastY-lastTUp.y) < 35 )
        {
          Wt.emit( id, {name: 'dblTap', eventObject: s}, x1, y1 );
          can.data( 'LastTouchUp', null );
          wasDoubleTap = true;
        }else if( Math.abs(touch.startX-touch.lastX) < 25 && Math.abs(touch.startY-touch.lastY) < 25 )
          can.data( 'LastTouchUp', {time:touch.lastTime, x: touch.lastX, y: touch.lastY} );
      }
      
      if( button === 4 )
      {
        tip.hide();
        //     var pagex = touch.lastX + can.offset().left;
        //     var pagey = touch.lastY + can.offset().top;
        //     Wt.emit( id, {name: 'rightClick', eventObject: s}, x1, y1, 0, pagex, pagey );
      }else if( !dragged && dt<200 && !wasDoubleTap )
      {
        Wt.emit( id, {name: 'userSingleClicked', eventObject: s}, x1, y1, 0, x1 + co.left, y1 + co.top );
        Wt.WT.OverlayShowPeakTip( id, can, cEl, x1, y1 );
        
        //Draw a line so it
        
        var cpad = can.data('chartPadding');
        if(!cpad) return;
        var yh = cEl.height - cpad.bottom - cpad.top;
        var cntxt = cEl.getContext("2d");
        if(cntxt.setLineDash)
          cntxt.setLineDash([1,2]);
        cntxt.beginPath();
        cntxt.strokeStyle = '#000000'; // black
        cntxt.moveTo(x1,cpad.top);
        cntxt.lineTo(x1,cpad.top+yh);
        cntxt.moveTo(x1-10,y1);
        cntxt.lineTo(x1+10,y1);
        cntxt.stroke();
        cntxt.closePath();
        
        if(cntxt.setLineDash)
          cntxt.setLineDash([]);
      }else if( dt > 200 && !wasDoubleTap )
      {
        if( can.data('HighBandwidth') && (!drawMode || !drawMode.highlight) )
        {
          Wt.emit( id, {name:'UserMouseUp'}, Math.round(x1), Math.round(y1), 0x4 );
        }else
        {
          var result = "" + dt + '&';
          if( drawMode && drawMode.highlight )
          result += x0+'&'+x1+'&'+y0+'&'+y1+'&'+co.left+'&'+co.top+'&'+button + '&0&0&0&0&0&0&0';
          else
          result += x0+'&'+x1+'&'+y0+'&'+y1+'&'+co.left+'&'+co.top+'&'+button + '&0&0&1&0&0&0&0';
          Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
        }
      }
    }else if( ntouches === 2 )
    {
      can.data( 'LastTouchUp', null );
      Wt.WT.DrawGammaLines('c'+id,false);
      var touch1 = touches[Object.keys(touches)[0]];
      var touch2 = touches[Object.keys(touches)[1]];
      var adx1 = Math.abs( touch1.startX - touch2.startX );
      var adx2 = Math.abs( touch1.lastX  - touch2.lastX );
      var ady1 = Math.abs( touch1.startY - touch2.startY );
      var ady2 = Math.abs( touch1.lastY  - touch2.lastY );
      var ddx = Math.abs( adx2 - adx1 );
      var ddy = Math.abs( ady2 - ady1 );
      
      var dt = Math.max(touch1.lastTime,touch2.lastTime) - Math.min(touch1.startTime,touch2.startTime);
      
      if( ddx >= ddy && ddx>20 )
      {
        if( adx1 < adx2 )
        { //zoom in, or if we are highlighting (e.g. time series chart), make equavelent of a shift-drag to the right
          if( touch1.startX > touch2.startX )
          touch2 = [touch1, touch1 = touch2][0];
          var result = "" + dt + '&';
          if( drawMode && drawMode.highlight )
          {
            result += Math.round(touch1.lastX) + '&' + Math.round(touch2.lastX)
            + '&' + Math.round(touch1.lastY) + '&' + Math.round(touch2.lastY)
            +'&'+co.left+'&'+co.top
            + '&1&0&0&0&0&0&1&0';
          }else
          {
            if( !can.data('HighBandwidth') )
              result += Math.round(touch1.startX) + '&' + Math.round(touch2.startX)
              + '&' + Math.round(touch1.startY) + '&' + Math.round(touch2.startY)
              +'&'+co.left+'&'+co.top
              + '&1&0&0&0&0&0&0&0';
          }
        
          if( result.length )
            Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
        }else
        {
          //zoom out
          var multx = adx1 / adx2;
          var cwidth = can.width();
          if( multx > 4 )      cwidth *= 0.049;
          else if( multx > 2 ) cwidth *= 0.029;
          else if( multx > 1 ) cwidth *= 0.019;
          
          var xlow = Math.round(0.5*(touch1.startX + touch2.startX) - 0.5*cwidth);
          var xhigh = Math.round(xlow + cwidth);
          var result = "" + dt + '&' + xhigh + '&' + xlow + '&' + touch1.startY + '&' + touch1.startY
          +'&'+co.left+'&'+co.top
          + '&1&0&0&0&0&0&0&0';
          if( !can.data('HighBandwidth') )
            Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
        }
      }else if( ddy>20 )
      {
        if( ady1 < ady2 )
        {
          if( touch1.startY > touch2.startY )
          touch2 = [touch1, touch1 = touch2][0];
          
          var result = "" + dt + '&' + Math.round(touch1.startX) + '&' + Math.round(touch2.startX)
          + '&' + Math.round(touch1.startY) + '&' + Math.round(touch2.startY)
          +'&'+co.left+'&'+co.top
          + '&1&0&0&0&0&0&0&0';
          Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
        }else //if( ady1 > ady2 )
        {
          var result = "" + dt + '&' + touch1.startX + '&' + touch2.startX
          + '&' + Math.round(cEl.height) + '&' + 0
          +'&'+co.left+'&'+co.top
          + '&1&0&0&0&0&0&0&0'; //button & keyCode & charCode & altKey & ctrlKey & metaKey & shiftKey & wheelDelta
          Wt.emit( id, {name: 'userDragged', eventObject: s}, result );
        }//if( ady1 < ady2 ) / else
      }//if( multx > multy ) / else
    }else
    {
      can.data( 'LastTouchUp', null );
      Wt.WT.DrawGammaLines('c'+id,false);
    }
    
    can.data('Touches', null);
  }catch(prob)
  {
    Wt.emit( id, {name: 'jsException', eventObject: s}, 'OverlayTouchEnd: ' + prob );
  }
  
}
);


WT_DECLARE_WT_MEMBER
(OverlayTouchChange, Wt::JavaScriptFunction, "OverlayTouchChange",
function(s,e)
{
  if( !e )
  {
    e = { changedTouches: [], touches: [], forceDraw: true };
  }else
  {
    e.preventDefault();
    e.stopPropagation();
  }
  
  try
  {
    var id = s.id;
    var cEl = this.getElement('c'+id);
    var can = $(cEl);
    
    if( !can || !cEl )
      throw 'OverlayTouchChange Error, no can or cEl';
    
    var touches = can.data('Touches');
    var cpadding = can.data('chartPadding');
    
    if( !touches || (e.changedTouches.length < 1 && !e.forceDraw) || !cpadding )
      return;
    
    var cbottom = cpadding.bottom;
    var ctop    = cpadding.top;
    var chartYHeight = cEl.height - cbottom - ctop;
    var cright  = cpadding.right;
    var cleft   = cpadding.left;
    var cWidth = can.width() - cright - cleft;
    var drawMode = can.data('drawMode');
    
    var tip = $('#' + id + '_tip');
    var highband = can.data('HighBandwidth');
    var serverside = can.data('ServerReferencePhotopeaks');
    var noSpecManip = can.data('NoSpecManip');
    
    var t, c;
    
    for( var i = 0; i < e.changedTouches.length; ++i )
    {
      c = e.changedTouches[i];
      if( c.identifier in touches )
      {
        t = touches[c.identifier];
        t.lastX = c.pageX - can.offset().left;
        t.lastY = c.pageY - can.offset().top;
        t.lastTime = (new Date()).getTime();
      }
    }
    
    var ntouches = Object.keys(touches).length;
//If we some touches that arent on, perhaps we should treat things a bit
//  differently...  Currently not working for some reason
/*    if( e.touches.length===1 && ntouches===2 && highband )
    {
      t = touches[Object.keys(touches)[0]];
      c = touches[Object.keys(touches)[1]];
      if( e.touches[0].identifier === Object.keys(touches)[0] )
        c = [t, t = c][0];
      t.lastX += c.lastDx;
      t.lastY += c.lastDY;
      t.lastDx = c.lastDx;
      t.lastDy = c.lastDy;
    }
*/
    
    var tip = $('#' + id + '_tip');
    tip.hide();
    
    var drawArrow = function(cntx,xstart,xend,y)
    {
      cntx.beginPath();
      cntx.strokeStyle = '#000000'; // black
      cntx.fillStyle = '#000000';
      cntx.moveTo(xstart, y);
      cntx.lineTo(xend, y);
      
      var mult = ((xstart-xend) > 0 ? 1 : -1);
      cntx.moveTo(xend,y);
      cntx.lineTo(xend + mult*10, y-4);
      cntx.lineTo(xend + mult*10, y+4);
      cntx.moveTo(xend,y);
      cntx.stroke();
      cntx.fill();
      cntx.closePath();
    };//drawArrow = function(y)
    
    
    var drawVerticalLine = function(cntx,x0,canv)
    {
      cntx.beginPath();
      cntx.strokeStyle = '#000000'; // black
      cntx.moveTo(x0, ctop);
      cntx.lineTo(x0, ctop+chartYHeight);
      cntx.stroke();
      cntx.closePath();
    };//drawVerticalLine(...)
    
    
    var drawHorizantalLine = function(cntx,y0,canv)
    {
      cntx.beginPath();
      cntx.strokeStyle = '#000000'; // black
      cntx.moveTo(cleft, y0);
      cntx.lineTo(cleft+cWidth, y0);
      cntx.stroke();
      cntx.closePath()
    };
    
    cEl.width = cEl.width;
    var context = cEl.getContext("2d");
    window.clearTimeout(can.data('prevCntrlMMTimer'));
    
    var rtt = can.data('rightTouchTimer');
    if( rtt && ntouches!==1 )
    {
      window.clearTimeout(rtt);
      can.data('rightTouchTimer',null);
    }
    
    if( !noSpecManip && Wt.WT.IsDeletePeakSwipe(touches) && !(drawMode && drawMode.highlight) )
    {
      var t1 = touches[Object.keys(touches)[0]];
      var t2 = touches[Object.keys(touches)[1]];
      var x0 = Math.min(t1.lastX,t2.lastX);
      var x1 = Math.max(t1.lastX,t2.lastX);
      
      context.beginPath();
      context.strokeStyle = '#000000'; // black
      context.fillStyle='rgba(61, 61, 61, 0.25)';
      context.fillRect(x0, ctop +5,x1-x0,cEl.height-ctop-cbottom-10);
      context.moveTo(x0, cEl.height- cbottom -5);
      context.lineTo(x0, ctop+5);
      context.moveTo(x1, cEl.height-cbottom-5);
      context.lineTo(x1, ctop+5);
      context.stroke();
      context.closePath();
      context.strokeText('Will Erase Peaks In Range', -45 + (x0+x1)/2.0, (ctop+5+cEl.height-cbottom-5)/2.0 );
      if( !serverside )
        Wt.WT.DrawGammaLines('c'+id,false);
      
      if(highband && !can.data('OrigZoom'))
      {
        window.clearTimeout(can.data('HBTimer'));
        Wt.emit(id, {name:'PinchZoom'}, Math.round(t1.startX), Math.round(t1.startX), Math.round(t2.startX), Math.round(t2.startX), 0 );
        can.data('OrigZoom',true);
      }
    }else if( !noSpecManip && Wt.WT.IsControlDragSwipe(touches) && !(drawMode && drawMode.highlight) )
    {
      var t1 = touches[Object.keys(touches)[0]];
      var t2 = touches[Object.keys(touches)[1]];
      var x0 = Math.min(t1.startX,t2.startX);
      var x1 = Math.max(t1.lastX,t2.lastX);
      if( x0 > x1 )
      x1 = [x0, x0 = x1][0];
      
      drawVerticalLine(context,x0,can);
      drawVerticalLine(context,x1,can);
      
      //Draw two blue dots on vertical lines at y+-20 pixels to show limits
      context.save();
      context.beginPath();
      context.strokeStyle = 'blue';
      context.fillStyle = 'blue';
      var yavrg = 0.5*(t1.startY+t2.startY);
      context.arc(x0, yavrg-20, 3, 0, 2 * Math.PI, false);
      context.arc(x0, yavrg+20, 3, 0, 2 * Math.PI, false);
      context.fill();
      
      var drawVArrow = function(cntx,x,ystart,yend)
      {
        cntx.beginPath();
        cntx.strokeStyle = '#000000'; // black
        cntx.fillStyle = '#000000';
        cntx.moveTo(x, ystart);
        cntx.lineTo(x, yend);
        var mult = ((ystart-yend) > 0) ? 1 : -1;
        cntx.moveTo(x,yend);
        cntx.lineTo(x-4, yend + mult*10);
        cntx.lineTo(x+4, yend + mult*10);
        cntx.moveTo(x,yend);
        cntx.stroke();
        cntx.fill();
        cntx.closePath();
      };//drawArrow = function(y)
      
      drawVArrow(context, x0-33, yavrg-40, yavrg-20-2 );
      drawVArrow(context, x0-33, yavrg+40, yavrg+20+2 );
      
      context.strokeStyle = '#ACACAC';
      context.strokeText('Keep fingers', x0-65, yavrg-4 );
      context.strokeText('level', x0-47, yavrg+8 );
      
      context.beginPath();
      context.arc(x1, yavrg-20, 3, 0, 2 * Math.PI, false);
      context.arc(x1, yavrg+20, 3, 0, 2 * Math.PI, false);
      context.fill();
      
      if( context.setLineDash )
        context.setLineDash([2,2]);
      context.beginPath();
      context.moveTo(x0, yavrg-20);
      context.lineTo(x1, yavrg-20);
      context.stroke();
      
      context.beginPath();
      context.moveTo(x0, yavrg+20);
      context.lineTo(x1, yavrg+20);
      context.stroke();
      
      context.restore();
      
      context.strokeText('Will search for peaks within', -45 + (x0+x1)/2.0, (ctop+5+cEl.height-cbottom-5)/2.0 );
      
      if( !serverside )
        Wt.WT.DrawGammaLines('c'+id,false);
      
      var ctrlLimit = can.data('ctrlKeyUpdateDt');
      if( !ctrlLimit )
      {
        Wt.WT.DrawContEst(id);
        Wt.emit( id, { name: 'cntrlMouseMove' }, Math.round(x0), Math.round(x1) );
      }else
      {
        var lcmm = can.data('prevCntrlMM');
        var now = (new Date()).getTime();
        var dt = now-lcmm;  //dt will be NaN if lcmm==null
        
        if( !lcmm || dt>=ctrlLimit )
        {
          Wt.emit( id, { name: 'cntrlMouseMove' }, Math.round(x0), Math.round(x1) );
          can.data('prevCntrlMM',now);
        }else
        {
          var pt = setTimeout( function(){Wt.WT.OverlayTouchChange(s,e);}, ctrlLimit-dt);
          can.data('prevCntrlMMTimer',pt);
        }
      }//if( !ctrlLimit ) / else
      
      if(highband && !can.data('OrigZoom'))
      {
        window.clearTimeout(can.data('HBTimer'));
        Wt.emit(id, {name:'PinchZoom'}, Math.round(t1.startX), Math.round(t1.startX), Math.round(t2.startX), Math.round(t2.startX), 0 );
        can.data('OrigZoom',true);
      }
    }else if( !noSpecManip && Wt.WT.IsAltShiftSwipe(touches) && !(drawMode && drawMode.highlight) )
    {
      var t1 = touches[Object.keys(touches)[0]];
      var t2 = touches[Object.keys(touches)[1]];
      var x0 = Math.min(t1.startX,t2.startX);
      var x1 = Math.max(t1.lastX,t2.lastX);
      
      context.beginPath();
      context.strokeStyle = '#000000'; // black
      context.fillStyle='rgba(255, 255, 0, 0.5)'; //0.6
      context.fillRect(x0, ctop +5,x1-x0,cEl.height-ctop-cbottom-10);
      context.moveTo(x0, cEl.height- cbottom -5);
      context.lineTo(x0, ctop+5);
      context.moveTo(x1, cEl.height-cbottom-5);
      context.lineTo(x1, ctop+5);
      context.stroke();
      context.closePath();
      context.strokeText('Will Count Gammas In Range', -45 + (x0+x1)/2.0, (ctop+5+cEl.height-cbottom-5)/2.0 );
      if( !serverside )
        Wt.WT.DrawGammaLines('c'+id,false);
      
      if(highband && !can.data('OrigZoom'))
      {
        window.clearTimeout(can.data('HBTimer'));
        Wt.emit(id, {name:'PinchZoom'}, Math.round(t1.startX), Math.round(t1.startX), Math.round(t2.startX), Math.round(t2.startX), 0 );
        can.data('OrigZoom',true);
      }
    }else if( ntouches===1 )
    {
      var xeqn;
      var touch = touches[Object.keys(touches)[0]];
      var x0 = Math.round( touch.startX );
      var y0 = Math.round( touch.startY );
      var x1 = Math.round( touch.lastX );
      var y1 = Math.round( touch.lastY );
      
      Wt.WT.OverlayUpdateMouseCoords( can, cEl, x1, y1 );
      if( !serverside )
        Wt.WT.DrawGammaLines('c'+id,false,x1);
        
      if( !noSpecManip )
        Wt.WT.DrawExpectedPeakConsequences( can, cEl, x1, y1 );
      
      if(highband)
        can.data('OrigZoom',false);
      
      if( (!highband || (Math.abs(x0-x1)>=15)) && rtt )
      {
        window.clearTimeout(rtt);
        can.data('rightTouchTimer',null);
      }
        
        
      if( drawMode && drawMode.highlight )
      {
        context.save();
        context.fillStyle='rgba(255, 255, 0, 0.605)';
        context.fillRect(x0,ctop,x1-x0,chartYHeight);
        drawVerticalLine(context,x0,can);
        drawVerticalLine(context,x1,can);
        context.restore();
        context.strokeText( 'Will use highlighted samples', x0+5, Math.round(ctop+0.5*chartYHeight+5) );
      }else
      {
        if( highband )
        {
          var callcount = 0;
          var pt = can.data('HBTimer');
          window.clearTimeout(pt);
          
          var doTheEmit = function()
          {
            var to;
            if( can.data('RenderWaiting') && callcount<33 )
            {
              ++callcount;
              to = setTimeout(doTheEmit, 33);
              can.data('HBTimer',to);
            }else
            {
              can.data('HBTimer',null);
              can.data('RenderWaiting',true);
              Wt.emit(id, {name:'UserMouseAltMove'}, Math.round(x1), Math.round(y1), touch.lastTime - touch.startTime );
            }
          };
            
          doTheEmit();
          
          if( !drawMode || !drawMode.highlight )
            return;
        }else
        {
          drawVerticalLine(context,x0,can);
          drawArrow( context, x0, x1, ctop + 0.05*chartYHeight );
          drawArrow( context, x0, x1, ctop + 0.35*chartYHeight );
          drawArrow( context, x0, x1, ctop + 0.65*chartYHeight );
          drawArrow( context, x0, x1, ctop + 0.95*chartYHeight );
          context.strokeText( 'Changing Displayed Energy Range', x0+5, Math.round(ctop+0.5*chartYHeight+5) );
        }
      }//if( button == 4 ) / else
    }else if( ntouches === 2 )
    {
      tip.hide();
      if( !serverside )
        Wt.WT.DrawGammaLines('c'+id,false);
        
      var touch1 = touches[Object.keys(touches)[0]];
      var touch2 = touches[Object.keys(touches)[1]];
      var adx1 = Math.abs( touch1.startX - touch2.startX );
      var adx2 = Math.abs( touch1.lastX  - touch2.lastX );
      var ady1 = Math.abs( touch1.startY - touch2.startY );
      var ady2 = Math.abs( touch1.lastY  - touch2.lastY );
      var ddx = Math.abs( adx2 - adx1 );
      var ddy = Math.abs( ady2 - ady1 );
      var areVertical = (adx2 > ady2);
      
      
      if( (highband && (ddx >= ddy || ddy < 20) && areVertical) || (!highband && ddx >= ddy && ddx>20) )
      {
        if( drawMode && drawMode.highlight )
        {
          var x0 = Math.round( touch1.lastX );
          var x1 = Math.round( touch2.lastX );
          context.save();
          context.fillStyle='rgba(255, 255, 0, 0.605)';
          context.fillRect(x0,ctop,x1-x0,chartYHeight);
          drawVerticalLine(context,x0,can);
          drawVerticalLine(context,x1,can);
          context.restore();
          context.strokeText( 'Will add/subtract highlighted samples', x0+5, Math.round(ctop+0.5*chartYHeight+5) );
        }else if( highband )
        {
          var t1x0 = Math.round(touch1.startX);
          var t1x1 = Math.round(touch1.lastX);
          var t2x0 = Math.round(touch2.startX);
          var t2x1 = Math.round(touch2.lastX);
          
          //Require at move of a few pixels, so this way the zoom wont change
          //  if the user really means to do a delete peak swipe, or define a
          //  ROI.  4 is barely enough to prevent accidental zooming, but 6
          //  is almost noticable when you zoom.
          if( ddx < 4 )
          {
//This commented out code would attempt to allow the chart to move left/right,
//  however its not working completely correct, and may not make a lot of sense
//  since this would interfere with the define of ROI.
//            var m = 0.5*(touch1.startX + touch2.startX);
//            if( t2x1 > t1x1 )  //t2 is to the right of t1
//            {
//              t1x1 = Math.round(m - 0.5*adx1);
//              t2x1 = Math.round(m + 0.5*adx1);
//            }else
//            {
//              t1x1 = Math.round(m + 0.5*adx1);
//              t2x1 = Math.round(m - 0.5*adx1);
//            }
            t1x1 = t1x0;
            t2x1 = t2x0;
          }else
          {
            var mult = (t1x1 < t2x1 && adx2 > adx1) ? 2 : -2;
            t1x1 += mult;
            t2x1 -= mult;
          }
          
          can.data('OrigZoom',false);
          var callcount = 0;
          var pt = can.data('HBTimer');
          window.clearTimeout(pt);
          
          var doTheEmit = function()
          {
            var to;
            if( can.data('RenderWaiting') && callcount<33 )
            {
              ++callcount;
              to = setTimeout(doTheEmit, 33);
              can.data('HBTimer',to);
            }else
            {
              can.data('HBTimer',null);
              can.data('RenderWaiting',true);
              Wt.emit(id, {name:'PinchZoom'}, t1x0, t1x1, t2x0, t2x1, 0 );
            }
          };
          
          doTheEmit();
        }else if( adx1 < adx2 )
        { //zoom in
          drawVerticalLine(context,touch1.startX,can);
          drawVerticalLine(context,touch2.startX,can);
          var xmin = Math.min(touch1.startX,touch2.startX);
          var xmax = Math.max(touch1.startX,touch2.startX);
          drawArrow( context, xmin, xmin-50, ctop + 0.25*chartYHeight );
          drawArrow( context, xmin, xmin-50, ctop + 0.75*chartYHeight );
          drawArrow( context, xmax, xmax+50, ctop + 0.25*chartYHeight );
          drawArrow( context, xmax, xmax+50, ctop + 0.75*chartYHeight );
          context.strokeText( 'Will Zoom In', Math.round(0.5*(touch1.startX + touch2.startX)-40), ctop+15 );
        }else
        { //zoom out
          drawVerticalLine(context,touch1.startX,can);
          drawVerticalLine(context,touch2.startX,can);
          drawArrow( context, touch1.startX, touch1.lastX, ctop + 0.25*chartYHeight );
          drawArrow( context, touch1.startX, touch1.lastX, ctop + 0.75*chartYHeight );
          drawArrow( context, touch2.startX, touch2.lastX, ctop + 0.25*chartYHeight );
          drawArrow( context, touch2.startX, touch2.lastX, ctop + 0.75*chartYHeight );
          
          var mult = adx1 / adx2;
          var zoommsg = "";
          if( mult > 4 )      zoommsg = 'all the way';
          else if( mult > 2 ) zoommsg = 'x4';
          else if( mult > 1 ) zoommsg = 'x2';
          context.strokeText( 'Will Zoom Out ' + zoommsg, Math.round(0.5*(touch1.startX + touch2.startX)-40), ctop+15 );
        }//if( adx1 < adx2 )
      }else if( ddx < ddy && ddy>20 && !(drawMode && drawMode.highlight) )
      {
        drawHorizantalLine(context,touch1.startY,can,cEl);
        drawHorizantalLine(context,touch2.startY,can,cEl);
        
        context.strokeText( (ady2 < ady1 ? 'Zoom-Out on Y-axis' : 'Zoom-In on Y-axis'),
                            -30 + cEl.width/2.0, Math.round(5 + 0.5*(touch1.startY+touch2.startY)) );
        
        if(highband && !can.data('OrigZoom'))
        {
          console.log( "Will ask to update range" );
          window.clearTimeout(can.data('HBTimer'));
          can.data('HBTimer',null);
          Wt.emit(id, {name:'PinchZoom'}, Math.round(touch1.startX), Math.round(touch1.startX), Math.round(touch2.startX), Math.round(touch2.startX), 0 );
          can.data('OrigZoom',true);
        }
      }else if( areVertical )
      {
        drawVerticalLine(context,touch1.startX,can);
        drawVerticalLine(context,touch2.startX,can);
      }
    }else
    {
//      if( !serverside )
      Wt.WT.DrawGammaLines('c'+id,false);
    }
  }catch(e)
  {
    Wt.emit( id, {name: 'jsException', eventObject: s}, 'OverlayTouchChange: ' + e );
  }
}
);


WT_DECLARE_WT_MEMBER
(OverlayOnClick, Wt::JavaScriptFunction, "OverlayOnClick",
function(s,e)
{
  var id = s.id;
  var can = $('#c'+id);
  var canElement = this.getElement('c'+id);
  if( !can || !canElement )
  {
    console.log('OverlayOnClick Error');
    return;
  }
  
  //this.cancelEvent(event, this.CancelPropagate);
  if( !can.data('mouseWasDrugged') )
  {
    var x = Math.round(e.pageX - can.offset().left);
    var y = Math.round(e.pageY - can.offset().top);
    var modifiers = 0x0;
    if(e.altKey)   modifiers |= 0x4;
    if(e.ctrlKey)  modifiers |= 0x2;
    if(e.metaKey)  modifiers |= 0x8;
    if(e.shiftKey) modifiers |= 0x1;
    
    var b = this.button(e);
    if(!b)
    {
      if (this.buttons & 1)      b = 1;
      else if (this.buttons & 2) b = 2;
      else if (this.buttons & 4) b = 4;
      else b = -1;
    }
    
    // Extra check for OSX, need to make Ctrl-left emulation as right click be
    // understood as a left click
    var mac=!!navigator.platform.match(/(Mac)/i);
    if (mac && e.ctrlKey && b===4)
    {
      b=1;
    }
    
    //Starting in Wt 3.3 we now get -1 from the above, so, well just assume its
    //  a left click (right clicks dont seem to get here).
    if( b===1 || b===-1 )
    {
      if( can.data('sumpeak') )
      {
        var xeqn = can.data('xeqn');
        if( xeqn )
        can.data('sumpeakclick',xeqn(x));
      }
      
      Wt.emit( id, {name: 'userSingleClicked', eventObject: s}, x, y, modifiers, e.pageX, e.pageY );
      canElement.width = canElement.width;
      Wt.WT.DrawGammaLines('c'+id,false,x);
      Wt.WT.OverlayUpdateMouseCoords(can,canElement,x,y);
      Wt.WT.DrawExpectedPeakConsequences(can,canElement,x,y);
    }
  }
}
);


WT_DECLARE_WT_MEMBER
(OverlayOnMouseDown, Wt::JavaScriptFunction, "OverlayOnMouseDown",
function(s,e)
{
  var id = s.id;
  var can = $('#c'+id);
  var canElement = this.getElement('c'+id);
  if( !can || !canElement )
  return;
  
  //  this.cancelEvent(e, this.CancelPropagate);
  canElement.width = canElement.width;
  var x = e.pageX - can.offset().left;
  var y = e.pageY - can.offset().top;
  can.data('startDragX', Math.round(x) );
  can.data('startDragY', Math.round(y) );
  can.data('startDragT', (new Date()).getTime() );
  can.data('mouseWasDrugged', false);
  //Update gamma lines we should have drawn
  
  if( can.data('HighBandwidth') )
  {
    var modifiers = 0x0;
    if(e.shiftKey) modifiers |= 0x1;
    if(e.ctrlKey)  modifiers |= 0x2;
    if(e.altKey)   modifiers |= 0x4;
    if(e.metaKey)  modifiers |= 0x8;
    Wt.emit( id, {name:'UserMouseDown'}, Math.round(x), Math.round(y), modifiers );
  }
  
  Wt.WT.DrawGammaLines('c'+id,false);
  
  try
  {
    var currentX = e.pageX - can.offset().left;
    var currentY = e.pageY - can.offset().top;
    Wt.WT.OverlayUpdateMouseCoords( can, canElement, currentX, currentY );
  }catch(error)
  {
  }
  
  if( e.ctrlKey )
  {
    
    var b = this.button(e);
    if(!b)
    {
      if (this.buttons & 1)      b = 1;
      else if (this.buttons & 2) b = 2;
      else if (this.buttons & 4) b = 4;
      else b = -1;
    }
    
    // Extra check for OSX, need to make Ctrl-left emulation as right click be
    // understood as a left click
    var mac=!!navigator.platform.match(/(Mac)/i);
    if (mac && b===4)
    {
      b=1;
    }
    
    if( b===1 )
    {
      can.data('cntrldwn', {x:Math.round(x), y:Math.round(y)} );
      Wt.emit( id, { name: 'cntrlMouseDown' }, Math.round(x) );
    }
  }
}
);

WT_DECLARE_WT_MEMBER
(OverlayOnMouseOver, Wt::JavaScriptFunction, "OverlayOnMouseOver",
function(s,e)
{
  var c = $('#c' + s.id );
  if( !c )
  return;
  //  this.cancelEvent(e, this.CancelPropagate);
  c.data('hasMouse', true);
}
);

WT_DECLARE_WT_MEMBER
(OverlayOnKeyPress, Wt::JavaScriptFunction, "OverlayOnKeyPress",
function(id,sender,e,signal)
{
  var can = $('#c'+id).eq(0);
  if( !can || !can.data('hasMouse') )
  return;
  var canElement = can.get(0);
  if( !canElement )
  return;
 
  //Key_Escape == 27
  if( (e.keyCode === 27) && (can.data('startDragX') !== null) )
  {
    canElement.width = canElement.width;
    can.data('startDragX',null);
    can.data('startDragY',null);
    can.data('startDragT',null);
    can.data('ContEst',null);
    
    canElement = this.getElement('c'+id);
    if( canElement )
    canElement.width = canElement.width;
    //Update gamma lines we should have drawn
    Wt.WT.DrawGammaLines('c'+id,false);
  }
  if(signal)
  {
    Wt.emit( id, {name: 'keyPressWhileMousedOver', event: e, eventObject: sender} );
    if(!e.metaKey)
    {
     e.preventDefault();
     e.stopPropagation();
    }
  }
}
);

WT_DECLARE_WT_MEMBER
(DrawExpectedPeakConsequences, Wt::JavaScriptFunction, "DrawExpectedPeakConsequences",
function( can, canElement, currentX, currentY )
{
  try
  {
    var txt;
    var xeqn = can.data('xeqn');
    var exeqn = can.data('exeqn');
    var xunitstr = can.data('xunit');
    var chartPadding = can.data('chartPadding');
    
    if( !chartPadding || !xeqn || !exeqn )
    return;
    
    var cbottom = chartPadding.bottom;
    var ctop    = chartPadding.top;
    var cright  = chartPadding.right;
    var cleft   = chartPadding.left;
    var lastx   = can.width() - cright;
    
    if( currentX < cleft || currentX > lastx )
    return;
    
    var context = canElement.getContext("2d");
    //should set the text color to blue or something
    //is there some way to figure out how many pixels long and high the text will be?
    
    var isNumber = function(o){return ! isNaN (o-0) && o != null;};
    
    var sumPeak   = can.data('sumpeak');
    var sumPeakClick = can.data('sumpeakclick');
    var compAngle = can.data('compangle');
    var compEdge  = can.data('compedge');
    var escPeaks  = can.data('escpeaks');
    
    var energy = xeqn( currentX );
    var canPairProduce = (energy > 1021.99782);
    
    if( !sumPeak && !isNumber(compAngle) && !compEdge && !(escPeaks && canPairProduce) )
    return;
    
    
    //If we're here, lets draw a line cooresponding to peaks current 'x' positons
    context.beginPath();
    var fontSize = 15;
    //    context.font = "normal " + fontSize + "px Arial";
    context.strokeStyle = '#000000'; // black
    context.moveTo(currentX, ctop);
    context.lineTo(currentX, canElement.height-cbottom);
    context.stroke();
    context.closePath();
    txt = "" + Math.round(10*energy)/10 + ' ' + xunitstr;
    context.strokeText( txt, currentX + 4, ctop+15+2*fontSize );
    
    if( sumPeak )
    {
      if( !sumPeakClick )
      {
        txt = 'Click to set sum peak first energy';
        context.save();
        context.strokeStyle = 'red';
        context.strokeText( txt, 0.5*(lastx+cleft)-50, ctop+0.25*can.height(), 100 );
        context.restore();
      }else
      {
        //should set a unique color
        var sumx = exeqn( energy + sumPeakClick );
        if( sumx > cleft && sumx < lastx )
        {
          context.beginPath();
          context.moveTo(sumx, ctop);
          context.lineTo(sumx, canElement.height-cbottom);
          context.stroke();
          context.closePath();
          txt = 'Sum Peak';
          context.strokeText( txt, sumx + 4, ctop+15+2*fontSize, 75 );
          txt = "" + Math.round(10*sumPeakClick)/10 + "+"
          + Math.round(10*energy)/10 + "="
          + Math.round(10*(energy+sumPeakClick))/10 + ' ' + xunitstr;
          context.strokeText( txt, sumx + 4, ctop+15+3*fontSize );
        }
        
        var anchorx = exeqn( sumPeakClick );
        if( anchorx > cleft && anchorx < lastx )
        {
          context.beginPath();
          context.moveTo(anchorx, ctop);
          context.lineTo(anchorx, canElement.height-cbottom);
          context.stroke();
          context.closePath();
          txt = "" + Math.round(10*sumPeakClick)/10+ ' ' + xunitstr;
          context.strokeText( txt, anchorx + 4, ctop+15+2*fontSize );
        }
        //var xval = xeqn( currentX );
      }
    }//if( sumPeak )
    
    
    
    if( escPeaks && canPairProduce )
    {
      var escapeE = energy - 510.99891;
      var escapex = exeqn( escapeE );
      if( escapex > cleft )
      {
        context.beginPath();
        context.moveTo(escapex, ctop);
        context.lineTo(escapex, canElement.height-cbottom);
        context.stroke();
        context.closePath();
        txt = 'Single Escape';
        context.strokeText( txt, escapex + 4, ctop+15 );
        txt = "" + Math.round(10*escapeE)/10 + ' ' + xunitstr;
        context.strokeText( txt, escapex + 4, ctop+15+fontSize );
      }
      
      escapeE = energy - 1021.99782;
      escapex = exeqn( escapeE );
      if( escapex > cleft )
      {
        context.beginPath();
        context.moveTo(escapex, ctop);
        context.lineTo(escapex, canElement.height-cbottom);
        context.stroke();
        context.closePath();
        txt = 'Double Escape';
        context.strokeText( txt, escapex + 4, ctop+15 );
        txt = "" + Math.round(10*escapeE)/10 + ' ' + xunitstr;
        context.strokeText( txt, escapex + 4, ctop+15+fontSize );
      }
    }//if( escPeaks && canPairProduce )
    
    
    if( isNumber(compAngle) )
    {
      var compAngleRad = compAngle*(3.14159265/180.0);
      var compEnergy = energy / (1 + ((energy/510.99891)*(1-Math.cos(compAngleRad))));
      
      compx = exeqn( compEnergy );
      if( compx > cleft )
      {
        context.beginPath();
        context.moveTo(compx, ctop);
        context.lineTo(compx, canElement.height-cbottom);
        context.stroke();
        context.closePath();
        txt = "" +  Math.round(10*compAngle)/10 + '\u00B0 Compton Peak';
        context.strokeText( txt, compx + 4, ctop+15 );
        txt = "" + Math.round(10*compEnergy)/10 + ' ' + xunitstr;
        context.strokeText( txt, compx + 4, ctop+15+fontSize );
      }
    }//if( isNumber(compAngle) )
    
    
    if( compEdge )
    {
      context.beginPath();
      var compedge = energy - (energy / (1 + (2*(energy/510.99891))));
      var edgex = exeqn( compedge );
      if( edgex > cleft )
      {
        context.moveTo(edgex, ctop);
        context.lineTo(edgex, canElement.height-cbottom);
        context.stroke();
        context.closePath();
        txt = 'Compton Edge';
        context.strokeText( txt, edgex + 4, ctop+15 );
        txt = "" + Math.round(10*compedge)/10 + ' ' + xunitstr;
        context.strokeText( txt, edgex + 4, ctop+15+fontSize );
      }
    }//if( compEdge )
  }catch(e)
  {
    console.log( 'Error drawing expected peak consequences' );
  }//try/catch
}
);

WT_DECLARE_WT_MEMBER
(DrawSearchEnergies, Wt::JavaScriptFunction, "DrawSearchEnergies",
function( id )
{
  var can = $('#'+id);
  var pad = can.data('chartPadding');
  var energies = can.data('SearchEnergies');
  var energyToPx = can.data('exeqn');
  var pxToEnergy = can.data('xeqn');

  if(!energies||!energyToPx||!pxToEnergy||!pad)
   return;

  var canEl = can.get(0);
  var cntx = canEl.getContext("2d");
  var minE = pxToEnergy(pad.left);
  var maxE = pxToEnergy(canEl.width-pad.right);

  var e, d, lpx, mpx, upx;

  cntx.save();
  cntx.strokeStyle='#4C4C4C';
  cntx.fillStyle='rgba(255, 204, 204, 0.5)';

  for(var i = 0; i < energies.length; ++i)
  {
    e = energies[i][0];
    if(e<minE||e>maxE)
      continue;

    d = energies[i][1];
    lpx = energyToPx(e-d);
    mpx = energyToPx(e);
    upx = energyToPx(e+d);

    cntx.fillRect(lpx,pad.top,upx-lpx,canEl.height-pad.bottom-pad.top);
    cntx.beginPath();
    cntx.moveTo(mpx,canEl.height-pad.bottom);
    cntx.lineTo(mpx,pad.top);
    cntx.stroke();
    cntx.closePath();
  }
  cntx.restore();
}
);

WT_DECLARE_WT_MEMBER
(DrawGammaLines, Wt::JavaScriptFunction, "DrawGammaLines",
function( id, eraseCan, mouseX )
{
  try
  {
    var can = $('#' + id);
    if( !can.length )
    return;
    
    var canEl = this.getElement(id);
    var energyToPx = can.data('exeqn');
    var pxToEnergy = can.data('xeqn');
    var noRefLineInfo = can.data('NoRefLineInfo');
    var serverside = can.data('ServerReferencePhotopeaks');
    
    if( !canEl || !energyToPx || !pxToEnergy )
      return;
    
    if( noRefLineInfo && serverside )
      return;
    
    if( eraseCan )
    canEl.width = canEl.width;
    Wt.WT.DrawSearchEnergies(id);
    
    var context = canEl.getContext("2d");
    
    var pad = can.data('chartPadding');
    if( !context || !pad )
    return;
    
    var cbottom = pad.bottom;
    var ctop    = pad.top;
    var cright  = pad.right;
    var cleft   = pad.left;
    var h = canEl.height-cbottom-ctop;
    var loglin = can.data('GammaLinesLogAndLin');
    var pxToCounts, countsToPx, minCounts, maxCounts, countrange;
    
    if( loglin )
    {
      pxToCounts = can.data('yeqn');
      countsToPx = can.data('cyeqn');
      if( !pxToCounts || !countsToPx )
      return;
      minCounts = pxToCounts(ctop+h);
      maxCounts = pxToCounts(ctop);
      countrange = maxCounts - minCounts;
    }
    
    var smallestDist = 9.0e20;
    var smallestDeltaPixel = 9.0e20;
    var nearestLine, nearestLineStyle, nearestLineParent, shielding, shieldingThickness, detector;
    
    var mousevalid = (typeof mouseX === "number");
    
    var doTheDraw = function( linedata )
    {
      try
      {
        if( !linedata )
        return;
        
        //Note that as of early 2013, only bleading edge versions of Chrome
        //  support dashed lines, so well add in a stub for other browsers
        if(!context.setLineDash)
        {
          context.setLineDash = function(){};
          //could also use mozDash or webkitLineDash...
        }
        
        var lines = linedata.lines;
        context.strokeStyle = linedata.color;
        
        var max_amp = -1.0;
        
        //Could do a slice in lines here so we only loop over the lines we need, rather
        //  than everyone
        var minE = pxToEnergy(cleft);
        var maxE = pxToEnergy(canEl.width-cright);
        
        var binarySearch = function( values, target )
        {
          //adapted from http://stackoverflow.com/questions/12369824/javascript-binary-search-insertion-preformance
          var l = 0,
          h = values.length - 1,
          m, comparison,
          comparator = function(a, b){return (a < b ? -1 : (a > b ? 1 : 0));};
          
          while( l <= h )
          {
            m = (l + h) >>> 1; /* Math.floor((l + h) / 2) */
            comparison = comparator(values[m].e, target);
            if (comparison < 0) {
              l = m + 1;
            } else if (comparison > 0) {
              h = m - 1;
            } else {
              return m;
            }
          }
          return l;
        };//var binarySearch
        
        var i,
            start = binarySearch( lines, minE ),
            end = binarySearch( lines, maxE );
        
        if( serverside )
        {
          max_amp = 1.0; 
        }else
        {
          for( i = start; i < end; ++i )
          {
            if( lines[i].e<minE || lines[i].e>maxE || lines[i].h<=0 )
              continue;
          
            if( loglin )
              max_amp = Math.max( minCounts + countrange*lines[i].h, max_amp );
            else
              max_amp = Math.max( lines[i].h, max_amp );
          }//for( var i in lines )
        }
        
        for( i = start; i < end; ++i )
        {
          var line = lines[i];
          if( line.e<minE || line.e>maxE || line.h<=0 )
            continue;
        
          var xpixel = Math.round(energyToPx(line.e));
          var d;
          if( loglin )
            d = Math.round( countsToPx(minCounts + countrange*(minCounts + countrange*line.h)/max_amp) );
          else
            d = Math.round( h*(1-line.h/max_amp) );
        
        
          if( !serverside )
          {
            context.beginPath();
        
          //if( line.decay === 'xray' )
            if( line.particle && line.particle !== 'gamma' )
              context.setLineDash([2,2]);
          
            if( loglin )
            {
              context.moveTo(xpixel,d);
            }else
            {
              if( (h-d) < 6 ) //Lets make sure the line is at least visible
              d = h-2;
              context.moveTo(xpixel,ctop+d);
            }
          
            context.lineTo(xpixel,canEl.height-cbottom);
            context.stroke();
            context.closePath();
            context.setLineDash([]);
          }//if( !serverside )
          
          if( mousevalid && (line.h>0.0) )
          {
            var dist = Math.abs(xpixel-mouseX);
            var scaleDist = (5.0+dist)/line.h;
            
            if( scaleDist<smallestDist && dist < 10)
            {
              nearestLine = line;
              nearestLineStyle = linedata.color;
              nearestLineParent = linedata.parent;
              smallestDist = scaleDist;
              smallestDeltaPixel = dist;
              shielding = linedata.shielding;
              shieldingThickness = linedata.shieldingThickness;
              detector = linedata.detector;
            }
          }//if( mouseX is valid coordinate )
        }//for( var i in lines )
      }catch(e)
      {
        console.log( 'doTheDraw: error drawing gamma lines: ' + e );
      }
    };//doTheDraw(...)
    
    var nuclines = can.data('persistedLines');
    for( var i in nuclines )
    doTheDraw( nuclines[i] );
    
    var currlines = can.data('gammalines');
    doTheDraw( currlines );
    
    if( nearestLine && smallestDeltaPixel < 20 && !noRefLineInfo )
    {
      var e = nearestLine.e;
      var x = energyToPx( e );
      var sf = nearestLine.h;
      var fontsize = 15;
      
      var txt, textdescrip;
      textdescrip = nearestLineParent + ', ' +  e + ' keV, rel. amp. ' + sf;
      
      if( nearestLine.decay )
      {
        if( nearestLine.particle === 'gamma' )
        {
          txt = nearestLine.decay;
          if( nearestLine.decay.indexOf('Capture') < 0 )
          txt += ' decay';
        }else if( nearestLine.particle === 'xray' )
        {
          textdescrip = nearestLine.el + ' x-ray, ' +  e + ' keV, I=' + sf;
        }else if( (typeof nearestLine.particle === 'string') && nearestLine.particle.length )
        {
          textdescrip = nearestLine.particle + ' from ' + nearestLineParent + ", " +  e + ' keV, I=' + sf;
        }
      }
      
      context.beginPath();
      context.strokeStyle = nearestLineStyle;
      context.fillStyle = nearestLineStyle;
      context.fillText( textdescrip, x+ 4, ctop+15 );
      if( txt )
        context.fillText( txt, x + 4, ctop+15+fontsize );
       
      if( nearestLine.particle === 'gamma' || nearestLine.particle === 'xray' )
      {
        var attTxt;
        if( shielding && shieldingThickness )
          attTxt = shieldingThickness + ' of ' + shielding;
        if( detector )
          attTxt = (attTxt ? (attTxt + ' with a ' + detector) : 'Assuming a ' + detector);
        if( attTxt )
          context.fillText( attTxt, x + 4, ctop + 15 + (txt?2:1)*fontsize );
        context.closePath();
      }
      
      //Now lets give the user an indication of which line it is
      context.beginPath();
      context.strokeStyle = "#FF0000";
      context.moveTo(x,canEl.height-cbottom+6);
      context.lineTo(x,canEl.height-cbottom);
      context.stroke();
      context.closePath();
    }//if( nearestLine )
  }catch(e)
  {
    console.log( 'Error drawing gamma lines: ' + e );
  }//try/catch
}//function( id, eraseCan )
);

WT_DECLARE_WT_MEMBER
(DrawContEst, Wt::JavaScriptFunction, "DrawContEst",
function( id )
{
  try
  {
    var canElement = Wt.WT.getElement('c'+id);
    var can = $('#c'+id);
    var coords = can.data('ContEst');
    
    if( !coords )
    return;
    
    var countsToYPix = can.data('cyeqn');
    var energyToXPix = can.data('exeqn');
    if( !countsToYPix || !energyToXPix )
    return;
    
    var x0 = energyToXPix(coords.x0);
    var y0 = countsToYPix(coords.y0);
    var x1 = energyToXPix(coords.x1);
    var y1 = countsToYPix(coords.y1);
    var dY = (coords.y1-coords.y0);
    
    var context = canElement.getContext("2d");
    context.save();
    context.beginPath();
    context.strokeStyle = "grey";
    context.lineWidth = 2;
    
    context.moveTo(x0, y0);
    for( var xpix=x0+5; xpix < x1; xpix += 5 )
    context.lineTo(xpix, countsToYPix(coords.y0+dY*(xpix-x0)/(x1-x0)));
    context.lineTo(x1, y1);
    
    context.stroke();
    context.lineWidth = 1;
    
    var lineAngle = Math.atan( (y1-y0)/(x1-x0) );
    var msg = 'approx continuum to use';
    var textw = context.measureText(msg).width*Math.cos(lineAngle);
    context.translate(0.5*(x0+x1-textw), y0-15);
    context.rotate( lineAngle );
    context.strokeText(msg, +5, -5 );
    context.restore();
  }catch(error)
  {
    console.log("Failed in DrawContEst: " + error );
  }
}
);



WT_DECLARE_WT_MEMBER
(OverlayUpdateMouseCoords, Wt::JavaScriptFunction, "OverlayUpdateMouseCoords",
function( can, canElement, currentX, currentY )
{
  try
  {
    var xeqn = can.data('xeqn');
    var yeqn = can.data('yeqn');
    var y2eqn = can.data('y2eqn');
    var xunitstr = can.data('xunit');
    
    var chartPadding = can.data('chartPadding');
    if( !chartPadding )
    return;
    
    var cbottom = chartPadding.bottom;
    var ctop    = chartPadding.top;
    var cright  = chartPadding.right;
    var cleft   = chartPadding.left;
    
    
    var context = canElement.getContext("2d");
    //should set the text color to blue or something
    //is there some way to figure out how many pixels long and high the text will be?
    
    var txt = [];
    if( yeqn )
    txt.push( 'y: ' + Math.round(10 * yeqn(currentY))/10 /*+ ' counts'*/ );
    
    if( y2eqn )
    txt.push( 'y2: ' + Math.round(10 * y2eqn(currentY))/10 /*+ ' counts'*/ );
    
    if( xeqn )
    {
      var txtstr = 'x: ' + Math.round(10 * xeqn(currentX))/10;
      if( xunitstr )
      txtstr = txtstr + ' ' + xunitstr;
      txt.push( txtstr );
    }//if( xeqn )
    
    var yoffset = 14.5*txt.length;
    
    //See http://tutorials.jenkov.com/html5-canvas/text.html for overview of positioning text
    var oldStrokeStyle = context.strokeStyle;
    var oldStrokeFont = context.font; //"10px sans-serif"
    var oldFillStyle = context.fillStyle;
    
    context.strokeStyle = '#00f'; // blue
    context.font        = '12px sans-serif'; //"16px Verdana";
    //context.textBaseline = "bottom";
                                  
    var maxWidth = 0;
    for( var i in txt )
    {
      var textMetrics = context.measureText(txt[i]);
      maxWidth = Math.max( maxWidth, textMetrics.width );
    }
                                  
    context.fillStyle='rgba(245, 245, 245, 0.8)';
    context.fillRect(canElement.width-cright-80-3, canElement.height-cbottom-yoffset-4, maxWidth+8, yoffset + 3);
    
    context.fillStyle   = '#00f';
    
    for( var i in txt )
    {
      context.fillText( txt[i], canElement.width-cright-80, canElement.height-cbottom-yoffset+8, 80 );
      yoffset -= 14.5;
    }//for( var i in txt )
                                  
    context.strokeStyle = oldStrokeStyle;
    context.fillStyle = oldFillStyle;
    context.font = oldStrokeFont;
  }catch(e)
  {
    console.log( 'Error updating mouse coordinates' );
  }//try/catch
}
);


WT_DECLARE_WT_MEMBER
(OverlayShowPeakTip, Wt::JavaScriptFunction, "OverlayShowPeakTip",
function( id, can, canElement, currentX, currentY )
{
  try
  {
    var tip = $('#'+id+'_tip');
    if(tip.length===0)
      tip = $('<div id="' + id + '_tip" class=\"peakinfo\"></div>').appendTo(can.parent());
 
 var setTipPos = function()
 {
 tip.show();
 var p=can.data('chartPadding');
 var t=can.height()-p.bottom-tip.height()-40;
 var l=can.width()-p.right-tip.width()-20;
 if(currentX>(l-20))
 {
   l=p.left+20;
   t=can.height()-p.bottom-tip.height()-5;
 }
 tip.offset({top:t,left:l });
 };
 
    var peaks = can.data('peaks');
    var dist = 99999.9;
    var nearpeak = -1;
    
    for( var i in peaks )
    {
      var peak = peaks[i];
      var d = Math.abs(currentX-peak.meanp);
      if( d < dist && (currentX > peak.xminp && currentX < peak.xmaxp) )
      {
        dist = d;
        nearpeak = i;
      }
    }
    
    if( nearpeak < 0 )
    {
      tip.hide();
      tip.data('curMean',null);
    }else if( peaks[nearpeak].mean === tip.data('curMean') )
    {
      setTipPos();
    }else
    {
      var peak = peaks[nearpeak];
      tip.data( 'curMean', peak.mean );
      var htmlTxt = ""
      +'<b>mean</b>: ' + peak.mean + ' keV<br>'
      +'<b>FWHM</b>: ' + (2.35482*peak.sigma).toFixed(2) + ' keV ('
      + (235.482*peak.sigma/peak.mean).toFixed(2) + '&#37;)<br>';
      if( peak.chi2 )
      htmlTxt += '<b>&chi;<sup>2</sup>/dof</b>: ' + peak.chi2 + '<br>';
      htmlTxt += '<b>peak area</b>: ' + peak.amp.toFixed(1) + "";
      
      if( peak.ampuncert )
      htmlTxt += '&plusmn;' + peak.ampuncert.toFixed(1) + "";
      
      htmlTxt += '<br><b>cont. area</b>: ' + peak.cont + '<br>' + "";
      
      if( peak.nuclide )
      htmlTxt += peak.nuclide + '<br>';
      
      tip.html( htmlTxt );
      //need to set the peak properties text hear
      //should transfer total area info, and continuum area info
      // to display as well
      setTipPos();
    }
  }catch(e)
  {
    //console.log( "Caught exception with peaks" )
  }
}
);

/*
WT_DECLARE_WT_MEMBER
(DrawZoomingOut, Wt::JavaScriptFunction, "DrawZoomingOut",
function(id)
{
  var can = $('#c'+id);
  var cEl = this.getElement('c'+id);
  if( !can || !cEl )
    return;
  
  var ctx = cEl.getContext("2d");
  var pad = can.data('chartPadding');
  if( !pad )
    return;
  
  ctx.save();
  
  var txt = 'Zooming Out';
  ctx.font = '16pt normal';
  var width = 60;
  try
  {
    var m = ctx.measureText(txt);
    width = m.width;
  }catch(e)
  {}
  
  
  function roundRect(x, y, w, h, rad)
  {
    ctx.beginPath();
    ctx.moveTo(x + rad, y);
    ctx.lineTo(x + w - rad, y);
    ctx.quadraticCurveTo(x + w, y, x + w, y + rad);
    ctx.lineTo(x + w, y + h - rad);
    ctx.quadraticCurveTo(x + w, y + h, x + w - rad, y + h);
    ctx.lineTo(x + rad, y + h);
    ctx.quadraticCurveTo(x, y + h, x, y + h - rad);
    ctx.lineTo(x, y + rad);
    ctx.quadraticCurveTo(x, y, x + rad, y);
    ctx.closePath();
//    ctx.stroke();
    ctx.fill();
  }
  
  var xpos = pad.left + 0.5*(cEl.width-pad.left-pad.right);
  var ypos = pad.top + 0.5*(cEl.height-pad.top-pad.bottom);
  
  ctx.fillStyle='rgba(245, 245, 245, 0.8)';
  roundRect( xpos -10 - 0.5*width, ypos-12, width + 20, 26, 4, true, true );
  
  ctx.strokeText(txt, xpos - 0.5*width, ypos+8 );
  
  ctx.restore();
});
*/
 
WT_DECLARE_WT_MEMBER
(OverlayOnMouseMove, Wt::JavaScriptFunction, "OverlayOnMouseMove",
function(sender,event)
{
  var id = sender.id;
  var can = $('#c'+id);
  var canElement = this.getElement('c'+id);
  if( !can || !canElement )
    return;
  
  var context = canElement.getContext("2d");
  
  var drawMode = can.data('drawMode');
  var chartPadding = can.data('chartPadding');
  if( !chartPadding || !drawMode )
  return;
  
  //    this.cancelEvent(event, this.CancelPropagate);
  
  var currentX = event.pageX - can.offset().left;
  var currentY = event.pageY - can.offset().top;
  var noSpecManip = can.data('NoSpecManip');
  
  
  //Update tooltip that gives peak information
  if( can.data('startDragX')===null )
  Wt.WT.OverlayShowPeakTip( id, can, canElement, currentX, currentY );
  
  
  //Clear the canvas
  canElement.width = canElement.width;
  
  
  var highlight = drawMode.highlight;
  var outline = drawMode.outline;
  var altShiftHighlight = drawMode.altShiftHighlight;
  
  var cbottom = chartPadding.bottom;
  var ctop    = chartPadding.top;
  var cright  = chartPadding.right;
  var cleft   = chartPadding.left;
  
  var startX = can.data('startDragX');
  var startY = can.data('startDragY');
  var t0 = can.data('startDragT');
  var t1 = (new Date()).getTime();
  
  var mouseButton = this.button(event);
  if(!mouseButton)
  {
    if (this.buttons & 1)
    mouseButton = 1;
    else if (this.buttons & 2)
    mouseButton = 2;
    else if (this.buttons & 4)
    mouseButton = 4;
    else
    mouseButton = -1;
  }


  if( can.data('HighBandwidth') && mouseButton===1 && !event.shiftKey && !event.ctrlKey )
  {
    var callcount = 0;  //If it takes long than a second, we're hosed anyway (I may take out callcount since its not really necassary)
    var pt = can.data('HBTimer');
    window.clearTimeout(pt);

    var evntname = 'UserMouseLeftMove';
    if( event.altKey )
      evntname = 'UserMouseAltMove';

    var doTheEmit = function()
    {
      if( can.data('RenderWaiting') && callcount<33 )
      {
        ++callcount;
        pt = setTimeout(doTheEmit, 33);
        can.data('HBTimer',pt);
      }else
      {
        can.data('RenderWaiting',true);
        Wt.emit(id, {name: evntname, eventObject: sender}, Math.round(currentX), Math.round(currentY), (t1-t0) );
        can.data('HBTimer',null);
      }
    };

    doTheEmit();
 
    if( event.altKey )
    {
      if( !highlight )
        return;
    }else
    {
      if( currentX>=startX )
        can.data('RenderWaiting',false);
//      else
//        Wt.WT.DrawZoomingOut(id);
      if( !highlight && currentX<startX )
      {
        can.data('mouseWasDrugged', true);
        return;
      }
    }
  }//if( do highbandwidth stuff )
  
  
  //Update mouse coordinate indicators
  Wt.WT.OverlayUpdateMouseCoords( can, canElement, currentX, currentY );
  
  //Update gamma lines we should have drawn
  Wt.WT.DrawGammaLines('c'+id,false,currentX);
  
  
  // Extra check for OSX, need to make Ctrl-left emulation as right click be
  // understood as a left click
  var mac=!!navigator.platform.match(/(Mac)/i);
  
  if (mac && event.ctrlKey && mouseButton===4)
  {
    mouseButton=1;
  }
  
  
  if( mouseButton <= 0 )
  Wt.WT.DrawExpectedPeakConsequences( can, canElement, currentX, currentY );
  
  if( startX===null || startY===null )
    return;
  
  can.data('mouseWasDrugged', true);
  
  var dx = currentX - startX;
  var dy = currentY - startY;
  var absDx = Math.abs(dx);
  var absDy = Math.abs(dy);
  
  var cntrldwn = can.data('cntrldwn');
  
  //if we've gotten to here, we are dragging or have the mouse down or somethign
  var tip = $('#' + id + '_tip');
  if( tip.length !== 0 )
  tip.hide();
  
  if( event.ctrlKey && event.altKey )
  return;
  
  if( event.ctrlKey && event.shiftKey )
  return;
  
  if( event.altKey && event.shiftKey && !altShiftHighlight && !highlight )
  return;
  
  if( mouseButton === 4 && (event.ctrlKey || event.altKey || event.shiftKey) )
  return;
  
  if( highlight && (mouseButton!==1 || event.ctrlKey) )
  return;
  
  if( event.ctrlKey && cntrldwn )
  {
    var ctrlLimit = can.data('ctrlKeyUpdateDt');
    if( !ctrlLimit )
    {
      Wt.WT.DrawContEst(id);
      Wt.emit( id, { name: 'cntrlMouseMove' }, cntrldwn.x, Math.round(currentX) );
    }else
    {
      var lcmm = can.data('prevCntrlMM');
      var now = (new Date()).getTime();
      var dt = now-lcmm;  //dt will be NaN if lcmm==null
      var pt = can.data('prevCntrlMMTimer');
      window.clearTimeout(pt);
      
      if( !lcmm || dt>=ctrlLimit )
      {
        Wt.emit( id, { name: 'cntrlMouseMove' }, cntrldwn.x, Math.round(currentX) );
        can.data('prevCntrlMM',now);
      }else
      {
        pt = setTimeout( function(){Wt.WT.OverlayOnMouseMove(sender,event);}, ctrlLimit-dt);
        can.data('prevCntrlMMTimer',pt);
      }
    }//if( !ctrlLimit ) / else
  }//if( event.ctrlKey )
  
  var shiftXRange = ((mouseButton===4 && !event.altKey && !event.shiftKey && !event.ctrlKey)
  || (mouseButton===1 && event.altKey && !event.shiftKey && !event.ctrlKeY && !highlight));
  
  if( shiftXRange )
  {
    var drawArrow = function(y)
    {
      context.beginPath();
      context.strokeStyle = '#000000'; // black
      context.fillStyle = '#000000';
      context.moveTo(startX, y);
      context.lineTo(currentX, y);
      
      var mult = -1;
      if( dx < 0 )
      mult = 1;
      
      context.moveTo(currentX,y);
      context.lineTo(currentX + mult*10, y-4);
      context.lineTo(currentX + mult*10, y+4);
      context.moveTo(currentX,y);
      context.stroke();
      context.fill();
      context.closePath();
    };
    
    var chartYHeight = canElement.height-cbottom-ctop;
    var ybottom = canElement.height-cbottom;
    //    context.moveTo(startX, 0.05*canElement.height);
    //    context.lineTo(startX, 0.95*canElement.height);
    context.beginPath();
    context.strokeStyle = '#000000'; // black
    context.moveTo(startX, ctop);
    context.lineTo(startX, ctop+chartYHeight);
    context.stroke();
    context.closePath();
    
    drawArrow( ctop + 0.05*chartYHeight );
    drawArrow( ctop + 0.35*chartYHeight );
    //        drawArrow( 0.5*canElement.height );
    drawArrow( ctop + 0.65*chartYHeight );
    drawArrow( ctop + 0.95*chartYHeight );
    
    if( !noSpecManip && mouseButton===4)
    {
      context.beginPath();
      context.strokeStyle = '#CCCCCC'; // grey
      context.moveTo(currentX, ctop);
      context.lineTo(currentX, ctop+chartYHeight);
      context.stroke();
      context.closePath();
      context.strokeStyle = '#000000'; // black
      var xeqn = can.data('xeqn');
      if( xeqn )
      {
        var currentE = Math.round(100*xeqn(currentX)+0.5)/100;
        var startE = Math.round(100*xeqn(startX)+0.5)/100;
        
        context.strokeText('Recalibrate data from ' + startE + ' to ' + currentE + ' keV', startX+5, 0.5*canElement.height);
      }else
      {
        context.strokeText('Recalibrate - move data in direction of arrows', startX+5, 0.5*canElement.height);
      }
    }else
    {
      context.strokeText('Changing Displayed Energy Range', startX+5, 0.5*canElement.height );
    }//if(mouseButton===4) / else
    
    return;
  }//if( shiftXRange )
  
  
  if( highlight || altShiftHighlight )
  context.fillStyle='rgba(255,255,0,0.605)';
  if( highlight && event.altKey )
    context.fillStyle='rgba(0,255,255,0.25)';
  
  // The following paints a green square on all of the canvas.
  // Used for debugging.
  //      context.fillStyle='rgba(0, 255, 0, 0.605)';
  //      context.fillRect(0, 0, 2000, 2000);
  
  var isAltShiftDragEvent = (event.shiftKey && event.altKey && altShiftHighlight && dx>0);
  
  if( absDx===absDy && !event.ctrlKey )
  {
    if( highlight )
    context.fillRect(startX,startY,dx, dy);
    
    if(outline)
    {
      context.beginPath();
      context.strokeStyle = "#544E4F";
      context.moveTo(startX, startY);
      context.lineTo(startX, currentY);
      context.lineTo(currentX, currentY);
      context.lineTo(currentX, startY);
      context.lineTo(startX, startY);
      context.stroke();
      context.closePath();
    }//if( outline )
    /*}else if( ((absDx > absDy) || event.ctrlKey) && !event.shiftKey ) 20121217 */
  }else if( (((absDx > absDy) || event.ctrlKey) && (!event.shiftKey || highlight)) || isAltShiftDragEvent )
  {
    if( highlight || isAltShiftDragEvent )
    context.fillRect(startX, ctop +5,dx,canElement.height-ctop-cbottom -5);
    
    if(outline && !event.shiftKey)
    {
      context.beginPath();
      context.strokeStyle = 'black';
      //context.strokeStyle = "#544E4F";
      context.moveTo(startX, canElement.height- cbottom );
      context.lineTo(startX, ctop );
      context.moveTo(currentX, canElement.height-cbottom );
      context.lineTo(currentX, ctop );
      context.stroke();
      context.closePath();
    }//if( outline )
    
    if( !event.ctrlKey && !highlight && !isAltShiftDragEvent && absDx>10 )
    {
      context.strokeStyle = 'black';
      
      if( dx>0 )
      context.strokeText('Zoom In', -25 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom)/2.0 );
      else
      {
        if( absDx < 0.02*canElement.width )
        context.strokeText('Zoom Out x2', -30 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom)/2.0 );
        else if( absDx < 0.04*canElement.width )
        context.strokeText('Zoom Out x4', -30 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom)/2.0 );
        else
        context.strokeText('Zoom Out Completely', -50 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom-5)/2.0 );
      }
    }
    
    if( !noSpecManip && isAltShiftDragEvent )
      context.strokeText('Count Gammas In Range', -30 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom)/2.0 );
    
    if( !noSpecManip && event.ctrlKey )
    {
      var countsToYPix = can.data('cyeqn');
      var energyToXPix = can.data('exeqn');
      
      if( cntrldwn && countsToYPix && energyToXPix )
      {
        // Now well draw some arrows to indicate something besides a zoom-in will happen
        context.strokeStyle = 'black';
        context.fillStyle = 'black';
        
        var xpixdown = cntrldwn.x;
        var yheight = 15 + canElement.height/2.0;
        
        context.beginPath();
        context.moveTo(xpixdown-35, yheight);
        context.lineTo(xpixdown, yheight);
        context.stroke();
        context.closePath();
        
        context.beginPath();
        context.moveTo(xpixdown,yheight);
        context.lineTo(xpixdown-10,yheight-4);
        context.lineTo(xpixdown-10,yheight+4);
        context.moveTo(xpixdown,yheight);
        context.fill();
        context.closePath();
        
        context.beginPath();
        context.moveTo(currentX+35, yheight);
        context.lineTo(currentX, yheight);
        context.stroke();
        context.closePath();
        
        context.beginPath();
        context.moveTo(currentX,yheight);
        context.lineTo(currentX+10,yheight-4);
        context.lineTo(currentX+10,yheight+4);
        context.moveTo(currentX,yheight);
        context.fill();
        context.closePath();
        
        yheight = canElement.height/4.0;
        context.strokeText('Will create peak inside', -50 + (xpixdown+currentX)/2.0, yheight);
      }
    }
    
  }else if( !event.shiftKey && !highlight )
  {
    if( highlight )
    context.fillRect(cleft,startY,canElement.width-cright-cleft, dy);
    
    if( outline )
    {
      context.beginPath();
      context.strokeStyle = "#544E4F";
      context.moveTo( cleft, startY);
      context.lineTo(canElement.width-cright, startY);
      context.moveTo( cleft, currentY);
      context.lineTo(canElement.width-cright, currentY);
      context.stroke();
      context.closePath();
    }//if( outline )
    
    if( absDy>10 )
    {
      if( dy > 0 )
      context.strokeText('Zoom-in on Y-axis', -30 + canElement.width/2.0, 5 + startY + 0.5*dy );
      else
      {
        if( absDy < 0.05*canElement.height )
        context.strokeText('Zoom-out on Y-axis x2', -30 + canElement.width/2.0, 5 + startY + 0.5*dy );
        else if( absDy < 0.075*canElement.height )
        context.strokeText('Zoom-out on Y-axis x4', -30 + canElement.width/2.0, 5 + startY + 0.5*dy );
        else
        context.strokeText('Zoom-out on Y-axis full', -30 + canElement.width/2.0, 5 + startY + 0.5*dy );
      }
    }
    
  }else if( !noSpecManip && event.shiftKey && !highlight )
  {
    context.beginPath();
    context.strokeStyle = '#000000'; // black
    context.fillStyle='rgba(61, 61, 61, 0.25)';
    context.fillRect(startX, ctop +5,dx,canElement.height-ctop-cbottom-10);
    context.moveTo(startX, canElement.height- cbottom -5);
    context.lineTo(startX, ctop+5);
    context.moveTo(currentX, canElement.height-cbottom-5);
    context.lineTo(currentX, ctop+5);
    context.stroke();
    context.closePath();
    
    if( !event.altKey )
    context.strokeText('Will Erase Peaks In Range', -45 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom-5)/2.0 );
    else
    context.strokeText('Will remove gamma count range', -45 + (startX+currentX)/2.0, (ctop+5+canElement.height-cbottom-5)/2.0 );
  }
}
);




