/* Note: this is at the same time valid JavaScript and C++. */


WT_DECLARE_WT_MEMBER
(LegendMoved, Wt::JavaScriptFunction, "LegendMoved",
function(legID)
{
  var spectrum = $('#'+legID);
  
  if( !spectrum )
  {
    console.log("LegendMoved: !spectrum");
    return;
  }
  
  var legid = spectrum.data('legid');
  if( !legid )
  {
    console.log("LegendMoved: !legid");
    return;
  }
  
  var legend = $('#' + legid);
  if( !legid )
  {
    console.log("LegendMoved: !legend");
    return;
  }
  
  var spec_offset = spectrum.offset();
  var leg_offset  = legend.offset();
  var spec_width  = spectrum.outerWidth();
  var leg_width   = legend.outerWidth();
 
  if( !spec_offset || !leg_offset || !spec_width || !leg_width )
  {
    console.log("LegendMoved: !spec_offset || !leg_offset || !spec_width || !leg_width");
    return;
  }
  spectrum.data( 'LTM', leg_offset.top - spec_offset.top );
  spectrum.data( 'LRM', spec_width + spec_offset.left - leg_offset.left - leg_width );
}
);



WT_DECLARE_WT_MEMBER
(AlignLegend, Wt::JavaScriptFunction, "AlignLegend",
function( chartid )
{
  if( !chartid )
  {
    console.log( 'Couldnt get chartid=' + chartid );
    return;
  }

  var spec = $('#'+chartid);
  if( !spec )
  {
    console.log( 'Couldnt get spec' );
    return;
  }
  
  var id = spec.data('legid');
  
  if( !id )
    return;
  
  var leg = $('#' + id );
  if( !leg || leg.length === 0 )
  {
    console.log( 'Couldnt get legend with id ' + id );
    return;
  }

  try{
    var getOffset = function( el ){var _x = 0;var _y = 0;while( el && !isNaN( el.offsetLeft ) && !isNaN( el.offsetTop ) ){_x += el.offsetLeft;_y += el.offsetTop;el = el.offsetParent;}return { top: _y, left: _x };};

    //Move the legend to the left top position, to give it a chance to fully
    //  render and expand its width to fit its contents before we ask for its
    //  width and stuff
    leg.offset({ top: 0, left: 0, bottom: "", right: "" });
    
    var toppx=0, leftpx=0;
    var hiddenParent = Wt.WT.isHidden(spec.get(0));

    if( !hiddenParent )
    {
      var specOffset = getOffset(spec.get(0));
      toppx = specOffset.top + spec.data('LTM');
      leftpx = specOffset.left + spec.outerWidth()
               - leg.outerWidth()-spec.data('LRM');
    }

    leg.offset({ top: toppx, left: leftpx, bottom: "", right: "" });

    if( hiddenParent )
      leg.hide();
    else
      leg.show();
  }catch(e)
  {
    console.log("Failed in AlignLegend: " + e);
  };
}
 );



WT_DECLARE_WT_MEMBER
(AlignOverlay, Wt::JavaScriptFunction, "AlignOverlay",
function( childId, parentId )
{
 // Adapted from http://stackoverflow.com/questions/442404/dynamically-retrieve-html-element-x-y-position-with-javascript
  var getOffset = function( el )
  {
   var _x = 0;
   var _y = 0;
   while( el && !isNaN( el.offsetLeft ) && !isNaN( el.offsetTop ) )
   {
       //_x += el.offsetLeft - el.scrollLeft;
       _x += el.offsetLeft;
       //_y += el.offsetTop - el.scrollTop;
       _y += el.offsetTop;
      el = el.offsetParent;
    }

    return { top: _y, left: _x };
  };//function getOffset( el )


  try {
    var child = $('#' + childId );
    var parent = $('#' + parentId );
    var parentEl = parent.get(0);
    var childEl = child.get(0);
    var can = $('#c'+childId);
    var childCan = can.get(0);
    var parentCan = $('#c'+parentId).get(0);

    //console.log('parentId is ' + parentId + ' childId is ' + childId);
    var scrollParentEl = null;
    var scrollParent = null;
    var scrollParentId = child.data('scrollParent');
    if( !scrollParentId || scrollParentId === "" )
      scrollParentId = parent.data('scrollParent');

 // Note that scrollParentId is really the id for the "outer div", not the
 // id for the widget represented by the scrollParent C++ class

    if( scrollParentId !== "" ) {
      scrollParent = $('#' + scrollParentId );
      scrollParentEl = scrollParent.get(0);
    }

    // parent.offset().top and parent.offset().left are something
    // See http://api.jquery.com/offset/
    // Sometimes parent.offset().top is 0 - but if it's non zero,
    // and less than about 180 or so, then we'll want to clip the top
    // off of canvas.
    // console.log("parent.offset().top is " + parent.offset().top);

    // We only copy the parent's offset for
    // the x axis.  We don't scroll the canvas
    // about the y axis, because this might
    // overlap with the tabs.
    // Instead, we will later on compensate for
    // the y-scrolled amount in the C++ callback function
    // (SpectrumChart::handleDrag()

    if( scrollParent )
    {
      var newOffset = parent.offset();

      // Only do this for the Anthony app
      //Align new y position with the original (unscrolled) y position of the chart
      newOffset.top = getOffset( parentEl ).top;
      child.offset( newOffset );
    }else
    {
      console.log('default behavior for lower chart childId ' + childId);
      child.offset( parent.offset() );
    }

    var w, h;
    if( parentEl && parentEl.wtWidth )
      w = parentEl.wtWidth;
    else
      w = parent.width();

    //console.log('parentEl.scrollTop is ' + parentEl.scrollTop);
    //console.log('childEl.scrollTop is ' + childEl.scrollTop);
    //console.log('child id is ' + childId);
    //console.log('parent id is ' + parentId);
    //console.log('scrollParentId is ' + scrollParentId);

     // This only make sense for the Anthony app
    var outerDivScrollAmount = 0;
    if ( scrollParentId && scrollParentId !== ""){
      //console.log('scrollParentEl.scrollTop is ' + scrollParentEl.scrollTop);
      outerDivScrollAmount = scrollParentEl.scrollTop;
    }

    if( parentEl && parentEl.wtHeight ){
      // This is where you want to set the height accordingly, depending
      // on how much the user has scrolled
      h = parentEl.wtHeight - outerDivScrollAmount;
      //console.log('Setting h to parentEl.wtHeight=' + h);
    } else {
      h = parent.height() - outerDivScrollAmount;
      //console.log('Setting h to parent.height()=' + h);
    }

    childEl.style.width =  w+'px';
    childEl.style.height = h+'px';
    
    childCan.setAttribute('width',w);
    childCan.setAttribute('height',h);

 //XXX - for some reason when mouse is clicked with the control key pressed, this
 //  AlignOverlay function gets called, which then messes up the dragging info
 //  if 'startDragX' and 'startDragY' are set to null.  So as a workaround
 //  these varables are intialized to null (which is apparently necessary)
 //  in CanvasForDragging::loadJs(), and following two lines are commented out
//     can.data('startDragX',null);
//     can.data('startDragY',null);

    childEl.wtWidth = w;
    childEl.wtHeight = h;
    //console.log('Set childEl.wtWidth=' + w + ', childEl.wtHeight=' + h );
     
    if( childEl.wtResize ){
      childEl.wtResize( childEl, w, h );
      //console.log('Called wtResize(' + w + ',' + h + ')');
    }
  }catch(error){
    if( console && console.log )
      console.log("Failed in AlignOverlay");
  }
}
);


WT_DECLARE_WT_MEMBER
(DrawnLegendTextMetric, Wt::JavaScriptFunction, "DrawnLegendTextMetric",
function(id)
{
  try {
    var el = this.getElement('c'+id);
    var cntxt = el.getContext("2d");
    cntxt.font = '12pt normal';
//    var m = cntxt.measureText("Foreground (2.6e+06 counts)");
    var m = cntxt.measureText("Real Time 99999.9 s");
    Wt.emit(id, {name: 'LegendTextMetric'}, m.width);
  }catch(e){
    if(console && console.log )
    console.log("Failed in DrawnLegendTextMetric for '" + id + "': " +e);
  }
}
);

