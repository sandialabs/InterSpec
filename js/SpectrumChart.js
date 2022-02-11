/* Note: this is at the same time valid JavaScript and C++. */



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


WT_DECLARE_WT_MEMBER
(CanvasToPngDownload, Wt::JavaScriptFunction, "CanvasToPngDownload",
function(elid,filename)
{
  try{
    var canvas = $('#' + elid).find('canvas')[0];
    if( canvas.length === 0 )
      throw 'Could not find canvas';
      
    if( filename.length === 0 )
      filename = "spectrum.png"; //ToDo: add in date/time or something.
    var dt = canvas.toDataURL('image/png');
    dt = dt.replace(/^data:image\\/[^;]*/, 'data:application/octet-stream');
    dt = dt.replace(/^data:application\\/octet-stream/, 'data:application/octet-stream;headers=Content-Disposition%3A%20attachment%3B%20filename='+filename);
    
    var link = document.createElement("a");
    link.download = filename;
    link.href = dt;
    link.target="_blank";
    link.click();
  }catch(e){
    console.log( 'Error saving PNG from spectrum: ' + e );
  }
}
);

