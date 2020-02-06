
DecayChainChart = function(elem, options) {
  var self = this;
  
  this.area = typeof elem === 'string' ? document.getElementById(elem) : elem;
  
  
  this.options = options || {};
  
  if( (typeof this.options.bottomPad) !== 'number' ) this.options.bottomPad = 0;
  if( (typeof this.options.leftPad) !== 'number' ) this.options.leftPad = 0;
  if( (typeof this.options.topPad) !== 'number' ) this.options.topPad = 0;
  if( (typeof this.options.rightPad) !== 'number' ) this.options.rightPad = 0;

  this.svg = d3.select(this.area).append("svg");
  
  this.chartbackground = this.svg.append("rect")
  .attr("id", "chartbackground"+this.area.id )
  .attr("class", "chartbackground" )
  .attr("clip-path", "url(#clip" + this.area.id + ")")
  .style("fill", "#ff0000")
  ;
  /*.attr("pointer-events", "all"); */
  
  this.chart = this.svg.append("g");
  this.decays = this.chart.append("g");
  this.xaxis = this.chart.append("g");
  this.yaxis = this.chart.append("g");
  
  this.clipPath = this.chart.append("svg:clipPath")
  .attr("id", "clip" + this.area.id )
  .append("svg:rect")
  .attr("x", 0)
  .attr("y", 0);

  const xaxislabel = '<tspan dy="-0.3em">Decreasing Atomic Mass</tspan>';
  const yaxislabel = '<tspan dy="1.0em">Increasing Atomic Number</tspan>';
  self.xaxistxt = self.xaxis.append("text").html(xaxislabel);
  self.yaxistxt = self.yaxis.append("text")
                      .html(yaxislabel)
                      .attr("transform", "rotate(-90)" )
                      .attr("y", 0);
                      
  self.xaxisline = self.xaxis.append("line")
    .attr("stroke-width", 2)
    .attr("stroke", "black");
  
  self.yaxisline = self.yaxis.append("line")
  .attr("stroke-width", 2)
  .attr("stroke", "black");
  
  /* this.svg = d3.select(self.chart).select('svg'); */
  
  /*
   Chart vs. Vis vs. Document Body
   ___________________________________________
   | _______________________________________ |
   ||                                       ||
   ||             decays <g>                ||
   ||                                       ||
   ||                                       ||
   ||_______________________________________||
   |                                         |
   |                svg <svg>                |
   |_________________________________________|
   */
  
  
  /*
  this.yAxisBody = this.decays.append("g")
  .attr("transform", "translate(0,0)");
  
  this.xAxisBody = this.decays.append("g")
  .attr("transform", "translate(0," + this.size.height + ")")
  ;
  */
  
  
  /* File out member variables (with dummy values) some functions expect to exist */
  this.data = [];
  this.fontSize = this.fontHeight = 8;
  
  //ToDo: add axis elements here, and just resize/position in redraw
  //this.decays.xaxis
  
  
  this.handleResize();
}


DecayChainChart.prototype.handleResize = function() {
  
  console.log( 'DecayChainChart.prototype.handleResize' );

  this.redraw();
}

DecayChainChart.prototype.setDecayData = function( data ){
  //console.log( 'data=', data );
  
  try
  {
    data = JSON.parse(data);
  }catch(e)
  {
    console.log( 'DecayChainChart.setDecayData: Error parsing data ', e );
    data = [];
  }
  
  if( !Array.isArray(data) )
    data = [];
  
  //Massage data so that each combination of AN and MN are only in the data once
  //  no matter how many isomeric states there are
  let filteredData = [], multiIsoData = [];
  
  //First, add the lowest isomeric state for each combination of {AN,MN}
  data.forEach( function(nuc){
    //Get number of entries with this AN and MN
    let numLessIso = data.reduce( function(acc, n){
      return acc + (n.atomicNumber==nuc.atomicNumber && n.massNumber==nuc.massNumber && n.iso < nuc.iso ? 1 : 0);
     }, 0);
    
    if( numLessIso == 0 )
      filteredData.push( nuc );
    else
      multiIsoData.push( nuc );
  });
  
  //Now, add the higher isomeric states to the entry already in filteredData
  multiIsoData.forEach( function(nuc){
    let val = filteredData.find( function(n){
      return (n.atomicNumber==nuc.atomicNumber && n.massNumber==nuc.massNumber);
    });
    if( val == undefined ){
      console.log( 'Unexpected error searching filteredData');
      return;
    }
    if( !val.additionalIsos )
      val.additionalIsos = [];
    val.additionalIsos.push( nuc );
  } );
  
  
  this.data = filteredData;
  this.redraw();
}


DecayChainChart.prototype.redraw = function() {
  let self = this;
  
  // cx and cy are total width and height of this widget
  const cx = this.area.clientWidth;
  const cy = this.area.clientHeight;
  
  this.svg.attr("width", cx);
  this.svg.attr("height", cy);
  this.chartbackground.attr("width", cx).attr("height", cy);
  
  //chartwidth and chartheight are the area of nuclide boxes plus the axis guides
  const chartwidth = Math.max(0, cx - this.options.leftPad - this.options.rightPad);
  const chartheight = Math.max(0, cy - this.options.topPad  - this.options.bottomPad);
  
  this.clipPath.attr("width", chartwidth).attr("height", chartheight);
  
  this.chart
      .attr("transform", "translate(" + this.options.leftPad + "," + this.options.topPad + ")");
  
  
  let atomicMasses = [], massNumbers = [];
  
  //Figure out AM and MN we will eventually plot
  this.data.forEach( function(nuc){
    if( atomicMasses.indexOf(nuc.massNumber) < 0 )
    atomicMasses.push( nuc.massNumber );
    if( massNumbers.indexOf(nuc.atomicNumber) < 0 )
    massNumbers.push( nuc.atomicNumber );
  } );
  
  atomicMasses.sort(function(a, b){return b-a});
  massNumbers.sort(function(a, b){return b-a});
  
  
  //Atomic mass will be on x-axis, and atomic number on y-axis
  

  //console.log( 'atomicMasses=', atomicMasses );
  //console.log( 'massNumbers=', massNumbers );
  
  //this.decays.attr("transform", "translate(" + yaxiswidth + ",0)");
  
  const w = chartwidth;
  const h = chartheight;
  const num_am_x = atomicMasses.length + 0.5;  //assume y-axis will take up half the width of an nuclide
  const num_mn_y = massNumbers.length + 1;  //assume x-axis will take up same height as a nuclide
  //console.log( 'w=', w, ', h=', h, ', num_am_x=', num_am_x, ', num_mn_y=', num_mn_y );
  
  //How to size everything is a work in progress - to say the least
  
  //We will draw a box for each nuclide, and place the nuclides text (ex.
  //  "U238") centered in each box.  However, we want the sizing to be adaptive,
  //  both how large the boxes/text are, and spacing between boxes.
  
  //Check how much space text will take up, and then assume this scales
  //  lineraly with font size (ToDo: is this reasonable?)
  //Right now we'll just loop over all elements and find the largest one; this
  //  can be improved if performance is an issue.
  let twenty_px_w = 1, twenty_px_h = 1;
  let labeltxt = function(nuc){
    /* return nuc.additionalIsos ? nuc.nuclide + "(m)" : nuc.nuclide; */
     let s = nuc.nuclide;
     if( !nuc.additionalIsos )
      return s;
    if( nuc.iso != 0 && s.includes('m') )
      s = s.substring(0,s.lastIndexOf('m'));
    return s + '<tspan dy="-0.5em" font-size="66%">(m)</tspan>';
  };
  
  this.data.forEach( function(nuc){
    let tmptxt = self.decays
                     .append("text")
                     .html( labeltxt(nuc) )
                     //.text( nuc.nuclide )
                     .attr("font-size", "20px");
    twenty_px_w = Math.max( twenty_px_w, tmptxt.node().getBoundingClientRect().width );
    twenty_px_h = Math.max( twenty_px_h, tmptxt.node().getBoundingClientRect().height );
    tmptxt.remove();
  } );
  
  //console.log( 'twenty_px_w=', twenty_px_w, ', twenty_px_h=', twenty_px_h );
  
  //We will use a minimum txt size of 5px, and a maximum of 105px.
  //At 20px we want the spacing to be {}
  
  //Define how much padding text should have inside each box.
  const txt_frac_x = 0.9, txt_frac_y = 0.8;
  const twenty_px_box_w = twenty_px_w / txt_frac_x;
  const twenty_px_box_h = twenty_px_h / txt_frac_y;
  
  //Now figure out what the x-sizing needs to be, and the y-sizing, and use
  //  the minimum of those to set the sizing
  
  //At a text size of 5px, we will have zero spacing between boxes, these next
  //  variables specify the fraction of box spacing between nuclides by the
  //  time we get to 105px size.
  let x_pad_mult_upper = 0.25;
  let y_pad_mult_upper = 0.25;
  
  const min_x_pad = 8, min_y_pad = 8;
  
  const min_font_size = 6;
  const min_font_frac = min_font_size / 20;
  const max_font_frac = 105 / 20;
  
  //We will vary the x spacing from 0 at 5px, to x_pad_mult_upper box height at 105px.
  //We will vary the y spacing from 0 at 5px, to y_pad_mult_upper box height at 105px.
  const width_fivepx = (min_x_pad + min_font_frac*twenty_px_box_w) * num_am_x;
  const width_onefivepx = (twenty_px_box_w * max_font_frac)*(num_am_x + x_pad_mult_upper*(num_am_x+1));
  const height_fivepx = (min_y_pad + min_font_frac*twenty_px_box_h) * num_mn_y;
  const height_onefivepx = (twenty_px_box_h * max_font_frac)*(num_mn_y + y_pad_mult_upper*(num_mn_y+1));
  
  //console.log( 'width_fivepx=', width_fivepx, ', width_onefivepx=', width_onefivepx );
  //console.log( 'height_fivepx=', height_fivepx, ', height_onefivepx=', height_onefivepx );
  
  //Get how much between 5px and 105px widths our current width is (may be negative or >1)
  let x_frac_between = (w - width_fivepx)/(width_onefivepx - width_fivepx);
  let y_frac_between = (h - height_fivepx)/(height_onefivepx - height_fivepx);
  
  //console.log( 'x_frac_between=', x_frac_between, ', y_frac_between=', y_frac_between );
  
  x_frac_between = Math.min( Math.max(x_frac_between,0), 1);
  y_frac_between = Math.min( Math.max(y_frac_between,0), 1);
  
  const font_size = min_font_size + 100*Math.min(x_frac_between,y_frac_between);
  const font_height = twenty_px_h * font_size / 20;
  
  this.fontSize = font_size;
  this.fontHeight = font_height;
  
  const box_w = twenty_px_box_w *(font_size/20);
  const box_h = twenty_px_box_h *(font_size/20);
  const pad_x = (w - num_am_x*box_w) / num_am_x;
  const pad_y = (h - num_mn_y*box_h) / num_mn_y;
  
  //console.log( 'font_size=', font_size );
  
  let rects = this.decays.selectAll("rect").data( this.data );
  rects.enter().append("rect");
  
  let rectXPos = function(d){
    const ind = atomicMasses.indexOf(d.massNumber);
    return (0.5 + ind)*pad_x + ind*box_w + 0.5*box_w; //the extra 0.5*box_w is room for y-axis
  };
  
  let rectYPos = function(d){
    const ind = massNumbers.indexOf( d.atomicNumber );
    return (0.5 + ind)*pad_y + ind*box_h;
  };
  
  rects.attr("x", rectXPos )
  .attr("y", rectYPos )
  .attr("width", box_w )
  .attr("height", box_h )
  .attr("rx", 0.08*box_w )
  .attr("stroke", "black")
  //.attr("fill", "none")
  .attr("fill", "#ffffff")  //add fill so mouse over will be active for entire rect area
  .attr("fill-opacity", "0")
  .attr("class", function (d) {
     return "DecayChartNuc "
            + d.nuclide
            + (d.isPrimaryNuc ? " primary" : "")
    } )
  .on("mouseover", function(d) {
    d3.selectAll('.DecayChartNuc').classed('MousedOver',false);
    d3.selectAll('.'+d.nuclide).classed('MousedOver',true);
  })
  .on("mouseout",  function(d) {
    d3.selectAll('.DecayChartNuc').classed('MousedOver',false);
  } )
  .on("click", self.nuclideClicked)
  .on("dblclick", function(d){ console.log( 'dblclick', d ); } ) //we need to call back to c++ here
  ;
  
  rects.exit().remove();
  
  
  let texts = this.decays.selectAll("text").data( this.data );
  
  texts.enter().append("text");
  
  let txtXPos = function(d){
    const ind = atomicMasses.indexOf(d.massNumber);
    return (0.5 + ind)*pad_x + (ind+0.5)*box_w + 0.5*box_w;  //the extra 0.5*box_w is room for y-axis
  };
  
  let txtYPos = function(d){
    const ind = massNumbers.indexOf( d.atomicNumber );
    return (0.5 + ind)*pad_y + ind*box_h + 0.35*font_height + 0.5*box_h;
  };
  
  texts.attr("x", txtXPos )
       .attr("y", txtYPos )
       .attr("class", function (d) {
         return "DecayChartNuc "
                 + d.nuclide
                 + (d.isPrimaryNuc ? " primary" : "")
                 + (d.additionalIsos ? " MultiIso" : "");
        } )
       .attr("style","text-anchor: middle")
       .on("mouseover", function(d) {
         d3.selectAll('.DecayChartNuc').classed('MousedOver',false);
         d3.selectAll('.'+d.nuclide).classed('MousedOver',true);
       })
       .on("mouseout",  function(d) {
         d3.selectAll('.DecayChartNuc').classed('MousedOver',false);
       } )
       .on("click", self.nuclideClicked)
       .html( labeltxt )
       .attr("font-family", "sans-serif")
       .attr("font-size", function(){ return font_size + "px";} )
       .attr("fill", "black");

  texts.exit().remove();
  
  
  //Set axis font size and reposition.
  self.yaxistxt.attr("font-size", Math.min(1.5*font_size,105));
  self.xaxistxt.attr("font-size", Math.min(1.5*font_size,105));
  
  const xtxth = self.xaxistxt.node().getBoundingClientRect().height;
  const xtxtw = self.xaxistxt.node().getBoundingClientRect().width;
  const ytxth = self.yaxistxt.node().getBoundingClientRect().height;
  const ytxtw = self.yaxistxt.node().getBoundingClientRect().width;
  
  self.xaxistxt.attr("x", 2*ytxtw )
      .attr("y", chartheight );
  self.yaxistxt.attr("x", -chartheight + 2*xtxth );
  
  const xaxisy = chartheight - 1.2*xtxth;
  const strokew = Math.min( Math.max(font_size/10, 1), 2.25 );
  
  self.xaxisline
    .attr("x1", 1.2*ytxtw)
    .attr("y1", xaxisy)
    .attr("x2", 1.2*ytxtw + 1.3*xtxtw)
    .attr("y2", xaxisy)
    .attr("stroke-width", strokew);
  
  self.yaxisline
    .attr("x1", 1.2*ytxtw)
    .attr("y1", xaxisy)
    .attr("x2", 1.2*ytxtw)
    .attr("y2", xaxisy - 1.15*ytxth)
    .attr("stroke-width", strokew);
  
  
}//redraw

DecayChainChart.prototype.nuclideClicked = function( d ) {
  console.log( 'nuclideClicked: ' + d.nuclide );
  //d is { nuclide: "Pa230", massNumber: 230, atomicNumber: 91, iso: 0, halfLive: "17.40 d", … }
  d3.selectAll('rect.DecayChartNuc').classed('selected',false);
  d3.selectAll('rect.DecayChartNuc.' + d.nuclide).classed('selected',true);
  //display text
}

