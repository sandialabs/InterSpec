/* Feature TODO list (created 20160220):
  -Fix intermitten issue of zooming in messes up (especially aver dragging starting from the y-axis title)
  [Done 20160218] -Add legend
  [Done 20160218: Christian] -Zoom in on y-axis
  [Done 20160815: Liza] -Draw Peaks
  -Add peak info display
  -Make it so x-axis binning is given seperately for each histogram
  [Done 20160310] -Add reference photopeak line drawing
  [Done 20160310] -Box that gives stats that follow mouse
  [Done 20160520] -x-axis tick labels can overlap (especially for higher enegy values zoomed in for a small range)
  [Done 20160520] -Resize plot area based on y-axis label sizes (but if actively zooming in with mouse wheel, keep centered at mouse)
  [Done 20160218] -scroll-wheel over y-axis adds remove default amount of padding above/bellow charts
  -scroll left-right (track pad two fingers) pan left/right
  [Done 20160218: Christian] -control (alt on mac) drag create duplicate charts that drag with mouse (for energy recalibration)
  -Add in a mock Wt.WT.emit(...) function and hook up signals
  [Done 20160218: Christian]Animate the click-and-drag-right zoom in effect
  -Some sliders to change scale factor on background and second spectrum (should appear only when something is clicked on to show them)
  -Background subtraction mode
  -Add statistical error bars
  [Done 20160815: Christian] -Add escape peak, compon scatter, etc. mouse interactions
  [Done 20160815: Christian] Add arrows to end of x axis to indicate if the chart extends to the right and left (e.g. can pan)
  [Done 20160521] -Add compact x-axis mode
  -Add widget to scale background and second spectrum
  [Done 20170318: Christian] -Add in (optional) x-axis slider widget (the one that shows the entire chart, and indicates what portion you are viewing and lets you drag to change this)
  [Done 20170218: Christian] -Add in InterSpec specific mouse functionality (drawing line for peak, deleting peaks, etc)
  [Done 20161216: Christian] -Add in ability to add in text to label features (very low priority)
  -Add ability to support polynomial binning with deviation pairs (instead of just lower edge energy)
  -Customize strokes, colors, lines, etc.
  [Done 20160520] -Dont let zoom in to less than one bin width
  [Done 20160521] -Fix x-axis labels that hang off end of plot area
  [Done 20160520] -Fix adding chart title and compensating for x-axis title
  -Customize mouse point to zoom-in/zoom-out where appropriate
  -Optimize frequency of rebinning of data (prevent extra rebinned data from being drawn)
  -Need some way to filter reference gamma lines to not draw insignificant lines.  Ex, Th232 gives ~900 dom elements, which can slow things down
*/

SpectrumChartD3 = function(elemid, options) {
  var self = this;
  this.chart = document.getElementById(elemid);

  this.cx = this.chart.clientWidth;
  this.cy = this.chart.clientHeight;
  this.options = options || {};

  this.options.yscale = "lin";
  this.options.gridx = false;
  this.options.gridy = false;
  this.options.compactXAxis = false; 
  this.options.adjustYAxisPadding = true;
  this.options.wheelScrollYAxis = true;
  this.options.showAnimation = false;
  this.options.showXAxisSliderChart = false;

  this.options.showUserLabels = (typeof options.showUserLabels == 'boolean') ? options.showUserLabels : false;
  this.options.showPeakLabels = (typeof options.showPeakLabels == 'boolean') ? options.showPeakLabels : false;
  this.options.showNuclideNames = (typeof options.showNuclideNames == 'boolean') ? options.showNuclideNames : false;
  this.options.showNuclideEnergies = (typeof options.showNuclideEnergies == 'boolean') ? options.showNuclideEnergies : false;
  
  this.options.showLegend = true;
  this.options.scaleBackgroundSecondary = (typeof options.scaleBackgroundSecondary == 'boolean') ? options.scaleBackgroundSecondary : false;

  this.options.refLineTopPad = 30;
  
  self.options.logYFracTop = 0.05;
  self.options.logYFracBottom = 0.025;
  self.options.linYFracTop = 0.1;
  self.options.linYFracBottom = 0.1;
  self.options.sqrtYFracTop = 0.1;
  self.options.sqrtYFracBottom = 0.1;

  this.options.showMouseStats = (typeof options.showMouseStats == 'boolean') ? options.showMouseStats : true;
  this.options.showComptonEdge = (typeof options.showComptonEdge == 'boolean') ? options.showComptonEdge : false;
  this.options.showComptonPeaks = (typeof options.showComptonPeaks == 'boolean') ? options.showComptonPeaks : false;
  this.options.showEscapePeaks = (typeof options.showEscapePeaks == 'boolean') ? options.showEscapePeaks : false;
  this.options.showSumPeaks = (typeof options.showSumPeaks == 'boolean') ? options.showSumPeaks : false;
  this.options.drawPeaks = (('drawPeaks' in options) && options.drawPeaks.constructor == Array) ? options.drawPeaks : [true,true,true];

  this.options.showXRangeArrows = (typeof options.showXRangeArrows == 'boolean') ? options.showXRangeArrows : true;

  this.options.maxScaleFactor = 10;

  this.padding = {
     "top":  5,
     "titlePad" : 5,
     "right":   10,
     "bottom": 5,
     "xTitlePad": 5,
     "left":     5,
     "labelPad": 5,
     "title":    23,
     "label":    8,
     "sliderChart":    8,
  };
  
  this.padding.leftComputed = this.padding.left + this.padding.title + this.padding.label + this.padding.labelPad;
  this.padding.topComputed = this.padding.top + this.padding.titlePad + 15;
  this.padding.bottomComputed = this.padding.bottom +  this.padding.xTitlePad + 15;
  
  this.size = {
    "width":  this.cx - this.padding.leftComputed - this.padding.right,
    "height": this.cy - this.padding.topComputed  - this.padding.bottomComputed,
    "sliderChartHeight": (this.cy - this.padding.topComputed  - this.padding.bottomComputed) / 10,
    "sliderChartWidth": this.cx - this.padding.leftComputed - this.padding.right - 30,
  };

  this.animation = {
    "duration": 200,
    "frames": 6
  };

  //pre IE9 support
  if (!Array.isArray) {
    Array.isArray = function(arg) {
      return Object.prototype.toString.call(arg) === '[object Array]';
    };
  }

  //When dragging the plot, both dragging_plot and zooming_plot will be
  //  true.  When only zooming (e.g. mouse wheel), then only zooming_plot
  //  will be true.
  this.dragging_plot = false;
  this.zooming_plot = false;

  this.refLines = [];

  // x-scale
  this.xScale = d3.scale.linear()
      .domain(this.options.xScaleDomain ? this.options.xScaleDomain : [0, 3000])
      .range([0, this.size.width]);

  // drag x-axis logic
  this.downx = Math.NaN;

  // y-scale (inverted domain)
  this.yScale = d3.scale.linear()
      .domain([0, 100])
      .nice()
      .range([0, this.size.height])
      .nice();

  this.downy = Math.NaN;

  // Finds distance between two points
  this.dist = function (a, b) {
    return Math.sqrt(Math.pow(a[0]-b[0],2) + Math.pow(a[1]-b[1],2));
  };


  this.getCountsForEnergy = function(spectrum, energy) {
    if (!self.rawData || !self.rawData.spectra || !spectrum || !spectrum.x)
      return -1;

    var channel, lowerchanval, counts = null;
    var spectrumIndex = self.rawData.spectra.indexOf(spectrum);

    if (spectrumIndex < 0)
      return -1;

    channel = d3.bisector(function(d){return d;}).right(spectrum.x, energy);

    if( spectrum.points && spectrum.points.length ){
      lowerchanval = d3.bisector(function(d){return d.x;}).left(spectrum.points,energy,1) - 1;
      counts = spectrum.points[lowerchanval].y;
    }
    return counts;
  };

  this.min_max_x_values = function() {
    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
      return [-1,-1];

    var min = null, max = null;
    self.rawData.spectra.forEach(function(spectrum) {
      if (!spectrum.x)
        return;

      if (min == null || spectrum.x[0] < min) min = spectrum.x[0];
      if (max == null || spectrum.x[spectrum.x.length-1] > max) max = spectrum.x[spectrum.x.length-1];
    });

    return [min,max];
  };

  this.displayed_raw_start = function(spectrum){
    if( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length )
      return -1;
    var xstart = self.xScale.domain()[0];
    if (!spectrum)
      spectrum = self.rawData.spectra[0]; // use foreground by default

   // switch to using:
   // var bisector = d3.bisector(function(d){return d.x;});
   // bisector.left(spectrum.x, xstart)

    var i = 0;
    while( i < spectrum.x.length && spectrum.x[i] < xstart )
      ++i;
    return i;
  };

  this.displayed_raw_end = function(spectrum){
    if( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
      return -1;
    var xend = self.xScale.domain()[1];
    if (!spectrum)
      spectrum = self.rawData.spectra[0]; // use foreground by default
    var i = spectrum.x.length - 1;
    while( i > 0 && spectrum.x[i] > xend )
      --i;
    return i + 1;
  };

  this.displayed_start = function(spectrum){
    if( !spectrum || !spectrum.points || !spectrum.points.length )
      return -1;
    var xstart = self.xScale.domain()[0];
    var i = 0;
    while( i < spectrum.points.length && spectrum.points[i].x <= xstart )
      ++i;
    return Math.max(i - 1,0);
  };

  this.displayed_end = function(spectrum){
    if( !spectrum || !spectrum.points || !spectrum.points.length )
      return -1;
    var xend = self.xScale.domain()[1];
    var i = spectrum.points.length - 1;
    while( i > 0 && spectrum.points[i].x >= xend )
      --i;
    return Math.min(i + 1, spectrum.points.length);
  };

  this.firstRaw = this.lastRaw = this.rebinFactor = -1;
  this.do_rebin();

  this.getYAxisDomain = function(){
    if (!this.rawData || !this.rawData.spectra || !this.rawData.spectra.length)
      return [3000,0];
    var self = this;
    var y0, y1;
    var minx = self.xScale.domain()[0], maxx = self.xScale.domain()[1];
    var foreground = this.rawData.spectra[0];
    var firstData = this.displayed_start(foreground);
    var lastData = this.displayed_end(foreground);

    if( firstData >= 0 ){
      y0 = y1 = foreground.points[firstData].y;

      this.rawData.spectra.forEach(function(spectrum) {
        firstData = self.displayed_start(spectrum);
        lastData = self.displayed_end(spectrum);

        for (var i = firstData; i < lastData; i++) {
          y0 = Math.min( y0, spectrum.points[i].y );
          y1 = Math.max( y1, spectrum.points[i].y );
        }
      });
    }else {
      y0 = 0;
      y1 = 3000;
    }

    if( y0 > y1 ) { y1 = [y0, y0 = y1][0]; }
    if( y0 == y1 ){ y0 -=1; y1 += 1; }


    if( this.options.yscale == "log" ) {
      //Specify the (approx) fraction of the chart that the scale should extend
      //  past where the data where hit.
      var yfractop = self.options.logYFracTop, yfracbottom = self.options.logYFracBottom;

      var y0Intitial = ((y0<=0.0) ? 0.1 : y0);
      var y1Intitial = ((y1<=0.0) ? 1.0 : y1);
      y1Intitial = ((y1Intitial<=y0Intitial) ? 1.1*y0Intitial : y1Intitial);

      var logY0 = Math.log10(y0Intitial);
      var logY1 = Math.log10(y1Intitial);

      var logLowerY = ((y0<=0.0) ? -1.0 : (logY0 - yfracbottom*(logY1-logY0)));
      var logUpperY = logY1 + yfractop*(logY1-logY0);

      var ylower = Math.pow( 10.0, logLowerY );
      var yupper = Math.pow( 10.0, logUpperY );

      y0 = ((y0<=0.0) ? 0.1 : ylower);
      y1 = ((y1<=0.0) ? 1.0 : yupper);
    } else if( this.options.yscale == "lin" )  {
      y0 = ((y0 <= 0.0) ? (1+self.options.linYFracBottom)*y0 : (1-self.options.linYFracBottom)*y0);
      y1 = (1 + self.options.linYFracTop)*y1;
    } else if( this.options.yscale == "sqrt" ) {
      y0 = ((y0 <= 0.0) ? 0.0 : (1-self.options.sqrtYFracBottom)*y0);
      y1 = (1+self.options.sqrtYFracTop)*y1;
    }

    return [y1,y0];
  };

  this.setYAxisDomain = function(){
    if( !isNaN(self.downy) )
      return;
    var yaxisDomain = self.getYAxisDomain(),
        y1 = yaxisDomain[0],
        y0 = yaxisDomain[1];
    this.yScale.domain([y1, y0]);
  };
  this.setYAxisDomain();


  // drag y-axis logic
  this.downy = Math.NaN;

  this.dragged = this.selected = null;


  this.vis = d3.select(this.chart).append("svg")
      .attr("width",  this.cx)
      .attr("height", this.cy)
      .append("g")
      .attr("transform", "translate(" + this.padding.leftComputed + "," + this.padding.topComputed + ")");

  this.plot = this.vis.append("rect")
      .attr("width", this.size.width)
      .attr("height", this.size.height)
      .attr("id", "chartarea"+this.chart.id )
      .style("fill", "#EEEEEE");
      //.attr("pointer-events", "all");

  /*
  Chart vs. Vis vs. Document Body
  ___________________________________________
  | _______________________________________ |
  ||  ___________________________________  ||
  || |                                   | ||
  || |           VIS (SPECTRUM)          | ||
  || |                                   | ||
  || |___________________________________| ||
  ||                                       ||
  ||                  CHART                ||
  ||_______________________________________||
  |                                         |
  |                                         |
  |                   BODY                  |
  |                                         |
  |                                         |
  |_________________________________________|
  */

  /* 
  NOTE: Restoring the zoom behavior only for touch events
        This does NOT affect the behavior for other mouse events.
        Restoring the zoom behavior allows pinch zooming/panning with touch interactions.
  */
  this.zoom = d3.behavior.zoom()
    .x(self.xScale)
    .y(self.yScale)
    .on("zoom", self.handleZoom())
    .on("zoomend", self.handleZoomEnd());

  /* Vis interactions */
  this.vis
    .call(this.zoom)
    .on("mousedown", self.handleVisMouseDown())
    .on("mouseup", self.handleVisMouseUp())
    .on("wheel", self.handleVisWheel())
    .on("touchstart", self.handleVisTouchStart())
    .on("touchmove", self.handleVisTouchMove())
    .on("touchend", self.handleVisTouchEnd());

  // Cancel the zooms for mouse events, we'll use our own implementation for these (however, keeping the touch zooming)
  this.vis
    .on("mousedown.zoom", null)
    .on("mousemove.zoom", null)
    .on("mouseup.zoom", null)
    .on("mouseover.zoom", null)
    .on("mouseout.zoom", null)
    .on("wheel.zoom", null)
    .on("click.zoom", null)
    .on("dblclick.zoom", null);

  d3.select(document.body)
    .on("mouseup", self.handleCancelAllMouseEvents());
  d3.select(window).on("mouseup", self.handleCancelAllMouseEvents())
    .on("mousemove", function() { if (self.sliderBoxDown || self.leftDragRegionDown || self.rightDragRegionDown || self.currentlyAdjustingSpectrumScale) { d3.event.preventDefault(); d3.event.stopPropagation(); }});



  /*  Chart Interactions */
  d3.select(this.chart)
    .on("mousemove", self.handleChartMouseMove())
    .on("mouseout", self.handleChartMouseOut())
    .on("mouseup", self.handleChartMouseUp())
    .on("wheel", self.handleChartWheel())
    .on("touchstart", self.handleChartTouchStart())
    .on("touchend", self.handleChartTouchEnd());

  /*
  To allow markers to be updated while mouse is outside of the chart, but still inside the visual.
  */
  this.yAxis = d3.svg.axis().scale(this.yScale)
   .orient("left")
   .innerTickSize(7)
   .outerTickSize(1)
   .ticks(0);

  this.yAxisBody = this.vis.append("g")
    .attr("class", "yaxis" )
    .attr("transform", "translate(0,0)")
    .call(this.yAxis);

  this.xAxis = d3.svg.axis().scale(this.xScale)
   .orient("bottom")
   .innerTickSize(7)
   .outerTickSize(1)
   .ticks(20, "f");

  this.xAxisBody = this.vis.append("g")
    .attr("class", "xaxis" )
    .attr("transform", "translate(0," + this.size.height + ")")
    .call(this.xAxis);

  this.vis.append("svg:clipPath")
    .attr("id", "clip" + this.chart.id )
    .append("svg:rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", this.size.width )
    .attr("height", this.size.height );

  self.addMouseInfoBox();

  self.peakVis = this.vis.append("g")
    .attr("class", "peakVis")
    .attr("transform","translate(0,0)")
    .attr("clip-path", "url(#clip" + this.chart.id + ")");

  //Make a <g> element to draw everything we want that follows the mouse around
  //  when we're displaying reference photopeaks.  If the mouse isnt close enough
  //  to a reference line, then this whole <g> will be hidden
  self.refLineInfo = this.vis.append("g")
    .attr("class", "refLineInfo")
    .style("display", "none");

  //Put the reference photopeak line text in its own <g> element so
  //  we can call getBBox() on it to get its extent to to decide where to
  //  position the text relative to the selected phtotopeak.
  self.refLineInfoTxt = self.refLineInfo.append("g");

  //Add the text to the <g>.  We will use tspan's to append each line of information
  self.refLineInfoTxt.append("text")
   .attr("x", 0)
   .attr("dy", "1em");

  //Add a small red circle on the x axis to help indicate which line the info is
  //  currently showing for.
  self.refLineInfo.append("circle")
    .attr("cx", 0)
    .attr("cy", self.size.height)
    .attr("r", 2)
    .style("fill","red");

  this.chartBody = this.vis.append("g")
    .attr("clip-path", "url(#clip" + this.chart.id + ")");


  // add Chart Title
  var title = this.options.title;
  this.options.title = null;
  this.setTitle( title, true );

  // Add the x-axis label
  if (this.options.xlabel) {
    this.xaxistitle = d3.select('svg').append("text")
        .attr("class", "xaxistitle")
        .text(this.options.xlabel)
        .attr("x", this.size.width/2)
        .attr("y", this.size.height )
        .attr("dy",29)
        .style("text-anchor","middle");
  }

  // Add y-axis label
  if (this.options.ylabel) {
    this.vis.append("g").append("text")
        .attr("class", "yaxistitle")
        .text(this.options.ylabel + (this.rebinFactor > 1 ? (" per " + this.rebinFactor + " Channels") : "") )
        .style("text-anchor","middle");
        //.attr("transform","translate(" + -40 + " " + this.size.height/2+") rotate(-90)");
  }
}

registerKeyboardHandler = function(callback) {
  var callback = callback;
  d3.select(window).on("keydown", callback);
}

//
// SpectrumChartD3 methods
//
SpectrumChartD3.prototype.getYScaleFactors = function() {
  var self = this;

  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  var result = [];
  for (var i = 0; i < self.rawData.spectra.length; ++i)
    result.push(self.rawData.spectra[i].yScaleFactor);

  // Return a copy of the y-scale factors
  return result.slice();
}

SpectrumChartD3.prototype.getSpectrumTitles = function() {
  var self = this;

  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  var result = [];
  self.rawData.spectra.forEach(function(spectrum) {
    if (spectrum.title)
      result.push(spectrum.title);
  });
  return result;
}

SpectrumChartD3.prototype.do_rebin = function() {
  var self = this;

  if( !this.rawData || !self.rawData.spectra || !self.rawData.spectra.length ) {
    return;
  }

  var foreground = this.rawData.spectra[0];

  this.rawData.spectra.forEach(function(spectrum, spectrumi) {
    var newRebin = 1;
    spectrum.points = [];

    var firstRaw = self.displayed_raw_start(spectrum);
    var lastRaw = self.displayed_raw_end(spectrum);

    var npoints = lastRaw - firstRaw;
    if( npoints > 1 && self.size.width > 2 ) {
      newRebin = Math.ceil( npoints / (self.size.width) );
    }

    if( newRebin != spectrum.rebinFactor || self.firstRaw !== firstRaw || self.lastRaw !== lastRaw ){
      spectrum.points = [];

      if( spectrum.rebinFactor != newRebin ){
        var txt = self.options.ylabel ? self.options.ylabel : "";
        if( newRebin !== 1 )
          txt += " per " + newRebin + " Channels"
        d3.select(".yaxistitle").text( txt );
      }

      spectrum.rebinFactor = newRebin;
      spectrum.firstRaw = firstRaw;
      spectrum.lastRaw = lastRaw;

      //Round firstRaw and lastRaw down and up to even multiples of newRebin
      firstRaw -= (firstRaw % newRebin);
      lastRaw += newRebin - (lastRaw % newRebin);
      if( firstRaw >= newRebin )
        firstRaw -= newRebin;
      if( lastRaw > spectrum.x.length )
        lastRaw = spectrum.x.length;


      //could do some optimizations here where we actually do a slightly larger
      //  range than displayed, so that next time we might not have to go back
      //  through the data to recompute things (it isnt clear if D3 will plot
      //  these datas, should check on this.)
      //Also, could _maybe_ use energy range, rather than indexes to track if we
      //  need to rebin the data or not...

      //console.log( "self.rebinFactor=" + self.rebinFactor );

      for( var i = firstRaw; i < lastRaw; i += newRebin ){
        var thisdata = { };
        if (i >= spectrum.x.length) {
          thisdata['x'] = spectrum.x[spectrum.x.length-1];
        }
        else
          thisdata['x'] = spectrum.x[i];

        if (spectrum.y.length > 0) {
          var key = 'y';
          thisdata[key] = 0;
          for( var j = 0; j < newRebin && i+j<spectrum.y.length; ++j )
            thisdata[key] += spectrum.y[i+j];
          thisdata[key] *= spectrum.yScaleFactor;
        }

        spectrum.points.push( thisdata );
      }
    }
  });
}

SpectrumChartD3.prototype.update = function() {
  var self = this;

  if (!this.rawData || !this.rawData.spectra || !this.rawData.spectra.length)
    return;

  for (var i = 0; i < this.rawData.spectra.length; ++i) {
    if (this['line'+i])
      this.vis.select("#spectrumline"+i).attr("d", this['line'+i](this.rawData.spectra[i].points));
  }


  /*
  if( self.xrange !== self.xScale.range() ) {
    console.log( "xrange changed changed" );
  }

  if( self.xdomain !== self.xScale.domain() ) {
    console.log( "xrange domain changed" );
  }

  if( self.yrange !== self.yScale.range() ) {
    console.log( "yrange range changed" );
  }

  if( self.ydomain !== self.yScale.domain() ) {
    console.log( "yrange domain changed" );
  }

  if( self.prevRebinFactor !== this.rebinFactor ) {
    console.log( "rebin factor changed" );
  }

  self.prevRebinFactor = this.rebinFactor;
  self.xrange = self.xScale.range();   //x height in px
  self.xdomain = self.xScale.domain(); //y in coordinate values
  self.yrange = self.yScale.range();   //x height in px
  self.ydomain = self.yScale.domain(); //y in coordinate values
  */

  if (d3.event && d3.event.keyCode) {
    d3.event.preventDefault();
    d3.event.stopPropagation();
  }
}

SpectrumChartD3.prototype.dataPointDrag = function() {
  var self = this;
  return function(d) {
    registerKeyboardHandler(self.keydown());
    document.onselectstart = function() { return false; };
    self.selected = self.dragged = d;
    self.update(false); // boolean set to false to indicate no animation needed

  }
}

SpectrumChartD3.prototype.setData = function( data, resetdomain ) {
  //need to make some consistency checks on data here
  //  -Has all necassary variables
  //  -Energy is monotonically increasing
  //  -All y's are the same length, and consistent with x.
  //  -No infs or nans.

  var self = this;

  if (data && data.spectra && data.spectra.length)
    for (var i = 0; i < data.spectra.length; ++i)
      if (this['line'+i])
        this.vis.selectAll("#spectrumline"+i).remove();

  this.vis.selectAll('path.line').remove();
  this.rawData = null;
  this.firstRaw = this.lastRaw = this.rebinFactor = -1; //force rebin calc

  try
  {
    if( !data || !data.spectra ) throw null;
    if( !Array.isArray(data.spectra) || data.spectra.length < 1 ) throw 'No spectrum-data specified';
    //check that x is same length or one longer than y.
    //Check that all specified y's are the same length
  }catch(e){
    if(e) console.log(e);
    this.redraw()();
    return;
  }

  this.rawData = data;

  for (var i = 0; i < this.rawData.spectra.length; ++i)
    this['line'+i] = null;
  
  this.do_rebin();  //this doesnt appear to be necassay, but JIC

  //Make it so the x-axis shows all the data
  if( resetdomain ) {
    var bounds = self.min_max_x_values();
    var minx = bounds[0], maxx = bounds[1];
    
    this.xScale.domain([minx, maxx]);
  }

  //reset the zoom, so first userzoom action wont behave wierdly
  // this.zoom.x(this.xScale);
  // this.zoom.y(this.yScale);

  // Hack: To properly choose the right set of points for the y-axis points
  function y(line) {
    return function(d) {
      var y = self.yScale(d['y']);
      if (isNaN(y)) y = 0;
      return y; 
    }
  }

  var maxYScaleFactor = 0.1;
  for (var i = 0; i < this.rawData.spectra.length; ++i)
    this.rawData.spectra[i].dataSum = 0;

  // Create the lines
  for (var i = 0; i < data.spectra.length; ++i) {
    var spectrumi = i;
    var spectrum = data.spectra[i];
    if (spectrum.y.length) {
      for (var j = 0; j < spectrum.y.length; ++j) {
        spectrum.dataSum += spectrum.y[j];
      }
      this['line' + i] = d3.svg.line()
        .interpolate("step-after")
        .x( function(d, pi) {
          return self.xScale(d.x);
        })
        .y( y(i) );

      this.chartBody.append("path")
        .attr("id", "spectrumline"+i)
        .attr("class", 'line')
        .attr("stroke", spectrum.lineColor ? spectrum.lineColor : 'black')
        .attr("d", this['line' + i](spectrum.points));

      if (spectrum.yScaleFactor)
        maxYScaleFactor = Math.max(spectrum.yScaleFactor, maxYScaleFactor);
    }
  }

  // + 10 to add at keast some padding for scale
  self.options.maxScaleFactor = maxYScaleFactor + 10;
  var maxsfinput;
  if (maxsfinput = document.getElementById("max-sf")) {
    maxsfinput.value = self.options.maxScaleFactor;
  }

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  // Update the spectrum drop down for adjusting y-scale factors
  if (currentsftitle = document.getElementById("current-sf-title")) {
    var titles = graph.getSpectrumTitles();
    currentsftitle.options.length = titles.length;
    titles.forEach(function(title,i) {
      if (i == 0)
        currentsftitle.options[i] = new Option("None", "", true, false);
      else
        currentsftitle.options[i] = new Option(title, title, false, false);
    });
  }

  this.addMouseInfoBox();

  this.updateLegend();
  this.drawScalerBackgroundSecondary();
  
  this.redraw()();
}


SpectrumChartD3.prototype.setReferenceLines = function( data ) {
  d3.selectAll("g.ref").remove();

  var default_colors = ["#0000FF","#006600", "#006666", "#0099FF","#9933FF", "#FF66FF", "#CC3333", "#FF6633","#FFFF99", "#CCFFCC", "#000CC", "#666666", "#003333"];

  if( !data ){
    this.refLines = null;
  } else {
    try {
      if( !Array.isArray(data) )
        throw "Input is not an array of reference lines";

      data.forEach( function(a,i){
         if( !a.color )
           a.color = default_colors[i%default_colors.length];
         if( !a.lines || !Array.isArray(a.lines) )
           throw "Reference lines does not contain an array of lines";
         a.lines.forEach( function(d){
           //{e:30.27,h:6.22e-05,particle:'xray',decay:'xray',el:'barium'}
           //particle in ["gamma", "xray", "beta", "alpha",   "positron", "electronCapture"];
           if( (typeof d.e !== "number") || (typeof d.h !== "number") || (typeof d.particle !== "string") )
             throw "Refernce line is invalid (" + JSON.stringify(d) + ")";
         });
       });

    //this.refLines = JSON.parse(JSON.stringify(data));  //creates deep copy, but then also have to go through and
    this.refLines = data;
    this.refLines.forEach( function(a,i){ a.lines.forEach( function(d){ d.parent = a; } ) } );
  }catch(e){
    this.refLines = null;
    console.log( "invalid input to setReferenceLines" );
  }
    this.refLines = data;
  }

  this.redraw()();
}

SpectrumChartD3.prototype.addMouseInfoBox = function(){
  if( this.mouseInfo )
    this.mouseInfo.remove();

  this.mouseInfo = this.vis.append("g")
                     .attr("class", "mouseInfo")
                     .style("display", "none")
                     .attr("transform","translate(" + this.size.width + "," + this.size.height + ")");

  this.mouseInfoBox = this.mouseInfo.append('rect')
               .attr("class", "mouseInfoBox")
               .attr('width', "12em")
               .attr('height', "2.5em")
               .attr('x', "-12.5em")
               .attr('y', "-3.1em");

  this.mouseInfo.append("g").append("text");
}


SpectrumChartD3.prototype.setLogY = function(){
  //To make the transition animated, see: http://bl.ocks.org/benjchristensen/2657838

  if( this.options.yscale === "log" )
    return;

  this.options.yscale = "log";

  this.yScale = d3.scale.log()
      .clamp(true)
      .domain([0, 100])
      .nice()
      .range([1, this.size.height])
      .nice();

  if( this.yGrid )
    this.yGrid.scale( this.yScale );

  this.redraw(this.options.showAnimation)();
}

SpectrumChartD3.prototype.setLinearY = function(){
  if( this.options.yscale === "lin" )
    return;

  this.options.yscale = "lin";
  this.yScale = d3.scale.linear()
      .domain([0, 100])
      .nice()
      .range([0, this.size.height])
      .nice();

  if( this.yGrid )
    this.yGrid.scale( this.yScale );

  this.redraw(this.options.showAnimation)();
}

SpectrumChartD3.prototype.setSqrtY = function(){
  if( this.options.yscale === "sqrt" )
    return;

  this.options.yscale = "sqrt";
  this.yScale = d3.scale.pow().exponent(0.5)
      .domain([0, 100])
      .range([0, this.size.height]);

  if( this.yGrid )
    this.yGrid.scale( this.yScale );

  this.redraw(this.options.showAnimation)();
}

SpectrumChartD3.prototype.setGridX = function( onstate ) {
  if( this.options.gridx == onstate )
    return;

  this.options.gridx = onstate;

  if( onstate ) {
    this.xGrid = d3.svg.axis().scale(this.xScale)
                   .orient("bottom")
                   .innerTickSize(-this.size.height)
                   .outerTickSize(0)
                   .tickFormat( "" )
                   .tickPadding(10)
                   .ticks( 20,"" );

    this.xGridBody = this.vis.insert("g", ".peakVis")
        .attr("width", this.size.width )
        .attr("height", this.size.height )
        .attr("class", "xgrid" )
        .attr("transform", "translate(0," + this.size.height + ")")
        .call( this.xGrid );
  } else {
    this.xGridBody.remove();
    this.xGrid = null;
    this.xGridBody = null;
  }

  this.redraw(true)();
}

SpectrumChartD3.prototype.setGridY = function( onstate ) {
  if( this.options.gridy == onstate )
    return;

  this.options.gridy = onstate;

  if( onstate ) {
    this.yGrid = d3.svg.axis().scale(this.yScale)
                   .orient("left")
                   .innerTickSize(-this.size.width)
                   .outerTickSize(0)
                   .tickFormat( "" )
                   .tickPadding(10);

    this.yGridBody = this.vis.insert("g", ".peakVis")
        .attr("width", this.size.width )
        .attr("height", this.size.height )
        .attr("class", "ygrid" )
        .attr("transform", "translate(0,0)")
        .call( this.yGrid );
  } else {
    this.yGridBody.remove();
    this.yGrid = null;
    this.yGridBody = null;
  }

  this.redraw()();
}


SpectrumChartD3.prototype.setAdjustYAxisPadding = function( adjust, pad ) {
  this.options.adjustYAxisPadding = Boolean(adjust);
  
  if( typeof pad === "number" )
    this.padding.left = pad;
    
  this.handleResize( false );
}

SpectrumChartD3.prototype.setCompactXAxis = function( compact ) {
  this.options.compactXAxis = Boolean(compact);
  
  //Might want to add ability to change xTitlePad here. 
  //this.padding.xTitlePad = 10;
    
  this.handleResize( false );
}

SpectrumChartD3.prototype.setShowLegend = function( show ) {
  this.options.showLegend = Boolean(show);
  this.updateLegend();
}


/* Mouse/touch interactions for the chart SVG */
SpectrumChartD3.prototype.handleChartMouseMove = function() {
  var self = this;

  return function() {
    self.mousemove()();

    // If no data is detected, then stop updating other mouse move parameters
    if(!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length || !isNaN(self.downx) || !isNaN(self.downy) || self.legdown || self.scalerdown) 
      return;

    // Prevent and stop default events from ocurring
    d3.event.preventDefault();
    d3.event.stopPropagation();

    // Get the absoulate minimum and maximum x-values valid in the data
    var minx, maxx, bounds;
    if( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length ){
      minx = 0;
      maxx = 3000;
    }else {
      bounds = self.min_max_x_values();
      minx = bounds[0];
      maxx = bounds[1];
    }
     
    // Set the mouse and chart parameters
    var m, x0px, x1px;

    m = d3.mouse(self.vis[0][0]);
    x0_min = self.xScale.range()[0];
    x1_max = self.xScale.range()[1];
    self.lastMouseMovePos = m;

    // Adjust the last mouse move position in case user starts draggin from out of bounds
    if (self.lastMouseMovePos[1] < 0)
      self.lastMouseMovePos[1] = 0;
    else if (self.lastMouseMovePos[1] > self.size.height)
      self.lastMouseMovePos[1] = self.size.height;

    // Do any other mousemove events that could be triggered by the left mouse drag
    self.handleMouseMoveSliderChart()();
    self.handleMouseMoveScaleFactorSlider()();


    if (d3.event.buttons === 1) {        // If left click being held down (left-click and drag)
      d3.select(document.body).attr("cursor", "move");

      // Holding the Shift-key and left-click dragging --> Delete Peaks Mode
      self.isDeletingPeaks = d3.event.shiftKey && !d3.event.altKey && !d3.event.ctrlKey && !d3.event.metaKey && !self.fittingPeak && !self.escapeKeyPressed;

      // Holding the Alt+Shift-key and left-click dragging --> Count Gammas Mode
      self.isCountingGammas = d3.event.altKey && d3.event.shiftKey && !d3.event.ctrlKey && !d3.event.metaKey && !self.fittingPeak && !self.escapeKeyPressed;

      // Holding the Alt+Ctrl-key + Left-click Dragging --> Recalibration Mode
      self.isRecalibrating = d3.event.altKey && d3.event.ctrlKey && !d3.event.metaKey && !d3.event.shiftKey && !self.fittingPeak && !self.escapeKeyPressed; 

      // Holding the Command-key + Left-click dragging --> Zoom-in Y Mode
      self.isZoomingInYAxis = d3.event.metaKey && !d3.event.altKey && !d3.event.ctrlKey && !d3.event.shiftKey && !self.fittingPeak && !self.escapeKeyPressed;

      // Holding the Alt-key + Left-click dragging ---> Undefined, maybe a future implementation?
      self.isUndefinedMouseAction = d3.event.altKey && !d3.event.ctrlKey && !d3.event.metaKey && !d3.event.shiftKey && !self.fittingPeak && !self.escapeKeyPressed;

      var isZoomingInXAxis = !d3.event.altKey && !d3.event.ctrlKey && !d3.event.metaKey && !d3.event.shiftKey && !self.fittingPeak && !self.escapeKeyPressed;

      // If erasing peaks
      if( self.fittingPeak ) {          // If fitting peaks
        self.handleMouseMovePeakFit();

      } else if (self.isDeletingPeaks) {   // If deleting peaks
        self.handleMouseMoveDeletePeak();

      } else if (self.isCountingGammas) {   // If counting gammas
        self.handleMouseMoveCountGammas();

      } else if (self.isRecalibrating) {      // If recalibrating the chart
        self.handleMouseMoveRecalibration();

      } else if (self.isZoomingInYAxis) {        // If zooming in y-axis
        self.handleMouseMoveZoomInY();

      } else if (self.isUndefinedMouseAction) {   // If undefined mouse action
        // Do nothing, save for future feature
        self.handleCancelAllMouseEvents()();

      } else if (isZoomingInXAxis) {    // If zooming in x-axis
        self.handleMouseMoveZoomInX();

      } else {
        self.handleCancelAllMouseEvents()();
      }

      return;

    } else if ( d3.event.buttons === 2 )  // Right Click Dragging: pans the chart left and right
      self.handlePanChart();
    

    self.updateFeatureMarkers(-1);
  }
}

SpectrumChartD3.prototype.handleChartMouseOut = function() {
  var self = this;

  return function () {
    if (!d3.event)
        return;

      if (!d3.select(d3.event.toElement)[0].parentNode || d3.event.toElement === document.body || d3.event.toElement.nodeName === "HTML" 
              || d3.event.toElement.nodeName === "DIV" || d3.event.toElement.offsetParent === document.body) {

        /* For debugging where the mouse is specifically out of */
        // if (!d3.select(d3.event.toElement)[0].parentNode)
        //   console.log("mouse out of the window");
        // else
        //   console.log("mouse out of the chart");

        // Cancel erasing peaks
        self.handleCancelMouseDeletePeak();

        // Cancel the right-click-drag action
        self.handleCancelMouseRecalibration();

        // Cancel the left-click-drag zoom action
        self.handleCancelMouseZoomInX();

        // Cancel the left-click-drag zooming in y axis action
        self.handleCancelMouseZoomInY();

        // Cancel count gammas
        self.handleCancelMouseCountGammas();
      }

      self.updateFeatureMarkers(-1);

      self.mousedOverRefLine = null;
      self.refLineInfo.style("display", "none");
      self.mouseInfo.style("display", "none");
  }
}

SpectrumChartD3.prototype.handleChartMouseUp = function() {
  var self = this;

  return function () {
    self.downx = Math.NaN;
    self.downy = Math.NaN;
    d3.select(document.body).style("cursor", "default");

    /* Here we can decide to either zoom-in or not after we let go of the left-mouse button from zooming in on the CHART.
       In this case, we are deciding to zoom-in if the above event occurs.

       If you want to do the opposite action, you call self.handleCancelMouseZoomInX()
     */

    self.handleMouseUpDeletePeak();

    self.handleMouseUpZoomInX();

    self.handleMouseUpZoomInY();

    self.handleMouseUpRecalibration();

    self.handleMouseUpCountGammas();

    self.lastMouseMovePos = null;
    self.sliderChartMouse = null;
  }
}

SpectrumChartD3.prototype.handleChartWheel = function () {
  var self = this;

  return function() {
    // Keep event from bubbling up any further
    if (d3.event) {
      d3.event.preventDefault();
      d3.event.stopPropagation();
    } else
      return;

    //Get mouse pixel x and y location
    var m = d3.mouse(self.vis[0][0]);

    // Handle y axis zooming if wheeling in y-axis
    if (m[0] < 0 && m[1] > 0 && m[0] < self.size.height) {
      self.handleYAxisWheel();
      return;
    }
  }
}

SpectrumChartD3.prototype.handleChartTouchStart = function() {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();
  }
}

SpectrumChartD3.prototype.handleChartTouchEnd = function() {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();

    self.sliderBoxDown = false;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;
    self.sliderChartTouch = null;
    self.savedSliderTouch = null;
  }
}


/* Mouse/Touch interactions for the vis svg */
SpectrumChartD3.prototype.handleVisMouseDown = function () {
  var self = this;
  // console.log("mousedown on plot!");

  return function () {
    self.dragging_plot = true;

    self.updateFeatureMarkers(null);

    self.mousedowntime = new Date();
    self.mousedownpos = d3.mouse(document.body);

    registerKeyboardHandler(self.keydown());

    if( !isNaN(self.downx) || !isNaN(self.downy) || self.legdown || self.scalerdown )
      return;

    // Cancel the default d3 event properties
    d3.event.preventDefault();
    d3.event.stopPropagation(); 

    var m = d3.mouse(self.vis[0][0]);

    if (d3.event.metaKey)
      self.zoomInYMouse = m;
    else
      self.zoomInYMouse = null;

    self.leftMouseDown = null;
    self.zoominmouse = self.deletePeaksMouse = self.countGammasMouse = self.recalibrationMousePos = null; 

    if( d3.event.buttons === 1 && m[0] >= 0 && m[0] < self.size.width && m[1] >= 0 && m[1] < self.size.height ) {    // if left click-and-drag and mouse is in bounds
      console.log("left mousedown");

      // Set cursor to move icon
      d3.select('body').style("cursor", "move");

      // Set the zoom in/erase mouse properties
      self.leftMouseDown = self.zoominmouse = self.deletePeaksMouse = self.countGammasMouse = self.recalibrationMousePos = m;
      self.peakFitMouseDown = m;
      self.origdomain = self.xScale.domain();
      self.zoominaltereddomain = false;
      self.zoominx0 = self.xScale.invert(m[0]);
      self.zoominy0 = self.xScale.invert(m[1]);

      self.recalibrationStartEnergy = [ self.xScale.invert(m[0]), self.xScale.invert(m[1]) ];
      self.isRecalibrating = false;

      // We are fitting peaks (if alt-key held)
      self.fittingPeak = d3.event.ctrlKey && !d3.event.altKey && !d3.event.metaKey && !d3.event.shiftKey && d3.event.keyCode !== 27;

      // Initially set the escape key flag false
      self.escapeKeyPressed = false;
      
      // Create the initial zoom box if we are not fitting peaks
      if( !self.fittingPeak ) {
        var zoomInXBox = d3.select("#zoomInXBox")
            zoomInXText = d3.select("#zoomInXText");

        zoomInXBox.remove();
        zoomInXText.remove();

        // Set the zoom-in box and display onto chart
        zoomInXBox = self.vis.append("rect")
          .attr("id", "zoomInXBox")
          .attr("class","leftbuttonzoombox")
          .attr("width", 1 )
          .attr("height", self.size.height)
          .attr("x", m[0])
          .attr("y", 0)
          .attr("pointer-events", "none");
      }

      self.updateFeatureMarkers(-1);
      
      self.zooming_plot = true;
      return false;

    } else if ( d3.event.buttons === 2 ) {    // listen to right-click mouse down event
      console.log("Right mouse down!");
      self.rightClickDown = d3.mouse(document.body);
      self.rightClickDrag = false;

      // Since this is the right-mouse button, we are not zooming in
      self.zooming_plot = false;
    }
  }
}

SpectrumChartD3.prototype.handleVisMouseUp = function () {
  var self = this;

  return function () {
    console.log("mouseup on vis!");

    if (!d3.event)
      return;

    // Set the client/page coordinates of the mouse
    var m = d3.mouse(self.vis[0][0])
        x = m[0],
        y = m[1],
        pageX = d3.event.pageX,
        pageY = d3.event.pageY,
        mouseEnergy = self.xScale.invert(x);

    // Handle any of default mouseup actions
    self.mouseup()();

    // We are not dragging the plot anymore
    self.dragging_plot = false;

    // Update feature marker positions
    self.updateFeatureMarkers(null);

    // If the slider chart is displayed and the user clicks on that, cancel the mouse click action
    if (y >= (self.size.height + 
      (self.xaxistitle != null && !d3.select(self.xaxistitle).empty() ? self.xaxistitle[0][0].clientHeight + 20 : 20) + 
      self.padding.sliderChart)) {
      return;
    }

    // Figure out right clicks
    if (d3.event.button === 2 && !self.rightClickDrag) {
      var foundPeakFromRightClick = false;
      var peakTitle;

      if (self.rawData.spectra) {

        self.rawData.spectra.forEach(function(spectra,specindex) {
          spectra.peaks.forEach( function(peak) {
            if (peak.lowerEnergy <= mouseEnergy && mouseEnergy <= peak.upperEnergy) {
              peakTitle = spectra.title ? spectra.title : "";
              foundPeakFromRightClick = true;

              console.log("Emit RIGHT CLICK (ON '" + peakTitle + "' PEAK) signal. (Peak = ", peak, ")", 
                          "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY );
              return;
            }
          });
        });
      }
      
      if (!foundPeakFromRightClick)
        console.log("Emit RIGHT CLICK (ON PLOT) signal!\nx = ", 
                    x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY);
      return;
    }
    
    // Figure out clicks and double clicks
    var nowtime = new Date();
    var clickDelay = 200;

    if (self.mousedownpos && self.dist(self.mousedownpos, d3.mouse(document.body)) < 5) {    // user clicked on screen
      if( nowtime - self.mousedowntime < clickDelay ) {

        if (self.lastClickEvent && nowtime - self.lastClickEvent < clickDelay) {    // check for double click
          if (self.mousewait) {
            window.clearTimeout(self.mousewait);
            self.mousewait = null;
          }
          console.log("Emit DOUBLE CLICK signal!", "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY);
        } else {
          self.mousewait = window.setTimeout((function(e) {
            self.updateFeatureMarkers(self.xScale.invert( x ));    // update the sum peak where user clicked

            return function() {
                      console.log( "Emit CLICK signal!", "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY );
                
                      self.unhighlightPeakFunc();
                      self.mousewait = null;
                  }
          })(d3.event), clickDelay);
        }
        self.lastClickEvent = new Date();
      }
    }

    if( !isNaN(self.downx) || !isNaN(self.downy) || self.legdown || self.scalerdown )
      return;

    // Handle fitting peaks (if needed)
    if (self.fittingPeak)
      self.handleMouseUpPeakFit();

    // Handle zooming in x-axis (if needed)
    if (self.zooming_plot) {
      self.handleMouseUpZoomInX();
    }

    // Handle deleting peaks (if needed)
    self.handleMouseUpDeletePeak();
    
    // Handle recalibration (if needed)
    self.handleMouseUpRecalibration();

    // Handle zooming in y-axis (if needed)
    self.handleMouseUpZoomInY();

    // HAndle counting gammas (if needed)
    self.handleMouseUpCountGammas();

    // Delete any other mouse actions going on
    self.leftMouseDown = null;
    self.zoominbox = null;
    self.zoominx0 = null;
    self.zoominy0 = null;
    self.zoominmouse = null;
    self.fittingPeak = null;
    self.escapeKeyPressed = false;

    // Set delete peaks mode off
    self.isDeletingPeaks = false;
    self.deletePeaksMouse = null;

    // Set the count gammas mode off
    self.countGammasMouse = null;

    // Set the recalibration mode off
    self.recalibrationMousePos = null;

    self.handleCancelMouseZoomInY();

    // Set the right click drag off
    self.rightClickDrag = false;

    // Cancel default d3 event properties
    d3.event.preventDefault();
    d3.event.stopPropagation();

    // Not zooming in anymore
    self.zooming_plot = false;

    // Not using x-axis slider anymore
    self.sliderBoxDown = false;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;
    self.sliderChartMouse = null;
    self.savedSliderMouse = null;
  }
}

SpectrumChartD3.prototype.handleVisWheel = function () {
  var self = this;

  return function () {
    var e = d3.event;

    //Keep event from bubbling up any further
    e.preventDefault();
    e.stopPropagation();

    //If the user is doing anything else, return
    //Note that if you do a two finger pinch on a mac book pro, you get e.ctrlKey==true and e.composed==true 
    if( !e || e.altKey || (e.ctrlKey && !e.composed) || e.shiftKey || e.metaKey || e.button != 0 || e.buttons != 0 ) {
     console.log( "Special condition with wheel, ignoring mousewheel" );
     console.log( e );
     return;
    }

    //Get mouse pixel x and y location
    var m = d3.mouse(self.vis[0][0]);

    //Make sure within chart area
    if( m[0] < 0 || m[0] > self.size.width || m[1] < 0 || m[1] > self.size.height ){

    // If wheeling in y-axis labels, zoom in the y-axis range
    if (m[0] < 0 && m[1] > 0 && m[0] < self.size.height) {
      self.handleYAxisWheel();
      return;
    }

     console.log( "Scroll outside of vis, ignoring mousewheel" );
     return;
    }  

    //If we are doing any other actions with the chart, then to bad.
    if( self.dragging_plot || self.zoominbox || self.fittingPeak ){
     console.log( "Plot is being dragged, zoomed, or peak fit, ignoring mousewheel" );
     return;
    }

    //Dont do anything if there is no data
    if( !self.rawData || !self.rawData.spectra || self.rawData.spectra.length < 2 ){
     console.log( "No data, ignoring mousewheel" );
     return;
    }

    var mindatax, maxdatax, bounds, foreground;
    foreground = self.rawData.spectra[0];
    bounds = self.min_max_x_values();
    mindatax = bounds[0];
    maxdatax = bounds[1];

    if ((self.xScale.domain()[0] == mindatax && e.deltaX < 0) || (self.xScale.domain()[1] == maxdatax && e.deltaX > 0)) {
      console.log("Cannot scroll past the minimum/maximum data from chart, ignoring mousewheel")
      return;
    }

    //Here we will set a timer so that if it is more than 1 second since the
    //  last wheel movement, we will consider the wheel event over, and start 
    //  fresh.
    //  This is just an example, and likely needs changed, or just removed.
    if( self.wheeltimer ){
     //This is not the first wheel event of the current user wheel action,
     //  lets clear the previous timeout (we'll reset a little bellow). 
     window.clearTimeout(self.wheeltimer);
    } else {
     //This is the first wheel event of this user wheel action, lets record
     //  initial mouse energy, counts, as well as the initial x-axis range.
     self.scroll_start_x = self.xScale.invert(m[0]);
     self.scroll_start_y = self.yScale.invert(m[1]);
     self.scroll_start_domain = self.xScale.domain();
     self.scroll_start_raw_channel = d3.bisector(function(d){return d;}).left(foreground.x, self.scroll_start_x);
     self.scroll_start_raw_channel = Math.max(0,self.scroll_start_raw_channel);
     self.scroll_start_raw_channel = Math.min(foreground.x.length-1,self.scroll_start_raw_channel);
     self.scroll_total_x = 0;
     self.scroll_total_y = 0;
    }

    //Function to clear out any variables assigned during scrolling, or finish
    //  up any actions that should be done 
    function wheelcleanup(e){
     console.log( "mousewheel, stopped" );
     self.wheeltimer = null;
     self.scroll_start_x = null;
     self.scroll_start_y = null;
     self.scroll_start_domain = null;
     self.scroll_start_raw_channel = null;
     self.scroll_total_x = null;
     self.scroll_total_y = null;
    }

    //Set a timeout to call wheelcleanup after a little time of not resieving 
    //  any user wheel actions.
    self.wheeltimer = window.setTimeout( wheelcleanup, 1000 );

    //scroll_total_x/y is the total number of scroll units the mouse has done
    //  since the user started doing this current mouse wheel.
    self.scroll_total_x += e.deltaX;
    self.scroll_total_y += e.deltaY;

    var MAX_SCROLL_TOTAL = 200;

    //Zoom all the way out by the time we get self.scroll_total_y = +200,
    //  or, zoom in to 3 bins by the time self.scroll_total_y = -200;
    self.scroll_total_y = Math.max( self.scroll_total_y, -MAX_SCROLL_TOTAL );
    self.scroll_total_y = Math.min( self.scroll_total_y, MAX_SCROLL_TOTAL );

    console.log("wheel on chart {" + self.scroll_total_x + "," + self.scroll_total_y + "}");

    var initial_range_x = self.scroll_start_domain[1] - self.scroll_start_domain[0];
    var terminal_range_x;
    if( self.scroll_total_y > 0 ){
      terminal_range_x = (maxdatax - mindatax);
    } else {
     //Find the bin one to the left of original mouse, and two to the right 
     var terminalmin = Math.max(0,self.scroll_start_raw_channel - 1);
     var terminalmax = Math.min(foreground.x.length-1,self.scroll_start_raw_channel + 2);
     terminalmin = foreground.x[terminalmin];
     terminalmax = foreground.x[terminalmax];
     terminal_range_x = terminalmax - terminalmin; 
    }

    var frac_y = Math.abs(self.scroll_total_y) / MAX_SCROLL_TOTAL;
    var new_range = initial_range_x + frac_y * (terminal_range_x - initial_range_x);  

    //Make it so the mouse is over the same energy as when the wheeling started
    var vis_mouse_frac = m[0] / self.size.width;  

    var new_x_min = self.scroll_start_x - (vis_mouse_frac * new_range);

    var new_x_max = new_x_min + new_range;


    //Now translate the chart left and right.  100 x wheel units is one initial 
    //  width left or right
    //  TODO: should probably make it so that on trackpads at least, there is
    //        is some threshold on the x-wheel before (like it has to be 
    //        greater than the y-wheel) before the panning is applied; this 
    //        would also imply treating the x-wheel using deltas on each event
    //        instead of using the cumulative totals like now.
    var mouse_dx_wheel = Math.min(initial_range_x,new_range) * (self.scroll_total_x / MAX_SCROLL_TOTAL);
    new_x_min += mouse_dx_wheel;
    new_x_max += mouse_dx_wheel;

    if( new_x_min < mindatax ){
     new_x_max += (mindatax - new_x_min);
     new_x_min = mindatax; 
    }

    if( new_x_max > maxdatax ){
     new_x_min = Math.max(mindatax,new_x_min-(new_x_max-maxdatax));
     new_x_max = maxdatax;
    }

    //Finally set the new x domain, and redraw the chart (which will take care
    //  of setting the y-domain).
    self.xScale.domain( [new_x_min,new_x_max] );
    self.redraw()();

    self.updateFeatureMarkers(-1);
  }
}

SpectrumChartD3.prototype.handleVisTouchStart = function() {
  var self = this;

  return function() {

      // Prevent default event actions from occurring (eg. zooming into page when trying to zoom into graph)
      d3.event.preventDefault();
      d3.event.stopPropagation();

      // Get the touches on the screen
      var t = d3.touches(self.vis[0][0]),
          touchHoldTimeInterval = 600;

      // Save the original zoom scale and translation
      self.savedZoomScale = self.zoom.scale();
      self.savedZoomTranslation = self.zoom.translate();

      // Represent where we initialized our touch start value
      self.touchStart = t;
      self.touchStartEvent = d3.event;
      self.touchPageStart = d3.touches(document.body).length === 1 ? [d3.event.pageX, d3.event.pageY] : null;

      if (t.length === 2) {
        self.countGammasStartTouches = self.createPeaksStartTouches = self.touchStart;
      }

      self.updateTouchesOnChart(self.touchStartEvent);

      // Boolean for the touch of a touch-hold signal
      self.touchHoldEmitted = false;

      self.touchHold = window.setTimeout((function(e) {
        var x = t[0][0],
            y = t[0][1],
            pageX = d3.event.pageX,
            pageY = d3.event.pageY;

        return function() {

          // Emit the tap signal, unhighlight any peaks that are highlighted
          if (self.touchStart && self.dist([pageX, pageY], self.touchPageStart) < 5 && !self.touchHoldEmitted) {
            console.log( "Emit TAP HOLD (RIGHT TAP) signal!", "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY );
            self.unhighlightPeakFunc();
            self.touchHoldEmitted = true;
          }

          // Clear the touch hold wait, the signal has already been emitted
          self.touchHold = null;
        }     
      })(d3.event), touchHoldTimeInterval);
    }
}

SpectrumChartD3.prototype.handleVisTouchMove = function() {
  var self = this;

  // Touch interaction helpers
  function isDeletePeakSwipe() {

    if (!self.touchesOnChart)
      return false;

    var keys = Object.keys(self.touchesOnChart);

    // Delete peak swipe = two-finger vertical swipe
    if (keys.length !== 2)
      return false;

    var maxDyDiff = 15,
        minDxDiff = 25;

    var t1 = self.touchesOnChart[keys[0]],
        t2 = self.touchesOnChart[keys[1]];

    var dy1 = t1.startY - t1.pageY,
        dy2 = t2.startY - t2.pageY,
        dyDiff = Math.abs(dy2 - dy1);

    if (dyDiff > maxDyDiff || Math.abs(t1.pageX - t2.pageX) < minDxDiff) {
      return false;
    }

    var dy = Math.min(dy1,dy2),
        dx1 = t1.startX - t1.pageX,
        dx2 = t2.startX - t2.pageX,
        dx = Math.abs(dx1 - dx2);

    return dy > dx && dy > maxDyDiff;
  }

  function isControlDragSwipe() {

    if (!self.touchesOnChart)
      return false;

    var keys = Object.keys(self.touchesOnChart);

    if (keys.length !== 2)
      return false;

    var t1 = self.touchesOnChart[keys[0]],
        t2 = self.touchesOnChart[keys[1]];

    if (t1.startX > t1.pageX || t2.startX > t2.pageX)
      return false;

    if( !isFinite(t1.startX) || !isFinite(t1.pageX)
    || !isFinite(t2.startX) || !isFinite(t2.pageX)
    || !isFinite(t1.startY) || !isFinite(t1.pageY)
    || !isFinite(t2.startY) || !isFinite(t2.pageY) )
      return false;

    var startdx = t1.startX - t2.startX;
        nowdx = t1.pageX - t2.pageX;
        yavrg = 0.5*(t1.startY+t2.startY);

    if( Math.abs(yavrg-t1.pageY) > 20 || 
        Math.abs(yavrg-t2.pageY) > 20 || 
        Math.abs(startdx-nowdx) > 20 ) 
      return false;

    return Math.abs(t1.pageX - t1.startX) > 30;
  }

  function isAltShiftSwipe() {
    var keys = Object.keys(self.touchesOnChart);

    if( keys.length !== 2 ) 
      return false;

    var t1 = self.touchesOnChart[keys[0]],
        t2 = self.touchesOnChart[keys[1]];

    if( Math.abs(t1.startX-t2.startX) > 20 || Math.abs(t1.pageX-t2.pageX) > 25 )
      return false;

    return ( (t1.pageX - t1.startX) > 30 );
  }

  function isZoomInYPinch() {
    if (!self.touchesOnChart)
      return false;

    var keys = Object.keys(self.touchesOnChart);

    if (keys.length !== 2)
      return false;

    var touch1 = self.touchesOnChart[keys[0]];
    var touch2 = self.touchesOnChart[keys[1]];
    var adx1 = Math.abs( touch1.startX - touch2.startX );
    var adx2 = Math.abs( touch1.pageX  - touch2.pageX );
    var ady1 = Math.abs( touch1.startY - touch2.startY );
    var ady2 = Math.abs( touch1.pageY  - touch2.pageY );
    var ddx = Math.abs( adx2 - adx1 );
    var ddy = Math.abs( ady2 - ady1 );
    var areVertical = (adx2 > ady2);

    return ddx < ddy && ddy>20
  }

  function deleteTouchLine() {
    // Delete the touch lines if they exist on the vis

    if (self.touchLineX) {
      self.touchLineX.remove();
      self.touchLineX = null;
    }

    if (self.touchLineY) {
      self.touchLineY.remove();
      self.touchLineY = null;
    }
  }

  return function() {

      // Prevent default event actions from occurring (eg. zooming into page when trying to zoom into graph)
      d3.event.preventDefault();
      d3.event.stopPropagation();

      // Nullify our touchstart position, we are now moving our touches
      self.touchStart = null;

      // Get the touches on the chart
      var t = d3.touches(self.vis[0][0]);

      if (t.length === 2) {
        self.deletePeaksTouches = t;
      }

      // Panning = one finger drag
      self.touchPan = t.length === 1;
      self.deletePeakSwipe = isDeletePeakSwipe() && !self.currentlyAdjustingSpectrumScale;
      self.controlDragSwipe = isControlDragSwipe() && !self.currentlyAdjustingSpectrumScale;
      self.altShiftSwipe = isAltShiftSwipe() && !self.currentlyAdjustingSpectrumScale;
      self.zoomInYPinch = isZoomInYPinch() && !self.currentlyAdjustingSpectrumScale;

      if (self.deletePeakSwipe) {
        self.handleTouchMoveDeletePeak();

      } else if (self.controlDragSwipe) {
        self.handleTouchMovePeakFit();

      } else if (self.altShiftSwipe) {
        self.handleTouchMoveCountGammas();

      } else if (self.zoomInYPinch) {
        self.handleTouchMoveZoomInY();

      } else if (self.currentlyAdjustingSpectrumScale) {
        self.handleTouchMoveAdjustSpectrumScale()();

      } else {
        self.handleCancelTouchCountGammas();
        self.handleCancelTouchDeletePeak();
        self.handleCancelTouchPeakFit();
        self.handleTouchCancelZoomInY();
      }

      /* 
      Clear the touch hold signal if:
        - More than one touches on chart
        - No touch positions on the page detected
        - We moved our touches by > 5 pixels
      */
      if (t.length > 1 || !self.touchPageStart || self.dist([d3.event.pageX, d3.event.pageY], self.touchPageStart) > 5) {
        if (self.touchHold) {
          window.clearTimeout(self.touchHold);
          self.touchHold = null;
        }
      }


      // Update mouse coordinates, feature markers on a touch pan action
      if (self.touchPan) {
        self.mousemove();
        self.updateMouseCoordText();
        self.updateFeatureMarkers(-1);
      }

      // Delete the touch line
      deleteTouchLine();

      // Update our map of touches on the chart
      self.updateTouchesOnChart(d3.event);

      self.lastTouches = t;
    }
}

SpectrumChartD3.prototype.handleVisTouchEnd = function() {
  var self = this;

  function updateTouchLine(touches) {

    // Touches is an error, abort function
    if (!touches)
      return;

    // Do not update the touch line (remove it) if there are more than two touches on the screen
    if (touches.length != 1) {
      deleteTouchLine();
      return;
    }

    // Set the coordinates of the touch
    var t = touches[0];

    // Create the x-value touch line, or update its coordinates if it already exists
    if (!self.touchLineX) {
      self.touchLineX = self.vis.append("line")
        .attr("class", "touchLine")
        .attr("x1", t[0])
        .attr("x2", t[0])
        .attr("y1", 0)
        .attr("y2", self.size.height);

    } else {
      self.touchLineX.attr("x1", t[0])
        .attr("x2", t[0])
    }

    // Create the y-value touch line, or update its coordinates if it already exists
    if (!self.touchLineY) {
      self.touchLineY = self.vis.append("line")
        .attr("class", "touchLine")
        .attr("x1", t[0]-10)
        .attr("x2", t[0]+10)
        .attr("y1", t[1])
        .attr("y2", t[1]);
    } else {
      self.touchLineY.attr("x1", t[0]-10)
        .attr("x2", t[0]+10)
        .attr("y1", t[1])
        .attr("y2", t[1]);
    }
  }
  function deleteTouchLine() {

    // Delete the touch lines if they exist on the vis
    if (self.touchLineX) {
      self.touchLineX.remove();
      self.touchLineX = null;
    }

    if (self.touchLineY) {
      self.touchLineY.remove();
      self.touchLineY = null;
    }
  }

  return function() {

      // Prevent default event actions from occurring (eg. zooming into page when trying to zoom into graph)
      d3.event.preventDefault();
      d3.event.stopPropagation();

      // Get the touches on the screen
      var t = d3.event.changedTouches;
      var visTouches = d3.touches(self.vis[0][0]);
      if (visTouches.length === 0) {
        self.touchesOnChart = null;
      }

      if (self.touchPan)
        console.log("touchend from pan!");

      else {
        // Detect tap/double tap signals
        if (t.length === 1 && self.touchStart) {

          // Get page, chart coordinates of event
          var x = self.touchStart[0][0],
              y = self.touchStart[0][1],
              pageX = d3.event.pageX,
              pageY = d3.event.pageY
              currentTapEvent = d3.event;

          // Set the double tap setting parameters
          var tapRadius = 35,                   // Radius area for where a double-tap is valid (anything outside this considered a single tap)
              doubleTapTimeInterval = 200;      // Time interval for double tap

          // Update the feature marker positions (argument added for sum peaks)
          self.updateFeatureMarkers(self.xScale.invert(x));

          // Update the touch line position
          updateTouchLine(self.touchStart);

          // Emit the proper TAP/DOUBLE-TAP signal
          if (self.touchPageStart && self.dist(self.touchPageStart, [pageX, pageY]) < tapRadius ) {

            if (self.lastTapEvent &&
                  self.lastTapEvent.timeStamp && currentTapEvent.timeStamp - self.lastTapEvent.timeStamp < doubleTapTimeInterval &&
                  self.dist([self.lastTapEvent.pageX, self.lastTapEvent.pageY], [pageX, pageY]) < tapRadius) {

              // Clear the single-tap wait if it exists, then emit the double tap signal
              if (self.tapWait) {
                window.clearTimeout(self.tapWait);
                self.tapWait = null;
              }

              // Emit the double-tap signal, clear any touch lines/highlighted peaks in chart
              console.log("Emit DOUBLE TAP signal!", "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY);
              deleteTouchLine();
              self.unhighlightPeakFunc();

            } else {

              // Create the single-tap wait emit action in case there is no more taps within double tap time interval
              self.tapWait = window.setTimeout((function(e) {

                // Move the feature markers to tapped coordinate
                self.updateFeatureMarkers(self.xScale.invert(x));    // update the sum peak where user clicked

                // Don't emit the tap signal if there was a tap-hold
                if (self.touchHoldEmitted)
                  return false;

                return function() {

                  // Emit the tap signal, unhighlight any peaks that are highlighted
                  console.log( "Emit TAP signal!", "\nx = ", x, ", y = ", y, ", pageX = ", pageX, ", pageY = ", pageY );
                  self.unhighlightPeakFunc();
                  var touchStartPosition = self.touchStart;

                  // Highlight peaks where tap position falls
                  if (self.pathEnergies) {
                    for (i = self.pathEnergies.length-1; i >= 0; i--) {
                      if (self.xScale.invert(x) >= self.pathEnergies[i].lower_energy && self.xScale.invert(x) <= self.pathEnergies[i].upper_energy) {
                        self.highlightPeakFunc(0, self.pathEnergies[i].path);
                        break;
                      }
                    }
                  }
                  // Clear the single-tap wait, the signal has already been emitted
                  self.tapWait = null;
                }     
              })(d3.event), doubleTapTimeInterval);
            
            }

            // Set last tap event to current one
            self.lastTapEvent = currentTapEvent;
          }

        } else    // Touch move detected, aborting tap signal
          self.updateFeatureMarkers(-1);
      }

      self.updateTouchesOnChart(d3.event);
      self.mousemove();
      self.updateMouseCoordText();

      self.handleTouchEndCountGammas();
      self.handleTouchEndDeletePeak();
      self.handleTouchEndPeakFit();
      self.handleTouchEndZoomInY();
      

      self.touchPan = false;
      self.touchZoom = false;
      self.touchStart = null;
      self.touchStart = null;
      self.touchHoldEmitted = false;

      self.deletePeakSwipe = false;
      self.controlDragSwipe = false;
      self.altShiftSwipe = false;
      self.zoomInYPinch = false;

      self.countGammasStartTouches = null;

      self.sliderBoxDown = false;
      self.leftDragRegionDown = false;
      self.rightDragRegionDown = false;
      self.sliderChartTouch = null;
      self.savedSliderTouch = null;
    };
}

// Pan chart is called when right-click dragging
SpectrumChartD3.prototype.handlePanChart = function () {
  var self = this;

  var minx, maxx, bounds;
  if( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length ){
    minx = 0;
    maxx = 3000;
  }else {
    bounds = self.min_max_x_values();
    minx = bounds[0];
    maxx = bounds[1];
  }

  // We are now right click dragging 
  self.rightClickDrag = true;

  // For some reason, using the mouse position for self.vis makes panning slightly buggy, so I'm using the document body coordinates instead
  var docMouse = d3.mouse(document.body);

  if (!docMouse || !self.rightClickDown)
    return;

  // Set the pan right / pan left booleans
  var panRight = docMouse[0] < self.rightClickDown[0],   // pan right if current mouse position to the left from previous mouse position
      panLeft = docMouse[0] > self.rightClickDown[0],    // pan left if current mouse position to the right from previous mouse position
      dx = 0;

  function newXDomain() {
    var currentX = docMouse[0],
        delta = currentX - self.rightClickDown[0],
        oldMin = self.xScale.range()[0],
        oldMax = self.xScale.range()[1];

    // Declare new x domain members
    var newXMin = oldMin, 
        newXMax = oldMax;

    newXMin = oldMin - delta;
    newXMax = oldMax - delta;

    newXMin = self.xScale.invert(newXMin);
    newXMax = self.xScale.invert(newXMax);

    // Make sure we don't go set the domain out of bounds from the data
    if( newXMin < minx ){
     newXMax += (minx - newXMin);
     newXMin = minx; 
    }
    if( newXMax > maxx ){
     newXMin = Math.max(minx,newXMin-(newXMax-maxx));
     newXMax = maxx;
    }
    return [newXMin, newXMax];
  }

  // Get the new x values
  var newXMin = newXDomain()[0],
      newXMax = newXDomain()[1];

  // Pan the chart
    self.xScale.domain( [newXMin,newXMax] );
    self.redraw()();

  // New mouse position set to current moue position
  self.rightClickDown = docMouse;
}


/* General mouse/key/touch interaction behaviors */
SpectrumChartD3.prototype.mousemove = function () {
  var self = this;

  return function() {
    //This function is called whenever a mouse movement occurs
    // console.log("mousemove function called from ", d3.event.type);

    var p = d3.mouse(self.vis[0][0]),
        t = d3.event.changedTouches;

    var energy = self.xScale.invert(p[0]),
        mousey = self.yScale.invert(p[1]);


    self.updateFeatureMarkers(-1);
    self.updateMouseCoordText();

    if( self.legdown ) {
      d3.event.preventDefault();
      d3.event.stopPropagation();

      var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX,
          y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY,
          calculated_x = d3.mouse(self.vis[0][0])[0]; // current mouse x position

      if ( calculated_x >= -self.padding.leftComputed && y >= 0 && 
           calculated_x <= self.cx && y <= self.cy ) {
        // console.log("change legend pos");
        var tx = (x - self.legdown.x) + self.legdown.x0;
        var ty = (y - self.legdown.y) + self.legdown.y0; 
        self.legend.attr("transform", "translate(" + tx + "," + ty + ")");
      }
    }

    if (self.scalerdown) {
      d3.event.preventDefault();
      d3.event.stopPropagation();

      var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX,
          y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY,
          calculated_x = d3.mouse(self.vis[0][0])[0]; // current mouse x position

      if ( calculated_x >= -self.padding.leftComputed && y >= 0 && 
           calculated_x <= self.cx && y <= self.cy ) {
        // console.log("change legend pos");
        var tx = (x - self.scalerdown.x) + self.scalerdown.x0;
        var ty = (y - self.scalerdown.y) + self.scalerdown.y0; 
        self.scalerWidget.attr("transform", "translate(" + tx + "," + ty + ")");
      }
    }

    if (self.adjustingBackgroundScale || self.adjustingSecondaryScale) {
      d3.event.preventDefault();
      d3.event.stopPropagation();
    }

    if( self.rawData && self.rawData.spectra && self.rawData.spectra.length ) {
      var foreground = self.rawData.spectra[0];
      var lowerchanval = d3.bisector(function(d){return d.x;}).left(foreground.points,energy,1) - 1;
      var counts = foreground.points[lowerchanval].y;
      var lefteenergy = foreground.points[lowerchanval].x;
      var rightenergy = foreground.points[lowerchanval+1>=foreground.points.length ? foreground.points.length-1 : lowerchanval+1].x;
      var midenergy = 0.5*(lefteenergy + rightenergy);
      var detchannel = lowerchanval*foreground.rebinFactor;

      //Here is where we would update some sort of display or whatever
      //console.log( "counts=" + counts + " in " + self.rebinFactor + " channels, starting at energy " + lefteenergy );

      // self.focus.attr("transform", "translate(" + self.xScale(midenergy) + "," + self.yScale(counts) + ")");
      // self.focus.style("display", null);
      // self.focus.select("text").text( "counts=" + counts + " in " + self.rebinFactor + " channels, starting at energy " + lefteenergy );

    } else {
      // self.focus.style("display", "none");
    }


    if (self.dragged) {
      //We make it here if a data point is dragged (which is not allowed)

      //self.dragged.y = self.yScale.invert(Math.max(0, Math.min(self.size.height, p[1])));
      self.update(false); // boolean set to false to indicate no animation needed
    };
    if (!isNaN(self.downx) && self.xScale.invert(p[0]) > 0) {       // make sure that xaxisDrag does not go lower than 0 (buggy behavior)

      // We make it here when a x-axis label is clicked on, and has been dragged a bit
      d3.select('body').style("cursor", "ew-resize");
      var rupx = self.xScale.invert(p[0]),
          xaxis1 = self.xScale.domain()[0],
          xaxis2 = self.xScale.domain()[1],
          xextent = xaxis2 - xaxis1;
      if (rupx != 0) {
        var changex, new_domain;
        changex = self.downx / rupx;
        new_domain = [xaxis1, xaxis1 + (xextent * changex)];

        if( self.rawData && self.rawData.spectra && self.rawData.spectra.length && new_domain[1] > self.rawData.spectra[0].x[self.rawData.spectra[0].x.length-1] )
          new_domain[1] = self.rawData.spectra[0].x[self.rawData.spectra[0].x.length-1];

        self.xScale.domain(new_domain);
        self.redraw()();
      }

      d3.event.preventDefault();
      d3.event.stopPropagation();
    };

    if (!isNaN(self.downy)) {
      d3.select('body').style("cursor", "ns-resize");
      var rupy = self.yScale.invert(p[1]),
          yaxis1 = self.yScale.domain()[1],
          yaxis2 = self.yScale.domain()[0],
          yextent = yaxis2 - yaxis1;
      if (rupy > 0) {
        var changey, new_domain;
        changey = self.downy / rupy;
        new_domain = [yaxis1 + (yextent * changey), yaxis1];

        self.yScale.domain(new_domain);
        self.redraw()();
      }

      d3.event.preventDefault();
      d3.event.stopPropagation();
    }
  }
}

SpectrumChartD3.prototype.mouseup = function () {
  var self = this;
  return function() {
    document.onselectstart = function() { return true; };
    d3.select('body').style("cursor", "auto");
    d3.select('body').style("cursor", "auto");
    if (!isNaN(self.downx)) {
      self.redraw()();
      self.downx = Math.NaN;
      //d3.event.preventDefault();
      //d3.event.stopPropagation();
    };
    if (!isNaN(self.downy)) {
      self.redraw()();
      self.downy = Math.NaN;
      //d3.event.preventDefault();
      //d3.event.stopPropagation();
    }
    if (self.dragged) {
      self.dragged = null
    }

    self.sliderBoxDown = false;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;
    self.sliderChartMouse = null;
    self.savedSliderMouse = null;
    self.currentlyAdjustingSpectrumScale = null;
  }
}

SpectrumChartD3.prototype.keydown = function () {
  var self = this;
  return function() {

    //if (!self.selected) return;

    switch (d3.event.keyCode) {
      case 27: { //escape
        self.escapeKeyPressed = true;

        if( self.fittingPeak ) {
          self.handleCancelMousePeakFit();
        }

        self.fittingPeak = false;
        self.handleCancelAllMouseEvents()();
        self.handleCancelAnimationZoom();
      }

      case 8: // backspace
      case 46: { // delete
        break;
      }
    }
  }
}

SpectrumChartD3.prototype.updateTouchesOnChart = function (touchEvent) {
  var self = this;

  // Don't do anything if touch event not recognized
  if (!touchEvent || !touchEvent.type.startsWith("touch"))
    return false;

  // Create dictionary of touches on the chart
  if (!self.touchesOnChart)
    self.touchesOnChart = {};

  // Add each touch start into map of touches on the chart
  var touch;
  for (i = 0; i < touchEvent.touches.length; i++) {
    // Get the touch 
    touch = touchEvent.touches[i];

    // Here we add a new attribute to each touch: the start coordinates of the touch
    //   If the touch already exists in our map, then just update the start coordinates of the new touch
    touch.startX = self.touchesOnChart[touch.identifier] ? self.touchesOnChart[touch.identifier].startX : touch.pageX;
    touch.startY = self.touchesOnChart[touch.identifier] ? self.touchesOnChart[touch.identifier].startY : touch.pageY;

    // Add/replace touch into dictionary
    self.touchesOnChart[touch.identifier] = touch;
  }

  // Delete any touches that are not on the screen anymore (read from 'touchend' event)
  if (touchEvent.type === "touchend") {
    for (i = 0; i < touchEvent.changedTouches.length; i++) {
      touch = touchEvent.changedTouches[i];

      if (self.touchesOnChart[touch.identifier])
        delete self.touchesOnChart[touch.identifier];
    }
  }
}


/* Handle specific interaction events */
SpectrumChartD3.prototype.handleZoom = function() {
  var self = this;

  function shouldZoomInX() {
    var keys = Object.keys(self.touchesOnChart);

    if (keys.length == 1) {
      var touch = self.touchesOnChart[keys[0]];
      var x0 = Math.round( touch.startX );
      var x1 = Math.round( touch.pageX );
      return Math.abs(x0-x1)>=15;
    }
    else if (keys.length == 0 || keys.length > 2)
      return false;

    var touch1 = self.touchesOnChart[keys[0]];
    var touch2 = self.touchesOnChart[keys[1]];
    var adx1 = Math.abs( touch1.startX - touch2.startX );
    var adx2 = Math.abs( touch1.pageX  - touch2.pageX );
    var ady1 = Math.abs( touch1.startY - touch2.startY );
    var ady2 = Math.abs( touch1.pageY  - touch2.pageY );
    var ddx = Math.abs( adx2 - adx1 );
    var ddy = Math.abs( ady2 - ady1 );
    var areVertical = (adx2 > ady2);

    return (ddx >= ddy && ddx>20);
  }

  return function() {
    // console.log( 'handleZoom' );
    var e = d3.event;
    var t = d3.touches(self.vis[0][0]);

    // Two touches are close together if distance between two touches < 70px
    //    We want to cancel the zoom if this is true to prevent redrawing while detecting swipes (may need to change/remove this)
    var areTwoTouchesCloseTogether = t.length === 2 && self.dist(t[0], t[1]) < 70; 

    if (e && e.sourceEvent && e.sourceEvent.type.startsWith("mouse")) {
      self.redraw()();
      return;
    }

    // Awesome hack: saving the zoom scale and translation allows us to prevent double-tap zooming!
    if (d3.event.sourceEvent === null || d3.event.sourceEvent.touches == 1 || areTwoTouchesCloseTogether ||
          self.deletePeakSwipe || self.controlDragSwipe || self.altShiftSwipe || self.zoomInYPinch || !shouldZoomInX()) {
      self.xScale.domain(self.xScale.domain());
      self.zoom.scale(self.savedZoomScale);
      self.zoom.translate(self.savedZoomTranslation);
      return false;
    }

    // Get chart domain values
    var xaxismin,
        xaxismax,
        bounds,
        x0 = self.xScale.domain()[0],
        x1 = self.xScale.domain()[1];

    // Set the proper min/max chart values w/respect to data
    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length) {
      xaxismin = 0;
      xaxismax = 3000;
    } else {
      bounds = self.min_max_x_values();
      xaxismin = bounds[0];
      xaxismax = bounds[1];
    }

    // Adjust domains properly to min/max chart values
    if( x0 < xaxismin ){
      x1 += xaxismin - x0;
    }else if( x1 > xaxismax ){
      x0 -= x1 - xaxismax;
    }
      
    // Set the x domain values
    x0 = Math.max( x0, xaxismin );
    x1 = Math.min( x1, xaxismax);

    // Update the x-domain to the new zoomed in values
    self.xScale.domain([x0,x1]);

    // Don't redraw (abort the function) if the domain did not change from last time
    //    Needed to add toPrecision(4) because domain would change (very) slightly sometimes 
    if (self.previousX0 && self.previousX1 && 
        (x0.toPrecision(4) === self.previousX0.toPrecision(4) && x1.toPrecision(4) === self.previousX1.toPrecision(4))) {
      return false;
    }

    // Set previous domain values to current ones
    self.previousX0 = x0;
    self.previousX1 = x1;

    // Redraw the chart
    self.redraw()();
  }
}

SpectrumChartD3.prototype.handleZoomEnd = function () {
  var self = this;

  return function() {

      // Reset the d3 zoom vector so next zoom action will be independant of this one
      self.zoom.x(self.xScale);
      self.zoom.y(self.yScale);
  };
}
  
SpectrumChartD3.prototype.handleYAxisWheel = function() {
  //This function doesnt have the best behaviour in the world, but its a start
  var self = this;
  // console.log("handleYAxisWheel");
  
  if( !d3.event )
    return;
    
  var m = d3.mouse(this.vis[0][0]);
  var t = d3.touches(this.vis[0][0]);

  if( m[0] > 0 && !t )
    return false;
  
  var wdelta;
  // wdelta = d3.event.wheelDelta == null ? d3.event.sourceEvent.wheelDelta : d3.event.wheelDelta;  // old way using source events
  wdelta = d3.event.deltaY ? d3.event.deltaY : d3.event.sourceEvent ? d3.event.sourceEvent.wheelDelta : 0;

  /* Implementation for touch interface */
  // if (!wdelta && wdelta !== 0) {
  //   var t1,t2,dx,dy;

  //   t1 = t[0];
  //   t2 = t[1];

  //   dy = Math.abs(t1[1] - t2[1]);

  //   wdelta = self.previous_dy - dy;
  // }

  
  var mult = 0;
  if( wdelta > 0 ){
    //zoom out
    mult = 0.02;
  } else if( wdelta != 0 ) {
    //zoom in
    mult = -0.02;
  }

  if( self.options.yscale == "log" ) {
    
    if( mult > 0 && self.options.logYFracBottom > 0.025 ){
      self.options.logYFracBottom -= mult;
    } else {
      self.options.logYFracTop += mult;
      self.options.logYFracTop = Math.min( self.options.logYFracTop, 10 );  
    }
    
    if( self.options.logYFracTop < 0 ){
      self.options.logYFracBottom += -self.options.logYFracTop;
      self.options.logYFracBottom = Math.min( self.options.logYFracBottom, 2.505 );
      console.log( 'self.options.logYFracBottom=' +  self.options.logYFracBottom );
      self.options.logYFracTop = 0;
    }
    
    
    self.options.logYFracTop = Math.max( self.options.logYFracTop, -0.95 );
    console.log( 'self.options.logYFracTop=' + self.options.logYFracTop );
    //
  } else if( self.options.yscale == "lin" ) {
    self.options.linYFracTop += mult;
    self.options.linYFracTop = Math.min( self.options.linYFracTop, 0.85 );
    self.options.linYFracTop = Math.max( self.options.linYFracTop, -0.95 );
    //self.options.linYFracBottom = 0.1;
  } else if( self.options.yscale == "sqrt" ) {
    self.options.sqrtYFracTop += mult;
    self.options.sqrtYFracTop = Math.min( self.options.sqrtYFracTop, 0.85 );
    self.options.sqrtYFracTop = Math.max( self.options.sqrtYFracTop, -0.95 );
    //self.options.sqrtYFracBottom = 0.1;  
  }            
            
  self.redraw()();
  
  if( d3.event.sourceEvent ){       
   d3.event.sourceEvent.preventDefault();
   d3.event.sourceEvent.stopPropagation();
  } else {
    d3.event.preventDefault();
    d3.event.stopPropagation();
  }
  return false;
}

SpectrumChartD3.prototype.handleResize = function( dontRedraw ) {
  var self = this;
  
  this.cx = this.chart.clientWidth;
  this.cy = this.chart.clientHeight;

  if (self.sliderChartPlot) {
    this.cy -= self.size.sliderChartHeight + self.padding.sliderChart + 
      (self.xaxistitle != null && !d3.select(self.xaxistitle).empty() ? self.xaxistitle[0][0].clientHeight + 20 : 20); 
  }
  
  var titleh = 0, xtitleh = 0, xlabelh = 7 + 22;
  if( this.options.title ) {
    d3.select('svg').selectAll(".title").forEach( function(t){
      titleh = t[0].getBBox().height;  
   });
  }
  
  if( this.xaxistitle )
    xtitleh = this.xaxistitle.node().getBBox().height;
  
  //Bellow isnt quite working, so stick with estimate of 22 above
  //self.xAxisBody.selectAll('g.tick').forEach( function(t){
      //console.log( t[0].getBBox().height );  //sometimes gives 22.21875, sometimes 4
   //});
  
  this.padding.topComputed = titleh + this.padding.top + (titleh > 0 ? this.padding.titlePad : 0);
  
  //Not sure where the -12 comes from, but things seem to work out...
  if( this.options.compactXAxis && xtitleh > 0 ) {
    var txth = 4 + this.padding.xTitlePad + xtitleh;
    this.padding.bottomComputed = this.padding.bottom + Math.max( txth, xlabelh-10 );
  } else {
    this.padding.bottomComputed = -12 + this.padding.bottom + xlabelh + (xtitleh > 0 ? this.padding.xTitlePad : 0) + xtitleh; 
  }
  
  
  // console.log("height beore")
  this.size.height = this.cy - this.padding.topComputed - this.padding.bottomComputed;
  this.size.sliderChartHeight = this.size.height / 10;
 
  this.calcLeftPadding( false );
 
  if (!this.origWidth)
    this.origWidth = this.size.width;
  this.size.width = this.cx - this.padding.leftComputed - this.padding.right;
  this.size.height = this.cy - this.padding.topComputed - this.padding.bottomComputed;

  this.xScale.range([0, this.size.width]);
  this.vis.attr("width",  this.size.width)
          .attr("height", this.size.height);
  this.vis.attr("transform", "translate(" + this.padding.leftComputed + "," + this.padding.topComputed + ")");
  this.plot.attr("width", this.size.width)
           .attr("height", this.size.height);

  d3.select(this.chart).select('svg')
      .attr("width", this.cx)
      .attr("height", this.cy)
      .attr("viewBox", "0 0 "+this.cx+" "+this.cy);

  this.vis.attr("width", this.cx )
          .attr("height", this.cy );

  d3.select("#chartarea"+this.chart.id)
          .attr("width", this.size.width )
          .attr("height", this.size.height );

  d3.select("#clip"+ this.chart.id).select("rect")
    .attr("width", this.size.width )
    .attr("height", this.size.height );

  //Fix the text position
  d3.select('svg').selectAll(".title")
      .attr("x", this.cx/2);

      // Christian: To keep the title from shifting during redraw, we translate the position of the label to stay in its place.
      // .attr("transform", "translate(" + ((this.size.width/2) - (this.origWidth/2)) + "," + 0 + ")");

  if( this.xaxistitle ){
    if( !this.options.compactXAxis ){
      this.xaxistitle.attr("x", this.cx/2);
      this.xaxistitle.attr("dy",xlabelh+this.padding.xTitlePad + (this.options.title ? 33 : 6)); // Christian: hard-coded 33 to account for padding between x-axis labels and slider chart

    } else {
      this.xaxistitle.attr("x", this.size.width - this.xaxistitle.node().getBBox().width + 30 );  //Not certain why this 30 is needed (I'm probably leaving out some element) but necassary to have label line up right
      this.xaxistitle.attr("dy",xtitleh + this.padding.xTitlePad + (this.options.title ? 31 : 4));
    }
    this.xaxistitle.attr("y", this.size.height);
  }

 // this.vis.selectAll(".yaxistitle")
    // .attr("transform","translate(" + -40 + " " + this.size.height/2+") rotate(-90)");

  this.xAxisBody.attr("width", this.size.width )
                .attr("height", this.size.height )
                .attr("transform", "translate(0," + this.size.height + ")")
                .call(this.xAxis);

  this.yAxisBody.attr("height", this.size.height )
                .call(this.yAxis);
  
  this.refLineInfo.select("circle").attr("cy", this.size.height);
  this.mouseInfo.attr("transform","translate(" + this.size.width + "," + this.size.height + ")");

  if( this.xGrid ) {
    this.xGrid.innerTickSize(-this.size.height)
    this.xGridBody.attr("width", this.size.width )
        .attr("height", this.size.height )
        .attr("transform", "translate(0," + this.size.height + ")")
        .call(this.xGrid);
  }

  if( this.yGrid ) {
    this.yGrid.innerTickSize(-this.size.width)
    this.yGridBody.attr("width", this.size.width )
        .attr("height", this.size.height )
        .attr("transform", "translate(0,0)")
        .call(this.yGrid);
  }
  
  //Make sure the legend stays visible
  if( this.legend ) {
    var trans = d3.transform(this.legend.attr("transform")).translate;
    var bb = this.legendBox.node().getBBox();
    var legx = ((trans[0]+bb.width) > this.cx) ? (this.cx-bb.width) : trans[0];
    var legy = ((trans[1]+bb.height) > this.cy) ? (this.cy-bb.height) : trans[1];
    this.legend.attr("transform", "translate(" + Math.max(0,legx) + "," + Math.max(legy,0) + ")" );
  }

  // Make sure the scaler stays visible
  if( this.scalerWidget ) {
    var trans = d3.transform(this.scalerWidget.attr("transform")).translate;
    var bb = this.scalerWidgetBox.node().getBBox();
    var scalerx = ((trans[0]+bb.width) > this.cx) ? (this.cx-bb.width) : trans[0];
    var scalery = ((trans[1]+bb.height) > this.cy) ? (this.cy-bb.height) : trans[1];
    this.scalerWidget.attr("transform", "translate(" + Math.max(0,scalerx) + "," + Math.max(scalery,0) + ")" );
  }

  if (!self.sliderChartWidthFactor) self.sliderChartWidthFactor = d3.select('svg').node().getBBox().width - self.size.sliderChartWidth;
  if (!self.sliderChartHeightFactor) self.sliderChartHeightFactor = self.size.sliderChartHeight / d3.select('svg')[0][0].clientHeight;
  self.size.sliderChartWidth = d3.select('svg').node().getBBox().width - self.sliderChartWidthFactor;
  self.size.sliderChartHeight = self.sliderChartHeightFactor * d3.select('svg')[0][0].clientHeight;
  if (this.options.showXAxisSliderChart) { self.drawXAxisSliderChart(); } 
  else                                   { self.cancelXAxisSliderChart(); }
  
  this.xScale.range([0, this.size.width]);
  this.yScale.range([0, this.size.height]);

  // this.zoom.x(this.xScale);
  // this.zoom.y(this.yScale);

  if( !dontRedraw ) {
    this.setData( this.rawData, false );
    this.redraw()();
  }
}



SpectrumChartD3.prototype.xticks = function() {

  var ticks = [];

  //it so
  //  the x axis labels (hopefully) always line up nicely where we want them
  //  e.g. kinda like multiple of 5, 10, 25, 50, 100, etc.
  var EPSILON = 1E-3;

  var rendermin = this.xScale.domain()[0];
  var rendermax = this.xScale.domain()[1];
  var range = rendermax - rendermin;
  
  //ndigstart makes up for larger numbers taking more pixels to render, so
  // it keeps numbers from overlapping.  Bellow uses 15 to kinda rpresent the
  //  width of numbers in pixels.
  var ndigstart = rendermin > 0 ? Math.floor(Math.log10(rendermin)) : 1;
  ndigstart = Math.max(1,Math.min(ndigstart,5));
  
  var nlabel = Math.floor( this.size.width / (50 + 15*(ndigstart-1)) );

  
  var renderInterval;// = range / 10.0;
  var n = Math.pow(10, Math.floor(Math.log10(range/nlabel)));
  var msd = range / (n * nlabel);


  if( isNaN(n) || !isFinite(n) || nlabel<=0 || n<=0.0 ) { //JIC
    console.log( "n=" + n + ", nlabel=" + nlabel + ", range=" + range );
    return ticks;
  }

  var subdashes = 0;

  //console.log( "msd=" + msd + ", nlabel=" + nlabel + ", ndigstart=" + ndigstart );

  if (msd < 1.5)
  {
    subdashes = 2;
    renderInterval = 0.5*n;
  }else if (msd < 3.3)
  {
    subdashes = 5;
    renderInterval = 0.5*n;
  }else if (msd < 7)
  {
    subdashes = 5;
    renderInterval = n;
  }else
  {
    subdashes = 10;
    renderInterval = n;
  }

  var biginterval = subdashes * renderInterval;
  var startEnergy = biginterval * Math.floor((rendermin + 1.0E-15) / biginterval);

  if( startEnergy < (rendermin-EPSILON*renderInterval) )
    startEnergy += biginterval;

  for( var i = -Math.floor(Math.floor(startEnergy-rendermin)/renderInterval); ; ++i)
  {
    if( i > 5000 )  //JIC
      break;

    var v = startEnergy + renderInterval * i;

    if( (v - rendermax) > EPSILON * renderInterval )
      break;

    var t = "";
    if( i>=0 && (i % subdashes == 0) )
      t += v;

    ticks.push( {value: v, major: (i % subdashes == 0), text: t } );
  }//for( intervals to draw ticks for )

  return ticks;
}

SpectrumChartD3.prototype.yticks = function() {
  var ticks = [];
  var EPSILON = 1.0E-3;

  var formatYNumber = function(v) {
    //poorly simulating "%.3g" in snprintf
    //SHould get rid of so many regexs, and shorten code (shouldnt there be a builtin function to print to "%.3"?)
    var t;
    if( v >= 1000 || v < 0.1 )
    {
      t = v.toPrecision(3);
      t = t.replace(/\.0+e/g, "e").replace(/\.0+$/g, "");
      if( t.indexOf('.') > 0 )
        t = t.replace(/0+e/g, "e");
    } else {
      t = v.toFixed( Math.max(0,2-Math.floor(Math.log10(v))) );
      if( t.indexOf('.') > 0 )
        t = t.replace(/0+$/g, "");
      t = t.replace(/\.$/g, "");
    }

    return t;
  }

  var renderymin = this.yScale.domain()[1],
      renderymax = this.yScale.domain()[0],
      heightpx = this.size.height;
  var range = renderymax - renderymin;

  if( this.options.yscale === "lin" )
  {
    //px_per_div: pixels between major and/or minor labels.
    var px_per_div = 50;

    //nlabel: approx number of major + minor labels we would like to have.
    var nlabel = heightpx / px_per_div;

    //renderInterval: Inverse of how many large labels to place between powers
    //  of 10 (1 is none, 0.5 is is one).
    var renderInterval;

    //n: approx how many labels will be used
    var n = Math.pow(10, Math.floor(Math.log10(range/nlabel)));

    //msd: approx how many sub dashes we would expect there to have to be
    //     to satisfy the spacing of labels we want with the given range.
    var msd = range / nlabel / n;

    if( isNaN(n) || !isFinite(n) || nlabel<=0 || n<=0.0 ) { //JIC
      console.log( "n=" + n + ", nlabel=" + nlabel + ", range=" + range );
      return ticks;
    }
      var subdashes = 0;

      if( msd < 1.5 )
      {
        subdashes = 2;
        renderInterval = 0.5*n;
      }else if (msd < 3.3)
      {
        subdashes = 5;
        renderInterval = 0.5*n;
      }else if (msd < 7)
      {
        subdashes = 5;
        renderInterval = n;
      }else
      {
        subdashes = 10;
        renderInterval = n;
      }

      var biginterval = subdashes * renderInterval;
      var starty = biginterval * Math.floor((renderymin + 1.0E-15) / biginterval);

      if( starty < (renderymin-EPSILON*renderInterval) )
        starty += biginterval;

      for( var i = -Math.floor(Math.floor(starty-renderymin)/renderInterval); ; ++i)
      {
        if( i > 500 )  //JIC
          break;

        var v = starty + renderInterval * i;

        if( (v - renderymax) > EPSILON * renderInterval )
          break;

        var t = "";
        if( i>=0 && ((i % subdashes) == 0) )
          t += formatYNumber(v);

        var len = ((i % subdashes) == 0) ? true : false;

        ticks.push( {value: v, major: len, text: t } );
      }//for( intervals to draw ticks for )
    }//case Chart::LinearScale:

    if( this.options.yscale === "log" )
    {
      //Get the power of 10 just bellow or equal to rendermin.
      var minpower = (renderymin > 0.0) ? Math.floor( Math.log10(renderymin) ) : -1;

      //Get the power of 10 just above or equal to renderymax.  If renderymax
      //  is less than or equal to 0, set power to be 0.
      var maxpower = (renderymax > 0.0) ? Math.ceil( Math.log10(renderymax) ): 0;

      //Adjust minpower and maxpower
      if( maxpower == minpower )
      {
        //Happens when renderymin==renderymax which is a power of 10
        ++maxpower;
        --minpower;
      }else if( maxpower > 2 && minpower < -1)
      {
        //We had a tiny value (possibly a fraction of a count), as well as a
        //  large value (>1000).
        minpower = -1;
      }else if( maxpower >= 0 && minpower < -1 && (maxpower-minpower) > 6 )
      {
        //we had a tiny power (1.0E-5), as well as one between 1 and 999,
        //  so we will only show the most significant decades
        minpower = maxpower - 5;
      }//if( minpower == maxpower ) / else / else


      //numdecades: number of decades the data covers, including the decade
      //  above and bellow the data.
      var numdecades = maxpower - minpower + 1;

      //minpxdecade: minimum number of pixels we need per decade.
      var minpxdecade = 25;

      //labeldelta: the number of decades between successive labeled large ticks
      //  each decade will have a large tick regardless
      var labeldelta = 1;

      //nticksperdecade: number of small+large ticks per decade
      var nticksperdecade = 10;


      if( Math.floor(heightpx / minpxdecade) < numdecades )
      {
        labeldelta = numdecades / (heightpx / minpxdecade);
        nticksperdecade = 1;
      }//if( (heightpx / minpxdecade) < numdecades )

      var t = null;
      var nticks = 0;
      var nmajorticks = 0;

      for( var decade = minpower; decade <= maxpower; ++decade )
      {
        var startcounts = Math.pow( 10.0, decade );
        var deltacounts = 10.0 * startcounts / nticksperdecade;
        var eps = deltacounts * EPSILON;

        if( (startcounts - renderymin) > -eps && (startcounts - renderymax) < eps )
        {
          t = ((decade%labeldelta)==0) ? formatYNumber(startcounts) : "";
          ++nticks;
          ++nmajorticks;
          ticks.push( {value: startcounts, major: true, text: t } );
        }//if( startcounts >= renderymin && startcounts <= renderymax )


        for( var i = 1; i < (nticksperdecade-1); ++i )
        {
          var y = startcounts + i*deltacounts;
          if( (y - renderymin) > -eps && (y - renderymax) < eps )
          {
            ++nticks;  
            ticks.push( {value: y, major: false, text: null } );
          }
        }//for( int i = 1; i < nticksperdecade; ++i )
      }//for( int decade = minpower; decade <= maxpower; ++decade )

      //If we have a decent number of (sub) labels, the user can orient
      //  themselves okay, so well get rid of the minor labels.
      if( (nticks > 8 || (heightpx/nticks) < 25 || nmajorticks > 1) && nmajorticks > 0 ) {
        for( var i = 0; i < ticks.length; ++i )
          if( !ticks[i].major )
            ticks[i].label = "";
      }
      
      if( nmajorticks < 1 && ticks.length ) {
        ticks[0].text = formatYNumber(ticks[0].value);
        ticks[ticks.length-1].text = formatYNumber(ticks[ticks.length-1].value);
      }
      
      if( ticks.length === 0 )
      {
        //cerr << "Forcing a single axis point in" << endl;
        var y = 0.5*(renderymin+renderymax);
        var t = formatYNumber(y);
        ticks.push( {value: y, major: true, text: t } );
      }
    }//case Chart::LogScale:

  //TODO: should to properly implement sqrt
  if( this.options.yscale === "sqrt" ){
    //Take the easy way out for now until I can get around to customizing
    this.yScale.copy().ticks(10)
        .forEach(function(e){ticks.push({value:e,major:true,text:formatYNumber(e)});});
  }


  return ticks;
}

SpectrumChartD3.prototype.updateMouseCoordText = function() {
  var self = this;

  //without d3.event
  if( !d3.event )
    return;
  if ( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length )
    return;

  var p = d3.mouse(self.vis[0][0]);
  // console.log("fix mouse text");

  if( !p ){
    p = d3.touch(self.vis[0][0]);
    p =  (p && p.length===1) ? p[0] : null;
  }

  if( !p ){
    self.mousedOverRefLine = null;
    self.refLineInfo.style("display", "none");
    self.mouseInfo.style("display", "none");
    return;
  }

  var energy = self.xScale.invert(p[0]);
  var y = self.yScale.invert(p[1]);

  //Find what channel this energy cooresponds
  var channel, lowerchanval, counts = null;
  var foreground = self.rawData.spectra[0];
  channel = (foreground.x && foreground.x.length) ? d3.bisector(function(d){return d;}).right(foreground.x, energy) : -1;
  if( foreground.points && foreground.points.length ){
    lowerchanval = d3.bisector(function(d){return d.x;}).left(foreground.points,energy,1) - 1;
    counts = foreground.points[lowerchanval].y; 
  }

  //console.log( "updateMouseCoordText: {" + energy + "," + y + "} and also fix self.refLines[0] usage");

  //Currently two majorish issues with the mouse text
  //  1) The counts for the foreground are all that is given. Should give for background/second
  //  2) Currently gives a channel (singular), but when there this.rebinFactor!=1, should give range
  //Also, formatting could be a bit better
  //Also need to make displaying this box an option
  //Could also add a blue dot in the data or something, along with lines to the axis, both features should be optional

  //Right now if mouse stats arent deesired we are just not showing them, but
  //  still updating them...

  if( self.options.showMouseStats ){
    self.mouseInfo.style("display", null );
    var mousetxt = self.mouseInfo.select("text");
    mousetxt.attr('dy', "-2em");
    mousetxt.selectAll("tspan").remove();

    var xmmsg = ""+(Math.round(10*energy)/10) + " keV";
    if( channel )
      xmmsg += ", chan: " + channel;
    var ymmsg = "";
    if( counts !== null )
        ymmsg += "counts: " + (Math.round(10*counts)/10) + (foreground.rebinFactor === 1 ? "" : ("/" + foreground.rebinFactor));
    ymmsg += (counts!==null?", ":"") + "y: " + (Math.round(10*y)/10);


    var ymmsglen = mousetxt.append('svg:tspan').attr('x', "-12em").attr('dy', "-1em")
                    .text( ymmsg ).node().getComputedTextLength();
    var xmmsglen = mousetxt.append('svg:tspan').attr('x', "-12em").attr('dy', "-1em")
            .text( xmmsg ).node().getComputedTextLength();

    //Resize the box to match the text size
    self.mouseInfoBox.attr('width', Math.max(xmmsglen,ymmsglen) + 9 );
  } else {
    self.mouseInfo.style("display", "none" );
  }

  var mindist = 9.0e20, nearestpx = 9.0e20;

  var nearestline = null;
  var h = self.size.height;
  var m = Math.min(h,self.options.refLineTopPad);
  var visy = Math.max(h-m,1);

  var reflines = d3.selectAll("g.ref");

  reflines[0].forEach( function(d,i){
    var yh = d.childNodes[0].y1.baseVal.value - d.childNodes[0].y2.baseVal.value;
    //var xpx = d.transform.baseVal[0].matrix.e;
    var xpx = self.xScale(d.__data__.e);
    var dpx = Math.abs(xpx - p[0]);

    //In principle, this check (or __data__.mousedover) shouldnt be necassary,
    //  but with out it sometimes lines will stay fat that arent supposed to.
    // Also, I think setting attr values is expensive, so only doing if necassary.
    if( d.__data__.mousedover && d !== self.mousedOverRefLine ){
      d3.select(d).attr("stroke-width",1).attr("dx", "-0.5" );
      d.__data__.mousedover = null;
    }

    var dist = dpx + dpx/(yh/visy);
    if( dist < mindist ) {
      mindist = dist;
      nearestline = d;
      nearestpx = dpx;
    }
  } );

  if( nearestpx > 10 ) {
    if( self.mousedOverRefLine ){
      d3.select(self.mousedOverRefLine).attr("stroke-width",1).attr("dx", "-0.5" );
      self.mousedOverRefLine.__data__.mousedover = null;
    }
    self.mousedOverRefLine = null;
    self.refLineInfo.style("display", "none");
    return;
  }

  if( self.mousedOverRefLine===nearestline )
    return;

  if( self.mousedOverRefLine ){
    d3.select(self.mousedOverRefLine).attr("stroke-width",1).attr("dx", "-0.5" );
    self.mousedOverRefLine.__data__.mousedover = null;
  }

  self.mousedOverRefLine = nearestline;

  var linedata = nearestline.__data__;
  linedata.mousedover = true;
  var e = linedata.e;
  var sf = linedata.h;
  var linepx = self.xScale(e);

  var txt, textdescrip, attTxt;
  var detector = linedata.parent.detector
  var shielding = linedata.parent.shielding;
  var shieldingThickness = linedata.parent.shieldingThickness;
  var nearestLineParent = linedata.parent.parent;


  //20160308: self.refLines[0] should now contain member variable particleSf
  //  that gives a map from the particle name to the scale factor applied to
  //  get sf. ex: var abs_br = self.refLines[0].particleSf[linedata.particle] * sf.

  textdescrip = (nearestLineParent ? (nearestLineParent + ', ') : "") +  e + ' keV, rel. amp. ' + sf;

  if( linedata.decay ) {
    if( linedata.particle === 'gamma' ) {
      txt = linedata.decay;
      if( linedata.decay.indexOf('Capture') < 0 )
        txt += ' decay';
    }else if( linedata.particle === 'xray' ) {
      textdescrip = linedata.el + ' x-ray, ' +  e + ' keV, I=' + sf;
    }else if( (typeof linedata.particle === 'string') && linedata.particle.length ){
      textdescrip = linedata.particle + ' from ' + nearestLineParent + ", " +  e + ' keV, I=' + sf;
    }
  }

  if( linedata.particle === 'gamma' || linedata.particle === 'xray' ) {
    if( shielding && shieldingThickness )
      attTxt = shieldingThickness + ' of ' + shielding;
    if( detector )
      attTxt = (attTxt ? (attTxt + ' with a ' + detector) : 'Assuming a ' + detector);
  }

  self.refLineInfo.style("display", null).attr("transform", "translate("+linepx+",0)" );
  var svgtxt = self.refLineInfoTxt.select("text")
                   .attr("dy", "1em")
                   .attr("fill", linedata.parent.color);
  svgtxt.selectAll("tspan").remove();
  svgtxt.append('svg:tspan').attr('x', 0).attr('dy', "1em").text( textdescrip );
  if( txt )
    svgtxt.append('svg:tspan').attr('x', 0).attr('dy', "1em").text( txt );
  if( attTxt )
    svgtxt.append('svg:tspan').attr('x', 0).attr('dy', "1em").text( attTxt );

  //Now detect if svgtxt is running off the right side of the chart, and if so
  //  put to the left of the line
  var tx = function(d) {
    var w = this.getBBox().width;
    return ((linepx+5+w)>self.size.width) ? ("translate("+(-5-w)+",0)") : "translate(5,0)";
  };
  self.refLineInfoTxt.attr("transform", tx );


  d3.select(nearestline).attr("stroke-width",2).select("line").attr("dx", "-1" );
}

SpectrumChartD3.prototype.updateFeatureMarkers = function(sumPeaksArgument) {
  var self = this;

  /* Christian: Just catching this weird error whenever user initially clicks the checkbox for one of the feature markers.
                I'm guessing it is thrown because the mouse does not exist for the vis element, so we don't update the feature markers.
   */
  try { d3.mouse(self.vis[0][0]) }                        
  catch (error) { console.log( "Error thrown: source event is null" ); return; }  // to catch and log certain error thrown in initial check of feature markers

  // Positional variables (for mouse and touch)
  var m = d3.mouse(self.vis[0][0]),
      t = d3.touches(self.vis[0][0]);

  // Adjust the mouse position accordingly to touch (because some of these functions use mouse position in touch devices)
  if (t.length > 0)
    m = t[0];

  // Chart coordinate values
  var energy = self.xScale.invert(m[0]),
      xmax = self.size.width,
      ymax = self.size.height;

  // Do not update feature markers if legend being dragged or energy value is undefined
  if ((t && self.legdown) || (t && self.scalerdown) || isNaN(energy))
    return;
  if (self.currentlyAdjustingSpectrumScale)
    return;

  var cursorIsOutOfBounds = (t && t.length > 0) ? (t[0][0] < 0 || t[0][0] > xmax) : (m[0] < 0  || m[0] > xmax || m[1] < 0 || m[1] > ymax);


  /* Mouse-edge Helpers: These are global helpers for feature markers that update the mouse edge position.
  */

  /* Mouse edge should be deleted if: 
      none of the scatter/escape peak options are unchecked 
      OR if cursor is out of bounds 
      OR if the user is currently dragging in the graph
  */
  function shouldDeleteMouseEdge() {
    return (!self.options.showComptonEdge && !self.options.showComptonPeaks && !self.options.showEscapePeaks && !self.options.showSumPeaks) || 
    cursorIsOutOfBounds || 
    self.dragging_plot;
  }
  function deleteMouseEdge( shouldDelete ) {
    if ( shouldDelete || shouldDeleteMouseEdge() )
    {
      if ( self.mouseEdge ) {
        self.mouseEdge.remove(); 
        self.mouseEdge = null;
      }
      if ( self.mouseEdgeText ) {
        self.mouseEdgeText.remove(); 
        self.mouseEdgeText = null;
      }
      return true;
    }
    return false;
  }
  function updateMouseEdge() {

    // If mouse edge could has been deleted, do not update the mouse edge
    if (deleteMouseEdge())
      return;

    // Update the mouse edge and corresponding text position 
    if ( self.mouseEdge ) {
        self.mouseEdge
          .attr("x1", m[0])
          .attr("x2", m[0]);
        self.mouseEdgeText
          .attr("x", m[0] + xmax/125 )
          .text( energy.toFixed(1) + " keV");
    } else {  
        // Create the mouse edge (and text next to it)
        self.mouseEdge = self.vis.append("line")
          .attr("class", "mouseLine")
          .attr("stroke-width", 2)
          .attr("stroke", "black")
          .attr("x1", m[0])
          .attr("x2", m[0])
          .attr("y1", 0)
          .attr("y2", self.size.height);

        self.mouseEdgeText = self.vis.append("text")
          .attr("class", "mouseLineText")
          .attr( "x", m[0] + xmax/125 )
          .attr( "y", self.size.height/4)
          .text( energy.toFixed(1) + " keV");
    }
  }


  function updateEscapePeaks() {
    // Calculations for the escape peak markers
    var singleEscapeEnergy = energy - 510.99891,
        singleEscapePix = self.xScale(singleEscapeEnergy),

        singleEscapeForwardEnergy = energy + 510.99891,
        singleEscapeForwardPix = self.xScale(singleEscapeForwardEnergy),

        doubleEscapeEnergy = energy - 1021.99782,
        doubleEscapePix = self.xScale(doubleEscapeEnergy),

        doubleEscapeForwardEnergy = energy + 1021.99782,
        doubleEscapeForwardPix = self.xScale(doubleEscapeForwardEnergy);

    // Deletes the marker for a single escape peak
    function deleteSingleEscape() {
      if( self.singleEscape ) {
        self.singleEscape.remove();
        self.singleEscape = null;
      }
      if ( self.singleEscapeText ) {
          self.singleEscapeText.remove();
          self.singleEscapeText = null;
      }
      if ( self.singleEscapeMeas ) {
          self.singleEscapeMeas.remove();
          self.singleEscapeMeas = null;
      }
    }

    // Deletes the marker for a double escape peak
    function deleteDoubleEscape() {
      if ( self.doubleEscape ) {
        self.doubleEscape.remove();
        self.doubleEscape = null;
      }
      if ( self.doubleEscapeText ) {
        self.doubleEscapeText.remove();
        self.doubleEscapeText = null;
      }
      if ( self.doubleEscapeMeas ) {
        self.doubleEscapeMeas.remove();
        self.doubleEscapeMeas = null;
      }
    }
    // Deletes the marker for a single forward escape peak
    function deleteSingleEscapeForward() {
      if ( self.singleEscapeForward ) {
        self.singleEscapeForward.remove();
        self.singleEscapeForward = null;
      }
      if ( self.singleEscapeForwardText ) {
          self.singleEscapeForwardText.remove();
          self.singleEscapeForwardText = null;
      }
      if ( self.singleEscapeForwardMeas ) {
          self.singleEscapeForwardMeas.remove();
          self.singleEscapeForwardMeas = null;
      }
    }

    // Deletes the marker for a double forward escape peak
    function deleteDoubleEscapeForward() {
      if ( self.doubleEscapeForward ) {
        self.doubleEscapeForward.remove();
        self.doubleEscapeForward = null;
      }
      if ( self.doubleEscapeForwardText ) {
          self.doubleEscapeForwardText.remove();
          self.doubleEscapeForwardText = null;
      }
      if ( self.doubleEscapeForwardMeas ) {
          self.doubleEscapeForwardMeas.remove();
          self.doubleEscapeForwardMeas = null;
      }
    }

    if (shouldDeleteMouseEdge())
      deleteMouseEdge(true);
    
    if( !self.options.showEscapePeaks || cursorIsOutOfBounds || self.dragging_plot ) {
      deleteSingleEscape();
      deleteDoubleEscape();
      deleteSingleEscapeForward();
      deleteDoubleEscapeForward();
      return;
    }

    var singleEscapeOutOfBounds = singleEscapePix < 0 || singleEscapePix > xmax,
        doubleEscapeOutOfBounds = doubleEscapePix < 0 || doubleEscapePix > xmax,
        singleEscapeForwardOutOfBounds = singleEscapeForwardPix < 0 || singleEscapeForwardPix > xmax,
        doubleEscapeForwardOutOfBounds = doubleEscapeForwardPix < 0 || doubleEscapeForwardPix > xmax;

    if ( doubleEscapeOutOfBounds ) {
      deleteDoubleEscape();

      if ( singleEscapeOutOfBounds )
        deleteSingleEscape();
    }

    if ( doubleEscapeForwardOutOfBounds ) {
      deleteDoubleEscapeForward();

      if ( singleEscapeForwardOutOfBounds )
        deleteSingleEscapeForward();
    }

    updateMouseEdge();

    if (!singleEscapeForwardOutOfBounds) {
      if( !self.singleEscapeForward && singleEscapeForwardEnergy >= 0 ) {  
        self.singleEscapeForward = self.vis.append("line")    // create single forward escape line
        .attr("class", "escapeLineForward")
        .attr("x1", singleEscapeForwardPix)
        .attr("x2", singleEscapeForwardPix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.singleEscapeForwardText = self.vis.append("text") // create Single Forward Escape label beside line
            .attr("class", "peakText")
            .attr( "x", singleEscapeForwardPix + xmax/200 )
            .attr( "y", self.size.height/5.3)
            .text( "Single Escape" );
      self.singleEscapeForwardMeas = self.vis.append("text") // Create measurement label besides line, under Single Escape label
            .attr("class", "peakText")
            .attr( "x", singleEscapeForwardPix + xmax/125 )
            .attr( "y", self.size.height/4.5)
            .text( singleEscapeForwardEnergy.toFixed(1) + " keV" );
      } else {
        if ( singleEscapeForwardEnergy < 0 && self.singleEscapeForward && self.singleEscapeForwardText && self.singleEscapeForwardMeas ) {
          deleteSingleEscapeForward();

        } else if ( self.singleEscapeForward ) {      // Move everything to where mouse is
          self.singleEscapeForward
            .attr("x1", singleEscapeForwardPix)
            .attr("x2", singleEscapeForwardPix);
          self.singleEscapeForwardText
            .attr("x", singleEscapeForwardPix + xmax/200 );
          self.singleEscapeForwardMeas
            .attr( "x", singleEscapeForwardPix + xmax/125 )
            .text( singleEscapeForwardEnergy.toFixed(1) + " keV" );
        }
      }
    }

    if (!doubleEscapeForwardOutOfBounds) {
      if( !self.doubleEscapeForward && doubleEscapeForwardEnergy >= 0 ) {  
        self.doubleEscapeForward = self.vis.append("line")    // create double forward escape line
        .attr("class", "escapeLineForward")
        .attr("x1", doubleEscapeForwardPix)
        .attr("x2", doubleEscapeForwardPix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.doubleEscapeForwardText = self.vis.append("text") // create double Forward Escape label beside line
            .attr("class", "peakText")
            .attr( "x", doubleEscapeForwardPix + xmax/200 )
            .attr( "y", self.size.height/5.3)
            .text( "Double Escape" );
      self.doubleEscapeForwardMeas = self.vis.append("text") // Create measurement label besides line, under double Escape label
            .attr("class", "peakText")
            .attr( "x", doubleEscapeForwardPix + xmax/125 )
            .attr( "y", self.size.height/4.5)
            .text( doubleEscapeForwardEnergy.toFixed(1) + " keV" );
      } else {
        if ( doubleEscapeForwardEnergy < 0 && self.doubleEscapeForward && self.doubleEscapeForwardText && self.doubleEscapeForwardMeas ) {
          deleteDoubleEscapeForward();

        } else if ( self.doubleEscapeForward ) {      // Move everything to where mouse is
          self.doubleEscapeForward
            .attr("x1", doubleEscapeForwardPix)
            .attr("x2", doubleEscapeForwardPix);
          self.doubleEscapeForwardText
            .attr("x", doubleEscapeForwardPix + xmax/200 );
          self.doubleEscapeForwardMeas
            .attr( "x", doubleEscapeForwardPix + xmax/125 )
            .text( doubleEscapeForwardEnergy.toFixed(1) + " keV" );
        }
      }
    }

    if ( singleEscapeOutOfBounds ) 
      return;

    // Draw single escape peak if not present in the grapy
    if( !self.singleEscape && singleEscapeEnergy >= 0 ) {  
      self.singleEscape = self.vis.append("line")    // create single escape line
      .attr("class", "peakLine")
      .attr("x1", singleEscapePix)
      .attr("x2", singleEscapePix)
      .attr("y1", 0)
      .attr("y2", self.size.height);
    self.singleEscapeText = self.vis.append("text") // create Single Escape label beside line
          .attr("class", "peakText")
          .attr( "x", singleEscapePix + xmax/200 )
          .attr( "y", self.size.height/5.3)
          .text( "Single Escape" );
    self.singleEscapeMeas = self.vis.append("text") // Create measurement label besides line, under Single Escape label
          .attr("class", "peakText")
          .attr( "x", singleEscapePix + xmax/125 )
          .attr( "y", self.size.height/4.5)
          .text( singleEscapeEnergy.toFixed(1) + " keV" );
    } else {
      if ( singleEscapeEnergy < 0 && self.singleEscape && self.singleEscapeText && self.singleEscapeMeas ) {
        self.singleEscape.remove();      // Delete lines and labels if out of bounds
        self.singleEscape = null;
        self.singleEscapeText.remove();
        self.singleEscapeText = null;
        self.singleEscapeMeas.remove();
        self.singleEscapeMeas = null;

      } else if ( self.singleEscape ) {      // Move everything to where mouse is
        self.singleEscape
          .attr("x1", singleEscapePix)
          .attr("x2", singleEscapePix);
        self.singleEscapeText
          .attr("x", singleEscapePix + xmax/200 );
        self.singleEscapeMeas
          .attr( "x", singleEscapePix + xmax/125 )
          .text( singleEscapeEnergy.toFixed(1) + " keV" );
      }
    }

    // Do not update the double escape peak marker anymore
    if (doubleEscapeOutOfBounds) 
      return;

      // Draw double escape peak if not present in the grapy
    if( !self.doubleEscape && doubleEscapeEnergy >= 0 ) {
      self.doubleEscape = self.vis.append("line")  // create double escape line
      .attr("class", "peakLine")
      .attr("stroke-width", 2)
      .attr("stroke", "black")
      .attr("x1", doubleEscapePix)
      .attr("x2", doubleEscapePix)
      .attr("y1", 0)
      .attr("y2", self.size.height);
    self.doubleEscapeText = self.vis.append("text") // create Double Escape label beside line
          .attr("class", "peakText")
          .attr( "x", doubleEscapePix + xmax/200 )
          .attr( "y", self.size.height/5.3)
          .text( "Double Escape" );
    self.doubleEscapeMeas = self.vis.append("text") // Create measurement label besides line, under Double Escape label
          .attr("class", "peakText")
          .attr( "x", doubleEscapePix + xmax/125 )
          .attr( "y", self.size.height/4.5)
          .text( doubleEscapeEnergy.toFixed(1) + " keV" );
    } else {
      if ( (doubleEscapeEnergy < 0) && self.doubleEscape && self.doubleEscapeText && self.doubleEscapeMeas ) {
        self.doubleEscape.remove();    // Delete lines and labels if out of bounds
      self.doubleEscape = null;
          self.doubleEscapeText.remove();
          self.doubleEscapeText = null;
          self.doubleEscapeMeas.remove();
          self.doubleEscapeMeas = null;

      } else if ( self.doubleEscape ) {    // Move everything to where mouse is

        self.doubleEscape
          .attr("x1", doubleEscapePix)
          .attr("x2", doubleEscapePix);
        self.doubleEscapeText
          .attr("x", doubleEscapePix + xmax/200 );
        self.doubleEscapeMeas
          .attr( "x", doubleEscapePix + xmax/125 )
          .text( doubleEscapeEnergy.toFixed(1) + " keV" );
      }
    }
  }

  function updateComptonPeaks() {

    var compAngleRad = self.comptonPeakAngle * (3.14159265/180.0)   // calculate radians of compton peak angle
    var comptonPeakEnergy = energy / (1 + ((energy/510.99891)*(1-Math.cos(compAngleRad)))); // get energy value from angle and current energy position
    var comptonPeakPix = self.xScale(comptonPeakEnergy);

    if (shouldDeleteMouseEdge())
      deleteMouseEdge(true);
    
    var comptonPeakOutOfBounds = comptonPeakPix < 0 || comptonPeakPix > xmax;

    // delete if compton peak option is turned off or cursor is out of the graph
    if( !self.options.showComptonPeaks || cursorIsOutOfBounds || comptonPeakOutOfBounds || self.dragging_plot ) {
      if( self.comptonPeak ) {
        self.comptonPeak.remove(); self.comptonPeak = null;
        // console.log( 'Should emit compton peak closed' );
      }
      if ( self.comptonPeakText ) {
        self.comptonPeakText.remove(); self.comptonPeakText = null;
      }
      if ( self.comptonPeakMeas ) {
        self.comptonPeakMeas.remove(); self.comptonPeakMeas = null;
      }

      if (t.length == 0) {        // if mouse movement (aka not touch movement), then update the mouse edge
        updateMouseEdge();
      }

      return;
    }

    if( !self.comptonPeak ) {
      // draw compton edge line here
      self.comptonPeak = self.vis.append("line")
        .attr("class", "peakLine")
        .attr("stroke-width", 2)
        .attr("stroke", "black")
        .attr("x1", comptonPeakPix)
        .attr("x2", comptonPeakPix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.comptonPeakText = self.vis.append("text")
        .attr("class", "peakText")
        .attr( "x", comptonPeakPix + xmax/200 )
        .attr( "y", self.size.height/10)
        .text( self.comptonPeakAngle + "° Compton Peak" );
      self.comptonPeakMeas = self.vis.append("text")
        .attr("class", "peakText")
        .attr( "x", comptonPeakPix + xmax/125 )
        .attr( "y", self.size.height/7.8)
        .text( comptonPeakEnergy.toFixed(1) + " keV" );
    } else {
      // update the compton peak edge line
      self.comptonPeak
        .attr("x1", comptonPeakPix)
        .attr("x2", comptonPeakPix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.comptonPeakText
        .attr("x", comptonPeakPix + xmax/200 )
        .attr( "y", self.size.height/10)
        .text( self.comptonPeakAngle + "° Compton Peak" );
      self.comptonPeakMeas
        .attr( "x", comptonPeakPix + xmax/125 )
        .attr( "y", self.size.height/7.8)
        .text( comptonPeakEnergy.toFixed(1) + " keV" );
    }
    updateMouseEdge();
  }

  function updateComptonEdge() {

    var compedge = energy - (energy / (1 + (2*(energy/510.99891))));
    var compEdgePix = self.xScale(compedge);
    
    if ( shouldDeleteMouseEdge() ) 
      deleteMouseEdge(true);

    var comptonEdgeOutOfBounds = compEdgePix < 0  || compEdgePix > xmax ;

    // delete if compton edge already if option is turned off or cursor is out of the graph
    if( !self.options.showComptonEdge || cursorIsOutOfBounds || comptonEdgeOutOfBounds || self.dragging_plot ) {
      if( self.comptonEdge ) {
        self.comptonEdge.remove(); 
        self.comptonEdge = null;
      }
      if ( self.comptonEdgeText ) {
        self.comptonEdgeText.remove(); 
        self.comptonEdgeText = null;
      }
      if ( self.comptonEdgeMeas ) {
          self.comptonEdgeMeas.remove(); 
          self.comptonEdgeMeas = null;
      }

      updateMouseEdge();
      return;
    }
    
    if( !self.comptonEdge ) {
      // draw compton edge line here
      self.comptonEdge = self.vis.append("line")
        .attr("class", "peakLine")
        .attr("stroke-width", 2)
        .attr("stroke", "black")
        .attr("x1", compEdgePix)
        .attr("x2", compEdgePix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.comptonEdgeText = self.vis.append("text")
        .attr("class", "peakText")
        .attr( "x", compEdgePix + xmax/200 )
        .attr( "y", self.size.height/22)
        .text( "Compton Edge" );
      self.comptonEdgeMeas = self.vis.append("text")
        .attr("class", "peakText")
        .attr( "x", compEdgePix + xmax/125 )
        .attr( "y", self.size.height/14)
        .text( compedge.toFixed(1) + " keV" );
    } else {
      self.comptonEdge
        .attr("x1", compEdgePix)
        .attr("x2", compEdgePix)
        .attr("y1", 0)
        .attr("y2", self.size.height);
      self.comptonEdgeText
        .attr("x", compEdgePix + xmax/200 )
        .attr("y", self.size.height/22);
      self.comptonEdgeMeas
        .attr( "x", compEdgePix + xmax/125 )
        .attr( "y", self.size.height/14)
        .text( compedge.toFixed(1) + " keV" );
    }
    updateMouseEdge();
  }

  function updateSumPeaks( clickedEnergy ) {

    function deleteClickedPeakMarker() {
      if( self.clickedPeak ) {
        self.clickedPeak.remove(); 
        self.clickedPeak = null;
      }
      if ( self.clickedPeakMeas ) {
        self.clickedPeakMeas.remove(); 
        self.clickedPeakMeas = null;
      }
    }

    function deleteSumPeakMarker() {
      if ( self.sumPeak ) {
        self.sumPeak.remove(); 
        self.sumPeak = null;
      }
      if ( self.sumPeakMeas ) {
        self.sumPeakMeas.remove(); 
        self.sumPeakMeas = null;
      }
      if ( self.sumPeakText ) {
          self.sumPeakText.remove(); 
          self.sumPeakText = null;
      }
    }

    function deleteLeftSumPeakMarker() {
      if ( self.leftSumPeak ) {
        self.leftSumPeak.remove(); 
        self.leftSumPeak = null;
      }
      if ( self.leftSumPeakMeas ) {
        self.leftSumPeakMeas.remove(); 
        self.leftSumPeakMeas = null;
      }
      if ( self.leftSumPeakText ) {
        self.leftSumPeakText.remove(); 
        self.leftSumPeakText = null;
      }
    }
    
    if ( shouldDeleteMouseEdge() ) 
      deleteMouseEdge(true);

    // delete if sum peak option is already turned off or cursor is out of the graph
    if( !self.options.showSumPeaks || cursorIsOutOfBounds || self.dragging_plot ) {
      // delete the sum peak corresponding help text 
      if ( self.sumPeakHelpText ) {
        self.sumPeakHelpText.remove(); 
        self.sumPeakHelpText = null;
      }

      if ( cursorIsOutOfBounds || self.dragging_plot ) {
        if ( self.clickedPeak ) {
          self.savedClickEnergy = self.xScale.invert( self.clickedPeak.attr("x1") );
        }
      } else if ( !self.options.showSumPeaks )
        self.savedClickEnergy = null;

      deleteClickedPeakMarker();
    }
    

    if( !self.options.showSumPeaks || cursorIsOutOfBounds || self.dragging_plot ) {
      deleteSumPeakMarker();
      deleteLeftSumPeakMarker();
      return;
    }

    if ( !self.options.showSumPeaks ) 
      self.savedClickEnergy = null;

    updateMouseEdge();

    var shouldUpdateClickedPeak = true,
        shouldUpdateSumPeak = true,
        shouldUpdateLeftSumPeak = true;

    if ( self.savedClickEnergy == null && clickedEnergy == null ) {
      shouldUpdateSumPeaks = shouldUpdateLeftSumPeak = shouldUpdateClickedPeak = false;
    }
    else if ( self.savedClickEnergy != null && (clickedEnergy == null || clickedEnergy < 0) )
      clickedEnergy = self.savedClickEnergy;
    else if ( clickedEnergy < 0 ) {
      if ( self.clickedPeak ) {
        self.savedClickEnergy = clickedEnergy = Number(self.clickedPeak.attr("energy")) ;
      }
      else shouldUpdateClickedPeak = false;
    }

    if ( !shouldUpdateClickedPeak && !shouldUpdateSumPeak && !shouldUpdateLeftSumPeak )
      return;

    var clickedEdgeOutOfBounds = false;
    if ( shouldUpdateClickedPeak ) {
      self.savedClickEnergy = null;

      var clickedEdgePix = self.xScale( clickedEnergy  );
      clickedEdgeOutOfBounds = clickedEdgePix < 0 || clickedEdgePix > xmax;

      if (clickedEdgeOutOfBounds) {
        self.savedClickEnergy = clickedEnergy;
        deleteClickedPeakMarker();
      }
      else if( !self.clickedPeak ) {
        // draw compton edge line here
        self.clickedPeak = self.vis.append("line")
          .attr("class", "peakLine")
          .attr("stroke-width", 2)
          .attr("stroke", "black")
          .attr("x1", clickedEdgePix)
          .attr("x2", clickedEdgePix)
          .attr("y1", 0)
          .attr("y2", self.size.height)
          .attr("energy", clickedEnergy);
        self.clickedPeakMeas = self.vis.append("text")
          .attr("class", "peakText")
          .attr( "x", clickedEdgePix + xmax/125 )
          .attr( "y", self.size.height/4)
          .text( clickedEnergy.toFixed(1) + " keV" );
      } else {

        self.clickedPeak
          .attr("x1", clickedEdgePix)
          .attr("x2", clickedEdgePix)
          .attr("y1", 0)
          .attr("y2", self.size.height)
          .attr("energy", clickedEnergy);
        self.clickedPeakMeas
          .attr( "x", clickedEdgePix + xmax/125 )
          .attr( "y", self.size.height/4)
          .text( clickedEnergy.toFixed(1) + " keV" );
      }
    }  

    if ( !self.clickedPeak && !self.savedClickEnergy ) { 
      shouldUpdateSumPeak = false;

      if ( !self.sumPeakHelpText ) {
        // create the sum peak help text
        self.sumPeakHelpText = self.vis.append("text")
          .attr("class", "peakText")
          .attr("fill", "red")
          .attr( "x", m[0] + xmax/125 )
          .attr( "y", self.size.height/3.5)
              .text( "Click to set sum peak first energy." );
      } else {
        // update the sum peak help text position
        self.sumPeakHelpText
          .attr( "x", m[0] + xmax/125 )
          .attr( "y", self.size.height/3.5)
      }

    } else {
        // delete sum peak help text
        if ( self.sumPeakHelpText ) {
          self.sumPeakHelpText.remove(); 
          self.sumPeakHelpText = null;
        }
    }

    if ( shouldUpdateLeftSumPeak && energy < clickedEnergy ) {
      var leftSumEnergy = clickedEnergy - energy,
          leftSumPix = self.xScale( leftSumEnergy  ),
          leftSumOutOfBounds = leftSumPix < 0 || leftSumPix > xmax;

      if( !self.leftSumPeak ) {
        // draw left-sum peak line here
        self.leftSumPeak = self.vis.append("line")
          .attr("class", "peakLine")
          .attr("stroke-width", 2)
          .attr("stroke", "black")
          .attr("x1", leftSumPix)
          .attr("x2", leftSumPix)
          .attr("y1", 0)
          .attr("y2", self.size.height);
        self.leftSumPeakText = self.vis.append("text")
          .attr("class", "peakText")
          .attr( "x", leftSumPix + xmax/200 )
          .attr( "y", self.size.height/3.4)
          .text( "Clicked Peak" );
        self.leftSumPeakMeas = self.vis.append("text")
          .attr("class", "peakText")
          .attr( "x", leftSumPix + xmax/125 )
          .attr( "y", self.size.height/3.0)
          .text( energy.toFixed(1) + "+" + leftSumEnergy.toFixed(1) + "=" + clickedEnergy.toFixed(1) + " keV" );
      } else {
        // update the left sum peak line here
        self.leftSumPeak
          .attr("x1", leftSumPix)
          .attr("x2", leftSumPix)
          .attr("y1", 0)
          .attr("y2", self.size.height); 
        self.leftSumPeakText
          .attr("class", "peakText")
          .attr( "x", leftSumPix + xmax/125 )
          .attr( "y", self.size.height/3.4);
        self.leftSumPeakMeas
          .attr( "x", leftSumPix + xmax/125 )
          .attr( "y", self.size.height/3.0)
          .text( energy.toFixed(1) + "+" + leftSumEnergy.toFixed(1) + "=" + clickedEnergy.toFixed(1) + " keV" );
        }

        if (leftSumOutOfBounds)
          deleteLeftSumPeakMarker();

    } else
      deleteLeftSumPeakMarker();


    if ( shouldUpdateSumPeak ) {
      var sumEnergy = energy + clickedEnergy,
          sumPix = self.xScale( sumEnergy  );

      var sumPeakOutOfBounds = sumPix < 0 || sumPix > xmax;
      if ( sumPeakOutOfBounds ) {
        deleteSumPeakMarker();
        return;
      }

      if (!clickedEnergy)
        return;

      if( !self.sumPeak ) {
        // draw sum peak line here
        self.sumPeak = self.vis.append("line")
          .attr("class", "peakLine")
          .attr("stroke-width", 2)
          .attr("stroke", "black")
          .attr("x1", sumPix)
          .attr("x2", sumPix)
          .attr("y1", 0)
          .attr("y2", self.size.height);
      self.sumPeakText = self.vis.append("text")
          .attr("class", "peakText")
          .attr( "x", sumPix + xmax/200 )
          .attr( "y", self.size.height/4)
          .text( "Sum Peak" );
      self.sumPeakMeas = self.vis.append("text")
          .attr("class", "peakText")
          .attr( "x", sumPix + xmax/125 )
          .attr( "y", self.size.height/3.7)
          .text( clickedEnergy.toFixed(1) + "+" + energy.toFixed(1) + "=" + sumEnergy.toFixed(1) + " keV" );
      } else {
        if (!self.sumPeakMeas)
          self.sumPeakMeas = self.vis.append("text")
                .attr("class", "peakText")
                .attr( "x", sumPix + xmax/125 )
                .attr( "y", self.size.height/3.7)
                .text( clickedEnergy.toFixed(1) + "+" + energy.toFixed(1) + "=" + sumEnergy.toFixed(1) + " keV" );

        self.sumPeak
          .attr("x1", sumPix)
          .attr("x2", sumPix)
          .attr("y1", 0)
          .attr("y2", self.size.height);
        self.sumPeakText
          .attr( "x", sumPix + xmax/125 )
          .attr( "y", self.size.height/4);
        self.sumPeakMeas
          .attr( "x", sumPix + xmax/125 )
          .attr( "y", self.size.height/3.7)
          .text( clickedEnergy.toFixed(1) + "+" + energy.toFixed(1) + "=" + sumEnergy.toFixed(1) + " keV" );
      }
    }
  }

  updateEscapePeaks();
  updateComptonPeaks();
  updateComptonEdge();
  updateSumPeaks(sumPeaksArgument);
}

SpectrumChartD3.prototype.rebinSpectrum = function(spectrumToBeAdjusted, linei) {
  var self = this;

  if( !this.rawData || !this.rawData.spectra || !this.rawData.spectra.length ) 
    return;

  // Check for the which corresponding spectrum line is the specified one to be rebinned
  if (linei == null || !spectrumToBeAdjusted) 
    return;

  var newRebin = 1;

  var firstRaw = self.displayed_raw_start(spectrumToBeAdjusted);
  var lastRaw = self.displayed_raw_end(spectrumToBeAdjusted);

  var npoints = lastRaw - firstRaw;
  if( npoints > 1 && this.size.width > 2 )
    newRebin = Math.ceil( npoints / this.size.width );

  // Since we can only adjust scale factors for background/secondary/others, we don't need to adjust y-axis title
  /*
  if( spectrumToBeAdjusted.rebinFactor != newRebin ){
    var txt = this.options.ylabel ? this.options.ylabel : "";
    if( newRebin !== 1 )
      txt += " per " + newRebin + " Channels"
    d3.select(".yaxistitle").text( txt );
  }
  */

  spectrumToBeAdjusted.rebinFactor = newRebin;
  spectrumToBeAdjusted.firstRaw = firstRaw;
  spectrumToBeAdjusted.lastRaw = lastRaw;

  //Round firstRaw and lastRaw down and up to even multiples of newRebin
  firstRaw -= (firstRaw % newRebin);
  lastRaw += newRebin - (lastRaw % newRebin);
  if( firstRaw >= newRebin )
    firstRaw -= newRebin;
  if( lastRaw > spectrumToBeAdjusted.x.length )
    lastRaw = spectrumToBeAdjusted.x.length;


  //could do some optimizations here where we actually do a slightly larger
  //  range than displayed, so that next time we might not have to go back
  //  through the data to recompute things (it isnt clear if D3 will plot
  //  these datas, should check on this.)
  //Also, could _maybe_ use energy range, rather than indexes to track if we
  //  need to rebin the data or not...

  //console.log( "self.rebinFactor=" + self.rebinFactor );
  var i = firstRaw;
  for( var pointi = 0; pointi < spectrumToBeAdjusted.points.length; pointi++ ){
    var thisdata = spectrumToBeAdjusted.points[pointi];

    if( spectrumToBeAdjusted.y.length > 0 ) {
      thisdata.y = 0;
      for( var j = 0; j < newRebin && i+j<spectrumToBeAdjusted.y.length; ++j )
        thisdata.y += spectrumToBeAdjusted.y[i+j];
      thisdata.y *= spectrumToBeAdjusted.yScaleFactor;
    }

    i += newRebin;
  }
}

SpectrumChartD3.prototype.updateLegend = function() {
  var self = this;
  
  if( !this.options.showLegend ) {
    if( this.legend ) {
      this.legend.remove();
      this.legend = null;
      this.legendBox = null;
      this.legBody = null;
      this.legendHeaderClose = null;
      var legendoption;
      if (legendoption = document.getElementById("legendoption"))
        legendoption.checked = false;
      console.log( 'Should emit legend closed' );
    }
    return;
  }
  
  if( !this.legend ) {

    function moveleg(){                 // move legend 
      if( self.legdown ) {
        // console.log(d3.event);
        d3.event.preventDefault();
        d3.event.stopPropagation();

        var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX;
        var y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY;

        var calculated_x = d3.mouse(self.vis[0][0])[0];

        if ( calculated_x >= -self.padding.leftComputed && y >= 0 && 
             calculated_x <= self.cx && y <= self.cy ) {
          console.log("change pos");
          var tx = (x - self.legdown.x) + self.legdown.x0;
          var ty = (y - self.legdown.y) + self.legdown.y0; 
          self.legend.attr("transform", "translate(" + tx + "," + ty + ")");
        }
      }
    }


    this.legend = d3.select(this.chart).select("svg").append("g")
                      .attr("class", "legend")
                      .attr("transform","translate(" + (this.cx - 120 - this.padding.right) + ","+ (this.padding.topComputed + 10) + ")");
    this.legendBox = this.legend.append('rect')
               .attr("class", "legendBack")
               .attr('width', "100")
               .attr('height', "1em")
               .attr( "rx", "5px")
               .attr( "ry", "5px");
    this.legBody = this.legend.append("g")
                       .attr("transform","translate(8,17)");
                       
    this.legendHeader = this.legend.append("g"); 
    this.legendHeader
               .style("display", "none")
               .append('rect')
               .attr("class", "legendHeader")
               .attr('width', "100px")
               .attr('height', "1.5em")
               .attr( "rx", "5px")
               .attr( "ry", "5px")
               .style("cursor", "pointer");
    
    //Add a close button to get rid of the legend
    this.legendHeaderClose = this.legendHeader.append('g').attr("transform","translate(4,4)");
    this.legendHeaderClose.append("rect")
            .attr("class", "d3closebut")
            .attr('height', "12")
            .attr('width', "12");
    this.legendHeaderClose.append("path")
        .attr("style", "stroke: white; stroke-width: 1.5px;" )
        .attr("d", "M 2,2 L 10,10 M 10,2 L 2,10");
    this.legendHeaderClose.on("click", function(){ self.options.showLegend = false; self.updateLegend(); } )
                          .on("touchend", function(){ self.options.showLegend = false; self.updateLegend(); } );   
               
    //this.legendHeader.append("text")
    //      .text("Legend")
    //      .attr("x", 62.5)
    //      .attr("y", "1.1em")
    //      .style("text-anchor","middle")
    //      .style("cursor", "pointer");
    
    this.legend.on("mouseover", function(d){if( !self.dragging_plot && !self.zooming_plot ) self.legendHeader.style("display", null);} )
      .on("mouseout", function(d){self.legendHeader.style("display", "none");} )
      .on("mousemove", moveleg)
      .on("touchmove", moveleg)
      .on("wheel", function(d){d3.event.preventDefault(); d3.event.stopPropagation();} );
    
    function mousedownleg(){
      console.log("mouse down on legend");
      if (d3.event.defaultPrevented) return;
      if( self.dragging_plot || self.zooming_plot ) return;
      d3.event.preventDefault();
      d3.event.stopPropagation();
      var trans = d3.transform(self.legend.attr("transform")).translate;

      var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX;
      var y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY;
      self.legdown = {x: x, y: y, x0: trans[0], y0: trans[1]};
    };
    
    this.legendHeader.on("mouseover", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.legend.attr("class", "legend activeLegend");} )
      .on("touchstart", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.legend.attr("class", "legend activeLegend"); } )
      .on("mouseout",  function(d) { if (self.legend) self.legend.attr("class", "legend"); } )
      .on("touchend",  function(d) { if (self.legend) self.legend.attr("class", "legend"); } )
      .on("mousedown.drag",  mousedownleg )
      .on("touchstart.drag",  mousedownleg )
      .on("touchend.drag", function() {self.legdown = null;})
      .on("mouseup.drag", function(){self.legdown = null;} )
      .on("mousemove.drag", moveleg)
      .on("mouseout.drag", moveleg);

    this.legend.on("touchstart", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.legendHeader.style("display", null); } )
      .on("touchstart.drag", mousedownleg)
      .on("touchend.drag",  function() 
      {
      if (!self.legend || !self.legendHeader) {
        self.legdown = null;
        return;
      }

      self.legend.attr("class", "legend"); 
      window.setTimeout(function() { self.legendHeader.style("display", "none"); }, 1500) 
      self.legdown = null; 
      });
  }
  
  var origtrans = d3.transform(this.legend.attr("transform")).translate;
  var bb = this.legend.node().getBBox();
  var fromRight = this.cx - origtrans[0] - this.legendBox.attr('width');
  
  this.legBody.selectAll("g").remove();

  var ypos = 0;
  this.rawData.spectra.forEach( function(spectrum,i){
    if( !spectrum || !spectrum.y.length )
      return;
    
    var sf = spectrum.yScaleFactor;
    var title = spectrum.title;
    var lt = spectrum.liveTime;
    var rt = spectrum.realTime;
    var neut = spectrum.neutrons;
    title = title ? title : ("Spectrum " + (i+1));
    var nsum = spectrum.dataSum;
    title += " (" + nsum.toFixed(nsum > 1000 ? 0 : 1) + " counts)"; 
   
    var thisentry = self.legBody.append("g")
      .attr("transform","translate(0," + ypos + ")");
    
    thisentry.append("path")
      .attr("class", "line" )
      .attr("stroke", spectrum.lineColor ? spectrum.lineColor : "black")
      .attr("d", "M0,-5 L12,-5");
    var thistxt = thisentry.append("text")
      .attr("class", "legentry")
      .attr( "x", 15 )
      .text(title);
    if( typeof lt === "number" )
      thistxt.append('svg:tspan').attr('x', "20").attr('y', thisentry.node().getBBox().height)
                    .text( "Live Time " + lt.toPrecision(4) );
    if( typeof rt === "number" )
      thistxt.append('svg:tspan').attr('x', "20").attr('y', thisentry.node().getBBox().height)
                    .text( "Real Time " + rt.toPrecision(4) );
    if( typeof neut === "number" )
      thistxt.append('svg:tspan').attr('x', "20").attr('y', thisentry.node().getBBox().height)
                    .text( "Neutron Count " + neut );
    if( typeof sf === "number" && sf != 1 )
      thistxt.append('svg:tspan').attr('x', "20").attr('y', thisentry.node().getBBox().height)
                    .text( "Scaled by " + sf.toPrecision(4) );
                    
    ypos += thisentry.node().getBBox().height + 5;
  });
                    
  //Resize the box to match the text size
  var w = this.legBody.node().getBBox().width + 15; 
  this.legendBox.attr('width', w );
  this.legendBox.attr('height', this.legBody.node().getBBox().height + 10 );
  this.legendHeaderClose.attr("transform","translate(" + (w-16) + ",4)");
  
  this.legendHeader.select('rect').attr('width', w );
  this.legendHeader.select('text').attr("x", w/2);
  
  //this.legendBox.attr('height', legtxt.node().getBBox().height + hh + 8 );

  //Set the transform so the space on the right of the legend stays the same
  this.legend.attr("transform", "translate(" + (this.cx - fromRight - w) + "," + origtrans[1] + ")" );
}


SpectrumChartD3.prototype.drawRefGammaLines = function() {
  //Drawing of the refrenece lines is super duper un-optimized!!!
  var self = this;

  if( !self.refLines.length || !self.refLines[0].lines  || !self.refLines[0].lines.length ) {
    //console.log( "no reference lines, " + self.refLines.length );
    return;
  }

  function getLinesInRange(xrange,lines) {
    var bisector = d3.bisector(function(d){return d.e;});
    var lindex = bisector.left( lines, xrange[0] );
    var rindex = bisector.right( lines, xrange[1] );
    return lines.slice(lindex,rindex).filter(function(d){return d.h > 1E-16;});
  }

  var lowerx = this.xScale.domain()[0], upperx = this.xScale.domain()[1];


  var reflines = [];
  self.refLines.forEach( function(input) {
    var lines = getLinesInRange(self.xScale.domain(),input.lines);
    input.maxVisibleAmp = d3.max(lines, function(d){return d.h;});  //same as lines[0].parent.maxVisibleAmp = ...
    reflines = reflines.concat( lines );
  });

  reflines.sort( function(l,r){ return ((l.e < r.e) ? -1 : (l.e===r.e ? 0 : 1)); } );

  var tx = function(d) { return "translate(" + self.xScale(d.e) + ",0)"; };
  var gy = self.vis.selectAll("g.ref")
            .data( reflines,function(d){return d.e;} )
            .attr("transform", tx)
            .attr("stroke-width",1);

  var gye = gy.enter().insert("g", "a")
    .attr("class", "ref")
    .attr("transform", tx);

  //color:'#0000FF',parent:'Ba133',age:'0.00 us',lines:[{e:30.27,h:6.22e-05,pa
  var stroke = function(d) { return d.parent.color; };

  var dashfunc = function(d){
    var particles = ["gamma", "xray", "beta", "alpha",   "positron", "electronCapture"];
    var dash      = [null,    ("3,3"),("1,1"),("3,2,1"), ("3,1"),    ("6,6") ];
    var index = particles.indexOf(d.particle);
    if( index < 0 ) { console.log( 'Invalid particle: ' + d.particle ); return null; }
    return (index > -1) ? dash[index] : null;
  };

  var h = self.size.height;
  var m = Math.min(h,self.options.refLineTopPad); //leave 20px margin at top of chart


  gye.append("line")
    .style("stroke-dasharray", dashfunc )
    .attr("stroke", stroke )
    // .attr("class", function(d,i){ return "refline " + i;} )
    .attr("y1", h )
    .attr("dx", "-0.5" )
    // .attr("y2", function(d){ console.log("y2="); return h - (h-m)*d.h/maxVisibleAmp; } )
    //.on("mouseover", function(d) { d3.select(this).style("stroke-width", "3");})
    //.on("mouseout",  function(d) { d3.select(this).style("stroke-width", "1");})
     ;

  // EXIT
  // Remove old elements as needed.
  gy.exit().remove();

  //Now update the height of all the lines.  If we did this in the gye.append("line")
  //  line above then the values for existing lines wouldnt be updated (only
  //  the new lines would have correct height)
  gy.select("line").attr("y2", function(d){ return Math.min(h - (h-m)*d.h/d.parent.maxVisibleAmp,h-2) ; } );
}

SpectrumChartD3.prototype.calcLeftPadding = function( updategeom ){
  
  if( !this.options.adjustYAxisPadding ) {
    this.padding.leftComputed = this.padding.left; 
    return;
  }
  
  var labels = this.vis.selectAll("g.y").selectAll("text");
  
  var labelw = 4;
  labels.forEach( function(label){
    labelw = Math.max( labelw, label.parentNode.getBBox().width );
  });
  
  var labelpad = this.padding.labelPad;
  var ticklen = 7; //hardcoded as major tick x1
  var title = d3.select(".yaxistitle");
  if( title )
    title.attr("transform","translate(-" + (labelw+labelpad) + " " + this.size.height/2+") rotate(-90)");
  
  var newleft = ticklen + labelw + labelpad + this.padding.left + 4;
  
  if( !updategeom ) {
    this.padding.leftComputed = newleft;
  } else if( Math.abs(this.padding.leftComputed - newleft) > 1.0 ) {
    this.padding.leftComputed = newleft;
    this.handleResize( true );
  } 
}


SpectrumChartD3.prototype.drawXTicks = function() {
  var self = this;
  var stroke = function(d) { return d ? "#ccc" : "#666"; };
  
  var xticks = self.xticks();
  var xtickvalues = xticks.map(function(d){return d.value;} );
  self.xAxis.tickValues( xtickvalues );

  // Regenerate x-ticks…
  self.xAxisBody.selectAll("g.tick").remove();

  self.xAxisBody.call(self.xAxis);
  
  self.xAxisBody.selectAll("text").style("cursor", "ew-resize")
    .on("mouseover", function(d) { d3.select(this).style("font-weight", "bold");})
    .on("mouseout",  function(d) { d3.select(this).style("font-weight", "normal");})
    .on("mousedown.drag",  self.xaxisDrag())
    .on("touchstart.drag", self.xaxisDrag());

  var xgticks = self.xAxisBody.selectAll('g.tick');
  var minorticks = xgticks.filter(function(d,i){ return xticks[i] && !xticks[i].major; } );
  minorticks.select('line').attr('y2', '4');
  minorticks.select('text').text("");

  //Check that the last tick doesnt go off the chart area.
  //  This could probably be accomplished MUCH more efficiently
  var majorticks = xgticks.filter(function(d,i){ return xticks[i] && xticks[i].major; } );
  if( this.options.compactXAxis ){
    //We have to check every tick to see if it overlaps with the title
    var xtitlex = this.xaxistitle.attr("x" );
    majorticks[0].forEach( function(tick){
      var txt = d3.select(tick).select('text')[0][0]; 
      if( (txt.getCTM().e + txt.getBBox().width) + 30 > xtitlex )  //Not sure why we need this 30, but its in positioning of the x-axis title too
        d3.select(txt)
        .text("");
    });
  } else {
    //We only need to check the last tick to see if it goes off the chart
    var lastmajor = majorticks[0].length ? majorticks[0][majorticks[0].length-1] : null; 
    if( lastmajor ) {
      lastmajor = d3.select(lastmajor).select('text')[0][0];
      if( (lastmajor.getCTM().e + lastmajor.getBBox().width) > this.cx )
        d3.select(lastmajor).text("");  
    }
  }
  //this.options.compactXAxis
  
  //Can calculate the width needed for the y-axis, and could adjust the plot area..
  //console.log( self.yAxisBody.node().getBBox().width );

  if( self.xGridBody ) {
    self.xGrid.tickValues( xtickvalues );
    self.xGridBody.remove();
    self.xGridBody = self.vis.insert("g", ".peakVis")
          .attr("width", self.size.width )
          .attr("height", self.size.height )
          .attr("class", "xgrid" )
          .attr("transform", "translate(0," + self.size.height + ")")
          .call( self.xGrid );
    self.xGridBody.selectAll('g.tick')
      .filter(function(d,i){ return !xticks[i].major; } )
      .attr("class","minorgrid");
  }
}


SpectrumChartD3.prototype.drawYTicks = function() {
  var self = this;
  var stroke = function(d) { return d ? "#ccc" : "#666"; };
  
  // Regenerate y-ticks…
    var ytick = self.yticks();
    var ytickvalues = ytick.map(function(d){return d.value;} );
    self.yScale.ticks(ytickvalues);

    if( self.yGridBody ) {
      self.yGrid.tickValues( ytickvalues );

      //Since the number of grid lines might change (but call(self.yGrid) expects the same number)
      //  we will remove, and re-add back in the grid... kinda a hack.
      //  Could probably go back to manually drawing the grid lined and get rid
      //  of yGrid and yGridBody...
      self.yGridBody.remove();
      self.yGridBody = self.vis.insert("g", ".peakVis")
        .attr("width", self.size.width )
        .attr("height", self.size.height )
        .attr("class", "ygrid" )
        .attr("transform", "translate(0,0)")
        .call( self.yGrid );
      self.yGridBody.selectAll('g.tick')
        .filter(function(d,i){return !ytick[i].major;} )
        .attr("class","minorgrid");
    }

    var ty = function(d) { return "translate(0," + self.yScale(d) + ")"; };
    var gy;

    gy = self.vis.selectAll("g.y")
      .data(ytickvalues, String)
      .attr("transform", ty);


    var fy = function(d,i){ return ytick[i].text; }; //self.yScale.tickFormat(10);
    gy.select("text").text( fy );

    var gye = gy.enter().insert("g", "a")
        .attr("class", "y")
        .attr("transform", ty)
        .attr("background-fill", "#FFEEB6");

    gye.append("line")
       .attr("stroke", stroke)
       .attr("class", "yaxistick")
       .attr("x1", function(d,i){ return ytick[i].major ? -7 : -4; } )
       .attr("x2", 0 );

    gye.append("text")
        .attr("class", "yaxislabel")
        .attr("x", -8.5)
        .attr("dy", ".35em")
        .attr("text-anchor", "end")
        .text(fy);
        // .style("cursor", "ns-resize")
        // .on("mouseover", function(d) { d3.select(this).style("font-weight", "bold");})
        // .on("mouseout",  function(d) { d3.select(this).style("font-weight", "normal");});
        // .on("mousedown.drag",  self.yaxisDrag())
        // .on("touchstart.drag", self.yaxisDrag());

    gy.exit().remove();
}

SpectrumChartD3.prototype.redraw = function() {
  var self = this;

  return function() {
    self.do_rebin();
    self.setYAxisDomain();

    var tx = function(d) { return "translate(" + self.xScale(d) + ",0)"; };
    var stroke = function(d) { return d ? "#ccc" : "#666"; };

    self.yAxisBody.call(self.yAxis);
    /*self.yAxisBody.selectAll("text")*/ d3.selectAll('.yaxislabel').style("cursor", "ns-resize")
        .on("mouseover", function(d) { d3.select(this).style("font-weight", "bold");})
        .on("mouseout",  function(d) { d3.select(this).style("font-weight", "normal");})
        .on("mousedown.drag",  self.yaxisDrag())
        .on("touchstart.drag", self.yaxisDrag());

    if (self.options.showXAxisSliderChart) { self.drawXAxisSliderChart(); } 
    else                                   { self.cancelXAxisSliderChart(); }

    self.drawYTicks();
    
    self.calcLeftPadding( true );
    
    self.drawXTicks();

    self.drawXAxisArrows();

    //The oringal code had the follwoing in, but I dont think is needed
    //self.plot.call(d3.behavior.zoom().x(self.xScale).y(self.yScale).on("zoom", self.redraw()));    
    
    self.drawPeaks();
    self.drawRefGammaLines();
    self.updateMouseCoordText();

    self.update();

    self.yAxisZoomedOutFully = true;
  }
}

SpectrumChartD3.prototype.drawXAxisArrows = function(show_arrow) {
  var self = this;

  if (self.options.showXRangeArrows) {
    var max_x;

    if (!self.xaxis_arrow) {
      self.xaxis_arrow = d3.select('.xaxis').append('svg:defs')
        .attr("id", "xaxisarrowdef")
        .append("svg:marker")
        .attr("id", "arrowhead")
        .attr('class', 'xaxisarrow')
        .attr("refX", 0)
        .attr("refY", 6)
        .attr("markerWidth", 9)
        .attr("markerHeight", 14)
        .attr("orient", 0)
        .append("path")
          .attr("d", "M2,2 L2,13 L8,7 L2,2")
          .style("stroke", "black");
    }

    max_x = self.min_max_x_values()[1];

    if (self.xScale.domain()[1] === max_x || (typeof show_arrow == 'boolean' && !show_arrow))            // should be a better way to determine if can still pan
      self.xAxisBody.select("path").attr("marker-end", null);
    else
      self.xAxisBody.select("path").attr("marker-end", "url(#arrowhead)");

  } else {
    if (self.xaxis_arrow) {
      self.xaxis_arrow.remove();
      self.xaxis_arrow = null;
    }

    d3.select("#xaxisarrowdef").remove();
  }
}

/* SpectrumChartD3 - Scale Factor Slider Methods */
SpectrumChartD3.prototype.drawScalerBackgroundSecondary = function() {
  var self = this;

  function hasYScaleFactors() {
    var result = false;
    self.rawData.spectra.forEach(function (spectrum) { 
      if (spectrum.yScaleFactor != null)
        result = true;
    });
    return result;
  }
  function allNegativeYScaleFactors() {
    var result = true;
    self.rawData.spectra.forEach(function (spectrum) { 
      if (spectrum.yScaleFactor && spectrum.yScaleFactor >= 0)
        result = false;
    });
    return result;
  }
  function onlyForegroundPresent() {
    var result = true;
    self.rawData.spectra.forEach(function (spectrum) { 
      if (spectrum && spectrum.title && spectrum.title.toUpperCase() != "FOREGROUND")
        result = false;
    });
    return result;
  }

  var onlyForegroundPresentInGraph = onlyForegroundPresent();
  var graphHasYScaleFactors = hasYScaleFactors();

  if (!self.rawData || !self.rawData.spectra || !graphHasYScaleFactors || onlyForegroundPresentInGraph || !self.options.scaleBackgroundSecondary) {
    if (self.scalerWidget) {
      self.scalerWidget.remove();
      self.scalerWidget = null;
      self.scalerWidgetBox = null;
      self.scalerWidgetHeader = null;
      self.scalerWidgetBody = null;
      self.scalerWidgetClose = null;
      self.scalerWidgetTitle = null;
      self.backgroundScale = null;
      var scaleroption;
      if (scaleroption = document.getElementById("scaleroption"))
        scaleroption.checked = false;
      console.log( 'Should emit scaler widget closed' );

      // Display proper error messages
      if (!self.rawData || !self.rawData.spectra)
        alert("No data specified for graph!");
      else if (!graphHasYScaleFactors)
        alert("No y-scale factors detected for any graph!");
      else if (onlyForegroundPresentInGraph)
        alert("Only chat foreground detected. Please specify a background/other chart for scaling y-values.");
    }
    return;
  }

  if (allNegativeYScaleFactors()) {
    if (self.scalerWidget) {
      self.scalerWidget.remove();
      self.scalerWidget = null;
      self.scalerWidgetBox = null;
      self.scalerWidgetHeader = null;
      self.scalerWidgetBody = null;
      self.scalerWidgetClose = null;
      self.scalerWidgetTitle = null;
      self.backgroundScale = null;
      console.log( 'Should emit scaler widget closed' );
    }
    return;
  }

  if (!self.scalerWidget) {

    function moveScalerWidget(){
      if( self.scalerdown && !self.adjustingBackgroundScale && !self.adjustingSecondaryScale ) {
        d3.event.preventDefault();
        d3.event.stopPropagation();

        var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX;
        var y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY;

        var calculated_x = d3.mouse(self.vis[0][0])[0];

        if ( calculated_x >= -self.padding.leftComputed && y >= 0 && 
             calculated_x <= self.cx && y <= self.cy ) {
          console.log("change pos");
          var tx = (x - self.scalerdown.x) + self.scalerdown.x0;
          var ty = (y - self.scalerdown.y) + self.scalerdown.y0; 
          self.scalerWidget.attr("transform", "translate(" + tx + "," + ty + ")");
        }
      }
    }

    self.scalerWidget = d3.select(this.chart).select("svg").append("g")
      .attr("class", "scalerwidget")
      .attr("transform","translate(" + (this.cx - (this.cx/2) - this.padding.right) + ","+ (this.padding.topComputed + 10) + ")");

    self.scalerWidgetBox = self.scalerWidget.append('rect')
      .attr("class", "scalerback")
      .attr('width', "100px")
      .attr('height', "5px")
      .attr( "rx", "5px")
      .attr( "ry", "5px");
    self.scalerWidgetBody = self.scalerWidget.append("g")
      .attr("transform","translate(8,17)");

    self.scalerWidgetHeader = self.scalerWidget.append("g"); 
    self.scalerWidgetHeader.style("display", "none")
      .append('rect')
      .attr("class", "scalerheader")
      .attr('width', "100px")
      .attr('height', "1.5em")
      .attr( "rx", "5px")
      .attr( "ry", "5px")
      .style("cursor", "pointer");
    
    //Add a close button to get rid of the scaler widget
    self.scalerWidgetClose = self.scalerWidgetHeader.append('g').attr("transform","translate(4,4)");
    self.scalerWidgetClose.append("rect")
      .attr("class", "d3closebut")
      .attr('height', "12")
      .attr('width', "12");
    self.scalerWidgetClose.append("path")
      .attr("style", "stroke: white; stroke-width: 1.5px;" )
      .attr("d", "M 2,2 L 10,10 M 10,2 L 2,10");
    self.scalerWidgetClose.on("click", function(){ self.options.scaleBackgroundSecondary = false; self.drawScalerBackgroundSecondary(); } )
                          .on("touchend", function(){ self.options.scaleBackgroundSecondary = false; self.drawScalerBackgroundSecondary(); } );   
    
    self.scalerWidget.on("mouseover", function(d){if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidgetHeader.style("display", null);} )
      .on("mouseout", function(d){self.scalerWidgetHeader.style("display", "none");} )
      .on("mousemove", moveScalerWidget)
      .on("touchmove", moveScalerWidget)
      .on("wheel", function(d){d3.event.preventDefault(); d3.event.stopPropagation();} );
    
    function mousedownscaler(){
      console.log("mouse down on scaler");
      if (d3.event.defaultPrevented) return;
      if( self.dragging_plot || self.zooming_plot ) return;
      d3.event.preventDefault();
      d3.event.stopPropagation();
      var trans = d3.transform(self.scalerWidget.attr("transform")).translate;

      var x = d3.event.x ? d3.event.x : d3.event.touches ?  d3.event.touches[0].clientX : d3.event.clientX;
      var y = d3.event.y ? d3.event.y : d3.event.touches ?  d3.event.touches[0].clientY : d3.event.clientY;
      self.scalerdown = {x: x, y: y, x0: trans[0], y0: trans[1]};
    }
    
    self.scalerWidgetHeader.on("mouseover", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidget.attr("class", "scalerwidget activescaler");} )
      .on("touchstart", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidget.attr("class", "scalerwidget activescaler"); } )
      .on("mouseout",  function(d) { self.scalerWidget.attr("class", "scalerwidget"); } )
      .on("touchend",  function(d) { if (self.scalerWidget) self.scalerWidget.attr("class", "scalerwidget"); } )
      .on("mousedown.drag",  mousedownscaler )
      .on("touchstart.drag",  mousedownscaler )
      .on("touchend.drag", function() {self.scalerdown = null;})
      .on("mouseup.drag", function(){self.scalerdown = null;} )
      .on("mousemove.drag", moveScalerWidget)
      .on("mouseout.drag", moveScalerWidget);

    self.scalerWidget.on("touchstart", function(d) { mousedownscaler(); if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidgetHeader.style("display", null); } )
      .on("touchend.drag",  function() 
      {
      if (!self.scalerWidget || !self.scalerWidgetHeader) {
        self.scalerdown = null;
        return; 
      }
      self.scalerWidget.attr("class", "scalerwidget"); 
      window.setTimeout(function() { self.scalerWidgetHeader.style("display", "none"); }, 1500) 
      self.scalerdown = null; 
      });
  }

  var origtrans = d3.transform(self.scalerWidget.attr("transform")).translate;
  var bb = self.scalerWidget.node().getBBox();
  var fromRight = self.cx - origtrans[0] - self.scalerWidgetBox.attr('width').replace("px","");

  self.scalerWidgetBody.selectAll("g").remove();

  function is_touch_device() {
    return 'ontouchstart' in window        // works on most browsers 
        || (navigator && navigator.maxTouchPoints);       // works on IE10/11 and Surface
  };

  var scalerWidth = 300;
  var toggleRadius = is_touch_device() ? 10 : 7;
  var xpos = 10;
  var ypos = 25;
  var scalerAdded = false;

  self.rawData.spectra.forEach(function(spectrum,i) {
    var spectrumScaleFactor = spectrum.yScaleFactor;
    var spectrumTitle = spectrum.title ? spectrum.title : "Spectrum " + i;

    if (i == 0 || spectrum.title == "Foreground")   // Don't add scaling functionality for foreground
      return;

    if (spectrumScaleFactor != null && spectrumScaleFactor >= 0) {
      spectrum.scale = d3.scale.linear()
        .domain([0, self.options.maxScaleFactor]) // TODO: Have global max scale factor
        .range([0, scalerWidth]);

      spectrum.scaleAxis = d3.svg.axis().scale(spectrum.scale)
        .orient("bottom")
        .innerTickSize(is_touch_device() ? 27 : 19)
        .outerTickSize(0)
        .ticks(5,'f');

      var spectrumSliderArea = self.scalerWidgetBody.append("g")
        .attr("id", spectrumTitle.replace(' ','-') + "SliderArea")
        .attr("transform","translate(" + xpos + "," + ypos + ")")
        .call(spectrum.scaleAxis);

      spectrum.sliderText = spectrumSliderArea.append("text")
        .attr("x", 0)
        .attr("y", 0)
        .attr("text-anchor", "start")
        .text(spectrumTitle + " Scale Factor: " + spectrumScaleFactor.toFixed(3));

      spectrum.sliderRect = spectrumSliderArea.append("rect")
        .attr("class", "scaleraxis")
        .attr("y", Number(spectrum.sliderText.attr("y")) + 10 + (is_touch_device() ? 5 : 0))
        .attr("rx", 5)
        .attr("ry", 5)
        .attr("width", scalerWidth)
        .attr("height", "1%");

      spectrum.sliderToggle = spectrumSliderArea.append("circle")
        .attr("class", "scalertoggle")
        .attr("cx", Math.min(spectrum.scale(spectrumScaleFactor), scalerWidth))
        .attr("cy", Number(spectrum.sliderRect.attr("y")) + toggleRadius/2)
        .attr("r", toggleRadius)
        .style("cursor", "pointer")
        .on("mousedown", function(){ self.currentlyAdjustingSpectrumScale = spectrum.title; })
        .on("mousemove", self.handleMouseMoveScaleFactorSlider())
        .on("mouseup", function(){ self.currentlyAdjustingSpectrumScale = null; })
        .on("touchstart", function(){ self.currentlyAdjustingSpectrumScale = spectrum.title; })
        .on("touchmove", self.handleMouseMoveScaleFactorSlider())
        .on("touchend", function(){ self.currentlyAdjustingSpectrumScale = null; });

      ypos += 60;
    }
  });

  //Resize the box to match the text size
  var w = self.scalerWidgetBody.node().getBBox().width + 15 + (2*xpos); 
  self.scalerWidgetHeader.attr("width", w)
  self.scalerWidgetBox.attr('width', w);
  self.scalerWidgetBox.attr('height', self.scalerWidgetBody.node().getBBox().height + 40 );
  self.scalerWidgetClose.attr("transform","translate(" + (w-16) + ",4)");

  if (!self.scalerWidgetTitle)
    self.scalerWidgetTitle = self.scalerWidget.append('text')
      .attr("x", w/2)
      .attr("text-anchor", "middle")
      .attr("transform", "translate(0," + 16 + ")")
      .attr("cursor", "pointer")
      .on("touchstart", function(d) { mousedownscaler();  if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidget.attr("class", "scalerwidget activescaler"); } )
      .on("mousemove", function(d) { if( !self.dragging_plot && !self.zooming_plot ) self.scalerWidget.attr("class", "scalerwidget activescaler"); } )
      .on("mouseout",  function(d) { self.scalerWidget.attr("class", "scalerwidget"); } )
      .on("touchend",  function(d) { self.scalerWidget.attr("class", "scalerwidget"); } )
      .on("mousedown.drag",  mousedownscaler )
      .on("touchstart.drag",  mousedownscaler )
      .on("touchend.drag", function() {self.scalerdown = null;})
      .on("mouseup.drag", function(){self.scalerdown = null;} )
      .on("mousemove.drag", moveScalerWidget)
      .text("Spectrum Y-Value Scaler")
  
  self.scalerWidgetHeader.select('rect').attr('width', w );
  self.scalerWidgetHeader.select('text').attr("x", w/2);

  //Set the transform so the space on the right of the legend stays the same
  self.scalerWidget.attr("transform", "translate(" + (this.cx - fromRight - w) + "," + origtrans[1] + ")" );
}

SpectrumChartD3.prototype.handleMouseMoveScaleFactorSlider = function() {
  var self = this;

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  /* Here is a modified, slightly more efficient version of redraw 
    specific to changing the scale factors for spectrums. 
  */
  function scaleFactorChangeRedraw(spectrum, linei) {
    self.updateLegend();
    self.rebinSpectrum(spectrum, linei);
    self.do_rebin();
    self.setYAxisDomain();
    self.drawYTicks();
    self.update();
    self.drawPeaks();
  }

  return function() {
    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length )
      return;
    if (!self.currentlyAdjustingSpectrumScale)
      return;

    // Check for the which corresponding spectrum line is the background
    var linei = null;
    var spectrum = null;
    for (var i = 0; i < self.rawData.spectra.length; ++i) {
      if (self.rawData.spectra[i].title.toUpperCase() == self.currentlyAdjustingSpectrumScale.toUpperCase()) {
        linei = i;
        spectrum = self.rawData.spectra[i];
        break;
      }
    }

    d3.event.preventDefault();
    d3.event.stopPropagation();

    if (linei == null || spectrum == null)
      return;

    d3.select(document.body).style("cursor", "pointer");

    var m = d3.mouse(spectrum.sliderRect[0][0]);
    var minrange = spectrum.scale.range()[0], maxrange = spectrum.scale.range()[1];
    var min = spectrum.scale.domain()[0];
    var max = spectrum.scale.domain()[1];

    var spectrumScaleFactor = Math.min(Math.max(min, spectrum.scale.invert(m[0])), max);
    var spectrumTitle = spectrum.title ? spectrum.title : "Spectrum " + i;

    spectrum.sliderToggle.attr("cx", spectrum.scale(spectrumScaleFactor));
    spectrum.sliderText.text(spectrumTitle + " Scale Factor: " + (needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed()));
    self.rawData.spectra[linei].yScaleFactor = spectrumScaleFactor;

    scaleFactorChangeRedraw(spectrum, linei);

    var currentsfinput;
    if (self.currentDropDownScaleFactorSpectrum == spectrumTitle && (currentsfinput = document.getElementById("current-sf")))
      currentsfinput.value = needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed();

    // Update the slider chart if needed
    if (self["sliderLine"+linei]) {
      var origdomain = self.xScale.domain();
      var origdomainrange = self.xScale.range();
      var origrange = self.yScale.domain();
      var bounds = self.min_max_x_values();
      var maxX = bounds[1];
      var minX = bounds[0];

      // Change the x and y-axis domain to the full range (for slider lines)
      self.xScale.domain([minX, maxX]);
      self.xScale.range([0, self.size.sliderChartWidth]);
      self.do_rebin();
      self.yScale.domain(self.getYAxisDomain());

      self.rawData.spectra.forEach(function(spec, speci) {
        if (self["sliderLine"+speci])
          self["sliderLine"+speci].attr("d", self["line"+speci](spec.points));
      })
      

      // Restore the original x and y-axis domain
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      self.do_rebin();
      self.yScale.domain(origrange);
    }
  }
}

SpectrumChartD3.prototype.handleTouchMoveAdjustSpectrumScale = function() {
  var self = this;

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  /* Here is a modified, slightly more efficient version of redraw 
    specific to changing the scale factors for spectrums. 
  */
  function scaleFactorChangeRedraw(spectrum, linei) {
    self.updateLegend();
    self.rebinSpectrum(spectrum, linei);
    self.do_rebin();
    self.setYAxisDomain();
    self.drawYTicks();
    self.update();
    self.drawPeaks();
  }

  return function() {
    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length )
      return;
    if (!self.currentlyAdjustingSpectrumScale)
      return;


    // Check for the which corresponding spectrum line is the background
    var linei = null;
    var spectrum = null;
    for (var i = 0; i < self.rawData.spectra.length; ++i) {
      if (self.rawData.spectra[i].title.toUpperCase() == self.currentlyAdjustingSpectrumScale.toUpperCase()) {
        linei = i;
        spectrum = self.rawData.spectra[i];
        break;
      }
    }

    if (linei == null || spectrum == null)
      return;

    d3.event.preventDefault();
    d3.event.stopPropagation();

    var t = d3.touches(spectrum.sliderRect[0][0]);;

    if (t.length !== 1)
      return;
    t = t[0];

    var minrange = spectrum.scale.range()[0], maxrange = spectrum.scale.range()[1];
    var min = spectrum.scale.domain()[0];
    var max = spectrum.scale.domain()[1];

    var spectrumScaleFactor = Math.min(Math.max(min, spectrum.scale.invert(m[0])), max);
    var spectrumTitle = spectrum.title ? spectrum.title : "Spectrum " + i;

    spectrum.sliderToggle.attr("cx", spectrum.scale(spectrumScaleFactor));
    spectrum.sliderText.text(spectrumTitle + " Scale Factor: " + (needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed()));
    self.rawData.spectra[linei].yScaleFactor = spectrumScaleFactor;

    scaleFactorChangeRedraw(spectrum, linei);

    var currentsfinput;
    if (self.currentDropDownScaleFactorSpectrum == spectrumTitle && (currentsfinput = document.getElementById("current-sf")))
      currentsfinput.value = needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed();

    // Update the slider chart if needed
    if (self["sliderLine"+linei]) {
      var origdomain = self.xScale.domain();
      var origdomainrange = self.xScale.range();
      var origrange = self.yScale.domain();
      var bounds = self.min_max_x_values();
      var maxX = bounds[1];
      var minX = bounds[0];

      // Change the x and y-axis domain to the full range (for slider lines)
      self.xScale.domain([minX, maxX]);
      self.xScale.range([0, self.size.sliderChartWidth]);
      self.do_rebin();
      self.yScale.domain(self.getYAxisDomain());

      self.rawData.spectra.forEach(function(spec, speci) {
        if (self["sliderLine"+speci])
          self["sliderLine"+speci].attr("d", self["line"+speci](spec.points));
      })

      // Restore the original x and y-axis domain
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      self.do_rebin();
      self.yScale.domain(origrange);
    }
  }
}


/* SpectrumChartD3 - X-axis Slider Methods */
SpectrumChartD3.prototype.drawXAxisSliderChart = function() {
  var self = this;

  // Cancel if the chart or raw data are not present
  if (!self.chart || d3.select(self.chart).empty() || !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length) {
    return;
  }
  // Cancel the action and clean up if the option for the slider chart is not checked
  if (!self.options.showXAxisSliderChart && self.origBBox && self.origHeight) {
    self.cancelXAxisSliderChart();
    return;
  }

  function drawSliderChartLines() {
    // Delete the data lines if they are present
    for (var i = 0; i < self.rawData.spectra.length; ++i)
      if (self['sliderLine' + i])
        self['sliderLine' + i].remove();

    for (var i = 0; i < self.rawData.spectra.length; ++i) {
      var spectrum = self.rawData.spectra[i];
      if (self['line'+i])
        self['sliderLine'+i] = self.sliderChartBody.append("path")
          .attr("id", 'sliderLine'+i)
          .attr("class", 'sline')
          .attr("stroke", spectrum.lineColor ? spectrum.lineColor : 'black')
          .attr("d", self['line'+i](spectrum.points))
          .attr("transform","scale(1," + (self.size.sliderChartHeight / self.size.height) + ")");
    }
  }

  function drawDragRegionLines() {
    d3.selectAll(".sliderDragRegionLine").remove();
    
    if (!self.sliderDragLeft || !self.sliderDragRight) {
      return;
    }
    var leftX = Number(self.sliderDragLeft.attr("x"));
    var leftY = Number(self.sliderDragLeft.attr("y"));
    var leftWidth =  Number(self.sliderDragLeft.attr("width"));
    var leftHeight = self.sliderDragLeft[0][0].height.baseVal.value;
    var rightX = Number(self.sliderDragRight.attr("x"));
    var rightY = Number(self.sliderDragRight.attr("y"));
    var rightWidth =  Number(self.sliderDragRight.attr("width"));
    var rightHeight = self.sliderDragRight[0][0].height.baseVal.value;

    var numberOfLines = 4;

    for (var i = 1; i < numberOfLines; i++) {
      self.sliderChart.append('line')
        .attr("class", "sliderDragRegionLine")
        .attr("x1", leftX + (i*leftWidth)/numberOfLines)
        .attr("x2", leftX + (i*leftWidth)/numberOfLines)
        .attr("y1", leftY + (leftHeight/4))
        .attr("y2", leftY + (3*leftHeight/4))
        .on("mousedown", self.handleMouseDownLeftSliderDrag())
        .on("mousemove", self.handleMouseMoveLeftSliderDrag(false))
        .on("mouseout", self.handleMouseOutLeftSliderDrag())
        .on("touchstart", self.handleTouchStartLeftSliderDrag())
        .on("touchmove", self.handleTouchMoveLeftSliderDrag(false));

      self.sliderChart.append('line')
        .attr("class", "sliderDragRegionLine")
        .attr("x1", rightX + (i*rightWidth)/numberOfLines)
        .attr("x2", rightX + (i*rightWidth)/numberOfLines)
        .attr("y1", rightY + (rightHeight/4))
        .attr("y2", rightY + (3*rightHeight/4))
        .on("mousedown", self.handleMouseDownRightSliderDrag())
        .on("mousemove", self.handleMouseMoveRightSliderDrag(false))
        .on("mouseout", self.handleMouseOutRightSliderDrag())
        .on("touchstart", self.handleTouchStartRightSliderDrag())
        .on("touchmove", self.handleTouchMoveRightSliderDrag(false));
    }
  }

  // Ensure that the slider chart height is always 1/10 of the chart height (for resizes)
  self.size.sliderChartHeight = self.size.height / 10;

  // Starting height of the chart (in px)
  var startingHeight = Number(d3.select(this.chart)[0][0].style.height.substring(0, d3.select(this.chart)[0][0].style.height.length - 2));

  // Extra padding calculated if x-axis title is present
  var extraPadding = self.xaxistitle != null && !d3.select(self.xaxistitle).empty() ? self.xaxistitle[0][0].clientHeight + 20 : 20;

  // Store the original bounding box and height of the chart without the x-axis slider chart
  if (!self.origBBox || !self.origHeight) {
    self.origBBox = d3.select(self.chart).select('svg').attr("viewBox");
    self.origHeight = d3.select(self.chart).select('svg').attr("height");
  }
  self.origBBox = d3.select(self.chart).select('svg').attr("viewBox");

  // Calculate the final chart height and position for storing the slider chart
  var height = (Number(self.origHeight) + self.size.sliderChartHeight + extraPadding + self.padding.sliderChart);

  // Set the chart svg viewBox and height accordingly to store the slider chart
  d3.select(self.chart).select('svg')
    .attr("viewBox", function() {
      var viewBox = self.origBBox.split(' ');
      return viewBox[0] + " " + viewBox[1] + " " + viewBox[2] + " " + height;
    })
    .attr("height", height);

  // Set the chart style height to have room for the slider chart
  if (startingHeight.toFixed() != (Number(self.origHeight) + self.size.sliderChartHeight + extraPadding + self.padding.sliderChart).toFixed()) {
    d3.select(this.chart)[0][0].style.height = height + "px";
  }

  // Store the original x and y-axis domain (we'll use these to draw the slider lines and position the slider box)
  var origdomain = self.xScale.domain();
  var origdomainrange = self.xScale.range();
  var origrange = self.yScale.domain();
  var bounds = self.min_max_x_values();
  var maxX = bounds[1];
  var minX = bounds[0];

  // Change the x and y-axis domain to the full range (for slider lines)
  self.xScale.domain([minX, maxX]);
  self.xScale.range([0, self.size.sliderChartWidth]);
  self.do_rebin();
  self.yScale.domain(self.getYAxisDomain());

  // Draw the elements for the slider chart
  if (d3.select("#sliderChart").empty()) {

    // G element of the slider chart
    self.sliderChart = d3.select("svg").append("g").attr("id", "sliderChart")
      .attr("transform", "translate(" + self.padding.leftComputed + "," + (480 + extraPadding + self.padding.sliderChart) + ")")
      // .on("mousemove", self.handleMouseMoveSliderChart());
      .on("touchstart", self.handleTouchStartSliderChart())
      .on("touchmove", self.handleTouchMoveSliderChart());

    // Plot area for data lines in slider chart
    self.sliderChartPlot = self.sliderChart.append("rect")
      .attr("id", "sliderchartarea"+self.chart.id )
      .attr("width", self.size.sliderChartWidth)
      .attr("height", self.size.sliderChartHeight)
      .style("fill", "#EEEEEE");

    // Chart body for slider (keeps the data lines)
    self.sliderChartBody = self.sliderChart.append("g")
      .attr("clip-path", "url(#sliderclip" + this.chart.id + ")");

    // Clip path for slider chart
    self.sliderChartClipPath = self.sliderChart.append('svg:clipPath')
        .attr("id", "sliderclip" + self.chart.id )
        .append("svg:rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", self.size.sliderChartWidth)
        .attr("height", self.size.sliderChartHeight);


    // For adding peaks into slider chart 
    // self.sliderPeakVis = d3.select("#sliderChart").append('g')
    //   .attr("id", "sliderPeakVis")
    //   .attr("class", "peakVis")
    //   .attr("transform", "translate(0,0)")
    //   .attr("clip-path", "url(#sliderclip" + this.chart.id + ")");

  } else {
    // Adjust the width, height, and transform for the slider chart elements (for resizing)
    self.sliderChartPlot.attr("width", self.size.sliderChartWidth)
      .attr("height", self.size.sliderChartHeight);
    self.sliderChartClipPath.attr("width", self.size.sliderChartWidth)
      .attr("height", self.size.sliderChartHeight);

    // self.sliderChart.attr("transform", "translate(0," + (self.size.height + extraPadding + self.padding.sliderChart) + ")");
  }

  // Commented out for adding peaks into slider chart sometime later
  // self.sliderPeakVis.selectAll('*').remove();
  // self.peakVis.select(function() {
  //   this.childNodes.forEach(function(path) {
  //     path = d3.select(path);
  //     self.sliderPeakVis.append('path')
  //       .attr("d", path.attr('d'))
  //       .attr("class", path.attr('class'))
  //       .attr("fill-opacity", path.attr('fill-opacity'))
  //       .attr("stroke-width", path.attr('stroke-width'))
  //       .attr("stroke", path.attr("stroke"))
  //       .attr("transform", "translate(0," + (self.size.height + extraPadding + self.padding.sliderChart) + ")");
  //   });
  // });

  // Add the slider draggable box and edges
  if (!self.sliderBox) {
    // Slider box
    self.sliderBox = self.sliderChart.append("rect")
      .attr("id", "sliderBox")
      .attr("class", "sliderBox")
      .attr("height", self.size.sliderChartHeight)
      .on("mousedown", self.handleMouseDownSliderBox())
      .on("touchstart", self.handleTouchStartSliderBox())
      .on("touchmove", self.handleTouchMoveSliderChart());

    // Left slider drag region
    self.sliderDragLeft = self.sliderChart.append("rect")
      .attr("id", "sliderDragLeft")
      .attr("class", "sliderDragRegion")
      .attr("rx", 2)
      .attr("ry", 2)
      .on("mousedown", self.handleMouseDownLeftSliderDrag())
      .on("mousemove", self.handleMouseMoveLeftSliderDrag(false))
      .on("mouseout", self.handleMouseOutLeftSliderDrag())
      .on("touchstart", self.handleTouchStartLeftSliderDrag())
      .on("touchmove", self.handleTouchMoveLeftSliderDrag(false));

    // Right slider drag region
    self.sliderDragRight = self.sliderChart.append("rect")
      .attr("id", "sliderDragRight")
      .attr("class", "sliderDragRegion")
      .attr("rx", 2)
      .attr("ry", 2)
      .on("mousedown", self.handleMouseDownRightSliderDrag())
      .on("mousemove", self.handleMouseMoveRightSliderDrag(false))
      .on("mouseout", self.handleMouseOutRightSliderDrag())
      .on("touchstart", self.handleTouchStartRightSliderDrag())
      .on("touchmove", self.handleTouchMoveRightSliderDrag(false));
  }

  var sliderBoxX = self.xScale(origdomain[0]);
  var sliderBoxWidth = self.xScale(origdomain[1]) - self.xScale(origdomain[0]);

  // Adjust the position of the slider box to the particular zoom region
  self.sliderBox.attr("x", sliderBoxX)
    .attr("width", sliderBoxWidth);
    // .on("mousemove", self.handleMouseMoveSliderChart());

  self.sliderDragLeft.attr("width", self.size.sliderChartWidth/100)
    .attr("height", self.size.sliderChartHeight/2.3);
  self.sliderDragRight.attr("width", self.size.sliderChartWidth/100)
    .attr("height", self.size.sliderChartHeight/2.3);

  self.sliderDragLeft.attr("x", sliderBoxX - Number(self.sliderDragLeft.attr("width"))/2)
    .attr("y", self.size.sliderChartHeight/2 - Number(self.sliderDragLeft.attr("height"))/2);

  self.sliderDragRight.attr("x", (sliderBoxX + sliderBoxWidth) - (Number(self.sliderDragRight.attr("width"))/2))
    .attr("y", self.size.sliderChartHeight/2 - Number(self.sliderDragRight.attr("height"))/2);

  drawSliderChartLines();
  drawDragRegionLines();

  // Restore the original x and y-axis domain
  self.xScale.domain(origdomain);
  self.xScale.range(origdomainrange);
  self.do_rebin();
  self.yScale.domain(origrange);
}

SpectrumChartD3.prototype.cancelXAxisSliderChart = function() {
  var self = this;
  self.size.sliderChartHeight = self.size.height / 10;

  var height = Number(d3.select(this.chart)[0][0].style.height.substring(0, d3.select(this.chart)[0][0].style.height.length - 2));

  if (!self.origBBox || !self.origHeight) {
    return;
  }

  if (self.sliderChart) {
    self.sliderChart.remove();
    self.sliderChartBody.remove();
    self.sliderChartPlot.remove();  
    self.sliderChartClipPath.remove();
    self.sliderBox.remove();
    d3.selectAll(".sliderDragRegion").remove();

    if (self.rawData && self.rawData.spectra && self.rawData.spectra.length)
      for (var i = 0; i < self.rawData.spectra.length; ++i)
        if (self['sliderLine'+i]) {
          self['sliderLine'+i].remove();
          self['sliderLine'+i] = null;
        }

    self.sliderChart = null;
    self.sliderChartPlot = null;
    self.sliderChartBody = null;
    self.sliderChartClipPath = null;
    self.sliderBox = null;
  }

  d3.select(self.chart)[0][0].style.height = self.origHeight + "px";
  d3.select(self.chart).select('svg').attr("viewBox", function() {
    var viewBoxAttrs = self.origBBox.split(' ');
    viewBoxAttrs[3] = self.origHeight.toString();
    return viewBoxAttrs.join(' ');
  });
  d3.select(self.chart).select('svg').attr("height", self.origHeight);

  self.origBBox = null;
  self.origHeight = null;

  self.sliderBoxDown = false;
  self.leftDragRegionDown = false;
  self.rightDragRegionDown = false;
  self.sliderChartMouse = null;
  self.savedSliderMouse = null;
}

SpectrumChartD3.prototype.handleMouseDownSliderBox = function() {
  var self = this;

  return function() {
    self.sliderBoxDown = true;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;

    // Initially set the escape key flag false
    self.escapeKeyPressed = false;    

    d3.select(document.body).style("cursor", "move");
  }
}

SpectrumChartD3.prototype.handleMouseMoveSliderChart = function() {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();
    // console.log("sliderboxdown = ", self.sliderBoxDown);
    // console.log("leftDragRegionDown = ", self.leftDragRegionDown);
    // console.log("rightDragRegionDown = ", self.rightDragRegionDown);

    if (self.leftDragRegionDown) {
      return self.handleMouseMoveLeftSliderDrag(true)();
    }

    if (self.rightDragRegionDown) {
      return self.handleMouseMoveRightSliderDrag(true)();
    }

    if (self.sliderBoxDown) {
      d3.select(document.body).style("cursor", "move");
      var m = d3.mouse(self.sliderChart[0][0]);
      var origdomain = self.xScale.domain();
      var origdomainrange = self.xScale.range();
      var bounds = self.min_max_x_values();
      var maxX = bounds[1];
      var minX = bounds[0];

      if (!self.sliderChartMouse) {
        self.sliderChartMouse = m;
      }

      var sliderBoxX = Number(self.sliderBox.attr("x"));
      var sliderBoxWidth = Number(self.sliderBox.attr("width"));
      var sliderDragRegionWidth = 3;
      var x = Math.min( self.size.sliderChartWidth - sliderBoxWidth, Math.max(0, sliderBoxX + (m[0] - self.sliderChartMouse[0])) );

      if ((sliderBoxX == 0 || sliderBoxX + sliderBoxWidth == self.size.sliderChartWidth)) {
        if (!self.savedSliderMouse)
          self.savedSliderMouse = m;
      }

      if (self.savedSliderMouse && m[0] != self.savedSliderMouse[0]) {
        if (sliderBoxX == 0 && m[0] < self.savedSliderMouse[0]) return;
        else if (sliderBoxX + sliderBoxWidth == self.size.sliderChartWidth && m[0] > self.savedSliderMouse[0]) return;
        else self.savedSliderMouse = null;
      }

      self.xScale.domain([minX, maxX]);
      self.xScale.range([0, self.size.sliderChartWidth]);
      self.sliderBox.attr("x", x);
      self.sliderDragLeft.attr("x", x);
      self.sliderDragRight.attr("x", x + sliderBoxWidth - sliderDragRegionWidth);

      origdomain = [ self.xScale.invert(x), self.xScale.invert(x + sliderBoxWidth) ];
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      self.redraw()();

      self.sliderChartMouse = m;
    }
  }
}

SpectrumChartD3.prototype.handleMouseDownLeftSliderDrag = function() {
  var self = this;

  return function() {
    self.leftDragRegionDown = true;

    // Initially set the escape key flag false
    self.escapeKeyPressed = false;    
  }
}

SpectrumChartD3.prototype.handleMouseMoveLeftSliderDrag = function(redraw) {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();

    if (self.sliderBoxDown) {
      return self.handleMouseMoveSliderChart()();
    }

    d3.select(document.body).style("cursor", "ew-resize");

    if (!self.leftDragRegionDown || !redraw) {
      return;
    }

    var m = d3.mouse(self.sliderChart[0][0]);
    var origdomain = self.xScale.domain();
    var origdomainrange = self.xScale.range();
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];
    var x = Math.max(m[0], 0);

    var sliderBoxX = self.xScale(origdomain[0]);
    var sliderBoxWidth = Number(self.sliderBox.attr("width"));
    var sliderDragRegionWidth = 3;
    var sliderDragPadding = 1;

    self.xScale.domain([minX, maxX]);
    self.xScale.range([0, self.size.sliderChartWidth]);

    if (m[0] > Number(self.sliderDragRight.attr("x") - sliderDragPadding)) {
      self.xScale.domain(origdomain);
      return;
    }

    self.sliderBox.attr("x", x);
    self.sliderDragLeft.attr("x", x);
    origdomain[0] = self.xScale.invert(x);

    self.xScale.domain(origdomain);
    self.xScale.range(origdomainrange);
    self.redraw()();
  }
}

SpectrumChartD3.prototype.handleMouseOutLeftSliderDrag = function() {
  var self = this;

  return function() {
    if (!self.leftDragRegionDown) {
      d3.select(document.body).style("cursor", "default");
    }
  }
}

SpectrumChartD3.prototype.handleMouseDownRightSliderDrag = function() {
  var self = this;

  return function() {
    self.rightDragRegionDown = true;

    // Initially set the escape key flag false
    self.escapeKeyPressed = false;    
  }
}

SpectrumChartD3.prototype.handleMouseMoveRightSliderDrag = function(redraw) {
  var self = this;

  return function() {

    if (self.sliderBoxDown) {
      return; //self.handleMouseMoveSliderChart()();
    }

    d3.event.preventDefault();
    d3.event.stopPropagation();
    d3.select('body').style("cursor", "ew-resize");

    if (!self.rightDragRegionDown || !redraw) {
      return;
    }

    var m = d3.mouse(self.sliderChart[0][0]);
    var origdomain = self.xScale.domain();
    var origdomainrange = self.xScale.range();
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];
    var x = Math.min(m[0], self.size.sliderChartWidth);

    var sliderBoxX = self.xScale(origdomain[0]);
    var sliderBoxWidth = Number(self.sliderBox.attr("width"));
    var sliderDragRegionWidth = 3;
    var sliderDragPadding = 1;

    self.xScale.domain([minX, maxX]);
    self.xScale.range([0, self.size.sliderChartWidth]);

    if (m[0] - sliderDragRegionWidth < Number(self.sliderDragLeft.attr("x")) + Number(self.sliderDragLeft.attr("width"))) {
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      return;
    }

    self.sliderBox.attr("width", Math.abs(x - Number(self.sliderDragRight.attr("x"))));
    self.sliderDragRight.attr("x", x - sliderDragRegionWidth);
    origdomain[1] = self.xScale.invert(x);

    self.xScale.domain(origdomain);
    self.xScale.range(origdomainrange);
    self.redraw()();
  }
}

SpectrumChartD3.prototype.handleMouseOutRightSliderDrag = function() {
  var self = this;

  return function() {
    if (!self.rightDragRegionDown) {
      d3.select(document.body).style("cursor", "default");
    }
  }
}

SpectrumChartD3.prototype.handleTouchStartSliderBox = function() {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();

    self.sliderBoxDown = true;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;

    var t = [d3.event.pageX, d3.event.pageY];
    var touchError = 25;
    var x1 = self.sliderBox[0][0].getBoundingClientRect().left;
    var x2 = self.sliderBox[0][0].getBoundingClientRect().right;

    if (d3.event.changedTouches.length !== 1)
      return;

    if (x1 - touchError <= t[0] && t[0] <= x1 + touchError) {
      self.sliderBoxDown = false;
      self.leftDragRegionDown = true;

    } else if (x2 - touchError <= t[0] && t[0] <= x2 + touchError) {
      self.sliderBoxDown = false;
      self.rightDragRegionDown = true;
    }
  }
}

SpectrumChartD3.prototype.handleTouchStartSliderChart = function() {
  var self = this;

  return function() {
    var t = [d3.event.pageX, d3.event.pageY];
    var touchError = 15;
    var x1 = self.sliderBox[0][0].getBoundingClientRect().left;
    var x2 = self.sliderBox[0][0].getBoundingClientRect().right;

    if (d3.event.changedTouches.length !== 1)
      return;

    if (x1 - touchError <= t[0] && t[0] <= x1 + touchError) {
      self.sliderBoxDown = false;
      self.leftDragRegionDown = true;

    } else if (x2 - touchError <= t[0] && t[0] <= x2 + touchError) {
      self.sliderBoxDown = false;
      self.rightDragRegionDown = true;
    }
  }
}

SpectrumChartD3.prototype.handleTouchMoveSliderChart = function() {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();

    if (self.leftDragRegionDown) {
      return self.handleTouchMoveLeftSliderDrag(true)();
    }

    if (self.rightDragRegionDown) {
      return self.handleTouchMoveRightSliderDrag(true)();
    }

    if (self.sliderBoxDown) {
      d3.select(document.body).style("cursor", "move");
      var t = d3.touches(self.sliderChart[0][0]);

      if (t.length !== 1)
        return;

      t = t[0];
      var origdomain = self.xScale.domain();
      var origdomainrange = self.xScale.range();
      var bounds = self.min_max_x_values();
      var maxX = bounds[1];
      var minX = bounds[0];

      if (!self.sliderChartTouch) {
        self.sliderChartTouch = t;
      }

      var sliderBoxX = Number(self.sliderBox.attr("x"));
      var sliderBoxWidth = Number(self.sliderBox.attr("width"));
      var sliderDragRegionWidth = 3;
      var x = Math.min( self.size.sliderChartWidth - sliderBoxWidth, Math.max(0, sliderBoxX + (t[0] - self.sliderChartTouch[0])) );

      if ((sliderBoxX == 0 || sliderBoxX + sliderBoxWidth == self.size.sliderChartWidth)) {
        if (!self.savedSliderTouch)
          self.savedSliderTouch = t;
      }

      if (self.savedSliderTouch && t[0] != self.savedSliderTouch[0]) {
        if (sliderBoxX == 0 && t[0] < self.savedSliderTouch[0]) return;
        else if (sliderBoxX + sliderBoxWidth == self.size.sliderChartWidth && t[0] > self.savedSliderTouch[0]) return;
        else self.savedSliderTouch = null;
      }

      self.xScale.domain([minX, maxX]);
      self.xScale.range([0, self.size.sliderChartWidth]);
      self.sliderBox.attr("x", x);
      self.sliderDragLeft.attr("x", x);
      self.sliderDragRight.attr("x", x + sliderBoxWidth - sliderDragRegionWidth);

      origdomain = [ self.xScale.invert(x), self.xScale.invert(x + sliderBoxWidth) ];
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      self.redraw()();

      self.sliderChartTouch = t;

      // IMPORTANT: To translate the current x-scale into the current zoom object
      self.zoom.x(self.xScale);
    }
  }
}

SpectrumChartD3.prototype.handleTouchStartLeftSliderDrag = function() {
  var self = this;

  return function() {
    self.leftDragRegionDown = true;
  }
}

SpectrumChartD3.prototype.handleTouchMoveLeftSliderDrag = function(redraw) {
  var self = this;

  return function() {
    d3.event.preventDefault();
    d3.event.stopPropagation();

    if (self.sliderBoxDown) {
      return;
    }

    if (!self.leftDragRegionDown) {
      return;
    }

    var t = d3.touches(self.sliderChart[0][0]);
    if (t.length !== 1)
      return;

    t = t[0];
    var origdomain = self.xScale.domain();
    var origdomainrange = self.xScale.range();
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];
    var x = Math.max(t[0], 0);

    var sliderBoxX = self.xScale(origdomain[0]);
    var sliderBoxWidth = Number(self.sliderBox.attr("width"));
    var sliderDragRegionWidth = 3;
    var sliderDragPadding = 1;

    self.xScale.domain([minX, maxX]);
    self.xScale.range([0, self.size.sliderChartWidth]);

    if (t[0] > Number(self.sliderDragRight.attr("x") - sliderDragPadding)) {
      self.xScale.domain(origdomain);
      return;
    }

    self.sliderBox.attr("x", x);
    self.sliderDragLeft.attr("x", x);
    origdomain[0] = self.xScale.invert(x);

    self.xScale.domain(origdomain);
    self.xScale.range(origdomainrange);
    self.redraw()();

    // IMPORTANT: To translate the current x-scale into the current zoom object
    self.zoom.x(self.xScale);
  }
}

SpectrumChartD3.prototype.handleTouchStartRightSliderDrag = function() {
  var self = this;

  return function() {
    self.rightDragRegionDown = true;
  }
}

SpectrumChartD3.prototype.handleTouchMoveRightSliderDrag = function(redraw) {
  var self = this;

  return function() {
    if (self.sliderBoxDown) {
      return;
    }

    d3.event.preventDefault();
    d3.event.stopPropagation();

    if (!self.rightDragRegionDown) {
      return;
    }

    var t = d3.touches(self.sliderChart[0][0]);
    if (t.length !== 1)
      return;
  
    t = t[0];
    var origdomain = self.xScale.domain();
    var origdomainrange = self.xScale.range();
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];
    var x = Math.min(t[0], self.size.sliderChartWidth);

    var sliderBoxX = self.xScale(origdomain[0]);
    var sliderBoxWidth = Number(self.sliderBox.attr("width"));
    var sliderDragRegionWidth = 3;
    var sliderDragPadding = 1;

    self.xScale.domain([minX, maxX]);
    self.xScale.range([0, self.size.sliderChartWidth]);

    if (t[0] - sliderDragRegionWidth < Number(self.sliderDragLeft.attr("x")) + Number(self.sliderDragLeft.attr("width"))) {
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      return;
    }

    self.sliderBox.attr("width", Math.abs(x - Number(self.sliderDragRight.attr("x"))));
    self.sliderDragRight.attr("x", x - sliderDragRegionWidth);
    origdomain[1] = self.xScale.invert(x);

    self.xScale.domain(origdomain);
    self.xScale.range(origdomainrange);
    self.redraw()();

    // IMPORTANT: To translate the current x-scale into the current zoom object
    self.zoom.x(self.xScale);
  }
}



/* SpectrumChartD3 - Animation Methods */
SpectrumChartD3.prototype.redrawZoomXAnimation = function(targetDomain) {
  var self = this;

  // Cancel animation if showAnimation option not checked
  if (!self.options.showAnimation) { return; }

  return function() {
    // Cancel the animation once reached desired target domain
    if (self.currentDomain == null || targetDomain == null || (self.currentDomain[0] == targetDomain[0] && self.currentDomain[1] == targetDomain[1])) {
      console.log("Time for animation = ", Math.floor(Date.now()) - self.startAnimationZoomTime, " ms");
      self.handleCancelAnimationZoom();
      return;
    }

    // Use fraction of time elapsed to calculate how far we will zoom in this frame
    var animationFractionTimeElapsed = Math.min( Math.max((Math.floor(Date.now()) - self.startAnimationZoomTime) / self.animation.duration), 1 );

    // Set x-axis domain to new values
    self.xScale.domain([ 
      Math.min( self.savedDomain[0] + (animationFractionTimeElapsed * (targetDomain[0] - self.savedDomain[0])), targetDomain[0] ),
      Math.max( self.savedDomain[1] - (animationFractionTimeElapsed * (self.savedDomain[1] - targetDomain[1])), targetDomain[1] )
    ]);
    self.currentDomain = self.xScale.domain();

    // Redraw and request a new animation frame
    self.redraw()();
    self.zoomAnimationID = requestAnimationFrame(self.redrawZoomXAnimation(targetDomain));
  }
}

SpectrumChartD3.prototype.redrawZoomInYAnimation = function(targetDomain,redraw) {
  var self = this;

  function roundTo3DecimalPlaces(num) { return Math.round(num * 1000) / 1000; }

  // Cancel animation if showAnimation option not checked
  if (!self.options.showAnimation) { return; }

  return function() {
    // Cancel the animation once reached desired target domain
    if (self.currentDomain == null || targetDomain ==  null || 
      (roundTo3DecimalPlaces(self.currentDomain[0]) == roundTo3DecimalPlaces(targetDomain[0]) && 
        roundTo3DecimalPlaces(self.currentDomain[1]) == roundTo3DecimalPlaces(targetDomain[1]))) {

      console.log("Time for animation = ", Math.floor(Date.now()) - self.startAnimationZoomTime, " ms");
      self.handleCancelAnimationZoom();
      return;
    }

    // Use fraction of time elapsed to calculate how far we will zoom in this frame
    var animationFractionTimeElapsed = Math.min( Math.max((Math.floor(Date.now()) - self.startAnimationZoomTime) / self.animation.duration), 1 );

    // Set y-axis domain to new values
    self.yScale.domain([ 
      Math.max( self.savedDomain[0] - (animationFractionTimeElapsed * (self.savedDomain[0] - targetDomain[0])), targetDomain[0] ),
      Math.min( self.savedDomain[1] + (animationFractionTimeElapsed * (targetDomain[1] - self.savedDomain[1])), targetDomain[1] )
    ]);
    self.currentDomain = self.yScale.domain();

    // Redraw and request a new animation frame
    redraw();
    self.zoomAnimationID = requestAnimationFrame(self.redrawZoomInYAnimation(targetDomain, redraw));
  }
}

SpectrumChartD3.prototype.redrawZoomOutYAnimation = function(targetDomain, redraw) {
  var self = this;

  function roundTo3DecimalPlaces(num) { return Math.round(num * 1000) / 1000; }

  // Cancel animation if showAnimation option not checked
  if (!self.options.showAnimation) { return; }

  return function() {
    // Cancel the animation once reached desired target domain
    if (self.currentDomain == null || targetDomain ==  null || 
      (roundTo3DecimalPlaces(self.currentDomain[0]) == roundTo3DecimalPlaces(targetDomain[0]) && 
        roundTo3DecimalPlaces(self.currentDomain[1]) == roundTo3DecimalPlaces(targetDomain[1]))) {

      console.log("Time for animation = ", Math.floor(Date.now()) - self.startAnimationZoomTime, " ms");
      self.handleCancelAnimationZoom();
      return;
    }

    // Use fraction of time elapsed to calculate how far we will zoom in this frame
    var animationFractionTimeElapsed = Math.min( Math.max((Math.floor(Date.now()) - self.startAnimationZoomTime) / self.animation.duration), 1 );

    // Set y-axis domain to new values
    self.yScale.domain([ 
      Math.min( self.savedDomain[0] + (animationFractionTimeElapsed * (targetDomain[0] - self.savedDomain[0])), targetDomain[0] ),
      Math.max( self.savedDomain[1] - (animationFractionTimeElapsed * (self.savedDomain[1] - targetDomain[1])), targetDomain[1] )
    ]);
    self.currentDomain = self.yScale.domain();

    // Redraw and request a new animation frame
    redraw();
    self.zoomAnimationID = requestAnimationFrame(self.redrawZoomOutYAnimation(targetDomain, redraw));
  }
}

SpectrumChartD3.prototype.handleCancelAnimationZoom = function() {
  var self = this;

  // Cancel the animation frames
  if (self.zoomAnimationID != null) {
    cancelAnimationFrame(self.zoomAnimationID);
  }

  // Set animation properties to null
  self.zoomAnimationID = null;
  self.currentDomain = null;
  self.savedDomain = null;
  self.startAnimationZoomTime = null;
}
 





/* Sets x-axis and y-axis drag for the chart. These are actions done by clicking and dragging
    one of the labels of the x-axis or y-axis.
*/

SpectrumChartD3.prototype.xaxisDrag = function() {
  var self = this;
  return function(d) {
    //This function is called once when you click on an x-axis label (which you can then start dragging it)
    //  And NOT when you click on the chart and drag it to pan
    console.log( 'xaxisDrag work' );
    document.onselectstart = function() { return false; };
    var p = d3.mouse(self.vis[0][0]);

    if (self.xScale.invert(p[0]) > 0)           // set self.downx equal to value of your mouse pos
      self.downx = self.xScale.invert(p[0]);
    else                                        // if mouse pos < 0, use 0.1 as the value (value of 0 is buggy)
      self.downx = 0.1;
    //console.log("p=" + p );
    //console.log("self.downx=" + self.downx);
    //console.log("self.vis[0][0]=" + self.vis[0][0]);
  }
}


SpectrumChartD3.prototype.yaxisDrag = function(d) {
  var self = this;
  return function(d) {
    console.log('yaxisDrag work');
    document.onselectstart = function() { return false; };
    var p = d3.mouse(self.vis[0][0]);
    self.downy = self.yScale.invert(p[1]);
  }
}


SpectrumChartD3.prototype.setShowUserLabels = function(d) {
  this.options.showUserLabels = d;
  this.redraw()();
}

SpectrumChartD3.prototype.setShowPeakLabels = function(d) {
  this.options.showPeakLabels = d;
  this.redraw()();
}

SpectrumChartD3.prototype.setShowNuclideNames = function(d) {
  this.options.showNuclideNames = d;
  this.redraw()();
}

SpectrumChartD3.prototype.setShowNuclideEnergies = function(d) {
  this.options.showNuclideEnergies = d;
  this.redraw()();
}



/* Set peak marker options for chart
    
    Currently have:
      - Compton Peaks (including angle)
      - Escape peaks
      - Sum peaks
      - Compton edges
*/
SpectrumChartD3.prototype.setComptonEdge = function(d) {
  this.options.showComptonEdge = d;
  if ( d ) {
    this.updateFeatureMarkers();
  }
}

SpectrumChartD3.prototype.setComptonPeakAngle = function(d) {
  var value = Number(d);
  if (value != NaN && 0 <= value && value <= 180 ) 
    this.comptonPeakAngle = value;
  else {
    this.comptonPeakAngle = 180;  // default angle is set to 180 degrees
  }
  this.updateFeatureMarkers();
}

SpectrumChartD3.prototype.setComptonPeaks = function(d) {
  this.options.showComptonPeaks = d;
  if ( !this.comptonPeakAngle )
    this.comptonPeakAngle = 180;
  if ( d ) {
    this.updateFeatureMarkers();
  }
}

SpectrumChartD3.prototype.setEscapePeaks = function(d) {
  this.options.showEscapePeaks = d;
  if ( d ) {
    this.updateFeatureMarkers();
  }
}

SpectrumChartD3.prototype.setSumPeaks = function(d) {
  this.options.showSumPeaks = d;
  if ( d ) {
    this.updateFeatureMarkers();
  }
}

// Sets the mouse stats view inside the chart
SpectrumChartD3.prototype.setShowMouseStats = function(d) {
  this.options.showMouseStats = d;
  this.mouseInfo.style("display", d ? null : "none");
  if( d )
    this.updateMouseCoordText();
}

// Sets whether or not peaks are highlighted
SpectrumChartD3.prototype.setShowPeaks = function(series,show) {
  this.options.drawPeaks[series] = show;
  this.redraw()();
}

// Sets the title of the graph
SpectrumChartD3.prototype.setTitle = function(title,dontRedraw) {
  
  var titleh = 0;
  if( (title == null || typeof title !== 'string') || title.length === 0 ){
    this.options.title = null;
    d3.select('#chartTitle').remove();
  } else {
    if( this.options.title )
      titleh = d3.select('svg').selectAll(".title").text( title ).node().getBBox().height;
    else
      titleh = d3.select('svg').append("text")
          .attr("id", "chartTitle")
          .attr("class", "title")
          .text(title)
          .attr("x", this.cx/2)
          .attr("dy", this.padding.title)
          .style("text-anchor","middle")
          .node().getBBox().height;
    this.options.title = title;
  }
  this.handleResize( dontRedraw ); 
}

SpectrumChartD3.prototype.setXRangeArrows = function(d) {
  var self = this;
  this.options.showXRangeArrows = d;
  this.drawXAxisArrows(d);
}

SpectrumChartD3.prototype.setShowXAxisSliderChart = function(d) {
  this.options.showXAxisSliderChart = d;
  this.drawXAxisSliderChart();
}

SpectrumChartD3.prototype.setShowAnimation = function(d) {
  this.options.showAnimation = d;
}

SpectrumChartD3.prototype.setAnimationDuration = function(d) {
  this.animation.duration = d;
}

SpectrumChartD3.prototype.setWheelScrollYAxis = function(d) {
  this.options.wheelScrollYAxis = d;
}

SpectrumChartD3.prototype.setDropDownSpectrumScaleFactor = function(d) {
  var self = this;
  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  self.currentDropDownScaleFactorSpectrum = d == "None" ? null : d;

  if (d != "None") {
    var valueChanged = false;
    self.rawData.spectra.forEach(function(spectrum) {
      if (spectrum.title && spectrum.title.toUpperCase() == self.currentDropDownScaleFactorSpectrum.toUpperCase() && d != "None" && !valueChanged) {
        valueChanged = true;

      var currentsfinput;
      var spectrumScaleFactor = spectrum.yScaleFactor;

      if (currentsfinput = document.getElementById("current-sf"))
        currentsfinput.value = needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed();
      return;
    }
    });
  }
}

SpectrumChartD3.prototype.setSpectrumScaleFactorWidget = function(d) {
  this.options.scaleBackgroundSecondary = d;
  this.drawScalerBackgroundSecondary();
}

SpectrumChartD3.prototype.setSpectrumScaleFactor = function(d) {
  var self = this;
  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;
  if (!self.currentlyAdjustingSpectrumScale && !self.currentDropDownScaleFactorSpectrum)
    return;

  var spectrumToCheck = self.currentlyAdjustingSpectrumScale ? self.currentlyAdjustingSpectrumScale : self.currentDropDownScaleFactorSpectrum;

  // Check for the which corresponding spectrum line is the to be adjusted
  var spectrumToBeAdjusted = null;
  var linei = null;
  self.rawData.spectra.forEach(function(spectrum,i) {
    if (spectrum.title && spectrum.title.toUpperCase() == spectrumToCheck.toUpperCase() && spectrumToBeAdjusted == null) {
      spectrumToBeAdjusted = spectrum;
      linei = i;
      return;
    }
  });

  if (spectrumToBeAdjusted == null || linei == null)
    return;
  if (!self['line'+linei] || self.vis.select("#spectrumline"+linei).empty())
    return;

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  /* Here is a modified, slightly more efficient version of redraw 
    specific to changing the scale factors for spectrums. 
  */
  function scaleFactorChangeRedraw(spectrum, linei) {
    self.updateLegend();
    self.rebinSpectrum(spectrum, linei);
    self.do_rebin();
    self.setYAxisDomain();
    self.drawYTicks();
    self.update();
    self.drawPeaks();
  }

  var spectrumScaleFactor = Math.min(Number(d), self.options.maxScaleFactor);
  spectrumToBeAdjusted.yScaleFactor = spectrumScaleFactor;
  scaleFactorChangeRedraw(spectrumToBeAdjusted, linei);

  if (spectrumToBeAdjusted.sliderText)
    spectrumToBeAdjusted.sliderText.text(spectrumToBeAdjusted.title + " Scale Factor: " + (needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed()));
  if (spectrumToBeAdjusted.sliderToggle && spectrumToBeAdjusted.scale)
    spectrumToBeAdjusted.sliderToggle.attr("cx", Math.max( Math.min(spectrumToBeAdjusted.scale(spectrumScaleFactor), spectrumToBeAdjusted.scale.range()[1]), spectrumToBeAdjusted.scale.range()[0] ));

  // Update the slider chart if needed
  if (self["sliderLine"+linei]) {
    var origdomain = self.xScale.domain();
    var origdomainrange = self.xScale.range();
    var origrange = self.yScale.domain();
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];

    // Change the x and y-axis domain to the full range (for slider lines)
    self.xScale.domain([minX, maxX]);
    self.xScale.range([0, self.size.sliderChartWidth]);
    self.do_rebin();
    self.yScale.domain(self.getYAxisDomain());

    self["sliderLine"+linei].attr("d", self["line"+linei](spectrumToBeAdjusted.points));

    // Restore the original x and y-axis domain
    self.xScale.domain(origdomain);
    self.xScale.range(origdomainrange);
    self.do_rebin();
    self.yScale.domain(origrange);
  }
}

SpectrumChartD3.prototype.setMaxScaleFactor = function(d) {
  var self = this;

  d = Math.max(Number(d), 0.1);

  if (self.options)
    this.options.maxScaleFactor = d;
  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  function needsDecimal(num) {
    return num % 1 != 0;
  }

  /* Here is a modified, slightly more efficient version of redraw 
    specific to changing the scale factors for spectrums. 
  */
  function scaleFactorChangeRedraw(spectrum, linei) {
    self.updateLegend();
    self.rebinSpectrum(spectrum, linei);
    self.do_rebin();
    self.setYAxisDomain();
    self.drawYTicks();
    self.update();
    self.drawPeaks();
  }

  self.rawData.spectra.forEach(function(spectrum,i) {
    // Don't scale the foreground
    if (i == 0)
      return;


    var spectrumScaleFactor = Math.min(spectrum.yScaleFactor, self.options.maxScaleFactor);
    var spectrumTitle = spectrum.title ? spectrum.title : "Spectrum " + i;
    spectrum.yScaleFactor = spectrumScaleFactor;

    var domain;
    if (!spectrum.scale || !spectrum.scaleAxis)
      return;

    domain = spectrum.scale.domain();
    spectrum.scale.domain([ domain[0], d ]);
    spectrum.scaleAxis.scale(spectrum.scale)
      .orient("bottom")
      .ticks(5,'f');
    d3.select("#" + spectrumTitle.replace(' ','-') + "SliderArea").call(spectrum.scaleAxis);

    if (spectrum.sliderText)
      spectrum.sliderText.text(spectrumTitle + " Scale Factor: " + (needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed()));
    if (spectrum.sliderToggle)
      spectrum.sliderToggle.attr("cx", Math.min(spectrum.scale(spectrum.yScaleFactor), spectrum.scale.range()[1]));

    scaleFactorChangeRedraw(spectrum, i);

    if (!self.vis.select("#spectrumline"+i).empty()) {
      self.vis.select("#spectrumline"+i).attr("d", self['line'+i](spectrum.points));
    }

    // Update the slider chart if needed
    if (self["sliderLine"+i]) {
      var origdomain = self.xScale.domain();
      var origdomainrange = self.xScale.range();
      var origrange = self.yScale.domain();
      var bounds = self.min_max_x_values();
      var maxX = bounds[1];
      var minX = bounds[0];

      // Change the x and y-axis domain to the full range (for slider lines)
      self.xScale.domain([minX, maxX]);
      self.xScale.range([0, self.size.sliderChartWidth]);
      self.do_rebin();
      self.yScale.domain(self.getYAxisDomain());

      self["sliderLine"+i].attr("d", self["line"+i](spectrum.points));

      // Restore the original x and y-axis domain
      self.xScale.domain(origdomain);
      self.xScale.range(origdomainrange);
      self.do_rebin();
      self.yScale.domain(origrange);
    }
  });
}


SpectrumChartD3.prototype.drawPeaks = function() {
  var self = this;

  self.peakVis.selectAll("*").remove();

  if( !this.rawData || !this.rawData.spectra ) 
    return;

  var minx = self.xScale.domain()[0], maxx = self.xScale.domain()[1];

  function offset_integral(roi,x0,x1){
    if( roi.type === 'NoOffset' )
      return 0.0;
    if( roi.type === 'External' ){
      console.log( 'External contrinuum not supported yet' );
      return 0.0;
    }

    x0 -= roi.referenceEnergy; x1 -= roi.referenceEnergy;
    var answer = 0.0;
    for( var i = 0; i < roi.coeffs.length; ++i )
      answer += (roi.coeffs[i]/(i+1)) * (Math.pow(x1,i+1) - Math.pow(x0,i+1));
    return Math.max( answer, 0.0 );
  }

  function roiPath(roi,points){
    var yl, yr;
    var lpx = self.xScale(roi.lowerEnergy), rpx = self.xScale(roi.upperEnergy);

    var paths = [];
    var bisector = d3.bisector(function(d){return d.x;});
    var xstartind = bisector.left( points, Math.max(roi.lowerEnergy,minx) );
    var xendind = bisector.right( points, Math.min(roi.upperEnergy,maxx) );
    if( xstartind>= (points.length-2) )
      return paths;
      
    if( xendind >= (points.length-2) )
      xendind = points.length - 2;

    //The continuum values used for the first and last bin of the ROI are fudged
    //  for now...  To be fixed

    //XXX - Need to check continuum type!!!

    var thisy = null, thisx = null, m, s, peak_area, cont_area;
    
    
    var firsty = offset_integral( roi, points[xstartind-(xstartind?1:0)].x, points[xstartind+(xstartind?0:1)].x );
    
    
    paths[0] = "M" + self.xScale(points[xstartind].x) + "," + self.yScale(firsty) + " L";
    for( var j = 0; j < roi.peaks.length; ++j )
      paths[j+1] = "";

    for( var i = xstartind; i < xendind; ++i ) {
      thisx = 0.5*(points[i].x + points[i+1].x);
      thisy = offset_integral( roi, points[i].x, points[i+1].x );
      paths[0] += " " + self.xScale(thisx) + "," + self.yScale(thisy);
      for( var j = 0; j < roi.peaks.length; ++j ) {
        m = roi.peaks[j].Centroid[0];
        s = roi.peaks[j].Width[0];
        if( thisx > (m - 5*s) && thisx < (m+5*s) ){
          if( !paths[j+1].length )
            paths[j+1] = "M" + self.xScale(thisx) + "," + self.yScale(thisy) + " L";
          else
            paths[j+1] += " " + self.xScale(thisx) + "," + self.yScale(thisy);
        }
      }
    }

    function erf(x) {
      //http://stackoverflow.com/questions/14846767/std-normal-cdf-normal-cdf-or-error-function
      var sign = (x >= 0) ? 1 : -1; // save the sign of x
      x = Math.abs(x);
      var a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027
          a5 = 1.061405429, p  = 0.3275911;
      var t = 1.0/(1.0 + p*x);
      var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
      return sign * y; // erf(-x) = -erf(x);
    }

    function gaus_integral( peak, x0, x1 ) {
      var peak_mean = peak.Centroid[0], peak_sigma = peak.Width[0], peak_amplitude = peak.Amplitude[0];
         //peak.LandauAmplitude, peak.LandauMode, peak.LandauSigma,
      var sqrt2 = 1.414213562373095;
      var erflowarg = (x0-peak_mean)/(sqrt2*peak_sigma);
      var erfhigharg = (x1-peak_mean)/(sqrt2*peak_sigma);
      return peak_amplitude * 0.5 * (erf(erfhigharg) - erf(erflowarg));
    };

    for( var i = xendind - 1; i >= xstartind; --i ) {
      peak_area = 0.0;
      thisx = 0.5*(points[i].x + points[i+1].x);
      cont_area = offset_integral( roi, points[i].x, points[i+1].x );

      roi.peaks.forEach( function(peak,peakn){
        if( peak.type !== 'GaussianDefined' ){
          console.log( 'Need to implement peak.type==' + peak.type );
          return;
        }
        if( peak.skewType !== 'NoSkew' )
          console.log( 'Need to implement peak skew type ' + peak.skewType );
        var area = gaus_integral( peak, points[i].x, points[i+1].x );
        peak_area += area;

        m = peak.Centroid[0];
        s = peak.Width[0];
        if( thisx > (m - 5*s) && thisx < (m+5*s) ){
          paths[peakn+1] += " " + self.xScale(thisx) + "," + self.yScale(cont_area + area);
        }
      });

      paths[0] += " " + self.xScale(thisx) + "," + self.yScale( peak_area + cont_area);
    }

    for( var i = 0; i < paths.length; ++i ) {
      var ind = paths[i].indexOf("L");
      if( ind > 2 )
        paths[i] += " " + paths[i].substr(1,ind-2);
    }
    //path += " " + self.xScale(points[xstartind].x) + "," + self.yScale(firsty);

    return paths;
  }//function roiPath(roi)

  function draw_roi(roi,specindex,spectrum) {
      /* roi:
      {
        type:'Linear',
        lowerEnergy:164.878,
        upperEnergy:203.311,
        referenceEnergy:164.878,
        coeffs:[49.2491,-0.365196],
        coeffUncerts:[2.88764,0.13389],
        fitForCoeff:[true,true],
        peaks: [
          {
            type:'GaussianDefined',
            skewType:'NoSkew',
            Centroid:[182.846,1.80707,true],
            Width:[7.17482,0.979233,true],
            Amplitude:[288.679,44.8641,true],
            LandauAmplitude:[0,-1,false],
            LandauMode:[0,-1,false],
            LandauSigma:[0,-1,false],
            Chi2:[0.573405,-1,false],
            forCalibration:true,
            forSourceFit:true,
            type:'NormalGamma',
            nuclide: {
              name: 'Mo99',
              decayParent:'Mo99',
              decayChild:'Tc99',
              DecayGammaEnergy:181.063
            }
          }]
        }
        */

        if( roi.type !== 'NoOffset' && roi.type !== 'Constant'
            && roi.type !== 'Linear' && roi.type !== 'Quardratic'
            && roi.type !==  'Cubic'  ){
          console.log( 'unrecognized roi.type: ' + roi.type );
          return;
        }

        if( roi.lowerEnergy > maxx || roi.upperEnergy < minx )
          return;

        if (!spectrum) {
          console.log("No spectrum specified to draw peaks");
          return;
        }

        var paths = roiPath( roi, spectrum.points, specindex );
        paths.forEach( function(p,num){
          if( num===1 && paths.length===2 )
            return;
          var path = self.peakVis.append("path")
                         .attr("d", p );

          var labels = [],
              labelAlreadyAdded = false;

          if (!self.peakLabelArray)                       /* This is a hack for gettng the correct highlighted peak to correspond with its label. */
            self.peakLabelArray = [];                     /* Declare an array of tuples that have a corresponding peak DOM with a label text*/
          if (!self.leftoverPeakLabels)
            self.leftoverPeakLabels = [];

          if (self.leftoverPeakLabels.length > 0) {                                                    // If there are still leftover labels not corresponding with a peak, then it
            self.peakLabelArray.push( { "path": path, "label": self.leftoverPeakLabels.shift() } );  // is the current peak's label
            labelAlreadyAdded = true;
          }
          else
            for (i = 0; i < roi.peaks.length; i++) {                                                      // Draw a label for each peak inside an ROI, put it in the label vector
              var peak = roi.peaks[i],
                  label = self.drawLabel(peak, path);
              if (label)
                labels.push(label);
            }

          if (labels.length > 0 && !labelAlreadyAdded)                                  // The first element inside the labels vector is the current peak's label
            self.peakLabelArray.push( { "path": path, "label": labels.shift() } );

          self.leftoverPeakLabels = self.leftoverPeakLabels.concat(labels);         // Any leftover peaks are put inside this array

          // Highlights a specified peak (darken the context)
          function highlightPeak(d, peak) {
            if (self.zooming_plot)
              return;

            if (self.highlightedPeak)
              unhighlightPeak(0, self.highlightedPeak);

            var thePeak = this;
            if (d3.event == null || d3.event.touches) {
              thePeak = peak;
              if (!thePeak || thePeak.attr("fill-opacity") != 0.6)
                return;
            }

            if (d3.select(thePeak).attr("fill-opacity") != 0.6)
              return;

            if (d3.event.buttons === 1 || d3.event.buttons === 2)
              return;

            if (!peak) {
              if (paths.length === 2) 
                d3.select(thePeak).attr("fill-opacity",0.8);
              else
                d3.select(thePeak).attr("stroke-width",2).attr("fill-opacity",0.8);

            } else {
              if (paths.length === 2) 
                peak.attr("fill-opacity",0.8);
              else
                peak.attr("stroke-width",2).attr("fill-opacity",0.8);
            }

            if (d3.select(thePeak)[0].parentNode) {    // if the 'this' pointer is the path for the peak, then declare the highlighted peak to be that
              self.highlightedPeak = thePeak;
            } else if (peak)
              self.highlightedPeak = peak;

            for (i = 0; i < self.peakLabelArray.length; i++)
              if (self.peakLabelArray[i].path === (peak ? peak : path)) {
                self.peakLabelArray[i].label.attr('stroke', 'black')
                  .attr("z-index", 100);
                self.highlightedLabel = self.peakLabelArray[i].label;
              }
          }

          // Unhighlights a specified peak (returns to default color)
          function unhighlightPeak(d, highlightedPeak) {
            if (self.highlightedLabel) {
              self.highlightedLabel.attr('stroke', 'none')
                    .attr("z-index", 0);
            }

            if (!highlightedPeak && !self.highlightedPeak) {
              self.highlightedPeak = null;
              self.highlightedLabel = null;
              return;
            }

            var peak = self.highlightedPeak;
            if (!highlightedPeak) {
              if (!Array.isArray(self.highlightedPeak)) 
                peak = d3.select(self.highlightedPeak);
              if (paths.length === 2)  peak.attr("fill-opacity",0.6);
              else                     peak.attr("stroke-width",1).attr("fill-opacity",0.6);
            } else {
              if (!Array.isArray(highlightedPeak)) 
                highlightedPeak = d3.select(highlightedPeak);
              if (paths.length === 2)  highlightedPeak.attr("fill-opacity",0.6);
              else                     highlightedPeak.attr("stroke-width",1).attr("fill-opacity",0.6);
            }


            self.highlightedPeak = null;
            self.highlightedLabel = null;
          }

          self.unhighlightPeakFunc = unhighlightPeak;
          self.highlightPeakFunc = highlightPeak; 

          function onRightClickOnPeak() {
            console.log("Emit RIGHT CLICK (ON PEAK) signal. (Peak roi = ", roi, ")");
          }

          self.pathEnergies.push( {"path": path, "lower_energy": roi.lowerEnergy, "upper_energy": roi.upperEnergy} );


          if( num === 0 ){
            path.attr("class", "roi")
              .style("fill",function(d){ return "blue"; } )
              .attr("fill-opacity",0.6);

            if( paths.length === 2 ){
              path.attr("stroke-width",1)
                .attr("stroke", "blue" )
                .on("mouseover", highlightPeak )
                .on("touchend", highlightPeak)
                .on("mouseout",  unhighlightPeak);
            }

          } else {
            path// .attr("class", "peak")
              .style("fill","blue" )
              .attr("fill-opacity",0.0)
              .attr("stroke-width",1)
              .attr("stroke", "blue" );

          }

          // For right-clicking on a peak
          path.on("contextmenu", onRightClickOnPeak);
        });
  };

  for (var i = 0; i < this.rawData.spectra.length; ++i) {
    var spectrumi = i;
    if (!this.rawData.spectra[i] || !this.rawData.spectra[i].peaks || !this.rawData.spectra[i].peaks.length || 
      (typeof self.options.drawPeaks[i] !== 'undefined' && !self.options.drawPeaks[i])) {
      continue;
    }

    this.rawData.spectra[i].peaks.forEach( function(roi){
      self.pathEnergies = [];
      draw_roi(roi,i,self.rawData.spectra[spectrumi]);
      
      //console.log( 'There are ' + (peaks ? peaks.length : 'null') + " peaks" );
    });
    
    if( self.current_fitting_peak )
      draw_roi( self.current_fitting_peak, 0, self.rawData.spectra[0] /* foreground */ );
  }
}

SpectrumChartD3.prototype.drawLabel = function(peak,path) {
  var self = this;

  if ( self.options.showUserLabels || self.options.showPeakLabels || self.options.showNuclideNames ) {
      var pathNode = path.node();                 // node for the peak element
      self.solver = new c.SimplexSolver();        // Initialize a default, new solver for each label creation

      var chartBox = d3.select("#chartarea"+self.chart.id)[0][0].getBBox();    // box coordinates for the chart area
      var chart = self.peakVis

      /*
      The path box coordinates are the box coordinates for the peak. These are used for positioning the label properly
      with respect to the peak.
      */
      var pathX = pathNode.getBBox().x,
          pathY = pathNode.getBBox().y,
          pathBoxHeight = pathNode.getBBox().height,
          pathBoxWidth = pathNode.getBBox().width,
          pathBoxMiddle = (pathX + pathBoxWidth) / 2;

      var peakEnergy = peak.Centroid[0].toFixed(2) + " keV";
      var nuclide = peak.nuclide;
      var userLabelText = peak.userLabel;

      var label, userLabel, peakEnergyLabel, nuclideNameLabel;

      // Main label DOM
      label = chart.append("text");
      label.attr('class', 'label')
        .attr("text-anchor", "start")
        .attr("x", pathX)
        .attr("y", pathY - 10 )
        .attr("energy", peakEnergy);

      if (self.options.showUserLabels) {                      // T-span element for the user label
        if (userLabelText) {                                  // Gets priority for top-most text in label
          userLabel = label.append("tspan");
          userLabel.attr("class", "userLabel")
            .attr("x", pathX)
            .attr("dy", 0)
            .text(userLabelText)
            .style('font-size', '7.5px');
        }
      }

      if (self.options.showPeakLabels) {                      // T-span element for the peak label
          peakEnergyLabel = label.append("tspan");            // Gets second in priority for top-most text in label
          peakEnergyLabel.attr("class", "peakEnergyLabel")
            .attr("x", pathX)
            .attr("dy", userLabel ?  "1em" : 0)              // If user label is not present, then this is the top-most element
            .text(peakEnergy)
            .style('font-size', '7.5px');

          // if (userLabel)
          //   peakEnergyLabel.attr("dx", -userLabel.node().getBBox().width/2);
      }

      if (self.options.showNuclideNames && nuclide) {         // T-span element for nuclide label
        nuclideNameLabel = label.append("tspan");             // Third in priority for top-most text in label
        nuclideNameLabel.attr("class", "nuclideLabel")
          .attr("x", pathX)
          .attr("dy", self.options.showUserLabels && self.options.showPeakLabels ? "1em" : self.options.showUserLabels || self.options.showPeakLabels ?  "1em" : 0)
          // .attr("dx", -label.node().getBBox().width/2)               
          .text(nuclide.name)
          .style('font-size', '7.5px');
      }

      if ( !userLabel && !peakEnergyLabel && !nuclideNameLabel) {      // Do not display an label if the user label, peak label, or nuclide name label are not displayed
        label.remove();                                                     // This means that although nuclide energy label option is selected, it is not displayed!
        return;
      }

      if (self.options.showNuclideEnergies && nuclide) {                    // Nuclide energy label displayed only if nuclide name labels are displayed!
        if (self.options.showNuclideNames)
          nuclideNameLabel.text(nuclide.name + ", " + nuclide.DecayGammaEnergy.toString() + " keV" );
      }

      // Inequality constants for constraints
      var greaterThan = c.GEQ,
          lessThan = c.LEQ;

      /*
      Strengths of constraints: These determine the prioirity of which constraints are given for a specific label.
      A required constraint means it must be followed, or else a Constraint Error is thrown.
      A strong constraint isn't required, but has priority over weaker constraints. (and so on for the weaker consraints.
      */
      var weak = c.Strength.weak,
          medium = c.Strength.medium,
          strong = c.Strength.strong,
          required = c.Strength.required;

      /*
      Create variable values for the Constraint Solver to handle. In particular, for centering each label at 
      the center of the peak, we want to capture the left-x value, the mid-point x-value, the width of the label,
      and the middle x-value for the path box of a peak.

      Initially, we provided the starting coordiantes for the labels and initialized their values in their variable creation.
      What we want to get in the end of adding all the constraints are different values for these positions so that the label placement
      is exactly where we want it to be.
      */
      var labelLeftX = new c.Variable ( { value : pathX } ),                                                       // Left coordinate for label
          labelRightY = new c.Variable ( { value : pathX + label.node().getBBox().width } ),                       // Right coordinate for label
          labelTopY = new c.Variable ( { value : label.node().getBBox().y } ),                                     // Top coordinate for label
          labelBottomY = new c.Variable ( { value : label.node().getBBox().y + label.node().getBBox().height } ),  // Bottom coordinate for label
          labelMid = new c.Variable( { value : (pathX + (label.node().getBBox().width) / 2) } ),                   // Mid-point coordinate for label (in terms of x)
          labelHeight = new c.Variable( { value : label.node().getBBox().height } ),                               // Height of label
          labelWidth = new c.Variable( { value : label.node().getBBox().width } ),                                 // Width of label
          peakMid = new c.Variable( { value : pathX + (pathBoxWidth / 2) } );                                      // Mid-point for the box for the peak element

      var cle;  // These are reusable variables for creating equations for constraints.

      /*
      Christian:
      The mid-point x-value of a peak and the width of the label will be staying constant throughout each added constraint.
      Add these constraints into the Constraint Solver.
      */
      self.solver.addConstraint(new c.StayConstraint( peakMid, required  ) );           // mid-point of the box for the peak element stays constant (you're not moving the peak!)
      self.solver.addConstraint(new c.StayConstraint( labelHeight, required  ) );       // keep height constant
      self.solver.addConstraint(new c.StayConstraint( labelWidth, required  ) );        // keep width constant
      self.solver.addConstraint(new c.StayConstraint( labelTopY, weak ) );             // try to keep initial y-position constant (right above the peak)



      /*
      Label assertions for the mid-point, width, and height values for a label. 
      These keep the internal size of a label consistent of what it currently is.
      */
      cle = c.plus(labelLeftX, c.divide(labelWidth, 2));                   
      var labelMidPointAssertion = new c.Equation( labelMid, cle );      // Left-x + (width/2) == mid-point

      cle = c.plus(labelLeftX, labelWidth);
      var labelWidthAssertion = new c.Equation( labelRightY, cle );      // left-x + width == right-x

      cle = c.plus(labelTopY, labelHeight);     
      var labelHeightAssertion = new c.Equation( labelBottomY, cle );    // top-y + height == bottom-y

      // var original_area = labelHeight.value * labelWidth.values

      self.solver.addConstraint( labelMidPointAssertion );           
      self.solver.addConstraint( labelWidthAssertion );
      self.solver.addConstraint( labelHeightAssertion );



      /*
      Label position assertions. These are added again to keep the size of the label consistent and accurate for collision detection with other elements.
      For example, the bottom coordinate of a label is always equal to the top_coordinate + height.
      Although many of these are obvious, it is important for consistency in collision detection.
      */
      self.solver.addConstraint( new c.Equation(labelTopY, c.minus(labelBottomY, labelHeight) ) );         // labelTopY = labelBottomY - label height
      self.solver.addConstraint( new c.Equation(labelBottomY, c.plus(labelTopY, labelHeight) ) );          // label bottom y = label top y + label height
      self.solver.addConstraint( new c.Equation(labelLeftX, c.minus(labelRightY, labelWidth) ) );          // label left x = label right x - label width
      self.solver.addConstraint( new c.Equation(labelRightY, c.plus(labelLeftX, labelWidth) ) );           // label right x = label left x + label width

      self.solver.addConstraint( new c.Equation(labelHeight, c.minus(labelBottomY, labelTopY) ) );         // label height = labelBottomY - labelTopY
      self.solver.addConstraint( new c.Equation(labelWidth, c.minus(labelRightY, labelLeftX) ) );          // label width = labelRightY - labelLeftX  
      self.solver.addConstraint( new c.Equation(labelMid, c.plus(labelLeftX, c.divide(labelWidth, 2) ) ) ); // label mid-point = lable_leftX + (labelWidth/2)

      self.solver.addConstraint( new c.Inequality( labelMid, greaterThan, labelLeftX ) );     // mid-point coordinate > left coordinate
      self.solver.addConstraint( new c.Inequality( labelRightY, greaterThan, labelLeftX) );  // right coordinate > left coordinate
      self.solver.addConstraint( new c.Inequality( labelBottomY, greaterThan, labelTopY ) ); // bottom coordinate > top coordinate

      // Label does not fall out of bounds - constraints
      self.solver.addConstraint( new c.Inequality(labelLeftX, greaterThan, 0) );                      // Left coordinate does not fall out of bounds to the left
      self.solver.addConstraint( new c.Inequality(labelTopY, greaterThan, self.padding.top*2) );      // Top coordinate does not fall out of bounds from the top
      self.solver.addConstraint( new c.Inequality( labelRightY, lessThan, chartBox.width ) );        // Right coordinate does not exceed right bounds of chart
      self.solver.addConstraint( new c.Inequality( labelBottomY, lessThan, chartBox.height ) );      // Bottom coordinate does not exceed bottom bounds of chart
      self.solver.addConstraint( new c.Inequality(labelBottomY, greaterThan, self.padding.top*2) );   // Bottom coordinates does not fall out of bounds from the top

      // self.solver.addConstraint( new c.Inequality( labelBottomY, lessThan, pathY, weak) );   // bottom y for label is above peak box


      /*
      To align the label properly with the peak, we try to keep the mid-point x-coordinate of a label
      aligned with the middle of the box for a peak element. To do this, we try to set the mid-point value
      for the label equal to the mid-point value for a peak's box.
      */
      var labelMidAtPeakMid = new c.Equation(labelMid, peakMid, strong, 5)
      self.solver.addConstraint( labelMidAtPeakMid );


      var numberOfSamePeaks = 1,                       // keeps track of the number of the same peak labels
          ypadding = 10,                                  // y-padding between labels
          xpadding = 20;                                  // x-padding between labels

      /* Returns true if another label overlaps the current label from the top. */
      function overlapFromTop(topY, otherLabel) {
        otherLabelSelect = d3.select(otherLabel);

        var otherLabelTopY = Number( ( otherLabelSelect.attr('y') == null ? otherLabel.getBBox().y : otherLabelSelect.attr('y') ) );
        var otherLabelBottomY = otherLabelTopY + otherLabel.getBBox().height;
        return otherLabelTopY < topY && otherLabelBottomY >= topY;
      }
      /* Returns true if another label overlaps the current label from the bottom. */
      function overlapFromBottom(topY, bottomY, otherLabel) {
        otherLabelSelect = d3.select(otherLabel);
        var otherLabelTopY = Number( ( otherLabelSelect.attr('y') == null ? otherLabel.getBBox().y : otherLabelSelect.attr('y') ) );
        return otherLabelTopY > topY && otherLabelTopY <= bottomY;
      }
      /* Returns true if another label overlaps the current label from the left side. */
      function overlapFromLeftSide(leftX, rightX, otherLabel){
        otherLabelSelect = d3.select(otherLabel);
        var otherLabelLeftX =  Number( ( otherLabelSelect.attr('x') == null ? otherLabel.getBBox().x : otherLabelSelect.attr('x') ) );
        var otherLabelRightY = otherLabelLeftX + otherLabel.getBBox().width;
        return otherLabelRightY >= leftX && otherLabelRightY <= rightX;
      }
      /* Returns true if another label overlaps the current label from the right side. */
      function overlapFromRightSide(leftX, rightX, otherLabel){
        otherLabelSelect = d3.select(otherLabel);
        var otherLabelLeftX =  Number( ( otherLabelSelect.attr('x') == null ? otherLabel.getBBox().x : otherLabelSelect.attr('x') ) );
        var otherLabelRightY = otherLabelLeftX + otherLabel.getBBox().width;
        return otherLabelLeftX <= rightX && otherLabelLeftX >= leftX;
      }
      /* Returns true if another label overlaps the current label directly. (Same left and top coordinates) */
      function overlapDirectly(leftX, topY, otherLabel){
        otherLabelSelect = d3.select(otherLabel);
        var otherLabelLeftX =  Number( ( otherLabelSelect.attr('x') == null ? otherLabel.getBBox().x : otherLabelSelect.attr('x') ) );
        var otherLabelTopY = Number( ( otherLabelSelect.attr('y') == null ? otherLabel.getBBox().y : otherLabelSelect.attr('y') ) );
        return otherLabelTopY == topY && otherLabelLeftX == leftX;
      }
      /* Returns true if another label overlaps the current label from anywhere. */
      function overlap( leftX, rightX, topY, bottomY, otherLabel  ) {
        var overlapFromSide = overlapFromLeftSide(leftX, rightX, otherLabel) || overlapFromRightSide(leftX, rightX, otherLabel);

        return (overlapFromTop(topY, otherLabel) && overlapFromSide) || 
               (overlapFromBottom(topY, bottomY, otherLabel) && overlapFromSide) ||
               (overlapFromSide && topY == Number(otherLabelSelect.attr('y'))) ||
               overlapDirectly(leftX, topY, otherLabel);
      }
      /* Adds a required inequality to the constraint solver. If an error is thrown and te required inequality cannot be added,
         then the inequality is changed to have a "strong" strength 
      */
      function addRequiredInequality( variable, inequality, value, otherLabel ) {
        try { self.solver.addConstraint( new c.Inequality( variable, inequality, value, required ), required ); }
        catch (e) { 
          console.log("Failed to make a constraint 'required', so making it 'strong'. Constraint was for ", peakEnergy, " on ", otherLabel.textContent);
          self.solver.addConstraint( new c.Inequality( variable, inequality, value, strong, 100 ), required ); 
        }
      }



      /*
      Main function called for fixing the overlaps for a label. To sum up how this function works:

        1. For each of the previous peaks drawn before this current label being drawn

          a.  If that box for a peak overlaps the current label from the side (but have the same top-coordinate values)
            i)  Move the current label to the right and down

          b.  Else if the box for a peak overlaps the current label from the top and to the side
            i)  Move the current label down

          c.  Else if the box for a peak overlaps the current label from the bottom and to the side
            i)  Move the current label up.

        2. For each of the previous labels drawn before this current label being drawn

          a.  If that other label overlaps the current label directly (same left and top coordiantes)
            i)  Move the current label down

          b.  Else if the other label overlaps the current label from the side, but the top coordinate values are the same
            i)  Move the current label to the right
            ii) Move the current label down

          c.  Else if the other label overlaps the current label from the side and the top
            i)  Move the current label down

          d.  Else if the other label overlaps the current label from the bottom and the side
            i)  Move the current label up
      */
      function fixOverlaps() {
        var overlappedLabels = [];

        // For overlapping peaks - not yet finsihed!
        self.peakVis.selectAll("path").each(function(d, i) {
          var peak = d3.select(this),
              peakBox = peak.node().getBBox();

          var peakLeftX =  peakBox.x,
              peakRightX = peakLeftX + peakBox.width,
              peakTopY = peakBox.y,
              peakBottomY = peakTopY + peakBox.height;

          var peakOverlapsFromTop = overlapFromTop(labelTopY.value, peak.node()),
              peakOverlapsFromBottom = overlapFromBottom(labelTopY.value, labelBottomY.value, peak.node()),
              peakOverlapsFromSide = overlapFromLeftSide( labelLeftX.value, labelRightY.value, peak.node() ) || 
                                        overlapFromRightSide(labelLeftX.value, labelRightY.value, peak.node()),
              overlapping = overlap(labelLeftX.value, labelRightY.value, labelTopY.value, labelBottomY.value, peak.node());

          if (overlapping) {
            // console.log(peakEnergy, " overlaps the peak ", peak.node());

            if (peakOverlapsFromSide && labelTopY.value == peakTopY) {
              // console.log(otherLabel.textContent, " is overlapping ", text, " from side, so moving ", text, " down.");
              addRequiredInequality( labelLeftX, greaterThan, peakRightX+xpadding, peak );                         // move label to the right
              addRequiredInequality( labelTopY, greaterThan, peakBottomY+ypadding, peak );                         // move label down
            }
            else if (peakOverlapsFromTop && peakOverlapsFromSide) {
              // console.log(otherLabel.textContent, " is overlapping ", text, " from top, so moving ", text, " down.");
              addRequiredInequality( labelTopY, greaterThan, peakBottomY+ypadding, peak );                         // move label down
            }
            else if (peakOverlapsFromBottom && peakOverlapsFromSide) {
              // console.log(otherLabel.textContent, " is overlapping ", text, " from bottom, so moving ", text, " up.");
              addRequiredInequality( labelBottomY, lessThan, peakTopY-ypadding, peak );                            // move label up
            }
          }
        });

        self.peakVis.selectAll("text").each(function(d, i){
          otherLabel = d3.select(this).node();

          if (otherLabel != label.node() && otherLabel.textContent == label.node().textContent) {                           // delete duplicate labels

                                                                                                                          // removel label if overlapping directly and same text content
            if ( !(self.options.showNuclideNames && !self.options.showPeakLabels && !self.options.showNuclideEnergies) || overlapDirectly( labelLeftX.value, labelTopY.value, otherLabel) ) {
              otherLabel.remove();                // If showing only nuclide name, don't delete duplicates

            }
            ++numberOfSamePeaks;
          }
          else if ( otherLabel != label.node() ) {
            otherLabelSelect = d3.select(otherLabel);

            var otherLabelLeftX =  Number(otherLabelSelect.attr('x')),
                otherLabelRightY = otherLabelLeftX + otherLabel.getBBox().width,
                otherLabelTopY = Number(otherLabelSelect.attr('y')),
                otherLabelBottomY = otherLabelTopY + otherLabel.getBBox().height;

            // Get booleans for overlap from other label
            var otherLabelOverlapsFromTop = overlapFromTop(labelTopY.value, otherLabel),
                otherLabelOverlapsFromBottom = overlapFromBottom(labelTopY.value, labelBottomY.value, otherLabel),
                otherLabelOverlapsFromSide = overlapFromLeftSide( labelLeftX.value, labelRightY.value, otherLabel ) || overlapFromRightSide(labelLeftX.value, labelRightY.value, otherLabel),
                directOverlap = overlapDirectly( labelLeftX.value, labelTopY.value, otherLabel),
                overlapping = overlap(labelLeftX.value, labelRightY.value, labelTopY.value, labelBottomY.value, otherLabel);

            if ( overlapping ) {
              overlappedLabels.push( otherLabel );

              if (directOverlap) {
                addRequiredInequality( labelTopY, greaterThan, otherLabelBottomY+ypadding, otherLabel );                         // move label down
              }
              else if (otherLabelOverlapsFromSide && labelTopY.value == otherLabelTopY) {
                // console.log(otherLabel.textContent, " is overlapping ", text, " from side, so moving ", text, " down.");
                addRequiredInequality( labelLeftX, greaterThan, otherLabelRightY+xpadding, otherLabel );                         // move label to the right
                addRequiredInequality( labelTopY, greaterThan, otherLabelBottomY+ypadding, otherLabel );                         // move label down
              }
              else if (otherLabelOverlapsFromTop && otherLabelOverlapsFromSide) {
                // console.log(otherLabel.textContent, " is overlapping ", text, " from top, so moving ", text, " down.");
                addRequiredInequality( labelTopY, greaterThan, otherLabelBottomY+ypadding, otherLabel );                         // move label down
              }
              else if (otherLabelOverlapsFromBottom && otherLabelOverlapsFromSide) {
                // console.log(otherLabel.textContent, " is overlapping ", text, " from bottom, so moving ", text, " up.");
                addRequiredInequality( labelBottomY, lessThan, otherLabelTopY-ypadding, otherLabel );                            // move label up
              }
              // Do something with the font size here
              /* I tried increasing the font size as much as it can without hitting the other labels, but it seems that it keeps the font at 1vw here. */
              // self.solver.addConstraint( new c.Inequality( labelHeight, greaterThan, labelHeight.value + 2, strong ) );       // label height > 1
              // self.solver.addConstraint( new c.Inequality( labelWidth, greaterThan, labelWidth.value + 2, strong ) ); 

            }
          }
          return overlappedLabels;
        });
      }

      var overlappedLabels = fixOverlaps();
      // fixOverlaps();          // can uncomment to run fixOverlaps twice, which may optimize how the labels overlap, with the sacrifice for performance

      /*
      Return label back to its original style.
      This means:
        - Label gets unbolded (if it was already bolded)
        - Label z-index goes back to default instead of being shown at the very top.
        - If an arrow/line is shown from the label to its corresponding peak, then that line/arrow is deleted.
      */
      function normalizeLabel(labelToNormalize) {                     // Return label back to original style on mouse-out.

            d3.select(labelToNormalize ? labelToNormalize : this)
              .style('z-index', 0)
              .style('cursor', 'default')
              .attr('stroke', 'none')
              .attr('font-weight', 'normal');

            // delete the pointer line from the label to the peak
            if (self.peakLabelLine) {
              self.peakLabelLine.remove();
              self.peakLabelLine = null;
            }
        }

      /*
      Highlight a selected label.
      This means:
        - Label becomes bold
        - Label's z-index goes to very top (so that the whole text is shown)
        - A line is drawn from the label to its corresponding peak
      */
      function highlightLabel() {
            if ( !self.dragging_plot ) {                          // Bold the label text and add a line (arrow) that points to the peak when moused over text.
              if ( self.labelToNormalize )
                normalizeLabel(self.labelToNormalize);

              d3.select(this)
                .style('cursor', 'default')
                .attr('stroke', 'black')
                .attr('font-weight', 'bold')
                .attr("z-index", 100);

              var x1, x2, y1, y2;

              if ( labelLeftX.value > peakMid.value ) { x1 = labelLeftX.value; }
              else if ( labelRightY.value < peakMid.value ) {  x1 = labelRightY.value; }
              else { x1 = labelMid.value; }

              // label heights are sometimes 1, mostly 18 --> look into the constraints why
              if ( labelTopY.value > pathY && labelBottomY.value > pathY ) {  y1 = labelTopY.value+(labelHeight.value - 10);  }    // if label is under peak
              else {  y1 = labelBottomY.value-labelHeight.value + (self.options.showNuclideNames && self.options.showNuclideEnergies ? 10 : 0);  }
                                                                      // To give some space between label and line

              x2 = pathX + pathBoxWidth/2; 
              y2 = pathY-2;

               // Here I am trying to draw an arrow marker for the line from the label to a peak
               if (!self.peakLabelArrow)
                self.peakLabelArrow = self.peakVis.append('svg:defs').append("svg:marker")
                              .attr("id", "triangle")
                              .attr('class', 'peaklabelarrow')
                              .attr('viewbox', "0 -5 10 10")
                              .attr("refX", 2.5)
                              .attr("refY", 2.5)
                              .attr("markerWidth", 30)
                              .attr("markerHeight", 30)
                              .attr("orient", "auto")
                              .append("path")
                                .attr("d", "M 0 0 20 6 0 12 3 6")
                                .attr("transform", "scale(0.4,0.4)")
                                .style("stroke", "black");

              self.peakLabelLine = self.peakVis.append('line')
                            .attr('class', 'peaklabelarrow')
                            .attr('x1', x1)
                            .attr('y1', y1)
                            .attr('x2', x2)
                            .attr('y2', y2)
                            .attr("marker-end", "url(#triangle)");


            } 

            self.labelToNormalize = this;
          }

      label
        .attr("x", labelLeftX.value)
        .attr("y", labelTopY.value)

        // Set the font-size here
        // .style("font-size", ((labelWidth.value * labelHeight.value) / original_area) + "vw")

        .on("mouseover", highlightLabel)
        .on("touchstart", highlightLabel)
        .on("mouseout",  normalizeLabel);

      if (userLabel)
        userLabel.attr("x", labelLeftX.value);
      if (peakEnergyLabel)
        peakEnergyLabel.attr("x", labelLeftX.value);
      if (nuclideNameLabel)
        nuclideNameLabel.attr("x", labelLeftX.value);
  }

    return label;

      /* Possible algorithm: Get list of all text nodes. Have a new array of grouped text nodes by overlapped text labels. This array will hold an array 
         of grouped nodes, and once this array has been created and filled, start doing stuff with it. Specifically, replace all those labels with something like
         "5 peaks from 100-200 kEV". Individual nodes will just be outputted normally ("500 kEV"). However, user may lose some information when zoomed all the
         way out of the graph. They would have to zoom-in further to see the peaks.
      */  
}



SpectrumChartD3.prototype.handleCancelAllMouseEvents = function() {
  var self = this;

  return function () {

    d3.select(document.body).style("cursor", "default");
    self.downx = Math.NaN;
    self.downy = Math.NaN;

    // Cancel recalibration
    self.handleCancelMouseRecalibration();

    // Cancel deleting peaks
    self.handleCancelMouseDeletePeak();

    // Cancel zooming in x-axis
    self.handleCancelMouseZoomInX();

    // Cancel zooming in y-axis
    self.handleCancelMouseZoomInY();

    // Cancel counting gammas
    self.handleCancelMouseCountGammas();

    // Canceling out all mouse events for chart
    self.zooming_plot = false;
    self.dragging_plot = false;
    self.recalibrationStartEnergy = null;

    // Cancel all kept mouse positions for mouse events
    self.leftMouseDown = null;
    self.zoominmouse = self.deletePeaksMouse = self.countGammasMouse = self.recalibrationMousePos = null;
    self.rightClickDown = null;

    // Cancel all legend interactions
    self.legdown = null;

    // Cancel scaler interactions
    self.scalerdown = null;

    // Cancel all slider drag region interactions
    self.sliderBoxDown = false;
    self.leftDragRegionDown = false;
    self.rightDragRegionDown = false;

    // Cancel all scaler widget interactions
    self.currentlyAdjustingSpectrumScale = null;
  }
}


//Liza: This function gets called whenever you press and hold the option key 
//      (alt on non-macs) then click and hold down with the left mouse button, 
//      and drag to the right.   
SpectrumChartD3.prototype.handleMouseMovePeakFit = function() {
  var self = this;

  console.log( "In handleMouseMovePeakFit + " + d3.mouse(self.vis[0][0])[0] );

  self.peakFitMouseMove = d3.mouse(self.vis[0][0]);
  self.peakFitTouchMove = d3.touches(self.vis[0][0]);

  var leftTouch, rightTouch;
  if (self.peakFitTouchMove.length > 0) {
    leftTouch = self.peakFitTouchMove[0][0] < self.peakFitTouchMove[1][0] ? self.peakFitTouchMove[0] : self.peakFitTouchMove[1];
    rightTouch = leftTouch === self.peakFitTouchMove[0] ? self.peakFitTouchMove[1] : self.peakFitTouchMove[0];
  }

  self.zooming_plot = false;

  // We don't want to redraw the reference lines if we are using touches
  if (self.peakFitTouchMove.length == 0)
    self.drawPeakFitReferenceLines();

  if (!self.zoominx0 && self.peakFitTouchMove.length > 0)
    self.zoominx0 = self.xScale.invert(leftTouch[0]);

  // Set the original X-domain if it does not exist
  if (!self.origdomain)
    self.origdomain = self.xScale.domain();
  
  //This next line is a hack to get around how D3 is managing the zoom, but we highkacked it for the peak fitting 
  self.xScale.domain( self.origdomain );
  
  //Get the initial energy the user 'clicked-down on'
  var x0_energy = self.zoominx0;
  // var y0_energy = self.zoominy0;
  
  //Get the current energy the mouse is at
  var x_energy = self.peakFitTouchMove.length == 0 ? self.xScale.invert(d3.mouse(self.vis[0][0])[0]) : self.xScale.invert(rightTouch[0]);
  //var y_counts = self.xScale.invert(mouse_pos_px[1]);
  
  console.log( "x0_energy=" + x0_energy + ", x_energy=" + x_energy );
  
  //Make sure x0_energy is less than x_energy
  if( x0_energy > x_energy )
    x_energy = [x0_energy, x0_energy = x_energy][0];
  
  var need_redraw_on_cancel = !!self.current_fitting_peak;
  self.current_fitting_peak = null;
  
  //If we arent displaying a spectrum, remove the current_fitting_peak and return 
  //  TODO: Have peak fitting for other spectra??
  if( !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length
      || self.rawData.spectra[0].y.length == 0 || this.rawData.spectra[0].y.length < 10 ) {
    if( need_redraw_on_cancel )     
      self.redraw()();
    return;
  }
  
  //Get rid of currently fit peak 
  self.current_fitting_peak = null;
  
  //Get the bins for the range we're fitting
  var foreground = self.rawData.spectra[0];
  var bisector = d3.bisector(function(d){return d;}); 
  var lbin = bisector.left(foreground.x, x0_energy);
  var rbin = bisector.right(foreground.x, x_energy);
  
  if( rbin < (foreground.x.length - 2) )
    rbin += 1;
  
  //Make sure we have at least a few bins  
  if( rbin < (lbin + 4) ){
    console.log( 'Not enough bins yet: ' + (rbin - lbin) );
    if( need_redraw_on_cancel )     
      self.redraw()();
    return;
  }
                  
  //Now get the actual energy range we are fitting (the lowere energy of left 
  //  bin, and upper energy of right bin)
  var lbin_lower_energy = foreground.x[lbin];
  var rbin_upper_energy = foreground.x[rbin + 1];  //FIXME: +1 will cause problem for last bin...
  
  //The polynomial continuum should be evalueated relative to a reference energy.
  //  Lets pick the middle of the ROI for this.
  var roi_ref_energy = 0.5*(rbin_upper_energy + lbin_lower_energy); 
  
  //Get the actual counts in the bins.  Note that energies is one longer than
  //  counts, so you can can access the end of the fitting range.
  var energies = foreground.x.slice(lbin,rbin+1); //FIXME: +1 dangerous on last bin!
  var counts = foreground.y.slice(lbin,rbin);
  
  //Here is where you would fit for your peaks.  As an example, I make a linear
  //  continuum based on first/last bin height (to provide an example of how 
  //  InterSpec defines the continuum as a density per keV, hence you have to 
  //  integrate over a bins width to get the bins continuum height), and then
  //  add a peak with an amplitude equal to the difference between the continuum
  //  area, and the actual data area. 

  
  //A function to integrate the area bellow the continuum between x0 and x1 
  //  for a given refrence energy
  function offset_eqn_integral(coeffs,refEnergy,x0,x1) {
    x0 -= refEnergy;
    x1 -= refEnergy;
    var answer = 0.0;
    for( var order = 0; order < coeffs.length; ++order ) {
      var exp = order + 1.0;
      answer += (coeffs[order]/exp) * (Math.pow(x1,exp) - Math.pow(x0,exp));
    }
    return answer;
  };
  
  //A function to get linear continuum coeficients using the hights of the first
  //  and last bin - this is only an example function!  do not use in real code,
  //  it can fail in many ways.
  //Function assumes xvals has one more element than the range to be fitted, and 
  //  yvals has exactly the number of bins to be fitted
  function linear_eqn_from_first_last_data( xvals, yvals, refenergy ) {
    var x1 = xvals[0];
    var dx1 = xvals[1] - x1;
    
    var x2 = xvals[xvals.length-2];
    var dx2 = xvals[xvals.length-1] - x2;
    
    x1 -= refenergy;
    x2 -= refenergy;
  
    var y1 = yvals[0];
    var y2 = yvals[yvals.length-1];

    var c = (2.0*x2*dx2 + dx2*dx2)/(2.0*x1*dx1 + dx1*dx1);
    var b = (y2-y1*c)/(dx2-dx1*c);
    var m = 2.0*(y1-b*dx1)/(2.0*x1*dx1+dx1*dx1);
    
    return [b,m];
  };

  function gaus_integral(mean, width, amp, x0, x1) {
      var peak_mean = mean, peak_sigma = width, peak_amplitude = amp;
         //peak.LandauAmplitude, peak.LandauMode, peak.LandauSigma,
      var sqrt2 = 1.414213562373095;
      var erflowarg = (x0-peak_mean)/(sqrt2*peak_sigma);
      var erfhigharg = (x1-peak_mean)/(sqrt2*peak_sigma);
      return peak_amplitude * 0.5 * (erf(erfhigharg) - erf(erflowarg));
  };

  function erf(x) {
      //http://stackoverflow.com/questions/14846767/std-normal-cdf-normal-cdf-or-error-function
      var sign = (x >= 0) ? 1 : -1; // save the sign of x
      x = Math.abs(x);
      var a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027
          a5 = 1.061405429, p  = 0.3275911;
      var t = 1.0/(1.0 + p*x);
      var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
      return sign * y; // erf(-x) = -erf(x);
  };

  var fit = peakFit.run(numeric.sub(energies, roi_ref_energy), counts);
  var params = fit.params;
  var polyd = fit.polyd;
  var peakchi2 = fit.chi2;
  var numPeaks = (params.length - (polyd + 1)) / 3;
  var sumdata = counts.reduce(function(a,b){return a+b;}, 0);
  var continuumCoefs = params.slice(0,polyd+1);
  var sumcontinuum = offset_eqn_integral(continuumCoefs,roi_ref_energy,lbin_lower_energy,rbin_upper_energy);

  var poly_type;
  switch (polyd) {
    case 0:
      poly_type = 'Constant';
    case 1:
      poly_type = 'Linear';
    case 2:
      poly_type = 'Quadratic';
    case 3:
    poly_type = 'Cubic';
  }

  self.current_fitting_peak =  { 
      type: poly_type,                       //options for you are are: 'NoOffset', 'Constant', 'Linear', 'Quardratic', 'Cubic'
      lowerEnergy: lbin_lower_energy,
      upperEnergy: rbin_upper_energy,
      referenceEnergy: roi_ref_energy,
      coeffs: continuumCoefs,  //Linear should have 2 coeficients, Quadratic 3, Cubic 4, constant 1
      coeffUncerts:[0,0],                  //you can put in actual values of uncertainties here
      fitForCoeff: Array.apply(null, Array(polyd + 1)).map(function(){return true}), //Same number of entries as coeffs; each entry should be true for now 
      peaks:[]
  };

  for (var i=polyd+1; i<params.length; i+=3) {
    var peakamp = Math.abs(sumdata - sumcontinuum);
    var gausmean = params[i+1] + roi_ref_energy;
    var peaksigma = Math.abs(params[i+2]);
    var gaus_sum = 0;
    for (var j=0; j<counts.length; j++) {
      var x0 = energies[j];
      var x1 = energies[j+1];
      gaus_sum += gaus_integral(gausmean, peaksigma, peakamp, x0, x1);
    }
    peakamp = Math.abs(gaus_sum - sumcontinuum);
    self.current_fitting_peak.peaks.push({
        type:'GaussianDefined',
        skewType:'NoSkew',
        Centroid:[gausmean,0,true],        //second element is uncertainty, third alsways true
        Width:[peaksigma,0.00431738,true], //second element is uncertainty, third alsways true
        Amplitude:[peakamp,0,true],        //second element is uncertainty, third alsways true
        LandauAmplitude:[0,-1,false],      //always this value
        LandauMode:[0,-1,false],           //always this value
        LandauSigma:[0,-1,false],          //always this value
        Chi2:[peakchi2,-1,false],          //second and third elemtn always 1 and false         
        forCalibration:false,              //always this value
        forSourceFit:false,                //always this value
        sourceType:'NotSpecified'          //always this value 
        }
    );
  }

  self.redraw()();
}

SpectrumChartD3.prototype.handleTouchMovePeakFit = function() {
  var self = this;

  if (!self.rawData || !self.rawData.spectra)
    return;

  // Clear the delete peaks mode
  self.handleCancelTouchDeletePeak();

  // Clear the count gammas mode
  self.handleCancelTouchCountGammas();

  // Cancel the zoom-in y mode
  self.handleTouchCancelZoomInY();


  var t = d3.touches(self.vis[0][0]);

  // Cancel the function if no two-finger swipes detected
  if (!t || t.length !== 2 || !self.createPeaksStartTouches) {
    self.handleCancelTouchPeakFit();
    return;
  }

  // Set the touch variables
  var leftStartTouch = self.createPeaksStartTouches[0][0] < self.createPeaksStartTouches[1][0] ? self.createPeaksStartTouches[0] : self.createPeaksStartTouches[1],
      rightStartTouch = leftStartTouch === self.createPeaksStartTouches[0] ? self.createPeaksStartTouches[1] : self.createPeaksStartTouches[0];

  var leftTouch = t[0][0] < t[1][0] ? t[0] : t[1],
      rightTouch = leftTouch === t[0] ? t[1] : t[0];

  if (leftTouch[0] < leftStartTouch[0]) {
    self.handleCancelTouchPeakFit();
    return;
  }

  self.isCreatingPeaks = true;
  self.fittingPeak = true;

  // Uncomment if you want to have peak fitting animation within touch
  // self.handleMouseMovePeakFit();

  // Set the length of the arrows
  var arrowLength = 25,
      arrowPadding = 7;

  // Finger touch pixel size = 57
  var touchPointRadius = 20; 

  // To keep track of some of the line objects being drawn
  var createPeakTouchCurrentLine,
      contEstLine, 
      contEstText = d3.select("#contEstText"),
      createPeakTouchText1 = d3.select("#createPeakTouchText1"),
      createPeakTouchText2 = d3.select("#createPeakTouchText2"),
      createPeakRightTouchPointTop = d3.select("#createPeakRightTouchPointTop"),
      createPeakRightTouchPointBottom = d3.select("#createPeakRightTouchPointBottom"),
      createPeakTouchLineTop = d3.select("#createPeakTouchLineTop"), 
      createPeakTouchLineBottom = d3.select("#createPeakTouchLineBottom"),
      createPeakTouchArrowLineTop = d3.select("#createPeakTouchArrowLineTop"),
      createPeakTouchArrowLineBottom = d3.select("#createPeakTouchArrowLineBottom");

  if (createPeakTouchArrowLineTop.empty()) {
    createPeakTouchArrowLineTop = self.vis.append("line")
      .attr("id", "createPeakTouchArrowLineTop")
      .attr("class", "createPeakArrowLine")
      .attr("x1", leftStartTouch[0] - 12)
      .attr("x2", leftStartTouch[0] - 12)
      .attr("y1", leftStartTouch[1] - touchPointRadius - 20)
      .attr("y2", leftStartTouch[1] - touchPointRadius - 10);

    var createPeakTouchTopArrow = self.vis.append('svg:defs').attr("id", "createPeakStartArrowDef")
                          .append("svg:marker")
                          .attr("id", "createPeakTouchTopArrow")
                          .attr('class', 'createPeakArrow')
                          .attr("refX", 2)
                          .attr("refY", 8)
                          .attr("markerWidth", 25)
                          .attr("markerHeight", 25)
                          .attr("orient", 90)
                          .append("path")
                            .attr("d", "M2,2 L2,13 L8,7 L2,2")
                            .style("stroke", "black");

    createPeakTouchArrowLineTop.attr("marker-end", "url(#createPeakTouchTopArrow)");
  }
  if (createPeakTouchArrowLineBottom.empty()) {
    createPeakTouchArrowLineBottom = self.vis.append("line")
      .attr("id", "createPeakTouchArrowLineBottom")
      .attr("class", "createPeakArrowLine")
      .attr("x1", leftStartTouch[0] - 12)
      .attr("x2", leftStartTouch[0] - 12)
      .attr("y1", leftStartTouch[1] + touchPointRadius + 10)
      .attr("y2", leftStartTouch[1] + touchPointRadius + 20);

    var createPeakTouchBottomArrow = self.vis.append('svg:defs').attr("id", "createPeakStartArrowDef")
                          .append("svg:marker")
                          .attr("id", "createPeakTouchBottomArrow")
                          .attr('class', 'createPeakArrow')
                          .attr("refX", 2)
                          .attr("refY", 8)
                          .attr("markerWidth", 25)
                          .attr("markerHeight", 25)
                          .attr("orient", 270)
                          .append("path")
                            .attr("d", "M2,2 L2,13 L8,7 L2,2")
                            .style("stroke", "black");

    createPeakTouchArrowLineBottom.attr("marker-start", "url(#createPeakTouchBottomArrow)");
  }


  // Create the left topmost and bottommost touch point
  if (d3.select("#createPeakTouchPoint1").empty()) {
    self.vis.append("circle")
      .attr("id", "createPeakTouchPoint1")
      .attr("class", "createPeakLeftTouchPoint")
      .attr("cx", leftStartTouch[0])
      .attr("cy", leftStartTouch[1] - touchPointRadius)
      .attr("r", 3);
  }
  if (d3.select("#createPeakTouchPoint2").empty()) {
    self.vis.append("circle")
      .attr("id", "createPeakTouchPoint2")
      .attr("class", "createPeakLeftTouchPoint")
      .attr("cx", leftStartTouch[0])
      .attr("cy", leftStartTouch[1] + touchPointRadius)
      .attr("r", 3);
  }

  // Create/update the right topmost and bottommost touch point
  if (createPeakRightTouchPointTop.empty()) {
    createPeakRightTouchPointTop = self.vis.append("circle")
      .attr("id", "createPeakRightTouchPointTop")
      .attr("class", "createPeakRightTouchPoint")
      .attr("cy", leftStartTouch[1] - touchPointRadius)
      .attr("r", 3);
  }
  createPeakRightTouchPointTop.attr("cx", rightTouch[0]);

  if (createPeakRightTouchPointBottom.empty()) {
    createPeakRightTouchPointBottom = self.vis.append("circle")
      .attr("id", "createPeakRightTouchPointBottom")
      .attr("class", "createPeakRightTouchPoint")
      .attr("cy", leftStartTouch[1] + touchPointRadius)
      .attr("r", 3);
  }
  createPeakRightTouchPointBottom.attr("cx", rightTouch[0]);

  // Create/update the touch lines
  if (createPeakTouchLineTop.empty()) {
    createPeakTouchLineTop = self.vis.append("line")
      .attr("id", "createPeakTouchLineTop")
      .attr("class", "createPeakTouchLine")
      .attr("x1", leftStartTouch[0])
      .attr("y1", leftStartTouch[1] - touchPointRadius)
      .attr("y2", leftStartTouch[1] - touchPointRadius);
  }
  createPeakTouchLineTop.attr("x2", rightTouch[0]);

  if (createPeakTouchLineBottom.empty()) {
    createPeakTouchLineBottom = self.vis.append("line")
      .attr("id", "createPeakTouchLineBottom")
      .attr("class", "createPeakTouchLine")
      .attr("x1", leftStartTouch[0])
      .attr("y1", leftStartTouch[1] + touchPointRadius)
      .attr("y2", leftStartTouch[1] + touchPointRadius);
  }
  createPeakTouchLineBottom.attr("x2", rightTouch[0]);

  // Create the leftmost starting point line 
  if (d3.select("#createPeakTouchStartLine").empty()) {
    self.vis.append("line")
      .attr("id", "createPeakTouchStartLine")
      .attr("class", "createPeakMouseLine")
      .attr("x1", leftStartTouch[0])
      .attr("x2", leftStartTouch[0])
      .attr("y1", 0)
      .attr("y2", self.size.height);
  }

  // Create/refer the rightmost current point line
  if (d3.select("#createPeakTouchCurrentLine").empty()) {
    createPeakTouchCurrentLine = self.vis.append("line")
      .attr("id", "createPeakTouchCurrentLine")
      .attr("class", "createPeakMouseLine")
      .attr("y1", 0)
      .attr("y2", self.size.height);

  } else 
    createPeakTouchCurrentLine = d3.select("#createPeakTouchCurrentLine");

  createPeakTouchCurrentLine.attr("x1", rightTouch[0])
    .attr("x2", rightTouch[0]);

  // Create the text for the touch level
  if (createPeakTouchText1.empty())
    createPeakTouchText1 = self.vis.append("text")
      .attr("id", "createPeakTouchText1")
      .attr("class", "createPeakTouchText")
      .attr("y", leftStartTouch[1] - 2)
      .text("Keep fingers");
  createPeakTouchText1.attr("x", leftStartTouch[0] - createPeakTouchText1[0][0].clientWidth - 5);

  if (createPeakTouchText2.empty())
    createPeakTouchText2 = self.vis.append("text")
      .attr("id", "createPeakTouchText2")
      .attr("class", "createPeakTouchText")
      .attr("y", leftStartTouch[1] + 8)
      .text("level");
  createPeakTouchText2.attr("x", leftStartTouch[0] - createPeakTouchText2[0][0].clientWidth - 25);

  // Create/refer to the text for the approx continuum
  if (contEstText.empty()) {
    contEstText = self.vis.append("text")
      .attr("id", "contEstText")
      .attr("class", "contEstLineText")
      .text("approx continuum to use");
  }

  // Get pixelated coordinates of mouse/starting positions of lines and coordinates
  var coordsX0 = self.xScale.invert(leftStartTouch[0]),
      coordsY0 = self.getCountsForEnergy(self.rawData.spectra[0], coordsX0),
      coordsX1 = self.xScale.invert(rightTouch[0]),
      coordsY1 = self.getCountsForEnergy(self.rawData.spectra[0], coordsX1),
      x0 = leftStartTouch[0],
      x1 = rightTouch[0],
      y0 = self.yScale( coordsY0 ),
      y1 = self.yScale( coordsY1 ),
      dy = coordsY1 - coordsY0,
      lineAngle = (180/Math.PI) * Math.atan( (y1-y0)/(x1-x0) );

  // Update the position and rotation of the continuum text using the pixelated coordinates
  contEstText.attr("x", x0 + (Math.abs(x0-x1)/2) - Number(contEstText[0][0].clientWidth)/2 )
    .attr("y", y0-15)
    .attr("transform", "rotate(" + (!isNaN(lineAngle) ? lineAngle : 0) + ", " + (x0 + ((x1-x0)/2)) + ", " + y0 + ")");


  // Create/refer to the reference text for creating a peak
  if (d3.select("#createPeakTouchText").empty())
    createPeakTouchText = self.vis.append("text")
      .attr("id", "createPeakTouchText")
      .attr("class", "mouseLineText")
      .attr("y", self.size.height/5)
      .text("Will create peak Inside");
  else
    createPeakTouchText = d3.select("#createPeakTouchText");

  // Move the create peaks text in the middle of the create peak mouse lines
  createPeakTouchText.attr("x", x0 + (Math.abs(x0-x1)/2) - Number(createPeakTouchText[0][0].clientWidth)/2 );


  // Remove the continuum estimation line (lines in this case, we create the effect using a numver of shorter 5px lines)
  d3.selectAll(".createPeakContEstLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Get the end coordinates of the first (5px) line in the estimated continuum
  var xpix = x0+5,
      ypix = self.yScale( coordsY0+dy*(xpix-x0)/(x1-x0) );

  // Stop updating if the y-coordinate could not be found for the line
  if (isNaN(ypix))
    return;

  // Create the first 5px line for the estimated continuum
  if (x1 > x0)
    contEstLine = self.vis.append("line")
      .attr("id", "contEstLine")
      .attr("class", "createPeakContEstLine")
      .attr("x1", x0)
      .attr("y1", y0)
      .attr("x2", xpix)
      .attr("y2", ypix);

  // Create and update the next 5px lines until you reach the end for the estimated continuum
  while (xpix < x1) {
    x0 = xpix;
    y0 = ypix;
    xpix += 5;
    ypix = self.yScale( coordsY0+dy*(xpix-leftStartTouch[0])/(x1-leftStartTouch[0]) );

    contEstLine = self.vis.append("line")
      .attr("id", "contEstLine")
      .attr("class", "createPeakContEstLine")
      .attr("x1", x0)
      .attr("y1", y0)
      .attr("x2", xpix)
      .attr("y2", ypix);
  }

  // Create the last 5px line for the estimated continuum line
  contEstLine = self.vis.append("line")
    .attr("id", "contEstLine")
    .attr("class", "createPeakContEstLine")
    .attr("x1", x0)
    .attr("y1", y0)
    .attr("x2", x1-2)
    .attr("y2", y1);
}

SpectrumChartD3.prototype.handleTouchEndPeakFit = function() {
  var self = this;

  var t = self.lastTouches;

  // Cancel the function if no two-finger swipes detected
  if (!t || t.length !== 2 || !self.createPeaksStartTouches) {
    self.handleCancelTouchPeakFit();
    return;
  }

  // Set the touch variables
  var leftStartTouch = self.createPeaksStartTouches[0][0] < self.createPeaksStartTouches[1][0] ? self.createPeaksStartTouches[0] : self.createPeaksStartTouches[1],
      rightStartTouch = leftStartTouch === self.createPeaksStartTouches[0] ? self.createPeaksStartTouches[1] : self.createPeaksStartTouches[0];

  var leftTouch = t[0][0] < t[1][0] ? t[0] : t[1],
      rightTouch = leftTouch === t[0] ? t[1] : t[0];

  // Emit the create peak signal
  if (self.controlDragSwipe && self.createPeaksStartTouches && self.lastTouches) {
    console.log("Emit CREATE PEAK signal from x0 = ", leftStartTouch[0], "(", self.xScale.invert(leftStartTouch[0]), 
      " kEV) to x1 = ", rightTouch[0], "(", self.xScale.invert(rightTouch[0]), 
      " kEV)");
  }

  // Delete all the create peak elements
  self.handleCancelTouchPeakFit();
}

SpectrumChartD3.prototype.handleCancelTouchPeakFit = function() {
  var self = this;

  // Delete the leftmost start line
  d3.select("#createPeakTouchStartLine").remove();
  
  // Delete the right most current mouse line
  d3.select("#createPeakTouchCurrentLine").remove(); 

  // Delete the arrow definitions pointing to the mouse lines
  d3.select("#createPeakStartArrowDef").remove();

  // Delete the estimated continuum text
  d3.select("#contEstText").remove();

  // Delete the arrows pointing to the mouse lines
  d3.selectAll(".createPeakArrow").forEach(function (arrows) {
    arrows.forEach(function(arrow) {
      arrow.remove();
    })
  });

  // Delete the arrow lines for the previous arrows
  d3.selectAll(".createPeakArrowLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Delete the touch points
  d3.selectAll(".createPeakLeftTouchPoint").forEach(function (points) {
    points.forEach(function(point) {
      point.remove();
    })
  });
  d3.selectAll(".createPeakRightTouchPoint").forEach(function (points) {
    points.forEach(function(point) {
      point.remove();
    })
  });

  // Delete the touch lines
  d3.selectAll(".createPeakTouchLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Delete the all the lines for the estimated continuum
  d3.selectAll(".createPeakContEstLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Delete the reference text for the create peak
  d3.select("#createPeakTouchText").remove();
  d3.selectAll(".createPeakTouchText").forEach(function (texts) {
    texts.forEach(function(text) {
      text.remove();
    })
  });

  self.controlDragSwipe = false;
  self.fittingPeak = null;
}

//Function called when use lets the mouse button up
SpectrumChartD3.prototype.handleMouseUpPeakFit = function() {
  console.log( 'Mouse up during peak fit' );

  this.fittingPeak = null;
  
  var keep_peaks = true;  //Do what you want here to decide if you should keep the peak(s) or not
  
  if( keep_peaks && this.current_fitting_peak ) {
   this.rawData.spectra[0].peaks.push( this.current_fitting_peak ); 
  }

  if (this.peakFitMouseMove && this.peakFitMouseDown) {
    console.log("Emit CREATE PEAK signal from x0 = ", this.peakFitMouseDown[0], "(", this.xScale.invert(this.peakFitMouseDown[0]), 
      " kEV) to x1 = ", this.peakFitMouseMove[0], "(", this.xScale.invert(this.peakFitMouseMove[0]), 
      " kEV)");
  }  
  
  this.current_fitting_peak = null;
  this.erasePeakFitReferenceLines();
  this.redraw()();
}

//Function called when you hit escape while fitting peak
SpectrumChartD3.prototype.handleCancelMousePeakFit = function() {
  console.log( 'Canceled peakfit' );
  this.fittingPeak = null;
  this.current_fitting_peak = null;
  this.erasePeakFitReferenceLines();
  this.redraw()();
}

SpectrumChartD3.prototype.drawPeakFitReferenceLines = function() {
  var self = this;

  if (!self.rawData || !self.rawData.spectra)
    return;

  // Cancel the function if no mouse position detected
  if (!self.lastMouseMovePos || !self.leftMouseDown)
    return;

  // Adjust the mouse move position with respect to the bounds of the vis
  if (self.lastMouseMovePos[0] < 0)
    self.lastMouseMovePos[0] = 0;
  else if (self.lastMouseMovePos[0] > self.size.width)
    self.lastMouseMovePos[0] = self.size.width;

  // Set the length of the arrows
  var arrowLength = 25,
      arrowPadding = 7;

  // To keep track of some of the line objects being drawn
  var createPeakMouseCurrentLine,
      startLineArrowLine, currentLineArrowLine,
      startLineArrow, currentLineArrow,
      createPeakMouseText,
      contEstLine, contEstText;

  // Create the leftmost starting point line 
  if (d3.select("#createPeakMouseStartLine").empty()) {
    self.vis.append("line")
      .attr("id", "createPeakMouseStartLine")
      .attr("class", "createPeakMouseLine")
      .attr("x1", self.leftMouseDown[0])
      .attr("x2", self.leftMouseDown[0])
      .attr("y1", 0)
      .attr("y2", self.size.height);
  }

  // Create/refer the rightmost current point line
  if (d3.select("#createPeakMouseCurrentLine").empty()) {
    createPeakMouseCurrentLine = self.vis.append("line")
      .attr("id", "createPeakMouseCurrentLine")
      .attr("class", "createPeakMouseLine")
      .attr("x1", self.lastMouseMovePos[0])
      .attr("x2", self.lastMouseMovePos[0])
      .attr("y1", 0)
      .attr("y2", self.size.height);

  } else 
    createPeakMouseCurrentLine = d3.select("#createPeakMouseCurrentLine");

  createPeakMouseCurrentLine.attr("x1", self.lastMouseMovePos[0])
    .attr("x2", self.lastMouseMovePos[0]);

  // Create the start arrow point of the leftmost line (staticly positioned)
  if (d3.select("#startLineArrowLine").empty()) {
    startLineArrowLine = self.vis.append("line")
      .attr("id", "startLineArrowLine")
      .attr("class", "createPeakArrowLine")
      .attr("x1", self.leftMouseDown[0] - arrowLength - arrowPadding)
      .attr("x2", self.leftMouseDown[0] - arrowPadding)
      .attr("y1", self.size.height / 2)
      .attr("y2", self.size.height / 2);

    startLineArrow = self.vis.append('svg:defs').attr("id", "createPeakStartArrowDef")
      .append("svg:marker")
      .attr("id", "createPeakStartArrow")
      .attr('class', 'createPeakArrow')
      .attr("refX", 2)
      .attr("refY", 7)
      .attr("markerWidth", 25)
      .attr("markerHeight", 25)
      .attr("orient", 0)
      .append("path")
        .attr("d", "M2,2 L2,13 L8,7 L2,2")
        .style("stroke", "black");

    startLineArrowLine.attr("marker-end", "url(#createPeakStartArrow)");
  }

  // Create the arrow pointing to the rightmost line, or refer to it if itdoesn't exist
  if (d3.select("#currentLineArrowLine").empty()) {
    currentLineArrowLine = self.vis.append("line")
      .attr("id", "currentLineArrowLine")
      .attr("class", "createPeakArrowLine")
      .attr("y1", self.size.height / 2)
      .attr("y2", self.size.height / 2);

    currentLineArrow = self.vis.append('svg:defs').attr("id", "createPeakCurrentArrowDef")
      .append("svg:marker")
      .attr("id", "createPeakCurrentArrow")
      .attr('class', 'createPeakArrow')
      .attr("refX", 2)
      .attr("refY", 7)
      .attr("markerWidth", 25)
      .attr("markerHeight", 25)
      .attr("orient", 180)
      .append("path")
        .attr("d", "M2,2 L2,13 L8,7 L2,2")
        .style("stroke", "black");              

    currentLineArrowLine.attr("marker-start", "url(#createPeakCurrentArrow)");
  }
  else
    currentLineArrowLine = d3.select("#currentLineArrowLine");

  // Update the position of the rightmost arrow
  currentLineArrowLine.attr("x1", self.lastMouseMovePos[0] + arrowPadding)
    .attr("x2", self.lastMouseMovePos[0] + arrowPadding + arrowLength);


  // Create/refer to the text for the approx continuum
  if (d3.select("#contEstText").empty()) {
    contEstText = self.vis.append("text")
      .attr("id", "contEstText")
      .attr("class", "contEstLineText")
      .text("approx continuum to use");

  } else
    contEstText = d3.select("#contEstText");


  // Get pixelated coordinates of mouse/starting positions of lines and coordinates
  var coordsX0 = self.xScale.invert(self.leftMouseDown[0]),
      coordsY0 = self.getCountsForEnergy(self.rawData.spectra[0], coordsX0),
      coordsX1 = self.xScale.invert(self.lastMouseMovePos[0]),
      coordsY1 = self.getCountsForEnergy(self.rawData.spectra[0], coordsX1),
      x0 = self.leftMouseDown[0],
      x1 = self.lastMouseMovePos[0],
      y0 = self.yScale( coordsY0 ),
      y1 = self.yScale( coordsY1 ),
      dy = coordsY1 - coordsY0,
      lineAngle = (180/Math.PI) * Math.atan( (y1-y0)/(x1-x0) );

  // Update the position and rotation of the continuum text using the pixelated coordinates
  contEstText.attr("x", Math.min(self.lastMouseMovePos[0], self.leftMouseDown[0]) + (Math.abs(self.leftMouseDown[0]-self.lastMouseMovePos[0])/2) - Number(contEstText.node().getBBox().width)/2 )
    .attr("y", y0-15)
    .attr("transform", "rotate(" + (!isNaN(lineAngle) ? lineAngle : 0) + ", " + (x0 + ((x1-x0)/2)) + ", " + y0 + ")");


  // Create/refer to the reference text for creating a peak
  if (d3.select("#createPeakMouseText").empty())
    createPeakMouseText = self.vis.append("text")
      .attr("id", "createPeakMouseText")
      .attr("class", "mouseLineText")
      .attr("y", self.size.height/5)
      .text("Will create peak Inside");
  else
    createPeakMouseText = d3.select("#createPeakMouseText");

  // Move the create peaks text in the middle of the create peak mouse lines
  createPeakMouseText.attr("x", Math.min(self.lastMouseMovePos[0], self.leftMouseDown[0]) + (Math.abs(self.leftMouseDown[0]-self.lastMouseMovePos[0])/2) - Number(createPeakMouseText.node().getBBox().width)/2 );


  // Remove the continuum estimation line (lines in this case, we create the effect using a numver of shorter 5px lines)
  d3.selectAll(".createPeakContEstLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Get the end coordinates of the first (5px) line in the estimated continuum
  var xpix = x0+5,
      ypix = self.yScale( coordsY0+dy*(xpix-x0)/(x1-x0) );

  // Stop updating if the y-coordinate could not be found for the line
  if (isNaN(ypix))
    return;

  // Create the first 5px line for the estimated continuum
  if (x1 > x0)
    contEstLine = self.vis.append("line")
      .attr("id", "contEstLine")
      .attr("class", "createPeakContEstLine")
      .attr("x1", x0)
      .attr("y1", y0)
      .attr("x2", xpix)
      .attr("y2", ypix);

  // Create and update the next 5px lines until you reach the end for the estimated continuum
  while (xpix < x1) {
    x0 = xpix;
    y0 = ypix;
    xpix += 5;
    ypix = self.yScale( coordsY0+dy*(xpix-self.leftMouseDown[0])/(x1-self.leftMouseDown[0]) );

    contEstLine = self.vis.append("line")
      .attr("id", "contEstLine")
      .attr("class", "createPeakContEstLine")
      .attr("x1", x0)
      .attr("y1", y0)
      .attr("x2", xpix)
      .attr("y2", ypix);
  }

  // Create the last 5px line for the estimated continuum line
  contEstLine = self.vis.append("line")
    .attr("id", "contEstLine")
    .attr("class", "createPeakContEstLine")
    .attr("x1", x0)
    .attr("y1", y0)
    .attr("x2", x1-2)
    .attr("y2", y1);
}

SpectrumChartD3.prototype.erasePeakFitReferenceLines = function() {
  var self = this;

  // Delete the leftmost start line
  d3.select("#createPeakMouseStartLine").remove();
  
  // Delete the right most current mouse line
  d3.select("#createPeakMouseCurrentLine").remove(); 

  // Delete the arrow definitions pointing to the mouse lines
  d3.select("#createPeakStartArrowDef").remove();

  // Delete the estimated continuum text
  d3.select("#contEstText").remove();

  // Delete the arrows pointing to the mouse lines
  d3.selectAll(".createPeakArrow").forEach(function (arrows) {
    arrows.forEach(function(arrow) {
      arrow.remove();
    })
  });

  // Delete the arrow lines for the previous arrows
  d3.selectAll(".createPeakArrowLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Delete the all the lines for the estimated continuum
  d3.selectAll(".createPeakContEstLine").forEach(function (lines) {
    lines.forEach(function(line) {
      line.remove();
    })
  });

  // Delete the reference text for the create peak
  d3.select("#createPeakMouseText").remove();
}


/* Zoom In X-Axis functions */
SpectrumChartD3.prototype.handleMouseMoveZoomInX = function () {
  var self = this;

  var zoomInXBox = d3.select("#zoomInXBox"),
      zoomInXText = d3.select("#zoomInXText");

  // Cancel erase peaks mode, we're zooming in
  self.handleCancelMouseDeletePeak();
  // Cancel recalibration mode, we're zooming in
  self.handleCancelMouseRecalibration();

  // Cancel the zooming in y mode
  self.handleCancelMouseZoomInY();

  // Cancel the count gammas mode
  self.handleCancelMouseCountGammas();

  if (!self.zoominmouse)  return;

  // Adjust the mouse move position with respect to the bounds of the vis
  if (self.lastMouseMovePos[0] < 0)
    self.lastMouseMovePos[0] = 0;
  else if (self.lastMouseMovePos[0] > self.size.width)
    self.lastMouseMovePos[0] = self.size.width;

  // We are now zooming in
  self.zooming_plot = true;

  // Restore the zoom-in box
  if (zoomInXBox.empty()) {
    zoomInXBox = self.vis.append("rect")
      .attr("id", "zoomInXBox")
      .attr("class","leftbuttonzoombox")
      .attr("width", 1 )
      .attr("height", self.size.height)
      .attr("x", self.lastMouseMovePos[0])
      .attr("y", 0)
      .attr("pointer-events", "none");
  }
  
  // If click-and-drag reaches out of bounds from the plot
  if (self.lastMouseMovePos[0] < 0 || self.lastMouseMovePos[0] > self.size.width) {
    // Do something if click and drag reaches out of bounds from plot
  }


  if( self.lastMouseMovePos[0] < self.zoominmouse[0] ) {      // If the mouse position is less than the zoombox starting position (zoom-out)

    // If we were animating a zoom, cancel that animation
    self.handleCancelAnimationZoom();

    // Remove the zoom-in x-axis text, we are zooming out
    zoomInXText.remove();

    self.zoominaltereddomain = true;
    zoomInXBox.attr("x", self.lastMouseMovePos[0])
      .attr("class","leftbuttonzoomoutbox")
      .attr("height", self.size.height)
      .attr("width", self.zoominmouse[0] - self.lastMouseMovePos[0] );

    //Do some zooming out
    var bounds = self.min_max_x_values();
    var xaxismin = bounds[0],
        xaxismax = bounds[1];

    var frac = 4 * (self.zoominmouse[0] - self.lastMouseMovePos[0]) / self.zoominmouse[0];

    if( !isFinite(frac) || isNaN(frac) || frac < 0.0 )
      frac = 0;

    var origdx = self.origdomain[1] - self.origdomain[0],
        maxdx = xaxismax - xaxismin,
        newdx = origdx + frac*(maxdx - origdx);

    var deltadx = newdx - origdx,
        newxmin = Math.max( self.origdomain[0] - 0.5*deltadx, xaxismin ),
        newxmax = Math.min( self.origdomain[1] + 0.5*deltadx, xaxismax );
    var olddomain = self.xScale.domain();
    var newdomain = [newxmin,newxmax];

    if (olddomain[0] != newdomain[0] || olddomain[1] != newdomain[1]) {
      self.xScale.domain( newdomain );
      self.redraw()();
    }

  } else {  
    // Set the zoombox (this gets called when current mouse position is greater than where the zoombox started) 
    if( self.zoominaltereddomain ) {
      var shouldRedraw = self.xScale.domain()[0] != self.origdomain[0] || self.xScale.domain()[1] != self.origdomain[1];
      self.xScale.domain( self.origdomain );
      self.zoominaltereddomain = false;
      zoomInXBox.attr("class","leftbuttonzoombox");
      if (shouldRedraw)
        self.redraw()();
    }

    // Update the zoomin-box x-position and width
    zoomInXBox.attr("x", self.zoominmouse[0])
      .attr("width", self.lastMouseMovePos[0] - self.zoominmouse[0] );

    if (self.lastMouseMovePos[0] - self.zoominmouse[0] > 7) {   // if zoom in box is at least 7px wide, update it
      if (zoomInXText.empty()) {
        zoomInXText = self.vis.append("text")
          .attr("id", "zoomInXText")
          .attr("class", "chartLineText")
          .attr("y", Number(zoomInXBox.attr("height"))/2)
          .text("Zoom In");
      }

      // keep zoom in label centered on the box
       zoomInXText.attr("x", Number(zoomInXBox.attr("x")) + (Number(zoomInXBox.attr("width"))/2) - (zoomInXText.node().getBBox().width/2) );
      //zoomInXText.attr("x", ((self.zoominmouse[0] + self.lastMouseMovePos[0])/2) - 30 );

    } else if (!zoomInXText.empty())       // delete if zoom in box not wide enough (eg. will not zoom)
      zoomInXText.remove();

  }
}

SpectrumChartD3.prototype.handleMouseUpZoomInX = function () {
  var self = this;

  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  if (self.zooming_plot) {
    var foreground = self.rawData.spectra[0];

    var lowerchanval = d3.bisector(function(d){return d.x0;}).left(foreground.points,self.xScale.invert(self.zoominmouse[0]),1) - 1;
    var higherchanval = d3.bisector(function(d){return d.x0;}).left(foreground.points,self.xScale.invert(self.lastMouseMovePos[0]),1) - 1;
    
    var yMinAtZoomRange = self.yScale(d3.min(foreground.points.slice(lowerchanval, higherchanval) , function(p) { return p.y; }));
    var yMaxAtZoomRange = self.yScale(d3.max(foreground.points.slice(lowerchanval, higherchanval) , function(p) { return p.y; }));

    // Get the current mouse position
    var m = d3.mouse(self.vis[0][0]);

    if( m[0] < self.zoominmouse[0] ) {
      //we were zooming out, nothing to do here

    } else {
      m[0] = m[0] < 0 ? 0 : m[0];

      var oldXScale = self.xScale.domain();

      //This is a big hack, because self.xScale.invert(m[0]) == self.zoominx0 otherwise.  I dont completely understand why!
      var x0 = self.xScale.invert(self.zoominmouse[0]),
          x1 = self.xScale.invert(m[0]);

      //require the mouse to have moved at least 6 pixels in x
      if( m[0] - self.zoominmouse[0] > 7 && self.rawData ) {
          
        var bounds = self.min_max_x_values();
        var mindatax = bounds[0], 
            maxdatax = bounds[1];
        
        //Make sure the new scale will span at least one bin , if not 
        //  make it just less than one bin, centered on that bin.
        var rawbi = d3.bisector(function(d){return d;}); 
        var lbin = rawbi.left(foreground.x, x0+(x1-x0));
        if( lbin > 1 && lbin < (foreground.x.length-1) && lbin === rawbi.left(foreground.x,x1+(x1-x0)) ) {
          //This doesnt actually work correctly, it doesnt actually center on the bin, meaning you're always a little off from where you want ...
          var corx0 = foreground.x[lbin-1];
          var corx1 = foreground.x[lbin];
          var p = 0.01*(corx1-corx0);
          corx0 += p;
          corx1 -= p;
          
          console.log( 'changing x0=' + x0 + ' and x1=' + x1 + " will set domain to " + (x0+(x1-x0)) + " through " + (x1+(x1-x0)) + " m[0]=" + m[0] + " self.zoominmouse[0]=" + self.zoominmouse[0] );
          self.xScale.domain( [Math.max(corx0+(x1-x0), mindatax),Math.min(corx1+(x1-x0), maxdatax)] );
        } else {
          self.xScale.domain([x0,x1]);
        }  

        // Update the peak markers
        self.updateFeatureMarkers(-1);
        
        if (self.options.showAnimation) {     // Draw zoom animations if option is checked

          // Cancel any current zoom animations
          self.handleCancelAnimationZoom();

          // Start new zoom animation
          self.currentDomain = self.savedDomain = oldXScale;
          self.zoomAnimationID = requestAnimationFrame(self.redrawZoomXAnimation([x0,x1]));
          self.startAnimationZoomTime = Math.floor(Date.now());

        } else { self.redraw()(); }   // Zoom animation unchecked; draw new x-axis range
      }
    }
  }


  self.handleCancelMouseZoomInX();

  self.zooming_plot = false;
}

SpectrumChartD3.prototype.handleCancelMouseZoomInX = function() {
  var self = this;

  var zoomInXBox = d3.select("#zoomInXBox"),
      zoomInXText = d3.select("#zoomInXText");

  // Delete zoom in box and text
  zoomInXText.remove();
  zoomInXBox.remove();

  // Not zooming in plot anymore
  self.zooming_plot = false;
}



/* Zoom In Y-Axis functions */
SpectrumChartD3.prototype.handleMouseMoveZoomInY = function () {
  var self = this;

  // Set the objects displayed for zooming in the y-axis
  var zoomInYBox = d3.select("#zoomInYBox"),
      zoomInYText = d3.select("#zoomInYText");

  // Get the mouse coordinates
  var m = d3.mouse(self.vis[0][0]);

  // Cancel the zooming mode
  self.handleCancelMouseZoomInX();

  // Cancel erase peaks mode, we're zooming in
  self.handleCancelMouseDeletePeak();

  // Cancel recalibration mode, we're zooming in
  self.handleCancelMouseRecalibration();

  // Cancel the count gammas mode
  self.handleCancelMouseCountGammas();

  // Now zooming in y mode
  self.zoomingYPlot = true;

  // Create the reference for the start of the zoom-in (Yaxis) box
  if (!self.zoomInYMouse)
    self.zoomInYMouse = self.lastMouseMovePos;

  // Adjust the zoom in y mouse in case user starts dragging from out of bounds
  if (self.zoomInYMouse[1] < 0)
    self.zoomInYMouse[1] = 0;
  else if (self.zoomInYMouse[1] > self.size.height)
    self.zoomInYMouse[1] = self.size.height;

  if (self.zoomInYMouse[0] < 0 || self.zoomInYMouse[0] > self.size.width) {
    self.handleCancelMouseZoomInY();
    return;
  }

  var height = self.lastMouseMovePos[1] - self.zoomInYMouse[1];

  // Set the zoom-y box
  if (zoomInYBox.empty()) {
    zoomInYBox = self.vis.append("rect")
      .attr("id", "zoomInYBox")
      .attr("class","leftbuttonzoombox")
      .attr("width", self.size.width )
      .attr("height", 0)
      .attr("x", 0)
      .attr("y", self.zoomInYMouse[1])
      .attr("pointer-events", "none");
  } else {
    if (height >= 0)  zoomInYBox.attr("height", height);
    else              zoomInYBox.attr("y", self.lastMouseMovePos[1])
                        .attr("height", Math.abs(height));
  }

  if (Math.abs(self.lastMouseMovePos[1] - self.zoomInYMouse[1]) > 10) {   // if zoom in box is at least 7px wide, update it
    if (zoomInYText.empty()) {
      zoomInYText = self.vis.append("text")
        .attr("id", "zoomInYText")
        .attr("class", "chartLineText");
      zoomInYText.attr("x", self.size.width/2 - (zoomInYText.node().getBBox().width/2));

    }

    // keep zoom in label centered on the box
    zoomInYText.attr("y", Number(zoomInYBox.attr("y")) + (Number(zoomInYBox.attr("height"))/2) );

    if (height > 0) {
      zoomInYText.text("Zoom-In on Y-axis");
      zoomInYBox.attr("class", "leftbuttonzoombox");
    }
    else {
      var zoomOutText = "Zoom-out on Y-axis";
      if (-height < 0.05*self.size.height)
        zoomOutText += " x2";
      else if (-height < 0.075*self.size.height)
        zoomOutText += " x4";
      else
        zoomOutText += " full";

      zoomInYText.text(zoomOutText);
      zoomInYBox.attr("class", "leftbuttonzoomoutboxy");
    }

  } else if (!zoomInYText.empty()) {       // delete if zoom in box not wide enough (eg. will not zoom)
    zoomInYText.remove();
  }
}

SpectrumChartD3.prototype.handleMouseUpZoomInY = function () {
  var self = this;

  var zoomInYBox = d3.select("#zoomInYBox"),
      zoomInYText = d3.select("#zoomInYText");

  /*
  Christian: This function is a similar method to the redraw function, except without certain actions
  taken out to only acount for redrawing the the zoomed chart by y-axis.
  */
  function redrawYAxis() {
    return function() {
      self.do_rebin();
      var tx = function(d) { return "translate(" + self.xScale(d) + ",0)"; };
      var stroke = function(d) { return d ? "#ccc" : "#666"; };

      self.yAxisBody.call(self.yAxis);
      self.yAxisBody.selectAll("text").style("cursor", "ns-resize")
          .on("mouseover", function(d) { d3.select(this).style("font-weight", "bold");})
          .on("mouseout",  function(d) { d3.select(this).style("font-weight", "normal");})
          .on("mousedown.drag",  self.yaxisDrag())
          .on("touchstart.drag", self.yaxisDrag());

      self.drawYTicks();
      
      self.calcLeftPadding( true );

      self.drawXTicks();

      self.drawXAxisArrows();
      
      self.drawPeaks();
      self.drawRefGammaLines();
      self.updateMouseCoordText();

      self.update(true);

      self.yAxisZoomedOutFully = false;
    }
  }

  if (zoomInYBox.empty()) {
    self.handleCancelMouseZoomInY();
    return;
  }

  // Set the y-values for where zoom-in occured
  var ypix1 = Math.min(Number(zoomInYBox.attr("y")), self.lastMouseMovePos[1]),
      ypix2 = Math.max(Number(zoomInYBox.attr("y")), self.lastMouseMovePos[1]),
      y1 = Math.min(self.yScale.invert(Number(zoomInYBox.attr("y"))), 
                    self.yScale.invert(self.lastMouseMovePos[1]) ),
      y2 = Math.max(self.yScale.invert(Number(zoomInYBox.attr("y"))), 
                    self.yScale.invert(self.lastMouseMovePos[1]) );

  var oldY1 = self.yScale.domain()[0];
  var oldY2 = self.yScale.domain()[1];

  if (self.zoomingYPlot) {
    if (y2 > y1 && Math.abs(ypix2-ypix1) > 10) {  // we are zooming in
      console.log("Emit Y AXIS ZOOM-IN signal!");

      if (self.options.showAnimation) {
        self.currentDomain = self.savedDomain = [oldY1, oldY2];

        if (self.zoomAnimationID != null) {
          cancelAnimationFrame(self.zoomAnimationID);
        }
        self.zoomAnimationID = requestAnimationFrame(self.redrawZoomInYAnimation([y2,y1], redrawYAxis()));
        self.startAnimationZoomTime = Math.floor(Date.now());

      } else { self.yScale.domain([y2,y1]); redrawYAxis()(); }

    } else {  // we are zooming out

      if (zoomInYText.empty()) {
        self.handleCancelMouseZoomInY();
        return;
      }

      // Zoom out completely if user dragged up a considerable amount
      if (zoomInYText.text().endsWith("full")) {

        // If zoom animations are checked, animate the zoom out
        if (self.options.showAnimation) {

          // Cancel previous zoom animation if there was one
          self.handleCancelAnimationZoom();

          // Start new zoom animation
          self.currentDomain = self.savedDomain = [oldY1, oldY2];
          self.zoomAnimationID = requestAnimationFrame(self.redrawZoomOutYAnimation([self.getYAxisDomain()[0],self.getYAxisDomain()[1]], redrawYAxis()));
          self.startAnimationZoomTime = Math.floor(Date.now());

        } else { self.redraw()(); }   // Redraw the chart normally if no zoom animations created

      } else { // Zoom out a portion of the full y-axis

        // This represents how much of the y-axis we will be zooming out
        var mult;

        // Get the mult value here
        // For some reason, setting the mult = 1 leads to buggy behavior, so here I set it to x2, x4 respectively
        if (zoomInYText.text().endsWith("x2"))  mult = 2; 
        else  mult = 4;

        // Get the old values of the current y-domain
        var oldY0 = self.yScale.domain()[0],
            oldY1 = self.yScale.domain()[1],
            oldRange = Math.abs(oldY1 - oldY0),
            centroid = oldY0 + 0.5*oldRange;

        // Values for the minimum and maximum y values
        var minY, maxY;

        if( !self.rawData || !self.rawData.y ){   // Manually set the min/max y values
          minY = 0.1;
          maxY = 3000;

        } else {                                  // Get the min, max y values from the data
          minY = self.getYAxisDomain()[1];
          maxY = self.getYAxisDomain()[0];
        }

        // Set the values for the new y-domain
        var newY0 = Math.max(centroid - mult*oldRange,
                             minY),
            newY1 = Math.min(newY0 + 2*mult*oldRange,
                             maxY);

        // Accounting for bounds checking, making sure new y domain does not reach beyond min/max values
        if (newY1 < 0 || newY1 > maxY)
          newY1 = maxY;
        if (newY0 > maxY || newY0 < 0)
          newY0 = minY;

        // Redraw the newly y-zoomed chart
        if (!(newY0 == minY && newY1 == maxY)) {

          // If zoom animations are checked, animate the zoom out
          if (self.options.showAnimation) {

            // Cancel previous zoom animation if there was one
            self.handleCancelAnimationZoom();

            // Start new zoom animation
            self.currentDomain = self.savedDomain = [oldY0, oldY1];
            self.zoomAnimationID = requestAnimationFrame(self.redrawZoomOutYAnimation([newY1, newY0], redrawYAxis()));
            self.startAnimationZoomTime = Math.floor(Date.now());

          } else { redrawYAxis()(); }   // Redraw the chart normally if no zoom animations created

        } else {  
          // If the new y domain reaches the min/max values and we are still trying to zoom out,
          //    then we attempt to zoom out again to capture the initial starting point of the chart
          self.redraw()();
        }
      }

    }
  }

  // Clean up objets from zooming in y-axis
  self.handleCancelMouseZoomInY();
}

SpectrumChartD3.prototype.handleCancelMouseZoomInY = function() {
  var self = this;

  // Set the objects displayed for zooming in the y-axis
  var zoomInYBox = d3.select("#zoomInYBox"),
      zoomInYText = d3.select("#zoomInYText");

  // Not zooming in y anymore
  self.zoomingYPlot = false;
  self.zoomInYMouse = null;

  // Delete zoom box and text
  zoomInYBox.remove();
  zoomInYText.remove();
}

SpectrumChartD3.prototype.handleTouchMoveZoomInY = function() {
  var self = this;

  if (!self.touchesOnChart)
    return;

  // Cancel delete peaks mode
  self.handleCancelTouchDeletePeak();

  // Cancel the create peaks mode
  self.handleCancelTouchPeakFit();

  // Cancel the count gammas mode
  self.handleCancelTouchCountGammas();

  var t = d3.touches(self.vis[0][0]);

  var keys = Object.keys(self.touchesOnChart);

  if (keys.length !== 2)
    return;
  if (t.length !== 2)
    return;

  var zoomInYTopLine = d3.select("#zoomInYTopLine");
  var zoomInYBottomLine = d3.select("#zoomInYBottomLine");
  var zoomInYText = d3.select("#zoomInYText");

  var touch1 = self.touchesOnChart[keys[0]];
  var touch2 = self.touchesOnChart[keys[1]];
  var adx1 = Math.abs( touch1.startX - touch2.startX );
  var adx2 = Math.abs( touch1.pageX  - touch2.pageX );
  var ady1 = Math.abs( touch1.startY - touch2.startY );
  var ady2 = Math.abs( touch1.pageY  - touch2.pageY );
  var ddx = Math.abs( adx2 - adx1 );
  var ddy = Math.abs( ady2 - ady1 );
  var areVertical = (adx2 > ady2);

  if (!touch1.visY)
    touch1.visY = t[0][1];
  if (!touch2.visY)
    touch2.visY = t[1][1];

  var topTouch = touch1.pageY < touch2.pageY ? touch1 : touch2;
  var bottomTouch = topTouch == touch1 ? touch2 : touch1;

  if (zoomInYTopLine.empty()) {
    zoomInYTopLine = self.vis.append("line")
      .attr("id", "zoomInYTopLine")
      .attr("class", "mouseLine")
      .attr("x1", 0)
      .attr("x2", self.size.width)
      .attr("y1", topTouch.visY)
      .attr("y2", topTouch.visY)
      .attr("count", self.yScale.invert(topTouch.visY));
  }

  if (zoomInYBottomLine.empty()) {
    zoomInYBottomLine = self.vis.append("line")
      .attr("id", "zoomInYBottomLine")
      .attr("class", "mouseLine")
      .attr("x1", 0)
      .attr("x2", self.size.width)
      .attr("y1", bottomTouch.visY)
      .attr("y2", bottomTouch.visY)
      .attr("count", self.yScale.invert(bottomTouch.visY));
  }

  if (zoomInYText.empty()) {
    zoomInYText = self.vis.append("text")
      .attr("id", "zoomInYText")
      .attr("class", "mouseLineText");
  }

  zoomInYText.text(function() {
    if (topTouch.visY > topTouch.startY && bottomTouch.visY < bottomTouch.startY)
      return "Zoom-Out on Y-axis";
    else if (topTouch.visY == topTouch.startY && bottomTouch.visY == bottomTouch.startY)
      return "";
    else
      return "Zoom-In on Y-axis";
  });

  zoomInYText.attr("x", self.size.width/2 - Number(zoomInYText.node().getBBox().width)/2)
    .attr("y", Number(topTouch.startY) + (bottomTouch.startY - topTouch.startY)/2);
}

SpectrumChartD3.prototype.handleTouchEndZoomInY = function() {
  var self = this;

  var zoomInYTopLine = d3.select("#zoomInYTopLine");
  var zoomInYBottomLine = d3.select("#zoomInYBottomLine");
  var zoomInYText = d3.select("#zoomInYText");

  if (zoomInYTopLine.empty() || zoomInYBottomLine.empty()) {
    self.handleTouchCancelZoomInY();
    return;
  }
  // Set the y-values for where zoom-in occured
  var ypix1 = Number(zoomInYTopLine.attr("y1")),
      ypix2 = Number(zoomInYBottomLine.attr("y1")),
      y1 = Number(zoomInYTopLine.attr("count")),
      y2 = Number(zoomInYBottomLine.attr("count"));

  console.log(y1, y2);
  /*
  This function is a similar method to the redraw function, except without certain actions
  taken out to only acount for redrawing the the zoomed chart by y-axis.
  */
  function redrawYAxis() {
    self.do_rebin();
    var tx = function(d) { return "translate(" + self.xScale(d) + ",0)"; };
    var stroke = function(d) { return d ? "#ccc" : "#666"; };

    self.yAxisBody.call(self.yAxis);
    self.yAxisBody.selectAll("text").style("cursor", "ns-resize")
        .on("mouseover", function(d) { d3.select(this).style("font-weight", "bold");})
        .on("mouseout",  function(d) { d3.select(this).style("font-weight", "normal");})
        .on("mousedown.drag",  self.yaxisDrag())
        .on("touchstart.drag", self.yaxisDrag());

    self.drawYTicks();
    
    self.calcLeftPadding( true );

    self.drawXTicks();

    self.drawXAxisArrows();
    
    self.drawPeaks();
    self.drawRefGammaLines();
    self.updateMouseCoordText();

    self.update();

    self.yAxisZoomedOutFully = false;
  }

  if (Math.abs(ypix2-ypix1) > 10) {  // we are zooming in
    if (zoomInYText.text().includes("Zoom-In")) {
      self.yScale.domain([y1,y2]);
      redrawYAxis();

    } else if (zoomInYText.text().includes("Zoom-Out")) {
        console.log("zooming out!");
        self.redraw()();
    }
  }
  self.handleTouchCancelZoomInY();
}

SpectrumChartD3.prototype.handleTouchCancelZoomInY = function() {
  var self = this;

  var zoomInYTopLine = d3.select("#zoomInYTopLine");
  var zoomInYBottomLine = d3.select("#zoomInYBottomLine");
  var zoomInYText = d3.select("#zoomInYText");

  zoomInYTopLine.remove();
  zoomInYBottomLine.remove();
  zoomInYText.remove();

  self.zoomInYPinch = null;
}


/* Energy Recalibration functions */
SpectrumChartD3.prototype.handleMouseMoveRecalibration = function() {
  var self = this;

  // Clear the zoom, we're recalibrating the chart
  self.handleCancelMouseZoomInX();

  // Cancel erase peaks mode, we're recalibrating the chart
  self.handleCancelMouseDeletePeak();

  // Cancel the zooming in y mode
  self.handleCancelMouseZoomInY();

  // Cancel the count gammas mode
  self.handleCancelMouseCountGammas();

  if (!self.recalibrationMousePos)
    return;
  if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length) 
    return;

  // Adjust the mouse move position with respect to the bounds of the vis
  if (self.lastMouseMovePos[0] < 0)
    self.lastMouseMovePos[0] = 0;
  else if (self.lastMouseMovePos[0] > self.size.width)
    self.lastMouseMovePos[0] = self.size.width;

  // Set the line objects to be referenced
  var recalibrationStartLine = d3.select("#recalibrationStartLine"),
      recalibrationText = d3.select("#recalibrationText"),
      recalibrationMousePosLines = d3.select("#recalibrationMousePosLines");

  var recalibrationG = d3.select("#recalibrationG");
  var recalibrationPeakVis = d3.select("#recalibrationPeakVis");

  // Set the line that symbolizes where user initially began right-click-and-drag
  if (recalibrationStartLine.empty()) {
    recalibrationStartLine = self.vis.append("line")
      .attr("id", "recalibrationStartLine")
      .attr("class", "mouseLine")
      .attr("x1", self.recalibrationMousePos[0])
      .attr("x2", self.recalibrationMousePos[0])
      .attr("y1", 0)
      .attr("y2", self.size.height);
  }

  if (recalibrationText.empty()) {                       // Right-click-and-drag text to say where recalibration ranges are
    recalibrationText = self.vis.append("text")
      .attr("id", "recalibrationText")
      .attr("class", "mouseLineText")
      .attr("x", self.recalibrationMousePos[0] + 5 /* As padding from the starting line */ )
      .attr("y", self.size.height/2)
      .text("Recalibrate data from " + self.xScale.invert(self.lastMouseMovePos[0]).toFixed(2) + " to " + self.xScale.invert(self.lastMouseMovePos[0]).toFixed(2) + " keV");

  } else   // Right-click-and-drag text already visible in vis, so just update it
    recalibrationText.text( "Recalibrate data from " + self.xScale.invert(recalibrationStartLine.attr("x1")).toFixed(2) + " to " + self.xScale.invert(self.lastMouseMovePos[0]).toFixed(2) + " keV" )
  

  // Draw the distance lines from the right-click-and-drag line to mouse position
  if (!self.recalibrationDistanceLines) {
    var distanceLine;
    self.recalibrationDistanceLines = [];

    var lineSpacing = self.size.height / 5;

    // Create 4 lines that draw from the right-click-and-drag initial starting point to where the mouse is
    for (i = 0; i < 4; i++) {
      distanceLine = self.vis.append("line")
                            .attr("class", "secondaryMouseLine")
                            .attr("x1", recalibrationStartLine.attr("x1"))
                            .attr("x2", self.lastMouseMovePos[0])
                            .attr("y1", lineSpacing)
                            .attr("y2", lineSpacing);

      // This is to add the the arrowhead end to the distance lines
      var arrow = self.vis.append('svg:defs')
        .attr("id", "recalibrationArrowDef")
        .append("svg:marker")
          .attr("id", "rightClickDragArrow")
          .attr('class', 'xaxisarrow')
          .attr("refX", 0)
          .attr("refY", 7)
          .attr("markerWidth", 25)
          .attr("markerHeight", 25)
          .attr("orient", 0)
          .append("path")
            .attr("d", "M2,2 L2,13 L8,7 L2,2")
            .style("stroke", "black");

      distanceLine.attr("marker-end", "url(#rightClickDragArrow)");

      self.recalibrationDistanceLines.push(distanceLine);
      lineSpacing += self.size.height / 5;
    }

  } else {    // Distance lines have already been drawn
    for (i = 0; i < self.recalibrationDistanceLines.length; i++) {
      var x2 = self.lastMouseMovePos[0];

      // If mouse position is to the left of the right-click-and-drag line
      if (self.lastMouseMovePos[0] < recalibrationStartLine.attr("x1")) {
        self.recalibrationDistanceLines[i].attr("x2", x2);     // Adjust the x-position of the line

        d3.selectAll("#rightClickDragArrow > path").each(function(){            // Flip the arrowhead to face towards the negative x-axis
          d3.select(this).attr("transform", "translate(8,14) rotate(180)");
        });
      }
      else {
        // To adjust for the length of the line with respect to its arrowhead
        if (self.recalibrationDistanceLines[i].attr("x2") > 0)
          x2 -= 8;  // Minus 8 to account for the arrow width connected to the line
        else
          x2 = 0;

        // Adjust the x-position of the line
        self.recalibrationDistanceLines[i].attr("x2", x2);

        // Un-flip the arrowhead (if it was flipped) back to pointing towards the positive x-axis
        d3.selectAll("#rightClickDragArrow > path").each(function(){
          d3.select(this).attr("transform", null);
        });
      }
    }
  }

  // Draw the line to represent the mouse position for recalibration
  if (recalibrationMousePosLines.empty())
    recalibrationMousePosLines = self.vis.append("line")
      .attr("id", "recalibrationMousePosLines")
      .attr("class", "lightMouseLine")
      .attr("y1", 0)
      .attr("y2", self.size.height);
   
  // Update the mouse position line for recalibration
  recalibrationMousePosLines.attr("x1", self.lastMouseMovePos[0])
    .attr("x2", self.lastMouseMovePos[0]);

  // Add the background for the recalibration animation
  if (recalibrationG.empty()) {
    recalibrationG = self.vis.append("g")
      .attr("id", "recalibrationG")
      .attr("clip-path", "url(#clip" + this.chart.id + ")");
  }

  // Add the foreground, background, and secondary lines for recalibration animation
  for (var i = 0; i < self.rawData.spectra.length; ++i) {
    var spectrum = self.rawData.spectra[i];
    var recalibrationLine = d3.select("#recalibrationLine"+i);

    if (recalibrationLine.empty() && self['line'+i]) {
      recalibrationLine = recalibrationG.append("path")
        .attr("id", "recalibrationLine"+i)
        .attr("class", "rline")
        .attr("stroke", spectrum.lineColor ? spectrum.lineColor : 'black')
        .attr("d", self['line'+i](spectrum.points));
    }
  }

  if (recalibrationPeakVis.empty() && self.peakVis) {
    recalibrationPeakVis = recalibrationG.append("g")
      .attr("id", "recalibrationPeakVis")
      .attr("class", "peakVis")
      .attr("transform","translate(0,0)")
      .attr("clip-path", "url(#clip" + this.chart.id + ")");

    self.peakVis.selectAll("path").each(function() {
      path = d3.select(this);
      recalibrationPeakVis.append("path")
        .attr("class", path.attr("class"))
        .attr("d", path.attr("d"))
        .attr("fill-opacity", 0.4)
        .style("fill", path.style("fill"));
    });
  }

  recalibrationPeakVis.attr("transform", "translate(" + (self.lastMouseMovePos[0] - self.recalibrationMousePos[0]) + ",0)");

  // Move the foreground, background, and secondary lines for recalibration animation with relation to mouse position
  for (var i = 0; i < self.rawData.spectra.length; ++i) {
    var recalibrationLine = d3.select("#recalibrationLine"+i);

    if (!recalibrationLine.empty())
      recalibrationLine.attr("transform", "translate(" + (self.lastMouseMovePos[0] - self.recalibrationMousePos[0]) + ",0)");
  }
}

SpectrumChartD3.prototype.handleMouseUpRecalibration = function() {
  var self = this;

  // Handle Right-click-and-drag (for recalibrating data)
  if (self.recalibrationMousePos) {

    // Emit the signal here
    if (self.isRecalibrating) {
      console.log("Emit RECALIBRATION SIGNAL from x0 = ", self.xScale(self.recalibrationStartEnergy[0]), 
        "(", self.recalibrationStartEnergy[0], " keV) to x1 = ", 
        self.lastMouseMovePos[0], "(", 
        self.xScale.invert(self.lastMouseMovePos[0]), " keV)");
    }

    self.handleCancelMouseRecalibration();
  }

  // User is not right-click-and-dragging any more, so set this to null
  if (self.isRecalibrating)
    self.recalibrationMousePos = null;
}

SpectrumChartD3.prototype.handleCancelMouseRecalibration = function() {
  var self = this;

  var recalibrationStartLine = d3.select("#recalibrationStartLine"),
      recalibrationText = d3.select("#recalibrationText"),
      recalibrationMousePosLines = d3.select("#recalibrationMousePosLines");

  var recalibrationG = d3.select("#recalibrationG");
  var recalibrationPeakVis = d3.select("#recalibrationPeakVis");

  // Remove the right-click-and-drag initial starting point line
  recalibrationStartLine.remove();

  // Remove the right-click-and-drag text
  recalibrationText.remove()

  // Remove the peak vis
  recalibrationPeakVis.remove();

  // Remove all the arrow defs created
  d3.selectAll("#recalibrationArrowDef").each(function(){
    this.remove();
  });

  // Remove the right-click-and-drag end-point line
  if (self.recalibrationDistanceLines) {
    for (i = 0; i < self.recalibrationDistanceLines.length; i++) {
      self.recalibrationDistanceLines[i].remove();
    }
    self.recalibrationDistanceLines = null;
  }

  // Remove the right-click-and-drag mouse line
  recalibrationMousePosLines.remove();

  self.isRecalibrating = false;

  // // User is not right-click-and-dragging any more, so set this to null
  if (self.isRecalibrating)
    self.recalibrationMousePos = null;

  recalibrationG.remove();
}



/* Delete Peak functions */
SpectrumChartD3.prototype.handleMouseMoveDeletePeak = function() {
  var self = this;

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  d3.event.preventDefault();
  d3.event.stopPropagation();

  // Cancel the zooming mode
  self.handleCancelMouseZoomInX();

  // Cancel the recalibration mode
  self.handleCancelMouseRecalibration();

  // Cancel the zooming in y mode
  self.handleCancelMouseZoomInY();

  // Cancel the count gammas mode
  self.handleCancelMouseCountGammas();


  if (!self.deletePeaksMouse)
    return;

  // Adjust the mouse move position with respect to the bounds of the vis
  if (self.lastMouseMovePos[0] < 0)
    self.lastMouseMovePos[0] = 0;
  else if (self.lastMouseMovePos[0] > self.size.width)
    self.lastMouseMovePos[0] = self.size.width;

  // Create the erase-peaks range box and text 
  if (deletePeaksBox.empty()) {
    deletePeaksBox = self.vis.append("rect")
      .attr("id", "deletePeaksBox")
      .attr("class","deletePeaksBox")
      .attr("width", Math.abs( self.deletePeaksMouse[0] - self.lastMouseMovePos[0] ))
      .attr("height", self.size.height)
      .attr("x", self.lastMouseMovePos[0] < self.deletePeaksMouse[0] ? self.lastMouseMovePos[0] : self.deletePeaksMouse[0])
      .attr("y", 0);

    deletePeaksText = self.vis.append("text")
      .attr("id", "deletePeaksText")
      .attr("class", "deletePeaksText")
      .attr("y", Number(deletePeaksBox.attr("height"))/2)
      .text("Will Erase Peaks In Range");

  } else {  // Erase-peaks range box has already been created, update it

    // Mouse position to the left of initial starting point of erase peaks box
    if (self.lastMouseMovePos[0] < self.deletePeaksMouse[0])
      deletePeaksBox.attr("x", self.lastMouseMovePos[0])

    // Adjust the width of the erase peaks box
    deletePeaksBox.attr("width", Math.abs( self.deletePeaksMouse[0] - self.lastMouseMovePos[0] ));
  }

  // Move the erase peaks text in the middle of the erase peaks range box
  deletePeaksText.attr("x", Number(deletePeaksBox.attr("x")) + (Number(deletePeaksBox.attr("width"))/2) - 40 );
}

SpectrumChartD3.prototype.handleMouseUpDeletePeak = function() {
  var self = this;

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  var deletePeaksRange;

  try {
    deletePeaksRange = [ 
      Math.min(self.xScale.invert(Number(deletePeaksBox.attr("x"))), self.xScale.invert(Number(deletePeaksBox.attr("x")) + Number(deletePeaksBox.attr("width")))), 
      Math.max(self.xScale.invert(Number(deletePeaksBox.attr("x"))), self.xScale.invert(Number(deletePeaksBox.attr("x")) + Number(deletePeaksBox.attr("width")))) 
      ];

    console.log("Emit ERASE PEAKS SIGNAL FROM ", deletePeaksRange[0], "keV to ", deletePeaksRange[1], " keV" );
  } catch (TypeError) { // For some reason, a type error is (seldom) returned when trying to access "x" attribute of deletePeaksBox, doesn't affect overall functionality though
    return;
  }

  self.handleCancelMouseDeletePeak();
}

SpectrumChartD3.prototype.handleCancelMouseDeletePeak = function() {
  var self = this;

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  // Delete the erase peaks box since we are not erasing peaks anymore
  deletePeaksBox.remove();

  // Delete the erase peaks text since we are not erasing peaks anymore
  deletePeaksText.remove();

  // We are not erasing peaks anymore
  self.isDeletingPeaks = false;
}

SpectrumChartD3.prototype.handleTouchMoveDeletePeak = function() {
  var self = this;

  var t = self.deletePeaksTouches;

  if (t.length !== 2) {
    self.handleCancelTouchDeletePeak();
    return;
  }

  // Cancel the count gammas mode
  self.handleCancelTouchCountGammas();

  // Cancel the create peaks mode
  self.handleCancelTouchPeakFit();

  // Cancel the zoom-in y mode
  self.handleTouchCancelZoomInY();

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  var leftTouch = t[0][0] < t[1][0] ? t[0] : t[1],
      rightTouch = leftTouch === t[0] ? t[1] : t[0];

  // Create the erase-peaks range box and text 
  if (deletePeaksBox.empty()) {
    deletePeaksBox = self.vis.append("rect")
      .attr("id", "deletePeaksBox")
      .attr("class","deletePeaksBox")
      .attr("height", self.size.height)
      .attr("y", 0);

    deletePeaksText = self.vis.append("text")
      .attr("id", "deletePeaksText")
      .attr("class", "deletePeaksText")
      .attr("y", Number(deletePeaksBox.attr("height"))/2)
      .text("Will Erase Peaks In Range");

  }

  // Adjust the positioning and width of the delete peaks box
  deletePeaksBox.attr("x", leftTouch[0])
    .attr("width", Math.abs( rightTouch[0] - leftTouch[0] ));


  // Move the erase peaks text in the middle of the erase peaks range box
  deletePeaksText.attr("x", Number(deletePeaksBox.attr("x")) + (Number(deletePeaksBox.attr("width"))/2) - Number(deletePeaksText.node().getBBox().width)/2 );
}

SpectrumChartD3.prototype.handleTouchEndDeletePeak = function() {
  var self = this;

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  var deletePeaksRange;

  try {
    deletePeaksRange = [ 
      Math.min(self.xScale.invert(Number(deletePeaksBox.attr("x"))), self.xScale.invert(Number(deletePeaksBox.attr("x")) + Number(deletePeaksBox.attr("width")))), 
      Math.max(self.xScale.invert(Number(deletePeaksBox.attr("x"))), self.xScale.invert(Number(deletePeaksBox.attr("x")) + Number(deletePeaksBox.attr("width")))) 
      ];

    console.log("Emit ERASE PEAKS SIGNAL FROM ", deletePeaksRange[0], "keV to ", deletePeaksRange[1], " keV" );
  } catch (TypeError) { // For some reason, a type error is (seldom) returned when trying to access "x" attribute of deletePeaksBox, doesn't affect overall functionality though
    return;
  }

  self.handleCancelTouchDeletePeak();
}

SpectrumChartD3.prototype.handleCancelTouchDeletePeak = function() {
  var self = this;

  var deletePeaksBox = d3.select("#deletePeaksBox"),
      deletePeaksText = d3.select("#deletePeaksText");

  // Delete the erase peaks box since we are not erasing peaks anymore
  deletePeaksBox.remove();

  // Delete the erase peaks text since we are not erasing peaks anymore
  deletePeaksText.remove();

  // We are not deleting peaks anymore, delete the touches
  self.deletePeaksTouches = null;
}


/* Count Gammas functions */
SpectrumChartD3.prototype.handleMouseMoveCountGammas = function() {
  var self = this

  d3.event.preventDefault();
  d3.event.stopPropagation();

  // Cancel the zooming mode
  self.handleCancelMouseZoomInX();

  // Cancel the recalibration mode
  self.handleCancelMouseRecalibration();

  // Cancel the zooming in y mode
  self.handleCancelMouseZoomInY();

  // Cancel the delete peak mode
  self.handleCancelMouseDeletePeak();

  if (!self.countGammasMouse || !self.lastMouseMovePos || !self.rawData || !self.rawData.spectra || !self.rawData.spectra.length)
    return;

  function maxChannelSpectrum() {
    var result;
    self.rawData.spectra.forEach(function(spectrum) {
      if (!spectrum.points || !spectrum.points.length)
        return;
      if (!result || spectrum.points[spectrum.points.length-1].x > result.points[result.points.length-1].x)
        result = spectrum;
    });
    return result;
  }

  function gammaIntegral(hist, lowerX, upperX) {
    var sum = 0.0;
    var countGammasBox = d3.select("#countGammasBox");
    var spectraTitles = self.getSpectrumTitles();

    for (var i = 0; i < spectraTitles.length; ++i)
      spectraTitles[i] = spectraTitles[i].toUpperCase();

    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length || countGammasBox.empty() || hist == null || spectraTitles.indexOf(hist.toUpperCase()) == -1)
      return sum;

    var spectrum;
    self.rawData.spectra.forEach(function(spec) {
      if (spec.title.toUpperCase() == hist.toUpperCase()) {
        spectrum = spec;
      }
    });

    if (!spectrum)
      return sum;

    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];

    lowerX = Math.min( maxX, Math.max(lowerX, minX) );
    upperX = Math.max( minX, Math.min(upperX, maxX) );

    if (lowerX == upperX)
      return sum;

    if (lowerX > upperX) {  // swap the two values
      var tmp = lowerX;
      lowerX = upperX;
      upperX = tmp;
    }

    var maxChannel = spectrum.points.length - 1;
    var lowerChannel = d3.bisector(function(d){return d.x;}).left(spectrum.points,lowerX,1) - 1;
    var upperChannel = d3.bisector(function(d){return d.x;}).left(spectrum.points,upperX,1) - 1;

    var lowerLowEdge = spectrum.points[lowerChannel].x;
    var lowerBinWidth = lowerChannel < spectrum.points.length - 1 ? spectrum.points[lowerChannel+1].x - spectrum.points[lowerChannel].x : 
                                                                spectrum.points[lowerChannel].x - spectrum.points[lowerChannel-1].x;
    var lowerUpEdge = lowerLowEdge + lowerBinWidth;

    if (lowerChannel === upperChannel) {
      var frac = (upperX - lowerX) / lowerBinWidth;
      console.log("lowerChannel == upper channel, counts = ", frac * spectrum.points[lowerChannel].y);
      return frac * spectrum.points[lowerChannel].y;
    }

    var fracLowBin = (lowerUpEdge - lowerX) / lowerBinWidth;
    sum += fracLowBin * spectrum.points[lowerChannel].y;

    var upperLowEdge = spectrum.points[upperChannel].x;
    var upperBinWidth = upperChannel < spectrum.points.length - 1 ? spectrum.points[upperChannel+1].x - spectrum.points[upperChannel].x : 
                                                                spectrum.points[upperChannel].x - spectrum.points[upperChannel-1].x;

    var fracUpBin = (upperX - upperLowEdge) / upperBinWidth;
    sum += fracUpBin * spectrum.points[upperChannel].y;

    for (var channel = lowerChannel + 1; channel < upperChannel; channel++) {
      sum += spectrum.points[channel].y;
    }

    return sum;
  }
  function needsDecimal(num) {
    return num % 1 != 0;
  }

  // Adjust the mouse move position with respect to the bounds of the vis
  if (self.lastMouseMovePos[0] < 0)
    self.lastMouseMovePos[0] = 0;
  else if (self.lastMouseMovePos[0] > self.size.width)
    self.lastMouseMovePos[0] = self.size.width;

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText"),
      countGammaRangeText = d3.select("#countGammaRangeText"),
      sigmaCount = d3.select("#sigmaCount");
  var countGammasRange = [ 
      Number(self.xScale.invert(Number(self.countGammasMouse[0]).toFixed(1))), 
      Number(self.xScale.invert(Number(self.lastMouseMovePos[0]).toFixed(1)))
  ];

  // Create the yellow box and text associated with it
  if (countGammasBox.empty()) { 
    countGammasBox = self.vis.append("rect")
      .attr("id", "countGammasBox")
      .attr("width", Math.abs( self.countGammasMouse[0] - self.lastMouseMovePos[0] ))
      .attr("height", self.size.height)
      .attr("x", self.lastMouseMovePos[0] < self.countGammasMouse[0] ? self.lastMouseMovePos[0] : self.countGammasMouse[0])
      .attr("y", 0);

    countGammasText = self.vis.append("text")
      .attr("id", "countGammasText")
      .attr("class", "countGammasText")
      .attr("y", Number(countGammasBox.attr("height"))/2)
      .text("Gamma Counts");

  } else {  // Adjust the width of the erase peaks box
    countGammasBox.attr("width", Math.abs( self.countGammasMouse[0] - self.lastMouseMovePos[0] ));
  }
  var ypos = Number(countGammasBox.attr("height"))/2 + 15;   // signifies the y-position of the text displayed

  // Mouse position to the left of initial starting point of count gammas box
  if (self.lastMouseMovePos[0] < self.countGammasMouse[0]) {
    countGammasBox.attr("class", "deletePeaksBox")
      .attr("x", self.lastMouseMovePos[0]);
    countGammasText.text("Will remove gamma count range");

    // Move the count gammas text in the middle of the count gammas box
    countGammasText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) - /*30*/ Number(countGammasText.node().getBBox().width)/2 );

    // Remove the counts text (we will not be counting gammas)
    self.getSpectrumTitles().forEach(function(spectrumTitle, i) {
      d3.select("#"+spectrumTitle.replace(' ','-')+"CountsText").remove();
      d3.select("#"+spectrumTitle.replace(' ','-')+"AsterickText").remove();
      d3.select("#asterickText"+i).remove();
      countGammaRangeText.remove();
    });
    sigmaCount.remove();
    return;

  } else {
    countGammasBox.attr("class", "countGammasBoxForward")
    countGammasText.text("Gamma Counts");

    if (countGammaRangeText.empty())
      countGammaRangeText = self.vis.append("text")
        .attr("id", "countGammaRangeText")
        .attr("class", "countGammasText")
        .attr("y", ypos);


    countGammaRangeText.text((needsDecimal(countGammasRange[0]) ? countGammasRange[0].toFixed(1) : countGammasRange[0].toFixed()) + 
                              " to " + 
                              (needsDecimal(countGammasRange[1]) ? countGammasRange[1].toFixed(1) : countGammasRange[1].toFixed()) + " keV");
    countGammaRangeText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) - /*30*/ Number(countGammaRangeText.node().getBBox().width)/2);
  }
  ypos += 15;

  // Move the count gammas text in the middle of the count gammas box
  countGammasText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) - /*30*/ Number(countGammasText.node().getBBox().width)/2 );

  // Display the count gammas text for all the spectrum
  var nforeground, nbackground;
  var scaleBackground, backgroundScaleFactor, backgroundSigma, backgroundScaleSigma;
  var nsigma = 0, isneg;
  var asterickText = "";
  var rightPadding = 50;
  var specialScaleSpectras = [];
  self.rawData.spectra.forEach(function(spectrum, i) {   // TODO: May need to change processing of
    if (!spectrum || spectrum.title == null)
      return;

    // Get information from the spectrum
    var spectrumCountsText = d3.select("#" + spectrum.title.replace(' ','-') + "CountsText");
    var spectrumScaleFactor = spectrum.yScaleFactor;
    var nspectrum = gammaIntegral(spectrum.title, countGammasRange[0], countGammasRange[1]);
    var spectrumGammaCount = Number((spectrumScaleFactor * nspectrum).toFixed(2));
    var countsText;

    // Save information for the foreground and background (for sigma comparison)
    if (spectrum.title.toUpperCase == "FOREGROUND" || i == 0)
      nforeground = nspectrum;
    else if (spectrum.title.toUpperCase() == "BACKGROUND") {
      nbackground = nspectrum;
      backgroundScaleFactor = spectrumScaleFactor;
      scaleBackground = nbackground * backgroundScaleFactor;
    }

    // Get the text to be displayed from the spectrum information
    if (spectrumScaleFactor != null && spectrumScaleFactor !== -1)
      countsText = spectrum.title + ": " + (needsDecimal(spectrumGammaCount) ? spectrumGammaCount.toFixed(2) : spectrumGammaCount.toFixed());
    if (spectrumScaleFactor != 1) {
      asterickText += "*";
      if (countsText)
        countsText += asterickText;
      specialScaleSpectras.push(asterickText + "scaled by " + 
        (needsDecimal(spectrumScaleFactor) ? spectrumScaleFactor.toFixed(3) : spectrumScaleFactor.toFixed()) + " from actual");
    }

    // Output the count gammas information to the chart
    if (countsText) {
      if (spectrumCountsText.empty())
        spectrumCountsText = self.vis.append("text")
          .attr("id", spectrum.title.replace(' ','-') + "CountsText")
          .attr("class", "countGammasText")
          .attr("y", ypos);
      spectrumCountsText.text(countsText); 
      spectrumCountsText.attr("x", Number(countGammasText.attr("x")) - rightPadding );
      ypos += 15;

    } else {
      spectrumCountsText.remove();
    }
  });

  // Get proper information for foreground-background sigma comparison
  if (nforeground && nbackground && backgroundScaleFactor) {
    backgroundSigma = Math.sqrt(nbackground);
    backgroundScaleSigma = backgroundScaleFactor * backgroundSigma;
    nsigma = backgroundScaleSigma == 0 ? 0 : (Number((Math.abs(nforeground - scaleBackground) / backgroundScaleSigma).toFixed(3)));
    isneg = scaleBackground > nforeground;
  }

  // Output foreground-background sigma information if it is available
  if (nsigma > 0) {
    if (sigmaCount.empty())
      sigmaCount = self.vis.append("text")
        .attr("id", "sigmaCount")
        .attr("class", "countGammasText")
        .attr("y", ypos);

    sigmaCount.attr("x", Number(countGammasText.attr("x")) - rightPadding + 10)
      .text("Foreground is " + (needsDecimal(nsigma) ? nsigma.toFixed(2) : nsigma.toFixed() ) + " σ " + (isneg ? "below" : "above") + " background.");
    ypos += 15;

  } else if (!sigmaCount.empty()) {
      sigmaCount.remove();
  }

  // Output all the corresponding asterick text with each spectrum (if scale factor != 1)
  specialScaleSpectras.forEach(function(string, i) {
    if (string == null || !string.length)
      return;
    var stringnode = d3.select("#asterickText" + i);

    if (stringnode.empty())
      stringnode = self.vis.append("text")
        .attr("id", "asterickText"+i)
        .attr("class", "countGammasText")
        .text(string);
    stringnode.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) + rightPadding/3)
      .attr("y", ypos);
    ypos += 15;
  });

  // Hide all the count gamma text when the mouse box is empty
  if (Number(countGammasBox.attr("width")) == 0)
    d3.selectAll(".countGammasText").attr("fill-opacity", 0);
  else
    d3.selectAll(".countGammasText").attr("fill-opacity", 1);
}

SpectrumChartD3.prototype.handleMouseUpCountGammas = function() {
  var self = this;

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText");

  var countGammasRange;

  try {
    countGammasRange = [ 
      Math.min(self.xScale.invert(Number(countGammasBox.attr("x"))), self.xScale.invert(Number(countGammasBox.attr("x")) + Number(countGammasBox.attr("width")))), 
      Math.max(self.xScale.invert(Number(countGammasBox.attr("x"))), self.xScale.invert(Number(countGammasBox.attr("x")) + Number(countGammasBox.attr("width")))) 
      ];

    if (self.lastMouseMovePos[0] < self.countGammasMouse[0]) {
      console.log("Emit REMOVE GAMMA COUNT SIGNAL FROM ", countGammasRange[0], "keV to ", countGammasRange[1], " keV" );
    } else {
      console.log("Emit COUNT GAMMAS SIGNAL FROM ", countGammasRange[0], "keV to ", countGammasRange[1], " keV" );
    }

    
  } catch (TypeError) { // For some reason, a type error is (seldom) returned when trying to access "x" attribute of countGammasBox, doesn't affect overall functionality though
    return;
  }

  self.handleCancelMouseCountGammas();
}

SpectrumChartD3.prototype.handleCancelMouseCountGammas = function() {
  var self = this;

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText");

  // Delete the count gammas box since we are not counting gammas anymore
  countGammasBox.remove();

  // Delete the count gamma texts since we are not counting gammas anymore
  self.vis.selectAll(".countGammasText").remove()

  // We are not erasing peaks anymore
  self.isCountingGammas = false;
}

SpectrumChartD3.prototype.handleTouchMoveCountGammas = function() {
  var self = this;

  // Cancel delete peaks mode
  self.handleCancelTouchDeletePeak();

  // Cancel the create peaks mode
  self.handleCancelTouchPeakFit();

  // Cancel the zoom-in y mode
  self.handleTouchCancelZoomInY();


  var t = d3.touches(self.vis[0][0]);

  if (t.length !== 2 || !self.countGammasStartTouches) {
    self.handleCancelTouchCountGammas();
    return;
  }

  function gammaIntegral(hist, lowerX, upperX) {
    var sum = 0.0;
    var countGammasBox = d3.select("#countGammasBox");
    var spectraTitles = self.getSpectrumTitles();

    for (var i = 0; i < spectraTitles.length; ++i)
      spectraTitles[i] = spectraTitles[i].toUpperCase();

    if (!self.rawData || !self.rawData.spectra || !self.rawData.spectra.length || !self.rawData.spectra[0].x || countGammasBox.empty() || hist == null || spectraTitles.indexOf(hist.toUpperCase()) == -1)
      return sum;

    var spectrum;
    self.rawData.spectra.forEach(function(spec) {
      if (spec.title.toUpperCase() == hist.toUpperCase()) {
        spectrum = spec;
      }
    });

    if (!spectrum)
      return sum;

    var y = 'y';
    var bounds = self.min_max_x_values();
    var maxX = bounds[1];
    var minX = bounds[0];

    lowerX = Math.min( maxX, Math.max(lowerX, minX) );
    upperX = Math.max( minX, Math.min(upperX, maxX) );

    if (lowerX == upperX)
      return sum;

    if (lowerX > upperX) {  // swap the two values
      var tmp = lowerX;
      lowerX = upperX;
      upperX = tmp;
    }

    var maxChannel = spectrum.points.length - 1;
    var lowerChannel = d3.bisector(function(d){return d.x;}).left(spectrum.points,lowerX,1) - 1;
    var upperChannel = d3.bisector(function(d){return d.x;}).left(spectrum.points,upperX,1) - 1;

    var lowerLowEdge = spectrum.points[lowerChannel].x;
    var lowerBinWidth = lowerChannel < spectrum.points.length - 1 ? spectrum.points[lowerChannel+1].x - spectrum.points[lowerChannel].x : 
                                                                spectrum.points[lowerChannel].x - spectrum.points[lowerChannel-1].x;
    var lowerUpEdge = lowerLowEdge + lowerBinWidth;

    if (lowerChannel === upperChannel) {
      var frac = (upperX - lowerX) / lowerBinWidth;
      console.log("lowerChannel == upper channel, counts = ", frac * spectrum.points[lowerChannel][y]);
      return frac * spectrum.points[lowerChannel][y];
    }

    var fracLowBin = (lowerUpEdge - lowerX) / lowerBinWidth;
    sum += fracLowBin * spectrum.points[lowerChannel][y];

    var upperLowEdge = spectrum.points[upperChannel].x;
    var upperBinWidth = upperChannel < spectrum.points.length - 1 ? spectrum.points[upperChannel+1].x - spectrum.points[upperChannel].x : 
                                                                spectrum.points[upperChannel].x - spectrum.points[upperChannel-1].x;

    var fracUpBin = (upperX - upperLowEdge) / upperBinWidth;
    sum += fracUpBin * spectrum.points[upperChannel][y];

    for (var channel = lowerChannel + 1; channel < upperChannel; channel++) {
      sum += spectrum.points[channel][y];
    }

    return sum;
  }
  function needsDecimal(num) {
    return num % 1 != 0;
  }

  var leftStartTouch = self.countGammasStartTouches[0][0] < self.countGammasStartTouches[1][0] ? self.countGammasStartTouches[0] : self.countGammasStartTouches[1],
      rightStartTouch = leftStartTouch === self.countGammasStartTouches[0] ? self.countGammasStartTouches[1] : self.countGammasStartTouches[0];

  var leftTouch = t[0][0] < t[1][0] ? t[0] : t[1],
      rightTouch = leftTouch === t[0] ? t[1] : t[0];

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText"),
      countGammaRangeText = d3.select("#countGammaRangeText"),
      foregroundCountsText = d3.select("#foregroundCountsText"),
      backgroundCountsText = d3.select("#backgroundCountsText"),
      secondaryCountsText = d3.select("#secondaryCountsText")
      foregroundAsterickText = d3.select("#foregroundAsterickText"),
      backgroundAsterickText = d3.select("#backgroundAsterickText"),
      secondaryAsterickText = d3.select("#secondaryAsterickText"),
      sigmaCount = d3.select("#sigmaCount");

  var countGammasRange = [ 
      self.xScale.invert(leftStartTouch[0]), 
      self.xScale.invert(rightTouch[0])
  ];

  // Create the yellow box associated with counting gammas
  if (countGammasBox.empty()) {
    countGammasBox = self.vis.append("rect")
      .attr("id", "countGammasBox")
      .attr("class", "countGammasBoxForward")
      .attr("width", Math.abs( leftStartTouch[0] - rightTouch[0] ))
      .attr("height", self.size.height)
      .attr("x", leftStartTouch[0])
      .attr("y", 0);

    countGammasText = self.vis.append("text")
      .attr("id", "countGammasText")
      .attr("class", "countGammasText")
      .attr("y", Number(countGammasBox.attr("height"))/2)
      .text("Count Gammas In Range");

  } else {
    // Adjust the width of the erase peaks box
    countGammasBox.attr("width", Math.abs( leftStartTouch[0] - rightTouch[0] ));
  }

  // Move the count gammas text in the middle of the count gammas box
  countGammasText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) - Number(countGammasText.node().getBBox().width)/2 );

  /*
  IMPORTANT: How raw data is structured: each attribute of raw data has 3 elements in the array which symbolize the data sets for
             the foreground, background, and secondary foreground respectively.
  
  ******* Scale factor in the form: [foregroundScaleFactor, backgroundScaleFactor, secondForegroundScaleFactor] ***********
  */
  var foregroundScaleFactor = self.rawData.spectra[0].yScaleFactor;;
  var backgroundScaleFactor = self.rawData.spectra[1].yScaleFactor;;
  var secondaryScaleFactor = self.rawData.spectra[2].yScaleFactor;;

  var nforeground = gammaIntegral("foreground", countGammasRange[0], countGammasRange[1]);
  var nbackground = gammaIntegral("background", countGammasRange[0], countGammasRange[1]);
  var nsecondary = gammaIntegral("secondary", countGammasRange[0], countGammasRange[1]);

  var foregroundGammaCount = Number((foregroundScaleFactor * nforeground).toFixed(2));
  var backgroundGammaCount = Number((backgroundScaleFactor * nbackground).toFixed(2));
  var secondaryGammaCount = Number((secondaryScaleFactor * nsecondary).toFixed(2));

  var scaleBackground = nbackground * backgroundScaleFactor;
  var backgroundSigma = Math.sqrt(nbackground);
  var backgroundScaleSigma = backgroundScaleFactor * backgroundSigma;
  var nsigma = backgroundScaleSigma == 0 ? 0 : (Number((Math.abs(nforeground - scaleBackground) / backgroundScaleSigma).toFixed(3)));
  var isneg = scaleBackground > nforeground;

  var foregroundText = null, backgroundText = null, secondaryText = null;

  var foregroundAsterick = foregroundScaleFactor != 1 ? "*" : "";
  var backgroundAsterick = (backgroundScaleFactor == 1 ? "" : foregroundScaleFactor != 1 ? "**" : "*");
  var secondaryAsterick = (secondaryScaleFactor == 1 ? "" : foregroundScaleFactor != 1 && backgroundScaleFactor != 1 ? "***" : foregroundScaleFactor == 1 && backgroundScaleFactor == 1 ? "*" : "**");

  var rightPadding = 50;

  if (foregroundScaleFactor != null && foregroundScaleFactor !== -1) {
    foregroundText = "Foreground: " + (needsDecimal(foregroundGammaCount) ? foregroundGammaCount.toFixed(2) : foregroundGammaCount.toFixed());
  }
  if (backgroundScaleFactor != null && backgroundScaleFactor !== -1) {
    backgroundText = "Background: " + (needsDecimal(backgroundGammaCount) ? backgroundGammaCount.toFixed(2) : backgroundGammaCount.toFixed());              
  }
  if (secondaryScaleFactor != null && secondaryScaleFactor !== -1) {
    secondaryText = "Secondary: " + (needsDecimal(secondaryGammaCount) ? secondaryGammaCount.toFixed(2) : secondaryGammaCount.toFixed());                 
  }

  if (foregroundText) {
    if (foregroundCountsText.empty()) {
      foregroundCountsText = self.vis.append("text")
        .attr("id", "foregroundCountsText")
        .attr("class", "countGammasText")
        .attr("y", Number(countGammasBox.attr("height"))/2 + 30);
    }
    foregroundCountsText.text(foregroundText + foregroundAsterick); 
    foregroundCountsText.attr("x", Number(countGammasText.attr("x")) - rightPadding );

  } else {
    foregroundCountsText.remove();
  }

  if (backgroundText) {
    if (backgroundCountsText.empty()) {
      backgroundCountsText = self.vis.append("text")
        .attr("id", "backgroundCountsText")
        .attr("class", "countGammasText")
        .attr("y", Number(countGammasBox.attr("height"))/2 + (!foregroundCountsText.empty() ? 45 : 30));  // Adjust height if foreground text not present
    }
    backgroundCountsText.text(backgroundText + backgroundAsterick); 
    backgroundCountsText.attr("x", Number(countGammasText.attr("x")) - rightPadding );

  } else {
    backgroundCountsText.remove();
  }

  if (secondaryText) {
    if (secondaryCountsText.empty()) {
      secondaryCountsText = self.vis.append("text")
        .attr("id", "secondaryCountsText")
        .attr("class", "countGammasText")                     // Adjust height if foreground/background text not present
        .attr("y", Number(countGammasBox.attr("height"))/2 + (!foregroundCountsText.empty() && !backgroundCountsText.empty() ? 60 : 
                                                              foregroundCountsText.empty() && backgroundCountsText.empty() ? 15 : 30));  
    }
    secondaryCountsText.text(secondaryText + secondaryAsterick); 
    secondaryCountsText.attr("x", Number(countGammasText.attr("x")) - rightPadding );

  } else {
    secondaryCountsText.remove();
  }

  if (nsigma > 0) {
    if (sigmaCount.empty()) 
      sigmaCount = self.vis.append("text")
        .attr("id", "sigmaCount")
        .attr("class", "countGammasText")
        .attr("y", (!secondaryCountsText.empty() ? Number(secondaryCountsText.attr("y")) + 15 : 
                    !backgroundCountsText.empty() ? Number(backgroundCountsText.attr("y")) + 15 :
                    Number(foregroundCountsText.attr("y")) + 15));  

    sigmaCount.attr("x", Number(countGammasText.attr("x")) - rightPadding + 10)
      .text("Foreground is " + (needsDecimal(nsigma) ? nsigma.toFixed(2) : nsigma.toFixed() ) + " σ " + (isneg ? "below" : "above") + " background.");

  } else if (!sigmaCount.empty()) {
      sigmaCount.remove();
  }

  // Add additional scale factor information for foreground, background, secondary
  if (!foregroundCountsText.empty() && foregroundAsterick != "" && foregroundAsterickText.empty()) {
    foregroundAsterickText = self.vis.append("text")
      .attr("id", "foregroundAsterickText")
      .attr("class", "countGammasText")
      .text(foregroundAsterick + "scaled by " + (needsDecimal(foregroundScaleFactor) ? foregroundScaleFactor.toFixed(3) : foregroundScaleFactor.toFixed()) + " from actual");
  }
  if (!backgroundCountsText.empty() && backgroundAsterick != "" && backgroundAsterickText.empty()) {
    backgroundAsterickText = self.vis.append("text")
      .attr("id", "backgroundAsterickText")
      .attr("class", "countGammasText")
      .text(backgroundAsterick + "scaled by " + (needsDecimal(backgroundScaleFactor) ? backgroundScaleFactor.toFixed(3) : backgroundScaleFactor.toFixed()) + " from actual");  
  }
  if (!secondaryCountsText.empty() && secondaryAsterick != "" && secondaryAsterickText.empty()) {
    secondaryAsterickText = self.vis.append("text")
      .attr("id", "secondaryAsterickText")
      .attr("class", "countGammasText")
      .text(secondaryAsterick + "scaled by " + (needsDecimal(secondaryScaleFactor) ? secondaryScaleFactor.toFixed(3) : secondaryScaleFactor.toFixed()) + " from actual");  
  }

  if (!foregroundAsterickText.empty())
    foregroundAsterickText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) + rightPadding/3)
      .attr("y", (!sigmaCount.empty() ? Number(sigmaCount.attr("y")) + 15 :
                  !secondaryCountsText.empty() ? Number(secondaryCountsText.attr("y")) + 15 : 
                  !backgroundCountsText.empty() ? Number(backgroundCountsText.attr("y")) + 15 :
                  Number(foregroundCountsText.attr("y")) + 15));  

  if (!backgroundAsterickText.empty())
    backgroundAsterickText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) + rightPadding/3)
      .attr("y", (!foregroundAsterickText.empty() ? Number(foregroundAsterickText.attr("y")) + 15 :
                  !sigmaCount.empty() ? Number(sigmaCount.attr("y")) + 15 :
                  !secondaryCountsText.empty() ? Number(secondaryCountsText.attr("y")) + 15 : 
                  !backgroundCountsText.empty() ? Number(backgroundCountsText.attr("y")) + 15 :
                  Number(foregroundCountsText.attr("y")) + 15));

  if (!secondaryAsterickText.empty())
    secondaryAsterickText.attr("x", Number(countGammasBox.attr("x")) + (Number(countGammasBox.attr("width"))/2) + rightPadding/3)
      .attr("y", (!backgroundAsterickText.empty() ? Number(backgroundAsterickText.attr("y")) + 15 :
                  !foregroundAsterickText.empty() ? Number(foregroundAsterickText.attr("y")) + 15 :
                  !sigmaCount.empty() ? Number(sigmaCount.attr("y")) + 15 :
                  !secondaryCountsText.empty() ? Number(secondaryCountsText.attr("y")) + 15 : 
                  !backgroundCountsText.empty() ? Number(backgroundCountsText.attr("y")) + 15 :
                  Number(foregroundCountsText.attr("y")) + 15));
}

SpectrumChartD3.prototype.handleTouchEndCountGammas = function() {
  var self = this;

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText");

  var countGammasRange;

  try {
    countGammasRange = [ 
      Math.min(self.xScale.invert(Number(countGammasBox.attr("x"))), self.xScale.invert(Number(countGammasBox.attr("x")) + Number(countGammasBox.attr("width")))), 
      Math.max(self.xScale.invert(Number(countGammasBox.attr("x"))), self.xScale.invert(Number(countGammasBox.attr("x")) + Number(countGammasBox.attr("width")))) 
      ];

    console.log("Emit COUNT GAMMAS SIGNAL FROM ", countGammasRange[0], "keV to ", countGammasRange[1], " keV" );

  } catch (TypeError) { // For some reason, a type error is (seldom) returned when trying to access "x" attribute of countGammasBox, doesn't affect overall functionality though
    return;
  }

  self.handleCancelTouchCountGammas();
}

SpectrumChartD3.prototype.handleCancelTouchCountGammas = function() {
  var self = this;

  var countGammasBox = d3.select("#countGammasBox"),
      countGammasText = d3.select("#countGammasText"),
      foregroundCountsText = d3.select("#foregroundCountsText"),
      backgroundCountsText = d3.select("#backgroundCountsText"),
      secondaryCountsText = d3.select("#secondaryCountsText");

  // Delete the count gammas box since we are not counting gammas anymore
  countGammasBox.remove();

  // Delete the count gamma texts since we are not counting gammas anymore
  d3.selectAll(".countGammasText").forEach(function (texts) {
    texts.forEach(function(text) {
      text.remove();
    })
  });
}





