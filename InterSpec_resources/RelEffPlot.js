/* Relative-efficiency chart.

 Shared by BOTH relative-activity tools - "Isotopics from peaks" (RelActManualGui) and
 "Isotopics by nuclides" (RelActAutoGui) - via RelEffChart, and it is also inlined into the HTML
 reports RelActAutoReport exports.  A change here therefore lands in all four places at once.

 Uncertainty convention: the shaded band around a fitted curve and the data-point error bars are
 both +-1 sigma (the band from the fitted curve's covariance, the bars statistical only).  Both
 tools' help files state this; keep them in step if it ever changes.
*/
RelEffPlot = function (elem,options) {
  this.chart = typeof elem === "string" ? document.getElementById(elem) : elem;

  this.options = options ? options : {};
  if( !this.options.margins )
    this.options.margins = {};
  if( typeof this.options.margins.top !== "number" )
    this.options.margins.top = 5;
  if( typeof this.options.margins.right !== "number" )
    this.options.margins.right = 5;
  if( typeof this.options.margins.bottom !== "number" )
    this.options.margins.bottom = 5;
  if( typeof this.options.margins.left !== "number" )
    this.options.margins.left = 5;
  if( typeof this.options.titleToAxisPadding !== "number" )
    this.options.titleToAxisPadding = 0;
  if( (typeof this.options.xAxisTitle !== "string") || (this.options.xAxisTitle.length === 0) )
    this.options.xAxisTitle = null;
  if( (typeof this.options.yAxisTitle !== "string") || (this.options.yAxisTitle.length === 0) )
    this.options.yAxisTitle = null;
      
  // Initialize custom colors with defaults
  this.customColors = ["#cfced2", "#a9ebdd", "#fd8273"];
      
  // Set the dimensions of the canvas / graph.
  // The container may not yet have been laid out (zero or sub-margin size) when the
  //  constructor runs - clamp so we never feed negative values to SVG attributes.
  //  A subsequent handleResize / setRelEffData call will set the real geometry.
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  const chartAreaWidth = Math.max( 0, parentWidth - this.options.margins.left - this.options.margins.right );
  const chartAreaHeight = Math.max( 0, parentHeight - this.options.margins.top - this.options.margins.bottom );

  // Setup the tooltip
  this.tooltip = d3.select(this.chart).append("div")
    .attr("class", "RelEffPlotTooltip")
    .style("opacity", 0);

  // Set the ranges
  this.xScale = d3.scale.linear().range([0, chartAreaWidth]);
  this.yScale = d3.scale.linear().range([chartAreaHeight, 0]);

  // Define the axes
  this.xAxis = d3.svg
    .axis()
    .scale(this.xScale)
    .orient("bottom")
    .tickFormat( x => `${x.toFixed(0)}` )
    .ticks(5);

  const self = this;
  this.yAxis = d3.svg
    .axis()
    .scale(this.yScale)
    .tickFormat( function(x){ return self.yTickLabel(x); } )
    .orient("left")
    .ticks(5);

  // Adds the svg canvas
  this.svg = d3.select(this.chart)
    .append("svg")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("class", "RelEffPlot");
  this.chartArea = this.svg  //This holds everything except the optional x and y axis titles
    .append("g")
    .attr("transform", "translate(" + this.options.margins.left + "," + this.options.margins.top + ")");
    
  // Create clipping path to prevent elements from being drawn outside the chart area
  this.chartArea.append("defs")
    .append("clipPath")
    .attr("id", "releffplot-clip-" + this.chart.id)
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", chartAreaWidth)
    .attr("height", chartAreaHeight);
    
  this.xScale.domain([0, 3000]);
  this.yScale.domain([0, 1]);

  // Setup mouse interactions for zoom functionality (before plotGroup so it's behind data points)
  this.setupMouseInteractions();

  // Create a clipped group for all plot elements (paths, circles, error bars)
  this.plotGroup = this.chartArea.append("g")
    .attr("clip-path", "url(#releffplot-clip-" + this.chart.id + ")");

  this.path_uncert = this.plotGroup.append("path")
    .attr("class", "RelEffPlotErrorBounds") //Applying in CSS doesnt seem to work for some reason
    .style("fill", "rgba(0, 0, 255, 0.1)")
    .style("stroke", "none")
    .style("pointer-events", "none")
    .style("stroke", "none")
    .style("pointer-events", "none");

  // Add the valueline path.
  this.path = this.plotGroup.append("path")    // Add the valueline path.
    .attr("class", "line")
    .style("pointer-events", "none");

  // Add the X Axis
  this.chartArea.append("g")
      .attr("class", "xAxis")
      .attr("transform", "translate(0," + chartAreaHeight + ")")
      .call(this.xAxis);

  // Add the Y Axis
  this.chartArea.append("g")
    .attr("class", "yAxis")
    .call(this.yAxis);

  if( this.options.xAxisTitle )
    this.setXAxisTitle( this.options.xAxisTitle, true );

  if( this.options.yAxisTitle )
    this.setYAxisTitle( this.options.yAxisTitle, true );

  // Initialize mouse interaction variables for zoom
  this.leftMouseDown = null;
  this.zoomAxis = null;  // 'x' or 'y' while a drag-zoom is in progress
  this.zooming = false;
}//RelEffPlot constructor


/* The y-axis tick label for value `x`.  Uses fixed notation with just enough decimals to tell adjacent ticks
   apart, unless compact scientific notation (written like the spectrum chart, e.g. "4e-5", "1e+4") makes the
   widest label on the axis shorter - the whole axis uses the same notation.  (We dont use d3's
   scale.tickFormat(), as it adds thousands separators.)
 */
RelEffPlot.prototype.yTickLabel = function( x ) {
  if( x == 0 )
    return "0";

  // d3 tick steps are 1, 2, or 5 x 10^n; the +0.01 guards against log10(1000) = 2.9999...
  const ticks = this.yScale.ticks( 5 );
  const step = (ticks.length > 1) ? Math.abs( ticks[1] - ticks[0] ) : Math.abs( x );
  if( !(step > 0) || !isFinite(step) )
    return x.toFixed( 2 );
  const stepExp = Math.floor( Math.log10(step) + 0.01 );

  const fixed = function( v ){
    return (v == 0) ? "0" : v.toFixed( Math.min( 12, Math.max( 0, -stepExp ) ) );
  };

  // Mantissa digits down to the decade of the tick step
  const sci = function( v ){
    if( v == 0 )
      return "0";
    const ndec = Math.floor( Math.log10( Math.abs(v) ) + 1.0E-9 ) - stepExp;
    return v.toExponential( Math.min( 12, Math.max( 0, ndec ) ) );
  };

  const widest = function( fmt ){
    return d3.max( ticks, function(t){ return fmt(t).length; } );
  };

  return (widest(sci) < widest(fixed)) ? sci(x) : fixed(x);
};//RelEffPlot.prototype.yTickLabel


RelEffPlot.prototype.setupMouseInteractions = function() {
  const self = this;

  // Add mouse interaction rect to the chart area - will be positioned later in layoutAndDraw
  this.mouseCapture = this.chartArea.append("rect")
    .attr("class", "mouse-capture")
    .style("fill", "none")
    .style("pointer-events", "all")
    .on("mousedown", function() {
      // Only a plain left-button drag starts a zoom; right/middle clicks (and ctrl-click on macOS,
      //  which opens a context menu whose mouseup we never see) would otherwise leave the zoom
      //  state machine stuck half-started.
      if (d3.event.button !== 0)
        return;

      // Remove any zoom box orphaned by a mouseup we never received (released outside the
      //  webview, context menu, etc.) - it sits above the data points and swallows events.
      self.chartArea.selectAll(".zoom-box, .zoom-text").remove();

      const mouse = d3.mouse(this);
      self.leftMouseDown = mouse;
      self.zoomStartX = self.xScale.invert(mouse[0]);

      // Like the spectrum chart, Cmd-drag zooms the y-axis; Alt-drag does too, as Windows/Linux have no
      //  usable Cmd key.  Captured once, so changing modifier keys mid-drag has no effect.
      self.zoomAxis = (d3.event.metaKey || d3.event.altKey) ? 'y' : 'x';
      d3.event.preventDefault();

      // Add document listeners for mouse move and up
      d3.select(document)
        .on("mousemove.releffzoom", function() { self.handleMouseMove(); })
        .on("mouseup.releffzoom", function() { self.handleMouseUp(); });
    })
    .on("mousemove", function() {
      if (!self.leftMouseDown) {
        // Could add tooltip functionality here in the future
      }
    });

  // On the whole SVG, so wheeling over a data point (which sits above the capture rect) works too
  this.svg.on("wheel.releffzoom", function() { self.handleWheel(); });
};//RelEffPlot.prototype.setupMouseInteractions


/* Like the spectrum chart, the mouse wheel zooms the x-axis about the energy under the mouse (wheel down
   zooms out), and a sideways scroll pans; the view stays within the data's energy range. */
RelEffPlot.prototype.handleWheel = function() {
  const e = d3.event;

  // Modifier keys are left for other uses, as in the spectrum chart (ctrl is a trackpad pinch, so is used)
  if( e.altKey || e.shiftKey || e.metaKey || this.leftMouseDown || !this.datasets )
    return;

  // Outside the plot area (axes, titles), leave the wheel to scroll the page
  const m = d3.mouse( this.chartArea.node() );
  const width = this.xScale.range()[1], height = this.yScale.range()[0];
  if( !(width > 0) || (m[0] < 0) || (m[0] > width) || (m[1] < 0) || (m[1] > height) )
    return;
  e.preventDefault();

  // Line or page based deltas (e.g., Firefox with a mouse wheel) to pixels
  const pxPerUnit = (e.deltaMode === 1) ? 16 : ((e.deltaMode === 2) ? height : 1);
  const dx = pxPerUnit * e.deltaX, dy = pxPerUnit * e.deltaY;

  const full = this.fullXDomain();
  const fullRange = full[1] - full[0];
  const domain = this.xScale.domain();
  const range = domain[1] - domain[0];

  // Pinch deltas are much smaller than wheel deltas, so zoom faster for them
  const zoomRate = e.ctrlKey ? 0.01 : 0.002;
  const newRange = Math.max( Math.min( 1, fullRange ), Math.min( fullRange, range * Math.exp( zoomRate * dy ) ) );

  // Keep the mouse over the same energy, then pan by the sideways scroll distance
  let lower = this.xScale.invert( m[0] ) - (m[0] / width) * newRange + (dx / width) * newRange;
  lower = Math.max( full[0], Math.min( lower, full[1] - newRange ) );

  if( (Math.abs(lower - domain[0]) > 1.0E-9*range) || (Math.abs(newRange - range) > 1.0E-9*range) )
    this.setXAxisRange( lower, lower + newRange );
};//RelEffPlot.prototype.handleWheel


/* The mouse position (chart-area coordinates) during a drag-zoom; for a y-zoom it is clamped to the plot
   area, so dragging past the top or bottom zooms to that edge. */
RelEffPlot.prototype.dragMousePos = function() {
  const mouse = d3.mouse(this.mouseCapture.node());
  if( this.zoomAxis === 'y' )
    mouse[1] = Math.max( 0, Math.min( mouse[1], this.yScale.range()[0] ) );
  return mouse;
};//RelEffPlot.prototype.dragMousePos


/* What releasing a drag-zoom from `start` to `current` would do: null if the drag is too short to act on,
   otherwise `{zoomIn, mult}`.  Drag right (x) or down (y) zooms in; drag left (x) zooms all the way out,
   while drag up (y) widens the y-range `mult` times (Infinity meaning all the way out), with the factor
   picked by drag distance like the spectrum chart (but in fixed pixels, as this chart is often short). */
RelEffPlot.prototype.dragZoomAction = function( start, current ) {
  if( this.zoomAxis === 'y' ){
    const dy = current[1] - start[1];
    if( Math.abs(dy) <= 10 )
      return null;
    if( dy > 0 )
      return { zoomIn: true, mult: 1 };
    return { zoomIn: false, mult: ((-dy < 20) ? 2 : ((-dy < 35) ? 4 : Infinity)) };
  }

  const dx = current[0] - start[0];
  return (Math.abs(dx) <= 5) ? null : { zoomIn: (dx > 0), mult: Infinity };
};//RelEffPlot.prototype.dragZoomAction


RelEffPlot.prototype.handleMouseMove = function() {
  if (!this.leftMouseDown) return;

  this.updateZoomBox( this.leftMouseDown, this.dragMousePos() );
};//RelEffPlot.prototype.handleMouseMove


RelEffPlot.prototype.handleMouseUp = function() {
  if (!this.leftMouseDown) return;

  const start = this.leftMouseDown;
  const current = this.dragMousePos();
  const action = this.dragZoomAction( start, current );

  // Remove zoom box and text
  this.chartArea.selectAll(".zoom-box, .zoom-text").remove();

  if( action && (this.zoomAxis === 'y') ){
    if( action.zoomIn ){
      const startY = this.yScale.invert( start[1] ), endY = this.yScale.invert( current[1] );
      this.setYAxisRange( Math.min(startY, endY), Math.max(startY, endY) );
    }else{
      this.zoomOutY( action.mult );
    }
  }else if( action ){
    if( action.zoomIn ){
      const startEnergy = this.xScale.invert( start[0] ), endEnergy = this.xScale.invert( current[0] );
      this.setXAxisRange( Math.min(startEnergy, endEnergy), Math.max(startEnergy, endEnergy) );
    }else{
      this.zoomOut();
    }
  }

  // Clean up
  this.leftMouseDown = null;
  this.zoomStartX = null;
  this.zoomAxis = null;

  // Remove document listeners
  d3.select(document)
    .on("mousemove.releffzoom", null)
    .on("mouseup.releffzoom", null);
};//RelEffPlot.prototype.handleMouseUp


/* Draws the band showing the range a drag-zoom will zoom to: vertical for an x-zoom, horizontal for a y-zoom. */
RelEffPlot.prototype.updateZoomBox = function( start, current ) {
  const chartAreaWidth = Math.max( 0, this.xScale.range()[1] );
  const chartAreaHeight = Math.max( 0, this.yScale.range()[0] ); // Get the height from y scale range
  const action = this.dragZoomAction( start, current );
  const isY = (this.zoomAxis === 'y');

  let zoomBox = this.chartArea.select(".zoom-box");
  let zoomText = this.chartArea.select(".zoom-text");

  if (zoomBox.empty()) {
    zoomBox = this.chartArea.append("rect")
      .attr("class", "zoom-box")
      .style("stroke-width", 1);

    zoomText = this.chartArea.append("text")
      .attr("class", "zoom-text chartLineText")
      .style("text-anchor", "middle")
      .style("dominant-baseline", "middle")
      .style("pointer-events", "none");
  }

  const x = isY ? 0 : Math.min( start[0], current[0] );
  const y = isY ? Math.min( start[1], current[1] ) : 0;
  const width = isY ? chartAreaWidth : Math.abs( current[0] - start[0] );
  const height = isY ? Math.abs( current[1] - start[1] ) : chartAreaHeight;

  // Apply appropriate CSS class and text
  const zoomingIn = isY ? (current[1] > start[1]) : (current[0] > start[0]);
  zoomBox.attr("class", "zoom-box " + (zoomingIn ? "leftbuttonzoombox" : "leftbuttonzoomoutboxy"));

  let label = zoomingIn ? "Zoom In" : "Zoom Out";
  if( isY && action && !action.zoomIn )
    label += isFinite(action.mult) ? (" x" + action.mult) : " full";

  zoomBox
    .attr("x", x)
    .attr("y", y)
    .attr("width", width)
    .attr("height", height);

  // Position text in center of zoom box; hide it if the box is too small (or wont do anything)
  zoomText
    .text( label )
    .attr("x", x + width / 2)
    .attr("y", y + height / 2)
    .style("visibility", (isY ? !!action : (width > 50)) ? "visible" : "hidden");
};//RelEffPlot.prototype.updateZoomBox


/* Sets the energy range, re-fitting the y-axis to what is visible within it. */
RelEffPlot.prototype.setXAxisRange = function(minEnergy, maxEnergy) {
  this.xScale.domain([minEnergy, maxEnergy]);
  this.yScale.domain( this.autoYDomain( this.xScale.domain() ) );

  this.chartArea.selectAll('.xAxis').call(this.xAxis);
  this.chartArea.selectAll('.yAxis').call(this.yAxis);

  this.relayoutOrRedraw();
};//RelEffPlot.prototype.setXAxisRange


RelEffPlot.prototype.setYAxisRange = function( minY, maxY ) {
  this.yScale.domain( [minY, maxY] );
  this.chartArea.selectAll('.yAxis').call(this.yAxis);

  this.relayoutOrRedraw();
};//RelEffPlot.prototype.setYAxisRange


RelEffPlot.prototype.zoomOut = function() {
  // Zoom out to show the full data extent
  const xDomain = this.fullXDomain();
  this.setXAxisRange( xDomain[0], xDomain[1] );
};//RelEffPlot.prototype.zoomOut


/* Widens the y-range `mult` times about its center, staying within the automatic y-range for the
   visible energies; a `mult` of Infinity zooms all the way out to that range. */
RelEffPlot.prototype.zoomOutY = function( mult ) {
  const full = this.autoYDomain( this.xScale.domain() );
  const domain = this.yScale.domain();
  const span = mult * (domain[1] - domain[0]);

  if( !isFinite(span) || (span >= (full[1] - full[0])) ){
    this.setYAxisRange( full[0], full[1] );
    return;
  }

  // Expand about the center, then shift back inside the full range if we went past either end
  let lower = 0.5*(domain[0] + domain[1]) - 0.5*span;
  lower = Math.max( full[0], Math.min( lower, full[1] - span ) );
  this.setYAxisRange( lower, lower + span );
};//RelEffPlot.prototype.zoomOutY


/* The energy range that shows every data point, with some padding. */
RelEffPlot.prototype.fullXDomain = function() {
  let min_x = Number.MAX_VALUE;
  let max_x = -Number.MAX_VALUE;

  (this.datasets || []).forEach(function(dataset) {
    if (dataset.data_vals && dataset.data_vals.length > 0) {
      min_x = Math.min(min_x, d3.min(dataset.data_vals, function(d) { return d.energy; }));
      max_x = Math.max(max_x, d3.max(dataset.data_vals, function(d) { return d.energy; }));
    }
  });

  // Handle edge case of no data
  if (min_x === Number.MAX_VALUE) min_x = 0;
  if (max_x === -Number.MAX_VALUE) max_x = 3000;

  const initial_x_range = max_x - min_x;
  let lower_padding = 0.025 * initial_x_range;
  if (Math.abs(initial_x_range) < 1.0) {
    lower_padding = 1.0;
  } else if (((min_x - lower_padding) < 50.0) || (lower_padding > 0.2 * min_x)) {
    lower_padding = 0.2 * min_x;
  }

  return [min_x - lower_padding, max_x + 0.025 * initial_x_range];
};//RelEffPlot.prototype.fullXDomain


/* The y-range that shows every data point, and the fit curves with their uncertainty bands, within the
   energy range `xDomain`, with some padding. */
RelEffPlot.prototype.autoYDomain = function( xDomain ) {
  let min_y = Number.MAX_VALUE;
  let max_y = -Number.MAX_VALUE;
  const addY = function( y ){
    if( isFinite(y) ){
      min_y = Math.min( min_y, y );
      max_y = Math.max( max_y, y );
    }
  };

  // Sample the curves at the same density they are drawn with
  const num_eqn_points = Math.max( 1, Math.floor( Math.abs(this.xScale.range()[1] - this.xScale.range()[0]) / 4 ) );

  (this.datasets || []).forEach(function(dataset) {
    (dataset.data_vals || []).forEach(function(d) {
      if( (d.energy >= xDomain[0]) && (d.energy <= xDomain[1]) )
        addY( d.eff );
    });

    if( dataset.fit_eqn ){
      for( let i = 0; i <= num_eqn_points; ++i ){
        const ene = xDomain[0] + i * (xDomain[1] - xDomain[0]) / num_eqn_points;
        const eff = dataset.fit_eqn(ene);
        const eff_uncert = dataset.fit_uncert_fcn ? dataset.fit_uncert_fcn(ene) : null;
        addY( eff - (eff_uncert ? eff_uncert : 0) );
        addY( eff + (eff_uncert ? eff_uncert : 0) );
      }
    }
  });

  // Handle edge case of no visible data
  if (min_y > max_y) {
    min_y = 0;
    max_y = 1;
  }

  const initial_y_range = (max_y > min_y) ? (max_y - min_y) : ((max_y !== 0) ? Math.abs(max_y) : 1);
  return [min_y - 0.1 * initial_y_range, max_y + 0.2 * initial_y_range];
};//RelEffPlot.prototype.autoYDomain


/* After an axis range changes, the y tick labels may be a different width than the left gutter was laid out
   for (e.g., "0.04" -> "0.0012", which would get clipped), so redo the layout, keeping the zoom, if so. */
RelEffPlot.prototype.relayoutOrRedraw = function() {
  if( !this.datasets )
    return;

  const yAxisNode = this.chartArea.select('.yAxis').node();
  const ytickw = yAxisNode ? yAxisNode.getBBox().width : 0;
  if( Math.abs( ytickw - this.layoutYTickWidth ) > 1 )
    this.layoutAndDraw( true );
  else
    this.redrawData();
};//RelEffPlot.prototype.relayoutOrRedraw


RelEffPlot.prototype.redrawData = function() {
  // This method redraws the paths and points with the current scales
  // We need to regenerate fit equation points for the current x-axis range
  if (!this.datasets || this.datasets.length === 0) return;
  
  const self = this;
  const xDomain = this.xScale.domain();
  const min_x = xDomain[0];
  const max_x = xDomain[1];
  const chartAreaWidth = this.xScale.range()[1] - this.xScale.range()[0];
  const num_eqn_points = chartAreaWidth / 4;
  
  // Regenerate fit equation points for current x-axis range
  this.datasets.forEach(function(dataset, index) {
    if (dataset.fit_eqn) {
      dataset.fit_eqn_points = [];
      for (let i = 0; i < num_eqn_points; ++i) {
        const ene = min_x + i * (max_x - min_x) / num_eqn_points;
        const eff = dataset.fit_eqn(ene);
        const eff_uncert = dataset.fit_uncert_fcn ? dataset.fit_uncert_fcn(ene) : null;
        
        dataset.fit_eqn_points.push({ 
          energy: ene, 
          eff: eff, 
          eff_uncert: eff_uncert 
        });
      }
    }
  });
  
  // Redraw paths with new data
  this.chartArea.selectAll("path.line").remove();
  this.chartArea.selectAll("path.RelEffPlotErrorBounds").remove();
  
  // Redraw fit lines and uncertainty bounds.
  //  Append into plotGroup (clipped, and BELOW the data-point circles and mouse-capture rect in
  //  hit-test order), and with pointer-events disabled on the line - otherwise the 2px curve
  //  stroke sits on top of the zoom-capture rect and the circles, silently eating the mousedown
  //  that starts a drag-zoom and the mouseover that shows a point's tooltip (this regressed in a
  //  merge once already; this is now the only place the curves are drawn).
  this.datasets.forEach(function(dataset, index) {
    if (dataset.fit_eqn && dataset.fit_eqn_points && dataset.fit_eqn_points.length > 0) {
      // Create path for fit line
      const valueline = d3.svg.line()
        .x(function(d) { return self.xScale(d.energy); })
        .y(function(d) { return self.yScale(d.eff); });

      const path = self.plotGroup.append("path")
        .attr("class", "line dataset-" + index)
        .attr("d", valueline(dataset.fit_eqn_points))
        .style("stroke", self.getDatasetColor(index))
        .style("pointer-events", "none");

      // Create path for uncertainty bounds if available
      if (dataset.fit_uncert_fcn) {
        const area = d3.svg.area()
          .x(function(d) { return self.xScale(d.energy); })
          .y0(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff - d.eff_uncert) : null; })  /* band is +-1 sigma, matching the data error bars */
          .y1(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff + d.eff_uncert) : null; });

        const path_uncert = self.plotGroup.append("path")
          .attr("class", "RelEffPlotErrorBounds dataset-" + index)
          .attr("d", area(dataset.fit_eqn_points))
          .style("fill", self.getDatasetColor(index, 0.1))
          .style("stroke", "none")
          .style("pointer-events", "none");
      }
    }
  });
  
  // Update data points positions
  this.chartArea.selectAll("circle")
    .attr("cx", function(d) { return self.xScale(d.energy); })
    .attr("cy", function(d) { 
      const val = self.yScale(d.eff);
      return isNaN(val) ? 0 : val;
    });
    
  // Update error bars positions
  this.chartArea.selectAll("line.errorbar")
    .attr('x1', function(d) { return self.xScale(d.energy); })
    .attr('x2', function(d) { return self.xScale(d.energy); })
    .attr('y1', function(d) {
      const val = self.yScale(d.eff + d.eff_uncert);
      return isNaN(val) ? 0 : val;
    })
    .attr('y2', function(d) {
      const val = self.yScale(d.eff - d.eff_uncert);
      return isNaN(val) ? 0 : val;
    });
};//RelEffPlot.prototype.redrawData


RelEffPlot.prototype.setMargins = function( margins ){
  if( !margins )
    return;
    
  if( typeof margins.top === "number" )
    this.options.margins.top = margins.top;
  if( typeof margins.right === "number" )
    this.options.margins.right = margins.right;
  if( typeof margins.bottom === "number" )
    this.options.margins.bottom = margins.bottom;
  if( typeof margins.left === "number" )
    this.options.margins.left = margins.left;
    
  this.handleResize();
}//RelEffPlot.prototype.setMargins


RelEffPlot.prototype.setXAxisTitle = function( title, dontCallResize ){
  if( !title || (title.length === 0) )
  {
    this.options.xAxisTitle = null;
    if( this.xaxistitle ){
      this.xaxistitle.remove();
      this.xaxistitle = null;
    }
    return;
  }//if( we are getting rid of the title )
  
  this.options.xAxisTitle = title;
  if( !this.xaxistitle ){
    this.xaxistitle = this.svg.append("text")
      .attr("class", "xaxistitle")
      .attr("x", 0.5*this.chart.clientWidth ) //roughly position title - will be updated when data set or resized
      .attr("y", this.chart.clientHeight - this.options.margins.bottom )
      .style("text-anchor", "middle")
      .style("dominant-baseline", "text-after-edge");
  }
  
  this.xaxistitle.text(title);
  
  if( !dontCallResize )
    this.handleResize();
}//RelEffPlot.prototype.setXAxisTitle


RelEffPlot.prototype.setYAxisTitle = function( title, dontCallResize ){
  if( !title || (title.length === 0) )
  {
    this.options.yAxisTitle = null;
    if( this.yaxistitle ){
      this.yaxistitle.remove();
      this.yaxistitle = null;
    }
    return;
  }//if( we are getting rid of the title )
  
  this.options.yAxisTitle = title;
  if( !this.yaxistitle ){
    this.yaxistitle = this.svg.append("text")
      .attr("class", "yaxistitle")
      .attr("transform", "rotate(-90)")
      .attr("y", -this.options.margins.left )
      .attr("x", -0.5*this.chart.clientHeight )
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .style("dominant-baseline", "text-after-edge");
  }
  
  this.yaxistitle.text(title);
  
  if( !dontCallResize )
    this.handleResize();
}//RelEffPlot.prototype.setYAxisTitle

/** Sets the data that will be plotted.
 
 @param datasets An array of datasets, each of which is an object with the following properties:
   - data_vals: An array of data points, each of which is an object with properties "energy", "mean", "counts", "counts_uncert", "eff", "eff_uncert", and "nuc_info".
   - fit_eqn: A function that takes a number as an argument (the x-value in the plot) and returns a y-value.
   - chi2_txt: A string that gets displayed on the chart.
   - fit_uncert_fcn: A function that gives the uncertainty on the solid line.
 
 */
RelEffPlot.prototype.setRelEffData = function (datasets) {
  // If not an array of dataset objects, convert to array format for backward compatibility
  if (!Array.isArray(datasets) || datasets.length === 0 || 
      (datasets.length > 0 && (!datasets[0].hasOwnProperty('data_vals') && !datasets[0].hasOwnProperty('fit_eqn')))) {
    // Handle backward compatibility - convert single dataset to array of datasets
    if (arguments.length >= 2) {
      const data_vals = arguments[0];
      const fit_eqn = arguments[1];
      const chi2_txt = arguments[2];
      const fit_uncert_fcn = arguments[3];
      
      datasets = [{ data_vals, fit_eqn, chi2_txt, fit_uncert_fcn }];
    } else {
      datasets = [];
    }
  }

  // Store the datasets
  this.datasets = datasets;

  // Extract chi2_txt values from all datasets into an array
  this.chi2TxtArray = datasets.map(dataset => 
    (typeof dataset.chi2_txt === "string" && dataset.chi2_txt.length > 0) ? dataset.chi2_txt : null
  ).filter(txt => txt !== null);
  
  // For backward compatibility
  this.chi2Txt = this.chi2TxtArray.length > 0 ? this.chi2TxtArray[0] : null;

  this.layoutAndDraw( false );
};//RelEffPlot.prototype.setRelEffData


/** Lays out the axes, titles, and plot area for the chart's current size and axis ranges, and (re)creates
 everything drawn in it.  If `keepZoom` is true the current axis ranges are kept, otherwise the chart is
 zoomed all the way out.
 */
RelEffPlot.prototype.layoutAndDraw = function( keepZoom ) {
  const self = this;
  const datasets = this.datasets;

  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;

  const titlePad = this.options.titleToAxisPadding;
  const xaxistitleBB = this.xaxistitle ? this.xaxistitle.node().getBBox() : null;
  const xtitleh = xaxistitleBB ? xaxistitleBB.height + titlePad : 0;

  const yaxistitleBB = this.yaxistitle ? this.yaxistitle.node().getBBox() : null;
  const ytitleh = yaxistitleBB ? yaxistitleBB.height + titlePad : 0;

  // Remove any existing chart info titles
  if (this.chartInfoTitles) {
    this.chartInfoTitles.forEach(title => {
      if (title) title.remove();
    });
  }
  
  // Initialize array for chart titles
  this.chartInfoTitles = [];
  let chi2_txt_height = 0;
  
  // Create and position text elements if we have any chi2_txt values
  if (this.chi2TxtArray.length > 0) {
    // First create temporary text elements to measure their widths
    const tempTexts = this.chi2TxtArray.map(txt => {
      return this.svg.append("text")
        .attr("class", "ChartInfoTitle temp")
        .text(txt)
        .style("visibility", "hidden"); // Hide it while measuring
    });
    
    // Measure all text elements to find the longest one
    const textWidths = tempTexts.map(txt => txt.node().getBBox().width);
    const maxTextWidth = Math.max(...textWidths) + 20; // Add 20px padding
    
    // Remove temporary elements
    tempTexts.forEach(txt => txt.remove());
    
    // Calculate starting position for the x-coordinate of all texts.
    //  The equation text is right-aligned, but a long one (e.g. a physical model) can be wider than the
    //  chart and run off the left edge - keep a gutter for the omitted-points icon, which lives in the
    //  top-left corner, so the two never overlap (see `updateOmittedPointsInfo`).
    const omittedIconGutter = datasets.some( d => (d.omitted_pts && d.omitted_pts.length) ) ? 20 : 0;
    const startX = Math.max( omittedIconGutter, parentWidth - maxTextWidth );
    
    // Create the actual text elements stacked vertically
    let currentY = 5;
    this.chi2TxtArray.forEach((txt, index) => {
      const textElem = this.svg.append("text")
        .attr("class", "ChartInfoTitle dataset-" + index)
        .attr("y", currentY)
        .attr("x", startX)
        .style("text-anchor", "start")  // Align text to start for consistent left edge
        .attr("dominant-baseline", "hanging") // Align text to the top of the text box
        .text(txt)
        .style("fill", this.getDatasetColor(index));
      
      this.chartInfoTitles.push(textElem);
      
      // Update Y position for next text element
      const textHeight = textElem.node().getBBox().height;
      currentY += textHeight + 2; // Add small spacing between lines
      
      // Track total height for layout calculations
      chi2_txt_height = currentY;
    });
  }
  
  // Use half the height for padding if there are any titles
  const chi2_txt_pad = chi2_txt_height > 0 ? chi2_txt_height + 5 : 0;

  // We will perform an initial estimate of the chart are - using the axis we have.
  //  Later on we will put current values into the axis, and then get the size, and
  //  update things
  let xticks = this.chartArea.selectAll('.xAxis');
  let xtickh = xticks.empty() ? 24.5 : d3.max(xticks[0], function (t) { return d3.select(t).node().getBBox().height; });

  let yticks = this.chartArea.selectAll('.yAxis');
  let ytickw = yticks.empty() ? 37 : Math.max(37, d3.max(yticks[0], function (t) { return d3.select(t).node().getBBox().width; })); //If we have no labels, we'll get like 7 px, so will require at least 37 px, which is the nominal value we should get

  let chartAreaWidth = Math.max( 0, parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw );
  let chartAreaHeight = Math.max( 0, parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad );

  // Pick the axis ranges before measuring the tick labels, so the gutter fits the labels actually shown
  //  (autoYDomain samples the curves at a density set by the plot width, so set the ranges first).
  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);

  if( !keepZoom ){
    this.xScale.domain( this.fullXDomain() );
    this.yScale.domain( this.autoYDomain( this.xScale.domain() ) );
  }

  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);

  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);

  // Now that we have x and y axis tick values (and their text), we will recompute axis sizes
  //  and correct chart area and positions using this information.  I dont think axis label
  //  texts should change an amount that effects things.
  xticks = this.chartArea.selectAll('.xAxis');
  xtickh = xticks.empty() ? 7 : d3.max(xticks[0], function(t){ return d3.select(t).node().getBBox().height; });
          
  yticks = this.chartArea.selectAll('.yAxis');
  ytickw = yticks.empty() ? 7 : d3.max(yticks[0], function(t){ return d3.select(t).node().getBBox().width; });
  this.layoutYTickWidth = ytickw;  // see relayoutOrRedraw

  // Clamp to non-negative - the container may not yet be laid out, and SVG attrs reject negative values.
  chartAreaWidth = Math.max( 0, parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw );
  chartAreaHeight = Math.max( 0, parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad );

  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);
       
  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);
        
  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);
               
  if (this.xaxistitle)
    this.xaxistitle
      .attr("x", this.options.margins.left + ytitleh + ytickw + 0.5*chartAreaWidth)
      .attr("y", parentHeight - this.options.margins.bottom);
         
  if (this.yaxistitle)
    this.yaxistitle
      .attr("y", this.options.margins.left)
      .attr("x", this.options.margins.top + chi2_txt_pad - 0.5*chartAreaHeight);
         
  this.chartArea.selectAll('.xAxis')
    .attr("transform", "translate(0," + chartAreaHeight + ")");
  this.chartArea
    .attr("transform", "translate(" + (this.options.margins.left + ytickw + ytitleh)
                       + "," + (this.options.margins.top + chi2_txt_pad) + ")");
    
  // Update clipping rectangle dimensions
  this.chartArea.select("#releffplot-clip-" + this.chart.id + " rect")
    .attr("width", chartAreaWidth)
    .attr("height", chartAreaHeight);
    
  // For some reason the circles position arent being updated on resize, so we'll just remove them first
  //  With later versions of D3 there is a .merge() function that is maybe relevant
  this.plotGroup.selectAll("circle").remove();
  this.plotGroup.selectAll("line.errorbar").remove();

  // Also drop any zoom box orphaned by a never-delivered mouseup - it paints above the data
  //  points and would otherwise block hover/drag in its slab until the next click.
  this.chartArea.selectAll(".zoom-box, .zoom-text").remove();
    
  // Calculate min and max counts across all datasets for radius scaling
  let minCounts = Number.MAX_VALUE;
  let maxCounts = -Number.MAX_VALUE;
  
  datasets.forEach(function(dataset) {
    if (dataset.data_vals && dataset.data_vals.length > 0) {
      dataset.data_vals.forEach(function(d) {
        if (d.counts > maxCounts) maxCounts = d.counts;
        if (d.counts < minCounts) minCounts = d.counts;
      });
    }
  });
  
  // Function to calculate radius based on counts (sqrt scaling, radius 1-4)
  function calculateRadius(counts) {
    if (minCounts === maxCounts) return 3; // Default radius if all counts are equal
    
    // Normalize counts to 0-1 range, then apply sqrt, then scale to 1-4
    const normalized = (counts - minCounts) / (maxCounts - minCounts);
    const sqrtNormalized = Math.sqrt(normalized);
    return 2 + sqrtNormalized * 3; // Scale from 2 to 5
  }

  // The color to draw a point with: the dominant nuclide's color, when a single nuclide accounts for more than
  //  half of the counts.  Returns null otherwise, so the ".multiiso"/".noiso" CSS rules take over.
  function pointNucColor(d) {
    if( !d.nuc_info || (d.nuc_info.length === 0) )
      return null;

    let max_contrib = 0, sum_contrib = 0, dominant_color = null;
    for( const el of d.nuc_info ) {
      const contrib = el.rel_act * el.br;
      sum_contrib += contrib;
      if( contrib > max_contrib ) {
        max_contrib = contrib;
        dominant_color = el.color;
      }
    }

    if( !dominant_color || (typeof dominant_color !== 'string') || (dominant_color.length === 0) )
      return null;

    return (max_contrib > 0.5*sum_contrib) ? dominant_color : null;
  }//pointNucColor

  // A point is in a "shared region" when another relative efficiency curve also has gammas at this energy.
  //  How the counts divide between the curves is not resolved by this peak itself (co-located gammas are
  //  degenerate locally); the division is fitted through the attenuation physics / curve shapes, constrained
  //  by the rest of the spectrum.  Such a point is faded (see `pointOpacity`) in proportion to this curve's
  //  share, so its opacity roughly shows how much this curve contributes at that energy.
  const blendedThreshold = 0.9;
  function isBlendedPoint(d) {
    return ((typeof d.blend_frac === 'number') && (d.blend_frac < blendedThreshold));
  }

  // For each energy shared by more than one rel. eff. curve, find the largest curve share there.  The curve
  //  contributing most is drawn fully opaque; the others are faded in proportion to their share of that maximum,
  //  so the point opacities roughly convey how much each curve contributes at that energy.  A shared cluster has
  //  the same `energy` on every curve's point (see ObsEff::energy), so rounding to the emitted precision keys them
  //  together.
  const maxBlendByEnergy = {};
  datasets.forEach(function(dataset) {
    const dv = dataset.data_vals;
    if( !dv )
      return;
    for( const d of dv ){
      const key = d.energy.toFixed(2);
      const frac = (typeof d.blend_frac === 'number') ? d.blend_frac : 1;
      if( !(key in maxBlendByEnergy) || (frac > maxBlendByEnergy[key]) )
        maxBlendByEnergy[key] = frac;
    }
  });

  // Opacity a point is drawn at: 1 for a sole-owner point, and for a shared point its share relative to the
  //  largest-contributing curve at that energy (floored so a surviving minor point stays visible/clickable).
  const opacityFloor = 0.25;
  function pointOpacity(d) {
    if( !isBlendedPoint(d) )
      return 1;
    const maxf = maxBlendByEnergy[d.energy.toFixed(2)] || 1;
    const rel = (maxf > 0) ? (d.blend_frac / maxf) : 1;
    return Math.max( opacityFloor, Math.min( 1, rel ) );
  }

  // Plot all data points across all datasets
  datasets.forEach(function(dataset, datasetIndex) {
    const data_vals = dataset.data_vals;
    
    // Everything below here requires data points, so if we dont have any, we can skip this
    if (!data_vals || data_vals.length === 0)
      return;
      
    // Add error bars
    let lines = self.plotGroup.selectAll("line.errorbar.dataset-" + datasetIndex)
      .data(data_vals);
    
    lines.enter()
      .append('line')
      .attr("class", "errorbar dataset-" + datasetIndex)
      .attr('x1', function(d) { return self.xScale(d.energy); })
      .attr('x2', function(d) { return self.xScale(d.energy); })
      .attr('y1', function(d) {
        const val = self.yScale(d.eff + d.eff_uncert);
        return isNaN(val) ? 0 : val;
      })
      .attr('y2', function(d) {
        const val = self.yScale(d.eff - d.eff_uncert);
        return isNaN(val) ? 0 : val;
      })
      .style("stroke", function(d) {
        // Fall back to dataset color if no valid color from nuc_info
        return pointNucColor(d) || self.getDatasetColor(datasetIndex);
      })
      .style("opacity", pointOpacity);

    // Add the data points
    self.plotGroup
      .selectAll("circle.dataset-" + datasetIndex)
      .data(data_vals)
      .enter()
      .append("circle")
      .attr("r", function(d) { return calculateRadius(d.counts); })
      .attr("cx", function (d) {
        return self.xScale(d.energy)
      })
      .attr("cy", function (d) {
        const val = self.yScale(d.eff);
        return isNaN(val) ? 0 : val;
      })
      .attr("class", function (d) {
        let baseClass = "dataset-" + datasetIndex + " ";
        if( isBlendedPoint(d) )
          baseClass += "blended ";

        // Guard nuc_info: if this callback throws for a single datum, d3 never attaches the
        //  mouseover/mouseout handlers to ANY circle (operators apply selection-wide in order),
        //  killing every tooltip on the chart.
        if (!d.nuc_info || (d.nuc_info.length === 0))
          return baseClass + "noiso";

        //Return the dominant nuclide, if no nuclide is over 50%, we'll use color defined by CSS "multiiso"
        //  TODO: use a gradient to indicate relative components
        let max_contrib = 0, sum_contrib = 0, max_nuc = null, dominant_color = null;
        for (const el of d.nuc_info) {
          const contrib = el.rel_act * el.br;
          sum_contrib += contrib;
          if (contrib > max_contrib) {
            max_contrib = contrib;
            max_nuc = el.nuc;
            dominant_color = el.color; // Get color from dominant nuclide
          }
        }

        if( !dominant_color )
          return baseClass + "noiso";

        if( max_contrib < 0.5*sum_contrib )
          return baseClass + "multiiso";
        
        if (max_nuc && (max_contrib > 0.5*sum_contrib)) {
          //remove problematic characters from max_nuc to create a valid CSS class name
          max_nuc = max_nuc.replace(/[^a-zA-Z0-9_-]/g, '');
          if (max_nuc.length === 0)
            return baseClass + "multiiso";
          
          const first_char = max_nuc.charAt(0);
          if ((first_char >= '0' && first_char <= '9') || (first_char === '-'))
            max_nuc = 'nuc-' + max_nuc;
          
          return baseClass + max_nuc;
        }
        
        return baseClass + "multiiso";
      })
      .style("fill", function(d) {
        // A null return leaves the fill to the ".RelEffPlot circle.noiso"/".multiiso" CSS rules
        if (d.nuc_info && d.nuc_info.length > 0)
          return pointNucColor(d);

        // Fall back to dataset color if no nuc_info at all
        return self.getDatasetColor(datasetIndex);
      })
      .style("stroke", null)
      // Shared points are faded in proportion to this curve's contribution (see `pointOpacity`), rather than
      //  drawn hollow, so the opacity roughly conveys how much this curve contributes at that energy.
      .style("opacity", pointOpacity)
      .on("mouseover", function (d, i) {
        self.tooltip.transition()
          .duration(200)
          .style("opacity", .9);

        // Bring a faded (shared) point up on hover so it can be inspected; use style opacity to match how the
        //  point is drawn (a `style` opacity would otherwise override an `attr` one).
        d3.select(this).transition()
          .duration('50')
          .style('opacity', Math.max( 0.85, pointOpacity(d) ))
          .attr("r", 6);

        const eqn_eff = dataset.fit_eqn ? dataset.fit_eqn(d.energy) : null;

        let txt = "";

        // With multiple rel. eff. curves, name the curve this point belongs to up-front, so the later
        //  "this curve" attribution wording (and marker coloring) is never ambiguous.
        if( datasets.length > 1 ){
          const curve_name = dataset.curve_name ? dataset.curve_name : ("Curve " + (datasetIndex+1));
          txt += "<div>Curve: " + curve_name + "</div>";
        }

        txt += "<div>Energy: " + (d.mean ? d.mean.toFixed(2) : d.energy.toFixed(2)) + " keV</div>"
          + "<div>Peak Area: " + d.counts.toFixed(1) + " &plusmn; " + d.counts_uncert.toFixed(1) + "</div>";


        // Add fit equation value if available for this dataset
        if (eqn_eff) {
          txt += "<div>RelEff Curve: " + eqn_eff.toPrecision(4);
          if( dataset.fit_uncert_fcn ){
            const crv_uncert = dataset.fit_uncert_fcn(d.energy);
            if( (crv_uncert !== null) && (crv_uncert > 0) )
              txt += " &plusmn; " + crv_uncert.toPrecision(4) + " (1&sigma;)";
          }
          txt += "</div>";
        }

        /* Data-point error bars/uncertainties are 1-sigma STATISTICAL only (counts_uncert/src_counts)
           - deliberately excluding any user "Add. Uncert.", which is a model-error term the fitted
           band accounts for, not a measurement error of the point. */
        txt += "<div>Measured RelEff: " + d.eff.toPrecision(4);
        if( d.eff_uncert && (d.eff_uncert > 0) ){
          txt += " &plusmn; " + d.eff_uncert.toPrecision(4) + " (1&sigma;)";
        }
        txt += "</div>";

        if( d.eff_uncert && (d.eff_uncert > 0) && eqn_eff ){
          const nsigma_off = Math.abs(d.eff - eqn_eff) / d.eff_uncert;
          txt += "<div>Peak fit " + nsigma_off.toPrecision(3) + " &sigma; from curve (peak uncert. only)</div>";
        }

        for (const el of d.nuc_info) {
          txt += "<div>&nbsp;&nbsp;" + el.nuc + ": br=" + el.br.toPrecision(4);
          if (el.rel_act)
            txt += ", RelAct=" + el.rel_act.toPrecision(4);
          txt += "</div>";
        }

        // When the problem has more than one rel. eff. curve, always show what fraction of this energy's summed
        //  counts the fit attributes to this curve - even when that is ~100% - so the split is never ambiguous.
        if( (datasets.length > 1) && (typeof d.blend_frac === 'number') ){
          txt += "<div class=\"RelEffPlotBlendNote\">Shared region: the fit attributes "
               + (100*d.blend_frac).toFixed(0)
               + "% of the summed counts at this energy to this curve.</div>";
        }

        self.tooltip.html(txt);

        // Make it so tooltip doesnt extend above/below/left/right of chart area.
        //  The tooltip is position:fixed, i.e. viewport coordinates - so use clientX/clientY
        //  (pageX/pageY differ by the document scroll), and clamp to the viewport: a tall tooltip
        //  (e.g. multi-curve shared-region notes with several nuclide lines) near the window top
        //  used to get a negative `top` and render entirely off-screen.
        const svg_location = d3.mouse(self.chartArea.node());
        const svg_bb = self.chartArea.node().getBoundingClientRect();
        const tt_bb = self.tooltip.node().getBoundingClientRect();
        const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
        const render_top = (svg_location[1] + tt_bb.height + 4 + 15 > svg_bb.height);
        let x_offset = render_right ? d3.event.clientX + 10 : d3.event.clientX - 10 - tt_bb.width;
        let y_offset = render_top ? d3.event.clientY - 10 - tt_bb.height : d3.event.clientY + 4;
        x_offset = Math.max( 2, Math.min( x_offset, window.innerWidth - tt_bb.width - 2 ) );
        y_offset = Math.max( 2, Math.min( y_offset, window.innerHeight - tt_bb.height - 2 ) );

        self.tooltip
          .style("left", x_offset + "px")
          .style("top", y_offset + "px");
      })
      .on('mouseout', function (d, i) {
        d3.select(this).transition()
          .duration('50')
          .style('opacity', pointOpacity(d))
          .attr("r", calculateRadius(d.counts));
        self.tooltip.transition()
          .duration(500)
          .style("opacity", 0);
      });
  });
  
  // Position the mouse capture rectangle to cover the chart area
  if (this.mouseCapture) {
    this.mouseCapture
      .attr("width", chartAreaWidth)
      .attr("height", chartAreaHeight);
  }

  this.updateOmittedPointsInfo();

  // Draws the fit curves (and uncertainty bands) for the current x-range
  this.redrawData();
};//RelEffPlot.prototype.layoutAndDraw


/* Shows a small info icon in the lower-left of the chart whenever some data points were intentionally left
   off it (see `RelActAutoSolution::show_obs_eff_point`), which on hover/tap lists the omitted energies and
   why each was left off.  Without this the points just silently disappear, and the user has no way to check
   that what was dropped is believable.
 */
RelEffPlot.prototype.updateOmittedPointsInfo = function () {
  const self = this;

  if( this.omittedInfoIcon ){
    this.omittedInfoIcon.remove();
    this.omittedInfoIcon = null;
  }
  if( this.omittedInfoPanel ){
    this.omittedInfoPanel.remove();
    this.omittedInfoPanel = null;
  }

  const datasets = this.datasets ? this.datasets : [];
  let num_omitted = 0;
  datasets.forEach( function(dataset){
    if( dataset.omitted_pts && dataset.omitted_pts.length )
      num_omitted += dataset.omitted_pts.length;
  } );

  if( num_omitted === 0 )
    return;

  // The icon/panel are absolutely positioned within the chart's container; make sure that container is
  //  actually a positioning ancestor, otherwise they escape to whatever dialog encloses the chart.
  if( window.getComputedStyle(this.chart).position === "static" )
    this.chart.style.position = "relative";

  const fmtCounts = function(v){
    return ((typeof v === 'number') && isFinite(v)) ? v.toPrecision(4) : "&mdash;";
  };

  let html = "<div class=\"RelEffPlotOmittedTitle\">Points not shown (" + num_omitted + ")</div>"
           + "<div class=\"RelEffPlotOmittedExplain\">These peaks were left off the chart because their"
           + " measured efficiency would not be meaningful:</div><table>"
           + "<tr class=\"RelEffPlotOmittedHeader\"><td>Energy</td><td>Source(s)</td>"
           + "<td>Fit counts</td><td>Expected</td><td>Why not shown</td></tr>";
  datasets.forEach( function(dataset, datasetIndex){
    if( !dataset.omitted_pts || !dataset.omitted_pts.length )
      return;

    if( datasets.length > 1 ){
      html += "<tr><td colspan=\"5\" class=\"RelEffPlotOmittedCurve\" style=\"color:"
            + self.getDatasetColor(datasetIndex) + "\">"
            + (dataset.curve_name ? dataset.curve_name : ("Curve " + (datasetIndex+1))) + "</td></tr>";
    }

    dataset.omitted_pts.forEach( function(pt){
      html += "<tr><td>" + pt.energy.toFixed(2) + " keV</td>"
            + "<td>" + (pt.srcs ? pt.srcs : "") + "</td>"
            + "<td>" + fmtCounts(pt.fit_counts) + "</td>"
            + "<td>" + fmtCounts(pt.expected_counts) + "</td>"
            + "<td>" + (pt.reason ? pt.reason : "") + "</td></tr>";
    } );
  } );
  html += "</table>";

  this.omittedInfoPanel = d3.select(this.chart).append("div")
    .attr("class", "RelEffPlotOmittedPanel")
    .style("opacity", 0)
    .style("visibility", "hidden")
    .style("pointer-events", "none")
    .html(html);

  // The panel is scrollable, so the user has to be able to move the mouse off the icon and onto it without
  //  it vanishing: hiding is delayed, and cancelled when the pointer enters either the icon or the panel.
  //  Clicking the icon pins the panel open (so a long list can be read without holding the mouse still).
  let hideTimer = null, pinned = false;

  const cancelHide = function(){
    if( hideTimer ){
      clearTimeout( hideTimer );
      hideTimer = null;
    }
  };

  const showPanel = function(){
    cancelHide();
    self.omittedInfoPanel.style("visibility", "visible").style("pointer-events", "auto")
      .transition().duration(150).style("opacity", 1);
  };

  const hidePanel = function(){
    cancelHide();
    pinned = false;
    // Stop blocking the chart immediately: if the fade-out transition is interrupted, the "end"
    //  callback never runs and the panel would stay visibility:visible at opacity 0 - an invisible
    //  near-full-chart overlay eating every hover and drag underneath.
    self.omittedInfoPanel.style("pointer-events", "none")
      .transition().duration(200).style("opacity", 0)
      .each("end", function(){
        if( !pinned )
          self.omittedInfoPanel.style("visibility", "hidden");
      });
  };

  const scheduleHide = function(){
    if( pinned )
      return;
    cancelHide();
    hideTimer = setTimeout( function(){ hideTimer = null; hidePanel(); }, 400 );
  };

  this.omittedInfoPanel
    .on("mouseover", cancelHide)
    .on("mouseout", scheduleHide);

  this.omittedInfoIcon = d3.select(this.chart).append("div")
    .attr("class", "RelEffPlotOmittedIcon")
    .attr("title", num_omitted + " data point(s) not shown - click to keep open")
    .text("i")
    .on("mouseover", showPanel)
    .on("mouseout", scheduleHide)
    .on("click", function(){
      if( pinned ){
        hidePanel();
      }else{
        pinned = true;
        showPanel();
      }
      d3.event.stopPropagation();
    });
};//RelEffPlot.prototype.updateOmittedPointsInfo

// Helper function to get colors for different datasets
RelEffPlot.prototype.getDatasetColor = function (index, alpha) {
  // If custom colors are provided for the first three datasets, use them
  if (index < this.customColors.length && this.customColors[index]) {
    const color = this.customColors[index];
    
    // If alpha is provided, return color with transparency
    if (alpha !== undefined) {
      // Handle both rgb and hex formats
      if (color.startsWith("rgb")) {
        const rgb = color.match(/\d+/g);
        return `rgba(${rgb[0]}, ${rgb[1]}, ${rgb[2]}, ${alpha})`;
      } else {
        // Convert hex to rgb
        const r = parseInt(color.slice(1, 3), 16);
        const g = parseInt(color.slice(3, 5), 16);
        const b = parseInt(color.slice(5, 7), 16);
        return `rgba(${r}, ${g}, ${b}, ${alpha})`;
      }
    }
    
    return color;
  }
  
  // Default color palette for additional datasets
  const colors = [
    "rgb(255, 127, 14)",  // orange
    "rgb(44, 160, 44)",   // green
    "rgb(214, 39, 40)",   // red
    "rgb(148, 103, 189)", // purple
    "rgb(140, 86, 75)",   // brown
    "rgb(227, 119, 194)", // pink
    "rgb(127, 127, 127)", // gray
    "rgb(188, 189, 34)",  // olive
    "rgb(23, 190, 207)"   // teal
  ];
  
  // Use modulo to cycle through colors if we have more datasets than colors
  // Adjust index to start from 0 for our fallback colors
  const adjustedIndex = index - this.customColors.length;
  const color = colors[adjustedIndex % colors.length];
  
  // If alpha is provided, return color with transparency
  if (alpha !== undefined) {
    // Extract RGB values and add alpha
    const rgb = color.match(/\d+/g);
    return `rgba(${rgb[0]}, ${rgb[1]}, ${rgb[2]}, ${alpha})`;
  }
  
  return color;
};

// Function to set custom colors from C++
RelEffPlot.prototype.setRelEffCurveColors = function (colors) {
  if (Array.isArray(colors) && colors.length > 0) {
    // Replace the first N colors where N is the minimum of colors.length and this.customColors.length
    for (let i = 0; i < Math.min(colors.length, this.customColors.length); i++) {
      if (colors[i] && typeof colors[i] === 'string') {
        this.customColors[i] = colors[i];
      }
    }
  }
  
  // If we have active datasets, redraw the chart to apply new colors
  if (this.datasets) {
    this.setRelEffData(this.datasets);
  }
};

RelEffPlot.prototype.handleResize = function () {
  if (this.datasets) {
    this.setRelEffData(this.datasets);
  } else if (this.data_vals) {
    // Extract chi2_txt from the array for backward compatibility
    const chi2_txt = this.chi2TxtArray && this.chi2TxtArray.length > 0 ? this.chi2TxtArray[0] : this.chi2Txt;
    
    // Backward compatibility for older code
    this.setRelEffData(this.data_vals, this.fit_eqn, chi2_txt, this.fit_uncert_fcn);
  }
  
  // Note: the mouse-capture rect is sized at the end of layoutAndDraw (which both branches above
  //  call) using the real plot-area size - axis tick/title space and the chi2-line padding
  //  subtracted.  A margins-only re-size here used to override that with a rect ~40px wider/taller
  //  than (and offset from) the plot area, so the capture geometry differed depending on whether
  //  the last call was setData or a resize.
};//RelEffPlot.prototype.handleResize
