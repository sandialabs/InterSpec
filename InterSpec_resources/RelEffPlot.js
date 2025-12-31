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
      
  // Set the dimensions of the canvas / graph
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  const chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right;
  const chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom;

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

  this.yAxis = d3.svg
    .axis()
    .scale(this.yScale)
    .tickFormat( function(x){
      if( x == 0 )
        return "0";
      if( x <= 1.0 )
        return x.toFixed(2);
      if( x <= 10.0 )
        return x.toFixed(1);
      return x.toFixed(0);
    } )
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
    .style("pointer-events", "none");
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
  this.zooming = false;
  
  // Setup mouse interactions for zoom functionality
  this.setupMouseInteractions();
}//RelEffPlot constructor


RelEffPlot.prototype.setupMouseInteractions = function() {
  const self = this;
  
  // Add mouse interaction rect to the chart area - will be positioned later in setRelEffData
  this.mouseCapture = this.chartArea.append("rect")
    .attr("class", "mouse-capture")
    .style("fill", "none")
    .style("pointer-events", "all")
    .on("mousedown", function() {
      const mouse = d3.mouse(this);
      self.leftMouseDown = mouse;
      self.zoomStartX = self.xScale.invert(mouse[0]);
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
};//RelEffPlot.prototype.setupMouseInteractions


RelEffPlot.prototype.handleMouseMove = function() {
  if (!this.leftMouseDown) return;
  
  const currentMouse = d3.mouse(this.mouseCapture.node());
  const startX = this.leftMouseDown[0];
  const currentX = currentMouse[0];
  
  // Determine if zooming in (left to right) or out (right to left)
  const zoomingIn = currentX > startX;
  
  // Update zoom box visual indicator
  this.updateZoomBox(startX, currentX, zoomingIn);
};//RelEffPlot.prototype.handleMouseMove


RelEffPlot.prototype.handleMouseUp = function() {
  if (!this.leftMouseDown) return;
  
  const currentMouse = d3.mouse(this.mouseCapture.node());
  const startX = this.leftMouseDown[0];
  const currentX = currentMouse[0];
  
  // Remove zoom box and text
  this.chartArea.select(".zoom-box").remove();
  this.chartArea.select(".zoom-text").remove();
  
  const startEnergy = this.xScale.invert(startX);
  const endEnergy = this.xScale.invert(currentX);
  
  // Determine if zooming in or out
  const zoomingIn = currentX > startX;
  
  if (Math.abs(currentX - startX) > 5) { // Minimum drag distance
    if (zoomingIn) {
      // Zoom in to the selected range
      this.setXAxisRange(Math.min(startEnergy, endEnergy), Math.max(startEnergy, endEnergy));
    } else {
      // Zoom out
      this.zoomOut();
    }
  }
  
  // Clean up
  this.leftMouseDown = null;
  this.zoomStartX = null;
  
  // Remove document listeners
  d3.select(document)
    .on("mousemove.releffzoom", null)
    .on("mouseup.releffzoom", null);
};//RelEffPlot.prototype.handleMouseUp


RelEffPlot.prototype.updateZoomBox = function(startX, currentX, zoomingIn) {
  const minX = Math.min(startX, currentX);
  const width = Math.abs(currentX - startX);
  
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
  
  // Apply appropriate CSS class and text
  if (zoomingIn) {
    zoomBox.attr("class", "zoom-box leftbuttonzoombox");
    zoomText.text("Zoom In");
  } else {
    zoomBox.attr("class", "zoom-box leftbuttonzoomoutboxy");
    zoomText.text("Zoom Out");
  }
  
  const chartAreaHeight = this.yScale.range()[0]; // Get the height from y scale range
  
  zoomBox
    .attr("x", minX)
    .attr("y", 0)
    .attr("width", width)
    .attr("height", chartAreaHeight);
    
  // Position text in center of zoom box
  zoomText
    .attr("x", minX + width / 2)
    .attr("y", chartAreaHeight / 2)
    .style("visibility", width > 50 ? "visible" : "hidden"); // Hide text if box too small
};//RelEffPlot.prototype.updateZoomBox


RelEffPlot.prototype.setXAxisRange = function(minEnergy, maxEnergy) {
  // Update x scale domain
  this.xScale.domain([minEnergy, maxEnergy]);
  
  // Update x axis
  this.chartArea.selectAll('.xAxis').call(this.xAxis);
  
  // Update y axis range based on visible data
  this.updateYAxisRange();
  
  // Redraw the data with new ranges
  if (this.datasets) {
    this.redrawData();
  }
};//RelEffPlot.prototype.setXAxisRange


RelEffPlot.prototype.zoomOut = function() {
  // Zoom out to show the full data extent
  if (!this.datasets || this.datasets.length === 0) {
    // Default range if no data
    this.setXAxisRange(0, 3000);
    return;
  }
  
  // Find min/max energy across all datasets
  let min_x = Number.MAX_VALUE;
  let max_x = Number.MIN_VALUE;
  
  this.datasets.forEach(function(dataset) {
    if (dataset.data_vals && dataset.data_vals.length > 0) {
      const data_min = d3.min(dataset.data_vals, function(d) { return d.energy; });
      const data_max = d3.max(dataset.data_vals, function(d) { return d.energy; });
      min_x = Math.min(min_x, data_min);
      max_x = Math.max(max_x, data_max);
    }
  });
  
  // Handle edge case of no data
  if (min_x === Number.MAX_VALUE) min_x = 0;
  if (max_x === Number.MIN_VALUE) max_x = 3000;
  
  // Apply the same padding logic as in setRelEffData
  const initial_x_range = max_x - min_x;
  let lower_padding = 0.025 * initial_x_range;
  if (Math.abs(initial_x_range) < 1.0) {
    lower_padding = 1.0;
  } else if (((min_x - lower_padding) < 50.0) || (lower_padding > 0.2 * min_x)) {
    lower_padding = 0.2 * min_x;
  }
  
  min_x -= lower_padding;
  max_x += 0.025 * initial_x_range;
  
  this.setXAxisRange(min_x, max_x);
};//RelEffPlot.prototype.zoomOut


RelEffPlot.prototype.updateYAxisRange = function() {
  if (!this.datasets || this.datasets.length === 0) return;
  
  const xDomain = this.xScale.domain();
  let min_y = Number.MAX_VALUE;
  let max_y = Number.MIN_VALUE;
  let hasVisibleData = false;
  
  // Find min/max y values for visible x range across all datasets
  this.datasets.forEach(function(dataset) {
    if (dataset.data_vals && dataset.data_vals.length > 0) {
      const visibleData = dataset.data_vals.filter(function(d) {
        return d.energy >= xDomain[0] && d.energy <= xDomain[1];
      });
      
      if (visibleData.length > 0) {
        hasVisibleData = true;
        const data_min_y = d3.min(visibleData, function(d) { return d.eff; });
        const data_max_y = d3.max(visibleData, function(d) { return d.eff; });
        min_y = Math.min(min_y, data_min_y);
        max_y = Math.max(max_y, data_max_y);
      }
    }
    
    // Also check fit equation points if they exist
    if (dataset.fit_eqn && dataset.fit_eqn_points) {
      const visibleFitData = dataset.fit_eqn_points.filter(function(d) {
        return d.energy >= xDomain[0] && d.energy <= xDomain[1];
      });
      
      if (visibleFitData.length > 0) {
        hasVisibleData = true;
        visibleFitData.forEach(function(d) {
          min_y = Math.min(min_y, d.eff - (d.eff_uncert ? 2 * d.eff_uncert : 0));
          max_y = Math.max(max_y, d.eff + (d.eff_uncert ? 2 * d.eff_uncert : 0));
        });
      }
    }
  });
  
  // Handle edge case of no visible data
  if (!hasVisibleData || min_y === Number.MAX_VALUE) {
    min_y = 0;
    max_y = 1;
  }
  
  // Apply padding
  const initial_y_range = max_y - min_y;
  min_y -= 0.1 * initial_y_range;
  max_y += 0.2 * initial_y_range;
  
  // Update y scale and axis
  this.yScale.domain([min_y, max_y]);
  this.chartArea.selectAll('.yAxis').call(this.yAxis);
};//RelEffPlot.prototype.updateYAxisRange


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
  
  // Redraw fit lines and uncertainty bounds
  this.datasets.forEach(function(dataset, index) {
    if (dataset.fit_eqn && dataset.fit_eqn_points && dataset.fit_eqn_points.length > 0) {
      // Create path for fit line
      const valueline = d3.svg.line()
        .x(function(d) { return self.xScale(d.energy); })
        .y(function(d) { return self.yScale(d.eff); });
        
      const path = self.chartArea.append("path")
        .attr("class", "line dataset-" + index)
        .attr("d", valueline(dataset.fit_eqn_points))
        .style("stroke", self.getDatasetColor(index));
        
      // Create path for uncertainty bounds if available
      if (dataset.fit_uncert_fcn) {
        const area = d3.svg.area()
          .x(function(d) { return self.xScale(d.energy); })
          .y0(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff - 2*d.eff_uncert) : null; })
          .y1(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff + 2*d.eff_uncert) : null; });
          
        const path_uncert = self.chartArea.append("path")
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
  if( typeof margins.right !== "number" )
    this.options.margins.right = margins.right;
  if( typeof margins.bottom !== "number" )
    this.options.margins.bottom = margins.bottom;
  if( typeof margins.left !== "number" )
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
  const self = this;

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
    
    // Calculate starting position for the x-coordinate of all texts
    const startX = parentWidth - maxTextWidth;
    
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

  let chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  let chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad;

  const num_eqn_points = chartAreaWidth / 4;

  // Find min/max values across all datasets
  let min_x = Number.MAX_VALUE;
  let max_x = Number.MIN_VALUE;
  let min_y = Number.MAX_VALUE;
  let max_y = Number.MIN_VALUE;
  
  // We'll collect all fit equation points here
  let all_fit_eqn_points = [];

  // Process all datasets to find limits and collect fit data
  datasets.forEach(function(dataset) {
    const data_vals = dataset.data_vals;
    const fit_eqn = dataset.fit_eqn;
    const fit_uncert_fcn = dataset.fit_uncert_fcn;
    
    if (data_vals && data_vals.length > 0) {
      min_x = Math.min(min_x, d3.min(data_vals, function (d) { return d.energy; }));
      max_x = Math.max(max_x, d3.max(data_vals, function (d) { return d.energy; }));
      min_y = Math.min(min_y, d3.min(data_vals, function (d) { return d.eff; }));
      max_y = Math.max(max_y, d3.max(data_vals, function (d) { return d.eff; }));
    }
    
    // Prepare fit equation points if available
    if (fit_eqn) {
      let fit_eqn_points = [];
      
      // We'll compute these later when we have the final min/max x values
      dataset.fit_eqn_points = []; // Store to use later
      all_fit_eqn_points.push({
        fit_eqn: fit_eqn,
        fit_uncert_fcn: fit_uncert_fcn,
        points: dataset.fit_eqn_points // Reference to array we'll fill later
      });
    }
  });
  
  // Handle edge case of no data
  if (min_x === Number.MAX_VALUE) min_x = 0;
  if (max_x === Number.MIN_VALUE) max_x = 3000;
  if (min_y === Number.MAX_VALUE) min_y = 0;
  if (max_y === Number.MIN_VALUE) max_y = 1;

  // Compute padding for axes
  const initial_x_range = max_x - min_x;
    
  let lower_padding = 0.025 * initial_x_range;
  if (Math.abs(initial_x_range) < 1.0) {
    lower_padding = 1.0;
  } else if (((min_x - lower_padding) < 50.0) || (lower_padding > 0.2 * min_x)) {
    lower_padding = 0.2 * min_x;
  }
    
  min_x -= lower_padding;
  max_x += 0.025 * initial_x_range;

  const initial_y_range = max_y - min_y;
  min_y -= 0.1 * initial_y_range;
  max_y += 0.2 * initial_y_range;

  // Now compute all fit equation points using final min/max x values
  all_fit_eqn_points.forEach(function(fit_data) {
    for (let i = 0; i < num_eqn_points; ++i) {
      const ene = min_x + i * (max_x - min_x) / num_eqn_points;
      const eff = fit_data.fit_eqn(ene);
      const eff_uncert = fit_data.fit_uncert_fcn ? fit_data.fit_uncert_fcn(ene) : null;
      
      fit_data.points.push({ 
        energy: ene, 
        eff: eff, 
        eff_uncert: eff_uncert 
      });
      
      // Update min/max y based on fit equations
      min_y = Math.min(min_y, eff - (eff_uncert ? 2 * eff_uncert : 0));
      max_y = Math.max(max_y, eff + (eff_uncert ? 2 * eff_uncert : 0));
    }
  });

  // Put initial estimate of areas into d3, so we can then get final sizes of axis and titles
  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);
  
  this.xScale.domain([min_x, max_x]);
  this.yScale.domain([min_y, max_y]);

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
          
  chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad;
    
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
    
  // Remove existing paths and create a group for each dataset
  this.plotGroup.selectAll("path.line").remove();
  this.plotGroup.selectAll("path.RelEffPlotErrorBounds").remove();
  
  // For each dataset with fit equation, create the fit line and uncertainty bounds
  all_fit_eqn_points.forEach(function(fit_data, index) {
    // Create path for fit line
    const valueline = d3.svg.line()
      .x(function(d) { return self.xScale(d.energy); })
      .y(function(d) { return self.yScale(d.eff); });
      
    const path = self.plotGroup.append("path")
      .attr("class", "line dataset-" + index)
      .attr("d", valueline(fit_data.points))
      .style("stroke", self.getDatasetColor(index))
      .style("pointer-events", "none");
      
    // Create path for uncertainty bounds if available
    if (fit_data.fit_uncert_fcn) {
      const area = d3.svg.area()
        .x(function(d) { return self.xScale(d.energy); })
        .y0(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff - 2*d.eff_uncert) : null; })
        .y1(function(d) { return d.eff_uncert !== null ? self.yScale(d.eff + 2*d.eff_uncert) : null; });
        
      const path_uncert = self.plotGroup.append("path")
        .attr("class", "RelEffPlotErrorBounds dataset-" + index)
        .attr("d", area(fit_data.points))
        .style("fill", self.getDatasetColor(index, 0.1))
        .style("stroke", "none")
        .style("pointer-events", "none");
    }
  });

  // For some reason the circles position arent being updated on resize, so we'll just remove them first
  //  With later versions of D3 there is a .merge() function that is maybe relevant
  this.plotGroup.selectAll("circle").remove();
  this.plotGroup.selectAll("line.errorbar").remove();
    
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
        // Check if there's a dominant nuclide with color
        if (d.nuc_info && d.nuc_info.length > 0) {
          let max_contrib = 0, dominant_color = null;
          let sum_contrib = 0;
          
          for (const el of d.nuc_info) {
            const contrib = el.rel_act * el.br;
            sum_contrib += contrib;
            if (contrib > max_contrib) {
              max_contrib = contrib;
              dominant_color = el.color;
            }
          }
          
          // Use the dominant nuclide's color if it's valid and contributes >50%
          if (dominant_color && typeof dominant_color === 'string' && 
              dominant_color.length > 0 && (max_contrib > 0.5 * sum_contrib)) {
            return dominant_color;
          }
        }
        
        // Fall back to dataset color if no valid color from nuc_info
        return self.getDatasetColor(datasetIndex);
      });
    
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
        
        if (d.nuc_info.length === 0)
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
        // Check if there's a dominant nuclide with color
        if (d.nuc_info && d.nuc_info.length > 0) {
          let max_contrib = 0, dominant_color = null;
          let sum_contrib = 0;
          
          for (const el of d.nuc_info) {
            const contrib = el.rel_act * el.br;
            sum_contrib += contrib;
            if (contrib > max_contrib) {
              max_contrib = contrib;
              dominant_color = el.color;
            }
          }

          if( !dominant_color )
            return null;  //We will rely on CSS of ".RelEffPlot circle.noiso" to handle this
          
          if( max_contrib < 0.5*sum_contrib )
            return null;  // We will rely on CSS of ".RelEffPlot circle.multiiso" to handle this
          
          // Use the dominant nuclide's color if it's valid and contributes >50%
          if (dominant_color && typeof dominant_color === 'string' && 
              dominant_color.length > 0 && (max_contrib > 0.5 * sum_contrib)) {
            return dominant_color;
          }
        }
        
        // Fall back to dataset color if no valid color from nuc_info
        return self.getDatasetColor(datasetIndex);
      })
      .on("mouseover", function (d, i) {
        self.tooltip.transition()
          .duration(200)
          .style("opacity", .9);

        d3.select(this).transition()
          .duration('50')
          .attr('opacity', '.85')
          .attr("r", 6);

        const eqn_eff = dataset.fit_eqn ? dataset.fit_eqn(d.energy) : null;

        let txt = "<div>Energy: " + (d.mean ? d.mean.toFixed(2) : d.energy.toFixed(2)) + " keV</div>"
          + "<div>Peak Area: " + d.counts.toFixed(1) + " &plusmn; " + d.counts_uncert.toFixed(1) + "</div>";


        // Add fit equation value if available for this dataset
        if (eqn_eff) {
          txt += "<div>RelEff Curve: " + eqn_eff.toPrecision(4);
          if( dataset.fit_uncert_fcn )
            txt += " &plusmn " + dataset.fit_uncert_fcn(d.energy);
          txt += "</div>";
        }

        txt += "<div>Measured RelEff: " + d.eff.toPrecision(4);
        if( d.eff_uncert && (d.eff_uncert > 0) ){
          txt += " &plusmn " + d.eff_uncert.toPrecision(4);
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
              
        self.tooltip.html(txt);

        // Make it so tooltip doesnt extend above/below/left/right of chart area
        const svg_location = d3.mouse(self.chartArea.node());
        const svg_bb = self.chartArea.node().getBoundingClientRect();
        const tt_bb = self.tooltip.node().getBoundingClientRect();
        const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
        const render_top = (svg_location[1] + tt_bb.height + 4 + 15 > svg_bb.height);
        const x_offset = render_right ? d3.event.pageX + 10 : d3.event.pageX - 10 - tt_bb.width;
        const y_offset = render_top ? d3.event.pageY - 10 - tt_bb.height : d3.event.pageY + 4;
          
        self.tooltip
          .style("left", x_offset + "px")
          .style("top", y_offset + "px");
      })
      .on('mouseout', function (d, i) {
        d3.select(this).transition()
          .duration('50')
          .attr('opacity', '1')
          .attr("r", calculateRadius(d.counts));
        self.tooltip.transition()
          .duration(500)
          .style("opacity", 0);
      });
  });
  
  // Auto-zoom out to show the full data extent when data is set
  this.zoomOut();
  
  // Position the mouse capture rectangle to cover the chart area
  if (this.mouseCapture) {
    this.mouseCapture
      .attr("width", chartAreaWidth)
      .attr("height", chartAreaHeight);
  }
};//RelEffPlot.prototype.setRelEffData

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
  
  // Update mouse capture rectangle size if it exists
  if (this.mouseCapture) {
    const parentWidth = this.chart.clientWidth;
    const parentHeight = this.chart.clientHeight;
    const chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right;
    const chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom;
    
    this.mouseCapture
      .attr("width", chartAreaWidth)
      .attr("height", chartAreaHeight);
  }
};//RelEffPlot.prototype.handleResize
