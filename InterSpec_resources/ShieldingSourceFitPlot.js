ShieldingSourceFitPlot = function (elem, options) {
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

  // Initialize display mode
  this.showChi = (typeof this.options.showChi === "boolean") ? this.options.showChi : true;

  // Initialize localized strings
  this.xAxisTitle = (typeof this.options.xAxisTitle === "string") ? this.options.xAxisTitle : "Energy (keV)";
  this.yAxisTitleChi = (typeof this.options.yAxisTitleChi === "string") ? this.options.yAxisTitleChi : "(observed-model)/uncert";
  this.yAxisTitleMult = (typeof this.options.yAxisTitleMult === "string") ? this.options.yAxisTitleMult : "Mult. of Model";
  this.tooltipChi = (typeof this.options.tooltipChi === "string") ? this.options.tooltipChi : "Chi-squared residual";
  this.tooltipMult = (typeof this.options.tooltipMult === "string") ? this.options.tooltipMult : "Multiple of model";

  // Set the dimensions of the canvas / graph
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  const chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right;
  const chartAreaHeight = Math.max(1, parentHeight - this.options.margins.top - this.options.margins.bottom);

  // Setup the tooltip
  this.tooltip = d3.select(this.chart).append("div")
    .attr("class", "ShieldingSourceFitPlotTooltip")
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
    .innerTickSize(7)
    .outerTickSize(1)
    .ticks(5);

  // Minor ticks axis - more ticks, no labels, shorter tick marks
  this.xAxisMinor = d3.svg
    .axis()
    .scale(this.xScale)
    .orient("bottom")
    .tickFormat("")
    .innerTickSize(3)
    .outerTickSize(0)
    .ticks(25);

  this.yAxis = d3.svg
    .axis()
    .scale(this.yScale)
    .tickFormat( function(x){
      if( x == 0 )
        return "0";
      if( Math.abs(x) < 0.01 )
        return x.toExponential(1);
      if( Math.abs(x) < 0.1 )
        return x.toFixed(3);
      if( Math.abs(x) <= 1.0 )
        return x.toFixed(2);
      if( Math.abs(x) <= 10.0 )
        return x.toFixed(1);
      return x.toFixed(0);
    } )
    .orient("left")
    .innerTickSize(7)
    .outerTickSize(1)
    .ticks(5);

  // Adds the svg canvas
  this.svg = d3.select(this.chart)
    .append("svg")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("class", "ShieldingSourceFitPlot");

  this.chartArea = this.svg  //This holds everything except the optional x and y axis titles
    .append("g")
    .attr("transform", "translate(" + this.options.margins.left + "," + this.options.margins.top + ")");

  // Create clipping path to prevent elements from being drawn outside the chart area
  this.chartArea.append("defs")
    .append("clipPath")
    .attr("id", "shieldingfitplot-clip-" + this.chart.id)
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", Math.max(0, chartAreaWidth))
    .attr("height", Math.max(0, chartAreaHeight));

  this.xScale.domain([0, 3000]);
  this.yScale.domain([-2, 2]);

  // Create a clipped group for all plot elements (paths, circles, error bars, reference lines)
  this.plotGroup = this.chartArea.append("g")
    .attr("clip-path", "url(#shieldingfitplot-clip-" + this.chart.id + ")");

  // Add the X Axis (major ticks)
  this.chartArea.append("g")
      .attr("class", "xAxis")
      .attr("transform", "translate(0," + chartAreaHeight + ")")
      .call(this.xAxis);

  // Add the X Axis minor ticks
  this.chartArea.append("g")
      .attr("class", "xAxisMinor")
      .attr("transform", "translate(0," + chartAreaHeight + ")")
      .call(this.xAxisMinor);

  // Add the Y Axis
  this.chartArea.append("g")
    .attr("class", "yAxis")
    .call(this.yAxis);

  // Add y-axis clickable area for mode toggling
  const self = this;
  this.yAxisClickArea = this.chartArea.append("rect")
    .attr("class", "yaxisarea")
    .attr("x", -60)
    .attr("y", 0)
    .attr("width", 60)
    .attr("height", chartAreaHeight)
    .style("fill", "transparent")
    .style("cursor", "pointer")
    .on("click", function() {
      self.toggleDisplayMode();
    });

  // Add axis titles if provided
  if( this.xAxisTitle && this.xAxisTitle.length > 0 )
    this.setXAxisTitle( this.xAxisTitle, true );

  const yAxisTitle = this.showChi ? this.yAxisTitleChi : this.yAxisTitleMult;
  if( yAxisTitle && yAxisTitle.length > 0 )
    this.setYAxisTitle( yAxisTitle, true );
};//ShieldingSourceFitPlot constructor


ShieldingSourceFitPlot.prototype.setLocalizations = function( locStrings ) {
  if( !locStrings )
    return;

  if( typeof locStrings.xAxisTitle === "string" )
    this.xAxisTitle = locStrings.xAxisTitle;
  if( typeof locStrings.yAxisTitleChi === "string" )
    this.yAxisTitleChi = locStrings.yAxisTitleChi;
  if( typeof locStrings.yAxisTitleMult === "string" )
    this.yAxisTitleMult = locStrings.yAxisTitleMult;
  if( typeof locStrings.tooltipChi === "string" )
    this.tooltipChi = locStrings.tooltipChi;
  if( typeof locStrings.tooltipMult === "string" )
    this.tooltipMult = locStrings.tooltipMult;

  // Update axis titles
  if( this.xaxistitle )
    this.xaxistitle.text( this.xAxisTitle );

  if( this.yaxistitle ) {
    const yAxisTitle = this.showChi ? this.yAxisTitleChi : this.yAxisTitleMult;
    this.yaxistitle.text( yAxisTitle );
  }
};//ShieldingSourceFitPlot.prototype.setLocalizations


ShieldingSourceFitPlot.prototype.setXAxisTitle = function( title, dontCallResize ){
  if( !title || (title.length === 0) )
  {
    this.xAxisTitle = null;
    if( this.xaxistitle ){
      this.xaxistitle.remove();
      this.xaxistitle = null;
    }
    return;
  }//if( we are getting rid of the title )

  this.xAxisTitle = title;
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
};//ShieldingSourceFitPlot.prototype.setXAxisTitle


ShieldingSourceFitPlot.prototype.setYAxisTitle = function( title, dontCallResize ){
  if( !title || (title.length === 0) )
  {
    if( this.yaxistitle ){
      this.yaxistitle.remove();
      this.yaxistitle = null;
    }
    return;
  }//if( we are getting rid of the title )

  const self = this;
  if( !this.yaxistitle ){
    this.yaxistitle = this.svg.append("text")
      .attr("class", "yaxistitle")
      .attr("transform", "rotate(-90)")
      .attr("y", -this.options.margins.left )
      .attr("x", -0.5*this.chart.clientHeight )
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .style("dominant-baseline", "text-after-edge")
      .style("cursor", "pointer")
      .on("click", function() {
        self.toggleDisplayMode();
      });
  }

  this.yaxistitle.text(title);

  if( !dontCallResize )
    this.handleResize();
};//ShieldingSourceFitPlot.prototype.setYAxisTitle


ShieldingSourceFitPlot.prototype.setMargins = function( margins ){
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
};//ShieldingSourceFitPlot.prototype.setMargins


ShieldingSourceFitPlot.prototype.setShowChi = function( showChi ) {
  this.showChi = showChi;

  // Update y-axis title
  const yAxisTitle = this.showChi ? this.yAxisTitleChi : this.yAxisTitleMult;
  this.setYAxisTitle( yAxisTitle, true );

  // Redraw the chart
  if( this.data )
    this.setData( this.data );
};//ShieldingSourceFitPlot.prototype.setShowChi


ShieldingSourceFitPlot.prototype.toggleDisplayMode = function() {
  this.setShowChi( !this.showChi );

  // Emit signal to C++ if Wt.emit exists
  if( typeof Wt !== 'undefined' && typeof Wt.emit === 'function' )
    Wt.emit( this.chart, 'displayModeChanged', this.showChi );
};//ShieldingSourceFitPlot.prototype.toggleDisplayMode


ShieldingSourceFitPlot.prototype.setData = function( data ) {
  const self = this;

  if( !data || !data.data_points || data.data_points.length === 0 ) {
    // Clear the chart
    this.plotGroup.selectAll("*").remove();
    this.data = null;
    return;
  }

  this.data = data;
  const data_points = data.data_points;

  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;

  const titlePad = this.options.titleToAxisPadding;
  const xaxistitleBB = this.xaxistitle ? this.xaxistitle.node().getBBox() : null;
  const xtitleh = xaxistitleBB ? xaxistitleBB.height + titlePad : 0;

  const yaxistitleBB = this.yaxistitle ? this.yaxistitle.node().getBBox() : null;
  const ytitleh = yaxistitleBB ? yaxistitleBB.height + titlePad : 0;

  // We will perform an initial estimate of the chart area
  let xticks = this.chartArea.selectAll('.xAxis');
  let xtickh = xticks.empty() ? 24.5 : d3.max(xticks[0], function (t) { return d3.select(t).node().getBBox().height; });

  let yticks = this.chartArea.selectAll('.yAxis');
  let ytickw = yticks.empty() ? 37 : Math.max(37, d3.max(yticks[0], function (t) { return d3.select(t).node().getBBox().width; }));

  let chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  let chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh;

  // Calculate data ranges
  const min_x = d3.min(data_points, d => d.energy);
  const max_x = d3.max(data_points, d => d.energy);

  let min_y, max_y;

  if( this.showChi ) {
    // Chi mode: include error bars and ensure y=0, y=-1, and y=1 are visible
    min_y = d3.min(data_points, d => d.chi);
    max_y = d3.max(data_points, d => d.chi);

    // Ensure at least [-1, 1] range for Chi mode
    min_y = Math.min(min_y, -1.0);
    max_y = Math.max(max_y, 1.0);

    // Ensure y=0 is visible
    min_y = Math.min(min_y, 0);
    max_y = Math.max(max_y, 0);
  } else {
    // Mult mode: include error bars and ensure y=0 and y=1 are visible
    min_y = d3.min(data_points, d => d.mult - d.mult_uncert);
    max_y = d3.max(data_points, d => d.mult + d.mult_uncert);

    // Ensure y=1 is visible
    min_y = Math.min(min_y, 1.0);
    max_y = Math.max(max_y, 1.0);

    // Ensure y=0 is visible
    min_y = Math.min(min_y, 0);
    max_y = Math.max(max_y, 0);
  }

  // Add padding
  const x_range = max_x - min_x;
  const y_range = max_y - min_y;

  const x_padding = 0.05 * x_range;
  const y_padding = 0.1 * y_range;

  this.xScale.domain([min_x - x_padding, max_x + x_padding]);
  this.yScale.domain([min_y - y_padding, max_y + y_padding]);

  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);

  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);

  this.chartArea.selectAll('.xAxisMinor')
    .call(this.xAxisMinor);

  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);

  // Now that we have x and y axis tick values, recompute sizes
  xticks = this.chartArea.selectAll('.xAxis');
  xtickh = xticks.empty() ? 7 : d3.max(xticks[0], function(t){ return d3.select(t).node().getBBox().height; });

  yticks = this.chartArea.selectAll('.yAxis');
  ytickw = yticks.empty() ? 7 : d3.max(yticks[0], function(t){ return d3.select(t).node().getBBox().width; });

  chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh;

  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);

  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);

  this.chartArea.selectAll('.xAxisMinor')
    .call(this.xAxisMinor);

  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);

  if( this.xaxistitle )
    this.xaxistitle
      .attr("x", this.options.margins.left + ytitleh + ytickw + 0.5*chartAreaWidth)
      .attr("y", parentHeight - this.options.margins.bottom);

  if( this.yaxistitle )
    this.yaxistitle
      .attr("y", this.options.margins.left)
      .attr("x", this.options.margins.top - 0.5*chartAreaHeight);

  this.chartArea.selectAll('.xAxis')
    .attr("transform", "translate(0," + chartAreaHeight + ")");
  this.chartArea.selectAll('.xAxisMinor')
    .attr("transform", "translate(0," + chartAreaHeight + ")");
  this.chartArea
    .attr("transform", "translate(" + (this.options.margins.left + ytickw + ytitleh)
                       + "," + this.options.margins.top + ")");

  // Update clipping rectangle dimensions
  this.chartArea.select("#shieldingfitplot-clip-" + this.chart.id + " rect")
    .attr("width", chartAreaWidth)
    .attr("height", chartAreaHeight);

  // Update y-axis click area
  if( this.yAxisClickArea )
    this.yAxisClickArea.attr("height", chartAreaHeight);

  // Clear previous plot elements
  this.plotGroup.selectAll("*").remove();

  // Add reference line (at y=0 for Chi mode, y=1 for Mult mode) - behind data points
  const refLineY = this.showChi ? 0.0 : 1.0;
  this.plotGroup.append("line")
    .attr("class", "refline")
    .attr("x1", this.xScale(min_x - x_padding))
    .attr("x2", this.xScale(max_x + x_padding))
    .attr("y1", this.yScale(refLineY))
    .attr("y2", this.yScale(refLineY))
    .style("stroke", "grey")
    .style("stroke-width", 1)
    .style("stroke-dasharray", "5,5");

  // Add error bars
  data_points.forEach(function(d) {
    if( !self.showChi && d.mult_uncert > 0 ) {
      self.plotGroup.append("line")
        .attr("class", "errorbar")
        .attr("x1", self.xScale(d.energy))
        .attr("x2", self.xScale(d.energy))
        .attr("y1", self.yScale(d.mult + d.mult_uncert))
        .attr("y2", self.yScale(d.mult - d.mult_uncert))
        .style("stroke", d.color);
    }
  });

  // Add data points
  this.plotGroup.selectAll("circle")
    .data(data_points)
    .enter()
    .append("circle")
    .attr("r", 4)
    .attr("cx", d => self.xScale(d.energy))
    .attr("cy", d => self.yScale(self.showChi ? d.chi : d.mult))
    .style("fill", d => d.color)
    .style("cursor", "pointer")
    .on("mouseover", function(d) {
      self.tooltip.transition()
        .duration(200)
        .style("opacity", .9);

      d3.select(this).transition()
        .duration('50')
        .attr('opacity', '.85')
        .attr("r", 6);

      let txt = "<div><b>" + d.nuclide + "</b> @ " + d.energy.toFixed(2) + " keV</div>";
      txt += "<div>Chi: " + d.chi.toFixed(3) + " σ</div>";
      txt += "<div>Mult: " + d.mult.toFixed(3);
      if( d.mult_uncert > 0 )
        txt += " ± " + d.mult_uncert.toFixed(3);
      txt += "</div>";

      if( d.sources && d.sources.length > 1 ) {
        txt += "<div><b>Contributions:</b></div>";
        for( const src of d.sources ) {
          txt += "<div>&nbsp;&nbsp;" + src.nuclide + ": " + (src.fraction * 100).toFixed(1) + "%</div>";
        }
      }

      self.tooltip.html(txt);

      // Position tooltip
      const svg_location = d3.mouse(self.chartArea.node());
      const svg_bb = self.chartArea.node().getBoundingClientRect();
      const tt_bb = self.tooltip.node().getBoundingClientRect();
      const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
      const render_top = (svg_location[1] + tt_bb.height + 4 + 15 > svg_bb.height);
      // Use clientX/clientY instead of pageX/pageY since tooltip uses position:fixed (viewport-relative)
      const x_offset = render_right ? d3.event.clientX + 10 : d3.event.clientX - 10 - tt_bb.width;
      const y_offset = render_top ? d3.event.clientY - 10 - tt_bb.height : d3.event.clientY + 4;

      self.tooltip
        .style("left", x_offset + "px")
        .style("top", y_offset + "px");
    })
    .on('mouseout', function(d) {
      d3.select(this).transition()
        .duration('50')
        .attr('opacity', '1')
        .attr("r", 4);
      self.tooltip.transition()
        .duration(500)
        .style("opacity", 0);
    })
    .on('click', function(d) {
      // Emit signal to C++ if Wt.emit exists
      if( typeof Wt !== 'undefined' && typeof Wt.emit === 'function' )
        Wt.emit( self.chart, 'dataPointClicked', d.energy );
    });

  // Calculate and display <dev> statistic
  let sum_chi_sq = 0;
  for( const d of data_points ) {
    sum_chi_sq += d.chi * d.chi;
  }
  const dev = Math.sqrt( sum_chi_sq / data_points.length );

  // Remove old dev text if it exists
  this.svg.selectAll(".devtext").remove();

  // Add dev text in upper right
  this.svg.append("text")
    .attr("class", "devtext")
    .attr("x", parentWidth - 10)
    .attr("y", 15)
    .style("text-anchor", "end")
    .style("font-size", "12px")
    .text("<dev>=" + dev.toFixed(2) + " σ");
};//ShieldingSourceFitPlot.prototype.setData


ShieldingSourceFitPlot.prototype.handleResize = function() {
  if( this.data )
    this.setData( this.data );
};//ShieldingSourceFitPlot.prototype.handleResize
