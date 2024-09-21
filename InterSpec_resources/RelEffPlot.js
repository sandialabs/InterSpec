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

  this.xScale.domain([0, 3000]);
  this.yScale.domain([0, 1]);

  // Add the valueline path.
  this.path = this.chartArea.append("path")    // Add the valueline path.
    .attr("class", "line");

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
}//RelEffPlot constructor


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


RelEffPlot.prototype.setRelEffData = function (data_vals, fit_eqn) {
  const self = this;

  if( !Array.isArray(data_vals) || (data_vals.length === 0) )
    data_vals = null;
    
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
    
  const titlePad = this.options.titleToAxisPadding;
  const xaxistitleBB = this.xaxistitle ? this.xaxistitle.node().getBBox() : null;
  const xtitleh = xaxistitleBB ? xaxistitleBB.height + titlePad : 0;
    
  const yaxistitleBB = this.yaxistitle ? this.yaxistitle.node().getBBox() : null;
  const ytitleh = yaxistitleBB ? yaxistitleBB.height + titlePad : 0;
    
  // We will perform an initial estimate of the chart are - using the axis we have.
  //  Later on we will put current values into the axis, and then get the size, and
  //  update things
  let xticks = this.chartArea.selectAll('.xAxis');
  let xtickh = xticks.empty() ? 24.5 : d3.max(xticks[0], function(t){ return d3.select(t).node().getBBox().height; });
    
  let yticks = this.chartArea.selectAll('.yAxis');
  let ytickw = yticks.empty() ? 37 : Math.max( 37, d3.max(yticks[0], function(t){ return d3.select(t).node().getBBox().width; }) ); //If we have no labels, we'll get like 7 px, so will require at least 37 px, which is the nominal value we should get
    
  let chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  let chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh;

  this.data_vals = data_vals;
  this.fit_eqn = fit_eqn;
    
  const num_eqn_points = chartAreaWidth / 4;

  // Scale the range of the data
  let min_x = data_vals ? d3.min(data_vals, function (d) { return d.energy; }) : 0;
  let max_x = data_vals ? d3.max(data_vals, function (d) { return d.energy; }) : 3000;
  const initial_x_range = max_x - min_x;
    
  let lower_padding = 0.025*initial_x_range;
  if( Math.abs(initial_x_range) < 1.0 )
    lower_padding = 1.0;
  else if( ((min_x - lower_padding) < 50.0) || (lower_padding > 0.2*min_x) )
    lower_padding = 0.2 * min_x;
    
  if( data_vals ){
    min_x -= lower_padding;
    max_x += 0.025 * initial_x_range;
  }

  let min_y = data_vals ? d3.min(data_vals, function (d) { return d.eff; }) : 0;
  let max_y = data_vals ? d3.max(data_vals, function (d) { return d.eff; }) : 1;
    
  let fit_eqn_points = [];
  if( fit_eqn ){
    for (let i = 0; i < num_eqn_points; ++i) {
      const ene = min_x + i * (max_x - min_x) / num_eqn_points;
      fit_eqn_points.push({ energy: ene, eff: fit_eqn(ene) });
    }
    min_y = Math.min(min_y, d3.min(fit_eqn_points, function (d) { return d.eff; }));
    max_y = Math.max(max_y, d3.max(fit_eqn_points, function (d) { return d.eff; }));
  }//if( fit_eqn )

  if( data_vals ){
    const initial_y_range = max_y - min_y;
    min_y -= 0.1 * initial_y_range;
    max_y += 0.2 * initial_y_range;
  }

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
  chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh;
    
  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);
       
  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);
        
  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);
               
  if( this.xaxistitle )
    this.xaxistitle
      .attr("x", this.options.margins.left + ytitleh + ytickw + 0.5*chartAreaWidth )
      .attr("y", parentHeight - this.options.margins.bottom );
         
  if( this.yaxistitle )
    this.yaxistitle
      .attr("y", this.options.margins.left )
      .attr("x",this.options.margins.top - 0.5*chartAreaHeight );
         
  this.chartArea.selectAll('.xAxis')
    .attr("transform", "translate(0," + chartAreaHeight + ")");
  this.chartArea
    .attr("transform", "translate(" + (this.options.margins.left + ytickw + ytitleh)
                       + "," + this.options.margins.top + ")");
    
  if( fit_eqn ){
    // Define the line
    const valueline = d3.svg.line()
      .x(function (d) { return self.xScale(d.energy); })
      .y(function (d) { return self.yScale(d.eff); });
    this.path
      .attr("d", valueline(fit_eqn_points));
  }else{
    this.path.attr("d", null);
  }

  // For some reason the circles position arent being updated on resize, so we'll just remove them first
  //  With later versions of D3 there is a .merge() function that is maybe relevant
  this.chartArea.selectAll("circle").remove();
  this.chartArea.selectAll("line.errorbar").remove();
    
  // Everything below here requires data points, so if we dont have any, we can skip this
  if( !data_vals )
    return;
      
  let lines = this.chartArea.selectAll("line.errorbar")
      .data(data_vals);
    
  lines.enter()
    .append('line')
    .attr("class", "errorbar")
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
    
  // Add the data points
  this.chartArea
    .selectAll("circle")
    .data(data_vals)
    .enter()
    .append("circle")
    .attr("r", 3)
    .attr("cx", function (d) {
      return self.xScale(d.energy)
    })
    .attr("cy", function (d) {
      const val = self.yScale(d.eff);
      return isNaN(val) ? 0 : val;
    })
    .attr("class", function (d) {
      if( d.nuc_info.length === 0 )
        return "noiso";

      //Return the dominant nuclide, if no nuclide is over 50%, we'll use "multiiso"
      //  TODO: use a gradient to indicate relative components
      let max_contrib = 0, sum_contrib = 0, max_nuc = null;
      for( const el of d.nuc_info ) {
        const contrib = el.rel_act * el.br;
        sum_contrib += contrib;
        if( contrib > max_contrib ){
          max_contrib = contrib;
          max_nuc = el.nuc;
        }
      }
          
      return (max_nuc && (max_contrib > 0.5*sum_contrib)) ? max_nuc : "multiiso";
    })
    .on("mouseover", function (d, i) {
      self.tooltip.transition()
        .duration(200)
        .style("opacity", .9);

      d3.select(this).transition()
        .duration('50')
        .attr('opacity', '.85')
        .attr("r", 6);

      let txt = "<div>Energy: " + (d.mean ? d.mean.toFixed(2) : d.energy.toFixed(2)) + " keV</div>"
        + "<div>Peak Area: " + d.counts.toFixed(1) + " &pm; " + d.counts_uncert.toFixed(1) + "</div>"
        + "<div>Measured RelEff: " + d.eff.toPrecision(5) + "</div>"
        + (fit_eqn ? "<div>RelEff Curve: " + fit_eqn(d.energy).toPrecision(5) + "</div>" : "");
      for (const el of d.nuc_info) {
        txt += "<div>&nbsp;&nbsp;" + el.nuc + ": br=" + el.br.toPrecision(4);
        if( el.rel_act )
          txt += ", RelAct=" + el.rel_act.toPrecision(4);
        txt += "</div>";
      }
            
      self.tooltip.html(txt);

      // Make it so tooltip doesnt extend above/below/left/right of chart area
      const svg_location = d3.mouse( self.chartArea.node() );
      const svg_bb = self.chartArea.node().getBoundingClientRect();
      const tt_bb = self.tooltip.node().getBoundingClientRect();
      const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
      const render_top = (svg_location[1] + tt_bb.height + 4 + 15 > svg_bb.height)
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
        .attr("r", 3);
      self.tooltip.transition()
        .duration(500)
        .style("opacity", 0);
    });
};//RelEffPlot.prototype.setRelEffData


RelEffPlot.prototype.handleResize = function () {
  this.setRelEffData(this.data_vals, this.fit_eqn);
};//RelEffPlot.prototype.handleResize
