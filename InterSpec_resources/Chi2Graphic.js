Chi2Graphic = function (elem,options) {
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
  if( (typeof this.options.yAxisChiTitle !== "string") || (this.options.yAxisChiTitle.length === 0) )
    this.options.yAxisChiTitle = null;
  if( (typeof this.options.yAxisScaleTitle !== "string") || (this.options.yAxisScaleTitle.length === 0) )
    this.options.yAxisScaleTitle = null;
    
  if( (this.options.displayType !== this.ChartType.Pull)
     && (this.options.displayType !== this.ChartType.Scale) )
    this.options.displayType = this.ChartType.Pull;
      
  // Set the dimensions of the canvas / graph
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  const chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right;
  const chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom;

  // Setup the tooltip
  this.tooltip = d3.select(this.chart).append("div")
    .attr("class", "Chi2GraphicTooltip")
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
    .attr("class", "Chi2Graphic");
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
      .attr("class", "xAxis " + (this.options.displayType === this.ChartType.Pull ? "Pull" : "Scale") )
      .attr("transform", "translate(0," + chartAreaHeight + ")")
      .call(this.xAxis);

  // Add the Y Axis
  this.chartArea.append("g")
    .attr("class", "yAxis")
    .call(this.yAxis);
      
  if( this.options.xAxisTitle )
    this.setXAxisTitle( this.options.xAxisTitle, true );
    
  this.setYAxisTitles( this.options.yAxisChiTitle, this.options.yAxisScaleTitle, true );
}//Chi2Graphic constructor

Chi2Graphic.prototype.ChartType = Object.freeze({
  Pull: 0,
  Scale: 1
} );

Chi2Graphic.prototype.setMargins = function( margins ){
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
}//Chi2Graphic.prototype.setMargins


Chi2Graphic.prototype.setXAxisTitle = function( title, dontCallResize ){
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
    this.xaxistitle = this.chartArea.append("text")
      .attr("class", "xaxistitle")
      .attr("x", 0.5*this.chart.clientWidth ) //roughly position title - will be updated when data set or resized
      .attr("y", this.chart.clientHeight - this.options.margins.bottom )
      .style("text-anchor", "middle");
  }
  
  this.xaxistitle.text(title);
  
  if( !dontCallResize )
    this.handleResize();
}//Chi2Graphic.prototype.setXAxisTitle


Chi2Graphic.prototype.setYAxisTitles = function( chi_title, scale_title, dontCallResize ){
  if( !chi_title || (chi_title.length === 0) )
    this.options.yAxisChiTitle = null;
  if( !scale_title || (scale_title.length === 0) )
    this.options.yAxisScaleTitle = null;
  
  const isPullChart = (this.ChartType.Pull === this.options.displayType);
  const title = isPullChart ? this.options.yAxisChiTitle : this.options.yAxisScaleTitle;
  
  if( !title ){
    this.yaxistitle.remove();
    this.yaxistitle = null;
    return;
  }
  
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
}//Chi2Graphic.prototype.setYAxisTitles

Chi2Graphic.prototype.setDisplayType = function( type ) {
  this.options.displayType = type;
  
  console.assert( (this.ChartType.Pull === type) || (this.ChartType.Scale === type),
    "Invalid value passed to Chi2Graphic.setDisplayType" );
  
  const isPullChart = (this.ChartType.Pull === type);
  this.chartArea.selectAll('.xAxis').attr("class", "xAxis " + (isPullChart ? "Pull" : "Scale") );
  
  this.setYAxisTitles( this.options.yAxisChiTitle, this.options.yAxisScaleTitle, true );
  
  this.chartArea.selectAll(".ScaleUnityLine").remove();
  if( !isPullChart ){
    this.chartArea
    .append('line')
    .attr("class", "ScaleUnityLine")
    .attr('x1', 0 )
    .attr('x2', 0 )
    .attr('y1', 0 )
    .attr('y2', 0 );
  }
  
  this.setData( this.data );
}//Chi2Graphic.prototype.setDisplayType


Chi2Graphic.prototype.setData = function( data ) {
  
  /*
   {
     "NDOF" : 1,
     "Points" : [
       {
         "color" : "rgb(0,185,255)",
         "energy" : 238.632003784179688,
         "nuclide" : "Th232",
         "numExpected" : 36187.425109183175664,
         "numObserved" : 38296.61328125,
         "numObservedUncert" : 204.2442626953125,
         "numSigmaOff" : 10.32679275409204,
         "observedOverExpected" : 1.058285113287366,
         "observedOverExpectedUncert" : 0.005644067299043
       },
       {
         "color" : "rgb(0,185,255)",
         "energy" : 338.32000732421875,
         "nuclide" : "Th232",
         "numExpected" : 7730.95941759361358,
         "numObserved" : 7781.0693359375,
         "numObservedUncert" : 97.467079162597656,
         "numSigmaOff" : 0.514121473367346,
         "observedOverExpected" : 1.006481720526155,
         "observedOverExpectedUncert" : 0.012607371724238
       },
       ...
     ]
   }
   
   */
  
  const self = this;
  this.data = data;
  
  self.tooltip.style("opacity", 0);
  
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  
  const isPullChart = (this.ChartType.Pull === this.options.displayType);
  
  const valid_data = (data && Array.isArray(data.Points) && data.Points.length
                             && (typeof data.Points[0].energy === "number" ));
  
  
  if( !this.chartInfoTitle ){
    this.chartInfoTitle = this.svg.append("text")
      .attr("class", "ChiInfoTitle")
      .attr("y",0);
  }
  
  let chi2Txt = "";
  if( valid_data ){
    let numsigma = 0.0, numpoints = 0;
    for( let i = 0; i < data.Points.length; ++i ){
      numsigma += Math.abs(data.Points[i].numSigmaOff);
      numpoints += 1;
    }
    chi2Txt = "<dev>=" + (numsigma/numpoints).toFixed(2) + "σ";
  }
  
  this.chartInfoTitle
    .attr("x", parentWidth - 15)
    .text(chi2Txt);
  
  let min_x = valid_data ? d3.min(data.Points, function (d) { return d.energy; }) : 0;
  let max_x = valid_data ? d3.max(data.Points, function (d) { return d.energy; }) : 3000;
  const initial_x_range = max_x - min_x;
    
  const yvaluefcn = function (d) { return isPullChart ? d.numSigmaOff : d.observedOverExpected; };
    
  let min_y = valid_data ? d3.min(data.Points, yvaluefcn ): (isPullChart ? -1 : 0);
  let max_y = valid_data ? d3.max(data.Points, yvaluefcn ): 1;
  
  const titlePad = this.options.titleToAxisPadding;
  const xaxistitleBB = this.xaxistitle ? this.xaxistitle.node().getBBox() : null;
  const xtitleh = isPullChart ? 0 : xaxistitleBB ? xaxistitleBB.height + titlePad : 0;
    
  const yaxistitleBB = this.yaxistitle ? this.yaxistitle.node().getBBox() : null;
  const ytitleh = yaxistitleBB ? yaxistitleBB.height + titlePad : 0;
  
  // TODO: Set this.yaxistitle to default font size, then if it is larger than the chart, and if so, resize the text, and if not
  
  // We will perform an initial estimate of the chart are - using the axis we have.
  //  Later on we will put current values into the axis, and then get the size, and
  //  update things
  let xticks = this.chartArea.selectAll('.xAxis');
  let xtickh = xticks.empty() ? 24.5 : d3.max(xticks[0], function(t){ return d3.select(t).node().getBBox().height; });
    
  let yticks = this.chartArea.selectAll('.yAxis');
  let ytickw = yticks.empty() ? 37 : Math.max( 37, d3.max(yticks[0], function(t){ return d3.select(t).node().getBBox().width; }) ); //If we have no labels, we'll get like 7 px, so will require at least 37 px, which is the nominal value we should get
    
  const chi2_txt_pad = this.chartInfoTitle ? 0.5*this.chartInfoTitle.node().getBBox().height : 0;
    
  let chartAreaWidth = parentWidth - this.options.margins.left - this.options.margins.right - ytitleh - ytickw;
  let chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad;

  let lower_padding = 0.025*initial_x_range;
  if( Math.abs(initial_x_range) < 1.0 )
    lower_padding = 1.0;
  else if( ((min_x - lower_padding) < 50.0) || (lower_padding > 0.2*min_x) )
    lower_padding = 0.2 * min_x;
    
  if( valid_data ){
    min_x -= lower_padding;
    max_x += 0.025 * initial_x_range;
  }
  
  //console.log( "Initial y :{", min_y, ",", max_y, "}" );
  if( isPullChart ){
    min_y = Math.min(-0.25, min_y);
    max_y = Math.max(0.25, max_y);
    if( valid_data ){
      const dy = Math.max( 0.5, max_y - min_y );
      min_y -= 0.15*dy;
      max_y += 0.15*dy;
    }
  }else{
    min_y = 0.0;
    if( valid_data ){
      max_y += 0.2 * max_y;
    }else{
      max_y = 2.0;
    }
    max_y = Math.max( max_y, 1.5 );
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
  chartAreaHeight = parentHeight - this.options.margins.top - this.options.margins.bottom - xtitleh - xtickh - chi2_txt_pad;
    
  this.xScale.range([0, chartAreaWidth]);
  this.yScale.range([chartAreaHeight, 0]);
       
  this.chartArea.selectAll('.xAxis')
    .call(this.xAxis);
        
  this.chartArea.selectAll('.yAxis')
    .call(this.yAxis);
   
  if( this.yaxistitle )
    this.yaxistitle
      .attr("y", this.options.margins.left )
      .attr("x",this.options.margins.top + chi2_txt_pad - 0.5*chartAreaHeight );
         
  this.chartArea.selectAll('.xAxis')
    .attr("transform", "translate(0," + chartAreaHeight + ")");
  this.chartArea
    .attr("transform", "translate(" + (this.options.margins.left + ytickw + ytitleh)
                       + "," + (this.options.margins.top + chi2_txt_pad) + ")");

  // For some reason the circles position arent being updated on resize, so we'll just remove them first
  //  With later versions of D3 there is a .merge() function that is maybe relevant
  this.chartArea.selectAll("circle").remove();
  this.chartArea.selectAll("line.errorbar").remove();
    
  if( !isPullChart ){
    const y1val = self.yScale(1);
    
    this.chartArea.selectAll(".ScaleUnityLine")
    .attr('x1', 0 )
    .attr('x2', chartAreaWidth )
    .attr('y1', y1val )
    .attr('y2', y1val );
  }
  
  const set_yaxis_info = function(){
    const xaxis_yloc = self.yScale( 0 );
    self.chartArea.selectAll('.xAxis').attr("transform", "translate(0," + xaxis_yloc + ")");
     
    if( self.xaxistitle )
      self.xaxistitle
        .attr("x", self.options.margins.left + ytitleh + ytickw + 0.5*chartAreaWidth )
        .attr("y", xaxis_yloc + xtickh + titlePad );
  };
  
  set_yaxis_info();
    
  // Everything below here requires data points, so if we dont have any, we can skip this
  if( !valid_data )
    return;
      
  let lines = this.chartArea.selectAll("line.errorbar")
      .data(data.Points);
    
  if( !isPullChart )
  {
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
  }//if( !isPullChart )
  
  
  // Add the data points
  this.chartArea
    .selectAll("circle")
    .data(data.Points)
    .enter()
    .append("circle")
    .attr("r", 5)
    .attr("fill", function (d) {
      return d.color ? d.color : "blue";
    })
    .attr("cx", function (d) {
      return self.xScale(d.energy)
    })
    .attr("cy", function (d) {
      const val = self.yScale( yvaluefcn(d) );
      return isNaN(val) ? 0 : val;
    })
    .attr("class", function (d) {
      return isPullChart ? "PullDataPoint" : "ScaleDataPoint";
    })
    .on("mouseover", function (d, i) {
      self.tooltip.transition()
        .duration(200)
        .style("opacity", .9);

      d3.select(this).transition()
        .duration('50')
        .attr('opacity', '.85')
        .attr("r", 7);
        
      let txt = "<div>Energy: " + d.energy.toFixed(2) + " keV</div>"
        + (d.numObserved ? "<div>Num Observed: " + d.numObserved.toFixed(1) + " &pm; " + d.numObservedUncert.toFixed(1) + "</div>" : "")
        + (d.numExpected ? "<div>Num Expected: " + d.numExpected.toFixed(1) + "</div>" : "")
        + (d.numBackground ? "<div>Num Background: " + d.numBackground.toFixed(1) + " &pm; " + d.numBackgroundUncert.toFixed(1) + "</div>" : "")
        + "<div>Num Sigma Off: " + d.numSigmaOff.toPrecision(5) + "</div>"
        + "<div>Obs/Exp: " + d.observedOverExpected.toPrecision(5) + " &pm; " + d.observedOverExpectedUncert.toPrecision(5) + "</div>"
        + "<div>Nuclide: " + d.nuclide + "</div>"
        + "</div>";
            
      self.tooltip.html(txt);

      // Make it so tooltip doesnt extend above/below/left/right of chart area
      const svg_location = d3.mouse( self.chartArea.node() );
      const svg_bb = self.chartArea.node().getBoundingClientRect();
      const tt_bb = self.tooltip.node().getBoundingClientRect();
      const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
      const render_top = ((svg_location[1] + tt_bb.height + 4 + 15) > svg_bb.height);
      
      const x_offset = render_right ? d3.event.pageX + 10 : d3.event.pageX - 10 - tt_bb.width;
      const y_offset = render_top
                        ? (((d3.event.pageY - 10 - tt_bb.height) < 0) ? 4 : (d3.event.pageY - 10 - tt_bb.height))
                        : ((d3.event.pageY + tt_bb.height + 4) > svg_bb.height) ? (svg_bb.height - tt_bb.height) : (d3.event.pageY + 4);
        
      self.tooltip
        .style("left", x_offset + "px")
        .style("top", y_offset + "px");
    })
    .on("mouseout", function (d, i) {
      d3.select(this).transition()
        .duration('50')
        .attr('opacity', '1')
        .attr("r", 5);
      self.tooltip.transition()
        .duration(500)
        .style("opacity", 0);
    });
};//Chi2Graphic.prototype.setData


Chi2Graphic.prototype.handleResize = function () {
  this.setData( this.data );
};//Chi2Graphic.prototype.handleResize
