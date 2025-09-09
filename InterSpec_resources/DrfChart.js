DrfChart = function (elem, options) {
  const self = this;
  
  this.chart = typeof elem === "string" ? document.getElementById(elem) : elem;
  this.options = options ? options : {};
  
  // Set default margins
  if (!this.options.margins)
    this.options.margins = {};
  if (typeof this.options.margins.top !== "number")
    this.options.margins.top = 5;
  if (typeof this.options.margins.right !== "number")
    this.options.margins.right = 50; // Extra space for right Y axis (when present)
  if (typeof this.options.margins.bottom !== "number")
    this.options.margins.bottom = 40;
  if (typeof this.options.margins.left !== "number")
    this.options.margins.left = 60;

  // Initialize chart dimensions
  this.updateDimensions();

  // Setup the tooltip
  this.tooltip = d3.select(this.chart).append("div")
    .attr("class", "DrfChartTooltip")
    .style("opacity", 0);

  // Create the SVG
  this.svg = d3.select(this.chart)
    .append("svg")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("class", "DrfChart");

  // Create the chart area
  this.chartArea = this.svg.append("g")
    .attr("transform", `translate(${this.options.margins.left}, ${this.options.margins.top})`);
    
  // Add chart area background rectangle (will be sized properly later)
  this.chartAreaBg = this.chartArea.append("rect")
    .attr("class", "chart-area")
    .attr("x", 0)
    .attr("y", 0);

  // Initialize scales
  this.xScale = d3.scale.linear()
    .domain([0, 3000])
    .range([0, this.chartAreaWidth]);
    
  this.efficiencyScale = d3.scale.linear()
    .domain([0, 1])
    .range([this.chartAreaHeight, 0]);
    
  this.fwhmScale = d3.scale.linear()
    .domain([0, 10])
    .range([this.chartAreaHeight, 0]);

  // Create axes
  this.xAxis = d3.svg.axis()
    .scale(this.xScale)
    .orient("bottom")
    .tickFormat(d => `${d.toFixed(0)}`);

  this.leftYAxis = d3.svg.axis()
    .scale(this.efficiencyScale)
    .orient("left")
    .tickFormat(d => {
      if (d === 0) return "0";
      if (d <= 1.0) return d.toFixed(2);
      if (d <= 10.0) return d.toFixed(1);
      return d.toFixed(0);
    });

  this.rightYAxis = d3.svg.axis()
    .scale(this.fwhmScale)
    .orient("right")
    .tickFormat(d => {
      if (d === 0) return "0";
      if (d <= 10.0) return d.toFixed(1);
      return d.toFixed(0);
    });

  // Add axes to chart
  this.xAxisGroup = this.chartArea.append("g")
    .attr("class", "x-axis")
    .attr("transform", `translate(0, ${this.chartAreaHeight})`)
    .call(this.xAxis);

  this.leftYAxisGroup = this.chartArea.append("g")
    .attr("class", "y-axis left")
    .call(this.leftYAxis);

  this.rightYAxisGroup = this.chartArea.append("g")
    .attr("class", "y-axis right")
    .attr("transform", `translate(${this.chartAreaWidth}, 0)`)
    .call(this.rightYAxis);

  // Add axis labels
  this.xAxisLabel = this.svg.append("text")
    .attr("class", "axis-label x-axis-label")
    .attr("text-anchor", "middle")
    .attr("x", this.options.margins.left + this.chartAreaWidth / 2)
    .attr("y", this.svg.node().getBoundingClientRect().height - 5)
    .text("Energy (keV)");

  this.leftYAxisLabel = this.svg.append("text")
    .attr("class", "axis-label y-axis-label left")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(-90)")
    .attr("x", -(this.options.margins.top + this.chartAreaHeight / 2))
    .attr("y", 15)
    .text("Efficiency");

  this.rightYAxisLabel = this.svg.append("text")
    .attr("class", "axis-label y-axis-label right")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(90)")
    .attr("x", this.options.margins.top + this.chartAreaHeight / 2)
    .attr("y", -(this.svg.node().getBoundingClientRect().width - 15))
    .text("FWHM (keV)");

  // Initialize line generators
  this.efficiencyLine = d3.svg.line()
    .x(d => this.xScale(d.energy))
    .y(d => this.efficiencyScale(d.efficiency));

  this.fwhmLine = d3.svg.line()
    .x(d => this.xScale(d.energy))
    .y(d => this.fwhmScale(d.fwhm));

  // Create clipping path for chart area
  this.chartArea.append("defs")
    .append("clipPath")
    .attr("id", "drfchart-clip-" + this.chart.id)
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);

  // Create group for clipped content
  this.plotGroup = this.chartArea.append("g")
    .attr("clip-path", "url(#drfchart-clip-" + this.chart.id + ")");

  // Add paths for lines
  this.efficiencyPath = this.plotGroup.append("path")
    .attr("class", "efficiency-line")
    .style("fill", "none")
    .style("stroke-width", "2px");

  this.fwhmPath = this.plotGroup.append("path")
    .attr("class", "fwhm-line") 
    .style("fill", "none")
    .style("stroke-width", "2px");

  // Initialize mouse interaction variables
  this.leftMouseDown = null;
  this.zooming = false;
  
  // Setup mouse interactions
  this.setupMouseInteractions();
  
  // Initialize data storage
  this.data = null;
};

DrfChart.prototype.updateDimensions = function() {
  const parentWidth = this.chart.clientWidth;
  const parentHeight = this.chart.clientHeight;
  this.chartAreaWidth = Math.max(0, parentWidth - this.options.margins.left - this.options.margins.right);
  this.chartAreaHeight = Math.max(0, parentHeight - this.options.margins.top - this.options.margins.bottom);
  
  // Update chart area background dimensions if it exists
  if (this.chartAreaBg) {
    this.chartAreaBg
      .attr("width", this.chartAreaWidth)
      .attr("height", this.chartAreaHeight);
  }
};

DrfChart.prototype.setupMouseInteractions = function() {
  const self = this;
  
  // Add mouse interaction to the chart area
  this.chartArea.append("rect")
    .attr("class", "mouse-capture")
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight)
    .style("fill", "none")
    .style("pointer-events", "all")
    .on("mousedown", function() {
      const mouse = d3.mouse(this);
      self.leftMouseDown = mouse;
      self.zoomStartX = self.xScale.invert(mouse[0]);
      d3.event.preventDefault();
      
      // Add document listeners for mouse move and up
      d3.select(document)
        .on("mousemove.drfzoom", function() { self.handleMouseMove(); })
        .on("mouseup.drfzoom", function() { self.handleMouseUp(); });
    })
    .on("mousemove", function() {
      if (!self.leftMouseDown) {
        self.updateTooltip(d3.mouse(this));
      }
    })
    .on("mouseout", function() {
      self.hideTooltip();
    });
};

DrfChart.prototype.handleMouseMove = function() {
  if (!this.leftMouseDown) return;
  
  const currentMouse = d3.mouse(this.chartArea.select(".mouse-capture").node());
  const startX = this.leftMouseDown[0];
  const currentX = currentMouse[0];
  
  // Determine if zooming in (left to right) or out (right to left)
  const zoomingIn = currentX > startX;
  
  // Update zoom box visual indicator
  this.updateZoomBox(startX, currentX, zoomingIn);
};

DrfChart.prototype.handleMouseUp = function() {
  if (!this.leftMouseDown) return;
  
  const currentMouse = d3.mouse(this.chartArea.select(".mouse-capture").node());
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
    .on("mousemove.drfzoom", null)
    .on("mouseup.drfzoom", null);
};

DrfChart.prototype.updateZoomBox = function(startX, currentX, zoomingIn) {
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
  
  zoomBox
    .attr("x", minX)
    .attr("y", 0)
    .attr("width", width)
    .attr("height", this.chartAreaHeight);
    
  // Position text in center of zoom box
  zoomText
    .attr("x", minX + width / 2)
    .attr("y", this.chartAreaHeight / 2)
    .style("visibility", width > 50 ? "visible" : "hidden"); // Hide text if box too small
};

DrfChart.prototype.setXAxisRange = function(minEnergy, maxEnergy) {
  // Update x scale domain
  this.xScale.domain([minEnergy, maxEnergy]);
  
  // Update x axis
  this.xAxisGroup.call(this.xAxis);
  
  // Update y axis ranges based on visible data
  this.updateYAxisRanges();
  
  // Redraw lines with new scales
  this.updateLines();
  
  // Notify C++ of range change if needed
  // this.WtEmit(this.chart.id, {name: 'xRangeChanged'}, minEnergy, maxEnergy);
};

DrfChart.prototype.zoomOut = function() {
  // Expand the current domain by 50%
  const currentDomain = this.xScale.domain();
  const center = (currentDomain[0] + currentDomain[1]) / 2;
  const halfWidth = (currentDomain[1] - currentDomain[0]) / 2;
  const newHalfWidth = halfWidth * 1.5;
  
  this.setXAxisRange(
    Math.max(0, center - newHalfWidth),
    center + newHalfWidth
  );
};

DrfChart.prototype.updateYAxisRanges = function() {
  if (!this.data || this.data.length === 0) return;
  
  const xDomain = this.xScale.domain();
  const visibleData = this.data.filter(d => 
    d.energy >= xDomain[0] && d.energy <= xDomain[1]
  );
  
  if (visibleData.length === 0) return;
  
  // Update efficiency scale
  const efficiencyExtent = d3.extent(visibleData, d => d.efficiency);
  if (efficiencyExtent[0] !== undefined && efficiencyExtent[1] !== undefined) {
    const efficiencyRange = efficiencyExtent[1] - efficiencyExtent[0];
    const efficiencyPadding = efficiencyRange * 0.1;
    this.efficiencyScale.domain([
      Math.max(0, efficiencyExtent[0] - efficiencyPadding),
      efficiencyExtent[1] + efficiencyPadding
    ]);
    this.leftYAxisGroup.call(this.leftYAxis);
  }
  
  // Update FWHM scale if we have FWHM data
  const fwhmData = visibleData.filter(d => d.fwhm !== undefined && d.fwhm !== null);
  if (fwhmData.length > 0) {
    const fwhmExtent = d3.extent(fwhmData, d => d.fwhm);
    if (fwhmExtent[0] !== undefined && fwhmExtent[1] !== undefined) {
      const fwhmRange = fwhmExtent[1] - fwhmExtent[0];
      const fwhmPadding = fwhmRange * 0.1;
      this.fwhmScale.domain([
        Math.max(0, fwhmExtent[0] - fwhmPadding),
        fwhmExtent[1] + fwhmPadding
      ]);
      this.rightYAxisGroup.call(this.rightYAxis).style("display", null);
      this.rightYAxisLabel.style("display", null);
    }
  } else {
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
  }
};

DrfChart.prototype.updateLines = function() {
  if (!this.data || this.data.length === 0) return;
  
  // Update efficiency line
  this.efficiencyPath.datum(this.data.filter(d => d.efficiency !== null && d.efficiency !== undefined))
    .attr("d", this.efficiencyLine);
    
  // Update FWHM line if we have FWHM data
  const fwhmData = this.data.filter(d => d.fwhm !== null && d.fwhm !== undefined);
  if (fwhmData.length > 0) {
    this.fwhmPath.datum(fwhmData)
      .attr("d", this.fwhmLine)
      .style("display", null);
  } else {
    this.fwhmPath.style("display", "none");
  }
};

DrfChart.prototype.updateTooltip = function(mouse) {
  if (!this.data || this.data.length === 0) return;
  
  const energy = this.xScale.invert(mouse[0]);
  
  // Find closest data point
  const bisector = d3.bisector(d => d.energy).left;
  const index = bisector(this.data, energy);
  const d0 = this.data[index - 1];
  const d1 = this.data[index];
  
  let closestPoint = null;
  if (!d0) closestPoint = d1;
  else if (!d1) closestPoint = d0;
  else closestPoint = energy - d0.energy > d1.energy - energy ? d1 : d0;
  
  if (!closestPoint) return;
  
  // Show tooltip
  let tooltipContent = `<div>Energy: ${energy.toFixed(1)} keV</div>`;
  tooltipContent += `<div>Efficiency: ${closestPoint.efficiency.toFixed(4)}</div>`;
  if (closestPoint.fwhm !== undefined && closestPoint.fwhm !== null) {
    tooltipContent += `<div>FWHM: ${closestPoint.fwhm.toFixed(2)} keV</div>`;
  }
  
  this.tooltip
    .html(tooltipContent)
    .style("left", (d3.event.pageX + 10) + "px")
    .style("top", (d3.event.pageY - 10) + "px")
    .transition()
    .duration(200)
    .style("opacity", 0.9);
};

DrfChart.prototype.hideTooltip = function() {
  this.tooltip
    .transition()
    .duration(500)
    .style("opacity", 0);
};

DrfChart.prototype.setData = function(data) {
  this.data = data;
  
  if (!data || data.length === 0) {
    this.efficiencyPath.style("display", "none");
    this.fwhmPath.style("display", "none");
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
    return;
  }
  
  // Check if we have FWHM data to determine layout
  const fwhmData = data.filter(d => d.fwhm !== undefined && d.fwhm !== null);
  const hasFwhmData = fwhmData.length > 0;
  
  // Adjust margins based on whether we have FWHM data
  const originalRightMargin = this.options.margins.right;
  this.options.margins.right = hasFwhmData ? 50 : 10;
  
  // If margin changed, update dimensions and scales
  if (originalRightMargin !== this.options.margins.right) {
    this.updateDimensions();
    this.xScale.range([0, this.chartAreaWidth]);
    this.efficiencyScale.range([this.chartAreaHeight, 0]);
    this.fwhmScale.range([this.chartAreaHeight, 0]);
    
    // Update chart area
    this.chartArea.attr("transform", `translate(${this.options.margins.left}, ${this.options.margins.top})`);
    
    // Update chart area background
    if (this.chartAreaBg) {
      this.chartAreaBg
        .attr("width", this.chartAreaWidth)
        .attr("height", this.chartAreaHeight);
    }
    
    // Update clipping path
    this.chartArea.select("#drfchart-clip-" + this.chart.id + " rect")
      .attr("width", this.chartAreaWidth)
      .attr("height", this.chartAreaHeight);
    
    // Update mouse capture rect
    this.chartArea.select(".mouse-capture")
      .attr("width", this.chartAreaWidth)
      .attr("height", this.chartAreaHeight);
    
    // Update x-axis position
    this.xAxisGroup.attr("transform", `translate(0, ${this.chartAreaHeight})`);
    
    // Update right y-axis position
    this.rightYAxisGroup.attr("transform", `translate(${this.chartAreaWidth}, 0)`);
    
    // Update axis labels
    const svgRect = this.svg.node().getBoundingClientRect();
    this.xAxisLabel
      .attr("x", this.options.margins.left + this.chartAreaWidth / 2);
    this.leftYAxisLabel
      .attr("x", -(this.options.margins.top + this.chartAreaHeight / 2));
    this.rightYAxisLabel
      .attr("x", this.options.margins.top + this.chartAreaHeight / 2)
      .attr("y", -(svgRect.width - 15));
  }
  
  // Update scales based on full data extent
  const energyExtent = d3.extent(data, d => d.energy);
  const efficiencyExtent = d3.extent(data, d => d.efficiency);
  
  this.xScale.domain(energyExtent);
  this.efficiencyScale.domain([
    Math.max(0, efficiencyExtent[0] - (efficiencyExtent[1] - efficiencyExtent[0]) * 0.1),
    efficiencyExtent[1] + (efficiencyExtent[1] - efficiencyExtent[0]) * 0.1
  ]);
  
  // Update FWHM scale if we have FWHM data
  if (hasFwhmData) {
    const fwhmExtent = d3.extent(fwhmData, d => d.fwhm);
    this.fwhmScale.domain([
      Math.max(0, fwhmExtent[0] - (fwhmExtent[1] - fwhmExtent[0]) * 0.1),
      fwhmExtent[1] + (fwhmExtent[1] - fwhmExtent[0]) * 0.1
    ]);
  }
  
  // Update axes
  this.xAxisGroup.call(this.xAxis);
  this.leftYAxisGroup.call(this.leftYAxis);
  if (hasFwhmData) {
    this.rightYAxisGroup.call(this.rightYAxis);
  }
  
  // Update lines
  this.updateLines();
  
  // Show/hide right axis and label based on whether we have FWHM data
  if (hasFwhmData) {
    this.rightYAxisGroup.style("display", null);
    this.rightYAxisLabel.style("display", null);
  } else {
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
  }
};

DrfChart.prototype.handleResize = function() {
  this.updateDimensions();
  
  // Update scales range
  this.xScale.range([0, this.chartAreaWidth]);
  this.efficiencyScale.range([this.chartAreaHeight, 0]);
  this.fwhmScale.range([this.chartAreaHeight, 0]);
  
  // Update chart area transform
  this.chartArea.attr("transform", `translate(${this.options.margins.left}, ${this.options.margins.top})`);
  
  // Update chart area background
  if (this.chartAreaBg) {
    this.chartAreaBg
      .attr("width", this.chartAreaWidth)
      .attr("height", this.chartAreaHeight);
  }
  
  // Update clipping path
  this.chartArea.select("#drfchart-clip-" + this.chart.id + " rect")
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);
    
  // Update mouse capture rect
  this.chartArea.select(".mouse-capture")
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);
  
  // Update axes
  this.xAxisGroup
    .attr("transform", `translate(0, ${this.chartAreaHeight})`)
    .call(this.xAxis);
    
  this.rightYAxisGroup
    .attr("transform", `translate(${this.chartAreaWidth}, 0)`)
    .call(this.rightYAxis);
    
  this.leftYAxisGroup.call(this.leftYAxis);
  
  // Update axis labels
  const svgRect = this.svg.node().getBoundingClientRect();
  this.xAxisLabel
    .attr("x", this.options.margins.left + this.chartAreaWidth / 2)
    .attr("y", svgRect.height - 5);
    
  this.leftYAxisLabel
    .attr("x", -(this.options.margins.top + this.chartAreaHeight / 2));
    
  this.rightYAxisLabel
    .attr("x", this.options.margins.top + this.chartAreaHeight / 2)
    .attr("y", -(svgRect.width - 15));
  
  // Update lines
  this.updateLines();
};

// Method to set x-axis range from C++
DrfChart.prototype.setXRange = function(minEnergy, maxEnergy) {
  this.setXAxisRange(minEnergy, maxEnergy);
};

// Method to get current x-axis range
DrfChart.prototype.getXRange = function() {
  return this.xScale.domain();
};

// Method to set line colors from C++ (kept for compatibility)
// Line colors are now handled entirely via CSS variables in DrfChart.css
DrfChart.prototype.setLineColors = function(efficiencyColor, fwhmColor) {
  // No-op: Colors are handled by CSS variables now
  // This method is kept for backward compatibility to prevent JS errors
};