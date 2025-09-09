// Helper function to calculate A_i for Akima interpolation
// This is equivalent to the calcA function in DetectorPeakResponse.cpp
function calcA(i, xy) {
  const n = xy.length;
  const x_i = xy[i].energy;
  const x_m1 = xy[i-1].energy;
  const x_m2 = xy[i-2].energy;
  const x_p1 = xy[i+1].energy;
  const x_p2 = xy[i+2].energy;
  const y_i = xy[i].efficiency;
  const y_m1 = xy[i-1].efficiency;
  const y_m2 = xy[i-2].efficiency;
  const y_p1 = xy[i+1].efficiency;
  const y_p2 = xy[i+2].efficiency;
  const p_i = (y_p1 - y_i) / (x_p1 - x_i);
  
  let A_i;
  if (i === 0) {
    A_i = p_i;
  } else if (i === 1) {
    const p_0 = (xy[1].efficiency - xy[0].efficiency) / (xy[1].energy - xy[0].energy);
    A_i = (p_0 + p_i) / 2.0;
  } else if (i === (n-1)) {
    A_i = p_i;
  } else if (i === (n-2)) {
    const p = (xy[n-2].efficiency - xy[n-3].efficiency) / (xy[n-2].energy - xy[n-3].energy);
    A_i = (p_i + p) / 2.0;
  } else {
    const p_1 = (y_p2 - y_p1) / (x_p2 - x_p1);
    const p_m1 = (y_i - y_m1) / (x_i - x_m1);
    const p_m2 = (y_m1 - y_m2) / (x_m1 - x_m2);
    
    const w1 = Math.abs(p_1 - p_i);
    const w2 = Math.abs(p_m1 - p_m2);
    if ((w1 + w2) === 0.0) {
      A_i = (p_m1 + p_i) / 2.0;
    } else {
      A_i = (w1 * p_m1 + w2 * p_i) / (w1 + w2);
    }
  }
  return A_i;
}

// Akima interpolation function - equivalent to DetectorPeakResponse::akimaInterpolate
function akimaInterpolate(z, xy) {
  // Find position using binary search
  let left = 0;
  let right = xy.length;
  while (left < right) {
    const mid = Math.floor((left + right) / 2);
    if (xy[mid].energy < z) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  
  const n = xy.length;
  
  if (left === 0) {
    return xy[0].efficiency;
  }
  if (left === n) {
    return xy[n-1].efficiency;
  }
  
  const i = left - 1;
  
  const x_i = xy[i].energy;
  const x_p1 = xy[i+1].energy;
  const y_i = xy[i].efficiency;
  const y_p1 = xy[i+1].efficiency;
  
  if ((n < 6) || (i < 2) || ((n-i) <= 3)) {
    // Use linear interpolation
    const d = (z - x_i) / (x_p1 - x_i);
    return y_i + d * (y_p1 - y_i);
  }
  
  const dx = z - x_i;
  const h_i = (x_p1 - x_i);
  const A_i = calcA(i, xy);
  const A_p1 = calcA(i+1, xy);
  const p_i = (y_p1 - y_i) / h_i;
  const c_i = (3.0 * p_i - 2.0 * A_i - A_p1) / h_i;
  const d_i = (A_p1 + A_i - 2.0 * p_i) / h_i / h_i;
  
  return y_i + dx * (A_i + dx * (c_i + dx * d_i));
}

// Exponential of log power series efficiency function - equivalent to DetectorPeakResponse::expOfLogPowerSeriesEfficiency
function expOfLogPowerSeriesEfficiency(energy, coefs) {
  let exparg = 0.0;
  const x = Math.log(energy);
  for (let i = 0; i < coefs.length; i++) {
    exparg += coefs[i] * Math.pow(x, i);
  }
  const answer = Math.exp(exparg);
  return (answer >= 0.0) ? answer : 0.0;
}

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
  this.m_minEnergy = 0;
  this.m_maxEnergy = 3000;
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
  // Update energy range
  this.m_minEnergy = minEnergy;
  this.m_maxEnergy = maxEnergy;
  
  // Update x scale domain
  this.xScale.domain([minEnergy, maxEnergy]);
  
  // Update x axis
  this.xAxisGroup.call(this.xAxis);
  
  // Update y axis ranges based on visible data
  this.updateYAxisRanges();
  
  // Update efficiency line
  this.updateEfficiencyLine();
  
  // Update FWHM line if available
  this.updateFwhmLine();
  
  // Notify C++ of range change if needed
  // this.WtEmit(this.chart.id, {name: 'xRangeChanged'}, minEnergy, maxEnergy);
};

DrfChart.prototype.zoomOut = function() {
  // Zoom all the way out to show the full data range
  let minEnergy = 50;
  let maxEnergy = 3000;
  
  // Get energy range from efficiency data if available
  if (this.efficiencyData && this.efficiencyData.energyExtent) {
    minEnergy = this.efficiencyData.energyExtent[0];
    maxEnergy = this.efficiencyData.energyExtent[1];
  } else if (this.m_minEnergy !== undefined && this.m_maxEnergy !== undefined) {
    minEnergy = this.m_minEnergy;
    maxEnergy = this.m_maxEnergy;
  }
  
  this.setXAxisRange(minEnergy, maxEnergy);
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
  if (this.fwhmData && this.fwhmData.length > 0) {
    const visibleFwhmData = this.fwhmData.filter(d => 
      d.energy >= xDomain[0] && d.energy <= xDomain[1]
    );
    if (visibleFwhmData.length > 0) {
      const fwhmExtent = d3.extent(visibleFwhmData, d => d.fwhm);
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
  } else {
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
  }
};


DrfChart.prototype.setEfficiencyData = function(efficiencyData) {
  this.efficiencyData = efficiencyData;
  
  // Set x range if energy extent is provided
  if (efficiencyData && efficiencyData.energyExtent) {
    this.setXAxisRange(efficiencyData.energyExtent[0], efficiencyData.energyExtent[1]);
  }
  
  this.updateEfficiencyLine();
  this.updateFwhmLine();
};

DrfChart.prototype.updateEfficiencyLine = function() {
  if (!this.efficiencyData) {
    this.efficiencyPath.style("display", "none");
    return;
  }
  
  // Generate efficiency data points based on the efficiency form
  let efficiencyPoints = [];
  
  // Common setup for both efficiency forms
  const minEnergy = this.m_minEnergy || 0;
  const maxEnergy = this.m_maxEnergy || 3000;
  const numPoints = Math.max(100, Math.min(600, Math.floor(this.chartAreaWidth / 2)));
  
  if (this.efficiencyData.form === 'kEnergyEfficiencyPairs') {
    // Use Akima interpolation to expand the pairs to more points
    for (let i = 0; i < numPoints; i++) {
      const energy = minEnergy + (i / (numPoints - 1)) * (maxEnergy - minEnergy);
      const energyInUnits = energy / this.efficiencyData.energyUnits;
      const efficiency = akimaInterpolate(energyInUnits, this.efficiencyData.pairs);
      
      // Skip invalid efficiency values
      if (!isFinite(efficiency) || efficiency < 0) {
        continue;
      }
      
      efficiencyPoints.push({ energy: energy, efficiency: efficiency });
    }
  } else if (this.efficiencyData.form === 'kExpOfLogPowerSeries') {
    // Generate points using the coefficients
    for (let i = 0; i < numPoints; i++) {
      const energy = minEnergy + (i / (numPoints - 1)) * (maxEnergy - minEnergy);
      const efficiency = this.calculateEfficiency(energy);
      
      // Skip invalid efficiency values
      if (!isFinite(efficiency) || efficiency < 0) {
        continue;
      }
      
      efficiencyPoints.push({ energy: energy, efficiency: efficiency });
    }
  }
  
  if (efficiencyPoints.length === 0) {
    this.efficiencyPath.style("display", "none");
    return;
  }
  
  // Update efficiency scale based on the generated points
  const efficiencyExtent = d3.extent(efficiencyPoints, d => d.efficiency);
  if (efficiencyExtent[0] !== undefined && efficiencyExtent[1] !== undefined) {
    const efficiencyRange = efficiencyExtent[1] - efficiencyExtent[0];
    const efficiencyPadding = efficiencyRange * 0.1;
    this.efficiencyScale.domain([
      Math.max(0, efficiencyExtent[0] - efficiencyPadding),
      efficiencyExtent[1] + efficiencyPadding
    ]);
    this.leftYAxisGroup.call(this.leftYAxis);
  }
  
  // Update the efficiency line
  this.efficiencyPath.datum(efficiencyPoints)
    .attr("d", this.efficiencyLine)
    .style("display", null);
};

DrfChart.prototype.setFwhmData = function(fwhmData) {
  this.fwhmData = fwhmData;
  this.updateFwhmLine();
};

DrfChart.prototype.updateFwhmLine = function() {
  if (!this.fwhmData || this.fwhmData.length === 0) {
    this.fwhmPath.style("display", "none");
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
    return;
  }
  
  // Filter out invalid FWHM values
  const validFwhmData = this.fwhmData.filter(d => 
    isFinite(d.fwhm) && d.fwhm >= 0 && isFinite(d.energy)
  );
  
  if (validFwhmData.length === 0) {
    this.fwhmPath.style("display", "none");
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
    return;
  }
  
  // Update FWHM line
  this.fwhmPath.datum(validFwhmData)
    .attr("d", this.fwhmLine)
    .style("display", null);
  
  // Update FWHM scale
  const fwhmExtent = d3.extent(validFwhmData, d => d.fwhm);
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
};

DrfChart.prototype.updateTooltip = function(mouse) {
  // Check if we have any data to show (efficiency or FWHM)
  if (!this.efficiencyData && (!this.fwhmData || this.fwhmData.length === 0)) return;
  
  const energy = this.xScale.invert(mouse[0]);
  
  // Calculate efficiency using the appropriate method
  let efficiency = null;
  if (this.efficiencyData) {
    if (this.efficiencyData.form === 'kEnergyEfficiencyPairs' || this.efficiencyData.form === 'kExpOfLogPowerSeries') {
      efficiency = this.calculateEfficiency(energy);
    } else if (this.efficiencyData.form === 'kFunctialEfficienyForm' && this.efficiencyData.pairs) {
      // For functional form, use the pre-generated pairs
      const bisector = d3.bisector(d => d.energy).left;
      const index = bisector(this.efficiencyData.pairs, energy);
      const d0 = this.efficiencyData.pairs[index - 1];
      const d1 = this.efficiencyData.pairs[index];
      
      let closestPoint = null;
      if (!d0) closestPoint = d1;
      else if (!d1) closestPoint = d0;
      else closestPoint = energy - d0.energy > d1.energy - energy ? d1 : d0;
      
      if (closestPoint) {
        efficiency = closestPoint.efficiency;
      }
    }
  }
  
  // Find FWHM data if available
  let fwhm = null;
  if (this.fwhmData && this.fwhmData.length > 0) {
    const bisector = d3.bisector(d => d.energy).left;
    const index = bisector(this.fwhmData, energy);
    const d0 = this.fwhmData[index - 1];
    const d1 = this.fwhmData[index];
    
    let closestPoint = null;
    if (!d0) closestPoint = d1;
    else if (!d1) closestPoint = d0;
    else closestPoint = energy - d0.energy > d1.energy - energy ? d1 : d0;
    
    if (closestPoint) {
      fwhm = closestPoint.fwhm;
    }
  }
  
  // Only show tooltip if we have at least one valid value
  if (efficiency === null && fwhm === null) return;
  
  // Show tooltip
  let tooltipContent = `<div>Energy: ${energy.toFixed(1)} keV</div>`;
  if (efficiency !== null) {
    const efficiencyStr = efficiency < 0.01 ? efficiency.toExponential(3) : efficiency.toFixed(4);
    tooltipContent += `<div>Efficiency: ${efficiencyStr}</div>`;
  }
  if (fwhm !== null) {
    tooltipContent += `<div>FWHM: ${fwhm.toFixed(2)} keV</div>`;
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

// Calculate efficiency based on efficiency form
DrfChart.prototype.calculateEfficiency = function(energy) {
  if (!this.efficiencyData) {
    return 0.0;
  }
  
  const energyInUnits = energy / this.efficiencyData.energyUnits;
  
  switch (this.efficiencyData.form) {
    case 'kEnergyEfficiencyPairs':
      return akimaInterpolate(energyInUnits, this.efficiencyData.pairs);
      
    case 'kExpOfLogPowerSeries':
      return expOfLogPowerSeriesEfficiency(energyInUnits, this.efficiencyData.coefficients);
      
    case 'kFunctialEfficienyForm':
      // For functional form, we still use the pre-calculated data points
      // since we can't evaluate the function directly in JavaScript
      return this.getEfficiencyFromData(energy);
      
    default:
      return 0.0;
  }
};

// Get efficiency from pre-calculated data points (for functional form or fallback)
DrfChart.prototype.getEfficiencyFromData = function(energy) {
  if (!this.data || this.data.length === 0) {
    return 0.0;
  }
  
  // Find closest data point
  const bisector = d3.bisector(d => d.energy).left;
  const index = bisector(this.data, energy);
  const d0 = this.data[index - 1];
  const d1 = this.data[index];
  
  let closestPoint = null;
  if (!d0) closestPoint = d1;
  else if (!d1) closestPoint = d0;
  else closestPoint = energy - d0.energy > d1.energy - energy ? d1 : d0;
  
  return closestPoint ? closestPoint.efficiency : 0.0;
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
  
  // Update efficiency line
  this.updateEfficiencyLine();
  
  // Update FWHM line
  this.updateFwhmLine();
};

// Method to set x-axis range from C++
DrfChart.prototype.setXRange = function(minEnergy, maxEnergy) {
  this.setXAxisRange(minEnergy, maxEnergy);
};

// Method to get current x-axis range
DrfChart.prototype.getXRange = function() {
  return this.xScale.domain();
};
