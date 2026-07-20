// Requires: DetectorPeakResponseJS.js

// Format a single value with `sigFigs` significant figures, choosing the
// shorter of fixed (toPrecision) and scientific (toExponential) notation,
// with ties broken in favour of scientific to avoid leading-zero clutter.
function drfChartFormatSigFigs(d, sigFigs) {
  if (d === 0) return "0";
  const absD = Math.abs(d);
  if (absD > 1.0) {
    // For values above 1, plain fixed notation is always compact enough.
    const decimalPlaces = Math.max(0, sigFigs - Math.floor(Math.log10(absD)) - 1);
    return d.toFixed(decimalPlaces);
  }
  const sci = d.toExponential(sigFigs - 1);
  const fix = d.toPrecision(sigFigs);
  return fix.length < sci.length ? fix : sci;
}

// Compute the [min, max] domain for the left (efficiency) Y-axis.
// When the efficiency is flat (min === max), use [0, value * 1.1] so the
// line is visible and the axis has a sensible scale.  Otherwise add 10%
// padding above and below, clamping the minimum to 0.
function drfChartEfficiencyDomain(extent) {
  const lo = extent[0], hi = extent[1];
  if (lo === hi) {
    // Flat line: anchor at zero, show 10% headroom above the value.
    return [0, hi > 0 ? hi * 1.1 : 1];
  }
  const pad = (hi - lo) * 0.1;
  return [Math.max(0, lo - pad), hi + pad];
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
    this.options.margins.left = 3; // px from SVG left edge to "Efficiency" title centre
  // computedLeft is derived dynamically from title + tick label widths; seed
  // it from margins.left so updateDimensions() has a valid value before the
  // first adjustLeftMargin() call.
  if (typeof this.options.margins.computedLeft !== "number")
    this.options.margins.computedLeft = 60;  // replaced by adjustLeftMargin() after first render

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
    .attr("transform", `translate(${this.options.margins.computedLeft}, ${this.options.margins.top})`);
    
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
    .tickFormat(d => drfChartFormatSigFigs(d, 2));

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
    .attr("x", this.options.margins.computedLeft + this.chartAreaWidth / 2)
    .attr("y", this.svg.node().getBoundingClientRect().height - 5)
    .text("Energy (keV)");

  this.leftYAxisLabel = this.svg.append("text")
    .attr("class", "axis-label y-axis-label left")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(-90)")
    .attr("x", -(this.options.margins.top + this.chartAreaHeight / 2))
    .attr("y", 15)   // will be recomputed by adjustLeftMargin()
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

  // Group for the per-angle response curves + uncertainty bands (drawn under
  // the flat efficiency/FWHM lines).
  this.angleGroup = this.plotGroup.append("g")
    .attr("class", "drf-angle-series");

  // Add paths for lines
  this.efficiencyPath = this.plotGroup.append("path")
    .attr("class", "efficiency-line")
    .style("fill", "none")
    .style("stroke-width", "2px");

  this.fwhmPath = this.plotGroup.append("path")
    .attr("class", "fwhm-line")
    .style("fill", "none")
    .style("stroke-width", "2px");

  // Legend for the per-angle curves (only shown when response angles are set).
  this.angleLegend = this.chartArea.append("g")
    .attr("class", "drf-angle-legend");

  // Per-angle response visualization state.
  this.responseSeries = null;     // parsed {angles, onAxisSolidAngleFraction, ...} or null
  this.effMode = "absolute";      // "absolute" (per emitted gamma) | "intrinsic"
  this.effIsLog = false;          // efficiency Y-axis currently logarithmic?
  // Colors for 0, 22.5, 45, 62.5, 90 degrees.
  this.angleColors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#d62728"];

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
  this.chartAreaWidth = Math.max(0, parentWidth - this.options.margins.computedLeft - this.options.margins.right);
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

  // Update the per-angle response curves
  this.updateAngleSeries();

  // Notify C++ of range change if needed
  // this.WtEmit(this.chart.id, {name: 'xRangeChanged'}, minEnergy, maxEnergy);
};

DrfChart.prototype.zoomOut = function() {
  // Zoom all the way out to show the full data range
  let minEnergy = 50;
  let maxEnergy = 3000;
  
  // Use detector extent if available
  if (this.detector) {
    const extent = this.detector.getEnergyExtent();
    minEnergy = extent[0];
    maxEnergy = extent[1];
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

  // Update efficiency scale (skipped in angle mode: updateAngleSeries owns the
  // logarithmic efficiency axis there).
  if (!this.responseSeries) {
    const efficiencyExtent = d3.extent(visibleData, d => d.efficiency);
    if (efficiencyExtent[0] !== undefined && efficiencyExtent[1] !== undefined) {
      this.efficiencyScale.domain( drfChartEfficiencyDomain(efficiencyExtent) );
      this.leftYAxisGroup.call(this.leftYAxis);
      this.adjustLeftMargin();
    }
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


DrfChart.prototype.setDetectorData = function(detectorData) {
  if (detectorData) {
    this.detector = new DetectorPeakResponseJS(detectorData);
  } else {
    this.detector = null;
  }
  
  // Set x range if energy extent is provided
  if (this.detector) {
    const energyExtent = this.detector.getEnergyExtent();
    this.setXAxisRange(energyExtent[0], energyExtent[1]);
  }
  
  this.updateEfficiencyLine();
  this.updateFwhmLine();
};


DrfChart.prototype.updateEfficiencyLine = function() {
  // In angle mode the flat curve is the angle-flat (far-field / infinite-plane)
  // intrinsic reference; it is only meaningful in the intrinsic view, and the
  // efficiency axis is owned by updateAngleSeries (do not reset it here).
  const angleMode = !!this.responseSeries;
  if (angleMode && this.effMode === "absolute") {
    this.efficiencyPath.style("display", "none");
    return;
  }

  if (!this.detector || !this.detector.hasEfficiency()) {
    this.efficiencyPath.style("display", "none");
    return;
  }

  // Generate efficiency data points
  let efficiencyPoints = [];

  // Get current x-axis domain
  const xDomain = this.xScale.domain();
  const minEnergy = xDomain[0];
  const maxEnergy = xDomain[1];
  const numPoints = Math.max(100, Math.min(600, Math.floor(this.chartAreaWidth / 2)));

  // Generate points using the detector class
  for (let i = 0; i < numPoints; i++) {
    const energy = minEnergy + (i / (numPoints - 1)) * (maxEnergy - minEnergy);
    const efficiency = this.detector.efficiency(energy);

    // Skip invalid efficiency values
    if (efficiency === null || !isFinite(efficiency) || efficiency < 0) {
      continue;
    }

    efficiencyPoints.push({ energy: energy, efficiency: efficiency });
  }

  if (efficiencyPoints.length === 0) {
    this.efficiencyPath.style("display", "none");
    return;
  }

  // Update efficiency scale based on the generated points (linear view only;
  // in angle mode updateAngleSeries owns the (log) efficiency domain).
  if (!angleMode) {
    const efficiencyExtent = d3.extent(efficiencyPoints, d => d.efficiency);
    if (efficiencyExtent[0] !== undefined && efficiencyExtent[1] !== undefined) {
      this.efficiencyScale.domain( drfChartEfficiencyDomain(efficiencyExtent) );
      this.leftYAxisGroup.call(this.leftYAxis);
      this.adjustLeftMargin();
    }
  }

  // Update the efficiency line (in intrinsic angle mode this is the far-field
  // reference; skip points that fall outside a log domain).
  if (angleMode)
    efficiencyPoints = efficiencyPoints.filter(d => d.efficiency > 0);
  this.efficiencyPath.datum(efficiencyPoints)
    .attr("d", this.efficiencyLine)
    .classed("drf-far-ref", angleMode)  // dashed far-field reference in angle mode
    .style("display", null);
};


// --- per-angle response visualization --------------------------------------

// Receives the JSON from DetectorPeakResponse::responseAngleSeriesJSON (or null
// to clear). Switches the efficiency Y-axis to logarithmic while angle curves
// are shown (efficiency spans decades), and restores linear otherwise.
DrfChart.prototype.setResponseSeries = function(series) {
  const hadSeries = !!this.responseSeries;
  this.responseSeries = (series && series.angles && series.angles.length) ? series : null;
  const haveSeries = !!this.responseSeries;

  if (haveSeries !== hadSeries) {
    // Swap the efficiency scale type; keep its pixel range.
    const range = this.efficiencyScale.range();
    this.efficiencyScale = (haveSeries ? d3.scale.log().clamp(true) : d3.scale.linear())
      .range(range);
    this.leftYAxis.scale(this.efficiencyScale);
    this.effIsLog = haveSeries;
    if (!haveSeries)
      this.leftYAxis.tickFormat(d => drfChartFormatSigFigs(d, 2));
  }

  this.updateEfficiencyLine();  // flat reference (intrinsic) / hidden (absolute)
  this.updateAngleSeries();
};


// Selects intrinsic (per-gamma-striking-face) vs absolute (per-emitted-gamma).
DrfChart.prototype.setEfficiencyMode = function(mode) {
  this.effMode = (mode === "intrinsic") ? "intrinsic" : "absolute";
  this.updateEfficiencyLine();
  this.updateAngleSeries();
};


// Maps a stored {energy,eff,sigma,flagged} pair to the plotted efficiency for
// the current mode (absolute as-sampled; intrinsic = absolute / on-axis solid
// angle fraction).
DrfChart.prototype.anglePlotValue = function(rawEff) {
  if (this.effMode !== "intrinsic")
    return rawEff;
  const f = this.responseSeries ? this.responseSeries.onAxisSolidAngleFraction : 0;
  return (f > 0) ? (rawEff / f) : rawEff;
};


DrfChart.prototype.updateAngleSeries = function() {
  const g = this.angleGroup;

  if (!this.responseSeries) {
    g.selectAll("*").remove();
    this.angleLegend.selectAll("*").remove();
    this.leftYAxisLabel.text("Efficiency");
    return;
  }

  const self = this;
  const xDomain = this.xScale.domain();

  // Build the plotted series (energy, val, lo, hi) within the visible x-range.
  const series = this.responseSeries.angles.map(function(a, idx) {
    const pts = [];
    for (let i = 0; i < a.pairs.length; i++) {
      const p = a.pairs[i];
      if (p.energy < xDomain[0] || p.energy > xDomain[1]) continue;
      const v = self.anglePlotValue(p.eff);
      if (!(v > 0) || !isFinite(v)) continue;
      const s = self.anglePlotValue(p.sigma);
      pts.push({ energy: p.energy, val: v,
                 lo: Math.max(v - s, v * 0.02), hi: v + s, flagged: p.flagged });
    }
    return { thetaDeg: a.thetaDeg, worstFlag: a.worstFlag,
             color: self.angleColors[idx % self.angleColors.length], pts: pts };
  }).filter(s => s.pts.length > 1);

  if (series.length === 0) {
    g.selectAll("*").remove();
    this.angleLegend.selectAll("*").remove();
    return;
  }

  // Efficiency (log) Y-domain over every plotted point (+ its upper band), and
  // the flat far-field reference when it is shown (intrinsic mode).
  let vmin = Infinity, vmax = -Infinity;
  series.forEach(s => s.pts.forEach(function(d){
    vmin = Math.min(vmin, d.val); vmax = Math.max(vmax, d.hi);
  }));
  if (this.effMode === "intrinsic" && this.detector && this.detector.hasEfficiency()) {
    for (let i = 0; i <= 20; i++) {
      const e = xDomain[0] + (i/20)*(xDomain[1]-xDomain[0]);
      const eff = this.detector.efficiency(e);
      if (eff > 0 && isFinite(eff)) { vmin = Math.min(vmin, eff); vmax = Math.max(vmax, eff); }
    }
  }
  if (!isFinite(vmin) || !isFinite(vmax) || vmax <= 0) {
    g.selectAll("*").remove(); this.angleLegend.selectAll("*").remove(); return;
  }
  vmin = Math.max(vmin, vmax * 1e-4);  // keep the log axis to ~4 decades
  this.efficiencyScale.domain([vmin, vmax]);
  this.leftYAxis.tickFormat(function(d){
    const l = Math.log10(d);
    return (Math.abs(l - Math.round(l)) < 1e-6) ? d.toExponential(0) : "";
  });
  this.leftYAxisGroup.call(this.leftYAxis);
  this.leftYAxisLabel.text(this.effMode === "intrinsic" ? "Intrinsic efficiency" : "Absolute efficiency");
  this.adjustLeftMargin();

  const line = d3.svg.line()
    .x(d => self.xScale(d.energy))
    .y(d => self.efficiencyScale(d.val));
  const band = d3.svg.area()
    .x(d => self.xScale(d.energy))
    .y0(d => self.efficiencyScale(Math.max(d.lo, self.efficiencyScale.domain()[0])))
    .y1(d => self.efficiencyScale(d.hi));

  // Uncertainty bands (accuracy cue), drawn under the lines.
  const bandSel = g.selectAll("path.drf-angle-band").data(series);
  bandSel.enter().append("path").attr("class", "drf-angle-band");
  bandSel.exit().remove();
  g.selectAll("path.drf-angle-band")
    .attr("d", d => band(d.pts))
    .style("fill", d => d.color)
    .style("stroke", "none")
    .style("opacity", 0.12);

  // Angle lines.
  const lineSel = g.selectAll("path.drf-angle-line").data(series);
  lineSel.enter().append("path").attr("class", "drf-angle-line")
    .style("fill", "none").style("stroke-width", "1.5px");
  lineSel.exit().remove();
  g.selectAll("path.drf-angle-line")
    .attr("d", d => line(d.pts))
    .style("stroke", d => d.color);

  this.updateEfficiencyLine();  // reflow the flat reference onto the new scale
  this.drawAngleLegend(series);
};


DrfChart.prototype.drawAngleLegend = function(series) {
  const lg = this.angleLegend;
  lg.selectAll("*").remove();

  const rows = series.map(s => ({
    color: s.color,
    text: s.thetaDeg + "°" + ((s.worstFlag && s.worstFlag !== "ok") ? " (" + s.worstFlag + ")" : "")
  }));
  rows.push({ color: null, text: (this.effMode === "intrinsic" ? "intrinsic eff." : "absolute eff. @ "
              + drfChartFormatSigFigs(this.responseSeries.distanceCm, 3) + " cm") });

  const x = this.chartAreaWidth - 150, y0 = 6, dy = 14;
  const g = lg.append("g").attr("transform", `translate(${x},${y0})`);
  rows.forEach(function(r, i){
    const row = g.append("g").attr("transform", `translate(0,${i*dy})`);
    if (r.color)
      row.append("line").attr("x1", 0).attr("x2", 16).attr("y1", 4).attr("y2", 4)
        .style("stroke", r.color).style("stroke-width", "2px");
    row.append("text").attr("x", 20).attr("y", 8)
      .attr("class", "drf-angle-legend-text").text(r.text);
  });
};


DrfChart.prototype.updateFwhmLine = function() {
  if (!this.detector || !this.detector.hasFwhm()) {
    this.fwhmPath.style("display", "none");
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
    return;
  }
  
  // Generate FWHM data points using the detector class
  let fwhmPoints = [];
  
  // Get current visible energy range (zoomed range)
  const xDomain = this.xScale.domain();
  const minEnergy = xDomain[0];
  const maxEnergy = xDomain[1];
  
  // Calculate number of points based on chart width
  const numPoints = Math.max(100, Math.min(600, Math.floor(this.chartAreaWidth / 2)));
  
  // Generate FWHM points using the detector class
  for (let i = 0; i < numPoints; i++) {
    const energy = minEnergy + (i / (numPoints - 1)) * (maxEnergy - minEnergy);
    const fwhm = this.detector.fwhm(energy);
    
    // Skip invalid FWHM values
    if (fwhm !== null && isFinite(fwhm) && fwhm >= 0 && fwhm < 9999.9) {
      fwhmPoints.push({ energy: energy, fwhm: fwhm });
    }
  }
  
  if (fwhmPoints.length === 0) {
    this.fwhmPath.style("display", "none");
    this.rightYAxisGroup.style("display", "none");
    this.rightYAxisLabel.style("display", "none");
    return;
  }
  
  // Update FWHM line
  this.fwhmPath.datum(fwhmPoints)
    .attr("d", this.fwhmLine)
    .style("display", null);
  
  // Update FWHM scale
  const fwhmExtent = d3.extent(fwhmPoints, d => d.fwhm);
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
  if (!this.detector) return;
  
  const energy = this.xScale.invert(mouse[0]);
  
  // Calculate efficiency using the detector class
  let efficiency = null;
  if (this.detector.hasEfficiency()) {
    efficiency = this.detector.efficiency(energy);
    
    // Validate efficiency value
    if (!isFinite(efficiency) || efficiency < 0) {
      efficiency = null;
    }
  }
  
  // Calculate FWHM using the detector class
  let fwhm = null;
  if (this.detector.hasFwhm()) {
    fwhm = this.detector.fwhm(energy);
    
    // Validate FWHM value
    if (!isFinite(fwhm) || fwhm < 0 || fwhm >= 9999.9) {
      fwhm = null;
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

DrfChart.prototype.handleResize = function() {
  this.updateDimensions();
  
  // Update scales range
  this.xScale.range([0, this.chartAreaWidth]);
  this.efficiencyScale.range([this.chartAreaHeight, 0]);
  this.fwhmScale.range([this.chartAreaHeight, 0]);
  
  // Update chart area transform
  this.chartArea.attr("transform", `translate(${this.options.margins.computedLeft}, ${this.options.margins.top})`);
  
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
  this.adjustLeftMargin();

  // Update axis labels
  const svgRect = this.svg.node().getBoundingClientRect();
  this.xAxisLabel
    .attr("x", this.options.margins.computedLeft + this.chartAreaWidth / 2)
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

  // Reflow the per-angle response curves onto the resized scales
  this.updateAngleSeries();
};

// Examine all rendered left Y-axis tick labels, find the minimum number of
// significant figures that makes every label string unique, reformat the
// labels in place, and return the maximum rendered text width.
DrfChart.prototype.reformatLeftYAxisLabels = function() {
  const tickNodes = this.leftYAxisGroup.selectAll("text");
  if (tickNodes.empty()) return 0;

  // In log (angle) mode the axis tickFormat already labels only the decade
  // ticks (blank elsewhere); leave the labels as-is and just measure width.
  if (this.effIsLog) {
    let maxWidth = 0;
    tickNodes.each(function() {
      const w = this.getBBox().width;
      if (w > maxWidth) maxWidth = w;
    });
    return maxWidth;
  }

  // Collect the numeric values D3 bound to each tick element.
  const values = [];
  tickNodes.each(function(d) { values.push(d); });

  // Find minimum sigFigs (1-6) that produces all-unique label strings.
  let sigFigs = 1;
  for (; sigFigs <= 6; sigFigs++) {
    const labels = values.map(v => drfChartFormatSigFigs(v, sigFigs));
    const unique = new Set(labels);
    if (unique.size === labels.length) break;
  }

  // Apply the chosen formatting and measure the widest result.
  let maxWidth = 0;
  tickNodes.each(function(d) {
    const text = drfChartFormatSigFigs(d, sigFigs);
    d3.select(this).text(text);
    const w = this.getBBox().width;
    if (w > maxWidth) maxWidth = w;
  });

  return maxWidth;
};

// Measure the widest left Y-axis tick label and adjust the left margin so
// labels never overlap the axis line or the "Efficiency" title.  Repositions
// everything that depends on the left margin.
DrfChart.prototype.adjustLeftMargin = function() {
  // Reformat labels with the minimum sig figs needed to keep them distinct,
  // and get back the widest rendered label width.
  const maxLabelWidth = this.reformatLeftYAxisLabels();

  // Measure the rendered height of the rotated "Efficiency" title; for a
  // rotated text element getBBox().height gives the font cap-height in px.
  // Fall back to 14 if the element isn't in the DOM yet.
  const titleNode = this.leftYAxisLabel.node();
  const titleHeight = titleNode ? Math.ceil(titleNode.getBBox().height) : 14;

  // Space needed from SVG left edge to the chart-area origin:
  //   margins.left – user-set gap from SVG left edge to title text
  //   titleHeight  – rendered height of the rotated "Efficiency" title
  //   maxLabelWidth – widest tick label
  //   tick length (6px) + inner gap (1px) = 1px
  const titleSlot = this.options.margins.left + titleHeight;
  const desiredLeft = titleSlot + Math.ceil(maxLabelWidth) + 7;

  if (desiredLeft === this.options.margins.computedLeft) return; // Nothing changed.

  this.options.margins.computedLeft = desiredLeft;

  // Recalculate chart area width and height.
  this.updateDimensions();

  // Move the chart area group.
  this.chartArea.attr("transform",
    `translate(${this.options.margins.computedLeft}, ${this.options.margins.top})`);

  // Update scales and dependent elements.
  this.xScale.range([0, this.chartAreaWidth]);
  this.efficiencyScale.range([this.chartAreaHeight, 0]);
  this.fwhmScale.range([this.chartAreaHeight, 0]);

  this.chartAreaBg
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);

  this.chartArea.select("#drfchart-clip-" + this.chart.id + " rect")
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);

  this.chartArea.select(".mouse-capture")
    .attr("width", this.chartAreaWidth)
    .attr("height", this.chartAreaHeight);

  this.xAxisGroup
    .attr("transform", `translate(0, ${this.chartAreaHeight})`)
    .call(this.xAxis);

  this.rightYAxisGroup
    .attr("transform", `translate(${this.chartAreaWidth}, 0)`)
    .call(this.rightYAxis);

  // Place the rotated title centre at margins.left + half its font height,
  // so the text is fully within the SVG and clear of the left edge.
  const svgRect = this.svg.node().getBoundingClientRect();
  this.xAxisLabel
    .attr("x", this.options.margins.computedLeft + this.chartAreaWidth / 2);

  this.leftYAxisLabel
    .attr("x", -(this.options.margins.top + this.chartAreaHeight / 2))
    .attr("y", this.options.margins.left + Math.round(titleHeight / 2));

  this.rightYAxisLabel
    .attr("x", this.options.margins.top + this.chartAreaHeight / 2)
    .attr("y", -(svgRect.width - 15));
};

// Method to set x-axis range from C++
DrfChart.prototype.setXRange = function(minEnergy, maxEnergy) {
  this.setXAxisRange(minEnergy, maxEnergy);
};

// Method to get current x-axis range
DrfChart.prototype.getXRange = function() {
  return this.xScale.domain();
};
