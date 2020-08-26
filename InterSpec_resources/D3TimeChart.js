/* InterSpec: an application to analyze spectral gamma radiation data.

Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.
For questions contact William Johnson via email at wcjohns@sandia.gov, or
alternative emails of interspec@sandia.gov.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/**
 * Constructor for SampleData Objects
 * @param realTime: value of a time measurement "endpoint" (x)
 * @param gammaCPS: gamma counts per second (y)
 * @param neutronCPS: neutron counts per second (y)
 */
DataPoint = function (time, gammaCPS, neutronCPS) {
  this.time = time; // realtime or livetime
  this.gammaCPS = gammaCPS;
  this.neutronCPS = neutronCPS;
};

DataPoint.prototype.setTime = function (time) {
  this.time = time;
};

DataPoint.prototype.setGammaCPS = function (gammaCPS) {
  this.gammaCPS = gammaCPS;
};

DataPoint.prototype.setNeutronCPS = function (neutronCPS) {
  this.neutronCPS = neutronCPS;
};

DetectorMetaData = function (gammaColor, neutronColor) {
  this.isGammaDetector = gammaColor ? true : false;
  this.isNeutronDetector = neutronColor ? true : false;
  this.gammaColor = gammaColor;
  this.neutronColor = neutronColor;
};

DetectorMetaData.prototype.setGammaColor = function (gammaColor) {
  this.gammaColor = gammaColor;
  this.isGammaDetector = true;
};

DetectorMetaData.prototype.setNeutronColor = function (neutronColor) {
  this.neutronColor = neutronColor;
  this.isNeutronDetector = true;
};

/**
 * Constructor for rudimentary 1-dimensional brush along x direction
 */
BrushX = function (scale, domain) {
  this.start = undefined;
  this.end = undefined;
  this.scale = scale;
  this.domain = domain;
};

/**
 * Sets the x-scale associated with this brush
 * @param {*} scale: D3 scale object
 */
BrushX.prototype.setScale = function (scale) {
  this.scale = scale;
  return this;
};

/**
 * Returns the x-scale associated with this brush
 */
BrushX.prototype.getScale = function () {
  return this.scale;
};

/**
 * Returns true if and only if the brush extent is empty.
 * When a brush is created, it is initially empty; the brush may also become empty with a single click on the background without moving, or if the extent is cleared.
 * A brush is considered empty if it has zero-width or zero-height. When the brush is empty, its extent is not strictly defined.
 */
BrushX.prototype.empty = function () {
  return !this.start || !this.end || this.start === this.end;
};

/**
 * Clears the extent, making the brush extent empty.
 */
BrushX.prototype.clear = function () {
  this.start = undefined;
  this.end = undefined;
};

/**
 * Sets the scaled start value of the
 * @param {Number} startCoord : d3 mouse x coordinate.
 */
BrushX.prototype.setStart = function (startCoord) {
  if (!this.scale) {
    console.log("Error: brush scale has not been set!");
  }
  var scaledCoord = this.scale.invert(startCoord);
  var domain = this.scale.domain();
  if (scaledCoord < domain[0]) {
    scaledCoord = domain[0];
  }
  if (scaledCoord > domain[1]) {
    scaledCoord = domain[1];
  }
  this.start = scaledCoord;
  return this;
};

BrushX.prototype.getStart = function () {
  return this.start;
};

BrushX.prototype.setEnd = function (endCoord) {
  if (!this.scale) {
    console.log("Error: brush scale has not been set!");
  }

  var scaledCoord = this.scale.invert(endCoord);
  var domain = this.scale.domain();

  if (scaledCoord < domain[0]) {
    scaledCoord = domain[0];
  }
  if (scaledCoord > domain[1]) {
    scaledCoord = domain[1];
  }
  this.end = scaledCoord;
  return this;
};

BrushX.prototype.getEnd = function () {
  return this.end;
};

BrushX.prototype.extent = function () {
  if (this.empty()) {
    return null;
  }

  return [this.start, this.end];
};

BrushX.prototype.getCenter = function () {
  if (this.empty()) {
    return null;
  }
  return (this.start + this.end) / 2;
};
/**
 * D3TimeChart object constructor.
 */
D3TimeChart = function (elem, options) {
  // this is called when the widget is loaded.
  var self = this;
  this.chart = typeof elem === "string" ? document.getElementById(elem) : elem;

  this.options = options || {};

  if (typeof this.options.xtitle !== "string") this.options.xtitle = "Time (s)";
  if (typeof this.options.y1title !== "string")
    this.options.y1title = "Gamma CPS";
  if (typeof this.options.y2title !== "string") this.options.y2title = "N CPS";
  if (typeof this.options.compactXAxis !== "boolean")
    this.options.compactXAxis = false;
  if (typeof this.options.gridx !== "boolean") this.options.gridx = false;
  if (typeof this.options.gridy !== "boolean") this.options.gridy = false;
  if (typeof this.options.chartLineWidth !== "number")
    this.options.chartLineWidth = 1;

  this.data = undefined; // formatted data
  this.selection = undefined;
  this.height = undefined;
  this.width = undefined;
  this.leadTime = 2.5;

  this.margin = {
    top: 50,
    right: 60,
    bottom: 50,
    left: 60,
  };

  this.sourceMap = [
    "IntrinsicActivity",
    "Calibration",
    "Background",
    "Foreground",
    "Unknown",
  ];

  this.svg = d3.select(this.chart).append("svg");
  this.linesG = this.svg.append("g").attr("class", "lines");
  this.axisBottomG = this.svg.append("g").attr("class", "axis");
  this.axisLeftG = this.svg.append("g").attr("class", "axis");
  this.axisRightG = this.svg.append("g").attr("class", "axis");

  this.rectG = this.svg.append("g").attr("class", "interaction_area");
  this.highlightRect = this.rectG.append("rect").attr("class", "selection");
  this.rect = this.rectG.append("rect"); //rectangle spanning the interactable area

  this.hoverToolTip = d3
    .select(this.chart)
    .append("div")
    .attr("id", "hover_info")
    .attr("class", "tooltip")
    .style("left", `${this.margin.left + 20}px`)
    .style("top", `${this.margin.top}px`);
}; //

/** Function to help emit callbacks back to C++
 
 @param elem: The element can be a DOM element, or the object ID of a WObject
 @param event: must be an object which indicates also the JavaScript event and event target
 @param args: array of args to pass into the Wt function
 */
D3TimeChart.prototype.WtEmit = function (elem, event) {
  if (!window.Wt) {
    console.warn(
      'Wt not found! Canceling "' + event.name + '" emit function...'
    );
    return;
  }

  //console.log( 'Emitting Wt event "' + ((event && event.name) ? event.name : 'null') + '", with ' + SpectrumChartD3.prototype.WtEmit.length + " arguments");

  // To support ES5 syntax in IE11, we replace spread operator with this
  var args = Array.prototype.slice.call(
    arguments,
    SpectrumChartD3.prototype.WtEmit.length
  );

  // Emit the function to Wt
  // Wt.emit( elem, event, ...args);  // ES6 syntax
  Wt.emit.apply(Wt, [elem, event].concat(args));
};

/**
 * Sets data members of the D3TimeChart object. Is called every time data is set in C++.
 * @param {Object} data
 */
D3TimeChart.prototype.setData = function (data) {
  //See the c++ function D3TimeChart::setData()
  console.log(data);
  if (!this.isValidData(data)) {
    console.log(
      "Runtime error in D3TimeChart.setData: Structure of data is not valid.\nDoing nothing..."
      );
    } else {
    this.rawData = data;
    var compressedData = this.compress(data, 30)
    console.log("here");
    console.log(compressedData);

    var formattedData = this.formatData(data);
    console.log(formattedData);

    this.data = [formattedData];
    console.log(this.data)


    // clear selection if there is any
    this.selection = null;


    // if height and width are set, may render directly.
    if (this.height && this.width) {
      this.render();
    }
  }
};

D3TimeChart.prototype.handleResize = function () {
  // This function is called when the Wt layout manager resizes the parent <div> element
  // Need to redraw everything
  // ...
  // Make sure to update the C++ code of the changed plotting size.

  console.log("Resized!");
  this.height = this.chart.clientHeight;
  this.width = this.chart.clientWidth;

  this.render();
};

/**
 * Renders the D3TimeChart.
 * @param {*} options : Optional argument to specify render options
 */
D3TimeChart.prototype.render = function (options) {
  if (!this.data) {
    console.log(
      "Runtime error in D3TimeChart.render: D3TimeChart data is not set.\nDoing nothing..."
    );
  } else if (!this.height || !this.width) {
    console.log(
      "Runtime error in D3TimeChart.render: dimensions of D3TimeChart div element are not set.\nDoing nothing..."
    );
    return;
  } else {
    console.log("Rendering...");

    var plotWidth = this.width - this.margin.left - this.margin.right;
    var plotHeight = this.height - this.margin.top - this.margin.bottom;

    var nPoints = this.data[0].sampleNumbers.length
    // check chart pixels vs data points. Compress if needed, and add to the data array.
    // where each data[i] is data compressed at level 2^i.
    if (plotWidth < nPoints) {
      var i = 1;
      for (var i = 1; i <= Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth))); i++) {
        // only compress if data doesn't already exist before
        if (this.data[i] == null) {
          console.log("Compressing!")
          this.data[i] = this.formatData(this.compress(this.rawData, Math.pow(2, i)))
        }
      }
    }
    console.log(this.data)
    this.compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));
    console.log(this.compressionIndex);

    // set dimensions of svg element and plot
    this.svg.attr("width", this.width);
    this.svg.attr("height", this.height);

    this.rect
      .attr("width", plotWidth)
      .attr("height", plotHeight)
      .attr("x", this.margin.left)
      .attr("y", this.margin.top)
      .attr("fill-opacity", 0);

    this.rect.on("dblclick", () => {
      this.handleDoubleClick();
    });

    // add a clipPath: everything outside of this area will not be drawn
    var clip = this.svg.select("#clip_th");

    if (!clip.empty()) {
      // update if already exists
      clip
        .select("rect")
        .attr("width", plotWidth)
        .attr("height", plotHeight)
        .attr("x", this.margin.left)
        .attr("y", this.margin.top - 2); // to account for stroke-width of path
    } else {
      this.svg
        .append("defs")
        .append("svg:clipPath")
        .attr("id", "clip_th")
        .append("svg:rect")
        .attr("width", plotWidth)
        .attr("height", plotHeight)
        .attr("x", this.margin.left)
        .attr("y", this.margin.top - 2);
    }

    this.linesG.attr("clip-path", "url(#clip_th)");

    // get scales
    var domains = this.selection
    ? {
      x: this.selection.domain,
      yGamma: this.data[this.selection.compressionIndex].domains.yGamma,
      yNeutron: this.data[this.selection.compressionIndex].domains.yNeutron,
    }
    : this.data[this.compressionIndex].domains;
    
    var scales = this.getScales(domains);

    var compressionIndex = this.selection ? this.selection.compressionIndex : this.compressionIndex;
    
    // add brush-highlight zooming
    var brush = new BrushX();

    brush.setScale(scales.xScale);

    var drag = d3.behavior
      .drag()
      .on("dragstart", () => {
        var coords = d3.mouse(this.rect.node());
        brush.setStart(coords[0]);
        this.mouseDownHighlight(coords[0]);
      })
      .on("drag", () => {
        brush.setEnd(d3.mouse(this.rect.node())[0]);
        var width =
          brush.getScale()(brush.getEnd()) - brush.getScale()(brush.getStart());
        this.mouseMoveHighlight(width);
      })
      .on("dragend", () => {
        this.handleBrush(brush);
        brush.clear();
        this.mouseUpHighlight();
      });

    this.rect.call(drag);

    this.updateChart(scales, compressionIndex, options);
  }
};

/**
 * Draws or updates plots and axes
 * @param {*} scales : Object of scales with properties xScale, yScaleGamma, yScaleNeutron
 * @param {*} options : Optional argument to specify render options
 */
D3TimeChart.prototype.updateChart = function (scales, compressionIndex, options) {
  var transitions = options && options.transitions ? true : false;

  var { xScale, yScaleGamma, yScaleNeutron } = scales;

  var HAS_GAMMA = true;
  var HAS_NEUTRON = false;

  // add/update hover interaction
  this.rect
    .on("mouseover", () => {
      this.showToolTip();
    })
    .on("mousemove", () => {
      var x = xScale.invert(d3.mouse(this.rect.node())[0]);

      var idx = this.findDataIndex(x, compressionIndex);
      var startTimeStamp = this.data[compressionIndex].startTimeStamps
        ? this.data[compressionIndex].startTimeStamps[idx]
        : null;
      var sourceType = this.data[compressionIndex].sourceTypes
        ? this.data[compressionIndex].sourceTypes[idx]
        : null;
      var data = [];
      for (var detName in this.data[compressionIndex].detectors) {
        var y = this.data[compressionIndex].detectors[detName].counts[idx * 2];

        data.push({
          detName: detName,
          gammaCPS: y.gammaCPS,
          neutronCPS: y.neutronCPS,
          startTimeStamp: startTimeStamp,
          sourceType: sourceType,
        });
      }
      var optargs = { sourceType: sourceType, startTimeStamp: startTimeStamp };
      this.updateToolTip(x, data, optargs);
    })
    .on("mouseout", () => {
      this.hideToolTip();
    });

  // plot data
  for (var detName in this.data[compressionIndex].detectors) {
    var { counts, meta } = this.data[compressionIndex].detectors[detName];

    var lineGamma = d3.svg
      .line()
      .x(function (d) {
        return xScale(d.time);
      })
      .y(function (d) {
        return yScaleGamma(d.gammaCPS);
      });

    var pathGamma = this.linesG.select(".det_" + detName + "_g");

    // if already drawn, just update
    if (!pathGamma.empty()) {
      if (transitions) {
        pathGamma.datum(counts).transition().duration(500).attr("d", lineGamma);
      } else {
        pathGamma.datum(counts).attr("d", lineGamma);
      }
    } else {
      this.linesG
        .append("path")
        .attr("class", "line det_" + detName + "_g")
        .datum(counts)
        .style("stroke", meta.gammaColor)
        .style("fill", "none")
        .attr("d", lineGamma);
    } // if (!pathGamma.empty())

    if (meta.isNeutronDetector) {
      HAS_NEUTRON = true;

      var lineNeutron = d3.svg
        .line()
        .x(function (d) {
          return xScale(d.time);
        })
        .y(function (d) {
          return yScaleNeutron(d.neutronCPS);
        });

      var pathNeutron = this.linesG.select(".det_" + detName + "_n");

      // if already drawn, just update
      if (!pathNeutron.empty()) {
        if (transitions) {
          pathNeutron
            .datum(counts)
            .transition()
            .duration(500)
            .attr("d", lineNeutron);
        } else {
          pathNeutron.datum(counts).attr("d", lineNeutron);
        }
      } else {
        this.linesG
          .append("path")
          .attr("class", "line det_" + detName + "_n")
          .attr("id", "det_" + detName + "_n")
          .datum(counts)
          .style("stroke", meta.neutronColor)
          .style("fill", "none")
          .attr("d", lineNeutron);
      } // if (!pathNeutron.empty())
    }
  }

  // set different tick counts for different viewport breakpoints
  var nSamples = this.data[compressionIndex].sampleNumbers.length;
  // if (this.width > 850) {
  //   nTicksX = nSamples;
  // } else if (this.width > 520) {
  //   nTicksX = Math.floor(nSamples / 2);
  // } else if (this.width > 280) {
  //   nTicksX = Math.floor(nSamples / 4);
  // } else {
  //   nTicksX = Math.floor(nSamples / 8);
  // }

  // plot axes and labels
  var xAxis = d3.svg.axis().scale(xScale);

  // update or create axis
  if (transitions) {
    this.axisBottomG
      .transition()
      .duration(500)
      .attr(
        "transform",
        "translate(0," + (this.height - this.margin.bottom) + ")"
      )
      .call(xAxis);
  } else {
    this.axisBottomG
      .attr(
        "transform",
        "translate(0," + (this.height - this.margin.bottom) + ")"
      )
      .call(xAxis);
  }

  var axisLabelX = this.svg.select("#th_label_x");

  // if already drawn, just update
  if (!axisLabelX.empty()) {
    // clear existing transforms
    axisLabelX.attr("transform", "none");

    // reposition
    axisLabelX.attr(
      "transform",
      `translate(${this.width / 2}, ${
        this.height -
        this.margin.bottom +
        this.axisBottomG.node().getBBox().height +
        15
      })`
    );
  } else {
    this.svg
      .append("text")
      .attr("class", "axis_label")
      .attr("id", "th_label_x")
      .attr(
        "transform",
        `translate(${this.width / 2}, ${
          this.height -
          this.margin.bottom +
          this.axisBottomG.node().getBBox().height +
          15
        })`
      )
      .style("text-anchor", "middle")
      .text(this.options.xtitle);
  } // if (!axisLabelX.empty())

  if (HAS_GAMMA) {
    var yAxisLeft = d3.svg
      .axis()
      .scale(yScaleGamma)
      .ticks(3)
      .orient("left")
      .tickFormat(d3.format(".1g"));

    // update or create axis
    this.axisLeftG
      .attr("transform", "translate(" + this.margin.left + ",0)")
      .call(yAxisLeft);

    var axisLabelY1 = this.svg.select("#th_label_y1");

    // if already drawn, just update
    if (!axisLabelY1.empty()) {
      // clear existing transforms
      axisLabelY1.attr("transform", "none");

      // reposition
      axisLabelY1.attr(
        "transform",
        `translate(${
          this.margin.left - this.axisLeftG.node().getBBox().width - 5
        }, ${this.height / 2}) rotate(-90)`
      );
    } else {
      this.svg
        .append("text")
        .attr("class", "axis_label")
        .attr("id", "th_label_y1")
        .attr(
          "transform",
          `translate(${
            this.margin.left - this.axisLeftG.node().getBBox().width - 5
          }, ${this.height / 2}) rotate(-90)`
        )
        .style("text-anchor", "middle")
        .text(this.options.y1title)
        .attr("font-size", "0.9em");
    } // if (!axisLabelY1.empty())
  } // if (HAS_GAMMA)

  if (HAS_NEUTRON) {
    var yAxisRight = d3.svg
      .axis()
      .scale(yScaleNeutron)
      .ticks(3)
      .orient("right");

    // create or update axis
    this.axisRightG
      .attr("transform", "translate(" + (this.width - this.margin.left) + ",0)")
      .call(yAxisRight);

    var axisLabelY2 = this.svg.select("#th_label_y2");

    // if already drawn, just update
    if (!axisLabelY2.empty()) {
      // clear existing transforms
      axisLabelY2.attr("transform", "none");

      // reposition
      axisLabelY2.attr(
        "transform",
        `translate(${
          this.width -
          this.margin.left +
          this.axisRightG.node().getBBox().width +
          10
        }, ${this.height / 2}) rotate(90)`
      );
    } else {
      this.svg
        .append("text")
        .attr("class", "axis_label")
        .attr("id", "th_label_y2")
        .attr(
          "transform",
          `translate(${
            this.width -
            this.margin.left +
            this.axisRightG.node().getBBox().width +
            10
          }, ${this.height / 2}) rotate(90)`
        )
        .style("text-anchor", "middle")
        .text(this.options.y2title);
    } // if (!axisLabelY2.empty())
  } // if (HAS_NEUTRON)
};

// HELPERS, CALLBACKS, AND EVENT HANDLERS //

/**
 * Basic data validation. Data should be object of form:
 * {
 *    realTimes: Numbers[N],
 *    sampleNumber: Numbers[N],
 *    gammaCounts: [{
 *                    detName: String,
 *                    color: String
 *                    counts: Numbers[N]
 *                  }
 *                  ...
 *                  ],
 *
 * <-- optional -->
 *    neutronCounts: [{
 *                    detName: String,
 *                    color: String
 *                    counts: Numbers[N]
 *                  },
 *                  ...
 *                  ],
 *
 *    occupancies: [{...}, {...}]
 * }
 *
 * @param {Object} data: data Object passed from Wt
 * @returns a boolean
 */
D3TimeChart.prototype.isValidData = function (data) {
  if (!data) {
    return false;
  }

  // Check properties
  if (
    !data.hasOwnProperty("sampleNumbers") ||
    !data.hasOwnProperty("realTimes") ||
    !data.hasOwnProperty("gammaCounts")
  ) {
    return false;
  }

  // Check data types
  if (
    !Array.isArray(data.sampleNumbers) ||
    !Array.isArray(data.realTimes) ||
    !Array.isArray(data.gammaCounts)
  ) {
    return false;
  }

  var nSamples = data.sampleNumbers.length;

  // Check matching lengths
  if (data.realTimes.length !== nSamples) {
    return false;
  }

  // check data of each gamma detector
  data.gammaCounts.forEach(function (detector) {
    // check  fields
    if (
      !detector.hasOwnProperty("color") ||
      !detector.hasOwnProperty("counts") ||
      detector.hasOwnProperty("detName")
    ) {
      return false;
    }

    // check length of counts array
    if (
      !detector.counts ||
      !Array.isArray(detector.counts) ||
      detector.counts.length !== nSamples
    ) {
      return false;
    }
  });

  //optional properties
  if (data.hasOwnProperty("neutronCounts")) {
    if (!Array.isArray(data.neutronCounts)) {
      return false;
    }

    data.neutronCounts.forEach(function (detector) {
      // check  fields
      if (
        !detector.hasOwnProperty("color") ||
        !detector.hasOwnProperty("counts") ||
        detector.hasOwnProperty("detName")
      ) {
        return false;
      }

      // check length of counts array
      if (
        !Array.isArray(detector.counts) ||
        detector.counts.length !== nSamples
      ) {
        return false;
      }
    });
  } // if (data.hasOwnProperty("neutronCounts"))

  return true;
};

/**
 * Formats data to several JSON objects for convenient plotting
 * returns object bundling sampleNumbers, realTimeIntervals, data domain, and detectors data.
 * e.g.
 * sampleNumbers: [1, 2, 3, ...],
 * realTimeIntervals: [[a,b], [c,d], ...],
 * meanIntervalTime: 0.500...
 * gpsCoordinates: [...], (optional)
 * occupied: [false, true, true, ...]
 * sourceTypes: [2, 3, 3, ..]
 * startTimeStamps: [...]
 * domains:
 * {
 *    x: [a, b],
 *    yGamma: [c, d]
 *    yNeutron: [e, f]
 * }
 *
 * detectors:
 * {
 *    det1: {
 *              meta: {
 *                        "isGammaDetector": true,
 *                        "isNeutronDetector": true,
 *                        "gammaColor": rgb(0,0,0),
 *                        "neutronColor": rgb(0,128,0)
 *                    },
 *              counts: [{DataPoint}, {DataPoint}, {DataPoint}, ...]
 *          },
 *
 *    det2: {...},
 *    ...
 * }
 *
 * @param {*} data
 */
D3TimeChart.prototype.formatData = function (data) {
  var nSamples = data.sampleNumbers.length;
  var detectors = {};
  var realTimeIntervals = this.getRealTimeIntervals(
    data.realTimes,
    data.sourceTypes
  );
  var meanIntervalTime = this.getMeanIntervalTime(
    data.realTimes,
    data.sourceTypes
  );

  // TODO: fix to use ES5
  // create occupied array:
  var occupied;
  if (data.hasOwnProperty("occupancies")) {
    var occupied = [];
    for (var i = 0; i < nSamples; i++) {
      var sampleNumber = (typeof data.sampleNumbers[i] === 'number' && isFinite(data.sampleNumbers[i])) ? data.sampleNumbers[i] : [...data.sampleNumbers[i]][0]
      occupied[i] = this.isOccupiedSample(
        sampleNumber,
        data.occupancies
      );
    }
  }

  // create startTimeStamps array:
  var startTimeStamps;
  if (
    data.hasOwnProperty("startTimeOffset") &&
    data.hasOwnProperty("startTimes")
  ) {
    startTimeStamps = data.startTimes.map(function (startTime) {
      if (startTime == null) {
        return null;
      }
      return data.startTimeOffset + startTime;
    });
  }

  var domains = this.getDomains(data);

  // format counts data
  data.gammaCounts.forEach(function (det) {
    // if new detector, initialize a new object and initialize metadata
    if (!detectors.hasOwnProperty(det.detName)) {
      detectors[det.detName] = {};
      detectors[det.detName].meta = new DetectorMetaData(det.color, undefined);
    } else {
      //detector already exists, so set gamma metadata
      detectors[det.detName].meta.setGammaColor(det.color);
    }
    // add gamma det livetimes if exists
    detectors[det.detName].gammaLiveTimes = det.liveTimes;

    // if no data already present for this detector, create a new array to start holding data for that detector and fill in the data.
    if (!detectors[det.detName].hasOwnProperty("data")) {
      detectors[det.detName].counts = [];

      for (var i = 0; i < nSamples; i++) {
        // use livetimes for cps if available; realtimes otherwise
        var cps = det.liveTimes
          ? det.counts[i] / det.liveTimes[i]
          : det.counts[i] / data.realTimes[i];

        // push line segment start
        detectors[det.detName].counts.push(
          new DataPoint(realTimeIntervals[i][0], cps, null)
        );
        // push line segment end
        detectors[det.detName].counts.push(
          new DataPoint(realTimeIntervals[i][1], cps, null)
        );
      }
    } else {
      // data is already present for this detector, so only set gamma CPS for each data point.
      for (var i = 0; i < nSamples; i++) {
        // use livetimes for cps if available; realtimes otherwise
        var cps = det.liveTimes
          ? det.counts[i] / det.liveTimes[i]
          : det.counts[i] / data.realTimes[i];
        detectors[det.detName].counts[2 * i].setGammaCPS(cps);
        detectors[det.detName].counts[2 * i + 1].setGammaCPS(cps);
      }
    }
  }); // data.gammaCounts.forEach

  if (data.hasOwnProperty("neutronCounts")) {
    data.neutronCounts.forEach(function (det) {
      // if new detector, initialize a new object and initialize metadata
      var detName = det.detName;
      if (!detectors.hasOwnProperty(det.detName)) {
        detectors[det.detName] = {};
        detectors[det.detName].meta = new DetectorMetaData(
          undefined,
          det.color
        );
      } else {
        // detector already exists, so set neutron metadata
        detectors[detName].meta.setNeutronColor(det.color);
      }
      // add neutron livetime if exists
      detectors[det.detName].neutronLiveTimes = det.liveTimes;

      // if no data already present for this detector, create a new array to start holding data for that detector and fill in the data.
      if (!detectors[det.detName].hasOwnProperty("counts")) {
        detectors[det.detName].counts = [];

        for (var i = 0; i < nSamples; i++) {
          // use livetimes for cps if available; realtimes otherwise
          var cps = det.liveTimes
            ? det.counts[i] / det.liveTimes[i]
            : det.counts[i] / data.realTimes[i];

          //push line segment start
          detectors[det.detName].counts.push(
            new DataPoint(realTimeIntervals[i][0], null, cps)
          );
          detectors[det.detName].counts.push(
            new DataPoint(realTimeIntervals[i][1], null, cps)
          );
        }
      } else {
        // data is already present for this detector, so only set neutron CPS for each data point.
        for (var i = 0; i < nSamples; i++) {
          // use livetimes for cps if available; realtimes otherwise
          var cps = det.liveTimes
            ? det.counts[i] / det.liveTimes[i]
            : det.counts[i] / data.realTimes[i];
          detectors[det.detName].counts[2 * i].setNeutronCPS(cps);
          detectors[det.detName].counts[2 * i + 1].setNeutronCPS(cps);
        }
      } // if (!detectors[det.detName].hasOwnProperty("data"))
    }); // data.neutronCounts.forEach
  } // if (data.hasOwnProperty("neutronCounts"))

  return {
    sampleNumbers: data.sampleNumbers,
    realTimeIntervals: realTimeIntervals,
    meanIntervalTime: meanIntervalTime,
    occupied: occupied,
    sourceTypes: data.sourceTypes,
    startTimeStamps: startTimeStamps,
    gpsCoordinates: data.gpsCoordinates,

    detectors: detectors,
    domains: domains,
  };
};

/**
 * Gets data domains.
 * @param {Object} data: raw data from Wt
 */
D3TimeChart.prototype.getDomains = function (data) {
  var realTimeIntervals = this.getRealTimeIntervals(
    data.realTimes,
    data.sourceTypes
  );

  var xMin = d3.min(realTimeIntervals, function (d) {
    return d[0];
  });
  var xMax = d3.max(realTimeIntervals, function (d) {
    return d[1];
  });

  var yMaxGamma = Number.MIN_SAFE_INTEGER;
  var yMaxNeutron = Number.MIN_SAFE_INTEGER;

  // Find max over gamma detectors
  var nSamples = data.sampleNumbers.length;

  for (var i = 0; i < data.gammaCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.gammaCounts[i].liveTimes
        ? data.gammaCounts[i].counts[j] / data.gammaCounts[i].liveTimes[j]
        : data.gammaCounts[i].counts[j] / data.realTimes[j];

      if (cps > yMaxGamma) {
        yMaxGamma = cps;
      }
    }
  }

  // Find max over neutron detectors
  for (var i = 0; i < data.neutronCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.neutronCounts[i].liveTimes
        ? data.neutronCounts[i].counts[j] / data.neutronCounts[i].liveTimes[j]
        : data.neutronCounts[i].counts[j] / data.realTimes[j];

      if (cps > yMaxNeutron) {
        yMaxNeutron = cps;
      }
    }
  }

  var gammaInterval = yMaxGamma > 0 ? [0, yMaxGamma] : undefined;
  var yNeutronInterval = yMaxNeutron > 0 ? [0, yMaxNeutron] : undefined;

  return {
    x: [xMin, xMax],
    yGamma: gammaInterval,
    yNeutron: yNeutronInterval,
  };
};

/**
 * Get scaling functions based on the data domain and element dimensions.
 * this.data, this.height, and this.width must be defined
 * Argument domains must be of form:
 * {
 *    x: [a, b],
 *    yGamma: [c, d],
 *    yNeutron: [e, f],
 * }
 */
D3TimeChart.prototype.getScales = function (domains) {
  if (!this.data) {
    console.log(
      "In D3TimeChart.getScales: domain not set.\nScales undefined..."
    );
    return {
      xScale: undefined,
      yScaleGamma: undefined,
      yScaleNeutron: undefined,
    };
  }
  var xScale = domains.x
    ? d3.scale
        .linear()
        .domain(domains.x)
        .range([this.margin.left, this.width - this.margin.right])
    : undefined;
  var yScaleGamma = domains.yGamma
    ? d3.scale
        .linear()
        .domain(domains.yGamma)
        .range([this.height - this.margin.bottom, this.margin.top])
    : undefined;
  var yScaleNeutron = domains.yNeutron
    ? d3.scale
        .linear()
        .domain(domains.yNeutron)
        .range([this.height - this.margin.bottom, this.margin.top])
    : undefined;

  return {
    xScale: xScale,
    yScaleGamma: yScaleGamma,
    yScaleNeutron: yScaleNeutron,
  };
};

/**
 * Computes sequential real time intervals for each data point from raw time segments.
 * @param {number[]} realTimes: realTimes data passed from Wt
 * @returns an array of length-two arrays which represent time intervals for individual samples.
 */
D3TimeChart.prototype.getRealTimeIntervals = function (realTimes, sourceTypes) {
  var realTimeIntervals = [];

  for (var i = 0; i < realTimes.length; i++) {
    if (sourceTypes && i === 0 && sourceTypes[i] === 2) {
      // center so background is not started at 0
      // to handle long lead-in:
      var leadTime = Math.max(-realTimes[i], -this.leadTime);
      realTimeIntervals[i] = [leadTime, 0];
    } else if (i === 0) {
      // first point is not background, so don't need to center
      realTimeIntervals[i] = [0, realTimes[i]];
    } else {
      var prevEndpoint = realTimeIntervals[i - 1][1];
      realTimeIntervals[i] = [prevEndpoint, prevEndpoint + realTimes[i]];
    }
  }

  return realTimeIntervals;
};

D3TimeChart.prototype.getMeanIntervalTime = function (realTimes, sourceTypes) {
  var n = 0;
  var acc = 0;
  for (var i = 0; i < realTimes.length; i++) {
    if (sourceTypes[i] == 3) {
      acc += realTimes[i];
      n += 1;
    }
  }
  return acc / n;
};

/**
 * Brush handler to  zoom into a selected domain
 * @param {} brush : d3 brush
 */
D3TimeChart.prototype.handleBrush = function (brush) {
  if (brush && !brush.empty() && brush.extent()[0] < brush.extent()[1]) {
    // find appropriate compression level for this range
    console.log(brush.extent());
    var leftIndex = this.findDataIndex(brush.extent()[0], 0)
    var rightIndex = this.findDataIndex(brush.extent()[1], 0)
    console.log([leftIndex, rightIndex])
    var nPoints = rightIndex - leftIndex + 1;
    var plotWidth = this.width - this.margin.left - this.margin.right;

    var compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));



    // set lower limit on extent size to 2 interval lengths
    var minExtentLeft = Math.max(
      brush.getCenter() - this.data[compressionIndex].meanIntervalTime,
      this.data[compressionIndex].domains.x[0]
    );
    var maxExtentRight = Math.min(
      brush.getCenter() + this.data[compressionIndex].meanIntervalTime,
      this.data[compressionIndex].domains.x[1]
    );

    var brushWidth = brush.extent()[1] - brush.extent()[0];

    var extent =
      brushWidth > 2 * this.data[compressionIndex].meanIntervalTime
        ? brush.extent()
        : [minExtentLeft, maxExtentRight];

    // update x-domain to the new domain
    this.selection ={domain: extent, compressionIndex: compressionIndex};

    var scales = this.getScales({
      x: this.selection.domain,
      yGamma: this.data[compressionIndex].domains.yGamma,
      yNeutron: this.data[compressionIndex].domains.yNeutron,
    });

    // brush.x(scales.xScale); // update xScale
    brush.setScale(scales.xScale);
    var transitions = compressionIndex === this.compressionIndex;
    this.updateChart(scales, compressionIndex, { transitions: transitions });

    // to clear the drawn rectangle.
    // this.brushG.call(brush.clear());
  }
};

/**
 * Double click handler to redraw chart with full (default) domain
 */
D3TimeChart.prototype.handleDoubleClick = function () {
  var transitions = this.selection.compressionIndex === this.compressionIndex;
  // clear selection
  this.selection = null;

  // redraw with full domain
  this.render({ transitions: transitions});
};

/**
 * Sets initial attributes of a highlight region
 * @param {Number} mouseX : x-coordinates of pointer in pixels relative to the containing element
 */
D3TimeChart.prototype.mouseDownHighlight = function (mouseX) {
  this.highlightRect
    .attr("height", this.height - this.margin.top - this.margin.bottom)
    .attr("x", mouseX)
    .attr("y", this.margin.top)
    .attr("width", 0);
};

/**
 * Updates width of highlight region on mouse move.
 * @param {Number} width : width in pixels
 */
D3TimeChart.prototype.mouseMoveHighlight = function (width) {
  if (width > 0) {
    this.highlightRect.attr("width", width);
  } else {
    this.highlightRect.attr("width", 0);
  }
};

/**
 * Clears rectangle dimensions
 */
D3TimeChart.prototype.mouseUpHighlight = function () {
  this.highlightRect.attr("height", 0).attr("width", 0);
};

D3TimeChart.prototype.showToolTip = function () {
  this.hoverToolTip.style("visibility", "visible");
};

/**
 * Handler to update tooltip display.
 * @param {*} time : real time of measurement (x-value)
 * @param {*} data : Array of data objects, where each object at minimum has fields: {detName, gammaCPS}. Optional fields: neutronCPS
 * @param {*} data : Object of optional keyword arguments. Accepted arguments include: startTimeStamp, sourceType
 */
D3TimeChart.prototype.updateToolTip = function (time, data, optargs) {
  this.hoverToolTip.html(this.createToolTipString(time, data, optargs));
};

D3TimeChart.prototype.hideToolTip = function () {
  this.hoverToolTip.style("visibility", "hidden");
};

D3TimeChart.prototype.createToolTipString = function (time, data, optargs) {
  var s =
    optargs.startTimeStamp != null
      ? `<div>${new Date(optargs.startTimeStamp).toUTCString()}</div>`
      : "";
  s +=
    optargs.sourceType != null
      ? `<div>Source: ${this.sourceMap[optargs.sourceType]}</div>`
      : "";
  s += `<div>Time: ${time.toPrecision(4)} s</div>`;

  // for each detector, give counts
  for (var i = 0; i < data.length; i++) {
    s += `<div>G CPS: ${data[i].gammaCPS.toPrecision(6)} (${
      data[i].detName
    })</div>`;

    if (data[i].neutronCPS != null) {
      // cps of 0 is still valid to display
      s += `<div>N CPS: ${data[i].neutronCPS.toPrecision(3)} (${
        data[i].detName
      })</div>`;
    }
  }

  return s;
};

/**
 * Performs search over data to find interval that the given time belongs to
 * Returns the index of the match in the data array
 * @param {} time : realTime of measurement (x-value)
 * @param {} compressionIndex : optional compression index to search over.
 */
D3TimeChart.prototype.findDataIndex = function (time, compressionIndex) {
  var cIdx = (compressionIndex != null) ? compressionIndex : this.compressionIndex;
  var highIdx = this.data[cIdx].realTimeIntervals.length - 1;
  var lowIdx = 0;

  while (lowIdx <= highIdx) {
    var midIdx = Math.floor((highIdx + lowIdx) / 2);
    var interval = this.data[cIdx].realTimeIntervals[midIdx];
    if (time >= interval[0] && time <= interval[1]) {
      return midIdx;
    } else if (time < interval[0]) {
      highIdx = midIdx - 1;
    } else if (time > interval[1]) {
      lowIdx = midIdx + 1;
    }
  }
  return -1;
};

D3TimeChart.prototype.isOccupiedSample = function (sampleNumber, occupancies) {
  // if occupancies.status is undefined, assume it is true...
  var status = occupancies.hasOwnProperty(status) ? occupancies.status : true;
  for (var i = 0; i < occupancies.length; i++) {
    if (
      sampleNumber >= occupancies[i].startSample &&
      sampleNumber <= occupancies[i].endSample
    ) {
      return !(true ^ status);
    }
  }
  return !(false ^ status);
};

// TODO: Fix to use ES5 (Set() is from ES6)
D3TimeChart.prototype.compress = function (data, n) {
  // Create output array template
  var out = {
    gammaCounts: [],
    realTimes: [],
    sampleNumbers: [],
  };

  // init flags
  var HAS_OCCUPANCY_DATA = false;
  var HAS_SOURCE_TYPE_DATA = false;
  var HAS_START_TIME_DATA = false;

  // Add optional fields:
  if (data.hasOwnProperty("neutronCounts")) {
    out.neutronCounts = [];
  }

  if (data.hasOwnProperty("occupancies")) {
    out.occupancies = data.occupancies;
    HAS_OCCUPANCY_DATA = true;
  }

  if (data.hasOwnProperty("sourceTypes")) {
    out.sourceTypes = [];
    HAS_SOURCE_TYPE_DATA = true;
  }

  if (
    data.hasOwnProperty("startTimeOffset") &&
    data.hasOwnProperty("startTimes")
  ) {
    out.startTimeOffset = data.startTimeOffset;
    out.startTimes = [];
    HAS_START_TIME_DATA = true;
  }

  // iterate over array with window of size n
  var length = data.sampleNumbers.length;

  // create occupied array:
  var occupied;
  if (HAS_OCCUPANCY_DATA) {
    var occupied = [];
    for (var i = 0; i < length; i++) {
      occupied[i] = this.isOccupiedSample(
        data.sampleNumbers[i],
        data.occupancies
      );
    }
  }

  // current index of output array
  var outIdx = 0;
  // current index of data array
  var i = 0;

  // data accumulators
  var detectorsAccumulator = {};

  while (i < length) {
    if (HAS_START_TIME_DATA) {
      out.startTimes[outIdx] = data.startTimes[i];
    }
    if (HAS_SOURCE_TYPE_DATA) {
      out.sourceTypes[outIdx] = data.sourceTypes[i];
    }
    var j = 0;
    while (j < n && i + j < length) {
      if (
        (!HAS_OCCUPANCY_DATA || occupied[i + j] === occupied[i]) &&
        (!HAS_SOURCE_TYPE_DATA ||
          data.sourceTypes[i + j] === data.sourceTypes[i])
      ) {
        // aggregation conditions: (occupancy data doesn't exist OR occupancies match) AND (source types data doesn't exist OR sourcetypes match).
        // aggregate the i + j'th realTime into the value of realTimesAccumulator[outIdx]. This will set the appropriate place inside the output array eventually.
        if (out.realTimes[outIdx] == null) {
          out.realTimes[outIdx] = 0;
        } // if (out.realTimes[outIdx] == null)
        out.realTimes[outIdx] += data.realTimes[i + j];

        // go into gammaCounts.
        // For each "detector" of gammaCounts...
        data.gammaCounts.forEach(function (detector) {
          // if the detectorsAccumulator doesn't already have that detector's name as a property, then add it to the detectorsAccumulator
          // And initialize its fields. The fields should contain:
          // gammaColor (string), gammaCounts (Number[]). Optional: gammaLiveTimes (Number[]), neutronColor (string), neutronCounts (Number[]), neutronLiveTimes(Number[])...
          if (!detectorsAccumulator.hasOwnProperty(detector.detName)) {
            detectorsAccumulator[detector.detName] = {};
          } // if (!detectorsAccumulator.hasOwnProperty(detector.detName))

          // if detector already exist but does not have gamma fields, then add them.
          if (
            !detectorsAccumulator[detector.detName].hasOwnProperty(
              "gammaCounts"
            )
          ) {
            detectorsAccumulator[detector.detName].gammaCounts = [];
          } // if (!detectorsAccumulator[detector.detName].hasOwnProperty("gammaCounts"))

          if (
            !detectorsAccumulator[detector.detName].hasOwnProperty("gammaColor")
          ) {
            detectorsAccumulator[detector.detName].gammaColor = detector.color;
          } // if (!detectorsAccumulator[detector.detName].hasOwnProperty("gammaColor"))

          if (
            detector.hasOwnProperty("liveTimes") &&
            !detectorsAccumulator[detector.detName].hasOwnProperty(
              "gammaLiveTimes"
            )
          ) {
            detectorsAccumulator[detector.detName].gammaLiveTimes = [];
          } // if (detector.hasOwnProperty("liveTimes") && !detectorsAccumulator[detector.detName].hasOwnProperty("gammaLiveTimes"))

          // gammaCounts[outIdx] will begin accumulating. So add to it the i + jth count.
          if (
            detectorsAccumulator[detector.detName].gammaCounts[outIdx] == null
          ) {
            detectorsAccumulator[detector.detName].gammaCounts[outIdx] = 0;
          } // if (detectorsAccumulator[detector.detName].gammaCounts[outIdx] == null)
          detectorsAccumulator[detector.detName].gammaCounts[outIdx] +=
            detector.counts[i + j];

          // if have liveTimes, it will also begin accumulating, so add to it the i + jth livetime.
          if (detector.hasOwnProperty("liveTimes")) {
            if (
              detectorsAccumulator[detector.detName].gammaLiveTimes[outIdx] ==
              null
            ) {
              detectorsAccumulator[detector.detName].gammaLiveTimes[outIdx] = 0;
            } // if (detectorsAccumulator[detector.detName].gammaLiveTimes[outIdx] == null)
            detectorsAccumulator[detector.detName].gammaLiveTimes[outIdx] +=
              detector.liveTimes[i + j];
          } // (detector.hasOwnProperty("liveTimes"))
        });

        // do same with neutronCounts
        if (data.hasOwnProperty("neutronCounts")) {
          data.neutronCounts.forEach(function (detector) {
            // if the detectorsAccumulator doesn't already have that detector's name as a property, then add it to the detectorsAccumulator
            // And initialize its fields. The fields should contain:
            // gammaColor (string), gammaCounts (Number[]). Optional: gammaLiveTimes (Number[]), neutronColor (string), neutronCounts (Number[]), neutronLiveTimes(Number[])...
            if (!detectorsAccumulator.hasOwnProperty(detector.detName)) {
              detectorsAccumulator[detector.detName] = {};
            }
            // if detector already exist but does not have neutron fields, then add them.
            if (
              !detectorsAccumulator[detector.detName].hasOwnProperty(
                "neutronCounts"
              )
            ) {
              detectorsAccumulator[detector.detName].neutronCounts = [];
            }

            if (
              !detectorsAccumulator[detector.detName].hasOwnProperty(
                "neutronColor"
              )
            ) {
              detectorsAccumulator[detector.detName].neutronColor =
                detector.color;
            }

            if (
              detector.hasOwnProperty("liveTimes") &&
              !detectorsAccumulator[detector.detName].hasOwnProperty(
                "neutronLiveTimes"
              )
            ) {
              detectorsAccumulator[detector.detName].neutronLiveTimes = [];
            }

            // neutronCounts[outIdx] will begin accumulating. So add to it the i + jth count.
            if (
              detectorsAccumulator[detector.detName].neutronCounts[outIdx] ==
              null
            ) {
              detectorsAccumulator[detector.detName].neutronCounts[outIdx] = 0;
            }
            detectorsAccumulator[detector.detName].neutronCounts[outIdx] +=
              detector.counts[i + j];

            // if have liveTimes, it will also begin accumulating, so add to it the i + jth livetime.
            if (detector.hasOwnProperty("liveTimes")) {
              if (
                detectorsAccumulator[detector.detName].neutronLiveTimes[
                  outIdx
                ] == null
              ) {
                detectorsAccumulator[detector.detName].neutronLiveTimes[
                  outIdx
                ] = 0;
              }
              detectorsAccumulator[detector.detName].neutronLiveTimes[outIdx] +=
                detector.liveTimes[i + j];
            }
          });
        } // if (data.hasOwnProperty("neutronCounts"))

        // add to sampleNumbers:
        if (out.sampleNumbers[outIdx] == null) {
          out.sampleNumbers[outIdx] = new Set();
        }
        out.sampleNumbers[outIdx].add(data.sampleNumbers[i + j]);

        // increment j
        j += 1;

        // if would reach exit condition, increment i and outIdx appropriately before exiting.
        if (j >= n || i + j >= length) {
          i = i + j;
          outIdx += 1;
        }
      } else {
        // the i + j'th data point cannot be aggregated with the rest... so increment i and outIdx appropriately and break out of the while-loop
        i = i + j;
        outIdx += 1;
        break;
      } // if ((!HAS_OCCUPANCY_DATA || occupied[i + j] === occupied[i]) && (!HAS_SOURCE_TYPE_DATA || data.sourceTypes[i + j] === data.sourceTypes[i]))
    } // while (j < n && i + j < length)
  } // while (i < length)

  // process for output
  for (var detName in detectorsAccumulator) {
    //gamma
    var gammaCount = {
      detName: detName,
      color: detectorsAccumulator[detName].gammaColor,
      counts: detectorsAccumulator[detName].gammaCounts,
    };
    // optional fields:
    if (detectorsAccumulator[detName].hasOwnProperty("gammaLiveTimes")) {
      gammaCount.liveTimes = detectorsAccumulator[detName].gammaLiveTimes;
    }

    out.gammaCounts.push(gammaCount);

    //neutron (optional)
    if (detectorsAccumulator[detName].hasOwnProperty("neutronCounts")) {
      var neutronCount = {
        detName: detName,
        color: detectorsAccumulator[detName].neutronColor,
        counts: detectorsAccumulator[detName].neutronCounts,
      };

      if (detectorsAccumulator[detName].hasOwnProperty("neutronLiveTimes")) {
        neutronCount.liveTimes = detectorsAccumulator[detName].neutronLiveTimes;
      }

      out.neutronCounts.push(neutronCount);
    } // if (detectorsAccumulator[detName].hasOwnProperty("neutronCounts"))
  } // for (var detName in detectorsAccumulator)

  return out;
};

// UNIMPLEMENTED

D3TimeChart.prototype.setHighlightRegions = function (regions) {
  //See the c++ function D3TimeChart::setHighlightRegionsToClient() for format of data
};

D3TimeChart.prototype.setXAxisTitle = function (title) {
  this.options.xtitle = title;
  //redraw x-title
};

D3TimeChart.prototype.setY1AxisTitle = function (title) {
  this.options.y1title = title;
  //redraw y1-title (e.g., gamma CPS axis title)
};

D3TimeChart.prototype.setY2AxisTitle = function () {
  this.options.y2title = title;
  //redraw y2-title (e.g., neutron CPS axis title)
};

D3TimeChart.prototype.setCompactXAxis = function (compact) {
  this.options.compactXAxis = compact;
  //Make x-zis title comapact or not
};

D3TimeChart.prototype.setGridX = function (show) {
  this.options.gridx = show;
  //add/remove horizantal grid lines
};

D3TimeChart.prototype.setGridY = function (show) {
  this.options.gridy = show;
  //add/remove vertical grid lines
};

D3TimeChart.prototype.setXAxisZoomSamples = function (firstsample, lastsample) {
  //Set it so only firstsample through lastsample sampels are visible
  //...
  //...

  //An example of using WtEmit:
  this.WtEmit(
    this.chart.id,
    { name: "timerangechange" },
    firstsample,
    lastsample
  );
};
