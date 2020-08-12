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
 * Constructor for Sample Objects
 * @param sampleNumber: integer sample number
 * @param realTimeInterval: Array of length 2 (i.e. [start_time, end_time])
 * @param gammaCPS: gamma counts per second
 * @param neutronCPS: neutron counts per second
 */
Sample = function (sampleNumber, realTimeInterval, gammaCPS, neutronCPS) {
  this.sampleNumber = sampleNumber;
  this.realTimeInterval = realTimeInterval;
  this.gammaCPS = gammaCPS;
  this.neutronCPS = neutronCPS;
};

/**
 * D3TimeChart object constructor.
 */
D3TimeChart = function (elem, options) {
  // this is called when the widget is loaded.
  var self = this;
  this.chart = typeof elem === "string" ? document.getElementById(elem) : elem;

  this.options = options || {};

  if (typeof this.options.xtitle !== "string") this.options.xlabel = "Time (s)";
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
  this.height = undefined;
  this.width = undefined;
  this.domains = undefined;
  this.leadTime = 2.5;

  this.margin = {
    top: 50,
    right: 50,
    bottom: 50,
    left: 50,
  };
  this.svg = d3.select(this.chart).append("svg");

  this.axisBottomG = this.svg.append("g").attr("class", "axis");
  this.axisLeftG = this.svg.append("g").attr("class", "axis");
  this.axisRightG = this.svg.append("g").attr("class", "axis");
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

D3TimeChart.prototype.handleResize = function () {
  // This function is called when the Wt layout manager resizes the parent <div> element
  // Need to redraw everything
  // ...
  // Make sure to update the C++ code of the changed plotting size.

  console.log("Resized!");
  console.log(
    "width: " +
      this.chart.clientWidth +
      ", " +
      "height: " +
      this.chart.clientHeight
  );
  this.height = this.chart.clientHeight;
  this.width = this.chart.clientWidth;

  this.render();
};

/**
 * Renders the D3TimeChart. Is called every time a resize occurs
 */
D3TimeChart.prototype.render = function () {
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

    // remove any existing drawn lines
    d3.selectAll(".line").remove();
    d3.selectAll(".axis_label").remove();

    // set dimensions of svg element and plot
    this.svg.attr("width", this.width);
    this.svg.attr("height", this.height);

    var { xScale, yScaleGamma, yScaleNeutron } = this.getScalers();

    var HAS_GAMMA = yScaleGamma ? true : false;
    var HAS_NEUTRON = yScaleNeutron ? true : false;

    // set different tick counts for different viewport breakpoints
    var nTicksX;
    if (this.width > 850) {
      nTicksX = this.data.nSamples;
    } else if (this.width > 520) {
      nTicksX = Math.floor(this.data.nSamples / 2);
    } else if (this.width > 280) {
      nTicksX = Math.floor(this.data.nSamples / 4);
    } else {
      nTicksX = Math.floor(this.data.nSamples / 8);
    }

    // plot axes
    var xAxis = d3.svg.axis().scale(xScale).ticks(nTicksX);
    this.axisBottomG
      .attr(
        "transform",
        "translate(0," + (this.height - this.margin.bottom) + ")"
      )
      .call(xAxis);

    if (HAS_GAMMA) {
      var yAxisLeft = d3.svg.axis().scale(yScaleGamma).ticks(3).orient("left");
      this.axisLeftG
        .attr("transform", "translate(" + this.margin.left + ",0)")
        .call(yAxisLeft);
    }

    if (HAS_NEUTRON) {
      var yAxisRight = d3.svg
        .axis()
        .scale(yScaleNeutron)
        .ticks(3)
        .orient("right");
      this.axisRightG
        .attr(
          "transform",
          "translate(" + (this.width - this.margin.left) + ",0)"
        )
        .call(yAxisRight);
    }

    // label axes
    this.svg
      .append("text")
      .attr("class", "axis_label")
      .attr(
        "transform",
        "translate(" +
          this.width / 2 +
          " ," +
          (this.height - this.margin.bottom + 30) +
          ")"
      )
      .style("text-anchor", "middle")
      .text("Real Time of Measurement (seconds)");

    this.svg
      .append("text")
      .attr("class", "axis_label")
      .attr("transform", `translate(15, ${this.height / 2}) rotate(-90) `)
      .style("text-anchor", "middle")
      .text("CPS");

    // plot data
    for (detector in this.data.samples) {
      var data = this.data.samples[detector];
      var lineGamma = "M";
      var lineNeutron = "M";

      data.forEach(function (d, i) {
        if (HAS_GAMMA) {
          var y0Gamma = yScaleGamma(d.gammaCPS);
        }
        if (HAS_NEUTRON) {
          var y0Neutron = yScaleNeutron(d.neutronCPS);
        }
        var x0 = xScale(d.realTimeInterval[1]);
        if (i === 0) {
          lineGamma += `${xScale(d.realTimeInterval[0])},${y0Gamma}H${x0}`;
          lineNeutron += `${xScale(d.realTimeInterval[0])},${y0Neutron}H${x0}`;
        } else {
          lineGamma += `H${x0}`;
          lineNeutron += `H${x0}`;
        }

        if (data[i + 1]) {
          if (HAS_GAMMA) {
            lineGamma += `V${yScaleGamma(data[i + 1].gammaCPS)}`;
          }
          if (HAS_NEUTRON) {
            lineNeutron += `V${yScaleNeutron(data[i + 1].neutronCPS)}`; //neturon
          }
        }
      });

      if (HAS_GAMMA) {
        this.svg
          .append("g")
          .attr("class", `line det_${detector.detName}`)
          .append("path")
          .style("stroke", this.data.detectors[detector].gammaColor)
          .style("fill", "none")
          .attr("d", lineGamma);
      }

      if (HAS_NEUTRON) {
        this.svg
          .append("g")
          .attr("class", `line det_${detector.detName}`)
          .append("path")
          .style("stroke", this.data.detectors[detector].neutronColor)
          .style("fill", "none")
          .attr("d", lineNeutron);
      }
    }
  }
};

/**
 * Get scaling functions based on the data domain and element dimensions.
 */
D3TimeChart.prototype.getScalers = function () {
  var xScale = d3.scale
    .linear()
    .domain(this.domains.x)
    .range([this.margin.left, this.width - this.margin.right]);
  var yScaleGamma = this.domains.yGamma
    ? d3.scale
        .linear()
        .domain(this.domains.yGamma)
        .range([this.height - this.margin.bottom, this.margin.top])
    : undefined;
  var yScaleNeutron = this.domains.yNeutron
    ? d3.scale
        .linear()
        .domain(this.domains.yNeutron)
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
D3TimeChart.prototype.getRealTimeIntervals = function (realTimes) {
  var realTimeIntervals = [];

  for (var i = 0; i < realTimes.length; i++) {
    if (i === 0) {
      // to handle long lead-in:
      var leadTime = Math.max(-realTimes[i], -this.leadTime);
      realTimeIntervals[i] = [leadTime, 0];
    } else {
      var prevEndpoint = realTimeIntervals[i - 1][1];
      realTimeIntervals[i] = [prevEndpoint, prevEndpoint + realTimes[i]];
    }
  }

  return realTimeIntervals;
};

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
    var formattedData = this.formatData(data);
    console.log(formattedData);

    this.data = formattedData;

    this.domains = this.getDomain(data);
    console.log("Domains:");
    console.log(this.domains);

    // if height and width are set, may render directly.
    if (this.height && this.width) {
      this.render();
    }
  }
};

/**
 * Gets data domain
 * @param {Object} data
 */
D3TimeChart.prototype.getDomain = function (data) {
  var realTimeIntervals = this.getRealTimeIntervals(data.realTimes);

  xMin = d3.min(realTimeIntervals, function (d) {
    return d[0];
  });
  xMax = d3.max(realTimeIntervals, function (d) {
    return d[1];
  });

  var yMaxGamma = Number.MIN_SAFE_INTEGER;
  var yMaxNeutron = Number.MIN_SAFE_INTEGER;

  // Find max over gamma detectors
  var nSamples = data.sampleNumbers.length;

  for (var i = 0; i < data.gammaCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.gammaCounts[i].counts[j] / data.realTimes[j];

      if (cps > yMaxGamma) {
        yMaxGamma = cps;
      }
    }
  }

  // Find max over neutron detectors
  for (var i = 0; i < data.neutronCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.neutronCounts[i].counts[j] / data.realTimes[j];

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
 * Formats data to several JSON objects for convenient plotting
 * e.g.
 * {
 *    detectors: { // for detector metadata
 *                  "det1": {
 *                            "isGammaDetector": true,
 *                            "isNeutronDetector": true,
 *                            "gammaColor": rgb(0,0,0),
 *                            "neutronColor": rgb(0,128,0)
 *                          },
 *                  "det2" : {...}
 *                  ...
 *               },
 *
 *    samples: {
 *                "det1": [{Sample}, {Sample}, {Sample}, ...],
 *                "det2": [{Sample}, {Sample}, {Sample}, ...],
 *                ...
 *             },
 *
 *    nSamples: 10
 * }
 * @param {*} data
 */
D3TimeChart.prototype.formatData = function (data) {
  var nSamples = data.sampleNumbers.length;
  var detectors = {};
  var samples = {};
  var realTimeIntervals = this.getRealTimeIntervals(data.realTimes);

  if (data.gammaCounts) {
    data.gammaCounts.forEach(function (d) {
      if (!detectors.hasOwnProperty(d.detName)) {
        // if new detector, initialize a new object
        detectors[d.detName] = {
          gammaColor: undefined,
          neutronColor: undefined,
          isGammaDetector: false,
          isNeutronDetector: false,
        };
      }
      detectors[d.detName].gammaColor = d.color;
      detectors[d.detName].isGammaDetector = true;

      // if no data already present for this detector, create a new array to start holding the data for that detector and fill in data.
      if (!samples.hasOwnProperty(d.detName)) {
        samples[d.detName] = [];

        for (var i = 0; i < nSamples; i++) {
          samples[d.detName].push(
            new Sample(
              data.sampleNumbers[i],
              realTimeIntervals[i],
              d.counts[i] / data.realTimes[i],
              null
            )
          );
        }
      } else {
        //data already present for this detector, so fill in the gamma counts
        for (var i = 0; i < nSamples; i++) {
          samples[d.detName][i].gammaCPS = d.counts[i] / data.realTimes[i];
        }
      } // if (!samples.hasOwnProperty(d.detName))
    }); // data.gammaCounts.forEach
  } // if (data.gammaCounts)

  if (data.neutronCounts) {
    data.neutronCounts.forEach(function (d) {
      if (!detectors.hasOwnProperty(d.detName)) {
        // if new detector, initialize a new object
        detectors[d.detName] = {
          gammaColor: undefined,
          neutronColor: undefined,
          isGammaDetector: false,
          isNeutronDetector: false,
        };
      }

      detectors[d.detName].neutronColor = d.color;
      detectors[d.detName].isNeutronDetector = true;

      if (!samples.hasOwnProperty(d.detName)) {
        samples[d.detName] = [];

        for (var i = 0; i < nSamples; i++) {
          samples[d.detName].push(
            new Sample(
              data.sampleNumbers[i],
              realTimeIntervals[i],
              null,
              d.counts[i] / data.realTimes[i]
            )
          );
        }
      } else {
        for (var i = 0; i < nSamples; i++) {
          samples[d.detName][i].neutronCPS = d.counts[i] / data.realTimes[i];
        }
      } // if (!samples.hasOwnProperty(d.detName))
    }); // data.neutronCounts.forEach
  } // if (data.neutronCounts)

  return { nSamples: nSamples, detectors: detectors, samples: samples };
};

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
