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

Sample = function (sampleNumber, cumRealTime, gammaCPS, neutronCPS) {
  this.sampleNumber = sampleNumber;
  this.cumRealTime = cumRealTime;
  this.gammaCPS = gammaCPS;
  this.neutronCPS = neutronCPS;
};

D3TimeChart = function (elem, options) {
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

  this.height = undefined;
  this.width = undefined;
  this.domains = undefined;

  this.svg = d3.select(this.chart).append("svg");

  this.axisBottomG = this.svg.append("g");
  this.axisLeftG = this.svg.append("g");
  this.lineG = this.svg.append("g");
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
  this.height = this.chart.clientWidth;
  this.width = this.chart.clientWidth;
};

D3TimeChart.prototype.render = function () {};

D3TimeChart.prototype.computeCumulativeRealTime = function (realTimes) {
  var cumRealTime = [];

  realTimes.reduce(function (acc, curr, i) {
    return (cumRealTime[i] = acc + curr);
  }, 0);

  return cumRealTime;
};

D3TimeChart.prototype.isValidData = function (data) {
  if (!data) {
    return false;
  }

  // check properties
  if (
    !data.hasOwnProperty("sampleNumbers") ||
    !data.hasOwnProperty("realTimes") ||
    !data.hasOwnProperty("gammaCounts") ||
    !data.hasOwnProperty("neutronCounts")
  ) {
    return false;
  }

  // Check data types
  if (
    !Array.isArray(data.gammaCounts) ||
    !Array.isArray(
      data.neutronCounts ||
        !Array.isArray(data.sampleNumbers) ||
        !Array.isArray(data.realTimes)
    )
  ) {
    return false;
  }

  var nSamples = data.sampleNumbers.length;

  // Check has at least one sample
  if (nSamples < 1) {
    return 0;
  }

  // Check matching lengths
  if (data.realTimes.length !== nSamples) {
    return false;
  }

  data.gammaCounts.forEach(function (detector) {
    if (
      !detector.counts ||
      !Array.isArray(detector.counts) ||
      detector.counts.length !== nSamples
    ) {
      return false;
    }
  });

  data.neutronCounts.forEach(function (detector) {
    if (
      !detector.counts ||
      !Array.isArray(detector.counts) ||
      detector.counts.length !== nSamples
    ) {
      return false;
    }
  });

  return true;
};

D3TimeChart.prototype.setData = function (data) {
  //See the c++ function D3TimeChart::setData()
  if (!this.isValidData(data)) {
    alert(
      "Error: Failed to set data-- structure of data is invalid. Aborting..."
    );
  } else {
    console.log(data);

    var formattedData = this.formatData(data);
    console.log(formattedData);

    this.data = formattedData;

    this.domains = this.getDomain(data);
  }
};

D3TimeChart.prototype.getDomain = function (data) {
  var xMax = Number.MIN_SAFE_INTEGER;
  var xMin = Number.MAX_SAFE_INTEGER;
  var yMax = Number.MIN_SAFE_INTEGER;
  var yMin = 0;

  // Find max
  // over gamma detectors
  var nSamples = data.sampleNumbers.length;

  for (var i = 0; i < data.gammaCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.gammaCounts[i].counts[j] / data.realTimes[j];

      if (cps > yMax) {
        yMax = cps;
      }
    }
  }

  // over neutron detectors
  for (var i = 0; i < data.neutronCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = data.neutronCounts[i].counts[j] / data.realTimes[j];

      if (cps > yMax) {
        yMax = cps;
      }
    }
  }

  var cumRealTime = this.computeCumulativeRealTime(data.realTimes);

  xMin = 0 - cumRealTime[0];
  xMax = cumRealTime[cumRealTime.length - 1] - cumRealTime[0];

  return {
    x: [xMin, xMax],
    y: [yMin, yMax],
  };
};

D3TimeChart.prototype.formatData = function (data) {
  var nSamples = data.sampleNumbers.length;
  var detectors = {};
  var samples = {};
  var cumRealTime = this.computeCumulativeRealTime(data.realTimes);

  if (data.gammaCounts) {
    data.gammaCounts.forEach(function (d) {
      if (!detectors.hasOwnProperty(d.detName)) {
        detectors[d.detName] = {
          gammaColor: undefined,
          neutronColor: undefined,
        };
      }
      detectors[d.detName].gammaColor = d.color;

      // if no data already present for this detector, create a new array to start holding the data for that detector and fill in data.
      if (!samples.hasOwnProperty(d.detName)) {
        samples[d.detName] = [];

        for (var i = 0; i < nSamples; i++) {
          samples[d.detName].push(
            new Sample(
              data.sampleNumbers[i],
              cumRealTime[i],
              d.counts[i] / data.realTimes[i],
              null
            )
          );
        }
      } else {
        // if (!samples.hasOwnProperty(d.detName))

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
        detectors[d.detName] = {
          gammaColor: undefined,
          neutronColor: undefined,
        };
      }

      detectors[d.detName].neutronColor = d.color;

      if (!samples.hasOwnProperty(d.detName)) {
        samples[d.detName] = [];

        for (var i = 0; i < nSamples; i++) {
          samples[d.detName].push(
            new Sample(
              data.sampleNumbers[i],
              cumRealTime[i],
              null,
              d.counts[i] / data.realTimes[i]
            )
          );
        }
      } else {
        // if (!samples.hasOwnProperty(d.detName))

        for (var i = 0; i < nSamples; i++) {
          samples[d.detName][i].neutronCPS = d.counts[i] / data.realTimes[i];
        }
      } // if (!samples.hasOwnProperty(d.detName))
    }); // data.neutronCounts.forEach
  } // if (data.neutronCounts)

  return { detectors: detectors, samples: samples };
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
