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
 * Constructor for custom error type, inherited from generic JS Error. Arguments are the same as for the generic JS Error objects.
 */
ValidationError = function (message, fileName, lineNumber) {
  var instance = new Error(message, fileName, lineNumber);
  Object.setPrototypeOf(instance, Object.getPrototypeOf(this));
  return instance;
};
ValidationError.prototype = Object.create(Error.prototype, {
  constructor: {
    value: Error,
    enumerable: false,
    writable: true,
    configurable: true,
  },
});
if (Object.setPrototypeOf) {
  Object.setPrototypeOf(ValidationError, Error);
} else {
  ValidationError.__proto__ = Error;
}

/**
 * Constructor for DataPoint objects. Represents a single point in the time history chart.
 * @param {Number} time: time of measurement (x)
 * @param {Number} gammaCPS: gamma counts per second (y)
 * @param {Number} neutronCPS: neutron counts per second (y)
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

/**
 * Constructor for DetectorMetaData objects.
 * @param {String} gammaColor: string representation of a color to use for the plotted gamma line
 * @param {String} neutronColor: string representation of a color to use for the plotted neutron line
 */
DetectorMetaData = function (gammaColor, neutronColor) {
  this.gammaColor = gammaColor;
  this.neutronColor = neutronColor;
};

DetectorMetaData.prototype.setGammaColor = function (gammaColor) {
  this.gammaColor = gammaColor;
};

DetectorMetaData.prototype.setNeutronColor = function (neutronColor) {
  this.neutronColor = neutronColor;
};

/**
 * Constructor for rudimentary 1-dimensional brush along x direction. Used to map mouse (drag) selection ranges to the data domain.
 * @param {Object} scale: a D3 scale object that sets the x-scale associated with this brush. The scale is used to map mouse positions into the data domain used to create the scale.
 */
BrushX = function (scale) {
  this.start = null;
  this.end = null;
  this.scale = scale;
};

/**
 * Sets or updates the x-scale associated with this brush. The scale is used to map mouse positions into data domain used to create the scale.
 * @param {Object} scale: D3 scale object
 */
BrushX.prototype.setScale = function (scale) {
  this.scale = scale;
  return this;
};

/**
 * Returns the x-scale associated with this brush. The scale is used to map mouse positions into the data domain used to create the scale.
 */
BrushX.prototype.getScale = function () {
  return this.scale;
};

/**
 * Returns true if and only if the brush extent is empty. A brush extent is empty if it has no start point or end point.
 * When a brush is created, it is initially empty; the brush may also become empty with a single click on the background without moving, or if the extent is cleared.
 * A brush is considered empty if it has zero-width. When the brush is empty, its extent is not strictly defined.
 */
BrushX.prototype.empty = function () {
  return !this.start || !this.end || this.start === this.end;
};

/**
 * Clears the extent, making the brush extent empty.
 */
BrushX.prototype.clear = function () {
  this.start = null;
  this.end = null;
};

/**
 * Sets the start value of the brush by mapping the mouse coordinate into into the data domain used to create the brush's scale.
 * @param {Number} startCoord : d3 mouse x coordinate.
 */
BrushX.prototype.setStart = function (startCoord) {
  if (!this.scale) {
    throw new Error("Brush scale has not been set!");
  }
  var scaledStart = this.scale.invert(startCoord);
  var domain = this.scale.domain();
  if (scaledStart < domain[0]) {
    scaledStart = domain[0];
  }
  if (scaledStart > domain[1]) {
    scaledStart = domain[1];
  }
  this.start = scaledStart;
};

BrushX.prototype.getStart = function () {
  return this.start;
};

/**
 * Sets the end value of the brush by mapping the mouse coordinate into into the data domain used to create the brush's scale.
 * @param {Number} endCoord : d3 mouse x coordinate.
 */
BrushX.prototype.setEnd = function (endCoord) {
  if (!this.scale) {
    throw new Error("Brush scale has not been set!");
  }

  var scaledEnd = this.scale.invert(endCoord);
  var domain = this.scale.domain();

  if (scaledEnd < domain[0]) {
    scaledEnd = domain[0];
  }
  if (scaledEnd > domain[1]) {
    scaledEnd = domain[1];
  }
  this.end = scaledEnd;
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

  /** OPTIONS */
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

  // minimum selection width option. Default value taken from D3SpectrumChart
  if (typeof this.options.minSelectionWidth !== "number")
    this.options.minSelectionWidth = 8;

  // option to use simplified gesture mode, which we currently define as mapping all drag gestures to zoom functionality.
  // current implementation of simplified gesture mode simply disables use of this.highlightOptions and treats any drag gesture as zoom, regardless of keyboard modifier.
  // To add extra functionality to other keys in simplified gesture mode, would need to use and make conditional changes to this.highlightOptions to add new key mappings under simplifiedGestureMode
  if (typeof this.options.useSimplifiedGestures !== "boolean")
    this.options.useSimplifiedGestures = false;

  // defines keyboard modifier keys and/or other metadata to use with drag gestures to achieve various highlight selection functionalities.
  // uses JS object to emulate behavior of Set object (e.g. {"none": true} instead of new Set(["none"])) to accommodate ES5
  this.highlightOptions = {
    foreground: {
      modifierKey: { none: true },
    },
    background: {
      modifierKey: { altKey: true },
    },
    zoom: {
      modifierKey: { ctrlKey: true, rightClick: true },
    },
  };

  // defines margin dimensions used in rendering of chart.
  this.margin = {
    top: 5,
    right: 60,
    bottom: 50,
    left: 60,
  };

  /** CONSTANTS */
  // map<integer, string> (i.e. array) of source integer code to string. 0 == IntrinsicActivity, 1== Calibration, 2 == Background, 3 == Foreground, 4 == Unknown
  this.SOURCE_MAP = Object.freeze([
    "IntrinsicActivity",
    "Calibration",
    "Background",
    "Foreground",
    "Unknown",
  ]);

  // colors used for highlight rectangles for various selection types.
  this.HIGHLIGHT_COLORS = Object.freeze({
    foreground: "rgb(255, 255, 0)",
    background: "rgb(0, 255, 255)",
    zoom: "rgb(0, 0, 0)",
  });

  /** COMPONENT STATE */
  // component state. When state changes, the component should usually respond by re-rendering/updating to reflect changes.
  this.state = {
    data: {
      formatted: null,
      raw: null,
      backgroundDuration: null,
      sampleNumberToIndexMap: null,
      unzoomedCompressionIndex: 0,
    },
    selection: null, // maybe would have been better to have named this "zoom": stores data related to zoom selection (e.g. data domain of magnified area, corresponding compression index to use for plotting this magnified data)
    regions: null,
    brush: new BrushX(),
    height: null,
    width: null,
  };

  /** SVG COMPONENT REFERENCES */
  this.svg = d3.select(this.chart).append("svg");
  this.verticalGridG = this.svg.append("g").attr("class", "grid");
  this.horizontalLeftGridG = this.svg.append("g").attr("class", "grid");
  this.horizontalRightGridG = this.svg.append("g").attr("class", "grid");
  this.linesG = this.svg.append("g").attr("class", "lines");
  this.axisBottomG = this.svg.append("g").attr("class", "axis xaxis");
  this.axisLeftG = this.svg.append("g").attr("class", "axis yaxis");
  this.axisRightG = this.svg.append("g").attr("class", "axis yaxis");

  this.highlightRegionsG = this.svg
    .append("g")
    .attr("class", "highlight_region");

  this.rectG = this.svg.append("g").attr("class", "interaction_area");
  this.highlightText = this.rectG
    .append("text")
    .attr("font-size", 11)
    .attr("font-weight", "bold");
  this.highlightRect = this.rectG.append("rect").attr("class", "selection");
  this.rect = this.rectG.append("rect"); //rectangle spanning the interactable area on plot area
  this.bottomAxisRect = this.rectG.append("rect"); // rectangle spanning additional interactable area on bottom axis area

  this.hoverToolTip = d3
    .select(this.chart)
    .append("div")
    .attr("id", "hover_info")
    .attr("class", "tooltip")
    .style("left", this.margin.left + 20 + "px")
    .style("top", this.margin.top + "px");

  this.occupancyLinesG = this.svg.append("g").attr("class", "occupancy_lines");

  /** MISC MEMBERS */
  this.escapeKeyPressed = false;
  this.shiftKeyHeld = false;
  this.highlightModifier = null; // holds the key pressed in conjunction with a highlight gesture to modify the action
  this.draggedForward = false;

  /** GLOBAL LISTENERS */
  // listeners to support esc canceling of highlighting and arrow-key panning
  document.addEventListener("keydown", function (evt) {
    evt = evt || window.event;
    if (evt.key === "Escape") {
      self.escapeKeyPressed = true;
      d3.select("body").style("cursor", "auto");
      self.mouseUpHighlight();
    } else if (evt.key === "ArrowLeft") {
      self.shiftSelection(-1);
    } else if (evt.key === "ArrowRight") {
      self.shiftSelection(1);
    }
  });
};

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
 * Sets this.state.data.raw, this.state.data.formatted, this.state.selection, and this.state.data.sampleNumberToIndexMap
 * @param {Object} rawData : raw data object sent from Wt
 */
D3TimeChart.prototype.setData = function (rawData) {
  try {
    //See the c++ function D3TimeChart::setData()
    if (!this.isValidRawData(rawData)) {
      throw new ValidationError("Structure of raw data is not valid.");
    }
    var formattedData = this.formatDataFromRaw(rawData);

    this.state.data.formatted = [formattedData];
    // console.log(this.state.data.formatted);

    // create inverted index of sample numbers  for fast lookup of array-indices from sample number keys
    var sampleToIndexMap = {};
    for (var i = 0; i < rawData.sampleNumbers.length; i++) {
      sampleToIndexMap[rawData.sampleNumbers[i]] = i;
    }

    // set other data members

    this.state.data.raw = rawData;
    // console.log(sampleToIndexMap);
    this.state.data.sampleNumberToIndexMap = sampleToIndexMap;

    // clear existing selection if there is any
    this.state.selection = null;

    // clear existing brush if there is any
    this.state.brush = new BrushX();

    // clear existing regions if there are any
    this.state.regions = null;

    // clear any existing plot lines drawn
    this.linesG.selectAll("path").remove();

    // if height and width are set, may initialize chart directly.
    if (this.state.height && this.state.width) {
      this.reinitializeChart();
    }
  } catch (err) {
    if (err instanceof ValidationError) {
      console.log(err.message + "\nDoing nothing...");
      // console.log(err); // log call stack to console
    } else {
      throw err;
    }
  }
};

D3TimeChart.prototype.handleResize = function () {
  // This function is called when the Wt layout manager resizes the parent <div> element
  // Need to redraw everything (incl size of svg element, )
  // ...
  // Make sure to update the C++ code of the changed plotting size.
  try {
    // console.log("Resized!");
    this.state.height = this.chart.clientHeight;
    this.state.width = this.chart.clientWidth;

    this.reinitializeChart();
  } catch (err) {
    if (err instanceof ValidationError) {
      console.log(err.message + "\nDoing nothing...");
      // console.log(err); // log call stack to console
    } else {
      throw err;
    }
  }
};

/**
 * Renders/re-initializes the D3TimeChart to be set up for plotting data and interaction. Compresses data for purposes of displaying, updates dimensions of svg element, updates dimensions of clip-path, and defines drag behavior over the figure
 * depends on data and chart size being already defined
 * @param {Object} options : Optional argument to specify render options
 */
D3TimeChart.prototype.reinitializeChart = function (options) {
  if (!this.state.data.formatted) {
    throw new ValidationError("D3TimeChart data is not set.");
  }
  if (!this.state.height || !this.state.width) {
    throw new ValidationError(
      "dimensions of D3TimeChart div element are not set."
    );
  }
  // console.log("Re-initializing...");

  this.margin.bottom = this.options.compactXAxis ? 25 : 50;

  var plotWidth = this.state.width - this.margin.left - this.margin.right;
  var plotHeight = this.state.height - this.margin.top - this.margin.bottom;

  var nPoints = this.state.data.formatted[0].sampleNumbers.length;

  // check chart pixels vs full set of  data points. Compress data if needed for purposes of rendering, and cache the compressed data inside the this.state.data.formatted array.
  // each data[i] is the data compressed at level 2^i.

  // impose lower plotWidth limit on compression calculations
  if (plotWidth > 10) {
    var compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));
    if (plotWidth < nPoints) {
      var i = 1;
      for (var i = 1; i <= compressionIndex; i++) {
        // only compress if data doesn't already exist before
        if (this.state.data.formatted[i] == null) {
          // console.log("Compressing!");
          this.state.data.formatted[i] = this.formatDataFromRaw(
            this.compress(this.state.data.raw, Math.pow(2, i))
          );
        }
      }
    } // if (plotWidth < nPoints)
    // console.log(this.state.data.formatted);
    this.state.data.unzoomedCompressionIndex = compressionIndex;
  }

  // set dimensions of svg element and plot
  this.svg.attr("width", this.state.width).attr("height", this.state.height);

  // set dimensions of interactable area.
  this.rect
    .attr("width", plotWidth)
    .attr("height", plotHeight)
    .attr("x", this.margin.left)
    .attr("y", this.margin.top)
    .attr("fill-opacity", 0);

  this.bottomAxisRect
    .attr("width", plotWidth)
    .attr("height", this.axisBottomG.node().getBBox().height)
    .attr("x", this.margin.left)
    .attr("y", this.margin.top + plotHeight)
    .attr("fill-opacity", 0);

  // add a clipPath: everything outside of this area will not be drawn
  var clip = this.svg.select("#clip_th");

  if (!clip.empty()) {
    // update if already exists
    clip
      .select("rect")
      .attr("width", plotWidth)
      .attr("height", plotHeight)
      .attr("x", this.margin.left)
      .attr("y", this.margin.top);
  } else {
    // else create new
    this.svg
      .append("defs")
      .append("svg:clipPath")
      .attr("id", "clip_th")
      .append("svg:rect")
      .attr("width", plotWidth)
      .attr("height", plotHeight)
      .attr("x", this.margin.left)
      .attr("y", this.margin.top);
  }

  this.linesG.attr("clip-path", "url(#clip_th)");
  this.highlightRegionsG.attr("clip-path", "url(#clip_th");
  this.occupancyLinesG.attr("clip-path", "url(#clip_th");

  // if have a selection, choose compression index to use for plotting based on number of points in the selection.
  if (this.state.selection) {
    var leftIndex = this.findDataIndex(this.state.selection.domain[0], 0);
    var rightIndex = this.findDataIndex(this.state.selection.domain[1], 0);
    var nPointsSelection = rightIndex - leftIndex + 1;

    this.state.selection.compressionIndex = Math.ceil(
      Math.log2(Math.ceil(nPointsSelection / plotWidth))
    );
  }

  // get scales
  var domains = this.state.selection
    ? {
        x: this.state.selection.domain,
        yGamma:
          this.state.data.formatted[this.state.selection.compressionIndex]
            .domains.yGamma,
        yNeutron:
          this.state.data.formatted[this.state.selection.compressionIndex]
            .domains.yNeutron,
      }
    : this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .domains;

  var scales = this.getScales(domains);

  var compressionIndex = this.state.selection
    ? this.state.selection.compressionIndex
    : this.state.data.unzoomedCompressionIndex;

  var brush = this.state.brush;

  // set scale of brush
  brush.setScale(scales.xScale);

  // DEFINE HANDLERS AND LISTENERS
  // TODO: To add additional touch functionality, add analogous touch gestures (to right click, alt key, ctrl key-- to perform zoom and background selection)

  /**
   * Handler for the initiation of a selection gesture.
   */
  var startSelection = () => {
    if (window.TouchEvent && d3.event instanceof TouchEvent) {
      d3.event.preventDefault();
      d3.event.stopPropagation();
    }

    if (this.escapeKeyPressed) this.escapeKeyPressed = false;

    var coords =
      window.TouchEvent &&
      (d3.event instanceof TouchEvent ||
        d3.event.sourceEvent instanceof TouchEvent)
        ? d3.touches(this.rect.node())[0]
        : d3.mouse(this.rect.node());
    // console.log(coords);
    brush.setStart(coords[0]);
    d3.select("body").style("cursor", "move");
    // console.log(d3.event.sourceEvent);

    // TODO: add analogous touch gestures to add additional touch functionality
    var TOUCH_ANALOGOUS_SHIFT = false;
    var TOUCH_ANALOGOUS_RIGHTCLICK = false;
    var TOUCH_ANALOGOUS_ALTKEY = false;
    var TOUCH_ANALOGOUS_CTRLCLICK = false;

    this.shiftKeyHeld =
      TOUCH_ANALOGOUS_SHIFT ||
      (d3.event.sourceEvent && d3.event.sourceEvent.shiftKey);

    if (
      TOUCH_ANALOGOUS_RIGHTCLICK ||
      (d3.event.type == "dragstart" &&
        d3.event.sourceEvent instanceof MouseEvent &&
        d3.event.sourceEvent.button == 2)
    ) {
      this.highlightModifier = "rightClick";
      this.mouseDownHighlight(coords[0], "rightClick");
    } else if (
      TOUCH_ANALOGOUS_ALTKEY ||
      (d3.event.type == "dragstart" &&
        d3.event.sourceEvent instanceof MouseEvent &&
        d3.event.sourceEvent.altKey)
    ) {
      this.highlightModifier = "altKey";
      this.mouseDownHighlight(coords[0], "altKey");
    } else if (
      TOUCH_ANALOGOUS_CTRLCLICK ||
      (d3.event.type == "dragstart" &&
        d3.event.sourceEvent instanceof MouseEvent &&
        d3.event.sourceEvent.ctrlKey)
    ) {
      this.highlightModifier = "ctrlKey";
      this.mouseDownHighlight(coords[0], "ctrlKey");
    } else {
      this.highlightModifier = "none";
      this.mouseDownHighlight(coords[0], "none");
    }
  };

  /**
   * Handler for the progression of a selection gesture.
   */
  var moveSelection = () => {
    if (window.TouchEvent && d3.event instanceof TouchEvent) {
      d3.event.preventDefault();
      d3.event.stopPropagation();
    }

    if (!this.escapeKeyPressed) {
      var coords =
        window.TouchEvent &&
        (d3.event instanceof TouchEvent ||
          d3.event.sourceEvent instanceof TouchEvent)
          ? d3.touches(this.rect.node())[0]
          : d3.mouse(this.rect.node());
      brush.setEnd(coords[0]);

      if (
        this.options.useSimplifiedGestures ||
        this.highlightModifier in this.highlightOptions.zoom.modifierKey
      ) {
        // if brush backward, call handler to handle zoom-out. Else, handle drawing the selection rectangle for zoom-in.
        if (brush.getEnd() < brush.getStart()) {
          this.handleDragBackZoom();
        } else {
          this.handleDragForwardZoom();
        }
      } else {
        // handle interactions other than zoom
        if (!this.options.useSimplifiedGestures) {
          // unnecessary check, but added to make it clear that if you wanted to add extra functionality to "simple gesture" mode, then you should handle things differently.
          // handle foreground or background selection

          var width =
            brush.getScale()(brush.getEnd()) -
            brush.getScale()(brush.getStart());
          this.mouseMoveHighlight(width);
        }
      }
    }
  };

  /**
   * Handler for the termination of a selection gesture.
   */

  var endSelection = () => {
    if (window.TouchEvent && d3.event instanceof TouchEvent) {
      d3.event.preventDefault();
      d3.event.stopPropagation();
    }
    if (brush.extent() != null) {
      if (this.escapeKeyPressed) {
        // clear selections and reset escape key
        brush.clear();
        this.escapeKeyPressed = false;
      } else {
        d3.select("body").style("cursor", "auto");
        // console.log(brush.extent());
        var lIdx = this.findDataIndex(brush.extent()[0], 0);
        var rIdx = this.findDataIndex(brush.extent()[1], 0);
        // console.log([
        //   this.state.data.formatted[0].sampleNumbers[lIdx],
        //   this.state.data.formatted[0].sampleNumbers[rIdx],
        // ]);
        if (
          this.options.useSimplifiedGestures ||
          this.highlightModifier in this.highlightOptions.zoom.modifierKey
        ) {
          this.handleBrushZoom();
        } else {
          // handle interactions other than zoom
          if (!this.options.useSimplifiedGestures) {
            // unnecessary check, but added to make it clear that if you wanted to add extra functionality to "simple gesture" mode, then you should handle things differently.
            // handle foreground or background selection

            var keyModifierMap = {
              altKey: 0x4,
              shiftKey: 0x1,
              none: 0x0,
            };
            this.WtEmit(
              this.chart.id,
              { name: "timedragged" },
              this.state.data.formatted[0].sampleNumbers[lIdx],
              this.state.data.formatted[0].sampleNumbers[rIdx],
              keyModifierMap[this.highlightModifier] |
                (keyModifierMap["shiftKey"] & this.shiftKeyHeld) // bitwise OR with the shift key modifier if held, 0 otherwise.
            );
          }
        }
      }
    }
    // clear
    this.mouseUpHighlight();
    brush.clear();
    this.highlightModifier = null;
    this.shiftKeyHeld = false;
  };

  // touch drag behavior
  this.rect
    .on("touchstart", startSelection)
    .on("touchmove", moveSelection)
    .on("touchend", endSelection);

  // mouse drag behavior
  var selectionDrag = d3.behavior
    .drag()
    .on("dragstart", startSelection)
    .on("drag", moveSelection)
    .on("dragend", endSelection);
  this.rect.call(selectionDrag);

  // pan drag behavior
  // initialize new variables for holding new selection and scales from panning
  var newSelection = null;
  var newScale = brush.getScale();
  var panDrag = d3.behavior
    .drag()
    .on("dragstart", () => {
      newSelection = this.state.selection;
      var coords = d3.mouse(this.rect.node());
      brush.setStart(coords[0]);
      d3.select("body").style("cursor", "ew-resize");
    })
    .on("drag", () => {
      if (newSelection) {
        brush.setEnd(d3.mouse(this.rect.node())[0]);
        var res = this.handleBrushPanSelection();
        if (res) {
          newSelection = res.newSelection;
          newScale = res.newScale;
        }
      }
    })
    .on("dragend", () => {
      // update selection, update scale
      if (newSelection) {
        this.state.selection = newSelection;
        brush.setScale(newScale.xScale);
        this.updateChart(newScale, newSelection.compressionIndex);
      }
      brush.clear();
    });

  this.bottomAxisRect.call(panDrag);

  this.bottomAxisRect
    .on("mouseover", () => {
      d3.select("body").style("cursor", "ew-resize");
    })
    .on("mouseout", () => {
      d3.select("body").style("cursor", "auto");
    });

  // mouse wheel behavior
  this.rect.node().onwheel = (evt) => {
    evt.preventDefault();
    this.handleMouseWheel(evt.deltaX, evt.deltaY, evt.x);
  };

  // Finally, initialization finished-- update chart
  this.updateChart(scales, compressionIndex, options);
};

/**
 * Draws or updates plots and axes
 * @param {Object} scales : Object of scales with properties xScale, yScaleGamma, yScaleNeutron
 * @param {Number} compressionIndex : Positive integer specifying data compression level to use for redrawing the chart. Compressed data are calculated during render() and cached inside this.state.data.formatted
 * @param {Object} options : Optional argument to specify additional render options if any
 */
D3TimeChart.prototype.updateChart = function (
  scales,
  compressionIndex,
  options
) {
  var xScale = scales.xScale;
  var yScaleGamma = scales.yScaleGamma;
  var yScaleNeutron = scales.yScaleNeutron;

  var HAS_GAMMA = true;
  var HAS_NEUTRON = false;

  var plotWidth = this.state.width - this.margin.left - this.margin.right;
  var plotHeight = this.state.height - this.margin.top - this.margin.bottom;

  // add/update hover interaction based on current scale (or zoom)
  this.rect
    .on("mouseover", () => {
      this.showToolTip();
    })
    .on("mousemove", () => {
      var x = xScale.invert(d3.mouse(this.rect.node())[0]);

      var idx = this.findDataIndex(x, compressionIndex);

      var startTimeStamp = this.state.data.formatted[compressionIndex]
        .startTimeStamps
        ? this.state.data.formatted[compressionIndex].startTimeStamps[idx]
        : null;

      var sourceType = this.state.data.formatted[compressionIndex].sourceTypes
        ? this.state.data.formatted[compressionIndex].sourceTypes[idx]
        : null;

      var sampleNumber =
        this.state.data.formatted[compressionIndex].sampleNumbers[idx];

      var tooltipData = [];
      for (var detName in this.state.data.formatted[compressionIndex]
        .detectors) {
        var y =
          this.state.data.formatted[compressionIndex].detectors[detName].counts[
            idx * 2
          ];

        tooltipData.push({
          detName: detName,
          gammaCPS: y.gammaCPS,
          neutronCPS: y.neutronCPS,
          startTimeStamp: startTimeStamp,
          sourceType: sourceType,
          sampleNumber: sampleNumber,
        });
      }
      var optargs = { sourceType: sourceType, startTimeStamp: startTimeStamp };
      this.updateToolTip(x, tooltipData, optargs);
    })
    .on("mouseout", () => {
      this.hideToolTip();
    });

  // plot data
  for (var detName in this.state.data.formatted[compressionIndex].detectors) {
    var counts =
      this.state.data.formatted[compressionIndex].detectors[detName].counts;

    // only use visible range if zoomed in, otherwise use full range.
    if (this.state.selection) {
      var lIdx = this.findDataIndex(
        this.state.selection.domain[0],
        compressionIndex
      );
      var rIdx = this.findDataIndex(
        this.state.selection.domain[1],
        compressionIndex
      );
      counts = counts.slice(lIdx * 2, (rIdx + 1) * 2);
    }

    var meta =
      this.state.data.formatted[compressionIndex].detectors[detName].meta;

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
      pathGamma.datum(counts).attr("d", lineGamma);
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
        pathNeutron.datum(counts).attr("d", lineNeutron);
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

  // update highlight region rendering
  if (this.state.regions && this.state.data.sampleNumberToIndexMap) {
    // console.log(this.state.regions);
    var chart = this;
    this.highlightRegionsG.selectAll("rect").each(function (d, i) {
      var startSample = chart.state.regions[i].startSample;
      var endSample = chart.state.regions[i].endSample;
      // console.log(endSample);
      var lIdx = chart.state.data.sampleNumberToIndexMap[startSample];
      var rIdx = chart.state.data.sampleNumberToIndexMap[endSample];
      // if (
      //   !(
      //     lIdx == null ||
      //     rIdx == null ||
      //     lIdx < 0 ||
      //     rIdx < 0 ||
      //     !Array.isArray(chart.state.data.formatted[0].realTimeIntervals) ||
      //     lIdx >= chart.state.data.formatted[0].realTimeIntervals.length ||
      //     rIdx >= chart.state.data.formatted[0].realTimeIntervals.length ||
      //     !Array.isArray(chart.state.data.formatted[0].realTimeIntervals[lIdx]) ||
      //     !Array.isArray(chart.state.data.formatted[0].realTimeIntervals[rIdx]) ||
      //     chart.state.data.formatted[0].realTimeIntervals[lIdx].length < 2 ||
      //     chart.state.data.formatted[0].realTimeIntervals[rIdx].length < 2
      //   )
      // ) {
      var startTime = chart.state.data.formatted[0].realTimeIntervals[lIdx][0];
      var endTime = chart.state.data.formatted[0].realTimeIntervals[rIdx][1];

      var lPixel = xScale(startTime);
      var rPixel = xScale(endTime);

      var highlightWidth = rPixel - lPixel > 2 ? rPixel - lPixel : 2;
      d3.select(this).attr("x", lPixel).attr("width", highlightWidth);
      // }
    });
  }

  // plot axes and labels

  var tickCount = 20;

  do {
    tickCount = Math.floor(tickCount / 2);

    // implement lower limit on tick count
    if (tickCount < 2) {
      break;
    }
    var xTickVals = xScale.ticks(tickCount);
    var xStart = xTickVals[0];
    var xStep = xTickVals[1] - xTickVals[0];
    var xTickCount = xTickVals.length;

    var chartDomain = this.state.selection
      ? this.state.selection.domain
      : this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
          .domains.x;
    var xTicksGenerated = this.generateTicks(
      xStart,
      xStep,
      xTickCount,
      10
    ).filter(
      (tick) => tick.value >= chartDomain[0] && tick.value <= chartDomain[1]
    );

    var xAxis = d3.svg
      .axis()
      .scale(xScale)
      .tickValues(xTicksGenerated.map((tick) => tick.value))
      .tickFormat(d3.format("~.6f"));

    // update or create axis
    this.axisBottomG
      .attr(
        "transform",
        "translate(0," + (this.state.height - this.margin.bottom) + ")"
      )
      .attr("id", "th_x-axis")
      .call(xAxis);

    // check whether there is any possibility for axis text overlap
    // if yes, then re-define axes with reduced (half) the current tick count
    var NUMBER_OF_TICKS_BETWEEN_MAJOR_TICKS = 5; // set this to the number of ticks between major ticks. Helps you get the next visible tick.

    var renderedTickCount = this.axisBottomG.selectAll("g.tick").size();
    var lastVisibleTick = this.axisBottomG.select(
      `g.tick:nth-child(${renderedTickCount})`
    );

    var secondToLastVisibleTick = this.axisBottomG.select(
      `g.tick:nth-child(${
        renderedTickCount - NUMBER_OF_TICKS_BETWEEN_MAJOR_TICKS
      })`
    );

    if (lastVisibleTick.empty() || secondToLastVisibleTick.empty()) {
      break;
    }

    var lastVisibleTickBound = lastVisibleTick.node().getBoundingClientRect();

    var secondToLastVisibleTickBound = secondToLastVisibleTick
      .node()
      .getBoundingClientRect();
  } while (secondToLastVisibleTickBound.right > lastVisibleTickBound.left);

  var xAxisArrowDefs = this.axisBottomG.select("#arrow_defs");

  if (xAxisArrowDefs.empty()) {
    xAxisArrowDefs = this.axisBottomG
      .append("svg:defs")
      .attr("id", "arrow_defs");

    xAxisArrowDefs
      .append("svg:marker")
      .attr("id", "right_arrow")
      .attr("class", "xaxisarrow")
      .attr("refX", -2)
      .attr("refY", 11)
      .attr("markerWidth", 9)
      .attr("markerHeight", 14)
      .attr("orient", 0)
      .append("path")
      .attr("d", "M0,0 L0,10 L5,5 L0,0")
      .attr("transform", "translate(0, 0)")
      .style("stroke", "black")
      .style("fill", "black");

    xAxisArrowDefs
      .append("svg:marker")
      .attr("id", "left_arrow")
      .attr("class", "xaxisarrow")
      .attr("refX", -2)
      .attr("refY", -1)
      .attr("markerWidth", 9)
      .attr("markerHeight", 14)
      .attr("orient", 180)
      .append("path")
      .attr("d", "M0,0 L0,10 L5,5 L0,0")
      .attr("transform", "translate(0, 0)")
      .style("stroke", "black")
      .style("fill", "black");
  }

  // if zoomed in, display zoom in marker.
  // Uses same zoom-in marker defined in SpectrumChartD3.js. Can define a different one in this file and use it if want to break the dependency.
  var axisBottomPath = this.axisBottomG.select("path");
  if (this.state.selection) {
    axisBottomPath.attr("marker-end", "url(#right_arrow)");
    axisBottomPath.attr("marker-start", "url(#left_arrow)");
  } else {
    axisBottomPath.attr("marker-end", null);
    axisBottomPath.attr("marker-start", null);
  }

  var axisLabelX = this.svg.select("#th_label_x");

  if (axisLabelX.empty()) {
    // create the element if not already created
    axisLabelX = this.svg
      .append("text")
      .attr("class", "xaxistitle")
      .attr("id", "th_label_x")
      .style("text-anchor", "middle")
      .text(this.options.xtitle);
  } // if (!axisLabelX.empty())

  // compute axis label translation to use
  var axisLabelXTranslation = this.options.compactXAxis
    ? "translate(" +
      (this.margin.left +
        this.state.width -
        axisLabelX.node().getBBox().width -
        20) +
      "," +
      (this.state.height -
        this.margin.bottom +
        this.axisBottomG.node().getBBox().height -
        3) +
      ")"
    : "translate(" +
      this.state.width / 2 +
      "," +
      (this.state.height -
        this.margin.bottom +
        this.axisBottomG.node().getBBox().height +
        15) +
      ")";

  axisLabelX.attr("transform", axisLabelXTranslation);

  var xLabelBoundingRect = axisLabelX.node().getBoundingClientRect();

  var USE_COMPACT_X_AXIS = this.options.compactXAxis;
  var axisBottomTicks = this.axisBottomG.selectAll("g.tick");

  // hide tick labels if they overlap with compact x-axis label
  axisBottomTicks.each(function () {
    var tickTransform = d3.transform(d3.select(this).attr("transform"));
    var tickText = d3.select(this).select("text");
    var tickBoundingRect = tickText.node().getBoundingClientRect();
    if (
      USE_COMPACT_X_AXIS &&
      tickBoundingRect.right > xLabelBoundingRect.left &&
      tickBoundingRect.left < xLabelBoundingRect.right
    ) {
      tickText.attr("visibility", "hidden");
    } else {
      tickText.attr("visibility", "visible");
    }
  });

  var dataBackgroundDuration = this.state.data.backgroundDuration;
  var firstTick = this.axisBottomG.select("g.tick:first-child");

  if (!firstTick.empty() && dataBackgroundDuration != null) {
    // add background duration and remove negative axis labels
    var leftMargin = this.margin.left;
    var axisBottomTicks = this.axisBottomG.selectAll("g.tick");
    axisBottomTicks.each(function () {
      var text = d3.select(this).select("text");
      var line = d3.select(this).select("line");
      if (parseFloat(text.node().textContent) < 0) {
        if (text.node().textContent === firstTick.node().textContent) {
          d3.select(this).attr("transform", "translate(" + leftMargin + ",0)");
          text.node().textContent = -dataBackgroundDuration;
        } else {
          d3.select(line.node()).attr("visibility", "hidden");
          d3.select(text.node()).attr("visibility", "hidden");
        }
      } // if (text.node().textContent < 0)
    }); // axisBottomTicks.each()
  } // if (dataBackgroundDuration != null)

  if (this.state.data.raw.occupancies) {
    var occupancies = this.state.data.raw.occupancies;
    if (
      this.occupancyLinesG.selectAll(".occupancy_line_group").size() !=
      occupancies.length
    ) {
      // if for any reason don't have all lines drawn, clear all existing lines and redraw
      this.occupancyLinesG.selectAll(".occupancy_line_group").remove();

      for (var i = 0; i < occupancies.length; i++) {
        var occupancyG = this.occupancyLinesG
          .append("g")
          .attr("class", "occupancy_line_group");
        var occupancyStartG = occupancyG
          .append("g")
          .attr("class", "occupancy_start_line_group");

        occupancyStartG
          .append("text")
          .text("occ. start")
          .style("fill", "#dedede");

        occupancyStartG
          .append("line")
          .attr("x1", 0)
          .attr("y1", this.margin.top)
          .attr("y2", this.state.height - this.margin.bottom)
          .attr("x2", 0)
          .style("stroke-width", 1)
          .style("stroke", "#dedede");

        var occupancyEndG = occupancyG
          .append("g")
          .attr("class", "occupancy_end_line_group");

        occupancyEndG.append("text").text("occ. end").style("fill", "#dedede");

        occupancyEndG
          .append("line")
          .attr("x1", 0)
          .attr("y1", this.margin.top)
          .attr("y2", this.state.height - this.margin.bottom)
          .attr("x2", 0)
          .style("stroke-width", 1)
          .style("stroke", "#dedede");
      }
    }

    // update positions
    var chart = this;
    this.occupancyLinesG
      .selectAll(".occupancy_line_group")
      .each(function (d, i) {
        var startSample = occupancies[i].startSample;
        var endSample = occupancies[i].endSample;
        var lIdx = chart.state.data.sampleNumberToIndexMap[startSample];
        var rIdx = chart.state.data.sampleNumberToIndexMap[endSample];

        var startTime =
          chart.state.data.formatted[0].realTimeIntervals[lIdx][0];
        var endTime = chart.state.data.formatted[0].realTimeIntervals[rIdx][1];

        var startLine = d3.select(this).select(".occupancy_start_line_group");
        var endLine = d3.select(this).select(".occupancy_end_line_group");

        startLine
          .attr("y1", chart.margin.top)
          .attr("y2", chart.state.height - chart.margin.bottom);

        startLine
          .select("text")
          .attr(
            "transform",
            "translate(2," +
              (chart.margin.top +
                startLine.select("text").node().getBBox().height) +
              ")"
          );
        startLine.attr("transform", "translate(" + xScale(startTime) + ",0)");

        endLine
          .select("line")
          .attr("y1", chart.margin.top)
          .attr("y2", chart.state.height - chart.margin.bottom);

        endLine
          .select("text")
          .attr(
            "transform",
            "translate(2," +
              (chart.margin.top +
                20 +
                endLine.select("text").node().getBBox().height) +
              ")"
          );
        endLine.attr("transform", "translate(" + xScale(endTime) + ",0)");

        // if end is the final data point, then hide the line to avoid overlap with axis line
        if (
          rIdx ===
          chart.state.data.formatted[0].realTimeIntervals.length - 1
        ) {
          endLine.attr("visibility", "hidden");
        }
      });
  } else {
    // if (this.state.data.raw.occupancies)
    this.occupancyLinesG.selectAll(".occupancy_line_group").remove();
  } // if (this.state.data.raw.occupancies)

  // format minor axis labels x-axis
  axisBottomTicks.each(function (d, i) {
    var tickText = d3.select(this).select("text");
    var tickLine = d3.select(this).select("line");
    if (!xTicksGenerated[i].isMajor) {
      tickText.attr("visibility", "hidden");
      tickLine.attr("y2", 4);
    }
  });

  // add grid lines
  if (this.options.gridx) {
    this.verticalGridG
      .attr(
        "transform",
        "translate(0," + (this.state.height - this.margin.bottom) + ")"
      )
      .attr("id", "th_x-grid")
      .call(xAxis.tickFormat("").tickSize(-plotHeight));

    var verticalGridLines = this.verticalGridG.selectAll("g.tick");
    verticalGridLines.each(function (d, i) {
      var tickLine = d3.select(this).select("line");
      if (!xTicksGenerated[i].isMajor) {
        tickLine.attr("stroke", "#f6f6f6");
      } else {
        tickLine.attr("stroke", "#d3d3d3");
      }
    });
  } else {
    this.verticalGridG.selectAll("*").remove();
  }

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

    // add grid lines
    if (this.options.gridy) {
      this.horizontalLeftGridG
        .attr("transform", "translate(" + this.margin.left + ",0)")
        .attr("id", "th_y1-grid")
        .call(yAxisLeft.tickFormat("").tickSize(-plotWidth));

      var horizontalLeftGridLines =
        this.horizontalLeftGridG.selectAll("g.tick");
      horizontalLeftGridLines.each(function (d, i) {
        var tickLine = d3.select(this).select("line");
        tickLine.attr("stroke", "#d3d3d3");
      });
    } else {
      this.horizontalLeftGridG.selectAll("*").remove();
    }

    var axisLabelY1 = this.svg.select("#th_label_y1");

    var axisLabelY1Text =
      compressionIndex == 0
        ? this.options.y1title
        : this.options.y1title +
          " per " +
          Math.pow(2, compressionIndex) +
          " Samples";

    // if already drawn, just update
    if (!axisLabelY1.empty()) {
      // reposition
      axisLabelY1
        .attr(
          "transform",
          "translate(" +
            (this.margin.left - this.axisLeftG.node().getBBox().width - 5) +
            "," +
            this.state.height / 2 +
            ") rotate(-90)"
        )
        .text(axisLabelY1Text);
    } else {
      // create new
      this.svg
        .append("text")
        .attr("class", "yaxistitle")
        .attr("id", "th_label_y1")
        .attr(
          "transform",
          "translate(" +
            (this.margin.left - this.axisLeftG.node().getBBox().width - 5) +
            "," +
            this.state.height / 2 +
            ") rotate(-90)"
        )
        .style("text-anchor", "middle")
        .text(axisLabelY1Text)
        .attr("font-size", "0.9em");
    } // if (!axisLabelY1.empty())

    // format minor axis labels:
  } // if (HAS_GAMMA)

  if (HAS_NEUTRON) {
    var yAxisRight = d3.svg
      .axis()
      .scale(yScaleNeutron)
      .ticks(3)
      .orient("right");

    // create or update axis
    this.axisRightG
      .attr(
        "transform",
        "translate(" + (this.state.width - this.margin.left) + ",0)"
      )
      .call(yAxisRight);

    // add grid lines
    if (this.options.gridy) {
      this.horizontalRightGridG
        .attr(
          "transform",
          "translate(" + (this.state.width - this.margin.left) + ",0)"
        )
        .attr("id", "th_y1-grid")
        .call(yAxisRight.tickFormat("").tickSize(-plotWidth));

      var horizontalRightGridLines =
        this.horizontalRightGridG.selectAll("g.tick");
      horizontalRightGridLines.each(function (d, i) {
        var tickLine = d3.select(this).select("line");
        tickLine.attr("stroke", "#f6f6f6");
      });
    } else {
      this.horizontalRightGridG.selectAll("*").remove();
    }

    var axisLabelY2 = this.svg.select("#th_label_y2");
    var axisLabelY2Text =
      compressionIndex == 0
        ? this.options.y2title
        : this.options.y2title +
          " per " +
          Math.pow(2, compressionIndex) +
          " Samples";

    // if already drawn, just update
    if (!axisLabelY2.empty()) {
      // reposition
      axisLabelY2
        .attr(
          "transform",
          "translate(" +
            (this.state.width -
              this.margin.left +
              this.axisRightG.node().getBBox().width +
              10) +
            "," +
            this.state.height / 2 +
            ") rotate(90)"
        )
        .text(axisLabelY2Text);
    } else {
      // create new
      this.svg
        .append("text")
        .attr("class", "yaxistitle")
        .attr("id", "th_label_y2")
        .attr(
          "transform",
          "translate(" +
            (this.state.width -
              this.margin.left +
              this.axisRightG.node().getBBox().width +
              10) +
            "," +
            this.state.height / 2 +
            ") rotate(90)"
        )
        .style("text-anchor", "middle")
        .text(axisLabelY2Text);
    } // if (!axisLabelY2.empty())
  } else {
    this.axisRightG.selectAll("*").remove();
    this.svg.select("#th_label_y2").remove();
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
 *                  },
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
 * @param {Object} rawData: data Object passed from Wt
 * @returns a boolean
 */
D3TimeChart.prototype.isValidRawData = function (rawData) {
  if (!rawData) {
    return false;
  }

  // Check properties
  if (
    !rawData.hasOwnProperty("sampleNumbers") ||
    !rawData.hasOwnProperty("realTimes") ||
    !rawData.hasOwnProperty("gammaCounts")
  ) {
    return false;
  }

  // Check data types
  if (
    !Array.isArray(rawData.sampleNumbers) ||
    !Array.isArray(rawData.realTimes) ||
    !Array.isArray(rawData.gammaCounts)
  ) {
    return false;
  }

  var nSamples = rawData.sampleNumbers.length;

  // Check matching lengths
  if (rawData.realTimes.length !== nSamples) {
    return false;
  }

  // check data of each gamma detector
  rawData.gammaCounts.forEach(function (detector) {
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
  if (rawData.hasOwnProperty("neutronCounts")) {
    if (!Array.isArray(rawData.neutronCounts)) {
      return false;
    }

    rawData.neutronCounts.forEach(function (detector) {
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
  } // if (rawData.hasOwnProperty("neutronCounts"))

  return true;
};

/**
 * Formats data to several JSON objects for convenient access and plotting. Data input assumed to be in same format as the raw data sent from Wt
 * returns object bundling sampleNumbers, realTimeIntervals, mean time of an interval, gps coordinates, occupancy status of intervals, source types of intervals, the timestamps at the start of each interval, data domain, and detectors data.
 * e.g.
 * {
 *   sampleNumbers: [1, 2, 3, ...],
 *   realTimeIntervals: [[a,b], [c,d], ...],
 *   meanIntervalTime: 0.500...
 *   gpsCoordinates: [...], (optional)
 *   occupied: [false, true, true, ...]
 *   sourceTypes: [2, 3, 3, ..]
 *   startTimeStamps: [...]
 *   domains:
 *   {
 *      x: [minTime, maxTime],
 *      yGamma: [minGamma, maxGamma]
 *      yNeutron: [minNeutron, yNeutron]
 *   }
 *
 *   detectors:
 *   {
 *      det1: {
 *                meta: {
 *                          "isGammaDetector": true,
 *                          "isNeutronDetector": true,
 *                          "gammaColor": rgb(0,0,0),
 *                          "neutronColor": rgb(0,128,0)
 *                      },
 *                counts: [{DataPoint}, {DataPoint}, {DataPoint}, ...],
 *                gammaLiveTimes: [481.999, 2.00028, ...] // livetimes used for computing CPS if available
 *                neutronLiveTimes: [481.671, 1.99889, ...]
 *            },
 *
 *      det2: {...},
 *      ...
 *   }
 * }
 *
 * @param {Object} rawData : data in same format as the raw data sent from Wt
 */
D3TimeChart.prototype.formatDataFromRaw = function (rawData) {
  var nSamples = rawData.sampleNumbers.length;
  var detectors = {};

  // get array of realTimes for intervals
  var realTimeIntervals = this.getRealTimeIntervals(
    rawData.realTimes,
    rawData.sourceTypes
  );

  // get average interval time
  var meanIntervalTime = this.getMeanIntervalTime(
    rawData.realTimes,
    rawData.sourceTypes
  );

  // format occupied array:
  var occupied;
  if (rawData.hasOwnProperty("occupancies")) {
    var occupied = [];
    for (var i = 0; i < nSamples; i++) {
      var sampleNumber =
        typeof rawData.sampleNumbers[i] === "number" &&
        isFinite(rawData.sampleNumbers[i])
          ? rawData.sampleNumbers[i]
          : parseInt(Object.keys(rawData.sampleNumbers[i])[0]);
      occupied[i] = this.isOccupiedSample(sampleNumber, rawData.occupancies);
    }
  }

  // format startTimeStamps array:
  var startTimeStamps;
  if (
    rawData.hasOwnProperty("startTimeOffset") &&
    rawData.hasOwnProperty("startTimes")
  ) {
    startTimeStamps = rawData.startTimes.map(function (startTime) {
      if (startTime == null) {
        return null;
      }
      return rawData.startTimeOffset + startTime;
    });
  }

  // get domains of the data (range of x values, range of yGamma values, range of yNeutron values)
  var domains = this.getDomainsFromRaw(rawData);

  // format detectors data (detector metadata, livetimes, and counts)
  rawData.gammaCounts.forEach(function (det) {
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
    if (!detectors[det.detName].hasOwnProperty("counts")) {
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

        // if cps > 0 ever, then is a gamma
        if (!detectors[det.detName].meta.isGammaDetector && cps > 0) {
          detectors[det.detName].meta.isGammaDetector = true;
        }
      }
    }
  }); // rawData.gammaCounts.forEach

  var HAS_NEUTRON = false;

  if (rawData.hasOwnProperty("neutronCounts")) {
    rawData.neutronCounts.forEach(function (det) {
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

          if (!detectors[det.detName].meta.isNeutronDetector && cps > 0) {
            detectors[det.detName].meta.isNeutronDetector = true;
            HAS_NEUTRON = true;
          }
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

          // if cps > 0 ever, then is a neutrondetector
          if (!detectors[det.detName].meta.isNeutronDetector && cps > 0) {
            detectors[det.detName].meta.isNeutronDetector = true;
            HAS_NEUTRON = true;
          }
        }
      } // if (!detectors[det.detName].hasOwnProperty("counts"))
    }); // rawData.neutronCounts.forEach
  } // if (rawData.hasOwnProperty("neutronCounts"))

  if (!HAS_NEUTRON) {
    // set margin
    this.margin.right = 10;
  } else {
    this.margin.right = 60;
  }

  return {
    sampleNumbers: rawData.sampleNumbers,
    realTimeIntervals: realTimeIntervals,
    meanIntervalTime: meanIntervalTime,
    occupied: occupied,
    sourceTypes: rawData.sourceTypes,
    startTimeStamps: startTimeStamps,
    gpsCoordinates: rawData.gpsCoordinates,

    detectors: detectors,
    domains: domains,
  };
};

/**
 * Gets data domains. Data assumed to be in same format as the raw data sent from Wt
 * @param {Object} rawData: data in same format as the raw data sent from Wt
 */
D3TimeChart.prototype.getDomainsFromRaw = function (rawData) {
  var realTimeIntervals = this.getRealTimeIntervals(
    rawData.realTimes,
    rawData.sourceTypes
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
  var nSamples = rawData.sampleNumbers.length;

  for (var i = 0; i < rawData.gammaCounts.length; i++) {
    for (var j = 0; j < nSamples; j++) {
      var cps = rawData.gammaCounts[i].liveTimes
        ? rawData.gammaCounts[i].counts[j] / rawData.gammaCounts[i].liveTimes[j]
        : rawData.gammaCounts[i].counts[j] / rawData.realTimes[j];

      if (cps > yMaxGamma) {
        yMaxGamma = cps;
      }
    }
  }

  // Find max over neutron detectors
  if (rawData.hasOwnProperty("neutronCounts")) {
    for (var i = 0; i < rawData.neutronCounts.length; i++) {
      for (var j = 0; j < nSamples; j++) {
        var cps = rawData.neutronCounts[i].liveTimes
          ? rawData.neutronCounts[i].counts[j] /
            rawData.neutronCounts[i].liveTimes[j]
          : rawData.neutronCounts[i].counts[j] / rawData.realTimes[j];

        if (cps > yMaxNeutron) {
          yMaxNeutron = cps;
        }
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
 * Get D3 linear scaling functions based on the data domain and chart element dimensions (height, width, margin).
 * this.state.data.formatted, this.state.height, and this.state.width must be defined for D3TimeChart
 *
 * @param {Object} domains: domains object, must be of form:
 * {
 *    x: [a, b],
 *    yGamma: [c, d],
 *    yNeutron: [e, f],
 * }
 */
D3TimeChart.prototype.getScales = function (domains) {
  var xScale = domains.x
    ? d3.scale
        .linear()
        .domain(domains.x)
        .range([this.margin.left, this.state.width - this.margin.right])
    : undefined;
  var yScaleGamma = domains.yGamma
    ? d3.scale
        .linear()
        .domain(domains.yGamma)
        .range([this.state.height - this.margin.bottom, this.margin.top])
    : undefined;
  var yScaleNeutron = domains.yNeutron
    ? d3.scale
        .linear()
        .domain(domains.yNeutron)
        .range([this.state.height - this.margin.bottom, this.margin.top])
    : undefined;

  return {
    xScale: xScale,
    yScaleGamma: yScaleGamma,
    yScaleNeutron: yScaleNeutron,
  };
};

/**
 * Computes array of sequential real time intervals for each data point from raw time segments, AND record background duration if there is a background
 * @param {Number[]} realTimes: realTimes array passed from Wt
 * @param {Number[]} sourceTypes: sourceTypes array passed from Wt
 * @returns an array of length-two arrays which represent time intervals for individual samples.
 */
D3TimeChart.prototype.getRealTimeIntervals = function (realTimes, sourceTypes) {
  var realTimeIntervals = [];

  for (var i = 0; i < realTimes.length; i++) {
    if (sourceTypes && i === 0 && sourceTypes[i] === 2) {
      // center so background is not started at 0
      // to reduce long lead-in for plotting, replace with a smaller value inside realTimeIntervals and then record the true backgroundDuration in this.state for future reference:
      var leadTime = Math.max(
        -realTimes[i],
        -Math.floor(realTimes.length * 0.12)
      );
      this.state.data.backgroundDuration = realTimes[i];
      realTimeIntervals[i] = [leadTime, 0];
    } else if (i === 0) {
      // first point is not background, so don't need to center
      this.state.data.backgroundDuration = null;
      realTimeIntervals[i] = [0, realTimes[i]];
    } else {
      var prevEndpoint = realTimeIntervals[i - 1][1];
      realTimeIntervals[i] = [prevEndpoint, prevEndpoint + realTimes[i]];
    }
  }

  return realTimeIntervals;
};

/**
 * Computes the mean interval time over all foreground sourcetype data
 * @param {Number[]} realTimes : realTimes array of raw data
 * @param {Number[]} sourceTypes : sourceTypes array of raw data
 * @returns a number which is the average interval time length over all foreground sourcetype data
 */
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
 * Function to handle shifting of the zoomed selection window.
 * @param {Number} n : shift amount
 */
D3TimeChart.prototype.shiftSelection = function (n) {
  if (this.state.selection) {
    var domain = this.state.selection.domain;
    var compressionIndex = this.state.selection.compressionIndex;

    var stepSize =
      this.state.data.formatted[compressionIndex].meanIntervalTime *
      n *
      Math.pow(2, compressionIndex);

    // compute shifted bound
    var rightBound = this.state.data.formatted[0].domains.x[1];
    var leftBound = this.state.data.formatted[0].domains.x[0];

    if (domain[1] + stepSize > rightBound) {
      var shiftAmount = rightBound - domain[1];
      domain[1] = rightBound;
      domain[0] += shiftAmount;
    } else if (domain[0] + stepSize < leftBound) {
      var shiftAmount = leftBound - domain[0];
      domain[0] = leftBound;
      domain[1] += shiftAmount;
    } else {
      domain[0] += stepSize;
      domain[1] += stepSize;
    }

    // update selection, update brush scale, update chart
    this.state.selection.domain = domain;

    var fullDomain = {
      x: domain,
      yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
      yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
    };
    var scale = this.getScales(fullDomain);
    this.state.brush.setScale(scale.xScale);
    this.updateChart(scale, compressionIndex);
  }
  ``;
};

/**
 * Function to handle mouse wheel for zooming and panning.
 * @param {Number} deltaX : integer deltaX value of the mouse scroll event
 * @param {Number} deltaY : integer deltaY value of the mouse scroll event
 * @param {Number} mouseX : integer x-coordinates of pointer in pixels relative to the containing element
 */
D3TimeChart.prototype.handleMouseWheel = function (deltaX, deltaY, mouseX) {
  var brush = this.state.brush;
  var xScale = brush.getScale();
  var focalPoint = xScale.invert(mouseX);

  var leftLimit =
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[0];
  var rightLimit =
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[1];

  var currentDomain = this.state.selection
    ? this.state.selection.domain
    : this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .domains.x;

  var minimumMeanIntervalTime = this.state.data.formatted[0].meanIntervalTime;

  // don't allow zoom in more than 2 mean interval lengths or zoom out if current domain is all the way zoomed out already.
  if (
    currentDomain[1] - currentDomain[0] <= minimumMeanIntervalTime * 2 &&
    deltaY < 0
  ) {
    return;
  } else if (
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[0] === currentDomain[0] &&
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[1] === currentDomain[1] &&
    deltaY > 0
  ) {
    return;
  }

  var zoomStepSize =
    0.001 * Math.exp(2, this.state.data.unzoomedCompressionIndex);

  let newLeftExtent;
  let newRightExtent;

  if (Math.abs(deltaY) >= Math.abs(deltaX)) {
    // if scroll vertical, zoom
    newLeftExtent = Math.max(
      currentDomain[0] -
        deltaY * Math.abs(currentDomain[0] - focalPoint) * zoomStepSize,
      leftLimit
    );
    newRightExtent = Math.min(
      currentDomain[1] +
        deltaY * Math.abs(currentDomain[1] - focalPoint) * zoomStepSize,
      rightLimit
    );
  } else {
    // if scroll horizontal, pan
    let panStepSize =
      deltaX *
      0.1 *
      this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .meanIntervalTime;

    if (currentDomain[0] + panStepSize < leftLimit) {
      panStepSize = leftLimit - currentDomain[0];
    } else if (currentDomain[1] + panStepSize > rightLimit) {
      panStepSize = rightLimit - currentDomain[1];
    }

    newLeftExtent = currentDomain[0] + panStepSize;
    newRightExtent = currentDomain[1] + panStepSize;
  }

  // if the window is smaller than 2x the mean interval time, then
  if (newRightExtent - newLeftExtent < 2 * minimumMeanIntervalTime) {
    var extentCenter = (newRightExtent + newLeftExtent) / 2;
    newLeftExtent = Math.max(extentCenter - minimumMeanIntervalTime, leftLimit);
    newRightExtent = Math.min(
      extentCenter + minimumMeanIntervalTime,
      rightLimit
    );
  }

  // compute new compression index to use
  var leftIndex = this.findDataIndex(newLeftExtent, 0);
  var rightIndex = this.findDataIndex(newRightExtent, 0);
  var nPoints = rightIndex - leftIndex + 1;
  var plotWidth = this.state.width - this.margin.left - this.margin.right;

  var compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));

  // obtain new scales, update brush, update selection
  var scales = this.getScales({
    x: [newLeftExtent, newRightExtent],
    yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
    yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
  });

  brush.setScale(scales.xScale);

  if (
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[0] === newLeftExtent &&
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[1] === newRightExtent
  ) {
    // completely zoomed out, so set selection to null
    this.state.selection = null;
  } else {
    this.state.selection = {
      domain: [newLeftExtent, newRightExtent],
      compressionIndex: compressionIndex,
    };
  }

  // update chart
  this.updateChart(scales, compressionIndex);
};

/**
 * Brush handler to zoom into a selected domain or zoom out if brushed back.
 * this.state.brush must be not null, and start and end properties must be set.
 */
D3TimeChart.prototype.handleBrushZoom = function () {
  var brush = this.state.brush;
  if (!brush || !brush.getStart() || !brush.getEnd()) {
    throw new Error("Brush start or end are undefined!");
  }
  // handle dragback (zoom out). Must update selection and brush scale upon mouseup event.
  // For dragback, handleDragBackZoom already handles chart updating--handleBrushZoom only updates selection and mouseup
  if (!this.draggedForward) {
    // if no selection (not already zoomed in), dragback does nothing.
    if (!this.state.selection) {
      return;
    }

    // calculate new extent if it is same as extent when all the way zoomed out, set selection to null.
    // naturalXScale is the xScale used when all the way zoomed out (using all data points).
    var naturalScale = this.getScales(
      this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .domains
    );
    var naturalXScale = naturalScale.xScale;

    var zoomOutAmount = -(
      naturalXScale.invert(brush.getScale()(brush.getEnd())) -
      naturalXScale.invert(brush.getScale()(brush.getStart()))
    );

    var newLeftExtent = Math.max(
      this.state.selection.domain[0] - zoomOutAmount,
      this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .domains.x[0]
    );
    var newRightExtent = Math.min(
      this.state.selection.domain[1] + zoomOutAmount,
      this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
        .domains.x[1]
    );
    var newExtent = [newLeftExtent, newRightExtent];

    // if extent is identical to the natural x-domain (i.e. when all the way zoomed out), then set selection to null. And set brush scale to the natural x-scale.
    if (
      newLeftExtent ===
        this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
          .domains.x[0] &&
      newRightExtent ===
        this.state.data.formatted[this.state.data.unzoomedCompressionIndex]
          .domains.x[1]
    ) {
      this.state.selection = null;
      // update brush scale
      brush.setScale(naturalXScale);
      this.updateChart(naturalScale, this.state.data.unzoomedCompressionIndex);
    } else {
      // Otherwise, update the selection to reflect the new domain and compression level, and update the scale of the brush.
      // compute new compression index to use
      var leftIndex = this.findDataIndex(newExtent[0], 0);
      var rightIndex = this.findDataIndex(newExtent[1], 0);
      var nPoints = rightIndex - leftIndex + 1;
      var plotWidth = this.state.width - this.margin.left - this.margin.right;

      var compressionIndex = Math.ceil(
        Math.log2(Math.ceil(nPoints / plotWidth))
      );

      // calculate new scale. Update brush
      var scales = this.getScales({
        x: newExtent,
        yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
        yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
      });
      this.state.selection = {
        domain: newExtent,
        compressionIndex: compressionIndex,
      };

      brush.setScale(scales.xScale);
      this.updateChart(scales, compressionIndex);
    } // if (newLeftExtent === this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains.x[0] && newRightExtent === this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains.x[1])
  } else {
    // handle drag forward (zoom-in)

    // if drawn window is too small, do nothing
    var width =
      brush.getScale()(brush.getEnd()) - brush.getScale()(brush.getStart());
    if (width < this.options.minSelectionWidth) {
      return;
    }

    // find appropriate compression level for this range
    var leftIndex = this.findDataIndex(brush.extent()[0], 0);
    var rightIndex = this.findDataIndex(brush.extent()[1], 0);
    var nPoints = rightIndex - leftIndex + 1;
    var plotWidth = this.state.width - this.margin.left - this.margin.right;

    var compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));

    // set lower limit on extent size to 2 interval lengths
    var minExtentLeft = Math.max(
      brush.getCenter() -
        this.state.data.formatted[compressionIndex].meanIntervalTime,
      this.state.data.formatted[compressionIndex].domains.x[0]
    );
    var maxExtentRight = Math.min(
      brush.getCenter() +
        this.state.data.formatted[compressionIndex].meanIntervalTime,
      this.state.data.formatted[compressionIndex].domains.x[1]
    );

    var brushWidth = brush.extent()[1] - brush.extent()[0];

    var extent =
      brushWidth >
      2 * this.state.data.formatted[compressionIndex].meanIntervalTime
        ? brush.extent()
        : [minExtentLeft, maxExtentRight];

    // update x-domain to the new domain

    // update selection
    this.state.selection = {
      domain: extent,
      compressionIndex: compressionIndex,
    };

    var scales = this.getScales({
      x: this.state.selection.domain,
      yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
      yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
    });

    // update brush scale
    brush.setScale(scales.xScale);
    this.updateChart(scales, compressionIndex);
  }
};

/**
 * Handler for drawing rectangle on brush forward.
 * this.state.brush must be not null, and start and end properties must be set.
 */
D3TimeChart.prototype.handleDragForwardZoom = function () {
  var brush = this.state.brush;
  if (!brush || !brush.getStart() || !brush.getEnd()) {
    throw new Error("Brush start or end are undefined!");
  }

  // if has been dragged backward prior to dragging forward, should set the flag
  if (!this.draggedForward) {
    // Set draggedforward flag. Flag is used to find conditions for redrawing chart at default position (prior to any drag-back zoomout occurring)
    this.draggedForward = true;

    //re-render the chart at default position if needed (prior to any drag-back zoomout occuring)
    if (this.state.selection) {
      var scales = this.getScales({
        x: this.state.selection.domain,
        yGamma:
          this.state.data.formatted[this.state.selection.compressionIndex]
            .domains.yGamma,
        yNeutron:
          this.state.data.formatted[this.state.selection.compressionIndex]
            .domains.yNeutron,
      });

      this.updateChart(scales, this.state.selection.compressionIndex);
    }
  }

  // compute width of drawn rectangle.
  var width =
    brush.getScale()(brush.getEnd()) - brush.getScale()(brush.getStart());

  this.mouseMoveHighlight(width);
};

/**
 * Handler for updating chart to zoom out on brush backward.
 * this.state.brush must be not null, and start and end properties must be set.
 */
D3TimeChart.prototype.handleDragBackZoom = function () {
  var brush = this.state.brush;
  if (!brush || !brush.getStart() || !brush.getEnd()) {
    throw new Error("Brush start or end are undefined!");
  }

  // clear rectangle if any drawn
  this.highlightRect.attr("width", 0);
  this.highlightText.attr("visibility", "hidden");

  // un-set draggedforward flag. Flag is used to find conditions for redrawing chart at default position (prior to any drag-back zoomout occurring)
  this.draggedForward = false;

  // if not zoomed in already (no selection), dragback does nothing.
  if (!this.state.selection) {
    return;
  }

  var naturalXScale = this.getScales(
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
  ).xScale;
  var zoomOutAmount = -(
    naturalXScale.invert(brush.getScale()(brush.getEnd())) -
    naturalXScale.invert(brush.getScale()(brush.getStart()))
  );

  // compute new extent
  var newLeftExtent = Math.max(
    this.state.selection.domain[0] - zoomOutAmount,
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[0]
  );
  var newRightExtent = Math.min(
    this.state.selection.domain[1] + zoomOutAmount,
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[1]
  );
  var newExtent = [newLeftExtent, newRightExtent];

  // compute new compression index to use
  var leftDataIndex = this.findDataIndex(newExtent[0], 0);
  var rightDataIndex = this.findDataIndex(newExtent[1], 0);
  var nPoints = rightDataIndex - leftDataIndex + 1;
  var plotWidth = this.state.width - this.margin.left - this.margin.right;

  var compressionIndex = Math.ceil(Math.log2(Math.ceil(nPoints / plotWidth)));

  var scales = this.getScales({
    x: newExtent,
    yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
    yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
  });

  this.updateChart(scales, compressionIndex);
};

/**
 * Sets initial attributes of a highlight region
 * @param {Number} mouseX : x-coordinates of pointer in pixels relative to the containing element
 */
D3TimeChart.prototype.mouseDownHighlight = function (mouseX, modifier) {
  if (this.options.useSimplifiedGestures) {
    this.highlightRect.attr("fill", this.HIGHLIGHT_COLORS.zoom);
    this.highlightText.text("Zoom in");
  } else {
    var foreground = this.highlightOptions.foreground;
    var background = this.highlightOptions.background;
    var zoom = this.highlightOptions.zoom;

    if (foreground && modifier in foreground.modifierKey) {
      this.highlightRect.attr("fill", this.HIGHLIGHT_COLORS.foreground);
      this.highlightText.text("Select foreground");
    } else if (background && modifier in background.modifierKey) {
      this.highlightRect.attr("fill", this.HIGHLIGHT_COLORS.background);
      this.highlightText.text("Select background");
    } else if (zoom && modifier in zoom.modifierKey) {
      this.highlightRect.attr("fill", this.HIGHLIGHT_COLORS.zoom);
      this.highlightText.text("Zoom in");
    }
  }

  this.highlightRect
    .attr("height", this.state.height - this.margin.top - this.margin.bottom)
    .attr("x", mouseX)
    .attr("y", this.margin.top)
    .attr("width", 1);

  this.highlightText
    .attr("x", mouseX - this.highlightText.node().getBBox().width / 2)
    .attr(
      "y",
      (this.state.height - this.margin.top - this.margin.bottom) / 2 +
        this.margin.top
    );
};

/**
 * Updates width of highlight region on mouse move.
 * @param {Number} width : width in pixels
 */
D3TimeChart.prototype.mouseMoveHighlight = function (width) {
  if (width > 0) {
    this.highlightRect.attr("width", width);
    if (width < this.options.minSelectionWidth) {
      this.highlightText.style("visibility", "hidden");
    } else {
      this.highlightText
        .attr("transform", `translate(${width / 2}, 0)`)
        .style("visibility", "visible");
    }
  } else {
    this.highlightRect.attr("width", 0);
  }
};

/**
 * Handles clearing of highlight rectangles on mouse-up.
 */
D3TimeChart.prototype.mouseUpHighlight = function () {
  this.highlightRect.attr("height", 0).attr("width", 0);
  this.highlightText.style("visibility", "hidden");
};

/**
 * Handles panning selection window by dragging on bottom axis.
 * this.state.brush must be not null, and start and end properties must be set.
 * @returns an object containing the new selection and new scale of the panned selection. The caller of this function must set the selection and brush scale outside on mouse-up event.
 */
D3TimeChart.prototype.handleBrushPanSelection = function () {
  // if not zoomed in already (no selection), drag does nothing.
  if (!this.state.selection) {
    return;
  }

  var brush = this.state.brush;

  if (!brush || !brush.getStart() || !brush.getEnd()) {
    throw new Error("Brush start or end are undefined!");
  }

  var domain = this.state.selection.domain;
  var compressionIndex = this.state.selection.compressionIndex;

  // compute new extent
  var panAmount =
    brush.getScale().invert(brush.getScale()(brush.getStart())) -
    brush.getScale().invert(brush.getScale()(brush.getEnd()));
  var rightBound =
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[1];
  var leftBound =
    this.state.data.formatted[this.state.data.unzoomedCompressionIndex].domains
      .x[0];

  // if will take you out of bound, do return and do nothing
  if (domain[0] + panAmount < leftBound || domain[1] + panAmount > rightBound) {
    return;
  }

  // compute new extent
  var newLeftExtent = domain[0] + panAmount;
  var newRightExtent = domain[1] + panAmount;
  var newExtent = [newLeftExtent, newRightExtent];

  // update chart
  var scales = this.getScales({
    x: newExtent,
    yGamma: this.state.data.formatted[compressionIndex].domains.yGamma,
    yNeutron: this.state.data.formatted[compressionIndex].domains.yNeutron,
  });

  this.updateChart(scales, compressionIndex);

  // return new selection and new scale so can later update when drag ends
  return {
    newSelection: {
      domain: newExtent,
      compressionIndex: compressionIndex,
    },
    newScale: scales,
  };
};

/**
 * Handles showing tooltip.
 */
D3TimeChart.prototype.showToolTip = function () {
  this.hoverToolTip.style("visibility", "visible");
};

/**
 * Handler to update tooltip display.
 * @param {Number} time : real time of measurement (x-axis)
 * @param {Object[]} data : Array of data objects, where each object at minimum has fields: {detName, gammaCPS}. Optional fields: neutronCPS
 * @param {Object} optargs : Object of optional keyword arguments. Accepted properties include: startTimeStamp, sourceType
 */
D3TimeChart.prototype.updateToolTip = function (time, data, optargs) {
  this.hoverToolTip.html(this.createToolTipString(time, data, optargs));
};

/**
 * Handler to hide tooltip.
 */
D3TimeChart.prototype.hideToolTip = function () {
  this.hoverToolTip.style("visibility", "hidden");
};

/**
 * Handler to create tooltip string from data.
 * @param {Number} time : real time of measurement (x-axis)
 * @param {Object[]} data : Array of data objects, where each object at minimum has fields: {detName, gammaCPS}. Optional fields: neutronCPS
 * @param {Object} optargs : Object of optional keyword arguments. Accepted properties include: startTimeStamp, sourceType
 * @returns the tooltip string
 */
D3TimeChart.prototype.createToolTipString = function (time, data, optargs) {
  var s =
    optargs.startTimeStamp != null
      ? "<div>" + new Date(optargs.startTimeStamp).toLocaleString() + "</div>"
      : "";

  // If want compression data in the tooltip, uncomment below
  // var compressionIndex = this.state.selection
  // ? this.state.selection.compressionIndex
  // : this.state.data.unzoomedCompressionIndex;

  // s += "<div>Data Compression Level: " + compressionIndex + "</div>";

  // // If want sourcetype data in the tooltip, uncomment below
  // s +=
  //   optargs.sourceType != null
  //     ? "<div>Source: " + this.SOURCE_MAP[optargs.sourceType] + "</div>"
  //     : "";

  s += "<div>Time: " + time.toPrecision(4) + " s</div>";

  // for each detector, give counts
  for (var i = 0; i < data.length; i++) {
    if (typeof data[i].sampleNumber === "number") {
      s += "<div>Sample Num: " + data[i].sampleNumber.toString();
    } else {
      var sampleNumbers = Object.keys(data[i].sampleNumber);
      if (sampleNumbers.length > 4) {
        s +=
          "<div>Sample Nums: " +
          sampleNumbers[0].toString() +
          "..." +
          sampleNumbers[sampleNumbers.length - 1].toString();
      } else {
        s += "<div>Sample Nums: " + sampleNumbers.toString();
      }
    }
    if (data[i].detName.length > 0) {
      s += " (" + data[i].detName + ")</div>";
    }

    s += "<div>G CPS: " + data[i].gammaCPS.toPrecision(6);

    if (data[i].detName.length > 0) {
      s += " (" + data[i].detName + ")</div>";
    }

    if (data[i].neutronCPS != null) {
      // cps of 0 is still valid to display
      s += "<div>N CPS: " + data[i].neutronCPS.toPrecision(3);

      if (data[i].detName.length > 0) {
        s += " (" + data[i].detName + ")</div>";
      }
    }
  }

  return s;
};

/**
 * Performs binary search over data to find interval that the given time belongs to
 * Returns the index of the match in the data array
 * @param {Number} time : realTime of measurement (x-value)
 * @param {Number} compressionIndex : optional compression index to search over. Default (no argument) is the current compressionIndex based on zoom amount.
 */
D3TimeChart.prototype.findDataIndex = function (time, compressionIndex) {
  var cIdx =
    compressionIndex != null
      ? compressionIndex
      : this.state.data.unzoomedCompressionIndex;
  var highIdx = this.state.data.formatted[cIdx].realTimeIntervals.length - 1;
  var lowIdx = 0;

  while (lowIdx <= highIdx) {
    var midIdx = Math.floor((highIdx + lowIdx) / 2);
    var interval = this.state.data.formatted[cIdx].realTimeIntervals[midIdx];
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

/**
 * Checks whether a given sample corresponds to a positive vehicle occupancy status.
 * @param {Number} sampleNumber : the integer sample number you are querying the occupancy status for
 * @param {Object[]} occupancies : raw occupancies data - an array of objects of form {color: [String], startSample: [Integer], endSample: [Integer]}
 */
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

/**
 * Compresses data by iteratively aggregating contiguous sets of n data points if they share the same attributes.
 * @param {Object} data : raw time-history data passed by Wt
 * @param {Number} n : a non-negative integer that specifies the compression level.
 * @returns the compressed data in the same JSON format as the raw (input) data passed by Wt.
 */
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

        // add to sampleNumbers. For ES5 compatibility, using object to mimic behavior of a set:
        if (out.sampleNumbers[outIdx] == null) {
          out.sampleNumbers[outIdx] = {};
        }
        out.sampleNumbers[outIdx][data.sampleNumbers[i + j]] = true;

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

/**
 * Highlights selected time intervals. Called by D3TimeChart::setHighlightRegionsToClient()
 * @param {Object[]} regions : array of objects. Objects are of form:
 * {
 *   startSample: int
 *   endSample: int
 *   fillColor: string
 *   type: string
 * }
 */
D3TimeChart.prototype.setHighlightRegions = function (regions) {
  if (
    !this.state.height ||
    !this.state.width ||
    this.state.height <= 0 ||
    this.state.width <= 0
  ) {
    return;
  }

  if (!regions || !Array.isArray(regions)) regions = [];

  if (
    this.state.data.formatted &&
    this.state.data.formatted.length &&
    this.state.data.sampleNumberToIndexMap
  ) {
    this.state.regions = regions;
    this.highlightRegionsG.selectAll("rect").remove();

    //See the c++ function D3TimeChart::setHighlightRegionsToClient() for format of data
    for (var i = 0; i < regions.length; i++) {
      // get index from sample number
      var startSample = regions[i].startSample;
      var endSample = regions[i].endSample;
      var fillColor = regions[i].fillColor;

      // Protect against invalid sample numbers specified in the regions
      if (
        !(startSample in this.state.data.sampleNumberToIndexMap) ||
        !(endSample in this.state.data.sampleNumberToIndexMap)
      )
        continue;

      var lIdx = this.state.data.sampleNumberToIndexMap[startSample];
      var rIdx = this.state.data.sampleNumberToIndexMap[endSample];

      if (
        !this.state.data.formatted[0].realTimeIntervals ||
        this.state.data.formatted[0].realTimeIntervals.length <= rIdx ||
        this.state.data.formatted[0].realTimeIntervals.length <= lIdx ||
        !Array.isArray(this.state.data.formatted[0].realTimeIntervals[lIdx]) ||
        !Array.isArray(this.state.data.formatted[0].realTimeIntervals[rIdx]) ||
        this.state.data.formatted[0].realTimeIntervals[lIdx].length < 2 ||
        this.state.data.formatted[0].realTimeIntervals[rIdx].length < 2
      ) {
        // don't draw anything
        continue;
      }

      // look up the corresponding time of the sample number using the index
      var startTime = this.state.data.formatted[0].realTimeIntervals[lIdx][0];
      var endTime = this.state.data.formatted[0].realTimeIntervals[rIdx][1];
      // draw a rectangle starting at the time and ending at the time with given height and fill color
      // console.log([startTime, endTime]);

      var scales = this.state.selection
        ? this.getScales({
            x: this.state.selection.domain,
            yGamma: this.state.data.formatted[0].domains.yGamma,
            yNeutron: this.state.data.formatted[0].domains.yNeutron,
          })
        : this.getScales(this.state.data.formatted[0].domains);
      var lPixel = scales.xScale(startTime);
      var rPixel = scales.xScale(endTime);

      var highlightWidth = rPixel - lPixel > 2 ? rPixel - lPixel : 2;
      this.highlightRegionsG
        .append("rect")
        .attr(
          "height",
          this.state.height - this.margin.top - this.margin.bottom
        )
        .attr("x", lPixel)
        .attr("y", this.margin.top)
        .attr("width", highlightWidth)
        .attr("fill", fillColor)
        .attr("fill-opacity", 0.5);
    }
  }
};

/**
 * Helper for generating major and minor chart tick values based on the initial values from the d3 automatic tick generator, which computes the tick values automatically from the scale used.
 * The reason we rely on the d3 automatic tick generator is because otherwise, very challenging to dynamically generate good tick values. d3 does a better job.
 * This function will divide the ticks generated from the d3 tick generator into minor divisions
 * @param {Number} start : starting tick value generated by d3 tick generator
 * @param {Number} step : step size between tick values generated by d3 tick generator
 * @param {Number} tickCount : number of ticks generated by d3 tick generator
 * @param {Number} minorDivisionsCount : number of times to divide original ticks by to form minor ticks. Must be even.
 * @returns an array of objects containing tick values and flag of whether the tick is major or minor
 */
D3TimeChart.prototype.generateTicks = function (
  start,
  step,
  tickCount,
  minorDivisionsCount
) {
  // must adjust to account for the space before and after the start and end tick values given by d3
  var tickArray = [];
  var adjustedTickCount = (tickCount - 1 + 2) * minorDivisionsCount;
  var adjustedStep = step / minorDivisionsCount;
  var adjustedStart = start - step;

  for (var i = 0; i < adjustedTickCount; i++) {
    // a hack to handle imprecisions of floating point. Keeps up to 5 decimal points if needed.
    var value = +(Math.round(adjustedStart + adjustedStep * i + "e+5") + "e-5");
    tickArray.push({
      value: value,
      isMajor: i % (minorDivisionsCount / 2) === 0,
    });
  }

  return tickArray;
};

// Unimplemented
D3TimeChart.prototype.setXAxisTitle = function (title) {
  this.options.xtitle = title;
  this.reinitializeChart();
  //redraw x-title
};

// Unimplemented
D3TimeChart.prototype.setY1AxisTitle = function (title) {
  this.options.y1title = title;
  this.reinitializeChart();
  //redraw y1-title (e.g., gamma CPS axis title)
};

// Unimplemented
D3TimeChart.prototype.setY2AxisTitle = function () {
  this.options.y2title = title;
  this.reinitializeChart();
  //redraw y2-title (e.g., neutron CPS axis title)
};

/**
 * Function to handle setting compact x axis. Called whenever toggle compact x-axis on/off. Triggers a resize.
 * @param {Boolean} compact : boolean value defining whether or not to use compact x-axis
 */
D3TimeChart.prototype.setCompactXAxis = function (compact) {
  //Make x-zis title comapact or not
  this.options.compactXAxis = compact;
  this.reinitializeChart();
};

D3TimeChart.prototype.setGridX = function (show) {
  this.options.gridx = show;
  this.reinitializeChart();
  //add/remove horizantal grid lines
};

D3TimeChart.prototype.setGridY = function (show) {
  this.options.gridy = show;
  this.reinitializeChart();
  //add/remove vertical grid lines
};

// Unimplemented
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
