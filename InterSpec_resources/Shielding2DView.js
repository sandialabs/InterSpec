/* Shielding2DView.js */

Shielding2DView = function(id, data, options) {
  var container = document.getElementById(id);
  if (!container) return;

  var self = this;
  this.container = container; // Store container reference
  
  // Store data
  this.data = null;
  
  // Options with defaults
  this.options = options || {};
  if (!this.options.padding) this.options.padding = {};
  if (typeof this.options.padding.top !== "number") this.options.padding.top = 40;
  if (typeof this.options.padding.right !== "number") this.options.padding.right = 40;
  if (typeof this.options.padding.bottom !== "number") this.options.padding.bottom = 40;
  if (typeof this.options.padding.left !== "number") this.options.padding.left = 40;
  
  // Display options
  this.showDetector = true;
  this.showLabels = true; // true = always, false = hover only
  this.pieAngle = 35; // Decreased to 35 degrees
  
  // Processed data (layers, bounds)
  this.processedLayers = [];
  this.processedBounds = null;
  
  // Create Controls
  var controls = document.createElement('div');
  controls.className = 'shielding-controls';
  
  var lblDet = document.createElement('label');
  var chkDet = document.createElement('input');
  chkDet.type = 'checkbox';
  chkDet.checked = this.showDetector;
  chkDet.onchange = function() { self.showDetector = this.checked; self.render(); };
  lblDet.appendChild(chkDet);
  lblDet.appendChild(document.createTextNode(' Show Detector'));
  controls.appendChild(lblDet);
  
  controls.appendChild(document.createElement('br'));

  var lblLbl = document.createElement('label');
  var chkLbl = document.createElement('input');
  chkLbl.type = 'checkbox';
  chkLbl.checked = this.showLabels;
  chkLbl.onchange = function() { self.showLabels = this.checked; self.render(); };
  lblLbl.appendChild(chkLbl);
  lblLbl.appendChild(document.createTextNode(' Always Show Labels'));
  controls.appendChild(lblLbl);

  container.appendChild(controls);

  // Create Tooltip
  this.tooltip = d3.select(container).append("div")
    .attr("class", "shielding-tooltip");

  // Create SVG
  this.svg = d3.select(container).append("svg")
    .attr("width", "100%")
    .attr("height", "100%")
    .style("background-color", "var(--d3spec-background-color, none)");

  this.g = this.svg.append("g");
  
  // Set initial data if provided
  if (data) {
    this.setData(data);
  }
};

// Set data and process it
Shielding2DView.prototype.setData = function(data) {
  this.data = data;
  this.processData();
  
  // Check if detector distance is more than 4x any outer dimension
  // If so, automatically uncheck "Show Detector"
  if (this.processedLayers && this.processedLayers.length > 0) {
      var maxOuterDim = 0;
      var geo = data.geometry;
      
      // Find the maximum outer dimension
      for (var i = 0; i < this.processedLayers.length; i++) {
          var layer = this.processedLayers[i];
          if (geo === "Spherical") {
              maxOuterDim = Math.max(maxOuterDim, layer.outerDims[0]);
          } else if (geo === "CylinderEndOn" || geo === "CylinderSideOn") {
              maxOuterDim = Math.max(maxOuterDim, layer.outerDims[0]); // radius
              maxOuterDim = Math.max(maxOuterDim, 2 * layer.outerDims[1]); // full length
          } else if (geo === "Rectangular") {
              maxOuterDim = Math.max(maxOuterDim, layer.outerDims[0]); // width
              maxOuterDim = Math.max(maxOuterDim, layer.outerDims[1]); // height
              maxOuterDim = Math.max(maxOuterDim, layer.outerDims[2]); // depth
          }
      }
      
      // Check if detector distance > 4x max outer dimension
      if (data.distance > 4 * maxOuterDim) {
          this.showDetector = false;
          // Update checkbox if container exists and checkbox is found
          if (this.container) {
              var chkDet = this.container.querySelector('input[type="checkbox"]');
              if (chkDet) {
                  chkDet.checked = false;
              }
          }
      }
  }
  
  this.render();
};

// Process data into layers and bounds
Shielding2DView.prototype.processData = function() {
  if (!this.data) return;
  
  var data = this.data;
  var self = this;
  
  // Initialize bounds
  var bounds = { minX: 0, maxX: 0, minY: 0, maxY: 0 };
  
  // Shielding Bounds
  var currentRad = 0;
  var currentCylRad = 0, currentCylHalfLen = 0;
  var currentW = 0, currentH = 0, currentD = 0;
  
  // Track layers for drawing
  var layers = [];
  
  for(var i=0; i<data.shieldings.length; i++) {
      var s = data.shieldings[i];
      var dims = s.dimensions;
      var geo = data.geometry;
      
      var layer = { 
          data: s, 
          index: i,
          isGeneric: (s.material.indexOf("Generic") === 0 && s.density === 0),
          path: "",
          innerDims: [],
          outerDims: []
      };
      
      if (geo === "Spherical") {
           // Wedge pointing UP (-Y)
           var innerR = currentRad;
           if (!layer.isGeneric) {
               currentRad += dims[0];
           }
           var outerR = layer.isGeneric ? innerR : currentRad; // Generic has zero thickness
           
           layer.innerDims = [innerR];
           layer.outerDims = [outerR];
           
           var rads = self.pieAngle * Math.PI / 180.0;
           
           var arc = d3.svg.arc()
              .innerRadius(innerR)
              .outerRadius(outerR)
              .startAngle(-rads/2)
              .endAngle(rads/2);
              
           layer.path = arc();
             layer.centroid = arc.centroid();
             // Connection point for label line: use centroid for spherical
             layer.connectionPoint = arc.centroid();
           
           // For spherical, we only draw from center (0) upward to -outerR
           bounds.minY = Math.min(bounds.minY, -outerR);
           bounds.maxY = Math.max(bounds.maxY, 0); // Center of sphere
           bounds.minX = Math.min(bounds.minX, -outerR * Math.sin(rads/2));
           bounds.maxX = Math.max(bounds.maxX, outerR * Math.sin(rads/2));
           
      } else if (geo === "CylinderEndOn") {
           // Rectangle. Z axis horizontal.
           var rThick = dims[0];
           var lThick = dims[1];
           
           var innerR = currentCylRad;
           var innerL = currentCylHalfLen;
           
           if (!layer.isGeneric) {
               currentCylRad += rThick;
               currentCylHalfLen += lThick;
           }
           
           var outerR = layer.isGeneric ? innerR : currentCylRad; // Generic has zero thickness
           var outerL = layer.isGeneric ? innerL : currentCylHalfLen;
           
           layer.innerDims = [innerR, innerL];
           layer.outerDims = [outerR, outerL];
           
           layer.shape = "rect";
           layer.outer = { x: -outerL, y: -outerR, w: 2*outerL, h: 2*outerR };
           
           // Path string (Outer minus Inner if applicable)
           layer.path = "M " + (-outerL) + " " + (-outerR) + 
                        " L " + (outerL) + " " + (-outerR) +
                        " L " + (outerL) + " " + (outerR) +
                        " L " + (-outerL) + " " + (outerR) + " Z";
           
           if (innerL > 0) {
               layer.path += " M " + (-innerL) + " " + (-innerR) +
                             " L " + (-innerL) + " " + (innerR) +
                             " L " + (innerL) + " " + (innerR) +
                             " L " + (innerL) + " " + (-innerR) + " Z";
           }
           
           layer.centroid = [0, outerR];
           // Connection point for label line: right edge, in top half of shell
           // Top half means negative Y (in SVG, negative Y is up)
           // Use a point between -outerR (top edge) and -innerR (top of inner edge if exists)
           var topY = -outerR + (outerR - innerR) * 0.25; // 25% down from top edge (in top half)
           layer.connectionPoint = [outerL, topY];
           
           bounds.maxX = Math.max(bounds.maxX, outerL);
           bounds.minX = Math.min(bounds.minX, -outerL);
           bounds.maxY = Math.max(bounds.maxY, outerR);
           bounds.minY = Math.min(bounds.minY, -outerR);
           
      } else if (geo === "CylinderSideOn") {
          // Circle.
          var rThick = dims[0];
          var lThick = dims[1]; // Length thickness
          var innerR = currentCylRad;
          var innerL = currentCylHalfLen;
          if (!layer.isGeneric) {
              currentCylRad += rThick;
              currentCylHalfLen += lThick;
          }
          var outerR = layer.isGeneric ? innerR : currentCylRad; // Generic has zero thickness
          var outerL = layer.isGeneric ? innerL : currentCylHalfLen;
          
          layer.innerDims = [innerR, innerL];
          layer.outerDims = [outerR, outerL];
          
          var arc = d3.svg.arc().innerRadius(innerR).outerRadius(outerR).startAngle(0).endAngle(2*Math.PI);
          layer.path = arc();
          layer.centroid = [0, outerR];
          // Connection point for label line: halfway through layer on upper half (top, at mid-radius)
          // In SVG, negative Y is up, so top of circle is at negative Y
          var midR = (innerR + outerR) / 2.0;
          layer.connectionPoint = [0, -midR]; // Top of circle at mid-radius (negative Y = up)
          
          bounds.maxX = Math.max(bounds.maxX, outerR);
          bounds.minX = Math.min(bounds.minX, -outerR);
          bounds.maxY = Math.max(bounds.maxY, outerR);
          bounds.minY = Math.min(bounds.minY, -outerR);
          
      } else if (geo === "Rectangular") {
          var wThick = dims[0]; var hThick = dims[1]; var dThick = dims[2];
          var innerW = currentW; var innerH = currentH; var innerD = currentD;
          if (!layer.isGeneric) {
              currentW += wThick; currentH += hThick; currentD += dThick;
          }
          var outerW = layer.isGeneric ? innerW : currentW; // Generic has zero thickness
          var outerH = layer.isGeneric ? innerH : currentH;
          var outerD = layer.isGeneric ? innerD : currentD;
          layer.innerDims = [innerW, innerH, innerD]; layer.outerDims = [outerW, outerH, outerD];
          
          layer.path = "M " + (-outerD) + " " + (-outerW) + 
                       " L " + (outerD) + " " + (-outerW) +
                       " L " + (outerD) + " " + (outerW) +
                       " L " + (-outerD) + " " + (outerW) + " Z";
          
          if (innerD > 0) {
               layer.path += " M " + (-innerD) + " " + (-innerW) +
                             " L " + (-innerD) + " " + (innerW) +
                             " L " + (innerD) + " " + (innerW) +
                             " L " + (innerD) + " " + (-innerW) + " Z";
          }
          layer.centroid = [0, outerW];
          // Connection point for label line: top edge, right side (outside line of shell)
          // Top edge is at y = -outerW (negative Y = up in SVG)
          layer.connectionPoint = [outerD, -outerW];
          
          bounds.maxX = Math.max(bounds.maxX, outerD);
          bounds.minX = Math.min(bounds.minX, -outerD);
          bounds.maxY = Math.max(bounds.maxY, outerW);
          bounds.minY = Math.min(bounds.minY, -outerW);
      }
      
      layers.push(layer);
  }
  
  // Store shielding-only bounds (without detector)
  this.processedLayers = layers;
  this.processedBoundsShieldingOnly = {
    minX: bounds.minX,
    maxX: bounds.maxX,
    minY: bounds.minY,
    maxY: bounds.maxY
  };
  
  // Store detector info for bounds calculation
  this.detectorInfo = {
    distance: data.distance,
    diameter: (data.detectorDiameter !== undefined) ? data.detectorDiameter : (3.0 * 25.4),
    geometry: data.geometry
  };
  
  // Calculate bounds with detector (will be recalculated in render() based on showDetector)
  this.processedBounds = this.calculateBounds(this.processedBoundsShieldingOnly, this.detectorInfo, this.showDetector);
};

// Calculate bounds including detector if needed
Shielding2DView.prototype.calculateBounds = function(shieldingBounds, detectorInfo, includeDetector) {
  var bounds = {
    minX: shieldingBounds.minX,
    maxX: shieldingBounds.maxX,
    minY: shieldingBounds.minY,
    maxY: shieldingBounds.maxY
  };
  
  if (includeDetector && detectorInfo) {
    var detDist = detectorInfo.distance;
    var detDiam = detectorInfo.diameter;
    var detRad = detDiam / 2.0;
    var detLen = detDiam; 
    
    if (detectorInfo.geometry === "Spherical") {
        // Top (-Y)
        bounds.minY = Math.min(bounds.minY, -detDist - detLen);
        bounds.maxX = Math.max(bounds.maxX, detRad);
        bounds.minX = Math.min(bounds.minX, -detRad);
    } else {
        // Right (+X)
        bounds.maxX = Math.max(bounds.maxX, detDist + detLen);
        bounds.maxY = Math.max(bounds.maxY, detRad);
        bounds.minY = Math.min(bounds.minY, -detRad);
    }
  }
  
  return bounds;
};

// Render the graphic
Shielding2DView.prototype.render = function() {
  if (!this.data || !this.processedLayers || !this.processedBoundsShieldingOnly) {
    this.processData();
  }
  
  // Recalculate bounds based on current showDetector state
  var bounds = this.calculateBounds(this.processedBoundsShieldingOnly, this.detectorInfo, this.showDetector);
  
  var self = this;
  var data = this.data;
  var layers = this.processedLayers;
  
  // Clear previous render
  this.g.selectAll("*").remove();
  
  // Get container dimensions
  var container = this.svg.node().parentElement;
  var width = container.clientWidth;
  var height = container.clientHeight;
  
  // Helper function to position tooltip within viewport
  function positionTooltip(tooltip, eventX, eventY, width, height) {
      var tooltipWidth = 200; // Approximate tooltip width
      var tooltipHeight = 150; // Approximate tooltip height
      var offset = 10;
      
      var left = eventX + offset;
      var top = eventY - offset;
      
      // Adjust if tooltip would go off right edge
      if (left + tooltipWidth > width) {
          left = eventX - tooltipWidth - offset;
      }
      
      // Adjust if tooltip would go off left edge
      if (left < 0) {
          left = offset;
      }
      
      // Adjust if tooltip would go off bottom edge
      if (top + tooltipHeight > height) {
          top = eventY - tooltipHeight - offset;
      }
      
      // Adjust if tooltip would go off top edge
      if (top < 0) {
          top = offset;
      }
      
      tooltip.style("left", left + "px")
             .style("top", top + "px");
  }
  
  function formatDensity(d) {
      if (d === 0) return "0.00";
      if (Math.abs(d) < 0.01) return d3.format(".2e")(d);
      return d3.format(".3g")(d);
  }
  
  function formatLength(length_mm, maxDecimals) {
      if (maxDecimals === undefined) maxDecimals = 2;
      if (length_mm === 0) return "0.00 mm";
      
      // Convert mm to base units for comparison
      var length_nm = length_mm * 1e6;
      var length_um = length_mm * 1e3;
      var length_cm = length_mm / 10;
      var length_m = length_mm / 1000;
      var length_km = length_mm / 1e6;
      
      if (length_nm < 1000) {
          return d3.format("." + maxDecimals + "f")(length_nm) + " nm";
      } else if (length_um < 1000) {
          return d3.format("." + maxDecimals + "f")(length_um) + " um";
      } else if (length_mm < 10) {
          return d3.format("." + maxDecimals + "f")(length_mm) + " mm";
      } else if (length_cm < 100) {
          return d3.format("." + maxDecimals + "f")(length_cm) + " cm";
      } else if (length_m < 1000) {
          return d3.format("." + maxDecimals + "f")(length_m) + " m";
      } else {
          return d3.format("." + maxDecimals + "f")(length_km) + " km";
      }
  }
  
  function formatMass(mass_g, maxDecimals) {
      if (maxDecimals === undefined) maxDecimals = 2;
      if (mass_g === 0) return "0.00 g";
      
      var mass_pg = mass_g * 1e12;
      var mass_ng = mass_g * 1e9;
      var mass_ug = mass_g * 1e6;
      var mass_mg = mass_g * 1e3;
      var mass_kg = mass_g / 1000;
      
      if (mass_g < 1e-9) {
          return d3.format("." + maxDecimals + "f")(mass_pg) + " pg";
      } else if (mass_g < 1e-6) {
          return d3.format("." + maxDecimals + "f")(mass_ng) + " ng";
      } else if (mass_g < 1e-3) {
          return d3.format("." + maxDecimals + "f")(mass_ug) + " ug";
      } else if (mass_g < 1.0) {
          return d3.format("." + maxDecimals + "f")(mass_mg) + " mg";
      } else if (mass_g < 1000) {
          return d3.format("." + maxDecimals + "f")(mass_g) + " g";
      } else {
          return d3.format("." + maxDecimals + "f")(mass_kg) + " kg";
      }
  }
  
  // Calculate Scale & Translation
  var padTop = self.options.padding.top;
  var padRight = self.options.padding.right;
  var padBottom = self.options.padding.bottom;
  var padLeft = self.options.padding.left;
  
  var availableW = width - padLeft - padRight;
  var availableH = height - padTop - padBottom;
  
  var contentW = bounds.maxX - bounds.minX;
  var contentH = bounds.maxY - bounds.minY;
  
  if (contentW <= 0) contentW = 1;
  if (contentH <= 0) contentH = 1;
  
  var scaleX = availableW / contentW;
  var scaleY = availableH / contentH;
  var scale = Math.min(scaleX, scaleY);
  
  if (scale <= 0) scale = 1;
  
  var cx = (bounds.minX + bounds.maxX) / 2;
  var cy = (bounds.minY + bounds.maxY) / 2;
  
  // Store scale for manual coordinate scaling
  this.currentScale = scale;
  this.currentTx = (width / 2) - cx * scale;
  this.currentTy = (height / 2) - cy * scale;
  
  // Helper functions to scale coordinates
  var scaleCoord = function(x) { return x * scale; };
  var scaleX = function(x) { return x * scale; };
  var scaleY = function(y) { return y * scale; };
  
  // Helper function to scale path coordinates in path strings
  function scalePathString(pathStr, scaleVal) {
      if (!pathStr) return "";
      // Scale all numbers in the path string
      // This regex matches numbers (including negative and decimals) in path commands
      return pathStr.replace(/([-\d.]+)/g, function(match) {
          return parseFloat(match) * scaleVal;
      });
  }
  
  // Apply only translation (no scale transform)
  this.g.attr("transform", "translate(" + this.currentTx + "," + this.currentTy + ")");
  
  // Draw Layers - separate material and generic layers
  var colorScale = d3.scale.category10();
  var materialLayers = layers.filter(function(d) { return !d.isGeneric; });
  var genericLayers = layers.filter(function(d) { return d.isGeneric; });
  
  // First pass: Draw material layers
  var materialLayerGroups = this.g.selectAll(".material-layer").data(materialLayers).enter()
      .append("g").attr("class", "layer material-layer");
      
  materialLayerGroups.append("path")
      .attr("class", function(d) { return d.isGeneric ? "generic-layer-path" : "shielding-layer-path"; })
      .attr("d", function(d) { 
          // For material layers, recreate path with scaled coordinates
          if (data.geometry === "Spherical") {
              // Recreate arc with scaled radii
              var rads = self.pieAngle * Math.PI / 180.0;
              var arc = d3.svg.arc()
                  .innerRadius(d.innerDims[0] * scale)
                  .outerRadius(d.outerDims[0] * scale)
                  .startAngle(-rads/2)
                  .endAngle(rads/2);
              return arc();
          } else if (data.geometry === "CylinderSideOn") {
              // Recreate arc with scaled radii
              var arc = d3.svg.arc()
                  .innerRadius(d.innerDims[0] * scale)
                  .outerRadius(d.outerDims[0] * scale)
                  .startAngle(0)
                  .endAngle(2*Math.PI);
              return arc();
          } else {
              // For rectangular paths, scale the path string coordinates
              return scalePathString(d.path, scale);
          }
      })
      .attr("fill", function(d, i) { 
          if (d.isGeneric) return "none";
          return colorScale(i); 
      })
      .attr("stroke", function(d) {
           if (d.isGeneric) return colorScale(d.index);
           return "var(--d3spec-axis-color, black)";
      })
      .style("stroke-width", function(d) {
          return (d.isGeneric ? 2 : 1) + "px"; // Fixed pixel width (scale already applied to group)
      })
      .attr("fill-opacity", function(d) {
          if (d.isGeneric) return 0;
          var dens = d.data.density; 
          var op = 0.2 + (dens / 15.0);
          if (op > 0.9) op = 0.9;
          return op;
      })
      .on("mouseover", function(d) {
           d3.select(this).attr("stroke", "orange").style("stroke-width", "3px"); // Fixed pixel width
           
           // Construct Tooltip
           var html = "<b>" + d.data.material + "</b><br/>";
           
           // For generic layers, show AN/AD instead of density, and skip mass
           if (d.isGeneric && d.data.atomicNumber !== undefined && d.data.arealDensity !== undefined) {
               html += "AN=" + d3.format(".1f")(d.data.atomicNumber) + ", AD=" + formatDensity(d.data.arealDensity) + " g/cm2<br/>";
           } else {
               html += "Density: " + formatDensity(d.data.density) + " g/cm3<br/>";
               html += "Mass: " + formatMass(d.data.mass) + "<br/>";
           }
           
           // Calculate and format thickness
           var thickness = [];
           if (d.innerDims.length > 0 && d.outerDims.length > 0) {
               for (var i = 0; i < Math.min(d.innerDims.length, d.outerDims.length); i++) {
                   thickness.push(d.outerDims[i] - d.innerDims[i]);
               }
               if (thickness.length > 0) {
                   html += "Thickness: " + thickness.map(v => formatLength(v, 1)).join(" x ") + "<br/>";
               }
           }
           
           // Format inner and outer dimensions with indicators
           var innerStr = "", outerStr = "";
           if (data.geometry === "Spherical") {
               if (d.innerDims.length > 0) {
                   innerStr = "Inner: " + formatLength(d.innerDims[0], 1);
               }
               outerStr = "Outer: " + formatLength(d.outerDims[0], 1);
           } else if (data.geometry === "CylinderEndOn" || data.geometry === "CylinderSideOn") {
               if (d.innerDims.length >= 2) {
                   var innerR = d.innerDims[0];
                   var innerFullL = 2 * d.innerDims[1];
                   innerStr = "Inner(RxL): " + formatLength(innerR, 1) + " x " + formatLength(innerFullL, 1);
               }
               var outerR = d.outerDims[0];
               var outerFullL = 2 * d.outerDims[1];
               outerStr = "Outer(RxL): " + formatLength(outerR, 1) + " x " + formatLength(outerFullL, 1);
           } else if (data.geometry === "Rectangular") {
               if (d.innerDims.length >= 3) {
                   innerStr = "Inner(WxHxD): " + formatLength(d.innerDims[0], 1) + " x " + formatLength(d.innerDims[1], 1) + " x " + formatLength(d.innerDims[2], 1);
               }
               outerStr = "Outer(WxHxD): " + formatLength(d.outerDims[0], 1) + " x " + formatLength(d.outerDims[1], 1) + " x " + formatLength(d.outerDims[2], 1);
           } else {
               // Fallback to generic format
               if (d.innerDims.length > 0) {
                   innerStr = "Inner: " + d.innerDims.map(v => formatLength(v, 1)).join(" x ");
               }
               outerStr = "Outer: " + d.outerDims.map(v => formatLength(v, 1)).join(" x ");
           }
           
           if (innerStr) {
               html += innerStr + "<br/>";
           }
           html += outerStr + "<br/>";
           
           if (d.data.sources && d.data.sources.length > 0) {
               html += "<hr/>Sources:<br/>" + d.data.sources.join("<br/>");
           }
           
           var eventX = d3.event.pageX - container.getBoundingClientRect().left;
           var eventY = d3.event.pageY - container.getBoundingClientRect().top;
           
           self.tooltip.html(html)
              .style("opacity", 1);
           positionTooltip(self.tooltip, eventX, eventY, width, height);
              
           self.g.select("#label-" + d.index).attr("font-weight", "bold").attr("fill", "orange");
      })
      .on("mouseout", function(d) {
           d3.select(this).attr("stroke", "var(--d3spec-axis-color, black)")
              .style("stroke-width", "1px"); // Fixed pixel width
           self.tooltip.style("opacity", 0);
           self.g.select("#label-" + d.index).attr("font-weight", "normal").attr("fill", "var(--d3spec-text-color, black)");
      });
  
  // Second pass: Draw generic layers on top as 2px lines
  var genericLayerGroups = this.g.selectAll(".generic-layer").data(genericLayers).enter()
      .append("g").attr("class", "layer generic-layer");
      
  genericLayerGroups.append("path")
      .attr("class", "generic-layer-path")
      .attr("d", function(d) { 
          // For generic layers, only show outer boundary as a line
          if (data.geometry === "Spherical") {
              // Just the outer arc
              var outerR = d.outerDims[0] * scale;
              var rads = self.pieAngle * Math.PI / 180.0;
              var arc = d3.svg.arc()
                 .innerRadius(outerR)
                 .outerRadius(outerR)
                 .startAngle(-rads/2)
                 .endAngle(rads/2);
              return arc();
          } else if (data.geometry === "CylinderSideOn") {
              // Just the outer circle
              var outerR = d.outerDims[0] * scale;
              var arc = d3.svg.arc()
                  .innerRadius(outerR)
                  .outerRadius(outerR)
                  .startAngle(0)
                  .endAngle(2*Math.PI);
              return arc();
          } else if (data.geometry === "CylinderEndOn") {
              // Just the outer rectangle outline
              var outerR = d.outerDims[0] * scale;
              var outerL = d.outerDims[1] * scale;
              return "M " + scaleX(-outerL) + " " + scaleY(-outerR) + 
                     " L " + scaleX(outerL) + " " + scaleY(-outerR) +
                     " L " + scaleX(outerL) + " " + scaleY(outerR) +
                     " L " + scaleX(-outerL) + " " + scaleY(outerR) + " Z";
          } else if (data.geometry === "Rectangular") {
              // Just the outer rectangle outline
              var outerW = d.outerDims[0] * scale;
              var outerD = d.outerDims[2] * scale;
              return "M " + scaleX(-outerD) + " " + scaleY(-outerW) + 
                     " L " + scaleX(outerD) + " " + scaleY(-outerW) +
                     " L " + scaleX(outerD) + " " + scaleY(outerW) +
                     " L " + scaleX(-outerD) + " " + scaleY(outerW) + " Z";
          }
      })
      .attr("fill", "none")
      .attr("stroke", function(d) { return colorScale(d.index); })
      .style("stroke-width", "2px")
      .on("mouseover", function(d) {
           d3.select(this).attr("stroke", "orange").style("stroke-width", "3px");
           
           // Construct Tooltip for generic layers
           var html = "<b>" + d.data.material + "</b><br/>";
           
           // For generic layers, show AN/AD instead of density, and skip mass
           if (d.data.atomicNumber !== undefined && d.data.arealDensity !== undefined) {
               html += "AN=" + d3.format(".1f")(d.data.atomicNumber) + ", AD=" + formatDensity(d.data.arealDensity) + " g/cm2<br/>";
           } else {
               html += "Density: " + formatDensity(d.data.density) + " g/cm3<br/>";
           }
           
           // Format inner and outer dimensions with indicators
           var innerStr = "", outerStr = "";
           if (data.geometry === "Spherical") {
               if (d.innerDims.length > 0) {
                   innerStr = "Inner: " + formatLength(d.innerDims[0], 1);
               }
               outerStr = "Outer: " + formatLength(d.outerDims[0], 1);
           } else if (data.geometry === "CylinderEndOn" || data.geometry === "CylinderSideOn") {
               if (d.innerDims.length >= 2) {
                   var innerR = d.innerDims[0];
                   var innerFullL = 2 * d.innerDims[1];
                   innerStr = "Inner(RxL): " + formatLength(innerR, 1) + " x " + formatLength(innerFullL, 1);
               }
               var outerR = d.outerDims[0];
               var outerFullL = 2 * d.outerDims[1];
               outerStr = "Outer(RxL): " + formatLength(outerR, 1) + " x " + formatLength(outerFullL, 1);
           } else if (data.geometry === "Rectangular") {
               if (d.innerDims.length >= 3) {
                   innerStr = "Inner(WxHxD): " + formatLength(d.innerDims[0], 1) + " x " + formatLength(d.innerDims[1], 1) + " x " + formatLength(d.innerDims[2], 1);
               }
               outerStr = "Outer(WxHxD): " + formatLength(d.outerDims[0], 1) + " x " + formatLength(d.outerDims[1], 1) + " x " + formatLength(d.outerDims[2], 1);
           }
           
           if (innerStr) {
               html += innerStr + "<br/>";
           }
           html += outerStr + "<br/>";
           
           if (d.data.sources && d.data.sources.length > 0) {
               html += "<hr/>Sources:<br/>" + d.data.sources.join("<br/>");
           }
           
           var eventX = d3.event.pageX - container.getBoundingClientRect().left;
           var eventY = d3.event.pageY - container.getBoundingClientRect().top;
           
           self.tooltip.html(html)
              .style("opacity", 1);
           positionTooltip(self.tooltip, eventX, eventY, width, height);
              
           self.g.select("#label-" + d.index).attr("font-weight", "bold").attr("fill", "orange");
      })
      .on("mouseout", function(d) {
           d3.select(this).attr("stroke", colorScale(d.index)).style("stroke-width", "2px");
           self.tooltip.style("opacity", 0);
           self.g.select("#label-" + d.index).attr("font-weight", "normal").attr("fill", "var(--d3spec-text-color, black)");
      });
  
  // Point Sources
  var pointSources = [];
  if (data.fitSources) {
      for(var i=0; i<data.fitSources.length; i++) {
          if (data.fitSources[i].type === "Point") {
              pointSources.push(data.fitSources[i]);
          }
      }
  }
  
  if (pointSources.length > 0) {
      // Draw a single circle for all point sources (all at origin)
      var sourceGroup = this.g.append("g").attr("class", "point-sources-group");
      
      var circle = sourceGroup.append("circle")
          .attr("class", "source-dot")
          .attr("cx", scaleX(0))
          .attr("cy", scaleY(0))
          .attr("r", 3) // Fixed pixel radius
          .attr("fill", "red")
          .attr("stroke", "black")
          .style("stroke-width", "1px") // Fixed pixel width
          .style("cursor", "pointer");
      
      // Build tooltip content for ALL point sources
      var tooltipHtml = "<b>Point Source";
      if (pointSources.length > 1) {
          tooltipHtml += "s (" + pointSources.length + ")";
      }
      tooltipHtml += "</b><br/>";
      
      pointSources.forEach(function(src, idx) {
          if (pointSources.length > 1) {
              tooltipHtml += "<hr/><b>Source " + (idx + 1) + ":</b><br/>";
          }
          
          if (src.nuclide) {
              tooltipHtml += "Nuclide: " + src.nuclide + "<br/>";
          }
          
          if (src.activityPretty !== undefined) {
              tooltipHtml += "Activity: " + src.activityPretty;
              if (src.activityUncert !== undefined && src.activityUncert > 0) {
                  // Format uncertainty - try to match the units of activityPretty
                  // For simplicity, just show as number with appropriate precision
                  var uncertStr = "";
                  if (src.activityUncert >= 1e9) {
                      uncertStr = d3.format(".2f")(src.activityUncert / 1e9) + " GBq";
                  } else if (src.activityUncert >= 1e6) {
                      uncertStr = d3.format(".2f")(src.activityUncert / 1e6) + " MBq";
                  } else if (src.activityUncert >= 1e3) {
                      uncertStr = d3.format(".2f")(src.activityUncert / 1e3) + " kBq";
                  } else {
                      uncertStr = d3.format(".2f")(src.activityUncert) + " Bq";
                  }
                  tooltipHtml += " ± " + uncertStr;
              }
              tooltipHtml += "<br/>";
          } else if (src.activity !== undefined) {
              // Fallback: format activity as Bq with appropriate units
              var actStr = "";
              if (src.activity >= 1e9) {
                  actStr = d3.format(".2f")(src.activity / 1e9) + " GBq";
              } else if (src.activity >= 1e6) {
                  actStr = d3.format(".2f")(src.activity / 1e6) + " MBq";
              } else if (src.activity >= 1e3) {
                  actStr = d3.format(".2f")(src.activity / 1e3) + " kBq";
              } else {
                  actStr = d3.format(".2f")(src.activity) + " Bq";
              }
              tooltipHtml += "Activity: " + actStr + "<br/>";
          }
          
          if (src.age !== undefined && src.age > 0) {
              // Age is typically in seconds, convert to readable format
              var ageYears = src.age / (365.25 * 24 * 3600);
              var ageDays = src.age / (24 * 3600);
              var ageHours = src.age / 3600;
              var ageStr = "";
              if (ageYears >= 1) {
                  ageStr = d3.format(".2f")(ageYears) + " years";
              } else if (ageDays >= 1) {
                  ageStr = d3.format(".2f")(ageDays) + " days";
              } else if (ageHours >= 1) {
                  ageStr = d3.format(".2f")(ageHours) + " hours";
              } else {
                  ageStr = d3.format(".2f")(src.age) + " seconds";
              }
              tooltipHtml += "Age: " + ageStr + "<br/>";
          }
      });
      
      circle.on("mouseover", function() {
           d3.select(this).attr("r", 5).attr("fill", "orange"); // Highlight on hover
           
           var eventX = d3.event.pageX - container.getBoundingClientRect().left;
           var eventY = d3.event.pageY - container.getBoundingClientRect().top;
           
           self.tooltip.html(tooltipHtml)
              .style("opacity", 1);
           positionTooltip(self.tooltip, eventX, eventY, width, height);
      })
      .on("mouseout", function() {
           d3.select(this).attr("r", 3).attr("fill", "red"); // Restore original
           self.tooltip.style("opacity", 0);
      });
  }
  
  // Calculate shielding-only maxX (without detector) for label positioning
  var shieldingMaxX = self.processedBoundsShieldingOnly.maxX;
  
  // Labels
  var labelData = [];
  layers.forEach(function(layer) {
      // Store the connection point Y coordinate in model coordinates
      var yPos = 0;
      if (layer.connectionPoint) {
          yPos = layer.connectionPoint[1]; // Y coordinate in model space
      } else if (data.geometry === "Spherical") {
          yPos = -layer.outerDims[0]; 
      } else if (data.geometry === "CylinderSideOn") {
          yPos = -layer.outerDims[0];
      } else {
          // For other geometries, use outer dimension
          yPos = layer.outerDims[0];
      }
      labelData.push({ layer: layer, yModel: yPos });
  });
  
  // Sort so outermost layer is at top
  // For spherical and cylinder side-on, outer layers have larger radius (more negative y)
  // For cylindrical end-on and rectangular, larger outerDims[0] means outer layer
  if (data.geometry === "Spherical" || data.geometry === "CylinderSideOn") {
      labelData.sort((a,b) => a.yModel - b.yModel); // More negative (outer) first
  } else {
      // For CylinderEndOn and Rectangular, sort by outer dimension (largest first = outermost)
      labelData.sort((a,b) => {
          var aOuter = a.layer.outerDims[0];
          var bOuter = b.layer.outerDims[0];
          return bOuter - aOuter; // Larger outer dimension = outer layer
      });
  }
  
  // Scale to screen coordinates
  labelData.forEach(function(ld) {
      ld.y = scaleY(ld.yModel); // Scale to screen coordinates
  });
  
  // For all geometries, ensure labels are ordered top to bottom (outermost at top)
  // and have proper spacing
  if (labelData.length > 1) {
      // Find the topmost and bottommost positions
      var allYs = labelData.map(function(ld) { return ld.y; });
      var topY = Math.min.apply(Math, allYs);
      var bottomY = Math.max.apply(Math, allYs);
      
      // For rectangular and CylinderEndOn, redistribute from top to bottom
      // For others, just ensure spacing
      if (data.geometry === "Rectangular" || data.geometry === "CylinderEndOn") {
          // Redistribute evenly from top to bottom
          for (var j = 0; j < labelData.length; j++) {
              var t = j / (labelData.length - 1);
              if (labelData.length === 1) t = 0;
              labelData[j].y = topY + (bottomY - topY) * t;
          }
      } else {
          // For other geometries, just ensure they're in order
          labelData.sort(function(a, b) { return a.y - b.y; });
      }
  }
  
  var labH = 50; // Preferred pixel spacing between labels (for three-line labels)
  var minSpacing = 10; // Minimum spacing between labels (must be at least this)
  for(var i=1; i<labelData.length; i++) {
      var requiredSpacing = Math.max(labH, minSpacing);
      
      // Calculate current spacing
      var currentSpacing = labelData[i].y - labelData[i-1].y;
      
      if (data.geometry === "Spherical" || data.geometry === "CylinderSideOn") {
          // For spherical and side-on, labels go downward (increasing Y, less negative)
          // Outer layer is at top (more negative Y), inner layers below (less negative Y)
          // After scaling: outer layer has smaller Y (top), inner layers have larger Y (bottom)
          if (currentSpacing < requiredSpacing) {
              labelData[i].y = labelData[i-1].y + requiredSpacing;
          }
      } else if (data.geometry === "CylinderEndOn" || data.geometry === "Rectangular") {
          // For end-on and rectangular, labels go from top (outermost) to bottom (innermost)
          // After scaling: labelData[0] should have smallest Y (top), labelData[i] should have larger Y (below)
          // Ensure labelData[i] is below labelData[i-1] with required spacing
          if (currentSpacing < requiredSpacing) {
              labelData[i].y = labelData[i-1].y + requiredSpacing;
          }
          // Also ensure it's not above the previous one
          if (labelData[i].y <= labelData[i-1].y) {
              labelData[i].y = labelData[i-1].y + requiredSpacing;
          }
      } else {
          // For other geometries, labels go upward (increasing Y, more positive)
          if (currentSpacing < requiredSpacing) {
              labelData[i].y = labelData[i-1].y + requiredSpacing;
          }
      }
  }
  
  labelData.forEach(function(ld) {
      var layer = ld.layer;
      var ly = ld.y; // Already scaled to screen coordinates
      
      // Position labels to the right of the shielding geometry
      // Use the maximum X of the shielding bounds (scaled) plus offset
      var lx = scaleX(shieldingMaxX) + 20; // Fixed 20px offset
      
      var lg = self.g.append("g")
          .attr("class", "label-group")
          .attr("id", "label-group-" + layer.index)
          .style("opacity", self.showLabels ? 1 : 0);
          
      // Use connection point for label line (top of shell) - scale coordinates
      var cx = 0, cy = 0;
      if (layer.connectionPoint) {
          cx = scaleX(layer.connectionPoint[0]);
          cy = scaleY(layer.connectionPoint[1]);
      } else if (layer.centroid) {
          cx = scaleX(layer.centroid[0]);
          cy = scaleY(layer.centroid[1]);
      }
      
      lg.append("polyline")
          .attr("class", "label-line")
          .attr("points", cx+","+cy + " " + (lx-5)+","+cy + " " + (lx-5)+","+ly + " " + lx+","+ly)
          .attr("fill", "none")
          .attr("stroke", "var(--d3spec-text-color, black)")
          .style("stroke-width", "1px"); // Fixed pixel width for UI elements
          
      // Create three-line label
      var labelText = lg.append("text")
          .attr("class", "fixed-size-text")
          .attr("id", "label-" + layer.index)
          .attr("x", lx + 2)
          .attr("y", ly)
          .attr("dy", "0em")
          .style("font-size", "12px") // Fixed pixel font size
          .attr("fill", "var(--d3spec-text-color, black)");
      
      // First line: material name
      labelText.append("tspan")
          .attr("x", lx + 2)
          .attr("dy", "0em")
          .text(layer.data.material);
      
      // Second line: density (or AN/AD for generic)
      var secondLineText = "";
      if (layer.isGeneric && layer.data.atomicNumber !== undefined && layer.data.arealDensity !== undefined) {
          // For generic layers, show AN and AD
          secondLineText = "AN=" + d3.format(".1f")(layer.data.atomicNumber) + ", AD=" + formatDensity(layer.data.arealDensity) + " g/cm2";
      } else {
          // For material layers, show density
          secondLineText = "\u03C1=" + formatDensity(layer.data.density) + " g/cm3";
      }
      labelText.append("tspan")
          .attr("x", lx + 2)
          .attr("dy", "1.2em")
          .text(secondLineText);
      
      // Third line: OD
      // OD is the primary outer dimension for each geometry
      var odStr = "";
      if (data.geometry === "Spherical") {
          odStr = "OD=" + formatLength(layer.outerDims[0], 1); // Outer radius
      } else if (data.geometry === "CylinderEndOn") {
          // outerDims[0] is radius, outerDims[1] is half-length, so full length is 2 * half-length
          var radius = layer.outerDims[0];
          var fullLength = 2 * layer.outerDims[1];
          odStr = "OD(RxL)=" + formatLength(radius, 1) + " x " + formatLength(fullLength, 1);
      } else if (data.geometry === "CylinderSideOn") {
          // outerDims[0] is radius, outerDims[1] is half-length, so full length is 2 * half-length
          var radius = layer.outerDims[0];
          var fullLength = 2 * layer.outerDims[1];
          odStr = "OD(RxL)=" + formatLength(radius, 1) + " x " + formatLength(fullLength, 1);
      } else if (data.geometry === "Rectangular") {
          // For rectangular, show all three dimensions: width x height x depth
          odStr = "OD(WxHxD)=" + formatLength(layer.outerDims[0], 1) + " x " + formatLength(layer.outerDims[1], 1) + " x " + formatLength(layer.outerDims[2], 1);
      }
      
      labelText.append("tspan")
          .attr("x", lx + 2)
          .attr("dy", "1.2em")
          .text(odStr);
  });
  
  if (!self.showLabels) {
      this.g.selectAll(".layer")
          .on("mouseover.label", function(d) { self.g.select("#label-group-" + d.index).style("opacity", 1); })
          .on("mouseout.label", function(d) { self.g.select("#label-group-" + d.index).style("opacity", 0); });
  }
  
  // Detector
  if (self.showDetector) {
      var detDist = data.distance;
      var detDiam = (data.detectorDiameter !== undefined) ? data.detectorDiameter : (3.0 * 25.4);
      var detRad = detDiam / 2.0;
      var detLen = detDiam; 
      var detGroup = this.g.append("g").attr("class", "detector-group");

      if (data.geometry === "Spherical") {
          // Top (-Y) - scale coordinates
          detGroup.append("rect")
              .attr("class", "detector-shape detector-rect")
              .attr("x", scaleX(-detRad))
              .attr("y", scaleY(-detDist - detLen))
              .attr("width", scaleX(detDiam))
              .attr("height", scaleY(detLen))
              .attr("fill", "#ccc")
              .attr("stroke", "var(--d3spec-axis-color, black)")
              .style("stroke-width", "1px"); // Fixed pixel width
              
          detGroup.append("text")
              .attr("x", scaleX(0))
              .attr("y", scaleY(-detDist - detLen/2))
              .attr("text-anchor", "middle")
              .attr("dy", "0.3em")
              .text("Detector")
              .style("font-size", Math.min(scaleX(detDiam/4), scaleY(detLen/4), 12) + "px") // Scale with box size, but ensure it fits
              .attr("fill", "black")
              .style("pointer-events", "none");
      } else {
          // Right (+X) - scale coordinates
          detGroup.append("rect")
              .attr("class", "detector-shape detector-rect")
              .attr("x", scaleX(detDist))
              .attr("y", scaleY(-detRad))
              .attr("width", scaleX(detLen))
              .attr("height", scaleY(detDiam))
              .attr("fill", "#ccc")
              .attr("stroke", "var(--d3spec-axis-color, black)")
              .style("stroke-width", "1px"); // Fixed pixel width

          detGroup.append("text")
              .attr("x", scaleX(detDist + detLen/2))
              .attr("y", scaleY(0))
              .attr("text-anchor", "middle")
              .attr("dy", "0.3em")
              .text("Detector")
              .style("font-size", Math.min(scaleX(detLen/4), scaleY(detDiam/4), 12) + "px") // Scale with box size, but ensure it fits
              .attr("fill", "black")
              .style("pointer-events", "none");
      }
      
      // Interaction
      detGroup.selectAll(".detector-shape")
          .on("mouseover", function() {
               d3.select(this).attr("stroke", "orange").style("stroke-width", "3px"); // Fixed pixel width
               
               // Find Outer Shielding Dimension
               var outerDim = 0;
               if (layers.length > 0) {
                   layers.forEach(function(l) {
                       if (data.geometry === "Spherical") outerDim = Math.max(outerDim, l.outerDims[0]);
                       else if (data.geometry === "CylinderSideOn") outerDim = Math.max(outerDim, l.outerDims[0]);
                       else if (data.geometry === "CylinderEndOn") outerDim = Math.max(outerDim, l.outerDims[1]);
                       else if (data.geometry === "Rectangular") outerDim = Math.max(outerDim, l.outerDims[1]);
                   });
               }
               
               var pFace = {x:0, y:0}, pCenter = {x:0, y:0}, pSurf = {x:0, y:0};
               
               if (data.geometry === "Spherical") {
                   pFace = {x:0, y:-detDist};
                   pCenter = {x:0, y:0};
                   pSurf = {x:0, y:-outerDim};
                   
                   self.g.append("line").attr("class", "guide-line")
                      .attr("x1", scaleX(-10)).attr("y1", scaleY(-outerDim))
                      .attr("x2", scaleX(10)).attr("y2", scaleY(-outerDim))
                      .attr("stroke", "orange").style("stroke-width", "2px"); // Fixed pixel width
                      
               } else {
                   pFace = {x:detDist, y:0};
                   pCenter = {x:0, y:0};
                   pSurf = {x:outerDim, y:0};
                   
                   self.g.append("line").attr("class", "guide-line")
                      .attr("x1", scaleX(outerDim)).attr("y1", scaleY(-10))
                      .attr("x2", scaleX(outerDim)).attr("y2", scaleY(10))
                      .attr("stroke", "orange").style("stroke-width", "2px"); // Fixed pixel width
               }
               
               self.g.append("line").attr("class", "guide-line")
                  .attr("x1", scaleX(pFace.x)).attr("y1", scaleY(pFace.y))
                  .attr("x2", scaleX(pCenter.x)).attr("y2", scaleY(pCenter.y))
                  .attr("stroke", "orange").style("stroke-width", "1.5px").style("stroke-dasharray", "5,5"); // Fixed pixel width
                  
               var distC = detDist; 
               var distS = detDist - outerDim;
               var txtOff = 15; // Fixed pixel offset
               var pTxtC = {}, pTxtS = {};
               
               if (data.geometry === "Spherical") {
                    pTxtC = { x: scaleX(txtOff), y: scaleY((pFace.y + pCenter.y)/2) };
                    pTxtS = { x: scaleX(-txtOff), y: scaleY((pFace.y + pSurf.y)/2) };
               } else {
                    pTxtC = { x: scaleX((pFace.x + pCenter.x)/2), y: scaleY(-txtOff) };
                    pTxtS = { x: scaleX((pFace.x + pSurf.x)/2), y: scaleY(txtOff) };
               }

               var anchorC = (data.geometry === "Spherical") ? "start" : "middle";
               var anchorS = (data.geometry === "Spherical") ? "end" : "middle";
               var dyC = (data.geometry === "Spherical") ? "0.3em" : "-0.3em";
               var dyS = (data.geometry === "Spherical") ? "0.3em" : "0.9em";

               self.g.append("text").attr("class", "guide-text")
                  .attr("x", pTxtC.x).attr("y", pTxtC.y)
                  .attr("text-anchor", anchorC).attr("dy", dyC)
                  .text("To center: " + formatLength(distC, 1))
                  .style("font-size", "11px").attr("fill", "orange"); // Fixed pixel font size
                  
               self.g.append("text").attr("class", "guide-text")
                  .attr("x", pTxtS.x).attr("y", pTxtS.y)
                  .attr("text-anchor", anchorS).attr("dy", dyS)
                  .text("To surface: " + formatLength(distS, 1))
                  .style("font-size", "11px").attr("fill", "orange"); // Fixed pixel font size
          })
          .on("mouseout", function() {
               d3.select(this).attr("stroke", "var(--d3spec-axis-color, black)").style("stroke-width", "1px"); // Fixed pixel width
               self.g.selectAll(".guide-line").remove();
               self.g.selectAll(".guide-text").remove();
          });
  }
  
  // Axis Indicators - position at lower left of viewing area
  var axisLen = 40; // Fixed pixel length
  var axOffset = 20; // Fixed pixel offset from edges
  // Position in screen coordinates (accounting for translation)
  var axX = -this.currentTx + axOffset;
  var axY = -this.currentTy + (height - axOffset);
  
  if (data.geometry === "CylinderSideOn") {
      // For side-on, show that length is perpendicular to screen
      var ag = this.g.append("g").attr("transform", "translate(" + axX + "," + axY + ")");
      
      // Draw a circle with a dot to indicate "into/out of screen"
      var circleRad = 8; // Fixed pixel radius
      ag.append("circle")
          .attr("class", "axis-line")
          .attr("cx", 0)
          .attr("cy", 0)
          .attr("r", circleRad)
          .attr("fill", "none")
          .attr("stroke", "var(--d3spec-axis-color, black)")
          .style("stroke-width", "2px"); // Fixed pixel width
      
      ag.append("circle")
          .attr("class", "axis-line")
          .attr("cx", 0)
          .attr("cy", 0)
          .attr("r", 2) // Fixed pixel radius
          .attr("fill", "var(--d3spec-axis-color, black)");
      
      ag.append("text")
          .attr("class", "fixed-size-text")
          .attr("x", circleRad + 4) // Fixed 4px offset
          .attr("y", 0)
          .attr("dy", "0.3em")
          .style("font-size", "10px") // Fixed pixel font size
          .attr("fill", "var(--d3spec-text-color, black)")
          .text("Length");
          
  } else if (data.geometry !== "Spherical") {
      var ag = this.g.append("g").attr("transform", "translate(" + axX + "," + axY + ")");
      
      ag.append("line")
          .attr("class", "axis-line")
          .attr("x1", 0).attr("y1", 0)
          .attr("x2", axisLen).attr("y2", 0)
          .attr("stroke", "var(--d3spec-axis-color, black)")
          .style("stroke-width", "2px"); // Fixed pixel width
          
      ag.append("line")
          .attr("class", "axis-line")
          .attr("x1", 0).attr("y1", 0)
          .attr("x2", 0).attr("y2", -axisLen)
          .attr("stroke", "var(--d3spec-axis-color, black)")
          .style("stroke-width", "2px"); // Fixed pixel width
          
      var hLabel = "", vLabel = "";
      if (data.geometry === "CylinderEndOn") {
          hLabel = "Length"; vLabel = "Radial"; 
      } else if (data.geometry === "Rectangular") {
          hLabel = "Depth"; vLabel = "Width"; 
      }
      
      var depthText = ag.append("text")
          .attr("class", "fixed-size-text")
          .attr("x", axisLen + 2) // Fixed 2px offset
          .attr("y", 0)
          .attr("dy", "0.3em")
          .style("font-size", "10px") // Fixed pixel font size
          .attr("fill", "var(--d3spec-text-color, black)")
          .text(hLabel);
          
      ag.append("text")
          .attr("class", "fixed-size-text")
          .attr("x", 0)
          .attr("y", -axisLen - 2) // Fixed 2px offset
          .attr("text-anchor", "middle")
          .style("font-size", "10px") // Fixed pixel font size
          .attr("fill", "var(--d3spec-text-color, black)")
          .text(vLabel);
      
      // For rectangular geometry, add Height indicator to the right
      if (data.geometry === "Rectangular") {
          // Calculate position 20px after the "Depth" text ends
          // The text starts at axisLen + 2, and we need to add the text width + 20px
          var textNode = depthText.node();
          var textWidth = textNode ? textNode.getBBox().width : 0;
          var depthTextEnd = axisLen + 2 + textWidth;
          var heightAxOffset = depthTextEnd + 20; // 20px after the "Depth" text ends
          var heightAg = this.g.append("g").attr("transform", "translate(" + (axX + heightAxOffset) + "," + axY + ")");
          
          // Draw a circle with a dot to indicate "into/out of screen" (Height direction)
          var circleRad = 8; // Fixed pixel radius
          heightAg.append("circle")
              .attr("class", "axis-line")
              .attr("cx", 0)
              .attr("cy", 0)
              .attr("r", circleRad)
              .attr("fill", "none")
              .attr("stroke", "var(--d3spec-axis-color, black)")
              .style("stroke-width", "2px"); // Fixed pixel width
          
          heightAg.append("circle")
              .attr("class", "axis-line")
              .attr("cx", 0)
              .attr("cy", 0)
              .attr("r", 2) // Fixed pixel radius
              .attr("fill", "var(--d3spec-axis-color, black)");
          
          heightAg.append("text")
              .attr("class", "fixed-size-text")
              .attr("x", circleRad + 4) // Fixed 4px offset
              .attr("y", 0)
              .attr("dy", "0.3em")
              .style("font-size", "10px") // Fixed pixel font size
              .attr("fill", "var(--d3spec-text-color, black)")
              .text("Height");
      }
  }
};

// Handle resize
Shielding2DView.prototype.handleResize = function() {
  if (this.data) {
    this.render();
  }
};
