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

// TODO:
//  - [ ] After drawing, or modifying drawings on the map, have buttons that let you load to foreground/background/secondary
//    - [ ] Could probably give some CPS, and live/real time stats from what is selected
//  - [ ] Could have CPS instead be energy ratio ... woudl need to tie into TimeChart for this
//  - [ ] Have option to show points for currently displayed foreground/background/secondary
//  - [ ] ability to specify lower and upper ranges for the color mapping (to highlight a line or other feature, for example), a log scale
        

LeafletRadMap = function (elem,options) {
  let self = this;
  this.parent = typeof elem === "string" ? document.getElementById(elem) : elem;

  if( !this.parent )
    throw new Error("LeafletRadMap: invalid parent identifier passed in");

  if( !this.parent.classList.contains("LeafletRadMap") )
    this.parent.classList.add("LeafletRadMap");
    
  this.mapDiv = document.createElement("div");
  this.mapDiv.classList.add("LeafletRadMapMap");
  this.parent.appendChild( this.mapDiv );

  this.gradientLegendDiv = document.createElement("div");
  this.gradientLegendDiv.classList.add("LeafletRadMapLegend");
  this.parent.appendChild( this.gradientLegendDiv );

  this.bottomDiv = document.createElement("div");
  this.bottomDiv.classList.add("LeafletRadMapBottomDiv");
  this.parent.appendChild( this.bottomDiv );
  
  this.btnsDiv = document.createElement("div");
  this.btnsDiv.classList.add("LeafletRadMapBtnsDiv");
  this.bottomDiv.appendChild( this.btnsDiv );

  this.copyBtnsDiv = document.createElement("div");
  this.copyBtnsDiv.classList.add("CopyMapBtnsDiv");
  this.bottomDiv.appendChild( this.copyBtnsDiv );

  this.markerLegendDiv = document.createElement("div");
  this.markerLegendDiv.classList.add("LeafletRadMapMarkerLegend");
  this.bottomDiv.appendChild( this.markerLegendDiv );
  

  this.options = options ? options : {};
  if( typeof this.options.heatMap !== "boolean" )
    this.options.heatMap = true;
  if( typeof this.options.minMarkerOpacity !== "number" )
    this.options.minMarkerOpacity = 0.1;
  if( typeof this.options.showHeatMap !== "boolean" )
    this.options.showHeatMap = false;
    
  // Text we might need to localize
  if( typeof this.options.foregroundTxt !== "string" )
    this.options.foregroundTxt = "Foreground";
  if( typeof this.options.backgroundTxt !== "string" )
    this.options.backgroundTxt = "Background";
  if( typeof this.options.secondaryTxt !== "string" )
    this.options.secondaryTxt = "Secondary";
  if( typeof this.options.gammaTxt !== "string" )
    this.options.gammaTxt = "Gamma";
  if( typeof this.options.displayedAsTxt !== "string" )
    this.options.displayedAsTxt = "Displayed as";
  if( typeof this.options.cpsText !== "string" )
    this.options.cpsText = "CPS";
  if( typeof this.options.loadTxt !== "string" )
    this.options.loadTxt = "Load";
  if( typeof this.options.measurementsAsTxt !== "string" )
    this.options.measurementsAsTxt = "measurements as";
  if( typeof this.options.loadTxtShort !== "string" )
    this.options.loadTxtShort = "Load as";
  if( typeof this.options.realTimeTxt !== "string" )
    this.options.realTimeTxt = "Real Time";
  if( typeof this.options.liveTimeTxt !== "string" )
    this.options.liveTimeTxt = "Live Time";
    
  this.map = L.map(this.mapDiv, {
    minZoom: 2,
    maxZoom: 23, // the heatmap plugin needs some help understanding the map's maxZoom
    editable: true,
    doubleClickZoom: false
  }).setView([37.67640130843344, -121.70667619429008], 15);  
  
  const basemapLayers = {
    Streets: L.esri.Vector.vectorBasemapLayer("ArcGIS:Streets", { apiKey: self.options.apiKey }).addTo( self.map ),
    Imagery: L.esri.Vector.vectorBasemapLayer("ArcGIS:Imagery", { apiKey: self.options.apiKey }),
    Topographic: L.esri.Vector.vectorBasemapLayer("ArcGIS:Topographic", { apiKey: self.options.apiKey }),
    Navigation: L.esri.Vector.vectorBasemapLayer("ArcGIS:Navigation", { apiKey: self.options.apiKey }),
    "Light Gray": L.esri.Vector.vectorBasemapLayer("ArcGIS:LightGray", { apiKey: self.options.apiKey }),
    "Dark gray": L.esri.Vector.vectorBasemapLayer("ArcGIS:DarkGray", { apiKey: self.options.apiKey }),
    "Streets Relief": L.esri.Vector.vectorBasemapLayer("ArcGIS:StreetsRelief", { apiKey: self.options.apiKey }),
  };
  
  L.control.layers(basemapLayers, null, { collapsed: true }).addTo( self.map );

  // Monitor tile loading - detect failures from proxy, no internet, etc.
  // The vectorBasemapLayer fires a "ready" event when the ArcGIS style is
  // fetched and MapLibre GL is initialized; on network failure this never fires.
  // However, the browser may cache the style, so _ready can be true even when
  // the server is unreachable - we must also check actual tile loading via
  // the MapLibre GL map's error/idle events.
  self._tilesLoaded = false;
  self._tileLoadFailedEmitted = false;

  function emitTileLoadFailed() {
    if( self._tileLoadFailedEmitted )
      return;
    self._tileLoadFailedEmitted = true;
    try {
      const interspecEl = document.querySelector( '.InterSpec' );
      if( interspecEl && window.Wt )
        Wt.emit( interspecEl.id, {name: 'tileLoadFailed'} );
    } catch( e ) {
      console.warn( "LeafletRadMap: error emitting tileLoadFailed:", e );
    }
  }

  function onLayerReady( layer ) {
    // Style loaded (possibly from cache). Hook into the MapLibre GL map to
    // check whether tiles actually load, or only the cached style was used.
    try {
      const glMap = layer._maplibreGL && layer._maplibreGL.getMaplibreMap
                  ? layer._maplibreGL.getMaplibreMap() : null;
      if( glMap ){
        let hadError = false;
        glMap.on( 'error', function() { hadError = true; } );
        glMap.once( 'idle', function() {
          if( !hadError ){
            self._tilesLoaded = true;
            if( self._tileLoadTimer )
              clearTimeout( self._tileLoadTimer );
          }
        } );
        return;
      }
    } catch( e ) {
      // Could not access MapLibre GL internals
    }

    // Fallback: if we cant access MapLibre GL internals, treat ready as loaded
    self._tilesLoaded = true;
    if( self._tileLoadTimer )
      clearTimeout( self._tileLoadTimer );
  }

  function startTileLoadMonitor( layer ) {
    if( self._tileLoadTimer )
      clearTimeout( self._tileLoadTimer );
    self._tilesLoaded = false;

    // If style already loaded (e.g., from browser cache), check tiles now
    if( layer._ready )
      onLayerReady( layer );

    // Also listen for future ready event (for non-cached case)
    layer.once( "ready", function() { onLayerReady( layer ); } );

    // Timeout: if tiles haven't confirmed loaded within 15s, warn user
    self._tileLoadTimer = setTimeout( function() {
      if( !self._tilesLoaded )
        emitTileLoadFailed();
    }, 15000 );
  }

  startTileLoadMonitor( basemapLayers.Streets );

  self.map.on( 'baselayerchange', function( e ) {
    startTileLoadMonitor( e.layer );
  } );

  // Add drawing tools
  this.map.pm.addControls({
    editControls: true,
    drawMarker: false,
    drawPolyline: false,
    drawCircleMarker: false,
    drawText: true,
    drawCircle: false,
    cutPolygon: true,
    removalMode: true,
  });
  
  L.control.scale({metric: true, imperial: false, position: 'bottomleft'}).addTo(this.map);
  
  L.control.measure({position: "bottomleft"}).addTo(this.map);
  
  this.markerGradientStops = {
    0.4: 'blue',
    0.6: 'cyan',
    0.7: 'lime',
    0.9: 'yellow',
    1.0: 'red'
  };
  this.grad_colors = this.defineGradientColors( self.markerGradientStops, self.options.minMarkerOpacity );
  
  // Setup a legend for the color gradient
  this.legendLowerCps = document.createElement("div");
  this.legendLowerCps.classList.add("LeglowVal");
  this.gradientLegendDiv.appendChild( this.legendLowerCps );

  this.legendCanvas = document.createElement('canvas');
  this.legendCanvas.classList.add("LegCanvas");
  this.gradientLegendDiv.appendChild( this.legendCanvas );

  this.legendUpperCps = document.createElement("div");
  this.legendUpperCps.classList.add("LegUpVal");
  this.gradientLegendDiv.appendChild( this.legendUpperCps );

  this.map.on('pm:create', (e) => {
    e.layer.on('pm:edit', (e) => {
      self.handleUserDrawingUpdate();
    } );
  
    self.handleUserDrawingUpdate();
  });
  
  this.map.on('pm:cut', (e) => {
    console.log( "geoman cut" );
    self.handleUserDrawingUpdate();
  } );
        
  this.map.on('pm:remove', (e) => {
    console.log( "geoman remove" );
    self.handleUserDrawingUpdate();
  } );
  
  this.map.on('pm:globalremovalmodetoggled', (e) => {
    console.log('pm:globalremovalmodetoggled');
    self.markerLayer.pm.disable();
  });


  // On initial map load, the `map.fitBounds(markerBounds)` is done before
  //  the map resizes, and so it will not be valid.  So we will mark when
  //  the autozoom (e.g., map.fitBounds()) should have set the size, using
  //  the `isAutoZoomed` member variable, and then when the map resizes, 
  //  if this variable is still set, we will re-auto-zoom.
  //  Any mouse down, or zooming will unset this variable - note however,
  //  that the zoomstart/zoomend callbacks will happen during data loading,
  //  so this unmarking happens even when we wouldnt like
  //  #TODO: Figure out how to set `this.isAutoZoomed` to false only when
  //         the user changes zoom, or pans.
  //         (but currently this mechanism works good-enough to get auto
  //          zoom to work when we first show map in InterSpec)
  this.isAutoZoomed = false;
  
  this.map.on('mousedown', (e) => { self.isAutoZoomed = false; });
  this.map.on('zoomstart', (e) => { self.isAutoZoomed = false; });
  this.map.on('zoomend',   (e) => { self.isAutoZoomed = false; });

  this.map.on('resize', (e) => { self.handleResize(); });

  // For debugging; it would be nice to display this somewhere maybe
  this.map.on('click', (e) => {
    console.log(`clicked at lat/lng: (${e.latlng.lat}, ${e.latlng.lng})` );
  });
}//LeafletRadMap constructor


LeafletRadMap.prototype.handleResize = function() {
  //console.log( 'handleResize, isAutoZoomed=', this.isAutoZoomed );

  // FIXME: this.isAutoZoomed will be set to false during normal data
  //        loading (unfortunately), not just when the the user changes
  //        things.

  if( !this.isAutoZoomed || !this.markerLayer )
    return;

  const markerBounds = this.markerLayer.getBounds();
  if( !markerBounds )
    return;

  const diagMeters = markerBounds.getNorthWest().distanceTo( markerBounds.getSouthEast() );
  
  this.map.invalidateSize();
  if( diagMeters > 10 ){
    this.map.fitBounds( markerBounds );
    //console.log( 'Set bounds to:', markerBounds );
  }else {
    this.map.setView( markerBounds.getCenter(), 18 );
  }
};//LeafletRadMap.prototype.handleResize = function() 


LeafletRadMap.prototype.WtEmit = function(elem, event) {
  if (!window.Wt) {
    console.warn('Wt not found! Canceling "' + event.name + '" emit function...');
    return;
  }
  
  // To support ES5 syntax in IE11, we replace spread operator with this
  let args = Array.prototype.slice.call(arguments, SpectrumChartD3.prototype.WtEmit.length);
  
  Wt.emit.apply(Wt, [elem, event].concat(args));
}


LeafletRadMap.prototype.updateMarkerLegend = function(){
  const self = this;
  
  let have_fore = false, have_back = false, have_second = false, have_other = false;
  for( const m of self.markers ){
    if( m.displayType === 0 )
      have_fore = true;
    else if( m.displayType === 1 )
      have_second = true;
    else if( m.displayType === 2 )
      have_back = true;
    else
      have_other = true;
    if( have_fore && have_second && have_back && have_other )
      break;
  }//for( const m of self.markers )
  
  self.markerLegendDiv.innerHTML = '';
  
  // if we have only one type of marker, we dont need a legend
  if( (0 + have_fore + have_back + have_second + have_other) <= 1 )
    return;
  
  if( have_fore )
    self.markerLegendDiv.innerHTML += '<div class="LegEntry"><div class="RadMapMarkerInner DispAsFore"></div><div>'
                                      + self.options.foregroundTxt + '</div></div>';
  if( have_back )
    self.markerLegendDiv.innerHTML += '<div class="LegEntry"><div class="RadMapMarkerInner DispAsBack"></div><div>'
                                      + self.options.backgroundTxt + '</div></div>';
  if( have_second )
    self.markerLegendDiv.innerHTML += '<div class="LegEntry"><div class="RadMapMarkerInner DispAsSecond"></div><div>'
                                      + self.options.secondaryTxt + '</div></div>';
  if( have_other )
    self.markerLegendDiv.innerHTML += '<div class="LegEntry"><div class="RadMapMarkerInner"></div><div>Not Displayed</div></div>';
  
  //self.markerLegendDiv.style.display = 'none';
}//updateMarkerLegend(...)

LeafletRadMap.prototype.handleUserDrawingUpdate = function(){
  const self = this;

  const wantedSamples = this.getSelectedSampleNumbers();
  
  if( wantedSamples.length === 0 ){
    self.btnsDiv.innerHTML = '<div class="NoSamplesSel"><p>No measurements selected</p></div>';
    self._appendCopyMapBtns();
    return;
  }

  function sendToInterSpec(evt, spectype, samples){
    // We attach the sigmal to InterSpec object, and not the actual map opbject because there is some issue with the
    //  signal being exposed when in an AuxWindow.
    self.WtEmit( document.querySelector('.InterSpec').id, {name: 'loadSamples', eventObject: evt}, samples, spectype );
    // If we want to go back to attaching to the `LeafletRadMap` instance, we would use:
    //  self.WtEmit( self.parent.id, {name: 'loadSamples', eventObject: evt}, samples, spectype );
  }

  self.btnsDiv.innerHTML = '';

  const msg = document.createElement('div');
  if( window.innerWidth < 675 )
    msg.innerHTML = self.options.loadTxtShort + ":";
  else
    msg.innerHTML = self.options.loadTxt + " " + wantedSamples.length + " " + self.options.measurementsAsTxt + ":";
  self.btnsDiv.appendChild( msg );

  const foreground = document.createElement('button');
  foreground.innerHTML = self.options.foregroundTxt;
  foreground.classList.add("Wt-btn");
  foreground.classList.add("with-label");
  foreground.addEventListener("click", function(evt){
    sendToInterSpec( evt, 'foreground', wantedSamples );
  });

  const background = document.createElement('button');
  background.innerHTML = self.options.backgroundTxt;
  background.classList.add("Wt-btn");
  background.classList.add("with-label");
  background.addEventListener("click", function(evt){
    sendToInterSpec( evt, 'background', wantedSamples );
  });

  const secondary = document.createElement('button');
  secondary.innerHTML = self.options.secondaryTxt;
  secondary.classList.add("Wt-btn");
  secondary.classList.add("with-label");
  secondary.addEventListener("click", function(evt){
    sendToInterSpec( evt, 'secondary', wantedSamples );
  });

  self.btnsDiv.appendChild( foreground );
  self.btnsDiv.appendChild( background );
  self.btnsDiv.appendChild( secondary );

  self._appendCopyMapBtns();
};//LeafletRadMap.prototype.handleUserDrawingUpdate


/** Returns the sample numbers of measurements currently within the drawn on polygons; or if no shapes
 * are drawn, returns all sample numbers.
 */
LeafletRadMap.prototype.getSelectedSampleNumbers = function(){
  const self = this;

  const layers = self.map.pm.getGeomanDrawLayers();
  const group = L.featureGroup();
  layers.forEach((layer)=>{
    group.addLayer(layer);
  });
    
  const shapeGJ = group.toGeoJSON();

  // Each entry in this `polygons` array will itself be an array of polygons.
  //  The first polygon of each entry will be the outline polygon.  The remaining 
  //  polygons are the exclude polygons
  //  polygons =[ [[outline poly],[exclude poly], [exclude poly]], [[outline poly]], ... ]
  const polygons = [];

  // Loop over each shape the user has drawn - note that if a user splits a shape using the cut tool, it will still only be one entry in shapeGJ.features for both shown shapes
  for( const feature of shapeGJ.features ) {
    if( feature.geometry.type == "Polygon" ){
      // If the user has draw a polygon, and optionally cut an area out, that didnt actually split the polygon, we will go here
      polygons.push( feature.geometry.coordinates );
    }else if( feature.geometry.type == "MultiPolygon" ){
      // If the user used the cut tool to split the polygon, we will go here
      for (const coords of feature.geometry.coordinates ) {
        polygons.push( coords );
      }
    } if( feature.geometry.type == "Point" ){
      // We get here if user added Text to the map
    } else {
      console.log( 'Unhandled geometry: ', feature );
    }
  }//for( const feature of shapeGJ.features )

  let contained_samples = [];
  if( !self.markerLayer )
    return contained_samples;

  // Now loop over all the markers, and find the ones in the polygons we want
  self.markerLayer.eachLayer(layer => {
    const latLng = layer.getLatLng();
    const coords = [latLng.lng,latLng.lat];

    if( polygons.length < 1 ){
      contained_samples.push( layer.sampleNumber );
      return;
    }
      

    let isInAnyOutline = false, isInAnyExclude = false;
    
    for( const coordinates of polygons ) {
      // I think 'feature' is an object the user has drawn.  Should be polygon, but may have excluded regions
      
      //console.log( 'feature.geometry.coordinates', feature.geometry.coordinates );
      const outline_polys = [], exclude_polys = [];
      if( coordinates.length > 0 )
        outline_polys.push( coordinates[0] );
        
      for( let i = 1; i < coordinates.length; ++i )
        exclude_polys.push( coordinates[i] );

      for (const poly of outline_polys ) {
        const isIn = d3.polygonContains( poly, coords );
        if( isIn ){
          isInAnyOutline = true;
          //console.log( 'Is in outline, ', layer.sampleNumber );
        }
      }

      for (const poly of exclude_polys ) {
        const isIn = d3.polygonContains( poly, coords );
        if( isIn ){
          isInAnyExclude = true;
          //console.log( 'Is in exclude' );
        }
      }
    }//for( const coordinates of polygons )

    if( isInAnyOutline && !isInAnyExclude )
      contained_samples.push( layer.sampleNumber );
  } );//self.markerLayer.eachLayer(layer => {

  return contained_samples;
};//LeafletRadMap.prototype.getSelectedSampleNumbers



/**  Takes in gradient stops (ex. {0.4: 'blue', 0.7: 'lime', 1.0: 'red'}), and returns an array
 with 256 entries that has the rgba components that map from 0 to 1. 
*/
LeafletRadMap.prototype.defineGradientColors = function(gradientStops, minOpacity){
  // This uses same methodology as https://github.com/mourner/simpleheat to create a gradient
  const self = this;

  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const gradient = ctx.createLinearGradient(0, 0, 0, 256);
  canvas.width = 1;
  canvas.height = 256;
  for( let i in gradientStops ) 
    gradient.addColorStop( i, gradientStops[i] );
  ctx.fillStyle = gradient;
  ctx.fillRect(0, 0, 1, 256);

  const gradData = ctx.getImageData(0, 0, 1, 256).data;
    
  let colors = []
  for( let i = 0; i < 256; ++i )
    colors.push( [gradData[4*i+0],gradData[4*i+1], gradData[4*i+2], minOpacity + ((i+1)/256.0)*(1-minOpacity)] );
    
  return colors;
};//LeafletRadMap.prototype.defineGradientColors


/** Returns a CSS string like 'rgba(124,232,5,0.32)' cooresponding to location on gradient of the input number.
@param frac Input number, in range of 0 to 1
@param grad_colors The array returned from defineGradientColors
*/
LeafletRadMap.prototype.gradColorForFrac = function( frac, grad_colors ){
  const index = Math.min( Math.max( Math.round(frac*grad_colors.length), 0), grad_colors.length - 1 );
  return "rgba(" + grad_colors[index][0] + "," + grad_colors[index][1] + "," + grad_colors[index][2] + "," + grad_colors[index][3] + ")";
};//LeafletRadMap.prototype.gradColorForFrac


LeafletRadMap.prototype.displayTypeClass = function( disp ) {
  if( disp === 0 )
    return "DispAsFore";
  if( disp === 1 )
    return "DispAsSecond";
  if( disp === 2 )
    return "DispAsBack";
  return "";
}//LeafletRadMap.prototype.setData(...)



LeafletRadMap.prototype.setData = function( data, dont_update_zoom ) {
  this.data = data;
  this.refresh( dont_update_zoom );
}//LeafletRadMap.prototype.setData(...)


LeafletRadMap.prototype.refresh = function( dont_update_zoom ){
  const self = this;

  const ndata = this.data && this.data.samples ? this.data.samples.length : 0;
 
  //Remove previous markers
  if( this.markerLayer ){
    this.markerLayer.remove();
    this.markerLayer = null;
  }

  // Remove all user drawn shapes
  this.map.eachLayer( function(layer){
    if (layer._path != null) {
      layer.remove();
    }
  });

  if( this.heatmap ){
    this.heatmap.remove();
    this.heatmap = null;
  }
  
  if( this.heatMapButton ){
    this.heatMapButton.remove();
    this.heatmap = null;
  }


  // We will hide drawing controls if we only have 0 or 1 unique GPS coordinates.
  //  Note that JS sets wont work if you insert arrays to them, as it will test
  //  for unique object references in that case, instead of array equality, so
  //  we will use seperate Sets for lat and lon
  const unique_lat = new Set(), unique_lon = new Set();

  const measurements = this.data ? this.data.samples : [];


  this.markers = [];
  this.heatMapData = [];
  this.maxGammaCps = 0;
  this.maxNeutCps = 0;
  this.minGammaCps = (measurements && (measurements.length > 0)) ? 1.0E15 : 0;
  this.minNeutronCps = (measurements && (measurements.length > 0)) ? 1.0E15 : 0;


  // First loop through to get the max gamma and neutron count rates, so we can set marker colors according to gradient
  for( sample of measurements )
  {
    if( isNaN(sample.gps[0]) || isNaN(sample.gps[1]) )
      continue;

    unique_lat.add( sample.gps[0] );
    unique_lon.add( sample.gps[1] );

    if( sample.nGDet > 0 )
    {
      const gamma_cps = sample.gSum / (sample.gLT > 0 ? sample.gLT : 1);
      self.maxGammaCps = Math.max( self.maxGammaCps, gamma_cps );
      self.minGammaCps = Math.min( self.minGammaCps, gamma_cps );
    }//if( sample.nGDet > 0 )

    if( sample.nNDet > 0 )
    {
      const neut_cps = sample.nSum / (sample.nRT > 0 ? sample.nRT : 1);
      self.maxNeutCps = Math.max( self.maxNeutCps, neut_cps );
      self.minNeutronCps = Math.min( self.minNeutronCps, neut_cps );
    }//if( sample.nNDet > 0 )
  }//for( sample of measurements )

  if( self.minGammaCps >= self.maxGammaCps )
    self.minGammaCps = 0;
  self.gammaCpsRange = (self.maxGammaCps > self.minGammaCps) ? (self.maxGammaCps - self.minGammaCps) : 1;

  if( self.minNeutronCps >= self.maxGammaCps )
    self.minNeutronCps = 0;
  self.neutCpsRange = (self.maxNeutCps > self.minNeutronCps) ? (self.maxNeutCps - self.minNeutronCps) : 1;


  for( sample of measurements )
  {
    const lat_lon = sample.gps;

    if( isNaN(lat_lon[0]) || isNaN(lat_lon[1]) )
      continue;

    let gamma_cps = -1, neut_cps = -1, rel_gamma_cps = 0;
    if( sample.nGDet > 0 )
    {
      const gamma_counts =  sample.gSum;
      const gamma_lt =  sample.gLT;
      gamma_cps = gamma_counts / (gamma_lt > 0 ? gamma_lt : 1);
      rel_gamma_rate = (gamma_cps - self.minGammaCps) / self.gammaCpsRange;
    }//if( sample.nGDet > 0 )

    if( sample.nNDet > 0 )
    {
      const neut_counts =  sample.nSum;
      const neut_lt =  sample.nRT; //Using neutron real-time, assuming no dead time
      neut_cps = neut_counts / (neut_lt > 0 ? neut_lt : 1);
    }//if( sample.nNDet > 0 )


    let marker = new L.marker( [lat_lon[0], lat_lon[1]], { 
      icon: L.divIcon({
        className: 'RadMapMarker',
        iconSize: [8, 8],
        html: '<div style="background-color: '
                + self.gradColorForFrac(rel_gamma_rate, self.grad_colors)
                + '" class="RadMapMarkerInner '
                + self.displayTypeClass(sample.disp)
                + '"/>'
      }),
      pmIgnore: true,  // without this geoman could delete the marker
      riseOnHover: true,
      draggable: false
    } );

    function onMarkerMouseover(m){
      let html = "";
      
      html += "<table class=\"MapMousePopup\">";
      html += `<tr><td>${self.options.gammaTxt} ${self.options.cpsText}</td><td>${m.gammaCps.toFixed(2)}</td></tr>`;
      if( typeof m.neutronCps === "number" )
        html += `<tr><td>Neutron ${self.options.cpsText}</td><td>${m.neutronCps.toFixed(4)}</td></tr>`;
      html += `<tr><td>${self.options.realTimeTxt}</td><td>${m.realTime.toFixed(3)}</td></tr>`;
      html += `<tr><td>${self.options.liveTimeTxt}</td><td>${m.gammaLiveTime.toFixed(3)}</td></tr>`;
      if( m.startTime )
      {
        const datetime = m.startTime.toISOString().slice(0, -1).split("T");
        const date = datetime.length ? datetime[0] : "";
        const time = datetime.length > 1 ? datetime[1] : "";
        html += `<tr><td>Date</td><td>${date}</td></tr>`;
        html += `<tr><td>Time</td><td>${time}</td></tr>`;
      }
      
      if( m.displayType >= 0 && m.displayType < 3 ){
        html += '<tr><td>' + self.options.displayedAsTxt + '</td><td>';
        if( m.displayType == 0 ) html += self.options.foregroundTxt;
        else if( m.displayType == 1 ) html += self.options.secondaryTxt;
        else if( m.displayType == 2 ) html += self.options.backgroundTxt;
        html += '</td></tr>';
      }

      if( m.sampleNumbers )
      {
        const val = (m.sampleNumbers.length > 14) ? (m.sampleNumbers.substr(0,11) + "...") : m.sampleNumbers;
        html += `<tr><td>Samples</td><td>${val}</td></tr>`;
      }else
      {
        html += `<tr><td>Sample</td><td>${m.sampleNumber}</td></tr>`;
      }

      html += `</table>`;

      m.setPopupContent( html );
      m.openPopup( m.getLatLng() );
    };//onMarkerMouseover function


    marker.bindPopup('').openPopup();
    marker.on('mouseover', function(){ onMarkerMouseover(marker); });

    marker.sampleNumber = sample.sample;
    if( sample.samples && ((typeof sample.samples) === "string") )
      marker.sampleNumbers = sample.samples;
    marker.realTime = sample.rt;
    marker.rateIntensity = rel_gamma_rate;
    marker.displayType = (typeof sample.disp === "number") ? sample.disp : 3;

    if( sample.nGDet > 0 )
    {
      marker.gammaCps = gamma_cps;
      marker.gammaLiveTime = sample.gLT;

      this.heatMapData.push( [lat_lon[0], lat_lon[1], rel_gamma_rate] );
    }//if( sample.nGDet > 0 )

    if( sample.nNDet > 0 )
    {
      marker.neutronCps = neut_cps;
    }//if( sample.nNDet > 0 )


    if( typeof sample.timeOffset === "number" )
      marker.startTime = new Date(sample.timeOffset + this.data.startTime);
    
    self.markers.push( marker );
  }//for( sample of data.samples )
  
  // Dont show drawing controls for 0 or 1 items.
  if( ((unique_lat.size > 1) || (unique_lon.size > 1)) != this.map.pm.controlsVisible() ){
    this.map.pm.toggleControls();
  }
    

  // Refer to https://github.com/Leaflet/Leaflet.heat for options
  const heatLayerOpts = {
    max: 1.0,
    radius: 35,
    blur: 25,
    minOpacity: self.options.minMarkerOpacity,
    gradient: this.markerGradientStops,
  };
  
  
  if( (unique_lat.size > 1) || (unique_lon.size > 1) ){
    self.heatMapButton = L.easyButton({
      leafletClasses: true,
  
      states: [{
        stateName: 'with-heatmap',
        icon:      'RemoveHeatmapBtn',
        title:     'Remove heatmap',
        onClick: function(btn, map) {
          if( self.heatmap )
            map.removeLayer(self.heatmap);
          self.options.showHeatMap = false;
          btn.state('no-heatmap'); 
          }
        }, {
          stateName: 'no-heatmap',
          icon:      'ShowHeatmapBtn',
          title:     'Show heatmap',
          onClick: function(btn, map) {
            if( self.heatMapData )
              self.heatmap = L.heatLayer(self.heatMapData, heatLayerOpts).addTo(self.map);
            self.options.showHeatMap = true;
            btn.state('with-heatmap'); 
          }
        }]
    });
  
    self.heatMapButton.addTo( self.map );

    if( self.options.showHeatMap ) {
      self.heatmap = L.heatLayer(self.heatMapData, heatLayerOpts).addTo(self.map);
      self.heatMapButton.state('with-heatmap');
    }else{
      self.heatMapButton.state('no-heatmap');
    }
  }//if( unique_coord.size > 1 ){

  if( self.markers.length ){
    self.map.invalidateSize();
    
    // See the 'iconCreateFunction' option form https://github.com/Leaflet/Leaflet.markercluster#other-options
    self.markerLayer = L.markerClusterGroup({
      removeOutsideVisibleBounds: false,
      maxClusterRadius: 5,
      singleMarkerMode: false,
      iconCreateFunction: function(cluster) {
        let maxClusterIntensity = 0, displayType = 3;
        for( m of cluster.getAllChildMarkers() ){
          maxClusterIntensity = Math.max( maxClusterIntensity, m.rateIntensity );
          displayType = Math.min( displayType, m.displayType )
        }
        
        return L.divIcon({
          className: 'RadMapMarker',
          iconSize: [22, 22],
          html: '<div style="background-color: '
                  + self.gradColorForFrac(maxClusterIntensity, self.grad_colors)
                  + '" class="RadMapMarkerMulti ' + self.displayTypeClass(displayType) + '">'
                  + '<div>' + cluster.getChildCount() + '</div>'
                + '</div>'
        });
      }
    });

    for( const m of self.markers )
      self.markerLayer.addLayer( m );
    self.map.addLayer(self.markerLayer);

    if( !dont_update_zoom ){
      const markerBounds = self.markerLayer.getBounds();
      const diagMeters = markerBounds.getNorthWest().distanceTo( markerBounds.getSouthEast() );
      
      if( diagMeters > 10 ){
        self.map.fitBounds( markerBounds );
        this.isAutoZoomed = true;
        //console.log( "Setting isAutoZoomed to true;");
      }else {
        this.isAutoZoomed = false;
        self.map.setView( markerBounds.getCenter(), 18 );
      }
    }//if( !dont_update_zoom )
  }else if( !dont_update_zoom ){
    this.isAutoZoomed = false;
    self.map.setView([37.67640130843344, -121.70667619429008], 15); 
  }

  // Show the bottom div whenever there are markers (for copy-URL buttons);
  //  gradient legend only shows when there are multiple markers.
  if( self.markers.length >= 1 ){
    self.bottomDiv.style.display = null;
  }else{
    self.bottomDiv.style.display = 'none';
  }

  if( (self.markers.length) > 1 && self.markerGradientStops ){
    self.gradientLegendDiv.style.display = null;

    let maxcps = self.maxGammaCps;
    let mincps = self.minGammaCps;
    if( maxcps < 1 )
      maxcps = maxcps.toFixed(3);
    else if( maxcps < 10 )
      maxcps = maxcps.toFixed(1);
    else
      maxcps = Math.round(maxcps)

    if( mincps < 1 )
      mincps = mincps.toFixed(3);
    else if( mincps < 10 )
      mincps = mincps.toFixed(1);
    else
      mincps = Math.round(mincps)


    self.legendUpperCps.innerHTML = maxcps + "<br />" + self.options.cpsText;
    self.legendLowerCps.innerHTML = mincps + "<br />" + self.options.cpsText;

    const width = self.legendCanvas.clientWidth;
    const height = self.legendCanvas.clientHeight;

    const ctx = self.legendCanvas.getContext('2d');
    const gradient = ctx.createLinearGradient(0, 0, 0, height);
    self.legendCanvas.width = width;
    self.legendCanvas.height = height;
    for( let i in self.markerGradientStops )
      gradient.addColorStop( i, self.markerGradientStops[i] );

    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, width, height);
  }else{
    self.gradientLegendDiv.style.display = 'none';
  }

  self.handleUserDrawingUpdate();
  self.updateMarkerLegend();
}//LeafletRadMap.prototype.refresh()


/** Populate the copyBtnsDiv with "Copy Google Maps URL" and "Copy Bing Maps URL" buttons.
 Only creates the buttons once; subsequent calls just show/hide the div.
 */
LeafletRadMap.prototype._appendCopyMapBtns = function() {
  const self = this;
  const hasData = self.data && self.data.samples && self.data.samples.length > 0;

  if( self.copyBtnsDiv.childElementCount === 0 && hasData )
  {
    const googleBtn = document.createElement( 'button' );
    googleBtn.innerHTML = self.options.copyGoogleMapsTxt || 'Copy Google Maps URL';
    googleBtn.classList.add( 'Wt-btn', 'LinkBtn' );
    googleBtn.addEventListener( 'click', function(){ self.copyGoogleMapsUrl(); } );

    const bingBtn = document.createElement( 'button' );
    bingBtn.innerHTML = self.options.copyBingMapsTxt || 'Copy Bing Maps URL';
    bingBtn.classList.add( 'Wt-btn', 'LinkBtn' );
    bingBtn.addEventListener( 'click', function(){ self.copyBingMapsUrl(); } );

    self.copyBtnsDiv.appendChild( googleBtn );
    self.copyBtnsDiv.appendChild( bingBtn );
  }

  self.copyBtnsDiv.style.display = hasData ? null : 'none';
};//LeafletRadMap.prototype._appendCopyMapBtns


/** Collect unique valid GPS points from data samples, preserving order. */
LeafletRadMap.prototype._getGpsPoints = function() {
  const samples = this.data ? this.data.samples : [];
  const points = [];
  const seen = new Set();

  for( let i = 0; i < samples.length; i++ )
  {
    const s = samples[i];
    if( !s.gps || isNaN(s.gps[0]) || isNaN(s.gps[1]) )
      continue;

    const key = s.gps[0].toFixed(6) + ',' + s.gps[1].toFixed(6);
    if( !seen.has(key) )
    {
      seen.add( key );
      points.push( [s.gps[0], s.gps[1]] );
    }
  }

  return points;
};//LeafletRadMap.prototype._getGpsPoints


/** Subsample an array of points to at most maxPoints, evenly spaced. */
LeafletRadMap.prototype._subsamplePoints = function( points, maxPoints ) {
  if( points.length <= maxPoints )
    return points;

  const result = [];
  const step = (points.length - 1) / (maxPoints - 1);

  for( let i = 0; i < maxPoints; i++ )
  {
    const idx = Math.round( i * step );
    result.push( points[idx] );
  }

  return result;
};//LeafletRadMap.prototype._subsamplePoints


/** Copy text to clipboard, showing a notification on success or failure. */
LeafletRadMap.prototype._copyToClipboard = function( text ) {
  const self = this;
  const successMsg = 'showMsg-info-' + (self.options.copyUrlSuccessTxt || 'Map URL copied to clipboard.');
  const failMsg = 'showMsg-error-' + (self.options.copyUrlFailTxt || 'Failed to copy map URL.');

  function notify( msg ) {
    if( window.Wt )
      Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'}, msg );
  }

  if( !navigator.clipboard )
  {
    let textArea = document.createElement( "textarea" );
    textArea.value = text;
    textArea.style.position = "fixed";
    textArea.style.top = "0";
    textArea.style.left = "0";
    textArea.style.width = "2em";
    textArea.style.height = "2em";
    textArea.style.opacity = "0";
    document.body.appendChild( textArea );
    textArea.focus();
    textArea.select();

    let successful = false;
    try
    {
      successful = document.execCommand( 'copy' );
    }catch( ex )
    {
      successful = false;
    }

    document.body.removeChild( textArea );
    notify( successful ? successMsg : failMsg );
    return;
  }//if( !navigator.clipboard )

  navigator.clipboard.writeText( text ).then( function() {
    notify( successMsg );
  }, function( err ) {
    console.warn( 'LeafletRadMap: Failed to copy to clipboard: ', err );
    notify( failMsg );
  });
};//LeafletRadMap.prototype._copyToClipboard


/** Build and copy a Google Maps URL from the current measurement GPS data. */
LeafletRadMap.prototype.copyGoogleMapsUrl = function() {
  const self = this;
  const points = self._getGpsPoints();

  if( points.length === 0 )
  {
    if( window.Wt )
      Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'},
        'showMsg-info-' + (self.options.copyUrlNoGpsTxt || 'No GPS data.') );
    return;
  }

  let url;
  if( points.length === 1 )
  {
    // Single point: center on location
    url = 'https://www.google.com/maps/@'
      + points[0][0].toFixed(6) + ',' + points[0][1].toFixed(6) + ',15z';
  }else
  {
    // Multiple points: directions mode with up to 10 waypoints
    const sampled = self._subsamplePoints( points, 10 );
    url = 'https://www.google.com/maps/dir';
    for( let i = 0; i < sampled.length; i++ )
      url += '/' + sampled[i][0].toFixed(6) + ',' + sampled[i][1].toFixed(6);
  }

  self._copyToClipboard( url );
};//LeafletRadMap.prototype.copyGoogleMapsUrl


/** Build and copy a Bing Maps URL from the current measurement GPS data.
 Uses the 'rtp' (route) parameter with 'pos' waypoints, since the 'sp=point'
 pushpin format is currently broken on the Bing Maps consumer website.
 */
LeafletRadMap.prototype.copyBingMapsUrl = function() {
  const self = this;
  const points = self._getGpsPoints();

  if( points.length === 0 )
  {
    if( window.Wt )
      Wt.emit( document.querySelector('.specviewer').id, {name:'miscSignal'},
        'showMsg-info-' + (self.options.copyUrlNoGpsTxt || 'No GPS data.') );
    return;
  }

  let url;
  if( points.length === 1 )
  {
    // Single point: center on location
    url = 'https://bing.com/maps/default.aspx?cp='
      + points[0][0].toFixed(6) + '~' + points[0][1].toFixed(6) + '&lvl=15';
  }else
  {
    // Multiple points: route waypoints, up to 25
    const sampled = self._subsamplePoints( points, 25 );
    let rtp = '';
    for( let i = 0; i < sampled.length; i++ )
    {
      if( i > 0 ) rtp += '~';
      rtp += 'pos.' + sampled[i][0].toFixed(6)
        + '_' + sampled[i][1].toFixed(6) + '_' + (i + 1);
    }
    url = 'https://bing.com/maps/default.aspx?rtp=' + rtp;
  }

  self._copyToClipboard( url );
};//LeafletRadMap.prototype.copyBingMapsUrl
