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

  this.legendDiv = document.createElement("div");
  this.legendDiv.classList.add("LeafletRadMapLegend");
  this.parent.appendChild( this.legendDiv );

  this.btnsDiv = document.createElement("div");
  this.btnsDiv.classList.add("LeafletRadMapBtnsDiv");
  this.parent.appendChild( this.btnsDiv );



  this.options = options ? options : {};
  if( typeof this.options.heatMap !== "boolean" )
    this.options.heatMap = true;
  if( typeof this.options.minMarkerOpacity !== "number" )
    this.options.minMarkerOpacity = 0.1;
  if( typeof this.options.showHeatMap !== "boolean" )
    this.options.showHeatMap = false;

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
  this.legendDiv.appendChild( this.legendLowerCps );

  this.legendCanvas = document.createElement('canvas');
  this.legendCanvas.classList.add("LegCanvas");
  this.legendDiv.appendChild( this.legendCanvas );

  this.legendUpperCps = document.createElement("div");
  this.legendUpperCps.classList.add("LegUpVal");
  this.legendDiv.appendChild( this.legendUpperCps );

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
}//LeafletRadMap constructor


LeafletRadMap.prototype.handleUserDrawingUpdate = function(){
  const self = this;

  const wantedSamples = this.getSelectedSampleNumbers();
  
  if( wantedSamples.length === 0 ){
    self.btnsDiv.innerHTML = '<div class="NoSamplesSel">No measurements selected</div>';
    return;
  }

  function sendToInterSpec(spectype, samples){
    console.log( 'Will send to InterSpec type:', spectype, ', samples:', samples );
  }

  self.btnsDiv.innerHTML = '';

  const msg = document.createElement('div');
  msg.innerHTML = "Load " + wantedSamples.length + " measurements as:";
  self.btnsDiv.appendChild( msg );

  const foreground = document.createElement('button');
  foreground.innerHTML = "Foreground";
  foreground.classList.add("Wt-btn");
  foreground.classList.add("with-label");
  foreground.addEventListener("click", function(){
    sendToInterSpec( 'foreground', wantedSamples );
  });

  const background = document.createElement('button');
  background.innerHTML = "Background";
  background.classList.add("Wt-btn");
  background.classList.add("with-label");
  background.addEventListener("click", function(){
    sendToInterSpec( 'background', wantedSamples );
  });

  const secondary = document.createElement('button');
  secondary.innerHTML = "Secondary";
  secondary.classList.add("Wt-btn");
  secondary.classList.add("with-label");
  secondary.addEventListener("click", function(){
    sendToInterSpec( 'secondary', wantedSamples );
  });

  self.btnsDiv.appendChild( foreground );
  self.btnsDiv.appendChild( background );
  self.btnsDiv.appendChild( secondary );
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
        //console.log('handleUserDrawingUpdate toGeoJSON feature poly', poly );
        const isIn = d3.polygonContains( poly, coords );
        if( isIn ){
          isInAnyOutline = true;
          //console.log( 'Is in outline, ', layer.sampleNumber );
        }
      }

      for (const poly of exclude_polys ) {
        //console.log('handleUserDrawingUpdate toGeoJSON feature poly', poly );
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


LeafletRadMap.prototype.setData = function( data ) {
  this.data = data;
  this.refresh();
}//LeafletRadMap.prototype.setData(...)


LeafletRadMap.prototype.refresh = function(){
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
  const unique_coord = new Set();

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

    unique_coord.add( sample.gps );

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
        className: 'MapCircleMarker',
        iconSize: [8, 8],
        html: '<div style="background-color: ' + self.gradColorForFrac(rel_gamma_rate, self.grad_colors) + '" class="MapCircleMarkerInner"/>'
      }),
      pmIgnore: true,  // without this geoman could delete the marker
      riseOnHover: true,
      draggable: false
    } );

    function onMarkerMouseover(m){
      let html = "";
      
      html += "<table>";
      html += `<tr><td>Gamma CPS</td><td>${m.gammaCps.toFixed(2)}</td></tr>`;
      if( typeof m.neutronCps === "number" )
        html += `<tr><td>Neutron CPS</td><td>${m.neutronCps.toFixed(4)}</td></tr>`;
      html += `<tr><td>Real Time</td><td>${m.realTime.toFixed(3)}</td></tr>`;
      html += `<tr><td>Live Time</td><td>${m.gammaLiveTime.toFixed(3)}</td></tr>`;
      if( m.startTime )
        html += `<tr><td>Time</td><td>${m.startTime.toISOString().slice(0, -1)}</td></tr>`;
      html += `<tr><td>Sample</td><td>${m.sampleNumber}</td></tr>`;
      html += `</table>`;

      m.setPopupContent( html );
      m.openPopup( m.getLatLng() );
    };//onMarkerMouseover function


    marker.bindPopup('').openPopup();
    marker.on('mouseover', function(){ onMarkerMouseover(marker); });

    marker.sampleNumber = sample.sample;
    marker.realTime = sample.rt;
    marker.rateIntensity = rel_gamma_rate;

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
  if( (unique_coord.size > 1) != this.map.pm.controlsVisible() ){
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
  
  
  if( unique_coord.size > 1 ){
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
    // See the 'iconCreateFunction' option form https://github.com/Leaflet/Leaflet.markercluster#other-options
    self.markerLayer = L.markerClusterGroup({
      removeOutsideVisibleBounds: false,
      maxClusterRadius: 5,
      singleMarkerMode: false,
      iconCreateFunction: function(cluster) {
        let maxClusterIntensity = 0;
        for( m of cluster.getAllChildMarkers() )
          maxClusterIntensity = Math.max( maxClusterIntensity, m.rateIntensity );
    
        return L.divIcon({
          className: 'MapCircleMarker',
          iconSize: [20, 20],
          html: '<div style="background-color: ' + self.gradColorForFrac(maxClusterIntensity, self.grad_colors) + '" class="MapCircleMarkerMulti">' + cluster.getChildCount() + '</div>'
        });
      }
    });

    for( const m of self.markers )
      self.markerLayer.addLayer( m );
    self.map.addLayer(self.markerLayer);

    const markerBounds = self.markerLayer.getBounds();
    const diagMeters = markerBounds.getNorthWest().distanceTo( markerBounds.getSouthEast() );
    if( diagMeters > 10 ){
      self.map.fitBounds( markerBounds );
    }else {
      self.map.setView( markerBounds.getCenter(), 18 )
    }
  }else{
    self.map.setView([37.67640130843344, -121.70667619429008], 15); 
  }

  if( (self.markers.length) > 1 && self.markerGradientStops ){
    self.legendDiv.style.display = null;
    self.btnsDiv.style.display = null;

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
      

    self.legendUpperCps.innerHTML = maxcps + "<br />CPS";
    self.legendLowerCps.innerHTML = mincps + "<br />CPS";

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
    self.legendDiv.style.display = 'none';
    self.btnsDiv.style.display = 'none';
  }

  self.handleUserDrawingUpdate();
}//LeafletRadMap.prototype.refresh()