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

D3TimeChart = function(elem, options) {
  var self = this;
  this.chart = typeof elem === 'string' ? document.getElementById(elem) : elem; 

  this.options = options || {}; 
  
  if( (typeof this.options.xtitle) !== 'string' ) this.options.xlabel = "Time (s)";
  if( (typeof this.options.y1title) !== 'string' ) this.options.y1title = "Gamma CPS";
  if( (typeof this.options.y2title) !== 'string' ) this.options.y2title = "N CPS";
  if( (typeof this.options.compactXAxis) !== 'boolean' ) this.options.compactXAxis = false;
  if( (typeof this.options.gridx) !== 'boolean' ) this.options.gridx = false;
  if( (typeof this.options.gridy) !== 'boolean' ) this.options.gridy = false;
  if( (typeof this.options.chartLineWidth) !== 'number' ) this.options.chartLineWidth = 1;
  

}//


/** Function to help emit callbacks back to C++
 
 @param elem: The element can be a DOM element, or the object ID of a WObject
 @param event: must be an object which indicates also the JavaScript event and event target
 @param args: array of args to pass into the Wt function
 */
D3TimeChart.prototype.WtEmit = function(elem, event) {
  if (!window.Wt) {
    console.warn('Wt not found! Canceling "' + event.name + '" emit function...');
    return;
  }

  //console.log( 'Emitting Wt event "' + ((event && event.name) ? event.name : 'null') + '", with ' + SpectrumChartD3.prototype.WtEmit.length + " arguments");

  // To support ES5 syntax in IE11, we replace spread operator with this
  var args = Array.prototype.slice.call(arguments, SpectrumChartD3.prototype.WtEmit.length);

  // Emit the function to Wt
  // Wt.emit( elem, event, ...args);  // ES6 syntax
  Wt.emit.apply(Wt, [elem, event].concat(args));
}

D3TimeChart.prototype.handleResize = function(){
  // This function is called when the Wt layout manager resizes the parent <div> element
  // Need to redraw everything
  // ...
  
  // Make sure to update the C++ code of the changed plotting size.
    
}

D3TimeChart.prototype.setData = function( data ){
  //See the c++ function D3TimeChart::setData()
}

D3TimeChart.prototype.setHighlightRegions = function( regions ){
  //See the c++ function D3TimeChart::setHighlightRegionsToClient() for format of data
}

D3TimeChart.prototype.setXAxisTitle = function( title ){
  this.options.xtitle = title;
  //redraw x-title
}

D3TimeChart.prototype.setY1AxisTitle = function( title ){
  this.options.y1title = title;
  //redraw y1-title (e.g., gamma CPS axis title)
}

D3TimeChart.prototype.setY2AxisTitle = function(){
  this.options.y2title = title;
  //redraw y2-title (e.g., neutron CPS axis title)
}

D3TimeChart.prototype.setCompactXAxis = function( compact ){
  this.options.compactXAxis = compact;
  //Make x-zis title comapact or not
}

D3TimeChart.prototype.setGridX = function( show ){
  this.options.gridx = show;
  //add/remove horizantal grid lines
}

D3TimeChart.prototype.setGridY = function( show ){
  this.options.gridy = show;
  //add/remove vertical grid lines
}

D3TimeChart.prototype.setXAxisZoomSamples = function( firstsample, lastsample ){
  //Set it so only firstsample through lastsample sampels are visible
  //...
  //...
  
  //An example of using WtEmit:
  this.WtEmit(this.chart.id, {name: 'timerangechange'}, firstsample, lastsample );
}


