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



MdaChi2Chart = function(elem, options) {
  let self = this;
  
  this.parent = typeof elem === 'string' ? document.getElementById(elem) : elem;
  
  this.margins = {top: 10, right: 5, bottom: 5, left: 5};
  
  const chartAreaWidth = Math.max(0,this.parent.clientWidth - this.margins.left - this.margins.right);
  const chartAreaHeight = Math.max(0,this.parent.clientHeight - this.margins.top - this.margins.bottom);
  
  // Set some simple example data for later when I forget format
  this.data = {
    data: [{x:0, y:0}, {x:1, y:1}],
    xtitle: "Activity (&mu;Ci)",
    ytitle: "Chi2",
    lineColor: "black",
    //axisColor: "black",
    chartBackgroundColor: "rgba(0,0,0,0)",
    //textColor: "black"
  };
  
  
  this.xscale = d3.scale.linear();
  this.yscale = d3.scale.linear();
  
  this.xaxis = d3.svg.axis()
    .scale(this.xscale)
    .orient("bottom").ticks(5);
  
  this.yaxis = d3.svg.axis()
    .scale(this.yscale)
    .orient("left").ticks(5);
  
  this.svg = d3.select(this.parent)
    .append("svg")
    .attr("width",this.parent.clientWidth)
    .attr("height",this.parent.clientHeight);
  
  this.rect = this.svg.append('rect') // outline for reference
    .attr({x: this.margins.left, y: this.margins.top,
       width: chartAreaWidth,
       height: chartAreaHeight,
       stroke: this.data.axisColor,
       'stroke-width': 0.5,
    fill: this.data.chartBackgroundColor }); // attributes in JS list format
  
  
  this.path = this.svg.append("path")
    .attr("class", "line")
    .style("fill", "none")
    .style( "stroke", this.data.lineColor )
    .style("stroke-width", 2);
  
  this.xaxisg = this.svg.append("g")         // Add the X Axis
    .attr("class", "xaxis")
    .call(this.xaxis)
    .attr( "style", "fill: none; stroke-width: 1;");

  this.yaxisg = this.svg.append("g")         // Add the Y Axis
    .attr("class", "yaxis")
    .call(this.yaxis)
    .attr( "style", "fill: none; stroke-width: 1;");
    
  // Add the text label for the x axis
  this.xtitle = this.svg.append("text")
    .attr("class", "xaxistitle")
    .style("text-anchor", "middle")
    .attr("dy", "-0.5em")
    .text("&nbsp;");
  
  // Add the text label for the Y axis
  this.ytitle = this.svg.append("text")
    .attr("transform", "rotate(-90)")
    .attr("class", "yaxistitle")
    .attr("y", 0)
    .attr("x",0)
    .attr("dy", "1em")
    .style("text-anchor", "middle")
    .text("&nbsp;");
  
  this.redraw();
}


MdaChi2Chart.prototype.redraw = function() {
  let self = this;
  
  // We'll make two passes at setting things.  Once to set approx axises and stuff, then get those
  //  dimensions, and use them to set final positions
  let leftMargin = this.margins.left;
  let bottomMargin = this.margins.bottom;
  let chartAreaWidth = Math.max(0,this.parent.clientWidth - leftMargin - this.margins.right);
  let chartAreaHeight = Math.max(0,this.parent.clientHeight - this.margins.top - bottomMargin);
  
  this.svg
    .attr("width",this.parent.clientWidth )
    .attr("height",this.parent.clientHeight );

  //console.log( "this.data.data: ", this.data.data );
  
  this.xscale
    .domain(d3.extent(this.data.data, function(d) {return d.x;}))
    .range([leftMargin, chartAreaWidth + leftMargin]);
  
  this.yscale
    .domain(d3.extent(this.data.data, function(d) {return d.y;}))
    .range([chartAreaHeight + this.margins.top, this.margins.top]);
  
  this.xaxis.scale( this.xscale );
  this.yaxis.scale( this.yscale );
    
  this.xaxisg.call(this.xaxis);
  this.yaxisg.call(this.yaxis);
  
  this.xtitle.text(this.data.xtitle);
  this.ytitle.text(this.data.ytitle);
    
    
  // Now get the actual axis and label sizes, and reposition everything correctly
  const xTitleHeight = this.xtitle.node().getBBox().height;
  const xAxisHeight = this.xaxisg.node().getBBox().height;
  
  const yTitleWidth = this.ytitle.node().getBBox().width;
  const yAxisWidth = this.yaxisg.node().getBBox().width;
  
  //console.log( 'xTitleHeight:', xTitleHeight, ', xAxisHeight:', xAxisHeight );
  //console.log( 'yTitleWidth:', yTitleWidth, ', yAxisWidth:', yAxisWidth );
  
  leftMargin = this.margins.left + yTitleWidth + yAxisWidth + 2; //2 px between labels and title
  bottomMargin = this.margins.bottom + xTitleHeight + xAxisHeight;
  
  chartAreaWidth = Math.max(0,this.parent.clientWidth - leftMargin - this.margins.right);
  chartAreaHeight = Math.max(0,this.parent.clientHeight - this.margins.top - bottomMargin);
  
  this.xscale
    .range([leftMargin, chartAreaWidth + leftMargin]);
  this.yscale
    .range([chartAreaHeight + this.margins.top, this.margins.top]);
  
  this.xaxis.scale( this.xscale );
  this.yaxis.scale( this.yscale );
  
  this.xtitle
    .attr("transform", "translate(" + (leftMargin + 0.5*chartAreaWidth) + "," + (chartAreaHeight + bottomMargin + this.margins.top) + ")");
  
  this.ytitle
    .attr("x",0 - ((chartAreaHeight + this.margins.top) / 2));
    
  this.xaxisg
    .attr("transform", "translate(0," + (chartAreaHeight + this.margins.top) + ")")
    .call(this.xaxis);
  
  this.yaxisg
      .attr("transform", "translate(" + leftMargin + ",0)")
      .call(this.yaxis);
      
  this.rect.attr({
        x: leftMargin,
        y: this.margins.top,
        width: chartAreaWidth,
        height: chartAreaHeight,
        fill: this.data.chartBackgroundColor
  });
      
  let valueline = d3.svg.line()
    .x(function(d) { return self.xscale(d.x);}) // apply the x scale to the x data
    .y(function(d) { return self.yscale(d.y);}) // apply the y scale to the y data
      
  this.path
    .attr("d", valueline(this.data.data))
    .style( "stroke", this.data.lineColor );
}


MdaChi2Chart.prototype.setData = function( data ){
  this.data = data;
  this.redraw();
}
