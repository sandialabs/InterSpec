RelEffPlot = function (elem) {
    let self = this;
    this.chart = typeof elem === "string" ? document.getElementById(elem) : elem;

    // Set the dimensions of the canvas / graph
    this.margins = { top: 30, right: 20, bottom: 30, left: 50 };
    const width = this.chart.clientWidth - this.margins.left - this.margins.right;
    const height = this.chart.clientHeight - this.margins.top - this.margins.bottom;

    // Setup the tooltip
    this.tooltip = d3.select("body").append("div")
        .attr("class", "RelEffPlotTooltip")
        .style("opacity", 0);

    // Set the ranges
    this.xScale = d3.scale.linear().range([0, width]);
    this.yScale = d3.scale.linear().range([height, 0]);

    // Define the axes
    this.xAxis = d3.svg
        .axis()
        .scale(this.xScale)
        .orient("bottom")
        .ticks(5);

    this.yAxis = d3.svg
        .axis()
        .scale(this.yScale)
        .orient("left")
        .ticks(5);

    // Adds the svg canvas
    this.svg = d3.select(this.chart)
        .append("svg")
        .attr("width", this.chart.clientWidth)
        .attr("height", this.chart.clientHeight)
        .attr("class", "RelEffPlot")
        .append("g")
        .attr("transform", "translate(" + this.margins.left + "," + this.margins.top + ")");

    this.xScale.domain([0, 3000]);
    this.yScale.domain([0, 1]);

    // Add the valueline path.
    this.path = this.svg.append("path")		// Add the valueline path.
        .attr("class", "line");

    // Add the X Axis
    this.svg.append("g")
        .attr("class", "xAxis")
        .attr("transform", "translate(0," + height + ")")
        .call(this.xAxis);

    // Add the Y Axis
    this.svg.append("g")
        .attr("class", "yAxis")
        .call(this.yAxis);
}//RelEffPlot constructor


RelEffPlot.prototype.setRelEffData = function (data_vals, fit_eqn) {
    const self = this;

    const width = this.chart.clientWidth - this.margins.left - this.margins.right;
    const height = this.chart.clientHeight - this.margins.top - this.margins.bottom;

    this.xScale
        .range([0, width]);
    this.yScale
        .range([height, 0]);
    this.svg
        .selectAll('.xAxis')
        .attr("transform", "translate(0," + height + ")");

    d3.select(this.chart).select("svg")
        .attr("width", this.chart.clientWidth)
        .attr("height", this.chart.clientHeight)

    this.data_vals = data_vals;
    this.fit_eqn = fit_eqn;

    if (!this.data_vals || !this.fit_eqn) {
        this.xScale.domain([0, 3000]);
        this.yScale.domain([0, 1]);
        this.svg.selectAll('.xAxis').call(this.xAxis);
        this.svg.selectAll('.yAxis').call(this.yAxis);
        this.path.attr("d", null);
        this.svg.selectAll("circle").remove();
        this.svg.selectAll("line.errorbar").remove();

        return;
    }//if( we dont have both data the rel eff function )


    
    
    const num_eqn_points = width / 4;

    // Scale the range of the data
    let min_x = d3.min(data_vals, function (d) { return d.energy; });
    let max_x = d3.max(data_vals, function (d) { return d.energy; });
    const initial_x_range = max_x - min_x;
    min_x -= ((min_x >=0 && min_x < (0.05*initial_x_range)) ? min_x : (0.05*initial_x_range));
    max_x += 0.05 * initial_x_range;


    let fit_eqn_points = [];
    for (let i = 0; i < num_eqn_points; ++i) {
        const ene = min_x + i * (max_x - min_x) / num_eqn_points;
        fit_eqn_points.push({ energy: ene, eff: fit_eqn(ene) });
    }

    let min_y = d3.min(data_vals, function (d) { return d.eff; });
    let max_y = d3.max(data_vals, function (d) { return d.eff; });
    min_y = Math.min(min_y, d3.min(fit_eqn_points, function (d) { return d.eff; }));
    max_y = Math.max(max_y, d3.max(fit_eqn_points, function (d) { return d.eff; }));

    const initial_y_range = max_y - min_y;
    min_y -= 0.1 * initial_y_range;
    max_y += 0.2 * initial_y_range;

    this.xScale.domain([min_x, max_x]);
    this.yScale.domain([min_y, max_y]);

    // Update the X Axis
    this.svg.selectAll('.xAxis')
        .call(this.xAxis);

    // Update the Y Axis
    this.svg.selectAll('.yAxis')
        .call(this.yAxis);

    // Define the line
    let valueline = d3.svg.line()
        .x(function (d) { return self.xScale(d.energy); })
        .y(function (d) { return self.yScale(d.eff); });

    this.path
        .attr("d", valueline(fit_eqn_points));


    // For some reason the circles position arent being updated on resize, so we'll just remove them first
    //  With later versions of D3 there is a .merge() function that is maybe relevant
    this.svg.selectAll("circle").remove();
    this.svg.selectAll("line.errorbar").remove();
    
    let lines = this.svg.selectAll("line.errorbar")
                .data(data_vals);
    
    lines.enter()
      .append('line')
      .attr("class", "errorbar")
      .attr('x1', function(d) { return self.xScale(d.energy); })
      .attr('x2', function(d) { return self.xScale(d.energy); })
      .attr('y1', function(d) { return self.yScale(d.eff + d.eff_uncert); })
      .attr('y2', function(d) { return self.yScale(d.eff - d.eff_uncert); });
    
    // Add the data points
    this.svg
        .selectAll("circle")
        .data(data_vals)
        .enter()
        .append("circle")
        .attr("r", 5)
        .attr("cx", function (d) {
            return self.xScale(d.energy)
        })
        .attr("cy", function (d) {
            return self.yScale(d.eff)
        })
        .attr("class", function (d) {
          
          if( d.nuc_info.length === 0 )
            return "noiso";

          if( d.nuc_info.length !== 1 )
            return "multiiso";
          
          const nuc = d.nuc_info[0].nuc;
          return nuc;
        })
        .on("mouseover", function (d, i) {

            self.tooltip.transition()
                .duration(200)
                .style("opacity", .9);

            d3.select(this).transition()
                .duration('50')
                .attr('opacity', '.85')
                .attr("r", 8);


            let txt = "<div>Energy: " + d.energy.toFixed(2) + " keV</div>"
                + "<div>Peak Area: " + d.counts.toFixed(1) + " &pm; " + d.counts_uncert.toFixed(1) + "</div>"
                + "<div>Measured RelEff: " + d.eff.toPrecision(5) + "</div>"
                + "<div>RelEff Curve: " + fit_eqn(d.energy).toPrecision(5) + "</div>";
            for (const el of d.nuc_info) {
                txt += "<div>&nbsp;&nbsp;" + el.nuc + ": br=" + el.br.toPrecision(4) + ", RelAct=" + el.rel_act.toPrecision(4) + "</div>"
            }
            
            self.tooltip.html(txt);

            // Make it so tooltip doesnt extend above/below/left/right of chart area
            const svg_location = d3.mouse( self.svg.node() );
            const svg_bb = self.svg.node().getBoundingClientRect();
            const tt_bb = self.tooltip.node().getBoundingClientRect();
            const render_right = ((svg_location[0] + tt_bb.width + 10 + 15) < svg_bb.width);
            const render_top = (svg_location[1] + tt_bb.height + 4 + 15 > svg_bb.height)
            const x_offset = render_right ? d3.event.pageX + 10 : d3.event.pageX - 10 - tt_bb.width;
            const y_offset = render_top ? d3.event.pageY - 10 - tt_bb.height : d3.event.pageY + 4;
            
            self.tooltip
              .style("left", x_offset + "px")
              .style("top", y_offset + "px");
        })
        .on('mouseout', function (d, i) {
            d3.select(this).transition()
                .duration('50')
                .attr('opacity', '1')
                .attr("r", 5);
            self.tooltip.transition()
                .duration(500)
                .style("opacity", 0);
        })
        ;
};//RelEffPlot.prototype.setRelEffData


RelEffPlot.prototype.handleResize = function () {
    this.setRelEffData(this.data_vals, this.fit_eqn);
};//RelEffPlot.prototype.handleResize