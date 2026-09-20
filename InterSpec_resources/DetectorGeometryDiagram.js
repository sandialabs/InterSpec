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

/* A side-elevation cross-section of a detector geometry (see DetectorGeometryDiagram.cpp,
   which does the geometry and hands over a list of r-z polycone regions).  d3 v3; ES2019-safe.

   Data: { box, zMin, zMax, rMax, knots: { z: [...], r: [...] }, labels: { front: "..." },
           regions: [ { id, kind, profile: [[z, rmin, rmax], ...], tip: { title, lines: [...] } } ] }
   Regions are listed outermost first; the crystal front face is z = 0, the source side negative z.

   The drawing scales with the container.  A common scale is chosen such that everything fits, but
   any gap between two knots thinner than `minPx` on screen is widened to `minPx` (the map from
   cm to pixels is piecewise linear through the knots), so a 0.5 mm window on a 10 cm detector is
   still a visible line.  Orientation follows the container: front at the top when it is taller
   than wide, front at the left otherwise.
 */
DetectorGeometryDiagram = function( elem, options )
{
  this.container = (typeof elem === 'string') ? document.getElementById(elem) : elem;

  // A stub element (a widget on a tab that has never been shown) resolves to null.  Returning here
  //  is what keeps a TypeError from aborting the rest of Wt's update block and taking unrelated
  //  widgets' event wiring with it.  The C++ side re-creates this object, with the geometry, once
  //  the real element replaces the stub.
  if( !this.container )
    return;

  // Wt re-emits every JavaScript member when it recreates a widget's DOM, so this constructor can
  //  run a second time on an element that already has our SVG and tooltip in it.  Clear them (and
  //  drop any previous observer) rather than stacking a second copy on top.
  if( this.roomObserver )
    this.roomObserver.disconnect();
  d3.select(this.container).selectAll('svg,div.DgdTooltip').remove();
  this.options = options || {};
  this.minPx = this.options.minPx || 2;
  // Below this much room beside the form, the drawing hides itself rather than being squeezed.
  this.minRoomPx = this.options.minRoomPx || 200;
  // The drawing takes its height from the form; below this it is a sliver, so hide it instead.
  this.minRoomHeightPx = this.options.minRoomHeightPx || 120;
  this.data = null;

  this.svg = d3.select(this.container).append('svg').attr('class', 'DgdSvg');
  this.g = this.svg.append('g');
  this.tooltip = d3.select(this.container).append('div').attr('class', 'DgdTooltip');

  // The room available is a property of the PARENT, not of this element: once hidden, this element
  // has no width of its own to measure and its own ResizeObserver stops firing.  So watch the row
  // the drawing shares with the form.
  const self = this;
  this.roomObserver = null;
  const parent = this.container.parentElement;
  if( parent && window.ResizeObserver )
  {
    this.roomObserver = new ResizeObserver( function(){ self.handleResize(); } );
    this.roomObserver.observe( parent );
  }
  this.updateRoom();
};

/* Width left for the drawing in the row it shares with the form: the row's width, less every other
   item in it and the gaps between them.  Valid whether or not the drawing is currently showing. */
DetectorGeometryDiagram.prototype.availableWidth = function()
{
  const parent = this.container.parentElement;
  if( !parent )
    return this.container.clientWidth;

  const style = window.getComputedStyle( parent );
  const gap = parseFloat( style.columnGap ) || parseFloat( style.gap ) || 0;

  let used = 0, gaps = 0;
  for( let i = 0; i < parent.children.length; ++i )
  {
    const sib = parent.children[i];
    if( sib === this.container )
      continue;
    if( window.getComputedStyle(sib).display === 'none' )
      continue;
    used += sib.offsetWidth;
    gaps += gap;
  }

  return parent.clientWidth - used - gaps;
};

/* Height available to the drawing: the row's, since this element stretches to it. */
DetectorGeometryDiagram.prototype.availableHeight = function()
{
  const parent = this.container.parentElement;
  return parent ? parent.clientHeight : this.container.clientHeight;
};

/* Shows or hides the drawing per #availableWidth / #availableHeight; returns whether it is showing. */
DetectorGeometryDiagram.prototype.updateRoom = function()
{
  const hide = ( this.availableWidth() < this.minRoomPx )
               || ( this.availableHeight() < this.minRoomHeightPx );
  const wasHidden = this.container.classList.contains( 'DgdNoRoom' );
  if( hide !== wasHidden )
  {
    if( hide )
      this.container.classList.add( 'DgdNoRoom' );
    else
      this.container.classList.remove( 'DgdNoRoom' );
    this.hideTooltip();
  }

  return !hide;
};

DetectorGeometryDiagram.prototype.setData = function( data )
{
  this.data = data;
  this.updateRoom();
  this.render();
};

DetectorGeometryDiagram.prototype.clear = function()
{
  this.data = null;
  this.g.selectAll('*').remove();
  this.hideTooltip();
};

DetectorGeometryDiagram.prototype.handleResize = function()
{
  // Re-decide whether there is room before drawing; `render()` no-ops while hidden (zero width).
  this.updateRoom();
  this.render();
};

DetectorGeometryDiagram.prototype.hideTooltip = function()
{
  this.tooltip.style('opacity', 0);
};

/* Sum of the on-screen lengths of the gaps between consecutive knots at scale s (px/cm), each gap
   at least minPx wide. */
function dgdInflatedExtent( knots, s, minPx )
{
  let e = 0;
  for( let i = 0; i + 1 < knots.length; ++i )
    e += Math.max( s*(knots[i+1] - knots[i]), minPx );
  return e;
}

/* A monotone piecewise-linear map cm -> px through the knots at scale s, thin gaps widened to
   minPx; values between knots interpolate within their gap, values outside continue at scale s. */
function dgdBuildMap( knots, s, minPx )
{
  const n = knots.length;
  const pos = [0];
  for( let i = 0; i + 1 < n; ++i )
    pos.push( pos[i] + Math.max( s*(knots[i+1] - knots[i]), minPx ) );

  return function( v ){
    if( n === 0 )
      return s*v;
    if( v <= knots[0] )
      return pos[0] + s*(v - knots[0]);
    if( v >= knots[n-1] )
      return pos[n-1] + s*(v - knots[n-1]);
    let lo = 0, hi = n - 1;
    while( hi - lo > 1 )
    {
      const mid = (lo + hi) >> 1;
      if( knots[mid] <= v )
        lo = mid;
      else
        hi = mid;
    }
    const span = knots[lo+1] - knots[lo];
    const f = (span > 0) ? (v - knots[lo]) / span : 0;
    return pos[lo] + f*(pos[lo+1] - pos[lo]);
  };
}

DetectorGeometryDiagram.prototype.render = function()
{
  const self = this;
  const w = this.container.clientWidth;
  const h = this.container.clientHeight;

  this.g.selectAll('*').remove();
  this.hideTooltip();

  const d = this.data;
  if( !d || !d.regions || !d.knots || (w < 60) || (h < 60) )
    return;

  this.svg.attr('width', w).attr('height', h);

  const minPx = this.minPx;
  const margin = 8;
  const labelBand = 14;   //room for the "front" label on the source side
  const vertical = (h >= 0.85*w);
  const availZ = (vertical ? h : w) - 2*margin - labelBand;
  const availR = (vertical ? w : h) - 2*margin;

  const zK = d.knots.z.slice().sort( (a,b) => a - b );
  const rK = d.knots.r.slice().sort( (a,b) => a - b );
  if( (zK.length < 2) || (rK.length < 2) || (availZ < 10) || (availR < 10) )
    return;

  // The largest common scale at which both the length and the (mirrored) width fit.
  const zSpan = zK[zK.length-1] - zK[0];
  const rSpan = rK[rK.length-1] - rK[0];
  let lo = 0, hi = 2*Math.max( availZ/Math.max(zSpan,1e-9), availR/Math.max(2*rSpan,1e-9) );
  for( let it = 0; it < 40; ++it )
  {
    const mid = 0.5*(lo + hi);
    if( (dgdInflatedExtent(zK, mid, minPx) <= availZ) && (2*dgdInflatedExtent(rK, mid, minPx) <= availR) )
      lo = mid;
    else
      hi = mid;
  }
  const s = Math.max( lo, 1e-6 );
  const mapZ = dgdBuildMap( zK, s, minPx );
  const mapR = dgdBuildMap( rK, s, minPx );
  const zExt = dgdInflatedExtent( zK, s, minPx );
  const rExt = dgdInflatedExtent( rK, s, minPx );

  // Screen placement: the drawing centred, the label band on the source side.
  let toXY, axisFrom, axisTo, labelPos, labelRotate;
  if( vertical )
  {
    const originY = 0.5*(h - (zExt + labelBand)) + labelBand;
    const centerX = 0.5*w;
    toXY = function( r, z ){
      const px = (r >= 0) ? mapR(r) : -mapR(-r);
      return [ centerX + px, originY + mapZ(z) ];
    };
    axisFrom = [centerX, originY - 4];
    axisTo = [centerX, originY + zExt + 4];
    labelPos = [centerX, originY - 4];
    labelRotate = 0;
  }else
  {
    const originX = 0.5*(w - (zExt + labelBand)) + labelBand;
    const centerY = 0.5*h;
    toXY = function( r, z ){
      const py = (r >= 0) ? mapR(r) : -mapR(-r);
      return [ originX + mapZ(z), centerY - py ];
    };
    axisFrom = [originX - 4, centerY];
    axisTo = [originX + zExt + 4, centerY];
    labelPos = [originX - 4, centerY];
    labelRotate = -90;
  }

  const fmt = function( p ){ return p[0].toFixed(2) + ',' + p[1].toFixed(2); };

  // Fills: one closed path per region, the full mirrored section (right outer, right inner, left
  //  inner, left outer).  A zero-width pinch along the axis cancels under nonzero fill.
  const regions = this.g.selectAll('g.DgdRegion').data( d.regions ).enter()
                    .append('g')
                    .attr('class', function(r){ return 'DgdRegion Dgd-' + r.kind; })
                    .attr('data-id', function(r){ return r.id; });

  regions.append('path').attr('class', 'DgdFill').attr('d', function( r ){
    const prof = r.profile;
    if( !prof || prof.length < 2 )
      return '';
    const pts = [];
    for( let i = 0; i < prof.length; ++i )
      pts.push( toXY( +prof[i][2], prof[i][0] ) );
    for( let i = prof.length - 1; i >= 0; --i )
      pts.push( toXY( +prof[i][1], prof[i][0] ) );
    for( let i = 0; i < prof.length; ++i )
      pts.push( toXY( -prof[i][1], prof[i][0] ) );
    for( let i = prof.length - 1; i >= 0; --i )
      pts.push( toXY( -prof[i][2], prof[i][0] ) );
    return 'M' + pts.map(fmt).join('L') + 'Z';
  });

  // Outlines, drawn after every fill so no region paints over a neighbour's edge: the outer edge
  //  on both sides, the front and back caps, and the inner (bore-side) edge only where it is off
  //  the axis.
  const outlines = this.g.append('g').attr('class', 'DgdOutline');
  const addLine = function( pts ){
    if( pts.length >= 2 )
      outlines.append('path').attr('d', 'M' + pts.map(fmt).join('L'));
  };
  d.regions.forEach( function( r ){
    const prof = r.profile;
    if( !prof || prof.length < 2 )
      return;
    for( const sign of [1, -1] )
    {
      addLine( prof.map( function(p){ return toXY( sign*p[2], p[0] ); } ) );
      // inner edge: runs of consecutive planes with rmin > 0
      let run = [];
      for( let i = 0; i < prof.length; ++i )
      {
        if( prof[i][1] > 0 )
          run.push( toXY( sign*prof[i][1], prof[i][0] ) );
        else
        {
          addLine( run );
          run = [];
        }
      }
      addLine( run );
    }
    // caps
    for( const idx of [0, prof.length - 1] )
    {
      const p = prof[idx];
      if( p[1] > 0 )
      {
        addLine( [ toXY(p[1], p[0]), toXY(p[2], p[0]) ] );
        addLine( [ toXY(-p[1], p[0]), toXY(-p[2], p[0]) ] );
      }else
      {
        addLine( [ toXY(-p[2], p[0]), toXY(p[2], p[0]) ] );
      }
    }
  });

  // Axis and the "front" label.
  this.g.append('path').attr('class', 'DgdAxis')
        .attr('d', 'M' + fmt(axisFrom) + 'L' + fmt(axisTo));
  const label = this.g.append('text').attr('class', 'DgdLabel')
        .attr('text-anchor', 'middle')
        .text( (d.labels && d.labels.front) ? d.labels.front : '' );
  if( labelRotate )
    label.attr('transform', 'translate(' + labelPos[0] + ',' + labelPos[1] + ') rotate(' + labelRotate + ')');
  else
    label.attr('x', labelPos[0]).attr('y', labelPos[1]);

  // Tooltips: title plus one line per fact, as text nodes (the material name is user text).
  const showTip = function( r ){
    const node = self.tooltip.node();
    while( node.firstChild )
      node.removeChild( node.firstChild );
    const b = document.createElement('b');
    b.appendChild( document.createTextNode( (r.tip && r.tip.title) ? r.tip.title : r.id ) );
    node.appendChild( b );
    const lines = (r.tip && r.tip.lines) ? r.tip.lines : [];
    for( let i = 0; i < lines.length; ++i )
    {
      node.appendChild( document.createElement('br') );
      node.appendChild( document.createTextNode( lines[i] ) );
    }
    self.tooltip.style('opacity', 1);
  };
  const moveTip = function(){
    const m = d3.mouse( self.container );
    const node = self.tooltip.node();
    const tw = node.offsetWidth || 200, th = node.offsetHeight || 80, off = 10;
    let left = m[0] + off, top = m[1] - off;
    if( left + tw > w ) left = m[0] - tw - off;
    if( left < 0 ) left = off;
    if( top + th > h ) top = m[1] - th - off;
    if( top < 0 ) top = off;
    self.tooltip.style('left', left + 'px').style('top', top + 'px');
  };
  regions.on('mouseover', function( r ){
      d3.select(this).classed('DgdHover', true);
      showTip( r );
      moveTip();
    })
    .on('mousemove', moveTip)
    .on('mouseout', function(){
      d3.select(this).classed('DgdHover', false);
      self.hideTooltip();
    })
    .on('touchstart', function( r ){
      showTip( r );
      moveTip();
    })
    .on('touchend', function(){ self.hideTooltip(); })
    .on('touchcancel', function(){ self.hideTooltip(); });
};
