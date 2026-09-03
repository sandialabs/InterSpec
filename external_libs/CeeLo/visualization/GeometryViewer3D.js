/* CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
 and X-ray detector efficiency - developed as part of InterSpec.

 Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
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
 * GeometryViewer3D - WebGL 3D geometry viewer for CeeLo GDML files.
 * No external dependencies. Pure WebGL + vanilla JS.
 *
 * Usage:
 *   var viewer = new GeometryViewer3D(document.getElementById('viewer'), {});
 *   viewer.loadGDML(gdmlXmlString);
 *   viewer.setSourcePosition(0, 0, -10);
 *
 * Compatible with future C++ Wt wrapper (follows RelEffPlot.js pattern).
 */

/* jshint esversion: 6 */
/* global DOMParser, XMLSerializer */

// ============================================================================
// Vector / Matrix math utilities (no external deps)
// ============================================================================

var GeometryViewer3D;

(function() {
"use strict";

function vec3_create(x, y, z) { return [x, y, z]; }
function vec3_sub(a, b) { return [a[0]-b[0], a[1]-b[1], a[2]-b[2]]; }
function vec3_add(a, b) { return [a[0]+b[0], a[1]+b[1], a[2]+b[2]]; }
function vec3_scale(v, s) { return [v[0]*s, v[1]*s, v[2]*s]; }
function vec3_dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function vec3_cross(a, b) {
  return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}
function vec3_length(v) { return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
function vec3_normalize(v) {
  var l = vec3_length(v);
  return l > 1e-12 ? [v[0]/l, v[1]/l, v[2]/l] : [0,0,0];
}
function vec3_lerp(a, b, t) {
  return [a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t, a[2]+(b[2]-a[2])*t];
}

// Column-major 4x4 matrix (OpenGL convention)
function mat4_identity() {
  return [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1];
}

function mat4_multiply(a, b) {
  var r = new Array(16);
  for (var c = 0; c < 4; c++) {
    for (var row = 0; row < 4; row++) {
      r[c*4+row] = a[row]*b[c*4] + a[4+row]*b[c*4+1] + a[8+row]*b[c*4+2] + a[12+row]*b[c*4+3];
    }
  }
  return r;
}

function mat4_perspective(fovY, aspect, near, far) {
  var f = 1.0 / Math.tan(fovY / 2.0);
  var nf = 1.0 / (near - far);
  return [f/aspect,0,0,0, 0,f,0,0, 0,0,(far+near)*nf,-1, 0,0,2*far*near*nf,0];
}

function mat4_lookAt(eye, center, up) {
  var f = vec3_normalize(vec3_sub(center, eye));
  var s = vec3_normalize(vec3_cross(f, up));
  var u = vec3_cross(s, f);
  return [
    s[0], u[0], -f[0], 0,
    s[1], u[1], -f[1], 0,
    s[2], u[2], -f[2], 0,
    -vec3_dot(s,eye), -vec3_dot(u,eye), vec3_dot(f,eye), 1
  ];
}

function mat4_translate(tx, ty, tz) {
  return [1,0,0,0, 0,1,0,0, 0,0,1,0, tx,ty,tz,1];
}

function mat4_invert(m) {
  var inv = new Array(16);
  inv[0]  =  m[5]*m[10]*m[15] - m[5]*m[11]*m[14] - m[9]*m[6]*m[15] + m[9]*m[7]*m[14] + m[13]*m[6]*m[11] - m[13]*m[7]*m[10];
  inv[4]  = -m[4]*m[10]*m[15] + m[4]*m[11]*m[14] + m[8]*m[6]*m[15] - m[8]*m[7]*m[14] - m[12]*m[6]*m[11] + m[12]*m[7]*m[10];
  inv[8]  =  m[4]*m[9]*m[15]  - m[4]*m[11]*m[13] - m[8]*m[5]*m[15] + m[8]*m[7]*m[13] + m[12]*m[5]*m[11] - m[12]*m[7]*m[9];
  inv[12] = -m[4]*m[9]*m[14]  + m[4]*m[10]*m[13] + m[8]*m[5]*m[14] - m[8]*m[6]*m[13] - m[12]*m[5]*m[10] + m[12]*m[6]*m[9];
  inv[1]  = -m[1]*m[10]*m[15] + m[1]*m[11]*m[14] + m[9]*m[2]*m[15] - m[9]*m[3]*m[14] - m[13]*m[2]*m[11] + m[13]*m[3]*m[10];
  inv[5]  =  m[0]*m[10]*m[15] - m[0]*m[11]*m[14] - m[8]*m[2]*m[15] + m[8]*m[3]*m[14] + m[12]*m[2]*m[11] - m[12]*m[3]*m[10];
  inv[9]  = -m[0]*m[9]*m[15]  + m[0]*m[11]*m[13] + m[8]*m[1]*m[15] - m[8]*m[3]*m[13] - m[12]*m[1]*m[11] + m[12]*m[3]*m[9];
  inv[13] =  m[0]*m[9]*m[14]  - m[0]*m[10]*m[13] - m[8]*m[1]*m[14] + m[8]*m[2]*m[13] + m[12]*m[1]*m[10] - m[12]*m[2]*m[9];
  inv[2]  =  m[1]*m[6]*m[15]  - m[1]*m[7]*m[14]  - m[5]*m[2]*m[15] + m[5]*m[3]*m[14] + m[13]*m[2]*m[7]  - m[13]*m[3]*m[6];
  inv[6]  = -m[0]*m[6]*m[15]  + m[0]*m[7]*m[14]  + m[4]*m[2]*m[15] - m[4]*m[3]*m[14] - m[12]*m[2]*m[7]  + m[12]*m[3]*m[6];
  inv[10] =  m[0]*m[5]*m[15]  - m[0]*m[7]*m[13]  - m[4]*m[1]*m[15] + m[4]*m[3]*m[13] + m[12]*m[1]*m[7]  - m[12]*m[3]*m[5];
  inv[14] = -m[0]*m[5]*m[14]  + m[0]*m[6]*m[13]  + m[4]*m[1]*m[14] - m[4]*m[2]*m[13] - m[12]*m[1]*m[6]  + m[12]*m[2]*m[5];
  inv[3]  = -m[1]*m[6]*m[11]  + m[1]*m[7]*m[10]  + m[5]*m[2]*m[11] - m[5]*m[3]*m[10] - m[9]*m[2]*m[7]   + m[9]*m[3]*m[6];
  inv[7]  =  m[0]*m[6]*m[11]  - m[0]*m[7]*m[10]  - m[4]*m[2]*m[11] + m[4]*m[3]*m[10] + m[8]*m[2]*m[7]   - m[8]*m[3]*m[6];
  inv[11] = -m[0]*m[5]*m[11]  + m[0]*m[7]*m[9]   + m[4]*m[1]*m[11] - m[4]*m[3]*m[9]  - m[8]*m[1]*m[7]   + m[8]*m[3]*m[5];
  inv[15] =  m[0]*m[5]*m[10]  - m[0]*m[6]*m[9]   - m[4]*m[1]*m[10] + m[4]*m[2]*m[9]  + m[8]*m[1]*m[6]   - m[8]*m[2]*m[5];

  var det = m[0]*inv[0] + m[1]*inv[4] + m[2]*inv[8] + m[3]*inv[12];
  if (Math.abs(det) < 1e-20) return mat4_identity();
  det = 1.0 / det;
  for (var i = 0; i < 16; i++) inv[i] *= det;
  return inv;
}

function mat3_normalFromMat4(m) {
  // Upper-left 3x3 inverse transpose
  var a00=m[0],a01=m[1],a02=m[2], a10=m[4],a11=m[5],a12=m[6], a20=m[8],a21=m[9],a22=m[10];
  var det = a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) + a02*(a10*a21-a11*a20);
  if (Math.abs(det) < 1e-20) return [1,0,0, 0,1,0, 0,0,1];
  det = 1.0/det;
  // Transpose of cofactor matrix
  return [
    (a11*a22-a12*a21)*det, (a12*a20-a10*a22)*det, (a10*a21-a11*a20)*det,
    (a02*a21-a01*a22)*det, (a00*a22-a02*a20)*det, (a01*a20-a00*a21)*det,
    (a01*a12-a02*a11)*det, (a02*a10-a00*a12)*det, (a00*a11-a01*a10)*det
  ];
}


// ============================================================================
// Mesh generation
// ============================================================================

function generateTubeMesh(rmin, rmax, halfZ, segments) {
  var positions = [], normals = [], indices = [];
  var i, a, ca, sa, idx;
  var TWO_PI = 2 * Math.PI;

  function pushVert(x,y,z, nx,ny,nz) {
    positions.push(x,y,z);
    normals.push(nx,ny,nz);
  }

  if (rmin < 1e-6) {
    // Solid cylinder
    // Outer side
    var baseIdx = 0;
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      ca = Math.cos(a); sa = Math.sin(a);
      pushVert(rmax*ca, rmax*sa, -halfZ, ca, sa, 0);
      pushVert(rmax*ca, rmax*sa, +halfZ, ca, sa, 0);
    }
    for (i = 0; i < segments; i++) {
      idx = baseIdx + i*2;
      indices.push(idx, idx+1, idx+2, idx+2, idx+1, idx+3);
    }

    // Top cap (z = +halfZ)
    baseIdx = positions.length / 3;
    pushVert(0, 0, +halfZ, 0, 0, 1); // center
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      pushVert(rmax*Math.cos(a), rmax*Math.sin(a), +halfZ, 0, 0, 1);
    }
    for (i = 0; i < segments; i++) {
      indices.push(baseIdx, baseIdx+1+i, baseIdx+2+i);
    }

    // Bottom cap (z = -halfZ)
    baseIdx = positions.length / 3;
    pushVert(0, 0, -halfZ, 0, 0, -1); // center
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      pushVert(rmax*Math.cos(a), rmax*Math.sin(a), -halfZ, 0, 0, -1);
    }
    for (i = 0; i < segments; i++) {
      indices.push(baseIdx, baseIdx+2+i, baseIdx+1+i);
    }
  } else {
    // Hollow tube (rmin > 0)
    // Outer side
    var base = 0;
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      ca = Math.cos(a); sa = Math.sin(a);
      pushVert(rmax*ca, rmax*sa, -halfZ, ca, sa, 0);
      pushVert(rmax*ca, rmax*sa, +halfZ, ca, sa, 0);
    }
    for (i = 0; i < segments; i++) {
      idx = base + i*2;
      indices.push(idx, idx+1, idx+2, idx+2, idx+1, idx+3);
    }

    // Inner side (normals inward)
    base = positions.length / 3;
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      ca = Math.cos(a); sa = Math.sin(a);
      pushVert(rmin*ca, rmin*sa, -halfZ, -ca, -sa, 0);
      pushVert(rmin*ca, rmin*sa, +halfZ, -ca, -sa, 0);
    }
    for (i = 0; i < segments; i++) {
      idx = base + i*2;
      indices.push(idx, idx+2, idx+1, idx+2, idx+3, idx+1);
    }

    // Top annular ring (z = +halfZ)
    base = positions.length / 3;
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      ca = Math.cos(a); sa = Math.sin(a);
      pushVert(rmax*ca, rmax*sa, +halfZ, 0, 0, 1);
      pushVert(rmin*ca, rmin*sa, +halfZ, 0, 0, 1);
    }
    for (i = 0; i < segments; i++) {
      idx = base + i*2;
      indices.push(idx, idx+2, idx+1, idx+1, idx+2, idx+3);
    }

    // Bottom annular ring (z = -halfZ)
    base = positions.length / 3;
    for (i = 0; i <= segments; i++) {
      a = (i / segments) * TWO_PI;
      ca = Math.cos(a); sa = Math.sin(a);
      pushVert(rmax*ca, rmax*sa, -halfZ, 0, 0, -1);
      pushVert(rmin*ca, rmin*sa, -halfZ, 0, 0, -1);
    }
    for (i = 0; i < segments; i++) {
      idx = base + i*2;
      indices.push(idx, idx+1, idx+2, idx+2, idx+1, idx+3);
    }
  }

  return {
    positions: new Float32Array(positions),
    normals: new Float32Array(normals),
    indices: new Uint32Array(indices)
  };
}

// Surface of revolution from a list of {z, rmin, rmax} planes -- a GDML
// <polycone>. CeeLo exports coaxial and/or bulletized crystals this way, so
// this is what draws the bore and the rounded front edge.
//
// Planes are given in the solid's own frame and may repeat a z value to make a
// step (that is how a bore mouth is expressed). Consecutive planes at the same
// z contribute no side wall, only the annular face between their radii.
function generatePolyconeMesh(planes, segments) {
  var positions = [], normals = [], indices = [];
  var TWO_PI = 2 * Math.PI;
  var i, k;

  // A polycone with no z-planes is malformed GDML; draw nothing rather than
  // throwing, which would blank the whole viewer.
  if (!planes || planes.length < 2) {
    return { positions: new Float32Array(0), normals: new Float32Array(0),
             indices: new Uint32Array(0) };
  }

  function pushVert(x, y, z, nx, ny, nz) {
    positions.push(x, y, z);
    normals.push(nx, ny, nz);
  }

  // Side walls (outer and inner), one quad strip per plane pair.
  function wall(rA, zA, rB, zB, outward) {
    if (rA < 1e-9 && rB < 1e-9) return;      // degenerate: on the axis
    if (Math.abs(zB - zA) < 1e-12 && Math.abs(rB - rA) < 1e-12) return;
    // Normal of the cone frustum in the meridian plane, rotated around z.
    var dz = zB - zA, dr = rB - rA;
    var len = Math.sqrt(dz * dz + dr * dr);
    if (len < 1e-12) return;
    var nr = (outward ? 1 : -1) * dz / len;
    var nz = (outward ? -1 : 1) * dr / len;
    var base = positions.length / 3;
    for (i = 0; i <= segments; i++) {
      var a = (i / segments) * TWO_PI, ca = Math.cos(a), sa = Math.sin(a);
      pushVert(rA * ca, rA * sa, zA, nr * ca, nr * sa, nz);
      pushVert(rB * ca, rB * sa, zB, nr * ca, nr * sa, nz);
    }
    for (i = 0; i < segments; i++) {
      var idx = base + i * 2;
      // Wind inner walls the other way so their front faces point into the
      // bore, matching generateTubeMesh's convention for a hollow tube.
      if (outward) indices.push(idx, idx + 1, idx + 2, idx + 2, idx + 1, idx + 3);
      else         indices.push(idx, idx + 2, idx + 1, idx + 2, idx + 3, idx + 1);
    }
  }

  // Annular face at constant z between two radii (end caps and bore steps).
  function annulus(rIn, rOut, z, nz) {
    if (rOut - rIn < 1e-9) return;
    var base = positions.length / 3;
    for (i = 0; i <= segments; i++) {
      var a = (i / segments) * TWO_PI, ca = Math.cos(a), sa = Math.sin(a);
      pushVert(rIn * ca, rIn * sa, z, 0, 0, nz);
      pushVert(rOut * ca, rOut * sa, z, 0, 0, nz);
    }
    for (i = 0; i < segments; i++) {
      var idx = base + i * 2;
      if (nz > 0) indices.push(idx, idx + 1, idx + 2, idx + 2, idx + 1, idx + 3);
      else        indices.push(idx, idx + 2, idx + 1, idx + 2, idx + 3, idx + 1);
    }
  }

  for (k = 0; k + 1 < planes.length; k++) {
    var p0 = planes[k], p1 = planes[k + 1];
    if (Math.abs(p1.z - p0.z) < 1e-12) {
      // A step: an annular face at constant z. The normal points away from the
      // material, which for an opening bore (rmin growing with z) is +z.
      if (Math.abs(p1.rmin - p0.rmin) > 1e-12) {
        annulus(Math.min(p0.rmin, p1.rmin), Math.max(p0.rmin, p1.rmin), p0.z,
                p1.rmin > p0.rmin ? +1 : -1);
      }
      if (Math.abs(p1.rmax - p0.rmax) > 1e-12) {
        annulus(Math.min(p0.rmax, p1.rmax), Math.max(p0.rmax, p1.rmax), p0.z,
                p1.rmax < p0.rmax ? +1 : -1);
      }
      continue;
    }
    wall(p0.rmax, p0.z, p1.rmax, p1.z, true);
    if (p0.rmin > 1e-9 || p1.rmin > 1e-9) {
      wall(p0.rmin, p0.z, p1.rmin, p1.z, false);
    }
  }

  // End caps.
  var first = planes[0], last = planes[planes.length - 1];
  annulus(first.rmin, first.rmax, first.z, -1);
  annulus(last.rmin, last.rmax, last.z, +1);

  return {
    positions: new Float32Array(positions),
    normals: new Float32Array(normals),
    indices: new Uint32Array(indices)
  };
}

function generateBoxMesh(halfX, halfY, halfZ) {
  // 6 faces, 4 verts each, 2 triangles each
  var p = [], n = [], idx = [];
  function face(v0,v1,v2,v3, nx,ny,nz) {
    var b = p.length / 3;
    p.push(v0[0],v0[1],v0[2]); n.push(nx,ny,nz);
    p.push(v1[0],v1[1],v1[2]); n.push(nx,ny,nz);
    p.push(v2[0],v2[1],v2[2]); n.push(nx,ny,nz);
    p.push(v3[0],v3[1],v3[2]); n.push(nx,ny,nz);
    idx.push(b,b+1,b+2, b+2,b+3,b);
  }
  var x=halfX, y=halfY, z=halfZ;
  face([+x,-y,-z],[+x,+y,-z],[+x,+y,+z],[+x,-y,+z], 1,0,0);  // +X
  face([-x,+y,-z],[-x,-y,-z],[-x,-y,+z],[-x,+y,+z], -1,0,0); // -X
  face([-x,+y,-z],[-x,+y,+z],[+x,+y,+z],[+x,+y,-z], 0,1,0);  // +Y
  face([-x,-y,+z],[-x,-y,-z],[+x,-y,-z],[+x,-y,+z], 0,-1,0); // -Y
  face([-x,-y,+z],[+x,-y,+z],[+x,+y,+z],[-x,+y,+z], 0,0,1);  // +Z
  face([+x,-y,-z],[-x,-y,-z],[-x,+y,-z],[+x,+y,-z], 0,0,-1); // -Z
  return {
    positions: new Float32Array(p),
    normals: new Float32Array(n),
    indices: new Uint32Array(idx)
  };
}

function generateSphereMesh(radius, stacks, slices) {
  var p = [], n = [], idx = [];
  var i, j, theta, phi, st, ct, sp, cp, x, y, z;
  for (i = 0; i <= stacks; i++) {
    theta = (i / stacks) * Math.PI;
    st = Math.sin(theta); ct = Math.cos(theta);
    for (j = 0; j <= slices; j++) {
      phi = (j / slices) * 2 * Math.PI;
      sp = Math.sin(phi); cp = Math.cos(phi);
      x = st * cp; y = st * sp; z = ct;
      p.push(radius*x, radius*y, radius*z);
      n.push(x, y, z);
    }
  }
  for (i = 0; i < stacks; i++) {
    for (j = 0; j < slices; j++) {
      var a = i*(slices+1)+j, b = a+slices+1;
      idx.push(a, b, a+1, a+1, b, b+1);
    }
  }
  return {
    positions: new Float32Array(p),
    normals: new Float32Array(n),
    indices: new Uint32Array(idx)
  };
}

// Generate axis lines geometry (returns {positions, colors} for GL_LINES)
function generateAxesMesh(length) {
  return {
    positions: new Float32Array([
      0,0,0, length,0,0,  // X
      0,0,0, 0,length,0,  // Y
      0,0,0, 0,0,length   // Z
    ]),
    colors: new Float32Array([
      1,0.3,0.3, 1,0.3,0.3,  // X = red
      0.3,1,0.3, 0.3,1,0.3,  // Y = green
      0.3,0.3,1, 0.3,0.3,1   // Z = blue
    ])
  };
}


// ============================================================================
// GDML Parser
// ============================================================================

function parseGDML(xmlString) {
  var parser = new DOMParser();
  var doc = parser.parseFromString(xmlString, "text/xml");
  var gdml = doc.querySelector("gdml");
  if (!gdml) throw new Error("Not a valid GDML document");

  // Parse <define> positions
  var defines = {};
  var defSection = gdml.querySelector("define");
  if (defSection) {
    var posEls = defSection.querySelectorAll("position");
    for (var i = 0; i < posEls.length; i++) {
      var pe = posEls[i];
      defines[pe.getAttribute("name")] = [
        parseFloat(pe.getAttribute("x") || "0"),
        parseFloat(pe.getAttribute("y") || "0"),
        parseFloat(pe.getAttribute("z") || "0")
      ];
    }
  }

  // Parse <materials>
  var materials = {};
  var matSection = gdml.querySelector("materials");
  if (matSection) {
    var matEls = matSection.querySelectorAll("material");
    for (var mi = 0; mi < matEls.length; mi++) {
      var mel = matEls[mi];
      var dEl = mel.querySelector("D");
      materials[mel.getAttribute("name")] = {
        name: mel.getAttribute("name"),
        density: dEl ? parseFloat(dEl.getAttribute("value")) : 0,
        state: mel.getAttribute("state") || "solid"
      };
    }
  }

  // Parse <solids>
  var solids = {};
  var solidSection = gdml.querySelector("solids");
  if (solidSection) {
    for (var si = 0; si < solidSection.children.length; si++) {
      var sel = solidSection.children[si];
      if (sel.nodeType !== 1) continue; // skip text nodes
      var sname = sel.getAttribute("name");
      var tag = sel.tagName.toLowerCase();

      if (tag === "tube") {
        solids[sname] = {
          type: "tube",
          name: sname,
          rmin: parseFloat(sel.getAttribute("rmin") || "0"),
          rmax: parseFloat(sel.getAttribute("rmax")),
          z: parseFloat(sel.getAttribute("z")),
          startphi: parseFloat(sel.getAttribute("startphi") || "0"),
          deltaphi: parseFloat(sel.getAttribute("deltaphi") || String(2*Math.PI))
        };
      } else if (tag === "polycone") {
        // Surface of revolution: coaxial and/or bulletized crystals.
        var planes = [];
        for (var zi = 0; zi < sel.children.length; zi++) {
          var zp = sel.children[zi];
          if (zp.nodeType !== 1 || zp.tagName.toLowerCase() !== "zplane") continue;
          planes.push({
            z: parseFloat(zp.getAttribute("z")),
            rmin: parseFloat(zp.getAttribute("rmin") || "0"),
            rmax: parseFloat(zp.getAttribute("rmax"))
          });
        }
        planes.sort(function (a, b) { return a.z - b.z; });
        solids[sname] = { type: "polycone", name: sname, planes: planes };
      } else if (tag === "box") {
        solids[sname] = {
          type: "box",
          name: sname,
          x: parseFloat(sel.getAttribute("x")),
          y: parseFloat(sel.getAttribute("y")),
          z: parseFloat(sel.getAttribute("z"))
        };
      } else if (tag === "sphere") {
        solids[sname] = {
          type: "sphere",
          name: sname,
          rmin: parseFloat(sel.getAttribute("rmin") || "0"),
          rmax: parseFloat(sel.getAttribute("rmax")),
          startphi: parseFloat(sel.getAttribute("startphi") || "0"),
          deltaphi: parseFloat(sel.getAttribute("deltaphi") || String(2*Math.PI)),
          starttheta: parseFloat(sel.getAttribute("starttheta") || "0"),
          deltatheta: parseFloat(sel.getAttribute("deltatheta") || String(Math.PI))
        };
      } else if (tag === "subtraction") {
        var firstRef = sel.querySelector("first").getAttribute("ref");
        var secondRef = sel.querySelector("second").getAttribute("ref");
        // Get position offset for second solid
        var subPos = [0, 0, 0];
        var posEl = sel.querySelector("position");
        var posRefEl = sel.querySelector("positionref");
        if (posEl) {
          subPos = [
            parseFloat(posEl.getAttribute("x") || "0"),
            parseFloat(posEl.getAttribute("y") || "0"),
            parseFloat(posEl.getAttribute("z") || "0")
          ];
        } else if (posRefEl) {
          var refName = posRefEl.getAttribute("ref");
          if (defines[refName]) subPos = defines[refName].slice();
        }
        solids[sname] = {
          type: "subtraction",
          name: sname,
          first: firstRef,
          second: secondRef,
          offset: subPos
        };
      }
    }
  }

  // Parse <structure>
  var volumes = {};
  var worldRef = null;
  var structSection = gdml.querySelector("structure");
  if (structSection) {
    var volEls = structSection.querySelectorAll("volume");
    for (var vi = 0; vi < volEls.length; vi++) {
      var vel = volEls[vi];
      var vname = vel.getAttribute("name");
      var matRef = vel.querySelector("materialref");
      var solRef = vel.querySelector("solidref");

      var physvols = [];
      var pvEls = vel.querySelectorAll("physvol");
      for (var pi = 0; pi < pvEls.length; pi++) {
        var pv = pvEls[pi];
        var pvName = pv.getAttribute("name") || "";
        var volRefEl = pv.querySelector("volumeref");
        var pvPos = [0, 0, 0];
        var pvPosEl = pv.querySelector("position");
        var pvPosRef = pv.querySelector("positionref");
        if (pvPosEl) {
          pvPos = [
            parseFloat(pvPosEl.getAttribute("x") || "0"),
            parseFloat(pvPosEl.getAttribute("y") || "0"),
            parseFloat(pvPosEl.getAttribute("z") || "0")
          ];
        } else if (pvPosRef) {
          var prName = pvPosRef.getAttribute("ref");
          if (defines[prName]) pvPos = defines[prName].slice();
        }
        physvols.push({
          name: pvName,
          volumeRef: volRefEl ? volRefEl.getAttribute("ref") : "",
          position: pvPos
        });
      }

      volumes[vname] = {
        name: vname,
        material: matRef ? matRef.getAttribute("ref") : "",
        solid: solRef ? solRef.getAttribute("ref") : "",
        physvols: physvols
      };
    }
  }

  // Find world volume
  var setupEl = gdml.querySelector("setup");
  if (setupEl) {
    var worldEl = setupEl.querySelector("world");
    if (worldEl) worldRef = worldEl.getAttribute("ref");
  }

  return { defines: defines, materials: materials, solids: solids, volumes: volumes, worldRef: worldRef };
}


// ============================================================================
// Subtraction decomposition -> tube primitives
// ============================================================================

function decomposeSubtraction(solids, subSolid) {
  var first = solids[subSolid.first];
  var second = solids[subSolid.second];
  if (!first || !second) return [];

  // Both should be tubes for our decomposition
  if (first.type !== "tube" || second.type !== "tube") {
    // Fallback: just render the outer solid
    return [{ solid: first, offsetZ: 0 }];
  }

  var outerR = first.rmax;
  var outerHalfZ = first.z / 2;
  var innerR = second.rmax;
  var innerHalfZ = second.z / 2;
  var dz = subSolid.offset[2]; // z-offset of inner relative to outer center
  var outerRmin = first.rmin;
  var innerRmin = second.rmin;

  // Inner cavity boundaries (in outer's local frame)
  var innerZBot = -innerHalfZ + dz;
  var innerZTop = +innerHalfZ + dz;
  var outerZBot = -outerHalfZ;
  var outerZTop = +outerHalfZ;

  var parts = [];

  if (innerR >= outerR - 0.001) {
    // Inner radius is same as outer: this is just a length difference (endcap)
    // Material at z < innerZBot or z > innerZTop
    if (innerZBot > outerZBot + 0.0001) {
      var h = innerZBot - outerZBot;
      parts.push({ type: "tube", rmin: outerRmin, rmax: outerR, halfZ: h/2,
                    offsetZ: (outerZBot + innerZBot) / 2 });
    }
    if (outerZTop > innerZTop + 0.0001) {
      var h2 = outerZTop - innerZTop;
      parts.push({ type: "tube", rmin: outerRmin, rmax: outerR, halfZ: h2/2,
                    offsetZ: (innerZTop + outerZTop) / 2 });
    }
    return parts;
  }

  // Cup or L-shape: inner removes a concentric region
  // Bottom region: from outerZBot to max(outerZBot, innerZBot) - full disk
  if (innerZBot > outerZBot + 0.0001) {
    var slabH = innerZBot - outerZBot;
    parts.push({ type: "tube", rmin: outerRmin, rmax: outerR, halfZ: slabH/2,
                 offsetZ: (outerZBot + innerZBot) / 2 });
  }

  // Side wall region: from max(outerZBot, innerZBot) to min(outerZTop, innerZTop)
  var wallBot = Math.max(outerZBot, innerZBot);
  var wallTop = Math.min(outerZTop, innerZTop);
  if (wallTop > wallBot + 0.0001) {
    var wallH = wallTop - wallBot;
    parts.push({ type: "tube", rmin: innerR, rmax: outerR, halfZ: wallH/2,
                 offsetZ: (wallBot + wallTop) / 2 });
  }

  // Top region: from min(outerZTop, innerZTop) to outerZTop - full disk
  if (outerZTop > innerZTop + 0.0001) {
    var capH = outerZTop - innerZTop;
    parts.push({ type: "tube", rmin: outerRmin, rmax: outerR, halfZ: capH/2,
                 offsetZ: (innerZTop + outerZTop) / 2 });
  }

  return parts;
}


// ============================================================================
// Volume computation
// ============================================================================

function computeSolidVolume(solids, solidName) {
  var s = solids[solidName];
  if (!s) return 0;
  if (s.type === "tube") {
    return Math.PI * (s.rmax*s.rmax - s.rmin*s.rmin) * s.z;
  } else if (s.type === "polycone") {
    // Exact for the conical frusta the planes describe.
    var v = 0;
    for (var i = 0; i + 1 < s.planes.length; i++) {
      var a = s.planes[i], b = s.planes[i+1];
      var dz = b.z - a.z;
      if (dz === 0) continue;
      var outer = a.rmax*a.rmax + a.rmax*b.rmax + b.rmax*b.rmax;
      var inner = a.rmin*a.rmin + a.rmin*b.rmin + b.rmin*b.rmin;
      v += Math.PI / 3 * dz * (outer - inner);
    }
    return v;
  } else if (s.type === "box") {
    return s.x * s.y * s.z;
  } else if (s.type === "sphere") {
    return (4/3) * Math.PI * (s.rmax*s.rmax*s.rmax - s.rmin*s.rmin*s.rmin);
  } else if (s.type === "subtraction") {
    return computeSolidVolume(solids, s.first) - computeSolidVolume(solids, s.second);
  }
  return 0;
}

function describeSolid(solids, solidName) {
  var s = solids[solidName];
  if (!s) return "";
  if (s.type === "polycone") {
    var pl = s.planes, n = pl.length;
    var zmin = pl[0].z, zmax = pl[n-1].z, rmax = 0, bore = 0, boreZ = zmax;
    for (var i = 0; i < n; i++) {
      if (pl[i].rmax > rmax) rmax = pl[i].rmax;
      if (pl[i].rmin > bore) bore = pl[i].rmin;
      if (pl[i].rmin > 1e-6 && pl[i].z < boreZ) boreZ = pl[i].z;
    }
    var d = "Polycone: R=" + rmax.toFixed(3) + " cm, L=" + (zmax - zmin).toFixed(3) + " cm";
    // The front edge is rounded when rmax at the front face is short of R.
    var frontR = pl[0].rmax;
    if (rmax - frontR > 1e-4) {
      d += ", bulletized r_b=" + (rmax - frontR).toFixed(3) + " cm";
    }
    if (bore > 1e-6) {
      d += ", bore R=" + bore.toFixed(3) + " cm x " + (zmax - boreZ).toFixed(3) + " cm deep";
    }
    return d + " (" + n + " z-planes)";
  } else if (s.type === "tube") {
    if (s.rmin < 0.001) {
      return "Cylinder: R=" + s.rmax.toFixed(3) + " cm, L=" + s.z.toFixed(3) + " cm";
    } else {
      var wallT = (s.rmax - s.rmin).toFixed(3);
      return "Tube: Rout=" + s.rmax.toFixed(3) + " cm, Rin=" + s.rmin.toFixed(3) + " cm, L=" + s.z.toFixed(3) + " cm (wall: " + wallT + " cm)";
    }
  } else if (s.type === "box") {
    return "Box: " + s.x.toFixed(3) + " x " + s.y.toFixed(3) + " x " + s.z.toFixed(3) + " cm";
  } else if (s.type === "sphere") {
    return "Sphere: R=" + s.rmax.toFixed(3) + " cm";
  } else if (s.type === "subtraction") {
    var first = solids[s.first];
    var second = solids[s.second];
    if (first && second && first.type === "tube" && second.type === "tube") {
      var sideT = first.rmax - second.rmax;
      var dz = s.offset[2] || 0;
      var frontT = (first.z - second.z) / 2 + dz;
      var backT = (first.z - second.z) / 2 - dz;
      var desc = "Shell: Rout=" + first.rmax.toFixed(3) + " cm, L=" + first.z.toFixed(3) + " cm";
      if (sideT > 0.0001) desc += ", side=" + sideT.toFixed(3) + " cm";
      if (frontT > 0.0001) desc += ", front=" + frontT.toFixed(3) + " cm";
      if (backT > 0.0001) desc += ", back=" + backT.toFixed(3) + " cm";
      return desc;
    }
    return describeSolid(solids, s.first) + " (shell)";
  }
  return "Unknown";
}


// ============================================================================
// Color scheme
// ============================================================================

var MATERIAL_COLORS = {
  'NaI':          [0.20, 0.70, 0.60],
  'LaBr3':        [0.30, 0.50, 0.80],
  'CZT':          [0.45, 0.30, 0.65],
  'Al':           [0.72, 0.72, 0.76],
  'Pb':           [0.40, 0.40, 0.44],
  'Water':        [0.30, 0.50, 0.90],
  'Soil':         [0.55, 0.40, 0.25],
  'Polyethylene': [0.82, 0.72, 0.22],
  'Pyrex':        [0.50, 0.65, 0.50],
  'Air':          null,
  'Vacuum':       null
};

function getColorForMaterial(matName) {
  if (MATERIAL_COLORS.hasOwnProperty(matName)) {
    return MATERIAL_COLORS[matName];
  }
  // Generate a deterministic color from name hash
  var hash = 0;
  for (var i = 0; i < matName.length; i++) {
    hash = ((hash << 5) - hash + matName.charCodeAt(i)) | 0;
  }
  return [
    0.3 + 0.5 * ((hash & 0xFF) / 255),
    0.3 + 0.5 * (((hash >> 8) & 0xFF) / 255),
    0.3 + 0.5 * (((hash >> 16) & 0xFF) / 255)
  ];
}

function getAlphaForVolume(volName, matName) {
  if (volName === "active_crystal") return 1.0;
  if (volName.indexOf("AttLV") === 0) return 0.45;
  if (volName.indexOf("SrcMaterial") === 0) return 0.30;
  if (volName.indexOf("SrcShield") === 0) return 0.40;
  // Default for unknown
  if (matName === "Vacuum" || matName === "Air") return 0;
  return 0.50;
}

function getLabelForVolume(volName, matName) {
  if (volName === "active_crystal") return matName + " Crystal (active)";
  if (volName.indexOf("AttLV") === 0) return matName + " Attenuator";
  if (volName === "SrcMaterialLV") return matName + " Source Material";
  if (volName.indexOf("SrcShieldOuter") === 0) return matName + " Outer Wall";
  if (volName.indexOf("SrcShieldBottom") === 0) return matName + " Bottom Wall";
  if (volName.indexOf("SrcShieldWellWall") === 0) return matName + " Well Wall";
  if (volName.indexOf("SrcShieldWellBottom") === 0) return matName + " Well Bottom";
  return volName + " (" + matName + ")";
}


// ============================================================================
// WebGL shader sources
// ============================================================================

var MAIN_VS = [
  "attribute vec3 aPosition;",
  "attribute vec3 aNormal;",
  "uniform mat4 uMVP;",
  "uniform mat4 uMV;",
  "uniform mat3 uNormalMat;",
  "varying vec3 vNormal;",
  "varying vec3 vPosition;",
  "void main() {",
  "  vNormal = normalize(uNormalMat * aNormal);",
  "  vPosition = (uMV * vec4(aPosition, 1.0)).xyz;",
  "  gl_Position = uMVP * vec4(aPosition, 1.0);",
  "}"
].join("\n");

var MAIN_FS = [
  "precision mediump float;",
  "varying vec3 vNormal;",
  "varying vec3 vPosition;",
  "uniform vec3 uColor;",
  "uniform float uAlpha;",
  "uniform float uEmissive;",
  "uniform vec3 uLightDir1;",
  "uniform vec3 uLightDir2;",
  "void main() {",
  "  vec3 N = normalize(vNormal);",
  "  if (!gl_FrontFacing) N = -N;",
  "  float diff1 = max(dot(N, uLightDir1), 0.0);",
  "  float diff2 = max(dot(N, uLightDir2), 0.0) * 0.4;",
  "  float ambient = 0.20;",
  "  vec3 viewDir = normalize(-vPosition);",
  "  vec3 halfDir = normalize(uLightDir1 + viewDir);",
  "  float spec = pow(max(dot(N, halfDir), 0.0), 48.0) * 0.35;",
  "  vec3 color = uColor * (ambient + diff1 * 0.65 + diff2) + vec3(spec) + vec3(uEmissive);",
  "  gl_FragColor = vec4(color, uAlpha);",
  "}"
].join("\n");

var PICK_VS = [
  "attribute vec3 aPosition;",
  "uniform mat4 uMVP;",
  "void main() {",
  "  gl_Position = uMVP * vec4(aPosition, 1.0);",
  "}"
].join("\n");

var PICK_FS = [
  "precision mediump float;",
  "uniform vec3 uPickColor;",
  "void main() {",
  "  gl_FragColor = vec4(uPickColor, 1.0);",
  "}"
].join("\n");

var LINE_VS = [
  "attribute vec3 aPosition;",
  "attribute vec3 aColor;",
  "uniform mat4 uMVP;",
  "varying vec3 vColor;",
  "void main() {",
  "  vColor = aColor;",
  "  gl_Position = uMVP * vec4(aPosition, 1.0);",
  "}"
].join("\n");

var LINE_FS = [
  "precision mediump float;",
  "varying vec3 vColor;",
  "void main() {",
  "  gl_FragColor = vec4(vColor, 1.0);",
  "}"
].join("\n");


// ============================================================================
// GeometryViewer3D constructor
// ============================================================================

GeometryViewer3D = function(elem, options) {
  this.container = typeof elem === "string" ? document.getElementById(elem) : elem;
  if (!this.container) throw new Error("GeometryViewer3D: invalid container element");

  this.options = options || {};
  this.options.backgroundColor = this.options.backgroundColor || [0.14, 0.14, 0.17, 1.0];
  this.options.showAxes = this.options.showAxes !== undefined ? this.options.showAxes : true;
  this.options.cylinderSegments = this.options.cylinderSegments || 48;
  this.options.sphereSegments = this.options.sphereSegments || 16;
  this.options.tooltipEnabled = this.options.tooltipEnabled !== undefined ? this.options.tooltipEnabled : true;

  // Create canvas
  this.canvas = document.createElement("canvas");
  this.canvas.style.width = "100%";
  this.canvas.style.height = "100%";
  this.canvas.style.display = "block";
  this.container.style.position = "relative";
  this.container.appendChild(this.canvas);

  // Create tooltip div
  this.tooltip = document.createElement("div");
  this.tooltip.style.cssText = "position:absolute;display:none;pointer-events:none;background:rgba(20,20,25,0.92);color:#e0e0e0;padding:8px 12px;border-radius:6px;font:13px/1.5 monospace;max-width:320px;z-index:10;border:1px solid rgba(255,255,255,0.15);box-shadow:0 4px 16px rgba(0,0,0,0.5);";
  this.container.appendChild(this.tooltip);

  // WebGL context
  this.gl = this.canvas.getContext("webgl", {
    antialias: true, alpha: false, depth: true, stencil: false,
    premultipliedAlpha: false
  });
  if (!this.gl) throw new Error("WebGL not supported");
  var gl = this.gl;

  // Compile shaders
  this.mainProgram = this._createProgram(MAIN_VS, MAIN_FS);
  this.pickProgram = this._createProgram(PICK_VS, PICK_FS);
  this.lineProgram = this._createProgram(LINE_VS, LINE_FS);

  // Get attribute/uniform locations
  this.mainLocs = {
    aPosition: gl.getAttribLocation(this.mainProgram, "aPosition"),
    aNormal:   gl.getAttribLocation(this.mainProgram, "aNormal"),
    uMVP:      gl.getUniformLocation(this.mainProgram, "uMVP"),
    uMV:       gl.getUniformLocation(this.mainProgram, "uMV"),
    uNormalMat:gl.getUniformLocation(this.mainProgram, "uNormalMat"),
    uColor:    gl.getUniformLocation(this.mainProgram, "uColor"),
    uAlpha:    gl.getUniformLocation(this.mainProgram, "uAlpha"),
    uEmissive: gl.getUniformLocation(this.mainProgram, "uEmissive"),
    uLightDir1:gl.getUniformLocation(this.mainProgram, "uLightDir1"),
    uLightDir2:gl.getUniformLocation(this.mainProgram, "uLightDir2")
  };
  this.pickLocs = {
    aPosition: gl.getAttribLocation(this.pickProgram, "aPosition"),
    uMVP:      gl.getUniformLocation(this.pickProgram, "uMVP"),
    uPickColor:gl.getUniformLocation(this.pickProgram, "uPickColor")
  };
  this.lineLocs = {
    aPosition: gl.getAttribLocation(this.lineProgram, "aPosition"),
    aColor:    gl.getAttribLocation(this.lineProgram, "aColor"),
    uMVP:      gl.getUniformLocation(this.lineProgram, "uMVP")
  };

  // Pick framebuffer
  this.pickFB = null;
  this.pickTex = null;
  this.pickDepth = null;
  this.pickWidth = 0;
  this.pickHeight = 0;

  // Scene
  this.renderables = [];    // {vao, vbo, nbo, ibo, numIndices, modelMatrix, color, alpha, label, volume, mass, dimStr, pickId, center}
  this.axesBuffers = null;
  this.sourceMarker = null; // renderable for source marker sphere
  this.sourcePosition = null;
  this.parsedGDML = null;
  this.distanceLine = null;  // {posBuffer, colBuffer, numVerts} for distance line overlay

  // Camera (free rotation via rotation matrix)
  this.camera = {
    rotMatrix: [1,0,0, 0,1,0, 0,0,1],  // 3x3 rotation matrix (column-major)
    distance: 25,
    target: [0, 0, 3],  // look-at point
    fov: 45 * Math.PI / 180,
    near: 0.1,
    far: 500
  };
  this.cameraInitialized = false;

  // Interaction state
  this.isDragging = false;
  this.isPanning = false;
  this.lastMouse = [0, 0];
  this.hoveredId = -1;
  this.dirty = true;
  this.pickDirty = true;
  this.mousePos = [0, 0];

  // Bind event handlers
  this._onMouseDown = this._handleMouseDown.bind(this);
  this._onMouseMove = this._handleMouseMove.bind(this);
  this._onMouseUp = this._handleMouseUp.bind(this);
  this._onDocMouseMove = this._handleDocMouseMove.bind(this);
  this._onDocMouseUp = this._handleDocMouseUp.bind(this);
  this._onWheel = this._handleWheel.bind(this);
  this._onContextMenu = function(e) { e.preventDefault(); };

  this.canvas.addEventListener("mousedown", this._onMouseDown);
  this.canvas.addEventListener("mousemove", this._onMouseMove);
  this.canvas.addEventListener("mouseup", this._onMouseUp);
  this.canvas.addEventListener("mouseleave", this._onMouseUp);
  this.canvas.addEventListener("wheel", this._onWheel, { passive: false });
  this.canvas.addEventListener("contextmenu", this._onContextMenu);

  // Touch events
  this._lastTouchDist = 0;
  this._lastTouchCenter = [0,0];
  this._onTouchStart = this._handleTouchStart.bind(this);
  this._onTouchMove = this._handleTouchMove.bind(this);
  this._onTouchEnd = this._handleTouchEnd.bind(this);
  this.canvas.addEventListener("touchstart", this._onTouchStart, {passive:false});
  this.canvas.addEventListener("touchmove", this._onTouchMove, {passive:false});
  this.canvas.addEventListener("touchend", this._onTouchEnd);

  // Resize observer
  this._resizeObserver = null;
  if (typeof ResizeObserver !== "undefined") {
    this._resizeObserver = new ResizeObserver(this.handleResize.bind(this));
    this._resizeObserver.observe(this.container);
  }

  // Initial size
  this.handleResize();

  // Create axes
  this._createAxes();

  // Render loop
  this._animFrame = null;
  this._renderLoop();
};

GeometryViewer3D.prototype._createProgram = function(vsSource, fsSource) {
  var gl = this.gl;
  function compile(type, src) {
    var s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
      console.error("Shader compile error:", gl.getShaderInfoLog(s));
    }
    return s;
  }
  var vs = compile(gl.VERTEX_SHADER, vsSource);
  var fs = compile(gl.FRAGMENT_SHADER, fsSource);
  var prog = gl.createProgram();
  gl.attachShader(prog, vs);
  gl.attachShader(prog, fs);
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
    console.error("Program link error:", gl.getProgramInfoLog(prog));
  }
  return prog;
};

GeometryViewer3D.prototype._createAxes = function() {
  var gl = this.gl;
  var mesh = generateAxesMesh(2.0);
  this.axesBuffers = {
    posBuffer: this._createBuffer(mesh.positions),
    colBuffer: this._createBuffer(mesh.colors),
    count: 6
  };
};

GeometryViewer3D.prototype._createBuffer = function(data) {
  var gl = this.gl;
  var buf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
  return buf;
};

GeometryViewer3D.prototype._createIndexBuffer = function(data) {
  var gl = this.gl;
  var buf = gl.createBuffer();
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buf);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, data, gl.STATIC_DRAW);
  return buf;
};


// ============================================================================
// Public API
// ============================================================================

GeometryViewer3D.prototype.loadGDML = function(xmlString) {
  // Save camera state if already initialized
  var savedCamera = this.cameraInitialized ? JSON.parse(JSON.stringify(this.camera)) : null;
  var savedBBox = this._sceneBBox ? this._sceneBBox.slice() : null;

  // Clean up old renderables
  this._cleanupRenderables();

  // Parse
  this.parsedGDML = parseGDML(xmlString);

  // Build scene
  this._buildScene();

  // Restore or reset camera
  if (savedCamera && savedBBox) {
    // Check if scene changed dramatically
    var newBBox = this._sceneBBox;
    var oldSize = Math.max(savedBBox[3]-savedBBox[0], savedBBox[4]-savedBBox[1], savedBBox[5]-savedBBox[2]);
    var newSize = Math.max(newBBox[3]-newBBox[0], newBBox[4]-newBBox[1], newBBox[5]-newBBox[2]);
    if (oldSize > 0 && (newSize / oldSize > 5 || newSize / oldSize < 0.2)) {
      this._resetCamera();
    } else {
      this.camera = savedCamera;
    }
  } else {
    this._resetCamera();
  }
  this.cameraInitialized = true;

  // Re-add source marker if set
  if (this.sourcePosition) {
    this._createSourceMarker(this.sourcePosition);
  }

  this.dirty = true;
  this.pickDirty = true;
};

GeometryViewer3D.prototype.setSourcePosition = function(x, y, z) {
  if (x === null || x === undefined) {
    this.sourcePosition = null;
    this.sourceMarker = null;
  } else {
    this.sourcePosition = [x, y, z];
    this._createSourceMarker(this.sourcePosition);
  }
  this.dirty = true;
};

GeometryViewer3D.prototype.resetCamera = function() {
  this._resetCamera();
  this.dirty = true;
};

GeometryViewer3D.prototype.setCameraPosition = function(x, y, z) {
  var t = this.camera.target;
  var d = vec3_sub([x,y,z], t);
  this.camera.distance = vec3_length(d);
  if (this.camera.distance < 0.01) this.camera.distance = 0.01;
  // Build rotation matrix from direction vector
  var fwd = vec3_normalize(d);  // camera forward = from target to eye
  // Choose an up hint that isn't parallel to fwd
  var upHint = (Math.abs(fwd[1]) < 0.9) ? [0,1,0] : [1,0,0];
  var right = vec3_normalize(vec3_cross(upHint, fwd));
  var up = vec3_cross(fwd, right);
  this.camera.rotMatrix = [right[0],right[1],right[2], up[0],up[1],up[2], fwd[0],fwd[1],fwd[2]];
  this.dirty = true;
};

GeometryViewer3D.prototype.getCameraPosition = function() {
  var eye = this._getEyePosition();
  return {x: eye[0], y: eye[1], z: eye[2]};
};

GeometryViewer3D.prototype.handleResize = function() {
  var rect = this.container.getBoundingClientRect();
  var dpr = window.devicePixelRatio || 1;
  var w = Math.max(1, Math.floor(rect.width * dpr));
  var h = Math.max(1, Math.floor(rect.height * dpr));
  if (this.canvas.width !== w || this.canvas.height !== h) {
    this.canvas.width = w;
    this.canvas.height = h;
    this.gl.viewport(0, 0, w, h);
    this.dirty = true;
    this.pickDirty = true;
  }
};

GeometryViewer3D.prototype.destroy = function() {
  if (this._animFrame) {
    cancelAnimationFrame(this._animFrame);
    this._animFrame = null;
  }
  if (this._resizeObserver) {
    this._resizeObserver.disconnect();
    this._resizeObserver = null;
  }
  this.canvas.removeEventListener("mousedown", this._onMouseDown);
  this.canvas.removeEventListener("mousemove", this._onMouseMove);
  this.canvas.removeEventListener("mouseup", this._onMouseUp);
  this.canvas.removeEventListener("mouseleave", this._onMouseUp);
  this.canvas.removeEventListener("wheel", this._onWheel);
  this.canvas.removeEventListener("contextmenu", this._onContextMenu);
  document.removeEventListener("mousemove", this._onDocMouseMove);
  document.removeEventListener("mouseup", this._onDocMouseUp);
  this.canvas.removeEventListener("touchstart", this._onTouchStart);
  this.canvas.removeEventListener("touchmove", this._onTouchMove);
  this.canvas.removeEventListener("touchend", this._onTouchEnd);
  this._cleanupRenderables();
  if (this.canvas.parentNode) this.canvas.parentNode.removeChild(this.canvas);
  if (this.tooltip.parentNode) this.tooltip.parentNode.removeChild(this.tooltip);
};


// ============================================================================
// Scene building
// ============================================================================

GeometryViewer3D.prototype._cleanupRenderables = function() {
  var gl = this.gl;
  for (var i = 0; i < this.renderables.length; i++) {
    var r = this.renderables[i];
    if (r.posBuffer) gl.deleteBuffer(r.posBuffer);
    if (r.normBuffer) gl.deleteBuffer(r.normBuffer);
    if (r.idxBuffer) gl.deleteBuffer(r.idxBuffer);
  }
  this.renderables = [];
  this.sourceMarker = null;
  this._updateDistanceLine(null);
};

GeometryViewer3D.prototype._buildScene = function() {
  var parsed = this.parsedGDML;
  if (!parsed || !parsed.worldRef) return;

  var worldVol = parsed.volumes[parsed.worldRef];
  if (!worldVol) return;

  var segs = this.options.cylinderSegments;
  var self = this;
  var pickIdCounter = 1;

  // Compute bounding box
  var bbMin = [Infinity, Infinity, Infinity];
  var bbMax = [-Infinity, -Infinity, -Infinity];

  worldVol.physvols.forEach(function(pv) {
    var vol = parsed.volumes[pv.volumeRef];
    if (!vol) return;
    var matName = vol.material;
    if (matName === "Vacuum" || matName === "Air") return;

    var color = getColorForMaterial(matName);
    if (!color) return;
    var alpha = getAlphaForVolume(vol.name, matName);
    var label = getLabelForVolume(vol.name, matName);
    var solidName = vol.solid;
    var solid = parsed.solids[solidName];
    if (!solid) return;

    var volume = computeSolidVolume(parsed.solids, solidName);
    var density = parsed.materials[matName] ? parsed.materials[matName].density : 0;
    var mass = volume * density;
    var dimStr = describeSolid(parsed.solids, solidName);

    // Generate mesh(es) for this physvol
    var meshes = self._generateMeshesForSolid(solid, parsed.solids, segs);

    meshes.forEach(function(meshInfo) {
      var mesh = meshInfo.mesh;
      var localOffsetZ = meshInfo.offsetZ || 0;

      // Model matrix: translate to physvol position + local offset
      var modelMat = mat4_translate(pv.position[0], pv.position[1], pv.position[2] + localOffsetZ);

      // Compute center for sorting
      var center = [pv.position[0], pv.position[1], pv.position[2] + localOffsetZ];

      // Update bounding box
      // Approximate with center and max extent
      var maxExtent = 0;
      for (var vi = 0; vi < mesh.positions.length; vi += 3) {
        var vx = mesh.positions[vi] + center[0];
        var vy = mesh.positions[vi+1] + center[1];
        var vz = mesh.positions[vi+2] + center[2];
        if (vx < bbMin[0]) bbMin[0] = vx;
        if (vy < bbMin[1]) bbMin[1] = vy;
        if (vz < bbMin[2]) bbMin[2] = vz;
        if (vx > bbMax[0]) bbMax[0] = vx;
        if (vy > bbMax[1]) bbMax[1] = vy;
        if (vz > bbMax[2]) bbMax[2] = vz;
      }

      var renderable = {
        posBuffer: self._createBuffer(mesh.positions),
        normBuffer: self._createBuffer(mesh.normals),
        idxBuffer: self._createIndexBuffer(mesh.indices),
        numIndices: mesh.indices.length,
        modelMatrix: modelMat,
        color: color,
        alpha: alpha,
        label: label,
        materialName: matName,
        density: density,
        volume: volume,
        mass: mass,
        dimStr: dimStr,
        pickId: pickIdCounter++,
        center: center,
        volName: vol.name
      };
      self.renderables.push(renderable);
    });
  });

  // Fallback bounding box
  if (bbMin[0] === Infinity) {
    bbMin = [-5,-5,-5];
    bbMax = [5,5,5];
  }
  this._sceneBBox = [bbMin[0], bbMin[1], bbMin[2], bbMax[0], bbMax[1], bbMax[2]];
};

GeometryViewer3D.prototype._generateMeshesForSolid = function(solid, allSolids, segs) {
  if (solid.type === "tube") {
    return [{ mesh: generateTubeMesh(solid.rmin, solid.rmax, solid.z/2, segs), offsetZ: 0 }];
  } else if (solid.type === "polycone") {
    return [{ mesh: generatePolyconeMesh(solid.planes, segs), offsetZ: 0 }];
  } else if (solid.type === "box") {
    return [{ mesh: generateBoxMesh(solid.x/2, solid.y/2, solid.z/2), offsetZ: 0 }];
  } else if (solid.type === "sphere") {
    return [{ mesh: generateSphereMesh(solid.rmax, this.options.sphereSegments, this.options.sphereSegments*2), offsetZ: 0 }];
  } else if (solid.type === "subtraction") {
    var parts = decomposeSubtraction(allSolids, solid);
    var result = [];
    for (var i = 0; i < parts.length; i++) {
      var p = parts[i];
      if (p.solid) {
        // Fallback: render the first solid directly
        var fm = this._generateMeshesForSolid(p.solid, allSolids, segs);
        for (var j = 0; j < fm.length; j++) {
          fm[j].offsetZ += p.offsetZ;
        }
        result = result.concat(fm);
      } else {
        // Decomposed tube part
        var mesh = generateTubeMesh(p.rmin, p.rmax, p.halfZ, segs);
        result.push({ mesh: mesh, offsetZ: p.offsetZ });
      }
    }
    return result;
  }
  return [];
};

GeometryViewer3D.prototype._createSourceMarker = function(pos) {
  var mesh = generateSphereMesh(0.3, 12, 24);
  this.sourceMarker = {
    posBuffer: this._createBuffer(mesh.positions),
    normBuffer: this._createBuffer(mesh.normals),
    idxBuffer: this._createIndexBuffer(mesh.indices),
    numIndices: mesh.indices.length,
    modelMatrix: mat4_translate(pos[0], pos[1], pos[2]),
    color: [1.0, 0.2, 0.2],
    alpha: 1.0,
    label: "Point Source",
    materialName: "",
    density: 0,
    volume: 0,
    mass: 0,
    dimStr: "Position: (" + pos[0].toFixed(1) + ", " + pos[1].toFixed(1) + ", " + pos[2].toFixed(1) + ") cm",
    pickId: 999,
    center: pos.slice(),
    volName: "source"
  };
};


// ============================================================================
// Distance measurement (source to detector)
// ============================================================================

GeometryViewer3D.prototype._getDetectorFaces = function() {
  if (!this.parsedGDML) return null;
  var parsed = this.parsedGDML;
  var worldVol = parsed.volumes[parsed.worldRef];
  if (!worldVol) return null;

  var crystalFaceCenter = null;
  var detGeomFaceZ = null;
  var detGeomFaceCenter = null;

  for (var i = 0; i < worldVol.physvols.length; i++) {
    var pv = worldVol.physvols[i];
    var vol = parsed.volumes[pv.volumeRef];
    if (!vol) continue;

    if (vol.name === "active_crystal") {
      var solid = parsed.solids[vol.solid];
      if (solid) {
        var halfZ = 0;
        if (solid.type === "tube") halfZ = solid.z / 2;
        else if (solid.type === "box") halfZ = solid.z / 2;
        else if (solid.type === "polycone") {
          halfZ = (solid.planes[solid.planes.length-1].z - solid.planes[0].z) / 2;
        }
        var fz = pv.position[2] - halfZ;
        crystalFaceCenter = [pv.position[0], pv.position[1], fz];
      }
    }

    if (vol.name.indexOf("AttLV") === 0) {
      var attSolid = parsed.solids[vol.solid];
      var attHalfZ = null;
      if (attSolid && attSolid.type === "subtraction") {
        var first = parsed.solids[attSolid.first];
        if (first && first.type === "tube") attHalfZ = first.z / 2;
      } else if (attSolid && attSolid.type === "tube") {
        attHalfZ = attSolid.z / 2;
      }
      if (attHalfZ !== null) {
        var afz = pv.position[2] - attHalfZ;
        if (detGeomFaceZ === null || afz < detGeomFaceZ) {
          detGeomFaceZ = afz;
          detGeomFaceCenter = [pv.position[0], pv.position[1], afz];
        }
      }
    }
  }

  if (!crystalFaceCenter) return null;

  return {
    crystalFaceCenter: crystalFaceCenter,
    detGeomFaceCenter: detGeomFaceCenter,
    hasAttenuator: detGeomFaceCenter !== null && detGeomFaceZ < crystalFaceCenter[2] - 0.001
  };
};

GeometryViewer3D.prototype._getSourceSurfaceClosestZ = function() {
  // Find the z-coordinate of the source material surface closest to crystal (on-axis)
  if (!this.parsedGDML) return null;
  var parsed = this.parsedGDML;
  var worldVol = parsed.volumes[parsed.worldRef];
  if (!worldVol) return null;

  for (var i = 0; i < worldVol.physvols.length; i++) {
    var pv = worldVol.physvols[i];
    var vol = parsed.volumes[pv.volumeRef];
    if (!vol || vol.name !== "SrcMaterialLV") continue;

    var solid = parsed.solids[vol.solid];
    if (!solid) return null;

    if (solid.type === "subtraction") {
      // Marinelli L-shape: on-axis closest source surface is at the well bottom
      var second = parsed.solids[solid.second];
      if (second && second.type === "tube") {
        return pv.position[2] + solid.offset[2] - second.z / 2;
      }
    } else if (solid.type === "tube") {
      // Simple tube source: top face closest to crystal
      return pv.position[2] + solid.z / 2;
    }
    return null;
  }
  return null;
};

GeometryViewer3D.prototype._getSourceCenter = function() {
  if (!this.parsedGDML) return null;
  var parsed = this.parsedGDML;
  var worldVol = parsed.volumes[parsed.worldRef];
  if (!worldVol) return null;

  for (var i = 0; i < worldVol.physvols.length; i++) {
    var pv = worldVol.physvols[i];
    var vol = parsed.volumes[pv.volumeRef];
    if (vol && vol.name === "SrcMaterialLV") {
      return pv.position.slice();
    }
  }
  return null;
};

GeometryViewer3D.prototype._isMarinelliSource = function() {
  // Check if source material uses a subtraction solid (Marinelli L-shape)
  if (!this.parsedGDML) return false;
  var parsed = this.parsedGDML;
  var worldVol = parsed.volumes[parsed.worldRef];
  if (!worldVol) return false;
  for (var i = 0; i < worldVol.physvols.length; i++) {
    var pv = worldVol.physvols[i];
    var vol = parsed.volumes[pv.volumeRef];
    if (vol && vol.name === "SrcMaterialLV") {
      var solid = parsed.solids[vol.solid];
      return solid && solid.type === "subtraction";
    }
  }
  return false;
};

GeometryViewer3D.prototype._computeDistanceInfo = function(hoveredObj) {
  var isPointSource = (hoveredObj.volName === "source");
  var isExtendedSource = (hoveredObj.volName === "SrcMaterialLV");
  if (!isPointSource && !isExtendedSource) return null;

  // Skip Marinelli beaker sources (geometry wraps around detector, distances not meaningful)
  if (isExtendedSource && this._isMarinelliSource()) return null;

  var faces = this._getDetectorFaces();
  if (!faces) return null;

  var result = {
    lineStart: null,
    lineEnd: faces.crystalFaceCenter,
    distances: []
  };

  var sourceCenter;
  if (isPointSource) {
    sourceCenter = this.sourcePosition ? this.sourcePosition.slice() : null;
  } else {
    sourceCenter = this._getSourceCenter();
  }
  if (!sourceCenter) return null;

  result.lineStart = sourceCenter;

  // Source center to crystal face
  var d1 = vec3_length(vec3_sub(sourceCenter, faces.crystalFaceCenter));
  result.distances.push({ label: "Source center \u2192 crystal face", value: d1 });

  // Source center to detector housing face
  if (faces.hasAttenuator) {
    var d2 = vec3_length(vec3_sub(sourceCenter, faces.detGeomFaceCenter));
    result.distances.push({ label: "Source center \u2192 housing face", value: d2 });
  }

  // Extended source: surface distances
  if (isExtendedSource) {
    var surfaceZ = this._getSourceSurfaceClosestZ();
    if (surfaceZ !== null) {
      var surfacePt = [faces.crystalFaceCenter[0], faces.crystalFaceCenter[1], surfaceZ];

      var d3 = vec3_length(vec3_sub(surfacePt, faces.crystalFaceCenter));
      result.distances.push({ label: "Source surface \u2192 crystal face", value: d3 });

      if (faces.hasAttenuator) {
        var d4 = vec3_length(vec3_sub(surfacePt, faces.detGeomFaceCenter));
        result.distances.push({ label: "Source surface \u2192 housing face", value: d4 });
      }
    }
  }

  return result;
};

GeometryViewer3D.prototype._updateDistanceLine = function(info) {
  var gl = this.gl;

  // Clean up old buffers
  if (this.distanceLine) {
    gl.deleteBuffer(this.distanceLine.posBuffer);
    gl.deleteBuffer(this.distanceLine.colBuffer);
    this.distanceLine = null;
    this.dirty = true;
  }

  if (!info || !info.lineStart || !info.lineEnd) return;

  var s = info.lineStart;
  var e = info.lineEnd;
  var positions = new Float32Array([s[0],s[1],s[2], e[0],e[1],e[2]]);
  var colors = new Float32Array([1.0,0.9,0.3, 1.0,0.9,0.3]);

  this.distanceLine = {
    posBuffer: this._createBuffer(positions),
    colBuffer: this._createBuffer(colors),
    numVerts: 2
  };
  this.dirty = true;
};

GeometryViewer3D.prototype._renderDistanceLine = function(vp) {
  var gl = this.gl;
  gl.useProgram(this.lineProgram);
  gl.uniformMatrix4fv(this.lineLocs.uMVP, false, vp);

  gl.bindBuffer(gl.ARRAY_BUFFER, this.distanceLine.posBuffer);
  gl.enableVertexAttribArray(this.lineLocs.aPosition);
  gl.vertexAttribPointer(this.lineLocs.aPosition, 3, gl.FLOAT, false, 0, 0);

  gl.bindBuffer(gl.ARRAY_BUFFER, this.distanceLine.colBuffer);
  gl.enableVertexAttribArray(this.lineLocs.aColor);
  gl.vertexAttribPointer(this.lineLocs.aColor, 3, gl.FLOAT, false, 0, 0);

  gl.lineWidth(2.0);
  gl.drawArrays(gl.LINES, 0, this.distanceLine.numVerts);
};


// ============================================================================
// 3x3 rotation helpers
// ============================================================================

function mat3_identity() { return [1,0,0, 0,1,0, 0,0,1]; }

function mat3_multiply(a, b) {
  // Column-major 3x3: a[col*3+row]
  var r = new Array(9);
  for (var c = 0; c < 3; c++) {
    for (var row = 0; row < 3; row++) {
      r[c*3+row] = a[row]*b[c*3] + a[3+row]*b[c*3+1] + a[6+row]*b[c*3+2];
    }
  }
  return r;
}

function mat3_fromAxisAngle(axis, angle) {
  var c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
  var x = axis[0], y = axis[1], z = axis[2];
  return [
    t*x*x + c,   t*x*y + s*z, t*x*z - s*y,
    t*x*y - s*z, t*y*y + c,   t*y*z + s*x,
    t*x*z + s*y, t*y*z - s*x, t*z*z + c
  ];
}

function mat3_transformVec(m, v) {
  return [
    m[0]*v[0] + m[3]*v[1] + m[6]*v[2],
    m[1]*v[0] + m[4]*v[1] + m[7]*v[2],
    m[2]*v[0] + m[5]*v[1] + m[8]*v[2]
  ];
}

// Re-orthogonalize a 3x3 rotation matrix to prevent drift
function mat3_orthogonalize(m) {
  // Gram-Schmidt on columns
  var c0 = vec3_normalize([m[0], m[1], m[2]]);
  var c1 = [m[3], m[4], m[5]];
  c1 = vec3_sub(c1, vec3_scale(c0, vec3_dot(c0, c1)));
  c1 = vec3_normalize(c1);
  var c2 = vec3_cross(c0, c1);
  return [c0[0],c0[1],c0[2], c1[0],c1[1],c1[2], c2[0],c2[1],c2[2]];
}


// ============================================================================
// Camera
// ============================================================================

GeometryViewer3D.prototype._getEyePosition = function() {
  var c = this.camera;
  // Camera looks along -Z in its local frame. Eye = target + R * [0,0,distance]
  var fwd = mat3_transformVec(c.rotMatrix, [0, 0, c.distance]);
  return vec3_add(c.target, fwd);
};

GeometryViewer3D.prototype._getCameraUp = function() {
  return mat3_transformVec(this.camera.rotMatrix, [0, 1, 0]);
};

GeometryViewer3D.prototype._getCameraRight = function() {
  return mat3_transformVec(this.camera.rotMatrix, [1, 0, 0]);
};

GeometryViewer3D.prototype._getViewMatrix = function() {
  var eye = this._getEyePosition();
  var up = this._getCameraUp();
  return mat4_lookAt(eye, this.camera.target, up);
};

GeometryViewer3D.prototype._getProjectionMatrix = function() {
  var aspect = this.canvas.width / Math.max(1, this.canvas.height);
  return mat4_perspective(this.camera.fov, aspect, this.camera.near, this.camera.far);
};

GeometryViewer3D.prototype._rotateCamera = function(dx, dy) {
  // Rotate around screen axes: dx -> rotate around screen up, dy -> rotate around screen right
  var sensitivity = 0.006;
  var up = this._getCameraUp();
  var right = this._getCameraRight();

  // Apply horizontal rotation (around screen up)
  if (Math.abs(dx) > 0.01) {
    var Rh = mat3_fromAxisAngle(up, -dx * sensitivity);
    this.camera.rotMatrix = mat3_multiply(Rh, this.camera.rotMatrix);
  }
  // Apply vertical rotation (around screen right)
  if (Math.abs(dy) > 0.01) {
    var Rv = mat3_fromAxisAngle(right, -dy * sensitivity);
    this.camera.rotMatrix = mat3_multiply(Rv, this.camera.rotMatrix);
  }
  // Re-orthogonalize periodically to prevent drift
  this.camera.rotMatrix = mat3_orthogonalize(this.camera.rotMatrix);
};

GeometryViewer3D.prototype._resetCamera = function() {
  if (!this._sceneBBox) {
    this.camera.target = [0, 0, 3];
    this.camera.distance = 25;
    // Default rotation: looking from behind-right, slightly above
    this._setDefaultRotation();
    return;
  }
  var bb = this._sceneBBox;

  // Include source position in bounding box
  if (this.sourcePosition) {
    var sp = this.sourcePosition;
    bb = [
      Math.min(bb[0], sp[0]-0.5), Math.min(bb[1], sp[1]-0.5), Math.min(bb[2], sp[2]-0.5),
      Math.max(bb[3], sp[0]+0.5), Math.max(bb[4], sp[1]+0.5), Math.max(bb[5], sp[2]+0.5)
    ];
  }

  var cx = (bb[0]+bb[3])/2, cy = (bb[1]+bb[4])/2, cz = (bb[2]+bb[5])/2;
  this.camera.target = [cx, cy, cz];

  var dx = bb[3]-bb[0], dy = bb[4]-bb[1], dz = bb[5]-bb[2];
  var maxDim = Math.max(dx, dy, dz);
  var halfFov = this.camera.fov / 2;
  this.camera.distance = (maxDim / 2) / Math.tan(halfFov) * 1.3;

  this._setDefaultRotation();
};

GeometryViewer3D.prototype._setDefaultRotation = function() {
  // Set camera to view from behind-right, slightly above the detector.
  // The detector extends along +Z, source is at -Z.
  // We want to look from a direction that shows both the source side (front face)
  // and the side of the detector.
  // Start with identity (looking along -Z), then rotate:
  //   - tilt up by ~25 degrees around X to see from slightly above
  //   - rotate left by ~25 degrees around Y to see from the side
  var Rx = mat3_fromAxisAngle([1,0,0], -0.40);  // tilt up
  var Ry = mat3_fromAxisAngle([0,1,0], 0.45);   // swing right
  this.camera.rotMatrix = mat3_multiply(Ry, Rx);
};


// ============================================================================
// Input handlers
// ============================================================================

GeometryViewer3D.prototype._handleMouseDown = function(e) {
  e.preventDefault();
  if (e.button === 0) {
    this.isDragging = true;
  } else if (e.button === 2 || e.button === 1) {
    this.isPanning = true;
  }
  this.lastMouse = [e.clientX, e.clientY];
  // Listen on document so drag continues even if cursor leaves canvas
  document.addEventListener("mousemove", this._onDocMouseMove);
  document.addEventListener("mouseup", this._onDocMouseUp);
};

GeometryViewer3D.prototype._handleMouseMove = function(e) {
  var rect = this.canvas.getBoundingClientRect();
  this.mousePos = [e.clientX - rect.left, e.clientY - rect.top];

  if (this.isDragging) {
    var dx = e.clientX - this.lastMouse[0];
    var dy = e.clientY - this.lastMouse[1];
    this._rotateCamera(dx, dy);
    this.lastMouse = [e.clientX, e.clientY];
    this.dirty = true;
    this.pickDirty = true;
  } else if (this.isPanning) {
    var dx2 = e.clientX - this.lastMouse[0];
    var dy2 = e.clientY - this.lastMouse[1];
    var panScale = this.camera.distance * 0.002;
    var right = this._getCameraRight();
    var up = this._getCameraUp();
    this.camera.target = vec3_add(this.camera.target, vec3_add(
      vec3_scale(right, -dx2 * panScale),
      vec3_scale(up, dy2 * panScale)
    ));
    this.lastMouse = [e.clientX, e.clientY];
    this.dirty = true;
    this.pickDirty = true;
  } else {
    // Hover pick only
    this.pickDirty = true;
  }
};

GeometryViewer3D.prototype._handleDocMouseMove = function(e) {
  if (this.isDragging) {
    var dx = e.clientX - this.lastMouse[0];
    var dy = e.clientY - this.lastMouse[1];
    this._rotateCamera(dx, dy);
    this.lastMouse = [e.clientX, e.clientY];
    this.dirty = true;
    this.pickDirty = true;
  } else if (this.isPanning) {
    var dx2 = e.clientX - this.lastMouse[0];
    var dy2 = e.clientY - this.lastMouse[1];
    var panScale = this.camera.distance * 0.002;
    var right = this._getCameraRight();
    var up = this._getCameraUp();
    this.camera.target = vec3_add(this.camera.target, vec3_add(
      vec3_scale(right, -dx2 * panScale),
      vec3_scale(up, dy2 * panScale)
    ));
    this.lastMouse = [e.clientX, e.clientY];
    this.dirty = true;
    this.pickDirty = true;
  }
};

GeometryViewer3D.prototype._handleDocMouseUp = function(e) {
  this.isDragging = false;
  this.isPanning = false;
  document.removeEventListener("mousemove", this._onDocMouseMove);
  document.removeEventListener("mouseup", this._onDocMouseUp);
};

GeometryViewer3D.prototype._handleMouseUp = function(e) {
  this.isDragging = false;
  this.isPanning = false;
  // Hide tooltip on mouse leave
  if (e && e.type === "mouseleave") {
    this.tooltip.style.display = "none";
    if (this.hoveredId !== -1) {
      this.hoveredId = -1;
      this.dirty = true;
    }
  }
};

GeometryViewer3D.prototype._handleWheel = function(e) {
  e.preventDefault();
  var factor = e.deltaY > 0 ? 1.08 : 0.92;
  this.camera.distance *= factor;
  this.camera.distance = Math.max(0.5, Math.min(500, this.camera.distance));
  this.dirty = true;
  this.pickDirty = true;
};

GeometryViewer3D.prototype._handleTouchStart = function(e) {
  e.preventDefault();
  if (e.touches.length === 1) {
    this.isDragging = true;
    this.lastMouse = [e.touches[0].clientX, e.touches[0].clientY];
  } else if (e.touches.length === 2) {
    this.isDragging = false;
    var dx = e.touches[1].clientX - e.touches[0].clientX;
    var dy = e.touches[1].clientY - e.touches[0].clientY;
    this._lastTouchDist = Math.sqrt(dx*dx + dy*dy);
    this._lastTouchCenter = [(e.touches[0].clientX + e.touches[1].clientX)/2,
                             (e.touches[0].clientY + e.touches[1].clientY)/2];
  }
};

GeometryViewer3D.prototype._handleTouchMove = function(e) {
  e.preventDefault();
  if (e.touches.length === 1 && this.isDragging) {
    var dx = e.touches[0].clientX - this.lastMouse[0];
    var dy = e.touches[0].clientY - this.lastMouse[1];
    this._rotateCamera(dx, dy);
    this.lastMouse = [e.touches[0].clientX, e.touches[0].clientY];
    this.dirty = true;
  } else if (e.touches.length === 2) {
    var dx2 = e.touches[1].clientX - e.touches[0].clientX;
    var dy2 = e.touches[1].clientY - e.touches[0].clientY;
    var dist = Math.sqrt(dx2*dx2 + dy2*dy2);
    if (this._lastTouchDist > 0) {
      var factor = this._lastTouchDist / dist;
      this.camera.distance *= factor;
      this.camera.distance = Math.max(0.5, Math.min(500, this.camera.distance));
      this.dirty = true;
    }
    this._lastTouchDist = dist;
  }
};

GeometryViewer3D.prototype._handleTouchEnd = function(e) {
  this.isDragging = false;
  this._lastTouchDist = 0;
};


// ============================================================================
// Picking
// ============================================================================

GeometryViewer3D.prototype._ensurePickFB = function() {
  var gl = this.gl;
  var w = Math.max(1, Math.floor(this.canvas.width / 2));
  var h = Math.max(1, Math.floor(this.canvas.height / 2));
  if (this.pickFB && this.pickWidth === w && this.pickHeight === h) return;

  if (this.pickFB) {
    gl.deleteFramebuffer(this.pickFB);
    gl.deleteTexture(this.pickTex);
    gl.deleteRenderbuffer(this.pickDepth);
  }

  this.pickWidth = w;
  this.pickHeight = h;
  this.pickTex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, this.pickTex);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, w, h, 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

  this.pickDepth = gl.createRenderbuffer();
  gl.bindRenderbuffer(gl.RENDERBUFFER, this.pickDepth);
  gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, w, h);

  this.pickFB = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFB);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.pickTex, 0);
  gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, this.pickDepth);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
};

GeometryViewer3D.prototype._renderPick = function() {
  var gl = this.gl;
  this._ensurePickFB();

  gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFB);
  gl.viewport(0, 0, this.pickWidth, this.pickHeight);
  gl.clearColor(0, 0, 0, 1);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.disable(gl.BLEND);
  gl.enable(gl.DEPTH_TEST);
  gl.depthMask(true);

  gl.useProgram(this.pickProgram);
  var proj = this._getProjectionMatrix();
  var view = this._getViewMatrix();
  var vp = mat4_multiply(proj, view);

  var allObjs = this.renderables.slice();
  if (this.sourceMarker) allObjs.push(this.sourceMarker);

  for (var i = 0; i < allObjs.length; i++) {
    var r = allObjs[i];
    if (r.alpha < 0.01) continue;

    var mvp = mat4_multiply(vp, r.modelMatrix);
    gl.uniformMatrix4fv(this.pickLocs.uMVP, false, mvp);

    // Encode pickId as RGB
    var id = r.pickId;
    gl.uniform3f(this.pickLocs.uPickColor,
      ((id >> 0) & 0xFF) / 255,
      ((id >> 8) & 0xFF) / 255,
      ((id >> 16) & 0xFF) / 255
    );

    gl.bindBuffer(gl.ARRAY_BUFFER, r.posBuffer);
    gl.enableVertexAttribArray(this.pickLocs.aPosition);
    gl.vertexAttribPointer(this.pickLocs.aPosition, 3, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, r.idxBuffer);
    gl.drawElements(gl.TRIANGLES, r.numIndices, gl.UNSIGNED_INT, 0);
  }

  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.viewport(0, 0, this.canvas.width, this.canvas.height);
};

GeometryViewer3D.prototype._doPick = function() {
  if (!this.pickDirty || this.isDragging || this.isPanning) return;
  this.pickDirty = false;

  // Check for OES_element_index_uint
  this._renderPick();

  var gl = this.gl;
  var dpr = window.devicePixelRatio || 1;
  var px = Math.floor(this.mousePos[0] * dpr / 2);
  var py = Math.floor((this.canvas.height / dpr - this.mousePos[1]) * dpr / 2);

  gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFB);
  var pixel = new Uint8Array(4);
  gl.readPixels(px, py, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixel);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);

  var id = pixel[0] + (pixel[1] << 8) + (pixel[2] << 16);
  if (id === 0) id = -1;

  if (id !== this.hoveredId) {
    this.hoveredId = id;
    this.dirty = true;
    this._updateTooltip();
  }
};

GeometryViewer3D.prototype._formatTooltipEntry = function(obj, isMain) {
  var html = '';
  var nameStyle = isMain ? 'font-weight:bold;margin-bottom:4px;color:#fff;' : 'font-weight:bold;margin-top:4px;color:#ccc;font-size:12px;';
  html += '<div style="' + nameStyle + '">' + obj.label + '</div>';
  var detailStyle = isMain ? '' : 'font-size:12px;color:#aaa;';
  if (obj.materialName) {
    html += '<div style="' + detailStyle + '">Material: ' + obj.materialName + ' (' + obj.density.toFixed(2) + ' g/cm\u00B3)</div>';
  }
  if (obj.dimStr) {
    html += '<div style="' + detailStyle + '">' + obj.dimStr + '</div>';
  }
  if (obj.volume > 0) {
    html += '<div style="' + detailStyle + '">Volume: ' + obj.volume.toFixed(2) + ' cm\u00B3, Mass: ' + obj.mass.toFixed(1) + ' g</div>';
  }
  return html;
};

GeometryViewer3D.prototype._findContainedObjects = function(outer) {
  // Find renderables whose center is geometrically inside the outer object's extent.
  // For shell/attenuator volumes, this means objects that the shell wraps around.
  // We use the parsed GDML solid info to get the outer object's radial and z extent,
  // then check if other objects' centers fall inside.
  var results = [];
  var seen = {};
  var oc = outer.center;

  // Get the outer object's bounding extent from its solid
  var outerExtent = this._getSolidExtent(outer);

  for (var i = 0; i < this.renderables.length; i++) {
    var r = this.renderables[i];
    if (r.pickId === outer.pickId) continue;
    if (r.volName === outer.volName) continue;
    if (seen[r.volName]) continue;
    if (r.volName === "source") continue;  // skip source marker

    var ic = r.center;
    // Check if this object's center is within the outer's bounding cylinder/box
    var relX = ic[0] - oc[0], relY = ic[1] - oc[1], relZ = ic[2] - oc[2];

    var inside = false;
    if (outerExtent.type === "cylinder") {
      var rDist = Math.sqrt(relX*relX + relY*relY);
      inside = (rDist < outerExtent.radius + 0.1) && (Math.abs(relZ) < outerExtent.halfZ + 0.1);
    } else if (outerExtent.type === "box") {
      inside = (Math.abs(relX) < outerExtent.halfX + 0.1) &&
               (Math.abs(relY) < outerExtent.halfY + 0.1) &&
               (Math.abs(relZ) < outerExtent.halfZ + 0.1);
    } else {
      // Fallback: distance-based
      var dist = Math.sqrt(relX*relX + relY*relY + relZ*relZ);
      inside = (dist < outerExtent.radius + 0.1);
    }

    if (inside) {
      seen[r.volName] = true;
      results.push(r);
    }
  }

  // Sort by volume (largest first)
  results.sort(function(a, b) { return b.volume - a.volume; });
  return results;
};

GeometryViewer3D.prototype._getSolidExtent = function(renderable) {
  // Get the bounding extent of a renderable's solid from parsed GDML
  if (!this.parsedGDML) return { type: "sphere", radius: 5 };
  var parsed = this.parsedGDML;
  var vol = parsed.volumes[renderable.volName];
  if (!vol) return { type: "sphere", radius: 5 };
  var solid = parsed.solids[vol.solid];
  if (!solid) return { type: "sphere", radius: 5 };

  if (solid.type === "tube") {
    return { type: "cylinder", radius: solid.rmax, halfZ: solid.z / 2 };
  } else if (solid.type === "polycone") {
    // Bounding cylinder: bulletization only removes material inside it.
    var pr = 0, pl = solid.planes;
    for (var pi = 0; pi < pl.length; pi++) if (pl[pi].rmax > pr) pr = pl[pi].rmax;
    return { type: "cylinder", radius: pr,
             halfZ: (pl[pl.length-1].z - pl[0].z) / 2 };
  } else if (solid.type === "box") {
    return { type: "box", halfX: solid.x / 2, halfY: solid.y / 2, halfZ: solid.z / 2 };
  } else if (solid.type === "subtraction") {
    var first = parsed.solids[solid.first];
    if (first && first.type === "tube") {
      return { type: "cylinder", radius: first.rmax, halfZ: first.z / 2 };
    }
  }
  return { type: "sphere", radius: 5 };
};

GeometryViewer3D.prototype._updateTooltip = function() {
  if (!this.options.tooltipEnabled) {
    this.tooltip.style.display = "none";
    this._updateDistanceLine(null);
    return;
  }

  if (this.hoveredId <= 0) {
    this.tooltip.style.display = "none";
    this._updateDistanceLine(null);
    return;
  }

  // Find the renderable with this pickId
  var allObjs = this.renderables.slice();
  if (this.sourceMarker) allObjs.push(this.sourceMarker);

  // Group by pickId (multiple parts of same object share pickId? No, they're unique)
  // But for subtraction-decomposed objects, multiple renderables share the same label/volume
  var found = null;
  for (var i = 0; i < allObjs.length; i++) {
    if (allObjs[i].pickId === this.hoveredId) {
      found = allObjs[i];
      break;
    }
  }

  if (!found) {
    this.tooltip.style.display = "none";
    this._updateDistanceLine(null);
    return;
  }

  // Build info for the hovered object
  var html = this._formatTooltipEntry(found, true);

  // Find contained/nested objects (objects whose center is inside the hovered object's bounding region)
  // Use a simple containment check: objects closer to detector axis with smaller extent
  var contained = this._findContainedObjects(found);
  if (contained.length > 0) {
    html += '<div style="border-top:1px solid rgba(255,255,255,0.15);margin-top:6px;padding-top:6px;font-size:12px;color:#aaa;">Contains:</div>';
    for (var ci = 0; ci < contained.length; ci++) {
      html += this._formatTooltipEntry(contained[ci], false);
    }
  }

  // Distance info for source volumes
  var distInfo = this._computeDistanceInfo(found);
  if (distInfo) {
    this._updateDistanceLine(distInfo);
    html += '<div style="border-top:1px solid rgba(255,255,255,0.15);margin-top:6px;padding-top:6px;">';
    html += '<div style="font-weight:bold;color:#ffee88;margin-bottom:2px;">Distances:</div>';
    for (var di = 0; di < distInfo.distances.length; di++) {
      var dd = distInfo.distances[di];
      html += '<div style="color:#ddd;">' + dd.label + ': ' + dd.value.toFixed(2) + ' cm</div>';
    }
    html += '</div>';
  } else {
    this._updateDistanceLine(null);
  }

  this.tooltip.innerHTML = html;
  this.tooltip.style.display = "block";

  // Position near mouse, clamped to container
  var rect = this.container.getBoundingClientRect();
  var tx = this.mousePos[0] + 16;
  var ty = this.mousePos[1] + 16;
  var tw = this.tooltip.offsetWidth;
  var th = this.tooltip.offsetHeight;
  if (tx + tw > rect.width - 10) tx = this.mousePos[0] - tw - 16;
  if (ty + th > rect.height - 10) ty = this.mousePos[1] - th - 16;
  this.tooltip.style.left = tx + "px";
  this.tooltip.style.top = ty + "px";
};


// ============================================================================
// Render
// ============================================================================

GeometryViewer3D.prototype._renderLoop = function() {
  var self = this;
  function loop() {
    self._animFrame = requestAnimationFrame(loop);
    self._doPick();
    if (!self.dirty) return;
    self.dirty = false;
    self._render();
  }
  loop();
};

GeometryViewer3D.prototype._render = function() {
  var gl = this.gl;
  var bg = this.options.backgroundColor;
  gl.clearColor(bg[0], bg[1], bg[2], bg[3]);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  // Enable uint32 indices
  var extUint = gl.getExtension("OES_element_index_uint");

  gl.enable(gl.DEPTH_TEST);
  gl.depthFunc(gl.LEQUAL);

  var proj = this._getProjectionMatrix();
  var view = this._getViewMatrix();
  var vp = mat4_multiply(proj, view);
  var eye = this._getEyePosition();

  // Light directions (in world space, normalized)
  var light1 = vec3_normalize([0.5, 0.7, 1.0]);
  var light2 = vec3_normalize([-0.6, -0.3, 0.5]);

  gl.useProgram(this.mainProgram);
  gl.uniform3fv(this.mainLocs.uLightDir1, light1);
  gl.uniform3fv(this.mainLocs.uLightDir2, light2);

  // Separate opaque and transparent
  var opaque = [], transparent = [];
  var allObjs = this.renderables.slice();
  if (this.sourceMarker) allObjs.push(this.sourceMarker);

  for (var i = 0; i < allObjs.length; i++) {
    var r = allObjs[i];
    if (r.alpha >= 0.99) {
      opaque.push(r);
    } else if (r.alpha > 0.01) {
      transparent.push(r);
    }
  }

  // Sort transparent by distance (back to front)
  transparent.sort(function(a, b) {
    var da = vec3_length(vec3_sub(a.center, eye));
    var db = vec3_length(vec3_sub(b.center, eye));
    return db - da;
  });

  // Render opaque
  gl.depthMask(true);
  gl.disable(gl.BLEND);
  this._renderObjects(opaque, vp, view, eye);

  // Render transparent
  gl.depthMask(false);
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  this._renderObjects(transparent, vp, view, eye);
  gl.depthMask(true);
  gl.disable(gl.BLEND);

  // Render axes
  if (this.options.showAxes) {
    this._renderAxes(vp);
  }

  // Render distance line (source to detector)
  if (this.distanceLine) {
    this._renderDistanceLine(vp);
  }
};

GeometryViewer3D.prototype._renderObjects = function(objects, vp, view, eye) {
  var gl = this.gl;
  var locs = this.mainLocs;

  for (var i = 0; i < objects.length; i++) {
    var r = objects[i];
    var mv = mat4_multiply(view, r.modelMatrix);
    var mvp = mat4_multiply(vp, r.modelMatrix);
    var normalMat = mat3_normalFromMat4(mv);

    gl.uniformMatrix4fv(locs.uMVP, false, mvp);
    gl.uniformMatrix4fv(locs.uMV, false, mv);
    gl.uniformMatrix3fv(locs.uNormalMat, false, normalMat);
    gl.uniform3fv(locs.uColor, r.color);
    gl.uniform1f(locs.uAlpha, r.alpha);
    gl.uniform1f(locs.uEmissive, r.pickId === this.hoveredId ? 0.18 : 0.0);

    gl.bindBuffer(gl.ARRAY_BUFFER, r.posBuffer);
    gl.enableVertexAttribArray(locs.aPosition);
    gl.vertexAttribPointer(locs.aPosition, 3, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, r.normBuffer);
    gl.enableVertexAttribArray(locs.aNormal);
    gl.vertexAttribPointer(locs.aNormal, 3, gl.FLOAT, false, 0, 0);

    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, r.idxBuffer);
    gl.drawElements(gl.TRIANGLES, r.numIndices, gl.UNSIGNED_INT, 0);
  }
};

GeometryViewer3D.prototype._renderAxes = function(vp) {
  var gl = this.gl;
  gl.useProgram(this.lineProgram);
  gl.uniformMatrix4fv(this.lineLocs.uMVP, false, vp);

  gl.bindBuffer(gl.ARRAY_BUFFER, this.axesBuffers.posBuffer);
  gl.enableVertexAttribArray(this.lineLocs.aPosition);
  gl.vertexAttribPointer(this.lineLocs.aPosition, 3, gl.FLOAT, false, 0, 0);

  gl.bindBuffer(gl.ARRAY_BUFFER, this.axesBuffers.colBuffer);
  gl.enableVertexAttribArray(this.lineLocs.aColor);
  gl.vertexAttribPointer(this.lineLocs.aColor, 3, gl.FLOAT, false, 0, 0);

  gl.lineWidth(2.0);
  gl.drawArrays(gl.LINES, 0, this.axesBuffers.count);
};


// End of IIFE
})();
