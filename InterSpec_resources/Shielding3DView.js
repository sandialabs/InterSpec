/* Shielding3DView.js - WebGL version */

Shielding3DView = function(id, data) {
  var container = document.getElementById(id);
  if (!container) return null;
  
  var self = this;
  this.container = container;
  this.data = null;
  this.meshes = [];
  this.gl = null;
  this.canvas = null;
  this.program = null;
  this.locs = null;

  // --- Constants & Units ---
  // data.distance and shielding dimensions are in base units (mm, since cm=10.0).
  // CM = 10.0.
  var CM = 10.0; 
  var INCH = 2.54 * CM;

  // --- WebGL Setup ---
  this.canvas = document.createElement('canvas');
  this.canvas.style.width = '100%';
  this.canvas.style.height = '100%';
  // Initial size
  this.canvas.width = container.clientWidth;
  this.canvas.height = container.clientHeight;
  container.innerHTML = '';
  container.appendChild(this.canvas);

  // ResizeObserver
  var resizeObserver = new ResizeObserver(function(entries) {
    for (var i = 0; i < entries.length; i++) {
      var entry = entries[i];
      if (entry.target === container) {
        var width = entry.contentRect.width;
        var height = entry.contentRect.height;
        
        // Update canvas size
        self.canvas.width = width;
        self.canvas.height = height;
        
        // Optionally force a render here if loop isn't running
      }
    }
  });
  resizeObserver.observe(container);

  this.gl = this.canvas.getContext('webgl') || this.canvas.getContext('experimental-webgl');
  if (!this.gl) {
    container.innerHTML = 'WebGL not supported';
    return null;
  }
  
  var gl = this.gl;
  var canvas = this.canvas;

  // --- Shaders ---
  var vertexShaderSrc = [
    'attribute vec3 aPosition;',
    'attribute vec3 aNormal;',
    'uniform mat4 uModelViewProjection;',
    'uniform mat4 uModel;',
    'uniform mat3 uNormalMatrix;',
    'varying vec3 vNormal;',
    'varying vec3 vPos;',
    'void main() {',
    '  vNormal = normalize(uNormalMatrix * aNormal);',
    '  vec4 pos = uModel * vec4(aPosition, 1.0);',
    '  vPos = pos.xyz;',
    '  gl_Position = uModelViewProjection * vec4(aPosition, 1.0);',
    '}'
  ].join('\n');

  var fragmentShaderSrc = [
    'precision mediump float;',
    'uniform vec4 uColor;',
    'uniform vec3 uLightDir;',
    'varying vec3 vNormal;',
    'varying vec3 vPos;',
    'void main() {',
    '  vec3 normal = normalize(vNormal);',
    '  float diff = max(dot(normal, uLightDir), 0.3);',
    '  gl_FragColor = vec4(uColor.rgb * diff, uColor.a);',
    '}'
  ].join('\n');

  function createShader(gl, type, source) {
    var shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      console.error('Shader error:', gl.getShaderInfoLog(shader));
      gl.deleteShader(shader);
      return null;
    }
    return shader;
  }

  var vertexShader = createShader(gl, gl.VERTEX_SHADER, vertexShaderSrc);
  var fragmentShader = createShader(gl, gl.FRAGMENT_SHADER, fragmentShaderSrc);
  var program = gl.createProgram();
  gl.attachShader(program, vertexShader);
  gl.attachShader(program, fragmentShader);
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    console.error('Program error:', gl.getProgramInfoLog(program));
    return;
  }

  this.program = program;
  this.locs = {
    aPosition: gl.getAttribLocation(program, 'aPosition'),
    aNormal: gl.getAttribLocation(program, 'aNormal'),
    uModelViewProjection: gl.getUniformLocation(program, 'uModelViewProjection'),
    uModel: gl.getUniformLocation(program, 'uModel'),
    uNormalMatrix: gl.getUniformLocation(program, 'uNormalMatrix'),
    uColor: gl.getUniformLocation(program, 'uColor'),
    uLightDir: gl.getUniformLocation(program, 'uLightDir')
  };
  
  var locs = this.locs;

  // --- Math Helpers ---
  var Mat4 = {
    create: function() { return new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]); },
    perspective: function(out, fovy, aspect, near, far) {
      var f = 1.0 / Math.tan(fovy / 2), nf = 1 / (near - far);
      out[0] = f / aspect; out[1] = 0; out[2] = 0; out[3] = 0;
      out[4] = 0; out[5] = f; out[6] = 0; out[7] = 0;
      out[8] = 0; out[9] = 0; out[10] = (far + near) * nf; out[11] = -1;
      out[12] = 0; out[13] = 0; out[14] = (2 * far * near) * nf; out[15] = 0;
      return out;
    },
    lookAt: function(out, eye, center, up) {
      var x0, x1, x2, y0, y1, y2, z0, z1, z2, len;
      var eyex = eye[0], eyey = eye[1], eyez = eye[2];
      var upx = up[0], upy = up[1], upz = up[2];
      var centerx = center[0], centery = center[1], centerz = center[2];
      z0 = eyex - centerx; z1 = eyey - centery; z2 = eyez - centerz;
      len = 1 / Math.sqrt(z0 * z0 + z1 * z1 + z2 * z2);
      z0 *= len; z1 *= len; z2 *= len;
      x0 = upy * z2 - upz * z1; x1 = upz * z0 - upx * z2; x2 = upx * z1 - upy * z0;
      len = Math.sqrt(x0 * x0 + x1 * x1 + x2 * x2);
      if (!len) { x0 = 0; x1 = 0; x2 = 0; } else { len = 1 / len; x0 *= len; x1 *= len; x2 *= len; }
      y0 = z1 * x2 - z2 * x1; y1 = z2 * x0 - z0 * x2; y2 = z0 * x1 - z1 * x0;
      len = Math.sqrt(y0 * y0 + y1 * y1 + y2 * y2);
      if (!len) { y0 = 0; y1 = 0; y2 = 0; } else { len = 1 / len; y0 *= len; y1 *= len; y2 *= len; }
      out[0] = x0; out[1] = y0; out[2] = z0; out[3] = 0;
      out[4] = x1; out[5] = y1; out[6] = z1; out[7] = 0;
      out[8] = x2; out[9] = y2; out[10] = z2; out[11] = 0;
      out[12] = -(x0 * eyex + x1 * eyey + x2 * eyez);
      out[13] = -(y0 * eyex + y1 * eyey + y2 * eyez);
      out[14] = -(z0 * eyex + z1 * eyey + z2 * eyez);
      out[15] = 1;
      return out;
    },
    multiply: function(out, a, b) {
      var a00=a[0], a01=a[1], a02=a[2], a03=a[3];
      var a10=a[4], a11=a[5], a12=a[6], a13=a[7];
      var a20=a[8], a21=a[9], a22=a[10], a23=a[11];
      var a30=a[12], a31=a[13], a32=a[14], a33=a[15];
      var b0=b[0], b1=b[1], b2=b[2], b3=b[3];
      out[0]=b0*a00+b1*a10+b2*a20+b3*a30; out[1]=b0*a01+b1*a11+b2*a21+b3*a31;
      out[2]=b0*a02+b1*a12+b2*a22+b3*a32; out[3]=b0*a03+b1*a13+b2*a23+b3*a33;
      b0=b[4]; b1=b[5]; b2=b[6]; b3=b[7];
      out[4]=b0*a00+b1*a10+b2*a20+b3*a30; out[5]=b0*a01+b1*a11+b2*a21+b3*a31;
      out[6]=b0*a02+b1*a12+b2*a22+b3*a32; out[7]=b0*a03+b1*a13+b2*a23+b3*a33;
      b0=b[8]; b1=b[9]; b2=b[10]; b3=b[11];
      out[8]=b0*a00+b1*a10+b2*a20+b3*a30; out[9]=b0*a01+b1*a11+b2*a21+b3*a31;
      out[10]=b0*a02+b1*a12+b2*a22+b3*a32; out[11]=b0*a03+b1*a13+b2*a23+b3*a33;
      b0=b[12]; b1=b[13]; b2=b[14]; b3=b[15];
      out[12]=b0*a00+b1*a10+b2*a20+b3*a30; out[13]=b0*a01+b1*a11+b2*a21+b3*a31;
      out[14]=b0*a02+b1*a12+b2*a22+b3*a32; out[15]=b0*a03+b1*a13+b2*a23+b3*a33;
      return out;
    },
    translate: function(out, a, v) {
      var x=v[0], y=v[1], z=v[2];
      if (a !== out) {
        out[0]=a[0]; out[1]=a[1]; out[2]=a[2]; out[3]=a[3];
        out[4]=a[4]; out[5]=a[5]; out[6]=a[6]; out[7]=a[7];
        out[8]=a[8]; out[9]=a[9]; out[10]=a[10]; out[11]=a[11];
      }
      out[12]=a[0]*x+a[4]*y+a[8]*z+a[12];
      out[13]=a[1]*x+a[5]*y+a[9]*z+a[13];
      out[14]=a[2]*x+a[6]*y+a[10]*z+a[14];
      out[15]=a[3]*x+a[7]*y+a[11]*z+a[15];
      return out;
    },
    rotate: function(out, a, rad, axis) {
      var x=axis[0], y=axis[1], z=axis[2];
      var len = Math.sqrt(x*x + y*y + z*z);
      if (len < 0.000001) return null;
      len = 1/len; x*=len; y*=len; z*=len;
      var s = Math.sin(rad), c = Math.cos(rad), t = 1-c;
      var a00=a[0], a01=a[1], a02=a[2], a03=a[3];
      var a10=a[4], a11=a[5], a12=a[6], a13=a[7];
      var a20=a[8], a21=a[9], a22=a[10], a23=a[11];
      var b00=x*x*t+c, b01=y*x*t+z*s, b02=z*x*t-y*s;
      var b10=x*y*t-z*s, b11=y*y*t+c, b12=z*y*t+x*s;
      var b20=x*z*t+y*s, b21=y*z*t-x*s, b22=z*z*t+c;
      out[0]=a00*b00+a10*b01+a20*b02; out[1]=a01*b00+a11*b01+a21*b02; out[2]=a02*b00+a12*b01+a22*b02; out[3]=a03*b00+a13*b01+a23*b02;
      out[4]=a00*b10+a10*b11+a20*b12; out[5]=a01*b10+a11*b11+a21*b12; out[6]=a02*b10+a12*b11+a22*b12; out[7]=a03*b10+a13*b11+a23*b12;
      out[8]=a00*b20+a10*b21+a20*b22; out[9]=a01*b20+a11*b21+a21*b22; out[10]=a02*b20+a12*b21+a22*b22; out[11]=a03*b20+a13*b21+a23*b22;
      if(a!==out) {out[12]=a[12]; out[13]=a[13]; out[14]=a[14]; out[15]=a[15];}
      return out;
    },
    scale: function(out, a, v) {
      var x=v[0], y=v[1], z=v[2];
      out[0]=a[0]*x; out[1]=a[1]*x; out[2]=a[2]*x; out[3]=a[3]*x;
      out[4]=a[4]*y; out[5]=a[5]*y; out[6]=a[6]*y; out[7]=a[7]*y;
      out[8]=a[8]*z; out[9]=a[9]*z; out[10]=a[10]*z; out[11]=a[11]*z;
      out[12]=a[12]; out[13]=a[13]; out[14]=a[14]; out[15]=a[15];
      return out;
    },
    normalFromMat4: function(out, a) {
      // Invert and transpose upper 3x3
      var a00=a[0], a01=a[1], a02=a[2], a03=a[3];
      var a10=a[4], a11=a[5], a12=a[6], a13=a[7];
      var a20=a[8], a21=a[9], a22=a[10], a23=a[11];
      var a30=a[12], a31=a[13], a32=a[14], a33=a[15];
      var b00=a00*a11-a01*a10, b01=a00*a12-a02*a10, b02=a00*a13-a03*a10;
      var b03=a01*a12-a02*a11, b04=a01*a13-a03*a11, b05=a02*a13-a03*a12;
      var b06=a20*a31-a21*a30, b07=a20*a32-a22*a30, b08=a20*a33-a23*a30;
      var b09=a21*a32-a22*a31, b10=a21*a33-a23*a31, b11=a22*a33-a23*a32;
      var det = b00*b11 - b01*b10 + b02*b09 + b03*b08 - b04*b07 + b05*b06;
      if (!det) return null;
      det = 1.0 / det;
      out[0] = (a11*b11 - a12*b10 + a13*b09) * det;
      out[1] = (a12*b08 - a10*b11 - a13*b07) * det;
      out[2] = (a10*b10 - a11*b08 + a13*b06) * det;
      out[3] = (a02*b10 - a01*b11 - a03*b09) * det;
      out[4] = (a00*b11 - a02*b08 + a03*b07) * det;
      out[5] = (a01*b08 - a00*b10 - a03*b06) * det;
      out[6] = (a31*b05 - a32*b04 + a33*b03) * det;
      out[7] = (a32*b02 - a30*b05 - a33*b01) * det;
      out[8] = (a30*b04 - a31*b02 + a33*b00) * det;
      return out;
    }
  };

  // --- Geometry Generators ---
  function createBufferInfo(gl, arrays) {
    var bufferInfo = {
      numElements: arrays.indices.length,
      indices: gl.createBuffer(),
      position: gl.createBuffer(),
      normal: gl.createBuffer()
    };
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bufferInfo.indices);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(arrays.indices), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, bufferInfo.position);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(arrays.position), gl.STATIC_DRAW);
    gl.bindBuffer(gl.ARRAY_BUFFER, bufferInfo.normal);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(arrays.normal), gl.STATIC_DRAW);
    return bufferInfo;
  }

  function createSphere(radius, latBands, longBands) {
    var positions = [];
    var normals = [];
    var indices = [];
    for (var lat = 0; lat <= latBands; lat++) {
      var theta = lat * Math.PI / latBands;
      var sinTheta = Math.sin(theta);
      var cosTheta = Math.cos(theta);
      for (var lon = 0; lon <= longBands; lon++) {
        var phi = lon * 2 * Math.PI / longBands;
        var sinPhi = Math.sin(phi);
        var cosPhi = Math.cos(phi);
        var x = cosPhi * sinTheta;
        var y = cosTheta;
        var z = sinPhi * sinTheta;
        normals.push(x, y, z);
        positions.push(radius * x, radius * y, radius * z);
      }
    }
    for (var lat = 0; lat < latBands; lat++) {
      for (var lon = 0; lon < longBands; lon++) {
        var first = (lat * (longBands + 1)) + lon;
        var second = first + longBands + 1;
        indices.push(first, second, first + 1);
        indices.push(second, second + 1, first + 1);
      }
    }
    return { position: positions, normal: normals, indices: indices };
  }

  function createCylinder(radius, length, radialSegments) {
    var positions = [];
    var normals = [];
    var indices = [];
    // Side
    for (var i = 0; i <= radialSegments; i++) {
      var theta = i * 2 * Math.PI / radialSegments;
      var x = Math.sin(theta);
      var y = Math.cos(theta);
      // Top edge
      normals.push(x, y, 0); positions.push(radius * x, radius * y, length / 2);
      // Bottom edge
      normals.push(x, y, 0); positions.push(radius * x, radius * y, -length / 2);
    }
    for (var i = 0; i < radialSegments; i++) {
        var idx = i * 2;
        indices.push(idx, idx + 1, idx + 2);
        indices.push(idx + 1, idx + 3, idx + 2);
    }
    // Caps (simple triangle fan approximation with center point)
    var topCenterIdx = positions.length / 3;
    positions.push(0, 0, length / 2); normals.push(0, 0, 1);
    var bottomCenterIdx = topCenterIdx + 1;
    positions.push(0, 0, -length / 2); normals.push(0, 0, -1);
    
    // Add cap vertices with correct normals
    var capStartIdx = bottomCenterIdx + 1;
    for (var i = 0; i <= radialSegments; i++) {
        var theta = i * 2 * Math.PI / radialSegments;
        var x = Math.sin(theta) * radius;
        var y = Math.cos(theta) * radius;
        positions.push(x, y, length/2); normals.push(0, 0, 1);
        positions.push(x, y, -length/2); normals.push(0, 0, -1);
    }
    for (var i = 0; i < radialSegments; i++) {
        var top1 = capStartIdx + i * 2;
        var top2 = top1 + 2;
        indices.push(topCenterIdx, top2, top1); // Top cap
        var bot1 = top1 + 1;
        var bot2 = top2 + 1;
        indices.push(bottomCenterIdx, bot1, bot2); // Bottom cap
    }
    return { position: positions, normal: normals, indices: indices };
  }

  function createBox(width, height, depth) {
    // 8 vertices, but normals differ per face, so 24 vertices
    var w = width/2, h = height/2, d = depth/2;
    var positions = [
      // Front
      -w, -h,  d,   w, -h,  d,   w,  h,  d,  -w,  h,  d,
      // Back
      -w, -h, -d,  -w,  h, -d,   w,  h, -d,   w, -h, -d,
      // Top
      -w,  h, -d,  -w,  h,  d,   w,  h,  d,   w,  h, -d,
      // Bottom
      -w, -h, -d,   w, -h, -d,   w, -h,  d,  -w, -h,  d,
      // Right
       w, -h, -d,   w,  h, -d,   w,  h,  d,   w, -h,  d,
      // Left
      -w, -h, -d,  -w, -h,  d,  -w,  h,  d,  -w,  h, -d
    ];
    var normals = [
      // Front
      0, 0, 1,  0, 0, 1,  0, 0, 1,  0, 0, 1,
      // Back
      0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1,
      // Top
      0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,
      // Bottom
      0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0,
      // Right
      1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,
      // Left
      -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0
    ];
    var indices = [
      0, 1, 2,      0, 2, 3,    // Front
      4, 5, 6,      4, 6, 7,    // Back
      8, 9, 10,     8, 10, 11,  // Top
      12, 13, 14,   12, 14, 15, // Bottom
      16, 17, 18,   16, 18, 19, // Right
      20, 21, 22,   20, 22, 23  // Left
    ];
    return { position: positions, normal: normals, indices: indices };
  }

  // --- Scene Objects Construction ---
  this.meshes = [];
  this.camRadius = 0;
  this.camTheta = Math.PI / 4;
  this.camPhi = Math.PI / 4;
  this.camTarget = [0, 0, 0];
  this.isDragging = false;
  this.lastMouseX = 0;
  this.lastMouseY = 0;
  this.button = -1;
  this.renderLoopId = null;
  
  // Initialize with data if provided
  if (data) {
    this.setData(data);
  }
  
  return this;
};

// Set data and rebuild meshes
Shielding3DView.prototype.setData = function(data) {
  this.data = data;
  
  // Clear existing meshes (but keep WebGL context)
  if (this.meshes) {
    // Clean up old buffers
    for (var i = 0; i < this.meshes.length; i++) {
      if (this.meshes[i].buffer) {
        if (this.meshes[i].buffer.position) this.gl.deleteBuffer(this.meshes[i].buffer.position);
        if (this.meshes[i].buffer.normal) this.gl.deleteBuffer(this.meshes[i].buffer.normal);
        if (this.meshes[i].buffer.indices) this.gl.deleteBuffer(this.meshes[i].buffer.indices);
      }
    }
  }
  
  this.meshes = [];
  
  // Remove old overlay if it exists
  var oldOverlay = this.container.querySelector('div[style*="position: absolute"]');
  if (oldOverlay && oldOverlay !== this.canvas) {
    oldOverlay.remove();
  }
  
  if (!data) return;
  
  var gl = this.gl;
  var container = this.container;
  var self = this;
  
  // Helper functions (same as in constructor)
  function createBufferInfo(gl, meshData) {
    var posBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(meshData.position), gl.STATIC_DRAW);
    
    var normBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, normBuf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(meshData.normal), gl.STATIC_DRAW);
    
    var idxBuf = gl.createBuffer();
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, idxBuf);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(meshData.indices), gl.STATIC_DRAW);
    
    return {
      position: posBuf,
      normal: normBuf,
      indices: idxBuf,
      numElements: meshData.indices.length
    };
  }
  
  function createSphere(radius, slices, stacks) {
    var vertices = [];
    var normals = [];
    var indices = [];
    for (var i = 0; i <= stacks; i++) {
      var theta = i * Math.PI / stacks;
      var sinTheta = Math.sin(theta);
      var cosTheta = Math.cos(theta);
      for (var j = 0; j <= slices; j++) {
        var phi = j * 2 * Math.PI / slices;
        var sinPhi = Math.sin(phi);
        var cosPhi = Math.cos(phi);
        var x = cosPhi * sinTheta;
        var y = cosTheta;
        var z = sinPhi * sinTheta;
        vertices.push(radius * x, radius * y, radius * z);
        normals.push(x, y, z);
      }
    }
    for (var i = 0; i < stacks; i++) {
      for (var j = 0; j < slices; j++) {
        var first = i * (slices + 1) + j;
        var second = first + slices + 1;
        indices.push(first, second, first + 1);
        indices.push(second, second + 1, first + 1);
      }
    }
    return { position: vertices, normal: normals, indices: indices };
  }
  
  function createCylinder(radius, height, segments) {
    var vertices = [];
    var normals = [];
    var indices = [];
    for (var i = 0; i <= segments; i++) {
      var angle = 2 * Math.PI * i / segments;
      var x = Math.cos(angle);
      var z = Math.sin(angle);
      vertices.push(radius * x, -height/2, radius * z);
      normals.push(x, 0, z);
      vertices.push(radius * x, height/2, radius * z);
      normals.push(x, 0, z);
    }
    for (var i = 0; i < segments; i++) {
      var base = i * 2;
      indices.push(base, base + 1, (base + 2) % (segments * 2 + 2));
      indices.push(base + 1, (base + 3) % (segments * 2 + 2), (base + 2) % (segments * 2 + 2));
    }
    // Top and bottom caps
    var centerBottom = segments * 2 + 2;
    var centerTop = centerBottom + 1;
    vertices.push(0, -height/2, 0);
    normals.push(0, -1, 0);
    vertices.push(0, height/2, 0);
    normals.push(0, 1, 0);
    for (var i = 0; i < segments; i++) {
      var base = i * 2;
      var next = ((i + 1) % segments) * 2;
      indices.push(centerBottom, next, base);
      indices.push(centerTop, base + 1, next + 1);
    }
    return { position: vertices, normal: normals, indices: indices };
  }
  
  function createBox(w, h, d) {
    var positions = [
      // Front
       w,  h,  d,   w, -h,  d,  -w, -h,  d,  -w,  h,  d,
      // Back
      -w,  h, -d,  -w, -h, -d,   w, -h, -d,   w,  h, -d,
      // Top
      -w,  h,  d,   w,  h,  d,   w,  h, -d,  -w,  h, -d,
      // Bottom
      -w, -h, -d,   w, -h, -d,   w, -h,  d,  -w, -h,  d,
      // Right
       w,  h, -d,   w,  h,  d,   w, -h,  d,   w, -h, -d,
      // Left
      -w, -h, -d,  -w, -h,  d,  -w,  h,  d,  -w,  h, -d
    ];
    var normals = [
      0, 0, 1,  0, 0, 1,  0, 0, 1,  0, 0, 1,
      0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1,
      0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,
      0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0,
      1, 0, 0,  1, 0, 0,  1, 0, 0,  1, 0, 0,
      -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0
    ];
    var indices = [
      0, 1, 2,      0, 2, 3,    // Front
      4, 5, 6,      4, 6, 7,    // Back
      8, 9, 10,     8, 10, 11,  // Top
      12, 13, 14,   12, 14, 15, // Bottom
      16, 17, 18,   16, 18, 19, // Right
      20, 21, 22,   20, 22, 23  // Left
    ];
    return { position: positions, normal: normals, indices: indices };
  }
  
  // Now create meshes (same logic as constructor)
  var meshes = this.meshes;
  
  // 1. Detector
  var detDistance = data.distance * CM; // Center of shielding to Face of Detector (data.distance is in mm, convert to raw units?)
  // Wait, data.distance in JSON is in mm (createJsonData divides by mm).
  // In JS, we want to render.
  // Let's stick to a unit scale.
  // If we assume 1 unit = 1 mm.
  // CM = 10.0 is weird if we want 1 unit = 1 mm.
  // Let's redefine units to be consistent with JSON input (mm).
  var MM = 1.0;
  var CM = 10.0 * MM;
  var INCH = 25.4 * MM;
  
  var detDistance = data.distance * MM;
  
  // Detector: Cylinder
  // Radius = diameter / 2
  var detDiam = (data.detectorDiameter !== undefined) ? (data.detectorDiameter * MM) : (3.0 * INCH);
  var detRad = detDiam / 2.0;
  
  // Length defaults to 3 inches or maybe should scale?
  // Let's keep length as 3 inches for now or equal to diameter if we want a 1:1 aspect ratio like typical NaI?
  // Actually typical detectors are often 1:1.
  var detLen = (data.detectorDiameter !== undefined) ? (data.detectorDiameter * MM) : (3.0 * INCH);
  
  // Position: Face at Z = detDistance. Cylinder center at Z = detDistance + detLen/2.
  var detZ = detDistance + detLen/2;
  
  // Add Cylinder mesh (Unit cylinder scaled)
  // We use generated buffers for specific sizes to avoid complex scaling normals issues if non-uniform
  meshes.push({
    buffer: createBufferInfo(gl, createCylinder(detRad, detLen, 32)),
    color: [0.7, 0.7, 0.7, 1.0], // Grey
    matrix: Mat4.translate(Mat4.create(), Mat4.create(), [0, 0, detZ]),
    transparent: false
  });

  // 2. Detector Body (Box 6")
  var boxDim = 6.0 * INCH;
  // Position: Attached to back of cylinder?
  // Back of cylinder is at Z = detDistance + detLen.
  // Box center at Z = detDistance + detLen + boxDim/2.
  // Box centered in X/Y?
  var boxZ = detDistance + detLen + boxDim/2;
  meshes.push({
    buffer: createBufferInfo(gl, createBox(boxDim, boxDim, boxDim)),
    color: [0.2, 0.2, 0.8, 1.0], // Blueish
    matrix: Mat4.translate(Mat4.create(), Mat4.create(), [0, 0, boxZ]),
    transparent: false
  });

  // 3. Handle (U shape on top)
  // Simple box handle
  var handW = 4 * INCH, handH = 1 * INCH, handD = 1 * INCH;
  var handZ = boxZ;
  var handY = boxDim/2 + handH/2;
  meshes.push({
    buffer: createBufferInfo(gl, createBox(handW, handH, handD)),
    color: [0.1, 0.1, 0.1, 1.0], // Dark
    matrix: Mat4.translate(Mat4.create(), Mat4.create(), [0, handY, handZ]),
    transparent: false
  });

  // 4. Shielding
  // Iterate data.shieldings.
  // Accumulate dimensions.
  var currentRad = 0;
  var currentW = 0, currentH = 0, currentD = 0;
  // For cylinders
  var currentCylRad = 0;
  var currentCylLen = 0;

  // Shielding is at 0,0,0.
  // Need to draw Inner to Outer.
  // data.shieldings is ordered inner to outer? Usually.
  
  for (var i = 0; i < data.shieldings.length; i++) {
    var s = data.shieldings[i];
    var dims = s.dimensions; // [d0, d1, d2]
    var geo = data.geometry; 
    
    var meshData = null;
    var color = [Math.random()*0.5+0.5, Math.random()*0.5+0.5, Math.random()*0.5+0.5, 0.4]; // Random transparent color
    if(s.material.indexOf("Lead") >= 0) color = [0.3, 0.3, 0.3, 0.5];
    if(s.material.indexOf("Iron") >= 0) color = [0.6, 0.3, 0.2, 0.5];
    if(s.material.indexOf("Copper") >= 0) color = [0.8, 0.5, 0.3, 0.5];
    if(s.material.indexOf("Aluminum") >= 0) color = [0.8, 0.8, 0.9, 0.5];
    
    if (geo === "Spherical") {
      currentRad += dims[0];
      meshData = createSphere(currentRad, 32, 32);
    } else if (geo === "CylinderEndOn" || geo === "CylinderSideOn") {
      // dims[0] = Radius Thickness, dims[1] = Length Thickness?
      // Assuming they are additive thicknesses based on typical code patterns
      currentCylRad += dims[0];
      currentCylLen += dims[1]; 
      meshData = createCylinder(currentCylRad, currentCylLen, 32);
      // Orientation: Cylinder along Z.
      // EndOn vs SideOn usually implies orientation relative to detector.
      // If EndOn, cylinder axis points to detector. Detector is at +Z. So Cylinder along Z is correct.
      // If SideOn, cylinder axis perpendicular to detector line. So Cylinder along Y or X.
      // Let's rotate 90 deg X if SideOn.
    } else if (geo === "Rectangular") {
      currentW += dims[0];
      currentH += dims[1];
      currentD += dims[2];
      meshData = createBox(currentW, currentH, currentD);
    }

    if (meshData) {
      var mat = Mat4.create();
      if (geo === "CylinderSideOn") {
         Mat4.rotate(mat, mat, Math.PI/2, [1, 0, 0]);
      }
      
      meshes.push({
        buffer: createBufferInfo(gl, meshData),
        color: color,
        matrix: mat,
        transparent: true
      });
    }
  }

  // 5. Point Sources
  // Check if there are point sources in data.fitSources
  var hasPointSource = false;
  if (data.fitSources) {
      for (var i = 0; i < data.fitSources.length; i++) {
          if (data.fitSources[i].type === "Point") {
              hasPointSource = true;
              break;
          }
      }
  }

  if (hasPointSource) {
      var srcRad = 2 * MM; // Small dot 2mm radius
      meshes.push({
        buffer: createBufferInfo(gl, createSphere(srcRad, 16, 16)),
        color: [1.0, 0.0, 0.0, 1.0], // Red
        matrix: Mat4.translate(Mat4.create(), Mat4.create(), [0, 0, 0]), // At origin
        transparent: false
      });
  }

  // --- Source Overlay ---
  if (data.fitSources && data.fitSources.length > 0) {
      var overlay = document.createElement('div');
      overlay.style.position = 'absolute';
      overlay.style.top = '10px';
      overlay.style.left = '10px';
      overlay.style.backgroundColor = 'rgba(255, 255, 255, 0.8)';
      overlay.style.padding = '10px';
      overlay.style.borderRadius = '5px';
      overlay.style.fontFamily = 'sans-serif';
      overlay.style.fontSize = '12px';
      overlay.style.pointerEvents = 'none'; // Let clicks pass through
      
      var html = "<strong>Sources:</strong><br/>";
      for (var i = 0; i < data.fitSources.length; i++) {
          var s = data.fitSources[i];
          html += s.nuclide + " (" + s.type + "): " + s.activityPretty;
          // if (s.activityUncert) html += " \u00B1 ...";
          html += "<br/>";
      }
      overlay.innerHTML = html;
      container.appendChild(overlay);
  }
  
  // Update camera target based on new data
  var detDistance = data.distance * MM;
  this.camTarget = [0, 0, detDistance/2];
  this.camRadius = detDistance * 1.5 + 10 * INCH;
  if (this.camRadius < 100) this.camRadius = 100;
};

  // --- Camera & Interaction ---
  // Initial zoom to show shielding and detector.
  // Extents: Shielding around 0. Detector around detDistance.
  // Center of view: (0 + detDistance)/2?
  // Radius of scene: detDistance + size.
  var camRadius = detDistance * 1.5 + 10 * INCH;
  if (camRadius < 100) camRadius = 100;
  
  var camTheta = Math.PI / 4; // 45 deg
  var camPhi = Math.PI / 4;
  var camTarget = [0, 0, detDistance/2];
  
  var isDragging = false;
  var lastMouseX = 0;
  var lastMouseY = 0;
  var button = -1;

  canvas.addEventListener('mousedown', function(e) {
    isDragging = true;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;
    button = e.button;
    e.preventDefault(); // Prevent text selection
  });

  window.addEventListener('mouseup', function() {
    isDragging = false;
  });

  canvas.addEventListener('mousemove', function(e) {
    if (!isDragging) return;
    var dx = e.clientX - lastMouseX;
    var dy = e.clientY - lastMouseY;
    lastMouseX = e.clientX;
    lastMouseY = e.clientY;

    if (button === 0) { // Left Button: Rotate
      camPhi -= dx * 0.01;
      camTheta -= dy * 0.01;
      // Clamp theta to avoid gimbal lock or flipping
      if (camTheta < 0.01) camTheta = 0.01;
      if (camTheta > Math.PI - 0.01) camTheta = Math.PI - 0.01;
    } else if (button === 2) { // Right Button: Pan
      // Move camTarget relative to camera view
      // Simple implementation: Pan in screen plane
      // We need Up and Right vectors of camera.
      // For now, just crude pan along world axes adjusted by phi
      var panSpeed = camRadius * 0.002; // Scale pan with distance
      // This is non-trivial without full view matrix basis access easily available here
      // Let's just modify camTarget Z/X/Y
      // Better: Just standard orbit pan logic roughly
      camTarget[0] -= dx * panSpeed * Math.cos(camPhi);
      camTarget[2] -= dx * panSpeed * Math.sin(camPhi);
      camTarget[1] += dy * panSpeed;
    }
  });
  
  canvas.addEventListener('contextmenu', function(e) { e.preventDefault(); return false; }); // Disable context menu

  canvas.addEventListener('wheel', function(e) {
    e.preventDefault();
    var delta = Math.sign(e.deltaY) * (camRadius * 0.1);
    camRadius += delta;
    if (camRadius < 1.0) camRadius = 1.0;
  });

  // --- Render Loop ---
  function render() {
    // Resize
    if (canvas.width !== canvas.clientWidth || canvas.height !== canvas.clientHeight) {
        canvas.width = canvas.clientWidth;
        canvas.height = canvas.clientHeight;
        gl.viewport(0, 0, canvas.width, canvas.height);
    }

    gl.clearColor(0.9, 0.9, 0.9, 1.0);
    gl.enable(gl.DEPTH_TEST);
    gl.enable(gl.CULL_FACE);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.enable(gl.BLEND);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    gl.useProgram(program);

    // Camera Matrix
    var eyeY = camTarget[1] + camRadius * Math.cos(camTheta);
    var horizRad = camRadius * Math.sin(camTheta);
    var eyeX = camTarget[0] + horizRad * Math.sin(camPhi);
    var eyeZ = camTarget[2] + horizRad * Math.cos(camPhi);
    
    var view = Mat4.lookAt(Mat4.create(), [eyeX, eyeY, eyeZ], camTarget, [0, 1, 0]);
    var projection = Mat4.perspective(Mat4.create(), Math.PI / 4, canvas.width / canvas.height, 1.0, 10000.0); // near 1, far 10000
    var viewProj = Mat4.multiply(Mat4.create(), projection, view);

    var lightDir = [1, 1, 1]; // normalized in shader usually
    var norm = Math.sqrt(3);
    gl.uniform3f(locs.uLightDir, 1/norm, 1/norm, 1/norm);

    // Draw Opaques First
    // Actually, we have a list. Let's separate them.
    var opaque = meshes.filter(function(m) { return !m.transparent; });
    var trans = meshes.filter(function(m) { return m.transparent; });
    
    // Sort transparent back-to-front?
    // Distance from camera to object center.
    trans.sort(function(a, b) {
       // Get position from matrix (col 3)
       var ax = a.matrix[12], ay = a.matrix[13], az = a.matrix[14];
       var bx = b.matrix[12], by = b.matrix[13], bz = b.matrix[14];
       var adist = (ax-eyeX)*(ax-eyeX) + (ay-eyeY)*(ay-eyeY) + (az-eyeZ)*(az-eyeZ);
       var bdist = (bx-eyeX)*(bx-eyeX) + (by-eyeY)*(by-eyeY) + (bz-eyeZ)*(bz-eyeZ);
       return bdist - adist; // Descending
    });

    var drawList = opaque.concat(trans);

    for (var i = 0; i < drawList.length; i++) {
      var mesh = drawList[i];
      var mvp = Mat4.multiply(Mat4.create(), viewProj, mesh.matrix);
      
      var normMat = new Float32Array(9);
      // Normal matrix is inverse transpose of model upper 3x3.
      // Since we only translate/rotate (and maybe uniform scale), model matrix upper 3x3 works if orthogonal.
      // Mat4.normalFromMat4 handles inverse transpose.
      Mat4.normalFromMat4(normMat, mesh.matrix);

      gl.uniformMatrix4fv(locs.uModelViewProjection, false, mvp);
      gl.uniformMatrix4fv(locs.uModel, false, mesh.matrix);
      gl.uniformMatrix3fv(locs.uNormalMatrix, false, normMat);
      gl.uniform4fv(locs.uColor, mesh.color);

      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, mesh.buffer.indices);
      gl.bindBuffer(gl.ARRAY_BUFFER, mesh.buffer.position);
      gl.vertexAttribPointer(locs.aPosition, 3, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(locs.aPosition);
      
      gl.bindBuffer(gl.ARRAY_BUFFER, mesh.buffer.normal);
      gl.vertexAttribPointer(locs.aNormal, 3, gl.FLOAT, false, 0, 0);
      gl.enableVertexAttribArray(locs.aNormal);

      if (mesh.transparent) {
          gl.depthMask(false); // Don't write depth for transparent
      } else {
          gl.depthMask(true);
      }

      gl.drawElements(gl.TRIANGLES, mesh.buffer.numElements, gl.UNSIGNED_SHORT, 0);
    }
    
    requestAnimationFrame(render);
  }
  requestAnimationFrame(render);
};
