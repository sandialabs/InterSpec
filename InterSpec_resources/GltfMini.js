/* GltfMini.js - just enough glTF 2.0 for InterSpec's hand-written WebGL views (see DrfAppearance):
 buffers embedded as data: URIs, triangle primitives with POSITION (+ optional NORMAL) and
 8/16/32-bit indices, node TRS/matrix transforms, and baseColorFactor materials.  No textures,
 skins, morph targets or animations.  ES5 only (Safari 13).

 GltfMini.parse(json) -> [{ name, extras, color: [r,g,b,a], position: Float32Array,
                            normal: Float32Array, indices: Uint16Array|Uint32Array }, ...]
 with positions and normals in the model's own frame (glTF metres), after each node's transform.
 A primitive without normals gets smooth ones, averaged from its triangles.
*/
var GltfMini = (function() {
  function decodeBase64(uri) {
    var b64 = uri.substring(uri.indexOf(',') + 1);
    var bin = atob(b64);
    var bytes = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    return bytes.buffer;
  }

  var TYPE_SIZE = { SCALAR: 1, VEC2: 2, VEC3: 3, VEC4: 4, MAT4: 16 };

  function readAccessor(doc, buffers, index) {
    var acc = doc.accessors[index];
    var view = doc.bufferViews[acc.bufferView];
    var n = acc.count * TYPE_SIZE[acc.type];
    var offset = (view.byteOffset || 0) + (acc.byteOffset || 0);
    var buf = buffers[view.buffer];
    switch (acc.componentType) {
      case 5126: return new Float32Array(buf.slice(offset, offset + 4*n));
      case 5125: return new Uint32Array(buf.slice(offset, offset + 4*n));
      case 5123: return new Uint16Array(buf.slice(offset, offset + 2*n));
      case 5121: return new Uint8Array(buf, offset, n);
    }
    throw new Error('GltfMini: unsupported component type ' + acc.componentType);
  }

  // Column-major 4x4 helpers
  function identity() { return [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]; }
  function mul(a, b) {
    var o = new Array(16);
    for (var c = 0; c < 4; c++)
      for (var r = 0; r < 4; r++) {
        var s = 0;
        for (var k = 0; k < 4; k++) s += a[k*4 + r] * b[c*4 + k];
        o[c*4 + r] = s;
      }
    return o;
  }
  function trs(node) {
    if (node.matrix) return node.matrix.slice();
    var t = node.translation || [0,0,0], q = node.rotation || [0,0,0,1], s = node.scale || [1,1,1];
    var x = q[0], y = q[1], z = q[2], w = q[3];
    return [
      (1 - 2*(y*y + z*z))*s[0], (2*(x*y + z*w))*s[0],     (2*(x*z - y*w))*s[0],     0,
      (2*(x*y - z*w))*s[1],     (1 - 2*(x*x + z*z))*s[1], (2*(y*z + x*w))*s[1],     0,
      (2*(x*z + y*w))*s[2],     (2*(y*z - x*w))*s[2],     (1 - 2*(x*x + y*y))*s[2], 0,
      t[0], t[1], t[2], 1 ];
  }

  // Smooth vertex normals for a primitive that has none: the area-weighted sum of its triangles'.
  function computeNormals(P, idx) {
    var N = new Float32Array(P.length);
    for (var t = 0; t + 2 < idx.length; t += 3) {
      var a = 3*idx[t], b = 3*idx[t+1], c = 3*idx[t+2];
      var ux = P[b] - P[a], uy = P[b+1] - P[a+1], uz = P[b+2] - P[a+2];
      var vx = P[c] - P[a], vy = P[c+1] - P[a+1], vz = P[c+2] - P[a+2];
      var nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
      [a, b, c].forEach(function(v) { N[v] += nx; N[v+1] += ny; N[v+2] += nz; });
    }
    return N;
  }

  function normalize(N) {
    for (var i = 0; i < N.length; i += 3) {
      var l = Math.sqrt(N[i]*N[i] + N[i+1]*N[i+1] + N[i+2]*N[i+2]) || 1;
      N[i] /= l; N[i+1] /= l; N[i+2] /= l;
    }
  }

  function parse(doc) {
    if (!doc.asset || String(doc.asset.version).charAt(0) !== '2')
      throw new Error('GltfMini: not glTF 2.0');
    var buffers = (doc.buffers || []).map(function(b) {
      if (!b.uri || b.uri.indexOf('data:') !== 0)
        throw new Error('GltfMini: only embedded (data: URI) buffers are supported');
      return decodeBase64(b.uri);
    });

    var out = [];
    var visited = {};
    function visit(nodeIndex, parentMatrix) {
      // Each node once: a file listing the same children over and over (or a cycle) could
      //  otherwise make this walk take forever.
      var node = doc.nodes[nodeIndex];
      if (!node || visited[nodeIndex]) return;
      visited[nodeIndex] = true;
      var m = mul(parentMatrix, trs(node));
      if (node.mesh !== undefined) {
        var mesh = doc.meshes[node.mesh];
        mesh.primitives.forEach(function(prim) {
          if ((prim.mode !== undefined && prim.mode !== 4) || (prim.attributes.POSITION === undefined))
            return;  // triangles only
          var pos = readAccessor(doc, buffers, prim.attributes.POSITION);
          var nrm = (prim.attributes.NORMAL !== undefined) ? readAccessor(doc, buffers, prim.attributes.NORMAL) : null;
          var idx = (prim.indices !== undefined) ? readAccessor(doc, buffers, prim.indices) : null;
          if (!idx) {
            idx = new Uint32Array(pos.length / 3);
            for (var k = 0; k < idx.length; k++) idx[k] = k;
          } else if (idx instanceof Uint8Array) {
            idx = new Uint16Array(idx);
          }
          var P = new Float32Array(pos.length), N = nrm ? new Float32Array(nrm.length) : null;
          for (var i = 0; i < pos.length; i += 3) {
            for (var r = 0; r < 3; r++) {
              P[i + r] = m[r]*pos[i] + m[4 + r]*pos[i+1] + m[8 + r]*pos[i+2] + m[12 + r];
              if (N) N[i + r] = m[r]*nrm[i] + m[4 + r]*nrm[i+1] + m[8 + r]*nrm[i+2];  // assumes no shear
            }
          }
          if (!N) N = computeNormals(P, idx);
          normalize(N);
          var mat = (prim.material !== undefined) ? doc.materials[prim.material] : null;
          var color = (mat && mat.pbrMetallicRoughness && mat.pbrMetallicRoughness.baseColorFactor)
                      ? mat.pbrMetallicRoughness.baseColorFactor : [0.7, 0.7, 0.7, 1.0];
          out.push({ name: node.name || mesh.name || '', extras: node.extras || {}, color: color,
                     position: P, normal: N, indices: idx });
        });
      }
      (node.children || []).forEach(function(c) { visit(c, m); });
    }

    var roots = [];
    if (doc.scenes && doc.scenes.length) {
      roots = doc.scenes[doc.scene || 0].nodes || [];
    } else {
      // No scene: every node that is nobody's child is a root.
      var isChild = {};
      (doc.nodes || []).forEach(function(n) { (n.children || []).forEach(function(c) { isChild[c] = true; }); });
      for (var ni = 0; ni < (doc.nodes || []).length; ni++)
        if (!isChild[ni]) roots.push(ni);
    }
    roots.forEach(function(n) { visit(n, identity()); });
    return out;
  }

  return { parse: parse };
})();
