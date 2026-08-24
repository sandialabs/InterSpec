# GDML geometry viewer

A dependency-free WebGL viewer for the GDML files CeeLo exports
(`EfficiencyCalculator::export_geant4_gdml`). Useful for eyeballing a detector +
attenuator + source-shield stack before committing a long simulation or GEANT4
reference run — a surprising number of geometry mistakes are obvious the moment
you look at them.

It is a standalone front-end asset: nothing in the C++ library builds, links, or
requires it, and it is not part of any CMake target.

## Use

```bash
# 1. Export a geometry from any example that calls export_geant4_gdml(), e.g.
cd build/examples
./benchmark_mc_configs --config 7 --energies 662 --events 20000   # -> detector_7.gdml

# 2. Open the test page in a browser and load that file.
open ../../visualization/test_geometry_viewer.html
```

`test_geometry_viewer.html` is a self-contained harness: it loads
`GeometryViewer3D.js` from the same directory, has a file picker for a `.gdml`,
and exposes the camera/slice controls.

## Embedding

```html
<script src="GeometryViewer3D.js"></script>
<script>
  var viewer = new GeometryViewer3D(document.getElementById('viewer'), {});
  viewer.loadGDML(gdmlXmlString);      // string contents of a .gdml file
  viewer.setSourcePosition(0, 0, -10); // cm, in the GDML frame
</script>
```

Pure WebGL and vanilla JavaScript — no three.js, no CDN, no build step. Supported
GDML solids are the ones CeeLo emits: boxes, tubes (including cup-shaped
attenuator subtraction solids), spheres, the Marinelli beaker's L-shaped source
region, and **polycones** — which is how a crystal with a bore and/or a
bulletized (rounded) front edge is exported, so those detectors draw with the
fillet and the bore tip in place. The tooltip reports the polycone's radius,
length, bulletizing radius and bore dimensions read back from the z-planes.

## Cross-checking against GEANT4

The viewer shows what *CeeLo exported*. To see what *GEANT4 built* from it — a
different question, and the one that catches export bugs — dump G4's own
polyhedron and render that:

```bash
ceelo_g4val detector_26.gdml \
    ../../tools/geant4_validation/macros/vis_crystal_vrml.mac out.csv --vis-batch
python3 ../../tools/geant4_validation/render_vrml.py g4_00.wrl \
    -o crystal.png --view side --clip x-      # --clip cuts away the near half
```
