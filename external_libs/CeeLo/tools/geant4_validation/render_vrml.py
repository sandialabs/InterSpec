#!/usr/bin/env python3
"""Render the VRML that GEANT4's VRML2FILE driver writes, to a PNG.

The mesh in a g4_NN.wrl is GEANT4's own polyhedron for each solid, so looking at
it answers "did GEANT4 build the shape the exporter intended?" -- a question the
GDML text cannot answer on its own. Pair it with
macros/vis_crystal_vrml.mac:

    ceelo_g4val detector_26.gdml macros/vis_crystal_vrml.mac out.csv --vis-batch
    python3 render_vrml.py g4_00.wrl -o crystal.png --view side

Dependency-free apart from Pillow (only for writing the PNG); the projection,
hidden-surface ordering and shading are done here.  Painter's algorithm, which
is exact enough for one convex-ish solid and needs no z-buffer.

Views: side (looking along -x, z to the right), front, back, iso.  Add --skip to
drop volumes by name substring (the world box is dropped by default).
"""
import argparse, math, re, sys


def parse_vrml(path):
    """Return [(name, points, faces)] for each Shape in the file."""
    txt = open(path).read()
    shapes = []
    # Each solid is preceded by a "#---------- SOLID: <name>" comment.
    chunks = re.split(r'#-+ SOLID: ', txt)[1:]
    for chunk in chunks:
        name = chunk.split('\n', 1)[0].strip()
        pm = re.search(r'point\s*\[(.*?)\]', chunk, re.S)
        im = re.search(r'coordIndex\s*\[(.*?)\]', chunk, re.S)
        if not pm or not im:
            continue
        pts = []
        for row in pm.group(1).replace('\n', ' ').split(','):
            row = row.strip()
            if not row:
                continue
            x, y, z = (float(v) for v in row.split())
            pts.append((x, y, z))
        idx = [int(v) for v in im.group(1).replace(',', ' ').split()]
        faces, cur = [], []
        for i in idx:
            if i == -1:
                if len(cur) >= 3:
                    faces.append(cur)
                cur = []
            else:
                cur.append(i)
        if len(cur) >= 3:
            faces.append(cur)
        shapes.append((name, pts, faces))
    return shapes


VIEWS = {
    # (eye direction the camera looks along, up vector)
    'side':  ((-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    'front': ((0.0, 0.0, 1.0),  (0.0, 1.0, 0.0)),
    'back':  ((0.0, 0.0, -1.0), (0.0, 1.0, 0.0)),
    'iso':   ((-0.6, -0.55, 0.58), (0.0, 0.0, 1.0)),
}


def normalize(v):
    n = math.sqrt(sum(c * c for c in v))
    return tuple(c / n for c in v) if n else v


def cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])


def clip_faces(shapes, clip):
    """Drop the half of every mesh nearest the camera, so the cut reveals the
    inside (bore, tip) instead of the outer skin.  `clip` is like 'y+', meaning
    keep facets whose centroid has y >= 0."""
    if not clip:
        return shapes
    axis = 'xyz'.index(clip[0])
    keep_positive = clip[1] == '+'
    out = []
    for name, pts, faces in shapes:
        kept = []
        for f in faces:
            c = sum(pts[i][axis] for i in f) / len(f)
            if (c >= 0) == keep_positive:
                kept.append(f)
        if kept:
            out.append((name, pts, kept))
    return out


def render(shapes, out_path, view='side', size=(1000, 750), margin=0.06,
           wire=False):
    from PIL import Image, ImageDraw

    fwd = normalize(VIEWS[view][0])          # direction the camera looks
    up0 = VIEWS[view][1]
    right = normalize(cross(fwd, up0))
    up = normalize(cross(right, fwd))

    def project(p):
        """Orthographic: (screen x, screen y, depth along view)."""
        return (sum(p[i] * right[i] for i in range(3)),
                sum(p[i] * up[i] for i in range(3)),
                sum(p[i] * fwd[i] for i in range(3)))

    tris = []
    for name, pts, faces in shapes:
        proj = [project(p) for p in pts]
        for f in faces:
            poly = [proj[i] for i in f]
            world = [pts[i] for i in f]
            # Face normal from the first three world vertices, for shading.
            a = tuple(world[1][i] - world[0][i] for i in range(3))
            b = tuple(world[2][i] - world[0][i] for i in range(3))
            n = normalize(cross(a, b))
            depth = sum(v[2] for v in poly) / len(poly)
            tris.append((depth, poly, n, name))

    if not tris:
        sys.exit('nothing to draw')

    xs = [v[0] for _, poly, _, _ in tris for v in poly]
    ys = [v[1] for _, poly, _, _ in tris for v in poly]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    W, H = size
    span = max(x1 - x0, y1 - y0) * (1 + 2 * margin)
    scale = min(W, H) / span
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2

    img = Image.new('RGB', (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    # Far faces first: for a single closed solid this is exactly right.
    tris.sort(key=lambda t: t[0], reverse=True)
    light = normalize((-0.4, -0.5, 0.75))

    for _, poly, n, name in tris:
        # Two-sided lighting: a cut view shows interior faces, whose normals
        # point away from the camera, and they must not go black.
        shade = 0.34 + 0.66 * abs(sum(n[i] * light[i] for i in range(3)))
        base = (90, 150, 215)
        col = tuple(int(min(255, c * shade)) for c in base)
        pix = [((v[0] - cx) * scale + W / 2, H / 2 - (v[1] - cy) * scale)
               for v in poly]
        # Outlines help on a coarse mesh but turn a 3000-facet polycone into
        # solid grey, so they are opt-in.
        draw.polygon(pix, fill=col, outline=(40, 55, 75) if wire else col)

    img.save(out_path)
    return len(tris)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('wrl')
    ap.add_argument('-o', '--out', default='geometry.png')
    ap.add_argument('--view', default='side', choices=sorted(VIEWS))
    ap.add_argument('--skip', action='append', default=['World'],
                    help='drop volumes whose name contains this (repeatable)')
    ap.add_argument('--only', default=None,
                    help='draw only volumes whose name contains this')
    ap.add_argument('--clip', default=None,
                    choices=['x+', 'x-', 'y+', 'y-', 'z+', 'z-'],
                    help='keep only facets on one side of a plane, so the cut '
                         'shows the interior (e.g. y+ with --view side)')
    ap.add_argument('--wire', action='store_true',
                    help='outline every facet (readable only on coarse meshes)')
    args = ap.parse_args()

    shapes = parse_vrml(args.wrl)
    kept = []
    for name, pts, faces in shapes:
        if args.only and args.only not in name:
            continue
        if any(s in name for s in args.skip):
            continue
        kept.append((name, pts, faces))

    kept = clip_faces(kept, args.clip)
    n = render(kept, args.out, view=args.view, wire=args.wire)
    print(f'{args.out}: {len(kept)} volume(s), {n} facets, view={args.view}')
    for name, pts, faces in kept:
        print(f'    {name}: {len(pts)} vertices, {len(faces)} facets')


if __name__ == '__main__':
    main()
