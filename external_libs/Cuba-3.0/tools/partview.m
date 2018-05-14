(*
	partview.m
		A partition viewer for Cuba results in Mathematica
		last modified 4 Feb 05 th
*)


BeginPackage["Cuba`"]

PartView::usage = "For a Cuba result obtained with Regions -> True,
PartView[result, dimx, dimy] displays the dimx-dimy plane of the
tessellation used in the integration."

Rect::usage = "Rect[{x1, y1}, {x2, y2}] is the graphics primitive used
by PartView to render a rectangle."

Begin["`PartView`"]

PartView[expr_, dimx_Integer, dimy_Integer] :=
Block[ {r, g, maxarea},
  r = Cases[expr,
    Region[ll_, ur_, ___] :> {ll[[{dimx, dimy}]], ur[[{dimx, dimy}]]},
    Infinity];
  maxarea = Times@@ (Max/@ #2 - Min/@ #1 &)@@ Transpose[r, {3, 1, 2}];
  (Show[#]; #)& @ Graphics[Apply[Rect, r, 1], AspectRatio -> 1]
]

maxarea = 1

Rect[{l_, d_}, {r_, u_}] := {
  { Hue[0, (ArcSin[1 - (r - l) (u - d)/maxarea]/(Pi/2))^2, 1],
    Rectangle[{l, d}, {r, u}] },
  { RGBColor[0, 0, 0],
    Line[{{l, d}, {r, d}, {r, u}, {l, u}, {l, d}}] }
}

End[]

EndPackage[]

