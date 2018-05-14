(* Test suite of Genz, used also by Sloan and Joe, and Novak and Ritter *)

seed = 4711

maxpoints = 150000

repeat = 20


(* Family 1: Oscillatory *)

f[1][x_, c_, w_] := Cos[2 Pi w[[1]] + c.x]


(* Family 2: Product peak *)

f[2][x_, c_, w_] := Times@@ MapThread[f2a, {x, c, w}]

f2a[xi_, ci_, wi_] := 1/(ci^-2 + (xi - wi)^2)


(* Family 3: Corner peak *)

f[3][x_, c_, w_] := (1 + c.x)^(-(Length[x] + 1))


(* Family 4: Gaussian *)

f[4][x_, c_, w_] := Exp[Plus@@ MapThread[f4a, {x, c, w}]]

f4a[xi_, ci_, wi_] := -ci^2 (xi - wi)^2


(* Family 5: Exponential *)

f[5][x_, c_, w_] := Exp[Plus@@ MapThread[f5a, {x, c, w}]]

f5a[xi_, ci_, wi_] := -ci Abs[xi - wi]


(* Family 6: Discontinuous *)

f[6][x_, c_, w_] := 0 /; x[[1]] > w[[1]] || x[[2]] > w[[2]]

f[6][x_, c_, w_] := Exp[c.x]


(* Novak & Ritter use
difficulty[fam_] := {9.00, 7.25, 1.85, 7.03, 2.04, 4.30}[[fam]]
*)

(* Sloan & Joe use
scale[dim_] := dim^Min[Max[.2 dim, 1], 2]

SetOptions[Interpolation, InterpolationOrder -> 2]

ifun[1] = Interpolation[{{5, 145.7}, {8, 354.0}, {10,  900.0}}];
ifun[2] = Interpolation[{{5, 261.0}, {8, 545.0}, {10, 1760.0}}];
ifun[3] = Interpolation[{{5, 433.0}, {8, 193.0}, {10,  185.0}}];
ifun[4] = Interpolation[{{5, 155.0}, {8, 382.0}, {10, 1230.0}}];
ifun[5] = Interpolation[{{5, 217.0}, {8, 674.0}, {10, 2040.0}}];
ifun[6] = Interpolation[{{5,  90.0}, {8, 240.0}, {10, 1470.0}}];

difficulty[fam_] := ifun[fam][ndim]/scale[ndim]
*)

difficulty[fam_] := {6.0, 18.0, 2.2, 15.2, 16.1, 16.4}[[fam]]

c[fam_] := Block[{r = w}, r difficulty[fam]/Plus@@ r]


Install["Vegas"]

Install["Suave"]

Install["Divonne"]

Install["Cuhre"]



SetAll[opt__] := (
  SetOptions[Vegas, opt];
  SetOptions[Suave, opt];
  SetOptions[Divonne, opt];
  SetOptions[Cuhre, opt];
  SetOptions[NIntegrate, opt];
)

SetAll[PrecisionGoal -> 3, MaxPoints -> maxpoints]

SetOptions[Divonne, Key1 -> -200]


def[f_][{x__}][{r__}] := (
  Attributes[idef] = {HoldAll};
  idef[int_, NIntegrate] := (
    int := Module[{count = 0, res},
             res = NIntegrate[f, r, EvaluationMonitor :> (++count)];
             {count, res}]
  ) /; $VersionNumber >= 5;
  idef[int_, Int_] := (
    int := Module[{count = 0, res},
             res = Int[(++count; f), r];
             {count, res}]
  );
  idef[vegas, Vegas];
  idef[suave, Suave];
  idef[divonne, Divonne];
  idef[cuhre, Cuhre];
  idef[nint, NIntegrate];
)

vars = Table[Unique["x"], {20}]

test[ndim_, fam_] :=
Block[ {w, xs = Take[vars, ndim]},
  w := Table[Random[], {ndim}];
  def[f[fam][xs, c[fam], w]][xs][{#, 0, 1}&/@ xs];
  {vegas, suave, divonne, cuhre, nint}
]


dotest[ndim_, from_:1, to_:6] :=
Block[ {dir = ToString[ndim]},
  If[ FileType[dir] =!= Directory, CreateDirectory[dir] ];
  Do[
    SeedRandom[seed];
    Put[ Table[test[ndim, fam], {repeat}],
         ToFileName[dir, "fam" <> ToString[fam]] ],
  {fam, from, to}]
]

