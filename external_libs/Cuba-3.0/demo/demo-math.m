Install["Vegas"]

Install["Suave"]

Install["Divonne"]

Install["Cuhre"]


test[n_] := {t[n, Vegas], t[n, Suave], t[n, Divonne], t[n, Cuhre]}

t[n_, int_] := int[f[n][x, y, z], {x,0,1}, {y,0,1}, {z,0,1}]


f[1][x_, y_, z_] := Sin[x] Cos[y] Exp[z]

f[2][x_, y_, z_] := 1/((x + y)^2 + .003) Cos[y] Exp[z]

f[3][x_, y_, z_] := 1/(3.75 - Cos[Pi x] - Cos[Pi y] - Cos[Pi z])

f[4][x_, y_, z_] := Abs[x^2 + y^2 + z^2 - .125]

f[5][x_, y_, z_] := Exp[-x^2 - y^2 - z^2]

f[6][x_, y_, z_] := 1/(1 - x y z + 10^-10)

f[7][x_, y_, z_] := Sqrt[Abs[x - y - z]]

f[8][x_, y_, z_] := Exp[-x y z]

f[9][x_, y_, z_] := x^2/(Cos[x + y + z + 1] + 5)

f[10][x_, y_, z_] := If[ x > .5, 1/Sqrt[x y z + 10^-5], Sqrt[x y z] ]

f[11][x_, y_, z_] := If[ x^2 + y^2 + z^2 < 1, 1, 0 ]

