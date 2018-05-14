:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Suave::usage =
	"Suave[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: NNew::usage = "NNew is an option of Suave.
	It specifies the number of new integrand evaluations in each subdivision."

:Evaluate: Flatness::usage = "Flatness is an option of Suave.
	It determines how prominently individual samples with a large fluctuation figure in the total fluctuation, which in turn determines how a region is split up.
	Explicitly, if F[i] is the individual fluctuation of sample i, the total fluctuation is computed as Sum[(1 + F[i])^p, {i, nsamples}]^(2/3/p), i.e. as the p-norm of the fluctuation vector to the power 2/3, where p is the number given by Flatness.
	Thus with increasing p, the fluctuation becomes more and more dominated by outliers, i.e. points with a large fluctuation.
	As suggested by the name Flatness, p should be chosen large for `flat' integrands and small for `volatile' integrands with high peaks.
	Note that since p appears in the exponent, one should not use too large values (say, no more than a few hundred) lest terms be truncated internally to prevent overflow."

:Evaluate: MinPoints::usage = "MinPoints is an option of Suave.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Suave.
	It can take the values Last or All which determine whether only the last (largest) or all sets of samples collected on a subregion over the iterations contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Suave.
	It can take the following values:
	False for Sobol quasi-random numbers (default),
	True or 0 for Mersenne Twister pseudo-random numbers,
	any other integer value n for Ranlux pseudo-random numbers of luxury level n."

:Evaluate: PseudoRandomSeed::usage = "PseudoRandomSeed is an option of Suave.
	It specifies the seed for the pseudo-random number generator."

:Evaluate: SharpEdges::usage = "SharpEdges is an option of Suave.
	It turns off smoothing of the importance function for integrands with sharp edges."

:Evaluate: Regions::usage = "Regions is an option of Suave.
	It specifies whether the regions into which the integration region has been cut are returned together with the integration results."

:Evaluate: Region::usage = "Region[ll, ur, res, df] describes a subregion:
	ll and ur are multidimensional equivalents of the region's lower left and upper right corner.
	res gives the integration results for the region in a list with entries of the form {integral, error, chi-square} for each component of the integrand.
	df is the number of degrees of freedom corresponding to the chi-square values in res."

:Evaluate: $Weight::usage = "$Weight is a global variable set by Suave during the evaluation of the integrand to the weight of the point being sampled."

:Evaluate: $Iteration::usage = "$Iteration is a global variable set by Suave during the evaluation of the integrand to the present iteration number."

:Evaluate: MapSample::usage = "MapSample is a function used to map the integrand over the points to be sampled."


:Evaluate: Begin["`Suave`"]

:Begin:
:Function: Suave
:Pattern: MLSuave[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, seed_,
  mineval_, maxeval_,
  nnew_, flatness_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, seed,
  mineval, maxeval,
  nnew, flatness}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer,
  Integer, Integer,
  Integer, Real}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Suave] = {HoldFirst}

:Evaluate: Options[Suave] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000, NNew -> 1000, Flatness -> 50,
	Verbose -> 1, Final -> Last,
	PseudoRandom -> False, PseudoRandomSeed -> 5489,
	SharpEdges -> False, Regions -> False, Compiled -> True}

:Evaluate: Suave[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}], ncomp,
	tags, vars, lower, range, jac, tmp, defs, intT,
	rel, abs, mineval, maxeval, nnew, flatness,
	verbose, final, level, seed, edges, regions, compiled,
	$Weight, $Iteration},
	  Message[Suave::optx, #, Suave]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Suave]];
	  {rel, abs, mineval, maxeval, nnew, flatness,
	    verbose, final, level, seed, edges, regions, compiled} =
	    tags /. {opt} /. Options[Suave];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, Thread[vars = defs], jac]];
	  intT = integrandT[f];
	  Block[#,
	    ncomp = Length[intT@@ RandomReal[1, ndim]];
	    MLSuave[ndim, ncomp, 10.^-rel, 10.^-abs,
	      Min[Max[verbose, 0], 3] +
	        If[final === Last, 4, 0] +
	        If[TrueQ[edges], 8, 0] +
	        If[TrueQ[regions], 128, 0] +
	        If[IntegerQ[level], 256 level, 0],
	      If[level =!= False && IntegerQ[seed], seed, 0],
	      mineval, maxeval,
	      nnew, flatness]
	  ]& @ vars
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Suave`t" <> ToString[n]]

:Evaluate: Attributes[foo] = {HoldAll}

:Evaluate: define[True, tmp_, defs_, jac_] := (
	TtoX := TtoX = Compile[tmp, defs];
	integrandT[f_] := Compile[tmp, eval[defs, Chop[f jac]//N],
	  {{_eval, _Real, 1}}] )

:Evaluate: define[_, tmp_, defs_, jac_] := (
	TtoX := TtoX = Function[tmp, defs];
	integrandT[f_] := Function[tmp, eval[defs, Chop[f jac]//N]] )

:Evaluate: eval[_, f_Real] = {f}

:Evaluate: eval[_, f:{__Real}] = f

:Evaluate: eval[x_, _] := (Message[Suave::badsample, ff, x]; {})

:Evaluate: sample[x_, w_, iter_] := (
	$Iteration = iter;
	Check[Flatten @ MapSample[
	  ($Weight = #[[1]]; intT@@ #[[2]])&,
	  Transpose[{w, Partition[x, ndim]}] ], {}] )

:Evaluate: MapSample = Map

:Evaluate: region[ll_, ur_, r___] := Region[TtoX@@ ll, TtoX@@ ur, r]

:Evaluate: Suave::badsample = "`` is not a real-valued function at ``."

:Evaluate: Suave::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Suave::badcomp = "Cannot integrate `` components."

:Evaluate: Suave::accuracy =
	"Desired accuracy was not reached within `` function evaluations on `` subregions."

:Evaluate: Suave::success = "Needed `` function evaluations on `` subregions."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Suave.tm
		Subregion-adaptive Vegas Monte Carlo integration
		by Thomas Hahn
		last modified 19 Dec 11 th
*/


#define SUAVE
#define ROUTINE "Suave"

#include "mathlink.h"
#include "decl.h"
#include "MSample.c"

/*********************************************************************/

static void Status(MLCONST char *msg, cint n1, cint n2)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 3);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Suave");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n1);
  MLPutInteger(stdlink, n2);
}

/*********************************************************************/

static inline void DoIntegrate(This *t)
{
  real integral[NCOMP], error[NCOMP], prob[NCOMP];
  cint fail = Integrate(t, integral, error, prob);

  if( fail < 0 ) {
    switch( fail ) {
    case -99:
      MLPutFunction(stdlink, "Abort", 0);
      return;
    case -1:
      Status("baddim", t->ndim, 0);
      break;
    case -2:
      Status("badcomp", t->ncomp, 0);
      break;
    }
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    Status(fail ? "accuracy" : "success", t->neval, t->nregions);
    MLPutFunction(stdlink, "Thread", 1);
    MLPutFunction(stdlink, "List", 3);
    MLPutRealList(stdlink, integral, t->ncomp);
    MLPutRealList(stdlink, error, t->ncomp);
    MLPutRealList(stdlink, prob, t->ncomp);
  }
}

/*********************************************************************/

void Suave(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness)
{
  This t;
  t.ndim = ndim;
  t.ncomp = ncomp;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = flags;
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.nnew = nnew;
  t.flatness = flatness;
  t.nregions = 0;
  t.neval = 0;

  DoIntegrate(&t);
  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

