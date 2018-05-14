:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Vegas::usage = "Vegas[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: NStart::usage = "NStart is an option of Vegas.
	It specifies the number of integrand evaluations per iteration to start with."

:Evaluate: NIncrease::usage = "NIncrease is an option of Vegas.
	It specifies the increase in the number of integrand evaluations per iteration."

:Evaluate: NBatch::usage = "NBatch is an option of Vegas.
	It specifies how many points are sent in one MathLink packet to be sampled by Mathematica."

:Evaluate: GridNo::usage = "GridNo is an option of Vegas.
	Vegas maintains an internal table in which it can memorize up to 10 grids, to be used on subsequent integrations.
	A GridNo between 1 and 10 selects the slot in this internal table.
	For other values the grid is initialized from scratch and discarded at the end of the integration."

:Evaluate: StateFile::usage = "StateFile is an option of Vegas.
	It specifies a file in which the internal state is stored after each iteration and from which it can be restored on a subsequent run.
	The state file is removed once the prescribed accuracy has been reached."

:Evaluate: MinPoints::usage = "MinPoints is an option of Vegas.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Vegas.
	It can take the values Last or All which determine whether only the last (largest) or all of the samples collected on a subregion over the iterations contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Vegas.
	It can take the following values:
	False for Sobol quasi-random numbers (default),
	True or 0 for Mersenne Twister pseudo-random numbers,
	any other integer value n for Ranlux pseudo-random numbers of luxury level n."

:Evaluate: PseudoRandomSeed::usage = "PseudoRandomSeed is an option of Vegas.
	It specifies the seed for the pseudo-random number generator."

:Evaluate: SharpEdges::usage = "SharpEdges is an option of Vegas.
	It turns off smoothing of the importance function for integrands with sharp edges."

:Evaluate: RetainStateFile::usage = "RetainStateFile is an option of Vegas.
	It determines whether a chosen state file is kept even if the integration terminates normally."

:Evaluate: $Weight::usage = "$Weight is a global variable set by Vegas during the evaluation of the integrand to the weight of the point being sampled."

:Evaluate: $Iteration::usage = "$Iteration is a global variable set by Suave during the evaluation of the integrand to the present iteration number."

:Evaluate: MapSample::usage = "MapSample is a function used to map the integrand over the points to be sampled."


:Evaluate: Begin["`Vegas`"]

:Begin:
:Function: Vegas
:Pattern: MLVegas[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, seed_,
  mineval_, maxeval_,
  nstart_, nincrease_, nbatch_,
  gridno_, statefile_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, seed,
  mineval, maxeval,
  nstart, nincrease, nbatch,
  gridno, statefile}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer,
  Integer, Integer,
  Integer, Integer, Integer,
  Integer, String}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Vegas] = {HoldFirst}

:Evaluate: Options[Vegas] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000,
	NStart -> 1000, NIncrease -> 500,
	NBatch -> 1000, GridNo -> 0, StateFile -> "",
	Verbose -> 1, Final -> All,
	PseudoRandom -> False, PseudoRandomSeed -> 5489,
	SharpEdges -> False, RetainStateFile -> False,
	Compiled -> True}

:Evaluate: Vegas[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}], ncomp,
	tags, vars, lower, range, jac, tmp, defs, intT,
	rel, abs, mineval, maxeval, nstart, nincrease, nbatch,
	gridno, verbose, final, level, seed, edges, retain,
	compiled, $Weight, $Iteration},
	  Message[Vegas::optx, #, Vegas]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Vegas]];
	  {rel, abs, mineval, maxeval, nstart, nincrease, nbatch,
	    gridno, state, verbose, final, level, seed, edges, retain,
	    compiled} = tags /. {opt} /. Options[Vegas];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, Thread[vars = defs], jac]];
	  intT = integrandT[f];
	  Block[#,
	    ncomp = Length[intT@@ RandomReal[1, ndim]];
	    MLVegas[ndim, ncomp, 10.^-rel, 10.^-abs,
	      Min[Max[verbose, 0], 3] +
	        If[final === Last, 4, 0] +
	        If[TrueQ[edges], 8, 0] +
	        If[TrueQ[retain], 16, 0] +
	        If[IntegerQ[level], 256 level, 0],
	      If[level =!= False && IntegerQ[seed], seed, 0],
	      mineval, maxeval,
	      nstart, nincrease, nbatch,
	      gridno, state]
	  ]& @ vars
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Vegas`t" <> ToString[n]]

:Evaluate: Attributes[foo] = {HoldAll}

:Evaluate: define[True, tmp_, defs_, jac_] :=
	integrandT[f_] := Compile[tmp, eval[defs, Chop[f jac]//N],
	  {{_eval, _Real, 1}}]

:Evaluate: define[_, tmp_, defs_, jac_] :=
	integrandT[f_] := Function[tmp, eval[defs, Chop[f jac]//N]]

:Evaluate: eval[_, f_Real] = {f}

:Evaluate: eval[_, f:{__Real}] = f

:Evaluate: eval[x_, _] := (Message[Vegas::badsample, ff, x]; {})

:Evaluate: sample[x_, w_, iter_] := (
	$Iteration = iter;
	Check[Flatten @ MapSample[
	  ($Weight = #[[1]]; intT@@ #[[2]])&,
	  Transpose[{w, Partition[x, ndim]}] ], {}] )

:Evaluate: MapSample = Map

:Evaluate: Vegas::badsample = "`` is not a real-valued function at ``."

:Evaluate: Vegas::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Vegas::badcomp = "Cannot integrate `` components."

:Evaluate: Vegas::accuracy =
	"Desired accuracy was not reached within `` function evaluations."

:Evaluate: Vegas::success = "Needed `` function evaluations."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Vegas.tm
		Vegas Monte Carlo integration
		by Thomas Hahn
		last modified 17 Apr 12 th
*/


#define VEGAS
#define ROUTINE "Vegas"

#include "mathlink.h"
#include "decl.h"
#include "MSample.c"

/*********************************************************************/

static void Status(MLCONST char *msg, cint n)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 2);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Vegas");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n);
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
      Status("baddim", t->ndim);
      break;
    case -2:
      Status("badcomp", t->ncomp);
      break;
    }
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    Status(fail ? "accuracy" : "success", t->neval);
    MLPutFunction(stdlink, "Thread", 1);
    MLPutFunction(stdlink, "List", 3);
    MLPutRealList(stdlink, integral, t->ncomp);
    MLPutRealList(stdlink, error, t->ncomp);
    MLPutRealList(stdlink, prob, t->ncomp);
  }
}

/*********************************************************************/

void Vegas(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease, cint nbatch,
  cint gridno, cchar *statefile)
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
  t.nstart = nstart;
  t.nincrease = nincrease;
  t.nbatch = nbatch;
  t.gridno = gridno;
  t.statefile = statefile;
  t.neval = 0;

  DoIntegrate(&t);
  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

