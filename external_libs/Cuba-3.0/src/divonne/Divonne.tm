:Evaluate: BeginPackage["Cuba`"]

:Evaluate: Divonne::usage =
	"Divonne[f, {x, xmin, xmax}..] computes a numerical approximation to the integral of the real scalar or vector function f.
	The output is a list with entries of the form {integral, error, chi-square probability} for each component of the integrand."

:Evaluate: Key1::usage = "Key1 is an option of Divonne.
	It determines sampling in the partitioning phase.\n
	Special cases:\n
	  Key1 = 7: use a degree-7 cubature rule,\n
	  Key1 = 9: use a degree-9 cubature rule,\n
	  Key1 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key1 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a random sample of n1 = Abs[Key1] points is used, where the sign of Key1 determines the type of sample:\n
	  Key1 > 0: use a Korobov quasi-random sample,\n
	  Key1 < 0: use a \"standard\" sample."

:Evaluate: Key2::usage = "Key2 is an option of Divonne.
	It determines sampling in the main integration phase.\n
	Special cases:\n
	  Key2 = 7: use a degree-7 cubature rule,\n
	  Key2 = 9: use a degree-9 cubature rule,\n
	  Key2 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key2 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a random sample is used, where the sign of Key2 determines the type of sample:\n
	  Key2 > 0: use a Korobov quasi-random sample,\n
	  Key2 < 0: use a \"standard\" sample,\n
	and n2 = Abs[Key2] determines the number of points:\n
	  n2 >= 40: sample n2 points,\n
	  n2 < 40: sample n2*nneed points, where nneed is the number of points needed to reach the prescribed accuracy, as estimated by Divonne from the results of the partitioning phase."

:Evaluate: Key3::usage = "Key3 is an option of Divonne.
	It sets the strategy for the refinement phase:\n
	  Key3 = 0: do not further treat the subregion,\n
	  Key3 = 1: split the subregion up once more,\n
	for other values the region is sampled a third time:\n
	  Key3 = 7: use a degree-7 cubature rule,\n
	  Key3 = 9: use a degree-9 cubature rule,\n
	  Key3 = 11: use a degree-11 cubature rule (available only in 3 dimensions),\n
	  Key3 = 13: use a degree-13 cubature rule (available only in 2 dimensions),\n
	otherwise a random sample is used, where the sign of Key3 determines the type of sample:\n
	  Key3 > 0: use a Korobov quasi-random sample,\n
	  Key3 < 0: use a \"standard\" sample,\n
	and n3 = Abs[Key3] determines the number of points:\n
	  n3 >= 40: sample n3 points,\n
	  n3 < 40: sample n3*nneed points, where nneed is the number of points needed to reach the prescribed accuracy, as estimated by Divonne from the results of the partitioning phase."

:Evaluate: MaxPass::usage = "MaxPass is an option of Divonne.
	It controls the partitioning termination.
	The partitioning phase is terminated when the estimated total number of integrand evaluations (partitioning plus main integration) does not decrease for MaxPass successive iterations."

:Evaluate: Border::usage = "Border is an option of Divonne.
	It specifies the width of the border of the integration region.
	Points falling into this border region are not sampled directly, but are extrapolated from two samples from the interior.
	The border width always refers to the unit hypercube, i.e. it is not rescaled if the integration region is not the unit hypercube."

:Evaluate: MaxChisq::usage = "MaxChisq is an option of Divonne.
	It specifies the maximum chi-square value a single subregion is allowed to have in the main integration phase.
	Regions which fail this chi-square test and whose sample averages differ by more than MinDeviation move on to the refinement phase."

:Evaluate: MinDeviation::usage = "MinDeviation is an option of Divonne.
	Regions which fail the chi-square test are not treated further if their sample averages differ by less than MinDeviation.
	MinDeviation is specified as the fraction of the requested error of the entire integral."

:Evaluate: Given::usage = "Given is an option of Divonne.
	It provides a list of points where the integrand might have peaks.
	Divonne will consider these points when partitioning the integration region."

:Evaluate: NExtra::usage = "NExtra is an option of Divonne.
	It specifies the maximum number of points that will be considered in the output of the PeakFinder function."

:Evaluate: PeakFinder::usage = "PeakFinder is an option of Divonne.
	It specifies the peak-finder function.
	This function is called whenever a region is up for subdivision and is supposed to point out possible peaks lying in the region, thus acting as the dynamic counterpart of the static list of points supplied with Given.
	It is invoked with two arguments, the multidimensional equivalents of the lower left and upper right corners of the region being investigated, and must return a (possibly empty) list of points."

:Evaluate: MinPoints::usage = "MinPoints is an option of Divonne.
	It specifies the minimum number of points to sample."

:Evaluate: Final::usage = "Final is an option of Divonne.
	It can take the values Last or All which determine whether only the last (largest) or all sets of samples collected on a subregion over the integration phases contribute to the final result."

:Evaluate: PseudoRandom::usage = "PseudoRandom is an option of Divonne.
	It can take the following values:
	False for Sobol quasi-random numbers (default),
	True or 0 for Mersenne Twister pseudo-random numbers,
	any other integer value n for Ranlux pseudo-random numbers of luxury level n."

:Evaluate: PseudoRandomSeed::usage = "PseudoRandomSeed is an option of Divonne.
	It specifies the seed for the pseudo-random number generator."

:Evaluate: Regions::usage = "Regions is an option of Divonne.
	It specifies whether the regions into which the integration region has been cut are returned together with the integration results."

:Evaluate: Region::usage = "Region[ll, ur, res, df] describes a subregion:
	ll and ur are multidimensional equivalents of the region's lower left and upper right corner.
	res gives the integration results for the region in a list with entries of the form {integral, error, chi-square} for each component of the integrand.
	df is the number of degrees of freedom corresponding to the chi-square values in res."

:Evaluate: $Phase::usage = "$Phase is a global variable set by Divonne during the evaluation of the integrand to the integration phase:\n
	0 = sampling of the points in xgiven,\n
	1 = partitioning phase,\n
	2 = main integration phase,\n
	3 = refinement phase."

:Evaluate: MapSample::usage = "MapSample is a function used to map the integrand over the points to be sampled."


:Evaluate: Begin["`Divonne`"]

:Begin:
:Function: Divonne
:Pattern: MLDivonne[ndim_, ncomp_,
  epsrel_, epsabs_, flags_, seed_,
  mineval_, maxeval_,
  key1_, key2_, key3_, maxpass_,
  border_, maxchisq_, mindeviation_,
  xgiven_, fgiven_, nextra_]
:Arguments: {ndim, ncomp,
  epsrel, epsabs, flags, seed,
  mineval, maxeval,
  key1, key2, key3, maxpass,
  border, maxchisq, mindeviation,
  xgiven, fgiven, nextra}
:ArgumentTypes: {Integer, Integer,
  Real, Real, Integer, Integer,
  Integer, Integer,
  Integer, Integer, Integer, Integer,
  Real, Real, Real,
  RealList, RealList, Integer}
:ReturnType: Manual
:End:

:Evaluate: Attributes[Divonne] = {HoldFirst}

:Evaluate: Options[Divonne] = {PrecisionGoal -> 3, AccuracyGoal -> 12,
	MinPoints -> 0, MaxPoints -> 50000,
	Key1 -> 47, Key2 -> 1, Key3 -> 1, MaxPass -> 5,
	Border -> 0, MaxChisq -> 10, MinDeviation -> .25,
	Given -> {}, NExtra -> 0, PeakFinder -> ({}&),
	Verbose -> 1, Final -> All,
	PseudoRandom -> False, PseudoRandomSeed -> 5489,
	Regions -> False, Compiled -> True}

:Evaluate: Divonne[f_, v:{_, _, _}.., opt___Rule] :=
	Block[ {ff = HoldForm[f], ndim = Length[{v}], ncomp,
	tags, vars, lower, range, jac, tmp, defs, intT, intX,
	rel, abs, mineval, maxeval, key1, key2, key3, maxpass, border,
	maxchisq, mindeviation,	given, nextra, peakfinder,
	final, verbose, level, seed, regions, compiled,
	$Phase},
	  Message[Divonne::optx, #, Divonne]&/@
	    Complement[First/@ {opt}, tags = First/@ Options[Divonne]];
	  {rel, abs, mineval, maxeval, key1, key2, key3, maxpass, border,
	    maxchisq, mindeviation, given, nextra, peakfinder,
	    verbose, final, level, seed, regions, compiled} =
	    tags /. {opt} /. Options[Divonne];
	  {vars, lower, range} = Transpose[{v}];
	  jac = Simplify[Times@@ (range -= lower)];
	  tmp = Array[tmpvar, ndim];
	  defs = Simplify[lower + range tmp];
	  Block[{Set}, define[compiled, tmp, vars, Thread[vars = defs], jac]];
	  intT = integrandT[f];
	  intX = integrandX[f];
	  Block[#,
	    ncomp = Length[intT@@ RandomReal[1, ndim]];
	    MLDivonne[ndim, ncomp, 10.^-rel, 10.^-abs,
	      Min[Max[verbose, 0], 3] +
	        If[final === Last, 4, 0] +
	        If[TrueQ[regions], 128, 0] +
	        If[IntegerQ[level], 256 level, 0],
	      If[level =!= False && IntegerQ[seed], seed, 0],
	      mineval, maxeval,
	      key1, key2, key3, maxpass,
	      N[border], N[maxchisq], N[mindeviation],
	      given, sample[given, 0, intX], nextra]
	  ]& @ vars
	]

:Evaluate: tmpvar[n_] := ToExpression["Cuba`Divonne`t" <> ToString[n]]

:Evaluate: Attributes[foo] = {HoldAll}

:Evaluate: define[True, tmp_, vars_, defs_, jac_] := (
	TtoX := TtoX = Compile[tmp, defs];
	integrandT[f_] := Compile[tmp, eval[defs, Chop[f jac]//N],
	  {{_eval, _Real, 1}}];
	integrandX[f_] := Compile[vars, eval[vars, Chop[f jac]//N],
	  {{_eval, _Real, 1}}] )

:Evaluate: define[_, tmp_, vars_, defs_, jac_] := (
	TtoX := TtoX = Function[tmp, defs];
	integrandT[f_] := Function[tmp, eval[defs, Chop[f jac]//N]];
	integrandX[f_] := Function[vars, eval[vars, Chop[f jac]//N]] )

:Evaluate: eval[_, f_Real] := {f}

:Evaluate: eval[_, f:{__Real}] := f

:Evaluate: eval[x_, _] := (Message[Divonne::badsample, ff, x]; {})

:Evaluate: sample[x_, p_, i_:intT] := (
	$Phase = p;
	Check[Flatten @ MapSample[i@@ # &, Partition[x, ndim]], {}] )

:Evaluate: MapSample = Map

:Evaluate: findpeak[b_, p_] := Check[Join[#, sample[#, p, intX]]& @
	N[Flatten[peakfinder@@ MapThread[TtoX, Partition[b, 2]]]], {}]

:Evaluate: region[ll_, ur_, r___] := Region[TtoX@@ ll, TtoX@@ ur, r]

:Evaluate: Divonne::badsample = "`` is not a real-valued function at ``."

:Evaluate: Divonne::baddim = "Cannot integrate in `` dimensions."

:Evaluate: Divonne::badcomp = "Cannot integrate `` components."

:Evaluate: Divonne::accuracy =
	"Desired accuracy was not reached within `` integrand evaluations on `` subregions.
	Estimate that MaxPoints needs to be increased by `` for this accuracy."

:Evaluate: Divonne::success = "Needed `` integrand evaluations on `` subregions."

:Evaluate: End[]

:Evaluate: EndPackage[]


/*
	Divonne.tm
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 19 Dec 11 th
*/


#define DIVONNE
#define ROUTINE "Divonne"

#include "mathlink.h"
#include "decl.h"
#include "MSample.c"

/*********************************************************************/

static void Status(MLCONST char *msg, cint n1, cint n2, cint n3)
{
  MLPutFunction(stdlink, "CompoundExpression", 2);
  MLPutFunction(stdlink, "Message", 4);
  MLPutFunction(stdlink, "MessageName", 2);
  MLPutSymbol(stdlink, "Divonne");
  MLPutString(stdlink, msg);
  MLPutInteger(stdlink, n1);
  MLPutInteger(stdlink, n2);
  MLPutInteger(stdlink, n3);
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
      Status("baddim", t->ndim, 0, 0);
      break;
    case -2:
      Status("badcomp", t->ncomp, 0, 0);
      break;
    }
    MLPutSymbol(stdlink, "$Failed");
  }
  else {
    Status(fail ? "accuracy" : "success", t->neval, t->nregions, fail);
    MLPutFunction(stdlink, "Thread", 1);
    MLPutFunction(stdlink, "List", 3);
    MLPutRealList(stdlink, integral, t->ncomp);
    MLPutRealList(stdlink, error, t->ncomp);
    MLPutRealList(stdlink, prob, t->ncomp);
  }
}

/*********************************************************************/

void Divonne(cint ndim, cint ncomp,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, cint maxpass,
  creal border, creal maxchisq, creal mindeviation,
  real *xgiven, clong nxgiven, real *fgiven, clong nfgiven,
  cnumber nextra)
{
  This t;
  t.ldxgiven = t.ndim = ndim;
  t.ncomp = ncomp;
  t.epsrel = epsrel;
  t.epsabs = epsabs;
  t.flags = flags;
  t.seed = seed;
  t.mineval = mineval;
  t.maxeval = maxeval;
  t.key1 = key1;
  t.key2 = key2;
  t.key3 = key3;
  t.maxpass = maxpass;
  t.border.upper = 1 - (t.border.lower = border);
  t.maxchisq = maxchisq;
  t.mindeviation = mindeviation;
  t.xgiven = xgiven;
  t.fgiven = fgiven;
  t.nextra = nextra;
  t.nregions = 0;
  t.neval = t.ngiven = nxgiven/ndim;

  DoIntegrate(&t);

  MLEndPacket(stdlink);
}

/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}

