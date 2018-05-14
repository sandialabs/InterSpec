/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 19 Dec 11 th
*/

#define DIVONNE
#define ROUTINE "Divonne"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata,
  creal epsrel, creal epsabs,
  cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  creal border, creal maxchisq, creal mindeviation,
  cnumber ngiven, ccount ldxgiven, real *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  t.ndim = ndim;
  t.ncomp = ncomp;
  t.integrand = integrand;
  t.userdata = userdata;
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
  t.ngiven = ngiven;
  t.xgiven = xgiven;
  t.ldxgiven = ldxgiven;
  t.nextra = nextra;
  t.peakfinder = peakfinder;
  t.nregions = 0;
  t.neval = 0;

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}

/*********************************************************************/

Extern void EXPORT(divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  creal *pborder, creal *pmaxchisq, creal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, real *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  This t;
  t.ndim = *pndim;
  t.ncomp = *pncomp;
  t.integrand = integrand;
  t.userdata = userdata;
  t.epsrel = *pepsrel;
  t.epsabs = *pepsabs;
  t.flags = *pflags;
  t.seed = *pseed;
  t.mineval = *pmineval;
  t.maxeval = *pmaxeval;
  t.key1 = *pkey1;
  t.key2 = *pkey2;
  t.key3 = *pkey3;
  t.maxpass = *pmaxpass;
  t.border.upper = 1 - (t.border.lower = *pborder);
  t.maxchisq = *pmaxchisq;
  t.mindeviation = *pmindeviation;
  t.ngiven = *pngiven;
  t.xgiven = xgiven;
  t.ldxgiven = *pldxgiven;
  t.nextra = *pnextra;
  t.peakfinder = peakfinder;
  t.nregions = 0;
  t.neval = 0;

  *pfail = Integrate(&t, integral, error, prob);
  *pnregions = t.nregions;
  *pneval = t.neval;
}

