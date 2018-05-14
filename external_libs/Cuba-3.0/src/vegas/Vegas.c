/*
	Vegas.c
		Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 19 Dec 11 th
*/


#define VEGAS
#define ROUTINE "Vegas"

#include "decl.h"
#include "CSample.c"

/*********************************************************************/

Extern void EXPORT(Vegas)(ccount ndim, ccount ncomp,
  Integrand integrand, void *userdata,
  creal epsrel, creal epsabs, cint flags, cint seed,
  cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease, cnumber nbatch,
  cint gridno, cchar *statefile,
  number *pneval, int *pfail,
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
  t.nstart = nstart;
  t.nincrease = nincrease;
  t.nbatch = nbatch;
  t.gridno = gridno;
  t.statefile = statefile;
  t.neval = 0;

  *pfail = Integrate(&t, integral, error, prob);
  *pneval = t.neval;
}

/*********************************************************************/

Extern void EXPORT(vegas)(ccount *pndim, ccount *pncomp,
  Integrand integrand, void *userdata,
  creal *pepsrel, creal *pepsabs, cint *pflags, cint *pseed,
  cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnstart, cnumber *pnincrease, 
  cnumber *pnbatch, cint *pgridno, cchar *statefile,
  number *pneval, int *pfail,
  real *integral, real *error, real *prob, int statefilelen)
{
  char *s = NULL;
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
  t.nstart = *pnstart;
  t.nincrease = *pnincrease;
  t.nbatch = *pnbatch;
  t.gridno = *pgridno;
  t.neval = 0;

  if( statefile ) {
	/* strip trailing spaces */
    while( statefilelen > 0 && statefile[statefilelen - 1] == ' ' )
      --statefilelen;
    if( statefilelen > 0 && (s = malloc(statefilelen + 1)) ) {
      memcpy(s, statefile, statefilelen);
      s[statefilelen] = 0;
    }
  }
  t.statefile = s;

  *pfail = Integrate(&t, integral, error, prob);
  *pneval = t.neval;

  free(s);
}

