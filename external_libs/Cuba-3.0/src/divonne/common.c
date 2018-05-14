/*
	common.c
		includes most of the modules
		this file is part of Divonne
		last modified 19 Dec 11 th
*/


#include "Random.c"
#include "ChiSquare.c"
#include "Rule.c"
#include "Sample.c"
#include "FindMinimum.c"
#include "Split.c"
#include "Explore.c"
#include "Iterate.c"

static inline bool BadDimension(cThis *t, ccount key)
{
  if( t->ndim > NDIM ) return true;
  if( IsSobol(key) ) return
    t->ndim < SOBOL_MINDIM || (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
  if( IsRule(key, t->ndim) ) return t->ndim < 1;
  return t->ndim < KOROBOV_MINDIM || t->ndim > KOROBOV_MAXDIM;
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

static inline void AllocGiven(This *t)
{
  real *xgiven = NULL, *fgiven = NULL;

  if( t->ngiven | t->nextra ) {
    cnumber nxgiven = t->ngiven*t->ndim;
    cnumber nxextra = t->nextra*t->ndim;
    cnumber nfgiven = t->ngiven*t->ncomp;
    cnumber nfextra = t->nextra*t->ncomp;

    Alloc(xgiven, nxgiven + nxextra + nfgiven + nfextra);
    t->xextra = xgiven + nxgiven;
    fgiven = t->xextra + nxextra;
    t->fextra = fgiven + nfgiven;

    if( nxgiven ) {
#ifdef MLVERSION
      Copy(xgiven, t->xgiven, nxgiven);
      Copy(fgiven, t->fgiven, nfgiven);
#else
      if( t->ldxgiven == t->ndim )
        Copy(xgiven, t->xgiven, nxgiven);
      else {
        number i;
        real *sgiven = t->xgiven, *dgiven = xgiven;
        for( i = 0; i < t->ngiven; ++i ) {
          Copy(dgiven, sgiven, t->ndim);
          sgiven += t->ldxgiven;
          dgiven += t->ndim;
        }
      }
      t->phase = 0;
      DoSample(t, t->ngiven, xgiven, fgiven);
#endif
    }
  }

  t->xgiven = xgiven;
  t->fgiven = fgiven;
}

