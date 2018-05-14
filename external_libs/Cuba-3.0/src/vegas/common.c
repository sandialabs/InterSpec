/*
	common.c
		Code common to Vegas.c and Vegas.tm
		this file is part of Vegas
		last modified 21 Dec 11 th
*/


#include "Random.c"
#include "ChiSquare.c"
#include "Grid.c"

static inline bool BadDimension(cThis *t)
{
  if( t->ndim > NDIM ) return true;
  return t->ndim < SOBOL_MINDIM ||
    (t->seed == 0 && t->ndim > SOBOL_MAXDIM);
}

static inline bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

