/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 19 Dec 11 th
*/


#include "ChiSquare.c"
#include "Rule.c"

#if ( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static inline
#endif
bool BadDimension(cThis *t)
{
  if( t->ndim > NDIM ) return true;
  return t->ndim < 2;
}

#if ( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static inline
#endif
bool BadComponent(cThis *t)
{
  if( t->ncomp > NCOMP ) return true;
  return t->ncomp < 1;
}

