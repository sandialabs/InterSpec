/*
	CSample.c
		the serial sampling routine
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 19 Dec 11 th
*/

#if( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static inline 
#endif
number SampleRaw(cThis *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  for( ; n; --n ) {
    if( t->integrand(&t->ndim, x, &t->ncomp, f, t->userdata
          VES_ONLY(, w++, &iter)
          DIV_ONLY(, &t->phase)) == ABORT ) return -1;
    x += t->ndim;
    f += t->ncomp;
  }
  return 0;
}

/*********************************************************************/
#if( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static inline 
#endif
void DoSampleSerial(This *t, cnumber n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  t->neval += n;
  if( SampleRaw(t, n, x, f VES_ONLY(, w, iter)) ) 
    longjmp(t->abort, -99);
}

/*********************************************************************/

#ifdef HAVE_FORK

static void DoSample(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter));
DIV_ONLY(static int Explore(This *t, cint iregion);)

#else

#define DoSample DoSampleSerial
#define Explore ExploreSerial
#define ForkCores(t)
#define WaitCores(t)

#endif

#ifdef DIVONNE
static inline count SampleExtra(This *t, cBounds *b)
{
  number n = t->nextra;
  t->peakfinder(&t->ndim, b, &n, t->xextra);
  DoSample(t, n, t->xextra, t->fextra);
  return n;
}
#endif

#include "common.c"

#ifdef HAVE_FORK
#include "Fork.c"
#endif

#include "Integrate.c"

