/*
	MSample.c
		the sampling routine for the
		Mathematica versions of the Cuba routines
		by Thomas Hahn
		last modified 19 Mar 12 th
*/


static void DoSample(This *t, cnumber n, real *x, real *f
  VES_ONLY(, real *w, ccount iter))
{
  real *mma_f;
  long mma_n;

  if( MLAbort ) longjmp(t->abort, -99);

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`" ROUTINE "`sample", 1 VES_ONLY(+2) DIV_ONLY(+1));
  MLPutRealList(stdlink, x, n*t->ndim);
  VES_ONLY(MLPutRealList(stdlink, w, n);
           MLPutInteger(stdlink, iter);)
  DIV_ONLY(MLPutInteger(stdlink, t->phase);)
  MLEndPacket(stdlink);

  MLNextPacket(stdlink);
  if( !MLGetRealList(stdlink, &mma_f, &mma_n) ) {
    MLClearError(stdlink);
    MLNewPacket(stdlink);
    longjmp(t->abort, -99);
  }
 
  t->neval += mma_n;

  if( mma_n != n*t->ncomp ) {
    MLDisownRealList(stdlink, mma_f, mma_n);
    longjmp(t->abort, -3);
  }
 
  Copy(f, mma_f, n*t->ncomp);
  MLDisownRealList(stdlink, mma_f, mma_n);
}

/*********************************************************************/

#ifdef DIVONNE
#define Explore ExploreSerial

static count SampleExtra(This *t, cBounds *b)
{
  count n, nget;
  real *mma_f;
  long mma_n;

  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, "Cuba`Divonne`findpeak", 2);
  MLPutRealList(stdlink, (real *)b, 2*t->ndim);
  MLPutInteger(stdlink, t->phase);
  MLEndPacket(stdlink);

  MLNextPacket(stdlink);
  if( !MLGetRealList(stdlink, &mma_f, &mma_n) ) {
    MLClearError(stdlink);
    MLNewPacket(stdlink);
    longjmp(t->abort, -99);
  }

  t->neval += nget = mma_n/(t->ndim + t->ncomp);

  n = IMin(nget, t->nextra);
  if( n ) {
    Copy(t->xextra, mma_f, n*t->ndim);
    Copy(t->fextra, mma_f + nget*t->ndim, n*t->ncomp);
  }

  MLDisownRealList(stdlink, mma_f, mma_n);

  return n;
}
#endif

/*********************************************************************/

#include "common.c"

#define ForkCores(t)
#define WaitCores(t)

#include "Integrate.c"

