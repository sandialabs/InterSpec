/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 21 Dec 11 th
*/


#include "stddecl.h"

#define Tag(x) ((x) | INT_MIN)
#define Untag(x) ((x) & INT_MAX)

typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;

typedef struct {
  real avg, err;
} PhaseResult;

typedef struct {
  real avg, spreadsq;
  real spread, secondspread;
  real nneed, maxerrsq, mindevsq;
  PhaseResult phase[2];
  int iregion;
} Totals;

typedef struct {
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;

typedef struct samples {
  real *x, *f;
  void (*sampler)(struct _this *t, ccount);
  cRule *rule;
  number n, neff;
  count coeff;
} Samples;

typedef const Samples cSamples;

typedef struct {
  real diff, err, spread;
} Errors;

typedef const Errors cErrors;

typedef int (*Integrand)(ccount *, creal *, ccount *, real *, void *, cint *);

typedef void (*PeakFinder)(ccount *, cBounds *, number *, real *);

typedef struct _this {
  count ndim, ncomp;
#ifndef MLVERSION
  Integrand integrand;
  void *userdata;
  PeakFinder peakfinder;
#ifdef HAVE_FORK
  int ncores, *child;
  int running, nchildren;
  fd_set children;
  real *frame;
  number nframe;
  SHM_ONLY(int shmid;)
#endif
#endif
  real epsrel, epsabs;
  int flags, seed;
  number mineval, maxeval;
  int key1, key2, key3;
  count maxpass;
  Bounds border;
  real maxchisq, mindeviation;
  number ngiven, nextra;
  real *xgiven, *fgiven;
  real *xextra, *fextra;
  count ldxgiven;
  count nregions;
  number neval, neval_opt, neval_cut;
  count phase;
  count selectedcomp, size;
  Samples samples[3];
  Totals *totals;
  Rule rule7, rule9, rule11, rule13;
  RNGState rng;
  void *voidregion;
  jmp_buf abort;
} This;

typedef const This cThis;

#define TYPEDEFREGION \
  typedef struct { \
    real avg, err, spread, chisq; \
    real fmin, fmax; \
    real xmin[NDIM], xmax[NDIM]; \
  } Result; \
  typedef const Result cResult; \
  typedef struct region { \
    int depth, next; \
    count isamples, cutcomp, xmajor; \
    real fmajor, fminor, vol; \
    Bounds bounds[NDIM]; \
    Result result[NCOMP]; \
  } Region

#define RegionPtr(n) (&((Region *)t->voidregion)[n])

#define CHUNKSIZE 4096

#define AllocRegions(t) \
  MemAlloc((t)->voidregion, ((t)->size = CHUNKSIZE)*sizeof(Region))

#define EnlargeRegions(t, n) if( (t)->nregions + n > (t)->size ) \
  ReAlloc((t)->voidregion, ((t)->size += CHUNKSIZE)*sizeof(Region))

#define SAMPLERDEFS \
  TYPEDEFREGION; \
  Region *region = RegionPtr(iregion); \
  cBounds *b = region->bounds; \
  Result *r = region->result; \
  cSamples *samples = &t->samples[region->isamples]; \
  real *x = samples->x, *f = samples->f; \
  cnumber n = samples->n

