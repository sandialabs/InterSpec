/*
	Sample.c
		most of what is related to sampling
		this file is part of Divonne
		last modified 17 Dec 11 th
*/


#define MARKMASK 0xfffffff
#define Marked(x) ((x) & ~MARKMASK)
#define Unmark(x) ((x) & MARKMASK)

#define EXTRAPOLATE_EPS (.25*t->border.lower)
/*#define EXTRAPOLATE_EPS 0x1p-26*/

/*********************************************************************/

static inline void SamplesIni(Samples *samples)
{
  samples->x = NULL;
}

/*********************************************************************/

static inline bool SamplesIniQ(cSamples *samples)
{
  return samples->x == NULL;
}

/*********************************************************************/

static inline void SamplesFree(cSamples *samples)
{
  free(samples->x);
}

/*********************************************************************/

static void SampleSobol(This *t, ccount iregion)
{
  SAMPLERDEFS;
  real avg[NCOMP], norm;
  number i;
  count dim, comp;

  for( i = 0; i < n; ++i ) {
    t->rng.getrandom(t, x);
    for( dim = 0; dim < t->ndim; ++x, ++dim )
      *x = b[dim].lower + *x*(b[dim].upper - b[dim].lower);
  }

  DoSample(t, n, samples->x, f);

  FCopy(avg, f);
  f += t->ncomp;
  for( i = 2; i < n; ++i )
    for( comp = 0; comp < t->ncomp; ++comp )
      avg[comp] += *f++;

  norm = region->vol/samples->neff;
  for( comp = 0; comp < t->ncomp; ++comp ) {
    r[comp].avg = norm*avg[comp];
    r[comp].err = 0;
  }
}

/*********************************************************************/

static void SampleKorobov(This *t, ccount iregion)
{
  SAMPLERDEFS;
  real *xlast = x + t->ndim, *flast = f + t->ncomp;
  real avg[NCOMP], norm;
  cnumber neff = samples->neff;
  number nextra = 0, i;
  real dist = 0;
  count dim, comp;

  for( i = 1; i < n; ++i ) {
    number c = i;
    for( dim = 0; dim < t->ndim; ++dim ) {
      creal dx = abs(2*c - neff)/(real)neff;
      *xlast++ = b[dim].lower + dx*(b[dim].upper - b[dim].lower);
      c = c*samples->coeff % neff;
    }
  }

  for( dim = 0; dim < t->ndim; ++dim ) {
    creal dx = (x[dim] = b[dim].upper) - t->border.upper;
    if( dx > 0 ) dist += Sq(dx);
  }

  if( dist > 0 ) {
    dist = sqrt(dist)/EXTRAPOLATE_EPS;
    for( dim = 0; dim < t->ndim; ++dim ) {
      real x2 = x[dim], dx = x2 - t->border.upper;
      if( dx > 0 ) {
        x[dim] = t->border.upper;
        x2 = t->border.upper - dx/dist;
      }
      xlast[dim] = x2;
    }
    nextra = 1;
  }

  DoSample(t, n + nextra, x, f);

  FCopy(avg, flast);
  flast += t->ncomp;
  for( i = 2; i < n; ++i )
    for( comp = 0; comp < t->ncomp; ++comp )
      avg[comp] += *flast++;

  if( nextra ) {
    for( comp = 0; comp < t->ncomp; ++comp )
      f[comp] += dist*(f[comp] - flast[comp]);
    for( dim = 0; dim < t->ndim; ++dim )
      x[dim] = b[dim].upper;
  }

  norm = region->vol/samples->neff;
  for( comp = 0; comp < t->ncomp; ++comp ) {
    r[comp].avg = norm*(avg[comp] + avg[comp] + f[comp]);
    r[comp].err = 0;
  }
}

/*********************************************************************/

#define IsSobol(k) NegQ(k)
#define IsRule(k, d) (k == 9 || k == 7 || (k == 11 && d == 3) || (k == 13 && d == 2))

/* The following coding is used for key1, key2, key3:
             0 = for key1, key2: use default,
                 for key3: do nothing,
             1 = for key3: split region again,
             7 = degree-7 cubature rule,
             9 = degree-9 cubature rule,
            11 = degree-11 cubature rule (only in 3 dims),
            13 = degree-13 cubature rule (only in 2 dims),
     -inf..-40 = absolute # of points, Sobol numbers,
       -39..-1 = multiplicator,        Sobol numbers,
         1..39 = multiplicator,        Korobov numbers,
       40..inf = absolute # of points, Korobov numbers.  */

static count SamplesLookup(This *t, Samples *samples, cint key,
  cnumber nwant, cnumber nmax, number nmin)
{
  number n;

  if( key == 13 && t->ndim == 2 ) {
    if( RuleIniQ(&t->rule13) ) Rule13Alloc(t);
    samples->rule = &t->rule13;
    samples->n = n = nmin = t->rule13.n;
    samples->sampler = SampleRule;
  }
  else if( key == 11 && t->ndim == 3 ) {
    if( RuleIniQ(&t->rule11) ) Rule11Alloc(t);
    samples->rule = &t->rule11;
    samples->n = n = nmin = t->rule11.n;
    samples->sampler = SampleRule;
  }
  else if( key == 9 ) {
    if( RuleIniQ(&t->rule9) ) Rule9Alloc(t);
    samples->rule = &t->rule9;
    samples->n = n = nmin = t->rule9.n;
    samples->sampler = SampleRule;
  }
  else if( key == 7 ) {
    if( RuleIniQ(&t->rule7) ) Rule7Alloc(t);
    samples->rule = &t->rule7;
    samples->n = n = nmin = t->rule7.n;
    samples->sampler = SampleRule;
  }
  else {
    n = Abs1(key);
    if( n < 40 ) n *= nwant;
    samples->sampler = (key < 0) ? SampleSobol :
      (n = n/2 + 1, SampleKorobov);
    samples->n = IMin(n, nmax);
  }

  samples->neff = samples->n;

  return IDim(n - nmax) | Marked(nmax - nmin);
}

/*********************************************************************/

static void SamplesAlloc(cThis *t, Samples *samples)
{
#define FIRST -INT_MAX
#define MarkLast(x) (x | Marked(INT_MAX))

#include "KorobovCoeff.c"

  number nx, nf;

  if( samples->sampler == SampleKorobov ) {
    enum { max = Elements(prime) - 2 };
    cint n = IMin(2*samples->n - 1, MAXPRIME);
    int i = Hash(n), p;
    count shift = 2 + NegQ(n - 1000);

    while( i = IMin(IDim(i), max),
           n > (p = prime[i + 1]) || n <= prime[i] ) {
      cint d = (n - Unmark(p)) >> ++shift;
      i += Min1(d);
    }

    samples->coeff = coeff[i][t->ndim - KOROBOV_MINDIM];
    samples->neff = p = Unmark(p);
    samples->n = p/2 + 1;
  }

  nx = t->ndim*(samples->n + 1);		/* need 1 for extrapolation */
  nf = t->ncomp*(samples->n + 1);

  Alloc(samples->x, nx + nf + t->ncomp + t->ncomp);
  samples->f = samples->x + nx;
}

/*********************************************************************/

static real Sample(This *t, creal *x0)
{
  real xtmp[2*NDIM], ftmp[2*NCOMP], *xlast = xtmp, f;
  real dist = 0;
  count dim, comp;
  number n = 1;

  for( dim = 0; dim < t->ndim; ++dim ) {
    creal x1 = *xlast++ = Min(Max(*x0++, 0.), 1.);
    real dx;
    if( (dx = x1 - t->border.lower) < 0 ||
        (dx = x1 - t->border.upper) > 0 ) dist += Sq(dx);
  }

  if( dist > 0 ) {
    dist = sqrt(dist)/EXTRAPOLATE_EPS;
    for( dim = 0; dim < t->ndim; ++dim ) {
      real x2 = xtmp[dim], dx, b;
      if( (dx = x2 - (b = t->border.lower)) < 0 ||
          (dx = x2 - (b = t->border.upper)) > 0 ) {
        xtmp[dim] = b;
        x2 = b - dx/dist;
      }
      *xlast++ = x2;
    }
    n = 2;
  }

  DoSample(t, n, xtmp, ftmp);

  comp = Untag(t->selectedcomp);
  f = ftmp[comp];
  if( n > 1 ) f += dist*(f - ftmp[comp + t->ncomp]);

  return Sign(t->selectedcomp)*f;
}

