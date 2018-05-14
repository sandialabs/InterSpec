/*
	Split.c
		determine optimal cuts for splitting a region
		this file is part of Divonne
		last modified 18 Dec 11 th
*/


#define BNDTOL .05
#define FRACT .5
#define BIG 1e10
#define SINGTOL 1e-4

#define LHSTOL .1
#define GAMMATOL .1

/* the next four macros must be in sync with the typedef of Bounds! */
#define Lower(d) (2*d)
#define Upper(d) (2*d + 1)
#define Dim(i) ((i) >> 1)
#define SignedDelta(i) (2*(i & 1) - 1)*delta[i]

typedef struct {
  count i;
  real save, delta;
  real f, df, fold;
  real lhs, row, sol;
} Cut;

/*********************************************************************/

static inline real Div(creal a, creal b)
{
  return (b != 0 && fabs(b) < BIG*fabs(a)) ? a/b : a;
}

/*********************************************************************/

static void SomeCut(This *t, Cut *cut, Bounds *b)
{
  count dim, maxdim;
  static count nextdim = 0;
  real xmid[NDIM], ymid, maxdev;

  for( dim = 0; dim < t->ndim; ++dim )
    xmid[dim] = .5*(b[dim].upper + b[dim].lower);
  ymid = Sample(t, xmid);

  maxdev = 0;
  maxdim = 0;
  for( dim = 0; dim < t->ndim; ++dim ) {
    real ylower, yupper, dev;
    creal x = xmid[dim];
    xmid[dim] = b[dim].lower;
    ylower = Sample(t, xmid);
    xmid[dim] = b[dim].upper;
    yupper = Sample(t, xmid);
    xmid[dim] = x;

    dev = fabs(ymid - .5*(ylower + yupper));
    if( dev >= maxdev ) {
      maxdev = dev;
      maxdim = dim;
    }
  }

  if( maxdev > 0 ) nextdim = 0;
  else maxdim = nextdim++ % t->ndim;

  cut->i = Upper(maxdim);
  cut->save = b[maxdim].upper;
  b[maxdim].upper = xmid[maxdim];
}

/*********************************************************************/

static inline real Volume(cThis *t, creal *delta)
{
  real vol = 1;
  count dim;
  for( dim = 0; dim < t->ndim; ++dim )
    vol *= delta[Lower(dim)] + delta[Upper(dim)];
  return vol;
}

/*********************************************************************/

static inline real SetupEqs(Cut *cut, ccount ncuts, real f)
{
  real sqsum = 0;
  Cut *c = &cut[ncuts];
  while( --c >= cut ) {
    sqsum += Sq(c->lhs = f - c->f);
    f = c->f;
  }
  return sqsum;
}

/*********************************************************************/

static inline void SolveEqs(Cut *cut, count ncuts,
  creal *delta, creal diff)
{
  real last = 0;
  real r = 1;
  Cut *c;

  for( c = cut; ; ++c ) {
    ccount dim = Dim(c->i);
    c->row = r -=
      Div(diff, (delta[Lower(dim)] + delta[Upper(dim)])*c->df);
    if( --ncuts == 0 ) break;
    last += r*c->lhs;
  }

  last = Div(c->lhs - last, r);

  for( ; c >= cut; last += (--c)->lhs ) {
    creal delmin = -(c->delta = delta[c->i]);
    creal delmax = FRACT*(delmin + c->save);
    c->sol = Div(last, c->df);
    if( c->sol > delmax ) c->sol = .75*delmax;
    if( c->sol < delmin ) c->sol = .75*delmin;
  }
}

/*********************************************************************/

static count FindCuts(This *t, Cut *cut, Bounds *bounds, creal vol,
  real *xmajor, creal fmajor, creal fdiff)
{
  cint sign = (fdiff < 0) ? -1 : 1;

  count ncuts = 0, icut;
  real delta[2*NDIM];
  real gamma, fgamma, lhssq;
  count dim, div;

  for( dim = 0; dim < t->ndim; ++dim ) {
    cBounds *b = &bounds[dim];
    creal xsave = xmajor[dim];
    real dist = b->upper - xsave;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      Cut *c = &cut[ncuts++];
      c->i = Upper(dim);
      c->save = dist;
      xmajor[dim] += dist *= FRACT;
      c->f = Sample(t, xmajor);
      xmajor[dim] = xsave;
    }
    delta[Upper(dim)] = dist;
  }

  for( dim = 0; dim < t->ndim; ++dim ) {
    cBounds *b = &bounds[dim];
    creal xsave = xmajor[dim];
    real dist = xsave - b->lower;
    if( dist >= BNDTOL*(b->upper - b->lower) ) {
      Cut *c = &cut[ncuts++];
      c->i = Lower(dim);
      c->save = dist;
      xmajor[dim] -= dist *= FRACT;
      c->f = Sample(t, xmajor);
      xmajor[dim] = xsave;
    }
    delta[Lower(dim)] = dist;
  }

  if( ncuts == 0 ) {
    SomeCut(t, cut, bounds);
    return 1;
  }

  for( ; ; ) {
    real mindiff = INFTY;
    Cut *mincut = cut;

    for( icut = 0; icut < ncuts; ++icut ) {
      Cut *c = &cut[icut];
      creal diff = fabs(fmajor - c->f);
      if( diff <= mindiff ) {
        mindiff = diff;
        mincut = c;
      }
    }

    gamma = Volume(t, delta)/vol;
    fgamma = fmajor + (gamma - 1)*fdiff;

    if( sign*(mincut->f - fgamma) < 0 ) break;

    if( --ncuts == 0 ) {
      SomeCut(t, cut, bounds);
      return 1;
    }

    delta[mincut->i] = mincut->save;
    memmove(mincut, mincut + 1, (char *)&cut[ncuts] - (char *)mincut);
  }

  for( icut = 0; icut < ncuts; ++icut ) {
    Cut *c = &cut[icut];
    c->fold = c->f;
    c->df = (c->f - fmajor)/delta[c->i];
  }

  lhssq = SetupEqs(cut, ncuts, fgamma);

repeat:
  SolveEqs(cut, ncuts, delta, gamma*fdiff);

  for( div = 1; div <= 16; div *= 4 ) {
    real gammanew, lhssqnew;

    for( icut = 0; icut < ncuts; ++icut ) {
      Cut *c = &cut[icut];
      real *x = &xmajor[Dim(c->i)];
      creal xsave = *x;
      delta[c->i] = c->delta + c->sol/div;
      *x += SignedDelta(c->i);
      c->f = Sample(t, xmajor);
      *x = xsave;
    }

    gammanew = Volume(t, delta)/vol;
    fgamma = fmajor + (gammanew - 1)*fdiff;
    lhssqnew = SetupEqs(cut, ncuts, fgamma);

    if( lhssqnew <= lhssq ) {
      real fmax;

      if( fabs(gammanew - gamma) < GAMMATOL*gamma ) break;
      gamma = gammanew;

      fmax = fabs(fgamma);
      for( icut = 0; icut < ncuts; ++icut ) {
        Cut *c = &cut[icut];
        creal dfmin = SINGTOL*c->df;
        creal sol = c->sol/div;
        real df = c->f - c->fold;
        df = (fabs(sol) < BIG*fabs(df)) ? df/sol : 1;
        c->df = (fabs(df) < fabs(dfmin)) ? dfmin : df;
        fmax = Max(fmax, fabs(c->f));
        c->fold = c->f;
      }

      if( lhssqnew < Sq((1 + fmax)*LHSTOL) ) break;
      lhssq = lhssqnew;
      goto repeat;
    }
  }

  for( icut = 0; icut < ncuts; ++icut ) {
    Cut *c = &cut[icut];
    real *b = (real *)bounds + c->i;
    c->save = *b;
    *b = xmajor[Dim(c->i)] + SignedDelta(c->i);
  }

  return ncuts;
}

/*********************************************************************/

static void Split(This *t, ccount iregion)
{
  TYPEDEFREGION;
  Region *region = RegionPtr(iregion);
  Cut cut[2*NDIM], *c;
  count ncuts, succ;
  int depth;
  real *b;

  t->selectedcomp = region->cutcomp;
  t->neval_cut -= t->neval;
  ncuts = FindCuts(t, cut, region->bounds, region->vol,
    (real *)region->result + region->xmajor, region->fmajor,
    region->fmajor - region->fminor);
  t->neval_cut += t->neval;

  depth = region->depth - ncuts;

  EnlargeRegions(t, ++ncuts);
  region = RegionPtr(iregion);
  region->depth = -ncuts;
  succ = iregion + region->next;
  region->next = t->nregions - iregion;
  b = (real *)region->bounds;

  region = RegionPtr(t->nregions);
  XCopy(region->bounds, b);
  region->depth = IDim(depth) + 1;
  region->next = 1;
  region->isamples = 0;
  for( c = cut; --ncuts; ++c ) {
    ccount ci = c->i;
    creal tmp = b[ci ^ 1];
    b[ci ^ 1] = b[ci];
    b[ci] = c->save;
    region = RegionPtr(++t->nregions);
    XCopy(region->bounds, b);
    region->depth = IDim(depth) + 1;
    region->next = 1;
    region->isamples = 0;
    ++depth;
    b[ci ^ 1] = tmp;
  }
  region->next = succ - t->nregions++;
}

