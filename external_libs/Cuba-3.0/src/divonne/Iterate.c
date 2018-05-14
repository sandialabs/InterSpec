/*
	Iterate.c
		recursion over regions
		this file is part of Divonne
		last modified 21 Dec 11 th
*/


static void Iterate(This *t, count iregion, cint depth, cint isamples,
  Totals *totals)
{
  TYPEDEFREGION;
  Region *parent, *region;
  typedef struct {
    real avg, err, spread, spreadsq;
  } Corr;
  Corr corr[NCOMP];
  count ireg, mreg = iregion;
  count comp, maxsplit;
  int last, idest, isrc;

  region = RegionPtr(iregion);
  region->depth = depth;
  region->next = -iregion - 1;
  if( isamples < 0 ) Split(t, iregion);
  else {
    region->isamples = isamples;
    ExploreSerial(t, iregion);
  }

  ireg = iregion + RegionPtr(iregion)->next;

  do {
    region = RegionPtr(ireg);
    if( region->depth > 0 ) {
      --region->depth;
FORK_ONLY(more:)
      ireg = Explore(t, ireg);
      if( ireg == -1 ) return;
      region = RegionPtr(ireg);
    }
    if( region->depth < 0 ) mreg = IMax(mreg, ireg);
    ireg += region->next;
  } while( ireg > 0 );

  FORK_ONLY(if( t->running ) goto more;)

  maxsplit = 1;
  for( ireg = mreg; ireg >= iregion; --ireg ) {
    parent = RegionPtr(ireg);
    maxsplit -= NegQ(parent->depth);
    if( parent->depth < 0 ) {
      count xreg;
      struct {
        count from, to;
      } todo[maxsplit], *tdmax = todo, *td;
      count nsplit = 0;
      real norm;

      memset(corr, 0, sizeof corr);

      tdmax->from = ireg + parent->next;
      tdmax->to = tdmax->from - parent->depth;
      ++tdmax;
      for( td = todo; td < tdmax; ++td ) {
        for( xreg = td->from; xreg < td->to; ++xreg ) {
          Region *region = RegionPtr(xreg);
          if( region->depth < 0 ) {
            tdmax->from = xreg + region->next;
            tdmax->to = tdmax->from - region->depth;
            ++tdmax;
          }
          else {
            ++nsplit;
            for( comp = 0; comp < t->ncomp; ++comp ) {
              cResult *r = &region->result[comp];
              Corr *c = &corr[comp];
              c->avg += r->avg;
              c->err += r->err;
              c->spread += Sq(r->spread);
            }
          }
        }
      }

      norm = 1./nsplit--;
      for( comp = 0; comp < t->ncomp; ++comp ) {
        Result *p = &parent->result[comp];
        Corr *c = &corr[comp];
        creal diff = fabs(p->avg - c->avg)*norm;
        c->avg = diff*norm*nsplit;
        c->err = (c->err == 0) ? 1 : 1 + diff/c->err;
        c->spread = (c->spread == 0) ? 1 : 1 + diff/sqrt(c->spread);
      }

      for( td = todo; td < tdmax; ++td )
        for( xreg = td->from; xreg < td->to; ++xreg ) {
          Region *region = RegionPtr(xreg);
          if( region->depth >= 0 ) {
            cnumber neff = t->samples[region->isamples].neff;
            for( comp = 0; comp < t->ncomp; ++comp ) {
              Result *r = &region->result[comp];
              Corr *c = &corr[comp];
              if( r->err > 0 ) r->err = r->err*c->err + c->avg;
              r->spread = r->spread*c->spread + c->avg*neff;
              c->spreadsq += Sq(r->spread);
            }
          }
        }
    }
  }

  if( totals )
    for( comp = 0; comp < t->ncomp; ++comp )
      totals[comp].spreadsq += corr[comp].spreadsq;

  for( last = -1, idest = isrc = iregion; iregion <= mreg; ++iregion ) {
    Region *region = RegionPtr(iregion);
    cint cur = NegQ(region->depth);
    switch( cur - last ) {
    case -1:
      memmove(RegionPtr(idest), RegionPtr(isrc),
        (iregion - isrc)*sizeof(Region));
      idest += iregion - isrc;
      break;
    case 1:
      isrc = iregion;
    }
    last = cur;
  }

  memmove(RegionPtr(idest), RegionPtr(iregion),
    (t->nregions - iregion)*sizeof(Region));
  t->nregions += idest - iregion;
}

