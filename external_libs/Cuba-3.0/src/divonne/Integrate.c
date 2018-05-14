/*
	Integrate.c
		partition the integration region until each region
		has approximately equal spread = 1/2 vol (max - min),
		then do a main integration over all regions
		this file is part of Divonne
		last modified 3 Jul 12 th
*/


#define INIDEPTH 3
#define DEPTH 5
#define POSTDEPTH 15

/*********************************************************************/

static int Integrate(This *t, real *integral, real *error, real *prob)
{
  TYPEDEFREGION;

  Totals totals[NCOMP];
  real nneed;
  count dim, comp, iter, pass = 0, err, iregion;
  number nwant, nmin = INT_MAX;
  ML_ONLY(number neff;)
  int fail;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Divonne input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  key1 %d\n  key2 %d\n  key3 %d\n  maxpass " COUNT "\n"
      "  border " REAL "\n  maxchisq " REAL "\n  mindeviation " REAL "\n"
      "  ngiven " NUMBER "\n  nextra " NUMBER,
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->key1, t->key2, t->key3, t->maxpass,
      t->border.lower, t->maxchisq, t->mindeviation,
      t->ngiven, t->nextra);
    Print(s);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t, t->key1) ||
      BadDimension(t, t->key2) ||   
      ((t->key3 & -2) && BadDimension(t, t->key3)) ) return -1;

  FORK_ONLY(t->nframe = 0;)
  t->neval_opt = t->neval_cut = 0;

  ForkCores(t);

  AllocGiven(t);

  AllocRegions(t);
  for( dim = 0; dim < t->ndim; ++dim ) {
    Bounds *b = &RegionPtr(0)->bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }
  t->nregions = 1;

  RuleIni(&t->rule7);
  RuleIni(&t->rule9);
  RuleIni(&t->rule11);
  RuleIni(&t->rule13);
  SamplesIni(&t->samples[0]);
  SamplesIni(&t->samples[1]);
  SamplesIni(&t->samples[2]);

  if( (fail = setjmp(t->abort)) ) goto abort;

  t->epsabs = Max(t->epsabs, NOTZERO);

  /* Step 1: partition the integration region */

  if( VERBOSE ) Print("Partitioning phase:");

  if( IsSobol(t->key1) || IsSobol(t->key2) || IsSobol(t->key3) )
    IniRandom(t);

  SamplesLookup(t, &t->samples[0], t->key1,
    (number)47, (number)INT_MAX, (number)0);
  SamplesAlloc(t, &t->samples[0]);

  t->totals = totals;
  Zap(totals);
  t->phase = 1;

  Iterate(t, 0, INIDEPTH, 0, NULL);

  for( iter = 1; ; ++iter ) {
    Totals *maxtot;
    count valid;

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Totals *tot = &totals[comp];
      tot->avg = tot->spreadsq = 0;
      tot->spread = tot->secondspread = -INFTY;
    }

    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];
        tot->avg += r->avg;
        tot->spreadsq += Sq(r->spread);
        if( r->spread > tot->spread ) {
          tot->secondspread = tot->spread;
          tot->spread = r->spread;
          tot->iregion = iregion;
        }
        else if( r->spread > tot->secondspread )
          tot->secondspread = r->spread;
      }
    }

    maxtot = totals;
    valid = 0;
    for( comp = 0; comp < t->ncomp; ++comp ) {
      Totals *tot = &totals[comp];
      integral[comp] = tot->avg;
      valid += tot->avg == tot->avg;
      if( tot->spreadsq > maxtot->spreadsq ) maxtot = tot;
      tot->spread = sqrt(tot->spreadsq);
      error[comp] = tot->spread/t->samples[0].neff;
    }

    if( VERBOSE ) {
      char s[128 + 64*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT " (pass " COUNT "):  " COUNT " regions\n"
        NUMBER7 " integrand evaluations so far,\n"
        NUMBER7 " in optimizing regions,\n"
        NUMBER7 " in finding cuts",
        iter, pass, t->nregions, t->neval, t->neval_opt, t->neval_cut);

      for( comp = 0; comp < t->ncomp; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL,
          comp + 1, integral[comp], error[comp]);

      Print(s);
    }

    if( valid == 0 ) goto abort;	/* all NaNs */

    if( t->neval > t->maxeval ) break;

    nneed = maxtot->spread/MaxErr(maxtot->avg);
    if( nneed < MAXPRIME ) {
      cnumber n = t->neval + t->nregions*(number)ceil(nneed);
      if( n < nmin ) {
        nmin = n;
        pass = 0;
      }
      else if( ++pass > t->maxpass && n >= t->mineval ) break;
    }

    Iterate(t, maxtot->iregion, DEPTH, -1, NULL);
  }

  /* Step 2: do a "full" integration on each region */

/* nneed = t->samples[0].neff + 1; */
  nneed = 2*t->samples[0].neff;
  for( comp = 0; comp < t->ncomp; ++comp ) {
    Totals *tot = &totals[comp];
    creal maxerr = MaxErr(tot->avg);
    tot->nneed = tot->spread/maxerr;
    nneed = Max(nneed, tot->nneed);
    tot->maxerrsq = Sq(maxerr);
    tot->mindevsq = tot->maxerrsq*Sq(t->mindeviation);
  }
  nwant = (number)Min(ceil(nneed), MARKMASK/40.);

  err = SamplesLookup(t, &t->samples[1], t->key2, nwant,
    (t->maxeval - t->neval)/t->nregions + 1, t->samples[0].n + 1);

  /* the number of points needed to reach the desired accuracy */
  fail = Unmark(err)*t->nregions;

  if( Marked(err) ) {
    if( VERBOSE ) Print("\nNot enough samples left for main integration.");
    for( comp = 0; comp < t->ncomp; ++comp )
      prob[comp] = -999;
    ML_ONLY(neff = t->samples[0].neff;)
  }
  else {
    bool can_adjust = (t->key3 == 1 && t->samples[1].sampler != SampleRule &&
      (t->key2 < 0 || t->samples[1].neff < MAXPRIME));
    count df, nlimit;

    SamplesAlloc(t, &t->samples[1]);

    if( VERBOSE ) {
      char s[128];
      sprintf(s, "\nMain integration on " COUNT
        " regions with " NUMBER " samples per region.",
        t->nregions, t->samples[1].neff);
      Print(s);
    }

    FClear(integral);
    FClear(error);
    FClear(prob);

    nlimit = t->maxeval - t->nregions*t->samples[1].n;
    df = 0;

#define CopyPhaseResults(f) \
  for( comp = 0; comp < t->ncomp; ++comp ) { \
    PhaseResult *p = &totals[comp].phase[f]; \
    cResult *r = &region->result[comp]; \
    p->avg = r->avg; \
    p->err = r->err; \
  }

#define Var2(f, res) Sq((res)->err ? (res)->err : r->spread/t->samples[f].neff)
#define Var(f) Var2(f, &tot->phase[f])

    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region;
      char s[64*NDIM + 256*NCOMP], *p = s;
      int todo;

refine:
      region = RegionPtr(iregion);
      CopyPhaseResults(0);
      t->phase = 2;
      region->isamples = 1;
      t->samples[1].sampler(t, iregion);
      CopyPhaseResults(1);

      if( can_adjust )
        for( comp = 0; comp < t->ncomp; ++comp )
          totals[comp].spreadsq -= Sq(region->result[comp].spread);

      nlimit += t->samples[1].n;
      todo = 0;

      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];

        if( t->neval < nlimit ) {
          creal avg2 = tot->phase[1].avg;
          creal diffsq = Sq(avg2 - tot->phase[0].avg);

          if( r->err*tot->nneed > r->spread ||
              diffsq > Max(t->maxchisq*(Var(0) + Var(1)), EPS*Sq(avg2)) ) {
            if( t->key3 && diffsq > tot->mindevsq ) {
              if( t->key3 == 1 ) {
                if( VERBOSE > 2 ) Print("\nSplit");
                t->phase = 1;
                Iterate(t, iregion, POSTDEPTH, 1, totals);

                if( can_adjust ) {
                  cnumber nnew = (tot->spreadsq/Sq(MARKMASK) > tot->maxerrsq) ?
                    MARKMASK :
                    (number)ceil(sqrt(tot->spreadsq/tot->maxerrsq));
                  if( nnew > nwant + nwant/64 ) {
                    ccount err = SamplesLookup(t, &t->samples[1], t->key2, nnew,
                      (t->maxeval - t->neval)/t->nregions + 1, t->samples[1].n);
                    fail += Unmark(err)*t->nregions;
                    nwant = nnew;
                    SamplesFree(&t->samples[1]);
                    SamplesAlloc(t, &t->samples[1]);

                    if( t->key2 > 0 && t->samples[1].neff >= MAXPRIME )
                      can_adjust = false;

                    if( VERBOSE > 2 ) {
                      char s[128];
                      sprintf(s, "Sampling remaining " COUNT
                        " regions with " NUMBER " points per region.",
                        t->nregions, t->samples[1].neff);
                      Print(s);
                    }
                  }
                }
                goto refine;
              }
              todo |= 3;
            }
            todo |= 1;
          }
        }
      }

      if( can_adjust ) {
        for( comp = 0; comp < t->ncomp; ++comp )
          totals[comp].maxerrsq -=
            Sq(region->result[comp].spread/t->samples[1].neff);
      }

      switch( todo ) {
      case 1:	/* get spread right */
        region->isamples = 1;
        ExploreSerial(t, iregion);
        break;

      case 3:	/* sample region again with more points */
        if( SamplesIniQ(&t->samples[2]) ) {
          SamplesLookup(t, &t->samples[2], t->key3,
            nwant, (number)INT_MAX, (number)0);
          SamplesAlloc(t, &t->samples[2]);
        }
        t->phase = 3;
        region->isamples = 2;
        t->samples[2].sampler(t, iregion);
        ExploreSerial(t, iregion);
        ++region->depth;	/* misused for df here */
        ++df;
      }

      if( VERBOSE > 2 ) {
        for( dim = 0; dim < t->ndim; ++dim ) {
          cBounds *b = &region->bounds[dim];
          p += sprintf(p,
            (dim == 0) ? "\nRegion (" REALF ") - (" REALF ")" :
                         "\n       (" REALF ") - (" REALF ")",
            b->lower, b->upper);
        }
      }

      for( comp = 0; comp < t->ncomp; ++comp ) {
        Result *r = &region->result[comp];
        Totals *tot = &totals[comp];

        creal x1 = tot->phase[0].avg;
        creal v1 = Var(0);
        creal x2 = tot->phase[1].avg;
        creal v2 = Var(1);
        creal r2 = v1 ? v2/v1 :
          Sq(t->samples[1].neff/(real)t->samples[0].neff);

        real norm = 1 + r2;
        real avg = x2 + r2*x1;
        real sigsq = v2;
        real chisq = Sq(x2 - x1);
        real chiden = v1 + v2;

        if( todo == 3 ) {
          creal x3 = r->avg;
          creal v3 = Var2(2, r);
          creal r3 = v2 ? v3/v2 :
            Sq(t->samples[2].neff/(real)t->samples[1].neff);

          norm = 1 + r3*norm;
          avg = x3 + r3*avg;
          sigsq = v3;
          chisq = v1*Sq(x3 - x2) + v2*Sq(x3 - x1) + v3*chisq;
          chiden = v1*v2 + v3*chiden;
        }

        avg = LAST ? r->avg : (sigsq *= norm = 1/norm, avg*norm);
        if( chisq > EPS ) chisq /= Max(chiden, NOTZERO);

        if( VERBOSE > 2 ) {
#define Out2(f, res) (res)->avg, r->spread/t->samples[f].neff, (res)->err
#define Out(f) Out2(f, &tot->phase[f])
          p += sprintf(p, "\n[" COUNT "] "
            REAL " +- " REAL "(" REAL ")\n    "
            REAL " +- " REAL "(" REAL ")", comp + 1, Out(0), Out(1));
          if( todo == 3 ) p += sprintf(p, "\n    "
            REAL " +- " REAL "(" REAL ")", Out2(2, r));
          p += sprintf(p, "  \tchisq " REAL, chisq);
        }

        integral[comp] += avg;
        error[comp] += sigsq;
        prob[comp] += chisq;

        r->avg = avg;
        r->spread = sqrt(sigsq);
        r->chisq = chisq;
      }

      if( VERBOSE > 2 ) Print(s);
    }

    for( comp = 0; comp < t->ncomp; ++comp )
      error[comp] = sqrt(error[comp]);

    df += t->nregions;

    if( VERBOSE > 2 ) {
      char s[16 + 128*NCOMP], *p = s;

      p += sprintf(p, "\nTotals:");

      for( comp = 0; comp < t->ncomp; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, integral[comp], error[comp], prob[comp], df);

      Print(s);
    }

    for( comp = 0; comp < t->ncomp; ++comp )
      prob[comp] = ChiSquare(prob[comp], df);

    ML_ONLY(neff = 1;)
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", t->nregions);
    for( iregion = 0; iregion < t->nregions; ++iregion ) {
      Region *region = RegionPtr(iregion);
      cBounds *b = region->bounds;
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < t->ndim; ++dim ) {
        lower[dim] = b[dim].lower;
        upper[dim] = b[dim].upper;
      }

      MLPutFunction(stdlink, "Cuba`Divonne`region", 4);

      MLPutRealList(stdlink, lower, t->ndim);
      MLPutRealList(stdlink, upper, t->ndim);

      MLPutFunction(stdlink, "List", t->ncomp);
      for( comp = 0; comp < t->ncomp; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->spread/neff, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->depth + 1);  /* misused for df */
    }
  }
#endif

abort:
  SamplesFree(&t->samples[2]);
  SamplesFree(&t->samples[1]);
  SamplesFree(&t->samples[0]);
  RuleFree(&t->rule13);
  RuleFree(&t->rule11);
  RuleFree(&t->rule9);
  RuleFree(&t->rule7);

  free(t->voidregion);
  free(t->xgiven);
  WaitCores(t);
  FORK_ONLY(FrameFree(t, ShmRm(t));)

  return fail;
}

