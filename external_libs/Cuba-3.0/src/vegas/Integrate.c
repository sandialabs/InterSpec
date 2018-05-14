/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 17 Apr 12 th
*/


static int Integrate(This *t, real *integral, real *error, real *prob)
{
  bin_t *bins;
  count dim, comp;
  int fail;
  struct {
    count niter;
    number nsamples, neval;
    Cumulants cumul[NCOMP];
    Grid grid[NDIM];
  } state;
  int statemsg = VERBOSE;
  struct stat st;

  if( VERBOSE > 1 ) {
    char s[512];
    sprintf(s, "Vegas input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  epsrel " REAL "\n  epsabs " REAL "\n"
      "  flags %d\n  seed %d\n"
      "  mineval " NUMBER "\n  maxeval " NUMBER "\n"
      "  nstart " NUMBER "\n  nincrease " NUMBER "\n"
      "  nbatch " NUMBER "\n  gridno %d\n"
      "  statefile \"%s\"",
      t->ndim, t->ncomp,
      t->epsrel, t->epsabs,
      t->flags, t->seed,
      t->mineval, t->maxeval,
      t->nstart, t->nincrease, t->nbatch,
      t->gridno, t->statefile);
    Print(s);
  }

  if( BadComponent(t) ) return -2;
  if( BadDimension(t) ) return -1;

  FrameAlloc(t, ShmRm(t));
  ForkCores(t);
  Alloc(bins, t->nbatch*t->ndim);

  if( (fail = setjmp(t->abort)) ) goto abort;

  IniRandom(t);

  if( t->statefile && *t->statefile == 0 ) t->statefile = NULL;

  if( t->statefile &&
      stat(t->statefile, &st) == 0 &&
      st.st_size == sizeof state && (st.st_mode & 0400) ) {
    cint h = open(t->statefile, O_RDONLY);
    read(h, &state, sizeof state);
    close(h);
    t->rng.skiprandom(t, t->neval = state.neval);

    if( VERBOSE ) {
      char s[256];
      sprintf(s, "\nRestoring state from %s.", t->statefile);
      Print(s);
    }
  }
  else {
    state.niter = 0;
    state.nsamples = t->nstart;
    Zap(state.cumul);
    GetGrid(t, state.grid);
  }

  /* main iteration loop */

  for( ; ; ) {
    number nsamples = state.nsamples;
    creal jacobian = 1./nsamples;
    Grid margsum[NCOMP][NDIM];

    Zap(margsum);

    for( ; nsamples > 0; nsamples -= t->nbatch ) {
      cnumber n = IMin(t->nbatch, nsamples);
      real *w = t->frame;
      real *x = w + n;
      real *f = x + n*t->ndim;
      real *lastf = f + n*t->ncomp;
      bin_t *bin = bins;

      while( x < f ) {
        real weight = jacobian;

        t->rng.getrandom(t, x);

        for( dim = 0; dim < t->ndim; ++dim ) {
          creal pos = *x*NBINS;
          ccount ipos = (count)pos;
          creal prev = (ipos == 0) ? 0 : state.grid[dim][ipos - 1];
          creal diff = state.grid[dim][ipos] - prev; 
          *x++ = prev + (pos - ipos)*diff;
          *bin++ = ipos;
          weight *= diff*NBINS;
        }

        *w++ = weight;
      }

      DoSample(t, n, w, f, t->frame, state.niter + 1);

      bin = bins;
      w = t->frame;

      while( f < lastf ) {
        creal weight = *w++;

        for( comp = 0; comp < t->ncomp; ++comp ) {
          real wfun = weight*(*f++);
          if( wfun ) {
            Cumulants *c = &state.cumul[comp];
            Grid *m = margsum[comp];

            c->sum += wfun;
            c->sqsum += wfun *= wfun;
            for( dim = 0; dim < t->ndim; ++dim )
              m[dim][bin[dim]] += wfun;
          }
        }

        bin += t->ndim;
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( comp = 0; comp < t->ncomp; ++comp ) {
      Cumulants *c = &state.cumul[comp];
      real avg, sigsq;
      real w = Weight(c->sum, c->sqsum, state.nsamples);

      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > MaxErr(c->avg));

      if( state.niter == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration " COUNT ":  " NUMBER " integrand evaluations so far",
        state.niter + 1, t->neval);

      for( comp = 0; comp < t->ncomp; ++comp ) {
        cCumulants *c = &state.cumul[comp];
        p += sprintf(p, "\n[" COUNT "] "
          REAL " +- " REAL "  \tchisq " REAL " (" COUNT " df)",
          comp + 1, c->avg, c->err, c->chisq, state.niter);
      }

      Print(s);
    }

    if( fail == 0 && t->neval >= t->mineval ) {
      if( t->statefile && KEEPFILE == 0 ) unlink(t->statefile);
      break;
    }

    if( t->neval >= t->maxeval && t->statefile == NULL ) break;

    if( t->ncomp == 1 )
      for( dim = 0; dim < t->ndim; ++dim )
        RefineGrid(t, state.grid[dim], margsum[0][dim]);
    else {
      for( dim = 0; dim < t->ndim; ++dim ) {
        Grid wmargsum;
        Zap(wmargsum);
        for( comp = 0; comp < t->ncomp; ++comp ) {
          real w = state.cumul[comp].avg;
          if( w != 0 ) {
            creal *m = margsum[comp][dim];
            count bin;
            w = 1/Sq(w);
            for( bin = 0; bin < NBINS; ++bin )
              wmargsum[bin] += w*m[bin];
          }
        }
        RefineGrid(t, state.grid[dim], wmargsum);
      }
    }

    ++state.niter;
    state.nsamples += t->nincrease;

    if( t->statefile ) {
      cint h = creat(t->statefile, 0666);
      if( h != -1 ) {
        state.neval = t->neval;
        write(h, &state, sizeof state);
        close(h);

        if( statemsg ) {
          char s[256];
          sprintf(s, "\nSaving state to %s.", t->statefile);
          Print(s);
          statemsg = false;
        }
      }
      if( t->neval >= t->maxeval ) break;
    }
  }

  for( comp = 0; comp < t->ncomp; ++comp ) {
    cCumulants *c = &state.cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, state.niter);
  }

abort:
  PutGrid(t, state.grid);
  free(bins);
  WaitCores(t);
  FrameFree(t);

  return fail;
}

