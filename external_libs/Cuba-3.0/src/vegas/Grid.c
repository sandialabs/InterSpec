/*
	Grid.c
		utility functions for the Vegas grid
		this file is part of Vegas
		last modified 13 Dec 11 th
*/


static inline void GetGrid(cThis *t, Grid *grid)
{
  count bin, dim;
  unsigned const int slot = abs(t->gridno) - 1;

  if( t->gridno < 0 && slot < MAXGRIDS ) griddim_[slot] = 0;

  if( slot < MAXGRIDS && gridptr_[slot] ) {
    if( griddim_[slot] == t->ndim ) {
      XCopy(grid, gridptr_[slot]);
      return;
    }
    free(gridptr_[slot]);
    gridptr_[slot] = NULL;
  }

  for( bin = 0; bin < NBINS; ++bin )
    grid[0][bin] = (bin + 1)/(real)NBINS;
  for( dim = 1; dim < t->ndim; ++dim )
    Copy(&grid[dim], &grid[0], 1);
}

/*********************************************************************/

static inline void PutGrid(cThis *t, Grid *grid)
{
  unsigned const int slot = abs(t->gridno) - 1;

  if( slot < MAXGRIDS ) {
    if( gridptr_[slot] == NULL ) Alloc(gridptr_[slot], t->ndim);
    griddim_[slot] = t->ndim;
    XCopy(gridptr_[slot], grid);
  }
}

/*********************************************************************/

static void RefineGrid(cThis *t, Grid grid, Grid margsum)
{
  real avgperbin, thisbin, newcur, delta;
  Grid imp, newgrid;
  int bin, newbin;

  /* smooth the f^2 value stored for each bin */
  real prev = margsum[0];
  real cur = margsum[1];
  real norm = margsum[0] = .5*(prev + cur);
  for( bin = 1; bin < NBINS - 1; ++bin ) {
    creal s = prev + cur;
    prev = cur;
    cur = margsum[bin + 1];
    norm += margsum[bin] = (s + cur)/3.;
  }
  norm += margsum[NBINS - 1] = .5*(prev + cur);

  if( norm == 0 ) return;
  norm = 1/norm;

  /* compute the importance function for each bin */
  avgperbin = 0;
  for( bin = 0; bin < NBINS; ++bin ) {
    real impfun = 0;
    if( margsum[bin] > 0 ) {
      creal r = margsum[bin]*norm;
      avgperbin += impfun = pow((r - 1)/log(r), 1.5);
    }
    imp[bin] = impfun;
  }
  avgperbin /= NBINS;

  /* redefine the size of each bin */
  cur = newcur = 0;
  thisbin = 0;
  bin = -1;
  for( newbin = 0; newbin < NBINS - 1; ++newbin ) {
    while( thisbin < avgperbin ) {
      thisbin += imp[++bin];
      prev = cur;
      cur = grid[bin];
    }
    thisbin -= avgperbin;
    delta = (cur - prev)*thisbin;
    newgrid[newbin] = SHARPEDGES ?
      cur - delta/imp[bin] :
      (newcur = Max(newcur,
        cur - 2*delta/(imp[bin] + imp[IDim(bin - 1)])));
  }
  Copy(grid, newgrid, NBINS - 1);
  grid[NBINS - 1] = 1;
}

