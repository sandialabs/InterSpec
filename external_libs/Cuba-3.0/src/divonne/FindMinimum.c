/*
	FindMinimum.c
		find minimum (maximum) of hyperrectangular region
		this file is part of Divonne
		last modified 22 Dec 11 th
*/


#define EPS 0x1p-52
#define RTEPS 0x1p-26
#define QEPS 0x1p-13

#define DELTA 0x1p-16
#define RTDELTA 0x1p-8
#define QDELTA 0x1p-4

/*
#define DELTA 1e-5
#define RTDELTA 3.1622776601683791e-3
#define QDELTA 5.6234132519034912e-2
*/

#define SUFTOL 8*QEPS*QDELTA
#define FTOL 5e-2
#define GTOL 1e-2

#define Hessian(i, j) hessian[(i)*t->ndim + j]

typedef struct { real dx, f; } Point;

/*********************************************************************/

static inline real Dot(ccount n, creal *a, creal *b)
{
  real sum = 0;
  count i;
  for( i = 0; i < n; ++i ) sum += a[i]*b[i];
  return sum;
}

/*********************************************************************/

static inline real Length(ccount n, creal *vec)
{
  return sqrt(Dot(n, vec, vec));
}

/*********************************************************************/

static inline void LinearSolve(cThis *t, ccount n, creal *hessian,
  creal *grad, real *p)
{
  int i, j;
  real dir;

  for( i = 0; i < n; ++i ) {
    dir = -grad[i];
    for( j = 0; j < i; ++j )
      dir -= Hessian(i, j)*p[j];
    p[i] = dir;
  }

  while( --i >= 0 ) {
    if( Hessian(i, i) <= 0 ) return;
    dir = p[i]/Hessian(i, i);
    for( j = i + 1; j < n; ++j )
      dir -= Hessian(j, i)*p[j];
    p[i] = dir;
  }
}

/*********************************************************************/

static void RenormalizeCholesky(cThis *t, ccount n, real *hessian,
  real *z, real alpha)
{
  count i, j;

  for( i = 0; i < n; ++i ) {
    creal dir = z[i];
    real beta = alpha*dir;
    real gamma = Hessian(i, i);
    real gammanew = Hessian(i, i) += beta*dir;

    if( i + 1 >= n || gammanew < 0 ||
        (gammanew < 1 && gamma > DBL_MAX*gammanew) ) return;

    gamma /= gammanew;
    beta /= gammanew;
    alpha *= gamma;

    if( gamma < .25 ) {
      for( j = i + 1; j < n; ++j ) {
        real delta = beta*z[j];
        z[j] -= dir*Hessian(j, i);
        Hessian(j, i) = Hessian(j, i)*gamma + delta;
      }
    }
    else {
      for( j = i + 1; j < n; ++j ) {
        z[j] -= dir*Hessian(j, i);
        Hessian(j, i) += beta*z[j];
      }
    }
  }
}

/*********************************************************************/

static void UpdateCholesky(cThis *t, ccount n, real *hessian,
  real *z, real *p)
{
  int i, j;
  real gamma = 0;

  for( i = 0; i < n; ++i ) {
    real dir = z[i];
    for( j = 0; j < i; ++j )
      dir -= Hessian(i, j)*p[j];
    p[i] = dir;
    gamma += Sq(dir)/Hessian(i, i);
  }
  gamma = Max(fabs(1 - gamma), EPS);

  while( --i >= 0 ) {
    creal dir = z[i] = p[i];
    real beta = dir/Hessian(i, i);
    creal gammanew = gamma + dir*beta;
    Hessian(i, i) *= gamma/gammanew;
    beta /= gamma;
    gamma = gammanew;
    for( j = i + 1; j < n; ++j ) {
      creal delta = beta*z[j];
      z[j] += dir*Hessian(j, i);
      Hessian(j, i) -= delta;
    }
  }
}

/*********************************************************************/

static inline void BFGS(cThis *t, ccount n, real *hessian,
  creal *gnew, creal *g, real *p, creal dx)
{
  real y[NDIM], c;
  count i, j;

  for( i = 0; i < n; ++i )
    y[i] = gnew[i] - g[i];
  c = dx*Dot(n, y, p);
  if( c < 1e-10 ) return;
  RenormalizeCholesky(t, n, hessian, y, 1/c);

  c = Dot(n, g, p);
  if( c >= 0 ) return;
  c = 1/sqrt(-c);
  for( i = 0; i < n; ++i )
    y[i] = c*g[i];
  UpdateCholesky(t, n, hessian, y, p);

  for( i = 0; i < n - 1; ++i )
    for( j = i + 1; j < n; ++j )
      Hessian(i, j) = Hessian(j, i);
}

/*********************************************************************/

static void Gradient(This *t, ccount nfree, ccount *ifree,
  cBounds *b, real *x, creal y, real *grad)
{
  count i;

  for( i = 0; i < nfree; ++i ) {
    ccount dim = Untag(ifree[i]);
    creal xd = x[dim];
    creal delta = (b[dim].upper - xd < DELTA) ? -DELTA : DELTA;
    x[dim] += delta;
    grad[i] = (Sample(t, x) - y)/delta;
    x[dim] = xd;
  }
}

/*********************************************************************/

static Point LineSearch(This *t, ccount nfree, ccount *ifree,
  creal *p, creal *xini, real fini, real *x,
  real step, creal range, creal grad,
  creal ftol, creal xtol, creal gtol)
{
  real tol = ftol, tol2 = tol + tol;
  Point cur = {0, fini};

  XCopy(x, xini);

  /* don't even try if
     a) we'd walk backwards,
     b) the range to explore is too small,
     c) the gradient is positive, i.e. we'd move uphill */

  if( step > 0 && range > tol2 && grad <= 0 ) {
    creal eps = RTEPS*fabs(range) + ftol;
    creal mingrad = -1e-4*grad, maxgrad = -gtol*grad;

    real end = range + eps;
    real maxstep = range - eps/(1 + RTEPS);

    Point min = cur, v = cur, w = cur;
    Point a = cur, b = {end, 0};
    real a1, b1 = end;

    /* distmin: distance along p from xini to the minimum,
       u: second-lowest point,
       v: third-lowest point,
       a, b: interval in which the minimum is sought. */

    real distmin = 0, dist, mid, q, r, s;
    count i;
    int shift;
    bool first;

    for( first = true; ; first = false ) {
      if( step >= maxstep ) {
        step = maxstep;
        maxstep = maxstep*(1 + .75*RTEPS) + .75*tol;
      }

      cur.dx = (fabs(step) >= tol) ? step : (step > 0) ? tol : -tol;
      dist = distmin + cur.dx;
      for( i = 0; i < nfree; ++i ) {
        ccount dim = ifree[i];
        x[dim] = xini[dim] + dist*p[i];
      }
      cur.f = Sample(t, x);

      if( cur.f <= min.f ) {
        v = w;
        w = min;
        min.f = cur.f;
        distmin = dist;

        /* shift everything to the new minimum position */
        maxstep -= cur.dx;
        v.dx -= cur.dx;
        w.dx -= cur.dx;
        a.dx -= cur.dx;
        b.dx -= cur.dx;
        if( cur.dx < 0 ) b = w;
        else a = w;

        tol = RTEPS*fabs(distmin) + ftol;
        tol2 = tol + tol;
      }
      else {
        if( cur.dx < 0 ) a = cur;
        else b = cur;
        if( cur.f <= w.f || w.dx == 0 ) v = w, w = cur;
        else if( cur.f <= v.f || v.dx == 0 || v.dx == w.dx ) v = cur;
      }

      if( distmin + b.dx <= xtol ) break;
      if( min.f < fini &&
          a.f - min.f <= fabs(a.dx)*maxgrad &&
          (fabs(distmin - range) > tol || maxstep < b.dx) ) break;

      mid = .5*(a.dx + b.dx);
      if( fabs(mid) <= tol2 - .5*(b.dx - a.dx) ) break;

      r = q = s = 0;
      if( fabs(end) > tol ) {
        if( first ) {
          creal s1 = w.dx*grad;
          creal s2 = w.f - min.f;
          s = (s1 - ((distmin == 0) ? 0 : 2*s2))*w.dx;
          q = 2*(s2 - s1);
        }
        else {
          creal s1 = w.dx*(v.f - min.f);
          creal s2 = v.dx*(w.f - min.f);
          s = s1*w.dx - s2*v.dx;
          q = 2*(s2 - s1);
        }
        if( q > 0 ) s = -s;
        q = fabs(q);
        r = end;
        if( step != b1 || b.dx <= maxstep ) end = step;
      }

      if( distmin == a.dx ) step = mid;
      else if( b.dx > maxstep ) step = (step < b.dx) ? -4*a.dx : maxstep;
      else {
        real num = a.dx, den = b.dx;
        if( fabs(b.dx) <= tol || (w.dx > 0 && fabs(a.dx) > tol) )
          num = b.dx, den = a.dx;
        num /= -den;
        step = (num < 1) ? .5*den*sqrt(num) : 5/11.*den*(.1 + 1/num);
      }

      if( step > 0 ) a1 = a.dx, b1 = step;
      else a1 = step, b1 = b.dx;
      if( fabs(s) < fabs(.5*q*r) && s > q*a1 && s < q*b1 ) {
        step = s/q;
        if( step - a.dx < tol2 || b.dx - step < tol2 )
          step = (mid > 0) ? tol : -tol;
      }
      else end = (mid > 0) ? b.dx : a.dx;
    }

    first = true;
    if( fabs(distmin - range) < tol ) {
      distmin = range;
      if( maxstep > b.dx ) first = false;
    }

    for( cur.dx = distmin, cur.f = min.f, shift = -1; ;
         cur.dx = Max(ldexp(distmin, shift), ftol), shift <<= 1 ) {
      for( i = 0; i < nfree; ++i ) {
        ccount dim = ifree[i];
        x[dim] = xini[dim] + cur.dx*p[i];
      }
      if( !first ) cur.f = Sample(t, x);

      if( cur.dx + b.dx <= xtol ) {
        cur.dx = 0;
        break;
      }
      if( fini - cur.f > cur.dx*mingrad ) break;
      if( cur.dx <= ftol ) {
        cur.dx = 0;
        break;
      }
      first = false;
    }
  }

  return cur;
}

/*********************************************************************/

static real LocalSearch(This *t, ccount nfree, ccount *ifree,
  cBounds *b, creal *x, creal fx, real *z)
{
  real delta, smax, sopp, spmax, snmax;
  real y[NDIM], fy, fz, ftest;
  real p[NDIM];
  int sign;
  count i;

  /* Choose a direction p along which to move away from the
     present x.  We choose the direction which leads farthest
     away from all borders. */

  smax = INFTY;
  for( i = 0; i < nfree; ++i ) {
    ccount dim = ifree[i];
    creal sp = b[dim].upper - x[dim];
    creal sn = x[dim] - b[dim].lower;
    if( sp < sn ) {
      smax = Min(smax, sn);
      p[i] = -1;
    }
    else {
      smax = Min(smax, sp);
      p[i] = 1;
    }
  }
  smax *= .9;

  /* Move along p until the integrand changes appreciably
     or we come close to a border. */

  XCopy(y, x);
  ftest = SUFTOL*(1 + fabs(fx));
  delta = RTDELTA/5;
  do {
    delta = Min(5*delta, smax);
    for( i = 0; i < nfree; ++i ) {
      ccount dim = ifree[i];
      y[dim] = x[dim] + delta*p[i];
    }
    fy = Sample(t, y);
    if( fabs(fy - fx) > ftest ) break;
  } while( delta != smax );

  /* Construct a second direction p' orthogonal to p, i.e. p.p' = 0.
     We let pairs of coordinates cancel in the dot product,
     i.e. we choose p'[0] = p[0], p'[1] = -p[1], etc.
     (It should really be 1/p and -1/p, but p consists of 1's and -1's.)
     For odd nfree, we let the last three components cancel by 
     choosing p'[nfree - 3] = p[nfree - 3],
              p'[nfree - 2] = -1/2 p[nfree - 2], and
              p'[nfree - 1] = -1/2 p[nfree - 1]. */

  sign = (nfree <= 1 && fy > fx) ? 1 : -1;
  spmax = snmax = INFTY;
  for( i = 0; i < nfree; ++i ) {
    ccount dim = ifree[i];
    real sp, sn;
    p[i] *= (nfree & 1 && nfree - i <= 2) ? -.5*sign : (sign = -sign);
    sp = (b[dim].upper - y[dim])/p[i];
    sn = (y[dim] - b[dim].lower)/p[i];
    if( p[i] > 0 ) {
      spmax = Min(spmax, sp);
      snmax = Min(snmax, sn);
    }
    else {
      spmax = Min(spmax, -sn);
      snmax = Min(snmax, -sp);
    }
  }
  smax = .9*spmax;
  sopp = .9*snmax;

  if( nfree > 1 && smax < snmax ) {
    real tmp = smax;
    smax = sopp;
    sopp = tmp;
    for( i = 0; i < nfree; ++i )
      p[i] = -p[i];
  }

  /* Move along p' until the integrand changes appreciably
     or we come close to a border. */

  XCopy(z, y);
  ftest = SUFTOL*(1 + fabs(fy));
  delta = RTDELTA/5;
  do {
    delta = Min(5*delta, smax);
    for( i = 0; i < nfree; ++i ) {
      ccount dim = ifree[i];
      z[dim] = y[dim] + delta*p[i];
    }
    fz = Sample(t, z);
    if( fabs(fz - fy) > ftest ) break;
  } while( delta != smax );

  if( fy != fz ) {
    real pleneps, grad, range, step;
    Point low;

    if( fy > fz ) {
      grad = (fz - fy)/delta;
      range = smax/.9;
      step = Min(delta + delta, smax);
    }
    else {
      grad = (fy - fz)/delta;
      range = sopp/.9 + delta;
      step = Min(delta + delta, sopp);
      XCopy(y, z);
      fy = fz;
      for( i = 0; i < nfree; ++i )
        p[i] = -p[i];
    }

    pleneps = Length(nfree, p) + RTEPS;
    low = LineSearch(t, nfree, ifree, p, y, fy, z, step, range, grad,
      RTEPS/pleneps, 0., RTEPS);
    fz = low.f;
  }

  if( fz != fx ) {
    real pleneps, grad, range, step;
    Point low;

    spmax = snmax = INFTY;
    for( i = 0; i < nfree; ++i ) {
      ccount dim = ifree[i];
      p[i] = z[dim] - x[dim];
      if( p[i] != 0 ) {
        creal sp = (b[dim].upper - x[dim])/p[i];
        creal sn = (x[dim] - b[dim].lower)/p[i];
        if( p[i] > 0 ) {
          spmax = Min(spmax, sp);
          snmax = Min(snmax, sn);
        }
        else {
          spmax = Min(spmax, -sn);
          snmax = Min(snmax, -sp);
        }
      }
    }

    grad = fz - fx;
    range = spmax;
    step = Min(.9*spmax, 2.);
    pleneps = Length(nfree, p) + RTEPS;
    if( fz > fx ) {
      delta = Min(.9*snmax, RTDELTA/pleneps);
      for( i = 0; i < nfree; ++i ) {
        ccount dim = ifree[i];
        z[dim] = x[dim] - delta*p[i];
      }
      fz = Sample(t, z);
      if( fz < fx ) {
        grad = (fz - fx)/delta;
        range = snmax;
        step = Min(.9*snmax, delta + delta);
        for( i = 0; i < nfree; ++i )
          p[i] = -p[i];
      }
      else if( delta < 1 ) grad = (fx - fz)/delta;
    }

    low = LineSearch(t, nfree, ifree, p, x, fx, z, step, range, grad,
      RTEPS/pleneps, 0., RTEPS);
    fz = low.f;
  }

  return fz;
}

/*********************************************************************/

static real FindMinimum(This *t, cBounds *b, real *xmin, real fmin)
{
  real hessian[NDIM*NDIM];
  real gfree[NDIM], p[NDIM];
  real tmp[NDIM], ftmp, fini = fmin;
  ccount maxeval = t->neval + 50*t->ndim;
  count nfree, nfix;
  count ifree[NDIM], ifix[NDIM];
  count dim, local;

  Zap(hessian);
  for( dim = 0; dim < t->ndim; ++dim )
    Hessian(dim, dim) = 1;

  /* Step 1: - classify the variables as "fixed" (sufficiently close
               to a border) and "free",
             - if the integrand is flat in the direction of the gradient
               w.r.t. the free dimensions, perform a local search. */

  for( local = 0; local < 2; ++local ) {
    bool resample = false;
    nfree = nfix = 0;
    for( dim = 0; dim < t->ndim; ++dim ) {
      if( xmin[dim] < b[dim].lower + (1 + fabs(b[dim].lower))*QEPS ) {
        xmin[dim] = b[dim].lower;
        ifix[nfix++] = dim;
        resample = true;
      }
      else if( xmin[dim] > b[dim].upper - (1 + fabs(b[dim].upper))*QEPS ) {
        xmin[dim] = b[dim].upper;
        ifix[nfix++] = Tag(dim);
        resample = true;
      }
      else ifree[nfree++] = dim;
    }

    if( resample ) fini = fmin = Sample(t, xmin);

    if( nfree == 0 ) goto releasebounds;

    Gradient(t, nfree, ifree, b, xmin, fmin, gfree);
    if( local || Length(nfree, gfree) > GTOL ) break;

    ftmp = LocalSearch(t, nfree, ifree, b, xmin, fmin, tmp);
    if( ftmp > fmin - (1 + fabs(fmin))*RTEPS )
      goto releasebounds;
    fmin = ftmp;
    XCopy(xmin, tmp);
  }

  while( t->neval <= maxeval ) {

    /* Step 2a: perform a quasi-Newton iteration on the free
                variables only. */

    if( nfree > 0 ) {
      real plen, pleneps;
      real minstep;
      count i, mini = 0, minfix = 0;
      Point low;

      LinearSolve(t, nfree, hessian, gfree, p);
      plen = Length(nfree, p);
      pleneps = plen + RTEPS;

      minstep = INFTY;
      for( i = 0; i < nfree; ++i ) {
        count dim = Untag(ifree[i]);
        if( fabs(p[i]) > EPS ) {
          real step;
          count fix;
          if( p[i] < 0 ) {
            step = (b[dim].lower - xmin[dim])/p[i];
            fix = dim;
          }
          else {
            step = (b[dim].upper - xmin[dim])/p[i];
            fix = Tag(dim);
          }
          if( step < minstep ) {
            minstep = step;
            mini = i;
            minfix = fix;
          }
        }
      }

      if( minstep*pleneps <= DELTA ) {
fixbound:
        ifix[nfix++] = minfix;

        if( mini < --nfree ) {
          creal diag = Hessian(mini, mini);

          Clear(tmp, mini);
          for( i = mini; i < nfree; ++i )
            tmp[i] = Hessian(i + 1, mini);

          for( i = mini; i < nfree; ++i ) {
            Move(&Hessian(i, 0), &Hessian(i + 1, 0), i);
            Hessian(i, i) = Hessian(i + 1, i + 1);
          }
          RenormalizeCholesky(t, nfree, hessian, tmp, diag);

          Move(&ifree[mini], &ifree[mini + 1], nfree - mini);
          Move(&gfree[mini], &gfree[mini + 1], nfree - mini);
        }
        continue;
      }

      low = LineSearch(t, nfree, ifree, p, xmin, fmin, tmp,
        Min(minstep, 1.), Min(minstep, 100.), Dot(nfree, gfree, p),
        RTEPS/pleneps, DELTA/pleneps, .2);

      if( low.dx > 0 ) {
        real fdiff;

        fmin = low.f;
        XCopy(xmin, tmp);

        Gradient(t, nfree, ifree, b, xmin, fmin, tmp);
        BFGS(t, nfree, hessian, tmp, gfree, p, low.dx);
        XCopy(gfree, tmp);

        if( fabs(low.dx - minstep) < QEPS*minstep ) goto fixbound;

        fdiff = fini - fmin;
        fini = fmin;
        if( fdiff > (1 + fabs(fmin))*FTOL ||
            low.dx*plen > (1 + Length(t->ndim, xmin))*FTOL ) continue;
      }
    }

    /* Step 2b: check whether freeing any fixed variable will lead
                to a reduction in f. */

releasebounds:
    if( nfix > 0 ) {
      real mingrad = INFTY;
      count i, mini = 0;
      bool repeat = false;

      Gradient(t, nfix, ifix, b, xmin, fmin, tmp);

      for( i = 0; i < nfix; ++i ) {
        creal grad = Sign(ifix[i])*tmp[i];
        if( grad < -RTEPS ) {
          repeat = true;
          if( grad < mingrad ) {
            mingrad = grad;
            mini = i;
          }
        }
      }

      if( repeat ) {
        gfree[nfree] = tmp[mini];
        ifree[nfree] = Untag(ifix[mini]);
        Clear(&Hessian(nfree, 0), nfree);
        Hessian(nfree, nfree) = 1;
        ++nfree;

        --nfix;
        Move(&ifix[mini], &ifix[mini + 1], nfix - mini);
        continue;
      }
    }

    break;
  }

  return fmin;
}

