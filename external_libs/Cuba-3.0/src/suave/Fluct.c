/*
	Fluct.c
		compute the fluctuation in the left and right half
		this file is part of Suave
		last modified 23 Aug 11 th
*/


#if defined(HAVE_LONG_DOUBLE) && defined(HAVE_POWL)

typedef long double realx;
#define XDBL_MAX_EXP LDBL_MAX_EXP
#define XDBL_MAX LDBL_MAX
#define powx powl
#define ldexpx ldexpl

#else

typedef double realx;
#define XDBL_MAX_EXP DBL_MAX_EXP
#define XDBL_MAX DBL_MAX
#define powx pow
#define ldexpx ldexp

#endif

typedef const realx crealx;

typedef struct {
  realx fluct;
  number n;
} Var;

/*********************************************************************/

static void Fluct(cThis *t, Var *var,
  cBounds *b, creal *w, number n, ccount comp, creal avg, creal err)
{
  creal *x = w + n;
  creal *f = x + n*t->ndim + comp;
  count nvar = 2*t->ndim;
  creal norm = 1/(err*Max(fabs(avg), err));
  creal flat = 2/3./t->flatness;
  crealx max = ldexpx(1., (int)((XDBL_MAX_EXP - 2)/t->flatness));

  Clear(var, nvar);

  while( n-- ) {
    count dim;
    crealx arg = 1 + fabs(*w++)*Sq(*f - avg)*norm;
    crealx ft = powx(arg < max ? arg : max, t->flatness);

    f += t->ncomp;

    for( dim = 0; dim < t->ndim; ++dim ) {
      Var *v = &var[2*dim + (*x++ >= b[dim].mid)];
      crealx f = v->fluct + ft;
      v->fluct = (f > XDBL_MAX/2) ? XDBL_MAX/2 : f;
      ++v->n;
    }
  }

  while( nvar-- ) {
    var->fluct = powx(var->fluct, flat);
    ++var;
  }
}

