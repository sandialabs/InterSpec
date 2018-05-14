/*
	ChiSquare.c
		the chi-square cdf
		after W.J. Kennedy and J.E. Gentle,
		Statistical computing, p. 116
		last modified 9 Feb 05 th
*/

#ifdef HAVE_ERF
#define Erf erf
#else
#include "Erf.c"
#endif

#if( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static inline 
#endif
real Normal(creal x)
{
  return .5*Erf(x/1.414213562373095048801689) + .5;
}

/*********************************************************************/

#if( !( defined(WIN32) || defined(UNDER_CE) || defined(_WIN32) || defined(WIN64) ) )
static
#endif
real ChiSquare(creal x, cint df)
{
  real y;

  if( df <= 0 ) return -999;

  if( x <= 0 ) return 0;
  if( x > 1000*df ) return 1;

  if( df > 1000 ) {
    if( x < 2 ) return 0;
    y = 2./(9*df);
    y = (pow(x/df, 1/3.) - (1 - y))/sqrt(y);
    if( y > 5 ) return 1;
    if( y < -18.8055 ) return 0;
    return Normal(y);
  }

  y = .5*x;

  if( df & 1 ) {
    creal sqrty = sqrt(y);
    real h = Erf(sqrty);
    count i;

    if( df == 1 ) return h;

    y = sqrty*exp(-y)/.8862269254527579825931;
    for( i = 3; i < df; i += 2 ) {
      h -= y;
      y *= x/i;
    }
    y = h - y;
  }
  else {
    real term = exp(-y), sum = term;
    count i;

    for( i = 1; i < df/2; ++i )
      sum += term *= y/i;
    y = 1 - sum;
  }

  return Max(0., y);
}

