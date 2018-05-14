/*
	Random.c
		quasi- and pseudo-random-number generation
		last modified 17 Dec 11 th
*/


/*
	PART 1: Sobol quasi-random-number generator
	adapted from ACM TOMS algorithm 659
*/

static void SobolGet(This *t, real *x)
{
  number seq = t->rng.sobol.seq++;
  count zerobit = 0, dim;

  while( seq & 1 ) {
    ++zerobit;
    seq >>= 1;
  }

  for( dim = 0; dim < t->ndim; ++dim ) {
    t->rng.sobol.prev[dim] ^= t->rng.sobol.v[dim][zerobit];
    x[dim] = t->rng.sobol.prev[dim]*t->rng.sobol.norm;
  }
}


static void SobolSkip(This *t, number n)
{
  while( n-- ) {
    number seq = t->rng.sobol.seq++;
    count zerobit = 0, dim;

    while( seq & 1 ) {
      ++zerobit;
      seq >>= 1;
    }

    for( dim = 0; dim < t->ndim; ++dim )
      t->rng.sobol.prev[dim] ^= t->rng.sobol.v[dim][zerobit];
  }
}


static inline void SobolIni(This *t)
{
  static number ini[9*40] = {
      3,   1,   0,   0,   0,   0,   0,   0,   0,
      7,   1,   1,   0,   0,   0,   0,   0,   0,
     11,   1,   3,   7,   0,   0,   0,   0,   0,
     13,   1,   1,   5,   0,   0,   0,   0,   0,
     19,   1,   3,   1,   1,   0,   0,   0,   0,
     25,   1,   1,   3,   7,   0,   0,   0,   0,
     37,   1,   3,   3,   9,   9,   0,   0,   0,
     59,   1,   3,   7,  13,   3,   0,   0,   0,
     47,   1,   1,   5,  11,  27,   0,   0,   0,
     61,   1,   3,   5,   1,  15,   0,   0,   0,
     55,   1,   1,   7,   3,  29,   0,   0,   0,
     41,   1,   3,   7,   7,  21,   0,   0,   0,
     67,   1,   1,   1,   9,  23,  37,   0,   0,
     97,   1,   3,   3,   5,  19,  33,   0,   0,
     91,   1,   1,   3,  13,  11,   7,   0,   0,
    109,   1,   1,   7,  13,  25,   5,   0,   0,
    103,   1,   3,   5,  11,   7,  11,   0,   0,
    115,   1,   1,   1,   3,  13,  39,   0,   0,
    131,   1,   3,   1,  15,  17,  63,  13,   0,
    193,   1,   1,   5,   5,   1,  27,  33,   0,
    137,   1,   3,   3,   3,  25,  17, 115,   0,
    145,   1,   1,   3,  15,  29,  15,  41,   0,
    143,   1,   3,   1,   7,   3,  23,  79,   0,
    241,   1,   3,   7,   9,  31,  29,  17,   0,
    157,   1,   1,   5,  13,  11,   3,  29,   0,
    185,   1,   3,   1,   9,   5,  21, 119,   0,
    167,   1,   1,   3,   1,  23,  13,  75,   0,
    229,   1,   3,   3,  11,  27,  31,  73,   0,
    171,   1,   1,   7,   7,  19,  25, 105,   0,
    213,   1,   3,   5,   5,  21,   9,   7,   0,
    191,   1,   1,   1,  15,   5,  49,  59,   0,
    253,   1,   1,   1,   1,   1,  33,  65,   0,
    203,   1,   3,   5,  15,  17,  19,  21,   0,
    211,   1,   1,   7,  11,  13,  29,   3,   0,
    239,   1,   3,   7,   5,   7,  11, 113,   0,
    247,   1,   1,   5,   3,  15,  19,  61,   0,
    285,   1,   3,   1,   1,   9,  27,  89,   7,
    369,   1,   1,   3,   7,  31,  15,  45,  23,
    299,   1,   3,   3,   9,   9,  25, 107,  39 };

  count dim, bit, nbits;
  number max, *pini = ini;
  cnumber nmax = 2*t->maxeval;

  for( nbits = 0, max = 1; max <= nmax; max <<= 1 ) ++nbits;
  t->rng.sobol.norm = 1./max;

  for( bit = 0; bit < nbits; ++bit )
    t->rng.sobol.v[0][bit] = (max >>= 1);

  for( dim = 1; dim < t->ndim; ++dim ) {
    number *pv = t->rng.sobol.v[dim], *pvv = pv;
    number powers = *pini++, j;
    int inibits = -1, bit;
    for( j = powers; j; j >>= 1 ) ++inibits;

    memcpy(pv, pini, inibits*sizeof(*pini));
    pini += 8;

    for( bit = inibits; bit < nbits; ++bit ) {
      number newv = *pvv, j = powers;
      int b;
      for( b = 0; b < inibits; ++b ) {
        if( j & 1 ) newv ^= pvv[b] << (inibits - b);
        j >>= 1;
      }
      pvv[inibits] = newv;
      ++pvv;
    }

    for( bit = 0; bit < nbits - 1; ++bit )
      pv[bit] <<= nbits - bit - 1;
  }

  t->rng.sobol.seq = 0;
  XClear(t->rng.sobol.prev);

  t->rng.getrandom = SobolGet;
  t->rng.skiprandom = SobolSkip;
}


/*
	PART 2: Mersenne Twister pseudo-random-number generator
	adapted from T. Nishimura's and M. Matsumoto's C code at
	http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*/

/* 32 or 53 random bits */
#define RANDOM_BITS 32


static inline state_t Twist(state_t a, state_t b)
{
  state_t mixbits = (a & 0x80000000) | (b & 0x7fffffff);
  state_t matrixA = (-(b & 1)) & 0x9908b0df;
  return (mixbits >> 1) ^ matrixA;
}


static inline void MersenneReload(state_t *state)
{
  state_t *s = state;
  int j;

  for( j = MERSENNE_N - MERSENNE_M + 1; --j; ++s )
    *s = s[MERSENNE_M] ^ Twist(s[0], s[1]);
  for( j = MERSENNE_M; --j; ++s )
    *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], s[1]);
  *s = s[MERSENNE_M - MERSENNE_N] ^ Twist(s[0], state[0]);
}


static inline state_t MersenneInt(state_t s)
{
  s ^= s >> 11;
  s ^= (s << 7) & 0x9d2c5680;
  s ^= (s << 15) & 0xefc60000;
  return s ^ (s >> 18);
}


static void MersenneGet(This *t, real *x)
{
  count next = t->rng.mersenne.next, dim;

  for( dim = 0; dim < t->ndim; ++dim ) {
#if RANDOM_BITS == 53
    state_t a, b;
#endif

    if( next >= MERSENNE_N ) {
      MersenneReload(t->rng.mersenne.state);
      next = 0;
    }

#if RANDOM_BITS == 53
    a = MersenneInt(t->rng.mersenne.state[next++]) >> 5;
    b = MersenneInt(t->rng.mersenne.state[next++]) >> 6;
    x[dim] = (67108864.*a + b)/9007199254740992.;
#else
    x[dim] = MersenneInt(t->rng.mersenne.state[next++])/4294967296.;
#endif
  }

  t->rng.mersenne.next = next;
}


static void MersenneSkip(This *t, number n)
{
#if RANDOM_BITS == 53
  n = 2*n*t->ndim + t->rng.mersenne.next;
#else
  n = n*t->ndim + t->rng.mersenne.next;
#endif
  t->rng.mersenne.next = n % MERSENNE_N;
  n /= MERSENNE_N;
  while( n-- ) MersenneReload(t->rng.mersenne.state);
}


static inline void MersenneIni(This *t)
{
  state_t seed = t->seed;
  state_t *next = t->rng.mersenne.state;
  count j;

  for( j = 1; j <= MERSENNE_N; ++j ) {
    *next++ = seed;
    seed = 0x6c078965*(seed ^ (seed >> 30)) + j;
    /* see Knuth TAOCP Vol 2, 3rd Ed, p. 106 for multiplier */
  }

  MersenneReload(t->rng.mersenne.state);
  t->rng.mersenne.next = 0;

  t->rng.getrandom = MersenneGet;
  t->rng.skiprandom = MersenneSkip;
}


/*
	PART 3: Ranlux subtract-and-borrow random-number generator 
	proposed by Marsaglia and Zaman, implemented by F. James with 
	the name RCARRY in 1991, and later improved by Martin Luescher 
	in 1993 to produce "Luxury Pseudorandom Numbers".
	Adapted from the CERNlib Fortran 77 code by F. James, 1993.

	The available luxury levels are:

	level 0  (p = 24): equivalent to the original RCARRY of Marsaglia
	         and Zaman, very long period, but fails many tests.
	level 1  (p = 48): considerable improvement in quality over level 0,
	         now passes the gap test, but still fails spectral test.
	level 2  (p = 97): passes all known tests, but theoretically still
	         defective.
	level 3  (p = 223): DEFAULT VALUE.  Any theoretically possible
	         correlations have very small chance of being observed.
	level 4  (p = 389): highest possible luxury, all 24 bits chaotic.
*/


static inline int RanluxInt(This *t, count n)
{
  int s = 0;

  while( n-- ) {
    s = t->rng.ranlux.state[t->rng.ranlux.j24] -
        t->rng.ranlux.state[t->rng.ranlux.i24] + t->rng.ranlux.carry;
    s += (t->rng.ranlux.carry = NegQ(s)) & (1 << 24);
    t->rng.ranlux.state[t->rng.ranlux.i24] = s;
    --t->rng.ranlux.i24;
    t->rng.ranlux.i24 += NegQ(t->rng.ranlux.i24) & 24;
    --t->rng.ranlux.j24;
    t->rng.ranlux.j24 += NegQ(t->rng.ranlux.j24) & 24;
  }

  return s;
}


static void RanluxGet(This *t, real *x)
{
/* The Generator proper: "Subtract-with-borrow",
   as proposed by Marsaglia and Zaman, FSU, March 1989 */

  count dim;

  for( dim = 0; dim < t->ndim; ++dim ) {
    cint nskip = (--t->rng.ranlux.n24 >= 0) ? 0 :
      (t->rng.ranlux.n24 = 24, t->rng.ranlux.nskip);
    cint s = RanluxInt(t, 1 + nskip);
    x[dim] = s*0x1p-24;
/* small numbers (with less than 12 significant bits) are "padded" */
    if( s < (1 << 12) )
      x[dim] += t->rng.ranlux.state[t->rng.ranlux.j24]*0x1p-48;
  }
}


static void RanluxSkip(This *t, cnumber n)
{
  RanluxInt(t, n + t->rng.ranlux.nskip*(n/24));
  t->rng.ranlux.n24 = 24 - n % 24;
}


static inline void RanluxIni(This *t)
{
  cint skip[] = {24, 48, 97, 223, 389,
    223, 223, 223, 223, 223, 223, 223, 223, 223, 223,
    223, 223, 223, 223, 223, 223, 223, 223, 223, 223};
  state_t seed = t->seed;
  state_t level = RNG;
  count i;

  if( level < sizeof skip ) level = skip[level];
  t->rng.ranlux.nskip = level - 24;

  t->rng.ranlux.i24 = 23;
  t->rng.ranlux.j24 = 9;
  t->rng.ranlux.n24 = 24;

  for( i = 0; i < 24; ++i ) {
    cint k = seed/53668;
    seed = 40014*(seed - k*53668) - k*12211;
    seed += NegQ(seed) & 2147483563;
    t->rng.ranlux.state[i] = seed & ((1 << 24) - 1);
  }

  t->rng.ranlux.carry = ~TrueQ(t->rng.ranlux.state[23]) & (1 << 24);

  t->rng.getrandom = RanluxGet;
  t->rng.skiprandom = RanluxSkip;
}


/*
	PART 4: User routines:

	- IniRandom sets up the random-number generator to produce a
	  sequence of at least n ndim-dimensional random vectors.

	- GetRandom retrieves one random vector.

	- SkipRandom skips over n random vectors.
*/

static inline void IniRandom(This *t)
{
  if( t->seed == 0 ) SobolIni(t);
  else if( RNG == 0 ) MersenneIni(t);
  else RanluxIni(t);
}

