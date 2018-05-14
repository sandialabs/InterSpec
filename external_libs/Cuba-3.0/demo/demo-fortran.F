* demo-fortran.F
* test program for the Cuba library
* last modified 27 Jun 12 th


	program CubaTest
	implicit none

	integer ndim, ncomp, last, seed, mineval, maxeval
	double precision epsrel, epsabs, userdata
	parameter (ndim = 3)
	parameter (ncomp = 1)
	parameter (userdata = 0)
	parameter (epsrel = 1D-3)
	parameter (epsabs = 1D-12)
	parameter (last = 4)
	parameter (seed = 0)
	parameter (mineval = 0)
	parameter (maxeval = 50000)

	integer nstart, nincrease, nbatch, gridno
	character*(*) statefile	
	parameter (nstart = 1000)
	parameter (nincrease = 500)
	parameter (nbatch = 1000)
	parameter (gridno = 0)
	parameter (statefile = "")

	integer nnew
	double precision flatness
	parameter (nnew = 1000)
	parameter (flatness = 25D0)

	integer key1, key2, key3, maxpass
	double precision border, maxchisq, mindeviation
	integer ngiven, ldxgiven, nextra
	parameter (key1 = 47)
	parameter (key2 = 1)
	parameter (key3 = 1)
	parameter (maxpass = 5)
	parameter (border = 0D0)
	parameter (maxchisq = 10D0)
	parameter (mindeviation = .25D0)
	parameter (ngiven = 0)
	parameter (ldxgiven = ndim)
	parameter (nextra = 0)

	integer key
	parameter (key = 0)

	external integrand

	double precision integral(ncomp), error(ncomp), prob(ncomp)
	integer verbose, nregions, neval, fail
	character*16 env

	integer c

	call getenv("CUBAVERBOSE", env)
	verbose = 2
	read(env, *, iostat=fail, end=999, err=999) verbose
999	continue

	print *, "-------------------- Vegas test --------------------"

	call vegas(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, nstart, nincrease, nbatch,
     &    gridno, statefile,
     &    neval, fail, integral, error, prob)

	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)

	print *, " "
	print *, "-------------------- Suave test --------------------"

	call suave(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose + last, seed,
     &    mineval, maxeval, nnew, flatness,
     &    nregions, neval, fail, integral, error, prob)

	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)

	print *, " "
	print *, "------------------- Divonne test -------------------"

	call divonne(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    nregions, neval, fail, integral, error, prob)

	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)

	print *, " "
	print *, "-------------------- Cuhre test --------------------"

	call cuhre(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose + last,
     &    mineval, maxeval, key,
     &    nregions, neval, fail, integral, error, prob)

	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
	end


************************************************************************

	integer function integrand(ndim, xx, ncomp, ff)
	implicit none
	integer ndim, ncomp
	double precision xx(*), ff(*)

#define x xx(1)
#define y xx(2)
#define z xx(3)
#define f ff(1)

#ifndef FUN
#define FUN 1
#endif

	double precision pi, rsq
	parameter (pi = 3.14159265358979323846D0)

	rsq = x**2 + y**2 + z**2

#if FUN == 1
	f = sin(x)*cos(y)*exp(z)
#elif FUN == 2
	f = 1/((x + y)**2 + .003D0)*cos(y)*exp(z)
#elif FUN == 3
	f = 1/(3.75D0 - cos(pi*x) - cos(pi*y) - cos(pi*z))
#elif FUN == 4
	f = abs(rsq - .125D0)
#elif FUN == 5
	f = exp(-rsq)
#elif FUN == 6
	f = 1/(1 - x*y*z + 1D-10)
#elif FUN == 7
	f = sqrt(abs(x - y - z))
#elif FUN == 8
	f = exp(-x*y*z)
#elif FUN == 9
	f = x**2/(cos(x + y + z + 1) + 5)
#elif FUN == 10
	if( x .gt. .5D0 ) then
	  f = 1/sqrt(x*y*z + 1D-5)
	else
	  f = sqrt(x*y*z)
	endif
#else
	if( rsq .lt. 1 ) then
	  f = 1
	else
	  f = 0
	endif
#endif

	integrand = 0
	end

