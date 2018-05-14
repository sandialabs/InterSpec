* cuba.F
* Fortran chooser for the Cuba routines
* last modified 3 Feb 05 th

#define VEGAS 1
#define SUAVE 2
#define DIVONNE 3
#define CUHRE 4


	subroutine Cuba(method, ndim, ncomp, integrand,
     &    integral, error, prob)
	implicit none
	integer method, ndim, ncomp
	external integrand
	double precision integral(*), error(*), prob(*)

	character*7 name(4)
	data name /"Vegas", "Suave", "Divonne", "Cuhre"/

	integer mineval, maxeval, verbose, last
	double precision epsrel, epsabs
	parameter (epsrel = 1D-3)
	parameter (epsabs = 1D-12)
	parameter (verbose = 2)
	parameter (last = 4)
	parameter (mineval = 0)
	parameter (maxeval = 50000)

	integer nstart, nincrease
	parameter (nstart = 1000)
	parameter (nincrease = 500)

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

	integer nregions, neval, fail


	if( method .eq. VEGAS ) then

	  call vegas(ndim, ncomp, integrand,
     &      epsrel, epsabs, verbose, mineval, maxeval,
     &      nstart, nincrease,
     &      neval, fail, integral, error, prob)
	  nregions = 1

	else if( method .eq. SUAVE ) then

	  call suave(ndim, ncomp, integrand,
     &      epsrel, epsabs, verbose + last, mineval, maxeval,
     &      nnew, flatness,
     &      nregions, neval, fail, integral, error, prob)

	else if( method .eq. DIVONNE ) then

	  call divonne(ndim, ncomp, integrand,
     &      epsrel, epsabs, verbose, mineval, maxeval,
     &      key1, key2, key3, maxpass,
     &      border, maxchisq, mindeviation,
     &      ngiven, ldxgiven, 0, nextra, 0,
     &      nregions, neval, fail, integral, error, prob)

	else if( method .eq. CUHRE ) then

	  call cuhre(ndim, ncomp, integrand,
     &      epsrel, epsabs, verbose + last, mineval, maxeval,
     &      key,
     &      nregions, neval, fail, integral, error, prob)

	else

	  print *, "invalid method ", method
	  return

	endif

	print *, "method   =", name(method)
	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(G20.12," +- ",G20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
	end

