/*
	WorkerIni.c
		set/run the init/exit functions for worker processes
		by Thomas Hahn
		last modified 6 Sep 12 th
*/


#include "stddecl.h"

extern workerini cubaini;

Extern void SUFFIX(cubasetinit)(subroutine f, void *arg)
{
  cubaini.initfun = f;
  cubaini.initarg = arg;
}


Extern void SUFFIX(cubasetexit)(subroutine f, void *arg)
{
  cubaini.exitfun = f;
  cubaini.exitarg = arg;
}


Extern void SUFFIX(cubaruninit)()
{
  if( cubaini.initfun ) cubaini.initfun(cubaini.initarg);
}


Extern void SUFFIX(cubarunexit)()
{
  if( cubaini.exitfun ) cubaini.exitfun(cubaini.exitarg);
}

