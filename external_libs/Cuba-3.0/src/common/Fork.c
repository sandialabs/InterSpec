/*
	Fork.c
		the parallel sampling routine
		for the C versions of the Cuba routines
		by Thomas Hahn
		last modified 6 Sep 12 th
*/

#define MINSLICE 10
#define MINCORES 1
//#define MINCORES 2

typedef struct {
  number n, m, i;
  VES_ONLY(count iter;)
  DIV_ONLY(int phase SHM_ONLY(, shmid);)
} Slice;

workerini cubaini;

#if defined(HAVE_SHMGET) && (defined(SUAVE) || defined(DIVONNE))
#define FRAMECOPY
#endif

#ifdef DEBUG
#define TERM_RED "\e[31m"
#define TERM_BLUE "\e[34m"
#define TERM_RESET "\e[0m\n"
#define MASTER(s, ...) \
fprintf(stderr, TERM_RED ROUTINE " master: " s TERM_RESET, __VA_ARGS__)
#define WORKER(s, ...) \
fprintf(stderr, TERM_BLUE ROUTINE " worker: " s TERM_RESET, __VA_ARGS__)
#else
#define MASTER(s, ...)
#define WORKER(s, ...)
#endif

/*********************************************************************/

#ifndef MSG_WAITALL
/* Windows */
#define MSG_WAITALL 0
#endif

static inline int readsock(int fd, void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = recv(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

static inline int writesock(int fd, const void *data, size_t n)
{
  ssize_t got;
  size_t remain = n;
  do got = send(fd, data, remain, MSG_WAITALL);
  while( got > 0 && (data += got, remain -= got) > 0 );
  return got;
}

/*********************************************************************/

static void DoSample(This *t, number n, creal *x, real *f
  VES_ONLY(, creal *w, ccount iter))
{
  cint ncores = IMin(t->ncores, n/MINSLICE);

  if( ncores < MINCORES ) DoSampleSerial(t, n, x, f VES_ONLY(, w, iter));
  else {
    Slice slice;
    int core, abort;
    char s[128];

    t->neval += n;

    slice.m = slice.n = (n + ncores - 1)/ncores;
    if( VERBOSE > 2 ) {
      sprintf(s, "sampling " NUMBER " points each on %d cores",
        slice.n, ncores);
      Print(s);
    }

    slice.i = 0;
    VES_ONLY(slice.iter = iter;)
    DIV_ONLY(slice.phase = t->phase;)

#ifdef DIVONNE
    if( n > t->nframe ) {
      FrameFree(t, ShmRm(t));
      t->nframe = n;
      FrameAlloc(t);
    }
    SHM_ONLY(slice.shmid = t->shmid;)
#endif

    SHM_ONLY(if( t->shmid != -1 ) {
      slice.m = n;
#ifdef FRAMECOPY
      VES_ONLY(Copy(t->frame, w, n);)
      Copy(t->frame + n*NW, x, n*t->ndim);
#endif
    })

    for( core = 0; core < ncores; ++core ) {
      cint fd = t->child[core];
      MASTER("sending " NUMBER " samples to core %d fd %d",
        slice.n, core, fd);
      writesock(fd, &slice, sizeof slice);
      SHM_ONLY(if( t->shmid == -1 )) {
        VES_ONLY(writesock(fd, w, slice.n*sizeof *w);
                 w += slice.n;)
        writesock(fd, x, slice.n*t->ndim*sizeof *x);
        x += slice.n*t->ndim;
      }
      slice.i += slice.n;
      n -= slice.n;
      slice.n = IMin(slice.n, n);
    }

    abort = 0;
    for( core = ncores; --core >= 0; ) {
      cint fd = t->child[core];
      MASTER("reading from core %d fd %d", core, fd);
      readsock(fd, &slice, sizeof slice);
      MASTER("reading " NUMBER " samples from core %d fd %d",
        slice.n, core, fd);
      if( slice.n == -1 ) abort = 1;
      else SHM_ONLY(if( t->shmid == -1 )) readsock(fd,
        f + slice.i*t->ncomp, slice.n*t->ncomp*sizeof *f);
    }
    if( abort ) longjmp(t->abort, -99);

#ifdef FRAMECOPY
    if( t->shmid != -1 )
      Copy(f, t->frame + slice.m*(NW + t->ndim), slice.m*t->ncomp);
#endif
  }
}

/*********************************************************************/

#ifdef DIVONNE

static inline int ReadyCore(cThis *t)
{
  int core;
  fd_set ready;

  memcpy(&ready, &t->children, sizeof ready);
  select(t->nchildren, &ready, NULL, NULL, NULL);

  for( core = 0; core < t->ncores; ++core )
    if( FD_ISSET(t->child[core], &ready) ) break;

  return core;
}

/*********************************************************************/

typedef struct {
  number neval, neval_opt, neval_cut;
  count nregions, iregion, retval;
} ExploreResult;

static int Explore(This *t, cint iregion)
{
  TYPEDEFREGION;
  Region *region;
  int ireg = iregion, core = t->running;

  if( t->ncores < MINCORES ) return ExploreSerial(t, iregion);

  if( t->running >= ((iregion < 0) ? 1 : t->ncores) ) {
    Totals totals[t->ncomp];
    cint fd = t->child[core = ReadyCore(t)];
    ExploreResult res;
    count comp, succ;

    --t->running;
    readsock(fd, &res, sizeof res);
    ireg = res.iregion;
    region = RegionPtr(ireg);
    succ = ireg + region->next;
    readsock(fd, region, sizeof(Region));
    if( --res.nregions > 0 ) {
      region->next = t->nregions - ireg;
      EnlargeRegions(t, res.nregions);
      readsock(fd, RegionPtr(t->nregions), res.nregions*sizeof(Region));
      t->nregions += res.nregions;
      RegionPtr(t->nregions-1)->next = succ - t->nregions + 1;
    }

    readsock(fd, totals, sizeof totals);
    for( comp = 0; comp < t->ncomp; ++comp )
      t->totals[comp].secondspread =
        Max(t->totals[comp].secondspread, totals[comp].secondspread);

    t->neval += res.neval;
    t->neval_opt += res.neval_opt;
    t->neval_cut += res.neval_cut;

    if( res.retval == -1 ) return -1;
  }

  if( iregion >= 0 ) {
    Slice slice;
    cint fd = t->child[core];
    slice.n = 0;
    slice.i = iregion;
    slice.phase = t->phase;
    region = RegionPtr(iregion);
    writesock(fd, &slice, sizeof slice);
    writesock(fd, &t->samples[region->isamples], sizeof(Samples));
    writesock(fd, region, sizeof *region);
    writesock(fd, t->totals, sizeof *t->totals);
    region->depth = 0;
    ++t->running;
  }

  return ireg;
}
#endif

/*********************************************************************/

static void DoChild(This *t, cint fd)
{
  Slice slice;

#ifdef DIVONNE
  TYPEDEFREGION;
  Totals totals[t->ncomp];
  ExploreResult res;

  t->totals = totals;
  t->ncores = 0;	/* no recursive forks */
  AllocRegions(t);
  SamplesIni(&t->samples[0]);
  t->samples[0].n = 0;
  SamplesIni(&t->samples[1]);
  t->samples[1].n = 0;
  SamplesIni(&t->samples[2]);
  t->samples[2].n = 0;
#endif

#ifdef SUAVE
  SHM_ONLY(if( t->shmid == -1 ))
    MemAlloc(t->frame, t->nframe*SAMPLESIZE);
#endif

  if( cubaini.initfun ) cubaini.initfun(cubaini.initarg);

  while( readsock(fd, &slice, sizeof slice) ) {
    number n = slice.n;
    DIV_ONLY(t->phase = slice.phase;)
    if( n > 0 ) {
      real VES_ONLY(*w,) *x, *f;
      WORKER("read " NUMBER " samples from fd %d", n, fd);

#ifdef DIVONNE
      if( n > t->nframe ) {
        FrameFree(t);
        t->nframe = n;
        SHM_ONLY(t->shmid = slice.shmid; ShmMap(t) else)
        MemAlloc(t->frame, t->nframe*SAMPLESIZE);
      }
#endif

      VES_ONLY(w = t->frame;)
      x = t->frame + slice.m*NW;
      f = x + slice.m*t->ndim;

      SHM_ONLY(if( t->shmid != -1 ) {
        VES_ONLY(w += slice.i;)
        x += slice.i*t->ndim;
        f += slice.i*t->ncomp;
      }
      else) {
        VES_ONLY(readsock(fd, w, n*sizeof *w);)
        readsock(fd, x, n*t->ndim*sizeof *x);
      }

      slice.n |= SampleRaw(t, n, x, f VES_ONLY(, w, slice.iter));
      WORKER("writing " NUMBER " samples to fd %d", n, fd);
      writesock(fd, &slice, sizeof slice);
      if( SHM_ONLY(t->shmid == -1 &&) slice.n != -1 )
        writesock(fd, f, slice.n*t->ncomp*sizeof *f);
    }
#ifdef DIVONNE
    else {
      Samples *samples, psamples;

      readsock(fd, &psamples, sizeof psamples);
      readsock(fd, RegionPtr(0), sizeof(Region));
      readsock(fd, totals, sizeof totals);
      t->nregions = 1;
      t->neval = t->neval_opt = t->neval_cut = 0;

      WORKER("read 1 region from fd %d", fd);

      samples = &t->samples[RegionPtr(0)->isamples];
      if( psamples.n != samples->n ) {
        SamplesFree(samples);
        *samples = psamples;
        SamplesAlloc(t, samples);
      }

      res.retval = ExploreSerial(t, 0);
      res.neval = t->neval;
      res.neval_opt = t->neval_opt;
      res.neval_cut = t->neval_cut;
      res.nregions = t->nregions;
      res.iregion = slice.i;
      WORKER("writing %d regions to fd %d", res.nregions, fd);
      writesock(fd, &res, sizeof res);
      writesock(fd, RegionPtr(0), t->nregions*sizeof(Region));
      writesock(fd, totals, sizeof totals);
    }
#endif
  }

  if( cubaini.exitfun ) cubaini.exitfun(cubaini.exitarg);

  exit(0);
}

/*********************************************************************/

#ifdef HAVE_GETLOADAVG
double cubaloadavg_;
#endif

static inline void ForkCores(This *t)
{
  int core;
  char s[128];
  cchar *env = getenv("CUBACORES");

  t->ncores = env ? atoi(env) : sysconf(_SC_NPROCESSORS_ONLN);
#ifdef HAVE_GETLOADAVG
  if( env == NULL || t->ncores < 0 ) {
    if( cubaloadavg_ < 0 ) getloadavg(&cubaloadavg_, 1);
    t->ncores = abs(t->ncores) - floor(cubaloadavg_);
  }
#endif

  DIV_ONLY(t->nchildren = t->running = 0;)

  if( t->ncores < MINCORES ) return;
  if( VERBOSE ) {
    sprintf(s, "using %d cores via "
#ifdef HAVE_SHMGET
      "shared memory",
#else
      "pipes",
#endif
      t->ncores);
    Print(s);
  }

  fflush(NULL);		/* make sure all buffers are flushed,
			   or else buffered content will be written
			   out multiply, at each child's exit(0) */

  Alloc(t->child, t->ncores);
  for( core = 0; core < t->ncores; ++core ) {
    int fd[2];
    pid_t pid;
    assert(
      socketpair(AF_LOCAL, SOCK_STREAM, 0, fd) != -1 &&
      (pid = fork()) != -1 );
    if( pid == 0 ) {
      close(fd[0]);
      DoChild(t, fd[1]);
    }
    MASTER("forked core %d pid %d pipe %d(master) -> %d(worker)",
      core, pid, fd[0], fd[1]);
    close(fd[1]);
    t->child[core] = fd[0];
    DIV_ONLY(FD_SET(fd[0], &t->children);
             t->nchildren = IMax(t->nchildren, fd[0] + 1);)
  }
}

/*********************************************************************/

static inline void WaitCores(cThis *t)
{
  if( t->ncores >= MINCORES ) {
    int core;
    pid_t pid;
    for( core = 0; core < t->ncores; ++core ) {
      MASTER("closing core %d fd %d", core, t->child[core]);
      close(t->child[core]);
    }
    free(t->child);
    for( core = 0; core < t->ncores; ++core ) {
      MASTER("waiting for core %d", core);
      wait(&pid);
      MASTER("core %d pid %d terminated", core, pid);
    }
  }
}

