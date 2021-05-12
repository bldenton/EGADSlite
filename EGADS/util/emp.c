/*
 *      EMP: Explicit Multithread Package
 *
 *      Copyright 2013-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif


#ifdef WIN32

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <process.h>


#ifndef __CUDA_ARCH__
static long
EMP_getseconds()
{
  SYSTEMTIME    systime;
  FILETIME      filtime;
  LARGE_INTEGER large;
  DWORDLONG     int64;

  GetSystemTime(&systime);
  SystemTimeToFileTime(&systime, &filtime);
  large.LowPart  = filtime.dwLowDateTime;
  large.HighPart = filtime.dwHighDateTime;
  int64          = large.QuadPart/10000000;
  return int64 & 0x7FFFFFFF;
}
#endif


/* initialize the EMP block */

__HOST_AND_DEVICE__
int EMP_Init(long *start)
{
#ifndef __CUDA_ARCH__
  int         np;
  char        *env;
  SYSTEM_INFO siSysInfo;

  GetSystemInfo(&siSysInfo);
  np = siSysInfo.dwNumberOfProcessors;

  env = getenv("EMPnumProc");
  if (env != NULL) {
    np = atoi(env);
    if (np < 1) np = 1;
  }
  if (start != NULL) *start = EMP_getseconds();
  return np;
#else
  return 1;
#endif
}


/* close the EMP block */

__HOST_AND_DEVICE__
long EMP_Done(long *start)
{
#ifndef __CUDA_ARCH__
  if (start == NULL) return 0L;
  return EMP_getseconds() - *start;
#else
  return 0L;
#endif
}


/* Waste a little time */

__HOST_AND_DEVICE__
void EMP_ThreadSpin()
{
#ifndef __CUDA_ARCH__
  Sleep(1);
#endif
}


/* Get the current Thread ID */

__HOST_AND_DEVICE__
long EMP_ThreadID()
{
#ifndef __CUDA_ARCH__
  return GetCurrentThreadId();
#else
  return 0L;
#endif
}


/* Spawn off another thread */

__HOST_AND_DEVICE__
HANDLE *EMP_ThreadCreate(void (*entry)(void *), void *arg)
{
#ifndef __CUDA_ARCH__
  HANDLE   *thread;
  unsigned  threadID;

  thread = (HANDLE *) malloc(sizeof(HANDLE));
  if (thread == NULL) return NULL;

  *thread = (HANDLE) _beginthreadex(NULL, 0, 
				    (unsigned int (__stdcall *)(void *)) entry,
                                    arg,  0, &threadID);
  return thread;
#else
  return NULL;
#endif
}


/* Ends the thread execution with the return code */

__HOST_AND_DEVICE__
void EMP_ThreadExit()
{
#ifndef __CUDA_ARCH__
  _endthreadex(0);
#endif
}


/* Wait for the thread to finish */

__HOST_AND_DEVICE__
void EMP_ThreadWait(HANDLE *thread)
{
#ifndef __CUDA_ARCH__
  if (WaitForSingleObject(*thread, INFINITE) == WAIT_FAILED)
    printf(" Warning: ThreadWait FAILED!\n");
#endif
}


/* Destroy the thread memory */

__HOST_AND_DEVICE__
void EMP_ThreadDestroy(HANDLE *thread)
{
#ifndef __CUDA_ARCH__
  CloseHandle(*thread);
  free(thread);
#endif
}


/* Create a lock (unlocked) */

__HOST_AND_DEVICE__
HANDLE *EMP_LockCreate()
{
#ifndef __CUDA_ARCH__
  HANDLE *mutex;

  mutex  = (HANDLE *) malloc(sizeof(HANDLE));
  if (mutex == NULL) return NULL;
  *mutex = CreateMutex(NULL, FALSE, NULL);
  if (*mutex == NULL) {
    printf(" ERROR: MUTEX not assigned (LockCreate)!\n");
    free(mutex);
    return NULL;
  }
  return mutex;
#else
  return NULL;
#endif
}


/* Destroy the lock memory */

__HOST_AND_DEVICE__
void EMP_LockDestroy(HANDLE *mutex)
{
#ifndef __CUDA_ARCH__
  ReleaseMutex(*mutex);
  free(mutex);
#endif
}


/* Set a lock (wait if already set) */

__HOST_AND_DEVICE__
void EMP_LockSet(HANDLE *mutex)
{
#ifndef __CUDA_ARCH__
  if (WaitForSingleObject(*mutex, INFINITE) == WAIT_FAILED)
    printf(" Warning: LockSet Wait FAILED!\n");
#endif
}


/* Gets the value of a lock (0-unset, 1-set) */

__HOST_AND_DEVICE__
int EMP_LockTest(HANDLE *mutex)
{
#ifndef __CUDA_ARCH__
  DWORD stat;

  stat = WaitForSingleObject(*mutex, 0L);
  if (stat == WAIT_FAILED) {
    return 1;
  } else {
    ReleaseMutex(*mutex);
    return 0;
  }
#else
  return 0;
#endif
}


/* Release a lock */

__HOST_AND_DEVICE__
void EMP_LockRelease(HANDLE *mutex)
{
#ifndef __CUDA_ARCH__
  if (ReleaseMutex(*mutex) == 0)
    printf(" Warning: LockRelease Unlock FAILED!\n");
#endif
}


#else


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <sys/time.h>


#ifndef __CUDA_ARCH__
static long
EMP_getseconds()
{
  struct timeval tv;

  gettimeofday(&tv, NULL);
  return tv.tv_sec;
}
#endif


/* Initialize the EMP block */

__HOST_AND_DEVICE__
int EMP_Init(/*@null@*/ long *start)
{
#ifndef __CUDA_ARCH__
  char *env;
  long nprocs = 1;

#ifdef _SC_NPROCESSORS_ONLN
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);
  if (nprocs < 1)
    printf("Could not determine number of CPUs online:\n%s\n", strerror(errno));
#else
  printf(stderr, "Could not determine number of CPUs");
#endif

  env = getenv("EMPnumProc");
  if (env != NULL) {
    nprocs = atoi(env);
    if (nprocs < 1) nprocs = 1;
  }
  if (start != NULL) *start = EMP_getseconds();
  return nprocs;
#else
  return 1;
#endif
}


/* Close the EMP block */

__HOST_AND_DEVICE__
long EMP_Done(long *start)
{
#ifndef __CUDA_ARCH__
  if (start == NULL) return 0L;
  return EMP_getseconds() - *start;
#else
  return 0L;
#endif
}


/* Waste a little time -- yeild */

__HOST_AND_DEVICE__
void EMP_ThreadSpin()
{
#ifndef __CUDA_ARCH__
  usleep(1000);
#endif
}


/* Get the current Thread ID */

__HOST_AND_DEVICE__
long EMP_ThreadID()
{
#ifndef __CUDA_ARCH__
  return (long) pthread_self();
#else
  return 0L;
#endif
}


/* Spawn off another thread */

__HOST_AND_DEVICE__
/*@null@*/ void *EMP_ThreadCreate(void (*entry)(void *), /*@null@*/ void *arg)
{
#ifndef __CUDA_ARCH__
  int            stat;
  pthread_t      *thread;
  pthread_attr_t attr;

  thread = (pthread_t *) malloc(sizeof(pthread_t));
  if (thread == NULL) return NULL;

  pthread_attr_init(&attr);
  stat = pthread_create(thread, &attr, (void * (*) (void *)) entry, arg);
  if (stat != 0) {
    free(thread);
    return NULL;
  }

  return (void *) thread;
#else
  return NULL;
#endif
}


/* Ends the thread execution */

__HOST_AND_DEVICE__
void EMP_ThreadExit()
{
#ifndef __CUDA_ARCH__
  pthread_exit(NULL);
#endif
}


/* Wait for the thread to finish */

__HOST_AND_DEVICE__
void EMP_ThreadWait(void *vthread)
{
#ifndef __CUDA_ARCH__
  int       stat;
  pthread_t *thread;

  thread = (pthread_t *) vthread;
  stat   = pthread_join(*thread, NULL);
  if (stat != 0)
    printf(" Threading ERROR: %d (ThreadWait)\n", stat);
#endif
}


/* Destroy the thread memory */

__HOST_AND_DEVICE__
void EMP_ThreadDestroy(/*@only@*/ void *vthread)
{
#ifndef __CUDA_ARCH__
#if defined(DARWIN) || defined(DARWIN64)
  pthread_t *thread;

  thread = (pthread_t *) vthread;
  pthread_detach(*thread);
#endif
  free(vthread);
#endif
}


/* Create a lock (unlocked) */

__HOST_AND_DEVICE__
/*@null@*/ void *EMP_LockCreate()
{
#ifndef __CUDA_ARCH__
  int             stat;
  pthread_mutex_t *mutex;

  mutex = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t));
  if (mutex != NULL) {
    stat = pthread_mutex_init(mutex, NULL);
    if (stat != 0) {
      printf(" Threading ERROR: %d (LockCreate)\n", stat);
      free(mutex);
      return NULL;
    }
  }
  return (void *) mutex;
#else
  return (void *) NULL;
#endif
}


/* Destroy the lock memory */

__HOST_AND_DEVICE__
void EMP_LockDestroy(/*@only@*/ void *vlock)
{
#ifndef __CUDA_ARCH__
  pthread_mutex_t *lock;

  lock = (pthread_mutex_t *) vlock;
  pthread_mutex_destroy(lock);
  free(lock);
#endif
}


/* Set a lock (wait if already set) */

__HOST_AND_DEVICE__
void EMP_LockSet(void *vlock)
{
#ifndef __CUDA_ARCH__
  int             stat;
  pthread_mutex_t *lock;

  lock = (pthread_mutex_t *) vlock;
  stat = pthread_mutex_lock(lock);
  if (stat != 0) printf(" Threading Warning: %d (LockSet)\n", stat);
#endif
}


/* Gets the value of a lock (0-unset, 1-set) */

__HOST_AND_DEVICE__
int EMP_LockTest(void *vlock)
{
#ifndef __CUDA_ARCH__
  int             stat;
  pthread_mutex_t *lock;

  lock = (pthread_mutex_t *) vlock;
  stat = pthread_mutex_trylock(lock);
  if (stat == EBUSY) {
    return 1;
  } else if (stat == 0) {
    pthread_mutex_unlock(lock);
  } else {
    printf(" Fatal Threading ERROR: %d (LockTest)\n", stat);
    exit(EXIT_FAILURE);
  }
#endif
  return 0;
}


/* Release a lock */

__HOST_AND_DEVICE__
void EMP_LockRelease(void *vlock)
{
#ifndef __CUDA_ARCH__
  pthread_mutex_t *lock;

  lock = (pthread_mutex_t *) vlock;
  pthread_mutex_unlock(lock);
#endif
}
#endif


/* structure to hold control info for EMP_for, EMP_sum, and EMP_min */

typedef struct {
  int    nindex;                    /* number of indices */
  int    index;                     /* current loop index */
  int    status;                    /* status flag */
  long   masterID;                  /* master thread ID */
  void   *mutex;                    /* the mutex or NULL for single thread */
  int    (*forFn)(int);             /* routine that does work in EMP_for */
  int    (*sumFn)(int, double *);   /* routine that does work in EMP_sum */
  double sum;                       /* global sum */
  int    (*minFn)(int, double *);   /* routine that does work for EMP_min */
  double min;                       /* global minimum */
  int    imin;                      /* index associated with minimum (or -1) */
} emp_T;


/* Inner routine used by EMP_for */

static void EMP_for_inner(void *Global)
{
  int   index, status;
  long  ID;
  emp_T *global;

  global = (emp_T *)Global;

  /* get our identifier */
  ID = EMP_ThreadID();

  /* loop as long as work remains and status==0 */
  while (1) {

    /* only one thread at a time here - controlled by a mutex */
    if (global->mutex != NULL ) EMP_LockSet(global->mutex);

    status        = global->status;
    index         = global->index;
    global->index = index + 1;

    if (global->mutex != NULL) EMP_LockRelease(global->mutex);

    /* check that status is still 0 */
    if (status != 0) break;

    /* check that there is really work to be done */
    if (index >= global->nindex) break;

    /* do the work */
    status = global->forFn(index);

    /* if an error was found, save status (inside a mutex) and return */
    if (status != 0) {
      if (global->mutex != NULL) EMP_LockSet(global->mutex);

      global->status = status;

      if (global->mutex != NULL) EMP_LockRelease(global->mutex);

      break;
    }
  }

  /* finished all work, so exit */
  if (ID != global->masterID) EMP_ThreadExit();
}


/* Simple parallel "for" loop */

int EMP_for(int maxproc, int nindex, int (*forFn)(int index))
{
  int   i, np;
  long  start;
  void  **threads = NULL;
  emp_T global;

  /* set up structure to hold parallelization info */
  global.nindex   = nindex;
  global.index    = 0;
  global.status   = 0;

  global.masterID = EMP_ThreadID();
  global.mutex    = NULL;

  global.forFn    = forFn;

  global.sumFn    = NULL;
  global.sum      = 0;

  global.minFn    = NULL;
  global.min      = 0;
  global.imin     = -1;

  /* np is limited by user-given number of processors */
  np = EMP_Init(&start);
  if (np > maxproc) np = maxproc;

  /* 
   * create the mutex to handle list synchronization and
   * get storage for the extra threads 
   */
  if (np > 1) {
    global.mutex = EMP_LockCreate();
    if (global.mutex == NULL) {
      np = 1;
    } else {
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(global.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and execute the inner function in them */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EMP_for_inner, &global);
      if (threads[i] == NULL) printf(" Error creating thread %d\n", i+1);
    }

  /* now run the inner routine from the original thread */
  EMP_for_inner(&global);

  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* cleanup */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);

  if (global.mutex != NULL) EMP_LockDestroy(global.mutex);
  if (threads      != NULL) free(threads);

  /* close the EMP package (for now) */
  EMP_Done(&start);

  return global.status;
}


/* Inner routine used by EMP_sum */

static void EMP_sum_inner(void *Global)
{
  int    index, status;
  long   ID;
  double sum;
  emp_T  *global;

  global = (emp_T *)Global;

  /* get our identifier */
  ID = EMP_ThreadID();

  /* loop as long as work remains and status==0 */
  while (1) {

    /* only one thread at a time here - controlled by a mutex */
    if (global->mutex != NULL) EMP_LockSet(global->mutex);

    status        = global->status;
    index         = global->index;
    global->index = index + 1;

    if (global->mutex != NULL) EMP_LockRelease(global->mutex);

    /* check that status is still 0 */
    if (status != 0) break;

    /* check that there is really work to be done */
    if (index >= global->nindex) break;

    /* do the work */
    status = global->sumFn(index, &sum);

    /* if an error was found, save status (inside a mutex) and return */
    if (status != 0) {
      if (global->mutex != NULL) EMP_LockSet(global->mutex);

      global->status = status;

      if (global->mutex != NULL) EMP_LockRelease(global->mutex);

      break;
    }

    /* add the sum to the global sum (inside a mutex) */
    if (global->mutex != NULL) EMP_LockSet(global->mutex);

    global->sum += sum;

    if (global->mutex != NULL) EMP_LockRelease(global->mutex);
  }

  /* finished all work, so exit */
  if (ID != global->masterID) EMP_ThreadExit();

}

/* Parallel loop to find a sum */

int EMP_sum(int maxproc, int nindex, int (*sumFn)(int index, double *sum),
            double *sum)
{
  int   i, np;
  long  start;
  void  **threads = NULL;
  emp_T global;

  /* set up structure to hold parallelization info */
  global.nindex   = nindex;
  global.index    = 0;
  global.status   = 0;

  global.masterID = EMP_ThreadID();
  global.mutex    = NULL;

  global.forFn    = NULL;

  global.sumFn    = sumFn;
  global.sum      = 0;

  global.minFn    = NULL;
  global.min      = 0;
  global.imin     = -1;

  /* np is limited by user-given number of processors */
  np = EMP_Init(&start);
  if (np > maxproc) np = maxproc;

  /* 
   * create the mutex to handle list synchronization and
   * get storage for the extra threads 
   */
  if (np > 1) {
    global.mutex = EMP_LockCreate();
    if (global.mutex == NULL) {
      np = 1;
    } else {
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(global.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and execute the inner routine in them */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EMP_sum_inner, &global);
      if (threads[i] == NULL) printf(" Error creating thread %d\n", i+1);
    }

  /* now run the inner routine from the original thread */
  EMP_sum_inner(&global);

  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* cleanup */
  if (threads != NULL) 
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);

  if (global.mutex != NULL) EMP_LockDestroy(global.mutex);
  if (threads != NULL) free(threads);

  /* close the EMP package (for now) */
  EMP_Done(&start);

  *sum = global.sum;
  return global.status;
}


#ifndef __clang_analyzer__
/* Inner routine used by EMP_min */

static void EMP_min_inner(void *Global)
{
  int    index, status;
  long   ID;
  double min;
  emp_T  *global;

  global = (emp_T *)Global;

  /* get our identifier */
  ID = EMP_ThreadID();

  /* loop as long as work remains and status==0 */
  while (1) {

    /* only one thread at a time here - controlled by a mutex */
    if (global->mutex != NULL) EMP_LockSet(global->mutex);

    status        = global->status;
    index         = global->index;
    min           = global->min;
    global->index = index + 1;

    if (global->mutex != NULL) EMP_LockRelease(global->mutex);

    /* check that status is still 0 */
    if (status != 0) break;

    /* check that there is really work to be done */
    if (index >= global->nindex) break;

    /* do the work */
    status = global->minFn(index, &min);

    /* if an error was found, save status (inside a mutex) and return */
    if (status != 0) {
      if (global->mutex != NULL) EMP_LockSet(global->mutex);

      global->status = status;

      if (global->mutex != NULL) EMP_LockRelease(global->mutex);

      break;
    }

    /* remember this min if small than global min (inside a mutex) */
    if (global->mutex != NULL) EMP_LockSet(global->mutex);

    if (min < global->min) {
        global->min  = min;
        global->imin = index;
    }

    if (global->mutex != NULL) EMP_LockRelease(global->mutex);
  }

  /* finished all work, so exit */
  if (ID != global->masterID) EMP_ThreadExit();

}


/* Parallel loop to find a minimum */

int EMP_min(int maxproc, int nindex, int (*minFn)(int index, double *min),
            double *min, int *imin)
{
  int   i, np;
  long  start;
  void  **threads = NULL;
  emp_T global;

  /* set up structure to hold parallelization info */
  global.nindex   = nindex;
  global.index    = 0;
  global.status   = 0;

  global.masterID = EMP_ThreadID();
  global.mutex    = NULL;

  global.forFn    = NULL;

  global.sumFn    = NULL;
  global.sum      = 0;

  global.minFn    = minFn;
  global.min      = *min;
  global.imin     = -1;

  /* np is limited by user-given number of processors */
  np = EMP_Init(&start);
  if (np > maxproc) np = maxproc;

  /* 
   * create the mutex to handle list synchronization and
   * get storage for the extra threads 
   */
  if (np > 1) {
    global.mutex = EMP_LockCreate();
    if (global.mutex == NULL) {
      np = 1;
    } else {
      threads = (void **) malloc((np-1)*sizeof(void *));
      if (threads == NULL) {
        EMP_LockDestroy(global.mutex);
        np = 1;
      }
    }
  }

  /* create the threads and execute the inner routine in them */
  if (threads != NULL)
    for (i = 0; i < np-1; i++) {
      threads[i] = EMP_ThreadCreate(EMP_min_inner, &global);
      if (threads[i] == NULL) printf(" Error creating thread %d\n", i+1);
    }

  /* now run the inner routine from the original thread */
  EMP_min_inner(&global);

  /* wait for all others to return */
  if (threads != NULL)
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadWait(threads[i]);

  /* cleanup */
  if (threads != NULL) 
    for (i = 0; i < np-1; i++)
      if (threads[i] != NULL) EMP_ThreadDestroy(threads[i]);

  if (global.mutex != NULL) EMP_LockDestroy(global.mutex);
  if (threads != NULL) free(threads);

  /* close the EMP package (for now) */
  EMP_Done(&start);

  *min  = global.min;
  *imin = global.imin;
  return global.status;
}
#endif
