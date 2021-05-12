/*
 *      EMP: Explicit Multithread Package
 *
 *             Function Header
 *
 *      Copyright 2013-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifndef _EMP_H_

#define _EMP_H_

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __DEVICE__
#undef __DEVICE__
#endif

#ifdef __CUDACC__
# define __HOST_AND_DEVICE__  __host__ __device__
# define __DEVICE__  __device__
#else
# define __HOST_AND_DEVICE__
# define __DEVICE__
#endif

#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus) || defined(WIN32)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef __cplusplus
extern "C" {
#endif

__HOST_AND_DEVICE__
           int   EMP_Init          __ProtoGlarp__(( /*@null@*/ long *start ));
__HOST_AND_DEVICE__
           long  EMP_Done          __ProtoGlarp__(( /*@null@*/ long *start ));

__HOST_AND_DEVICE__
/*@null@*/ void *EMP_ThreadCreate  __ProtoGlarp__(( void (*entry)(void *),
                                                    /*@null@*/ void *arg ));
__HOST_AND_DEVICE__
           void  EMP_ThreadExit    __ProtoGlarp__((  ));
__HOST_AND_DEVICE__
           void  EMP_ThreadWait    __ProtoGlarp__(( void *thread ));
__HOST_AND_DEVICE__
           void  EMP_ThreadSpin    __ProtoGlarp__((  ));
__HOST_AND_DEVICE__
           long  EMP_ThreadID      __ProtoGlarp__((  ));
__HOST_AND_DEVICE__
           void  EMP_ThreadDestroy __ProtoGlarp__(( /*@only@*/ void *thread ));

__HOST_AND_DEVICE__
/*@null@*/ void *EMP_LockCreate    __ProtoGlarp__(( ));
__HOST_AND_DEVICE__
           void  EMP_LockSet       __ProtoGlarp__(( void *lock ));
__HOST_AND_DEVICE__
           int   EMP_LockTest      __ProtoGlarp__(( void *lock ));
__HOST_AND_DEVICE__
           void  EMP_LockRelease   __ProtoGlarp__(( void *lock ));
__HOST_AND_DEVICE__
           void  EMP_LockDestroy   __ProtoGlarp__(( /*@only@*/ void *lock ));

__HOST_AND_DEVICE__
           int   EMP_for           __ProtoGlarp__(( int maxproc, int nindex,
                                                    int (*forFn)(int index) ));
__HOST_AND_DEVICE__
           int   EMP_sum           __ProtoGlarp__(( int maxproc, int nindex,
                                                    int (*sumFn)(int index,
                                                                 double *sum),
                                                    double *sum ));
__HOST_AND_DEVICE__
           int   EMP_min           __ProtoGlarp__(( int maxproc, int nindex,
                                                    int (*minFn)(int index,
                                                                 double *min),
                                                    double *min, int *imin ));

#ifdef __cplusplus
}
#endif

#endif  /*_EMP_H_*/

