/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Memory Handling Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egadsTypes.h"

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __DEVICE__
#undef __DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#define __DEVICE__ extern "C" __device__
#else
#define __HOST_AND_DEVICE__
#define __DEVICE__
#endif


__HOST_AND_DEVICE__ /*@null@*/ /*@out@*/ /*@only@*/ void *
EG_alloc(size_t nbytes)
{
  if (nbytes == 0) return NULL;
  return malloc(nbytes);
}


__HOST_AND_DEVICE__ /*@null@*/ /*@only@*/ void *
EG_calloc(size_t nele, size_t size)
{
#ifndef __CUDA_ARCH__
  if (nele*size == 0) return NULL;
  return calloc(nele, size);
#else
  void *ptr = EG_alloc(nele*size);
  return memset(ptr, 0, nele*size);
#endif
}


__HOST_AND_DEVICE__ /*@null@*/ /*@only@*/ void *
EG_reall(/*@null@*/ /*@only@*/ /*@returned@*/ void *ptr, size_t nbytes)
{
#ifndef __CUDA_ARCH__
  if (nbytes == 0) return NULL;
  return realloc(ptr, nbytes);
#else
  return NULL;
#endif
}


__HOST_AND_DEVICE__ void
EG_free(/*@null@*/ /*@only@*/ void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}


__HOST_AND_DEVICE__ /*@null@*/ /*@only@*/ char *
EG_strdup(/*@null@*/ const char *str)
{
#ifndef __CUDA_ARCH__
  int  i, len;
  char *dup;

  if (str == NULL) return NULL;

  len = strlen(str) + 1;
  dup = (char *) EG_alloc(len*sizeof(char));
  if (dup != NULL)
    for (i = 0; i < len; i++) dup[i] = str[i];

  return dup;
#else
  return NULL;
#endif
}

