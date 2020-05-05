#ifndef CUDA_UTIL_H
#define CUDA_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __NVCC__
# include <cuda.h>
# include <cuda_runtime.h>
#endif

#ifndef __CUDA_SAFE_CALL
#define __CUDA_SAFE_CALL(call)                                                 \
do {                                                                           \
  cudaError_t cuda_error = call;                                               \
  if( cudaSuccess != cuda_error ) {                                            \
    fprintf(stderr,"CUDA Error: %s, in '%s', line %i\n",                       \
            cudaGetErrorString(cuda_error),__FILE__,__LINE__);                 \
  }                                                                            \
} while(0)
#endif

#ifndef __CUDA_ERROR_CHECK
#define __CUDA_ERROR_CHECK                                                     \
do {                                                                           \
  cudaError err = cudaGetLastError();                                          \
  if ( cudaSuccess != err ) {                                                  \
    printf("cudaCheckError() failed at %s:%i : %s\n",__FILE__,__LINE__,        \
           cudaGetErrorString( err ) );                                        \
  }                                                                            \
                                                                               \
  /* More careful checking. However, this will affect performance. */          \
  /* Comment away if needed.                                       */          \
  err = cudaDeviceSynchronize();                                               \
  if( cudaSuccess != err ) {                                                   \
    fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",      \
             __FILE__, __LINE__, cudaGetErrorString( err ) );                  \
  }                                                                            \
} while(0)
#endif

#ifdef __NVCC__
# define __HOST_AND_DEVICE__ __host__ __device__
# define __DEVICE__ __device__
#else
# define __HOST_AND_DEVICE__
# define __DEVICE__
#endif

/**
 * The following macros are for convenience to perform general CUDA operations
 **/

#if defined(__NVCC__) && !defined(__CUDA_ARCH__)

# define DEVICE_MALLOC(devptraddr, nmemb)                                      \
do {                                                                           \
  *(devptraddr) = NULL;                                                        \
  __CUDA_SAFE_CALL(cudaMalloc((void**)(devptraddr), nmemb));                   \
} while(0)

# define DEVICE_REALLOC(devdstaddr, devsrc, type, old_nmemb, nmemb)            \
do {                                                                           \
  type *temp;                                                                  \
  DEVICE_MALLOC((void **)&temp, (nmemb)*sizeof(type));                         \
  if (temp != NULL) {                                                          \
    MEMCPY_DEVICE_TO_DEVICE(temp, devsrc, (old_nmemb)*sizeof(type));           \
    DEVICE_FREE(devsrc);                                                       \
  }                                                                            \
  *(devdstaddr) = temp;                                                        \
} while(0)

# define DEVICE_FREE(devptr) __CUDA_SAFE_CALL(cudaFree(devptr))

# define MEMCPY_HOST_TO_DEVICE(devptr, hostptr, nmemb) __CUDA_SAFE_CALL(cudaMemcpy(devptr, hostptr, nmemb, cudaMemcpyHostToDevice))

# define MEMCPY_DEVICE_TO_HOST(hostptr, devptr, nmemb) __CUDA_SAFE_CALL(cudaMemcpy(hostptr, devptr, nmemb, cudaMemcpyDeviceToHost))

# define MEMCPY_DEVICE_TO_DEVICE(dst, src, nmemb) __CUDA_SAFE_CALL(cudaMemcpy(dst, src, nmemb, cudaMemcpyDeviceToDevice))

# define MEMCPY_HOST_TO_HOST(dst, src, nmemb) __CUDA_SAFE_CALL(cudaMemcpy(dst, src, nmemb, cudaMemcpyHostToHost))

# define STRDUP_HOST_TO_DEVICE(devptraddr, hostptr)                            \
do {                                                                           \
  char*  string_d;                                                             \
  size_t strcnt = strlen(hostptr) + 1;                                         \
  DEVICE_MALLOC((void **)&string_d, strcnt*sizeof(char));                      \
  MEMCPY_HOST_TO_HOST(devptraddr, &string_d, sizeof(char*));                   \
  if ( string_d != NULL ) {                                                    \
    MEMCPY_HOST_TO_DEVICE(string_d, hostptr, strcnt*sizeof(char));             \
  }                                                                            \
} while(0)

# define STRLEN_DEVICE(devptraddr, N)                                          \
do {                                                                           \
  size_t i=0;                                                                  \
  char*  string_d;                                                             \
  char   c;                                                                    \
  MEMCPY_DEVICE_TO_HOST(&string_d, devptraddr, sizeof(char*));                 \
  do {                                                                         \
    MEMCPY_DEVICE_TO_HOST(&(c), &(string_d[i]), sizeof(char));                 \
  } while(c != '\0' && ++i);                                                   \
  N = i;                                                                       \
} while(0)

# define STRNDUP_HOST_TO_DEVICE(devptraddr, hostptr, N)                        \
do {                                                                           \
  char*  string_d;                                                             \
  DEVICE_MALLOC((void **)&string_d, N*sizeof(char));                           \
  MEMCPY_HOST_TO_HOST(devptraddr, &string_d, sizeof(char*));                   \
  if ( string_d != NULL ) {                                                    \
    MEMCPY_HOST_TO_DEVICE(string_d, hostptr, N*sizeof(char));                  \
  }                                                                            \
} while(0)

# define STRNCPY_DEVICE_TO_HOST(hostptr, devptraddr, N)                        \
do {                                                                           \
  char*  string_d=*devptraddr;                                                 \
  size_t i=0;                                                                  \
  do {                                                                         \
    MEMCPY_DEVICE_TO_HOST(&(hostptr[i]), &((*devptraddr)[i]), sizeof(char));   \
  } while(hostptr[i] != '\0' && ++i <= N);                                     \
} while(0)

#else

# define DEVICE_MALLOC(devptraddr, nmemb)  *(devptraddr) = malloc(nmemb)

# define DEVICE_REALLOC(devdstaddr, devsrc, type, old_nmemb, nmemb) *(devdstaddr) = (type*) realloc(devsrc,(nmemb)*sizeof(type))

# define DEVICE_FREE(devptr) free(devptr)

# define MEMCPY_HOST_TO_DEVICE(devptr, hostptr, nmemb) memcpy(devptr, hostptr, nmemb)

# define MEMCPY_DEVICE_TO_HOST(hostptr, devptr, nmemb) memcpy(hostptr, devptr, nmemb)

# define MEMCPY_DEVICE_TO_DEVICE(dst, src, nmemb) memcpy(dst, src, nmemb)

# define STRDUP_HOST_TO_DEVICE(devptraddr, hostptr) *(devptraddr) = strdup(hostptr)

# define STRNDUP_HOST_TO_DEVICE(devptraddr, hostptr, N)                        \
do {                                                                           \
  *(devptraddr) = strndup(hostptr, N);                                         \
  if (*(devptraddr) == NULL) status = EGADS_MALLOC;                            \
} while(0)

# define STRNCPY_DEVICE_TO_HOST(hostptr, devptraddr, N) strncpy(hostptr, *(devptraddr), N)

#endif

#endif
