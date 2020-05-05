#ifndef _EGADS_STRING_H_
#define _EGADS_STRING_H_

#include <string.h>

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __DEVICE__
#undef __DEVICE__
#endif

#ifdef __CUDACC__
# define __HOST_AND_DEVICE__ __host__ __device__
# define __DEVICE__ __device__
#else
# define __HOST_AND_DEVICE__
# define __DEVICE__
#endif

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__ __HOST_AND_DEVICE__
#else
#define __ProtoExt__ __HOST_AND_DEVICE__ extern
#endif

__ProtoExt__ char *EG_strcpy( char *dest, const char *src );

__ProtoExt__ char *EG_strncpy( char *dest, const char *src, size_t n );

__ProtoExt__ char *EG_strcat( char *dest, const char *src );

__ProtoExt__ int   EG_strncmp( const char *str_a, const char *str_b, size_t n );

#ifdef __cplusplus
}
#endif

#endif

