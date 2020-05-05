#include <stdlib.h>
#include <string.h>

#include "liteString.h"


/*@-mayaliasunique@*/
__HOST_AND_DEVICE__ char *
EG_strcpy(char *dest, const char *src)
{
#ifndef __CUDA_ARCH__
  return strcpy(dest, src);
#else
  int i = 0;
  do {
    dest[i] = src[i];
  } while (src[i++] != '\0');
  return dest;
#endif
}


__HOST_AND_DEVICE__ char *
EG_strncpy(char *dest, const char *src, size_t n)
{
#ifndef __CUDA_ARCH__
  return strncpy(dest, src, n);
#else
  int i = 0;
  do {
    dest[i] = src[i];
  } while (src[i] != '\0' && ++i <= n);
  return dest;
#endif
}


__HOST_AND_DEVICE__ char *
EG_strcat(char *dest, const char *src)
{
#ifndef __CUDA_ARCH__
  return strcat(dest, src);
#else
  int i = 0;
  while (dest[i] != '\0') i++;
  EG_strcpy(dest+i, src);
  return dest;
#endif
}


__HOST_AND_DEVICE__ int
EG_strncmp(const char *str_a, const char *str_b, size_t n)
{
#ifndef __CUDA_ARCH__
  return strncmp(str_a, str_b, n);
#else
  int match = '\0';
  size_t i = 0;
  size_t done = 0;
  while ((i < n) && (match == 0) && !done) {
    if ((str_a[i] == '\0') || (str_b[i] == '\0')) {
      done = 1;
    } else if (str_a[i] != str_b[i]) {
      match = i+1;
      if (((int)str_a[i] - (int)str_b[i]) < 0) match = 0 - (i + 1);
    }
    ++i;
  }
  return match;
#endif
}
/*@+mayaliasunique@*/
