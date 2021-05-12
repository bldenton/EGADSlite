#ifndef __UVMAP_LIB_H__

#define __UVMAP_LIB_H__

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <errno.h>
#include <sys/stat.h>

#ifdef _WIN32
  #include <process.h>
  #include <direct.h>
  #include <io.h> 
  #define snprintf _snprintf
  typedef __int64 LONG_LONG_int;
#else
  #include <libgen.h>
  #include <unistd.h>
  #include <sys/utsname.h>
  typedef int64_t LONG_LONG_int;
#endif

#ifndef __UG_TYPEDEF_H__
  #ifdef _UG_LONG_LONG_INT_
    #define _UG_INT64_	1
    typedef LONG_LONG_int INT_;
  #else
    typedef int INT_;
  #endif

  typedef INT_ INT_2D[2];
  typedef INT_ INT_3D[3];
  typedef double DOUBLE_2D[2];
  typedef double DOUBLE_3D[3];
#endif

#ifndef MAX
  #define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
  #define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef NINT
  #define NINT(x) ((INT_) (((x) >= 0.0) ? floor((x)+0.5) : -floor(0.5-(x))))
#endif

#ifndef EGADS_EXTRAPOL
  #define EGADS_EXTRAPOL -37
#endif
#ifndef EGADS_UVMAP
  #define EGADS_UVMAP -35
#endif
#ifndef EGADS_READERR
  #define EGADS_READERR -32
#endif
#ifndef EGADS_WRITERR
  #define EGADS_WRITERR -19
#endif
#ifndef EGADS_MALLOC
  #define EGADS_MALLOC -4
#endif
#ifndef EGADS_NOTFOUND
  #define EGADS_NOTFOUND -1
#endif
#ifndef EGADS_SUCCESS
  #define EGADS_SUCCESS 0
#endif

#include "UVMAP_LIB_INC.h"

#endif
