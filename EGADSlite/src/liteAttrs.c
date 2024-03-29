/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Attribute Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef __CUDACC__
#include "egadsString.h"
#else
#include <string.h>
#endif
#include "egadsTypes.h"
#include "egadsInternals_lite.h"

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



__HOST_AND_DEVICE__ int
EGlite_attributeNumSeq(const egObject *obj, const char *name, int *num)
{
  int     i, length, find = -1;
  egAttrs *attrs;

  *num = 0;
  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (name == NULL)              return EGADS_NONAME;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_SUCCESS;
  
  length = strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] == ' ') return EGADS_SEQUERR;
  
  for (i = 0; i < attrs->nseqs; i++)
    if (strcmp(attrs->seqs[i].root,name) == 0) {
      find = i;
      break;
    }
  if (find == -1) return EGADS_SUCCESS;

  *num = attrs->seqs[find].nSeq;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_attributeNum(const egObject *obj, int *num)
{
  egAttrs *attrs;

  *num = 0;
  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_SUCCESS;

  *num = attrs->nattrs;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_attributeGet(const egObject *obj, int index, const char **name, 
                int *atype, int *len, /*@null@*/ const int **ints, 
                /*@null@*/ const double **reals, /*@null@*/ const char **str)
{
  egAttrs *attrs;

  *name  = NULL;
  *atype = 0;
  if (ints  != NULL) *ints  = NULL;
  if (reals != NULL) *reals = NULL;
  if (str   != NULL) *str   = NULL;
  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) {
    printf(" EGADS Error: NULL Attributes (EGlite_attributeGet)!\n");
    return EGADS_INDEXERR;
  }
  if ((index < 1) || (index > attrs->nattrs)) {
    printf(" EGADS Error: Index Error %d [1-%d] (EGlite_attributeGet)!\n",
           index, attrs->nattrs);
    return EGADS_INDEXERR;
  }

  *name  = attrs->attrs[index-1].name;
  *atype = attrs->attrs[index-1].type;
  *len   = attrs->attrs[index-1].length;
  if (*atype == ATTRINT) {
    if (ints != NULL) 
      if (*len <= 1) {
        *ints = &attrs->attrs[index-1].vals.integer;
      } else {
        *ints =  attrs->attrs[index-1].vals.integers;
      }
  } else if (*atype == ATTRREAL) {
    if (reals != NULL)
      if (*len <= 1) {
        *reals = &attrs->attrs[index-1].vals.real;
      } else {
        *reals =  attrs->attrs[index-1].vals.reals;
      }
  } else if (*atype == ATTRCSYS) {
    *len = attrs->attrs[index-1].length - 12;
    if (reals != NULL) *reals = attrs->attrs[index-1].vals.reals;
    if (*len < 0) *len = 0;
  } else {
    if (str != NULL) *str = attrs->attrs[index-1].vals.string;
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_attributeRet(const egObject *obj, const char *name, int *atype, 
                int *len, /*@null@*/ const int **ints, 
                          /*@null@*/ const double **reals, 
                          /*@null@*/ const char **str)
{
  int     i, index;
  egAttrs *attrs;

  *atype = 0;
  *len   = 0;
  if (ints  != NULL) *ints  = NULL;
  if (reals != NULL) *reals = NULL;
  if (str   != NULL) *str   = NULL;
  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (name == NULL)              return EGADS_NONAME;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_NOTFOUND;

  index = -1;
  for (i = 0; i < attrs->nattrs; i++)
#ifdef __CUDACC__
    if (EGlite_strncmp(attrs->attrs[i].name, name, 256) == 0) {
#else
    if (strcmp(attrs->attrs[i].name, name) == 0) {
#endif
      index = i;
      break;
    }
  if (index == -1) return EGADS_NOTFOUND;

  *atype = attrs->attrs[index].type;
  *len   = attrs->attrs[index].length;
  if (*atype == ATTRINT) {
    if (ints != NULL) 
      if (*len <= 1) {
        *ints = &attrs->attrs[index].vals.integer;
      } else {
        *ints =  attrs->attrs[index].vals.integers;
      }
  } else if (*atype == ATTRREAL) {
    if (reals != NULL)
      if (*len <= 1) {
        *reals = &attrs->attrs[index].vals.real;
      } else {
        *reals =  attrs->attrs[index].vals.reals;
      }
  } else if (*atype == ATTRCSYS) {
    *len = attrs->attrs[index].length - 12;
    if (reals != NULL) *reals = attrs->attrs[index].vals.reals;
    if (*len < 0) *len = 0;
  } else {
    if (str != NULL) *str = attrs->attrs[index].vals.string;
  }
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_attributeRetSeq(const egObject *obj, const char *name, int index, int *atype,
                   int *len, /*@null@*/ const int **ints,
                             /*@null@*/ const double **reals,
                             /*@null@*/ const char **str)
{
  int     i, length, stat, find = -1;
  char    *fullname;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (name == NULL)              return EGADS_NONAME;
  if (index <= 0)                return EGADS_INDEXERR;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_NOTFOUND;
  
  length = strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] == ' ') return EGADS_SEQUERR;
  
  for (i = 0; i < attrs->nseqs; i++)
    if (strcmp(attrs->seqs[i].root,name) == 0) {
      find = i;
      break;
    }
  if (find == -1) {
    if (index != 1) return EGADS_INDEXERR;
    return EGlite_attributeRet(obj, name, atype, len, ints, reals, str);
  }
  if (index > attrs->seqs[find].nSeq) {
    printf(" EGADS Error: Index %d [1-%d] (EGlite_attributeRetSeq)!\n",
           index, attrs->seqs[find].nSeq);
    return EGADS_INDEXERR;
  }
  
  fullname = (char *) EGlite_alloc((length+8)*sizeof(char));
  if (fullname == NULL) return EGADS_MALLOC;
  snprintf(fullname, length+8, "%s %d", name, index);
  stat = EGlite_attributeRet(obj, fullname, atype, len, ints, reals, str);
  EGlite_free(fullname);
  
  return stat;
}
