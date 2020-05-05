#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteClasses.h"
#include "liteDevice.h"


#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif


#ifdef __cplusplus
extern "C" {
#endif

__HOST_AND_DEVICE__ extern int EG_evaluate(const egObject *geom,
                                           const double *param, double *ev);
__HOST_AND_DEVICE__ extern int EG_attributeRet(const egObject *obj,
                                               const char *name, int *atype,
                                               int *len,
                                               /*@null@*/ const int **ints,
                                               /*@null@*/ const double **reals,
                                               /*@null@*/ const char **str);
__HOST_AND_DEVICE__ extern int EG_addStrAttr(egObject *obj, const char *name,
                                             const char *str);

int EG_attributeRetDev(const egObject *obj_d, const char *name, int *atype,
                       int *len, int **ints, double **reals, char **str);
int EG_addStrAttrDev(egObject *obj_d, const char *name, const char *str);
int EG_evaluateDev(const egObject *geom_d, const double *param, double *ev);

typedef struct {
  char         *name;
  int          atype;
  int          len;
  const int    *ints;
  const double *reals;
  char         *str;
  int          status;
} attribute_t;

typedef struct {
  double param[2];
  double ev[18];
  int    nev;
  int    status;
} evaluate_t;

#ifdef __cplusplus
}
#endif


__device__
size_t EG_strlen(const char *s)
{
  size_t i;
  for (i = 0; ; ++i)
    if (s[i] == '\0') break;
  return i-1;
}


__device__
size_t EG_strcmp(const char *s1, const char *s2)
{
  size_t i = 0;
  do {
    if (s1[i] != s2[i]) {
      if (s1[i] < s2[i]) return (size_t)-1;
      return 1;
    }
  } while (s1[i++] != '\0');
  return 0;
}


__device__
static int EG_addStrAttr_d(egObject *obj, const char *name, const char *str)
{
  int     i, length, find = -1;
  egAttr  *attr;
  egAttrs *attrs;
  
  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  
  if ((name == NULL) || (str == NULL)) {
      printf(" EGADS Internal: NULL Name/Value (EG_addStrAttr)!\n");
    return EGADS_NONAME;
  }
  length = EG_strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] <= ' ') {
      length = 0;
      break;
    }
  if (length == 0) {
    printf(" EGADS Internal: BAD Name (EG_addStrAttr)!\n");
    return EGADS_INDEXERR;
  }
  attrs = (egAttrs *) obj->attrs;
  
  if (attrs != NULL)
    for (i = 0; i < attrs->nattrs; i++)
      if (EG_strcmp(attrs->attrs[i].name,name) == 0) {
        find = i;
        break;
      }
  
  if ((find != -1) && (attrs != NULL)) {
    
    /* an existing attribute -- reset the values */
    
    if (attrs->attrs[find].type == ATTRINT) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.integers);
    } else if (attrs->attrs[find].type == ATTRREAL) {
      if (attrs->attrs[find].length != 1)
        EG_free(attrs->attrs[find].vals.reals);
    } else if (attrs->attrs[find].type == ATTRCSYS) {
      EG_free(attrs->attrs[find].vals.reals);
    } else if (attrs->attrs[find].type == ATTRSTRING) {
      EG_free(attrs->attrs[find].vals.string);
    }
    
  } else {
    
    if (attrs == NULL) {
      attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
      if (attrs == NULL) {
        printf(" EGADS Internal: Attrs MALLOC for %s (EG_addStrAttr)!\n",
               name);
        return EGADS_MALLOC;
      }
      attrs->nattrs = 0;
      attrs->attrs  = NULL;
      attrs->nseqs  = 0;
      attrs->seqs   = NULL;
      obj->attrs    = attrs;
    }
    if (attrs->attrs == NULL) {
      attr = (egAttr *) EG_alloc((attrs->nattrs+1)*sizeof(egAttr));
    } else {
      attr = (egAttr *) EG_reall(attrs->attrs,
                                 (attrs->nattrs+1)*sizeof(egAttr));
    }
    if (attr == NULL) {
      printf(" EGADS Internal: Attr MALLOC for %s (EG_addStrAttr)!\n",
             name);
      return EGADS_MALLOC;
    }
    attrs->attrs = attr;
    find = attrs->nattrs;
    attrs->attrs[find].vals.string = NULL;
    attrs->attrs[find].name        = EG_strdup(name);
    if (attrs->attrs[find].name == NULL) return EGADS_MALLOC;
    attrs->nattrs += 1;
  }
  
  attrs->attrs[find].type        = ATTRSTRING;
  attrs->attrs[find].length      = 0;
  attrs->attrs[find].vals.string = EG_strdup(str);
  if (attrs->attrs[find].vals.string != NULL)
    attrs->attrs[find].length = EG_strlen(attrs->attrs[find].vals.string);

  return EGADS_SUCCESS;
}

__global__
void attributeRetKernel(const egObject *obj_d, attribute_t* attr_d)
{
  attr_d->status = EG_attributeRet(obj_d, attr_d->name, &attr_d->atype,
                                   &attr_d->len, &attr_d->ints,
                                   &attr_d->reals, (const char **)&attr_d->str);
}


__global__
void addStrAttrKernel(egObject *obj_d, attribute_t* attr_d)
{
  attr_d->status = EG_addStrAttr_d(obj_d, attr_d->name, attr_d->str);
}


__global__
void evaluateKernel(const egObject *geom_d, evaluate_t* eval_d)
{
  eval_d->status = EG_evaluate(geom_d, eval_d->param, eval_d->ev);

  switch(geom_d->oclass) {
    case(PCURVE):
      eval_d->nev = 6;
      break;
    case(CURVE):
    case(EDGE):
      eval_d->nev = 9;
      break;
    case(SURFACE):
    case(FACE):
      eval_d->nev = 18;
      break;
    default:
      eval_d->nev = 3;
  }
}


int EG_attributeRetDev(const egObject *obj_d, const char *name, int *atype,
                       int *len, int **ints, double **reals, char **str)
{
  attribute_t *attr_d;
  attribute_t attr_, *attr_h = &attr_;
  int         status;

  DEVICE_MALLOC((void**)&attr_d, sizeof(attribute_t));
  EG_SET_STR(&(attr_d->name), name);

  attributeRetKernel<<<1, 1>>>(obj_d, attr_d);
  __CUDA_ERROR_CHECK;

  MEMCPY_DEVICE_TO_HOST(attr_h, attr_d, sizeof(attribute_t));

  *atype = attr_h->atype;
  *len = attr_h->len;
  if (*atype == ATTRINT) {
    if (attr_h->ints != NULL)
      if (attr_h->len <= 1) {
        *ints = (int *) EG_alloc(attr_h->len*sizeof(int));
        MEMCPY_DEVICE_TO_HOST(*ints, attr_d->ints, attr_h->len*sizeof(int));
      } else {
        MEMCPY_DEVICE_TO_HOST(*ints, attr_d->ints, sizeof(int));
      }
  } else if (*atype == ATTRREAL) {
    if (reals != NULL)
      if (attr_h->len <= 1) {
        *reals = (double *) EG_alloc(attr_h->len*sizeof(double));
        MEMCPY_DEVICE_TO_HOST(*reals, attr_d->reals, attr_h->len*sizeof(double));
      } else {
        MEMCPY_DEVICE_TO_HOST(*reals, attr_d->reals, sizeof(double));
      }
  } else if (*atype == ATTRCSYS) {
    *len = attr_h->len - 12;
    if (reals != NULL) {
      *reals = (double *) EG_alloc(attr_h->len*sizeof(double));
      MEMCPY_DEVICE_TO_HOST(*reals, attr_d->reals, attr_h->len*sizeof(double));
    }
    if (*len < 0) *len = 0;
  } else if (*atype == ATTRSTRING) {
    if (str != NULL) {
      *str = (char *) EG_alloc((attr_h->len+1)*sizeof(char));
      *str[0] = '\0';
      if( 0 < attr_h->len )
        MEMCPY_DEVICE_TO_HOST(*str, attr_h->str, (attr_h->len+1)*sizeof(char));
    }
  }
  status = attr_h->status;

  EG_FREE(attr_h->name);
  EG_FREE(attr_d);

  return status;
}


int EG_addStrAttrDev(egObject *obj_d, const char *name, const char *str)
{
  attribute_t *attr_d;
  attribute_t attr_, *attr_h = &attr_;
  int         status;

  DEVICE_MALLOC((void**)&attr_d, sizeof(attribute_t));
  EG_SET_STR(&(attr_d->name), name);
  EG_SET_STR(&(attr_d->str), str);

  addStrAttrKernel<<<1, 1>>>(obj_d, attr_d);
  __CUDA_ERROR_CHECK;

  MEMCPY_DEVICE_TO_HOST(attr_h, attr_d, sizeof(attribute_t));
  status = attr_h->status;

  EG_FREE(attr_h->str);
  EG_FREE(attr_h->name);
  EG_FREE(attr_d);

  return status;
}


int EG_evaluateDev(const egObject *geom_d, const double *param, double *ev)
{
  evaluate_t *eval_d;
  evaluate_t eval_, *eval_h = &eval_;
  short      oclass;
  int        nparam;
  int        i;
  int        status;

  EG_GET(&oclass, &(geom_d->oclass), sizeof(short), 1);
  switch(oclass) {
    case(PCURVE):
    case(CURVE):
    case(EDGE):
      nparam = 1;
      break;
    default:
      nparam = 2;
  }

  DEVICE_MALLOC((void**)&eval_d, sizeof(evaluate_t));
  EG_COPY(&(eval_d->param), param, sizeof(double), nparam);

  evaluateKernel<<<1, 1>>>(geom_d, eval_d);
  __CUDA_ERROR_CHECK;

  MEMCPY_DEVICE_TO_HOST(eval_h, eval_d, sizeof(evaluate_t));
  for( i=0; i<eval_h->nev; ++i) ev[i] = eval_h->ev[i];
  status = eval_h->status;

  EG_FREE(eval_d);

  return status;
}
