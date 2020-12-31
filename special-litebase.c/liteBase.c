/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Base Object Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteClasses.h"
#include "emp.h"
#include "liteDevice.h"

#if !defined(WIN32) && !defined(__CYGWIN__)
#include <execinfo.h>
#endif

#define STRING(a)       #a
#define STR(a)          STRING(a)

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __PROTO_H_AND_D__
#undef __PROTO_H_AND_D__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#define __PROTO_H_AND_D__   extern "C" __host__ __device__
  extern "C" int EG_open(egObject **context);
  extern "C" __host__ __device__ int EG_close(egObject *context);
  extern "C" int EG_loadModel(egObject *context, /*@unused@*/ int bflg,
                              const char *name, egObject **model);
#else
#define __HOST_AND_DEVICE__
#define __PROTO_H_AND_D__ extern
#endif

__PROTO_H_AND_D__ int  EG_importModel( egObject *context, const size_t nbytes,
                                       const char *stream, egObject **model );
__PROTO_H_AND_D__ int  EG_exactInit( );


static const char *EGADSprop[2] = {STR(EGADSPROP),
                  "\nEGADSprop: Copyright 2011-2019 MIT. All Rights Reserved."};


#ifndef __CUDA_ARCH__
static void
EG_traceback()
{
#if !defined(WIN32) && !defined(__CYGWIN__)
  int    i;
  void   *array[100];
  size_t size;
  char   **strings;
  
  size = backtrace(array, 100);
  strings = backtrace_symbols (array, size);
  i = size;
  printf ("\nObtained %d stack frames:\n", i);
  for (i = 0; i < size; i++)
    printf ("%s\n", strings[i]);
  free (strings);
  printf("\n");
#endif
}
#endif


__HOST_AND_DEVICE__ void
EG_revision(int *major, int *minor, const char **OCCrev)
{
  *major  = EGADSMAJOR;
  *minor  = EGADSMINOR;
/*@-observertrans@*/
  *OCCrev = "(lite version)";
/*@+observertrans@*/
}


/*@-nullret@*/
__HOST_AND_DEVICE__ static int
EG_freeBlind(egObject *object)
{
  liteGeometry *lgeom;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;
  egObject     object_, *object_h = &object_;

  if  (object == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(object_h, object);
  if  (object_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object_h->oclass < PCURVE) ||
      (object_h->oclass > MODEL))      return EGADS_NOTTOPO;
  if (object_h->blind == NULL)         return EGADS_SUCCESS;
  
  if (object_h->oclass <= SURFACE) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    lgeom = (liteGeometry *) object_h->blind;
    EG_GET_GEOM(lgeom_h, lgeom);
    if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
    EG_FREE(lgeom_h->data);
  } else if (object_h->oclass == LOOP) {
    liteLoop lloop_, *lloop_h = &lloop_;
    lloop = (liteLoop *) object_h->blind;
    EG_GET_LOOP(lloop_h, lloop);
    EG_FREE(lloop_h->edges);
    EG_FREE(lloop_h->senses);
  } else if (object_h->oclass == FACE) {
    liteFace lface_, *lface_h = &lface_;
    lface = (liteFace *) object_h->blind;
    EG_GET_FACE(lface_h, lface);
    EG_FREE(lface_h->loops);
    EG_FREE(lface_h->senses);
  } else if (object_h->oclass == SHELL) {
    liteShell lshell_, *lshell_h = &lshell_;
    lshell = (liteShell *) object_h->blind;
    EG_GET_SHELL(lshell_h, lshell);
    EG_FREE(lshell_h->faces);
  } else if (object_h->oclass == BODY) {
    liteBody lbody_, *lbody_h = &lbody_;
    lbody = (liteBody *) object_h->blind;
    EG_GET_BODY(lbody_h, lbody);
    EG_FREE(lbody_h->pcurves.objs);
    EG_FREE(lbody_h->curves.objs);
    EG_FREE(lbody_h->surfaces.objs);
    EG_FREE(lbody_h->nodes.objs);
    EG_FREE(lbody_h->edges.objs);
    EG_FREE(lbody_h->loops.objs);
    EG_FREE(lbody_h->faces.objs);
    EG_FREE(lbody_h->shells.objs);
    EG_FREE(lbody_h->senses);
  } else if (object_h->oclass == MODEL) {
    liteModel lmodel_, *lmodel_h = &lmodel_;
    lmodel = (liteModel *) object_h->blind;
    EG_GET_MODEL(lmodel_h, lmodel);
    EG_FREE(lmodel_h->bodies);
  }
  
  EG_FREE(object_h->blind);
  object_h->blind = NULL;
  EG_SET_OBJECT(&object, object_h);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_makeObject(/*@null@*/ egObject *context, egObject **obj)
{
  int      outLevel;
  egObject *object, *prev;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egCntxt  *cntx;
  egObject context_, *context_h = &context_;
  egObject object_, *object_h = &object_;
  
  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context_h->blind;
  if (cntx == NULL)                    return EGADS_NODATA;
  EG_GET_CNTXT(cntx_h, cntx);
  outLevel = cntx_h->outLevel;
  if (cntx_h->mutex != NULL) EMP_LockSet(cntx_h->mutex);

  /* any objects in the pool? */
  object = cntx_h->pool;
  if (object == NULL) {
    EG_NEW_OBJECT(&object);
    if (object == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc on Object (EG_makeObject)!\n");
      if (cntx_h->mutex != NULL) EMP_LockRelease(cntx_h->mutex);
      return EGADS_MALLOC;
    }
  } else {
    EG_GET_OBJECT(object_h, object);
    EG_SET_OBJECT_PTR(&(cntx->pool), &(object_h->next));
    object_h->prev = NULL;
  }

  prev                  = cntx_h->last;
  object_h->magicnumber = MAGIC;
  object_h->oclass      = NIL;
  object_h->mtype       = 0;
  object_h->tref        = NULL;
  object_h->attrs       = NULL;
  object_h->blind       = NULL;
  object_h->topObj      = context;
  object_h->prev        = prev;
  object_h->next        = NULL;
  EG_SET_OBJECT(&object, object_h);
  EG_SET_OBJECT_PTR(&(prev->next), &object);

  *obj = object;
  EG_SET_OBJECT_PTR(&(cntx->last), obj);
  if (cntx_h->mutex != NULL) EMP_LockRelease(cntx_h->mutex);
  return EGADS_SUCCESS;
}


int
EG_open(egObject **context)
{
  int      i, status = EGADS_SUCCESS;
  egObject *object;
  egCntxt  *cntx;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egObject object_, *object_h = &object_;
  void     *mutex;
  
  EG_NEW_CNTXT(&cntx);
  if (cntx == NULL) return EGADS_MALLOC;
  EG_NEW_OBJECT(&object);
  if (object == NULL) {
    EG_FREE_CNTXT(cntx);
    return EGADS_MALLOC;
  }
  for (i = 0; i < MTESSPARAM; i++) cntx_h->tess[i] = 0.0;
  cntx_h->outLevel  = 1;
  cntx_h->signature = (char **) EGADSprop;
  cntx_h->usrPtr    = NULL;
  cntx_h->threadID  = EMP_ThreadID();
  cntx_h->mutex     = EMP_LockCreate();
  cntx_h->pool      = NULL;
  cntx_h->last      = object;
  if (cntx_h->mutex == NULL)
    printf(" EMP Error: mutex creation = NULL (EG_open)!\n");
  EG_SET_CNTXT(cntx, cntx_h);
  if (EGADS_SUCCESS != status) {
    printf(" Incomplete context creation (EG_open)!\n");
    EG_FREE_CNTXT(cntx);
    return status;
  }
  
  object_h->magicnumber = MAGIC;
  object_h->oclass      = CONTXT;
  object_h->mtype       = 1;                /* lite version */
  object_h->tref        = NULL;             /* not used here */
  object_h->attrs       = NULL;
  object_h->blind       = cntx;
  object_h->topObj      = NULL;             /* our single model */
  object_h->prev        = NULL;
  object_h->next        = NULL;
  EG_SET_OBJECT(&object, object_h);

  EG_exactInit();
  *context = object;
  mutex = cntx_h->mutex;
  EG_GET_CNTXT(cntx_h, cntx);
  if (cntx_h->mutex == NULL && mutex != NULL) {	  /* cntx is on device */
    EMP_LockRelease(mutex);
    EMP_LockDestroy(mutex);
  }
  return EGADS_SUCCESS;
}
/*@+nullret@*/


__HOST_AND_DEVICE__ int
EG_referenceObject(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_referenceTopObj(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ /*@kept@*/ /*@null@*/ egObject *
EG_context(const egObject *obj)
{
  int      cnt;
  egObject *object, *topObj;
  egObject obj_, *obj_h = &obj_;
  egObject object_, *object_h = &object_;
  
  if (obj == NULL) {
    printf(" EGADS Internal: EG_context called with NULL!\n");
    return NULL;
  }
  EG_GET_OBJECT(obj_h, obj);
  if (obj_h->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context Object NOT an ego!\n");
    return NULL;
  }
  if (obj_h->oclass == CONTXT) return (egObject *) obj;
  
  object = obj_h->topObj;
  if (object == NULL) {
    printf(" EGADS Internal: EG_context topObj is NULL!\n");
    return NULL;
  }
  EG_GET_OBJECT(object_h, object);
  if (object_h->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context topObj NOT an ego!\n");
    return NULL;
  }
  if (object_h->oclass == CONTXT) return object;
  
  cnt = 0;
  do {
    egObject topObj_, *topObj_h = &topObj_;
    topObj = object_h->topObj;
    if (topObj == NULL) {
      printf(" EGADS Internal: %d EG_context contents of topObj is NULL!\n",
             cnt);
      return NULL;
    }
    EG_GET_OBJECT(topObj_h, topObj);
    if (topObj_h->magicnumber != MAGIC) {
      printf(" EGADS Internal: %d EG_context contents of topObj NOT an ego!\n",
             cnt);
      return NULL;
    }
    if (topObj_h->oclass == CONTXT) return topObj;
    object = topObj;
    EG_GET_OBJECT(object_h, object);
    cnt++;
  } while (object != NULL);
  
  printf(" EGADS Internal: Cannot find context -- depth = %d!\n", cnt);
#ifndef __CUDA_ARCH__
  EG_traceback();
#endif
  return NULL;
}


__HOST_AND_DEVICE__ int
EG_sameThread(const egObject *obj)
{
  egObject *context;
#ifndef __CUDA_ARCH__
  egCntxt  *cntxt;
#endif
  
  if (obj == NULL)               return 1;
  if (obj->magicnumber != MAGIC) return 1;
  context = EG_context(obj);
  if (context == NULL)           return 1;
  
#ifndef __CUDA_ARCH__
  cntxt = (egCntxt *) context->blind;
  if (cntxt->threadID == EMP_ThreadID()) return 0;
#endif
  return 1;
}


__HOST_AND_DEVICE__ int
EG_outLevel(const egObject *obj)
{
  egObject *context;
  egCntxt  *cntxt;
  
  if (obj == NULL)               return 0;
  if (obj->magicnumber != MAGIC) return 0;
  context = EG_context(obj);
  if (context == NULL)           return 0;
  
  cntxt = (egCntxt *) context->blind;
  return cntxt->outLevel;
}


__HOST_AND_DEVICE__ int
EG_setOutLevel(egObject *context, int outLevel)
{
  int     old;
  egCntxt *cntx;
  
  if  (context == NULL)                 return EGADS_NULLOBJ;
  if  (context->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if  (context->oclass != CONTXT)       return EGADS_NOTCNTX;
  if ((outLevel < 0) || (outLevel > 3)) return EGADS_RANGERR;
  cntx = (egCntxt *) context->blind;
  if  (cntx == NULL)                    return EGADS_NODATA;
  old            = cntx->outLevel;
  cntx->outLevel = outLevel;
  
  return old;
}


__HOST_AND_DEVICE__ int
EG_setTessParam(egObject *context, int iParam, double value, double *oldValue)
{
  egCntxt *cntx;
  
  *oldValue = 0.0;
  if  (context == NULL)                      return EGADS_NULLOBJ;
  if  (context->magicnumber != MAGIC)        return EGADS_NOTOBJ;
  if  (context->oclass != CONTXT)            return EGADS_NOTCNTX;
  if ((iParam < 1) || (iParam > MTESSPARAM)) return EGADS_RANGERR;
  cntx = (egCntxt *) context->blind;
  if  (cntx == NULL)                          return EGADS_NODATA;
  *oldValue            = cntx->tess[iParam-1];
  cntx->tess[iParam-1] = value;
  
  return EGADS_SUCCESS;
}


int
EG_loadModel(egObject *context, /*@unused@*/ int bflg, const char *name,
             egObject **model)
{
  int      status;
  size_t   nbytes, ntest;
  char     *stream;
  FILE     *fp;
  egObject context_, *context_h = &context_;

  *model = NULL;
  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;

  if (name != NULL) {
    fp = fopen(name, "rb");
    if (fp == NULL) return EGADS_NOTFOUND;

    fseek(fp, 0, SEEK_END);
    nbytes = ftell(fp);
    rewind(fp);

    stream = (char *) EG_alloc(nbytes+1);
    if (stream == NULL) return EGADS_MALLOC;

    ntest = fread(stream, sizeof(char), nbytes, fp);
    if (ntest != nbytes) {
      printf(" EGADSlite Error: Stream expected to be %zd long but is %zd!\n",
             nbytes, ntest);
      EG_free(stream);
      fclose(fp);
      return EGADS_NOLOAD;
    }

    fclose(fp);

    status = EG_importModel(context, nbytes, stream, model);

    EG_free(stream);
    return status;
  }

  EG_GET_OBJECT_PTR(model, &(context->topObj));
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_deleteObject(egObject *object)
{
  int      i, j;
  egObject *pobj, *nobj, *context;
  egCntxt  *cntx;
  egTessel *tess;
  egObject object_, *object_h = &object_;
  void     *nil = NULL;

  if (object == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(object_h, object);
  if (object_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object_h->oclass == EMPTY)      return EGADS_EMPTY;
  if (object_h->oclass == REFERENCE)  return EGADS_REFERCE;

  if (object_h->oclass == TESSELLATION) {
    egObject context_, *context_h = &context_;
    egCntxt  cntx_, *cntx_h = &cntx_;

    context = EG_context(object);
    if (context == NULL)              return EGADS_NOTCNTX;
    EG_GET_OBJECT(context_h, context);
    cntx = (egCntxt *) context_h->blind;
    if (cntx == NULL)                 return EGADS_NODATA;
    EG_GET_CNTXT(cntx_h, cntx);
    if (cntx_h->mutex != NULL) EMP_LockSet(cntx->mutex);
    tess = (egTessel *) object_h->blind;
    if (tess != NULL) {
      egTessel tess_, *tess_h = &tess_;
      EG_GET_TESSEL(tess_h, tess);
      if (tess_h->xyzs != NULL) EG_FREE(tess_h->xyzs);
      if (tess_h->tess1d != NULL) {
        egTess1D tess1d_, *tess1d_h = &tess1d_;
        for (i = 0; i < tess_h->nEdge; i++) {
          EG_GET_TESS1D(tess1d_h, &(tess_h->tess1d[i]));
          if (tess1d_.faces[0].faces != NULL)
            EG_FREE(tess1d_.faces[0].faces);
          if (tess1d_.faces[1].faces != NULL)
            EG_FREE(tess1d_.faces[1].faces);
          if (tess1d_.faces[0].tric  != NULL)
            EG_FREE(tess1d_.faces[0].tric);
          if (tess1d_.faces[1].tric  != NULL)
            EG_FREE(tess1d_.faces[1].tric);
          if (tess1d_.xyz    != NULL)
            EG_FREE(tess1d_.xyz);
          if (tess1d_.t      != NULL)
            EG_FREE(tess1d_.t);
          if (tess1d_.global != NULL)
            EG_FREE(tess1d_.global);
        }
        EG_FREE(tess_h->tess1d);
      }
      if (tess_h->tess2d != NULL) {
        egTess2D tess2d_, *tess2d_h = &tess2d_;
        for (i = 0; i < 2*tess_h->nFace; i++) {
          EG_GET_TESS2D(tess2d_h, &(tess_h->tess2d[i]));
          if (tess2d_.mKnots != NULL)
            EG_deleteObject(tess2d_.mKnots);
          if (tess2d_.xyz    != NULL)
            EG_FREE(tess2d_.xyz);
          if (tess2d_.uv     != NULL)
            EG_FREE(tess2d_.uv);
          if (tess2d_.global != NULL)
            EG_FREE(tess2d_.global);
          if (tess2d_.ptype  != NULL)
            EG_FREE(tess2d_.ptype);
          if (tess2d_.pindex != NULL)
            EG_FREE(tess2d_.pindex);
          if (tess2d_.bary   != NULL)
            EG_FREE(tess2d_.bary);
          if (tess2d_.frame != NULL)
            EG_FREE(tess2d_.frame);
          if (tess2d_.frlps != NULL)
            EG_FREE(tess2d_.frlps);
          if (tess2d_.tris   != NULL)
            EG_FREE(tess2d_.tris);
          if (tess2d_.tric   != NULL)
            EG_FREE(tess2d_.tric);
          if (tess2d_.patch  != NULL) {
            egPatch  patch_, *patch_h = &patch_;
            for (j = 0; j < tess2d_.npatch; j++) {
              EG_GET_PATCH(patch_h, &(tess2d_.patch[j]));
              if (patch_.ipts != NULL)
                EG_FREE(patch_.ipts);
              if (patch_.bounds != NULL)
                EG_FREE(patch_.bounds);
            }
            EG_FREE(tess2d_.patch);
          }
        }
        EG_FREE(tess_h->tess2d);
      }
      if (tess_h->globals != NULL) EG_FREE(tess_h->globals);
      EG_FREE(tess);
      object_h->oclass = EMPTY;
      object_h->blind  = NULL;
/*@-nullret@*/
      EG_SET_OBJECT(&object, object_h);
/*@+nullret@*/
      
      /* patch up the lists & put the object in the pool */
      pobj = object_h->prev;          /* always have a previous -- context! */
      nobj = object_h->next;
      if (nobj == NULL) {
        if (object != cntx_h->last)
          printf(" EGADS Info: Context Last NOT Object Next w/ NULL!\n");
        EG_SET_OBJECT_PTR(&(cntx->last), &pobj);
      } else {
        EG_SET_OBJECT_PTR(&(nobj->prev), &pobj);
      }
      if (pobj == NULL) {
        printf(" EGADS Info: PrevObj is NULL (EG_destroyObject)!\n");
      } else {
        EG_SET_OBJECT_PTR(&(pobj->next), &nobj);
      }
      EG_SET_OBJECT_PTR(&(object->prev), &nil);
      EG_SET_OBJECT_PTR(&(object->next), &(cntx_h->pool));
      EG_SET_OBJECT_PTR(&(cntx->pool), &object);
      if (cntx_h->mutex != NULL) EMP_LockRelease(cntx_h->mutex);
    }
    return EGADS_SUCCESS;
  }

  /* report other deletes? */
  
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_getInfo(const egObject *object, int *oclass, int *mtype, egObject **top,
           egObject **prev, egObject **next)
{
  egObject object_, *object_h = &object_;

  if (object == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(object_h, object);
  if (object_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object_h->oclass == EMPTY)      return EGADS_EMPTY;

  *oclass = object_h->oclass;
  *mtype  = object_h->mtype;
  *top    = object_h->topObj;
  *prev   = object_h->prev;
  *next   = object_h->next;

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_close(egObject *context)
{
  int      i, status;
  egAttrs  *attrs;
  egObject *obj, *next, *last;
  egObject obj_, *obj_h = &obj_;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egCntxt  *cntx;
  egObject context_, *context_h = &context_;

  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context_h->blind;
  if (cntx == NULL)                    return EGADS_NODATA;
  EG_GET_CNTXT(cntx_h, cntx);

  /* delete tessellation objects */
  
  obj  = context_h->next;
  last = NULL;
  while (obj != NULL) {
    EG_GET_OBJECT_PTR(&next, &(obj->next));
    EG_GET_OBJECT(obj_h, obj);
    if (obj_h->oclass == TESSELLATION)
      if (EG_deleteObject(obj) == EGADS_SUCCESS) {
        obj = last;
        if (obj == NULL) {
          next = context_h->next;
        } else {
          next = obj_h->next;
        }
      }
    last = obj;
    obj  = next;
  }

  /* delete all objects */
  if (cntx_h->mutex != NULL) EMP_LockSet(cntx_h->mutex);
  
  obj = context_h->next;
  while (obj != NULL) {
    EG_GET_OBJECT(obj_h, obj);
    next = obj_h->next;
    if ((obj_h->oclass >= PCURVE) && (obj_h->oclass <= MODEL)) {
      status = EG_freeBlind(obj);
      if (status != EGADS_SUCCESS)
        printf(" EGADS Info: freeBlind = %d in Cleanup (EG_close)!\n", status);
    }
    attrs = (egAttrs *) obj_h->attrs;
    if (attrs != NULL) {
      egAttrs attrs_, *attrs_h = &attrs_;
      egAttr  attr_, *attr_h = &attr_;
      EG_GET_ATTRS(attrs_h, attrs);
      for (i = 0; i < attrs_h->nattrs; i++) {
        EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
        EG_FREE(attr_.name);
        if (attr_.type == ATTRINT) {
          if (attr_.length > 1) EG_FREE(attr_.vals.integers);
        } else if ((attr_.type == ATTRREAL) ||
                   (attr_.type == ATTRCSYS)) {
          if (attr_.length > 1) EG_FREE(attr_.vals.reals);
        } else if (attr_.type == ATTRSTRING) {
          EG_FREE(attr_.vals.string);
        }
      }
      EG_FREE(attrs_h->attrs);
      EG_FREE(attrs);
    }
    EG_FREE(obj);
    obj = next;
  }
  
  /* clean up the pool */
  
  obj = cntx_h->pool;
  while (obj != NULL) {
    EG_GET_OBJECT(obj_h, obj);
    if (obj_h->magicnumber != MAGIC) {
      printf(" EGADS Info: Found BAD Object in Cleanup (EG_close)!\n");
      printf("             Class = %d\n", obj_h->oclass);
      break;
    }
    next = obj_h->next;
    EG_FREE(obj);
    obj = next;
  }

  context_h->magicnumber = 0;
  context_h->oclass      = EMPTY;
  EG_FREE(context);
  if (cntx_h->mutex != NULL) EMP_LockRelease(cntx_h->mutex);
  if (cntx_h->mutex != NULL) EMP_LockDestroy(cntx_h->mutex);
  EG_FREE(cntx);
  
  return EGADS_SUCCESS;
}
