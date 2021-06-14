/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Lite Base Object Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "egadsTypes_lite.h"
#include "egadsInternals_lite.h"
#include "liteClasses.h"
#include "emp_lite.h"
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
  extern "C" int EGlite_open(egObject **context);
  extern "C" __host__ __device__ int EGlite_close(egObject *context);
  extern "C" int EGlite_loadModel(egObject *context, /*@unused@*/ int bflg,
                              const char *name, egObject **model);
#else
#define __HOST_AND_DEVICE__
#define __PROTO_H_AND_D__ extern
#endif

__PROTO_H_AND_D__ int  EGlite_importModel( egObject *context, const size_t nbytes,
                                       const char *stream, egObject **model );
__PROTO_H_AND_D__ int  EGlite_exactInit( );


static const char *EGADSprop[2] = {STR(EGADSPROP),
                  "\nEGADSprop: Copyright 2011-2019 MIT. All Rights Reserved."};


#ifndef __CUDA_ARCH__
static void
EGlite_traceback()
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
EGlite_revision(int *major, int *minor, const char **OCCrev)
{
  *major  = EGADSMAJOR;
  *minor  = EGADSMINOR;
/*@-observertrans@*/
  *OCCrev = "(lite version)";
/*@+observertrans@*/
}


/*@-nullret@*/
__HOST_AND_DEVICE__ static int
EGlite_freeBlind(egObject *object)
{
  liteGeometry *lgeom;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;
  egObject     object_, *object_h = &object_;

  if  (object == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(object_h, object);
  if  (object_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object_h->oclass < PCURVE) ||
      (object_h->oclass > MODEL))      return EGADS_NOTTOPO;
  if (object_h->blind == NULL)         return EGADS_SUCCESS;
  
  if (object_h->oclass <= SURFACE) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    lgeom = (liteGeometry *) object_h->blind;
    EGlite_GET_GEOM(lgeom_h, lgeom);
    if (lgeom_h->header != NULL) EGlite_FREE(lgeom_h->header);
    EGlite_FREE(lgeom_h->data);
  } else if (object_h->oclass == LOOP) {
    liteLoop lloop_, *lloop_h = &lloop_;
    lloop = (liteLoop *) object_h->blind;
    EGlite_GET_LOOP(lloop_h, lloop);
    EGlite_FREE(lloop_h->edges);
    EGlite_FREE(lloop_h->senses);
  } else if (object_h->oclass == FACE) {
    liteFace lface_, *lface_h = &lface_;
    lface = (liteFace *) object_h->blind;
    EGlite_GET_FACE(lface_h, lface);
    EGlite_FREE(lface_h->loops);
    EGlite_FREE(lface_h->senses);
  } else if (object_h->oclass == SHELL) {
    liteShell lshell_, *lshell_h = &lshell_;
    lshell = (liteShell *) object_h->blind;
    EGlite_GET_SHELL(lshell_h, lshell);
    EGlite_FREE(lshell_h->faces);
  } else if (object_h->oclass == BODY) {
    liteBody lbody_, *lbody_h = &lbody_;
    lbody = (liteBody *) object_h->blind;
    EGlite_GET_BODY(lbody_h, lbody);
    EGlite_FREE(lbody_h->pcurves.objs);
    EGlite_FREE(lbody_h->curves.objs);
    EGlite_FREE(lbody_h->surfaces.objs);
    EGlite_FREE(lbody_h->nodes.objs);
    EGlite_FREE(lbody_h->edges.objs);
    EGlite_FREE(lbody_h->loops.objs);
    EGlite_FREE(lbody_h->faces.objs);
    EGlite_FREE(lbody_h->shells.objs);
    EGlite_FREE(lbody_h->senses);
  } else if (object_h->oclass == MODEL) {
    liteModel lmodel_, *lmodel_h = &lmodel_;
    lmodel = (liteModel *) object_h->blind;
    EGlite_GET_MODEL(lmodel_h, lmodel);
    EGlite_FREE(lmodel_h->bodies);
  }
  
  EGlite_FREE(object_h->blind);
  object_h->blind = NULL;
  EGlite_SET_OBJECT(&object, object_h);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_makeObject(/*@null@*/ egObject *context, egObject **obj)
{
  int      outLevel;
  egObject *object, *prev;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egCntxt  *cntx;
  egObject context_, *context_h = &context_;
  egObject object_, *object_h = &object_;
  
  if (context == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context_h->blind;
  if (cntx == NULL)                    return EGADS_NODATA;
  EGlite_GET_CNTXT(cntx_h, cntx);
  outLevel = cntx_h->outLevel;
  if (cntx_h->mutex != NULL) EMP_LITE_LockSet(cntx_h->mutex);

  /* any objects in the pool? */
  object = cntx_h->pool;
  if (object == NULL) {
    EGlite_NEW_OBJECT(&object);
    if (object == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc on Object (EGlite_makeObject)!\n");
      if (cntx_h->mutex != NULL) EMP_LITE_LockRelease(cntx_h->mutex);
      return EGADS_MALLOC;
    }
  } else {
    EGlite_GET_OBJECT(object_h, object);
    EGlite_SET_OBJECT_PTR(&(cntx->pool), &(object_h->next));
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
  EGlite_SET_OBJECT(&object, object_h);
  EGlite_SET_OBJECT_PTR(&(prev->next), &object);

  *obj = object;
  EGlite_SET_OBJECT_PTR(&(cntx->last), obj);
  if (cntx_h->mutex != NULL) EMP_LITE_LockRelease(cntx_h->mutex);
  return EGADS_SUCCESS;
}


int
EGlite_open(egObject **context)
{
  int      i, status = EGADS_SUCCESS;
  egObject *object;
  egCntxt  *cntx;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egObject object_, *object_h = &object_;
  void     *mutex;
  
  EGlite_NEW_CNTXT(&cntx);
  if (cntx == NULL) return EGADS_MALLOC;
  EGlite_NEW_OBJECT(&object);
  if (object == NULL) {
    EGlite_FREE_CNTXT(cntx);
    return EGADS_MALLOC;
  }
  for (i = 0; i < MTESSPARAM; i++) cntx_h->tess[i] = 0.0;
  cntx_h->outLevel  = 1;
  cntx_h->signature = (char **) EGADSprop;
  cntx_h->usrPtr    = NULL;
  cntx_h->threadID  = EMP_LITE_ThreadID();
  cntx_h->mutex     = EMP_LITE_LockCreate();
  cntx_h->pool      = NULL;
  cntx_h->last      = object;
  if (cntx_h->mutex == NULL)
    printf(" EMP Error: mutex creation = NULL (EGlite_open)!\n");
  EGlite_SET_CNTXT(cntx, cntx_h);
  if (EGADS_SUCCESS != status) {
    printf(" Incomplete context creation (EGlite_open)!\n");
    EGlite_FREE_CNTXT(cntx);
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
  EGlite_SET_OBJECT(&object, object_h);

  EGlite_exactInit();
  *context = object;
  mutex = cntx_h->mutex;
  EGlite_GET_CNTXT(cntx_h, cntx);
  if (cntx_h->mutex == NULL && mutex != NULL) {	  /* cntx is on device */
    EMP_LITE_LockRelease(mutex);
    EMP_LITE_LockDestroy(mutex);
  }
  return EGADS_SUCCESS;
}
/*@+nullret@*/


__HOST_AND_DEVICE__ int
EGlite_referenceObject(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_referenceTopObj(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ /*@kept@*/ /*@null@*/ egObject *
EGlite_context(const egObject *obj)
{
  int      cnt;
  egObject *object, *topObj;
  egObject obj_, *obj_h = &obj_;
  egObject object_, *object_h = &object_;
  
  if (obj == NULL) {
    printf(" EGADS Internal: EGlite_context called with NULL!\n");
    return NULL;
  }
  EGlite_GET_OBJECT(obj_h, obj);
  if (obj_h->magicnumber != MAGIC) {
    printf(" EGADS Internal: EGlite_context Object NOT an ego!\n");
    return NULL;
  }
  if (obj_h->oclass == CONTXT) return (egObject *) obj;
  
  object = obj_h->topObj;
  if (object == NULL) {
    printf(" EGADS Internal: EGlite_context topObj is NULL!\n");
    return NULL;
  }
  EGlite_GET_OBJECT(object_h, object);
  if (object_h->magicnumber != MAGIC) {
    printf(" EGADS Internal: EGlite_context topObj NOT an ego!\n");
    return NULL;
  }
  if (object_h->oclass == CONTXT) return object;
  
  cnt = 0;
  do {
    egObject topObj_, *topObj_h = &topObj_;
    topObj = object_h->topObj;
    if (topObj == NULL) {
      printf(" EGADS Internal: %d EGlite_context contents of topObj is NULL!\n",
             cnt);
      return NULL;
    }
    EGlite_GET_OBJECT(topObj_h, topObj);
    if (topObj_h->magicnumber != MAGIC) {
      printf(" EGADS Internal: %d EGlite_context contents of topObj NOT an ego!\n",
             cnt);
      return NULL;
    }
    if (topObj_h->oclass == CONTXT) return topObj;
    object = topObj;
    EGlite_GET_OBJECT(object_h, object);
    cnt++;
  } while (object != NULL);
  
  printf(" EGADS Internal: Cannot find context -- depth = %d!\n", cnt);
#ifndef __CUDA_ARCH__
  EGlite_traceback();
#endif
  return NULL;
}


__HOST_AND_DEVICE__ int
EGlite_sameThread(const egObject *obj)
{
  egObject *context;
#ifndef __CUDA_ARCH__
  egCntxt  *cntxt;
#endif
  
  if (obj == NULL)               return 1;
  if (obj->magicnumber != MAGIC) return 1;
  context = EGlite_context(obj);
  if (context == NULL)           return 1;
  
#ifndef __CUDA_ARCH__
  cntxt = (egCntxt *) context->blind;
  if (cntxt->threadID == EMP_LITE_ThreadID()) return 0;
#endif
  return 1;
}


__HOST_AND_DEVICE__ int
EGlite_outLevel(const egObject *obj)
{
  egObject *context;
  egCntxt  *cntxt;
  
  if (obj == NULL)               return 0;
  if (obj->magicnumber != MAGIC) return 0;
  context = EGlite_context(obj);
  if (context == NULL)           return 0;
  
  cntxt = (egCntxt *) context->blind;
  return cntxt->outLevel;
}


__HOST_AND_DEVICE__ int
EGlite_setOutLevel(egObject *context, int outLevel)
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
EGlite_setTessParam(egObject *context, int iParam, double value, double *oldValue)
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
EGlite_loadModel(egObject *context, /*@unused@*/ int bflg, const char *name,
             egObject **model)
{
  int      status;
  size_t   nbytes, ntest;
  char     *stream;
  FILE     *fp;
  egObject context_, *context_h = &context_;

  *model = NULL;
  if (context == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;

  if (name != NULL) {
    fp = fopen(name, "rb");
    if (fp == NULL) return EGADS_NOTFOUND;

    fseek(fp, 0, SEEK_END);
    nbytes = ftell(fp);
    rewind(fp);

    stream = (char *) EGlite_alloc(nbytes+1);
    if (stream == NULL) return EGADS_MALLOC;

    ntest = fread(stream, sizeof(char), nbytes, fp);
    if (ntest != nbytes) {
      printf(" EGADSlite Error: Stream expected to be %zd long but is %zd!\n",
             nbytes, ntest);
      EGlite_free(stream);
      fclose(fp);
      return EGADS_NOLOAD;
    }

    fclose(fp);

    status = EGlite_importModel(context, nbytes, stream, model);

    EGlite_free(stream);
    return status;
  }

  EGlite_GET_OBJECT_PTR(model, &(context->topObj));
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_deleteObject(egObject *object)
{
  int      i, j;
  egObject *pobj, *nobj, *context;
  egCntxt  *cntx;
  egTessel *tess;
  egObject object_, *object_h = &object_;
  void     *nil = NULL;

  if (object == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(object_h, object);
  if (object_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object_h->oclass == EMPTY)      return EGADS_EMPTY;
  if (object_h->oclass == REFERENCE)  return EGADS_REFERCE;

  if (object_h->oclass == TESSELLATION) {
    egObject context_, *context_h = &context_;
    egCntxt  cntx_, *cntx_h = &cntx_;

    context = EGlite_context(object);
    if (context == NULL)              return EGADS_NOTCNTX;
    EGlite_GET_OBJECT(context_h, context);
    cntx = (egCntxt *) context_h->blind;
    if (cntx == NULL)                 return EGADS_NODATA;
    EGlite_GET_CNTXT(cntx_h, cntx);
    if (cntx_h->mutex != NULL) EMP_LITE_LockSet(cntx->mutex);
    tess = (egTessel *) object_h->blind;
    if (tess != NULL) {
      egTessel tess_, *tess_h = &tess_;
      EGlite_GET_TESSEL(tess_h, tess);
      if (tess_h->xyzs != NULL) EGlite_FREE(tess_h->xyzs);
      if (tess_h->tess1d != NULL) {
        egTess1D tess1d_, *tess1d_h = &tess1d_;
        for (i = 0; i < tess_h->nEdge; i++) {
          EGlite_GET_TESS1D(tess1d_h, &(tess_h->tess1d[i]));
          if (tess1d_.faces[0].faces != NULL)
            EGlite_FREE(tess1d_.faces[0].faces);
          if (tess1d_.faces[1].faces != NULL)
            EGlite_FREE(tess1d_.faces[1].faces);
          if (tess1d_.faces[0].tric  != NULL)
            EGlite_FREE(tess1d_.faces[0].tric);
          if (tess1d_.faces[1].tric  != NULL)
            EGlite_FREE(tess1d_.faces[1].tric);
          if (tess1d_.xyz    != NULL)
            EGlite_FREE(tess1d_.xyz);
          if (tess1d_.t      != NULL)
            EGlite_FREE(tess1d_.t);
          if (tess1d_.global != NULL)
            EGlite_FREE(tess1d_.global);
        }
        EGlite_FREE(tess_h->tess1d);
      }
      if (tess_h->tess2d != NULL) {
        egTess2D tess2d_, *tess2d_h = &tess2d_;
        for (i = 0; i < 2*tess_h->nFace; i++) {
          EGlite_GET_TESS2D(tess2d_h, &(tess_h->tess2d[i]));
          if (tess2d_.mKnots != NULL)
            EGlite_deleteObject(tess2d_.mKnots);
          if (tess2d_.xyz    != NULL)
            EGlite_FREE(tess2d_.xyz);
          if (tess2d_.uv     != NULL)
            EGlite_FREE(tess2d_.uv);
          if (tess2d_.global != NULL)
            EGlite_FREE(tess2d_.global);
          if (tess2d_.ptype  != NULL)
            EGlite_FREE(tess2d_.ptype);
          if (tess2d_.pindex != NULL)
            EGlite_FREE(tess2d_.pindex);
          if (tess2d_.bary   != NULL)
            EGlite_FREE(tess2d_.bary);
          if (tess2d_.frame != NULL)
            EGlite_FREE(tess2d_.frame);
          if (tess2d_.frlps != NULL)
            EGlite_FREE(tess2d_.frlps);
          if (tess2d_.tris   != NULL)
            EGlite_FREE(tess2d_.tris);
          if (tess2d_.tric   != NULL)
            EGlite_FREE(tess2d_.tric);
          if (tess2d_.patch  != NULL) {
            egPatch  patch_, *patch_h = &patch_;
            for (j = 0; j < tess2d_.npatch; j++) {
              EGlite_GET_PATCH(patch_h, &(tess2d_.patch[j]));
              if (patch_.ipts != NULL)
                EGlite_FREE(patch_.ipts);
              if (patch_.bounds != NULL)
                EGlite_FREE(patch_.bounds);
            }
            EGlite_FREE(tess2d_.patch);
          }
        }
        EGlite_FREE(tess_h->tess2d);
      }
      if (tess_h->globals != NULL) EGlite_FREE(tess_h->globals);
      EGlite_FREE(tess);
      object_h->oclass = EMPTY;
      object_h->blind  = NULL;
/*@-nullret@*/
      EGlite_SET_OBJECT(&object, object_h);
/*@+nullret@*/
      
      /* patch up the lists & put the object in the pool */
      pobj = object_h->prev;          /* always have a previous -- context! */
      nobj = object_h->next;
      if (nobj == NULL) {
        if (object != cntx_h->last)
          printf(" EGADS Info: Context Last NOT Object Next w/ NULL!\n");
        EGlite_SET_OBJECT_PTR(&(cntx->last), &pobj);
      } else {
        EGlite_SET_OBJECT_PTR(&(nobj->prev), &pobj);
      }
      if (pobj == NULL) {
        printf(" EGADS Info: PrevObj is NULL (EGlite_destroyObject)!\n");
      } else {
        EGlite_SET_OBJECT_PTR(&(pobj->next), &nobj);
      }
      EGlite_SET_OBJECT_PTR(&(object->prev), &nil);
      EGlite_SET_OBJECT_PTR(&(object->next), &(cntx_h->pool));
      EGlite_SET_OBJECT_PTR(&(cntx->pool), &object);
      if (cntx_h->mutex != NULL) EMP_LITE_LockRelease(cntx_h->mutex);
    }
    return EGADS_SUCCESS;
  }

  /* report other deletes? */
  
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EGlite_getInfo(const egObject *object, int *oclass, int *mtype, egObject **top,
           egObject **prev, egObject **next)
{
  egObject object_, *object_h = &object_;

  if (object == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(object_h, object);
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
EGlite_close(egObject *context)
{
  int      i, status;
  egAttrs  *attrs;
  egObject *obj, *next, *last;
  egObject obj_, *obj_h = &obj_;
  egCntxt  cntx_, *cntx_h = &cntx_;
  egCntxt  *cntx;
  egObject context_, *context_h = &context_;

  if (context == NULL)                 return EGADS_NULLOBJ;
  EGlite_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context_h->blind;
  if (cntx == NULL)                    return EGADS_NODATA;
  EGlite_GET_CNTXT(cntx_h, cntx);

  /* delete tessellation objects */
  
  obj  = context_h->next;
  last = NULL;
  while (obj != NULL) {
    EGlite_GET_OBJECT_PTR(&next, &(obj->next));
    EGlite_GET_OBJECT(obj_h, obj);
    if (obj_h->oclass == TESSELLATION)
      if (EGlite_deleteObject(obj) == EGADS_SUCCESS) {
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
  if (cntx_h->mutex != NULL) EMP_LITE_LockSet(cntx_h->mutex);
  
  obj = context_h->next;
  while (obj != NULL) {
    EGlite_GET_OBJECT(obj_h, obj);
    next = obj_h->next;
    if ((obj_h->oclass >= PCURVE) && (obj_h->oclass <= MODEL)) {
      status = EGlite_freeBlind(obj);
      if (status != EGADS_SUCCESS)
        printf(" EGADS Info: freeBlind = %d in Cleanup (EGlite_close)!\n", status);
    }
    attrs = (egAttrs *) obj_h->attrs;
    if (attrs != NULL) {
      egAttrs attrs_, *attrs_h = &attrs_;
      egAttr  attr_, *attr_h = &attr_;
      EGlite_GET_ATTRS(attrs_h, attrs);
      for (i = 0; i < attrs_h->nattrs; i++) {
        EGlite_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
        EGlite_FREE(attr_.name);
        if (attr_.type == ATTRINT) {
          if (attr_.length > 1) EGlite_FREE(attr_.vals.integers);
        } else if ((attr_.type == ATTRREAL) ||
                   (attr_.type == ATTRCSYS)) {
          if (attr_.length > 1) EGlite_FREE(attr_.vals.reals);
        } else if (attr_.type == ATTRSTRING) {
          EGlite_FREE(attr_.vals.string);
        }
      }
      EGlite_FREE(attrs_h->attrs);
      EGlite_FREE(attrs);
    }
    EGlite_FREE(obj);
    obj = next;
  }
  
  /* clean up the pool */
  
  obj = cntx_h->pool;
  while (obj != NULL) {
    EGlite_GET_OBJECT(obj_h, obj);
    if (obj_h->magicnumber != MAGIC) {
      printf(" EGADS Info: Found BAD Object in Cleanup (EGlite_close)!\n");
      printf("             Class = %d\n", obj_h->oclass);
      break;
    }
    next = obj_h->next;
    EGlite_FREE(obj);
    obj = next;
  }

  context_h->magicnumber = 0;
  context_h->oclass      = EMPTY;
  EGlite_FREE(context);
  if (cntx_h->mutex != NULL) EMP_LITE_LockRelease(cntx_h->mutex);
  if (cntx_h->mutex != NULL) EMP_LITE_LockDestroy(cntx_h->mutex);
  EGlite_FREE(cntx);
  
  return EGADS_SUCCESS;
}
