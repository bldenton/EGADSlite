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

#if !defined(WIN32) && !defined(__CYGWIN__)
#include <execinfo.h>
#endif

#define STRING(a)       #a
#define STR(a)          STRING(a)


extern int  EG_importModel(egObject *context, const size_t nbytes,
                           const char *stream, egObject **model);
extern void EG_exactInit( );

static char *EGADSprop[2] = {STR(EGADSPROP),
                             "\nEGADSprop: Copyright 2011-2020 MIT. All Rights Reserved."};


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


static int
EG_freeBlind(egObject *object)
{
  liteGeometry *lgeom;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;

  if  (object == NULL)               return EGADS_NULLOBJ;
  if  (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((object->oclass < PCURVE) ||
      (object->oclass > MODEL))      return EGADS_NOTTOPO;
  if (object->blind == NULL)         return EGADS_SUCCESS;
  
  if (object->oclass <= SURFACE) {
    lgeom = (liteGeometry *) object->blind;
    if (lgeom->header != NULL) EG_free(lgeom->header);
    EG_free(lgeom->data);
  } else if (object->oclass == LOOP) {
    lloop = (liteLoop *) object->blind;
    EG_free(lloop->edges);
    EG_free(lloop->senses);
  } else if (object->oclass == FACE) {
    lface = (liteFace *) object->blind;
    EG_free(lface->loops);
    EG_free(lface->senses);
  } else if (object->oclass == SHELL) {
    lshell = (liteShell *) object->blind;
    EG_free(lshell->faces);
  } else if (object->oclass == BODY) {
    lbody = (liteBody *) object->blind;
    EG_free(lbody->pcurves.objs);
    EG_free(lbody->curves.objs);
    EG_free(lbody->surfaces.objs);
    EG_free(lbody->nodes.objs);
    EG_free(lbody->edges.objs);
    EG_free(lbody->loops.objs);
    EG_free(lbody->faces.objs);
    EG_free(lbody->shells.objs);
    EG_free(lbody->senses);
  } else if (object->oclass == MODEL) {
    lmodel = (liteModel *) object->blind;
    EG_free(lmodel->bodies);
  }
  
  EG_free(object->blind);
  object->blind = NULL;
  return EGADS_SUCCESS;
}


int
EG_makeObject(/*@null@*/ egObject *context, egObject **obj)
{
  int      outLevel;
  egObject *object, *prev;
  egCntxt  *cntx;
  
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;
  outLevel = cntx->outLevel;
  if (cntx->mutex != NULL) EMP_LockSet(cntx->mutex);
  
  /* any objects in the pool? */
  object = cntx->pool;
  if (object == NULL) {
    object = (egObject *) EG_alloc(sizeof(egObject));
    if (object == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Malloc on Object (EG_makeObject)!\n");
      if (cntx->mutex != NULL) EMP_LockRelease(cntx->mutex);
      return EGADS_MALLOC;
    }
  } else {
    cntx->pool   = object->next;
    object->prev = NULL;
  }
  
  prev                = cntx->last;
  object->magicnumber = MAGIC;
  object->oclass      = NIL;
  object->mtype       = 0;
  object->tref        = NULL;
  object->attrs       = NULL;
  object->blind       = NULL;
  object->topObj      = context;
  object->prev        = prev;
  object->next        = NULL;
  prev->next          = object;
  
  *obj = object;
  cntx->last = *obj;
  if (cntx->mutex != NULL) EMP_LockRelease(cntx->mutex);
  return EGADS_SUCCESS;
}


void
EG_revision(int *major, int *minor, char **OCCrev)
{
  *major  = EGADSMAJOR;
  *minor  = EGADSMINOR;
  *OCCrev = NULL;
}


int
EG_open(egObject **context)
{
  int      i;
  egObject *object;
  egCntxt  *cntx;
  
  cntx   = (egCntxt *) EG_alloc(sizeof(egCntxt));
  if (cntx == NULL) return EGADS_MALLOC;
  object = (egObject *) EG_alloc(sizeof(egObject));
  if (object == NULL) {
    EG_free(cntx);
    return EGADS_MALLOC;
  }
  for (i = 0; i < MTESSPARAM; i++) cntx->tess[i] = 0.0;
  cntx->outLevel  = 1;
  cntx->signature = EGADSprop;
  cntx->usrPtr    = NULL;
  cntx->threadID  = EMP_ThreadID();
  cntx->mutex     = EMP_LockCreate();
  cntx->pool      = NULL;
  cntx->last      = object;
  if (cntx->mutex == NULL)
    printf(" EMP Error: mutex creation = NULL (EG_open)!\n");
  
  object->magicnumber = MAGIC;
  object->oclass      = CONTXT;
  object->mtype       = 1;                /* lite version */
  object->tref        = NULL;             /* not used here */
  object->attrs       = NULL;
  object->blind       = cntx;
  object->topObj      = NULL;             /* our single model */
  object->prev        = NULL;
  object->next        = NULL;

  EG_exactInit();
  *context = object;
  return EGADS_SUCCESS;
}


int
EG_referenceObject(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


int
EG_referenceTopObj(/*@unused@*/ egObject *object,
                   /*@unused@*/ /*@null@*/ const egObject *ref)
{
  return EGADS_SUCCESS;
}


/*@kept@*/ /*@null@*/ egObject *
EG_context(const egObject *obj)
{
  int      cnt;
  egObject *object, *topObj;
  
  if (obj == NULL) {
    printf(" EGADS Internal: EG_context called with NULL!\n");
    return NULL;
  }
  if (obj->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context Object NOT an ego!\n");
    return NULL;
  }
  if (obj->oclass == CONTXT) return (egObject *) obj;
  
  object = obj->topObj;
  if (object == NULL) {
    printf(" EGADS Internal: EG_context topObj is NULL!\n");
    return NULL;
  }
  if (object->magicnumber != MAGIC) {
    printf(" EGADS Internal: EG_context topObj NOT an ego!\n");
    return NULL;
  }
  if (object->oclass == CONTXT) return object;
  
  cnt = 0;
  do {
    topObj = object->topObj;
    if (topObj == NULL) {
      printf(" EGADS Internal: %d EG_context contents of topObj is NULL!\n",
             cnt);
      return NULL;
    }
    if (topObj->magicnumber != MAGIC) {
      printf(" EGADS Internal: %d EG_context contents of topObj NOT an ego!\n",
             cnt);
      return NULL;
    }
    if (topObj->oclass == CONTXT) return topObj;
    object = topObj;
    cnt++;
  } while (object != NULL);
  
  printf(" EGADS Internal: Cannot find context -- depth = %d!\n", cnt);
  EG_traceback();
  return NULL;
}


int
EG_sameThread(const egObject *obj)
{
  egObject *context;
  egCntxt  *cntxt;
  
  if (obj == NULL)               return 1;
  if (obj->magicnumber != MAGIC) return 1;
  context = EG_context(obj);
  if (context == NULL)           return 1;
  
  cntxt = (egCntxt *) context->blind;
  if (cntxt->threadID == EMP_ThreadID()) return 0;
  return 1;
}


int
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


int
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


int
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
  int    status;
  size_t nbytes, ntest;
  char   *stream;
  FILE   *fp;

  *model = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;

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

  *model = context->topObj;
  return EGADS_SUCCESS;
}


int
EG_deleteObject(egObject *object)
{
  int      i, j;
  egObject *pobj, *nobj, *context;
  egCntxt  *cntx;
  egTessel *tess;
  
  if (object == NULL)               return EGADS_NULLOBJ;
  if (object->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (object->oclass == EMPTY)      return EGADS_EMPTY;
  if (object->oclass == REFERENCE)  return EGADS_REFERCE;

  if (object->oclass == TESSELLATION) {
    context = EG_context(object);
    if (context == NULL)            return EGADS_NOTCNTX;
    cntx = (egCntxt *) context->blind;
    if (cntx == NULL)               return EGADS_NODATA;
    if (cntx->mutex != NULL) EMP_LockSet(cntx->mutex);
    tess = (egTessel *) object->blind;
    if (tess != NULL) {
      if (tess->xyzs != NULL) EG_free(tess->xyzs);
      if (tess->tess1d != NULL) {
        for (i = 0; i < tess->nEdge; i++) {
          if (tess->tess1d[i].faces[0].faces != NULL)
            EG_free(tess->tess1d[i].faces[0].faces);
          if (tess->tess1d[i].faces[1].faces != NULL)
            EG_free(tess->tess1d[i].faces[1].faces);
          if (tess->tess1d[i].faces[0].tric  != NULL)
            EG_free(tess->tess1d[i].faces[0].tric);
          if (tess->tess1d[i].faces[1].tric  != NULL)
            EG_free(tess->tess1d[i].faces[1].tric);
          if (tess->tess1d[i].xyz    != NULL)
            EG_free(tess->tess1d[i].xyz);
          if (tess->tess1d[i].t      != NULL)
            EG_free(tess->tess1d[i].t);
          if (tess->tess1d[i].global != NULL)
            EG_free(tess->tess1d[i].global);
        }
        EG_free(tess->tess1d);
      }
      if (tess->tess2d != NULL) {
        for (i = 0; i < 2*tess->nFace; i++) {
          if (tess->tess2d[i].mKnots != NULL)
            EG_deleteObject(tess->tess2d[i].mKnots);
          if (tess->tess2d[i].xyz    != NULL)
            EG_free(tess->tess2d[i].xyz);
          if (tess->tess2d[i].uv     != NULL)
            EG_free(tess->tess2d[i].uv);
          if (tess->tess2d[i].global != NULL)
            EG_free(tess->tess2d[i].global);
          if (tess->tess2d[i].ptype  != NULL)
            EG_free(tess->tess2d[i].ptype);
          if (tess->tess2d[i].pindex != NULL)
            EG_free(tess->tess2d[i].pindex);
          if (tess->tess2d[i].bary   != NULL)
            EG_free(tess->tess2d[i].bary);
          if (tess->tess2d[i].frame != NULL)
            EG_free(tess->tess2d[i].frame);
          if (tess->tess2d[i].frlps != NULL)
            EG_free(tess->tess2d[i].frlps);
          if (tess->tess2d[i].tris   != NULL)
            EG_free(tess->tess2d[i].tris);
          if (tess->tess2d[i].tric   != NULL)
            EG_free(tess->tess2d[i].tric);
          if (tess->tess2d[i].patch  != NULL) {
            for (j = 0; j < tess->tess2d[i].npatch; j++) {
              if (tess->tess2d[i].patch[j].ipts != NULL)
                EG_free(tess->tess2d[i].patch[j].ipts);
              if (tess->tess2d[i].patch[j].bounds != NULL)
                EG_free(tess->tess2d[i].patch[j].bounds);
            }
            EG_free(tess->tess2d[i].patch);
          }
        }
        EG_free(tess->tess2d);
      }
      if (tess->globals != NULL) EG_free(tess->globals);
      EG_free(tess);
      object->oclass = EMPTY;
      object->blind  = NULL;
      
      /* patch up the lists & put the object in the pool */
      pobj = object->prev;          /* always have a previous -- context! */
      nobj = object->next;
      if (nobj == NULL) {
        if (object != cntx->last)
          printf(" EGADS Info: Context Last NOT Object Next w/ NULL!\n");
        cntx->last = pobj;
      } else {
        nobj->prev = pobj;
      }
      if (pobj == NULL) {
        printf(" EGADS Info: PrevObj is NULL (EG_destroyObject)!\n");
      } else {
        pobj->next = nobj;
      }
      object->prev = NULL;
      object->next = cntx->pool;
      cntx->pool   = object;
      if (cntx->mutex != NULL) EMP_LockRelease(cntx->mutex);
    }
    return EGADS_SUCCESS;
  }

  /* report other deletes? */
  
  
  return EGADS_SUCCESS;
}


int
EG_close(egObject *context)
{
  int      i, status;
  egAttrs  *attrs;
  egObject *obj, *next, *last;
  egCntxt  *cntx;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  cntx = (egCntxt *) context->blind;
  if (cntx == NULL)                  return EGADS_NODATA;

  /* delete tessellation objects */
  
  obj  = context->next;
  last = NULL;
  while (obj != NULL) {
    next = obj->next;
    if (obj->oclass == TESSELLATION)
      if (EG_deleteObject(obj) == EGADS_SUCCESS) {
        obj = last;
        if (obj == NULL) {
          next = context->next;
        } else {
          next = obj->next;
        }
      }
    last = obj;
    obj  = next;
  }

  /* delete all objects */
  if (cntx->mutex != NULL) EMP_LockSet(cntx->mutex);
  
  obj = context->next;
  while (obj != NULL) {
    next = obj->next;
    if ((obj->oclass >= PCURVE) && (obj->oclass <= MODEL)) {
      status = EG_freeBlind(obj);
      if (status != EGADS_SUCCESS)
        printf(" EGADS Info: freeBlind = %d in Cleanup (EG_close)!\n", status);
    }
    attrs = (egAttrs *) obj->attrs;
    if (attrs != NULL) {
      for (i = 0; i < attrs->nattrs; i++) {
        EG_free(attrs->attrs[i].name);
        if (attrs->attrs[i].type == ATTRINT) {
          if (attrs->attrs[i].length > 1) EG_free(attrs->attrs[i].vals.integers);
        } else if ((attrs->attrs[i].type == ATTRREAL) ||
                   (attrs->attrs[i].type == ATTRCSYS)) {
          if (attrs->attrs[i].length > 1) EG_free(attrs->attrs[i].vals.reals);
        } else if (attrs->attrs[i].type == ATTRSTRING) {
          EG_free(attrs->attrs[i].vals.string);
        }
      }
      EG_free(attrs->attrs);
      EG_free(attrs);
    }
    EG_free(obj);
    obj = next;
  }
  
  /* clean up the pool */
  
  obj = cntx->pool;
  while (obj != NULL) {
    if (obj->magicnumber != MAGIC) {
      printf(" EGADS Info: Found BAD Object in Cleanup (EG_close)!\n");
      printf("             Class = %d\n", obj->oclass);
      break;
    }
    next = obj->next;
    EG_free(obj);
    obj = next;
  }

  context->magicnumber = 0;
  context->oclass      = EMPTY;
  EG_free(context);
  if (cntx->mutex != NULL) EMP_LockRelease(cntx->mutex);
  if (cntx->mutex != NULL) EMP_LockDestroy(cntx->mutex);
  EG_free(cntx);
  
  return EGADS_SUCCESS;
}
