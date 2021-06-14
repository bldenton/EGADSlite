/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             FORTRAN Bindings for Base & Effective Topo Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INT8 unsigned long long

#include "egadsTypes.h"
#include "egadsInternals.h"


  extern void EG_revision(int *major, int *minor, const char **OCCrev);
  extern int  EG_open(egObject **context);
  extern int  EG_updateThread(egObject *context);
  extern int  EG_deleteObject(egObject *object);
  extern int  EG_makeTransform(egObject *context, const double *xform, 
                               egObject **oform);
  extern int  EG_getTransformation(const egObject *oform, double *xform);
  extern int  EG_setOutLevel(egObject *context, int outLevel);
  extern int  EG_setFixedKnots(egObject *context, int fixed);
  extern int  EG_setFullAttrs(egObject *context, int full);
  extern int  EG_setTessParam(egObject *context, int iParam, double value,
                             double *oldValue);
  extern int  EG_getContext(egObject *object, egObject **context);
  extern int  EG_getInfo(const egObject *object, int *oclass, int *mtype, 
                         egObject **top, egObject **prev, egObject **next);
  extern int  EG_copyObject(const egObject *object, /*@null@*/ void *ptr,
                            egObject **copy);
  extern int  EG_flipObject(const egObject *object, egObject **copy);
  extern int  EG_close(egObject *context);
  extern int  EG_initEBody(egObject *tess, double angle, egObject **EBody);
  extern int  EG_finishEBody(egObject *EBody);
  extern int  EG_makeEFace(egObject *EBody, int nFace, /*@null@*/
                           const egObject **Faces, egObject **EFace);
  extern int  EG_makeAttrEFaces(egObject *EBody, const char *attrName,
                                int *nEFace, egObject ***EFaces);
  extern int  EG_effectiveMap(egObject *EObject, double *eparam,
                              egObject **Object, double *param);
  extern int  EG_effectiveEdgeList(egObject *EEdge, int *nedges,
                                   egObject ***edges, int **senses,
                                   double **tstart);



void EG_c2f(/*@null@*/ const char *string, char *name, int nameLen)
{
  int  i, len;

  for (i = 0; i < nameLen; i++) name[i] = ' ';
  if (string == NULL) return;
  len = strlen(string);
  if (len > nameLen) len = nameLen;
  for (i = 0; i < len; i++) name[i] = string[i];
}


/*@null@*/ char *EG_f2c(const char *name, int nameLen)
{
  char *string;
  int  i, len;

  for (len = nameLen; len >= 1; len--) if (name[len-1] != ' ') break;
  string = (char *) EG_alloc((len+1)*sizeof(char));
  if (string != NULL) {
    for (i = 0; i < len; i++) string[i] = name[i];
    string[len] = 0;
  }
  return string;
}


void
#ifdef WIN32
IG_REVISION (int *major, int *minor, char *str, int strLen)
#else
ig_revision_(int *major, int *minor, char *str, int strLen)
#endif
{
  const char *ostr;
  
  EG_revision(major, minor, &ostr);
  EG_c2f(ostr, str, strLen);
}


int
#ifdef WIN32
IG_OPEN (INT8 *cntxt)
#else
ig_open_(INT8 *cntxt)
#endif
{
  int      stat;
  egObject *context;

  *cntxt = 0;
  stat   = EG_open(&context);
  if (stat == EGADS_SUCCESS) *cntxt = (INT8) context;
  return stat;
}


int
#ifdef WIN32
IG_UPDATETHREAD (INT8 *obj)
#else
ig_updatethread_(INT8 *obj)
#endif
{
  egObject *object;
  
  object = (egObject *) *obj;
  return EG_updateThread(object);
}


int
#ifdef WIN32
IG_DELETEOBJECT (INT8 *obj)
#else
ig_deleteobject_(INT8 *obj)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_deleteObject(object);
}


int
#ifdef WIN32
IG_MAKETRANSFORM (INT8 *cntxt, double *xform, INT8 *ofrm)
#else
ig_maketransform_(INT8 *cntxt, double *xform, INT8 *ofrm)
#endif
{
  int      stat;
  egObject *context, *oform;

  *ofrm   = 0;
  context = (egObject *) *cntxt;
  stat    = EG_makeTransform(context, xform, &oform);
  if (stat == EGADS_SUCCESS) *ofrm = (INT8) oform;
  return stat;
}


int
#ifdef WIN32
IG_GETTRANSFORM (INT8 *ofrm, double *xform)
#else
ig_gettransform_(INT8 *ofrm, double *xform)
#endif
{
  egObject *oform;

  oform = (egObject *) *ofrm;
  return EG_getTransformation(oform, xform);
}


int
#ifdef WIN32
IG_SETOUTLEVEL (INT8 *cntxt, int *out)
#else
ig_setoutlevel_(INT8 *cntxt, int *out)
#endif
{
  egObject *context;

  context = (egObject *) *cntxt;
  return EG_setOutLevel(context, *out);
}


int
#ifdef WIN32
IG_SETFIXEDKNOTS (INT8 *cntxt, int *out)
#else
ig_setfixedknots_(INT8 *cntxt, int *out)
#endif
{
  egObject *context;

  context = (egObject *) *cntxt;
  return EG_setFixedKnots(context, *out);
}


int
#ifdef WIN32
IG_SETFULLATTRS (INT8 *cntxt, int *out)
#else
ig_setfullattrs_(INT8 *cntxt, int *out)
#endif
{
  egObject *context;

  context = (egObject *) *cntxt;
  return EG_setFullAttrs(context, *out);
}


int
#ifdef WIN32
IG_SETTESSPARAM (INT8 *cntxt, int *iparam, double *val, double *oldval)
#else
ig_settessparam_(INT8 *cntxt, int *iparam, double *val, double *oldval)
#endif
{
  egObject *context;
  
  context = (egObject *) *cntxt;
  return EG_setTessParam(context, *iparam, *val, oldval);
}


int
#ifdef WIN32
IG_GETCONTEXT (INT8 *obj, INT8 *cntxt)
#else
ig_getcontext_(INT8 *obj, INT8 *cntxt)
#endif
{
  int      stat;
  egObject *object, *context;
  
  *cntxt = 0;
  object = (egObject *) *obj;
  stat   = EG_getContext(object, &context);
  if (stat == EGADS_SUCCESS) *cntxt = (INT8) context;
  return stat;
}


int
#ifdef WIN32
IG_GETINFO (INT8 *obj, int *oclass, int *mtype, INT8 *top, INT8 *prv, INT8 *nxt)
#else
ig_getinfo_(INT8 *obj, int *oclass, int *mtype, INT8 *top, INT8 *prv, INT8 *nxt)
#endif
{
  int      stat;
  egObject *object, *topObj, *prev, *next;
  
  *prv   = *nxt = 0;
  object = (egObject *) *obj;
  stat   = EG_getInfo(object, oclass, mtype, &topObj, &prev, &next);
  if (stat == EGADS_SUCCESS) {
    *top = (INT8) topObj;
    *prv = (INT8) prev;
    *nxt = (INT8) next;
  }
  return stat;
}


int
#ifdef WIN32
IG_COPYOBJECT (INT8 *obj, INT8 *ofrm, INT8 *cp)
#else
ig_copyobject_(INT8 *obj, INT8 *ofrm, INT8 *cp)
#endif
{
  int      stat;
  egObject *object, *xform, *copy;

  *cp    = 0;
  object = (egObject *) *obj;
  xform  = (egObject *) *ofrm;
  stat   = EG_copyObject(object, xform, &copy);
  if (stat == EGADS_SUCCESS) *cp = (INT8) copy;
  return stat;
}


int
#ifdef WIN32
IG_FLIPOBJECT (INT8 *obj, INT8 *cp)
#else
ig_flipobject_(INT8 *obj, INT8 *cp)
#endif
{
  int      stat;
  egObject *object, *copy;

  *cp    = 0;
  object = (egObject *) *obj;
  stat   = EG_flipObject(object, &copy);
  if (stat == EGADS_SUCCESS) *cp = (INT8) copy;
  return stat;
}


int
#ifdef WIN32
IG_CLOSE (INT8 *obj)
#else
ig_close_(INT8 *obj)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_close(object);
}


int
#ifdef WIN32
IG_INITEBODY (INT8 *tess, double *angle, INT8 *ebody)
#else
ig_initebody_(INT8 *tess, double *angle, INT8 *ebody)
#endif
{
  int      stat;
  egObject *object, *body;
  
  *ebody = 0;
  object = (egObject *) *tess;
  stat   = EG_initEBody(object, *angle, &body);
  if (stat == EGADS_SUCCESS) *ebody = (INT8) body;
  return stat;
}


int
#ifdef WIN32
IG_FINISHEBODY (INT8 *obj)
#else
ig_finishebody_(INT8 *obj)
#endif
{
  egObject *object;

  object = (egObject *) *obj;
  return EG_finishEBody(object);
}


int
#ifdef WIN32
IG_MAKEEFACE (INT8 *ebody, int *nface, INT8 *faces, INT8 *eface)
#else
ig_makeeface_(INT8 *ebody, int *nface, INT8 *faces, INT8 *eface)
#endif
{
  int      i, stat;
  egObject *object, *oface, **objs = NULL;
  
  *eface = 0;
  object = (egObject *) *ebody;
  if (*nface > 0) {
    objs = (egObject **) EG_alloc(*nface*sizeof(egObject *));
    if (objs == NULL) return EGADS_MALLOC;
    for (i = 0; i < *nface; i++)
      objs[i] = (egObject *) faces[i];
  }
  stat = EG_makeEFace(object, *nface, (const egObject **) objs, &oface);
  if (objs != NULL) EG_free(objs);
  if (stat == EGADS_SUCCESS) *eface = (INT8) oface;
  return stat;
}


int
#ifdef WIN32
IG_MAKEATTREFACES (INT8 *ebody, const char *attr, int *nface, INT8 **efaces,
                   int attrLen)
#else
ig_makeattrefaces_(INT8 *ebody, const char *attr, int *nface, INT8 **efaces,
                   int attrLen)
#endif
{
  int      i, stat;
  char     *fattr;
  INT8     *cobjs;
  egObject *object, **objs;
  
  *efaces = 0;
  *nface  = 0;
  object  = (egObject *) *ebody;
  fattr   = EG_f2c(attr, attrLen);
  if (fattr == NULL) return EGADS_MALLOC;
  stat = EG_makeAttrEFaces(object, fattr, nface, &objs);
  EG_free(fattr);
  if (stat != EGADS_SUCCESS) return stat;
  EG_free(fattr);
  *efaces = cobjs = (INT8 *) EG_alloc(*nface*sizeof(INT8));
  if (cobjs == NULL) {
    EG_free(objs);
    return EGADS_SUCCESS;
  }
  for (i = 0; i < *nface; i++) cobjs[i] = (INT8) objs[i];
  EG_free(objs);
  
  return EGADS_SUCCESS;
}


int
#ifdef WIN32
IG_EFFECTIVEMAP (INT8 *eobj, double *eparam, INT8 *obj, double *param)
#else
ig_effectivemap_(INT8 *eobj, double *eparam, INT8 *obj, double *param)
#endif
{
  int      stat;
  egObject *eobject, *object;
  
  *obj    = 0;
  eobject = (egObject *) *eobj;
  stat    = EG_effectiveMap(eobject, eparam, &object, param);
  if (stat == EGADS_SUCCESS) *obj = (INT8) object;
  return stat;
}


int
#ifdef WIN32
IG_EFFECTIVEEDGELIST (INT8 *eobj, int *nedges, INT8 **edges, int **senses,
                      double **tstart)
#else
ig_effectiveedgelist_(INT8 *eobj, int *nedges, INT8 **edges, int **senses,
                      double **tstart)
#endif
{
  int      i, stat;
  INT8     *cobjs;
  egObject *eobject, **edgos;
  
  *edges  = NULL;
  *senses = NULL;
  *tstart = NULL;
  eobject = (egObject *) *eobj;
  stat    = EG_effectiveEdgeList(eobject, nedges, &edgos, senses, tstart);
  if (stat != EGADS_SUCCESS) return stat;
  
  if (*nedges > 0) {
    *edges = cobjs = (INT8 *) EG_alloc(*nedges*sizeof(INT8));
    if (cobjs == NULL) {
      EG_free(edgos);
      EG_free(senses);
      EG_free(tstart);
      *senses = NULL;
      *tstart = NULL;
      return EGADS_MALLOC;
    }
    for (i = 0; i < *nedges; i++) cobjs[i] = (INT8) edgos[i];
  }
  EG_free(edgos);
  
  return EGADS_SUCCESS;
}
