/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Import a Model from EGADS (via a string)
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteClasses.h"


//#define DEBUG


  extern int EG_close(egObject *context);
  extern int EG_evaluate(const egObject *geom, const double *param, double *ev);
  extern int EG_attributeRet(const egObject *obj, const char *name, int *atype,
                             int *len, /*@null@*/ const int **ints,
                             /*@null@*/ const double **reals,
                             /*@null@*/ const char **str);

typedef struct {
  void    *data;
  size_t  ptr;
  size_t  size;
  int     swap;
} stream_T;



static void
swap(void *buffer, size_t size)
{
  char *buf, save;
  
  buf = buffer;
  if (size == 2) {
    save   = buf[1];
    buf[1] = buf[0];
    buf[0] = save;
  } else if (size == 4) {
    save   = buf[3];
    buf[3] = buf[0];
    buf[0] = save;
    save   = buf[2];
    buf[2] = buf[1];
    buf[1] = save;
  } else {
    save   = buf[7];
    buf[7] = buf[0];
    buf[0] = save;
    save   = buf[6];
    buf[6] = buf[1];
    buf[1] = save;
    save   = buf[5];
    buf[5] = buf[2];
    buf[2] = save;
    save   = buf[4];
    buf[4] = buf[3];
    buf[3] = save;
  }
}


static int
Fread(void *data, size_t size, int nitems, stream_T *stream)
{
  int  i;
  char *buf;
  
  buf = data;
  memcpy(data, &(((char *) stream->data)[stream->ptr]), size*nitems);
  if ((size != sizeof(char)) && (stream->swap == 1))
    for (i = 0; i < nitems; i++) swap(&buf[i*size], size);
  stream->ptr += size*nitems;

  return nitems;
}


static int
EG_addStrAttr(egObject *obj, const char *name, const char *str)
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
  length = strlen(name);
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
      if (strcmp(attrs->attrs[i].name,name) == 0) {
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
    attrs->attrs[find].length = strlen(attrs->attrs[find].vals.string);

  return EGADS_SUCCESS;
}


static void
EG_freeAttrs(egAttrs **attrx)
{
  int     i;
  egAttrs *attrs;
  
  attrs = *attrx;
  if (attrs == NULL) return;
  
  *attrx = NULL;
  /* remove any attributes */
  for (i = 0; i < attrs->nattrs; i++) {
    if (attrs->attrs[i].name != NULL) EG_free(attrs->attrs[i].name);
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


static int
EG_readString(stream_T *fp, char **string)
{
  int    len;
  size_t n;
  
  *string = NULL;
  n = Fread(&len, sizeof(int), 1, fp);
  if (n   != 1) return EGADS_READERR;
  if (len <  0) return EGADS_READERR;
  if (len == 0) return EGADS_SUCCESS;
  
  *string = (char *) EG_alloc(len*sizeof(char));
  if (*string == NULL) return EGADS_MALLOC;
  
  n = Fread(*string, sizeof(char), len, fp);
  if (n != len) {
    EG_free(*string);
    *string = NULL;
    return EGADS_READERR;
  }
  
  return EGADS_SUCCESS;
}


static int
EG_readAttrs(stream_T *fp, egAttrs **attrx)
{
  int     nattr, i, stat;
  size_t  n;
  egAttr  *attr;
  egAttrs *attrs;
  
  *attrx = NULL;
  n = Fread(&nattr, sizeof(int), 1, fp);
  if (n     != 1) return EGADS_READERR;
  if (nattr == 0) return EGADS_SUCCESS;
  
  attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
  if (attrs == NULL) return EGADS_MALLOC;
  attr  = (egAttr *)  EG_alloc(nattr*sizeof(egAttr));
  if (attr == NULL) {
    EG_free(attrs);
    return EGADS_MALLOC;
  }
  attrs->nattrs = nattr;
  attrs->attrs  = attr;
  for (i = 0; i < nattr; i++) {
    attr[i].name   = NULL;
    attr[i].length = 1;
    attr[i].type   = ATTRINT;
  }
  
/*@-mustfreefresh@*/
  for (i = 0; i < nattr; i++) {
    n = Fread(&attr[i].type,   sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    n = Fread(&attr[i].length, sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    stat = EG_readString(fp, &attr[i].name);
    if (stat != EGADS_SUCCESS) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    if (attr[i].type == ATTRINT) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = Fread(&attr[i].vals.integer, sizeof(int),              1, fp);
      } else if (attr[i].length > 1) {
        attr[i].vals.integers = (int *) EG_alloc(attr[i].length*sizeof(int));
        if (attr[i].vals.integers == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread(attr[i].vals.integers, sizeof(int), attr[i].length, fp);
      }
      if (n != attr[i].length) {
        EG_freeAttrs(&attrs);
        return EGADS_READERR;
      }
    } else if ((attr[i].type == ATTRREAL) || (attr[i].type == ATTRCSYS)) {
      n = attr[i].length;
      if (attr[i].length == 1) {
        n = Fread(&attr[i].vals.real, sizeof(double),              1, fp);
      } else if (attr[i].length > 1) {
        attr[i].vals.reals = (double *) EG_alloc(attr[i].length*sizeof(double));
        if (attr[i].vals.reals == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread(attr[i].vals.reals, sizeof(double), attr[i].length, fp);
      }
      if (n != attr[i].length) return EGADS_READERR;
    } else {
      stat = EG_readString(fp, &attr[i].vals.string);
      if (stat != EGADS_SUCCESS) return EGADS_READERR;
    }
  }
/*@+mustfreefresh@*/
  
  *attrx = attrs;
  return EGADS_SUCCESS;
}


static int
EG_readGeometry(liteGeometry *lgeom, int *iref, stream_T *fp)
{
  int n, nhead, ndata;
  
  *iref         = 0;
  lgeom->ref    = NULL;
  lgeom->header = NULL;
  lgeom->data   = NULL;
  n = Fread(iref,   sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&nhead, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ndata, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  
  if (nhead != 0) {
    lgeom->header = (int *) EG_alloc(nhead*sizeof(int));
    if (lgeom->header == NULL) return EGADS_MALLOC;
    n = Fread(lgeom->header, sizeof(int), nhead, fp);
    if (n != nhead) return EGADS_READERR;
  }
  lgeom->data = (double *) EG_alloc(ndata*sizeof(double));
  if (lgeom->data == NULL) return EGADS_MALLOC;
  n = Fread(lgeom->data, sizeof(double), ndata, fp);
  if (n != ndata) return EGADS_READERR;
  
  return EGADS_SUCCESS;
}


static int
EG_readBody(egObject *context, egObject *mobject, int bindex, stream_T *fp)
{
  int          i, j, m, n, stat, mtype, iref, ntypes[8];
  double       t, d, x0[2], x1[2], data[6];
//#ifdef DEBUG
  int          atype, alen;
  const int    *ints;
  const double *reals;
  const char   *str;
//#endif
  egObject     *obj, *bobj, *pcobj;
  liteGeometry *lgeom;
  liteNode     *lnode;
  liteEdge     *ledge;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;

  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (mobject == NULL)               return EGADS_NULLOBJ;
  if (mobject->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject->oclass != MODEL)      return EGADS_NOTMODEL;
  lmodel = (liteModel *) mobject->blind;

  n = Fread(&mtype, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(ntypes, sizeof(int), 8, fp);
  if (n != 8) return EGADS_READERR;
  lbody = (liteBody *) EG_alloc(sizeof(liteBody));
  if (lbody == NULL) return EGADS_MALLOC;
  lbody->pcurves.objs   = NULL;
  lbody->curves.objs    = NULL;
  lbody->surfaces.objs  = NULL;
  lbody->nodes.objs     = NULL;
  lbody->edges.objs     = NULL;
  lbody->loops.objs     = NULL;
  lbody->faces.objs     = NULL;
  lbody->shells.objs    = NULL;
  lbody->senses         = NULL;
  lbody->pcurves.nobjs  = 0;
  lbody->curves.nobjs   = 0;
  lbody->surfaces.nobjs = 0;
  lbody->nodes.nobjs    = 0;
  lbody->edges.nobjs    = 0;
  lbody->loops.nobjs    = 0;
  lbody->faces.nobjs    = 0;
  lbody->shells.nobjs   = 0;
#ifdef DEBUG
  printf(" Reading Body #%d: %d %d %d %d %d %d %d %d\n", bindex+1, ntypes[0],
         ntypes[1], ntypes[2], ntypes[3], ntypes[4], ntypes[5], ntypes[6],
         ntypes[7]);
#endif

  stat = EG_makeObject(context, &bobj);
  if (stat != EGADS_SUCCESS) {
    EG_free(lbody);
    return stat;
  }
  bobj->oclass = BODY;
  bobj->mtype  = mtype;
  bobj->blind  = lbody;
  lmodel->bodies[bindex] = bobj;
  
  /* make all of the objects */
  if (ntypes[0] > 0) {
    lbody->pcurves.objs = (egObject **) EG_alloc(ntypes[0]*sizeof(egObject *));
    if (lbody->pcurves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[0]; i++) lbody->pcurves.objs[i] = NULL;
    lbody->pcurves.nobjs = ntypes[0];
    for (i = 0; i < ntypes[0]; i++) {
      stat = EG_makeObject(context, &lbody->pcurves.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[1] > 0) {
    lbody->curves.objs = (egObject **) EG_alloc(ntypes[1]*sizeof(egObject *));
    if (lbody->curves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[1]; i++) lbody->curves.objs[i] = NULL;
    lbody->curves.nobjs = ntypes[1];
    for (i = 0; i < ntypes[1]; i++) {
      stat = EG_makeObject(context, &lbody->curves.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[2] > 0) {
    lbody->surfaces.objs = (egObject **) EG_alloc(ntypes[2]*sizeof(egObject *));
    if (lbody->surfaces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[2]; i++) lbody->surfaces.objs[i] = NULL;
    lbody->surfaces.nobjs = ntypes[2];
    for (i = 0; i < ntypes[2]; i++) {
      stat = EG_makeObject(context, &lbody->surfaces.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[3] > 0) {
    lbody->nodes.objs = (egObject **) EG_alloc(ntypes[3]*sizeof(egObject *));
    if (lbody->nodes.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[3]; i++) lbody->nodes.objs[i] = NULL;
    lbody->nodes.nobjs = ntypes[3];
    for (i = 0; i < ntypes[3]; i++) {
      stat = EG_makeObject(context, &lbody->nodes.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[4] > 0) {
    lbody->edges.objs = (egObject **) EG_alloc(ntypes[4]*sizeof(egObject *));
    if (lbody->edges.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[4]; i++) lbody->edges.objs[i] = NULL;
    lbody->edges.nobjs = ntypes[4];
    for (i = 0; i < ntypes[4]; i++) {
      stat = EG_makeObject(context, &lbody->edges.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[5] > 0) {
    lbody->loops.objs = (egObject **) EG_alloc(ntypes[5]*sizeof(egObject *));
    if (lbody->loops.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[5]; i++) lbody->loops.objs[i] = NULL;
    lbody->loops.nobjs = ntypes[5];
    for (i = 0; i < ntypes[5]; i++) {
      stat = EG_makeObject(context, &lbody->loops.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[6] > 0) {
    lbody->faces.objs = (egObject **) EG_alloc(ntypes[6]*sizeof(egObject *));
    if (lbody->faces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[6]; i++) lbody->faces.objs[i] = NULL;
    lbody->faces.nobjs = ntypes[6];
    for (i = 0; i < ntypes[6]; i++) {
      stat = EG_makeObject(context, &lbody->faces.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }
  if (ntypes[7] > 0) {
    lbody->shells.objs = (egObject **) EG_alloc(ntypes[7]*sizeof(egObject *));
    if (lbody->shells.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[7]; i++) lbody->shells.objs[i] = NULL;
    lbody->shells.nobjs = ntypes[7];
    for (i = 0; i < ntypes[7]; i++) {
      stat = EG_makeObject(context, &lbody->shells.objs[i]);
      if (stat != EGADS_SUCCESS) return stat;
    }
  }

#ifdef DEBUG
  printf(" Reading %d PCurves...\n", ntypes[0]);
#endif
  /* pcurves */
  if (lbody->pcurves.objs != NULL)
    for (i = 0; i < lbody->pcurves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lgeom = (liteGeometry *) EG_alloc(sizeof(liteGeometry));
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        if (lgeom->header != NULL) EG_free(lgeom->header);
        if (lgeom->data   != NULL) EG_free(lgeom->data);
        EG_free(lgeom);
        return stat;
      }
      if (iref != 0) lgeom->ref = lbody->pcurves.objs[iref-1];
      obj = lbody->pcurves.objs[i];
      obj->oclass = PCURVE;
      obj->mtype  = mtype;
      obj->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      if (mtype != BSPLINE) continue;
      for (n = 2; n < lgeom->header[2]; n++) {
        m     = lgeom->header[3] + 2*n - 4;
        x0[0] = lgeom->data[m+2] - lgeom->data[m  ];
        x0[1] = lgeom->data[m+3] - lgeom->data[m+1];
        d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
        if (d != 0.0) {
          x0[0] /= d;
          x0[1] /= d;
        }
        x1[0] = lgeom->data[m+4] - lgeom->data[m+2];
        x1[1] = lgeom->data[m+5] - lgeom->data[m+3];
        d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
        if (d != 0.0) {
          x1[0] /= d;
          x1[1] /= d;
        }
        d = x0[0]*x1[0] + x0[1]*x1[1];
        if (d < -0.95) {
#ifdef DEBUG
          printf(" EGADS Info: PCurve %d dot flip at %d/%d (%lf) -- %lf %lf!\n",
                 i, n-2, lgeom->header[2]-2, d,
                 lgeom->data[0], lgeom->data[lgeom->header[3]-1]);
#endif
          stat = EG_addStrAttr(obj, ".Bad", "CPrev");
          if (stat != EGADS_SUCCESS)
            printf("             EG_addStrAttr CPrev= %d\n", stat);
        }
      }
    }
  
  /* curves */
#ifdef DEBUG
  printf(" Reading %d Curves...\n", ntypes[1]);
#endif
  if (lbody->curves.objs != NULL)
    for (i = 0; i < lbody->curves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lgeom = (liteGeometry *) EG_alloc(sizeof(liteGeometry));
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        if (lgeom->header != NULL) EG_free(lgeom->header);
        if (lgeom->data   != NULL) EG_free(lgeom->data);
        EG_free(lgeom);
        return stat;
      }
      if (iref != 0) lgeom->ref = lbody->curves.objs[iref-1];
      obj = lbody->curves.objs[i];
      obj->oclass = CURVE;
      obj->mtype  = mtype;
      obj->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }
  
  /* surfaces */
#ifdef DEBUG
  printf(" Reading %d Surfaces...\n", ntypes[2]);
#endif
  if (lbody->surfaces.objs != NULL)
    for (i = 0; i < lbody->surfaces.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lgeom = (liteGeometry *) EG_alloc(sizeof(liteGeometry));
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        if (lgeom->header != NULL) EG_free(lgeom->header);
        if (lgeom->data   != NULL) EG_free(lgeom->data);
        EG_free(lgeom);
        return stat;
      }
/*@-nullderef@*/
      if (iref < 0) lgeom->ref = lbody->curves.objs[-iref-1];
/*@+nullderef@*/
      if (iref > 0) lgeom->ref = lbody->surfaces.objs[iref-1];
      obj = lbody->surfaces.objs[i];
      obj->oclass = SURFACE;
      obj->mtype  = mtype;
      obj->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }
  
  /* nodes */
#ifdef DEBUG
  printf(" Reading %d Nodes...\n", ntypes[3]);
#endif
  if (lbody->nodes.objs != NULL)
    for (i = 0; i < lbody->nodes.nobjs; i++) {
      lnode = (liteNode *) EG_alloc(sizeof(liteNode));
      if (lnode == NULL) return EGADS_MALLOC;
      n = Fread(lnode->xyz,  sizeof(double), 3, fp);
      if (n != 3) {
        EG_free(lnode);
        return EGADS_READERR;
      }
      n = Fread(&lnode->tol, sizeof(double), 1, fp);
      if (n != 1) {
        EG_free(lnode);
        return EGADS_READERR;
      }
      obj = lbody->nodes.objs[i];
      obj->oclass = NODE;
      obj->mtype  = 0;
      obj->blind  = lnode;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }
  
  /* edges */
#ifdef DEBUG
  printf(" Reading %d Edges...\n", ntypes[4]);
#endif
  if (lbody->edges.objs != NULL)
    for (i = 0; i < lbody->edges.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      ledge = (liteEdge *) EG_alloc(sizeof(liteEdge));
      if (ledge == NULL) return EGADS_MALLOC;
      ledge->curve    = NULL;
      ledge->nodes[0] = NULL;
      ledge->nodes[1] = NULL;
  
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) ledge->curve = lbody->curves.objs[iref-1];
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) ledge->nodes[0] = lbody->nodes.objs[iref-1];
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) ledge->nodes[1] = lbody->nodes.objs[iref-1];
/*@+nullderef@*/
     
      n = Fread(ledge->trange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_free(ledge);
        return EGADS_READERR;
      }
      n = Fread(ledge->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_free(ledge);
        return EGADS_READERR;
      }
      n = Fread(&ledge->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_free(ledge);
        return EGADS_READERR;
      }

      obj = lbody->edges.objs[i];
      obj->oclass = EDGE;
      obj->mtype  = mtype;
      obj->blind  = ledge;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }

  /* loops */
#ifdef DEBUG
  printf(" Reading %d Loops...\n", ntypes[5]);
#endif
  if (lbody->loops.objs != NULL)
    for (i = 0; i < lbody->loops.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lloop = (liteLoop *) EG_alloc(sizeof(liteLoop));
      if (lloop == NULL) return EGADS_MALLOC;
      lloop->nedges  = m;
      lloop->surface = NULL;
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(lloop);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
        lloop->surface = lbody->surfaces.objs[iref-1];
/*@+nullderef@*/
        m *= 2;
      }
      n = Fread(lloop->bbox,  sizeof(double), 6, fp);
      if (n != 6) {
        EG_free(lloop);
        return EGADS_READERR;
      }
      lloop->senses = (int *) EG_alloc(lloop->nedges*sizeof(int));
      if (lloop->senses == NULL) {
        EG_free(lloop);
        return EGADS_MALLOC;
      }
      n = Fread(lloop->senses,  sizeof(int), lloop->nedges, fp);
      if (n != lloop->nedges) {
        EG_free(lloop->senses);
        EG_free(lloop);
        return EGADS_READERR;
      }
      lloop->edges = (egObject **) EG_alloc(m*sizeof(egObject *));
      if (lloop->edges == NULL) {
        EG_free(lloop->senses);
        EG_free(lloop);
        return EGADS_MALLOC;
      }
      for (j = 0; j < m; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(lloop->edges);
          EG_free(lloop->senses);
          EG_free(lloop);
          return EGADS_READERR;
        }
        if (j < lloop->nedges) {
/*@-nullderef@*/
          lloop->edges[j] = lbody->edges.objs[iref-1];
/*@+nullderef@*/
        } else {
/*@-nullderef@*/
          lloop->edges[j] = lbody->pcurves.objs[iref-1];
/*@+nullderef@*/
        }
      }
      if (lloop->surface != NULL) {
        for (n = 0; n < lloop->nedges; n++) {
          ledge = lloop->edges[n]->blind;
          pcobj = lloop->edges[lloop->nedges+n];
          stat  = EG_attributeRet(pcobj, ".Bad", &atype, &alen,
                                  &ints, &reals, &str);
          if (stat != EGADS_SUCCESS) continue;
          EG_evaluate(pcobj, &ledge->trange[0], data);
          d     = sqrt(data[2]*data[2] + data[3]*data[3]);
          x0[0] = x0[1] = 0.0;
          if (d != 0.0) {
            x0[0] = data[2]/d;
            x0[1] = data[3]/d;
          }
          for (j = 1; j < 1000; j++) {
            t    = ledge->trange[0]+j*(ledge->trange[1]-ledge->trange[0])/999.;
            EG_evaluate(pcobj, &t, data);
            d     = sqrt(data[2]*data[2] + data[3]*data[3]);
            x1[0] = x1[1] = 0.0;
            if (d != 0.0) {
              x1[0] = data[2]/d;
              x1[1] = data[3]/d;
            }
            if (x0[0]*x1[0] + x0[1]*x1[1] < -0.95) {
              stat = EG_addStrAttr(pcobj, ".Bad", "fold");
              if (stat != EGADS_SUCCESS)
                printf(" EGADS Info: EG_addStrAttr fold= %d\n", stat);
            }
            x0[0] = x1[0];
            x0[1] = x1[1];
          }
        }
      }
      
      obj = lbody->loops.objs[i];
      obj->oclass = LOOP;
      obj->mtype  = mtype;
      obj->blind  = lloop;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }
  
  /* faces */
#ifdef DEBUG
  printf(" Reading %d Faces...\n", ntypes[6]);
#endif
  if (lbody->faces.objs != NULL)
    for (i = 0; i < lbody->faces.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lface = (liteFace *) EG_alloc(sizeof(liteFace));
      if (lface == NULL) return EGADS_MALLOC;
      lface->nloops  = m;
      lface->surface = NULL;
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_free(lface);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
        lface->surface = lbody->surfaces.objs[iref-1];
#endif
/*@+nullderef@*/
      }
      n = Fread(lface->urange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_free(lface);
        return EGADS_READERR;
      }
      n = Fread(lface->vrange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_free(lface);
        return EGADS_READERR;
      }
      n = Fread(lface->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_free(lface);
        return EGADS_READERR;
      }
      n = Fread(&lface->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_free(lface);
        return EGADS_READERR;
      }
      lface->senses = (int *) EG_alloc(lface->nloops*sizeof(int));
      if (lface->senses == NULL) {
        EG_free(lface);
        return EGADS_READERR;
      }
      n = Fread(lface->senses,  sizeof(int), lface->nloops, fp);
      if (n != lface->nloops) {
        EG_free(lface->senses);
        EG_free(lface);
        return EGADS_READERR;
      }
      lface->loops = (egObject **) EG_alloc(lface->nloops*sizeof(egObject *));
      if (lface->loops == NULL) {
        EG_free(lface->senses);
        EG_free(lface);
        return EGADS_MALLOC;
      }
      for (j = 0; j < lface->nloops; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(lface->loops);
          EG_free(lface->senses);
          EG_free(lface);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        lface->loops[j] = lbody->loops.objs[iref-1];
#endif
/*@+nullderef@*/
      }

      obj = lbody->faces.objs[i];
      obj->oclass = FACE;
      obj->mtype  = mtype;
      obj->blind  = lface;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
//#ifdef DEBUG
      for (j = 0; j < lface->nloops; j++) {
        lloop = lface->loops[j]->blind;
        if (lloop->surface == NULL) continue;
        for (n = 0; n < lloop->nedges; n++) {
          stat = EG_attributeRet(lloop->edges[lloop->nedges+n], ".Bad", &atype,
                                 &alen, &ints, &reals, &str);
          if ((stat == EGADS_SUCCESS) && (atype == ATTRSTRING))
            if (strcmp(str, "fold") == 0)
              printf(" EGADS Info: Body %d Face %d Loop#%d/Edge#%d Bad PCurve -- %s!\n",
                     bindex+1, i+1, j+1, n+1, str);
        }
      }
//#endif
    }

  /* shells */
#ifdef DEBUG
  printf(" Reading %d Shells...\n", ntypes[7]);
#endif
  if (lbody->shells.objs != NULL)
    for (i = 0; i < lbody->shells.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      lshell = (liteShell *) EG_alloc(sizeof(liteShell));
      if (lshell == NULL) return EGADS_MALLOC;
      lshell->nfaces = m;
      n = Fread(lshell->bbox, sizeof(double), 6, fp);
      if (n != 6) {
        EG_free(lshell);
        return EGADS_READERR;
      }
      lshell->faces = (egObject **) EG_alloc(m*sizeof(egObject *));
      if (lshell->faces == NULL) {
        EG_free(lshell);
        return EGADS_MALLOC;
      }
      for (j = 0; j < m; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(lshell->faces);
          EG_free(lshell);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        lshell->faces[j] = lbody->faces.objs[iref-1];
#endif
/*@+nullderef@*/
      }

      obj = lbody->shells.objs[i];
      obj->oclass = SHELL;
      obj->mtype  = mtype;
      obj->blind  = lshell;
      stat = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
      if (stat != EGADS_SUCCESS) return stat;
    }
  
  /* finish off the body */
  
  if (lbody->shells.nobjs != 0) {
    lbody->senses = (int *) EG_alloc(lbody->shells.nobjs*sizeof(int));
    if (lbody->senses == NULL) return EGADS_MALLOC;
    n = Fread(lbody->senses, sizeof(int), lbody->shells.nobjs, fp);
    if (n != lbody->shells.nobjs) return EGADS_READERR;
  }
  n = Fread(lbody->bbox, sizeof(double), 6, fp);
  if (n != 6) return EGADS_READERR;
  stat = EG_readAttrs(fp, (egAttrs **) &bobj->attrs);
  if (stat != EGADS_SUCCESS) return stat;

  return EGADS_SUCCESS;
}


int
EG_importModel(egObject *context, const size_t nbytes, const char *stream,
               egObject **model)
{
  int       i, n, rev[2];
  liteModel *lmodel;
  egObject  *obj;
  stream_T  myStream;
  stream_T  *fp = &myStream;
  
  *model = NULL;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (context->topObj != NULL)       return EGADS_EXISTS;

  fp->size = nbytes;
  fp->ptr  = 0;
  fp->data = (void *) stream;
  fp->swap = 0;

  /* get header */
  n = Fread(&i,             sizeof(int),    1, fp);
  if (n != 1) {
    return EGADS_READERR;
  }
  if (i != MAGIC) {
    swap(&i, sizeof(int));
    if (i != MAGIC) {
      printf(" EGADS Error: Not a EGADS Lite file!\n");
      return EGADS_READERR;
    }
    fp->swap = 1;
  }
  n = Fread(rev,            sizeof(int),    2, fp);
  if (n != 2) {
    return EGADS_READERR;
  }
  if (rev[0] != 1) {
    printf(" EGADS Error: EGADS Lite file revision = %d %d!\n", rev[0], rev[1]);
    return EGADS_READERR;
  }
  
  lmodel = (liteModel *) EG_alloc(sizeof(liteModel));
  if (lmodel == NULL) {
    printf(" EGADS Error: Malloc of Model!\n");
    return EGADS_MALLOC;
  }
  n = Fread(lmodel->bbox,   sizeof(double), 6, fp);
  if (n != 6) {
    EG_free(lmodel);
    return EGADS_READERR;
  }
  n = Fread(&lmodel->nbody, sizeof(int),    1, fp);
  if (n != 1) {
    EG_free(lmodel);
    return EGADS_READERR;
  }
  lmodel->bodies = NULL;
  if (lmodel->nbody > 0) {
    lmodel->bodies = (egObject **) EG_alloc(lmodel->nbody*sizeof(egObject *));
    if (lmodel->bodies == NULL) {
      printf(" EGADS Error: Malloc of %d Bodies!\n", lmodel->nbody);
      EG_free(lmodel);
      return EGADS_MALLOC;
    }
    for (i = 0; i < lmodel->nbody; i++) lmodel->bodies[i] = NULL;
  }
  i = EG_makeObject(context, &obj);
  if (i != EGADS_SUCCESS) {
    EG_free(lmodel->bodies);
    printf(" EGADS Error: makeObject on Model = %d!\n", i);
    EG_free(lmodel);
    return i;
  }
  obj->oclass = MODEL;
  obj->mtype  = 0;
  obj->blind  = lmodel;
  i = EG_readAttrs(fp, (egAttrs **) &obj->attrs);
  if (i != EGADS_SUCCESS) {
    EG_close(context);
    return i;
  }

  /* get all of the bodies */
  for (n = 0; n < lmodel->nbody; n++) {
    i = EG_readBody(context, obj, n, fp);
    if (i == EGADS_SUCCESS) continue;
    /* errorred out -- cleanup */
    EG_close(context);
    return i;
  }

  *model = context->topObj = obj;

  return EGADS_SUCCESS;
}
