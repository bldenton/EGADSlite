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
#include "liteDevice.h"


/* #define DEBUG    */
/* #define FULLATTR */


  extern int EG_close(egObject *context);
#ifdef __NVCC__
  extern int EG_evaluateDev(const egObject *geom_d, const double *param,
                            double *ev);
  extern int EG_attributeRetDev(const egObject *obj_d, const char *name,
                                int *atype, int *len, int **ints,
                                double **reals, char **str);
  extern int EG_addStrAttrDev(egObject *obj_d, const char *name,
                              const char *str);
#else
  extern int EG_evaluate(const egObject *geom, const double *param, double *ev);
  extern int EG_attributeRet(const egObject *obj, const char *name, int *atype,
                             int *len, /*@null@*/ const int **ints,
                                       /*@null@*/ const double **reals,
                                       /*@null@*/ const char **str);
#endif


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
  
  buf = (char *) buffer;
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
  
  buf = (char *) data;
  memcpy(data, &(((char *) stream->data)[stream->ptr]), size*nitems);
  if ((size != sizeof(char)) && (stream->swap == 1))
    for (i = 0; i < nitems; i++) swap(&buf[i*size], size);
  stream->ptr += size*nitems;

  return nitems;
}


#ifdef FULLATTR
static void
EG_attrBuildSeq(egAttrs *attrs)
{
  int       i, j, l, n, snum, *hit, nospace = 0, nseqs = 0;
  char      *root, *newname;
  egAttr    *attr;
  egAttrSeq *seqs = NULL, *tmp;
  egAttrs   attrs_, *attrs_h = &attrs_;

  EG_GET_ATTRS(attrs_h, attrs);
  
  hit = (int *) EG_alloc(attrs_h->nattrs*sizeof(int));
  if (hit == NULL) {
    printf(" EGADS Internal: Malloc on %d attributes!\n", attrs_h->nattrs);
    return;
  }
  for (i = 0; i < attrs_h->nattrs; i++) hit[i] = 0;
  
  /* build the sequence structure */
  for (i = 0; i < attrs_h->nattrs; i++) {
    if (hit[i] != 0) continue;
    if (attrs->attrs[i].name == NULL) continue;
    l = strlen(attrs->attrs[i].name);
    for (n = j = 0; j < l; j++)
      if (attrs->attrs[i].name[j] == 32) n++;
    if (n == 0) {
      /* only the first can have no space */
      nospace = 1;
    } else if (n >  1) {
      printf(" EGADS Internal: More than a single space (%d) in an Attr name!\n",
             n);
      continue;
    }
    /* make the root name */
    n    = l;
    root = attrs->attrs[i].name;
    if (nospace == 0) {
      root = EG_strdup(attrs->attrs[i].name);
      if (root == NULL) {
        printf(" EGADS Internal: Null root on %s!\n", attrs->attrs[i].name);
        EG_free(hit);
        return;
      }
      for (n = 0; n < l; n++)
        if (root[n] == 32) {
          root[n] = 0;
          break;
        }
    }
    
    /* count the members */
    snum = hit[i] = 1;
    for (j = i+1; j < attrs->nattrs; j++) {
      if (hit[j] != 0) continue;
      if (attrs->attrs[j].name == NULL) continue;
      if (strlen(attrs->attrs[j].name) < n) continue;
      if (attrs->attrs[j].name[n] == 32)
        if (strncmp(attrs->attrs[i].name, attrs->attrs[j].name, n-1) == 0)
          snum++;
    }
    
    /* only a single member */
    if (snum == 1) {
      if (nospace == 1) continue;
      /* remove existing seq number */
      EG_free(attrs->attrs[i].name);
      attrs->attrs[i].name = root;
      continue;
    }
    
    /* build the sequence */
    if (nseqs == 0) {
      seqs = (egAttrSeq *) EG_alloc(sizeof(egAttrSeq));
      if (seqs == NULL) {
        EG_free(root);
        printf(" EGADS Internal: Malloc on Base Sequence!\n");
        continue;
      }
    } else {
      tmp = (egAttrSeq *) EG_reall(seqs, (nseqs+1)*sizeof(egAttrSeq));
      if (tmp == NULL) {
        EG_free(root);
        printf(" EGADS Internal: Malloc on %d Sequence!\n", nseqs+1);
        continue;
      }
      seqs = tmp;
    }
    seqs[nseqs].attrSeq = (int *) EG_alloc(snum*sizeof(int));
    if (seqs[nseqs].attrSeq == NULL) {
      EG_free(root);
      printf(" EGADS Internal: Malloc on %d Attr Seq Pointers!\n", snum);
      continue;
    }
    if (nospace == 1) root = EG_strdup(attrs->attrs[i].name);
    seqs[nseqs].nSeq = snum;
    seqs[nseqs].root = root;

    /* load the sequence */
    seqs[nseqs].attrSeq[0] = i;
    snum = 1;
    for (j = i+1; j < attrs->nattrs; j++) {
      if (hit[j] != 0) continue;
      if (attrs->attrs[j].name == NULL) continue;
      if (strlen(attrs->attrs[j].name) < n) continue;
      if (attrs->attrs[j].name[n] == 32)
        if (strncmp(attrs->attrs[i].name, attrs->attrs[j].name, n-1) == 0) {
          seqs[nseqs].attrSeq[snum] = j;
          snum++;
          hit[j] = 1;
        }
    }
    
    /* check/correct the sequence numbers */
    for (j = 0; j < snum; j++) {
      attr = &attrs->attrs[seqs[nseqs].attrSeq[j]];
      if ((j != 0) || (nospace == 0)) {
        l    = 0;
        sscanf(&attr->name[n], "%d", &l);
#ifdef DEBUG
        if (l == j+1) printf(" seq = %d, oldname = %s\n", j+1, attr->name);
#endif
        if (l == j+1) continue;
      }
      newname = (char *) EG_alloc((n+8)*sizeof(char));
      if (newname == NULL) {
        printf(" EGADS Internal: Malloc on name %s!\n", root);
        continue;
      }
      snprintf(newname, n+8, "%s %d", root, j+1);
      EG_free(attr->name);
      attr->name = newname;
#ifdef DEBUG
      printf(" seq = %d, newname = %s\n", j+1, newname);
#endif
    }
#ifdef DEBUG
    printf(" seq %d: root = %s,  snum = %d\n",
           nseqs, seqs[nseqs].root, seqs[nseqs].nSeq);
    for (j = 0; j < seqs[nseqs].nSeq; j++) {
      attr = &attrs->attrs[seqs[nseqs].attrSeq[j]];
      printf(" %d: %d  %s\n", j+1, seqs[nseqs].attrSeq[j], attr->name);
    }
#endif
    nseqs++;
  }
  EG_free(hit);

  attrs_h->nseqs = nseqs;
  attrs_h->seqs  = seqs;
  EG_SET_ATTRS(attrs, attrs_h);
}
#endif


static int
EG_addStrAttr(egObject *obj, const char *name, const char *str)
{
  int     i, length, find = -1;
  char    *temp;
  egAttr  *attr;
  egAttrs *attrs;
  egObject obj_, *obj_h = &obj_;
  egAttrs  attrs_, *attrs_h = &attrs_;
  egAttr   attr_, *attr_h = &attr_;
  
  if (obj == NULL)               return EGADS_NULLOBJ;
  EG_GET_OBJECT(obj_h, obj);
  if (obj_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj_h->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj_h->oclass == NIL)        return EGADS_EMPTY;
  if (obj_h->oclass == REFERENCE)  return EGADS_REFERCE;
  
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
  attrs = (egAttrs *) obj_h->attrs;
  
  if (attrs != NULL) {
    EG_NEW(&temp, char, length+2);
    if (temp != NULL) {
      EG_GET_ATTRS(attrs_h, attrs);
      for (i = 0; i < attrs_h->nattrs; i++) {
        EG_GET_STRN(temp, &(attrs_h->attrs[i].name), length+1);
        if (strcmp(temp,name) == 0) {
          find = i;
          break;
        }
      }
      EG_FREE(temp);
    }
  }
  
  if ((find != -1) && (attrs != NULL)) {
    
    /* an existing attribute -- reset the values */
    
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));
    if (attr_h->type == ATTRINT) {
      if (attr_h->length != 1)
        EG_FREE(attr_h->vals.integers);
    } else if (attr_h->type == ATTRREAL) {
      if (attr_h->length != 1)
        EG_FREE(attr_h->vals.reals);
    } else if (attr_h->type == ATTRCSYS) {
      EG_FREE(attr_h->vals.reals);
    } else if (attr_h->type == ATTRSTRING) {
      EG_FREE(attr_h->vals.string);
    }
    
  } else {
    
    if (attrs == NULL) {
      EG_NEW(&attrs, egAttrs, 1);
      if (attrs == NULL) {
        printf(" EGADS Internal: Attrs MALLOC for %s (EG_addStrAttr)!\n",
               name);
        return EGADS_MALLOC;
      }
      attrs_h->nattrs = 0;
      attrs_h->attrs  = NULL;
      attrs_h->nseqs  = 0;
      attrs_h->seqs   = NULL;
/*@-nullret@*/
      EG_SET_ATTRS(attrs, attrs_h);
/*@+nullret@*/
      EG_SET_OBJECT_PTR(&(obj->attrs), &attrs);
    }
    EG_GET_ATTRS(attrs_h, attrs);
    if (attrs_h->attrs == NULL) {
      EG_NEW(&attr, egAttr, (attrs_h->nattrs+1));
    } else {
      EG_REALLOC(&attr, attrs_h->attrs, egAttr, attrs_h->nattrs,
                 (attrs_h->nattrs+1));
    }
    if (attr == NULL) {
      printf(" EGADS Internal: Attr MALLOC for %s (EG_addStrAttr)!\n",
             name);
      return EGADS_MALLOC;
    }
    attrs_h->attrs = attr;
    find = attrs_h->nattrs;
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));
    attr_h->vals.string = NULL;
    EG_SET_STR(&(attr_h->name), name);
    if (attr_h->name == NULL) return EGADS_MALLOC;
    EG_SET_ATTR(&(attrs_h->attrs[find]), attr_h);
    attrs_h->nattrs += 1;
  }
  
  EG_GET_ATTR(attr_h, &(attrs_h->attrs[find]));

  attr_h->type        = ATTRSTRING;
  attr_h->length      = 0;
  EG_SET_STR(&(attr_h->vals.string), str);
  if (attr_h->vals.string != NULL) {
    attr_h->length = strlen(str);
  }
  EG_SET_ATTR(&(attrs_h->attrs[find]), attr_h);

  EG_SET_ATTRS(attrs, attrs_h);

  return EGADS_SUCCESS;
}


static void
EG_freeAttrs(egAttrs **attrx)
{
  int     i;
  egAttrs *attrs;
  egAttrs attrs_, *attrs_h = &attrs_;
  egAttr  attr_, *attr_h = &attr_;
  egAttrs *nil = NULL;
  
  attrs = *attrx;
  if (attrs == NULL) return;
  EG_GET_ATTRS(attrs_h, *attrx);
  
  EG_SET_ATTRS_PTR(attrx, nil);
  /* remove any attributes */
  for (i = 0; i < attrs_h->nattrs; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
    if (attr_.name != NULL) EG_FREE(attr_.name);
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


static int
EG_readString(stream_T *fp, char **string)
{
  int    len;
  size_t n;
  char   *string_h;
  int    status = EGADS_SUCCESS;
  
  *string = string_h = NULL;
  n = Fread(&len, sizeof(int), 1, fp);
  if (n   != 1) return EGADS_READERR;
  if (len <  0) return EGADS_READERR;
  if (len == 0) return EGADS_SUCCESS;
  
  string_h = (char *) EG_alloc(len*sizeof(char));
  if (string_h == NULL) return EGADS_MALLOC;
  
  n = Fread(string_h, sizeof(char), len, fp);
  if (n != len) {
    EG_free(string_h);
    *string = NULL;
    return EGADS_READERR;
  }
  EG_SET_STR(&(string[0]), string_h);

  EG_free(string_h);
  
  return status;
}


static int
EG_readAttrs(stream_T *fp, egAttrs **attrx)
{
  int     nattr, i, status;
#ifdef FULLATTR
  int     len, j, nspace = 0;
#endif
  size_t  n;
  egAttr  *attr;
  egAttrs *attrs;
  egAttrs attrs_, *attrs_h = &attrs_;
  egAttr  attr_, *attr_h = &attr_;
  void    *temp;
  
  *attrx = NULL;
  n = Fread(&nattr, sizeof(int), 1, fp);
  if (n     != 1) return EGADS_READERR;
  if (nattr == 0) return EGADS_SUCCESS;
  
  EG_NEW(&attrs, egAttrs, 1);
  if (attrs == NULL) return EGADS_MALLOC;
  EG_NEW(&attr, egAttr, nattr);
  if (attr == NULL) {
    EG_FREE(attrs);
    return EGADS_MALLOC;
  }
  attrs_h->nattrs = nattr;
  attrs_h->attrs  = attr;
  EG_SET_ATTRS(attrs, attrs_h);
  for (i = 0; i < nattr; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));
    attr_.name   = NULL;
    attr_.length = 1;
    attr_.type   = ATTRINT;
    EG_SET_ATTR(&(attrs_h->attrs[i]), attr_h);
  }
  
/*@-mustfreefresh@*/
  for (i = 0; i < nattr; i++) {
    EG_GET_ATTR(attr_h, &(attrs_h->attrs[i]));

    n = Fread(&attr_.type,   sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    n = Fread(&attr_.length, sizeof(int), 1, fp);
    if (n != 1) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
    status = EG_readString(fp, &attr_.name);
    if (status != EGADS_SUCCESS) {
      EG_freeAttrs(&attrs);
      return EGADS_READERR;
    }
#ifdef FULLATTR
    if (attr_.name != NULL) {
      len = strlen(attr_.name);
      for (j = 0; j < len; j++)
        if (attr_.name[j] == 32) nspace++;
    }
#endif
    if (attr_.type == ATTRINT) {
      n = attr_.length;
      if (attr_.length == 1) {
        n = Fread(&attr_.vals.integer, sizeof(int), 1, fp);
      } else if (attr_.length > 1) {
        EG_NEW(&attr_.vals.integers, int, attr_.length);
        if (attr_.vals.integers == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        temp = EG_alloc(attr_.length*sizeof(int));
        if (temp == NULL) {
          EG_FREE(attr_.vals.integers);
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread((int *) temp, sizeof(int), attr_.length, fp);
        EG_COPY(attr_.vals.integers, temp, int, attr_.length);
        EG_free(temp);
      }
      if (n != attr_.length) {
        EG_FREE(attr_.vals.integers);
        EG_freeAttrs(&attrs);
        return EGADS_READERR;
      }
    } else if ((attr_.type == ATTRREAL) || (attr_.type == ATTRCSYS)) {
      n = attr_.length;
      if (attr_.length == 1) {
        n = Fread(&attr_.vals.real, sizeof(double), 1, fp);
      } else if (attr_.length > 1) {
        EG_NEW(&attr_.vals.reals, double, attr_.length);
        if (attr_.vals.reals == NULL) {
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        temp = EG_alloc(attr_.length*sizeof(double));
        if (temp == NULL) {
          EG_FREE(attr_.vals.reals);
          EG_freeAttrs(&attrs);
          return EGADS_MALLOC;
        }
        n = Fread((double *) temp, sizeof(double), attr_.length, fp);
        EG_COPY(attr_.vals.reals, temp, double, attr_.length);
        EG_free(temp);
      }
      if (n != attr_.length) {
        EG_FREE(attr_.vals.reals);
        EG_freeAttrs(&attrs);
        return EGADS_READERR;
      }
    } else {
      status = EG_readString(fp, &attr_.vals.string);
      if (status != EGADS_SUCCESS) return EGADS_READERR;
    }
    EG_SET_ATTR(&(attrs_h->attrs[i]), attr_h);
  }
/*@+mustfreefresh@*/

#ifdef FULLATTR
  /* sequences exist! */
  if (nspace != 0) EG_attrBuildSeq(attrs);
#endif
  
  *attrx = attrs;
  return EGADS_SUCCESS;
}


static int
EG_readGeometry(liteGeometry *lgeom, int *iref, stream_T *fp)
{
  int          n, nhead, ndata;
  liteGeometry lgeom_, *lgeom_h = &lgeom_;
  void         *temp;
  
  EG_GET_GEOM(lgeom_h, lgeom);

  *iref           = 0;
  lgeom_h->ref    = NULL;
  lgeom_h->header = NULL;
  lgeom_h->data   = NULL;
/*@-nullret@*/
  EG_SET_GEOM(lgeom, lgeom_h);
/*@+nullret@*/
  n = Fread(iref,   sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&nhead, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(&ndata, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  
  if (nhead != 0) {
    EG_NEW(&(lgeom_h->header), int, nhead);
    EG_COPY(&(lgeom->header), &(lgeom_h->header), int *, 1);
    if (lgeom_h->header == NULL) return EGADS_MALLOC;
    temp = EG_alloc(nhead*sizeof(int));
    if (temp == NULL) return EGADS_MALLOC;
    n = Fread((int *) temp, sizeof(int), nhead, fp);
    EG_COPY(lgeom_h->header, temp, int, nhead);
    EG_free(temp);
    if (n != nhead) return EGADS_READERR;
  }
  EG_NEW(&(lgeom_h->data), double, ndata);
  EG_COPY(&(lgeom->data), &(lgeom_h->data), double *, 1);
  if (lgeom_h->data == NULL) return EGADS_MALLOC;
  temp = EG_alloc(ndata*sizeof(double));
  if (temp == NULL) return EGADS_MALLOC;
  n = Fread((double *) temp, sizeof(double), ndata, fp);
  EG_COPY(lgeom_h->data, temp, double, ndata);
  EG_free(temp);
  if (n != ndata) return EGADS_READERR;
  
  return EGADS_SUCCESS;
}


static int
EG_readBody(egObject *context, egObject *mobject, int bindex, stream_T *fp)
{
  int          i, j, m, n, stat, mtype, iref, ntypes[8];
  double       t, d, x0[2], x1[2], data[6];
  int          atype, alen;
  const int    *ints;
  const double *reals;
  const char   *str;
  egObject     *obj, *bobj, *pcobj;
  liteGeometry *lgeom;
  liteNode     *lnode;
  liteEdge     *ledge;
  liteLoop     *lloop;
  liteFace     *lface;
  liteShell    *lshell;
  liteBody     *lbody;
  liteModel    *lmodel;
  egObject     context_, *context_h = &context_;
  egObject     mobject_, *mobject_h = &mobject_;
  egObject     bobj_, *bobj_h = &bobj_;
  liteModel    lmodel_, *lmodel_h = &lmodel_;
  liteBody     lbody_, *lbody_h = &lbody_;
  void         *nil = NULL;

  if (context == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(context_h, context);
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (mobject == NULL)                 return EGADS_NULLOBJ;
  EG_GET_OBJECT(mobject_h, mobject);
  if (mobject_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (mobject_h->oclass != MODEL)      return EGADS_NOTMODEL;
  lmodel = (liteModel *) mobject_h->blind;
  EG_GET_MODEL(lmodel_h, lmodel);

  n = Fread(&mtype, sizeof(int), 1, fp);
  if (n != 1) return EGADS_READERR;
  n = Fread(ntypes, sizeof(int), 8, fp);
  if (n != 8) return EGADS_READERR;
  EG_NEW(&lbody, liteBody, 1);
  if (lbody == NULL) return EGADS_MALLOC;
  EG_GET_BODY(lbody_h, lbody);
  lbody_h->pcurves.objs   = NULL;
  lbody_h->curves.objs    = NULL;
  lbody_h->surfaces.objs  = NULL;
  lbody_h->nodes.objs     = NULL;
  lbody_h->edges.objs     = NULL;
  lbody_h->loops.objs     = NULL;
  lbody_h->faces.objs     = NULL;
  lbody_h->shells.objs    = NULL;
  lbody_h->senses         = NULL;
  lbody_h->pcurves.nobjs  = 0;
  lbody_h->curves.nobjs   = 0;
  lbody_h->surfaces.nobjs = 0;
  lbody_h->nodes.nobjs    = 0;
  lbody_h->edges.nobjs    = 0;
  lbody_h->loops.nobjs    = 0;
  lbody_h->faces.nobjs    = 0;
  lbody_h->shells.nobjs   = 0;

#ifdef DEBUG
  printf(" Reading Body #%d: %d %d %d %d %d %d %d %d\n", bindex+1, ntypes[0],
         ntypes[1], ntypes[2], ntypes[3], ntypes[4], ntypes[5], ntypes[6],
         ntypes[7]);
#endif

  stat = EG_makeObject(context, &bobj);
  if (stat != EGADS_SUCCESS) {
    EG_FREE(lbody);
    return stat;
  }
  EG_GET_OBJECT(bobj_h, bobj);
  bobj_h->oclass = BODY;
  bobj_h->mtype  = mtype;
  bobj_h->blind  = lbody;
  EG_SET_OBJECT(&bobj, bobj_h);
  EG_SET_OBJECT_PTR(&(lmodel_h->bodies[bindex]), &bobj);
  
  /* make all of the objects */
  if (ntypes[0] > 0) {
    EG_NEW(&(lbody_h->pcurves.objs), egObject *, ntypes[0]);
    if (lbody_h->pcurves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[0]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->pcurves.objs[i]), &nil);
    lbody_h->pcurves.nobjs = ntypes[0];
    for (i = 0; i < ntypes[0]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->pcurves.objs[i]), &obj);
    }
  }
  if (ntypes[1] > 0) {
    EG_NEW(&(lbody_h->curves.objs), egObject *, ntypes[1]);
    if (lbody_h->curves.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[1]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->curves.objs[i]), &nil);
    lbody_h->curves.nobjs = ntypes[1];
    for (i = 0; i < ntypes[1]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->curves.objs[i]), &obj);
    }
  }
  if (ntypes[2] > 0) {
    EG_NEW(&(lbody_h->surfaces.objs), egObject *, ntypes[2]);
    if (lbody_h->surfaces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[2]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->surfaces.objs[i]), &nil);
    lbody_h->surfaces.nobjs = ntypes[2];
    for (i = 0; i < ntypes[2]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->surfaces.objs[i]), &obj);
    }
  }
  if (ntypes[3] > 0) {
    EG_NEW(&(lbody_h->nodes.objs), egObject *, ntypes[3]);
    if (lbody_h->nodes.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[3]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->nodes.objs[i]), &nil);
    lbody_h->nodes.nobjs = ntypes[3];
    for (i = 0; i < ntypes[3]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->nodes.objs[i]), &obj);
    }
  }
  if (ntypes[4] > 0) {
    EG_NEW(&(lbody_h->edges.objs), egObject *, ntypes[4]);
    if (lbody_h->edges.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[4]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->edges.objs[i]), &nil);
    lbody_h->edges.nobjs = ntypes[4];
    for (i = 0; i < ntypes[4]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->edges.objs[i]), &obj);
    }
  }
  if (ntypes[5] > 0) {
    EG_NEW(&(lbody_h->loops.objs), egObject *, ntypes[5]);
    if (lbody_h->loops.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[5]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->loops.objs[i]), &nil);
    lbody_h->loops.nobjs = ntypes[5];
    for (i = 0; i < ntypes[5]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->loops.objs[i]), &obj);
    }
  }
  if (ntypes[6] > 0) {
    EG_NEW(&(lbody_h->faces.objs), egObject *, ntypes[6]);
    if (lbody_h->faces.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[6]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->faces.objs[i]), &nil);
    lbody_h->faces.nobjs = ntypes[6];
    for (i = 0; i < ntypes[6]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->faces.objs[i]), &obj);
    }
  }
  if (ntypes[7] > 0) {
    EG_NEW(&(lbody_h->shells.objs), egObject *, ntypes[7]);
    if (lbody_h->shells.objs == NULL) return EGADS_SUCCESS;
    for (i = 0; i < ntypes[7]; i++)
      EG_SET_OBJECT_PTR(&(lbody_h->shells.objs[i]), &nil);
    lbody_h->shells.nobjs = ntypes[7];
    for (i = 0; i < ntypes[7]; i++) {
      stat = EG_makeObject(context, &obj);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT_PTR(&(lbody_h->shells.objs[i]), &obj);
    }
  }
/*@-nullret@*/
  EG_SET_BODY(lbody, lbody_h);
/*@+nullret@*/

#ifdef DEBUG
  printf(" Reading %d PCurves...\n", ntypes[0]);
#endif
  /* pcurves */
  if (lbody_h->pcurves.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    int header_h[4];
    double data_h[6];
    for (i = 0; i < lbody_h->pcurves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
      EG_GET_GEOM(lgeom_h, lgeom);
      if (iref != 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->pcurves.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->pcurves.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = PCURVE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
      if (mtype != BSPLINE) continue;
      EG_GET_HEADERS_AT(header_h, lgeom_h->header, 4);
      for (n = 2; n < header_h[2]; n++) {
        m     = header_h[3] + 2*n - 4;
        EG_GET_DATA_AT(data_h, &(lgeom_h->data[m]), 6);
        x0[0] = data_h[  2] - data_h[  0];
        x0[1] = data_h[  3] - data_h[  1];
        d     = sqrt(x0[0]*x0[0] + x0[1]*x0[1]);
        if (d != 0.0) {
          x0[0] /= d;
          x0[1] /= d;
        }
        x1[0] = data_h[  4] - data_h[  2];
        x1[1] = data_h[  5] - data_h[  3];
        d     = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
        if (d != 0.0) {
          x1[0] /= d;
          x1[1] /= d;
        }
        d = x0[0]*x1[0] + x0[1]*x1[1];
        if (d < -0.95) {
#ifdef DEBUG
          double last_data_h;
          EG_GET_DATA_AT(&last_data_h, &(lgeom_h->data[header_h[3]-1]), 1);
          printf(" EGADS Info: PCurve %d dot flip at %d/%d (%lf) -- %lf %lf!\n",
                 i, n-2, header_h[2]-2, d,
                 data_h[0], last_data_h);
#endif
          stat = EG_addStrAttr(obj, ".Bad", "CPrev");
          if (stat != EGADS_SUCCESS)
            printf("             EG_addStrAttr CPrev= %d\n", stat);
        }
      }
    }
  }

  /* curves */
#ifdef DEBUG
  printf(" Reading %d Curves...\n", ntypes[1]);
#endif
  if (lbody_h->curves.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->curves.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
      if (iref != 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->curves.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->curves.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = CURVE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }
  
  /* surfaces */
#ifdef DEBUG
  printf(" Reading %d Surfaces...\n", ntypes[2]);
#endif
  if (lbody_h->surfaces.objs != NULL) {
    liteGeometry lgeom_, *lgeom_h = &lgeom_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->surfaces.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lgeom, liteGeometry, 1);
      if (lgeom == NULL) return EGADS_MALLOC;
      stat = EG_readGeometry(lgeom, &iref, fp);
      if (stat != EGADS_SUCCESS) {
        EG_GET_GEOM(lgeom_h, lgeom);
        if (lgeom_h->header != NULL) EG_FREE(lgeom_h->header);
        if (lgeom_h->data   != NULL) EG_FREE(lgeom_h->data);
        EG_FREE(lgeom);
        return stat;
      }
/*@-nullderef@*/
      if (iref < 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->curves.objs[-iref-1]));
      }
/*@+nullderef@*/
      if (iref > 0) {
        EG_SET_OBJECT_PTR(&(lgeom->ref), &(lbody_h->surfaces.objs[iref-1]));
      }
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->surfaces.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = SURFACE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lgeom;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }
  
  /* nodes */
#ifdef DEBUG
  printf(" Reading %d Nodes...\n", ntypes[3]);
#endif
  if (lbody_h->nodes.objs != NULL) {
    liteNode lnode_, *lnode_h = &lnode_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->nodes.nobjs; i++) {
      EG_NEW(&lnode, liteNode, 1);
      if (lnode == NULL) return EGADS_MALLOC;
      n = Fread(lnode_h->xyz,  sizeof(double), 3, fp);
      if (n != 3) {
        EG_FREE(lnode);
        return EGADS_READERR;
      }
      n = Fread(&lnode_h->tol, sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(lnode);
        return EGADS_READERR;
      }
      EG_SET_NODE(lnode, lnode_h);
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->nodes.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = NODE;
      obj_h->mtype  = 0;
      obj_h->blind  = lnode;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }
  
  /* edges */
#ifdef DEBUG
  printf(" Reading %d Edges...\n", ntypes[4]);
#endif
  if (lbody_h->edges.objs != NULL) {
    liteEdge ledge_, *ledge_h = &ledge_;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->edges.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&ledge, liteEdge, 1);
      if (ledge == NULL) return EGADS_MALLOC;
      EG_GET_EDGE(ledge_h, ledge);
      ledge_h->curve    = NULL;
      ledge_h->nodes[0] = NULL;
      ledge_h->nodes[1] = NULL;
  
      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->curve),
                                       &(lbody_h->curves.objs[iref-1]));
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->nodes[0]),
                                       &(lbody_h->nodes.objs[iref-1]));
/*@+nullderef@*/
      n = Fread(&iref, sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullderef@*/
      if (iref != 0) EG_GET_OBJECT_PTR(&(ledge_h->nodes[1]),
                                       &(lbody_h->nodes.objs[iref-1]));
/*@+nullderef@*/
     
      n = Fread(ledge_h->trange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
      n = Fread(ledge_h->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
      n = Fread(&ledge_h->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(ledge);
        return EGADS_READERR;
      }
/*@-nullret@*/
      EG_SET_EDGE(ledge, ledge_h);
/*@+nullret@*/

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->edges.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = EDGE;
      obj_h->mtype  = mtype;
      obj_h->blind  = ledge;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }

  /* loops */
#ifdef DEBUG
  printf(" Reading %d Loops...\n", ntypes[5]);
#endif
  if (lbody_h->loops.objs != NULL) {
    liteLoop lloop_, *lloop_h = &lloop_;
    int      *itemp;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->loops.nobjs; i++) {
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lloop, liteLoop, 1);
      if (lloop == NULL) return EGADS_MALLOC;
      lloop_h->nedges  = m;
      lloop_h->surface = NULL;
      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(lloop);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
        EG_GET_OBJECT_PTR(&(lloop_h->surface),
                          &(lbody_h->surfaces.objs[iref-1]));
/*@+nullderef@*/
        m *= 2;
      }
      n = Fread(lloop_h->bbox,  sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lloop);
        return EGADS_READERR;
      }
      if (lloop_h->nedges != 0) {
        egObject **otemp;
        EG_NEW(&(lloop_h->senses), int, lloop_h->nedges);
        if (lloop_h->senses == NULL) {
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }

        itemp = (int *) EG_alloc(lloop_h->nedges*sizeof(int));
        if (itemp == NULL) {
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        n = Fread((int *) itemp,  sizeof(int), lloop_h->nedges, fp);
        if (n != lloop_h->nedges) {
          EG_free(itemp);
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_READERR;
        }
        EG_COPY(lloop_h->senses, itemp, int, lloop_h->nedges);
        EG_free(itemp);

        EG_NEW(&(lloop_h->edges), egObject *, m);
        if (lloop_h->edges == NULL) {
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        otemp = (egObject **) EG_alloc(m*sizeof(egObject *));
        if (otemp == NULL) {
          EG_FREE(lloop_h->edges);
          EG_FREE(lloop_h->senses);
          EG_FREE(lloop);
          return EGADS_MALLOC;
        }
        for (j = 0; j < m; j++) {
          n = Fread(&iref, sizeof(int), 1, fp);
          if (n != 1) {
            EG_free(otemp);
            EG_FREE(lloop_h->edges);
            EG_FREE(lloop_h->senses);
            EG_FREE(lloop);
            return EGADS_READERR;
          }
          if (j < lloop_h->nedges) {
/*@-nullderef@*/
            EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->edges.objs[iref-1]));
/*@+nullderef@*/
          } else {
/*@-nullderef@*/
            EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->pcurves.objs[iref-1]));
/*@+nullderef@*/
          }
        }
        EG_COPY(lloop_h->edges, otemp, egObject *, m);
        if (lloop_h->surface != NULL) {
          liteEdge ledge_, *ledge_h = &ledge_;
          for (n = 0; n < lloop_h->nedges; n++) {
            EG_GET_OBJECT_PTR(&ledge, &(otemp[n]->blind));
            EG_GET_EDGE(ledge_h, ledge);
            pcobj = otemp[lloop_h->nedges+n];
#ifdef __NVCC__
            stat  = EG_attributeRetDev(pcobj, ".Bad", &atype, &alen,
                                       NULL, NULL, NULL);
            if (stat != EGADS_SUCCESS) continue;
            EG_evaluateDev(pcobj, &ledge_h->trange[0], data);
#else
            stat  = EG_attributeRet(pcobj, ".Bad", &atype, &alen,
                                    &ints, &reals, &str);
            if (stat != EGADS_SUCCESS) continue;
            EG_evaluate(pcobj, &ledge_h->trange[0], data);
#endif
            d     = sqrt(data[2]*data[2] + data[3]*data[3]);
            x0[0] = x0[1] = 0.0;
            if (d != 0.0) {
              x0[0] = data[2]/d;
              x0[1] = data[3]/d;
            }
            for (j = 1; j < 1000; j++) {
              t = ledge_h->trange[0]+j*(ledge_h->trange[1]-ledge_h->trange[0])/999.;
#ifdef __NVCC__
              EG_evaluateDev(pcobj, &t, data);
#else
              EG_evaluate(pcobj, &t, data);
#endif
              d = sqrt(data[2]*data[2] + data[3]*data[3]);
              x1[0] = x1[1] = 0.0;
              if (d != 0.0) {
                x1[0] = data[2]/d;
                x1[1] = data[3]/d;
              }
              if (x0[0]*x1[0] + x0[1]*x1[1] < -0.95) {
#ifdef __NVCC__
                stat = EG_addStrAttrDev(pcobj, ".Bad", "fold");
#else
                stat = EG_addStrAttr(pcobj, ".Bad", "fold");
#endif
                if (stat != EGADS_SUCCESS)
                  printf(" EGADS Info: EG_addStrAttr fold= %d\n", stat);
              }
              x0[0] = x1[0];
              x0[1] = x1[1];
            }
          }
        }
        EG_free(otemp);
      }
/*@-nullret@*/
      EG_SET_LOOP(lloop, lloop_h);
/*@+nullret@*/
      
      EG_GET_OBJECT_PTR(&obj, &(lbody_h->loops.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = LOOP;
      obj_h->mtype  = mtype;
      obj_h->blind  = lloop;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }
  
  /* faces */
#ifdef DEBUG
  printf(" Reading %d Faces...\n", ntypes[6]);
#endif
  if (lbody_h->faces.objs != NULL) {
    liteFace lface_, *lface_h = &lface_;
#ifdef DEBUG
    liteLoop lloop_, *lloop_h = &lloop_;
    liteLoop *lloop_d;
#endif
    int      *itemp;
    egObject obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->faces.nobjs; i++) {
      egObject **otemp;
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lface, liteFace, 1);
      if (lface == NULL) return EGADS_MALLOC;
      lface_h->nloops  = m;
      lface_h->surface = NULL;
      n = Fread(&iref,  sizeof(int), 1, fp);
      if (n != 1) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      if (iref != 0) {
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(lface_h->surface),
                          &(lbody_h->surfaces.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
      n = Fread(lface_h->urange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(lface_h->vrange, sizeof(double), 2, fp);
      if (n != 2) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(lface_h->bbox,   sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      n = Fread(&lface_h->tol,   sizeof(double), 1, fp);
      if (n != 1) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      EG_NEW(&(lface_h->senses), int, lface_h->nloops);
      if (lface_h->senses == NULL) {
        EG_FREE(lface);
        return EGADS_READERR;
      }
      itemp = (int *) EG_alloc(lface_h->nloops*sizeof(int));
      if (itemp == NULL) {
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_MALLOC;
      }
      n = Fread((int *) itemp,  sizeof(int), lface_h->nloops, fp);
      if (n != lface_h->nloops) {
        EG_free(itemp);
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_READERR;
      }
      EG_COPY(lface_h->senses, itemp, int, lface_h->nloops);
      EG_free(itemp);
      EG_NEW(&(lface_h->loops), egObject *, lface_h->nloops);
      if (lface_h->loops == NULL) {
        EG_FREE(lface_h->senses);
        EG_FREE(lface);
        return EGADS_MALLOC;
      }
      otemp = (egObject **) EG_alloc(lface_h->nloops*sizeof(egObject *));
      for (j = 0; j < lface_h->nloops; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(otemp);
          EG_FREE(lface_h->loops);
          EG_FREE(lface_h->senses);
          EG_FREE(lface);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->loops.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
/*@-nullpass@*/
      EG_COPY(lface_h->loops, otemp, egObject *, lface_h->nloops);
/*@+nullpass@*/
      EG_free(otemp);

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->faces.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = FACE;
      obj_h->mtype  = mtype;
      obj_h->blind  = lface;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
#ifdef DEBUG
      for (j = 0; j < lface_h->nloops; j++) {
        EG_GET_OBJECT_PTR(&(obj), &(lface_h->loops[j]));
        EG_GET_OBJECT_PTR(&(obj), &(obj->blind));
        lloop_d = (liteLoop *) obj;
        EG_GET_LOOP(lloop_h, lloop_d);
        if (lloop_h->surface == NULL) continue;
        for (n = 0; n < lloop_h->nedges; n++) {
          EG_GET_OBJECT_PTR(&(obj), &(lloop_h->edges[lloop_h->nedges+n]));
#ifdef __NVCC__
          stat  = EG_attributeRetDev(obj, ".Bad", &atype, &alen,
                                     NULL, NULL, &str);
#else
          stat = EG_attributeRet(obj, ".Bad", &atype,
                                 &alen, &ints, &reals, &str);
#endif
          if ((stat == EGADS_SUCCESS) && (atype == ATTRSTRING)) {
            if (strcmp(str, "fold") == 0)
              printf(" EGADS Info: Body %d Face %d Loop#%d/Edge#%d Bad PCurve -- %s!\n",
                     bindex+1, i+1, j+1, n+1, str);
#ifdef __NVCC__
            EG_free(str);
#endif
          }
        }
      }
#endif
/*@-nullret@*/
      EG_SET_FACE(lface, lface_h);
/*@+nullret@*/
    }
  }

  /* shells */
#ifdef DEBUG
  printf(" Reading %d Shells...\n", ntypes[7]);
#endif
  if (lbody_h->shells.objs != NULL) {
    liteShell lshell_, *lshell_h = &lshell_;
    egObject  obj_, *obj_h = &obj_;
    for (i = 0; i < lbody_h->shells.nobjs; i++) {
      egObject  **otemp;
      n = Fread(&mtype, sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      n = Fread(&m,     sizeof(int), 1, fp);
      if (n != 1) return EGADS_READERR;
      EG_NEW(&lshell, liteShell, 1);
      if (lshell == NULL) return EGADS_MALLOC;
      lshell_h->nfaces = m;
      n = Fread(lshell_h->bbox, sizeof(double), 6, fp);
      if (n != 6) {
        EG_FREE(lshell);
        return EGADS_READERR;
      }
      EG_NEW(&(lshell_h->faces), egObject *, m);
      if (lshell_h->faces == NULL) {
        EG_FREE(lshell);
        return EGADS_MALLOC;
      }
      otemp = (egObject **) EG_alloc(m*sizeof(egObject *));
      if (otemp == NULL) {
        EG_FREE(lshell_h->faces);
        EG_FREE(lshell);
        return EGADS_MALLOC;
      }
      for (j = 0; j < m; j++) {
        n = Fread(&iref, sizeof(int), 1, fp);
        if (n != 1) {
          EG_free(otemp);
          EG_FREE(lshell_h->faces);
          EG_FREE(lshell);
          return EGADS_READERR;
        }
/*@-nullderef@*/
#ifndef __clang_analyzer__
        EG_GET_OBJECT_PTR(&(otemp[j]), &(lbody_h->faces.objs[iref-1]));
#endif
/*@+nullderef@*/
      }
      EG_COPY(lshell_h->faces, otemp, egObject *, m);
      EG_free(otemp);
      EG_SET_SHELL(lshell, lshell_h);

      EG_GET_OBJECT_PTR(&obj, &(lbody_h->shells.objs[i]));
      EG_GET_OBJECT(obj_h, obj);
      obj_h->oclass = SHELL;
      obj_h->mtype  = mtype;
      obj_h->blind  = lshell;
      stat = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
      if (stat != EGADS_SUCCESS) return stat;
      EG_SET_OBJECT(&obj, obj_h);
    }
  }
  
  /* finish off the body */
  
  if (lbody_h->shells.nobjs != 0) {
    int *itemp;
    EG_NEW(&(lbody_h->senses), int, lbody_h->shells.nobjs);
    if (lbody_h->senses == NULL) return EGADS_MALLOC;
    itemp = (int *) EG_alloc(lbody_h->shells.nobjs*sizeof(int));
    if (itemp == NULL) {
      EG_FREE(lbody_h->senses);
      return EGADS_MALLOC;
    }
    n = Fread(itemp, sizeof(int), lbody_h->shells.nobjs, fp);
    if (n != lbody_h->shells.nobjs) {
      EG_free(itemp);
      EG_FREE(lbody_h->senses);
      return EGADS_READERR;
    }
    EG_COPY(lbody_h->senses, itemp, int, lbody_h->shells.nobjs);
    EG_free(itemp);
  }
  n = Fread(lbody_h->bbox, sizeof(double), 6, fp);
  if (n != 6) return EGADS_READERR;
  EG_GET_OBJECT(bobj_h, bobj);
  stat = EG_readAttrs(fp, (egAttrs **) &bobj_h->attrs);
  if (stat != EGADS_SUCCESS) return stat;
  EG_SET_OBJECT(&bobj, bobj_h);
/*@-nullret@*/
  EG_SET_BODY(lbody, lbody_h);
/*@+nullret@*/

  return EGADS_SUCCESS;
}


static int
EG_freeLiteModel(/*@null@*/ liteModel *model)
{
  egObject *bodies_d = NULL;

  if (NULL != model) {
    EG_COPY(&(bodies_d), &(model->bodies), egObject *, 1);
    if (NULL != bodies_d) EG_FREE(bodies_d);
    EG_FREE(model);
  }
  return EGADS_SUCCESS;
}


static int
EG_allocLiteModel(int nbody, double bbox[6], liteModel **model)
{
  int       i;
  liteModel obj_, *obj_h = &obj_;
  liteModel *obj = NULL;
  void      *nil = NULL;

  obj_h->nbody = nbody;
  obj_h->bbox[0] = bbox[0]; obj_h->bbox[1] = bbox[1]; obj_h->bbox[2] = bbox[2];
  obj_h->bbox[3] = bbox[3]; obj_h->bbox[4] = bbox[4]; obj_h->bbox[5] = bbox[5];
  obj_h->bodies = NULL;
  if (nbody > 0) {
    EG_NEW((void **)&(obj_h->bodies), egObject *, obj_h->nbody);
    if (obj_h->bodies == NULL) goto modelCleanup;
  }

  EG_NEW((void **)&(obj), liteModel, 1);
  if (obj == NULL) goto modelCleanup;
  EG_COPY(obj, obj_h, liteModel, 1);
  if (obj == NULL) goto modelCleanup;

  for (i = 0; i < obj_h->nbody; ++i) {
    egObject *ptr;
/*@-nullderef@*/
    EG_COPY(&(obj_h->bodies[i]), &(nil), void *, 1);
    EG_GET_OBJECT_PTR(&ptr, &(obj_h->bodies[i]));
/*@+nullderef@*/
    if (ptr != NULL) goto modelCleanup;
  }

  *model = obj;
  return EGADS_SUCCESS;

modelCleanup:
/*@-dependenttrans@*/
  EG_FREE(obj);
/*@+dependenttrans@*/
  EG_FREE(obj_h->bodies);
  return EGADS_MALLOC;
}


int
EG_importModel(egObject *context, const size_t nbytes, const char *stream,
               egObject **model)
{
  int       i, n, rev[2];
  liteModel *lmodel = NULL;
  egObject  obj_, *obj_h = &obj_;
  egObject  *obj;
  egObject  context_;
  egObject  *context_h = &context_;
  stream_T  myStream;
  stream_T  *fp = &myStream;
  int       nbody;
  double    bbox[6];
  
  *model = NULL;

  EG_GET_OBJECT(context_h, context);

  if (context_h == NULL)               return EGADS_NULLOBJ;
  if (context_h->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context_h->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (context_h->topObj != NULL)       return EGADS_EXISTS;

  fp->size = nbytes;
  fp->ptr  = 0;
  fp->data = (void *) stream;
  fp->swap = 0;

  /* get header */
  n = Fread(&i,     sizeof(int),    1, fp);
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
  n = Fread(rev,    sizeof(int),    2, fp);
  if (n != 2) {
    return EGADS_READERR;
  }
  if (rev[0] != 1) {
    printf(" EGADS Error: EGADS Lite file revision = %d %d!\n", rev[0], rev[1]);
    return EGADS_READERR;
  }

  n = Fread(bbox,   sizeof(double), 6, fp);
  if (n != 6) {
    return EGADS_READERR;
  }
  n = Fread(&nbody, sizeof(int),    1, fp);
  if (n != 1) {
    return EGADS_READERR;
  }
  EG_allocLiteModel(nbody, bbox, &lmodel);

  i = EG_makeObject(context, &obj);
  if (i != EGADS_SUCCESS) {
    EG_freeLiteModel(lmodel);
    printf(" EGADS Error: makeObject on Model = %d!\n", i);
    return i;
  }
  EG_GET_OBJECT(obj_h, obj);
  obj_h->oclass = MODEL;
  obj_h->mtype  = 0;
  obj_h->blind  = lmodel;
  i = EG_readAttrs(fp, (egAttrs **) &obj_h->attrs);
  if (i != EGADS_SUCCESS) {
    EG_close(context);
    return i;
  }
/*@-nullret@*/
  EG_SET_OBJECT(&obj, obj_h);
/*@+nullret@*/

  /* get all of the bodies */
  for (n = 0; n < nbody; n++) {
    i = EG_readBody(context, obj, n, fp);
    if (i == EGADS_SUCCESS) continue;
    /* errorred out -- cleanup */
    EG_close(context);
    return i;
  }

  EG_SET_OBJECT_PTR(&(context->topObj), &obj);
  *model = obj;

  return EGADS_SUCCESS;
}
