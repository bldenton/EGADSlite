/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Attribute Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "egadsTypes.h"
#include "egadsInternals.h"

#ifdef WIN32
#define snprintf _snprintf
#endif


#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])



extern int EG_fullAttrs( const egObject *obj );
extern int EG_evaluate( const egObject *geom, /*@null@*/ const double *param,
                        double *results );
extern int EG_getTopology( const egObject *topo, egObject **geom, int *oclass,
                           int *type, /*@null@*/ double *limits, int *nChildren,
                           egObject ***children, int **senses );



int
EG_attributePrint(const egObject *obj)
{
  int     i, j;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_SUCCESS;
  
  printf("\n Attributes:\n");
  for (i = 0; i < attrs->nattrs; i++) {
    printf("    %s: ", attrs->attrs[i].name);
    if (attrs->attrs[i].type == ATTRINT) {
      if (attrs->attrs[i].length <= 1) {
        printf("%d\n", attrs->attrs[i].vals.integer);
      } else {
        for (j = 0; j < attrs->attrs[i].length; j++)
          printf("%d ", attrs->attrs[i].vals.integers[j]);
        printf("\n");
      }
    } else if ((attrs->attrs[i].type == ATTRREAL) ||
               (attrs->attrs[i].type == ATTRCSYS)) {
      if (attrs->attrs[i].length <= 1) {
        printf("%lf\n", attrs->attrs[i].vals.real);
      } else {
        for (j = 0; j < attrs->attrs[i].length; j++)
          printf("%lf ", attrs->attrs[i].vals.reals[j]);
        printf("\n");
      }
    } else if (attrs->attrs[i].type == ATTRPTR) {
#if defined(WIN32) && defined(_OCC64)
      printf("%llx\n", (long long) attrs->attrs[i].vals.string);
#else
      printf("%lx\n", (long) attrs->attrs[i].vals.string);
#endif
    } else {
      printf("%s\n", attrs->attrs[i].vals.string);
    }
  }

  return EGADS_SUCCESS;
}


void
EG_attrBuildSeq(egAttrs *attrs)
{
  int       i, j, l, n, snum, *hit, nospace = 0, nseqs = 0;
  char      *root, *newname;
  egAttr    *attr;
  egAttrSeq *seqs = NULL, *tmp;
  
  /* cleanup what might have been attached */
  for (i = 0; i < attrs->nseqs; i++) {
    EG_free(attrs->seqs[i].root);
    EG_free(attrs->seqs[i].attrSeq);
  }
  if (attrs->seqs != NULL) EG_free(attrs->seqs);
  attrs->seqs  = NULL;
  attrs->nseqs = 0;
  if (attrs->nattrs == 0) return;
  
  hit = (int *) EG_alloc(attrs->nattrs*sizeof(int));
  if (hit == NULL) {
    printf(" EGADS Internal: Malloc on %d attributes!\n", attrs->nattrs);
    return;
  }
  for (i = 0; i < attrs->nattrs; i++) hit[i] = 0;
  
  /* build the sequence structure */
  for (i = 0; i < attrs->nattrs; i++) {
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

  attrs->nseqs = nseqs;
  attrs->seqs  = seqs;
}


static int
EG_constructCSys(egObject *obj, int outLevel, int len, const double *reals,
                 double *csys)
{
  int    stat, idir, flip, oclass, mtype, n, *senses;
  double val, pos[3], dir1[3], dir2[3], dir3[3], result[18];
  egObject *geom, **children;
  
  pos[0]  = pos[1]  = pos[2]  = 0.0;
  dir1[0] = dir1[1] = dir1[2] = 0.0;
  dir2[0] = dir2[1] = dir2[2] = 0.0;
  if (len == 9) {
    pos[0]  = reals[0];
    pos[1]  = reals[1];
    pos[2]  = reals[2];
    dir1[0] = reals[3];
    dir1[1] = reals[4];
    dir1[2] = reals[5];
    dir2[0] = reals[6];
    dir2[1] = reals[7];
    dir2[2] = reals[8];
  } else if (len == 6) {
    if ((obj->oclass == FACE) || (obj->oclass == SURFACE)) {
      stat = EG_evaluate(obj, reals, result);
      if (stat != EGADS_SUCCESS) return stat;
      pos[0]  = result[0];
      pos[1]  = result[1];
      pos[2]  = result[2];
      dir2[0] = result[3];
      dir2[1] = result[4];
      dir2[2] = result[5];
      dir3[0] = result[6];
      dir3[1] = result[7];
      dir3[2] = result[8];
      CROSS(dir1, dir2, dir3);
      if ((obj->oclass == FACE) && (obj->mtype == SREVERSE)) {
        dir1[0] = -dir1[0];
        dir1[1] = -dir1[1];
        dir1[2] = -dir1[2];
      }
      dir1[0] *= reals[2];
      dir1[1] *= reals[2];
      dir1[2] *= reals[2];
      dir2[0]  = reals[3];
      dir2[1]  = reals[4];
      dir2[2]  = reals[5];
    } else if (obj->oclass == NODE) {
      stat = EG_getTopology(obj, &geom, &oclass, &mtype, pos,
                            &n, &children, &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_getTopology = %d (EG_attributeAdd)!\n", stat);
        return stat;
      }
      dir1[0] = reals[0];
      dir1[1] = reals[1];
      dir1[2] = reals[2];
      dir2[0] = reals[3];
      dir2[1] = reals[4];
      dir2[2] = reals[5];
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: CSys with len = %d (EG_attributeAdd)!\n", len);
      return EGADS_ATTRERR;
    }
  } else if (len == 5) {
    if ((obj->oclass == EDGE) || (obj->oclass == CURVE)) {
      stat = EG_evaluate(obj, reals, result);
      if (stat != EGADS_SUCCESS) return stat;
      pos[0]  = result[0];
      pos[1]  = result[1];
      pos[2]  = result[2];
      dir1[0] = result[3]*reals[1];
      dir1[1] = result[4]*reals[1];
      dir1[2] = result[5]*reals[1];
      dir2[0] = reals[2];
      dir2[1] = reals[3];
      dir2[2] = reals[4];
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: CSys with len = %d (EG_attributeAdd)!\n", len);
      return EGADS_ATTRERR;
    }
  } else if (len == 3) {
    if ((obj->oclass == FACE) || (obj->oclass == SURFACE)) {
      val  = reals[2];
      flip = 1;
      if (val < 0.0) flip = -1;
      val *= flip;
      idir = val + 0.0001;
      if ((idir < 1) || (idir > 4)) {
        if (outLevel > 0)
          printf(" EGADS Error: CSys with dir = %lf (EG_attributeAdd)!\n",
                 reals[2]);
        return EGADS_RANGERR;
      }
      stat = EG_evaluate(obj, reals, result);
      if (stat != EGADS_SUCCESS) return stat;
      pos[0]  = result[0];
      pos[1]  = result[1];
      pos[2]  = result[2];
      dir2[0] = result[3];
      dir2[1] = result[4];
      dir2[2] = result[5];
      dir3[0] = result[6];
      dir3[1] = result[7];
      dir3[2] = result[8];
      CROSS(dir1, dir2, dir3);
      if ((obj->oclass == FACE) && (obj->mtype == SREVERSE)) {
        dir1[0] = -dir1[0];
        dir1[1] = -dir1[1];
        dir1[2] = -dir1[2];
      }
      dir1[0] *= flip;
      dir1[1] *= flip;
      dir1[2] *= flip;
      if (idir == 1) {
        dir2[0] =  result[3];
        dir2[1] =  result[4];
        dir2[2] =  result[5];
      } else if (idir == 2) {
        dir2[0] =  result[6];
        dir2[1] =  result[7];
        dir2[2] =  result[8];
      } else if (idir == 3) {
        dir2[0] = -result[3];
        dir2[1] = -result[4];
        dir2[2] = -result[5];
      } else {
        dir2[0] = -result[6];
        dir2[1] = -result[7];
        dir2[2] = -result[8];
      }
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: CSys with len = %d (EG_attributeAdd)!\n", len);
      return EGADS_ATTRERR;
    }
  } else {
    if (outLevel > 0)
      printf(" EGADS Error: CSys with len = %d (EG_attributeAdd)!\n", len);
    return EGADS_ATTRERR;
  }
  
  /* normalize the directions */
  val = sqrt(DOT(dir1, dir1));
  if (val == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: CSys Dir1 is zero (EG_attributeAdd)!\n");
    return EGADS_DEGEN;
  }
  dir1[0] /= val;
  dir1[1] /= val;
  dir1[2] /= val;
  val = sqrt(DOT(dir2, dir2));
  if (val == 0.0) {
    if (outLevel > 0)
      printf(" EGADS Error: CSys Dir2 is zero (EG_attributeAdd)!\n");
    return EGADS_DEGEN;
  }
  dir2[0] /= val;
  dir2[1] /= val;
  dir2[2] /= val;
  
  /* check for orthogonality */
  if (fabs(DOT(dir1, dir2)) > 1.e-5) {
    if (outLevel > 0)
      printf(" EGADS Error: CSys Dirs Not Orthongonal (EG_attributeAdd)!\n");
    return EGADS_NOTORTHO;
  }
  
  CROSS(dir3, dir1, dir2);
  csys[ 0] = pos[0];
  csys[ 1] = pos[1];
  csys[ 2] = pos[2];
  csys[ 3] = dir1[0];
  csys[ 4] = dir1[1];
  csys[ 5] = dir1[2];
  csys[ 6] = dir2[0];
  csys[ 7] = dir2[1];
  csys[ 8] = dir2[2];
  csys[ 9] = dir3[0];
  csys[10] = dir3[1];
  csys[11] = dir3[2];
  return EGADS_SUCCESS;
}


static void
EG_adjustCSys(egObject *obj, int len, const double *xform, double *data)
{
  int    i;
  double csys[12], val;
  
  /* transform the CSys */
  for (i = 0; i < 12; i++) csys[i] = data[len+i];
  data[len   ] = csys[0]*xform[ 0] + csys[ 1]*xform[ 1] +
                 csys[2]*xform[ 2] +          xform[ 3];
  data[len+ 1] = csys[0]*xform[ 4] + csys[ 1]*xform[ 5] +
                 csys[2]*xform[ 6] +          xform[ 7];
  data[len+ 2] = csys[0]*xform[ 8] + csys[ 1]*xform[ 9] +
                 csys[2]*xform[10] +          xform[11];
  data[len+ 3] = csys[3]*xform[ 0] + csys[ 4]*xform[ 1] + csys[ 5]*xform[ 2];
  data[len+ 4] = csys[3]*xform[ 4] + csys[ 4]*xform[ 5] + csys[ 5]*xform[ 6];
  data[len+ 5] = csys[3]*xform[ 8] + csys[ 4]*xform[ 9] + csys[ 5]*xform[10];
  val          = sqrt(data[len+ 3]*data[len+ 3] + data[len+ 4]*data[len+ 4] +
                      data[len+ 5]*data[len+ 5]);
  if (val != 0.0) {
    data[len+ 3] /= val;
    data[len+ 4] /= val;
    data[len+ 5] /= val;
  }
  data[len+ 6] = csys[6]*xform[ 0] + csys[ 7]*xform[ 1] + csys[ 8]*xform[ 2];
  data[len+ 7] = csys[6]*xform[ 4] + csys[ 7]*xform[ 5] + csys[ 8]*xform[ 6];
  data[len+ 8] = csys[6]*xform[ 8] + csys[ 7]*xform[ 9] + csys[ 8]*xform[10];
  val          = sqrt(data[len+ 6]*data[len+ 6] + data[len+ 7]*data[len+ 7] +
                      data[len+ 8]*data[len+ 8]);
  if (val != 0.0) {
    data[len+ 6] /= val;
    data[len+ 7] /= val;
    data[len+ 8] /= val;
  }
  data[len+ 9] = csys[9]*xform[ 0] + csys[10]*xform[ 1] + csys[11]*xform[ 2];
  data[len+10] = csys[9]*xform[ 4] + csys[10]*xform[ 5] + csys[11]*xform[ 6];
  data[len+11] = csys[9]*xform[ 8] + csys[10]*xform[ 9] + csys[11]*xform[10];
  val          = sqrt(data[len+ 9]*data[len+ 9] + data[len+10]*data[len+10] +
                      data[len+11]*data[len+11]);
  if (val != 0.0) {
    data[len+ 9] /= val;
    data[len+10] /= val;
    data[len+11] /= val;
  }

  /* transform the original data */
  if (len == 9) {
    for (i = 0; i < 9; i++) data[i] = data[len+i];
  } else if (len == 6) {
    if (obj->oclass == NODE) {
      for (i = 0; i < 6; i++) data[i] = data[len+i+3];
    } else {
      for (i = 0; i < 3; i++) data[i+3] = data[len+i+6];
    }
  } else if (len == 5) {
    for (i = 0; i < 3; i++) data[i+2] = data[len+i+6];
  } else if (len != 3) {
    printf(" EGADS Internal: EG_adjustCSys with len = %d!\n", len);
  }

}


int
EG_attributeNumSeq(const egObject *obj, const char *name, int *num)
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


static int
EG_attributeMerge(int flg, egObject *obj, const char *name, int atype, int lenx,
                  /*@null@*/ const int  *ints, /*@null@*/ const double *reals,
                  /*@null@*/ const char *str)
{
  int     i, stat, len, length, outLevel, find = -1;
  double  csys[12];
  egAttr  *attr;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (EG_sameThread(obj))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(obj);

  /* special CSYS copy flag */
  len = lenx;
  if ((atype == ATTRCSYS) && (lenx < 0)) len = -lenx;

  if (name == NULL) {
    if (outLevel > 0) 
      printf(" EGADS Error: NULL Name (EG_attributeAdd)!\n");
    return EGADS_NONAME;
  }
  length = strlen(name);
  if (flg == 0)
    for (i = 0; i < length; i++)
      if (name[i] <= ' ') {
        length = 0;
        break;
      }
  if (length == 0) {
    if (outLevel > 0) 
      printf(" EGADS Error: BAD Name (EG_attributeAdd)!\n");
    return EGADS_INDEXERR;
  }
  if ((atype != ATTRINT)  && (atype != ATTRREAL) && (atype != ATTRSTRING) &&
      (atype != ATTRCSYS) && (atype != ATTRPTR)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Bad Attr Type (%d) for %s (EG_attributeAdd)!\n",
             atype, name);
    return EGADS_INDEXERR;
  }
  if (((atype == ATTRINT)    && (ints  == NULL)) ||
      ((atype == ATTRREAL)   && (reals == NULL)) ||
      ((atype == ATTRCSYS)   && (reals == NULL)) ||
      ((atype == ATTRSTRING) && (str   == NULL))) {
    if (outLevel > 0) 
      printf(" EGADS Error: NULL data for %s  type = %d (EG_attributeAdd)!\n",
             name, atype);
    return EGADS_NODATA;
  }
  if ((len <= 0) && (atype != ATTRSTRING) && (atype != ATTRPTR)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Bad Attr Length (%d) for %s (EG_attributeAdd)!\n",
             len, name);
    return EGADS_INDEXERR;
  }
  
  attrs = (egAttrs *) obj->attrs;
  if ((attrs != NULL) && (flg == 0)) {
    for (i = 0; i < attrs->nseqs; i++)
      if (strcmp(attrs->seqs[i].root,name) == 0) {
        find = i;
        break;
      }
    if (find != -1) {
      if (outLevel > 0)
        printf(" EGADS Error: Sequenced Attribute for %s (EG_attributeAdd)!\n",
               name);
      return EGADS_SEQUERR;
    }
  }
  
  /* get the CSys */
  if ((atype == ATTRCSYS) && (reals != NULL) && (lenx > 0)) {
    stat = EG_constructCSys(obj, outLevel, len, reals, csys);
    if (stat != EGADS_SUCCESS) return stat;
  }

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
        if (outLevel > 0) 
          printf(" EGADS Error: Attrs MALLOC for %s (EG_attributeAdd)!\n",
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
        if (outLevel > 0) 
          printf(" EGADS Error: Attr MALLOC for %s (EG_attributeAdd)!\n",
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

  attrs->attrs[find].type   = atype;
  attrs->attrs[find].length = len;
  if (atype == ATTRINT) {
    if (ints != NULL)
      if (len == 1) {
        attrs->attrs[find].vals.integer = *ints;
      } else {
        attrs->attrs[find].vals.integers = (int *) EG_alloc(len*sizeof(int));
        if (attrs->attrs[find].vals.integers == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < len; i++)
            attrs->attrs[find].vals.integers[i] = ints[i];
        }
      }
  } else if (atype == ATTRREAL) {
    if (reals != NULL)
      if (len == 1) {
        attrs->attrs[find].vals.real = *reals;
      } else {
        attrs->attrs[find].vals.reals = (double *) EG_alloc(len*sizeof(double));
        if (attrs->attrs[find].vals.reals == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < len; i++)
            attrs->attrs[find].vals.reals[i] = reals[i];
        }
      }
  } else if (atype == ATTRCSYS) {
    if (reals != NULL) {
      if (lenx < 0) {
        attrs->attrs[find].vals.reals = (double *) EG_alloc(len*sizeof(double));
        if (attrs->attrs[find].vals.reals == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < len; i++)
            attrs->attrs[find].vals.reals[i] = reals[i];
        }
      } else {
        attrs->attrs[find].length += 12;
        attrs->attrs[find].vals.reals = (double *) EG_alloc((len+12)*sizeof(double));
        if (attrs->attrs[find].vals.reals == NULL) {
          attrs->attrs[find].length = 0;
        } else {
          for (i = 0; i < len; i++)
            attrs->attrs[find].vals.reals[i] = reals[i];
          /* fill in the actual CSys */
          for (i = 0; i < 12; i++)
            attrs->attrs[find].vals.reals[len+i] = csys[i];
        }
      }
    }
  } else if (atype == ATTRPTR) {
    attrs->attrs[find].length = -1;
    attrs->attrs[find].vals.string = (char *) str;
  } else {
    attrs->attrs[find].length = 0;
    attrs->attrs[find].vals.string = EG_strdup(str);
    if (attrs->attrs[find].vals.string != NULL) {
      attrs->attrs[find].length = strlen(attrs->attrs[find].vals.string);
    } else {
      if (outLevel > 0)
        printf(" EGADS Info: Attrs EG_strdup NULL for %s (EG_attributeAdd)!\n",
               name);
    }
  }
  return EGADS_SUCCESS;
}


int
EG_attributeAdd(egObject *obj, const char *name, int atype, int len,
                /*@null@*/ const int  *ints, /*@null@*/ const double *reals,
                /*@null@*/ const char *str)
{
  return EG_attributeMerge(0, obj, name, atype, len, ints, reals, str);
}


static int
EG_attrSameValue(egAttr *attr, int atype, int len, /*@null@*/ const int  *ints,
                 /*@null@*/ const double *reals, /*@null@*/ const char *str)
{
  int i;

  if (attr->type != atype) return EGADS_OUTSIDE;
  
  if (atype == ATTRINT) {
    if (len != attr->length) return EGADS_OUTSIDE;
    if (ints == NULL) return EGADS_OUTSIDE;
    if (len == 1) {
      if (attr->vals.integer != *ints) return EGADS_OUTSIDE;
    } else {
      for (i = 0; i < len; i++)
        if (attr->vals.integers[i] != ints[i]) return EGADS_OUTSIDE;
    }

  } else if (atype == ATTRREAL) {
    if (len != attr->length) return EGADS_OUTSIDE;
    if (reals == NULL) return EGADS_OUTSIDE;
    if (len == 1) {
      if (attr->vals.real != *reals) return EGADS_OUTSIDE;
    } else {
      for (i = 0; i < len; i++)
        if (attr->vals.reals[i] != reals[i]) return EGADS_OUTSIDE;
    }
  } else if (atype == ATTRCSYS) {
    if (len != attr->length-12) return EGADS_OUTSIDE;
    if (reals == NULL) return EGADS_OUTSIDE;
    for (i = 0; i < len; i++)
      if (attr->vals.reals[i] != reals[i]) return EGADS_OUTSIDE;
  } else if (atype == ATTRPTR) {
    if (str == NULL) return EGADS_OUTSIDE;
    if (attr->vals.string != (char *) str) return EGADS_OUTSIDE;
  } else {
    if (str == NULL) return EGADS_OUTSIDE;
    if (attr->vals.string == NULL) return EGADS_OUTSIDE;
    if (strcmp(attr->vals.string,str) != 0) return EGADS_OUTSIDE;
  }
  
  return EGADS_SUCCESS;
}


int
EG_attributeAddSeq(egObject *obj, const char *name, int atype, int len,
                   /*@null@*/ const int  *ints, /*@null@*/ const double *reals,
                   /*@null@*/ const char *str)
{
  int     i, stat, length, outLevel, seq, find = -1;
  char    *newname;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (EG_sameThread(obj))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(obj);

  if (name == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL Name (EG_attributeAddSeq)!\n");
    return EGADS_NONAME;
  }
  length = strlen(name);
  for (i = 0; i < length; i++)
    if (name[i] <= ' ') {
      length = 0;
      break;
    }
  if (length == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: BAD Name (EG_attributeAddSeq)!\n");
    return EGADS_INDEXERR;
  }
  if ((atype != ATTRINT)  && (atype != ATTRREAL) && (atype != ATTRSTRING) &&
      (atype != ATTRCSYS) && (atype != ATTRPTR)) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Attr Type (%d) for %s (EG_attributeAddSeq)!\n",
             atype, name);
    return EGADS_INDEXERR;
  }
  if (((atype == ATTRINT)    && (ints  == NULL)) ||
      ((atype == ATTRREAL)   && (reals == NULL)) ||
      ((atype == ATTRCSYS)   && (reals == NULL)) ||
      ((atype == ATTRSTRING) && (str   == NULL))) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL data for %s  type = %d (EG_attributeAddSeq)!\n",
             name, atype);
    return EGADS_NODATA;
  }
  if ((len <= 0) && ((atype == ATTRINT) || (atype == ATTRREAL))) {
    if (outLevel > 0)
      printf(" EGADS Error: Bad Attr Length (%d) for %s (EG_attributeAddSeq)!\n",
             len, name);
    return EGADS_INDEXERR;
  }
  attrs = (egAttrs *) obj->attrs;
  
  /* what is the sequence status? */
  if (attrs != NULL)
    for (i = 0; i < attrs->nseqs; i++)
      if (strcmp(attrs->seqs[i].root,name) == 0) {
        find = i;
        break;
      }
  
  if ((find == -1) || (attrs == NULL)) {
    /* no sequence */
    if (attrs != NULL)
      for (i = 0; i < attrs->nattrs; i++)
        if (strcmp(attrs->attrs[i].name,name) == 0) {
          find = i;
          break;
        }
    if ((find == -1) || (attrs == NULL))
      return EG_attributeMerge(1, obj, name, atype, len, ints, reals, str);
  
    /* do we have the same value? */
    stat = EG_attrSameValue(&attrs->attrs[find], atype, len, ints, reals, str);
    if (stat == 0) return EGADS_SUCCESS;
    newname = (char *) EG_alloc((length+8)*sizeof(char));
    if (newname == NULL) return EGADS_MALLOC;
    seq = 2;
    snprintf(newname, length+8, "%s %d", name, seq);
    stat = EG_attributeMerge(1, obj, newname, atype, len, ints, reals, str);
    EG_free(newname);
    if (stat != EGADS_SUCCESS) return stat;
    
    EG_attrBuildSeq(obj->attrs);
    return seq;
  }
  
  /* same values as something in the sequence? */
  for (i = 0; i < attrs->seqs[find].nSeq; i++) {
    stat = EG_attrSameValue(&attrs->attrs[attrs->seqs[find].attrSeq[i]],
                            atype, len, ints, reals, str);
    if (stat == 0) return EGADS_SUCCESS;
  }
  
  /* add to an existing sequence */
  newname = (char *) EG_alloc((length+8)*sizeof(char));
  if (newname == NULL) return EGADS_MALLOC;
  seq = attrs->seqs[find].nSeq+1;
  snprintf(newname, length+8, "%s %d", name, seq);
  stat = EG_attributeMerge(1, obj, newname, atype, len, ints, reals, str);
  EG_free(newname);
  if (stat != EGADS_SUCCESS) return stat;
  
  EG_attrBuildSeq(obj->attrs);
  return seq;
}


int
EG_attributeDel(egObject *obj, /*@null@*/ const char *name)
{
  int     i, j, k, outLevel, find = -1;
  char    *cptr;
  egAttrs *attrs;

  if (obj == NULL)               return EGADS_NULLOBJ;
  if (obj->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (obj->oclass == EMPTY)      return EGADS_EMPTY;
  if (obj->oclass == NIL)        return EGADS_EMPTY;
  if (obj->oclass == REFERENCE)  return EGADS_REFERCE;
  if (EG_sameThread(obj))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(obj);

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_SUCCESS;

  if (name == NULL) {

    /* delete all attributes associated with the object */
    obj->attrs = NULL;
    for (i = 0; i < attrs->nseqs; i++) {
      EG_free(attrs->seqs[i].root);
      EG_free(attrs->seqs[i].attrSeq);
    }
    if (attrs->seqs != NULL) EG_free(attrs->seqs);
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

  } else {

    /* are we sequenced? */
    
    for (i = 0; i < attrs->nseqs; i++)
      if (strcmp(attrs->seqs[i].root,name) == 0) {
        find = i;
        break;
      }
    if (find != -1) {
      /* delete all attribues in the sequence */
      for (i = attrs->seqs[find].nSeq-1; i >= 0; i--) {
        j = attrs->seqs[find].attrSeq[i];
        EG_free(attrs->attrs[j].name);
        if (attrs->attrs[j].type == ATTRINT) {
          if (attrs->attrs[j].length > 1)
            EG_free(attrs->attrs[j].vals.integers);
        } else if ((attrs->attrs[j].type == ATTRREAL) ||
                   (attrs->attrs[j].type == ATTRCSYS)) {
          if (attrs->attrs[j].length > 1)
            EG_free(attrs->attrs[j].vals.reals);
        } else if (attrs->attrs[j].type == ATTRSTRING) {
          EG_free(attrs->attrs[j].vals.string);
        }
        for (k = j+1; k < attrs->nattrs; k++)
          attrs->attrs[k-1] = attrs->attrs[k];
        attrs->nattrs -= 1;
      }
      EG_attrBuildSeq(attrs);
      return EGADS_SUCCESS;
    }
    
    /* delete the named attribute */
    for (i = 0; i < attrs->nattrs; i++)
      if (strcmp(attrs->attrs[i].name,name) == 0) {
        find = i;
        break;
      }
    if (find == -1) {
      if (outLevel > 0) 
        printf(" EGADS Error: No Attribute -> %s (EG_attributeDel)!\n",
               name);
      return EGADS_NOTFOUND;
    }
    cptr = attrs->attrs[find].name;
    if (attrs->attrs[find].type == ATTRINT) {
      if (attrs->attrs[find].length > 1) 
        EG_free(attrs->attrs[find].vals.integers);
    } else if ((attrs->attrs[find].type == ATTRREAL) ||
               (attrs->attrs[find].type == ATTRCSYS)) {
      if (attrs->attrs[find].length > 1) 
        EG_free(attrs->attrs[find].vals.reals);
    } else if (attrs->attrs[find].type == ATTRSTRING) {
      EG_free(attrs->attrs[find].vals.string);
    }
    for (i = find+1; i < attrs->nattrs; i++)
      attrs->attrs[i-1] = attrs->attrs[i];
    attrs->nattrs -= 1;

    /* are we sequenced? */
    j = strlen(name);
    for (i = 0; i < j; i++)
      if (name[i] == 32) {
        EG_attrBuildSeq(attrs);
        break;
      }
    EG_free(cptr);
  }

  return EGADS_SUCCESS;
}


int
EG_attributeNum(const egObject *obj, int *num)
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


int
EG_attributeGet(const egObject *obj, int index, const char **name, 
                int *atype, int *len, /*@null@*/ const int **ints, 
                /*@null@*/ const double **reals, /*@null@*/ const char **str)
{
  int     outLevel;
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
  outLevel = EG_outLevel(obj);

  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) {
    if (outLevel > 0) 
      printf(" EGADS Error: NULL Attributes (EG_attributeGet)!\n");
    return EGADS_INDEXERR;
  }
  if ((index < 1) || (index > attrs->nattrs)) {
    if (outLevel > 0) 
      printf(" EGADS Error: Index Error %d [1-%d] (EG_attributeGet)!\n",
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


int
EG_attributeRet(const egObject *obj, const char *name, int *atype, 
                int *len, /*@null@*/ const int **ints, 
                          /*@null@*/ const double **reals, 
                          /*@null@*/ const char **str)
{
  int     i, outLevel, index;
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
  outLevel = EG_outLevel(obj);

  if (name == NULL) {
    if (outLevel > 0) 
      printf(" EGADS Error: NULL Name (EG_attributeRet)!\n");
    return EGADS_NONAME;
  }
  attrs = (egAttrs *) obj->attrs;
  if (attrs == NULL) return EGADS_NOTFOUND;

  index = -1;
  for (i = 0; i < attrs->nattrs; i++)
    if (strcmp(attrs->attrs[i].name, name) == 0) {
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


int
EG_attributeRetSeq(const egObject *obj, const char *name, int index, int *atype,
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
    return EG_attributeRet(obj, name, atype, len, ints, reals, str);
  }
  if (index > attrs->seqs[find].nSeq) {
    printf(" EGADS Error: Index %d [1-%d] (EG_attributeRetSeq)!\n",
           index, attrs->seqs[find].nSeq);
    return EGADS_INDEXERR;
  }
  
  fullname = (char *) EG_alloc((length+8)*sizeof(char));
  if (fullname == NULL) return EGADS_MALLOC;
  snprintf(fullname, length+8, "%s %d", name, index);
  stat = EG_attributeRet(obj, fullname, atype, len, ints, reals, str);
  EG_free(fullname);
  
  return stat;
}


int
EG_attributeXDup(const egObject *src, /*@null@*/ const double *xform,
                       egObject *dst)
{
  int     i, j, k, l, n, stat, len, outLevel, fullAttr, freer, *ints;
  double  *reals;
  char    *str, *name;
  egAttr  *attr;
  egAttrs *sattrs, *dattrs;
  
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->oclass == EMPTY)      return EGADS_EMPTY;
  if (src->oclass == NIL)        return EGADS_EMPTY;
  if (src->oclass == REFERENCE)  return EGADS_REFERCE;
  outLevel = EG_outLevel(src);
  fullAttr = EG_fullAttrs(src);
  
  if (dst == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst (EG_attributeDup)!\n");
    return EGADS_NULLOBJ;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attributeDup)!\n");
    return EGADS_NOTOBJ;
  }
  
  /* remove any current attributes */
  if (fullAttr == 0) {
    dattrs = dst->attrs;
    if (dattrs != NULL) {
      dst->attrs = NULL;
      for (i = 0; i < dattrs->nseqs; i++) {
        EG_free(dattrs->seqs[i].root);
        EG_free(dattrs->seqs[i].attrSeq);
      }
      if (dattrs->seqs != NULL) EG_free(dattrs->seqs);
      dattrs->nseqs = 0;
      dattrs->seqs  = NULL;
      for (i = 0; i < dattrs->nattrs; i++) {
        EG_free(dattrs->attrs[i].name);
        if (dattrs->attrs[i].type == ATTRINT) {
          if (dattrs->attrs[i].length > 1) EG_free(dattrs->attrs[i].vals.integers);
        } else if ((dattrs->attrs[i].type == ATTRREAL) ||
                   (dattrs->attrs[i].type == ATTRCSYS)) {
          if (dattrs->attrs[i].length > 1) EG_free(dattrs->attrs[i].vals.reals);
        } else if (dattrs->attrs[i].type == ATTRSTRING) {
          EG_free(dattrs->attrs[i].vals.string);
        }
      }
      EG_free(dattrs->attrs);
      EG_free(dattrs);
    }
  }
  
  sattrs = src->attrs;
  if (sattrs == NULL) return EGADS_SUCCESS;
  for (n = i = 0; i < sattrs->nattrs; i++)
    if (sattrs->attrs[i].type != ATTRPTR) n++;
  if (n == 0) return EGADS_SUCCESS;
  
  dattrs = dst->attrs;
  if (dattrs == NULL) {

    /* copy the attributes */
    dattrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
    if (dattrs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: dst Malloc (EG_attributeDup)!\n");
      return EGADS_MALLOC;
    }
    dattrs->nattrs = 0;
    dattrs->attrs  = NULL;
    dattrs->nseqs  = 0;
    dattrs->seqs   = NULL;
    dst->attrs     = dattrs;
    attr           = (egAttr *) EG_alloc(n*sizeof(egAttr));
    if (attr == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: dst attr Malloc (EG_attributeDup)!\n");
      return EGADS_MALLOC;
    }
    for (k = i = 0; i < sattrs->nattrs; i++) {
      if (sattrs->attrs[i].type == ATTRPTR) continue;
      attr[k].name   = EG_strdup(sattrs->attrs[i].name);
      attr[k].length = sattrs->attrs[i].length;
      attr[k].type   = sattrs->attrs[i].type;
      if (attr[k].type == ATTRINT) {
        if (attr[k].length <= 1) {
          attr[k].vals.integer = sattrs->attrs[i].vals.integer;
        } else {
          attr[k].vals.integers = (int *) EG_alloc(attr[k].length*sizeof(int));
          if (attr[k].vals.integers == NULL) {
            attr[k].length = 0;
          } else {
            for (j = 0; j < attr[k].length; j++)
              attr[k].vals.integers[j] = sattrs->attrs[i].vals.integers[j];
          }
        }
      } else if (attr[k].type == ATTRREAL) {
        if (attr[k].length <= 1) {
          attr[k].vals.real = sattrs->attrs[i].vals.real;
        } else {
          attr[k].vals.reals = (double *) EG_alloc(attr[k].length*sizeof(double));
          if (attr[k].vals.reals == NULL) {
            attr[k].length = 0;
          } else {
            for (j = 0; j < attr[i].length; j++)
              attr[k].vals.reals[j] = sattrs->attrs[i].vals.reals[j];
          }
        }
      } else if (attr[k].type == ATTRCSYS) {
        attr[k].vals.reals = (double *) EG_alloc(attr[k].length*sizeof(double));
        if (attr[k].vals.reals == NULL) {
          attr[k].length = 0;
        } else {
          for (j = 0; j < attr[k].length; j++)
            attr[k].vals.reals[j] = sattrs->attrs[i].vals.reals[j];
          if (xform != NULL)
            EG_adjustCSys(dst, attr[k].length-12, xform, attr[k].vals.reals);
        }
      } else {
        attr[k].vals.string = EG_strdup(sattrs->attrs[i].vals.string);
      }
      k++;
    }
    dattrs->nattrs = n;
    dattrs->attrs  = attr;
    
  } else {
    
    /* merge the attributes */
    for (i = 0; i < sattrs->nattrs; i++) {
      if (sattrs->attrs[i].type == ATTRPTR) continue;
      freer = 0;
      len   = sattrs->attrs[i].length;
      ints  = NULL;
      reals = NULL;
      str   = NULL;
      if (sattrs->attrs[i].type == ATTRINT) {
        if (len <= 1) {
          ints = &sattrs->attrs[i].vals.integer;
        } else {
          ints = sattrs->attrs[i].vals.integers;
        }
      } else if (sattrs->attrs[i].type == ATTRREAL) {
        if (len <= 1) {
          reals = &sattrs->attrs[i].vals.real;
        } else {
          reals = sattrs->attrs[i].vals.reals;
        }
      } else if (sattrs->attrs[i].type == ATTRCSYS) {
        reals = sattrs->attrs[i].vals.reals;
        if (xform != NULL) {
          reals = (double *) EG_alloc(len*sizeof(double));
          if (reals == NULL) {
            if (outLevel > 0)
              printf(" EGADS Warning: Malloc of = %d CSYS (EG_attributeDup)!\n",
                     len);
            continue;
          }
          for (j = 0; j < len; j++) reals[j] = sattrs->attrs[i].vals.reals[j];
          EG_adjustCSys(dst, len-12, xform, reals);
          freer = 1;
        }
        /* flag an ATTRCSYS copy */
        len = -len;
      } else {
        str = sattrs->attrs[i].vals.string;
      }
      name = sattrs->attrs[i].name;
      l    = strlen(name);
      for (j = 0; j < l; j++)
        if (name[j] == 32) break;
      if (j != l) {
        name    = EG_strdup(sattrs->attrs[i].name);
        if (name == NULL)  {
          if (outLevel > 0)
            printf(" EGADS Warning: Malloc of = %s (EG_attributeDup)!\n",
                   sattrs->attrs[i].name);
          if (freer == 1) EG_free(reals);
          continue;
        }
        name[j] = 0;
      }
      stat = EG_attributeAddSeq(dst, name, sattrs->attrs[i].type, len,
                                ints, reals, str);
      if (stat < EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: EG_attributeAddSeq = %d (EG_attributeDup)!\n",
                 stat);
      if (name != sattrs->attrs[i].name) EG_free(name);
      if (freer == 1) EG_free(reals);
    }
  }
  EG_attrBuildSeq(dattrs);
  
  return EGADS_SUCCESS;
}


int
EG_attributeDup(const egObject *src, egObject *dst)
{
  if (EG_sameThread(src)) return EGADS_CNTXTHRD;
  if (EG_sameThread(dst)) return EGADS_CNTXTHRD;
  return EG_attributeXDup(src, NULL, dst);
}


int
EG_attributeCommon(const egObject *src, egObject *dst)
{
  int          i, j, hit, stat, outLevel, snum, dnum, stype, slen, dtype, dlen;
  double       reals[3];
  const char   *name, *dstr, *sstr;
  const int    *dints,  *sints;
  const double *dreals, *sreals;
  
  if (src == NULL)               return EGADS_NULLOBJ;
  if (src->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (src->oclass == EMPTY)      return EGADS_EMPTY;
  if (src->oclass == NIL)        return EGADS_EMPTY;
  if (src->oclass == REFERENCE)  return EGADS_REFERCE;
  if (EG_sameThread(src))        return EGADS_CNTXTHRD;
  if (EG_sameThread(dst))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(src);
  
  if (dst == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst (EG_attributeCommon)!\n");
    return EGADS_NULLOBJ;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attributeCommon)!\n");
    return EGADS_NOTOBJ;
  }
  stat = EG_attributeNum(src, &snum);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_attributeNum src = %d (EG_attributeCommon)!\n",
             stat);
    return stat;
  }
  stat = EG_attributeNum(dst, &dnum);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_attributeNum dst = %d (EG_attributeCommon)!\n",
             stat);
    return stat;
  }
  if ((dnum == 0) && (snum == 0)) return EGADS_SUCCESS;
  
  /* handle EGADS attributes first */
  hit  = 0;
  stat = EG_attributeRet(dst, ".tParams", &dtype, &dlen, &dints, &dreals, &dstr);
  if (stat == EGADS_SUCCESS) {
    if ((dtype != ATTRREAL) || (dlen != 3)) {
      EG_attributeDel(dst, ".tParams");
      hit = 1;
    } else {
      stat = EG_attributeRet(src, ".tParams", &stype, &slen,
                             &sints, &sreals, &sstr);
      if (stat == EGADS_SUCCESS) {
        if ((stype == ATTRREAL) && (slen == 3)) {
          if ((sreals[0] != dreals[0]) || (sreals[1] != dreals[1]) ||
              (sreals[2] != dreals[2])) {
            reals[0] = sreals[0];
            if (dreals[0] < reals[0]) reals[0] = dreals[0];
            reals[1] = sreals[1];
            if (dreals[1] < reals[1]) reals[1] = dreals[1];
            reals[2] = sreals[2];
            if (dreals[2] < reals[2]) reals[2] = dreals[2];
            stat = EG_attributeAdd(dst, ".tParams", dtype, dlen, NULL, reals,
                                   NULL);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_attributeAdd tPs = %d (EG_attributeCommon)!\n",
                       stat);
              return stat;
            }
          }
        }
      }
    }
  } else {
    hit = 1;
  }
  if (hit == 1) {
    stat = EG_attributeRet(src, ".tParams", &stype, &slen,
                           &sints, &sreals, &sstr);
    if (stat == EGADS_SUCCESS)
      if ((stype == ATTRREAL) && (slen == 3)) {
        stat = EG_attributeAdd(dst, ".tParams", stype, slen, NULL, sreals,
                               NULL);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_attributeAdd TPs = %d (EG_attributeCommon)!\n",
                   stat);
          return stat;
        }
      }
  }
  
  hit  = 0;
  stat = EG_attributeRet(dst, ".tParam", &dtype, &dlen, &dints, &dreals, &dstr);
  if (stat == EGADS_SUCCESS) {
    if ((dtype != ATTRREAL) || (dlen != 3)) {
      EG_attributeDel(dst, ".tParam");
      hit = 1;
    } else {
      stat = EG_attributeRet(src, ".tParam", &stype, &slen,
                             &sints, &sreals, &sstr);
      if (stat == EGADS_SUCCESS) {
        if ((stype == ATTRREAL) && (slen == 3)) {
          if ((sreals[0] != dreals[0]) || (sreals[1] != dreals[1]) ||
              (sreals[2] != dreals[2])) {
            reals[0] = sreals[0];
            if (dreals[0] < reals[0]) reals[0] = dreals[0];
            reals[1] = sreals[1];
            if (dreals[1] < reals[1]) reals[1] = dreals[1];
            reals[2] = sreals[2];
            if (dreals[2] < reals[2]) reals[2] = dreals[2];
            stat = EG_attributeAdd(dst, ".tParam", dtype, dlen, NULL, reals,
                                   NULL);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: EG_attributeAdd tP = %d (EG_attributeCommon)!\n",
                       stat);
              return stat;
            }
          }
        }
      }
    }
  } else {
    hit = 1;
  }
  if (hit == 1) {
    stat = EG_attributeRet(src, ".tParam", &stype, &slen,
                           &sints, &sreals, &sstr);
    if (stat == EGADS_SUCCESS)
      if ((stype == ATTRREAL) && (slen == 3)) {
        stat = EG_attributeAdd(dst, ".tParam", stype, slen, NULL, sreals,
                               NULL);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_attributeAdd TP = %d (EG_attributeCommon)!\n",
                   stat);
          return stat;
        }
      }
  }
  
  stat = EG_attributeRet(dst, ".nPos", &dtype, &dlen, &dints, &dreals, &dstr);
  if (stat == EGADS_SUCCESS) {
    if ((dtype != ATTRINT) || (dlen != 1)) {
      EG_attributeDel(dst, ".nPos");
    } else {
      stat = EG_attributeRet(src, ".nPos", &stype, &slen,
                             &sints, &sreals, &sstr);
      if (stat == EGADS_SUCCESS) {
        if ((stype == ATTRINT) && (slen == 1)) {
          hit  = dints[0] + sints[0] + 1;
          stat = EG_attributeAdd(dst, ".nPos", dtype, dlen, &hit, reals, NULL);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: EG_attributeAdd nP = %d (EG_attributeCommon)!\n",
                     stat);
            return stat;
          }
        } else {
          EG_attributeDel(dst, ".nPos");
        }
      } else {
        EG_attributeDel(dst, ".nPos");
      }
    }
  }

  /* remove all other EGADS attributes */
  do {
    hit  = 0;
    stat = EG_attributeNum(dst, &dnum);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_attributeNum Dst = %d (EG_attributeCommon)!\n",
               stat);
      return stat;
    }
    for (i = 1; i <= dnum; i++) {
      stat = EG_attributeGet(dst, i, &name, &dtype, &dlen,
                             &dints, &dreals, &dstr);
      if (stat != EGADS_SUCCESS)        continue;
      if (name[0] != '.')               continue;
      if (strcmp(name,".tParams") == 0) continue;
      if (strcmp(name,".tParam")  == 0) continue;
      if (strcmp(name,".nPos")    == 0) continue;
      stat = EG_attributeDel(dst, name);
      if (stat != EGADS_SUCCESS) continue;
      hit = 1;
      break;
    }
  } while (hit == 1);
  
  /* find common attributes -- delete others */
  do {
    hit  = 0;
    stat = EG_attributeNum(dst, &dnum);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_attributeNum DSt = %d (EG_attributeCommon)!\n",
               stat);
      return stat;
    }
    for (i = 1; i <= dnum; i++) {
      stat = EG_attributeGet(dst, i, &name, &dtype, &dlen,
                             &dints, &dreals, &dstr);
      if (stat != EGADS_SUCCESS)      continue;
      if (name[0] == '.')             continue;
      stat = EG_attributeRet(src, name, &stype, &slen, &sints, &sreals, &sstr);
      if (stat  != EGADS_SUCCESS)     goto bail;
      if (stype != dtype)             goto bail;
      if (stype == ATTRINT) {
        if (slen != dlen)             goto bail;
        if ((sints == NULL) ||
            (dints == NULL))          goto bail;
        for (j = 0; j < slen; j++)
          if (dints[j] != sints[j])   goto bail;
      } else if ((stype == ATTRREAL) || (stype == ATTRCSYS)) {
        if (slen != dlen)             goto bail;
        if ((sreals == NULL) ||
            (dreals == NULL))         goto bail;
        for (j = 0; j < slen; j++)
          if (dreals[j] != sreals[j]) goto bail;
      } else if (stype == ATTRPTR) {
        if (dstr != sstr)             goto bail;
      } else {
        if ((dstr == NULL) ||
            (sstr == NULL))           goto bail;
        if (strcmp(dstr,sstr) != 0)   goto bail;
      }
      continue;
    bail:
      stat = EG_attributeDel(dst, name);
      if (stat != EGADS_SUCCESS) continue;
      hit = 1;
      break;
    }
  } while (hit == 1);
  
  return EGADS_SUCCESS;
}
