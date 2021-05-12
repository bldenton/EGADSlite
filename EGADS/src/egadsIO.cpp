/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Load & Save Functions
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

//#define STEPATTRS

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "egadsClasses.h"
#include <IGESControl_Controller.hxx>
#include <IGESBasic_Name.hxx>
#ifdef STEPATTRS
#include <StepBasic_Product.hxx>
#include <StepBasic_ProductDefinition.hxx>
#include <StepBasic_ProductDefinitionFormation.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <StepRepr_NextAssemblyUsageOccurrence.hxx>
#include <StepRepr_RepresentationItem.hxx>
#endif
#include <Interface_Graph.hxx>
#include <Interface_EntityIterator.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>


class egadsLabel
{
public:
  TopoDS_Shape shape;
  char        *shapeName;
};


#ifdef WIN32
#define snprintf _snprintf
#endif


#define INTERIM

#define UVTOL    1.e-4


  extern "C" void EG_revision( int *major, int *minor, char **OCCrev );
  extern "C" void EG_initOCC( );
  extern "C" int  EG_destroyTopology( egObject *topo );
  extern "C" int  EG_fullAttrs( const egObject *obj );
  extern "C" void EG_attrBuildSeq( egAttrs *attrs );
  extern "C" void EG_readAttrs( egObject *obj, int nattr, FILE *fp );
  extern "C" void EG_writeAttr( egAttrs *attrs, FILE *fp );
  extern "C" int  EG_writeNumAttr( egAttrs *attrs );
  extern "C" int  EG_dereferenceTopObj( egObject *object,
                                        /*@null@*/ const egObject *ref );

  extern "C" int  EG_loadModel( egObject *context, int bflg, const char *name, 
                                egObject **model );
  extern "C" int  EG_saveModel( const egObject *model, const char *name );
  extern "C" int  EG_saveTess( egObject *tess, const char *name );
  extern "C" int  EG_loadTess( egObject *body, const char *name,
                               egObject **tess );

  extern "C" int  EG_attributeAdd( egObject *obj, const char *name, int type,
                                   int len, /*@null@*/ const int    *ints,
                                            /*@null@*/ const double *reals,
                                            /*@null@*/ const char   *str );
  extern "C" int  EG_statusTessBody( egObject *tess, egObject **body,
                                     int *state, int *npts );
  extern "C" int  EG_getBodyTopos( const egObject *body, egObject *src,
                                   int oclass, int *ntopo, egObject ***topos );
  extern "C" int  EG_getTessEdge( const egObject *tess, int indx, int *len,
                                  const double **xyz, const double **t );
  extern "C" int  EG_getTessFace( const egObject *tess, int indx, int *len,
                                  const double **xyz, const double **uv,
                                  const int **ptype, const int **pindex,
                                  int *ntri, const int **tris, const int **tric );
  extern "C" int  EG_initTessBody( egObject *object, egObject **tess );
  extern "C" int  EG_setTessEdge( const egObject *tess, int index, int len,
                                  const double *xyz, const double *t );
  extern "C" int  EG_setTessFace( const egObject *tess, int index, int len,
                                  const double *xyz, const double *uv,
                                  int ntri, const int *tris );
  extern "C" int  EG_isSame ( const egObject *obj1, const egObject *obj2 );
  extern "C" int  EG_invEvaluate( const egObject *obj, double *xyz,
                                  double *param, double *results );
  extern "C" int  EG_writeEBody( const egObject *EBody, FILE *fp );
  extern "C" int  EG_readEBody( FILE *fp, egObject *body, egObject **EBody );

  extern     void EG_splitPeriodics( egadsBody *body );
  extern     void EG_splitMultiplicity( egadsBody *body, int outLevel );
  extern     int  EG_traverseBody( egObject *context, int i, egObject *bobj, 
                                   egObject *topObj, egadsBody *body,
                                   int *nerr );


void
EG_revision(int *major, int *minor, char **OCCrev)
{
#ifdef INTERIM
  static char OCCrevStr[42];
  
  snprintf(OCCrevStr, 41, "Interim Release with OpenCASCADE %d.%d.%d",
           OCC_VERSION_MAJOR, OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
#else
  static char OCCrevStr[26];
  
  snprintf(OCCrevStr, 25, "with OpenCASCADE %d.%d.%d", OCC_VERSION_MAJOR,
           OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
#endif
  *major  = EGADSMAJOR;
  *minor  = EGADSMINOR;
  *OCCrev = OCCrevStr;
}


void
EG_attriBodyTrav(const egObject *obj, const TopoDS_Shape& shape, egadsBody *pbody)
{
  if (obj->blind == NULL) return;

  if (obj->oclass == NODE) {
  
    int index = pbody->nodes.map.FindIndex(shape);
    if (index != 0)
      EG_attributeDup(obj, pbody->nodes.objs[index-1]);
    else
      printf(" EGADS Internal: Dropping Node attributes!\n");
  
  } else if (obj->oclass == EDGE) {
  
    egadsEdge *pedge = (egadsEdge *) obj->blind;
    int index = pbody->edges.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->edges.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Edge attributes!\n");
    }

    TopoDS_Vertex V1, V2;
    TopExp::Vertices(TopoDS::Edge(shape), V1, V2);

    EG_attriBodyTrav(pedge->nodes[0], V1, pbody);
    if (obj->mtype == TWONODE)
      EG_attriBodyTrav(pedge->nodes[1], V2, pbody);
  
  } else if (obj->oclass == LOOP) {
  
    egadsLoop *ploop = (egadsLoop *) obj->blind;
    int index = pbody->loops.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->loops.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Loop attributes!\n");
    }

    BRepTools_WireExplorer Exp(ploop->loop);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(ploop->edges[i], Exp.Current(), pbody);
    }
  
  } else if (obj->oclass == FACE) {
  
    egadsFace *pface = (egadsFace *) obj->blind;
    int index = pbody->faces.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->faces.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Face attributes!\n");
    }

    TopExp_Explorer Exp(pface->face, TopAbs_WIRE);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(pface->loops[i], Exp.Current(), pbody);
    }
  
  } else if (obj->oclass == SHELL) {
  
    egadsShell *pshell = (egadsShell *) obj->blind;
    int index = pbody->shells.map.FindIndex(shape);
    if (index != 0) {
      EG_attributeDup(obj, pbody->shells.objs[index-1]);
    } else {
      printf(" EGADS Internal: Dropping Shell attributes!\n");
    }

    TopExp_Explorer Exp(pshell->shell, TopAbs_FACE);
    for (int i = 0; Exp.More(); Exp.Next(), i++) {
      EG_attriBodyTrav(pshell->faces[i], Exp.Current(), pbody);
    }
  }
}


int
EG_attriBodyDup(const egObject *src, egObject *dst)
{
  int          i, j, n, stat, nents, nattr;
  double       tmin, tmax, trange[2], xyz[3];
  egObject     *aobj, *dobj, *geom;
  egAttrs      *attrs;
  TopoDS_Shape shape;
  
  if ((src == NULL) || (dst == NULL)) return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)     return EGADS_NOTOBJ;
  if  (src->oclass < NODE)            return EGADS_NOTTOPO;
  if  (src->blind == NULL)            return EGADS_NODATA;
  int outLevel = EG_outLevel(src);
  int fullAttr = EG_fullAttrs(src);
  
  if (src->oclass == MODEL) {
    egadsModel *pmdl = (egadsModel *) src->blind;
    for (i = 0; i < pmdl->nbody; i++) {
      stat = EG_attriBodyDup(pmdl->bodies[i], dst);
      if (stat != EGADS_SUCCESS) return stat;
    }
    return EGADS_SUCCESS;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attriBodyDup)!\n");
    return EGADS_NOTOBJ;
  }
  if (dst->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not a BODY (EG_attriBodyDup)!\n");
    return EGADS_NOTBODY;
  }
  if (dst->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst pointer (EG_attriBodyDup)!\n");
    return EGADS_NODATA;
  }
  egadsBody *pbody = (egadsBody *) dst->blind;
  
  if (src->oclass == BODY) {
  
    // use hashed body data on both ends
    egadsBody *pbods = (egadsBody *) src->blind;
    shape = pbody->shape;
    if (shape.IsSame(pbods->shape)) EG_attributeDup(src, dst);
    
    nents = pbods->shells.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->shells.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->shells.map(i+1);
      j     = pbody->shells.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->shells.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }

    nents = pbods->faces.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->faces.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->faces.map(i+1);
      j     = pbody->faces.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->faces.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }
    
    nents = pbods->loops.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->loops.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->loops.map(i+1);
      j     = pbody->loops.map.FindIndex(shape);
      if (j == 0) continue;                     // not in the dst body
      dobj = pbody->loops.objs[j-1];
      EG_attributeDup(aobj, dobj);
    }

    nents = pbods->edges.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->edges.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->edges.map(i+1);
      j     = pbody->edges.map.FindIndex(shape);
      if (j == 0) {
        if (fullAttr == 0) continue;
        if (aobj->mtype == DEGENERATE) continue;
        egadsEdge *pedge = (egadsEdge *) aobj->blind;
        trange[0] = pedge->trange[0];
        trange[1] = pedge->trange[1];
        geom      = pedge->curve;
        n = pbody->edges.map.Extent();
        for (j = 1; j <= n; j++) {
          dobj = pbody->edges.objs[j-1];
          if (dobj->mtype == DEGENERATE) continue;
          if (EG_isSame(aobj, dobj) == 0) {
            if (dobj->mtype == ONENODE) {
              if (aobj->mtype == ONENODE) EG_attributeDup(aobj, dobj);
              continue;
            }
            egadsEdge *pedgd = (egadsEdge *) dobj->blind;
            egadsNode *pnod0 = (egadsNode *) pedgd->nodes[0]->blind;
            egadsNode *pnod1 = (egadsNode *) pedgd->nodes[1]->blind;
            stat = EG_invEvaluate(geom, pnod0->xyz, &tmin, xyz);
            if (stat != EGADS_SUCCESS) continue;
            stat = EG_invEvaluate(geom, pnod1->xyz, &tmax, xyz);
            if (stat != EGADS_SUCCESS) continue;
            if ((tmin+UVTOL < trange[0]) || (tmax-UVTOL > trange[1])) continue;
            EG_attributeDup(aobj, dobj);
          }
        }
      } else {
        dobj = pbody->edges.objs[j-1];
        EG_attributeDup(aobj, dobj);
      }
    }
  
    nents = pbods->nodes.map.Extent();
    for (i = 0; i < nents; i++) {
      aobj = pbods->nodes.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = attrs->nattrs;
      if (nattr <= 0) continue;
      shape = pbods->nodes.map(i+1);
      j     = pbody->nodes.map.FindIndex(shape);
      if (j == 0) {
        if (fullAttr == 0) continue;
        n = pbody->nodes.map.Extent();
        for (j = 1; j <= n; j++) {
          dobj = pbody->nodes.objs[j-1];
          if (EG_isSame(aobj, dobj) == 0) {
            EG_attributeDup(aobj, dobj);
            break;
          }
        }
      } else {
        dobj = pbody->nodes.objs[j-1];
        EG_attributeDup(aobj, dobj);
      }
    }

  } else {
  
    // traverse the source to find objects with attributes
    if (((dst->mtype == SOLIDBODY) || (dst->mtype == SHEETBODY)) &&
         (src->oclass == SHELL)) {
      egadsShell *pshell = (egadsShell *) src->blind;
      TopExp_Explorer Exp(pbody->shape, TopAbs_SHELL);
      bool found = false;
      for (; Exp.More(); Exp.Next()) {
        if (pshell->shell.IsSame(Exp.Current())) {
          EG_attriBodyTrav(src, Exp.Current(), pbody);
          found = true;
          break;
        }
      }
      if (!found)
        printf(" EGADS Internal: Dropping Shell and sub-shape attributes!\n");
    } else {
      EG_attriBodyTrav(src, pbody->shape, pbody);
    }
  }
  
  return EGADS_SUCCESS;
}


int
EG_attriBodyCopy(const egObject *src, /*@null@*/ double *xform, egObject *dst)
{
  int      i, nents, nattr;
  egObject *aobj, *dobj;
  egAttrs  *attrs;
  
  if ((src == NULL) || (dst == NULL)) return EGADS_NULLOBJ;
  if  (src->magicnumber != MAGIC)     return EGADS_NOTOBJ;
  if  (src->oclass < NODE)            return EGADS_NOTTOPO;
  if  (src->blind == NULL)            return EGADS_NODATA;
  int outLevel = EG_outLevel(src);
  
  if (src->oclass == MODEL) {
    if (outLevel > 0)
      printf(" EGADS Error: src MODEL not supported (EG_attriBodyCopy)!\n");
    return EGADS_NOTMODEL;
  }
  if (dst->magicnumber != MAGIC) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not an EGO (EG_attriBodyCopy)!\n");
    return EGADS_NOTOBJ;
  }
  if (src->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: src not a BODY (EG_attriBodyCopy)!\n");
    return EGADS_NOTBODY;
  }
  if (dst->oclass != BODY) {
    if (outLevel > 0)
      printf(" EGADS Error: dst not a BODY (EG_attriBodyCopy)!\n");
    return EGADS_NOTBODY;
  }
  if (src->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL src pointer (EG_attriBodyCopy)!\n");
    return EGADS_NODATA;
  }
  if (dst->blind == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: NULL dst pointer (EG_attriBodyCopy)!\n");
    return EGADS_NODATA;
  }
  EG_attributeXDup(src, xform, dst);
  egadsBody *pbods = (egadsBody *) src->blind;
  egadsBody *pbody = (egadsBody *) dst->blind;
  
  nents = pbods->shells.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->shells.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->shells.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->faces.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->faces.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->faces.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }
    
  nents = pbods->loops.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->loops.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->loops.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }

  nents = pbods->edges.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->edges.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->edges.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }
  
  nents = pbods->nodes.map.Extent();
  for (i = 0; i < nents; i++) {
    aobj = pbods->nodes.objs[i];
    if (aobj->attrs == NULL) continue;
    attrs = (egAttrs *) aobj->attrs;
    nattr = attrs->nattrs;
    if (nattr <= 0) continue;
    dobj = pbody->nodes.objs[i];
    EG_attributeXDup(aobj, xform, dobj);
  }
  
  return EGADS_SUCCESS;
}


void
EG_readAttrs(egObject *obj, int nattr, FILE *fp)
{
  int     i, j, n, type, namlen, len, ival, nseq, *ivec = NULL;
  char    *name, cval, *string = NULL;
  double  rval, *rvec = NULL;
  egAttrs *attrs = NULL;
  egAttr  *attr  = NULL;
  
  attr = (egAttr *) EG_alloc(nattr*sizeof(egAttr));
  if (attr != NULL) {
    attrs = (egAttrs *) EG_alloc(sizeof(egAttrs));
    if (attrs == NULL) {
      EG_free(attr);
      attr = NULL;
    }
  }
  
  for (nseq = n = i = 0; i < nattr; i++) {
    j = fscanf(fp, "%d %d %d", &type, &namlen, &len);
    if (j != 3) break;
    name = NULL;
    if ((attrs != NULL) && (namlen != 0))
      name = (char *) EG_alloc((namlen+1)*sizeof(char));
    if (name != NULL) {
      fscanf(fp, "%s", name);
      for (j = 0; j < namlen; j++)
        if (name[j] == 127) {
          nseq++;
          name[j] = 32;
          break;
        }
    } else {
      for (j = 0; j < namlen; j++) fscanf(fp, "%c", &cval);
    }
    if (type == ATTRINT) {
      if (len == 1) {
        fscanf(fp, "%d", &ival);
      } else {
        ivec = NULL;
        if ((name != NULL) && (len != 0))
          ivec = (int *) EG_alloc(len*sizeof(int));
        if (ivec == NULL) {
          for (j = 0; j < len; j++) fscanf(fp, "%d", &ival);
          if (name != NULL) {
            EG_free(name);
            name = NULL;
          }
        } else {
          for (j = 0; j < len; j++) fscanf(fp, "%d", &ivec[j]);
        }
      }
    } else if ((type == ATTRREAL) || (type == ATTRCSYS)) {
      if (len == 1) {
        fscanf(fp, "%lf", &rval);
      } else {
        rvec = NULL;
        if ((name != NULL) && (len != 0))
          rvec = (double *) EG_alloc(len*sizeof(double));
        if (rvec == NULL) {
          for (j = 0; j < len; j++) fscanf(fp, "%lf", &rval);
          if (name != NULL) {
            EG_free(name);
            name = NULL;
          }
        } else {
          for (j = 0; j < len; j++) fscanf(fp, "%lf", &rvec[j]);
        }
      } 
    } else {
      do {
        j = fscanf(fp, "%c", &cval);
        if (j == 0) break;
      } while (cval != '#');
      string = NULL;
      if (name != NULL)
        string = (char *) EG_alloc((len+1)*sizeof(char));
      if (string != NULL) {
        for (j = 0; j < len; j++) fscanf(fp, "%c", &string[j]);
        string[len] = 0;
      } else {
        for (j = 0; j < len; j++) fscanf(fp, "%c", &cval);
        if (name != NULL) {
          EG_free(name);
          name = NULL;
        }
      }
    }

    if (name != NULL) {
      attr[n].name   = name;
      attr[n].type   = type;
      attr[n].length = len;
      if (type == ATTRINT) {
        if (len == 1) {
          attr[n].vals.integer  = ival;
        } else {
          attr[n].vals.integers = ivec;
        }
      } else if ((type == ATTRREAL) || (type == ATTRCSYS)) {
        if (len == 1) {
          attr[n].vals.real  = rval;
        } else {
          attr[n].vals.reals = rvec;
        }
      } else {
        attr[n].vals.string = string;
      }
      n++;
    }
  }
  
  if (attrs != NULL) {
    attrs->nattrs = n;
    attrs->attrs  = attr;
    attrs->nseqs  = 0;
    attrs->seqs   = NULL;
    if (nseq != 0) EG_attrBuildSeq(attrs);
    obj->attrs    = attrs;
  }
}


void
EG_initOCC()
{
#if CASVER >= 690
  int  np;
  char *env;

  env = getenv("EMPnumProc");
  np  = 2;
  if (env != NULL) np = atoi(env);
  if (np > 1) BOPAlgo_Options::SetParallelMode(Standard_True);
#endif
  if (((OCC_VERSION_MAJOR       != 7) || (OCC_VERSION_MINOR != 3) ||
       (OCC_VERSION_MAINTENANCE == 0)) &&
      ((OCC_VERSION_MAJOR       != 7) || (OCC_VERSION_MINOR != 4) ||
       (OCC_VERSION_MAINTENANCE == 0)))
    printf(" EGADS WARNING: OpenCASCADE %d.%d.%d NOT an Authorized Release!\n",
           OCC_VERSION_MAJOR, OCC_VERSION_MINOR, OCC_VERSION_MAINTENANCE);
}


static int
EG_readTess(FILE *fp, egObject *body, egObject **tess)
{
  int     i, j, ir, status, nnode, nedge, nface, n[3], len, ntri, nattr;
  int     *ptype, *pindex, *tris, *tric;
  double  *xyz, *param;
  egObject *obj;
  
  *tess  = NULL;
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }
  
  ir = fscanf(fp, "%d %d %d", &n[0], &n[1], &n[2]);
  if (ir != 3) {
    printf(" EGADS Error: Header with only %d words (EG_readTess)!\n", ir);
    return EGADS_INDEXERR;
  }
  if ((nnode != n[0]) || (nedge != n[1]) || (nface != n[2])) {
    printf(" EGADS Error: Count mismatch %d %d  %d %d  %d %d (EG_readTess)!\n",
           nnode, n[0], nedge, n[1], nface, n[2]);
    return EGADS_INDEXERR;
  }
  
  /* initialize the Tessellation Object */
  status = EG_initTessBody(body, tess);
  if (status != EGADS_SUCCESS) return status;
  EG_dereferenceTopObj(body, *tess);
  
  /* do the Edges */
  for (i = 0; i < nedge; i++) {
    len = 0;
    fscanf(fp, "%d", &len);
    if (len == 0) continue;
    xyz   = (double *) malloc(3*len*sizeof(double));
    param = (double *) malloc(  len*sizeof(double));
    if ((xyz == NULL) || (param == NULL)) {
      printf(" EGADS Error: malloc on Edge %d -- len = %d (EG_readTess)!\n",
               i+1, len);
      if (xyz   != NULL) free(xyz);
      if (param != NULL) free(param);
      EG_deleteObject(*tess);
      *tess = NULL;
      return EGADS_MALLOC;
    }
    for (j = 0; j < len; j++) {
      ir = fscanf(fp, "%le %le %le %le", &xyz[3*j  ], &xyz[3*j+1], &xyz[3*j+2],
                  &param[j]);
      if (ir != 4) {
        printf(" EGADS Error: %d/%d Read got %d out of 4 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    status = EG_setTessEdge(*tess, i+1, len, xyz, param);
    free(xyz);
    free(param);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessEdge %d = %d (EG_readTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      return status;
    }
  }
  
  /* do the Faces */
  for (i = 0; i < nface; i++) {
    len = ntri = 0;
    fscanf(fp, "%d %d", &len, &ntri);
    if ((len == 0) || (ntri == 0)) continue;
    xyz    = (double *) malloc(3*len*sizeof(double));
    param  = (double *) malloc(2*len*sizeof(double));
    ptype  = (int *)    malloc(  len* sizeof(int));
    pindex = (int *)    malloc(  len* sizeof(int));
    tris   = (int *)    malloc(3*ntri*sizeof(int));
    tric   = (int *)    malloc(3*ntri*sizeof(int));
    if ((xyz    == NULL) || (param == NULL) || (ptype == NULL) ||
        (pindex == NULL) || (tris  == NULL) || (tric  == NULL)) {
      printf(" EGADS Error: malloc on Face %d -- lens = %d %d (EG_readTess)!\n",
             i+1, len, ntri);
      if (xyz    != NULL) free(xyz);
      if (param  != NULL) free(param);
      if (ptype  != NULL) free(ptype);
      if (pindex != NULL) free(pindex);
      if (tris   != NULL) free(tris);
      if (tric   != NULL) free(tric);
      EG_deleteObject(*tess);
      *tess = NULL;
      return EGADS_MALLOC;
    }
    for (j = 0; j < len; j++) {
      ir = fscanf(fp, "%le %le %le %le %le %d %d", &xyz[3*j], &xyz[3*j+1],
                  &xyz[3*j+2], &param[2*j], &param[2*j+1], &ptype[j], &pindex[j]);
      if (ir != 7) {
        printf(" EGADS Error: %d/%d Read got %d out of 7 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        free(ptype);
        free(pindex);
        free(tris);
        free(tric);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    for (j = 0; j < ntri; j++) {
      ir = fscanf(fp, "%d %d %d %d %d %d", &tris[3*j], &tris[3*j+1],
                  &tris[3*j+2], &tric[3*j], &tric[3*j+1], &tric[3*j+2]);
      if (ir != 6) {
        printf(" EGADS Error: %d/%d Read got %d out of 6 (EG_readTess)!\n",
               j+1, len, ir);
        free(xyz);
        free(param);
        free(ptype);
        free(pindex);
        free(tris);
        free(tric);
        EG_deleteObject(*tess);
        *tess = NULL;
        return EGADS_READERR;
      }
    }
    status = EG_setTessFace(*tess, i+1, len, xyz, param, ntri, tris);
    free(xyz);
    free(param);
    free(ptype);
    free(pindex);
    free(tris);
    free(tric);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Warning: EG_setTessFace %d = %d (EG_readTess)!\n",
             i+1, status);
/*    EG_deleteObject(*tess);
      *tess = NULL;
      return status; */
    }
  }
  
  /* close up the open tessellation */
  status = EG_statusTessBody(*tess, &obj, &i, &len);
  if (status == EGADS_OUTSIDE) {
    printf(" EGADS Warning: Tessellation Object is incomplete (EG_readTess)!\n");
    egTessel *btess = (egTessel *) (*tess)->blind;
    btess->done = 1;
  } else if (status != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_statusTessBody = %d (EG_readTess)!\n", status);
    EG_deleteObject(*tess);
    *tess = NULL;
    return status;
  }
  if ((status != EGADS_OUTSIDE) && (i != 1)) {
    printf(" EGADS Warning: Tessellation Object is %d (EG_readTess)!\n", i);
/*  EG_deleteObject(*tess);
    *tess = NULL;
    return EGADS_TESSTATE;  */
    egTessel *btess = (egTessel *) (*tess)->blind;
    btess->done = 1;
  }
  
  /* attach the attributes */
  fscanf(fp, "%d\n", &nattr);
  if (nattr != 0) EG_readAttrs(*tess, nattr, fp);
  
  return EGADS_SUCCESS;
}


int
EG_loadModel(egObject *context, int bflg, const char *name, egObject **model)
{
  int          i, j, stat, outLevel, len, nattr, nerr, hite, hitf, nbs, egads;
  int          oclass, ibody;
  double       scale = 1.0;
  const char   *units;
  egObject     *omodel, *aobj;
  TopoDS_Shape source;
  egadsModel   *mshape = NULL;
  egadsLabel   *labels = NULL;
  FILE         *fp;
  
  *model = NULL;
  egads  = 0;
  nbs    = 0;
  if (context == NULL)               return EGADS_NULLOBJ;
  if (context->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (context->oclass != CONTXT)     return EGADS_NOTCNTX;
  if (EG_sameThread(context))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(context);
  
  if (name == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: NULL Filename (EG_loadModel)!\n");
    return EGADS_NONAME;
  }
  
  /* does file exist? */

  fp = fopen(name, "r");
  if (fp == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Not Found (EG_loadModel)!\n", name);
    return EGADS_NOTFOUND;
  }
  fclose(fp);
  
  /* find extension */
  
  len = strlen(name);
  for (i = len-1; i > 0; i--)
    if (name[i] == '.') break;
  if (i == 0) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Extension in %s (EG_loadModel)!\n", name);
    return EGADS_NODATA;
  }
  
  if ((strcasecmp(&name[i],".step") == 0) || 
      (strcasecmp(&name[i],".stp") == 0)) {

    /* STEP files */

    STEPControl_Reader aReader;
    IFSelect_ReturnStatus status = aReader.ReadFile(name);
    if (status != IFSelect_RetDone) {
      if (outLevel > 0)
        printf(" EGADS Error: STEP Read of %s = %d (EG_loadModel)!\n", 
               name, status);
      return EGADS_NOLOAD;
    }

    // inspect the root transfers
    if (outLevel > 2)
      aReader.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);

    if ((bflg&16) != 0) egads = -1;
    if ((bflg&4)  == 0) {
      TColStd_SequenceOfAsciiString unitLength, unitAngle, solidAngle;
      aReader.FileUnits(unitLength, unitAngle, solidAngle);
      if (unitLength.Length() >= 1) {
        if (unitLength.Length() > 1)
          printf(" EGADS Info: # unitLengths = %d\n", unitLength.Length());
        units = unitLength(1).ToCString();
        if ((strcasecmp(units,"inch")   == 0) ||
            (strcasecmp(units,"inches") == 0)) {
          scale = 1.0/25.4;
          printf("  STEP Reader Info: Using %s\n", units);
        } else if ((strcasecmp(units,"foot") == 0) ||
                   (strcasecmp(units,"feet") == 0)) {
          scale = 1.0/304.8;
          printf("  STEP Reader Info: Using %s\n", units);
        } else if ((strcasecmp(units,"metre")  == 0) ||
                   (strcasecmp(units,"meter")  == 0) ||
                   (strcasecmp(units,"meters") == 0)) {
          scale = 1.0/1000.0;
          printf("  STEP Reader Info: Using %s\n", units);
        } else if ((strcasecmp(units,"centimetre")  == 0) ||
                   (strcasecmp(units,"centimeter")  == 0) ||
                   (strcasecmp(units,"centimeters") == 0)) {
          scale = 1.0/100.0;
          printf("  STEP Reader Info: Using %s\n", units);
        } else if ((strcasecmp(units,"millimetre")  == 0) ||
                   (strcasecmp(units,"millimeter")  == 0) ||
                   (strcasecmp(units,"millimeters") == 0)) {
          printf("  STEP Reader Info: Using %s\n", units);
        } else {
          printf(" EGADS STEP Info: Cannot convert %s -- using millimeters!\n",
                 units);
        }
      }
    }

    int nroot = aReader.NbRootsForTransfer();
    if (outLevel > 1)
      printf(" EGADS Info: %s Entries = %d\n", name, nroot);

    for (i = 1; i <= nroot; i++) {
      Standard_Boolean ok = aReader.TransferRoot(i);
      if ((!ok) && (outLevel > 0))
        printf(" EGADS Warning: Transfer %d/%d is not OK!\n", i, nroot);
    }

    nbs = aReader.NbShapes();
    if (nbs <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: %s has No Shapes (EG_loadModel)!\n", 
               name);
      return EGADS_NOLOAD;
    }
    if (outLevel > 1)    
      printf(" EGADS Info: %s has %d Shape(s)\n", name, nbs);

#ifdef STEPATTRS
    const Handle(XSControl_WorkSession)& workSession = aReader.WS();
    const Handle(XSControl_TransferReader)& TR = workSession->TransferReader();
    Handle(Standard_Type) tNAUO = STANDARD_TYPE(
                                       StepRepr_NextAssemblyUsageOccurrence);
    Handle(Standard_Type) tPD   = STANDARD_TYPE(StepBasic_ProductDefinition);
#endif

    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 1; i <= nbs; i++) {
      TopoDS_Shape aShape = aReader.Shape(i);
#ifdef STEPATTRS
      Handle(Standard_Transient) ent = TR->EntityFromShapeResult(aShape, 3);
      if (!ent.IsNull()) {
        if (ent->DynamicType() == tNAUO) {
          Handle(StepRepr_NextAssemblyUsageOccurrence) NAUO =
              Handle(StepRepr_NextAssemblyUsageOccurrence)::DownCast(ent);
          if (!NAUO.IsNull()) {
            Interface_EntityIterator subs = workSession->Graph().Sharings(NAUO);
            for (subs.Start(); subs.More(); subs.Next()) {
              Handle(StepRepr_ProductDefinitionShape) PDS =
                Handle(StepRepr_ProductDefinitionShape)::DownCast(subs.Value());
              if (!PDS.IsNull()) {
                Handle(StepBasic_ProductDefinitionRelationship) PDR =
                              PDS->Definition().ProductDefinitionRelationship();
                if (!PDR.IsNull()) {
                  if (PDR->HasDescription() && PDR->Description()->Length() > 0) {
                    printf(" NAUOa = %s\n", PDR->Description()->ToCString());
                  } else if (PDR->Name()->Length() > 0) {
                    printf(" NAUOb = %s\n", PDR->Name()->ToCString());
                  } else {
                    printf(" NAUOc = %s\n", PDR->Id()->ToCString());
                  }
                }
              }
            }
          }
        } else if (ent->DynamicType() == tPD) {
          Handle(StepBasic_ProductDefinition) PD =
              Handle(StepBasic_ProductDefinition)::DownCast(ent);
          if (!PD.IsNull()) {
            Handle(StepBasic_Product) Prod = PD->Formation()->OfProduct();
            if (Prod->Name()->UsefullLength() > 0) {
              printf(" PDa   = %s\n", Prod->Name()->ToCString());
            } else {
              printf(" PDb   = %s\n", Prod->Id()->ToCString());
            }
          }
        }
      }
#endif
      if ((bflg&8) != 0) {
        ShapeUpgrade_UnifySameDomain uShape(aShape, Standard_True,
                                            Standard_True, Standard_True);
//      uShape.SetLinearTolerance(100.0*Precision::Confusion());
//      uShape.SetAngularTolerance(10.0*Precision::Angular());
        uShape.Build();
        aShape = uShape.Shape();
      }
      if (scale != 1.0) {
        gp_Trsf form = gp_Trsf();
        form.SetValues(scale, 0.0,   0.0,   0.0,
                       0.0,   scale, 0.0,   0.0,
                       0.0,   0.0,   scale, 0.0);
        BRepBuilderAPI_Transform xForm(aShape, form, Standard_True);
        if (!xForm.IsDone()) {
          printf(" EGADS Warning: Can't scale Body %d (EG_loadModel)!\n", i);
        } else {
          aShape = xForm.ModifiedShape(aShape);
        }
      }
      builder3D.Add(compound, aShape);
    }
    source = compound;
    
    if (outLevel > 0) {
      hite = hitf = 0;
      TopTools_IndexedMapOfShape MapE, MapF;
      TopExp::MapShapes(source, TopAbs_EDGE, MapE);
      for (i = 1; i <= MapE.Extent(); i++) {
        Standard_Real t1, t2;
        TopoDS_Shape shape = MapE(i);
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        Handle(Geom_BSplineCurve) hBSpline =
          Handle(Geom_BSplineCurve)::DownCast(hCurve);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hite++;
      }
      TopExp::MapShapes(source, TopAbs_FACE, MapF);
      for (int i = 1; i <= MapF.Extent(); i++) {
        TopoDS_Shape shape = MapF(i);
        TopoDS_Face  Face  = TopoDS::Face(shape);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
        Handle(Geom_BSplineSurface) hBSpline =
          Handle(Geom_BSplineSurface)::DownCast(hSurface);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
      }
      if (hite+hitf != 0)
        printf(" EGADS Info: Import Model has %d C0 Edges & %d C0 Faces!\n",
               hite, hitf);
    }
    
  } else if ((strcasecmp(&name[i],".iges") == 0) || 
             (strcasecmp(&name[i],".igs") == 0)) {
             
    /* IGES files */
    
    IGESControl_Reader iReader;
    Standard_Integer stats = iReader.ReadFile(name);
    if (stats != IFSelect_RetDone) {
      if (outLevel > 0)
        printf(" EGADS Error: IGES Read of %s = %d (EG_loadModel)!\n", 
               name, stats);
      return EGADS_NOLOAD;
    }
    iReader.TransferRoots();
    if ((bflg&16) != 0) egads = -1;

    nbs = iReader.NbShapes();
    if (nbs <= 0) {
      if (outLevel > 0)
        printf(" EGADS Error: %s has No Shapes (EG_loadModel)!\n", 
               name);
      return EGADS_NOLOAD;
    }
    if (outLevel > 1)    
      printf(" EGADS Info: %s has %d Shape(s)\n", name, nbs);
    
    const Handle(XSControl_WorkSession)& workSession = iReader.WS();
    const Handle(XSControl_TransferReader)& transferReader =
          workSession->TransferReader();
    labels = new egadsLabel[nbs];

    TopoDS_Compound compound;
    BRep_Builder    builder3D;
    builder3D.MakeCompound(compound);
    for (i = 1; i <= nbs; i++) {
      TopoDS_Shape aShape = iReader.Shape(i);
      
      labels[i-1].shapeName = NULL;
      const Handle(IGESData_IGESEntity)& shapeEntity =
            Handle(IGESData_IGESEntity)::DownCast(
                  transferReader->EntityFromShapeResult(aShape, -1));
      if (shapeEntity->HasName())
        labels[i-1].shapeName = EG_strdup(shapeEntity->NameValue()->ToCString());
      
      if ((bflg&8) != 0) {
        ShapeUpgrade_UnifySameDomain uShape(aShape);
        uShape.Build();
        aShape = uShape.Shape();
      }
      labels[i-1].shape = aShape;
      builder3D.Add(compound, aShape);
    }
    source = compound;
    
    if (outLevel > 0) {
      hite = hitf = 0;
      TopTools_IndexedMapOfShape MapE, MapF;
      TopExp::MapShapes(source, TopAbs_EDGE, MapE);
      for (i = 1; i <= MapE.Extent(); i++) {
        Standard_Real t1, t2;
        TopoDS_Shape shape = MapE(i);
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        Handle(Geom_BSplineCurve) hBSpline =
          Handle(Geom_BSplineCurve)::DownCast(hCurve);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hite++;
      }
      TopExp::MapShapes(source, TopAbs_FACE, MapF);
      for (int i = 1; i <= MapF.Extent(); i++) {
        TopoDS_Shape shape = MapF(i);
        TopoDS_Face  Face  = TopoDS::Face(shape);
        Handle(Geom_Surface) hSurface = BRep_Tool::Surface(Face);
        Handle(Geom_BSplineSurface) hBSpline =
          Handle(Geom_BSplineSurface)::DownCast(hSurface);
        if (hBSpline.IsNull()) continue;
        if (hBSpline->Continuity() == GeomAbs_C0)  hitf++;
      }
      if (hite+hitf != 0)
        printf(" EGADS Info: Import Model has %d C0 Edges & %d C0 Faces!\n",
               hite, hitf);
    }

  } else if ((strcasecmp(&name[i],".brep") == 0) ||
             (strcasecmp(&name[i],".egads") == 0)) {
  
    /* Native OCC file */
    if (strcasecmp(&name[i],".egads") == 0) egads = 1;

    BRep_Builder builder;
    if (!BRepTools::Read(source, name, builder)) {
      if (outLevel > 0)
        printf(" EGADS Warning: Read Error on %s (EG_loadModel)!\n", name);
      return EGADS_NOLOAD;
    }

  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Extension in %s Not Supported (EG_loadModel)!\n",
             name);
    return EGADS_NODATA;
  }
  
  int nWire  = 0;
  int nFace  = 0;
  int nSheet = 0;
  int nSolid = 0;
  
  TopExp_Explorer Exp;
  for (Exp.Init(source, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) nWire++;
  for (Exp.Init(source, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) nFace++;
  for (Exp.Init(source, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) nSheet++;
  for (Exp.Init(source, TopAbs_SOLID); Exp.More(); Exp.Next()) nSolid++;

  if (outLevel > 1)
    printf("\n EGADS Info: %s has %d Solids, %d Sheets, %d Faces and %d Wires\n",
           name, nSolid, nSheet, nFace, nWire);
           
  int nBody = nWire+nFace+nSheet+nSolid;
  /* can we promote Edges to Wires? */
  if ((nBody == 0) || (egads == -1)) {
    double tol = 1.e-7;
/*
    for (Exp.Init(source, TopAbs_VERTEX); Exp.More(); Exp.Next()) {
      TopoDS_Vertex Vert = TopoDS::Vertex(Exp.Current());
      double tolv = BRep_Tool::Tolerance(Vert);
      if (tolv > tol) tol = tolv;
    }
 */
    int nEdge, add = 0;
    j = 0;
    for (Exp.Init(source, TopAbs_EDGE, TopAbs_WIRE); Exp.More(); Exp.Next()) j++;
    if (j != 0) {
      nEdge = j;
      TopoDS_Edge *Edges = new TopoDS_Edge[nEdge];
      /* remove small Edges */
      j = 0;
      for (Exp.Init(source, TopAbs_EDGE, TopAbs_WIRE); Exp.More(); Exp.Next()) {
        TopoDS_Shape shape = Exp.Current();
        TopoDS_Edge  Edge  = TopoDS::Edge(shape);
        Standard_Real t1, t2;
        Handle(Geom_Curve) hCurve = BRep_Tool::Curve(Edge, t1, t2);
        GeomAdaptor_Curve AC(hCurve);
        if (GCPnts_AbscissaPoint::Length(AC, t1, t2) < tol) continue;
        Edges[j] = Edge;
        j++;
      }
      nEdge = j;
      TopoDS_Compound compound;
      BRep_Builder    builder3D;
      builder3D.MakeCompound(compound);
      if (nBody != 0) builder3D.Add(compound, source);
      while (j != 0) {
        for (i = 0; i < nEdge; i++)
          if (!Edges[i].IsNull()) {
            BRepBuilderAPI_MakeWire MW;
            try {
              MW.Add(Edges[i]);
              TopoDS_Vertex V1, V2;
              TopExp::Vertices(Edges[i], V2, V1, Standard_True);
              if (outLevel > 1) printf(" start Edge %d:", i+1);
              Edges[i].Nullify();
              if (V2.IsSame(V1)) {
                if (outLevel > 1) printf("\n");
                break;
              }
              gp_Pnt pv1 = BRep_Tool::Pnt(V1);
              gp_Pnt pv2 = BRep_Tool::Pnt(V2);
              int hit;
              do {
                hit = 0;
                for (int k = 0; k < nEdge; k++)
                  if (!Edges[k].IsNull()) {
                    TopExp::Vertices(Edges[k], V2, V1, Standard_True);
                    gp_Pnt pv = BRep_Tool::Pnt(V1);
                    if (pv.Distance(pv1) < tol) {
                      if (outLevel > 1) printf(" --%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv1 = pv;
                      hit++;
                      break;
                    }
                    if (pv.Distance(pv2) < tol) {
                      if (outLevel > 1) printf(" -+%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv2 = pv;
                      hit++;
                      break;
                    }
                    pv = BRep_Tool::Pnt(V2);
                    if (pv.Distance(pv1) < tol) {
                      if (outLevel > 1) printf(" +-%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv1 = pv;
                      hit++;
                      break;
                    }
                    if (pv.Distance(pv2) < tol) {
                      if (outLevel > 1) printf(" ++%d", k+1);
                      MW.Add(Edges[k]);
                      Edges[k].Nullify();
                      pv2 = pv;
                      hit++;
                      break;
                    }
                  }
              } while (hit != 0);
              if (outLevel > 1) printf("\n");
            }
            catch (const Standard_Failure& e) {
              printf(" EGADS Warning: Cannot Add Edge in Wire!\n");
              printf("                %s\n", e.GetMessageString());
              continue;
            }
            catch (...) {
              printf(" EGADS Warning: Cannot Add Edge in Wire!\n");
              continue;
            }
            if (MW.Error()) {
              if (outLevel > 0)
                printf(" EGADS Error: Problem with Imported Edge!\n");
              continue;
            }
            if (!MW.IsDone()) {
              if (outLevel > 0)
                printf(" EGADS Error: Problem with Loop Conversion!\n");
              continue;
            }
            TopoDS_Wire wire = MW.Wire();
            builder3D.Add(compound, wire);
            nWire++;
            nBody++;
            add++;
          }
        j = 0;
        for (i = 0; i < nEdge; i++) if (!Edges[i].IsNull()) j++;
      }
      source = compound;
      delete [] Edges;
    }
    if (add != 0)
      if (outLevel > 0)
        printf(" EGADS Info: %d Edges converted to %d WireBodies (EG_loadModel)!\n",
               nEdge, add);
  }
  if (nBody == 0) {
    source.Nullify();
    if (outLevel > 0)
      printf(" EGADS Warning: Nothing found in %s (EG_loadModel)!\n", name);
    return EGADS_NODATA;
  }
  
  mshape              = new egadsModel;
  mshape->shape       = source;
  mshape->nobjs       = nBody;
  mshape->nbody       = nBody;
  mshape->bbox.filled = 0;
  mshape->bodies = new egObject*[nBody];
  for (i = 0; i < nBody; i++) {
    stat = EG_makeObject(context, &mshape->bodies[i]);
    if (stat != EGADS_SUCCESS) {
      for (int j = 0; j < i; j++) {
        egObject  *obj   = mshape->bodies[j];
        egadsBody *pbody = (egadsBody *) obj->blind;
        delete pbody;
        EG_deleteObject(mshape->bodies[j]);
      }
      delete [] mshape->bodies;
      delete mshape;
      return stat;
    }
    egObject  *pobj    = mshape->bodies[i];
    egadsBody *pbody   = new egadsBody;
    pbody->nodes.objs  = NULL;
    pbody->edges.objs  = NULL;
    pbody->loops.objs  = NULL;
    pbody->faces.objs  = NULL;
    pbody->shells.objs = NULL;
    pbody->senses      = NULL;
    pbody->bbox.filled = 0;
    pbody->massFill    = 0;
    pobj->blind        = pbody;
  }
  
  i = 0;
  for (Exp.Init(mshape->shape, TopAbs_WIRE,  TopAbs_FACE);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer wCheck(pbody->shape);
      if (!wCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid WireBody!\n", i);
        } else {
          BRepCheck_Analyzer wfCheck(fixedShape);
          if (!wfCheck.IsValid()) {
            if (outLevel > 0)
              printf(" EGADS Warning: Fixed Body %d is an Invalid WireBody!\n",
                     i);
          } else {
            pbody->shape = fixedShape;
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_FACE,  TopAbs_SHELL);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer fCheck(pbody->shape);
      if (!fCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid FaceBody!\n", i);
        } else {
          BRepCheck_Analyzer ffCheck(fixedShape);
          if (!ffCheck.IsValid()) {
            if (outLevel > 0)
              printf(" EGADS Warning: Fixed Body %d is an Invalid FaceBody!\n",
                     i);
          } else {
            pbody->shape = fixedShape;
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_SHELL, TopAbs_SOLID);
       Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer sCheck(pbody->shape);
      if (!sCheck.IsValid()) {
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid SheetBody!\n", i);
        } else {
          if (fixedShape.ShapeType() != TopAbs_SHELL) {
            if (outLevel > 0)
              printf(" EGADS Reject: Fixed Body %d is No longer a SheetBody!\n",
                     i);
          } else {
            BRepCheck_Analyzer sfCheck(fixedShape);
            if (!sfCheck.IsValid()) {
              if (outLevel > 0)
                printf(" EGADS Warning: Fixed Body %d is an Invalid SheetBody!\n",
                       i);
            } else {
              pbody->shape = fixedShape;
            }
          }
        }
      }
    }
  }
  for (Exp.Init(mshape->shape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
    egObject  *obj   = mshape->bodies[i++];
    egadsBody *pbody = (egadsBody *) obj->blind;
    pbody->shape     = Exp.Current();
    if (egads == 0) {
      BRepCheck_Analyzer sCheck(pbody->shape);
      if (!sCheck.IsValid()) {
/*      ShapeFix_ShapeTolerance sTol;
        sTol.SetTolerance(pbody->shape, 1.e-4, TopAbs_SHELL);  */
        Handle_ShapeFix_Shape sfs = new ShapeFix_Shape(pbody->shape);
        sfs->Perform();
        TopoDS_Shape fixedShape = sfs->Shape();
        if (fixedShape.IsNull()) {
          if (outLevel > 0)
            printf(" EGADS Warning: Body %d is an Invalid SolidBody!\n", i);
        } else {
          if (fixedShape.ShapeType() != TopAbs_SHELL) {
            if (outLevel > 0)
              printf(" EGADS Reject: Fixed Body %d is No longer a SolidBody!\n",
                     i);
          } else {
            BRepCheck_Analyzer sfCheck(fixedShape);
            if (!sfCheck.IsValid()) {
              if (outLevel > 0)
                printf(" EGADS Warning: Fixed Body %d is an Invalid SolidBody!\n",
                       i);
            } else {
              pbody->shape = fixedShape;
            }
          }
        }
      }
    }
  }
  
  stat = EG_makeObject(context, &omodel);
  if (stat != EGADS_SUCCESS) {
    source.Nullify();
    for (i = 0; i < nBody; i++) {
      egObject  *obj   = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) obj->blind;
      delete pbody;
      EG_deleteObject(mshape->bodies[i]);
    }
    delete [] mshape->bodies;
    delete mshape;
    return stat;
  }
  omodel->oclass = MODEL;
  omodel->blind  = mshape;
  EG_referenceObject(omodel, context);
  
  for (j = i = 0; i < nBody; i++) {
    egObject  *pobj  = mshape->bodies[i];
    egadsBody *pbody = (egadsBody *) pobj->blind;
    pobj->topObj     = omodel;
    if (((bflg&1) == 0) && (egads == 0)) EG_splitPeriodics(pbody);
    if (((bflg&2) != 0) && (egads == 0)) EG_splitMultiplicity(pbody, outLevel);
    stat = EG_traverseBody(context, i, pobj, omodel, pbody, &nerr);
    if (stat != EGADS_SUCCESS) {
      if (egads == 1) {
        mshape->nbody = i;
        EG_destroyTopology(omodel);
        return stat;
      }
      if (outLevel > 0)
        printf(" EGADS Warning: Body  %d parse = %d!\n", i+1, stat);
      continue;
    }
    mshape->bodies[j] = mshape->bodies[i];
    j++;
  }
  if (j != mshape->nbody) {
    if (outLevel > 0)
      printf(" EGADS Inform:  %d Bodies not included!\n\n", mshape->nbody-j);
    mshape->nbody = j;
    if (j == 0) {
      EG_destroyTopology(omodel);
      return EGADS_TOPOERR;
    }
  }
  *model = omodel;
  
  /* possibly assign attributes from IGES/STEP read */
  if (labels != NULL) {
    for (i = 0; i < nbs; i++) {
      if (labels[i].shapeName == NULL) continue;
      char *value = labels[i].shapeName;
/*    printf(" Shape %2d: %s\n", i+1, value);  */
      for (int ibody = 0; ibody < mshape->nbody; ibody++) {
        egObject  *pobj   = mshape->bodies[ibody];
        egadsBody *pbody  = (egadsBody *) pobj->blind;
        if (pbody->shape == labels[i].shape)
          EG_attributeAdd(pobj, "Name", ATTRSTRING, 1, NULL, NULL, value);
        j = pbody->nodes.map.FindIndex(labels[i].shape);
        if (j != 0) EG_attributeAdd(pbody->nodes.objs[j-1],  "Name", ATTRSTRING,
                                    1, NULL, NULL, value);
        j = pbody->edges.map.FindIndex(labels[i].shape);
        if (j != 0) EG_attributeAdd(pbody->edges.objs[j-1],  "Name", ATTRSTRING,
                                    1, NULL, NULL, value);
        j = pbody->loops.map.FindIndex(labels[i].shape);
        if (j != 0) EG_attributeAdd(pbody->loops.objs[j-1],  "Name", ATTRSTRING,
                                    1, NULL, NULL, value);
        j = pbody->faces.map.FindIndex(labels[i].shape);
        if (j != 0) EG_attributeAdd(pbody->faces.objs[j-1],  "Name", ATTRSTRING,
                                    1, NULL, NULL, value);
        j = pbody->shells.map.FindIndex(labels[i].shape);
        if (j != 0) EG_attributeAdd(pbody->shells.objs[j-1], "Name", ATTRSTRING,
                                    1, NULL, NULL, value);
      }
      EG_free(value);
    }
    delete [] labels;
  }
  if (egads != 1) return EGADS_SUCCESS;

  /* get the attributes from the EGADS files */
  
  fp = fopen(name, "r");
  if (fp == NULL) {
    printf(" EGADS Info: Cannot reOpen %s (EG_loadModel)!\n", name);
    return EGADS_SUCCESS;
  }
  char line[81];
  for (;;) {
    line[0] = line[1] = ' ';
    if (fgets(line, 81, fp) == NULL) break;
    if ((line[0] == '#') && (line[1] == '#')) break;
  }
  
  // got the header
  if ((line[0] == '#') && (line[1] == '#')) {
    if (outLevel > 1) printf(" Header = %s\n", line);
    // get number of model attributes
    fscanf(fp, "%d", &nattr);
    if (nattr != 0) EG_readAttrs(omodel, nattr, fp);
    for (i = 0; i < nBody; i++) {
      int otype,  oindex;
      int rsolid, rshell, rface, rloop, redge, rnode;
      int nsolid, nshell, nface, nloop, nedge, nnode;

      fscanf(fp, " %d %d %d %d %d %d %d", &rsolid, &rshell, 
             &rface, &rloop, &redge, &rnode, &nattr);
      if (outLevel > 2)
        printf(" read = %d %d %d %d %d %d %d\n", rsolid, rshell, 
               rface, rloop, redge, rnode, nattr);
      egObject  *pobj  = mshape->bodies[i];
      egadsBody *pbody = (egadsBody *) pobj->blind;
      nnode  = pbody->nodes.map.Extent();
      nedge  = pbody->edges.map.Extent();
      nloop  = pbody->loops.map.Extent();
      nface  = pbody->faces.map.Extent();
      nshell = pbody->shells.map.Extent();
      nsolid = 0;
      if (pobj->mtype == SOLIDBODY) nsolid = 1;
      if ((nnode != rnode) || (nedge  != redge)  || (nloop  != rloop) ||
          (nface != rface) || (nshell != rshell) || (nsolid != rsolid)) {
        printf(" EGADS Info: %d %d, %d %d, %d %d, %d %d, %d %d, %d %d",
               nnode, rnode, nedge,  redge,  nloop,  rloop,
               nface, rface, nshell, rshell, nsolid, rsolid);
        printf("  MisMatch on Attributes (EG_loadModel)!\n");
        fclose(fp);
        return EGADS_SUCCESS;
      }
      // got the correct body -- transfer the attributes
      if (nattr != 0) EG_readAttrs(pobj, nattr, fp);
      for (;;)  {
        j = fscanf(fp, "%d %d %d\n", &otype, &oindex, &nattr);
        if (outLevel > 2)
          printf(" %d:  attr header = %d %d %d\n", 
                 j, otype, oindex, nattr);
        if (j     != 3) break;
        if (otype == 0) break;
        if (otype == 1) {
          aobj = pbody->shells.objs[oindex];
        } else if (otype == 2) {
          aobj = pbody->faces.objs[oindex];
        } else if (otype == 3) {
          aobj = pbody->loops.objs[oindex];
        } else if (otype == 4) {
          aobj = pbody->edges.objs[oindex];
        } else {
          aobj = pbody->nodes.objs[oindex];
        }
        EG_readAttrs(aobj, nattr, fp);
      }
    }
    
    /* get the ancillary objects from the EGADS files */
    line[0] = line[1] = ' ';
    fgets(line, 81, fp);
    if ((line[0] == '#') && (line[1] == '#')) {
      j = fscanf(fp, "%hd", &omodel->mtype);
      if ((j != 1) || (omodel->mtype < mshape->nbody)) {
        printf(" EGADS Info: Ext failure in %s  %d %d (EG_loadModel)!\n",
               name, j, omodel->mtype);
        omodel->mtype = 0;
        fclose(fp);
        return EGADS_SUCCESS;
      }
      egObject** bodies = new egObject*[omodel->mtype];
      for (j = 0; j < mshape->nbody; j++) bodies[j] = mshape->bodies[j];
      for (j = mshape->nbody; j < omodel->mtype; j++) bodies[j] = NULL;
      delete [] mshape->bodies;
      mshape->nobjs  = omodel->mtype;
      mshape->bodies = bodies;
      for (j = mshape->nbody; j < omodel->mtype; j++) {
        i = fscanf(fp, "%d %d", &oclass, &ibody);
        if (i != 2) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d (EG_loadModel)!\n",
                 name, i);
          break;
        }
        if (ibody > j) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d %d (EG_loadModel)!\n",
                 name, ibody, j);
          break;
        }
        if (oclass == TESSELLATION) {
          stat = EG_readTess(fp,  mshape->bodies[ibody-1], &mshape->bodies[j]);
          EG_referenceObject(mshape->bodies[ibody-1], mshape->bodies[j]);
        } else {
          stat = EG_readEBody(fp, mshape->bodies[ibody-1], &mshape->bodies[j]);
        }
        if (stat != EGADS_SUCCESS) {
          omodel->mtype = j;
          mshape->nobjs = omodel->mtype;
          printf(" EGADS Info: Ext read failure in %s  %d %d %d (EG_loadModel)!\n",
                 name, oclass, ibody, stat);
          break;
        }
        EG_referenceObject(mshape->bodies[j], omodel);
        EG_removeCntxtRef(mshape->bodies[j]);
        mshape->bodies[j]->topObj = omodel;
      }
    }
  } else {
    printf(" EGADS Info: EGADS Header not found in %s (EG_loadModel)!\n", 
           name);
  }
  
  fclose(fp);
  return EGADS_SUCCESS;
}


void
EG_writeAttr(egAttrs *attrs, FILE *fp)
{
  char c;
  int  namln;
  
  int    nattr = attrs->nattrs;
  egAttr *attr = attrs->attrs;
  for (int i = 0; i < nattr; i++) {
    if (attr[i].type == ATTRPTR) continue;
    namln = 0;
    if (attr[i].name != NULL) namln = strlen(attr[i].name);
    fprintf(fp, "%d %d %d\n", attr[i].type, namln, attr[i].length);
    if (namln != 0) {
      for (int j = 0; j < namln; j++) {
        c = attr[i].name[j];
        if (c == 32) c = 127;
        fprintf(fp, "%c", c);
      }
      fprintf(fp, "\n");
    }
    if (attr[i].type == ATTRINT) {
      if (attr[i].length == 1) {
        fprintf(fp, "%d\n", attr[i].vals.integer);
      } else {
        for (int j = 0; j < attr[i].length; j++)
          fprintf(fp, "%d ", attr[i].vals.integers[j]);
        fprintf(fp, "\n");
      }
    } else if ((attr[i].type == ATTRREAL) || (attr[i].type == ATTRCSYS)) {
      if (attr[i].length == 1) {
        fprintf(fp, "%19.12le\n", attr[i].vals.real);
      } else {
        for (int j = 0; j < attr[i].length; j++)
          fprintf(fp, "%19.12le ", attr[i].vals.reals[j]);
        fprintf(fp, "\n");
      }    
    } else if (attr[i].type == ATTRSTRING) {
      if (attr[i].length != 0) {
        fprintf(fp, "#%s\n", attr[i].vals.string);
      } else {
        /* we should never get here! */
        fprintf(fp, "#\n");
      }
    }
  }
}


int
EG_writeNumAttr(egAttrs *attrs)
{
  int num = 0;
  for (int i = 0; i < attrs->nattrs; i++) {
    if (attrs->attrs[i].type == ATTRPTR) continue;
    num++;
  }
  
  return num;
}


static void
EG_writeAttrs(const egObject *obj, FILE *fp)
{
  int     i, nsolid, nshell, nface, nloop, nedge, nnode, nattr = 0;
  egAttrs *attrs;
  
  attrs = (egAttrs *) obj->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;
  
  if (obj->oclass == MODEL) {
  
    fprintf(fp, "%d\n", nattr);
    if (nattr != 0) EG_writeAttr(attrs, fp);
    
  } else {

    egadsBody *pbody = (egadsBody *) obj->blind;
    nnode  = pbody->nodes.map.Extent();
    nedge  = pbody->edges.map.Extent();
    nloop  = pbody->loops.map.Extent();
    nface  = pbody->faces.map.Extent();
    nshell = pbody->shells.map.Extent();
    nsolid = 0;
    if (obj->mtype == SOLIDBODY) nsolid = 1;
    fprintf(fp, "  %d  %d  %d  %d  %d  %d  %d\n", nsolid, nshell, nface,  
            nloop, nedge, nnode, nattr);
    if (nattr != 0) EG_writeAttr(attrs, fp);
    
    for (i = 0; i < nshell; i++) {
      egObject *aobj = pbody->shells.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    1 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }

    for (i = 0; i < nface; i++) {
      egObject *aobj = pbody->faces.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    2 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }
    
    for (i = 0; i < nloop; i++) {
      egObject *aobj = pbody->loops.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    3 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }
        
    for (i = 0; i < nedge; i++) {
      egObject *aobj = pbody->edges.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    4 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }
        
    for (int i = 0; i < nnode; i++) {
      egObject *aobj = pbody->nodes.objs[i];
      if (aobj->attrs == NULL) continue;
      attrs = (egAttrs *) aobj->attrs;
      nattr = EG_writeNumAttr(attrs);
      if (nattr <= 0) continue;
      fprintf(fp, "    5 %d %d\n", i, nattr);
      EG_writeAttr(attrs, fp);
    }
    fprintf(fp, "    0 0 0\n");

  }

}


static int
EG_writeTess(const egObject *tess, FILE *fp)
{
  int          j, status, len;
  int          ntri, iedge, iface, nnode, nedge, nface, nattr = 0;
  const double *pxyz  = NULL, *puv    = NULL, *pt    = NULL;
  const int    *ptype = NULL, *pindex = NULL, *ptris = NULL, *ptric = NULL;
  egAttrs      *attrs;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  if (tess->blind == NULL)          return EGADS_NODATA;

  // get the body from tessellation
  egTessel  *btess = (egTessel *) tess->blind;
  egObject  *body  = btess->src;

  // get the sizes
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  if (body->oclass == EBODY) {
    status = EG_getBodyTopos(body, NULL, EEDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, EFACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  } else {
    status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
    if (status != EGADS_SUCCESS) return status;
    status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
    if (status != EGADS_SUCCESS) return status;
  }

  // write the number of nodes, edges, and faces
  fprintf(fp, " %d %d %d\n", nnode, nedge, nface);

  // write out the edge tessellation
  for (iedge = 0; iedge < nedge; iedge++) {
    status = EG_getTessEdge(tess, iedge+1, &len, &pxyz, &pt);
    if (status != EGADS_SUCCESS) return status;
    fprintf(fp, " %d\n", len);
    if (len == 0) continue;
    for (j = 0; j < len; j++)
      fprintf(fp, "%19.12le %19.12le %19.12le %19.12le\n", pxyz[3*j],
              pxyz[3*j+1], pxyz[3*j+2], pt[j]);
  }

  // write out face tessellations
  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess, iface+1, &len, &pxyz, &puv, &ptype, &pindex,
                            &ntri, &ptris, &ptric);
    if (status != EGADS_SUCCESS) return status;
    fprintf(fp, " %d %d\n", len, ntri);
    if ((len == 0) || (ntri == 0)) continue;
    for (j = 0; j < len; j++)
      fprintf(fp, "%19.12le %19.12le %19.12le %19.12le %19.12le %d %d\n",
              pxyz[3*j], pxyz[3*j+1], pxyz[3*j+2], puv[2*j], puv[2*j+1],
              ptype[j], pindex[j]);
    for (j = 0; j < ntri; j++)
      fprintf(fp, "%d %d %d %d %d %d\n", ptris[3*j], ptris[3*j+1], ptris[3*j+2],
              ptric[3*j], ptric[3*j+1], ptric[3*j+2]);
  }

  // write out the tessellation attributes
  attrs = (egAttrs *) tess->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;
  fprintf(fp, "%d\n", nattr);
  if (nattr != 0) EG_writeAttr(attrs, fp);

  return EGADS_SUCCESS;
}


int
EG_saveModel(const egObject *model, const char *name)
{
  int            i, j, n, len, outLevel, ibody, nbody, stat;
  TopoDS_Shape   wshape;
  const egObject **objs;
  FILE           *fp;
  
  if  (model == NULL)               return EGADS_NULLOBJ;
  if  (model->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((model->oclass != MODEL) &&
      (model->oclass != BODY))      return EGADS_NOTMODEL;
  outLevel = EG_outLevel(model);

  if (name == NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: NULL Filename (EG_saveModel)!\n");
    return EGADS_NONAME;
  }
  
  /* does file exist? */

  fp = fopen(name, "r");
  if (fp != NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Exists (EG_saveModel)!\n", name);
    fclose(fp);
    return EGADS_EXISTS;
  }
  
  /* find extension */
  
  len = strlen(name);
  for (i = len-1; i > 0; i--)
    if (name[i] == '.') break;
  if (i == 0) {
    if (outLevel > 0)
      printf(" EGADS Warning: No Extension in %s (EG_saveModel)!\n", name);
    return EGADS_NODATA;
  }
  
  if (model->oclass == MODEL) {
    egadsModel *mshape = (egadsModel *) model->blind;
    wshape             = mshape->shape;
    nbody              = mshape->nbody;
    objs               = (const egObject **) mshape->bodies;
  } else {
    egadsBody *mshape  = (egadsBody *) model->blind;
    wshape             = mshape->shape;
    nbody              = 1;
    objs               = &model;
  }
  
  if ((strcasecmp(&name[i],".step") == 0) || 
      (strcasecmp(&name[i],".stp") == 0)) {

    /* STEP files */
    
    STEPControl_Writer aWriter;
    TopExp_Explorer Exp;
    const STEPControl_StepModelType aVal = STEPControl_AsIs;
    for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);
    for (Exp.Init(wshape, TopAbs_SOLID);
         Exp.More(); Exp.Next()) aWriter.Transfer(Exp.Current(), aVal);	
  
    if (!aWriter.Write(name)) {
      printf(" EGADS Warning: STEP Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }
    
  } else if ((strcasecmp(&name[i],".iges") == 0) || 
             (strcasecmp(&name[i],".igs") == 0)) {
             
    /* IGES files */

    try {
      IGESControl_Controller::Init();
      IGESControl_Writer iWrite;
      TopExp_Explorer Exp;
      for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());
      for (Exp.Init(wshape, TopAbs_SOLID);
           Exp.More(); Exp.Next()) iWrite.AddShape(Exp.Current());

      iWrite.ComputeModel();
      if (!iWrite.Write(name)) {
        printf(" EGADS Warning: IGES Write Error (EG_saveModel)!\n");
        return EGADS_WRITERR;
      }
    }
    catch (...)
    {
      printf(" EGADS Warning: Internal IGES Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }

  } else if ((strcasecmp(&name[i],".brep") == 0) ||
             (strcasecmp(&name[i],".egads") == 0)) {
  
    /* Native OCC file or our filetype */

    if (!BRepTools::Write(wshape, name)) {
      printf(" EGADS Warning: OCC Write Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }
    if (strcasecmp(&name[i],".brep") == 0) return EGADS_SUCCESS;
    
    /* append the attributes -- output in the read order */
    
    FILE *fp = fopen(name, "a");
    if (fp == NULL) {
      printf(" EGADS Warning: EGADS Open Error (EG_saveModel)!\n");
      return EGADS_WRITERR;
    }
    fprintf(fp, "\n##EGADS HEADER FILE-REV 1 ##\n");
    /* write model attributes */
    if (model->oclass == MODEL) {
      EG_writeAttrs(model, fp);
    } else {
      fprintf(fp, "0\n");
    }
    TopExp_Explorer Exp;
    for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    for (Exp.Init(wshape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
      TopoDS_Shape shape = Exp.Current();
      for (i = 0; i < nbody; i++) {
        const egObject *obj   = objs[i];
        egadsBody      *pbody = (egadsBody *) obj->blind;
        if (shape.IsSame(pbody->shape)) {
          EG_writeAttrs(obj, fp);
          break;
        }
      }
    }
    /* add non-Body types */
    if ((model->oclass == MODEL) && (model->mtype > nbody)) {
      fprintf(fp, "##EGADS HEADER FILE-EXT 1 ##\n");
      fprintf(fp, "%hd\n", model->mtype);
      for (j = nbody; j < model->mtype; j++) {
        const egObject *obj = objs[j];
        TopoDS_Shape bshape;
        if (obj->oclass == TESSELLATION) {
          egTessel  *btess = (egTessel *) obj->blind;
          egObject  *src   = btess->src;
          if (src->oclass == EBODY) {
            for (i = nbody; i < model->mtype; i++) {
              const egObject *bod = objs[i];
              if (bod == src) break;
            }
            if (i == model->mtype) {
              printf(" EGADS Internal: Ancillary tess %d -- cannot find EBody!\n",
                     j+1);
              break;
            }
            fprintf(fp, "%d %d\n", obj->oclass, i+1);
            stat = EG_writeTess(obj, fp);
            if (stat != EGADS_SUCCESS) {
              printf(" EGADS Error: Ancillary tessellation %d -- status = %d\n",
                     j+1, stat);
              break;
            }
            continue;
          }
          egadsBody *pbody = (egadsBody *) src->blind;
          bshape           = pbody->shape;
        } else {
          egEBody   *ebody = (egEBody *) obj->blind;
          egObject  *ref   = ebody->ref;
          egadsBody *pbody = (egadsBody *) ref->blind;
          bshape           = pbody->shape;
        }
        n = ibody = 0;
        for (Exp.Init(wshape, TopAbs_WIRE,  TopAbs_FACE); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_FACE,  TopAbs_SHELL); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_SHELL, TopAbs_SOLID); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        for (Exp.Init(wshape, TopAbs_SOLID); Exp.More(); Exp.Next()) {
          TopoDS_Shape shape = Exp.Current();
          for (i = 0; i < nbody; i++) {
            const egObject *obj   = objs[i];
            egadsBody      *pbody = (egadsBody *) obj->blind;
            if (shape.IsSame(pbody->shape)) {
              n++;
              if (shape.IsSame(bshape)) ibody = n;
              break;
            }
          }
        }
        if ((n != nbody) || (ibody == 0)) {
          printf(" EGADS Internal: Ancillary objects -- n = %d [%d], ib = %d\n",
                 n, nbody, ibody);
          break;
        }
        fprintf(fp, "%d %d\n", obj->oclass, ibody);
        if (obj->oclass == TESSELLATION) {
          stat = EG_writeTess(obj, fp);
        } else {
          stat = EG_writeEBody(obj, fp);
        }
        if (stat != EGADS_SUCCESS) {
          printf(" EGADS Error: Ancillary objects %d -- status = %d\n",
                 j+1, stat);
          break;
        }
      }
    }
    fclose(fp);

  } else {
    if (outLevel > 0)
      printf(" EGADS Warning: Extension in %s Not Supported (EG_saveModel)!\n",
             name);
    return EGADS_NODATA;
  }

  return EGADS_SUCCESS;
}


int
EG_saveTess(egObject *tess, const char *name)
{
  int          status, len,   outLevel, stat;
  int          npts,   ntri,  iedge,    iface;
  int          nnode,  nedge, nface,    nattr = 0;
  const double *pxyz  = NULL, *puv    = NULL, *pt    = NULL;
  const int    *ptype = NULL, *pindex = NULL, *ptris = NULL, *ptric = NULL;
  egObject     *body;
  egAttrs      *attrs;

  if (tess == NULL)                 return EGADS_NULLOBJ;
  if (tess->magicnumber != MAGIC)   return EGADS_NOTOBJ;
  if (tess->oclass != TESSELLATION) return EGADS_NOTTESS;
  outLevel = EG_outLevel(tess);

  /* does file exist? */

  FILE *fp = fopen(name, "r");
  if (fp != NULL) {
    if (outLevel > 0)
      printf(" EGADS Warning: File %s Exists (EG_saveTess)!\n", name);
    status = EGADS_EXISTS;
    goto cleanup;
  }

  // get the body from tessellation
  status = EG_statusTessBody(tess, &body, &stat, &npts);
  if (status != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Warning: Tessellation is Open (EG_saveTess)!\n");
    status = EGADS_TESSTATE;
    goto cleanup;
  }

  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;
  status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) goto cleanup;

  fp = fopen(name, "w");
  if (fp == NULL) {
    printf(" EGADS Warning: File %s Open Error (EG_saveTess)!\n", name);
    status = EGADS_WRITERR;
    goto cleanup;
  }

  // write the number of nodes, edges, and faces
  fwrite(&nnode, sizeof(int), 1, fp);
  fwrite(&nedge, sizeof(int), 1, fp);
  fwrite(&nface, sizeof(int), 1, fp);

  // write out the edge tessellation
  for (iedge = 0; iedge < nedge; iedge++) {

    status = EG_getTessEdge(tess, iedge+1, &len, &pxyz, &pt);
    if (status != EGADS_SUCCESS) goto cleanup;

    fwrite(&len, sizeof(int),        1, fp);
    if (len == 0) continue;
    fwrite(pxyz, sizeof(double), 3*len, fp);
    fwrite(pt,   sizeof(double),   len, fp);
    if (ferror(fp)) {
      printf ("EGADS Warning: File %s Write Error (EG_saveTess)!\n", name);
      status = EGADS_WRITERR;
      goto cleanup;
    }
  }

  // write out face tessellations
  for (iface = 0; iface < nface; iface++) {
    status = EG_getTessFace(tess, iface+1, &len, &pxyz, &puv, &ptype, &pindex,
                            &ntri, &ptris, &ptric);
    if (status != EGADS_SUCCESS) goto cleanup;

    fwrite(&len,   sizeof(int),         1, fp);
    fwrite(&ntri,  sizeof(int),         1, fp);
    if ((len == 0) || (ntri == 0)) continue;
    fwrite(pxyz,   sizeof(double),  3*len, fp);
    fwrite(puv,    sizeof(double),  2*len, fp);
    fwrite(ptype,  sizeof(int),       len, fp);
    fwrite(pindex, sizeof(int),       len, fp);
    fwrite(ptris,  sizeof(int),    3*ntri, fp);
    fwrite(ptric,  sizeof(int),    3*ntri, fp);
    if (ferror(fp)) {
      printf ("EGADS Warning: File %s Write Error (EG_saveTess)!\n", name);
      status = EGADS_WRITERR;
      goto cleanup;
    }
  }

  // write out the tessellation attributes
  attrs = (egAttrs *) tess->attrs;
  if (attrs != NULL) nattr = attrs->nattrs;
  fprintf(fp, "%d\n", nattr);
  if (nattr != 0) EG_writeAttr(attrs, fp);

  status = EGADS_SUCCESS;

cleanup:
  fclose(fp);

  return status;
}


int
EG_loadTess(egObject *body, const char *name, egObject **tess)
{
  int      i, ir, status, nnode, nedge, nface, n[3], len, ntri, nattr;
  int      *ptype, *pindex, *tris, *tric;
  double   *xyz, *param;
  egObject *obj;
  FILE     *fp;
  
  *tess  = NULL;
  status = EG_getBodyTopos(body, NULL, NODE, &nnode, NULL);
  if (status != EGADS_SUCCESS) return status;
  status = EG_getBodyTopos(body, NULL, EDGE, &nedge, NULL);
  if (status != EGADS_SUCCESS) return status;
  status = EG_getBodyTopos(body, NULL, FACE, &nface, NULL);
  if (status != EGADS_SUCCESS) return status;
  
  fp = fopen(name, "r");
  if (fp == NULL) {
    printf(" EGADS Error: File %s Does not Exist (EG_loadTess)!\n", name);
    return EGADS_EXISTS;
  }
  
  ir = fread(n, sizeof(int), 3, fp);
  if (ir != 3) {
    printf(" EGADS Error: Header with only %d words (EG_loadTess)!\n", ir);
    fclose(fp);
    return EGADS_INDEXERR;
  }
  if ((nnode != n[0]) || (nedge != n[1]) || (nface != n[2])) {
    printf(" EGADS Error: Count mismatch %d %d  %d %d  %d %d (EG_loadTess)!\n",
           nnode, n[0], nedge, n[1], nface, n[2]);
    fclose(fp);
    return EGADS_INDEXERR;
  }
  
  /* initialize the Tessellation Object */
  status = EG_initTessBody(body, tess);
  if (status != EGADS_SUCCESS) {
    fclose(fp);
    return status;
  }
  
  /* do the Edges */
  for (i = 0; i < nedge; i++) {
    len = 0;
    fread(&len, sizeof(int), 1, fp);
    if (len == 0) continue;
    xyz   = (double *) malloc(3*len*sizeof(double));
    param = (double *) malloc(  len*sizeof(double));
    if ((xyz == NULL) || (param == NULL)) {
      printf(" EGADS Error: malloc on Edge %d -- len = %d (EG_loadTess)!\n",
               i+1, len);
      if (xyz   != NULL) free(xyz);
      if (param != NULL) free(param);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return EGADS_MALLOC;
    }
    fread(xyz,   sizeof(double), 3*len, fp);
    fread(param, sizeof(double),   len, fp);
    status = EG_setTessEdge(*tess, i+1, len, xyz, param);
    free(xyz);
    free(param);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessEdge %d = %d (EG_loadTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return status;
    }
  }
  
  /* do the Faces */
  for (i = 0; i < nface; i++) {
    len = ntri = 0;
    fread(&len,  sizeof(int), 1, fp);
    fread(&ntri, sizeof(int), 1, fp);
    if ((len == 0) || (ntri == 0)) continue;
    xyz    = (double *) malloc(3*len*sizeof(double));
    param  = (double *) malloc(2*len*sizeof(double));
    ptype  = (int *)    malloc(  len* sizeof(int));
    pindex = (int *)    malloc(  len* sizeof(int));
    tris   = (int *)    malloc(3*ntri*sizeof(int));
    tric   = (int *)    malloc(3*ntri*sizeof(int));
    if ((xyz    == NULL) || (param == NULL) || (ptype == NULL) ||
        (pindex == NULL) || (tris  == NULL) || (tric  == NULL)) {
      printf(" EGADS Error: malloc on Face %d -- lens = %d %d (EG_loadTess)!\n",
             i+1, len, ntri);
      if (xyz    != NULL) free(xyz);
      if (param  != NULL) free(param);
      if (ptype  != NULL) free(ptype);
      if (pindex != NULL) free(pindex);
      if (tris   != NULL) free(tris);
      if (tric   != NULL) free(tric);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return EGADS_MALLOC;
    }
    fread(xyz,    sizeof(double), 3*len,  fp);
    fread(param,  sizeof(double), 2*len,  fp);
    fread(ptype,  sizeof(int),      len,  fp);
    fread(pindex, sizeof(int),      len,  fp);
    fread(tris,   sizeof(int),    3*ntri, fp);
    fread(tric,   sizeof(int),    3*ntri, fp);
    status = EG_setTessFace(*tess, i+1, len, xyz, param, ntri, tris);
    free(xyz);
    free(param);
    free(ptype);
    free(pindex);
    free(tris);
    free(tric);
    if (status != EGADS_SUCCESS) {
      printf(" EGADS Error: EG_setTessFace %d = %d (EG_loadTess)!\n",
             i+1, status);
      EG_deleteObject(*tess);
      *tess = NULL;
      fclose(fp);
      return status;
    }
  }
  
  /* close up the open tessellation */
  status = EG_statusTessBody(*tess, &obj, &i, &len);
  if (status != EGADS_SUCCESS) {
    printf(" EGADS Error: EG_statusTessBody = %d (EG_loadTess)!\n", status);
    EG_deleteObject(*tess);
    *tess = NULL;
    fclose(fp);
    return status;
  }
  if (i != 1) {
    printf(" EGADS Error: New Tessellation Object is Open (EG_loadTess)!\n");
    EG_deleteObject(*tess);
    *tess = NULL;
    fclose(fp);
    return EGADS_TESSTATE;
  }
  
  /* attach the attributes */
  fscanf(fp, "%d\n", &nattr);
  if (nattr != 0) EG_readAttrs(*tess, nattr, fp);
  fclose(fp);
  
  return EGADS_SUCCESS;
}
