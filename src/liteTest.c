/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Lite Read & Tester
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"

#ifdef WIN32
#define LONG long long
#else
#define LONG long
#endif

static void
attrOut(int level, ego object)
{
  int          i, j, k, stat, nattr, atype, alen;
  const int    *ints;
  const double *reals;
  const char   *name, *str;
  LONG         pointer;

  stat = EG_attributeNum(object, &nattr);
  if (stat  != EGADS_SUCCESS) return;
  if (nattr == 0)             return;
  
  for (i = 1; i <= nattr; i++) {
    stat = EG_attributeGet(object, i, &name, &atype, &alen, &ints, &reals, &str);
    if (stat != EGADS_SUCCESS) continue;
    for (k = 0; k < 2*level+2; k++) printf(" ");
    printf("Attr: %s = ", name);
    if (atype == ATTRINT) {
      for (j = 0; j < alen; j++) printf("%d ", ints[j]);
    } else if ((atype == ATTRREAL) || (atype == ATTRCSYS)) {
      for (j = 0; j < alen; j++) printf("%lf ", reals[j]);
    } else if (atype == ATTRPTR) {
      pointer = (LONG) str;
#ifdef WIN32
      printf("%llx", pointer);
#else
      printf("%lx", pointer);
#endif
    } else {
      printf("%s", str);
    }
    printf("\n");
  }
}


static void
parseOut(int level, ego object, /*@null@*/ ego body, int sense)
{
  int    i, stat, oclass, mtype, nobjs, periodic, index, *senses, *ivec;
  ego    geom, *objs;
  double limits[4], bbox[6], *rvec;
  LONG   pointer;
  static char *classType[27] = {"CONTEXT", "TRANSFORM", "TESSELLATION",
                                "NIL", "EMPTY", "REFERENCE", "", "",
                                "", "", "PCURVE", "CURVE", "SURFACE", "", 
                                "", "", "", "", "", "", "NODE",
                                "EGDE", "LOOP", "FACE", "SHELL",
                                "BODY", "MODEL"};
  static char *curvType[9] = {"Line", "Circle", "Ellipse", "Parabola",
                              "Hyperbola", "Trimmed", "Bezier", "BSpline", 
                              "Offset"};
  static char *surfType[11] = {"Plane", "Spherical", "Cylinder", "Revolution",
                               "Toroidal", "Trimmed" , "Bezier", "BSpline", 
                               "Offset", "Conical", "Extrusion"};
  
  pointer = (LONG) object;
  oclass  = object->oclass;
  mtype   = object->mtype;
  
  /* geometry */
  if ((oclass >= PCURVE) && (oclass <= SURFACE)) {
    stat = EG_getGeometry(object, &oclass, &mtype, &geom, &ivec, &rvec);
    if (stat != EGADS_SUCCESS) {
      printf(" parseOut: %d EG_getGeometry return = %d\n", level, stat);
      return;
    }
    stat = EG_getRange(object, limits, &periodic);

    /* output most data except the axes info */
    if (oclass != SURFACE) {

      for (i = 0; i < 2*level; i++) printf(" ");
#ifdef WIN32
      printf("%s %llx  range = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1], periodic);
#else
      printf("%s %lx  range = %le %le  per = %d\n", 
             classType[oclass], pointer, limits[0], limits[1], periodic);
#endif

      for (i = 0; i < 2*level+2; i++) printf(" ");
      if (oclass == PCURVE) {
        switch (mtype) {
          case CIRCLE:
            printf("%s  radius = %lf\n", curvType[mtype-1], rvec[6]);
            break;
          case ELLIPSE:
            printf("%s  major = %lf, minor = %lf\n", 
                   curvType[mtype-1], rvec[6], rvec[7]);
            break;
          case PARABOLA:
            printf("%s  focus = %lf\n", curvType[mtype-1], rvec[6]);
            break;
          case HYPERBOLA:
            printf("%s  major = %lf, minor = %lf\n", 
                   curvType[mtype-1], rvec[6], rvec[7]);
            break;
          case TRIMMED:
            printf("%s  first = %lf, last = %lf\n", 
                   curvType[mtype-1], rvec[0], rvec[1]);
            break;
          case BEZIER:
            printf("%s  flags = %x, degree = %d, #CPs = %d\n", 
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2]);
            break;       
          case BSPLINE:
            printf("%s  flags = %x, degree = %d, #CPs = %d, #knots = %d\n", 
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
/*          for (i = 0; i < 2*level+2; i++) printf(" ");
            printf("knots =");
            for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
            printf("\n");  */
            break;
          case OFFSET:
            printf("%s  offset = %lf\n", curvType[mtype-1], rvec[0]);
            break;
          case 0:
            printf("unknown curve type!\n");
            break;
          default:
            printf("%s   %lf %lf   %lf %lf\n", curvType[mtype-1],
                   rvec[0], rvec[1], rvec[2], rvec[3]);
        }
      } else {
        switch (mtype) {
          case CIRCLE:
            printf("%s  radius = %lf\n", curvType[mtype-1], rvec[9]);
            break;
          case ELLIPSE:
            printf("%s  major = %lf, minor = %lf\n", 
                   curvType[mtype-1], rvec[9], rvec[10]);
            break;
          case PARABOLA:
            printf("%s  focus = %lf\n", curvType[mtype-1], rvec[9]);
            break;
          case HYPERBOLA:
            printf("%s  major = %lf, minor = %lf\n", 
                   curvType[mtype-1], rvec[9], rvec[10]);
            break;
          case TRIMMED:
            printf("%s  first = %lf, last = %lf\n", 
                   curvType[mtype-1], rvec[0], rvec[1]);
            break;
          case BEZIER:
            printf("%s  flags = %x, degree = %d, #CPs = %d\n", 
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2]);
            break;       
          case BSPLINE:
            printf("%s  flags = %x, degree = %d, #CPs = %d, #knots = %d\n", 
                   curvType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
/*          for (i = 0; i < 2*level+2; i++) printf(" ");
            printf("knots =");
            for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
            printf("\n");  */
            break;
          case OFFSET:
            printf("%s  offset = %lf\n", curvType[mtype-1], rvec[3]);
            break;
          case 0:
            printf("unknown curve type!\n");
            break;
          default:
            printf("%s\n", curvType[mtype-1]);
        }
      }

    } else {
    
      for (i = 0; i < 2*level; i++) printf(" ");
#ifdef WIN32
      printf("%s %llx  Urange = %le %le  Vrange = %le %le  per = %d\n",
             classType[oclass], pointer, limits[0], limits[1],
                                         limits[2], limits[3], periodic);
#else
      printf("%s %lx  Urange = %le %le  Vrange = %le %le  per = %d\n", 
             classType[oclass], pointer, limits[0], limits[1], 
                                         limits[2], limits[3], periodic);
#endif

      for (i = 0; i < 2*level+2; i++) printf(" ");
      switch (mtype) {
        case SPHERICAL:
          printf("%s  radius = %lf\n", surfType[mtype-1], rvec[9]);
          break;
        case CONICAL:
          printf("%s  angle = %lf, radius = %lf\n", 
                 surfType[mtype-1], rvec[12], rvec[13]);
          printf("    rvec = %lf %lf %lf   %lf %lf %lf  \n", 
                 rvec[0], rvec[1], rvec[2], rvec[3], rvec[4],  rvec[5]);
          printf("           %lf %lf %lf   %lf %lf %lf  \n", 
                 rvec[6], rvec[7], rvec[8], rvec[9], rvec[10], rvec[11]);
          break;
        case CYLINDRICAL:
          printf("%s  radius = %lf\n", surfType[mtype-1], rvec[12]);
          break;
        case TOROIDAL:
          printf("%s  major = %lf, minor = %lf\n", 
                 surfType[mtype-1], rvec[12], rvec[13]);
          break;
        case BEZIER:
          printf("%s  flags = %x, U deg = %d #CPs = %d, V deg = %d #CPs = %d\n", 
                 surfType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3], ivec[4]);
          break;
        case BSPLINE:
          printf("%s  flags = %x, U deg = %d #CPs = %d #knots = %d ",
                 surfType[mtype-1], ivec[0], ivec[1], ivec[2], ivec[3]);
          printf(" V deg = %d #CPs = %d #knots = %d\n",
                 ivec[4], ivec[5], ivec[6]);
/*        for (i = 0; i < 2*level+2; i++) printf(" ");
          printf("Uknots =");
          for (i = 0; i < ivec[3]; i++) printf(" %lf", rvec[i]);
          for (i = 0; i < 2*level+2; i++) printf(" ");
          printf("\nVknots =");
          for (i = 0; i < ivec[6]; i++) printf(" %lf", rvec[i+ivec[3]]);
          printf("\n"); */
          break;
        case TRIMMED:
          printf("%s  U trim = %lf %lf, V trim = %lf %lf\n", 
                 surfType[mtype-1], rvec[0], rvec[1], rvec[2], rvec[3]);
          break;
        case OFFSET:
          printf("%s  offset = %lf\n", surfType[mtype-1], rvec[0]);
          break;
        case 0:
          printf("unknown surface type!\n");
          break;
        default:
          printf("%s\n", surfType[mtype-1]);
      }
    }

    if (ivec != NULL) EG_free(ivec);
    if (rvec != NULL) EG_free(rvec);
    attrOut(level, object);
    if (geom != NULL) parseOut(level+1, geom, body, 0);
    return;
  }

  /* output class and pointer data */

  for (i = 0; i < 2*level; i++) printf(" ");
  index = 0;
  if ((oclass >= NODE) && (oclass < BODY))
    if (body != NULL) index = EG_indexBodyTopo(body, object);
#ifdef WIN32
  if (sense == 0) {
    printf("%s %llx %d\n", classType[oclass], pointer, index);
  } else {
    printf("%s %llx %d  sense = %d\n", classType[oclass], pointer, index, sense);
  }
#else
  if (sense == 0) {
    printf("%s %lx %d\n", classType[oclass], pointer, index);
  } else {
    printf("%s %lx %d  sense = %d\n", classType[oclass], pointer, index, sense);
  }
#endif

  /* topology*/
  if ((oclass >= NODE) && (oclass <= MODEL)) {
    stat = EG_getTopology(object, &geom, &oclass, &mtype, limits, &nobjs,
                          &objs, &senses);
    if (stat != EGADS_SUCCESS) {
      printf(" parseOut: %d EG_getTopology return = %d\n", level, stat);
      return;
    }
    if (oclass == NODE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      printf("XYZ = %lf %lf %lf\n", limits[0], limits[1], limits[2]);
    } else if (oclass == EDGE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      if (mtype == DEGENERATE) {
        printf("tRange = %lf %lf -- Degenerate!\n", limits[0], limits[1]);
      } else {
        printf("tRange = %lf %lf\n", limits[0], limits[1]);
      }
    } else if (oclass == FACE) {
      for (i = 0; i < 2*level+2; i++) printf(" ");
      printf("uRange = %lf %lf, vRange = %lf %lf\n", limits[0], limits[1], 
                                                     limits[2], limits[3]);
    }
    if (oclass != NODE) {
      stat = EG_getBoundingBox(object, bbox);
      if (stat == EGADS_SUCCESS) {
        for (i = 0; i < 2*level+2; i++) printf(" ");
        printf("BBox = %lf %lf %lf  %lf %lf %lf\n", bbox[0], bbox[1], bbox[2],
               bbox[3], bbox[4], bbox[5]);
      }
    }
    attrOut(level, object);
    
    if ((geom != NULL) && (mtype != DEGENERATE))
      parseOut(level+1, geom, body, 0);
    if (senses == NULL) {
      if (oclass == MODEL) {
        for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], objs[i], 0);
      } else {
        for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], body, 0);
      }
    } else {
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i], body, senses[i]);
    }
    if ((geom != NULL) && (oclass == LOOP))
      for (i = 0; i < nobjs; i++) parseOut(level+1, objs[i+nobjs], body, 0);
    
  }

}


int main(int argc, char *argv[])
{
  int i, j, k, n, nn, stat, oclass, mtype, nbodies, *senses;
  ego context, model, geom, *bodies, *objs, *nobjs;
  
  if (argc != 2) {
    printf(" Usage: liteTest liteFile\n\n");
    exit(EXIT_FAILURE);
  }
  /* initialize */
  printf(" EG_open          = %d\n", EG_open(&context));
  printf(" EG_loadModel     = %d  %s\n", EG_loadModel(context, 0, argv[1],
                                                      &model), argv[1]);
  
  /* test bodyTopo functions */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies,
                        &bodies, &senses);
  if (stat == EGADS_SUCCESS)
    for (i = 0; i < nbodies; i++) {
      stat = EG_getBodyTopos(bodies[i], NULL, NODE, &n, &objs);
      if (stat != EGADS_SUCCESS) {
        printf(" ERROR: %d Node EG_getBodyTopos = %d!\n", i, stat);
        continue;
      }
      for (j = 0; j < n; j++) {
        k = EG_indexBodyTopo(bodies[i], objs[j]);
        if (k != j+1) printf("  Node Index = %d but should be %d!\n", k, j+1);
      }
      EG_free(objs);
      stat = EG_getBodyTopos(bodies[i], NULL, EDGE, &n, &objs);
      if (stat != EGADS_SUCCESS) {
        printf(" ERROR: %d Edge EG_getBodyTopos = %d!\n", i, stat);
        continue;
      }
      for (j = 0; j < n; j++) {
        k = EG_indexBodyTopo(bodies[i], objs[j]);
        if (k != j+1) printf("  Edge Index = %d but should be %d!\n", k, j+1);
      }
      EG_free(objs);
      stat = EG_getBodyTopos(bodies[i], NULL, LOOP, &n, &objs);
      if (stat != EGADS_SUCCESS) {
        printf(" ERROR: %d Loop EG_getBodyTopos = %d!\n", i, stat);
        continue;
      }
      for (j = 0; j < n; j++) {
        k = EG_indexBodyTopo(bodies[i], objs[j]);
        if (k != j+1) printf("  Loop Index = %d but should be %d!\n", k, j+1);
        stat = EG_getBodyTopos(bodies[i], objs[j], NODE, &nn, &nobjs);
        if (stat != EGADS_SUCCESS) {
          printf("  Loop Index = %d GetNodes = %d\n", j+1, stat);
          continue;
        }
        printf("     loop %d has %d Nodes!\n", j, nn);
        EG_free(nobjs);
        stat = EG_getBodyTopos(bodies[i], objs[j], FACE, &nn, &nobjs);
        if (stat != EGADS_SUCCESS) {
          printf("  Loop Index = %d GetFaces = %d\n", j+1, stat);
          continue;
        }
        printf("     loop %d has %d Faces!\n", j, nn);
        EG_free(nobjs);
      }
      EG_free(objs);
      stat = EG_getBodyTopos(bodies[i], NULL, FACE, &n, &objs);
      if (stat != EGADS_SUCCESS) {
        printf(" ERROR: Face %d EG_getBodyTopos = %d!\n", i, stat);
        continue;
      }
      for (j = 0; j < n; j++) {
        k = EG_indexBodyTopo(bodies[i], objs[j]);
        if (k != j+1) printf("  Face Index = %d but should be %d!\n", k, j+1);
      }
      if (objs != NULL) EG_free(objs);
      stat = EG_getBodyTopos(bodies[i], NULL, SHELL, &n, &objs);
      if (stat != EGADS_SUCCESS) {
        printf(" ERROR: Shell %d EG_getBodyTopos = %d!\n", i, stat);
        continue;
      }
      if (objs != NULL) {
        for (j = 0; j < n; j++) {
          k = EG_indexBodyTopo(bodies[i], objs[j]);
          if (k != j+1) printf("  Shell Index = %d but should be %d!\n", k, j+1);
        }
        EG_free(objs);
      }
    }
  printf(" \n");
  
  /* output the entire model structure */
  parseOut(0, model, NULL, 0);
  printf(" \n");

  printf(" EG_close         = %d\n", EG_close(context));
  return 0;
}
