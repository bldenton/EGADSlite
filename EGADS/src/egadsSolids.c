/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Solid Primitives Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egads.h"
#include "egads_dot.h"
#include "egadsStack.h"

#if defined(_MSC_VER) && (_MSC_VER < 1900)
#define __func__  __FUNCTION__
#endif

//#define PCURVE_SENSITIVITY

#define PI               3.1415926535897931159979635

#define CROSS(a,b,c)  a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                      a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                      a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define CROSS_DOT(a_dot,b,b_dot,c,c_dot) \
  a_dot[0] = (b_dot[1]*c[2]) + (b[1]*c_dot[2]) - (b_dot[2]*c[1]) - (b[2]*c_dot[1]);\
  a_dot[1] = (b_dot[2]*c[0]) + (b[2]*c_dot[0]) - (b_dot[0]*c[2]) - (b[0]*c_dot[2]);\
  a_dot[2] = (b_dot[0]*c[1]) + (b[0]*c_dot[1]) - (b_dot[1]*c[0]) - (b[1]*c_dot[0])

#define DOT(a,b)      (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define MAX(a,b)      ((a) > (b) ? (a) : (b))
#define MIN(a,b)      ((a) < (b) ? (a) : (b))

extern int EG_outLevel( const ego object );


/* checks if an ego is of the expected class and type */
static int
EG_expectClassType(int oclass, int mtype, int eclass, int etype)
{
  static 
  const char *classType[27] = {"CONTEXT", "TRANSFORM", "TESSELLATION",
                               "NIL", "EMPTY", "REFERENCE", "", "",
                               "", "", "PCURVE", "CURVE", "SURFACE", "",
                               "", "", "", "", "", "", "NODE",
                               "EGDE", "LOOP", "FACE", "SHELL",
                               "BODY", "MODEL"};
  static 
  const char *bodyType[9] = {"", "", "", "", "",
                             "WIREBODY", "FACEBODY", "SHEETBODY", "SOLIDBODY"};
  static 
  const char *topoType[5] = {"ONENODE", "TWONODE", "OPEN", "CLOSED", "DEGENERATE"};
  static 
  const char *curvType[9] = {"LINE", "CIRCLE", "ELLIPSE", "PARABOLA",
                             "HYPERBOLA", "TRIMMED", "BEZIER", "BSPLINE",
                             "OFFSET"};
  static 
  const char *surfType[11] = {"PLANE", "SPHERICAL", "CYLINDER", "REVOLUTION",
                              "TOROIDAL", "TRIMMED" , "BEZIER", "BSPLINE",
                              "OFFSET", "CONICAL", "EXTRUSION"};
  const char **objType, **eType;

  if ((oclass != eclass) || (mtype != etype)) {
    if (oclass == BODY) {
      objType = bodyType;
    } else if ( (oclass == EDGE) || (oclass == LOOP) || (oclass == SHELL) ) {
      objType = topoType;
    } else if (oclass == SURFACE) {
      objType = surfType;
    } else if ((oclass == CURVE) || (oclass == PCURVE)) {
      objType = curvType;
    } else {
      printf(" EGADS Error: Unexpected oclass %d (%s)!\n",
             oclass, __func__);
      return EGADS_GEOMERR;
    }
    if (eclass == BODY) {
      eType = bodyType;
    } else if ( (eclass == EDGE) || (eclass == LOOP) || (eclass == SHELL) ) {
      eType = topoType;
    } else if (eclass == SURFACE) {
      eType = surfType;
    } else if ((eclass == CURVE) || (eclass == PCURVE)) {
      eType = curvType;
    } else {
      printf(" EGADS Error: Unexpected eclass %d (%s)!\n",
             eclass, __func__);
      return EGADS_GEOMERR;
    }
    printf(" EGADS Error: Object %s %s is not expected %s %s (%s)!\n",
           classType[oclass], objType[mtype-1], classType[eclass],
           eType[etype-1], __func__);
    return EGADS_GEOMERR;
  }

  return EGADS_SUCCESS;
}


static void
EG_Ax2(double *dirz, double *dirx, double *diry)
{
  double   A, Aabs, B, Babs, C, Cabs;

  /* taken from gp_Ax2 constructor
   *
   * This logic introduces discontinuities for finite differencing though...
   */
  A = dirz[0];
  B = dirz[1];
  C = dirz[2];
  Aabs = A;
  if (Aabs < 0) Aabs = - Aabs;
  Babs = B;
  if (Babs < 0) Babs = - Babs;
  Cabs = C;
  if (Cabs < 0) Cabs = - Cabs;

  if      ( Babs <= Aabs && Babs <= Cabs) {
    if (Aabs > Cabs) { dirx[0] = -C; dirx[1] = 0.; dirx[2] =  A; }
    else             { dirx[0] =  C; dirx[1] = 0.; dirx[2] = -A; }
  }
  else if ( Aabs <= Babs && Aabs <= Cabs) {
    if (Babs > Cabs) { dirx[0] = 0.; dirx[1] = -C; dirx[2] =  B; }
    else             { dirx[0] = 0.; dirx[1] =  C; dirx[2] = -B; }
  }
  else {
    if (Aabs > Babs) { dirx[0] = -B; dirx[1] =  A; dirx[2] = 0.; }
    else             { dirx[0] =  B; dirx[1] = -A; dirx[2] = 0.; }
  }

  CROSS(diry, dirz, dirx);
}


static void
EG_Ax2_dot(double *dirz, double *dirz_dot,
           double *dirx, double *dirx_dot,
           double *diry, double *diry_dot)
{
  double A, Aabs, B, Babs, C, Cabs;
  double A_dot, B_dot, C_dot;

  /* taken from gp_Ax2 constructor and linearized */
  A = dirz[0];
  B = dirz[1];
  C = dirz[2];

  A_dot = dirz_dot[0];
  B_dot = dirz_dot[1];
  C_dot = dirz_dot[2];

  Aabs = A;
  if (Aabs < 0) Aabs = -Aabs;
  Babs = B;
  if (Babs < 0) Babs = -Babs;
  Cabs = C;
  if (Cabs < 0) Cabs = -Cabs;

  if      (Babs <= Aabs && Babs <= Cabs) {
    if (Aabs > Cabs) {
      dirx[0]     = -C;     dirx[1]     = 0.; dirx[2]     =  A;
      dirx_dot[0] = -C_dot; dirx_dot[1] = 0.; dirx_dot[2] =  A_dot;
    } else {
      dirx[0]     =  C;     dirx[1]     = 0.; dirx[2]     = -A;
      dirx_dot[0] =  C_dot; dirx_dot[1] = 0.; dirx_dot[2] = -A_dot;
    }
  }
  else if (Aabs <= Babs && Aabs <= Cabs) {
    if (Babs > Cabs) {
      dirx[0]     = 0.; dirx[1]     = -C;     dirx[2]     =  B;
      dirx_dot[0] = 0.; dirx_dot[1] = -C_dot; dirx_dot[2] =  B_dot;
    } else {
      dirx[0]     = 0.; dirx[1]     =  C;     dirx[2]     = -B;
      dirx_dot[0] = 0.; dirx_dot[1] =  C_dot; dirx_dot[2] = -B_dot;
    }
  }
  else {
    if (Aabs > Babs) {
      dirx[0]     = -B;     dirx[1]     =  A;     dirx[2]     = 0.;
      dirx_dot[0] = -B_dot; dirx_dot[1] =  A_dot; dirx_dot[2] = 0.;
    } else {
      dirx[0]     =  B;     dirx[1]     = -A;     dirx[2]     = 0.;
      dirx_dot[0] =  B_dot; dirx_dot[1] = -A_dot; dirx_dot[2] = 0.;
    }
  }

  CROSS(diry, dirz, dirx);
  CROSS_DOT(diry_dot, dirz,dirz_dot, dirx,dirx_dot);
}


/* nodes that define each edge of a BOX */
static const int iedgeNodesBox[12][2] = {{1,0},
                                         {0,2},
                                         {3,2},
                                         {1,3},
                                         {5,4},
                                         {4,6},
                                         {7,6},
                                         {5,7},
                                         {1,5},
                                         {0,4},
                                         {3,7},
                                         {2,6}};


/* edges that make up each face in a BOX */
static const int ifaceEdgesBox[6][4] = {{ 0, 1, 2, 3},   /* x-min */
                                        { 4, 7, 6, 5},   /* x-max */
                                        { 8, 4, 9, 0},   /* y-min */
                                        {10, 2,11, 6},   /* y-max */
                                        { 3,10, 7, 8},   /* z-min */
                                        { 1, 9, 5,11}};  /* z-max */


int
EG_makeSolidBox(egObject *context, const double *data, egObject **body)
{
  /* all entities are consistent with OCC construction,
   * but the final body edge/node numbering is not consistent
   */

  /* extract box data */
  double x  = data[0];
  double y  = data[1];
  double z  = data[2];
  double dx = data[3];
  double dy = data[4];
  double dz = data[5];

  /* node coordinates */
  double rnode[8][3] = {{x,    y,    z+dz },
                        {x,    y,    z    },
                        {x,    y+dy, z+dz },
                        {x,    y+dy, z    },
                        {x+dx, y,    z+dz },
                        {x+dx, y,    z    },
                        {x+dx, y+dy, z+dz },
                        {x+dx, y+dy, z    }};

  /* surface data definitions */
  double rsurf[6][9] = { {x,    y,    z,     0, 0, 1,  0,-1, 0},    /* x-min */
                         {x+dx, y,    z,     0, 0, 1,  0,-1, 0},    /* x-max */
                         {x,    y,    z,     0, 0, 1,  1, 0, 0},    /* y-min */
                         {x,    y+dy, z,     0, 0, 1,  1, 0, 0},    /* y-max */
                         {x,    y,    z,     1, 0, 0,  0, 1, 0},    /* z-min */
                         {x,    y,    z+dz,  1, 0, 0,  0, 1, 0} };  /* z-max */

  /* senses of each modulus 2 face */
  int fsens[2] = {SREVERSE, SFORWARD};

  /* edge senses in each modulus 2 loop */
  int esens[2][4] = {{SFORWARD, SFORWARD, SREVERSE, SREVERSE},
                     {SREVERSE, SFORWARD, SFORWARD, SREVERSE}};

  /* loop sense for each face */
  int lsens[1] = {SFORWARD};

  int      stat = EGADS_SUCCESS;
  int      i, j, outLevel;
  double   xyz0[3], xyz1[3], rdata[6], tdata[2];
  ego      nodes[8], edges[12], faces[6];
  ego      ref, line, surf, enodes[2], ledges[4], loop, shell;
  objStack stack;

  outLevel = EG_outLevel(context);

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create nodes */
  for (i = 0; i < 8; i++) {
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rnode[i], 0, NULL, NULL, &nodes[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* create lines and edges */
  for (i = 0; i < 12; i++) {

    enodes[0] = nodes[iedgeNodesBox[i][0]];
    enodes[1] = nodes[iedgeNodesBox[i][1]];

    for (j = 0; j < 3; j++) {
      xyz0[j] = rnode[iedgeNodesBox[i][0]][j];
      xyz1[j] = rnode[iedgeNodesBox[i][1]][j];
    }

    rdata[0] = xyz0[0];
    rdata[1] = xyz0[1];
    rdata[2] = xyz0[2];
    rdata[3] = xyz1[0] - xyz0[0];
    rdata[4] = xyz1[1] - xyz0[1];
    rdata[5] = xyz1[2] - xyz0[2];

    stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL, rdata, &line);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, line);
    if (stat != EGADS_SUCCESS) goto cleanup;

    tdata[0] = 0;
    tdata[1] = sqrt(rdata[3]*rdata[3] + rdata[4]*rdata[4] + rdata[5]*rdata[5]);

    stat = EG_makeTopology(context, line, EDGE, TWONODE, tdata,
                           2, enodes, NULL, &edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* create surfaces, loops, and faces */
  for (i = 0; i < 6; i++) {

    stat = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL, rsurf[i], &surf);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, surf);
    if (stat != EGADS_SUCCESS) goto cleanup;

    for (j = 0; j < 4; j++)
      ledges[j] = edges[ifaceEdgesBox[i][j]];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL,
                           4, ledges, esens[i%2], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, surf, FACE, fsens[i%2], NULL,
                           1, &loop, lsens, &faces[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* make the shell */
  stat = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL,
                         6, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* make the final body */
  stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                         &shell, NULL, body);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (%s)!\n", i, __func__);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  return stat;
}


int
EG_makeSolidBox_dot(const double *data, const double *data_dot, egObject *body)
{
  /* extract box data */
  double x  = data[0];
  double y  = data[1];
  double z  = data[2];
  double dx = data[3];
  double dy = data[4];
  double dz = data[5];

  double x_dot  = data_dot[0];
  double y_dot  = data_dot[1];
  double z_dot  = data_dot[2];
  double dx_dot = data_dot[3];
  double dy_dot = data_dot[4];
  double dz_dot = data_dot[5];

  /* node coordinates */
  double rnode[8][3] = {{x,    y,    z+dz},
                        {x,    y,    z   },
                        {x,    y+dy, z+dz},
                        {x,    y+dy, z   },
                        {x+dx, y,    z+dz},
                        {x+dx, y,    z   },
                        {x+dx, y+dy, z+dz},
                        {x+dx, y+dy, z   }};

  double rnode_dot[8][3] = {{x_dot,        y_dot,        z_dot+dz_dot},
                            {x_dot,        y_dot,        z_dot       },
                            {x_dot,        y_dot+dy_dot, z_dot+dz_dot},
                            {x_dot,        y_dot+dy_dot, z_dot       },
                            {x_dot+dx_dot, y_dot,        z_dot+dz_dot},
                            {x_dot+dx_dot, y_dot,        z_dot       },
                            {x_dot+dx_dot, y_dot+dy_dot, z_dot+dz_dot},
                            {x_dot+dx_dot, y_dot+dy_dot, z_dot       }};

  /* surface data definitions */
  double rsurf[6][9] = { {x,    y,    z,      0, 0, 1,  0,-1, 0},    /* x-min */
                         {x+dx, y,    z,      0, 0, 1,  0,-1, 0},    /* x-max */
                         {x,    y,    z,      0, 0, 1,  1, 0, 0},    /* y-min */
                         {x,    y+dy, z,      0, 0, 1,  1, 0, 0},    /* y-max */
                         {x,    y,    z,      1, 0, 0,  0, 1, 0},    /* z-min */
                         {x,    y,    z+dz,   1, 0, 0,  0, 1, 0} };  /* z-max */

  double rsurf_dot[6][9] =
    { {x_dot,        y_dot,        z_dot       , 0, 0, 0,  0, 0, 0},   /* x-min */
      {x_dot+dx_dot, y_dot,        z_dot       , 0, 0, 0,  0, 0, 0},   /* x-max */
      {x_dot,        y_dot,        z_dot       , 0, 0, 0,  0, 0, 0},   /* y-min */
      {x_dot,        y_dot+dy_dot, z_dot       , 0, 0, 0,  0, 0, 0},   /* y-max */
      {x_dot,        y_dot,        z_dot       , 0, 0, 0,  0, 0, 0},   /* z-min */
      {x_dot,        y_dot,        z_dot+dz_dot, 0, 0, 0,  0, 0, 0} }; /* z-max */

  int      stat = EGADS_SUCCESS;
  int      i, j, oclass, mtype, outLevel, *senses;
  int      nshell, nface, nloop, nedge, nnode;
  double   xyz0[3], xyz0_dot[3], xyz1[3], xyz1_dot[3], rdata[6], rdata_dot[6];
  ego      nodes[8], edges[12], *faces, *shells;
  ego      ref, line, plane, *enodes, *ledges, *loops;

  outLevel = EG_outLevel(body);


  /* get the Shell from the Body */
  stat = EG_getTopology(body, &ref, &oclass, &mtype,
                        rdata, &nshell, &shells, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, BODY, SOLIDBODY);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (nshell != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: body must have one SHELL (%s)!\n", __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }


  /* get the Faces from the Shell */
  stat = EG_getTopology(shells[0], &ref, &oclass, &mtype,
                        rdata, &nface, &faces, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, SHELL, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (nface != 6) {
    if (outLevel > 0)
      printf(" EGADS Error: SHELL (with %d Faces) must have 6 Faces (%s)!\n",
             nface, __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* set surfaces sensitivity and extract edges from the loop in each face */
  for (i = 0; i < 6; i++) {

    /* get the Loop from the Face */
    stat = EG_getTopology(faces[i], &plane, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nloop != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
               i+1, nloop, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set the surface sensitivity */
    stat = EG_setGeometry_dot(plane, SURFACE, PLANE, NULL, rsurf[i], rsurf_dot[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;


    /* get the Edges from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if (nedge != 4) {
      if (outLevel > 0)
        printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
               i+1, nedge, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* save off the edges in the known order */
    for (j = 0; j < 4; j++) {
      edges[ifaceEdgesBox[i][j]] = ledges[j];
    }
  }

  /* extract all the nodes from the edges and set line sensitivity */
  for (i = 0; i < 12; i++) {

    /* get the Nodes from the Edge */
    stat = EG_getTopology(edges[i], &line, &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract the nodes */
    nodes[iedgeNodesBox[i][0]] = enodes[0];
    nodes[iedgeNodesBox[i][1]] = enodes[1];

    for (j = 0; j < 3; j++) {
      xyz0[j] = rnode[iedgeNodesBox[i][0]][j];
      xyz1[j] = rnode[iedgeNodesBox[i][1]][j];

      xyz0_dot[j] = rnode_dot[iedgeNodesBox[i][0]][j];
      xyz1_dot[j] = rnode_dot[iedgeNodesBox[i][1]][j];
    }

    rdata[0] = xyz0[0];
    rdata[1] = xyz0[1];
    rdata[2] = xyz0[2];
    rdata[3] = xyz1[0] - xyz0[0];
    rdata[4] = xyz1[1] - xyz0[1];
    rdata[5] = xyz1[2] - xyz0[2];

    rdata_dot[0] = xyz0_dot[0];
    rdata_dot[1] = xyz0_dot[1];
    rdata_dot[2] = xyz0_dot[2];
    rdata_dot[3] = xyz1_dot[0] - xyz0_dot[0];
    rdata_dot[4] = xyz1_dot[1] - xyz0_dot[1];
    rdata_dot[5] = xyz1_dot[2] - xyz0_dot[2];

    /* set the line sensitivity */
    stat = EG_setGeometry_dot(line, CURVE, LINE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* set node sensitivity */
  for (i = 0; i < 8; i++) {
    stat = EG_setGeometry_dot(nodes[i], NODE, 0, NULL, rnode[i], rnode_dot[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  return stat;
}


int
EG_makeSolidSphere(egObject *context, int stypx, const double *data,
                   egObject **body)
{
  /* center and radius of the sphere */
  double xcent[3] = {data[0], data[1], data[2]};
  double r        = data[3];

  int      stat = EGADS_SUCCESS;
  int      i, outLevel, nface, lsens[1] = {SFORWARD};
  int      esens[4] = {SREVERSE, SREVERSE, SFORWARD, SFORWARD};
  double   rdata[10], tdata[2];
  ego      sphere, trimmed, circles[2], nodes[2], edges[8], loop, faces[2];
  ego      edge, shell, ref;
  objStack stack;

  outLevel = EG_outLevel(context);

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Sphere */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = 1;        /* x-axis */
  rdata[4] = 0;
  rdata[5] = 0;
  rdata[6] = 0;        /* y-axis */
  rdata[7] = 1;
  rdata[8] = 0;
  rdata[9] = r;        /* radius */
  stat = EG_makeGeometry(context, SURFACE, SPHERICAL, NULL, NULL,
                         rdata, &sphere);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, sphere);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Circle curve for the u == 0 Edge */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = 1;        /* x-axis (consistent with sphere) */
  rdata[4] = 0;
  rdata[5] = 0;
  rdata[6] = 0;        /* y-axis (cross of sphere x-axis and y-axis) */
  rdata[7] = 0;
  rdata[8] = 1;
  rdata[9] = r;        /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Nodes for the Edge(s) */
  rdata[0] = xcent[0];
  rdata[1] = xcent[1];
  rdata[2] = xcent[2] - r;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  rdata[0] = xcent[0];
  rdata[1] = xcent[1];
  rdata[2] = xcent[2] + r;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    /* create the Circle curve for the u == PI Edge */
    rdata[0] = xcent[0]; /* center */
    rdata[1] = xcent[1];
    rdata[2] = xcent[2];
    rdata[3] = -1;       /* x-axis (reversed so sphere v is consistent with circle t) */
    rdata[4] =  0;
    rdata[5] =  0;
    rdata[6] =  0;       /* y-axis (cross of sphere x-axis and y-axis) */
    rdata[7] =  0;
    rdata[8] =  1;
    rdata[9] =  r;       /* radius */
    stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           rdata, &circles[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, circles[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Degenerate Edges on the Node */
    tdata[0] = PI;
    tdata[1] = 2.*PI;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[1], NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Edge on the second Circle */
    tdata[0] = -PI/2.;
    tdata[1] =  PI/2.;

    stat = EG_makeTopology(context, circles[1], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Degenerate Edges on the Node */
    tdata[0] = PI;
    tdata[1] = 2.*PI;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[0], NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the first Circle */
    tdata[0] = -PI/2.;
    tdata[1] =  PI/2.;

    stat = EG_makeTopology(context, circles[0], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* P-curve */
    rdata[0] = 0.;    rdata[1] =  PI/2;   /* v == PI/2 VMAX */
    rdata[2] = 1.;    rdata[3] =  0.  ;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = PI;    rdata[1] = 0.;      /* u == PI UMIN   */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 0.;    rdata[1] = -PI/2.;  /* v == -PI/2 VMIN */
    rdata[2] = 1.;    rdata[3] =  0.   ;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 2*PI;    rdata[1] = 0.;    /* u == 2*PI UMAX  */
    rdata[2] =   0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* trim the surface */
    rdata[0] =  PI;
    rdata[1] =  2.*PI;
    rdata[2] = -PI/2.;
    rdata[3] =  PI/2.;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, sphere, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, edges, esens, &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* re-make Degenerate Edges on the Nodes */
    tdata[0] = 0;
    tdata[1] = PI;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[1], NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[0], NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* swap the edges */
    edge = edges[3];
    edges[3] = edges[1];
    edges[1] = edge;

    /* create new P-curves */
    rdata[0] = 0;     rdata[1] = 0.;    /* u == 0 UMIN       */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = PI;    rdata[1] = 0.;    /* u == PI UMAX   */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* trim the surface */
    rdata[0] =  0.0;
    rdata[1] =  PI;
    rdata[2] = -PI/2.;
    rdata[3] =  PI/2.;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, sphere, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, edges, esens, &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 2;
    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

  } else {

    /* make the Degenerate Edges on the first Node */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[1], NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the u == 0 Circle */
    tdata[0] = -PI/2.;
    tdata[1] =  PI/2.;

    stat = EG_makeTopology(context, circles[0], EDGE, TWONODE,
                           tdata, 2, nodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Degenerate Edges on the second Node */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, &nodes[0], NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    edges[3] = edges[1]; /* repeat the circle edge */

    /* create P-curves */
    rdata[0] = 0.;    rdata[1] =  PI/2.; /* v == PI/2 VMAX */
    rdata[2] = 1.;    rdata[3] =  0.   ;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 0.;    rdata[1] = 0.;     /* u == 0 UMIN       */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 0.;    rdata[1] = -PI/2;  /* v == -PI/2 VMIN */
    rdata[2] = 1.;    rdata[3] =  0.  ;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 2*PI;  rdata[1] = 0.;    /* u == 2*PI UMAX   */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, sphere, LOOP, CLOSED,
                           NULL, 4, edges, esens, &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 1;
    stat = EG_makeTopology(context, sphere, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

  }

  /* make the shell */
  stat = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL,
                         nface, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* make the final body */
  stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                         &shell, NULL, body);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (%s)!\n", i, __func__);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  return stat;
}


int
EG_makeSolidSphere_dot(int stypx, const double *data, const double *data_dot,
                       egObject *body)
{
  /* center and radius of the sphere */
  double xcent[3] = {data[0], data[1], data[2]};
  double r        = data[3];

  double xcent_dot[3] = {data_dot[0], data_dot[1], data_dot[2]};
  double r_dot        = data_dot[3];

  int      stat = EGADS_SUCCESS;
#ifdef PCURVE_SENSITIVITY
  int      j, *ivec=NULL;
#endif
  int      i, oclass, mtype, outLevel, *senses;
  int      nshell, nface, nloop, nedge, nnode;
  double   rdata[10], rdata_dot[10], *rvec=NULL;
  ego      *faces, *shells;
  ego      ref, circle, trimmed, sphere, *enodes, *ledges, *loops;

  outLevel = EG_outLevel(body);

  /* get the Shell from the Body */
  stat = EG_getTopology(body, &ref, &oclass, &mtype,
                        rdata, &nshell, &shells, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, BODY, SOLIDBODY);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (nshell != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: body must have one SHELL (%s)!\n", __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* get the Faces from the Shell */
  stat = EG_getTopology(shells[0], &ref, &oclass, &mtype,
                        rdata, &nface, &faces, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, SHELL, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    if ( nface != 2 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 2 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

  } else {

    if ( nface != 1 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 1 Face (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

  }

  /* set surfaces sensitivity and extract edges from the loop in each face */
  for (i = 0; i < nface; i++) {

    /* get the Loop from the Face */
    stat = EG_getTopology(faces[i], &trimmed, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nloop != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
               i+1, nloop, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* get the PCurves from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    
#ifdef PCURVE_SENSITIVITY
    /* set PCURVE sensitivities */
    rdata_dot[0] = rdata_dot[1] = rdata_dot[2] = rdata_dot[3] = 0;
    for (j = 0; j < 4; j++) {
      /* generally calling EG_getGeometry to get rvec for EG_setGeometry_dot is not correct,
       * but here it is ok because the direction was specified as a unit vector */
      stat = EG_getGeometry(ledges[4+j], &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, rvec, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;
      EG_free(ivec); ivec=NULL;
      EG_free(rvec); rvec=NULL;
    }
#endif
    
    if (stypx > 0) { /* split periodic */

      /* set trim the surface sensitivity */
      if (i == 0) {
        rdata[0] =  PI;
        rdata[1] =  2.*PI;
        rdata[2] = -PI/2.;
        rdata[3] =  PI/2.;
      } else {
        rdata[0] =  0.0;
        rdata[1] =  PI;
        rdata[2] = -PI/2.;
        rdata[3] =  PI/2.;
      }

      rdata_dot[0] = 0;
      rdata_dot[1] = 0;
      rdata_dot[2] = 0;
      rdata_dot[3] = 0;

      stat = EG_setGeometry_dot(trimmed, SURFACE, TRIMMED, NULL, rdata, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Sphere from the Trimmed surface */
      stat = EG_getGeometry(trimmed, &oclass, &mtype, &sphere, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;

    } else {
      sphere = trimmed;
    }

    /* create the Sphere data */
    rdata[0] = xcent[0]; /* center */
    rdata[1] = xcent[1];
    rdata[2] = xcent[2];
    rdata[3] = 1;        /* x-axis */
    rdata[4] = 0;
    rdata[5] = 0;
    rdata[6] = 0;        /* y-axis */
    rdata[7] = 1;
    rdata[8] = 0;
    rdata[9] = r;        /* radius */

    rdata_dot[0] = xcent_dot[0]; /* center */
    rdata_dot[1] = xcent_dot[1];
    rdata_dot[2] = xcent_dot[2];
    rdata_dot[3] = 0;            /* x-axis */
    rdata_dot[4] = 0;
    rdata_dot[5] = 0;
    rdata_dot[6] = 0;            /* y-axis */
    rdata_dot[7] = 0;
    rdata_dot[8] = 0;
    rdata_dot[9] = r_dot;        /* radius */

    /* set the Sphere surface sensitivity */
    stat = EG_setGeometry_dot(sphere, SURFACE, SPHERE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* get the Edges from the Loop */
  stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                        rdata, &nedge, &ledges, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (nedge != 4) {
    if (outLevel > 0)
      printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
             i+1, nedge, __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* get the Nodes from the Edge with the u == 0 Circle */
  stat = EG_getTopology(ledges[1], &circle, &oclass, &mtype,
                        rdata, &nnode, &enodes, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the u == 0 Circle data */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = 1;        /* x-axis (consistent with sphere) */
  rdata[4] = 0;
  rdata[5] = 0;
  rdata[6] = 0;        /* y-axis (cross of sphere x-axis and y-axis) */
  rdata[7] = 0;
  rdata[8] = 1;
  rdata[9] = r;        /* radius */

  rdata_dot[0] = xcent_dot[0]; /* center */
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2];
  rdata_dot[3] = 0;            /* x-axis (consistent with sphere) */
  rdata_dot[4] = 0;
  rdata_dot[5] = 0;
  rdata_dot[6] = 0;            /* y-axis (cross of sphere x-axis and y-axis) */
  rdata_dot[7] = 0;
  rdata_dot[8] = 0;
  rdata_dot[9] = r_dot;        /* radius */

  /* set the Circle sensitivity */
  stat = EG_setGeometry_dot(circle, CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* set the Nodes sensitivities */
  rdata[0] = xcent[0];
  rdata[1] = xcent[1];
  rdata[2] = xcent[2] - r;

  rdata_dot[0] = xcent_dot[0];
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2] - r_dot;
  stat = EG_setGeometry_dot(enodes[0], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  rdata[0] = xcent[0];
  rdata[1] = xcent[1];
  rdata[2] = xcent[2] + r;

  rdata_dot[0] = xcent_dot[0];
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2] + r_dot;
  stat = EG_setGeometry_dot(enodes[1], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    /* get the Nodes from the Edge with the u == PI Circle */
    stat = EG_getTopology(ledges[3], &circle, &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the u == PI Circle data */
    rdata[0] = xcent[0]; /* center */
    rdata[1] = xcent[1];
    rdata[2] = xcent[2];
    rdata[3] = -1;       /* x-axis (reversed so sphere v is consistent with circle t) */
    rdata[4] =  0;
    rdata[5] =  0;
    rdata[6] =  0;       /* y-axis (cross of sphere x-axis and y-axis) */
    rdata[7] =  0;
    rdata[8] =  1;
    rdata[9] = r;        /* radius */

    rdata_dot[0] = xcent_dot[0]; /* center */
    rdata_dot[1] = xcent_dot[1];
    rdata_dot[2] = xcent_dot[2];
    rdata_dot[3] = 0;            /* x-axis */
    rdata_dot[4] = 0;
    rdata_dot[5] = 0;
    rdata_dot[6] = 0;            /* y-axis */
    rdata_dot[7] = 0;
    rdata_dot[8] = 0;
    rdata_dot[9] = r_dot;        /* radius */

    /* set the Circle sensitivity */
    stat = EG_setGeometry_dot(circle, CURVE, CIRCLE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

  }

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }
#ifdef PCURVE_SENSITIVITY
  EG_free(ivec);
#endif
  EG_free(rvec);

  return stat;
}


/* nodes that define each edge of a split CONE */
static const int iedgeNodesCone[5][2] = {{1,2},
                                         {0,1},
                                         {0,-1},
                                         {0,2},
                                         {2,1}};


/* edges that make up each face in a split CONE */
static const int ifaceEdgesCone[3][4] = {{ 0, 1, 2, 3},
                                         { 2, 1, 4, 3},
                                         { 0, 4, -1, -1}};


int
EG_makeSolidCone(egObject *context, int stypx, const double *data,
                 egObject **body)
{
  /* vertex, center of base, and radius */
  double xvert[3] = {data[0], data[1], data[2]};
  double xbase[3] = {data[3], data[4], data[5]};
  double r        = data[6];

  /* orderings taken from OCC resulting body */

  int esens[2][4] = { {SREVERSE, SREVERSE, SFORWARD, SFORWARD},
                      {SFORWARD, SFORWARD, SREVERSE, SREVERSE} };

  double pcurve[2][4][4] = { { /* Face 0 */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {   0,    0, 0, 1}, /* UMIN */
                              {   0,    0, 1, 0}, /* VMIN */
                              {  PI,    0, 0, 1}, /* UMAX */
                             },
                             { /* Face 1 */
                              {   0,    0, 1, 0}, /* VMIN */
                              {2*PI,    0, 0, 1}, /* UMAX */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {  PI,    0, 0, 1}, /* UMIN */
                             } };

  int      stat = EGADS_SUCCESS;
  int      i, j, outLevel, nface, oclass, mtype, lsens[1] = {SFORWARD}, *ivec=NULL;
  double   rdata[14], tdata[2], height, angle, vmax;
  double   dirx[3], diry[3], dirz[3], xyz[3], *rvec=NULL;
  ego      cone, plane, trimmed, circle, lines[2], nodes[3], enodes[2];
  ego      edges[5], ledges[8], loop, faces[3], shell, ref;
  objStack stack;

  outLevel = EG_outLevel(context);

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  dirz[0] = xbase[0]-xvert[0];
  dirz[1] = xbase[1]-xvert[1];
  dirz[2] = xbase[2]-xvert[2];

  height = sqrt( DOT(dirz,dirz) );

  if (height == 0.0) {
    printf(" EGADS Error: Zero height (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  dirz[0] /= height;
  dirz[1] /= height;
  dirz[2] /= height;

  EG_Ax2(dirz, dirx, diry);

  angle = atan(r/height);
  vmax = r/sin(angle);

  /* set vmax for the P-curves */
  pcurve[0][0][1] = vmax;
  pcurve[1][2][1] = vmax;

  /* create the Cone */
  rdata[ 0] = xvert[0]; /* center */
  rdata[ 1] = xvert[1];
  rdata[ 2] = xvert[2];
  rdata[ 3] = dirx[0];  /* x-axis */
  rdata[ 4] = dirx[1];
  rdata[ 5] = dirx[2];
  rdata[ 6] = diry[0];  /* y-axis */
  rdata[ 7] = diry[1];
  rdata[ 8] = diry[2];
  rdata[ 9] = dirz[0];  /* z-axis */
  rdata[10] = dirz[1];
  rdata[11] = dirz[2];
  rdata[12] = angle;    /* angle  */
  rdata[13] = 0;        /* radius */
  stat = EG_makeGeometry(context, SURFACE, CONICAL, NULL, NULL,
                         rdata, &cone);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, cone);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_getGeometry(cone, &oclass, &mtype, &ref, &ivec, &rvec);
  if ((stat != EGADS_SUCCESS) || (rvec == NULL)) goto cleanup;

  /* get the axes for the cone */
  dirx[0] = rvec[ 3];
  dirx[1] = rvec[ 4];
  dirx[2] = rvec[ 5];
  diry[0] = rvec[ 6];
  diry[1] = rvec[ 7];
  diry[2] = rvec[ 8];
  dirz[0] = rvec[ 9];
  dirz[1] = rvec[10];
  dirz[2] = rvec[11];

  /* create the Plane for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  stat = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                         rdata, &plane);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, plane);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Circle curve for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circle);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circle);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the Vertex */
  rdata[0] = xvert[0];
  rdata[1] = xvert[1];
  rdata[2] = xvert[2];
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the Circle */
  xyz[0] = rdata[0] = xbase[0] + dirx[0]*r;
  xyz[1] = rdata[1] = xbase[1] + dirx[1]*r;
  xyz[2] = rdata[2] = xbase[2] + dirx[2]*r;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Line from the vertex to the base */
  rdata[0] = xvert[0];
  rdata[1] = xvert[1];
  rdata[2] = xvert[2];
  rdata[3] = xyz[0] - xvert[0];
  rdata[4] = xyz[1] - xvert[1];
  rdata[5] = xyz[2] - xvert[2];
  stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                         rdata, &lines[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, lines[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    /* create the 2nd Node on the Circle */
    xyz[0] = rdata[0] = xbase[0] - dirx[0]*r;
    xyz[1] = rdata[1] = xbase[1] - dirx[1]*r;
    xyz[2] = rdata[2] = xbase[2] - dirx[2]*r;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the 2nd Line from the vertex to the base */
    rdata[0] = xvert[0];
    rdata[1] = xvert[1];
    rdata[2] = xvert[2];
    rdata[3] = xyz[0] - xvert[0];
    rdata[4] = xyz[1] - xvert[1];
    rdata[5] = xyz[2] - xvert[2];
    stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           rdata, &lines[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, lines[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;


    /* make 1st Edge on the Circle */
    tdata[0] = 0;
    tdata[1] = PI;

    enodes[0] = nodes[iedgeNodesCone[0][0]];
    enodes[1] = nodes[iedgeNodesCone[0][1]];

    stat = EG_makeTopology(context, circle, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the first Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[iedgeNodesCone[1][0]];
    enodes[1] = nodes[iedgeNodesCone[1][1]];

    stat = EG_makeTopology(context, lines[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Degenerate Edges on the vertex Node */
    tdata[0] = 0;
    tdata[1] = PI;

    enodes[0] = nodes[iedgeNodesCone[2][0]];

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, enodes, NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the second Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[iedgeNodesCone[3][0]];
    enodes[1] = nodes[iedgeNodesCone[3][1]];

    stat = EG_makeTopology(context, lines[1], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make 2nd Edge on the Circle */
    tdata[0] =   PI;
    tdata[1] = 2*PI;

    enodes[0] = nodes[iedgeNodesCone[4][0]];
    enodes[1] = nodes[iedgeNodesCone[4][1]];

    stat = EG_makeTopology(context, circle, EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[4]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Faces on the Cone */

    /* Edge and P-curves for Face 1 loop */
    for (j = 0; j < 4; j++) {
      ledges[j] = edges[ifaceEdgesCone[0][j]];

      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[0][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* trim the surface */
    rdata[0] = 0;
    rdata[1] = PI;
    rdata[2] = 0;
    rdata[3] = vmax;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, cone, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* re-make the Degenerate Edges on the vertex Node */
    tdata[0] = PI;
    tdata[1] = 2*PI;

    enodes[0] = nodes[0];

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, enodes, NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* Edge and P-curves for Face 2 loop */
     for (j = 0; j < 4; j++) {
      ledges[j] = edges[ifaceEdgesCone[1][j]];

      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[1][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* trim the surface */
    rdata[0] = PI;
    rdata[1] = 2*PI;
    rdata[2] = 0;
    rdata[3] = vmax;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, cone, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the second Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the base Face */
    //esens[1][0] == SFORWARD;
    //esens[1][1] == SFORWARD;

    ledges[0] = edges[ifaceEdgesCone[2][0]];
    ledges[1] = edges[ifaceEdgesCone[2][1]];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, plane, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 3;

  } else {

    /* make Edge on the Circle */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[1];

    stat = EG_makeTopology(context, circle, EDGE, ONENODE,
                           tdata, 1, enodes, NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the first Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[0];
    enodes[1] = nodes[1];

    stat = EG_makeTopology(context, lines[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Degenerate Edges on the vertex Node */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[0];

    stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                           tdata, 1, enodes, NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;


    /* make Faces on the Cone */
    ledges[0] = edges[0];
    ledges[1] = edges[1];
    ledges[2] = edges[2];
    ledges[3] = edges[1];

    /* modify the UMAX P-curve to the complete periodic */
    pcurve[0][3][0] = 2*PI;

    /* P-curves */
    for (j = 0; j < 4; j++) {
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[0][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, cone, LOOP, CLOSED,
                           NULL, 4, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, cone, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the base Face */
    // esens[1][0] = SFORWARD;

    ledges[0] = edges[0];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, plane, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 2;
  }

  /* make the shell */
  stat = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL,
                         nface, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* make the final body */
  stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                         &shell, NULL, body);
  if (stat != EGADS_SUCCESS) goto cleanup;


cleanup:
  EG_free(ivec);
  EG_free(rvec);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (%s)!\n", i, __func__);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  return stat;
}


int
EG_makeSolidCone_dot(int stypx, const double *data, const double *data_dot,
                       egObject *body)
{
  /* vertex, center of base, and radius */
  double xvert[3] = {data[0], data[1], data[2]};
  double xbase[3] = {data[3], data[4], data[5]};
  double r        = data[6];

  double xvert_dot[3] = {data_dot[0], data_dot[1], data_dot[2]};
  double xbase_dot[3] = {data_dot[3], data_dot[4], data_dot[5]};
  double r_dot        = data_dot[6];

#ifdef PCURVE_SENSITIVITY
  double pcurve[2][4][4] = { { /* Face 0 */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {   0,    0, 0, 1}, /* UMIN */
                              {   0,    0, 1, 0}, /* VMIN */
                              {  PI,    0, 0, 1}, /* UMAX */
                             },
                             { /* Face 1 */
                              {   0,    0, 1, 0}, /* VMIN */
                              {2*PI,    0, 0, 1}, /* UMAX */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {  PI,    0, 0, 1}, /* UMIN */
                             } };


  double pcurve_dot[2][4][4] = { { /* Face 0 */
                                  {   0, 0000, 0, 0}, /* VMAX */
                                  {   0,    0, 0, 0}, /* UMIN */
                                  {   0,    0, 0, 0}, /* VMIN */
                                  {   0,    0, 0, 0}, /* UMAX */
                                 },
                                 { /* Face 1 */
                                  {   0,    0, 0, 0}, /* VMIN */
                                  {   0,    0, 0, 0}, /* UMAX */
                                  {   0, 0000, 0, 0}, /* VMAX */
                                  {   0,    0, 0, 0}, /* UMIN */
                                 } };
#endif

  int    stat = EGADS_SUCCESS;
#ifdef PCURVE_SENSITIVITY
  int    j, *ivec=NULL;
#endif
  int    i, oclass, mtype, outLevel, *senses;
  int    nshell, nface, nloop, nedge, nnode;
  double dirx[3], dirx_dot[3], diry[3], diry_dot[3], dirz[3], dirz_dot[3];
  double rdata[14], rdata_dot[14], height, height_dot, height2;
  double angle, angle_dot, vmax, vmax_dot, xyz[3], xyz_dot[3];
  double rcone[14], rcone_dot[14], *rvec=NULL, *rvec_dot=NULL;
  ego    *faces, *shells, edges[5], nodes[3];
  ego    ref, circle, lines[2], trimmed, cone, plane, *enodes, *ledges, *loops;

  outLevel = EG_outLevel(body);

  dirz[0] = xbase[0]-xvert[0];
  dirz[1] = xbase[1]-xvert[1];
  dirz[2] = xbase[2]-xvert[2];

  dirz_dot[0] = xbase_dot[0]-xvert_dot[0];
  dirz_dot[1] = xbase_dot[1]-xvert_dot[1];
  dirz_dot[2] = xbase_dot[2]-xvert_dot[2];

  height = sqrt( DOT(dirz,dirz) );

  if (height == 0.0) {
    printf(" EGADS Error: Zero height (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  height_dot = (dirz[0]*dirz_dot[0] + dirz[1]*dirz_dot[1] +
                dirz[2]*dirz_dot[2])/height;

  height2 = height*height;

  dirz_dot[0] = dirz_dot[0]/height - dirz[0]*height_dot/height2;
  dirz_dot[1] = dirz_dot[1]/height - dirz[1]*height_dot/height2;
  dirz_dot[2] = dirz_dot[2]/height - dirz[2]*height_dot/height2;

  dirz[0] = dirz[0]/height;
  dirz[1] = dirz[1]/height;
  dirz[2] = dirz[2]/height;

  EG_Ax2_dot(dirz, dirz_dot,
             dirx, dirx_dot,
             diry, diry_dot);

  angle = atan(r/height);
  angle_dot = (-r*height_dot/height2 + r_dot/height)/(1 + pow(r/height,2));

  vmax = r/sin(angle);
  vmax_dot = r_dot/sin(angle) - angle_dot/tan(angle) * r/sin(angle) ;

#ifdef PCURVE_SENSITIVITY
  /* set vmax for the P-curves */
  pcurve[0][0][1] = vmax;
  pcurve[1][2][1] = vmax;

  pcurve_dot[0][0][1] = vmax_dot;
  pcurve_dot[1][2][1] = vmax_dot;
#endif

  /* the Cone data */
  rcone[ 0] = xvert[0]; /* center */
  rcone[ 1] = xvert[1];
  rcone[ 2] = xvert[2];
  rcone[ 3] = dirx[0];  /* x-axis */
  rcone[ 4] = dirx[1];
  rcone[ 5] = dirx[2];
  rcone[ 6] = diry[0];  /* y-axis */
  rcone[ 7] = diry[1];
  rcone[ 8] = diry[2];
  rcone[ 9] = dirz[0];  /* z-axis */
  rcone[10] = dirz[1];
  rcone[11] = dirz[2];
  rcone[12] = angle;    /* angle  */
  rcone[13] = 0;        /* radius */

  rcone_dot[ 0] = xvert_dot[0]; /* center */
  rcone_dot[ 1] = xvert_dot[1];
  rcone_dot[ 2] = xvert_dot[2];
  rcone_dot[ 3] = dirx_dot[0];  /* x-axis */
  rcone_dot[ 4] = dirx_dot[1];
  rcone_dot[ 5] = dirx_dot[2];
  rcone_dot[ 6] = diry_dot[0];  /* y-axis */
  rcone_dot[ 7] = diry_dot[1];
  rcone_dot[ 8] = diry_dot[2];
  rcone_dot[ 9] = dirz_dot[0];  /* z-axis */
  rcone_dot[10] = dirz_dot[1];
  rcone_dot[11] = dirz_dot[2];
  rcone_dot[12] = angle_dot;    /* angle  */
  rcone_dot[13] = 0;        /* radius */

  /* get the Shell from the Body */
  stat = EG_getTopology(body, &ref, &oclass, &mtype,
                        rdata, &nshell, &shells, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, BODY, SOLIDBODY);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if ( nshell != 1 ) {
    if (outLevel > 0)
      printf(" EGADS Error: body must have one SHELL (%s)!\n", __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* get the Faces from the Shell */
  stat = EG_getTopology(shells[0], &ref, &oclass, &mtype,
                        rdata, &nface, &faces, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_expectClassType(oclass, mtype, SHELL, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    if ( nface != 3 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 3 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set surfaces sensitivity and extract edges from the loop in the cone faces */
    for (i = 0; i < 2; i++) {

      /* get the Loop from the Face */
      stat = EG_getTopology(faces[i], &trimmed, &oclass, &mtype,
                            rdata, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nloop != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
                 i+1, nloop, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      /* set trim the surface sensitivity */
      if (i == 0) {
        rdata[0] = 0;
        rdata[1] = PI;
        rdata[2] = 0;
        rdata[3] = vmax;
      } else {
        rdata[0] = PI;
        rdata[1] = 2*PI;
        rdata[2] = 0;
        rdata[3] = vmax;
      }

      rdata_dot[0] = 0;
      rdata_dot[1] = 0;
      rdata_dot[2] = 0;
      rdata_dot[3] = vmax_dot;

      stat = EG_setGeometry_dot(trimmed, SURFACE, TRIMMED, NULL, rdata, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Sphere from the Trimmed surface */
      stat = EG_getGeometry(trimmed, &oclass, &mtype, &cone, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* set the cone sensitivities */
      stat = EG_setGeometry_dot(cone, SURFACE, CONICAL, NULL, rcone, rcone_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Edges from the Loop */
      stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                            rdata, &nedge, &ledges, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nedge != 4) {
        if (outLevel > 0)
          printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
                 i+1, nedge, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      edges[ifaceEdgesCone[i][0]] = ledges[0];
      edges[ifaceEdgesCone[i][1]] = ledges[1];
      edges[ifaceEdgesCone[i][2]] = ledges[2];
      edges[ifaceEdgesCone[i][3]] = ledges[3];

#ifdef PCURVE_SENSITIVITY
      /* set PCURVE sensitivities */
      for (j = 0; j < 4; j++) {
        stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, pcurve[i][j], pcurve_dot[i][j]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
#endif
    }

    /* extract all the nodes from the edges */
    for (i = 0; i < 5; i++) {

      /* get the Nodes from the Edge */
      stat = EG_getTopology(edges[i], &ref, &oclass, &mtype,
                            rdata, &nnode, &enodes, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (iedgeNodesCone[i][1] >= 0) {

        stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* extract the nodes */
        nodes[iedgeNodesCone[i][0]] = enodes[0];
        nodes[iedgeNodesCone[i][1]] = enodes[1];
      } else {

        stat = EG_expectClassType(oclass, mtype, EDGE, DEGENERATE);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* extract the node */
        nodes[iedgeNodesCone[i][0]] = enodes[0];
      }
    }

    stat = EG_getGeometry_dot(cone, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the cone */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    /* get the Plane from the 3rd Face */
    stat = EG_getTopology(faces[2], &plane, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract geometry from edges */
    stat = EG_getTopology(edges[0], &circle, &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[1], &lines[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[3], &lines[1], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* the 2nd Node on the Circle */
    xyz[0] = rdata[0] = xbase[0] - dirx[0]*r;
    xyz[1] = rdata[1] = xbase[1] - dirx[1]*r;
    xyz[2] = rdata[2] = xbase[2] - dirx[2]*r;

    xyz_dot[0] = rdata_dot[0] = xbase_dot[0] - dirx_dot[0]*r - dirx[0]*r_dot;
    xyz_dot[1] = rdata_dot[1] = xbase_dot[1] - dirx_dot[1]*r - dirx[1]*r_dot;
    xyz_dot[2] = rdata_dot[2] = xbase_dot[2] - dirx_dot[2]*r - dirx[2]*r_dot;

    stat = EG_setGeometry_dot(nodes[2], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* the 2nd Line from the vertex to the base */
    rdata[0] = xvert[0];
    rdata[1] = xvert[1];
    rdata[2] = xvert[2];
    rdata[3] = xyz[0] - xvert[0];
    rdata[4] = xyz[1] - xvert[1];
    rdata[5] = xyz[2] - xvert[2];

    rdata_dot[0] = xvert_dot[0];
    rdata_dot[1] = xvert_dot[1];
    rdata_dot[2] = xvert_dot[2];
    rdata_dot[3] = xyz_dot[0] - xvert_dot[0];
    rdata_dot[4] = xyz_dot[1] - xvert_dot[1];
    rdata_dot[5] = xyz_dot[2] - xvert_dot[2];

    stat = EG_setGeometry_dot(lines[1], CURVE, LINE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

  } else {

    if ( nface != 2 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 2 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* get the Loop from the Face */
    stat = EG_getTopology(faces[0], &cone, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nloop != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
               1, nloop, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* get the PCurves from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
    /* modify the UMAX P-curve to the complete periodic */
    pcurve[0][3][0] = 2*PI;

    /* set PCURVE sensitivities */
    for (j = 0; j < 4; j++) {
      stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, pcurve[0][j], pcurve_dot[0][j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }
#endif

    /* set the cone sensitivities */
    stat = EG_setGeometry_dot(cone, SURFACE, CONICAL, NULL, rcone, rcone_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_getGeometry_dot(cone, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the cone */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    /* get the Edges from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if (nedge != 4) {
      if (outLevel > 0)
        printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
               1, nedge, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    edges[0] = ledges[0];
    edges[1] = ledges[1];
    edges[2] = ledges[2];
    edges[3] = ledges[3];

    /* get the Line and Nodes from the Edge */
    stat = EG_getTopology(edges[1], &lines[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract the nodes */
    nodes[0] = enodes[0];
    nodes[1] = enodes[1];

    /* get the Plane from the 2rd Face */
    stat = EG_getTopology(faces[1], &plane, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* get the circle */
    stat = EG_getTopology(edges[0], &circle, &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* the Plane for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0]; /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];

  rdata_dot[0] = xbase_dot[0]; /* center */
  rdata_dot[1] = xbase_dot[1];
  rdata_dot[2] = xbase_dot[2];
  rdata_dot[3] = dirx_dot[0]; /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];

  stat = EG_setGeometry_dot(plane, SURFACE, PLANE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Circle curve for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */

  rdata_dot[0] = xbase_dot[0]; /* center */
  rdata_dot[1] = xbase_dot[1];
  rdata_dot[2] = xbase_dot[2];
  rdata_dot[3] = dirx_dot[0];  /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];
  rdata_dot[9] = r_dot;        /* radius */

  stat = EG_setGeometry_dot(circle, CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Node on the Vertex */
  rdata[0] = xvert[0];
  rdata[1] = xvert[1];
  rdata[2] = xvert[2];

  rdata_dot[0] = xvert_dot[0];
  rdata_dot[1] = xvert_dot[1];
  rdata_dot[2] = xvert_dot[2];

  stat = EG_setGeometry_dot(nodes[0], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the first Node on the Circle */
  xyz[0] = rdata[0] = xbase[0] + dirx[0]*r;
  xyz[1] = rdata[1] = xbase[1] + dirx[1]*r;
  xyz[2] = rdata[2] = xbase[2] + dirx[2]*r;

  xyz_dot[0] = rdata_dot[0] = xbase_dot[0] + dirx_dot[0]*r + dirx[0]*r_dot;
  xyz_dot[1] = rdata_dot[1] = xbase_dot[1] + dirx_dot[1]*r + dirx[1]*r_dot;
  xyz_dot[2] = rdata_dot[2] = xbase_dot[2] + dirx_dot[2]*r + dirx[2]*r_dot;

  stat = EG_setGeometry_dot(nodes[1], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Line from the vertex to the base */
  rdata[0] = xvert[0];
  rdata[1] = xvert[1];
  rdata[2] = xvert[2];
  rdata[3] = xyz[0] - xvert[0];
  rdata[4] = xyz[1] - xvert[1];
  rdata[5] = xyz[2] - xvert[2];

  rdata_dot[0] = xvert_dot[0];
  rdata_dot[1] = xvert_dot[1];
  rdata_dot[2] = xvert_dot[2];
  rdata_dot[3] = xyz_dot[0] - xvert_dot[0];
  rdata_dot[4] = xyz_dot[1] - xvert_dot[1];
  rdata_dot[5] = xyz_dot[2] - xvert_dot[2];

  stat = EG_setGeometry_dot(lines[0], CURVE, LINE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
#ifdef PCURVE_SENSITIVITY
  EG_free(ivec);
#endif
  EG_free(rvec);
  EG_free(rvec_dot);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  return stat;
}


/* nodes that define each edge of a split CYLINDER */
static const int iedgeNodesCyl[6][2] = {{1,3},   /* 0 - top Circle */
                                        {0,1},   /* 1 - u == 0 Line */
                                        {0,2},   /* 2 - base Circle */
                                        {2,3},   /* 3 - u == PI Line */
                                        {3,1},   /* 4 - top Circle */
                                        {2,0}};  /* 5 - base Circle */


/* edges that make up each face in a split CYLINDER */
static const int ifaceEdgesCyl[4][4] = {{ 0, 1,  2,  3 },
                                        { 5, 1,  4,  3 },
                                        { 2, 5, -1, -1 },
                                        { 0, 4, -1, -1 }};


int
EG_makeSolidCylinder(egObject *context, int stypx, const double *data,
                     egObject **body)
{
  /* center of base, center of top, and radius */
  double xbase[3] = {data[0], data[1], data[2]};
  double xcent[3] = {data[3], data[4], data[5]};
  double r        =  data[6];

  /* orderings taken from OCC resulting body */

  int esens[2][4] = { {SREVERSE, SREVERSE, SFORWARD, SFORWARD},
                      {SFORWARD, SFORWARD, SREVERSE, SREVERSE} };

  double pcurve[2][4][4] = { { /* Face 0 */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {   0,    0, 0, 1}, /* UMIN */
                              {   0,    0, 1, 0}, /* VMIN */
                              {  PI,    0, 0, 1}, /* UMAX */
                             },
                             { /* Face 1 */
                              {   0,    0, 1, 0}, /* VMIN */
                              {2*PI,    0, 0, 1}, /* UMAX */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {  PI,    0, 0, 1}, /* UMIN */
                             } };

  int      stat = EGADS_SUCCESS;
  int      i, j, outLevel, nface, oclass, mtype, lsens[1] = {SFORWARD}, *ivec=NULL;
  double   rdata[13], tdata[2], height, vmax;
  double   dirx[3], diry[3], dirz[3], xyz0[3], xyz1[3], *rvec=NULL;
  ego      cylinder, planes[2], trimmed, circles[2], lines[2], nodes[4];
  ego      enodes[2], edges[6], ledges[8], loop, faces[4], shell, ref;
  objStack stack;

  outLevel = EG_outLevel(context);

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  dirz[0] = xcent[0]-xbase[0];
  dirz[1] = xcent[1]-xbase[1];
  dirz[2] = xcent[2]-xbase[2];

  height = sqrt( DOT(dirz,dirz) );
  if (height == 0.0) {
    printf(" EGADS Error: Zero height (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  dirz[0] /= height;
  dirz[1] /= height;
  dirz[2] /= height;

  EG_Ax2(dirz, dirx, diry);

  vmax = height;

  /* set vmax for the P-curves */
  pcurve[0][0][1] = vmax;
  pcurve[1][2][1] = vmax;

  /* create the Cylinder */
  rdata[ 0] = xbase[0]; /* center */
  rdata[ 1] = xbase[1];
  rdata[ 2] = xbase[2];
  rdata[ 3] = dirx[0];  /* x-axis */
  rdata[ 4] = dirx[1];
  rdata[ 5] = dirx[2];
  rdata[ 6] = diry[0];  /* y-axis */
  rdata[ 7] = diry[1];
  rdata[ 8] = diry[2];
  rdata[ 9] = dirz[0];  /* z-axis */
  rdata[10] = dirz[1];
  rdata[11] = dirz[2];
  rdata[12] = r;        /* radius */
  stat = EG_makeGeometry(context, SURFACE, CYLINDRICAL, NULL, NULL,
                         rdata, &cylinder);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, cylinder);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_getGeometry(cylinder, &oclass, &mtype, &ref, &ivec, &rvec);
  if ((stat != EGADS_SUCCESS) || (rvec == NULL)) goto cleanup;

  /* get the axes for the cylinder */
  dirx[0] = rvec[ 3];
  dirx[1] = rvec[ 4];
  dirx[2] = rvec[ 5];
  diry[0] = rvec[ 6];
  diry[1] = rvec[ 7];
  diry[2] = rvec[ 8];
  dirz[0] = rvec[ 9];
  dirz[1] = rvec[10];
  dirz[2] = rvec[11];

  /* create the Plane for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0]; /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  stat = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                         rdata, &planes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, planes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Plane for the top */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  stat = EG_makeGeometry(context, SURFACE, PLANE, NULL, NULL,
                         rdata, &planes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, planes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Circle curve for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Circle curve for the top */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circles[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circles[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the base Circle */
  xyz0[0] = rdata[0] = xbase[0] + dirx[0]*r;
  xyz0[1] = rdata[1] = xbase[1] + dirx[1]*r;
  xyz0[2] = rdata[2] = xbase[2] + dirx[2]*r;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the Node on the top Circle */
  xyz1[0] = rdata[0] = xcent[0] + dirx[0]*r;
  xyz1[1] = rdata[1] = xcent[1] + dirx[1]*r;
  xyz1[2] = rdata[2] = xcent[2] + dirx[2]*r;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the u == 0 Line from the base to top */
  rdata[0] = xyz0[0];
  rdata[1] = xyz0[1];
  rdata[2] = xyz0[2];
  rdata[3] = xyz1[0] - xyz0[0];
  rdata[4] = xyz1[1] - xyz0[1];
  rdata[5] = xyz1[2] - xyz0[2];
  stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                         rdata, &lines[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, lines[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    /* create the 2nd Node on the base Circle */
    xyz0[0] = rdata[0] = xbase[0] - dirx[0]*r;
    xyz0[1] = rdata[1] = xbase[1] - dirx[1]*r;
    xyz0[2] = rdata[2] = xbase[2] - dirx[2]*r;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the 2nd Node on the top Circle */
    xyz1[0] = rdata[0] = xcent[0] - dirx[0]*r;
    xyz1[1] = rdata[1] = xcent[1] - dirx[1]*r;
    xyz1[2] = rdata[2] = xcent[2] - dirx[2]*r;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the u == PI Line from the base to the top */
    rdata[0] = xyz0[0];
    rdata[1] = xyz0[1];
    rdata[2] = xyz0[2];
    rdata[3] = xyz1[0] - xyz0[0];
    rdata[4] = xyz1[1] - xyz0[1];
    rdata[5] = xyz1[2] - xyz0[2];
    stat = EG_makeGeometry(context, CURVE, LINE, NULL, NULL,
                           rdata, &lines[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, lines[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make 1st Edge on the top Circle */
    tdata[0] = 0;
    tdata[1] = PI;

    enodes[0] = nodes[iedgeNodesCyl[0][0]];
    enodes[1] = nodes[iedgeNodesCyl[0][1]];

    stat = EG_makeTopology(context, circles[1], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the u == 0 Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[iedgeNodesCyl[1][0]];
    enodes[1] = nodes[iedgeNodesCyl[1][1]];

    stat = EG_makeTopology(context, lines[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edges on the base Circle */
    tdata[0] = 0;
    tdata[1] = PI;

    enodes[0] = nodes[iedgeNodesCyl[2][0]];
    enodes[1] = nodes[iedgeNodesCyl[2][1]];

    stat = EG_makeTopology(context, circles[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Edge on the u == PI Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[iedgeNodesCyl[3][0]];
    enodes[1] = nodes[iedgeNodesCyl[3][1]];

    stat = EG_makeTopology(context, lines[1], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make 2nd Edge on the top Circle */
    tdata[0] =   PI;
    tdata[1] = 2*PI;

    enodes[0] = nodes[iedgeNodesCyl[4][0]];
    enodes[1] = nodes[iedgeNodesCyl[4][1]];

    stat = EG_makeTopology(context, circles[1], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[4]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[4]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make 2nd Edge on the base Circle */
    tdata[0] =   PI;
    tdata[1] = 2*PI;

    enodes[0] = nodes[iedgeNodesCyl[5][0]];
    enodes[1] = nodes[iedgeNodesCyl[5][1]];

    stat = EG_makeTopology(context, circles[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[5]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[5]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Faces on the Cylinder */

    /* Edge and P-curves for Face 1 loop */
    for (j = 0; j < 4; j++) {
      ledges[j] = edges[ifaceEdgesCyl[0][j]];

      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[0][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* trim the surface */
    rdata[0] = 0;
    rdata[1] = PI;
    rdata[2] = 0;
    rdata[3] = vmax;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, cylinder, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the first Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* Edge and P-curves for Face 2 loop */
    for (j = 0; j < 4; j++) {
      ledges[j] = edges[ifaceEdgesCyl[1][j]];

      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[1][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* trim the surface */
    rdata[0] = PI;
    rdata[1] = 2*PI;
    rdata[2] = 0;
    rdata[3] = vmax;
    stat = EG_makeGeometry(context, SURFACE, TRIMMED, cylinder, NULL,
                           rdata, &trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, trimmed);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the second Face */
    stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                           NULL, 4, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the top Face */
    //esens[1][0] == SFORWARD;
    //esens[1][1] == SFORWARD;

    ledges[0] = edges[ifaceEdgesCyl[3][0]];
    ledges[1] = edges[ifaceEdgesCyl[3][1]];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, planes[1], FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the base Face */
    //esens[0][0] == SREVERSE;
    //esens[0][1] == SREVERSE;

    ledges[0] = edges[ifaceEdgesCyl[2][0]];
    ledges[1] = edges[ifaceEdgesCyl[2][1]];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 2, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, planes[0], FACE, SREVERSE,
                           NULL, 1, &loop, lsens, &faces[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 4;

  } else {

    /* make Edge on the top Circle */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[1];

    stat = EG_makeTopology(context, circles[1], EDGE, ONENODE,
                           tdata, 1, enodes, NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the u == 0 Edge on the first Line */
    tdata[0] = 0;
    tdata[1] = vmax;

    enodes[0] = nodes[0];
    enodes[1] = nodes[1];

    stat = EG_makeTopology(context, lines[0], EDGE, TWONODE,
                           tdata, 2, enodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Edge on the base Circle */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[0];

    stat = EG_makeTopology(context, circles[0], EDGE, ONENODE,
                           tdata, 1, enodes, NULL, &edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Faces on the Cylinder */
    ledges[0] = edges[0];
    ledges[1] = edges[1];
    ledges[2] = edges[2];
    ledges[3] = edges[1];

    /* modify the UMAX P-curve to the complete periodic */
    pcurve[0][3][0] = 2*PI;

    /* P-curves */
    for (j = 0; j < 4; j++) {
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[0][j], &ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, ledges[4+j]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* make the Loop and the cylinder Face */
    stat = EG_makeTopology(context, cylinder, LOOP, CLOSED,
                           NULL, 4, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, cylinder, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the top Face */
    // esens[1][0] = SFORWARD;

    ledges[0] = edges[0];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, planes[1], FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the base Face */
    // esens[0][0] == SREVERSE;

    ledges[0] = edges[2];

    stat = EG_makeTopology(context, NULL, LOOP, CLOSED,
                           NULL, 1, ledges, esens[0], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, planes[0], FACE, SREVERSE,
                           NULL, 1, &loop, lsens, &faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 3;
  }

  /* make the shell */
  stat = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL,
                         nface, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* make the final body */
  stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                         &shell, NULL, body);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
  EG_free(ivec);
  EG_free(rvec);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (%s)!\n", i, __func__);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  return stat;
}


int
EG_makeSolidCylinder_dot(int stypx, const double *data, const double *data_dot,
                       egObject *body)
{
  /* center of base, center of top, and radius */
  double xbase[3] = {data[0], data[1], data[2]};
  double xcent[3] = {data[3], data[4], data[5]};
  double r        =  data[6];

  double xbase_dot[3] = {data_dot[0], data_dot[1], data_dot[2]};
  double xcent_dot[3] = {data_dot[3], data_dot[4], data_dot[5]};
  double r_dot        = data_dot[6];

#ifdef PCURVE_SENSITIVITY
  double pcurve[2][4][4] = { { /* Face 0 */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {   0,    0, 0, 1}, /* UMIN */
                              {   0,    0, 1, 0}, /* VMIN */
                              {  PI,    0, 0, 1}, /* UMAX */
                             },
                             { /* Face 1 */
                              {   0,    0, 1, 0}, /* VMIN */
                              {2*PI,    0, 0, 1}, /* UMAX */
                              {   0, 0000, 1, 0}, /* VMAX */
                              {  PI,    0, 0, 1}, /* UMIN */
                             } };


  double pcurve_dot[2][4][4] = { { /* Face 0 */
                                  {   0, 0000, 0, 0}, /* VMAX */
                                  {   0,    0, 0, 0}, /* UMIN */
                                  {   0,    0, 0, 0}, /* VMIN */
                                  {   0,    0, 0, 0}, /* UMAX */
                                 },
                                 { /* Face 1 */
                                  {   0,    0, 0, 0}, /* VMIN */
                                  {   0,    0, 0, 0}, /* UMAX */
                                  {   0, 0000, 0, 0}, /* VMAX */
                                  {   0,    0, 0, 0}, /* UMIN */
                                 } };
#endif

  int    stat = EGADS_SUCCESS;
#ifdef PCURVE_SENSITIVITY
  int    j, *ivec=NULL;
#endif
  int    i, oclass, mtype, outLevel, *senses;
  int    nshell, nface, nloop, nedge, nnode;
  double dirx[3], dirx_dot[3], diry[3], diry_dot[3], dirz[3], dirz_dot[3];
  double rdata[14], rdata_dot[14], height, height_dot, height2;
  double vmax, vmax_dot, xyz0[3], xyz0_dot[3], xyz1[3], xyz1_dot[3];
  double rcyl[13], rcyl_dot[13], *rvec=NULL, *rvec_dot=NULL;
  ego    *faces, *shells, *enodes, *ledges, *loops, edges[6], nodes[4];
  ego    ref, circles[2], lines[2], trimmed, cylinder, planes[2];

  outLevel = EG_outLevel(body);

  dirz[0] = xcent[0]-xbase[0];
  dirz[1] = xcent[1]-xbase[1];
  dirz[2] = xcent[2]-xbase[2];

  dirz_dot[0] = xcent_dot[0]-xbase_dot[0];
  dirz_dot[1] = xcent_dot[1]-xbase_dot[1];
  dirz_dot[2] = xcent_dot[2]-xbase_dot[2];

  height = sqrt( DOT(dirz,dirz) );
  if (height == 0.0) {
    printf(" EGADS Error: Zero height (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  height_dot = (dirz[0]*dirz_dot[0] + dirz[1]*dirz_dot[1] +
                dirz[2]*dirz_dot[2])/height;

  height2 = height*height;

  dirz_dot[0] = dirz_dot[0]/height - dirz[0]*height_dot/height2;
  dirz_dot[1] = dirz_dot[1]/height - dirz[1]*height_dot/height2;
  dirz_dot[2] = dirz_dot[2]/height - dirz[2]*height_dot/height2;

  dirz[0] = dirz[0]/height;
  dirz[1] = dirz[1]/height;
  dirz[2] = dirz[2]/height;

  EG_Ax2_dot(dirz, dirz_dot,
             dirx, dirx_dot,
             diry, diry_dot);

  vmax     = height;
  vmax_dot = height_dot;

#ifdef PCURVE_SENSITIVITY
  /* set vmax for the P-curves */
  pcurve[0][0][1] = vmax;
  pcurve[1][2][1] = vmax;

  pcurve_dot[0][0][1] = vmax_dot;
  pcurve_dot[1][2][1] = vmax_dot;
#endif

  /* the Cylinder data */
  rcyl[ 0] = xbase[0]; /* center */
  rcyl[ 1] = xbase[1];
  rcyl[ 2] = xbase[2];
  rcyl[ 3] = dirx[0];  /* x-axis */
  rcyl[ 4] = dirx[1];
  rcyl[ 5] = dirx[2];
  rcyl[ 6] = diry[0];  /* y-axis */
  rcyl[ 7] = diry[1];
  rcyl[ 8] = diry[2];
  rcyl[ 9] = dirz[0];  /* z-axis */
  rcyl[10] = dirz[1];
  rcyl[11] = dirz[2];
  rcyl[12] = r;        /* radius */

  rcyl_dot[ 0] = xbase_dot[0]; /* center */
  rcyl_dot[ 1] = xbase_dot[1];
  rcyl_dot[ 2] = xbase_dot[2];
  rcyl_dot[ 3] = dirx_dot[0];  /* x-axis */
  rcyl_dot[ 4] = dirx_dot[1];
  rcyl_dot[ 5] = dirx_dot[2];
  rcyl_dot[ 6] = diry_dot[0];  /* y-axis */
  rcyl_dot[ 7] = diry_dot[1];
  rcyl_dot[ 8] = diry_dot[2];
  rcyl_dot[ 9] = dirz_dot[0];  /* z-axis */
  rcyl_dot[10] = dirz_dot[1];
  rcyl_dot[11] = dirz_dot[2];
  rcyl_dot[12] = r_dot;        /* radius */

  /* get the Shell from the Body */
  stat = EG_getTopology(body, &ref, &oclass, &mtype,
                        rdata, &nshell, &shells, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, BODY, SOLIDBODY);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if ( nshell != 1 ) {
    if (outLevel > 0)
      printf(" EGADS Error: body must have one SHELL (%s)!\n", __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* get the Faces from the Shell */
  stat = EG_getTopology(shells[0], &ref, &oclass, &mtype,
                        rdata, &nface, &faces, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, SHELL, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    if (nface != 4) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 4 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set surfaces sensitivity and extract edges from the loop in the cylinder faces */
    for (i = 0; i < 2; i++) {

      /* get the Loop from the Face */
      stat = EG_getTopology(faces[i], &trimmed, &oclass, &mtype,
                            rdata, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (nloop != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
                 i+1, nloop, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      /* set trim the surface sensitivity */
      if (i == 0) {
        rdata[0] = 0;
        rdata[1] = PI;
        rdata[2] = 0;
        rdata[3] = vmax;
      } else {
        rdata[0] = PI;
        rdata[1] = 2*PI;
        rdata[2] = 0;
        rdata[3] = vmax;
      }

      rdata_dot[0] = 0;
      rdata_dot[1] = 0;
      rdata_dot[2] = 0;
      rdata_dot[3] = vmax_dot;

      stat = EG_setGeometry_dot(trimmed, SURFACE, TRIMMED, NULL, rdata, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Cylinder from the Trimmed surface */
      stat = EG_getGeometry(trimmed, &oclass, &mtype, &cylinder, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* set the cone sensitivities */
      stat = EG_setGeometry_dot(cylinder, SURFACE, CYLINDRICAL, NULL, rcyl, rcyl_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Edges from the Loop */
      stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                            rdata, &nedge, &ledges, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nedge != 4) {
        if (outLevel > 0)
          printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
                 i+1, nedge, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      edges[ifaceEdgesCyl[i][0]] = ledges[0];
      edges[ifaceEdgesCyl[i][1]] = ledges[1];
      edges[ifaceEdgesCyl[i][2]] = ledges[2];
      edges[ifaceEdgesCyl[i][3]] = ledges[3];

#ifdef PCURVE_SENSITIVITY
      /* set PCURVE sensitivities */
      for (j = 0; j < 4; j++) {
        stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, pcurve[i][j], pcurve_dot[i][j]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
#endif
    }

    /* extract all the nodes from the edges */
    for (i = 0; i < 6; i++) {

      /* get the Nodes from the Edge */
      stat = EG_getTopology(edges[i], &ref, &oclass, &mtype,
                            rdata, &nnode, &enodes, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* extract the nodes */
      nodes[iedgeNodesCyl[i][0]] = enodes[0];
      nodes[iedgeNodesCyl[i][1]] = enodes[1];
    }

    stat = EG_getGeometry_dot(cylinder, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the cylinder */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    EG_free(rvec);     rvec     = NULL;
    EG_free(rvec_dot); rvec_dot = NULL;

    /* get the Planes from the last Faces */
    stat = EG_getTopology(faces[2], &planes[1], &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(faces[3], &planes[0], &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract geometry from edges */
    stat = EG_getTopology(edges[0], &circles[1], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[1], &lines[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[2], &circles[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[3], &lines[1], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* the 2nd Node on the base Circle */
    xyz0[0] = rdata[0] = xbase[0] - dirx[0]*r;
    xyz0[1] = rdata[1] = xbase[1] - dirx[1]*r;
    xyz0[2] = rdata[2] = xbase[2] - dirx[2]*r;

    xyz0_dot[0] = rdata_dot[0] = xbase_dot[0] - dirx_dot[0]*r - dirx[0]*r_dot;
    xyz0_dot[1] = rdata_dot[1] = xbase_dot[1] - dirx_dot[1]*r - dirx[1]*r_dot;
    xyz0_dot[2] = rdata_dot[2] = xbase_dot[2] - dirx_dot[2]*r - dirx[2]*r_dot;

    stat = EG_setGeometry_dot(nodes[2], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* the 2nd Node on the top Circle */
    xyz1[0] = rdata[0] = xcent[0] - dirx[0]*r;
    xyz1[1] = rdata[1] = xcent[1] - dirx[1]*r;
    xyz1[2] = rdata[2] = xcent[2] - dirx[2]*r;

    xyz1_dot[0] = rdata_dot[0] = xcent_dot[0] - dirx_dot[0]*r - dirx[0]*r_dot;
    xyz1_dot[1] = rdata_dot[1] = xcent_dot[1] - dirx_dot[1]*r - dirx[1]*r_dot;
    xyz1_dot[2] = rdata_dot[2] = xcent_dot[2] - dirx_dot[2]*r - dirx[2]*r_dot;

    stat = EG_setGeometry_dot(nodes[3], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* the u == PI Line */
    rdata[0] = xyz0[0];
    rdata[1] = xyz0[1];
    rdata[2] = xyz0[2];
    rdata[3] = xyz1[0] - xyz0[0];
    rdata[4] = xyz1[1] - xyz0[1];
    rdata[5] = xyz1[2] - xyz0[2];

    rdata_dot[0] = xyz0_dot[0];
    rdata_dot[1] = xyz0_dot[1];
    rdata_dot[2] = xyz0_dot[2];
    rdata_dot[3] = xyz1_dot[0] - xyz0_dot[0];
    rdata_dot[4] = xyz1_dot[1] - xyz0_dot[1];
    rdata_dot[5] = xyz1_dot[2] - xyz0_dot[2];

    stat = EG_setGeometry_dot(lines[1], CURVE, LINE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

  } else {

    if (nface != 3) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 3 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* get the Loop from the Face */
    stat = EG_getTopology(faces[0], &cylinder, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nloop != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
               1, nloop, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set the cylinder sensitivities */
    stat = EG_setGeometry_dot(cylinder, SURFACE, CYLINDRICAL, NULL, rcyl, rcyl_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_getGeometry_dot(cylinder, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the cylinder */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    EG_free(rvec);     rvec     = NULL;
    EG_free(rvec_dot); rvec_dot = NULL;

    /* get the Edges from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if (nedge != 4) {
      if (outLevel > 0)
        printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
               1, nedge, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    edges[0] = ledges[0];
    edges[1] = ledges[1];
    edges[2] = ledges[2];
    edges[3] = ledges[3];

#ifdef PCURVE_SENSITIVITY
    /* set PCURVE sensitivities */
    rdata_dot[0] = rdata_dot[1] = rdata_dot[2] = rdata_dot[3] = 0;
    for (j = 0; j < 4; j++) {
      /* generally calling EG_getGeometry to get rvec for EG_setGeometry_dot is not correct,
       * but here it is ok because the direction was specified as a unit vector */
      stat = EG_getGeometry(ledges[4+j], &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, rvec, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;
      EG_free(ivec); ivec=NULL;
      EG_free(rvec); rvec=NULL;
    }
#endif

    /* get the Line and Nodes from the Edge */
    stat = EG_getTopology(edges[1], &lines[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract the nodes */
    nodes[0] = enodes[0];
    nodes[1] = enodes[1];

    /* get the Planes from the last Faces */
    stat = EG_getTopology(faces[1], &planes[1], &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(faces[2], &planes[0], &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* get the Circles from Edges */
    stat = EG_getTopology(edges[0], &circles[1], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, ONENODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(edges[2], &circles[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, ONENODE);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* the Plane for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0]; /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];

  rdata_dot[0] = xbase_dot[0]; /* center */
  rdata_dot[1] = xbase_dot[1];
  rdata_dot[2] = xbase_dot[2];
  rdata_dot[3] = dirx_dot[0]; /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];

  stat = EG_setGeometry_dot(planes[0], SURFACE, PLANE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Plane for the top */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0]; /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];

  rdata_dot[0] = xcent_dot[0]; /* center */
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2];
  rdata_dot[3] = dirx_dot[0]; /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];

  stat = EG_setGeometry_dot(planes[1], SURFACE, PLANE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Circle curve for the base */
  rdata[0] = xbase[0]; /* center */
  rdata[1] = xbase[1];
  rdata[2] = xbase[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */

  rdata_dot[0] = xbase_dot[0]; /* center */
  rdata_dot[1] = xbase_dot[1];
  rdata_dot[2] = xbase_dot[2];
  rdata_dot[3] = dirx_dot[0];  /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];
  rdata_dot[9] = r_dot;        /* radius */

  stat = EG_setGeometry_dot(circles[0], CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Circle curve for the top */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = r;        /* radius */

  rdata_dot[0] = xcent_dot[0]; /* center */
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2];
  rdata_dot[3] = dirx_dot[0];  /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];
  rdata_dot[9] = r_dot;        /* radius */

  stat = EG_setGeometry_dot(circles[1], CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Node on the base Circle */
  xyz0[0] = rdata[0] = xbase[0] + dirx[0]*r;
  xyz0[1] = rdata[1] = xbase[1] + dirx[1]*r;
  xyz0[2] = rdata[2] = xbase[2] + dirx[2]*r;

  xyz0_dot[0] = rdata_dot[0] = xbase_dot[0] + dirx_dot[0]*r + dirx[0]*r_dot;
  xyz0_dot[1] = rdata_dot[1] = xbase_dot[1] + dirx_dot[1]*r + dirx[1]*r_dot;
  xyz0_dot[2] = rdata_dot[2] = xbase_dot[2] + dirx_dot[2]*r + dirx[2]*r_dot;

  stat = EG_setGeometry_dot(nodes[0], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the first Node on the Circle */
  xyz1[0] = rdata[0] = xcent[0] + dirx[0]*r;
  xyz1[1] = rdata[1] = xcent[1] + dirx[1]*r;
  xyz1[2] = rdata[2] = xcent[2] + dirx[2]*r;

  xyz1_dot[0] = rdata_dot[0] = xcent_dot[0] + dirx_dot[0]*r + dirx[0]*r_dot;
  xyz1_dot[1] = rdata_dot[1] = xcent_dot[1] + dirx_dot[1]*r + dirx[1]*r_dot;
  xyz1_dot[2] = rdata_dot[2] = xcent_dot[2] + dirx_dot[2]*r + dirx[2]*r_dot;

  stat = EG_setGeometry_dot(nodes[1], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the u == 0 Line */
  rdata[0] = xyz0[0];
  rdata[1] = xyz0[1];
  rdata[2] = xyz0[2];
  rdata[3] = xyz1[0] - xyz0[0];
  rdata[4] = xyz1[1] - xyz0[1];
  rdata[5] = xyz1[2] - xyz0[2];

  rdata_dot[0] = xyz0_dot[0];
  rdata_dot[1] = xyz0_dot[1];
  rdata_dot[2] = xyz0_dot[2];
  rdata_dot[3] = xyz1_dot[0] - xyz0_dot[0];
  rdata_dot[4] = xyz1_dot[1] - xyz0_dot[1];
  rdata_dot[5] = xyz1_dot[2] - xyz0_dot[2];

  stat = EG_setGeometry_dot(lines[0], CURVE, LINE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
  EG_free(rvec);
  EG_free(rvec_dot);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  return stat;
}


/* curve associated with each edge of a split TORUS */
static const int iedgeCurveTorus[8] = {0,  /* 0 - u == 0 */
                                       0,  /* 1 */
                                       1,  /* 2 - v == 0 */
                                       1,  /* 3 */
                                       2,  /* 4 - u == PI */
                                       2,  /* 5 */
                                       3,  /* 6 - v == PI */
                                       3}; /* 7 */


/* nodes that define each edge of a split TORUS */
static const int iedgeNodesTorus[8][2] = {{0,1},   /* 0 - u == 0  v 0-PI   */
                                          {1,0},   /* 1 - u == 0  v PI-2PI */
                                          {0,2},   /* 2 - v == 0  u 0-PI   */
                                          {2,0},   /* 3 - v == 0  u PI-2PI */
                                          {2,3},   /* 4 - u == PI v 0-PI   */
                                          {3,2},   /* 5 - u == PI v PI-2PI */
                                          {1,3},   /* 6 - v == PI u 0-PI   */
                                          {3,1}};  /* 7 - v == PI u PI-2PI */


/* edges that make up each face in a split TORUS */
static const int ifaceEdgesTorus[4][4] = {{ 3, 0, 7, 4},   /* Face 0 */
                                          { 6, 0, 2, 4},   /* Face 1 */
                                          { 2, 1, 6, 5},   /* Face 2 */
                                          { 7, 1, 3, 5}};  /* Face 3 */


/* trmming parameter spaces for split TORUS */
static const double ftrimTorus[4][4] = { {PI, 2*PI,  0,   PI},   /* Face 0 */
                                         { 0,   PI,  0,   PI},   /* Face 1 */
                                         { 0,   PI, PI, 2*PI},   /* Face 2 */
                                         {PI, 2*PI, PI, 2*PI} }; /* Face 3 */


int
EG_makeSolidTorus(egObject *context, int stypx, const double *data,
                  egObject **body)
{
  /* center of base, center of top, major radius, and minor radius */
  double xcent[3] = {data[0], data[1], data[2]};
  double dirz[3]  = {data[3], data[4], data[5]};
  double majr     = data[6];
  double minr     = data[7];

  /* orderings taken from OCC resulting body */

  int esens[4][4] = { {SFORWARD, SFORWARD, SREVERSE, SREVERSE},   /* Face 0 */
                      {SREVERSE, SREVERSE, SFORWARD, SFORWARD},   /* Face 1 */
                      {SREVERSE, SREVERSE, SFORWARD, SFORWARD},   /* Face 2 */
                      {SFORWARD, SFORWARD, SREVERSE, SREVERSE} }; /* Face 3 */

  const double pcurve[4][4][4] = { { /* Face 0 */
                                    {   0,    0, 1, 0}, /* VMIN */
                                    {2*PI,    0, 0, 1}, /* UMAX */
                                    {   0,   PI, 1, 0}, /* VMAX */
                                    {  PI,    0, 0, 1}, /* UMIN */
                                   },
                                   { /* Face 1 */
                                    {   0,   PI, 1, 0}, /* VMAX */
                                    {   0,    0, 0, 1}, /* UMIN */
                                    {   0,    0, 1, 0}, /* VMIN */
                                    {  PI,    0, 0, 1}, /* UMAX */
                                   },
                                   { /* Face 2 */
                                    {   0, 2*PI, 1, 0}, /* VMAX */
                                    {   0,    0, 0, 1}, /* UMIN */
                                    {   0,   PI, 1, 0}, /* VMIN */
                                    {  PI,    0, 0, 1}, /* UMAX */
                                   },
                                   { /* Face 3 */
                                    {   0,   PI, 1, 0}, /* VMIN */
                                    {2*PI,    0, 0, 1}, /* UMAX */
                                    {   0, 2*PI, 1, 0}, /* VMAX */
                                    {  PI,    0, 0, 1}, /* UMIN */
                                   } };

  int      stat = EGADS_SUCCESS;
  int      i, j, outLevel, nface, oclass, mtype, *ivec=NULL;
  int      lsens[1] = {SFORWARD};
  double   rdata[14], tdata[2], R, r, height;
  double   dirx[3], diry[3], *rvec=NULL;
  ego      torus, trimmed, circles[4], nodes[4], enodes[2], edges[8], ledges[8];
  ego      loop, faces[4], shell, ref;
  objStack stack;

  outLevel = EG_outLevel(context);

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  R = majr + minr;
  r = majr - minr;

  height = sqrt( DOT(dirz,dirz) );
  if (height == 0.0) {
    printf(" EGADS Error: Zero magnitude axis (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  dirz[0] /= height;
  dirz[1] /= height;
  dirz[2] /= height;

  EG_Ax2(dirz, dirx, diry);

  /* create the Torus */
  rdata[ 0] = xcent[0]; /* center */
  rdata[ 1] = xcent[1];
  rdata[ 2] = xcent[2];
  rdata[ 3] = dirx[0];  /* x-axis */
  rdata[ 4] = dirx[1];
  rdata[ 5] = dirx[2];
  rdata[ 6] = diry[0];  /* y-axis */
  rdata[ 7] = diry[1];
  rdata[ 8] = diry[2];
  rdata[ 9] = dirz[0];  /* z-axis */
  rdata[10] = dirz[1];
  rdata[11] = dirz[2];
  rdata[12] = majr;     /* major radius */
  rdata[13] = minr;     /* minor radius */
  stat = EG_makeGeometry(context, SURFACE, TOROIDAL, NULL, NULL,
                         rdata, &torus);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, torus);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_getGeometry(torus, &oclass, &mtype, &ref, &ivec, &rvec);
  if ((stat != EGADS_SUCCESS) || (rvec == NULL)) goto cleanup;

  /* get the axes for the cylinder */
  dirx[0] = rvec[ 3];
  dirx[1] = rvec[ 4];
  dirx[2] = rvec[ 5];
  diry[0] = rvec[ 6];
  diry[1] = rvec[ 7];
  diry[2] = rvec[ 8];
  dirz[0] = rvec[ 9];
  dirz[1] = rvec[10];
  dirz[2] = rvec[11];

  /* create the minor u == 0 Circle curve */
  rdata[0] = xcent[0] + dirx[0]*majr; /* center */
  rdata[1] = xcent[1] + dirx[1]*majr;
  rdata[2] = xcent[2] + dirx[2]*majr;
  rdata[3] = dirx[0];   /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = dirz[0];   /* y-axis */
  rdata[7] = dirz[1];
  rdata[8] = dirz[2];
  rdata[9] = minr;      /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circles[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* create the major v == 0 Circle curve */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = R;        /* radius */
  stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                         rdata, &circles[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, circles[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;


  /* create the outer Node on both Circles */
  rdata[0] = xcent[0] + dirx[0]*R;
  rdata[1] = xcent[1] + dirx[1]*R;
  rdata[2] = xcent[2] + dirx[2]*R;
  stat = EG_makeTopology(context, NULL, NODE, 0,
                         rdata, 0, NULL, NULL, &nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, nodes[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    /* create the minor u == PI Circle curve */
    rdata[0] = xcent[0] - dirx[0]*majr; /* center */
    rdata[1] = xcent[1] - dirx[1]*majr;
    rdata[2] = xcent[2] - dirx[2]*majr;
    rdata[3] = -dirx[0];   /* x-axis */
    rdata[4] = -dirx[1];
    rdata[5] = -dirx[2];
    rdata[6] = dirz[0];    /* y-axis */
    rdata[7] = dirz[1];
    rdata[8] = dirz[2];
    rdata[9] = minr;       /* radius */
    stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           rdata, &circles[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, circles[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the major v == PI Circle curve */
    rdata[0] = xcent[0]; /* center */
    rdata[1] = xcent[1];
    rdata[2] = xcent[2];
    rdata[3] = dirx[0];  /* x-axis */
    rdata[4] = dirx[1];
    rdata[5] = dirx[2];
    rdata[6] = diry[0];  /* y-axis */
    rdata[7] = diry[1];
    rdata[8] = diry[2];
    rdata[9] = r;        /* radius */
    stat = EG_makeGeometry(context, CURVE, CIRCLE, NULL, NULL,
                           rdata, &circles[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, circles[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the inner Node on the u == 0 Circle */
    rdata[0] = xcent[0] + dirx[0]*r;
    rdata[1] = xcent[1] + dirx[1]*r;
    rdata[2] = xcent[2] + dirx[2]*r;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the outer Node on the u == PI Circle */
    rdata[0] = xcent[0] - dirx[0]*R;
    rdata[1] = xcent[1] - dirx[1]*R;
    rdata[2] = xcent[2] - dirx[2]*R;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the inner Node on the u == PI Circle */
    rdata[0] = xcent[0] - dirx[0]*r;
    rdata[1] = xcent[1] - dirx[1]*r;
    rdata[2] = xcent[2] - dirx[2]*r;
    stat = EG_makeTopology(context, NULL, NODE, 0,
                           rdata, 0, NULL, NULL, &nodes[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, nodes[3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make all the Edges */
    for (i = 0; i < 8; i++) {
      if (i % 2 == 0) {
        tdata[0] = 0;
        tdata[1] = PI;
      } else {
        tdata[0] = PI;
        tdata[1] = 2*PI;
      }

      enodes[0] = nodes[iedgeNodesTorus[i][0]];
      enodes[1] = nodes[iedgeNodesTorus[i][1]];

      stat = EG_makeTopology(context, circles[iedgeCurveTorus[i]], EDGE, TWONODE,
                             tdata, 2, enodes, NULL, &edges[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, edges[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    /* make all the Faces */
    for (i = 0; i < 4; i++) {

      for (j = 0; j < 4; j++) {
        ledges[j] = edges[ifaceEdgesTorus[i][j]];

        /* P-curve */
        stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, pcurve[i][j], &ledges[4+j]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, ledges[4+j]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      /* trim the surface */
      stat = EG_makeGeometry(context, SURFACE, TRIMMED, torus, NULL,
                             ftrimTorus[i], &trimmed);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, trimmed);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the Loop and the Face */
      stat = EG_makeTopology(context, trimmed, LOOP, CLOSED,
                             NULL, 4, ledges, esens[i], &loop);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, loop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_makeTopology(context, trimmed, FACE, SFORWARD,
                             NULL, 1, &loop, lsens, &faces[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, faces[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    nface = 4;

  } else {

    /* make Edge on the u == 0 Circle */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[0];

    stat = EG_makeTopology(context, circles[0], EDGE, ONENODE,
                           tdata, 1, enodes, NULL, &edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Edge on the v == 0 Circle */
    tdata[0] = 0;
    tdata[1] = 2*PI;

    enodes[0] = nodes[0];

    stat = EG_makeTopology(context, circles[1], EDGE, ONENODE,
                           tdata, 1, enodes, NULL, &edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, edges[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make Faces on the Cylinder */
    ledges[0] = edges[1];
    ledges[1] = edges[0];
    ledges[2] = edges[1];
    ledges[3] = edges[0];

    /* P-curve */
    rdata[0] = 0.;    rdata[1] =  2*PI;   /* v == 2*PI VMAX */
    rdata[2] = 1.;    rdata[3] =  0.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &ledges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, ledges[4+0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 0.;    rdata[1] = 0.;      /* u == 0 UMIN   */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &ledges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, ledges[4+1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 0.;    rdata[1] = 0.;      /* v == 0 VMIN */
    rdata[2] = 1.;    rdata[3] = 0.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &ledges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, ledges[4+2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    rdata[0] = 2*PI;  rdata[1] = 0.;      /* u == 2*PI UMAX  */
    rdata[2] = 0.;    rdata[3] = 1.;
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, rdata, &ledges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, ledges[4+3]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* make the Loop and the cylinder Face */
    stat = EG_makeTopology(context, torus, LOOP, CLOSED,
                           NULL, 4, ledges, esens[1], &loop);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_makeTopology(context, torus, FACE, SFORWARD,
                           NULL, 1, &loop, lsens, &faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, faces[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    nface = 1;
  }

  /* make the shell */
  stat = EG_makeTopology(context, NULL, SHELL, CLOSED, NULL,
                         nface, faces, NULL, &shell);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(&stack, shell);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* make the final body */
  stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                         &shell, NULL, body);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
  EG_free(ivec);
  EG_free(rvec);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (%s)!\n", i, __func__);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  return stat;
}


int
EG_makeSolidTorus_dot(int stypx, const double *data, const double *data_dot,
                       egObject *body)
{
  /* center of base, center of top, major radius, and minor radius */
  double xcent[3] = {data[0], data[1], data[2]};
  double dirz[3]  = {data[3], data[4], data[5]};
  double majr     = data[6];
  double minr     = data[7];

  double xcent_dot[3] = {data_dot[0], data_dot[1], data_dot[2]};
  double dirz_dot[3]  = {data_dot[3], data_dot[4], data_dot[5]};
  double majr_dot     = data_dot[6];
  double minr_dot     = data_dot[7];

  int    stat = EGADS_SUCCESS;
#ifdef PCURVE_SENSITIVITY
  int      j, *ivec=NULL;
#endif
  int    i, oclass, mtype, outLevel, *senses;
  int    nshell, nface, nloop, nedge, nnode;
  double dirx[3], dirx_dot[3], diry[3], diry_dot[3], r, r_dot, R, R_dot;
  double rdata[14], rdata_dot[14], height, height_dot, height2;
  double rtorus[14], rtorus_dot[14], *rvec=NULL, *rvec_dot=NULL;
  ego    *faces, *shells, edges[8], nodes[4];
  ego    ref, circles[4], trimmed, torus, *enodes, *ledges, *loops;

  outLevel = EG_outLevel(body);

  R = majr + minr;
  r = majr - minr;

  R_dot = majr_dot + minr_dot;
  r_dot = majr_dot - minr_dot;

  height = sqrt( DOT(dirz,dirz) );

  if (height == 0.0) {
    printf(" EGADS Error: Zero height (%s)!\n", __func__);
    return EGADS_GEOMERR;
  }

  height_dot = (dirz[0]*dirz_dot[0] + dirz[1]*dirz_dot[1] +
                dirz[2]*dirz_dot[2])/height;

  height2 = height*height;

  dirz_dot[0] = dirz_dot[0]/height - dirz[0]*height_dot/height2;
  dirz_dot[1] = dirz_dot[1]/height - dirz[1]*height_dot/height2;
  dirz_dot[2] = dirz_dot[2]/height - dirz[2]*height_dot/height2;

  dirz[0] = dirz[0]/height;
  dirz[1] = dirz[1]/height;
  dirz[2] = dirz[2]/height;

  EG_Ax2_dot(dirz, dirz_dot,
             dirx, dirx_dot,
             diry, diry_dot);

  /* the Torus data */
  rtorus[ 0] = xcent[0]; /* center */
  rtorus[ 1] = xcent[1];
  rtorus[ 2] = xcent[2];
  rtorus[ 3] = dirx[0];  /* x-axis */
  rtorus[ 4] = dirx[1];
  rtorus[ 5] = dirx[2];
  rtorus[ 6] = diry[0];  /* y-axis */
  rtorus[ 7] = diry[1];
  rtorus[ 8] = diry[2];
  rtorus[ 9] = dirz[0];  /* z-axis */
  rtorus[10] = dirz[1];
  rtorus[11] = dirz[2];
  rtorus[12] = majr;     /* major radius */
  rtorus[13] = minr;     /* minor radius */

  rtorus_dot[ 0] = xcent_dot[0]; /* center */
  rtorus_dot[ 1] = xcent_dot[1];
  rtorus_dot[ 2] = xcent_dot[2];
  rtorus_dot[ 3] = dirx_dot[0];  /* x-axis */
  rtorus_dot[ 4] = dirx_dot[1];
  rtorus_dot[ 5] = dirx_dot[2];
  rtorus_dot[ 6] = diry_dot[0];  /* y-axis */
  rtorus_dot[ 7] = diry_dot[1];
  rtorus_dot[ 8] = diry_dot[2];
  rtorus_dot[ 9] = dirz_dot[0];  /* z-axis */
  rtorus_dot[10] = dirz_dot[1];
  rtorus_dot[11] = dirz_dot[2];
  rtorus_dot[12] = majr_dot;     /* major radius */
  rtorus_dot[13] = minr_dot;     /* minor radius */

  /* get the Shell from the Body */
  stat = EG_getTopology(body, &ref, &oclass, &mtype,
                        rdata, &nshell, &shells, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, BODY, SOLIDBODY);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if ( nshell != 1 ) {
    if (outLevel > 0)
      printf(" EGADS Error: body must have one SHELL (%s)!\n", __func__);
    stat = EGADS_TOPOERR;
    goto cleanup;
  }

  /* get the Faces from the Shell */
  stat = EG_getTopology(shells[0], &ref, &oclass, &mtype,
                        rdata, &nface, &faces, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_expectClassType(oclass, mtype, SHELL, CLOSED);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (stypx > 0) { /* split periodic */

    if ( nface != 4 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 4 Faces (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set surfaces sensitivity and extract edges from the loop in the torus faces */
    for (i = 0; i < 4; i++) {

      /* get the Loop from the Face */
      stat = EG_getTopology(faces[i], &trimmed, &oclass, &mtype,
                            rdata, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (nloop != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
                 i+1, nloop, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      /* set trim the surface sensitivity */
      rdata_dot[0] = 0;
      rdata_dot[1] = 0;
      rdata_dot[2] = 0;
      rdata_dot[3] = 0;

      stat = EG_setGeometry_dot(trimmed, SURFACE, TRIMMED, NULL, ftrimTorus[i], rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Sphere from the Trimmed surface */
      stat = EG_getGeometry(trimmed, &oclass, &mtype, &torus, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* set the torus sensitivities */
      stat = EG_setGeometry_dot(torus, SURFACE, TOROIDAL, NULL, rtorus, rtorus_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the Edges from the Loop */
      stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                            rdata, &nedge, &ledges, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nedge != 4) {
        if (outLevel > 0)
          printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
                 i+1, nedge, __func__);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      edges[ifaceEdgesTorus[i][0]] = ledges[0];
      edges[ifaceEdgesTorus[i][1]] = ledges[1];
      edges[ifaceEdgesTorus[i][2]] = ledges[2];
      edges[ifaceEdgesTorus[i][3]] = ledges[3];

#ifdef PCURVE_SENSITIVITY
      /* set PCURVE sensitivities */
      rdata_dot[0] = rdata_dot[1] = rdata_dot[2] = rdata_dot[3] = 0;
      for (j = 0; j < 4; j++) {
        /* generally calling EG_getGeometry to get rvec for EG_setGeometry_dot is not correct,
         * but here it is ok because the direction was specified as a unit vector */
        stat = EG_getGeometry(ledges[4+j], &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, rvec, rdata_dot);
        if (stat != EGADS_SUCCESS) goto cleanup;
        EG_free(ivec); ivec=NULL;
        EG_free(rvec); rvec=NULL;
      }
#endif
    }

    /* extract all the nodes from the edges */
    for (i = 0; i < 8; i++) {

      /* get the Nodes from the Edge */
      stat = EG_getTopology(edges[i], &circles[iedgeCurveTorus[i]], &oclass, &mtype,
                            rdata, &nnode, &enodes, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_expectClassType(oclass, mtype, EDGE, TWONODE);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* extract the nodes */
      nodes[iedgeNodesTorus[i][0]] = enodes[0];
      nodes[iedgeNodesTorus[i][1]] = enodes[1];
    }

    stat = EG_getGeometry_dot(torus, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the torus */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    /* the minor u == PI Circle curve */
    rdata[0] = xcent[0] - dirx[0]*majr; /* center */
    rdata[1] = xcent[1] - dirx[1]*majr;
    rdata[2] = xcent[2] - dirx[2]*majr;
    rdata[3] = -dirx[0];   /* x-axis */
    rdata[4] = -dirx[1];
    rdata[5] = -dirx[2];
    rdata[6] = dirz[0];   /* y-axis */
    rdata[7] = dirz[1];
    rdata[8] = dirz[2];
    rdata[9] = minr;       /* radius */

    rdata_dot[0] = xcent_dot[0] - dirx_dot[0]*majr - dirx[0]*majr_dot; /* center */
    rdata_dot[1] = xcent_dot[1] - dirx_dot[1]*majr - dirx[1]*majr_dot;
    rdata_dot[2] = xcent_dot[2] - dirx_dot[2]*majr - dirx[2]*majr_dot;
    rdata_dot[3] = -dirx_dot[0];   /* x-axis */
    rdata_dot[4] = -dirx_dot[1];
    rdata_dot[5] = -dirx_dot[2];
    rdata_dot[6] = dirz_dot[0];   /* y-axis */
    rdata_dot[7] = dirz_dot[1];
    rdata_dot[8] = dirz_dot[2];
    rdata_dot[9] = minr_dot;       /* radius */

    stat = EG_setGeometry_dot(circles[2], CURVE, CIRCLE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the major v == PI Circle curve */
    rdata[0] = xcent[0]; /* center */
    rdata[1] = xcent[1];
    rdata[2] = xcent[2];
    rdata[3] = dirx[0];  /* x-axis */
    rdata[4] = dirx[1];
    rdata[5] = dirx[2];
    rdata[6] = diry[0];  /* y-axis */
    rdata[7] = diry[1];
    rdata[8] = diry[2];
    rdata[9] = r;        /* radius */

    rdata_dot[0] = xcent_dot[0]; /* center */
    rdata_dot[1] = xcent_dot[1];
    rdata_dot[2] = xcent_dot[2];
    rdata_dot[3] = dirx_dot[0];  /* x-axis */
    rdata_dot[4] = dirx_dot[1];
    rdata_dot[5] = dirx_dot[2];
    rdata_dot[6] = diry_dot[0];  /* y-axis */
    rdata_dot[7] = diry_dot[1];
    rdata_dot[8] = diry_dot[2];
    rdata_dot[9] = r_dot;        /* radius */

    stat = EG_setGeometry_dot(circles[3], CURVE, CIRCLE, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the inner Node on the u == 0 Circle */
    rdata[0] = xcent[0] + dirx[0]*r;
    rdata[1] = xcent[1] + dirx[1]*r;
    rdata[2] = xcent[2] + dirx[2]*r;

    rdata_dot[0] = xcent_dot[0] + dirx_dot[0]*r + dirx[0]*r_dot;
    rdata_dot[1] = xcent_dot[1] + dirx_dot[1]*r + dirx[1]*r_dot;
    rdata_dot[2] = xcent_dot[2] + dirx_dot[2]*r + dirx[2]*r_dot;

    stat = EG_setGeometry_dot(nodes[1], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the outer Node on the u == PI Circle */
    rdata[0] = xcent[0] - dirx[0]*R;
    rdata[1] = xcent[1] - dirx[1]*R;
    rdata[2] = xcent[2] - dirx[2]*R;

    rdata_dot[0] = xcent_dot[0] - dirx_dot[0]*R - dirx[0]*R_dot;
    rdata_dot[1] = xcent_dot[1] - dirx_dot[1]*R - dirx[1]*R_dot;
    rdata_dot[2] = xcent_dot[2] - dirx_dot[2]*R - dirx[2]*R_dot;

    stat = EG_setGeometry_dot(nodes[2], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* create the inner Node on the u == PI Circle */
    rdata[0] = xcent[0] - dirx[0]*r;
    rdata[1] = xcent[1] - dirx[1]*r;
    rdata[2] = xcent[2] - dirx[2]*r;

    rdata_dot[0] = xcent_dot[0] - dirx_dot[0]*r - dirx[0]*r_dot;
    rdata_dot[1] = xcent_dot[1] - dirx_dot[1]*r - dirx[1]*r_dot;
    rdata_dot[2] = xcent_dot[2] - dirx_dot[2]*r - dirx[2]*r_dot;

    stat = EG_setGeometry_dot(nodes[3], NODE, 0, NULL, rdata, rdata_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;

  } else {

    if ( nface != 1 ) {
      if (outLevel > 0)
        printf(" EGADS Error: SHELL (with %d Faces) must have 1 Face (%s)!\n",
               nface, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* get the Loop from the Face */
    stat = EG_getTopology(faces[0], &torus, &oclass, &mtype,
                          rdata, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nloop != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: Face %d (with %d Loops) must have 1 Loop (%s)!\n",
               1, nloop, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    /* set the torus sensitivities */
    stat = EG_setGeometry_dot(torus, SURFACE, TOROIDAL, NULL, rtorus, rtorus_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_getGeometry_dot(torus, &rvec, &rvec_dot);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if ((rvec == NULL) || (rvec_dot == NULL)) goto cleanup;

    /* get the axes for the cylinder */
    dirx[0] = rvec[ 3];
    dirx[1] = rvec[ 4];
    dirx[2] = rvec[ 5];
    diry[0] = rvec[ 6];
    diry[1] = rvec[ 7];
    diry[2] = rvec[ 8];
    dirz[0] = rvec[ 9];
    dirz[1] = rvec[10];
    dirz[2] = rvec[11];

    dirx_dot[0] = rvec_dot[ 3];
    dirx_dot[1] = rvec_dot[ 4];
    dirx_dot[2] = rvec_dot[ 5];
    diry_dot[0] = rvec_dot[ 6];
    diry_dot[1] = rvec_dot[ 7];
    diry_dot[2] = rvec_dot[ 8];
    dirz_dot[0] = rvec_dot[ 9];
    dirz_dot[1] = rvec_dot[10];
    dirz_dot[2] = rvec_dot[11];

    /* get the Edges from the Loop */
    stat = EG_getTopology(loops[0], &ref, &oclass, &mtype,
                          rdata, &nedge, &ledges, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, LOOP, CLOSED);
    if (stat != EGADS_SUCCESS) goto cleanup;
    if (nedge != 4) {
      if (outLevel > 0)
        printf(" EGADS Error: Loops for Face %d (with %d Edges) must have 4 Edges (%s)!\n",
               1, nedge, __func__);
      stat = EGADS_TOPOERR;
      goto cleanup;
    }

    edges[1] = ledges[0];
    edges[0] = ledges[1];
    edges[1] = ledges[2];
    edges[0] = ledges[3];

#ifdef PCURVE_SENSITIVITY
    /* set PCURVE sensitivities */
    rdata_dot[0] = rdata_dot[1] = rdata_dot[2] = rdata_dot[3] = 0;
    for (j = 0; j < 4; j++) {
      /* generally calling EG_getGeometry to get rvec for EG_setGeometry_dot is not correct,
       * but here it is ok because the direction was specified as a unit vector */
      stat = EG_getGeometry(ledges[4+j], &oclass, &mtype, &ref, &ivec, &rvec);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_setGeometry_dot(ledges[4+j], PCURVE, LINE, NULL, rvec, rdata_dot);
      if (stat != EGADS_SUCCESS) goto cleanup;
      EG_free(ivec); ivec=NULL;
      EG_free(rvec); rvec=NULL;
    }
#endif

    /* get the Circles and Node from the Edge */
    stat = EG_getTopology(edges[0], &circles[0], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, ONENODE);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* extract the node */
    nodes[0] = enodes[0];

    /* get the Circles from Edge */
    stat = EG_getTopology(edges[1], &circles[1], &oclass, &mtype,
                          rdata, &nnode, &enodes, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_expectClassType(oclass, mtype, EDGE, ONENODE);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* the minor u == 0 Circle curve */
  rdata[0] = xcent[0] + dirx[0]*majr; /* center */
  rdata[1] = xcent[1] + dirx[1]*majr;
  rdata[2] = xcent[2] + dirx[2]*majr;
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = dirz[0];  /* y-axis */
  rdata[7] = dirz[1];
  rdata[8] = dirz[2];
  rdata[9] = minr;     /* radius */

  rdata_dot[0] = xcent_dot[0] + dirx_dot[0]*majr + dirx[0]*majr_dot; /* center */
  rdata_dot[1] = xcent_dot[1] + dirx_dot[1]*majr + dirx[1]*majr_dot;
  rdata_dot[2] = xcent_dot[2] + dirx_dot[2]*majr + dirx[2]*majr_dot;
  rdata_dot[3] = dirx_dot[0];  /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = dirz_dot[0];  /* y-axis */
  rdata_dot[7] = dirz_dot[1];
  rdata_dot[8] = dirz_dot[2];
  rdata_dot[9] = minr_dot;     /* radius */

  stat = EG_setGeometry_dot(circles[0], CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the major v == 0 Circle curve */
  rdata[0] = xcent[0]; /* center */
  rdata[1] = xcent[1];
  rdata[2] = xcent[2];
  rdata[3] = dirx[0];  /* x-axis */
  rdata[4] = dirx[1];
  rdata[5] = dirx[2];
  rdata[6] = diry[0];  /* y-axis */
  rdata[7] = diry[1];
  rdata[8] = diry[2];
  rdata[9] = R;        /* radius */

  rdata_dot[0] = xcent_dot[0]; /* center */
  rdata_dot[1] = xcent_dot[1];
  rdata_dot[2] = xcent_dot[2];
  rdata_dot[3] = dirx_dot[0];  /* x-axis */
  rdata_dot[4] = dirx_dot[1];
  rdata_dot[5] = dirx_dot[2];
  rdata_dot[6] = diry_dot[0];  /* y-axis */
  rdata_dot[7] = diry_dot[1];
  rdata_dot[8] = diry_dot[2];
  rdata_dot[9] = R_dot;        /* radius */

  stat = EG_setGeometry_dot(circles[1], CURVE, CIRCLE, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* the Node on the base Circle */
  rdata[0] = xcent[0] + dirx[0]*R;
  rdata[1] = xcent[1] + dirx[1]*R;
  rdata[2] = xcent[2] + dirx[2]*R;

  rdata_dot[0] = xcent_dot[0] + dirx_dot[0]*R + dirx[0]*R_dot;
  rdata_dot[1] = xcent_dot[1] + dirx_dot[1]*R + dirx[1]*R_dot;
  rdata_dot[2] = xcent_dot[2] + dirx_dot[2]*R + dirx[2]*R_dot;

  stat = EG_setGeometry_dot(nodes[0], NODE, 0, NULL, rdata, rdata_dot);
  if (stat != EGADS_SUCCESS) goto cleanup;

cleanup:
#ifdef PCURVE_SENSITIVITY
  EG_free(ivec);
#endif
  EG_free(rvec);
  EG_free(rvec_dot);

  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (%s)!\n", stat, __func__);
  }

  return stat;
}
