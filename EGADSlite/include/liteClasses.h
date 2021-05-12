#ifndef EGADSLCLASSES_H
#define EGADSLCLASSES_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADSlite Object Header
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


typedef struct {
  egObject *ref;                  /* reference object or NULL */
  int      *header;
  double   *data;
} liteGeometry;


typedef struct {
  double xyz[3];
  double tol;
} liteNode;


typedef struct {
  egObject *curve;                /* curve object */
  egObject *nodes[2];             /* pointer to ego nodes */
  double   trange[2];
  double   bbox[6];               /* bounding box */
  double   tol;
} liteEdge;


typedef struct {
  egObject *surface;              /* associated non-planar surface
                                     will have pcurves after edges (nonNULL) */
  int      nedges;                /* number of edges */
  egObject **edges;               /* edge objects (*2 if surface is nonNULL) */
  int      *senses;               /* sense for each edge */
  double   bbox[6];               /* bounding box */
} liteLoop;


typedef struct {
  egObject *surface;              /* surface object */
  int      nloops;                /* number of loops */
  egObject **loops;               /* loop objects */
  int      *senses;               /* outer/inner for each loop */
  double   urange[2];
  double   vrange[2];
  double   bbox[6];               /* bounding box */
  double   tol;
} liteFace;


typedef struct {
  int       nfaces;               /* number of faces */
  egObject **faces;               /* face objects */
  double   bbox[6];               /* bounding box */
} liteShell;


typedef struct {
  int       nobjs;                /* number in the map */
  egObject **objs;                /* vector of egos for map */
} liteMap;


typedef struct {
  liteMap pcurves;
  liteMap curves;
  liteMap surfaces;
  liteMap nodes;
  liteMap edges;
  liteMap loops;
  liteMap faces;
  liteMap shells;
  int     *senses;                /* shell outer/inner (solids) */
  double   bbox[6];               /* bounding box */
} liteBody;


typedef struct {
  int      nbody;                 /* number of bodies */
  egObject **bodies;              /* vector of pointers to bodies */
  double   bbox[6];               /* bounding box */
} liteModel;

#endif
