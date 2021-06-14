#ifndef EGADSCLASSES_H
#define EGADSCLASSES_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             C++/OpenCASCADE Object Header
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egadsOCC.h"
#include "Surreal/SurrealS.h"

class egadsPCurve
{
public:
  Handle(Geom2d_Curve) handle;
  egObject             *ref;
  int                  topFlg;
  int                  dataLen;
  int                  *header;
  double               *data;
  SurrealS<1>          *data_dot;
  double               trange[2];
};


class egadsCurve
{
public:
  Handle(Geom_Curve) handle;
  egObject           *ref;
  int                topFlg;
  int                dataLen;
  int                *header;
  double             *data;
  SurrealS<1>        *data_dot;
  double             trange[2];
};


class egadsSurface
{
public:
  Handle(Geom_Surface) handle;
  egObject             *ref;
  int                  topFlg;
  int                  dataLen;
  int                  *header;
  double               *data;
  SurrealS<1>          *data_dot;
  double               urange[2];
  double               vrange[2];
};


class egadsBox
{
public:
  int    filled;
  double box[6];
};


class egadsNode
{
public:
  TopoDS_Vertex node;
  double        xyz[3];
  egadsBox      bbox;
  int           filled;
  SurrealS<1>   xyz_dot[3];
};


class egadsEdge
{
public:
  TopoDS_Edge edge;
  egObject    *curve;                   // curve object
  egObject    *nodes[2];                // pointer to ego nodes
  int         topFlg;
  double      trange[2];
  egadsBox    bbox;
  int         filled;
  SurrealS<1> trange_dot[2];
};


class egadsLoop
{
public:
  TopoDS_Wire loop;
  egObject    *surface;                 // associated non-planar surface
                                        // will have pcurves after edges (nonNULL)
  int         nedges;                   // number of edges
  egObject    **edges;                  // edge objects (*2 if surface is nonNULL)
  int         *senses;                  // sense for each edge
  int         topFlg;
  egadsBox    bbox;
};


class egadsFace
{
public:
  TopoDS_Face face;
  egObject    *surface;                 // surface object
  int         nloops;                   // number of loops
  egObject    **loops;                  // loop objects
  int         *senses;                  // outer/inner for each loop
  int         topFlg;
  double      urange[2];
  double      vrange[2];
  egadsBox    bbox;
};


class egadsShell
{
public:
  TopoDS_Shell shell;
  int          nfaces;                  // number of faces
  egObject    **faces;                  // face objects
  int         topFlg;
  egadsBox    bbox;
};


class egadsMap
{
public:
  TopTools_IndexedMapOfShape map;
  egObject                   **objs;    // vector of egos w/ map
};


class egadsBody
{
public:
  TopoDS_Shape shape;                   // OCC topology
  egadsMap     nodes;
  egadsMap     edges;
  egadsMap     loops;
  egadsMap     faces;
  egadsMap     shells;
  int          *senses;                 // shell outer/inner (solids)
  egadsBox     bbox;
  int          massFill;
  double       massProp[14];
};


class egadsModel
{
public:
  TopoDS_Shape shape;                   // OCC topology
  int          nobjs;                   // number of total egObjects
  int          nbody;                   // number of bodies
  egObject     **bodies;                // vector of pointers to egObjects
  egadsBox     bbox;
};

#endif
