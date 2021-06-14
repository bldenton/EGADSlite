#ifndef EGADSTYPES_H
#define EGADSTYPES_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             General Object Header
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egadsErrors.h"


#define EGADSMAJOR     1
#define EGADSMINOR    19
#define EGADSPROP     EGADSprop: Revision 1.19

#define MAGIC      98789
#define MTESSPARAM     2

/* OBJECT CLASSES */

#define CONTXT         0
#define TRANSFORM      1
#define TESSELLATION   2
#define NIL            3        /* allocated but not assigned */
#define EMPTY          4        /* in the pool */
#define REFERENCE      5
#define PCURVE        10
#define CURVE         11
#define SURFACE       12
#define NODE          20
#define EDGE          21
#define LOOP          22
#define FACE          23
#define SHELL         24
#define BODY          25
#define MODEL         26
#define EEDGE         31
#define ELOOPX        32        /* ELOOP conflicts with errno.h */
#define EFACE         33
#define ESHELL        34
#define EBODY         35


/* MEMBER TYPES */

  /* PCURVES & CURVES */
#define LINE           1
#define CIRCLE         2
#define ELLIPSE        3
#define PARABOLA       4
#define HYPERBOLA      5
#define TRIMMED        6
#define BEZIER         7
#define BSPLINE        8
#define OFFSET         9

  /* SURFACES */
#define PLANE          1
#define SPHERICAL      2
#define CYLINDRICAL    3
#define REVOLUTION     4 
#define TOROIDAL       5
#define CONICAL       10
#define EXTRUSION     11

  /* TOPOLOGY */
#define SREVERSE      -1
#define NOMTYPE        0
#define SFORWARD       1
#define ONENODE        1
#define TWONODE        2
#define OPEN           3
#define CLOSED         4
#define DEGENERATE     5
#define WIREBODY       6
#define FACEBODY       7
#define SHEETBODY      8
#define SOLIDBODY      9


/* ATTRIBUTE TYPES */

#define ATTRINT        1
#define ATTRREAL       2
#define ATTRSTRING     3
#define ATTRCSYS      12
#define ATTRPTR       13


/* SOLID/GENERAL BOOLEAN OPERATIONS */

#define SUBTRACTION    1
#define INTERSECTION   2
#define FUSION         3
#define SPLITTER       4


/* SOLID BODY TYPES */

#define BOX            1
#define SPHERE         2
#define CONE           3
#define CYLINDER       4
#define TORUS          5


/* ISOCLINE TYPES */

#define UISO	       0
#define VISO           1


/* FACE SOURCE TYPES */

#define NODEOFF        1
#define EDGEOFF        2
#define FACEDUP        3
#define FACECUT        4
#define FACEOFF        5


typedef struct {
  char     *name;               /* Attribute Name */
  int      type;                /* Attribute Type */
  int      length;              /* number of values */
  union {
    int    integer;             /* single int -- length == 1 */
    int    *integers;           /* multiple ints */
    double real;                /* single double -- length == 1 */
    double *reals;              /* mutiple doubles */
    char   *string;             /* character string (no single char) */
  } vals;
} egAttr;


typedef struct {
  char *root;                   /* root name for the sequenced attribute */
  int  nSeq;                    /* number in the sequence (2 and on) */
  int  *attrSeq;                /* vector of ordered sequenced attributes */
} egAttrSeq;


typedef struct {
  int       nattrs;             /* number of attributes */
  egAttr    *attrs;             /* the attributes */
  int       nseqs;              /* number of sequenced attributes */
  egAttrSeq *seqs;              /* the sequenced attributes */
} egAttrs;


typedef struct egObject {
  int     magicnumber;          /* must be set to validate the object */
  short   oclass;               /* object Class */
  short   mtype;                /* member Type */
  void   *attrs;                /* object Attributes or Reference */
  void   *blind;		/* blind pointer to object data */
  struct egObject *topObj;      /* top of the hierarchy or context (if top) */
  struct egObject *tref;        /* threaded list of references */
  struct egObject *prev;        /* back pointer */
  struct egObject *next;        /* forward pointer */
} egObject;
typedef struct egObject* ego;


typedef struct {
  int      outLevel;		/* output level for messages
                                   0 none, 1 minimal, 2 verbose, 3 debug */
  int      fixedKnots;          /* always use evenly spaced knots */
  int      fullAttrs;           /* persist Edge and Node Attributes */
  double   tess[MTESSPARAM];    /* global tessellation parameters */
  char     **signature;
  void     *usrPtr;
  long     threadID;            /* the OS' thread identifier */
  void     *mutex;              /* this thread's mutex */
  egObject *pool;               /* available object structures for use */
  egObject *last;               /* the last object in the list */
} egCntxt;


typedef struct {
  int index;                    /* face index or last for more than 1 */
  int nface;                    /* number of faces (when more than 1) */
  int *faces;                   /* the face indices when multiples */
  int *tric;                    /* connection into tris for each face */
} egFconn;


typedef struct {
  int    tri;                   /* triangle index */
  double w[2];                  /* barycentric coordinates (weights for 0/1) */
} egBary;


typedef struct {
  egObject *obj;                /* edge object */
  int      nodes[2];            /* node indices */
  egFconn  faces[2];            /* minus and plus face connectivity */
  double   *xyz;                /* point coordinates */
  double   *t;                  /* parameter values */
  int      *global;             /* global indices */
  int      npts;                /* number of points */
} egTess1D;


typedef struct {
  int *ipts;                    /* index for point (nu*nv) */
  int *bounds;                  /* bound index (2*nu + 2*nv) */
  int nu;
  int nv;
} egPatch;


typedef struct {
  egObject *mKnots;             /* mapped surface if exists */
  double   *xyz;                /* point coordinates */
  double   *uv;                 /* parameter values */
  int      *global;             /* global indices */
  int      *ptype;              /* point type */
  int      *pindex;             /* point index */
  egBary   *bary;               /* barycentric coordinates in the frame */
  int      *frame;              /* initial triangulation */
  int      *frlps;              /* initial triangulation loop counts */
  int      *tris;               /* triangle indices */
  int      *tric;               /* triangle neighbors */
  egPatch  *patch;              /* patches */
  int      npts;                /* number of points */
  int      nframe;              /* number of tris in the initial frame */
  int      nfrlps;              /* number of loops in the initial frame */
  int      ntris;               /* number of tris */
  int      npatch;              /* number of quad patches */
  int      tfi;                 /* tfi flag (0 - not, 1 - tfi used) */
} egTess2D;


typedef struct {
  egObject *src;                /* source of the tessellation */
  double   *xyzs;               /* storage for geom */
  egTess1D *tess1d;             /* Edge tessellations */
  egTess2D *tess2d;             /* Face tessellations (tris then quads) */
  int      *globals;            /* global definitions */
  double   params[6];           /* suite of parameters used */
  double   tparam[MTESSPARAM];
  int      nGlobal;             /* number of Global vertices */
  int      nEdge;               /* number of Edge tessellations */
  int      nFace;               /* number of Face tessellations */
  int      nu;                  /* number of us for surface / ts for curve */
  int      nv;                  /* number of vs for surface tessellation */
  int      done;                /* is Object complete? */
} egTessel;


typedef struct {
  int       nobjs;              /* number in the map */
  egObject **objs;              /* vector of egos for map */
} egEMap;


typedef struct {
  egObject *edge;               /* Bounding Edge */
  int      curve;               /* -1 internal, 0 - line, 1 - curve */
  int      npts;                /* number of verts */
  double   *ts;                 /* the t values (npts in length) */
  double   dstart[3];           /* displacement in xyz first Node */
  double   dend[3];             /* displacement in xyz last Node */
} egEdVert;


typedef struct {
  int      iedge;               /* the Edge Index */
  int      sense;               /* sense use for Edge */
  egObject *nstart;             /* Node object @ beginning */
  double   tstart;              /* t for beginning of segment */
  double   tend;                /* t for end of segment */
} egEEseg;


typedef struct {
  egEdVert *sedges;             /* source Edge structure */
  int      nsegs;               /* number of Edge segments */
  egEEseg  *segs;               /* Edge segments */
  double   trange[2];
  egObject *nodes[2];           /* pointer to ego Nodes */
} egEEdge;


typedef struct {
  egObject *edge;               /* Edge object */
  int      sense;               /* sense use for Edge */
  int      npts;                /* number of discrete points */
  int      *iuv;                /* index into uvmap for Edge ts */
} egEdgeUV;


typedef struct {
  egEMap   eedges;              /* list of EEdges (in order) */
  int      *senses;             /* sense for each EEdge */
  double   area;                /* area in UVmap */
  int      nedge;               /* length of edgeUVs */
  egEdgeUV *edgeUVs;            /* the UVs for each Edge in the ELoop */
} egELoop;


typedef struct {
  egObject *face;               /* Face object */
  double   tol;                 /* max displacement magnitude */
  int      start;               /* offset for start in larger triangulation */
  int      nuvs;                /* length of uvs (x 2) */
  int      ndeflect;            /* length of deflect (x 3) */
  int      ntris;               /* number of triangles (single Face) */
  int      *uvtris;             /* tri indices into uvs */
  int      *uvtric;             /* neighbors (NULL until needed) */
  double   *uvs;                /* UVs for the triangle vertices -- Face Tess */
  double   *deflect;            /* displacement in xyz for frame indices */
} egEPatch;


typedef struct {
  egEdVert *sedges;             /* source Edge structure */
  int      npatch;              /* number of Faces */
  egEPatch *patches;            /* the Face(s) and discrete data */
  egEMap   eloops;              /* list of ELoops */
  int      *senses;
  int      *trmap;              /* triangle remapping -- can be null */
  void     *uvmap;              /* UVmap structure */
  double   range[4];
  int      last;                /* last triangle -- single Face */
} egEFace;


typedef struct {
  egEMap efaces;                /* list of EFaces */
} egEShell;


typedef struct {
  egObject *ref;                /* source of the EBody (Body Tess or Body) */
  egEMap   eedges;
  egEMap   eloops;
  egEMap   efaces;
  egEMap   eshells;
  int      *senses;
  double   angle;               /* Open Edge Node removal criteria */
  int      done;                /* is EBody complete? */
  int      nedge;               /* the number of source Edges */
  egEdVert *edges;              /* the source Edge discretizations */
} egEBody;

#endif
