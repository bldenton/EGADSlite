#include <math.h>
#include <string.h>

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* shifted erf = 0 at x = -1 and 1 at x = 1
           min is better than max (eg angles) a < b
           max is better than min (sizes)     b < a    */
#define ERF(a,b,x)        (0.5*(1.-erf(3.554147*(x-(0.5*(a+b)))/(0.5*(a-b)))))

#define CROSS(a,b,c)      c[0] = (a[1]*b[2]) - (a[2]*b[1]);\
                          c[1] = (a[2]*b[0]) - (a[0]*b[2]);\
                          c[2] = (a[0]*b[1]) - (a[1]*b[0])
#define DOT(a,b)          (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define DOT4(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define PI                  3.1415926535897931159979635
#define PIEPS               3.2
#define qEPS                1.e-14
#define EPS08               1.E-08
#define ANGCUT              2.9  // ~ 170 DEG
#define ISQUAD              0
#define ISVERTEX            1
#define SWAP                0
#define COLLAPSE            1
#define SPLIT               2
#define SWAPCOLLAPSE        0
#define DOUBLECOLLAPSE      1
#define SWAPDOUBLECOLLAPSE  2
#define SWAPSPLIT           3
#define DOUBLESPLIT         4
#define SWAPDOUBLESPLIT     5
#define DOUBLESWAP          6
#define QA0                 0
#define QA1                 1
#define QA2                 200
#define QA3                 40000
#define QACB                8000000  /* crosses domain bounds: SUPER INVALID */


typedef struct{
  int verts[4], qadj[4], id;
} Quad;


typedef struct{
  int    *verts, *quads, type; /* -1 interior
                                   0 its vertices are linked to bounds
                                   1 links to bounds directly  */
  int    nV, nQ;               /* nV = n + 1 =  origin(1) + peaks (n)  */
  int    *idxV, *idxQ, *area;
  double *angle, *ratio;
} vStar;


typedef struct {
  int    fID, oriQ, oriV, plotcount, totQ, totV, pp, sizeV, sizeQ, *qIdx,
         *qAdj, **valence, *vInv, *vType, *remQ, *remV, invsteps, regBd, *degen;
  ego    face;
  double range[4],  *xyzs, *uvs, minArea, maxArea, avArea, *bdAng, fin;
  vStar  **star;
} meshMap;


typedef struct {
  ego     tess,   *faces;
  int     nedges, nfaces;
  meshMap **qm;
} bodyQuad;


typedef struct{
  int verts[6], vals[6];
  int q[2];
} quadGroup;

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif


__HOST_AND_DEVICE__ int  EG_outLevel(const egObject *object);

__HOST_AND_DEVICE__ int  EG_createMeshMap(bodyQuad *bodydata);
__HOST_AND_DEVICE__ int  EG_meshRegularization(meshMap *qm);
__HOST_AND_DEVICE__ int  EG_makeQuadTess(bodyQuad bodydata, ego *quadTess);
__HOST_AND_DEVICE__ void EG_destroyMeshMap(bodyQuad *bodydata);
