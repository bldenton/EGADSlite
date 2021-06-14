/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Internal UVmap Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "UVMAP_LIB.h"

//#define SMOOTHUV

FILE *UG_Output_File = NULL;

extern /*@null@*/ /*@out@*/ /*@only@*/ void *EG_alloc( size_t nbytes );
extern /*@null@*/ /*@only@*/ void *EG_reall( /*@null@*/ /*@only@*/ 
                                  /*@returned@*/ void *ptr, size_t nbytes );
extern void EG_free( /*@null@*/ /*@only@*/ void *pointer );

extern int  EG_inTriExact( double *t1, double *t2, double *t3, double *p,
                           double *w );



/*@null@*/ /*@out@*/ /*@only@*/ static void *
EG_uvmapMalloc(INT_ *err_flag, size_t size)
{
  void *ptr;

  *err_flag = EGADS_SUCCESS;
  ptr = EG_alloc(size);
  if (ptr == NULL) *err_flag = EGADS_MALLOC;
  return ptr;
}


/*@null@*/ /*@only@*/ static void *
EG_uvmapRealloc(INT_ *err_flag, void *ptr, size_t size)
{
  void *tmp;

  *err_flag = EGADS_SUCCESS;
  tmp = EG_reall(ptr, size);
  if (tmp == NULL) *err_flag = EGADS_MALLOC;
  return tmp;
}


static void
EG_triRemap(int *trmap, int itri, int flag, int *verts, double *ws)
{
  int    i1, i2, i3, tris[3];
  double w[3];
  
  if (trmap         == NULL) return;
  if (trmap[itri-1] == 0)    return;
  
  tris[0]  = verts[0];
  tris[1]  = verts[1];
  tris[2]  = verts[2];
  w[0]     = ws[0];
  w[1]     = ws[1];
  w[2]     = ws[2];
  
  i1       =  trmap[itri-1]       & 3;
  i2       = (trmap[itri-1] >> 2) & 3;
  i3       = (trmap[itri-1] >> 4) & 3;
  
  if (flag == 1) {
    /* EGADS tri is source -- uvmap is destination */
    verts[0] = tris[i1-1];
    verts[1] = tris[i2-1];
    verts[2] = tris[i3-1];
    ws[0]    = w[i1-1];
    ws[1]    = w[i2-1];
    ws[2]    = w[i3-1];
  } else {
    /* uvmap is source -- EGADS tri is destination */
    verts[i1-1] = tris[0];
    verts[i2-1] = tris[1];
    verts[i3-1] = tris[2];
    ws[i1-1]    = w[0];
    ws[i2-1]    = w[1];
    ws[i3-1]    = w[2];
  }
}


#ifdef SMOOTHUV
static int
EG_uvmapSmooth(uvmap_struct *uvstruct, double *xyzs)
{
  int    i, i1, i2, i3, iter, *verts;
  double d, delta, maxd, uv[2], *sums;
  
  /* find bounding vertices */
  verts = (int *) EG_alloc(uvstruct->nnode*sizeof(int));
  if (verts == NULL) return EGADS_MALLOC;
  for (i = 0; i < uvstruct->nnode; i++) verts[i] = 0;
  
  for (i = 1; i <= uvstruct->nbface; i++) {
    i1   = uvstruct->inibf[i][0];
    i2   = uvstruct->inibf[i][1];
    i3   = uvstruct->inibf[i][2];
    if (uvstruct->ibfibf[i][0] < 0) {
      verts[i2-1]++;
      verts[i3-1]++;
    }
    if (uvstruct->ibfibf[i][1] < 0) {
      verts[i1-1]++;
      verts[i3-1]++;
    }
    if (uvstruct->ibfibf[i][2] < 0) {
      verts[i1-1]++;
      verts[i2-1]++;
    }
  }

  /* smooth the interior */
  sums = (double *) EG_alloc(3*uvstruct->nnode*sizeof(double));
  if (sums == NULL) {
    EG_free(verts);
    return EGADS_MALLOC;
  }
  
  /* iterate */
  for (iter = 0; iter < 100000; iter++) {
    for (i = 0; i < 3*uvstruct->nnode; i++) sums[i] = 0.0;
    
    for (i  = 0; i < uvstruct->nbface; i++) {
      i1    = uvstruct->inibf[i+1][0];
      i2    = uvstruct->inibf[i+1][1];
      i3    = uvstruct->inibf[i+1][2];
      uv[0] = uvstruct->u[i1][0] + uvstruct->u[i2][0] + uvstruct->u[i3][0];
      uv[1] = uvstruct->u[i1][1] + uvstruct->u[i2][1] + uvstruct->u[i3][1];
      sums[3*i1-3] += uv[0] - 3.0*uvstruct->u[i1][0];
      sums[3*i1-2] += uv[1] - 3.0*uvstruct->u[i1][1];
      sums[3*i1-1] += 2.0;
      sums[3*i2-3] += uv[0] - 3.0*uvstruct->u[i2][0];
      sums[3*i2-2] += uv[1] - 3.0*uvstruct->u[i2][1];
      sums[3*i2-1] += 2.0;
      sums[3*i3-3] += uv[0] - 3.0*uvstruct->u[i3][0];
      sums[3*i3-2] += uv[1] - 3.0*uvstruct->u[i3][1];
      sums[3*i3-1] += 2.0;
    }
  
    i1   = 0;
    maxd = delta = 0.0;
    for (i = 0; i < uvstruct->nnode; i++) {
      if (verts[i]    != 0)   continue;
      if (sums[3*i+2] == 0.0) continue;
      uv[0] = 0.5*sums[3*i  ]/sums[3*i+2];
      uv[1] = 0.5*sums[3*i+1]/sums[3*i+2];
      d     = sqrt(uv[0]*uv[0] + uv[1]*uv[1]);
      if (d > maxd) maxd = d;
      delta += d;
      i1++;

      uvstruct->u[i+1][0] += uv[0];
      uvstruct->u[i+1][1] += uv[1];
    }
    if (delta/i1 < 1.e-13) break;

  }
  printf(" EGADS Info: iter %d  maxUVdelta = %le  RMS = %le\n",
         iter+1, maxd, delta/i1);
  
  /* check UV areas */
  for (i = 1; i <= uvstruct->nbface; i++) {
    i1   = uvstruct->inibf[i][0];
    i2   = uvstruct->inibf[i][1];
    i3   = uvstruct->inibf[i][2];
    d    = ((uvstruct->u[i1][0]-uvstruct->u[i3][0])*
            (uvstruct->u[i2][1]-uvstruct->u[i3][1]) -
            (uvstruct->u[i1][1]-uvstruct->u[i3][1])*
            (uvstruct->u[i2][0]-uvstruct->u[i3][0]));
    if (d <= 0.0)
      printf(" EGADS Error: Folded UV space @ tri %d = %lf\n", i, d);
  }

  EG_free(sums);
  EG_free(verts);
  return EGADS_SUCCESS;
}
#endif


/*               ********** Exposed Entry Points **********               */

void
EG_uvmapInit()
{
  uvmap_register_ext_free(EG_free);
  uvmap_register_ext_malloc(EG_uvmapMalloc);
  uvmap_register_ext_realloc(EG_uvmapRealloc);
}


/* generate the UVmap from a connected set of triangles
 *
 * ntri  = number of input triangles
 * tris  = input triangles (3*ntri in length)
 *         each 3 indices (1-bias) create a right-handed triangle
 * itris = the source Face index for the triangle (ntri in length)
 * nvrt  = number of vertices supporting the triangles
 * vrts  = the input 3D coordinates (3*nvrt in Length)
 *         each triad contains the XYZs for the vertex
 * trmap = the returned pointer to a triangle mapping (can be null)
 * uvmap = the returned pointer to the internal uvmap structure
 *
 */

int
EG_uvmapMake(int ntri, int *tris, int *itris, int nvrt, double *vrts,
             double *range, int **trmap, void **uvmap)
{
  int          i, i0, i1, i2, n, stat, *map;
  double       *uv = NULL;
  uvmap_struct *uvstruct;
  
  *uvmap = NULL;
  *trmap = NULL;
  stat   = EG_uvmapGen(1, ntri, nvrt, 1, 0, itris, tris, vrts, &uv, uvmap);
  if ((stat != EGADS_SUCCESS) || (uv == NULL)) {
    printf(" EGADS Error: EG_uvmap_gen = %d (EG_uvmapMake)!\n", stat);
    return stat;
  }
  
  range[0] = range[1] = uv[0];
  range[2] = range[3] = uv[1];
  for (i = 1; i < nvrt; i++) {
    if (uv[2*i  ] < range[0]) range[0] = uv[2*i  ];
    if (uv[2*i  ] > range[1]) range[1] = uv[2*i  ];
    if (uv[2*i+1] < range[2]) range[2] = uv[2*i+1];
    if (uv[2*i+1] > range[3]) range[3] = uv[2*i+1];
  }
  EG_free(uv);
  uvstruct = (uvmap_struct *) *uvmap;
  if (uvstruct == NULL) {
    printf(" EGADS Error: NULL uvmap structure (EG_uvmapMake)!\n");
    return EGADS_UVMAP;
  }
  /* for some reason this is not set! */
  uvstruct->nnode = nvrt;

#ifdef SMOOTHUV
  /* smooth our new UVs */
  stat = EG_uvmapSmooth(uvstruct, vrts);
  if (stat != EGADS_SUCCESS)
    printf(" EGADS Warning: uvmapSmooth = %d (EG_uvmapMake)!\n", stat);
#endif

  /* see if we have had an index re-ordering */
  n = 0;
  for (i = 1; i <= ntri; i++) {
    if ((uvstruct->inibf[i][0] == tris[3*i-3]) ||
        (uvstruct->inibf[i][1] == tris[3*i-2]) ||
        (uvstruct->inibf[i][2] == tris[3*i-1])) continue;
    n++;
  }
  if (n != 0) {
    printf(" EGADS Info: %d of %d triangles reordered (EG_uvmapMake)!\n",
           n, ntri);
    map = (int *) EG_alloc(ntri*sizeof(int));
    if (map == NULL) {
      printf(" EGADS Error: Malloc of %d triangle indices (EG_uvmapMake)!\n",
             ntri);
      uvmap_struct_free(uvmap);
      return EGADS_MALLOC;
    }
    for (i = 0; i <  ntri; i++) map[i] = 0;
    for (i = 1; i <= ntri; i++) {
      if ((uvstruct->inibf[i][0] == tris[3*i-3]) ||
          (uvstruct->inibf[i][1] == tris[3*i-2]) ||
          (uvstruct->inibf[i][2] == tris[3*i-1])) continue;
      i0 = i1 = i2 = 0;
      if  (uvstruct->inibf[i][0] == tris[3*i-3]) i0 = 1;
      if  (uvstruct->inibf[i][0] == tris[3*i-2]) i0 = 2;
      if  (uvstruct->inibf[i][0] == tris[3*i-1]) i0 = 3;
      if  (uvstruct->inibf[i][1] == tris[3*i-3]) i1 = 1;
      if  (uvstruct->inibf[i][1] == tris[3*i-2]) i1 = 2;
      if  (uvstruct->inibf[i][1] == tris[3*i-1]) i1 = 3;
      if  (uvstruct->inibf[i][2] == tris[3*i-3]) i2 = 1;
      if  (uvstruct->inibf[i][2] == tris[3*i-2]) i2 = 2;
      if  (uvstruct->inibf[i][2] == tris[3*i-1]) i2 = 3;
      if ((i0 == 0) || (i1 == 0) || (i2 == 0)) {
        printf(" EGADS Error: Tri mapping %d %d %d  %d %d %d (EG_uvmapMake)!\n",
               tris[3*i-3], tris[3*i-2], tris[3*i-1], uvstruct->inibf[i][0],
               uvstruct->inibf[i][1], uvstruct->inibf[i][2]);
        uvmap_struct_free(uvmap);
        return EGADS_UVMAP;
      }
      map[i-1] = (i2 << 4) + (i1 << 2) + i0;
    }
    *trmap = map;
  }

  return EGADS_SUCCESS;
}


/* return triangle containing the input UV
 *
 * uvmap = pointer to the internal uvmap structure
 * trmap = pointer to triangle map -- can be null
 * uv    = the input target UV (2 in length)
 * fID   = the returned Face ID
 * itri  = the returned index (1-bias) into tris for found location
 * verts = the 3 vertex indices for the triangle
 * ws    = the weights in the triangle for the vertices (3 in len)
 */

int
EG_uvmapLocate(void *uvmap, int *trmap, double *uv, int *fID, int *itri,
               int *verts, double *ws)
{
  int          i, i1, i2, i3, stat, cls = 0;
  double       w[3], neg = 0.0;
  uvmap_struct *uvstruct;
  
  verts[0] = verts[1] = verts[2] = 0;
  ws[0]    = ws[1]    = ws[2]    = 0.0;
  stat = EG_uvmapFindUV(1, uv, uvmap, fID, itri, verts, ws);
  EG_triRemap(trmap, *itri, 0, verts, ws);
  if (stat != EGADS_NOTFOUND) return stat;
  
  /* need to extrapolate */
  uvstruct = (uvmap_struct *) uvmap;
  for (i = 1; i <= uvstruct->nbface; i++) {
    i1   = uvstruct->inibf[i][0];
    i2   = uvstruct->inibf[i][1];
    i3   = uvstruct->inibf[i][2];
    stat = EG_inTriExact(uvstruct->u[i1], uvstruct->u[i2], uvstruct->u[i3],
                         uv, w);
    if (stat == EGADS_SUCCESS) {
      *fID     = uvstruct->idibf[i];
      *itri    = i;
      ws[0]    = w[0];
      ws[1]    = w[1];
      ws[2]    = w[2];
      verts[0] = i1;
      verts[1] = i2;
      verts[2] = i3;
      EG_triRemap(trmap, i, 0, verts, ws);
/*    printf(" EGADS Info: Exact -> Ws = %le %le %le (EG_uvmapLocate)!\n",
             ws[0], ws[1], ws[2]);  */
      return EGADS_SUCCESS;
    }
    if (w[1] < w[0]) w[0] = w[1];
    if (w[2] < w[0]) w[0] = w[2];
    if (cls == 0) {
      cls = i;
      neg = w[0];
    } else {
      if (w[0] > neg) {
        cls = i;
        neg = w[0];
      }
    }
  }
  if (cls == 0) return EGADS_NOTFOUND;
  
  /* extrapolate */
  *fID  = uvstruct->idibf[cls];
  *itri = cls;
  i1    = uvstruct->inibf[cls][0];
  i2    = uvstruct->inibf[cls][1];
  i3    = uvstruct->inibf[cls][2];
  EG_inTriExact(uvstruct->u[i1], uvstruct->u[i2], uvstruct->u[i3], uv, ws);
  verts[0] = i1;
  verts[1] = i2;
  verts[2] = i3;
  EG_triRemap(trmap, cls, 0, verts, ws);
/*
  if (neg < -1.e-4) {
    printf(" EGADS Info: Extrapolate -> Ws = %le %le %le (EG_uvmapLocate)!\n",
           ws[0], ws[1], ws[2]);
    return EGADS_EXTRAPOL;
  }
 */
  return EGADS_SUCCESS;
}


/* return uvmap UV given index
 *
 * uvmap = pointer to the internal uvmap structure
 * index = index into uvmap
 * uv    = uv in uvmap
 */

void
EG_getUVmap(void *uvmap, int index, double *uv)
{
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  uv[0]    = uvstruct->u[index][0];
  uv[1]    = uvstruct->u[index][1];
}


/* return uvmap UV given uv and fID
 *
 * uvmap = pointer to the internal uvmap structure
 * trmap = pointer to triangle map -- can be null
 * fuv   = the input Face UV (2 in length)
 * fuvs  = uvs for the Face
 * tris  = tris indices for fuvs
 * tbeg  = triangle start index
 * tend  = triangle end index
 * uv    = uv in uvmap
 */

int
EG_uv2UVmap(void *uvmap, int *trmap, double *fuv, double *fuvs, int *tris,
            int tbeg, int tend, double *uv)
{
  int          i, j, stat, i1, i2, i3, verts[3], cls = 0;
  double       w[3], neg;
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  for (i = tbeg; i <= tend; i++) {
    j    = i - tbeg;
    i1   = tris[3*j  ] - 1;
    i2   = tris[3*j+1] - 1;
    i3   = tris[3*j+2] - 1;
    stat = EG_inTriExact(&fuvs[2*i1], &fuvs[2*i2], &fuvs[2*i3], fuv, w);
    if (stat == EGADS_SUCCESS) {
      verts[0] = uvstruct->inibf[i][0];
      verts[1] = uvstruct->inibf[i][1];
      verts[2] = uvstruct->inibf[i][2];
      EG_triRemap(trmap, i, 1, verts, w);
      i1    = verts[0];
      i2    = verts[1];
      i3    = verts[2];
      uv[0] = w[0]*uvstruct->u[i1][0] + w[1]*uvstruct->u[i2][0] +
              w[2]*uvstruct->u[i3][0];
      uv[1] = w[0]*uvstruct->u[i1][1] + w[1]*uvstruct->u[i2][1] +
              w[2]*uvstruct->u[i3][1];
      return EGADS_SUCCESS;
    }
    if (w[1] < w[0]) w[0] = w[1];
    if (w[2] < w[0]) w[0] = w[2];
    if (cls == 0) {
      cls = j;
      neg = w[0];
    } else {
      if (w[0] > neg) {
        cls = j;
        neg = w[0];
      }
    }
  }
  if (cls == 0) return EGADS_NOTFOUND;

  /* extrapolate */
  i1 = tris[3*cls  ] - 1;
  i2 = tris[3*cls+1] - 1;
  i3 = tris[3*cls+2] - 1;
  EG_inTriExact(&fuvs[2*i1], &fuvs[2*i2], &fuvs[2*i3], fuv, w);
  verts[0] = uvstruct->inibf[cls+tbeg][0];
  verts[1] = uvstruct->inibf[cls+tbeg][1];
  verts[2] = uvstruct->inibf[cls+tbeg][2];
  EG_triRemap(trmap, cls+tbeg, 1, verts, w);
  i1    = verts[0];
  i2    = verts[1];
  i3    = verts[2];
  uv[0] = w[0]*uvstruct->u[i1][0] + w[1]*uvstruct->u[i2][0] +
          w[2]*uvstruct->u[i3][0];
  uv[1] = w[0]*uvstruct->u[i1][1] + w[1]*uvstruct->u[i2][1] +
          w[2]*uvstruct->u[i3][1];
  
  return EGADS_SUCCESS;
}


int
EG_uvmapWrite(void *uvmap, int *trmap, FILE *fp)
{
  int          i, trmp = 0, msrch = 0;
  uvmap_struct *uvstruct;
  
  uvstruct = (uvmap_struct *) uvmap;
  if (uvstruct->mdef  != 1) return EGADS_UVMAP;
  if (uvstruct->ndef  != 1) return EGADS_UVMAP;
  if (uvstruct->msrch != NULL) msrch = 1;
  if (trmap           != NULL) trmp  = 1;
  fprintf(fp, "%d %d %d %d %d %d\n", (int) uvstruct->isrch,
          (int) uvstruct->ibface, (int) uvstruct->nbface, (int) uvstruct->nnode,
          msrch, trmp);
  
  for (i = 1; i <= uvstruct->nbface; i++) {
    fprintf(fp, "%d ", uvstruct->idibf[i]);
    if ((i%20 == 0) && (i != uvstruct->nbface)) fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  
  for (i = 1; i <= uvstruct->nbface; i++)
    fprintf(fp, "%d %d %d\n", uvstruct->inibf[i][0], uvstruct->inibf[i][1],
            uvstruct->inibf[i][2]);
  
  for (i = 1; i <= uvstruct->nbface; i++)
    fprintf(fp, "%d %d %d\n", uvstruct->ibfibf[i][0], uvstruct->ibfibf[i][1],
            uvstruct->ibfibf[i][2]);
        
  for (i = 1; i <= uvstruct->nnode; i++)
    fprintf(fp, "%19.12le %19.12le\n", uvstruct->u[i][0], uvstruct->u[i][1]);
  
  if (uvstruct->msrch != NULL) {
    for (i = 1; i <= uvstruct->nbface; i++) {
      fprintf(fp, "%d ", uvstruct->msrch[i]);
      if ((i%20 == 0) && (i != uvstruct->nbface)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
  if (trmap != NULL) {
    for (i = 1; i <= uvstruct->nbface; i++) {
      fprintf(fp, "%d ", trmap[i-1]);
      if ((i%20 == 0) && (i != uvstruct->nbface)) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
  }
    
  return EGADS_SUCCESS;
}


int
EG_uvmapRead(FILE *fp, double *range, void **uvmap, int **trmap)
{
  int          i, cnt, msrch, trmp, *map = NULL, stat = EGADS_SUCCESS;
  INT_         err = 0;
  uvmap_struct *uvstruct;
  
  *uvmap   = NULL;
  *trmap   = NULL;
  uvstruct = (uvmap_struct *) uvmap_malloc(&err, sizeof(uvmap_struct));
  if ((err != 0) || (uvstruct == NULL)) {
    printf(" EGADS Error: Failed to allocate UVmap (EG_uvmapRead)!\n");
    return EGADS_MALLOC;
  }
  uvstruct->ndef   = 1;
  uvstruct->mdef   = 1;
  uvstruct->idef   = 1;
  uvstruct->isrch  = 0;
  uvstruct->ibface = 0;
  uvstruct->nbface = 0;
  uvstruct->nnode  = 0;

  uvstruct->idibf  = NULL;
  uvstruct->msrch  = NULL;
  uvstruct->inibf  = NULL;
  uvstruct->ibfibf = NULL;
  uvstruct->u      = NULL;

  cnt = fscanf(fp, "%d %d %d %d %d %d", &uvstruct->isrch,  &uvstruct->ibface,
               &uvstruct->nbface, &uvstruct->nnode, &msrch, &trmp);
  if (cnt != 6) {
    printf(" EGADS Error: Header read len = %d (EG_uvmapRead)!\n", cnt);
    stat = EGADS_READERR;
    goto malloc;
  }
  
  uvstruct->idibf = (INT_ *) uvmap_malloc(&err,
                                          (uvstruct->nbface+1)*sizeof(INT_));
  if ((err != 0) || (uvstruct->idibf == NULL)) {
    printf(" EGADS Error: malloc %d id (EG_uvmapRead)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto malloc;
  }
  uvstruct->idibf[0] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    cnt = fscanf(fp, "%d", &uvstruct->idibf[i]);
    if (cnt != 1) {
      printf(" EGADS Error: idibf read len = %d/%d (EG_uvmapRead)!\n",
             i, uvstruct->nbface);
      stat = EGADS_READERR;
      goto malloc;
    }
  }
  
  uvstruct->inibf = (INT_3D *) uvmap_malloc(&err,
                                          (uvstruct->nbface+1)*sizeof(INT_3D));
  if ((err != 0) || (uvstruct->inibf == NULL)) {
    printf(" EGADS Error: malloc %d in (EG_uvmapRead)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto malloc;
  }
  uvstruct->inibf[0][0] = uvstruct->inibf[0][1] = uvstruct->inibf[0][2] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    cnt = fscanf(fp, "%d %d %d", &uvstruct->inibf[i][0], &uvstruct->inibf[i][1],
                 &uvstruct->inibf[i][2]);
    if (cnt != 3) {
      printf(" EGADS Error: inibf read len = %d/%d (EG_uvmapRead)!\n",
             i, uvstruct->nbface);
      stat = EGADS_READERR;
      goto malloc;
    }
  }
  
  uvstruct->ibfibf = (INT_3D *) uvmap_malloc(&err,
                                           (uvstruct->nbface+1)*sizeof(INT_3D));
  if ((err != 0) || (uvstruct->ibfibf == NULL)) {
    printf(" EGADS Error: malloc %d ibf (EG_uvmapRead)!\n", uvstruct->nbface);
    stat = EGADS_MALLOC;
    goto malloc;
  }
  uvstruct->ibfibf[0][0] = uvstruct->ibfibf[0][1] = uvstruct->ibfibf[0][2] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    cnt = fscanf(fp, "%d %d %d", &uvstruct->ibfibf[i][0],
                 &uvstruct->ibfibf[i][1], &uvstruct->ibfibf[i][2]);
    if (cnt != 3) {
      printf(" EGADS Error: ibfibf read len = %d/%d (EG_uvmapRead)!\n",
             i, uvstruct->nbface);
      stat = EGADS_READERR;
      goto malloc;
    }
  }

  uvstruct->u = (DOUBLE_2D *) uvmap_malloc(&err,
                                         (uvstruct->nnode+1)*sizeof(DOUBLE_2D));
  if ((err != 0) || (uvstruct->u == NULL)) {
    printf(" EGADS Error: malloc %d u (EG_uvmapRead)!\n", uvstruct->nnode);
    stat = EGADS_MALLOC;
    goto malloc;
  }
  uvstruct->u[0][0] = uvstruct->u[0][1] = 0.0;
  for (i = 1; i <= uvstruct->nnode; i++) {
    cnt = fscanf(fp, "%lf %lf", &uvstruct->u[i][0], &uvstruct->u[i][1]);
    if (cnt != 2) {
      printf(" EGADS Error: u read len = %d/%d (EG_uvmapRead)!\n",
             i, uvstruct->nnode);
      stat = EGADS_READERR;
      goto malloc;
    }
  }

  range[0] = range[1] = uvstruct->u[1][0];
  range[2] = range[3] = uvstruct->u[1][1];
  for (i = 2; i <= uvstruct->nnode; i++) {
    if (uvstruct->u[i][0] < range[0]) range[0] = uvstruct->u[i][0];
    if (uvstruct->u[i][0] > range[1]) range[1] = uvstruct->u[i][0];
    if (uvstruct->u[i][1] < range[2]) range[2] = uvstruct->u[i][1];
    if (uvstruct->u[i][1] > range[3]) range[3] = uvstruct->u[i][1];
  }
  
  if (msrch == 1) {
    uvstruct->msrch = (INT_ *) uvmap_malloc(&err,
                                            (uvstruct->nbface+1)*sizeof(INT_));
    if ((err != 0) || (uvstruct->msrch == NULL)) {
      printf(" EGADS Error: malloc %d msrch (EG_uvmapRead)!\n",
             uvstruct->nbface);
      stat = EGADS_MALLOC;
      goto malloc;
    }
    uvstruct->msrch[0] = 0;
    for (i = 1; i <= uvstruct->nbface; i++) {
      cnt = fscanf(fp, "%d", &uvstruct->msrch[i]);
      if (cnt != 1) {
        printf(" EGADS Error: msrch read len = %d/%d (EG_uvmapRead)!\n",
               i, uvstruct->nbface);
        stat = EGADS_READERR;
        goto malloc;
      }
    }
  }

  if (trmp == 1) {
    map = (int *) EG_alloc(uvstruct->nbface*sizeof(int));
    if (map == NULL) {
      printf(" EGADS Error: malloc %d trmap (EG_uvmapRead)!\n",
             uvstruct->nbface);
      stat = EGADS_MALLOC;
      goto malloc;
    }
    for (i = 0; i < uvstruct->nbface; i++) {
      cnt = fscanf(fp, "%d", &map[i]);
      if (cnt != 1) {
        printf(" EGADS Error: trmap read len = %d/%d (EG_uvmapRead)!\n",
               i, uvstruct->nbface);
        stat = EGADS_READERR;
        goto malloc;
      }
    }
    *trmap = map;
  }

  *uvmap  = uvstruct;
  return EGADS_SUCCESS;
  
malloc:
/*@+dependenttrans@*/
  uvmap_free(uvstruct->idibf);
/*@-nullpass@*/
  uvmap_free(uvstruct->msrch);
/*@+nullpass@*/
  uvmap_free(uvstruct->inibf);
  uvmap_free(uvstruct->ibfibf);
  uvmap_free(uvstruct->u);
  uvmap_free(uvstruct);
  return stat;
}


int
EG_uvmapCopy(void *uvsrc, int *trsrc, void **uvmap, int **trmap)
{
  int          i, *map = NULL;
  INT_         err = 0;
  uvmap_struct *uvstruct, *sstruct;
  
  *uvmap   = NULL;
  *trmap   = NULL;
  sstruct  = (uvmap_struct *) uvsrc;
  uvstruct = (uvmap_struct *) uvmap_malloc(&err, sizeof(uvmap_struct));
  if ((err != 0) || (uvstruct == NULL)) {
    printf(" EGADS Error: Failed to allocate UVmap (EG_uvmapRead)!\n");
    return EGADS_MALLOC;
  }
  uvstruct->ndef   = 1;
  uvstruct->mdef   = 1;
  uvstruct->idef   = 1;
  uvstruct->isrch  = sstruct->isrch;
  uvstruct->ibface = sstruct->ibface;
  uvstruct->nbface = sstruct->nbface;
  uvstruct->nnode  = sstruct->nnode;

  uvstruct->idibf  = NULL;
  uvstruct->msrch  = NULL;
  uvstruct->inibf  = NULL;
  uvstruct->ibfibf = NULL;
  uvstruct->u      = NULL;
  
  uvstruct->idibf = (INT_ *) uvmap_malloc(&err,
                                          (uvstruct->nbface+1)*sizeof(INT_));
  if ((err != 0) || (uvstruct->idibf == NULL)) {
    printf(" EGADS Error: malloc %d id (EG_uvmapCopy)!\n", uvstruct->nbface);
    goto malloc;
  }
  uvstruct->idibf[0] = 0;
  for (i = 1; i <= uvstruct->nbface; i++)
    uvstruct->idibf[i] = sstruct->idibf[i];
  
  uvstruct->inibf = (INT_3D *) uvmap_malloc(&err,
                                           (uvstruct->nbface+1)*sizeof(INT_3D));
  if ((err != 0) || (uvstruct->inibf == NULL)) {
    printf(" EGADS Error: malloc %d in (EG_uvmapCopy)!\n", uvstruct->nbface);
    goto malloc;
  }
  uvstruct->inibf[0][0] = uvstruct->inibf[0][1] = uvstruct->inibf[0][2] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    uvstruct->inibf[i][0] = sstruct->inibf[i][0];
    uvstruct->inibf[i][1] = sstruct->inibf[i][1];
    uvstruct->inibf[i][2] = sstruct->inibf[i][2];
  }
  
  uvstruct->ibfibf = (INT_3D *) uvmap_malloc(&err,
                                           (uvstruct->nbface+1)*sizeof(INT_3D));
  if ((err != 0) || (uvstruct->ibfibf == NULL)) {
    printf(" EGADS Error: malloc %d ibf (EG_uvmapCopy)!\n", uvstruct->nbface);
    goto malloc;
  }
  uvstruct->ibfibf[0][0] = uvstruct->ibfibf[0][1] = uvstruct->ibfibf[0][2] = 0;
  for (i = 1; i <= uvstruct->nbface; i++) {
    uvstruct->ibfibf[i][0] = sstruct->ibfibf[i][0];
    uvstruct->ibfibf[i][1] = sstruct->ibfibf[i][1];
    uvstruct->ibfibf[i][2] = sstruct->ibfibf[i][2];
  }

  uvstruct->u = (DOUBLE_2D *) uvmap_malloc(&err,
                                         (uvstruct->nnode+1)*sizeof(DOUBLE_2D));
  if ((err != 0) || (uvstruct->u == NULL)) {
    printf(" EGADS Error: malloc %d u (EG_uvmapCopy)!\n", uvstruct->nnode);
    goto malloc;
  }
  uvstruct->u[0][0] = uvstruct->u[0][1] = 0.0;
  for (i = 1; i <= uvstruct->nnode; i++) {
    uvstruct->u[i][0] = sstruct->u[i][0];
    uvstruct->u[i][1] = sstruct->u[i][1];
  }
  
  if (sstruct->msrch != NULL) {
    uvstruct->msrch = (INT_ *) uvmap_malloc(&err,
                                            (uvstruct->nbface+1)*sizeof(INT_));
    if ((err != 0) || (uvstruct->msrch == NULL)) {
      printf(" EGADS Error: malloc %d msrch (EG_uvmapCopy)!\n",
             uvstruct->nbface);
      goto malloc;
    }
    uvstruct->msrch[0] = 0;
    for (i = 1; i <= uvstruct->nbface; i++)
      uvstruct->msrch[i] = sstruct->msrch[i];
  }

  if (trsrc != NULL) {
    map = (int *) EG_alloc(uvstruct->nbface*sizeof(int));
    if (map == NULL) {
      printf(" EGADS Error: malloc %d trmap (EG_uvmapCopy)!\n",
             uvstruct->nbface);
      goto malloc;
    }
    for (i = 0; i < uvstruct->nbface; i++) map[i] = trsrc[i];
    *trmap = map;
  }

  *uvmap  = uvstruct;
  return EGADS_SUCCESS;
  
malloc:
/*@+dependenttrans@*/
  uvmap_free(uvstruct->idibf);
/*@-nullpass@*/
  uvmap_free(uvstruct->msrch);
/*@+nullpass@*/
  uvmap_free(uvstruct->inibf);
  uvmap_free(uvstruct->ibfibf);
  uvmap_free(uvstruct->u);
  uvmap_free(uvstruct);
  return EGADS_MALLOC;
}
