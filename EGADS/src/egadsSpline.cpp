/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Blend, Rule and Sculpting Functions
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "egads.h"
#include "egads_dot.h"
#include "egadsClasses.h"
#include "egadsStack.h"
#include "Surreal/SurrealS.h"

#define EGADS_SPLINE_VELS

#ifdef EGADS_SPLINE_VELS
#include "egadsSplineVels.h"
#endif

//#define PCURVE_SENSITIVITY

//#define SPLINE_TECPLOT_DEBUG
//#define BLEND_FPE
//#define BLEND_SPLIT_CONSTRUCTION
//#define MAXIMUM_RATIO_WARNING
//#define MATCH_TIP_SNOR_NNOR (talk to Elaine about this one...)
//#define DUMP_SECTIONS

/* Ordering of edges in loops for spline faces, must
 * walk counter-clockwise around the uv-rectangle.
 * Used in EG_splineGeom and EG_splineGeom_dot.
 */
#define UMIN 0
#define VMIN 1
#define UMAX 2
#define VMAX 3

#if !defined(WIN32) && defined(BLEND_FPE)
// floating point exceptions
  #include <fenv.h>

#if defined(__APPLE__) && defined(__MACH__)

// This fixes the missing feenbableexcept on MacOS.  The fix reuses in verbatim a code available on
// https://github.com/ArduPilot/ardupilot/blob/master/libraries/AP_Common/missing/fenv.h.

// Public domain polyfill for feenableexcept on OS X
// http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c

int feenableexcept(unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
  // previous masks
  unsigned int old_excepts;

  if (fegetenv(&fenv))
    return -1;

  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return fesetenv(&fenv) ? -1 : old_excepts;
}
#endif
#endif

#define PI               3.1415926535897931159979635

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define MAX(a,b)         ((a) > (b) ? (a) : (b))
#define MIN(a,b)         ((a) < (b) ? (a) : (b))

/* number of knots to define half of a tip treatment */
#define NTIP 11

/* number of samples to construct tip TE surface */
#define NSUBTE 1
#define NTIPTE (NTIP+1 + NSUBTE*NTIP)

#define INFTY 1e200 /* consistent with OCC */

extern "C" /*@kept@*/ /*@null@*/ egObject *EG_context( const ego object );
extern "C" int EG_sameThread( const ego object );
extern "C" int EG_outLevel( const ego object );
extern "C" int EG_fixedKnots( const ego object );
extern "C" int EG_isPlanar( const ego object );
extern "C" int EG_getPlane( const ego object, ego *plane );
extern "C" int EG_spline2dAppx( ego context,     int    endc,
                                /*@null@*/ const double *uknot,
                                /*@null@*/ const double *vknot,
                                /*@null@*/ const int    *vdata,
                                /*@null@*/ const double *wesT,
                                /*@null@*/ const double *easT,
                                /*@null@*/ const double *south,
                                /*@null@*/       double *snor,
                                /*@null@*/ const double *north,
                                /*@null@*/       double *nnor, int imax, int jmax,
                                           const double *xyz,  double tol,
                                ego *esurf );

extern "C" int EG_spline1dEval( int *ivec, double *rdata, double t,
                                double *point );
template<int N, class T>
int EG_spline1dEval( int *ivec, SurrealS<N> *data, T& t, SurrealS<N> *point );
extern "C" int EG_spline1dDeriv( int *ivec, double *rdata, int der, double t,
                                 double *deriv );
template<int N, class T>
int EG_spline1dDeriv( int *ivec, SurrealS<N> *data, int der, T& t,
                      SurrealS<N> *point );
template<class T>
int
EG_spline1dFit(int endx, int imaxx, const T *xyz, const T *kn,
               double tol, int *ivec, T **rdata);

template<int N, class T>
int EG_spline2dEval( int *ivec, SurrealS<N> *data, const T *uv,
                     SurrealS<N> *point );
extern "C" int EG_spline2dDeriv( int *ivec, double *data, int der,
                                 const double *uv, double *deriv );
template<int N, class T>
int EG_spline2dDeriv( int *ivec, SurrealS<N> *data, int der, const T *uv,
                      SurrealS<N> *deriv );

template<class T>
int
EG_spline2dAppr(int endc, int imaxx, int jmaxx, const T *xyz,
                /*@null@*/ const T   *uknot, /*@null@*/ const T *vknot,
                /*@null@*/ const int *vdata,
                /*@null@*/ const T   *wesT,  /*@null@*/ const T *easT,
                /*@null@*/ const T   *south, /*@null@*/       T *snor,
                /*@null@*/ const T   *north, /*@null@*/       T *nnor,
                double tol, int *header, T **rdata);
template<class T>
int
EG_isoCurve( const int *header2d, const T *data2d,
             const int ik, const int jk, int *header, T **data );
template<class T>
int
EG_subSpline1d( const int *fheader, const T *fdata,
                const int i1, const int iN,
                int *sheader, T **sdata );
template<class T>
int
EG_subSpline2d( const int *fheader, const T *fdata,
                const int i1, const int iN,
                const int j1, const int jN, int *sheader, T **sdata );

#ifdef EGADS_SPLINE_VELS
extern "C" void EG_getGeometryLen( const egObject *geom, int *nivec, int *nrvec );
#endif

namespace // private to this file (no-name namespace)
{

template<class T>
struct egSequ
{
  int ncp;                   /* number of control points */
  int nave;                  /* number used to average positions */
  T   *knots;                /* the knot positions */
};

template<class T>
struct egSpline
{
  int header[7];             /* spline header information */
  T   *data;                 /* spline data */

  egSpline() : data(NULL)
  {
    for (int i = 0; i < 7; i++) header[i] = 0;
  }
  ~egSpline() { EG_free(data); }

  void clear()
  {
    for (int i = 0; i < 7; i++) header[i] = 0;
    EG_free(data); data = NULL;
  }

  egSpline& operator=(const egSpline& spline)
  {
    EG_free(data); data = NULL;
    for (int i = 0; i < 7; i++) header[i] = spline.header[i];

    int icp   = header[2];
    int iknot = header[3];
    int jcp   = header[5];
    int jknot = header[6];

    int len = 0;
    if (jcp > 0)
      len = iknot+jknot+3*icp*jcp;
    else
      len = iknot + 3*icp;

    data = (T*)EG_alloc(len*sizeof(T));
    for (int i = 0; i < len; i++) data[i] = spline.data[i];

    return *this;
  }
};

#if defined(DEBUG) || defined(SPLINE_TECPLOT_DEBUG)
static double value(double val)
{
  return val;
}

template <int N>
static double value(SurrealS<N> valS)
{
  return valS.value();
}
#endif

#ifdef SPLINE_TECPLOT_DEBUG
extern "C"
int
EG_spline2dEval(int *ivec, double *data, const double *uv, double *point);

template<class T>
static void
EG_tecplotSpline(int *header, T *data, const char *filename)
{
  T uv[2], point[9];
  int ideg  = header[1];
  int icp   = header[2];
  int iknot = header[3]-2*ideg;
  int jdeg  = header[4];
  int jcp   = header[5];
  int jknot = header[6]-2*jdeg;

  int imax  = iknot + ideg*(iknot-1);
  int jmax  = jknot + jdeg*(jknot-1);

  T *knotu =  data+ideg;
  T *knotv = &data[header[3]+jdeg];
  T *cp    = &data[header[3]+header[6]];

  printf( "Writing %s...\n", filename );
  FILE* fp = fopen( filename, "w" );

  fprintf( fp, "\"\"\n" );
  fprintf( fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"Xu\", \"Yu\", \"Zu\", \"Xv\", \"Yv\", \"Zv\", \"u\", \"v\"\n" );
  fprintf( fp, "ZONE T=\"%s\", I=%d, J=%d\n", filename, imax, jmax );

  for (int j = 0; j < jknot-1; j++) {
    for (int jj = (j == 0 ? 0 : 1); jj <= jdeg+1; jj++) {

      uv[1] = knotv[j]*(1.-jj/double(jdeg+1)) + knotv[j+1]*(jj/double(jdeg+1));

      for (int i = 0; i < iknot-1; i++) {
        for (int ii = (i == 0 ? 0 : 1); ii <= ideg+1; ii++) {

          uv[0] = knotu[i]*(1.-ii/double(ideg+1)) + knotu[i+1]*(ii/double(ideg+1));

          EG_spline2dDeriv(header, data, 1, uv, point);

          for (int k = 0; k < 9; k++)
            fprintf( fp, "%22.15e " , value(point[k]) );

          fprintf( fp, "%22.15e " , value(uv[0]) );
          fprintf( fp, "%22.15e\n", value(uv[1]) );
        }
      }
    }
  }

  fprintf( fp, "\"\"\n" );
  fprintf( fp, "ZONE T=\"%s cp\", I=%d, J=%d\n", filename, icp, jcp );

  for (int j = 0; j < jcp; j++) {
    for (int i = 0; i < icp; i++) {
      fprintf( fp, "%22.15e ", value(cp[3*(i+j*icp)  ]) );
      fprintf( fp, "%22.15e ", value(cp[3*(i+j*icp)+1]) );
      fprintf( fp, "%22.15e ", value(cp[3*(i+j*icp)+2]) );
      for (int k = 0; k < 8; k++)
        fprintf( fp, "0 " );
      fprintf( fp, "\n" );
    }
  }

  fclose(fp);
}

static void
EG_tecplotSplinePoints(int imax, int jmax, const SurrealS<1> *xyz,
                       const char *filename) {}

static void
EG_tecplotSplinePoints(int imax, int jmax, const double *xyzs,
                       const char *filename)
{
  printf( "Writing %s...\n", filename );
  FILE* fp = fopen( filename, "w" );

  fprintf( fp, "\"\"\n" );
  fprintf( fp, "VARIABLES = \"X\", \"Y\", \"Z\"\n" );
  fprintf( fp, "ZONE T=\"%s\", I=%d, J=%d\n", filename, imax, jmax );

  for (int j = 0; j < jmax; j++) {
    for (int i = 0; i < imax; i++) {
      for (int k = 0; k < 3; k++)
        fprintf( fp, "%22.15e ", xyzs[3*(i+j*imax)+k]);
      fprintf( fp, "\n" );
    }
  }

  fclose(fp);
}
#endif


#ifdef DUMP_SECTIONS
static void
EG_dumpSections(int nsec, const ego *secs, const char *filename)
{
  int    stat = EGADS_SUCCESS;
  int    sense[1] = {SFORWARD};
  double tdata[2];
  ego *bodies, model, eedge, eloop, context;

  context = EG_context(secs[0]);

  bodies = (ego*)EG_alloc(nsec*sizeof(ego));
  for (int i = 0; i < nsec; i++) bodies[i] = NULL;

  for (int i = 0; i < nsec; i++) {
    if (secs[i]->oclass == NODE) {
      tdata[0] = 0;
      tdata[1] = 1;
      stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, tdata, 1,
                             (ego *) &secs[i], sense, &eedge);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, 1, &eedge,
                             sense, &eloop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &eloop,
                             NULL, &bodies[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      EG_deleteObject(eloop);
      EG_deleteObject(eedge);
    } else if (secs[i]->oclass == LOOP) {

      stat = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1,
                             (ego *) &secs[i], NULL, &bodies[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;

    } else if (secs[i]->oclass == FACE) {

      stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1,
                             (ego *) &secs[i], NULL, &bodies[i]);
      if (stat != EGADS_SUCCESS) goto cleanup;

    }
  }

  stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, nsec, bodies,
                         NULL, &model);
  if (stat != EGADS_SUCCESS) goto cleanup;

  printf(" EGADS Debug: Writing %s\n", filename);
  stat = EG_saveModel(model, filename);
  EG_deleteObject(model);

cleanup:
  for (int i = 0; i < nsec; i++)
    EG_deleteObject(bodies[i]);

  if (stat != EGADS_SUCCESS) {
    printf(" EGADS Warning: Failed to write %s (EG_dumpSections)!\n", filename);
  }
}
#endif


static int
EG_findTE(int nstripe, ego *edges, int *te)
{
  int    stat = EGADS_SUCCESS;
  int    j, nnode, oclass, mtype, *senses;
  double mlen, nlen, limits[4];
  ego    geom, *nodes;

  *te = -1;
  if (nstripe == 2) return EGADS_SUCCESS;

  /* find the te Edge */
  if (nstripe == 3) {
    mlen = 1.e100;
    for (j = 0; j < nstripe; j++) {
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, limits, &nnode,
                            &nodes, &senses);
      if (stat != EGADS_SUCCESS) return stat;

      stat = EG_arcLength(edges[j], limits[0], limits[1], &nlen);
      if (stat != EGADS_SUCCESS) return stat;
      if (nlen > mlen) continue;
      mlen = nlen;
      *te   = j;
    }
  }

  return EGADS_SUCCESS;
}


/* knot sequence utility functions */

template<class T>
static int
EG_allocSeq(int nstripe, egSequ<T> **sequ)
{
  int    i;
  egSequ<T> *seq;

  *sequ = NULL;
  seq   = (egSequ<T> *) EG_alloc(nstripe*sizeof(egSequ<T>));
  if (seq == NULL) return EGADS_MALLOC;
  for (i = 0; i < nstripe; i++) {
    seq[i].ncp   = 0;
    seq[i].nave  = 0;
    seq[i].knots = NULL;
  }

  *sequ = seq;
  return EGADS_SUCCESS;
}


template<class T>
static void
EG_freeSeq(int nstripe, /*@only@*/ egSequ<T> *seq)
{
  int i;

  for (i = 0; i < nstripe; i++)
    if (seq[i].knots != NULL) EG_free(seq[i].knots);

  EG_free(seq);
}


template<class T>
static int
EG_mergeSeq(int stripe, egSequ<T> *seq, int n, T *nomulti)
{
  int i;
  T   dt;

  /*  for (i = 0; i < n; i++) printf(" %d/%d: %lf\n", i+1, n, nomulti[i]);  */
  /* keep existing if larger */
  if (seq[stripe].ncp > n) {
    return EGADS_SUCCESS;
  }

  /* same number of knots -- sum with previous */
  if (seq[stripe].ncp == n) {
    if (seq[stripe].nave == 0) { /* keep equi-spaced */
      return EGADS_SUCCESS;
    }
    for (i = 0; i < n; i++) {
      seq[stripe].knots[i] *= seq[stripe].nave;
      seq[stripe].knots[i] += (nomulti[i  ] - nomulti[0])/
                              (nomulti[n-1] - nomulti[0]);
    }
    dt = seq[stripe].knots[n-1];
    for (i = 0; i < n; i++) seq[stripe].knots[i] /= dt;
    seq[stripe].nave++;
    return EGADS_SUCCESS;
  }

  /* greater sequence -- use it */
  if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
  seq[stripe].knots = (T *) EG_alloc(n*sizeof(T));
  if (seq[stripe].knots == NULL) {
    return EGADS_MALLOC;
  }
  for (i = 0; i < n; i++)
    seq[stripe].knots[i] = nomulti[i] - nomulti[0];
  dt = seq[stripe].knots[n-1];
  for (i = 0; i < n; i++) seq[stripe].knots[i] /= dt;
  seq[stripe].nave = 1;
  seq[stripe].ncp  = n;

  return EGADS_SUCCESS;
}


template<class T>
static int
EG_setSeq(int stripe, egSequ<T> *seq, int num, int sens,
          T *range, /*@null@*/ int *iinfo, /*@null@*/ T *rinfo)
{
  int stat = EGADS_SUCCESS;
  int i, nk, n, ndeg;
  T   dt, *nomulti, *tmp;

  if ((iinfo == NULL) || (rinfo == NULL)) {
    /* not a spline */
    if (seq[stripe].ncp >= num) return EGADS_SUCCESS;
    if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
    seq[stripe].knots = (T *) EG_alloc(num*sizeof(T));
    if (seq[stripe].knots == NULL) return EGADS_MALLOC;
    /* equi-spaced sampling */
    dt = 1.0/(num-1.0);
    for (i = 0; i < num; i++) seq[stripe].knots[i] = i*dt;
    seq[stripe].nave = 0;
    seq[stripe].ncp  = num;
  } else {
    ndeg = iinfo[1];
    nk   = iinfo[3] - 2*ndeg;
    nomulti = (T *) EG_alloc(nk*sizeof(T));
    if (nomulti == NULL) return EGADS_MALLOC;

    n = 0;
    nomulti[n++] = range[0];
    for (i = 0; i < nk; i++) {

      if (rinfo[i+ndeg] < range[0]) continue;
      if (rinfo[i+ndeg] > range[1]) break;

      if (n == 1) {
        if (rinfo[i+ndeg] - range[0] > 1e-4)
          nomulti[n++] = rinfo[i+ndeg];
      } else {
        nomulti[n] = rinfo[i+ndeg];
        if (nomulti[n] != nomulti[n-1]) n++;
      }
    }
    if (range[1] - nomulti[n-1] > 1e-4) {
      nomulti[n++] = range[1];
    } else {
      nomulti[n-1] = range[1];
    }

    if ((seq[stripe].ncp < num) && (n < num)) {
      EG_free(nomulti);
      if (seq[stripe].knots != NULL) EG_free(seq[stripe].knots);
      seq[stripe].knots = (T *) EG_alloc(num*sizeof(T));
      if (seq[stripe].knots == NULL) return EGADS_MALLOC;
      /* equi-spaced sampling */
      dt = 1.0/(num-1.0);
      for (i = 0; i < num; i++) seq[stripe].knots[i] = i*dt;
      seq[stripe].nave = 0;
      seq[stripe].ncp  = num;
      return EGADS_SUCCESS;
    } else if ( sens == -1 ) {

      tmp = (T *) EG_alloc(n*sizeof(T));
      if (tmp == NULL) return EGADS_MALLOC;
      for (i = 0; i < n; i++)
        tmp[i] = nomulti[i];
      
      /* reverse the spacing for a negative sense */
      range[1] = nomulti[n-1];
      for (i = 0; i < n; i++)
        nomulti[i] = range[1] - tmp[n-1-i];
      EG_free(tmp);
    }

    /* merge the current sequence with possibly existing one */
    stat = EG_mergeSeq(stripe, seq, n, nomulti);
    EG_free(nomulti);
  }

  return stat;
}


#ifdef EGADS_SPLINE_VELS
static void
setTrange(double *trange, double *trangeD, double *trangeD_dot)
{
  trange[0] = trangeD[0];
  trange[1] = trangeD[1];
}

static void
setTrange(SurrealS<1> *trange, double *trangeD, double *trangeD_dot)
{
  trange[0].value() = trangeD[0]; trange[0].deriv() = trangeD_dot[0];
  trange[1].value() = trangeD[1]; trange[1].deriv() = trangeD_dot[1];
}

static void
setRinfo(int len, double **rvec, double *rvecD, double *rvecD_dot)
{
  (*rvec) = (double*)EG_alloc(len*sizeof(double));
  for (int i = 0; i < len; i++)
    (*rvec)[i] = rvecD[i];
}

static void
setRinfo(int len, SurrealS<1> **rvec, double *rvecD, double *rvecD_dot)
{
  (*rvec) = (SurrealS<1>*)EG_alloc(len*sizeof(SurrealS<1>));
  for (int i = 0; i < len; i++) {
    (*rvec)[i].value() = rvecD[i];
    (*rvec)[i].deriv() = rvecD_dot[i];
  }
}
#endif


template<class T>
static int
EG_setSequence(
#ifdef EGADS_SPLINE_VELS
    const egadsSplineVels *vels,
#endif
    int nsec, const ego *secs, int te, int nstripe,
               egSequ<T> **ncp_out)
{
  int    i, j, k, n, jj, outLevel, stat = EGADS_SUCCESS;
  int    oclass, mtype, nnode, *senses=NULL, *iinfo=NULL, *itmp=NULL;
  double data[18];
#ifdef EGADS_SPLINE_VELS
  double trangeD[2], trangeD_dot[2], *rinfoD, *rinfoD_dot;
  ego    top, prev, next;
#endif
  T      *rinfo=NULL, trange[2];
  ego    loop=NULL, ref, geom;
  ego    *chldrn=NULL, *edges=NULL, *nodes=NULL;
  egSequ<T> *ncp=NULL;

  *ncp_out = NULL;

  outLevel = EG_outLevel(secs[0]);

  /* get the number of sample points per Edge */
  stat = EG_allocSeq(nstripe, &ncp);
  if (stat !=  EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Edges (EG_setSequence)!\n", nstripe);
    goto cleanup;
  }
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      if (mtype == WIREBODY) {
        loop = chldrn[0];
      } else if (mtype == FACEBODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &n, &chldrn,
                              &senses);
        loop = chldrn[0];
      }
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (j = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nnode,
                            &nodes, &itmp);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTopo = %d (EG_setSequence)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      if (mtype == DEGENERATE) continue;

#ifdef EGADS_SPLINE_VELS

      if (vels != NULL && vels->velocityOfRange != NULL) {

        stat = (*(vels->velocityOfRange))(vels->usrData, secs, i, edges[jj],
                                          trangeD, trangeD_dot);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d velocityOfRange = %d (EG_setSequence)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
        setTrange(trange, trangeD, trangeD_dot);

      } else if (EG_hasGeometry_dot(edges[jj]) == EGADS_SUCCESS) {

        stat = EG_getRange(edges[jj], trange, &k);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_setSequence)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }

      } else {

        stat = EG_getRange(edges[jj], trangeD, &k);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_setSequence)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
        trange[0] = trangeD[0];
        trange[1] = trangeD[1];
      }
#else
      stat = EG_getRange(edges[jj], trange, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d getRange = %d (EG_ruled)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
#endif


      mtype = -1;
      do {
        if ((mtype == BEZIER) || (mtype == BSPLINE)) {
          EG_free(iinfo);
          iinfo = NULL;
        }
        geom = ref;
#ifdef EGADS_SPLINE_VELS
        stat = EG_getInfo(geom, &oclass, &mtype, &top, &prev, &next);
        if (stat != EGADS_SUCCESS) goto cleanup;

        if (mtype == BSPLINE && vels != NULL && vels->velocityOfBspline != NULL) {
          stat = (*(vels->velocityOfBspline))(vels->usrData, secs, i, edges[jj], geom,
                                              &iinfo, &rinfoD, &rinfoD_dot);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d velocityOfBspline = %d (EG_setSequence)!\n",
                     i+1, jj+1, stat);
            goto cleanup;
          }
          stat = EG_setGeometry_dot(geom, oclass, mtype, iinfo, rinfoD, rinfoD_dot);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d EG_setGeometry_dot = %d (EG_setSequence)!\n",
                     i+1, jj+1, stat);
            goto cleanup;
          }

          int ilen, len;
          EG_getGeometryLen(geom, &ilen, &len);
          setRinfo(len, &rinfo, rinfoD, rinfoD_dot);
          EG_free(rinfoD);     rinfoD = NULL;
          EG_free(rinfoD_dot); rinfoD_dot = NULL;
        } else if (EG_hasGeometry_dot(geom) == EGADS_SUCCESS) {
          stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rinfo);
        } else {
          double *rvec;
          int ilen, len;
          stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rvec);
          if (stat != EGADS_SUCCESS) goto cleanup;
          EG_getGeometryLen(geom, &ilen, &len);
          rinfo = (T*)EG_alloc(len*sizeof(T));
          for (k = 0; k < len; k++) rinfo[k] = rvec[k];
          EG_free(rvec);
        }
#else
        stat = EG_getGeometry(geom, &oclass, &mtype, &ref, &iinfo, &rinfo);
#endif
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d getGeom = %d (EG_setSequence)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }
        if (mtype != BSPLINE) { EG_free(rinfo); rinfo = NULL; }
      } while ((mtype == TRIMMED) || (mtype == OFFSET));

      stat = EGADS_NOTFOUND;
      if (mtype == LINE) {
        stat = EG_setSeq<T>(j, ncp, 3, senses[jj], trange, NULL, NULL);
      } else if (mtype == CIRCLE) {
/*      This fixes windtunnel6 and rule8 for OCC 7.4
        k = 36.*(value(trange[1])-value(trange[0]))/PI;
        if (k < 4) k = 4;
        stat = EG_setSeq<T>(j, ncp,  k, trange, NULL, NULL);
 */
        stat = EG_setSeq<T>(j, ncp, 12, senses[jj], trange, NULL, NULL);
      } else if (mtype == ELLIPSE) {
/*      This fixes windtunnel6 and rule8 for OCC 7.4
        k = 36.*(value(trange[1])-value(trange[0]))/PI;
        if (k < 4) k = 4;
        stat = EG_setSeq<T>(j, ncp,  k, trange, NULL, NULL);
*/
        stat = EG_setSeq<T>(j, ncp, 12, senses[jj], trange, NULL, NULL);
      } else if (mtype == PARABOLA) {
        stat = EG_setSeq<T>(j, ncp, 12, senses[jj], trange, NULL, NULL);
      } else if (mtype == HYPERBOLA) {
        stat = EG_setSeq<T>(j, ncp, 12, senses[jj], trange, NULL, NULL);
      } else if (mtype == BEZIER) {
        k    = iinfo[2];
        if (k < 12) k = 12;
        stat = EG_setSeq<T>(j, ncp, k, senses[jj], trange, NULL, NULL);
        EG_free(iinfo); iinfo = NULL;
      } else if (mtype == BSPLINE) {
        stat = EG_setSeq<T>(j, ncp, 12, senses[jj], trange, iinfo, rinfo);
        EG_free(iinfo); iinfo = NULL;
        EG_free(rinfo); rinfo = NULL;
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d Edge %d setSeq = %d (EG_setSequence)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      j++;
    }
  }

  if (te > 0) {

    /* make sure sampling sequences are identical */
    if (nstripe == 2 || te == 2) {
      i = 0; j = 1;
    } else if (te == 0) {
      i = 1; j = 2;
    } else {
      i = 0; j = 2;
    }

    stat = EG_mergeSeq(i, ncp, ncp[j].ncp, ncp[j].knots);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_mergeSeq(j, ncp, ncp[i].ncp, ncp[i].knots);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat == EGADS_SUCCESS)
    *ncp_out = ncp;
  else
    EG_freeSeq(nstripe, ncp);

  return stat;
}


template<class T>
static void
EG_checkDirs(ego edge, T t, T *pnt, T *dir2)
{
  int stat;
  T   dot, data[9], tan[3], dir[3];

  stat = EG_evaluate(edge, &t, data);
  if (stat != EGADS_SUCCESS) return;
  dir[0]   = (data[0]-pnt[0]);
  dir[1]   = (data[1]-pnt[1]);
  dir[2]   = (data[2]-pnt[2]);
  dot      = sqrt(DOT(dir, dir));
  if (dot != 0.0) {
    dir[0] /= dot;
    dir[1] /= dot;
    dir[2] /= dot;
  }
  tan[0]    = dir2[0];
  tan[1]    = dir2[1];
  tan[2]    = dir2[2];
  dot       = sqrt(DOT(dir2, dir2));
  if (dot != 0.0) {
    tan[0] /= dot;
    tan[1] /= dot;
    tan[2] /= dot;
  }
  if (DOT(dir, tan) < 0.0) {
//  printf(" checkDir flip = %lf!\n", DOT(dir, tan));
    dir2[0] = -dir2[0];
    dir2[1] = -dir2[1];
    dir2[2] = -dir2[2];
  }
}


template<class T>
static int
EG_wingTipSpline(const int outLevel, const T &ratio, ego sect, ego seci,
                 double v, egSpline<T> **surfs, egSpline<T> *tipsurf)
{
  int    stat, i, j, k, kk, n, oclass, mtype, nFace, nLoop, nEdge, nNode, nKnotu, nKnotv;
  int    *info, *senses, *lsenses=NULL, open, deg, iknot, nsub;
  double *gdata=NULL;
  double nlen, tol, limits[4], norm[3], snorm[3], xu[3], xv[3];
  T      t, s, uv[2], alen, blen, b0norm, x0[3], x1[3], a0[3], b0[3];
  T      *data, *xyz=NULL, *uKnots=NULL, *vKnots=NULL;
  T      *west=NULL, *east=NULL, *snor=NULL, *nnor=NULL, results[18];
  ego    geom, rGeom, loop;
  ego    *faces, *loops, *edges, *nodes;

  if (sect->oclass != FACE && sect->mtype != FACEBODY) {
    if (outLevel > 0)
      printf(" EGADS Error: Section is not a Face or FACEBODY (EG_wingTip)!\n");
    return EGADS_GEOMERR;
  }

  /* get the geometry for the tip section */
  stat = EG_getTopology(sect, &geom, &oclass, &mtype, limits, &nLoop,
                        &loops, &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getTopology = %d for Section (EG_wingTip)!\n",
             stat);
    return EGADS_TOPOERR;
  }
  if (oclass == BODY) {
    stat = EG_getTopology(loops[0], &geom, &oclass, &mtype, limits, &nLoop,
                          &loops, &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getTopology = %d for Section (EG_wingTip)!\n",
               stat);
      return EGADS_TOPOERR;
    }
  }


  /* look at the loops */
  if (nLoop != 1) {
    if (outLevel > 0)
      printf(" EGADS Error: Section has %d Loops (EG_wingTip)!\n",
             nLoop);
    return EGADS_TOPOERR;
  }

  /* get the plane normal */
  stat = EG_getGeometry(geom, &oclass, &mtype, &rGeom, &info, &gdata);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getGeometry = %d (EG_wingTip)!\n",
             stat);
    goto cleanup;
  }
  if (mtype != PLANE) {
    if (outLevel > 0)
      printf(" EGADS Error: Section Face is not Planar (EG_wingTip)!\n");
    goto cleanup;
  }
  xu[0]    = gdata[3];
  xu[1]    = gdata[4];
  xu[2]    = gdata[5];
  xv[0]    = gdata[6];
  xv[1]    = gdata[7];
  xv[2]    = gdata[8];
  CROSS(norm, xu, xv);
  nlen     = sqrt(DOT(norm, norm));
  norm[0] /= nlen;
  norm[1] /= nlen;
  norm[2] /= nlen;

  EG_free(info);  info  = NULL;
  EG_free(gdata); gdata = NULL;

  /* extract the first node from the tip section */
  stat = EG_getTopology(loops[0], &geom, &oclass, &mtype, limits, &nEdge,
                        &edges, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_getTopology(edges[0], &geom, &oclass, &mtype, limits, &nNode,
                        &nodes, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_evaluate(nodes[0], NULL, xu);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* get the loop for the interior section */
  if (seci->oclass == BODY) {
    stat = EG_getTopology(seci, &geom, &oclass, &mtype, limits, &nFace,
                          &faces, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    stat = EG_getTopology(faces[0], &geom, &oclass, &mtype, limits, &nLoop,
                          &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    loop = loops[0];
  } else  if (seci->oclass == FACE) {
    stat = EG_getTopology(seci, &geom, &oclass, &mtype, limits, &nLoop,
                          &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    loop = loops[0];
  } else {
    loop = seci;
  }

  stat = EG_getTopology(loop, &geom, &oclass, &mtype, limits, &nEdge,
                        &edges, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_getTopology(edges[0], &geom, &oclass, &mtype, limits, &nNode,
                        &nodes, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  stat = EG_evaluate(nodes[0], NULL, xv);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* check if the tip section plane normal needs to be flipped */
  snorm[0] = xu[0] - xv[0];
  snorm[1] = xu[1] - xv[1];
  snorm[2] = xu[2] - xv[2];
  nlen     = sqrt(DOT(snorm, snorm));
  snorm[0] /= nlen;
  snorm[1] /= nlen;
  snorm[2] /= nlen;

  if (DOT(norm, snorm) < 0.0) {
    norm[0] = -norm[0];
    norm[1] = -norm[1];
    norm[2] = -norm[2];
  }

  /* make the grid of points. nKnotu must be odd. */
  nKnotu = 2*NTIP+1;
#ifdef MATCH_TIP_SNOR_NNOR
  nsub = 0;
#else
  nsub = 3;
#endif

  /* get the Knots */
  tol    = 1e-7;
  deg    = surfs[0]->header[1];
  iknot  = surfs[0]->header[3]-2*deg;
  data   = surfs[0]->data;
  nKnotv = iknot + nsub*(iknot-1);
  uKnots = (T *) EG_alloc(nKnotu*sizeof(T));
  vKnots = (T *) EG_alloc(nKnotv*sizeof(T));
  if ((uKnots == NULL) || (vKnots == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Knot Alloc %d %d (EG_wingTip)!\n", nKnotu, nKnotv);
    stat = EGADS_MALLOC;
    goto cleanup;
  }

#ifdef MATCH_TIP_SNOR_NNOR
  snor = (T *) EG_alloc(3*nKnotu*sizeof(T));
  nnor = (T *) EG_alloc(3*nKnotu*sizeof(T));
  if ((snor == NULL) || (nnor == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Knot Alloc %d (EG_wingTip)!\n",
             info[3]);
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for (k = 0; k < 3*nKnotu; k++) snor[k] = nnor[k] = 0.0;
#endif

  // equally spaced in u to match the trailing edge face
  for (k = 0; k < nKnotu; k++) uKnots[k] = k/double(nKnotu-1);

  //for (k = 0; k < nKnotv; k++) vKnots[k] = data[k+3];
  n = 0;
  vKnots[n++] = data[deg];
  for (k = 0; k < iknot-1; k++)
    for (kk = 1; kk <= nsub+1; kk++)
      vKnots[n++] = data[k  +deg]*(1.-kk/double(nsub+1)) +
                    data[k+1+deg]*(kk/double(nsub+1));

  xyz = (T *) EG_alloc((nKnotu+2)*3*nKnotv*sizeof(T));
  if (xyz == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: %d x 3 Allocate (EG_wingTip)!\n",
             nKnotv);
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  west = &xyz[(nKnotu  )*3*nKnotv];
  east = &xyz[(nKnotu+1)*3*nKnotv];

  n = 0;
  for (j = 0; j < 2; j++) {

    for (k = 0; k < nKnotv; k++) {
      t = vKnots[k];
      if (n != 0) t = 1.0 - vKnots[k];

      uv[0] = t;
      uv[1] = v;
      stat = EG_spline2dDeriv(surfs[j]->header, surfs[j]->data, 1, uv, results);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: %d/%d Eval = %d (EG_wingTip)!\n",
                 k, nKnotv, stat);
        goto cleanup;
      }

      xyz[(k*nKnotu+n)*3  ] = results[0];
      xyz[(k*nKnotu+n)*3+1] = results[1];
      xyz[(k*nKnotu+n)*3+2] = results[2];

#ifdef MATCH_TIP_SNOR_NNOR
      if (k == 0 && j == 0) {
        snor[3*n  ] =  results[3];
        snor[3*n+1] =  results[4];
        snor[3*n+2] =  results[5];
      } else if (k == nKnotv-1 && j == 0) {
        nnor[3*n  ] = -results[3];
        nnor[3*n+1] = -results[4];
        nnor[3*n+2] = -results[5];
      } else if (k == 0 && j == 1) {
        snor[3*n  ] = -results[3];
        snor[3*n+1] = -results[4];
        snor[3*n+2] = -results[5];
      } else if (k == nKnotv-1 && j == 1) {
        nnor[3*n  ] =  results[3];
        nnor[3*n+1] =  results[4];
        nnor[3*n+2] =  results[5];
      }
#endif

      x1[0] = results[6];
      x1[1] = results[7];
      x1[2] = results[8];
      alen  = sqrt(DOT(x1,x1));
      if (alen != 0.0) {
        x1[0] /= alen;
        x1[1] /= alen;
        x1[2] /= alen;
      }
      alen  = 1.0;
      if (n == 0) {
        if (DOT(x1, norm) < 0.0) alen = -1.0;
        west[3*k  ] = alen*x1[0];
        west[3*k+1] = alen*x1[1];
        west[3*k+2] = alen*x1[2];
      } else {
        if (DOT(x1, norm) > 0.0) alen = -1.0;
        east[3*k  ] = alen*x1[0];
        east[3*k+1] = alen*x1[1];
        east[3*k+2] = alen*x1[2];
      }
    }
    n = nKnotu-1;
  }

#ifdef MAXIMUM_RATIO_WARNING
  double w[3], e[3], A[4], b[2], denom, st[2], X0[3], X1[3], Xm[3], dX[3], r;
  double ratioMax = value(ratio);
  n = nKnotu-1;
  for (k = 0; k < nKnotv; k++) {

    X0[0] = value(xyz[(k*nKnotu  )*3  ]);
    X0[1] = value(xyz[(k*nKnotu  )*3+1]);
    X0[2] = value(xyz[(k*nKnotu  )*3+2]);
    X1[0] = value(xyz[(k*nKnotu+n)*3  ]);
    X1[1] = value(xyz[(k*nKnotu+n)*3+1]);
    X1[2] = value(xyz[(k*nKnotu+n)*3+2]);
    w[0]  = value(-west[3*k  ]);
    w[1]  = value(-west[3*k+1]);
    w[2]  = value(-west[3*k+2]);
    e[0]  = value(-east[3*k  ]);
    e[1]  = value(-east[3*k+1]);
    e[2]  = value(-east[3*k+2]);

    /* J  = [[w[0], e[0]],
             [w[1], e[1]],
             [w[2], e[2]]]
       solve  A  = (J' * J)
          and b  =  J' * (x0 - x1)
       where  J' = transpose(J) */
    A[0] = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
    A[1] = w[0]*e[0] + w[1]*e[1] + w[2]*e[2];
    A[2] = A[1];
    A[3] = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];

    b[0] = w[0]*(X0[0] - X1[0]) + w[1]*(X0[1] - X1[1]) + w[2]*(X0[2] - X1[2]);
    b[1] = e[0]*(X0[0] - X1[0]) + e[1]*(X0[1] - X1[1]) + e[2]*(X0[2] - X1[2]);

    /* solve A*delta = b using Cramer's rule */
    denom    =  A[0]*A[3] - A[2]*A[1];
    if (denom == 0.0) continue; /* this means the vectors are parallel,
                                   which is not a problem */
    st[0] = (b[0]*A[3] - b[1]*A[1])/denom;
    st[1] = (A[0]*b[1] - A[2]*b[0])/denom;

    /* if both are positive, the vectors point inwards and could cross over */
    if ((st[0] > 0) && (st[1] > 0)) {

      Xm[0] = 0.5*(X0[0] + X1[0]);
      Xm[1] = 0.5*(X0[1] + X1[1]);
      Xm[2] = 0.5*(X0[2] + X1[2]);

      dX[0] = X0[0] - X1[0];
      dX[1] = X0[1] - X1[1];
      dX[2] = X0[2] - X1[2];

      r = 0.5*sqrt(DOT(dX,dX));

      X0[0] = X0[0] - w[0]*st[0];
      X0[1] = X0[1] - w[1]*st[0];
      X0[2] = X0[2] - w[2]*st[0];
      X1[0] = X1[0] + e[0]*st[1];
      X1[1] = X1[1] + e[1]*st[1];
      X1[2] = X1[2] + e[2]*st[1];

      dX[0] = X0[0] - Xm[0];
      dX[1] = X0[1] - Xm[1];
      dX[2] = X0[2] - Xm[2];

      ratioMax = MIN( ratioMax, sqrt(DOT(dX,dX))/r );

      dX[0] = X1[0] - Xm[0];
      dX[1] = X1[1] - Xm[1];
      dX[2] = X1[2] - Xm[2];

      ratioMax = MIN( ratioMax, sqrt(DOT(dX,dX))/r );
    }
  }

  if ((ratio > ratioMax) && (ratioMax > 0.0)) {
    printf("EGADS Warning: Ratio %lf > %lf might give cusped or folded shape (EG_wingTip)!\n",
           value(ratio), ratioMax);
  }
#endif

  n = nKnotu-1;
  for (open = k = 0; k < nKnotv; k++) {
    x0[0] = 0.5*(xyz[(k*nKnotu  )*3  ] + xyz[(k*nKnotu+n)*3  ]);
    x0[1] = 0.5*(xyz[(k*nKnotu  )*3+1] + xyz[(k*nKnotu+n)*3+1]);
    x0[2] = 0.5*(xyz[(k*nKnotu  )*3+2] + xyz[(k*nKnotu+n)*3+2]);
    a0[0] = 0.5*(xyz[(k*nKnotu  )*3  ] - xyz[(k*nKnotu+n)*3  ]);
    a0[1] = 0.5*(xyz[(k*nKnotu  )*3+1] - xyz[(k*nKnotu+n)*3+1]);
    a0[2] = 0.5*(xyz[(k*nKnotu  )*3+2] - xyz[(k*nKnotu+n)*3+2]);
    alen  = sqrt(DOT(a0,a0));
    if (alen != 0.0) {
      a0[0] /= alen;
      a0[1] /= alen;
      a0[2] /= alen;
    }

    if ((k == 0)        && (alen > tol)) open = -1;
    if ((k == nKnotv-1) && (alen > tol)) open =  1;
    if ((k == 0)        && (open != -1)) alen = 0.0;
    if ((k == nKnotv-1) && (open !=  1)) alen = 0.0;
    b0[0] = 0.5*(west[3*k  ] - east[3*k  ]);
    b0[1] = 0.5*(west[3*k+1] - east[3*k+1]);
    b0[2] = 0.5*(west[3*k+2] - east[3*k+2]);
    b0norm = sqrt(DOT(b0, b0));
    if (b0norm == 0.0) b0norm = 1.0;
    blen = alen*ratio/b0norm;

    /* create the tip ellipse */
    for (i = 1; i < nKnotu-1; i++) {
      s = i/double(nKnotu-1);
      xyz[(k*nKnotu+i)*3  ] = x0[0] + alen*a0[0]*cos(PI*s) +
                              blen*((1-s)*west[3*k  ] - s*east[3*k  ])*sin(PI*s);
      xyz[(k*nKnotu+i)*3+1] = x0[1] + alen*a0[1]*cos(PI*s) +
                              blen*((1-s)*west[3*k+1] - s*east[3*k+1])*sin(PI*s);
      xyz[(k*nKnotu+i)*3+2] = x0[2] + alen*a0[2]*cos(PI*s) +
                              blen*((1-s)*west[3*k+2] - s*east[3*k+2])*sin(PI*s);
    }

    /* scale the tangents using the magnitude of the tangent from an ellipse */
    west[3*k  ] *=  PI*blen;
    west[3*k+1] *=  PI*blen;
    west[3*k+2] *=  PI*blen;
    east[3*k  ] *= -PI*blen;
    east[3*k+1] *= -PI*blen;
    east[3*k+2] *= -PI*blen;
  }

#ifdef SPLINE_TECPLOT_DEBUG
  {
    char filename[42];
    sprintf(filename, "wingtipVectors_%d.dat", (int)(v+1));
    printf( "Writing %s...\n", filename );
    FILE* fp = fopen( filename, "w" );

    fprintf( fp, "\"\"\n" );
    fprintf( fp, "VARIABLES = \"X\", \"Y\", \"Z\", \"nx\", \"ny\", \"nz\"\n" );
    fprintf( fp, "ZONE T=\"%s\", I=%d\n", filename, 2*nKnotv );

    for (int i = 0; i < nKnotv; i++) {
      for (int k = 0; k < 3; k++)
        fprintf( fp, "%22.15e ", value(xyz[3*(i*nKnotu  )+k]));
      for (int k = 0; k < 3; k++)
        fprintf( fp, "%22.15e ", value(west[3*i+k]));
      fprintf( fp, "\n" );
    }
    for (int i = nKnotv-1; i >= 0; i--) {
      for (int k = 0; k < 3; k++)
        fprintf( fp, "%22.15e ", value(xyz[3*(i*nKnotu+n)+k]));
      for (int k = 0; k < 3; k++)
        fprintf( fp, "%22.15e ", value(east[3*i+k]));
      fprintf( fp, "\n" );
    }

    fclose(fp);

    sprintf(filename, "wingtipPoints_%d.dat", (int)(v+1));
    EG_tecplotSplinePoints(nKnotu, nKnotv, xyz, filename);
  }
#endif

  /* fit the points */
  stat = EG_spline2dAppr<T>(0, nKnotu, nKnotv, xyz, uKnots, vKnots, NULL,
                            west, east, NULL, snor, NULL, nnor,
                            1.e-8, tipsurf->header, &tipsurf->data);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_spline2dAppr = %d (EG_wingTip)!\n",
             stat);
    goto cleanup;
  }

  stat = EGADS_SUCCESS;

cleanup:
  EG_free(lsenses);
  EG_free(info);
  EG_free(xyz);
  EG_free(gdata);
  EG_free(snor);
  EG_free(nnor);
  EG_free(uKnots);
  EG_free(vKnots);

  return stat;
}


static int
EG_getSecNodes(int outLevel, int inode, ego sec, int nstripe, ego *nodes)
{
  int    stat;
  int    nnode, nedge, nloop, nface, oclass, mtype, *senses, *sens, j, uclosed;
  double limits[4];
  ego    *node = NULL, *edges, *loops, *faces, eref, loop=NULL;

  if (!(sec->oclass == BODY && (sec->mtype == WIREBODY || sec->mtype == FACEBODY)) &&
      sec->oclass != FACE && sec->oclass != LOOP && sec->oclass != NODE) {
    printf("EGADS Error: getSecNodes: sec is not a FaceBody, WireBody, Face, Loop or Node!\n");
    return EGADS_GEOMERR;
  }

  if (sec->oclass == NODE) {
    for (j = 0; j < nstripe+1; j++)
      nodes[j*inode] = sec;
    return EGADS_SUCCESS;
  } else if (sec->oclass == BODY) {
    if (sec->mtype == WIREBODY) {
      // get the loop
      stat = EG_getTopology(sec, &eref, &oclass, &mtype,
                            limits, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      loop = loops[0];
    } else if (sec->mtype == FACEBODY) {

      // get the face
      stat = EG_getTopology(sec, &eref, &oclass, &mtype,
                            limits, &nface, &faces, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      // get the loop
      stat = EG_getTopology(faces[0], &eref, &oclass, &mtype,
                            limits, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      loop = loops[0];
    }
  } else if (sec->oclass == FACE) {
    // get the loop
    stat = EG_getTopology(sec, &eref, &oclass, &mtype,
                          limits, &nloop, &loops, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    loop = loops[0];
  } else {
    loop = sec;
  }

  // get edges
  stat = EG_getTopology(loop, &eref, &oclass, &mtype,
                        limits, &nedge, &edges, &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (nedge != nstripe) {
    stat = EGADS_GEOMERR;
    goto cleanup;
  }

  uclosed = 1;
  if (mtype == OPEN) uclosed = 0;

  for (j = 0; j < nedge; j++) {
    stat = EG_getTopology(edges[j], &eref, &oclass, &mtype,
                          limits, &nnode, &node, &sens);
    nodes[j*inode] = NULL;
    if ((mtype == DEGENERATE) || (node == NULL)) continue;
    if (senses[j] == SFORWARD || mtype == ONENODE)
      nodes[j*inode] = node[0];
    else
      nodes[j*inode] = node[1];
  }

  nodes[nedge*inode] = NULL;
  if (node != NULL)
    if (uclosed == 0) {
      /* set the last node if open */
      if (senses[nedge-1] == SFORWARD || mtype == ONENODE)
        nodes[nedge*inode] = node[1];
      else
        nodes[nedge*inode] = node[0];
    } else {
      /* repeat the node for a uclosed loop */
      nodes[nedge*inode] = nodes[0];
    }

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf("EGADS Errror: EG_getSecNodes = %d!\n", stat);
  }
  return stat;
}


template<class T>
static int
EG_loft2spline(const int outLevel, const T *vknot, const int *vdata,
               const T **drnd, egSequ<T> *seq, int jmax, const T *xyz,
               double tol, int *header, T **data)
{
  int     status = EGADS_SUCCESS;
  int     i, j, imax, freeNT = 0, freeST = 0;
  T       dist, dy;
  T       x0[3], x1[3], nnor[3], snor[3], *norT, *souT;
  const T *north, *south;

  imax = seq->ncp;

  /* check for degenerate sides */
  north = south = NULL;
  norT  = souT  = NULL;
  i  = 0;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
    if (outLevel > 0)
      printf(" EGADS Error: Imin (west) degenerate (EG_ loft2spline)!\n");
    return EGADS_DEGEN;
  }
  j  = 0;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  loft2spline: Jmin (south) degenerate!\n");
#endif
    if (drnd[1] != NULL) {
      south = drnd[1];
      x0[0] = drnd[1][1];
      x0[1] = drnd[1][2];
      x0[2] = drnd[1][3];
      x1[0] = drnd[1][5];
      x1[1] = drnd[1][6];
      x1[2] = drnd[1][7];
      CROSS(snor, x0, x1);
      dist  = DOT(snor, snor);
      if ((dist == 0.0) || (DOT(x0, x1) > tol) ||
          (fabs(1.0-DOT(x0, x0)) > tol) ||
          (fabs(1.0-DOT(x1, x1)) > tol)) {
        if (outLevel > 0)
          printf(" EGADS Error: BAD South Axes (EG_loft2spline)!\n");
        return EGADS_NOTORTHO;
      }
      dist   = 1.0/sqrt(dist);
      snor[0] *= dist;
      snor[1] *= dist;
      snor[2] *= dist;
      souT     = snor;
#ifdef DEBUG
      printf("               with normal = %lf %lf %lf!\n",
             value(snor[0]), value(snor[1]), value(snor[2]));
#endif
    }
  } else {
    if (drnd[1] != NULL)
      if (drnd[1][0] != 0.0) {
        souT = (T*)EG_alloc(3*imax*sizeof(T));
        if (souT == NULL) { status = EGADS_MALLOC; goto cleanup; }
        freeST = 1;

        for (i = 0; i < imax; i++) {
          souT[3*i+0] = drnd[1][1]*drnd[1][0];
          souT[3*i+1] = drnd[1][2]*drnd[1][0];
          souT[3*i+2] = drnd[1][3]*drnd[1][0];
        }
      }
  }
  i  = imax-1;
  dy = 0.0;
  for (j = 1; j < jmax; j++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i)+(j-1)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i)+(j-1)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i)+(j-1)*imax)+2]));
  if (dy <= tol) {
    if (outLevel > 0)
      printf(" EGADS Error: Imax (east) degenerate (EG_ loft2spline)!\n");
    return EGADS_DEGEN;
  }
  j  = jmax-1;
  dy = 0.0;
  for (i = 1; i < imax; i++)
    dy += sqrt((xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ])*
               (xyz[3*((i)+(j)*imax)  ]-xyz[3*((i-1)+(j)*imax)  ]) +
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1])*
               (xyz[3*((i)+(j)*imax)+1]-xyz[3*((i-1)+(j)*imax)+1]) +
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2])*
               (xyz[3*((i)+(j)*imax)+2]-xyz[3*((i-1)+(j)*imax)+2]));
  if (dy <= tol) {
#ifdef DEBUG
    printf("  loft2spline: Jmax (north) degenerate!\n");
#endif
    if (drnd[3] != NULL) {
      north = drnd[3];
      x0[0] = drnd[3][1];
      x0[1] = drnd[3][2];
      x0[2] = drnd[3][3];
      x1[0] = drnd[3][5];
      x1[1] = drnd[3][6];
      x1[2] = drnd[3][7];
      CROSS(nnor, x0, x1);
      dist  = DOT(nnor, nnor);
      if ((dist == 0.0) || (DOT(x0, x1) > tol)||
          (fabs(1.0-DOT(x0, x0)) > tol) ||
          (fabs(1.0-DOT(x1, x1)) > tol)) {
        if (outLevel > 0)
          printf(" EGADS Error: BAD North Axes (EG_spline2d)!\n");
        return EGADS_NOTORTHO;
      }
      dist   = 1.0/sqrt(dist);
      nnor[0] *= dist;
      nnor[1] *= dist;
      nnor[2] *= dist;
      norT     = nnor;
#ifdef DEBUG
      printf("               with normal = %lf %lf %lf!\n",
             value(nnor[0]), value(nnor[1]), value(nnor[2]));
#endif
    }
  } else {
    if (drnd[3] != NULL)
      if (drnd[3][0] != 0.0) {
        norT = (T*)EG_alloc(3*imax*sizeof(T));
        if (norT == NULL) { status = EGADS_MALLOC; goto cleanup; }
        freeNT = 1;

        for (i = 0; i < imax; i++) {
          norT[3*i+0] = drnd[3][1]*drnd[3][0];
          norT[3*i+1] = drnd[3][2]*drnd[3][0];
          norT[3*i+2] = drnd[3][3]*drnd[3][0];
        }
      }
  }

  if (imax == 3) {
    status = EG_spline2dAppr<T>(1, imax, jmax, xyz, seq->knots, vknot, vdata,
                                NULL, NULL, south, souT, north, norT,
                                tol, header, data);
  } else {
    /* blend2xy (x and/or y is "b") fails with 2 */
    status = EG_spline2dAppr<T>(1, imax, jmax, xyz, seq->knots, vknot, vdata,
                                drnd[0], drnd[2], south, souT, north, norT,
                                tol, header, data);

  }

cleanup:

  if (freeST == 1) EG_free(souT);
  if (freeNT == 1) EG_free(norT);

  return status;
}


#ifdef EGADS_SPLINE_VELS
static int
EG_secSplinePointsVels(const egadsSplineVels *vels,
                       int outLevel, int lsec, int nsec, const ego *secs, int j,
                       egSequ<double> *ncp, int *planar, double *xyzs,
                       double* t1, double* tN)
{
  printf("EGADS Internal Error: EG_secSplinePointsVels called with double!\n");
  return EGADS_GEOMERR;
}


static int
EG_secSplinePointsVels(const egadsSplineVels *vels,
                       int outLevel, int lsec, int nsec, const ego *secs, int j,
                       egSequ< SurrealS<1> > *ncp, int *planar,
                       SurrealS<1> *xyzs, SurrealS<1> *t1, SurrealS<1> *tN)
{
  int         stat = EGADS_SUCCESS;
  int         i, jj, k, nn, npt, nchld, nedge, oclass, mtype;
  int         *senses, *iinfo=NULL;
  double      data[18], trangeD[2], trangeD_dot[2];
  double      xyz_dot[3], *ts=NULL, *ts_dot=NULL, *xs=NULL, *xs_dot=NULL;
  double      xyz[3], tbeg[3], tbeg_dot[3], tend[3], tend_dot[3];
  SurrealS<1> t=0, dt, point[6], v1[3], v2[3], trange[2];
  ego         loop=NULL, ref, *chldrn, *edges;


  npt = *planar = 0;
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
//    stat = EG_evaluate(secs[i], NULL, point);
      stat = (*(vels->velocityOfNode))(vels->usrData, secs, lsec+i, secs[i],
                                       NULL, xyz, xyz_dot);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Node %d eval = %d (EG_secSplinePointsVels)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      point[0].value() = xyz[0]; point[0].deriv() = xyz_dot[0];
      point[1].value() = xyz[1]; point[1].deriv() = xyz_dot[1];
      point[2].value() = xyz[2]; point[2].deriv() = xyz_dot[2];

      /* set the node sensitivity */
      stat = EG_setGeometry_dot(secs[i], NODE, 0, NULL, point);
      if (stat != EGADS_SUCCESS) goto cleanup;

      for (k = 0; k < ncp[j].ncp; k++, npt++) {
        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
      t1[3*i  ] = t1[3*i+1] = t1[3*i+2] = 0.0;
      tN[3*i  ] = tN[3*i+1] = tN[3*i+2] = 0.0;
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
      if (ref->mtype != PLANE) *planar = 1;
    } else if (oclass == BODY) {
      if (mtype == WIREBODY) {
        loop = chldrn[0];
      } else if (mtype == FACEBODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                              &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
        loop = chldrn[0];
      }
      if (EG_isPlanar(loop) != EGADS_SUCCESS) *planar = 1;
    } else {
      loop = secs[i];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) *planar = 1;
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nedge, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (nn = jj = 0; jj < nedge; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_secSplinePointsVels)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      if (mtype == DEGENERATE) continue;
      //stat = EG_evaluate(chldrn[0], NULL, v1);
      xyz_dot[0] = xyz_dot[1] = xyz_dot[2] = 0.;
      stat = (*(vels->velocityOfNode))(vels->usrData, secs, lsec+i, chldrn[0],
                                       edges[jj], xyz, xyz_dot);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d EDGE %d eval n1 = %d (EG_secSplinePointsVels)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      v1[0].value() = xyz[0]; v1[0].deriv() = xyz_dot[0];
      v1[1].value() = xyz[1]; v1[1].deriv() = xyz_dot[1];
      v1[2].value() = xyz[2]; v1[2].deriv() = xyz_dot[2];

      /* set the node sensitivity */
      stat = EG_setGeometry_dot(chldrn[0], NODE, 0, NULL, v1);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (nchld == 1) {
        v2[0] = v1[0];
        v2[1] = v1[1];
        v2[2] = v1[2];
      } else {
        //stat = EG_evaluate(chldrn[1], NULL, v2);
        stat = (*(vels->velocityOfNode))(vels->usrData, secs, lsec+i, chldrn[1],
                                         edges[jj], xyz, xyz_dot);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d EDGE %d getTOPOn2 = %d (EG_secSplinePointsVels)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }
        v2[0].value() = xyz[0]; v2[0].deriv() = xyz_dot[0];
        v2[1].value() = xyz[1]; v2[1].deriv() = xyz_dot[1];
        v2[2].value() = xyz[2]; v2[2].deriv() = xyz_dot[2];

        /* set the node sensitivity */
        stat = EG_setGeometry_dot(chldrn[1], NODE, 0, NULL, v2);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (nn == j) {
        if (vels->velocityOfRange != NULL) {
          stat = (*(vels->velocityOfRange))(vels->usrData, secs, lsec+i, edges[jj],
                                            trangeD, trangeD_dot);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d velocityOfRange = %d (EG_secKnotPointsVels)!\n",
                     i+1, jj+1, stat);
            goto cleanup;
          }
        } else if (EG_hasGeometry_dot(edges[jj]) == EGADS_SUCCESS) {
          stat = EG_getRange_dot(edges[jj], trangeD, trangeD_dot, &k);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_secKnotPointsVels)!\n",
                     i+1, jj+1, stat);
            goto cleanup;
          }
        } else {
          stat = EG_getRange(edges[jj], trangeD, &k);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_secKnotPointsVels)!\n",
                     i+1, jj+1, stat);
            goto cleanup;
          }
          trangeD_dot[0] = 0;
          trangeD_dot[1] = 0;
        }
        trange[0].value() = trangeD[0]; trange[0].deriv() = trangeD_dot[0];
        trange[1].value() = trangeD[1]; trange[1].deriv() = trangeD_dot[1];

        ts     = (double*)EG_alloc(  ncp[j].ncp*sizeof(double));
        ts_dot = (double*)EG_alloc(  ncp[j].ncp*sizeof(double));
        xs     = (double*)EG_alloc(3*ncp[j].ncp*sizeof(double));
        xs_dot = (double*)EG_alloc(3*ncp[j].ncp*sizeof(double));

        dt = trange[1] - trange[0];
        for (k = 0; k < ncp[j].ncp; k++) {
          if (senses[jj] == 1) {
            t = trange[0] + ncp[j].knots[k]*dt;
          } else {
            t = trange[1] - ncp[j].knots[k]*dt;
          }
          ts[k]     = t.value();
          ts_dot[k] = t.deriv();
        }

        stat = (*(vels->velocityOfEdge))(vels->usrData, secs, lsec+i, edges[jj],
                                         ncp[j].ncp,
                                         ts, ts_dot,
                                         xs, xs_dot,
                                         tbeg, tbeg_dot,
                                         tend, tend_dot);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d Edge %d velocityOfEdge = %d (EG_secSplinePointsVels)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }

        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          /* needed for EG_checkDirs */
          if (senses[jj] == 1) {
            t = trange[0] + ncp[j].knots[k]*dt;
          } else {
            t = trange[1] - ncp[j].knots[k]*dt;
          }

          point[0].value() = xs[3*k+0]; point[0].deriv() = xs_dot[3*k+0];
          point[1].value() = xs[3*k+1]; point[1].deriv() = xs_dot[3*k+1];
          point[2].value() = xs[3*k+2]; point[2].deriv() = xs_dot[3*k+2];

          if (k == 0) {
            point[3].value() = tbeg[0]; point[3].deriv() = tbeg_dot[0];
            point[4].value() = tbeg[1]; point[4].deriv() = tbeg_dot[1];
            point[5].value() = tbeg[2]; point[5].deriv() = tbeg_dot[2];
          } else if (k == ncp[j].ncp-1) {
            point[3].value() = tend[0]; point[3].deriv() = tend_dot[0];
            point[4].value() = tend[1]; point[4].deriv() = tend_dot[1];
            point[5].value() = tend[2]; point[5].deriv() = tend_dot[2];
          }

          xyzs[3*npt  ] = point[0];
          xyzs[3*npt+1] = point[1];
          xyzs[3*npt+2] = point[2];
          if (k == 0) {
            if (senses[jj] == 1) {
              t1[3*i  ]     =  point[3]*dt;
              t1[3*i+1]     =  point[4]*dt;
              t1[3*i+2]     =  point[5]*dt;
              xyzs[3*npt  ] =  v1[0];
              xyzs[3*npt+1] =  v1[1];
              xyzs[3*npt+2] =  v1[2];
              EG_checkDirs(edges[jj], t+0.001*dt, point, &t1[3*i]);
            } else {
              t1[3*i  ]     = -point[3]*dt;
              t1[3*i+1]     = -point[4]*dt;
              t1[3*i+2]     = -point[5]*dt;
              xyzs[3*npt  ] =  v2[0];
              xyzs[3*npt+1] =  v2[1];
              xyzs[3*npt+2] =  v2[2];
              EG_checkDirs(edges[jj], t-0.001*dt, point, &t1[3*i]);
            }
          } else if (k == ncp[j].ncp-1) {
            if (senses[jj] == 1) {
              tN[3*i  ]     = -point[3]*dt;
              tN[3*i+1]     = -point[4]*dt;
              tN[3*i+2]     = -point[5]*dt;
              xyzs[3*npt  ] =  v2[0];
              xyzs[3*npt+1] =  v2[1];
              xyzs[3*npt+2] =  v2[2];
              EG_checkDirs(edges[jj], t-0.001*dt, point, &tN[3*i]);
            } else {
              tN[3*i  ]     =  point[3]*dt;
              tN[3*i+1]     =  point[4]*dt;
              tN[3*i+2]     =  point[5]*dt;
              xyzs[3*npt  ] =  v1[0];
              xyzs[3*npt+1] =  v1[1];
              xyzs[3*npt+2] =  v1[2];
              EG_checkDirs(edges[jj], t+0.001*dt, point, &tN[3*i]);
            }
          }
        }
        EG_free(ts);     ts     = NULL;
        EG_free(ts_dot); ts_dot = NULL;
        EG_free(xs);     xs     = NULL;
        EG_free(xs_dot); xs_dot = NULL;
        break;
      }
      nn++;
    }
  }

  stat = EGADS_SUCCESS;

cleanup:
  EG_free(ts);     ts     = NULL;
  EG_free(ts_dot); ts_dot = NULL;
  EG_free(xs);     xs     = NULL;
  EG_free(xs_dot); xs_dot = NULL;

  return stat;
}
#endif


template<class T>
static int
EG_secSplinePoints(
#ifdef EGADS_SPLINE_VELS
    const egadsSplineVels *vels,
#endif
    int outLevel, int lsec, int nsec, const ego *secs, int j,
                   egSequ<T> *ncp, int *planar, T *xyzs, T *t1, T *tN)
{
  int    stat = EGADS_SUCCESS;
  int    i, jj, k, nn, npt, nchld, nedge, oclass, mtype;
  int    *senses, *iinfo=NULL;
  double data[18];
  T      t, point[18], v1[3], v2[3], trange[2], dt;
  ego    loop=NULL, ref, *chldrn, *edges;

#ifdef EGADS_SPLINE_VELS
  if (vels != NULL) {
    return EG_secSplinePointsVels(vels, outLevel, lsec, nsec, secs, j, ncp,
                                  planar, xyzs, t1, tN);
  }
#endif

  npt = *planar = 0;
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      stat = EG_evaluate(secs[i], NULL, point);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Node %d eval = %d (EG_secSplinePoints)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      for (k = 0; k < ncp[j].ncp; k++, npt++) {
        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
      t1[3*i  ] = t1[3*i+1] = t1[3*i+2] = 0.0;
      tN[3*i  ] = tN[3*i+1] = tN[3*i+2] = 0.0;
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
      if (ref->mtype != PLANE) *planar = 1;
    } else if (oclass == BODY) {
      if (mtype == WIREBODY) {
        loop = chldrn[0];
      } else if (mtype == FACEBODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                              &senses);
        loop = chldrn[0];
      }
      if (EG_isPlanar(loop) != EGADS_SUCCESS) *planar = 1;
    } else {
      loop = secs[i];
      if (EG_isPlanar(loop) != EGADS_SUCCESS) *planar = 1;
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nedge, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (nn = jj = 0; jj < nedge; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d getTOPO = %d (EG_secSplinePoints)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      if (mtype == DEGENERATE) continue;
      stat = EG_evaluate(chldrn[0], NULL, v1);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Sec %d EDGE %d eval n1 = %d (EG_secSplinePoints)!\n",
                 i+1, j+1, stat);
        goto cleanup;
      }
      if (nchld == 1) {
        v2[0] = v1[0];
        v2[1] = v1[1];
        v2[2] = v1[2];
      } else {
        stat = EG_evaluate(chldrn[1], NULL, v2);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d EDGE %d getTOPOn2 = %d (EG_secSplinePoints)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }
      }
      if (nn == j) {
        stat = EG_getRange(edges[jj], trange, &k);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d getRange = %d (EG_secSplinePoints)!\n",
                   i+1, j+1, stat);
          goto cleanup;
        }
        dt = trange[1] - trange[0];
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          if (senses[jj] == 1) {
            t = trange[0] + ncp[j].knots[k]*dt;
          } else {
            t = trange[1] - ncp[j].knots[k]*dt;
          }
          stat = EG_evaluate(edges[jj], &t, point);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Section %d Edge %d eval = %d (EG_secSplinePoints)!\n",
                     i+1, j+1, stat);
            goto cleanup;
          }
          xyzs[3*npt  ] = point[0];
          xyzs[3*npt+1] = point[1];
          xyzs[3*npt+2] = point[2];

          if (k == 0) {
            if (senses[jj] == 1) {
              t1[3*i  ]     =  point[3]*dt;
              t1[3*i+1]     =  point[4]*dt;
              t1[3*i+2]     =  point[5]*dt;
              xyzs[3*npt  ] =  v1[0];
              xyzs[3*npt+1] =  v1[1];
              xyzs[3*npt+2] =  v1[2];
              EG_checkDirs(edges[jj], t+0.001*dt, point, &t1[3*i]);
            } else {
              t1[3*i  ]     = -point[3]*dt;
              t1[3*i+1]     = -point[4]*dt;
              t1[3*i+2]     = -point[5]*dt;
              xyzs[3*npt  ] =  v2[0];
              xyzs[3*npt+1] =  v2[1];
              xyzs[3*npt+2] =  v2[2];
              EG_checkDirs(edges[jj], t-0.001*dt, point, &t1[3*i]);
            }
          } else if (k == ncp[j].ncp-1) {
            if (senses[jj] == 1) {
              tN[3*i  ]     = -point[3]*dt;
              tN[3*i+1]     = -point[4]*dt;
              tN[3*i+2]     = -point[5]*dt;
              xyzs[3*npt  ] =  v2[0];
              xyzs[3*npt+1] =  v2[1];
              xyzs[3*npt+2] =  v2[2];
              EG_checkDirs(edges[jj], t-0.001*dt, point, &tN[3*i]);
            } else {
              tN[3*i  ]     =  point[3]*dt;
              tN[3*i+1]     =  point[4]*dt;
              tN[3*i+2]     =  point[5]*dt;
              xyzs[3*npt  ] =  v1[0];
              xyzs[3*npt+1] =  v1[1];
              xyzs[3*npt+2] =  v1[2];
              EG_checkDirs(edges[jj], t+0.001*dt, point, &tN[3*i]);
            }
          }
        }
        break;
      }
      nn++;
    }
  }

  stat = EGADS_SUCCESS;

cleanup:
  return stat;
}


#ifdef EGADS_SPLINE_VELS
static int
EG_secKnotPointsVels(const egadsSplineVels *vels,
                     int outLevel, int nsec, const ego *secs,
                     int nstripe, int uclosed, double *xyzs)
{
  printf("EGADS Internal Error: EG_secKnotPointsVels called with double!\n");
  return EGADS_GEOMERR;
}


static int
EG_secKnotPointsVels(const egadsSplineVels *vels,
                     int outLevel, int nsec, const ego *secs,
                     int nstripe, int uclosed, SurrealS<1> *xyzs)
{
  int         stat = EGADS_SUCCESS;
  int         i, k, n, jj, n0, n1, nchld, npt, oclass, mtype;
  int         *senses, *iinfo=NULL;
  double      data[18], trangeD[2], trangeD_dot[2], xyz[3], xyz_dot[3];
  double      ts[4], ts_dot[4], xs[3*4], xs_dot[3*4];
  double      tbeg[3], tbeg_dot[3], tend[3], tend_dot[3];
  SurrealS<1> point[18], t, trange[2];
  ego         ref, loop=NULL, *edges, *chldrn;

  for (npt = i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
//    stat = EG_evaluate(secs[i], NULL, point);
      stat = (*(vels->velocityOfNode))(vels->usrData, secs, i, secs[i], NULL, xyz,
                                       xyz_dot);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Node eval = %d (EG_secKnotPointsVels)!\n",
                 i+1, stat);
        goto cleanup;
      }
      point[0].value() = xyz[0]; point[0].deriv() = xyz_dot[0];
      point[1].value() = xyz[1]; point[1].deriv() = xyz_dot[1];
      point[2].value() = xyz[2]; point[2].deriv() = xyz_dot[2];

      for (k = 0; k < 3*nstripe+1-uclosed; k++, npt++) {
        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      if (mtype == WIREBODY) {
        loop = chldrn[0];
      } else if (mtype == FACEBODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                              &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
        loop = chldrn[0];
      }
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (n1 = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d GetTOPO = %d (EG_secKnotPointsVels)!\n",
                 i+1, jj+1, stat);
        goto cleanup;
      }
      if (mtype == DEGENERATE) continue;
      if (vels->velocityOfRange != NULL) {
        stat = (*(vels->velocityOfRange))(vels->usrData, secs, i, edges[jj],
                                          trangeD, trangeD_dot);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d velocityOfRange = %d (EG_secKnotPointsVels)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
      } else if (EG_hasGeometry_dot(edges[jj]) == EGADS_SUCCESS) {
        stat = EG_getRange_dot(edges[jj], trangeD, trangeD_dot, &k);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_secKnotPointsVels)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
      } else {
        stat = EG_getRange(edges[jj], trangeD, &k);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Sec %d Edge %d GetRange = %d (EG_secKnotPointsVels)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
        trangeD_dot[0] = 0;
        trangeD_dot[1] = 0;
      }
      trange[0].value() = trangeD[0]; trange[0].deriv() = trangeD_dot[0];
      trange[1].value() = trangeD[1]; trange[1].deriv() = trangeD_dot[1];

      n0 = 1;
      if ((uclosed == 0) && (n1 == 0)) {
        n1++;
        n0 = 0;
      }

      for (k = 0; k < 4; k++) {
        if (senses[jj] == 1) {
          t = trange[0] + k*(trange[1] - trange[0])/3.0;
        } else {
          t = trange[1] - k*(trange[1] - trange[0])/3.0;
        }
        ts[k]     = t.value();
        ts_dot[k] = t.deriv();
      }

      stat = (*(vels->velocityOfEdge))(vels->usrData, secs, i, edges[jj],
                                       4,
                                       ts, ts_dot,
                                       xs, xs_dot,
                                       tbeg, tbeg_dot,
                                       tend, tend_dot);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d velocityOfEdge = %d (EG_secKnotPointsVels)!\n",
                 i+1, jj+1, stat);
        goto cleanup;
      }

      for (k = n0; k < 4; k++, npt++) {
        point[0].value() = xs[3*k+0]; point[0].deriv() = xs_dot[3*k+0];
        point[1].value() = xs[3*k+1]; point[1].deriv() = xs_dot[3*k+1];
        point[2].value() = xs[3*k+2]; point[2].deriv() = xs_dot[3*k+2];

        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
    }
  }

cleanup:
  return stat;
}
#endif


template<class T>
static int
EG_secKnotPoints(
#ifdef EGADS_SPLINE_VELS
                 const egadsSplineVels *vels,
#endif
                 int outLevel, int nsec, const ego *secs,
                 int nstripe, int uclosed, T *xyzs)
{
  int    stat = EGADS_SUCCESS;
  int    i, k, n, jj, n0, n1, nchld, npt, oclass, mtype;
  int    *senses, *iinfo=NULL;
  double data[18];
  T      point[18], trange[2], t;
  ego    ref, loop=NULL, *edges, *chldrn;

#ifdef EGADS_SPLINE_VELS
  if (vels != NULL) {
    return EG_secKnotPointsVels(vels, outLevel, nsec, secs, nstripe, uclosed, xyzs);
  }
#endif

  for (npt = i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    if (oclass == NODE) {
      stat = EG_evaluate(secs[i], NULL, point);
      for (k = 0; k < 3*nstripe+1-uclosed; k++, npt++) {
        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
      continue;
    } else if (oclass == FACE) {
      loop = chldrn[0];
    } else if (oclass == BODY) {
      if (mtype == WIREBODY) {
        loop = chldrn[0];
      } else if (mtype == FACEBODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                              &senses);
        loop = chldrn[0];
      }
    } else {
      loop = secs[i];
    }
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) continue;
    for (n1 = jj = 0; jj < n; jj++) {
      stat = EG_getTopology(edges[jj], &ref, &oclass, &mtype, data, &nchld,
                            &chldrn, &iinfo);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d EDGE %d GetTOPO = %d (EG_secKnotPoints)!\n",
                 i+1, jj+1, stat);
        goto cleanup;
      }
      if (mtype == DEGENERATE) continue;
      stat = EG_getRange(edges[jj], trange, &k);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d GetRange = %d (EG_secKnotPoints)!\n",
                 i+1, jj+1, stat);
        goto cleanup;
      }
      n0 = 1;
      if ((uclosed == 0) && (n1 == 0)) {
        n1++;
        n0 = 0;
      }
      for (k = n0; k < 4; k++, npt++) {
        if (senses[jj] == 1) {
          t = trange[0] + k*(trange[1] - trange[0])/3.0;
        } else {
          t = trange[1] - k*(trange[1] - trange[0])/3.0;
        }
        stat = EG_evaluate(edges[jj], &t, point);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d Edge %d eval = %d (EG_secKnotPoints)!\n",
                   i+1, jj+1, stat);
          goto cleanup;
        }
        xyzs[3*npt  ] = point[0];
        xyzs[3*npt+1] = point[1];
        xyzs[3*npt+2] = point[2];
      }
    }
  }

cleanup:
  return stat;
}

static int
EG_getSecEdge(int outLevel, const ego sec, int stripe, ego *secEdge, int *sens)
{
  int    stat = EGADS_SUCCESS;
  int    nchld, oclass, mtype;
  int    *senses;
  double data[18];
  ego    ref, loop=NULL, *edges, *chldrn;

  stat = EG_getTopology(sec, &ref, &oclass, &mtype, data, &nchld, &chldrn,
                        &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (oclass == NODE) {
    *secEdge = sec;
    *sens = 0;
    return EGADS_SUCCESS;
  } else if (oclass == FACE) {
    loop = chldrn[0];
  } else if (oclass == BODY) {
    if (mtype == WIREBODY) {
      loop = chldrn[0];
    } else if (mtype == FACEBODY) {
      stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nchld, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      loop = chldrn[0];
    }
  } else {
    loop = sec;
  }
  stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nchld, &edges,
                        &senses);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: GetTOPO = %d (EG_getSecEdge)!\n",
             stat);
    goto cleanup;
  }

  *secEdge = edges[stripe];
  *sens = senses[stripe];

cleanup:
  return stat;
}

static int
EG_secEquivEdge(int outLevel, const ego *secs, int stripe)
{
  int    stat = EGADS_SUCCESS;
  int    i, sens;
  ego    secEdges[2];

  for (i = 0; i < 2; i++) {
    stat = EG_getSecEdge(outLevel, secs[i], stripe, &secEdges[i], &sens);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (secEdges[0]->oclass != secEdges[1]->oclass) return EGADS_OUTSIDE;

  /* check if the section Edges are equivalent */
  stat = EG_isEquivalent(secEdges[0], secEdges[1]);

cleanup:
  return stat;
}

#ifdef EGADS_SPLINE_VELS
static int
EG_nodeVels(const egadsSplineVels *vels,
            int outLevel, int lsec, const ego *secs, const ego node,
            double *xyzs)
{
  return EG_evaluate(node, NULL, xyzs);
}

static int
EG_nodeVels(const egadsSplineVels *vels,
            int outLevel, int lsec, const ego *secs, const ego node,
            SurrealS<1> *xyzs)
{
  int         stat = EGADS_SUCCESS;
  double      xyz[3], xyz_dot[3];
  SurrealS<1> point[3];

  if (vels == NULL)
    return EG_evaluate(node, NULL, xyzs);

  stat = (*(vels->velocityOfNode))(vels->usrData, secs, lsec, node,
      NULL, xyz, xyz_dot);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Section %d Node eval = %d (EG_nodeVels)!\n",
             lsec+1, stat);
    goto cleanup;
  }
  point[0].value() = xyz[0]; point[0].deriv() = xyz_dot[0];
  point[1].value() = xyz[1]; point[1].deriv() = xyz_dot[1];
  point[2].value() = xyz[2]; point[2].deriv() = xyz_dot[2];

  /* set the node sensitivity */
  stat = EG_setGeometry_dot(node, NODE, 0, NULL, point);
  if (stat != EGADS_SUCCESS) goto cleanup;

  xyzs[0] = point[0];
  xyzs[1] = point[1];
  xyzs[2] = point[2];

cleanup:
  return stat;
}
#endif


template<class T>
static int
EG_blendSpline(
#ifdef EGADS_SPLINE_VELS
               const egadsSplineVels *vels,
#endif
               int nsex, const ego *secs, T *rc1, T *rcN, int begRC, int endRC,
               int nstripe, int uclosed, int vclosed, int planar,
               int *nsecC0_out, ego **secsC0_out,
               ego **nodes_out, egSpline<T> **curvsU_out,
               egSpline<T> **curvsV_out, egSpline<T> **surfs_out,
               int &tip, int &te, egSpline<T> *tipsurfs, egSpline<T> *tipcurvs)
{
  int    stat = EGADS_SUCCESS;
  int    i, j, k, n, ii, jj, kk, outLevel, oclass, mtype, npt, isrf;
  int    inode, fixed, icrvU, icrvV, nsecC0=0;
  int    nnode=0, ncurvU=0, ncurvV, nsurf=0, n0, *senses=NULL;
  int    nsec, ntip, rite, left, nKnotv, deg, iknotc, *vdata=NULL;
  int    *vdataTE=NULL, *aggC0=NULL;
  double data[18];
  T      *t1, *tN, *vknot=NULL, *vknotC0=NULL, *vknotTE=NULL, max2, mknot;
  T      knotc[2*NTIPTE-1], *snor=NULL, *nnor=NULL, *xyzs=NULL;
  T      dx, dy, dz, dv, dv0, dv1, t, pointL[18], pointU[18], dX[3];
  T      d2xdt2L, d2ydt2L, d2zdt2L, d2sdt2L, d2xdt2R, d2ydt2R, d2zdt2R, d2sdt2R;
  ego    loop=NULL, ref;
  ego    *chldrn=NULL, *edges=NULL, *nodes=NULL, *secsC0=NULL;
  egSequ<T> *ncp=NULL;
  egSpline<T> *surfs=NULL, *curvsV=NULL, *curvsU=NULL, *neigbr[2]={NULL,NULL};
#ifdef BLEND_SPLIT_CONSTRUCTION
  int nC0, iC0;
#else
  egSpline<T> *blendsurfs=NULL;
#endif
  const T *sides[4]={NULL,NULL,NULL,NULL};

  *nsecC0_out = 0;
  *secsC0_out = NULL;
  *nodes_out  = NULL;
  *curvsU_out = NULL;
  *curvsV_out = NULL;
  *surfs_out  = NULL;

  outLevel = EG_outLevel(secs[0]);
  fixed    = EG_fixedKnots(secs[0]);

  nsec = nsex;
  if (nsec < 0) nsec = -nsec;

  /* look for tip treatment */
  ntip = 0;
  tip = 0;
  if ((secs[     0]->oclass == FACE || secs[     0]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    tip += 1;
    ntip++;
  }
  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    tip += 2;
    ntip++;
  }

  if ((ntip > 0) && !(nstripe == 3 || nstripe == 2) ) {
    if (outLevel > 0)
      printf(" EGADS Error: Section has %d Edges (may only have 2 or 3) (EG_blend)!\n",
             nstripe);
    return EGADS_GEOMERR;
  }

  /* set the knots in the loft direction */
  vknot   = (T *)   EG_alloc(nsec*sizeof(T));
  vdata   = (int *) EG_alloc(nsec*sizeof(int));
  vknotC0 = (T *)   EG_alloc(nsec*sizeof(T));
  aggC0   = (int *) EG_alloc(2*nsec*sizeof(int));
  xyzs    = (T *)   EG_alloc(3*(3*nstripe+1-uclosed)*nsec*sizeof(T));
  secsC0  = (ego *) EG_alloc(nsec*sizeof(ego));
  if ((vknot == NULL)   || (xyzs == NULL)  || (vdata == NULL)   ||
      (vknotC0 == NULL) || (aggC0 == NULL) || (secsC0 == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for (EG_blend)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i <   nsec; i++) secsC0[i] = NULL;
  for (i = 0; i < 2*nsec; i++) aggC0[i]  = 0;
  
#ifdef EGADS_SPLINE_VELS
  /* Perform one redundant call to make sure the section is populated with dot information (if available) */
  if (vels != NULL) {
    stat = EG_secKnotPoints(vels, outLevel, nsec, secs, nstripe, uclosed, xyzs);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }
#endif
#ifdef EGADS_SPLINE_VELS
  stat = EG_secKnotPoints(vels, outLevel, nsec, secs, nstripe, uclosed, xyzs);
#else
  stat = EG_secKnotPoints(outLevel, nsec, secs, nstripe, uclosed, xyzs);
#endif
  if (stat != EGADS_SUCCESS) goto cleanup;

  n0 = 3*nstripe+1-uclosed;
  /* arc-length spaced */
  for (i = 0; i < nsec; i++) vknot[i] = 0.0;
  dz = n0;
  dv = 0.0;
  for (i = 0; i < n0; i++) {
    dy = 0.0;
    for (j = 1; j < nsec; j++) {
      dy += sqrt((xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ])*
                 (xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ]) +
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1])*
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1]) +
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2])*
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2]));
    }
    if (dy == 0.0) {
      dz -= 1.0;
      continue;
    }
    dv += dy;
    dx = 0.0;
    for (j = 1; j < nsec; j++) {
      dx += sqrt((xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ])*
                 (xyzs[3*(i+j*n0)  ]-xyzs[3*(i+(j-1)*n0)  ]) +
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1])*
                 (xyzs[3*(i+j*n0)+1]-xyzs[3*(i+(j-1)*n0)+1]) +
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2])*
                 (xyzs[3*(i+j*n0)+2]-xyzs[3*(i+(j-1)*n0)+2]))/dy;
      vknot[j] += dx;
    }
  }
  for (n = j = 1; j < nsec; j++)
    if (vknot[j] < vknot[j-1]) n++;
  if ((n == 1) && (fixed == 0)) {
    for (j = 0; j < nsec; j++) vknot[j] /= dz;
    dv /= dz;
  } else {
    /* equally spaced */
    for (jj = j = 0; j < nsec-1; j++)
      if (secs[j] != secs[j+1]) jj++;
    for (k = j = 0; j < nsec-1; j++) {
      dy       = k;
      vknot[j] = dy/jj;
      if (secs[j] != secs[j+1]) k++;
    }
    vknot[nsec-1] = 1.;
    dv = 1.0;
  }

  /* always starts with the first section */
  nsecC0    = 1;
  secsC0[0] = secs[0];
#ifdef BLEND_SPLIT_CONSTRUCTION
  aggC0[0]  = 0;
#else
  aggC0[0]  = 1;
#endif
  vdata[0] = vdata[nsec-1] = 1;
  for (j = 1; j < nsec-1; j++) {
    if ((j < nsec-3) && (vknot[j] == vknot[j+1]) && (vknot[j] == vknot[j+2])) {
      if ((j < nsec-4) && (vknot[j] == vknot[j+3])) {
        printf("EG_blend: Multiplicity > 3 at %d\n", j);
        goto cleanup;
      }
#ifdef DEBUG
      printf("repeated vknot at j=%d, %d, and %d\n", j, j+1, j+2);
#endif
      if (nsex > 0) {
#ifdef BLEND_SPLIT_CONSTRUCTION
        aggC0[2*nsecC0-1] = j+1;
        aggC0[2*nsecC0  ] = j+2;
#else
        aggC0[2*nsecC0-1] = j+2;
        aggC0[2*nsecC0  ] = j+3;
#endif
        secsC0[nsecC0++ ] = secs[j];
        vdata[j++] = 1;
        vdata[j++] = 1;
        vdata[j  ] = 1;
      } else {
        vdata[j++] = -3;
        vdata[j++] =  1;
        vdata[j  ] = +3;
      }
    } else if (vknot[j] == vknot[j+1]) {
#ifdef DEBUG
      printf("repeated vknot at j=%d and %d\n", j, j+1);
#endif
      vdata[j++] = -2;
      vdata[j  ] = +2;
    } else {
      vdata[j  ] =  1;
    }
#ifdef DEBUG
  printf("after pass 1\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d   %lf\n", j, vdata[j], value(vknot[j]));
#endif
  }

  /* and cap of the last section */
#ifdef BLEND_SPLIT_CONSTRUCTION
  aggC0[2*nsecC0-1] = nsec;
#else
  aggC0[2*nsecC0-1] = nsec+1;
#endif
  secsC0[nsecC0++ ] = secs[nsec-1];

  /* for all multiplicity 2 knots, determine where the flat spot is */
  for (j = 1; j < nsec-1; j++)
    if ((vdata[j] == -2) && (vdata[j+1] == +2)) {
      if (j < 2) {
        /* flat on left */
        vdata[j+1] = 1;
      } else if (j > nsec-4) {
        /* flat on right */
        vdata[j  ] = 1;
      } else if ((vknot[j-1] == vknot[j-2]) && (vknot[j+3] == vknot[j+2])) {
        /* flat on both sides */
        printf("EG_blend: C1 -- Too few sections on either side of %d and %d\n",
               j, j+1);
        stat = EGADS_DEGEN;
        goto cleanup;
      } else if (vknot[j-1] == vknot[j-2]) {
        /* flat on left because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on left at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j+1] = 1;
      } else if (vknot[j+3] == vknot[j+2]) {
        /* flat on right because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on right at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j  ] = 1;
      } else {
        /* find second derivatives on both sides */
        for (left = rite = i = 0; i < n0; i++) {
          d2xdt2L = ((xyzs[3*(i+(j  )*n0)  ]-
                      xyzs[3*(i+(j-1)*n0)  ])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)  ]-
                      xyzs[3*(i+(j-2)*n0)  ])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);
          d2ydt2L = ((xyzs[3*(i+(j  )*n0)+1]-
                      xyzs[3*(i+(j-1)*n0)+1])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)+1]-
                      xyzs[3*(i+(j-2)*n0)+1])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);
          d2zdt2L = ((xyzs[3*(i+(j  )*n0)+2]-
                      xyzs[3*(i+(j-1)*n0)+2])/(vknot[j  ]-vknot[j-1])
                  -  (xyzs[3*(i+(j-1)*n0)+2]-
                      xyzs[3*(i+(j-2)*n0)+2])/(vknot[j-1]-vknot[j-2]))
                  / (vknot[j  ]-vknot[j-2]);

          d2xdt2R = ((xyzs[3*(i+(j+3)*n0)  ]-
                      xyzs[3*(i+(j+2)*n0)  ])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)  ]-
                      xyzs[3*(i+(j+1)*n0)  ])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);
          d2ydt2R = ((xyzs[3*(i+(j+3)*n0)+1]-
                      xyzs[3*(i+(j+2)*n0)+1])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)+1]-
                      xyzs[3*(i+(j+1)*n0)+1])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);
          d2zdt2R = ((xyzs[3*(i+(j+3)*n0)+2]-
                      xyzs[3*(i+(j+2)*n0)+2])/(vknot[j+3]-vknot[j+2])
                  -  (xyzs[3*(i+(j+2)*n0)+2]-
                      xyzs[3*(i+(j+1)*n0)+2])/(vknot[j+2]-vknot[j+1]))
                  / (vknot[j+3]-vknot[j+1]);

          d2sdt2L = d2xdt2L*d2xdt2L + d2ydt2L*d2ydt2L + d2zdt2L*d2zdt2L;
          d2sdt2R = d2xdt2R*d2xdt2R + d2ydt2R*d2ydt2R + d2zdt2R*d2zdt2R;
          max2    = fabs(d2sdt2L);
          if (max2 < fabs(d2sdt2R)) max2 = fabs(d2sdt2R);
          if (max2 == 0.0) {
            left++;
          } else if (fabs(d2sdt2L-d2sdt2R)/max2 < 1.e-6) {
            left++;
          } else if (d2sdt2L > d2sdt2R) {
            rite++;
          } else {
            left++;
          }
        }
        if (rite > left) {
          /* flat on right because left curvature is smaller */
#ifdef DEBUG
          printf("flat on right %d %d\n", left, rite);
#endif
          vdata[j  ] = 1;
        } else {
          /* flat on left because right curvture is smaller  */
#ifdef DEBUG
          printf("flat on left %d %d\n", left, rite);
#endif
          vdata[j+1] = 1;
        }
      }
      j++;
    }
#ifdef DEBUG
  printf("after pass 2\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d\n", j, vdata[j]);
#endif
  EG_free(xyzs); xyzs = NULL;

  if (ntip > 0) {
    if (nstripe == 3) {

      for (npt = i = 0; i < nsec; i++) {
        stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                              &senses);
        if (stat != EGADS_SUCCESS) continue;
        if (oclass == NODE) {
          continue;
        } else if (oclass == FACE) {
          loop = chldrn[0];
        } else if (oclass == BODY) {
          if (mtype == WIREBODY) {
            loop = chldrn[0];
          } else if (mtype == FACEBODY) {
            stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &n, &chldrn,
                                  &senses);
            loop = chldrn[0];
          }
        } else {
          loop = secs[i];
        }
        stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &n, &edges,
                              &senses);
        if (stat != EGADS_SUCCESS) continue;

        stat = EG_findTE(nstripe, edges, &te);
        if (stat != EGADS_SUCCESS) goto cleanup;
        break;
      }
    } else {
      te = 2;
    }
  } else {
    te = -1;
  }

  /* get the sampling points for each stripe */
#ifdef EGADS_SPLINE_VELS
  stat = EG_setSequence(vels, nsec, secs, te, nstripe, &ncp);
#else
  stat = EG_setSequence(nsec, secs, te, nstripe, &ncp);
#endif
  if (stat != EGADS_SUCCESS) goto cleanup;

  inode = nsecC0;
  icrvU = nsecC0-1;
  icrvV = nsecC0;
  isrf = nsecC0-1;

  nnode  = (nstripe+1)*nsecC0;
  ncurvU = (nstripe+1)*icrvU;
  ncurvV = nstripe*icrvV;
  nsurf  = nstripe*isrf;

  /* make a surface for each patch/stripe */
  nodes  = (ego*)EG_alloc( nnode *sizeof(ego));
  curvsU = new egSpline<T>[ncurvU];
  curvsV = new egSpline<T>[ncurvV];
  surfs  = new egSpline<T>[nsurf];
  if ((nodes == NULL) || (curvsU == NULL) || (curvsV == NULL) ||
      (surfs == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation (EG_blend)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }

#ifdef BLEND_SPLIT_CONSTRUCTION
  for (iC0 = 0; iC0 < nsecC0-1; iC0++) {
    nC0 = aggC0[2*iC0+1]-aggC0[2*iC0];

    /* normalize the knots based on the C0 sections */
    dy = vknot[aggC0[2*iC0+1]-1]-vknot[aggC0[2*iC0]];
    for (i = 0; i < nC0; i++) {
      vknotC0[i] = (vknot[i+aggC0[2*iC0]] - vknot[aggC0[2*iC0]])/dy;
    }

    for (j = 0; j < nstripe; j++) {
      xyzs = (T *) EG_alloc(3*(ncp[j].ncp+2)*nC0*sizeof(T));
      if (xyzs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Allocation for %dX%d points (EG_blend)!\n",
                 ncp[j].ncp, nC0);
        stat = EGADS_MALLOC;
        goto cleanup;
      }
      t1 = &xyzs[3* ncp[j].ncp   *nC0];
      tN = &xyzs[3*(ncp[j].ncp+1)*nC0];

      /* get the spline points between the two C0 sections for the current stripe */
      stat = EG_secSplinePoints(outLevel, nC0, &secs[aggC0[2*iC0]], j, ncp,
                                &planar, xyzs, t1, tN);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* index now based on the C0 sections */
      i = iC0;

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSplinePoints_%d.dat", j);
        EG_tecplotSplinePoints(ncp[j].ncp, nC0, xyzs, filename);
      }
#endif

      if (nC0 == 2 && iC0 > 0 && iC0 < nsecC0-2 ) {
        /* single ruled spline */
        if ((planar == 1) || (ncp[j].ncp == 3)) {
          stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                    1.e-8, surfs[i+j*isrf].header,
                                    &surfs[i+j*isrf].data);
        } else {
          stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                    NULL, t1, tN, NULL, NULL, NULL, NULL,
                                    1.e-8, surfs[i+j*isrf].header,
                                    &surfs[i+j*isrf].data);
        }
      } else {

        if (planar == 0) {
          sides[0] = t1;
          sides[2] = tN;
        } else {
          sides[0] = NULL;
          sides[2] = NULL;
        }

        sides[1] = iC0 == 0 ? rc1 : NULL;
        sides[3] = iC0 == nsecC0-2 ? rcN : NULL;

        /* get the BSpline surface */
        stat = EG_loft2spline(outLevel, vknotC0, &vdata[aggC0[2*iC0]],
                              sides, &ncp[j], nC0, xyzs, 1.e-8,
                              surfs[i+j*isrf].header, &surfs[i+j*isrf].data);
      }
      EG_free(xyzs); xyzs = NULL;
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Strip %d splined = %d (EG_blend)!\n", j+1, stat);
        goto cleanup;
      }

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSpline_%dx%d.dat", i,j);
        EG_tecplotSpline(surfs[i+j*isrf].header, surfs[i+j*isrf].data, filename);
      }
#endif

      /* get the u = 0 curveU[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                         0, -1,
                         curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      if (j == nstripe-1) {
        if (uclosed == 1) {
          /* uclosed loop so first and last curves are the same */
        } else {
          /* get the u = 1 curveU[i,j+1] */
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             surfs[i+j*isrf].header[2]-1, -1,
                             curvsU[i+(j+1)*icrvU].header,
                             &curvsU[i+(j+1)*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }

      /* get the v = 0 curveV[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                         -1, 0,
                         curvsV[i+j*icrvV].header, &curvsV[i+j*icrvV].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      /* get the v = 1 curveV[i+1,j] on the last section */
      if ( i == nsecC0-2 ) {
        stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                           -1, surfs[i+j*isrf].header[5]-1,
                           curvsV[(i+1)+j*icrvV].header,
                           &curvsV[(i+1)+j*icrvV].data);

        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d isoCurve+ = %d (EG_blend)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
      }
    }
  }

#else

  blendsurfs = new egSpline<T>[nstripe];
  if (blendsurfs == NULL) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Surfaces (EG_blend)!\n", nstripe);
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  for (j = 0; j < nstripe; j++) {
    if (j == te) continue;
    xyzs = (T *) EG_alloc(3*(ncp[j].ncp+2)*nsec*sizeof(T));
    if (xyzs == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation for %dX%d points (EG_blend)!\n",
               ncp[j].ncp, nsec);
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    t1 = &xyzs[3* ncp[j].ncp   *nsec];
    tN = &xyzs[3*(ncp[j].ncp+1)*nsec];

    if (planar == 0) {
      sides[0] = t1;
      sides[2] = tN;
    } else {
      sides[0] = NULL;
      sides[2] = NULL;
    }

    sides[1] = begRC > 0 ? rc1 : NULL;
    sides[3] = endRC > 0 ? rcN : NULL;

    /* get the spline points across all sections for the current stripe */
#ifdef EGADS_SPLINE_VELS
    stat = EG_secSplinePoints(vels, outLevel, 0, nsec, secs, j, ncp, &n, xyzs, t1, tN);
#else
    stat = EG_secSplinePoints(outLevel, 0, nsec, secs, j, ncp, &n, xyzs, t1, tN);
#endif
    if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef SPLINE_TECPLOT_DEBUG
    {
      char filename[42];
      sprintf(filename, "blendSplinePoints_%d.dat", j);
      EG_tecplotSplinePoints(ncp[j].ncp, nsec, xyzs, filename);
    }
#endif


    /* get the BSpline surface */
    stat = EG_loft2spline(outLevel, vknot, vdata,
                          sides, &ncp[j], nsec, xyzs, 1.e-8,
                          blendsurfs[j].header, &blendsurfs[j].data);
    EG_free(xyzs); xyzs = NULL;
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Strip %d splined = %d (EG_blend)!\n", j+1, stat);
      goto cleanup;
    }

#ifdef SPLINE_TECPLOT_DEBUG
    {
      char filename[42];
      sprintf(filename, "blendSpline_%d.dat", j);
      EG_tecplotSpline(blendsurfs[j].header, blendsurfs[j].data, filename);
    }
#endif
  }

  if (ntip > 0) {

    if (nstripe == 2 || te == 2) {
      neigbr[0] = &blendsurfs[0];
      neigbr[1] = &blendsurfs[1];
    } else if (te == 1) {
      neigbr[0] = &blendsurfs[2];
      neigbr[1] = &blendsurfs[0];
    } else {
      neigbr[0] = &blendsurfs[1];
      neigbr[1] = &blendsurfs[2];
    }

    if (tip & 1) {
      stat = EG_wingTipSpline<T>(outLevel, rc1[1], secs[0], secs[1],
                                 0.0, neigbr, &tipsurfs[0]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    if (tip & 2) {
      stat = EG_wingTipSpline<T>(outLevel, rcN[1], secs[nsec-1], secs[nsec-2],
                                 1.0, neigbr, &tipsurfs[1]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    if (nstripe == 3) {

      j = te;

      if (tip & 1) {
        /* get the v = 0 on the first tip */
        stat = EG_isoCurve(tipsurfs[0].header, tipsurfs[0].data,
                           -1, 0,
                           tipcurvs[0].header, &tipcurvs[0].data);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      if (tip & 2) {
        /* get the v = 0 on the last tip */
        stat = EG_isoCurve(tipsurfs[1].header, tipsurfs[1].data,
                           -1, 0,
                           tipcurvs[1].header, &tipcurvs[1].data);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      nKnotv  = nsec + ntip*NTIPTE;
      xyzs    = (T *) EG_alloc(3*(ncp[j].ncp+2)*nKnotv*sizeof(T));
      vknotTE = (T *) EG_alloc(nKnotv*sizeof(T));
      vdataTE = (int *) EG_alloc(nKnotv*sizeof(int));
      if ((xyzs == NULL) ||
          (vknotTE == NULL) ||
          (vdataTE == NULL)) {
        if (outLevel > 0)
          printf(" EGADS Error: Allocation (EG_blend)!\n");
        stat = EGADS_MALLOC;
        goto cleanup;
      }
      t1 = &xyzs[3* ncp[j].ncp   *nKnotv];
      tN = &xyzs[3*(ncp[j].ncp+1)*nKnotv];

      dv0 = 0.0;
      dv1 = 0.0;
      ii = tip & 1 ? NTIPTE : 0;
      for (i = 0; i < ii; i++) {
        vknotTE[i] = 0.0;
        vdataTE[i] = 1;
      }
      for (i = 0; i < nsec; i++) {
         vdataTE[ii+i] = vdata[i];
       }
      for (i = nsec+ii; i < nKnotv; i++) {
        vknotTE[i] = 0.0;
        vdataTE[i] = 1;
      }

      /* get the spline points across all sections for the current stripe */
#ifdef EGADS_SPLINE_VELS
      stat = EG_secSplinePoints(vels, outLevel, 0, nsec, secs, j, ncp, &n,
                                &xyzs[3*ncp[j].ncp*ii], &t1[3*ii], &tN[3*ii]);
#else
      stat = EG_secSplinePoints(outLevel, 0, nsec, secs, j, ncp, &n,
                                &xyzs[3*ncp[j].ncp*ii], &t1[3*ii], &tN[3*ii]);
#endif
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (tip & 1) {

        snor = (T *) EG_alloc(3*ncp[j].ncp*sizeof(T));
        if (snor == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Allocation (EG_blend)!\n");
          stat = EGADS_MALLOC;
          goto cleanup;
        }
        for (i = 0; i < 3*ncp[j].ncp; i++) snor[i] = 0.0;

        /* repeat points to create C1 section */
        jj = ncp[j].ncp*ii;
        npt = ncp[j].ncp*(ii-1);
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          xyzs[3*npt  ] = xyzs[3*(jj+k)  ];
          xyzs[3*npt+1] = xyzs[3*(jj+k)+1];
          xyzs[3*npt+2] = xyzs[3*(jj+k)+2];
        }

        t1[3*(ii-1)  ] = t1[3*(ii)  ];
        t1[3*(ii-1)+1] = t1[3*(ii)+1];
        t1[3*(ii-1)+2] = t1[3*(ii)+2];

        tN[3*(ii-1)  ] = tN[3*(ii)  ];
        tN[3*(ii-1)+1] = tN[3*(ii)+1];
        tN[3*(ii-1)+2] = tN[3*(ii)+2];

        vknotTE[ii-1] = vknot[0]*dv;
        vdataTE[ii  ] = +2;

        /* get the midpoint knot value for the TE tip curve */
        mknot = tipcurvs[0].data[tipcurvs[0].header[3]/2];
        npt = 0;

        deg    = tipcurvs[0].header[1];
        iknotc = tipcurvs[0].header[3]-2*deg;
        n = 0;
        knotc[n++] = tipcurvs[0].data[deg];
        for (k = 0; k < iknotc-1; k++)
          for (kk = 1; kk <= NSUBTE+1; kk++)
            knotc[n++] = tipcurvs[0].data[k+deg]*(1.-kk/double(NSUBTE+1)) +
                         tipcurvs[0].data[k+1+deg]*(kk/double(NSUBTE+1));

        /* create section points */
        for (i = 0; i < NTIPTE; i++) {
          t = knotc[NTIPTE-1+i];
          stat = EG_spline1dDeriv(tipcurvs[0].header, tipcurvs[0].data, 1, t,
                                  pointL);
          if (stat != EGADS_SUCCESS) goto cleanup;

          t = knotc[NTIPTE-1-i];
          stat = EG_spline1dDeriv(tipcurvs[0].header, tipcurvs[0].data, 1, t,
                                  pointU);
          if (stat != EGADS_SUCCESS) goto cleanup;

          dX[0] = pointU[0] - pointL[0];
          dX[1] = pointU[1] - pointL[1];
          dX[2] = pointU[2] - pointL[2];

          for (k = 0; k < ncp[j].ncp; k++, npt++) {
            t = ncp[j].knots[k];
            xyzs[3*npt  ] = pointL[0] + t*dX[0];
            xyzs[3*npt+1] = pointL[1] + t*dX[1];
            xyzs[3*npt+2] = pointL[2] + t*dX[2];
          }

          t1[3*i  ] =  dX[0];
          t1[3*i+1] =  dX[1];
          t1[3*i+2] =  dX[2];

          tN[3*i  ] = -dX[0];
          tN[3*i+1] = -dX[1];
          tN[3*i+2] = -dX[2];

          if (i == 0) {
            snor[0] = pointL[3];
            snor[1] = pointL[4];
            snor[2] = pointL[5];

            snor[3*(ncp[j].ncp-1)+0] = -pointL[3];
            snor[3*(ncp[j].ncp-1)+1] = -pointL[4];
            snor[3*(ncp[j].ncp-1)+2] = -pointL[5];
          }
        }

        t = knotc[2*NTIPTE-2];
        stat = EG_spline1dDeriv(tipcurvs[0].header, tipcurvs[0].data, 1, t,
                                pointL);
        if (stat != EGADS_SUCCESS) goto cleanup;

        t = knotc[0];
        stat = EG_spline1dDeriv(tipcurvs[0].header, tipcurvs[0].data, 1, t,
                                pointU);
        if (stat != EGADS_SUCCESS) goto cleanup;

        dX[0] = pointU[0] - pointL[0];
        dX[1] = pointU[1] - pointL[1];
        dX[2] = pointU[2] - pointL[2];

        /* compute the normal distances for the tip */
        dv0 = sqrt(DOT(dX,dX));
        dv0 *= 0.5*rc1[1];

        /* set the v-knots for the extended section */
        for (i = 0; i < ii; i++)
          vknotTE[i] = dv0*(knotc[NTIPTE-1+i] - mknot)/(1.0-mknot);
      }

      for (i = 0; i < nsec; i++) {
        vknotTE[ii+i] = vknot[i]*dv+dv0;
      }
      dv += dv0;

      if (tip & 2) {

        nnor = (T *) EG_alloc(3*ncp[j].ncp*sizeof(T));
        if (nnor == NULL) {
          if (outLevel > 0)
            printf(" EGADS Error: Allocation (EG_blend)!\n");
          stat = EGADS_MALLOC;
          goto cleanup;
        }
        for (i = 0; i < 3*ncp[j].ncp; i++) nnor[i] = 0.0;

        npt = ncp[j].ncp*(ii+nsec-1);

        /* repeat points to create C1 section */
        jj = npt;
        npt += ncp[j].ncp;
        for (k = 0; k < ncp[j].ncp; k++, npt++) {
          xyzs[3*npt  ] = xyzs[3*(jj+k)  ];
          xyzs[3*npt+1] = xyzs[3*(jj+k)+1];
          xyzs[3*npt+2] = xyzs[3*(jj+k)+2];
        }

        t1[3*(ii+nsec)  ] = t1[3*(ii+nsec-1)  ];
        t1[3*(ii+nsec)+1] = t1[3*(ii+nsec-1)+1];
        t1[3*(ii+nsec)+2] = t1[3*(ii+nsec-1)+2];

        tN[3*(ii+nsec)  ] = tN[3*(ii+nsec-1)  ];
        tN[3*(ii+nsec)+1] = tN[3*(ii+nsec-1)+1];
        tN[3*(ii+nsec)+2] = tN[3*(ii+nsec-1)+2];

        vknotTE[ii+nsec] = vknot[nsec-1]*dv;
        vdataTE[ii+nsec-1] = -2;

        /* get the midpoint knot value for the TE tip curve */
        mknot = tipcurvs[1].data[tipcurvs[1].header[3]/2];

        deg    = tipcurvs[1].header[1];
        iknotc = tipcurvs[1].header[3]-2*deg;

        n = 0;
        knotc[n++] = tipcurvs[1].data[deg];
        for (k = 0; k < iknotc-1; k++)
          for (kk = 1; kk <= NSUBTE+1; kk++)
            knotc[n++] = tipcurvs[1].data[k+deg]*(1.-kk/double(NSUBTE+1)) +
                         tipcurvs[1].data[k+1+deg]*(kk/double(NSUBTE+1));

        /* create section points */
        for (i = 1; i < NTIPTE; i++) {
          t = knotc[2*NTIPTE-2-i];
          stat = EG_spline1dDeriv(tipcurvs[1].header, tipcurvs[1].data, 1, t,
                                  pointL);
          if (stat != EGADS_SUCCESS) goto cleanup;

          t = knotc[i];
          stat = EG_spline1dDeriv(tipcurvs[1].header, tipcurvs[1].data, 1, t,
                                  pointU);
          if (stat != EGADS_SUCCESS) goto cleanup;

          dX[0] = pointU[0] - pointL[0];
          dX[1] = pointU[1] - pointL[1];
          dX[2] = pointU[2] - pointL[2];

          for (k = 0; k < ncp[j].ncp; k++, npt++) {
            t = ncp[j].knots[k];
            xyzs[3*npt  ] = pointL[0] + t*dX[0];
            xyzs[3*npt+1] = pointL[1] + t*dX[1];
            xyzs[3*npt+2] = pointL[2] + t*dX[2];
          }

          t1[3*(ii+nsec+i)  ] =  dX[0];
          t1[3*(ii+nsec+i)+1] =  dX[1];
          t1[3*(ii+nsec+i)+2] =  dX[2];

          tN[3*(ii+nsec+i)  ] = -dX[0];
          tN[3*(ii+nsec+i)+1] = -dX[1];
          tN[3*(ii+nsec+i)+2] = -dX[2];

          if (i == NTIPTE-1) {
            nnor[0] = pointL[3];
            nnor[1] = pointL[4];
            nnor[2] = pointL[5];

            nnor[3*(ncp[j].ncp-1)+0] = -pointL[3];
            nnor[3*(ncp[j].ncp-1)+1] = -pointL[4];
            nnor[3*(ncp[j].ncp-1)+2] = -pointL[5];
          }
        }

        t = knotc[2*NTIPTE-2];
        stat = EG_spline1dDeriv(tipcurvs[1].header, tipcurvs[1].data, 1, t,
                                pointL);
        if (stat != EGADS_SUCCESS) goto cleanup;

        t = knotc[0];
        stat = EG_spline1dDeriv(tipcurvs[1].header, tipcurvs[1].data, 1, t,
                                pointU);
        if (stat != EGADS_SUCCESS) goto cleanup;

        dX[0] = pointU[0] - pointL[0];
        dX[1] = pointU[1] - pointL[1];
        dX[2] = pointU[2] - pointL[2];

        /* compute the normal distances for the tip */
        dv1 = sqrt(DOT(dX,dX));
        dv1 *= 0.5*rcN[1];

        for (i = ii+nsec+1; i < nKnotv; i++)
          vknotTE[i] = dv1*knotc[i-(ii+nsec)]/mknot+dv;
        dv += dv1;
      }

      /* normalize the knots */
      dx = vknotTE[0];
      dy = vknotTE[nKnotv-1]-vknotTE[0];
      for (i = 0; i < nKnotv; i++) {
        vknotTE[i] = (vknotTE[i] - dx)/dy;
      }

      if (tip & 1) {
        snor[0] *= dv/(2*dv0);
        snor[1] *= dv/(2*dv0);
        snor[2] *= dv/(2*dv0);

        snor[3*(ncp[j].ncp-1)+0] *= dv/(2*dv0);
        snor[3*(ncp[j].ncp-1)+1] *= dv/(2*dv0);
        snor[3*(ncp[j].ncp-1)+2] *= dv/(2*dv0);
      }

      if (tip & 2) {
        nnor[0] *= dv/(2*dv1);
        nnor[1] *= dv/(2*dv1);
        nnor[2] *= dv/(2*dv1);

        nnor[3*(ncp[j].ncp-1)+0] *= dv/(2*dv1);
        nnor[3*(ncp[j].ncp-1)+1] *= dv/(2*dv1);
        nnor[3*(ncp[j].ncp-1)+2] *= dv/(2*dv1);
      }

#ifdef SPLINE_TECPLOT_DEBUG
      EG_tecplotSplinePoints(ncp[j].ncp, nKnotv, xyzs, "pointsTE.dat");
#endif

      /* get the BSpline surface */
      stat = EG_spline2dAppr<T>(1, ncp[j].ncp, nKnotv, xyzs, NULL, vknotTE,
                                vdataTE, t1, tN, NULL, snor, NULL, nnor,
                                1.e-8, blendsurfs[j].header, &blendsurfs[j].data);
      EG_free(xyzs); xyzs = NULL;
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: EG_spline2dAppr = %d (EG_wingTip)!\n",
                 stat);
        goto cleanup;
      }

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSpline_%d.dat", j);
        EG_tecplotSpline(blendsurfs[j].header, blendsurfs[j].data, filename);
      }
#endif
    }
  }


  for (j = 0; j < nstripe; j++) {
    for (i = 0; i < nsecC0-1; i++) {

      if ((j == te) && (tip & 1)) {
        if ( i > 0 ) aggC0[2*i] += NTIPTE;
        aggC0[2*i+1] += NTIPTE;
      }
      if ((j == te) && (i == nsecC0-2) && (tip & 2)) {
        aggC0[2*i+1] += NTIPTE;
      }

      /* get the sub-BSpline surface */
      stat = EG_subSpline2d(blendsurfs[j].header, blendsurfs[j].data,
                            1, blendsurfs[j].header[2]-1,
                            aggC0[2*i], aggC0[2*i+1],
                            surfs[i+j*isrf].header, &surfs[i+j*isrf].data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Strip %d splined = %d (EG_blend)!\n",
                 j+1, stat);
        goto cleanup;
      }

      if ((j == te) && (tip & 1)) {
        if ( i > 0 ) aggC0[2*i] -= NTIPTE;
        aggC0[2*i+1] -= NTIPTE;
      }
      if ((j == te) && (i == nsecC0-2) && (tip & 2)) {
        aggC0[2*i+1] -= NTIPTE;
      }


#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSubSpline_%dx%d.dat", i,j);
        EG_tecplotSpline(surfs[i+j*isrf].header, surfs[i+j*isrf].data, filename);
      }
#endif

      /* get the u = 0 curveU[i,j] */
      if (j == te) {
        /* don't get the u = 0 on the TE surface */
        if (j > 0) {
          /* get the u = 1 on the previous surface */
          stat = EG_isoCurve(surfs[i+(j-1)*isrf].header, surfs[i+(j-1)*isrf].data,
                             surfs[i+(j-1)*isrf].header[2]-1, -1,
                             curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        } else {
          /* get the u = 1 on the last surface */
          stat = EG_isoCurve(surfs[i+(nstripe-1)*isrf].header,
                             surfs[i+(nstripe-1)*isrf].data,
                             surfs[i+(nstripe-1)*isrf].header[2]-1, -1,
                             curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      } else {
        stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                           0, -1,
                           curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
      }

      if (j == nstripe-1) {
        if (uclosed == 1) {
          /* uclosed loop so first and last curves are the same */
        } else {
          /* get the u = 1 curveU[i,j+1] */
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             surfs[i+j*isrf].header[2]-1, -1,
                             curvsU[i+(j+1)*icrvU].header,
                             &curvsU[i+(j+1)*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }

      /* get the v = 0 curveV[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data, -1, 0,
                         curvsV[i+j*icrvV].header, &curvsV[i+j*icrvV].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      /* get the v = 1 curveV[i+1,j] on the last section */
      if ( i == nsecC0-2 ) {
        if (vclosed == 1) {
          /* vclosed loop so first and last curves are the same */
        } else {
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             -1, surfs[i+j*isrf].header[5]-1,
                             curvsV[(i+1)+j*icrvV].header,
                             &curvsV[(i+1)+j*icrvV].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve+ = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }
    }
  }
#endif

  for (i = 0; i < nsecC0; i++) {
    /* get the nodes of this section */
    stat = EG_getSecNodes( outLevel, inode, secsC0[i], nstripe, &nodes[i] );
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (vclosed == 1) {
    for (j = 0; j < nstripe+1; j++)
      nodes[nsecC0-1+j*inode] = nodes[0+j*inode];
  }

  stat = EGADS_SUCCESS;

cleanup:
  /* clean up all of our temps */
  EG_free(aggC0);
  EG_free(vknot);
  EG_free(vknotC0);
  EG_free(vknotTE);
  EG_free(vdata);
  EG_free(vdataTE);
  EG_free(xyzs);
  EG_free(snor);
  EG_free(nnor);
#ifndef BLEND_SPLIT_CONSTRUCTION
  delete [] blendsurfs;
#endif
  if (ncp != NULL) EG_freeSeq(nstripe, ncp);
  if (stat != EGADS_SUCCESS) {
    EG_free(secsC0);
    EG_free(nodes);
    delete [] surfs;
    delete [] curvsU;
    delete [] curvsV;
  } else {
    *nsecC0_out = nsecC0;
    *secsC0_out = secsC0;
    *nodes_out = nodes;
    *curvsU_out = curvsU;
    *curvsV_out = curvsV;
    *surfs_out = surfs;
  }

  return stat;
}


template<class T>
static int
EG_ruledSpline(
#ifdef EGADS_SPLINE_VELS
               const egadsSplineVels *vels,
#endif
               int nsec, const ego *secs, int nstripe, int uclosed, int vclosed,
               ego **nodes_out, egSpline<T> **curvsU_out,
               egSpline<T> **curvsV_out, egSpline<T> **surfs_out)
{
  int    i, j, outLevel, stat, isrf, inode, icrvU, icrvV;
  int    nnode=0, ncurvU=0, ncurvV, nsurf=0, planar;
  T      t1[6], tN[6], *xyzs=NULL;
  ego    *nodes=NULL;
  egSequ<T> *ncp=NULL;
  egSpline<T> *surfs=NULL, *curvsV=NULL, *curvsU=NULL;

  *nodes_out  = NULL;
  *curvsU_out = NULL;
  *curvsV_out = NULL;
  *surfs_out  = NULL;

  outLevel = EG_outLevel(secs[0]);

  /* get the sampling points for each stripe */
#ifdef EGADS_SPLINE_VELS
  stat = EG_setSequence(vels, nsec, secs, 0, nstripe, &ncp);
#else
  stat = EG_setSequence(nsec, secs, 0, nstripe, &ncp);
#endif
  if (stat != EGADS_SUCCESS) goto cleanup;

  inode = nsec;
  icrvU = nsec-1;
  icrvV = nsec;
  isrf = nsec-1;

  nnode  = (nstripe+1)*nsec;
  ncurvU = (nstripe+1)*icrvU;
  ncurvV = nstripe*icrvV;
  nsurf  = nstripe*isrf;

  /* make a nodes, curves, and surfaces */
  nodes = (ego*)EG_alloc( nnode *sizeof(ego));
  curvsU = new egSpline<T>[ncurvU];
  curvsV = new egSpline<T>[ncurvV];
  surfs = new egSpline<T>[nsurf];
  if ((nodes == NULL) || (curvsU == NULL) || (curvsV == NULL) ||
      (surfs == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation (EG_ruled)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  for (j = 0; j < nnode; j++) nodes[j] = NULL;

  for (i = 0; i < nsec; i++) {
    /* get the nodes of this section */
    stat = EG_getSecNodes( outLevel, inode, secs[i], nstripe, &nodes[i] );
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* propagate equivalent nodes */
  for (j = 0; j < nstripe+1; j++) {
    for (i = 0; i < nsec-1; i++) {
      if (EG_isEquivalent(nodes[i+j*inode], nodes[i+1+j*inode]) == EGADS_SUCCESS) {
        nodes[i+1+j*inode] = nodes[i+j*inode];
      }
    }
  }

  if (vclosed == 1) {
    for (j = 0; j < nstripe+1; j++)
      nodes[0+j*inode] = nodes[nsec-1+j*inode];
  }


  for (j = 0; j < nstripe; j++) {
    xyzs = (T *) EG_alloc(6*ncp[j].ncp*sizeof(T));
    if (xyzs == NULL)  {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation for %dX%d points (EG_ruled)!\n",
               ncp[j].ncp, nsec);
      stat = EGADS_MALLOC;
      goto cleanup;
    }

    /* first construct all surfaces (skipping degenerate ones) */
    for (i = 0; i < nsec-1; i++) {

      /* skip multipicty of edges */
      stat = EG_secEquivEdge(outLevel, &secs[i], j);
      if (stat < EGADS_SUCCESS) goto cleanup;
      if (stat == EGADS_SUCCESS) continue;

      /* get the spline points between pairs of sections for the current stripe */
#ifdef EGADS_SPLINE_VELS
      stat = EG_secSplinePoints(vels, outLevel, i, 2, &secs[i], j, ncp, &planar, xyzs,
                                t1, tN);
#else
      stat = EG_secSplinePoints(outLevel, i, 2, &secs[i], j, ncp, &planar, xyzs,
                                t1, tN);
#endif
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the BSpline surface */
      if ((planar == 1) || (ncp[j].ncp == 3)) {
        stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                  1.e-8, surfs[i+j*isrf].header,
                                  &surfs[i+j*isrf].data);
      } else {
        stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                  NULL, t1, tN, NULL, NULL, NULL, NULL,
                                  1.e-8, surfs[i+j*isrf].header,
                                  &surfs[i+j*isrf].data);
      }
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d spline2d = %d (EG_ruled)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "ruleSpline_%dx%d.dat", i,j);
        EG_tecplotSpline(surfs[i+j*isrf].header, surfs[i+j*isrf].data, filename);
      }
#endif
    }

    /* construct the curves from the surfaces */
    for (i = 0; i < nsec-1; i++) {

      if (surfs[i+j*isrf].data == NULL) continue;

      if ((nodes[i+j*inode] == nodes[i+1+j*inode]) ||
          (nodes[i+j*inode] == nodes[0  +j*inode] && i == nsec-2 && vclosed == 1)) {
        /* degenerate or periodic and degenerate */
      } else {
        /* get the u = 0 curveU[i,j] */
        stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                           0, -1,
                           curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_ruled)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
      }

      if (j == nstripe-1) {
        if (uclosed == 1 || nodes[i+(j+1)*inode] == nodes[i+1+(j+1)*inode]) {
          /* uclosed loop so first and last curves are the same, or the end is degenerate */
        } else {
          /* get the u = 1 curveU[i,j+1] */
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             surfs[i+j*isrf].header[2]-1, -1,
                             curvsU[i+(j+1)*icrvU].header,
                             &curvsU[i+(j+1)*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_ruled)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }

      /* get the v = 0 curveV[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                         -1, 0,
                         curvsV[i+j*icrvV].header, &curvsV[i+j*icrvV].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_ruled)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      /* get the v = 1 curveV[i+1,j] on the last section, or if the next section is NULL */
      if (i == nsec-2 || surfs[MIN((i+1),nsec-2)+j*isrf].data == NULL) {
        if (vclosed == 1 && i == nsec-2) {
          /* first and last section is the same to create uclosed geometry */
        } else {
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             -1, surfs[i+j*isrf].header[5]-1,
                             curvsV[(i+1)+j*icrvV].header,
                             &curvsV[(i+1)+j*icrvV].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve+ = %d (EG_ruled)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }
    }
    EG_free(xyzs); xyzs = NULL;
  }

  stat = EGADS_SUCCESS;

cleanup:
  /* clean up all of our temps */
  EG_free(xyzs);
  if (ncp != NULL) EG_freeSeq(nstripe, ncp);
  if (stat != EGADS_SUCCESS) {
    EG_free(nodes);
    delete [] surfs;
    delete [] curvsU;
    delete [] curvsV;
  } else {
    *nodes_out  = nodes;
    *curvsU_out = curvsU;
    *curvsV_out = curvsV;
    *surfs_out  = surfs;
  }

  return stat;
}


static int
EG_splineGeom(int nsec, const ego *secs, int nstripe, int uclosed, int vclosed,
              int ncap, int te, int tip, ego *nodes, egSpline<double> *splcurvsU,
              egSpline<double> *splcurvsV, egSpline<double> *splsurfs,
              int *ncurvV_out, ego **edgesU_out, ego **edgesV_out,
              int *nface_out, ego **faces_out, objStack *stack)
{
  int    i, j, n, outLevel, stat, degenSurf = 0;
  int    ncurvU=0, ncurvV=0, nface=0, inode, icrvU, icrvV, isrf;
  int    esens[4], lsens[1] = {SFORWARD};
  double data[18], ts[2];
  ego    context, loop, ref, nds[2], edges[8];
  ego    pcurvs[4]={NULL,NULL,NULL,NULL};
  ego    *curvsU=NULL, *curvsV=NULL, *surfs=NULL;
  ego    *edgesU=NULL, *edgesV=NULL, *faces=NULL;

  *ncurvV_out = 0;
  *edgesU_out = NULL;
  *edgesV_out = NULL;
  *nface_out  = 0;
  *faces_out  = NULL;

  /* Geometric entities are laid out as follows
   *
   *
   *                           curvsV[i+1,j]
   *  sec[i+1]  nodes[i+1,j]--------------------- nodes[i+1,j+1]
   *              |                                 |
   *              |                                 |
   *              |                                 |
   *  curvsU[i,j] |   v        surfs[i,j]           | curvsU[i,j+1]
   *              |   ^                             |
   *              |   |                             |
   *              |   +--> u                        |
   *              |                                 |
   *  sec[i]    nodes[i,j]----------------------- nodes[i,j+1]
   *                           curvsV[i,j]
   *
   *                           stripe[j]
   */

  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL) return EGADS_NOTCNTX;

  /* setup strides */
  inode = nsec;
  icrvU = nsec-1;
  icrvV = nsec;
  isrf  = nsec-1;

  ncurvU = (nstripe+1)*icrvU;
  ncurvV = nstripe*icrvV;
  nface  = nstripe*isrf + ncap;

  curvsU = (ego *) EG_alloc(ncurvU*sizeof(ego));
  curvsV = (ego *) EG_alloc(ncurvV*sizeof(ego));
  edgesU = (ego *) EG_alloc(ncurvU*sizeof(ego));
  edgesV = (ego *) EG_alloc(ncurvV*sizeof(ego));
  surfs =  (ego *) EG_alloc((nface-ncap)*sizeof(ego));
  faces =  (ego *) EG_alloc(nface*sizeof(ego));
  if ((curvsU == NULL) ||
      (curvsV == NULL) ||
      (edgesU == NULL) ||
      (edgesV == NULL) ||
      (surfs == NULL) ||
      (faces == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for %d Surfaces (EG_splineGeom)!\n",
             nface);
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  for (i = 0; i < ncurvU; i++) {
    curvsU[i] = NULL;
    edgesU[i] = NULL;
  }
  for (i = 0; i < ncurvV; i++) {
    curvsV[i] = NULL;
    edgesV[i] = NULL;
  }
  for (i = 0; i < nface-ncap; i++) surfs[i] = NULL;
  for (i = 0; i < nface     ; i++) faces[i] = NULL;

  ts[0] = 0.0;
  ts[1] = 1.0;

  /* construct ego geometry from spline information */
  for (j = 0; j < nstripe; j++) {
    degenSurf = 0;
    for (i = 0; i < nsec-1; i++) {

      if (splsurfs[i+j*isrf].data != NULL) {

        stat = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                               splsurfs[i+j*isrf].header, splsurfs[i+j*isrf].data,
                               &surfs[i+j*isrf]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Surface %d/%d EG_makeGeometry = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
#ifdef SPLINE_TECPLOT_DEBUG
          char filename[42];
          sprintf(filename, "spline_%dx%d_fail.dat", i,j);
          EG_tecplotSpline(splsurfs[i+j*isrf].header, splsurfs[i+j*isrf].data,
                           filename);
#endif
          goto cleanup;
        }
        stat = EG_stackPush(stack, surfs[i+j*isrf]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      /*-------------------------*/
      /*      curvsU[i,j]        */
      /*-------------------------*/

      if (splcurvsU[i+j*icrvU].data == NULL) {

        /* degenerate EDGE */
        n = 1;
        nds[0] = nodes[i+  j*inode];

        /* make the edge on the u-constant curve */
        stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, n,
                               nds, NULL, &edgesU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Degenerate Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
        stat = EG_stackPush(stack, edgesU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) goto cleanup;

      } else {

        /* make the curve that spans the sections */
        stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                               splcurvsU[i+j*icrvU].header,
                               splcurvsU[i+j*icrvU].data, &curvsU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Curve %d/%d EG_makeGeometry = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
        stat = EG_stackPush(stack, curvsU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        n = 2;
        nds[0] = nodes[i+  j*inode];
        nds[1] = nodes[i+1+j*inode];

        /* make the edge on the u-constant curve */
        stat = EG_makeTopology(context, curvsU[i+j*isrf], EDGE, TWONODE, ts, n,
                               nds, NULL, &edgesU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
        stat = EG_stackPush(stack, edgesU[i+j*icrvU]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      /*-------------------------*/
      /*      curvsU[i,j+1]      */
      /*-------------------------*/

      if (j == nstripe-1) {
        if (uclosed == 1) {
          /* uclosed loop so first and last curves are the same */
          edgesU[i+(j+1)*icrvU] = edgesU[i+0*icrvU];
        } else {

          if (splcurvsU[i+(j+1)*icrvU].data == NULL) {

            /* degenerate EDGE */
            n = 1;
            nds[0] = nodes[i  +(j+1)*inode];

            /* make the edge on the curve */
            stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE,
                                   ts, n, nds, NULL, &edgesU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Degnerate Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                       j+1, i+1, stat);
              goto cleanup;
            }
            stat = EG_stackPush(stack, edgesU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) goto cleanup;

          } else {
            /* make the last curve that spans the sections */
            stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                   splcurvsU[i+(j+1)*icrvU].header,
                                   splcurvsU[i+(j+1)*icrvU].data,
                                   &curvsU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Edge %d/%d EG_makeGeometry = %d (EG_splineGeom)!\n",
                       j+1, i+1, stat);
              goto cleanup;
            }
            stat = EG_stackPush(stack, curvsU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) goto cleanup;

            n = 2;
            nds[0] = nodes[i  +(j+1)*inode];
            nds[1] = nodes[i+1+(j+1)*inode];

            /* make the edge on the curve */
            stat = EG_makeTopology(context, curvsU[i+(j+1)*icrvU], EDGE, TWONODE,
                                   ts, n, nds, NULL, &edgesU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                       j+1, i+1, stat);
              goto cleanup;
            }
            stat = EG_stackPush(stack, edgesU[i+(j+1)*icrvU]);
            if (stat != EGADS_SUCCESS) goto cleanup;
          }
        }
      }

      if (splsurfs[i+j*isrf].data == NULL) {
        /* degenerate surface: copy the edge from the previous surface */
        edgesV[(i+1)+j*icrvV] = edgesV[i+j*icrvV];
        if (vclosed && i == nsec-2) {
          edgesV[0+j*icrvV] = edgesV[i+j*icrvV];
        }
        degenSurf++;
        continue;
      }

      /*-------------------------*/
      /*      curvsV[i,j]        */
      /*-------------------------*/

      /* don't make v constant curves on the TE surface tips */
      if ( !((j == te) && (i == 0) && (tip & 1)) ) {
        /* make the curve on the current section/stripe */
        stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                               splcurvsV[i+j*icrvV].header,
                               splcurvsV[i+j*icrvV].data,
                                 &curvsV[i+j*icrvV]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d EG_makeGeometry = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
        stat = EG_stackPush(stack, curvsV[i+j*icrvV]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        if (secs[i]->oclass == NODE) {
          n = 1;
          nds[0] = nodes[i+j*inode];

          /* make the DEGENERATE edge on the node */
          stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, n,
                                 nds, NULL, &edgesV[i+j*icrvV]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        } else if (nstripe == 1 && uclosed == 1) {
          n = 1;
          nds[0] = nodes[i+ j   *inode];

          /* make the ONENODE edge on the curve */
          stat = EG_makeTopology(context, curvsV[i+j*icrvV], EDGE, ONENODE, ts,
                                 n, nds, NULL, &edgesV[i+j*icrvV]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        } else {
          n = 2;
          nds[0] = nodes[i+ j   *inode];
          nds[1] = nodes[i+(j+1)*inode];

          /* make the TWONODE edge on the curve */
          stat = EG_makeTopology(context, curvsV[i+j*icrvV], EDGE, TWONODE, ts,
                                 n, nds, NULL, &edgesV[i+j*icrvV]);
          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
        stat = EG_stackPush(stack, edgesV[i+j*icrvV]);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      /*-------------------------*/
      /*      curvsV[i+1,j]      */
      /*-------------------------*/

      if ( !((j == te) && (i == nsec-2) && (tip & 2)) ) {
        /* make the curve on the last section/stripe */
        if (i == nsec-2 || splsurfs[MIN((i+1),nsec-2)+j*isrf].data == NULL) {
          if (vclosed == 1 && i == nsec-2) {
            /* first and last sections are repeated to close the geometry */
            edgesV[(i+1)+j*icrvV] = edgesV[0+j*icrvV];
          } else {
            stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                                   splcurvsV[(i+1)+j*icrvV].header,
                                   splcurvsV[(i+1)+j*icrvV].data,
                                   &curvsV[(i+1)+j*icrvV]);
            if (stat != EGADS_SUCCESS) {
              if (outLevel > 0)
                printf(" EGADS Error: Edge %d/%d EG_makeGeometry = %d (EG_splineGeom)!\n",
                       j+1, i+1, stat);
              goto cleanup;
            }
            stat = EG_stackPush(stack, curvsV[(i+1)+j*icrvV]);
            if (stat != EGADS_SUCCESS) goto cleanup;

            if (secs[i+1]->oclass == NODE) {
              n = 1;
              nds[0] = nodes[i+1+j*inode];

              /* make the DEGENERATE edge on the node */
              stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, n,
                                     nds, NULL, &edgesV[(i+1)+j*icrvV]);
              if (stat != EGADS_SUCCESS) {
                if (outLevel > 0)
                  printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                         j+1, i+1, stat);
                goto cleanup;
              }
            } else if (nstripe == 1 && uclosed == 1) {
              n = 1;
              nds[0] = nodes[i+1+ j   *inode];

              /* make the ONENODE edge on the curve */
              stat = EG_makeTopology(context, curvsV[(i+1)+j*icrvV], EDGE, ONENODE,
                                     ts, n, nds, NULL, &edgesV[(i+1)+j*icrvV]);
              if (stat != EGADS_SUCCESS) {
                if (outLevel > 0)
                  printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                         j+1, i+1, stat);
                goto cleanup;
              }
            } else {
              n = 2;
              nds[0] = nodes[i+1+ j   *inode];
              nds[1] = nodes[i+1+(j+1)*inode];

              /* make the TWONODE edge on the curve */
              stat = EG_makeTopology(context, curvsV[(i+1)+j*icrvV], EDGE, TWONODE,
                                     ts, n, nds, NULL, &edgesV[(i+1)+j*icrvV]);
              if (stat != EGADS_SUCCESS) {
                if (outLevel > 0)
                  printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                         j+1, i+1, stat);
                goto cleanup;
              }
            }
            stat = EG_stackPush(stack, edgesV[(i+1)+j*icrvV]);
            if (stat != EGADS_SUCCESS) goto cleanup;
          }
        }
      }
    }

    if (degenSurf > 0) {
      if (degenSurf == nsec-1) {
        i = nsec-2;
        int esens;
        /* all surfaces are degenerate, use the edge from the first section */
        stat = EG_getSecEdge(outLevel, secs[0], j, &ref, &esens);
        if (stat != EGADS_SUCCESS) goto cleanup;
        if (ref->oclass != EDGE) {
          stat = EGADS_TOPOERR;
          if (outLevel > 0)
            printf(" EGADS Error: Stripe %d section 1 is not and Edge (EG_splineGeom)!\n",
                   j+1);
          goto cleanup;
        }
        egadsEdge *pedge = (egadsEdge *) ref->blind;

        if (esens >= SFORWARD) {
          nds[0] = nodes[i+1+ j   *inode];
          nds[1] = nodes[i+1+(j+1)*inode];
        } else {
          nds[0] = nodes[i+1+(j+1)*inode];
          nds[1] = nodes[i+1+ j   *inode];
        }
        n = 2;

        /* make the TWONODE edge on the curve */
        stat = EG_makeTopology(context, pedge->curve, EDGE, TWONODE,
                               pedge->trange, n, nds, NULL, &edgesV[(i+1)+j*icrvV]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
        stat = EG_stackPush(stack, edgesV[(i+1)+j*icrvV]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        
        if (esens <= SREVERSE) {
          ego edge = edgesV[(i+1)+j*icrvV];
          stat = EG_flipObject(edge, &edgesV[(i+1)+j*icrvV]);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }
      }

      /* copy over edgesV for degenerate surfaces */
      for (i = nsec-2; i >= 0; i--) {
        if (splsurfs[i+j*isrf].data == NULL) {
          edgesV[i+j*icrvV] = edgesV[(i+1)+j*icrvV];
        }
      }
    }
  }


  data[0] = 0.; data[1] = 0.; /* u == 0 UMIN */
  data[2] = 0.; data[3] = 1.;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(stack, pcurvs[0]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.; data[1] = 0.; /* v == 0 VMIN */
  data[2] = 1.; data[3] = 0.;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(stack, pcurvs[1]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  data[0] = 1.; data[1] = 0.; /* u == 1 UMAX */
  data[2] = 0.; data[3] = 1.;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[2]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(stack, pcurvs[2]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  data[0] = 0.; data[1] = 1.; /* v == 1 VMAX */
  data[2] = 1.; data[3] = 0.;
  stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[3]);
  if (stat != EGADS_SUCCESS) goto cleanup;
  stat = EG_stackPush(stack, pcurvs[3]);
  if (stat != EGADS_SUCCESS) goto cleanup;

  esens[UMIN] = SREVERSE;
  esens[VMIN] = SFORWARD;
  esens[UMAX] = SFORWARD;
  esens[VMAX] = SREVERSE;

  edges[UMIN+4] = pcurvs[0];
  edges[VMIN+4] = pcurvs[1];
  edges[UMAX+4] = pcurvs[2];
  edges[VMAX+4] = pcurvs[3];

  /* construct faces */
  for (j = 0; j < nstripe; j++) {
    for (i = 0; i < nsec-1; i++) {

      if (surfs[i+j*isrf] == NULL) continue;

      /* skip trailing edge face here */
      if ((j == te) && (i == 0     ) && (tip & 1)) continue;
      if ((j == te) && (i == nsec-2) && (tip & 2)) continue;

      edges[UMIN] = edgesU[i+  j   *icrvU];
      edges[VMIN] = edgesV[i+  j   *icrvV];
      edges[UMAX] = edgesU[i+ (j+1)*icrvU];
      edges[VMAX] = edgesV[i+1+j   *icrvV];
      ref         = surfs[i+j*isrf];

      /* make the loop from the edges */
      stat = EG_makeTopology(context, ref, LOOP, CLOSED, NULL, 4, edges, esens,
                             &loop);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Loop %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }
      stat = EG_stackPush(stack, loop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the face from the loop */
      stat = EG_makeTopology(context, ref, FACE, SFORWARD, NULL, 1, &loop, lsens,
                             &faces[i+j*isrf]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d/%d EG_makeTopology = %d (EG_splineGeom)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }
      stat = EG_stackPush(stack, faces[i+j*isrf]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }
  }

  /* set outputs */
  *ncurvV_out = ncurvV;
  *edgesU_out = edgesU; edgesU = NULL;
  *edgesV_out = edgesV; edgesV = NULL;
  *nface_out  = nface;
  *faces_out  = faces;  faces = NULL;

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_splineGeom)!\n", stat);
  }

  /* clean up all of our temps */
  EG_free(faces);
  EG_free(surfs);
  EG_free(edgesU);
  EG_free(edgesV);
  EG_free(curvsU);
  EG_free(curvsV);

  return stat;
}


static int
EG_splineGeom_dot(ego body, int nsec, const ego *secs, int nstripe,
                  int uclosed, int vclosed, int ncap, int te, int tip, ego *nodes,
                  egSpline< SurrealS<1> > *splcurvsU,
                  egSpline< SurrealS<1> > *splcurvsV,
                  egSpline< SurrealS<1> > *splsurfs)
{
#ifdef PCURVE_SENSITIVITY
                         /* point,  direction */
  double rpcurve[4][4] = { { 0, 0,  0, 1},   /* UMIN */
                           { 0, 0,  1, 0},   /* VMIN */
                           { 1, 0,  0, 1},   /* UMAX */
                           { 0, 1,  1, 0} }; /* VMAX */

  double rpcurve_dot[4] = {0,0,0,0}; /* P-curves do not use paramters */
#endif

  int    i, j, k, outLevel, stat, oclass, mtype, fsens;
  int    nchldrn, nedge, nface, inode, icrvU, icrvV, isrf;
  int    *senses, *esens, iedge[4];
  SurrealS<1> ts[2] = {0, 1};
  double data[18];
  ego    surf, curv, loop, shell, ref;
  ego    *chldrn, *nds, *edges, *faces;

  outLevel = EG_outLevel(body);

  /* clear any old sensitivities */
  stat = EG_setGeometry_dot(body, BODY, 0, NULL, NULL, NULL);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getTopology = %d (EG_splineGeom_dot)!\n", stat);
    goto cleanup;
  }

  /* get the shell */
  stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                        &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;

  if (mtype == FACEBODY) {
    /* the shell is a face for a FACEBODY */
    nface = 1;
    faces = chldrn;
  } else {
    /* get the faces from the shell */
    shell = chldrn[0];
    stat = EG_getTopology(shell, &ref, &oclass, &mtype, data, &nface, &faces,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* setup strides */
  inode = nsec;
  icrvU = nsec-1;
  icrvV = nsec;
  isrf  = nsec-1;

  if ( nface != nstripe*isrf + ncap ) {
    if (outLevel > 0)
      printf(" EGADS Error: Face Count = %d is not %d (EG_splineGeom_dot)!\n",
             nface, nstripe*isrf + ncap);
    return EGADS_TOPOERR;
  }

  for (j = 0; j < nstripe; j++) {
    for (i = 0; i < nsec-1; i++) {

      /* get the surface and loop for the face */
      stat = EG_getTopology(faces[i+j*isrf], &ref, &oclass, &fsens, data,
                            &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nchldrn != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d Loop Count = %d is not 1 (EG_splineGeom_dot)!\n",
                 i+j*isrf, nchldrn);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }
      surf = ref;
      loop = chldrn[0];

      /* set the sensitivity for the surface */
      stat = EG_setGeometry_dot(surf, SURFACE, BSPLINE,
                                splsurfs[i+j*isrf].header,
                                splsurfs[i+j*isrf].data);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* get the edges for the loop */
      stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nedge, &edges,
                            &esens);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if ( !(((j == te) && (i == 0     ) && (tip & 1)) ||
             ((j == te) && (i == nsec-2) && (tip & 2))) &&
           nedge != 4 ) {
        if (outLevel > 0)
          printf(" EGADS Error: Face %d Edge Count = %d is not 4 (EG_splineGeom_dot)!\n",
                 i+j*isrf, nedge);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      iedge[0] = iedge[1] = iedge[2] = iedge[3] = -1;
      for (k = 0; k < nedge; k++) {
        /* re-discover the edge order in the loop */
        egadsEdge *pedge = (egadsEdge *) edges[k]->blind;

        ego node0 = pedge->nodes[0];
        ego node1 = pedge->nodes[1];

        if ((EG_isEquivalent(node0, nodes[ i   +j*inode]) == EGADS_SUCCESS) &&
            (EG_isEquivalent(node1, nodes[(i+1)+j*inode]) == EGADS_SUCCESS)) {
          /* use sens when VMIN and VMAX are degenerate */
          if (esens[k] == SREVERSE*fsens) iedge[UMIN] = k;
          if (esens[k] == SFORWARD*fsens) iedge[UMAX] = k;
        } else if ((EG_isEquivalent(node0, nodes[i+ j   *inode]) == EGADS_SUCCESS) &&
                   (EG_isEquivalent(node1, nodes[i+(j+1)*inode]) == EGADS_SUCCESS)) {
          iedge[VMIN] = k;
        } else if ((EG_isEquivalent(node0, nodes[ i   +(j+1)*inode]) == EGADS_SUCCESS) &&
                   (EG_isEquivalent(node1, nodes[(i+1)+(j+1)*inode]) == EGADS_SUCCESS)) {
          iedge[UMAX] = k;
        } else if ((EG_isEquivalent(node0, nodes[(i+1)+ j   *inode]) == EGADS_SUCCESS) &&
                   (EG_isEquivalent(node1, nodes[(i+1)+(j+1)*inode]) == EGADS_SUCCESS)) {
          iedge[VMAX] = k;
        }
      }

      /* check that all edges were found, excluding modified trailing edges */
      for (k = 0; k < 4; k++) {
        if (iedge[k] != -1) continue;
        if ( (j == te) && (i == 0     ) && (tip & 1) && (k == VMIN) ) continue;
        if ( (j == te) && (i == nsec-2) && (tip & 2) && (k == VMAX) ) continue;

        printf(" EGADS Error: Sections for sensitivities mismatch (EG_splineGeom_dot)!\n");
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

#ifdef PCURVE_SENSITIVITY
      /* set P-curve sensitivities to zero */
      for (k = 0; k < 4; k++) {
        if ( (j == te) && (i == 0     ) && (tip & 1) ) continue;
        if ( (j == te) && (i == nsec-2) && (tip & 2) ) continue;

        stat = EG_setGeometry_dot(edges[4+iedge[k]], PCURVE, LINE, NULL,
                                  rpcurve[k], rpcurve_dot);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
#endif

      /* get the curv and nodes for the u = 0 edge */
      stat = EG_getTopology(edges[iedge[UMIN]], &curv, &oclass, &mtype, data,
                            &nchldrn, &nds, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* set the sensitivity for the curve */
      stat = EG_setGeometry_dot(curv, CURVE, BSPLINE,
                                splcurvsU[i+j*icrvU].header,
                                splcurvsU[i+j*icrvU].data);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* set the sensitivity for the t-range of the edge */
      stat = EG_setRange_dot(edges[iedge[UMIN]], EDGE, ts);
      if (stat != EGADS_SUCCESS) goto cleanup;

      if (j == nstripe-1 && uclosed == 0) {
        /* get the curv for the u = 1 edge */
        stat = EG_getTopology(edges[iedge[UMAX]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity for the curv */
        stat = EG_setGeometry_dot(curv, CURVE, BSPLINE,
                                  splcurvsU[i+(j+1)*icrvU].header,
                                  splcurvsU[i+(j+1)*icrvU].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity for the t-range of the edge */
        stat = EG_setRange_dot(edges[iedge[UMAX]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      /* don't make v constant curves on the TE surface tips */
      if ( !((j == te) && (i == 0) && (tip & 1)) ) {
        /* get the curv for the v = 0 edge */
        stat = EG_getTopology(edges[iedge[VMIN]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* copy the node sensitivity in the lower left corner (or just once
           if the section is a NODE) */
        if (secs[i]->oclass != NODE || j == 0) {
          stat = EG_copyGeometry_dot(nodes[i+j*inode], NULL, NULL, nds[0]);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        /* copy the last node on the section if it's open */
        if (j == nstripe-1 && uclosed == 0) {
          stat = EG_copyGeometry_dot(nodes[i+(j+1)*inode], NULL, NULL, nds[1]);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        if (mtype != DEGENERATE) {
          /* set the sensitivity for the curv */
          stat = EG_setGeometry_dot(curv, CURVE, BSPLINE,
                                    splcurvsV[i+j*icrvV].header,
                                    splcurvsV[i+j*icrvV].data);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        /* set the sensitivity for the t-range of the edge */
        stat = EG_setRange_dot(edges[iedge[VMIN]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }

      if ( !((j == te) && (i == nsec-2) && (tip & 2)) ) {
        if (i == nsec-2 && vclosed == 0) {
          /* get the curv for the v = 1 edge on the last section */
          stat = EG_getTopology(edges[iedge[VMAX]], &curv, &oclass, &mtype, data,
                                &nchldrn, &nds, &senses);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* copy the node sensitivity in the top left corner (or just once
             if the section is a NODE) */
          if (secs[i+1]->oclass != NODE || j == 0) {
            stat = EG_copyGeometry_dot(nodes[i+1+j*inode], NULL, NULL, nds[0]);
            if (stat != EGADS_SUCCESS) goto cleanup;
          }

          /* copy the last node on the section if it's open */
          if (j == nstripe-1 && uclosed == 0) {
            stat = EG_copyGeometry_dot(nodes[i+1+(j+1)*inode], NULL, NULL,
                                       nds[1]);
            if (stat != EGADS_SUCCESS) goto cleanup;
          }

          if (mtype != DEGENERATE) {
            /* set the sensitivity for the curv */
            stat = EG_setGeometry_dot(curv, CURVE, BSPLINE,
                                      splcurvsV[i+1+j*icrvV].header,
                                      splcurvsV[i+1+j*icrvV].data);
            if (stat != EGADS_SUCCESS) goto cleanup;
          }

          /* set the sensitivity for the t-range of the edge */
          stat = EG_setRange_dot(edges[iedge[VMAX]], EDGE, ts);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }
      }
    }
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: status = %d (EG_splineGeom_dot)!\n", stat);
  }

  return stat;
}


template<class T>
static int
EG_blendNodeSpline(
#ifdef EGADS_SPLINE_VELS
               const egadsSplineVels *vels,
#endif
               int nsex, const ego *secs, int vclosed,
               int *nsecC0_out, ego **secsC0_out,
               ego **nodes_out, egSpline<T> **curvsU_out)
{
  int    stat = EGADS_SUCCESS;
  int    i, j, k, n, jj, outLevel;
  int    fixed, nsecC0=0, ncurvU=0;
  int    nsec, rite, left, *vdata=NULL;
  int    *vdataTE=NULL, *aggC0=NULL;
  T      *xyzs=NULL, *vknot=NULL, *vknotC0=NULL, *vknotTE=NULL, max2;
  T      dx, dy;
  T      d2xdt2L, d2ydt2L, d2zdt2L, d2sdt2L, d2xdt2R, d2ydt2R, d2zdt2R, d2sdt2R;
  ego    *nodes=NULL, *secsC0=NULL;
  egSpline<T> *curvsU=NULL;
#ifdef BLEND_SPLIT_CONSTRUCTION
  int nC0, iC0;
#else
  egSpline<T> blendcurve;
#endif

  *nsecC0_out = 0;
  *secsC0_out = NULL;
  *nodes_out  = NULL;
  *curvsU_out = NULL;

  outLevel = EG_outLevel(secs[0]);
  fixed    = EG_fixedKnots(secs[0]);

  nsec = nsex;
  if (nsec < 0) nsec = -nsec;

  /* set the knots in the loft direction */
  vknot   = (T *)   EG_alloc(  nsec*sizeof(T));
  vdata   = (int *) EG_alloc(  nsec*sizeof(int));
  vknotC0 = (T *)   EG_alloc(  nsec*sizeof(T));
  aggC0   = (int *) EG_alloc(2*nsec*sizeof(int));
  xyzs    = (T *)   EG_alloc(3*nsec*sizeof(T));
  secsC0  = (ego *) EG_alloc(  nsec*sizeof(ego));
  nodes   = (ego *) EG_alloc(  nsec*sizeof(ego));
  if ((vknot == NULL)   || (xyzs == NULL)  || (vdata == NULL)   ||
      (vknotC0 == NULL) || (aggC0 == NULL) || (secsC0 == NULL)  ||
      (nodes == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation for (EG_blend)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for (i = 0; i <   nsec; i++) secsC0[i] = NULL;
  for (i = 0; i < 2*nsec; i++) aggC0[i]  = 0;

  for (i = 0; i < nsec; i++) {
    /* get the node of this section */
    stat = EG_getSecNodes( outLevel, 1, secs[i], 0, &nodes[i] );
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (vclosed == 1) {
    nodes[nsec-1] = nodes[0];
  }

  /* get the spline points across all sections */
  for (i = 0; i < nsec; i++) {
#ifdef EGADS_SPLINE_VELS
    stat = EG_nodeVels(vels, outLevel, i, secs, nodes[i], &xyzs[3*i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
#else
    stat = EG_evaluat(nodes[j  ], NULL, &xyzs[3*i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
#endif
  }

  /* arc-length spaced */
  for (i = 0; i < nsec; i++) vknot[i] = 0.0;

  dy = 0.0;
  for (j = 1; j < nsec; j++) {
    dy += sqrt((xyzs[3*j  ]-xyzs[3*(j-1)  ])*
               (xyzs[3*j  ]-xyzs[3*(j-1)  ]) +
               (xyzs[3*j+1]-xyzs[3*(j-1)+1])*
               (xyzs[3*j+1]-xyzs[3*(j-1)+1]) +
               (xyzs[3*j+2]-xyzs[3*(j-1)+2])*
               (xyzs[3*j+2]-xyzs[3*(j-1)+2]));
  }
  dx = 0.0;
  for (j = 1; j < nsec; j++) {
    dx += sqrt((xyzs[3*j  ]-xyzs[3*(j-1)  ])*
               (xyzs[3*j  ]-xyzs[3*(j-1)  ]) +
               (xyzs[3*j+1]-xyzs[3*(j-1)+1])*
               (xyzs[3*j+1]-xyzs[3*(j-1)+1]) +
               (xyzs[3*j+2]-xyzs[3*(j-1)+2])*
               (xyzs[3*j+2]-xyzs[3*(j-1)+2]))/dy;
    vknot[j] += dx;
  }

  for (n = j = 1; j < nsec; j++)
    if (vknot[j] < vknot[j-1]) n++;

  if ((n == 1) && (fixed == 0)) {
    for (j = 0; j < nsec; j++) vknot[j] /= dx;
  } else {
    /* equally spaced */
    for (jj = j = 0; j < nsec-1; j++)
      if (secs[j] != secs[j+1]) jj++;
    for (k = j = 0; j < nsec-1; j++) {
      dy       = k;
      vknot[j] = dy/jj;
      if (secs[j] != secs[j+1]) k++;
    }
    vknot[nsec-1] = 1.;
  }

  /* always starts with the first section */
  nsecC0    = 1;
  secsC0[0] = secs[0];
#ifdef BLEND_SPLIT_CONSTRUCTION
  aggC0[0]  = 0;
#else
  aggC0[0]  = 1;
#endif
  vdata[0] = vdata[nsec-1] = 1;
  for (j = 1; j < nsec-1; j++) {
    if ((j < nsec-3) && (vknot[j] == vknot[j+1]) && (vknot[j] == vknot[j+2])) {
      if ((j < nsec-4) && (vknot[j] == vknot[j+3])) {
        printf("EG_blend: Multiplicity > 3 at %d\n", j);
        goto cleanup;
      }
#ifdef DEBUG
      printf("repeated vknot at j=%d, %d, and %d\n", j, j+1, j+2);
#endif
      if (nsex > 0) {
#ifdef BLEND_SPLIT_CONSTRUCTION
        aggC0[2*nsecC0-1] = j+1;
        aggC0[2*nsecC0  ] = j+2;
#else
        aggC0[2*nsecC0-1] = j+2;
        aggC0[2*nsecC0  ] = j+3;
#endif
        secsC0[nsecC0++ ] = secs[j];
        vdata[j++] = 1;
        vdata[j++] = 1;
        vdata[j  ] = 1;
      } else {
        vdata[j++] = -3;
        vdata[j++] =  1;
        vdata[j  ] = +3;
      }
    } else if (vknot[j] == vknot[j+1]) {
#ifdef DEBUG
      printf("repeated vknot at j=%d and %d\n", j, j+1);
#endif
      vdata[j++] = -2;
      vdata[j  ] = +2;
    } else {
      vdata[j  ] =  1;
    }
#ifdef DEBUG
  printf("after pass 1\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d   %lf\n", j, vdata[j], value(vknot[j]));
#endif
  }

  /* and cap of the last section */
#ifdef BLEND_SPLIT_CONSTRUCTION
  aggC0[2*nsecC0-1] = nsec;
#else
  aggC0[2*nsecC0-1] = nsec+1;
#endif
  secsC0[nsecC0++ ] = secs[nsec-1];

  /* for all multiplicity 2 knots, determine where the flat spot is */
  for (j = 1; j < nsec-1; j++)
    if ((vdata[j] == -2) && (vdata[j+1] == +2)) {
      if (j < 2) {
        /* flat on left */
        vdata[j+1] = 1;
      } else if (j > nsec-4) {
        /* flat on right */
        vdata[j  ] = 1;
      } else if ((vknot[j-1] == vknot[j-2]) && (vknot[j+3] == vknot[j+2])) {
        /* flat on both sides */
        printf("EG_blend: C1 -- Too few sections on either side of %d and %d\n",
               j, j+1);
        stat = EGADS_DEGEN;
        goto cleanup;
      } else if (vknot[j-1] == vknot[j-2]) {
        /* flat on left because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on left at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j+1] = 1;
      } else if (vknot[j+3] == vknot[j+2]) {
        /* flat on right because of another multiplicity 2 too close */
#ifdef DEBUG
        printf("flat on right at %d and %d due to close data\n", j, j+1);
#endif
        vdata[j  ] = 1;
      } else {
        /* find second derivatives on both sides */
        left = rite = 0;
        d2xdt2L = ((xyzs[3*(j  )  ]-
                    xyzs[3*(j-1)  ])/(vknot[j  ]-vknot[j-1])
                -  (xyzs[3*(j-1)  ]-
                    xyzs[3*(j-2)  ])/(vknot[j-1]-vknot[j-2]))
                / (vknot[j  ]-vknot[j-2]);
        d2ydt2L = ((xyzs[3*(j  )+1]-
                    xyzs[3*(j-1)+1])/(vknot[j  ]-vknot[j-1])
                -  (xyzs[3*(j-1)+1]-
                    xyzs[3*(j-2)+1])/(vknot[j-1]-vknot[j-2]))
                / (vknot[j  ]-vknot[j-2]);
        d2zdt2L = ((xyzs[3*(j  )+2]-
                    xyzs[3*(j-1)+2])/(vknot[j  ]-vknot[j-1])
                -  (xyzs[3*(j-1)+2]-
                    xyzs[3*(j-2)+2])/(vknot[j-1]-vknot[j-2]))
                / (vknot[j  ]-vknot[j-2]);

        d2xdt2R = ((xyzs[3*(j+3)  ]-
                    xyzs[3*(j+2)  ])/(vknot[j+3]-vknot[j+2])
                -  (xyzs[3*(j+2)  ]-
                    xyzs[3*(j+1)  ])/(vknot[j+2]-vknot[j+1]))
                / (vknot[j+3]-vknot[j+1]);
        d2ydt2R = ((xyzs[3*(j+3)+1]-
                    xyzs[3*(j+2)+1])/(vknot[j+3]-vknot[j+2])
                -  (xyzs[3*(j+2)+1]-
                    xyzs[3*(j+1)+1])/(vknot[j+2]-vknot[j+1]))
                / (vknot[j+3]-vknot[j+1]);
        d2zdt2R = ((xyzs[3*(j+3)+2]-
                    xyzs[3*(j+2)+2])/(vknot[j+3]-vknot[j+2])
                -  (xyzs[3*(j+2)+2]-
                    xyzs[3*(j+1)+2])/(vknot[j+2]-vknot[j+1]))
                / (vknot[j+3]-vknot[j+1]);

        d2sdt2L = d2xdt2L*d2xdt2L + d2ydt2L*d2ydt2L + d2zdt2L*d2zdt2L;
        d2sdt2R = d2xdt2R*d2xdt2R + d2ydt2R*d2ydt2R + d2zdt2R*d2zdt2R;
        max2    = fabs(d2sdt2L);
        if (max2 < fabs(d2sdt2R)) max2 = fabs(d2sdt2R);
        if (max2 == 0.0) {
          left++;
        } else if (fabs(d2sdt2L-d2sdt2R)/max2 < 1.e-6) {
          left++;
        } else if (d2sdt2L > d2sdt2R) {
          rite++;
        } else {
          left++;
        }

        if (rite > left) {
          /* flat on right because left curvature is smaller */
#ifdef DEBUG
          printf("flat on right %d %d\n", left, rite);
#endif
          vdata[j  ] = 1;
        } else {
          /* flat on left because right curvture is smaller  */
#ifdef DEBUG
          printf("flat on left %d %d\n", left, rite);
#endif
          vdata[j+1] = 1;
        }
      }
      j++;
    }
#ifdef DEBUG
  printf("after pass 2\n");
  for (j = 0; j < nsec; j++)
    printf("vdata[%2d]=%2d\n", j, vdata[j]);
#endif


#ifdef BLEND_SPLIT_CONSTRUCTION
  for (iC0 = 0; iC0 < nsecC0-1; iC0++) {
    nC0 = aggC0[2*iC0+1]-aggC0[2*iC0];

    /* normalize the knots based on the C0 sections */
    dy = vknot[aggC0[2*iC0+1]-1]-vknot[aggC0[2*iC0]];
    for (i = 0; i < nC0; i++) {
      vknotC0[i] = (vknot[i+aggC0[2*iC0]] - vknot[aggC0[2*iC0]])/dy;
    }

    for (j = 0; j < nstripe; j++) {
      xyzs = (T *) EG_alloc(3*(ncp[j].ncp+2)*nC0*sizeof(T));
      if (xyzs == NULL) {
        if (outLevel > 0)
          printf(" EGADS Error: Allocation for %dX%d points (EG_blend)!\n",
                 ncp[j].ncp, nC0);
        stat = EGADS_MALLOC;
        goto cleanup;
      }
      t1 = &xyzs[3* ncp[j].ncp   *nC0];
      tN = &xyzs[3*(ncp[j].ncp+1)*nC0];

      /* get the spline points between the two C0 sections for the current stripe */
      stat = EG_secSplinePoints(outLevel, nC0, &secs[aggC0[2*iC0]], j, ncp,
                                &planar, xyzs, t1, tN);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* index now based on the C0 sections */
      i = iC0;

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSplinePoints_%d.dat", j);
        EG_tecplotSplinePoints(ncp[j].ncp, nC0, xyzs, filename);
      }
#endif

      if (nC0 == 2 && iC0 > 0 && iC0 < nsecC0-2 ) {
        /* single ruled spline */
        if ((planar == 1) || (ncp[j].ncp == 3)) {
          stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                    1.e-8, surfs[i+j*isrf].header,
                                    &surfs[i+j*isrf].data);
        } else {
          stat = EG_spline2dAppr<T>(1, ncp[j].ncp, 2, xyzs, ncp[j].knots, NULL,
                                    NULL, t1, tN, NULL, NULL, NULL, NULL,
                                    1.e-8, surfs[i+j*isrf].header,
                                    &surfs[i+j*isrf].data);
        }
      } else {

        if (planar == 0) {
          sides[0] = t1;
          sides[2] = tN;
        } else {
          sides[0] = NULL;
          sides[2] = NULL;
        }

        sides[1] = iC0 == 0 ? rc1 : NULL;
        sides[3] = iC0 == nsecC0-2 ? rcN : NULL;

        /* get the BSpline surface */
        stat = EG_loft2spline(outLevel, vknotC0, &vdata[aggC0[2*iC0]],
                              sides, &ncp[j], nC0, xyzs, 1.e-8,
                              surfs[i+j*isrf].header, &surfs[i+j*isrf].data);
      }
      EG_free(xyzs); xyzs = NULL;
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Strip %d splined = %d (EG_blend)!\n", j+1, stat);
        goto cleanup;
      }

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "blendSpline_%dx%d.dat", i,j);
        EG_tecplotSpline(surfs[i+j*isrf].header, surfs[i+j*isrf].data, filename);
      }
#endif

      /* get the u = 0 curveU[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                         0, -1,
                         curvsU[i+j*icrvU].header, &curvsU[i+j*icrvU].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      if (j == nstripe-1) {
        if (uclosed == 1) {
          /* uclosed loop so first and last curves are the same */
        } else {
          /* get the u = 1 curveU[i,j+1] */
          stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                             surfs[i+j*isrf].header[2]-1, -1,
                             curvsU[i+(j+1)*icrvU].header,
                             &curvsU[i+(j+1)*icrvU].data);

          if (stat != EGADS_SUCCESS) {
            if (outLevel > 0)
              printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                     j+1, i+1, stat);
            goto cleanup;
          }
        }
      }

      /* get the v = 0 curveV[i,j] */
      stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                         -1, 0,
                         curvsV[i+j*icrvV].header, &curvsV[i+j*icrvV].data);

      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Edge %d/%d isoCurve- = %d (EG_blend)!\n",
                 j+1, i+1, stat);
        goto cleanup;
      }

      /* get the v = 1 curveV[i+1,j] on the last section */
      if ( i == nsecC0-2 ) {
        stat = EG_isoCurve(surfs[i+j*isrf].header, surfs[i+j*isrf].data,
                           -1, surfs[i+j*isrf].header[5]-1,
                           curvsV[(i+1)+j*icrvV].header,
                           &curvsV[(i+1)+j*icrvV].data);

        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Edge %d/%d isoCurve+ = %d (EG_blend)!\n",
                   j+1, i+1, stat);
          goto cleanup;
        }
      }
    }
  }

#else

  /* get the BSpline curve */
  stat = EG_spline1dFit( 1, nsec, xyzs,
                         vknot, 1.e-8, blendcurve.header,
                         &blendcurve.data );
  EG_free(xyzs); xyzs = NULL;
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: Curve spline = %d (EG_blend)!\n", stat);
    goto cleanup;
  }


  ncurvU = nsecC0-1;

  nodes  = (ego*)EG_reall( nodes, nsecC0 *sizeof(ego));
  curvsU = new egSpline<T>[ncurvU];
  if ((nodes == NULL) || (curvsU == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation (EG_blend)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  /* get nodes for C0 sections */
  for (i = 0; i < nsecC0; i++) {
    /* get the node of this section */
    stat = EG_getSecNodes( outLevel, 1, secsC0[i], 0, &nodes[i] );
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (vclosed == 1) {
    nodes[nsecC0-1] = nodes[0];
  }

  /* make curves */
  for (i = 0; i < nsecC0-1; i++) {

    /* get the sub-BSpline curve */
    stat = EG_subSpline1d(blendcurve.header, blendcurve.data,
                          aggC0[2*i], aggC0[2*i+1],
                          curvsU[i].header, &curvsU[i].data);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: subSpline1d = %d (EG_blend)!\n",
               stat);
      goto cleanup;
    }
  }
#endif

  stat = EGADS_SUCCESS;

cleanup:
  /* clean up all of our temps */
  EG_free(aggC0);
  EG_free(vknot);
  EG_free(vknotC0);
  EG_free(vknotTE);
  EG_free(vdata);
  EG_free(vdataTE);
  EG_free(xyzs);

  if (stat != EGADS_SUCCESS) {
    EG_free(secsC0);
    EG_free(nodes);
    delete [] curvsU;
  } else {
    *nsecC0_out = nsecC0;
    *secsC0_out = secsC0;
    *nodes_out = nodes;
    *curvsU_out = curvsU;
  }

  return stat;
}


template<class T>
static int
EG_ruledNodeSpline(
#ifdef EGADS_SPLINE_VELS
               const egadsSplineVels *vels,
#endif
               int nsec, const ego *secs, int vclosed,
               ego **nodes_out, egSpline<T> **curvsU_out)
{
  int    i, outLevel, stat;
  int    nnode=0, ncurvU=0;
  ego    *nodes=NULL;
  egSpline<T> *curvsU=NULL;

  *nodes_out  = NULL;
  *curvsU_out = NULL;

  outLevel = EG_outLevel(secs[0]);

  nnode  = nsec;
  ncurvU = nsec-1;

  /* make a nodes, curves, and surfaces */
  nodes = (ego*)EG_alloc( nnode *sizeof(ego));
  curvsU = new egSpline<T>[ncurvU];
  if ((nodes == NULL) || (curvsU == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation (EG_ruled)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  for (i = 0; i < nsec; i++) {
    /* get the node of this section */
    stat = EG_getSecNodes( outLevel, 1, secs[i], 0, &nodes[i] );
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  if (vclosed == 1) {
    nodes[0] = nodes[nsec-1];
  }

  /* first construct all surfaces (skipping degenerate ones) */
  for (i = 0; i < nsec-1; i++) {

    /* skip multipicty of nodes */
    stat = EG_isEquivalent(nodes[i], nodes[i+1]);
    if (stat < EGADS_SUCCESS) goto cleanup;
    if (stat == EGADS_SUCCESS) {
      nodes[i+1] = nodes[i];
      continue;
    }

    /* create the spline data */
    curvsU[i].header[0] = 0; /* Bit flag */
    curvsU[i].header[1] = 1; /* Degree */
    curvsU[i].header[2] = 2; /* nCP */
    curvsU[i].header[3] = 4; /* nKnots */

    curvsU[i].data = (T*)EG_alloc((4+3*2)*sizeof(T));
    if (curvsU[i].data == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation (EG_ruled)!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }

    curvsU[i].data[0] = 0; /* first knot points */
    curvsU[i].data[1] = 0;

    curvsU[i].data[2] = 1; /* last knot points */
    curvsU[i].data[3] = 1;

    /* evaluate the node coordinates */
#ifdef EGADS_SPLINE_VELS
    stat = EG_nodeVels(vels, outLevel, i  , secs, nodes[i  ], curvsU[i].data+4);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_nodeVels(vels, outLevel, i+1, secs, nodes[i+1], curvsU[i].data+7);
    if (stat != EGADS_SUCCESS) goto cleanup;
#else
    stat = EG_evaluat(nodes[i  ], NULL, curvsU[i].data+4);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_evaluat(nodes[i+1], NULL, curvsU[i].data+7);
    if (stat != EGADS_SUCCESS) goto cleanup;
#endif
  }

  stat = EGADS_SUCCESS;

cleanup:
  /* clean up all of our temps */
  if (stat != EGADS_SUCCESS) {
    EG_free(nodes);
    delete [] curvsU;
  } else {
    *nodes_out  = nodes;
    *curvsU_out = curvsU;
  }

  return stat;
}


static int
EG_splineNodeGeom(int nsec, const ego *secs, int vclosed,
                  ego *nodes, egSpline<double> *splcurvsU,
                  ego **edgesU_out, objStack *stack)
{
  int    i, n, outLevel, stat;
  int    ncurvU=0;
  double ts[2];
  ego    context, nds[2];
  ego    *curvsU=NULL;
  ego    *edgesU=NULL;

  *edgesU_out = NULL;

  /* Geometric entities are laid out as follows
   *
   *  sec[i+1]  nodes[i+1]
   *              |
   *              |
   *              |
   *  curvsU[i]   |   v
   *              |   ^
   *              |   |
   *              |   +--> u
   *              |
   *  sec[i]    nodes[i]
   *
   */

  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL) return EGADS_NOTCNTX;

  ncurvU = nsec-1;

  curvsU = (ego *) EG_alloc(ncurvU*sizeof(ego));
  edgesU = (ego *) EG_alloc(ncurvU*sizeof(ego));
  if ((curvsU == NULL) ||
      (edgesU == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: Allocation (EG_splineNodeGeom)!\n");
    stat = EGADS_MALLOC;
    goto cleanup;
  }

  for (i = 0; i < ncurvU; i++) {
    curvsU[i] = NULL;
    edgesU[i] = NULL;
  }

  ts[0] = 0.0;
  ts[1] = 1.0;

  /* construct ego geometry from spline information */
  for (i = 0; i < nsec-1; i++) {

    /*-------------------------*/
    /*      curvsU[i]          */
    /*-------------------------*/

    /* make the curve that spans the sections */
    stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                           splcurvsU[i].header,
                           splcurvsU[i].data, &curvsU[i]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Curve %d EG_makeGeometry = %d (EG_splineNodeGeom)!\n",
               i+1, stat);
      goto cleanup;
    }
    stat = EG_stackPush(stack, curvsU[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    n = 2;
    nds[0] = nodes[i  ];
    nds[1] = nodes[i+1];

    /* make the edge on the u-constant curve */
    stat = EG_makeTopology(context, curvsU[i], EDGE, TWONODE, ts, n,
                           nds, NULL, &edgesU[i]);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Edge %d EG_makeTopology = %d (EG_splineNodeGeom)!\n",
               i+1, stat);
      goto cleanup;
    }
    stat = EG_stackPush(stack, edgesU[i]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    ncurvU++;
  }

  /* set outputs */
  *edgesU_out = edgesU; edgesU = NULL;

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_splineNodeGeom)!\n", stat);
  }

  /* clean up all of our temps */
  EG_free(edgesU);
  EG_free(curvsU);

  return stat;
}


static int
EG_splineNodeGeom_dot(ego body, int nsec, const ego *secs,
                      ego *nodes, egSpline< SurrealS<1> > *splcurvsU)
{
  int    i, outLevel, stat, oclass, mtype;
  int    nchldrn, nedge, icurvU;
  int    *senses;
  SurrealS<1> ts[2] = {0, 1};
  double data[4];
  ego    curv, loop, ref;
  ego    *chldrn, *nds, *edges;

  outLevel = EG_outLevel(body);

  /* clear any old sensitivities */
  stat = EG_setGeometry_dot(body, BODY, 0, NULL, NULL, NULL);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: getTopology = %d (EG_splineNodeGeom_dot)!\n", stat);
    goto cleanup;
  }

  /* get the loop */
  stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                        &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;
  if (mtype != WIREBODY) {
    if (outLevel > 0)
      printf(" EGADS Error: body is not a WIREBODY (EG_splineNodeGeom_dot)!\n");
    goto cleanup;
  }

  /* get the edges from the loop */
  loop = chldrn[0];
  stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &nedge, &edges,
                        &senses);
  if (stat != EGADS_SUCCESS) goto cleanup;


  icurvU = 0;
  for (i = 0; i < nsec-1; i++) {

    if (splcurvsU[i].data == NULL) continue;

    /* get the curv and nodes for the u = 0 edge */
    stat = EG_getTopology(edges[icurvU], &curv, &oclass, &mtype, data,
                          &nchldrn, &nds, &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity for the curve */
    stat = EG_setGeometry_dot(curv, CURVE, BSPLINE,
                              splcurvsU[icurvU].header,
                              splcurvsU[icurvU].data);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* set the sensitivity for the t-range of the edge */
    stat = EG_setRange_dot(edges[icurvU], EDGE, ts);
    if (stat != EGADS_SUCCESS) goto cleanup;

    /* copy over node sensitivities */
    stat = EG_copyGeometry_dot(nodes[i], NULL, NULL, nds[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (icurvU == nedge-1) {
      stat = EG_copyGeometry_dot(nodes[i+1], NULL, NULL, nds[1]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    icurvU++;
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: status = %d (EG_splineNodeGeom_dot)!\n", stat);
  }

  return stat;
}


static int
EG_checkSections(int nsec, const ego *secs, const char *func,
                 int &nstripe, int &ncap, int &bodyType, int &planar,
                 int &uclosed, int &vclosed, int &allNodes)
{
  int    i, j, k, jj, n, stat, oclass, mtype, outLevel;
  int    *senses=NULL;
  double data[18];
  ego    context, loop, ref, geom, loopBeg=NULL, loopEnd=NULL;
  ego    *chldrn, *edges;

  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);

  /* check for excess multiplicity */
  for (n = i = 0; i < nsec-1; i++) {
    if (EG_isEquivalent(secs[i], secs[i+1]) != EGADS_SUCCESS) {
      n = 0;
    }
    n++;
    if (n > 3) {
      printf(" EGADS Error:: Section %d has multiplicity > 3 (%s)\n", i, func);
      return EGADS_TOPOCNT;
    }
    if (n > 2 && (i == 0 || i == nsec-2)) {
      printf(" EGADS Error:: Section %d has multiplicity > 2 (%s)\n", i, func);
      return EGADS_TOPOCNT;
    }
  }

  /* look at the input and check to see if OK */
  uclosed = 1;
  bodyType = SOLIDBODY;
  ncap = nstripe = vclosed = planar = 0;
  allNodes = 1;

  /* first check if all sections are Nodes */
  for (i = 0; i < nsec; i++) {
    if (secs[i] == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: NULL Section Object %d (%s)!\n", i+1, func);
      return EGADS_NULLOBJ;
    }
    if (secs[i]->magicnumber != MAGIC) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is not an EGO (%s)!\n", i+1, func);
      return EGADS_NOTOBJ;
    }
    if (secs[i]->blind == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d has no data (%s)!\n", i+1, func);
      return EGADS_NODATA;
    }
    if (EG_context(secs[i]) != context) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Context Mismatch (%s)!\n", i+1, func);
      return EGADS_MIXCNTX;
    }
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d getTopology = %d (%s)!\n",
               i+1, stat, func);
      return stat;
    }

    loop = NULL;
    if (oclass == FACE) {
      allNodes = 0;
    } else if (oclass == LOOP) {
      loop = secs[i];
    } else if (oclass == BODY) {
      if (secs[i]->mtype == WIREBODY) {
        loop = chldrn[0];
        if ((i == 0) || (i == nsec-1)) bodyType = SHEETBODY;
      } else {
        allNodes = 0;
      }
    } else if (oclass != NODE){
      allNodes = 0;
    }

    if (loop == NULL) continue;
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &jj, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Loop getTopology = %d (%s)!\n",
               i+1, stat, func);
      return stat;
    }
    if (mtype != DEGENERATE) allNodes = 0;
  }

  if (allNodes == 1) {
    /* check if closed in v-direction */
    if (EG_isEquivalent(secs[0], secs[nsec-1]) == EGADS_SUCCESS) {
      vclosed = 1;
      ncap = 0;
    }
    bodyType = WIREBODY;
    return EGADS_SUCCESS;
  }

  /* at least one section has edges, check for consistency */
  for (i = 0; i < nsec; i++) {
    stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d getTopology = %d (%s)!\n",
               i+1, stat, func);
      return stat;
    }
    if (oclass == NODE) {
      if ((i != 0) && (i != nsec-1)) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Node and not Bound (%s)!\n",
                 i+1, func);
        return EGADS_NOTTOPO;
      }
      loop = NULL;
    } else if (oclass == FACE) {
      if (n != 1) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Face with %d Loops (%s)!\n",
                 i+1, n, func);
        return EGADS_TOPOERR;
      }
      loop = chldrn[0];
      if ((i == 0) || (i == nsec-1)) ncap++;
    } else if (oclass == BODY) {
      if (secs[i]->mtype == WIREBODY) {
        loop = chldrn[0];
        if ((i == 0) || (i == nsec-1)) bodyType = SHEETBODY;
      } else if (secs[i]->mtype == FACEBODY) {

        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &n, &chldrn,
                              &senses);
        if (n != 1 || stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: Section %d is FaceBody with %d Loops (%s)!\n",
                   i+1, n, func);
          return EGADS_TOPOERR;
        }
        loop = chldrn[0];
        if ((i == 0) || (i == nsec-1)) ncap++;
      } else {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d is Not a WireBody or FaceBody (%s)!\n",
                 i+1, func);
        return EGADS_NOTTOPO;
      }
    } else if (oclass == LOOP) {
      loop = secs[i];
      if ((i == 0) || (i == nsec-1)) bodyType = SHEETBODY;
    } else {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d is Not a Loop (%s)!\n", i+1, func);
      return EGADS_NOTTOPO;
    }

    if (loop == NULL) continue;
    if (i == 0     ) loopBeg = loop;
    if (i == nsec-1) loopEnd = loop;
    if (EG_isPlanar(loop) != EGADS_SUCCESS) planar = 1;
    stat = EG_getTopology(loop, &ref, &oclass, &mtype, data, &jj, &edges,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: Section %d Loop getTopology = %d (%s)!\n",
               i+1, stat, func);
      return stat;
    }
    if (mtype == OPEN) { uclosed = 0; bodyType = SHEETBODY; }
    for (n = j = 0; j < jj; j++) {
      stat = EG_getTopology(edges[j], &geom, &oclass, &mtype, data, &k, &chldrn,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d Edge %d getTopo = %d (%s)!\n",
                 i+1, j+1, stat, func);
        return stat;
      }
      if (mtype != DEGENERATE) n++;
    }
    if (nstripe == 0) {
      nstripe = n;
    } else {
      if (n != nstripe) {
        if (outLevel > 0)
          printf(" EGADS Error: Section %d has %d Edges -- prev = %d (%s)!\n",
                 i+1, n, nstripe, func);
        return EGADS_TOPOERR;
      }
    }
  }
  if (nstripe == 0) {
    if (outLevel > 0)
      printf(" EGADS Error: Zero Edges found (%s)!\n", func);
    return EGADS_TOPOERR;
  }
  if (uclosed == 0) ncap = 0;

  /* check if closed in v-direction */
  if ((loopBeg != NULL) && (loopEnd != NULL) &&
      (EG_isEquivalent(loopBeg, loopEnd) == EGADS_SUCCESS)) {
    vclosed = 1;
    ncap = 0;
  }
  if ((uclosed == 1) && (vclosed == 1)) bodyType = SOLIDBODY;

  return EGADS_SUCCESS;
}


static int
EG_uniqueSections(int nsec_in, const ego *secs_in,
                  int &nsec, ego **secs)
{
  int i;
  *secs = (ego*)EG_alloc(nsec_in*sizeof(ego));
  if (*secs == NULL) return EGADS_MALLOC;

  /* get the unique sections */
  nsec = 0;
  for (i = 0; i < nsec_in-1; i++) {
    if (EG_isEquivalent(secs_in[i], secs_in[i+1]) == EGADS_SUCCESS) continue;
    (*secs)[nsec++] = secs_in[i];
  }
  if (EG_isEquivalent(secs_in[nsec_in-2], secs_in[nsec_in-1]) != EGADS_SUCCESS)
    (*secs)[nsec++] = secs_in[nsec_in-1];

  return EGADS_SUCCESS;
}

} // namespace


extern "C" int
EG_blend(int nsex, const ego *secs, /*@null@*/ double *rc1,
                                    /*@null@*/ double *rcN, ego *result)
{
  typedef double T;
  int    i, j, k, n, m, mm, outLevel, stat, nstripe, ncap, cap, ntip;
  int    oclass, mtype, attrint, te, deg, iknotc, header[4], vknotN, inode=0;
  int    bodyType, uclosed, vclosed, allNodes, tip, nface=0, icrvU, icrvV, ncurvV=0, isrf=0;
  int    *senses, planar, *csens=NULL, tsens[8], lsens[1] = {SFORWARD};
  int    nsec, nsecC0, nedge, begRC, endRC, tips[2], tipsec[2], bstrp[2];
  double data[MAX(18,NTIPTE+2+2*NTIPTE)], ts[2], mknot, vknot[2], *sknotv;
  double knotc[2*NTIPTE-1], rc1S[8], rcNS[8], *rc1s=NULL, *rcNs=NULL, len;
  ego    context, surf, loop, ref, shell, body;
  ego    tedges[16]={NULL}, pcurvs[5], tnodes[3], tipcurv, tipnodes[2], tipedges[4];
  ego    *chldrn=NULL, *secsC0=NULL, *nodes=NULL, *edgesU=NULL;
  ego    *edgesV=NULL, *faces=NULL, *cedges=NULL;
  objStack    stack;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;
  egSpline<T> tipsurfs[2], tipcurvs[2], curv;

  nsec    = nsex;
  if (nsec < 0) nsec = -nsec;
  *result = NULL;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL)               return EGADS_NOTCNTX;

#if !defined(WIN32) && defined(BLEND_FPE)
  feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
#endif

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec, secs, "EG_blend",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;


  /* special case for all nodes */
  if (allNodes == 1) {

#ifdef EGADS_SPLINE_VELS
    stat = EG_blendNodeSpline(NULL, nsex, secs, vclosed,
                              &nsecC0, &secsC0, &nodes, &splcurvsU);
#else
    stat = EG_blendNodeSpline(nsex, secs, vclosed,
                              &nsecC0, &secsC0, &nodes, &splcurvsU);
#endif
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: blendNodeSpline = %d (EG_blend)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom(nsecC0, secsC0, vclosed, nodes, splcurvsU,
                             &edgesU, &stack);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom = %d (EG_blend)!\n", stat);
      goto cleanup;
    }

    senses = (int*)EG_alloc((nsecC0-1)*sizeof(int));
    if (senses == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation (EG_blend)!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nsecC0-1; i++) senses[i] = SFORWARD;

    stat = EG_makeTopology(context, NULL, LOOP, vclosed == 1 ? CLOSED : OPEN, NULL,
                           nsecC0-1, edgesU, senses, &loop);
    EG_free(senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo Loop = %d (EG_blend)!\n", stat);
      goto cleanup;
    }
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    body = NULL;
    stat = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1,
                           &loop, NULL, &body);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo WireBody = %d (EG_blend)!\n", stat);
      if (body != NULL) EG_deleteObject(body);
      goto cleanup;
    }

    stat = EGADS_SUCCESS;
    *result = body;
    goto cleanup;
  }

  /* sections with edges */

  /* check end conditions */
  begRC = endRC  = -1;
  if (rc1 != NULL) {
    i = 0;
    if (secs[i]->oclass == NODE) begRC = 1; /* nose treatment */
    else if (rc1[0] != 0       ) begRC = 2; /* tangencey */
  }
  if (rcN != NULL) {
    i = nsec-1;
    if (secs[i]->oclass == NODE) endRC = 1; /* nose treatment */
    else if (rcN[0] != 0       ) endRC = 2; /* tangencey */
  }
  if (((begRC == 1) || (endRC == 1)) && (nsec <= 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 2 for Nose Treatment (EG_blend)!\n");
    return EGADS_GEOMERR;
  }
  if ((begRC == 1) && (endRC == 1) && (nsec <= 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 3 for 2Nose Treatment (EG_blend)!\n");
    return EGADS_GEOMERR;
  }

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec, secs, "EG_blend_secs.egads");
#endif

  /* extract nose treatment */
  if (begRC == 1) {
    for (i = 0; i < 8; i++) {
      rc1S[i] = rc1[i];
    }

    /* normalize */
    len = sqrt(rc1S[1]*rc1S[1] + rc1S[2]*rc1S[2] + rc1S[3]*rc1S[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first south axis (EG_blend)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[1] /= len;
    rc1S[2] /= len;
    rc1S[3] /= len;

    len = sqrt(rc1S[5]*rc1S[5] + rc1S[6]*rc1S[6] + rc1S[7]*rc1S[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second south axis (EG_blend)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[5] /= len;
    rc1S[6] /= len;
    rc1S[7] /= len;

    rc1s = rc1S;
  }
  if (endRC == 1) {
    for (i = 0; i < 8; i++) {
      rcNS[i] = rcN[i];
    }

    /* normalize */
    len = sqrt(rcNS[1]*rcNS[1] + rcNS[2]*rcNS[2] + rcNS[3]*rcNS[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first north axis (EG_blend)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[1] /= len;
    rcNS[2] /= len;
    rcNS[3] /= len;

    len = sqrt(rcNS[5]*rcNS[5] + rcNS[6]*rcNS[6] + rcNS[7]*rcNS[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second north axis (EG_blend)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[5] /= len;
    rcNS[6] /= len;
    rcNS[7] /= len;

    rcNs = rcNS;
  }

  /* extract tangency treatment */
  if (begRC == 2) {
    for (i = 0; i < 4; i++) {
      rc1S[i] = rc1[i];
    }
    rc1s = rc1S;
  }
  if (endRC == 2) {
    for (i = 0; i < 4; i++) {
      rcNS[i] = rcN[i];
    }
    rcNs = rcNS;
  }

  /* extract tip treatment */
  if ((secs[     0]->oclass == FACE || secs[     0]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rc1S[i] = rc1[i];
    }
    rc1s = rc1S;
  }
  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rcNS[i] = rcN[i];
    }
    rcNs = rcNS;
  }

#ifdef EGADS_SPLINE_VELS
  stat = EG_blendSpline<T>(NULL, nsex, secs, rc1s, rcNs, begRC, endRC,
                           nstripe, uclosed, vclosed, planar, &nsecC0, &secsC0,
                           &nodes, &splcurvsU, &splcurvsV, &splsurfs,
                           tip, te, tipsurfs, tipcurvs);
#else
  stat = EG_blendSpline<T>(nsex, secs, rc1s, rcNs, begRC, endRC,
                           nstripe, uclosed, vclosed, planar, &nsecC0, &secsC0,
                           &nodes, &splcurvsU, &splcurvsV, &splsurfs,
                           tip, te, tipsurfs, tipcurvs);
#endif
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: blendSpline = %d (EG_blend)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom(nsecC0, secsC0, nstripe, uclosed, vclosed, ncap, te, tip,
                       nodes, splcurvsU, splcurvsV, splsurfs,
                       &ncurvV, &edgesU, &edgesV, &nface, &faces, &stack);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: splineGeom = %d (EG_blend)!\n", stat);
    goto cleanup;
  }

  inode = nsecC0;
  icrvU = nsecC0-1;
  icrvV = nsecC0;
  isrf  = nsecC0-1;

  /* add cap faces if necessary */
  if (ncap > 0) {

    cedges = (ego *) EG_alloc(2*nstripe*sizeof(ego));
    csens  = (int *) EG_alloc(  nstripe*sizeof(int));
    if ((cedges == NULL) || (csens == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocating edges for loop!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for (j = 0; j < 2*nstripe; j++) cedges[j] = NULL;
    for (j = 0; j <   nstripe; j++) csens[j]  = SFORWARD;

    cap = 0;
    for (k = 0; k < 2; k++) {
      if ((tip & 1) && (k == 0)) { cap++; continue; }
      if ((tip & 2) && (k == 1)) { cap++; continue; }
      i = (k == 0) ? 0 : nsecC0-1;
      stat = EG_getTopology(secsC0[i], &ref, &oclass, &mtype, data,
                            &n, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (oclass == BODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data,
                              &n, &chldrn, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (oclass != FACE) continue;

      /* grab the spline edges on the cap */
      for (j = 0; j < nstripe; j++) cedges[j] = edgesV[i+j*icrvV];

      if (ref->mtype == PLANE) {

        /* make the loop (no need for surface with planar) */
        stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, nstripe,
                               cedges, csens, &loop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Cap %d Planar Loop = %d (EG_blend)!\n",
                   cap, stat);
          goto cleanup;
        }

      } else {

        /* make PCurves for the loop */
        for (j = 0; j < nstripe; j++) {
          /* construct pcurves on the surface */
          double eps = 1.e-7;
          EG_tolerance(cedges[j], &eps);
          if (eps < 1.e-7) eps = 1.e-7;
          for (int jj = 0; jj < 2; jj++) {
            stat = EG_otherCurve(ref, cedges[j], eps, &cedges[j+nstripe]);
            if (stat != EGADS_SUCCESS) {
              if (jj == 0) {
                eps *= 100.0;
                continue;
              }
              if (outLevel > 0)
                printf(" EGADS Error: %d otherCurve First = %d (EG_blend)!\n",
                       j, stat);
              goto cleanup;
            }
          }
          stat = EG_stackPush(&stack, cedges[j+nstripe]);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        /* make the loop with the original surface */
        stat = EG_makeTopology(context, ref, LOOP, CLOSED, NULL, nstripe,
                               cedges, csens, &loop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Cap %d Loop = %d (EG_blend)!\n",
                   cap, stat);
          goto cleanup;
        }
      }
      stat = EG_stackPush(&stack, loop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the face */
      stat = EG_makeTopology(context, ref, FACE, mtype,
                             NULL, 1, &loop, lsens, &faces[nface-ncap+cap]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo Cap %d Face = %d (EG_blend)!\n",
                 cap, stat);
        goto cleanup;
      }
      stat = EG_stackPush(&stack, faces[nface-ncap+cap]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      EG_attributeDup(secsC0[i], faces[nface-ncap+cap]);

      attrint = -(k+3);
      stat    = EG_attributeAdd(faces[nface-ncap+cap], ".blendStrip", ATTRINT,
                                1, &attrint, NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
                 attrint, stat);
      cap++;
    }
  }

  /* extract special wing tip treatment */
  cap  = 0;
  ntip = 0;
  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    bstrp[0]  = -1;
    tips[0]   = nface-ncap+cap;
    tipsec[0] = 0;
    cap++;
    ntip++;
  }

  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) &&
      (planar == 0) && (ntip == 0)) cap++;

  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    bstrp[1]  = -2;
    tips[1]   = nface-ncap+cap;
    tipsec[1] = nsecC0-1;
    cap++;
    ntip++;
  }

  /* modify cap Faces with wingTips */
  if ((ntip != 0) && (cedges != NULL)) {

    /* construct the pcurves for the wing tip surfaces */
    tsens[0] = SFORWARD;
    tsens[1] = SFORWARD;
    tsens[2] = SFORWARD;
    tsens[3] = SFORWARD;
    tsens[4] = SREVERSE;
    tsens[5] = SREVERSE;

    data[0] =  0.0; data[1] =  0.0;
    data[2] =  0.0; data[3] =  1.0; /* u == 0 */
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, pcurvs[0]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    data[0] =  0.0; data[1] =  1.0;
    data[2] =  1.0; data[3] =  0.0; /* v == 1 */
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, pcurvs[1]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    data[0] =  1.0; data[1] =  1.0;
    data[2] =  0.0; data[3] = -1.0; /* u == 1 */
    stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;
    stat = EG_stackPush(&stack, pcurvs[2]);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (nstripe == 2) {
      data[0] =  1.0; data[1] =  0.0;
      data[2] = -1.0; data[3] =  0.0; /* v == 0 */
      stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data, &pcurvs[3]);
      if (stat != EGADS_SUCCESS) goto cleanup;
      stat = EG_stackPush(&stack, pcurvs[3]);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    cap = 0;
    for (k = 0; k < 2; k++) {
      if (!(tip & 1) && (k == 0)) { cap++; continue; }
      if (!(tip & 2) && (k == 1)) { cap++; continue; }
      i = (k == 0) ? 0 : nsecC0-1;

      for (j = 0; j < nstripe; j++) {
        cedges[j] = edgesV[tipsec[k]+j*icrvV];
      }

#ifdef SPLINE_TECPLOT_DEBUG
      {
        char filename[42];
        sprintf(filename, "wingtip_%d.dat", k+1);
        EG_tecplotSpline(tipsurfs[k].header, tipsurfs[k].data, filename);
      }
#endif

      stat = EG_makeGeometry(context, SURFACE, BSPLINE, NULL, tipsurfs[k].header,
                             tipsurfs[k].data, &surf);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Tip Surface %d EG_makeGeometry = %d (EG_blend)!\n",
                 k, stat);
        goto cleanup;
      }
      stat = EG_stackPush(&stack, surf);
      if (stat != EGADS_SUCCESS) goto cleanup;

      ts[0] = 0.0; ts[1] = 1.0;
      /* set the pcurves */
      if (nstripe == 2) {
        nedge = 4;

        tedges[0] = cedges[0];

        stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1,
                               &nodes[tipsec[k]+1*inode], NULL, &tedges[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tedges[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        tedges[2] = cedges[1];

        stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1,
                               &nodes[tipsec[k]+2*inode], NULL, &tedges[3]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tedges[3]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        tedges[4] = pcurvs[0];
        tedges[5] = pcurvs[1];
        tedges[6] = pcurvs[2];
        tedges[7] = pcurvs[3];
      } else {
        nedge = 5;

        /* grab edges depending on which one is the te edge */
        int js[3][2] = {{1,2}, {2,0}, {0,1}};
        n = 0; j = js[te][0];
        tedges[n]   = cedges[j];
        tedges[5+n] = pcurvs[n];

        n = 1; j = js[te][0];
        stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1,
                               &nodes[tipsec[k]+(j+1)*inode], NULL,
                               &tedges[n]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tedges[n]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        tedges[5+n] = pcurvs[n];

        n = 2; j = js[te][1];
        tedges[n]   = cedges[j];
        tedges[5+n] = pcurvs[n];

        /* get the midpoint knot value for the TE tip curve */
        mknot = tipcurvs[k].data[tipcurvs[k].header[3]/2];

        /* make the tip node */
        stat = EG_spline1dEval(tipcurvs[k].header, tipcurvs[k].data, mknot,
                               data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_makeTopology(context, NULL, NODE, 0, data, 0, NULL, NULL,
                               &tipnodes[k]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tipnodes[k]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        tnodes[0] = nodes[tipsec[k]+(te+1)*inode];
        tnodes[1] = tipnodes[k];
        tnodes[2] = nodes[tipsec[k]+ te   *inode];

        /* make the tip curve and edges */
        stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL,
                               tipcurvs[k].header, tipcurvs[k].data, &tipcurv);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tipcurv);
        if (stat != EGADS_SUCCESS) goto cleanup;

        ts[0] = 0.0;
        ts[1] = mknot;
        stat = EG_makeTopology(context, tipcurv, EDGE, TWONODE, ts, 2, tnodes,
                               NULL, &tipedges[2*k]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tipedges[2*k]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        header[0] = 0;
        header[1] = 1;
        header[2] = 2;
        header[3] = 4;

        data[0] = data[1] = 0.0;
        data[2] = data[3] = mknot; /* v == 0 */
        data[4] = 0.0;   data[5] = 0.0;
        data[6] = mknot; data[7] = 0.0;
        stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                               &pcurvs[3]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, pcurvs[3]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        tedges[4]   = tipedges[2*k];
        tedges[5+4] = pcurvs[3];

        ts[0] = mknot;
        ts[1] = 1.0;
        stat = EG_makeTopology(context, tipcurv, EDGE, TWONODE, ts, 2,
                               &tnodes[1], NULL, &tipedges[2*k+1]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tipedges[2*k+1]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        data[0] = data[1] = mknot;
        data[2] = data[3] = 1.0; /* v == 0 */
        data[4] = mknot; data[5] = 0.0;
        data[6] = 1.0;   data[7] = 0.0;
        stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                               &pcurvs[4]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, pcurvs[4]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        tedges[3]   = tipedges[2*k+1];
        tedges[5+3] = pcurvs[4];
      }

      /* make the loop with the new surface */
      stat = EG_makeTopology(context, surf, LOOP, CLOSED, NULL, nedge,
                             tedges, tsens, &loop);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo Tip %d Loop = %d (EG_blend)!\n",
                 k, stat);
        goto cleanup;
      }
      stat = EG_stackPush(&stack, loop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the tip treatment face
       * SERVERSE here because the loop is constructed in reverse in EG_wingTipSpline
       */
      stat = EG_makeTopology(context, surf, FACE, SREVERSE,
                             NULL, 1, &loop, lsens, &faces[tips[k]]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo Tip %d Face = %d (EG_blend)!\n",
                 k, stat);
        goto cleanup;
      }
      stat = EG_stackPush(&stack, faces[tips[k]]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      EG_attributeDup(secsC0[i], faces[tips[k]]);

      stat = EG_attributeAdd(faces[tips[k]], ".blendStrip", ATTRINT, 1,
                             &bstrp[k], NULL, NULL);
      if (stat != EGADS_SUCCESS)
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
                 attrint, stat);
    }

    /* construct the trailing edge face */
    if (nstripe == 3) {
      j = te;
      for (i = 0; i < nsecC0-1; i += nsecC0-2) {

        /* make the surface */
        stat = EG_makeGeometry(context, SURFACE, BSPLINE, NULL,
                               splsurfs[i+j*isrf].header,
                               splsurfs[i+j*isrf].data, &surf);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, surf);
        if (stat != EGADS_SUCCESS) goto cleanup;

        curv.clear();
        stat = EG_isoCurve(splsurfs[i+j*isrf].header, splsurfs[i+j*isrf].data,
                           0, -1,
                           curv.header, &curv.data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* get the vknots of the surface, and the offset for the 2nd tip */
        sknotv = &curv.data[curv.header[1]];
        vknotN = curv.header[3] - 2*curv.header[1] - NTIPTE;

        vknot[0] = 0.;
        vknot[1] = 1.;
        if ((i == 0       ) && (tip & 1))
          vknot[0] = curv.data[curv.header[1] + NTIPTE];
        if ((i == nsecC0-2) && (tip & 2))
          vknot[1] = curv.data[curv.header[1] + vknotN];

        nedge = 2;
        if ((i == 0       ) && (tip & 1)) nedge += 3; else nedge += 1;
        if ((i == nsecC0-2) && (tip & 2)) nedge += 3; else nedge += 1;

        k = 0;
        tsens[k] = SREVERSE;
        tedges[k] = edgesU[i+j*icrvU];

        ts[0] = 0.0;
        ts[1] = 1.0;

        header[0] = 0;
        header[1] = 1;
        header[2] = 2;
        header[3] = 4;

        data[0] = data[1] = 0.0;
        data[2] = data[3] = 1.0; /* u == 0 */
        data[4] = 0.0; data[5] = vknot[0];
        data[6] = 0.0; data[7] = vknot[1];
        stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                               &tedges[k+nedge]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tedges[k+nedge]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        k++;

        if ((i == 0) && (tip & 1)) {
          /* surface next to tip treatment at v == 0 */
          tsens[k+0] = SREVERSE;
          tsens[k+1] = SFORWARD;
          tsens[k+2] = SREVERSE;

          tedges[k+0] = tipedges[1];
//        tedges[k+1] DEGENERATE
          tedges[k+2] = tipedges[0];

          stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1,
                                 &tipnodes[0], NULL, &tedges[k+1]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+1]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* get the midpoint knot value for the TE tip curve */
          mknot = tipcurvs[0].data[tipcurvs[0].header[3]/2];

          deg    = tipcurvs[0].header[1];
          iknotc = tipcurvs[0].header[3]-2*deg;
          n = 0;
          knotc[n++] = tipcurvs[0].data[deg];
          for (m = 0; m < iknotc-1; m++)
            for (mm = 1; mm <= NSUBTE+1; mm++)
              knotc[n++] = tipcurvs[0].data[m+deg]*(1.-mm/double(NSUBTE+1)) +
                           tipcurvs[0].data[m+1+deg]*(mm/double(NSUBTE+1));

          header[0] = 0;
          header[1] = 1;
          header[2] = NTIPTE;
          header[3] = NTIPTE+2;

          /* u == 0 */
          data[0       ] = mknot;
          data[NTIPTE+1] = 1.0;
          for (n = 0; n < NTIPTE; n++) {
            data[1+n] = knotc[NTIPTE-1+n];      /* sampling of the tip curve to
                                                   make the TE surface */
            data[NTIPTE+2 + 2*n  ] = 0.;
            data[NTIPTE+2 + 2*n+1] = sknotv[n]; /* knots from the surface */
          }
          stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                                 &tedges[k+0+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+0+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* degenerate EDGE */
          data[0] = 0.0; data[1] =  0.0;
          data[2] = 1.0; data[3] =  0.0; /* v == 0 */
          stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                                 &tedges[k+1+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+1+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* u == 1 */
          data[0       ] = 0.0;
          data[NTIPTE+1] = mknot;
          for (n = 0; n < NTIPTE; n++) {
            data[1+n] = knotc[n];     /* sampling of the tip curve to make the
                                         TE surface (reversed) */
            data[NTIPTE+2 + 2*n  ] = 1.;
            data[NTIPTE+2 + 2*n+1] = sknotv[NTIPTE-1-n]; /* knots from the
                                                          surface (reversed) */
          }
          stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                                 &tedges[k+2+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+2+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          k += 3;
        } else {

          /* no tip treatment at v == 0 */
          tsens[k]  = SFORWARD;
          tedges[k] = edgesV[i+j*icrvV];

          data[0] = 0.; data[1] = 0.; /* v == 0 */
          data[2] = 1.; data[3] = 0.;
          stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                                 &tedges[k+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          k += 1;
        }

        tsens[k]  = SFORWARD;
        tedges[k] = edgesU[i+(j+1)*icrvU];

        header[0] = 0;
        header[1] = 1;
        header[2] = 2;
        header[3] = 4;

        data[0] = data[1] = 0.0;
        data[2] = data[3] = 1.0; /* u == 1 */
        data[4] = 1.; data[5] = vknot[0];
        data[6] = 1.; data[7] = vknot[1];
        stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                               &tedges[k+nedge]);
        if (stat != EGADS_SUCCESS) goto cleanup;
        stat = EG_stackPush(&stack, tedges[k+nedge]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        k += 1;

        if ((i == nsecC0-2) && (tip & 2)) {

          tsens[k+0] = SFORWARD;
          tsens[k+1] = SFORWARD;
          tsens[k+2] = SFORWARD;

          tedges[k+0] = tipedges[2];
//        tedges[k+1] DEGENERATE
          tedges[k+2] = tipedges[3];

          stat = EG_makeTopology(context, NULL, EDGE, DEGENERATE, ts, 1,
                                 &tipnodes[1], NULL, &tedges[k+1]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+1]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* get the midpoint knot value for the TE tip curve */
          mknot = tipcurvs[1].data[tipcurvs[1].header[3]/2];

          deg    = tipcurvs[1].header[1];
          iknotc = tipcurvs[1].header[3]-2*deg;
          n = 0;
          knotc[n++] = tipcurvs[1].data[deg];
          for (m = 0; m < iknotc-1; m++)
            for (mm = 1; mm <= NSUBTE+1; mm++)
              knotc[n++] = tipcurvs[1].data[m+deg]*(1.-mm/double(NSUBTE+1)) +
                           tipcurvs[1].data[m+1+deg]*(mm/double(NSUBTE+1));

          header[0] = 0;
          header[1] = 1;
          header[2] = NTIPTE;
          header[3] = NTIPTE+2;

          /* u == 1 */
          data[0       ] = 0.0;
          data[NTIPTE+1] = mknot;
          for (n = 0; n < NTIPTE; n++) {
            data[1+n] = knotc[n];         /* sampling of the tip curve to make
                                             the TE surface */
            data[NTIPTE+2 + 2*n  ] = 1.;
            data[NTIPTE+2 + 2*n+1] = sknotv[n+vknotN]; /* knots from the surface */
          }
          stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                                 &tedges[k+0+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+0+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* degenerate EDGE */
          data[0] =  1.0; data[1] =  1.0;
          data[2] = -1.0; data[3] =  0.0; /* v == 1 */
          stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                                 &tedges[k+1+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+1+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* u == 0 */
          data[0       ] = mknot;
          data[NTIPTE+1] = 1.0;
          for (n = 0; n < NTIPTE; n++) {
            data[1+n] = knotc[NTIPTE-1+n]; /* sampling of the tip curve to make
                                              the TE surface (reversed) */
            data[NTIPTE+2 + 2*n  ] = 0.0;
            /* knots from the surface (reversed) */
            data[NTIPTE+2 + 2*n+1] = sknotv[NTIPTE-1-n+vknotN];

          }
          stat = EG_makeGeometry(context, PCURVE, BSPLINE, NULL, header, data,
                                 &tedges[k+2+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+2+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          k += 3;
        } else {

          /* no tip treatment at v == 1 */
          tsens[k]  = SREVERSE;
          tedges[k] = edgesV[i+1+j*icrvV];

          data[0] = 0.; data[1] = 1.; /* v == 1 */
          data[2] = 1.; data[3] = 0.;
          stat = EG_makeGeometry(context, PCURVE, LINE, NULL, NULL, data,
                                 &tedges[k+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;
          stat = EG_stackPush(&stack, tedges[k+nedge]);
          if (stat != EGADS_SUCCESS) goto cleanup;

          k += 1;
        }

        /* make the loop with the surface */
        stat = EG_makeTopology(context, surf, LOOP, CLOSED, NULL, nedge,
                               tedges, tsens, &loop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Tip %d Loop = %d (EG_blend)!\n",
                   k, stat);
          goto cleanup;
        }
        stat = EG_stackPush(&stack, loop);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* make the trailing edge face */
        stat = EG_makeTopology(context, surf, FACE, SFORWARD,
                               NULL, 1, &loop, lsens, &faces[i+j*isrf]);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Tip %d Face = %d (EG_blend)!\n",
                   k, stat);
          goto cleanup;
        }
        stat = EG_stackPush(&stack, faces[i+j*isrf]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        if (nsecC0-2 == 0) break;
      }
    }
  }

  /* set face attributes */
  for (j = 0; j < nstripe; j++) {
    for (i = 0; i < nsecC0-1; i++) {

      attrint = j+1;
      stat    = EG_attributeAdd(faces[i+j*isrf], ".blendStrip", ATTRINT, 1,
                                &attrint, NULL, NULL);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Warning: Strip %d blendStrip = %d (EG_blend)!\n",
                 attrint, stat);
      }
    }
  }

  if (nface == 1 && bodyType == SHEETBODY) {
    stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1,
                           &faces[0], NULL, &body);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo FaceBody = %d (EG_blend)!\n", stat);
      if (body != NULL) EG_deleteObject(body);
      goto cleanup;
    }
  } else {

    /* put all the faces together */
    stat = EG_makeTopology(context, NULL, SHELL, bodyType == SOLIDBODY ? CLOSED : OPEN,
        NULL, nface, faces, NULL, &shell);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo Shell = %d (EG_blend)!\n", stat);
      goto cleanup;
    }
    stat = EG_stackPush(&stack, shell);
    if (stat != EGADS_SUCCESS) goto cleanup;

    body = NULL;
    if (bodyType == SOLIDBODY) {
      stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                             &shell, NULL, &body);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo SolidBody = %d (EG_blend)!\n", stat);
        if (body != NULL) EG_deleteObject(body);
        goto cleanup;
      }
    } else {
      stat = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1,
                             &shell, NULL, &body);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo SheetBody = %d (EG_blend)!\n", stat);
        if (body != NULL) EG_deleteObject(body);
        goto cleanup;
      }
    }
  }

  stat = EGADS_SUCCESS;
  *result = body;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_blend)!\n", stat);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (EG_blend)!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  EG_free(faces);
  EG_free(edgesU);
  EG_free(edgesV);

  EG_free(nodes);
  EG_free(secsC0);

  EG_free(cedges);
  EG_free(csens);
  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}


extern "C" int
EG_blend_dot(ego body, int nsex, const ego *secs,
             /*@null@*/ double *rc1, /*@null@*/ double *rc1_dot,
             /*@null@*/ double *rcN, /*@null@*/ double *rcN_dot)
{
  typedef SurrealS<1> T;
  int    i, j, k, kk, outLevel, stat, nstripe, ncap, cap, oclass, mtype;
  int    bodyType, uclosed, vclosed, allNodes, tip, te, isrf, inode, nface=0, nloop, nedge;
  int    *senses, planar, nchldrn=0, ntip, data_dot;
  int    nsec, nsecC0, begRC, endRC, tips[2]={0,0}, tipsec[2]={0,0}, iedge[2], found;
  double data[18], trange[2];
  T      rc1S[8], rcNS[8], *rc1s=NULL, *rcNs=NULL, result[3], mknot, ts[2], len;
  ego    context, surf, sec_surf, cap_surf, curv, ref, shell, *tedges, *nds;
  ego    *chldrn, *secsC0=NULL, *nodes=NULL, *edges, *faces=NULL, *loops;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;
  egSpline<T> tipsurfs[2], tipcurvs[4];
#ifdef PCURVE_SENSITIVITY
  int    ilen, rlen, *ivec=NULL, header[4];
  double *rvec=NULL
  T      sdata[8], *svec=NULL
#endif

  nsec = nsex;
  if (nsec < 0) nsec = -nsec;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL)               return EGADS_NOTCNTX;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec, secs, "EG_blend_dot",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* check data_dot */
  data_dot = 1;
  for (i = 0; i < nsec; i++) {
    if (EG_hasGeometry_dot(secs[i]) == EGADS_NOTFOUND) {
      data_dot = 0;
      if (outLevel > 0)
        printf(" EGADS Error: Section %d without data_dot (EG_blend_dot)!\n",
               i+1);
    }
  }

  if (data_dot == 0) return EGADS_NODATA;

  /* special case for all nodes */
  if (allNodes == 1) {

#ifdef EGADS_SPLINE_VELS
    stat = EG_blendNodeSpline(NULL, nsex, secs, vclosed,
                              &nsecC0, &secsC0, &nodes, &splcurvsU);
#else
    stat = EG_blendNodeSpline(nsex, secs, vclosed, &nodes, &splcurvsU);
#endif
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: ruledNodeSpline = %d (EG_blend_dot)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom_dot(body, nsecC0, secsC0, nodes, splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom_dot = %d (EG_blend_dot)!\n", stat);
      goto cleanup;
    }

    stat = EGADS_SUCCESS;
    goto cleanup;
  }

  /* sections with edges */

  /* check end conditions */
  begRC = endRC  = -1;
  if (rc1 != NULL) {
    i = 0;
    if (secs[i]->oclass == NODE) begRC = 1; /* nose treatment */
    else if (rc1[0] != 0       ) begRC = 2; /* tangencey */
  }
  if (rcN != NULL) {
    i = nsec-1;
    if (secs[i]->oclass == NODE) endRC = 1; /* nose treatment */
    else if (rcN[0] != 0       ) endRC = 2; /* tangencey */
  }
  if (((begRC == 1) || (endRC == 1)) && (nsec <= 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 2 for Nose Treatment (EG_blend_dot)!\n");
    return EGADS_GEOMERR;
  }
  if ((begRC == 1) && (endRC == 1) && (nsec <= 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 3 for 2Nose Treatment (EG_blend_dot)!\n");
    return EGADS_GEOMERR;
  }

  if ((rc1 != NULL) && (rc1_dot == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: rc1 != NULL but rc1_dot == NULL (EG_blend_dot)!\n");
    return EGADS_NODATA;
  }

  if ((rcN != NULL) && (rcN_dot == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: rcN != NULL but rcN_dot == NULL (EG_blend_dot)!\n");
    return EGADS_NODATA;
  }

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec, secs, "EG_blend_dot_secs.egads");
#endif

  /* extract nose treatment */
  if (begRC == 1) {
    for (i = 0; i < 8; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }

    /* normalize */
    len = sqrt(rc1S[1]*rc1S[1] + rc1S[2]*rc1S[2] + rc1S[3]*rc1S[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first south axis (EG_blend_dot)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[1] /= len;
    rc1S[2] /= len;
    rc1S[3] /= len;

    len = sqrt(rc1S[5]*rc1S[5] + rc1S[6]*rc1S[6] + rc1S[7]*rc1S[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second south axis (EG_blend_dot)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[5] /= len;
    rc1S[6] /= len;
    rc1S[7] /= len;

    rc1s = rc1S;
  }
  if (endRC == 1) {
    for (i = 0; i < 8; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }

    /* normalize */
    len = sqrt(rcNS[1]*rcNS[1] + rcNS[2]*rcNS[2] + rcNS[3]*rcNS[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first north axis (EG_blend_dot)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[1] /= len;
    rcNS[2] /= len;
    rcNS[3] /= len;

    len = sqrt(rcNS[5]*rcNS[5] + rcNS[6]*rcNS[6] + rcNS[7]*rcNS[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second north axis (EG_blend_dot)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[5] /= len;
    rcNS[6] /= len;
    rcNS[7] /= len;

    rcNs = rcNS;
  }

  /* extract tangency treatment */
  if (begRC == 2) {
    for (i = 0; i < 4; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }
    rc1s = rc1S;
  }
  if (endRC == 2) {
    for (i = 0; i < 4; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }
    rcNs = rcNS;
  }

  /* extract tip treatment */
  if ((secs[     0]->oclass == FACE || secs[     0]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }
    rc1s = rc1S;
  }
  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }
    rcNs = rcNS;
  }

#ifdef EGADS_SPLINE_VELS
  stat = EG_blendSpline<T>(NULL, nsex, secs, rc1s, rcNs, begRC, endRC,
                           nstripe, uclosed, vclosed, planar, &nsecC0, &secsC0,
                           &nodes, &splcurvsU, &splcurvsV, &splsurfs,
                           tip, te, tipsurfs, tipcurvs);
#else
  stat = EG_blendSpline<T>(nsex, secs, rc1s, rcNs, begRC, endRC,
                           nstripe, uclosed, vclosed, planar, &nsecC0, &secsC0,
                           &nodes, &splcurvsU, &splcurvsV, &splsurfs,
                           tip, te, tipsurfs, tipcurvs);
#endif
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: blendSpline = %d (EG_blend_dot)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom_dot(body, nsecC0, secsC0, nstripe, uclosed, vclosed,
                           ncap, te, tip,
                           nodes, splcurvsU, splcurvsV, splsurfs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: splineGeom = %d (EG_blend_dot)!\n", stat);
    goto cleanup;
  }

  inode = nsecC0;
  isrf  = nsecC0-1;

  if (ncap > 0) {
    /* get the shell */
    stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    shell = chldrn[0];

    /* get the faces */
    stat = EG_getTopology(shell, &ref, &oclass, &mtype, data, &nface, &faces,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* extract special wing tip treatment */
  cap  = 0;
  ntip = 0;
  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    tips[0]   = nface-ncap+cap;
    tipsec[0] = 0;
    cap++;
    ntip++;
  }

  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) && (ntip == 0)) cap++;

  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    tips[1]   = nface-ncap+cap;
    tipsec[1] = nsecC0-1;
    cap++;
    ntip++;
  }

  if (ncap > 0) {
    cap = 0;
    for (k = 0; k < 2; k++) {
      if ((tip & 1) && (k == 0)) continue;
      if ((tip & 2) && (k == 1)) continue;
      i = (k == 0) ? 0 : nsecC0-1;
      stat = EG_getTopology(secsC0[i], &sec_surf, &oclass, &mtype, data,
                            &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (oclass == BODY) {
        stat = EG_getTopology(chldrn[0], &sec_surf, &oclass, &mtype, data,
                              &nchldrn, &chldrn, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (oclass != FACE) continue;

      stat = EG_getTopology(faces[nface-ncap+cap], &cap_surf, &oclass, &mtype,
                            data, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_copyGeometry_dot(sec_surf, NULL, NULL, cap_surf);
      if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
      if (sec_surf->mtype != PLANE) {
        /* get the edges for the loop */
        stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, data, &nedge,
                              &edges, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set P-curve sensitivities to zero (though it's not really correct
         * as we need EG_otherCurve_dot...) */
        for (jj = 0; jj < nedge; jj++) {
          stat = EG_getGeometry(edges[nedge+jj], &oclass, &mtype, &ref, &ivec,
                                &rvec);
          if (stat != EGADS_SUCCESS) goto cleanup;

          EG_getGeometryLen(edges[nedge+jj], &ilen, &rlen);
          svec = (T*)EG_alloc(rlen*sizeof(T));
          for (j = 0; j < rlen; j++) svec[j] = rvec[j];

          stat = EG_setGeometry_dot(edges[nedge+jj], oclass, mtype, ivec, svec);
          if (stat != EGADS_SUCCESS) goto cleanup;

          EG_free(ivec); ivec=NULL;
          EG_free(rvec); rvec=NULL;
          EG_free(svec); svec=NULL;
        }
      }
#endif

      cap++;
    }
  }

  /* set sensitivities for cap Faces with wingTips */
  if ((ntip != 0) && (faces != NULL)) {

    for (k = 0; k < 2; k++) {
      if (!(tip & 1) && (k == 0)) continue;
      if (!(tip & 2) && (k == 1)) continue;

      /* get the surface and loop */
      stat = EG_getTopology(faces[tips[k]], &surf, &oclass, &mtype, data,
                            &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nloop != 1) {
        printf(" EGADS Error: Wing Tip Face %d Loop = %d (should be 1) (EG_blend_dot)!\n",
               tips[k], nloop);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      stat = EG_setGeometry_dot(surf, SURFACE, BSPLINE, tipsurfs[k].header,
                                tipsurfs[k].data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: setGeometry_dot Tip Face %d status %d(EG_blend_dot)!\n",
                 tips[k], stat);
        goto cleanup;
      }

      /* get the edges */
      stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, data, &nedge,
                            &tedges, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
      /* set P-curve sensitivities to zero which is correct */
      for (kk = 0; kk < nedge; kk++) {
        stat = EG_getGeometry(tedges[nedge+kk], &oclass, &mtype, &ref, &ivec,
                              &rvec);
        if (stat != EGADS_SUCCESS) goto cleanup;

        if (mtype != LINE) continue;

        EG_getGeometryLen(tedges[nedge+kk], &ilen, &rlen);
        svec = (T*)EG_alloc(rlen*sizeof(T));
        for (jj = 0; jj < rlen; jj++) svec[jj] = rvec[jj];

        stat = EG_setGeometry_dot(tedges[nedge+kk], PCURVE, LINE, NULL, svec);
        if (stat != EGADS_SUCCESS) goto cleanup;

        EG_free(ivec); ivec=NULL;
        EG_free(rvec); rvec=NULL;
        EG_free(svec); svec=NULL;
      }
#endif

      /* set the sensitivity for the t-range of the DEGENERATE edges */
      for (kk = 0; kk < nedge; kk++) {
        if (tedges[kk]->mtype == DEGENERATE) {
          ts[0] = 0.0;
          ts[1] = 1.0;
          stat = EG_setRange_dot(tedges[kk], EDGE, ts);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }
      }

      if (nstripe == 2) {

        /* check the edge count */
        if (nedge != 4) {
          printf(" EGADS Error: Wing Tip Face %d nEdge = %d (should be 4) (EG_blend_dot)!\n",
                 tips[k], nedge);
          stat = EGADS_TOPOERR;
          goto cleanup;
        }

      } else {

        /* check the edge count */
        if (nedge != 5) {
          printf(" EGADS Error: Wing Tip Face %d nEdge = %d (should be 5) (EG_blend_dot)!\n",
                 tips[k], nedge);
          stat = EGADS_TOPOERR;
          goto cleanup;
        }

        /* set the sensitivity of the TE tip curve(s) */

        iedge[0] = iedge[1] = -1;
        for (kk = 0; kk < nedge; kk++) {
          if (tedges[kk]->mtype == DEGENERATE) continue;

          /* re-discover the edge order in the loop */
          egadsEdge *pedge = (egadsEdge *) tedges[kk]->blind;

          ego node0 = pedge->nodes[0];
          ego node1 = pedge->nodes[1];

          /* exclude edges on the section */
          found = 0;
          for (j = 0; j < nstripe; j++) {
            if ((EG_isEquivalent(node0, nodes[tipsec[k]+ j   *inode]) == EGADS_SUCCESS) &&
                (EG_isEquivalent(node1, nodes[tipsec[k]+(j+1)*inode]) == EGADS_SUCCESS) ){
              found = 1;
            }
          }
          if (found == 1)  continue;

          if (EG_isEquivalent(node1, nodes[tipsec[k]+ te   *inode]) == EGADS_SUCCESS) {
            iedge[0] = kk;
          } else if (EG_isEquivalent(node0, nodes[tipsec[k]+(te+1)*inode]) == EGADS_SUCCESS) {
            iedge[1] = kk;
          }
        }

        /* check that all edges were found, excluding modified trailing edges */
        for (kk = 0; kk < 2; kk++) {
          if (iedge[kk] != -1) continue;

          printf(" EGADS Internal Error: Failed to determine loop ordering (EG_blend_dot)!\n");
          stat = EGADS_TOPOERR;
          goto cleanup;
        }

        /*** first Edge ***/
        stat = EG_getTopology(tedges[iedge[0]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(curv, CURVE, BSPLINE, tipcurvs[k].header,
                                  tipcurvs[k].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* get the midpoint knot value for the TE tip curve */
        mknot = tipcurvs[k].data[tipcurvs[k].header[3]/2];

        /* tip node sensitivity */
        stat = EG_spline1dEval(tipcurvs[k].header, tipcurvs[k].data, mknot,
                               result);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(nds[0], NODE, 0, NULL, result);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_copyGeometry_dot(nodes[tipsec[k]+ te   *inode], NULL, NULL,
                                   nds[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
        /* p-curve sensitivity */
        header[0] = 0;
        header[1] = 1;
        header[2] = 2;
        header[3] = 4;

        sdata[0] = sdata[1] = mknot;
        sdata[2] = sdata[3] = 1.0; /* v == 0 */
        sdata[4] = mknot; sdata[5] = 0.0;
        sdata[6] = 1.0;   sdata[7] = 0.0;
        stat = EG_setGeometry_dot(tedges[iedge[0]+nedge], PCURVE, BSPLINE,
                                  header, sdata);
        if (stat != EGADS_SUCCESS) goto cleanup;
#endif

        /* set the sensitivity for the t-range of the edge */
        ts[0] = mknot;
        ts[1] = 1.0;
        stat = EG_setRange_dot(tedges[iedge[0]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;


        /*** the 2nd Edge ***/
        stat = EG_getTopology(tedges[iedge[1]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(curv, CURVE, BSPLINE, tipcurvs[k].header,
                                  tipcurvs[k].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_copyGeometry_dot(nodes[tipsec[k]+(te+1)*inode], NULL, NULL,
                                   nds[0]);
        if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
        /* p-curve sensitivity */
        header[0] = 0;
        header[1] = 1;
        header[2] = 2;
        header[3] = 4;

        sdata[0] = sdata[1] = 0.0;
        sdata[2] = sdata[3] = mknot; /* v == 0 */
        sdata[4] = 0.0;   sdata[5] = 0.0;
        sdata[6] = mknot; sdata[7] = 0.0;
        stat = EG_setGeometry_dot(tedges[iedge[1]+nedge], PCURVE, BSPLINE,
                                  header, sdata);
        if (stat != EGADS_SUCCESS) goto cleanup;
#endif
        /* set the sensitivity for the t-range of the edge */
        ts[0] = 0.0;
        ts[1] = mknot;
        stat = EG_setRange_dot(tedges[iedge[1]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
    }

    /* set the trailing edge face sensitivities */
    if (nstripe == 3) {
      j = te;
      for (i = 0; i < nsecC0-1; i += nsecC0-2) {

        stat = EG_getTopology(faces[i+j*isrf], &surf, &oclass, &mtype, data,
                              &nloop, &loops, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(surf, SURFACE, BSPLINE,
                                  splsurfs[i+j*isrf].header,
                                  splsurfs[i+j*isrf].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* get the edges for the loop */
        stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, data, &nedge,
                              &edges, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
#error "Implement TE p-curve sensitvitiy"
#endif

        for (kk = 0; kk < nedge; kk++) {
          stat = EG_getRange(edges[kk], trange, &nloop);
          if (stat != EGADS_SUCCESS) goto cleanup;

          /* set the sensitivity for the t-range of the edge */
          ts[0] = trange[0];
          ts[1] = trange[1];
          stat = EG_setRange_dot(edges[kk], EDGE, ts);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        if (nsecC0-2 == 0) break;
      }
    }
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_blend_dot)!\n", stat);
  }

  /* clean up all of our temps */
  EG_free(nodes);
  EG_free(secsC0);
#ifdef PCURVE_SENSITIVITY
  EG_free(ivec);
  EG_free(rvec);
  EG_free(svec);
#endif

  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}


extern "C" int
EG_ruled(int nsec_in, const ego *secs_in, ego *result)
{
  typedef double T;
  int    i, j, k, n, outLevel, stat, nstripe, oclass, mtype, bodyType, planar;
  int    nsec, inode, icrvV, ncurvV=0, nface=0, degenclosed;
  int    uclosed, vclosed, allNodes, cap, ncap, *senses=NULL, *csens=NULL, lsens[1] = {SFORWARD};
  double data[18];
  ego    context, loop, ref, shell, body;
  ego    *secs=NULL, *chldrn;
  ego    *nodes=NULL, *edgesU=NULL, *edgesV=NULL, *faces=NULL, *cedges=NULL;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;
  objStack    stack;

  *result = NULL;
  if (nsec_in <= 1)                     return EGADS_EMPTY;
  if (secs_in == NULL)                  return EGADS_NULLOBJ;
  if (secs_in[0] == NULL)               return EGADS_NULLOBJ;
  if (secs_in[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs_in[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs_in[0]);
  context  = EG_context(secs_in[0]);
  if (context == NULL)                  return EGADS_NOTCNTX;

  /* create stack for gracefully cleaning up objects */
  stat  = EG_stackInit(&stack);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec_in, secs_in, "EG_ruled",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec, secs, "EG_ruled_secs.egads");
#endif

  stat = EG_uniqueSections(nsec_in, secs_in, nsec, &secs);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* special case for all nodes */
  if (allNodes == 1) {

#ifdef EGADS_SPLINE_VELS
    stat = EG_ruledNodeSpline(NULL, nsec, secs, vclosed, &nodes, &splcurvsU);
#else
    stat = EG_ruledNodeSpline(nsec, secs, vclosed, &nodes, &splcurvsU);
#endif
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: ruledNodeSpline = %d (EG_ruled)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom(nsec, secs, vclosed, nodes, splcurvsU,
                             &edgesU, &stack);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom = %d (EG_ruled)!\n", stat);
      goto cleanup;
    }

    senses = (int*)EG_alloc((nsec-1)*sizeof(int));
    if (senses == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocation (EG_ruled)!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for (i = 0; i < nsec-1; i++) senses[i] = SFORWARD;
    
    stat = EG_makeTopology(context, NULL, LOOP, vclosed == 1 ? CLOSED : OPEN, NULL,
                           nsec-1, edgesU, senses, &loop);
    EG_free(senses);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo Loop = %d (EG_ruled)!\n", stat);
      goto cleanup;
    }
    stat = EG_stackPush(&stack, loop);
    if (stat != EGADS_SUCCESS) goto cleanup;

    body = NULL;
    stat = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1,
                           &loop, NULL, &body);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo WireBody = %d (EG_ruled)!\n", stat);
      if (body != NULL) EG_deleteObject(body);
      goto cleanup;
    }

    stat = EGADS_SUCCESS;
    *result = body;
    goto cleanup;
  }

  /* sections with edges */
#ifdef EGADS_SPLINE_VELS
  stat = EG_ruledSpline(NULL, nsec, secs, nstripe, uclosed, vclosed, &nodes, &splcurvsU,
                        &splcurvsV, &splsurfs);
#else
  stat = EG_ruledSpline(nsec, secs, nstripe, uclosed, vclosed, &nodes, &splcurvsU,
                        &splcurvsV, &splsurfs);
#endif
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: ruledSpline = %d (EG_ruled)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom(nsec, secs, nstripe, uclosed, vclosed, ncap, -1, 0,
                       nodes, splcurvsU, splcurvsV, splsurfs,
                       &ncurvV, &edgesU, &edgesV, &nface, &faces, &stack);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: splineGeom = %d (EG_ruled)!\n", stat);
    goto cleanup;
  }

  inode = nsec;
  icrvV = nsec;

  /* add cap faces if necessary */
  if (ncap > 0) {

    cedges = (ego *) EG_alloc(2*nstripe*sizeof(ego));
    csens  = (int *) EG_alloc(  nstripe*sizeof(int));
    if ((cedges == NULL) || (csens == NULL)) {
      if (outLevel > 0)
        printf(" EGADS Error: Allocating edges for loop!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for (j = 0; j < 2*nstripe; j++) cedges[j] = NULL;
    for (j = 0; j <   nstripe; j++) csens[j]  = SFORWARD;

    cap = 0;
    for (k = 0; k < 2; k++) {
      i = (k == 0) ? 0 : nsec-1;
      stat = EG_getTopology(secs[i], &ref, &oclass, &mtype, data, &n,
                            &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (oclass == BODY) {
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &n,
                              &chldrn, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (oclass != FACE) continue;

      /* grab the spline edges on the cap */
      for (j = 0; j < nstripe; j++) cedges[j] = edgesV[i+j*icrvV];

      if (ref->mtype == PLANE) {

        /* make the loop (no need for surface with planar) */
        stat = EG_makeTopology(context, NULL, LOOP, CLOSED, NULL, nstripe,
                               cedges, csens, &loop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Cap %d Planar Loop = %d (EG_ruled)!\n",
                   cap, stat);
          goto cleanup;
        }

      } else {

        /* make PCurves for the loop */
        for (j = 0; j < nstripe; j++) {
          /* construct pcurves on the surface */
          double eps = 1.e-7;
          EG_tolerance(cedges[j], &eps);
          if (eps < 1.e-7) eps = 1.e-7;
          for (int jj = 0; jj < 2; jj++) {
            stat = EG_otherCurve(ref, cedges[j], eps, &cedges[j+nstripe]);
            if (stat != EGADS_SUCCESS) {
              if (jj == 0) {
                eps *= 100.0;
                continue;
              }
              if (outLevel > 0)
                printf(" EGADS Error: %d otherCurve First = %d (EG_ruled)!\n",
                       j, stat);
              goto cleanup;
            }
          }
          stat = EG_stackPush(&stack, cedges[j+nstripe]);
          if (stat != EGADS_SUCCESS) goto cleanup;
        }

        /* make the loop with the original surface */
        stat = EG_makeTopology(context, ref, LOOP, CLOSED, NULL, nstripe,
                               cedges, csens, &loop);
        if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: makeTopo Cap %d Loop = %d (EG_ruled)!\n",
                   cap, stat);
          goto cleanup;
        }
        for (j = nstripe; j < 2*nstripe; j++) cedges[j] = NULL;
      }
      stat = EG_stackPush(&stack, loop);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* make the face */
      stat = EG_makeTopology(context, ref, FACE, mtype,
                             NULL, 1, &loop, lsens, &faces[nface-ncap+cap]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo Cap %d Face = %d (EG_ruled)!\n",
                 cap, stat);
        goto cleanup;
      }
      stat = EG_stackPush(&stack, faces[nface-ncap+cap]);
      if (stat != EGADS_SUCCESS) goto cleanup;

      EG_attributeDup(secs[i], faces[nface-ncap+cap]);
      cap++;
    }
  }

  body = NULL;
  if (nface == 1 && bodyType == SHEETBODY) {
    stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1,
                           &faces[0], NULL, &body);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo FaceBody = %d (EG_ruled)!\n", stat);
      if (body != NULL) EG_deleteObject(body);
      goto cleanup;
    }
  } else {

    /* remove any degenerate (NULL) faces */
    j = 0;
    for (i = 0; i < nface; i++) {
      if (faces[i] == NULL) continue;
      faces[j++] = faces[i];
    }
    nface = j;

    /* check if degeneracies creates a closed solid */
    degenclosed = 1;
    for (i = 0; i < nsec-1; i++) {
      j = 0;
      if (nodes[i+1+j*inode] != nodes[i+j*inode]) degenclosed = 0;
      j = nstripe;
      if (nodes[i+1+j*inode] != nodes[i+j*inode]) degenclosed = 0;
    }
    if (vclosed && degenclosed) bodyType = SOLIDBODY;

    /* put all the faces together */
    stat = EG_makeTopology(context, NULL, SHELL, bodyType == SOLIDBODY ? CLOSED : OPEN,
                           NULL, nface, faces, NULL, &shell);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeTopo Shell = %d (EG_ruled)!\n", stat);
      goto cleanup;
    }
    stat = EG_stackPush(&stack, shell);
    if (stat != EGADS_SUCCESS) goto cleanup;

    if (bodyType == SOLIDBODY) {
      stat = EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, 1,
                             &shell, NULL, &body);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo SolidBody = %d (EG_ruled)!\n", stat);
        if (body != NULL) EG_deleteObject(body);
        goto cleanup;
      }
    } else {
      stat = EG_makeTopology(context, NULL, BODY, SHEETBODY, NULL, 1,
                             &shell, NULL, &body);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: makeTopo SheetBody = %d (EG_ruled)!\n", stat);
        if (body != NULL) EG_deleteObject(body);
        goto cleanup;
      }
    }
  }

  stat = EGADS_SUCCESS;
  *result = body;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_ruled)!\n", stat);
  }

  /* clean up all of our temps */
  EG_stackPop(&stack, &ref);
  while (ref != NULL) {
    i = EG_deleteObject(ref);
    if ((i != EGADS_SUCCESS) && (outLevel > 0))
      printf(" EGADS Internal: EG_deleteObject = %d (EG_ruled)!\n", i);
    EG_stackPop(&stack, &ref);
  }
  EG_stackFree(&stack);

  EG_free(secs);
  EG_free(faces);
  EG_free(edgesU);
  EG_free(edgesV);
  EG_free(nodes);

  EG_free(cedges);
  EG_free(csens);
  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}


extern "C" int
EG_ruled_dot(ego body, int nsec_in, const ego *secs_in)
{
  typedef SurrealS<1> T;
  int    i, k, outLevel, stat, nstripe, oclass, mtype, bodyType, planar;
  int    nsec, uclosed, vclosed, allNodes, data_dot, ncap, cap, nchldrn, nface, *senses;
  double data[18];
  ego    context, ref, shell, sec_surf, cap_surf;
  ego    *secs=NULL, *chldrn, *nodes=NULL, *faces;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;
#ifdef PCURVE_SENSITIVITY
  int    ilen, rlen, *ivec=NULL, nedge, header[4];
  double *rvec=NULL
  T      *svec=NULL;
#endif

  if (nsec_in <= 1)                     return EGADS_EMPTY;
  if (secs_in == NULL)                  return EGADS_NULLOBJ;
  if (secs_in[0] == NULL)               return EGADS_NULLOBJ;
  if (secs_in[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs_in[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs_in[0]);
  context  = EG_context(secs_in[0]);
  if (context == NULL)                  return EGADS_NOTCNTX;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec_in, secs_in, "EG_ruled_dot",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* check input has data_dot */
  data_dot = 1;
  for (i = 0; i < nsec_in; i++) {
    if (EG_hasGeometry_dot(secs_in[i]) == EGADS_NOTFOUND) {
      data_dot = 0;
      if (outLevel > 0)
        printf(" EGADS Error: Section %d no data_dot (EG_ruled_dot)!\n",
               i+1);
    }
  }

  if (data_dot == 0) return EGADS_NODATA;

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec_in, secs_in, "EG_ruled_dot_secs.egads");
#endif

  stat = EG_uniqueSections(nsec_in, secs_in, nsec, &secs);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* special case for all nodes */
  if (allNodes == 1) {

#ifdef EGADS_SPLINE_VELS
    stat = EG_ruledNodeSpline(NULL, nsec, secs, vclosed, &nodes, &splcurvsU);
#else
    stat = EG_ruledNodeSpline(nsec, secs, vclosed, &nodes, &splcurvsU);
#endif
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: ruledNodeSpline = %d (EG_ruled_dot)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom_dot(body, nsec, secs, nodes, splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom_dot = %d (EG_ruled_dot)!\n", stat);
      goto cleanup;
    }

    stat = EGADS_SUCCESS;
    goto cleanup;
  }

  /* sections with edges */
#ifdef EGADS_SPLINE_VELS
  stat = EG_ruledSpline(NULL, nsec, secs, nstripe, uclosed, vclosed, &nodes, &splcurvsU,
                        &splcurvsV, &splsurfs);
#else
  stat = EG_ruledSpline(nsec, secs, nstripe, uclosed, vclosed, &nodes, &splcurvsU,
                        &splcurvsV, &splsurfs);
#endif
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: ruledSpline = %d (EG_ruled_dot)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom_dot(body, nsec, secs, nstripe, uclosed, vclosed, ncap,
                           -1, 0, nodes, splcurvsU, splcurvsV, splsurfs);
  if (stat != EGADS_SUCCESS) {
     if (outLevel > 0)
       printf(" EGADS Error: splineGeom_dot = %d (EG_ruled_dot)!\n", stat);
     goto cleanup;
   }

  /* copy over cap face sensitivities (curves and nodes are already set) */
  if (ncap > 0) {

    /* get the shell */
    stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    shell = chldrn[0];

    /* get the faces */
    stat = EG_getTopology(shell, &ref, &oclass, &mtype, data, &nface, &faces,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    cap = 0;
    for (k = 0; k < 2; k++) {
      i = (k == 0) ? 0 : nsec-1;
      stat = EG_getTopology(secs[i], &sec_surf, &oclass, &mtype, data,
                            &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (oclass == BODY) {
        stat = EG_getTopology(chldrn[0], &sec_surf, &oclass, &mtype, data,
                              &nchldrn, &chldrn, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (oclass != FACE) continue;

      stat = EG_getTopology(faces[nface-ncap+cap], &cap_surf, &oclass, &mtype,
                            data, &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_copyGeometry_dot(sec_surf, NULL, NULL, cap_surf);
      if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef PCURVE_SENSITIVITY
      if (sec_surf->mtype != PLANE) {
        /* get the edges for the loop */
        stat = EG_getTopology(chldrn[0], &ref, &oclass, &mtype, data, &nedge,
                              &edges, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set P-curve sensitivities to zero (though it's not really correct
         * as we need EG_otherCurve_dot...) */
        for (jj = 0; jj < nedge; jj++) {
          stat = EG_getGeometry(edges[nedge+jj], &oclass, &mtype, &ref, &ivec,
                                &rvec);
          if (stat != EGADS_SUCCESS) goto cleanup;

          EG_getGeometryLen(edges[nedge+jj], &ilen, &rlen);
          svec = (T*)EG_alloc(rlen*sizeof(T));
          for (j = 0; j < rlen; j++) svec[j] = rvec[j];

          stat = EG_setGeometry_dot(edges[nedge+jj], oclass, mtype, ivec, svec);
          if (stat != EGADS_SUCCESS) goto cleanup;

          EG_free(ivec); ivec=NULL;
          EG_free(rvec); rvec=NULL;
          EG_free(svec); svec=NULL;
        }
      }
#endif

      cap++;
    }
  }

  stat = EGADS_SUCCESS;

cleanup:

#ifdef PCURVE_SENSITIVITY
  EG_free(ivec);
  EG_free(rvec);
  EG_free(svec);
#endif

  EG_free(secs);
  EG_free(nodes);
  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}


#ifdef EGADS_SPLINE_VELS
extern "C"
int EG_blend_vels( int nsex, const ego *secs,
                   /*@null@*/ const double *rc1,
                   /*@null@*/ const double *rc1_dot,
                   /*@null@*/ const double *rcN,
                   /*@null@*/ const double *rcN_dot,
                   egadsSplineVels *vels, ego body )
{
  typedef SurrealS<1> T;
  int    i, j, k, jj, kk, outLevel, stat, nstripe, ncap, cap, oclass, mtype, bodyType;
  int    uclosed, vclosed, allNodes, tip, te, isrf, inode, nface=0, nloop, nedge;
  int    *senses, planar, nchldrn=0, ntip, *ivec=NULL, ilen, rlen;
  int    nsec, nsecC0, begRC, endRC, tips[2]={0,0}, tipsec[2]={0,0}, iedge[2], found;
  double data[18], *rvec=NULL;
  T      rc1S[8], rcNS[8], *rc1s=NULL, *rcNs=NULL, result[3], mknot, ts[2], len;
  T      *rinfo=NULL;
  ego    context, surf, curv, ref, shell, sec_surf, cap_surf;
  ego    *chldrn, *secsC0=NULL, *nodes=NULL, *faces=NULL;
  ego    *loops, *tedges, *nds;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;
  egSpline<T> tipsurfs[2], tipcurvs[4];

  nsec    = nsex;
  if (nsec < 0) nsec = -nsec;
  if (nsec <= 1)                     return EGADS_EMPTY;
  if (secs == NULL)                  return EGADS_NULLOBJ;
  if (vels == NULL)                  return EGADS_NULLOBJ;
  if (secs[0] == NULL)               return EGADS_NULLOBJ;
  if (secs[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs[0]);
  context  = EG_context(secs[0]);
  if (context == NULL)               return EGADS_NOTCNTX;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec, secs, "EG_blend_vels",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* special case for all nodes */
  if (allNodes == 1) {

    stat = EG_blendNodeSpline(vels, nsex, secs, vclosed,
                              &nsecC0, &secsC0, &nodes, &splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: ruledNodeSpline = %d (EG_blend_vels)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom_dot(body, nsecC0, secsC0, nodes, splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom_dot = %d (EG_blend_vels)!\n", stat);
      goto cleanup;
    }

    /* remove any temporary sensitivities from the sections */
    for (i = 0; i < nsec; i++) {
      stat = EG_setGeometry_dot(secs[i], 0, 0, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    stat = EGADS_SUCCESS;
    goto cleanup;
  }

  /* sections with edges */

  /* check end conditions but not data_dot*/
  begRC = endRC  = -1;
  if (rc1 != NULL) {
    i = 0;
    if (secs[i]->oclass == NODE) begRC = 1; /* nose treatment */
    else if (rc1[0] != 0       ) begRC = 2; /* tangencey */
  }
  if (rcN != NULL) {
    i = nsec-1;
    if (secs[i]->oclass == NODE) endRC = 1; /* nose treatment */
    else if (rcN[0] != 0       ) endRC = 2; /* tangencey */
  }
  if (((begRC == 1) || (endRC == 1)) && (nsec <= 2)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 2 for Nose Treatment (EG_blend_dot)!\n");
    return EGADS_GEOMERR;
  }
  if ((begRC == 1) && (endRC == 1) && (nsec <= 3)) {
    if (outLevel > 0)
      printf(" EGADS Error: nsec must be > 3 for 2Nose Treatment (EG_blend_dot)!\n");
    return EGADS_GEOMERR;
  }

  if ((rc1 != NULL) && (rc1_dot == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: rc1 != NULL but rc1_dot == NULL (EG_blend_dot)!\n");
    return EGADS_NODATA;
  }

  if ((rcN != NULL) && (rcN_dot == NULL)) {
    if (outLevel > 0)
      printf(" EGADS Error: rcN != NULL but rcN_dot == NULL (EG_blend_dot)!\n");
    return EGADS_NODATA;
  }

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec, secs, "EG_blend_vels_secs.egads");
#endif

  /* extract nose treatment */
  if (begRC == 1) {
    for (i = 0; i < 8; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }

    /* normalize */
    len = sqrt(rc1S[1]*rc1S[1] + rc1S[2]*rc1S[2] + rc1S[3]*rc1S[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first south axis (EG_blend_vels)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[1] /= len;
    rc1S[2] /= len;
    rc1S[3] /= len;

    len = sqrt(rc1S[5]*rc1S[5] + rc1S[6]*rc1S[6] + rc1S[7]*rc1S[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second south axis (EG_blend_vels)!\n");
      return EGADS_GEOMERR;
    }
    rc1S[5] /= len;
    rc1S[6] /= len;
    rc1S[7] /= len;

    rc1s = rc1S;
  }
  if (endRC == 1) {
    for (i = 0; i < 8; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }

    /* normalize */
    len = sqrt(rcNS[1]*rcNS[1] + rcNS[2]*rcNS[2] + rcNS[3]*rcNS[3]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude first north axis (EG_blend_vels)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[1] /= len;
    rcNS[2] /= len;
    rcNS[3] /= len;

    len = sqrt(rcNS[5]*rcNS[5] + rcNS[6]*rcNS[6] + rcNS[7]*rcNS[7]);
    if (len == 0) {
      if (outLevel > 0)
        printf(" EGADS Error: Zero magnitude second north axis (EG_blend_vels)!\n");
      return EGADS_GEOMERR;
    }
    rcNS[5] /= len;
    rcNS[6] /= len;
    rcNS[7] /= len;

    rcNs = rcNS;
  }

  /* extract tangency treatment */
  if (begRC == 2) {
    for (i = 0; i < 4; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }
    rc1s = rc1S;
  }
  if (endRC == 2) {
    for (i = 0; i < 4; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }
    rcNs = rcNS;
  }

  /* extract tip treatment */
  if ((secs[     0]->oclass == FACE || secs[     0]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rc1S[i].value() = rc1[i];
      rc1S[i].deriv() = rc1_dot[i];
    }
    rc1s = rc1S;
  }
  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (planar == 0) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    for (i = 0; i < 2; i++) {
      rcNS[i].value() = rcN[i];
      rcNS[i].deriv() = rcN_dot[i];
    }
    rcNs = rcNS;
  }

  stat = EG_blendSpline<T>(vels, nsex, secs, rc1s, rcNs, begRC, endRC,
                           nstripe, uclosed, vclosed, planar, &nsecC0, &secsC0,
                           &nodes, &splcurvsU, &splcurvsV, &splsurfs,
                           tip, te, tipsurfs, tipcurvs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: blendSpline = %d (EG_blend_vels)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom_dot(body, nsecC0, secsC0, nstripe, uclosed, vclosed, ncap,
                           te, tip, nodes, splcurvsU, splcurvsV, splsurfs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: splineGeom = %d (EG_blend_vels)!\n", stat);
    goto cleanup;
  }

  inode = nsecC0;
  isrf  = nsecC0-1;

  if (ncap > 0) {
    /* get the shell */
    stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    shell = chldrn[0];

    /* get the faces */
    stat = EG_getTopology(shell, &ref, &oclass, &mtype, data, &nface, &faces,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  /* extract special wing tip treatment */
  cap  = 0;
  ntip = 0;
  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) &&
      (rc1 != NULL) && (rc1[0] == 0.0)) {
    tips[0]   = nface-ncap+cap;
    tipsec[0] = 0;
    cap++;
    ntip++;
  }

  if ((secs[0]->oclass == FACE || secs[0]->mtype == FACEBODY) && (ntip == 0)) cap++;

  if ((secs[nsec-1]->oclass == FACE || secs[nsec-1]->mtype == FACEBODY) &&
      (rcN != NULL) && (rcN[0] == 0.0)) {
    tips[1]   = nface-ncap+cap;
    tipsec[1] = nsecC0-1;
    cap++;
    ntip++;
  }

  if (ncap > 0) {
    cap = 0;
    for (k = 0; k < 2; k++) {
      if ((tip & 1) && (k == 0)) continue;
      if ((tip & 2) && (k == 1)) continue;
      i = (k == 0) ? 0 : nsecC0-1;
      stat = EG_getTopology(secsC0[i], &sec_surf, &oclass, &mtype, data,
                            &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (oclass == BODY) {
        stat = EG_getTopology(chldrn[0], &sec_surf, &oclass, &mtype, data,
                              &nchldrn, &chldrn, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
      if (oclass != FACE) continue;

      stat = EG_getTopology(faces[nface-ncap+cap], &cap_surf, &oclass, &mtype,
                            data, &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* Cap surface sensitivity is optional for EG_blend_vels */
      if (EG_hasGeometry_dot(sec_surf) == EGADS_SUCCESS) {
        stat = EG_copyGeometry_dot(sec_surf, NULL, NULL, cap_surf);
        if (stat != EGADS_SUCCESS) goto cleanup;
      } else {
        /* set sensitivity to INFTY so EG_hasGeometry_dot  and
                                       EG_copyGeometry_dot work */
        stat = EG_getGeometry(cap_surf, &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) goto cleanup;
        EG_getGeometryLen(cap_surf, &ilen, &rlen);
        rinfo = (T*)EG_alloc(rlen*sizeof(T));
        for (jj = 0; jj < rlen; jj++) {
          rinfo[jj].value() = rvec[jj];
          rinfo[jj].deriv() = INFTY;
        }
        stat = EG_setGeometry_dot(cap_surf, oclass, mtype, ivec, rinfo);
        if (stat != EGADS_SUCCESS) goto cleanup;

        EG_free(ivec);  ivec=NULL;
        EG_free(rvec);  rvec=NULL;
        EG_free(rinfo); rinfo=NULL;
      }

      cap++;
    }
  }

  /* set sensitivities for cap Faces with wingTips */
  if ((ntip != 0) && (faces != NULL)) {

    for (k = 0; k < 2; k++) {
      if (!(tip & 1) && (k == 0)) continue;
      if (!(tip & 2) && (k == 1)) continue;

      /* get the surface and loop */
      stat = EG_getTopology(faces[tips[k]], &surf, &oclass, &mtype, data,
                            &nloop, &loops, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;
      if (nloop != 1) {
        printf(" EGADS Error: Wing Tip Face %d Loop = %d (should be 1) (EG_blend_vels)!\n",
               tips[k], nloop);
        stat = EGADS_TOPOERR;
        goto cleanup;
      }

      stat = EG_setGeometry_dot(surf, SURFACE, BSPLINE, tipsurfs[k].header,
                                tipsurfs[k].data);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: setGeometry_dot Tip Face %d status %d(EG_blend_vels)!\n",
                 tips[k], stat);
        goto cleanup;
      }

      if (nstripe == 2) {

        /* check the edge count */
        stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, data, &nedge,
                              &tedges, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
        if (nedge != 4) {
          printf(" EGADS Error: Wing Tip Face %d nEdge = %d (should be 4) (EG_blend_vels)!\n",
                 tips[k], nedge);
          stat = EGADS_TOPOERR;
          goto cleanup;
        }
      } else {

        /* get the edges */
        stat = EG_getTopology(loops[0], &ref, &oclass, &mtype, data, &nedge,
                              &tedges, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;
        if (nedge != 5) {
          printf(" EGADS Error: Wing Tip Face %d nEdge = %d (should be 5) (EG_blend_vels)!\n",
                 tips[k], nedge);
          stat = EGADS_TOPOERR;
          goto cleanup;
        }

        /* set the sensitivity of the TE tip curve(s) */

        iedge[0] = iedge[1] = -1;
        for (kk = 0; kk < nedge; kk++) {
          if (tedges[kk]->mtype == DEGENERATE) continue;

          /* re-discover the edge order in the loop */
          egadsEdge *pedge = (egadsEdge *) tedges[kk]->blind;

          ego node0 = pedge->nodes[0];
          ego node1 = pedge->nodes[1];

          /* exclude edges on the section */
          found = 0;
          for (j = 0; j < nstripe; j++) {
            if ((EG_isEquivalent(node0, nodes[tipsec[k]+ j   *inode]) == EGADS_SUCCESS) &&
                (EG_isEquivalent(node1, nodes[tipsec[k]+(j+1)*inode]) == EGADS_SUCCESS)) {
              found = 1;
            }
          }
          if (found == 1) continue;

          if (EG_isEquivalent(node1, nodes[tipsec[k]+ te   *inode]) == EGADS_SUCCESS) {
            iedge[0] = kk;
          } else if (EG_isEquivalent(node0, nodes[tipsec[k]+(te+1)*inode]) == EGADS_SUCCESS) {
            iedge[1] = kk;
          }
        }

        /* check that all edges were found, excluding modified trailing edges */
        for (kk = 0; kk < 2; kk++) {
          if (iedge[kk] != -1) continue;

          printf(" EGADS Internal Error: Failed to determine loop ordering (EG_blend_vels)!\n");
          stat = EGADS_TOPOERR;
          goto cleanup;
        }


        /*** the 1st edge ***/
        stat = EG_getTopology(tedges[iedge[0]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(curv, CURVE, BSPLINE, tipcurvs[k].header,
                                  tipcurvs[k].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* get the midpoint knot value for the TE tip curve */
        mknot = tipcurvs[k].data[tipcurvs[k].header[3]/2];

        /* tip node sensitivity */
        stat = EG_spline1dEval(tipcurvs[k].header, tipcurvs[k].data, mknot,
                               result);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(nds[0], NODE, 0, NULL, result);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_copyGeometry_dot(nodes[tipsec[k]+ te   *inode], NULL, NULL,
                                   nds[1]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity for the t-range of the edge */
        ts[0] = mknot;
        ts[1] = 1.0;
        stat = EG_setRange_dot(tedges[iedge[0]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;


        /*** the 2nd edge ***/
        stat = EG_getTopology(tedges[iedge[1]], &curv, &oclass, &mtype, data,
                              &nchldrn, &nds, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(curv, CURVE, BSPLINE, tipcurvs[k].header,
                                  tipcurvs[k].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_copyGeometry_dot(nodes[tipsec[k]+(te+1)*inode], NULL, NULL,
                                   nds[0]);
        if (stat != EGADS_SUCCESS) goto cleanup;

        /* set the sensitivity for the t-range of the edge */
        ts[0] = 0.0;
        ts[1] = mknot;
        stat = EG_setRange_dot(tedges[iedge[1]], EDGE, ts);
        if (stat != EGADS_SUCCESS) goto cleanup;
      }
    }

    /* set the trailing edge face sensitivities */
    if (nstripe == 3) {
      j = te;
      for (i = 0; i < nsecC0-1; i += nsecC0-2) {

        stat = EG_getTopology(faces[i+j*isrf], &surf, &oclass, &mtype, data,
                              &nloop, &loops, &senses);
        if (stat != EGADS_SUCCESS) goto cleanup;

        stat = EG_setGeometry_dot(surf, SURFACE, BSPLINE,
                                  splsurfs[i+j*isrf].header,
                                  splsurfs[i+j*isrf].data);
        if (stat != EGADS_SUCCESS) goto cleanup;

        if (nsecC0-2 == 0) break;
      }
    }
  }

  /* remove any temporary sensitivities from the sections */
  for (i = 0; i < nsec; i++) {
    stat = EG_setGeometry_dot(secs[i], 0, 0, NULL, NULL);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_blend_vels)!\n", stat);
  }

  /* clean up all of our temps */
  EG_free(nodes);
  EG_free(secsC0);
  EG_free(ivec);
  EG_free(rvec);
  EG_free(rinfo);

  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}


extern "C"
int EG_ruled_vels(int nsec_in, const ego *secs_in, egadsSplineVels *vels, ego body)
{
  typedef SurrealS<1> T;
  int    i, k, jj, outLevel, stat, nstripe, oclass, mtype, bodyType, planar, rlen;
  int    nsec=0, uclosed, vclosed, allNodes, ncap, cap, nchldrn, nface, *senses, *ivec=NULL, ilen;
  double data[18], *rvec=NULL;
  T      *rinfo=NULL;
  ego    *secs=NULL, context, ref, shell, sec_surf, cap_surf;
  ego    *chldrn, *nodes=NULL, *faces;
  egSpline<T> *splsurfs=NULL, *splcurvsU=NULL, *splcurvsV=NULL;

  if (nsec_in <= 1)                     return EGADS_EMPTY;
  if (secs_in == NULL)                  return EGADS_NULLOBJ;
  if (vels == NULL)                     return EGADS_NULLOBJ;
  if (secs_in[0] == NULL)               return EGADS_NULLOBJ;
  if (secs_in[0]->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (EG_sameThread(secs_in[0]))        return EGADS_CNTXTHRD;
  outLevel = EG_outLevel(secs_in[0]);
  context  = EG_context(secs_in[0]);
  if (context == NULL)                  return EGADS_NOTCNTX;

  /* look at the input and check to see if OK */
  stat = EG_checkSections(nsec_in, secs_in, "EG_ruled_vels",
                          nstripe, ncap, bodyType, planar,
                          uclosed, vclosed, allNodes);
  if (stat != EGADS_SUCCESS) goto cleanup;

#ifdef DUMP_SECTIONS
  EG_dumpSections(nsec_in, secs_in, "EG_ruled_vels_secs.egads");
#endif

  stat = EG_uniqueSections(nsec_in, secs_in, nsec, &secs);
  if (stat != EGADS_SUCCESS) goto cleanup;

  /* special case for all nodes */
  if (allNodes == 1) {

    stat = EG_ruledNodeSpline(vels, nsec, secs, vclosed, &nodes, &splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: ruledNodeSpline = %d (EG_ruled_vels)!\n", stat);
      goto cleanup;
    }

    stat = EG_splineNodeGeom_dot(body, nsec, secs, nodes, splcurvsU);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: splineNodeGeom_dot = %d (EG_ruled_vels)!\n", stat);
      goto cleanup;
    }

    /* remove any temporary sensitivities from the sections */
    for (i = 0; i < nsec; i++) {
      stat = EG_setGeometry_dot(secs[i], 0, 0, NULL, NULL);
      if (stat != EGADS_SUCCESS) goto cleanup;
    }

    stat = EGADS_SUCCESS;
    goto cleanup;
  }

  /* sections with edges */
  stat = EG_ruledSpline(vels, nsec, secs, nstripe, uclosed, vclosed, &nodes, &splcurvsU,
                        &splcurvsV, &splsurfs);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: ruledSpline = %d (EG_ruled_vels)!\n", stat);
    goto cleanup;
  }

  stat = EG_splineGeom_dot(body, nsec, secs, nstripe, uclosed, vclosed, ncap,
                           -1, 0, nodes, splcurvsU, splcurvsV, splsurfs);
  if (stat != EGADS_SUCCESS) {
     if (outLevel > 0)
       printf(" EGADS Error: splineGeom_dot = %d (EG_ruled_vels)!\n", stat);
     goto cleanup;
   }

  /* copy over cap face sensitivities (curves and nodes are already set) */
  if (ncap > 0) {

    /* get the shell */
    stat = EG_getTopology(body, &ref, &oclass, &mtype, data, &nchldrn, &chldrn,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;
    shell = chldrn[0];

    /* get the faces */
    stat = EG_getTopology(shell, &ref, &oclass, &mtype, data, &nface, &faces,
                          &senses);
    if (stat != EGADS_SUCCESS) goto cleanup;

    cap = 0;
    for (k = 0; k < 2; k++) {
      i = (k == 0) ? 0 : nsec-1;
      if ((secs[i]->oclass != FACE)) continue;

      stat = EG_getTopology(secs[i], &sec_surf, &oclass, &mtype, data,
                            &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      stat = EG_getTopology(faces[nface-ncap+cap], &cap_surf, &oclass, &mtype,
                            data, &nchldrn, &chldrn, &senses);
      if (stat != EGADS_SUCCESS) goto cleanup;

      /* Cap surface sensitivity is optional for EG_ruled_vels */
      if (EG_hasGeometry_dot(sec_surf) == EGADS_SUCCESS) {
        stat = EG_copyGeometry_dot(sec_surf, NULL, NULL, cap_surf);
        if (stat != EGADS_SUCCESS) goto cleanup;
      } else {
        /* set sensitivity to INFTY so EG_hasGeometry_dot  and
                                       EG_copyGeometry_dot work */
        stat = EG_getGeometry(cap_surf, &oclass, &mtype, &ref, &ivec, &rvec);
        if (stat != EGADS_SUCCESS) goto cleanup;
        EG_getGeometryLen(cap_surf, &ilen, &rlen);
        rinfo = (T*)EG_alloc(rlen*sizeof(T));
        for (jj = 0; jj < rlen; jj++) {
          rinfo[jj].value() = rvec[jj];
          rinfo[jj].deriv() = INFTY;
        }
        stat = EG_setGeometry_dot(cap_surf, oclass, mtype, ivec, rinfo);
        if (stat != EGADS_SUCCESS) goto cleanup;

        EG_free(ivec);  ivec=NULL;
        EG_free(rvec);  rvec=NULL;
        EG_free(rinfo); rinfo=NULL;
      }

      cap++;
    }
  }

  /* remove any temporary sensitivities from the sections */
  for (i = 0; i < nsec; i++) {
    stat = EG_setGeometry_dot(secs[i], 0, 0, NULL, NULL);
    if (stat != EGADS_SUCCESS) goto cleanup;
  }

  stat = EGADS_SUCCESS;

cleanup:
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: exit with status = %d (EG_ruled_vels)!\n", stat);
  }

  /* clean up all of our temps */
  EG_free(secs);
  EG_free(nodes);
  EG_free(ivec);
  EG_free(rvec);
  EG_free(rinfo);

  delete [] splsurfs;
  delete [] splcurvsU;
  delete [] splcurvsV;

  return stat;
}
#endif
