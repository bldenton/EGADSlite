/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADSlite Geometry Evaluation (non-recursive) Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef WIN32
#define false 0
#define true  (!false)
typedef int bool;
#else
#include <stdbool.h>
#endif

#include "egadsTypes.h"
#include "egadsInternals.h"
#include "liteString.h"
#include "liteClasses.h"


typedef struct {
  int    uk;
  int    vk;
  double uv[2];
  double dist2;
} liteIndex;


#define PARAMACC         1.0e-4         /* parameter accuracy */
#define MAXDEG          25
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#define CROSS(a,b,c)       a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                           a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                           a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif
#ifdef __PROTO_H_AND_D__
#undef __PROTO_H_AND_D__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#define __PROTO_H_AND_D__   extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#define __PROTO_H_AND_D__ extern
#endif


__PROTO_H_AND_D__ void EG_EvaluateQuotientRule(  int dim, int der_count,
                                                 int v_stride, double *v );
__PROTO_H_AND_D__ void EG_EvaluateQuotientRule2( int dim, int der_count,
                                                 int v_stride, double *v );
__PROTO_H_AND_D__ bool EG_Bezier1DRat( int dim, int order, int cv_stride,
                                       const double *cv, double t0, double t1,
                                       int der_count, double t, int v_stride,
                                       double *v );

__PROTO_H_AND_D__ int  EG_attributeRet( const egObject *obj, const char *name,
                                        int *atype, int *len,
                                        /*@null@*/ const int **ints,
                                        /*@null@*/ const double **reals,
                                        /*@null@*/ const char **str );
__PROTO_H_AND_D__ int  EG_getRange( const egObject *geom, double *range,
                                    int *periodic );
__PROTO_H_AND_D__ int  EG_evaluate( const egObject *geom, const double *param,
                                    double *result );
__PROTO_H_AND_D__ int  EG_invEvaluate( const egObject *geom, double *xyz,
                                       double *param, double *result );
__PROTO_H_AND_D__ int  EG_getTopology( const egObject *topo, egObject **geom,
                                       int *oclass, int *mtype,
                                       /*@null@*/ double *limits,
                                       int *nChildren,
                                       egObject ***children, int **senses );
__PROTO_H_AND_D__ int  EG_getEdgeUV( const egObject *face,
                                     const egObject *edge,
                                     int sense, double t, double *result );
__PROTO_H_AND_D__ int  EG_getBodyTopos( const egObject *body,
                                        /*@null@*/ egObject *src,
                                        int oclass, int *ntopo,
                                        /*@null@*/ egObject ***topos );
__PROTO_H_AND_D__ int  EG_getBody( const egObject *object, egObject **body );



__HOST_AND_DEVICE__ static int
FindSpan(int nKnots, int degree, double u, double *U)
{
  int n, low, mid, high;
  
  if (u <= U[degree]) return degree;
  n = nKnots - degree - 1;
  if (u >= U[n]) return n-1;
  
  low  = degree;
  high = n;
  mid  = (low+high)/2;
  while ((u < U[mid]) || (u >= U[mid+1])) {
    if (u < U[mid]) {
      high = mid;
    } else {
      low  = mid;
    }
    mid = (low+high)/2;
  }
  
  return mid;
}


__HOST_AND_DEVICE__ static void
DersBasisFuns(int i, int p, double u, double *knot, int der, double **ders)
{
  int    j, k, j1, j2, r, s1, s2, rk, pk;
#ifndef __CUDA_ARCH__
  double d, saved, temp, ndu[MAXDEG+1][MAXDEG+1];
  double a[2][MAXDEG+1], left[MAXDEG+1], right[MAXDEG+1];
#else
  double d, saved, temp, **ndu;
  double *a[2], *left, *right;

  ndu = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (ndu == NULL) {
    printf("Malloc error\n");
    return;
  }
  for (j = 0; j < MAXDEG+1; ++j) {
    ndu[j] = (double *) malloc((MAXDEG+1)*sizeof(double));
    if (ndu[j] == NULL) {
      printf("Malloc error %d\n",j);
      return;
    }
  }
  for (j = 0; j < 2; ++j) {
    a[j] = (double*)malloc((MAXDEG+1)*sizeof(double));
    if (a[j] == NULL) {
      printf("Malloc error %d\n",j);
      return;
    }
  }
  left = (double *) malloc((MAXDEG+1)*sizeof(double));
  if (left == NULL) {
    printf("Malloc error\n");
    return;
  }
  right = (double *) malloc((MAXDEG+1)*sizeof(double));
  if (right == NULL) {
    printf("Malloc error\n");
    return;
  }
#endif
  
  ndu[0][0] = 1.0;
  for (j = 1; j <= p; j++) {
    left[j]  = u - knot[i+1-j];
    right[j] = knot[i+j] - u;
    saved = 0.0;
    for (r = 0; r < j; r++) {
      ndu[j][r] = right[r+1] + left[j-r];
      temp      = ndu[r][j-1]/ndu[j][r];
      ndu[r][j] = saved + right[r+1]*temp;
      saved     = left[j-r]*temp;
    }
    ndu[j][j] = saved;
  }
  
  for (j = 0; j <= p; j++) ders[0][j] = ndu[j][p];  /* basis function */
  
  /* compute derivatives */
  for (r = 0; r <= p; r++ ) {
    s1      = 0;
    s2      = 1;
    a[0][0] = 1.0;
    /* compute k'th derivative */
    for (k = 1; k <= der; k++ ) {
      d  = 0.0;
      rk = r - k;
      pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0]/ndu[pk+1][rk];
        d        = a[s2][0]*ndu[rk][pk];
      }
      j1 = rk >= -1  ? 1   : -rk;
      j2 = (r-1<=pk) ? k-1 : p-r;
      for (j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
        d       +=  a[s2][j]*ndu[rk+j][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
        d       +=  a[s2][k]*ndu[r][pk];
      }
      ders[k][r] = d;
      /* switch rows */
      j  = s1;
      s1 = s2;
      s2 = j;
    }
  }
  
  r = p;
  for (k = 1; k <= der; k++) {
    for (j = 0; j <= p; j++ ) ders[k][j] *= r;
    r *= p - k;
  }
  
#ifdef __CUDA_ARCH__
  free(right);
  free(left);
  for (j = 0; j < 2;        ++j) free(a[j]);
  for (j = 0; j < MAXDEG+1; ++j) free(ndu[j]);
  free(ndu);
#endif

}

 
__HOST_AND_DEVICE__ static int
EG_splinePCDeriv(int *ivec, double *data, double t, double *deriv)
{
  int    der = 2;
  int    i, j, k, degree, nKnots, span, dt;
  double x, y, wsum, v[9];  /* note: v is sized for der <= 2! */
#ifndef __CUDA_ARCH__
  double Nders[MAXDEG+1][MAXDEG+1], *Nder[MAXDEG+1], *CP, *w;
#else
  double **Nders, **Nder, *CP, *w;

  Nders = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nders == NULL) {
    printf("Malloc Nders error\n");
    return EGADS_MALLOC;
  }
  for (j = 0; j < MAXDEG+1; ++j) {
    Nders[j] = (double *) malloc((MAXDEG+1)*sizeof(double));
    if (Nders[j] == NULL) {
      printf("Malloc Nders[%d] error\n",j);
      return EGADS_MALLOC;
    }
  }
  Nder = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nder == NULL) {
    printf("Malloc Nder error\n");
    return EGADS_MALLOC;
  }
#endif
  
  for (k = 0; k <= der; k++) deriv[2*k  ] = deriv[2*k+1] = 0.0;
  
  degree = ivec[1];
  nKnots = ivec[3];
  CP     = data + ivec[3];
  w      = CP + 2*ivec[2];
  dt     = MIN(der, degree);
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_splinePCDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree >= MAXDEG) {
    printf(" EG_splinePCDeriv: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }
  for (i = 0; i <= degree; i++) Nder[i] = &Nders[i][0];
  
  span = FindSpan(nKnots, degree, t, data);
  DersBasisFuns(span,     degree, t, data, dt, Nder);
  
  if (ivec[0] == 0) {
    
    for (k = 0; k <= dt; k++)
      for (j = 0; j <= degree; j++) {
        i             = span-degree+j;
        deriv[2*k  ] += Nders[k][j]*CP[2*i  ];
        deriv[2*k+1] += Nders[k][j]*CP[2*i+1];
      }
    
  } else {

    for (k = 0; k <= dt; k++) {
      x = y = wsum = 0.0;
      for (j = 0; j <= degree; j++) {
        i     = span-degree+j;
        x    += Nders[k][j]*w[i]*CP[2*i  ];
        y    += Nders[k][j]*w[i]*CP[2*i+1];
        wsum += Nders[k][j]*w[i];
      }
      v[3*k  ] = x;
      v[3*k+1] = y;
      v[3*k+2] = wsum;
    }
    EG_EvaluateQuotientRule(2, der, 3, v);
    for (k = 0; k <= der; k++) {
      deriv[2*k  ] = v[3*k  ];
      deriv[2*k+1] = v[3*k+1];
    }

  }
  
#ifdef __CUDA_ARCH__
  free(Nder);
  for (j = 0; j < MAXDEG+1; ++j) free(Nders[j]);
  free(Nders);
#endif
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_spline1dDeriv(int *ivec, double *data, double t, double *deriv)
{
  int    der = 2;
  int    i, j, k, degree, nKnots, span, dt;
  double x, y, z, wsum, v[12];  /* note: v is sized for der <= 2! */
#ifndef __CUDA_ARCH__
  double Nders[MAXDEG+1][MAXDEG+1], *Nder[MAXDEG+1], *CP, *w;
#else
  double **Nders, **Nder, *CP, *w;

  Nders = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nders == NULL) {
    printf("Malloc Nders error\n");
    return EGADS_MALLOC;
  }
  for (j = 0; j < MAXDEG+1; ++j) {
    Nders[j] = (double *) malloc((MAXDEG+1)*sizeof(double));
    if (Nders[j] == NULL) {
      printf("Malloc Nders[%d] error\n",j);
      return EGADS_MALLOC;
    }
  }
  Nder = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nder == NULL) {
    printf("Malloc Nder error\n");
    return EGADS_MALLOC;
  }
#endif
  
  for (k = 0; k <= der; k++) deriv[3*k  ] = deriv[3*k+1] = deriv[3*k+2] = 0.0;
  
  degree = ivec[1];
  nKnots = ivec[3];
  CP     = data + ivec[3];
  w      = CP + 3*ivec[2];
  dt     = MIN(der, degree);
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_spline1dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree >= MAXDEG) {
    printf(" EG_spline1dDeriv: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }
  for (i = 0; i <= degree; i++) Nder[i] = &Nders[i][0];
  
  span = FindSpan(nKnots, degree, t, data);
  DersBasisFuns(span,     degree, t, data, dt, Nder);
  
  if (ivec[0] == 0) {
    
    for (k = 0; k <= dt; k++)
      for (j = 0; j <= degree; j++) {
        i             = span-degree+j;
        deriv[3*k  ] += Nders[k][j]*CP[3*i  ];
        deriv[3*k+1] += Nders[k][j]*CP[3*i+1];
        deriv[3*k+2] += Nders[k][j]*CP[3*i+2];
      }
    
  } else {
    
    for (k = 0; k <= dt; k++) {
      x = y = z = wsum = 0.0;
      for (j = 0; j <= degree; j++) {
        i     = span-degree+j;
        x    += Nders[k][j]*w[i]*CP[3*i  ];
        y    += Nders[k][j]*w[i]*CP[3*i+1];
        z    += Nders[k][j]*w[i]*CP[3*i+2];
        wsum += Nders[k][j]*w[i];
      }
      v[4*k  ] = x;
      v[4*k+1] = y;
      v[4*k+2] = z;
      v[4*k+3] = wsum;
    }
    EG_EvaluateQuotientRule(3, der, 4, v);
    for (k = 0; k <= der; k++) {
      deriv[3*k  ] = v[4*k  ];
      deriv[3*k+1] = v[4*k+1];
      deriv[3*k+2] = v[4*k+2];
    }

  }

#ifdef __CUDA_ARCH__
  free(Nder);
  for (j = 0; j < MAXDEG+1; ++j) free(Nders[j]);
  free(Nders);
#endif
  return EGADS_SUCCESS;
}

 
__HOST_AND_DEVICE__ static int
EG_spline2dDeriv(int *ivec, double *data, const double *uv, double *deriv)
{
  int    der = 2;
  int    i, j, k, l, m, s, degu, degv, nKu, nKv, nCPu, spanu, spanv, du, dv;
  double v[24];  /* note: v is sized for der <= 2! */
#ifndef __CUDA_ARCH__
  double *Kv, *CP, *w, *NderU[MAXDEG+1], *NderV[MAXDEG+1];
  double Nu[MAXDEG+1][MAXDEG+1], Nv[MAXDEG+1][MAXDEG+1], temp[4*MAXDEG];
#else
  double *Kv, *CP, *w, **NderU, **NderV;
  double **Nu, **Nv, *temp;

  NderU = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (NderU == NULL) {
    printf("Malloc NderU error\n");
    return EGADS_MALLOC;
  }
  NderV = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (NderV == NULL) {
    printf("Malloc NderV error\n");
    return EGADS_MALLOC;
  }
  Nu = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nu == NULL) {
    printf("Malloc Nu error\n");
    return EGADS_MALLOC;
  }
  for (j = 0; j < MAXDEG+1; ++j) {
    Nu[j] = (double *) malloc((MAXDEG+1)*sizeof(double));
    if (Nu[j] == NULL) {
      printf("Malloc Nu[%d] error\n",j);
      return EGADS_MALLOC;
    }
  }
  Nv = (double **) malloc((MAXDEG+1)*sizeof(double *));
  if (Nv == NULL) {
    printf("Malloc Nv error\n");
    return EGADS_MALLOC;
  }
  for (j = 0; j < MAXDEG+1; ++j) {
    Nv[j] = (double *) malloc((MAXDEG+1)*sizeof(double));
    if (Nv[j] == NULL) {
      printf("Malloc Nv[%d] error\n",j);
      return EGADS_MALLOC;
    }
  }
  temp = (double *) malloc(4*MAXDEG*sizeof(double));
  if (temp == NULL) {
    printf("Malloc temp error\n");
    return EGADS_MALLOC;
  }
#endif
  
  degu = ivec[1];
  nCPu = ivec[2];
  nKu  = ivec[3];
  degv = ivec[4];
  nKv  = ivec[6];
  Kv   = data + ivec[3];
  CP   = data + ivec[3] + ivec[6];
  w    = CP   + 3*ivec[2]*ivec[5];
  du   = MIN(der, degu);
  dv   = MIN(der, degv);
  for (m = l = 0; l <= der; l++)
    for (k = 0; k <= der-l; k++, m++) {
      deriv[3*m  ] = deriv[3*m+1] = deriv[3*m+2] = 0.0;
          v[4*m  ] =     v[4*m+1] =     v[4*m+2] = v[4*m+3] = 0.0;
    }
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_spline2dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degu >= MAXDEG) {
    printf(" EG_spline2dDeriv: degreeU %d >= %d!\n", degu, MAXDEG);
    return EGADS_CONSTERR;
  }
  if (degv >= MAXDEG) {
    printf(" EG_spline2dDeriv: degreeV %d >= %d!\n", degv, MAXDEG);
    return EGADS_CONSTERR;
  }
  for (i = 0; i <= degu; i++) NderU[i] = &Nu[i][0];
  for (i = 0; i <= degv; i++) NderV[i] = &Nv[i][0];
  
  spanu = FindSpan(nKu, degu, uv[0], data);
  DersBasisFuns(spanu,  degu, uv[0], data, du, NderU);
  spanv = FindSpan(nKv, degv, uv[1], Kv);
  DersBasisFuns(spanv,  degv, uv[1], Kv, dv, NderV);

  if (ivec[0] == 0) {
    
    for (m = l = 0; l <= dv; l++)
      for (k = 0; k <= der-l; k++, m++) {
        if (k > du) continue;
        for (s = 0; s <= degv; s++) {
          temp[3*s] = temp[3*s+1] = temp[3*s+2] = 0.0;
          for (j = 0; j <= degu; j++) {
            i            = spanu-degu+j + nCPu*(spanv-degv+s);
            temp[3*s  ] += Nu[k][j]*CP[3*i  ];
            temp[3*s+1] += Nu[k][j]*CP[3*i+1];
            temp[3*s+2] += Nu[k][j]*CP[3*i+2];
          }
        }
        for (s = 0; s <= degv;  s++) {
          deriv[3*m  ] += Nv[l][s]*temp[3*s  ];
          deriv[3*m+1] += Nv[l][s]*temp[3*s+1];
          deriv[3*m+2] += Nv[l][s]*temp[3*s+2];
        }
      }
    /* reorder to match EGADS (der = 2) */
    temp[0]   = deriv[ 6];
    temp[1]   = deriv[ 7];
    temp[2]   = deriv[ 8];
    deriv[ 6] = deriv[ 9];
    deriv[ 7] = deriv[10];
    deriv[ 8] = deriv[11];
    deriv[ 9] = temp[0];
    deriv[10] = temp[1];
    deriv[11] = temp[2];
    
  } else {

    for (m = l = 0; l <= dv; l++)
      for (k = 0; k <= der-l; k++, m++) {
        if (k > du) continue;
        for (s = 0; s <= degv; s++) {
          temp[4*s  ] = temp[4*s+1] = temp[4*s+2] = temp[4*s+3] = 0.0;
          for (j = 0; j <= degu; j++) {
            i            = spanu-degu+j + nCPu*(spanv-degv+s);
            temp[4*s  ] += Nu[k][j]*w[i]*CP[3*i  ];
            temp[4*s+1] += Nu[k][j]*w[i]*CP[3*i+1];
            temp[4*s+2] += Nu[k][j]*w[i]*CP[3*i+2];
            temp[4*s+3] += Nu[k][j]*w[i];
          }
        }
        for (s = 0; s <= degv;  s++) {
          v[4*m  ] += Nv[l][s]*temp[4*s  ];
          v[4*m+1] += Nv[l][s]*temp[4*s+1];
          v[4*m+2] += Nv[l][s]*temp[4*s+2];
          v[4*m+3] += Nv[l][s]*temp[4*s+3];
        }
      }
    /* reorder to match EGADS (der = 2) */
    temp[0] = v[ 8];
    temp[1] = v[ 9];
    temp[2] = v[10];
    temp[3] = v[11];
    v[ 8]   = v[12];
    v[ 9]   = v[13];
    v[10]   = v[14];
    v[11]   = v[15];
    v[12]   = temp[0];
    v[13]   = temp[1];
    v[14]   = temp[2];
    v[15]   = temp[3];
    EG_EvaluateQuotientRule2(3, der, 4, v);
    for (m = l = 0; l <= der; l++)
      for (k = 0; k <= der-l; k++, m++) {
        deriv[3*m  ] = v[4*m  ];
        deriv[3*m+1] = v[4*m+1];
        deriv[3*m+2] = v[4*m+2];
      }

  }
  
#ifdef __CUDA_ARCH__
  free(temp);
  for (j = 0; j < MAXDEG+1; ++j) free(NderV[j]);
  free(NderV);
  for (j = 0; j < MAXDEG+1; ++j) free(NderU[j]);
  free(NderU);
  free(Nv);
  free(Nu);
#endif
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_bezierPCDeriv(int *ivec, double *data, double t, double *deriv)
{
  bool   stat;
  int    i, k, degree;
#ifndef __CUDA_ARCH__
  double Q[2*MAXDEG+2], S[2*MAXDEG], T[2*MAXDEG-2], CP[3*MAXDEG+3], *w;
#else
  double *Q, *S, *T, *CP, *w;

  Q = (double *) malloc((2*MAXDEG+2)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return EGADS_MALLOC;
  }
  S = (double *) malloc(2*MAXDEG*sizeof(double));
  if (S == NULL) {
    printf("Malloc S error\n");
    return EGADS_MALLOC;
  }
  T = (double *) malloc((2*MAXDEG-2)*sizeof(double));
  if (T == NULL) {
    printf("Malloc T error\n");
    return EGADS_MALLOC;
  }
  CP = (double *) malloc((3*MAXDEG+3)*sizeof(double));
  if (CP == NULL) {
    printf("Malloc CP error\n");
    return EGADS_MALLOC;
  }
#endif
  
  for (k = 0; k < 6; k++) deriv[k] = 0.0;
  
  degree = ivec[2]-1;
  w      = data + 2*ivec[2];
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_bezierPCDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree <= 0) {
    printf(" EG_bezierPCDeriv: degree %d <= 0!\n", degree);
    return EGADS_CONSTERR;
  }

  if (degree >= MAXDEG) {
    printf(" EG_bezierPCDeriv: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }
  
  if (ivec[0] == 0) {
    
    for (i = 0; i <= degree; i++) {
      Q[2*i  ] = data[2*i  ];
      Q[2*i+1] = data[2*i+1];
    }
    
    for (k = 1; k <= degree; k++)
      for (i = 0; i <= degree-k; i++) {
        Q[2*i  ] = (1.0-t)*Q[2*i  ] + t*Q[2*i+2];
        Q[2*i+1] = (1.0-t)*Q[2*i+1] + t*Q[2*i+3];
      }
    deriv[0] = Q[0];
    deriv[1] = Q[1];
    
    for (i = 0; i < degree; i++) {
      Q[2*i  ] = data[2*i  ];
      Q[2*i+1] = data[2*i+1];
      S[2*i  ] = data[2*i+2];
      S[2*i+1] = data[2*i+3];
    }
    
    for (k = 1; k < degree; k++)
      for (i = 0; i < degree-k; i++) {
        Q[2*i  ] = (1.0-t)*Q[2*i  ] + t*Q[2*i+2];
        Q[2*i+1] = (1.0-t)*Q[2*i+1] + t*Q[2*i+3];
        S[2*i  ] = (1.0-t)*S[2*i  ] + t*S[2*i+2];
        S[2*i+1] = (1.0-t)*S[2*i+1] + t*S[2*i+3];
      }
    deriv[2] = degree*(S[0] - Q[0]);
    deriv[3] = degree*(S[1] - Q[1]);
    
    if (degree > 1) {
      for (i = 0; i < degree-1; i++) {
        Q[2*i  ] = data[2*i  ];
        Q[2*i+1] = data[2*i+1];
        S[2*i  ] = data[2*i+2];
        S[2*i+1] = data[2*i+3];
        T[2*i  ] = data[2*i+4];
        T[2*i+1] = data[2*i+5];
      }
      
      for (k = 1; k < degree-1; k++)
        for (i = 0; i < degree-k-1; i++) {
          Q[2*i  ] = (1.0-t)*Q[2*i  ] + t*Q[2*i+2];
          Q[2*i+1] = (1.0-t)*Q[2*i+1] + t*Q[2*i+3];
          S[2*i  ] = (1.0-t)*S[2*i  ] + t*S[2*i+2];
          S[2*i+1] = (1.0-t)*S[2*i+1] + t*S[2*i+3];
          T[2*i  ] = (1.0-t)*T[2*i  ] + t*T[2*i+2];
          T[2*i+1] = (1.0-t)*T[2*i+1] + t*T[2*i+3];
        }
      deriv[4] = degree*(degree-1)*(T[0] - 2.0*S[0] + Q[0]);
      deriv[5] = degree*(degree-1)*(T[1] - 2.0*S[1] + Q[1]);
    }
    
  } else {
    
    for (i = 0; i <= degree; i++) {
      CP[3*i  ] = data[2*i  ]*w[i];
      CP[3*i+1] = data[2*i+1]*w[i];
      CP[3*i+2] = w[i];
    }
    stat = EG_Bezier1DRat(2, ivec[2], 3, CP, 0.0, 1.0, 2, t, 2, deriv);
    if (!stat) return EGADS_GEOMERR;
    
  }
  
#ifdef __CUDA_ARCH__
  free(CP);
  free(T);
  free(S);
  free(Q);
#endif
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static void
deCasteljau(double *P, int n, double u, double *C)
{
  int    i, k;
#ifndef __CUDA_ARCH__
  double Q[3*MAXDEG+3];
#else
  double *Q;

  Q = (double *) malloc((3*MAXDEG+3)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return;
  }
#endif
  
  Q[0] = Q[1] = Q[2] = 0.0;
  for (i = 0; i <= n; i++) {
    Q[3*i  ] = P[3*i  ];
    Q[3*i+1] = P[3*i+1];
    Q[3*i+2] = P[3*i+2];
  }
  
  for (k = 1; k <= n; k++)
    for (i = 0; i <= n-k; i++) {
      Q[3*i  ] = (1.0-u)*Q[3*i  ] + u*Q[3*i+3];
      Q[3*i+1] = (1.0-u)*Q[3*i+1] + u*Q[3*i+4];
      Q[3*i+2] = (1.0-u)*Q[3*i+2] + u*Q[3*i+5];
    }
  
  C[0] = Q[0];
  C[1] = Q[1];
  C[2] = Q[2];

#ifdef __CUDA_ARCH__
  free(Q);
#endif
}


__HOST_AND_DEVICE__ static void
deCasteljauD1(double *P, int n, double u, double *C)
{
  int    i, k;
#ifndef __CUDA_ARCH__
  double Q[3*MAXDEG], S[3*MAXDEG];
#else
  double *Q, *S;

  Q = (double *) malloc((3*MAXDEG)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return;
  }
  S = (double *) malloc((3*MAXDEG)*sizeof(double));
  if (S == NULL) {
    printf("Malloc S error\n");
    return;
  }
#endif
  
  Q[0] = Q[1] = Q[2] = 0.0;
  S[0] = S[1] = S[2] = 0.0;
  for (i = 0; i < n; i++) {
    Q[3*i  ] = P[3*i  ];
    Q[3*i+1] = P[3*i+1];
    Q[3*i+2] = P[3*i+2];
    S[3*i  ] = P[3*i+3];
    S[3*i+1] = P[3*i+4];
    S[3*i+2] = P[3*i+5];
  }
  
  for (k = 1; k < n; k++)
    for (i = 0; i < n-k; i++) {
      Q[3*i  ] = (1.0-u)*Q[3*i  ] + u*Q[3*i+3];
      Q[3*i+1] = (1.0-u)*Q[3*i+1] + u*Q[3*i+4];
      Q[3*i+2] = (1.0-u)*Q[3*i+2] + u*Q[3*i+5];
      S[3*i  ] = (1.0-u)*S[3*i  ] + u*S[3*i+3];
      S[3*i+1] = (1.0-u)*S[3*i+1] + u*S[3*i+4];
      S[3*i+2] = (1.0-u)*S[3*i+2] + u*S[3*i+5];
    }
  
  C[0] = n*(S[0] - Q[0]);
  C[1] = n*(S[1] - Q[1]);
  C[2] = n*(S[2] - Q[2]);

#ifdef __CUDA_ARCH__
  free(S);
  free(Q);
#endif
}


__HOST_AND_DEVICE__ static void
deCasteljauD2(double *P, int n, double u, double *C)
{
  int    i, k;
#ifndef __CUDA_ARCH__
  double Q[3*MAXDEG-3], S[3*MAXDEG-3], T[3*MAXDEG-3];
#else
  double *Q, *S, *T;

  Q = (double *) malloc((3*MAXDEG-3)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return;
  }
  S = (double *) malloc((3*MAXDEG-3)*sizeof(double));
  if (S == NULL) {
    printf("Malloc S error\n");
    return;
  }
  T = (double *) malloc((3*MAXDEG-3)*sizeof(double));
  if (T == NULL) {
    printf("Malloc T error\n");
    return;
  }
#endif
  
  C[0] = C[1] = C[2] = 0.0;
  if (n <= 1) return;
  
  for (i = 0; i < n-1; i++) {
    Q[3*i  ] = P[3*i  ];
    Q[3*i+1] = P[3*i+1];
    Q[3*i+2] = P[3*i+2];
    S[3*i  ] = P[3*i+3];
    S[3*i+1] = P[3*i+4];
    S[3*i+2] = P[3*i+5];
    T[3*i  ] = P[3*i+6];
    T[3*i+1] = P[3*i+7];
    T[3*i+2] = P[3*i+8];
  }
  
  for (k = 1; k < n-1; k++)
    for (i = 0; i < n-k-1; i++) {
      Q[3*i  ] = (1.0-u)*Q[3*i  ] + u*Q[3*i+3];
      Q[3*i+1] = (1.0-u)*Q[3*i+1] + u*Q[3*i+4];
      Q[3*i+2] = (1.0-u)*Q[3*i+2] + u*Q[3*i+5];
      S[3*i  ] = (1.0-u)*S[3*i  ] + u*S[3*i+3];
      S[3*i+1] = (1.0-u)*S[3*i+1] + u*S[3*i+4];
      S[3*i+2] = (1.0-u)*S[3*i+2] + u*S[3*i+5];
      T[3*i  ] = (1.0-u)*T[3*i  ] + u*T[3*i+3];
      T[3*i+1] = (1.0-u)*T[3*i+1] + u*T[3*i+4];
      T[3*i+2] = (1.0-u)*T[3*i+2] + u*T[3*i+5];
    }
  
  C[0] = n*(n-1)*(T[0] - 2.0*S[0] + Q[0]);
  C[1] = n*(n-1)*(T[1] - 2.0*S[1] + Q[1]);
  C[2] = n*(n-1)*(T[2] - 2.0*S[2] + Q[2]);

#ifdef __CUDA_ARCH__
  free(T);
  free(S);
  free(Q);
#endif
}


__HOST_AND_DEVICE__ static int
EG_bezier1dDeriv(int *ivec, double *data, double t, double *deriv)
{
  bool   stat;
  int    i, k, degree;
#ifndef __CUDA_ARCH__
  double Q[3*MAXDEG+3], S[3*MAXDEG], T[3*MAXDEG-3], CP[4*MAXDEG+4], *w;
#else
  double *Q, *S, *T, *CP, *w;

  Q = (double *) malloc((3*MAXDEG+3)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return EGADS_MALLOC;
  }
  S = (double *) malloc((3*MAXDEG)*sizeof(double));
  if (S == NULL) {
    printf("Malloc S error\n");
    return EGADS_MALLOC;
  }
  T = (double *) malloc((3*MAXDEG-3)*sizeof(double));
  if (T == NULL) {
    printf("Malloc T error\n");
    return EGADS_MALLOC;
  }
  CP = (double *) malloc((4*MAXDEG+4)*sizeof(double));
  if (CP == NULL) {
    printf("Malloc CP error\n");
    return EGADS_MALLOC;
  }
#endif
  
  for (k = 0; k < 9; k++) deriv[k] = 0.0;
  
  degree = ivec[2]-1;
  w      = data + 3*ivec[2];
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_bezier1dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if (degree <= 0) {
    printf(" EG_bezier1dDeriv: degree %d <= 0!\n", degree);
    return EGADS_CONSTERR;
  }
  if (degree >= MAXDEG) {
    printf(" EG_bezier1dDeriv: degree %d >= %d!\n", degree, MAXDEG);
    return EGADS_CONSTERR;
  }
  
  if (ivec[0] == 0) {
    
    for (i = 0; i <= degree; i++) {
      Q[3*i  ] = data[3*i  ];
      Q[3*i+1] = data[3*i+1];
      Q[3*i+2] = data[3*i+2];
    }
    
    for (k = 1; k <= degree; k++)
      for (i = 0; i <= degree-k; i++) {
        Q[3*i  ] = (1.0-t)*Q[3*i  ] + t*Q[3*i+3];
        Q[3*i+1] = (1.0-t)*Q[3*i+1] + t*Q[3*i+4];
        Q[3*i+2] = (1.0-t)*Q[3*i+2] + t*Q[3*i+5];
      }
    deriv[0] = Q[0];
    deriv[1] = Q[1];
    deriv[2] = Q[2];
    
    for (i = 0; i < degree; i++) {
      Q[3*i  ] = data[3*i  ];
      Q[3*i+1] = data[3*i+1];
      Q[3*i+2] = data[3*i+2];
      S[3*i  ] = data[3*i+3];
      S[3*i+1] = data[3*i+4];
      S[3*i+2] = data[3*i+5];
    }
    
    for (k = 1; k < degree; k++)
      for (i = 0; i < degree-k; i++) {
        Q[3*i  ] = (1.0-t)*Q[3*i  ] + t*Q[3*i+3];
        Q[3*i+1] = (1.0-t)*Q[3*i+1] + t*Q[3*i+4];
        Q[3*i+2] = (1.0-t)*Q[3*i+2] + t*Q[3*i+5];
        S[3*i  ] = (1.0-t)*S[3*i  ] + t*S[3*i+3];
        S[3*i+1] = (1.0-t)*S[3*i+1] + t*S[3*i+4];
        S[3*i+2] = (1.0-t)*S[3*i+2] + t*S[3*i+5];
      }
    deriv[3] = degree*(S[0] - Q[0]);
    deriv[4] = degree*(S[1] - Q[1]);
    deriv[5] = degree*(S[2] - Q[2]);
    
    if (degree > 1) {
      for (i = 0; i < degree-1; i++) {
        Q[3*i  ] = data[3*i  ];
        Q[3*i+1] = data[3*i+1];
        Q[3*i+2] = data[3*i+2];
        S[3*i  ] = data[3*i+3];
        S[3*i+1] = data[3*i+4];
        S[3*i+2] = data[3*i+5];
        T[3*i  ] = data[3*i+6];
        T[3*i+1] = data[3*i+7];
        T[3*i+2] = data[3*i+8];
      }
      
      for (k = 1; k < degree-1; k++)
        for (i = 0; i < degree-k-1; i++) {
          Q[3*i  ] = (1.0-t)*Q[3*i  ] + t*Q[3*i+3];
          Q[3*i+1] = (1.0-t)*Q[3*i+1] + t*Q[3*i+4];
          Q[3*i+2] = (1.0-t)*Q[3*i+2] + t*Q[3*i+5];
          S[3*i  ] = (1.0-t)*S[3*i  ] + t*S[3*i+3];
          S[3*i+1] = (1.0-t)*S[3*i+1] + t*S[3*i+4];
          S[3*i+2] = (1.0-t)*S[3*i+2] + t*S[3*i+5];
          T[3*i  ] = (1.0-t)*T[3*i  ] + t*T[3*i+3];
          T[3*i+1] = (1.0-t)*T[3*i+1] + t*T[3*i+4];
          T[3*i+2] = (1.0-t)*T[3*i+2] + t*T[3*i+5];
        }
      deriv[6] = degree*(degree-1)*(T[0] - 2.0*S[0] + Q[0]);
      deriv[7] = degree*(degree-1)*(T[1] - 2.0*S[1] + Q[1]);
      deriv[8] = degree*(degree-1)*(T[2] - 2.0*S[2] + Q[2]);
    }
    
  } else {

    for (i = 0; i <= degree; i++) {
      CP[4*i  ] = data[3*i  ]*w[i];
      CP[4*i+1] = data[3*i+1]*w[i];
      CP[4*i+2] = data[3*i+2]*w[i];
      CP[4*i+3] = w[i];
    }
    stat = EG_Bezier1DRat(3, ivec[2], 4, CP, 0.0, 1.0, 2, t, 3, deriv);
    if (!stat) return EGADS_GEOMERR;

  }
  
#ifdef __CUDA_ARCH__
  free(CP);
  free(T);
  free(S);
  free(Q);
#endif
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_bezier2dDeriv(int *ivec, double *data, const double *uv, double *deriv)
{
  int    i, j, n, m, header[7];
#ifndef __CUDA_ARCH__
  double Q[3*MAXDEG+3], P[3*MAXDEG+3], D[2*MAXDEG+4+4*(MAXDEG+1)*(MAXDEG+1)];
#else
  double *Q, *P, *D;

  Q = (double *) malloc((3*MAXDEG+3)*sizeof(double));
  if (Q == NULL) {
    printf("Malloc Q error\n");
    return EGADS_MALLOC;
  }
  P = (double *) malloc((3*MAXDEG+3)*sizeof(double));
  if (P == NULL) {
    printf("Malloc P error\n");
    return EGADS_MALLOC;
  }
  D = (double *) malloc((2*MAXDEG+4+4*(MAXDEG+1)*(MAXDEG+1))*sizeof(double));
  if (D == NULL) {
    printf("Malloc D error\n");
    return EGADS_MALLOC;
  }
#endif
  
  for (j = 0; j < 18; j++) deriv[j] = 0.0;
  
  n = ivec[2]-1;
  m = ivec[4]-1;
  if ((ivec[0] != 0) && (ivec[0] != 2)) {
    printf(" EG_bezier2dDeriv: flag = %d!\n", ivec[0]);
    return EGADS_GEOMERR;
  }
  if ((n >= MAXDEG) || (m >= MAXDEG)) {
    printf(" EG_bezier2dDeriv: degree %d or %d >= %d!\n", n, m, MAXDEG);
    return EGADS_CONSTERR;
  }
  
  if (ivec[0] == 0) {
    
    if (n <= m) {
      
      for (j = 0; j <= m; j++)
        deCasteljau(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljau(Q, m, uv[1], deriv);
      
      for (j = 0; j <= m; j++)
        deCasteljauD1(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljau(Q, m, uv[1], &deriv[3]);
      
      for (j = 0; j <= m; j++)
        deCasteljau(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljauD1(Q, m, uv[1], &deriv[6]);

      for (j = 0; j <= m; j++)
        deCasteljauD2(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljau(Q, m, uv[1], &deriv[9]);
      
      for (j = 0; j <= m; j++)
        deCasteljauD1(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljauD1(Q, m, uv[1], &deriv[12]);
      
      for (j = 0; j <= m; j++)
        deCasteljau(&data[3*j*ivec[2]], n, uv[0], &Q[3*j]);
      deCasteljauD2(Q, m, uv[1], &deriv[15]);
      
    } else {
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljau(P, m, uv[1], &Q[3*i]);
      }
      deCasteljau(Q, n, uv[0], deriv);
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljau(P, m, uv[1], &Q[3*i]);
      }
      deCasteljauD1(Q, n, uv[0], &deriv[3]);
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljauD1(P, m, uv[1], &Q[3*i]);
      }
      deCasteljau(Q, n, uv[0], &deriv[6]);
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljau(P, m, uv[1], &Q[3*i]);
      }
      deCasteljauD2(Q, n, uv[0], &deriv[9]);
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljauD1(P, m, uv[1], &Q[3*i]);
      }
      deCasteljauD1(Q, n, uv[0], &deriv[12]);
      
      for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
          P[3*j  ] = data[3*(i+j*ivec[2])  ];
          P[3*j+1] = data[3*(i+j*ivec[2])+1];
          P[3*j+2] = data[3*(i+j*ivec[2])+2];
        }
        deCasteljauD2(P, m, uv[1], &Q[3*i]);
      }
      deCasteljau(Q, n, uv[0], &deriv[15]);
    }
    
  } else {
    
    /* convert to BSpline and evaluate */
    header[0] =   ivec[0];
    header[1] =   ivec[1];
    header[2] =   ivec[2];
    header[3] = 2*ivec[2];
    header[4] =   ivec[3];
    header[5] =   ivec[4];
    header[6] = 2*ivec[4];
    for (n = i = 0; i < ivec[2]; i++, n++) D[n] = 0.0;
    for (    i = 0; i < ivec[2]; i++, n++) D[n] = 1.0;
    for (    i = 0; i < ivec[4]; i++, n++) D[n] = 0.0;
    for (    i = 0; i < ivec[4]; i++, n++) D[n] = 1.0;
    m = 3*ivec[2]*ivec[4];
    if ((ivec[0]&2) != 0) m += ivec[2]*ivec[4];
    for (    i = 0; i < m; i++, n++) D[n] = data[i];
    return EG_spline2dDeriv(header, D, uv, deriv);
    
  }
  
#ifdef __CUDA_ARCH__
  free(D);
  free(P);
  free(P);
#endif
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static void
EG_rotatePC(double *axes, double *data, double *result)
{
  result[0] = data[0]*axes[2] + data[1]*axes[4] + axes[0];
  result[1] = data[0]*axes[3] + data[1]*axes[5] + axes[1];
  result[2] = data[2]*axes[2] + data[3]*axes[4];
  result[3] = data[2]*axes[3] + data[3]*axes[5];
  result[4] = data[4]*axes[2] + data[5]*axes[4];
  result[5] = data[4]*axes[3] + data[5]*axes[5];
}


__HOST_AND_DEVICE__ static void
EG_rotate2D(double *axes, double *data, double *result)
{
  double zaxis[3], *xaxis, *yaxis;
  
  xaxis = &axes[3];
  yaxis = &axes[6];
  CROSS(zaxis, xaxis, yaxis);
  result[0] = data[0]*axes[3] + data[1]*axes[6] + data[2]*zaxis[0] + axes[0];
  result[1] = data[0]*axes[4] + data[1]*axes[7] + data[2]*zaxis[1] + axes[1];
  result[2] = data[0]*axes[5] + data[1]*axes[8] + data[2]*zaxis[2] + axes[2];
  result[3] = data[3]*axes[3] + data[4]*axes[6] + data[5]*zaxis[0];
  result[4] = data[3]*axes[4] + data[4]*axes[7] + data[5]*zaxis[1];
  result[5] = data[3]*axes[5] + data[4]*axes[8] + data[5]*zaxis[2];
  result[6] = data[6]*axes[3] + data[7]*axes[6] + data[8]*zaxis[0];
  result[7] = data[6]*axes[4] + data[7]*axes[7] + data[8]*zaxis[1];
  result[8] = data[6]*axes[5] + data[7]*axes[8] + data[8]*zaxis[2];
}


__HOST_AND_DEVICE__ static void
EG_rotate3D(double *axes, double *data, double *result)
{
  result[ 0] = data[ 0]*axes[3] + data[ 1]*axes[6] + data[ 2]*axes[ 9] + axes[0];
  result[ 1] = data[ 0]*axes[4] + data[ 1]*axes[7] + data[ 2]*axes[10] + axes[1];
  result[ 2] = data[ 0]*axes[5] + data[ 1]*axes[8] + data[ 2]*axes[11] + axes[2];
  result[ 3] = data[ 3]*axes[3] + data[ 4]*axes[6] + data[ 5]*axes[ 9];
  result[ 4] = data[ 3]*axes[4] + data[ 4]*axes[7] + data[ 5]*axes[10];
  result[ 5] = data[ 3]*axes[5] + data[ 4]*axes[8] + data[ 5]*axes[11];
  result[ 6] = data[ 6]*axes[3] + data[ 7]*axes[6] + data[ 8]*axes[ 9];
  result[ 7] = data[ 6]*axes[4] + data[ 7]*axes[7] + data[ 8]*axes[10];
  result[ 8] = data[ 6]*axes[5] + data[ 7]*axes[8] + data[ 8]*axes[11];
  result[ 9] = data[ 9]*axes[3] + data[10]*axes[6] + data[11]*axes[ 9];
  result[10] = data[ 9]*axes[4] + data[10]*axes[7] + data[11]*axes[10];
  result[11] = data[ 9]*axes[5] + data[10]*axes[8] + data[11]*axes[11];
  result[12] = data[12]*axes[3] + data[13]*axes[6] + data[14]*axes[ 9];
  result[13] = data[12]*axes[4] + data[13]*axes[7] + data[14]*axes[10];
  result[14] = data[12]*axes[5] + data[13]*axes[8] + data[14]*axes[11];
  result[15] = data[15]*axes[3] + data[16]*axes[6] + data[17]*axes[ 9];
  result[16] = data[15]*axes[4] + data[16]*axes[7] + data[17]*axes[10];
  result[17] = data[15]*axes[5] + data[16]*axes[8] + data[17]*axes[11];
}


__HOST_AND_DEVICE__ static int
EG_evaluateGeomX(const egObject *geomx, const double *param, double *result)
{
  int            i, stat;
  double         data[18], norm[3], axes[15], dist, radius;
  double         *tan, *tanv, *gdata;
  liteGeometry   *lgeom;
  const egObject *geom;

  geom = geomx;

recurse:
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  
  stat = EGADS_NOTFOUND;

  if (geom->oclass == PCURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    switch (geom->mtype) {
      case LINE:
        result[0] = *param*gdata[2] + gdata[0];
        result[1] = *param*gdata[3] + gdata[1];
        result[2] = gdata[2];
        result[3] = gdata[3];
        result[4] = result[5] = 0.0;
        stat      = EGADS_SUCCESS;
        break;
        
      case CIRCLE:
        data[0] =  gdata[6]*cos(*param);
        data[1] =  gdata[6]*sin(*param);
        data[2] = -gdata[6]*sin(*param);
        data[3] =  gdata[6]*cos(*param);
        data[4] = -gdata[6]*cos(*param);
        data[5] = -gdata[6]*sin(*param);
        EG_rotatePC(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case ELLIPSE:
        data[0] =  gdata[6]*cos(*param);
        data[1] =  gdata[7]*sin(*param);
        data[2] = -gdata[6]*sin(*param);
        data[3] =  gdata[7]*cos(*param);
        data[4] = -gdata[6]*cos(*param);
        data[5] = -gdata[7]*sin(*param);
        EG_rotatePC(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case PARABOLA:
        if (gdata[6] == 0.0) return EGADS_GEOMERR;
        data[1] = *param;
        data[0] = 0.25*data[1]*data[1]/gdata[6];
        data[2] = 0.5*data[1]/gdata[6];
        data[3] = 1.0;
        data[4] = 0.5/gdata[6];
        data[5] = 0.0;
        EG_rotatePC(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case HYPERBOLA:
        data[0] = gdata[6]*cosh(*param);
        data[1] = gdata[7]*sinh(*param);
        data[2] = gdata[6]*sinh(*param);
        data[3] = gdata[7]*cosh(*param);
        data[4] = gdata[6]*cosh(*param);
        data[5] = gdata[7]*sinh(*param);
        EG_rotatePC(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case TRIMMED:
        geom = lgeom->ref;
        goto recurse;
        
      case BEZIER:
        stat = EG_bezierPCDeriv(lgeom->header, gdata, param[0], result);
        break;
        
      case BSPLINE:
        stat = EG_splinePCDeriv(lgeom->header, gdata, param[0], result);
        break;
        
      case OFFSET:
        /* should be handled before getting here! */
        stat = EGADS_REFERCE;
        break;
    }
    
  } else if (geom->oclass == CURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    switch (geom->mtype) {
      case LINE:
        result[0] = *param*gdata[3] + gdata[0];
        result[1] = *param*gdata[4] + gdata[1];
        result[2] = *param*gdata[5] + gdata[2];
        result[3] = gdata[3];
        result[4] = gdata[4];
        result[5] = gdata[5];
        result[6] = result[7] = result[8] =  0.0;
        stat      = EGADS_SUCCESS;
        break;
        
      case CIRCLE:
        data[0] =  gdata[9]*cos(*param);
        data[1] =  gdata[9]*sin(*param);
        data[2] =  0.0;
        data[3] = -gdata[9]*sin(*param);
        data[4] =  gdata[9]*cos(*param);
        data[5] =  0.0;
        data[6] = -gdata[9]*cos(*param);
        data[7] = -gdata[9]*sin(*param);
        data[8] =  0.0;
        EG_rotate2D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case ELLIPSE:
        data[0] =  gdata[ 9]*cos(*param);
        data[1] =  gdata[10]*sin(*param);
        data[2] =  0.0;
        data[3] = -gdata[ 9]*sin(*param);
        data[4] =  gdata[10]*cos(*param);
        data[5] =  0.0;
        data[6] = -gdata[ 9]*cos(*param);
        data[7] = -gdata[10]*sin(*param);
        data[8] =  0.0;
        EG_rotate2D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case PARABOLA:
        if (gdata[9] == 0.0) return EGADS_GEOMERR;
        data[1] = *param;
        data[0] = 0.25*data[1]*data[1]/gdata[9];
        data[2] = 0.0;
        data[3] = 0.5*data[1]/gdata[9];
        data[4] = 1.0;
        data[5] = 0.0;
        data[6] = 0.5/gdata[9];
        data[7] = 0.0;
        data[8] = 0.0;
        EG_rotate2D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case HYPERBOLA:
        data[0] = gdata[ 9]*cosh(*param);
        data[1] = gdata[10]*sinh(*param);
        data[2] = 0.0;
        data[3] = gdata[ 9]*sinh(*param);
        data[4] = gdata[10]*cosh(*param);
        data[5] = 0.0;
        data[6] = gdata[ 9]*cosh(*param);
        data[7] = gdata[10]*sinh(*param);
        data[8] = 0.0;
        EG_rotate2D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case TRIMMED:
        geom = lgeom->ref;
        goto recurse;
        
      case BEZIER:
        stat = EG_bezier1dDeriv(lgeom->header, gdata, param[0], result);
        break;
        
      case BSPLINE:
        stat = EG_spline1dDeriv(lgeom->header, gdata, param[0], result);
        break;
        
      case OFFSET:
        /* should be handled before getting here! */
        stat = EGADS_REFERCE;
        break;
    }
    
  } else {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    switch (geom->mtype) {
      case PLANE:
        result[ 0] = param[0]*gdata[3] + param[1]*gdata[6] + gdata[0];
        result[ 1] = param[0]*gdata[4] + param[1]*gdata[7] + gdata[1];
        result[ 2] = param[0]*gdata[5] + param[1]*gdata[8] + gdata[2];
        result[ 3] = gdata[3];
        result[ 4] = gdata[4];
        result[ 5] = gdata[5];
        result[ 6] = gdata[6];
        result[ 7] = gdata[7];
        result[ 8] = gdata[8];
        result[ 9] = result[10] = result[11] = 0.0;
        result[12] = result[13] = result[14] = 0.0;
        result[15] = result[16] = result[17] = 0.0;
        stat      = EGADS_SUCCESS;
        break;
        
      case SPHERICAL:
        radius =  gdata[9];
        tan    = &gdata[3];
        tanv   = &gdata[6];
        CROSS(norm, tan, tanv);
        dist   = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
        if (gdata[9] < 0.0) {
          dist   = -dist;
          radius = -radius;
        }
        if (dist != 0.0) {
          norm[0] /= dist;
          norm[1] /= dist;
          norm[2] /= dist;
        }
        for (i = 0; i < 9; i++) axes[i]   = gdata[i];
        for (i = 0; i < 3; i++) axes[i+9] = norm[i];
        data[ 0] =  radius*cos(param[0])*cos(param[1]);
        data[ 1] =  radius*sin(param[0])*cos(param[1]);
        data[ 2] =  radius*sin(param[1]);
        data[ 3] = -radius*sin(param[0])*cos(param[1]);
        data[ 4] =  radius*cos(param[0])*cos(param[1]);
        data[ 5] =  0.0;
        data[ 6] = -radius*cos(param[0])*sin(param[1]);
        data[ 7] = -radius*sin(param[0])*sin(param[1]);
        data[ 8] =  radius*cos(param[1]);
        data[ 9] = -radius*cos(param[0])*cos(param[1]);
        data[10] = -radius*sin(param[0])*cos(param[1]);
        data[11] =  0.0;
        data[12] =  radius*sin(param[0])*sin(param[1]);
        data[13] = -radius*cos(param[0])*sin(param[1]);
        data[14] =  0.0;
        data[15] = -radius*cos(param[0])*cos(param[1]);
        data[16] = -radius*sin(param[0])*cos(param[1]);
        data[17] = -radius*sin(param[1]);
        EG_rotate3D(axes, data, result);
        stat     = EGADS_SUCCESS;
        break;
        
      case CONICAL:
        radius   =  gdata[13] + param[1]*sin(gdata[12]);
        data[ 0] =  radius*cos(param[0]);
        data[ 1] =  radius*sin(param[0]);
        data[ 2] =  param[1]*cos(gdata[12]);
        data[ 3] = -radius*sin(param[0]);
        data[ 4] =  radius*cos(param[0]);
        data[ 5] =  0.0;
        data[ 6] =  sin(gdata[12])*cos(param[0]);
        data[ 7] =  sin(gdata[12])*sin(param[0]);
        data[ 8] =  cos(gdata[12]);
        data[ 9] = -radius*cos(param[0]);
        data[10] = -radius*sin(param[0]);
        data[11] =  0.0;
        data[12] = -sin(gdata[12])*sin(param[0]);
        data[13] =  sin(gdata[12])*cos(param[0]);
        data[14] =  0.0;
        data[15] = data[16] = data[17] = 0.0;
        EG_rotate3D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case CYLINDRICAL:
        data[ 0] =  gdata[12]*cos(param[0]);
        data[ 1] =  gdata[12]*sin(param[0]);
        data[ 2] =  param[1];
        data[ 3] = -gdata[12]*sin(param[0]);
        data[ 4] =  gdata[12]*cos(param[0]);
        data[ 5] =  0.0;
        data[ 6] =  0.0;
        data[ 7] =  0.0;
        data[ 8] =  1.0;
        data[ 9] = -gdata[12]*cos(param[0]);
        data[10] = -gdata[12]*sin(param[0]);
        data[11] =  0.0;
        data[12] = data[13] = data[14] = 0.0;
        data[15] = data[16] = data[17] = 0.0;
        EG_rotate3D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case TOROIDAL:
        radius   =  gdata[12] + gdata[13]*cos(param[1]);
        data[ 0] =  radius*cos(param[0]);
        data[ 1] =  radius*sin(param[0]);
        data[ 2] =  gdata[13]*sin(param[1]);
        data[ 3] = -radius*sin(param[0]);
        data[ 4] =  radius*cos(param[0]);
        data[ 5] =  0.0;
        data[ 6] =  -gdata[13]*sin(param[1])*cos(param[0]);
        data[ 7] =  -gdata[13]*sin(param[1])*sin(param[0]);
        data[ 8] =   gdata[13]*cos(param[1]);
        data[ 9] = -radius*cos(param[0]);
        data[10] = -radius*sin(param[0]);
        data[11] =  0.0;
        data[12] =  gdata[13]*sin(param[1])*sin(param[0]);
        data[13] = -gdata[13]*sin(param[1])*cos(param[0]);
        data[14] =  0.0;
        data[15] = -gdata[13]*cos(param[1])*cos(param[0]);
        data[16] = -gdata[13]*cos(param[1])*sin(param[0]);
        data[17] = -gdata[13]*sin(param[1]);
        EG_rotate3D(gdata, data, result);
        stat    = EGADS_SUCCESS;
        break;
        
      case REVOLUTION:
        /* should be handled before getting here! */
        stat = EGADS_REFERCE;
        break;
        
      case EXTRUSION:
        /* should be handled before getting here! */
        stat = EGADS_REFERCE;
        break;
        
      case TRIMMED:
        geom = lgeom->ref;
        goto recurse;
        
      case BEZIER:
        stat = EG_bezier2dDeriv(lgeom->header, gdata, param, result);
        break;
        
      case BSPLINE:
        stat = EG_spline2dDeriv(lgeom->header, gdata, param, result);
        break;
        
      case OFFSET:
        /* should be handled before getting here! */
        stat = EGADS_REFERCE;
        break;
    }
    
  }
  
  return stat;
}


/*
 * copy of EG_evaluateGeom & EG_eval3deriv
 *     used by EG_eval3deriv to avoid recursion
 */

__HOST_AND_DEVICE__ static int
EG_eval3derivZ(const egObject *geom, const double *param, double *result)
{
  int    stat;
  double step = 1.e-7;
  double uv[2], uminus[18], uplus[18], vminus[18], vplus[18];
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;

  uv[1] = 0.0;
  if (geom->oclass == PCURVE) {
    
    uv[0] = *param - step;
    stat  = EG_evaluateGeomX(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = *param + step;
    stat  = EG_evaluateGeomX(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = (uplus[4] - uminus[4])/(2.0*step);
    result[1] = (uplus[5] - uminus[5])/(2.0*step);
    
  } else if (geom->oclass == CURVE) {
    
    uv[0] = *param - step;
    stat  = EG_evaluateGeomX(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = *param + step;
    stat  = EG_evaluateGeomX(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = (uplus[6] - uminus[6])/(2.0*step);
    result[1] = (uplus[7] - uminus[7])/(2.0*step);
    result[2] = (uplus[8] - uminus[8])/(2.0*step);

  } else {
    
    uv[0] = param[0] - step;
    uv[1] = param[1];
    stat  = EG_evaluateGeomX(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = param[0] + step;
    stat  = EG_evaluateGeomX(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = param[0];
    uv[1] = param[1] - step;
    stat  = EG_evaluateGeomX(geom, uv, vminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[1] = param[1] + step;
    stat  = EG_evaluateGeomX(geom, uv, vplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[ 0] = (uplus[ 9] - uminus[ 9])/(2.0*step);
    result[ 1] = (uplus[10] - uminus[10])/(2.0*step);
    result[ 2] = (uplus[11] - uminus[11])/(2.0*step);
    result[ 3] = (uplus[12] - uminus[12])/(2.0*step);
    result[ 4] = (uplus[13] - uminus[13])/(2.0*step);
    result[ 5] = (uplus[14] - uminus[14])/(2.0*step);
    result[ 6] = (vplus[12] - vminus[12])/(2.0*step);
    result[ 7] = (vplus[13] - vminus[13])/(2.0*step);
    result[ 8] = (vplus[14] - vminus[14])/(2.0*step);
    result[ 9] = (vplus[15] - vminus[15])/(2.0*step);
    result[10] = (vplus[16] - vminus[16])/(2.0*step);
    result[11] = (vplus[17] - vminus[17])/(2.0*step);
    
  }
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_evaluateGeomZ(const egObject *geom, const double *param, double *result)
{
  int          stat;
  double       data[18], norm[3], axes[15], d, dist, dum[3], radius, der3[12];
  double       *tan, *tanv, ru, *gdata;
  liteGeometry *lgeom;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;

  if (geom->oclass == PCURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3derivZ(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      norm[0] =  result[3];
      norm[1] = -result[2];
      dist    = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
      if (dist != 0.0) {
        result[0] += norm[0]*gdata[0]/dist;
        result[1] += norm[1]*gdata[0]/dist;
        axes[0]    =  result[5];
        axes[1]    = -result[4];
        radius     = (norm[0]*axes[0] + norm[1]*axes[1])/(dist*dist*dist);
        result[2] += axes[0]*gdata[0]/dist - norm[0]*radius;
        result[3] += axes[1]*gdata[0]/dist - norm[1]*radius;
        /* second derivatives -- need 3rd from pcurve! */
        axes[2]    =  der3[1];
        axes[3]    = -der3[0];
        axes[4]    = axes[2]*gdata[0]/dist;
        axes[5]    = axes[3]*gdata[0]/dist;
        axes[4]   -= 2.0*axes[0]*gdata[3]*radius;
        axes[5]   -= 2.0*axes[1]*gdata[3]*radius;
        d          = (norm[0]*axes[2] + norm[1]*axes[3] +
                      axes[0]*axes[0] + axes[1]*axes[1])/(dist*dist*dist);
        axes[4]   += norm[0]*gdata[0]*(3.0*radius*radius*dist - d);
        axes[5]   += norm[1]*gdata[0]*(3.0*radius*radius*dist - d);
        result[4] += axes[4];
        result[5] += axes[5];
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  } else if (geom->oclass == CURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3derivZ(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      tan  = &result[3];
      CROSS(norm, tan, gdata);
      dist = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
      if (dist != 0.0) {
        result[0] += norm[0]*gdata[3]/dist;
        result[1] += norm[1]*gdata[3]/dist;
        result[2] += norm[2]*gdata[3]/dist;
        tanv       = &result[6];
        CROSS(axes, tanv, gdata);
        radius     = norm[0]*axes[0] + norm[1]*axes[1] + norm[2]*axes[2];
        radius    /= dist*dist*dist;
        result[3] += axes[0]*gdata[3]/dist - norm[0]*radius;
        result[4] += axes[1]*gdata[3]/dist - norm[1]*radius;
        result[5] += axes[2]*gdata[3]/dist - norm[2]*radius;
        /* second derivatives -- need 3rd from curve! */
#ifndef __clang_analyzer__
        tanv       = &axes[3];
        CROSS(tanv, der3, gdata);
#else
        axes[3]    = axes[4] = axes[5] = 0.0;
#endif
        axes[6]    = axes[3]*gdata[3]/dist;
        axes[7]    = axes[4]*gdata[3]/dist;
        axes[8]    = axes[5]*gdata[3]/dist;
        axes[6]   -= 2.0*axes[0]*gdata[3]*radius;
        axes[7]   -= 2.0*axes[1]*gdata[3]*radius;
        axes[8]   -= 2.0*axes[2]*gdata[3]*radius;
        d          = norm[0]*axes[3] + norm[1]*axes[4] + norm[2]*axes[5] +
                     axes[0]*axes[0] + axes[1]*axes[1] + axes[2]*axes[2];
        d         /= dist*dist*dist;
        axes[6]   += norm[0]*gdata[3]*(3.0*radius*radius*dist - d);
        axes[7]   += norm[1]*gdata[3]*(3.0*radius*radius*dist - d);
        axes[8]   += norm[2]*gdata[3]*(3.0*radius*radius*dist - d);
        result[6] += axes[6];
        result[7] += axes[7];
        result[8] += axes[8];
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  } else {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == REVOLUTION) {
      stat = EG_evaluateGeomX(lgeom->ref, &param[1], data);
      if (stat != EGADS_SUCCESS) return stat;
      tan        = &gdata[3];
      tanv       = &data[3];
      CROSS(norm, tan, tanv);
      axes[6]    = norm[0];
      axes[7]    = norm[1];
      axes[8]    = norm[2];
      tanv       = &data[6];
      CROSS(norm, tan, tanv);
      axes[9]    = norm[0];
      axes[10]   = norm[1];
      axes[11]   = norm[2];
      axes[0]    = data[0] - gdata[0];
      axes[1]    = data[1] - gdata[1];
      axes[2]    = data[2] - gdata[2];
      CROSS(norm, tan, axes);
      axes[3]    = norm[0];
      axes[4]    = norm[1];
      axes[5]    = norm[2];
      norm[0]   *= sin(param[0]);
      norm[1]   *= sin(param[0]);
      norm[2]   *= sin(param[0]);
      dist       = (gdata[3]*axes[0] + gdata[4]*axes[1] +
                    gdata[5]*axes[2])*(1.0-cos(param[0]));
      norm[0]   += gdata[3]*dist;
      norm[1]   += gdata[4]*dist;
      norm[2]   += gdata[5]*dist;
      result[ 0] = axes[0]*cos(param[0]) + norm[0] + gdata[0];
      result[ 1] = axes[1]*cos(param[0]) + norm[1] + gdata[1];
      result[ 2] = axes[2]*cos(param[0]) + norm[2] + gdata[2];
      
      dist       = (gdata[3]*axes[0] + gdata[4]*axes[1] + gdata[5]*axes[2]);
      result[ 3] = axes[3]*cos(param[0]) +
                   (gdata[3]*dist  - axes[0])*sin(param[0]);
      result[ 4] = axes[4]*cos(param[0]) +
                   (gdata[4]*dist  - axes[1])*sin(param[0]);
      result[ 5] = axes[5]*cos(param[0]) +
                   (gdata[5]*dist  - axes[2])*sin(param[0]);
      
      result[ 9] = -axes[3]*sin(param[0]) +
                    (gdata[3]*dist  - axes[0])*cos(param[0]);
      result[10] = -axes[4]*sin(param[0]) +
                    (gdata[4]*dist  - axes[1])*cos(param[0]);
      result[11] = -axes[5]*sin(param[0]) +
                    (gdata[5]*dist  - axes[2])*cos(param[0]);
      
      dist       = (gdata[3]*data[3] + gdata[4]*data[4] + gdata[5]*data[5]);
      result[ 6] = data[3]*cos(param[0]) + axes[6]*sin(param[0]) +
                   gdata[3]*(1.0-cos(param[0]))*dist;
      result[ 7] = data[4]*cos(param[0]) + axes[7]*sin(param[0]) +
                   gdata[4]*(1.0-cos(param[0]))*dist;
      result[ 8] = data[5]*cos(param[0]) + axes[8]*sin(param[0]) +
                   gdata[5]*(1.0-cos(param[0]))*dist;
      
      result[12] = axes[6]*cos(param[0]) +
                   (gdata[3]*dist - data[3])*sin(param[0]);
      result[13] = axes[7]*cos(param[0]) +
                   (gdata[4]*dist - data[4])*sin(param[0]);
      result[14] = axes[8]*cos(param[0]) +
                   (gdata[5]*dist - data[5])*sin(param[0]);
      
      dist       = (gdata[3]*data[6] + gdata[4]*data[7] + gdata[5]*data[8]);
      result[15] = data[6]*cos(param[0]) + axes[ 9]*sin(param[0]) +
                   gdata[3]*(1.0-cos(param[0]))*dist;
      result[16] = data[7]*cos(param[0]) + axes[10]*sin(param[0]) +
                   gdata[4]*(1.0-cos(param[0]))*dist;
      result[17] = data[8]*cos(param[0]) + axes[11]*sin(param[0]) +
                   gdata[5]*(1.0-cos(param[0]))*dist;
      return EGADS_SUCCESS;
        
    } else if (geom->mtype == EXTRUSION) {
      stat = EG_evaluateGeomX(lgeom->ref, param, data);
      if (stat != EGADS_SUCCESS) return stat;
      result[ 0] = data[0] + param[1]*gdata[0];
      result[ 1] = data[1] + param[1]*gdata[1];
      result[ 2] = data[2] + param[1]*gdata[2];
      result[ 3] = data[3];
      result[ 4] = data[4];
      result[ 5] = data[5];
      
      result[ 6] = gdata[0];
      result[ 7] = gdata[1];
      result[ 8] = gdata[2];
      result[ 9] = data[6];
      result[10] = data[7];
      result[11] = data[8];
      result[12] = result[13] = result[14] = 0.0;
      result[15] = result[16] = result[17] = 0.0;
      return EGADS_SUCCESS;
    
    } else if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3derivZ(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      tan  = &result[3];
      tanv = &result[6];
      CROSS(norm, tan, tanv);
      dist = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
      if (dist != 0.0) {
        result[ 0] += norm[0]*gdata[0]/dist;
        result[ 1] += norm[1]*gdata[0]/dist;
        result[ 2] += norm[2]*gdata[0]/dist;

        tan         = &result[9];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[0]     = dum[0];
        axes[1]     = dum[1];
        axes[2]     = dum[2];
        tan         = &result[ 3];
        tanv        = &result[12];
        CROSS(dum, tan, tanv);
        axes[0]    += dum[0];
        axes[1]    += dum[1];
        axes[2]    += dum[2];
        tan         = &result[12];
        tanv        = &result[ 6];
        CROSS(dum, tan, tanv);
        axes[3]     = dum[0];
        axes[4]     = dum[1];
        axes[5]     = dum[2];
        tan         = &result[ 3];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[3]    += dum[0];
        axes[4]    += dum[1];
        axes[5]    += dum[2];
                  
        /* second derivatives -- need 3rd from surface! */
#ifdef __clang_analyzer__
        der3[0]     = der3[1] = der3[2] = der3[3] = der3[ 4] = der3[ 5] = 0;
        der3[6]     = der3[7] = der3[8] = der3[9] = der3[10] = der3[11] = 0;
#endif
        tan         = &der3[0];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[6]     = dum[0];
        axes[7]     = dum[1];
        axes[8]     = dum[2];
        tan         = &result[3];
        tanv        = &der3[3];
        CROSS(dum, tan, tanv);
        axes[6]    += dum[0];
        axes[7]    += dum[1];
        axes[8]    += dum[2];
        tan         = &result[9];
        tanv        = &result[12];
        CROSS(dum, tan, tanv);
        axes[6]    += 2.0*dum[0];
        axes[7]    += 2.0*dum[1];
        axes[8]    += 2.0*dum[2];
                  
        tan         = &der3[6];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[ 9]    = dum[0];
        axes[10]    = dum[1];
        axes[11]    = dum[2];
        tan         = &result[3];
        tanv        = &der3[9];
        CROSS(dum, tan, tanv);
        axes[ 9]   += dum[0];
        axes[10]   += dum[1];
        axes[11]   += dum[2];
        tan         = &result[12];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[ 9]   += 2.0*dum[0];
        axes[10]   += 2.0*dum[1];
        axes[11]   += 2.0*dum[2];
                  
        tan         = &result[12];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[12]    = dum[0];
        axes[13]    = dum[1];
        axes[14]    = dum[2];
        tan         = &result[3];
        tanv        = &der3[6];
        CROSS(dum, tan, tanv);
        axes[12]   += dum[0];
        axes[13]   += dum[1];
        axes[14]   += dum[2];
        tan         = &result[9];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[12]   += 2.0*dum[0];
        axes[13]   += 2.0*dum[1];
        axes[14]   += 2.0*dum[2];
                  
        radius      = norm[0]*axes[0] + norm[1]*axes[1] + norm[2]*axes[2];
        radius     *= gdata[0]/(dist*dist*dist);
        result[ 3] += axes[0]*gdata[0]/dist - norm[0]*radius;
        result[ 4] += axes[1]*gdata[0]/dist - norm[1]*radius;
        result[ 5] += axes[2]*gdata[0]/dist - norm[2]*radius;
        d           = norm[0]*axes[6] + norm[1]*axes[7] + norm[2]*axes[8];
        d          += axes[0]*axes[0] + axes[1]*axes[1] + axes[2]*axes[2];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*radius*dist/gdata[0] - d;
        result[ 9] += axes[6]*gdata[0]/dist - 2.0*axes[0]*radius + norm[0]*d;
        result[10] += axes[7]*gdata[0]/dist - 2.0*axes[1]*radius + norm[1]*d;
        result[11] += axes[8]*gdata[0]/dist - 2.0*axes[2]*radius + norm[2]*d;
        ru          = radius;
                  
        radius      = norm[0]*axes[3] + norm[1]*axes[4] + norm[2]*axes[5];
        radius     *= gdata[0]/(dist*dist*dist);
        result[ 6] += axes[3]*gdata[0]/dist - norm[0]*radius;
        result[ 7] += axes[4]*gdata[0]/dist - norm[1]*radius;
        result[ 8] += axes[5]*gdata[0]/dist - norm[2]*radius;
        d           = norm[0]*axes[9] + norm[1]*axes[10] + norm[2]*axes[11];
        d          += axes[3]*axes[3] + axes[4]*axes[ 4] + axes[5]*axes[ 5];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*radius*dist/gdata[0] - d;
        result[15] += axes[ 9]*gdata[0]/dist - 2.0*axes[3]*radius + norm[0]*d;
        result[16] += axes[10]*gdata[0]/dist - 2.0*axes[4]*radius + norm[1]*d;
        result[17] += axes[11]*gdata[0]/dist - 2.0*axes[5]*radius + norm[2]*d;
                  
        d           = axes[3]*axes[ 0] + axes[4]*axes[ 1] + axes[5]*axes[ 2];
        d          += norm[0]*axes[12] + norm[1]*axes[13] + norm[2]*axes[14];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*ru*dist/gdata[0] - d;
        result[12] += axes[12]*gdata[0]/dist - axes[0]*radius -
                      axes[ 3]*ru + norm[0]*d;
        result[13] += axes[13]*gdata[0]/dist - axes[1]*radius -
                      axes[ 4]*ru + norm[1]*d;
        result[14] += axes[14]*gdata[0]/dist - axes[2]*radius -
                      axes[ 5]*ru + norm[2]*d;
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  }
  
  /* should never get here! */
  return EGADS_SUCCESS;
}


/* 
 * 3rd derivatives are needed for only the Offset evaulations
 *     for now use finite differences -- later should be made analytic!
 */

__HOST_AND_DEVICE__ static int
EG_eval3deriv(const egObject *geom, const double *param, double *result)
{
  int    stat;
  double step = 1.e-7;
  double uv[2], uminus[18], uplus[18], vminus[18], vplus[18];
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) &&
      (geom->oclass != CURVE)  && (geom->oclass != SURFACE))
                                   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;

  uv[1] = 0.0;
  if (geom->oclass == PCURVE) {
    
    uv[0] = *param - step;
    stat  = EG_evaluateGeomZ(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = *param + step;
    stat  = EG_evaluateGeomZ(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = (uplus[4] - uminus[4])/(2.0*step);
    result[1] = (uplus[5] - uminus[5])/(2.0*step);
    
  } else if (geom->oclass == CURVE) {
    
    uv[0] = *param - step;
    stat  = EG_evaluateGeomZ(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = *param + step;
    stat  = EG_evaluateGeomZ(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = (uplus[6] - uminus[6])/(2.0*step);
    result[1] = (uplus[7] - uminus[7])/(2.0*step);
    result[2] = (uplus[8] - uminus[8])/(2.0*step);

  } else {
    
    uv[0] = param[0] - step;
    uv[1] = param[1];
    stat  = EG_evaluateGeomZ(geom, uv, uminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = param[0] + step;
    stat  = EG_evaluateGeomZ(geom, uv, uplus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[0] = param[0];
    uv[1] = param[1] - step;
    stat  = EG_evaluateGeomZ(geom, uv, vminus);
    if (stat != EGADS_SUCCESS) return stat;
    uv[1] = param[1] + step;
    stat  = EG_evaluateGeomZ(geom, uv, vplus);
    if (stat != EGADS_SUCCESS) return stat;
    result[ 0] = (uplus[ 9] - uminus[ 9])/(2.0*step);
    result[ 1] = (uplus[10] - uminus[10])/(2.0*step);
    result[ 2] = (uplus[11] - uminus[11])/(2.0*step);
    result[ 3] = (uplus[12] - uminus[12])/(2.0*step);
    result[ 4] = (uplus[13] - uminus[13])/(2.0*step);
    result[ 5] = (uplus[14] - uminus[14])/(2.0*step);
    result[ 6] = (vplus[12] - vminus[12])/(2.0*step);
    result[ 7] = (vplus[13] - vminus[13])/(2.0*step);
    result[ 8] = (vplus[14] - vminus[14])/(2.0*step);
    result[ 9] = (vplus[15] - vminus[15])/(2.0*step);
    result[10] = (vplus[16] - vminus[16])/(2.0*step);
    result[11] = (vplus[17] - vminus[17])/(2.0*step);
    
  }
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_evaluateGeom(const egObject *geom, const double *param, double *result)
{
  int          stat;
  double       data[18], norm[3], axes[15], d, dist, dum[3], radius, der3[12];
  double       *tan, *tanv, ru, *gdata;
  liteGeometry *lgeom;

  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;

  if (geom->oclass == PCURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3deriv(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      norm[0] =  result[3];
      norm[1] = -result[2];
      dist    = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
      if (dist != 0.0) {
        result[0] += norm[0]*gdata[0]/dist;
        result[1] += norm[1]*gdata[0]/dist;
        axes[0]    =  result[5];
        axes[1]    = -result[4];
        radius     = (norm[0]*axes[0] + norm[1]*axes[1])/(dist*dist*dist);
        result[2] += axes[0]*gdata[0]/dist - norm[0]*radius;
        result[3] += axes[1]*gdata[0]/dist - norm[1]*radius;
        /* second derivatives -- need 3rd from pcurve! */
        axes[2]    =  der3[1];
        axes[3]    = -der3[0];
        axes[4]    = axes[2]*gdata[0]/dist;
        axes[5]    = axes[3]*gdata[0]/dist;
        axes[4]   -= 2.0*axes[0]*gdata[3]*radius;
        axes[5]   -= 2.0*axes[1]*gdata[3]*radius;
        d          = (norm[0]*axes[2] + norm[1]*axes[3] +
                      axes[0]*axes[0] + axes[1]*axes[1])/(dist*dist*dist);
        axes[4]   += norm[0]*gdata[0]*(3.0*radius*radius*dist - d);
        axes[5]   += norm[1]*gdata[0]*(3.0*radius*radius*dist - d);
        result[4] += axes[4];
        result[5] += axes[5];
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  } else if (geom->oclass == CURVE) {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3deriv(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      tan  = &result[3];
      CROSS(norm, tan, gdata);
      dist = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
      if (dist != 0.0) {
        result[0] += norm[0]*gdata[3]/dist;
        result[1] += norm[1]*gdata[3]/dist;
        result[2] += norm[2]*gdata[3]/dist;
        tanv       = &result[6];
        CROSS(axes, tanv, gdata);
        radius     = norm[0]*axes[0] + norm[1]*axes[1] + norm[2]*axes[2];
        radius    /= dist*dist*dist;
        result[3] += axes[0]*gdata[3]/dist - norm[0]*radius;
        result[4] += axes[1]*gdata[3]/dist - norm[1]*radius;
        result[5] += axes[2]*gdata[3]/dist - norm[2]*radius;
        /* second derivatives -- need 3rd from curve! */
#ifndef __clang_analyzer__
        tanv       = &axes[3];
        CROSS(tanv, der3, gdata);
#else
        axes[3]    = axes[4] = axes[5] = 0.0;
#endif
        axes[6]    = axes[3]*gdata[3]/dist;
        axes[7]    = axes[4]*gdata[3]/dist;
        axes[8]    = axes[5]*gdata[3]/dist;
        axes[6]   -= 2.0*axes[0]*gdata[3]*radius;
        axes[7]   -= 2.0*axes[1]*gdata[3]*radius;
        axes[8]   -= 2.0*axes[2]*gdata[3]*radius;
        d          = norm[0]*axes[3] + norm[1]*axes[4] + norm[2]*axes[5] +
                     axes[0]*axes[0] + axes[1]*axes[1] + axes[2]*axes[2];
        d         /= dist*dist*dist;
        axes[6]   += norm[0]*gdata[3]*(3.0*radius*radius*dist - d);
        axes[7]   += norm[1]*gdata[3]*(3.0*radius*radius*dist - d);
        axes[8]   += norm[2]*gdata[3]*(3.0*radius*radius*dist - d);
        result[6] += axes[6];
        result[7] += axes[7];
        result[8] += axes[8];
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  } else {
    lgeom = (liteGeometry *) geom->blind;
    gdata = lgeom->data;
    if (gdata == NULL) return EGADS_NODATA;
    if (geom->mtype == REVOLUTION) {
      stat = EG_evaluateGeomX(lgeom->ref, &param[1], data);
      if (stat != EGADS_SUCCESS) return stat;
      tan        = &gdata[3];
      tanv       = &data[3];
      CROSS(norm, tan, tanv);
      axes[6]    = norm[0];
      axes[7]    = norm[1];
      axes[8]    = norm[2];
      tanv       = &data[6];
      CROSS(norm, tan, tanv);
      axes[9]    = norm[0];
      axes[10]   = norm[1];
      axes[11]   = norm[2];
      axes[0]    = data[0] - gdata[0];
      axes[1]    = data[1] - gdata[1];
      axes[2]    = data[2] - gdata[2];
      CROSS(norm, tan, axes);
      axes[3]    = norm[0];
      axes[4]    = norm[1];
      axes[5]    = norm[2];
      norm[0]   *= sin(param[0]);
      norm[1]   *= sin(param[0]);
      norm[2]   *= sin(param[0]);
      dist       = (gdata[3]*axes[0] + gdata[4]*axes[1] +
                    gdata[5]*axes[2])*(1.0-cos(param[0]));
      norm[0]   += gdata[3]*dist;
      norm[1]   += gdata[4]*dist;
      norm[2]   += gdata[5]*dist;
      result[ 0] = axes[0]*cos(param[0]) + norm[0] + gdata[0];
      result[ 1] = axes[1]*cos(param[0]) + norm[1] + gdata[1];
      result[ 2] = axes[2]*cos(param[0]) + norm[2] + gdata[2];
      
      dist       = (gdata[3]*axes[0] + gdata[4]*axes[1] + gdata[5]*axes[2]);
      result[ 3] = axes[3]*cos(param[0]) +
                   (gdata[3]*dist  - axes[0])*sin(param[0]);
      result[ 4] = axes[4]*cos(param[0]) +
                   (gdata[4]*dist  - axes[1])*sin(param[0]);
      result[ 5] = axes[5]*cos(param[0]) +
                   (gdata[5]*dist  - axes[2])*sin(param[0]);
      
      result[ 9] = -axes[3]*sin(param[0]) +
                    (gdata[3]*dist  - axes[0])*cos(param[0]);
      result[10] = -axes[4]*sin(param[0]) +
                    (gdata[4]*dist  - axes[1])*cos(param[0]);
      result[11] = -axes[5]*sin(param[0]) +
                    (gdata[5]*dist  - axes[2])*cos(param[0]);
      
      dist       = (gdata[3]*data[3] + gdata[4]*data[4] + gdata[5]*data[5]);
      result[ 6] = data[3]*cos(param[0]) + axes[6]*sin(param[0]) +
                   gdata[3]*(1.0-cos(param[0]))*dist;
      result[ 7] = data[4]*cos(param[0]) + axes[7]*sin(param[0]) +
                   gdata[4]*(1.0-cos(param[0]))*dist;
      result[ 8] = data[5]*cos(param[0]) + axes[8]*sin(param[0]) +
                   gdata[5]*(1.0-cos(param[0]))*dist;
      
      result[12] = axes[6]*cos(param[0]) +
                   (gdata[3]*dist - data[3])*sin(param[0]);
      result[13] = axes[7]*cos(param[0]) +
                   (gdata[4]*dist - data[4])*sin(param[0]);
      result[14] = axes[8]*cos(param[0]) +
                   (gdata[5]*dist - data[5])*sin(param[0]);
      
      dist       = (gdata[3]*data[6] + gdata[4]*data[7] + gdata[5]*data[8]);
      result[15] = data[6]*cos(param[0]) + axes[ 9]*sin(param[0]) +
                   gdata[3]*(1.0-cos(param[0]))*dist;
      result[16] = data[7]*cos(param[0]) + axes[10]*sin(param[0]) +
                   gdata[4]*(1.0-cos(param[0]))*dist;
      result[17] = data[8]*cos(param[0]) + axes[11]*sin(param[0]) +
                   gdata[5]*(1.0-cos(param[0]))*dist;
      return EGADS_SUCCESS;
        
    } else if (geom->mtype == EXTRUSION) {
      stat = EG_evaluateGeomX(lgeom->ref, param, data);
      if (stat != EGADS_SUCCESS) return stat;
      result[ 0] = data[0] + param[1]*gdata[0];
      result[ 1] = data[1] + param[1]*gdata[1];
      result[ 2] = data[2] + param[1]*gdata[2];
      result[ 3] = data[3];
      result[ 4] = data[4];
      result[ 5] = data[5];
      
      result[ 6] = gdata[0];
      result[ 7] = gdata[1];
      result[ 8] = gdata[2];
      result[ 9] = data[6];
      result[10] = data[7];
      result[11] = data[8];
      result[12] = result[13] = result[14] = 0.0;
      result[15] = result[16] = result[17] = 0.0;
      return EGADS_SUCCESS;
    
    } else if (geom->mtype == OFFSET) {
      stat = EG_evaluateGeomX(lgeom->ref, param, result);
      if (stat != EGADS_SUCCESS) return stat;
      stat = EG_eval3deriv(lgeom->ref, param, der3);
      if (stat != EGADS_SUCCESS) return stat;
      tan  = &result[3];
      tanv = &result[6];
      CROSS(norm, tan, tanv);
      dist = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
      if (dist != 0.0) {
        result[ 0] += norm[0]*gdata[0]/dist;
        result[ 1] += norm[1]*gdata[0]/dist;
        result[ 2] += norm[2]*gdata[0]/dist;

        tan         = &result[9];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[0]     = dum[0];
        axes[1]     = dum[1];
        axes[2]     = dum[2];
        tan         = &result[ 3];
        tanv        = &result[12];
        CROSS(dum, tan, tanv);
        axes[0]    += dum[0];
        axes[1]    += dum[1];
        axes[2]    += dum[2];
        tan         = &result[12];
        tanv        = &result[ 6];
        CROSS(dum, tan, tanv);
        axes[3]     = dum[0];
        axes[4]     = dum[1];
        axes[5]     = dum[2];
        tan         = &result[ 3];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[3]    += dum[0];
        axes[4]    += dum[1];
        axes[5]    += dum[2];
                  
        /* second derivatives -- need 3rd from surface! */
#ifdef __clang_analyzer__
        der3[0]     = der3[1] = der3[2] = der3[3] = der3[ 4] = der3[ 5] = 0;
        der3[6]     = der3[7] = der3[8] = der3[9] = der3[10] = der3[11] = 0;
#endif
        tan         = &der3[0];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[6]     = dum[0];
        axes[7]     = dum[1];
        axes[8]     = dum[2];
        tan         = &result[3];
        tanv        = &der3[3];
        CROSS(dum, tan, tanv);
        axes[6]    += dum[0];
        axes[7]    += dum[1];
        axes[8]    += dum[2];
        tan         = &result[9];
        tanv        = &result[12];
        CROSS(dum, tan, tanv);
        axes[6]    += 2.0*dum[0];
        axes[7]    += 2.0*dum[1];
        axes[8]    += 2.0*dum[2];
                  
        tan         = &der3[6];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[ 9]    = dum[0];
        axes[10]    = dum[1];
        axes[11]    = dum[2];
        tan         = &result[3];
        tanv        = &der3[9];
        CROSS(dum, tan, tanv);
        axes[ 9]   += dum[0];
        axes[10]   += dum[1];
        axes[11]   += dum[2];
        tan         = &result[12];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[ 9]   += 2.0*dum[0];
        axes[10]   += 2.0*dum[1];
        axes[11]   += 2.0*dum[2];
                  
        tan         = &result[12];
        tanv        = &result[6];
        CROSS(dum, tan, tanv);
        axes[12]    = dum[0];
        axes[13]    = dum[1];
        axes[14]    = dum[2];
        tan         = &result[3];
        tanv        = &der3[6];
        CROSS(dum, tan, tanv);
        axes[12]   += dum[0];
        axes[13]   += dum[1];
        axes[14]   += dum[2];
        tan         = &result[9];
        tanv        = &result[15];
        CROSS(dum, tan, tanv);
        axes[12]   += 2.0*dum[0];
        axes[13]   += 2.0*dum[1];
        axes[14]   += 2.0*dum[2];
                  
        radius      = norm[0]*axes[0] + norm[1]*axes[1] + norm[2]*axes[2];
        radius     *= gdata[0]/(dist*dist*dist);
        result[ 3] += axes[0]*gdata[0]/dist - norm[0]*radius;
        result[ 4] += axes[1]*gdata[0]/dist - norm[1]*radius;
        result[ 5] += axes[2]*gdata[0]/dist - norm[2]*radius;
        d           = norm[0]*axes[6] + norm[1]*axes[7] + norm[2]*axes[8];
        d          += axes[0]*axes[0] + axes[1]*axes[1] + axes[2]*axes[2];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*radius*dist/gdata[0] - d;
        result[ 9] += axes[6]*gdata[0]/dist - 2.0*axes[0]*radius + norm[0]*d;
        result[10] += axes[7]*gdata[0]/dist - 2.0*axes[1]*radius + norm[1]*d;
        result[11] += axes[8]*gdata[0]/dist - 2.0*axes[2]*radius + norm[2]*d;
        ru          = radius;
                  
        radius      = norm[0]*axes[3] + norm[1]*axes[4] + norm[2]*axes[5];
        radius     *= gdata[0]/(dist*dist*dist);
        result[ 6] += axes[3]*gdata[0]/dist - norm[0]*radius;
        result[ 7] += axes[4]*gdata[0]/dist - norm[1]*radius;
        result[ 8] += axes[5]*gdata[0]/dist - norm[2]*radius;
        d           = norm[0]*axes[9] + norm[1]*axes[10] + norm[2]*axes[11];
        d          += axes[3]*axes[3] + axes[4]*axes[ 4] + axes[5]*axes[ 5];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*radius*dist/gdata[0] - d;
        result[15] += axes[ 9]*gdata[0]/dist - 2.0*axes[3]*radius + norm[0]*d;
        result[16] += axes[10]*gdata[0]/dist - 2.0*axes[4]*radius + norm[1]*d;
        result[17] += axes[11]*gdata[0]/dist - 2.0*axes[5]*radius + norm[2]*d;
                  
        d           = axes[3]*axes[ 0] + axes[4]*axes[ 1] + axes[5]*axes[ 2];
        d          += norm[0]*axes[12] + norm[1]*axes[13] + norm[2]*axes[14];
        d          *= gdata[0]/(dist*dist*dist);
        d           = 3.0*radius*ru*dist/gdata[0] - d;
        result[12] += axes[12]*gdata[0]/dist - axes[0]*radius -
                      axes[ 3]*ru + norm[0]*d;
        result[13] += axes[13]*gdata[0]/dist - axes[1]*radius -
                      axes[ 4]*ru + norm[1]*d;
        result[14] += axes[14]*gdata[0]/dist - axes[2]*radius -
                      axes[ 5]*ru + norm[2]*d;
      }
      return EGADS_SUCCESS;
    } else {
      return EG_evaluateGeomX(geom, param, result);
    }
    
  }
  
  /* should never get here! */
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnPCurve(const egObject *geom, const double *coor, double *range,
                   double *t, double *uv)
{
  int    i, j, stat;
  double a, b = 0.0, dt = 0.0, tt, dist, ldist = 0.0, pw[2], result[6];
  
  /* netwon-raphson from picked position */
  for (i = 0; i < 40; i++) {
    
    /* line search */
    for (j = 0; j < 20; j++) {
      
      if (j == 19) dt = 0.0;
      tt = *t + dt;
      stat = EG_evaluateGeom(geom, &tt, result);
      if (stat != EGADS_SUCCESS) return stat;
      pw[0] = result[0] - coor[0];
      pw[1] = result[1] - coor[1];
      dist  = sqrt(pw[0]*pw[0] + pw[1]*pw[1]);
      if (dist < 1.e-08) break;
      b     = -(pw[0]*result[2] + pw[1]*result[3]);
      
      if ((i == 0) || (dist < ldist)) {
        break;
      } else {
        dt /= 2.0;
      }
    }
    
    *t = tt;
    if (j                 ==     20) break;  /* line search failed */
    if (dist               < 1.e-08) break;  /* converged! */
    if (fabs(ldist - dist) < 1.e-08) break;
    
    a  = (result[2]*result[2] + result[3]*result[3]) +
         (    pw[0]*result[4] +     pw[1]*result[5]);
    if (a == 0.0) break;
    dt = b/a;
    ldist = dist;

    /* limit dt within the valid range */
    tt = *t + dt;
    if (tt < range[0]) dt = range[0] - *t;
    if (tt > range[1]) dt = range[1] - *t;
  }
/*
  if (i == 40 || j == 20)
    printf(" EGADS Info: %d NearestOnPCurve not Converged!\n", geom->mtype);  */
  
  uv[0] = result[0];
  uv[1] = result[1];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnCurve(const egObject *geom, const double *coor, double *range,
                  double *t, double *xyz)
{
  int    i, j, stat;
  double a, b = 0.0, dt = 0.0, tt, dist, ldist = 0.0, pw[3], result[9];
  
  /* netwon-raphson from picked position */
  for (i = 0; i < 40; i++) {
    
    /* line search */
    for (j = 0; j < 20; j++) {
      
      if (j == 19) dt = 0.0;
      tt = *t + dt;
      stat = EG_evaluateGeom(geom, &tt, result);
      if (stat != EGADS_SUCCESS) return stat;
      pw[0] = result[0] - coor[0];
      pw[1] = result[1] - coor[1];
      pw[2] = result[2] - coor[2];
      dist  = sqrt(pw[0]*pw[0] + pw[1]*pw[1] + pw[2]*pw[2]);
      if (dist < 1.e-08) break;
      b     = -(pw[0]*result[3] + pw[1]*result[4] + pw[2]*result[5]);
      
      if ((i == 0) || (dist < ldist)) {
        break;
      } else {
        dt /= 2.0;
      }
    }
    
    *t = tt;
    if (j                 ==     20) break;  /* line search failed */
    if (dist               < 1.e-08) break;  /* converged! */
    if (fabs(ldist - dist) < 1.e-08) break;

    a  = (result[3]*result[3] + result[4]*result[4] + result[5]*result[5]) +
         (    pw[0]*result[6] +     pw[1]*result[7] +     pw[2]*result[8]);
    if (a == 0.0) break;
    dt = b/a;
    ldist = dist;

    /* limit dt within the valid range */
    tt = *t + dt;
    if (tt < range[0]) dt = range[0] - *t;
    if (tt > range[1]) dt = range[1] - *t;
  }
/*
  if (i == 40 || j == 20)
    printf(" EGADS Info: %d NearestOnCurve not Converged!\n", geom->mtype);  */
  
  xyz[0] = result[0];
  xyz[1] = result[1];
  xyz[2] = result[2];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnSurface(const egObject *geom, const double *point, double *uv,
                    double *coor)
{
  int    count, j, stat, per;
  double a00, a10, a11, det, dist, ldist, dx[3], uvs[2], result[18];
  double range[4], b0 = 0.0, b1 = 0.0, delta[2] = {0.0, 0.0};
  
  stat = EG_getRange(geom, range, &per);
  if (stat != EGADS_SUCCESS) return stat;

  /* newton iteration */
  ldist = uvs[0] = uvs[1] = 0.0;
  for (count = 0; count < 40; count++) {

    /* line search */
    for (j = 0; j < 20; j++) {

      if (j == 19) delta[0] = delta[1] = 0.0;
      uvs[0] = uv[0] + delta[0];
      uvs[1] = uv[1] + delta[1];

      stat = EG_evaluateGeom(geom, uvs, result);
      if (stat != EGADS_SUCCESS) return stat;
      dx[0] = result[0] - point[0];
      dx[1] = result[1] - point[1];
      dx[2] = result[2] - point[2];
      dist  = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      if (dist < 1.e-8) break;
      b0  =    -dx[0]*result[ 3] -     dx[1]*result[ 4] -     dx[2]*result[ 5];
      b1  =    -dx[0]*result[ 6] -     dx[1]*result[ 7] -     dx[2]*result[ 8];
      if ((count == 0) || (dist < ldist)) {
        break;
      } else {
        delta[0] /= 2.0;
        delta[1] /= 2.0;
      }
    }

    /* update the solution */
    uv[0] = uvs[0];
    uv[1] = uvs[1];

    if (j               ==    20) break;  /* line search failed */
    if (fabs(dist-ldist) < 1.e-8) break;
    if (dist             < 1.e-8) break;

    a00 = result[3]*result[ 3] + result[4]*result[ 4] + result[5]*result[ 5] +
              dx[0]*result[ 9] +     dx[1]*result[10] +     dx[2]*result[11];
    a10 = result[3]*result[ 6] + result[4]*result[ 7] + result[5]*result[ 8] +
              dx[0]*result[12] +     dx[1]*result[13] +     dx[2]*result[14];
    a11 = result[6]*result[ 6] + result[7]*result[ 7] + result[8]*result[ 8] +
              dx[0]*result[15] +     dx[1]*result[16] +     dx[2]*result[17];

    det = a00*a11 - a10*a10;
    if (det == 0.0) {
      coor[0] = result[0];
      coor[1] = result[1];
      coor[2] = result[2];
      return EGADS_DEGEN;
    }
    det      = 1.0/det;
    delta[0] = det*(b0*a11 - b1*a10);
    delta[1] = det*(b1*a00 - b0*a10);

    /* limit delta within the valid uv range */
    uvs[0]   = uv[0] + delta[0];
    uvs[1]   = uv[1] + delta[1];
    if (uvs[0] < range[0]) delta[0] = range[0] - uv[0];
    if (uvs[0] > range[1]) delta[0] = range[1] - uv[0];
    if (uvs[1] < range[2]) delta[1] = range[2] - uv[1];
    if (uvs[1] > range[3]) delta[1] = range[3] - uv[1];

    /* save off the last distance */
    ldist  = dist;

    /* stop if delta is very small */
    if (sqrt(delta[0]*delta[0] + delta[1]*delta[1]) < 1.e-12) break;
/*  printf("   %d: %lf %lf   %le\n", count, uv[0], uv[1], ldist);  */
  }

  coor[0] = result[0];
  coor[1] = result[1];
  coor[2] = result[2];
  if (j     == 20) return EGADS_EMPTY;
  if (count == 40) return EGADS_EMPTY;

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnPCurveLM(const egObject *geom, const double *point, double *range,
                     double *t, double *coor)
{
  int    status, istep, i;
  double lambda, ts, tt, result[6], A, b, delta;
  double obj_old, obj_tmp, data[6];
  
  /* initialize Levenberg-Marquardt */
  lambda = 1.0;
  
  ts     = *t;
  status = EG_evaluateGeom(geom, &ts, result);
  if (status != EGADS_SUCCESS) return status;
  coor[0] = result[0];
  coor[1] = result[1];
  
  obj_old = (point[0] - result[0])*(point[0] - result[0]) +
            (point[1] - result[1])*(point[1] - result[1]);
  if (obj_old < 1.e-14) return EGADS_SUCCESS;
  
  /* take Levenberg-Marquardt steps */
  for (istep = 0; istep < 250; istep++) {
    
    /* set up A  = (J' * J + lambda * diag(J' * J))
          and b  =  J' * (point - result)
       where  J' = transpose(J) */
    A = (result[2]*result[2] + result[3]*result[3]) * (1.0 + lambda);
    if (A == 0.0) break;
    b = result[2]*(point[0] - result[0]) + result[3]*(point[1] - result[1]);
    
    /* find a temp uv and associated objective function */
    delta = b/A;
    tt    = ts + delta;
    if (tt < range[0]) tt = range[0];
    if (tt > range[1]) tt = range[1];
    status = EG_evaluateGeom(geom, &tt, data);
    if (status != EGADS_SUCCESS) return status;
    
    obj_tmp = (point[0] - data[0])*(point[0] - data[0]) +
              (point[1] - data[1])*(point[1] - data[1]);
    
    /* if step is better, accept it and halve lambda (making it more Newton) */
    if (obj_tmp < obj_old) {
      obj_old = obj_tmp;
      ts      = tt;
      for (i  = 0; i < 4; i++) result[i] = data[i];
      lambda /= 2.0;
      if (lambda  < 1.e-14) lambda = 1.e-14;
      if (obj_old < 1.e-14) break;
    } else {
      /* otherwise, reject it and double lambda (more gradient-descent-like) */
      lambda *= 2.0;
    }
    if ((tt == range[0]) || (tt == range[1])) break;
    
    /* stop if delta is very small */
    if (fabs(delta) < 1.e-12) break;
  }
  
  /* return the solution */
  *t      = ts;
  coor[0] = result[0];
  coor[1] = result[1];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnCurveLM(const egObject *geom, const double *point, double *range,
                    double *t, double *coor)
{
  int    status, istep, i;
  double lambda, ts, tt, result[9], A, b, delta;
  double obj_old, obj_tmp, data[9];
  
  /* initialize Levenberg-Marquardt */
  lambda = 1.0;
  
  ts     = *t;
  status = EG_evaluateGeom(geom, &ts, result);
  if (status != EGADS_SUCCESS) return status;
  coor[0] = result[0];
  coor[1] = result[1];
  coor[2] = result[2];
  
  obj_old = (point[0] - result[0])*(point[0] - result[0]) +
            (point[1] - result[1])*(point[1] - result[1]) +
            (point[2] - result[2])*(point[2] - result[2]);
  if (obj_old < 1.e-14) return EGADS_SUCCESS;
  
  /* take Levenberg-Marquardt steps */
  for (istep = 0; istep < 250; istep++) {
    
    /* set up A  = (J' * J + lambda * diag(J' * J))
          and b  =  J' * (point - result)
       where  J' = transpose(J) */
    A = (result[3]*result[3] + result[4]*result[4] + result[5]*result[5]) *
        (1.0 + lambda);
    if (A == 0.0) break;
    b = result[3]*(point[0] - result[0]) + result[4]*(point[1] - result[1]) +
        result[5]*(point[2] - result[2]);
    
    /* find a temp uv and associated objective function */
    delta = b/A;
    tt    = ts + delta;
    if (tt < range[0]) tt = range[0];
    if (tt > range[1]) tt = range[1];
    status = EG_evaluateGeom(geom, &tt, data);
    if (status != EGADS_SUCCESS) return status;
    
    obj_tmp = (point[0] - data[0])*(point[0] - data[0]) +
              (point[1] - data[1])*(point[1] - data[1]) +
              (point[2] - data[2])*(point[2] - data[2]);
    
    /* if step is better, accept it and halve lambda (making it more Newton) */
    if (obj_tmp < obj_old) {
      obj_old = obj_tmp;
      ts      = tt;
      for (i  = 0; i < 6; i++) result[i] = data[i];
      lambda /= 2.0;
      if (lambda  < 1.e-14) lambda = 1.e-14;
      if (obj_old < 1.e-14) break;
    } else {
      /* otherwise, reject it and double lambda (more gradient-descent-like) */
      lambda *= 2.0;
    }
    if ((tt == range[0]) || (tt == range[1])) break;
    
    /* stop if delta is very small */
    if (fabs(delta) < 1.e-12) break;
  }
  
  /* return the solution */
  *t      = ts;
  coor[0] = result[0];
  coor[1] = result[1];
  coor[2] = result[2];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_nearestOnSurfaceLM(const egObject *geom, const double *point, double *uv,
                      double *coor)
{
  int    status, istep, i, per;
  double lambda, uvs[2], uvt[2], result[18], A[4], b[2], denom, delta[2], range[4];
  double obj_old, obj_tmp, data[18];
  
  /* initialize Levenberg-Marquardt */
  lambda = 1.0;
  
  uvs[0] = uv[0];
  uvs[1] = uv[1];
  status = EG_evaluateGeom(geom, uvs, result);
  if (status != EGADS_SUCCESS) return status;

  obj_old = (point[0] - result[0])*(point[0] - result[0]) +
            (point[1] - result[1])*(point[1] - result[1]) +
            (point[2] - result[2])*(point[2] - result[2]);
  if (obj_old < 1.e-14) return EGADS_SUCCESS;
  
  status = EG_getRange(geom, range, &per);
  if (status != EGADS_SUCCESS) return status;

  /* take Levenberg-Marquardt steps */
  for (istep = 0; istep < 250; istep++) {
    
    /* set up A  = (J' * J + lambda * diag(J' * J))
          and b  =  J' * (point - result)
       where  J' = transpose(J) */
    A[0] = (result[3]*result[3] + result[4]*result[4] + result[5]*result[5]) *
           (1.0 + lambda);
    A[1] =  result[3]*result[6] + result[4]*result[7] + result[5]*result[8];
    A[2] = A[1];
    A[3] = (result[6]*result[6] + result[7]*result[7] + result[8]*result[8]) *
           (1.0 + lambda);
    
    b[0] = result[3]*(point[0] - result[0]) + result[4]*(point[1] - result[1]) +
           result[5]*(point[2] - result[2]);
    b[1] = result[6]*(point[0] - result[0]) + result[7]*(point[1] - result[1]) +
           result[8]*(point[2] - result[2]);
    
    /* solve A*delta = b using Cramer's rule */
    denom    =  A[0]*A[3] - A[2]*A[1];
    if (denom == 0.0) break;
    delta[0] = (b[0]*A[3] - b[1]*A[1])/denom;
    delta[1] = (A[0]*b[1] - A[2]*b[0])/denom;
    
    /* find a temp uv and associated objective function */
    uvt[0]   = uvs[0] + delta[0];
    uvt[1]   = uvs[1] + delta[1];
    if (uvt[0] < range[0]) { uvt[0] = range[0]; delta[0] = uvt[0] - uvs[0]; }
    if (uvt[0] > range[1]) { uvt[0] = range[1]; delta[0] = uvt[0] - uvs[0]; }
    if (uvt[1] < range[2]) { uvt[1] = range[2]; delta[1] = uvt[1] - uvs[1]; }
    if (uvt[1] > range[3]) { uvt[1] = range[3]; delta[1] = uvt[1] - uvs[1]; }
    status   = EG_evaluateGeom(geom, uvt, data);
    if (status != EGADS_SUCCESS) return status;
    
    obj_tmp = (point[0] - data[0])*(point[0] - data[0]) +
              (point[1] - data[1])*(point[1] - data[1]) +
              (point[2] - data[2])*(point[2] - data[2]);
    
    /* if step is better, accept it and halve lambda (making it more Newton) */
    if (obj_tmp < obj_old) {
      obj_old = obj_tmp;
      uvs[0]  = uvt[0];
      uvs[1]  = uvt[1];
      for (i  = 0; i < 9; i++) result[i] = data[i];
      lambda /= 2.0;
      if (lambda  < 1.e-14) lambda = 1.e-14;
      if (obj_old < 1.e-14) break;
    } else {
      /* otherwise, reject it and double lambda (more gradient-descent-like) */
      lambda *= 2.0;
    }
    
    /* stop if delta is very small */
    if (sqrt(delta[0]*delta[0] + delta[1]*delta[1]) < 1.e-12) break;
  }
  
  /* return the solution */
  uv[0]   = uvs[0];
  uv[1]   = uvs[1];
  coor[0] = result[0];
  coor[1] = result[1];
  coor[2] = result[2];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static void
EG_orderCandidates(liteIndex *cand, double dist2, double *uv, int i, int j)
{
  if (dist2 >= cand[3].dist2) return;
  
  if (dist2 < cand[0].dist2) {
    cand[3]       = cand[2];
    cand[2]       = cand[1];
    cand[1]       = cand[0];
    cand[0].uk    = i;
    cand[0].vk    = j;
    cand[0].uv[0] = uv[0];
    if (j != 0) cand[0].uv[1] = uv[1];
    cand[0].dist2 = dist2;
  } else if (dist2 < cand[1].dist2) {
    cand[3]       = cand[2];
    cand[2]       = cand[1];
    cand[1].uk    = i;
    cand[1].vk    = j;
    cand[1].uv[0] = uv[0];
    if (j != 0) cand[1].uv[1] = uv[1];
    cand[1].dist2 = dist2;
  } else if (dist2 < cand[2].dist2) {
    cand[3]       = cand[2];
    cand[2].uk    = i;
    cand[2].vk    = j;
    cand[2].uv[0] = uv[0];
    if (j != 0) cand[2].uv[1] = uv[1];
    cand[2].dist2 = dist2;
  } else {
    cand[3].uk    = i;
    cand[3].vk    = j;
    cand[3].uv[0] = uv[0];
    if (j != 0) cand[3].uv[1] = uv[1];
    cand[3].dist2 = dist2;
  }
}


__HOST_AND_DEVICE__ int
EG_invEvaGeomLimits(const egObject *geomx, /*@null@*/ const double *limits,
                    const double *xyz, double *param, double toler,
                    double *result)
{
  int            i, j, ii, iii, jjj, k, stat, per, atype, ulen, vlen, cnt;
  int            jDiv, iDiv;
  double         a, b, tx, tt, period, tol, coord[3], srange[4] = {0.,0.,0.,0.};
  double         pt[3], uvs[2], uvx[2], range[4], data[18];
  double         *urats = NULL, *vrats = NULL;
  liteGeometry   *lgeom, *lref;
  liteIndex      cand[4];
  const int      *ints;
  const double   *reals;
  const char     *str;
  const egObject *geom;
  static double  ratios[5]  = {0.02, 0.25, 0.5,  0.75, 0.98};
  static double  finrat[10] = {0.02, 0.1,  0.2,  0.3,  0.4,
                               0.5,  0.6,  0.7,  0.8,  0.98};
  static double  xfinrt[20] = {0.02, 0.05, 0.1,  0.15, 0.2,
                               0.25, 0.3,  0.35, 0.4,  0.45,
                               0.5,  0.55, 0.6,  0.65, 0.7,
                               0.75, 0.8,  0.85, 0.9,  0.98};
  static double  percrv[11] = {0./8., 0.010, 1./8., 2./8., 3./8., 4./8.,
                               5./8., 6./8., 7./8., 0.990, 8./8.};

  geom = geomx;
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE) &&
      (geom->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  
  stat = EG_getRange(geom, range, &per);
  if (stat != EGADS_SUCCESS) return stat;
  srange[0] = range[0];
  srange[1] = range[1];
  if (geom->oclass == SURFACE) {
    srange[2] = range[2];
    srange[3] = range[3];
  }
  
  /* do we re-limit? */
  if (limits != NULL) {
    range[0] = limits[0];
    range[1] = limits[1];
    if (geom->oclass == SURFACE) {
      range[2] = limits[2];
      range[3] = limits[3];
    }
  }
  
  if (geom->oclass == PCURVE) {
    
    while (geom->mtype == TRIMMED) {
      lgeom = (liteGeometry *) geom->blind;
      geom  = lgeom->ref;
    }
    
    /* find good starting point */
    b = 0.0;
    if (geom->mtype == BEZIER) {
      for (i = 0; i < 20; i++) {
        tx   = (1.0-xfinrt[i])*range[0] + xfinrt[i]*range[1];
        stat = EG_evaluateGeom(geom, &tx, data);
        if (stat != EGADS_SUCCESS) return stat;
        a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
            (data[1]-xyz[1])*(data[1]-xyz[1]);
        if (i == 0) {
          *param = tx;
          b      = a;
        } else {
          if (a < b) {
            *param = tx;
            b      = a;
          }
        }
      }
    } else if (geom->mtype == BSPLINE) {
      lgeom         = (liteGeometry *) geom->blind;
      cand[0].dist2 = cand[1].dist2 = cand[2].dist2 = cand[3].dist2 = 1.e308;
      cand[0].uk    = cand[1].uk    = cand[2].uk    = cand[3].uk    = 0;
      cand[0].vk    = cand[1].vk    = cand[2].vk    = cand[3].vk    = 0;
      cand[0].uv[0] = cand[1].uv[0] = cand[2].uv[0] = cand[3].uv[0] = 0;
      cand[0].uv[1] = cand[1].uv[1] = cand[2].uv[1] = cand[3].uv[1] = 0;
      k    = lgeom->header[1];
      stat = EG_attributeRet(geom, ".Bad", &atype, &ulen, &ints, &reals, &str);
      if ((stat == EGADS_SUCCESS) && (atype == ATTRSTRING))
        if (EG_strncmp(str, "fold", 4) == 0) k *= 2;
      for (j = i = 1; i < lgeom->header[3]; i++) {
        if (lgeom->data[i-1] <  range[0])       continue;
        if (lgeom->data[i-1] == lgeom->data[i]) continue;
        tx = range[1];
        for (ii = 1; ii <= k; ii++) {
          tx = lgeom->data[i-1] + ii*(lgeom->data[i] - lgeom->data[i-1])/(k+1);
          if (tx > range[1]) break;
          stat = EG_evaluateGeom(geom, &tx, data);
          if (stat != EGADS_SUCCESS) return stat;
          a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
              (data[1]-xyz[1])*(data[1]-xyz[1]);
          EG_orderCandidates(cand, a, &tx, i, 0);
          if (j == 1) {
            *param = tx;
            b      = a;
          } else {
            if (a < b) {
              *param = tx;
              b      = a;
            }
          }
          j++;
        }
        if (tx > range[1]) break;
      }
      if (j < 20) {
        for (i = 0; i < 20; i++) {
          tx   = (1.0-xfinrt[i])*range[0] + xfinrt[i]*range[1];
          stat = EG_evaluateGeom(geom, &tx, data);
          if (stat != EGADS_SUCCESS) return stat;
          a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
              (data[1]-xyz[1])*(data[1]-xyz[1]);
          EG_orderCandidates(cand, a, &tx, i+1, 0);
          if (j == 1) {
            *param = tx;
            b      = a;
          } else {
            if (a < b) {
              *param = tx;
              b      = a;
            }
          }
          j++;
        }
      }
    } else {
      if (geom->mtype == BEZIER) {
        urats = xfinrt;
        ulen = 20;
      } else if (geom->mtype == LINE) {
        urats = ratios;
        ulen = 1;
      } else if ((geom->mtype == CIRCLE) || (geom->mtype == ELLIPSE)) {
        urats = percrv;
        ulen = 11;
      } else {
        urats = finrat;
        ulen = 10;
      }
      for (i = 0; i < ulen; i++) {
        tx   = (1.0-urats[i])*range[0] + urats[i]*range[1];
        stat = EG_evaluateGeom(geom, &tx, data);
        if (stat != EGADS_SUCCESS) return stat;
        a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
            (data[1]-xyz[1])*(data[1]-xyz[1]);
        if (i == 0) {
          *param = tx;
          b      = a;
        } else {
          if (a < b) {
            *param = tx;
            b      = a;
          }
        }
      }
    }
    tx   = *param;
    stat = EG_evaluateGeom(geom, &tx, data);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = data[0];
    result[1] = data[1];
    a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
        (data[1]-xyz[1])*(data[1]-xyz[1]);
    if (geom->mtype == BSPLINE) {
      for (i = 0; i < 2; i++) {
        if (cand[i].uk == 0) continue;
        tt   = cand[i].uv[0];
        stat = EG_nearestOnPCurve(geom, xyz, range, &tt, coord);
        if (stat != EGADS_SUCCESS) return stat;
        b = (coord[0]-xyz[0])*(coord[0]-xyz[0]) +
            (coord[1]-xyz[1])*(coord[1]-xyz[1]);
        if (b < a) {
          a         = b;
          tx        = *param;
          data[0]   = result[0];
          data[1]   = result[1];
          *param    = tt;
          result[0] = coord[0];
          result[1] = coord[1];
        }
      }
      EG_nearestOnPCurveLM(geom, xyz, range, param, result);
    } else {
      stat = EG_nearestOnPCurve(geom, xyz, range, param, result);
      if (stat != EGADS_SUCCESS) return stat;
    }
    b = (result[0]-xyz[0])*(result[0]-xyz[0]) +
        (result[1]-xyz[1])*(result[1]-xyz[1]);
    if (b > a) {
      printf(" EGADS Info: %d NearestOnP diverge %le vs %le\n",
             geom->mtype, sqrt(a), sqrt(b));
      *param    = tx;
      result[0] = data[0];
      result[1] = data[1];
    }
    
  } else if (geom->oclass == CURVE) {

    while (geom->mtype == TRIMMED) {
      lgeom = (liteGeometry *) geom->blind;
      geom  = lgeom->ref;
    }
    
    /* find good starting point */
    b = 0.0;
    if (geom->mtype == BSPLINE) {
      lgeom         = (liteGeometry *) geom->blind;
      cand[0].dist2 = cand[1].dist2 = cand[2].dist2 = cand[3].dist2 = 1.e308;
      cand[0].uk    = cand[1].uk    = cand[2].uk    = cand[3].uk    = 0;
      cand[0].vk    = cand[1].vk    = cand[2].vk    = cand[3].vk    = 0;
      cand[0].uv[0] = cand[1].uv[0] = cand[2].uv[0] = cand[3].uv[0] = 0;
      cand[0].uv[1] = cand[1].uv[1] = cand[2].uv[1] = cand[3].uv[1] = 0;
      k = lgeom->header[1];
      for (j = i = 1; i < lgeom->header[3]; i++) {
        if (lgeom->data[i-1] <  range[0])       continue;
        if (lgeom->data[i-1] == lgeom->data[i]) continue;
        tx = range[1];
        for (ii = 1; ii <= k; ii++) {
          tx = lgeom->data[i-1] + ii*(lgeom->data[i] - lgeom->data[i-1])/(k+1);
          if (tx > range[1]) break;
          stat = EG_evaluateGeom(geom, &tx, data);
          if (stat != EGADS_SUCCESS) return stat;
          a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
              (data[1]-xyz[1])*(data[1]-xyz[1]) +
              (data[2]-xyz[2])*(data[2]-xyz[2]);
          EG_orderCandidates(cand, a, &tx, i, 0);
          if (j == 1) {
            *param = tx;
            b      = a;
          } else {
            if (a < b) {
              *param = tx;
              b      = a;
            }
          }
          j++;
        }
        if (tx > range[1]) break;
      }
      if (j < 20) {
        for (i = 0; i < 20; i++) {
          tx   = (1.0-xfinrt[i])*range[0] + xfinrt[i]*range[1];
          stat = EG_evaluateGeom(geom, &tx, data);
          if (stat != EGADS_SUCCESS) return stat;
          a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
              (data[1]-xyz[1])*(data[1]-xyz[1]) +
              (data[2]-xyz[2])*(data[2]-xyz[2]);
          EG_orderCandidates(cand, a, &tx, i+1, 0);
          if (j == 1) {
            *param = tx;
            b      = a;
          } else {
            if (a < b) {
              *param = tx;
              b      = a;
            }
          }
          j++;
        }
      }
    } else {
      if (geom->mtype == BEZIER) {
        urats = xfinrt;
        ulen = 20;
      } else if (geom->mtype == LINE) {
        urats = ratios;
        ulen = 1;
      } else if ((geom->mtype == CIRCLE) || (geom->mtype == ELLIPSE)) {
        urats = percrv;
        ulen = 11;
      } else {
        urats = finrat;
        ulen = 10;
      }
      for (i = 0; i < ulen; i++) {
        tx   = (1.0-urats[i])*range[0] + urats[i]*range[1];
        stat = EG_evaluateGeom(geom, &tx, data);
        if (stat != EGADS_SUCCESS) return stat;
        a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
            (data[1]-xyz[1])*(data[1]-xyz[1]) +
            (data[2]-xyz[2])*(data[2]-xyz[2]);
        if (i == 0) {
          *param = tx;
          b      = a;
        } else {
          if (a < b) {
            *param = tx;
            b      = a;
          }
        }
      }
    }
    tx   = *param;
    stat = EG_evaluateGeom(geom, &tx, data);
    if (stat != EGADS_SUCCESS) return stat;
    result[0] = data[0];
    result[1] = data[1];
    result[2] = data[2];
    a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
        (data[1]-xyz[1])*(data[1]-xyz[1]) +
        (data[2]-xyz[2])*(data[2]-xyz[2]);
    if (geom->mtype == BSPLINE) {
      for (i = 0; i < 2; i++) {
        if (cand[i].uk == 0) continue;
        tt   = cand[i].uv[0];
        stat = EG_nearestOnCurve(geom, xyz, range, &tt, coord);
        if (stat != EGADS_SUCCESS) return stat;
        b = (coord[0]-xyz[0])*(coord[0]-xyz[0]) +
            (coord[1]-xyz[1])*(coord[1]-xyz[1]) +
            (coord[2]-xyz[2])*(coord[2]-xyz[2]);
        if (b < a) {
          a         = b;
          tx        = *param;
          data[0]   = result[0];
          data[1]   = result[1];
          data[2]   = result[2];
          *param    = tt;
          result[0] = coord[0];
          result[1] = coord[1];
          result[2] = coord[2];
        }
      }
      EG_nearestOnCurveLM(geom, xyz, range, param, result);
    } else {
      stat = EG_nearestOnCurve(geom, xyz, range, param, result);
      if (stat != EGADS_SUCCESS) return stat;
    }
    b = (result[0]-xyz[0])*(result[0]-xyz[0]) +
        (result[1]-xyz[1])*(result[1]-xyz[1]) +
        (result[2]-xyz[2])*(result[2]-xyz[2]);
    if (b > a) {
      printf(" EGADS Info: %d NearestOnC diverge %le vs %le\n",
             geom->mtype, sqrt(a), sqrt(b));
      *param    = tx;
      result[0] = data[0];
      result[1] = data[1];
      result[2] = data[2];
    }

  } else {

    tol = toler;
    if (tol == 0.0) tol = 1.e-8;

    while (geom->mtype == TRIMMED) {
      lgeom = (liteGeometry *) geom->blind;
      geom  = lgeom->ref;
    }

    /* do different things based on surface type */
    b = 1.e308;
    if (geom->mtype == BSPLINE) {
      lgeom         = (liteGeometry *) geom->blind;
      cand[0].dist2 = cand[1].dist2 = cand[2].dist2 = cand[3].dist2 = 1.e308;
      cand[0].uk    = cand[1].uk    = cand[2].uk    = cand[3].uk    = 0;
      cand[0].vk    = cand[1].vk    = cand[2].vk    = cand[3].vk    = 0;
      cand[0].uv[0] = cand[1].uv[0] = cand[2].uv[0] = cand[3].uv[0] = 0;
      cand[0].uv[1] = cand[1].uv[1] = cand[2].uv[1] = cand[3].uv[1] = 0;
      cnt  = 0;
      iDiv = lgeom->header[3]-2*lgeom->header[1];
      jDiv = lgeom->header[6]-2*lgeom->header[4];
      if ((iDiv > 2) && (jDiv > 2)) {
        for (j = lgeom->header[4]+1;
             j < lgeom->header[6]-lgeom->header[4]; j++) {
          if (lgeom->data[j+lgeom->header[3]  ] ==
              lgeom->data[j+lgeom->header[3]-1]) continue;
          uvs[1] = 0.5*(lgeom->data[j+lgeom->header[3]  ] +
                        lgeom->data[j+lgeom->header[3]-1]);
          if (uvs[1] < range[2]) continue;
          if (uvs[1] > range[3]) break;
          for (i = lgeom->header[1]+1;
               i < lgeom->header[3]-lgeom->header[1]; i++) {
            if (lgeom->data[i] == lgeom->data[i-1]) continue;
            uvs[0] = 0.5*(lgeom->data[i] + lgeom->data[i-1]);
            if (uvs[0] < range[0]) continue;
            if (uvs[0] > range[1]) break;
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            EG_orderCandidates(cand, a, uvs, i, j);
            cnt++;
          }
        }
      } else if ((iDiv > 2) && (jDiv <= 2)) {
        for (j = 1; j < 4; j++) {
          uvs[1] = range[2] + j*(range[3]-range[2])/4.0;
          for (i = lgeom->header[1]+1;
               i < lgeom->header[3]-lgeom->header[1]; i++) {
            if (lgeom->data[i] == lgeom->data[i-1]) continue;
            uvs[0] = 0.5*(lgeom->data[i] + lgeom->data[i-1]);
            if (uvs[0] < range[0]) continue;
            if (uvs[0] > range[1]) break;
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            EG_orderCandidates(cand, a, uvs, i, j);
            cnt++;
          }
        }
      } else if ((iDiv <= 2) && (jDiv > 2)) {
        for (j = lgeom->header[4]+1;
             j < lgeom->header[6]-lgeom->header[4]; j++) {
          if (lgeom->data[j+lgeom->header[3]  ] ==
              lgeom->data[j+lgeom->header[3]-1]) continue;
          uvs[1] = 0.5*(lgeom->data[j+lgeom->header[3]  ] +
                        lgeom->data[j+lgeom->header[3]-1]);
          if (uvs[1] < range[2]) continue;
          if (uvs[1] > range[3]) break;
          for (i = 1; i < 4; i++) {
            uvs[0] = range[0] + i*(range[1]-range[0])/4.0;
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            EG_orderCandidates(cand, a, uvs, i, j);
            cnt++;
          }
        }
      }
      /* was the sampling enough? */
      if ( (cnt == 0) ||
          ((cnt < 2) && ((lgeom->header[1] <= 2) || (lgeom->header[4] <= 2)))) {
        /* no -- just do a fine UV grid */
        for (j = 0; j < 10; j++) {
          uvs[1] = (1.0-finrat[j])*range[2] + finrat[j]*range[3];
          for (i = 0; i < 10; i++) {
            uvs[0] = (1.0-finrat[i])*range[0] + finrat[i]*range[1];
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            if (a < b) {
              b         = a;
              param[0]  = uvs[0];
              param[1]  = uvs[1];
              result[0] = data[0];
              result[1] = data[1];
              result[2] = data[2];
            }
          }
        }
        uvs[0] = param[0];
        uvs[1] = param[1];
        stat   = EG_nearestOnSurface(geom, xyz, param, pt);
        if (stat == EGADS_SUCCESS) {
          a    = (pt[0]-xyz[0])*(pt[0]-xyz[0]) +
                 (pt[1]-xyz[1])*(pt[1]-xyz[1]) +
                 (pt[2]-xyz[2])*(pt[2]-xyz[2]);
          if (b < a) {
            param[0]  = uvs[0];
            param[1]  = uvs[1];
          } else {
            result[0] = pt[0];
            result[1] = pt[1];
            result[2] = pt[2];
          }
        } else {
          param[0]  = uvs[0];
          param[1]  = uvs[1];
        }
      } else {
        /* yes -- subsample based on order */
        for (k = 0; k < 4; k++) {
          i = cand[k].uk;
          j = cand[k].vk;
          if (iDiv <= 2) i = lgeom->header[1]+1;
          if (jDiv <= 2) j = lgeom->header[4]+1;
          if ((i == 0) || (j == 0)) continue;
          for (jjj = 1; jjj <= lgeom->header[4]; jjj++) {
            uvs[1] =      lgeom->data[j+lgeom->header[3]-1] +
                     jjj*(lgeom->data[j+lgeom->header[3]  ] -
                          lgeom->data[j+lgeom->header[3]-1])/(lgeom->header[4]+1);
            for (iii = 1; iii <= lgeom->header[1]; iii++) {
              uvs[0] =      lgeom->data[i-1] +
                       iii*(lgeom->data[i]-lgeom->data[i-1])/(lgeom->header[1]+1);
              stat   = EG_evaluateGeom(geom, uvs, data);
              if (stat != EGADS_SUCCESS) continue;
              a = tx = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                       (data[1]-xyz[1])*(data[1]-xyz[1]) +
                       (data[2]-xyz[2])*(data[2]-xyz[2]);
              uvx[0] = uvs[0];
              uvx[1] = uvs[1];
              stat   = EG_nearestOnSurface(geom, xyz, uvx, pt);
              if (stat == EGADS_SUCCESS)
                if ((uvx[0] >= srange[0]) && (uvx[0] <= srange[1]) &&
                    (uvx[1] >= srange[2]) && (uvx[1] <= srange[3]))
                  tx = (pt[0]-xyz[0])*(pt[0]-xyz[0]) +
                       (pt[1]-xyz[1])*(pt[1]-xyz[1]) +
                       (pt[2]-xyz[2])*(pt[2]-xyz[2]);
              if (tx < a) {
                if (tx < b) {
                  b         = tx;
                  param[0]  = uvx[0];
                  param[1]  = uvx[1];
                  result[0] = pt[0];
                  result[1] = pt[1];
                  result[2] = pt[2];
                }
              } else {
                if (a < b) {
                  b         = a;
                  param[0]  = uvs[0];
                  param[1]  = uvs[1];
                  result[0] = data[0];
                  result[1] = data[1];
                  result[2] = data[2];
                }
              }
              if (b < tol*tol) break;
            }
            if (b < tol*tol) break;
          }
          if (b < tol*tol) break;
        }
      }
      /* sometimes the second derivative puts us in bad places! */
      EG_nearestOnSurfaceLM(geom, xyz, param, result);
    } else if (geom->mtype == EXTRUSION) {
      lgeom = (liteGeometry *) geom->blind;
      if (lgeom->ref->mtype == BSPLINE) {
        lref = (liteGeometry *) lgeom->ref->blind;
        ii   = lref->header[1];
        ii  *= lref->header[3] - 2*lref->header[1];
        if (ii < 6) ii = 6;
        for (j = 0; j < 10; j++) {
          uvs[1] = (1.0-finrat[j])*range[2] + finrat[j]*range[3];
          for (i = 0; i < ii; i++) {
            uvs[0] = range[0] + i*(range[1]-range[0])/(ii-1);
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            if (a < b) {
              b        = a;
              param[0] = uvs[0];
              param[1] = uvs[1];
            }
          }
        }
        uvs[0] = param[0];
        uvs[1] = param[1];
        stat   = EG_nearestOnSurface(geom, xyz, param, result);
        a      = (result[0]-xyz[0])*(result[0]-xyz[0]) +
                 (result[1]-xyz[1])*(result[1]-xyz[1]) +
                 (result[2]-xyz[2])*(result[2]-xyz[2]);
        if ((b < a) || (stat != EGADS_SUCCESS)) {
/*        printf(" Info: NEAREST Fails -- Seed point closer    stat = %d\n",
                 stat);  */
          EG_evaluateGeom(geom, uvs, data);
          param[0]  = uvs[0];
          param[1]  = uvs[1];
          result[0] = data[0];
          result[1] = data[1];
          result[2] = data[2];
        }
        /* sometimes the second derivative puts us in bad places! */
        EG_nearestOnSurfaceLM(geom, xyz, param, result);
      } else {
        for (j = 0; j < 20; j++) {
          uvs[1] = (1.0-xfinrt[j])*range[2] + xfinrt[j]*range[3];
          for (i = 0; i < 20; i++) {
            uvs[0] = (1.0-xfinrt[i])*range[0] + xfinrt[i]*range[1];
            stat   = EG_evaluateGeom(geom, uvs, data);
            if (stat != EGADS_SUCCESS) continue;
            a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
                (data[1]-xyz[1])*(data[1]-xyz[1]) +
                (data[2]-xyz[2])*(data[2]-xyz[2]);
            if (a < b) {
              b        = a;
              param[0] = uvs[0];
              param[1] = uvs[1];
            }
          }
        }
        uvs[0] = param[0];
        uvs[1] = param[1];
        stat   = EG_nearestOnSurface(geom, xyz, param, result);
        a      = (result[0]-xyz[0])*(result[0]-xyz[0]) +
                 (result[1]-xyz[1])*(result[1]-xyz[1]) +
                 (result[2]-xyz[2])*(result[2]-xyz[2]);
        if (b < a) {
          printf(" EGADS Info: NearestOn EXTRUSION fails -- Seed closer stat = %d\n",
                 stat);
          EG_evaluateGeom(geom, uvs, data);
          param[0]  = uvs[0];
          param[1]  = uvs[1];
          result[0] = data[0];
          result[1] = data[1];
          result[2] = data[2];
        }
      }
    } else {

      if ((geom->mtype == BEZIER) ||
          (geom->mtype == OFFSET)) {
        urats = xfinrt;
        ulen = 20;
        vrats = xfinrt;
        vlen = 20;
      } else if (geom->mtype == REVOLUTION) {
        urats = percrv;
        ulen = 11;
        vrats = xfinrt;
        vlen = 20;
      } else if (geom->mtype == PLANE) {
        urats = ratios;
        ulen = 1;
        vrats = ratios;
        vlen = 1;
      } else if ((geom->mtype == SPHERICAL) ||
                 (geom->mtype == TOROIDAL)) {
        urats = percrv;
        ulen = 11;
        vrats = percrv;
        vlen = 11;
      } else if ((geom->mtype == CYLINDRICAL) ||
                 (geom->mtype == CONICAL)) {
        urats = percrv;
        ulen = 11;
        vrats = finrat;
        vlen = 10;
      } else {
        urats = finrat;
        ulen = 10;
        vrats = finrat;
        vlen = 10;
      }

      for (j = 0; j < vlen; j++) {
        uvs[1] = (1.0-vrats[j])*range[2] + vrats[j]*range[3];
        for (i = 0; i < ulen; i++) {
          uvs[0] = (1.0-urats[i])*range[0] + urats[i]*range[1];
          stat   = EG_evaluateGeom(geom, uvs, data);
          if (stat != EGADS_SUCCESS) continue;
          a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
              (data[1]-xyz[1])*(data[1]-xyz[1]) +
              (data[2]-xyz[2])*(data[2]-xyz[2]);
          if (a < b) {
            b        = a;
            param[0] = uvs[0];
            param[1] = uvs[1];
          }
        }
      }
      uvs[0] = param[0];
      uvs[1] = param[1];
      stat   = EG_nearestOnSurface(geom, xyz, param, result);
      a      = (result[0]-xyz[0])*(result[0]-xyz[0]) +
               (result[1]-xyz[1])*(result[1]-xyz[1]) +
               (result[2]-xyz[2])*(result[2]-xyz[2]);
      if (b < a) {
        printf(" EGADS Info: NearestOn %d %dx%d fails -- Seed closer stat = %d\n",
               geom->mtype, ulen, vlen, stat);
        EG_evaluateGeom(geom, uvs, data);
        param[0]  = uvs[0];
        param[1]  = uvs[1];
        result[0] = data[0];
        result[1] = data[1];
        result[2] = data[2];
      }
    }

    if ((per&1) != 0) {
      period = srange[1] - srange[0];
      if ((param[0]+PARAMACC < srange[0]) || (param[0]-PARAMACC > srange[1]))
        if (param[0]+PARAMACC < srange[0]) {
          if (param[0]+period-PARAMACC < srange[1]) param[0] += period;
        } else {
          if (param[0]-period+PARAMACC > srange[0]) param[0] -= period;
        }
    }
    if ((per&2) != 0) {
      period = srange[3] - srange[2];
      if ((param[1]+PARAMACC < srange[2]) || (param[1]-PARAMACC > srange[3]))
        if (param[1]+PARAMACC < srange[2]) {
          if (param[1]+period-PARAMACC < srange[3]) param[1] += period;
        } else {
          if (param[1]-period+PARAMACC > srange[2]) param[1] -= period;
        }
    }
    
  }
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_invEvaluateGeomGuess(const egObject *geom, /*@null@*/ const double *limits,
                        double *xyz, double *param, double *result)
{
  int    stat, per;
  double range[2];
  
  if  (geom == NULL)               return EGADS_NULLOBJ;
  if  (geom->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if ((geom->oclass != PCURVE) && (geom->oclass != CURVE)  &&
      (geom->oclass != SURFACE))   return EGADS_NOTGEOM;
  if  (geom->blind == NULL)        return EGADS_NODATA;
  
  if (geom->oclass == PCURVE) {

    stat = EG_getRange(geom, range, &per);
    if (stat != EGADS_SUCCESS) return stat;
    stat = EG_nearestOnPCurve(geom, xyz, range, param, result);
    
  } else if (geom->oclass == CURVE) {
    
    if (limits == NULL) {
      stat = EG_getRange(geom, range, &per);
      if (stat != EGADS_SUCCESS) return stat;
    } else {
      range[0] = limits[0];
      range[1] = limits[1];
    }
    stat = EG_nearestOnCurve(geom, xyz, range, param, result);
    
  } else {
    
    /* Surfaces */
    stat = EG_nearestOnSurface(geom, xyz, param, result);
    if ((stat == EGADS_EMPTY) || (stat == EGADS_DEGEN)) stat = EGADS_SUCCESS;
    
  }

  return stat;
}


__HOST_AND_DEVICE__ static int
EG_inFace3D(const egObject *face, const double *uv, /*@null@*/ double *p,
            /*@null@*/ double *uvx)
{
  int      i, j, stat;
  double   a, b, tx, param[2], pt[3], xyz[3], tang[3], edata[9], data[18];
  double   tol, uvtol, ttol, nrm[3], norm[3], tng[3];
  egObject *loop, *edge, *surface, *node;
  liteEdge *pedge;
  liteLoop *ploop;
  liteFace *pface;
  
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  pface   = (liteFace *) face->blind;
  surface = pface->surface;
  if (surface == NULL)            return EGADS_NULLOBJ;

  uvtol = pface->urange[1] - pface->urange[0];
  if (pface->vrange[1] - pface->vrange[0] > uvtol)
    uvtol = pface->vrange[1] - pface->vrange[0];
  uvtol *= 1.e-10;
  if (p == NULL) {
    if (uv[0]+uvtol < pface->urange[0]) return EGADS_OUTSIDE;
    if (uv[0]-uvtol > pface->urange[1]) return EGADS_OUTSIDE;
    if (uv[1]+uvtol < pface->vrange[0]) return EGADS_OUTSIDE;
    if (uv[1]-uvtol > pface->vrange[1]) return EGADS_OUTSIDE;
  }
  
  /* get closest edge point to our UV and store away the tangent */
  node     = NULL;
  b        = 1.e308;
  pt[0]    = pt[1]    = pt[2]   = 0.0;
  tang[0]  = tang[1]  = tang[2] = 0.0;
  param[0] = param[1] = 0.0;
  stat     = EG_evaluate(surface, uv, data);
  if (stat != EGADS_SUCCESS) return stat;
  for (i = 0; i < pface->nloops; i++) {
    loop  = pface->loops[i];
    if (loop == NULL) continue;
    ploop = (liteLoop *) loop->blind;
    if (ploop == NULL) continue;
    for (j = 0; j < ploop->nedges; j++) {
      edge = ploop->edges[j];
      if (edge == NULL) continue;
      if (edge->mtype == DEGENERATE) continue;
      pedge = (liteEdge *) edge->blind;
      tol   = pedge->tol;
      stat  = EG_invEvaluate(edge, data, &tx, xyz);
      if (stat != EGADS_SUCCESS) return stat;
      a = (data[0]-xyz[0])*(data[0]-xyz[0]) +
          (data[1]-xyz[1])*(data[1]-xyz[1]) +
          (data[2]-xyz[2])*(data[2]-xyz[2]);
      if (a < b) {
        b     = a;
        node  = NULL;
        pt[0] = xyz[0];
        pt[1] = xyz[1];
        pt[2] = xyz[2];
        stat  = EG_getEdgeUV(face, edge, ploop->senses[j], tx, param);
        if (stat != EGADS_SUCCESS) return stat;
        if (b < tol*tol) {
          if (p != NULL) {
            p[0]   = pt[0];
            p[1]   = pt[1];
            p[2]   = pt[2];
          }
          if (uvx != NULL) {
            uvx[0] = param[0];
            uvx[1] = param[1];
          }
          return EGADS_SUCCESS;
        }
        stat = EG_evaluate(edge, &tx, edata);
        if (stat != EGADS_SUCCESS) return stat;
        tang[0] = ploop->senses[j]*edata[3];
        tang[1] = ploop->senses[j]*edata[4];
        tang[2] = ploop->senses[j]*edata[5];
        ttol    = (pedge->trange[1] - pedge->trange[0])*1.e-9;
        if (tx <= pedge->trange[0]+ttol) node = pedge->nodes[0];
        if (tx >= pedge->trange[1]-ttol) node = pedge->nodes[1];
      }
    }
  }
  if (b > 1.e307) return EGADS_NOTFOUND;
  
  /* correct tangent (average of both endpoints) for hitting a Node */
  if (node != NULL) {
    tang[0] = tang[1] = tang[2] = 0.0;
    for (i = 0; i < pface->nloops; i++) {
      loop  = pface->loops[i];
      if (loop == NULL) continue;
      ploop = (liteLoop *) loop->blind;
      if (ploop == NULL) continue;
      for (j = 0; j < ploop->nedges; j++) {
        edge = ploop->edges[j];
        if (edge == NULL) continue;
        if (edge->mtype == DEGENERATE) continue;
        pedge = (liteEdge *) edge->blind;
        if (pedge->nodes[0] == node) {
          stat = EG_evaluate(edge, &pedge->trange[0], edata);
          if (stat != EGADS_SUCCESS) return stat;
          tng[0] = ploop->senses[j]*edata[3];
          tng[1] = ploop->senses[j]*edata[4];
          tng[2] = ploop->senses[j]*edata[5];
          b      = sqrt(tng[0]*tng[0] + tng[1]*tng[1] + tng[2]*tng[2]);
          if (b != 0.0) {
            tang[0] += tng[0]/b;
            tang[1] += tng[1]/b;
            tang[2] += tng[2]/b;
          }
        }
        if (pedge->nodes[1] == node) {
          stat = EG_evaluate(edge, &pedge->trange[1], edata);
          if (stat != EGADS_SUCCESS) return stat;
          tng[0] = ploop->senses[j]*edata[3];
          tng[1] = ploop->senses[j]*edata[4];
          tng[2] = ploop->senses[j]*edata[5];
          b      = sqrt(tng[0]*tng[0] + tng[1]*tng[1] + tng[2]*tng[2]);
          if (b != 0.0) {
            tang[0] += tng[0]/b;
            tang[1] += tng[1]/b;
            tang[2] += tng[2]/b;
          }
        }
      }
    }
  }
  
  /* got nearest point to an edge -- look at what side we are on */
  if (p != NULL) {
    p[0]   = pt[0];
    p[1]   = pt[1];
    p[2]   = pt[2];
  }
  if (uvx != NULL) {
    uvx[0] = param[0];
    uvx[1] = param[1];
  }
  pt[0] -= data[0];
  pt[1] -= data[1];
  pt[2] -= data[2];
  a      = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);
  if (a < pface->tol) return EGADS_SUCCESS;
  
  pt[0]   /= a;
  pt[1]   /= a;
  pt[2]   /= a;
  b        = sqrt(tang[0]*tang[0] + tang[1]*tang[1] + tang[2]*tang[2]);
  tang[0] /= b;
  tang[1] /= b;
  tang[2] /= b;
  CROSS(nrm, pt, tang);
  stat     = EG_evaluate(surface, param, data);
  if (stat != EGADS_SUCCESS) return stat;
  pt[0]    = data[3];
  pt[1]    = data[4];
  pt[2]    = data[5];
  tang[0]  = data[6];
  tang[1]  = data[7];
  tang[2]  = data[8];
  CROSS(norm, pt, tang);
  if (face->mtype == SREVERSE) {
    norm[0] = -norm[0];
    norm[1] = -norm[1];
    norm[2] = -norm[2];
  }
  b = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
  if (b == 0.0) return EGADS_DEGEN;
  norm[0] /= b;
  norm[1] /= b;
  norm[2] /= b;

  if (norm[0]*nrm[0]+norm[1]*nrm[1]+norm[2]*nrm[2] < 0.0) return EGADS_OUTSIDE;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_inFaceX(const egObject *face, const double *uva, /*@null@*/ double *pt,
           /*@null@*/ double *uvx)
{
  int      i, j, stat, sense, sper;
  double   dist, d, dd, t, ts, utol, vtol, ttol, etol, period;
  double   uv[2], srange[4], uvs[2], xyz[3], data[18], eval[9];
  egObject *loop, *edge, *surface, *pcurve, *pcurv, *node;
  liteNode *pnode;
  liteEdge *pedge;
  liteLoop *ploop;
  liteFace *pface;
  
  if (face == NULL)               return EGADS_NULLOBJ;
  if (face->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (face->oclass != FACE)       return EGADS_NOTTOPO;
  if (face->blind == NULL)        return EGADS_NODATA;
  pface   = (liteFace *) face->blind;
  surface = pface->surface;
  if (surface == NULL)            return EGADS_NULLOBJ;
  if (surface->mtype == PLANE)    return EG_inFace3D(face, uva, pt, uvx);
  
  /* only works when there are PCurves */
  utol  = (pface->urange[1] - pface->urange[0])*1.e-9;
  vtol  = (pface->vrange[1] - pface->vrange[0])*1.e-9;
  uv[0] = uva[0];
  uv[1] = uva[1];
  stat  = EG_getRange(surface, srange, &sper);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_getRange PC = %d in EG_inFace!\n", stat);
    return stat;
  }
  /* is uv in range? */
  if ((sper&1) != 0) {
    period = srange[1] - srange[0];
    if ((uv[0]+utol < pface->urange[0]) || (uv[0]-utol > pface->urange[1]))
      if (uv[0]+utol < pface->urange[0]) {
        while (uv[0]+period-utol < pface->urange[1]) uv[0] += period;
      } else {
        while (uv[0]-period+utol > pface->urange[0]) uv[0] -= period;
      }
  }
  if ((sper&2) != 0) {
    period = srange[3] - srange[2];
    if ((uv[1]+vtol < pface->vrange[0]) || (uv[1]-vtol > pface->vrange[1]))
      if (uv[1]+vtol < pface->vrange[0]) {
        while (uv[1]+period-vtol < pface->vrange[1]) uv[1] += period;
      } else {
        while (uv[1]-period+vtol > pface->vrange[0]) uv[1] -= period;
      }
  }
  if (pt == NULL) {
    if (uv[0]+utol < pface->urange[0]) return EGADS_OUTSIDE;
    if (uv[0]-utol > pface->urange[1]) return EGADS_OUTSIDE;
    if (uv[1]+vtol < pface->vrange[0]) return EGADS_OUTSIDE;
    if (uv[1]-vtol > pface->vrange[1]) return EGADS_OUTSIDE;
  }
  stat = EG_evaluate(surface, uv, data);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_evaluate S = %d in EG_inFace!\n", stat);
    return stat;
  }
  xyz[0] = data[0];
  xyz[1] = data[1];
  xyz[2] = data[2];
  
  /* get closest edge point to our UV and store away the pcurve/t */
  dist  = 1.e300;
  pcurv = node = NULL;
  ts    = 0.0;
  sense = 1;
  for (i = 0; i < pface->nloops; i++) {
    loop  = pface->loops[i];
    if (loop == NULL) continue;
    ploop = (liteLoop *) loop->blind;
    if (ploop == NULL) continue;
    for (j = 0; j < ploop->nedges; j++) {
      edge = ploop->edges[j];
      if (edge == NULL) continue;
      if (edge->blind == NULL) continue;
      if (edge->mtype == DEGENERATE) continue;
      pedge  = (liteEdge *) edge->blind;
      etol   = pedge->tol;
      pcurve = ploop->edges[j+ploop->nedges];
      if (pcurve == NULL) continue;
      stat   = EG_invEvaGeomLimits(pcurve, pedge->trange, uv, &ts, 0.0, uvs);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_invEvaluate = %d in EG_inFace!\n", stat);
        return stat;
      }
      /* are we as close as the relative UV tolerance? */
      if ((fabs(uv[0]-uvs[0]) < utol) && (fabs(uv[1]-uvs[1]) < vtol)) {
        if (pt != NULL) {
          stat = EG_evaluate(edge, &ts, eval);
          if (stat != EGADS_SUCCESS) {
            printf(" EG_evaluate e = %d in EG_inFace!\n", stat);
            return stat;
          }
          pt[0]  = eval[0];
          pt[1]  = eval[1];
          pt[2]  = eval[2];
        }
        if (uvx != NULL) {
          uvx[0] = uvs[0];
          uvx[1] = uvs[1];
        }
        return EGADS_SUCCESS;
      }
      d = (uv[0]-uvs[0])*(uv[0]-uvs[0]) + (uv[1]-uvs[1])*(uv[1]-uvs[1]);
      if (d >= dist) continue;
      /* are we within the Edge tolerance? */
      stat = EG_evaluate(edge, &ts, eval);
      if (stat != EGADS_SUCCESS) {
        printf(" EG_evaluate E = %d in EG_inFace!\n", stat);
        return stat;
      }
      dd = sqrt((xyz[0]-eval[0])*(xyz[0]-eval[0]) +
                (xyz[1]-eval[1])*(xyz[1]-eval[1]) +
                (xyz[2]-eval[2])*(xyz[2]-eval[2]));
      if (dd <= etol) {
        if (pt != NULL) {
          pt[0]  = eval[0];
          pt[1]  = eval[1];
          pt[2]  = eval[2];
        }
        if (uvx != NULL) {
          uvx[0] = uvs[0];
          uvx[1] = uvs[1];
        }
        return EGADS_SUCCESS;
      }
      if (d < dist) {
        dist  = d;
        pcurv = pcurve;
        t     = ts;
        sense = ploop->senses[j];
        node  = NULL;
        /* hit Node (can't find internal Edge segment) */
        ttol  = (pedge->trange[1] - pedge->trange[0])*1.e-9;
        if ((ts <= pedge->trange[0]+ttol) || (ts >= pedge->trange[1]-ttol))
          if (pedge->nodes[0] != pedge->nodes[1]) {
            if (ts <= pedge->trange[0]+ttol) node = pedge->nodes[0];
            if (ts >= pedge->trange[1]-ttol) node = pedge->nodes[1];
          }
      }
    }
  }
  if (pcurv == NULL) {
    printf(" No Edge Found in EG_inFace!\n");
    return EGADS_NOTFOUND;
  }

  /* do we need to pivot around a Node (find min Edge dot)? */
  if (node != NULL) {
    pnode = (liteNode *) node->blind;
    if (pnode == NULL) {
      printf(" NULL Node in EG_inFace!\n");
      return EGADS_NODATA;
    }
/*  printf(" Nearest on Node @ (%lf,%lf,%lf)\n", pnode->xyz[0], pnode->xyz[1],
           pnode->xyz[2]);  */
    dist = 1.0;
    for (i = 0; i < pface->nloops; i++) {
      loop  = pface->loops[i];
      if (loop == NULL) continue;
      ploop = (liteLoop *) loop->blind;
      if (ploop == NULL) continue;
      for (j = 0; j < ploop->nedges; j++) {
        edge = ploop->edges[j];
        if (edge == NULL) continue;
        if (edge->blind == NULL) continue;
        pedge = (liteEdge *) edge->blind;
        if (pedge->nodes[0] == node) {
          pcurve = ploop->edges[j+ploop->nedges];
          if (pcurve != NULL) {
            ts   = pedge->trange[0];
            stat = EG_evaluate(pcurve, &ts, data);
            if (stat != EGADS_SUCCESS) {
              printf("  EGADS Internal: EG_evaluate = %d in EG_inFace!\n", stat);
              return stat;
            }
            /* dot product between tangent & direction to point */
            uvs[0] = uv[0] - data[0];
            uvs[1] = uv[1] - data[1];
            d = sqrt(uvs[0]*uvs[0] + uvs[1]*uvs[1]);
            if (d != 0.0) {
              uvs[0] /= d;
              uvs[1] /= d;
            }
            d = sqrt(data[2]*data[2] + data[3]*data[3]);
            if (d != 0.0) {
              data[2] /= d;
              data[3] /= d;
            }
            d = uvs[0]*data[2] + uvs[1]*data[3];
/*          printf("   hit first node for %d/%d  %lf\n", j+1, i+1, d);  */
            if (fabs(d) < dist) {
              dist  = fabs(d);
              pcurv = pcurve;
              t     = ts;
              sense = ploop->senses[j];
            }
          }
        }
        if (pedge->nodes[1] == node) {
          pcurve = ploop->edges[j+ploop->nedges];
          if (pcurve != NULL) {
            ts   = pedge->trange[1];
            stat = EG_evaluate(pcurve, &ts, data);
            if (stat != EGADS_SUCCESS) {
              printf("  EGADS Internal: EG_evaluate = %d in EG_inFace!\n", stat);
              return stat;
            }
            /* dot product between tangent & direction to point */
            uvs[0] = uv[0] - data[0];
            uvs[1] = uv[1] - data[1];
            d = sqrt(uvs[0]*uvs[0] + uvs[1]*uvs[1]);
            if (d != 0.0) {
              uvs[0] /= d;
              uvs[1] /= d;
            }
            d = sqrt(data[2]*data[2] + data[3]*data[3]);
            if (d != 0.0) {
              data[2] /= d;
              data[3] /= d;
            }
            d = uvs[0]*data[2] + uvs[1]*data[3];
/*          printf("   hit last  node for %d/%d  %lf\n", j+1, i+1, d);  */
            if (fabs(d) < dist) {
              dist  = fabs(d);
              pcurv = pcurve;
              t     = ts;
              sense = ploop->senses[j];
            }
          }
        }
      }
    }
  }
  
  /* what side of the Edge are we on? */
  stat = EG_evaluate(pcurv, &t, data);
  if (stat != EGADS_SUCCESS) {
    printf("  EGADS Internal: EG_evaluate = %d in EG_inFace!\n", stat);
    return stat;
  }
  uvs[0] = uv[0] - data[0];
  uvs[1] = uv[1] - data[1];
  d = sqrt(uvs[0]*uvs[0] + uvs[1]*uvs[1]);
  if (d != 0.0) {
    uvs[0] /= d;
    uvs[1] /= d;
  }
  d = sqrt(data[2]*data[2] + data[3]*data[3]);
  if (d != 0.0) {
    data[2] /= d;
    data[3] /= d;
  }
  d = uvs[0]*data[3] - uvs[1]*data[2];
  if (d == 0.0) {
    printf(" EGADS Internal: Cross = 0 (EG_inFace)!\n");
    printf("                 %le %le  %le %le\n",
           uvs[0], uvs[1], data[2], data[3]);
  }

  if (d*sense*face->mtype > 0.0) {
/*  printf("  d = %le  sense = %d  face or = %d\n", d, sense, face->mtype);
    printf("          uv = %lf %lf  on Edge (%lf) = %lf %lf  %d\n",
           uv[0], uv[1], t, data[0], data[1], pcurv->mtype);
    printf("          uvbox = %lf %lf  %lf %lf\n", pface->urange[0],
           pface->urange[1], pface->vrange[0], pface->vrange[1]);
    printf("                  %lf %lf  %lf %lf  %d\n",
           srange[0], srange[1], srange[2], srange[3], sper);
    for (i = 0; i < pface->nloops; i++) {
      loop  = pface->loops[i];
      if (loop == NULL) continue;
      ploop = (liteLoop *) loop->blind;
      if (ploop == NULL) continue;
      for (j = 0; j < ploop->nedges; j++) {
        edge = ploop->edges[j];
        if (edge == NULL) continue;
        if (edge->blind == NULL) continue;
        pedge  = (liteEdge *) edge->blind;
        pcurve = ploop->edges[j+ploop->nedges];
        if (pcurve != pcurv) continue;
        for (stat = 0; stat < 100; stat++) {
          ts = pedge->trange[0] + stat*(pedge->trange[1] - pedge->trange[0])/99.;
          EG_evaluate(pcurv, &ts, eval);
          printf("       %lf: %lf %lf  %lf %lf\n",
                 ts, eval[0], eval[1], eval[2], eval[3]);
        }
      }
    }  */
    if (uvx != NULL) {
      uvx[0] = data[0];
      uvx[1] = data[1];
    }
    if (pt != NULL) {
      if (node != NULL) {
        pnode = (liteNode *) node->blind;
        if (pnode != NULL) {
          pt[0] = pnode->xyz[0];
          pt[1] = pnode->xyz[1];
          pt[2] = pnode->xyz[2];
          return EGADS_OUTSIDE;
        }
      }
      uvs[0] = data[0];
      uvs[1] = data[1];
      stat   = EG_evaluate(surface, uvs, data);
      if (stat != EGADS_SUCCESS)
        printf(" EGADS Internal: eval = %d!\n", stat);
      pt[0]  = data[0];
      pt[1]  = data[1];
      pt[2]  = data[2];
    }
    return EGADS_OUTSIDE;
  }
  
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_getPerpDir(double *norm, egObject *edge, int sense, double t, double *dir)
{
  int    stat;
  double dist, result[9], *dt;
  
  dt   = &result[3];
  stat = EG_evaluate(edge, &t, result);
  if (stat != EGADS_SUCCESS) return stat;
  CROSS(dir, norm, dt);
  dist = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
  if (sense == SFORWARD) dist = -dist;
  if (dist != 0.0) {
    dir[0] /= dist;
    dir[1] /= dist;
    dir[2] /= dist;
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ int
EG_getWindingAngle(egObject *edge, double t, double *angle)
{
  int      i, j, k, n, stat, cx, nface, oclass, mtype, nloop, sense, *senses;
  double   dist, ang, uv[2], lims[4], result[18], *du, *dv;
  double   norms[2][3], dir[2][3], n1[3], n2[3], x0[2], x1[2], p0[2], p1[2];
  egObject *body, *ref, *sur, **faces, **loops, **edges;
  
  *angle = -90.0;
  if (edge == NULL)               return EGADS_NULLOBJ;
  if (edge->magicnumber != MAGIC) return EGADS_NOTOBJ;
  if (edge->oclass != EDGE)       return EGADS_NOTTOPO;
  if (edge->blind == NULL)        return EGADS_NODATA;
  if (edge->mtype == DEGENERATE)  return EGADS_DEGEN;
  
  stat = EG_getBody(edge, &body);
  if (stat != EGADS_SUCCESS)      return stat;
  if (body == NULL)               return EGADS_NODATA;
  if (body->mtype != SOLIDBODY)   return EGADS_NOTBODY;
  
  stat = EG_getBodyTopos(body, edge, FACE, &nface, &faces);
  if (stat != EGADS_SUCCESS) return stat;
  if (nface == 1) {
    *angle = 180.0;
    EG_free(faces);
    return EGADS_SUCCESS;
  } else if (nface != 2) {
    EG_free(faces);
    return EGADS_TOPOCNT;
  }
  
  du = &result[3];
  dv = &result[6];
  for (i = 0; i < nface; i++) {
    /* find the sense */
    stat = EG_getTopology(faces[i], &sur, &oclass, &mtype, lims, &nloop, &loops,
                          &senses);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      return stat;
    }
    for (sense = j = 0; j < nloop; j++) {
      stat = EG_getTopology(loops[j], &ref, &oclass, &mtype, lims, &n, &edges,
                            &senses);
      if (stat != EGADS_SUCCESS) {
        EG_free(faces);
        return stat;
      }
      for (k = 0; k < n; k++) {
        if (edges[k] != edge) continue;
        sense = senses[k];
        break;
      }
      if (sense != 0) break;
    }
    /* get the normals */
    stat = EG_getEdgeUV(faces[i], edge, sense, t, uv);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      return stat;
    }
    stat = EG_evaluate(faces[i], uv, result);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      return stat;
    }
    CROSS(norms[i], du, dv);
    dist = sqrt(norms[i][0]*norms[i][0] + norms[i][1]*norms[i][1] +
                norms[i][2]*norms[i][2]);
    if (faces[i]->mtype == SREVERSE) dist = -dist;
    if (dist != 0.0) {
      norms[i][0] /= dist;
      norms[i][1] /= dist;
      norms[i][2] /= dist;
    }
    /* get the perpendicular to the Edge */
    stat = EG_getPerpDir(norms[i], edge, sense, t, dir[i]);
    if (stat != EGADS_SUCCESS) {
      EG_free(faces);
      return stat;
    }
  }
  EG_free(faces);
  
  /* project on to 2D by cross products of normals */
  CROSS(n2, norms[0], norms[1]);
  dist = sqrt(DOT(n2, n2));
  if (dist < 1.e-7) {
                                   *angle =   0.0;
    if (DOT(dir[0], dir[1]) < 0.0) *angle = 180.0;
    return EGADS_SUCCESS;
  }
  n2[0] /= dist;
  n2[1] /= dist;
  n2[2] /= dist;
  CROSS(n1, norms[0], n2);
  
  x0[0] = DOT(norms[0], norms[0]);
  x0[1] = DOT(n1,       norms[0]);
  x1[0] = DOT(norms[0], norms[1]);
  x1[1] = DOT(n1,       norms[1]);
  p0[0] = DOT(norms[0], dir[0]);
  p0[1] = DOT(n1,       dir[0]);
  p1[0] = DOT(norms[0], dir[1]);
  p1[1] = DOT(n1,       dir[1]);
  
  /* are we convex? */
                                     cx = 0;
  if (x0[0]*p1[0]+x0[1]*p1[1] > 0.0) cx = 1;
  if (x1[0]*p0[0]+x1[1]*p0[1] > 0.0) cx = 1;
  dist = DOT(dir[0], dir[1]);
  if ((dist < -1.0) || (dist > 1.0))
    printf(" EG_getWindingAngle: dot = %20.12lf\n", dist);
  if (dist < -1.0) dist = -1.0;
  if (dist >  1.0) dist =  1.0;
  ang  = acos(dist);
  if (cx == 1) ang = 2.0*PI - ang;

  *angle = 180.0*ang/PI;
  return EGADS_SUCCESS;
}
