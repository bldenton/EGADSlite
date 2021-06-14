#ifndef BSPLINEFIT_H
#define BSPLINEFIT_H

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif

__ProtoExt__ int EG_spline1dEval( int *ivec, double *rdata, double t,
                                  double *point );
__ProtoExt__ int EG_spline1dEval_dot( int *ivec, const double *rdata,
                                      const double *rdata_dot, double t,
                                      double *point, double *point_dot );
__ProtoExt__ int EG_spline1dDeriv( int *ivec, double *rdata, int der, double t,
                                   double *deriv );
__ProtoExt__ int EG_spline1dDeriv_dot( int *ivec, double *rdata,
                                       double *rdata_dot, int der, double t,
                                       double *deriv, double *deriv_dot );

__ProtoExt__ int EG_spline1dFit( int endx, int imaxx, const double *xyz,
                                 const double *kn, double tol, int *header,
                                 double **rdata );
__ProtoExt__ int EG_spline1dFit_dot( int endx, int imaxx,
                                     const double *xyz, const double *xyz_dot,
                                     const double *kn, const double *kn_dot,
                                     double tol, int *header,
                                     double **rdata, double **rdata_dot );

__ProtoExt__ int EG_spline1dTan( int imaxx, const double *t1, const double *xyz,
                                 const double *tn, const double *kn, double tol,
                                 int *ivec, double **rdata );
__ProtoExt__ int EG_spline1dTan_dot( int imaxx,
                                     const double *t1,  const double *t1_dot,
                                     const double *xyz, const double *xyz_dot,
                                     const double *tn,  const double *tn_dot,
                                     const double *kn,  const double *kn_dot,
                                     double tol, int *ivec,
                                     double **rdata,    double **rdata_dot );

__ProtoExt__ int EG_spline2dEval( int *ivec, double *data, const double *uv,
                                  double *point );
__ProtoExt__ int EG_spline2dDeriv( int *ivec, double *data, int der,
                                   const double *uv, double *deriv );
  
__ProtoExt__ int EG_spline2dAprx( int e, int im, int jm, const double *xyz,
                                  const double *uknot,   const double *vknot,
                                  const int    *vdata,
                                  const double *wesT,    const double *easT,
                                  const double *south,         double *snor,
                                  const double *north,         double *nnor,
                                  double tol, int *header,     double **rdata );

#ifdef __cplusplus
} /* extern "C" */

#include "Surreal/SurrealS.h"

/* Surreal interface for computing derivatives of splines directly */

template<int N, class T>
int EG_spline1dEval(int *ivec, SurrealS<N> *data, T& t, SurrealS<N> *point);

template<int N, class T>
int EG_spline1dDeriv(int *ivec, SurrealS<N> *data, int der, T& t,
                     SurrealS<N> *point);

template<class T>
int EG_spline1dFit(int endx, int imaxx, const T *xyz,
                   const T *kn, double tol, int *ivec,
                   T **rdata);

template<class T>
int EG_spline1dTan(int imaxx, const T *t1, const T *xyz, const T *tn,
                   const T *kn, double tol, int *ivec,
                   T **rdata);

template<int N, class T>
int EG_spline2dEval(int *ivec, SurrealS<N> *data, const T *uv,
                    SurrealS<N> *point);

template<int N, class T>
int EG_spline2dDeriv(int *ivec, SurrealS<N> *data, int der, const T *uv,
                     SurrealS<N> *deriv);

template<int N>
int EG_spline2dAprx(int endc, int imax, int jmax, const SurrealS<N> *xyz,
                    const SurrealS<N> *uknot,     const SurrealS<N> *vknot,
                    const int *vdata,
                    const SurrealS<N> *wesT,      const SurrealS<N> *easT,
                    const SurrealS<N> *south,           SurrealS<N> *snor,
                    const SurrealS<N> *north,           SurrealS<N> *nnor,
                    double tol, int *header,            SurrealS<N> **rdata);

#ifdef WIN32
/* all explicit instantiations need to be exposed here for WIN32 */
template __declspec( dllimport )
         int EG_spline1dEval<1, double>(int *, SurrealS<1> *, double&,
                                        SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline1dEval<1, SurrealS<1>>(int *, SurrealS<1> *, SurrealS<1>&,
                                             SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline1dFit< double >(int, int, const double *,
                                      const double *,
                                      double, int *, double **);
template __declspec( dllimport )
         int EG_spline1dFit< SurrealS<1> >(int, int, const SurrealS<1> *,
                                           const SurrealS<1> *,
                                           double, int *, SurrealS<1> **);

template __declspec( dllimport )
         int EG_spline1dTan< SurrealS<1> >(int imaxx, const SurrealS<1> *t1,
                                           const SurrealS<1> *xyz,
                                           const SurrealS<1> *tn,
                                           const SurrealS<1> *kn, double tol,
                                           int *ivec, SurrealS<1> **rdata);

template __declspec( dllimport )
         int EG_spline2dEval<1, double>(int *, SurrealS<1> *, const double *,
                                        SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline2dEval<1, SurrealS<1>>(int *, SurrealS<1> *,
                                             const SurrealS<1> *, SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline2dDeriv<1, double>(int *, SurrealS<1> *, int,
                                         const double *, SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline2dDeriv<1, SurrealS<1>>(int *, SurrealS<1> *, int,
                                              const SurrealS<1> *, SurrealS<1> *);

template __declspec( dllimport )
         int EG_spline2dAprx<1>(int endc, int imax, int jmax,
                                const SurrealS<1> *, const SurrealS<1> *,
                                const SurrealS<1> *, const int *,
                                const SurrealS<1> *, const SurrealS<1> *,
                                const SurrealS<1> *,       SurrealS<1> *,
                                const SurrealS<1> *,       SurrealS<1> *,
                                double tol, int *header,   SurrealS<1> **);
#endif

#endif

#endif
