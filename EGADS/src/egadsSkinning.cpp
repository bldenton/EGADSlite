/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Skin a series of BSpline Curves
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include "egads.h"
#include "egads_dot.h"

#include <stdio.h>
#include <math.h>

extern "C" int EG_sameThread( const ego object );
extern "C" int EG_outLevel( const ego object );

#define EPS 1.E-08
#define PI  3.1415926535897931159979635


/* Computes Bsplines values for a given knot sequence */
template<class T>
static T
OneBasisFun (int p, int m, T U[], int i, T u)
{
  int j, k;
  T   saved, Uleft, Uright, temp, *N;

  // Special cases
  if ( ((i == 0) && (u == U[0])) ||  ((i == m-p-1) && (u == U[m])) ) return 1.0;
  if ( (u < U[i]) || (u >= U[i+p+1]) ) return 0.0;
  N = (T *) EG_alloc((p+1) *sizeof(T));
  if (N == NULL) return EGADS_MALLOC;

  for (j = 0; j <= p; ++j) {
    if ( ( u >= U[i+j] ) && ( u < U[i+j+1] ) ) {
      N[j] = 1.0;
    } else {
      N[j] = 0.0;
    }
  }
  for ( k = 1; k <= p; ++k ) {
    if ( N[0] == 0.0 ) {
      saved = 0.0;
    } else {
      saved = ( (u-U[i]) * N[0] )/( U[i+k]-U[i] );
    }
    for ( j = 0; j < p-k+1; ++j ) {
      Uleft  = U[i+j+1];
      Uright = U[i+j+k+1];
      if (N[j+1] == 0.0) {
        N[j]  = saved;
        saved = 0.0;
      } else {
        temp  = N[j+1]/( Uright - Uleft );
        N[j]  = saved + ( Uright - u )*temp;
        saved = ( u - Uleft )*temp;
      }
    }
  }
  saved = N[0];
  EG_free(N);
  return saved;
}


template<class T>
static int
matsol(T    A[],           /* (in)  matrix to be solved (stored rowwise) */
                           /* (out) upper-triangular form of matrix */
       T    b[],           /* (in)  right hand side */
                           /* (out) right-hand side after swapping */
       int  n,             /* (in)  size of matrix */
       T    x[])           /* (out) solution of A*x=b */
{

  int ir, jc, kc, imax;
  T   amax, swap, fact;
  double EPS12=1.E-12;

  /* --------------------------------------------------------------- */

  /* reduce each column of A */
  for (kc = 0; kc < n; kc++) {

    /* find pivot element */
    imax = kc;
    amax = fabs(A[kc*n+kc]);

    for (ir = kc+1; ir < n; ir++) {
      if (fabs(A[ir*n+kc]) > amax) {
        imax = ir;
        amax = fabs(A[ir*n+kc]);
      }
    }

    /* check for possibly-singular matrix (ie, near-zero pivot) */
    if (amax < EPS12) return EGADS_DEGEN;

    /* if diagonal is not pivot, swap rows in A and b */
    if (imax != kc) {
      for (jc = 0; jc < n; jc++) {
        swap         = A[kc  *n+jc];
        A[kc  *n+jc] = A[imax*n+jc];
        A[imax*n+jc] = swap;
      }

      swap    = b[kc  ];
      b[kc  ] = b[imax];
      b[imax] = swap;
    }

    /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
    for (ir = kc+1; ir < n; ir++) {
      fact = A[ir*n+kc] / A[kc*n+kc];

      for (jc = kc+1; jc < n; jc++) A[ir*n+jc] -= fact * A[kc*n+jc];

      b[ir] -= fact * b[kc];

      A[ir*n+kc] = 0;
    }
  }

  /* back-substitution pass */
  x[n-1] = b[n-1] / A[(n-1)*n+(n-1)];
  for (jc = n-2; jc >= 0; jc--) {
    x[jc] = b[jc];
    for (kc = jc+1; kc < n; kc++) x[jc] -= A[jc*n+kc] * x[kc];
    x[jc] /= A[jc*n+jc];
  }
  return EGADS_SUCCESS;
}


// This function is a duplicate. Copied from egasSplineFit.cpp
template<class T>
static int
FindSpan(int n, int p, const T& u, T *U)
{
  int low, mid, high;

  if ( fabs ( u - U[n+1] ) < EPS) return n;
  low  = p; high = n+1;
  mid  = ( low + high )/2;
  while ( ( u < U[mid] ) || ( u >= U[mid+1] ) ) {
    if ( u < U[mid] ) {
      high = mid;
    } else {
      low  = mid;
    }
    mid = ( low + high )/2;
  }
  return mid;
}


/* Finds the knot multiplicity */
template<class T>
static int
getMultiplicity(T *vector, int idx, int size)
{
  int j, m;

  m = 1;
  for ( j = idx+1; j < size; ++j )
    if ( fabs ( vector[idx] - vector[j] ) < EPS ) ++m;

  return m;
}


/* Merges two knot vectors (with highest multiplicity) and stores the resulting
   vector in kOut (kOut gets overwritten!);
    - kOut should have at least dim = dim_kIn+dim_kOut
    - *dim_kOut saves the ACTUAL length of the merged vector
*/
template<class T>
static int
mergeKnotVectors(int dim_kIn, T *kIn, int *dim_kOut, T *kOut)
{
  int i, j, m, idx, mu, mv, maxm, dimAux;
  T *auxK = NULL;

  if ( *dim_kOut == 0 ) {
    if ( dim_kIn == 0 )  return EGADS_EMPTY;
    *dim_kOut = dim_kIn;
    for (i = 0 ; i< dim_kIn; ++i) kOut[i] = kIn[i];
    return EGADS_SUCCESS;
  }

  dimAux = *dim_kOut;
  auxK   = (T *) EG_alloc(dimAux*sizeof(T));
  if (auxK == NULL ) return EGADS_MALLOC;
  // Copy data from kOut because it will overwrite it
  for ( i = 0; i < dimAux ; ++ i) auxK[i] = kOut[i];
  i = 0 ; j = 0 ; idx = 0 ;
  while ( i < dim_kIn && j < dimAux ) {
    if ( fabs ( kIn[i] - auxK[j] ) < EPS ) { // equal knots->get highest multiplicity
      mu   = getMultiplicity(kIn, i, dim_kIn);
      mv   = getMultiplicity(auxK, j, dimAux);
      maxm = mu;
      if (mv > mu) maxm = mv;
        for ( m = 0; m < maxm; ++m)  {
          kOut[idx] = kIn[i];
          ++idx;
        }
    } else if ( kIn[i] > auxK[j] ) { // insert new knot
      	mv = getMultiplicity (auxK, j, dimAux);
        mu = 0;
        for ( m = 0; m < mv; ++m) {
          kOut[idx] = auxK[j];
          ++idx;
        }
    } else {  // Copy old knot
        mu = getMultiplicity (kIn, i, dim_kIn);
        mv = 0;
        for (m = 0; m < mu; ++m) {
          kOut[idx] = kIn[i];
          ++idx;
        }
    }
    i += mu ;
    j += mv ;
  }
  while( i < dim_kIn ) {
    mu = getMultiplicity (kIn, i, dim_kIn);
    for (m = 0; m < mu; ++m) {
        kOut[idx] = kIn[i];
        ++idx;
    }
    i += mu;
  }
  while( j < dimAux ) {
    mv = getMultiplicity (auxK, j, dimAux);
    for (m = 0; m < mv; ++m) {
        kOut[idx] = auxK[j];
        ++idx;
    }
    j += mv;
  }
  *dim_kOut = idx;

  EG_free(auxK);
  return  EGADS_SUCCESS;
}


/* Compares two knot sequences and collects the missing knots for refinement
	- stores new knots or repeated knots (if higher multiplicity).
	- *Udiff pointer should have minimum size = dimUnew
	- Assumes knot sequences are ordered!
*/
template<class T>
static int
findMissingKnots (int outLevel, int dimUold,   T *Uold, int dimUnew, T *Unew,
                  int *dimUdiff, T *Udiff)
{
  int i, j, k, m, mu_old, mu_new, mu_old0, mu_new0, idx;

  mu_old0 = getMultiplicity (Uold, 0, dimUold);
  mu_new0 = getMultiplicity (Unew, 0, dimUnew);
  if (fabs ( Unew[0] - Uold[0] ) > EPS  ) {
    if (outLevel > 0)
      printf(" EG_skinning: First knot and last should be the same\n");
    return EGADS_RANGERR;
  }
  if (mu_old0 != mu_new0) {
    if (outLevel > 0)
      printf(" EG_skinning: (Assuming clamped knots) curves degree must be equal\n");
    return EGADS_GEOMERR;
  }
  i = mu_old0 ; j = mu_new0 ; idx = 0;
  while ( j < dimUnew ) {
    if ( fabs ( Uold[i] - Unew[j] ) < EPS )	{ // if they are the same. check multiplicity
      mu_old = getMultiplicity(Uold, i, dimUold);
      mu_new = getMultiplicity(Unew, j, dimUnew);
      if ( mu_new > mu_old )	{
        for ( k = 0; k < mu_new - mu_old; ++k )	{
          Udiff[idx] = Unew[k];
          ++idx;
        }
      }
      i += mu_old;
      j += mu_new;
    } else if ( Uold[i] > Unew[j] )  {   // insert new knot
      mu_new = getMultiplicity(Unew, j, dimUnew);
      for ( m = 0; m < mu_new ; ++m) {
        Udiff[idx] = Unew[j];
        ++idx;
      }
      j += mu_new;
    } else {  // do nothing
      if ( i < dimUold ) {
        mu_old = getMultiplicity(Uold, i, dimUold);
        i += mu_old;
      }
    }
  }
  *dimUdiff = idx;

  return EGADS_SUCCESS;
}


/* Knot Refinement: Stores new control points (**newPts) and newKnots (**newKnots)
   NewPts collects all the new control points. It should be previously allocated
          and at least dim = (nP+1)*3
   NewKnots collects the new knot sequence. It should be previously allocated
          and at least dim = m+r +2
   ACTUAL lengths are stored in pointers *dimNewKnot and *dimNewCps
*/
template<class T>
static int
RefineKnotVectCurve (int n, int p, T *splineData, int r, T *insertKnots,
                     int *dimNewCps, T *newPts,int *dimNewKnot,
                     T *newKnots)
{
  int d, i, j, k, l, m, a, b, nP, offset, ind;
  T   *splineKnots = NULL, alfa;

  m  = n + p + 1;  // Tot knots in Spline data( acually number is m+1 )
  nP = n + r + 1; // num of new ctrl points ( acually number is nP+1 )
  if ( n < 0  || p < 0 || r < 0 ) return EGADS_INDEXERR;
  if ( (splineData == NULL ) || (insertKnots == NULL ) ) return EGADS_NULLOBJ;
  splineKnots = (T *) EG_alloc( (m+1) *sizeof(T));
  if (splineKnots == NULL ) return EGADS_MALLOC;

  for ( i = 0; i <= m; ++i)  splineKnots[i] = splineData[i];
  a = FindSpan(n, p, insertKnots[0], splineKnots);
  b = FindSpan(n, p, insertKnots[r], splineKnots);
  b++;
  offset = m + 1 ;
  for ( j = 0; j <= a-p; ++j ) {
    for ( d = 0; d < 3; ++d )
      newPts[j*3 + d] = splineData[j*3 +d+offset];  // m is the data offset
  }
  for ( j = b-1; j <= n; ++j) {
    for( d = 0; d < 3; ++d)
      newPts[ (j+r+1)*3 + d ] = splineData[ j*3 +d+offset ];
  }
  for ( j = 0  ; j <= a; ++j )        newKnots[j] = splineKnots[j];
  for ( j = b+p; j <= m; ++j )    newKnots[j+r+1] = splineKnots[j];
  i = b + p - 1;  k = b + p + r;
  for ( j = r; j >= 0; --j ) {
    while ( ( insertKnots[j] <= splineKnots[i] ) && ( i > a ) ) {
      for ( d = 0; d < 3; ++d )
        newPts[ (k-p-1)*3 + d ] = splineData[ (i-p-1)*3 +d+offset ];
      newKnots[k] = splineKnots[i];
      --k;
      --i;
    }
    for ( d = 0; d < 3; ++d )   newPts[ (k-p-1)*3 + d ] = newPts[ (k-p)*3 +d ];
    for ( l = 1; l <= p; ++l ) {
      ind = k - p + l;
      alfa = newKnots[k+l] - insertKnots[j];
      if ( fabs(alfa) == 0. ) {
        for ( d = 0; d < 3; ++d )
          newPts[ (ind -1)*3+d ] = newPts[ ind*3+d ];
      } else {
        alfa /= ( newKnots[k+l] - splineKnots[ i-p+l ] );
        for ( d = 0; d < 3; ++d )
          newPts[ (ind-1)*3+d ] = alfa*newPts[ (ind-1)*3+d ] +(1.0-alfa)*newPts[ ind*3+d ];
      }
    }
    newKnots[k] = insertKnots[j];
    --k;
  }
  *dimNewCps  = nP + 1;
  *dimNewKnot = (m+1) + (r+1);
  EG_free(splineKnots);

  return EGADS_SUCCESS;
}

/*  Assuming all curves have equal degree, makes the curves compatible:
	- creates a common (unique) knot vector.
	- creates new control points for each curve.
*/
template<class T>
static int
makeCurvesCompatible(int nC, ego *splineCurves, int *nP, T **controlPoints,
                     T **knotVector, int *dimKnotVector, int *degree)
{
  int c = 0, k = 0, p = 0, d = 0, oclass, ctype, stat, offset, dimMergedKnots;
  int dimKnotIsert = 0, dimKnotCheck = 0, nPts = 0,  deg = 0, knotSum = 0;
  int **splineInfo = NULL;
  T   **splineData = NULL, *cp = NULL, *mergedKnots = NULL;
  T   *knotInsert = NULL,  *knotCheck = NULL, *newCtrlPts = NULL;
  ego geom;

  int outLevel = EG_outLevel(splineCurves[0]);

  *controlPoints = NULL;
  *knotVector    = NULL;
  if ( nC <= 0 ) return EGADS_INDEXERR;
  splineData = (T **) EG_alloc(nC*sizeof(T *));
  if ( splineData == NULL ) return  EGADS_MALLOC;
  splineInfo = (int    **) EG_alloc(nC*sizeof(int *));
  if ( splineInfo == NULL )  {
    EG_free(splineData);
    return  EGADS_MALLOC;
  }

  for ( c = 0; c < nC; ++c ) {
    stat = EG_getGeometry(splineCurves[c], &oclass, &ctype, &geom,
                          &splineInfo[c], &splineData[c]);
    if ( stat != EGADS_SUCCESS ) {
      if (outLevel > 0)
        printf(" EG_getGeometry from curve %d  = %d\n", c, stat);
      goto bail;
    }
    if ( oclass != CURVE ) {
      if (outLevel > 0)
        printf(" EG_skinning: Curve %d not a Curve!\n", c);
      stat = EGADS_NOTGEOM;
      goto bail;
    }
    if ( ctype != BSPLINE ) {
      if (outLevel > 0)
        printf(" EG_skinning: Curve %d not a BSpline!\n", c);
      stat = EGADS_NOTGEOM;
      goto bail;
    }
    if ( ( splineInfo[c][0] & 2 ) != 0 ) {
      if (outLevel > 0)
        printf(" EG_skinning: Curve %d is Rational!\n", c);
      stat = EGADS_NOTGEOM;
      goto bail;
    }
  }
  // Check Degree: must be the same (ATM)
  deg = splineInfo[0][1];
  for (c = 1; c < nC; ++c) {
    if (deg != splineInfo[c][1])    {
      if (outLevel > 0)
        printf(" EG_skinning: Curves need to be of same degree!\n");
      stat = EGADS_NOTGEOM;
      goto bail;
    }
  }
  nPts           = splineInfo[0][2];
  dimMergedKnots = splineInfo[0][3];
  stat = 0;
  // Check if knot sequences are different
  for ( c = 1; c < nC; ++c ) {
    if ( nPts != splineInfo[c][2] ) {
      stat = 1;
      break;
    }
    if ( dimMergedKnots != splineInfo[c][3] ) {
      stat = 1;
      break;
    }
    for ( k = 0; k < dimMergedKnots; ++k ) {
      if ( splineData[c][k] != splineData[0][k] ) {
        stat = 1;
        break;
      }
    }
    if ( stat == 1 ) break;
  }
  if (stat == 1) {   // Create Common Knot Sequence
    knotSum = 0;
    for ( c = 0 ; c < nC ; ++c )
      knotSum += splineInfo[c][3]; // maximum size= all knots are different
    mergedKnots = (T *) EG_alloc(knotSum *sizeof(T));
    if ( mergedKnots == NULL ) goto bail;
    dimMergedKnots = 0;
    for ( c = 0 ; c < nC ; ++c ) {
      stat = mergeKnotVectors( splineInfo[c][3], splineData[c],
                               &dimMergedKnots, mergedKnots);
      if ( stat != EGADS_SUCCESS)  {
        if (outLevel > 0)
          printf(" EG_skinning: mergeKnotVectors for curve %d   = %d\n", c, stat);
        goto bail;
      }
    }
    knotSum = 0;
    for ( c = 0 ; c < nC ; ++c )
      knotSum += splineInfo[c][3]; // maximum size= all knots are different
    cp         = (T*) EG_alloc(knotSum*3*nC  *sizeof(T));  // max size of control points
    knotInsert = (T*) EG_alloc(dimMergedKnots*sizeof(T));
    knotCheck  = (T*) EG_alloc(knotSum 	     *sizeof(T));
    newCtrlPts = (T*) EG_alloc(knotSum*3     *sizeof(T));
    if ( cp == NULL  || knotInsert == NULL || newCtrlPts == NULL ||
        knotCheck == NULL )  {
      stat = EGADS_MALLOC;
      goto bail;
    }
    dimKnotIsert = 0;
    for ( c = 0; c < nC; ++c ) {
      stat = findMissingKnots(outLevel, splineInfo[c][3], splineData[c], dimMergedKnots, mergedKnots,
                              &dimKnotIsert, knotInsert);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EG_skinning: findMissingKnots %d = %d\n", c, stat);
        goto bail;
      }
      if ( dimKnotIsert > 0 ) {
        // there are new knots in this sequence: collect knots and control points.
        stat = RefineKnotVectCurve (splineInfo[c][2]-1, deg, splineData[c],
                                    dimKnotIsert-1, knotInsert, &nPts, newCtrlPts,
                                    &dimKnotCheck, knotCheck);
        if (stat != EGADS_SUCCESS) goto bail;
        // knotcheck is a vector which should coincide with the mergedKnot sequence that we found
        if (dimKnotCheck != dimMergedKnots)	{
          if (outLevel > 0)
            printf(" EG_skinning: Knot Refinement has the wrong dimension!!\n");
          stat = EGADS_GEOMERR;
          goto bail;
        }
        for( k = 0; k< dimKnotCheck; ++k) {
          if (fabs( knotCheck[k] - mergedKnots[k] ) > EPS) {
            if (outLevel > 0)
              printf(" EG_skinning: Knot Refinement has the wrong sequence!!\n");
            stat = EGADS_GEOMERR;
            goto bail;
          }
        }
      } else {
        // no new knots for this curve so copy directly from spline data
        nPts   = splineInfo[c][2];
        offset = splineInfo[c][3];
        for ( p = 0; p < nPts; ++p ) {
          for ( d = 0 ; d < 3 ; ++ d)
            newCtrlPts[ p*3 + d ] = splineData[c][ offset+d+p*3 ];
        }
      }  // add curve control points to global pointer cp
      for ( p = 0; p < nPts ; ++p ) {
        for ( d = 0; d < 3; ++d )
          cp[ d + p*3 + c*nPts*3 ] = newCtrlPts[ d + p*3 ];
      }
    }
  } else {
    // All knots are the same. Copy directly from spline data
    deg            = splineInfo[0][1];
    nPts           = splineInfo[0][2];
    dimMergedKnots = splineInfo[0][3];
    mergedKnots    = (T*) EG_alloc(dimMergedKnots*sizeof(T));
    cp             = (T*) EG_alloc(nC*nPts*3     *sizeof(T));  // max size of control points
    if ( mergedKnots == NULL || cp == NULL ) {
      stat = EGADS_MALLOC;
      goto bail;
    }
    for ( k = 0; k < dimMergedKnots; ++k ) mergedKnots[k] = splineData[0][k];
    for ( c = 0; c < nC; ++c ) {
      offset = dimMergedKnots;
      for ( p = 0; p < nPts ; ++p ) {
        for ( d = 0; d < 3; ++d )
          cp[ d + p*3 + c*nPts*3 ] = splineData[c][offset+d];
        offset+=3;
      }
    }
  }
  // Copy all data to input variables & leave
  *degree        = deg;  // Assuming that one day degree elevation function will be implemented, *degree will update
  *dimKnotVector = dimMergedKnots;
  *knotVector    = (T *) EG_alloc( dimMergedKnots*sizeof(T));
  if ( *knotVector == NULL ) {
    stat = EGADS_MALLOC;
    goto bail;
  }
  for ( k = 0; k < dimMergedKnots ; ++k )
    (*knotVector)[k] = mergedKnots[k];
  *controlPoints = (T *) EG_alloc(nPts*nC*3*sizeof(T));
  if( *controlPoints == NULL ) {
    stat = EGADS_MALLOC;
    goto bail ;
  }
  *nP = nPts ;
  for (c = 0; c < nC; ++c) {
    for ( p = 0;p < nPts ; ++p) {
      for (d = 0 ; d < 3 ; ++ d)
        (*controlPoints)[ d + p*3 + c*nPts*3 ] = cp[ d + p*3 + c*nPts*3 ];
    }
  }

bail:
  for (c = 0; c < nC; ++c) {
    EG_free(splineInfo[c]);
    EG_free(splineData[c]);
  }
  EG_free(splineData);
  EG_free(splineInfo);
  EG_free(mergedKnots);
  EG_free(cp);
  EG_free(knotCheck);
  EG_free(knotInsert);
  EG_free(newCtrlPts);

  return stat;
}

template<class T>
int makeSkinningGeom(int nC, ego *sectionCurves, int skinning_degree,
                     int *splineInfo, T **splineData)
{
  T   *P_cross, *Uknots; // ptr to control points and commont knot vector
  int c = 0, d = 0, i = 0, j = 0, n=0, it = 0, stat, dimVknots = 0;
  int nP = 0, dimUknots = 0, degree = 0;
  int blend = 0, dataLength = 0, offset = 0;
  T   sumVal, diam, dist, locCord;
  T   *A = NULL, *A_aux = NULL, *net[3] , *b = NULL;
  T   *v_param = NULL , *vKnots = NULL;

  int outLevel = EG_outLevel(sectionCurves[0]);

  (*splineData) = NULL;

  stat = makeCurvesCompatible<T>(nC, sectionCurves, &nP, &P_cross, &Uknots,
                              &dimUknots, &degree);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeCurvesCompatible = %d (EG_skinning)!\n", stat);
    return stat;
  }

  // Get parameters for skinning:
  v_param = (T *) EG_alloc(nC *sizeof(T));
  if ( v_param == NULL ) return EGADS_MALLOC;
  v_param[0]    = 0.0;
  v_param[nC-1] = 1.0;
  for ( c = 1; c < nC-1; ++c ) {
    sumVal = 0.0;
    for ( n = 0; n < nP; ++n ) {
      diam = 0.0;
      for ( i = 1; i < nC; ++i ) {
        dist=0.0;
        for ( d = 0; d < 3; ++d )
          dist += pow(P_cross[d+n*3+i*3*nP]-P_cross[d+n*3+(i-1)*3*nP], 2);
        diam += sqrt(dist);
      }
      locCord = 0.0;
      for( d = 0; d < 3; ++d )
        locCord += pow(P_cross[d+n*3+c*3*nP]-P_cross[d+n*3+(c-1)*3*nP], 2);
      sumVal += sqrt( locCord )/diam;
    }
    v_param[c] = v_param[c-1] + ( 1./(double)nP )*sumVal;
  }
  dimVknots = nC + skinning_degree +1;
  vKnots    = (T *) EG_alloc(dimVknots *sizeof(T));
  if ( vKnots == NULL ) {
    EG_free(v_param);
    EG_free(P_cross);
    EG_free(Uknots);
    return EGADS_MALLOC;
  }
  for ( i = 0; i <= skinning_degree; ++i ) {
    vKnots[i] = 0.0;
    vKnots[ dimVknots-1-skinning_degree+i ] = 1.0;
  }
  for (j = 1; j < nC-skinning_degree; ++j ) {
    vKnots[ j+skinning_degree ] = 0.0;
    for ( i = j; i <= j+skinning_degree-1; ++i )
      vKnots[j+skinning_degree] += v_param[i];
    vKnots[skinning_degree+j]  /= (double) skinning_degree;
  }
  // FIND SURFACE CONTROL NET ----> INTERPOLATING CROSS CURVES + SOLVING LINEAR SYSTEMS
  A      = (T *) EG_alloc(nC*nC*sizeof(T));
  A_aux  = (T *) EG_alloc(nC*nC*sizeof(T));
  net[0] = (T *) EG_alloc(nP*nC*sizeof(T));
  net[1] = (T *) EG_alloc(nP*nC*sizeof(T));
  net[2] = (T *) EG_alloc(nP*nC*sizeof(T));
  b      = (T *) EG_alloc(   nC*sizeof(T));
  if ( (A == NULL) || (A_aux == NULL) || (net[0] == NULL) || (net[1] == NULL) ||
       (net[2] == NULL) || (b == NULL) ) {
    if (b      != NULL) EG_free(b);
    if (net[2] != NULL) EG_free(net[2]);
    if (net[1] != NULL) EG_free(net[1]);
    if (net[0] != NULL) EG_free(net[0]);
    if (A_aux  != NULL) EG_free(A_aux);
    if (A      != NULL) EG_free(A);
    EG_free(vKnots);
    EG_free(v_param);
    EG_free(P_cross);
    EG_free(Uknots);
    return EGADS_MALLOC;
  }
  for ( it = i = 0; i < nC; ++i )
    for ( j = 0; j < nC; ++j, ++it )
      A[it] = OneBasisFun(skinning_degree, dimVknots-1, vKnots, j, v_param[i]);

  for ( blend = 0; blend < nP; ++blend ) {
    for ( d = 0; d < 3; ++d ) {
      for ( n = 0; n < nC;    ++n ) b[n] = P_cross[ d+blend*3+n*3*nP ];
      for ( n = 0; n < nC*nC; ++n ) A_aux[n] = A[n];
      // Finds the Control Points such that the curve interpolates the original control points.
      stat = matsol(A_aux, b, nC, &net[d][blend*nC]);
      if (stat != EGADS_SUCCESS) {
        if (outLevel > 0)
          printf(" EGADS Error: Solving Linear System = %d!!\n", stat);
        EG_free(b);
        EG_free(net[2]);
        EG_free(net[1]);
        EG_free(net[0]);
        EG_free(A_aux);
        EG_free(A);
        EG_free(vKnots);
        EG_free(v_param);
        EG_free(P_cross);
        EG_free(Uknots);
        return stat;
      }
    }
  }
  EG_free(b);
  EG_free(A_aux);
  EG_free(A);
  EG_free(v_param);
  EG_free(P_cross);

  /* SURFACE BSPLINES */
  splineInfo[0] = 0;
  splineInfo[1] = degree;
  splineInfo[2] = nP;
  splineInfo[3] = dimUknots; // TOTAL KNOTS IN CROSS SECTIONS
  splineInfo[4] = skinning_degree;
  splineInfo[5] = nC;
  splineInfo[6] = dimVknots; // TOTAL KNOTS IN BLENDER FUNCTIONS
  dataLength    = splineInfo[3] + nC*nP*3 + splineInfo[6];
  *splineData   = (T *) EG_alloc(dataLength *sizeof(T));
  if (*splineData == NULL ) {
    EG_free(net[2]);
    EG_free(net[1]);
    EG_free(net[0]);
    EG_free(vKnots);
    EG_free(Uknots);
    return EGADS_MALLOC;
  }

  /* KNOT SEQUENCE  (cross-sections direction) */
  for ( i = 0; i < dimUknots; ++i ) (*splineData)[i] = Uknots[i];
  EG_free(Uknots);
  /* KNOT SEQUENCE  (skinning direction) */
  for ( i = 0; i < dimVknots; ++i ) (*splineData)[dimUknots+i] = vKnots[i];
  EG_free(vKnots);
  /* CONTOL NET */
  offset = dimVknots + dimUknots;
  for (c = 0; c < nC; ++c)
    for (n = 0; n < nP; ++n) {
      for (d = 0; d < 3; d++) (*splineData)[offset+d] = net[d][n*nC+c];
      offset+=3;
    }
  EG_free(net[2]);
  EG_free(net[1]);
  EG_free(net[0]);

  return EGADS_SUCCESS;
}


extern "C"
#ifdef STANDALONE
static
#endif
int EG_skinning(int nC, ego *sectionCurves, int skinning_degree, ego *surface)
{
  int    stat, i, outLevel, data_dot=1;
  int    splineInfo[7];
  ego    context;

  *surface = NULL;
  if (nC <= 0)                                return EGADS_RANGERR;
  for (i = 0; i < nC; i++) {
    if (sectionCurves[i] == NULL)               return EGADS_NULLOBJ;
    if (sectionCurves[i]->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if (sectionCurves[i]->oclass != CURVE)      return EGADS_NOTGEOM;
    if (EG_sameThread(sectionCurves[i]))        return EGADS_CNTXTHRD;
    if (EG_hasGeometry_dot(sectionCurves[i]) != EGADS_SUCCESS) data_dot=0;
  }
  stat = EG_getContext(sectionCurves[0], &context);
  if (stat != EGADS_SUCCESS) return stat;
  outLevel = EG_outLevel(sectionCurves[0]);

  if (skinning_degree > nC-1) {  // max spline degree allowed in v-direction is nC
    if (outLevel > 0)
      printf(" EGADS Warning: Reducing skinning degree to %d from %d (EG_skinning)\n",
             nC-1, skinning_degree);
    skinning_degree = nC-1;
  }

  if (data_dot == 1) {
    SurrealS<1> *splineData = NULL;

    stat = makeSkinningGeom< SurrealS<1> >(nC, sectionCurves, skinning_degree,
                                           splineInfo, &splineData);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeSkinningGeom = %d (EG_skinning)!\n", stat);
      return stat;
    }

    stat = EG_makeGeometry(context, SURFACE, BSPLINE,
                           NULL, splineInfo, splineData, surface);
    EG_free(splineData);
    if ( stat != EGADS_SUCCESS ) {
      if (outLevel > 0)
        printf(" EG_skinning: EG_makeGeometry SPLINE SURFACE = %d\n", stat);
      return stat;
    }
  } else {
    double *splineData = NULL;

    stat = makeSkinningGeom<double>(nC, sectionCurves, skinning_degree,
                                    splineInfo, &splineData);
    if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: makeSkinningGeom = %d (EG_skinning)!\n", stat);
      return stat;
    }

    stat = EG_makeGeometry(context, SURFACE, BSPLINE,
                           NULL, splineInfo, splineData, surface);
    EG_free(splineData);
    if ( stat != EGADS_SUCCESS ) {
      if (outLevel > 0)
        printf(" EG_skinning: EG_makeGeometry SPLINE SURFACE = %d\n", stat);
      return stat;
    }
  }
  return EGADS_SUCCESS;
}


extern "C"
#ifdef STANDALONE
static
#endif
int EG_skinning_dot(ego surface, int nC, ego *sectionCurves)
{
  int    stat, outLevel, i, skinning_degree;
  int    splineInfo[7], oclass, mtype, *ivec=NULL;
  double *rvec=NULL;
  SurrealS<1> *splineData = NULL;
  ego refGeom;

  if (nC <= 0)                                  return EGADS_RANGERR;
  for (i = 0; i < nC; i++) {
    if (sectionCurves[i] == NULL)               return EGADS_NULLOBJ;
    if (sectionCurves[i]->magicnumber != MAGIC) return EGADS_NOTOBJ;
    if (sectionCurves[i]->oclass != CURVE)      return EGADS_NOTGEOM;
    if (EG_sameThread(sectionCurves[i]))        return EGADS_CNTXTHRD;
    if (EG_hasGeometry_dot(sectionCurves[i]) != EGADS_SUCCESS) {
      outLevel = EG_outLevel(sectionCurves[0]);
      if (outLevel > 0)
        printf(" EGADS Error: Curve[%d] does not have sensitivities (EG_skinning_dot)!\n", i);
      return EGADS_NODATA;
    }
  }
  outLevel = EG_outLevel(sectionCurves[0]);


  stat = EG_getGeometry(surface, &oclass, &mtype, &refGeom, &ivec, &rvec);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_getGeometry = %d (EG_skinning_dot)!\n", stat);
    return stat;
  }

  skinning_degree = ivec[4];
  EG_free(ivec); ivec=NULL;
  EG_free(rvec); rvec=NULL;

  stat = makeSkinningGeom< SurrealS<1> >(nC, sectionCurves, skinning_degree,
                                  splineInfo, &splineData);
  if (stat != EGADS_SUCCESS) {
    if (outLevel > 0)
      printf(" EGADS Error: makeSkinningGeom = %d (EG_skinning_dot)!\n", stat);
    return stat;
  }

  stat = EG_setGeometry_dot(surface, SURFACE, BSPLINE,
                            splineInfo, splineData);
  EG_free(splineData);
  if ( stat != EGADS_SUCCESS ) {
    if (outLevel > 0)
      printf(" EGADS Error: EG_setGeometry_dot = %d (EG_skinning_dot)!\n", stat);
    return stat;
  }
  return EGADS_SUCCESS;
}


#ifdef STANDALONE
// This function "manufactures" cross-section curves (sine-waves) for testing
static int
createCrossSections(ego context, ego *splineCurves, int nC, int degree)
{
  int    c, i, j, k, n, stat, nCP, Nknots, offset;
  int    splineHeader[4];
  double t, **data;

  data = (double **) EG_alloc(nC*sizeof(double*));
  if (data == NULL) return EGADS_MALLOC;
  splineHeader[0] = 0;
  splineHeader[1] = degree;
  // Control Points
  for (c = 0; c < nC; ++c) {
    // Uniform knots (clamped sequence)
    nCP             = degree+1+nC-c; // each spline curve a different knot seq.
    Nknots          = nCP + degree+1;
    splineHeader[2] = nCP;
    splineHeader[3] = Nknots;
    data[c]         = (double*) EG_alloc((Nknots + nCP*3) *sizeof(double));
    if(data[c] == NULL ) {
      EG_free(data);
      return EGADS_MALLOC;
    }
    for (k = 0; k <= degree; ++k) {
      data[c][k] = 0.0;
      data[c][Nknots-1-degree+k] = 1.0;
    }
    for (j = 1; j < nCP-degree; ++j) {
      data[c][j+degree] = 0.0;
      for (i = j; i <= j+degree-1; ++i)
        data[c][j+degree] = (double)j/(nCP-degree);
    }
    offset = Nknots;
    for (n = 0; n < nCP; ++n) {
      t = (double)n/(double)nCP;
      data[c][offset]   = t;                //x-coord
      data[c][offset+1] = sin(4*PI*t);      //y-coord
      data[c][offset+2] = (double)c/(double)nC;
      offset	       +=3;
    }
    stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL, splineHeader, data[c],
                           &splineCurves[c]);
    if (stat != EGADS_SUCCESS) {
      for (i = 0; i < nC; ++i) EG_free(data[i]);
      EG_free(data);
      return stat;
    }
  }
  for (i = 0; i < nC; ++i) EG_free(data[i]);
  EG_free(data);

  return EGADS_SUCCESS;
}


int main(int argc, char *argv[])
{
  int c, stat;
  int nC              = 9; // NUMBER OF CURVES
  int degree          = 3; // Splines Degree
  int skinning_degree = 3; // nC - 2;
  int sense;
  ego *sectionEdges = NULL , *loop= NULL, objs[2], nodes[2];
  ego face, *body = NULL, model ;
  double t, range[2], limits[4], xyz[9];
  ego context, *sectionCurves, surface;
  if ( argc != 4 && argc != 1 ) {
    printf("\n Usage: makeSkinning #Curves Spline-degree Surface-degree\n\n");
    return 1;
  }
  if (argc == 4) {
    nC              = atoi(argv[1]);
    degree          = atoi(argv[2]);
    skinning_degree = atoi(argv[3]);
  }
  sectionCurves = (ego *) EG_alloc(nC*sizeof(ego));
  if ( sectionCurves == NULL ) goto bail;
  printf(" EG_open            = %d\n", EG_open(&context));
  stat = createCrossSections(context, sectionCurves, nC, degree);
  if ( stat != EGADS_SUCCESS ) {
    printf(" CREATE CROSS SECTIONS             = %d\n", stat);
    goto bail;
  }
  stat = EG_skinning(nC, sectionCurves, skinning_degree, &surface);
  if(stat != EGADS_SUCCESS) {
    printf(" EG_surfaceSkinning     = %d\n", stat);
    goto bail;
  }
  if ( stat == EGADS_SUCCESS ) {
    // SURFACE LIMITS: TAKE FIRST AND LAST MEMBERS IN KNOT SEQUEnCE.
    limits[0] = 0.0;
    limits[1] = 1.0;
    limits[2] = 0.0;
    limits[3] = 1.0;
    stat =  EG_makeFace( surface, SFORWARD, limits, &face);
    if ( stat != EGADS_SUCCESS ) {
      printf(" EG_makeFace  SPLINE face     = %d\n", stat);
      goto bail;
    }
    /* NOW WE CONSTRUCT EDGES MADE OF THE BSPLINE CROSS SECTION CURVES.
       IF SKINNING WORKED, SURFACE SHOULD INTERSECT THESE CURVES */
    sectionEdges = (ego *) EG_alloc(nC *sizeof(ego));
    loop = (ego *) EG_alloc(nC *sizeof(ego));
    body = (ego *) EG_alloc( (nC+1) *sizeof(ego));
    if ( sectionEdges == NULL || body == NULL || loop == NULL ) {
      printf(" NULL POINTERS IN MAIN\n ");
      goto bail;
    }
    for ( c = 0; c < nC; ++c ) {
      t    = 0.0;
      stat = EG_evaluate(sectionCurves[c], &t, xyz);
      if(stat != EGADS_SUCCESS) {
        printf(" EG_evaluate on curve %d,  NODE 0  = %d\n", c, stat);
        goto bail;
      }
      stat = EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL,
                             &nodes[0]);
      if(stat != EGADS_SUCCESS) {
        printf(" EG_makeTopology on curve %d,  NODE 0  = %d\n",c,stat);
        goto bail;
      }
      t    = 1.0;
      stat = EG_evaluate(sectionCurves[c], &t, xyz);
      if(stat != EGADS_SUCCESS) {
        printf(" EG_evaluate on curve %d,  NODE 1  = %d\n", c, stat);
        goto bail;
      }
      stat= EG_makeTopology(context, NULL, NODE, 0, xyz, 0, NULL, NULL,
                            &nodes[1]);
      if(stat != EGADS_SUCCESS) {
        printf(" EG_makeTopology on curve %d,  NODE 1 = %d\n", c, stat);
        goto bail;
      }
      objs[0]  = nodes[0];
      objs[1]  = nodes[1];
      range[0] = limits[0];
      range[1] = limits[1];
      stat = EG_makeTopology(context, sectionCurves[c], EDGE,
                             TWONODE, range, 2, objs, NULL, &sectionEdges[c]);
      if( stat != EGADS_SUCCESS ) {
        printf(" EG_makeTopology for spline curve  %d  =  %d\n",c,stat);
        goto bail;
      }
      sense = 1;
      stat  = EG_makeTopology(context, NULL, LOOP, OPEN, NULL, 1,
                               &sectionEdges[c], &sense, &loop[c]);
      if( stat != EGADS_SUCCESS ) {
        printf(" EG_makeTopology for spline Loop %d  =  %d\n", c, stat);
        goto bail;
      }
      stat = EG_makeTopology(context, NULL, BODY, WIREBODY, NULL, 1, &loop[c],
                             NULL, &body[c]);
      if(stat != EGADS_SUCCESS) {
        printf(" EG_makeTopology for WIREBODY   %d  =  %d\n",c,stat);
        goto bail;
      }
    }
    stat = EG_makeTopology(context, NULL, BODY, FACEBODY, NULL, 1, &face,
                           NULL, &body[nC]);
    if ( stat != EGADS_SUCCESS ) {
      printf(" EG_makeTopology for FACE BODY  =  %d\n", stat);
      goto bail;
    }
    stat = EG_makeTopology(context, NULL, MODEL, 0, NULL, nC+1, body, NULL,
                           &model);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeTopology for MODEL  = %d\n", stat);
      goto bail;
    }
    printf(" EG_saveModel       = %d\n", EG_saveModel(model, "skinning.egads"));

bail:
    if (stat != EGADS_SUCCESS) printf("Error Skinning Process: %d\n", stat);
  }
  EG_free(loop);
  EG_free(body);
  EG_free(sectionEdges);
  EG_free(sectionCurves);
  EG_close(context);
  return 0;
}
#endif
