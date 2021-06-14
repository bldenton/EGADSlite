#include <stdlib.h>
#include <string.h>


#ifdef LITE
#ifdef WIN32
#define false 0
#define true  (!false)
typedef int bool;
#else
#include <stdbool.h>
#endif
#define TEMPLATE
#define DOUBLE double
#else
#include "Surreal/SurrealS.h"
#define TEMPLATE template<class TT>
#define DOUBLE TT
#endif

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif


/*
 * These routines were taken from the OpenNURBS package with the date of
 * 2013-07-11, see <http://www.openNURBS.org>.
 */


__HOST_AND_DEVICE__ static double
BinomialCoefficient(int i, int j)
{
#define MAX_HALF_N 26
  static const double bc[((MAX_HALF_N-2)*(MAX_HALF_N-1))/2 + MAX_HALF_N - 2] =
  {15.0, 20.0, 28.0, 56.0, 70.0, 45.0, 120.0, 210.0, 252.0, 66.0,
    220.0, 495.0, 792.0, 924.0, 91.0, 364.0, 1001.0, 2002.0, 3003.0,
    3432.0, 120.0, 560.0, 1820.0, 4368.0, 8008.0, 11440.0, 12870.0,
    153.0, 816.0, 3060.0, 8568.0, 18564.0, 31824.0, 43758.0, 48620.0,
    190.0, 1140.0, 4845.0, 15504.0, 38760.0, 77520.0, 125970.0,
    167960.0, 184756.0, 231.0, 1540.0, 7315.0, 26334.0, 74613.0,
    170544.0, 319770.0, 497420.0, 646646.0, 705432.0, 276.0, 2024.0,
    10626.0, 42504.0, 134596.0, 346104.0, 735471.0, 1307504.0,
    1961256.0, 2496144.0, 2704156.0, 325.0, 2600.0, 14950.0, 65780.0,
    230230.0, 657800.0, 1562275.0, 3124550.0, 5311735.0, 7726160.0,
    9657700.0, 10400600.0, 378.0, 3276.0, 20475.0, 98280.0, 376740.0,
    1184040.0, 3108105.0, 6906900.0, 13123110.0, 21474180.0,
    30421755.0, 37442160.0, 40116600.0, 435.0, 4060.0, 27405.0,
    142506.0, 593775.0, 2035800.0, 5852925.0, 14307150.0, 30045015.0,
    54627300.0, 86493225.0, 119759850.0, 145422675.0, 155117520.0,
    496.0, 4960.0, 35960.0, 201376.0, 906192.0, 3365856.0,
    10518300.0, 28048800.0, 64512240.0, 129024480.0, 225792840.0,
    347373600.0, 471435600.0, 565722720.0, 601080390.0, 561.0,
    5984.0, 46376.0, 278256.0, 1344904.0, 5379616.0, 18156204.0,
    52451256.0, 131128140.0, 286097760.0, 548354040.0, 927983760.0,
    1391975640.0, 1855967520.0, 2203961430.0, 2333606220.0, 630.0,
    7140.0, 58905.0, 376992.0, 1947792.0, 8347680.0, 30260340.0,
    94143280.0, 254186856.0, 600805296.0, 1251677700.0, 2310789600.0,
    3796297200.0, 5567902560.0, 7307872110.0, 8597496600.0,
    9075135300.0, 703.0, 8436.0, 73815.0, 501942.0, 2760681.0,
    12620256.0, 48903492.0, 163011640.0, 472733756.0, 1203322288.0,
    2707475148.0, 5414950296.0, 9669554100.0, 15471286560.0,
    22239974430.0, 28781143380.0, 33578000610.0, 35345263800.0,
    780.0, 9880.0, 91390.0, 658008.0, 3838380.0, 18643560.0,
    76904685.0, 273438880.0, 847660528.0, 2311801440.0, 5586853480.0,
    12033222880.0, 23206929840.0, 40225345056.0, 62852101650.0,
    88732378800.0, 113380261800.0, 131282408400.0, 137846528820.0,
    861.0, 11480.0, 111930.0, 850668.0, 5245786.0, 26978328.0,
    118030185.0, 445891810.0, 1471442973.0, 4280561376.0,
    11058116888.0, 25518731280.0, 52860229080.0, 98672427616.0,
    166509721602.0, 254661927156.0, 353697121050.0, 446775310800.0,
    513791607420.0, 538257874440.0, 946.0, 13244.0, 135751.0,
    1086008.0, 7059052.0, 38320568.0, 177232627.0, 708930508.0,
    2481256778.0, 7669339132.0, 21090682613.0, 51915526432.0,
    114955808528.0, 229911617056.0, 416714805914.0, 686353797976.0,
    1029530696964.0, 1408831480056.0, 1761039350070.0,
    2012616400080.0, 2104098963720.0, 1035.0, 15180.0, 163185.0,
    1370754.0, 9366819.0, 53524680.0, 260932815.0, 1101716330.0,
    4076350421.0, 13340783196.0, 38910617655.0, 101766230790.0,
    239877544005.0, 511738760544.0, 991493848554.0, 1749695026860.0,
    2818953098830.0, 4154246671960.0, 5608233007146.0,
    6943526580276.0, 7890371113950.0, 8233430727600.0, 1128.0,
    17296.0, 194580.0, 1712304.0, 12271512.0, 73629072.0,
    377348994.0, 1677106640.0, 6540715896.0, 22595200368.0,
    69668534468.0, 192928249296.0, 482320623240.0, 1093260079344.0,
    2254848913647.0, 4244421484512.0, 7309837001104.0,
    11541847896480.0, 16735679449896.0, 22314239266528.0,
    27385657281648.0, 30957699535776.0, 32247603683100.0, 1225.0,
    19600.0, 230300.0, 2118760.0, 15890700.0, 99884400.0,
    536878650.0, 2505433700.0, 10272278170.0, 37353738800.0,
    121399651100.0, 354860518600.0, 937845656300.0, 2250829575120.0,
    4923689695575.0, 9847379391150.0, 18053528883775.0,
    30405943383200.0, 47129212243960.0, 67327446062800.0,
    88749815264600.0, 108043253365600.0, 121548660036300.0,
    126410606437752.0, 1326.0, 22100.0, 270725.0, 2598960.0,
    20358520.0, 133784560.0, 752538150.0, 3679075400.0,
    15820024220.0, 60403728840.0, 206379406870.0, 635013559600.0,
    1768966344600.0, 4481381406320.0, 10363194502115.0,
    21945588357420.0, 42671977361650.0, 76360380541900.0,
    125994627894135.0, 191991813933920.0, 270533919634160.0,
    352870329957600.0, 426384982032100.0, 477551179875952.0,
    495918532948104.0};

  int n, half_n, bc_i;

  if (i  < 0 || j  < 0) return  0.0;
  if (0 == i || 0 == j) return  1.0;
  n = i+j;
  if (1 == i || 1 == j) return (double)n;
  if (4 == n)           return  6.0;
  if (5 == n)           return 10.0;

  if (n%2) return BinomialCoefficient(i-1,j)+BinomialCoefficient(i,j-1);

  half_n = n >> 1;
  if (half_n > MAX_HALF_N)
           return BinomialCoefficient(i-1,j)+BinomialCoefficient(i,j-1);

  if (i > half_n) i = n - i;
 /* at this point we have n even,
  * MAX_HALF_N*2 >= n >= 6 and 1 < i <= n/2
  * and we grab the answer from the bc[] table.
  */
  half_n -= 2;
  bc_i = ((half_n*(half_n+1))>>1) + i - 3;
  return bc[bc_i];

#undef MAX_HALF_N
}


/*
Description:
  Use the quotient rule to compute derivatives of a one parameter
  rational function F(t) = X(t)/W(t), where W is a scalar
  and F and X are vectors of dimension dim.
Parameters:
  dim       - [in]
  der_count - [in] number of derivative (>=0)
  v_stride  - [in] (>= dim+1)
  v         - [in/out]
    v[] is an array of length (der_count+1)*v_stride.
    The input v[] array contains  derivatives of the numerator and
    denominator	functions in the order (X, W), (Xt, Wt), (Xtt, Wtt), ...
    In general, the (dim+1) coordinates of the d-th derivative
    are in (v[n],...,v[n+dim]) where n = d*v_stride.
    In the output v[] array the derivatives of X are replaced with
    the derivatives of F and the derivatives of W are divided by
    w = v[dim].
*/

__HOST_AND_DEVICE__ TEMPLATE
void EG_EvaluateQuotientRule(int dim, int der_count, int v_stride, DOUBLE *v)
{
  /*
    The quotient rule says the n-th derivative is

     (n)      (n)          (n)             (n-1)    (1)              (1)    (n-1)
    f  (t) = x   (t)  -  (w  (t)*f(t) + n*w    (t)*f  (t) + ... + n*w  (t)*f    (t))
             ---------------------------------------------------------------------
                                            w(t)

                                                                 (i)   (j)
      (The missing summands look like  BinomialCoefficient(i,j)*w   * f   )
  */
#ifndef __clang_analyzer__
  DOUBLE wt, w2, *f, *x, *w;
  int    i, j, n, df;

  wt = v[dim];
  wt = 1.0/wt;
  i  = (der_count+1)*v_stride;
  x  = v;
  while (i--) *x++ *= wt;

  if (der_count) {
    // 1st derivative - faster special case
    f  = v;            // f  =  func(t)
    x  = v + v_stride; // x  =  numerator'(t)/w
    wt = -x[dim];      // wt = -denominator'(t)/w
    j  = dim; while (j--) *x++ +=  wt* *f++;
    if (der_count > 1) {
      // 2nd derivative - faster special case
      f = v + v_stride;
      x = f + v_stride;
      // v = func(t), f = func'(t), x = numerator''(t)/w,
      // * wt = -2*denominator'(t)/w, w2 = denominator''(t)/w
      wt *= 2.0;
      w2  = -x[dim];
      j   = dim; while (j--) *x++ += w2* *v++ + wt* *f++;
      if (der_count > 2) {
        df = v_stride-dim;
        // higher derivatives use slower loop
        v -= dim;
        x  = v + v_stride*2;
        for (n = 3; n <= der_count; n++) {
          // computing n-th derivative
          f  = v;
          x += v_stride; // x = numerator^(n)/weight
          w  = v + n*v_stride + dim;
          for (i = 0; i < n; i++) {
            // f = value of i-th derivative
            // w = ((n-i)-th derivative of denominator)/weight
            wt = -BinomialCoefficient(n-i,i) * *w;
            w -= v_stride;
            j  = dim; while (j--) *x++ += *f++ * wt;
            x -= dim;
            f += df;
          }
        }
      }
    }
  }
#endif
}


/*
Description:
  Use the quotient rule to compute partial derivatives of a two parameter
  rational function F(s,t) = X(s,t)/W(s,t), where W is a scalar
  and F and X are vectors of dimension dim.
Parameters:
  dim       - [in]
  der_count - [in] number of derivative (>=0)
  v_stride  - [in] (>= dim+1)
  v         - [in/out]
    v[] is an array of length (der_count+2)*(der_count+1)*v_stride.
    The input array contains derivatives of the numerator and denominator
		functions in the order X, W, Xs, Ws, Xt, Wt, Xss, Wss, Xst, Wst, Xtt, Wtt, ...
    In general, the (i,j)-th derivatives are in the (dim+1) entries of v[]
		v[k], ..., answer[k+dim], where	k = ((i+j)*(i+j+1)/2 + j)*v_stride.
    In the output v[] array the derivatives of X are replaced with
    the derivatives of F and the derivatives of W are divided by
    w = v[dim].
*/

__HOST_AND_DEVICE__ TEMPLATE
void EG_EvaluateQuotientRule2(int dim, int der_count, int v_stride, DOUBLE *v)
{
  DOUBLE F, Fs, Ft, ws, wt, wss, wtt, wst, *f, *x;
  int    i, j, n, q, ii, jj, Fn;

  // comment notation:
  //  X = value of numerator
  //  W = value of denominator
  //  F = X/W
  //  Xs = partial w.r.t. 1rst parameter
  //  Xt = partial w.r.t. 2nd parameter
  //  ...
  //

  // divide everything by the weight
  F = v[dim];
  F = 1.0/F;
  if (v_stride > dim+1) {
    i = ((der_count+1)*(der_count+2)>>1);
    x = v;
    j = dim+1;
    q = v_stride-j;
    while (i--) {
      jj = j;
      while (jj--) *x++ *= F;
      x += q;
    }
  } else {
    i = (((der_count+1)*(der_count+2))>>1)*v_stride;
    x = v;
    while (i--) *x++ *= F;
  }

  if (der_count) {
                              // first derivatives
    f  = v;                   // f = F
    x  = v + v_stride;        // x = Xs/w, x[v_stride] = Xt/w
    ws = -x[dim];             // ws = -Ws/w
    wt = -x[dim+v_stride];    // wt = -Wt/w
    j  = dim;
    while (j--) {
      F            = *f++;
      *x          += ws*F;
      x[v_stride] += wt*F;
      x++;
    }

    if (der_count > 1) {
                              // 2nd derivatives
      f  += (v_stride-dim);   // f = Fs, f[cvdim] = Ft
      x   = v + 3*v_stride;   // x = Xss, x[v_stride] = Xst, x[2*v_stride] = Xtt
      wss = -x[dim];          // wss = -wss/W
      wst = -x[v_stride+dim]; // wst = -Wst/W
      n   = 2*v_stride;
      wtt = -x[n+dim];        // wtt = -Wtt/w
      j   = dim;
      while (j--) {
        F            = *v++;
        Ft           = f[v_stride];
        Fs           = *f++;
        *x          += wss*F + 2.0*ws*Fs;     // Dss
        x[v_stride] += wst*F + wt*Fs + ws*Ft; // Dst
        x[n]        += wtt*F + 2.0*wt*Ft;     // Dtt
        x++;
      }

      if (der_count > 2) {
        // general loop for higher derivatives
        v -= dim;               // restore v pointer to input value
        f  = v + 6*v_stride;    // f = Xsss
        for (n = 3; n <= der_count; n++) {
          for (j = 0; j <= n; j++ ) {
            // f = Ds^i Dt^j
            i = n-j;
            for (ii = 0; ii <= i; ii++) {
              ws = BinomialCoefficient(ii,i-ii);
              for (jj = ii?0:1; jj <= j; jj++) {
                q  = ii+jj;
                Fn = ((q*(q+1))/2 + jj)*v_stride+dim;
                // wt = -(i choose ii)*(j choose jj)*W(ii,jj)
                wt = -ws*BinomialCoefficient(jj,j-jj)*v[Fn];
                q  = n-q;
                Fn = ((q*(q+1))/2 + j-jj)*v_stride; // X(i-ii,j-jj) = v[Fn]
                for (q = 0; q < dim; q++) f[q] += wt*v[Fn+q];
              }
            }
            f += v_stride;
          }
        }
      }
    }
  }

}


__HOST_AND_DEVICE__ TEMPLATE
static void
EG_IncreaseBezierDegree(int dim, int order, int cv_stride, DOUBLE* cv)

{
  double a0, a1, d, c0, c1;
  int    j;
  DOUBLE *newcv   = cv;
  const int cvdim = dim+1;
  const int dcv   = cv_stride - cvdim;

  j = cv_stride*order;
  newcv += j;
#ifdef LITE
  memcpy(newcv, newcv-cv_stride, cvdim*sizeof(*newcv));
#else
  for (int ii = 0; ii < cvdim; ii++) newcv[ii] = (newcv-cv_stride)[ii];
#endif
  newcv -= (dcv+1);
  cv = newcv - cv_stride;
  a0 = order;
  a1 = 0.0;
  d  = 1.0/a0;
  while (--order) {
    a0 -= 1.0;
    a1 += 1.0;
    c0  = d*a0;
    c1  = d*a1;
    j   = cvdim;
    while (j--) {
      *newcv = c0 * *cv + c1 * *newcv;
      cv--;
      newcv--;
    }
    cv    -= dcv;
    newcv -= dcv;
  }
}


__HOST_AND_DEVICE__ TEMPLATE
static bool
EG_RemoveBezierSingAt0(int dim, int order, int cv_stride, DOUBLE* cv)
{
  const int cvdim = dim+1;
  int j,k,ord0;

  ord0 = order;
  while(cv[dim] == 0.0) {
    order--;
    if (order < 2)
      return false;
    j = dim;
    while(j--) {
      if (cv[j] != 0.0)
        return false;
    }
    for (j=0;  j<order;  j++) {
      for (k=0;  k<cvdim;  k++)
        cv[j*cv_stride+k] = (order*cv[(j+1)*cv_stride+k])/(j+1);
    }
  }
  while (order < ord0)
    EG_IncreaseBezierDegree(dim, order++, cv_stride, cv);
  return true;
}


/*****************************************************************************
Evaluate a Rational Bezier

INPUT:
  dim
    (>= 1) dimension of Bezier's range
  order
    (>= 2) (order = degree+1)
  cv
    array of (dim+1)*order doubles that define the Bezier's control vertices.
  t0, t1  (t0 != t1)
    Bezier's domain.   Mathematically, Beziers have domain [0,1].  In practice
    Beziers are frequently evaluated at (t-t0)/(t1-t0) and the chain
    rule is used to evaluate the derivatives.  This function is the most
    efficient place to apply the chain rule.
  t
    Evaluation parameter
  der_count
    (>= 0)  number of derivatives to evaluate
  answer
     answer[i] is NULL or points to an array of dim doubles.
OUTPUT:
    0: successful
   -1: failure - rational function had nonzero numerator and zero
                 denominator
  answer
     bez(t)   = answer[0]
     bez'(t)  = answer[1]
              ...
        (n)
     bez  (t) = answer[n]
COMMENTS:
  Use de Casteljau's algorithm.  Rational fuctions with removable singularities
  (like x^2/x) are efficiently and correctly handled.
*****************************************************************************/

__HOST_AND_DEVICE__ TEMPLATE
bool EG_Bezier1DRat(int dim,              // dimension
                    int order,            // order
                    int cv_stride,        // cv_stride >= dim+1
                    const DOUBLE *cv,     // cv[order*cv_stride] array
                    double t0, double t1, // domain
                    int der_count,        // number of derivatives to compute
                    DOUBLE t,             // evaluation parameter
                    int v_stride,         // v_stride (>=dimension)
                    DOUBLE *v           ) // v[(der_count+1)*v_stride] array
{
  unsigned char stack_buffer[4*64*sizeof(DOUBLE)];
  DOUBLE delta_t;
  DOUBLE alpha0;
  DOUBLE alpha1;
  DOUBLE *cv0, *cv1;
  int i, j, k;
  DOUBLE *CV, *tmp;
  void* free_me = 0;
  const int degree = order-1;
  const int cvdim  = dim+1;
  size_t sizeofCV;

  if (cv_stride < cvdim) cv_stride = cvdim;

#ifdef LITE
  memset(v, 0, v_stride*(der_count+1)*sizeof(*v));
#else
  for (i = 0; i < v_stride*(der_count+1); i++) v[i] = 0.;
#endif

  i = order*cvdim;
  j = 0;
  if (der_count > degree) j = (der_count-degree)*cvdim;

  sizeofCV = (i+j)*sizeof(*CV);
  CV = (DOUBLE *)( (sizeofCV <= sizeof(stack_buffer)) ? stack_buffer : (free_me=malloc(sizeofCV)) );
#ifdef LITE
  if (j) memset(CV+i, 0, j*sizeof(*CV));
#else
  if (j) for (int ii = 0; ii < j; ii++) CV[ii+i] = 0.;
#endif
  cv0 = CV;
  if (t0 == t || (t <= 0.5*(t0+t1) && t != t1)) {
    for (i = 0; i < order; i++) {
#ifdef LITE
      memcpy(cv0, cv, cvdim*sizeof(*cv0));
#else
      for (int ii = 0; ii < cvdim; ii++) cv0[ii] = cv[ii];
#endif
      cv0 += cvdim;
      cv  += cv_stride;
    }
    cv     -= (cv_stride*order);
    delta_t = 1.0/(t1 - t);
    alpha1  = 1.0/(t1-t0);
    alpha0  = (t1-t)*alpha1;
    alpha1 *= t-t0;
  } else {
    cv += (cv_stride*order);
    k   = order;
    while( k--) {
      cv  -= cv_stride;
#ifdef LITE
      memcpy( cv0, cv, cvdim*sizeof(*cv0) );
#else
      for (int ii = 0; ii < cvdim; ii++) cv0[ii] = cv[ii];
#endif
      cv0 += cvdim;
    }
    delta_t = 1.0/(t0 - t);
    alpha0  = 1.0/(t1-t0);
    alpha1  = (t1-t)*alpha0;
    alpha0 *= t-t0;
  }

  /* deCasteljau (from the right) */
  if (alpha1 != 0.0) {
    j = order; while (--j) {
      cv0 = CV;
      cv1 = cv0 + cvdim;
      i   = j;
      while (i--) {
        k = cvdim;
        while (k--) {
#ifndef __clang_analyzer__
          *cv0 = *cv0 * alpha0 + *cv1 * alpha1;
#else
          *cv0 =        alpha0 +        alpha1;
#endif
          cv0++;
          cv1++;
        }
      }
    }
  }

  /* check for removable singularity */
#ifndef __clang_analyzer__
  if (CV[dim] == 0.0) {
    if ( !EG_RemoveBezierSingAt0(dim,order,cvdim,CV) ) {
      if (free_me) free(free_me);
      return false;
    }
  }
#endif

  /* Lee (from the right) */
  if (der_count) {
    tmp=CV;
    alpha0 = order;
    j   = (der_count>=order) ? order : der_count+1;
    CV += cvdim*j; while(--j) {
      alpha0 -= 1.0; cv1 = CV; cv0 = cv1-cvdim;
      i       = j;
      while(i--) {
        alpha1 = alpha0 * delta_t;
        k=cvdim; while(k--) {
          cv0--;
          cv1--;
#ifndef __clang_analyzer__
          *cv1 = alpha1*(*cv1 - *cv0);
#else
          *cv1 = alpha1;
#endif
        }
      }
    }
    CV=tmp;
  }

  if (order == 2) {
    j = cv_stride;
    for (i = 0; i < cvdim; i++, j++)
      if (cv[i] == cv[j]) CV[i] = cv[i];
  }

  EG_EvaluateQuotientRule(dim, der_count, cvdim, CV);

  for (i = 0; i <= der_count; i++) {
#ifdef LITE
    memcpy(v, CV, dim*sizeof(*v));
#else
#ifndef __clang_analyzer__
    for (int ii = 0; ii < dim; ii++) v[ii] = CV[ii];
#endif
#endif
    v  += v_stride;
    CV += cvdim;
  }

  if (free_me) free(free_me);

  return true;
}

#ifndef LITE
template void EG_EvaluateQuotientRule( int dim, int der_count, int v_stride,
                                       SurrealS<1> *v );
template void EG_EvaluateQuotientRule( int dim, int der_count, int v_stride,
                                       double *v );

template void EG_EvaluateQuotientRule2( int dim, int der_count, int v_stride,
                                        SurrealS<1> *v );
template void EG_EvaluateQuotientRule2( int dim, int der_count, int v_stride,
                                        double *v );

template bool EG_Bezier1DRat( int dim, int order, int cv_stride,
                              const SurrealS<1> *cv, double t0, double t1,
                              int der_count, SurrealS<1> t, int v_stride,
                              SurrealS<1> *v );
template bool EG_Bezier1DRat( int dim, int order, int cv_stride,
                              const double *cv, double t0, double t1,
                              int der_count, double t, int v_stride,
                              double *v );
#endif

