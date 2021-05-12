#include "egads.h"
#include "regQuads.h"

//#define DEBUG
//#define DEBUG2
//#define REPORT

#ifdef DEBUG
 #define REPORT
 #define DEBUG2
#endif

#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ extern "C" __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif


#ifdef _MSC_VER
#if _MSC_VER <= 1700
double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	  double t, y;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
#endif
#endif

#ifdef REPORT
  #include <time.h>
#endif


__HOST_AND_DEVICE__ extern void
EG_getSidepoint(const ego face, double fac, const double *uvm,
                const double *uvp, /*@null@*/ const double *uvl,
                /*@null@*/ const double *uvr, double *uv);


/* BASE-LEVEL FUNCTION */
__HOST_AND_DEVICE__ static int inList(int n, int *list, int p)
{
  int i;

  if (n == 0) return -1;
  for (i = 0; i < n; i++)
    if (list[i] == p) return i;
  return -1;
}


/* IO FUNCTIONS */
__HOST_AND_DEVICE__ static void meshCount(meshMap *qm, int *nI, int *nV, int *nQ)
{
  int i, qSum, vSum = 0, vSum2 = 0;

  for (i = qSum = 0; i < qm->totQ; i++)
    if (qm->qIdx[4 * i] != -2) qSum++;
  for (i = vSum = 0; i < qm->totV; i++) {
      if (qm->vType[i] == -2) continue;
      vSum2++;
      if (qm->valence[i][1] != 4) vSum++;
  }
  *nQ = qSum;
  *nI = vSum;
  *nV = vSum2;
}


#ifdef REPORT
static void meshStats(meshMap *qm, int sweep)
{
  int  i,len, mv, *valptr = NULL;
  char buffer[200];
  FILE *fout = NULL;

  if (sweep == 0)
    snprintf(buffer,200, "OriMeshStats_face_%d_.txt", qm->fID);
  else if (sweep == 10)
    snprintf(buffer,200, "PreMeshStats_face_%d_.txt", qm->fID);
  else
    snprintf(buffer,200, "RegMeshStats_face_%d_.txt", qm->fID);
  fout = fopen(buffer, "w");
  if (fout == NULL) {
      printf(" Error writing in file %s\n", buffer);
      return;
  }
  len = qm->totV;
  for (mv = i = 0; i < len; i++) {
    if (qm->vType[i] == -2) continue;
    mv = MAX(mv, qm->valence[i][1]);
  }
  valptr = (int*)EG_alloc(2 * (mv + 1) * sizeof(int));
  if (valptr == NULL) {
    fclose(fout);
    printf(" meshStats VALENCE POINTER IS NULL!!\n");
    return;
  }
  for (i = 0 ; i < mv; i++) {
      valptr[mv + i] = 0;valptr[i] = 0;
  }
  for (i = 0; i < len; ++i) {
      if (qm->vType[i] == -2) continue;
      if (qm->vType[i] == -1) ++valptr[     qm->valence[i][1]];
      else                    ++valptr[mv + qm->valence[i][1]];
  }
  printf(" Writing in %s \n", buffer);
  fprintf(fout,"---- TOTAL VERTICES %d TOTAL QUADS %d --------------\n",
          qm->totV - qm->remV[0], qm->totQ - qm->remQ[0]);
  fprintf(fout," INTERIOR VERTICES\n");
  for (i = 0; i < mv; ++i) {
      if (valptr[i]  > 0)
        fprintf(fout," VALENCE %d = %d VERTICES\n", i, valptr[i]);
  }
  fprintf(fout," BOUNDARY VERTICES\n");
  for (i = 0; i < mv; ++i) {
      if (valptr[2 * i]  > 0)
        fprintf(fout," VALENCE %d = %d VERTICES\n", i, valptr[2 * i]);
  }
  meshCount(qm, &i, &len, &mv);
  fprintf(fout,"\n\nMesh has %d QUADS, %d / %d irregular vertices (%.2f %% ) ============\n",
          mv, i, len, (double) i * 100.0 / (double) len);
  EG_free(valptr);
  fclose(fout);
}


static void gnuData(meshMap *qm, /*@null@*/ char *name)
{
  int    i, j, k, v, d;
  double eval[18*2], average[2], dist;
  char   name2[100];
  FILE   *fout = NULL;

  if (qm->fID == 0) return;
  if (name == NULL) {
      snprintf(name2, 100, "FACE_%d_M_%d",qm->fID, qm->plotcount++);
      printf(" Writing in File %s  \n", name2);
      fout = fopen(name2, "w");
  } else {
      printf(" Writing in File %s  \n", name);
      fout = fopen(name, "w");
  }
  if (fout == NULL) return;

  for (i = 0; i < qm->totQ; ++i) {
      if (qm->qIdx[4 * i]                == -2 ||
          qm->vType[qm->qIdx[4 * i] - 1] == -2) continue;
      for (k = 0; k < 4; ++k) {
          v  = qm->qIdx[4*i + k] - 1;
          fprintf(fout, "%lf %lf %lf %d %lf %lf \n", qm->xyzs[3*v  ],
                  qm->xyzs[3*v+1], qm->xyzs[3*v+2], v + 1,
                  qm->uvs [2*v  ], qm->uvs [2*v+1]);
          if (qm->vType[v] != -1) continue;
          dist = 0.0;
          EG_evaluate(qm->face, &qm->uvs[2*v], eval);
          for (d = 0; d < 3; ++d)
            dist += (eval[d] - qm->xyzs[3*v + d]) *
            (eval[d] - qm->xyzs[3*v + d]);
          dist = sqrt(dist);
          if (dist > EPS08) {
              printf(" DIST = %11.2e  IN QUAD %d  VERTEX %d. UVs and xyzs are mismatched.  UV %lf  %lf \n",
                     dist,i+1, v+1, qm->uvs[2*v], qm->uvs[2*v + 1]);
              for (d = 0; d < 3; ++d)
                printf(" %lf  != %lf \t", eval[0], qm->xyzs[3*v + d]);
          }
      }
      v = qm->qIdx[4*i] - 1;
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n",  qm->xyzs[3*v  ],
              qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1,
              qm->uvs [2*v    ], qm->uvs [2*v + 1]);
      fprintf(fout,"\n\n");
      average[0] = 0.0; average[1] = 0.0;
      eval[18] = eval[19] = eval[20] = 0.0;
      for (j = 0; j < 4; j++) {
          k = qm->qIdx[4 * i + j] - 1;
          average[0] += 0.25 * qm->uvs[2 * k    ];
          average[1] += 0.25 * qm->uvs[2 * k + 1];
          eval[18]   += 0.25 * qm->xyzs[3 * k];
          eval[19]   += 0.25 * qm->xyzs[3 * k +1];
          eval[20]   += 0.25 * qm->xyzs[3 * k +2];
      }
      EG_invEvaluateGuess(qm->face, &eval[18], average, eval);
      fprintf(fout, "%lf %lf %lf %d %lf %lf\n", eval[0], eval[1], eval[2],
              i + 1, average[0], average[1]);
      fprintf(fout,"\n\n");
  }
  fclose(fout);
  if (name) {
      snprintf(name2, 100, "Edges_%s", name);
      printf(" Writing in File %s  \n", name2);
      fout = fopen(name2, "w");
      if (fout == NULL ) return;
      for (i = 0; i < qm->totQ; i++) {
          if (qm->qIdx[4 * i] == -2) continue;
          for (j = 0; j < 4; j++) {
              if (qm->qAdj[4 * i + j] == -1 ) {
                  fprintf(fout, "%lf %lf %lf %d\n",
                          qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1)    ],
                          qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 1],
                          qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 2], qm->qIdx[4 * i + j]);
                  fprintf(fout, "%lf %lf %lf %d\n\n\n",
                          qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1)    ],
                          qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 1],
                          qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 2], qm->qIdx[4 * i + (j + 1)%4]);
              }
          }
      }
      fclose(fout);
  }
}


static void wvsData(meshMap *qm, /*@null@*/ char *buffer)
{
  int  i, j, k, n;
  FILE *fil = NULL;
  char aux[100];

  if (qm->fID == 0 ) return;
  if (buffer == NULL) {
      snprintf(aux, 100, "wvsFace_%d_M_%d", qm->fID, qm->plotcount++);
      fil = fopen(aux, "w");
  } else {
      fil = fopen (buffer, "w");
      printf(" Writing in File %s  \n", buffer);
  }
  if (fil == NULL) {
      printf(" I couldn't write in wvsData file !!!!\n ");
      return ;
  }
  fprintf(fil, "%d %d\n", qm->totV, qm->totQ);
  for (i = 0; i < qm->totV; i++) {
      fprintf(fil, "%lf %lf %lf \n", qm->xyzs[3 * i    ],
              qm->xyzs[3 * i + 1], qm->xyzs[3 * i + 2]);
  }
  fprintf(fil,"\n");
  for (n = i = 0; i < qm->totQ; i++) {
      fprintf(fil, "%d %d %d %d\n",  qm->qIdx[4 * i    ],
              qm->qIdx[4 * i + 1], qm->qIdx[4 * i + 2],
              qm->qIdx[4 * i + 3]);
      if (qm->qAdj[4 * i    ] == -1 ) n++;
      if (qm->qAdj[4 * i + 1] == -1 ) n++;
      if (qm->qAdj[4 * i + 2] == -1 ) n++;
      if (qm->qAdj[4 * i + 3] == -1 ) n++;
  }
  fprintf(fil,"%d\n", n );
  for (k = i = 0; i < qm->totQ; i++) {
      for (j = 0; j < 4; j++) {
          if (qm->qAdj[4 * i + j] == -1 ) {
              fprintf(fil, "%lf %lf %lf\n",
                      qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1)    ],
                      qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 1],
                      qm->xyzs[3 * (qm->qIdx[4 * i + j]- 1) + 2]);
              fprintf(fil, "%lf %lf %lf\n",
                      qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1)    ],
                      qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 1],
                      qm->xyzs[3 * (qm->qIdx[4 * i + (j + 1)%4]- 1) + 2]);
              k++;
          }
      }
  }
  if (k != n)
    printf(" wvsData Error filling the edges segments!! %d != %d\n", k, n);
  fclose(fil);
  return ;
}

#endif


/* DEBUG FUNCTIONS */
#ifdef DEBUG

__HOST_AND_DEVICE__ static void printQuad (meshMap *qm, int id)
{
  int    i, v = 0 , val = 0;
  double uv[3], uva[2], xyz[18];

  uv[0] = uv[1] = uv[2] = 0.0;
  --id;
  printf("---------------------------------------\n");
  printf(" QUAD %d HAS VERTICES ", id+1);
  for (i = 0; i < 4; ++i) {
      v   = qm->qIdx[4 * id + i];
      if (v <= 0) printf(" DEL %d %d \n ", v, 0);
      else {
          val = qm->valence[v - 1][1];
          if (qm->vType[v -1] > 0) printf(" %d ( val = %d Type %d BD ANG %lf)  ",
                                          v, val, qm->vType[v-1], qm->bdAng[v-1]);
          else printf(" %d ( val = %d Type %d)  ", v, val, qm->vType[v-1]);
      }
  }
  printf("\t AND ADJACENT QUADS ");
  for (i = 0; i < 4; ++i) printf(" %d ", qm->qAdj[4*id + i]);
  printf("\n\n Coordinates \n");
  for (i = 0; i <= 4; i++) {
      v = qm->qIdx[4 * id + i%4] - 1;
      printf("%lf %lf %lf %d %lf %lf \n",
             qm->xyzs[3 * v    ], qm->xyzs[3 * v + 1],
             qm->xyzs[3 * v + 2], v + 1,
             qm->uvs [2 * v + 1], qm->uvs [2 * v    ]);
      if (i < 4) {
          uv[0] += 0.25 * qm->xyzs[3 * v    ];
          uv[1] += 0.25 * qm->xyzs[3 * v + 1];
          uv[2] += 0.25 * qm->xyzs[3 * v + 2];
      }
  }
  i = EG_invEvaluate(qm->face, uv, uva, xyz);
  printf("\n\n\n %d  %lf %lf %lf %d %lf %lf \n",
         i, xyz[0], xyz[1], xyz[2], id + 1, uv[0], uv[1]);
  printf(" \n\n---------------------------------------\n");
}


__HOST_AND_DEVICE__ static void printQuadGroup(meshMap *qm, quadGroup qg)
{
  int i;

  for (i = 0; i < 2; i++) {
      if (qg.q[i] < 0 || qg.q[i] > qm->totQ ) {
          printf(" Quad %d is out of bounds !!! \n", qg.q[i] );
          return;
      }
      printQuad(qm, qg.q[i]);
  }
  for (i = 0; i < 6; i++) {
      printf(" QV (%d) = %d ", i, qg.verts[i] );
      if (qg.verts[i] < 0) printf(" NO VAL !! NEGATIVE! \n");
      if (qg.verts[i] ==0 || qg.verts[i] > qm->totV) {
          printf(" QUAD HAS VERTEX %d OUT OF BOUNDS!!\n ", qg.verts[i] );
      }
      else printf("type %d valence %d \n", qm->vType[qg.verts[i] - 1],
                  qg.vals[i]);
  }
}


__HOST_AND_DEVICE__ static int checkMesh(meshMap *qm)
{
  int i, j, k, quad, val, v1, v2;

  for (i = 0; i < qm->totV; i++) {
      if (qm->vType[i] == -2) continue;
      val = qm->valence[i][2];
      if (val <= 0 || val > qm->totV) {
          printf(" checkMesh vertex %d valence out of bounds %d (max %d ) !!\n",
                 i + 1, qm->valence[i][2], qm->totV);
          return EGADS_INDEXERR;
      }
      quad = qm->valence[i][0];
      if (quad <= 0 || quad > qm->totQ) {
          printf("In checkMesh quad for %d is out of bounds!!! %d > %d \n ",
                 i + 1, quad, qm->totQ);
          return EGADS_INDEXERR;
      }
      for (j = 0; j < val; ++j) {
          v2 = qm->valence[i][3 + j] - 1;
          for (k = 0; k < qm->valence[v2][2]; k++)
            if (qm->valence[v2][3 + k] == i + 1) {
                k = -1;
                break;
            }
          if (k != -1) {
              printf(" checkMesh Vertex %d has assigned %d as link but %d doesn't point at %d\n",
                     i + 1, v2 + 1, v2 + 1,  i + 1);
              for (k = 0; k < qm->valence[i][2]; k++)
                printf(" V1 %d ==> %d \n ", k, qm->valence[i][3+k]);
              printf("\n\n");
              for (k = 0; k < qm->valence[v2][2]; k++)
                printf(" V1 %d ==> %d \n ", k, qm->valence[v2][3+k]);
              return EGADS_INDEXERR;
          }
      }
  }
  for (i = 0; i < qm->totQ; i++) {
      if (qm->qIdx[4 * i] == -2) continue;
      for (j = 0; j < 4; ++j) {
          v1 = qm->qAdj[4 * i  + j] - 1;
          if (v1 < 0) continue;
          val = -1;
          for (k = 0; k < 4; ++k) {
              if (qm->qAdj[4 * v1 + k] == i + 1) {
                  val = 1;
                  break;
              }
          }
          if (val == -1) {
              printf (" checkMesh quads %d and %d don't point at each other\n",
                       i +1, v1 + 1);
              printQuad(qm, i + 1);
              printQuad(qm,v1 + 1);
              return EGADS_INDEXERR;
          }
      }
  }
  return EGADS_SUCCESS;
}
#endif


__HOST_AND_DEVICE__ static void updateVertex(meshMap *qm, int vID, double *uv)
{
  double eval[18];

  if (EG_evaluate(qm->face, uv, eval) != EGADS_SUCCESS ||
      qm->vType[vID -1]   != -1) return;
  qm->uvs [2 * (vID - 1)    ] = uv  [0];
  qm->uvs [2 * (vID - 1) + 1] = uv  [1];
  qm->xyzs[3 * (vID - 1)    ] = eval[0];
  qm->xyzs[3 * (vID - 1) + 1] = eval[1];
  qm->xyzs[3 * (vID - 1) + 2] = eval[2];
}


__HOST_AND_DEVICE__ static int EG_regBoundQuad(meshMap *qm, int qID)
{
  int iA, iB, iC, iD, vA, vB, vC, vD;

  if (qm->regBd == 0) return -1;
  iA = qm->qIdx [4 * (qID -1)   ] -1;
  iB = qm->qIdx [4 * (qID -1) +1] -1;
  iC = qm->qIdx [4 * (qID -1) +2] -1;
  iD = qm->qIdx [4 * (qID -1) +3] -1;
  if ((qm->vType[iA] == -1 && qm->vType[iB] == -1 &&
      qm->vType[iC] == -1 && qm->vType[iD] == -1) ||
      qm->vType[iA] * qm->valence[iA][2]   == -2 ||
      qm->vType[iB] * qm->valence[iB][2]   == -2 ||
      qm->vType[iC] * qm->valence[iC][2]   == -2 ||
      qm->vType[iD] * qm->valence[iD][2]   == -2)
    return -1;
  vA = vB = vC = vD = 4;
  if (qm->vType[iA] > 0) vA = qm->valence[iA][1];
  if (qm->vType[iB] > 0) vB = qm->valence[iB][1];
  if (qm->vType[iC] > 0) vC = qm->valence[iC][1];
  if (qm->vType[iD] > 0) vD = qm->valence[iD][1];
  if (vA * vB * vC * vD == 4 * 4 * 4 * 4) return 1;
  else return 0;
}


__HOST_AND_DEVICE__ static int
EG_adjQtoPair(meshMap *qm, int qID, int v1, int v2, int *adj)
{
  int i, aux = -1;

  adj[0] = -1; adj[1] = -1;
  for (i = 0; i < 4; ++i) {
      if (qm->qIdx[4*(qID - 1) + i] == v1) adj[0] = i;
      if (qm->qIdx[4*(qID - 1) + i] == v2) aux    = i;
      if (aux != -1 && adj[0] != -1) break;
  }
  if (aux == -1 || adj[0] == -1) return EGADS_SUCCESS;
  if      (abs(adj[0] - aux) == 3) adj[0] = 3;
  else if (aux < adj[0]          ) adj[0] = aux;
  adj[1] = qm->qAdj[4*(qID - 1) + adj[0]];
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int EG_quadVertIdx(meshMap *qm, int q, int v)
{
  int i = 0;

  for (i = 0; i < 4; i++)
    if (qm->qIdx[4 * (q - 1) + i] == v) return i;
  return -1;
}


__HOST_AND_DEVICE__ static void EG_freeStar(vStar **star)
{
  if ((*star) == NULL) return;
  EG_free((*star)->verts);
  EG_free((*star)->quads);
  EG_free((*star)->idxQ);
  EG_free((*star)->idxV);
  EG_free((*star)->area);
  EG_free((*star)->angle);
  EG_free((*star)->ratio);
  EG_free (*star);
  *star = NULL;
}


__HOST_AND_DEVICE__ static int EG_buildStar(meshMap *qm, vStar **star, int vID)
{
  int i = 0, id0 = -1, q = 0, auxV, auxQ, v = 0, quadID, prevQuad,
      it = 0, it2 = 0, adj[2], *vertex = NULL, *quads = NULL;
  int qLoop[8] = {0, 1, 2, 3, 0, 1, 2, 3};

  vertex = (int *) EG_alloc(qm->totV * sizeof(int));
  quads  = (int *) EG_alloc(qm->totQ * sizeof(int));
  if (vertex == NULL || quads == NULL) {
      printf(" EG_buildStar MALLOC at quads & verts!!\n ");
      if (vertex != NULL) EG_free(vertex);
      if (quads  != NULL) EG_free(quads);
      return EGADS_MALLOC;
  }
  // quads are -1 bias
  quadID      = qm->valence[vID - 1][0] - 1;
  vertex[v++] = vID;
  it          = 0;
  do {
      id0 = EG_quadVertIdx(qm, quadID + 1, vID);
      if (id0 < 0) {
          printf(" In EG_buildStar id for Vert %d in Quad %d is %d !!\n",
                 vID, quadID + 1, id0);
#ifdef DEBUG
          printQuad(qm, quadID + 1);
#endif
          EG_free(vertex);
          EG_free(quads);
          return EGADS_INDEXERR;
      }
      for (i = 1; i <= 2; ++i)
        vertex[v++] = qm->qIdx[4 * quadID + qLoop[id0 + i]];
      quads[q++] = quadID + 1;
      prevQuad   = quadID;
      quadID     = qm->qAdj[4 * prevQuad + qLoop[id0 + 3 ]] - 1;
      if (quadID  < 0) { //make a "ghost" quad
          auxQ        = prevQuad;
          vertex[v++] = qm->qIdx[4*auxQ + qLoop[id0 + 3]];
          auxV        = qm->qIdx[4*auxQ + qLoop[id0 + 1]];
          it2         = 0;
          do {
              i    = EG_adjQtoPair(qm, auxQ + 1, vID, auxV, adj);
              if (adj[1] == -1 || i != EGADS_SUCCESS) break;
              auxQ = adj[1] - 1;
              i    = EG_quadVertIdx(qm, auxQ + 1, vID);
              if (i < 0) {
                  printf(" In buildStar vertex Id %d in quad %d is %d\n",
                         vID, quadID + 1, i);
#ifdef DEBUG
                  printQuad(qm, quadID + 1);
#endif
                  EG_free(vertex);
                  EG_free(quads);
                  return EGADS_INDEXERR;
              }
              auxV = qm->qIdx[4 * auxQ + qLoop[i + 1]];
              it2++;
              if (it2 > 200) {
                  printf(" stuck in interior loop of build star!!!!!!!!\n");
                  EG_free(quads);
                  EG_free(vertex);
                  return EGADS_RANGERR;
              }
          } while (adj[1] != - 1);
          quads[q++]  = -1;
          vertex[v++] = -1;
          quadID      = auxQ;
      }
      if (quadID < 0) {
          printf(" I am stuck in build star. Pointing a NULL quad \n");
          EG_free(vertex);
          EG_free(quads);
          return EGADS_INDEXERR;
      }
      it++;
      if (it > 200) {
          printf(" EG_buildStar:: stuck in outer loop of build star!!!!!!!!\n");
          EG_free(vertex);
          EG_free(quads);
          return EGADS_RANGERR;
      }
  } while (quadID + 1 != quads[0]);

  if (*star != NULL) EG_freeStar(&(*star));
  *star = (vStar *) EG_alloc(sizeof(vStar));
  if ((*star) == NULL) {
      EG_free(vertex);
      EG_free(quads);
      return EGADS_MALLOC;
  }
  (*star)->nQ    = q;
  (*star)->nV    = v;
  (*star)->type  = -1;
  (*star)->verts = (int    *) EG_alloc(    v * sizeof(int));
  (*star)->quads = (int    *) EG_alloc(    q * sizeof(int));
  (*star)->idxV  = (int    *) EG_alloc(2 * v * sizeof(int));
  (*star)->idxQ  = (int    *) EG_alloc(2 * q * sizeof(int));
  (*star)->area  = (int    *) EG_alloc(    q * sizeof(int));
  (*star)->ratio = (double *) EG_alloc(    q * sizeof(double));
  (*star)->angle = (double *) EG_alloc(    q * sizeof(double));
  if ((*star)->verts == NULL || (*star)->quads == NULL ||
      (*star)->idxV  == NULL || (*star)->idxQ  == NULL ||
      (*star)->area  == NULL || (*star)->ratio == NULL ||
      (*star)->angle == NULL ) {
      EG_free((*star));
      EG_free(vertex);
      EG_free( quads);
      return EGADS_MALLOC;
  }
  for (i = 0; i < q; ++i) {
      (*star)->quads[i    ] = quads[i];
      (*star)->idxQ [i    ] = i;
      (*star)->idxQ [q + i] = i;
  }
  for (i = 0; i < v; ++i) {
      (*star)->verts[i    ] = vertex[i];
      (*star)->idxV [i    ] = i;
      (*star)->idxV [v + i] = i + 1;
  }
  EG_free (vertex);
  EG_free (quads);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int EG_setValence(meshMap *qm, int vID)
{
  int   i, i4, ip, in, *iv = NULL, j, si, n;

  if (vID <= 0 || vID > qm->totV) {
        printf(" VERTEX %d is out of bounds !\n ", vID);
        return EGADS_INDEXERR;
  }
  si = vID - 1;
  i  = EG_buildStar(qm, &qm->star[si], vID);
  if (i != EGADS_SUCCESS || qm->star[si] == NULL) {
      printf(" In setValence for vertex %d stat from EG_buildStar = %d\n",
             vID, i);
      return i;
  }
  if (qm->valence[si][2] < qm->star[si]->nQ) {
      qm->valence[si] = (int *) EG_reall(qm->valence[si],
                                         (3 + qm->star[si]->nQ) * sizeof(int));
      if (qm->valence[si] == NULL) {
          printf(" VALENCE is NULL %d \n ", si + 1 );
          return EGADS_MALLOC;
      }
  }
  n  =  qm->star[si]->nQ;
  iv = (int *) EG_alloc(3 * n * sizeof(int));
  if (iv == NULL) {
    printf(" EG_setValence iv pointer is NULL \n");
    return EGADS_MALLOC;
  }
  qm->valence[si][2] = n; // real valence
  qm->valence[si][1] = n;
  if (qm->vType[si] == 2 ||
      qm->vType[si] == 0) qm->valence[si][1] = n + 2; //adjusted based on vertex type, eg, a corner is regular when val = 2
  else if (qm->vType[si] == 3) qm->valence[si][1] = n + 1;
  else if (qm->vType[si] >= 5) qm->valence[si][1] = n - (qm->vType[si] - 4);
  for (in = ip = i4 = i = 0; i < n; i++) {
      if (qm->star[si]->verts[2 * i + 1] == -1) continue;
      if      (qm->vType[qm->star[si]->verts[2 * i + 1] - 1] >= 4)
        iv[2 * n + (i4++)] = qm->star[si]->verts[2 * i + 1];
      else if (qm->vType[qm->star[si]->verts[2 * i + 1] - 1] >= 2)
        iv[    n + (ip++)] = qm->star[si]->verts[2 * i + 1];
      else
        iv[          in++] = qm->star[si]->verts[2 * i + 1];
  }
  for (j = i = 0; i < i4; i++)
    qm->valence[si][3 + (j++)] = iv[2 * n + i];
  for (i = 0; i < ip; i++)
    qm->valence[si][3 + (j++)] = iv[    n + i];
  for (i = 0; i < in; i++)
    qm->valence[si][3 + (j++)] = iv[i];
  EG_free(iv);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_backupQuads(meshMap *qm, int *nq, int *qlist, Quad **quad)
{
  int   i, j, q, v, qcount, *qaux = NULL, k, add = 0;

  qaux      = (int *) EG_alloc(qm->totQ * sizeof(int));
  if (qaux == NULL) return EGADS_MALLOC;
  for (qcount = q = 0; q < *nq; q++) {
      if (qlist[q] == -1) continue;
      for (i = 0; i < 4; i++) {
          v    = qm->qIdx[4 * (qlist[q] - 1) + i] - 1;
          if (qm->star[v] == NULL) {
              EG_free (qaux);
              printf("Star for vertex %d is NULL !!\n ", v+ 1);
              return EGADS_MALLOC;
          }
          for (j = 0; j < qm->star[v]->nQ; j++) {
              if (qm->star[v]->quads[j] == -1) continue;
              for (add = k = 0; k < qcount; k++) {
                  if (qaux[k] == qm->star[v]->quads[j]) {
                      add = 1;
                      break;
                  }
              }
              if (add == 0 || qcount == 0) qaux[qcount++] = qm->star[v]->quads[j];
          }
      }
  }
  (*quad)      = (Quad *) EG_alloc(qcount * sizeof(Quad));
  if ((*quad) == NULL) {
      EG_free(qaux);
      return EGADS_MALLOC;
  }
  for (q = 0; q < qcount; q++) {
      (*quad)[q].id = qaux[q];
      for (j = 0; j < 4; j++) {
          (*quad)[q].qadj [j] = qm->qAdj[4 * (qaux[q] - 1) + j];
          (*quad)[q].verts[j] = qm->qIdx[4 * (qaux[q] - 1) + j];
      }
  }
  EG_free(qaux);
  *nq = qcount;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int EG_restoreQuads(meshMap *qm, Quad *quad, int nq)
{
  int i, j, *vid = NULL, k, stat;

  vid      = (int *) EG_alloc(qm->totV * sizeof(int));
  if (vid == NULL) return EGADS_MALLOC;
  for (k  = i = 0; i < nq; i++) {
      if (quad[i].id == -1) continue;
      for (j = 0; j < 4; j++) {
          qm->qAdj[4 * (quad[i].id - 1) + j] = quad[i].qadj [j];
          qm->qIdx[4 * (quad[i].id - 1) + j] = quad[i].verts[j];
          if (k == qm->totV) continue;
          if (inList(k, vid, quad[i].verts[j]) == -1) {
              qm->valence[quad[i].verts[j] - 1][0] = quad[i].id;
              vid[k++] = quad[i].verts[j];
          }
      }
  }
  for (i   = 0; i < k; i++) {
      stat = EG_setValence(qm, vid[i]);
      if (stat != EGADS_SUCCESS) {
          printf(" EG_restoreQuads failed at setting valence for %d \n",
                 vid[i]);
          EG_free(vid);
          return stat;
      }
  }
  EG_free(vid);
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static void
EG_centroid(meshMap *qm, int n, int *list, double *uvOUT, int usequv)
{
  int i, sumuv = 0,  sumxyz = 0;
  double  xyz0[18], xyz1[18], uv[2];

  if (usequv == 1) {
      uv[0] = uvOUT[0];
      uv[1] = uvOUT[1];
      if (EG_evaluate(qm->face, uv, xyz0) != EGADS_SUCCESS ) return ;
      sumuv++; sumxyz++;
  } else {
      uv  [0] =  uv [1] = 0.0;
      xyz0[0] = xyz0[1] = xyz0[2] = 0.0;
  }
  for (i = 0; i < n; i++) {
      if (qm->vType  [list[i] -1] *
          qm->valence[list[i] -1][2] == -2) continue;
      xyz0[0] += qm->xyzs[3 * (list[i] - 1)    ];
      xyz0[1] += qm->xyzs[3 * (list[i] - 1) + 1];
      xyz0[2] += qm->xyzs[3 * (list[i] - 1) + 2];
      sumxyz++;
      if (qm->vType [list[i]-1]  > 0 &&
          qm-> degen[list[i]-1] == 1) {
#ifdef DEBUG
              printf(" EG_centroid SKIP %d --> DEGENERATE\n",list[i]);
#endif
              continue;
      }
      uv[0]  += qm->uvs [2 * (list[i] - 1)    ];
      uv[1]  += qm->uvs [2 * (list[i] - 1) + 1];
      sumuv++;
  }
  if (sumxyz > 1 ) {
    xyz0[0] /= (double) sumxyz;
    xyz0[1] /= (double) sumxyz;
    xyz0[2] /= (double) sumxyz;
  }
  if (sumuv > 1 ) {
    uv [0]  /= (double) sumuv;
    uv [1]  /= (double) sumuv;
  }
  uvOUT[0] = uv[0];
  uvOUT[1] = uv[1];
  i        = EG_invEvaluateGuess(qm->face, xyz0, uv, xyz1);
  if (i   != EGADS_SUCCESS  || uv[0] < qm->range[0] ||
       uv[0] > qm->range[1] || uv[1] < qm->range[2] ||
       uv[1] > qm->range[3] ) {
#ifdef DEBUG
      printf(" FACE %d OUT OF RANGE %lf %lf %lf %lf\n", qm->fID, qm->range[0],
             qm->range[1], qm->range[2], qm->range[3]);
      printf("  %lf  %lf  %lf  %lf   %lf  GUESS\n", uvOUT[0], uvOUT[1], xyz0[0],
             xyz1[1], xyz0[2]);
      printf("  %lf  %lf  %lf  %lf   %lf  RES1  \n", uv[0], uv[1], xyz1[0],
             xyz1[1], xyz1[2]);
      printf("  %lf  %lf  %lf  %lf   %lf  RES1  \n", uv[0], uv[1], xyz1[0],
             xyz1[1], xyz1[2]);
      xyz1[0] = xyz1[1] = xyz1[2];
      i       = EG_invEvaluate(qm->face, xyz0, uv, xyz1);
      printf("  %d  %lf  %lf  %lf  %lf   %lf  RES2  \n", i, uv[0], uv[1],
             xyz1[0], xyz1[1], xyz1[2]);
#endif
  }
  uvOUT[0] = uv[0];
  uvOUT[1] = uv[1];
}


__HOST_AND_DEVICE__ static int
EG_angAtNPlane(meshMap *qm, int vC, int v1, int v2, double *normal, double *angle)
{
  int    stat, sing = 0, t1, t2, tc;
  double cross[3], pro1[3], pro2[3], r[2], xyz[18], v01[3],
         v02[3], norm1, norm2, c, uv[2], dot;

  *angle = 0.0;
  t1 = MIN(qm-> vType[v1 - 1], 0);
  t2 = MIN(qm-> vType[v2 - 1], 0);
  tc = MIN(qm-> vType[vC - 1], 0);
  if (tc == 0 && qm->degen[vC - 1] == 1) tc++;
  if (t1 == 0 && qm->degen[v1 - 1] == 1) t1++;
  if (t2 == 0 && qm->degen[v2 - 1] == 1) t2++;
  stat      = EG_evaluate(qm->face, &qm->uvs[2 * (vC - 1)], xyz);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_angAtNPlane EG_evaluate %d --> %d!!\n", vC, stat);
      return stat;
  }
  norm1     = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
  norm2     = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
  norm1     = sqrt(norm1 * norm2);
  // check that it is not singular point. Boundary vertices use also normal plane at a point inside.
  if (norm1 < qEPS || tc == 1) {
    sing = 1;
    uv[0] = uv[1] = 0.0;
    if (t1 <= 0) {
        uv[0] = qm->uvs[2 * (v1 - 1)    ];
        uv[1] = qm->uvs[2 * (v1 - 1) + 1];
    }
    if (t2 <= 0) {
        uv[0] = 0.5 * (uv[0] + qm->uvs[2 * (v2 - 1)    ]);
        uv[1] = 0.5 * (uv[1] + qm->uvs[2 * (v2 - 1) + 1]);
    }
    stat = EG_evaluate(qm->face, uv, xyz);
    if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
    printf(" INSIDE NORM1 < EPS EG_angAtNPlane EG_evaluate %d --> %d!!\n", vC, stat);
#endif
      return stat;
    }
    if (stat != EGADS_SUCCESS) return stat;
    norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
    norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
    norm1 = sqrt(norm1 * norm2);
    if (norm1 < qEPS) {
#ifdef DEBUG  // NOSE CASE HAS DEGENERATE QUAD EVEN AFTER EG_CENTROID. IT IS ENTERING HERE SOMETIMES
    // NEED TO FIX IT AT SOME POINT !!!!
    printf("EG_angleAtVnormalPlane VERT %d TYPE %d LINK TO %d %d "
           " du %lf %lf %lf dv %lf %lf %lf\n",
            vC, qm->vType[vC -1], v1, v2,
            xyz[3], xyz[4], xyz[5],
            xyz[6], xyz[7], xyz[8]);
#endif
        return EGADS_GEOMERR;
    }
  }
  v01[0] = xyz[3]; v01[1] = xyz[4]; v01[2] = xyz[5];
  v02[0] = xyz[6]; v02[1] = xyz[7]; v02[2] = xyz[8];
  if (qm->face->mtype == SREVERSE) {
      CROSS(v02, v01, normal);
  } else {
      CROSS(v01, v02, normal);
  }
  norm1      = sqrt(normal[0] *normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
  normal[0] /= norm1;
  normal[1] /= norm1;
  normal[2] /= norm1;
  normal[3] = xyz[0]; normal[4] = xyz[1]; normal[5] = xyz[2];
  c         = normal[0] * xyz[0] + normal[1] * xyz[1] + normal[2] * xyz[2];
  r[0]      = r[1] = 0.0;
  if (sing == 0) {
    dot = normal[0] * qm->xyzs[3 * (v1 - 1)    ] +
          normal[1] * qm->xyzs[3 * (v1 - 1) + 1] +
          normal[2] * qm->xyzs[3 * (v1 - 1) + 2] ;
    r[0] = (c - dot);
    dot = normal[0] * qm->xyzs[3 * (v2 - 1)    ] +
          normal[1] * qm->xyzs[3 * (v2 - 1) + 1] +
          normal[2] * qm->xyzs[3 * (v2 - 1) + 2] ;
    r[1] = (c - dot);
  }
  pro1[0] = qm->xyzs[3 * (v1 - 1)    ] + r[0] * normal[0];
  pro1[1] = qm->xyzs[3 * (v1 - 1) + 1] + r[0] * normal[1];
  pro1[2] = qm->xyzs[3 * (v1 - 1) + 2] + r[0] * normal[2];
	pro2[0] = qm->xyzs[3 * (v2 - 1)    ] + r[1] * normal[0];
  pro2[1] = qm->xyzs[3 * (v2 - 1) + 1] + r[1] * normal[1];
  pro2[2] = qm->xyzs[3 * (v2 - 1) + 2] + r[1] * normal[2];
  // re-center at vertex
  v01[0] = pro1[0] - qm->xyzs[3 * (vC - 1)    ];
  v01[1] = pro1[1] - qm->xyzs[3 * (vC - 1) + 1];
  v01[2] = pro1[2] - qm->xyzs[3 * (vC - 1) + 2];
  v02[0] = pro2[0] - qm->xyzs[3 * (vC - 1)    ];
  v02[1] = pro2[1] - qm->xyzs[3 * (vC - 1) + 1];
  v02[2] = pro2[2] - qm->xyzs[3 * (vC - 1) + 2];
  norm1  = sqrt( DOT(v01, v01) * DOT(v02, v02) );
  CROSS (v01, v02, cross);
  // projection gives length zero vectors; take just the "angle" as it is
  if (norm1 < qEPS) dot = 1.0;
  else              dot = DOT(v01, v02) / norm1;
  if (fabs(dot - 1.0) < EPS08) {
    *angle = 0.0;
    /*
    printf(" DIFFERENCE BETWEEN COORDINATES \n");
    printf(" %lf %lf %lf \n",qm->xyzs[3 * (v1 - 1)    ], qm->xyzs[3 * (v1 - 1)+1], qm->xyzs[3 * (v1 - 1) +2]);
    printf("%lf %lf %lf \n", proj1[0], proj1[1], proj1[2]);
    printf(" %lf %lf %lf \n",qm->xyzs[3 * (v2 - 1)    ], qm->xyzs[3 * (v2 - 1)+1], qm->xyzs[3 * (v2 - 1) +2]);
    printf("%lf %lf %lf \n", proj2[0], proj2[1], proj2[2]);
    printf(" ANGLE 0\n");*/
  }
  else if (fabs(dot + 1.0) < EPS08) *angle = PI;
  else                              *angle = acos(dot);
  norm1 = sqrt(DOT(cross, cross));
  if (norm1 > 1.e-10) {
       cross[0] /= norm1;
       cross[1] /= norm1;
       cross[2] /= norm1;
  }
  if (DOT(normal, cross) < 0.0) *angle =  2.0 * PI - *angle;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_angAtBdVert(meshMap *qm, int v, int *links, double *size)
{
  int i, j, k;
  double dist = 0.0, normal[6];

  *size = 0.0;
  if (qm->vType[v - 1]  <    0) return EGADS_INDEXERR;
  if (qm->star [v - 1] == NULL) {
      printf(" Looking at corners: buildstar %d is NULL \n ", v);
      return EGADS_MALLOC;
  }
  for (links[0] = links[1] = k = i = 0; i < qm->star[v-1]->nQ; i++) {
      j = qm->star[v-1]->verts[2 * i + 1] - 1;
      if (j >=0 && qm->vType[j] != -1) k++;
      if (qm->star[v-1]->quads[i] == -1) {
          links[1] = qm->star[v-1]->verts[                    2 * i + 1];
          links[0] = qm->star[v-1]->verts[qm->star[v-1]->idxV[2 * i + 3]];
      }
  }
  if (k >= 3) {
      *size = PI * 0.5; // boundary vertex is connected to more than two bounds (?)
      return EGADS_SUCCESS;
  } else if (k != 2) {
      printf(" EG_angAtBdVert:: vertex %d is at surface bounds and connected only to another boundary vertex !!\n ", v);
      return EGADS_GEOMERR;
  } else {
	  dist =  (qm->xyzs[3 * (links[1] - 1)]     - qm->xyzs[3 * (links[0] - 1)]    ) *
	          (qm->xyzs[3 * (links[1] - 1)]     - qm->xyzs[3 * (links[0] - 1)]    ) ;
	  dist += (qm->xyzs[3 * (links[1] - 1) + 1] - qm->xyzs[3 * (links[0] - 1) + 1]) *
	          (qm->xyzs[3 * (links[1] - 1) + 1] - qm->xyzs[3 * (links[0] - 1) + 1]) ;
      dist += (qm->xyzs[3 * (links[1] - 1) + 2] - qm->xyzs[3 * (links[0] - 1) + 2]) *
	          (qm->xyzs[3 * (links[1] - 1) + 2] - qm->xyzs[3 * (links[0] - 1) + 2]) ;
	  if (dist < 1.e-10) {
		  *size = 2.0 * PI;
#ifdef DEBUG
		  printf(" FACE %d V %d IS PERIODIC. ASSIGN ANGLE 2pi\n", qm->fID, v);
#endif
		  return EGADS_SUCCESS;
	  }
    return EG_angAtNPlane(qm, v, links[0], links[1], normal, &(*size));
  }
}


__HOST_AND_DEVICE__ static int EG_nValenceCount(meshMap *qm, int q, int n)
{
  int i, count, val;

  if (q <= 0 || q > qm->totQ ) {
      printf("EG_nValenceCount for a wrong index!! %d != [%d,%d]\n", q, 1, qm->totQ);
      return -1;
  }
  for (count = i = 0; i < 4; i++) {
      val    = qm->valence[qm->qIdx[4 * (q - 1) + i] - 1][1];
      if (val < 0) return val;
      if (n == val ||
         (n >    5 && val > 5)) count++;
  }
  return count;
}


__HOST_AND_DEVICE__ static double
EG_quadSize(meshMap *qm, int qID, /*@null@*/ double *tris)
{
  int iA, iB, iC, iD;
  double cross[3], vAB[3], vAC[3], vAD[3], qs[4];

  // Quad vertices
  iA  = qm->qIdx[4*(qID - 1)    ] - 1;
  iB  = qm->qIdx[4*(qID - 1) + 1] - 1;
  iC  = qm->qIdx[4*(qID - 1) + 2] - 1;
  iD  = qm->qIdx[4*(qID - 1) + 3] - 1;
  // first config ABC, ACD
  vAB[0] = qm->xyzs[3 * iB    ] - qm->xyzs[3 * iA    ];
  vAB[1] = qm->xyzs[3 * iB + 1] - qm->xyzs[3 * iA + 1];
  vAB[2] = qm->xyzs[3 * iB + 2] - qm->xyzs[3 * iA + 2];
  vAC[0] = qm->xyzs[3 * iC    ] - qm->xyzs[3 * iA    ];
  vAC[1] = qm->xyzs[3 * iC + 1] - qm->xyzs[3 * iA + 1];
  vAC[2] = qm->xyzs[3 * iC + 2] - qm->xyzs[3 * iA + 2];
  vAD[0] = qm->xyzs[3 * iD    ] - qm->xyzs[3 * iA    ];
  vAD[1] = qm->xyzs[3 * iD + 1] - qm->xyzs[3 * iA + 1];
  vAD[2] = qm->xyzs[3 * iD + 2] - qm->xyzs[3 * iA + 2];
  CROSS(vAB, vAC, cross);
  qs[0] = 0.5 * sqrt(DOT(cross, cross)) ;
  CROSS(vAC, vAD, cross);
  qs[1] = 0.5 * sqrt(DOT(cross, cross)) ;
  // second config BCD, BAD
  vAB[0] = qm->xyzs[3 * iC    ] - qm->xyzs[3 * iB    ];
  vAB[1] = qm->xyzs[3 * iC + 1] - qm->xyzs[3 * iB + 1];
  vAB[2] = qm->xyzs[3 * iC + 2] - qm->xyzs[3 * iB + 2];
  vAC[0] = qm->xyzs[3 * iD    ] - qm->xyzs[3 * iB    ];
  vAC[1] = qm->xyzs[3 * iD + 1] - qm->xyzs[3 * iB + 1];
  vAC[2] = qm->xyzs[3 * iD + 2] - qm->xyzs[3 * iB + 2];
  vAD[0] = qm->xyzs[3 * iA    ] - qm->xyzs[3 * iB    ];
  vAD[1] = qm->xyzs[3 * iA + 1] - qm->xyzs[3 * iB + 1];
  vAD[2] = qm->xyzs[3 * iA + 2] - qm->xyzs[3 * iB + 2];
  CROSS(vAB, vAC, cross);
  qs[2]  = 0.5 * sqrt(DOT(cross, cross)) ;
  CROSS(vAC, vAD, cross);
  qs[3]  = 0.5 * sqrt(DOT(cross, cross)) ;
  if (tris) { // used in area to get the triangle sizes ABC, BCD, CDA, DAB
      tris[0] = qs[0];
      tris[1] = qs[2];
      tris[2] = qs[1];
      tris[3] = qs[3];
  }
  return MIN(fabs(qs[0] + qs[1]), fabs(qs[2] + qs[3]));
}


/*
 * Returns <0 EGADS_ERROR
 * = 0         (QA0) area is good
 * = 1       (QA1) area has obtuse angles
 * = 200     (QA2) corner has obtuse angle or quad is very DELTA like
 * = 40000   (QA3) selfintersects
 * = 8000000 (QACB) crosses the domain boundary!
 * ratio = ABC / ACD (ideal =1 so triangle split is equal forming parallelogram)
 */

__HOST_AND_DEVICE__ static int EG_vertexArea(meshMap *qm, int vID)
  {
   int s, i, k, k1, vA, vB, vC, vD, count, doublet, qV[4], ori[4], bv[4], lr[2];
   int selfint = 1, cw, s1, s2, totArea = 0, ta, tb, td, iv = -1, bvp;
   double pV[12], c, cross[3],vn[6], vAB[3], vAC[3], vAD[3], tria[8], vr[2], t;
   double ang[4], qa[2], dot, pen = 1.0, ma = 0.0, maiv = 0.0, mba = 0.0;

   for (s = 0; s < qm->star[vID-1]->nQ; s++) {
     qm->star[vID-1]->area [s] = QA0;
     qm->star[vID-1]->ratio[s] = 1.0;
     qm->star[vID-1]->angle[s] = PI * 0.5;
     if(qm->star[vID-1]->quads[s] == -1) continue;
     (void) EG_quadSize(qm, qm->star[vID-1]->quads[s], tria);
 #ifdef DEBUG2
     printf(" TRIAREAS %lf %lf  %lf  %lf  \n",
            tria[0], tria[2], tria[1], tria[3]);
 #endif
     qV[0] = qm->qIdx[4 * (qm->star[vID-1]->quads[s] - 1 )    ] - 1;
     qV[1] = qm->qIdx[4 * (qm->star[vID-1]->quads[s] - 1 ) + 1] - 1;
     qV[2] = qm->qIdx[4 * (qm->star[vID-1]->quads[s] - 1 ) + 2] - 1;
     qV[3] = qm->qIdx[4 * (qm->star[vID-1]->quads[s] - 1 ) + 3] - 1;
     // project coordinates
     selfint =  1;
     iv      = -1;
     mba     = ma = 0.0;
     bvp     = 0;
     for (doublet = count = k = 0; k < 4; ++k) {
       ori[k]   =  1;
       ang[k]   =  0.0;
       bv [k]   = -1;
       if (qV[k] + 1 == vID) iv = k;
       if (qm->valence[qV[k]][2] * qm->vType[qV[k]] == -2) {
           ori[k]  = 1;
           doublet = 1;
       } else {
         vA =  k;
         vB = (k + 1)%4;
         vC = (k + 2)%4;
         vD = (k + 3)%4;
         ta = MIN(qm-> vType[qV[vA]], 0);
         tb = MIN(qm-> vType[qV[vB]], 0);
         td = MIN(qm-> vType[qV[vD]], 0);
         if (ta == 0 && qm->degen[qV[vA]] == 1) ta++;
         if (tb == 0 && qm->degen[qV[vB]] == 1) tb++;
         if (td == 0 && qm->degen[qV[vD]] == 1) td++;
         if (EG_angAtNPlane(qm, qV[vA] + 1, qV[vB] + 1, qV[vD] + 1,
                            vn, &ang[k]) != EGADS_SUCCESS ) {
                              ang[k] = PI * 0.5;
                              ori[k] = 1;
                              continue;
                           }
         if (ta < 1 && tb < 1 && td < 1) {
          c = vn[0] * vn[3] + vn[1] * vn[4] + vn[2] * vn[5];
          for (i = 0; i < 4; i++) {
            k1  = qV[(k + i)%4];
            dot = vn[0] * qm->xyzs[3 * k1    ]+
                  vn[1] * qm->xyzs[3 * k1 + 1]+
                  vn[2] * qm->xyzs[3 * k1 + 2];
            t   = (c - dot);
              pV[3 * i    ] = qm->xyzs[3 * k1    ] + t * vn[0];
              pV[3 * i + 1] = qm->xyzs[3 * k1 + 1] + t * vn[1];
              pV[3 * i + 2] = qm->xyzs[3 * k1 + 2] + t * vn[2];
          }
          vAB[0] = pV[3] - pV[0]; vAB[1] = pV[4]  - pV[1]; vAB[2] = pV[5]  - pV[2];
          vAC[0] = pV[6] - pV[0]; vAC[1] = pV[7]  - pV[1]; vAC[2] = pV[8]  - pV[2];
          vAD[0] = pV[9] - pV[0]; vAD[1] = pV[10] - pV[1]; vAD[2] = pV[11] - pV[2];
        } else {
          vAB[0] = qm->xyzs[3 *qV[vB]    ] - qm->xyzs[3 *qV[vA]    ];
          vAB[1] = qm->xyzs[3 *qV[vB] + 1] - qm->xyzs[3 *qV[vA] + 1];
          vAB[2] = qm->xyzs[3 *qV[vB] + 2] - qm->xyzs[3 *qV[vA] + 2];
          vAC[0] = qm->xyzs[3 *qV[vC]    ] - qm->xyzs[3 *qV[vA]    ];
          vAC[1] = qm->xyzs[3 *qV[vC] + 1] - qm->xyzs[3 *qV[vA] + 1];
          vAC[2] = qm->xyzs[3 *qV[vC] + 2] - qm->xyzs[3 *qV[vA] + 2];
          vAD[0] = qm->xyzs[3 *qV[vD]    ] - qm->xyzs[3 *qV[vA]    ];
          vAD[1] = qm->xyzs[3 *qV[vD] + 1] - qm->xyzs[3 *qV[vA] + 1];
          vAD[2] = qm->xyzs[3 *qV[vD] + 2] - qm->xyzs[3 *qV[vA] + 2];
        }
        CROSS(vAB, vAC, cross);
        if (DOT(vn, cross) < qEPS) ori[k] = -1;
        else {
             CROSS(vAC, vAD, cross);
             if (DOT(vn, cross) < qEPS) ori[k] = -1;
        }
        ma = MAX(ma, ang[k]);
        if ((qm->vType[qV[k]] > 0 && qm->bdAng[qV[k]] > PIEPS) ||
            (qm->vType[qm->valence[qV[k]][3] - 1] > 0 &&
             qm->bdAng[qm->valence[qV[k]][3] - 1] > PIEPS)) mba = MAX(mba, ang[k]);
             if ( qm->vType[qV[k]] != -1 ||
                 (qm->vType[qm->valence[qV[k]][3] -1] > 0 &&
                  qm->bdAng[qm->valence[qV[k]][3] -1] > PIEPS)) {
                   bv[k] = k;
                   if (qm->vType[qV[bv[k]]] == -1)  bvp = MAX (bvp, 1);
                   else {
                     bvp = MAX(qm->vType[qV[k]], bvp);
                     if (qm->bdAng[qV[k]] > PIEPS ) {//&& fabs(qm->bdAng[qV[k]] - 2.0 * PI) > qEPS) {
                    // IF QV[k] IS BOUNDARY VERTEX CHECK IF IT CROSSES DOMAIN
                    for (lr[0] = lr[1] = i = 0; i < qm->star[qV[k]]->nQ; i++) {
                    if ( qm->star[qV[k]]->quads[i] == -1) {
                        lr[1] = qm->star[qV[k]]->verts[2 * i + 1]                        - 1;
                        lr[0] = qm->star[qV[k]]->verts[qm->star[qV[k]]->idxV[2 * i + 3]] - 1;
                        break;
                      }
                    }
                    if (lr[1] != qV[(k +1)%4] && lr[0] != qV[(k +3)%4]) {
                      for (k1 = 1; k1 <= 2; k1++) {
                        vB = (k + k1    )%4;
                        vC = (k + k1 + 1)%4;
                        if (qV[vB] == lr[0] || qV[vB] == lr[1] ||
                            qV[vC] == lr[0] || qV[vC] == lr[1] )continue;
                        vAB[0] = qm->xyzs[3 * qV[vB]    ] - qm->xyzs[3 * qV[k]    ];
                        vAB[1] = qm->xyzs[3 * qV[vB] + 1] - qm->xyzs[3 * qV[k] + 1];
                        vAB[2] = qm->xyzs[3 * qV[vB] + 2] - qm->xyzs[3 * qV[k] + 2];
                        vAC[0] = qm->xyzs[3 * qV[vC]    ] - qm->xyzs[3 * qV[k]    ];
                        vAC[1] = qm->xyzs[3 * qV[vC] + 1] - qm->xyzs[3 * qV[k] + 1];
                        vAC[2] = qm->xyzs[3 * qV[vC] + 2] - qm->xyzs[3 * qV[k] + 2];
                        vAD[0] = qm->xyzs[3 *  lr[0]    ] - qm->xyzs[3 * qV[k]    ];
                        vAD[1] = qm->xyzs[3 *  lr[0] + 1] - qm->xyzs[3 * qV[k] + 1];
                        vAD[2] = qm->xyzs[3 *  lr[0] + 2] - qm->xyzs[3 * qV[k] + 2];
                        CROSS(vAB, vAC, cross);
                        cw = 0;
                        if (DOT(vn, cross) < 0.0 ) cw = 1;
                        // check that BC and CD dont cross the boundary
                        CROSS(vAB, vAD, cross);
                        s1 = s2 = 0;
                        if (DOT(vn, cross) < 0.0) s1 = 1;
                        CROSS(vAD, vAC, cross);
                        if (DOT(vn, cross) < 0.0) s2 = 1;
                        if ((cw == 0 && s1 == 0 && s2 == 0) ||
                            (cw == 1 && s1 == 1 && s2 == 1)) {
                                qm->star[vID-1]->area[s] = QACB;
#ifdef DEBUG
                          printf(" Alert! QUAD %d LINK %d %d CROSSES THE DOMAIN BOUNDS ! REJECT\n ",
                                  qm->star[vID-1]->quads[s], qV[vB] + 1, qV[vC] + 1);
                          printf(" LINK AB AL %d %d  %d %d = %d\n ",
                                  qV[k] + 1, qV[vB] + 1, qV[k] + 1,  lr[0] + 1, s1);
                          printf(" LINK AL AC %d %d  %d %d = %d\n ",
                                  qV[k] + 1,  lr[0] + 1, qV[k] + 1, qV[vC] + 1, s2);
                          printQuad(qm, qm->star[vID-1]->quads[s]);
#endif
                        }
                      }
                    }
                  }
                }
            }
         }
         if (ori[k] == -1) continue;
         if (selfint != 2) selfint = 0;
         count++;
       }
       qm->star[vID-1]->angle[s] = ang[iv];
       if (qm->star[vID-1]->area[s] == QA0 ) {
           if (selfint == 1 || count  < 2) qm->star[vID-1]->area[s] = QA3;
           else if (         count   != 4) qm->star[vID-1]->area[s] = QA1;
       }
       // now look for special cases: sharp interior corners, boundary quads...
       vr[0] = vr[1] = 0.0;
       qa[0] = tria[0] + tria[2];
       qa[1] = tria[1] + tria[3];
       maiv  = MAX(MAX(ang[iv], ang[(iv+1)%4]), ang[(iv+3)%4]);
       if (doublet == 0 && qa[0] > qEPS && qa[1] > qEPS &&
           qm->star[vID-1]->area[s] <= QA1) {
           k  =  iv    % 2; // id of corresponding diagonal split
           k1 = (iv + 3)%4; // id of triangle with centre iv
           if(tria[k] > qEPS && tria[k + 2] > qEPS)
             vr[0] = MIN(tria[k] / tria[k + 2], tria[k + 2] / tria[k]);
           if(tria[k1] > qEPS && tria[(k1 + 2)%4] > qEPS)
             vr[1] = MIN(tria[k1] / tria[(k1 + 2)%4], tria[(k1 + 2)%4] / tria[k1]);
           if (qm->star[vID-1]->area[s] == QA0) {
               pen = 1.0 - ERF(ANGCUT, 0.5 * ( PI + ANGCUT), maiv); // penalty kicks when angle > ANGCUT
 #ifdef DEBUG
               printf(" Angle %lf CUT %lf PENALTY %lf \n", maiv, ANGCUT, pen);
 #endif
               if (ang[(iv + 2) %4] > ANGCUT &&
                            ang[iv] < ANGCUT) qm->star[vID-1]->ratio[s] = vr[0];
               else qm->star[vID-1]->ratio[s] = MIN(vr[0], vr[1]);
           } else {
               vA = 0;
               if (qa[0] > qa[1] ) vA = 1;
               if (tria[k] + tria [k + 2] > qa[vA]) // MIN AREA = qa[vA]
                 qm->star[vID-1]->ratio[s] = -vr[0];
               else if (tria[(k1 + 2)%4] > qEPS &&
                        tria[k1] / tria[(k1 + 2)%4] < 1.0)
                 qm->star[vID-1]->ratio[s] = -vr[1];
               else
                 qm->star[vID-1]->ratio[s] = vr[0];
           }
       }
       if (selfint == 1 && bvp >= 4) qm->star[vID-1]->area[s] = QACB;
       if (qm->star[vID-1]->area[s] >= QA3)
         qm->star[vID-1]->ratio[s] = -qm->star[vID-1]->angle[s];
       if (doublet == 1) {
           if (qm->star[vID-1]->area [s] < QA3) {
               qm->star[vID-1]->area [s] = QA0;
               qm->star[vID-1]->ratio[s] = 1.0;
           }
       } else if (qm->star[vID-1]->area[s] < QA3 &&
                 ((bvp  > 1 && mba  >= PIEPS) ||
                  (bvp == 1 && maiv >= PIEPS))) {
 #ifdef DEBUG2
           printf(" bvp %d mba %lf maiv %lf\n", bvp, mba, maiv);
 #endif
           qm->star[vID-1]->area[s] = bvp * QA2;
       }
       if (qm->star[vID-1]->area [s] ==  QA1 &&
           qm->star[vID-1]->ratio[s] <  -0.5) qm->star[vID-1]->area[s] *= 2;
 #ifdef DEBUG2
       if (qm->star[vID-1]->area[s] != QA0) {
           printf("\n\n STAR %d SIGN AREA %d VAL %lf %lf -->RATIO %lf\n",
           vID, qm->star[vID-1]->area[s], qa[0],qa[1],
           qm->star[vID-1]->ratio[s]);
           if (qm->star[vID-1]->area[s] % QA2 == 0 )
               printf(" ************ SHARP CORNER %d = %lf ***************\n", bvp, ma);
           else if (qm->star[vID-1]->area[s] == QA3)
               printf(" ************ SELFINTERSECTS        ***************\n");
           else if (qm->star[vID-1]->area[s] == QACB)
               printf(" ************ CROSSES DOMAIN BOUNDS ***************\n");
           else if (qm->star[vID-1]->area[s] % QA1 == 0)
               printf(" ************ OBTUSE QUAD           ***************\n");
       }
 #endif
       qm->star[vID-1]->ratio[s] *= pen;
       if (qm->star[vID-1]->area [s] == QA0 && qm->star[vID-1]->ratio[s] < 1.e-05) {
           qm->star[vID-1]->area [s] = QA1;
           qm->star[vID-1]->ratio[s] = 0.0;
       }
 #ifdef DEBUG2
       printf("\n==========  star area %d Q %d BP %d AT V CENTRE %lf=======================\n",
       vID, qm->star[vID-1]->quads[s],  bvp,  qm->star[vID-1]->angle[s]);
       printf(" AREAS %d  penalty %lf ratio %lf\n",
              qm->star[vID-1]->area[s], pen, qm->star[vID-1]->ratio[s]);
       printf(" ratios TRIS %lf %lf = %lf  %lf %lf = %lf  -- \n ",
              tria[0], tria[2], vr[0], tria[1], tria[3], vr[1]);
       printf(" Angles (%d) = %lf (%d) = %lf (%d) = %lf (%d) = %lf \n ",
       qV[0] + 1, ang[0], qV[1] + 1, ang[1],
       qV[2] + 1, ang[2], qV[3] + 1, ang[3]);
       printf("\n=================================\n");
 #endif
       totArea  += qm->star[vID-1]->area[s];
   }
   return totArea;
 }


__HOST_AND_DEVICE__ static int
EG_placeVertex(meshMap *qm, int vID, int pass, int full)
{
  int q, v, j, i, k, *vl = NULL, ta[2], v1, v2, v3, corner, vla = 0,
      la = -1, lb, nt, block = 0, round, doublet = 0, bt = -1, i0;
  double uv[4], res[6], angopt =0.0, angpen, errbd, ptol = 0.25, auxd, normal[6];
  int *qInfo = NULL;
#ifdef DEBUG2
  char buffer[100];
  int d1, d2, update = 0;
  double p[2], pos[18];
  FILE *fout;
#endif

  v = vID - 1;
  if (qm->vType[v] != -1 || qm->star[v] == NULL) {
#ifdef DEBUG
      printf(" VERTEX %d is boundary star %p!\n ", vID, qm->star[v]);
#endif
      return QA0;
  }
#ifdef DEBUG2
  snprintf(buffer, 100,"PLACEStart_%d_%d", qm->plotcount, vID);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (d1 = 0; d1 < qm->star[v]->nQ; d1++) {
          p[0] = p[1] = 0.0;
          for (d2 = 0; d2 <= 4; d2++) {
              k = qm->qIdx[ 4 * (qm->star[v]->quads[d1] -1) + d2%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*k  ],
              qm->xyzs[3*k + 1], qm->xyzs[3*k + 2], k + 1 );
              if ( d2 == 4 ) break;
              p[0] += 0.25 * qm->uvs[2 * k    ];
              p[1] += 0.25 * qm->uvs[2 * k + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, p, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
          qm->star[v]->quads[d1]);
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
  if (qm->valence[v][2] == 2) {
      i     = qm->valence[v][3] - 1;
      j     = qm->valence[v][4] - 1;
      uv[0] = 0.5 * (qm->uvs[2 * i    ] + qm->uvs[2 * j    ]);
      uv[1] = 0.5 * (qm->uvs[2 * i + 1] + qm->uvs[2 * j + 1]);
      updateVertex(qm, vID, uv);
      return QA0;
  }
  qInfo  = (int*) EG_alloc(2 * qm->star[v]->nQ * sizeof(int));
  vl     = (int*) EG_alloc(    qm->star[v]->nV * sizeof(int));
  if (qInfo == NULL || vl == NULL) return EGADS_MALLOC;
  corner = k = -1;
  bt     = qm->star[v]->type;
  if (bt >= 0) {
      if(qm->vType[qm->valence[v][3]-1] > 0 &&
         qm->bdAng[qm->valence[v][3]-1] > PIEPS) k = 1;
      else k = 0;
  }
  for (i = 0; i < qm->star[v]->nQ; i++) {
      qInfo[i] = k;
      vl   [i] = 0;
      v1       = qm->star[v]->verts[qm->star[v]->idxV[2 *i + 1]] - 1;
      v2       = qm->star[v]->verts[qm->star[v]->idxV[2 *i + 2]] - 1;
      v3       = qm->star[v]->verts[qm->star[v]->idxV[2 *i + 3]] - 1;
      if (qm->vType[v1] > 0) {
          if (qm->bdAng[v1] > PIEPS) {
              qInfo[i] = 3;
              if (la == -1) la = 2 * i + 1;
          }
      }
      if (qm->vType[v2] > 0) {
          if (qm->bdAng[v2] > PIEPS) qInfo[i] = 3;
          else                       qInfo[i] = MAX(qInfo[i], 2);
          if (qm->vType[v2] == 2)      corner = i;
      }
      if (qm->vType[v3] > 0) {
          if (qm->bdAng[v3] > PIEPS) qInfo[i] = 3;
          else                       qInfo[i] = MAX(qInfo[i], 2);
      }
      if (    qm->vType[qm->valence[v1][3] -1] > 0 ) {
          if (qm->bdAng[qm->valence[v1][3] -1] > PIEPS)
               qInfo[i] = MAX(qInfo[i], 1);
          else qInfo[i] = MAX(qInfo[i], 0);
      }
      if (    qm->vType[qm->valence[v2][3] -1] > 0 ) {
          if (qm->bdAng[qm->valence[v2][3] -1] > PIEPS)
               qInfo[i] = MAX(qInfo[i], 1);
          else qInfo[i] = MAX(qInfo[i], 0);
      }
      if (    qm->vType[qm->valence[v3][3] -1] > 0 ) {
          if (qm->bdAng[qm->valence[v3][3] -1] > PIEPS)
               qInfo[i] = MAX(qInfo[i], 1);
          else qInfo[i] = MAX(qInfo[i], 0);
      }
      j = qm->valence[v ][1] * qm->valence[v1][1] *
          qm->valence[v2][1] * qm->valence[v3][1];
      qInfo[qm->star[v]->nQ + i] = 0;
      if      (j%3 == 0) qInfo[qm->star[v]->nQ + i] = -1;
      else if (j   == 4 * 4 * 4 * 4 ||
               (qm->star[v]->nQ != 4 && bt == -1)) qInfo[qm->star[v]->nQ + i] = 1;
#ifdef DEBUG2
      printf(" BT %d QUAD %d FROM CENTRE %d IS TYPE %d VAL4 %d \n",
      bt, qm->star[v]->quads[i], v + 1, qInfo[i],  qInfo[qm->star[v]->nQ + i]);
#endif
  }
  uv [0] = qm->uvs[2 * v    ];
  uv [1] = qm->uvs[2 * v + 1];
  angopt = 0.0;
  k      = 0;
  if      (la     != -1) {
    k = la;
    vla = qm->star[v]->verts[la] - 1;
  }
  else if (corner >=  0) k = 2 * corner + 2;
  if (k > 0) {
      angopt  = qm->bdAng[qm->star[v]->verts[k] - 1];
      angopt /= (double)(qm->valence[qm->star[v]->verts[k]-1][2] - 1);
#ifdef DEBUG2
      printf(" ANGLE AT BOUNDARY VERTEX %d is %lf VALENCE %d \t",
      qm->star[v]->verts[k], angopt, qm->valence[qm->star[v]->verts[k]-1][2]);
#endif
  }
  block   = 0;
  nt      = qm->star[v]->nV;
  if (bt == 3 && full == 0) full = 1;
  ta[0]   = ta[1] = QA0;
  for (i0 = round = 0; round < 2; round++) {
      for (q = 0; q <= nt; q++) {
          if (i0 == 1) {
              if (q > 0 && (block > 0 || qm-> star[v]->nQ == 3 )) break;
#ifdef DEBUG2
              printf ("\n\n PLACEVERT ROUND %d la %d q %d / %d PASS %d block %d\n ",
              round, la, q, nt, pass, block);
#endif
              updateVertex(qm, vID, uv);
              if (q == 0) {
                uv[2] = uv[0];
                uv[3] = uv[1];
                EG_centroid (qm, qm->star[v]->nQ, &qm->valence[v][3], &uv[2], (round+1)%2);
                if (la != -1 && ta[0] == QA0 && qm->bdAng[vla] > 3.9 &&
                       qm->degen[vla] == 0) {
                    if (qm->star[v]-> nQ != 4 || full != 2) continue;
                    uv[2] = 0.5 * uv[0] + 0.5 * uv[2];
                    uv[3] = 0.5 * uv[1] + 0.5 * uv[3];
                  }
              } else if (round == 0) {
                  if ((    qm->star[v]->nQ %2 == 0 && q > qm->star[v]->nQ / 2) ||
                      (q > qm->star[v]->nQ )) break;
                  vl[0]  = qm->star[v]->verts[                  2 * (q - 1) + 1 ] - 1;
                  vl[1]  = qm->star[v]->verts[qm->star[v]->idxV[2 * (q - 1) + 5]] - 1;
                  j      = qm->star[v]->verts[qm->star[v]->idxV[2 * (q - 1) + 7]] - 1;
                  if ( j == vl[0] ||
                  qm->vType[vl[0]] * qm->valence[vl[0]][2] == -2 ||
                  qm->vType[vl[1]] * qm->valence[vl[1]][2] == -2 ) continue;
                  if (la != -1 ) {
                      v1  = qm->star[v]->idxV[la + 2];
                      v2  = qm->star[v]->idxV[la + (qm->star[v]->nV -1)- 2];
                      if ((vl[0] == qm->star[v]->verts[v1] - 1||
                           vl[0] == qm->star[v]->verts[v2] - 1) &&
                          (vl[1] == qm->star[v]->verts[v1] - 1  ||
                           vl[1] == qm->star[v]->verts[v2] - 1)) continue;
                  }
                  if ((qm->vType[vl[0]] > 0 && qm->degen[vl[0]] == 1) ||
                      (qm->vType[vl[1]] > 0 && qm->degen[vl[1]] == 1)) {
                        vl[0]++; vl[1]++;
                        EG_centroid(qm, 2, vl, &uv[2], 0);
                   } else
                      EG_getSidepoint(qm->face, 0.5, &qm->uvs[2 * vl[0]],
                                     &qm->uvs[2 * vl[1]], NULL, NULL, &uv[2]);
              } else if (la != -1) {
                  lb = qm->star[v]->idxV [la + q];
                  if( (abs(lb - la) < 3 || abs(lb - la) > (nt - 1) - 3 ) ||
                      (ta[0] == QA0 && q % 2 == 1))  continue;
                  vl[0] = qm->star[v]->verts[la] - 1;
                  vl[1] = qm->star[v]->verts[lb] - 1;
                  if ((qm->vType[vl[0]] > 0 && qm->degen[vl[0]] == 1) ||
                      (qm->vType[vl[1]] > 0 && qm->degen[vl[1]] == 1)) {
                        vl[0]++; vl[1]++;
                        EG_centroid(qm, 2, vl, &uv[2], 0);
                   }
                   else EG_getSidepoint(qm->face, 0.5, &qm->uvs[2 * vl[0]],
                                       &qm->uvs[2 * vl[1]], NULL, NULL, &uv[2]);
              } else {
                  if (q - 1 >= qm->star[v]->nQ || ta[0] == QA0 ) break;
                  vl[0]  = vID;
                  i      = 1;
                  for (k = 0; k < qm->star[v]->nQ; k++) {
                      j  = qm->star[v]->verts[2 * k + 1] -1;
                      if (k == q - 1 || k == q || qm->vType[j] >   0 ||
                          qm->vType[j] * qm->valence[j][2] == -2 ) continue;
                      vl[i++] = j + 1;
                  }
                  if (i == 2) continue;
                  EG_centroid(qm, i, vl, &uv[2], 0);
#ifdef DEBUG2
                  printf(" CENTROID EXCLUDING QUAD %d \n", qm->star[v]->quads[q - 1]);
#endif
              }
              updateVertex(qm, vID, &uv[2]);
          }
          ta[i0]          = EG_vertexArea(qm, vID);
          angpen          = 0.0;
          res[3 * i0    ] = 1.0;
          res[3 * i0 + 1] = 1.0;
          res[3 * i0 + 2] = 1.0;
          for (j = 0 ;j < qm->star[v]->nQ; j++ ) {
#ifdef DEBUG2
              printf(" BT %d QINFO %d =  %d  QTYPE %d \n",
                       bt, qm->star[v]->quads[j], qInfo[j], qInfo[j + qm->star[v]->nQ]);
#endif
              if (qm->star[v]->ratio[j] < 0.0)
                  res[3 * i0] = MIN (res[3 * i0], qm->star[v]->ratio[j]);
              else {
                  errbd  = -1.0;
                  angpen = MAX(angpen, qm->star[v]->angle[j]);
                  if (corner == j) {
                      errbd  = fabs(angopt - qm->star[v]->angle[corner]);
#ifdef DEBUG2
                      printf(" CORNER --> ANGOPT %lf AT BOUND %lf ERRBD %lf\n",
                              angopt, qm->star[v]->angle[corner], errbd);
#endif
                 } else if (2 * j + 1 == la && qm->degen[vla] == 0 &&
                   EG_angAtNPlane (qm, vla + 1, qm->star[v]->verts[la + 1], v + 1,
                                   normal, &auxd) == EGADS_SUCCESS) {
                                     errbd = fabs(angopt - auxd);
#ifdef DEBUG2
  printf (" LA = 2j+1 %d %d %d --> ANGOPT %lf ANG %lf\t",
            vla + 1, qm->star[v]->verts[la + 1], v + 1, angopt, auxd);
#endif
                 } else if (qm->star[v]->idxV[2 * j + 3] == la &&
                           qm->degen[vla] == 0) {
                      v1 = qm->star[v]->idxV[la + (qm->star[v]->nV -1) - 1];
                      if (EG_angAtNPlane(qm, vla + 1, v + 1, qm->star[v]->verts[v1],
                                         normal, &auxd) == EGADS_SUCCESS) {
                      errbd = fabs(angopt - auxd);
#ifdef DEBUG2
printf("LA = 2j + 3 %d %d %d --> ANGOPT %lf ANG %lf\t",
        vla + 1,vID, qm->star[v]->verts[v1], angopt, auxd);
#endif

                    }
                } else if (qm->star[v]->nQ == 4 && qm->star[v]->angle[j] < ANGCUT &&
                           qm->star[v]->angle[(2+j)%4] < ANGCUT ) {
                             errbd = fabs(qm->star[v]->angle[j] - qm->star[v]->angle[(2 + j)%4]);
#ifdef DEBUG2
printf(" OPP ANGLES --> %lf %lf --> %lf \n",
        qm->star[v]->angle[j], qm->star[v]->angle[(2 + j)%4],
        fabs(qm->star[v]->angle[j]- qm->star[v]->angle[(2 + j)%4]));
#endif
                }
                auxd = fabs(1.0 - ERF(ANGCUT,PI, angpen));
#ifdef DEBUG2
printf(" DOMINANT ANGLE %lf MAX ANGLE %lf\n", errbd, auxd);
#endif
                  if (errbd > 0.0) errbd = MIN (auxd, (PI - errbd) / PI);
                  else             errbd = auxd;
                  if  (qInfo[qm->star[v]->nQ + j] == -1 && qm->star[v]->area[j] == QA0) continue;
                  if (errbd < EPS08) {
                      res[3 *i0 + 1] = MIN (res[3 *i0 + 1], -angpen);
                      res[3 *i0 + 2] = MIN (res[3 *i0 + 2], -angpen);
                  } else if  (qInfo[j] >= bt ||
                              (bt < 3   && qm->star[v]-> nQ == 4 &&
                                 qInfo[qm->star[v]->nQ + j] == 1))
                              res[3 *i0 + 1] = MIN (res[3 *i0 + 1], errbd * qm->star[v]->ratio[j]);
                  else res[3 *i0 + 2] = MIN (res[3 *i0 + 2], errbd * qm->star[v]->ratio[j]);
                  /*else if  (qInfo[j] >= bt || (qm->star[v]-> nQ == 4 && qInfo[qm->star[v]->nQ + j] == 1))
                       res[3 *i0 + 1] = MIN (res[3 *i0 + 1], errbd * qm->star[v]->ratio[j]);
                  else res[3 *i0 + 2] = MIN (res[3 *i0 + 2], errbd * qm->star[v]->ratio[j]);*/
#ifdef DEBUG2
  printf(" j %d AREA %d RATIO %lf ERRBD %lf res %lf %lf %lf\n",
           j, qm->star[v]->area[j], qm->star[v]->ratio[j], errbd,
           res[3*i0], res[3 * i0 + 1], res[3 * i0 + 2]);
#endif
              }
          }
          if (i0 == 0) {
              i0 = 1;
              q--;
              continue;
          }
#ifdef DEBUG2
          if (q == 0 ) printf(" CENTROID %d\n", round);
          else printf(" \n\n q - %d / %d MID POINTS BETWEEN %d %d \n ",q ,nt,  vl[0], vl[1]);
          printf(" NOW    AREA %d NEG %lf PRIORITY %lf OTHER %lf\n", ta[1], res[3], res[4], res[5]);
          printf(" STORED AREA %d NEG %lf PRIORITY %lf OTHER %lf\n", ta[0], res[0], res[1], res[2]);
          snprintf(buffer, 100,"PLACE_%d_%d_%d_%d", qm->plotcount, vID, q, round);
          printf("Writing in %s\n ", buffer);
          fout = fopen(buffer,"w");
          if (fout != NULL ) {
              for (d1  = 0; d1 < qm->star[v]->nQ; d1++ ) {
                  p[0] = p[1] = 0.0;
                  for (d2 = 0; d2 <= 4; d2++) {
                      k   = qm->qIdx[ 4 * (qm->star[v]->quads[d1] -1) + d2%4 ] - 1;
                      fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*k  ],
                      qm->xyzs[3*k + 1], qm->xyzs[3*k + 2], k + 1 );
                      if (d2 == 4) break;
                      p[0] += 0.25 * qm->uvs[2 * k    ];
                      p[1] += 0.25 * qm->uvs[2 * k + 1];
                  }
                  fprintf(fout,"\n\n");
                  EG_evaluate(qm->face, p, pos);
                  fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
                  qm->star[v]->quads[d1]);
                  fprintf(fout,"\n\n");
              }
              fclose (fout);
          }
#endif
          // compare areas and see if we should update
          k = 0;
          /*if (pass == -1) {
              if (la == -1 && (ta[1] < QA2 || (qm->pp == 1 && ta[1] <= QA2))) {
	            k = 2;
                block = 1;
	          }
          }*/
          if (pass == -1) {
              if ((bt == 2 &&  ta[1] == QA0) ||
                  (bt  < 2 && (ta[1]   < QA2 || (qm->pp == 1 && ta[1] <= QA2)))) {
	            k = 2;
                block = 1;
	          }
          }
          else if (ta[0] <  ta[1]) continue;
          else if (ta[1] <  ta[0])  k = 1;
          else if (ta[1] != QA0) {
               if (ta[1] >= QA3 && la != -1 )continue;
               else if (res[3] <    0.0) {
                    if (res[3] > res[0]) k = 1;
              }else if (res[0] < 0.0 || res[4] + res[5] > res[1] + res[2]) k = 1;
          } else {
#ifdef DEBUG2
              printf(" LOOKING CAREFULLY :: PASS %d  FULL %d \n", pass, full);
#endif
              if (fabs (res[1] - 1.0) < 1.e-10 &&
			      fabs (res[2] - 1.0) < 1.e-10 &&
				  fabs (res[4] - 1.0) < 1.e-10 &&
				  fabs (res[5] - 1.0) < 1.e-10 ) k = 2;
			  else if (q == 0 && res[4] >= ptol && res[5] >= ptol &&
					  (bt < 3 || full   != 1 )) k = 2;
              if (k == 0) {
                  if (res[4] < ptol) {
                      if (res[1] > res[4]) continue;
                      k = 1;
                  } else {
                      if      (res[1] < ptol) k = 1;
                      else if (res[5] < ptol) {
                           if (res[2] > res[5]) continue;
                           k = 1;
                      }
					  else if (fabs(res[1] - res[4]) < qEPS) {
                        //if (res[2] > res[5]) continue;
						//k = 1;
						if (res[2] < res[5]) k = 1;
						else continue;
					  } else if (fabs(res[2] - res[5]) < qEPS) {
                        if(res[1] < res[4]) k = 1;
						else continue;
					  }
					  else if (res[2] < ptol || res[1] < ptol || res[4] > res[1]) k = 1;
                  }
              }
          }
          if (k == 0) {
            updateVertex(qm, vID, uv);
            continue;
          }
#ifdef DEBUG2
          update = 1;
          printf(" k = %d pass %d q %d r %d\n", k, pass, q , round);
          printf(" NEW  AREA %d RATIOS NEG %lf PIV %lf SEC %lf\n", ta[1], res[3], res[4], res[5]);
          printf(" BEST AREA %d RATIOS NEG %lf PIV %lf SEC %lf\n", ta[0], res[0], res[1], res[2]);
#endif
          uv [0] = uv [2];
          uv [1] = uv [3];
          res[0] = res[3];
          res[1] = res[4];
          res[2] = res[5];
          ta [0] = ta [1];
          if (ta[0] != QA0) continue;
          if (doublet == 1 || k == 2 || q == 0 || full == 10) block = 1;
      }
      if (block > 0) break;
  }
  updateVertex(qm, vID, uv);
#ifdef DEBUG2
  printf(" UPDATING AREAS ???? %d \n", update);
  snprintf(buffer, 100,"PLACE_%d_%d", qm->plotcount, vID);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (d1  = 0; d1 < qm->star[v]->nQ; d1++ ) {
          p[0] = p[1] = 0.0;
          for (d2 = 0; d2 <= 4; d2++) {
              k   = qm->qIdx[ 4 * (qm->star[v]->quads[d1] -1) + d2%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*k  ],
              qm->xyzs[3*k + 1], qm->xyzs[3*k + 2], k + 1 );
              if (d2 == 4) break;
              p[0] += 0.25 * qm->uvs[2 * k    ];
              p[1] += 0.25 * qm->uvs[2 * k + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, p, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
          qm->star[v]->quads[d1]);
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
  EG_free(vl);
  EG_free(qInfo);
  return ta[0];
}


__HOST_AND_DEVICE__ int EG_createMeshMap(bodyQuad *bodydata)
{
  int          f, stat = 0, j, q, k, kk, kOK, len, iA, iB, iC;
  int          ntri, nquad, e4[4] ;
  const int    *tris, *tric, *ptype, *pindex;
  double       angle, xyz[18], norm1, norm2, u01[3], u02[3], v01[3], v02[3], uC[3], vC[3];
  const double *xyzs, *uvs;
  int          qV[6]    = { 0, 1, 2, 5, 0, 1};
  int          qLoop[5] = { 0, 1, 2, 3, 0   };
  ego          tess;
  egTessel     *btess;

  tess  = bodydata->tess;
  btess = (egTessel *) tess->blind;

  bodydata->qm = (meshMap **) EG_alloc(bodydata->nfaces * sizeof(meshMap *));
  if (bodydata->qm == NULL) {
    printf(" EG_createMeshMap Failed to allocate quad map pointer!\n");
    return  EGADS_MALLOC;
  }
  for (f = 0; f < bodydata->nfaces; ++f) {
      bodydata->qm[f] = (meshMap *) EG_alloc(sizeof(meshMap));
      if (bodydata->qm[f] == NULL) {
          printf("EG_createMeshMap face %d VOID ptr\n ", f + 1);
          continue;
      }
      bodydata->qm[f]->fID        = f + 1;
      bodydata->qm[f]->plotcount  = 0;
      /* Edges associated to face */
      stat = EG_getTessFace(bodydata->tess, f + 1, &len,
                            &xyzs, &uvs, &ptype, &pindex, &ntri,
                            &tris, &tric);
      if (stat != EGADS_SUCCESS) {
          printf("EG_createMeshMap :: EG_getTessFace %d = %d !!\n",f + 1, stat);
          bodydata->qm[f]->fID   = 0;
      }
      // CHECK TESSELATION IS WELL CONSTRUCTED
      for (q = 0; q < ntri; q++) {
         iA     = tris[3 * q     ] - 1;
         iB     = tris[3 * q  + 1] - 1;
         iC     = tris[3 * q  + 2] - 1;

         u01[0] = xyzs[3 * iB    ] - xyzs[3 * iA    ];
         u01[1] = xyzs[3 * iB + 1] - xyzs[3 * iA + 1];
         u01[2] = xyzs[3 * iB + 2] - xyzs[3 * iA + 2];

         u02[0] = xyzs[3 * iC    ] - xyzs[3 * iA    ];
         u02[1] = xyzs[3 * iC + 1] - xyzs[3 * iA + 1];
         u02[2] = xyzs[3 * iC + 2] - xyzs[3 * iA + 2];
         CROSS(u01, u02, uC);
         norm1   = sqrt(uC[0] * uC[0] + uC[1] * uC[1] + uC[2] * uC[2]);
         if ( norm1 == 0. ) {
            printf(" \n EG_createMeshMap FACE %d TRI %d, zero area !!!\n",
            f +1, q + 1 );
            bodydata->qm[f]->fID       = 0;
         }
         uC[0] /= norm1;
         uC[1] /= norm1;
         uC[2] /= norm1;
         for (kk = 0; kk < 3; kk++) {
           k      = tric[3 * q + kk] - 1;
           if (k < 0 || k < q) continue;
           iA     = tris[3 * k     ] - 1;
           iB     = tris[3 * k  + 1] - 1;
           iC     = tris[3 * k  + 2] - 1;

           v01[0] = xyzs[3 * iB    ] - xyzs[3 * iA    ];
           v01[1] = xyzs[3 * iB + 1] - xyzs[3 * iA + 1];
           v01[2] = xyzs[3 * iB + 2] - xyzs[3 * iA + 2];

           v02[0] = xyzs[3 * iC    ] - xyzs[3 * iA    ];
           v02[1] = xyzs[3 * iC + 1] - xyzs[3 * iA + 1];
           v02[2] = xyzs[3 * iC + 2] - xyzs[3 * iA + 2];
           CROSS(v01, v02, vC);
           norm1  = sqrt(vC[0] * vC[0] + vC[1] * vC[1] + vC[2] * vC[2]);
           if ( norm1 == 0. ) {
              printf(" \n EG_createMeshMap FACE %d TRI %d, zero area !!!\n",
              f +1, k + 1 );
              bodydata->qm[f]->fID       = 0;
           }
           vC[0] /= norm1;
           vC[1] /= norm1;
           vC[2] /= norm1;

           if (DOT(uC, vC) < -EPS08 ) {
              printf(" \n EG_createMeshMap FACE %d TRI %d, %d %.12f wrong OR !!!\n",
              f +1, q + 1, k + 1, DOT(uC, vC) );
              bodydata->qm[f]->fID       = 0;
           }
        }
      }
      nquad                   = (int) ntri/2;
      bodydata->qm[f]->oriV   = len;
      bodydata->qm[f]->oriQ   = nquad;
      bodydata->qm[f]->sizeV  = 2 * len;
      bodydata->qm[f]->sizeQ  = 2 * nquad;
      bodydata->qm[f]->totV   = len;
      bodydata->qm[f]->totQ   = nquad;
      bodydata->qm[f]->vInv   = NULL;
      bodydata->qm[f]->regBd  = 1;
      if( btess->tess2d[f].tfi == 1 || bodydata->qm[f]->fID == 0) {
          bodydata->qm[f]->xyzs    = NULL;
          bodydata->qm[f]->uvs     = NULL;
          bodydata->qm[f]->vType   = NULL;
          bodydata->qm[f]->qIdx    = NULL;
          bodydata->qm[f]->qAdj    = NULL;
          bodydata->qm[f]->remQ    = NULL;
          bodydata->qm[f]->remV    = NULL;
          bodydata->qm[f]->valence = NULL;
          bodydata->qm[f]->star    = NULL;
          bodydata->qm[f]->bdAng   = NULL;
          bodydata->qm[f]->degen   = NULL;
          bodydata->qm[f]->fID     = 0;
          continue;
      }
      bodydata->qm[f]->xyzs    = (double *) EG_alloc(3*(2 * len  )*sizeof(double));
      bodydata->qm[f]->uvs     = (double *) EG_alloc(2*(2 * len  )*sizeof(double));
      bodydata->qm[f]->vType   = (int    *) EG_alloc(  (2 * len  )*sizeof(  int ));
      bodydata->qm[f]->qIdx    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->qAdj    = (int    *) EG_alloc(4*(2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->remQ    = (int    *) EG_alloc(  (2 * nquad)*sizeof(  int ));
      bodydata->qm[f]->remV    = (int    *) EG_alloc(  (2 * len  )*sizeof(  int ));
      bodydata->qm[f]->valence = (int   **) EG_alloc(  (2 * len  )*sizeof(  int*));
      bodydata->qm[f]->star    = (vStar **) EG_alloc(  (2 * len  )*sizeof(vStar*));
      if (bodydata->qm[f]->qIdx  == NULL || bodydata->qm[f]->qAdj    == NULL ||
          bodydata->qm[f]->xyzs  == NULL || bodydata->qm[f]->uvs     == NULL ||
          bodydata->qm[f]->vType == NULL || bodydata->qm[f]->remQ    == NULL ||
          bodydata->qm[f]->remV  == NULL || bodydata->qm[f]->valence == NULL ||
          bodydata->qm[f]->star  == NULL) {
          bodydata->qm[f]->fID = 0;
          continue;
      }
      bodydata->qm[f]->remQ[0]  = 0;
      bodydata->qm[f]->remV[0]  = 0;
      bodydata->qm[f]->invsteps = 0;
      for (j = 0; j < 2 * len; j++) {
          bodydata->qm[f]->star[j]    = NULL;
          bodydata->qm[f]->valence[j] = (int *) EG_alloc(3 * sizeof(int));
          if (bodydata->qm[f]->valence[j] == NULL) {
              bodydata->qm[f]->fID = 0;
              continue;
          }
          bodydata->qm[f]->valence[j][2] = 0;
      }
      if (bodydata->qm[f]->fID == 0) continue;
      bodydata->qm[f]->face = bodydata->faces[f];
      stat = EG_getRange(bodydata->qm[f]->face, bodydata->qm[f]->range, &j);
      if (stat != EGADS_SUCCESS)
          printf(" EG_createMeshMap :: getRange = %d !! \n", stat);
      bodydata->qm[f]->range[0] -= 1.e-4;
      bodydata->qm[f]->range[1] += 1.e-4;
      bodydata->qm[f]->range[2] -= 1.e-4;
      bodydata->qm[f]->range[3] += 1.e-4;
      for (k = j = 0; j < len; j++) {
          bodydata->qm[f]->vType  [j]         = ptype[j];
          if (ptype[j] >= 0 ) k++;
          bodydata->qm[f]->uvs    [2 * j    ] = uvs[2 * j    ];
          bodydata->qm[f]->uvs    [2 * j + 1] = uvs[2 * j + 1];
          if (uvs[2 * j    ] < bodydata->qm[f]->range[0] ||
              uvs[2 * j    ] > bodydata->qm[f]->range[1] ||
              uvs[2 * j + 1] < bodydata->qm[f]->range[2] ||
              uvs[2 * j + 1] > bodydata->qm[f]->range[3] ) {
              printf(" EG_createMeshMap :: vertex %d = %f  %f  is out of range !! \n ",
                     j + 1, uvs[2 * j], uvs[2 * j + 1]);
              printf(" Range FACE %d --- > %f  %f  %f  %f\n ", bodydata->qm[f]->fID,
                     bodydata->qm[f]->range[0],bodydata->qm[f]->range[1],
                     bodydata->qm[f]->range[2],bodydata->qm[f]->range[3]);
              bodydata->qm[f]->fID = 0;
              break;
          }
          bodydata->qm[f]->xyzs[3 * j    ] = xyzs[3 * j    ];
          bodydata->qm[f]->xyzs[3 * j + 1] = xyzs[3 * j + 1];
          bodydata->qm[f]->xyzs[3 * j + 2] = xyzs[3 * j + 2];
      }
      bodydata->qm[f]->bdAng = (double *) EG_alloc(k * sizeof(double));
      bodydata->qm[f]->degen = (int    *) EG_alloc(k * sizeof(   int));
      if (bodydata->qm[f]->bdAng == NULL ||
          bodydata->qm[f]->degen == NULL){
          bodydata->qm[f]->fID = 0;
          continue;
      }
      for (j = 0; j < len; j++) {
          if(bodydata->qm[f]->vType[j] == -1) break;
          bodydata->qm[f]->degen[j] = 0;
          stat  = EG_evaluate(bodydata->qm[f]->face, &bodydata->qm[f]->uvs[2*j],
                              xyz);
          if (stat != EGADS_SUCCESS)
              printf(" EG_createMeshMap :: evaluate = %d !! \n", stat);
          norm1 = xyz[3] * xyz[3] + xyz[4] * xyz[4] + xyz[5] * xyz[5];
          norm2 = xyz[6] * xyz[6] + xyz[7] * xyz[7] + xyz[8] * xyz[8];
          if (sqrt(norm1) < qEPS || sqrt(norm2) < qEPS) {
#ifdef DEBUG
              printf("\n\n FACE %d V %d is degenerate\n", f + 1, j + 1);
              printf(" DU %lf %lf %lf DV %lf %lf %lf\n",
                     xyz[3], xyz[4], xyz[5], xyz[6], xyz[7], xyz[8]);
#endif
              bodydata->qm[f]->degen[j] = 1;
          }
      }
      for (j = 0; j < nquad; j++) {
          for (k = 0; k < 4; ++k)
            bodydata->qm[f]->qIdx[4*j + k] = tris[6*j + qV[k+1]];
      }
      for (j = 0; j < nquad; j++)
        {
          kk  = 0;
          kOK = 0;
          q   = 0;
          while (q < nquad) {
              if (q == j) {
                  if (j == nquad - 1 && kk < 4) {
                      bodydata->qm[f]->qAdj[4*j + kk++] = -1;
                      if (kk == 4) break;
                      q = 0;
                  }
                  else q++;
              }
              if  (q == nquad) break;
              for (k = 0; k < 4; ++k) {
                  if ((bodydata->qm[f]->qIdx[4*j + qLoop[kk    ]] ==
                      bodydata->qm[f]->qIdx[4*q + qLoop[k     ]] ||
                      bodydata->qm[f]->qIdx[4*j + qLoop[kk    ]] ==
                          bodydata->qm[f]->qIdx[4*q + qLoop[k  + 1]]) &&
                      (bodydata->qm[f]->qIdx[4*j + qLoop[kk + 1]] ==
                          bodydata->qm[f]->qIdx[4*q + qLoop[k     ]] ||
                          bodydata->qm[f]->qIdx[4*j + qLoop[kk + 1]] ==
                              bodydata->qm[f]->qIdx[4*q + qLoop[k  + 1]])) {
                      bodydata->qm[f]->qAdj[4*j + kk++] = q + 1;
                      q   = -1;
                      kOK =  1;
                      k   =  4;
                      if (kk == 4) q = nquad;
                  }
              }
              if ((kOK == 0) && (q >= nquad -1)) {
                  bodydata->qm[f]->qAdj[4*j + kk++] = -1;
                  q = -1;
                  if (kk == 4) break;
              }
              else kOK = 0;
              q++;
          }
          if (kOK == 0 && kk < 4) {
              while (kk < 4)
                bodydata->qm[f]->qAdj[4*j + kk++] = -1;
          }
        }
      for (j = 0; j < nquad; j++) {
          for (q = 0; q < 4; ++q)
            bodydata->qm[f]->valence[bodydata->qm[f]->qIdx[4 * j + q] - 1][0] = j + 1;
      }
      for (j = 0; j < len; j++) {
          stat = EG_setValence (bodydata->qm[f], j + 1);
          if (stat != EGADS_SUCCESS) {
              printf("In EG_createMeshMap :: set valence at %d is %d!!\n ",
                     j + 1, stat);
              bodydata->qm[f]->fID = 0;
              break;
          }
      }
      for (e4[0]= e4[1] = j = 0; j < len; j++) {
          if (bodydata->qm[f]->vType[j] == -1) break;
          stat     = EG_angAtBdVert(bodydata->qm[f], j + 1, e4, &angle);
          if (stat != EGADS_SUCCESS || angle < EPS08 ) {
#ifdef REPORT
              printf(" FACE %d EG_angAtBdVert %d angle %f\n ",f + 1, stat, angle);
              printf(" Vertices: %d %d %d \n ", j+ 1, e4[0], e4[1]);
#endif
              if (bodydata->qm[f]->vType[j] == 0) bodydata->qm[f]->vType[j] = 2;
          }
          else if (angle < 0.75 * PI) bodydata->qm[f]->vType[j] = 2;
          else if (angle < 1.25 * PI) bodydata->qm[f]->vType[j] = 3;
          else if (angle < 1.75 * PI) bodydata->qm[f]->vType[j] = 4;
          else                        bodydata->qm[f]->vType[j] = 5;
		  if (bodydata->qm[f]->valence[j][2] < bodydata->qm[f]->vType[j])
              bodydata->qm[f]->vType[j] = bodydata->qm[f]->valence[j][2];
		  if (bodydata->qm[f]->vType[j] == 2)
                   bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] + 2;
          else if (bodydata->qm[f]->vType[j] == 3)
                   bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] + 1;
          else if (bodydata->qm[f]->vType[j] >= 5)
                   bodydata->qm[f]->valence[j][1] = bodydata->qm[f]->valence[j][2] - 1;
          bodydata->qm[f]->bdAng[j] = angle;
      }
  }
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_makeValidMesh(meshMap *qm, int nP, /*@null@*/ int *pList, int fullReg)
{
  int    si, v, q, i, ii, j, k, kv, it = 0, itMax, sum = 0, *qlist = NULL, kq = 0;
  int    stat = EGADS_SUCCESS, *mv = NULL, *area = NULL, fr, pass;
  double *uvxyz = NULL;
#ifdef DEBUG2
  double pos[18], uv[2];
  FILE   *fout = NULL;
  char   buffer[100];

  gnuData(qm, NULL);
#endif

  mv      = (int   *)EG_alloc(qm->totV * sizeof(int));
  qlist   = (int   *)EG_alloc(qm->totQ * sizeof(int));
  if (mv == NULL || qlist == NULL) return EGADS_MALLOC;
  if (fullReg == 0) { // move around only affected vertices
      if (nP  == 0 || pList == NULL) goto cleanup;
      for (kv = j = 0; j < nP; j++) {
          si = pList[j] - 1;
          if (qm->vType[si] == -2) continue;
          if (qm->vType[si] == -1 &&
              inList (kv, mv, si) == -1 ) mv[kv++] = si;
          if( qm->star[si] == NULL) {
              printf(" STAR %d is NULL !!!\n ", si + 1);
              stat = EGADS_MALLOC;
              goto cleanup;
          }
          for (i = 0; i < qm->star[si]->nQ; i++) {
              if (qm->star[si]->quads[i] == -1) continue;
              v = qm->star[si]->verts[2 * i + 1] - 1;
              if (qm->vType[v] == -1 && inList (kv, mv, v) == -1 ) mv[kv++] = v;
              v = qm->star[si]->verts[2 * i + 2] - 1;
              if (qm->vType[v] == -1 && inList (kv, mv, v) == -1 ) mv[kv++] = v;
          }
      }
      if ( kv == 0 ) {
          EG_free(mv);
          EG_free(qlist);
          return EGADS_SUCCESS;
      }
      itMax = 5;
      uvxyz = (double*)EG_alloc(5 * kv * sizeof(double));
      if ( uvxyz == NULL ) {
          EG_free(mv);
          EG_free(qlist);
          return EGADS_MALLOC;
      }
      for (j = 0; j < kv; j++) {
          uvxyz[5 * j    ] = qm->uvs [2 * mv[j]    ];
          uvxyz[5 * j + 1] = qm->uvs [2 * mv[j] + 1];
          uvxyz[5 * j + 2] = qm->xyzs[3 * mv[j]    ];
          uvxyz[5 * j + 3] = qm->xyzs[3 * mv[j] + 1];
          uvxyz[5 * j + 4] = qm->xyzs[3 * mv[j] + 2];
      }
      for ( kq = k = 0 ; k < kv; k++) {
          for (i = 0; i < qm->star[mv[k]]->nQ; i++) {
              if (kq == 0 || inList(kq, qlist, qm->star[mv[k]]->quads[i]) == -1)
                qlist[kq++] = qm->star[mv[k]]->quads[i];
          }
      }
  } else {
      for (kv = i = 0 ; i < qm->totV; i++ )
          if ( qm->vType[i] == -1) mv[kv++] = i;
      for (kq = i = 0 ; i < qm->totQ; i++ )
          if (qm->vType[qm->qIdx[4*i]-1] != -2) qlist[kq++] = i + 1;
      itMax = 40;
  }
  area = (int *)EG_alloc(kv * sizeof(int));
  if (area == NULL ) {
      EG_free(mv);
      EG_free(uvxyz);
      EG_free(qlist);
      return EGADS_MALLOC;
  }
  for (i = 0; i < kv; i++) {
      if (qm->star[mv[i]] == NULL ) {
          printf(" Star in v %d is nULL \n ", i + 1 );
          stat = EGADS_MALLOC;
          goto cleanup;
      }
  }
#ifdef DEBUG2
  printf(" EG_makeValid mesh for %d regular  ---- \n ", fullReg);
  snprintf(buffer, 100,"MKSA_%d", qm->plotcount++);
  printf("Writing in %s\n ", buffer);
  fout = fopen(buffer,"w");
  if (fout != NULL ) {
      for (i = 0 ; i < kq; i++ ) {
          uv[0] = uv[1] = 0.0;
          for (j = 0; j <= 4; j++) {
              v = qm->qIdx[ 4 * (qlist[i] - 1) + j%4 ] - 1;
              fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*v  ],
              qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1 );
              if ( j == 4 ) break;
              uv[0] += 0.25 * qm->uvs[2 * v    ];
              uv[1] += 0.25 * qm->uvs[2 * v + 1];
          }
          fprintf(fout,"\n\n");
          EG_evaluate(qm->face, uv, pos);
          fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
          qlist[i] );
          fprintf(fout,"\n\n");
      }
      fclose (fout);
  }
#endif
#ifndef __clang_analyzer__
  pass = -1;
#endif
  fr   = fullReg;
  for (k = 0 ; k < kv; k++) {
      qm->star[mv[k]]->type = -1;
      for (i = 0; i < qm->star[mv[k]]->nV; ++i) {
          si = qm->star[mv[k]]->verts[i] - 1;
          if (qm->vType[si] >= 0) {
              if (qm->bdAng[si] > PIEPS) qm->star[mv[k]]->type = 3;
              else                       qm->star[mv[k]]->type = MAX (qm->star[mv[k]]->type, 2);
          } else {
              for (ii = 0 ; ii < qm->star[si]->nQ; ii++) {
                  v   = qm->star[si]->verts[2 * ii + 1] - 1;
                  if (qm->vType[v] >= 0 ) {
                      if (qm->bdAng[v] > PIEPS) qm->star[mv[k]]->type = MAX (qm->star[mv[k]]->type, 1);
                      else                      qm->star[mv[k]]->type = MAX (qm->star[mv[k]]->type, 0);
                  }
              }
          }
      }
  }
  for (it = 0 ; it < itMax; it++) {
      if (it % 2 == 0 && fullReg == 1) pass = -1;
	  else pass = 1;
      if (fullReg == 1 && it == itMax - 3) fr = 2;
      if (it == itMax -1 ) fr++;
      for (sum = k = 0 ; k < kv; k++) {
          if (qm->star[mv[k]] ->type <= 0) continue;
          area[k] = EG_placeVertex(qm, mv[k] + 1, pass, fr);
          sum     = MAX(sum,area[k]);
#ifdef DEBUG2
          if (area[k] >= QA1) printf("VERTEX %d has area %d\n",
                                     mv[k] + 1, area[k]);
#endif
      }
      for (k = 0 ; k < kv; k++) {
          if (qm->star[mv[k]] ->type > 0 ) continue;
          area[k] = EG_placeVertex(qm, mv[k] + 1, pass, fr);
           sum    = MAX(sum,area[k]);
#ifdef DEBUG2
          if (area[k] >= QA1) printf("VERTEX %d has area %d\n",
                                     mv[k] + 1, area[k]);
#endif
      }
#ifdef DEBUG2
      printf(" MAKEVALID IT %d ===================  SUM %d  \n ", it, sum );
      snprintf(buffer, 100,"MKLA_%d", qm->plotcount++);
      printf("Writing in %s\n ", buffer);
      fout       = fopen(buffer,"w");
      if (fout != NULL ) {
          for (i = 0 ; i < kq; i++ ) {
              uv[0] = uv[1] = 0.0;
              for (j = 0; j <= 4; j++) {
                  v = qm->qIdx[ 4 * (qlist[i] - 1) + j%4 ] - 1;
                  fprintf(fout, "%lf %lf %lf %d\n",  qm->xyzs[3*v  ],
                  qm->xyzs[3*v + 1], qm->xyzs[3*v + 2], v + 1 );
                  if ( j == 4 ) break;
                  uv[0] += 0.25 * qm->uvs[2 * v    ];
                  uv[1] += 0.25 * qm->uvs[2 * v + 1];
              }
              fprintf(fout,"\n\n");
              EG_evaluate(qm->face, uv, pos);
              fprintf(fout, "%lf %lf %lf %d\n", pos[0], pos[1], pos[2],
              qlist[i] );
              fprintf(fout,"\n\n");
          }
          fclose (fout);
      }
      printf(" MAKEVALID ROUND %d CHECK ALL QUADS \n", it);
      if (qm->vInv) {
          for ( i = 0 ; i < qm->vInv[0]; i++ )
              printf(" INVALID VERTEX %d --> %d\n", i + 1, qm->vInv[i+1]);
      }
#endif
      if (sum == QA0 && fullReg == 0) break;
  }
  //if (sum  < QA2 && fullReg == 0) sum = QA0;
  if (sum  < QA2) sum = QA0;  // allowing area 100, large angles but still valid. The normal plane angle function is not 100% reliable
  if (sum != QA0) {
      if (fullReg == 0) {
          stat = EGADS_GEOMERR;
#ifdef DEBUG2
          if (sum < QA2)
              printf(" EG MAKE VALID Mesh has obtuse angles\n");
          else
              printf("EG MAKE VALID  Mesh has invalid quads\n");
          snprintf(buffer, 100,"face_%d_InvalidMesh_%d",
          qm->fID, qm->invsteps);
          gnuData(qm , buffer);
#endif
          qm->invsteps++;
          if (uvxyz) {
              for (j = 0; j < kv; j++) {
                  if (qm->vType[mv[j]] != -1) continue;
                  qm->uvs  [2 * mv[j]    ] = uvxyz[5 * j    ];
                  qm->uvs  [2 * mv[j] + 1] = uvxyz[5 * j + 1];
                  qm->xyzs [3 * mv[j]    ] = uvxyz[5 * j + 2];
                  qm->xyzs [3 * mv[j] + 1] = uvxyz[5 * j + 3];
                  qm->xyzs [3 * mv[j] + 2] = uvxyz[5 * j + 4];
              }
          }
      }
      else stat = EGADS_GEOMERR;
  }
  for (sum = q = 0; q < kv; q++) {
      k   = -1;
      if (qm->vInv) k = inList(qm->vInv[0], &qm->vInv[1], mv[q] + 1);
      j = area[q];
      if ((fullReg > 0 && j != QA0) ||
      (k == -1     && ((qm->pp == 0 && j >= QA2) ||
      (qm->pp == 1 && j >  QA2)))) {
#ifdef DEBUG
          printf(" Q = %d ---> AREA %d   SUM %d \n", mv[q] + 1, area[q], sum);
#endif
          sum++;
      }
  }
  if (qm->vInv != NULL) {
      for (k = 0 ; k < qm->vInv[0]; k++)  {
          if (qm->vInv[1 + k] <= 0 || qm->vInv[1 + k] > qm->totV ||
              qm->star [qm->vInv[1 + k] -1] == NULL ||
              qm->vType[qm->vInv[1 + k] -1] == -2) {
              qm->vInv[1 + k] = -2;
              continue;
          }
          if (EG_vertexArea(qm, qm->vInv[1 + k]) < QA2)
              qm->vInv[1 + k] = -2; // vInv has all the initial invalid quads.If fixed, remove from list
      }
  }
#ifdef DEBUG2
  printf(" EG_makeValidMesh stat %d\n", stat);
  gnuData(qm, NULL);
#ifndef __CUDA_ARCH__
  if (qm->plotcount > 50000) fprintf(stderr," SOMHI !\n");
#else
  if (qm->plotcount > 50000) printf(" SOMHI !\n");
#endif
#endif
#ifdef REPORT
  if ( fullReg > 0 )
    printf(" Face %d -->FINAL MESH %d\n", qm->fID, stat);
#endif
  cleanup:
  EG_free(uvxyz);
  EG_free(mv);
  EG_free(area);
  EG_free(qlist);
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_createQuadGroup(meshMap *qm, quadGroup *qg, int q0, int q1)
{
  int i, j, k, ids[8], piv = 0, aux, vaux[6], vc[4], c1 = 0, c2 = 0;

  if (q0 < 0 || q1 < 0 || q0 > qm->totQ || q1 > qm->totQ ||
      qm->qIdx[4 * ( q0 - 1 )] == -2 ||
      qm->qIdx[4 * ( q1 - 1 )] == -2 ) {
#ifdef DEBUG
      printf(" EG_commonVerts wrong quad pair %d %d\n", q0, q1);
      if (q0 > 0)
        printQuad(qm, q0);
      if (q1 > 0)
        printQuad(qm, q1);
#endif
      return EGADS_INDEXERR;

  }
  qg->q[0] = q0;
  qg->q[1] = q1;
  if (q0 == q1) {
      printf("EG_createQuadGroup :: Quads %d are the same!!!\n", q0);
      return EGADS_INDEXERR;
  }
  for (k = i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
          if (qm->qIdx[4 * (q0 - 1) + i] == qm->qIdx[4 * (q1 - 1) + j] &&
              qm->qIdx[4 * (q0 - 1) + i] != -1) {
              ids[2 *k    ] = i;
              ids[2 *k + 1] = j;
              vc [   k++  ] = qm->qIdx[4 * (q1 - 1) + j];
          }
      }
      if ( k == 2) break;
  }
  if (k < 2) {
      printf("EG_createQuadGroup: Quad %d %d are not adjacent!!\n", q0, q1);
#ifdef DEBUG
      printQuad(qm, q0);
      printQuad(qm, q1);
#endif
      return EGADS_INDEXERR;
  }
  if (qm->valence[vc[1] -1][1] >
  qm->valence[vc[0] -1][1]) {
      ids[0] = ids[2];
      ids[1] = ids[3];
  }
  piv    = 0;
  if (qm->qAdj[4 * (qg->q[0] - 1) + ids[0]] == qg->q[1]) piv = 1;
  qg->verts[0] = qm->qIdx[4 * (qg->q[piv] - 1) +  ids[piv]];
  qg->verts[1] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 1)%4];
  qg->verts[2] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 2)%4];
  qg->verts[3] = qm->qIdx[4 * (qg->q[piv] - 1) + (ids[piv] + 3)%4];
  aux = (piv + 1)%2;
  qg->verts[4] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 2)%4];
  qg->verts[5] = qm->qIdx[4 * (qg->q[aux] - 1) + (ids[aux] + 3)%4];
  if (piv == 1) {
      aux      = qg->q[0];
      qg->q[0] = qg->q[1];
      qg->q[1] = aux;
  }
  if (qm->valence[qg->verts[0] -1][1] ==
      qm->valence[qg->verts[3] -1][1]) {
      if (qm->valence[qg->verts[1] -1][1] != 4) c1++;
      if (qm->valence[qg->verts[5] -1][1] != 4) c1++;
      if (qm->valence[qg->verts[2] -1][1] != 4) c2++;
      if (qm->valence[qg->verts[4] -1][1] != 4) c2++;
      if (c2 > c1) {
          vaux[0]      = qg->verts[0];
          vaux[1]      = qg->verts[1];
          vaux[2]      = qg->verts[2];
          qg->verts[0] = qg->verts[3];
          qg->verts[1] = qg->verts[4];
          qg->verts[2] = qg->verts[5];
          qg->verts[3] = vaux[0];
          qg->verts[4] = vaux[1];
          qg->verts[5] = vaux[2];
          aux          = qg->q[0];
          qg->q[0]     = qg->q[1];
          qg->q[1]     = aux;
      }
  }
  qg->vals[0] = qm->valence[qg->verts[0] -1][1];
  qg->vals[1] = qm->valence[qg->verts[1] -1][1];
  qg->vals[2] = qm->valence[qg->verts[2] -1][1];
  qg->vals[3] = qm->valence[qg->verts[3] -1][1];
  qg->vals[4] = qm->valence[qg->verts[4] -1][1];
  qg->vals[5] = qm->valence[qg->verts[5] -1][1];
  return EGADS_SUCCESS;
}


/* Assuming qID is collapsing through v */
__HOST_AND_DEVICE__ static int EG_validCollapse(meshMap *qm, int qID, int v)
{
  int j, k, kk, id, link[2], aux, aux2, type, val;
  if (v <= 0 || v > qm->totV) return EGADS_INDEXERR;
  id     = EG_quadVertIdx(qm, qID, v);
  aux    = qm->qIdx[4 * (qID - 1) + (id + 2)%4] - 1;
  if (qm->vType[aux] != -1) {
      for (j = 0; j < qm->valence[v - 1][2]; j++) {
          k  = qm->valence[v - 1][3 + j] -1;
          if (qm->vType[k] == -1) continue;
          if (EG_quadVertIdx(qm, qID, k + 1) < 0) {
#ifdef DEBUG
              printf(" Vertex %d is bounds and not in quad %d. "
                  "This will produce a new artificial edge, ie boundary vertices from "
                  "different edges will be linked. Not at all desirable \n ",
                     k + 1, qID );
#endif
              return 0;
          }
      }
  }
  link[0]  = qm->qIdx[4 * (qID - 1) + (id + 1)%4];
  link[1]  = qm->qIdx[4 * (qID - 1) + (id + 3)%4];
  for (j   = 0; j < 2; j++) {
      type = qm->vType  [link[j] - 1];
      val  = qm->valence[link[j] - 1][2];
      if (type > 0 ) {
          if ( val <= type) {
#ifdef DEBUG
              printf(" link %d has valence %d type %d !\n ", link[j], val, type);
#endif
              return 0;
          }else continue;
      }
      else if (val != 3) continue;
      if (qm->vType[link[(j+1)%2] -1] >=4 ) return 0;
      for (k  = 0; k < 3; k++) {
          aux = qm->valence[link[j] - 1][3 + k] - 1;
          if (EG_quadVertIdx(qm, qID,  aux + 1) >= 0) continue;
          if (qm->valence[aux][2] == qm->vType[aux]) {
#ifdef DEBUG
              printf(" Chain Effect Doublets %d -->%d %d \n",
                     aux + 1, qm->valence[aux][2], qm->vType[aux]);
#endif
              return 0;
          }
          if (qm->valence[aux][2] != 3) continue;
          for (kk = 0; kk < 3; kk++) {
              aux2 = qm->valence[aux][3 + kk] - 1;
              if (aux2 + 1 == link[j]) continue;
              if (qm->vType[aux] >= 4 ||
                  qm->vType[aux] == 2 ||
                  qm->valence[aux2][2] == qm->vType[aux] ||
                  (qm->valence[aux2][2] == 3 && qm->vType[aux] == -1)) {
#ifdef DEBUG
                  printf(" KK Chain Effect Doublets %d -->%d %d \n",
                         aux2 + 1, qm->valence[aux][2], qm->vType[aux]);
#endif
                  return 0;
              }
          }
      }
  }
  return 1;
}


/* Assuming we will break link v1-v2 */
__HOST_AND_DEVICE__ static int EG_validSwap(meshMap *qm, int v1, int v2)
{
  int i, vs[2], type, val;

  vs[0]    = v1 - 1; vs[1] = v2 - 1;
  for (i   = 0; i < 2; i++) {
      type = qm->vType  [vs[i]];
      val  = qm->valence[vs[i]][2];
      if ((type > 0 && val <= type) || (type == -1 && val <= 3)) return 0;
  }
  return 1;
}


__HOST_AND_DEVICE__ static int
EG_swappingOperation(meshMap *qm, quadGroup qg, int swap, int *activity)
{
  int   nq, stat, i0, i1, i, j, k, adj;
  int   *list, qID[2], adjQmap[6];
  double  uv[4] ;
  Quad *quad = NULL;

  swap      = swap%3;
  *activity = 0;
  if (swap == 0) {
#ifdef DEBUG
      printf(" swapping throu 0-3 will result in the same pair!! \n ");
#endif
      return EGADS_INDEXERR;
  }
  qID[0]     = qg.q[0]; qID[1] = qg.q[1];
  nq         = 2;
  stat       = EG_backupQuads(qm, &nq, qID, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf(" In EG_swappingOperation :: EG_backupQuad %d !!\n ", stat);
      return stat;
  }
  i0 = EG_quadVertIdx(qm, qID[0], qg.verts[0]); // centre
  i1 = EG_quadVertIdx(qm, qID[1], qg.verts[3]); // opposite
  // Modify Quads and get adj map
  qm->qIdx[4 * (qID[0] - 1)] = qg.verts[swap];
  qm->qIdx[4 * (qID[1] - 1)] = qg.verts[swap];
  for (i = 0; i < 3; i++) {
      adjQmap[i    ] = qm->qAdj[4 * (qID[0] - 1) + (i + i0)%4];
      adjQmap[i + 3] = qm->qAdj[4 * (qID[1] - 1) + (i + i1)%4];
      qm->qIdx[4 * (qID[0] - 1) + i + 1] = qg.verts[(swap + 3 + i)%6];
      qm->qIdx[4 * (qID[1] - 1) + i + 1] = qg.verts[(swap + 1 + i)  ];
  }
  qm->qAdj[4 * (qID[0] - 1)    ] = qID[1];
  qm->qAdj[4 * (qID[1] - 1) + 3] = qID[0];
  for (i = 0; i < 3; i++) {
      adj = adjQmap[(3 + i + swap)%6] - 1;
      qm->qAdj[4 * (qID[0] - 1) + i + 1] = adj + 1;
      if (adj >= 0) {
          for (j = 0; j < 4; j++) if (qm->qAdj[4 * adj + j] == qID[1])
            qm->qAdj[4 * adj + j] = qID[0];
      }
      adj = adjQmap[i + swap] -1;
      qm->qAdj[4 * (qID[1] - 1) + i] = adj + 1;
      if (adj >= 0) {
          for (j = 0; j < 4; j++) if (qm->qAdj[4 * adj + j] == qID[0])
            qm->qAdj[4 * adj + j] = qID[1];
      }
  }
  for (i = 0; i < 4; i++) {
      j = qm->qIdx[4 * (qID[0] - 1) + i] - 1;
      qm->valence[j][0] = qID[0];
      j = qm->qIdx[4 * (qID[1] - 1) + i] - 1;
      qm->valence[j][0] = qID[1];
  }

  list = (int *) EG_alloc (qm->totV * sizeof(int));
  if (list == NULL) {
      EG_free(quad);
      return EGADS_MALLOC;
  }
  for (i0 = i = 0; i < 4; i++) {
      if (inList(i0, list, qg.verts[i] ) == -1) list[i0++] = qg.verts[i];
      for (j = 1; j < qm->star[qg.verts[i] - 1]->nV; j++) {
          k  = qm->star[qg.verts[i] - 1]->verts[j];
          if ( k == -1) continue;
          if (inList(i0, list, k ) == -1 ) list[i0++] = k;
      }
  }
  for (i = 0; i < i0; i++) {
      stat = EG_setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
#ifndef __CUDA_ARCH__
          fprintf(stderr," EG_swappingOperation :: ERROR SETTING VALENCE %d \n ",
                  list[i]);
#else
          printf(" EG_swappingOperation :: ERROR SETTING VALENCE %d \n ",
                 list[i]);
#endif
          EG_free(list);
          return stat;
      }
  }
  EG_free(list);
  *activity = 1;
  i     = EG_vertexArea(qm,qg.verts[0]);
  j     = EG_vertexArea(qm,qg.verts[3]);
  uv[0] = qm->uvs[2*(qg.verts[0]-1)    ];
  uv[1] = qm->uvs[2*(qg.verts[0]-1) + 1];
  uv[2] = qm->uvs[2*(qg.verts[3]-1)    ];
  uv[3] = qm->uvs[2*(qg.verts[3]-1) + 1];
  if (i != QACB && j != QACB) {
      if (i != QA0 || j != QA0) {
          if (qm->vType[qg.verts[0]-1] == -1 &&
              qm->vType[qm->valence[qg.verts[0]-1][3]-1] < 4) {
              EG_centroid(qm, qm->valence[qg.verts[0] -1][2],
                          &qm->valence[qg.verts[0] -1][3],
                          &qm->uvs[2 * (qg.verts[0] - 1)], 1);
              updateVertex(qm, qg.verts[0], &qm->uvs[2 * (qg.verts[0] - 1)] );
          }
          if (qm->vType[qg.verts[3]-1] == -1 &&
              qm->vType[qm->valence[qg.verts[3]-1][3]-1] < 4) {
              EG_centroid(qm, qm->valence[qg.verts[3] -1][2],
                          &qm->valence[qg.verts[3] -1][3],
                          &qm->uvs[2 * (qg.verts[3] - 1)], 1);
              updateVertex(qm, qg.verts[3], &qm->uvs[2 * (qg.verts[3] - 1)] );
          }
          i0 = EG_vertexArea(qm, qg.verts[0]);
          i1 = EG_vertexArea(qm, qg.verts[3]);
          if (i0 + i1 > i + j) {
              if(qm->vType[qg.verts[0]-1])
                updateVertex(qm, qg.verts[0], uv);
              if(qm->vType[qg.verts[3]-1])
                updateVertex(qm, qg.verts[3], &uv[2]);
          }
      }
      if (EG_makeValidMesh(qm, 6, qg.verts, 0) == EGADS_SUCCESS) {
          EG_free(quad);
          return EGADS_SUCCESS;
      }
  } else {
#ifdef DEBUG
    printf(" EG_swappingOperation : reject swap. It crosses the bounds %d  %d  \n ",
           i, j);
    gnuData(qm, NULL);
#endif
  }
#ifdef DEBUG
  printf(" EG_swapping operation didn't work well: restore quads\n ");
#endif
  if(qm->vType[qg.verts[0]-1])
    updateVertex(qm, qg.verts[0], uv);
  if(qm->vType[qg.verts[3]-1])
    updateVertex(qm, qg.verts[3], &uv[2]);
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_swappingOperation :: restore quads %d !!\n ", stat);
      return stat;
  }
  *activity = 0;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_splittingOperation(meshMap *qm, int vC, int vL, int vR, int *activity)
{
  int   qIdx[4], modQ[4], verts[4], adj[2], poly[4], q, newQ, i, j, stat;
  int   id0 = -1, id1 = -1, dist, links[4], vals[4], addedV = 0, nq, si, *list = NULL, n;
  double uv[4];

  Quad  *quad = NULL;

  si          = vC - 1;
  *activity   = 0;
  if (qm->star[si] == NULL ) {
      printf(" In splittingOperation build star %d is NULL !!\n", vC);
      return EGADS_MALLOC;
  }
  nq        = qm->star[si]->nQ;
  stat      = EG_backupQuads(qm, &nq, qm->star[si]->quads, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      EG_free(quad);
      printf("In EG_splittingOperation :: EG_backupQuads %d !!\n", stat);
      return stat;
  }
  if (qm->remQ[0] > 0) {
      poly[3] = qm->remV[qm->remV[0]--];
      newQ    = qm->remQ[qm->remQ[0]--];
  } else {
      poly[3] = qm->totV + 1;
      newQ    = qm->totQ + 1;
      if (poly[3] >  qm->sizeV) {
#ifdef DEBUG
          printf(" We have duplicated the number of initial vertices. This is too much. \n");
#endif
          EG_free(quad);
          return EGADS_INDEXERR;
      }
      ++qm->totV; ++qm->totQ;
      addedV = 1;
  }
  qm->vType[poly[3] - 1] = -1;
  if (qm->vType[vC - 1] != -1) {
      for (q = 0; q < qm->star[si]->nQ; q++)
        if (qm->star[si]->quads[q] == -1) break;
      for (i = 0; i < 2; i++) {
          if (i == 1) dist = (qm->star[si]->nV - 1) - 4;
          else        dist = 4;
          id0 = qm->star[si]->idxQ[q + i];
          links[2 * i    ] = qm->star[si]->verts[qm->star[si]->idxV[2 * id0 + 1       ]];
          links[2 * i + 1] = qm->star[si]->verts[qm->star[si]->idxV[2 * id0 + 1 + dist]];
          vals [2 * i    ] = qm->valence [links[2*i    ] - 1][1];
          vals [2 * i + 1] = qm->valence [links[2*i + 1] - 1][1];
      }
/*    dist = 4;
      i    = 0;    */
      if (vals[0] * vals[1] > vals[2] * vals[3]) {
          vL = links[3];
          vR = links[2];
      } else {
          vL = links[0];
          vR = links[1];
      }
  }
  id0 = - 1; id1 = -1;
  for (j = 0; j < qm->star[si]->nQ; j++) {
      if (qm->star[si]->verts[2 * j + 1] == vL) id0 = j;
      if (qm->star[si]->verts[2 * j + 1] == vR) id1 = j;
  }
  if (id0 == -1 || id1 == -1) {
#ifdef DEBUG
    printf("EG_splittingOperation  Problem finding vertices %d !!\n ",
           EGADS_INDEXERR);
    printf(" vL %d vR %d \n ", vL, vR);
    for (j = 0; j < qm->star[si]->nV; j++)
      printf(" j %d =%d\n ", j, qm->star[si]->verts[j]);
#endif
    EG_free(quad);
    return EGADS_INDEXERR;
  }
  // update stars
  poly [0] = qm->star[si]->verts[0];
  poly [1] = qm->star[si]->verts[2 * id0 + 1];
  poly [2] = qm->star[si]->verts[2 * id1 + 1];
  uv[2]    = qm->uvs[2 * (poly[0] -1)    ];
  uv[3]    = qm->uvs[2 * (poly[0] -1) + 1];
  qIdx [0] = id0;
  qIdx [1] = qm->star[si]->idxQ[id1 + qm->star[si]->nQ - 1];
  qIdx [2] = id1;
  qIdx [3] = qm->star[si]->idxQ[id0 + qm->star[si]->nQ - 1];
  verts[0] = poly[1];
  verts[1] = poly[2];
  verts[2] = poly[2];
  verts[3] = poly[1];
  qm->qIdx[4 * (newQ - 1)    ] = poly[1];
  qm->qIdx[4 * (newQ - 1) + 1] = poly[0];
  qm->qIdx[4 * (newQ - 1) + 2] = poly[2];
  qm->qIdx[4 * (newQ - 1) + 3] = poly[3];
  for (i = 0; i < 4; ++i) {
      modQ[i] = qm->star[si]->quads[qIdx[i]];
      qm->qAdj[4 * (newQ - 1) + i] = modQ[i];
      qm->valence[qm->qIdx[4*(newQ - 1)+i]-1][0] = newQ;
      if (modQ[i] == -1) continue;
      (void) EG_adjQtoPair(qm, modQ[i], poly[0], verts[i], adj);
      qm->qAdj[4 * (modQ[i] - 1) + adj[0]] = newQ;
  }
  j = qIdx[2];
  q = qm->star[si]->quads[j++];
  while (q != qm->star[si]->quads[qIdx[0]]) {
      for (i = 0; i < 4; ++i)
        if (qm->qIdx[4 * (q - 1) + i] == poly[0])
          qm->qIdx[4 * (q - 1) + i] = poly[3];
      q = qm->star[si]->quads[qm->star[si]->idxQ[j++]];
      if (q == -1) {
          if (qm->star[si]->quads[qIdx[0]] == -1) break;
          else q = qm->star[si]->quads[qm->star[si]->idxQ[j++]];
      }
  }
  list = (int *) EG_alloc(qm->totV * sizeof(int));
  if (list == NULL) {
      EG_free(quad);
      return EGADS_MALLOC;
  }
  EG_setValence(qm, poly[3]);
  for (n = i = 0; i < 4; i++) {
      if (inList(n, list, poly[i]) == -1) list[n++] = poly[i];
      for (j = 1; j < qm->star[poly[i] - 1]->nV; j++) {
          q  = qm->star[poly[i] - 1]->verts[j];
          if ( q == -1 ) continue;
          if (inList(n, list, q ) == -1 ) list[n++] = q;
      }
  }
  for (i = 0 ; i < n; i++) {
      stat = EG_setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
#ifndef __CUDA_ARCH__
          fprintf(stderr," EG_splittingOperation :: ERROR SETTING VALENCE %d \n ",
                  poly[i]);
#else
          printf(" EG_splittingOperation :: ERROR SETTING VALENCE %d \n ",
                 poly[i]);
#endif
          EG_free(quad);
          EG_free(list);
          return stat;
      }
  }
  EG_free(list);
  si = poly[3] - 1;
  if (qm->star[si] == NULL) {
#ifndef __CUDA_ARCH__
      fprintf(stderr," EG_splittingOperation star for new vertex %d is NULL\n",
              poly[3]);
#else
      printf(" EG_splittingOperation star for new vertex %d is NULL\n",
             poly[3]);
#endif
      return EGADS_MALLOC;
  }
  EG_centroid(qm, qm->star[si]->nQ, &qm->valence[si][3], uv, 0);
  updateVertex(qm, poly[3], uv);
  if (qm->vType[poly[0]-1] == -1) {
      EG_centroid(qm, qm->star[poly[0]-1]->nQ, &qm->valence[poly[0]-1][3], uv, 0);
      updateVertex(qm, poly[0], uv);
  }
#ifdef DEBUG
  printf(" SPLIT VERTEX %d INTRO %d\n", poly[0], poly[3]);
  gnuData(qm, NULL);
#endif
  id1 = 0;
  id0 = EG_placeVertex(qm, poly[3], 1, 0);
  if (qm->vType[poly[0]-1] == -1)
    id1 = EG_placeVertex(qm, poly[0], 1, 0);
  if (id0 >= QA2) {
      if (id0 >= QA3 ) id0 = 1;
      else {
          for ( id0 = i = 0 ; i < qm->star[poly[3] -1]->nQ; i++)
            if (qm->vType[qm->star[poly[3] -1]->verts[2 * i + 1]-1] >= 4) id0 = 1;
      }
  } else id0 = 0;
  if (id1 >= QA2) {
      if (id1 >= QA3 ) id1 = 1;
      else {
          for (id1 = i = 0 ; i < qm->star[poly[0] -1]->nQ; i++)
            if (qm->vType[qm->star[poly[0] -1]->verts[2 * i + 1]-1] >= 4) id1 = 1;
      }
  } else id1 = 0;
  if (id0 == 0 && id1 == 0) {
      *activity = 1;
#ifdef DEBUG
    printf(" --------SPLIT THRU  AREA %d  %d  %d  MAX %d\n",
           id0, id1, poly[0], QA2);
    gnuData(qm, NULL);
#endif
    if (EG_makeValidMesh(qm, 4, poly, 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
    }
  }

#ifdef DEBUG
  printf("EG_splittingOperation mesh valid failed: restore quads\n");
  gnuData(qm, NULL);
#endif
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_splittingOperation :: EG_restoreQuads %d !!\n ", stat);
      return stat;
  }
  if (addedV == 1) {
      --qm->totV;
      --qm->totQ;
  } else {
      qm->remV[++qm->remV[0]] = poly[3];
      qm->remQ[++qm->remQ[0]] = newQ;
  }
  qm->vType[poly[3] - 1] = -2;
  for (i = 0; i < 4; i++)
    qm->qIdx[4 * (newQ - 1) + i] = -2;
  *activity = 0;
  if (qm->vType[poly[0] -1] == -1 ) {
      qm->uvs[2 * (poly[0] -1)   ] = uv[2] ;
      qm->uvs[2 * (poly[0] -1) +1] = uv[3] ;
      updateVertex(qm, poly[0], &uv[2]);
  }
#ifdef DEBUG
  gnuData(qm, NULL);
#endif
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_mergeVertices(meshMap *qm, int qC, int centre, int *activity)
{
  int    stat, i, j, q, adjq, adjPair[2], auxQ, oldQ[8], nq, doublet = 0;
  int    piv[4] = {1, 0, 3, 2}, n, *list = NULL;
  double uv[2], uvxyz[10];
  Quad   *quad = NULL;

  nq        = 1;
  stat      = EG_backupQuads(qm, &nq, &qC, &quad);
  if (stat != EGADS_SUCCESS || quad == NULL) {
      printf(" In EG_mergeVertices stat in EG_backpQuads %d !!\n ", stat);
      return stat;
  }
  EG_centroid(qm, 4, &qm->qIdx[4 * (qC - 1)], uv, 0);
  *activity  = 0;
  i          = EG_quadVertIdx(qm, qC, centre);
  if ( qm->vType[centre - 1] != -1 ) {
      i = (i+ 2)%4;
      centre = qm->qIdx[4 * (qC - 1) + i];
  }
#ifdef DEBUG
  printf(" EG_mergeVertices collapsing %d through %d = %d \n ",
         qC, centre, i);
  gnuData(qm, NULL);
#endif

  if (qm->vType[centre - 1] != -1 || i < 0) {
      printf(" EG_mergeVertices collapse through %d is type %d and id in quad %d->%d !!\n",
             centre, qm->vType[centre - 1], qC, i);
#ifdef DEBUG
      printQuad(qm, qC);
#endif
      return EGADS_INDEXERR;
  }

  if (qm->valence[centre-1][2] == 2) doublet = 1;
  for (q = 0; q < 4; ++q) {
      oldQ[q    ] = qm->qIdx[4 * (qC - 1) + (q + i)%4];
      oldQ[q + 4] = qm->qAdj[4 * (qC - 1) + (q + i)%4];
  }
  if (qm->vType[oldQ[0]-1] == -1 && qm->vType[oldQ[2]-1] == -1 &&
      qm->vType[qm->valence[oldQ[0]-1][3] - 1] > 0 &&
      qm->bdAng[qm->valence[oldQ[0]-1][3] - 1] > PIEPS) {
      i = (i + 2)%4;
/*    centre = qm->qIdx[4 * (qC - 1) + i];    */
      for (q = 0; q < 4; ++q) {
          oldQ[q    ] = qm->qIdx[4 * (qC - 1) + (q + i)%4];
          oldQ[q + 4] = qm->qAdj[4 * (qC - 1) + (q + i)%4];
      }
  }
  uvxyz[0] = qm->uvs [2 * (oldQ[0] - 1)    ];
  uvxyz[1] = qm->uvs [2 * (oldQ[0] - 1) + 1];
  uvxyz[2] = qm->xyzs[3 * (oldQ[0] - 1)    ];
  uvxyz[3] = qm->xyzs[3 * (oldQ[0] - 1) + 1];
  uvxyz[4] = qm->xyzs[3 * (oldQ[0] - 1) + 2];
  uvxyz[5] = qm->uvs [2 * (oldQ[2] - 1)    ];
  uvxyz[6] = qm->uvs [2 * (oldQ[2] - 1) + 1];
  uvxyz[7] = qm->xyzs[3 * (oldQ[2] - 1)    ];
  uvxyz[8] = qm->xyzs[3 * (oldQ[2] - 1) + 1];
  uvxyz[9] = qm->xyzs[3 * (oldQ[2] - 1) + 2];
  // save old quads vertices and adjacents to update map correctly
  for (i = 0; i < 4; i++) {
      q = oldQ[4 + i];
      if (q < 0) continue;
      adjq = oldQ[4 + piv[i]]; // - 1 bias
      stat = EG_adjQtoPair(qm, q, oldQ[i], oldQ[(i + 1)%4], adjPair);
      if (stat != EGADS_SUCCESS || adjPair[1] != qC) {
          printf("EG_mergeVertices stat %d  adj %d to pair %d %d is %d \n ",
                 stat, q ,oldQ[i], oldQ[(i + 1)%4], adjPair[1]);
          return  EGADS_INDEXERR;
      }
      qm->qAdj[4 * (q - 1) + adjPair[0]] = adjq;
  }
  // Eliminate vertex oldQ[0] from all the quads and its valences
  for (i = 0; i < qm->star[oldQ[0] - 1]->nQ; ++i) {
      q = qm->star[oldQ[0] - 1]->quads[i];
      if      (q == -1) continue; // ghost quad
      else if (q == qC) {  // eliminate quad qC
          for (j = 0; j < 4; ++ j) {
              qm->qIdx[4 * (qC - 1) + j] = -2;
              qm->qAdj[4 * (qC - 1) + j] = -2;
          }
      } else {
          for (j = 0; j < 4; ++ j)
            if (qm->qIdx[4 * (q - 1) + j] == oldQ[0])
              qm->qIdx[4 * (q - 1) + j] = oldQ[2];
      }
  }
  // Point all the collapsed quad vertices to a valid quad
  for (i = 1; i < 4; i++) {
      if (qm->valence[oldQ[i] - 1][0] != qC) continue;
      for (q = 0; q < 4; q++) {
          auxQ = oldQ[4 + q];
          if (auxQ < 0) continue;
          if (EG_quadVertIdx(qm, auxQ, oldQ[i]) >= 0) {
              qm->valence[oldQ[i] - 1][0] = auxQ;
              break;
          }
      }
      if (qm->valence[oldQ[i] - 1][0] == qC) {
          printf(" EG_mergeVertices opposite vertex still points at collapsed quad %d !!\n",
                 EGADS_INDEXERR);
          EG_free(list);
          return EGADS_INDEXERR;
      }
  }
  EG_setValence(qm, oldQ[2]);
  EG_setValence(qm, oldQ[1]);
  EG_setValence(qm, oldQ[3]);
  list = (int *) EG_alloc(qm->totV * sizeof(int));
  if (list == NULL) {
      EG_free(quad);
      return EGADS_MALLOC;
  }
  n = 0;
  for (i = 1; i < 4; i++) {
      for (j = 1; j < qm->star[oldQ[i] - 1]->nV; j++) {
          q  = qm->star[oldQ[i] - 1]->verts[j];
          if ( q == -1 ) continue;
          if (inList(n, list, q ) == -1 ) list[n++] = q;
      }
  }
  for ( i = 0 ; i < n; i++ ) {
      stat = EG_setValence(qm, list[i]);
      if (stat != EGADS_SUCCESS) {
          printf(" EG_mergeVertices error setting valence %d \n ", list[i]);
          EG_free(quad);
          EG_free(list);
          return stat;
      }
  }
  EG_free(list);
  if (qm->vType[oldQ[3] - 1] >= 4 ||
      qm->vType[oldQ[1] - 1] >= 4 ) doublet = 0;
  if (qm->vType[oldQ[2] - 1] == -1 && doublet == 0 &&
      qm->vType[qm->valence[oldQ[2]-1][3] - 1] < 4 ) {
      updateVertex(qm, oldQ[2], uv);
      EG_placeVertex(qm, oldQ[2], 1,  0);
  }
  qm->vType[oldQ[0] - 1]  = -2; // -2 = removed
  qm->remQ[++qm->remQ[0]] = qC;
  qm->remV[++qm->remV[0]] = oldQ[0];
  *activity  = 1;
#ifdef DEBUG
  printf(" BEFORE CALLING MAKEVALID \n ");
  gnuData(qm, NULL);
#endif

  if (doublet == 1) {
#ifdef DEBUG
      printf(" COLLAPSING DOUBLET\n ");
#endif
      EG_free(quad);
      if (qm->vType[oldQ[1] -1] == -1)
        EG_placeVertex(qm, oldQ[1],1, 0);
      if (qm->vType[oldQ[3] -1] == -1)
        EG_placeVertex(qm, oldQ[3], 1, 0);
      return EGADS_SUCCESS;
  }
  if (EG_makeValidMesh(qm, 3, &oldQ[1], 0) == EGADS_SUCCESS) {
      EG_free(quad);
      return EGADS_SUCCESS;
  }
#ifdef DEBUG
  printf("EG_mergeVertices collapse quad %d thru vertex %d -->%d !!\n",
         qC, oldQ[0], stat);
#endif
  qm->vType[oldQ[0] - 1]  = -1;
  stat = EG_restoreQuads(qm, quad, nq);
  EG_free(quad);
  if (stat != EGADS_SUCCESS) {
      printf("EG_mergeVertices :: EG_restoreQuads stat %d \n", stat);
      return stat;
  }
  i = -1;
  if (qm->remQ[0] > 0)
    i = inList(qm->remQ[0], &qm->remQ[1], qC);
  if (i != -1) {
      for (j = i + 1; j < qm->remQ[0] - 1; j++) {
          qm->remQ[j] = qm->remQ[j + 1];
          qm->remV[j] = qm->remV[j + 1];
      }
  }
  --qm->remQ[0];
  --qm->remV[0];
  if (qm->vType[oldQ[0]-1] == -1 ) {
      qm->uvs [2 * (oldQ[0] - 1)    ] = uvxyz[0];
      qm->uvs [2 * (oldQ[0] - 1) + 1] = uvxyz[1];
      qm->xyzs[3 * (oldQ[0] - 1)    ] = uvxyz[2];
      qm->xyzs[3 * (oldQ[0] - 1) + 1] = uvxyz[3];
      qm->xyzs[3 * (oldQ[0] - 1) + 2] = uvxyz[4];
  }
  if (qm->vType[oldQ[2]-1] == -1 ) {
      qm->uvs [2 * (oldQ[2] - 1)    ] = uvxyz[5];
      qm->uvs [2 * (oldQ[2] - 1) + 1] = uvxyz[6];
      qm->xyzs[3 * (oldQ[2] - 1)    ] = uvxyz[7];
      qm->xyzs[3 * (oldQ[2] - 1) + 1] = uvxyz[8];
      qm->xyzs[3 * (oldQ[2] - 1) + 2] = uvxyz[9];
  }
  *activity = 0;
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int EG_swap(meshMap *qm, int qIn, int *activity)
{
  int       stat, q, swap = 0, bq1 = 0, bq2 = 0, corner = 0;
  quadGroup qg;
  *activity  = 0;
  if (qm->qIdx[4 * (qIn - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_swap for quad %d is void \n ", qIn );
#endif
      return EGADS_SUCCESS;
  }
  if (qm->regBd == 1 && EG_regBoundQuad(qm, qIn) == 0) bq1 = 1;
  for (swap = q = 0; q < 4; q++) {
      qg.q[0]        = qIn;
      qg.q[1]        = qm->qAdj[4 * (qIn - 1) + q];
      if (qg.q[1]   == -1) continue;
      if (qm->regBd ==  1 && EG_regBoundQuad(qm, qg.q[1]) == 0) bq2 = 1;
      stat           = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat      != EGADS_SUCCESS) {
#ifdef DEBUG
          printf(" EG_swap->EG_createQuadGroup  %d !! \n", stat);
          printQuadGroup(qm, qg);
#endif
          return stat;
      }
#ifdef DEBUG
      printf(" SWAP GROUP %d %d\n", qg.q[0], qg.q[1]);
      printQuadGroup(qm, qg);
#endif


      if ( EG_validSwap(qm, qg.verts[0], qg.verts[3]) != 1) continue;
      corner = 0;
      if ((bq1 == 1 || bq2 == 1) &&
          ((qm->vType  [qg.verts[0] -1]    > 0 &&
              qm->valence[qg.verts[0] -1][1] > 4 ) ||
              (qm->vType  [qg.verts[3] -1]    > 0 &&
                  qm->valence[qg.verts[3] -1][1] > 4 ))) {
#ifdef DEBUG
          printQuadGroup (qm, qg);
#endif

          if (qm->vType [qg.verts[0] -1] == 2 ||
              qm->vType [qg.verts[3] -1] == 2 ) {
              corner = 1;
              if (qg.vals[1] * qg.vals[4] <=
                  qg.vals[2] * qg.vals[5] ) swap = 2;
              else swap = 1;
          } else {
              if      (qm->vType [qg.verts[1] -1] *
                  qm->vType [qg.verts[4] -1] == 1 ) swap = 1;
              else if (qm->vType [qg.verts[2] -1] *
                  qm->vType [qg.verts[5] -1] == 1 ) swap = 2;
          }
      } else {
          if      (qg.vals[0] <= 4 ) continue;
          if      (qg.vals[1] * qg.vals[4] == 9) swap = 1;
          else if (qg.vals[2] * qg.vals[5] == 9) swap = 2;
          else if (qg.vals[3] >= 5) {
              if  (    qg.vals[1] * qg.vals[4] == 12) swap = 1;
              else if (qg.vals[2] * qg.vals[5] == 12) swap = 2;
          } else if (qg.vals[3] ==4 ) {
              if (qg.vals[1] == 3 && qg.vals[2] == 3) {
                  if (qm->vType[qg.verts[4] - 1] > 0 ) {
                      if (qm->vType[qg.verts[5]-1] > 0 ) continue;
                      else swap = 2;
                  } else if (qm->vType[qg.verts[5] - 1] > 0 ) {
                      if (qm->vType[qg.verts[4] - 1] > 0 ) continue;
                      else swap = 1;
                  }
                  else if (qg.vals[4] <= qg.vals[5]) swap = 1;
                  else swap = 2;
              }
              if (qg.vals[4] == 3 && qg.vals[5] == 3) {
                  if (qm->vType[qg.verts[1] - 1] > 0 ) {
                      if (qm->vType[qg.verts[2]-1] > 0 ) continue;
                      else swap = 2;
                  } else if (qm->vType[qg.verts[2] - 1] > 0 ) {
                      if (qm->vType[qg.verts[1] - 1] > 0 ) continue;
                      else swap = 1;
                  }
                  else if (qg.vals[1] <= qg.vals[2]) swap = 1;
                  else swap = 2;
              }
          }
      }
      if (swap != 0) {
          if (corner == 0 && qm->regBd == 1 &&
              (qm->vType[qg.verts[swap        ]-1] != -1 ||
                  qm->vType[qg.verts[(swap + 3)%6]-1] != -1 )) swap = 0;
          else break;
      }
  }
  if (swap  == 0) return EGADS_SUCCESS;
  return EG_swappingOperation(qm, qg, swap, &(*activity));
}


__HOST_AND_DEVICE__ static int
EG_doubleSwap(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int piv5 = -1, piv3 = -1, q5, i, adjPiv5, stat, adj[2], swap = 0, vopp3;

  *activity = 0;
  if      (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  else if (qg.vals[0] == 4 && forcing == 0) return EGADS_SUCCESS;
  else if (qg.vals[0] == 5) {
    if  (qg.vals[2] * qg. vals[4] == 15) {
      piv3 = 4;
      if (qg.vals[2] == 3) piv3 = 2;
      swap = 0;
    }
    else if (qg.vals[1] == 3 || qg.vals[5] == 3) {
      piv3 = 1; if (qg.vals[1] != 3) piv3 = 5;
      swap = 1;
    }
  }
  else if (forcing == 1) {
    if  (qg.vals[1] != 3 && qg.vals[5] != 3) return EGADS_SUCCESS;
    piv3 = 1; if (qg.vals[1] != 3) piv3 = 5;
    if (qg.vals[1] != 5 && qg.vals[5] != 5) swap = 1;
  }
  if (piv3 == -1) return EGADS_SUCCESS;
  vopp3 = (piv3 + 3)%6;
  if (qm->valence[qg.verts[vopp3] -1][1] != 4) return EGADS_SUCCESS;
  if (swap == 0) {
    piv5 = (vopp3 + 1)%6;
    if (piv5%3 == 0) piv5 = (vopp3 + 5)%6;
  }
  else {
    piv5 = (vopp3 + 1)%6;
    if (piv5%3 != 0) piv5 = (vopp3 + 5)%6;
  }
  q5   = 0;
  if (vopp3 > 3) q5 = 1;
  stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[vopp3], qg.verts[piv5], adj);
#ifdef DEBUG
  printf(" rDOUBLE SWAP :: ADJ of  %d PAIR %d %d IS %d \n",
         qg.q[q5], qg.verts[vopp3], qg.verts[piv5], adj[1]);
#endif

  if (stat != EGADS_SUCCESS || adj[1] == -1 ||
      (qm->regBd == 1 && EG_regBoundQuad(qm, adj[1]) == 1)) {
#ifdef DEBUG
    if (stat != EGADS_SUCCESS)
      printf(" EG_doubleSwap adjToPair %d -->%d !!\n", adj[1], stat);
#endif
    return stat;
  }
  i       = EG_quadVertIdx(qm, adj[1], qg.verts[vopp3]);
  adjPiv5 = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (adjPiv5 == qg.verts[piv5])
    adjPiv5 = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if (      swap == 0 && qm->valence[adjPiv5 -1][1] > 4) return EGADS_SUCCESS;
  else if ((swap == 1 && qm->valence[adjPiv5 -1][1] < 5) ||
           (qm->vType[adjPiv5 - 1] != -1 &&
            qm->vType[qg.verts[vopp3] - 1] != -1)) return EGADS_SUCCESS;
  piv5      = qg.verts[0];
  stat      = EG_swappingOperation(qm, qg, piv3, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
    printf(" EG_doubleSwap: at first swap activity %d, stat -->%d !!\n ",
           *activity, stat);
#endif
    return stat;
  }
  stat = EG_swap(qm, adj[1], &i);
  if (i == 0 && stat == EGADS_SUCCESS) {
    (void) EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
    for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == piv5) break;
    stat = EG_swappingOperation(qm, qg, swap, &i);
    if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
    else {
#ifdef DEBUG
      printf("In EG_doubleSwap undoing first swap activity %d stat %d \n!!",
             i, stat);
#endif
      return stat;
    }
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_collapse(meshMap *qm, int qID, int *activity, int forcing, int cID)
{
  int i, vC, vO, v5[5], v3[5], qb[4], stat = EGADS_SUCCESS, i3 = 0,
      val[2], links[3], act = 0;

  *activity = 0;
  if (qID <=0 || qID > qm -> totQ || qm->qIdx[4 * (qID - 1)] == -2 ||
      (qm->regBd == 1 && EG_regBoundQuad (qm, qID) == 1)) {
#ifdef DEBUG
      printf(" EG_collapse for quad %d is empty or touches bounds\n ", qID);
#endif
      return EGADS_SUCCESS;
  }
  for (i = 0; i < 4; i++) qb[i] = qm->qIdx[4 * (qID - 1) + i];
  for (v3[0] = v5[0] = i = 0; i < 4; i++) {
      val[0]   = qm->valence[qb[i] -1][1];
      if      (val[0] == 3) v3[++v3[0]] = i;
      else if (val[0] >= 5) v5[++v5[0]] = i;
      else if (val[0] == 2) {
          stat = EG_mergeVertices(qm, qID, qb[i], &(*activity));
          if (stat != EGADS_SUCCESS || *activity == 0) {
              printf("EG_collapse :: I can't fail whilst removing doublets !!!!!!!! \n ");
#ifdef DEBUG
              gnuData(qm, NULL);
              printQuad (qm, qID );
#endif
              return EGADS_GEOMERR;
          }
          break;
      }
  }
  if (forcing == 0 && *activity == 0) {
      if (v3[0] + v5[0] <= 2 || v3[0] * v5[0] == 0) return EGADS_SUCCESS;
      if (v5[0] == 1 || (v5[0] == 2 && abs(v5[1] - v5[2])%2 == 0)) {
          vC   = qb[ (v5[1] + 1)%4];
          if (EG_validCollapse(qm, qID, vC) == 1) {
              stat      = EG_mergeVertices(qm, qID, vC, &(*activity));
              if (stat != EGADS_SUCCESS) {
                  printf(" EG_collapse :: EG_mergeVertices %d activity %d \n",
                         stat, *activity);
                  return stat;
              }
          }
      }
  }
  if (*activity == 0 && forcing > 0) {
      if (cID > 0) {
          if(EG_validCollapse(qm, qID, cID) == 1)
            stat = EG_mergeVertices(qm, qID, cID, &(*activity));
      } else {
          if (EG_nValenceCount(qm, qID, 3) == 0 && forcing == 1) return EGADS_SUCCESS;
          i = 0; i3 = 0;
          while (i < 2) {
              vC = qm->qIdx[4 * (qID - 1) +  i3       ];
              vO = qm->qIdx[4 * (qID - 1) + (i3 + 2)%4];
              if ( qm->vType  [vC - 1]    == -1 && ( forcing == 3 ||
                  (qm->valence[vO - 1][2] <  10 && qm->valence[vC - 1][2] ==  3)) )
                {
                  links[0] = qm->qIdx[4 * (qID - 1) + (i3 + 1)%4];
                  links[1] = qm->qIdx[4 * (qID - 1) + (i3 + 3)%4];
                  val  [0] = qm->valence[links[0] - 1][1];
                  val  [1] = qm->valence[links[1] - 1][1];
                  if ((qm->valence[links[0] - 1][1] >= 5 ||
                       qm->valence[links[1] - 1][1] >= 5 || i == 1) &&
                      EG_validCollapse(qm, qID, vC) == 1) {
                      stat      = EG_mergeVertices(qm, qID, vC, &(*activity));
                      if (stat != EGADS_SUCCESS) {
                          printf("EG_collapse :: forcing EG_mergeVertices quad %d !! \n",
                                 stat);
#ifdef DEBUG
                          printQuad(qm, qID);
#endif
                          return stat;
                      }
                      if (*activity > 0) break;
                  }
                }
              i3++;
              if (i3 == 4) {
                  i++;
                  i3 = 0;
              }
          }
      }
  }
  if (stat  == EGADS_SUCCESS && *activity == 0) return stat;
  for (i = 0; i < 4; i++) {
      links[0] = qb[i] - 1;
      if (qm->vType[links[0]] == -2) continue;
      links[1] = qm->valence[links[0]][3] - 1;
      links[2] = qm->valence[links[0]][4] - 1;
      for (i3  = 0; i3 < 3; i3++) {
          if (qm->vType[links[i3]] != -1 ||
              qm->valence[links[i3]][2] != 2) continue;
          stat = EG_collapse(qm, qm->valence[links[i3]][0], &act, 0, 0);
          if (stat != EGADS_SUCCESS || act == 0) {
              printf(" EG_collapse removing doublet in doublet loop : V %d stat %d act %d !!!!\n ",
                     links[i3] + 1, stat, act );
              return EGADS_GEOMERR;
          }
      }
  }
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_swapCollapse(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int  stat, i, i3 = -1, q5, qC, v03, vL5adj, vOpp3, swap = 0, adj[2], v0;

  *activity = 0;
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
#ifdef DEBUG
  printf(" SwAP COLLAPSE \n");
  printQuadGroup(qm, qg);
#endif
  if (qg.vals[1] * qg.vals[2] == 15) {
      i3 = 1;
      if    (qg.vals[i3] == 5) i3 = 2;
  //    forcing = 0;
  } else if (qg.vals[4] * qg.vals[5] == 15) {
      i3 = 4;
      if (qg.vals[i3] == 5) i3 = 5;
    //  forcing = 0;
  } else if (qg.vals[0] * qg.vals[3] == 20) {
      if      (qg.vals[1] == 3) i3 = 1;
      else if (qg.vals[5] == 3) i3 = 5;
      else return EGADS_SUCCESS;
      if (i3 == 1) {
          q5   = 0;
          stat = EG_adjQtoPair(qm, qg.q[0], qg.verts[2], qg.verts[3], adj);
      }
      else {
          q5   = 1;
          stat = EG_adjQtoPair(qm, qg.q[1], qg.verts[4], qg.verts[3], adj);
      }
      if (stat != EGADS_SUCCESS || adj[1] == -1) {
          if (stat != EGADS_SUCCESS)
            printf("EG_swapCollapse centre 4-5 adj to pair %d !!\n", stat);
          return stat;
      }
      stat = EG_createQuadGroup(qm, &qg, qg.q[q5], adj[1]);
      if (stat != EGADS_SUCCESS) {
          printf("Inside EG_swapCollapse :: EG_createGroup %d !!\n ", stat);
          return stat;
      }
      return EG_swapCollapse(qm, qg, forcing, &(*activity));
  }
  else if (forcing == 1) {
      if (qg.vals[1] * qg.vals[5] == 25 || qg.vals[2] * qg.vals[4] == 25) {
          i3 = 2; if (qg.vals[1] != 5) i3 = 1;
      }
      else if (qg.vals[1] * qg.vals[4] == 15) {
          i3 = 1; if (qg.vals[1] != 3) i3 = 4;
      }
      else if (qg.vals[2] * qg.vals[5] == 15) {
          i3 = 2; if (qg.vals[2] != 3) i3 = 5;
      }
      else if (qg.vals[1] * qg.vals[5] == 9 || qg.vals[2] * qg.vals[4] == 9) {
          v0        = qg.verts[0];
          if (qm->regBd == 0 || (
              qm->vType[qg.verts[1] -1] == -1 &&
              qm->vType[qg.verts[4] -1] == -1)) {
              stat      = EG_swappingOperation(qm, qg, 1, &(*activity));
              if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
                  printf("forcing swapcollapse:: after swapping %d activity %d !!\\n ",
                         stat, *activity);
#endif
                  return stat;
              }
              qC = qg.q[0];
              if (EG_nValenceCount(qm, qC, 3) < 2) qC = qg.q[1];
              stat = EG_collapse(qm, qC, &(*activity), 1, 0);
              if (stat  != EGADS_SUCCESS || *activity == 0) {
                  (void) EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
                  for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
                  stat = EG_swappingOperation(qm, qg, swap, &i);
                  if (stat  == EGADS_SUCCESS && i == 1) *activity = 0;
              }
              return stat;
          }
      } else return EGADS_SUCCESS;
  } else return EGADS_SUCCESS;
  vOpp3 = (i3 + 3)%6;
  q5    = 0; if (vOpp3 > 3) q5 = 1;
  qC    = qg.q[(q5 + 1)%2];
  if (EG_validCollapse(qm, qC, qg.verts[i3]) == 0) return EGADS_SUCCESS;
  v03   = (vOpp3 + 1)%6;
  if (v03%3 != 0) v03 = (vOpp3 + 5)%6;
  if (qg.vals[vOpp3] == 3) return EGADS_SUCCESS;
  stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[vOpp3], qg.verts[v03], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) return stat;
  i      = EG_quadVertIdx(qm, adj[1], qg.verts[vOpp3]);
  vL5adj = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
  if (  vL5adj == qg.verts[v03]) vL5adj = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
  if ((forcing == 0 && qg.vals[vOpp3] <= 4 && qm->valence[vL5adj - 1][1] >= 4) ||
      (forcing == 1 && qm->valence[vL5adj - 1][1] > 4)) return EGADS_SUCCESS;
  stat      = EG_createQuadGroup(qm, &qg, qg.q[q5], adj[1]);
#ifdef DEBUG
  printf(" SWAP COLLAPSE SEEMS TO WORK. CREATING SWAP GROUP FORCING %d \n",
         forcing);
      printQuadGroup(qm, qg);
#endif
  if (stat != EGADS_SUCCESS) {
#ifdef DEBUG
      printf("EG_swapCollapse before swap: EG_createQuadGroup -->%d !!\n",
             stat);
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  for (swap = 0; swap < 6; swap++)
    if (qg.verts[swap] == vL5adj) break;
  v0        = qg.verts[0];
  if (qm->regBd == 1 && (
      EG_regBoundQuad(qm, qC)             == 1 ||
      qm->vType[qg.verts[ swap     ] -1] != -1 ||
      qm->vType[qg.verts[(swap+3)%6] -1] != -1)) return EGADS_SUCCESS;
#ifdef DEBUG
  printf(" SWAP IS VALID -- SWAP THRU %d \n ", swap);
#endif
  stat      = EG_swappingOperation(qm, qg, swap, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("EG_swapCollapse after swapping %d activity %d !!\n",
             swap, *activity);
#endif
      return stat;
  }
  stat = EG_collapse(qm, qC, &(*activity), 1, 0);
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf("EG_swapCollapse after collapsing %d activity %d  !!\n",
             stat, *activity);
#endif
      (void) EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
      stat = EG_swappingOperation(qm, qg, swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_doubleCollapse(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int i, stat, i3, i5, q3;

  *activity  = 0;
  if (qm->regBd == 1 && (
      EG_regBoundQuad(qm, qg.q[0]) == 1 ||
      EG_regBoundQuad(qm, qg.q[1]) == 1 )) return EGADS_SUCCESS;
  if (qg.vals[0] * qg.vals[3] == 16) {
      for (i5 = 0; i5 < 2; i5++) {
          for (i = i3 = 0; i3 < 6; i3++) {
              if (qg.vals[i3] != 3) continue;
              if (i5 == 1 || qg.vals[(i3 + 1)%6] >= 5 ||
                  qg.vals[(i3 + 5)%6] >= 5) {
                  i  =  1;
                  i5 = -1;
                  break;
              }
          }
          if (i5 == -1) break;
      }
      if (i == 0) return EGADS_SUCCESS;

      i5 = (i3 + 1)%6;
      if (i5%3 == 0) i5 = (i3 + 5)%6;
      if (forcing == 0 && (qg.vals[i5] < 5 ||
          ((qg.vals[(i3 + 3)%6] >= 5 && qg.vals[(i5 + 3)%6] != 3) ||
              qg.vals[(i3 + 3)%6] * qg.vals[(i5 + 3)%6] == 16))) return EGADS_SUCCESS;
      if (forcing == 0 && (qg.vals[i5] < 5 && qg.vals[(i5 + 3)%6] < 5 &&
          qg.vals[(i3 + 3)%6] == 4)) return EGADS_SUCCESS;
      else if (forcing == 1) {
          if (     qg.vals[(i3 + 3)%6] == 3) i5 = (i3 + 3)%6;
          else if (qg.vals[(i3 + 2)%6] >= 5) i5 = (i3 + 2)%6;
          else if (qg.vals[(i3 + 4)%6] >= 5) i5 = (i3 + 4)%6;
          else return EGADS_SUCCESS;
      }
  }
  else if (qg.vals[0] * qg.vals[3] == 12) {
      i3 = 0;
      if (qg.vals[3] == 3) i3 = 3;
      i5 = (i3 + 1)%6;
      if (qg.vals[i5] < 5) i5 = (i3 + 5)%6;
      if (qg.vals[i5] < 5 || qg.vals[(i5 + 3)%6] == 4) return EGADS_SUCCESS;
  }
  else return EGADS_SUCCESS;
  q3 = 0;
  if (forcing == 1 && qg.vals[0] * qg.vals[3] == 16) {
      if (i3 >= 3) q3 = 1;
  } else {
      if (i5 >= 3) q3 = 1;
  }
  for (i = 0; i < 2; i++) {
      if (qm->qIdx[4 * (qg.q[(q3 + i)%2] -1)] == -2) continue;
      stat        = EG_collapse(qm, qg.q[(q3 + i)%2], &i3, 1, 0);
      if (stat   == EGADS_SUCCESS) {
          if (i3 == 0) break;
          (*activity)++;
      }
      if (stat != EGADS_SUCCESS) {
          printf(" EG_doubleCollapse after EG_forceColapse quad %d->%d !!\n ",
                 qg.q[(q3 + i)%2], stat);
          return stat;
      }
  }
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_swapDoubleCollapse(meshMap *qm, quadGroup qg, int *activity)
{
  int k, swap = 1, id, j, stat, v0;

  *activity = 0;
  if (qg.vals[0] * qg.vals[3] != 20 ||
      qg.vals[2] * qg.vals[4] !=  9 ||
      EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) swap = 2;
  if (qg.vals[swap] * qg.vals[(swap+3)%6] != 12) return EGADS_SUCCESS;
  j = 2;
  if (swap == 2) j = 4;
  if (qm->vType[qm->valence[qg.verts[3] -1][3] - 1] > 0 &&
      qm->vType[qm->valence[qg.verts[j] -1][3] - 1] > 0 ) {
#ifdef DEBUG
      gnuData(qm, NULL);
      printf("SWAPPING NEAR BOUND WONT ALLOW IT. \n");
      printQuadGroup(qm, qg);
      printf(" V3 %d - %d and V ADJ %d -> %d\n",
             qg.verts[3], qm->valence[qg.verts[3] -1][3],
             qg.verts[j], qm->valence[qg.verts[j] -1][3]);
#endif
      return EGADS_SUCCESS;
  }
  v0 = qg.verts[0];
#ifdef DEBUG
  printf("SWAP DOUBLE COLLAPSE SWAP THRU %d \n ", swap);
  printQuadGroup(qm, qg);
#endif
  stat      = EG_swappingOperation(qm, qg, swap, &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" In EG_swapDoubleCollapse: Swapping went %d activity %d \n",
             stat, *activity);
#endif
      return stat;
  }
  *activity = 1;
  for (k = 0; k < 2; k++) {
      id = EG_quadVertIdx(qm, qg.q[k], qg.verts[3]);
      if (id >= 0) break;
  }
  if (id < 0) {
      printf(" I can't find vertex %d in quads!!! \n", qg.verts[3]);
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return EGADS_INDEXERR;
  }
#ifdef DEBUG
  printf(" NOW COLLAPSE %d \n", qg.q[k]);
#endif
  stat      = EG_collapse(qm, qg.q[k], &(*activity), 1, 0);
  if (stat == EGADS_SUCCESS && *activity > 0) return stat;
  if (*activity == 0 ) {
      j = qg.q[k] - 1;
#ifdef DEBUG
  printf(" UNDO SWAP \n");
#endif
      if ( qm->qIdx[4 * j] == -2 ) {
          printf(" QUAD IS VOID BUT ACTIVITY 0 \n ");
          return EGADS_INDEXERR;
      }
  }
  stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
  if (stat != EGADS_SUCCESS) {
      printf("EG_swapDoubleCollapse after invalid force collapse :: EG_createQuadGroup is %d \n ",
             stat);
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
  stat = EG_swappingOperation(qm, qg, swap, &j);
  if (stat == EGADS_SUCCESS && j > 0 ) *activity = 0;
  else {
#ifdef DEBUG
      printf(" In EG swapDoubleCollapse: I failed to collapse after swapping! s = %d act = %d \n ",
             stat, j);
#endif
  }
  return stat;
}


__HOST_AND_DEVICE__ static int EG_split(meshMap *qm, int qID, int *activity)
{
  int si, v, q, id0 = 0, i, stat, dist = 0, validSplit = 0;
  int poly[3], val[3], id6[2], links[4];

  *activity   = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_split for quad %d is -2\n ", qID );
#endif
      return EGADS_SUCCESS;
  }
  for (v = 0; v < 4; v++) {
      poly[0]   = qm->qIdx[4 * (qID - 1) + v];
      val [0]   = qm->valence[poly[0] - 1][1];
      if (val[0] < 5 || qm->vType[poly[0] - 1] >= 4) continue;
      si = poly[0] - 1;
      if (qm->star[si] == NULL ) {
          printf("EG_split star at vertex %d is NULL !!\n ", si + 1);
          return EGADS_MALLOC;
      }
      id6[0]  = -1;
      id6[1]  = -1;
      for (validSplit = q = 0; q < qm->star[si]->nQ; ++q) {
          if (qm->star[si]->quads[q] == - 1) continue;
          id0     = 2 * q + 1;
          poly[1] = qm->star[si]->verts[id0];
          val [1] = qm->valence[poly[1] -1][1];
          if (qm->regBd == 1 && qm->vType[poly[1] -1] != -1 ) continue;
          for (i = 0; i <= qm->valence[poly[0] - 1][2] - 4; ++i) {
              dist    = 4 + 2*i;
              poly[2] = qm->star[si]->verts[qm->star[si]->idxV[id0 + dist]];
              if (poly[2] < 0 || (
                  qm->regBd == 1 && qm->vType[poly[2] -1] != -1 )) continue; // ghost vertex
              val[2]  = qm->valence[poly[2] -1][1];
              if (val[1] == 3 && val[2] == 3 &&
                  qm->vType[poly[0] - 1] == -1) {
                  if ((qm->star[si]->nQ == 6 && dist == 6) || qm->star[si]->nQ != 6) {
                      validSplit = 1;
                      break;
                  }
              }
              else if (val[1] <= 3 && val[2] <= 4 && id6[0] == -1) id6[0] = id0;
              else if (val[2] <= 3 && val[1] <= 4 && id6[0] == -1)
                id6[0] = qm->star[si]->idxV[id0 + dist];
              else if (val[1] <= 4 && val[2] <= 4 && id6[0] == -1) id6[1] = id0;
#ifdef DEBUG
              printf(" USING VALS V %d = %d V %d = %d VALIDSPLIT %d id6 %d  %d\n",
                     poly[1], val[1], poly[2], val[2], validSplit, id6[0], id6[1]);
#endif
          }
          if (validSplit == 1) break;
      }
      if (validSplit == 0 && val[0] >= 6 && (id6[0] != -1 || id6[1] != -1)) {
          validSplit  = 1;
          /* boundary vertex: Special split since regular = 3 ( not val 4)  */
          if (qm->vType[poly[0] - 1] == 3) {
              for (q = 0; q < qm->star[si]->nQ; q++)
                if (qm->star[si]->quads[q] == -1) break;
              id6[0] = - 1; id6[1] = -1;
              for (i = 0; i < 2; i++) {
                  id0 = qm->star[si]->idxQ[q + i];
                  if (i == 1) dist = (qm->star[si]->nV - 1) - 4;
                  else        dist = 4 + 2 * i;
                  links[2 * i    ] = 2 * id0 + 1;
                  links[2 * i + 1] = qm->star[si]->idxV[2 * id0 + 1 + dist];
                  if (  qm->valence[qm->star[si]->verts[links[2*i     ]] - 1][1] == 4) {
                      if (qm->valence[qm->star[si]->verts[links[2*i + 1 ]] - 1][1] == 3 &&
                          id6[0] == -1) id6[0] = i;
                      else if (qm->valence[qm->star[si]->verts[links[2*i + 1]] -1][1] == 4  &&
                          id6[1] == -1) id6[1] = i;
                  }
              }
              dist    = 4;
              if      (id6[0] != -1) id0 = links[3 * id6[0]];
              else if (id6[1] != -1) id0 = links[3 * id6[1]];
              else  validSplit = 0;
          } else {
              dist   = 6;
              id0    = id6[0]; if (id0 < 0) id0 = id6[1];
              val[0] = 500;
              if (qm->regBd == 1) {
#ifdef DEBUG
                printf(" CHECK THAT WE DON't INTRODUCE A BOUNDARY BACK %d\n",
                       si + 1);
#endif
                for (id0 = i = 0 ; i < qm->star[si]->nQ; i++) {
                  id6[0] = qm->star[si]->idxV[2* i + 1       ];
                  id6[1] = qm->star[si]->idxV[2* i + 1 + dist];
                  if (qm->vType[qm->star[si]->verts[id6[0]] - 1] > 0 ||
                      qm->vType[qm->star[si]->verts[id6[1]] - 1] > 0 ) continue;
                  val[1] = qm->valence[qm->star[si]->verts[id6[0]] - 1][2];
                  val[2] = qm->valence[qm->star[si]->verts[id6[1]] - 1][2];
                  if (val[0] > val[1] * val[2]) {
                    id0    = id6[0];
                    val[0] = val[1] * val[2];
                  }
                }
                if (id0 == 0) validSplit = 0;
              }
          }
      }
    if (validSplit == 1)  {
      stat = EG_splittingOperation(qm, qm->star[si]->verts[0], qm->star[si]->verts[id0],
                                   qm->star[si]->verts[qm->star[si]->idxV[id0 + dist]],
                                   &(*activity));
      return stat;
    }
  }
  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_doubleSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int i, j, stat, piv[2] = {1, 5};
  *activity = 0;
  if (qg.vals[1] != 3) {
      piv[0] = 5; piv[1] = 1;
  }
  if (qg.vals[piv[0]] != 3) return EGADS_SUCCESS;
  if ((qm->vType[qg.verts[0]      - 1] >= 0 &&
      qm->vType[qg.verts[piv[1]] - 1] >= 0) ||
      qm->vType[qg.verts[0]      - 1] >= 4) return EGADS_SUCCESS;
  if ((forcing == 0 && (qg.vals[0] < 5 || qg.vals[1] * qg.vals[5]  != 15)) ||
      (forcing == 1 &&  qg.vals[0] * qg.vals[piv[1]] <= 16)) return EGADS_SUCCESS;
  if (qg.vals[0] < 5 || qg.vals[1] * qg.vals[5]  != 15) return EGADS_SUCCESS;

#ifdef DEBUG
  printf(" DOUBLE SPLIT FUNCTION FORCING %d\n", forcing);
  printQuadGroup(qm, qg);
#endif

  if  (qm->regBd == 1 &&
       (qm->vType[qg.verts[piv[0]] - 1] > 0 ||
        qm->vType[qg.verts[piv[1]] - 1] > 0 ) ) return EGADS_SUCCESS;
  stat      = EG_splittingOperation(qm, qg.verts[0], qg.verts[piv[0]],
                                    qg.verts[piv[1]], &(*activity));
  if (stat != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
    printf(" EG_doubleSplit: force 1st split through %d - %d -->%d activity %d !!\n ",
           qg.verts[0], qg.verts[piv[0]], stat, *activity);
#endif
    return stat;
  }
  for (j    = 0; j < 2; j++)
    if (EG_quadVertIdx(qm, qg.q[j], qg.verts[piv[1]]) >= 0) break;
  stat      = EG_split(qm, qg.q[j], &i);
  if (stat == EGADS_SUCCESS && i == 0) {
    j     = qm->valence[qg.verts[0] - 1][0];
    if ( forcing == 0) {
      stat  = EG_collapse(qm, j, &i, 1, 0);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
    }
    return stat;
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_swapDoubleSplit(meshMap *qm, quadGroup qg, int *activity)
{
  int       i5, q, i55, i3, val3, i0, v30[2], i,  stat, adj[2], q0[2], v0;
  quadGroup sqg;

  *activity = 0;
  for (q = i3 = 0; i3 < 6; i3++)
    if (qg.vals[i3] == 3 && i3%3 != 0) {
        q = 1;
        break;
    }
  if (q == 0) return EGADS_SUCCESS;
  i5  = (i3 + 3)%6;
  if (qg.vals[i5] < 5) return EGADS_SUCCESS;
  i55 = (i5 + 1)%6;
  if (qg.vals[i55] < 5) {
      i55 = (i5 + 5)%6;
      if (qg.vals[i55] < 5) return EGADS_SUCCESS;
  }
  q = 0;
  if (EG_quadVertIdx(qm, qg.q[q], qg.verts[i5]) < 0) q = 1;
  stat = EG_adjQtoPair(qm, qg.q[q], qg.verts[i5], qg.verts[i55], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) {
      if (stat != EGADS_SUCCESS)
        printf("In EG_swapDoubleSplit: adjacent to pair %d -->!!\n ", stat);
      return stat;
  }
  if ((i3 + 1)%6 == 0 || (i3 + 5)%6 == 0) i0 = 0;
  else                                    i0 = 3;
  v30[0] = qg.verts[i0];
  v30[1] = qg.verts[(i0 + 3)%6];
  q0[0]  = qg.q[0];
  q0[1]  = qg.q[1];
  val3   = qg.verts[i3];
  stat   = EG_createQuadGroup(qm, &qg, qg.q[q], adj[1]);
  if (stat != EGADS_SUCCESS) {
      printf("Inside EG_swapDoubleSplit: before swapping EG_createQuadGroup stat %d\n ",
             stat);
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  for (i0 = 0; i0 < 6; i0++) if (qg.verts[i0] == v30[0]) break;
  v0            = qg.verts[0];
  if(qm->regBd == 1 && (qm->vType[qg.verts[i0] -1] > 0 ||
      qm->vType[qg.verts[(i0+3)%6] -1] > 0 )) return EGADS_SUCCESS;
  stat          = EG_swappingOperation(qm, qg, i0, &(*activity));
  if (stat     != EGADS_SUCCESS || *activity == 0) {
      if (stat != EGADS_SUCCESS)
        printf(" EG_swapDoubleSplit error at swap: %d activity %d!!\n ",
               stat, *activity);
      return stat;
  }
  i = 0;
  if (EG_quadVertIdx(qm, q0[0], val3) < 0) i = 1;
  stat = EG_adjQtoPair(qm, q0[i], v30[0], v30[1], adj);
  if (stat != EGADS_SUCCESS || adj[1] == -1) {
      if (stat != EGADS_SUCCESS)
        printf(" EG_swapDoubleSplit after swapping adjacent to pair %d !!\n ",
               stat);
      return stat;
  }
  stat = EG_createQuadGroup(qm, &sqg, q0[i], adj[1]);
  if (stat != EGADS_SUCCESS) {
      printf("Inside EG_swapDoubleSplit: before splitting EG_createQuadGroup stat %d\n ",
             stat);
#ifdef DEBUG
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  stat   = EG_doubleSplit(qm, sqg, 0, &i);
  if (i == 0 && stat == EGADS_SUCCESS) {
      stat      = EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      if (stat != EGADS_SUCCESS) {
          printf(" In EG_swapDoubleSplit: recover swap create group %d\n ",
                 stat);
#ifdef DEBUG
          printQuadGroup(qm, qg);
#endif
          return EGADS_SUCCESS;
      }
      for (i = 0; i < 6; i++) if (qg.verts[i] == v0) break;
      stat = EG_swappingOperation(qm, qg, i, &q);
      if (stat == EGADS_SUCCESS && q == 1) *activity = 0;
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_swapSplit(meshMap *qm, quadGroup qg, int forcing, int *activity)
{
  int   stat = EGADS_SUCCESS, i, j, i3 = -1, i5 = -1, v3opp = -1,
      q5, vL5, vL5adj, swap = 0, adj[2], v0, si;

#ifdef DEBUG
  printf(" SwAP Split \n");
  printQuadGroup(qm, qg);
#endif
  *activity = 0;
  if  (qg.vals[0] * qg.vals[3] != 20 ||
      EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0) return EGADS_SUCCESS;
  if      (qg.vals[1] == 3) i3 = 1;
  else if (qg.vals[5] == 3) i3 = 5;
  for (i = 1; i < 6; i++)
    if (qg.vals[i] == 5) i5 = i;
  if (i3 != -1 && (i5 == -1 || (i5 != -1 && i5 == (i3 + 3)%6))) {
      v3opp   = (i3 +3)%6;
      i5      = v3opp;
      forcing = 0;
  }
  else if (i3 == -1 && forcing == 1 && (i5 == -1 || i5 == 2 || i5 == 4)) {
      if ( i5 == -1) v3opp = -1;
      else           v3opp = i5;
  }
  else return EGADS_SUCCESS;

  if (v3opp == -1) {
      for (i = 0; i < 2; i++) {
          j  = 2 + 2 * i;
          if (i == 0) vL5 = 1;
          else        vL5 = 5;
          if (qm->valence[qg.verts[j] -1][1]  == 3) continue;
          stat = EG_adjQtoPair(qm, qg.q[i], qg.verts[j],
                               qg.verts[vL5], adj);
          if (stat != EGADS_SUCCESS || adj[1] == -1) continue;
          q5 = EG_quadVertIdx(qm, adj[1], qg.verts[j]);
          vL5adj = qm->qIdx[4 * (adj[1] - 1) + (q5 + 1)%4];
          if (vL5adj == qg.verts[vL5])
            vL5adj = qm->qIdx[4 * (adj[1] - 1) + (q5 + 3)%4];
          if (qm->valence[vL5adj - 1][1] == 3) {
              i5   = j;
              swap = j;
              break;
          }
      }
  } else {
      vL5 = (v3opp + 1)%6;
      if (vL5%3 == 0) vL5 = (v3opp + 5)%6;
      q5 = 0;
      if (EG_quadVertIdx(qm, qg.q[q5], qg.verts[vL5]) < 0) q5 = 1;
      stat = EG_adjQtoPair(qm, qg.q[q5], qg.verts[v3opp],
                           qg.verts[vL5], adj);
      if (stat     != EGADS_SUCCESS || adj[1] == -1) {
          if (stat != EGADS_SUCCESS)
            printf("EG_swapSplit: EG_adjQtoPair from quad %d is %d\n!!",
                   qg.q[q5], stat);
          return stat;
      }
      i           = EG_quadVertIdx(qm, adj[1], qg.verts[v3opp]);
      vL5adj      = qm->qIdx[4 * (adj[1] - 1) + (i + 1)%4];
      if (vL5adj == qg.verts[vL5])
        vL5adj  = qm->qIdx[4 * (adj[1] - 1) + (i + 3)%4];
      if (i3 != -1 && (qg.vals[v3opp] == 5 || qm->valence[vL5adj - 1][1] == 3))
        swap = i3;
      else if (forcing == 1 && (qg.vals[v3opp] == 5 || qm->valence[vL5adj - 1][1] == 3))
        swap = v3opp;
  }
  if (swap%3 == 0 || (qm->vType[qg.verts [swap       ] - 1] > 0 &&
      qm->vType[qg.verts[(swap + 3)%6] - 1] > 0) ||
      (qm->regBd == 1 &&
          (qm->vType[qg.verts [swap       ] - 1] > 0 ||
              qm->vType[qg.verts[(swap + 3)%6] - 1] > 0)) ) return EGADS_SUCCESS;
  stat         = EG_swappingOperation(qm, qg, swap, &(*activity));
  v0           = qg. verts[0];
  if (stat    != EGADS_SUCCESS || *activity == 0) {
#ifdef DEBUG
      printf(" In swapSplit thru %d : EG_swappingOperation went %d activity %d !!\n ",
             swap, stat, *activity);
      printQuadGroup(qm, qg);
#endif
      return stat;
  }
  si = qg.verts[i5] - 1;
  if ( qm->star[si] == NULL ) {
      printf("In swapSplit star at %d is NULL !!\n", qg.verts[i5]);
      return EGADS_MALLOC;
  }
  for (i = 0; i < qm->star[si]->nQ; i++)
    if (qm->star[si]->verts[2 * i + 1] == qg.verts[3]) break;
  if (qm->star[si]->verts[2 * i + 1] != qg.verts[3]) {
      return EGADS_INDEXERR;
  }
  adj[0] = qm->star[si]->verts[qm->star[si]->idxV[2 * i + 1 + 6]];
  if (qm->star[si]->nQ == 5) {
      adj[1] = qm->star[si]->verts[qm->star[si]->idxV[2 * i + 1 + 4]];
      if (qm->valence[adj[1] - 1][1] < qm->valence[adj[0] - 1][1])
        adj[0] = adj[1];
  }
#ifdef DEBUG
  printf("In swapSplit call splitting op using %d %d %d !!\n",
         qg.verts[i5], qg.verts[3], adj[0]);
#endif
  if (qm->regBd ==  0 ||
      (qm->vType[qg.verts[3] -1] == -1 &&
          qm->vType[adj[0] -     1] == -1 ))
    stat     = EG_splittingOperation(qm, qg.verts[i5], qg.verts[3], adj[0], &i);
  if (stat    != EGADS_SUCCESS || i == 0) {
#ifdef DEBUG
      printf("In swapSplit splittingOperation -->%d !!\n", stat);
#endif
      (void) EG_createQuadGroup(qm, &qg, qg.q[0], qg.q[1]);
      for (swap = 0; swap < 6; swap++) if (qg.verts[swap] == v0) break;
      stat = EG_swappingOperation(qm, qg,swap, &i);
      if (stat == EGADS_SUCCESS && i == 1) *activity = 0;
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_basicOperation(meshMap *qm, int qID, int type, int *activity)
{
  int stat = EGADS_SUCCESS;
  switch (type) {
    case SWAP:
      stat = EG_swap(qm, qID, &(*activity));
      break;
    case COLLAPSE:
      stat = EG_collapse(qm, qID, &(*activity), 0, 0);
      break;
    case SPLIT:
      stat = EG_split(qm, qID, &(*activity));
      break;
  }
  if (stat != EGADS_SUCCESS && *activity == 0) return EGADS_SUCCESS;
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_composeOperation(meshMap *qm, quadGroup qg, int type, int forcing,
                    int *activity)
{
  int stat = EGADS_SUCCESS;
  switch (type) {
    case DOUBLESWAP:
      stat = EG_doubleSwap(qm, qg, forcing, &(*activity));
      break;
    case SWAPCOLLAPSE:
      stat = EG_swapCollapse(qm, qg,forcing, &(*activity));
      break;
    case DOUBLECOLLAPSE:
      stat = EG_doubleCollapse(qm, qg, forcing,  &(*activity));
      break;
    case SWAPDOUBLECOLLAPSE:
      stat = EG_swapDoubleCollapse(qm, qg, &(*activity));
      break;
    case SWAPSPLIT:
      stat = EG_swapSplit(qm, qg, forcing, &(*activity));
      break;
    case DOUBLESPLIT:
      stat = EG_doubleSplit(qm, qg, forcing, &(*activity));
      break;
    case SWAPDOUBLESPLIT:
      stat = EG_swapDoubleSplit(qm, qg, &(*activity));
      break;
  }
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_cleanQuad(meshMap *qm, int qID, int useAdj, int transfer, int forcing,
             int *activity)
{
  int stat, i, q, qadj, act = 0, aux;
  int opBasic[3] = {COLLAPSE, SWAP, SPLIT};
  int opComp [7] = {SWAPCOLLAPSE, DOUBLECOLLAPSE, SWAPDOUBLECOLLAPSE, DOUBLESWAP,
                    SWAPSPLIT,    DOUBLESPLIT,    SWAPDOUBLESPLIT};
  quadGroup qg;

  *activity = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
    printf(" EG_cleanQuad %d is empty\n", qID);
#endif
    return EGADS_SUCCESS;
  }

  i = qm->valence[qm->qIdx[4* (qID - 1)    ] - 1][1] *
      qm->valence[qm->qIdx[4* (qID - 1) + 1] - 1][1] *
      qm->valence[qm->qIdx[4* (qID - 1) + 2] - 1][1] *
      qm->valence[qm->qIdx[4* (qID - 1) + 3] - 1][1];
  if (i == 4 * 4 * 4 * 4) return EGADS_SUCCESS;
  if (transfer <= 0) {
      if (qm->valence[qm->qIdx[4* (qID - 1)    ] - 1][1] > 5 ||
          qm->valence[qm->qIdx[4* (qID - 1) + 1] - 1][1] > 5 ||
          qm->valence[qm->qIdx[4* (qID - 1) + 2] - 1][1] > 5 ||
          qm->valence[qm->qIdx[4* (qID - 1) + 3] - 1][1] > 5 ) {
          opBasic[1] = SPLIT;
          opBasic[2] = SWAP;
      }
      for (i = 0; i < 3; i++) {
          stat       = EG_basicOperation(qm, qID, opBasic[i], &act);
          *activity += act;
          if (stat != EGADS_SUCCESS) {
              printf(" EG_cleanQuad Basic operation %d !!\n ", stat);
              return stat;
          }
          if (act > 0 || qID < 0 || qm->qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
      }
  }
  if (useAdj == 0)  return EGADS_SUCCESS;
  if (qID < 0 || qm->qIdx[4 * (qID - 1)] == -2) {
#ifdef DEBUG
      printf(" EG_cleanQuad %d is empty !!\n", qID);
#endif
      return EGADS_SUCCESS;
  }

  aux = 0;
  for (i = 0; i < 7; i++) {
#ifdef DEBUG
      printf(" Q %d FORCING %d OPERATION NO %d --> %d \n", qID, aux, i, opComp[i]);
#endif
      if (transfer == 1 && opComp[i] == DOUBLESWAP) continue;
      for (q   = 0; q < 4; q++) {
          if (qID < 0 || qm->qIdx[4 * (qID - 1)] == -2 )return EGADS_SUCCESS;
          qadj = qm->qAdj[4 * (qID - 1) + q];
          if (qadj == -1 ) continue;
          stat = EG_createQuadGroup(qm, &qg, qID, qadj);
          if (stat != EGADS_SUCCESS) {
              printf("Inside EG_cleanQuad: EG_createQuadGroup stat -->%d !!\n ",
                     stat);
#ifdef DEBUG
              printQuadGroup(qm, qg);
#endif
              return stat;
          }
          stat       = EG_composeOperation(qm, qg, opComp[i], aux, &act);
          *activity += act;
          if (stat != EGADS_SUCCESS) {
              printf("In EG_cleanQuads restor mesh after EG_composeOperation stat -->%d !!\n ",
                     stat);
              return stat;
          }
          if (act > 0 || qID < 0 || qm->qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
      }
      if (i == 6 && forcing == 1 && aux == 0) {
          aux = 1;
          i  = -1;
      }
  }
  if (forcing == 1 && EG_nValenceCount(qm, qID, 3) > 1) {
#ifdef DEBUG
      printf("EG_cleanQuad FORCE COLLAPSE AROUND QUAD %d \n ", qID);
#endif

      stat = EG_collapse(qm, qID, &(*activity), 1, 0);
      if (stat != EGADS_SUCCESS || *activity > 0) {
          if (stat != EGADS_SUCCESS)
            printf("EG_cleanQuad force Collapse went %d!1\n ", stat);
          return stat;
      }
  }

  return EGADS_SUCCESS;
}


__HOST_AND_DEVICE__ static int
EG_transferValences(meshMap *qm, int *qID, int try5533, int *transfering,
                    int *activity)
{
  int       i, j, q, adj, swap = 0, stat, min, qAux[2], bdq, qv35;
  quadGroup qg;

  *activity = 0;
  if (qID[0] < 0 || qm->qIdx[4 * (qID[0] - 1)] == -2 ) {
#ifdef DEBUG
      printf(" EG_transferValences for quad %d is out of bounds or quad is regular!!\n ",
             qID[0]);
#endif
      return EGADS_SUCCESS;
  }
  bdq = EG_regBoundQuad(qm, qID[0]);
  if (*transfering == 0) {
      stat          = EG_cleanQuad(qm, qID[0], 1, 0, try5533, &(*activity));
      if (stat     != EGADS_SUCCESS || *activity > 0) {
          qID[0]    = -1;
          if (stat != EGADS_SUCCESS)
              printf("EG_transferValences: EG_cleanQuad %d -->%d!!\n",
              qID[0], stat);
          return stat;
      }
  }
  if (qm->regBd == 1 && bdq == 1) {
#ifdef DEBUG
      printf(" EG_transferValences for quad %d is out of bounds or quad is regular!!\n ",
             qID[0]);
#endif
      return EGADS_SUCCESS;
  }
  j = qm -> valence[qm->qIdx[4 * (qID[0] - 1)    ] - 1][1] *
      qm -> valence[qm->qIdx[4 * (qID[0] - 1) + 1] - 1][1] *
      qm -> valence[qm->qIdx[4 * (qID[0] - 1) + 2] - 1][1] *
      qm -> valence[qm->qIdx[4 * (qID[0] - 1) + 3] - 1][1];
  qv35 = 0;
  if (j % 15 == 0) {
      qv35 = 1;
      if (j % 9 == 0 || j % 25 == 0) qv35 = 2;
  }
  if (qID[0] == -1 || qm->qIdx[4 * (qID[0] - 1)] == -2 ||
           j ==  4 * 4 * 4 * 4) return EGADS_SUCCESS;
  for (q   = 0; q < 4; q++) {
      swap = 0;
      adj  = qm->qAdj[4 * (qID[0] - 1) + q];
      if (adj == -1 || adj == qID[1] ||
         (qm->regBd == 1 && bdq == -1 &&  EG_regBoundQuad(qm, adj) ==  1 )) {
          if (*transfering == 0 && qv35 == 1 ) {
              EG_collapse(qm, qID[0], &(*activity), 1, 0);
              if (stat     != EGADS_SUCCESS || *activity > 0) {
                  qID[0]    = -1;
                  if (stat != EGADS_SUCCESS)
                      printf("EG_transferValences: EG_cleanQuad %d -->%d!!\n",
                      qID[0], stat);
                  return stat;
              }
          }
          continue;
      }
      if (*transfering == 0) {
          stat      = EG_cleanQuad(qm, adj, 1, 0, try5533, &(*activity));
          if (stat != EGADS_SUCCESS || *activity > 0) {
              if (stat != EGADS_SUCCESS)
                  printf(" EG_TransferValence EG_cleanQuad -->%d !!\n ", stat);
              qID[0] = -1;
              return stat;
          }
      }
      if (qm->qAdj[4 * (qID[0] - 1) + q] != adj) continue;
      stat = EG_createQuadGroup(qm, &qg, qID[0], adj);
      if (stat != EGADS_SUCCESS) {
          printf(" Inside EG_transferValences EG_createQuadGroup %d !!\n", stat);
#ifdef DEBUG
          printQuadGroup(qm, qg);
#endif
          return stat;
      }
      if (*transfering == 0 && qg.vals[0] * qg.vals[3] == 15 && qID[1] != -1) {
          j = 0;
          if (qg.q[j] == qID[0]) j = 1;
          qID[0] = qg.q[j];
          qID[1] = qg.q[(j + 1)%2];
          for (i = 0 ; i < 4; i++) {
               j = qm->qIdx [4 * (qID[0] - 1 ) + i ] - 1;
              if (qm->vType[qm->valence[j][3] - 1] > 0) {
                  stat      = EG_collapse(qm, qID[0], &(*activity), 1, 0);
                  if (stat != EGADS_SUCCESS)
                    printf("EG_transferValences: forceCollapse gave %d !!\n",
                           stat);
                  qID[0]    = qID[1];
                  qID[1]    = -1;
                  return stat;
              }
          }
          return EG_transferValences(qm, qID, 0, &(*transfering), &(*activity));
      }
      if     (qg.vals[0] * qg.vals[3] >= 20) min = 12;
      else if (try5533 == 1 && *transfering == 0) {
          if      (qg.vals[0] * qg.vals[3] >= 25) min = 16;
          else if (qg.vals[0] * qg.vals[3] == 16) min =  9;
          else continue;
      } else continue;
      if (     qg.vals[1] * qg.vals[4] <= min) swap = 1;
      else if (qg.vals[2] * qg.vals[5] <= min) swap = 2;
      if (EG_validSwap(qm, qg.verts[0], qg.verts[3]) == 0 || swap == 0) continue;
      stat          = EG_swappingOperation(qm, qg, swap, &(*activity));
      if (stat     != EGADS_SUCCESS || *activity == 0) {
          if (stat != EGADS_SUCCESS)
#ifdef DEBUG
          printf(" EG_transferValences EG_swappingOperation stat %d activity %d !!\n",
                 stat, *activity);
#endif
          return stat;
      }
      i        = 0;
      if (min == 9 || min == 16 || *transfering == 2) *transfering = 1;
      else                         *transfering = 0;
      if (qID[1] > 0 && qm->qIdx[4 * (qID[1] - 1)] != -1 &&
          qg.q[0] == qID[1]) i = 1;
      else if (EG_nValenceCount(qm,qg.q[0], 5) == 2 &&
               EG_nValenceCount(qm,qg.q[1], 5) >= 1 &&
               EG_nValenceCount(qm,qg.q[1], 3) == 1) i = 1;
      else if (EG_nValenceCount(qm,qg.q[1], 4) <
               EG_nValenceCount(qm,qg.q[0], 4)) i = 1;
      qID[0] = qg.q[ i       ];
      qID[1] = qg.q[(i + 1)%2];
#ifdef DEBUG
      printf(" ASSIGNED QUAS %d %d TRANSFER %d\n", qID[0], qID[1], *transfering);
#endif
      if (*transfering == 1) {
          for (j = 0; j < 2; j++) {
              i         = 1;
              qAux[0]   = qID[ j       ];
              qAux[1]   = qID[(j + 1)%2];
              stat      = EG_transferValences(qm, qAux, 0, &i, &min);
              if (stat != EGADS_SUCCESS) {
                  printf(" EG_transferValences: separating valences after forcing %d!!\n",
                  stat);
                  return stat;
              }
          }
          *transfering = 0;
          qID[0]       = qAux[0];
          qID[1]       = qAux[1];
      }
      if (*activity > 0) break;
  }
  return stat;
}


__HOST_AND_DEVICE__ void EG_destroyMeshMap(bodyQuad *bodydata)
{
  int i, j;

  if (bodydata->qm == NULL) return;
  for (i = 0; i < bodydata->nfaces; ++i) {
      if (bodydata->qm[i]) {
          for (j = 0; j < bodydata->qm[i]->sizeV; ++j) {
              if (bodydata->qm[i]->valence)
                EG_free(bodydata->qm[i]->valence[j]);
              if (bodydata->qm[i]->star)
                EG_freeStar(&bodydata->qm[i]->star[j]);
          }
          EG_free(bodydata->qm[i]->valence);
          EG_free(bodydata->qm[i]->qIdx);
          EG_free(bodydata->qm[i]->qAdj);
          EG_free(bodydata->qm[i]->uvs);
          EG_free(bodydata->qm[i]->xyzs);
          EG_free(bodydata->qm[i]->remQ);
          EG_free(bodydata->qm[i]->remV);
          EG_free(bodydata->qm[i]->vType);
          EG_free(bodydata->qm[i]->star);
          EG_free(bodydata->qm[i]->vInv);
          EG_free(bodydata->qm[i]->bdAng);
          EG_free(bodydata->qm[i]->degen);
          EG_free(bodydata->qm[i]);
      }
  }
  EG_free(bodydata->qm);
  bodydata->qm = NULL;
}


__HOST_AND_DEVICE__ static int resizeQm(meshMap *qm)
{
  int    stat = EGADS_SUCCESS, nV, vRem, nQ, qRem, i, j, k;
  int    *vpiv = NULL, *qpiv = NULL;
  int    *verts = NULL, *adjs = NULL, *type = NULL, bV, bQ;
  double *uvxyz = NULL;

  if(qm->fID == 0) return EGADS_SUCCESS;
  vpiv = (int *) EG_alloc(qm->totV * sizeof(int));
  qpiv = (int *) EG_alloc(qm->totQ * sizeof(int));
  if (vpiv == NULL || qpiv == NULL) return EGADS_MALLOC;
  for (nV = vRem = i = 0; i < qm->totV; i++) {
      if (qm->vType[i] != -2) vpiv[i] = nV++;
      else {
          vpiv[i] = -2;
          vRem++;
      }
  }
  for (nQ = qRem = i = 0; i < qm->totQ; i++) {
      if (qm->qIdx[4*i] != -2) qpiv[i] = nQ++;
      else {
          qpiv[i] = -2;
          qRem++;
      }
  }
  if (vRem != qRem) {
      printf(" In resizeQm: I have %d removed vertices and %d quads!! they should match!!!!!\n ",
             vRem, qRem);
      stat = EGADS_INDEXERR;
      goto cleanup;
  }
  if (vRem == 0) {
      EG_free(vpiv);
      EG_free(qpiv);
      return EGADS_SUCCESS;
  }
  bV    = qm->totV;
  bQ    = qm->totQ;
  uvxyz = (double *) EG_alloc(5 * bV *sizeof(double));
  type  = (int    *) EG_alloc(    bV *sizeof(int));
  verts = (int    *) EG_alloc(4 * bQ *sizeof(int));
  adjs  = (int    *) EG_alloc(4 * bQ *sizeof(int));
  if (uvxyz == NULL || type == NULL || verts == NULL || adjs == NULL) {
      stat = EGADS_MALLOC;
      goto cleanup;
  }
  for (j = 0; j < bQ; ++j) {
      for (k = 0; k < 4; ++k) {
          verts[4*j + k] = qm->qIdx[4*j + k];
          adjs [4*j + k] = qm->qAdj[4*j + k];
      }
  }
  for (j = 0; j < bV; ++j) {
      type [  j    ] = qm->vType[  j    ];
      uvxyz[5*j    ] = qm->uvs  [2*j    ];
      uvxyz[5*j + 1] = qm->uvs  [2*j + 1];
      uvxyz[5*j + 2] = qm->xyzs [3*j    ];
      uvxyz[5*j + 3] = qm->xyzs [3*j + 1];
      uvxyz[5*j + 4] = qm->xyzs [3*j + 2];
  }
  qm->totV    = nV;
  qm->totQ    = nQ;
  qm->remQ[0] = 0;
  qm->remV[0] = 0;
  for (i = 0; i < bQ; i++) {
      if (qpiv[i] == -2) continue;
      for (k = 0; k < 4; ++k) {
          qm->qIdx[4*qpiv[i] + k ] = vpiv[verts[4*i + k] - 1] + 1;
          if (adjs[4*i + k] == -1) qm->qAdj[4*qpiv[i] + k] = - 1;
          else                     qm->qAdj[4*qpiv[i] + k] = qpiv[adjs[4*i + k] - 1] + 1;
          qm->valence[qm->qIdx[4*qpiv[i] + k] - 1][0] = qpiv[i] + 1;
      }
  }
  for (i = 0; i < bV; i++) {
      if (vpiv[i] == -2) continue;
      j = vpiv[i];
      qm->vType  [j      ] = type [  i    ];
      qm->uvs    [2*j    ] = uvxyz[5*i    ];
      qm->uvs    [2*j + 1] = uvxyz[5*i + 1];
      qm->xyzs   [3*j    ] = uvxyz[5*i + 2];
      qm->xyzs   [3*j + 1] = uvxyz[5*i + 3];
      qm->xyzs   [3*j + 2] = uvxyz[5*i + 4];
      EG_setValence (qm, j + 1);
  }
  cleanup:
  EG_free(vpiv);
  EG_free(qpiv);
  EG_free(uvxyz);
  EG_free(type);
  EG_free(verts);
  EG_free(adjs);
  return stat;
}


__HOST_AND_DEVICE__ static int
EG_cleanNeighborhood(meshMap *qm, int qID,  int transfer, int *activity)
{
  int si, i, act, stat = 0, count = 0, j, v[4];

  *activity = 0;
  if (qm->qIdx[4 * (qID - 1)] == -2) return EGADS_SUCCESS;
  for (i = 0; i < 4; i++) v[i]  = qm->qIdx[4 * (qID - 1) + i];
  stat      = EG_cleanQuad(qm, qID, 1, transfer, 0, &act);
  if (stat != EGADS_SUCCESS) return stat;
  *activity += act;
  for (i = 0; i < 4; i++) {
      if (qm->vType[v[i] - 1] == -2) continue;
      si = v[i] - 1;
      if ( qm->star[si] == NULL ) {
          printf(" EG_cleanNeighborhood star at %d is NULL !!\n", v[i]);
          return EGADS_MALLOC;
      }
      for (j = 0; j < qm->star[si]->nQ; j++) {
          if (qm->star[si]->quads[j] == -1) continue;
          if (EG_cleanQuad(qm, qm->star[si]->quads[j], 1, transfer, 0, &act) != EGADS_SUCCESS) {
              return stat;
          }
          count += act;
      }
  }
  *activity += count;
#ifdef DEBUG
  return checkMesh(qm);
#else
  return EGADS_SUCCESS;
#endif
}


#ifdef REPORT
__HOST_AND_DEVICE__ static void meshProperties(meshMap *qm, int sweep)
{
  int j, k, i, v[4], quarts[8], n, Nx = 6;
  double minT, maxT, minR, maxR, uv0[2], uv1[2], uv3[2], p[18];
  FILE *fil =  NULL;
  char buffer[100];

  if (sweep == 0)
    snprintf(buffer,100, "InitialMeshProperties_%d.txt", qm->fID);
  else
    snprintf(buffer,100, "finalMeshProperties_%d.txt", qm->fID);
  fil = fopen(buffer,"w");
  if (fil == NULL ) {
      printf(" meshProperties Error writing in file %s\n", buffer);
      return;
  }
  minT = 2.0 * PI; maxT = 0.0;
  minR = 1.0; maxR = 0.0;
  quarts[0] = quarts[1] = quarts[2] = quarts[3] = 0;
  quarts[4] = quarts[5] = quarts[6] = quarts[7] = 0;
  for (k = 0 ; k < qm->totV; k++) {
      if (qm->vType[k] == -2) continue;
      EG_vertexArea(qm, k + 1);
      for ( j = 0; j < qm->star[k]->nQ; j++) {
          if (qm->star[k]->quads[j] == -1) continue;
          minT = MIN(minT, qm->star[k]->angle[j]);
          maxT = MAX(maxT, qm->star[k]->angle[j]);
          minR = MIN(minR, qm->star[k]->ratio[j]);
          maxR = MAX(maxR, qm->star[k]->ratio[j]);
          if     (qm->star[k]->angle[j] < (20. /180. * PI)) quarts[0]++;
          else if(qm->star[k]->angle[j] < (80. /180. * PI)) quarts[1]++;
          else if(qm->star[k]->angle[j] < (160./180. * PI)) quarts[2]++;
          else                                                quarts[3]++;
          if     (qm->star[k]->ratio[j] < 0.1 ) quarts[4]++;
          else if(qm->star[k]->ratio[j] < 0.25) quarts[5]++;
          else if(qm->star[k]->ratio[j] < 0.6 ) quarts[6]++;
		  else if(qm->star[k]->ratio[j] < 0.75) quarts[6]++;
          else                                  quarts[7]++;
      }
  }
  printf("\n\n FACE %d ---> MIN MAX ANGLES (%lf %lf)   MIN MAX RATIOS (%lf  %lf )\n",
         qm->fID, minR, maxT, minR, maxR);
  printf("Vertices by angles : \n"
      " [0, 15] = %d [15,40] = %d [40, 155] = %d [155, 180] = %d  \n",
      quarts[0] ,quarts[1] ,quarts[2],quarts[3]);
  printf("Vertices by triangle aspect ratio : \n"
      " [0, 0.1] = %d [0.1, 0.25] = %d [0.26,0.6] %d [0.6,1] = %d\n",
      quarts[5] ,quarts[5] ,quarts[6], quarts[7]);
  if (fil) {
      fprintf(fil, "\n\n #FACE %d ---> MIN MAX ANGLES (%lf %lf)   MIN MAX RATIOS (%lf  %lf )\n",
              qm->fID, minR, maxT, minR, maxR);
      fprintf(fil, "#Vertices by angles : \n"
              " #[0, 15] = %d [15,40] = %d [40, 155] = %d [155, 180] = %d  \n",
              quarts[0] ,quarts[1] ,quarts[2],quarts[3]);
      fprintf(fil, "#Vertices by triangle aspect ratio : \n"
              " #[0, 0.1] = %d [0.1, 0.25] = %d [0.26,0.6] %d [0.6,1] = %d\n",
              quarts[5] ,quarts[5] ,quarts[6], quarts[7]);
  }
  for (k = 0 ; k < qm->totQ; k++) {
      minT = 2.0 * PI; maxT = 0.0;
      minR = 1.0;      maxR = 0.0;
      v[0] = v[1] = v[2] = v[3] = 0;
      for (i = 0 ; i < 4; i++) {
          v[i] = qm->qIdx[4 * k + i] -1;
          if (qm->vType[v[i]] == -2) break;
          for (j = 0 ; j < qm->star[v[i]]->nQ; j++)
            if (qm->star[v[i]]->quads[j] == k + 1) break;
          minT = MIN(minT, qm->star[v[i]]->angle[j]);
          maxT = MAX(maxT, qm->star[v[i]]->angle[j]);
          minR = MIN(minR, qm->star[v[i]]->ratio[j]);
          maxR = MAX(maxR, qm->star[v[i]]->ratio[j]);
      }
      uv0[0] = qm->uvs[2 * v[0]    ];
      uv0[1] = qm->uvs[2 * v[0] + 1];
      uv1[0] = qm->uvs[2 * v[1]    ];
      uv1[1] = qm->uvs[2 * v[1] + 1];
      for (n = 0; n <= Nx; n++) {
          uv3[0] = uv0[0] + (qm->uvs[2 * v[3]    ] - qm->uvs[2 * v[0]    ]) * n / Nx;
          uv3[1] = uv0[1] + (qm->uvs[2 * v[3] + 1] - qm->uvs[2 * v[0] + 1]) * n / Nx;
          if (EG_evaluate(qm->face, uv3, p) != EGADS_SUCCESS) continue;
          fprintf(fil,"%f %f %f %f %f %f %f %d \n",
                  p[0], p[1], p[2], minT* 180.0 / PI,
                  maxT* 180.0 / PI, minR, maxR, k + 1);
          uv3[0] = uv1[0] + (qm->uvs[2 * v[2]    ] - qm->uvs[2 * v[1]    ]) * n / Nx;
          uv3[1] = uv1[1] + (qm->uvs[2 * v[2] + 1] - qm->uvs[2 * v[1] + 1]) * n / Nx;
          if (EG_evaluate(qm->face, uv3, p) != EGADS_SUCCESS) continue;
          fprintf(fil,"%f %f %f %f %f %f %f %d\n\n\n",
                  p[0], p[1], p[2], minT* 180.0 / PI,
                  maxT* 180.0 / PI, minR, maxR, k +1);
      }
      uv1[0] = qm->uvs[2 * v[3]    ];
      uv1[1] = qm->uvs[2 * v[3] + 1];
      for (n = 0; n <= Nx; n++) {
          uv3[0] = uv0[0] + (qm->uvs[2 * v[1]    ] - qm->uvs[2 * v[0]    ]) * n / Nx;
          uv3[1] = uv0[1] + (qm->uvs[2 * v[1] + 1] - qm->uvs[2 * v[0] + 1]) * n / Nx;
          if (EG_evaluate(qm->face, uv3, p) != EGADS_SUCCESS) continue;
          fprintf(fil,"%f %f %f %f %f %f %f %d\n",
                  p[0], p[1], p[2], minT* 180.0 / PI,
                  maxT* 180.0 / PI, minR, maxR, k + 1);
          uv3[0] = uv1[0] + (qm->uvs[2 * v[2]    ] - qm->uvs[2 * v[3]    ]) * n / Nx;
          uv3[1] = uv1[1] + (qm->uvs[2 * v[2] + 1] - qm->uvs[2 * v[3] + 1]) * n / Nx;
          if (EG_evaluate(qm->face, uv3, p) != EGADS_SUCCESS) continue;
          fprintf(fil,"%f %f %f %f %f %f %f %d\n\n\n",
                  p[0], p[1], p[2], minT* 180.0 / PI,
                  maxT* 180.0 / PI, minR, maxR, k + 1);
      }
      if (fil) fprintf(fil,"\n\n");
  }
  if (fil) fclose(fil);
  minT = 2.0 * PI; maxT = 0.0;
  minR = 1.0; maxR = 0.0;
  quarts[0] = quarts[1] = quarts[2] = quarts[3] = 0;
  quarts[4] = quarts[5] = quarts[6] = quarts[7] = 0;
  for (k = 0 ; k < qm->totV; k++) {
      if (qm->vType[k] < 0) continue;
      for ( j = 0; j < qm->star[k]->nQ; j++) {
          if (qm->star[k]->quads[j] == -1) continue;
          minT = MIN(minT, qm->star[k]->angle[j]);
          maxT = MAX(maxT, qm->star[k]->angle[j]);
          minR = MIN(minR, qm->star[k]->ratio[j]);
          maxR = MAX(maxR, qm->star[k]->ratio[j]);
          if     (qm->star[k]->angle[j] < (15. /180. * PI)) quarts[0]++;
          else if(qm->star[k]->angle[j] < (40. /180. * PI)) quarts[1]++;
          else if(qm->star[k]->angle[j] < (155./180. * PI)) quarts[2]++;
          else                                                quarts[3]++;
          if     (qm->star[k]->ratio[j] < 0.1 ) quarts[4]++;
          else if(qm->star[k]->angle[j] < 0.25) quarts[5]++;
          else if(qm->star[k]->angle[j] < 0.6 ) quarts[6]++;
          else                                    quarts[7]++;
      }
  }
  minT *= 180.0 / PI;
  maxT *= 180.0 / PI;
  printf("AT THE BOUNDARY \nVertices by angles : \n"
         " [0, 15] = %d [15,40] = %d [40, 155] %d [155,180] = %d \n",
         quarts[0], quarts[1], quarts[2], quarts[3]);
  printf("Vertices by triangle aspect ratio : \n"
         " [0, 0.1] = %d [0.1, 0.25] = %d [0.26,0.6] %d [0.6,1] = %d\n",
         quarts[5], quarts[5], quarts[6], quarts[7]);
  printf("minT = %lf  maxT = %lf\n", minT, maxT);
}
#endif


__HOST_AND_DEVICE__ int EG_meshRegularization(meshMap *qm)
{
  int    i, j, k, s, q, stat = EGADS_SUCCESS,  ni = 0, n0 = 0;
  int    ITMAX, it = 0, activity = 0, totActivity = 0, loopact;
  int    iV0, iV, qPair[4], prevPair[2], totV, totV0, vQ0, vQ, transfer = 0, *skipQuad = NULL, sq;
  double minArea, maxArea, minArea0, maxArea0, avArea0, avArea, *qArea = NULL;

  if(qm      == NULL) return EGADS_SUCCESS;
  if(qm->fID == 0   ) return EGADS_SUCCESS;
#ifdef REPORT
  char   buffer[100];
  clock_t time0, time1, ttime;
  int min;
  float secs;
  printf("__________________________IO FILES ______________________________\n");
  snprintf(buffer,100, "gnuInit_%d.txt", qm->fID);
  gnuData(qm, buffer);
  snprintf(buffer,100, "wvsInit_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100, "initMeshStats_%d.txt", qm->fID);
  meshStats(qm, 0);
  meshProperties(qm,0);
  qm->invsteps = 0;
  printf("*****************************************************************\n");
#endif


  // GET RANGE FOR EACH POINT
  meshCount(qm, &iV0, &totV0, &vQ0);
  // PREPROCESS: Reduce the number of quads by forcing collapses everywhere
#ifdef REPORT
  time0 = clock();
#endif
#ifdef DEBUG
  printf(" IMPOSING REG BOUNDS %d\n", qm->regBd);
#endif
  skipQuad = (int *)    EG_alloc (qm->totQ * sizeof(int   ));
  qArea    = (double *) EG_alloc (qm->totQ * sizeof(double));
  if (skipQuad == NULL || qArea == NULL) {
      qm->fID = 0;
      return EGADS_MALLOC;
  }
  ITMAX = qm->totQ;
  qm->pp = 1;
  minArea0 = 10000000.00;
  maxArea0 = avArea0 = 0.0;
  minArea  = 10000000.00;
  maxArea  = avArea = 0.0;
  for (n0 = sq = loopact = it = 0 ; it < ITMAX; it++) {
      minArea = 10000000.00;
      maxArea = 0.0;
      avArea  = 0.0;
      for (j = i = 0; i < qm->totQ; i++) {
          qArea[i] = 0.0;
          if (qm->qIdx[ 4 * i ] == -2) continue;
          qArea[i] = EG_quadSize(qm, i + 1, NULL);
          if (it == 0) {
              minArea0 = MIN(qArea[i], minArea0);
              maxArea0 = MAX(qArea[i], maxArea0);
              avArea0 += qArea[i];
          }
          minArea = MIN(qArea[i], minArea);
          maxArea = MAX(qArea[i], maxArea);
          avArea += qArea[i];
          j++;
      }
      if (it ==0) avArea0 /= (double)j;
      avArea /= (double)j;
      if (maxArea >= 1.5 * maxArea0 || minArea >= avArea0) {
#ifdef DEBUG
         printf(" LEAVE BECAUSE %lf > %lf OR %lf > %lf \n",
                maxArea, 1.5 * maxArea0, avArea, 1.5 * avArea0);
#endif
break;
      }
#ifdef DEBUG
      if ( it % 100 == 0 )
        printf(" PREPROCESS %d  / %d :: Areas Min max Av %1.2e %1.2e %1.2e ori %1.2e %1.2e %1.2e inv steps %d act %d\n",
               it, ITMAX, minArea, maxArea, avArea, minArea0, maxArea0, avArea0
               ,qm->invsteps, loopact);
#endif

      if (n0 == 0) {
          n0 = 1; // only first round is allowed to have invalid quads (self int)
          for (ni = i = 0; i < qm->totV; i++) {
              if(qm->vType[i] != -1) continue;
              if (EG_vertexArea(qm, i + 1) > QA2) ni++;
          }
          if (ni > 0 ) {
              //#ifdef REPORT
              printf("\n\n Attention! face %d has initally %d / %d INVALID vertex positions \n",
                     qm->fID, ni, qm->totV);
              //#endif
              qm->vInv = (int *) EG_alloc ((ni + 1) * sizeof(int));
              if (qm->vInv == NULL) {
                  EG_free (skipQuad);
                  EG_free (qArea   );
                  qm->fID = 0;
                  return EGADS_MALLOC;
              }
              for (s = i = 0; i < qm->totV; i++) {
                  if (qm->vType[i] != -1) continue;
                  for (j = 0; j < qm->star[i]->nQ; j++) {
                      if (s >= ni) {
#ifdef DEBUG
                          printf(" OUT OF BOUNDS !!! %d > %d \n", s, ni);
#endif
                          break;
                      }
                      if (qm->star[i]->area[j] > QA2 &&
                              (s == 0 || inList(s, &qm->vInv[1], i + 1) == -1)) {
                          qm->vInv[++s] = i + 1;
                          break;
                      }
                  }
                  if (s >= ni) break;
              }
              qm->vInv[0] = ni;
          }
      }
      if (qm->vInv != NULL && qm->vInv[0] > 0) {
          for (j = 0; j < qm->vInv[0]; j++) {
              k = qm->vInv[1 + j] - 1;
              if (k < 0 || k >= qm->totV)  continue;
              for (i = 0; i < qm->star[k]->nQ; i++) {
                  if(qm->star[k]->quads[i] == -1) continue;
                  stat  = EG_collapse(qm, qm->star[k]->quads[i], &activity, 3, 0);
                  if (stat != EGADS_SUCCESS) {
                      printf(" EG_meshRegularization :: preprocess EG_collapse for quad %d -->%d!!\n ",
                             qm->vInv[1 + j], stat);
                      EG_free (skipQuad);
                      EG_free (qArea   );
                      EG_free(qm->vInv);
                      qm->fID = 0;
                      return stat;
                  }
                  if (activity > 0) {
                      loopact += activity;
                      if (qm->vType[k] == -2 ) qm->vInv[1 + j] = -2;
                      else {
                          for (s = i = 0; i < qm->star[k]->nQ; i++) {
                              if(qm->star[k]->quads[i] == -1) continue;
                              s += qm->star[k]->area[i];
                          }
                          if (s <= QA2) qm->vInv[1 + j] = -2;
                      }
                      if (qm->vInv[1 + j] == -2) break;
                  }
              }
          }
      }
      for (k = i = 0; i < qm->totQ; i++) {
          if (qm->qIdx[4 * i] == -2 || inList (sq, skipQuad, i + 1) != -1) continue;
          if (qArea[i] <= avArea ) {
              stat  = EG_collapse(qm, i+1, &activity, 2, 0);
              if (stat != EGADS_SUCCESS) {
                  printf(" EG_meshRegularization :: preprocess EG_collapse for quad %d -->%d!!\n ",
                         k, stat);
              }
              loopact += activity;
              if (activity == 0) {
                  if (sq <= qm->totQ -1) skipQuad[sq++] = i + 1;
              }
              else break;
          }
      }
      if (it < ITMAX / 2 && sq == qm->totQ) sq = 0 ;
  }
#ifdef REPORT
#ifndef __CUDA_ARCH__
  fprintf(stderr, "RATIO MIN MAX AREA %lf %lf  av %lf \n",
          minArea, maxArea, avArea);
#else
  printf("RATIO MIN MAX AREA %lf %lf  av %lf \n", minArea, maxArea, avArea);
#endif
#endif
  qm->minArea = minArea0;
  qm->maxArea = maxArea0;
  qm->avArea  =  avArea0;
#ifdef REPORT
#ifndef __CUDA_ARCH__
  fprintf(stderr, " MIN MAX AREA %lf %lf AVERAGE %lf RATIO MIN / AV %lf \n",
#else
  printf(" MIN MAX AREA %lf %lf AVERAGE %lf RATIO MIN / AV %lf \n",
#endif
         qm ->minArea, qm -> maxArea, qm ->avArea, qm->minArea/qm->avArea);
#endif
  if ( qm ->minArea / qm -> avArea > 0.5 ) {
#ifdef REPORT
      printf(" DEACTIVATING SMALL AREA: %lf  %lf ratio %lf\n",
             qm->minArea, qm->avArea, qm->minArea / qm->avArea);
      qm->minArea = qm->maxArea = qm->avArea = 0.0;
#endif
  }
  EG_free (skipQuad);
  EG_free (qArea   );
  qm->pp = 0;
#ifdef REPORT
  time1 = clock();
  ttime = (double)(time1 - time0) / CLOCKS_PER_SEC;
  min   = floor(ttime) / 60;
  secs  = ttime -  (double)min * 60.0;
  printf("\n\n");
  printf("****************** PREPROCESSING: %d collapses ******************\n", loopact);
  printf(" Min Max areas %1.2e %1.2e vs %1.2e %1.2e\n", minArea, maxArea, minArea0, maxArea0);
  meshCount(qm, &iV, &totV, &vQ);
  printf(" Original mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
         vQ0, iV0, totV0, (double) iV0 * 100.0 / (double)totV0);
  printf(" Current  mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
         vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Invalid steps  %d\t CPU time: %d mins and %f sec\n",
         qm->invsteps, min, secs);
  printf("__________________________IO FILES ______________________________\n");
  snprintf(buffer,100, "gnuPreprocess_%d.txt", qm->fID);
  gnuData(qm, buffer);
  snprintf(buffer,100, "wvsPreprocess_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100, "preprocessMeshStats_%d.txt", qm->fID);
  meshStats(qm, 10);
  qm->invsteps = 0;
  printf("*****************************************************************\n");
  time0 = clock();
#endif
  k = 0;
  if (qm->vInv != NULL && qm->vInv[0] > 0 ) {
      for (j = 0 ; j < qm->vInv[0]; j++ ) {
          if ( qm->vInv[ 1 + j ] == -2 ) continue;
          k++;
      }
      EG_free(qm->vInv);
      qm->vInv = NULL;
  }
  if (k == 1) {
#ifndef __CUDA_ARCH__
      fprintf(stderr,"After Preprocessing face %d there are still %d invalid quads from the original tessellation "
              "that couldn't be eliminated. Let's see what happens....!!!! \n ",
              qm->fID, k);
#else
      printf("After Preprocessing face %d there are still %d invalid quads from the original tessellation "
             "that couldn't be eliminated. Let's see what happens....!!!! \n ",
             qm->fID, k);
#endif
  }
  stat = resizeQm(qm);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_meshRegularization resizing after preprocess %d !!\n", stat);
      qm->fID = 0;
      return stat;
  }
  ITMAX = 100;
  for (it = 0; it < ITMAX; it++) { //deactivae regular bound:: bring it back during transfer
      totActivity = 0;
#ifdef DEBUG
      printf("CHECKIIN IT %d / %d \n", it, ITMAX);
      gnuData(qm, NULL);
#endif
      for (q = 0; q < qm->totQ; q++) {
          if (qm->qIdx[4 * q] == -2) continue; //can be a deleted quad
          stat = EG_cleanNeighborhood(qm, q + 1, 0, &activity);
          if (stat != EGADS_SUCCESS) {
              printf(" In EG_cleanMesh: EG_CleanNeighborhood for quad %d -->%d!!\n ",
                     q + 1, stat);
              qm->fID = 0;
              return stat;
          }
          totActivity += activity;
          if (it == -1) {
              stat = resizeQm(qm);
              if (stat != EGADS_SUCCESS)
                printf(" In EG_cleanMesh: resizeQm = %d\n", stat);
            //printf(" TRANSFER VALENCES  IV %d vQ %d totV %d \n ", iV ,vQ, totV);
            //snprintf(buffer,100, "r1_%d_%d.txt", qm->fID, q);
            //gnuData(qm, buffer);
            //snprintf(buffer,100, "wvsr1_%d_%d.txt", qm->fID,q);
            //wvsData(qm, buffer);
          }
      }
      meshCount(qm, &iV, &totV, &vQ);
      if (totActivity == 0 || iV <= 2) break;
  }
#ifdef REPORT
  (void) resizeQm(qm);
  printf(" TRANSFER VALENCES  IV %d vQ %d totV %d \n ", iV ,vQ, totV);
  snprintf(buffer,100, "notransfer_%d.txt", qm->fID);
  gnuData(qm, buffer);
  snprintf(buffer,100, "wvsnotransfer_%d.txt", qm->fID);
  wvsData(qm, buffer);
#endif
  if (iV > 2) {
      ITMAX = 100;
      for (k  = 0; k <= 1; k++) {
#ifdef DEBUG
          if (k == 0) printf("Translating valences around mesh round 1\n ");
          else        printf("Translating valences around mesh round 2 :: forcing \n ");
#endif
          if (qm->regBd == 1) {
              for (q = 0; q < qm->totQ; q++) {
                  if (qm->qIdx[4 * q] == -2) continue; //can be a deleted quad
                  j = i = 0;
                  if (EG_regBoundQuad(qm, q + 1) != 0) {
                      for (i = 0; i < 4; i++) {
                          if(qm->valence[qm->qIdx[4 * q + i] - 1][2] == 3) {
                              j = 1;
                              break;
                          }
                      }
                      if (j == 1) {
                        i = qm->qIdx[4 * q + i];
                        if (EG_quadSize(qm, q + 1, NULL) <= 2.0 * qm->minArea) j = 1;
                        else j = 0;
                      }
                  } else {
                      j = -1;
                      stat = EG_swap(qm, q + 1, &activity);
                      if ( stat != EGADS_SUCCESS || activity > 0) {
#ifdef DEBUG
                          printf(" We have swapped quads  %d \n", q + 1);
#endif
                          continue;
                      }
                      if      (qm->vType  [qm->qIdx[4 * q    ] -1]    > 0 &&
                               qm->valence[qm->qIdx[4 * q    ] -1][1] > 4) j = 0;
                      else if (qm->vType  [qm->qIdx[4 * q + 1] -1]    > 0 &&
                               qm->valence[qm->qIdx[4 * q + 1] -1][1] > 4) j = 1;
                      else if (qm->vType  [qm->qIdx[4 * q + 2] -1] > 0 &&
                               qm->valence[qm->qIdx[4 * q + 2] -1][1] > 4) j = 2;
                      else if (qm->vType  [qm->qIdx[4 * q + 3] -1] > 0 &&
                               qm->valence[qm->qIdx[4 * q + 3] -1][1] > 4) j = 3;
                      if (j >= 0) {
                          if ((qm->qAdj[4 * q    ] == -1 || qm->qAdj[4 * q + 1] == -1 ||
                               qm->qAdj[4 * q + 2] == -1 || qm->qAdj[4 * q + 3] == -1 ) &&
                                  (qm->vType[qm->qIdx[4 * q    ]-1] != 2 &&
                                   qm->vType[qm->qIdx[4 * q + 1]-1] != 2 &&
                                   qm->vType[qm->qIdx[4 * q + 2]-1] != 2 &&
                                   qm->vType[qm->qIdx[4 * q + 3]-1] != 2)) continue;
                          i = qm->qIdx[4 * q + (j + 1)%4];
                          if (qm->vType[i - 1] != -1) i = qm->qIdx[4 * q + (j + 3)%4];
                          j = 1;
                      }
                  }
                  if (j == 1 && EG_validCollapse(qm, q + 1, i) == 1) {
#ifdef DEBUG
                      printf(" FORCE COLLAPSE: EITHER TOO SMALL OR IRREGULARBOUNDARY \n");
#endif
                      stat         = EG_collapse(qm, q + 1, &activity, 1, i);
                      totActivity += activity;
                      if (stat != EGADS_SUCCESS) {
                          printf(" EG_meshRegularization :: preprocess EG_collapse for quad %d -->%d!!\n ",
                                 q + 1, stat);
                          qm->fID = 0;
                          return stat;
                      }
                  }
              }
          }
          for ( it = 0 ; it < ITMAX; it++) {
              totActivity  =  0;
              prevPair[0]  = -1;
              prevPair[1]  = -1;
              for (q = 0; q < qm->totQ; q++) {
                  transfer = 0;
                  if (qm->qIdx[4 * q] == -2 || q + 1 == prevPair[0] ||
                          q + 1 == prevPair[1]) continue;
#ifdef DEBUG
                  printf("VERTEX TRANSLATE %d \n", q + 1);
#endif
                  qPair[0]  = q + 1;
                  qPair[1]  =    -1;
                  stat      = EG_transferValences(qm, qPair, k, &transfer, &activity);
                  if (stat     != EGADS_SUCCESS) {
                      printf("EG_meshRegularization:: EG_transferValences %d\n", stat);
                      qm->fID = 0;
                      return stat;
                  }
                  totActivity += activity;
                  if (activity == 0 || qPair[0] < 0 ||
                          qm->qIdx[4*(qPair[0]-1)] == -2) continue;
                  for (j = 0; j < 20; j++) {
                      activity = 0;
                      for (i  = 0; i < 2; i++) {
                          if (qPair[i] < 0 || qm->qIdx[4*(qPair[i]-1)] == -2) continue;
                          stat = EG_cleanNeighborhood(qm, qPair[i], transfer, &loopact);
                          if (stat != EGADS_SUCCESS) {
                              printf("EG_meshRegularization:: EG_cleanNeighborhood in transferValences %d\n", stat);
                              qm->fID = 0;
                              return stat;
                          }
                          activity +=loopact;
                      }
                      if (activity > 0) break;
                      stat      = EG_transferValences(qm, qPair, k,
                                                      &transfer, &activity);
                      if (stat != EGADS_SUCCESS) {
                          printf("EG_meshRegularization:: EG_transferValences after j loop %d\n", stat);
                          qm->fID = 0;
                          return stat;
                      }
                      if (activity == 0 || qPair[0] < 0 ||
                              qm->qIdx[4*(qPair[0]-1)] == -2) break;
                  }
                  if (iV <= 2) break;
                  prevPair[0] = qPair[0];
                  prevPair[1] = qPair[1];
              }
              if (iV <= 2) break;
              for (q = 0; q < qm->totQ; q++) {
                  if (qm->qIdx[4 * q] == -2) continue;
                  stat = EG_cleanQuad(qm, q + 1, 1, 0, 0, &i);
                  if (stat != EGADS_SUCCESS) {
                      printf(" In EG_meshRegularization clean quad %d !!\n ", stat);
                      qm->fID = 0;
                      return stat;
                  }
                  totActivity += i;
              }
              meshCount(qm, &iV, &totV, &vQ);
              if (iV <= 2 || totActivity == 0 ) break;
          }
      }
      if ( k == 1 ) qm->minArea = qm->avArea = 0.0;
  }
#ifdef DEBUG
  printf(" MESH BEFORE FINAL ROUND \n");
  gnuData(qm, NULL);
#endif
  for (q = 0; q < qm->totQ; q++) {
      if (qm->qIdx[4 * q] == -2) continue; //can be a deleted quad
      if ((qm->vType[qm->qIdx[4 * q    ]-1] * qm->valence[qm->qIdx[4 * q    ]-1][2] == -3 &&
           qm->vType[qm->qIdx[4 * q + 2]-1] * qm->valence[qm->qIdx[4 * q + 2]-1][2] == -3) ||
          (qm->vType[qm->qIdx[4 * q + 1]-1] * qm->valence[qm->qIdx[4 * q + 1]-1][2] == -3 &&
           qm->vType[qm->qIdx[4 * q + 3]-1] * qm->valence[qm->qIdx[4 * q + 3]-1][2] == -3) )
           stat = EG_collapse (qm, q + 1, &i, 1, -1);
      else stat = EG_cleanNeighborhood(qm, q + 1, 0, &activity);
      if (stat    != EGADS_SUCCESS) {
          printf(" In EG_cleanMesh: EG_CleanNeighborhood for quad %d -->%d!!\n ",
                 q + 1, stat);
          qm->fID = 0;
          return stat;
      }
      totActivity += activity;
  }
  for (q = 0; q < qm->totQ; q++) {
      if (qm->qIdx[4 * q] == -2) continue; //can be a deleted quad
      stat         = EG_cleanQuad(qm, q + 1, 0, 0, 0, &activity);
      if (stat    != EGADS_SUCCESS) {
          printf(" In EG_cleanMesh: EG_CleanNeighborhood for quad %d -->%d!!\n ",
                 q + 1, stat);
          qm->fID = 0;
          return stat;
      }
      totActivity += activity;
  }
  stat      = resizeQm(qm);
  if (stat != EGADS_SUCCESS) {
      printf("EG_meshRegularization final resize %d !!\n", stat);
      qm->fID = 0;
      return stat;
  }
#ifdef REPORT
  time1 = clock();
  ttime = (double)(time1 - time0) / CLOCKS_PER_SEC;
  min   = floor(ttime) / 60;
  secs  = ttime -  (double)min * 60.0;
  meshCount(qm, &iV, &totV, &vQ);
  printf("\n\n");
  printf("********************* REGULARIZATION RESULTS *********************\n");
  meshCount(qm, &iV, &totV, &vQ);
  printf(" Original mesh: %d quads and %d / %d irregular vertices (%.2f %%) \n",
         vQ0, iV0, totV0, (double) iV0 * 100.0 / (double)totV0);
  printf(" Final    mesh: %d quads and %d / %d irregular vertices (%.2f %%)\n",
         vQ, iV, totV, (double) iV * 100.0 / (double) totV);
  printf(" Invalid steps  %d\t CPU time: %d mins and %f sec\n",
         qm->invsteps, min, secs);
  printf("__________________________ IO FILES _______________________________\n");
  snprintf(buffer,100, "regularizationMeshStats_%d.txt", qm->fID);
  meshStats(qm, 10);
  snprintf(buffer,100, "wvsRegular_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100,"gnuRegular_%d.txt", qm->fID);
  gnuData(qm, buffer);
  printf("*******************************************************************\n");
#endif
  stat = EG_makeValidMesh(qm, 0, NULL, 1);
#ifdef REPORT
  printf("__________________________IO FILES ______________________________\n");
  snprintf(buffer,100, "gnuFin_%d.txt", qm->fID);
  gnuData(qm, buffer);
  snprintf(buffer,100, "wvsFin_%d.txt", qm->fID);
  wvsData(qm, buffer);
  snprintf(buffer,100, "finMeshStats_%d.txt", qm->fID);
  meshStats(qm, 100);
  meshProperties(qm,100);
  qm->invsteps = 0;
  printf("*****************************************************************\n");
  if (stat != EGADS_SUCCESS) printf("Deactivating face %d stat %d\n",
                                    qm->fID, stat);
#endif
  if (stat != EGADS_SUCCESS) qm->fID = 0;
  return  stat;
}


__HOST_AND_DEVICE__ int EG_makeQuadTess(bodyQuad bodydata, ego *quadTess)
{
  int          i, j, npts, nt, stat, outLevel, nedges, mQ, *tris;
  const int    *ptype, *pindex, *trs, *trc;
  const double *xyzs, *ts, *uvs;
  ego          tess, newTess, obj, *edges;
  egTessel     *btess, *ntess;
#ifndef LITE
  int          *iq;
#endif

  *quadTess = NULL;
  tess      = bodydata.tess;
  btess     = (egTessel *) tess->blind;
  outLevel  = EG_outLevel(tess);

  /* initialize the new tessellation object */
  stat = EG_initTessBody(btess->src, &newTess);
  if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_initTessBody = %d (EG_makeQuadTess)!\n", stat);
      return stat;
  }
  ntess = (egTessel *) newTess->blind;
  if (btess->src->oclass == EBODY) {
    stat  = EG_getBodyTopos(btess->src, NULL, EEDGE, &nedges, &edges);
  } else {
    stat  = EG_getBodyTopos(btess->src, NULL,  EDGE, &nedges, &edges);
  }
  if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_getBodyTopos E = %d (EG_makeQuadTess)!\n",
               stat);
      EG_deleteObject(newTess);
      return stat;
  }

  /* copy the Edges */
  for (i = 0; i < nedges; i++) {
      if (edges[i]->mtype == DEGENERATE) continue;
      stat = EG_getTessEdge(tess, i+1, &npts, &xyzs, &ts);
      if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_getTessEdge %d = %d (EG_makeQuadTess)!\n",
                   i+1, stat);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
      }
      if (npts == 0) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_getTessEdge %d -- no points (EG_makeQuadTess)!\n",
                   i+1);
          EG_free(edges);
          EG_deleteObject(newTess);
          return EGADS_INDEXERR;
      }

      stat = EG_setTessEdge(newTess, i+1, npts, xyzs, ts);
      if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_setTessEdge %d = %d (EG_makeQuadTess)!\n",
                   i+1, stat);
          EG_free(edges);
          EG_deleteObject(newTess);
          return stat;
      }
  }
  EG_free(edges);

  /* make the Face tessellations */
  for (mQ = i = 0; i < bodydata.nfaces; i++) {
      if (bodydata.qm[i]) {
        if (bodydata.qm[i]->totQ > mQ) mQ = bodydata.qm[i]->totQ;
        if ((btess->tess2d[i].tfi != 1) && (bodydata.qm[i]->fID != 0)) continue;
      }
      stat = EG_getTessFace(tess, i+1, &npts, &xyzs, &uvs, &ptype, &pindex,
                            &nt, &trs, &trc);
      if ((stat != EGADS_SUCCESS) || (nt == 0)) {
          if (outLevel > 0)
            if (stat != EGADS_SUCCESS) {
                printf(" EGADS Error: EG_getTessFace %d = %d (EG_quadTess)!\n",
                       i+1, stat);
            } else {
                printf(" EGADS Error: Face %d has no tessellation (EG_quadTess)!\n",
                       i+1);
            }
          EG_deleteObject(newTess);
          return stat;
      }
      stat = EG_setTessFace(newTess, i+1, npts, xyzs, uvs, nt, trs);
      if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            EG_deleteObject(newTess);
          return stat;
      }
      ntess->tess2d[i].tfi = 1;
  }
#ifdef DEBUG
  printf(" Max Quads = %d\n", mQ);
#endif
  tris = (int *) EG_alloc(6*mQ*sizeof(int));
  if (tris == NULL) {
      if (outLevel > 0)
        printf(" EGADS Error: Quad Allocation = %d (EG_makeQuadTess)!\n", mQ);
      EG_deleteObject(newTess);
      return EGADS_MALLOC;
  }

  /* fill in quads as triangle pairs */
  for (i = 0; i < bodydata.nfaces; i++) {
      if (bodydata.qm[i]        == NULL) continue;
      if ((btess->tess2d[i].tfi == 1) || (bodydata.qm[i]->fID == 0)) continue;


      for (j = 0; j < bodydata.qm[i]->totQ; j++) {
          tris[6*j  ] = bodydata.qm[i]->qIdx[4*j  ];
          tris[6*j+1] = bodydata.qm[i]->qIdx[4*j+1];
          tris[6*j+2] = bodydata.qm[i]->qIdx[4*j+2];
          tris[6*j+3] = bodydata.qm[i]->qIdx[4*j  ];
          tris[6*j+4] = bodydata.qm[i]->qIdx[4*j+2];
          tris[6*j+5] = bodydata.qm[i]->qIdx[4*j+3];
      }

      stat = EG_setTessFace(newTess, i+1, bodydata.qm[i]->totV,
                            bodydata.qm[i]->xyzs,
                            bodydata.qm[i]->uvs,
                            bodydata.qm[i]->totQ*2, tris);
      if (stat != EGADS_SUCCESS) {
          if (outLevel > 0)
            printf(" EGADS Error: EG_setTessFace %d = %d (EG_makeQuadTess)!\n",
                   i+1, stat);
          EG_free(tris);
          EG_deleteObject(newTess);
          return stat;
      }
  }
  EG_free(tris);

  /* close up the open tessellation */
  stat = EG_statusTessBody(newTess, &obj, &i, &npts);
  if (stat != EGADS_SUCCESS) {
      if (outLevel > 0)
        printf(" EGADS Error: EG_statusTessBody = %d (EG_makeQuadTess)!\n", stat);
      EG_deleteObject(newTess);
      return stat;
  }
  if (i != 1) {
      if (outLevel > 0)
        printf(" EGADS Error: New Tessellation Object is Open (EG_makeQuadTess)!\n");
      EG_deleteObject(newTess);
      return EGADS_TESSTATE;
  }

#ifndef LITE
  iq = (int *) EG_alloc(bodydata.nfaces*sizeof(int));
  if (iq != NULL) {
    for (i = 0; i < bodydata.nfaces; i++) iq[i] = ntess->tess2d[i].ntris/2;
    stat = EG_attributeAdd(newTess, ".mixed", ATTRINT, bodydata.nfaces,
                           iq, NULL, NULL);
    if (stat != EGADS_SUCCESS)
      printf(" EGADS Warning: EG_attributeAdd M = %d (EG_makeQuadTess)!\n",
             stat);
    EG_free(iq);
  }
  stat = EG_attributeAdd(newTess, ".tessType", ATTRSTRING, 4,
                         NULL, NULL, "Quad");
  if (stat != EGADS_SUCCESS)
    if (outLevel > 0)
      printf(" EGADS Warning: EG_attributeAdd Q = %d (EG_makeQuadTess)!\n",
             stat);
#endif

  *quadTess = newTess;
  return EGADS_SUCCESS;
}


#ifdef STANDALONE
int main (int argc, char *argv[])
{
  int          stat = 0,  f , i, j, iBody, oclass, mtype, nbody;
  int          atype, alen, *senses,FACECHOICE = -1;
  const int    *ints;
  float        arg;
  double       box[6], size, params[3];
  const double *reals;
  const char   *OCCrev, *string;
  ego          context, tess, model, geom, *bodies, *dum, newTess;
  bodyQuad     *bodydata;
#ifdef REPORT
  clock_t      start_t, end_t, total_t;
  double       time, frac;
  int min;
  start_t = clock();
#endif

  if (argc != 3 && argc != 2 && argc != 5 && argc != 6 && argc != 7) {
      printf("\n Usage: %s = (1) filename, [Tesselation:  (2) angle (3) maxlen and chord (4)] [facechoice] \n\n", argv[0]);
      printf(" You have entered : \n");
      for (i = 0; i < argc; i++) printf(" argv[%d] = %s\n ", i, argv[i]);
      return 1;
  }
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  stat = EG_open(&context);
  if (stat != EGADS_SUCCESS) return 1;
  stat = EG_loadModel(context, 0, argv[1], &model);
  printf(" EG_loadModel      = %d\n", stat);
  printf(" EG_getBoundingBox = %d\n", EG_getBoundingBox(model, box));
  printf("       BoundingBox = %lf %lf %lf\n", box[0], box[1], box[2]);
  printf("                     %lf %lf %lf\n", box[3], box[4], box[5]);
  printf(" \n");
  size = box[3]-box[0];
  if (size < box[4]-box[1]) size = box[4]-box[1];
  if (size < box[5]-box[2]) size = box[5]-box[2];


  if (argc > 3) {
      sscanf(argv[2], "%f", &arg);
      params[2] = arg;
      sscanf(argv[3], "%f", &arg);
      params[0] = arg;
      sscanf(argv[4], "%f", &arg);
      params[1] = arg;
      printf(" Using angle = %lf,  relSide = %lf,  relSag = %lf\n",
             params[2], params[0], params[1]);

  } else {
    params[0] =  0.025;
    params[1] =  0.001;
    params[2] = 15.0;
  }
  params[0] *= size;
  params[1] *= size;
  if (argc == 3) sscanf(argv[2], "%d", &FACECHOICE);
  if (argc >= 6) sscanf(argv[5], "%d", &FACECHOICE);
  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
                        &bodies, &senses);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_getTopology = %d\n", stat);
      return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  bodydata = (bodyQuad *) malloc(nbody*sizeof(bodyQuad));
  if (bodydata == NULL) {
      printf(" MALLOC Error on Body storage!\n");
      return 1;
  }

  /* fill our structure a body at at time */
  for (iBody = 0; iBody < nbody; iBody++) {
      /* FORCE QUADS FOR DEBUG MODE: IF WANTED uncomment the 4 lines */
    /*stat = EG_attributeAdd(bodies[iBody], ".qParams",
                           ATTRSTRING, 4, NULL, NULL, "off");
    if (stat != EGADS_SUCCESS)
        printf(" Body %d: attributeAdd = %d\n", iBody, stat);*/
      EG_getTopology(bodies[iBody], &geom, &oclass,
                     &mtype, NULL, &j, &dum, &senses);
      if (mtype == WIREBODY) {
          printf(" Body %d: Type = WireBody\n",  iBody+1);
      } else if (mtype == FACEBODY) {
          printf(" Body %d: Type = FaceBody\n",  iBody+1);
      } else if (mtype == SHEETBODY) {
          printf(" Body %d: Type = SheetBody\n", iBody+1);
      } else {
          printf(" Body %d: Type = SolidBody\n", iBody+1);
      }
      stat = EG_getBodyTopos(bodies[iBody], NULL, FACE,
                             &bodydata[iBody].nfaces, &bodydata[iBody].faces);
      i    = EG_getBodyTopos(bodies[iBody], NULL, EDGE,
                             &bodydata[iBody].nedges, NULL);
      if ((stat != EGADS_SUCCESS) || (i != EGADS_SUCCESS)) {
          printf(" EG_getBodyTopos Face = %d\n", stat);
          printf(" EG_getBodyTopos Edge = %d\n", i);
          continue;
      }
      stat = EG_makeTessBody(bodies[iBody], params, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
          printf(" EG_makeTessBody %d = %d\n", iBody, stat);
          continue;
      }
      tess = bodydata[iBody].tess;
      /* disable regularization in EGADS */
      stat = EG_attributeAdd(tess, ".qRegular", ATTRSTRING, 3, NULL, NULL, "Off");
      if (stat != EGADS_SUCCESS)
        printf(" EG_attributeAdd qRegular %d = %d\n", iBody, stat);
      stat = EG_quadTess(tess, &bodydata[iBody].tess);
      if (stat != EGADS_SUCCESS) {
          printf(" EG_quadTess %d = %d  -- reverting...\n", iBody, stat);
          bodydata[iBody].tess = tess;
          continue;
      }
      EG_deleteObject(tess);
  }
  for (iBody = 0; iBody < nbody; iBody++) {
      stat = EG_attributeRet(bodydata[iBody].tess, ".tessType", &atype,
                             &alen, &ints, &reals, &string);
      if (stat != EGADS_SUCCESS) {
          printf(" Tessellation is NOT Quadded!\n");
          continue;
      }
      if (atype != ATTRSTRING) {
          printf(" Tessellation Flag is the Wrong Type!\n");
          continue;
      }
      if (strcmp(string, "Quad") != 0) {
          printf(" Tessellation Flag is Not Quad = %s!\n", string);
          continue;
      }
      if (FACECHOICE > bodydata[iBody].nfaces) FACECHOICE = -1;
      if (EG_createMeshMap(&bodydata[iBody]) == EGADS_SUCCESS) {
        for (f = 0; f < bodydata[iBody].nfaces; ++f) {
            if (FACECHOICE >= 0) f = FACECHOICE - 1;
            stat      = EG_meshRegularization(bodydata[iBody].qm[f]);
            if (stat != EGADS_SUCCESS) {
#ifndef __CUDA_ARCH__
                fprintf(stderr, " EG_meshRegularization face %d / %d = %d DEACTIVATE !!!\n ",
                        f + 1, bodydata[iBody].nfaces,  stat);
#else
                printf(" EG_meshRegularization face %d / %d = %d DEACTIVATE !!!\n ",
                        f + 1, bodydata[iBody].nfaces,  stat);
#endif
                bodydata[iBody].qm[f]->fID = 0;
            }
            if (FACECHOICE >= 0) break;
        }
      }
      if (FACECHOICE == -1) {
          stat = EG_makeQuadTess(bodydata[iBody], &newTess);
          if (stat != EGADS_SUCCESS) {
#ifndef __CUDA_ARCH__
              fprintf(stderr, " EG_makeQuadTess = %d!\n ", stat);
#else
              printf(" EG_makeQuadTess = %d!\n ", stat);
#endif
              goto cleanup;
          }
          stat = EG_deleteObject(newTess);
          if (stat != EGADS_SUCCESS) {
#ifndef __CUDA_ARCH__
              fprintf(stderr, " EG_deleteObject = %d!\n ", stat);
#else
              printf(" EG_deleteObject = %d!\n ", stat);
#endif
              goto cleanup;
          }
      }
  }
  cleanup:
  if (stat != EGADS_SUCCESS ) printf(" EG_main stat in cleanup %d !!\n", stat);
  for (iBody = 0; iBody < nbody; iBody++) {
      EG_destroyMeshMap(&bodydata[iBody]);
      EG_free(bodydata[iBody].faces);
      EG_deleteObject(bodydata[iBody].tess);
  }
  EG_free(bodydata);
#ifdef REPORT
  end_t   = clock();
  total_t = end_t - start_t;
  time    = (double) total_t / CLOCKS_PER_SEC;
  min     = floor(time) / 60;
  frac    = time -  min * 60;
#ifndef __CUDA_ARCH__
  fprintf(stderr, "Total time taken by CPU: %d minutes and %f seconds\n",
          min, frac);
#else
  printf("Total time taken by CPU: %d minutes and %f seconds\n",
         min, frac);
#endif
#endif
  EG_deleteObject(model);
  EG_close(context);
  return 0;
}
#endif
