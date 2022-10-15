#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include <petsc/private/hashmapi.h>

#ifdef PETSC_HAVE_EGADS
#include <egads_lite.h>

PETSC_INTERN PetscErrorCode DMPlex_EGADSlite_EDGE_XYZtoUV_Internal(const PetscScalar[], ego, const PetscScalar[], const PetscInt, const PetscInt, PetscScalar[]);
PETSC_INTERN PetscErrorCode DMPlex_EGADSlite_FACE_XYZtoUV_Internal(const PetscScalar[], ego, const PetscScalar[], const PetscInt, const PetscInt, PetscScalar[]);

PetscErrorCode DMPlex_EGADSlite_EDGE_XYZtoUV_Internal(const PetscScalar coords[], ego obj, const PetscScalar range[], const PetscInt v, const PetscInt dE, PetscScalar paramsV[])
{
  PetscInt       loopCntr = 0;
  PetscScalar    dx, dy, dz, lambda, tolr, obj_old, obj_tmp;
  PetscScalar    delta, A, b;
  PetscScalar    ts[2], tt[2], eval[18], data[18];
  //PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  ts[0] = (range[0] + range[1]) / 2.;
  
  while (tolr >= 1.0e-12) {
    PetscCall(EGlite_evaluate(obj, ts, eval));
    dx = coords[v*dE+0] - eval[0];
    dy = coords[v*dE+1] - eval[1];
    dz = coords[v*dE+2] - eval[2];
    obj_old = dx*dx + dy*dy + dz*dz;
	
    if (obj_old < 1.0E-14) {tolr = obj_old; break;}
	
    A = (eval[3]*eval[3] + eval[4]*eval[4] + eval[5]*eval[5]) * (1.0 + lambda);
    if (A == 0.0) {PetscCall(PetscPrintf(PETSC_COMM_SELF, "A = 0.0 \n")); break;}
    b = eval[3]*dx + eval[4]*dy + eval[5]*dz;    

    /* Solve A*delta = b */
    delta = b/A;
	
    /* Find a temp (u,v) and associated objective function */
    tt[0] = ts[0] + delta;
    if (tt[0] < range[0]) {tt[0] = range[0]; delta = tt[0] - ts[0];}
    if (tt[0] > range[1]) {tt[0] = range[1]; delta = tt[0] - ts[0];}
	
    PetscCall(EGlite_evaluate(obj, tt, data));
	
    obj_tmp = (coords[v*dE+0] - data[0])*(coords[v*dE+0] - data[0]) +
              (coords[v*dE+1] - data[1])*(coords[v*dE+1] - data[1]) +
              (coords[v*dE+2] - data[2])*(coords[v*dE+2] - data[2]);
	
    /* If step is better, accept it and halve lambda (making it more Newton-like) */
    if (obj_tmp < obj_old) {
      obj_old = obj_tmp;
      ts[0] = tt[0];
      for (int jj = 0; jj < 18; ++jj) eval[jj] = data[jj];
      lambda /= 2.0;
      if (lambda < 1.0E-14) lambda = 1.0E-14;
      if (obj_old < 1.0E-14) {tolr = obj_old; break;}
    } else {
      /* Otherwise reject it and double lambda (making it more gradient-descent like) */
      lambda *= 2.0;
    }
	
    if ((tt[0] == range[0]) || (tt[0] == range[1])) break;
    if (fabs(delta) < 1.0E-12) {tolr = obj_old; break;}
	
    tolr = obj_old;
	
    loopCntr += 1;
    if (loopCntr > 30) break;
  }
  paramsV[v*3+0] = ts[0];
  paramsV[v*3+1] = 0.;
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlex_EGADSlite_FACE_XYZtoUV_Internal(const PetscScalar coords[], ego obj, const PetscScalar range[], const PetscInt v, const PetscInt dE, PetscScalar paramsV[])
{
  PetscInt       loopCntr = 0;
  PetscScalar    dx, dy, dz, lambda, tolr, denom, obj_old, obj_tmp;
  PetscScalar    uvs[2], uvt[2], delta[2], A[4], b[2], eval[18], data[18];
  //PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  uvs[0] = (range[0] + range[1]) / 2.;
  uvs[1] = (range[2] + range[3]) / 2.;
  
  while (tolr >= 1.0e-12) {
    PetscCall(EGlite_evaluate(obj, uvs, eval));
    dx = coords[v*dE+0] - eval[0];
    dy = coords[v*dE+1] - eval[1];
    dz = coords[v*dE+2] - eval[2];
    obj_old = dx*dx + dy*dy + dz*dz;
	
    if (obj_old < 1.0E-14) {tolr = obj_old; break;}
	
    A[0] = (eval[3]*eval[3] + eval[4]*eval[4] + eval[5]*eval[5]) * (1.0 + lambda);
    A[1] =  eval[3]*eval[6] + eval[4]*eval[7] + eval[5]*eval[8];
    A[2] = A[1];
    A[3] = (eval[6]*eval[6] + eval[7]*eval[7] + eval[8]*eval[8]) * (1.0 + lambda);
	
    b[0] = eval[3]*dx + eval[4]*dy + eval[5]*dz;
    b[1] = eval[6]*dx + eval[7]*dy + eval[8]*dz;
	
    /* Solve A*delta = b using Cramer's Rule */
    denom = A[0]*A[3] - A[2]*A[1];
    if (denom == 0.0) {PetscCall(PetscPrintf(PETSC_COMM_SELF, "denom = 0.0 \n"));}
    delta[0] = (b[0]*A[3] - b[1]*A[1]) / denom;
    delta[1] = (A[0]*b[1] - A[2]*b[0]) / denom;
	
    /* Find a temp (u,v) and associated objective function */
    uvt[0] = uvs[0] + delta[0];
    uvt[1] = uvs[1] + delta[1];
    
    if (uvt[0] < range[0]) {uvt[0] = range[0]; delta[0] = uvt[0] - uvs[0];}
    if (uvt[0] > range[1]) {uvt[0] = range[1]; delta[0] = uvt[0] - uvs[0];}
    if (uvt[1] < range[2]) {uvt[1] = range[2]; delta[1] = uvt[1] - uvs[1];}
    if (uvt[1] > range[3]) {uvt[1] = range[3]; delta[1] = uvt[1] - uvs[1];}
    
    PetscCall(EGlite_evaluate(obj, uvt, data));
    
    obj_tmp = (coords[v*dE+0] - data[0])*(coords[v*dE+0] - data[0]) +
              (coords[v*dE+1] - data[1])*(coords[v*dE+1] - data[1]) +
              (coords[v*dE+2] - data[2])*(coords[v*dE+2] - data[2]);
	
    /* If step is better, accept it and halve lambda (making it more Newton-like) */
    if (obj_tmp < obj_old) {
      obj_old = obj_tmp;
      uvs[0] = uvt[0];
      uvs[1] = uvt[1];
      for (int jj = 0; jj < 18; ++jj) eval[jj] = data[jj];
      lambda /= 2.0;
      if (lambda < 1.0E-14) lambda = 1.0E-14;
      if (obj_old < 1.0E-14) {tolr = obj_old; break;}
    } else {
      /* Otherwise reject it and double lambda (making it more gradient-descent like) */
      lambda *= 2.0;
    }
	
    if (sqrt(delta[0]*delta[0] + delta[1]*delta[1]) < 1.0E-12) {tolr = obj_old; break;}
	
    tolr = obj_old;
	
    loopCntr += 1;
    if (loopCntr > 30) break;
  }
  paramsV[v*3+0] = uvs[0];
  paramsV[v*3+1] = uvs[1];
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlexSnapToGeomModel_EGADSLite_Internal(DM dm, PetscInt p, ego model, PetscInt bodyID, PetscInt faceID, PetscInt edgeID, const PetscScalar mcoords[], PetscScalar gcoords[])
{
  DM             cdm;
  ego           *bodies;
  ego            geom, body, obj;
  /* result has to hold derviatives, along with the value */
  double         params[3], result[18], paramsV[16*3], range[4];
  int            Nb, oclass, mtype, *senses, peri;
  Vec            coordinatesLocal;
  PetscScalar   *coords = NULL;
  PetscInt       Nv, v, Np = 0, pm;
  PetscInt       dE, d;
  PetscReal      pTolr = 1.0e-14;
  //PetscErrorCode ierr;

  PetscFunctionBeginHot;
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetCoordinateDim(dm, &dE));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinatesLocal));
  PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  if (bodyID >= Nb) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];

  if      (edgeID >= 0) {PetscCall(EGlite_objectBodyTopo(body, EDGE, edgeID, &obj)); Np = 1;}
  else if (faceID >= 0) {PetscCall(EGlite_objectBodyTopo(body, FACE, faceID, &obj)); Np = 2;}
  else {
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }

  /* Calculate parameters (t or u,v) for vertices */
  PetscCall(DMPlexVecGetClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords));
  Nv  /= dE;
  if (Nv == 1) {
    PetscCall(DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords));
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }
  if (Nv > 16) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Cannot handle %D coordinates associated to point %D", Nv, p);

  /* Correct EGADSlite 2pi bug when calculating nearest point on Periodic Surfaces */
  PetscCall(EGlite_getRange(obj, range, &peri));

  for (v = 0; v < Nv; ++v) {
    if (edgeID > 0) {PetscCall(DMPlex_EGADSlite_EDGE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV));}
    else            {PetscCall(DMPlex_EGADSlite_FACE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV));}
  }
  PetscCall(DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords));
  /* Calculate parameters (t or u,v) for new vertex at edge midpoint */
  for (pm = 0; pm < Np; ++pm) {
    params[pm] = 0.;
    for (v = 0; v < Nv; ++v) params[pm] += paramsV[v*3+pm];
    params[pm] /= Nv;
  }
  if ((params[0] + pTolr < range[0]) || (params[0] - pTolr > range[1])) {SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on t or u", p);}
  if (Np > 1 && ((params[1] + pTolr < range[2]) || (params[1] - pTolr > range[3]))) {SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on v", p);}
  /* Put coordinates for new vertex in result[] */
  PetscCall(EGlite_evaluate(obj, params, result));
  for (d = 0; d < dE; ++d) gcoords[d] = result[d];
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexEGADSliteDestroy_Private(void *context)
{
  if (context) EGlite_close((ego) context);
  return 0;
}

static PetscErrorCode DMPlexCreateEGADSlite_Internal(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  PetscInt       cStart, cEnd, c;
  /* EGADSLite variables */
  ego            geom, *bodies, *objs, *nobjs, *mobjs, *lobjs;
  int            oclass, mtype, nbodies, *senses;
  int            b;
  /* PETSc variables */
  DM             dm;
  PetscHMapI     edgeMap = NULL;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, maxCorners = 0, numVertices = 0, newVertices = 0, numEdges = 0, numCells = 0, newCells = 0, numQuads = 0, cOff = 0, fOff = 0;
  PetscInt      *cells  = NULL, *cone = NULL;
  PetscReal     *coords = NULL;
  PetscMPIInt    rank;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscCall(MPI_Comm_rank(comm, &rank));
  if (!rank) {
    const PetscInt debug = 0;

    /* ---------------------------------------------------------------------------------------------------
    Generate Petsc Plex
      Get all Nodes in model, record coordinates in a correctly formatted array
      Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
      We need to uniformly refine the initial geometry to guarantee a valid mesh
    */

    /* Calculate cell and vertex sizes */
    PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));
    PetscCall(PetscHMapICreate(&edgeMap));
    numEdges = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l, Nv, v;

      PetscCall(EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int Ner  = 0, Ne, e, Nc;

        PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int Nv, id;
          PetscHashIter iter;
          PetscBool     found;

          PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          id   = EGlite_indexBodyTopo(body, edge);
          PetscCall(PetscHMapIFind(edgeMap, id-1, &iter, &found));
          if (!found) {PetscCall(PetscHMapISet(edgeMap, id-1, numEdges++));}
          ++Ner;
        }
        if (Ner == 2)      {Nc = 2;}
        else if (Ner == 3) {Nc = 4;}
        else if (Ner == 4) {Nc = 8; ++numQuads;}
        else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Cannot support loop with %d edges", Ner);
        numCells += Nc;
        newCells += Nc-1;
        maxCorners = PetscMax(Ner*2+1, maxCorners);
      }
      PetscCall(EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      for (v = 0; v < Nv; ++v) {
        ego vertex = nobjs[v];

        id = EGlite_indexBodyTopo(body, vertex);
        /* TODO: Instead of assuming contiguous ids, we could use a hash table */
        numVertices = PetscMax(id, numVertices);
      }
      EGlite_free(lobjs);
      EGlite_free(nobjs);
    }
    PetscCall(PetscHMapIGetSize(edgeMap, &numEdges));
    newVertices  = numEdges + numQuads;
    numVertices += newVertices;

    dim        = 2; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    cdim       = 3; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    numCorners = 3; /* Split cells into triangles */
    PetscCall(PetscMalloc3(numVertices*cdim, &coords, numCells*numCorners, &cells, maxCorners, &cone));

    /* Get vertex coordinates */
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nv, v;

      PetscCall(EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      for (v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    dummy;

        PetscCall(EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
        id   = EGlite_indexBodyTopo(body, vertex);
        coords[(id-1)*cdim+0] = limits[0];
        coords[(id-1)*cdim+1] = limits[1];
        coords[(id-1)*cdim+2] = limits[2];
      }
      EGlite_free(nobjs);
    }
    PetscCall(PetscHMapIClear(edgeMap));
    fOff     = numVertices - newVertices + numEdges;
    numEdges = 0;
    numQuads = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int Nl, l;

      PetscCall(EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e;

        lid  = EGlite_indexBodyTopo(body, loop);
        PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
        for (e = 0; e < Ne; ++e) {
          ego       edge = objs[e];
          int       eid, Nv;
          PetscHashIter iter;
          PetscBool     found;

          PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          ++Ner;
          eid  = EGlite_indexBodyTopo(body, edge);
          PetscCall(PetscHMapIFind(edgeMap, eid-1, &iter, &found));
          if (!found) {
            PetscInt v = numVertices - newVertices + numEdges;
            double range[4], params[3] = {0., 0., 0.}, result[18];
            int    periodic[2];

            PetscCall(PetscHMapISet(edgeMap, eid-1, numEdges++));
            PetscCall(EGlite_getRange(edge, range, periodic));
            params[0] = 0.5*(range[0] + range[1]);
            PetscCall(EGlite_evaluate(edge, params, result));
            coords[v*cdim+0] = result[0];
            coords[v*cdim+1] = result[1];
            coords[v*cdim+2] = result[2];
          }
        }
        if (Ner == 4) {
          PetscInt v = fOff + numQuads++;
          ego     *fobjs, face;
          double   range[4], params[3] = {0., 0., 0.}, result[18];
          int      Nf, fid, periodic[2];

          PetscCall(EGlite_getBodyTopos(body, loop, FACE, &Nf, &fobjs));
          face = fobjs[0];
          fid  = EGlite_indexBodyTopo(body, face);
          if (Nf != 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Loop %d has %d faces, instead of 1 (%d)", lid-1, Nf, fid);
          PetscCall(EGlite_getRange(face, range, periodic));
          params[0] = 0.5*(range[0] + range[1]);
          params[1] = 0.5*(range[2] + range[3]);
          PetscCall(EGlite_evaluate(face, params, result));
          coords[v*cdim+0] = result[0];
          coords[v*cdim+1] = result[1];
          coords[v*cdim+2] = result[2];
        }
      }
    }
    if (numEdges + numQuads != newVertices) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Number of new vertices %D != %D previous count", numEdges + numQuads, newVertices);

    /* Get cell vertices by traversing loops */
    numQuads = 0;
    cOff     = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l;

      PetscCall(EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e, nc = 0, c, Nt, t;

        lid  = EGlite_indexBodyTopo(body, loop);
        PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int points[3];
          int eid, Nv, v, tmp;

          eid  = EGlite_indexBodyTopo(body, edge);
          PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          else                     ++Ner;
          if (Nv != 2) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Edge %d has %d vertices != 2", eid, Nv);

          for (v = 0; v < Nv; ++v) {
            ego vertex = nobjs[v];

            id = EGlite_indexBodyTopo(body, vertex);
            points[v*2] = id-1;
          }
          {
            PetscInt edgeNum;

            PetscCall(PetscHMapIGet(edgeMap, eid-1, &edgeNum));
            points[1] = numVertices - newVertices + edgeNum;
          }
          /* EGADS loops are not oriented, but seem to be in order, so we must piece them together */
          if (!nc) {
            for (v = 0; v < Nv+1; ++v) cone[nc++] = points[v];
          } else {
            if (cone[nc-1] == points[0])      {cone[nc++] = points[1]; if (cone[0] != points[2]) cone[nc++] = points[2];}
            else if (cone[nc-1] == points[2]) {cone[nc++] = points[1]; if (cone[0] != points[0]) cone[nc++] = points[0];}
            else if (cone[nc-3] == points[0]) {tmp = cone[nc-3]; cone[nc-3] = cone[nc-1]; cone[nc-1] = tmp; cone[nc++] = points[1]; if (cone[0] != points[2]) cone[nc++] = points[2];}
            else if (cone[nc-3] == points[2]) {tmp = cone[nc-3]; cone[nc-3] = cone[nc-1]; cone[nc-1] = tmp; cone[nc++] = points[1]; if (cone[0] != points[0]) cone[nc++] = points[0];}
            else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Edge %d does not match its predecessor", eid);
          }
        }
        if (nc != 2*Ner)     SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Number of corners %D != %D", nc, 2*Ner);
        if (Ner == 4) {cone[nc++] = numVertices - newVertices + numEdges + numQuads++;}
        if (nc > maxCorners) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Number of corners %D > %D max", nc, maxCorners);
        /* Triangulate the loop */
        switch (Ner) {
          case 2: /* Bi-Segment -> 2 triangles */
            Nt = 2;
            cells[cOff*numCorners+0] = cone[0];
            cells[cOff*numCorners+1] = cone[1];
            cells[cOff*numCorners+2] = cone[2];
            ++cOff;
            cells[cOff*numCorners+0] = cone[0];
            cells[cOff*numCorners+1] = cone[2];
            cells[cOff*numCorners+2] = cone[3];
            ++cOff;
            break;
          case 3: /* Triangle   -> 4 triangles */
            Nt = 4;
            cells[cOff*numCorners+0] = cone[0];
            cells[cOff*numCorners+1] = cone[1];
            cells[cOff*numCorners+2] = cone[5];
            ++cOff;
            cells[cOff*numCorners+0] = cone[1];
            cells[cOff*numCorners+1] = cone[2];
            cells[cOff*numCorners+2] = cone[3];
            ++cOff;
            cells[cOff*numCorners+0] = cone[5];
            cells[cOff*numCorners+1] = cone[3];
            cells[cOff*numCorners+2] = cone[4];
            ++cOff;
            cells[cOff*numCorners+0] = cone[1];
            cells[cOff*numCorners+1] = cone[3];
            cells[cOff*numCorners+2] = cone[5];
            ++cOff;
            break;
          case 4: /* Quad       -> 8 triangles */
            Nt = 8;
            cells[cOff*numCorners+0] = cone[0];
            cells[cOff*numCorners+1] = cone[1];
            cells[cOff*numCorners+2] = cone[7];
            ++cOff;
            cells[cOff*numCorners+0] = cone[1];
            cells[cOff*numCorners+1] = cone[2];
            cells[cOff*numCorners+2] = cone[3];
            ++cOff;
            cells[cOff*numCorners+0] = cone[3];
            cells[cOff*numCorners+1] = cone[4];
            cells[cOff*numCorners+2] = cone[5];
            ++cOff;
            cells[cOff*numCorners+0] = cone[5];
            cells[cOff*numCorners+1] = cone[6];
            cells[cOff*numCorners+2] = cone[7];
            ++cOff;
            cells[cOff*numCorners+0] = cone[8];
            cells[cOff*numCorners+1] = cone[1];
            cells[cOff*numCorners+2] = cone[3];
            ++cOff;
            cells[cOff*numCorners+0] = cone[8];
            cells[cOff*numCorners+1] = cone[3];
            cells[cOff*numCorners+2] = cone[5];
            ++cOff;
            cells[cOff*numCorners+0] = cone[8];
            cells[cOff*numCorners+1] = cone[5];
            cells[cOff*numCorners+2] = cone[7];
            ++cOff;
            cells[cOff*numCorners+0] = cone[8];
            cells[cOff*numCorners+1] = cone[7];
            cells[cOff*numCorners+2] = cone[1];
            ++cOff;
            break;
          default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d edges, which we do not support", lid, Ner);
        }
        if (debug) {
          for (t = 0; t < Nt; ++t) {
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "  LOOP Corner NODEs Triangle %D (", t));
            for (c = 0; c < numCorners; ++c) {
              if (c > 0) {PetscCall(PetscPrintf(PETSC_COMM_SELF, ", "));}
              PetscCall(PetscPrintf(PETSC_COMM_SELF, "%D", cells[(cOff-Nt+t)*numCorners+c]));
            }
            PetscCall(PetscPrintf(PETSC_COMM_SELF, ")\n"));
          }
        }
      }
      EGlite_free(lobjs);
    }
  }
  if (cOff != numCells) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Count of total cells %D != %D previous count", cOff, numCells);
  PetscCall(DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numVertices, numCorners, PETSC_TRUE, cells, cdim, coords, &dm));
  PetscCall(PetscFree3(coords, cells, cone));
  PetscCall(PetscInfo(dm, " Total Number of Unique Cells    = %D (%D)\n", numCells, newCells));
  PetscCall(PetscInfo(dm, " Total Number of Unique Vertices = %D (%D)\n", numVertices, newVertices));
  /* Embed EGADS model in DM */
  {
    PetscContainer modelObj, contextObj;

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &modelObj));
    PetscCall(PetscContainerSetPointer(modelObj, model));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Context", (PetscObject) contextObj));
    PetscCall(PetscContainerDestroy(&contextObj));
  }
  /* Label points */
  PetscCall(DMCreateLabel(dm, "EGADS Body ID"));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Face ID"));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Edge ID"));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Vertex ID"));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));
  cOff = 0;
  for (b = 0; b < nbodies; ++b) {
    ego body = bodies[b];
    int id, Nl, l;

    PetscCall(EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
    for (l = 0; l < Nl; ++l) {
      ego  loop = lobjs[l];
      ego *fobjs;
      int  lid, Nf, fid, Ner = 0, Ne, e, Nt = 0, t;

      lid  = EGlite_indexBodyTopo(body, loop);
      PetscCall(EGlite_getBodyTopos(body, loop, FACE, &Nf, &fobjs));
      if (Nf > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);
      fid  = EGlite_indexBodyTopo(body, fobjs[0]);
      EGlite_free(fobjs);
      PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
      for (e = 0; e < Ne; ++e) {
        ego             edge = objs[e];
        int             eid, Nv, v;
        PetscInt        points[3], support[2], numEdges, edgeNum;
        const PetscInt *edges;

        eid  = EGlite_indexBodyTopo(body, edge);
        PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
        if (mtype == DEGENERATE) continue;
        else                     ++Ner;
        for (v = 0; v < Nv; ++v) {
          ego vertex = nobjs[v];

          id   = EGlite_indexBodyTopo(body, vertex);
          PetscCall(DMLabelSetValue(edgeLabel, numCells + id-1, eid));
          points[v*2] = numCells + id-1;
        }
        PetscCall(PetscHMapIGet(edgeMap, eid-1, &edgeNum));
        points[1] = numCells + numVertices - newVertices + edgeNum;

        PetscCall(DMLabelSetValue(edgeLabel, points[1], eid));
        support[0] = points[0];
        support[1] = points[1];
        PetscCall(DMPlexGetJoin(dm, 2, support, &numEdges, &edges));
        if (numEdges != 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Vertices (%D, %D) should only bound 1 edge, not %D", support[0], support[1], numEdges);
        PetscCall(DMLabelSetValue(edgeLabel, edges[0], eid));
        PetscCall(DMPlexRestoreJoin(dm, 2, support, &numEdges, &edges));
        support[0] = points[1];
        support[1] = points[2];
        PetscCall(DMPlexGetJoin(dm, 2, support, &numEdges, &edges));
        if (numEdges != 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Vertices (%D, %D) should only bound 1 edge, not %D", support[0], support[1], numEdges);
        PetscCall(DMLabelSetValue(edgeLabel, edges[0], eid));
        PetscCall(DMPlexRestoreJoin(dm, 2, support, &numEdges, &edges));
      }
      switch (Ner) {
        case 2: Nt = 2;break;
        case 3: Nt = 4;break;
        case 4: Nt = 8;break;
        default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Loop with %d edges is unsupported", Ner);
      }
      for (t = 0; t < Nt; ++t) {
        PetscCall(DMLabelSetValue(bodyLabel, cOff+t, b));
        PetscCall(DMLabelSetValue(faceLabel, cOff+t, fid));
      }
      cOff += Nt;
    }
    EGlite_free(lobjs);
  }
  PetscCall(PetscHMapIDestroy(&edgeMap));
  PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
  for (c = cStart; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt  clSize, cl, bval, fval;

    PetscCall(DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &clSize, &closure));
    PetscCall(DMLabelGetValue(bodyLabel, c, &bval));
    PetscCall(DMLabelGetValue(faceLabel, c, &fval));
    for (cl = 0; cl < clSize*2; cl += 2) {
      PetscCall(DMLabelSetValue(bodyLabel, closure[cl], bval));
      PetscCall(DMLabelSetValue(faceLabel, closure[cl], fval));
    }
    PetscCall(DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &clSize, &closure));
  }
  *newdm = dm;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexEGADSLitePrintModel_Internal(ego model)
{
  ego            geom, *bodies, *nobjs, *mobjs, *lobjs, *shobjs, *fobjs, *eobjs;
  int            oclass, mtype, *senses, *shsenses, *fsenses, *lsenses, *esenses;
  int            Nb, b;
  //PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* test bodyTopo functions */
  PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", Nb));

  for (b = 0; b < Nb; ++b) {
    ego body = bodies[b];
    int id, sh, Nsh, f, Nf, l, Nl, e, Ne, v, Nv;
	
    /* List Topology of Bodies */
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "   BODY %d TOPOLOGY SUMMARY \n", b));

    /* Output Basic Model Topology */
    PetscCall(EGlite_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of SHELLS: %d \n", Nsh));

    PetscCall(EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of FACES: %d \n", Nf));

    PetscCall(EGlite_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of LOOPS: %d \n", Nl));

    PetscCall(EGlite_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of EDGES: %d \n", Ne));

    PetscCall(EGlite_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of NODES: %d \n", Nv));
	
    EGlite_free(shobjs);
    EGlite_free(fobjs);
    EGlite_free(lobjs);
    EGlite_free(eobjs);
    EGlite_free(nobjs);
	
    /* List Topology of Bodies */
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      TOPOLOGY DETAILS \n", b));

    /* Get SHELL info which associated with the current BODY */
    PetscCall(EGlite_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses));
	
    for (sh = 0; sh < Nsh; ++sh) {
      ego shell   = shobjs[sh];
      int shsense = shsenses[sh];
      
      id   = EGlite_indexBodyTopo(body, shell);
      PetscCall(PetscPrintf(PETSC_COMM_SELF, "         SHELL ID: %d :: sense = %d\n", id, shsense));
      
      /* Get FACE infor associated with current SHELL */
      PetscCall(EGlite_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses));
      
      for (f = 0; f < Nf; ++f) {
        ego face   = fobjs[f];

        id   = EGlite_indexBodyTopo(body, face);
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "           FACE ID: %d \n", id));
          
        /* Get LOOP info associated with current FACE */
        PetscCall(EGlite_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses));
        
        for (l = 0; l < Nl; ++l) {
          ego loop   = lobjs[l];
          int lsense = lsenses[l];

          id   = EGlite_indexBodyTopo(body, loop);
          PetscCall(PetscPrintf(PETSC_COMM_SELF, "             LOOP ID: %d :: sense = %d\n", id, lsense));

          /* Get EDGE info associated with the current LOOP */
          PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &esenses));
          for (e = 0; e < Ne; ++e) {
            ego    edge      = eobjs[e];
            ego    topRef, prev, next;
            int    esense    = esenses[e];
            double range[4]  = {0., 0., 0., 0.};
            int    peri;

            id   = EGlite_indexBodyTopo(body, edge); //CHKERRQ(ierr);
            PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "               EDGE ID: %d :: sense = %d\n", id, esense));
                  
            if (mtype == DEGENERATE) { PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 EDGE %d is DEGENERATE \n", id)); }
            PetscCall(EGlite_getRange(edge, range, &peri));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 Peri = %d :: Range = %lf, %lf, %lf, %lf \n", peri, range[0], range[1], range[2], range[3]));

            /* Get NODE info associated with the current EDGE */
            PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
                  
            for (v = 0; v < Nv; ++v) {
              ego    vertex = nobjs[v];
              double limits[4];
              int    dummy;
    
              PetscCall(EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
              id = EGlite_indexBodyTopo(body, vertex);
              PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 NODE ID: %d \n", id));
              PetscCall(PetscPrintf(PETSC_COMM_SELF, "                    (x, y, z) = (%lf, %lf, %lf) \n", limits[0], limits[1], limits[2]));
            }
          }
        }
      }
    }
  }
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n\n"));
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCreateEGADSlite(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, vertexLabel;
  // EGADSLite variables
  ego             geom, *bodies, *mobjs, *fobjs, *lobjs, *eobjs, *nobjs;
  ego             topRef, prev, next;
  int             oclass, mtype, nbodies, *senses, *lSenses, *eSenses;
  int             b;
  // PETSc variables
  DM              dm;
  PetscHMapI      edgeMap = NULL, bodyIndexMap = NULL, bodyVertexMap = NULL, bodyEdgeMap = NULL, bodyFaceMap = NULL, bodyEdgeGlobalMap = NULL;
  PetscInt        dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numEdges = 0, numFaces = 0, numCells = 0, edgeCntr = 0;
  PetscInt        cellCntr = 0, numPoints = 0;
  PetscInt        *cells  = NULL;
  const PetscInt  *cone = NULL;
  PetscReal       *coords = NULL;
  PetscMPIInt      rank;
  //PetscErrorCode   ierr;

  PetscFunctionBeginUser;
  PetscCall(MPI_Comm_rank(comm, &rank));
  if (!rank) {
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex
    //  Get all Nodes in model, record coordinates in a correctly formatted array
    //  Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
    //  We need to uniformly refine the initial geometry to guarantee a valid mesh

    // Caluculate cell and vertex sizes
    PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));

    PetscCall(PetscHMapICreate(&edgeMap));
    PetscCall(PetscHMapICreate(&bodyIndexMap));
    PetscCall(PetscHMapICreate(&bodyVertexMap));
    PetscCall(PetscHMapICreate(&bodyEdgeMap));
    PetscCall(PetscHMapICreate(&bodyEdgeGlobalMap));
    PetscCall(PetscHMapICreate(&bodyFaceMap));

    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      int             Nf, Ne, Nv;
      PetscHashIter   BIiter, BViter, BEiter, BEGiter, BFiter, EMiter;
      PetscBool       BIfound, BVfound, BEfound, BEGfound, BFfound, EMfound;

      PetscCall(PetscHMapIFind(bodyIndexMap, b, &BIiter, &BIfound));
      PetscCall(PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound));
      PetscCall(PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound));
      PetscCall(PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound));
      PetscCall(PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound));

      if (!BIfound)  {PetscCall(PetscHMapISet(bodyIndexMap, b, numFaces + numEdges + numVertices));}
      if (!BVfound)  {PetscCall(PetscHMapISet(bodyVertexMap, b, numVertices));}
      if (!BEfound)  {PetscCall(PetscHMapISet(bodyEdgeMap, b, numEdges));}
      if (!BEGfound) {PetscCall(PetscHMapISet(bodyEdgeGlobalMap, b, edgeCntr));}
      if (!BFfound)  {PetscCall(PetscHMapISet(bodyFaceMap, b, numFaces));}

      PetscCall(EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      PetscCall(EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));
      PetscCall(EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      EGlite_free(fobjs);
      EGlite_free(eobjs);
      EGlite_free(nobjs);

      // Remove DEGENERATE EDGES from Edge count
      PetscCall(EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));
      int Netemp = 0;
      for (int e = 0; e < Ne; ++e) {
        ego     edge = eobjs[e];
        int     eid;

        PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
        eid = EGlite_indexBodyTopo(body, edge);

        PetscCall(PetscHMapIFind(edgeMap, edgeCntr + eid - 1, &EMiter, &EMfound));
        if (mtype == DEGENERATE) {
          if (!EMfound) {PetscCall(PetscHMapISet(edgeMap, edgeCntr + eid - 1, -1));}
        } else {
          ++Netemp;
          if (!EMfound) {PetscCall(PetscHMapISet(edgeMap, edgeCntr + eid - 1, Netemp));}
        }
      }
      edgeCntr    += Ne;
      EGlite_free(eobjs);

      // Determine Number of Cells
      PetscCall(EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      for (int f = 0; f < Nf; ++f) {
        ego     face = fobjs[f];
        int     edgeTemp = 0;

        PetscCall(EGlite_getBodyTopos(body, face, EDGE, &Ne, &eobjs));
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];

          PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
          if (mtype != DEGENERATE) {++edgeTemp;}
        }
        numCells += (2 * edgeTemp);
        EGlite_free(eobjs);
      }
      EGlite_free(fobjs);

      numFaces    += Nf;
      numEdges    += Netemp;
      numVertices += Nv;
    }

    // Set up basic DMPlex parameters
    dim        = 2;     // Assumes 3D Models :: Need to handle 2D modles in the future
    cdim       = 3;     // Assumes 3D Models :: Need to update to handle 2D modles in future
    numCorners = 3;     // Split Faces into triangles
    numPoints  = numVertices + numEdges + numFaces;   // total number of coordinate points

    PetscCall(PetscMalloc2(numPoints*cdim, &coords, numCells*numCorners, &cells));

    // Get Vertex Coordinates and Set up Cells
    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      int             Nf, Ne, Nv;
      PetscInt        bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
      PetscHashIter   BViter, BEiter, BEGiter, BFiter, EMiter;
      PetscBool       BVfound, BEfound, BEGfound, BFfound, EMfound;

      // Vertices on Current Body
      PetscCall(EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));

      PetscCall(PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound));
      if (!BVfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);
      PetscCall(PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart));

      for (int v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    id, dummy;

        PetscCall(EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
        id = EGlite_indexBodyTopo(body, vertex);
	
        coords[(bodyVertexIndexStart + id - 1)*cdim + 0] = limits[0];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 1] = limits[1];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 2] = limits[2];
      }
      EGlite_free(nobjs);

      // Edge Midpoint Vertices on Current Body
      PetscCall(EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));

      PetscCall(PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound));
      if (!BEfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);
      PetscCall(PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart));

      PetscCall(PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound));
      if (!BEGfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);
      PetscCall(PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart));

      for (int e = 0; e < Ne; ++e) {
        ego          edge = eobjs[e];
        double       range[2], avgt[1], cntrPnt[9];
        int          eid, eOffset;
        int          periodic;

        PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
        if (mtype == DEGENERATE) {continue;}

        eid = EGlite_indexBodyTopo(body, edge);

        // get relative offset from globalEdgeID Vector
        PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
        if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d not found in edgeMap", bodyEdgeGlobalIndexStart + eid - 1);
        PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

        PetscCall(EGlite_getRange(edge, range, &periodic));
        avgt[0] = (range[0] + range[1]) /  2.;

        PetscCall(EGlite_evaluate(edge, avgt, cntrPnt));
	  
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 2] = cntrPnt[2];
      }
      EGlite_free(eobjs);

      // Face Midpoint Vertices on Current Body
      PetscCall(EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));

      PetscCall(PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound));
      if (!BFfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);
      PetscCall(PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart));

      for (int f = 0; f < Nf; ++f) {
        ego       face = fobjs[f];
        double    range[4], avgUV[2], cntrPnt[18];
        int       peri, id;

        id = EGlite_indexBodyTopo(body, face);
        PetscCall(EGlite_getRange(face, range, &peri));

        avgUV[0] = (range[0] + range[1]) / 2.;
        avgUV[1] = (range[2] + range[3]) / 2.;
        PetscCall(EGlite_evaluate(face, avgUV, cntrPnt));

        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 2] = cntrPnt[2];
      }
      EGlite_free(fobjs);

      // Define Cells :: Note - This could be incorporated in the Face Midpoint Vertices Loop but was kept separate for clarity
      PetscCall(EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      for (int f = 0; f < Nf; ++f) {
        ego      face = fobjs[f];
        int      fID, midFaceID, midPntID, startID, endID, Nl;

        fID = EGlite_indexBodyTopo(body, face);
        midFaceID = numVertices + numEdges + bodyFaceIndexStart + fID - 1;
        // Must Traverse Loop to ensure we have all necessary information like the sense (+/- 1) of the edges.
        // TODO :: Only handles single loop faces (No holes). The choices for handling multiloop faces are:
        //            1) Use the DMPlexCreateEGADSFromFile() with the -dm_plex_egads_with_tess = 1 option.
        //               This will use a default EGADS tessellation as an initial surface mesh.
        //            2) Create the initial surface mesh via a 2D mesher :: Currently not availble (?future?)
        //               May I suggest the XXXX as a starting point?

        PetscCall(EGlite_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lSenses));

        if (Nl > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face has %d Loops. Can only handle Faces with 1 Loop. Please use --dm_plex_egads_with_tess = 1 Option", Nl);
        for (int l = 0; l < Nl; ++l) {
          ego      loop = lobjs[l];

          PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses));
          for (int e = 0; e < Ne; ++e) {
            ego     edge = eobjs[e];
            int     eid, eOffset;

            PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
            eid = EGlite_indexBodyTopo(body, edge);
            if (mtype == DEGENERATE) { continue; }

            // get relative offset from globalEdgeID Vector
            PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
            if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
            PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

            midPntID = numVertices + bodyEdgeIndexStart + eOffset - 1;

            PetscCall(EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));

            if (eSenses[e] > 0) { startID = EGlite_indexBodyTopo(body, nobjs[0]); endID = EGlite_indexBodyTopo(body, nobjs[1]); }
            else { startID = EGlite_indexBodyTopo(body, nobjs[1]); endID = EGlite_indexBodyTopo(body, nobjs[0]); }

            // Define 2 Cells per Edge with correct orientation
            cells[cellCntr*numCorners + 0] = midFaceID;
            cells[cellCntr*numCorners + 1] = bodyVertexIndexStart + startID - 1;
            cells[cellCntr*numCorners + 2] = midPntID;

            cells[cellCntr*numCorners + 3] = midFaceID;
            cells[cellCntr*numCorners + 4] = midPntID;
            cells[cellCntr*numCorners + 5] = bodyVertexIndexStart + endID - 1;

            cellCntr = cellCntr + 2;
          }
        }
      }
      EGlite_free(fobjs);
    }
  }

  // Generate DMPlex
  PetscCall(DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm));
  PetscCall(PetscFree2(coords, cells));
  PetscCall(PetscInfo(dm, " Total Number of Unique Cells    = %D \n", numCells));
  PetscCall(PetscInfo(dm, " Total Number of Unique Vertices = %D \n", numVertices));

  // Embed EGADS model in DM
  {
    PetscContainer modelObj, contextObj;

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &modelObj));
    PetscCall(PetscContainerSetPointer(modelObj, model));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Context", (PetscObject) contextObj));
    PetscCall(PetscContainerDestroy(&contextObj));
  }
  // Label points
  PetscInt   nStart, nEnd;

  PetscCall(DMCreateLabel(dm, "EGADS Body ID"));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Face ID"));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Edge ID"));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Vertex ID"));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));

  PetscCall(DMPlexGetHeightStratum(dm, 2, &nStart, &nEnd));

  cellCntr = 0;
  for (b = 0; b < nbodies; ++b) {
    ego             body = bodies[b];
    int             Nv, Ne, Nf;
    PetscInt        bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
    PetscHashIter   BViter, BEiter, BEGiter, BFiter, EMiter;
    PetscBool       BVfound, BEfound, BEGfound, BFfound, EMfound;

    PetscCall(PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound));
    if (!BVfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);
    PetscCall(PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart));

    PetscCall(PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound));
    if (!BEfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);
    PetscCall(PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart));

    PetscCall(PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound));
    if (!BFfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);
    PetscCall(PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart));

    PetscCall(PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound));
    if (!BEGfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);
    PetscCall(PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart));

    PetscCall(EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
    for (int f = 0; f < Nf; ++f) {
      ego   face = fobjs[f];
      int   fID, Nl;

      fID  = EGlite_indexBodyTopo(body, face);

      PetscCall(EGlite_getBodyTopos(body, face, LOOP, &Nl, &lobjs));
      for (int l = 0; l < Nl; ++l) {
        ego  loop = lobjs[l];
        int  lid;

        lid  = EGlite_indexBodyTopo(body, loop);
        if (Nl > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);

        PetscCall(EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses));
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];
          int     eid, eOffset;

          // Skip DEGENERATE Edges
          PetscCall(EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
          if (mtype == DEGENERATE) {continue;}
          eid = EGlite_indexBodyTopo(body, edge);

          // get relative offset from globalEdgeID Vector
          PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
          if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
          PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

          PetscCall(EGlite_getBodyTopos(body, edge, NODE, &Nv, &nobjs));
          for (int v = 0; v < Nv; ++v){
            ego vertex = nobjs[v];
            int vID;

            vID = EGlite_indexBodyTopo(body, vertex);
            PetscCall(DMLabelSetValue(bodyLabel, nStart + bodyVertexIndexStart + vID - 1, b));
            PetscCall(DMLabelSetValue(vertexLabel, nStart + bodyVertexIndexStart + vID - 1, vID));
          }
          EGlite_free(nobjs);

          PetscCall(DMLabelSetValue(bodyLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, b));
          PetscCall(DMLabelSetValue(edgeLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, eid));

          // Define Cell faces
          for (int jj = 0; jj < 2; ++jj){
            PetscCall(DMLabelSetValue(bodyLabel, cellCntr, b));
            PetscCall(DMLabelSetValue(faceLabel, cellCntr, fID));
            PetscCall(DMPlexGetCone(dm, cellCntr, &cone));

            PetscCall(DMLabelSetValue(bodyLabel, cone[0], b));
            PetscCall(DMLabelSetValue(faceLabel, cone[0], fID));

            PetscCall(DMLabelSetValue(bodyLabel, cone[1], b));
            PetscCall(DMLabelSetValue(edgeLabel, cone[1], eid));

            PetscCall(DMLabelSetValue(bodyLabel, cone[2], b));
            PetscCall(DMLabelSetValue(faceLabel, cone[2], fID));

            cellCntr = cellCntr + 1;
          }
        }
      }
      EGlite_free(lobjs);

      PetscCall(DMLabelSetValue(bodyLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, b));
      PetscCall(DMLabelSetValue(faceLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, fID));
    }
    EGlite_free(fobjs);
  }

  PetscCall(PetscHMapIDestroy(&edgeMap));
  PetscCall(PetscHMapIDestroy(&bodyIndexMap));
  PetscCall(PetscHMapIDestroy(&bodyVertexMap));
  PetscCall(PetscHMapIDestroy(&bodyEdgeMap));
  PetscCall(PetscHMapIDestroy(&bodyEdgeGlobalMap));
  PetscCall(PetscHMapIDestroy(&bodyFaceMap));

  *newdm = dm;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCreateEGADSlite_Tess_Internal(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel              bodyLabel, faceLabel, edgeLabel, vertexLabel;
  /* EGADSLite variables */
  ego                  geom, *bodies, *fobjs;
  int                  b, oclass, mtype, nbodies, *senses;
  int                  totalNumTris = 0, totalNumPoints = 0;
  double               boundBox[6] = {0., 0., 0., 0., 0., 0.}, tessSize;
  /* PETSc variables */
  DM                   dm;
  PetscHMapI           pointIndexStartMap = NULL, triIndexStartMap = NULL, pTypeLabelMap = NULL, pIndexLabelMap = NULL;
  PetscHMapI           pBodyIndexLabelMap = NULL, triFaceIDLabelMap = NULL, triBodyIDLabelMap = NULL;
  PetscInt             dim = -1, cdim = -1, numCorners = 0, counter = 0;
  PetscInt            *cells  = NULL;
  const PetscInt      *cone = NULL;
  PetscReal           *coords = NULL;
  PetscMPIInt          rank;
  //PetscErrorCode       ierr;

  PetscFunctionBeginUser;
  PetscCall(MPI_Comm_rank(comm, &rank));
  if (!rank) {
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex from EGADSlite created Tessellation of geometry
    // ---------------------------------------------------------------------------------------------------

    // Caluculate cell and vertex sizes
    PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));

    PetscCall(PetscHMapICreate(&pointIndexStartMap));
    PetscCall(PetscHMapICreate(&triIndexStartMap));
    PetscCall(PetscHMapICreate(&pTypeLabelMap));
    PetscCall(PetscHMapICreate(&pIndexLabelMap));
    PetscCall(PetscHMapICreate(&pBodyIndexLabelMap));
    PetscCall(PetscHMapICreate(&triFaceIDLabelMap));
    PetscCall(PetscHMapICreate(&triBodyIDLabelMap));

    /* Create Tessellation of Bodies */
    ego tessArray[nbodies];

    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      double          params[3] = {0.0, 0.0, 0.0};    // Parameters for Tessellation
      int             Nf, bodyNumPoints = 0, bodyNumTris = 0;
      PetscHashIter   PISiter, TISiter;
      PetscBool       PISfound, TISfound;

      /* Store Start Index for each Body's Point and Tris */
      PetscCall(PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound));
      PetscCall(PetscHMapIFind(triIndexStartMap, b, &TISiter, &TISfound));

      if (!PISfound)  {PetscCall(PetscHMapISet(pointIndexStartMap, b, totalNumPoints));}
      if (!TISfound)  {PetscCall(PetscHMapISet(triIndexStartMap, b, totalNumTris));}

      /* Calculate Tessellation parameters based on Bounding Box */
      /* Get Bounding Box Dimensions of the BODY */
      PetscCall(EGlite_getBoundingBox(body, boundBox));
      tessSize = boundBox[3] - boundBox[0];
      if (tessSize < boundBox[4] - boundBox[1]) tessSize = boundBox[4] - boundBox[1];
      if (tessSize < boundBox[5] - boundBox[2]) tessSize = boundBox[5] - boundBox[2];

      // TODO :: May want to give users tessellation parameter options //
      params[0] = 0.0250 * tessSize;
      params[1] = 0.0075 * tessSize;
      params[2] = 15.0;

      PetscCall(EGlite_makeTessBody(body, params, &tessArray[b]));

      PetscCall(EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));

      for (int f = 0; f < Nf; ++f) {
        ego             face = fobjs[f];
        int             len, fID, ntris;
        const int      *ptype, *pindex, *ptris, *ptric;
        const double   *pxyz, *puv;

        // Get Face ID //
        fID = EGlite_indexBodyTopo(body, face);

        // Checkout the Surface Tessellation //
        PetscCall(EGlite_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric));

        // Determine total number of triangle cells in the tessellation //
        bodyNumTris += (int) ntris;

        // Check out the point index and coordinate //
        for (int p = 0; p < len; ++p) {
          int global;

          PetscCall(EGlite_localToGlobal(tessArray[b], fID, p+1, &global));

          // Determine the total number of points in the tessellation //
          bodyNumPoints = PetscMax(bodyNumPoints, global);
        }
      }
      EGlite_free(fobjs);

      totalNumPoints += bodyNumPoints;
      totalNumTris += bodyNumTris;
    }
  //}  - Original End of (!rank)

    dim = 2;
    cdim = 3;
    numCorners = 3;
    //PetscInt counter = 0;

    /* NEED TO DEFINE MATRICES/VECTORS TO STORE GEOM REFERENCE DATA   */
    /* Fill in below and use to define DMLabels after DMPlex creation */
    PetscCall(PetscMalloc2(totalNumPoints*cdim, &coords, totalNumTris*numCorners, &cells));

    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      int             Nf;
      PetscInt        pointIndexStart;
      PetscHashIter   PISiter;
      PetscBool       PISfound;

      PetscCall(PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound));
      if (!PISfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in pointIndexStartMap", b);
      PetscCall(PetscHMapIGet(pointIndexStartMap, b, &pointIndexStart));

      PetscCall(EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));

      for (int f = 0; f < Nf; ++f) {
        /* Get Face Object */
        ego              face = fobjs[f];
        int              len, fID, ntris;
        const int       *ptype, *pindex, *ptris, *ptric;
        const double    *pxyz, *puv;

        /* Get Face ID */
        fID = EGlite_indexBodyTopo(body, face);

        /* Checkout the Surface Tessellation */
        PetscCall(EGlite_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric));

        /* Check out the point index and coordinate */
        for (int p = 0; p < len; ++p) {
          int              global;
          PetscHashIter    PTLiter, PILiter, PBLiter;
          PetscBool        PTLfound, PILfound, PBLfound;

          PetscCall(EGlite_localToGlobal(tessArray[b], fID, p+1, &global));

          /* Set the coordinates array for DAG */
          coords[((global-1+pointIndexStart)*3) + 0] = pxyz[(p*3)+0];
          coords[((global-1+pointIndexStart)*3) + 1] = pxyz[(p*3)+1];
          coords[((global-1+pointIndexStart)*3) + 2] = pxyz[(p*3)+2];

          /* Store Geometry Label Information for DMLabel assignment later */
          PetscCall(PetscHMapIFind(pTypeLabelMap, global-1+pointIndexStart, &PTLiter, &PTLfound));
          PetscCall(PetscHMapIFind(pIndexLabelMap, global-1+pointIndexStart, &PILiter, &PILfound));
          PetscCall(PetscHMapIFind(pBodyIndexLabelMap, global-1+pointIndexStart, &PBLiter, &PBLfound));

          if (!PTLfound)  {PetscCall(PetscHMapISet(pTypeLabelMap, global-1+pointIndexStart, ptype[p]));}
          if (!PILfound)  {PetscCall(PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, pindex[p]));}
          if (!PBLfound)  {PetscCall(PetscHMapISet(pBodyIndexLabelMap, global-1+pointIndexStart, b));}

          if (ptype[p] < 0) { PetscCall(PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, fID));}
        }

        for (int t = 0; t < (int) ntris; ++t){
          int             global, globalA, globalB;
          PetscHashIter   TFLiter, TBLiter;
          PetscBool       TFLfound, TBLfound;

          PetscCall(EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 0], &global));
          cells[(counter*3) +0] = global-1+pointIndexStart;

          PetscCall(EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 1], &globalA));
          cells[(counter*3) +1] = globalA-1+pointIndexStart;

          PetscCall(EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 2], &globalB));
          cells[(counter*3) +2] = globalB-1+pointIndexStart;

          PetscCall(PetscHMapIFind(triFaceIDLabelMap, counter, &TFLiter, &TFLfound));
          PetscCall(PetscHMapIFind(triBodyIDLabelMap, counter, &TBLiter, &TBLfound));

          if (!TFLfound)  {PetscCall(PetscHMapISet(triFaceIDLabelMap, counter, fID));}
          if (!TBLfound)  {PetscCall(PetscHMapISet(triBodyIDLabelMap, counter, b));}

          counter += 1;
        }
      }
      EGlite_free(fobjs);
    }
  }

  //Build DMPlex
  PetscCall(DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, totalNumTris, totalNumPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm));
  PetscCall(PetscFree2(coords, cells));

  // Embed EGADS model in DM
  {
    PetscContainer modelObj, contextObj;

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &modelObj));
    PetscCall(PetscContainerSetPointer(modelObj, model));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADSLite Context", (PetscObject) contextObj));
    PetscCall(PetscContainerDestroy(&contextObj));
  }

  // Label Points
  PetscCall(DMCreateLabel(dm, "EGADS Body ID"));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Face ID"));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Edge ID"));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMCreateLabel(dm, "EGADS Vertex ID"));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));

   /* Get Number of DAG Nodes at each level */
  int   fStart, fEnd, eStart, eEnd, nStart, nEnd;

  PetscCall(DMPlexGetHeightStratum(dm, 0, &fStart, &fEnd));
  PetscCall(DMPlexGetHeightStratum(dm, 1, &eStart, &eEnd));
  PetscCall(DMPlexGetHeightStratum(dm, 2, &nStart, &nEnd));

  /* Set DMLabels for NODES */
  for (int n = nStart; n < nEnd; ++n) {
    int             pTypeVal, pIndexVal, pBodyVal;
    PetscHashIter   PTLiter, PILiter, PBLiter;
    PetscBool       PTLfound, PILfound, PBLfound;

    //Converted to Hash Tables
    PetscCall(PetscHMapIFind(pTypeLabelMap, n - nStart, &PTLiter, &PTLfound));
    if (!PTLfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pTypeLabelMap", n);
    PetscCall(PetscHMapIGet(pTypeLabelMap, n - nStart, &pTypeVal));

    PetscCall(PetscHMapIFind(pIndexLabelMap, n - nStart, &PILiter, &PILfound));
    if (!PILfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pIndexLabelMap", n);
    PetscCall(PetscHMapIGet(pIndexLabelMap, n - nStart, &pIndexVal));

    PetscCall(PetscHMapIFind(pBodyIndexLabelMap, n - nStart, &PBLiter, &PBLfound));
    if (!PBLfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pBodyLabelMap", n);
    PetscCall(PetscHMapIGet(pBodyIndexLabelMap, n - nStart, &pBodyVal));

    PetscCall(DMLabelSetValue(bodyLabel, n, pBodyVal));
    if (pTypeVal == 0) {PetscCall(DMLabelSetValue(vertexLabel, n, pIndexVal));}
    if (pTypeVal >  0) {PetscCall(DMLabelSetValue(edgeLabel, n, pIndexVal));}
    if (pTypeVal <  0) {PetscCall(DMLabelSetValue(faceLabel, n, pIndexVal));}
  }

  /* Set DMLabels for Edges - Based on the DMLabels of the EDGE's NODES */
  for (int e = eStart; e < eEnd; ++e) {
    int    bodyID_0, vertexID_0, vertexID_1, edgeID_0, edgeID_1, faceID_0, faceID_1;

    PetscCall(DMPlexGetCone(dm, e, &cone));
    PetscCall(DMLabelGetValue(bodyLabel, cone[0], &bodyID_0));    // Do I need to check the other end?
    PetscCall(DMLabelGetValue(vertexLabel, cone[0], &vertexID_0));
    PetscCall(DMLabelGetValue(vertexLabel, cone[1], &vertexID_1));
    PetscCall(DMLabelGetValue(edgeLabel, cone[0], &edgeID_0));
    PetscCall(DMLabelGetValue(edgeLabel, cone[1], &edgeID_1));
    PetscCall(DMLabelGetValue(faceLabel, cone[0], &faceID_0));
    PetscCall(DMLabelGetValue(faceLabel, cone[1], &faceID_1));

    PetscCall(DMLabelSetValue(bodyLabel, e, bodyID_0));

    if (edgeID_0 == edgeID_1) { PetscCall(DMLabelSetValue(edgeLabel, e, edgeID_0)); }
    else if (vertexID_0 > 0 && edgeID_1 > 0) { PetscCall(DMLabelSetValue(edgeLabel, e, edgeID_1)); }
    else if (vertexID_1 > 0 && edgeID_0 > 0) { PetscCall(DMLabelSetValue(edgeLabel, e, edgeID_0)); }
    else { /* Do Nothing */ }
  }

  /* Set DMLabels for Cells */
  for (int f = fStart; f < fEnd; ++f){
    int             edgeID_0;
    PetscInt        triBodyVal, triFaceVal;
    PetscHashIter   TFLiter, TBLiter;
    PetscBool       TFLfound, TBLfound;

    // Convert to Hash Table
    PetscCall(PetscHMapIFind(triFaceIDLabelMap, f - fStart, &TFLiter, &TFLfound));
    if (!TFLfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triFaceIDLabelMap", f);
    PetscCall(PetscHMapIGet(triFaceIDLabelMap, f - fStart, &triFaceVal));

    PetscCall(PetscHMapIFind(triBodyIDLabelMap, f - fStart, &TBLiter, &TBLfound));
    if (!TBLfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triBodyIDLabelMap", f);
    PetscCall(PetscHMapIGet(triBodyIDLabelMap, f - fStart, &triBodyVal));

    PetscCall(DMLabelSetValue(bodyLabel, f, triBodyVal));
    PetscCall(DMLabelSetValue(faceLabel, f, triFaceVal));

    /* Finish Labeling previously unlabeled DMPlex Edges - Assumes Triangular Cell (3 Edges Max) */
    PetscCall(DMPlexGetCone(dm, f, &cone));

    for (int jj = 0; jj < 3; ++jj) {
      PetscCall(DMLabelGetValue(edgeLabel, cone[jj], &edgeID_0));

      if (edgeID_0 < 0) {
        PetscCall(DMLabelSetValue(bodyLabel, cone[jj], triBodyVal));
        PetscCall(DMLabelSetValue(faceLabel, cone[jj], triFaceVal));
      }
    }
  }

  *newdm = dm;
  PetscFunctionReturn(0);
}
#endif

/*@
  DMPlexInflateToEGADSliteGeomModel - Snaps the vertex coordinates of a DMPlex object representing the mesh to its geometry if some vertices depart from the model. This usually happens with non-conforming refinement.

  Collective on dm

  Input Parameter:
. dm - The uninflated DM object representing the mesh

  Output Parameter:
. dm - The inflated DM object representing the mesh

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
PetscErrorCode DMPlexInflateToEGADSliteGeomModel(DM dm)
{
#if defined(PETSC_HAVE_EGADS)
  /* EGADS Variables */
  ego            model, geom, body, face, edge;
  ego           *bodies;
  int            Nb, oclass, mtype, *senses;
  double         result[3];
  /* PETSc Variables */
  DM             cdm;
  PetscContainer modelObj;
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  Vec            coordinates;
  PetscScalar   *coords;
  PetscInt       bodyID, faceID, edgeID, vertexID;
  PetscInt       cdim, d, vStart, vEnd, v;
  //PetscErrorCode ierr;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_EGADS)
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSLite Model", (PetscObject *) &modelObj));
  if (!modelObj) PetscFunctionReturn(0);
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinates));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));

  PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
  PetscCall(EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));

  PetscCall(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));
  PetscCall(VecGetArrayWrite(coordinates, &coords));
  for (v = vStart; v < vEnd; ++v) {
    PetscScalar *vcoords;

    PetscCall(DMLabelGetValue(bodyLabel, v, &bodyID));
    PetscCall(DMLabelGetValue(faceLabel, v, &faceID));
    PetscCall(DMLabelGetValue(edgeLabel, v, &edgeID));
    PetscCall(DMLabelGetValue(vertexLabel, v, &vertexID));

    if (bodyID >= Nb) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
    body = bodies[bodyID];

    PetscCall(DMPlexPointLocalRef(cdm, v, coords, (void *) &vcoords));
    if (edgeID > 0) {
      /* Snap to EDGE at nearest location */
      double params[1];
      PetscCall(EGlite_objectBodyTopo(body, EDGE, edgeID, &edge));
      PetscCall(EGlite_invEvaluate(edge, vcoords, params, result)); // Get (x,y,z) of nearest point on EDGE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (faceID > 0) {
      /* Snap to FACE at nearest location */
      double params[2];
      PetscCall(EGlite_objectBodyTopo(body, FACE, faceID, &face));
      PetscCall(EGlite_invEvaluate(face, vcoords, params, result)); // Get (x,y,z) of nearest point on FACE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    }
  }
  PetscCall(VecRestoreArrayWrite(coordinates, &coords));
  /* Clear out global coordinates */
  PetscCall(VecDestroy(&dm->coordinates));
#endif
  PetscFunctionReturn(0);
}


/*@C
  DMPlexCreateEGADSLiteFromFile - Create a DMPlex mesh from an EGADSLite file.

  Collective

  Input Parameters:
+ comm     - The MPI communicator
- filename - The name of the EGADSLite file

  Output Parameter:
. dm       - The DM object representing the mesh

  Level: beginner

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSFromFile()
@*/
PetscErrorCode DMPlexCreateEGADSLiteFromFile(MPI_Comm comm, const char filename[], DM *dm)
{
  PetscMPIInt    rank;
#if defined(PETSC_HAVE_EGADS)
  ego            context= NULL, model = NULL;
#endif
  PetscBool      printModel = PETSC_FALSE, tessModel = PETSC_FALSE, newModel = PETSC_FALSE;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidCharPointer(filename, 2);
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_print_model", &printModel, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_tess_model", &tessModel, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_new_model", &newModel, NULL));
  PetscCall(MPI_Comm_rank(comm, &rank));
#if defined(PETSC_HAVE_EGADS)
  if (!rank) {

    PetscCall(EGlite_open(&context));
    PetscCall(EGlite_loadModel(context, 0, filename, &model));
    if (printModel) {PetscCall(DMPlexEGADSLitePrintModel_Internal(model));}

  }
  if (tessModel)     {PetscCall(DMPlexCreateEGADSlite_Tess_Internal(comm, context, model, dm));}
  else if (newModel) {PetscCall(DMPlexCreateEGADSlite_Internal(comm, context, model, dm));}
  else               {PetscCall(DMPlexCreateEGADSlite(comm, context, model, dm));}
  PetscFunctionReturn(0);
#else
  SETERRQ(comm, PETSC_ERR_SUP, "This method requires EGADSLite support. Reconfigure using --download-egads");
#endif
}
