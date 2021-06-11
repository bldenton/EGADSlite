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
  PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  ts[0] = (range[0] + range[1]) / 2.;
  
  while (tolr >= 1.0e-12) {
    ierr = EGlite_evaluate(obj, ts, eval); CHKERRQ(ierr);
    dx = coords[v*dE+0] - eval[0];
	dy = coords[v*dE+1] - eval[1];
	dz = coords[v*dE+2] - eval[2];
	obj_old = dx*dx + dy*dy + dz*dz;
	
	if (obj_old < 1.0E-14) {tolr = obj_old; break;}
	
	A = (eval[3]*eval[3] + eval[4]*eval[4] + eval[5]*eval[5]) * (1.0 + lambda);
	if (A == 0.0) {ierr = PetscPrintf(PETSC_COMM_SELF, "A = 0.0 \n"); CHKERRQ(ierr); break;}
	b = eval[3]*dx + eval[4]*dy + eval[5]*dz;    

	/* Solve A*delta = b */
	delta = b/A;
	
	/* Find a temp (u,v) and associated objective function */
	tt[0] = ts[0] + delta;
	if (tt[0] < range[0]) {tt[0] = range[0]; delta = tt[0] - ts[0];}
	if (tt[0] > range[1]) {tt[0] = range[1]; delta = tt[0] - ts[0];}
	
	ierr = EGlite_evaluate(obj, tt, data); CHKERRQ(ierr);
	
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
	}
	else {
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
  PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  uvs[0] = (range[0] + range[1]) / 2.;
  uvs[1] = (range[2] + range[3]) / 2.;
  
  while (tolr >= 1.0e-12) {
    ierr = EGlite_evaluate(obj, uvs, eval); CHKERRQ(ierr);
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
	if (denom == 0.0) {ierr = PetscPrintf(PETSC_COMM_SELF, "denom = 0.0 \n"); CHKERRQ(ierr);}
	delta[0] = (b[0]*A[3] - b[1]*A[1]) / denom;
	delta[1] = (A[0]*b[1] - A[2]*b[0]) / denom;
	
	/* Find a temp (u,v) and associated objective function */
	uvt[0] = uvs[0] + delta[0];
	uvt[1] = uvs[1] + delta[1];
	
	if (uvt[0] < range[0]) {uvt[0] = range[0]; delta[0] = uvt[0] - uvs[0];}
	if (uvt[0] > range[1]) {uvt[0] = range[1]; delta[0] = uvt[0] - uvs[0];}
	if (uvt[1] < range[2]) {uvt[1] = range[2]; delta[1] = uvt[1] - uvs[1];}
	if (uvt[1] > range[3]) {uvt[1] = range[3]; delta[1] = uvt[1] - uvs[1];}
	
	ierr = EGlite_evaluate(obj, uvt, data); CHKERRQ(ierr);
	
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
	}
	else {
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
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetCoordinateDim(dm, &dE);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinatesLocal);CHKERRQ(ierr);
  ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];

  if      (edgeID >= 0) {ierr = EGlite_objectBodyTopo(body, EDGE, edgeID, &obj);CHKERRQ(ierr); Np = 1;}
  else if (faceID >= 0) {ierr = EGlite_objectBodyTopo(body, FACE, faceID, &obj);CHKERRQ(ierr); Np = 2;}
  else {
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }

  /* Calculate parameters (t or u,v) for vertices */
  ierr = DMPlexVecGetClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
  Nv  /= dE;
  if (Nv == 1) {
    ierr = DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }
  if (Nv > 16) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Cannot handle %D coordinates associated to point %D", Nv, p);

  /* Correct EGADSlite 2pi bug when calculating nearest point on Periodic Surfaces */
  ierr = EGlite_getRange(obj, range, &peri);CHKERRQ(ierr);

  for (v = 0; v < Nv; ++v) {
	if (edgeID > 0) {ierr = DMPlex_EGADSlite_EDGE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV); CHKERRQ(ierr);}
    else            {ierr = DMPlex_EGADSlite_FACE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV); CHKERRQ(ierr);}
  }
  ierr = DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
  /* Calculate parameters (t or u,v) for new vertex at edge midpoint */
  for (pm = 0; pm < Np; ++pm) {
    params[pm] = 0.;
    for (v = 0; v < Nv; ++v) params[pm] += paramsV[v*3+pm];
    params[pm] /= Nv;
  }
  if ((params[0] + pTolr < range[0]) || (params[0] - pTolr > range[1])) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on t or u", p);}
  if (Np > 1 && ((params[1] + pTolr < range[2]) || (params[1] - pTolr > range[3]))) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on v", p);}
  /* Put coordinates for new vertex in result[] */
  ierr = EGlite_evaluate(obj, params, result);CHKERRQ(ierr);
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
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
  if (!rank) {
    const PetscInt debug = 0;

    /* ---------------------------------------------------------------------------------------------------
    Generate Petsc Plex
      Get all Nodes in model, record coordinates in a correctly formatted array
      Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
      We need to uniformly refine the initial geometry to guarantee a valid mesh
    */

    /* Calculate cell and vertex sizes */
    ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&edgeMap);CHKERRQ(ierr);
    numEdges = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l, Nv, v;

      ierr = EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int Ner  = 0, Ne, e, Nc;

        ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);
        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int Nv, id;
          PetscHashIter iter;
          PetscBool     found;

          ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);
          if (mtype == DEGENERATE) continue;
          id   = EGlite_indexBodyTopo(body, edge);CHKERRQ(ierr);
          ierr = PetscHMapIFind(edgeMap, id-1, &iter, &found);CHKERRQ(ierr);
          if (!found) {ierr = PetscHMapISet(edgeMap, id-1, numEdges++);CHKERRQ(ierr);}
          ++Ner;
        }
        if (Ner == 2)      {Nc = 2;}
        else if (Ner == 3) {Nc = 4;}
        else if (Ner == 4) {Nc = 8; ++numQuads;}
        else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Cannot support loop with %d edges", Ner);
        numCells += Nc;
        newCells += Nc-1;
        maxCorners = PetscMax(Ner*2+1, maxCorners);
      }
      ierr = EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);
      for (v = 0; v < Nv; ++v) {
        ego vertex = nobjs[v];

        id = EGlite_indexBodyTopo(body, vertex);
        /* TODO: Instead of assuming contiguous ids, we could use a hash table */
        numVertices = PetscMax(id, numVertices);
      }
      EGlite_free(lobjs);
      EGlite_free(nobjs);
    }
    ierr = PetscHMapIGetSize(edgeMap, &numEdges);CHKERRQ(ierr);
    newVertices  = numEdges + numQuads;
    numVertices += newVertices;

    dim        = 2; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    cdim       = 3; /* Assume 3D Models :: Need to update to handle 2D Models in the future */
    numCorners = 3; /* Split cells into triangles */
    ierr = PetscMalloc3(numVertices*cdim, &coords, numCells*numCorners, &cells, maxCorners, &cone);CHKERRQ(ierr);

    /* Get vertex coordinates */
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nv, v;

      ierr = EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);
      for (v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    dummy;

        ierr = EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
        id   = EGlite_indexBodyTopo(body, vertex);CHKERRQ(ierr);
        coords[(id-1)*cdim+0] = limits[0];
        coords[(id-1)*cdim+1] = limits[1];
        coords[(id-1)*cdim+2] = limits[2];
      }
      EGlite_free(nobjs);
    }
    ierr = PetscHMapIClear(edgeMap);CHKERRQ(ierr);
    fOff     = numVertices - newVertices + numEdges;
    numEdges = 0;
    numQuads = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int Nl, l;

      ierr = EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e;

        lid  = EGlite_indexBodyTopo(body, loop);CHKERRQ(ierr);
        ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);
        for (e = 0; e < Ne; ++e) {
          ego       edge = objs[e];
          int       eid, Nv;
          PetscHashIter iter;
          PetscBool     found;

          ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);
          if (mtype == DEGENERATE) continue;
          ++Ner;
          eid  = EGlite_indexBodyTopo(body, edge);CHKERRQ(ierr);
          ierr = PetscHMapIFind(edgeMap, eid-1, &iter, &found);CHKERRQ(ierr);
          if (!found) {
            PetscInt v = numVertices - newVertices + numEdges;
            double range[4], params[3] = {0., 0., 0.}, result[18];
            int    periodic[2];

            ierr = PetscHMapISet(edgeMap, eid-1, numEdges++);CHKERRQ(ierr);
            ierr = EGlite_getRange(edge, range, periodic);CHKERRQ(ierr);
            params[0] = 0.5*(range[0] + range[1]);
            ierr = EGlite_evaluate(edge, params, result);CHKERRQ(ierr);
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

          ierr = EGlite_getBodyTopos(body, loop, FACE, &Nf, &fobjs);CHKERRQ(ierr);
          face = fobjs[0];
          fid  = EGlite_indexBodyTopo(body, face);CHKERRQ(ierr);
          if (Nf != 1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Loop %d has %d faces, instead of 1 (%d)", lid-1, Nf, fid);
          ierr = EGlite_getRange(face, range, periodic);CHKERRQ(ierr);
          params[0] = 0.5*(range[0] + range[1]);
          params[1] = 0.5*(range[2] + range[3]);
          ierr = EGlite_evaluate(face, params, result);CHKERRQ(ierr);
          coords[v*cdim+0] = result[0];
          coords[v*cdim+1] = result[1];
          coords[v*cdim+2] = result[2];
        }
      }
    }
    if (numEdges + numQuads != newVertices) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Number of new vertices %D != %D previous count", numEdges + numQuads, newVertices);

    /* Get cell vertices by traversing loops */
    numQuads = 0;
    cOff     = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l;

      ierr = EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e, nc = 0, c, Nt, t;

        lid  = EGlite_indexBodyTopo(body, loop);CHKERRQ(ierr);
        ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int points[3];
          int eid, Nv, v, tmp;

          eid  = EGlite_indexBodyTopo(body, edge);
          ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);
          if (mtype == DEGENERATE) continue;
          else                     ++Ner;
          if (Nv != 2) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Edge %d has %d vertices != 2", eid, Nv);

          for (v = 0; v < Nv; ++v) {
            ego vertex = nobjs[v];

            id = EGlite_indexBodyTopo(body, vertex);
            points[v*2] = id-1;
          }
          {
            PetscInt edgeNum;

            ierr = PetscHMapIGet(edgeMap, eid-1, &edgeNum);CHKERRQ(ierr);
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
            else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Edge %d does not match its predecessor", eid);
          }
        }
        if (nc != 2*Ner)     SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Number of corners %D != %D", nc, 2*Ner);
        if (Ner == 4) {cone[nc++] = numVertices - newVertices + numEdges + numQuads++;}
        if (nc > maxCorners) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Number of corners %D > %D max", nc, maxCorners);
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
          default: SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d edges, which we do not support", lid, Ner);
        }
        if (debug) {
          for (t = 0; t < Nt; ++t) {
            ierr = PetscPrintf(PETSC_COMM_SELF, "  LOOP Corner NODEs Triangle %D (", t);CHKERRQ(ierr);
            for (c = 0; c < numCorners; ++c) {
              if (c > 0) {ierr = PetscPrintf(PETSC_COMM_SELF, ", ");CHKERRQ(ierr);}
              ierr = PetscPrintf(PETSC_COMM_SELF, "%D", cells[(cOff-Nt+t)*numCorners+c]);CHKERRQ(ierr);
            }
            ierr = PetscPrintf(PETSC_COMM_SELF, ")\n");CHKERRQ(ierr);
          }
        }
      }
      EGlite_free(lobjs);
    }
  }
  if (cOff != numCells) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Count of total cells %D != %D previous count", cOff, numCells);
  ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numVertices, numCorners, PETSC_TRUE, cells, cdim, coords, &dm);CHKERRQ(ierr);
  ierr = PetscFree3(coords, cells, cone);CHKERRQ(ierr);
  ierr = PetscInfo2(dm, " Total Number of Unique Cells    = %D (%D)\n", numCells, newCells);CHKERRQ(ierr);
  ierr = PetscInfo2(dm, " Total Number of Unique Vertices = %D (%D)\n", numVertices, newVertices);CHKERRQ(ierr);
  /* Embed EGADS model in DM */
  {
    PetscContainer modelObj, contextObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADSLite Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &contextObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(contextObj, context);CHKERRQ(ierr);
    ierr = PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADSLite Context", (PetscObject) contextObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&contextObj);CHKERRQ(ierr);
  }
  /* Label points */
  ierr = DMCreateLabel(dm, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  cOff = 0;
  for (b = 0; b < nbodies; ++b) {
    ego body = bodies[b];
    int id, Nl, l;

    ierr = EGlite_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
    for (l = 0; l < Nl; ++l) {
      ego  loop = lobjs[l];
      ego *fobjs;
      int  lid, Nf, fid, Ner = 0, Ne, e, Nt = 0, t;

      lid  = EGlite_indexBodyTopo(body, loop);CHKERRQ(ierr);
      ierr = EGlite_getBodyTopos(body, loop, FACE, &Nf, &fobjs);CHKERRQ(ierr);
      if (Nf > 1) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);CHKERRQ(ierr);
      fid  = EGlite_indexBodyTopo(body, fobjs[0]);CHKERRQ(ierr);
      EGlite_free(fobjs);
      ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);
      for (e = 0; e < Ne; ++e) {
        ego             edge = objs[e];
        int             eid, Nv, v;
        PetscInt        points[3], support[2], numEdges, edgeNum;
        const PetscInt *edges;

        eid  = EGlite_indexBodyTopo(body, edge);
        ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);
        if (mtype == DEGENERATE) continue;
        else                     ++Ner;
        for (v = 0; v < Nv; ++v) {
          ego vertex = nobjs[v];

          id   = EGlite_indexBodyTopo(body, vertex);
          ierr = DMLabelSetValue(edgeLabel, numCells + id-1, eid);CHKERRQ(ierr);
          points[v*2] = numCells + id-1;
        }
        ierr = PetscHMapIGet(edgeMap, eid-1, &edgeNum);CHKERRQ(ierr);
        points[1] = numCells + numVertices - newVertices + edgeNum;

        ierr = DMLabelSetValue(edgeLabel, points[1], eid);CHKERRQ(ierr);
        support[0] = points[0];
        support[1] = points[1];
        ierr = DMPlexGetJoin(dm, 2, support, &numEdges, &edges);CHKERRQ(ierr);
        if (numEdges != 1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Vertices (%D, %D) should only bound 1 edge, not %D", support[0], support[1], numEdges);
        ierr = DMLabelSetValue(edgeLabel, edges[0], eid);CHKERRQ(ierr);
        ierr = DMPlexRestoreJoin(dm, 2, support, &numEdges, &edges);CHKERRQ(ierr);
        support[0] = points[1];
        support[1] = points[2];
        ierr = DMPlexGetJoin(dm, 2, support, &numEdges, &edges);CHKERRQ(ierr);
        if (numEdges != 1) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Vertices (%D, %D) should only bound 1 edge, not %D", support[0], support[1], numEdges);
        ierr = DMLabelSetValue(edgeLabel, edges[0], eid);CHKERRQ(ierr);
        ierr = DMPlexRestoreJoin(dm, 2, support, &numEdges, &edges);CHKERRQ(ierr);
      }
      switch (Ner) {
        case 2: Nt = 2;break;
        case 3: Nt = 4;break;
        case 4: Nt = 8;break;
        default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Loop with %d edges is unsupported", Ner);
      }
      for (t = 0; t < Nt; ++t) {
        ierr = DMLabelSetValue(bodyLabel, cOff+t, b);CHKERRQ(ierr);
        ierr = DMLabelSetValue(faceLabel, cOff+t, fid);CHKERRQ(ierr);
      }
      cOff += Nt;
    }
    EGlite_free(lobjs);
  }
  ierr = PetscHMapIDestroy(&edgeMap);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  for (c = cStart; c < cEnd; ++c) {
    PetscInt *closure = NULL;
    PetscInt  clSize, cl, bval, fval;

    ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
    ierr = DMLabelGetValue(bodyLabel, c, &bval);CHKERRQ(ierr);
    ierr = DMLabelGetValue(faceLabel, c, &fval);CHKERRQ(ierr);
    for (cl = 0; cl < clSize*2; cl += 2) {
      ierr = DMLabelSetValue(bodyLabel, closure[cl], bval);CHKERRQ(ierr);
      ierr = DMLabelSetValue(faceLabel, closure[cl], fval);CHKERRQ(ierr);
    }
    ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
  }
  *newdm = dm;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexEGADSLitePrintModel_Internal(ego model)
{
  ego            geom, *bodies, *nobjs, *mobjs, *lobjs, *shobjs, *fobjs, *eobjs;
  int            oclass, mtype, *senses, *shsenses, *fsenses, *lsenses, *esenses;
  int            Nb, b;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* test bodyTopo functions */
  ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", Nb);CHKERRQ(ierr);

  for (b = 0; b < Nb; ++b) {
    ego body = bodies[b];
    int id, sh, Nsh, f, Nf, l, Nl, e, Ne, v, Nv;
	
	/* List Topology of Bodies */
	ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "   BODY %d TOPOLOGY SUMMARY \n", b);CHKERRQ(ierr);

    /* Output Basic Model Topology */
    ierr = EGlite_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "      Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "      Number of FACES: %d \n", Nf);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "      Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "      Number of EDGES: %d \n", Ne);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "      Number of NODES: %d \n", Nv);CHKERRQ(ierr);
	
	EGlite_free(shobjs);
	EGlite_free(fobjs);
	EGlite_free(lobjs);
	EGlite_free(eobjs);
	EGlite_free(nobjs);
	
	/* List Topology of Bodies */
	ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "      TOPOLOGY DETAILS \n", b);CHKERRQ(ierr);

    /* Get SHELL info which associated with the current BODY */
    ierr = EGlite_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses);CHKERRQ(ierr);
	
	for (sh = 0; sh < Nsh; ++sh) {
	  ego shell   = shobjs[sh];
	  int shsense = shsenses[sh];
	  
	  id   = EGlite_indexBodyTopo(body, shell);
	  ierr = PetscPrintf(PETSC_COMM_SELF, "         SHELL ID: %d :: sense = %d\n", id, shsense);CHKERRQ(ierr);
	  
	  /* Get FACE infor associated with current SHELL */
	  ierr = EGlite_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses);CHKERRQ(ierr);
	  
	  for (f = 0; f < Nf; ++f) {
		ego face   = fobjs[f];

        id   = EGlite_indexBodyTopo(body, face);
        ierr = PetscPrintf(PETSC_COMM_SELF, "           FACE ID: %d \n", id);CHKERRQ(ierr);
          
		/* Get LOOP info associated with current FACE */
        ierr = EGlite_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses);CHKERRQ(ierr);
        
        for (l = 0; l < Nl; ++l) {
          ego loop   = lobjs[l];
		  int lsense = lsenses[l];

          id   = EGlite_indexBodyTopo(body, loop);
          ierr = PetscPrintf(PETSC_COMM_SELF, "             LOOP ID: %d :: sense = %d\n", id, lsense);CHKERRQ(ierr);

          /* Get EDGE info associated with the current LOOP */
          ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &esenses);CHKERRQ(ierr);
          for (e = 0; e < Ne; ++e) {
            ego    edge      = eobjs[e];
			ego    topRef, prev, next;
			int    esense    = esenses[e];
            double range[4]  = {0., 0., 0., 0.};
            int    peri;

		    id   = EGlite_indexBodyTopo(body, edge); //CHKERRQ(ierr);
			ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "               EDGE ID: %d :: sense = %d\n", id, esense); CHKERRQ(ierr);
                  
			if (mtype == DEGENERATE) { ierr = PetscPrintf(PETSC_COMM_SELF, "                 EDGE %d is DEGENERATE \n", id);CHKERRQ(ierr); }
            ierr = EGlite_getRange(edge, range, &peri);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "                 Peri = %d :: Range = %lf, %lf, %lf, %lf \n", peri, range[0], range[1], range[2], range[3]);

            /* Get NODE info associated with the current EDGE */
            ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);
                  
            for (v = 0; v < Nv; ++v) {
              ego    vertex = nobjs[v];
              double limits[4];
              int    dummy;
    
              ierr = EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
              id = EGlite_indexBodyTopo(body, vertex);
              ierr = PetscPrintf(PETSC_COMM_SELF, "                 NODE ID: %d \n", id);CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_SELF, "                    (x, y, z) = (%lf, %lf, %lf) \n", limits[0], limits[1], limits[2]);
            }
		  }
		}
	  }
    }
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n\n");
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
  PetscErrorCode   ierr;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
  if (!rank) {
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex
    //  Get all Nodes in model, record coordinates in a correctly formatted array
    //  Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
    //  We need to uniformly refine the initial geometry to guarantee a valid mesh

    // Caluculate cell and vertex sizes
    ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);

    ierr = PetscHMapICreate(&edgeMap);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&bodyIndexMap);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&bodyVertexMap);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&bodyEdgeMap);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&bodyEdgeGlobalMap);CHKERRQ(ierr);
    ierr = PetscHMapICreate(&bodyFaceMap);CHKERRQ(ierr);

    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      int             Nf, Ne, Nv;
      PetscHashIter   BIiter, BViter, BEiter, BEGiter, BFiter, EMiter;
      PetscBool       BIfound, BVfound, BEfound, BEGfound, BFfound, EMfound;

      ierr = PetscHMapIFind(bodyIndexMap, b, &BIiter, &BIfound);CHKERRQ(ierr);
      ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound);CHKERRQ(ierr);
      ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound);CHKERRQ(ierr);
      ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound);CHKERRQ(ierr);
      ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound);CHKERRQ(ierr);

      if (!BIfound)  {ierr = PetscHMapISet(bodyIndexMap, b, numFaces + numEdges + numVertices);CHKERRQ(ierr);}
      if (!BVfound)  {ierr = PetscHMapISet(bodyVertexMap, b, numVertices);CHKERRQ(ierr);}
      if (!BEfound)  {ierr = PetscHMapISet(bodyEdgeMap, b, numEdges);CHKERRQ(ierr);}
      if (!BEGfound) {ierr = PetscHMapISet(bodyEdgeGlobalMap, b, edgeCntr);CHKERRQ(ierr);}
      if (!BFfound)  {ierr = PetscHMapISet(bodyFaceMap, b, numFaces);CHKERRQ(ierr);}

      ierr = EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs);CHKERRQ(ierr);
      ierr = EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs);CHKERRQ(ierr);
      ierr = EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);
      EGlite_free(fobjs);
      EGlite_free(eobjs);
      EGlite_free(nobjs);

      // Remove DEGENERATE EDGES from Edge count
      ierr = EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs);CHKERRQ(ierr);
      int Netemp = 0;
      for (int e = 0; e < Ne; ++e) {
        ego     edge = eobjs[e];
        int     eid;

        ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next);CHKERRQ(ierr);
        eid = EGlite_indexBodyTopo(body, edge);CHKERRQ(ierr);

        ierr = PetscHMapIFind(edgeMap, edgeCntr + eid - 1, &EMiter, &EMfound);CHKERRQ(ierr);
        if (mtype == DEGENERATE) {
          if (!EMfound) {ierr = PetscHMapISet(edgeMap, edgeCntr + eid - 1, -1);CHKERRQ(ierr);}
        }
        else {
          ++Netemp;
          if (!EMfound) {ierr = PetscHMapISet(edgeMap, edgeCntr + eid - 1, Netemp);CHKERRQ(ierr);}
        }
      }
	  edgeCntr    += Ne;
      EGlite_free(eobjs);

      // Determine Number of Cells
      ierr = EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs);CHKERRQ(ierr);
      for (int f = 0; f < Nf; ++f) {
        ego     face = fobjs[f];
        int     edgeTemp = 0;

        ierr = EGlite_getBodyTopos(body, face, EDGE, &Ne, &eobjs);CHKERRQ(ierr);
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];

          ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next);CHKERRQ(ierr);
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

    ierr = PetscMalloc2(numPoints*cdim, &coords, numCells*numCorners, &cells);CHKERRQ(ierr);

    // Get Vertex Coordinates and Set up Cells
    for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
      int             Nf, Ne, Nv;
      PetscInt        bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
      PetscHashIter   BViter, BEiter, BEGiter, BFiter, EMiter;
      PetscBool       BVfound, BEfound, BEGfound, BFfound, EMfound;

      // Vertices on Current Body
      ierr = EGlite_getBodyTopos(body, NULL, NODE, &Nv, &nobjs);CHKERRQ(ierr);

      ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound);CHKERRQ(ierr);
      if (!BVfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);
      ierr = PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart);CHKERRQ(ierr);

      for (int v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    id, dummy;

        ierr = EGlite_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
        id = EGlite_indexBodyTopo(body, vertex);CHKERRQ(ierr);
	
        coords[(bodyVertexIndexStart + id - 1)*cdim + 0] = limits[0];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 1] = limits[1];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 2] = limits[2];
      }
      EGlite_free(nobjs);

      // Edge Midpoint Vertices on Current Body
      ierr = EGlite_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs);CHKERRQ(ierr);

      ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound);CHKERRQ(ierr);
      if (!BEfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);
      ierr = PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart);CHKERRQ(ierr);

      ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound);CHKERRQ(ierr);
      if (!BEGfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);
      ierr = PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart);CHKERRQ(ierr);

      for (int e = 0; e < Ne; ++e) {
        ego          edge = eobjs[e];
        double       range[2], avgt[1], cntrPnt[9];
        int          eid, eOffset;
        int          periodic;

        ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next);CHKERRQ(ierr);
        if (mtype == DEGENERATE) {continue;}

        eid = EGlite_indexBodyTopo(body, edge);CHKERRQ(ierr);

        // get relative offset from globalEdgeID Vector
        ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound);CHKERRQ(ierr);
        if (!EMfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d not found in edgeMap", bodyEdgeGlobalIndexStart + eid - 1);
        ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset);CHKERRQ(ierr);

        ierr = EGlite_getRange(edge, range, &periodic);CHKERRQ(ierr);
        avgt[0] = (range[0] + range[1]) /  2.;

        ierr = EGlite_evaluate(edge, avgt, cntrPnt);CHKERRQ(ierr);
	  
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 2] = cntrPnt[2];
      }
      EGlite_free(eobjs);

      // Face Midpoint Vertices on Current Body
      ierr = EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs);CHKERRQ(ierr);

      ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound);CHKERRQ(ierr);
      if (!BFfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);
      ierr = PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart);CHKERRQ(ierr);

      for (int f = 0; f < Nf; ++f) {
        ego       face = fobjs[f];
        double    range[4], avgUV[2], cntrPnt[18];
        int       peri, id;

        id = EGlite_indexBodyTopo(body, face);
        ierr = EGlite_getRange(face, range, &peri);CHKERRQ(ierr);

        avgUV[0] = (range[0] + range[1]) / 2.;
        avgUV[1] = (range[2] + range[3]) / 2.;
        ierr = EGlite_evaluate(face, avgUV, cntrPnt);CHKERRQ(ierr);

        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 2] = cntrPnt[2];
      }
      EGlite_free(fobjs);

      // Define Cells :: Note - This could be incorporated in the Face Midpoint Vertices Loop but was kept separate for clarity
      ierr = EGlite_getBodyTopos(body, NULL, FACE, &Nf, &fobjs);CHKERRQ(ierr);
      for (int f = 0; f < Nf; ++f) {
        ego      face = fobjs[f];
        int      fID, midFaceID, midPntID, startID, endID, Nl;

        fID = EGlite_indexBodyTopo(body, face);CHKERRQ(ierr);
        midFaceID = numVertices + numEdges + bodyFaceIndexStart + fID - 1;
        // Must Traverse Loop to ensure we have all necessary information like the sense (+/- 1) of the edges.
        // TODO :: Only handles single loop faces (No holes). The choices for handling multiloop faces are:
        //            1) Use the DMPlexCreateEGADSFromFile() with the -dm_plex_egads_with_tess = 1 option.
        //               This will use a default EGADS tessellation as an initial surface mesh.
        //            2) Create the initial surface mesh via a 2D mesher :: Currently not availble (?future?)
        //               May I suggest the XXXX as a starting point?

        ierr = EGlite_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lSenses);CHKERRQ(ierr);

        if (Nl > 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face has %d Loops. Can only handle Faces with 1 Loop. Please use --dm_plex_egads_with_tess = 1 Option", Nl);
        for (int l = 0; l < Nl; ++l) {
          ego      loop = lobjs[l];

          ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses);CHKERRQ(ierr);
          for (int e = 0; e < Ne; ++e) {
            ego     edge = eobjs[e];
            int     eid, eOffset;

            ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next);CHKERRQ(ierr);
            eid = EGlite_indexBodyTopo(body, edge);
            if (mtype == DEGENERATE) { continue; }

            // get relative offset from globalEdgeID Vector
            ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound);CHKERRQ(ierr);
            if (!EMfound) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
            ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset);CHKERRQ(ierr);

            midPntID = numVertices + bodyEdgeIndexStart + eOffset - 1;

            ierr = EGlite_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);

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
  ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm);CHKERRQ(ierr);
  ierr = PetscFree2(coords, cells);CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Cells    = %D \n", numCells);CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Vertices = %D \n", numVertices);CHKERRQ(ierr);

  // Embed EGADS model in DM
  {
    PetscContainer modelObj, contextObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &contextObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(contextObj, context);CHKERRQ(ierr);
    ierr = PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADS Context", (PetscObject) contextObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&contextObj);CHKERRQ(ierr);
  }
  // Label points
  PetscInt   nStart, nEnd;

  ierr = DMCreateLabel(dm, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);

  ierr = DMPlexGetHeightStratum(dm, 2, &nStart, &nEnd);CHKERRQ(ierr);

  cellCntr = 0;
  for (b = 0; b < nbodies; ++b) {
    ego             body = bodies[b];
    int             Nv, Ne, Nf;
    PetscInt        bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
    PetscHashIter   BViter, BEiter, BEGiter, BFiter, EMiter;
    PetscBool       BVfound, BEfound, BEGfound, BFfound, EMfound;

    ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound);CHKERRQ(ierr);
    if (!BVfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);
    ierr = PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart);CHKERRQ(ierr);

    ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound);CHKERRQ(ierr);
    if (!BEfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);
    ierr = PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart);CHKERRQ(ierr);

    ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound);CHKERRQ(ierr);
    if (!BFfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);
    ierr = PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart);CHKERRQ(ierr);

    ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound);CHKERRQ(ierr);
    if (!BEGfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);
    ierr = PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
    for (int f = 0; f < Nf; ++f) {
      ego   face = fobjs[f];
      int   fID, Nl;

      fID  = EGlite_indexBodyTopo(body, face);CHKERRQ(ierr);

      ierr = EGlite_getBodyTopos(body, face, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
      for (int l = 0; l < Nl; ++l) {
        ego  loop = lobjs[l];
        int  lid;

        lid  = EGlite_indexBodyTopo(body, loop);CHKERRQ(ierr);
        if (Nl > 1) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);CHKERRQ(ierr);

        ierr = EGlite_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses);CHKERRQ(ierr);
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];
          int     eid, eOffset;

          // Skip DEGENERATE Edges
          ierr = EGlite_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next);CHKERRQ(ierr);
          if (mtype == DEGENERATE) {continue;}
          eid = EGlite_indexBodyTopo(body, edge);CHKERRQ(ierr);

          // get relative offset from globalEdgeID Vector
          ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound);CHKERRQ(ierr);
          if (!EMfound) SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
          ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset);CHKERRQ(ierr);

          ierr = EGlite_getBodyTopos(body, edge, NODE, &Nv, &nobjs);CHKERRQ(ierr);
          for (int v = 0; v < Nv; ++v){
            ego vertex = nobjs[v];
            int vID;

            vID = EGlite_indexBodyTopo(body, vertex);CHKERRQ(ierr);
            ierr = DMLabelSetValue(bodyLabel, nStart + bodyVertexIndexStart + vID - 1, b);CHKERRQ(ierr);
            ierr = DMLabelSetValue(vertexLabel, nStart + bodyVertexIndexStart + vID - 1, vID);CHKERRQ(ierr);
          }
          EGlite_free(nobjs);

          ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, b);CHKERRQ(ierr);
          ierr = DMLabelSetValue(edgeLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, eid);CHKERRQ(ierr);

          // Define Cell faces
          for (int jj = 0; jj < 2; ++jj){
            ierr = DMLabelSetValue(bodyLabel, cellCntr, b);CHKERRQ(ierr);
            ierr = DMLabelSetValue(faceLabel, cellCntr, fID);CHKERRQ(ierr);
            ierr = DMPlexGetCone(dm, cellCntr, &cone);CHKERRQ(ierr);

            ierr = DMLabelSetValue(bodyLabel, cone[0], b);CHKERRQ(ierr);
            ierr = DMLabelSetValue(faceLabel, cone[0], fID);CHKERRQ(ierr);

            ierr = DMLabelSetValue(bodyLabel, cone[1], b);CHKERRQ(ierr);
            ierr = DMLabelSetValue(edgeLabel, cone[1], eid);CHKERRQ(ierr);

            ierr = DMLabelSetValue(bodyLabel, cone[2], b);CHKERRQ(ierr);
            ierr = DMLabelSetValue(faceLabel, cone[2], fID);CHKERRQ(ierr);

            cellCntr = cellCntr + 1;
          }
        }
      }
      EGlite_free(lobjs);

      ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, b);CHKERRQ(ierr);
      ierr = DMLabelSetValue(faceLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, fID);CHKERRQ(ierr);
    }
    EGlite_free(fobjs);
  }

  ierr = PetscHMapIDestroy(&edgeMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyIndexMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyVertexMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeGlobalMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyFaceMap);CHKERRQ(ierr);

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
  PetscErrorCode       ierr;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
  if (!rank) {
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex from EGADSlite created Tessellation of geometry
    // ---------------------------------------------------------------------------------------------------

  // Caluculate cell and vertex sizes
  ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);

  ierr = PetscHMapICreate(&pointIndexStartMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&triIndexStartMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&pTypeLabelMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&pIndexLabelMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&pBodyIndexLabelMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&triFaceIDLabelMap);CHKERRQ(ierr);
  ierr = PetscHMapICreate(&triBodyIDLabelMap);CHKERRQ(ierr);

  /* Create Tessellation of Bodies */
  ego tessArray[nbodies];

  for (b = 0; b < nbodies; ++b) {
    ego             body = bodies[b];
    double          params[3] = {0.0, 0.0, 0.0};    // Parameters for Tessellation
    int             Nf, bodyNumPoints = 0, bodyNumTris = 0;
    PetscHashIter   PISiter, TISiter;
    PetscBool       PISfound, TISfound;

    /* Store Start Index for each Body's Point and Tris */
    ierr = PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound);CHKERRQ(ierr);
    ierr = PetscHMapIFind(triIndexStartMap, b, &TISiter, &TISfound);CHKERRQ(ierr);

    if (!PISfound)  {ierr = PetscHMapISet(pointIndexStartMap, b, totalNumPoints);CHKERRQ(ierr);}
    if (!TISfound)  {ierr = PetscHMapISet(triIndexStartMap, b, totalNumTris);CHKERRQ(ierr);}

    /* Calculate Tessellation parameters based on Bounding Box */
    /* Get Bounding Box Dimensions of the BODY */
    ierr = EGlite_getBoundingBox(body, boundBox);
    tessSize = boundBox[3] - boundBox[0];
    if (tessSize < boundBox[4] - boundBox[1]) tessSize = boundBox[4] - boundBox[1];
    if (tessSize < boundBox[5] - boundBox[2]) tessSize = boundBox[5] - boundBox[2];

    // TODO :: May want to give users tessellation parameter options //
    params[0] = 0.0250 * tessSize;
    params[1] = 0.0075 * tessSize;
    params[2] = 15.0;

    ierr = EGlite_makeTessBody(body, params, &tessArray[b]);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);

    for (int f = 0; f < Nf; ++f) {
      ego             face = fobjs[f];
      int             len, fID, ntris;
      const int      *ptype, *pindex, *ptris, *ptric;
      const double   *pxyz, *puv;

      // Get Face ID //
      fID = EGlite_indexBodyTopo(body, face);

      // Checkout the Surface Tessellation //
      ierr = EGlite_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric);CHKERRQ(ierr);

      // Determine total number of triangle cells in the tessellation //
      bodyNumTris += (int) ntris;

      // Check out the point index and coordinate //
      for (int p = 0; p < len; ++p) {
        int global;

        ierr = EGlite_localToGlobal(tessArray[b], fID, p+1, &global);

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
  ierr = PetscMalloc2(totalNumPoints*cdim, &coords, totalNumTris*numCorners, &cells);CHKERRQ(ierr);

  for (b = 0; b < nbodies; ++b) {
    ego             body = bodies[b];
    int             Nf;
    PetscInt        pointIndexStart;
    PetscHashIter   PISiter;
    PetscBool       PISfound;

    ierr = PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound);CHKERRQ(ierr);
    if (!PISfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in pointIndexStartMap", b);
    ierr = PetscHMapIGet(pointIndexStartMap, b, &pointIndexStart);CHKERRQ(ierr);

    ierr = EGlite_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);

    for (int f = 0; f < Nf; ++f) {
      /* Get Face Object */
      ego              face = fobjs[f];
      int              len, fID, ntris;
      const int       *ptype, *pindex, *ptris, *ptric;
      const double    *pxyz, *puv;

      /* Get Face ID */
      fID = EGlite_indexBodyTopo(body, face);

      /* Checkout the Surface Tessellation */
      ierr = EGlite_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric);CHKERRQ(ierr);

      /* Check out the point index and coordinate */
      for (int p = 0; p < len; ++p) {
        int              global;
        PetscHashIter    PTLiter, PILiter, PBLiter;
        PetscBool        PTLfound, PILfound, PBLfound;

        ierr = EGlite_localToGlobal(tessArray[b], fID, p+1, &global);

        /* Set the coordinates array for DAG */
        coords[((global-1+pointIndexStart)*3) + 0] = pxyz[(p*3)+0];
        coords[((global-1+pointIndexStart)*3) + 1] = pxyz[(p*3)+1];
        coords[((global-1+pointIndexStart)*3) + 2] = pxyz[(p*3)+2];

        /* Store Geometry Label Information for DMLabel assignment later */
        ierr = PetscHMapIFind(pTypeLabelMap, global-1+pointIndexStart, &PTLiter, &PTLfound);CHKERRQ(ierr);
        ierr = PetscHMapIFind(pIndexLabelMap, global-1+pointIndexStart, &PILiter, &PILfound);CHKERRQ(ierr);
        ierr = PetscHMapIFind(pBodyIndexLabelMap, global-1+pointIndexStart, &PBLiter, &PBLfound);CHKERRQ(ierr);

        if (!PTLfound)  {ierr = PetscHMapISet(pTypeLabelMap, global-1+pointIndexStart, ptype[p]);CHKERRQ(ierr);}
        if (!PILfound)  {ierr = PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, pindex[p]);CHKERRQ(ierr);}
        if (!PBLfound)  {ierr = PetscHMapISet(pBodyIndexLabelMap, global-1+pointIndexStart, b);CHKERRQ(ierr);}

        if (ptype[p] < 0) { ierr = PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, fID);CHKERRQ(ierr);}
      }

      for (int t = 0; t < (int) ntris; ++t){
        int             global, globalA, globalB;
        PetscHashIter   TFLiter, TBLiter;
        PetscBool       TFLfound, TBLfound;

        ierr = EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 0], &global);
        cells[(counter*3) +0] = global-1+pointIndexStart;

        ierr = EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 1], &globalA);
        cells[(counter*3) +1] = globalA-1+pointIndexStart;

        ierr = EGlite_localToGlobal(tessArray[b], fID, ptris[(t*3) + 2], &globalB);
        cells[(counter*3) +2] = globalB-1+pointIndexStart;

        ierr = PetscHMapIFind(triFaceIDLabelMap, counter, &TFLiter, &TFLfound);CHKERRQ(ierr);
        ierr = PetscHMapIFind(triBodyIDLabelMap, counter, &TBLiter, &TBLfound);CHKERRQ(ierr);

        if (!TFLfound)  {ierr = PetscHMapISet(triFaceIDLabelMap, counter, fID);CHKERRQ(ierr);}
        if (!TBLfound)  {ierr = PetscHMapISet(triBodyIDLabelMap, counter, b);CHKERRQ(ierr);}

        counter += 1;
      }
    }
    EGlite_free(fobjs);
  }
}

  //Build DMPlex
  ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, totalNumTris, totalNumPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm);CHKERRQ(ierr);
  ierr = PetscFree2(coords, cells);CHKERRQ(ierr);

  // Embed EGADS model in DM
  {
    PetscContainer modelObj, contextObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &contextObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(contextObj, context);CHKERRQ(ierr);
    ierr = PetscContainerSetUserDestroy(contextObj, DMPlexEGADSliteDestroy_Private);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADS Context", (PetscObject) contextObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&contextObj);CHKERRQ(ierr);
  }

  // Label Points
  ierr = DMCreateLabel(dm, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);

   /* Get Number of DAG Nodes at each level */
  int   fStart, fEnd, eStart, eEnd, nStart, nEnd;

  ierr = DMPlexGetHeightStratum(dm, 0, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 2, &nStart, &nEnd);CHKERRQ(ierr);

  /* Set DMLabels for NODES */
  for (int n = nStart; n < nEnd; ++n) {
    int             pTypeVal, pIndexVal, pBodyVal;
    PetscHashIter   PTLiter, PILiter, PBLiter;
    PetscBool       PTLfound, PILfound, PBLfound;

    //Converted to Hash Tables
    ierr = PetscHMapIFind(pTypeLabelMap, n - nStart, &PTLiter, &PTLfound);CHKERRQ(ierr);
    if (!PTLfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pTypeLabelMap", n);
    ierr = PetscHMapIGet(pTypeLabelMap, n - nStart, &pTypeVal);CHKERRQ(ierr);

    ierr = PetscHMapIFind(pIndexLabelMap, n - nStart, &PILiter, &PILfound);CHKERRQ(ierr);
    if (!PILfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pIndexLabelMap", n);
    ierr = PetscHMapIGet(pIndexLabelMap, n - nStart, &pIndexVal);CHKERRQ(ierr);

    ierr = PetscHMapIFind(pBodyIndexLabelMap, n - nStart, &PBLiter, &PBLfound);CHKERRQ(ierr);
    if (!PBLfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pBodyLabelMap", n);
    ierr = PetscHMapIGet(pBodyIndexLabelMap, n - nStart, &pBodyVal);CHKERRQ(ierr);

    ierr = DMLabelSetValue(bodyLabel, n, pBodyVal);CHKERRQ(ierr);
    if (pTypeVal == 0) {ierr = DMLabelSetValue(vertexLabel, n, pIndexVal);CHKERRQ(ierr);}
    if (pTypeVal >  0) {ierr = DMLabelSetValue(edgeLabel, n, pIndexVal);CHKERRQ(ierr);}
    if (pTypeVal <  0) {ierr = DMLabelSetValue(faceLabel, n, pIndexVal);CHKERRQ(ierr);}
  }

  /* Set DMLabels for Edges - Based on the DMLabels of the EDGE's NODES */
  for (int e = eStart; e < eEnd; ++e) {
    int    bodyID_0, vertexID_0, vertexID_1, edgeID_0, edgeID_1, faceID_0, faceID_1;

    ierr = DMPlexGetCone(dm, e, &cone);CHKERRQ(ierr);
    ierr = DMLabelGetValue(bodyLabel, cone[0], &bodyID_0);CHKERRQ(ierr);    // Do I need to check the other end?
    ierr = DMLabelGetValue(vertexLabel, cone[0], &vertexID_0);CHKERRQ(ierr);
    ierr = DMLabelGetValue(vertexLabel, cone[1], &vertexID_1);CHKERRQ(ierr);
    ierr = DMLabelGetValue(edgeLabel, cone[0], &edgeID_0);CHKERRQ(ierr);
    ierr = DMLabelGetValue(edgeLabel, cone[1], &edgeID_1);CHKERRQ(ierr);
    ierr = DMLabelGetValue(faceLabel, cone[0], &faceID_0);CHKERRQ(ierr);
    ierr = DMLabelGetValue(faceLabel, cone[1], &faceID_1);CHKERRQ(ierr);

    ierr = DMLabelSetValue(bodyLabel, e, bodyID_0);CHKERRQ(ierr);

    if (edgeID_0 == edgeID_1) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_0);CHKERRQ(ierr); }
    else if (vertexID_0 > 0 && edgeID_1 > 0) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_1);CHKERRQ(ierr); }
    else if (vertexID_1 > 0 && edgeID_0 > 0) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_0);CHKERRQ(ierr); }
    else { /* Do Nothing */ }
  }

  /* Set DMLabels for Cells */
  for (int f = fStart; f < fEnd; ++f){
    int             edgeID_0;
    PetscInt        triBodyVal, triFaceVal;
    PetscHashIter   TFLiter, TBLiter;
    PetscBool       TFLfound, TBLfound;

    // Convert to Hash Table
    ierr = PetscHMapIFind(triFaceIDLabelMap, f - fStart, &TFLiter, &TFLfound);CHKERRQ(ierr);
    if (!TFLfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triFaceIDLabelMap", f);
    ierr = PetscHMapIGet(triFaceIDLabelMap, f - fStart, &triFaceVal);CHKERRQ(ierr);

    ierr = PetscHMapIFind(triBodyIDLabelMap, f - fStart, &TBLiter, &TBLfound);CHKERRQ(ierr);
    if (!TBLfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triBodyIDLabelMap", f);
    ierr = PetscHMapIGet(triBodyIDLabelMap, f - fStart, &triBodyVal);CHKERRQ(ierr);

    ierr = DMLabelSetValue(bodyLabel, f, triBodyVal);CHKERRQ(ierr);
    ierr = DMLabelSetValue(faceLabel, f, triFaceVal);CHKERRQ(ierr);

    /* Finish Labeling previously unlabeled DMPlex Edges - Assumes Triangular Cell (3 Edges Max) */
    ierr = DMPlexGetCone(dm, f, &cone);CHKERRQ(ierr);

    for (int jj = 0; jj < 3; ++jj) {
      ierr = DMLabelGetValue(edgeLabel, cone[jj], &edgeID_0);CHKERRQ(ierr);

      if (edgeID_0 < 0) {
        ierr = DMLabelSetValue(bodyLabel, cone[jj], triBodyVal);CHKERRQ(ierr);
        ierr = DMLabelSetValue(faceLabel, cone[jj], triFaceVal);CHKERRQ(ierr);
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
  PetscErrorCode ierr;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_EGADS)
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  if (!modelObj) PetscFunctionReturn(0);
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);

  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EGlite_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);

  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = VecGetArrayWrite(coordinates, &coords);CHKERRQ(ierr);
  for (v = vStart; v < vEnd; ++v) {
    PetscScalar *vcoords;

    ierr = DMLabelGetValue(bodyLabel, v, &bodyID);CHKERRQ(ierr);
    ierr = DMLabelGetValue(faceLabel, v, &faceID);CHKERRQ(ierr);
    ierr = DMLabelGetValue(edgeLabel, v, &edgeID);CHKERRQ(ierr);
    ierr = DMLabelGetValue(vertexLabel, v, &vertexID);CHKERRQ(ierr);

    if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
    body = bodies[bodyID];

    ierr = DMPlexPointLocalRef(cdm, v, coords, (void *) &vcoords);CHKERRQ(ierr);
    if (edgeID > 0) {
      /* Snap to EDGE at nearest location */
      double params[1];
      ierr = EGlite_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
      ierr = EGlite_invEvaluate(edge, vcoords, params, result);CHKERRQ(ierr); // Get (x,y,z) of nearest point on EDGE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (faceID > 0) {
      /* Snap to FACE at nearest location */
      double params[2];
      ierr = EGlite_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);
      ierr = EGlite_invEvaluate(face, vcoords, params, result);CHKERRQ(ierr); // Get (x,y,z) of nearest point on FACE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    }
  }
  ierr = VecRestoreArrayWrite(coordinates, &coords);CHKERRQ(ierr);
  /* Clear out global coordinates */
  ierr = VecDestroy(&dm->coordinates);CHKERRQ(ierr);
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
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidCharPointer(filename, 2);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_print_model", &printModel, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_tess_model", &tessModel, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egadslite_new_model", &newModel, NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);
#if defined(PETSC_HAVE_EGADS)
  if (!rank) {

    ierr = EGlite_open(&context);CHKERRQ(ierr);
    ierr = EGlite_loadModel(context, 0, filename, &model);CHKERRQ(ierr);
    if (printModel) {ierr = DMPlexEGADSLitePrintModel_Internal(model);CHKERRQ(ierr);}

  }
  if (tessModel)     {ierr = DMPlexCreateEGADSlite_Tess_Internal(comm, context, model, dm);CHKERRQ(ierr);}
  else if (newModel) {ierr = DMPlexCreateEGADSlite_Internal(comm, context, model, dm);CHKERRQ(ierr);}
  else               {ierr = DMPlexCreateEGADSlite(comm, context, model, dm);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
#else
  SETERRQ(comm, PETSC_ERR_SUP, "This method requires EGADSLite support. Reconfigure using --download-egads");
#endif
}
