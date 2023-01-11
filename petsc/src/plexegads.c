#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include <petsc/private/hashmapi.h>

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

/* We need to understand how to natively parse STEP files. There seems to be only one open source implementation of
   the STEP parser contained in the OpenCASCADE package. It is enough to make a strong man weep:

     https://github.com/tpaviot/oce/tree/master/src/STEPControl

   The STEP, and inner EXPRESS, formats are ISO standards, so they are documented

     https://stackoverflow.com/questions/26774037/documentation-or-specification-for-step-and-stp-files
     http://stepmod.sourceforge.net/express_model_spec/

   but again it seems that there has been a deliberate effort at obfuscation, probably to raise the bar for entrants.
*/

#ifdef PETSC_HAVE_EGADS
PETSC_INTERN PetscErrorCode DMPlexSnapToGeomModel_EGADS_Internal(DM, PetscInt, ego, PetscInt, PetscInt, PetscInt, const PetscScalar[], PetscScalar[]);
PETSC_INTERN PetscErrorCode DMPlexSnapToGeomModel_EGADSLite_Internal(DM, PetscInt, ego, PetscInt, PetscInt, PetscInt, const PetscScalar[], PetscScalar[]);
PETSC_INTERN PetscErrorCode DMPlexInflateToEGADSliteGeomModel(DM);
PETSC_INTERN PetscErrorCode DMPlex_EGADS_EDGE_XYZtoUV_Internal(const PetscScalar[], ego, const PetscScalar[], const PetscInt, const PetscInt, PetscScalar[]);
PETSC_INTERN PetscErrorCode DMPlex_EGADS_FACE_XYZtoUV_Internal(const PetscScalar[], ego, const PetscScalar[], const PetscInt, const PetscInt, PetscScalar[]);


PetscErrorCode DMPlex_EGADS_GeomDecode_Internal(const PetscInt geomClass, const PetscInt geomType, char **retClass, char **retType)
{
  PetscFunctionBeginHot

  /* EGADS Object Type */
  if (geomClass == CONTXT) {*retClass = (char*)"CONTEXT";}
  if (geomClass == TRANSFORM) {*retClass = (char*)"TRANSFORM";}
  if (geomClass == TESSELLATION) {*retClass = (char*)"TESSELLATION";}
  if (geomClass == NIL) {*retClass = (char*)"NIL";}
  if (geomClass == EMPTY) {*retClass = (char*)"EMPTY";}
  if (geomClass == REFERENCE) {*retClass = (char*)"REFERENCE";}
  if (geomClass == PCURVE) {*retClass = (char*)"PCURVE";}
  if (geomClass == CURVE) {*retClass = (char*)"CURVE";}
  if (geomClass == SURFACE) {*retClass = (char*)"SURFACE";}
  if (geomClass == NODE) {*retClass = (char*)"NODE";}
  if (geomClass == EDGE) {*retClass = (char*)"EDGE";}
  if (geomClass == LOOP) {*retClass = (char*)"LOOP";}
  if (geomClass == FACE) {*retClass = (char*)"FACE";}
  if (geomClass == SHELL) {*retClass = (char*)"SHELL";}
  if (geomClass == BODY) {*retClass = (char*)"BODY";}
  if (geomClass == MODEL) {*retClass = (char*)"MODEL";}
	
  /* PCURVES & CURVES */
  if (geomClass == PCURVE || geomClass == CURVE) {
    if (geomType == LINE) {*retType = (char*)"LINE";}
	  if (geomType == CIRCLE) {*retType = (char*)"CIRCLE";}
	  if (geomType == ELLIPSE) {*retType = (char*)"ELLIPSE";}
	  if (geomType == PARABOLA) {*retType = (char*)"PARABOLA";}
	  if (geomType == HYPERBOLA) {*retType = (char*)"HYPERBOLA";}
	  if (geomType == TRIMMED) {*retType = (char*)"TRIMMED";}
	  if (geomType == BEZIER) {*retType = (char*)"BEZIER";}
	  if (geomType == BSPLINE) {*retType = (char*)"BSPLINE";}
	  if (geomType == OFFSET) {*retType = (char*)"OFFSET";}
  }

  /* SURFACE */
  if (geomClass == SURFACE) {
    if (geomType == PLANE) {*retType = (char*)"PLANE";}
	  if (geomType == SPHERICAL) {*retType = (char*)"SPHERICAL";}
	  if (geomType == CYLINDRICAL) {*retType = (char*)"CYLINDRICAL";}
	  if (geomType == REVOLUTION) {*retType = (char*)"REVOLUTION";}
	  if (geomType == TOROIDAL) {*retType = (char*)"TOROIDAL";}
	  if (geomType == CONICAL) {*retType = (char*)"CONICAL";}
	  if (geomType == EXTRUSION) {*retType = (char*)"EXTRUSION";}
	  if (geomType == BEZIER) {*retType = (char*)"BEZIER";}
	  if (geomType == BSPLINE) {*retType = (char*)"BSPLINE";}
  }
  
  /* TOPOLOGY */
  if (geomClass == NODE || geomClass == EDGE || geomClass == LOOP ||
      geomClass == FACE || geomClass == SHELL || geomClass == BODY ||
	    geomClass == MODEL) {
    if (geomType == SREVERSE) {*retType = (char*)"SREVERSE";}
	  if (geomType == NOMTYPE) {*retType = (char*)"NOMTYPE";}
	  if (geomType == SFORWARD && geomClass == FACE) {*retType = (char*)"SFORWARD";}
	  if (geomType == ONENODE && geomClass == EDGE) {*retType = (char*)"ONENODE";}
	  if (geomType == TWONODE) {*retType = (char*)"TWONODE";}
	  if (geomType == OPEN) {*retType = (char*)"OPEN";}
	  if (geomType == CLOSED) {*retType = (char*)"CLOSED";}
	  if (geomType == DEGENERATE) {*retType = (char*)"DEGENERATE";}
	  if (geomType == WIREBODY) {*retType = (char*)"WIREBODY";}
	  if (geomType == FACEBODY) {*retType = (char*)"FACEBODY";}
	  if (geomType == SHEETBODY) {*retType = (char*)"SHEETBODY";}
	  if (geomType == SOLIDBODY) {*retType = (char*)"SOLIDBODY";}
  }
	
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlex_EGADS_EDGE_XYZtoUV_Internal(const PetscScalar coords[], ego obj, const PetscScalar range[], const PetscInt v, const PetscInt dE, PetscScalar paramsV[])
{
  PetscInt       loopCntr = 0;
  PetscScalar    dx, dy, dz, lambda, tolr, obj_old, obj_tmp, target;
  PetscScalar    delta, A, b;
  PetscScalar    ts[2], tt[2], eval[18], data[18];
  //PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  target = 1.0E-20;
  ts[0] = (range[0] + range[1]) / 2.;
  
  while (tolr >= target) {
    PetscCall(EG_evaluate(obj, ts, eval));
    dx = coords[v*dE+0] - eval[0];
	  dy = coords[v*dE+1] - eval[1];
	  dz = coords[v*dE+2] - eval[2];
	  obj_old = dx*dx + dy*dy + dz*dz;
	
	  if (obj_old < target) {tolr = obj_old; break;}
	
	  A = (eval[3]*eval[3] + eval[4]*eval[4] + eval[5]*eval[5]) * (1.0 + lambda);
	  if (A == 0.0) {PetscCall(PetscPrintf(PETSC_COMM_SELF, "A = 0.0 \n")); break;}
	  b = eval[3]*dx + eval[4]*dy + eval[5]*dz;    

	  /* Solve A*delta = b */
	  delta = b/A;
	
	  /* Find a temp (u,v) and associated objective function */
	  tt[0] = ts[0] + delta;
	  if (tt[0] < range[0]) {tt[0] = range[0]; delta = tt[0] - ts[0];}
	  if (tt[0] > range[1]) {tt[0] = range[1]; delta = tt[0] - ts[0];}
	
	  PetscCall(EG_evaluate(obj, tt, data));
	
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
	    if (obj_old < target) {tolr = obj_old; break;}
	  }
	  else {
	    /* Otherwise reject it and double lambda (making it more gradient-descent like) */
	    lambda *= 2.0;
	  }
	
	  if ((tt[0] == range[0]) || (tt[0] == range[1])) break;
	  if (fabs(delta) < target) {tolr = obj_old; break;}
	
	  tolr = obj_old;
	
	  loopCntr += 1;
	  if (loopCntr > 100) break;
  }
  paramsV[v*3+0] = ts[0];
  paramsV[v*3+1] = 0.;
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlex_EGADS_FACE_XYZtoUV_Internal(const PetscScalar coords[], ego obj, const PetscScalar range[], const PetscInt v, const PetscInt dE, PetscScalar paramsV[])
{
  PetscInt       loopCntr = 0;
  PetscScalar    dx, dy, dz, lambda, tolr, denom, obj_old, obj_tmp, target;
  PetscScalar    uvs[2], uvt[2], delta[2], A[4], b[2], eval[18], data[18];
  //PetscErrorCode ierr;
  
  PetscFunctionBeginHot
  /* Initialize Levenberg-Marquardt parameters */
  lambda = 1.0;
  tolr = 1.0;
  target = 1.0E-20;
  uvs[0] = (range[0] + range[1]) / 2.;
  uvs[1] = (range[2] + range[3]) / 2.;
  
  while (tolr >= target) {
    PetscCall(EG_evaluate(obj, uvs, eval));
    dx = coords[v*dE+0] - eval[0];
	  dy = coords[v*dE+1] - eval[1];
	  dz = coords[v*dE+2] - eval[2];
	  obj_old = dx*dx + dy*dy + dz*dz;
	
	  if (obj_old < target) {tolr = obj_old; break;}
	
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
	
	  PetscCall(EG_evaluate(obj, uvt, data));
	
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
	    if (obj_old < target) {tolr = obj_old; break;}
	  }
	  else {
	    /* Otherwise reject it and double lambda (making it more gradient-descent like) */
	    lambda *= 2.0;
	  }
	
	  if (sqrt(delta[0]*delta[0] + delta[1]*delta[1]) < target) {tolr = obj_old; break;}
	
	  tolr = obj_old;
	
	  loopCntr += 1;
	  if (loopCntr > 100) break;
  }
  paramsV[v*3+0] = uvs[0];
  paramsV[v*3+1] = uvs[1];
  
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlexSnapToGeomModel_EGADS_Internal(DM dm, PetscInt p, ego model, PetscInt bodyID, PetscInt faceID, PetscInt edgeID, const PetscScalar mcoords[], PetscScalar gcoords[])
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
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  if (bodyID >= Nb) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  //PetscCheckFalse(bodyID >= Nb, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);    // Call SETERRQ() if true
  body = bodies[bodyID];

  if      (edgeID >= 0) {PetscCall(EG_objectBodyTopo(body, EDGE, edgeID, &obj)); Np = 1;}
  else if (faceID >= 0) {PetscCall(EG_objectBodyTopo(body, FACE, faceID, &obj)); Np = 2;}
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
  //PetscCheckFalse(Nv > 16, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Cannot handle %D coordinates associated to point %D", Nv, p);  // Call SETERRQ() if true

  /* Correct EGADSlite 2pi bug when calculating nearest point on Periodic Surfaces */
  PetscCall(EG_getRange(obj, range, &peri));

  for (v = 0; v < Nv; ++v) {
	  if (edgeID > 0) {PetscCall(DMPlex_EGADS_EDGE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV));}
    else            {PetscCall(DMPlex_EGADS_FACE_XYZtoUV_Internal(coords, obj, range, v, dE, paramsV));}
  }
  PetscCall(DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords));
  /* Calculate parameters (t or u,v) for new vertex at edge midpoint */
  for (pm = 0; pm < Np; ++pm) {
    params[pm] = 0.;
    for (v = 0; v < Nv; ++v) params[pm] += paramsV[v*3+pm];
    params[pm] /= Nv;
  }
  if ((params[0] + pTolr < range[0]) || (params[0] - pTolr > range[1])) {SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on t or u", p);}
  //PetscCheckFalse((params[0] + pTolr < range[0]) || (params[0] - pTolr > range[1]), PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on t or u", p);  // Call SETERRQ() if True
  if (Np > 1 && ((params[1] + pTolr < range[2]) || (params[1] - pTolr > range[3]))) {SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on v", p);}
  //PetscCheckFalse(Np > 1 && ((params[1] + pTolr < range[2]) || (params[1] - pTolr > range[3])), PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Point %D had bad interpolation on v", p);    // Call SETERRQ() if Ture
  /* Put coordinates for new vertex in result[] */
  PetscCall(EG_evaluate(obj, params, result));
  for (d = 0; d < dE; ++d) gcoords[d] = result[d];
  PetscFunctionReturn(0);
}
#endif

/*@
  DMPlexSnapToGeomModel - Given a coordinate point 'mcoords' on the mesh point 'p', return the closest coordinate point 'gcoords' on the geometry model associated with that point.

  Not collective

  Input Parameters:
+ dm      - The DMPlex object
. p       - The mesh point
- mcoords - A coordinate point lying on the mesh point

  Output Parameter:
. gcoords - The closest coordinate point on the geometry model associated with 'p' to the given point

  Note: Returns the original coordinates if no geometry model is found. Right now the only supported geometry model is EGADS.

  Level: intermediate

.seealso: DMRefine(), DMPlexCreate(), DMPlexSetRefinementUniform()
@*/
// Orig :: PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar mcoords[], PetscScalar gcoords[])
PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, PetscInt dE, const PetscScalar mcoords[], PetscScalar gcoords[])
{
  // Orig :: PetscInt       dE, d;
  PetscInt       d;

  PetscFunctionBeginHot;
  // Orig :: PetscCall(DMGetCoordinateDim(dm, &dE));
#ifdef PETSC_HAVE_EGADS
  {
    DM_Plex       *plex = (DM_Plex *) dm->data;
    DMLabel        bodyLabel, faceLabel, edgeLabel;
    PetscInt       bodyID, faceID, edgeID;
    PetscContainer modelObj;
    ego            model;
    PetscBool      islite = PETSC_FALSE;
	  PetscBool      dontSnapToGeom = PETSC_FALSE;
	  //PetscErrorCode ierr;

    /* Added New Option to Skip SnapToGeometry Updates */
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_without_snap_to_geom", &dontSnapToGeom, NULL));
    if (dontSnapToGeom == PETSC_TRUE) {
      for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
      PetscFunctionReturn(0);
    }

    PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
    PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
    PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
    if (!bodyLabel || !faceLabel || !edgeLabel || plex->ignoreModel) {
      for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
      PetscFunctionReturn(0);
    }
    PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
    if (!modelObj) {
      PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSlite Model", (PetscObject *) &modelObj));
      islite = PETSC_TRUE;
    }
    if (!modelObj) {
      for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
      PetscFunctionReturn(0);
    }
    PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
    PetscCall(DMLabelGetValue(bodyLabel, p, &bodyID));
    PetscCall(DMLabelGetValue(faceLabel, p, &faceID));
    PetscCall(DMLabelGetValue(edgeLabel, p, &edgeID));
    /* Allows for "Connective" Plex Edges present in models with multiple non-touching Entities */
    if (bodyID < 0) {
      for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
      PetscFunctionReturn(0);
    }
    if (islite) {PetscCall(DMPlexSnapToGeomModel_EGADSLite_Internal(dm, p, model, bodyID, faceID, edgeID, mcoords, gcoords));}
    else        {PetscCall(DMPlexSnapToGeomModel_EGADS_Internal(dm, p, model, bodyID, faceID, edgeID, mcoords, gcoords));}
  }
#else
  for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
#endif
  PetscFunctionReturn(0);
}


#if defined(PETSC_HAVE_EGADS)
static PetscErrorCode DMPlexEGADSPrintModel_Internal(ego model)
{
  ego            geom, *bodies, *nobjs, *mobjs, *lobjs, *shobjs, *fobjs, *eobjs;
  int            oclass, mtype, *senses, *shsenses, *fsenses, *lsenses, *esenses;
  int            Nb, b;
  //PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* test bodyTopo functions */
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", Nb));

  for (b = 0; b < Nb; ++b) {
    ego body = bodies[b];
    int id, sh, Nsh, f, Nf, l, Nl, e, Ne, v, Nv;
	
	  /* List Topology of Bodies */
	  PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));
	  PetscCall(PetscPrintf(PETSC_COMM_SELF, "   BODY %d TOPOLOGY SUMMARY \n", b));

    /* Output Basic Model Topology */
    PetscCall(EG_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of SHELLS: %d \n", Nsh));

    PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of FACES: %d \n", Nf));

    PetscCall(EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of LOOPS: %d \n", Nl));

    PetscCall(EG_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of EDGES: %d \n", Ne));

    PetscCall(EG_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "      Number of NODES: %d \n", Nv));
	
	  EG_free(shobjs);
	  EG_free(fobjs);
	  EG_free(lobjs);
	  EG_free(eobjs);
	  EG_free(nobjs);
	
	  /* List Topology of Bodies */
	  PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));
	  PetscCall(PetscPrintf(PETSC_COMM_SELF, "      TOPOLOGY DETAILS \n", b));

    /* Get SHELL info which associated with the current BODY */
    PetscCall(EG_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses));
	
	  for (sh = 0; sh < Nsh; ++sh) {
	    ego shell   = shobjs[sh];
	    int shsense = shsenses[sh];
	  
	    id   = EG_indexBodyTopo(body, shell);
	    PetscCall(PetscPrintf(PETSC_COMM_SELF, "         SHELL ID: %d :: sense = %d\n", id, shsense));
	  
	    /* Get FACE infor associated with current SHELL */
	    PetscCall(EG_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses));
	  
	    for (f = 0; f < Nf; ++f) {
		    ego face   = fobjs[f];
		    ego gRef, gPrev, gNext;
		    int goclass, gmtype, *gpinfo;
		    double *gprv;
		    char *gClass=(char*)"", *gType=(char*)"";
		    double fdata[4];
		    ego fRef, fPrev, fNext;
		    int foclass, fmtype;

        id   = EG_indexBodyTopo(body, face);
          
		    /* Get LOOP info associated with current FACE */
        PetscCall(EG_getTopology(face, &geom, &oclass, &mtype, fdata, &Nl, &lobjs, &lsenses));
        PetscCall(EG_getInfo(face, &foclass, &fmtype, &fRef, &fPrev, &fNext));
		    PetscCall(EG_getGeometry(geom, &goclass, &gmtype, &gRef, &gpinfo, &gprv));
		    PetscCall(EG_getInfo(geom, &goclass, &gmtype, &gRef, &gPrev, &gNext));
		    PetscCall(DMPlex_EGADS_GeomDecode_Internal(goclass, gmtype, &gClass, &gType));
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "           FACE ID: %d :: sense = %d \n", id, fmtype));
		    PetscCall(PetscPrintf(PETSC_COMM_SELF, "             GEOMETRY CLASS: %s \n", gClass));
		    PetscCall(PetscPrintf(PETSC_COMM_SELF, "             GEOMETRY TYPE:  %s \n\n", gType));
		    PetscCall(PetscPrintf(PETSC_COMM_SELF, "             RANGE (umin, umax) = (%f, %f) \n", fdata[0], fdata[1]));
		    PetscCall(PetscPrintf(PETSC_COMM_SELF, "                   (vmin, vmax) = (%f, %f) \n\n", fdata[2], fdata[3]));
		
        for (l = 0; l < Nl; ++l) {
          ego loop   = lobjs[l];
          int lsense = lsenses[l];

          id   = EG_indexBodyTopo(body, loop);
          PetscCall(PetscPrintf(PETSC_COMM_SELF, "             LOOP ID: %d :: sense = %d\n", id, lsense));

          /* Get EDGE info associated with the current LOOP */
          PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &esenses));
          for (e = 0; e < Ne; ++e) {
            ego    edge      = eobjs[e];
            ego    topRef, prev, next;
            int    esense    = esenses[e];
            double range[4]  = {0., 0., 0., 0.};
            int    peri;
            ego gEref, gEprev, gEnext;
            int gEoclass, gEmtype;
            char *gEclass=(char*)"", *gEtype=(char*)"";
			
            PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
            PetscCall(EG_getInfo(geom, &gEoclass, &gEmtype, &gEref, &gEprev, &gEnext));
            PetscCall(DMPlex_EGADS_GeomDecode_Internal(gEoclass, gEmtype, &gEclass, &gEtype));
			
            id   = EG_indexBodyTopo(body, edge);
            PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "               EDGE ID: %d :: sense = %d\n", id, esense));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 GEOMETRY CLASS: %s \n", gEclass));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 GEOMETRY TYPE:  %s \n", gEtype));
                  
            if (mtype == DEGENERATE) { PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 EDGE %d is DEGENERATE \n", id)); }
            PetscCall(EG_getRange(edge, range, &peri));
            PetscCall(PetscPrintf(PETSC_COMM_SELF, "                 Peri = %d :: Range = %lf, %lf, %lf, %lf \n", peri, range[0], range[1], range[2], range[3]));

            /* Get NODE info associated with the current EDGE */
            PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
                  
            for (v = 0; v < Nv; ++v) {
              ego    vertex = nobjs[v];
              double limits[4];
              int    dummy;
    
              PetscCall(EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
              id = EG_indexBodyTopo(body, vertex);
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

static PetscErrorCode DMPlexEGADSDestroy_Private(void *context)
{
  if (context) EG_close((ego) context);
  return 0;
}

static PetscErrorCode DMPlexCreateEGADS_Internal(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  PetscInt       cStart, cEnd, c;
  /* EGADS variables */
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
    PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));
    PetscCall(PetscHMapICreate(&edgeMap));
    numEdges = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nl, l, Nv, v;

      PetscCall(EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int Ner  = 0, Ne, e, Nc;

        PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int Nv, id;
          PetscHashIter iter;
          PetscBool     found;

          PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          id   = EG_indexBodyTopo(body, edge);
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
      PetscCall(EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      for (v = 0; v < Nv; ++v) {
        ego vertex = nobjs[v];

        id = EG_indexBodyTopo(body, vertex);
        /* TODO: Instead of assuming contiguous ids, we could use a hash table */
        numVertices = PetscMax(id, numVertices);
      }
      EG_free(lobjs);
      EG_free(nobjs);
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

      PetscCall(EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      for (v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    dummy;

        PetscCall(EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
        id   = EG_indexBodyTopo(body, vertex);
        coords[(id-1)*cdim+0] = limits[0];
        coords[(id-1)*cdim+1] = limits[1];
        coords[(id-1)*cdim+2] = limits[2];
      }
      EG_free(nobjs);
    }
    PetscCall(PetscHMapIClear(edgeMap));
    fOff     = numVertices - newVertices + numEdges;
    numEdges = 0;
    numQuads = 0;
    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int Nl, l;

      PetscCall(EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e;

        lid  = EG_indexBodyTopo(body, loop);
        PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
        for (e = 0; e < Ne; ++e) {
          ego       edge = objs[e];
          int       eid, Nv;
          PetscHashIter iter;
          PetscBool     found;

          PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          ++Ner;
          eid  = EG_indexBodyTopo(body, edge);
          PetscCall(PetscHMapIFind(edgeMap, eid-1, &iter, &found));
          if (!found) {
            PetscInt v = numVertices - newVertices + numEdges;
            double range[4], params[3] = {0., 0., 0.}, result[18];
            int    periodic[2];

            PetscCall(PetscHMapISet(edgeMap, eid-1, numEdges++));
            PetscCall(EG_getRange(edge, range, periodic));
            params[0] = 0.5*(range[0] + range[1]);
            PetscCall(EG_evaluate(edge, params, result));
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

          PetscCall(EG_getBodyTopos(body, loop, FACE, &Nf, &fobjs));
          face = fobjs[0];
          fid  = EG_indexBodyTopo(body, face);
          if (Nf != 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Loop %d has %d faces, instead of 1 (%d)", lid-1, Nf, fid);
          PetscCall(EG_getRange(face, range, periodic));
          params[0] = 0.5*(range[0] + range[1]);
          params[1] = 0.5*(range[2] + range[3]);
          PetscCall(EG_evaluate(face, params, result));
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

      PetscCall(EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
        int lid, Ner = 0, Ne, e, nc = 0, c, Nt, t;

        lid  = EG_indexBodyTopo(body, loop);
        PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];
          int points[3];
          int eid, Nv, v, tmp;

          eid  = EG_indexBodyTopo(body, edge);
          PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
          if (mtype == DEGENERATE) continue;
          else                     ++Ner;
          if (Nv != 2) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Edge %d has %d vertices != 2", eid, Nv);

          for (v = 0; v < Nv; ++v) {
            ego vertex = nobjs[v];

            id = EG_indexBodyTopo(body, vertex);
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
      EG_free(lobjs);
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
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Context", (PetscObject) contextObj));
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

    PetscCall(EG_getBodyTopos(body, NULL, LOOP, &Nl, &lobjs));
    for (l = 0; l < Nl; ++l) {
      ego  loop = lobjs[l];
      ego *fobjs;
      int  lid, Nf, fid, Ner = 0, Ne, e, Nt = 0, t;

      lid  = EG_indexBodyTopo(body, loop);
      PetscCall(EG_getBodyTopos(body, loop, FACE, &Nf, &fobjs));
      if (Nf > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);
      fid  = EG_indexBodyTopo(body, fobjs[0]);
      EG_free(fobjs);
      PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses));
      for (e = 0; e < Ne; ++e) {
        ego             edge = objs[e];
        int             eid, Nv, v;
        PetscInt        points[3], support[2], numEdges, edgeNum;
        const PetscInt *edges;

        eid  = EG_indexBodyTopo(body, edge);
        PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));
        if (mtype == DEGENERATE) continue;
        else                     ++Ner;
        for (v = 0; v < Nv; ++v) {
          ego vertex = nobjs[v];

          id   = EG_indexBodyTopo(body, vertex);
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
    EG_free(lobjs);
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


static PetscErrorCode DMPlexCreateEGADS(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, vertexLabel;
  // EGADS variables
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
    PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));

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

      PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      PetscCall(EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));
      PetscCall(EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));
      EG_free(fobjs);
      EG_free(eobjs);
      EG_free(nobjs);

      // Remove DEGENERATE EDGES from Edge count
      PetscCall(EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));
      int Netemp = 0;
      for (int e = 0; e < Ne; ++e) {
        ego     edge = eobjs[e];
        int     eid;

        PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
        eid = EG_indexBodyTopo(body, edge);

        PetscCall(PetscHMapIFind(edgeMap, edgeCntr + eid - 1, &EMiter, &EMfound));
        if (mtype == DEGENERATE) {
          if (!EMfound) {PetscCall(PetscHMapISet(edgeMap, edgeCntr + eid - 1, -1));}
        }
        else {
          ++Netemp;
          if (!EMfound) {PetscCall(PetscHMapISet(edgeMap, edgeCntr + eid - 1, Netemp));}
        }
      }
	    edgeCntr += Ne;
      EG_free(eobjs);

      // Determine Number of Cells
      PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      for (int f = 0; f < Nf; ++f) {
        ego     face = fobjs[f];
        int     edgeTemp = 0;

        PetscCall(EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs));
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];

          PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
          if (mtype != DEGENERATE) {++edgeTemp;}
        }
        numCells += (2 * edgeTemp);
        EG_free(eobjs);
      }
      EG_free(fobjs);

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
      PetscCall(EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs));

      PetscCall(PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound));
      if (!BVfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);
      PetscCall(PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart));

      for (int v = 0; v < Nv; ++v) {
        ego    vertex = nobjs[v];
        double limits[4];
        int    id, dummy;

        PetscCall(EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses));
        id = EG_indexBodyTopo(body, vertex);
	
        coords[(bodyVertexIndexStart + id - 1)*cdim + 0] = limits[0];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 1] = limits[1];
        coords[(bodyVertexIndexStart + id - 1)*cdim + 2] = limits[2];
      }
      EG_free(nobjs);

      // Edge Midpoint Vertices on Current Body
      PetscCall(EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs));

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

        PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
        if (mtype == DEGENERATE) {continue;}

        eid = EG_indexBodyTopo(body, edge);

        // get relative offset from globalEdgeID Vector
        PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
        if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d not found in edgeMap", bodyEdgeGlobalIndexStart + eid - 1);
        PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

        PetscCall(EG_getRange(edge, range, &periodic));
        avgt[0] = (range[0] + range[1]) /  2.;

        PetscCall(EG_evaluate(edge, avgt, cntrPnt));
	  
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 2] = cntrPnt[2];
      }
      EG_free(eobjs);

      // Face Midpoint Vertices on Current Body
      PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));

      PetscCall(PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound));
      if (!BFfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);
      PetscCall(PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart));

      for (int f = 0; f < Nf; ++f) {
        ego       face = fobjs[f];
        double    range[4], avgUV[2], cntrPnt[18];
        int       peri, id;

        id = EG_indexBodyTopo(body, face);
        PetscCall(EG_getRange(face, range, &peri));

        avgUV[0] = (range[0] + range[1]) / 2.;
        avgUV[1] = (range[2] + range[3]) / 2.;
        PetscCall(EG_evaluate(face, avgUV, cntrPnt));

        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 1] = cntrPnt[1];
        coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 2] = cntrPnt[2];
      }
      EG_free(fobjs);

      // Define Cells :: Note - This could be incorporated in the Face Midpoint Vertices Loop but was kept separate for clarity
      PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
      for (int f = 0; f < Nf; ++f) {
        ego      face = fobjs[f];
        int      fID, midFaceID, midPntID, startID, endID, Nl;

        fID = EG_indexBodyTopo(body, face);
        midFaceID = numVertices + numEdges + bodyFaceIndexStart + fID - 1;
        // Must Traverse Loop to ensure we have all necessary information like the sense (+/- 1) of the edges.
        // TODO :: Only handles single loop faces (No holes). The choices for handling multiloop faces are:
        //            1) Use the DMPlexCreateEGADSFromFile() with the -dm_plex_egads_with_tess = 1 option.
        //               This will use a default EGADS tessellation as an initial surface mesh.
        //            2) Create the initial surface mesh via a 2D mesher :: Currently not availble (?future?)
        //               May I suggest the XXXX as a starting point?

        PetscCall(EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lSenses));

        if (Nl > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face has %d Loops. Can only handle Faces with 1 Loop. Please use --dm_plex_egads_with_tess = 1 Option", Nl);
        for (int l = 0; l < Nl; ++l) {
          ego      loop = lobjs[l];

          PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses));
          for (int e = 0; e < Ne; ++e) {
            ego     edge = eobjs[e];
            int     eid, eOffset;

            PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
            eid = EG_indexBodyTopo(body, edge);
            if (mtype == DEGENERATE) { continue; }

            // get relative offset from globalEdgeID Vector
            PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
            if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
            PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

            midPntID = numVertices + bodyEdgeIndexStart + eOffset - 1;

            PetscCall(EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses));

            if (eSenses[e] > 0) { startID = EG_indexBodyTopo(body, nobjs[0]); endID = EG_indexBodyTopo(body, nobjs[1]); }
            else { startID = EG_indexBodyTopo(body, nobjs[1]); endID = EG_indexBodyTopo(body, nobjs[0]); }

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
      EG_free(fobjs);
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
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Context", (PetscObject) contextObj));
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

    PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
    for (int f = 0; f < Nf; ++f) {
      ego   face = fobjs[f];
      int   fID, Nl;

      fID  = EG_indexBodyTopo(body, face);

      PetscCall(EG_getBodyTopos(body, face, LOOP, &Nl, &lobjs));
      for (int l = 0; l < Nl; ++l) {
        ego  loop = lobjs[l];
        int  lid;

        lid  = EG_indexBodyTopo(body, loop);
        if (Nl > 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);

        PetscCall(EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses));
        for (int e = 0; e < Ne; ++e) {
          ego     edge = eobjs[e];
          int     eid, eOffset;

          // Skip DEGENERATE Edges
          PetscCall(EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next));
          if (mtype == DEGENERATE) {continue;}
          eid = EG_indexBodyTopo(body, edge);

          // get relative offset from globalEdgeID Vector
          PetscCall(PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound));
          if (!EMfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);
          PetscCall(PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset));

          PetscCall(EG_getBodyTopos(body, edge, NODE, &Nv, &nobjs));
          for (int v = 0; v < Nv; ++v){
            ego vertex = nobjs[v];
            int vID;

            vID = EG_indexBodyTopo(body, vertex);
            PetscCall(DMLabelSetValue(bodyLabel, nStart + bodyVertexIndexStart + vID - 1, b));
            PetscCall(DMLabelSetValue(vertexLabel, nStart + bodyVertexIndexStart + vID - 1, vID));
          }
          EG_free(nobjs);

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
      EG_free(lobjs);

      PetscCall(DMLabelSetValue(bodyLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, b));
      PetscCall(DMLabelSetValue(faceLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, fID));
    }
    EG_free(fobjs);
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


static PetscErrorCode DMPlexCreateEGADS_Tess_Internal(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel              bodyLabel, faceLabel, edgeLabel, vertexLabel;
  /* EGADS variables */
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
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses));

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
    PetscCall(EG_getBoundingBox(body, boundBox));
    tessSize = boundBox[3] - boundBox[0];
    if (tessSize < boundBox[4] - boundBox[1]) tessSize = boundBox[4] - boundBox[1];
    if (tessSize < boundBox[5] - boundBox[2]) tessSize = boundBox[5] - boundBox[2];

    // TODO :: May want to give users tessellation parameter options //
    params[0] = 0.0250 * tessSize;
    params[1] = 0.0075 * tessSize;
    params[2] = 15.0;

    PetscCall(EG_makeTessBody(body, params, &tessArray[b]));

    PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));

    for (int f = 0; f < Nf; ++f) {
      ego             face = fobjs[f];
      int             len, fID, ntris;
      const int      *ptype, *pindex, *ptris, *ptric;
      const double   *pxyz, *puv;

      // Get Face ID //
      fID = EG_indexBodyTopo(body, face);

      // Checkout the Surface Tessellation //
      PetscCall(EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric));

      // Determine total number of triangle cells in the tessellation //
      bodyNumTris += (int) ntris;

      // Check out the point index and coordinate //
      for (int p = 0; p < len; ++p) {
        int global;

        PetscCall(EG_localToGlobal(tessArray[b], fID, p+1, &global));

        // Determine the total number of points in the tessellation //
        bodyNumPoints = PetscMax(bodyNumPoints, global);
      }
    }
    EG_free(fobjs);

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

    PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));

    for (int f = 0; f < Nf; ++f) {
      /* Get Face Object */
      ego              face = fobjs[f];
      int              len, fID, ntris;
      const int       *ptype, *pindex, *ptris, *ptric;
      const double    *pxyz, *puv;

      /* Get Face ID */
      fID = EG_indexBodyTopo(body, face);

      /* Checkout the Surface Tessellation */
      PetscCall(EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric));

      /* Check out the point index and coordinate */
      for (int p = 0; p < len; ++p) {
        int              global;
        PetscHashIter    PTLiter, PILiter, PBLiter;
        PetscBool        PTLfound, PILfound, PBLfound;

        PetscCall(EG_localToGlobal(tessArray[b], fID, p+1, &global));

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

        PetscCall(EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 0], &global));
        cells[(counter*3) +0] = global-1+pointIndexStart;

        PetscCall(EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 1], &globalA));
        cells[(counter*3) +1] = globalA-1+pointIndexStart;

        PetscCall(EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 2], &globalB));
        cells[(counter*3) +2] = globalB-1+pointIndexStart;

        PetscCall(PetscHMapIFind(triFaceIDLabelMap, counter, &TFLiter, &TFLfound));
        PetscCall(PetscHMapIFind(triBodyIDLabelMap, counter, &TBLiter, &TBLfound));

        if (!TFLfound)  {PetscCall(PetscHMapISet(triFaceIDLabelMap, counter, fID));}
        if (!TBLfound)  {PetscCall(PetscHMapISet(triBodyIDLabelMap, counter, b));}

        counter += 1;
      }
    }
    EG_free(fobjs);
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
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj));
    PetscCall(PetscContainerDestroy(&modelObj));

    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &contextObj));
    PetscCall(PetscContainerSetPointer(contextObj, context));
    PetscCall(PetscContainerSetUserDestroy(contextObj, DMPlexEGADSDestroy_Private));
    PetscCall(PetscObjectCompose((PetscObject) dm, "EGADS Context", (PetscObject) contextObj));
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
  DMPlexInflateToEGADSGeomModel - Snaps the vertex coordinates of a DMPlex object representing the mesh to its geometry if some vertices depart from the model. This usually happens with non-conforming refinement.

  Collective on dm

  Input Parameter:
. dm - The uninflated DM object representing the mesh

  Output Parameter:
. dm - The inflated DM object representing the mesh

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
PetscErrorCode DMPlexInflateToEGADSGeomModel(DM dm)
{
#if defined(PETSC_HAVE_EGADS)
  /* EGADS Variables */
  ego            model, geom, body, face, edge, vertex;
  ego           *bodies;
  int            Nb, oclass, mtype, *senses;
  double         result[4];
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
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
  if (!modelObj) PetscFunctionReturn(0);
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinates));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));

  PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));

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
    if (vertexID > 0) {
      PetscCall(EG_objectBodyTopo(body, NODE, vertexID, &vertex));
      PetscCall(EG_evaluate(vertex, NULL, result));
      //PetscCall(EG_getTopology(vertex, &geom, &oclass, &mtype, result, NULL, NULL, NULL));
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (edgeID > 0) {
      /* Snap to EDGE at nearest location */
      double params[1];
      PetscCall(EG_objectBodyTopo(body, EDGE, edgeID, &edge));
      PetscCall(EG_invEvaluate(edge, vcoords, params, result)); // Get (x,y,z) of nearest point on EDGE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (faceID > 0) {
      /* Snap to FACE at nearest location */
      double params[2];
      PetscCall(EG_objectBodyTopo(body, FACE, faceID, &face));
      PetscCall(EG_invEvaluate(face, vcoords, params, result)); // Get (x,y,z) of nearest point on FACE
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    }
  }
  PetscCall(VecRestoreArrayWrite(coordinates, &coords));
  /* Clear out global coordinates */
  PetscCall(VecDestroy(&dm->coordinates));
#endif
  PetscFunctionReturn(0);
}


/*@
  DMPlexInflateToGeomModel - Snaps the vertex coordinates of a DMPlex object representing the mesh to its geometry if some vertices depart from the model. This usually happens with non-conforming refinement.

  Collective on dm

  Input Parameter:
. dm - The uninflated DM object representing the mesh

  Output Parameter:
. dm - The inflated DM object representing the mesh

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
PetscErrorCode DMPlexInflateToGeomModel(DM dm)
{
#if defined(PETSC_HAVE_EGADS)
  /* PETSc Variables */
  PetscContainer modelObj;
  PetscBool      islite = PETSC_FALSE;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_EGADS)
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
  if (!modelObj) {
    PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSlite Model", (PetscObject *) &modelObj));
    islite = PETSC_TRUE;
  }
  
  if (modelObj) {
    if (islite) {PetscCall(DMPlexInflateToEGADSliteGeomModel(dm));}
    else        {PetscCall(DMPlexInflateToEGADSGeomModel(dm));}
  }
#endif
  PetscFunctionReturn(0);
}


PetscErrorCode ConvertEGADSModelToAllBSplines(ego model)
{
  ego	           context = NULL, geom, *bodies, *fobjs;
  int            oclass, mtype;
  int	          *senses;
  int	           Nb, Nf;
  //PetscErrorCode ierr;

  // Get the number of bodies and body objects in the model
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  //PetscCall(PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", Nb));
  
  // Get all Faces on the body    <-- Only working with 1 body at the moment.
  ego body = bodies[0];
  PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
  ego newFaces[Nf];

  // Convert the 1st Face to a BSpline Geometry
  for (int ii = 0; ii < Nf; ++ii) {
    ego     face = fobjs[ii];		// ii was 0 when looking at a single face
    ego     gRef, gPrev, gNext, *lobjs;
    int     goclass, gmtype, *gpinfo;
	  int     Nl, *lsenses;		// orig --> id
    double *gprv;
    char   *gClass=(char*)"", *gType=(char*)"";
  
    //id   = EG_indexBodyTopo(body, face);			// get FACE ID
    PetscCall(EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses));	// Get FACES Geometry object (geom_
    PetscCall(EG_getGeometry(geom, &goclass, &gmtype, &gRef, &gpinfo, &gprv));				// Get geometry object info
    PetscCall(EG_getInfo(geom, &goclass, &gmtype, &gRef, &gPrev, &gNext));					// Get geometry info

    DMPlex_EGADS_GeomDecode_Internal(goclass, gmtype, &gClass, &gType);   // Decode Geometry integers
    
    // Convert current FACE to a BSpline Surface
    ego     bspline;
    ego     bRef, bPrev, bNext;
    int     boclass, bmtype, *bpinfo;
    double *bprv;
    char   *bClass=(char*)"", *bType=(char*)"";

    PetscCall(EG_convertToBSpline(face, &bspline));
  
    PetscCall(EG_getGeometry(bspline, &boclass, &bmtype, &bRef, &bpinfo, &bprv));				// Get geometry object info
    PetscCall(EG_getInfo(bspline, &boclass, &bmtype, &bRef, &bPrev, &bNext));					// Get geometry info
  
    DMPlex_EGADS_GeomDecode_Internal(boclass, bmtype, &bClass, &bType);   // Decode Geometry integers

    // Get Context from FACE
    context = NULL;
    PetscCall(EG_getContext(face, &context));
	
	  // Silent WARNING Regarding OPENCASCADE 7.5
	  EG_setOutLevel(context, 0);
  
    ego newgeom;
    PetscCall(EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bpinfo, bprv, &newgeom));
  
    // Create new FACE based on new SURFACE geometry
    double data[4];
    int    periodic;
    PetscCall(EG_getRange(newgeom, data, &periodic));

    ego newface;
    PetscCall(EG_makeFace(newgeom, SFORWARD, data, &newface));
	
	  // Reinstate WARNING Regarding OPENCASCADE 7.5
	  PetscCall(EG_setOutLevel(context, 1));
	
	  newFaces[ii] = newface;
  }
  
  // Sew New Faces together to get a new model
  ego newmodel;
  PetscCall(EG_sewFaces(Nf, newFaces, 0.0, 0, &newmodel));
  *model = *newmodel;
  PetscFunctionReturn(0);
}


/*@C
  DMPlexCreateEGADSFromFile - Create a DMPlex mesh from an EGADS, IGES, or STEP file.

  Collective

  Input Parameters:
+ comm     - The MPI communicator
- filename - The name of the EGADS, IGES, or STEP file

  Output Parameter:
. dm       - The DM object representing the mesh

  Level: beginner

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSliteFromFile()
@*/
PetscErrorCode DMPlexCreateEGADSFromFile(MPI_Comm comm, const char filename[], DM *dm)
{
  PetscMPIInt    rank;
#if defined(PETSC_HAVE_EGADS)
  ego            context= NULL, model = NULL;
#endif
  PetscBool      printModel = PETSC_FALSE, tessModel = PETSC_FALSE, newModel = PETSC_FALSE;
  PetscBool      shapeOpt = PETSC_FALSE;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidCharPointer(filename, 2);
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_print_model", &printModel, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_tess_model", &tessModel, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_new_model", &newModel, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_shape_opt", &shapeOpt, NULL));
  PetscCall(MPI_Comm_rank(comm, &rank));
#if defined(PETSC_HAVE_EGADS)
  if (!rank) {
	  
    PetscCall(EG_open(&context));
    PetscCall(EG_loadModel(context, 0, filename, &model));
    if (shapeOpt)   {PetscCall(ConvertEGADSModelToAllBSplines(model));}
    if (printModel) {PetscCall(DMPlexEGADSPrintModel_Internal(model));}

  }
  if (tessModel)     {PetscCall(DMPlexCreateEGADS_Tess_Internal(comm, context, model, dm));}
  else if (newModel) {PetscCall(DMPlexCreateEGADS_Internal(comm, context, model, dm));}
  else               {PetscCall(DMPlexCreateEGADS(comm, context, model, dm));}
  PetscFunctionReturn(0);
#else
  SETERRQ(comm, PETSC_ERR_SUP, "This method requires EGADS support. Reconfigure using --download-egads");
#endif
}


/*@C
  DMPlex_Surface_Grad - Exposes the Geometry's Control Points and Weights and Calculates the Mesh Topology Boundary Nodes Gradient
                        with respect the associated geometry's Control Points and Weights.

  Collective

  Input Parameters:
. dm      - The DM object representing the mesh with PetscContainer containing an EGADS geometry model


  Output Parameter:
. dm       - The DM object representing the mesh with PetscContainers containing the EGADS geometry model,
             Array-Hash Table Geometry Control Point Pair, Array-Hash Table Geometry Weights Pair and
			 Matrix-Hash Table Surface Gradient Pair

  Level: intermediate

.seealso:
@*/
PetscErrorCode DMPlex_Surface_Grad(DM dm)
{
	ego 			      model, geom, *bodies, *fobjs;
	PetscContainer 	modelObj;
	int            	oclass, mtype, *senses;
	int            	Nb, Nf;
	PetscHMapI      faceCntrlPtRow_Start = NULL, faceCPWeightsRow_Start = NULL;
	PetscHMapI      pointSurfGradRow_Start = NULL;
	Mat             pointSurfGrad;
  IS			        faceLabelValues, edgeLabelValues, vertexLabelValues;
	PetscInt	      faceLabelSize, edgeLabelSize, vertexLabelSize;
	//PetscErrorCode 	ierr;
	
	PetscFunctionBegin;
	PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
	if (!modelObj) {
	  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSlite Model", (PetscObject *) &modelObj));
	}

	// Get attached EGADS model (pointer)
	PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
	
	// Get the bodies in the model
	PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
	ego body = bodies[0];		// Only operate on 1st body. Model should only have 1 body.
	
	// Get the total number of FACEs in the model
	PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
	
	// Get the total number of points and IDs in the DMPlex with a "EGADS Face Label"
	// This will provide the total number of DMPlex points on the boundary of the geometry
	PetscCall(DMGetLabelIdIS(dm, "EGADS Face ID", &faceLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Face ID", &faceLabelSize));
	
	PetscCall(DMGetLabelIdIS(dm, "EGADS Edge ID", &edgeLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Edge ID", &edgeLabelSize));
	
	PetscCall(DMGetLabelIdIS(dm, "EGADS Vertex ID", &vertexLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Vertex ID", &vertexLabelSize));
	
	const PetscInt	*faceIndices, *edgeIndices, *vertexIndices;
	PetscCall(ISGetIndices(faceLabelValues, &faceIndices));
	PetscCall(ISGetIndices(edgeLabelValues, &edgeIndices));
	PetscCall(ISGetIndices(vertexLabelValues, &vertexIndices));
		
	// Get the points associated with each FACE, EDGE and VERTEX label in the DM
	PetscInt	totalNumPoints = 0;
	for (int ii = 0; ii < faceLabelSize; ++ii) {
	  // Cycle through FACE labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Face ID", faceIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(faceLabelValues, &faceIndices));
	PetscCall(ISDestroy(&faceLabelValues));
	
	for (int ii = 0; ii < edgeLabelSize; ++ii) {
	  // Cycle Through EDGE Labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Edge ID", edgeIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(edgeLabelValues, &edgeIndices));
	PetscCall(ISDestroy(&edgeLabelValues));
	
	for (int ii = 0; ii < vertexLabelSize; ++ii) {
	  // Cycle Through VERTEX Labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Vertex ID", vertexIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(vertexLabelValues, &vertexIndices));
	PetscCall(ISDestroy(&vertexLabelValues));

	int      maxNumCPs = 0;
	int      totalNumCPs = 0;
	ego 	   bRef, bPrev, bNext, fgeom, *lobjs;
	//ego		   bspline;
	int 	   id, boclass, bmtype, *bpinfo;
	int		   foclass, fmtype, Nl, *lsenses;
	double 	*bprv;
	double	 fdata[4];
	
	// Create Hash Tables
	PetscInt   cntr = 0, wcntr = 0;
	PetscCall(PetscHMapICreate(&faceCntrlPtRow_Start));
	PetscCall(PetscHMapICreate(&faceCPWeightsRow_Start));
	
	for (int ii = 0; ii < Nf; ++ii) {
	  // Need to get the maximum number of Control Points defining the FACEs
	  ego		face = fobjs[ii];
	  //char 	*bClass=(char*)"", *bType=(char*)"";
	  int   maxNumCPs_temp;
  
	  id   = EG_indexBodyTopo(body, face);
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	  PetscCall(EG_getInfo(fgeom, &boclass, &bmtype, &bRef, &bPrev, &bNext));
	  maxNumCPs_temp = bpinfo[2] * bpinfo[5];
	  totalNumCPs += bpinfo[2] * bpinfo[5];
	  
	  if (maxNumCPs_temp > maxNumCPs) {maxNumCPs = maxNumCPs_temp;}
	}

  PetscInt    *cpCoordDataLengthPtr, *wDataLengthPtr;
	PetscInt     cpCoordDataLength = 3 * totalNumCPs;
	PetscInt     wDataLength = totalNumCPs;
	cpCoordDataLengthPtr = &cpCoordDataLength;
	wDataLengthPtr = &wDataLength;
	//PetscScalar cntrlPtCoords[cpCoordDataLength];
	//PetscScalar cntrlPtWeights[wDataLength];
	PetscScalar *cntrlPtCoords, *cntrlPtWeights;
	PetscMalloc1(cpCoordDataLength, &cntrlPtCoords);
	PetscMalloc1(wDataLength, &cntrlPtWeights);
	for (int ii = 0; ii < Nf; ++ii) {
	  // Need to Populate Control Point Coordinates and Weight Vectors
	  ego		        face = fobjs[ii];
	  //char 	*bClass=(char*)"", *bType=(char*)"";
	  //int     maxNumCPs_temp;
	  PetscHashIter hashKeyIter, wHashKeyIter;
	  PetscBool   	hashKeyFound, wHashKeyFound;
  
	  id   = EG_indexBodyTopo(body, face);
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	  PetscCall(EG_getInfo(fgeom, &boclass, &bmtype, &bRef, &bPrev, &bNext));
		
	  // Store Face ID to 1st Row of Control Point Vector
	  PetscCall(PetscHMapIFind(faceCntrlPtRow_Start, id, &hashKeyIter, &hashKeyFound));
		
	  if (!hashKeyFound)  {
        PetscCall(PetscHMapISet(faceCntrlPtRow_Start, id, cntr));
	  }
		
	  int offsetCoord = bpinfo[3] + bpinfo[6];
	  for (int jj = 0; jj < 3 * bpinfo[2] * bpinfo[5]; ++jj) {
	    cntrlPtCoords[cntr] = bprv[offsetCoord + jj];
	    cntr += 1;
	  }
		
	  // Store Face ID to 1st Row of Control Point Weight Vector
	  PetscCall(PetscHMapIFind(faceCPWeightsRow_Start, id, &wHashKeyIter, &wHashKeyFound));
		
	  if (!wHashKeyFound)  {
	    PetscCall(PetscHMapISet(faceCPWeightsRow_Start, id, wcntr));
	  }
		
	  int offsetWeight = bpinfo[3] + bpinfo[6] + (3 * bpinfo[2] * bpinfo[5]);
	  for (int jj = 0; jj < bpinfo[2] *bpinfo[5]; ++jj) {
	    cntrlPtWeights[wcntr] = bprv[offsetWeight + jj];
		  wcntr += 1;
	  }
	}
	
	
	// Attach Control Point and Weight Data to DM
	{
	  PetscContainer cpOrgObj, cpCoordObj, cpCoordLengthObj;
	  PetscContainer wOrgObj, wValObj, wDataLengthObj;

	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpOrgObj));
	  PetscCall(PetscContainerSetPointer(cpOrgObj, faceCntrlPtRow_Start));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Hash Table", (PetscObject) cpOrgObj));
	  PetscCall(PetscContainerDestroy(&cpOrgObj));
	  
	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpCoordObj));
	  PetscCall(PetscContainerSetPointer(cpCoordObj, cntrlPtCoords));
	  //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Coordinates", (PetscObject) cpCoordObj));
	  PetscCall(PetscContainerDestroy(&cpCoordObj));
		
	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpCoordLengthObj));
	  PetscCall(PetscContainerSetPointer(cpCoordLengthObj, cpCoordDataLengthPtr));
	  //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Coordinate Data Length", (PetscObject) cpCoordLengthObj));
	  PetscCall(PetscContainerDestroy(&cpCoordLengthObj));

		
	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wOrgObj));
	  PetscCall(PetscContainerSetPointer(wOrgObj, faceCPWeightsRow_Start));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weights Hash Table", (PetscObject) wOrgObj));
	  PetscCall(PetscContainerDestroy(&wOrgObj));

	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wValObj));
	  PetscCall(PetscContainerSetPointer(wValObj, cntrlPtWeights));
	  //PetscCall(PetscContainerSetUserDestroy(wValObj, DMPlexEGADSDestroy_Private));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight Data", (PetscObject) wValObj));
	  PetscCall(PetscContainerDestroy(&wValObj));

	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wDataLengthObj));
	  PetscCall(PetscContainerSetPointer(wDataLengthObj, wDataLengthPtr));
	  //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight Data Length", (PetscObject) wDataLengthObj));
	  PetscCall(PetscContainerDestroy(&wDataLengthObj));	

	}
	
	// Define Matrix to store  Surface Gradient information dx_i/dCPj_i
	PetscInt		    gcntr = 0;
	const PetscInt  rowSize = 3 * maxNumCPs * totalNumPoints;
	const PetscInt  colSize = 4 * Nf;
	//PetscScalar		surfaceGrad[rowSize][colSize];
	
	// Create Point Surface Gradient Matrix
	MatCreate(PETSC_COMM_WORLD, &pointSurfGrad);
  MatSetSizes(pointSurfGrad,PETSC_DECIDE, PETSC_DECIDE, rowSize, colSize);
  MatSetType(pointSurfGrad, MATAIJ);
  MatSetUp(pointSurfGrad);
	
	// Create Hash Table to store Point's stare row in surfaceGrad[][]
	PetscCall(PetscHMapICreate(&pointSurfGradRow_Start));
	
	// Get Coordinates for the DMPlex point
	DM 			       cdm;
	PetscInt	     dE, Nv;
	Vec            coordinatesLocal;
	PetscScalar   *coords = NULL;
	PetscCall(DMGetCoordinateDM(dm, &cdm));
	PetscCall(DMGetCoordinateDim(dm, &dE));
	PetscCall(DMGetCoordinatesLocal(dm, &coordinatesLocal));
	
	// CYCLE THROUGH FACEs
	for (int ii = 0; ii < Nf; ++ii) {
	  ego		           face = fobjs[ii];
	  ego             *eobjs, *nobjs;
	  PetscInt	       fid, Ne, Nn;	
	  DMLabel		       faceLabel, edgeLabel, nodeLabel;
	  PetscHMapI       currFaceUniquePoints = NULL;
	  IS			         facePoints, edgePoints, nodePoints;
	  const PetscInt	*fIndices, *eIndices, *nIndices;
	  PetscInt         fSize, eSize, nSize;
	  PetscHashIter    fHashKeyIter, eHashKeyIter, nHashKeyIter, pHashKeyIter;
	  PetscBool   	   fHashKeyFound, eHashKeyFound, nHashKeyFound, pHashKeyFound;
	  PetscInt	       cfCntr = 0;
		
	  // Get Geometry Object for the Current FACE
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	
	  // Get all EDGE and NODE objects attached to the current FACE
	  PetscCall(EG_getBodyTopos(body, face, EDGE,  &Ne, &eobjs));
	  PetscCall(EG_getBodyTopos(body, face, NODE,  &Nn, &nobjs));
		
	  // Get all DMPlex Points that have DMLabel "EGADS Face ID" and store them in a Hash Table for later use
	  fid  = EG_indexBodyTopo(body, face);
	  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
	  PetscCall(DMLabelGetStratumIS(faceLabel, fid, &facePoints));
	  PetscCall(ISGetIndices(facePoints, &fIndices));
	  PetscCall(ISGetSize(facePoints, &fSize));
		
	  PetscCall(PetscHMapICreate(&currFaceUniquePoints));
	  
	  for (int jj = 0; jj < fSize; ++jj) {
	    PetscCall(PetscHMapIFind(currFaceUniquePoints, fIndices[jj], &fHashKeyIter, &fHashKeyFound));
			
		  if (!fHashKeyFound) {
		    PetscCall(PetscHMapISet(currFaceUniquePoints, fIndices[jj], cfCntr));
	      cfCntr += 1;
		  }
			
		  PetscCall(PetscHMapIFind(pointSurfGradRow_Start, fIndices[jj], &pHashKeyIter, &pHashKeyFound));
			
		  if (!pHashKeyFound) {
		    PetscCall(PetscHMapISet(pointSurfGradRow_Start, fIndices[jj], gcntr));
		    gcntr += 3 * maxNumCPs;
		  }
	  }
	  PetscCall(ISRestoreIndices(facePoints, &fIndices));
	  PetscCall(ISDestroy(&facePoints));
		
	  // Get all DMPlex Points that have DMLable "EGADS Edge ID" attached to the current FACE and store them in a Hash Table for later use.
	  for (int jj = 0; jj < Ne; ++jj) {
      ego			edge = eobjs[jj];
		  PetscBool	containLabelValue;
			
		  id = EG_indexBodyTopo(body, edge);
		  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
		  PetscCall(DMLabelHasValue(edgeLabel, id, &containLabelValue));
			
		  if (containLabelValue) {
		    PetscCall(DMLabelGetStratumIS(edgeLabel, id, &edgePoints));
		    PetscCall(ISGetIndices(edgePoints, &eIndices));
	      PetscCall(ISGetSize(edgePoints, &eSize));
			
		    for (int kk = 0; kk < eSize; ++kk) {
		      PetscCall(PetscHMapIFind(currFaceUniquePoints, eIndices[kk], &eHashKeyIter, &eHashKeyFound));
			
			    if (!eHashKeyFound) {
			      PetscCall(PetscHMapISet(currFaceUniquePoints, eIndices[kk], cfCntr));
			      cfCntr += 1;
			    }
					
			    PetscCall(PetscHMapIFind(pointSurfGradRow_Start, eIndices[kk], &pHashKeyIter, &pHashKeyFound));
			
			    if (!pHashKeyFound) {
			      PetscCall(PetscHMapISet(pointSurfGradRow_Start, eIndices[kk], gcntr));
			      gcntr += 3 * maxNumCPs;
			    }
		    }
		    PetscCall(ISRestoreIndices(edgePoints, &eIndices));
		    PetscCall(ISDestroy(&edgePoints));
		  }
	  }
		
	  // Get all DMPlex Points that have DMLabel "EGADS Vertex ID" attached to the current FACE and store them in a Hash Table for later use.
	  for (int jj = 0; jj < Nn; ++jj) {
	    ego		node = nobjs[jj];
		  id = EG_indexBodyTopo(body, node);
		  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &nodeLabel));
		  PetscCall(DMLabelGetStratumIS(nodeLabel, id, &nodePoints));
		  PetscCall(ISGetIndices(nodePoints, &nIndices));
		  PetscCall(ISGetSize(nodePoints, &nSize));
			
		  for (int kk = 0; kk < nSize; ++kk) {
		    PetscCall(PetscHMapIFind(currFaceUniquePoints, nIndices[kk], &nHashKeyIter, &nHashKeyFound));
			
		    if (!nHashKeyFound) {
		      PetscCall(PetscHMapISet(currFaceUniquePoints, nIndices[kk], cfCntr));
			    cfCntr += 1;
		    }
				
		    PetscCall(PetscHMapIFind(pointSurfGradRow_Start, nIndices[kk], &pHashKeyIter, &pHashKeyFound));	
		    if (!pHashKeyFound) {
			    PetscCall(PetscHMapISet(pointSurfGradRow_Start, nIndices[kk], gcntr));
			    gcntr += 3 * maxNumCPs;
		    }
		  }
		  PetscCall(ISRestoreIndices(nodePoints, &nIndices));
		  PetscCall(ISDestroy(&nodePoints));
	  }
	
	  // Get the Total Number of entries in the Hash Table
	  PetscInt	currFaceUPSize;
	  PetscCall(PetscHMapIGetSize(currFaceUniquePoints, &currFaceUPSize));
		
	  // Get Keys
	  PetscInt	currFaceUPKeys[currFaceUPSize], off = 0;
	  PetscCall(PetscHMapIGetKeys(currFaceUniquePoints, &off, currFaceUPKeys));
		
	  // Cycle through all points on the current FACE
	  for (int jj = 0; jj < currFaceUPSize; ++jj) {
	    PetscInt	currPointID = currFaceUPKeys[jj];
		  PetscCall(DMPlexVecGetClosure(cdm, NULL, coordinatesLocal, currPointID, &Nv, &coords));

		  // Get UV position of FACE
		  double	params[2], range[4], eval[18]; //, eval2[18], paramsV[4], result[18];
		  int	    peri; 
		  PetscCall(EG_getRange(face, range, &peri));
		  PetscCall(DMPlex_EGADS_FACE_XYZtoUV_Internal(coords, face, range, 0, dE, params));	
		  PetscCall(EG_evaluate(face, params, eval));	
	
		  // Make a new SURFACE Geometry by changing the location of the Control Points
		  int    prvSize = bpinfo[3] + bpinfo[6] + (4 * bpinfo[2]*bpinfo[5]);
		  double nbprv[prvSize];
			
		  // Cycle through each Control Point
		  double  deltaCoord = 1.0E-4;
		  int     offset = bpinfo[3] + bpinfo[6];
		  int     wOffset = offset + (3 * bpinfo[2] * bpinfo[5]);
		  for (int ii = 0; ii < bpinfo[2]*bpinfo[5]; ++ii){
		    // Cycle through each direction (x, then y, then z)
		    for (int kk = 0; kk < 4; ++kk) {
		      // Reinitialize nbprv[] values because we only want to change one value at a time
		      for (int mm = 0; mm < prvSize; ++mm) {
			      nbprv[mm] = bprv[mm];
			    }
					
			    if (kk == 0) {	        //X
			      nbprv[offset + 0] = bprv[offset + 0] + deltaCoord;
			      nbprv[offset + 1] = bprv[offset + 1];
			      nbprv[offset + 2] = bprv[offset + 2];
			    } else if (kk == 1) {	//Y
			      nbprv[offset + 0] = bprv[offset + 0];
			      nbprv[offset + 1] = bprv[offset + 1] + deltaCoord;
			      nbprv[offset + 2] = bprv[offset + 2];
			    } else if (kk == 2) {	//Z
			      nbprv[offset + 0] = bprv[offset + 0];
			      nbprv[offset + 1] = bprv[offset + 1];
			      nbprv[offset + 2] = bprv[offset + 2] + deltaCoord;
			    } else if (kk == 3) {	// Weights
			      nbprv[wOffset + ii] = bprv[wOffset + ii] + deltaCoord;
			    } else {
			     // currently do nothing
			    }		
					
			    // Create New Surface Based on New Control Points or Weights
			    ego newgeom, context;
			    PetscCall(EG_setOutLevel(context, 0));
			    PetscCall(EG_open(&context));
			    PetscCall(EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bpinfo, nbprv, &newgeom));
			    PetscCall(EG_setOutLevel(context, 1));
					
			    // Evaluate new (x, y, z) Point Position based on new Surface Definition
			    double newCoords[18]; //, newParams[2];
			    //double newRange[4];
			    PetscCall(EG_getRange(newgeom, range, &peri));
			    PetscCall(DMPlex_EGADS_FACE_XYZtoUV_Internal(coords, newgeom, range, 0, dE, params));
			    PetscCall(EG_evaluate(newgeom, params, newCoords));

			    // Now Calculate the Surface Gradient for the change in x-component Control Point
			    PetscScalar dxdCx = (newCoords[0] - coords[0]) / deltaCoord;
			    PetscScalar dxdCy = (newCoords[1] - coords[1]) / deltaCoord;
			    PetscScalar dxdCz = (newCoords[2] - coords[2]) / deltaCoord;
					
			    // Store Gradient Information in surfaceGrad[][] Matrix
			    PetscInt   startRow;
			    PetscCall(PetscHMapIGet(pointSurfGradRow_Start, currPointID, &startRow));
/*
			    surfaceGrad[startRow + (ii * 3) + 0][((fid - 1) * 4) + kk] = dxdCx;
			    surfaceGrad[startRow + (ii * 3) + 1][((fid - 1) * 4) + kk] = dxdCy;
			    surfaceGrad[startRow + (ii * 3) + 2][((fid - 1) * 4) + kk] = dxdCz;
*/		
			    // Store Results in Petsc Matrix
			    PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 0, ((fid - 1) * 4) + kk, dxdCx, INSERT_VALUES));
			    PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 1, ((fid - 1) * 4) + kk, dxdCy, INSERT_VALUES));
			    PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 2, ((fid - 1) * 4) + kk, dxdCz, INSERT_VALUES));
		    }
		   offset += 3;
		  }
		  PetscCall(DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, currPointID, &Nv, &coords));
	  }
	}
	
	// Assemble Point Surface Grad Matrix
	MatAssemblyBegin(pointSurfGrad, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pointSurfGrad, MAT_FINAL_ASSEMBLY);
	
	// Attach Surface Gradient Hash Table and Matrix to DM
	{
	  PetscContainer surfGradOrgObj, surfGradObj;

	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &surfGradOrgObj));
	  PetscCall(PetscContainerSetPointer(surfGradOrgObj, pointSurfGradRow_Start));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Gradient Hash Table", (PetscObject) surfGradOrgObj));
	  PetscCall(PetscContainerDestroy(&surfGradOrgObj));

	  PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &surfGradObj));
	  PetscCall(PetscContainerSetPointer(surfGradObj, pointSurfGrad));
	  //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
	  PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Gradient Matrix", (PetscObject) surfGradObj));
	  PetscCall(PetscContainerDestroy(&surfGradObj));
	}
	EG_free(fobjs);
	PetscFunctionReturn(0);
}


/*@C
  DMPlexGeomDataAndGrads - Exposes Control Points and Control Point Weights defining the underlying geometry allowing user manipulation of the geometry.
                           Calculates the DM Point location, surface area and volume gradients wrt to Control Point and Control Point Weights using Finite
                           Difference (small perturbation of Control Point coordinates or Control Point Weight value).

  Collective

  Input Parameters:
. dm            - The DM object representing the mesh with PetscContainer containing an EGADS geometry model
. fullGeomGrad  - PetscBool flag. Determines how the Surface Area and Volume Gradients wrt to Control Points and Control Point Weights are calculated.
                      PETSC_FALSE :: Surface Area Gradient wrt Control Points and Control Point Weights are calculated using the change in the local
                                     FACE changes (not the entire body). Volume Gradients are not calculated. Faster computations.
                      PETSC_TRUE  :: Surface Area Gradietn wrt to Control Points and Control Point Weights are calculated using the change observed in
                                     the entire solid body. Volume Gradients are calculated. Slower computation due to the need to generate a new solid
                                     body geometry for every Control Point and Control Point Weight change.


  Output Parameter:
. dm       - The updated DM object representing the mesh with PetscContainers containing the Control Point, Control Point Weight and Gradient Data.

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSliteFromFile(), DMPlexModifyEGADSGeomModel()
@*/
PetscErrorCode DMPlexGeomDataAndGrads(DM dm, PetscBool fullGeomGrad) 
{
	ego 			      model, geom, *bodies, *fobjs;
	PetscContainer 	modelObj;
	int            	oclass, mtype, *senses;
	int            	Nb, Nf;
	PetscHMapI      faceCntrlPtRow_Start = NULL, faceCPWeightsRow_Start = NULL;
	PetscHMapI      pointSurfGradRow_Start = NULL;
	Mat             pointSurfGrad, cpEquiv;
  IS			        faceLabelValues, edgeLabelValues, vertexLabelValues;
	PetscInt	      faceLabelSize, edgeLabelSize, vertexLabelSize;
	//PetscErrorCode 	ierr;
	
	//PetscFunctionBegin;
	PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
	if (!modelObj) {
	  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSlite Model", (PetscObject *) &modelObj));
	}
  
  // Throw Error is DM does not have an attached EGADS geometry model
  if (!modelObj) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DM does not have required EGADS Geometry Model attached. Please generate DM with a Geometry Model attached.");

	// Get attached EGADS model (pointer)
	PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
	
	// Get the bodies in the model
	PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
	ego body = bodies[0];		// Only operate on 1st body. Model should only have 1 body.
	
	// Get the total number of FACEs in the model
	PetscCall(EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs));
	
	// Get the total number of points and IDs in the DMPlex with a "EGADS Face Label"
	// This will provide the total number of DMPlex points on the boundary of the geometry
	PetscCall(DMGetLabelIdIS(dm, "EGADS Face ID", &faceLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Face ID", &faceLabelSize));
	
	PetscCall(DMGetLabelIdIS(dm, "EGADS Edge ID", &edgeLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Edge ID", &edgeLabelSize));
	
	PetscCall(DMGetLabelIdIS(dm, "EGADS Vertex ID", &vertexLabelValues));
	PetscCall(DMGetLabelSize(dm, "EGADS Vertex ID", &vertexLabelSize));
	
	const PetscInt	*faceIndices, *edgeIndices, *vertexIndices;
	PetscCall(ISGetIndices(faceLabelValues, &faceIndices));
	PetscCall(ISGetIndices(edgeLabelValues, &edgeIndices));
	PetscCall(ISGetIndices(vertexLabelValues, &vertexIndices));
		
	// Get the points associated with each FACE, EDGE and VERTEX label in the DM
	PetscInt	totalNumPoints = 0;
	for (int ii = 0; ii < faceLabelSize; ++ii) {
	  // Cycle through FACE labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Face ID", faceIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(faceLabelValues, &faceIndices));
	PetscCall(ISDestroy(&faceLabelValues));
	
	for (int ii = 0; ii < edgeLabelSize; ++ii) {
	  // Cycle Through EDGE Labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Edge ID", edgeIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(edgeLabelValues, &edgeIndices));
	PetscCall(ISDestroy(&edgeLabelValues));
	
	for (int ii = 0; ii < vertexLabelSize; ++ii) {
	  // Cycle Through VERTEX Labels
	  PetscInt	size;
	  PetscCall(DMGetStratumSize(dm, "EGADS Vertex ID", vertexIndices[ii], &size));
	  totalNumPoints += size;
	}
	PetscCall(ISRestoreIndices(vertexLabelValues, &vertexIndices));
	PetscCall(ISDestroy(&vertexLabelValues));

	int     maxNumCPs = 0;
	int     totalNumCPs = 0;
	ego 	  bRef, bPrev, bNext, fgeom, *lobjs;
	int 	  id, boclass, bmtype, *bpinfo;
	int		  foclass, fmtype, Nl, *lsenses;
	double *bprv;
	double	fdata[4];
	
	// Create Hash Tables
	PetscInt   cntr = 0, wcntr = 0, vcntr = 0;
	PetscCall(PetscHMapICreate(&faceCntrlPtRow_Start));
	PetscCall(PetscHMapICreate(&faceCPWeightsRow_Start));
	
	for (int ii = 0; ii < Nf; ++ii) {
	  // Need to get the maximum number of Control Points defining the FACEs
	  ego		face = fobjs[ii];
	  int   maxNumCPs_temp;
  
	  id   = EG_indexBodyTopo(body, face);
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	  PetscCall(EG_getInfo(fgeom, &boclass, &bmtype, &bRef, &bPrev, &bNext));
	  maxNumCPs_temp = bpinfo[2] * bpinfo[5];
	  totalNumCPs += bpinfo[2] * bpinfo[5];
	  
	  if (maxNumCPs_temp > maxNumCPs) {maxNumCPs = maxNumCPs_temp;}
	}

  PetscInt    *cpCoordDataLengthPtr, *wDataLengthPtr;
	PetscInt     cpCoordDataLength = 3 * totalNumCPs;
	PetscInt     wDataLength = totalNumCPs;
	cpCoordDataLengthPtr = &cpCoordDataLength;
	wDataLengthPtr = &wDataLength;
  
	PetscScalar *cntrlPtCoords, *cntrlPtWeights;
	PetscMalloc1(cpCoordDataLength, &cntrlPtCoords);
	PetscMalloc1(wDataLength, &cntrlPtWeights);

  // For dSA/dCPi 
  PetscScalar *gradSACP, *gradSAW, *gradVCP, *gradVW;
  PetscMalloc1(cpCoordDataLength, &gradSACP);
  PetscMalloc1(wDataLength, &gradSAW);
  PetscMalloc1(cpCoordDataLength, &gradVCP);
  PetscMalloc1(wDataLength, &gradVW);
  
  // Control Point - Vertex/Edge/Face Relationship
  PetscInt  *cp_vertex, *cp_edge, *cp_face;
  PetscInt  *w_vertex, *w_edge, *w_face;
  PetscMalloc3(totalNumCPs, &cp_vertex, totalNumCPs, &cp_edge, totalNumCPs, &cp_face);
  PetscMalloc3(wDataLength, &w_vertex, wDataLength, &w_edge, wDataLength, &w_face);

	for (int ii = 0; ii < Nf; ++ii) {
	  // Need to Populate Control Point Coordinates and Weight Vectors
	  ego		          face = fobjs[ii];
    ego            *vobjs, *eobjs;
    int             offsetCoord, offsetWeight;
    PetscInt        Nv, Ne, wRowStart = 0;
	  PetscHashIter   hashKeyIter, wHashKeyIter;
	  PetscBool   	  hashKeyFound, wHashKeyFound;
  
	  id   = EG_indexBodyTopo(body, face);
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	  PetscCall(EG_getInfo(fgeom, &boclass, &bmtype, &bRef, &bPrev, &bNext));
    PetscCall(EG_getBodyTopos(body, face, NODE,  &Nv, &vobjs));
		
	  // Store Face ID to 1st Row of Control Point Vector
	  PetscCall(PetscHMapIFind(faceCntrlPtRow_Start, id, &hashKeyIter, &hashKeyFound));
		
	  if (!hashKeyFound)  {PetscCall(PetscHMapISet(faceCntrlPtRow_Start, id, cntr));}
		
	  offsetCoord = bpinfo[3] + bpinfo[6];
	  for (int jj = 0; jj < 3 * bpinfo[2] * bpinfo[5]; ++jj) {
	    cntrlPtCoords[cntr] = bprv[offsetCoord + jj];
	    cntr += 1;
	  }
		
	  // Store Face ID to 1st Row of Control Point Weight Vector
	  PetscCall(PetscHMapIFind(faceCPWeightsRow_Start, id, &wHashKeyIter, &wHashKeyFound));
		
	  if (!wHashKeyFound)  {
	    PetscCall(PetscHMapISet(faceCPWeightsRow_Start, id, wcntr));
      wRowStart = wcntr;
	  }
    
	  offsetWeight = bpinfo[3] + bpinfo[6] + (3 * bpinfo[2] * bpinfo[5]);
	  for (int jj = 0; jj < bpinfo[2] * bpinfo[5]; ++jj) {
	    cntrlPtWeights[wcntr] = bprv[offsetWeight + jj];
      cp_face[wcntr] = id;
      w_face[wcntr] = id;
      wcntr += 1;
	  }
    
    // Associate Control Points with Vertix IDs
    PetscScalar xcp, ycp, zcp;
    offsetCoord = bpinfo[3] + bpinfo[6];
    for (int jj = 0; jj < 3 * bpinfo[2] * bpinfo[5]; jj+=3) {
      xcp = bprv[offsetCoord + jj + 0];
      ycp = bprv[offsetCoord + jj + 1];
      zcp = bprv[offsetCoord + jj + 2];
      
      //Initialize Control Point and Weight to Vertex ID relationship to -1
      cp_vertex[vcntr] = -1;
      w_vertex[vcntr] = -1;
      cp_edge[vcntr] = -1;
      w_edge[vcntr] = -1;
      
      for (int kk = 0; kk < Nv; ++kk) {
        int       vid;
        double    vCoords[3];
        PetscScalar vDelta;
        ego vertex = vobjs[kk];
        vid   = EG_indexBodyTopo(body, vertex);
        PetscCall(EG_evaluate(vertex, NULL, vCoords));
        vDelta = sqrt(pow(vCoords[0] - xcp, 2.0) + pow(vCoords[1] - ycp, 2.0) + pow(vCoords[2] - zcp, 2.0));

        if (vDelta < 1.0E-15) {
          cp_vertex[vcntr] = vid;
          w_vertex[vcntr] = vid;
        }
      }
      vcntr += 1;
    }
    EG_free(vobjs);
    
    // Associate Control Points with Edge IDs
    PetscCall(EG_getBodyTopos(body, face, EDGE,  &Ne, &eobjs));
    
    int  cpV1, cpV2;
    int  minID, maxID;
    
    // Along vmin axis
    minID = wRowStart;
    maxID = wRowStart + (bpinfo[2] - 1);
    cpV1 = cp_vertex[minID];
    cpV2 = cp_vertex[maxID];
    for (int jj = 0; jj < Ne; ++jj) {
      ego   edge = eobjs[jj];
      ego   egeom, *nobjs;
      int   eoclass, emtype, Nn, *nsenses;
      int   n1ID, n2ID, eid;
      
      eid  = EG_indexBodyTopo(body, edge);
      PetscCall(EG_getTopology(edge, &egeom, &eoclass, &emtype, NULL, &Nn, &nobjs, &nsenses));
      
      if (emtype != DEGENERATE) {
        // Get IDs for current Edge's End Vertices
        n1ID   = EG_indexBodyTopo(body, nobjs[0]);
        n2ID   = EG_indexBodyTopo(body, nobjs[1]);
        
        if ((cpV1 == n1ID || cpV1 == n2ID) && (cpV2 == n1ID || cpV2 == n2ID)) {
          for (int kk = minID + 1; kk < maxID; ++kk) {            
            cp_edge[kk] = eid;
            w_edge[kk] = eid;
          }
        }
      }
    }
    
    // Along vmax axis
    minID = wRowStart + (bpinfo[2] * (bpinfo[5] - 1));
    maxID = wRowStart + (bpinfo[2] * bpinfo[5] - 1);
    
    cpV1 = cp_vertex[minID];
    cpV2 = cp_vertex[maxID];
    for (int jj = 0; jj < Ne; ++jj) {
      ego   edge = eobjs[jj];
      ego   egeom, *nobjs;
      int   eoclass, emtype, Nn, *nsenses;
      int   n1ID, n2ID, eid;
      
      eid   = EG_indexBodyTopo(body, edge);
      PetscCall(EG_getTopology(edge, &egeom, &eoclass, &emtype, NULL, &Nn, &nobjs, &nsenses));
      
      if (emtype != DEGENERATE) {
        // Get IDs for current Edge's End Vertices
        n1ID   = EG_indexBodyTopo(body, nobjs[0]);
        n2ID   = EG_indexBodyTopo(body, nobjs[1]);
        
        if ((cpV1 == n1ID || cpV1 == n2ID) && (cpV2 == n1ID || cpV2 == n2ID)) {
          for (int kk = minID + 1; kk < maxID - 1; ++kk) {            
            cp_edge[kk] = eid;
            w_edge[kk] = eid;
          }
        }
      }
    }
    
    // Along umin axis
    minID = wRowStart;
    maxID = wRowStart + (bpinfo[2] * (bpinfo[5] - 1));
    
    cpV1 = cp_vertex[minID];
    cpV2 = cp_vertex[maxID];
    for (int jj = 0; jj < Ne; ++jj) {
      ego   edge = eobjs[jj];
      ego   egeom, *nobjs;
      int   eoclass, emtype, Nn, *nsenses;
      int   n1ID, n2ID, eid;
      
      eid   = EG_indexBodyTopo(body, edge);
      PetscCall(EG_getTopology(edge, &egeom, &eoclass, &emtype, NULL, &Nn, &nobjs, &nsenses));
      
      if (emtype != DEGENERATE) {
        // Get IDs for current Edge's End Vertices
        n1ID   = EG_indexBodyTopo(body, nobjs[0]);
        n2ID   = EG_indexBodyTopo(body, nobjs[1]);
        
        if ((cpV1 == n1ID || cpV1 == n2ID) && (cpV2 == n1ID || cpV2 == n2ID)) {
          for (int kk = minID + bpinfo[2]; kk < maxID; kk+=bpinfo[2]) {            
            cp_edge[kk] = eid;
            w_edge[kk] = eid;
          }
        }
      }
    }
    
    // Along umax axis
    minID = wRowStart + (bpinfo[2] - 1);
    maxID = wRowStart + (bpinfo[2] * bpinfo[5]) - 1;
    cpV1 = cp_vertex[minID];
    cpV2 = cp_vertex[maxID];
    for (int jj = 0; jj < Ne; ++jj) {
      ego   edge = eobjs[jj];
      ego   egeom, *nobjs;
      int   eoclass, emtype, Nn, *nsenses;
      int   n1ID, n2ID, eid;
      
      eid  = EG_indexBodyTopo(body, edge);
      PetscCall(EG_getTopology(edge, &egeom, &eoclass, &emtype, NULL, &Nn, &nobjs, &nsenses));
      
      if (emtype != DEGENERATE) {
        // Get IDs for current Edge's End Vertices
        n1ID   = EG_indexBodyTopo(body, nobjs[0]);
        n2ID   = EG_indexBodyTopo(body, nobjs[1]);
        
        if ((cpV1 == n1ID || cpV1 == n2ID) && (cpV2 == n1ID || cpV2 == n2ID)) {
          for (int kk = minID + bpinfo[2]; kk < maxID; kk+=bpinfo[2]) {            
            cp_edge[kk] = eid;
            w_edge[kk] = eid;
          }
        }
      }
    }
    EG_free(eobjs);
	}
    
  // Determine Control Point Equivalance Matrix relating Control Points between Surfaces
  //     Note: The Weights will also be tied together in the same manner
  //           Also can use the Weight Hash Table for Row Start ID of each Face
	const PetscInt  cpRowSize = totalNumCPs;
	const PetscInt  cpColSize = cpRowSize;
  PetscInt       *maxNumRelatePtr;
  PetscInt        maxNumRelate = 0;
	
	// Create Point Surface Gradient Matrix
	MatCreate(PETSC_COMM_WORLD, &cpEquiv);
  MatSetSizes(cpEquiv,PETSC_DECIDE, PETSC_DECIDE, cpRowSize, cpColSize);
  MatSetType(cpEquiv, MATAIJ);
  MatSetUp(cpEquiv);
  
  for (int ii = 0; ii < totalNumCPs; ++ii) {
    PetscScalar x1, y1, z1;
    PetscInt    maxRelateTemp = 0;
    x1 = cntrlPtCoords[(3 * ii) + 0];
    y1 = cntrlPtCoords[(3 * ii) + 1];
    z1 = cntrlPtCoords[(3 * ii) + 2];
    
    for (int jj = 0; jj < totalNumCPs; ++jj) {
      PetscScalar x2, y2, z2;
      PetscScalar cpDelta, eqFactor;
      x2 = cntrlPtCoords[(3 * jj) + 0];
      y2 = cntrlPtCoords[(3 * jj) + 1];
      z2 = cntrlPtCoords[(3 * jj) + 2];
      
      cpDelta = sqrt(pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0));
      if (cpDelta < 1.0E-15) {
        eqFactor = 1.0;
        maxRelateTemp += 1;
      } else {
        eqFactor = 0.0;
      }

      // Store Results in Petsc Matrix
      PetscCall(MatSetValue(cpEquiv, ii, jj, eqFactor, INSERT_VALUES));
    }
    if (maxRelateTemp > maxNumRelate) {maxNumRelate = maxRelateTemp;}
  }
  maxNumRelatePtr = &maxNumRelate;
  
  // Assemble Point Surface Grad Matrix
	MatAssemblyBegin(cpEquiv, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(cpEquiv, MAT_FINAL_ASSEMBLY);
  	
	// Attach Control Point and Weight Data to DM
	{
	  PetscContainer cpOrgObj, cpCoordObj, cpCoordLengthObj;
	  PetscContainer wOrgObj, wValObj, wDataLengthObj;
    PetscContainer cpEquivObj, cp_faceObj, cp_edgeObj, cp_vertexObj;
    PetscContainer w_faceObj, w_edgeObj, w_vertexObj;
    PetscContainer maxNumRelateObj;
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Hash Table", (PetscObject *) &cpOrgObj));
    if (!cpOrgObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpOrgObj));
      PetscCall(PetscContainerSetPointer(cpOrgObj, faceCntrlPtRow_Start));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Hash Table", (PetscObject) cpOrgObj));
      PetscCall(PetscContainerDestroy(&cpOrgObj));
    } else {
      PetscCall(PetscContainerSetPointer(cpOrgObj, faceCntrlPtRow_Start));
    }
	  
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Coordinates", (PetscObject *) &cpCoordObj));
    if (!cpCoordObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpCoordObj));
      PetscCall(PetscContainerSetPointer(cpCoordObj, cntrlPtCoords));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Coordinates", (PetscObject) cpCoordObj));
      PetscCall(PetscContainerDestroy(&cpCoordObj));
    } else {
      PetscCall(PetscContainerSetPointer(cpCoordObj, cntrlPtCoords));
    }
		
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Coordinate Data Length", (PetscObject *) &cpCoordLengthObj));
    if (!cpCoordLengthObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpCoordLengthObj));
      PetscCall(PetscContainerSetPointer(cpCoordLengthObj, cpCoordDataLengthPtr));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Coordinate Data Length", (PetscObject) cpCoordLengthObj));
      PetscCall(PetscContainerDestroy(&cpCoordLengthObj));
    } else {
      PetscCall(PetscContainerSetPointer(cpCoordLengthObj, cpCoordDataLengthPtr));
    }

    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weights Hash Table", (PetscObject *) &wOrgObj));
    if (!wOrgObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wOrgObj));
      PetscCall(PetscContainerSetPointer(wOrgObj, faceCPWeightsRow_Start));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weights Hash Table", (PetscObject) wOrgObj));
      PetscCall(PetscContainerDestroy(&wOrgObj));
    } else {
      PetscCall(PetscContainerSetPointer(wOrgObj, faceCPWeightsRow_Start));
    }

    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weight Data", (PetscObject *) &wValObj));
    if (!wValObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wValObj));
      PetscCall(PetscContainerSetPointer(wValObj, cntrlPtWeights));
      //PetscCall(PetscContainerSetUserDestroy(wValObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight Data", (PetscObject) wValObj));
      PetscCall(PetscContainerDestroy(&wValObj));
    } else {
      PetscCall(PetscContainerSetPointer(wValObj, cntrlPtWeights));
    }

    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weight Data Length", (PetscObject *) &wDataLengthObj));
    if (!wDataLengthObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &wDataLengthObj));
      PetscCall(PetscContainerSetPointer(wDataLengthObj, wDataLengthPtr));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight Data Length", (PetscObject) wDataLengthObj));
      PetscCall(PetscContainerDestroy(&wDataLengthObj));	
    } else {
      PetscCall(PetscContainerSetPointer(wDataLengthObj, wDataLengthPtr));
    }

    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Equivalancy Matrix", (PetscObject *) &cpEquivObj));
    if (!cpEquivObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cpEquivObj));
      PetscCall(PetscContainerSetPointer(cpEquivObj, cpEquiv));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Equivalancy Matrix", (PetscObject) cpEquivObj));
      PetscCall(PetscContainerDestroy(&cpEquivObj));
    } else {
      PetscCall(PetscContainerSetPointer(cpEquivObj, cpEquiv));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Maximum Number Control Point Equivalency", (PetscObject *) &maxNumRelateObj));
    if (!maxNumRelateObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &maxNumRelateObj));
      PetscCall(PetscContainerSetPointer(maxNumRelateObj, maxNumRelatePtr));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Maximum Number Control Point Equivalency", (PetscObject) maxNumRelateObj));
      PetscCall(PetscContainerDestroy(&maxNumRelateObj));
    } else {
      PetscCall(PetscContainerSetPointer(maxNumRelateObj, maxNumRelatePtr));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point - Face Map", (PetscObject *) &cp_faceObj));
    if (!cp_faceObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cp_faceObj));
      PetscCall(PetscContainerSetPointer(cp_faceObj, cp_face));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point - Face Map", (PetscObject) cp_faceObj));
      PetscCall(PetscContainerDestroy(&cp_faceObj));
    } else {
      PetscCall(PetscContainerSetPointer(cp_faceObj, cp_face));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weight - Face Map", (PetscObject *) &w_faceObj));
    if (!w_faceObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &w_faceObj));
      PetscCall(PetscContainerSetPointer(w_faceObj, w_face));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight - Face Map", (PetscObject) w_faceObj));
      PetscCall(PetscContainerDestroy(&w_faceObj));
    } else {
      PetscCall(PetscContainerSetPointer(w_faceObj, w_face));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point - Edge Map", (PetscObject *) &cp_edgeObj));
    if (!cp_edgeObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cp_edgeObj));
      PetscCall(PetscContainerSetPointer(cp_edgeObj, cp_edge));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point - Edge Map", (PetscObject) cp_edgeObj));
      PetscCall(PetscContainerDestroy(&cp_edgeObj));
    } else {
      PetscCall(PetscContainerSetPointer(cp_edgeObj, cp_edge));
    }

    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weight - Edge Map", (PetscObject *) &w_edgeObj));
    if (!w_edgeObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &w_edgeObj));
      PetscCall(PetscContainerSetPointer(w_edgeObj, w_edge));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight - Edge Map", (PetscObject) w_edgeObj));
      PetscCall(PetscContainerDestroy(&w_edgeObj));
    } else {
      PetscCall(PetscContainerSetPointer(w_edgeObj, w_edge));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point - Vertex Map", (PetscObject *) &cp_vertexObj));
    if (!cp_vertexObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cp_vertexObj));
      PetscCall(PetscContainerSetPointer(cp_vertexObj, cp_vertex));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point - Vertex Map", (PetscObject) cp_vertexObj));
      PetscCall(PetscContainerDestroy(&cp_vertexObj));
    } else {
      PetscCall(PetscContainerSetPointer(cp_vertexObj, cp_vertex));
    }
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weight - Vertex Map", (PetscObject *) &w_vertexObj));
    if (!w_vertexObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &w_vertexObj));
      PetscCall(PetscContainerSetPointer(w_vertexObj, w_vertex));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Control Point Weight - Vertex Map", (PetscObject) w_vertexObj));
      PetscCall(PetscContainerDestroy(&w_vertexObj));
    } else {
      PetscCall(PetscContainerSetPointer(cp_vertexObj, w_vertexObj));
    } 
	}
	
	// Define Matrix to store  Geometry Gradient information dGeom_i/dCPj_i
	PetscInt		    gcntr = 0;
	const PetscInt  rowSize = 3 * maxNumCPs * totalNumPoints;
	const PetscInt  colSize = 4 * Nf;
	
	// Create Point Surface Gradient Matrix
	MatCreate(PETSC_COMM_WORLD, &pointSurfGrad);
  MatSetSizes(pointSurfGrad,PETSC_DECIDE, PETSC_DECIDE, rowSize, colSize);
  MatSetType(pointSurfGrad, MATAIJ);
  MatSetUp(pointSurfGrad);
	
	// Create Hash Table to store Point's stare row in surfaceGrad[][]
	PetscCall(PetscHMapICreate(&pointSurfGradRow_Start));
	
	// Get Coordinates for the DMPlex point
	DM 			       cdm;
	PetscInt	     dE, Nv;
	Vec            coordinatesLocal;
	PetscScalar   *coords = NULL;
  
	PetscCall(DMGetCoordinateDM(dm, &cdm));
	PetscCall(DMGetCoordinateDim(dm, &dE));
	PetscCall(DMGetCoordinatesLocal(dm, &coordinatesLocal));
	
	// CYCLE THROUGH FACEs
  PetscScalar  maxGrad = 0.;
	for (int ii = 0; ii < Nf; ++ii) {
	  ego		          face = fobjs[ii];
	  ego            *eobjs, *nobjs;
	  PetscInt	      fid, Ne, Nn;	
	  DMLabel		      faceLabel, edgeLabel, nodeLabel;
	  PetscHMapI      currFaceUniquePoints = NULL;
	  IS			        facePoints, edgePoints, nodePoints;
	  const PetscInt *fIndices, *eIndices, *nIndices;
	  PetscInt        fSize, eSize, nSize;
	  PetscHashIter   fHashKeyIter, eHashKeyIter, nHashKeyIter, pHashKeyIter;
	  PetscBool   	  fHashKeyFound, eHashKeyFound, nHashKeyFound, pHashKeyFound;
	  PetscInt	      cfCntr = 0;
		
	  // Get Geometry Object for the Current FACE
	  PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
	  PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
	
	  // Get all EDGE and NODE objects attached to the current FACE
	  PetscCall(EG_getBodyTopos(body, face, EDGE,  &Ne, &eobjs));
	  PetscCall(EG_getBodyTopos(body, face, NODE,  &Nn, &nobjs));
		
	  // Get all DMPlex Points that have DMLabel "EGADS Face ID" and store them in a Hash Table for later use
	  fid   = EG_indexBodyTopo(body, face);
	  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
	  PetscCall(DMLabelGetStratumIS(faceLabel, fid, &facePoints));
	  PetscCall(ISGetIndices(facePoints, &fIndices));
	  PetscCall(ISGetSize(facePoints, &fSize));
		
	  PetscCall(PetscHMapICreate(&currFaceUniquePoints));
	  
	  for (int jj = 0; jj < fSize; ++jj) {
	    PetscCall(PetscHMapIFind(currFaceUniquePoints, fIndices[jj], &fHashKeyIter, &fHashKeyFound));
			
      if (!fHashKeyFound) {
        PetscCall(PetscHMapISet(currFaceUniquePoints, fIndices[jj], cfCntr));
	      cfCntr += 1;
      }
			
      PetscCall(PetscHMapIFind(pointSurfGradRow_Start, fIndices[jj], &pHashKeyIter, &pHashKeyFound));
			
      if (!pHashKeyFound) {
        PetscCall(PetscHMapISet(pointSurfGradRow_Start, fIndices[jj], gcntr));
        gcntr += 3 * maxNumCPs;
      }
	  }
	  PetscCall(ISRestoreIndices(facePoints, &fIndices));
	  PetscCall(ISDestroy(&facePoints));
		
	  // Get all DMPlex Points that have DMLable "EGADS Edge ID" attached to the current FACE and store them in a Hash Table for later use.
	  for (int jj = 0; jj < Ne; ++jj) {
      ego			edge = eobjs[jj];
		  PetscBool	containLabelValue;
			
      id = EG_indexBodyTopo(body, edge);
      PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
      PetscCall(DMLabelHasValue(edgeLabel, id, &containLabelValue));
			
      if (containLabelValue) {
        PetscCall(DMLabelGetStratumIS(edgeLabel, id, &edgePoints));
        PetscCall(ISGetIndices(edgePoints, &eIndices));
        PetscCall(ISGetSize(edgePoints, &eSize));
			
        for (int kk = 0; kk < eSize; ++kk) {
          PetscCall(PetscHMapIFind(currFaceUniquePoints, eIndices[kk], &eHashKeyIter, &eHashKeyFound));
			
          if (!eHashKeyFound) {
            PetscCall(PetscHMapISet(currFaceUniquePoints, eIndices[kk], cfCntr));
            cfCntr += 1;
          }
					
          PetscCall(PetscHMapIFind(pointSurfGradRow_Start, eIndices[kk], &pHashKeyIter, &pHashKeyFound));
			
          if (!pHashKeyFound) {
            PetscCall(PetscHMapISet(pointSurfGradRow_Start, eIndices[kk], gcntr));
            gcntr += 3 * maxNumCPs;
          }
        }
        PetscCall(ISRestoreIndices(edgePoints, &eIndices));
        PetscCall(ISDestroy(&edgePoints));
      }
	  }
		
	  // Get all DMPlex Points that have DMLabel "EGADS Vertex ID" attached to the current FACE and store them in a Hash Table for later use.
	  for (int jj = 0; jj < Nn; ++jj) {
	    ego		node = nobjs[jj];
		  id = EG_indexBodyTopo(body, node);
		  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &nodeLabel));
		  PetscCall(DMLabelGetStratumIS(nodeLabel, id, &nodePoints));
		  PetscCall(ISGetIndices(nodePoints, &nIndices));
		  PetscCall(ISGetSize(nodePoints, &nSize));
			
		  for (int kk = 0; kk < nSize; ++kk) {
		    PetscCall(PetscHMapIFind(currFaceUniquePoints, nIndices[kk], &nHashKeyIter, &nHashKeyFound));
			
		    if (!nHashKeyFound) {
		      PetscCall(PetscHMapISet(currFaceUniquePoints, nIndices[kk], cfCntr));
			    cfCntr += 1;
		    }
				
		    PetscCall(PetscHMapIFind(pointSurfGradRow_Start, nIndices[kk], &pHashKeyIter, &pHashKeyFound));	
		    if (!pHashKeyFound) {
			    PetscCall(PetscHMapISet(pointSurfGradRow_Start, nIndices[kk], gcntr));
			    gcntr += 3 * maxNumCPs;
		    }
		  }
		  PetscCall(ISRestoreIndices(nodePoints, &nIndices));
		  PetscCall(ISDestroy(&nodePoints));
	  }
	
	  // Get the Total Number of entries in the Hash Table
	  PetscInt	currFaceUPSize;
	  PetscCall(PetscHMapIGetSize(currFaceUniquePoints, &currFaceUPSize));
		
	  // Get Keys
	  PetscInt	currFaceUPKeys[currFaceUPSize], off = 0;
	  PetscCall(PetscHMapIGetKeys(currFaceUniquePoints, &off, currFaceUPKeys));
		
    // Get Current Face Surface Area
    PetscScalar  fSA, faceData[14];
    PetscCall(EG_getMassProperties(face, faceData));
    fSA = faceData[1];
    
    // Get Start Row in cpEquiv Matrix
    PetscHashIter   Witer;
    PetscBool       Wfound;
    PetscInt  faceWStartRow;
    PetscCall(PetscHMapIFind(faceCPWeightsRow_Start, fid, &Witer, &Wfound));
    if (!Wfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Weights Hash Table");
    PetscCall(PetscHMapIGet(faceCPWeightsRow_Start, fid, &faceWStartRow));

	  // Cycle through all points on the current FACE
	  for (int jj = 0; jj < currFaceUPSize; ++jj) {
	    PetscInt	currPointID = currFaceUPKeys[jj];
      PetscCall(DMPlexVecGetClosure(cdm, NULL, coordinatesLocal, currPointID, &Nv, &coords));

      // Get UV position of FACE
      double  params[2], range[4], eval[18];
      int			peri;
      PetscCall(EG_getRange(face, range, &peri));
      PetscCall(DMPlex_EGADS_FACE_XYZtoUV_Internal(coords, face, range, 0, dE, params));	
      PetscCall(EG_evaluate(face, params, eval));	
  
      // Make a new SURFACE Geometry by changing the location of the Control Points
      int    prvSize = bpinfo[3] + bpinfo[6] + (4 * bpinfo[2]*bpinfo[5]);
      double nbprv[prvSize];
			
      // Cycle through each Control Point
      double  denomNew, denomOld;
      double  deltaCoord = 1.0E-4;
      int     offset = bpinfo[3] + bpinfo[6];
      int     wOffset = offset + (3 * bpinfo[2] * bpinfo[5]);
      for (int ii = 0; ii < bpinfo[2]*bpinfo[5]; ++ii){
        // Cycle through each direction (x, then y, then z)
        if (jj == 0) {
          // Get the Number Control Points that are the same as the current points
          //    We are looking for repeated Control Points
          PetscInt    commonCPcntr = 0;
          for (int mm = 0; mm < bpinfo[2]*bpinfo[5]; ++mm) {
            PetscScalar   matValue;
            PetscCall(MatGetValue(cpEquiv, faceWStartRow + ii, faceWStartRow + mm, &matValue));
          
            if (matValue > 0.0) {
              commonCPcntr += 1;
            }
          }
        }
        
        for (int kk = 0; kk < 4; ++kk) {
          // Reinitialize nbprv[] values because we only want to change one value at a time
          for (int mm = 0; mm < prvSize; ++mm) {
            nbprv[mm] = bprv[mm];
          }          
          
          if (kk == 0) {	      //X
            nbprv[offset + 0] = bprv[offset + 0] + deltaCoord;
            nbprv[offset + 1] = bprv[offset + 1];
            nbprv[offset + 2] = bprv[offset + 2];
            denomNew = nbprv[offset + 0];
            denomOld = bprv[offset + 0];
          } else if (kk == 1) {	//Y
            nbprv[offset + 0] = bprv[offset + 0];
            nbprv[offset + 1] = bprv[offset + 1] + deltaCoord;
            nbprv[offset + 2] = bprv[offset + 2];
            denomNew = nbprv[offset + 1];
            denomOld = bprv[offset + 1];
          } else if (kk == 2) {	//Z
            nbprv[offset + 0] = bprv[offset + 0];
            nbprv[offset + 1] = bprv[offset + 1];
            nbprv[offset + 2] = bprv[offset + 2] + deltaCoord;
            denomNew = nbprv[offset + 2];
            denomOld = bprv[offset + 2];
          } else if (kk == 3) {	// Weights
            nbprv[wOffset + ii] = bprv[wOffset + ii] + deltaCoord;
            denomNew = nbprv[wOffset + ii];
            denomOld = bprv[wOffset + ii];
          } else {
            // currently do nothing
          }	
                   
          // Create New Surface Based on New Control Points or Weights
          ego newgeom, context;
          PetscCall(EG_getContext(face, &context));
          PetscCall(EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bpinfo, nbprv, &newgeom));
					
          // Evaluate new (x, y, z) Point Position based on new Surface Definition
          double newCoords[18];
          PetscCall(EG_getRange(newgeom, range, &peri));
          PetscCall(DMPlex_EGADS_FACE_XYZtoUV_Internal(coords, newgeom, range, 0, dE, params));
          PetscCall(EG_evaluate(newgeom, params, newCoords));
          
          // Calculate Surface Area Gradients wrt Control Points and Weights using the local discrete FACE only
          //      NOTE 1: Will not provide Volume Gradient wrt to Control Points and Weights.
          //      NOTE 2: This is faster than below where an entire new solid geometry is created for each
          //              Control Point and Weight gradient
          if (!fullGeomGrad) {
            // Create new FACE based on new SURFACE geometry
            if (jj == 0) {    // only for 1st DMPlex Point because we only per CP or Weight
              double newFaceRange[4];
              int newFacePeri;
              PetscCall(EG_getRange(newgeom, newFaceRange, &newFacePeri));
            
              ego newface;
              PetscCall(EG_makeFace(newgeom, SFORWARD, newFaceRange, &newface));
            
              // Get New Face Surface Area
              PetscScalar  newfSA, newFaceData[14];
              PetscCall(EG_getMassProperties(newface, newFaceData));
              newfSA = newFaceData[1];
            
              // Update Control Points
              PetscHashIter   CPiter, Witer;
              PetscBool       CPfound, Wfound;
              PetscInt        faceCPStartRow, faceWStartRow;
              
              PetscScalar  dSAdCPi;
              dSAdCPi = (newfSA - fSA) / (denomNew - denomOld);
              
              if (kk < 3) {
                PetscCall(PetscHMapIFind(faceCntrlPtRow_Start, fid, &CPiter, &CPfound));
                if (!CPfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Hash Table");
                PetscCall(PetscHMapIGet(faceCntrlPtRow_Start, fid, &faceCPStartRow));

                gradSACP[faceCPStartRow + (ii * 3) + kk] = dSAdCPi;

                if (fabs(dSAdCPi) > maxGrad) maxGrad = fabs(dSAdCPi);

              } else if (kk == 3) {
                PetscCall(PetscHMapIFind(faceCPWeightsRow_Start, fid, &Witer, &Wfound));
                if (!Wfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Hash Table");
                PetscCall(PetscHMapIGet(faceCPWeightsRow_Start, fid, &faceWStartRow));
                
                gradSAW[faceWStartRow + ii] = dSAdCPi;

              } else {
                // Do Nothing
              }
            }
          }

          // Now Calculate the Surface Gradient for the change in x-component Control Point
          PetscScalar dxdCx = (newCoords[0] - coords[0]) / deltaCoord;
          PetscScalar dxdCy = (newCoords[1] - coords[1]) / deltaCoord;
          PetscScalar dxdCz = (newCoords[2] - coords[2]) / deltaCoord;
					
          // Store Gradient Information in surfaceGrad[][] Matrix
          PetscInt   startRow;
          PetscCall(PetscHMapIGet(pointSurfGradRow_Start, currPointID, &startRow));
	
          // Store Results in Petsc Matrix
          PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 0, ((fid - 1) * 4) + kk, dxdCx, INSERT_VALUES));
          PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 1, ((fid - 1) * 4) + kk, dxdCy, INSERT_VALUES));
          PetscCall(MatSetValue(pointSurfGrad, startRow + (ii * 3) + 2, ((fid - 1) * 4) + kk, dxdCz, INSERT_VALUES));
        }
        offset += 3;
      }
      PetscCall(DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, currPointID, &Nv, &coords));
	  }
	}
  
	// Assemble Point Surface Grad Matrix
	MatAssemblyBegin(pointSurfGrad, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(pointSurfGrad, MAT_FINAL_ASSEMBLY);
  
  if (fullGeomGrad) {
    // Calculate Surface Area and Volume Control Point and Control Point Weight Gradients
    //    Note: This is much slower than above due to a new solid geometry being created for
    //          each change in Control Point and Control Point Weight. However, this method
    //          will provide the Volume Gradient.
    
    // Get Current Face Surface Area
    PetscScalar  bodyVol, bodySA, bodyData[14];
    PetscCall(EG_getMassProperties(body, bodyData));
    bodyVol = bodyData[0];
    bodySA = bodyData[1];
    
    // Cycle through Control Points
    for (int ii = 0; ii < totalNumCPs; ++ii) {    // ii should also be the row in cpEquiv for the Control Point
      // Cycle through X, Y, Z, W changes
      for (int jj = 0; jj < 4; ++jj) {
        // Cycle Through Faces
        double  denomNew = 0.0, denomOld = 0.0;
        double  deltaCoord = 1.0E-4;
        ego     newFaces[Nf];
        for (int kk = 0; kk < Nf; ++kk) {
          ego         face;
          PetscInt    currFID = kk + 1;
          PetscCall(EG_objectBodyTopo(body, FACE, currFID, &face));    // Get Current FACE
      
          // Get Geometry Object for the Current FACE
          PetscCall(EG_getTopology(face, &fgeom, &foclass, &fmtype, fdata, &Nl, &lobjs, &lsenses));
          PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
          
          // Make a new SURFACE Geometry by changing the location of the Control Points
          int    prvSize = bpinfo[3] + bpinfo[6] + (4 * bpinfo[2]*bpinfo[5]);
          double nbprv[prvSize];
        
          // Reinitialize nbprv[] values because we only want to change one value at a time
          for (int mm = 0; mm < prvSize; ++mm) {
            nbprv[mm] = bprv[mm];
          }
          
          // Get Control Point Row and Column Start for cpEquiv
          PetscHashIter   Witer;
          PetscBool       Wfound;
          PetscInt        faceWStartRow;
          PetscCall(PetscHMapIFind(faceCPWeightsRow_Start, currFID, &Witer, &Wfound));
          if (!Wfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Weights Hash Table");
          PetscCall(PetscHMapIGet(faceCPWeightsRow_Start, currFID, &faceWStartRow));
          
          // Modify the Current Control Point on this FACE and All Other FACES
          // IMPORTANT!!! If you do not move all identical Control Points on other FACES
          //              you will not generate a solid body. You will generate a set of
          //              disconnected surfaces that have gap(s) between them.
          int     offset = bpinfo[3] + bpinfo[6];
          int     wOffset = offset + (3 * bpinfo[2] * bpinfo[5]);
          for (int mm = 0; mm < bpinfo[2]*bpinfo[5]; ++mm) {
            PetscScalar   matValue;
            PetscCall(MatGetValue(cpEquiv, ii, faceWStartRow + mm, &matValue));
              
            if (matValue > 0.0) {
              if (jj == 0) {	      //X
                nbprv[offset + (3 * mm) + 0] = bprv[offset + (3 * mm) + 0] + deltaCoord;
                nbprv[offset + (3 * mm) + 1] = bprv[offset + (3 * mm) + 1];
                nbprv[offset + (3 * mm) + 2] = bprv[offset + (3 * mm) + 2];
                denomNew = nbprv[offset + (3 * mm) + 0];
                denomOld = bprv[offset + (3 * mm) + 0];
              } else if (jj == 1) {	//Y
                nbprv[offset + (3 * mm) + 0] = bprv[offset + (3 * mm) + 0];
                nbprv[offset + (3 * mm) + 1] = bprv[offset + (3 * mm) + 1] + deltaCoord;
                nbprv[offset + (3 * mm) + 2] = bprv[offset + (3 * mm) + 2];
                denomNew = nbprv[offset + (3 * mm) + 1];
                denomOld = bprv[offset + (3 * mm) + 1];
              } else if (jj == 2) {	//Z
                nbprv[offset + (3 * mm) + 0] = bprv[offset + (3 * mm) + 0];
                nbprv[offset + (3 * mm) + 1] = bprv[offset + (3 * mm) + 1];
                nbprv[offset + (3 * mm) + 2] = bprv[offset + (3 * mm) + 2] + deltaCoord;
                denomNew = nbprv[offset + (3 * mm) + 2];
                denomOld = bprv[offset + (3 * mm) + 2];
              } else if (jj == 3) {	// Weights
                nbprv[wOffset + mm] = bprv[wOffset + mm] + deltaCoord;
                denomNew = nbprv[wOffset + mm];
                denomOld = bprv[wOffset + mm];
              } else {
                // currently do nothing
              }
            }
          }
          
          // Create New Surface Based on New Control Points or Weights
          ego newgeom, context;
          PetscCall(EG_getContext(face, &context));
          PetscCall(EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bpinfo, nbprv, &newgeom));
          
          // Create New FACE based on modified geometry
          double newFaceRange[4];
          int    newFacePeri;
          PetscCall(EG_getRange(newgeom, newFaceRange, &newFacePeri));
            
          ego newface;
          PetscCall(EG_makeFace(newgeom, SFORWARD, newFaceRange, &newface));
          
          // store new face for later assembly
          newFaces[kk] = newface;
        }
        
        // X-WANT TO BUILD THE NEW GEOMETRY, X-GET NEW SA AND PERFORM dSA/dCPi CALCS HERE <---
        // Sew New Faces together to get a new model
        ego newmodel;
        PetscCall(EG_sewFaces(Nf, newFaces, 0.0, 0, &newmodel));
   
        // Get Surface Area and Volume of New/Updated Solid Body
        PetscScalar newData[14];
        PetscCall(EG_getTopology(newmodel, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
        ego nbody = bodies[0];
        PetscCall(EG_getMassProperties(nbody, newData));
      
        PetscScalar  dSAdCPi, dVdCPi;
        PetscScalar  nbodyVol = newData[0], nbodySA = newData[1];
        
        // Calculate Gradients wrt to Control Points and Control Points Weights depending on jj value
        dSAdCPi = (nbodySA - bodySA) / (denomNew - denomOld);
        dVdCPi = (nbodyVol - bodyVol) / (denomNew - denomOld);
              
        if (jj < 3) {
          // Gradienst wrt to Control Points
          gradSACP[(ii * 3) + jj] = dSAdCPi;
          gradVCP[(ii * 3) + jj] = dVdCPi;
        } else if (jj == 3) {
          // Gradients wrt to Control Point Weights
          gradSAW[ii] = dSAdCPi;
          gradVW[ii] = dVdCPi;
        } else {
          // Do Nothing
        } 
      }
    }
  }  
	
	// Attach Surface Gradient Hash Table and Matrix to DM
	{
	  PetscContainer surfGradOrgObj, surfGradObj;
    PetscContainer gradSACPObj, gradSAWObj;
    PetscContainer gradVCPObj, gradVWObj;

    PetscCall(PetscObjectQuery((PetscObject) dm, "Surface Gradient Hash Table", (PetscObject *) &surfGradOrgObj));
    if (!surfGradOrgObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &surfGradOrgObj));
      PetscCall(PetscContainerSetPointer(surfGradOrgObj, pointSurfGradRow_Start));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Gradient Hash Table", (PetscObject) surfGradOrgObj));
      PetscCall(PetscContainerDestroy(&surfGradOrgObj));
    } else {
      PetscCall(PetscContainerSetPointer(surfGradOrgObj, pointSurfGradRow_Start));
    } 

    PetscCall(PetscObjectQuery((PetscObject) dm, "Surface Gradient Matrix", (PetscObject *) &surfGradObj));
    if (!surfGradObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &surfGradObj));
      PetscCall(PetscContainerSetPointer(surfGradObj, pointSurfGrad));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Gradient Matrix", (PetscObject) surfGradObj));
      PetscCall(PetscContainerDestroy(&surfGradObj));
    } else {
      PetscCall(PetscContainerSetPointer(surfGradObj, pointSurfGrad));
    } 
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Surface Area Control Point Gradient", (PetscObject *) &gradSACPObj));
    if (!gradSACPObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &gradSACPObj));
      PetscCall(PetscContainerSetPointer(gradSACPObj, gradSACP));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Area Control Point Gradient", (PetscObject) gradSACPObj));
      PetscCall(PetscContainerDestroy(&gradSACPObj));
    } else {
      PetscCall(PetscContainerSetPointer(gradSACPObj, gradSACP));
    } 
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Surface Area Weights Gradient", (PetscObject *) &gradSAWObj));
    if (!gradSAWObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &gradSAWObj));
      PetscCall(PetscContainerSetPointer(gradSAWObj, gradSAW));
      //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Surface Area Weights Gradient", (PetscObject) gradSAWObj));
      PetscCall(PetscContainerDestroy(&gradSAWObj));
    } else {
      PetscCall(PetscContainerSetPointer(gradSAWObj, gradSAW));
    } 

    if (fullGeomGrad) {
      PetscCall(PetscObjectQuery((PetscObject) dm, "Volume Control Point Gradient", (PetscObject *) &gradVCPObj));
      if (!gradVCPObj) {
        PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &gradVCPObj));
        PetscCall(PetscContainerSetPointer(gradVCPObj, gradVCP));
        //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
        PetscCall(PetscObjectCompose((PetscObject) dm, "Volume Control Point Gradient", (PetscObject) gradVCPObj));
        PetscCall(PetscContainerDestroy(&gradVCPObj));
      } else {
        PetscCall(PetscContainerSetPointer(gradVCPObj, gradVCP));
      } 
      
      PetscCall(PetscObjectQuery((PetscObject) dm, "Volume Weights Gradient", (PetscObject *) &gradVWObj));
      if (!gradVWObj) {
        PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &gradVWObj));
        PetscCall(PetscContainerSetPointer(gradVWObj, gradVW));
        //PetscCall(PetscContainerSetUserDestroy(cpCoordObj, DMPlexEGADSDestroy_Private));
        PetscCall(PetscObjectCompose((PetscObject) dm, "Volume Weights Gradient", (PetscObject) gradVWObj));
        PetscCall(PetscContainerDestroy(&gradVWObj));
      } else {
        PetscCall(PetscContainerSetPointer(gradVWObj, gradVW));
      } 
    }
	}
	EG_free(fobjs);
	PetscFunctionReturn(0);
}


/*@C
  DMPlexModifyEGADSGeomModel - Generates a new EGADS geometry model based in user provided Control Points
                               and Control Points Weights. Optionally, the function will inflate the DM 
                               to the new geometry and save the new geometry to a file.

  Collective

  Input Parameters:
. dm            - The DM object representing the mesh with PetscContainer containing an EGADS geometry model
. comm          - MPI_Comm object
. newCP[]       - C Array of [x, y, z] New/Updated Control Point Coordinates defining the geometry (See DMPlexGeomDataAndGrads() for format)
. newW[]        - C Array of New/Updated Control Point Weights associated with the Control Points defining the new geometry (See DMPlexGemGrads() for format)
. autoInflate   - PetscBool Flag denoting if the user would like to inflate the DM points to the new geometry.
. saveGeom      - PetscBool Flag denoting if the user would iike to save the new geometry to a file.
. stpName       - Char Array indicating the name of the file to save the new geometry to. Extension must be included and will denote type of file written.
                      *.stp or *.step = STEP File
                      *.igs or *.iges = IGES File
                              *.egads = EGADS File
                          *.egadslite = EGADSlite File
                               *.brep = BRep File (OpenCASCADE File)


  Output Parameter:
. dm       - The updated DM object representing the mesh with PetscContainers containing the updated/modified geometry

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSLiteFromFile(), DMPlexGeomDataAndGrads()
@*/
PetscErrorCode DMPlexModifyEGADSGeomModel(DM dm, MPI_Comm comm, PetscScalar newCP[], PetscScalar newW[], PetscBool autoInflate, PetscBool saveGeom, char *stpName) 
{ 
  /* EGADS/EGADSlite variables */
  ego            context, model, geom, *bodies, *lobjs, *fobjs;
  int            oclass, mtype, *senses, *lsenses;
  int            Nb, Nf, Nl, id;
  /* PETSc variables */
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  PetscContainer modelObj, cpHashTableObj, wHashTableObj;
  PetscHMapI     cpHashTable = NULL, wHashTable = NULL;
  //PetscErrorCode ierr;

  PetscFunctionBegin;
  // Look to see if DM has a Container with either a EGADS or EGADS Model
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
  if (!modelObj) {
    PetscCall(PetscObjectQuery((PetscObject) dm, "EGADSlite Model", (PetscObject *) &modelObj));
  }
  if (!modelObj) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DM does not have a EGADS Geometry Model attached to it!");

  // Get attached EGADS model (pointer)
  PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
  
  // Look to see if DM has Container for Geometry Control Point Data
  PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Hash Table", (PetscObject *) &cpHashTableObj));  
  PetscCall(PetscObjectQuery((PetscObject) dm, "Control Point Weights Hash Table", (PetscObject *) &wHashTableObj));
  
  if (!cpHashTableObj || !wHashTableObj) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "DM does not have required Geometry Data attached! Please run DMPlexGeomDataAndGrads() Function first.");
  
  // Get attached EGADS model Control Point and Weights Hash Tables and Data Arrays (pointer)
  PetscCall(PetscContainerGetPointer(cpHashTableObj, (void **) &cpHashTable));
  PetscCall(PetscContainerGetPointer(wHashTableObj, (void **) &wHashTable));
  
  // Get the number of bodies and body objects in the model
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));
  
  // Get all Faces on the body
  ego body = bodies[0];
  PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
  ego newFaces[Nf];
 
  // Update Control Point and Weight definitions for each surface
  for (int jj = 0; jj < Nf; ++jj) {
    ego     face = fobjs[jj];
    ego     bRef, bPrev, bNext;
    ego     fgeom;
    int     offset;
    int     boclass, bmtype, *bpinfo;
    double *bprv;

    // Get FACE ID and other Geometry Data
    id   = EG_indexBodyTopo(body, face);
    PetscCall(EG_getTopology(face, &fgeom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses));
    PetscCall(EG_getGeometry(fgeom, &boclass, &bmtype, &bRef, &bpinfo, &bprv));
    PetscCall(EG_getInfo(fgeom, &boclass, &bmtype, &bRef, &bPrev, &bNext));

    // Update Control Points
    PetscHashIter   CPiter, Witer;
    PetscBool       CPfound, Wfound;
    PetscInt        faceCPStartRow, faceWStartRow;
	
    PetscCall(PetscHMapIFind(cpHashTable, id, &CPiter, &CPfound));
    if (!CPfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Hash Table");
    PetscCall(PetscHMapIGet(cpHashTable, id, &faceCPStartRow));
	
    PetscCall(PetscHMapIFind(wHashTable, id, &Witer, &Wfound));
    if (!Wfound) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "FACE ID not found in Control Point Weights Hash Table");
    PetscCall(PetscHMapIGet(wHashTable, id, &faceWStartRow));
	  
    // UPDATE CONTROL POINTS Locations
    offset = bpinfo[3] + bpinfo[6];
    for (int ii = 0; ii < 3 * bpinfo[2] * bpinfo[5]; ++ii){
      bprv[offset + ii] = newCP[faceCPStartRow + ii];
    }
	
    // UPDATE CONTROL POINT WEIGHTS
    offset = bpinfo[3] + bpinfo[6] + 3 * bpinfo[2] * bpinfo[5];
    for (int ii = 0; ii < bpinfo[2] * bpinfo[5]; ++ii){
      bprv[offset + ii] = newW[faceWStartRow + ii];
    }
  
    // Get Context from FACE
    context = NULL;
    PetscCall(EG_getContext(face, &context));
	
    // Create New Surface
    ego newgeom;
    PetscCall(EG_makeGeometry(context, SURFACE, BSPLINE, NULL, bpinfo, bprv, &newgeom));
  
    // Create new FACE based on new SURFACE geometry
    double data[4];
    int    periodic;
    PetscCall(EG_getRange(newgeom, data, &periodic));

    ego newface;
    PetscCall(EG_makeFace(newgeom, SFORWARD, data, &newface));
    newFaces[jj] = newface;
  }
  EG_free(fobjs);
  
  // Sew New Faces together to get a new model
  ego newmodel;
  PetscCall(EG_sewFaces(Nf, newFaces, 0.0, 0, &newmodel));
    
  // Get the total number of NODEs on the original geometry. (This will be the same for the new geometry)
  int  totalNumNode;
  ego *nobjTotal;
  PetscCall(EG_getBodyTopos(body, NULL, NODE, &totalNumNode, &nobjTotal));
  EG_free(nobjTotal);
  
  // Initialize vector to store equivalent NODE indices between the 2 goemetrys
  // FORMAT :: vector index is the Original Geometry's NODE ID, the vector Value is the New Geometry's NODE ID
  int nodeIDEquiv[totalNumNode+1];
  
  // Now we need to Map the NODE and EDGE IDs from each Model
  PetscCall(EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs));
  
  // New CAD
  ego *newbodies, newgeomtest, *nfobjs;
  int  nNf, newNb, newoclass, newmtype, *newsenses;
  PetscCall(EG_getTopology(newmodel, &newgeomtest, &newoclass, &newmtype, NULL, &newNb, &newbodies, &newsenses));
  
  ego newbody = newbodies[0];
  PetscCall(EG_getBodyTopos(newbody, NULL, FACE, &nNf, &nfobjs));
  
  if (newNb > 1 ) {PetscCall(PetscPrintf(PETSC_COMM_SELF, "  ERROR :: newNb > 1 || newNb = %d \n", newNb));}
  
  // Find Equivalent Nodes
  for (int ii = 0; ii < Nf; ++ii) {
    double  fdata[4];
    int     peri;
    
    // Get Current FACE [u, v] Ranges
    PetscCall(EG_getRange(fobjs[ii], fdata, &peri));
  
    // Equate NODE IDs between 2 FACEs by working through (u, v) limits of FACE
    for (int jj = 0; jj < 2; ++jj){
      for (int kk = 2; kk < 4; ++kk) {
        double params[2] = {fdata[jj], fdata[kk]};
        double eval[18];
        PetscCall(EG_evaluate(fobjs[ii], params, eval));
  
        // Original Body
        ego *nobjsOrigFace;
        int  origNn;
        PetscCall(EG_getBodyTopos(body, fobjs[ii], NODE, &origNn, &nobjsOrigFace));
  
        double minVal = 1.0E10;
        double evalCheck[18];
        int    equivOrigNodeID = -1;
        for (int mm = 0; mm < origNn; ++mm) {
          double delta = 1.0E10;
          PetscCall(EG_evaluate(nobjsOrigFace[mm], NULL, evalCheck));
          delta = sqrt(pow(evalCheck[0] - eval[0], 2.0) + pow(evalCheck[1] - eval[1], 2.0) + pow(evalCheck[2] - eval[2], 2.0));
	
          if (delta < minVal) {
            equivOrigNodeID = EG_indexBodyTopo(body, nobjsOrigFace[mm]);
            minVal = delta;
          }
        }
        EG_free(nobjsOrigFace);
        
        // New Body
        ego *nobjsNewFace;
        int  newNn;
        PetscCall(EG_getBodyTopos(newbody, nfobjs[ii], NODE, &newNn, &nobjsNewFace));
  
        minVal = 1.0E10;
        int equivNewNodeID = -1;
        for (int mm = 0; mm < newNn; ++mm) {
          double delta = 1.0E10;
          PetscCall(EG_evaluate(nobjsNewFace[mm], NULL, evalCheck));
          delta = sqrt(pow(evalCheck[0] - eval[0], 2.0) + pow(evalCheck[1] - eval[1], 2.0) + pow(evalCheck[2] - eval[2], 2.0));
	
          if (delta < minVal) {
            equivNewNodeID = EG_indexBodyTopo(newbody, nobjsNewFace[mm]);
            minVal = delta;
          }
        }
        EG_free(nobjsNewFace);

        // Store equivalent NODE IDs
        nodeIDEquiv[equivOrigNodeID] = equivNewNodeID;
      }
    }
  }

  // Find Equivalent EDGEs
  //   Get total number of EDGEs on Original Geometry
  int  totalNumEdge;
  ego *eobjsOrig;
  PetscCall(EG_getBodyTopos(body, NULL, EDGE, &totalNumEdge, &eobjsOrig));
  EG_free(eobjsOrig);
  
  //   Get total number of EDGEs on New Geometry
  int  totalNumEdgeNew;
  ego *eobjsNew;
  PetscCall(EG_getBodyTopos(newbody, NULL, EDGE, &totalNumEdgeNew, &eobjsNew));
  EG_free(eobjsNew);
  
  // Initialize EDGE ID equivalent vector
  // FORMAT :: vector index is the Original Geometry's EDGE ID, the vector Value is the New Geometry's EDGE ID
  int edgeIDEquiv[totalNumEdge + 1];
  
  // Find Equivalent EDGEs
  for (int ii = 0; ii < Nf; ++ii) {
	  // Get Original Geometry EDGE's NODEs
	  int numOrigEdge, numNewEdge;
	  PetscCall(EG_getBodyTopos(body, fobjs[ii], EDGE, &numOrigEdge, &eobjsOrig));
	  PetscCall(EG_getBodyTopos(newbody, nfobjs[ii], EDGE, &numNewEdge, &eobjsNew));
	  
	  // new loop below
	  for (int nn = 0; nn < numOrigEdge; ++nn) {
      ego origEdge = eobjsOrig[nn];
      ego geomEdgeOrig, *nobjsOrig;
      int oclassEdgeOrig, mtypeEdgeOrig;
      int NnOrig, *nsensesEdgeOrig;
      PetscCall(EG_getTopology(origEdge, &geomEdgeOrig, &oclassEdgeOrig, &mtypeEdgeOrig, NULL, &NnOrig, &nobjsOrig, &nsensesEdgeOrig)); 
	  
      PetscBool isSame = PETSC_FALSE;
      for (int jj = 0; jj < numNewEdge; ++jj) {
        ego newEdge = eobjsNew[jj];
        ego geomEdgeNew, *nobjsNew;
        int oclassEdgeNew, mtypeEdgeNew;
        int NnNew, *nsensesEdgeNew;
        PetscCall(EG_getTopology(newEdge, &geomEdgeNew, &oclassEdgeNew, &mtypeEdgeNew, NULL, &NnNew, &nobjsNew, &nsensesEdgeNew)); 
		
        if (mtypeEdgeOrig == mtypeEdgeNew) {
          // Only operate if the EDGE types are the same
          for (int kk = 0; kk < NnNew; ++kk) {
            int nodeIDOrigGeom, nodeIDNewGeom;
            nodeIDOrigGeom = EG_indexBodyTopo(body, nobjsOrig[kk]);
            nodeIDNewGeom = EG_indexBodyTopo(newbody, nobjsNew[kk]);
            
            if (nodeIDNewGeom == nodeIDEquiv[nodeIDOrigGeom]) {
              isSame = PETSC_TRUE;
            } else {
              isSame = PETSC_FALSE;
              kk = NnNew;  // skip ahead because first NODE failed test and order is important
            }
          }
			
          if (isSame == PETSC_TRUE) {
            int edgeIDOrig, edgeIDNew;
            edgeIDOrig = EG_indexBodyTopo(body, origEdge);
            edgeIDNew = EG_indexBodyTopo(newbody, newEdge);
            edgeIDEquiv[edgeIDOrig] = edgeIDNew;
            jj = numNewEdge;
          }
        }
      }
	  }
    EG_free(eobjsOrig);
    EG_free(eobjsNew);	  
  }
  EG_free(fobjs);
  EG_free(nfobjs); 

  // Modify labels to point to the IDs on the new Geometry
  IS	isNodeID, isEdgeID;
  
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));
  
  PetscCall(ISCreateGeneral(comm, totalNumNode+1, nodeIDEquiv, PETSC_COPY_VALUES, &isNodeID));
  PetscCall(ISCreateGeneral(comm, totalNumEdge+1, edgeIDEquiv, PETSC_COPY_VALUES, &isEdgeID));
  PetscCall(DMLabelPermuteValues(vertexLabel, isNodeID));
  PetscCall(DMLabelPermuteValues(edgeLabel, isEdgeID));
  PetscCall(ISDestroy(&isNodeID));
  PetscCall(ISDestroy(&isEdgeID));

  // Attempt to point to the new geometry
  PetscCall(PetscContainerSetPointer(modelObj, newmodel));
  
  // save updated model to file
  if (saveGeom == PETSC_TRUE && stpName != NULL) {PetscCall(EG_saveModel(newmodel, stpName));}
  
  // Inflate Mesh to EGADS Model
  if (autoInflate == PETSC_TRUE) {PetscCall(DMPlexInflateToEGADSGeomModel_tuv(dm));}
  
  PetscFunctionReturn(0);
}


/*@C
  DMPlexGetEGADSGeomModel_tuv - Gets the [t] (EDGES) and [u, v] (FACES) geometry parameters of DM points 
                                that are associated geometry relationships. Requires a DM with a EGADS
                                model attached.

  Collective

  Input Parameters:
. dm      - The DM object representing the mesh with PetscContainer containing an EGADS geometry model


  Output Parameter:
. dm       - The DM object representing the mesh with PetscContainers containing the associated [t] and
             [u, v] parameter data for associated DM points.

  Level: intermediate

.seealso:  DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSLiteFromFile(), DMPlexGeomDataAndGrads()
@*/
PetscErrorCode DMPlexGetEGADSGeomModel_tuv(DM dm)
{
#if defined(PETSC_HAVE_EGADS)
  /* EGADS Variables */
  ego            model, geom, body, face, edge;
  ego           *bodies;
  int            Nb, oclass, mtype, *senses;
  double         result[4];
  /* PETSc Variables */
  DM             cdm;
  PetscContainer modelObj;
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  Vec            coordinates;
  PetscScalar   *coords;
  PetscInt       bodyID, faceID, edgeID, vertexID;
  PetscInt       cdim, vStart, vEnd, v;
  //PetscErrorCode ierr;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_EGADS)
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
  if (!modelObj) PetscFunctionReturn(0);
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinates));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));

  PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));

  PetscCall(DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd));
  PetscCall(VecGetArrayWrite(coordinates, &coords));
  
  // Define t, u, v arrays to be stored in a PetscContainer after populated
  PetscScalar  *t_point, *u_point, *v_point;
  PetscMalloc3(vEnd - vStart, &t_point, vEnd - vStart, &u_point, vEnd - vStart, &v_point);
  
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
      PetscCall(EG_objectBodyTopo(body, EDGE, edgeID, &edge));
      PetscCall(EG_invEvaluate(edge, vcoords, params, result)); // Get (t) of nearest point on EDGE
      t_point[v - vStart] = params[0];
      u_point[v - vStart] = 0.0;
      v_point[v - vStart] = 0.0;
      //for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (faceID > 0) {
      /* Snap to FACE at nearest location */
      double params[2];
      PetscCall(EG_objectBodyTopo(body, FACE, faceID, &face));
      PetscCall(EG_invEvaluate(face, vcoords, params, result)); // Get (x,y,z) of nearest point on FACE
      t_point[v - vStart] = 0.0;
      u_point[v - vStart] = params[0];
      v_point[v - vStart] = params[1];
      //for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else {
      t_point[v - vStart] = 0.0;
      u_point[v - vStart] = 0.0;
      v_point[v - vStart] = 0.0;
    }
  }
  PetscCall(VecRestoreArrayWrite(coordinates, &coords));
  /* Clear out global coordinates */
  PetscCall(VecDestroy(&dm->coordinates));
  
  /* Store in PetscContainters */
  {
    PetscContainer t_pointObj, u_pointObj, v_pointObj;

    PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Edge t Parameter", (PetscObject *) &t_pointObj));
    if (!t_pointObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &t_pointObj));
      PetscCall(PetscContainerSetPointer(t_pointObj, t_point));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Point - Edge t Parameter", (PetscObject) t_pointObj));
      PetscCall(PetscContainerDestroy(&t_pointObj));
    } else {
      PetscCall(PetscContainerSetPointer(t_pointObj, t_point));
    } 
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Face u Parameter", (PetscObject *) &u_pointObj));
    if (!u_pointObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &u_pointObj));
      PetscCall(PetscContainerSetPointer(u_pointObj, u_point));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Point - Face u Parameter", (PetscObject) u_pointObj));
      PetscCall(PetscContainerDestroy(&u_pointObj));
    } else {
      PetscCall(PetscContainerSetPointer(u_pointObj, u_point));
    } 
    
    PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Face v Parameter", (PetscObject *) &v_pointObj));
    if (!v_pointObj) {
      PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &v_pointObj));
      PetscCall(PetscContainerSetPointer(v_pointObj, v_point));
      PetscCall(PetscObjectCompose((PetscObject) dm, "Point - Face v Parameter", (PetscObject) v_pointObj));
      PetscCall(PetscContainerDestroy(&v_pointObj));
    } else {
      PetscCall(PetscContainerSetPointer(v_pointObj, v_point));
    } 
  }
#endif
  PetscFunctionReturn(0);
}


/*@C
  DMPlexInflateToEGADSGeomModel_tuv - Inflates the DM to the associated underlying geometry using the
                                      [t] {EDGES) and [u, v] (FACES} associated parameters. Requires a
                                      DM with an EGADS model attached and a previous call to 
                                      DMPlexGetEGADSGeomModel_tuv()

  Collective

  Input Parameters:
. dm      - The DM object representing the mesh with PetscContainer containing an EGADS geometry model


  Output Parameter:
. dm       - The updated DM object inflated to the associated underlying geometry. This updates the [x, y, z]
             coordinates of DM points associated with geometry.

  Level: intermediate

.seealso:  DMPLEX, DMCreate(), DMPlexCreateEGADS(), DMPlexCreateEGADSLiteFromFile(), DMPlexGeomDataAndGrads(), DMPlexGetEGADSGeomModel_tuv()
@*/
PetscErrorCode DMPlexInflateToEGADSGeomModel_tuv(DM dm)
{
#if defined(PETSC_HAVE_EGADS)
  /* EGADS Variables */
  ego            model, geom, body, face, edge, vertex;
  ego           *bodies;
  int            Nb, oclass, mtype, *senses;
  double         result[18], params[2];
  /* PETSc Variables */
  DM             cdm;
  PetscContainer modelObj;
  PetscContainer t_pointObj, u_pointObj, v_pointObj;
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  Vec            coordinates;
  PetscScalar   *coords;
  PetscScalar   *t_point, *u_point, *v_point;
  PetscInt       bodyID, faceID, edgeID, vertexID;
  PetscInt       cdim, d, vStart, vEnd, v;
  //PetscErrorCode ierr;
#endif

  PetscFunctionBegin;
#if defined(PETSC_HAVE_EGADS)
  PetscCall(PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj));
  PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Edge t Parameter", (PetscObject *) &t_pointObj));
  PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Face u Parameter", (PetscObject *) &u_pointObj));
  PetscCall(PetscObjectQuery((PetscObject) dm, "Point - Face v Parameter", (PetscObject *) &v_pointObj));
  
  if (!modelObj) PetscFunctionReturn(0);
  if (!t_pointObj) PetscFunctionReturn(0);
  if (!u_pointObj) PetscFunctionReturn(0);
  if (!v_pointObj) PetscFunctionReturn(0);
  
  PetscCall(DMGetCoordinateDim(dm, &cdim));
  PetscCall(DMGetCoordinateDM(dm, &cdm));
  PetscCall(DMGetCoordinatesLocal(dm, &coordinates));
  PetscCall(DMGetLabel(dm, "EGADS Body ID", &bodyLabel));
  PetscCall(DMGetLabel(dm, "EGADS Face ID", &faceLabel));
  PetscCall(DMGetLabel(dm, "EGADS Edge ID", &edgeLabel));
  PetscCall(DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel));
  
  PetscCall(PetscContainerGetPointer(t_pointObj, (void **) &t_point));
  PetscCall(PetscContainerGetPointer(u_pointObj, (void **) &u_point));
  PetscCall(PetscContainerGetPointer(v_pointObj, (void **) &v_point));

  PetscCall(PetscContainerGetPointer(modelObj, (void **) &model));
  PetscCall(EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses));

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
    if (vertexID > 0) {
      /* Snap to Vertices */
      PetscCall(EG_objectBodyTopo(body, NODE, vertexID, &vertex));
      PetscCall(EG_evaluate(vertex, NULL, result));
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (edgeID > 0) {
      /* Snap to EDGE */
      params[0] = t_point[v - vStart];
      PetscCall(EG_objectBodyTopo(body, EDGE, edgeID, &edge));
      PetscCall(EG_evaluate(edge, params, result));
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    } else if (faceID > 0) {
      /* Snap to FACE */
      params[0] = u_point[v - vStart];
      params[1] = v_point[v - vStart];
      PetscCall(EG_objectBodyTopo(body, FACE, faceID, &face));
      PetscCall(EG_evaluate(face, params, result));
      for (d = 0; d < cdim; ++d) vcoords[d] = result[d];
    }
  }
  PetscCall(VecRestoreArrayWrite(coordinates, &coords));
  /* Clear out global coordinates */
  PetscCall(VecDestroy(&dm->coordinates));
#endif
  PetscFunctionReturn(0);
}
