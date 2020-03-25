#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel         bodyLabel, faceLabel, edgeLabel;
  PetscContainer  modelObj;
  PetscInt        cdim, d, bodyID, faceID, edgeID;
  PetscInt        vecSize, vStart, vEnd, eConeSize, ni;
  PetscInt        ix[3];
  const PetscInt *eCone = NULL;
  PetscScalar     xyzA[3], xyzB[3];
  Vec             vCoords;
  ego            *bodies;
  ego             model, geom, body, face, edge;
  double          point[3], params[3], result[3];
  double          resultA[3], resultB[3];
  double          eval[18];
  int             Nb, oclass, mtype, *senses, pntCntr;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];
  
  //for (d = 0; d < cdim; ++d) point[d] = initCoords[d];  // Don't need
  
  /* Get Vector with all DMPlex Node Coordinates */
  ierr = VecCreate(PETSC_COMM_WORLD, &vCoords); CHKERRQ(ierr);
  ierr = VecSetType(vCoords, VECSTANDARD); CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &vCoords); CHKERRQ(ierr);    
  ierr = VecGetSize(vCoords, &vecSize); CHKERRQ(ierr);
    
  ierr = VecAssemblyBegin(vCoords); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vCoords); CHKERRQ(ierr);
    
  /* Get DMPlex Start ID for Vertices */
  ierr = DMPlexGetHeightStratum(dm, 2, &vStart, &vEnd); CHKERRQ(ierr); CHKERRQ(ierr);
 
  /* Get ID's of Vertex at each end of the Edge */
  ierr = DMPlexGetConeSize(dm, p, &eConeSize); CHKERRQ(ierr);
  ierr = DMPlexGetCone(dm, p, &eCone); CHKERRQ(ierr);
 
  /* Get Coordinates of Vertices defining the 1st Edge on Face 1 */
  ni = cdim;     // Assume 3D Space ni = 3, need to update to reference cdim
  pntCntr = 0;
  for (int ii = (eCone[0]-vStart)*cdim; ii < (eCone[0]-vStart)*cdim+cdim; ++ii){
    ix[pntCntr] = ii;    //was ix[ii]
    pntCntr = pntCntr + 1;
  }
  ierr = VecGetValues(vCoords, ni, ix, xyzA); CHKERRQ(ierr);    // Vertex A coordinates (x,y,z)

  pntCntr = 0;
  for (int ii = (eCone[1]-vStart)*cdim; ii < (eCone[1]-vStart)*cdim+cdim; ++ii){
    ix[pntCntr] = ii;
    pntCntr = pntCntr + 1;
  }
  ierr = VecGetValues(vCoords, ni, ix, xyzB); CHKERRQ(ierr);  // Vertex B coordinates (x,y,z)
    
  /* Snap new Vertex to Geometry */
  if (edgeID >= 0) {
    /* Snap to EDGE Geometry */
    double   paramsA[1], paramsB[1], paramsC[1];
    double   paramsNew[1];
    ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    ierr = EG_invEvaluate(edge, xyzA, paramsA, resultA); CHKERRQ(ierr);   // t- parameter for Vertex A
    ierr = EG_invEvaluate(edge, xyzB, paramsB, resultB); CHKERRQ(ierr);   // t- parameter for Vertex B
    double range[2];
    int    peri;
    ierr = EG_getRange(edge, range, &peri); CHKERRQ(ierr);
    
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;          // t- parameter for New Vertex at midpoint
    
    if ((paramsNew[0] >= range[0]) && (paramsNew[0] <= range[1])){
      ierr = EG_evaluate(edge, paramsNew, eval); CHKERRQ(ierr);  // (x,y,z) coordinates for New Vertex at midpoint
    } else {
      ierr = EG_invEvaluate(edge, initCoords, paramsC, eval); CHKERRQ(ierr);
    }
    
    //if (p == 238){
    ierr = PetscPrintf(PETSC_COMM_SELF, "  -- EDGE Refinement -- \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", edgeID); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsA[0]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsB[0]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsNew[0]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    range(tstart, tend) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
    //}
    
    //ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    //ierr = EG_invEvaluate(edge, point, params, eval);
    
    
  } else {
    /* Snap to FACE Geometry searching along DMPlex normal */
    double   paramsA[2], paramsB[2], paramsC[2];
    double   paramsNew[2];
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);    // Corresponding FACE ID for DMPlex edge
    
    ierr = EG_invEvaluate(face, xyzA, paramsA, resultA); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex A
    ierr = EG_invEvaluate(face, xyzB, paramsB, resultB); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex B
    
    /* Get (u,v) range for FACE */
    double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
    
    /* Calculate (u, v) parameters for refined point */
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;
    paramsNew[1] = (paramsA[1] + paramsB[1]) / 2.;
    
    if ((paramsNew[0] >= range[0] && paramsNew[0] <= range[1]) && (paramsNew[1] >= range[2] && paramsNew[1] <= range[3])){
      ierr = EG_evaluate(face, paramsNew, eval); CHKERRQ(ierr);
    } else {
      ierr = EG_invEvaluate(edge, initCoords, paramsC, eval); CHKERRQ(ierr);
    }
    
    ierr = PetscPrintf(PETSC_COMM_SELF, "  -- FACE Refinement -- \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", edgeID); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsNew[0], paramsNew[1]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    range(ustart, uend) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "         (vstart, vend) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
    
    ///* Calculate new (x, y, z) coordinates based on (u, v) parameter space */
    //ierr = EG_evaluate(face, paramsNew, eval); CHKERRQ(ierr);
  }
  for (d = 0; d < cdim; ++d) coords[d] = eval[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
