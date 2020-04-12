#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, tLabel, uLabel, vLabel;
  PetscContainer  modelObj;
  PetscInt        cdim, d, bodyID, faceID, edgeID;
  PetscScalar     tParamA, uParamA, vParamA, tParamB, uParamB, vParamB;
  PetscInt        vecSize, vStart, vEnd, eConeSize, ni;
  PetscInt        ix[3];
  const PetscInt *eCone = NULL;
  PetscScalar     xyzA[3], xyzB[3];
  Vec             vCoords;
  ego            *bodies;
  ego             model, geom, body, face, edge;
  double          point[3], params[3], result[3];
  double          resultA[3], resultB[3];
  double          eval[18], evalC[18];
  double          errTolr = 1.0e-4;
  int             Nb, oclass, mtype, *senses, pntCntr;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS t param", &tLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS u param", &uLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS v param", &vLabel);CHKERRQ(ierr);
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
  
  /* -------------------------- */
  /*   Trial Debug              */
  /* -------------------------- */
  //ierr = DMLabelGetValue(tLabel, eCone[0], &tParamA);CHKERRQ(ierr);
  //ierr = DMLabelGetValue(uLabel, eCone[0], &uParamA);CHKERRQ(ierr);
  //ierr = DMLabelGetValue(vLabel, eCone[0], &vParamA);CHKERRQ(ierr);
  //ierr = DMLabelGetValue(tLabel, eCone[1], &tParamB);CHKERRQ(ierr);
  //ierr = DMLabelGetValue(uLabel, eCone[1], &uParamB);CHKERRQ(ierr);
  //ierr = DMLabelGetValue(vLabel, eCone[1], &vParamB);CHKERRQ(ierr);
  //
  //if (uParamA != NAN) {
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "      uParamA = %lf \n", uParamA); CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "      vParamA = %lf \n", vParamA); CHKERRQ(ierr);
  //}
  //
  //if (uParamB != NAN) {
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "      uParamB = %lf \n", uParamB); CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_SELF, "      vParamB = %lf \n", vParamB); CHKERRQ(ierr);
  //}
  /* End of Trial Debug */
 
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
    
    /* Debug */
    //if (p == 238){
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  -- EDGE Refinement -- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", edgeID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsA[0]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsB[0]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsNew[0]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    range(tstart, tend) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
    //}
    
    //ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    //ierr = EG_invEvaluate(edge, point, params, eval);
    
    
  } else {
    /* Snap to FACE Geometry searching along DMPlex normal */
    double   paramsA[2], paramsB[2], paramsC[2];
    double   paramsNew[2];
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);    // Corresponding FACE ID for DMPlex edge
    
    /* Check tolerance */
    //double  tol;
    //ierr = EG_getTolerance(face, &tol); CHKERRQ(ierr);
    //tol = 1.0e-10;
    //ierr = EG_setTolerance(face, 1.0e-10); CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_SELF, "    tol = %.6e \n", tol); CHKERRQ(ierr);
    //ierr = EG_getTolerance(face, &tol); CHKERRQ(ierr);
    /* ------- */
    
    // Instead of doing EG_invEvaluate. Do we want to use t- and (u,v) DMlabels?
    ierr = EG_invEvaluate(face, xyzA, &paramsA, &resultA); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex A
    ierr = EG_invEvaluate(face, xyzB, &paramsB, &resultB); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex B
    
    // Check if we get the same xyz Location back -- if not we have to correct
    ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
    double delta;
    double dx, dy, dz, du, dv;
    delta = 0.;
    
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]))/xyzA[ii];
    }
    delta = sqrt(delta);
    
    int  xyzCheck = 0;
    if (delta > errTolr){
        xyzCheck = 1;      
        
    //  paramsA[0] = paramsA[0] / 2.;
    //  paramsA[1] = paramsA[1] / 2.;
    //  ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
    //  for (int jj=0; jj<3; ++jj){
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    PointA doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    //    
    //    dx = eval[0] - xyzA[0];
    //    dy = eval[1] - xyzA[1];
    //    dv = (eval[3]*dy - eval[4]*dx)/(eval[3]*eval[7] - eval[4]*eval[6]);
    //    du = (dx - eval[6]*dv)/eval[3];
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    du = %.6e \n", du); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    dv = %.6e \n", dv); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[3] = %.6e \n", eval[3]); CHKERRQ(ierr);
    //    ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[4] = %.6e \n", eval[4]); CHKERRQ(ierr);
    //    
    //    paramsA[0] = paramsA[0] - du;
    //    paramsA[1] = paramsA[1] - dv;
    //    ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
    //    
    //    delta = 0.;
    //    for (int ii=0; ii<3; ++ii){
    //      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]));
    //    }
    //    delta = sqrt(delta);
    //  }
    }
    
    ierr = EG_evaluate(face, paramsB, eval); CHKERRQ(ierr);
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]))/xyzB[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){
        xyzCheck = 1;
        
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    PointB doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    }   
    // End Check 
    
    /* Get (u,v) range for FACE */
    double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
    
    /* Calculate (u, v) parameters for refined point */
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;
    paramsNew[1] = (paramsA[1] + paramsB[1]) / 2.;
    
    int paramTest = 0;
    double xMax, xMin, yMax, yMin, zMax, zMin;
    if ((paramsNew[0] >= range[0] && paramsNew[0] <= range[1]) && (paramsNew[1] >= range[2] && paramsNew[1] <= range[3])){
      ierr = EG_evaluate(face, paramsNew, &eval); CHKERRQ(ierr);
      
      if (xyzCheck == 1){
        for (int ii = 0; ii <3; ++ii){
          eval[ii] = initCoords[ii];
        }
      }
      
      //// 04-07-2020 :: Try to fix Refinement error //
      //ierr = PetscPrintf(PETSC_COMM_SELF, "  A  eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
      //paramTest = 1;
      //if (xyzA[0] >= xyzB[0]) {
      //  xMax = xyzA[0];
      //  xMin = xyzB[0];
      //} else {
      //  xMax = xyzB[0];
      //  xMin = xyzA[0];
      //}
      //ierr = PetscPrintf(PETSC_COMM_SELF, "  xMax, xMin = (%lf, %lf) \n", xMax, xMin);CHKERRQ(ierr);
      //
      //if (xyzA[1] >= xyzB[1]) {
      //  yMax = xyzA[1];
      //  yMin = xyzB[1];
      //} else {
      //  yMax = xyzB[1];
      //  yMin = xyzA[1];
      //}
      //ierr = PetscPrintf(PETSC_COMM_SELF, "  yMax, yMin = (%lf, %lf) \n", yMax, yMin);CHKERRQ(ierr);
      //
      //if (xyzA[2] >= xyzB[2]) {
      //  zMax = xyzA[2];
      //  zMin = xyzB[2];
      //} else {
      //  zMax = xyzB[2];
      //  zMin = xyzA[2];
      //}
      //ierr = PetscPrintf(PETSC_COMM_SELF, "  zMax, zMin = (%lf, %lf) \n", zMax, zMin);CHKERRQ(ierr);
      //
      //if ((eval[0] >= xMin && eval[0] <= xMax) && (eval[1] >= yMin && eval[1] <= yMax) && (eval[2] >= zMin && eval[2] <= zMax)){
      //  //eval[0] = (xyzA[0] + xyzB[0]) / 2.;
      //  //eval[1] = (xyzA[1] + xyzB[1]) / 2.;
      //  //eval[2] = (xyzA[2] + xyzB[2]) / 2.;
      //  //ierr = PetscPrintf(PETSC_COMM_SELF, "  B  eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
      //} else {
      //  // Do Nothing
      //  eval[0] = (xyzA[0] + xyzB[0]) / 2.;
      //  eval[1] = (xyzA[1] + xyzB[1]) / 2.;
      //  eval[2] = (xyzA[2] + xyzB[2]) / 2.;
      //  ierr = PetscPrintf(PETSC_COMM_SELF, "  B  eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
      //  ierr = EG_invEvaluate(face, eval, paramsC, eval); CHKERRQ(ierr);
      //}
      ////   END Change //
      
    } else {
      ierr = EG_invEvaluate(face, initCoords, &paramsC, &eval); CHKERRQ(ierr);
    }
    
    //if ((eCone[0] == vStart + 53 || eCone[1] == vStart + 53)) { // && (eCone[0] == vStart - 1 + 519 || eCone[1] == vStart - 1 + 519)){
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  --------------------- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  -- FACE Refinement -- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "  --------------------- \n"); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", edgeID); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    vStart = %d \n", vStart); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsNew[0], paramsNew[1]); CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "    initCoords (x, y, z) = (%lf, %lf, %lf) \n", initCoords[0], initCoords[1], initCoords[2]);CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "    evalC (x, y, z) = (%lf, %lf, %lf) \n", evalC[0], evalC[1], evalC[2]);CHKERRQ(ierr);
    //  //ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsC[0], paramsC[1]); CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "    range(ustart, uend) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
    //  ierr = PetscPrintf(PETSC_COMM_SELF, "         (vstart, vend) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
    //}
    
    ///* Calculate new (x, y, z) coordinates based on (u, v) parameter space */
    //ierr = EG_evaluate(face, paramsNew, eval); CHKERRQ(ierr);
  }
  for (d = 0; d < cdim; ++d) coords[d] = eval[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
