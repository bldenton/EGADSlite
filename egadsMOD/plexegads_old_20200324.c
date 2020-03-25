#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel        bodyLabel, faceLabel, edgeLabel;
  PetscContainer modelObj;
  PetscInt       cdim, d, bodyID, faceID, edgeID;
  ego           *bodies;
  ego            model, geom, body, face, edge;
  double         point[3], params[3], result[3];
  int            Nb, oclass, mtype, *senses;
  PetscErrorCode ierr;

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
  for (d = 0; d < cdim; ++d) point[d] = initCoords[d];
  if (edgeID >= 0) {
    /* Snap to EDGE Geometry */
    ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    ierr = EG_invEvaluate(edge, point, params, result);
  } else {
    /* Snap to FACE Geometry searching along DMPlex normal */
    ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------------- \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "--- IN FACE REFINEMENT --- \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "-------------------------- \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);    // Corresponding FACE ID for DMPlex edge
    
    /* Get Vector with all DMPlex Node Coordinates */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Local Vertices Coordinates \n");CHKERRQ(ierr);
    Vec    vCoords;
    PetscInt vecSize;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &vCoords); CHKERRQ(ierr);
    ierr = VecSetType(vCoords, VECSTANDARD); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &vCoords); CHKERRQ(ierr);
    
    ierr = VecGetSize(vCoords, &vecSize); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " vecSize = %d \n", vecSize); CHKERRQ(ierr);  // DEBUG Statement
    
    ierr = VecAssemblyBegin(vCoords); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vCoords); CHKERRQ(ierr);
    
    /* Get the 2 DMPlex Faces connected to the DMPlex Edge */
    PetscInt pSupportSize;
    //PetscInt cAconeSize, cBconeSize;
    //PetscInt cAeAconeSize, cAeBconeSize;
    //PetscInt cBeAconeSize, cBeBconeSize;
    const PetscInt *pSupport = NULL;
    //const PetscInt *cAcone = NULL, *cAconeOrient = NULL;
    //const PetscInt *cBcone = NULL, *cBconeOrient = NULL;
    //const PetscInt *cAeAcone = NULL, *cAeBcone = NULL;
    //const PetscInt *cBeAcone = NULL, *cBeBcone = NULL;
    
    ierr = DMPlexGetSupportSize(dm, p, &pSupportSize); CHKERRQ(ierr);
    ierr = DMPlexGetSupport(dm, p, &pSupport); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Face A ID = %d :: Face B ID = %d \n", pSupport[0], pSupport[1]);CHKERRQ(ierr);
    
    /* Get Normals from Faces attached to Edge */
    PetscReal *vol = NULL;
    PetscReal centroidA[3], normalA[3];
    PetscReal centroidB[3], normalB[3];
    ierr = DMPlexComputeCellGeometryFVM(dm, pSupport[0], vol, centroidA, normalA); CHKERRQ(ierr);
    ierr = DMPlexComputeCellGeometryFVM(dm, pSupport[1], vol, centroidB, normalB); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  FACE %d \n", pSupport[0]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    centroidA (x, y, z) = (%lf, %lf, %lf) \n", centroidA[0], centroidA[1], centroidA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    normalA (x, y, z) = (%lf, %lf, %lf) \n", normalA[0], normalA[1], normalA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  FACE %d \n", pSupport[1]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    centroidB (x, y, z) = (%lf, %lf, %lf) \n", centroidB[0], centroidB[1], centroidB[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    normalB (x, y, z) = (%lf, %lf, %lf) \n", normalB[0], normalB[1], normalB[2]); CHKERRQ(ierr);
    
    /*   Now Calculate the Average Normal of the 2 Faces connected to the Edge and Normalize it */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Calculate normAvg \n");CHKERRQ(ierr);
    double normAvg[3];
    double len;
    for (int ii = 0; ii < 3; ++ii){
      normAvg[ii] = normalA[ii] + normalB[ii];
    }
    ierr = PetscPrintf(PETSC_COMM_SELF, "  normAvg (x, y, z) = (%lf, %lf, %lf) \n", normAvg[0], normAvg[1], normAvg[2]);CHKERRQ(ierr);
    
    len = sqrt(normAvg[0]*normAvg[0] + normAvg[1]*normAvg[1] + normAvg[2]*normAvg[2]);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  normAvg_len = %lf \n", len);CHKERRQ(ierr);
    
    for (int ii = 0; ii < 3; ++ii){
      normAvg[ii] = normAvg[ii] / len;
    }
    ierr = PetscPrintf(PETSC_COMM_SELF, "  normAvg/len (x, y, z) = (%lf, %lf, %lf) \n", normAvg[0], normAvg[1], normAvg[2]);CHKERRQ(ierr);
    
    /*   Save Initial guess Point on Edge. Use as base point  */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Defining Initial Guess Point (initPoint) \n");CHKERRQ(ierr);
    double initPoint[3], vecA[3];
    for (int ii = 0; ii < 3; ++ii){
      initPoint[ii] = point[ii];
    }
    ierr = PetscPrintf(PETSC_COMM_SELF, "  initPoint (x, y, z) = (%lf, %lf, %lf) \n", initPoint[0], initPoint[1], initPoint[2]);CHKERRQ(ierr);
    
    /*   Perform Multiple Inverse Evaluations moving along the normal vector  */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Perform Multiple Inverse Evaluations \n");CHKERRQ(ierr);
    ierr = EG_invEvaluate(face, point, params, result); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  result (x, y, z) = (%lf, %lf, %lf) \n", result[0], result[1], result[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", params[0], params[1]); CHKERRQ(ierr);
    for (int jj = 0; jj < 10; ++jj){
      // Calculate vector from initial point to solution point
      for (int ii = 0; ii < 3; ++ii){
        vecA[ii] = result[ii] - initPoint[ii];
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  vecA (x, y, z) = (%lf, %lf, %lf) \n", vecA[0], vecA[1], vecA[2]); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "  initPoint (x, y, z) = (%lf, %lf, %lf) \n", initPoint[0], initPoint[1], initPoint[2]); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "  initCoords (x, y, z) = (%lf, %lf, %lf) \n", initCoords[0], initCoords[1], initCoords[2]);CHKERRQ(ierr);
      
      // Dot Product between vector and normal
      len = 0.;
      for (int ii = 0; ii < 3; ++ii){
        len = len + vecA[ii]*normAvg[ii];
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  len = %lf \n", len);CHKERRQ(ierr);
      
      // Set new guess point and Inverse Evaluate
      for (int ii = 0; ii < 3; ++ii){
        point[ii] = abs(len)*normAvg[ii] + initPoint[ii]; // initPoint[] = result[] originally
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  point (x, y, z) = (%lf, %lf, %lf) \n", point[0], point[1], point[2]);CHKERRQ(ierr);
      ierr = EG_invEvaluate(face, point, params, result); CHKERRQ(ierr); 
      ierr = PetscPrintf(PETSC_COMM_SELF, "  result (x, y, z) = (%lf, %lf, %lf) \n", result[0], result[1], result[2]);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", params[0], params[1]); CHKERRQ(ierr);
    }
    
    /* -------------------------------------------------- */    
    /*  Operate on the (u, v) parameter space of the Face */
    /* -------------------------------------------------- */
    PetscInt eConeSize;
    const PetscInt *eCone = NULL;
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertex for DMPlex Edge \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, p, &eConeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, p, &eCone); CHKERRQ(ierr);          // Now we have Vertex IDs for the Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d has Vertices %d and %d \n", p, eCone[0], eCone[1]); CHKERRQ(ierr);
    
    /*   Get DMPlex Start ID for Vertices */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices Start ID \n"); CHKERRQ(ierr);
    PetscInt vStart, vEnd;
    ierr = DMPlexGetHeightStratum(dm, 2, &vStart, &vEnd); CHKERRQ(ierr); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  vStart = %d \n", vStart); CHKERRQ(ierr);
    
    /*   Get Coordinates of Vertices defining the 1st Edge on Face 1 */
    PetscInt ix[3], ni;
    PetscScalar xyzA[3], xyzB[3];
    ni = 3;     // Assume 3D Space ni = 3
    //double vecA[3], vecB[3];
    //double lenA;
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Getting Vertex A Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    ix[] start = %d \n", eCone[0]-vStart);CHKERRQ(ierr);
    int pntCntr = 0;
    for (int ii = (eCone[0]-vStart)*3; ii < (eCone[0]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;    //was ix[ii]
      //ierr = VecGetValues(coords, ni, ix, &xyz); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge  was &xyzA
      //xyzA[pntCntr] = xyz;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(vCoords, ni, ix, xyzA); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]); CHKERRQ(ierr);
    
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Getting Vertex B Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    ix[] start = %d \n", eCone[1]-vStart);CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (eCone[1]-vStart)*3; ii < (eCone[1]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;    //was ix[ii]
      //ierr = VecGetValues(coords, ni, ix, &xyz); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge  was &xyzA
      //xyzA[pntCntr] = xyz;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(vCoords, ni, ix, xyzB); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]); CHKERRQ(ierr);
    
    /*   Get Parameters for Vertex A of Edge */
    double    paramsA[3], resultA[3];
    ierr = EG_invEvaluate(face, xyzA, paramsA, resultA); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  resultA (x, y, z) = (%lf, %lf, %lf) \n", resultA[0], resultA[1], resultA[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "             (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
    
    /*   Get Parameters for Vertex B of Edge */
    double    paramsB[3], resultB[3];
    ierr = EG_invEvaluate(face, xyzB, paramsB, resultB); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  resultB (x, y, z) = (%lf, %lf, %lf) \n", resultB[0], resultB[1], resultB[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "             (u, v) = (%lf, %lf) \n", paramsB[0], paramsB[1]); CHKERRQ(ierr);
    
    /* Calculate (u, v) parameters for refined point */
    double paramsNew[2];
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;
    paramsNew[1] = (paramsA[1] + paramsB[1]) / 2.;
    
    /* Calculate new (x, y, z) coordinates based on (u, v) parameter space */
    double eval[18];
    ierr = EG_evaluate(face, paramsNew, eval); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
    
    for (int ii = 0; ii < 3; ++ii){
      result[ii] = eval[ii];
    }
    
    /* Indicate End of Refinement for this Edge */
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Completed \n"); CHKERRQ(ierr);
  }
  for (d = 0; d < cdim; ++d) coords[d] = result[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
