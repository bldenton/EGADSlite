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
    ierr = PetscPrintf(PETSC_COMM_SELF, "--- IN FACE REFINEMENT --- \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    point p = %d \n", p); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    Face ID = %d \n", faceID); CHKERRQ(ierr);
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);    // Corresponding FACE ID for DMPlex edge
    
    /* Get Vector with all DMPlex Node Coordinates */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Local Vertices Coordinates \n");CHKERRQ(ierr);
    Vec    coords;
    PetscInt vecSize;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &coords); CHKERRQ(ierr);
    ierr = VecSetType(coords, VECSTANDARD); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(dm, &coords); CHKERRQ(ierr);
    
    ierr = VecGetSize(coords, &vecSize); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " vecSize = %d \n", vecSize); CHKERRQ(ierr);  // DEBUG Statement
    
    ierr = VecAssemblyBegin(coords); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(coords); CHKERRQ(ierr);
    
    /* Get the 2 DMPlex Faces connected to the DMPlex Edge */
    PetscInt pSupportSize;
    PetscInt cAconeSize, cBconeSize;
    PetscInt cAeAconeSize, cAeBconeSize;
    PetscInt cBeAconeSize, cBeBconeSize;
    const PetscInt *pSupport = NULL;
    const PetscInt *cAcone = NULL, *cAconeOrient = NULL;
    const PetscInt *cBcone = NULL, *cBconeOrient = NULL;
    const PetscInt *cAeAcone = NULL, *cAeBcone = NULL;
    const PetscInt *cBeAcone = NULL, *cBeBcone = NULL;
    
    ierr = DMPlexGetSupportSize(dm, p, &pSupportSize); CHKERRQ(ierr);
    ierr = DMPlexGetSupport(dm, p, &pSupport); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Face ID 1 = %d :: Face ID 2 = %d \n", pSupport[0], pSupport[1]);CHKERRQ(ierr);
    
    /* Get Normals from Faces attached to Edge */
    PetscReal *vol = NULL;
    PetscReal centroidA[3], normalA[3];
    PetscReal centroidB[3], normalB[3];
    ierr = DMPlexComputeCellGeometryFVM(dm, pSupport[0], vol, centroidA, normalA); CHKERRQ(ierr);
    ierr = DMPlexComputeCellGeometryFVM(dm, pSupport[1], vol, centroidB, normalB); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  FACE A \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    centroidA (x, y, z) = (%lf, %lf, %lf) \n", centroidA[0], centroidA[1], centroidA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    normalA (x, y, z) = (%lf, %lf, %lf) \n", normalA[0], normalA[1], normalA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  FACE B \n"); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    centroidB (x, y, z) = (%lf, %lf, %lf) \n", centroidB[0], centroidB[1], centroidB[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "    normalB (x, y, z) = (%lf, %lf, %lf) \n", normalB[0], normalB[1], normalB[2]); CHKERRQ(ierr);
    
    /* Get All DMPlex Edges attached to the 2 DMPlex Faces */
    /*   1st Cell */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting DMPlex Edges for Cell 1 \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, pSupport[0], &cAconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, pSupport[0], &cAcone); CHKERRQ(ierr);          // Now we have all the DMPlex Edge IDs for the face
    ierr = DMPlexGetConeOrientation(dm, pSupport[0], &cAconeOrient); CHKERRQ(ierr);  // Now we have the Orientation for the DMPlex Egdes for the face
    /*   2nd Cell */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting DMPlex Edges for Cell 2 \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, pSupport[1], &cBconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, pSupport[1], &cBcone); CHKERRQ(ierr);          // Now we have all the DMPlex Edge IDs for the face
    ierr = DMPlexGetConeOrientation(dm, pSupport[1], &cBconeOrient); CHKERRQ(ierr);  // Now we have the Orientation for the DMPlex Egdes for the face
    
    /* Get Nodes for 2 DMPlex Edges attached to each of the 2 DMPlex Faces */
    /*   1st Face */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices for Cell 1, Edge A \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, cAcone[0], &cAeAconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, cAcone[0], &cAeAcone); CHKERRQ(ierr);          // Now we have the Vertices on Edge A of the 1st Face
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d has Vertices %d and %d \n", cAcone[0], cAeAcone[0], cAeAcone[1]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices for Cell 1, Edge B \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, cAcone[1], &cAeBconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, cAcone[1], &cAeBcone); CHKERRQ(ierr);          // Now we have the Vertices on Edge B of the 1st Face
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d has Vertices %d and %d \n", cAcone[1], cAeBcone[0], cAeBcone[1]);CHKERRQ(ierr);
    /*   2nd Face */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices for Cell 2, Edge A \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, cBcone[0], &cBeAconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, cBcone[0], &cBeAcone); CHKERRQ(ierr);          // Now we have the Vertices on Edge A of the 2nd Face
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d has Vertices %d and %d \n", cBcone[0], cBeAcone[0], cBeAcone[1]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices for Cell 2, Edge B \n");CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(dm, cBcone[1], &cBeBconeSize); CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, cBcone[1], &cBeBcone); CHKERRQ(ierr);          // Now we have the Vertices on Edge B of the 2nd Face
    ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d has Vertices %d and %d \n", cBcone[1], cBeBcone[0], cBeBcone[1]);CHKERRQ(ierr);
    /*   Get DMPlex Start ID for Vertices */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Vertices Start ID \n");CHKERRQ(ierr);
    PetscInt vStart, vEnd;
    ierr = DMPlexGetHeightStratum(dm, 2, &vStart, &vEnd); CHKERRQ(ierr);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " vStart = %d \n", vStart);CHKERRQ(ierr);
    /*   Get Coordinates of Vertices defining the 1st Edge on Face 1 */
    PetscInt ix[3], ni;  //ix[3]
    PetscScalar xyzA[3], xyzB[3];  // xyzA[3], xyzB[3]
    ni = 3;     // Assume 3D Space ni = 3
    double vecA[3], vecB[3], normA[3];
    double lenA;
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 1, Edge A, Vertex A Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  ix[] start = %d \n", cAeAcone[0]-vStart);CHKERRQ(ierr);
    int pntCntr = 0;
    for (int ii = (cAeAcone[0]-vStart)*3; ii < (cAeAcone[0]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;    //was ix[ii]
      //ierr = VecGetValues(coords, ni, ix, &xyz); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge  was &xyzA
      //xyzA[pntCntr] = xyz;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzA); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 1, Edge A, Vertex B Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  ix[] start = %d \n", cAeAcone[1]-vStart);CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cAeAcone[1]-vStart)*3; ii < (cAeAcone[1]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzB); CHKERRQ(ierr);  // Get coordinates of 2st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]); CHKERRQ(ierr);
    /*   Calculate vector along 1st Face Edge */
    //if (cAconeOrient[0] < 0){
    //  for (int ii = 0; ii < 3; ++ii){
    //    vecA[ii] = xyzA[ii] - xyzB[ii];
    //  }
    //} else {
      for (int ii = 0; ii < 3; ++ii){
        vecA[ii] = xyzB[ii] - xyzA[ii];
      }
    //}
    /*   Get Coordinates of Vertices defining the 2nd Edge on Face 1 */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 1, Edge B, Vertex A Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  ix[] start = %d \n", cAeBcone[0]-vStart);CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cAeBcone[0]-vStart)*3; ii < (cAeBcone[0]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzA); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 1, Edge B, Vertex B Coordinates \n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  ix[] start = %d \n", cAeBcone[1]-vStart);CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cAeBcone[1]-vStart)*3; ii < (cAeBcone[1]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzB); CHKERRQ(ierr);  // Get coordinates of 2st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]); CHKERRQ(ierr);
    /*   Calculate vector along 1st Face Edge */
    //if (cAconeOrient[1] < 0){
    //  for (int ii = 0; ii < 3; ++ii){
    //    vecB[ii] = xyzA[ii] - xyzB[ii];
    //  }
    //} else {
      for (int ii = 0; ii < 3; ++ii){
        vecB[ii] = xyzB[ii] - xyzA[ii];
      }
    //}    
    /*   Calculate Normal of Face 1 using the Cross Product vecA x vecB */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Calculating Face 1 Normal Vector \n");CHKERRQ(ierr);
    normA[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
    normA[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
    normA[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
    lenA = sqrt(normA[0]*normA[0] + normA[1]*normA[1] + normA[2]*normA[2]);    // length of Normal for Surface A
    ierr = PetscPrintf(PETSC_COMM_SELF, "  normA (x, y, z) = (%lf, %lf, %lf) \n", normA[0], normA[1], normA[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  lenA = %lf \n", lenA);CHKERRQ(ierr);
    
    
    /*   Get Coordinates of Vertices defining the 1st Edge on Face 2 */
    //PetscInt ix[3], ni;
    //PetscScalar xyzA[3], xyzB[3];
    //ni = 3;     // Assume 3D Space
    //double vecA[3], vecB[3], normA[3];
    double normB[3];
    double lenB;
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 2, Edge A, Vertex A Coordinates \n");CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cBeAcone[0]-vStart)*3; ii < (cBeAcone[0]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzA); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]); CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 2, Edge A, Vertex B Coordinates \n");CHKERRQ(ierr); 
    pntCntr = 0; 
    for (int ii = (cBeAcone[1]-vStart)*3; ii < (cBeAcone[1]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzB); CHKERRQ(ierr);  // Get coordinates of 2st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]); CHKERRQ(ierr);
    /*   Calculate vector along 1st Face Edge */
    //if (cBconeOrient[0] < 0){
    //  for (int ii = 0; ii < 3; ++ii){
    //    vecA[ii] = xyzA[ii] - xyzB[ii];
    //  }
    //} else {
      for (int ii = 0; ii < 3; ++ii){
        vecA[ii] = xyzB[ii] - xyzA[ii];
      }
    //}
    /*   Get Coordinates of Vertices defining the 2nd Edge on Face 2 */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 2, Edge B, Vertex A Coordinates \n");CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cBeBcone[0]-vStart)*3; ii < (cBeBcone[0]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzA); CHKERRQ(ierr);  // Get coordinates of 1st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "Getting Face 2, Edge B, Vertex B Coordinates \n");CHKERRQ(ierr);
    pntCntr = 0;
    for (int ii = (cBeBcone[1]-vStart)*3; ii < (cBeBcone[1]-vStart)*3+3; ++ii){    //*** NEED TO UPDATE replace *3 with *coordDim  ***//
      ix[pntCntr] = ii;
      pntCntr = pntCntr + 1;
    }
    ierr = VecGetValues(coords, ni, ix, xyzB); CHKERRQ(ierr);  // Get coordinates of 2st point on Edge
    ierr = PetscPrintf(PETSC_COMM_SELF, "  xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]); CHKERRQ(ierr);
    /*   Calculate vector along 1st Face Edge */
    //if (cBconeOrient[1] < 0){
    //  for (int ii = 0; ii < 3; ++ii){
    //    vecB[ii] = xyzA[ii] - xyzB[ii];
    //  }
    //} else {
      for (int ii = 0; ii < 3; ++ii){
        vecB[ii] = xyzB[ii] - xyzA[ii];
      }
    //}    
    /*   Calculate Normal of Face 1 using the Cross Product vecA x vecB */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Calculating Face 2 Normal \n");CHKERRQ(ierr);
    normB[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
    normB[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
    normB[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
    lenB = sqrt(normB[0]*normB[0] + normB[1]*normB[1] + normB[2]*normB[2]);    // length of Normal for Surface B
    ierr = PetscPrintf(PETSC_COMM_SELF, "  normB (x, y, z) = (%lf, %lf, %lf) \n", normB[0], normB[1], normB[2]);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  lenB = %lf \n", lenB);CHKERRQ(ierr);
    
    /*   Now Calculate the Average Normal of the 2 Faces connected to the Edge and Normalize it */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Calculate normAvg \n");CHKERRQ(ierr);
    double normAvg[3];
    double len;
    for (int ii = 0; ii < 3; ++ii){
      normAvg[ii] = normA[ii] + normB[ii];
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
    double initPoint[3];
    for (int ii = 0; ii < 3; ++ii){
      initPoint[ii] = point[ii];
    }
    ierr = PetscPrintf(PETSC_COMM_SELF, "  initPoint (x, y, z) = (%lf, %lf, %lf) \n", initPoint[0], initPoint[1], initPoint[2]);CHKERRQ(ierr);
    
    /*   Perform Multiple Inverse Evaluations moving along the normal vector  */
    ierr = PetscPrintf(PETSC_COMM_SELF, "Perform Multiple Inverse Evaluations \n");CHKERRQ(ierr);
    ierr = EG_invEvaluate(face, point, params, result); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "  result (x, y, z) = (%lf, %lf, %lf) \n", result[0], result[1], result[2]);CHKERRQ(ierr);
    for (int jj = 0; jj < 3; ++jj){
      // Calculate vector from initial point to solution point
      for (int ii = 0; ii < 3; ++ii){
        vecA[ii] = result[ii] - initPoint[ii];
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  vecA (x, y, z) = (%lf, %lf, %lf) \n", vecA[0], vecA[1], vecA[2]);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "  initPoint (x, y, z) = (%lf, %lf, %lf) \n", initPoint[0], initPoint[1], initPoint[2]);CHKERRQ(ierr);
      
      // Dot Product between vector and normal
      len = 0.;
      for (int ii = 0; ii < 3; ++ii){
        len = len + vecA[ii]*normAvg[ii];
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  len = %lf \n", len);CHKERRQ(ierr);
      
      // Set new guess point and Inverse Evaluate
      for (int ii = 0; ii < 3; ++ii){
        point[ii] = len*normAvg[ii];
      }
      ierr = PetscPrintf(PETSC_COMM_SELF, "  point (x, y, z) = (%lf, %lf, %lf) \n", point[0], point[1], point[2]);CHKERRQ(ierr);
      ierr = EG_invEvaluate(face, point, params, result); CHKERRQ(ierr); 
      ierr = PetscPrintf(PETSC_COMM_SELF, "  result (x, y, z) = (%lf, %lf, %lf) \n", result[0], result[1], result[2]);CHKERRQ(ierr);
    }
        
    /* Perform Inverse Evaulation for point Snapped to Geometry - Original */
    //ierr = EG_invEvaluate(face, point, params, result);
  }
  for (d = 0; d < cdim; ++d) coords[d] = result[d];
  ierr = PetscPrintf(PETSC_COMM_SELF, "  Completed \n"); CHKERRQ(ierr);
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
