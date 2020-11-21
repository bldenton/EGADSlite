#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, vertexLabel;
  PetscContainer  modelObj;
  PetscInt        cdim, d, bodyID, faceID, edgeID, vertexID;
  PetscInt        vecSize, vStart, vEnd, eConeSize, ni;
  PetscInt        ix[3];
  const PetscInt *eCone = NULL;
  PetscScalar     xyz[3], xyzA[3], xyzB[3];
  Vec             vCoords;
  ego            *bodies;
  ego             model, geom, body, face, edge;
  double          resultA[3], resultB[3];
  double          range[4];
  int             peri;
  double          eval[18];
  double          errTolr = 1.0e-4;
  double          delta;
  int             Nb, oclass, mtype, *senses, pntCntr;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n Start DMPlexSnapToGeomModel \n");CHKERRQ(ierr);
  
  /* Output p ID for DEBUG */
  ierr = PetscPrintf(PETSC_COMM_SELF, "  Current DMPlex edge || p = %d \n", p); CHKERRQ(ierr);
  
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &vertexID);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  if (bodyID < 0 || bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);		// Added bodyID < 0 condition
  body = bodies[bodyID];

  
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
  
  /* DEBUG Statements */
  ierr = PetscPrintf(PETSC_COMM_SELF, "-- p = %d \n", p); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      faceID = %d \n", faceID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      edgeID = %d \n", edgeID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      vertexID = %d \n", vertexID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(vertexLabel, eCone[0], &vertexID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, eCone[0], &edgeID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, eCone[0], &faceID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      faceID = %d \n", faceID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      edgeID = %d \n", edgeID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      vertexID = %d \n", vertexID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(vertexLabel, eCone[1], &vertexID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, eCone[1], &edgeID); CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, eCone[1], &faceID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "   eCone[1] = %d \n", eCone[1]); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      faceID = %d \n", faceID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      edgeID = %d \n", edgeID); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "      vertexID = %d \n", vertexID); CHKERRQ(ierr);
  
  ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &vertexID);CHKERRQ(ierr);
  /* */
 
  /* Get Coordinates of Vertices defining the 1st Edge on Face 1 */
  ni = cdim;
  pntCntr = 0;
  for (int ii = (eCone[0]-vStart)*cdim; ii < (eCone[0]-vStart)*cdim+cdim; ++ii){
    ix[pntCntr] = ii;
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
  if (bodyID >= 0 && edgeID >= 0) {
    /* Snap to EDGE Geometry */
    double   paramsA[1], paramsB[1];
    double   paramsNew[1];
    ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    ierr = EG_invEvaluate(edge, xyzA, paramsA, resultA); CHKERRQ(ierr);   // t- parameter for Vertex A
    ierr = EG_invEvaluate(edge, xyzB, paramsB, resultB); CHKERRQ(ierr);   // t- parameter for Vertex B
    
    /* Get FACE range */
    ierr = EG_getRange(edge, range, &peri); CHKERRQ(ierr);
    
    /* Check to make sure we get the rigth t parameter */
    ierr = EG_evaluate(edge, paramsA, eval);
    
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]))/xyzA[ii];
    }
    delta = sqrt(abs(delta));
    
    if (delta > errTolr){
      //ierr = PetscPrintf(PETSC_COMM_SELF, " ERROR :: EDGE Point A doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    } 
    
    /* Check t parameter for Point B */
    ierr = EG_evaluate(edge, paramsB, eval);
    
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]))/xyzB[ii];
    }
    delta = sqrt(abs(delta));
    
    if (delta > errTolr){
      //ierr = PetscPrintf(PETSC_COMM_SELF, " ERROR :: EDGE Point B doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    }

    /* Calculate t parameter for New Vertex at Midpoint */    
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;          // t- parameter for New Vertex at midpoint
    
    if ((paramsNew[0] >= range[0]) && (paramsNew[0] <= range[1])){
      ierr = EG_evaluate(edge, paramsNew, eval); CHKERRQ(ierr);  // (x,y,z) coordinates for New Vertex at midpoint
    } else {
      ierr = PetscPrintf(PETSC_COMM_SELF, "  ERROR :: New Vertex is outside of EDGE range || p = %d \n", p); CHKERRQ(ierr);
	  for (int ii = 0; ii < cdim; ++ii) xyz[ii] = initCoords[ii];
      ierr = EG_invEvaluate(edge, xyz, paramsNew, eval); CHKERRQ(ierr);
    }
  
  } else if (bodyID >= 0 && faceID > 0) {
    /* Snap to FACE Geometry searching along DMPlex normal */
    double   paramsA[2], paramsB[2];
    double   paramsNew[2];
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);    // Corresponding FACE ID for DMPlex edge
    
    /* Get (u,v) parameters for Vertex A & Vertex B */
    ierr = EG_invEvaluate(face, xyzA, paramsA, resultA); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex A
    ierr = EG_invEvaluate(face, xyzB, paramsB, resultB); CHKERRQ(ierr);    // Get (u,v) parameters for Vertex B
    
    /* Check if we get the same xyz Location back for Point A */
    ierr = EG_evaluate(face, paramsA, eval); CHKERRQ(ierr);

    /* Get FACE range */
    ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
    
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzA[ii]) * (eval[ii]-xyzA[ii]))/xyzA[ii];
    }
    delta = sqrt(delta);
    
    if (delta > errTolr){
      //ierr = PetscPrintf(PETSC_COMM_SELF, " ERROR :: Point A doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    }

    /* Check if we get the same xyz location back for Point B */  
    ierr = EG_evaluate(face, paramsB, eval); CHKERRQ(ierr);
    delta = 0.;
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((eval[ii]-xyzB[ii]) * (eval[ii]-xyzB[ii]))/xyzB[ii];
    }
    delta = sqrt(delta);

    if (delta > errTolr){
      //ierr = PetscPrintf(PETSC_COMM_SELF, " ERROR :: Point B doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
    }
    
    /* Get (u,v) range for FACE */
    ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
    
    /* Calculate (u, v) parameters for refined point */
    paramsNew[0] = (paramsA[0] + paramsB[0]) / 2.;
    paramsNew[1] = (paramsA[1] + paramsB[1]) / 2.;
    
    if ((paramsNew[0] >= range[0] && paramsNew[0] <= range[1]) && (paramsNew[1] >= range[2] && paramsNew[1] <= range[3])){
      ierr = EG_evaluate(face, paramsNew, eval); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_SELF, " ERROR :: FACE has to use invEvaluate() || p = %d \n", p); CHKERRQ(ierr);
	  for (int ii = 0; ii < cdim; ++ii) xyz[ii] = initCoords[ii];
      ierr = EG_invEvaluate(face, xyz, paramsNew, eval); CHKERRQ(ierr);
    }
  } else {
	/* Interior Nodes - Pass through coordinates */
	for (d = 0; d < cdim; ++d) eval[d] = initCoords[d];
  }
  for (d = 0; d < cdim; ++d) coords[d] = eval[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode DMPlexInflateToGeomModel(DM dm)
{
	/* Variables for Function */
	PetscErrorCode  	ierr;
	DMLabel         	bodyLabel, faceLabel, edgeLabel, vertexLabel;
	PetscContainer  	modelObj;
	PetscInt        	cdim, bodyID, faceID, edgeID, vertexID;
	PetscInt        	vecSize, vStart, vEnd, p;
	ego             	model, geom, body, face, edge;
	ego				   *bodies;
	int             	Nb, oclass, mtype, *senses, pntCntr;
	PetscInt        	ix[3];
	double          	xyz[3];
	PetscScalar         result[3];
	Vec             	vCoords, inflateCoords;

	PetscFunctionBegin;
	/* Get Coordinate Dimension of DMPlex */
	ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_SELF, "\n Start DMPlexInflateToGeomModel \n");CHKERRQ(ierr);
	
#ifdef PETSC_HAVE_EGADS
	ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
	ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
	ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
	ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
	
	ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
	ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);

	/* Get Vector with all DMPlex Node Coordinates */
	ierr = VecCreate(PETSC_COMM_WORLD, &vCoords); CHKERRQ(ierr);
	ierr = VecSetType(vCoords, VECSTANDARD); CHKERRQ(ierr);
	ierr = DMGetCoordinatesLocal(dm, &vCoords); CHKERRQ(ierr);    
	ierr = VecGetSize(vCoords, &vecSize); CHKERRQ(ierr);
    
	ierr = VecAssemblyBegin(vCoords); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(vCoords); CHKERRQ(ierr);
	//ierr = PetscPrintf(PETSC_COMM_SELF, "\n Set up Original Coordinate Vector \n");CHKERRQ(ierr);
	
	/* Setup Vector for Inflated Coordinates */
	ierr = VecCreate(PETSC_COMM_WORLD, &inflateCoords); CHKERRQ(ierr);
	ierr = VecSetType(inflateCoords, VECSTANDARD); CHKERRQ(ierr);
	ierr = VecSetSizes(inflateCoords, PETSC_DECIDE, vecSize); CHKERRQ(ierr);
    
	/* Get DMPlex Start ID for Vertices */
	ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); CHKERRQ(ierr); CHKERRQ(ierr);
	
	for (p = vStart; p < vEnd; ++p) {
		ierr = PetscPrintf(PETSC_COMM_SELF, "  p = %d \n", p);CHKERRQ(ierr);	
		
		ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(vertexLabel, p, &vertexID);CHKERRQ(ierr);
		
		ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyID = %d \n", bodyID);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "    faceID = %d \n", faceID);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeID = %d \n", edgeID);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "    vertexID = %d \n", vertexID);CHKERRQ(ierr);
		
		if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);		
		body = bodies[bodyID];
		
		/* Get Original Coordinates stored in the DMPlex */
		pntCntr = 0;
		for (int ii = (p-vStart)*cdim; ii < (p-vStart)*cdim+cdim; ++ii){
			ix[pntCntr] = ii;
			pntCntr = pntCntr + 1;
		}
		ierr = VecGetValues(vCoords, cdim, ix, xyz); CHKERRQ(ierr);    // Vertex coordinates (x,y,z)
		ierr = PetscPrintf(PETSC_COMM_SELF, "    xyz[x, y, z]  = [%lf, %lf, %lf] \n", xyz[0], xyz[1], xyz[2]);CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_SELF, "  Original Coordinates Retrieved \n");CHKERRQ(ierr);
				
		if (edgeID > 0) {
			/* Snap to EDGE at nearest location */
			double   params[1];
			ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge); CHKERRQ(ierr);
			ierr = EG_invEvaluate(edge, xyz, params, result); CHKERRQ(ierr);   // Get (x,y,z) of nearest point on EDGE
		} else if (faceID > 0) {
			/* Snap to FACE at nearest location */
			double   params[2];
			ierr = EG_objectBodyTopo(body, FACE, faceID, &face); CHKERRQ(ierr);    // Corresponding FACE for FaceID
			ierr = EG_invEvaluate(face, xyz, params, result); CHKERRQ(ierr);    // Get (x,y,z) of nearest point on FACE
		} else {
			for ( int ii = 0; ii < 3; ++ii) {
				result[ii] = xyz[ii];
			}
		}
		ierr = PetscPrintf(PETSC_COMM_SELF, "    result[x, y, z]  = [%lf, %lf, %lf] \n", result[0], result[1], result[2]);CHKERRQ(ierr);
		ierr = VecSetValues(inflateCoords, cdim, ix, result, INSERT_VALUES); CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_SELF, "  Inflated Coordinates Stored \n");CHKERRQ(ierr);
	}
	
	/* Prepare inflateCoords vector for DMPlex insertion */
	ierr = VecAssemblyBegin(inflateCoords); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(inflateCoords); CHKERRQ(ierr);
	
	/* Insert Inflated Coordinates into DMPlex */
	ierr = DMSetCoordinatesLocal(dm, inflateCoords); CHKERRQ(ierr);
	//ierr = PetscPrintf(PETSC_COMM_SELF, "  dmMeah Coordinates Updated \n");CHKERRQ(ierr);
	ierr = VecDestroy(&inflateCoords); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "  inflateCoords Destroyed \n");CHKERRQ(ierr);
#else
	// Do Nothing
#endif  
	PetscFunctionReturn(0);
}