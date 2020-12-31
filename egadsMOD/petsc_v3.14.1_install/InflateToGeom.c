
/*@C
  DMPlexInflateToGeomModel - Inflates the cTetgen/Tetgen generated Steiner points of 3D Plex to the EGAD CAD model.
                             This is useful when the initial PLC is not retained by cTetgen/Tetgen.
                             cTetgen/Tetgen only understands the PLC and does not utilize the DM Plex's embedded EGADS data.
                             

  Collective

  Input Parameters:
+ dm     - Uninflated DM Plex with EGADS Embedded Model

  Output Parameter:
. dm     - Updated DM with points on EGADS CAD Geometry

  Level: beginner

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
PetscErrorCode DMPlexInflateToGeomModel(DM dm)
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
	
	/* Setup Vector for Inflated Coordinates */
	ierr = VecCreate(PETSC_COMM_WORLD, &inflateCoords); CHKERRQ(ierr);
	ierr = VecSetType(inflateCoords, VECSTANDARD); CHKERRQ(ierr);
	ierr = VecSetSizes(inflateCoords, PETSC_DECIDE, vecSize); CHKERRQ(ierr);
    
	/* Get DMPlex Start ID for Vertices */
	ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd); CHKERRQ(ierr); CHKERRQ(ierr);
	
	for (p = vStart; p < vEnd; ++p) {
		ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(vertexLabel, p, &vertexID);CHKERRQ(ierr);
		
		if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);		
		body = bodies[bodyID];
		
		/* Get Original Coordinates stored in the DMPlex */
		pntCntr = 0;
		for (int ii = (p-vStart)*cdim; ii < (p-vStart)*cdim+cdim; ++ii){
			ix[pntCntr] = ii;
			pntCntr = pntCntr + 1;
		}
		ierr = VecGetValues(vCoords, cdim, ix, xyz); CHKERRQ(ierr);
				
		if (edgeID > 0) {
			/* Snap to EDGE at nearest location */
			double   params[1];
			ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge); CHKERRQ(ierr);
			ierr = EG_invEvaluate(edge, xyz, params, result); CHKERRQ(ierr);
		} else if (faceID > 0) {
			/* Snap to FACE at nearest location */
			double   params[2];
			ierr = EG_objectBodyTopo(body, FACE, faceID, &face); CHKERRQ(ierr);
			ierr = EG_invEvaluate(face, xyz, params, result); CHKERRQ(ierr);
		} else { 
		  for ( int ii = 0; ii < 3; ++ii) {
		    result[ii] = xyz[ii];}
		  }
		ierr = VecSetValues(inflateCoords, cdim, ix, result, INSERT_VALUES); CHKERRQ(ierr);
	}
	
	/* Prepare inflateCoords vector for DMPlex insertion */
	ierr = VecAssemblyBegin(inflateCoords); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(inflateCoords); CHKERRQ(ierr);
	
	/* Insert Inflated Coordinates into DMPlex */
	ierr = DMSetCoordinatesLocal(dm, inflateCoords); CHKERRQ(ierr);
	ierr = VecDestroy(&inflateCoords); CHKERRQ(ierr);
#else
	// Do Nothing
#endif  
	PetscFunctionReturn(0);
}

