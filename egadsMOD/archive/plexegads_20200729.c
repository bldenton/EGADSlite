#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

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
PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar mcoords[], PetscScalar gcoords[])
{
#ifdef PETSC_HAVE_EGADS
  DM             cdm;
  DMLabel        bodyLabel, faceLabel, edgeLabel;
  PetscContainer modelObj;
  PetscInt       bodyID, faceID, edgeID;
  ego           *bodies;
  ego            model, geom, body, obj;
  /* result has to hold derviatives, along with the value */
  double         params[3], result[18], paramsV[16*3], resultV[16*3];
  int            Nb, oclass, mtype, *senses;
  Vec            coordinatesLocal;
  PetscScalar   *coords = NULL;
  PetscInt       Nv, v, Np = 0, pm;
  PetscReal      PARAMACC = 1.0e-4;         /* parameter accuracy */
#endif
  PetscInt       dE, d;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &dE);CHKERRQ(ierr);
#ifdef PETSC_HAVE_EGADS
  ierr = DMGetLabel(dm, "EGADS Body ID",   &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID",   &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID",   &edgeLabel);CHKERRQ(ierr);
  if (!bodyLabel || !faceLabel || !edgeLabel) {
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinatesLocal);CHKERRQ(ierr);
  ierr = DMLabelGetValue(bodyLabel,   p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel,   p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel,   p, &edgeID);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  
  if (bodyID < 0) {PetscFunctionReturn(0);}  // This accounts for "connective" edges if there are multiple bodies
  //if (bodyID < 0 || bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];
  
  //ierr = PetscPrintf(PETSC_COMM_SELF, "  (p, bodyID, EdgeID, FaceID) = ( %d, %d, %d, %d) \n", p, bodyID, edgeID, faceID);CHKERRQ(ierr);
  
  if (edgeID >= 0)      {ierr = EG_objectBodyTopo(body, EDGE, edgeID, &obj);CHKERRQ(ierr); Np = 1;}
  else if (faceID >= 0) {ierr = EG_objectBodyTopo(body, FACE, faceID, &obj);CHKERRQ(ierr); Np = 2;}
  //else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Point %D is not in edge or face label for EGADS", p);
  else {for (d = 0; d <dE; ++d) gcoords[d] = mcoords[d]; PetscFunctionReturn(0); }
  /* Calculate parameters (t or u,v) for vertices */
  ierr = DMPlexVecGetClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
  Nv  /= dE;
  if (Nv == 1) {
    ierr = DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
    for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
    PetscFunctionReturn(0);
  }
  if (Nv > 16) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Cannot handle %D coordinates associated to point %D", Nv, p);
  // NEW CODE
  double initCoordsA[3], initCoordsB[3]; //, initCoordsAA[3], initCoordsBB[3];
  double initUVa[2] = {0., 0.}, initUVb[2] = {0., 0.};
  
  double range[4]; // Edge [tmin, tmax 0 0 ] :: Face [umin, umax, vmin, vmax]
  double initCoords[3];
  int    peri;
  ierr = EG_getRange(obj, range, &peri); CHKERRQ(ierr);
  // END NEW CODE
  for (v = 0; v < Nv; ++v) {
	  ierr = EG_invEvaluate(obj, &coords[v*dE], &paramsV[v*3], &resultV[v*3]);CHKERRQ(ierr);
	  // NEW CODE
	  if ( v == 0) { 
		initCoordsA[0] = coords[v*dE+0];
		initCoordsA[1] = coords[v*dE+1];
		initCoordsA[2] = coords[v*dE+2];
		initUVa[0] = paramsV[v*3+0];
		initUVa[1] = paramsV[v*3+1];
		//ierr = EG_evaluate(obj, initUVa, initCoordsAA);CHKERRQ(ierr);
	  }
	  if ( v == 1) { 
		initCoordsB[0] = coords[v*dE+0];
		initCoordsB[1] = coords[v*dE+1];
		initCoordsB[2] = coords[v*dE+2];
		initUVb[0] = paramsV[v*3+0];
		initUVb[1] = paramsV[v*3+1];
		//ierr = EG_evaluate(obj, initUVb, initCoordsBB);CHKERRQ(ierr);
	  }
	  
	  if (peri > 0){
		// u-parameter fix
		if (paramsV[v*3+0] + PARAMACC < range[0]){
			paramsV[v*3+0] += 2. * PETSC_PI;
		} else if (paramsV[v*3+0] - PARAMACC > range[1]) {
			paramsV[v*3+0] -= 2. * PETSC_PI;
		} else {
			// Do Nothing
		}
		
		if (peri > 1 ) {
			// v-parameter fix
			if (paramsV[v*3+1] + PARAMACC < range[2]){
				paramsV[v*3+1] += 2. * PETSC_PI;
			} else if (paramsV[v*3+1] - PARAMACC > range[3]) {
				paramsV[v*3+1] -= 2. * PETSC_PI;
			} else {
				// Do Nothing
			} 
		}
	  }
	  
	  // END NEW CODE
  }
  ierr = DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
  /* Calculate parameters (t or u,v) for new vertex at edge midpoint */
  for (pm = 0; pm < Np; ++pm) {
    params[pm] = 0.;
    for (v = 0; v < Nv; ++v) {params[pm] += paramsV[v*3+pm];}
    params[pm] /= Nv;
  }
  // TODO Check
    //Adouble range[4]; // Edge [tmin, tmax 0 0 ] :: Face [umin, umax, vmin, vmax]
	//Adouble initCoords[3];
    //Aint    peri;
    //Aierr = EG_getRange(obj, range, &peri); CHKERRQ(ierr);
	
	if (edgeID > 0 ) {
		if ((params[0] < range[0]) || (params[0] > range[1])) {
			//for (d = 0; d < dE; ++d) {
			//	initCoords[d]= 0.;
			//	for (v = 0; v < Nv; ++v) { initCoords[d] += coords[v*dE+d]; }
			//	coords[d] /= Nv;
			//}
			//ierr = EG_invEvaluate(obj, initCoords, params, result); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "WARNING :: invEvaluate() used on EGADS EDGE %d for DMPlex Node %d \n", edgeID, p); CHKERRQ(ierr);
		} else {
			//ierr = EG_evaluate(obj, params, result);CHKERRQ(ierr);
		}
	} else if (faceID > 0) {
		if ((params[0] < range[0]) || (params[0] > range[1]) || (params[1] < range[2]) || (params[1] > range[3])) {
			//for (d = 0; d < dE; ++d) {
			//	initCoords[d]= 0.;
			//	for (v = 0; v < Nv; ++v) { initCoords[d] += coords[v*dE+d]; }
			//	coords[d] /= Nv;
			//}
			//ierr = EG_invEvaluate(obj, initCoords, params, result); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "WARNING :: invEvaluate() used on EGADS FACE %d for DMPlex Node %d \n", faceID, p); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "              initCoordsA = [%lf, %lf, %lf] \n", initCoordsA[0], initCoordsA[1], initCoordsA[2]); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "                  initUVa = [%lf, %lf] \n", initUVa[0], initUVa[1]); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "              initCoordsB = [%lf, %lf, %lf] \n", initCoordsB[0], initCoordsB[1], initCoordsB[2]); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "                  initUVb = [%lf, %lf] \n", initUVb[0], initUVb[1]); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "                   params = [%lf, %lf] \n", params[0], params[1]); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "                    Range = [%lf, %lf, %lf, %lf] \n", range[0], range[1], range[2], range[3]); CHKERRQ(ierr);
			//ierr = PetscPrintf(PETSC_COMM_SELF, "           initCoordsAA = [%lf, %lf, %lf] \n", initCoordsAA[0], initCoordsAA[1], initCoordsAA[2]); CHKERRQ(ierr);
			//ierr = PetscPrintf(PETSC_COMM_SELF, "           initCoordsBB = [%lf, %lf, %lf] \n", initCoordsBB[0], initCoordsBB[1], initCoordsBB[2]); CHKERRQ(ierr);
			
			//
			/*
			if (peri > 0){
				// u-parameter fix
				if (params[0] + PARAMACC < range[0]){
					params[0] += 2. * PETSC_PI;
				} else if (params[0] - PARAMACC > range[1]) {
					params[0] -= 2. * PETSC_PI;
				} else {
					// Do Nothing
				}
				
				if (peri > 1 ) {
				// v-parameter fix
					if (params[1] + PARAMACC < range[2]){
						params[1] += 2. * PETSC_PI;
					} else if (params[1] - PARAMACC > range[3]) {
						params[1] -= 2. * PETSC_PI;
					} else {
						// Do Nothing
					} 
				}
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "                                  paramsMOD = [%lf, %lf] \n", params[0], params[1]); CHKERRQ(ierr);
			*/
			//
		} else {
			//ierr = EG_evaluate(obj, params, result);CHKERRQ(ierr);
		}
	} else {
		// Do Nothing
	}
  //---------------------------------------
  /* Put coordinates for new vertex in result[] */
  ierr = EG_evaluate(obj, params, result);CHKERRQ(ierr);
  for (d = 0; d < dE; ++d) gcoords[d] = result[d];
#else
  for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
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
		//ierr = PetscPrintf(PETSC_COMM_SELF, "  p = %d \n", p);CHKERRQ(ierr);	
		
		ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
		ierr = DMLabelGetValue(vertexLabel, p, &vertexID);CHKERRQ(ierr);
		
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyID = %d \n", bodyID);CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    faceID = %d \n", faceID);CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeID = %d \n", edgeID);CHKERRQ(ierr);
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    vertexID = %d \n", vertexID);CHKERRQ(ierr);
		
		if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);		
		body = bodies[bodyID];
		
		/* Get Original Coordinates stored in the DMPlex */
		pntCntr = 0;
		for (int ii = (p-vStart)*cdim; ii < (p-vStart)*cdim+cdim; ++ii){
			ix[pntCntr] = ii;
			pntCntr = pntCntr + 1;
		}
		ierr = VecGetValues(vCoords, cdim, ix, xyz); CHKERRQ(ierr);    // Vertex coordinates (x,y,z)
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    xyz[x, y, z]  = [%lf, %lf, %lf] \n", xyz[0], xyz[1], xyz[2]);CHKERRQ(ierr);
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
		//ierr = PetscPrintf(PETSC_COMM_SELF, "    result[x, y, z]  = [%lf, %lf, %lf] \n", result[0], result[1], result[2]);CHKERRQ(ierr);
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