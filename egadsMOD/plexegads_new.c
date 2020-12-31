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
  double paramsUV[2] = {0., 0.};
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
	  paramsUV[0] = paramsV[v*3+0];
	  paramsUV[1] = paramsV[v*3+1];
	  
	  /* Remove next line to get rid of parameter correction */
	  //Bierr = DMPlexGeomParameterCorrect(obj, dE, edgeID, faceID, paramsUV, coords); CHKERRQ(ierr);
	  /*---------------------------------------------------- */
	  
	  paramsV[v*3+0] = paramsUV[0];
	  paramsV[v*3+1] = paramsUV[1];
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

PETSC_EXTERN PetscErrorCode DMPlexGeomParameterCorrect(ego obj, PetscInt cdm, PetscInt edgeID, PetscInt faceID, PetscScalar params[], const PetscScalar coords[])
{
	/* Variables for Function */
	PetscErrorCode  	ierr;
	DMLabel         	bodyLabel, faceLabel, edgeLabel, vertexLabel;
	PetscContainer  	modelObj;
	//PetscInt        	cdim, bodyID, faceID, edgeID, vertexID;
	PetscInt        	vecSize, vStart, vEnd, p;
	ego             	model, geom, body, face, edge;
	ego				   *bodies;
	int             	Nb, oclass, mtype, *senses, pntCntr;
	PetscInt        	ix[3];
	double          	xyz[3];
	PetscScalar         result[18];
	PetscScalar         delta, errTolr = 1.0e-4;
	Vec             	vCoords, inflateCoords;

	PetscFunctionBegin;
	/* Get Coordinate Dimension of DMPlex */
	//ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);
	
	ierr = PetscPrintf(PETSC_COMM_SELF, "\n Start DMPlexGeomParameterCorrect \n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "       cdm = %d \n", cdm);CHKERRQ(ierr);
#ifdef PETSC_HAVE_EGADS

if (edgeID > 0 ) {
    ierr = EG_evaluate(obj, params, result);
	
	double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(obj, range, &peri); CHKERRQ(ierr);
    
    delta = 0.;
    for (int ii=0; ii<cdm; ++ii){
      delta = delta + ((result[ii]-coords[ii]) * (result[ii]-coords[ii]));
    }
    delta = sqrt(abs(delta));
    
    if (delta > errTolr){
		double dx, dy, dz, dt;
		double dtdx, dtdy, dtdz;
       
       for (int jj=0; jj<3; ++jj){
        //ierr = PetscPrintf(PETSC_COMM_SELF, " EDGE Point B doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  EDGE :: p = %d \n", p); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzB (x, y, z) = (%lf, %lf, %lf) \n", xyzB[0], xyzB[1], xyzB[2]);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "               (t) = (%lf) \n", paramsB[0]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (tmin, tmax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    range (vmin, vmax) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeID);CHKERRQ(ierr);
      
        dx = result[0] - coords[0];
		if (cdm > 1) { dy = result[1] - coords[1];}
        if (cdm > 2) { dz = result[2] - coords[2];}
		if (cdm > 3) { ierr = PetscPrintf(PETSC_COMM_SELF, " Error :: > 3 coordinate dimension not supportted! \n");CHKERRQ(ierr); }
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        //dv = (eval[3]*dy - eval[4]*dx)/(eval[3]*eval[7] - eval[4]*eval[6]);
        //du = (dx - eval[6]*dv)/eval[3];
        
        if (result[3] != 0.){
          dtdx = 1. / result[3];
        } else {
          dtdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdx = %lf \n", dtdx); CHKERRQ(ierr);
        
		if (cdm > 1) {
			if (result[4] != 0.){
				dtdy = 1. / result[4];
			} else {
				dtdy = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdy = %lf \n", dtdy); CHKERRQ(ierr);
		}
        
		if (cdm > 2) {
			if (result[5] != 0.){
				dtdz = 1. / result[5];
			} else {
				dtdz = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dtdz = %lf \n", dtdz); CHKERRQ(ierr);
		}
       
        //dudx = eval[3];
        //dudy = eval[4];
        //dudz = eval[5];
        
        dt = dtdx*dx;
		if (cdm > 1) { dt += dtdy*dy; }
		if (cdm > 2) { dt += dtdz*dz; }
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dt = %.6e \n", dt); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[3] = %.6e \n", eval[3]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[4] = %.6e \n", eval[4]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[5] = %.6e \n", eval[5]); CHKERRQ(ierr);
        
        if (params[0] - dt < range[0]){
          //paramsB[0] = (paramsB[0] + range[0]) / 2.;
          //paramsB[0] = paramsB[0]; // - (dt/abs(dt))*1.e-4;
          params[0] = range[0];
        } else if (params[0] - dt > range[1]) {
          //paramsB[0] = (paramsB[0] + range[1]) / 2.;
          //paramsB[0] = paramsB[0]; // - (dt/abs(dt))*1.e-4;
          params[0] = range[1];		// was range[0]
        } else {
          params[0] = params[0] - dt;
        }
        
        //paramsA[0] = paramsA[0] - du;
        //paramsA[1] = paramsA[1] - dv;
		
		//***params[0] = params[0] - dt;
		ierr = EG_evaluate(obj, params, &result); CHKERRQ(ierr);
        
        //**ierr = EG_evaluate(edge, paramsB, &eval); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii = 0; ii < cdm; ++ii){
          delta = delta + ((result[ii]-coords[ii]) * (result[ii]-coords[ii]));
        }
        delta = sqrt(delta);
       }
     }    
     // End of Check
}
else if (faceID > 0) {
    // Check if we get the same xyz Location back -- if not we have to correct
    ierr = EG_evaluate(obj, params, &result); CHKERRQ(ierr);
    //double delta;
    double dx = 0., dy = 0., dz = 0., du = 0., dv = 0.;
    double dudx = 0., dudy = 0., dudz = 0., dvdx = 0., dvdy = 0., dvdz = 0.;
    
    double range[4];    // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(obj, range, &peri); CHKERRQ(ierr);
    
    delta = 0.;  
    for (int ii=0; ii<3; ++ii){
      delta = delta + ((result[ii]-coords[ii]) * (result[ii]-coords[ii]));
    }
    delta = sqrt(delta);
    
    //int  xyzCheck = 0;
    if (delta > errTolr){
      //xyzCheck = 1;      
      
      // Trial for better Evaluation of (u,v) parameters  
      //paramsA[0] = paramsA[0] / 2.;
      //paramsA[1] = paramsA[1] / 2.;
      //ierr = EG_evaluate(face, paramsA, &eval); CHKERRQ(ierr);
      //int jj = 0;
      //while (delta > errTolr) {
      for (int jj = 0; jj < 3; ++jj){
        //ierr = PetscPrintf(PETSC_COMM_SELF, "Point A doesn't match DMPlex Coordinates \n"); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  FACE :: p = %d \n", p); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "  jj = %d \n", jj); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    delta = %.6e \n", delta); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eCone[0] = %d \n", eCone[0]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    xyzA (x, y, z) = (%lf, %lf, %lf) \n", xyzA[0], xyzA[1], xyzA[2]);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "            (u, v) = (%lf, %lf) \n", paramsA[0], paramsA[1]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval (x, y, z) = (%lf, %lf, %lf) \n", eval[0], eval[1], eval[2]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (umin, umax) = (%lf, %lf) \n", range[0], range[1]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    range (vmin, vmax) = (%lf, %lf) \n", range[2], range[3]);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    peri = %d \n", peri);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    faceId = %d \n", faceID);CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    edgeId = %d \n", edgeID);CHKERRQ(ierr);
        
        dx = result[0] - coords[0];
        if (cdm > 1) { dy = result[1] - coords[1]; }
        if (cdm > 2) { dz = result[2] - coords[2]; }
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dx = %lf \n", dx); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dy = %lf \n", dy); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dz = %lf \n", dz); CHKERRQ(ierr);
        
        //dv = (eval[3]*dy - eval[4]*dx)/(eval[3]*eval[7] - eval[4]*eval[6]);
        //du = (dx - eval[6]*dv)/eval[3];
        
        if (result[3] != 0.){
          dudx = 1. / result[3];
        } else {
          dudx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dudx = %lf \n", dudx); CHKERRQ(ierr);
		
		if (result[6] != 0.){
          dvdx = 1. / result[6];
        } else {
          dvdx = 0.;
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdx = %lf \n", dvdx); CHKERRQ(ierr);
        
		if ( cdm > 1) {
			if (result[4] != 0.){
				dudy = 1. / result[4];
			} else {
				dudy = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dudy = %lf \n", dudy); CHKERRQ(ierr);

			if (result[7] != 0.){
				dvdy = 1. / result[7];
			} else {
				dvdy = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdy = %lf \n", dvdy); CHKERRQ(ierr);
		}
        
		if ( cdm > 2) {
			if (result[5] != 0.){
				dudz = 1. / result[5];
			} else {
				dudz = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dudz = %lf \n", dudz); CHKERRQ(ierr);

			if (result[8] != 0.){
				dvdz = 1. / result[8];
			} else {
				dvdz = 0.;
			}
			ierr = PetscPrintf(PETSC_COMM_SELF, "    dvdz = %lf \n", dvdz); CHKERRQ(ierr);
		}
       
        //dudx = eval[3];
        //dudy = eval[4];
        //dudz = eval[5];
        //dvdx = eval[6];
        //dvdy = eval[7];
        //dvdz = eval[8];
        
        du = dudx*dx;
        dv = dvdx*dx;
		if ( cdm > 1) { du += dudy*dy; dv += dvdy*dy; }
		if ( cdm > 2) { du += dudz*dz; dv += dvdz*dz; }
        
        ierr = PetscPrintf(PETSC_COMM_SELF, "    du = %.6e \n", du); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_SELF, "    dv = %.6e \n", dv); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[3] = %.6e \n", eval[3]); CHKERRQ(ierr);
        //ierr = PetscPrintf(PETSC_COMM_SELF, "    eval[4] = %.6e \n", eval[4]); CHKERRQ(ierr);
        
		// --- Modified these calculation ---
        if (params[0] - du < range[0]){
          //paramsA[0] = (paramsA[0] + range[0]) / 2.;
          //paramsA[0] = paramsA[0]; // - (du/abs(du))*1.e-4;
          params[0] = range[0];
        } else if (params[0] - du > range[1]) {
          //paramsA[0] = (paramsA[0] + range[1]) / 2.;
          //paramsA[0] = paramsA[0]; // - (du/abs(du))*1.e-4;
          params[0] = range[1];		// was range[0]
        } else {
          params[0] = params[0] - du;
        }
        
        if (params[1] - dv < range[2]){
          //paramsA[1] = (paramsA[1] + range[2]) / 2.;
          //paramsA[1] = paramsA[1]; // - (dv/abs(dv))*1.e-4;
          params[1] = range[2];
        } else if (params[1] - dv > range[3]) {
          //paramsA[1] = (paramsA[1] + range[3]) / 2.;
          //paramsA[1] = paramsA[1]; // - (dv/abs(dv))*1.e-4;
          params[1] = range[3];
        } else {
          params[1] = params[1] - dv;
        }
        //----------------------
		
        //params[0] = params[0] - du;
        //params[1] = params[1] - dv;
        
        ierr = EG_evaluate(obj, params, &result); CHKERRQ(ierr);
        
        delta = 0.;
        for (int ii = 0; ii < cdm; ++ii){
          delta = delta + ((result[ii]-coords[ii]) * (result[ii]-coords[ii]));
        }
        delta = sqrt(delta);
        
        //jj = ++jj;
      } // End of Trial
    }
}
else {
	// Do Nothing
}

#else
	// Do Nothing
#endif  
	PetscFunctionReturn(0);
}
