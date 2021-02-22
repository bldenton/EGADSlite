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
  double         params[3], result[18], paramsV[16*3], resultV[16*3], range[4];
  int            Nb, oclass, mtype, *senses;
  Vec            coordinatesLocal;
  PetscScalar   *coords = NULL;
  PetscReal      PARAMACC = 1.0e-4;
  PetscInt       Nv, v, Np = 0, pm, peri;
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
  
  /* Allows for "Connective" Plex Edges present in models with multiple non-touching Entities */
  if (bodyID < 0) {PetscFunctionReturn(0);}
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];

  if (edgeID >= 0)      {ierr = EG_objectBodyTopo(body, EDGE, edgeID, &obj);CHKERRQ(ierr); Np = 1;}
  else if (faceID >= 0) {ierr = EG_objectBodyTopo(body, FACE, faceID, &obj);CHKERRQ(ierr); Np = 2;}
  else {for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d]; PetscFunctionReturn(0); }

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
  ierr = EG_getRange(obj, range, &peri); CHKERRQ(ierr);
  for (v = 0; v < Nv; ++v) {
	  ierr = EG_invEvaluate(obj, &coords[v*dE], &paramsV[v*3], &resultV[v*3]);CHKERRQ(ierr);
	  if (peri > 0) {
		  if (paramsV[v*3+0] + PARAMACC < range[0]) { paramsV[v*3+0] += 2. * PETSC_PI; }
		  else if (paramsV[v*3+0] - PARAMACC > range[1]) { paramsV[v*3+0] -= 2. * PETSC_PI; }
		  else { /* Do Nothing */}
	  }
	  if (peri > 1) {
		  if (paramsV[v*3+1] + PARAMACC < range[2]) { paramsV[v*3+1] += 2. * PETSC_PI; }
		  else if (paramsV[v*3+1] - PARAMACC > range[3]) { paramsV[v*3+1] -= 2. * PETSC_PI; }
		  else { /* Do Nothing */}
	  }	  
  }
  ierr = DMPlexVecRestoreClosure(cdm, NULL, coordinatesLocal, p, &Nv, &coords);CHKERRQ(ierr);
  /* Calculate parameters (t or u,v) for new vertex at edge midpoint */
  for (pm = 0; pm < Np; ++pm) {
    params[pm] = 0.;
    for (v = 0; v < Nv; ++v) {params[pm] += paramsV[v*3+pm];}
    params[pm] /= Nv;
  }
  /* TODO Check -- MAY NOT NEED BECAUSE OF PARAMETER CHECK ABOVE
    double range[4]; // [umin, umax, vmin, vmax]
    int    peri;
    ierr = EG_getRange(face, range, &peri);CHKERRQ(ierr);
    if ((paramsNew[0] < range[0]) || (paramsNew[0] > range[1]) || (paramsNew[1] < range[2]) || (paramsNew[1] > range[3])) SETERRQ();
  */
  /* Put coordinates for new vertex in result[] */
  ierr = EG_evaluate(obj, params, result);CHKERRQ(ierr);
  for (d = 0; d < dE; ++d) gcoords[d] = result[d];
#else
  for (d = 0; d < dE; ++d) gcoords[d] = mcoords[d];
#endif
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_EGADS)
static PetscErrorCode DMPlexEGADSPrintModel(ego model)
{
  ego            geom, *bodies, *nobjs, *mobjs, *lobjs, *shobjs, *fobjs, *eobjs;
  int            oclass, mtype, *senses, *shsenses, *fsenses, *lsenses, *esenses;
  int            Nb, b;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* test bodyTopo functions */
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", Nb);CHKERRQ(ierr);

  for (b = 0; b < Nb; ++b) {
    ego body = bodies[b];
    int id, sh, Nsh, f, Nf, l, Nl, e, Ne, v, Nv;

    /* Output Basic Model Topology */
    ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

    ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);

    ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

    ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &eobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of EDGES: %d \n", Ne);CHKERRQ(ierr);

    ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &nobjs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);
	
	EG_free(shobjs);
	EG_free(fobjs);
	EG_free(lobjs);
	EG_free(eobjs);
	EG_free(nobjs);
	
	/* List Topology of Bodies */
	ierr = PetscPrintf(PETSC_COMM_SELF, "\n"); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "   Body %d Topology :: \n", b);CHKERRQ(ierr);

    /* Get SHELL info which associated with the current BODY */
    ierr = EG_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses);CHKERRQ(ierr);
	
	for (sh = 0; sh < Nsh; ++sh) {
	  ego shell   = shobjs[sh];
	  int shsense = shsenses[sh];
	  
	  id   = EG_indexBodyTopo(body, shell);
	  ierr = PetscPrintf(PETSC_COMM_SELF, "      SHELL ID: %d :: sense = %d\n", id, shsense);CHKERRQ(ierr);
	  
	  /* Get FACE infor associated with current SHELL */
	  ierr = EG_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses);CHKERRQ(ierr);
	  
	  for (f = 0; f < Nf; ++f) {
		ego face   = fobjs[f];

        id   = EG_indexBodyTopo(body, face);
        ierr = PetscPrintf(PETSC_COMM_SELF, "        FACE ID: %d \n", id);CHKERRQ(ierr);
          
		/* Get LOOP info associated with current FACE */
        ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lsenses);CHKERRQ(ierr);
        
        for (l = 0; l < Nl; ++l) {
          ego loop   = lobjs[l];
		  int lsense = lsenses[l];

          id   = EG_indexBodyTopo(body, loop);
          ierr = PetscPrintf(PETSC_COMM_SELF, "          LOOP ID: %d :: sense = %d\n", id, lsense);CHKERRQ(ierr);

          /* Get EDGE info associated with the current LOOP */
          ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &esenses);CHKERRQ(ierr);
          for (e = 0; e < Ne; ++e) {
            ego    edge      = eobjs[e];
			ego    topRef, prev, next;
			int    esense    = esenses[e];
            double range[4]  = {0., 0., 0., 0.};
            int    peri;

		    id   = EG_indexBodyTopo(body, edge); //CHKERRQ(ierr);
			ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "            EDGE ID: %d :: sense = %d\n", id, esense); CHKERRQ(ierr);
                  
			if (mtype == DEGENERATE) { ierr = PetscPrintf(PETSC_COMM_SELF, "              EDGE %d is DEGENERATE \n", id);CHKERRQ(ierr); }
            ierr = EG_getRange(edge, range, &peri);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_SELF, "              Peri = %d :: Range = %lf, %lf, %lf, %lf \n", peri, range[0], range[1], range[2], range[3]);

            /* Get NODE info associated with the current EDGE */
            ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses);CHKERRQ(ierr);
                  
            for (v = 0; v < Nv; ++v) {
              ego    vertex = nobjs[v];
              double limits[4];
              int    dummy;
    
              ierr = EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses);CHKERRQ(ierr);
              id = EG_indexBodyTopo(body, vertex);
              ierr = PetscPrintf(PETSC_COMM_SELF, "              NODE ID: %d \n", id);CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_SELF, "                 (x, y, z) = (%lf, %lf, %lf) \n", limits[0], limits[1], limits[2]);
            }
		  }
		}
	  }
    }
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n\n");
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexEGADSDestroy_Private(void *context)
{
  if (context) EG_close((ego) context);
  return 0;
}

/*@C
  DMPlexCreateEGADS - Creates a DMPlex object representing a triangular surface mesh of the geometry and attaches a PetscContainer with the EGADS/EGADSlite geometry model and context. The DMPlex object includes DMLabels referencing the relevant DMPlex nodes to Geometry data. 

  Collective

  Input Parameters:
  comm         - Petsc MPI_Comm Object
  ego context  - EGADS/EGADSlite context object
  ego model    - EGADS/EGADSlite model object

  Output Parameter:
. newdm        - The DM object with triangular surface mesh, model, context and appropriate labels.

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
static PetscErrorCode DMPlexCreateEGADS(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel         bodyLabel, faceLabel, edgeLabel, vertexLabel;
  // EGADS/EGADSLite variables 
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
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) { 
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex
    //  Get all Nodes in model, record coordinates in a correctly formatted array
    //  Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
    //  We need to uniformly refine the initial geometry to guarantee a valid mesh
	
	// Caluculate cell and vertex sizes 
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses); CHKERRQ(ierr);
	
    ierr = PetscHMapICreate(&edgeMap); CHKERRQ(ierr);	
	ierr = PetscHMapICreate(&bodyIndexMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyVertexMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyEdgeMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyEdgeGlobalMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyFaceMap); CHKERRQ(ierr);
	
	for (b = 0; b < nbodies; ++b) {
      ego             body = bodies[b];
	  int             Nf, Ne, Nv;
	  PetscHashIter   BIiter, BViter, BEiter, BEGiter, BFiter, EMiter;
	  PetscBool       BIfound, BVfound, BEfound, BEGfound, BFfound, EMfound;
	  
	  ierr = PetscHMapIFind(bodyIndexMap, b, &BIiter, &BIfound); CHKERRQ(ierr);
	  ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound); CHKERRQ(ierr);
	  ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound); CHKERRQ(ierr);
	  ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound); CHKERRQ(ierr);
	  ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound); CHKERRQ(ierr);
	  
	  if (!BIfound)  {ierr = PetscHMapISet(bodyIndexMap, b, numFaces + numEdges + numVertices); CHKERRQ(ierr);}
	  if (!BVfound)  {ierr = PetscHMapISet(bodyVertexMap, b, numVertices); CHKERRQ(ierr);}
	  if (!BEfound)  {ierr = PetscHMapISet(bodyEdgeMap, b, numEdges); CHKERRQ(ierr);}
	  if (!BEGfound) {ierr = PetscHMapISet(bodyEdgeGlobalMap, b, edgeCntr); CHKERRQ(ierr);}
	  if (!BFfound)  {ierr = PetscHMapISet(bodyFaceMap, b, numFaces); CHKERRQ(ierr);}
	  	  
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);	  
	  ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  ierr = EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs); CHKERRQ(ierr);
	  EG_free(fobjs);
	  EG_free(eobjs);
	  EG_free(nobjs);
	  
	  // Remove DEGENERATE EDGES from Edge count 
	  ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  int Netemp = 0;
	  for (int e = 0; e < Ne; ++e) {
		  ego     edge = eobjs[e];
		  int     eid;
		  
		  ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		  eid = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);
		  
		  ierr = PetscHMapIFind(edgeMap, edgeCntr + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);		  
		  if (mtype == DEGENERATE) {
		    if (!EMfound) {ierr = PetscHMapISet(edgeMap, edgeCntr + eid - 1, -1); CHKERRQ(ierr);}
		  }
		  else {
			++Netemp;
		    if (!EMfound) {ierr = PetscHMapISet(edgeMap, edgeCntr + eid - 1, Netemp); CHKERRQ(ierr);}
		  }
	  }
	  EG_free(eobjs);
	  
	  // Determine Number of Cells 
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);
	  for (int f = 0; f < Nf; ++f) {
        ego     face = fobjs[f];
		int     edgeTemp = 0;
	  			
	  	ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);			
	  	for (int e = 0; e < Ne; ++e) {
	  	  ego     edge = eobjs[e];
	  	  
	  	  ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
	  	  if (mtype != DEGENERATE) {++edgeTemp;}
	  	}
	  	numCells += (2 * edgeTemp);
	  	EG_free(eobjs);
	  }
	  EG_free(fobjs);
	  
	  numFaces    += Nf;
	  numEdges    += Netemp;
	  numVertices += Nv;
	  edgeCntr    += Ne;
	}
	
	// Set up basic DMPlex parameters 
	dim        = 2;		// Assumes 3D Models :: Need to handle 2D modles in the future
	cdim       = 3;     // Assumes 3D Models :: Need to update to handle 2D modles in future
	numCorners = 3;     // Split Faces into triangles 
    numPoints  = numVertices + numEdges + numFaces;   // total number of coordinate points

	ierr = PetscMalloc2(numPoints*cdim, &coords, numCells*numCorners, &cells); CHKERRQ(ierr);
	
	// Get Vertex Coordinates and Set up Cells 
	for (b = 0; b < nbodies; ++b) {
	  ego             body = bodies[b];
	  int             Nf, Ne, Nv; 
	  PetscInt        bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
	  PetscHashIter   BViter, BEiter, BEGiter, BFiter, EMiter;
	  PetscBool       BVfound, BEfound, BEGfound, BFfound, EMfound;
	  
	  // Vertices on Current Body 
	  ierr = EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound); CHKERRQ(ierr);
	  if (!BVfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);}	  
	  ierr = PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart); CHKERRQ(ierr);
	  
	  for (int v = 0; v < Nv; ++v) {
	    ego    vertex = nobjs[v];
		double limits[4];
		int    id, dummy;
		
		ierr = EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses); CHKERRQ(ierr);
		id = EG_indexBodyTopo(body, vertex); CHKERRQ(ierr);

		coords[(bodyVertexIndexStart + id - 1)*cdim + 0] = limits[0];
		coords[(bodyVertexIndexStart + id - 1)*cdim + 1] = limits[1];
		coords[(bodyVertexIndexStart + id - 1)*cdim + 2] = limits[2];
	  }
	  EG_free(nobjs);
	  
	  // Edge Midpoint Vertices on Current Body 
	  ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound); CHKERRQ(ierr);
	  if (!BEfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);}	  
	  ierr = PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound); CHKERRQ(ierr);
	  if (!BEGfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);}	  
	  ierr = PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart); CHKERRQ(ierr);
	  
	  for (int e = 0; e < Ne; ++e) {
	    ego          edge = eobjs[e];
		double       range[2], avgt[1], cntrPnt[9];
		int          eid, eOffset;
		int          periodic;
		
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {continue;}
		
		eid = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);
		
		// get relative offset from globalEdgeID Vector 
		ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);
	    if (!EMfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d not found in edgeMap", bodyEdgeGlobalIndexStart + eid - 1);}	  
	    ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset); CHKERRQ(ierr);
		
		ierr = EG_getRange(edge, range, &periodic); CHKERRQ(ierr);
		avgt[0] = (range[0] + range[1]) /  2.;
		
		ierr = EG_evaluate(edge, avgt, cntrPnt); CHKERRQ(ierr);
		coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 1] = cntrPnt[1];
		coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 2] = cntrPnt[2];
	  }
	  EG_free(eobjs);

	  // Face Midpoint Vertices on Current Body 
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound); CHKERRQ(ierr);
	  if (!BFfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b); 
	  ierr = PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart); CHKERRQ(ierr);
	  
	  for (int f = 0; f < Nf; ++f) {
		ego       face = fobjs[f];
		double    range[4], avgUV[2], cntrPnt[18];
		int       peri, id;
		
		id = EG_indexBodyTopo(body, face); 
		ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
		
		avgUV[0] = (range[0] + range[1]) / 2.;
		avgUV[1] = (range[2] + range[3]) / 2.;
		ierr = EG_evaluate(face, avgUV, cntrPnt); CHKERRQ(ierr);

		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 0] = cntrPnt[0];
		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 1] = cntrPnt[1];
		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 2] = cntrPnt[2];
	  }
	  EG_free(fobjs);
	  
	  // Define Cells :: Note - This could be incorporated in the Face Midpoint Vertices Loop but was kept separate for clarity 
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);
	  for (int f = 0; f < Nf; ++f) {
		ego      face = fobjs[f];
		int      fID, midFaceID, midPntID, startID, endID, Nl;
		
		fID = EG_indexBodyTopo(body, face); CHKERRQ(ierr);
		midFaceID = numVertices + numEdges + bodyFaceIndexStart + fID - 1;
		// Must Traverse Loop to ensure we have all necessary information like the sense (+/- 1) of the edges. 
		// TODO :: Only handles single loop faces (No holes). The choices for handling multiloop faces are:    
		//            1) Use the DMPlexCreateEGADSFromFile() with the -dm_plex_egads_with_tess = 1 option.     
		//               This will use a default EGADS tessellation as an initial surface mesh.                
		//            2) Create the initial surface mesh via a 2D mesher :: Currently not availble (?future?)  
		//               May I suggest the XXXX as a starting point?                                           
		
		ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &lobjs, &lSenses); CHKERRQ(ierr);
		
	    if (Nl > 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face has %d Loops. Can only handle Faces with 1 Loop. Please use --dm_plex_egads_with_tess = 1 Option", Nl);
		for (int l = 0; l < Nl; ++l) {
          ego      loop = lobjs[l];
		  
          ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses); CHKERRQ(ierr);
		  for (int e = 0; e < Ne; ++e) {
		    ego     edge = eobjs[e];
		    int     eid, eOffset;
			
		    ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
			eid = EG_indexBodyTopo(body, edge);
		    if (mtype == DEGENERATE) { continue; }
			
		    // get relative offset from globalEdgeID Vector 
		    ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);
	        if (!EMfound) {SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);}	  
	        ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset); CHKERRQ(ierr);
			
			midPntID = numVertices + bodyEdgeIndexStart + eOffset - 1;
		  
		    ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses); CHKERRQ(ierr);
			
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
  CHKMEMQ;
  // Generate DMPlex   
  ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm); CHKERRQ(ierr);
  ierr = PetscFree2(coords, cells); CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Cells    = %D \n", numCells); CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Vertices = %D \n", numVertices); CHKERRQ(ierr);
  
  // Embed EGADS model in DM 
  {
    PetscContainer modelObj, contextObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dm, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &contextObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(contextObj, context);CHKERRQ(ierr);
    ierr = PetscContainerSetUserDestroy(contextObj, DMPlexEGADSDestroy_Private);CHKERRQ(ierr);
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
	  
	ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound); CHKERRQ(ierr);
	if (!BVfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);}	  
	ierr = PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart); CHKERRQ(ierr);
	
	ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound); CHKERRQ(ierr);
	if (!BEfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);}	  
	ierr = PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart); CHKERRQ(ierr);
	
    ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound); CHKERRQ(ierr);
	if (!BFfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b);}	  
	ierr = PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart); CHKERRQ(ierr);
	
    ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound); CHKERRQ(ierr);
    if (!BEGfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);}	  
    ierr = PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart); CHKERRQ(ierr);
	
	ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs); CHKERRQ(ierr);
	for (int f = 0; f < Nf; ++f) {
	  ego   face = fobjs[f];
      int   fID, Nl;
	  
	  fID  = EG_indexBodyTopo(body, face); CHKERRQ(ierr);
	  
	  ierr = EG_getBodyTopos(body, face, LOOP, &Nl, &lobjs);CHKERRQ(ierr);
	  for (int l = 0; l < Nl; ++l) {
        ego  loop = lobjs[l];
		int  lid;
		
		lid  = EG_indexBodyTopo(body, loop); CHKERRQ(ierr);
	    if (Nl > 1) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_SUP, "Loop %d has %d > 1 faces, which is not supported", lid, Nf);CHKERRQ(ierr);
	  
		ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses);CHKERRQ(ierr);
		for (int e = 0; e < Ne; ++e) {
		  ego     edge = eobjs[e];
		  int     eid, eOffset;
		
		  // Skip DEGENERATE Edges
		  ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		  if (mtype == DEGENERATE) {continue;}
		  eid = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);
		  
		  // get relative offset from globalEdgeID Vector
		  ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);
		  if (!EMfound) {SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);}	  
		  ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset); CHKERRQ(ierr);
				
		  ierr = EG_getBodyTopos(body, edge, NODE, &Nv, &nobjs); CHKERRQ(ierr);		
		  for (int v = 0; v < Nv; ++v){
		    ego vertex = nobjs[v];
		    int vID;
		  
		    vID = EG_indexBodyTopo(body, vertex);CHKERRQ(ierr);
		    ierr = DMLabelSetValue(bodyLabel, nStart + bodyVertexIndexStart + vID - 1, b); CHKERRQ(ierr);
		    ierr = DMLabelSetValue(vertexLabel, nStart + bodyVertexIndexStart + vID - 1, vID); CHKERRQ(ierr);
		  }
		  EG_free(nobjs);

		  ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, b); CHKERRQ(ierr);
		  ierr = DMLabelSetValue(edgeLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, eid); CHKERRQ(ierr);
			
		  // Define Cell faces
		  for (int jj = 0; jj < 2; ++jj){
		    ierr = DMLabelSetValue(bodyLabel, cellCntr, b); CHKERRQ(ierr);
		    ierr = DMLabelSetValue(faceLabel, cellCntr, fID);CHKERRQ(ierr);
		    ierr = DMPlexGetCone(dm, cellCntr, &cone); CHKERRQ(ierr);
		 
		    ierr = DMLabelSetValue(bodyLabel, cone[0], b); CHKERRQ(ierr);
		    ierr = DMLabelSetValue(faceLabel, cone[0], fID);CHKERRQ(ierr);
			  
		    ierr = DMLabelSetValue(bodyLabel, cone[1], b); CHKERRQ(ierr);
		    ierr = DMLabelSetValue(edgeLabel, cone[1], eid);CHKERRQ(ierr);
			  
		   ierr = DMLabelSetValue(bodyLabel, cone[2], b);CHKERRQ(ierr);
		   ierr = DMLabelSetValue(faceLabel, cone[2], fID);CHKERRQ(ierr);

		   cellCntr = cellCntr + 1;
		  }
		}
	  }
	  EG_free(lobjs);

	  ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, b);CHKERRQ(ierr);
	  ierr = DMLabelSetValue(faceLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, fID);CHKERRQ(ierr);
	}
	EG_free(fobjs);
  }
 
  ierr = PetscHMapIDestroy(&edgeMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyIndexMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyVertexMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeGlobalMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyFaceMap); CHKERRQ(ierr);
  
  *newdm = dm;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCreateEGADSTess(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
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
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) { 
    // ---------------------------------------------------------------------------------------------------
    // Generate Petsc Plex from EGADSlite created Tessellation of geometry
    // ---------------------------------------------------------------------------------------------------
	
	// Caluculate cell and vertex sizes 
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses); CHKERRQ(ierr);
	
    ierr = PetscHMapICreate(&pointIndexStartMap); CHKERRQ(ierr);	
	ierr = PetscHMapICreate(&triIndexStartMap); CHKERRQ(ierr);	
	ierr = PetscHMapICreate(&pTypeLabelMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&pIndexLabelMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&pBodyIndexLabelMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&triFaceIDLabelMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&triBodyIDLabelMap); CHKERRQ(ierr);

    /* Create Tessellation of Bodies */
    ego   tessArray[nbodies];
  
    for ( b = 0; b < nbodies; ++b) {
	  ego             body = bodies[b];
	  double          params[3] = {0.0, 0.0, 0.0};		// Parameters for Tessellation	  
	  int             Nf, bodyNumPoints = 0, bodyNumTris = 0;
	  PetscHashIter   PISiter, TISiter;
	  PetscBool       PISfound, TISfound;
	  
	  /* Store Start Index for each Body's Point and Tris */
	  ierr = PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound); CHKERRQ(ierr);
	  ierr = PetscHMapIFind(triIndexStartMap, b, &TISiter, &TISfound); CHKERRQ(ierr);
	  
	  if (!PISfound)  {ierr = PetscHMapISet(pointIndexStartMap, b, totalNumPoints); CHKERRQ(ierr);}
	  if (!TISfound)  {ierr = PetscHMapISet(triIndexStartMap, b, totalNumTris); CHKERRQ(ierr);}	  
	  
	  /* Calculate Tessellation parameters based on Bounding Box */
	  /* Get Bounding Box Dimensions of the BODY */
	  ierr = EG_getBoundingBox(body, boundBox);
	  tessSize = boundBox[3] - boundBox[0];
	  if ( tessSize < boundBox[4] - boundBox[1]) tessSize = boundBox[4] - boundBox[1];
	  if ( tessSize < boundBox[5] - boundBox[2]) tessSize = boundBox[5] - boundBox[2];
	  
	  // TODO :: May want to give users tessellation parameter options //
	  params[0] = 0.0250 * tessSize;
	  params[1] = 0.0075 * tessSize;
	  params[2] = 15.0;

	  ierr = EG_makeTessBody(body, params, &tessArray[b]); CHKERRQ(ierr);
	  
	  ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
	  	  
	  for (int f = 0; f < Nf; ++f) {
	    ego             face = fobjs[f];
		int             len, fID, ntris;
		const int      *ptype, *pindex, *ptris, *ptric;
		const double   *pxyz, *puv;
		  
		// Get Face ID //
		fID = EG_indexBodyTopo(body, face);
		  
		// Checkout the Surface Tessellation //
		ierr = EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric); CHKERRQ(ierr);
		  
		// Determine total number of triangle cells in the tessellation //
		bodyNumTris += (int) ntris;
		  
		// Check out the point index and coordinate //
		for (int p = 0; p < len; ++p) {
	      int global;
		  
		  ierr = EG_localToGlobal(tessArray[b], fID, p+1, &global);

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
  ierr = PetscMalloc2(totalNumPoints*cdim, &coords, totalNumTris*numCorners, &cells);CHKERRQ(ierr);
  
  for ( b = 0; b < nbodies; ++b) {
	ego             body = bodies[b];
	int             Nf;
	PetscInt        pointIndexStart; 
	PetscHashIter   PISiter;
	PetscBool       PISfound;
	
	ierr = PetscHMapIFind(pointIndexStartMap, b, &PISiter, &PISfound); CHKERRQ(ierr);
	if (!PISfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in pointIndexStartMap", b);}	  
	ierr = PetscHMapIGet(pointIndexStartMap, b, &pointIndexStart); CHKERRQ(ierr);
	
	ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
	  	  
	for (int f = 0; f < Nf; ++f) {
	  /* Get Face Object */
	  ego              face = fobjs[f];  
	  int              len, fID, ntris; 
	  const int       *ptype, *pindex, *ptris, *ptric;
	  const double    *pxyz, *puv;
		  
	  /* Get Face ID */
	  fID = EG_indexBodyTopo(body, face);
		  
	  /* Checkout the Surface Tessellation */
	  ierr = EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric); CHKERRQ(ierr);
		  
	  /* Check out the point index and coordinate */
	  for (int p = 0; p < len; ++p) {
		int              global;
		PetscHashIter    PTLiter, PILiter, PBLiter;
		PetscBool        PTLfound, PILfound, PBLfound;
		
		ierr = EG_localToGlobal(tessArray[b], fID, p+1, &global);
			
		/* Set the coordinates array for DAG */
		coords[((global-1+pointIndexStart)*3) + 0] = pxyz[(p*3)+0];
		coords[((global-1+pointIndexStart)*3) + 1] = pxyz[(p*3)+1];
		coords[((global-1+pointIndexStart)*3) + 2] = pxyz[(p*3)+2];
			
		/* Store Geometry Label Information for DMLabel assignment later */
		ierr = PetscHMapIFind(pTypeLabelMap, global-1+pointIndexStart, &PTLiter, &PTLfound); CHKERRQ(ierr);
        ierr = PetscHMapIFind(pIndexLabelMap, global-1+pointIndexStart, &PILiter, &PILfound); CHKERRQ(ierr);
        ierr = PetscHMapIFind(pBodyIndexLabelMap, global-1+pointIndexStart, &PBLiter, &PBLfound); CHKERRQ(ierr);

        if (!PTLfound)  {ierr = PetscHMapISet(pTypeLabelMap, global-1+pointIndexStart, ptype[p]); CHKERRQ(ierr);}	
        if (!PILfound)  {ierr = PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, pindex[p]); CHKERRQ(ierr);}	
        if (!PBLfound)  {ierr = PetscHMapISet(pBodyIndexLabelMap, global-1+pointIndexStart, b); CHKERRQ(ierr);}			  

		if (ptype[p] < 0) { ierr = PetscHMapISet(pIndexLabelMap, global-1+pointIndexStart, fID); CHKERRQ(ierr);}
	  }

	  for (int t = 0; t < (int) ntris; ++t){
		int             global, globalA, globalB;
		PetscHashIter   TFLiter, TBLiter;
	    PetscBool       TFLfound, TBLfound;
		  
		ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 0], &global);
		cells[(counter*3) +0] = global-1+pointIndexStart;

		ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 1], &globalA);
		cells[(counter*3) +1] = globalA-1+pointIndexStart;
		
		ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 2], &globalB);
		cells[(counter*3) +2] = globalB-1+pointIndexStart;
		  
		ierr = PetscHMapIFind(triFaceIDLabelMap, counter, &TFLiter, &TFLfound); CHKERRQ(ierr);
        ierr = PetscHMapIFind(triBodyIDLabelMap, counter, &TBLiter, &TBLfound); CHKERRQ(ierr);
		  
		if (!TFLfound)  {ierr = PetscHMapISet(triFaceIDLabelMap, counter, fID); CHKERRQ(ierr);}	
        if (!TBLfound)  {ierr = PetscHMapISet(triBodyIDLabelMap, counter, b); CHKERRQ(ierr);}

		counter += 1;
	  }
	}
	EG_free(fobjs);
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
    ierr = PetscContainerSetUserDestroy(contextObj, DMPlexEGADSDestroy_Private);CHKERRQ(ierr);
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
  for (int n = nStart; n < nEnd; ++n){
	int             pTypeVal, pIndexVal, pBodyVal;
    PetscHashIter   PTLiter, PILiter, PBLiter;
	PetscBool       PTLfound, PILfound, PBLfound;	
	  
	//Converted to Hash Tables
	ierr = PetscHMapIFind(pTypeLabelMap, n - nStart, &PTLiter, &PTLfound); CHKERRQ(ierr);
	if (!PTLfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pTypeLabelMap", n);}	  
	ierr = PetscHMapIGet(pTypeLabelMap, n - nStart, &pTypeVal); CHKERRQ(ierr);
	
	ierr = PetscHMapIFind(pIndexLabelMap, n - nStart, &PILiter, &PILfound); CHKERRQ(ierr);
	if (!PILfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pIndexLabelMap", n);}	  
	ierr = PetscHMapIGet(pIndexLabelMap, n - nStart, &pIndexVal); CHKERRQ(ierr);
	
	ierr = PetscHMapIFind(pBodyIndexLabelMap, n - nStart, &PBLiter, &PBLfound); CHKERRQ(ierr);
	if (!PBLfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in pBodyLabelMap", n);}	  
	ierr = PetscHMapIGet(pBodyIndexLabelMap, n - nStart, &pBodyVal); CHKERRQ(ierr);
	
	ierr = DMLabelSetValue(bodyLabel, n, pBodyVal); CHKERRQ(ierr);
	if(pTypeVal == 0) { ierr = DMLabelSetValue(vertexLabel, n, pIndexVal); CHKERRQ(ierr);}
	if(pTypeVal > 0) { ierr = DMLabelSetValue(edgeLabel, n, pIndexVal); CHKERRQ(ierr);}
	if(pTypeVal < 0) { ierr = DMLabelSetValue(faceLabel, n, pIndexVal); CHKERRQ(ierr);}
	  
  }
  
  /* Set DMLabels for Edges - Based on the DMLabels of the EDGE's NODES */
  for (int e = eStart; e < eEnd; ++e) {
	int    bodyID_0, vertexID_0, vertexID_1, edgeID_0, edgeID_1, faceID_0, faceID_1;
	  
	ierr = DMPlexGetCone(dm, e, &cone); CHKERRQ(ierr);
	ierr = DMLabelGetValue(bodyLabel, cone[0], &bodyID_0); CHKERRQ(ierr);		// Do I need to check the other end?
	ierr = DMLabelGetValue(vertexLabel, cone[0], &vertexID_0); CHKERRQ(ierr);
	ierr = DMLabelGetValue(vertexLabel, cone[1], &vertexID_1); CHKERRQ(ierr);
	ierr = DMLabelGetValue(edgeLabel, cone[0], &edgeID_0); CHKERRQ(ierr);
	ierr = DMLabelGetValue(edgeLabel, cone[1], &edgeID_1); CHKERRQ(ierr);
	ierr = DMLabelGetValue(faceLabel, cone[0], &faceID_0); CHKERRQ(ierr);
	ierr = DMLabelGetValue(faceLabel, cone[1], &faceID_1); CHKERRQ(ierr);
	  
	ierr = DMLabelSetValue(bodyLabel, e, bodyID_0); CHKERRQ(ierr);
	  
	if (edgeID_0 == edgeID_1) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_0); CHKERRQ(ierr); }
	else if (vertexID_0 > 0 && edgeID_1 > 0) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_1); CHKERRQ(ierr); }
	else if (vertexID_1 > 0 && edgeID_0 > 0) { ierr = DMLabelSetValue(edgeLabel, e, edgeID_0); CHKERRQ(ierr); }
	else { /* Do Nothing */ }
  }
  
  /* Set DMLabels for Cells */
  for (int f = fStart; f < fEnd; ++f){
	int             edgeID_0;
	PetscInt        triBodyVal, triFaceVal;
	PetscHashIter   TFLiter, TBLiter;
	PetscBool       TFLfound, TBLfound;
	
    // Convert to Hash Table
	ierr = PetscHMapIFind(triFaceIDLabelMap, f - fStart, &TFLiter, &TFLfound); CHKERRQ(ierr);
	if (!TFLfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triFaceIDLabelMap", f);}	  
	ierr = PetscHMapIGet(triFaceIDLabelMap, f - fStart, &triFaceVal); CHKERRQ(ierr);
	
	ierr = PetscHMapIFind(triBodyIDLabelMap, f - fStart, &TBLiter, &TBLfound); CHKERRQ(ierr);
	if (!TBLfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "DAG Point %d not found in triBodyIDLabelMap", f);}
    ierr = PetscHMapIGet(triBodyIDLabelMap, f - fStart, &triBodyVal); CHKERRQ(ierr);
	
	ierr = DMLabelSetValue(bodyLabel, f, triBodyVal); CHKERRQ(ierr);
	ierr = DMLabelSetValue(faceLabel, f, triFaceVal); CHKERRQ(ierr);
	  
	/* Finish Labeling previously unlabeled DMPlex Edges - Assumes Triangular Cell (3 Edges Max) */
	ierr = DMPlexGetCone(dm, f, &cone); CHKERRQ(ierr);
	  
	for (int jj = 0; jj < 3; ++jj) {
	  ierr = DMLabelGetValue(edgeLabel, cone[jj], &edgeID_0); CHKERRQ(ierr);
	  
	  if (edgeID_0 < 0) {
		ierr = DMLabelSetValue(bodyLabel, cone[jj], triBodyVal); CHKERRQ(ierr);
	  	ierr = DMLabelSetValue(faceLabel, cone[jj], triFaceVal); CHKERRQ(ierr);	
	  }
	}
  }
  
  *newdm = dm;
  PetscFunctionReturn(0);
}

/*@C
  DMPlexInflateToGeomModel - Snaps the vertix coordinates of a refined DMPlex object representing the mesh to its geometry if it was not performed during the refinement process

  Collective

  Input Parameters:
  dm       - The uninflated DM object representing the mesh

  Output Parameter:
. dm       - The inflated DM object representing the mesh

  Level: intermediate

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/

PETSC_EXTERN PetscErrorCode DMPlexInflateToGeomModel(DM dm)
{
	/* Variables for Function */
	//PETSC Variables
	PetscErrorCode  	ierr;
	DMLabel         	bodyLabel, faceLabel, edgeLabel, vertexLabel;
	PetscContainer  	modelObj;
	PetscInt        	cdim, bodyID, faceID, edgeID, vertexID;
	PetscInt        	vecSize, vStart, vEnd, p;
	PetscInt        	ix[3];
	double          	xyz[3];
	Vec             	vCoords, inflateCoords;
	// EGADS/EGADSlite Variables for Function
    ego             	model, geom, body, face, edge;
	ego				   *bodies;
	int             	Nb, oclass, mtype, *senses, pntCntr;
	PetscScalar         result[3];
	
	PetscFunctionBegin;
	/* Get Coordinate Dimension of DMPlex */
	ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);
	
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
		ierr = VecGetValues(vCoords, cdim, ix, xyz); CHKERRQ(ierr);    // Vertex coordinates (x,y,z)
				
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
		ierr = VecSetValues(inflateCoords, cdim, ix, result, INSERT_VALUES); CHKERRQ(ierr);
	}
	
	/* Prepare inflateCoords vector for DMPlex insertion */
	ierr = VecAssemblyBegin(inflateCoords); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(inflateCoords); CHKERRQ(ierr);
	
	/* Insert Inflated Coordinates into DMPlex */
	ierr = DMSetCoordinatesLocal(dm, inflateCoords); CHKERRQ(ierr);
	ierr = VecDestroy(&inflateCoords); CHKERRQ(ierr);
  
	PetscFunctionReturn(0);
}
#endif

/*@C
  DMPlexCreateEGADSFromFile - Create a DMPlex mesh from an EGADSLite file.

  Collective

  Input Parameters:
+ comm     - The MPI communicator
- filename - The name of the ExodusII file

  Output Parameter:
. dm       - The DM object representing the mesh

  Level: beginner

.seealso: DMPLEX, DMCreate(), DMPlexCreateEGADS()
@*/
PetscErrorCode DMPlexCreateEGADSFromFile(MPI_Comm comm, const char filename[], DM *dm)
{
  PetscMPIInt    rank;
#if defined(PETSC_HAVE_EGADS)
  ego            context= NULL, model = NULL;
#endif
  PetscBool      printModel = PETSC_FALSE, tessModel = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidCharPointer(filename, 2);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_print_model", &printModel, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_with_tess", &tessModel, NULL); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
#if defined(PETSC_HAVE_EGADS)
  if (!rank) {

    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, filename, &model);CHKERRQ(ierr);
    if (printModel) {ierr = DMPlexEGADSPrintModel(model);CHKERRQ(ierr);}

  }
  
  if (tessModel) {ierr = DMPlexCreateEGADSTess(comm, context, model, dm); CHKERRQ(ierr);}
  else {ierr = DMPlexCreateEGADS(comm, context, model, dm); CHKERRQ(ierr);}

#else
  SETERRQ(comm, PETSC_ERR_SUP, "This method requires EGADSLite support. Reconfigure using --download-egads");
#endif
  PetscFunctionReturn(0);
}
