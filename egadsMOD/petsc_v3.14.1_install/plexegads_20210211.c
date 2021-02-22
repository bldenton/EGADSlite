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
  ego            geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *shobjs, *fobjs, *eobjs;
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

    ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &objs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);
	
	EG_free(shobjs);
	EG_free(fobjs);
	EG_free(lobjs);
	EG_free(eobjs);
	EG_free(objs);
	
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
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexEGADSDestroy_Private(void *context)
{
  if (context) EG_close((ego) context);
  return 0;
}

static PetscErrorCode DMPlexCreateEGADS(MPI_Comm comm, ego context, ego model, DM *newdm)
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel;
  /* EGADSLite variables */
  ego            geom, *bodies, *mobjs, *objs, *fobjs, *eobjs, *nobjs;
  int            oclass, mtype, nbodies, *senses;
  int            b;
  /* PETSc variables */
  DM             dm;
  PetscHMapI     edgeMap = NULL, bodyIndexMap = NULL, bodyVertexMap = NULL, bodyEdgeMap = NULL, bodyFaceMap = NULL, bodyEdgeGlobalMap = NULL;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numEdges = 0, numFaces = 0, numCells = 0, edgeCntr = 0;
  PetscInt       cellCntr = 0, numPoints = 0; 
  PetscInt      *cells  = NULL;
  PetscReal     *coords = NULL;
  PetscMPIInt    rank;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) { 
    /* ---------------------------------------------------------------------------------------------------
    Generate Petsc Plex
      Get all Nodes in model, record coordinates in a correctly formatted array
      Cycle through bodies, cycle through loops, recorde NODE IDs in a correctly formatted array
      We need to uniformly refine the initial geometry to guarantee a valid mesh
    */
	
	/* Caluculate cell and vertex sizes */
	ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses); CHKERRQ(ierr);
	
    ierr = PetscHMapICreate(&edgeMap); CHKERRQ(ierr);	
	ierr = PetscHMapICreate(&bodyIndexMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyVertexMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyEdgeMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyEdgeGlobalMap); CHKERRQ(ierr);
	ierr = PetscHMapICreate(&bodyFaceMap); CHKERRQ(ierr);
	
	for (b = 0; b < nbodies; ++b) {
      ego           body = bodies[b];
	  int           Nf, Ne, Nv;
	  PetscHashIter BIiter, BViter, BEiter, BEGiter, BFiter, EMiter;
	  PetscBool     BIfound, BVfound, BEfound, BEGfound, BFfound, EMfound;
	  
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
	  
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[ Nf = %d :: Ne = %d :: Nv = %d]\n", Nf, Ne, Nv); CHKERRQ(ierr);
	  
	  /* Remove DEGENERATE EDGES from Edge count */
	  ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  int Netemp = 0;
	  for (int e = 0; e < Ne; ++e) {
		  ego     edge = eobjs[e];
		  ego     topRef, prev, next;
		  int     eid;
		  
		  ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		  eid = EG_indexBodyTopo(body, edge); //CHKERRQ(ierr);
		  
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
	  
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[ line 336 ]\n"); CHKERRQ(ierr);
	  
	  /* Determine Number of Cells */
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_SELF, "  Nf = %d \n", Nf); CHKERRQ(ierr);
	  for (int f = 0; f < Nf; ++f) {
        ego face = fobjs[f];
	  	
	  	ierr = PetscPrintf(PETSC_COMM_SELF, "  f = %d \n", f); CHKERRQ(ierr);		
	  	ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);	
        //ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &senses);CHKERRQ(ierr);		
        ierr = PetscPrintf(PETSC_COMM_SELF, "  Ne = %d \n", Ne); CHKERRQ(ierr);
	  
	  	int edgeTemp = 0;
	  	for (int e = 0; e < Ne; ++e) {
	  	  ego   edge = eobjs[e];
	  	  ego   topRef, prev, next;
	  	  
	  	  ierr = PetscPrintf(PETSC_COMM_SELF, "  e = %d \n", e); CHKERRQ(ierr);
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
	ierr = PetscPrintf(PETSC_COMM_SELF, "  NumFaces = %d \n", numFaces);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "  numEdges = %d \n", numEdges);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "  numVertices = %d \n", numVertices);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "  edgeCntr = %d \n", edgeCntr);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF, "  numCells = %d \n", numCells);CHKERRQ(ierr);
	
	/* Set up basic DMPlex parameters */
	dim        = 2;		/* Assumes 3D Models :: Need to handle 2D modles in the future */
	cdim       = 3;     /* Assumes 3D Models :: Need to update to handle 2D modles in future */
	numCorners = 3;     /* Split Faces into triangles */
    numPoints  = numVertices + numEdges + numFaces;   /* total number of coordinate points */

	ierr = PetscMalloc2(numPoints*cdim, &coords, numCells*numCorners, &cells); CHKERRQ(ierr);
	
	/* Get Vertex Coordinates and Set up Cells */
	for (b = 0; b < nbodies; ++b) {
	  ego           body = bodies[b];
	  int           Nf, Ne, Nv; 
	  PetscInt      bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
	  PetscHashIter BViter, BEiter, BEGiter, BFiter, EMiter;
	  PetscBool     BVfound, BEfound, BEGfound, BFfound, EMfound;
	  
	  /* Vertices on Current Body */
	  ierr = EG_getBodyTopos(body, NULL, NODE, &Nv, &nobjs); //CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyVertexMap, b, &BViter, &BVfound); CHKERRQ(ierr);
	  if (!BVfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyVertexMap", b);}	  
	  ierr = PetscHMapIGet(bodyVertexMap, b, &bodyVertexIndexStart); CHKERRQ(ierr);
	  
	  for (int v = 0; v < Nv; ++v) {
	    ego    vertex = nobjs[v];
		double limits[4];
		int    id, dummy;
		
		ierr = EG_getTopology(vertex, &geom, &oclass, &mtype, limits, &dummy, &mobjs, &senses); CHKERRQ(ierr);
		id = EG_indexBodyTopo(body, vertex); //CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  coords(%d) filled \n", bodyVertexIndexStart + id - 1); CHKERRQ(ierr);
		coords[(bodyVertexIndexStart + id - 1)*cdim + 0] = limits[0];
		coords[(bodyVertexIndexStart + id - 1)*cdim + 1] = limits[1];
		coords[(bodyVertexIndexStart + id - 1)*cdim + 2] = limits[2];
	  }
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(nobjs) line 413\n"); CHKERRQ(ierr);
	  EG_free(nobjs);
	  //EG_free(mobjs);
	  //EG_free(senses);
	  
	  /* Edge Midpoint Vertices on Current Body */
	  ierr = EG_getBodyTopos(body, NULL, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyEdgeMap, b, &BEiter, &BEfound); CHKERRQ(ierr);
	  if (!BEfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeMap", b);}	  
	  ierr = PetscHMapIGet(bodyEdgeMap, b, &bodyEdgeIndexStart); CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyEdgeGlobalMap, b, &BEGiter, &BEGfound); CHKERRQ(ierr);
	  if (!BEGfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyEdgeGlobalMap", b);}	  
	  ierr = PetscHMapIGet(bodyEdgeGlobalMap, b, &bodyEdgeGlobalIndexStart); CHKERRQ(ierr);
	  
	  for (int e = 0; e < Ne; ++e) {
	    ego          edge = eobjs[e];
		ego          topRef, prev, next;
		double       range[2], avgt[1], cntrPnt[9];
		int          eid, eOffset;
		int          periodic;
		
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {continue;}
		
		eid = EG_indexBodyTopo(body, edge); // CHKERRQ(ierr);
		
		/* get relative offset from globalEdgeID Vector */
		ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);
	    if (!EMfound) {SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d not found in edgeMap", bodyEdgeGlobalIndexStart + eid - 1);}	  
	    ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset); CHKERRQ(ierr);
		
		ierr = EG_getRange(edge, range, &periodic); CHKERRQ(ierr);
		avgt[0] = (range[0] + range[1]) /  2.;
		
		ierr = EG_evaluate(edge, avgt, cntrPnt); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  coords(%d) filled \n", numVertices + bodyEdgeIndexStart + eOffset - 1);CHKERRQ(ierr);
		coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 0] = cntrPnt[0];
        coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 1] = cntrPnt[1];
		coords[(numVertices + bodyEdgeIndexStart + eOffset - 1)*cdim + 2] = cntrPnt[2];
	  }
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(eobjs) line 451\n");CHKERRQ(ierr);
	  EG_free(eobjs);

	  /* Face Midpoint Vertices on Current Body */
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); //CHKERRQ(ierr);
	  
	  ierr = PetscHMapIFind(bodyFaceMap, b, &BFiter, &BFfound); CHKERRQ(ierr);
	  if (!BFfound) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "Body %d not found in bodyFaceMap", b); 
	  ierr = PetscHMapIGet(bodyFaceMap, b, &bodyFaceIndexStart); CHKERRQ(ierr);
	  
	  for (int f = 0; f < Nf; ++f) {
		ego     face = fobjs[f];
		double  range[4], avgUV[2], cntrPnt[9];
		int     peri, id;
		
		id = EG_indexBodyTopo(body, face); 
		ierr = EG_getRange(face, range, &peri); CHKERRQ(ierr);
		
		avgUV[0] = (range[0] + range[1]) / 2.;
		avgUV[1] = (range[2] + range[3]) / 2.;
		ierr = EG_evaluate(face, avgUV, cntrPnt); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  coords(%d) filled \n", numVertices + numEdges + bodyFaceIndexStart + id - 1); CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "     [x, y, z) = [%lf, %lf, %lf] \n", cntrPnt[0], cntrPnt[1], cntrPnt[2]); CHKERRQ(ierr);
		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 0] = cntrPnt[0];
		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 1] = cntrPnt[1];
		coords[(numVertices + numEdges + bodyFaceIndexStart + id - 1)*cdim + 2] = cntrPnt[2];
	  }
	  EG_free(fobjs);
	  
	  /* Define Cells :: Note - This could be incorporated in the Face Midpoint Vertices Loop but was kept separate for clarity */
	  ierr = EG_getBodyTopos(body, NULL, FACE, &Nf, &fobjs); 
	  for (int f = 0; f < Nf; ++f) {
		ego           face = fobjs[f];
		//ego          *lobjs;
		int           fID, midFaceID, midPntID, startID, endID, Nl;
		int          *lSenses;
		
		fID = EG_indexBodyTopo(body, face); //CHKERRQ(ierr);
		midFaceID = numVertices + numEdges + bodyFaceIndexStart + fID - 1;
		ierr = PetscPrintf(PETSC_COMM_SELF, "  fID = %d \n", fID);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  Nf = %d \n", Nf);CHKERRQ(ierr);
		/* Must Traverse Loop to ensure we have all necessary information like the sense (+/- 1) of the edges. */
		/* TODO :: Only handles single loop faces (No holes). The choices for handling multiloop faces are:    */
		/*            1) Use the DMPlexCreateEGADSFromFile() with the -dm_plex_egads_with_tess = 1 option.     */
		/*               This will use a default EGADS tessellation as an initial surface mesh.                */
		/*            2) Create the initial surface mesh via a 2D mesher :: Currently not availble (?future?)  */
		/*               May I suggest the XXXX as a starting point?                                           */
		
		ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &objs, &lSenses); //CHKERRQ(ierr);
		//ierr = EG_getBodyTopos(body, face, LOOP, &Nl, &objs); //CHKERRQ(ierr);
		
	    if (Nl > 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face has %d Loops. Can only handle Faces with 1 Loop. Please use --dm_plex_egads_with_tess = 1 Option", Nl);
		ierr = PetscPrintf(PETSC_COMM_SELF, "  Nl = %d \n", Nl);CHKERRQ(ierr);
		for (int l = 0; l < Nl; ++l) {
          ego      loop = objs[l];
		  int      lid, *eSenses;
		  lid = EG_indexBodyTopo(body, loop);
		  
          ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &eobjs, &eSenses); //CHKERRQ(ierr);
		  //ierr = EG_getBodyTopos(body, loop, EDGE, &Ne, &eobjs);
		  
		  ierr = PetscPrintf(PETSC_COMM_SELF, "  l = %d \n", l); //CHKERRQ(ierr);
		  ierr = PetscPrintf(PETSC_COMM_SELF, "  lid = %d \n", lid); //CHKERRQ(ierr);
		  ierr = PetscPrintf(PETSC_COMM_SELF, "  Ne = %d \n", Ne); //CHKERRQ(ierr);
		  for (int e = 0; e < Ne; ++e) {
		    ego   edge = eobjs[e];
		    ego   topRef, prev, next;
			ego   *nobjs;
		    int   eid, eOffset;
		    
			ierr = PetscPrintf(PETSC_COMM_SELF, "    e = %d \n", e); //CHKERRQ(ierr);
			
		    ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); //CHKERRQ(ierr);
			eid = EG_indexBodyTopo(body, edge);	// WHY IS THIS GIVING A NEGATIVE NUMBER ON THE 3rd LOOP??????
		    if (mtype == DEGENERATE) {
				ierr = PetscPrintf(PETSC_COMM_SELF, "  Edge %d is DEGENERATE \n", eid); //CHKERRQ(ierr);
				continue;
			}
		    //eid = EG_indexBodyTopo(body, edge);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    eid = %d \n", eid); //CHKERRQ(ierr);
			
		    /* get relative offset from globalEdgeID Vector */
		    ierr = PetscHMapIFind(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &EMiter, &EMfound); CHKERRQ(ierr);
	        if (!EMfound) {SETERRQ3(PETSC_COMM_SELF, PETSC_ERR_SUP, "Edge %d of Body %d not found in edgeMap. Global Edge ID :: %d", eid, b, bodyEdgeGlobalIndexStart + eid - 1);}	  
	        ierr = PetscHMapIGet(edgeMap, bodyEdgeGlobalIndexStart + eid - 1, &eOffset); CHKERRQ(ierr);
			
			ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyEdgeGlobalIndexStart = %d \n", bodyEdgeGlobalIndexStart); CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    eOffset = %d \n", eOffset); //CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    eSenses = %d \n", eSenses[e]); //CHKERRQ(ierr);
			
			midPntID = numVertices + bodyEdgeIndexStart + eOffset - 1;
		  
		    ierr = EG_getTopology(edge, &geom, &oclass, &mtype, NULL, &Nv, &nobjs, &senses); CHKERRQ(ierr);
			
		    if (eSenses[e] > 0) { startID = EG_indexBodyTopo(body, nobjs[0]); endID = EG_indexBodyTopo(body, nobjs[1]); }
		    else { startID = EG_indexBodyTopo(body, nobjs[1]); endID = EG_indexBodyTopo(body, nobjs[0]); }
			
			/* Define 2 Cells per Edge with correct orientation */
			ierr = PetscPrintf(PETSC_COMM_SELF, "  cells(%d) filled \n", cellCntr*numCorners);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    Edge ID = %d \n", eid);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    startID = %d \n", startID);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    endID = %d \n", endID);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    midFaceID = %d \n", midFaceID);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyVertexIndexStart = %d \n", bodyVertexIndexStart);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyVertexIndexStart + startID - 1 = %d \n", bodyVertexIndexStart + startID - 1);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    midPntID = %d \n", midPntID);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF, "    bodyVertexIndexStart + endID - 1 = %d \n", bodyVertexIndexStart + endID - 1);CHKERRQ(ierr);
			
			cells[cellCntr*numCorners + 0] = midFaceID;
			cells[cellCntr*numCorners + 1] = bodyVertexIndexStart + startID - 1;
			cells[cellCntr*numCorners + 2] = midPntID;
			
			cells[cellCntr*numCorners + 3] = midFaceID;
			cells[cellCntr*numCorners + 4] = midPntID;
			cells[cellCntr*numCorners + 5] = bodyVertexIndexStart + endID - 1;
			
			cellCntr = cellCntr + 2;
			
			ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(nobjs) line 532\n");CHKERRQ(ierr);
			//EG_free(nobjs);
			//EG_free(lobjs);
		  }
		  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(eobjs) line 535\n");CHKERRQ(ierr);
		  EG_free(eobjs);
		}
		EG_free(objs);
	  }
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(fobjs) line 538\n");CHKERRQ(ierr);
	  EG_free(fobjs);
	}
  }
  CHKMEMQ;
  /* Generate DMPlex */
  ierr = PetscPrintf(PETSC_COMM_SELF, "dim = %d\n", dim);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "numCells = %d\n", numCells);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "numPoints = %d\n", numPoints);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "numCorners = %d\n", numCorners);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "cdim = %d\n", cdim);CHKERRQ(ierr);
  
  ierr = DMPlexCreateFromCellListPetsc(PETSC_COMM_WORLD, dim, numCells, numPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dm); CHKERRQ(ierr);
  ierr = PetscFree2(coords, cells); CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Cells    = %D \n", numCells); CHKERRQ(ierr);
  ierr = PetscInfo1(dm, " Total Number of Unique Vertices = %D \n", numVertices); CHKERRQ(ierr);
  
  /* Embed EGADS model in DM */
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
  /* Label points */
  PetscInt  fStart, fEnd, eStart, eEnd, nStart, nEnd;
 
  ierr = DMCreateLabel(dm, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dm, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr);
  
  ierr = DMPlexGetHeightStratum(dm, 0, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 2, &nStart, &nEnd);CHKERRQ(ierr);
  
  cellCntr = 0;
  for (b = 0; b < nbodies; ++b) {
    ego           body = bodies[b];
	ego          *fobjs;
	PetscInt      Nv, Ne, Nf, bodyVertexIndexStart, bodyEdgeIndexStart, bodyEdgeGlobalIndexStart, bodyFaceIndexStart;
	PetscHashIter BViter, BEiter, BEGiter, BFiter, EMiter;
	PetscBool     BVfound, BEfound, BEGfound, BFfound, EMfound;
	  
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
	  ego  *eobjs;
      int   fID;
	  
	  fID  = EG_indexBodyTopo(body, face); CHKERRQ(ierr);
      ierr = EG_getBodyTopos(body, face, EDGE, &Ne, &eobjs); CHKERRQ(ierr);
	  for (int e = 0; e < Ne; ++e) {
        ego           edge = eobjs[e];
        ego          *nobjs;
		ego           topRef, prev, next;
        int           eid, eOffset;
		
		/* Skip DEGENERATE Edges */
		ierr = EG_getInfo(edge, &oclass, &mtype, &topRef, &prev, &next); CHKERRQ(ierr);
		if (mtype == DEGENERATE) {continue;}
		eid = EG_indexBodyTopo(body, edge); CHKERRQ(ierr);
		
		/* get relative offset from globalEdgeID Vector */
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
		ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(nobjs) line 634\n");CHKERRQ(ierr);
		//EG_free(nobjs);

		ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, b);CHKERRQ(ierr);
		ierr = DMLabelSetValue(edgeLabel, nStart + numVertices + bodyEdgeIndexStart + eOffset - 1, eid);CHKERRQ(ierr);
		
		// Can I define Cell faces HERE?  -- Let's try
		for (int jj = 0; jj < 2; ++jj){
          const PetscInt   *cone = NULL;
		  
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
	  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(eobjs) line 660\n");CHKERRQ(ierr);
	  //EG_free(eobjs);

	  ierr = DMLabelSetValue(bodyLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, b);CHKERRQ(ierr);
	  ierr = DMLabelSetValue(faceLabel, nStart + numVertices + numEdges + bodyFaceIndexStart + fID - 1, fID);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(fobjs) line 666\n");CHKERRQ(ierr);
	//EG_free(fobjs);
  }
  ierr = PetscPrintf(PETSC_COMM_SELF, "[AFTER EG_free(bodies) line 669\n");CHKERRQ(ierr);
  //EG_free(bodies);
  
  ierr = PetscHMapIDestroy(&edgeMap);CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyIndexMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyVertexMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyEdgeGlobalMap); CHKERRQ(ierr);
  ierr = PetscHMapIDestroy(&bodyFaceMap); CHKERRQ(ierr);
  
  *newdm = dm;
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
  PetscBool      printModel = PETSC_FALSE; //, tessModel = PETSC_FALSE;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidCharPointer(filename, 2);
  ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_print_model", &printModel, NULL);CHKERRQ(ierr);
  //ierr = PetscOptionsGetBool(NULL, NULL, "-dm_plex_egads_with_tess", &tessModel, NULL); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
#if defined(PETSC_HAVE_EGADS)
  if (!rank) {

    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, filename, &model);CHKERRQ(ierr);
    if (printModel) {ierr = DMPlexEGADSPrintModel(model);CHKERRQ(ierr);}

  }
  ierr = DMPlexCreateEGADS(comm, context, model, dm);CHKERRQ(ierr);
#else
  SETERRQ(comm, PETSC_ERR_SUP, "This method requires EGADSLite support. Reconfigure using --download-egads");
#endif
  PetscFunctionReturn(0);
}
