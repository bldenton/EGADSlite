static const char help[] = "Test of EGADSLite CAD functionality";

#include <petscdmplex.h>

#include <egads.h>
#include <petsc.h>

void surfArea(DM dm);

typedef struct {
  char filename[PETSC_MAX_PATH_LEN];
} AppCtx;

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  options->filename[0] = '\0';

  ierr = PetscOptionsBegin(comm, "", "EGADSPlex Problem Options", "EGADSLite");CHKERRQ(ierr);
  ierr = PetscOptionsString("-filename", "The EGADSLite file", "ex9.c", options->filename, options->filename, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[])
{
  DMLabel        bodyLabel, faceLabel, edgeLabel, vertexLabel; // markerLabel;
  PetscInt       cStart, cEnd;
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs, *fobjs, *eobjs, *shobjs;
  ego            *tess;
  int            oclass, mtype, nbodies, *senses, *bsenses, *shsenses, *fsenses, *lsenses, *esenses;
  int            b, maxNumEdges;
  double         boundBox[6], tessSize;
  /* PETSc variables */
  DM             dm, dmNozzle, dmMesh;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0, depth = 0;
  PetscInt       numFaces = 0, numEdges = 0, numPoints = 0;
  PetscInt      *cells  = NULL, *coneOrient = NULL, *cone = NULL, *coneSize = NULL;
  PetscReal     *coords = NULL;
  MPI_Comm       comm;
  PetscMPIInt    rank;
  AppCtx         ctx;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = ProcessOptions(comm, &ctx);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    /* Open EGADs file and load EGADs model data */
    ierr = EG_open(&context);CHKERRQ(ierr);
    ierr = EG_loadModel(context, 0, ctx.filename, &model);CHKERRQ(ierr);

    /* test bodyTopo functions */
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &bsenses);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", nbodies);CHKERRQ(ierr);
	
	/* Initialize Variable to determine the maximum number of EDGES in any BODY */
	maxNumEdges = 0;

    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nsh, sh, Nf, f, Nl, l, Ne, e, Nv, v;
	  
	  /* Get Bounding Box Dimensions of the BODY */
	  ierr = EG_getBoundingBox(body, boundBox);

      /* Output Basic Model Topology */
      ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &shobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of EDGES: %d \n", Ne);CHKERRQ(ierr);
	  
	  /* Check to see if current BODY has more EDGEs than the previous BODY */
	  maxNumEdges = PetscMax(maxNumEdges, Ne);

      ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);
      
      /* Get SHELL info which associated with the current BODY */
      ierr = EG_getTopology(body, &geom, &oclass, &mtype, NULL, &Nsh, &shobjs, &shsenses);CHKERRQ(ierr);

      for (sh = 0; sh < Nsh; ++sh){
        ego shell = shobjs[sh];
        int shsense = shsenses[sh];
        
        id   = EG_indexBodyTopo(body, shell);
        ierr = PetscPrintf(PETSC_COMM_SELF, "      SHELL ID: %d :: sense = %d\n", id, shsense);CHKERRQ(ierr);
        
        /* Get FACE info which associated with the current SHELL */
        ierr = EG_getTopology(shell, &geom, &oclass, &mtype, NULL, &Nf, &fobjs, &fsenses);CHKERRQ(ierr);
        
        for (f = 0; f < Nf; ++f){
          ego face = fobjs[f];
          int fsense = fsenses[f];
          
          id   = EG_indexBodyTopo(body, face);
          ierr = PetscPrintf(PETSC_COMM_SELF, "        FACE ID: %d \n", id);CHKERRQ(ierr);
          
          /* Get LOOP info which associated with the current FACE */
          ierr = EG_getTopology(face, &geom, &oclass, &mtype, NULL, &Nl, &objs, &lsenses);CHKERRQ(ierr);
          
          for (l = 0; l < Nl; ++l) {
            ego loop = objs[l];
            int lsense = lsenses[l];
  
            id   = EG_indexBodyTopo(body, loop);
            ierr = PetscPrintf(PETSC_COMM_SELF, "          LOOP ID: %d :: sense = %d\n", id, lsense);CHKERRQ(ierr);
    
            /* Get EDGE info which associated with the current LOOP */
            ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &esenses);CHKERRQ(ierr);
			
            for (e = 0; e < Ne; ++e) {
              ego edge = objs[e];
  			  int esense = esenses[e];
    
              id = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
              ierr = PetscPrintf(PETSC_COMM_SELF, "            EDGE ID: %d :: sense = %d\n", id, esense);CHKERRQ(ierr);
    
              double range[4] = {0., 0., 0., 0.};
              int    peri;
    
              ierr = EG_getRange(objs[e], range, &peri);
              ierr = PetscPrintf(PETSC_COMM_SELF, "              Range = %lf, %lf, %lf, %lf \n", range[0], range[1], range[2], range[3]);
			  
              /* Get NODE info which associated with the current EDGE */
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
  }
  
  /* Create Tessellation of Bodies */
  ego *tessArray[nbodies];
  int totalNumTris = 0, totalNumPoints = 0;
  int pointIndexStart[nbodies], triIndexStart[nbodies];
  
  for ( b = 0; b < nbodies; ++b) {
	  ego body = bodies[b];
	  double params[3] = {0.0, 0.0, 0.0};		// Parameters for Tessellation
	  
	  int f, Nf;
	  int bodyNumPoints = 0, bodyNumTris = 0;
	  
	  /* Store Start Index for each Body's Point and Tris */
	  pointIndexStart[b] = totalNumPoints;
	  triIndexStart[b] = totalNumTris;
	  
	  /* Calculate Tessellation parameters based on Bounding Box */
	  tessSize = boundBox[3] - boundBox[0];
	  if ( tessSize < boundBox[4] - boundBox[1]) tessSize = boundBox[4] - boundBox[1];
	  if ( tessSize < boundBox[5] - boundBox[2]) tessSize = boundBox[5] - boundBox[2];
	  
	  params[0] = 0.025 * tessSize;
	  params[1] = 0.001 * tessSize;
	  params[2] = 15.0;
	  
	  //ierr = EG_makeTessBody(body, params, &tess);
	  ierr = EG_makeTessBody(body, params, &tessArray[b]); CHKERRQ(ierr);
	  
	  ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
	  	  
	  for ( f = 0; f < Nf; ++f) {
		  /* Get Face Object */
		  ego face = fobjs[f];
		  
		  int len, *ptype, *pindex, *ntris, *ptris, *ptric;
		  double *pxyz, *puv;
		  
		  /* Get Face ID */
		  int fID = EG_indexBodyTopo(body, face);
		  
		  /* Checkout the Surface Tessellation */
		  //ierr = EG_getTessFace(tess, fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric);
		  ierr = EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric); CHKERRQ(ierr);
		  
		  /* Determine total number of triangle cells in the tessellation */
		  bodyNumTris = bodyNumTris + (int) ntris;
		  
		  /* Check out the point index and coordinate */
		  for (int p = 0; p < len; ++p) {
			int global;
			//ierr = EG_localToGlobal(tess, fID, p+1, &global);
			ierr = EG_localToGlobal(tessArray[b], fID, p+1, &global);

			/* Determine the total number of points in the tessellation */
			bodyNumPoints = PetscMax(bodyNumPoints, global);
		  }
	  }
	  totalNumPoints = totalNumPoints + bodyNumPoints;
	  totalNumTris = totalNumTris + bodyNumTris;		  
  }
  
  dim = 2;
  cdim = 3;
  numCorners = 3;
  PetscInt counter = 0;
  
  /* NEED TO DEFINE MATRICES/VECTORS TO STORE GEOM REFERENCE DATA   */
  /* Fill in below and use to define DMLabels after DMPlex creation */
  int pointLabels[totalNumPoints][3];			// ptype, pindex, bodyIndex
  int triLabels[totalNumTris][2];				// fID, bodyIndex
  
  ierr = PetscMalloc2(totalNumPoints*cdim, &coords, totalNumTris*numCorners, &cells);CHKERRQ(ierr);
  
  for ( b = 0; b < nbodies; ++b) {
	  ego body = bodies[b];
	  
	  int f, Nf;
	  
	  ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
	  	  
	  for ( f = 0; f < Nf; ++f) {
		  /* Get Face Object */
		  ego face = fobjs[f];
		  
		  int len, *ptype, *pindex, *ntris, *ptris, *ptric;
		  double *pxyz, *puv;
		  
		  /* Get Face ID */
		  int fID = EG_indexBodyTopo(body, face);
		  
		  /* Checkout the Surface Tessellation */
		  //ierr = EG_getTessFace(tess, fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric); CHKERRQ(ierr);
		  ierr = EG_getTessFace(tessArray[b], fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric); CHKERRQ(ierr);
		  
		  /* Check out the point index and coordinate */
		  for (int p = 0; p < len; ++p) {
			int global;
			//ierr = EG_localToGlobal(tess, fID, p+1, &global);
			ierr = EG_localToGlobal(tessArray[b], fID, p+1, &global);
			
			/* Set the coordinates array for DAG */
			coords[((global-1+pointIndexStart[b])*3) + 0] = pxyz[(p*3)+0];
			coords[((global-1+pointIndexStart[b])*3) + 1] = pxyz[(p*3)+1];
			coords[((global-1+pointIndexStart[b])*3) + 2] = pxyz[(p*3)+2];
			
			/* Store Geometry Label Information for DMLabel assignment later */
			pointLabels[global-1+pointIndexStart[b]][0] = ptype[p];
			pointLabels[global-1+pointIndexStart[b]][1] = pindex[p];
			pointLabels[global-1+pointIndexStart[b]][2] = b;
			
			if (ptype[p] < 0) { pointLabels[global-1+pointIndexStart[b]][1] = fID; }
		  }
		  
		  for (int t = 0; t < (int) ntris; ++t){
			  int global, globalA, globalB;
			  //ierr = EG_localToGlobal(tess, fID, ptris[(t*3) + 0], &global);
			  ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 0], &global);
			  cells[(counter*3) +0] = global-1+pointIndexStart[b];
			  //ierr = EG_localToGlobal(tess, fID, ptris[(t*3) + 1], &globalA);
			  ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 1], &globalA);
			  cells[(counter*3) +1] = globalA-1+pointIndexStart[b];
			  //ierr = EG_localToGlobal(tess, fID, ptris[(t*3) + 2], &globalB);
			  ierr = EG_localToGlobal(tessArray[b], fID, ptris[(t*3) + 2], &globalB);
			  cells[(counter*3) +2] = globalB-1+pointIndexStart[b];
			  			  
			  triLabels[counter][0] = fID;
			  triLabels[counter][1] = b;
			  
			  counter = counter + 1;
		  }
	  }
  } 

  //Build DMPlex   
  ierr = DMPlexCreateFromCellList(PETSC_COMM_WORLD, dim, totalNumTris, totalNumPoints, numCorners, PETSC_TRUE, cells, cdim, coords, &dmNozzle);CHKERRQ(ierr); 
  ierr = PetscFree2(coords, cells);CHKERRQ(ierr);

  /* Attached the EGADS Model to the Surface Plex */
  {
    PetscContainer modelObj;

    ierr = PetscContainerCreate(PETSC_COMM_SELF, &modelObj);CHKERRQ(ierr);
    ierr = PetscContainerSetPointer(modelObj, model);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject) dmNozzle, "EGADS Model", (PetscObject) modelObj);CHKERRQ(ierr);
    ierr = PetscContainerDestroy(&modelObj);CHKERRQ(ierr);
  }
  
  ierr = DMCreateLabel(dmNozzle, "EGADS Body ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Face ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Edge ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMCreateLabel(dmNozzle, "EGADS Vertex ID");CHKERRQ(ierr);
  ierr = DMGetLabel(dmNozzle, "EGADS Vertex ID", &vertexLabel);CHKERRQ(ierr); 
  
   /* Get Number of DAG Nodes at each level */
  int fDAGlevel, eDAGlevel, nDAGlevel;
  int fStart, fEnd, eStart, eEnd, nStart, nEnd;
  
  ierr = DMPlexGetHeightStratum(dmNozzle, 0, &fStart, &fEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dmNozzle, 1, &eStart, &eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dmNozzle, 2, &nStart, &nEnd);CHKERRQ(ierr);
  
  /* Set DMLabels for NODES */
  for (int ii = nStart; ii < nEnd; ++ii){
	  ierr = DMLabelSetValue(bodyLabel, ii, pointLabels[ii-nStart][2]); CHKERRQ(ierr);
	  if(pointLabels[ii-nStart][0] == 0) { ierr = DMLabelSetValue(vertexLabel, ii, pointLabels[ii-nStart][1]); CHKERRQ(ierr);}
	  if(pointLabels[ii-nStart][0] > 0) { ierr = DMLabelSetValue(edgeLabel, ii, pointLabels[ii-nStart][1]); CHKERRQ(ierr);}
	  if(pointLabels[ii-nStart][0] < 0) { ierr = DMLabelSetValue(faceLabel, ii, pointLabels[ii-nStart][1]); CHKERRQ(ierr);}
  }
  
  /* Set DMLabels for Edges - Based on the DMLabels of the EDGE's NODES */
  for (int ii = eStart; ii < eEnd; ++ii) {
	  int bodyID_0, vertexID_0, vertexID_1, edgeID_0, edgeID_1, faceID_0, faceID_1;
	  
	  ierr = DMPlexGetCone(dmNozzle, ii, &cone); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(bodyLabel, cone[0], &bodyID_0); CHKERRQ(ierr);		// Do I need to check the other end?
	  ierr = DMLabelGetValue(vertexLabel, cone[0], &vertexID_0); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(vertexLabel, cone[1], &vertexID_1); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(edgeLabel, cone[0], &edgeID_0); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(edgeLabel, cone[1], &edgeID_1); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(faceLabel, cone[0], &faceID_0); CHKERRQ(ierr);
	  ierr = DMLabelGetValue(faceLabel, cone[1], &faceID_1); CHKERRQ(ierr);
	  
	  ierr = DMLabelSetValue(bodyLabel, ii, bodyID_0); CHKERRQ(ierr);
	  
	  if (edgeID_0 == edgeID_1){
		  ierr = DMLabelSetValue(edgeLabel, ii, edgeID_0); CHKERRQ(ierr);
	  } else if (vertexID_0 > 0 && edgeID_1 > 0){
		  ierr = DMLabelSetValue(edgeLabel, ii, edgeID_1); CHKERRQ(ierr);
	  } else if (vertexID_1 > 0 && edgeID_0 > 0){
		  ierr = DMLabelSetValue(edgeLabel, ii, edgeID_0); CHKERRQ(ierr);
	  } else {
		  // Do Nothing
	  }
  }
  
  /* Set DMLabels for Cells */
  for (int ii = fStart; ii < fEnd; ++ii){
	  int edgeID_0;
	  
	  ierr = DMLabelSetValue(bodyLabel, ii, triLabels[ii-fStart][1]); CHKERRQ(ierr);
	  ierr = DMLabelSetValue(faceLabel, ii, triLabels[ii-fStart][0]); CHKERRQ(ierr);
	  
	  /* Finish Labeling previously unlabeled DMPlex Edges - Assumes Triangular Cell (3 Edges Max) */
	  ierr = DMPlexGetCone(dmNozzle, ii, &cone); CHKERRQ(ierr);
	  
	  for (int jj = 0; jj < 3; ++jj) {
	  	ierr = DMLabelGetValue(edgeLabel, cone[jj], &edgeID_0); CHKERRQ(ierr);
	  
	  	if (edgeID_0 < 0) {
	  		ierr = DMLabelSetValue(bodyLabel, cone[jj], triLabels[ii-fStart][1]); CHKERRQ(ierr);
	  		ierr = DMLabelSetValue(faceLabel, cone[jj], triLabels[ii-fStart][0]); CHKERRQ(ierr);	
		}
	  }
  }

  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNozzle \n");CHKERRQ(ierr);
  ierr = DMView(dmNozzle, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  ierr = DMLabelView(bodyLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(faceLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(edgeLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = DMLabelView(vertexLabel, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);  
  
  ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr); 
  
  // Remove when not using Tetgen
  ierr = DMPlexGenerate(dmNozzle, "tetgen", PETSC_TRUE, &dmMesh); CHKERRQ(ierr);

  /* Output State of DMLabels for dmMesh after Volumetric Mesh generated */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  ierr = DMView(dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  /* Output Pre-refinement Volumetric Mesh */
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view3");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  surfArea(dmMesh);
  
  /* Refine Volumetric Mesh (dmMesh) */
  // Petsc Refinement
  // 1st time
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 1st Refinement \n");CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view4");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  surfArea(dmMesh);
  
  // 2nd Time
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 2nd Refinement \n");CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view5");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  surfArea(dmMesh);
  
  // 3rd Time
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh Created Trying 3rd Refinement \n");CHKERRQ(ierr);
  ierr = DMSetFromOptions(dmMesh);CHKERRQ(ierr);    // Check Snap_to_Geometry on Volumetric Mesh
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view6");CHKERRQ(ierr);
  
  // Calculate Surface Area
  //BsurfArea(dmNozzle);
  surfArea(dmMesh);
  
  /* Destry DMPlexes to free memory */
  //ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  ierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);		// Strange error

  /* Close EGADSlite file */
  ierr = EG_close(context);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  test:
    suffix: sphere_0
    args: -filename ${wPETSC_DIR}/share/petsc/datafiles/meshes/unit_sphere.egadslite -dm_view ::ascii_info_detail

TEST*/


void surfArea(DM dm) {
  DMLabel        bodyLabel, faceLabel;
  double         surfaceArea = 0., volume = 0.;
  PetscReal      vol, centroid[3], normal[3];
  PetscInt       dim, cStart, cEnd, fStart, fEnd;
  PetscInt       bodyID, faceID;
  PetscErrorCode ierr;
  
  ierr = DMGetDimension(dm, &dim); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "    dim = %d \n", dim);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Body ID",   &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID",   &faceLabel);CHKERRQ(ierr);
  
  if ( dim == 2 ) {
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
    for (int ii = cStart; ii < cEnd; ++ii) {
		ierr = DMLabelGetValue(faceLabel, ii, &faceID);CHKERRQ(ierr);
		if ( faceID >= 0) {
	      ierr = DMPlexComputeCellGeometryFVM(dm, ii, &vol, &centroid, &normal); CHKERRQ(ierr);
	      surfaceArea += vol;
		}
	}
  }

  if ( dim == 3 ) {
	  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
	  for (int ii = fStart; ii < fEnd; ++ii) {
		ierr = DMLabelGetValue(faceLabel, ii, &faceID);CHKERRQ(ierr);
		if ( faceID >= 0 ) {
		  ierr = DMPlexComputeCellGeometryFVM(dm, ii, &vol, &centroid, &normal); CHKERRQ(ierr);
		  surfaceArea += vol;
		}
	  }
	  
	  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
	  for (int ii = cStart; ii < cEnd; ++ii) {
		ierr = DMLabelGetValue(bodyLabel, ii, &bodyID);CHKERRQ(ierr);
		if ( bodyID >= 0 ) {
		  ierr = DMPlexComputeCellGeometryFVM(dm, ii, &vol, &centroid, &normal); CHKERRQ(ierr);
		  volume += vol;
		}
	  }
  }
  
  if ( dim == 2 ) {
	  ierr = PetscPrintf(PETSC_COMM_SELF, "    Surface Area = %.6e \n", surfaceArea);CHKERRQ(ierr);
  } else if ( dim == 3 ) {
	  ierr = PetscPrintf(PETSC_COMM_SELF, "    Volume = %.6e \n", volume);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_SELF, "    Surface Area = %.6e \n", surfaceArea);CHKERRQ(ierr);  
  } else {
	  // Do Nothing
  }	
}