static const char help[] = "Test of EGADSLite CAD functionality";

#include <petscdmplex.h>

#include <egads.h>
#include <petsc.h>

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
  /* PETSc variables */
  DM             dm, dmNozzle, dmMesh;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0, depth = 0;
  PetscInt       numFaces = 0, numEdges = 0, numPoints = 0;
  PetscInt      *cells  = NULL, *coneOrient = NULL, *cones = NULL, *coneSize = NULL;
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
  for ( b = 0; b < nbodies; ++b) {
	  ego body = bodies[b];
	  double params[3] = {0.0, 0.0, 0.0};		// Parameters for Tessellation
	  
	  int f, Nf;
	  
	  ierr = EG_makeTessBody(body, params, &tess); 
	  ierr = PetscPrintf(PETSC_COMM_SELF, "ierr = %d \n", ierr);CHKERRQ(ierr);
	  CHKERRQ(ierr);
	  
	  ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &fobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);
	  
	  for ( f = 0; f < Nf; ++f) {
		  /* Get Face Object */
		  ego face = fobjs[f];
		  
		  int len, *ptype, *pindex, *ntris, *ptris, *ptric;
		  double *pxyz, *puv;
		  
		  /* Get Face ID */
		  int fID = EG_indexBodyTopo(body, face);
		  
		  /* Checkout the Surface Tessellation */
		  ierr = EG_getTessFace(tess, fID, &len, &pxyz, &puv, &ptype, &pindex, &ntris, &ptris, &ptric);
		  ierr = PetscPrintf(PETSC_COMM_SELF, "ierr = %d \n", ierr);CHKERRQ(ierr);
		  CHKERRQ(ierr);
		  
		  /* See what we have from the tessellation */
		  ierr = PetscPrintf(PETSC_COMM_SELF, "  (fID, len, ntris) = (%d, %d, %d) \n", fID, len, ntris);CHKERRQ(ierr);
		  
		  /* Check out the point index and coordinate */
		  for (int p = 0; p < len; ++p) {
			int global;
			//ierr = EG_localToGlobal(tess, fID, pindex[p], &global);
			ierr = PetscPrintf(PETSC_COMM_SELF, "  pIndex (type, local, global) :: (x, y, z) = (%d, %d, %d) :: (%lf, %lf, %lf) \n", ptype[p], pindex[p], global, pxyz[(p*3)+0], pxyz[(p*3)+1], pxyz[(p*3)+2]);CHKERRQ(ierr);
		  }
	  }
  }

  //ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  //ierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);		// Strange error

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
