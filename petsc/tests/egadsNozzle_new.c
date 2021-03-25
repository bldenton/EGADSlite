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
  /* PETSc variables */
  DM             dm, dmNozzle = NULL, dmMesh = NULL;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0, depth = 0;
  PetscInt       numFaces = 0, numEdges = 0, numPoints = 0;
  //PetscInt      *cells  = NULL, *coneOrient = NULL, *cones = NULL, *coneSize = NULL;
  //PetscReal     *coords = NULL, PARAMACC = 1.0e-4;
  MPI_Comm       comm;
  PetscMPIInt    rank;
  AppCtx         ctx;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = ProcessOptions(comm, &ctx);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  if (!rank) {
    ierr = DMPlexCreateEGADSFromFile(comm, ctx.filename, &dmNozzle); CHKERRQ(ierr);
  }

  /* View Current State of Plex */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmNozzle \n");CHKERRQ(ierr);
  ierr = DMView(dmNozzle, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  // Refines Surface Mesh per option -dm_refine
  ierr = DMSetFromOptions(dmNozzle); CHKERRQ(ierr);
  ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr);
  
  // Generate Volumetric Mesh using TETGEN - Remove when not using Tetgen
  ierr = DMPlexGenerate(dmNozzle, "tetgen", PETSC_TRUE, &dmMesh); CHKERRQ(ierr);  
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view2");CHKERRQ(ierr);
  
  /* Inflate Mesh to EGADS Geometry */
  ierr = DMPlexInflateToGeomModel(dmMesh); CHKERRQ(ierr);
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view4");CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF, "\n Inflated dmMesh \n");CHKERRQ(ierr);
 
  /* Output State of DMLabels for dmMesh after Volumetric Mesh generated */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  ierr = DMView(dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);
  
  // Moved above - ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view2");CHKERRQ(ierr);
  //Bierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view3");CHKERRQ(ierr);
  
  // Refine Volumetric Mesh per option -dm_refine
  ierr = DMSetFromOptions(dmMesh); CHKERRQ(ierr);
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view3");CHKERRQ(ierr); 
  
  // Calculate 2D :: Surface Area;  3D :: Volume & Surface Area
  //BsurfArea(dmNozzle);
  //surfArea(dmMesh);
 
  //ierr = DMDestroy(&dm);CHKERRQ(ierr);
  //ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  //Bierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);		// Strange error

  /* Close EGADSlite file */
  //ierr = EG_close(context);CHKERRQ(ierr);
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
