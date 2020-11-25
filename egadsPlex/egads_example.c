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
  /* PETSc variables */
  DM             dmNozzle = NULL, dmMesh = NULL;
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
  
  /* Set Plex from Options and Save to file*/
  ierr = DMSetFromOptions(dmNozzle); CHKERRQ(ierr);
  ierr = DMViewFromOptions(dmNozzle, NULL, "-dm_view");CHKERRQ(ierr);;
  
  /* Generate Volumetric Mesh using TETGEN - Remove when not using Tetgen */
  ierr = DMPlexGenerate(dmNozzle, "tetgen", PETSC_TRUE, &dmMesh); CHKERRQ(ierr);
  
  /* Output State of Plex after Volumetric Mesh generated */
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n dmMesh \n");CHKERRQ(ierr);
  ierr = DMView(dmMesh, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF, "\n");CHKERRQ(ierr);  
  
  /* Save Volumetric to file */
  ierr = DMViewFromOptions(dmMesh, NULL, "-dm_view2");CHKERRQ(ierr);
 
  /* Destroy DMPlexes before exiting program */
  ierr = DMDestroy(&dmMesh);CHKERRQ(ierr);
  ierr = DMDestroy(&dmNozzle);CHKERRQ(ierr);
  
  /* Close EGADSlite file */
  ierr = PetscFinalize();
  return ierr;
}
