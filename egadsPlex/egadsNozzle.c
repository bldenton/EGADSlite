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
  DMLabel        bodyLabel, faceLabel, edgeLabel;
  PetscInt       cStart, cEnd, c;
  /* EGADSLite variables */
  ego            context, model, geom, *bodies, *objs, *nobjs, *mobjs, *lobjs;
  int            oclass, mtype, nbodies, *senses;
  int            b;
  /* PETSc variables */
  DM             dm;
  PetscInt       dim = -1, cdim = -1, numCorners = 0, numVertices = 0, numCells = 0;
  PetscInt      *cells  = NULL;
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
    ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbodies, &bodies, &senses);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF, " Number of BODIES (nbodies): %d \n", nbodies);CHKERRQ(ierr);

    for (b = 0; b < nbodies; ++b) {
      ego body = bodies[b];
      int id, Nsh, Nf, Nl, l, Ne, e, Nv, v;

      /* Output Basic Model Topology */
      ierr = EG_getBodyTopos(body, NULL, SHELL, &Nsh, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of SHELLS: %d \n", Nsh);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, FACE,  &Nf, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of FACES: %d \n", Nf);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, LOOP,  &Nl, &lobjs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of LOOPS: %d \n", Nl);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, EDGE,  &Ne, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of EDGES: %d \n", Ne);CHKERRQ(ierr);

      ierr = EG_getBodyTopos(body, NULL, NODE,  &Nv, &objs);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_SELF, "   Number of NODES: %d \n", Nv);CHKERRQ(ierr);

      for (l = 0; l < Nl; ++l) {
        ego loop = lobjs[l];
\
        id   = EG_indexBodyTopo(body, loop);
        ierr = PetscPrintf(PETSC_COMM_SELF, "          LOOP ID: %d\n", id);CHKERRQ(ierr);

        /* Get EDGE info which associated with the current LOOP */
        ierr = EG_getTopology(loop, &geom, &oclass, &mtype, NULL, &Ne, &objs, &senses);CHKERRQ(ierr);

        for (e = 0; e < Ne; ++e) {
          ego edge = objs[e];

          id = EG_indexBodyTopo(body, edge);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_SELF, "            EDGE ID: %d\n", id);CHKERRQ(ierr);

          double range[4] = {0., 0., 0., 0.};
          double point[3] = {0., 0., 0.};
          int    peri;

          ierr = EG_getRange(objs[e], range, &peri);
          ierr = PetscPrintf(PETSC_COMM_SELF, "            Range = %lf, %lf, %lf, %lf \n", range[0], range[1], range[2], range[3]);

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

            point[0] = point[0] + limits[0];
            point[1] = point[1] + limits[1];
            point[2] = point[2] + limits[2];
          }

#if 0
          point[0] = point[0]/2.;
          point[1] = point[1]/2.;
          point[2] = point[2]/2.;

          double trange[2];

          trange[0] = 0.;
          trange[1] = 0.;
          double *params[4];
          double *result[3];
          double *xyzresult[9];
          double t=0.;

          ierr = EG_nearestOnCurve(objs[e], point, range, t, &xyzresult);
          ierr = PetscPrintf(PETSC_COMM_SELF, " (t1, t2) = (%lf, %lf) \n", params[0], params[1]);
          ierr = PetscPrintf(PETSC_COMM_SELF, " (x, y, z) = (%lf, %lf, %lf) \n", xyzresult[0], xyzresult[1], xyzresult[2]);
#endif
        }
      }
    }
  }

//  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
//
//  ierr = DMViewFromOptions(dm, NULL, "-dm_view");CHKERRQ(ierr);
//  ierr = DMDestroy(&dm);CHKERRQ(ierr);

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
