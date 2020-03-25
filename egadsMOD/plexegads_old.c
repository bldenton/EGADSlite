#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/

#ifdef PETSC_HAVE_EGADS
#include <egads.h>
#endif

PetscErrorCode DMPlexSnapToGeomModel(DM dm, PetscInt p, const PetscScalar initCoords[], PetscScalar coords[])
{
  DMLabel        bodyLabel, faceLabel, edgeLabel;
  PetscContainer modelObj;
  PetscInt       cdim, d, bodyID, faceID, edgeID;
  ego           *bodies;
  ego            model, geom, body, face, edge;
  double         point[3], params[3], result[3];
  int            Nb, oclass, mtype, *senses;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDim(dm, &cdim);CHKERRQ(ierr);

#ifdef PETSC_HAVE_EGADS
  ierr = DMGetLabel(dm, "EGADS Body ID", &bodyLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Face ID", &faceLabel);CHKERRQ(ierr);
  ierr = DMGetLabel(dm, "EGADS Edge ID", &edgeLabel);CHKERRQ(ierr);
  ierr = DMLabelGetValue(bodyLabel, p, &bodyID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(faceLabel, p, &faceID);CHKERRQ(ierr);
  ierr = DMLabelGetValue(edgeLabel, p, &edgeID);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject) dm, "EGADS Model", (PetscObject *) &modelObj);CHKERRQ(ierr);
  ierr = PetscContainerGetPointer(modelObj, (void **) &model);CHKERRQ(ierr);
  ierr = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &Nb, &bodies, &senses);CHKERRQ(ierr);
  if (bodyID >= Nb) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Body %D is not in [0, %d)", bodyID, Nb);
  body = bodies[bodyID];
  for (d = 0; d < cdim; ++d) point[d] = initCoords[d];
  if (edgeID >= 0) {
    ierr = EG_objectBodyTopo(body, EDGE, edgeID, &edge);CHKERRQ(ierr);
    ierr = EG_invEvaluate(edge, point, params, result);
  } else {
    ierr = EG_objectBodyTopo(body, FACE, faceID, &face);CHKERRQ(ierr);
    ierr = EG_invEvaluate(face, point, params, result);
  }
  for (d = 0; d < cdim; ++d) coords[d] = result[d];
#else
  for (d = 0; d < cdim; ++d) coords[d] = initCoords[d];
#endif
  PetscFunctionReturn(0);
}
