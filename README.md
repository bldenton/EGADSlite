# EGADSlite
Project is intended to integrate the capabilities of EGADSlite into Petsc. 

When used with Petsc (option -with-egads = 1), Petsc will read a .egadslite CAD file and generate a 2D DM Plex with the CAD topology embedded.
The embedded information allows the Plex (mesh) to "snap to the geometry" the Plex is representing upon refinement. This technology allows the
Plex (mesh) to better approximate the domain and/or boundary conditions as it is refined.

The initially created 2D DM Plex can be used to generate a 3D Plex (mesh) via Petsc's implementation of cTetgen or Tetgen. In this case, the CAD
topology is transferred from the 2D Plex (mesh) to the 3D Plex (mesh) and gains the "snap to geometry" capability. In addition, boundary conditions
can be assigned via geometric association.

Additional Capabilies include:
     1) Inflate to Geometery :: When 3D Plexes (meshes) are generated via cTetgen or Tetgen with Steiner points, it must be "inflated" so the
                                Steiner points lie on the associated CAD surfaces.
     2) Print Model Topology :: Print Model Topology to screen
     
