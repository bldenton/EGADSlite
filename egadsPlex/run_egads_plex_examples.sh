#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make test_egads_examples

# Run code with all options
#  Sphere Mesh Creation
./test_egads_examples -filename ../examples/unit_sphere.egadslite -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_sphere_surf.h5 -dm_view2 hdf5:mesh_sphere_vol.h5

#  Cylinder Mesh Creation
./test_egads_examples -filename ../examples/cylinder_example.egadslite -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_cylinder_surf.h5 -dm_view2 hdf5:mesh_cylinder_vol.h5

#  Cone Mesh Creation
./test_egads_examples -filename ../examples/cone_example.egadslite -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_cone_surf.h5 -dm_view2 hdf5:mesh_cone_vol.h5

# Complex Development Files -- ADVANCED FUTURE CHECK
#./test_egads_examples -filename ../examples/complex.egadslite -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_complex.h5 -dm_view2 hdf5:mesh_vol_complex.h5 -dm_view3 hdf5:mesh_vol_complex.h5 -dm_view4 hdf5:mesh_vol_complex_r1.h5 -dm_view5 hdf5:mesh_vol_complex_r2.h5 -dm_view6 hdf5:mesh_vol_complex_r3.h5 -dm_view7 hdf5:mesh_vol_complex_r4.h5 -dm_view8 hdf5:mesh_vol_complex_r5.h5 -dm_refine 1

# Nozzle Development Files
./test_egads_examples -filename ../examples/Nozzle_fromStep.egadslite -dm_refine 0 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_nozzle_surf.h5 -dm_view2 hdf5:mesh_nozzle_vol.h5

# -- Save for reference to generate .xmf files from .h5 files -- */
#${PETSC_DIR}/lib/petsc/bin/petsc_gen_xdmf.py mesh_nozzle.h5

