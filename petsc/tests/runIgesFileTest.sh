#!/bin/bash

# Compile Code
echo "---------------------------"
echo "     Recompiling Code"
echo "---------------------------"
make clean
make testNozzleIges

# Nozzle Mesh Creation using Universal Mesh Approach
./testNozzleIges -filename ../../examples/Nozzle_example.igs -dm_plex_egads_without_snap_to_geom 0 -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_uni_nozzle.h5 -dm_view2 hdf5:mesh_uni_vol_nozzle.h5 -dm_view3 hdf5:mesh_uni_vol_nozzle_postRefine.h5 -dm_view4 hdf5:mesh_uni_vol_nozzle_inflate.h5

# Nozzle Mesh Creation using Tessellation
#./testNozzleIges -filename ../examples/Nozzle_example.igs -dm_plex_egads_with_tess 1 -dm_refine 1 -dm_plex_egads_print_model 1 -dm_view hdf5:mesh_tess_nozzle.h5 -dm_view2 hdf5:mesh_tess_vol_nozzle.h5 -dm_view3 hdf5:mesh_tess_vol_nozzle_postRefine.h5 -dm_view4 hdf5:mesh_vol_r1_RR.h5 -dm_view5 hdf5:mesh_vol_r2_RR.h5 -dm_view6 hdf5:mesh_vol_r3_RR.h5 -dm_view7 hdf5:mesh_vol_r4_RR.h5 -dm_view8 hdf5:mesh_vol_r5_RR.h5


#${PETSC_DIR}/lib/petsc/bin/petsc_gen_xdmf.py mesh_nozzle.h5

