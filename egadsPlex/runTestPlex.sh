#!/bin/bash

# Compile Code
make clean
make

# Run code with all options
./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh.h5
./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r1.h5 -dm_refine 1
./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r5.h5 -dm_refine 5
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r10.h5 -dm_refine 10
#mpiexec -n 32 ./testPlex -filename ../examples/unit_sphere.egadslite -dm_view hdf5:mesh_r25.h5 -dm_refine 25

