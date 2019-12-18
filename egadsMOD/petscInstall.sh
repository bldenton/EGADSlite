#!/bin/bash

# Unload petsc module if loaded
module unload petsc

# Go to build directory
cd /kif1/data/shared/software/builddir

# Remove current unpacked folder
rm -rf petsc-3.12.0

# Unpack tarball
tar -xzf $UBCESLAB_SWENV_PREFIX/sourcesdir/petsc/petsc-lite-3.12.0.tar.gz

# Copy patch for hdf5
cp /kif1/data/shared/software/patch/hdf5.patch /kif1/data/shared/software/builddir/petsc-3.12.0/hdf5.patch

# Apply the patch
cd /kif1/data/shared/software/builddir/petsc-3.12.0
patch -p1 < hdf5.patch

# Copy Modified plexrefine
cp /kif1/data/shared/software/egadsMOD/plexrefine.c /kif1/data/shared/software/builddir/petsc-3.12.0/src/dm/impls/plex/plexrefine.c

# Copy Modified petscvariables
#cp /kif1/data/shared/software/egadsMOD/petscvariables /kif1/data/shared/software/builddir/petsc-3.12.0/lib/petsc/conf/petscvariables

# Install the modified petsc code 
cd /eng/home/bldenton/ubaceslab/build_scripts
./build_petsc_egads_version.sh

