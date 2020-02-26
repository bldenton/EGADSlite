#!/bin/bash

# Unload petsc module if loaded
module unload petsc

# Go to build directory
cd $UBCESLAB_SWENV_PREFIX/builddir

# Remove current unpacked folder
rm -rf petsc-3.12.4

# Unpack tarball
tar -xzf $UBCESLAB_SWENV_PREFIX/sourcesdir/petsc/petsc-lite-3.12.4.tar.gz

# Petsc v3.12.4 doesn't need the patch. Only Petsc v 3.12.0 does
## Copy patch for hdf5
#cp /kif1/data/shared/software/patch/hdf5.patch /kif1/data/shared/software/builddir/petsc-3.12.0/hdf5.patch
#
## Apply the patch
#cd /kif1/data/shared/software/builddir/petsc-3.12.0
#patch -p1 < hdf5.patch

# Copy Modified plexrefine
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexrefine.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexrefine.c

# Copy Modified petscvariables
#cp /kif1/data/shared/software/egadsMOD/petscvariables /kif1/data/shared/software/builddir/petsc-3.12.0/lib/petsc/conf/petscvariables

# Copy egads.py
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/egads.py $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/config/BuildSystem/config/packages

# Copy makefile with EGADS library
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/makefile $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/makefile

# Copy plexegads.c to correct location
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c

# Install the modified petsc code 
cd ~/ubaceslab/build_scripts
./build_petsc_egads_version.sh

