#!/bin/bash
echo "Updating plexegads.c in Petsc"
# Copy plexegads.c to correct location
cp /mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/egadsMOD/plexegads.c $UBCESLAB_SWENV_PREFIX/builddir/petsc-3.12.4/src/dm/impls/plex/plexegads.c

# Update Petsc
cd /mnt/c/Users/Brandon/software/builddir/petsc-3.12.4
make -f ./gmakefile
make install

cd $goEGADSlite
