#!/bin/bash
UPDATE_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/ESP_v1.18/EngSketchPad/src/EGADS
DESTIN_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite/EGADS

### Update EGADS Source Code ###
echo ""
echo "[Updating EGADS version for Petsc]"

echo "    [Updating EGADS /include/ folder]"
cp -rf $UPDATE_DIR/include $DESTIN_DIR

echo "    [Updating EGADS /src/ folder]"
cp -rf $UPDATE_DIR/src $DESTIN_DIR

echo "    [Updating EGADS /util/ folder]"
cp -rf $UPDATE_DIR/util $DESTIN_DIR

echo "[COMPLETED EGADS Codebase Update]"
echo ""
