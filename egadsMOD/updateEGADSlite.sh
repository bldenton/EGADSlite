#!/bin/bash
UPDATE_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/ESP_v1.18/EngSketchPad/src/EGADS
UPDATE_DIR_ESP=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/ESP_v1.18/EngSketchPad
DESTIN_DIR=/mnt/c/Users/Brandon/Documents/School/Dissertation/Software/EGADSlite

echo "Updating EGADSlite version for Petsc"
echo ""

# Update include folder
cp $UPDATE_DIR/include/egads.h $DESTIN_DIR/include/egads.h
cp $UPDATE_DIR/include/egadsErrors.h $DESTIN_DIR/include/egadsErrors.h
cp $UPDATE_DIR/src/egadsInternals.h $DESTIN_DIR/include/egadsInternals.h
cp $UPDATE_DIR/src/egadsTris.h $DESTIN_DIR/include/egadsTris.h
cp $UPDATE_DIR_ESP/include/egadsTypes.h $DESTIN_DIR/include/egadsTypes.h
cp $UPDATE_DIR/include/emp.h $DESTIN_DIR/include/emp.h
cp $UPDATE_DIR/lite/liteClasses.h $DESTIN_DIR/include/liteClasses.h
cp $UPDATE_DIR/lite/cudaUtil.h $DESTIN_DIR/include/cudaUtil.h
cp $UPDATE_DIR/lite/liteDevice.h $DESTIN_DIR/include/liteDevice.h
cp $UPDATE_DIR/lite/liteString.h $DESTIN_DIR/include/liteString.h
cp $UPDATE_DIR/util/regQuads.h $DESTIN_DIR/include/regQuads.h

# Update src folder
cp $UPDATE_DIR/lite/Makefile $DESTIN_DIR/src/Makefile
cp $UPDATE_DIR/lite/NMakefile $DESTIN_DIR/src/NMakefile
cp $UPDATE_DIR/src/egadsQuads.c $DESTIN_DIR/src/egadsQuads.c
cp $UPDATE_DIR/src/egadsRobust.c $DESTIN_DIR/src/egadsRobust.c
cp $UPDATE_DIR/src/egadsTess.c $DESTIN_DIR/src/egadsTess.c
cp $UPDATE_DIR/src/egadsTessInp.c $DESTIN_DIR/src/egadsTessInp.c
cp $UPDATE_DIR/src/egadsTris.c $DESTIN_DIR/src/egadsTris.c
cp $UPDATE_DIR/lite/egadslite.def $DESTIN_DIR/src/egadslite.def
cp $UPDATE_DIR/lite/egadslite.rc $DESTIN_DIR/src/egadslite.rc
cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/src/emp.c
cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/src/evaluate.c
cp $UPDATE_DIR/lite/liteAttrs.c $DESTIN_DIR/src/liteAttrs.c
cp $UPDATE_DIR/lite/liteBase.c $DESTIN_DIR/src/liteBase.c
cp $UPDATE_DIR/lite/liteGeom.c $DESTIN_DIR/src/liteGeom.c
cp $UPDATE_DIR/lite/liteImport.c $DESTIN_DIR/src/liteImport.c
cp $UPDATE_DIR/lite/liteMemory.c $DESTIN_DIR/src/liteMemory.c
cp $UPDATE_DIR/lite/liteTess.mak $DESTIN_DIR/src/liteTess.mak
cp $UPDATE_DIR/lite/liteTess.make $DESTIN_DIR/src/liteTess.make
cp $UPDATE_DIR/lite/liteTest.c $DESTIN_DIR/src/liteTest.c
cp $UPDATE_DIR/lite/liteTest.mak $DESTIN_DIR/src/liteTest.mak
cp $UPDATE_DIR/lite/liteTest.make $DESTIN_DIR/src/liteTest.make
cp $UPDATE_DIR/lite/liteTopo.c $DESTIN_DIR/src/liteTopo.c
cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/src/rational.c
cp $UPDATE_DIR/lite/NR.make $DESTIN_DIR/src/NR.make
cp $UPDATE_DIR/lite/liteDevice.c $DESTIN_DIR/src/liteDevice.c
cp $UPDATE_DIR/lite/liteString.c $DESTIN_DIR/src/liteString.c
cp $UPDATE_DIR/lite/liteTessNR.make $DESTIN_DIR/src/liteTessNR.make

# Previously in util folder -- now moved to src folder for Petsc EGADSlite v1.18
cp $UPDATE_DIR/util/README.txt $DESTIN_DIR/src/README.txt
cp $UPDATE_DIR/util/SurrealD1_btest.cpp $DESTIN_DIR/src/SurrealD1_btest.cpp
cp $UPDATE_DIR/util/SurrealD4_btest.cpp $DESTIN_DIR/src/SurrealD4_btest.cpp
cp $UPDATE_DIR/util/SurrealS1_btest.cpp $DESTIN_DIR/src/SurrealS1_btest.cpp
cp $UPDATE_DIR/util/SurrealS4_btest.cpp $DESTIN_DIR/src/SurrealS4_btest.cpp
cp $UPDATE_DIR/util/ThreadTest.c $DESTIN_DIR/src/ThreadTest.c
cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/src/emp.c
cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/src/evaluate.c
cp $UPDATE_DIR/util/limitTessBody.c $DESTIN_DIR/src/limitTessBody.c
cp $UPDATE_DIR/util/limits.mak $DESTIN_DIR/src/limits.mak
cp $UPDATE_DIR/util/limits.make $DESTIN_DIR/src/limits.make
cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/src/rational.c
cp $UPDATE_DIR/util/retessFaces.c $DESTIN_DIR/src/retessFaces.c
cp $UPDATE_DIR/util/regQuads.h $DESTIN_DIR/src/regQuads.h
cp $UPDATE_DIR/util/regQuads.c $DESTIN_DIR/src/regQuads.c

# Update util folder - All files have been moved to either include or src or Petsc integration
# They continue to be updated here as backup
cp $UPDATE_DIR/util/README.txt $DESTIN_DIR/util/README.txt
cp $UPDATE_DIR/util/SurrealD1_btest.cpp $DESTIN_DIR/util/SurrealD1_btest.cpp
cp $UPDATE_DIR/util/SurrealD4_btest.cpp $DESTIN_DIR/util/SurrealD4_btest.cpp
cp $UPDATE_DIR/util/SurrealS1_btest.cpp $DESTIN_DIR/util/SurrealS1_btest.cpp
cp $UPDATE_DIR/util/SurrealS4_btest.cpp $DESTIN_DIR/util/SurrealS4_btest.cpp
cp $UPDATE_DIR/util/ThreadTest.c $DESTIN_DIR/util/ThreadTest.c
cp $UPDATE_DIR/util/emp.c $DESTIN_DIR/util/emp.c
cp $UPDATE_DIR/util/evaluate.c $DESTIN_DIR/util/evaluate.c
cp $UPDATE_DIR/util/limitTessBody.c $DESTIN_DIR/util/limitTessBody.c
cp $UPDATE_DIR/util/limits.mak $DESTIN_DIR/util/limits.mak
cp $UPDATE_DIR/util/limits.make $DESTIN_DIR/util/limits.make
cp $UPDATE_DIR/util/rational.c $DESTIN_DIR/util/rational.c
cp $UPDATE_DIR/util/retessFaces.c $DESTIN_DIR/util/retessFaces.c
cp $UPDATE_DIR/util/regQuads.h $DESTIN_DIR/util/regQuads.h
cp $UPDATE_DIR/util/regQuads.c $DESTIN_DIR/util/regQuads.c