#EGADS_LIBDIR = /PETSc3/cad/EGADSlite/src
#EGADS_LIBDIR2= ../src
#EGADS_LIBDIR2=$ESP_ROOT/lib
#CFLAGS = -I../include

# Note: I replaced -legadslite with -legads

testPlexEx45: src_ts_tutorials_ex45_EGADS.o
	${CLINKER} -o $@ $< -legadslite ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -L$(EGADS_LIBDIR2)

clean::
	${RM} testPlexEx45 src_ts_tutorials_ex45_EGADS.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
