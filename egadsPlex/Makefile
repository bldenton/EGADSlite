EGADS_LIBDIR = /PETSc3/cad/EGADSlite/src
CFLAGS = -I../include

testPlex: egadsPlex.o
	${CLINKER} -o $@ $< ${PETSC_LIB} ${CC_LINKER_SLFLAG}${EGADS_LIBDIR} -L${EGADS_LIBDIR} -legadslite

clean::
	${RM} testPlex egadsPlex.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
