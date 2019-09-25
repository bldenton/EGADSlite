LIBS = \
  -legadslite -L../src -I../include -lpetsc -L$PETSC_DIR/lib -I$PETSC_DIR/include 

CC = mpicc

all: testPlex

testPlex: egadPlex.c
	$(CC) $? $(LIBS) -o $@

clean:
	rm -f testPlex egadsPlex.o
