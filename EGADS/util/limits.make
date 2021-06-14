#
include ../include/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif

default:	$(TDIR)/limits $(ODIR)/limitTessBody.o

$(TDIR)/limits:	$(ODIR)/limits.o $(LDIR)/$(SHLIB)
	$(CXX) -o $(TDIR)/limits $(ODIR)/limits.o -L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/limits.o:	limitTessBody.c ../src/egadsTris.h ../include/egads.h \
			../include/egadsTypes.h ../include/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I../src -DSTANDALONE \
		limitTessBody.c -o $(ODIR)/limits.o

$(ODIR)/limitTessBody.o:	limitTessBody.c ../src/egadsTris.h \
				../include/egads.h ../include/egadsTypes.h \
				../include/egadsErrors.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I../src limitTessBody.c \
		-o $(ODIR)/limitTessBody.o

clean:
	-rm $(ODIR)/limitTessBody.o $(ODIR)/limits.o

cleanall:	clean
	-rm $(TDIR)/limits
