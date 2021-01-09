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

$(ODIR)/limits.o:	limitTessBody_lite.c ../src/egadsTris_lite.h ../include/egads_lite.h \
			../include/egadsTypes_lite.h ../include/egadsErrors_lite.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I../src -DSTANDALONE \
		limitTessBody_lite.c -o $(ODIR)/limits.o

$(ODIR)/limitTessBody.o:	limitTessBody_lite.c ../src/egadsTris_lite.h \
				../include/egads_lite.h ../include/egadsTypes_lite.h \
				../include/egadsErrors_lite.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I../src limitTessBody_lite.c \
		-o $(ODIR)/limitTessBody.o

clean:
	-rm $(ODIR)/limitTessBody.o $(ODIR)/limits.o

cleanall:	clean
	-rm $(TDIR)/limits
