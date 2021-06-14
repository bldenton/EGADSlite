#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
ifdef ESP_BLOC
ODIR = $(ESP_BLOC)/obj
TDIR = $(ESP_BLOC)/test
else
ODIR = .
TDIR = $(ESP_ROOT)/bin
endif


$(TDIR)/vHOtess:	$(ODIR)/vHOtess.o $(ODIR)/egadsHOtess.o
	$(CXX) -o $(TDIR)/vHOtess $(ODIR)/vHOtess.o $(ODIR)/egadsHOtess.o \
		-L$(LDIR) -lwsserver -legads -lpthread -lz $(RPATH) -lm

$(ODIR)/vHOtess.o:	vHOtess.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) vHOtess.c -o $(ODIR)/vHOtess.o

$(ODIR)/egadsHOtess.o:	egadsHOtess.c $(IDIR)/egads.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) egadsHOtess.c \
		-o $(ODIR)/egadsHOtess.o

clean:
	-rm $(ODIR)/vHOtess.o $(ODIR)/egadsHOtess.o

cleanall:	clean
	-rm $(TDIR)/vHOtess
