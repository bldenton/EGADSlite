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

$(TDIR)/triServer:	$(ODIR)/triServer.o
	$(CXX) -o $(TDIR)/triServer $(ODIR)/triServer.o \
		-L$(LDIR) -lwsserver -lpthread -lz $(RPATH) -lm

$(ODIR)/triServer.o:	triServer.c $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) triServer.c -o $(ODIR)/triServer.o

clean:
	-rm $(ODIR)/triServer.o

cleanall:	clean
	-rm $(TDIR)/triServer
