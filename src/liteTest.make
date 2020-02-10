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


$(TDIR)/liteTest:	$(ODIR)/liteTest.o
	$(CC) -o $(TDIR)/liteTest $(ODIR)/liteTest.o -L$(LDIR) -legadslite \
		$(RPATH) -lpthread -lm

$(ODIR)/liteTest.o:	liteTest.c
	$(CC) -c $(COPTS) $(DEFINE) -I$(IDIR) liteTest.c -o $(ODIR)/liteTest.o

clean:
	-rm $(ODIR)/liteTest.o

cleanall:	clean
	-rm $(TDIR)/liteTest
