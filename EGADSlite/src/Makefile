#
include ../include/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
IDIR  = $(ESP_ROOT)/include
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

# allow for clang static analyzer
ifneq (,$(findstring analyzer,$(CC)))
DEFINE += -DDEBUG
endif
ifneq (,$(findstring intercept,$(CC)))
DEFINE += -DDEBUG
endif

VPATH = $(ODIR)

OBJS  = liteBase.o liteMemory.o liteGeom.o liteTopo.o liteAttrs.o liteImport.o \
	liteString.o

ifeq ($(ESP_ARCH),DARWIN64)
default:	$(LDIR)/libegadslite.dylib $(LDIR)/libegadslitestatic.a
else
default:	$(LDIR)/libegadslite.so    $(LDIR)/libegadslitestatic.a
endif

$(LDIR)/libegadslite.so:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslite.so
	rm $(LDIR)/libegadslite.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined \
		-o $(LDIR)/libegadslite.so $(OBJS) liteTess.o liteTris.o \
		liteQuads.o liteTessInp.o egadsRobust.o emp.o evaLite.o \
		ratLite.o liteRegQuads.o -lpthread -lm )

$(LDIR)/libegadslite.dylib:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslite.dylib
	rm $(LDIR)/libegadslite.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegadslite.dylib \
                $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o \
		egadsRobust.o emp.o evaLite.o ratLite.o liteRegQuads.o \
		-undefined error -install_name '@rpath/libegadslite.dylib' \
                -current_version $(EGREV) )

$(LDIR)/libegadslitestatic.a:	$(OBJS) liteTess.o liteTris.o liteQuads.o \
				liteTessInp.o egadsRobust.o emp.o evaLite.o \
				ratLite.o liteRegQuads.o
	touch $(LDIR)/libegadslitestatic.a
	rm $(LDIR)/libegadslitestatic.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libegadslitestatic.a $(OBJS) \
                liteTess.o liteTris.o liteQuads.o liteTessInp.o egadsRobust.o \
		emp.o evaLite.o ratLite.o liteRegQuads.o ; $(RANLB) )

$(OBJS): %.o:	%.c ../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
		../include/egadsTypes_lite.h ../include/liteClasses.h ../include/liteDevice.h
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I. -I../src $< \
		-o $(ODIR)/$@

liteTess.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTess_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTess_lite.c -o $(ODIR)/liteTess.o

liteTris.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTris_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                ../src/egadsTris_lite.c -o $(ODIR)/liteTris.o

liteQuads.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsQuads_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -I../include -I. -I../src \
                ../src/egadsQuads_lite.c -o $(ODIR)/liteQuads.o

liteTessInp.o:	../include/egadsErrors_lite.h ../include/egadsInternals_lite.h \
                ../include/egadsTypes_lite.h ../include/liteClasses.h ../src/egadsTessInp_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
                -I../src ../src/egadsTessInp_lite.c -o $(ODIR)/liteTessInp.o

emp.o:		../src/emp_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../src/emp_lite.c \
		-o $(ODIR)/emp.o

evaLite.o:	../src/evaluate_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I. -I../src \
		../src/evaluate_lite.c -o $(ODIR)/evaLite.o

ratLite.o:	../src/rational_lite.c
	$(CC) -c $(COPTS) $(DEFINE) -DLITE ../src/rational_lite.c \
		-o $(ODIR)/ratLite.o

liteRegQuads.o:	../src/regQuads_lite.c ../include/regQuads_lite.h
	$(CC) -c $(COPTS) $(DEFINE) -DLITE -I../include -I../src \
		../src/regQuads_lite.c -o $(ODIR)/liteRegQuads.o

egadsRobust.o:	../src/egadsRobust_lite.c
	$(CC) -c $(COPTS) $(DEFINE) ../src/egadsRobust_lite.c \
		-o $(ODIR)/egadsRobust.o

# Static analysis
SCANDIR=/tmp/scanEGADSlite
SCANEXCLUDE=
include ../include/STANALYZER.make

lint:
	$(LINT) -I$(IDIR) -DLITE -I../src -I. liteAttrs.c liteGeom.c liteMemory.c liteBase.c liteImport.c liteString.c liteTest.c liteTopo.c ../src/evaluate_lite.c ../src/rational_lite.c ../src/regQuads_lite.c -evalorder -exportlocal -shiftimplementation -aliasunique -retvalother -uniondef -formatcode

clean:
	-(cd $(ODIR); rm $(OBJS) liteTess.o liteTris.o liteQuads.o liteTessInp.o evaLite.o ratLite.o liteRegQuads.o )

cleanall:	clean
	touch $(LDIR)/libegadslite.dylib $(LDIR)/libegadslite.so
	-rm $(LDIR)/libegadslite*
