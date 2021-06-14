#
ifndef ESP_ROOT
$(error ESP_ROOT must be set -- Please fix the environment...)
endif
#
include ../include/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
IDIR  = $(ESP_ROOT)/include
SIDIR = $(IDIR)/Surreal
ifdef ESP_BLOC
ODIR  = $(ESP_BLOC)/obj
else
ODIR  = .
endif

# allow for clang static analyzer
ifneq (,$(findstring analyzer,$(CC)))
CPPOPT += -UNo_Exception -UNDEBUG
DEFINE += -DDEBUG
endif
ifneq (,$(findstring intercept,$(CC)))
CPPOPT += -UNo_Exception -UNDEBUG
DEFINE += -DDEBUG
endif

VPATH = $(ODIR)

OBJSP = egadsGeom.o egadsHLevel.o    egadsTopo.o     egadsCopy.o     egadsSpline.o \
	egadsIO.o   egadsSplineFit.o egadsSkinning.o egadsTessSens.o BRepLib_FuseEdges.o
OBJS  = egadsBase.o egadsMemory.o egadsAttrs.o  egadsTess.o   egadsTessInp.o \
	egadsTris.o egadsQuads.o  egadsFit.o    egadsRobust.o egadsSBO.o \
	prmCfit.o   prmGrid.o     prmUV.o       egadsExport.o  \
	egadsEffect.o egadsSolids.o
FOBJS = fgadsBase.o fgadsMemory.o fgadsAttrs.o  fgadsTess.o   fgadsHLevel.o \
	fgadsGeom.o fgadsTopo.o
SRINC = $(SIDIR)/SurrealD.h $(SIDIR)/SurrealD_Lazy.h $(SIDIR)/SurrealD_Trad.h \
        $(SIDIR)/SurrealS.h $(SIDIR)/SurrealS_Lazy.h $(SIDIR)/SurrealS_Trad.h \
	$(SIDIR)/always_inline.h

default:	includes uts $(LDIR)/$(SHLIB) \
		$(LDIR)/libegadstatic.a $(LDIR)/libfgads.a $(LDIR)/libemp.a
	@echo " *** EGADS Build Complete! ***"

$(LDIR)/libegads.so:	$(OBJSP) $(OBJS) emp.o evaluate.o rational.o \
			regQuads.o $(LDIR)/libuvmap.a
	touch $(LDIR)/libegads.so
	rm $(LDIR)/libegads.so
	(cd $(ODIR); $(CC) -shared -Wl,-no-undefined -o $(LDIR)/libegads.so \
		$(OBJSP) $(OBJS) emp.o evaluate.o rational.o regQuads.o \
		$(LIBPATH) -luvmap $(LIBS) $(CPPSLB) $(RPATH) )

$(LDIR)/libegads.dylib:	$(OBJSP) $(OBJS) emp.o evaluate.o rational.o \
			regQuads.o $(LDIR)/libuvmap.a
	touch $(LDIR)/libegads.dylib
	rm $(LDIR)/libegads.dylib
	(cd $(ODIR); $(CC) -dynamiclib -o $(LDIR)/libegads.dylib $(OBJSP) \
		$(OBJS) emp.o evaluate.o rational.o regQuads.o \
		$(RPATH) $(LIBPATH) -luvmap $(LIBS) $(CPPSLB) -undefined error \
		-install_name '@rpath/libegads.dylib' \
		-compatibility_version $(CASREV) \
		-current_version $(EGREV) )

# CYGWIN
$(LDIR)/libegads.dll:	$(OBJSP) $(OBJS) emp.o evaluate.o rational.o \
			regQuads.o $(LDIR)/libuvmap.a
	touch $(LDIR)/libegads.dll
	touch $(LDIR)/libegads.dll.a
	rm $(LDIR)/libegads.dll
	rm $(LDIR)/libegads.dll.a
	(cd $(ODIR); $(CC) -shared -o $(LDIR)/libegads.dll \
		-Wl,--out-implib,$(LDIR)/libegads.dll.a \
		$(OBJSP) $(OBJS) emp.o evaluate.o rational.o regQuads.o \
		$(LIBPATH) -luvmap $(LIBS) $(CPPSLB) )

$(LDIR)/libegadstatic.a:	$(OBJSP) $(OBJS) emp.o evaluate.o rational.o \
				regQuads.o $(LDIR)/libuvmap.a
	touch $(LDIR)/libegadstatic.a
	rm $(LDIR)/libegadstatic.a
	cp $(LDIR)/libuvmap.a $(LDIR)/libegadstatic.a
	(cd $(ODIR); \
		ar $(LOPTS) $(LDIR)/libegadstatic.a $(OBJS) $(OBJSP) emp.o \
		evaluate.o rational.o regQuads.o; $(RANLB) )

$(LDIR)/libuvmap.a:
	$(MAKE) -C ../util

$(LDIR)/libfgads.a:		$(FOBJS)
	touch $(LDIR)/libfgads.a
	rm $(LDIR)/libfgads.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libfgads.a $(FOBJS); $(RANLB) )

$(LDIR)/libemp.a:		emp.o
	touch $(LDIR)/libemp.a
	rm $(LDIR)/libemp.a
	(cd $(ODIR); ar $(LOPTS) $(LDIR)/libemp.a emp.o; $(RANLB) )

$(ODIR)/SurrealD1:	../util/SurrealD1_btest.cpp
	$(CXX) -o $(ODIR)/SurrealD1 -O -I. ../util/SurrealD1_btest.cpp
	$(ODIR)/SurrealD1

$(ODIR)/SurrealD4:	../util/SurrealD4_btest.cpp
	$(CXX) -o $(ODIR)/SurrealD4 -O -I. ../util/SurrealD4_btest.cpp
	$(ODIR)/SurrealD4

$(ODIR)/SurrealS1:	../util/SurrealS1_btest.cpp
	$(CXX) -o $(ODIR)/SurrealS1 -O -I. ../util/SurrealS1_btest.cpp
	$(ODIR)/SurrealS1

$(ODIR)/SurrealS4:	../util/SurrealS4_btest.cpp
	$(CXX) -o $(ODIR)/SurrealS4 -O -I. ../util/SurrealS4_btest.cpp
	$(ODIR)/SurrealS4

uts:	$(ODIR)/SurrealD1 $(ODIR)/SurrealD4 $(ODIR)/SurrealS1 $(ODIR)/SurrealS4
	$(MAKE) -C ../util

$(OBJS): %.o:	%.c ../include/egadsErrors.h egadsInternals.h egadsStack.h \
		../include/egadsTypes.h egadsTris.h
	$(CC) -c $(COPTS) $(DEFINE) $(XDEF) -I../include -I../util \
		$< -o $(ODIR)/$@

$(FOBJS): %.o:	%.c ../include/egadsErrors.h egadsInternals.h egadsStack.h \
		../include/egadsTypes.h 
	$(CC) -c $(COPTS) $(DEFINE) $(XDEF) -I../include $< -o $(ODIR)/$@

$(OBJSP): %.o:	%.cpp ../include/egadsErrors.h egadsOCC.h egadsInternals.h \
		../include/egadsTypes.h ../include/egads_dot.h \
		egadsClasses.h egadsSplineFit.h egadsStack.h Surreal/SurrealS.h
	$(CXX) -c $(CPPOPT) $(DEFINE) $(XDEF) $(INCS) -I../include \
		-I. $< -o $(ODIR)/$@

emp.o:		../util/emp.c
	$(CC) -c $(COPTS) $(DEFINE) -D$(ESP_ARCH) -I../include ../util/emp.c \
		-o $(ODIR)/emp.o

evaluate.o:	../util/evaluate.c
	$(CXX) -x c++ -c $(CPPOPT) $(DEFINE) $(INCS) -I. -I../include \
		../util/evaluate.c -o $(ODIR)/evaluate.o

rational.o:	../util/rational.c
	$(CXX) -x c++ -c $(CPPOPT) $(DEFINE) $(INCS) -I. -I../include \
		../util/rational.c -o $(ODIR)/rational.o

regQuads.o:	../util/regQuads.c ../util/regQuads.h
	$(CC) -c $(COPTS) $(DEFINE) $(INCS) -I../include -I../util \
		../util/regQuads.c -o $(ODIR)/regQuads.o

includes: 	$(IDIR)/$(ESP_ARCH) $(IDIR)/STANALYZER.make \
		$(IDIR)/egads.h $(IDIR)/egads_dot.h $(IDIR)/egads.inc \
		$(IDIR)/egadsErrors.h $(IDIR)/egadsTypes.h $(IDIR)/emp.h \
		$(IDIR)/prm.h $(IDIR)/wsserver.h $(IDIR)/wsserver.inc \
		$(IDIR)/wsss.h $(IDIR)/egadsTris.h $(IDIR)/egadsSplineFit.h \
		$(IDIR)/egadsSplineVels.h $(IDIR)/egadsf90.inc $(SRINC)

$(IDIR)/$(ESP_ARCH):	../include/$(ESP_ARCH)
	touch $(IDIR)/$(ESP_ARCH)
	rm $(IDIR)/$(ESP_ARCH)
	cp -p ../include/$(ESP_ARCH) $(IDIR)/$(ESP_ARCH)

$(IDIR)/STANALYZER.make:	../include/STANALYZER.make
	touch $(IDIR)/STANALYZER.make
	rm $(IDIR)/STANALYZER.make
	cp -p ../include/STANALYZER.make $(IDIR)/STANALYZER.make

$(IDIR)/egads.h:	../include/egads.h
	touch $(IDIR)/egads.h
	rm $(IDIR)/egads.h
	cp -p ../include/egads.h $(IDIR)/egads.h

$(IDIR)/egads_dot.h:	../include/egads_dot.h
	touch $(IDIR)/egads_dot.h
	rm $(IDIR)/egads_dot.h
	cp -p ../include/egads_dot.h $(IDIR)/egads_dot.h

$(IDIR)/egads.inc:	../include/egads.inc
	touch $(IDIR)/egads.inc
	rm $(IDIR)/egads.inc
	cp -p ../include/egads.inc $(IDIR)/egads.inc

$(IDIR)/egadsErrors.h:	../include/egadsErrors.h
	touch $(IDIR)/egadsErrors.h
	rm $(IDIR)/egadsErrors.h
	cp -p ../include/egadsErrors.h $(IDIR)/egadsErrors.h

$(IDIR)/egadsTypes.h:	../include/egadsTypes.h
	touch $(IDIR)/egadsTypes.h
	rm $(IDIR)/egadsTypes.h
	cp -p ../include/egadsTypes.h $(IDIR)/egadsTypes.h

$(IDIR)/emp.h:		../include/emp.h
	touch $(IDIR)/emp.h
	rm $(IDIR)/emp.h
	cp -p ../include/emp.h $(IDIR)/emp.h

$(IDIR)/prm.h:		../include/prm.h
	touch $(IDIR)/prm.h
	rm $(IDIR)/prm.h
	cp -p ../include/prm.h $(IDIR)/prm.h

$(IDIR)/wsserver.h:	../include/wsserver.h
	touch $(IDIR)/wsserver.h
	rm $(IDIR)/wsserver.h
	cp -p ../include/wsserver.h $(IDIR)/wsserver.h

$(IDIR)/wsserver.inc:	../include/wsserver.inc
	touch $(IDIR)/wsserver.inc
	rm $(IDIR)/wsserver.inc
	cp -p ../include/wsserver.inc $(IDIR)/wsserver.inc

$(IDIR)/wsss.h:		../include/wsss.h
	touch $(IDIR)/wsss.h
	rm $(IDIR)/wsss.h
	cp -p ../include/wsss.h $(IDIR)/wsss.h

$(IDIR)/egadsTris.h:		egadsTris.h
	touch $(IDIR)/egadsTris.h
	rm $(IDIR)/egadsTris.h
	cp -p egadsTris.h $(IDIR)/egadsTris.h

$(IDIR)/egadsSplineFit.h:	egadsSplineFit.h
	touch $(IDIR)/egadsSplineFit.h
	rm $(IDIR)/egadsSplineFit.h
	cp -p egadsSplineFit.h $(IDIR)/egadsSplineFit.h

$(IDIR)/egadsSplineVels.h:	egadsSplineVels.h
	touch $(IDIR)/egadsSplineVels.h
	rm $(IDIR)/egadsSplineVels.h
	cp -p egadsSplineVels.h $(IDIR)/egadsSplineVels.h

$(IDIR)/egadsf90.inc:	../include/egadsf90.inc
	touch $(IDIR)/egadsf90.inc
	rm $(IDIR)/egadsf90.inc
	cp -p ../include/egadsf90.inc $(IDIR)/egadsf90.inc

$(SRINC):	Surreal/SurrealD.h Surreal/SurrealS.h Surreal/always_inline.h \
		Surreal/SurrealD_Lazy.h Surreal/SurrealD_Trad.h \
		Surreal/SurrealS_Lazy.h Surreal/SurrealS_Trad.h 
	mkdir -p $(SIDIR)
	cp Surreal/$(notdir $@) $@

# Static analysis
SCANDIR=/tmp/scanEGADS
SCANEXCLUDE=
include ../include/STANALYZER.make

lint:
	@echo EGADS library:
	$(LINT) -I../include -I../util egadsBase.c egadsAttrs.c egadsMemory.c egadsTess.c egadsTessInp.c egadsTris.c egadsQuads.c egadsFit.c egadsRobust.c egadsSBO.c ../util/regQuads.c egadsExport.c egadsSolids.c egadsEffect.c -evalorder -exportlocal -incompletetype -boolops -predboolothers -dependenttrans -nestcomment
	@echo " "
	@echo FORTRAN bindings:
	$(LINT) -I../include fgadsBase.c fgadsTopo.c fgadsGeom.c fgadsAttrs.c fgadsMemory.c fgadsTess.c fgadsHLevel.c

clean:
	-(cd ../util; make clean)
	-(cd $(ODIR); rm $(FOBJS) $(OBJSP) $(OBJS) emp.o evaluate.o rational.o \
		regQuads.o SurrealD1 SurrealD4 SurrealS1 SurrealS4)

cleanall:	clean
	-(cd ../util; make cleanall)
	-rm $(LDIR)/libegads* $(LDIR)/libfgads.a $(LDIR)/libemp.a
	-rm $(IDIR)/$(ESP_ARCH) $(IDIR)/egads.h $(IDIR)/egads_dot.h
	-rm $(IDIR)/egadsErrors.h $(IDIR)/egadsTypes.h $(IDIR)/emp.h 
	-rm $(IDIR)/prm.h $(IDIR)/wsserver.h $(IDIR)/wsserver.inc $(IDIR)/wsss.h
	-rm $(IDIR)/egadsTris.h $(IDIR)/egadsSplineFit.h $(IDIR)/egads.inc
	-rm $(IDIR)/egadsSplineVels.h $(IDIR)/egadsf90.inc $(SRINC)
	-rmdir $(SIDIR)
