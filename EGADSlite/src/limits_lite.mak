#
!include ..\include\$(ESP_ARCH).$(MSVC)
SDIR = $(MAKEDIR)
IDIR = $(ESP_ROOT)\include
LDIR = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!ELSE
ODIR = .
TDIR = $(ESP_BLOC)\bin
!ENDIF

default:	$(TDIR)\limits.exe $(ODIR)\limitTessBody.obj

$(TDIR)\limits.exe:	$(ODIR)\limits.obj $(LDIR)\egads.lib
	cl /Fe$(TDIR)\limits.exe $(ODIR)\limits.obj $(LIBPTH) egads.lib
	$(MCOMP) /manifest $(TDIR)\limits.exe.manifest \
		/outputresource:$(TDIR)\limits.exe;1

$(ODIR)\limits.obj:	limitTessBody_lite.c $(IDIR)\egads_lite.h $(IDIR)\egadsTypes_lite.h \
		$(IDIR)\egadsErrors_lite.h $(SDIR)\..\src\egadsTris_lite.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(SDIR)\..\src -DSTANDALONE \
		limitTessBody_lite.c /Fo$(ODIR)\limits.obj

$(ODIR)\limitTessBody.obj:	limitTessBody_lite.c $(IDIR)\egads_lite.h $(IDIR)\egadsTypes_lite.h \
		$(IDIR)\egadsErrors_lite.h $(SDIR)\..\src\egadsTris_lite.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(SDIR)\..\src limitTessBody_lite.c \
		/Fo$(ODIR)\limitTessBody.obj

clean:
	-del $(ODIR)\limitTessBody.obj $(ODIR)\limits.obj 

cleanall:	clean
	-del $(TDIR)\limits.exe
