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

$(ODIR)\limits.obj:	limitTessBody.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h $(SDIR)\..\src\egadsTris.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(SDIR)\..\src -DSTANDALONE \
		limitTessBody.c /Fo$(ODIR)\limits.obj

$(ODIR)\limitTessBody.obj:	limitTessBody.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
		$(IDIR)\egadsErrors.h $(SDIR)\..\src\egadsTris.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(SDIR)\..\src limitTessBody.c \
		/Fo$(ODIR)\limitTessBody.obj

clean:
	-del $(ODIR)\limitTessBody.obj $(ODIR)\limits.obj 

cleanall:	clean
	-del $(TDIR)\limits.exe
