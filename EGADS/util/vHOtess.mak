#
IDIR  = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
SDIR  = $(MAKEDIR)
UDIR  = $(SDIR)\..\util
LDIR  = $(ESP_ROOT)\lib
!IFDEF ESP_BLOC
ODIR  = $(ESP_BLOC)\obj
TDIR  = $(ESP_BLOC)\test
!ELSE
ODIR  = .
TDIR  = $(ESP_ROOT)\bin
!ENDIF

$(TDIR)\vHOtess.exe:	$(ODIR)\vHOtess.obj $(ODIR)\egadsHOtess.obj \
			$(LDIR)\egads.lib $(LDIR)\wsserver.lib $(LDIR)\z.lib
	cl /Fe$(TDIR)\vHOtess.exe $(ODIR)\vHOtess.obj $(ODIR)\egadsHOtess.obj \
		/link /LIBPATH:$(LDIR) wsserver.lib z.lib egads.lib ws2_32.lib
	$(MCOMP) /manifest $(TDIR)\vHOtess.exe.manifest \
		/outputresource:$(TDIR)\vHOtess.exe;1

$(ODIR)\vHOtess.obj:	vHOtess.c $(IDIR)\egads.h $(IDIR)\egadsTypes.h \
			$(IDIR)\egadsErrors.h $(IDIR)\wsss.h $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers vHOtess.c \
		/Fo$(ODIR)\vHOtess.obj

$(ODIR)\egadsHOtess.obj:	$(UDIR)\egadsHOtess.c $(IDIR)\egads.h \
				$(IDIR)\egadsTypes.h $(IDIR)\egadsErrors.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) $(UDIR)\egadsHOtess.c \
		/Fo$(ODIR)\egadsHOtess.obj

clean:
	-del $(ODIR)\egadsHOtess.obj $(ODIR)\vHOtess.obj

cleanall:	clean
	-del $(TDIR)\vHOtess.exe $(TDIR)\vHOtess.exe.manifest
