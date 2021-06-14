#
IDIR = $(ESP_ROOT)\include
!include $(IDIR)\$(ESP_ARCH).$(MSVC)
LDIR = $(ESP_ROOT)\lib
!ifdef ESP_BLOC
ODIR = $(ESP_BLOC)\obj
TDIR = $(ESP_BLOC)\test
!else
ODIR = .
TDIR = $(ESP_ROOT)\bin
!endif

$(TDIR)\triServer.exe:	$(ODIR)\triServer.obj
	cl /Fe$(TDIR)\triServer.exe $(ODIR)\triServer.obj \
		$(LIBPTH) wsserver.lib z.lib ws2_32.lib
	$(MCOMP) /manifest $(TDIR)\triServer.exe.manifest \
		/outputresource:$(TDIR)\triServer.exe;1

$(ODIR)\triServer.obj:	triServer.c $(IDIR)\wsserver.h
	cl /c $(COPTS) $(DEFINE) -I$(IDIR) -I$(IDIR)\winhelpers triServer.c \
		/Fo$(ODIR)\triServer.obj

clean:
	-del $(ODIR)\triServer.obj

cleanall:	clean
	-del $(TDIR)\triServer.exe
