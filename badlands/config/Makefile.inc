CFGFILE = $(TOP)/config/Makefile.cfg.osx
include $(CFGFILE)

EXEC=bin/$(EXE)

BUILDDIR=$(TOP)/build
OBJDIR=$(BUILDDIR)/objects
MODDIR=$(BUILDDIR)/modules
LIBDIR=$(BUILDDIR)/lib

FFLAGS= $(FCFLAGS) -I$(MODDIR)
	
LDFLAGS= -L$(LIBDIR)/

