# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015
#
# Description : Makefile for Ocean Circulation Model
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

TOPDIR = $(shell pwd)

include $(CONFFILE)

TARGET = OCEANCIRC
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

SRCDIRs = $(TOPDIR)/SwanSrc

FSOURCES = Circulation.f90 OceanModel.f90 WaveModel.f90

FOBJS=$(FSOURCES:.f90=.o)

.PHONY: dist module dust clobber

dist : module

module :
	@echo "*************************************************"
	@echo
	@echo "Build OceanCirc modules and library."
	@echo
	@echo "*************************************************"
	for i in $(SRCDIRs) ; do   \
    	( cd $$i ; make dist) ;       \
	done
	make $(TARGETLIB)

$(TARGETLIB): $(FOBJS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(BADLANDS_F) -c $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) $< -o $@

dust:
	rm -f $(SRCDIRs)/*.o $(SRCDIRs)/*.mod $(SRCDIRs)/*.original
	rm -f $(SRCDIRs)/*~ $(SRCDIRs)/*_genmod*
	rm -fv *.bak *~ *.o *.mod *.original

clobber: dust
