# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015 
#
# Description : Makefile for SpmModel
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = SPM
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = Hydrology.f90 VisualiseSurface.f90 VisualiseDrainage.f90 WatershedDelineation.f90 \
	SpmHillslope.f90 Geomorphology.f90 

FOBJS=$(FSOURCES:.f90=.o)

.PHONY: dist module dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build SPM modules and library."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB): $(COBJS) $(FOBJS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(BADLANDS_F) -c $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) $< -o $@ 
	
dust:
	rm -fv *.bak *~ *.o *.mod *.original
	
clobber: dust
