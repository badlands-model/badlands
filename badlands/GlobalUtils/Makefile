# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015 
#
# Description : Makefile for GlobalUtils
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

FSTG = ./Utility

TARGET = GlobalUtils
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = $(FSTG)/StringUtility.f90 GlobalClass.f90 TopologyClass.f90 \
	Sorting.f90 HydrologyClass.f90 ForcesClass.f90 MeshPartition.f90

FOBJS=$(FSOURCES:.f90=.o)

.PHONY: dist module dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build GlobalUtils modules and library."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB): $(COBJS) $(FOBJS)
	cp *.o $(FSTG)/*.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(BADLANDS_F) -c $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) $< -o $@ 
	
dust:
	rm -fv *.bak *~ *.o *.mod *.original $(FSTG)/*.o $(FSTG)/*.original 
	
clobber: dust
