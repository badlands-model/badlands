# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015
#
# Description : Makefile for Ice Sheet Model
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

TARGET = ICE
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = IceModel.f90 

FOBJS=$(FSOURCES:.f90=.o)

.PHONY: dist module dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build IceSheet modules and library."
	@echo	
	@echo "*************************************************"
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
	rm -fv *.bak *~ *.o *.mod *.original
	
clobber: dust
