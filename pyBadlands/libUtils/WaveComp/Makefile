# Fortran optimisation flags
F77FLAGS = -O3 -m64 -w -fno-second-underscore

FCFLAGS=-O3 -funroll-loops --param max-unroll-times=2 -cpp \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none \
		-shared -fPIC

#F77FLAGS = -m64 -w -fno-second-underscore -shared  -fPIC -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -ffpe-trap=zero,overflow,underflow
#FCFLAGS= -m64 -w -fno-second-underscore -shared  -fPIC -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0  -ffpe-trap=zero,overflow,underflow

TOPDIR = $(shell pwd)

BUILDDIR=$(TOPDIR)/build
OBJDIR=$(BUILDDIR)/objects
MODDIR=$(BUILDDIR)/modules
LIBDIR=$(BUILDDIR)/lib

# FC compiler to use
F90=mpif90

FFLAGS= $(FCFLAGS) -I$(MODDIR)

TOPDIR = $(shell pwd)

TARGET = WaveParam
TARGETLIB = $(LIBDIR)/lib$(TARGET).so

SRCDIRs = $(TOPDIR)/SwanSrc
SOURCE1 = WaveClass.f90
SOURCE2 = SimWaves.f90

OBJS=./*.o

FOB1 = $(SOURCE1:.f90=.o)
FOB2 = $(SOURCE2:.f90=.o)

.PHONY: dist module dust clobber

dist :
	$(F90) WaveClass.f90 -c $(FFLAGS)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)
	@mkdir -p $(LIBDIR)
	cp *.o $(OBJDIR)
	cp *.mod $(MODDIR)
	@echo
	@echo "*************************************************"
	@echo "Build modules and library for Ocean Waves"
	@echo
	@echo "*************************************************"
	for i in $(SRCDIRs) ; do   \
    	( cd $$i ; make dist) ;       \
	done
	make module

module :
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
			if [ -f $$f ]; then \
				mv $$f $(MODDIR);       \
		fi; 		\
	done
	make $(TARGETLIB)

$(TARGETLIB):
	$(F90) SimWaves.f90 -c $(FFLAGS)
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
	ar -rcs $(TARGETLIB) $(OBJDIR)/*.o

%.o : %.f90
	$(F90) $(FFLAGS) -c $< -o $@

dust :
	rm -f $(SRCDIRs)/*.o $(SRCDIRs)/*.mod $(SRCDIRs)/*.original
	rm -f $(SRCDIRs)/*~ $(SRCDIRs)/*_genmod*
	rm -fv *.bak *~ *.o *.mod *.original *_genmod*

clobber : dust
	rm -rfv $(BUILDDIR)
