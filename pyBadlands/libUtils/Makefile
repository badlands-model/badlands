##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##

## Makefile used to wrap C/Fortran functions used in Badlands

.PHONY: all clean

DIR=$(shell pwd)

F90 = gfortran
F2PY = f2py
FFLAGS = -shared -O2 -fPIC
LIBDIR = .
INCDIR = .

# SWAN compilation
TOP=$(shell pwd)
SWAN = simswan
DIRMODS = WaveComp
BUILDDIR = $(TOP)/WaveComp/build
OBJDIR = $(BUILDDIR)/objects
MODDIR = $(BUILDDIR)/modules
LIBDIR2 = $(BUILDDIR)/lib
FMPI = mpif90
F2PYMPI = f2py --fcompiler=gfortran
# Fortran optimisation flags
FCFLAGS = -O3 -funroll-loops --param max-unroll-times=2 -cpp \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none

#FCFLAGS= -m64 -w -fno-second-underscore -shared  -fPIC -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0  -ffpe-trap=zero,overflow,underflow

FFLAGS2 = $(FCFLAGS) -I$(MODDIR)
F2PY_FLAGS = --f90exec=${FMPI} -I$(MODDIR)
LDFLAGS = -L$(LIBDIR2)/

# Depression filling algorithm using Planchon and Darboux
PD_FMOD = PDclass
PD_LMOD = PDalgo

# Flow parameters computation
FLOW_LMOD = FLOWalgo

# Fast loop computation
LOOP_LMOD = FASTloop

# Flow network stack computation
FNET_FMOD = FLOWstack
FNET_LMOD = FLWnetwork

# Finite Volume discretisation computation
FV_FMOD = FVclass
FV_LMOD = FVframe

# Orographic precipitation computation
OR_FMOD = ORclass
OR_LMOD = ORmodel

all: waves ${PD_FMOD}.so ${PD_LMOD}.so ${FLOW_LMOD}.so ${LOOP_LMOD}.so ${FNET_FMOD}.so ${FNET_LMOD}.so ${FV_FMOD}.so ${FV_LMOD}.so ${OR_FMOD}.so ${OR_LMOD}.so sfd.so

	@echo ""
	@echo ""
	@echo "---------------------------------------------"
	@echo ""
	@echo ""
	@echo "Set the library path in your ~/.bashrc file"
	@echo "Add the following to LD_LIBRARY_PATH and PYTHON_PATH: "
	@echo ${DIR}
	@echo ""
	@echo "Alternatively use a symbolic link"
	@echo "ln -s "${DIR}" ~/.python"
	@echo ""
	@echo ""
	@echo "---------------------------------------------"
	@echo ""
	@echo ""

waves:
	@echo "---------------------------------------------"
	@echo ""
	@echo ""
	@echo "Ocean Plugin Makefile"
	@echo "Build Swan libraries: "
	for i in $(DIRMODS) ; do   \
				( cd $$i ; make dist) ;       \
	done
	@$(if $(wildcard SimSwan.so),rm -f *.so,)
	make libs
	@echo "*************************************************"
	@echo "Build Badlands other fortran modules"
	@echo
	@echo "*************************************************"

libs :	${SWAN}.so

${SWAN}.so: ${SWAN}.f90
	${F2PYMPI} ${F2PY_FLAGS} $(LDFLAGS) -lWaveParam -c -m ${SWAN} $<

%.o : %.f90
	$(FMPI) $(FFLAGS2) -c $< -o $@

${PD_FMOD}.so: ${PD_FMOD}.f90
	${F90} ${FFLAGS} -o lib${PD_FMOD}.so $<

${PD_LMOD}.so: ${PD_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -L${LIBDIR} -I${INCDIR} -l${PD_FMOD} -m ${PD_LMOD} $<

${FLOW_LMOD}.so: ${FLOW_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -m ${FLOW_LMOD} $<

${LOOP_LMOD}.so: ${LOOP_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -m ${LOOP_LMOD} $<

${FNET_FMOD}.so: ${FNET_FMOD}.f90
	${F90} ${FFLAGS} -o lib${FNET_FMOD}.so $<

${FNET_LMOD}.so: ${FNET_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -L${LIBDIR} -I${INCDIR} -l${FNET_FMOD} -m ${FNET_LMOD} $<

${FV_FMOD}.so: ${FV_FMOD}.f90
	${F90} ${FFLAGS} -o lib${FV_FMOD}.so $<

${FV_LMOD}.so: ${FV_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -L${LIBDIR} -I${INCDIR} -l${FV_FMOD} -m ${FV_LMOD} $<

${OR_FMOD}.so: ${OR_FMOD}.f90
	${F90} ${FFLAGS} -o lib${OR_FMOD}.so $<

${OR_LMOD}.so: ${OR_LMOD}.f90
	${F2PY} ${F2PY_FLAGS} -c -L${LIBDIR} -I${INCDIR} -l${OR_FMOD} -m ${OR_LMOD} $<

sfd.so: sfd.c sfd.pyf
	${F2PY} ${F2PY_FLAGS} -c -m sfd sfd.c sfd.pyf

clean:
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make clobber) ;   \
	done
	rm -fv *~ *.bak *.o *.mod *.original *.so
	rm -rfv $(BUILDDIR) *.dSYM
