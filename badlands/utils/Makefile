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
FFLAGS = -shared -O3 -fPIC

FMOD1 = classpd
FMOD2 = classfv
FMOD3 = classoro

all:
	${F90} ${FFLAGS} -c ${FMOD1}.f90
	${F90} ${FFLAGS} -c ${FMOD2}.f90
	${F90} ${FFLAGS} -c ${FMOD3}.f90

clean:
	rm -fv *~ *.bak *.o *.mod *.original *.so
	rm -rfv *.dSYM
