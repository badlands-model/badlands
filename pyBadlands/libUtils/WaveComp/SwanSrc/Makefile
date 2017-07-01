# Fortran optimisation flags
F77FLAGS = -O3 -m64 -w -fno-second-underscore -shared  -fPIC
FCFLAGS=-O3 -funroll-loops --param max-unroll-times=2 -cpp \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none \
		-shared -fPIC


BUILDDIR=../build
OBJDIR=$(BUILDDIR)/objects
MODDIR=$(BUILDDIR)/modules
LIBDIR=$(BUILDDIR)/lib

# FC compiler to use
F90=mpif90

FFLAGS= $(FCFLAGS) -I$(MODDIR)

LDFLAGS= -L$(LIBDIR)/



TOPDIR=$(TOP)

SWAN_OBJS = swmod1.o swmod2.o m_constants.o m_fileio.o serv_xnl4v5.o \
	mod_xnl4v5.o SwanGriddata.o SwanGridobjects.o SwanCompdata.o SwanCoupler.o SwanpyReef.o \
	swanpre1.o swanpre2.o swancom1.o swancom2.o swancom3.o swancom4.o \
	swancom5.o swanout1.o swanout2.o swanser.o swanparll.o SwanReadGrid.o \
	SwanReadADCGrid.o SwanReadTriangleGrid.o SwanReadEasymeshGrid.o SwanInitCompGrid.o \
	SwanCheckGrid.o SwanCreateEdges.o SwanGridTopology.o SwanGridVert.o SwanGridCell.o \
	SwanGridFace.o SwanPrintGridInfo.o SwanFindPoint.o SwanPointinMesh.o SwanBpntlist.o \
	SwanPrepComp.o SwanVertlist.o SwanCompUnstruc.o SwanDispParm.o SwanPropvelX.o \
	SwanSweepSel.o SwanPropvelS.o SwanTranspAc.o SwanTranspX.o SwanDiffPar.o \
	SwanGSECorr.o SwanInterpolatePoint.o SwanInterpolateAc.o SwanInterpolateOutput.o \
	SwanConvAccur.o SwanConvStopc.o SwanThreadBounds.o SwanFindObstacles.o \
	SwanCrossObstacle.o SwanComputeForce.o SwanIntgratSpc.o SwanBndStruc.o \
	SwanReadfort18.o SwanSumOverNodes.o SwanMinOverNodes.o SwanMaxOverNodes.o \
	ocpids.o ocpcre.o ocpmix.o

HCAT_OBJS = swanhcat.o

.SUFFIXES: .f .f90

.PHONY: dist

dist:
	make swan
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		cp $$f $(MODDIR);       \
		fi; 		\
	done

swan: $(SWAN_OBJS)

.f.o:
	$(F90) $< -c $(F77FLAGS)

.f90.o:
	$(F90) $< -c $(FFLAGS)
