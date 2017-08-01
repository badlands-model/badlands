# Fortran optimisation flags
F77FLAGS = -O3 -m64 -w -fno-second-underscore -shared  -fPIC
FCFLAGS=-O3 -funroll-loops --param max-unroll-times=2 -cpp \
    -ffree-form -ftree-vectorize -ffast-math -lstdc++ -ffree-line-length-none \
		-shared -fPIC

#F77FLAGS = -m64 -w -fno-second-underscore -shared  -fPIC -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -ffpe-trap=zero,overflow,underflow
#FCFLAGS= -m64 -w -fno-second-underscore -shared  -fPIC -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0  -ffpe-trap=zero,overflow,underflow

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
	make mod
	@echo "*************************************************"
	@echo "Build SWAN fortran modules"
	@echo
	@echo "*************************************************"
	cp *.o $(OBJDIR)
	cp *.mod $(MODDIR)

mod:
	$(F90) swmod1.f -c $(F77FLAGS)
	$(F90) swmod2.f -c $(F77FLAGS)
	$(F90) m_constants.f90 -c $(FFLAGS)
	$(F90) m_fileio.f90 -c $(FFLAGS)
	$(F90) serv_xnl4v5.f90 -c $(FFLAGS)
	$(F90) mod_xnl4v5.f90 -c $(FFLAGS)
	$(F90) SwanGriddata.f90 -c $(FFLAGS)
	$(F90) SwanGridobjects.f90 -c $(FFLAGS)
	$(F90) SwanCompdata.f90 -c $(FFLAGS)
	$(F90) SwanCoupler.f90 -c $(FFLAGS)
	$(F90) SwanpyReef.f -c $(F77FLAGS)
	$(F90) swanpre1.f -c $(F77FLAGS)
	$(F90) swanpre2.f -c $(F77FLAGS)
	$(F90) swancom1.f -c $(F77FLAGS)
	$(F90) swancom2.f -c $(F77FLAGS)
	$(F90) swancom3.f -c $(F77FLAGS)
	$(F90) swancom4.f -c $(F77FLAGS)
	$(F90) swancom5.f -c $(F77FLAGS)
	$(F90) swanout1.f -c $(F77FLAGS)
	$(F90) swanout2.f -c $(F77FLAGS)
	$(F90) swanser.f -c $(F77FLAGS)
	$(F90) swanparll.f -c $(F77FLAGS)
	$(F90) SwanReadGrid.f90 -c $(FFLAGS)
	$(F90) SwanReadADCGrid.f90 -c $(FFLAGS)
	$(F90) SwanReadTriangleGrid.f90 -c $(FFLAGS)
	$(F90) SwanReadEasymeshGrid.f90 -c $(FFLAGS)
	$(F90) SwanInitCompGrid.f90 -c $(FFLAGS)
	$(F90) SwanCheckGrid.f90 -c $(FFLAGS)
	$(F90) SwanCreateEdges.f90 -c $(FFLAGS)
	$(F90) SwanGridTopology.f90 -c $(FFLAGS)
	$(F90) SwanGridVert.f90 -c $(FFLAGS)
	$(F90) SwanGridCell.f90 -c $(FFLAGS)
	$(F90) SwanGridFace.f90 -c $(FFLAGS)
	$(F90) SwanPrintGridInfo.f90 -c $(FFLAGS)
	$(F90) SwanFindPoint.f90 -c $(FFLAGS)
	$(F90) SwanPointinMesh.f90 -c $(FFLAGS)
	$(F90) SwanBpntlist.f90 -c $(FFLAGS)
	$(F90) SwanGridCell.f90 -c $(FFLAGS)
	$(F90) SwanPrepComp.f90 -c $(FFLAGS)
	$(F90) SwanVertlist.f90 -c $(FFLAGS)
	$(F90) SwanCompUnstruc.f90 -c $(FFLAGS)
	$(F90) SwanDispParm.f90 -c $(FFLAGS)
	$(F90) SwanPropvelX.f90 -c $(FFLAGS)
	$(F90) SwanSweepSel.f90 -c $(FFLAGS)
	$(F90) SwanPropvelS.f90 -c $(FFLAGS)
	$(F90) SwanTranspAc.f90 -c $(FFLAGS)
	$(F90) SwanTranspX.f90 -c $(FFLAGS)
	$(F90) SwanDiffPar.f90 -c $(FFLAGS)
	$(F90) SwanGSECorr.f90 -c $(FFLAGS)
	$(F90) SwanInterpolatePoint.f90 -c $(FFLAGS)
	$(F90) SwanInterpolateAc.f90 -c $(FFLAGS)
	$(F90) SwanInterpolateOutput.f90 -c $(FFLAGS)
	$(F90) SwanConvAccur.f90 -c $(FFLAGS)
	$(F90) SwanConvStopc.f90 -c $(FFLAGS)
	$(F90) SwanThreadBounds.f90 -c $(FFLAGS)
	$(F90) SwanFindObstacles.f90 -c $(FFLAGS)
	$(F90) SwanCrossObstacle.f90 -c $(FFLAGS)
	$(F90) SwanComputeForce.f90 -c $(FFLAGS)
	$(F90) SwanIntgratSpc.f90 -c $(FFLAGS)
	$(F90) SwanBndStruc.f90 -c $(FFLAGS)
	$(F90) SwanReadfort18.f90 -c $(FFLAGS)
	$(F90) SwanSumOverNodes.f90 -c $(FFLAGS)
	$(F90) SwanMinOverNodes.f90 -c $(FFLAGS)
	$(F90) SwanMaxOverNodes.f90 -c $(FFLAGS)
	$(F90) ocpids.f -c $(F77FLAGS)
	$(F90) ocpcre.f -c $(F77FLAGS)
	$(F90) ocpmix.f -c $(F77FLAGS)

swmod1.o : swmod1.f
	$(F90) $< -c $(F77FLAGS)

swmod2.o : swmod2.f
	$(F90) $< -c $(F77FLAGS)

m_constants.o : m_constants.f90
	$(F90) $< -c $(FFLAGS)

m_fileio.o : m_fileio.f90
	$(F90) $< -c $(FFLAGS)

serv_xnl4v5.o : serv_xnl4v5.f90
	$(F90) $< -c $(FFLAGS)

SwanGriddata.o : SwanGriddata.f90
	$(F90) $< -c $(FFLAGS)

SwanGridobjects.o : SwanGridobjects.f90
	$(F90) $< -c $(FFLAGS)

SwanCompdata.o : SwanCompdata.f90
	$(F90) $< -c $(FFLAGS)

SwanCoupler.o : SwanCoupler.f90
	$(F90) $< -c $(FFLAGS)

SwanpyReef.o : SwanpyReef.f
	$(F90) $< -c $(F77FLAGS)

swanpre1.o : swanpre1.f
	$(F90) $< -c $(F77FLAGS)

swanpre2.o : swanpre2.f
	$(F90) $< -c $(F77FLAGS)

swancom1.o : swancom1.f
	$(F90) $< -c $(F77FLAGS)

swancom2.o : swancom2.f
	$(F90) $< -c $(F77FLAGS)

swancom3.o : swancom3.f
	$(F90) $< -c $(F77FLAGS)

swancom4.o : swancom4.f
	$(F90) $< -c $(F77FLAGS)

swancom5.o : swancom5.f
	$(F90) $< -c $(F77FLAGS)

swanout1.o : swanout1.f
	$(F90) $< -c $(F77FLAGS)

swanout2.o : swanout2.f
	$(F90) $< -c $(F77FLAGS)

swanser.o : swanser.f
	$(F90) $< -c $(F77FLAGS)

swanparll.o : swanparll.f
	$(F90) $< -c $(F77FLAGS)

SwanReadGrid.o : SwanReadGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanReadADCGrid.o : SwanReadADCGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanReadTriangleGrid.o : SwanReadTriangleGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanReadEasymeshGrid.o : SwanReadEasymeshGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanInitCompGrid.o : SwanInitCompGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanCheckGrid.o : SwanCheckGrid.f90
	$(F90) $< -c $(FFLAGS)

SwanCreateEdges.o : SwanCreateEdges.f90
	$(F90) $< -c $(FFLAGS)

SwanGridTopology.o : SwanGridTopology.f90
	$(F90) $< -c $(FFLAGS)

SwanGridVert.o : SwanGridVert.f90
	$(F90) $< -c $(FFLAGS)

SwanGridCell.o : SwanGridCell.f90
	$(F90) $< -c $(FFLAGS)

SwanGridFace.o : SwanGridFace.f90
	$(F90) $< -c $(FFLAGS)

SwanPrintGridInfo.o : SwanPrintGridInfo.f90
	$(F90) $< -c $(FFLAGS)

SwanFindPoint.o : SwanFindPoint.f90
	$(F90) $< -c $(FFLAGS)

SwanPointinMesh.o : SwanPointinMesh.f90
	$(F90) $< -c $(FFLAGS)

SwanBpntlist.o : SwanBpntlist.f90
	$(F90) $< -c $(FFLAGS)

SwanPrepComp.o : SwanPrepComp.f90
	$(F90) $< -c $(FFLAGS)

SwanVertlist.o : SwanVertlist.f90
	$(F90) $< -c $(FFLAGS)

SwanCompUnstruc.o : SwanCompUnstruc.f90
	$(F90) $< -c $(FFLAGS)

SwanDispParm.o : SwanDispParm.f90
	$(F90) $< -c $(FFLAGS)

SwanPropvelX.o : SwanPropvelX.f90
	$(F90) $< -c $(FFLAGS)

SwanSweepSel.o : SwanSweepSel.f90
	$(F90) $< -c $(FFLAGS)

SwanPropvelS.o : SwanPropvelS.f90
	$(F90) $< -c $(FFLAGS)

SwanTranspAc.o : SwanTranspAc.f90
	$(F90) $< -c $(FFLAGS)

SwanTranspX.o : SwanTranspX.f90
	$(F90) $< -c $(FFLAGS)

SwanDiffPar.o : SwanDiffPar.f90
	$(F90) $< -c $(FFLAGS)

SwanGSECorr.o : SwanGSECorr.f90
	$(F90) $< -c $(FFLAGS)

SwanInterpolatePoint.o : SwanInterpolatePoint.f90
	$(F90) $< -c $(FFLAGS)

SwanInterpolateAc.o : SwanInterpolateAc.f90
	$(F90) $< -c $(FFLAGS)

SwanInterpolateOutput.o : SwanInterpolateOutput.f90
	$(F90) $< -c $(FFLAGS)

SwanConvAccur.o : SwanConvAccur.f90
	$(F90) $< -c $(FFLAGS)

SwanConvStopc.o : SwanConvStopc.f90
	$(F90) $< -c $(FFLAGS)

SwanThreadBounds.o : SwanThreadBounds.f90
	$(F90) $< -c $(FFLAGS)

SwanFindObstacles.o : SwanFindObstacles.f90
	$(F90) $< -c $(FFLAGS)

SwanCrossObstacle.o : SwanCrossObstacle.f90
	$(F90) $< -c $(FFLAGS)

SwanComputeForce.o : SwanComputeForce.f90
	$(F90) $< -c $(FFLAGS)

SwanIntgratSpc.o : SwanIntgratSpc.f90
	$(F90) $< -c $(FFLAGS)

SwanBndStruc.o : SwanBndStruc.f90
	$(F90) $< -c $(FFLAGS)

SwanReadfort18.o : SwanReadfort18.f90
	$(F90) $< -c $(FFLAGS)

SwanSumOverNodes.o : SwanSumOverNodes.f90
	$(F90) $< -c $(FFLAGS)

SwanMinOverNodes.o : SwanMinOverNodes.f90
	$(F90) $< -c $(FFLAGS)

SwanMaxOverNodes.o : SwanMaxOverNodes.f90
	$(F90) $< -c $(FFLAGS)

ocpids.o : ocpids.f
	$(F90) $< -c $(F77FLAGS)

ocpcre.o : ocpcre.f
	$(F90) $< -c $(F77FLAGS)

ocpmix.o : ocpmix.f
	$(F90) $< -c $(F77FLAGS)

.f.o:
	$(F90) $< -c $(F77FLAGS)

.f90.o:
	$(F90) $< -c $(FFLAGS)
