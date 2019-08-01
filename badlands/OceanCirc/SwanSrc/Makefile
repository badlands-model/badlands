UP=../..

TOP=$(shell pwd)/$(UP)

CONFFILE= $(TOP)/config/Makefile.inc
include $(CONFFILE)

TOPDIR=$(TOP)


F77FLAGS = -O3 -m64 -w -fno-second-underscore
#F77FLAGS= -O0 -g -fbounds-check -fbacktrace -finit-real=zero\
    -ffpe-trap=zero,overflow,invalid -fno-common\
    -Wtabs -Wuninitialized

SWAN_OBJS = swmod1.o swmod2.o m_constants.o m_fileio.o serv_xnl4v5.o \
	mod_xnl4v5.o SwanGriddata.o SwanGridobjects.o SwanCompdata.o SwanCoupler.o SwanSGFM.o \
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
	#make hconcat
	cp *.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		cp $$f $(MODDIR);       \
		fi; 		\
	done

hconcat: $(HCAT_OBJS)

swan: $(SWAN_OBJS)

.f.o:
	$(BADLANDS_F) $< -c $(F77FLAGS)

.f90.o:
	$(BADLANDS_F) $< -c $(FFLAGS)
