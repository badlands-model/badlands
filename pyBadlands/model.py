import time
import numpy as np
import mpi4py.MPI as mpi

from pyBadlands import (diffLinear, elevationTIN, flowNetwork, forceSim,
                        FVmethod, isoFlex, partitionTIN, raster2TIN,
                        visualiseFlow, visualiseTIN, xmlParser)

# profiling support
import cProfile
import os
import pstats
import StringIO


class Model(object):
    """State object for the pyBadlands model."""

    def __init__(self):
        """
        Constructor.
        """

        # simulation state
        self.tNow = 0.
        self.outputStep = 0
        self.applyDisp = False
        self.simStarted = False

        self._rank = mpi.COMM_WORLD.rank
        self._size = mpi.COMM_WORLD.size
        self._comm = mpi.COMM_WORLD

    def load_xml(self, filename, verbose=False):
        """
        Load an XML configuration file.
        """

        # only the first node should create a unique output dir
        self.input = xmlParser.xmlParser(filename, makeUniqueOutputDir=(self._rank == 0))
        self.tNow = self.input.tStart

        # sync the chosen output dir to all nodes
        self.input.outDir = self._comm.bcast(self.input.outDir, root=0)

        # seed the random number generator consistently on all nodes
        seed = None
        if self._rank == 0:
            # limit to max uint32
            seed = np.random.mtrand.RandomState().tomaxint() % 0xFFFFFFFF
        seed = self._comm.bcast(seed, root=0)
        np.random.seed(seed)

        self.load_dem(self.input.demfile, verbose)

    def load_dem(self, filename, verbose=False):
        """
        Load a DEM input file.
        """

        # 1. Get DEM regular grid and create Badlands TIN.

        recGrid = raster2TIN.raster2TIN(filename, areaDelFactor=self.input.Afactor)

        self.fixIDs = recGrid.boundsPt + recGrid.edgesPt

        force = forceSim.forceSim(self.input.seafile, self.input.seapos,
                                  self.input.rainMap, self.input.rainTime,
                                  self.input.rainVal, self.input.orographic,
                                  self.input.rbgd, self.input.rmin, self.input.rmax ,
                                  self.input.windx, self.input.windy, self.input.tauc,
                                  self.input.tauf, self.input.nm, self.input.cw,
                                  self.input.hw, self.input.ortime, self.input.tectFile,
                                  self.input.tectTime, recGrid.regX,
                                  recGrid.regY, self.input.tDisplay)

        if self.input.disp3d:
            force.time3d = self.input.time3d
            if self.input.merge3d == 0. or self.input.merge3d > recGrid.resEdges:
                force.merge3d = self.input.Afactor * recGrid.resEdges * 0.5
            else:
                force.merge3d = self.input.merge3d

        # 2. Partition the TIN
        walltime = time.clock()

        FVmesh = FVmethod.FVmethod(recGrid.tinMesh['vertices'],
                                   recGrid.tinMesh['triangles'],
                                   recGrid.tinMesh['edges'])

        # Perform partitioning by equivalent domain splitting
        partitionIDs, RowProc, ColProc = partitionTIN.simple(recGrid.tinMesh['vertices'][:, 0],
                                                             recGrid.tinMesh['vertices'][:, 1])
        FVmesh.partIDs = partitionIDs

        # Get each partition global node ID
        inGIDs = np.where(partitionIDs == self._rank)[0]

        if self._rank == 0 and self._size > 1 and verbose:
            print " - partition TIN amongst processors ", time.clock() - walltime

        # 3. Build Finite Volume discretisation

        # Define overlapping partitions
        allIDs, localTIN = partitionTIN.overlap(recGrid.tinMesh['vertices'][:, 0],
                                                recGrid.tinMesh['vertices'][:, 1],
                                                RowProc, ColProc,
                                                2 * recGrid.resEdges,
                                                verbose)

        # Set parameters of the finite volume mesh
        tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'],
                                  localTIN['edges'])

        walltime = time.clock()

        # Define Finite Volume parameters
        totPts = len(recGrid.tinMesh['vertices'][:, 0])
        FVmesh.neighbours = np.zeros((totPts, 20), dtype=np.int32)
        FVmesh.neighbours.fill(-2)
        FVmesh.edge_length = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.vor_edges = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.control_volumes = np.zeros(totPts, dtype=np.float)

        # Compute Finite Volume parameters
        tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, allIDs, totPts,
                                                  recGrid.resEdges * self.input.Afactor, verbose)

        FVmesh.neighbours[tGIDs, :tMesh.maxNgbh] = tNgbh
        FVmesh.edge_length[tGIDs, :tMesh.maxNgbh] = tEdgs
        FVmesh.vor_edges[tGIDs, :tMesh.maxNgbh] = tVors
        FVmesh.control_volumes[tGIDs] = tVols

        if self._rank == 0 and verbose:
            print " - FV mesh ", time.clock() - walltime

        # 4. Interpolate elevation
        walltime = time.clock()

        inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == self._rank)[0]
        inIDs += recGrid.boundsPt

        local_elev = np.zeros(totPts)
        local_elev.fill(-1.e6)

        if self.input.restart:
            local_cum = np.zeros(totPts)
            local_cum.fill(-1.e6)
            if self.input.flexure:
                local_cumflex = np.zeros(totPts)
                local_cumflex.fill(-1.e6)
                local_elev[inIDs],local_cum[inIDs],local_cumflex[inIDs] = recGrid.load_hdf5_flex(self.input.rfolder,
                                                                        self.input.rstep,FVmesh.node_coords[inIDs, :2])
            else:
                local_elev[inIDs],local_cum[inIDs] = recGrid.load_hdf5(self.input.rfolder,self.input.rstep,
                                                                        FVmesh.node_coords[inIDs, :2])
            self._comm.Allreduce(mpi.IN_PLACE, local_elev, op=mpi.MAX)
            self._comm.Allreduce(mpi.IN_PLACE, local_cum, op=mpi.MAX)
            self.cumdiff = local_cum
            self.cumdiff[:recGrid.boundsPt] = 0.
            if self.input.flexure:
                self._comm.Allreduce(mpi.IN_PLACE, local_cumflex, op=mpi.MAX)
                self.cumflex = local_cumflex
                self.cumflex[:recGrid.boundsPt] = 0.
        else:
            local_elev[inIDs] = elevationTIN.getElevation(recGrid.regX, recGrid.regY,
                                                          recGrid.regZ, FVmesh.node_coords[inIDs, :2])
            self._comm.Allreduce(mpi.IN_PLACE, local_elev, op=mpi.MAX)
            self.cumdiff = np.zeros(totPts)
            if self.input.flexure:
                self.cumflex = np.zeros(totPts)

        elevation = elevationTIN.update_border_elevation(local_elev, FVmesh.neighbours,
                                                         FVmesh.edge_length, recGrid.boundsPt,
                                                         btype=self.input.btype)
        # Set default to no rain
        force.update_force_TIN(FVmesh.node_coords[:,:2])
        self.rain = np.zeros(totPts, dtype=float)

        # Define variables
        self.hillslope = diffLinear()
        self.hillslope.CDaerial = self.input.CDa
        self.hillslope.CDmarine = self.input.CDm

        self.flow = flowNetwork()
        self.flow.erodibility = self.input.SPLero
        self.flow.m = self.input.SPLm
        self.flow.n = self.input.SPLn
        self.flow.mindt = self.input.minDT
        self.flow.bedrock = self.input.bedrock
        self.flow.alluvial = self.input.alluvial
        self.flow.esmooth = self.input.esmooth
        self.flow.dsmooth = self.input.dsmooth
        self.flow.xycoords = FVmesh.node_coords[:,:2]
        self.flow.spl = self.input.spl
        self.flow.capacity = self.input.capacity
        self.flow.filter = self.input.filter
        self.flow.depo = self.input.depo

        if self._rank == 0 and verbose:
            print " - interpolate elevation on grid ", time.clock() - walltime

        # Flexural isostasy initialisation
        if self.input.flexure:
            walltime = time.clock()
            nx = self.input.fnx
            ny = self.input.fny
            if nx == 0:
                nx = recGrid.nx
            if ny == 0:
                ny = recGrid.ny
            if self.input.elasticH is None:
                elasticT = str(self.input.elasticGrid)
            else:
                elasticT = self.input.elasticH

            self.flex = isoFlex.isoFlex()
            self.flex.buildGrid(nx, ny, self.input.youngMod, self.input.dmantle, self.input.dsediment,
                elasticT, self.input.flexbounds, FVmesh.node_coords[:,:2])

            self.tinFlex = np.zeros(totPts, dtype=float)
            force.getSea(self.tNow)
            self.tinFlex = self.flex.get_flexure(elevation, self.cumdiff, force.sealevel,
                                recGrid.boundsPt, initFlex=True)
            self.tinFlex = force.disp_border(self.tinFlex, FVmesh.neighbours, FVmesh.edge_length, recGrid.boundsPt)
            self.cumflex += self.tinFlex
            if self._rank == 0:
                print "   - Compute flexural isostasy ", time.clock() - walltime

        # Save state for subsequent calls
        # TODO: there is a lot of stuff here. Can we reduce it?
        self.recGrid = recGrid
        self.elevation = elevation
        self.FVmesh = FVmesh
        self.allIDs = allIDs
        self.inIDs = inIDs
        self.inGIDs = inGIDs
        self.tMesh = tMesh
        self.force = force
        self.totPts = totPts

    def rebuildMesh(self, verbose=False):
        """
        Build TIN after 3D displacements.
        """

        # Build the Finite Volume representation
        self.fixIDs = self.recGrid.boundsPt + self.recGrid.edgesPt
        walltime = time.clock()

        FVmesh = FVmethod.FVmethod(self.recGrid.tinMesh['vertices'],
                                   self.recGrid.tinMesh['triangles'],
                                   self.recGrid.tinMesh['edges'])

        # Perform partitioning by equivalent domain splitting
        partitionIDs, RowProc, ColProc = partitionTIN.simple(self.recGrid.tinMesh['vertices'][:, 0],
                                                             self.recGrid.tinMesh['vertices'][:, 1])
        FVmesh.partIDs = partitionIDs

        # Get each partition global node ID
        inGIDs = np.where(partitionIDs == self._rank)[0]

        if self._rank == 0 and verbose:
            print " - partition TIN amongst processors ", time.clock() - walltime

        # Define overlapping partitions
        allIDs, localTIN = partitionTIN.overlap(self.recGrid.tinMesh['vertices'][:, 0],
                                                self.recGrid.tinMesh['vertices'][:, 1],
                                                RowProc, ColProc,
                                                2 * self.recGrid.resEdges,
                                                verbose)

        # Set parameters of the finite volume mesh
        tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'],
                                  localTIN['edges'])

        walltime = time.clock()

        # Define Finite Volume parameters
        totPts = len(self.recGrid.tinMesh['vertices'][:, 0])
        FVmesh.neighbours = np.zeros((totPts, 20), dtype=np.int32)
        FVmesh.neighbours.fill(-2)
        FVmesh.edge_length = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.vor_edges = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.control_volumes = np.zeros(totPts, dtype=np.float)

        # Compute Finite Volume parameters
        tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, allIDs, totPts,
                                                  self.recGrid.resEdges * self.input.Afactor,
                                                  verbose)

        FVmesh.neighbours[tGIDs, :tMesh.maxNgbh] = tNgbh
        FVmesh.edge_length[tGIDs, :tMesh.maxNgbh] = tEdgs
        FVmesh.vor_edges[tGIDs, :tMesh.maxNgbh] = tVors
        FVmesh.control_volumes[tGIDs] = tVols

        if self._rank == 0 and verbose:
            print " - FV mesh ", time.clock() - walltime

        inIDs = np.where(FVmesh.partIDs[self.recGrid.boundsPt:] == self._rank)[0]
        inIDs += self.recGrid.boundsPt

        # Set default with no rain
        self.force.update_force_TIN(FVmesh.node_coords[:,:2])
        self.rain = np.zeros(totPts, dtype=float)
        self.rain[inIDs] = self.force.get_Rain(self.tNow,self.elevation, inIDs)

        # Update flexural isostasy
        if self.input.flexure:
            self.tinFlex = np.zeros(totPts, dtype=float)
            self.flex.update_flexure_parameters(FVmesh.node_coords[:,:2])

        # Save state for subsequent calls
        # TODO: there is a lot of stuff here. Can we reduce it?
        self.flow.xycoords = FVmesh.node_coords[:, :2]
        self.FVmesh = FVmesh
        self.allIDs = allIDs
        self.inIDs = inIDs
        self.inGIDs = inGIDs
        self.tMesh = tMesh
        self.totPts = totPts


    def compute_flow(self, verbose=False):
        """
        Compute flows and update the model.
        """

        Flow_time = time.clock()

        # 1. Perform pit filling
        walltime = time.clock()

        # Update sea-level
        self.force.getSea(self.tNow)

        self.fillH = None

        if self.input.depo == 0 or self.input.capacity or self.input.filter:
            self.flow.maxdep = 0.
            self.flow.maxh = 0.
            if self.flow.esmooth is None and self.input.filter:
                self.flow.esmooth = 0
            if self.flow.dsmooth is None and self.input.filter:
                self.flow.dsmooth = 0
        else:
            self.fillH = elevationTIN.pit_filling_PD(self.elevation, self.FVmesh.neighbours,
                                                 self.recGrid.boundsPt,
                                                 self.force.sealevel - self.input.sealimit,
                                                 self.input.fillmax)

        if self._rank == 0 and verbose and self.input.spl and not self.input.filter:
            print " -   depression-less algorithm PD with stack", time.clock() - walltime

        # 2. Compute stream network
        walltime = time.clock()
        ngbhs = self.FVmesh.neighbours[self.allIDs, :]
        edges = self.FVmesh.vor_edges[self.allIDs, :]
        distances = self.FVmesh.edge_length[self.allIDs, :]
        self.flow.SFD_receivers(self.fillH, self.elevation, ngbhs, edges, distances,
                                self.allIDs, self.force.sealevel - self.input.sealimit)

        if self._rank == 0 and verbose:
            print " -   compute receivers parallel ", time.clock() - walltime

        # Distribute evenly local minimas to processors
        walltime = time.clock()
        self.flow.localbase = np.array_split(self.flow.base, self._size)[self._rank]
        self.flow.ordered_node_array()
        if self._rank == 0 and verbose:
            print " -   compute stack order locally ", time.clock() - walltime

        walltime = time.clock()
        stackNbs = self._comm.allgather(len(self.flow.localstack))
        globalstack = np.zeros(sum(stackNbs),dtype=self.flow.localstack.dtype)
        self._comm.Allgatherv(sendbuf=[self.flow.localstack, mpi.INT],
                              recvbuf=[globalstack, (stackNbs, None), mpi.INT])
        self.flow.stack = globalstack

        if self._rank == 0 and verbose:
            print " -   send stack order globally ", time.clock() - walltime

        # 3. Compute discharge
        walltime = time.clock()
        self.flow.compute_flow(self.FVmesh.control_volumes, self.rain)
        if self._rank == 0 and verbose:
            print " -   compute discharge ", time.clock() - walltime

    def compute_flux(self, tEnd, verbose=False):
        """
        Compute sediment fluxes.
        """

        # 1. Compute CFL condition
        walltime = time.clock()
        if self.input.Hillslope:
            self.hillslope.dt_pstability(self.FVmesh.edge_length[self.inGIDs, :self.tMesh.maxNgbh])
        else:
            self.hillslope.CFL = tEnd - self.tNow

        if self.input.depo == 1:
            if self.input.filter:
                self.flow.CFL = self.input.maxDT
            elif self.input.spl:
                self.flow.dt_fstability(self.fillH, self.inGIDs)
            else:
                self.flow.dt_fstability(self.elevation, self.inGIDs)
        else:
            if self.input.filter:
                self.flow.CFL = self.input.maxDT
            else:
                self.flow.dt_fstability(self.elevation, self.inGIDs)

        CFLtime = min(self.flow.CFL, self.hillslope.CFL)
        CFLtime = max(self.input.minDT, CFLtime)
        if self._rank == 0 and verbose:
            print " -   Get CFL time step ", time.clock() - walltime

        # 2. Compute sediment fluxes
        # Initial cumulative elevation change
        walltime = time.clock()
        diff_flux = self.hillslope.sedflux(self.flow.diff_flux, self.force.sealevel, self.elevation,
                                           self.FVmesh.control_volumes)
        if self._rank == 0 and verbose:
            print " -   Get hillslope fluxes ", time.clock() - walltime

        walltime = time.clock()
        xyMin = [self.recGrid.regX.min(), self.recGrid.regY.min()]
        xyMax = [self.recGrid.regX.max(), self.recGrid.regY.max()]
        tstep, sedrate = self.flow.compute_sedflux(self.FVmesh.control_volumes, self.elevation, self.fillH,
                                            xyMin, xyMax,diff_flux, CFLtime, self.force.sealevel, self.cumdiff)
        if self._rank == 0 and verbose:
            print " -   Get stream fluxes ", time.clock() - walltime

        # Update surface parameters and time
        timestep = min(tstep, tEnd - self.tNow)
        diff = sedrate * timestep

        if self.input.filter:
            smthdiff = self.flow.gaussian_filter(diff)
            smthdiff[:self.recGrid.boundsPt] = 0.
            self.elevation += smthdiff
            self.cumdiff += smthdiff
        else:
            self.elevation += diff
            self.cumdiff += diff

        if self.applyDisp:
            self.elevation += self.disp * timestep

        self.tNow += timestep

        if self._rank == 0 and verbose:
            print " - Flow computation ", time.clock() - Flow_time

    def write_output(self, outDir, step):
        """
        Write HDF5 output.
        """

        out_time = time.clock()
        visXlim = np.zeros(2)
        visYlim = np.zeros(2)
        visXlim[0] = self.recGrid.rectX.min()
        visXlim[1] = self.recGrid.rectX.max()
        visYlim[0] = self.recGrid.rectY.min()
        visYlim[1] = self.recGrid.rectY.max()

        # Done when TIN has been built/rebuilt
        outPts, outCells = visualiseTIN.output_cellsIDs(self.allIDs, self.inIDs, visXlim, visYlim,
                                        self.FVmesh.node_coords[:, :2], self.tMesh.cells)
        tcells = np.zeros(self._size)
        tcells[self._rank] = len(outCells)
        self._comm.Allreduce(mpi.IN_PLACE, tcells, op=mpi.MAX)
        tnodes = np.zeros(self._size)
        tnodes[self._rank] = len(self.allIDs)
        self._comm.Allreduce(mpi.IN_PLACE, tnodes, op=mpi.MAX)

        # Done for every visualisation step
        flowIDs, polylines = visualiseFlow.output_Polylines(outPts, self.flow.receivers[outPts],
                        visXlim, visYlim, self.FVmesh.node_coords[:, :2])
        fnodes = np.zeros(self._size)
        fnodes[self._rank] = len(flowIDs)
        self._comm.Allreduce(mpi.IN_PLACE, fnodes, op=mpi.MAX)
        fline = np.zeros(self._size)
        fline[self._rank] = len(polylines[:, 0])
        self._comm.Allreduce(mpi.IN_PLACE, fline, op=mpi.MAX)

        # Compute flow parameters
        self.flow.compute_parameters()

        # Write HDF5 files
        if self.input.flexure:
            visualiseTIN.write_hdf5_flexure(self.input.outDir, self.input.th5file, step, self.tMesh.node_coords[:,:2],
                                    self.elevation[self.allIDs], self.rain[self.allIDs], self.flow.discharge[self.allIDs],
                                    self.cumdiff[self.allIDs], self.cumflex[self.allIDs], outCells, self._rank, self.input.oroRain)
        else:
            visualiseTIN.write_hdf5(self.input.outDir, self.input.th5file, step, self.tMesh.node_coords[:,:2],
                                    self.elevation[self.allIDs], self.rain[self.allIDs], self.flow.discharge[self.allIDs],
                                    self.cumdiff[self.allIDs], outCells, self._rank, self.input.oroRain)
        visualiseFlow.write_hdf5(self.input.outDir, self.input.fh5file, step, self.FVmesh.node_coords[flowIDs, :2],
                                 self.elevation[flowIDs], self.flow.discharge[flowIDs], self.flow.chi[flowIDs],
                                 self.flow.basinID[flowIDs], polylines, self._rank)

        # Combine HDF5 files and write time series
        if self._rank == 0:
            visualiseTIN.write_xmf(self.input.outDir, self.input.txmffile, self.input.txdmffile, step, self.tNow,
                                   tcells, tnodes, self.input.th5file, self.force.sealevel, self._size,
                                   self.input.flexure, self.input.oroRain)
            visualiseFlow.write_xmf(self.input.outDir, self.input.fxmffile, self.input.fxdmffile,
                                    step, self.tNow, fline, fnodes, self.input.fh5file, self._size)
            print "   - Writing outputs (%0.02f seconds; tNow = %s)" % (time.clock() - out_time, self.tNow)

    def run_to_time(self, tEnd, profile=False):
        """
        Run the simulation to a specified point in time (tEnd).

        If profile is True, dump cProfile output to /tmp.
        """

        if profile:
            pid = os.getpid()
            pr = cProfile.Profile()
            pr.enable()

        if not self.simStarted:
            # Anything in here will be executed once at the start of time
            self.force.next_rain = self.force.T_rain[0, 0]
            self.force.next_disp = self.force.T_disp[0, 0]
            self.force.next_display = self.input.tStart
            self.exitTime = self.input.tEnd
            if self.input.flexure:
                self.force.next_flexure = self.input.tStart + self.input.ftime
            else:
                self.force.next_flexure = self.exitTime + self.input.tDisplay
            self.simStarted = True

        last_time = time.clock()
        last_output = time.clock()
        while self.tNow < tEnd:
            diff = time.clock() - last_time

            # At most, display output every 5 seconds
            if time.clock() - last_output >= 5.0:
                print 'tNow = %s (step took %0.02f seconds)' % (self.tNow, diff)
                last_output = time.clock()
            last_time = time.clock()

            # Load Rain regime
            if self.force.next_rain <= self.tNow and self.force.next_rain < self.input.tEnd:
                self.rain = np.zeros(self.totPts, dtype=float)
                self.rain[self.inIDs] = self.force.get_Rain(self.tNow, self.elevation, self.inIDs)
                self._comm.Allreduce(mpi.IN_PLACE, self.rain, op=mpi.MAX)

            # Load Tectonic Grid
            if not self.input.disp3d:
                if self.force.next_disp <= self.tNow and self.force.next_disp < self.input.tEnd:
                    ldisp = np.zeros(self.totPts, dtype=float)
                    ldisp.fill(-1.e6)
                    ldisp[self.inIDs] = self.force.load_Tecto_map(self.tNow,self.inIDs)
                    self._comm.Allreduce(mpi.IN_PLACE, ldisp, op=mpi.MAX)
                    self.disp = self.force.disp_border(ldisp, self.FVmesh.neighbours,
                                                       self.FVmesh.edge_length, self.recGrid.boundsPt)
                    self.applyDisp = True
            else:
                if self.force.next_disp <= self.tNow and self.force.next_disp < self.input.tEnd:
                    updateMesh = self.force.load_Disp_map(self.tNow, self.FVmesh.node_coords[:, :2], self.inIDs)
                    if updateMesh:
                        self.force.dispZ = self.force.disp_border(self.force.dispZ, self.FVmesh.neighbours,
                                           self.FVmesh.edge_length, self.recGrid.boundsPt)
                        if self.input.flexure:
                            self.recGrid.tinMesh, self.elevation, self.cumdiff, self.cumflex = self.force.apply_XY_dispacements_flexure(
                                self.recGrid.areaDel, self.fixIDs, self.elevation, self.cumdiff, self.cumflex)
                        else:
                            self.recGrid.tinMesh, self.elevation, self.cumdiff = self.force.apply_XY_dispacements(
                                self.recGrid.areaDel, self.fixIDs, self.elevation, self.cumdiff)
                        self.rebuildMesh()

            # Run the simulation for a bit
            self.compute_flow(verbose=False)

            # Isostatic flexure
            if self.tNow >= self.force.next_flexure:
                flextime = time.clock()
                self.force.getSea(self.tNow)
                self.tinFlex = self.flex.get_flexure(self.elevation, self.cumdiff,
                            self.force.sealevel,self.recGrid.boundsPt,initFlex=False)
                # Get border values
                self.tinFlex = self.force.disp_border(self.tinFlex, self.FVmesh.neighbours,
                                                      self.FVmesh.edge_length, self.recGrid.boundsPt)
                # Update flexural parameters
                self.elevation += self.tinFlex
                self.cumflex += self.tinFlex
                # Update next flexure time
                self.force.next_flexure += self.input.ftime
                if self._rank == 0:
                    print "   - Compute flexural isostasy ", time.clock() - flextime

            if self.tNow >= self.force.next_display:
                # time to write output
                self.write_output(outDir=self.input.outDir, step=self.outputStep)
                # Update next display time
                self.force.next_display += self.input.tDisplay
                self.outputStep += 1
                last_output = time.clock()

            tStop = min([self.force.next_display, self.force.next_flexure, tEnd, self.force.next_disp, self.force.next_rain])
            self.compute_flux(tEnd=tStop, verbose=False)

        diff = time.clock() - last_time
        print 'tNow = %s (%0.02f seconds)' % (self.tNow, diff)
        # Isostatic flexure
        if self.input.flexure:
            flextime = time.clock()
            self.force.getSea(self.tNow)
            self.tinFlex = self.flex.get_flexure(self.elevation, self.cumdiff,
                        self.force.sealevel,self.recGrid.boundsPt,initFlex=False)
            # Get border values
            self.tinFlex = self.force.disp_border(self.tinFlex, self.FVmesh.neighbours,
                                                  self.FVmesh.edge_length, self.recGrid.boundsPt)
            # Update flexural parameters
            self.elevation += self.tinFlex
            self.cumflex += self.tinFlex
            # Update next flexure time
            self.force.next_flexure += self.input.ftime
            if self._rank == 0:
                print "   - Compute flexural isostasy ", time.clock() - flextime
        # Output
        self.write_output(outDir=self.input.outDir, step=self.outputStep)
        self.force.next_display += self.input.tDisplay
        self.outputStep += 1

        if profile:
            pr.disable()
            s = StringIO.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.dump_stats('/tmp/profile-%d' % pid)

    def ncpus(self):
        """Return the number of CPUs used to generate the results."""
        return 1
