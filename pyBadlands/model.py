# TODO: see how well this works when you need to mix ipyparallel and mpi-hidden code. hopefully %%px will work automatically! do a demo notebook to show this, perhaps

import numpy as np
import time

from pyBadlands import (diffLinear, elevationTIN, flowNetwork, forceSim,
                        FVmethod, partitionTIN, raster2TIN, visualiseFlow,
                        visualiseTIN, xmlParser)


class Model(object):
    """State object for the pyBadlands model."""

    def __init__(self):
        """Constructor."""

        # simulation state
        self.tNow = 0.
        self.outputStep = 0
        self.simStarted = False

    def load_xml(self, filename):
        """Load an XML configuration file."""
        self.input = xmlParser.xmlParser(filename)
        self.load_dem(self.input.demfile)

    def load_dem(self, filename):
        """Load a DEM input file."""
        # 1. Get DEM regular grid and create Badlands TIN.

        # single processor for now
        rank = 0
        size = 1

        # TODO: remove this and the outputDir stuff from raster2TIN
        import tempfile
        tempdir = tempfile.mkdtemp()

        recGrid = raster2TIN.raster2TIN(filename, outputDir=tempdir, areaDelFactor=self.input.Afactor)

        self.fixIDs = recGrid.boundsPt + recGrid.edgesPt

        force = forceSim.forceSim(self.input.seafile, self.input.seapos,
                                  self.input.rainMap, self.input.rainTime,
                                  self.input.rainVal, self.input.tectFile,
                                  self.input.tectTime, recGrid.regX,
                                  recGrid.regY, self.input.tDisplay)

        if self.input.disp3d:
            force.time3d = self.input.time3d
            if self.input.merge3d == 0. or self.input.merge3d > recGrid.resEdges:
                force.merge3d = self.input.Afactor * recGrid.resEdges * 0.5

        tinData = np.column_stack((recGrid.tinMesh['vertices'][:, 0],
                                   recGrid.tinMesh['vertices'][:, 1]))

        # TODO: see if you can store just one of recGrid/tinData

        # 2. Partition the TIN
        walltime = time.clock()

        FVmesh = FVmethod.FVmethod(recGrid.tinMesh['vertices'],
                                   recGrid.tinMesh['triangles'],
                                   recGrid.tinMesh['edges'])

        # Perform partitioning by equivalent domain splitting
        rowProc = 1
        colProc = size
        partitionIDs, RowProc, ColProc = partitionTIN.simple(tinData[:, 0],
                                                             tinData[:, 1])
        FVmesh.partIDs = partitionIDs

        # Get each partition global node ID
        inGIDs = np.where(partitionIDs == rank)[0]

        if rank == 0 and size > 1:
            print " - partition TIN amongst processors ", time.clock() - walltime

        # 3. Build Finite Volume discretisation

        # Define overlapping partitions
        allIDs, localTIN = partitionTIN.overlap(tinData[:, 0], tinData[:, 1],
                                                rowProc, colProc,
                                                2 * recGrid.resEdges)

        # Set parameters of the finite volume mesh
        tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'],
                                  localTIN['edges'])

        walltime = time.clock()

        # Define Finite Volume parameters
        totPts = len(tinData[:, 0])
        FVmesh.neighbours = np.zeros((totPts, 20), dtype=np.int)
        FVmesh.neighbours.fill(-2)
        FVmesh.edge_length = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.vor_edges = np.zeros((totPts, 20), dtype=np.float)
        FVmesh.control_volumes = np.zeros(totPts, dtype=np.float)

        # Compute Finite Volume parameters
        tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, allIDs, totPts, recGrid.resEdges * self.input.Afactor)

        FVmesh.neighbours[tGIDs, :tMesh.maxNgbh] = tNgbh
        FVmesh.edge_length[tGIDs, :tMesh.maxNgbh] = tEdgs
        FVmesh.vor_edges[tGIDs, :tMesh.maxNgbh] = tVors
        FVmesh.control_volumes[tGIDs] = tVols

        if rank == 0:
            print " - FV mesh ", time.clock() - walltime

        # 4. Interpolate elevation
        walltime = time.clock()

        inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == rank)[0]
        inIDs += recGrid.boundsPt

        local_elev = np.zeros(totPts)
        local_elev.fill(-1.e6)
        local_elev[inIDs] = elevationTIN.getElevation(recGrid.regX, recGrid.regY, recGrid.regZ, FVmesh.node_coords[inIDs, :2])
        # comm.Allreduce(MPI.IN_PLACE, local_elev, op=MPI.MAX)  # not needed for single threaded

        elevation = elevationTIN.update_border_elevation(local_elev, FVmesh.neighbours, FVmesh.edge_length, recGrid.boundsPt, btype=self.input.btype)

        # set default of no rain
        self.rain = np.zeros(totPts, dtype=float)

        # Define variables
        self.cumdiff = np.zeros(totPts)
        self.hillslope = diffLinear()
        self.hillslope.CDaerial = self.input.CDa
        self.hillslope.CDmarine = self.input.CDm

        self.flow = flowNetwork()
        self.flow.erodibility = self.input.SPLero
        self.flow.m = self.input.SPLm
        self.flow.n = self.input.SPLn
        self.flow.mindt = self.input.minDT

        if rank == 0:
            print " - interpolate elevation on grid ", time.clock() - walltime

        # save state for subsequent calls
        # TODO: there is a lot of stuff here. Can we reduce it?
        self.tinData = tinData
        self.recGrid = recGrid
        self.elevation = elevation
        self.FVmesh = FVmesh
        self.allIDs = allIDs
        self.inIDs = inIDs
        self.inGIDs = inGIDs
        self.tMesh = tMesh
        self.force = force
        self.totPts = totPts

    def compute_flow(self, verbose=False):
        """
        Compute flows and update the model.

        tEnd: simulated time at which to stop computation
        """
        # single core
        size = 1
        rank = 0

        Flow_time = time.clock()

        # 1. Perform pit filling
        walltime = time.clock()

        # Update sea-level
        self.force.getSea(self.tNow)
        self.fillH = elevationTIN.pit_filling_PD(self.elevation, self.FVmesh.neighbours, self.recGrid.boundsPt, self.force.sealevel - self.input.sealimit, self.input.fillmax)
        if rank == 0 and verbose:
            print " -   depression-less algorithm PD with stack", time.clock() - walltime

        # 2. Compute stream network
        walltime = time.clock()
        ngbhs = self.FVmesh.neighbours[self.allIDs, :]
        edges = self.FVmesh.vor_edges[self.allIDs, :]
        distances = self.FVmesh.edge_length[self.allIDs, :]
        self.flow.SFD_receivers(self.fillH, self.elevation, ngbhs, edges, distances, self.allIDs, self.force.sealevel - self.input.sealimit)
        if rank == 0 and verbose:
            print " -   compute receivers parallel ", time.clock() - walltime

        # Distribute evenly local minimas to processors
        walltime = time.clock()
        self.flow.localbase = np.array_split(self.flow.base, size)[rank]
        self.flow.ordered_node_array()
        if rank == 0 and verbose:
            print " -   compute stack order locally ", time.clock() - walltime

        walltime = time.clock()
        # single threaded
        '''
        stackNbs = comm.allgather(len(flow.localstack))
        globalstack = np.zeros(sum(stackNbs),dtype=flow.localstack.dtype)
        comm.Allgatherv(sendbuf=[flow.localstack, MPI.INT],
                     recvbuf=[globalstack, (stackNbs, None), MPI.INT])
        flow.stack = globalstack
        '''
        self.flow.stack = self.flow.localstack
        if rank == 0 and verbose:
            print " -   send stack order globally ", time.clock() - walltime

        # 3. Compute discharge
        walltime = time.clock()
        self.flow.compute_flow(self.FVmesh.control_volumes, self.rain, True)
        if rank == 0 and verbose:
            print " -   compute discharge ", time.clock() - walltime

    def compute_flux(self, tEnd, verbose=False):
        # single core
        size = 1
        rank = 0

        self.exitTime = min([self.force.next_display, self.force.next_disp, self.force.next_rain])
        
        # 4. Compute CFL condition
        walltime = time.clock()
        self.hillslope.dt_pstability(self.FVmesh.edge_length[self.inGIDs, :self.tMesh.maxNgbh])

        self.flow.dt_fstability(self.FVmesh.node_coords[:, :2], self.fillH, self.inGIDs)

        CFLtime = min(self.flow.CFL, self.hillslope.CFL)
        CFLtime = max(self.input.minDT, CFLtime)
        if rank == 0 and verbose:
            print " -   Get CFL time step ", time.clock() - walltime

        # 5. Compute sediment fluxes
        # Initial cumulative elevation change
        walltime = time.clock()

        # FIXME: not sure if diff_flux should be updating self.flow.diff_flux; not shown in notebook
        diff_flux = self.hillslope.sedflux(self.flow.diff_flux, self.force.sealevel, self.elevation, self.FVmesh.control_volumes)
        if rank == 0 and verbose:
            print " -   Get hillslope fluxes ", time.clock() - walltime

        walltime = time.clock()
        xyMin = [self.recGrid.regX.min(), self.recGrid.regY.min()]
        xyMax = [self.recGrid.regX.max(), self.recGrid.regY.max()]
            
        tstep, sedrate = self.flow.compute_sedflux(self.FVmesh.control_volumes, self.elevation, self.fillH,
                         self.FVmesh.node_coords[:, :2], xyMin, xyMax, diff_flux, CFLtime, self.force.sealevel)
        if rank == 0 and verbose:
            print " -   Get stream fluxes ", time.clock() - walltime

        # Update surface
        timestep = min([tstep, tEnd, self.exitTime - self.tNow])
        diff = sedrate * timestep
        tID = np.where(self.fillH-self.elevation>0)[0]
        
        self.elevation += diff
        self.cumdiff += diff
        self.tNow += timestep
        
    def write_output(self, outDir, step):
        """
        Write HDF5 output.

        outDir: directory to write the output to
        step: TODO confirm what this is
        """

        # single core
        size = 1
        rank = 0

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
        tcells = np.zeros(size)
        tcells[rank] = len(outCells)
        # comm.Allreduce(MPI.IN_PLACE,tcells,op=MPI.MAX)
        tnodes = np.zeros(size)
        tnodes[rank] = len(self.allIDs)
        # comm.Allreduce(MPI.IN_PLACE,tnodes,op=MPI.MAX)

        # Done for every visualisation step
        flowIDs, polylines = visualiseFlow.output_Polylines(outPts, self.flow.receivers[outPts],
                        visXlim, visYlim, self.FVmesh.node_coords[:, :2])
        fnodes = np.zeros(size)
        fnodes[rank] = len(flowIDs)
        # comm.Allreduce(MPI.IN_PLACE,fnodes,op=MPI.MAX)
        fline = np.zeros(size)
        fline[rank] = len(polylines[:, 0])
        # comm.Allreduce(MPI.IN_PLACE,fline,op=MPI.MAX)

        # Compute flow parameters
        self.flow.compute_parameters(self.FVmesh.node_coords[:, :2])

        # Write HDF5 files
        visualiseTIN.write_hdf5(self.input.outDir, self.input.th5file, step, self.tMesh.node_coords[:,:2],
                                self.elevation[self.allIDs], self.flow.discharge[self.allIDs], self.cumdiff[self.allIDs],
                                outCells, rank) 
        visualiseFlow.write_hdf5(self.input.outDir, self.input.fh5file, step, self.FVmesh.node_coords[flowIDs, :2],
                                 self.elevation[flowIDs], self.flow.discharge[flowIDs], self.flow.chi[flowIDs],
                                 self.flow.basinID[flowIDs], polylines, rank)

        # Combine HDF5 files and write time series
        if rank == 0:
            visualiseTIN.write_xmf(self.input.outDir, self.input.txmffile, self.input.txdmffile,
                                   step, self.tNow, tcells, tnodes, self.input.th5file, size)
            visualiseFlow.write_xmf(self.input.outDir, self.input.fxmffile, self.input.fxdmffile,
                                    step, self.tNow, fline, fnodes, self.input.fh5file, size)
            print "   - Writing outputs (%0.02f seconds; tNow = %s)" % (time.clock() - out_time, self.tNow)

    def run_to_time(self, tEnd):
        """Run the simulation to a specified point in time (tEnd)."""
        
        if not self.simStarted:
            self.force.next_rain = self.force.T_rain[0, 0]
            self.force.next_disp = self.force.T_disp[0, 0]
            self.force.next_display = self.input.tStart
            self.exitTime = self.input.tEnd
            self.simStarted = True

        last_time = time.clock()
        last_output = time.clock()
        while self.tNow < tEnd:
            diff = time.clock() - last_time

            # at most, display output every 5 seconds
            if time.clock() - last_output >= 5.0:
                print 'tNow = %s (step took %0.02f seconds)' % (self.tNow, diff)
                last_output = time.clock()
            last_time = time.clock()

            # Load Rain Map
            if self.force.next_rain <= self.tNow and self.force.next_rain < self.input.tEnd:
                self.rain = np.zeros(self.totPts, dtype=float)
                self.rain[self.inIDs] = self.force.load_Rain_map(self.tNow, self.FVmesh.node_coords[self.inIDs, :2])
                # comm.Allreduce(MPI.IN_PLACE, self.rain, op=MPI.MAX)

            # Load Tectonic Grid
            if not self.input.disp3d:
                if self.force.next_disp <= self.tNow and self.force.next_disp < self.input.tEnd:
                    ldisp = np.zeros(self.totPts, dtype=float)
                    ldisp.fill(-1.e6)
                    ldisp[self.inIDs] = self.force.load_Tecto_map(self.tNow, self.FVmesh.node_coords[self.inIDs, :2])
                    # comm.Allreduce(MPI.IN_PLACE, ldisp, op=MPI.MAX)
                    self.disp = self.force.disp_border(ldisp, self.FVmesh.neighbours, self.FVmesh.edge_length, self.recGrid.boundsPt)
            else:
                if self.force.next_disp <= self.tNow and self.force.next_disp < self.input.tEnd:
                    self.force.load_Disp_map(self.tNow, self.FVmesh.node_coords[:, :2], self.inIDs)
                    self.force.disp_border(self.force.dispZ, self.FVmesh.neighbours, self.FVmesh.edge_length, self.recGrid.boundsPt)
                    self.recGrid.tinMesh, self.elevation, self.cumdiff = self.force.apply_XY_dispacements(self.recGrid.areaDel, self.fixIDs, self.elevation, self.cumdiff)

            # run the simulation for a bit
            self.compute_flow(verbose=False)
            
            if self.tNow >= self.force.next_display:
                # time to write output
                self.write_output(outDir=self.input.outDir, step=self.outputStep)
                self.force.next_display += self.input.tDisplay
                self.outputStep += 1
                last_output = time.clock()
                
            tStop = min([self.force.next_display, tEnd, self.force.next_disp, self.force.next_rain])
            self.compute_flux(tEnd=tStop, verbose=False)
        
        diff = time.clock() - last_time
        print 'tNow = %s (%0.02f seconds)' % (self.tNow, diff)
        self.write_output(outDir=self.input.outDir, step=self.outputStep)
        self.force.next_display += self.input.tDisplay
        self.outputStep += 1
