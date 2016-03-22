# TODO: see how well this works when you need to mix ipyparallel and mpi-hidden code. hopefully %%px will work automatically! do a demo notebook to show this, perhaps

import numpy as np
import time

from pyBadlands import (diffLinear, elevationTIN, flowNetwork, FVmethod, partitionTIN,
                        raster2TIN, visualiseFlow, visualiseTIN)


class Model(object):
    """State object for the pyBadlands model."""

    def __init__(self):
        """Constructor."""
        # set default parameters

        # Planchon & Darboux filling thickness limit [m]
        self.fillmax = 1.
        # Precipitation rate [m/a]
        self.rainVal = 1.
        # Sea level position [m]
        self.seapos = -3450.
        # Limit flow network computation based
        # on relative sea level position [m]
        self.sealimit = 100.
        # Diffusion coefficients [m2/a]
        self.CDa = 0.5   # aerial
        self.CDm = 0.5   # marine
        # Stream power law coefficients
        self.SPLm = 0.5
        self.SPLn = 1.
        # Stream power law erodibility [m^(1-2*SPLm)/a]
        self.SPLero = 1.e-5

        # simulation time parameters
        self.tNow = 0.
        self.minDT = 1.

    def load_dem(self, filename):
        """Load a DEM input file."""
        # 1. Get DEM regular grid and create Badlands TIN.

        # single processor for now
        rank = 0
        size = 1

        # TODO: remove this and the outputDir stuff from raster2TIN
        import tempfile
        tempdir = tempfile.mkdtemp()

        recGrid = raster2TIN.raster2TIN(filename, outputDir=tempdir)
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

        if rank == 0:
            print " - partition TIN amongst processors ",time.clock() - walltime

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
        tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, allIDs,
                                                               totPts)

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
        local_elev[inIDs] = elevationTIN.getElevation(recGrid.regX, recGrid.regY, recGrid.regZ, FVmesh.node_coords[inIDs,:2])
        # comm.Allreduce(MPI.IN_PLACE, local_elev, op=MPI.MAX)  # not needed for single threaded

        elevation = elevationTIN.update_border_elevation(local_elev, FVmesh.neighbours, FVmesh.edge_length, recGrid.boundsPt, btype='slope')

        # set default of no rain
        self.rain = np.zeros(totPts, dtype=float)

        self.sealevel = self.seapos

        # Define variables
        self.cumdiff = np.zeros(totPts)
        self.hillslope = diffLinear()
        self.hillslope.CDaerial = self.CDa
        self.hillslope.CDmarine = self.CDm
        self.hillslope.dt_pstability(FVmesh.edge_length[inGIDs, :tMesh.maxNgbh])

        # FIXME: if the user changes the flow parameters (e.g. model.SPLero = something) they won't be updated correctly until load_dem is rerun
        self.flow = flowNetwork()
        self.flow.erodibility = self.SPLero
        self.flow.m = self.SPLm
        self.flow.n = self.SPLn
        self.flow.mindt = self.minDT

        if rank == 0:
            print " - interpolate elevation on grid ", time.clock() - walltime

        # tinData = np.column_stack((tinData, elevation))

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

    def compute_flow(self, tEnd):
        """
        Compute flows and update the model.

        tEnd: simulated time at which to stop computation
        """
        # single core
        size = 1
        rank = 0

        while self.tNow <= tEnd:
            Flow_time = time.clock()

            # 1. Perform pit filling
            walltime = time.clock()

            fillH = elevationTIN.pit_filling_PD(self.elevation, self.FVmesh.neighbours, self.recGrid.boundsPt, self.sealevel - self.sealimit, self.fillmax, 0.01)
            if rank == 0:
                print " - depression-less algorithm PD with stack", time.clock() - walltime

            # 2. Compute stream network
            walltime = time.clock()
            ngbhs = self.FVmesh.neighbours[self.allIDs, :]
            edges = self.FVmesh.vor_edges[self.allIDs, :]
            distances = self.FVmesh.edge_length[self.allIDs, :]
            self.flow.SFD_receivers(fillH, self.elevation, ngbhs, edges, distances, self.allIDs, self.sealevel - self.sealimit)
            if rank == 0:
                print " - compute receivers parallel ", time.clock() - walltime

            # Distribute evenly local minimas to processors
            walltime = time.clock()
            self.flow.localbase = np.array_split(self.flow.base, size)[rank]
            self.flow.ordered_node_array()
            if rank == 0:
                print " - compute stack order locally ", time.clock() - walltime

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
            if rank == 0:
                print " - send stack order globally ", time.clock() - walltime

            # 3. Compute discharge
            walltime = time.clock()
            self.flow.compute_flow(self.FVmesh.control_volumes, self.rain, True)
            if rank == 0:
                print " - compute discharge ", time.clock() - walltime

            # 4. Compute CFL condition
            walltime = time.clock()

            self.flow.dt_fstability(self.FVmesh.node_coords[:, :2], fillH, self.inGIDs)

            CFLtime = min(self.flow.CFL, self.hillslope.CFL)
            CFLtime = max(self.minDT, CFLtime)
            if rank == 0:
                print " - Get CFL time step ", time.clock() - walltime

            # 5. Compute sediment fluxes
            walltime = time.clock()
            # FIXME: not sure if diff_flux should be updating self.flow.diff_flux; not shown in notebook
            diff_flux = self.hillslope.sedflux(self.flow.diff_flux, self.sealevel, self.elevation, self.FVmesh.control_volumes)
            if rank == 0:
                print " - Get hillslope fluxes ", time.clock() - walltime

            walltime = time.clock()
            xyMin = [self.recGrid.regX.min(), self.recGrid.regY.min()]
            xyMax = [self.recGrid.regX.max(), self.recGrid.regY.max()]
            tstep, sedrate = self.flow.compute_sedflux(self.FVmesh.control_volumes, self.elevation, fillH,
                             self.FVmesh.node_coords[:, :2], xyMin, xyMax, diff_flux, CFLtime, self.sealevel, True)
            if rank == 0:
                print " - Get stream fluxes ", time.clock() - walltime

            # Update surface
            print("tstep = %s" % tstep)
            # FIXME: taking a guess here; I think that tstep is how long we simulated for the sediment flux, so I'm going to add that to the current time
            diff = sedrate * tstep
            self.elevation = self.elevation + diff
            self.cumdiff += diff
            self.tNow += tstep
            if rank == 0:
                print " - Flow computation ", time.clock() - Flow_time

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
        if not hasattr(self, 'tcells'):
            self.outPts, self.outCells = visualiseTIN.output_cellsIDs(self.allIDs, self.inIDs, visXlim, visYlim, self.FVmesh.node_coords[:, :2], self.tMesh.cells)
            self.tcells = np.zeros(size)
            self.tcells[rank] = len(self.outCells)
            # comm.Allreduce(MPI.IN_PLACE,tcells,op=MPI.MAX)
            self.tnodes = np.zeros(size)
            self.tnodes[rank] = len(self.allIDs)
            # comm.Allreduce(MPI.IN_PLACE,tnodes,op=MPI.MAX)

        # Done for every visualisation step
        flowIDs, polylines = visualiseFlow.output_Polylines(self.outPts, self.flow.receivers[self.outPts], visXlim, visYlim, self.FVmesh.node_coords[:, :2])
        fnodes = np.zeros(size)
        fnodes[rank] = len(flowIDs)
        # comm.Allreduce(MPI.IN_PLACE,fnodes,op=MPI.MAX)
        fline = np.zeros(size)
        fline[rank] = len(polylines[:, 0])
        # comm.Allreduce(MPI.IN_PLACE,fline,op=MPI.MAX)

        # Compute flow parameters
        self.flow.compute_parameters(self.FVmesh.node_coords[:, :2])

        th5file = 'tin.time'
        txmffile = 'tin.time'
        txdmffile = 'tin.series.xdmf'
        fh5file = 'flow.time'
        fxmffile = 'flow.time'
        fxdmffile = 'flow.series.xdmf'

        # Write HDF5 files
        visualiseTIN.write_hdf5(outDir, th5file, step, self.tMesh.node_coords[:, :2], self.elevation[self.allIDs],
                                self.flow.discharge[self.allIDs], self.cumdiff[self.allIDs], self.outCells, rank)
        visualiseFlow.write_hdf5(outDir, fh5file, step, self.FVmesh.node_coords[flowIDs, :2], self.elevation[flowIDs],
                                self.flow.discharge[flowIDs], self.flow.chi[flowIDs], self.flow.basinID[flowIDs], polylines, rank)

        # Combine HDF5 files and write time series
        if rank == 0:
            visualiseTIN.write_xmf(outDir, txmffile, txdmffile, step, self.tNow, self.tcells, self.tnodes, th5file, size)
            visualiseFlow.write_xmf(outDir, fxmffile, fxdmffile, step, self.tNow, fline, fnodes, fh5file, size)
            print " - Writing outputs ", time.clock() - out_time
