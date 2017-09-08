##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file defines the functions used to build Badlands meshes and surface grids.
"""

import time
import numpy as np
import mpi4py.MPI as mpi

from pyBadlands import (partitionTIN, FVmethod, elevationTIN, raster2TIN, oceanDyn,
                        eroMesh, strataMesh, isoFlex, stratiWedge, carbMesh, forceSim)

def construct_mesh(input, filename, verbose=False):
    """
    The following function is taking parsed values from the XML to:
        - build model grids & meshes,
        - initialise Finite Volume discretisation,
        - define the partitioning when parallelisation is enable.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD

    cumflex = None
    flex = None
    wave = None
    tinFlex = None
    strata = None
    mapero = None

    # Get DEM regular grid and create Badlands TIN.
    recGrid = raster2TIN.raster2TIN(filename, areaDelFactor=input.Afactor)

    fixIDs = recGrid.boundsPt + recGrid.edgesPt

    force = forceSim.forceSim(input.seafile, input.seapos, input.rainMap,
                input.rainTime, input.rainVal, input.orographic, input.orographiclin,
                input.rbgd, input.rmin, input.rmax, input.rzmax, input.windx,
                input.windy, input.tauc, input.tauf, input.nm,
                input.cw, input.hw, input.ortime, input.tectFile,
                input.tectTime, recGrid.regX, recGrid.regY, input.riverPos,
                input.riverTime, input.riverQws, input.riverRck, input.riverNb,
                input.rockNb, input.tDisplay)

    if input.disp3d:
        force.time3d = input.time3d
        if input.merge3d == 0. or input.merge3d > recGrid.resEdges:
            force.merge3d = input.Afactor * recGrid.resEdges * 0.5
        else:
            force.merge3d = input.merge3d

    # Partition the TIN
    walltime = time.clock()
    FVmesh = FVmethod.FVmethod(recGrid.tinMesh['vertices'], recGrid.tinMesh['triangles'],
                                recGrid.tinMesh['edges'])

    # Perform partitioning by equivalent domain splitting
    partitionIDs, RowProc, ColProc = partitionTIN.simple(recGrid.tinMesh['vertices'][:, 0],
                                                         recGrid.tinMesh['vertices'][:, 1])
    FVmesh.partIDs = partitionIDs

    # Get each partition global node ID
    inGIDs = np.where(partitionIDs == rank)[0]

    # Build Finite Volume discretisation
    # Define overlapping partitions
    lGIDs, localTIN = partitionTIN.overlap(recGrid.tinMesh['vertices'][:, 0], recGrid.tinMesh['vertices'][:, 1],
                                            RowProc, ColProc, 2*recGrid.resEdges, verbose)

    # Set parameters of the finite volume mesh
    tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'], localTIN['edges'])

    if rank == 0 and size > 1 and verbose:
        print " - partition TIN amongst processors and create local TINs", time.clock() - walltime

    # Define Finite Volume parameters
    walltime = time.clock()
    totPts = len(recGrid.tinMesh['vertices'][:, 0])
    FVmesh.neighbours = np.zeros((totPts, 20), dtype=np.int32, order='F')
    FVmesh.neighbours.fill(-2)
    FVmesh.edge_length = np.zeros((totPts, 20), dtype=np.float, order='F')
    FVmesh.vor_edges = np.zeros((totPts, 20), dtype=np.float, order='F')
    FVmesh.control_volumes = np.zeros(totPts, dtype=np.float)

    # Compute Finite Volume parameters
    tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, lGIDs, totPts,
                                                  recGrid.resEdges*input.Afactor, verbose)

    FVmesh.neighbours[tGIDs,:tMesh.maxNgbh] = tNgbh
    FVmesh.edge_length[tGIDs,:tMesh.maxNgbh] = tEdgs
    FVmesh.vor_edges[tGIDs,:tMesh.maxNgbh] = tVors
    FVmesh.control_volumes[tGIDs] = tVols

    if rank == 0 and verbose:
        print " - FV mesh ", time.clock() - walltime

    # Define TIN parameters
    if input.flexure:
        elevation, cumdiff, cumhill, cumflex, inIDs, parentIDs = _define_TINparams(totPts, input, FVmesh, recGrid, verbose)
    else:
        elevation, cumdiff, cumhill, inIDs, parentIDs = _define_TINparams(totPts, input, FVmesh, recGrid, verbose)

    # Build stratigraphic and erodibility meshes
    if input.laytime > 0 and input.erolays >= 0:
        strata, mapero = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, rank, verbose)
    elif input.laytime > 0:
        strata = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, rank, verbose)
    elif input.erolays >= 0:
        mapero = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, rank, verbose)

    # Set default to no rain
    force.update_force_TIN(FVmesh.node_coords[:,:2])

    # Flexural isostasy initialisation
    if input.flexure:
        flex, tinFlex, cumflex = _init_flexure(FVmesh, input, recGrid, force, elevation,
                                                cumdiff, cumflex, totPts, rank, verbose)

    # Wave grid initialisation
    if input.waveOn:
        wave = _init_wave(input, recGrid, force, rank, verbose)

    # Stratigraphic TIN initialisation
    if input.rockNb > 0:
        layNb = int((input.tEnd - input.tStart)/input.laytime)+2
        bPts = recGrid.boundsPt
        ePts = recGrid.edgesPt
        if input.restart:
            straTIN = stratiWedge.stratiWedge(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                            ePts, input.layersData, input.actlay, input.outDir, input.strath5file,
                            input.rockNb, recGrid.regX, recGrid.regY, elevation, input.rockCk, cumdiff,
                            input.rfolder, input.rstep)
        else:
            straTIN = stratiWedge.stratiWedge(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                                    ePts, input.layersData, input.actlay, input.outDir, input.strath5file,
                                    input.rockNb, recGrid.regX, recGrid.regY, elevation, input.rockCk)
    else:
        straTIN = None
        
    # Stratigraphic grid in case of carbonate and/or pelagic growth functions
    if input.carbonate or input.pelagic:
        layNb = int((input.tEnd - input.tStart)/input.tDisplay)+2
        bPts = recGrid.boundsPt
        ePts = recGrid.edgesPt
        if input.restart:
            carbTIN = carbMesh.carbMesh(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                            ePts, input.layersData, input.outDir, input.strath5file,
                            recGrid.regX, recGrid.regY, elevation,
                            input.rfolder, input.rstep)
        else:
            carbTIN = carbMesh.carbMesh(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                            ePts, input.layersData, input.outDir, input.strath5file,
                            recGrid.regX, recGrid.regY, elevation)
    else:
        carbTIN = None

    return recGrid, FVmesh, force, tMesh, lGIDs, fixIDs, \
        inIDs, parentIDs, inGIDs, totPts, elevation, cumdiff, \
        cumhill, cumflex, strata, mapero, tinFlex, flex, wave, \
        straTIN, carbTIN

def reconstruct_mesh(recGrid, input, verbose=False):
    """
    The following function is used after 3D displacements to:
        - rebuild model grids & meshes,
        - reinitialise Finite Volume discretisation,
        - redefine the partitioning when parallelisation is enable.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD

    walltime = time.clock()
    FVmesh = FVmethod.FVmethod(recGrid.tinMesh['vertices'], recGrid.tinMesh['triangles'],
                               recGrid.tinMesh['edges'])

    # Perform partitioning by equivalent domain splitting
    partitionIDs, RowProc, ColProc = partitionTIN.simple(recGrid.tinMesh['vertices'][:, 0],
                                                         recGrid.tinMesh['vertices'][:, 1])
    FVmesh.partIDs = partitionIDs

    # Get each partition global node ID
    inGIDs = np.where(partitionIDs == rank)[0]

    if rank == 0 and verbose:
        print " - partition TIN amongst processors ", time.clock() - walltime

    # Define overlapping partitions
    walltime = time.clock()
    lGIDs, localTIN = partitionTIN.overlap(recGrid.tinMesh['vertices'][:, 0],
                                           recGrid.tinMesh['vertices'][:, 1],
                                           RowProc, ColProc, 2*recGrid.resEdges,
                                           verbose)

    # Set parameters of the finite volume mesh
    tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'], localTIN['edges'])

    # Define Finite Volume parameters
    totPts = len(recGrid.tinMesh['vertices'][:, 0])
    FVmesh.neighbours = np.zeros((totPts, 20), dtype=np.int32, order='F')
    FVmesh.neighbours.fill(-2)
    FVmesh.edge_length = np.zeros((totPts, 20), dtype=np.float, order='F')
    FVmesh.vor_edges = np.zeros((totPts, 20), dtype=np.float, order='F')
    FVmesh.control_volumes = np.zeros(totPts, dtype=np.float)

    # Compute Finite Volume parameters
    tGIDs, tNgbh, tEdgs, tVors, tVols = tMesh.construct_FV(inGIDs, lGIDs, totPts,
                                            recGrid.resEdges*input.Afactor, verbose)

    FVmesh.neighbours[tGIDs,:tMesh.maxNgbh] = tNgbh
    FVmesh.edge_length[tGIDs,:tMesh.maxNgbh] = tEdgs
    FVmesh.vor_edges[tGIDs,:tMesh.maxNgbh] = tVors
    FVmesh.control_volumes[tGIDs] = tVols

    if rank == 0 and verbose:
        print " - reconstructed FV mesh ", time.clock() - walltime

    inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == rank)[0]
    inIDs += recGrid.boundsPt
    elevationTIN.assign_parameter_pit(FVmesh.neighbours, FVmesh.control_volumes, input.diffnb,
                                      input.diffprop, recGrid.boundsPt, input.fillmax)

    return FVmesh, tMesh, lGIDs, inIDs, inGIDs, totPts

def _define_TINparams(totPts, input, FVmesh, recGrid, verbose=False):
    """
    This function is defining the main values declared on the TIN.
    """

    # Initialise MPI communications
    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    walltime = time.clock()

    inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == rank)[0]
    inIDs += recGrid.boundsPt

    local_elev = np.zeros(totPts)
    local_elev.fill(-1.e6)

    # In case of a restart read values from HDF5 files
    if input.restart:
        local_cum = np.zeros(totPts)
        local_cum.fill(-1.e6)
        local_hill = np.zeros(totPts)
        local_hill.fill(-1.e6)
        if input.flexure:
            local_cumflex = np.zeros(totPts)
            local_cumflex.fill(-1.e6)
            local_elev[inIDs],local_cum[inIDs],local_hill[inIDs], local_cumflex[inIDs] = recGrid.load_hdf5_flex(input.rfolder,
                                                                input.rstep,FVmesh.node_coords[inIDs, :2])
        else:
            local_elev[inIDs],local_cum[inIDs],local_hill[inIDs] = recGrid.load_hdf5(input.rfolder,input.rstep,
                                                                FVmesh.node_coords[inIDs, :2])
        comm.Allreduce(mpi.IN_PLACE, local_elev, op=mpi.MAX)
        comm.Allreduce(mpi.IN_PLACE, local_cum, op=mpi.MAX)
        comm.Allreduce(mpi.IN_PLACE, local_hill, op=mpi.MAX)
        # Get cumulative erosion/deposition values
        cumdiff = local_cum
        cumdiff[:recGrid.boundsPt] = 0.
        cumhill = local_hill
        cumhill[:recGrid.boundsPt] = 0.
        if input.flexure:
            # Get cumulative flexural values
            comm.Allreduce(mpi.IN_PLACE, local_cumflex, op=mpi.MAX)
            cumflex = local_cumflex
            cumflex[:recGrid.boundsPt] = 0.

    # Otherwise interpolate elevation from DEM to TIN
    else:
        local_elev[inIDs] = elevationTIN.getElevation(recGrid.regX, recGrid.regY,
                                            recGrid.regZ, FVmesh.node_coords[inIDs, :2])
        comm.Allreduce(mpi.IN_PLACE, local_elev, op=mpi.MAX)
        # Initialise TIN parameters
        cumdiff = np.zeros(totPts)
        cumhill = np.zeros(totPts)
        if input.flexure:
            cumflex = np.zeros(totPts)

    # Assign boundary values
    elevation, parentIDs = elevationTIN.update_border_elevation(local_elev, FVmesh.neighbours,
                                FVmesh.edge_length, recGrid.boundsPt, btype=input.btype)

    # Define pit filling algorithm
    elevationTIN.assign_parameter_pit(FVmesh.neighbours, FVmesh.control_volumes, input.diffnb,
                                      input.diffprop, recGrid.boundsPt, input.fillmax)

    if rank == 0 and verbose:
        print " - define paramters on TIN grid ", time.clock() - walltime

    if input.flexure:
        return elevation, cumdiff, cumhill, cumflex, inIDs, parentIDs
    else:
        return elevation, cumdiff, cumhill, inIDs, parentIDs

def _build_strateroMesh(input, FVmesh, recGrid, cumdiff, rank, verbose=False):
    """
    This function is creating the stratigraphic mesh and the erodibility maps
    in cases where these functions are turned on.
    """

    # Build stratigraphic mesh
    if input.laytime > 0:

        walltime = time.clock()
        sdx = input.stratdx
        if sdx == 0:
            sdx = recGrid.rectX[1] - recGrid.rectX[0]
        bbX = [recGrid.rectX.min(),recGrid.rectX.max()]
        bbY = [recGrid.rectY.min(),recGrid.rectY.max()]
        layNb = int((input.tEnd - input.tStart)/input.laytime)+2
        strata = None
        if input.restart:
            strata = strataMesh.strataMesh(sdx, bbX, bbY, layNb, FVmesh.node_coords[:, :2],
                                input.outDir, input.sh5file, cumdiff, input.rfolder, input.rstep)
        else:
            strata = strataMesh.strataMesh(sdx, bbX, bbY, layNb, FVmesh.node_coords[:, :2],
                                input.outDir, input.sh5file)
        if rank == 0 and verbose:
            print " - create stratigraphic regions ", time.clock() - walltime

    # Define pre-existing erodibility maps
    if input.erolays >= 0:

        walltime = time.clock()
        bPts = recGrid.boundsPt
        if input.restart:
            mapero = eroMesh.eroMesh(input.erolays, input.eroMap, input.eroVal, input.SPLero,
                                    input.thickMap, input.thickVal, FVmesh.node_coords[:, :2],
                                    recGrid.regX, recGrid.regY, bPts, recGrid.edgesPt, input.outDir,
                                    rfolder=input.rfolder, rstep=input.rstep)
        else:
            mapero = eroMesh.eroMesh(input.erolays, input.eroMap, input.eroVal, input.SPLero,
                                     input.thickMap, input.thickVal, FVmesh.node_coords[:, :2],
                                     recGrid.regX, recGrid.regY, bPts, recGrid.edgesPt, input.outDir,
                                     rfolder=None, rstep=0)

        if rank == 0 and verbose:
            print " - create erodibility mesh ", time.clock() - walltime

    if input.laytime > 0 and input.erolays >= 0:
        return strata, mapero
    elif input.laytime > 0:
        return strata
    else:
        return mapero

def _init_flexure(FVmesh, input, recGrid, force, elevation, cumdiff,
                  cumflex, totPts, rank, verbose=False):
    """
    Initialise flexural isostasy.
    """

    # Initialise flexure parameters for gFlex library.
    walltime = time.clock()

    nx = input.fnx
    ny = input.fny
    if nx == 0:
        nx = recGrid.nx
    if ny == 0:
        ny = recGrid.ny
    if input.elasticH is None:
        elasticT = str(input.elasticGrid)
    else:
        elasticT = input.elasticH

    flex = isoFlex.isoFlex()
    flex.buildGrid(nx, ny, input.youngMod, input.dmantle, input.dsediment,
                elasticT, input.flexbounds, FVmesh.node_coords[:,:2])

    tinFlex = np.zeros(totPts, dtype=float)
    force.getSea(input.tStart)
    tinFlex = flex.get_flexure(elevation, cumdiff, force.sealevel,
                               recGrid.boundsPt, initFlex=True)
    tinFlex = force.disp_border(tinFlex, FVmesh.neighbours, FVmesh.edge_length, recGrid.boundsPt)
    cumflex += tinFlex
    if rank == 0 and verbose:
        print "   - Initialise flexural isostasy ", time.clock() - walltime

    return flex, tinFlex, cumflex

def _init_wave(input, recGrid, force, rank, verbose=False):
    """
    Initialise wave mesh.
    """

    # Initialise flexure parameters for gFlex library.
    walltime = time.clock()

    minX = recGrid.rectX.min()-input.resW
    maxX = recGrid.rectX.max()+input.resW
    minY = recGrid.rectY.min()-input.resW
    maxY = recGrid.rectY.max()+input.resW
    waveX = np.arange(minX,maxX+input.resW,input.resW)
    waveY = np.arange(minY,maxY+input.resW,input.resW)

    wave = oceanDyn.oceanDyn(input.resW,waveX,waveY)
    force.getSea(input.tStart)

    if rank == 0 and verbose:
        print "   - Initialise wave grid ", time.clock() - walltime

    return wave
