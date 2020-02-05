##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file defines the functions used to build **badlands** meshes and surface grids.
"""

import os
import time
import numpy as np
from scipy.interpolate import griddata

if 'READTHEDOCS' not in os.environ:
    from badlands import (partitionTIN, FVmethod, elevationTIN, raster2TIN, waveSed,
                            eroMesh, strataMesh, isoFlex, stratiWedge, carbMesh, forceSim)


def construct_mesh(input, filename, verbose=False):
    """
    The following function is taking parsed values from the XML to:

    * build model grids & meshes,
    * initialise Finite Volume discretisation,
    * define the partitioning when parallelisation is enable.

    Args:
        input: class containing XML input file parameters.
        filename: (str) this is a string containing the path to the regular grid file.
        verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).

    Returns
    -------
    recGrid
        class describing the regular grid characteristics.
    FVmesh
        class describing the finite volume mesh.
    force
        class describing the forcing parameters.
    tMesh
        class describing the TIN mesh.
    lGIDs
        numpy 1D array containing the node indices.
    fixIDs
        numpy 1D array containing the fixed node indices.
    inGIDs
        numpy 1D array containing the node indices inside the mesh.
    totPts
        total number of points in the mesh.
    elevation
        numpy array containing the elevations for the domain.
    cumdiff
        cumulative total erosion/deposition changes
    cumhill
        cumulative hillslope erosion/deposition changes
    cumfail
        cumulative failure induced erosion/deposition changes
    cumflex
        cumlative changes induced by flexural isostasy
    strata
        stratigraphic class parameters
    mapero
        underlying erodibility map characteristics
    tinFlex
        class describing the flexural TIN mesh.
    flex
        class describing the flexural isostasy functions.
    wave
        class describing the wave functions.
    straTIN
        class describing the stratigraphic TIN mesh.
    carbTIN
        class describing the carbonate TIN mesh.
    """

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
                input.rockNb, input.tDisplay, input.carbValSp1, input.carbValSp2,
                input.carbTime)

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
    inGIDs = np.where(partitionIDs == 0)[0]

    # Build Finite Volume discretisation
    # Define overlapping partitions
    lGIDs, localTIN = partitionTIN.overlap(recGrid.tinMesh['vertices'][:, 0], recGrid.tinMesh['vertices'][:, 1],
                                            RowProc, ColProc, 2*recGrid.resEdges, verbose)

    # Set parameters of the finite volume mesh
    tMesh = FVmethod.FVmethod(localTIN['vertices'], localTIN['triangles'], localTIN['edges'])

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

    if verbose:
        print(" - FV mesh ", time.clock() - walltime)

    # Define TIN parameters
    if input.flexure:
        elevation, cumdiff, cumhill, cumfail, cumflex, inIDs, parentIDs = _define_TINparams(totPts, input, FVmesh, recGrid, verbose)
    else:
        elevation, cumdiff, cumhill, cumfail, inIDs, parentIDs = _define_TINparams(totPts, input, FVmesh, recGrid, verbose)

    # Build stratigraphic and erodibility meshes
    if ((input.laytime and input.laytime > 0) and
       (input.erolays and input.erolays >= 0)):
        strata, mapero = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, verbose)
    elif (input.laytime and input.laytime > 0):
        strata = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, verbose)
    elif (input.erolays and input.erolays >= 0):
        mapero = _build_strateroMesh(input, FVmesh, recGrid, cumdiff, verbose)

    # Set default to no rain
    force.update_force_TIN(FVmesh.node_coords[:,:2])

    # Flexural isostasy initialisation
    if input.flexure:
        flex, tinFlex, cumflex = _init_flexure(FVmesh, input, recGrid, force, elevation,
                                                cumdiff, cumflex, totPts, verbose)

    # Wavesed grid initialisation
    if input.waveSed:
        ref_elev = get_reference_elevation(input,recGrid,elevation)
        wave = _init_wavesed(input,ref_elev, recGrid, force, verbose)
        wave.build_tree(FVmesh.node_coords[:,:2])

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
    if input.carbonate:
        layNb = int((input.tEnd - input.tStart)/input.tDisplay)+2
        bPts = recGrid.boundsPt
        ePts = recGrid.edgesPt
        if input.carbonate2:
            nbSed = 3
        else:
            nbSed = 2

        if input.restart:
            carbTIN = carbMesh.carbMesh(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                            ePts, input.layersData, input.outDir, input.strath5file, input.baseMap, nbSed,
                            recGrid.regX, recGrid.regY, elevation, input.rfolder, input.rstep)
        else:
            carbTIN = carbMesh.carbMesh(layNb, input.initlayers, FVmesh.node_coords[:, :2], bPts,
                            ePts, input.layersData, input.outDir, input.strath5file, input.baseMap, nbSed,
                            recGrid.regX, recGrid.regY, elevation)
    else:
        carbTIN = None

    return recGrid, FVmesh, force, tMesh, lGIDs, fixIDs, \
        inIDs, parentIDs, inGIDs, totPts, elevation, cumdiff, \
        cumhill, cumfail, cumflex, strata, mapero, tinFlex, flex, wave, \
        straTIN, carbTIN

def reconstruct_mesh(recGrid, input, verbose=False):
    """
    The following function is used after 3D displacements to:

    * rebuild model grids & meshes,
    * reinitialise Finite Volume discretisation,
    * redefine the partitioning when parallelisation is enable.

    Args:
        recGrid: class describing the regular grid characteristics.
        input: class containing XML input file parameters.
        verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).

    Returns
    -------
    FVmesh
        class describing the finite volume mesh.
    tMesh
        class describing the TIN mesh.
    lGIDs
        numpy 1D array containing the node indices.
    inIDs
        numpy 1D array containing the local node indices inside the mesh.
    inGIDs
        numpy 1D array containing the node indices inside the mesh.
    totPts
        total number of points in the mesh.
    """

    walltime = time.clock()
    FVmesh = FVmethod.FVmethod(recGrid.tinMesh['vertices'], recGrid.tinMesh['triangles'],
                               recGrid.tinMesh['edges'])

    # Perform partitioning by equivalent domain splitting
    partitionIDs, RowProc, ColProc = partitionTIN.simple(recGrid.tinMesh['vertices'][:, 0],
                                                         recGrid.tinMesh['vertices'][:, 1])
    FVmesh.partIDs = partitionIDs

    # Get each partition global node ID
    inGIDs = np.where(partitionIDs == 0)[0]

    if verbose:
        print(" - partition TIN amongst processors ", time.clock() - walltime)

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

    if verbose:
        print(" - reconstructed FV mesh ", time.clock() - walltime)

    inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == 0)[0]
    inIDs += recGrid.boundsPt
    elevationTIN.assign_parameter_pit(FVmesh.neighbours, FVmesh.control_volumes, input.diffnb,
                                      input.diffprop, input.propa, input.propb, recGrid.boundsPt,
                                      input.fillmax)

    return FVmesh, tMesh, lGIDs, inIDs, inGIDs, totPts

def _define_TINparams(totPts, input, FVmesh, recGrid, verbose=False):
    """
    This function is defining the main values declared on the TIN.
    """

    walltime = time.clock()

    inIDs = np.where(FVmesh.partIDs[recGrid.boundsPt:] == 0)[0]
    inIDs += recGrid.boundsPt

    local_elev = np.zeros(totPts)
    local_elev.fill(-1.e6)

    # In case of a restart read values from HDF5 files
    if input.restart:
        local_cum = np.zeros(totPts)
        local_cum.fill(-1.e6)
        local_hill = np.zeros(totPts)
        local_hill.fill(-1.e6)
        local_fail = np.zeros(totPts)
        local_fail.fill(-1.e6)
        if input.flexure:
            local_cumflex = np.zeros(totPts)
            local_cumflex.fill(-1.e6)
            local_elev[inIDs],local_cum[inIDs],local_hill[inIDs], local_fail[inIDs], local_cumflex[inIDs] = recGrid.load_hdf5_flex(input.rfolder,
                                                                input.rstep,FVmesh.node_coords[inIDs, :2])
        else:
            local_elev[inIDs],local_cum[inIDs],local_hill[inIDs], local_fail[inIDs] = recGrid.load_hdf5(input.rfolder,input.rstep,
                                                                FVmesh.node_coords[inIDs, :2])

        # Get cumulative erosion/deposition values
        cumdiff = local_cum
        cumdiff[:recGrid.boundsPt] = 0.
        cumhill = local_hill
        cumhill[:recGrid.boundsPt] = 0.
        cumfail = local_hill
        cumfail[:recGrid.boundsPt] = 0.
        if input.flexure:
            # Get cumulative flexural values
            cumflex = local_cumflex
            cumflex[:recGrid.boundsPt] = 0.

    # Otherwise interpolate elevation from DEM to TIN
    else:
        local_elev[inIDs] = elevationTIN.getElevation(recGrid.regX, recGrid.regY,
                                            recGrid.regZ, FVmesh.node_coords[inIDs, :2])
        # Initialise TIN parameters
        cumdiff = np.zeros(totPts)
        cumhill = np.zeros(totPts)
        cumfail = np.zeros(totPts)
        if input.flexure:
            cumflex = np.zeros(totPts)

    # Assign boundary values
    elevation, parentIDs = elevationTIN.update_border_elevation(local_elev, FVmesh.neighbours,
                                FVmesh.edge_length, recGrid.boundsPt, btype=input.btype)

    # Define pit filling algorithm
    elevationTIN.assign_parameter_pit(FVmesh.neighbours, FVmesh.control_volumes, input.diffnb,
                                      input.diffprop, input.propa, input.propb,
                                      recGrid.boundsPt, input.fillmax)

    if verbose:
        print(" - define paramters on TIN grid ", time.clock() - walltime)

    if input.flexure:
        return elevation, cumdiff, cumhill, cumfail, cumflex, inIDs, parentIDs
    else:
        return elevation, cumdiff, cumhill, cumfail, inIDs, parentIDs

def _build_strateroMesh(input, FVmesh, recGrid, cumdiff, verbose=False):
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
                                input.outDir, input.sh5file, input.poro0, input.poroC, cumdiff, input.rfolder, input.rstep)
        else:
            strata = strataMesh.strataMesh(sdx, bbX, bbY, layNb, FVmesh.node_coords[:, :2],
                                input.outDir, input.sh5file, input.poro0, input.poroC)
        if verbose:
            print(" - create stratigraphic regions ", time.clock() - walltime)

    # Define pre-existing erodibility maps
    if input.erolays and input.erolays >= 0:

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

        if verbose:
            print(" - create erodibility mesh ", time.clock() - walltime)

    if ((input.laytime and input.laytime > 0) and
       (input.erolays and input.erolays >= 0)):
        return strata, mapero
    elif (input.laytime and input.laytime > 0):
        return strata
    else:
        return mapero

def _init_flexure(FVmesh, input, recGrid, force, elevation, cumdiff,
                  cumflex, totPts, verbose=False):
    """
    Initialise flexural isostasy.
    """

    # Initialise flexure parameters for gFlex library.
    walltime = time.clock()

    nx = input.fnx
    ny = input.fny
    elasticT2 = None
    if nx == 0:
        nx = recGrid.nx
    if ny == 0:
        ny = recGrid.ny
    if input.elasticH is not None:
        elasticT = input.elasticH
    elif input.elasticGrid is not None:
        elasticT = str(input.elasticGrid)
    else:
        elasticT = input.elasticA1
        elasticT2 = input.elasticA2

    flex = isoFlex.isoFlex()
    flex.buildGrid(nx, ny, input.youngMod, input.dmantle, input.dsediment,
                elasticT, elasticT2, input.flexbounds, FVmesh.node_coords[:,:2], input.ftime)

    tinFlex = np.zeros(totPts, dtype=float)
    ref_elev = get_reference_elevation(input,recGrid,elevation)
    force.getSea(input.tStart,input.udw,ref_elev)
    tinFlex = flex.get_flexure(elevation, cumdiff, force.sealevel,
                               recGrid.boundsPt, initFlex=True)
    tinFlex = force.disp_border(tinFlex, FVmesh.neighbours, FVmesh.edge_length, recGrid.boundsPt)
    cumflex += tinFlex
    if verbose:
        print("   - Initialise flexural isostasy ", time.clock() - walltime)

    return flex, tinFlex, cumflex

def _init_wavesed(input, z0,recGrid, force, verbose=False):
    """
    Initialise wavesed mesh.
    """

    # Initialise wavesed parameters.
    walltime = time.clock()

    wave = waveSed.waveSed(input, recGrid, Ce=input.wCe, Cd=50.)
    force.getSea(input.tStart,input.udw,z0)

    if verbose:
        print("   - Initialise wavesed grid ", time.clock() - walltime)

    return wave


def get_reference_elevation(input, recGrid, elevation):
    """
    The following function define the elevation from the TIN to a regular grid...

    Args:
        input: class containing XML input file parameters.
        recGrid: class describing the regular grid characteristics.
        elevation: TIN elevation mesh.

    Returns:
        - ref_elev - interpolated elevation on the regular grid
    """
    if input.searef:
        x_ref, y_ref = input.searef
        pts = recGrid.tinMesh["vertices"]
        ref_elev = griddata(points=pts, values=elevation, xi=[x_ref, y_ref], method="nearest")
    else:
        ref_elev = 0.

    return ref_elev
