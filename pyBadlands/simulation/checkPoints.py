##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file is used to checkpoint and is the main entry to write simulation output.
"""

import time
import numpy as np

from pyBadlands import (visualiseFlow, visualiseTIN, eroMesh)

def write_checkpoints(input, recGrid, lGIDs, inIDs, tNow, FVmesh, \
                      tMesh, force, flow, rain, elevation, fillH, \
                      cumdiff, cumhill, cumfail, wavediff, step, prop, \
                      mapero=None, cumflex=None):
    """
    Create the checkpoint files (used for HDF5 output).
    """

    deepb = input.deepbasin

    if (input.erolays and input.erolays >= 0):
        eroOn = True
    else:
        eroOn = False
    out_time = time.clock()
    visXlim = np.zeros(2)
    visYlim = np.zeros(2)
    visXlim[0] = recGrid.rectX.min()
    visXlim[1] = recGrid.rectX.max()
    visYlim[0] = recGrid.rectY.min()
    visYlim[1] = recGrid.rectY.max()

    # Done when TIN has been built/rebuilt
    if FVmesh.outPts is None and FVmesh.outCells is None:
        FVmesh.outPts, FVmesh.outCells = visualiseTIN.output_cellsIDs(lGIDs, inIDs,
                                            visXlim, visYlim, FVmesh.node_coords[:, :2],
                                            tMesh.cells)
    tcells = np.zeros(1)
    tcells[0] = len(FVmesh.outCells)
    tnodes = np.zeros(1)
    tnodes[0] = len(lGIDs)

    # Done for every visualisation step
    flowIDs, polylines = visualiseFlow.output_Polylines(FVmesh.outPts, flow.receivers[FVmesh.outPts],
                                                        visXlim, visYlim, FVmesh.node_coords[:, :2])
    fnodes = np.zeros(1)
    fnodes[0] = len(flowIDs)
    fline = np.zeros(1)
    fline[0] = len(polylines[:, 0])

    # Compute flow parameters
    if deepb >= 5000.:
        deepb = force.sealevel
    flow.view_receivers(fillH, elevation, FVmesh.neighbours, FVmesh.vor_edges,
                        FVmesh.edge_length, lGIDs, deepb) #force.sealevel)
    flow.compute_parameters()
    visdis = np.copy(flow.discharge)
    seaIDs = np.where(elevation<deepb)[0] #force.sealevel)[0]
    if len(seaIDs)>0:
        visdis[seaIDs] = 1.
        flow.basinID[seaIDs] = -1
    visdis[visdis<1.] = 1.

    rockOn = False
    if input.carbonate or input.pelagic:
        rockOn = True

    # Write HDF5 files
    if input.waveSed and tNow > input.tStart:
        waveOn = True
        meanH = force.meanH[lGIDs]
        meanS = force.meanS[lGIDs]
        wdiff = wavediff[lGIDs]
    else:
        waveOn = False
        meanH = None
        meanS = None
        wdiff = None

    if input.flexure:
        visualiseTIN.write_hdf5_flexure(input.outDir, input.th5file, step, tMesh.node_coords[:,:2],
                                    elevation[lGIDs], rain[lGIDs], visdis[lGIDs], cumdiff[lGIDs],
                                    cumhill[lGIDs], cumfail[lGIDs], cumflex[lGIDs], FVmesh.outCells, input.oroRain,
                                    eroOn, flow.erodibility[lGIDs], FVmesh.control_volumes[lGIDs],
                                    waveOn, meanH, meanS, wdiff, rockOn, prop[lGIDs,:])
    else:
        visualiseTIN.write_hdf5(input.outDir, input.th5file, step, tMesh.node_coords[:,:2],
                                elevation[lGIDs], rain[lGIDs], visdis[lGIDs], cumdiff[lGIDs],
                                cumhill[lGIDs], cumfail[lGIDs], FVmesh.outCells, input.oroRain, eroOn,
                                flow.erodibility[lGIDs], FVmesh.control_volumes[lGIDs],
                                waveOn, meanH, meanS, wdiff, rockOn,
                                prop[lGIDs,:])

    if flow.sedload is not None:
            if flow.flowdensity is not None:
                visualiseFlow.write_hdf5(input.outDir, input.fh5file, step, FVmesh.node_coords[flowIDs, :2], elevation[flowIDs],
                                     visdis[flowIDs], flow.chi[flowIDs], flow.sedload[flowIDs], flow.flowdensity[flowIDs], flow.basinID[flowIDs], polylines)
            else:
                zeros = np.zeros(len(flowIDs),dtype=float)
                visualiseFlow.write_hdf5(input.outDir, input.fh5file, step, FVmesh.node_coords[flowIDs, :2], elevation[flowIDs],
                                     visdis[flowIDs], flow.chi[flowIDs], flow.sedload[flowIDs], zeros, flow.basinID[flowIDs], polylines)
    else:
        zeros = np.zeros(len(flowIDs),dtype=float)
        if flow.flowdensity is not None:
            visualiseFlow.write_hdf5(input.outDir, input.fh5file, step, FVmesh.node_coords[flowIDs, :2], elevation[flowIDs],
                                 visdis[flowIDs], flow.chi[flowIDs], zeros, flow.flowdensity[flowIDs], flow.basinID[flowIDs], polylines)
        else:
            visualiseFlow.write_hdf5(input.outDir, input.fh5file, step, FVmesh.node_coords[flowIDs, :2], elevation[flowIDs],
                                 visdis[flowIDs], flow.chi[flowIDs], zeros, zeros, flow.basinID[flowIDs], polylines)

    # Combine HDF5 files and write time series
    visualiseTIN.write_xmf(input.outDir, input.txmffile, input.txdmffile, step, tNow, tcells,
                           tnodes, input.th5file, force.sealevel, input.flexure,
                           input.oroRain, eroOn, waveOn, rockOn, prop.shape[1])

    visualiseFlow.write_xmf(input.outDir, input.fxmffile, input.fxdmffile, step, tNow,
                            fline, fnodes, input.fh5file)

    print("   - Writing outputs (%0.02f seconds; tNow = %s)" % (time.clock() - out_time, tNow))

    # Record erodibility maps
    if (input.erolays and input.erolays >= 0):
        mapero.write_hdf5_erolay(step)
