##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file saves simulation information to allow for checkpointing. It not only stores
some data for restarting a new simulation but also is the main entry to write the
simulation output.

**Badlands** uses `hdf5`_ library for outputs generation.
To read the **hdf5** files, the code creates a *XML schema* (**.xmf** files) describing
how the data is stored in the **hdf5** files.

Then `xdmf`_ files provides support to visualise temporal evolution of the output with
applications such as **Paraview** or Visit.

.. image:: img/output.png
   :scale: 50 %
   :alt: Badlands output folder
   :align: center

Above is a typical structure that you will find in **badlands** output folder:

- **h5** folder contains the **hdf5** data, all the information computed by the model are stored in these files. You will have at least the *tin* (surface) and *flow* (stream network) dataset and also the *sed* (stratigraphy) data if the stratal structure is computed in your simulation.
- **xmf** folder contains the XML files used to read the **hdf5** files contained in the **h5** folder.
- the **.xml** input file used to build this specific model.
- two **.xdmf** files for the surface (**tin_series.xdmf**) and the flow network (**flow_series.xdmf**) that read the **xmf** files through time.

"""

import time
import numpy as np
import os
if 'READTHEDOCS' not in os.environ:
    from badlands import (visualiseFlow, visualiseTIN, eroMesh)

def write_checkpoints(input, recGrid, lGIDs, inIDs, tNow, FVmesh, \
                      tMesh, force, flow, rain, elevation, fillH, \
                      cumdiff, cumhill, cumfail, wavediff, step, prop, \
                      mapero=None, cumflex=None):
    """
    Create the checkpoint files (used for HDF5 output).

    Args:
        input: class containing XML input file parameters.
        recGrid: class describing the regular grid characteristics.
        lGIDs: global nodes indices.
        inIDs: local nodes indices.
        tNow: current time step
        FVmesh: finite volume mesh
        tMesh: TIN mesh
        force:  forcing parameters
        flow: flow parameters
        rain: rain value
        elevation: elevation mesh
        fillH: filled elevation mesh
        cumdiff: cumulative total erosion/deposition changes
        cumhill: cumulative hillslope erosion/deposition changes
        cumfail: cumulative failure induced erosion/deposition changes
        wavediff: cumulative wave induced erosion/deposition
        step: iterative step for output
        prop: proportion of each sediment type
        mapero: imposed erodibility maps
        cumflex: cumlative changes induced by flexural isostasy
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
    flow.compute_parameters(elevation,force.sealevel)
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
                                    waveOn, meanH, meanS, wdiff, rockOn, prop[lGIDs,:], force.sealevel)
    else:
        visualiseTIN.write_hdf5(input.outDir, input.th5file, step, tMesh.node_coords[:,:2],
                                elevation[lGIDs], rain[lGIDs], visdis[lGIDs], cumdiff[lGIDs],
                                cumhill[lGIDs], cumfail[lGIDs], FVmesh.outCells, input.oroRain, eroOn,
                                flow.erodibility[lGIDs], FVmesh.control_volumes[lGIDs],
                                waveOn, meanH, meanS, wdiff, rockOn,
                                prop[lGIDs,:], force.sealevel)

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

    return
