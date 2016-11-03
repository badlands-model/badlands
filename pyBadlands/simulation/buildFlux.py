##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This file is the main entry point to compute flow network and associated sedimentary fluxes.
"""

import time
import numpy as np
import mpi4py.MPI as mpi

from pyBadlands import (elevationTIN)

def streamflow(input, FVmesh, recGrid, force, hillslope, flow, elevation, \
                 lGIDs, rain, tNow, verbose=False):
    """
    Compute flow network.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD

    flow.Flow_time = time.clock()

    # Update sea-level
    walltime = time.clock()
    force.getSea(tNow)
    fillH = None

    # Update river input
    force.getRivers(tNow)
    riverrain = rain+force.rivQw

    if input.depo == 0 or input.capacity or input.filter:
        flow.maxdep = 0.
        flow.maxh = 0.
        if flow.esmooth is None and input.filter:
            flow.esmooth = 0
        if flow.dsmooth is None and input.filter:
            flow.dsmooth = 0
        # Build an initial depression-less surface at start time if required
        if input.tStart == tNow and input.nopit == 1:
            sea_lvl =  force.sealevel - input.sealimit
            elevation = elevationTIN.pit_stack_PD(elevation,sea_lvl,input.nopit)
            fillH = elevation
    else:
        # fillH = elevationTIN.pit_filling_PD(elevation, FVmesh.neighbours,
        #                            recGrid.boundsPt, force.sealevel-input.sealimit
        #                            input.fillmax)
        sea_lvl =  force.sealevel - input.sealimit
        #fillH = elevationTIN.pit_stack_PD(elevation,sea_lvl)
        # Build an initial depression-less surface at start time if required
        if input.tStart == tNow and input.nopit == 1 :
            fillH = elevationTIN.pit_stack_PD(elevation,sea_lvl,input.nopit)
            elevation = fillH
        else:
            fillH = elevationTIN.pit_stack_PD(elevation,sea_lvl,0)

    if rank == 0 and verbose and input.spl and not input.filter:
        print " -   depression-less algorithm PD with stack", time.clock() - walltime

    # Compute stream network
    walltime = time.clock()
    if input.nHillslope:
        flow.SFD_nreceivers(hillslope.Sc, fillH, elevation, FVmesh.neighbours,
                            FVmesh.vor_edges, FVmesh.edge_length,
                            lGIDs, force.sealevel-input.sealimit)
    else:
        flow.SFD_receivers(fillH, elevation, FVmesh.neighbours,
                           FVmesh.vor_edges, FVmesh.edge_length,
                           lGIDs, force.sealevel-input.sealimit)

    if rank == 0 and verbose:
        print " -   compute receivers parallel ", time.clock() - walltime

    # Distribute evenly local minimas to processors
    walltime = time.clock()
    flow.localbase = np.array_split(flow.base, size)[rank]
    flow.ordered_node_array()
    if rank == 0 and verbose:
        print " -   compute stack order locally ", time.clock() - walltime

    walltime = time.clock()
    stackNbs = comm.allgather(len(flow.localstack))
    globalstack = np.zeros(sum(stackNbs), dtype=flow.localstack.dtype)
    comm.Allgatherv(sendbuf=[flow.localstack, mpi.INT],
                    recvbuf=[globalstack, (stackNbs, None), mpi.INT])
    flow.stack = globalstack
    if rank == 0 and verbose:
        print " -   send stack order globally ", time.clock() - walltime

    # Compute discharge
    walltime = time.clock()
    flow.compute_flow(FVmesh.control_volumes, riverrain)
    if rank == 0 and verbose:
        print " -   compute discharge ", time.clock() - walltime

    return fillH, elevation

def sediment_flux(input, recGrid, hillslope, FVmesh, tMesh, flow, force, applyDisp, \
                  mapero, cumdiff, fillH, disp, inGIDs, elevation, tNow, tEnd, verbose=False):
    """
    Compute sediment fluxes.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD

    # Compute CFL condition
    walltime = time.clock()
    if input.Hillslope:
        hillslope.dt_stability(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
    elif input.nHillslope:
        hillslope.dt_stability(flow.diff_cfl)
    else:
        hillslope.CFL = tEnd-tNow

    if input.depo == 1:
        if input.filter:
            flow.CFL = input.maxDT
        elif input.spl:
            flow.dt_stability(fillH, inGIDs)
        else:
            flow.dt_stability(elevation, inGIDs)
    else:
        if input.filter:
            flow.CFL = input.maxDT
        else:
            flow.dt_stability(elevation, inGIDs)

    CFLtime = min(flow.CFL, hillslope.CFL)
    CFLtime = min(CFLtime,tEnd-tNow)
    CFLtime = max(input.minDT, CFLtime)
    CFLtime = min(input.maxiDT, CFLtime)

    if rank == 0 and verbose:
        print " -   Get CFL time step ", time.clock() - walltime

    # Compute sediment fluxes
    # Initial cumulative elevation change
    walltime = time.clock()
    diff_flux = hillslope.sedflux(flow.diff_flux, force.sealevel, elevation, FVmesh.control_volumes)
    if rank == 0 and verbose:
        print " -   Get hillslope fluxes ", time.clock() - walltime

    walltime = time.clock()
    xyMin = [recGrid.regX.min(), recGrid.regY.min()]
    xyMax = [recGrid.regX.max(), recGrid.regY.max()]
    id = np.where(force.rivQs>0)
    tmp = force.rivQs[id]
    tstep, sedrate = flow.compute_sedflux(FVmesh.control_volumes, elevation, fillH, xyMin, xyMax,
                                          diff_flux, CFLtime, force.rivQs, force.sealevel, cumdiff,
                                          input.perc_dep, input.slp_cr)
    if rank == 0 and verbose:
        print " -   Get stream fluxes ", time.clock() - walltime

    # Update surface parameters and time
    timestep = min(tstep, tEnd-tNow)
    diff = sedrate * timestep

    if input.filter:
        smthdiff = flow.gaussian_filter(diff)
        smthdiff[:recGrid.boundsPt] = 0.
        elevation += smthdiff
        cumdiff += smthdiff
    else:
        elevation += diff
        cumdiff += diff

    if applyDisp:
        elevation += disp * timestep

    # Update erodibility values
    if input.erolays >= 0:
        mapero.getErodibility(diff)
        flow.erodibility = mapero.erodibility

    tNow += timestep

    if rank == 0 and verbose:
        print " - Flow computation ", time.clock() - flow.Flow_time

    return tNow,elevation,cumdiff
