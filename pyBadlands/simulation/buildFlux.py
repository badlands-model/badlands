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
from matplotlib import path

from pyBadlands import (elevationTIN)

def streamflow(input, FVmesh, recGrid, force, hillslope, flow, elevation, \
                 lGIDs, rain, tNow, verbose=False):
    """
    Compute flow network.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD

    # Update sea-level
    walltime = time.clock()
    force.getSea(tNow)
    fillH = None

    # Update river input
    force.getRivers(tNow)
    riverrain = rain+force.rivQw

    # Build an initial depression-less surface at start time if required
    if input.tStart == tNow and input.nopit == 1 :
        fillH = elevationTIN.pit_stack_PD(elevation,input.nopit,force.sealevel)
        elevation = fillH
    else:
        fillH = elevationTIN.pit_stack_PD(elevation,0,force.sealevel)

    if rank == 0 and verbose and input.spl:
        print " -   depression-less algorithm PD with stack", time.clock() - walltime

    # Compute stream network
    walltime = time.clock()
    flow.SFD_receivers(fillH, elevation, FVmesh.neighbours,
                       FVmesh.vor_edges, FVmesh.edge_length,
                       lGIDs)

    if rank == 0 and verbose:
        print " -   compute receivers parallel ", time.clock() - walltime

    # Distribute evenly local minimas to processors on filled surface
    walltime = time.clock()
    flow.localbase = np.array_split(flow.base, size)[rank]
    flow.ordered_node_array_filled()
    if rank == 0 and verbose:
        print " -   compute stack order locally for filled surface", time.clock() - walltime

    walltime = time.clock()
    stackNbs = comm.allgather(len(flow.localstack))
    globalstack = np.zeros(sum(stackNbs), dtype=flow.localstack.dtype)
    comm.Allgatherv(sendbuf=[flow.localstack, mpi.INT],
                    recvbuf=[globalstack, (stackNbs, None), mpi.INT])
    flow.stack = globalstack
    if rank == 0 and verbose:
        print " -   send stack order for filled surface globally ", time.clock() - walltime

    # Distribute evenly local minimas on real surface
    walltime = time.clock()
    flow.localbase1 = np.array_split(flow.base1, size)[rank]
    flow.ordered_node_array_elev()
    if rank == 0 and verbose:
        print " -   compute stack order locally for real surface", time.clock() - walltime

    walltime = time.clock()
    stackNbs1 = comm.allgather(len(flow.localstack1))
    globalstack1 = np.zeros(sum(stackNbs1), dtype=flow.localstack1.dtype)
    comm.Allgatherv(sendbuf=[flow.localstack1, mpi.INT],
                    recvbuf=[globalstack1, (stackNbs1, None), mpi.INT])
    flow.stack1 = globalstack1
    if rank == 0 and verbose:
        print " -   send stack order for real surface globally ", time.clock() - walltime

    # Compute a unique ID for each local depression and their downstream draining nodes
    flow.compute_parameters_depression(fillH,elevation,FVmesh.control_volumes,force.sealevel)

    # Compute discharge
    walltime = time.clock()
    flow.compute_flow(elevation, FVmesh.control_volumes, riverrain)
    if rank == 0 and verbose:
        print " -   compute discharge ", time.clock() - walltime

    return fillH, elevation

def sediment_flux(input, recGrid, hillslope, FVmesh, tMesh, flow, force, rain, lGIDs, applyDisp, straTIN, \
                  mapero, cumdiff, cumhill, fillH, disp, inGIDs, elevation, tNow, tEnd, verbose=False):
    """
    Compute sediment fluxes.
    """

    rank = mpi.COMM_WORLD.rank
    size = mpi.COMM_WORLD.size
    comm = mpi.COMM_WORLD
    flow_time = time.clock()

    # Get active layer
    if straTIN is not None:
        walltime = time.clock()
        flow.activelay[flow.activelay<1.] = 1.
        flow.activelay[flow.activelay>straTIN.activeh] = straTIN.activeh
        straTIN.get_active_layer(flow.activelay)
        activelay = straTIN.alayR
        flow.straTIN = 1
        # Set the average erodibility based on rock types in the
        # active layer
        th = np.sum(activelay,axis=1).reshape(len(elevation),1)
        flow.erodibility = np.sum(straTIN.rockCk*activelay/th,axis=1)
        eroCk = straTIN.rockCk
        if rank == 0 and verbose:
            print " -   Get active layer ", time.clock() - walltime
    else:
        activelay = None
        eroCk = 0.

    # Find border/inside nodes
    ids = np.arange(len(FVmesh.control_volumes))
    tmp1 = np.where(FVmesh.control_volumes>0.)[0]
    xyMin = [recGrid.regX.min()-1., recGrid.regY.min()-1.]
    xyMax = [recGrid.regX.max()+1., recGrid.regY.max()+1.]
    domain = path.Path([(xyMin[0],xyMin[1]),(xyMax[0],xyMin[1]), (xyMax[0],xyMax[1]), (xyMin[0],xyMax[1])])
    tmp2 = domain.contains_points(flow.xycoords)
    insideIDs = np.intersect1d(tmp1,ids[tmp2])
    borders = np.zeros(len(FVmesh.control_volumes),dtype=int)
    borders[insideIDs] = 1
    outsideIDs = np.where(borders==0)[0]

    # Compute CFL condition
    walltime = time.clock()
    if input.Hillslope:
        hillslope.dt_stability(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
        hillslope.dt_stability_ms(FVmesh.edge_length[inGIDs,:tMesh.maxNgbh])
    else:
        hillslope.CFL = tEnd-tNow
    flow.dt_stability(fillH, inGIDs)
    CFLtime = min(flow.CFL, hillslope.CFL)
    if CFLtime>1.:
        CFLtime = float(round(CFLtime-0.5,0))
    if rank == 0 and verbose:
        print 'CFL for hillslope and flow ',hillslope.CFL,flow.CFL,CFLtime
    CFLtime = min(CFLtime, tEnd - tNow)
    CFLtime = max(input.minDT, CFLtime)
    CFLtime = min(input.maxDT, CFLtime)
    if rank == 0 and verbose:
        print " -   Get CFL time step ", time.clock() - walltime

    # Compute sediment fluxes
    # Initial cumulative elevation change
    walltime = time.clock()
    ids = np.where(force.rivQs>0)
    tmp = force.rivQs[ids]

    timestep, sedchange, erosion, deposition = flow.compute_sedflux(FVmesh.control_volumes, elevation, rain, fillH,
                                          borders, domain, CFLtime, activelay, eroCk, force.rivQs, force.sealevel,
                                          input.perc_dep, input.slp_cr, FVmesh.neighbours, verbose)
    if rank == 0 and verbose:
        print " -   Get stream fluxes ", time.clock() - walltime
    ed = np.sum(sedchange,axis=1)
    elevation += ed
    cumdiff += ed

    # Compute marine sediment diffusion
    if hillslope.CDriver > 0.:
        walltime = time.clock()

        # Initialise marine sediments diffusion array
        it = 0
        diffstep = timestep
        frac = np.zeros(deposition.shape, order='F')
        diffsed = np.zeros(len(elevation))
        seaIDs = np.zeros(len(elevation),dtype=int)

        # Perform river related sediment diffusion
        while diffstep > 0. and it < 100:
            sumdep = np.sum(deposition,axis=1)
            tmp = np.where(np.logical_and(sumdep>0.,elevation<force.sealevel))[0]
            if it == 0:
                diffsed[tmp] = sumdep[tmp]
            # Only diffuse sediment pile which are above 1 m thick
            tmpIDs = np.where(diffsed>1.)[0]
            seaIDs[tmpIDs] = 1
            # Compute marine fluxes
            flow.compute_marine_diffusion(elevation, borders, seaIDs, FVmesh.neighbours,
                           FVmesh.vor_edges, FVmesh.edge_length, lGIDs)
            diffcoeff = hillslope.sedfluxmarine(force.sealevel, elevation,
                                                 FVmesh.control_volumes)
            diffcoeff[outsideIDs] = 0.
            diffmarine = diffcoeff*flow.diff_flux

            # Define maximum time step
            maxstep = min(hillslope.CFLms,diffstep)
            # Check for excessive erosion thicknesses
            limitIDs = np.where(np.logical_and(-diffmarine*maxstep>sumdep,diffmarine<0.))[0]
            if len(limitIDs) > 0:
                mindt = -sumdep[limitIDs]/diffmarine[limitIDs]
                maxstep = max(mindt.min(),input.minDT)
            diffstep -= maxstep

            # Distribute rock based on their respective proportions in the deposited columns
            if straTIN is not None:
                # Get proportion of each rock contained in the marine deposits
                frac.fill(0.)
                frac[tmp,:] = deposition[tmp,:]/sumdep[tmp].reshape(len(tmp),1)
                # Compute multi-rock diffusion
                sedpropflux = flow.compute_sediment_diffusion(elevation, borders, seaIDs, frac,
                                    diffcoeff*maxstep, FVmesh.neighbours, FVmesh.vor_edges,
                                    FVmesh.edge_length, lGIDs)
                # Update deposition
                deposition += sedpropflux
                deposition[deposition<0] = 0.

            # Update elevation, erosion/deposition and remaining river sediment to diffuse
            elevation += diffmarine*maxstep
            cumdiff += diffmarine*maxstep
            diffsed += diffmarine*maxstep
            if straTIN is None:
                deposition = diffsed.reshape(len(diffsed),1)
            seaIDs[tmpIDs] = 0
            it += 1
        if rank == 0 and verbose:
            print " -   Get river sediment marine fluxes ", time.clock() - walltime

    difftype = 0
    if straTIN is not None:
        walltime = time.clock()
        straTIN.update_erosion(erosion)
        straTIN.update_deposition(deposition)
        if rank == 0 and verbose:
            print " -   Update erosion layer ", time.clock() - walltime
        difftype = 1

    # Compute hillslope processes
    walltime = time.clock()
    flow.compute_hillslope_diffusion(elevation, borders, FVmesh.neighbours,
                       FVmesh.vor_edges, FVmesh.edge_length,lGIDs, difftype)
    diffcoeff = hillslope.sedflux(force.sealevel, elevation, FVmesh.control_volumes)
    cdiff = diffcoeff * flow.diff_flux * timestep
    if straTIN is None:
        diff_flux = np.zeros(len(cdiff))
        diff_flux[insideIDs] = cdiff[insideIDs]
        diff = diff_flux
        if input.btype == 'outlet':
            diff[insideIDs[0]] = 0.
        elevation += diff
    else:
        # Extract the thickness that will be eroded due to diffusion
        maxlayh = -cdiff
        maxlayh[maxlayh<1.] = 1.
        straTIN.get_active_layer(maxlayh)
        difflay = straTIN.alayR
        sumdep = np.sum(difflay,axis=1)
        tmp = np.where(sumdep>0.)[0]
        frac[tmp,:] = difflay[tmp,:]/sumdep[tmp].reshape(len(tmp),1)
        # Compute multi-rock diffusion
        tIDs = np.ones(len(elevation))
        sedpropflux = flow.compute_sediment_diffusion(elevation, borders, tIDs, frac,
                            diffcoeff*timestep, FVmesh.neighbours, FVmesh.vor_edges,
                            FVmesh.edge_length,lGIDs)
        if input.btype == 'outlet':
            sedpropflux[insideIDs[0],:] = 0.
        # Update elevation
        diff = np.sum(sedpropflux,axis=1)
        elevation += diff
        # Update active layer
        erosion.fill(0.)
        tmp = np.where(diff<0.)
        erosion[tmp,:] = sedpropflux[tmp,:]
        straTIN.update_erosion(erosion)
        deposition.fill(0.)
        tmp = np.where(diff>0.)
        deposition[tmp,:] = sedpropflux[tmp,:]
        straTIN.update_deposition(deposition)

    if input.btype == 'slope':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]-0.1
    elif input.btype == 'flat':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]
    elif input.btype == 'wall':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]+100.
    elif input.btype == 'outlet':
        elevation[1:len(flow.parentIDs)] = elevation[flow.parentIDs[1:]]+100.
    elif input.btype == 'wall1':
        elevation[:len(flow.parentIDs)] = elevation[flow.parentIDs]-0.1
        elevation[:recGrid.nx+1] = elevation[flow.parentIDs[:recGrid.nx+1]]+100.
    cumdiff += diff
    cumhill += diff
    if rank == 0 and verbose:
        print " -   Get hillslope fluxes ", time.clock() - walltime

    if applyDisp:
        elevation += disp * timestep

    # Update erodibility values
    if input.erolays >= 0:
        mapero.getErodibility(diff)
        flow.erodibility = mapero.erodibility

    tNow += timestep

    if rank == 0 and verbose:
        print " - Flow computation ", time.clock() - flow_time

    return tNow,elevation,cumdiff,cumhill
