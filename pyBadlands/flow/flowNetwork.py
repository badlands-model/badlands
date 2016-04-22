##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to Badlands stream network computation.
"""

import math
import numpy
import warnings
import mpi4py.MPI as mpi
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter

from pyBadlands.libUtils import SFDalgo as SFD
from pyBadlands.libUtils import FLOWalgo
from pyBadlands.libUtils import FLWnetwork
from pyBadlands.libUtils import PDalgo

class flowNetwork:
    """
    Class for handling flow network computation based on Braun & Willett's
    algorithm.
    """

    def __init__(self):
        """
        Initialization.
        """

        self.xycoords = None
        self.base = None
        self.localbase = None
        self.receivers = None
        self.arrdonor = None
        self.delta = None
        self.donors = None
        self.localstack = None
        self.stack = None
        self.partFlow = None
        self.maxdonors = 0
        self.CFL = None
        self.erodibility = None
        self.m = None
        self.n = None
        self.mindt = None
        self.alluvial = 0.
        self.bedrock = 0.
        self.esmooth = 0.
        self.dsmooth = 0.
        self.spl = False
        self.capacity = False
        self.filter = False
        self.depo = 0

        self.discharge = None
        self.localsedflux = None
        self.maxh = None
        self.maxdep = None
        self.diff_flux = None
        self.chi = None
        self.basinID = None

        self.xgrid = None
        self.xgrid = None
        self.xi = None
        self.yi = None
        self.xyi = None

        return

    def SFD_receivers(self, fillH, elev, neighbours, edges, distances, globalIDs, sea):
        """
        Single Flow Direction function computes downslope flow directions by inspecting the neighborhood
        elevations around each node. The SFD method assigns a unique flow direction towards the steepest
        downslope neighbor.

        Parameters
        ----------
        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable : neighbours
            Numpy integer-type array with the neighbourhood IDs.

        variable : edges
            Numpy real-type array with the voronoi edges length for each neighbours of local IDs.

        variable : distances
            Numpy real-type array with the distances between each connection in the local TIN.

        variable: globalIDs
            Numpy integer-type array containing for local nodes their global IDs.

        variable : sea
            Current elevation of sea level.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Call the SFD function from libUtils
        if self.depo == 0 or self.capacity or self.filter:
            base, receivers, diff_flux = SFD.sfdcompute.directions_base(elev, \
                neighbours, edges, distances, globalIDs, sea)

            # Send local base level globally
            comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)
            bpos = numpy.where(base >= 0)[0]
            self.base = numpy.random.shuffle(bpos)
            
            # Send local receivers globally
            comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local diffusion flux globally
            comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux

        else:
            base, receivers, maxh, maxdep, diff_flux = SFD.sfdcompute.directions(fillH, elev, \
                neighbours, edges, distances, globalIDs, sea)

            # Send local base level globally
            comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)
            bpos = numpy.where(base >= 0)[0]
            self.base = numpy.random.shuffle(bpos)

            # Send local receivers globally
            comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local maximum deposition globally
            comm.Allreduce(mpi.IN_PLACE,maxh,op=mpi.MAX)
            self.maxh = maxh

            # Send local maximum deposition globally
            comm.Allreduce(mpi.IN_PLACE,maxdep,op=mpi.MAX)
            self.maxdep = maxdep

            # Send local diffusion flux globally
            comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux

        return

    def _donors_number_array(self):
        """
        Creates an array containing the number of donors for each node.
        """

        self.arrdonor = None
        numPts = len(self.receivers)
        self.arrdonor = numpy.zeros(numPts, dtype=int)
        maxID = numpy.max(self.receivers)
        self.arrdonor[:(maxID+1)] = numpy.bincount(self.receivers)
        self.maxdonors = self.arrdonor.max()

        return

    def _delta_array(self):
        """
        Creates the "delta" array, which is a list containing, for each
        node, the array index where that node's donor list begins.
        """

        self.delta = None
        nbdonors = len(self.arrdonor)
        self.delta = numpy.zeros( nbdonors+1, dtype=int)
        self.delta.fill(nbdonors)
        self.delta[-2::-1] -= numpy.cumsum(self.arrdonor[::-1])

        return

    def ordered_node_array(self):
        """
        Creates an array of node IDs that is arranged in order from downstream
        to upstream.
        """

        # Build donors array for each node
        self._donors_number_array()

        # Create the delta array
        self._delta_array()

        # Using libUtils stack create the ordered node array
        self.donors,lstcks = FLWnetwork.fstack.build(self.localbase,self.receivers,self.delta)

        # Create local stack
        stids = numpy.where(lstcks > -1 )[0]
        self.localstack = lstcks[stids]

        return

    def compute_flow(self, Acell, rain, parallel=False):
        """
        Calculates the drainage area and water discharge at each node.

        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        variable : rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).

        variable : parallel
            Boolean to inform if the discharge algorithm need to be run in parallel (True) or serial (False).
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        numPts = len(Acell)

        self.discharge = numpy.zeros(numPts, dtype=float)
        self.discharge[self.stack] = Acell[self.stack]*rain[self.stack]

        # Compute discharge using libUtils
        if(parallel):
            discharge = FLOWalgo.flowcompute.discharge(self.localstack,self.receivers,self.discharge)
            comm.Allreduce(mpi.IN_PLACE,discharge,op=mpi.MAX)
        else:
            discharge = FLOWalgo.flowcompute.discharge(self.stack,self.receivers,self.discharge)

        self.discharge = discharge

        return

    def compute_parameters(self):
        """
        Calculates the catchment IDs and the Chi parameter (Willett 2014).
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Get basin starting IDs for each local partition
        cumbase = numpy.zeros(size+1)
        for i in range(size):
            cumbase[i+1] = len(numpy.array_split(self.base, size)[i])+cumbase[i]+1

        # Compute discharge using libUtils
        splexp = self.m / self.n
        chi, basinID = FLOWalgo.flowcompute.parameters(self.localstack,self.receivers,
                                               self.discharge,self.xycoords,splexp,cumbase[rank])
        comm.Allreduce(mpi.IN_PLACE,chi,op=mpi.MAX)
        comm.Allreduce(mpi.IN_PLACE,basinID,op=mpi.MAX)

        self.chi = chi
        self.basinID = basinID

        return

    def compute_sedflux(self, Acell, elev, fillH, xymin, xymax, diff_flux, dt, sealevel, cumdiff):
        """
        Calculates the sediment flux at each node.

        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        variable : xymin
            Numpy array containing the minimal XY coordinates of TIN grid (excuding border cells).

        variable : xymax
            Numpy array containing the maximal XY coordinates of TIN grid (excuding border cells).

        variable : diff_flux
            Numpy arrays representing the fluxes due to linear diffusion equation.

        variable : dt
            Real value corresponding to the maximal stability time step.

        variable : sealevel
            Real value giving the sea-level height at considered time step.

        variable : cumdiff
            Numpy array containing the cumulative deposit thicknesses.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Compute sediment flux using libUtils
        if(size > 1):
            # Purely erosive case
            if self.spl and self.depo == 0:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_ero_only(self.localstack,self.receivers, \
                                      self.xycoords,xymin,xymax,self.discharge,elev, \
                                      diff_flux,self.erodibility,self.m,self.n,sealevel,dt)

            # Stream power law and mass is not conserved
            elif self.spl and self.filter:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_nocapacity_quick(self.localstack,self.receivers, \
                         self.xycoords,Acell,xymin,xymax,self.discharge,elev,diff_flux,self.erodibility, \
                         self.m,self.n,sealevel,dt)

            # Stream power law
            elif self.spl:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_nocapacity(self.localstack,self.receivers,self.xycoords, \
                         Acell,xymin,xymax,self.maxh,self.maxdep,self.discharge,fillH,elev,diff_flux, \
                         self.erodibility,self.m,self.n,sealevel,dt)

            # River carrying capacity case
            else:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_capacity(self.localstack,self.receivers,self.xycoords,\
                         Acell,xymin,xymax,self.discharge,elev,diff_flux,cumdiff,self.erodibility, \
                         self.m,self.n,self.bedrock,self.alluvial,sealevel,dt)

            timestep = numpy.zeros(1)
            timestep[0] = newdt
            comm.Allreduce(mpi.IN_PLACE,timestep,op=mpi.MIN)
            newdt = timestep[0]
            comm.Allreduce(mpi.IN_PLACE,sedflux,op=mpi.MAX)
            tempIDs = numpy.where(sedflux < -9.5e5)
            sedflux[tempIDs] = 0.
            newdt = max(self.mindt,newdt)
            sedrate = sedflux
        else:
            # Purely erosive case
            if self.spl and self.depo == 0:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_ero_only(self.localstack,self.receivers, \
                                      self.xycoords,xymin,xymax,self.discharge,elev, \
                                      diff_flux,self.erodibility,self.m,self.n,sealevel,dt)

            # Stream power law and mass is not conserved
            elif self.spl and self.filter:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_nocapacity_quick(self.localstack,self.receivers, \
                         self.xycoords,Acell,xymin,xymax,self.discharge,elev,diff_flux,self.erodibility, \
                         self.m,self.n,sealevel,dt)

            # Stream power law
            elif self.spl:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_nocapacity(self.localstack,self.receivers,self.xycoords, \
                         Acell,xymin,xymax,self.maxh,self.maxdep,self.discharge,fillH,elev,diff_flux, \
                         self.erodibility,self.m,self.n,sealevel,dt)

            # River carrying capacity case
            else:
                sedflux, newdt = FLOWalgo.flowcompute.sedflux_capacity(self.localstack,self.receivers,self.xycoords,\
                         Acell,xymin,xymax,self.discharge,elev,diff_flux,cumdiff,self.erodibility, \
                         self.m,self.n,self.bedrock,self.alluvial,sealevel,dt)

            tempIDs = numpy.where(sedflux < -9.5e5)
            sedflux[tempIDs] = 0.
            newdt = max(self.mindt,newdt)
            sedrate = sedflux

        return newdt,sedrate

    def gaussian_filter(self, diff):
        """
        Gaussian filter operation used to smooth erosion and deposition
        thicknesses for large simulation time steps. Using this operation
        implies that the resulting simulation is not conserving mass.

        Parameters
        ----------
        variable : diff
            Numpy arrays containing the erosion and deposition thicknesses.
        """

        if self.xgrid is None:
            dx = self.xycoords[1,0] - self.xycoords[0,0]
            xmin, xmax = min(self.xycoords[:,0]), max(self.xycoords[:,0])
            ymin, ymax = min(self.xycoords[:,1]), max(self.xycoords[:,1])
            self.xgrid = numpy.arange(xmin,xmax+dx,dx)
            self.ygrid = numpy.arange(ymin,ymax+dx,dx)
            self.xi, self.yi = numpy.meshgrid(self.xgrid, self.ygrid)
            self.xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]

        depZ = numpy.copy(diff)
        depZ = depZ.clip(0.)

        eroZ = numpy.copy(diff)
        eroZ = eroZ.clip(max=0.)

        tree = cKDTree(self.xycoords[:,:2])
        distances, indices = tree.query(self.xyi, k=3)
        if len(depZ[indices].shape) == 3:
            zd_vals = depZ[indices][:,:,0]
            ze_vals = eroZ[indices][:,:,0]
        else:
            zd_vals = depZ[indices]
            ze_vals = eroZ[indices]
        zdi = numpy.average(zd_vals,weights=(1./distances), axis=1)
        zei = numpy.average(ze_vals,weights=(1./distances), axis=1)

        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            zdi[onIDs] = depZ[indices[onIDs,0]]
            zei[onIDs] = eroZ[indices[onIDs,0]]

        depzi = numpy.reshape(zdi,(len(self.ygrid),len(self.xgrid)))
        erozi = numpy.reshape(zei,(len(self.ygrid),len(self.xgrid)))

        smthDep = gaussian_filter(depzi, sigma=self.dsmooth)
        smthEro = gaussian_filter(erozi, sigma=self.esmooth)

        rgi_dep = RegularGridInterpolator((self.ygrid, self.xgrid), smthDep)
        zdepsmth = rgi_dep((self.xycoords[:,1],self.xycoords[:,0]))
        rgi_ero = RegularGridInterpolator((self.ygrid, self.xgrid), smthEro)
        zerosmth = rgi_ero((self.xycoords[:,1],self.xycoords[:,0]))

        return zdepsmth + zerosmth

    def dt_fstability(self, elev, locIDs):
        """
        This pyfortran function computes the maximal timestep to ensure computation stability
        of the flow processes. This CFL-like condition is computed using erodibility
        coefficients, discharges plus elevations and distances between TIN nodes and
        their respective reveivers.

        Parameters
        ----------
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable: locIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Compute the local value for time stability
        dt = FLOWalgo.flowcompute.flowcfl(locIDs,self.receivers,self.xycoords,elev, \
                                      self.discharge,self.erodibility,self.m,self.n)

        # Global mimimum value for diffusion stability
        CFL = numpy.zeros(1)
        CFL[0] = dt
        comm.Allreduce(mpi.IN_PLACE,CFL,op=mpi.MIN)
        self.CFL = CFL[0]

    def dt_pstability(self, elev, locIDs):
        """
        This pure python function computes the maximal timestep to ensure computation stability
        of the flow processes. This CFL-like condition is computed using erodibility
        coefficients, discharges plus elevations and distances between TIN nodes and
        their respective reveivers.

        Parameters
        ----------
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable: locIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Compute elevation difference and distance between receiver and donor
        dz0 = elev[locIDs] - elev[self.receivers[locIDs]]
        dx = self.xycoords[locIDs,0] - self.xycoords[self.receivers[locIDs],0]
        dy = self.xycoords[locIDs,1] - self.xycoords[self.receivers[locIDs],1]
        dist0 = numpy.sqrt(dx**2 + dy**2)
        dshg0 = self.discharge[locIDs]

        # Limit to positive values of dz and discharge
        tmpID = numpy.where(dz0 > 0.)
        tmpdz = dz0[tmpID]
        tmpdist = dist0[tmpID]
        tmpdisc = dshg0[tmpID]
        tmpID2 = numpy.where( tmpdisc > 0.)
        discharge = tmpdisc[tmpID2]
        dist = tmpdist[tmpID2]
        dz = tmpdz[tmpID2]

        # Compute the local value for time stability
        CFL = numpy.zeros(1)
        CFL[0] = numpy.amin(dist / (self.erodibility * (discharge)**(self.m) * \
                ( dz / dist )**(self.n-1.) ))

        # Global mimimum value for diffusion stability
        comm.Allreduce(mpi.IN_PLACE,CFL,op=mpi.MIN)
        self.CFL = CFL[0]

    ### THIS IS NOT USED
    def distribute_deposits(self, dt, elev, GIDs, Ngbs, Acell, sedrate, sea):
        """
        Distribute sediment on base-level over the surroundings nodes in each catchment.

        Parameters
        ----------
        variable : dt
            Real value corresponding to the maximal stability time step.

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable: GIDs
            Numpy integer-type array containing for local nodes their global IDs.

        variable : Ngbs
            Numpy array containing ID of neighbor nodes.

        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        variable : sedrate
            Numpy float-type array containing sediment fluxes values (m/a).

        variable : sea
            Sea level elevation (m).
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Get basin starting IDs for each local partition
        cumbase = numpy.zeros(size+1)
        for i in range(size):
            cumbase[i+1] = len(numpy.array_split(self.base, size)[i])+cumbase[i]+1

        # Compute basin IDs using libUtils for all nodes in the local partition
        basinID = FLOWalgo.flowcompute.base_ids(self.localstack,self.receivers,cumbase[rank])
        comm.Allreduce(mpi.IN_PLACE,basinID,op=mpi.MAX)

        # Look for nodes at the edges of each basin
        tmpID = FLOWalgo.flowcompute.basin_edges(GIDs,Ngbs[GIDs],basinID)
        eIDs = numpy.where(tmpID >= 0)[0]
        edgeID = tmpID[eIDs]

        # Find baselevel nodes under deposition and above sea-level
        depID = numpy.where(sedrate[GIDs]*dt > 0.1)[0]
        tmpBase = numpy.intersect1d(self.localbase,depID)
        tmpbID = numpy.where(elev[tmpBase] > sea)[0]
        depBase = tmpBase[tmpbID]

        # Store temporary elevation array
        tmpElev = numpy.copy( elev )
        lIDs = numpy.empty(len(elev), dtype=int)

        # Define for all local enclosed basins the maximum volume of deposition
        newdt = dt
        timePD = numpy.zeros(len(depBase))
        hBasin = []
        idBasin = []
        for i in range(len(depBase)):

            timePD[i] = -1.
            # Find edges nodes IDs for the given basin
            basinEdges = numpy.where(basinID[edgeID] == basinID[depBase[i]])[0]
            edgID = numpy.where(edgeID[basinEdges] != depBase[i])[0]
            edgeList = edgeID[basinEdges][edgID]

            # Find edges node ID with minimal elevation: it is the limit before
            # overfilling of the basin
            edgeMinID = numpy.argmin(tmpElev[edgeList]) # + sedrate[edgeList] * newdt)
            overspillID = edgeList[edgeMinID]
            overspillH = tmpElev[edgeList[edgeMinID]] #+ sedrate[edgeList[edgeMinID]] * newdt

            if sedrate[depBase[i]] * newdt > overspillH - elev[depBase[i]]:
                # Find all points in the considered enclosed basin
                baseIDs = numpy.where(basinID == basinID[depBase[i]])[0]
                tmpElev[baseIDs] += sedrate[baseIDs] * newdt
                tmpElev[depBase[i]] = elev[depBase[i]]
                if sedrate[overspillID] * newdt < 0:
                    tmpElev[overspillID] = overspillH

                # Find all points below overfilling elevation
                depIDs = numpy.where(tmpElev[baseIDs] < overspillH)[0]
                if sedrate[overspillID] * newdt < 0:
                    tmpElev[overspillID] += sedrate[overspillID] * newdt
                depositIDs = baseIDs[depIDs]

                # Find edges for this set of points it will be considered as the boundary nodes
                # for pit filling
                basinNgbh = numpy.unique(Ngbs[depositIDs].flatten())
                nIDs = numpy.where(basinNgbh >= 0)[0]
                Neighbors = numpy.unique(basinNgbh[nIDs])
                boundPts = numpy.setdiff1d(Neighbors,depositIDs)

                if not (numpy.in1d(overspillID,boundPts).all()):
                    raise ValueError('Problem defining catchment based pit filling algorithm.')

                # Merge both boundary points and potential depositional points
                allID = numpy.concatenate((boundPts, depositIDs), axis=0)

                #tmp = numpy.sort(allID, axis=None)
                #if len(tmp[numpy.diff(tmp) == 0]) > 0:
                #    raise ValueError('Duplicate points found in the depositional points.')

                # Compute sediment volume to distribute from the base level point
                depVol = Acell[depBase[i]] * sedrate[depBase[i]] * dt

                # Compute volume to reach spillover
                hDep = overspillH - tmpElev[allID]
                hDep = numpy.clip(hDep,0,overspillH - tmpElev[depBase[i]])
                maxVol = numpy.sum(hDep * Acell[allID])

                if depVol > maxVol:
                    # Perform Planchon & Darboux pit filling algorithm on the given set of points
                    epsilon = 0.001
                    lIDs.fill(-1)
                    lIDs[allID] = numpy.arange(len(allID))
                    hMax = PDalgo.pdcompute.basin_filling(allID, lIDs, tmpElev,
                                                    Ngbs[allID], epsilon, len(boundPts))
                    idBasin.append(allID)
                    hBasin.append(hMax)
                    # Get the corresponding time step based on sedimentary flux
                    hDep = hMax - tmpElev[allID]
                    maxVol = numpy.sum(hDep * Acell[allID])
                    if maxVol > 0:
                        timePD[i] = Acell[depBase[i]] * sedrate[depBase[i]] / maxVol
                        newdt = min(newdt,timePD[i])
                        newdt = max(self.mindt,newdt)
                    else:
                        raise ValueError('Maximum volume of the basin is null.')
                else:
                    idBasin.append(allID)
                    hBasin.append(hDep)
            else:
                idBasin.append([-1])
                hBasin.append([0.])

        # Force the time step to ensure stability
        timestep = numpy.zeros(1)
        timestep[0] = newdt
        comm.Allreduce(mpi.IN_PLACE,timestep,op=mpi.MIN)
        newdt = timestep[0]

        # Update erosion deposition changes with the new timestep
        distH = numpy.zeros(len(elev))
        diff = sedrate * newdt

        # Get the elevation changes using stability criteria
        tmpElev = elev + diff

        # Define for all enclosed basins the deposition thicknesses
        for i in range(len(depBase)):
            if idBasin[i][0] >= 0:
                allID = idBasin[i]
                tmpElev[depBase[i]] = elev[depBase[i]]

                # Compute sediment volume to distribute from the base level point
                depVol = Acell[depBase[i]] * diff[depBase[i]]

                # Check if the requested time step is sufficient for overfilling
                if timePD[i] <= newdt:
                    # Thickness up to overfilling
                    hMax = hBasin[i]

                    # Get the deposition volume to reach overfilling
                    hDep = hMax - tmpElev[allID]
                    maxVol = numpy.sum(hDep * Acell[allID])

                    # Ensure sediment volume conservation
                    if depVol == maxVol:
                        distH[allID] = hDep
                    elif depVol > maxVol:
                        hplus = (depVol - maxVol) / numpy.sum(Acell[allID])
                        distH[allID] = hDep + hplus
                    else:
                        distH[allID] = hDep * depVol / maxVol
                else:
                    # Sort basin elevation from basin baselevel upward
                    tmpID = numpy.argsort(tmpElev[allID])
                    idSort = allID[tmpID]
                    #if idSort[0] != depBase[i]:
                    #    raise ValueError('The lowest point is not the base-level one.')

                    # Fill recursively until deposition volume is reached
                    hFill = PDalgo.pdcompute.fill_recursive(depVol, Acell[idSort], tmpElev[idSort])
                    distH[idSort] = hFill

        # Export local changes globally
        diff = (distH + tmpElev) - elev
        comm.Allreduce(mpi.IN_PLACE,diff,op=mpi.MAX)

        #diff = sedrate * newdt
        #newdt = dt

        return newdt, diff
