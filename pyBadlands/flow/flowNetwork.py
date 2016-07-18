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
import pyBadlands.libUtils.sfd as sfd
from pyBadlands.libUtils import FLOWalgo
from pyBadlands.libUtils import FLWnetwork

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
        self.diff_cfl = None
        self.chi = None
        self.basinID = None

        self.xgrid = None
        self.xgrid = None
        self.xi = None
        self.yi = None
        self.xyi = None

        self._comm = mpi.COMM_WORLD
        self._rank = self._comm.Get_rank()
        self._size = self._comm.Get_size()

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
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        variable : distances
            Numpy real-type array with the distances between each connection in the TIN.

        variable: globalIDs
            Numpy integer-type array containing for local nodes their global IDs.

        variable : sea
            Current elevation of sea level.
        """

        # Call the SFD function from libUtils
        if self.depo == 0 or self.capacity or self.filter:
            base, receivers, diff_flux = sfd.directions_base(elev, neighbours, edges, distances, globalIDs, sea)

            # Send local base level globally
            self._comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)

            bpos = numpy.where(base >= 0)[0]
            self.base = base[bpos]
            numpy.random.shuffle(self.base)
            # Send local receivers globally
            self._comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local diffusion flux globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux
        else:
            base, receivers, maxh, maxdep, diff_flux = sfd.directions(fillH, elev, neighbours, edges, distances, globalIDs, sea)

            # Send local base level globally
            self._comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)
            bpos = numpy.where(base >= 0)[0]
            self.base = base[bpos]
            numpy.random.shuffle(self.base)

            # Send local receivers globally
            self._comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local maximum height globally
            self._comm.Allreduce(mpi.IN_PLACE,maxh,op=mpi.MAX)
            self.maxh = maxh

            # Send local maximum deposition globally
            self._comm.Allreduce(mpi.IN_PLACE,maxdep,op=mpi.MAX)
            self.maxdep = maxdep

            # Send local diffusion flux globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux

    def SFD_nreceivers(self, Sc, fillH, elev, neighbours, edges, distances, globalIDs, sea):
        """
        Single Flow Direction function computes downslope flow directions by inspecting the neighborhood
        elevations around each node. The SFD method assigns a unique flow direction towards the steepest
        downslope neighbor. In addition it compute the hillslope non-linear diffusion

        Parameters
        ----------
        variable : Sc
            Critical slope for non-linear diffusion.

        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable : neighbours
            Numpy integer-type array with the neighbourhood IDs.

        variable : edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        variable : distances
            Numpy real-type array with the distances between each connection in the TIN.

        variable: globalIDs
            Numpy integer-type array containing for local nodes their global IDs.

        variable : sea
            Current elevation of sea level.
        """

        # Call the SFD function from libUtils
        if self.depo == 0 or self.capacity or self.filter:
            base, receivers, diff_flux, diff_cfl = SFD.sfdcompute.directions_base_nl(elev, \
                neighbours, edges, distances, globalIDs, sea, Sc)

            # Send local base level globally
            self._comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)

            bpos = numpy.where(base >= 0)[0]
            self.base = base[bpos]
            numpy.random.shuffle(self.base)
            # Send local receivers globally
            self._comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local diffusion flux globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux

            # Send local diffusion CFL condition globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_cfl,op=mpi.MIN)
            self.diff_cfl = diff_cfl
        else:
            base, receivers, maxh, maxdep, diff_flux, diff_cfl = SFD.sfdcompute.directions_nl(fillH, \
                elev, neighbours, edges, distances, globalIDs, sea, Sc)

            # Send local base level globally
            self._comm.Allreduce(mpi.IN_PLACE,base,op=mpi.MAX)
            bpos = numpy.where(base >= 0)[0]
            self.base = base[bpos]
            numpy.random.shuffle(self.base)

            # Send local receivers globally
            self._comm.Allreduce(mpi.IN_PLACE,receivers,op=mpi.MAX)
            self.receivers = receivers

            # Send local maximum height globally
            self._comm.Allreduce(mpi.IN_PLACE,maxh,op=mpi.MAX)
            self.maxh = maxh

            # Send local maximum deposition globally
            self._comm.Allreduce(mpi.IN_PLACE,maxdep,op=mpi.MAX)
            self.maxdep = maxdep

            # Send local diffusion flux globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
            self.diff_flux = diff_flux

            # Send local diffusion CFL condition globally
            self._comm.Allreduce(mpi.IN_PLACE,diff_cfl,op=mpi.MIN)
            self.diff_cfl = diff_cfl

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

    def compute_flow(self, Acell, rain):
        """
        Calculates the drainage area and water discharge at each node.

        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        variable : rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).
        """

        numPts = len(Acell)

        self.discharge = numpy.zeros(numPts, dtype=float)
        self.discharge[self.stack] = Acell[self.stack] * rain[self.stack]

        # Compute discharge using libUtils
        self.discharge = FLOWalgo.flowcompute.discharge(self.localstack, self.receivers, self.discharge)
        self._comm.Allreduce(mpi.IN_PLACE, self.discharge, op=mpi.MAX)

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

        # Parallel case
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

        # Serial case
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

        K = 3

        if self.xgrid is None:
            dx = self.xycoords[1,0] - self.xycoords[0,0]
            xmin, xmax = min(self.xycoords[:,0]), max(self.xycoords[:,0])
            ymin, ymax = min(self.xycoords[:,1]), max(self.xycoords[:,1])
            self.xgrid = numpy.arange(xmin,xmax+dx,dx)
            self.ygrid = numpy.arange(ymin,ymax+dx,dx)
            self.xi, self.yi = numpy.meshgrid(self.xgrid, self.ygrid)

            # querying the cKDTree later becomes a bottleneck, so distribute the xyi array across all MPI nodes
            xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]
            splits = numpy.array_split(xyi, self._size)
            self.split_lengths = numpy.array(map(len, splits)) * K
            self.localxyi = splits[self._rank]
            self.query_shape = (xyi.shape[0], K)

        depZ = numpy.copy(diff)
        depZ = depZ.clip(0.)

        eroZ = numpy.copy(diff)
        eroZ = eroZ.clip(max=0.)

        tree = cKDTree(self.xycoords[:,:2])

        # Querying the KDTree is rather slow, so we split it across MPI nodes
        # FIXME: the Allgatherv fails if we don't flatten the array first - why?
        nelems = self.query_shape[0] * self.query_shape[1]
        indices = numpy.empty(self.query_shape, dtype=numpy.int64)
        localdistances, localindices = tree.query(self.localxyi, k=K)

        distances_flat = numpy.empty(nelems, dtype=numpy.float64)
        self._comm.Allgatherv(numpy.ravel(localdistances), [distances_flat, (self.split_lengths, None)])

        indices_flat = numpy.empty(nelems, dtype=numpy.int64)
        self._comm.Allgatherv(numpy.ravel(localindices), [indices_flat, (self.split_lengths, None)])

        distances = distances_flat.reshape(self.query_shape)
        indices = indices_flat.reshape(self.query_shape)

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

    def dt_stability(self, elev, locIDs):
        """
        This function computes the maximal timestep to ensure computation stability
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
