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
import time
import numpy
import warnings
import mpi4py.MPI as mpi
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter
from matplotlib import path

import h5py
import pandas as pd
import xml.etree.ElementTree as ETO

import os
if 'READTHEDOCS' not in os.environ:
    import pyBadlands.libUtils.sfd as sfd
    from pyBadlands.libUtils import PDalgo
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
        self.base1 = None
        self.localbase = None
        self.localbase1 = None
        self.receivers = None
        self.receivers1 = None
        self.arrdonor = None
        self.delta = None
        self.donors = None
        self.donors1 = None
        self.localstack = None
        self.localstack1 = None
        self.stack = None
        self.stack1 = None
        self.partFlow = None
        self.maxdonors = 0
        self.CFL = None
        self.erodibility = None
        self.m = None
        self.n = None
        self.mindt = None
        self.spl = False
        self.depo = 0

        self.discharge = None
        self.localsedflux = None
        self.maxh = None
        self.maxdep = None
        self.diff_flux = None
        self.diff_cfl = None
        self.chi = None
        self.basinID = None
        self.pitID = None
        self.pitVolume = None
        self.pitDrain = None
        self.allDrain = None

        self.xgrid = None
        self.xgrid = None
        self.xi = None
        self.yi = None
        self.xyi = None
        self.parentIDs = None

        self._comm = mpi.COMM_WORLD
        self._rank = self._comm.Get_rank()
        self._size = self._comm.Get_size()

    def compute_hillslope_diffusion(self, elev, neighbours, edges, distances, globalIDs):
        """
        Perform hillslope evolution based on diffusion processes.

        Parameters:
            elev: Numpy arrays containing the elevation of the TIN nodes.
            neighbours: Numpy integer-type array with the neighbourhood IDs.
            edges: Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.
            distances: Numpy real-type array with the distances between each connection in the TIN.
            globalIDs: Numpy integer-type array containing for local nodes their global IDs.
        """
        diff_flux = sfd.diffusion(elev, neighbours, edges, distances, globalIDs)

        # Send local diffusion flux globally
        self._comm.Allreduce(mpi.IN_PLACE,diff_flux,op=mpi.MAX)
        self.diff_flux = diff_flux

    def SFD_receivers(self, fillH, elev, neighbours, edges, distances, globalIDs):
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
        """

        # Call the SFD function from libUtils
        # Get the directions from true surface
        base1, receivers1 = sfd.directions_base(elev, neighbours, edges, distances, globalIDs)

        # Send local base level globally
        self._comm.Allreduce(mpi.IN_PLACE,base1,op=mpi.MAX)
        bpos = numpy.where(base1 >= 0)[0]
        self.base1 = base1[bpos]

        # Send local receivers globally
        self._comm.Allreduce(mpi.IN_PLACE,receivers1,op=mpi.MAX)
        self.receivers1 = receivers1

        # Get the directions from filled surface
        base, receivers, maxh, maxdep = sfd.directions(fillH, elev, neighbours, edges, distances, globalIDs)

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

    def _donors_number_array1(self):
        """
        Creates an array containing the number of donors for each node.
        """

        self.arrdonor = None
        numPts = len(self.receivers1)
        self.arrdonor = numpy.zeros(numPts, dtype=int)
        maxID = numpy.max(self.receivers1)
        self.arrdonor[:(maxID+1)] = numpy.bincount(self.receivers1)
        self.maxdonors = self.arrdonor.max()

        return

    def _delta_array1(self):
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

    def ordered_node_array_filled(self):
        """
        Creates an array of node IDs that is arranged in order from downstream
        to upstream for filled surface.
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

    def ordered_node_array_elev(self):
        """
        Creates an array of node IDs that is arranged in order from downstream
        to upstream for real surface.
        """

        # Build donors array for each node
        self._donors_number_array1()

        # Create the delta array
        self._delta_array1()

        # Using libUtils stack create the ordered node array
        self.donors1,lstcks = FLWnetwork.fstack.build(self.localbase1,self.receivers1,self.delta)

        # Create local stack
        stids = numpy.where(lstcks > -1 )[0]
        self.localstack1 = lstcks[stids]

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

    def compute_parameters_depression(self, fillH, elev, Acell, sealevel, debug=False):
        """
        Calculates each depression maximum deposition volume and its downstream draining node.

        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        variable : sealevel
            Real value giving the sea-level height at considered time step.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Compute pit ID and volume using libUtils
        pitID, pitVolume = FLOWalgo.flowcompute.basinparameters(self.localstack1,self.receivers1,
                                                                elev,fillH,Acell)
        comm.Allreduce(mpi.IN_PLACE,pitID,op=mpi.MAX)
        comm.Allreduce(mpi.IN_PLACE,pitVolume,op=mpi.MAX)

        self.pitID = pitID
        self.pitVolume = pitVolume

        # Find the depression node IDs
        pIDs = numpy.where(self.pitVolume>=0.)[0]

        if(len(pIDs)>0):
            xyidd = numpy.where(self.pitVolume==self.pitVolume.max())[0]
            # Order the pits based on filled elevation from top to bottom
            orderPits = numpy.argsort(fillH[pIDs])[::-1]

            # Find the depression or edge, marine point where a given pit is draining
            self.pitDrain = FLOWalgo.flowcompute.basindrainage(orderPits,self.pitID,self.receivers,pIDs,
                                                          fillH,sealevel)
            self.allDrain = FLOWalgo.flowcompute.basindrainageall(orderPits,self.pitID,self.receivers,pIDs)

            # Debugging plotting function
            #debug = True
            if debug:
                self.visualise_draining_path(pIDs, elev, self.pitDrain, fillH, 'drain')
                self.visualise_draining_path(pIDs, elev, self.allDrain, fillH, 'alldrain')
        else:
            self.pitDrain = -numpy.ones(len(pitID))
            self.allDrain = -numpy.ones(len(pitID))

        return

    def compute_sedflux(self, Acell, elev, fillH, xymin, xymax, dt, rivqs, sealevel,
        cumdiff, perc_dep, slp_cr, sigma, verbose=False):
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

        variable : dt
            Real value corresponding to the maximal stability time step.

        variable : rivqs
            Numpy arrays representing the sediment fluxes from rivers.

        variable : sealevel
            Real value giving the sea-level height at considered time step.

        variable : cumdiff
            Numpy array containing the cumulative deposit thicknesses.

        variable : slp_cr
            Critical slope used to force aerial deposition for alluvial plain.

        variable : perc_dep
            Maximum percentage of deposition at any given time interval.

        variable : sigma
            Marine sedimentation gaussian filter parameter.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        newdt = numpy.copy(dt)
        sedflux = numpy.zeros(len(elev))
        verbose = False

        # Compute sediment flux using libUtils
        # Stream power law
        if self.spl:
            if rank==0 and verbose:
                time0 = time.clock()
                time1 = time.clock()

            # Find border/inside nodes
            ids = numpy.arange(len(Acell))
            tmp1 = numpy.where(Acell>0.)[0]
            domain = path.Path([(xymin[0],xymin[1]),(xymax[0],xymin[1]), (xymax[0],xymax[1]), (xymin[0],xymax[1])])
            tmp2 = domain.contains_points(self.xycoords)
            insideIDs = numpy.intersect1d(tmp1,ids[tmp2])
            borders = numpy.zeros(len(Acell),dtype=int)
            borders[insideIDs] = 1

            cdepo, cero = FLOWalgo.flowcompute.streampower(self.localstack,self.receivers,self.pitID,self.pitVolume, \
                     self.pitDrain,self.xycoords,Acell,self.maxh,self.maxdep,self.discharge,fillH,elev,rivqs, \
                     self.erodibility,self.m,self.n,perc_dep,slp_cr,sealevel,newdt,borders)
            comm.Allreduce(mpi.IN_PLACE,cdepo,op=mpi.MAX)
            comm.Allreduce(mpi.IN_PLACE,cero,op=mpi.MIN)

            if self.depo == 0:
                volChange = cero
            else:
                volChange = cdepo+cero
            if rank==0 and verbose:
                print "   - Compute sediment volumetric flux ", time.clock() - time1
                time1 = time.clock()

            # Find overfilling catchments
            ids = numpy.where(numpy.logical_and(volChange>self.pitVolume,self.pitVolume>0.))[0]

            # Check if there are some internally drained depressions within the computational domain?
            intID = numpy.where(numpy.logical_and(self.allDrain == self.pitID,self.pitID>=0))[0]
            if len(intID)>0 and len(ids)>0:
                search = domain.contains_points(self.xycoords[intID])
                # For all these closed basins find the ones overfilled
                if len(search) > 0:
                    overfilled = numpy.intersect1d(intID[search],ids)
                    # Limit the time step to restrict deposition in these basins
                    if len(overfilled) > 0:
                        # Compute the percentage of overfilling
                        percOver = self.pitVolume[overfilled]/volChange[overfilled]
                        newdt = dt*percOver.min()

            if newdt>1.:
                newdt = float(round(newdt-0.5,0))
            newdt = max(self.mindt,newdt)

            if rank==0 and verbose:
                print "   - Compute depressions connectivity ", time.clock() - time1
                time1 = time.clock()
            if newdt < dt:
                cdepo, cero = FLOWalgo.flowcompute.streampower(self.localstack,self.receivers,self.pitID,self.pitVolume, \
                    self.pitDrain,self.xycoords,Acell,self.maxh,self.maxdep,self.discharge,fillH,elev,rivqs, \
                    self.erodibility,self.m,self.n,perc_dep,slp_cr,sealevel,newdt,borders)
                comm.Allreduce(mpi.IN_PLACE,cdepo,op=mpi.MAX)
                comm.Allreduce(mpi.IN_PLACE,cero,op=mpi.MIN)
                volChange = cdepo+cero
                if rank==0 and verbose:
                    print "   - Compute volumetric fluxes with updated dt ", time.clock() - time1
                    time1 = time.clock()

                # Ensure no overfilling remains
                ids = numpy.where(numpy.logical_and(volChange>self.pitVolume,self.pitVolume>0.))[0]
                search = domain.contains_points(self.xycoords[intID])
                if (len(search)>0) and (len(ids)>0):
                    overfilled = numpy.intersect1d(intID[search],ids)
                    if len(overfilled) > 0:
                        print 'WARNING: overfilling persists after time-step limitation.',len(overfilled)
                    #assert len(overfilled) == 0, 'WARNING: overfilling persists after time-step limitation.'

            # Compute erosion
            ero = numpy.zeros(len(cero))
            ero[insideIDs] = cero[insideIDs]
            erosion = numpy.zeros(len(ero))
            erosion[insideIDs] = ero[insideIDs]/Acell[insideIDs]
            if rank==0 and verbose:
                print "   - Compute erosion ", time.clock() - time1
                time1 = time.clock()

            # Compute deposition
            if self.depo == 0:
                # Purely erosive case
                deposition = numpy.zeros(len(cdepo))
            else:
                depo = numpy.zeros(len(cdepo))
                depo[insideIDs] = cdepo[insideIDs]
                deposition = numpy.zeros(len(depo))

                tmp = numpy.where(elev>sealevel)[0]
                landIDs = numpy.intersect1d(tmp,insideIDs)

                # Compute alluvial plain deposition
                depID = numpy.where(numpy.logical_and(fillH==elev,depo>0.))[0]
                plainID = numpy.intersect1d(depID,landIDs)
                if len(plainID) > 0:
                    deposition[plainID] = depo[plainID]/Acell[plainID]
                    depo[plainID] = 0.
                    if rank==0 and verbose:
                        print "   - Compute plain deposition ", time.clock() - time1
                        time1 = time.clock()

                # Compute land pit deposition
                pitIDs = numpy.where(numpy.logical_and(elev>sealevel,fillH>sealevel))[0]
                volIDs = numpy.where(self.pitVolume>0.)[0]
                tmpIDs = numpy.intersect1d(volIDs,pitIDs)
                depID = numpy.where(depo>0.)[0]
                landIDs = numpy.intersect1d(depID,tmpIDs)
                if len(landIDs) > 0:
                    perc = numpy.zeros(len(depo))
                    # Get the percentage to deposit
                    perc[landIDs] = depo[landIDs]/self.pitVolume[landIDs]
                    tmp = numpy.where(perc>1)[0]
                    overfilled = numpy.intersect1d(tmp,insideIDs)
                    if len(overfilled) > 0:
                        print 'WARNING: overfilling persists during land pit deposition.',len(overfilled)
                    #assert len(overfilled) == 0, 'WARNING: overfilling persists during land pit deposition.'

                    for p in range(len(landIDs)):
                        tmp = numpy.where(self.pitID==landIDs[p])[0]
                        perc[tmp] = perc[landIDs[p]]
                    perc[perc>1.] = 1.
                    tmp = numpy.where(perc>0.)
                    deposition[tmp] = (fillH[tmp]-elev[tmp])*perc[tmp]
                    depo[landIDs] = 0.
                    if rank==0 and verbose:
                        print "   - Compute land pit deposition ", time.clock() - time1
                        time1 = time.clock()

                # Compute water deposition
                tmp = numpy.where(numpy.logical_and(depo>0.,elev<=sealevel))[0]
                seaIDs = numpy.intersect1d(tmp,insideIDs)
                if len(seaIDs) > 0:
                    seavol = numpy.zeros(len(depo))
                    seavol[seaIDs] = depo[seaIDs]
                    # Distribute marine sediments based on angle of repose
                    seadep = PDalgo.pdstack.marine_sed(elev, seavol, borders, sealevel)
                    if rank==0 and verbose:
                        print "   - Compute marine deposition ", time.clock() - time1
                        time1 = time.clock()
                    if sigma > 0.:
                        smthdep = self.gaussian_diffusion(seadep,sigma)
                        frac = numpy.sum(seadep*Acell)/numpy.sum(smthdep*Acell)
                        deposition += smthdep*frac
                    else:
                        deposition += seadep
                    depo[seaIDs] = 0.
                    if rank==0 and verbose:
                        print "   - Smooth marine deposition ", time.clock() - time1
                        time1 = time.clock()

                # Is there some remaining deposits?
                tmp = numpy.where(depo>0)[0]
                remainIDs = numpy.intersect1d(tmp,insideIDs)
                if len(remainIDs)>0:
                    if rank == 0:
                        print 'WARNING: forced deposition is performed during this timestep.',len(remainIDs)
                    deposition[remainIDs] = depo[remainIDs]/Acell[remainIDs]

            # Define erosion/deposition changes
            sedflux[insideIDs] = erosion[insideIDs]+deposition[insideIDs]
            if rank==0 and verbose:
                print "   - Total sediment flux time ", time.clock() - time0

        return newdt,sedflux

    def gaussian_diffusion(self, diff, dsmooth):
        """
        Gaussian filter operation used to smooth diffusion related deposition
        thicknesses.

        Parameters
        ----------
        variable : diff
            Numpy arrays containing the deposition thicknesses.

        variable : dsmooth
            Smoothing parameter.
        """

        if self.xgrid is None:
            dx = self.xycoords[1,0] - self.xycoords[0,0]
            xmin, xmax = min(self.xycoords[:,0]), max(self.xycoords[:,0])
            ymin, ymax = min(self.xycoords[:,1]), max(self.xycoords[:,1])
            self.xgrid = numpy.arange(xmin,xmax+dx,dx)
            self.ygrid = numpy.arange(ymin,ymax+dx,dx)
            self.xi, self.yi = numpy.meshgrid(self.xgrid, self.ygrid)

            # Querying the cKDTree later becomes a bottleneck, so distribute the xyi array across all MPI nodes
            xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]
            splits = numpy.array_split(xyi, self._size)
            self.split_lengths = numpy.array(map(len, splits)) * 3
            self.localxyi = splits[self._rank]
            self.query_shape = (xyi.shape[0], 3)

        depZ = numpy.copy(diff)
        tree = cKDTree(self.xycoords[:,:2])

        # Querying the KDTree is rather slow, so we split it across MPI nodes
        nelems = self.query_shape[0] * self.query_shape[1]
        indices = numpy.empty(self.query_shape, dtype=numpy.int64)
        localdistances, localindices = tree.query(self.localxyi, k=3)

        distances_flat = numpy.empty(nelems, dtype=numpy.float64)
        self._comm.Allgatherv(numpy.ravel(localdistances), [distances_flat, (self.split_lengths, None)])

        indices_flat = numpy.empty(nelems, dtype=numpy.int64)
        self._comm.Allgatherv(numpy.ravel(localindices), [indices_flat, (self.split_lengths, None)])

        distances = distances_flat.reshape(self.query_shape)
        indices = indices_flat.reshape(self.query_shape)

        if len(depZ[indices].shape) == 3:
            zd_vals = depZ[indices][:,:,0]
        else:
            zd_vals = depZ[indices]
        with numpy.errstate(divide='ignore'):
            zdi = numpy.average(zd_vals,weights=(1./distances), axis=1)

        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            zdi[onIDs] = depZ[indices[onIDs,0]]

        depzi = numpy.reshape(zdi,(len(self.ygrid),len(self.xgrid)))

        smthDep = gaussian_filter(depzi, sigma=dsmooth)

        rgi_dep = RegularGridInterpolator((self.ygrid, self.xgrid), smthDep)
        zdepsmth = rgi_dep((self.xycoords[:,1],self.xycoords[:,0]))

        return zdepsmth

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

        return

    def visualise_draining_path(self, pIDs, elev, drain, fillH, filename):
        """
        Debugging function used to plot draining pathway between depressions.

        Parameters
        ----------
        variable : pIDs
            Numpy array containing all depression node IDs.

        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.

        variable : drain
            Numpy arrays containing the draining pit ID.

        variable : filename
            Name of the output file.
        """

        n1 = len(pIDs)
        coords = numpy.zeros((n1*2,3))
        coords[:n1,0] = self.xycoords[pIDs,0]
        coords[:n1,1] = self.xycoords[pIDs,1]
        coords[:n1,2] = elev[pIDs]
        coords[n1:,0] = self.xycoords[drain[pIDs],0]
        coords[n1:,1] = self.xycoords[drain[pIDs],1]
        coords[n1:,2] = elev[drain[pIDs]]
        connect = numpy.zeros((n1,2), dtype=int)
        connect[:,0] = numpy.arange(n1)
        connect[:,1] = numpy.arange(n1)+n1
        h5file = filename+'.hdf5'
        with h5py.File(h5file, "w") as f:
            # Write node coordinates and elevation
            f.create_dataset('coords',shape=(n1*2,3), dtype='float32', compression='gzip')
            f["coords"][:,:3] = coords

            f.create_dataset('connect',shape=(n1,2), dtype='int32', compression='gzip')
            f["connect"][:,:2] = connect

        xmf_file = filename+'.xmf'
        f= open(str(xmf_file),'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(' <Domain>\n')
        f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
        f.write('      <Time Type="Single" Value="0"/>\n')
        f.write('      <Grid Name="Block.0">\n')
        f.write('         <Topology Type="Polyline" NodesPerElement="2" ')
        f.write('NumberOfElements="%d" BaseOffset="0">\n'%n1)
        f.write('          <DataItem Format="HDF" DataType="Int" ')
        f.write('Dimensions="%d 2">%s:/connect</DataItem>\n'%(n1,h5file))
        f.write('         </Topology>\n')
        f.write('         <Geometry Type="XYZ">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 3">%s:/coords</DataItem>\n'%(n1*2,h5file))
        f.write('         </Geometry>\n')
        f.write('      </Grid>\n')
        f.write('    </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()

        df = pd.DataFrame({'X':self.xycoords[pIDs,0],'Y':self.xycoords[pIDs,1],'Z':elev[pIDs],'V':self.pitVolume[pIDs],
                             'ID':self.pitID[pIDs],'Drain':self.pitDrain[pIDs]})
        df.to_csv(filename+'vol.csv',columns=['X', 'Y', 'Z', 'V', 'ID','Drain'], sep=',', index=False)

        # df = pd.DataFrame({'X':self.xycoords[:,0],'Y':self.xycoords[:,1],'Z':elev[:],'W':fillH[:]})
        # df.to_csv(filename+'water.csv',columns=['X', 'Y', 'Z', 'W'], sep=',', index=False)

        return
