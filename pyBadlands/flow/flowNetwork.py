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

    def __init__(self, input):
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
        self.mindt = None
        self.spl = False
        self.depo = 0

        self.discharge = None
        self.localsedflux = None
        self.maxh = None
        self.maxdep = None
        self.diff_cfl = None
        self.chi = None
        self.basinID = None
        self.pitID = None
        self.pitVolume = None
        self.pitDrain = None
        self.allDrain = None

        self.xgrid = None
        self.ygrid = None
        self.xi = None
        self.yi = None
        self.xyi = None
        self.Afactor = None
        self.parentIDs = None
        self.dx = None
        self.distances = None
        self.indices = None
        self.onIDs = None

        self.mp = input.mp
        self.mt = input.mt
        self.nt = input.nt
        self.kt = input.kt
        self.kw = input.kw
        self.b = input.b
        self.deepb = input.deepbasin
        self.critdens = input.denscrit
        self.flowdensity = None
        self.sedload = None
        self.straTIN = 0
        self.activelay = None

        self.borders = None
        self.domain = None
        self.insideIDs = None
        self.outsideIDs = None
        self.borders2 = None
        self.insideIDs2 = None
        self.outsideIDs2 = None

        FLOWalgo.flowcompute.eroparams(input.incisiontype,input.SPLm,input.SPLn,input.mt,
                                       input.nt,input.kt,input.kw,input.b,input.bedslptype)

    def compute_hillslope_diffusion(self, elev, neighbours, edges, distances, globalIDs, type, Sc):
        """
        Perform hillslope evolution based on diffusion processes.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.

        type
            Flag to compute the diffusion when multiple rocks are used.
        """

        if type == 0:
            if Sc > 0.:
                tSc = numpy.zeros(1)
                tSc[0] = Sc
                diff_flux = sfd.diffusionnl(tSc, elev, self.borders2, neighbours, edges, distances, globalIDs)
            else:
                diff_flux = sfd.diffusion(elev, self.borders2, neighbours, edges, distances, globalIDs)
        else:
            diff_flux = sfd.diffusionero(elev, self.borders2, neighbours, edges, distances, globalIDs)

        return diff_flux

    def compute_marine_diffusion(self, elev, depoH, neighbours, edges, distances, coeff,
                                 globalIDs, seal, maxth, tstep):
        """
        Perform river transported marine sediments diffusion.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        dep
            Numpy arrays flagging the deposited nodes.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        diff_flux, ndt = FLOWalgo.flowcompute.diffmarine(elev, self.borders, depoH, neighbours, edges,
                                           distances, coeff, globalIDs, seal, maxth, tstep)

        # Send local diffusion flux globally
        mindt = numpy.array(ndt)

        return diff_flux,mindt

    def compute_failure_diffusion(self, elev, depoH, neighbours, edges, distances, coeff,
                                 globalIDs, maxth, tstep):
        """
        Perform slope failure transported sediments diffusion.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        dep
            Numpy arrays flagging the deposited nodes.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        diff_flux, ndt = FLOWalgo.flowcompute.difffailure(elev, self.borders, depoH, neighbours, edges,
                                           distances, coeff, globalIDs, maxth, tstep)

        # Send local diffusion flux globally
        mindt = numpy.array(ndt)

        return diff_flux,mindt

    def compute_failure(self, elev, sfail):
        """
        Perform erosion induced by slope failure.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        sfail
            Critical slope to initiate slope failure.
        """

        erosion = FLOWalgo.flowcompute.slumpero(self.localstack,self.receivers,self.xycoords, \
                             elev,sfail,self.borders)

        return erosion

    def compute_sediment_marine(self, elev, dep, sdep, coeff, neighbours, seal, maxth,
                                edges, distances, globalIDs):
        """
        Perform sediment diffusion for multiple rock types.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        dep
            Numpy arrays containing the rock deposition.

        coeff
            Numpy arrays containing the coefficient value for the diffusion algorithm.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        diff_prop, diff_flux = FLOWalgo.flowcompute.diffsedmarine(elev, self.borders, dep, sdep,
                                           seal, maxth, coeff, neighbours, edges, distances, globalIDs)

        return diff_prop,diff_flux

    def compute_sediment_hillslope(self, elev, difflay, coeff, neighbours,
                                    edges, layh, distances, globalIDs):
        """
        Perform sediment diffusion for multiple rock types.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        difflay
            Numpy arrays containing the rock type fractions in the active layer.

        coeff
            Numpy arrays containing the coefficient value for the diffusion algorithm.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        layh
            Numpy arrays containing the thickness of the active layer.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        sumdiff, ero, depo = FLOWalgo.flowcompute.diffsedhillslope(elev, self.borders, difflay,
                                           layh, coeff, neighbours, edges, distances, globalIDs)

        return sumdiff,ero,depo

    def SFD_receivers(self, fillH, elev, neighbours, edges, distances, globalIDs):
        """
        Single Flow Direction function computes downslope flow directions by inspecting the neighborhood
        elevations around each node. The SFD method assigns a unique flow direction towards the steepest
        downslope neighbor.

        Parameters
        ----------
        fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        elev
            Numpy arrays containing the elevation of the TIN nodes.

        neighbours
            Numpy integer-type array with the neighbourhood IDs.

        edges
            Numpy real-type array with the voronoi edges length for each neighbours of the TIN nodes.

        distances
            Numpy real-type array with the distances between each connection in the TIN.

        globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """

        # Call the SFD function from libUtils
        # Get the directions from true surface
        base1, receivers1 = sfd.directions_base(elev, neighbours, edges, distances, globalIDs)

        # Send local base level globally
        bpos = numpy.where(base1 >= 0)[0]
        self.base1 = base1[bpos]

        # Send local receivers globally
        self.receivers1 = receivers1

        # Get the directions from filled surface
        base, receivers, maxh, maxdep = sfd.directions(fillH, elev, neighbours, edges, distances, globalIDs)

        # Send local base level globally
        bpos = numpy.where(base >= 0)[0]
        self.base = base[bpos]
        numpy.random.shuffle(self.base)

        # Send local receivers globally
        self.receivers = receivers

        # Send local maximum height globally
        self.maxh = maxh

        # Send local maximum deposition globally
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

    def compute_flow(self, sealevel, elev, Acell, rain):
        """
        Calculates the drainage area and water discharge at each node.

        Parameters
        ----------
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        Acell
            Numpy float-type array containing the voronoi area for each nodes (in m^2)

        rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).
        """

        numPts = len(Acell)

        self.discharge = numpy.zeros(numPts, dtype=float)
        self.discharge[self.stack] = Acell[self.stack] * rain[self.stack]

        # Compute discharge using libUtils
        self.discharge, self.activelay = FLOWalgo.flowcompute.discharge(sealevel, self.localstack, self.receivers,
                                                        elev, self.discharge)

    def view_receivers(self, fillH, elev, neighbours, edges, distances, globalIDs, sea):
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
        base, receivers = sfd.dirview(fillH, elev, neighbours, edges, distances, globalIDs, sea)

        # Send local base level globally
        bpos = numpy.where(base >= 0)[0]
        self.base = base[bpos]
        numpy.random.shuffle(self.base)

        # Send local receivers globally
        self.receivers = receivers
        self.localbase = self.base
        self.ordered_node_array_filled()

        globalstack = self.localstack
        self.stack = globalstack

        return

    def compute_parameters(self):
        """
        Calculates the catchment IDs and the Chi parameter (Willett 2014).
        """

        # Get basin starting IDs for each local partition
        cumbase = numpy.zeros(2)
        for i in range(1):
            cumbase[i+1] = len(numpy.array_split(self.base, 1)[i])+cumbase[i]+1

        # Compute discharge using libUtils
        chi, basinID = FLOWalgo.flowcompute.parameters(self.localstack,self.receivers,
                                               self.discharge,self.xycoords,cumbase[0])

        self.chi = chi
        self.basinID = basinID

    def compute_parameters_depression(self, fillH, elev, Acell, sealevel, debug=False):
        """
        Calculates each depression maximum deposition volume and its downstream draining node.

        Parameters
        ----------

        fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.

        elev
            Numpy arrays containing the elevation of the TIN nodes.

        Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)

        sealevel
            Real value giving the sea-level height at considered time step.
        """

        # Compute pit ID and volume using libUtils
        pitID, pitVolume = FLOWalgo.flowcompute.basinparameters(self.localstack1,self.receivers1,
                                                                elev,fillH,Acell)
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

    def compute_sedflux(self, Acell, elev, rain, fillH, dt, actlay, rockCk, rivqs,
        sealevel, perc_dep, slp_cr, ngbh, verbose=False):
        """
        Calculates the sediment flux at each node.

        Parameters
        ----------
        Acell
            Numpy float-type array containing the voronoi area for each nodes (in m^2)

        elev
            Numpy arrays containing the elevation of the TIN nodes.

        rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).

        fillH
            Numpy array containing the lake elevations.

        dt
            Real value corresponding to the maximal stability time step.

        actlay
            Active layer composition.

        rockCk
            Rock erodibility values.

        rivqs
            Numpy arrays representing the sediment fluxes from rivers.

        sealevel
            Real value giving the sea-level height at considered time step.

        slp_cr
            Critical slope used to force aerial deposition for alluvial plain.

        perc_dep
            Maximum percentage of deposition at any given time interval.
        """

        check = False

        newdt = numpy.copy(dt)

        if actlay is None:
            sedflux = numpy.zeros((len(elev),1))
        else:
            sedflux = numpy.zeros((len(elev),len(rockCk)))

        # Compute sediment flux using libUtils
        # Stream power law
        if self.spl:
            if verbose:
                time0 = time.clock()
                time1 = time.clock()

            # Find border/inside nodes
            if self.mp>0.:
                if self.straTIN == 1:
                    rp = numpy.power(rain,self.mp).reshape((len(elev),1))
                    eroCoeff = rockCk * rp
                else:
                    eroCoeff = self.erodibility*numpy.power(rain,self.mp)
                    eroCoeff.reshape((len(elev),1))
            else:
                if self.straTIN == 1:
                    eroCoeff = numpy.tile(rockCk,(len(elev),1))
                else:
                    eroCoeff = self.erodibility.reshape((len(elev),1))
            if actlay is None:
                actlay = numpy.zeros((len(elev),1))

            cdepo, cero, sedload, flowdensity = FLOWalgo.flowcompute.streampower(self.critdens, self.localstack,self.receivers,self.pitID, \
                     self.pitVolume,self.pitDrain,self.xycoords,Acell,self.maxh,self.maxdep,self.discharge,fillH, \
                     elev,rivqs,eroCoeff,actlay,perc_dep,slp_cr,sealevel,sealevel+self.deepb,newdt,self.borders)

            if self.depo == 0:
                volChange = cero
            else:
                volChange = cdepo+cero
            if verbose:
                print("   - Compute sediment volumetric flux ", time.clock() - time1)
                time1 = time.clock()

            # Find overfilling catchments
            tmpChange,id1,id2,nb1,nb2 = FLOWalgo.flowcompute.getid1(volChange,self.pitVolume,self.allDrain,self.pitID)

            # Check if there are some internally drained depressions within the computational domain?
            if nb1>0 and nb2>0:
                ids = id1[:nb1]
                intID = id2[:nb2]
                search = self.domain.contains_points(self.xycoords[intID])
                # For all these closed basins find the ones overfilled
                if len(search) > 0:
                    overfilled = numpy.intersect1d(intID[search],ids)
                    # Limit the time step to restrict deposition in these basins
                    if len(overfilled) > 0:
                        # Compute the percentage of overfilling
                        percOver = self.pitVolume[overfilled]/tmpChange[overfilled]
                        newdt = dt*percOver.min()

            if newdt>1.:
                newdt = float(round(newdt-0.5,0))
            newdt = max(self.mindt,newdt)

            if verbose:
                print("   - Compute depressions connectivity ", time.clock() - time1)
                time1 = time.clock()

            if newdt < dt:
                cdepo, cero, sedload, flowdensity = FLOWalgo.flowcompute.streampower(self.critdens, self.localstack,self.receivers,self.pitID, \
                        self.pitVolume,self.pitDrain,self.xycoords,Acell,self.maxh,self.maxdep,self.discharge,fillH, \
                        elev,rivqs,eroCoeff,actlay,perc_dep,slp_cr,sealevel,sealevel+self.deepb,newdt,self.borders)
                volChange = cdepo+cero
                if verbose:
                    print("   - Compute volumetric fluxes with updated dt ", time.clock() - time1)
                    time1 = time.clock()

                if check:
                    # Ensure no overfilling remains
                    tmpChange = numpy.sum(volChange,axis=1)
                    ids = numpy.where(numpy.logical_and(tmpChange>self.pitVolume,self.pitVolume>0.))[0]
                    search = self.domain.contains_points(self.xycoords[intID])
                    if (len(search)>0) and (len(ids)>0):
                        overfilled = numpy.intersect1d(intID[search],ids)
                        if len(overfilled) > 0:
                            print('WARNING: overfilling persists after time-step limitation.',len(overfilled))

            # Update river sediment load in kg/s
            sedld = numpy.sum(sedload,axis=1)
            self.sedload = sedld/(newdt*3.154e7)
            den = flowdensity/1000.
            self.flowdensity = den

            # Compute erosion
            erosion = numpy.zeros(cero.shape)
            erosion[self.insideIDs,:] = cero[self.insideIDs,:]/Acell[self.insideIDs].reshape(len(self.insideIDs),1)
            if verbose:
                print("   - Compute erosion ", time.clock() - time1)
                time1 = time.clock()

            # Compute deposition
            if self.depo == 0:
                # Purely erosive case
                deposition = numpy.zeros(cdepo.shape)
            else:
                depo = numpy.zeros(cdepo.shape)
                depo[self.insideIDs,:] = cdepo[self.insideIDs,:]
                deposition = numpy.zeros(depo.shape)

                # Compute alluvial plain deposition
                plainid,landid,seaid,perc,nplain,nland,nsea,ndepo = FLOWalgo.flowcompute.getids(fillH,elev,depo,
                                                                    self.pitVolume,sealevel)
                depo = ndepo

                if nplain > 0:
                    plainID = plainid[:nplain]
                    deposition[plainID,:] = depo[plainID,:]/Acell[plainID].reshape(len(plainID),1)
                    depo[plainID,:] = 0.
                    if verbose:
                        print("   - Compute plain deposition ", time.clock() - time1)
                        time1 = time.clock()

                # Compute land pit deposition
                if nland > 0:
                    landIDs = landid[:nland]
                    for p in range(len(landIDs)):
                        tmp = numpy.where(self.pitID==landIDs[p])[0]

                        if len(tmp) == 1:
                            deposition[tmp,:] = (fillH[tmp]-elev[tmp])*perc[landIDs[p],:]
                        else:
                            deposition[tmp,:] = (fillH[tmp]-elev[tmp]).reshape(len(tmp),1)*perc[landIDs[p],:]
                        tmpd = numpy.sum(deposition[tmp,:]*Acell[tmp].reshape(len(Acell[tmp]),1))
                        dfrac = numpy.sum(depo[landIDs[p],:])/tmpd
                        deposition[tmp,:] *= dfrac

                    depo[landIDs,:] = 0.
                    if verbose:
                        print("   - Compute land pit deposition ", time.clock() - time1)
                        time1 = time.clock()

                # Compute water deposition
                if nsea > 0:
                    # Distribute marine sediments based on angle of repose
                    seaIDs = seaid[:nsea]
                    seavol = numpy.zeros(depo.shape)
                    seavol[seaIDs,:] = depo[seaIDs,:]
                    seadep = PDalgo.pdstack.marine_distribution(elev, seavol, sealevel, self.borders, seaIDs)
                    deposition += seadep
                    depo[seaIDs,:] = 0.
                    if verbose:
                        print("   - Compute marine deposition ", time.clock() - time1)
                        time1 = time.clock()

                # Is there some remaining deposits?
                if numpy.any(depo):
                    print('WARNING: forced deposition is performed during this timestep.')
                    deposition[self.insideIDs,:] += depo[self.insideIDs,:]/Acell[self.insideIDs].reshape(len(self.insideIDs),1)

            # Define erosion/deposition changes
            sedflux[self.insideIDs,:] = erosion[self.insideIDs,:]+deposition[self.insideIDs,:]

            if verbose:
                print("   - Total sediment flux time ", time.clock() - time0)

        return newdt,sedflux,erosion,deposition

    def gaussian_diffusion(self, diff, dsmooth):
        """
        Gaussian filter operation used to smooth diffusion related deposition
        thicknesses.

        Parameters
        ----------
        diff
            Numpy arrays containing the deposition thicknesses.

        dsmooth
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
            splits = numpy.array_split(xyi, 1)
            split_lengths = numpy.array(list(map(len, splits))) * 3
            localxyi = splits[0]
            query_shape = (xyi.shape[0], 3)

            # Build Tree
            tree = cKDTree(self.xycoords[:,:2])

            # Querying the KDTree is rather slow, so we split it across MPI nodes
            nelems = query_shape[0] * query_shape[1]
            indices = numpy.empty(query_shape, dtype=numpy.int64)
            localdistances, localindices = tree.query(localxyi, k=3)

            self.distances = localdistances
            self.indices = localindices
            self.onIDs = numpy.where(self.distances[:,0] == 0)[0]

        depZ = numpy.copy(diff)

        if len(depZ[self.indices].shape) == 3:
            zd_vals = depZ[self.indices][:,:,0]
        else:
            zd_vals = depZ[self.indices]

        with numpy.errstate(divide='ignore'):
            zdi = numpy.average(zd_vals,weights=(1./self.distances), axis=1)

        if len(self.onIDs) > 0:
            zdi[self.onIDs] = depZ[self.indices[self.onIDs,0]]

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
        elev
            Numpy arrays containing the elevation of the TIN nodes.

        locIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """

        # Compute the local value for time stability
        dt = FLOWalgo.flowcompute.flowcfl(locIDs,self.receivers,self.xycoords,elev, \
                                      self.discharge,self.erodibility)

        # Global mimimum value for diffusion stability
        self.CFL = dt

        return

    def visualise_draining_path(self, pIDs, elev, drain, fillH, filename):
        """
        Debugging function used to plot draining pathway between depressions.

        Parameters
        ----------
        pIDs
            Numpy array containing all depression node IDs.

        elev
            Numpy arrays containing the elevation of the TIN nodes.

        drain
            Numpy arrays containing the draining pit ID.

        filename
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
