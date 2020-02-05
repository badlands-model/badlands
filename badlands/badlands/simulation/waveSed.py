##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Regional scale model of wave propogation and associated sediment transport.

Important:
    The wave model is based on **Airy wave theory** and takes into account wave refraction based on
    **Huygen's principle**.

Airy wave theory:

.. image:: img/airy.png
   :scale: 90 %
   :alt: airy wave theory
   :align: center

Huygen's principle:

.. image:: img/huygens.jpg
   :scale: 110 %
   :alt: Huygens theory
   :align: center

The sediment entrainment is computed from wave shear stress and transport according to both
wave direction and longshore transport. Deposition is dependent of shear stress and diffusion.

Note:
    The model is intended to quickly simulate the average impact of wave induced sediment transport
    at large scale and over geological time period.

"""

import os
import math
import time
import errno
import numpy as np
import pandas as pd

from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import interpn
from scipy.ndimage.filters import gaussian_filter

from matplotlib.path import Path
from collections import OrderedDict

if 'READTHEDOCS' not in os.environ:
    from random import *
    from badlands import waveseds as ocean

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class waveSed:
    """
    Class for building wave based on **linear wave theory**.


    Initialisation function.

    Args:
        input: class containing XML input file parameters.
        recGrid: class describing the regular grid characteristics.
        Ce: sediment entrainment coefficient [Default is 1.].
        Cd: sediment diffusion coefficient [Default is 30.].
    """
    def __init__(self, input, recGrid, Ce, Cd):

        # Gravity [L/T2]
        self.grav = 9.81
        # Sea water density [M/L3]
        self.rhow = 1027
        # Sediment density [M/L3]
        self.rhos = 2650
        # Porosity
        self.poro = 0.4
        # Bottom friction coefficient
        self.fric = 0.032
        # Kinematic viscosity water (20C) [m2/s]
        self.nu = 1.004*1.e-6

        self.dia = input.d50
        self.tsteps = input.tsteps
        self.dsteps = input.dsteps
        self.Ce = Ce
        self.Cd = Cd
        # Maximum eroded thickness for each time step
        self.Ero = input.wEro*input.tWave


        # Non-dimensional diameter
        self.ds = self.dia*np.power(self.grav*(self.rhos/self.rhow-1)/(self.nu*self.nu),1./3.)

        # Van Rijn formula
        if self.ds <= 4.:
            self.tau_cr = 0.24*np.power(self.ds,-1.)
        elif self.ds<= 10.:
            self.tau_cr = 0.14*np.power(self.ds,-0.64)
        elif self.ds<= 20.:
            self.tau_cr = 0.04*np.power(self.ds,-0.1)
        elif self.ds<= 150.:
            self.tau_cr = 0.013*np.power(self.ds,0.29)
        else:
            self.tau_cr = 0.055

        self.tau_cr = self.tau_cr*self.grav*self.dia*(self.rhos-self.rhow)
        self.wavebase = input.waveBase
        self.resfac = input.resW

        minX = recGrid.rectX.min()-input.resW
        maxX = recGrid.rectX.max()+input.resW
        minY = recGrid.rectY.min()-input.resW
        maxY = recGrid.rectY.max()+input.resW
        self.dx = float(input.resW)
        self.regX = np.arange(minX,maxX+self.dx,self.dx)
        self.regY = np.arange(minY,maxY+self.dx,self.dx)

        self.nx = len(self.regX)
        self.ny = len(self.regY)
        self.regZ = np.zeros((self.nx,self.ny),order='F')

        self.xi, self.yi = np.meshgrid(self.regX, self.regY)
        self.XY = np.column_stack((self.xi.flatten(),self.yi.flatten()))
        self.wtree2 = cKDTree(self.XY)

        self.sealvl = 0.
        self.inland = None
        self.depth = None
        self.sear = None
        self.seac = None
        self.landc = None
        self.landr = None
        self.transX = None
        self.transY = None
        self.erodep = None
        self.waveS = None
        self.waveH = None
        self.dists = None
        self.inds = None
        self.wxyTIN = None
        self.wtree = None
        self.dists2 = None
        self.inds2 = None

        self.regularlayer = None

        return

    def build_tree(self, xyTIN):
        """
        Update wave mesh.

        Args:
            xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        """

        # Update TIN grid kdtree for interpolation
        self.wxyTIN = xyTIN
        self.wtree = cKDTree(xyTIN)
        tindx = xyTIN[1,0] - xyTIN[0,0]
        schpts = max(int(self.dx*self.dx/(tindx*tindx)),4)
        self.dists, self.inds = self.wtree.query(self.XY, k=schpts)

        self.dists2, self.inds2 = self.wtree2.query(self.wxyTIN, k=4)

        return

    def compute_wavesed(self,tNow,input,force,elev,actlay):
        """
        This function computes **wave evolution** and wave-induced **sedimentary changes**.

        Args:
            tNow: current simulation time.
            input: class containing XML input file parameters.
            force: class containing wave forcing parameters.
            elev: elevation of the TIN.
            actlay: active layer from TIN.

        Returns
        -------
        tED
            numpy array containing wave induced erosion/deposition changes.
        nactlay
            numpy array containing the updated active layer rock-type content.
        """

        self.sealvl = force.sealevel

        self._findland(elev, actlay, self.sealvl)

        inside = 0
        if actlay is not None:
            nactlay = np.copy(actlay)
        else:
            nactlay = None

        for w in range(input.waveNb):
            if tNow >= input.waveTime[w,0] and tNow < input.waveTime[w,1]:
                for clim in range(input.climNb[w]):
                    perc = input.wavePerc[w][clim]
                    direction = input.waveWd[w][clim]
                    height = input.waveWh[w][clim]

                    # Define wave source direction
                    source = self._wavesource(direction)

                    # Compute wave parameters for given condition
                    self._cmptwaves(source, h0=height, sigma=1.)
                    t1 = time.clock()

                    # Compute sediment transport
                    self._cmptsed(perc,sigma=1.)

                    if clim > 0:
                        avedz += self.erodep
                        avewH += self.waveH*perc
                        avewS += self.waveS*perc
                    else:
                        inside = 1
                        avedz = self.erodep
                        avewH = self.waveH*perc
                        avewS = self.waveS*perc

        # Interpolate to TIN nodes
        if inside > 0:
            # Interpolate wave mesh information on TIN
            h = avewH.flatten('F')
            s = avewS.flatten('F')
            ed = avedz.flatten('F')
            h_vals = h[self.inds2]
            s_vals = s[self.inds2]
            ed_vals = ed[self.inds2]

            with np.errstate(invalid='ignore'):
                force.meanH = np.average(h_vals,weights=(1./self.dists2), axis=1)
                force.meanS = np.average(s_vals,weights=(1./self.dists2), axis=1)
                tED = np.average(ed_vals,weights=(1./self.dists2), axis=1)

            onIDs = np.where(self.dists2[:,0] == 0)[0]
            if len(onIDs) > 0:
                force.meanH[onIDs] = h[self.inds2[onIDs,0]]
                force.meanS[onIDs] = s[self.inds2[onIDs,0]]
                tED[onIDs] = ed[self.inds2[onIDs,0]]

            if actlay is not None:
                al = self.regularlayer.flatten('F')
                al_vals = al[self.inds2]
                with np.errstate(invalid='ignore'):
                    tal = np.average(al_vals,weights=(1./self.dists2), axis=1)
                if len(onIDs) > 0:
                    tal[onIDs] = al[self.inds2[onIDs,0]]
                tal[tal<0.] = 0.
                nactlay[:,0] = tal
            # force.meanH = interpn( (self.regX, self.regY), avewH, (self.wxyTIN), method='linear')
            # force.meanS = interpn( (self.regX, self.regY), avewS, (self.wxyTIN), method='linear')
            # tED = interpn( (self.regX, self.regY), avedz, (self.wxyTIN), method='linear')
        else:
            force.meanH = np.zeros(len(self.wxyTIN))
            force.meanS = np.zeros(len(self.wxyTIN))
            tED = np.zeros(len(self.wxyTIN))

        return tED, nactlay

    def _wavesource(self, dir=0.):
        """
        This function defines wave source boundary conditions from input directions.

        Args:
            dir: wave direction from input condition.

        Returns:
            - src - numpy array containing the imposed wave boundary conditions.
        """

        src = np.zeros(self.regZ.shape)

        src.fill(-2)
        # East
        if dir == 0:
            src[-1,:] = 0
        # North
        elif dir == 90:
            src[:,-1] = 0
        # West
        elif dir == 180:
            src[0,:] = 0
        # South
        elif dir == 270:
            src[:,0] = 0
        # North-East
        elif dir > 0 and dir < 90:
            src[-1,-1] = 0
        # North-West
        elif dir > 90 and dir < 180:
            src[0,-1] = 0
        # South-West
        elif dir > 180 and dir < 270:
            src[0,0] = 0
        # South-East
        elif dir > 270:
            src[-1,0] = 0

        src[self.landr,self.landc] = -2

        return src

    def _findland(self, elev, actlay, lvl=0.):
        """
        This function computes the land IDs as well as the lake IDs.

        Args:
            elev: elevation of the TIN.
            actlay: active layer from TIN.
            lvl: sea-level position.
        """

        self.sealvl = lvl

        # Interpolate TIN elevation on wave mesh
        z_vals = elev[self.inds]
        with np.errstate(invalid='ignore'):
            nz = np.average(z_vals,weights=(1./self.dists), axis=1)
        onIDs = np.where(self.dists[:,0] == 0)[0]
        if len(onIDs) > 0:
            nz[onIDs] = elev[self.inds[onIDs,0]]
        self.regZ = np.reshape(nz,(self.nx,self.ny),order='F')

        # Interpolate TIN active layer composition on wave mesh
        if actlay is not None:
            a_vals = actlay[self.inds,0]
            with np.errstate(invalid='ignore'):
                na = np.average(a_vals,weights=(1./self.dists), axis=1)
            if len(onIDs) > 0:
                na[onIDs] = actlay[self.inds[onIDs,0],0]
            self.regularlayer = np.reshape(na,(self.nx,self.ny),order='F')

        # Specify land/sea areas
        tmpxi = self.xi
        tmpXY = self.XY

        self.inland = np.ones(self.regZ.shape)
        self.depth = self.sealvl - self.regZ
        self.sear,self.seac = np.where(self.depth>0)
        self.inland[self.sear,self.seac]=0
        self.landr,self.landc = np.where(self.depth<=0)

        return

    def _cmptwaves(self, src=None, h0=0., sigma=1., shadow=0, shoalC=0.99):
        """
        Waves are transformed from deep to shallow water assuming shore-parallel depth contours. The
        orientation of wave fronts is determine by wave celerity and refraction due to depth variations
        and travel time in the domain is calculated from Huygen's principle.

        Args:
            src: position of wave boundary condition.
            h0: wave height value along the boundary.
            sigma: smoothing coefficient.
            shadow: considering shadow effect (1) or no shadow (0) [default is 0].
            shoalC: coefficent at attenuation in shoaling region [default is 0.99].
        """

        waveC,waveL,travel,self.waveH = ocean.airymodel(self.dx,
                                                        shoalC,h0,self.depth,
                                                        src,self.inland,shadow)

        lkr,lkc = np.where(np.logical_and(travel<0,self.depth>0))
        self.waveH[lkr,lkc] *= 0.05
        self.waveH[self.waveH>h0*1.25] = h0*1.25

        #self.waveU = np.copy(travel)
        #waveC = gaussian_filter(waveC, sigma=sigma)
        waveL = gaussian_filter(waveL, sigma=sigma)
        waveL[self.landr,self.landc] = 0.

        # Breaking wave height
        Hb = np.zeros(waveC.shape)
        # McCowan (1894)
        Hb[self.sear,self.seac] = 0.78*self.depth[self.sear,self.seac]

        # Wave height [L]
        self.waveH = gaussian_filter(self.waveH, sigma=sigma)
        self.waveH[self.landr,self.landc] = 0.
        breakr,breakc = np.where(self.waveH>Hb)
        self.waveH[breakr,breakc] = Hb[breakr,breakc]

        # Wave direction [radians]
        travel[travel<0.] = travel.max()+10.
        gradx,grady = np.gradient(travel,edge_order=2)
        waveD = np.arctan2(grady, gradx)
        waveD = waveD%(np.pi*2)

        # Wave maximum orbital velocity [L/T]
        waveU = np.zeros(waveC.shape)
        tmp3 = np.sqrt(self.grav/self.depth[self.sear,self.seac])
        waveU[self.sear,self.seac] = 0.5*self.waveH[self.sear,self.seac]*tmp3
        waveU = gaussian_filter(waveU, sigma=sigma)
        tmpr1,tmpc1 = np.where(self.depth>self.wavebase)
        waveU[tmpr1,tmpc1] = 0.

        # Bathymetric contour angle
        gradx,grady = np.gradient(self.regZ,edge_order=2)
        cDir = np.arctan2(grady, gradx)+np.pi/2.
        cDir = cDir%(np.pi*2)

        # Wave transport direction
        self.transpX = np.cos(waveD)
        self.transpY = np.sin(waveD)

        # Longshore drift contour
        tr,tc = np.where(abs(waveD-cDir)>0.5*np.pi)
        cDir[tr,tc] = cDir[tr,tc]+np.pi

        # Sediment transport direction
        tmpr,tmpc = np.where(np.logical_and(self.depth>0.,self.depth<self.wavebase*0.5))
        self.transpX[tmpr,tmpc] = np.cos(cDir[tmpr,tmpc])
        self.transpY[tmpr,tmpc] = np.sin(cDir[tmpr,tmpc])
        self.transpX[self.landr,self.landc] = 0.
        self.transpY[self.landr,self.landc] = 0.

        # Wave period [T]
        waveT = np.zeros(waveC.shape)
        tmpr1,tmpc1 = np.where(self.depth>0.5)
        waveT[tmpr1,tmpc1] = waveL[tmpr1,tmpc1]/np.sqrt(self.grav*self.depth[tmpr1,tmpc1])

        # Friction factor
        fric = np.zeros(waveC.shape)
        kbb = 2.*np.pi*self.dia/12.
        R = waveU*waveT/(2.*np.pi*kbb)
        tmpr2,tmpc2 = np.where(R>0.)
        fric[tmpr2,tmpc2] = 1.39*np.power(R[tmpr2,tmpc2],-0.52)

        # Shear stress (N/m2)
        self.waveS = 0.5*self.rhow*fric*np.power(waveU,2.)
        self.waveS = gaussian_filter(self.waveS, sigma=sigma)
        self.waveS[self.landr,self.landc] = 0.

        return

    def _cmptsed(self, perc, sigma=1.):
        """
        Compute wave induced sedimentation (erosion/deposition).

        Args:
            perc: percentage of activity for the considered wave climate
            sigma: smoothing coefficient.
        """

        # Thickness of entrained sediment [L]
        self.waveS[self.waveS<1.e-5] = 0.
        r,c=np.where(self.waveS>0.)
        Hent = np.zeros(self.waveS.shape)
        #Hent[r,c] = self.Ce*np.log(self.waveS[r,c]/self.tau_cr)*perc
        Hent[r,c] = -self.Ce*np.log(np.sqrt(np.power(self.tau_cr/self.waveS[r,c],2)))*perc
        Hent[Hent<0.] = 0.
        r,c=np.where(np.logical_and(Hent>0.,Hent>0.25*self.depth))
        Hent[r,c] = 0.25*self.depth[r,c]
        if sigma>0:
            Hent = gaussian_filter(Hent, sigma=sigma)
        Hent[self.landr,self.landc] = 0.
        Hent[Hent>self.Ero] = self.Ero

        if np.max(Hent)==0.:
            self.erodep = np.zeros(self.waveS.shape)
            return

        # Limit erosion thickness based on active layer composition
        if self.regularlayer is not None:
            tr1,tc1 = np.where(Hent>self.regularlayer)
            if len(tr1)>0:
                Hent[tr1,tc1] = self.regularlayer[tr1,tc1]

        # Proportion of transport in X,Y direction
        tot = np.abs(self.transpX)+np.abs(self.transpY)
        tr,tc = np.where(tot>0)
        tX = np.zeros(self.waveS.shape)
        tY = np.zeros(self.waveS.shape)
        tX[tr,tc] = self.transpX[tr,tc]/tot[tr,tc]
        tY[tr,tc] = self.transpY[tr,tc]/tot[tr,tc]

        # Compute sediment transport
        wdz,distw = ocean.wavtransport(self.tsteps,self.depth,Hent,tX,tY)

        # Diffuse marine coefficient
        area = self.dx*self.dx
        CFL = float(area*area/(4.*self.Cd*area))
        Cdiff = self.Cd/area

        # Perform wave related sediment diffusion
        nelev = -self.depth+wdz-Hent

        # Compute maximum marine fluxes and maximum timestep to avoid excessive diffusion erosion
        ndz = ocean.wavdiffusion(nelev, wdz, Cdiff, self.Ero, CFL, self.dsteps)

        # Distribute sediment
        if sigma>0.:
            val = gaussian_filter(ndz+distw, sigma=sigma)
            totval = np.sum(val)
            if totval>0.:
                frac = np.sum(ndz+distw)/totval
            else:
                frac = 1.
            val = frac*val
        else:
            val = ndz+distw

        dz = val-Hent

        r,c = np.where(np.logical_and(dz>0,self.depth<-2.))
        dz[r,c] = 0.

        if self.regularlayer is not None:
            self.regularlayer += dz
            r2,c2 = np.where(self.regularlayer<0.)
            self.regularlayer[r2,c2] = 0.
            dz[r2,c2] = 0.

        self.erodep = dz

        return
