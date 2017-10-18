##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Regional scale model of wave propogation and associated sediment transport. The wave model
is based on Airy wave theory and takes into account wave refraction based on Huygen's principle.
The sediment entrainment is computed from wave shear stress and transport according to both
wave direction and longshore transport. Deposition is dependent of shear stress and diffusion.
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
from matplotlib import _cntr as cntr

import pyBadlands.libUtils.WAVEsed as ocean

from random import *

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class waveSed():
    """
    Class for building wave based on linear wave theory.
    """

    def __init__(self, input, recGrid, Ce, Cd):
        """
        Initialization function.

        Parameters
        ----------
        variable : filename
            Bathymetry file to load.

        variable : wavebase
            Maximum depth for wave influence (m) [Default is 10].

        variable: resfac
            Requested resolution factor for wave and sediment transport computation.

        variable: dia
            Sediment average diameter (m).

        variable: Ce
            Sediment entrainment coefficient [Default is 1.].

        variable: Cd
            Sediment diffusion coefficient [Default is 30.].
        """

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
        self.Ero = input.wEro

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
            self.tau_cr = 0.045

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

        Parameters
        ----------
        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m)
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

        self.sealvl = force.sealevel

        # t1 = time.clock()
        self.findland(elev, actlay, self.sealvl)
        # print 'findland (s)',time.clock()-t1
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
                    source = self.wavesource(direction)

                    # Compute wave parameters for given condition
                    self.cmptwaves(source, h0=height, sigma=1.)
                    # print 'cmptwaves (s)',time.clock()-t1
                    t1 = time.clock()

                    # Compute sediment transport
                    self.cmptsed(perc,sigma=1.)
                    # print 'cmptsed (s)',time.clock()-t1

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

        return tED,nactlay

    def wavesource(self, dir=0.):
        """
        This function defines wave source boundary conditions from input directions.

        Parameters
        ----------

        variable: dir
            Wave direction from input condition.
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

    def compute_shoreline(self, mlen=0.):
        """
        This function computes the shoreline position for a given sea-level.

        Parameters
        ----------

        variable: mlen
            Minimum island perimeter length to consider.
        """

        c = cntr.Cntr(self.xi, self.yi, self.regZ.T)
        contour = c.trace(self.sealvl)

        nseg = len(contour) // 2
        contours, codes = contour[:nseg], contour[nseg:]
        contourList = []
        start = True

        # Loop through each contour
        for c in range(len(contours)):
            tmpts =  contours[c]
            closed = False
            if tmpts[0,0] == tmpts[-1,0] and tmpts[0,1] == tmpts[-1,1]:
                closed = True

            # Remove duplicate points
            unique = OrderedDict()
            for p in zip(tmpts[:,0], tmpts[:,1]):
                unique.setdefault(p[:2], p)
            pts = np.asarray(unique.values())

            if closed:
                cpts = np.zeros((len(pts)+1,2), order='F')
                cpts[0:len(pts),0:2] = pts
                cpts[-1,0:2] = pts[0,0:2]

                # Get contour length
                arr = cpts
                val = (arr[:-1,:] - arr[1:,:]).ravel()
                dist = val.reshape((arr.shape[0]-1,2))
                lgth = np.sum(np.sqrt(np.sum(dist**2, axis=1)))
            else:
                lgth = 1.e8
                cpts = pts

            if len(cpts) > 2 and lgth > mlen:
                contourList.append(cpts)
                if start:
                    contourPts = cpts
                    start = False
                else:
                    contourPts = np.concatenate((contourPts,cpts))

        return contourPts, contourList

    def findland(self, elev, actlay, lvl=0.):
        """
        This function computes the land IDs as well as the lake IDs.

        Parameters
        ----------
        variable: elev
            Elevation of the TIN.

        variable: lvl
            Sea-level position.
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

        # xy, xylist = self.compute_shoreline(0.)
        #self.regZ.dump("elev.dat")

        # Find lakes and assign IDs as land
        # t11 = time.clock()
        # for p in range(len(xylist)):
        #     if xylist[p][0,0] == xylist[p][-1,0] and xylist[p][0,1] == xylist[p][-1,1]:
        #         mpath = Path( xylist[p] )
        #         mask_flat = mpath.contains_points(self.XY)
        #         ar = mask_flat.reshape(self.xi.shape).astype(int)
        #         tc,tr = np.where(np.logical_and(ar.T==1,self.regZ<self.sealvl))
        #         self.inland[tc,tr] = 1
        # print 'lakes (s)',time.clock()-t11,len(xylist)
        return

    def cmptwaves(self, src=None, h0=0., sigma=1., shadow=0, shoalC=0.99):
        """
        Waves are transformed from deep to shallow water assuming shore-parallel depth contours. The
        orientation of wave fronts is determine by wave celerity and refraction due to depth variations
        and travel time in the domain is calculated from Huygen's principle.

        Parameters
        ----------
        variable: src
            Position of wave boundary condition.

        variable: h0
            Wave height value along the boundary.

        variable: sigma
            Smoothing coefficient.

        variable: shadow
            Considering shadow effect (1) or no shadow (0) [default is 0].

        variable: shoalC
            Coefficent at attenuation in shoaling region [default is 0.99].
        """

        waveC,waveL,travel,self.waveH = ocean.wavesed.airymodel(self.dx,
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
        R = waveU*waveT/(2.*np.pi*2.5*self.dia)
        tmpr2,tmpc2 = np.where(R>0.)
        fric[tmpr2,tmpc2] = 0.237*np.power(R[tmpr2,tmpc2],-0.52)

        # Shear stress (N/m2)
        self.waveS = 0.5*self.rhow*fric*np.power(waveU,2.)
        self.waveS = gaussian_filter(self.waveS, sigma=sigma)
        self.waveS[self.landr,self.landc] = 0.

        return

    def cmptsed(self, perc, sigma=1.):
        """
        Compute wave induced sedimentation (erosion/deposition).

        Parameters
        ----------
        variable: sigma
            Smoothing coefficient.
        """

        # Thickness of entrained sediment [L]
        self.waveS[self.waveS<1.e-5] = 0.
        r,c=np.where(self.waveS>0.)
        Hent = np.zeros(self.waveS.shape)
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
        wdz,distw = ocean.wavesed.transport(self.tsteps,self.depth,Hent,tX,tY)

        # Diffuse marine coefficient
        area = self.dx*self.dx
        CFL = float(area*area/(4.*self.Cd*area))
        Cdiff = self.Cd/area

        # Perform wave related sediment diffusion
        nelev = -self.depth+wdz-Hent

        # Compute maximum marine fluxes and maximum timestep to avoid excessive diffusion erosion
        ndz = ocean.wavesed.diffusion(nelev, wdz, Cdiff, self.Ero, CFL, self.dsteps)

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
        #self.erodep.dump("erodep.dat")

        return
