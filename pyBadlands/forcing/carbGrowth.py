##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines several functions used to define carbonate evolution in Badlands
simulation based on 3 forcing parameters: depth, wave and sedimentation rate.
"""
import os
import time
import numpy
import pandas
import triangle
import mpi4py.MPI as mpi
from pyBadlands.libUtils import ORmodel
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import OrderedDict
from matplotlib import _cntr as cntr

class carbGrowth:
    """
    This class defines external carbonate growth parameters.
    """
    def __init__(self, input=None, regX=None, regY=None, boundsPt=None):

        self.regX = regX
        self.regY = regY
        self.boundsPt = boundsPt
        self.tXY = None

        self.growth = input.carbGrowth
        self.depthfile = input.carbDepth
        self.sedfile = input.carbSed
        self.wavefile = input.carbWave

        self.depthval = None
        self.depthfct = None
        self.depthFunc = None
        self.sedval = None
        self.sedfct = None
        self.sedFunc = None
        self.waveval = None
        self.wavefct = None
        self.waveFunc = None

        self.sedgrowth = None
        self.depthgrowth = None
        self.wavegrowth = None

        self.mlen = input.islandPerim
        self.mdist = input.coastdist
        self.Afactor = input.Afactor
        self.tree = None
        self.dx = None
        self.nx = None
        self.ny = None
        self.xi = None
        self.yi = None
        self.distances = None
        self.indices = None

        self.tinBase1 = None
        self.baseMap = input.baseMap

        if self.depthfile != None:
            self._build_depth_function()
        if self.sedfile != None:
            self._build_sed_function()
        if self.wavefile != None:
            self._build_wave_function()

    def build_basement(self,tXY):
        """
        Using Pandas library to read the basement map file and define consolidated and
        soft sediment region.
        """
        self.tXY = tXY
        self.tinBase1 = numpy.ones(len(tXY))

        # Read basement file
        Bmap = pandas.read_csv(str(self.baseMap), sep=r'\s+', engine='c',
                    header=None, na_filter=False, dtype=numpy.float, low_memory=False)

        rectBase = numpy.reshape(Bmap.values,(len(self.regX), len(self.regY)),order='F')
        self.tinBase1[self.boundsPt:] = interpolate.interpn( (self.regX, self.regY), rectBase,
                                                    tXY[self.boundsPt:,:], method='linear')

        return

    def _build_depth_function(self):
        """
        Using Pandas library to read the depth control file and define depth interpolation
        function based on Scipy 1D linear function.
        """

        # Read depth control file
        depthdata = pandas.read_csv(self.depthfile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.depthval = numpy.zeros(len(depthdata.values[:,0])+2)
        self.depthfct = numpy.zeros(len(self.depthval))

        self.depthval[1:-1] = depthdata.values[:,0]
        self.depthfct[1:-1] = depthdata.values[:,1]
        self.depthval[0] = -1.0e7
        self.depthfct[0] = self.depthfct[1]
        self.depthval[-1] = 1.e7
        self.depthfct[-1] = self.depthfct[-2]

        self.depthFunc = interpolate.interp1d(self.depthval, self.depthfct, kind='linear')

    def _build_sed_function(self):
        """
        Using Pandas library to read the sedimentation control file and define sedimentation interpolation
        function based on Scipy 1D linear function.
        """

        # Read sedimentation rate file
        seddata = pandas.read_csv(self.sedfile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.sedval = numpy.zeros(len(seddata.values[:,0])+2)
        self.sedfct = numpy.zeros(len(self.sedval))

        self.sedval[1:-1] = seddata.values[:,0]
        self.sedfct[1:-1] = seddata.values[:,1]
        self.sedval[0] = -1.0e7
        self.sedfct[0] = self.sedfct[1]
        self.sedval[-1] = 1.e7
        self.sedfct[-1] = self.sedfct[-2]

        self.sedFunc = interpolate.interp1d(self.sedval, self.sedfct, kind='linear')

    def _build_wave_function(self):
        """
        Using Pandas library to read the wave control file and define wave interpolation
        function based on Scipy 1D linear function.
        """

        # Read wave control file
        wavedata = pandas.read_csv(self.wavefile, sep=r'\s+', engine='c',
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)

        self.waveval = numpy.zeros(len(wavedata.values[:,0])+2)
        self.wavefct = numpy.zeros(len(self.waveval))

        self.waveval[1:-1] = wavedata.values[:,0]
        self.wavefct[1:-1] = wavedata.values[:,1]
        self.waveval[0] = -1.0e7
        self.wavefct[0] = self.wavefct[1]
        self.waveval[-1] = 1.e7
        self.wavefct[-1] = self.wavefct[-2]

        self.waveFunc = interpolate.interp1d(self.waveval, self.wavefct, kind='linear')

    def _getWaveFct(self, wavefield):
        """
        Computes for a given wave field to carbonate wave dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing wave height.
        """

        if self.wavefile == None:
            self.wavegrowth = numpy.ones(len(wavefield))
        else:
            self.wavegrowth = self.waveFunc(wavefield)

        return

    def _getSedFct(self, sedfield):
        """
        Computes for a given sedimentation rate dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing sedimentation rate.
        """

        if self.sedfile == None:
            self.sedgrowth = numpy.ones(len(sedfield))
        else:
            self.sedgrowth = self.sedFunc(sedfield)

        return

    def _getDepthFct(self, depthfield):
        """
        Computes for a given depth field to carbonate wave dependent growth function.

        Parameters
        ----------
        depthfield : numpy array containing depth.
        """

        if self.depthfile == None:
            self.depthgrowth = numpy.ones(len(depthfield))
        else:
            self.depthgrowth = self.depthFunc(-depthfield)

        return

    def computeShoreline(self,z,lvl=0.):
        """
        This function computes the shoreline position for a given sea-level.
        Parameters
        ----------
        variable: z
            Mesh relative elevation to sea-level.
        variable: lvl
            Water level defined in the input.
        """

        c = cntr.Cntr(self.xi, self.yi, z)
        contour = c.trace(lvl)

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
            pts = numpy.asarray(unique.values())

            if closed:
                cpts = numpy.zeros((len(pts)+1,2), order='F')
                cpts[0:len(pts),0:2] = pts
                cpts[-1,0:2] = pts[0,0:2]

                # Get contour length
                arr = cpts
                val = (arr[:-1,:] - arr[1:,:]).ravel()
                dist = val.reshape((arr.shape[0]-1,2))
                lgth = numpy.sum(numpy.sqrt(numpy.sum(dist**2, axis=1)))
            else:
                lgth = 1.e8
                cpts = pts

            if len(cpts) > 2 and lgth > self.mlen:
                contourList.append(cpts)
                if start:
                    contourPts = cpts
                    start = False
                else:
                    contourPts = numpy.concatenate((contourPts,cpts))

        return contourPts

    def oceanIDs(self,xy,depthfield):

        tree = cKDTree(xy)
        distances, indices = tree.query(self.tXY, k=1)
        seaIDs = numpy.where(numpy.logical_and(distances[:]>=self.mdist,depthfield<=0.))[0]

        return seaIDs

    def buildReg(self,tXY):
        """
        Build regular grid for shoreline contour calculation.
        """

        self.tXY = tXY
        self.tree = cKDTree(self.tXY)
        self.dx = (self.tXY[1,0] - self.tXY[0,0])*self.Afactor

        if self.nx is None:
            self.nx = int((self.tXY[:,0].max() - self.tXY[:,1].min())/self.dx+1)
            self.ny = int((self.tXY[:,1].max() - self.tXY[:,1].min())/self.dx+1)
            xi = numpy.linspace(self.tXY[:,0].min(), self.tXY[:,0].max(), self.nx)
            yi = numpy.linspace(self.tXY[:,1].min(), self.tXY[:,1].max(), self.ny)
            self.xi, self.yi = numpy.meshgrid(xi, yi)
            xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]

        self.distances, self.indices = self.tree.query(xyi, k=3)

        return

    def getDistanceShore(self,depthfield):
        """
        Computes IDs of nodes at a given distance from shoreline.
        """

        if len(depthfield[self.indices].shape) == 3:
            z_vals = depthfield[self.indices][:,:,0]
        else:
            z_vals = depthfield[self.indices]

        zi = numpy.average(z_vals,weights=(1./self.distances), axis=1)
        onIDs = numpy.where(self.distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            zi[onIDs] = depthfield[self.indices[onIDs,0]]

        z = numpy.reshape(zi,(self.ny, self.nx))

        xy = self.computeShoreline(z)

        seaIDs = self.oceanIDs(xy, depthfield)

        return seaIDs

    def computeCarbonate(self, wavefield, sedfield, depthfield, dt):
        """
        Computes carbonate growth.
        """

        if self.mdist == 0.:
            if self.baseMap is not None:
                seaIds = numpy.where(numpy.logical_and(self.tinBase1==0,depthfield<0.))[0]
            else:
                seaIds = numpy.where(depthfield<0.)[0]
        else:
            seaIds = self.getDistanceShore(depthfield)

        growth = numpy.zeros(len(depthfield))
        growth.fill(1.1e6)

        # Get each controlling function values
        if self.depthfile != None:
            self._getDepthFct(depthfield)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.depthgrowth[seaIds])
        if self.sedfile != None:
            self._getSedFct(sedfield)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.sedgrowth[seaIds])
        if self.wavefile != None:
            self._getWaveFct(wavefield)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.wavegrowth[seaIds])
        growth[growth>1.e6] = 0.

        # Average growth function limitation
        val = self.growth*growth*dt
        val[val<0.] = 0.
        val[seaIds] = numpy.minimum(val[seaIds],-depthfield[seaIds]*0.98)
        #tids = numpy.where(numpy.logical_and(depthfield<-20.,val>0.))[0]

        return val
