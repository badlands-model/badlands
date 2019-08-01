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
from pyBadlands.libUtils import ORmodel
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import OrderedDict
from legacycontour import _cntr as cntr

class carbGrowth:
    """
    This class defines external carbonate growth parameters.
    """
    def __init__(self, input=None, regX=None, regY=None, tinBase=None):

        self.regX = regX
        self.regY = regY
        self.tXY = None

        self.depthfile = input.carbDepth
        self.sedfile = input.carbSed
        self.wavefile = input.carbWave

        self.depthfile2 = input.carbDepth2
        self.sedfile2 = input.carbSed2
        self.wavefile2 = input.carbWave2

        self.depthval = None
        self.depthfct = None
        self.depthFunc = None
        self.sedval = None
        self.sedfct = None
        self.sedFunc = None
        self.waveval = None
        self.wavefct = None
        self.waveFunc = None

        self.depthval2 = None
        self.depthfct2 = None
        self.depthFunc2 = None
        self.sedval2 = None
        self.sedfct2 = None
        self.sedFunc2 = None
        self.waveval2 = None
        self.wavefct2 = None
        self.waveFunc2 = None

        self.carbonate = input.carbonate
        self.sedgrowth = None
        self.depthgrowth = None
        self.wavegrowth = None

        self.carbonate2 = input.carbonate2
        self.sedgrowth2 = None
        self.depthgrowth2 = None
        self.wavegrowth2 = None

        self.mlen = input.islandPerim
        self.mdist = input.coastdist

        self.mlen2 = input.islandPerim2
        self.mdist2 = input.coastdist2

        self.Afactor = input.Afactor
        self.tree = None
        self.dx = None
        self.nx = None
        self.ny = None
        self.xi = None
        self.yi = None
        self.distances = None
        self.indices = None

        self.tinBase = tinBase

        if self.depthfile != None:
            self._build_depth_function(1)
        if self.sedfile != None:
            self._build_sed_function(1)
        if self.wavefile != None:
            self._build_wave_function(1)

        if self.depthfile2 != None:
            self._build_depth_function(2)
        if self.sedfile2 != None:
            self._build_sed_function(2)
        if self.wavefile2 != None:
            self._build_wave_function(2)

    def _build_depth_function(self,id):
        """
        Using Pandas library to read the depth control file and define depth interpolation
        function based on Scipy 1D linear function.
        """

        # Read depth control file
        if id == 1:
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

        if id == 2:
            depthdata = pandas.read_csv(self.depthfile2, sep=r'\s+', engine='c',
                                   header=None, na_filter=False,
                                   dtype=numpy.float, low_memory=False)

            self.depthval2 = numpy.zeros(len(depthdata.values[:,0])+2)
            self.depthfct2 = numpy.zeros(len(self.depthval))

            self.depthval2[1:-1] = depthdata.values[:,0]
            self.depthfct2[1:-1] = depthdata.values[:,1]
            self.depthval2[0] = -1.0e7
            self.depthfct2[0] = self.depthfct2[1]
            self.depthval2[-1] = 1.e7
            self.depthfct2[-1] = self.depthfct2[-2]

            self.depthFunc2 = interpolate.interp1d(self.depthval2, self.depthfct2, kind='linear')

    def _build_sed_function(self,id):
        """
        Using Pandas library to read the sedimentation control file and define sedimentation interpolation
        function based on Scipy 1D linear function.
        """

        # Read sedimentation rate file
        if id == 1:
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

        if id == 2:
            seddata = pandas.read_csv(self.sedfile2, sep=r'\s+', engine='c',
                                   header=None, na_filter=False,
                                   dtype=numpy.float, low_memory=False)

            self.sedval2 = numpy.zeros(len(seddata.values[:,0])+2)
            self.sedfct2 = numpy.zeros(len(self.sedval))

            self.sedval2[1:-1] = seddata.values[:,0]
            self.sedfct2[1:-1] = seddata.values[:,1]
            self.sedval2[0] = -1.0e7
            self.sedfct2[0] = self.sedfct2[1]
            self.sedval2[-1] = 1.e7
            self.sedfct2[-1] = self.sedfct2[-2]

            self.sedFunc2 = interpolate.interp1d(self.sedval2, self.sedfct2, kind='linear')

    def _build_wave_function(self, id):
        """
        Using Pandas library to read the wave control file and define wave interpolation
        function based on Scipy 1D linear function.
        """

        # Read wave control file
        if id == 1:
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

        if id == 2:
            wavedata = pandas.read_csv(self.wavefile2, sep=r'\s+', engine='c',
                                   header=None, na_filter=False,
                                   dtype=numpy.float, low_memory=False)

            self.waveval2 = numpy.zeros(len(wavedata.values[:,0])+2)
            self.wavefct2 = numpy.zeros(len(self.waveval2))

            self.waveval2[1:-1] = wavedata.values[:,0]
            self.wavefct2[1:-1] = wavedata.values[:,1]
            self.waveval2[0] = -1.0e7
            self.wavefct2[0] = self.wavefct2[1]
            self.waveval2[-1] = 1.e7
            self.wavefct2[-1] = self.wavefct2[-2]

            self.waveFunc2 = interpolate.interp1d(self.waveval2, self.wavefct2, kind='linear')

    def _getWaveFct(self, wavefield, id):
        """
        Computes for a given wave field to carbonate wave dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing wave height.
        """

        if id == 1:
            if self.wavefile == None:
                self.wavegrowth = numpy.ones(len(wavefield))
            else:
                self.wavegrowth = self.waveFunc(wavefield)

        if id == 2:
            if self.wavefile2 == None:
                self.wavegrowth2 = numpy.ones(len(wavefield))
            else:
                self.wavegrowth2 = self.waveFunc2(wavefield)

        return

    def _getSedFct(self, sedfield, id):
        """
        Computes for a given sedimentation rate dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing sedimentation rate.
        """

        if id == 1:
            if self.sedfile == None:
                self.sedgrowth = numpy.ones(len(sedfield))
            else:
                self.sedgrowth = self.sedFunc(sedfield)

        if id == 2:
            if self.sedfile2 == None:
                self.sedgrowth2 = numpy.ones(len(sedfield))
            else:
                self.sedgrowth2 = self.sedFunc2(sedfield)

        return

    def _getDepthFct(self, depthfield, id):
        """
        Computes for a given depth field to carbonate wave dependent growth function.

        Parameters
        ----------
        depthfield : numpy array containing depth.
        """

        if id == 1:
            if self.depthfile == None:
                self.depthgrowth = numpy.ones(len(depthfield))
            else:
                self.depthgrowth = self.depthFunc(-depthfield)

        if id == 2:
            if self.depthfile2 == None:
                self.depthgrowth2 = numpy.ones(len(depthfield))
            else:
                self.depthgrowth2 = self.depthFunc2(-depthfield)

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
            pts = numpy.asarray(list(unique.values()))

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

    def computeCarbonate(self, wavefield, sedfield, depthfield, growthsp1, growthsp2, dt):
        """
        Computes carbonate growth.
        """

        if self.mdist == 0.:
            if self.tinBase is not None:
                tmpids = numpy.where(depthfield<0.)[0]
                seaIds = numpy.where(numpy.logical_and(self.tinBase==0,depthfield<0.))[0]
            else:
                seaIds = numpy.where(depthfield<0.)[0]
        else:
            seaIds = self.getDistanceShore(depthfield)

        growth = numpy.zeros(len(depthfield))
        growth.fill(1.1e6)

        if self.carbonate2:
            growth2 = numpy.zeros(len(depthfield))
            growth2.fill(1.1e6)

        # Get each controlling function values
        if self.depthfile != None:
            self._getDepthFct(depthfield,1)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.depthgrowth[seaIds])
        if self.sedfile != None:
            self._getSedFct(sedfield,1)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.sedgrowth[seaIds])
        if self.wavefile != None:
            self._getWaveFct(wavefield,1)
            growth[seaIds] = numpy.minimum(growth[seaIds],self.wavegrowth[seaIds])
        growth[growth>1.e6] = 0.

        if self.carbonate2:
            if self.depthfile2 != None:
                self._getDepthFct(depthfield,2)
                growth2[seaIds] = numpy.minimum(growth2[seaIds],self.depthgrowth2[seaIds])
            if self.sedfile2 != None:
                self._getSedFct(sedfield,2)
                growth2[seaIds] = numpy.minimum(growth2[seaIds],self.sedgrowth2[seaIds])
            if self.wavefile2 != None:
                self._getWaveFct(wavefield,2)
                growth2[seaIds] = numpy.minimum(growth2[seaIds],self.wavegrowth2[seaIds])
            growth2[growth2>1.e6] = 0.

        # Average growth function limitation
        val = growthsp1*growth*dt
        val[val<0.] = 0.
        val[seaIds] = numpy.minimum(val[seaIds],-depthfield[seaIds]*0.9)
        tmpid = numpy.where(numpy.logical_and(val==val.max(),val>0))[0]
        if self.carbonate2:
            val2 = growthsp2*growth2*dt
            val2[val2<0.] = 0.
            val2[seaIds] = numpy.minimum(val2[seaIds],-depthfield[seaIds]*0.9)
        else:
            val2 = None

        return val, val2
