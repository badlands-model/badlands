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
import numpy
import pandas
import triangle
import mpi4py.MPI as mpi
from pyBadlands.libUtils import ORmodel
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.spatial import cKDTree

class carbGrowth:
    """
    This class defines external carbonate growth parameters.
    """
    def __init__(self, input=None):

        self.growth = input.carbGrowth
        self.ero = input.carbEro
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

        if self.depthfile != None:
            self._build_depth_function()
        if self.sedfile != None:
            self._build_sed_function()
        if self.wavefile != None:
            self._build_wave_function()

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

    def getWaveFct(self, wavefield):
        """
        Computes for a given wave field to carbonate wave dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing wave height.
        """

        if self.wavefile == None:
            self.wavegrowth = np.ones(len(wavefield))
        else:
            self.wavegrowth = self.waveFunc(wavefield)

        return

    def getSedFct(self, sedfield):
        """
        Computes for a given sedimentation rate dependent growth function.

        Parameters
        ----------
        wavefield : numpy array containing sedimentation rate.
        """

        if self.sedfile == None:
            self.sedgrowth = np.ones(len(sedfield))
        else:
            self.sedgrowth = self.sedFunc(sedfield)

        return

    def getDepthFct(self, depthfield):
        """
        Computes for a given depth field to carbonate wave dependent growth function.

        Parameters
        ----------
        depthfield : numpy array containing depth.
        """

        if self.depthfile == None:
            self.depthgrowth = np.ones(len(depthfield))
        else:
            self.depthgrowth = self.depthFunc(-depthfield)

        return

    def computeCarbonate(self, wavefield, sedfield, depthfield):
        """
        Computes carbonate growth.
        """

        # Get each controlling function values
        self.getSedFct(sedfield)
        self.getWaveFct(wavefield)
        self.getDepthFct(depthfield)

        # Average growth function limitation
        growth = self.growth*(self.depthgrowth + self.sedgrowth + self.wavegrowth)/3.

        
        return
