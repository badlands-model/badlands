##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines pelagic evolution in **badlands** simulation based on forcing parameter: **depth**.

"""
import os
import numpy
import pandas

from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.spatial import cKDTree

class pelagicGrowth:
    """
    This class defines external pelagic growth parameters.

    Args:
        input: class containing XML input file parameters.
    """
    def __init__(self, input=None):

        self.growth = input.pelGrowth
        self.depthfile = input.pelDepth

        self.depthval = None
        self.depthfct = None
        self.depthFunc = None

        self.depthgrowth = None

        if self.depthfile != None:
            self._build_depth_function()

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

        return

    def _getDepthFct(self, depthfield):
        """
        Computes for a given depth field the pelagic growth function.

        Parameters
        ----------
        depthfield : numpy array containing depth.
        """

        if self.depthfile == None:
            self.depthgrowth = numpy.zeros(len(depthfield))
        else:
            self.depthgrowth = self.depthFunc(-depthfield)

        return

    def computePelagic(self, depthfield, dt):
        """
        Computes pelagic growth.

        Args:
            depthfield : numpy array containing depth.
            dt: pelagic growth time step in years.

        Returns:
            - growth - numpy array containing the growth (in metres) of pelagic.
        """

        # Get each controlling function values
        self._getDepthFct(depthfield)

        # Average growth function limitation
        growth = self.growth*self.depthgrowth*dt
        growth[growth<0.] = 0.

        return growth
