##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to Badlands hillslope computation
based on linear diffusion.
"""

import math
import numpy
import warnings
from pyBadlands.libUtils import FLOWalgo
import mpi4py.MPI as MPI

class diffLinear:
    """
    Class for handling hillslope computation using a linear diffusion equation.
    """

    def __init__(self):
        '''Initialization.'''
        self.CDaerial = None
        self.CDmarine = None
        self.CFL = None

    def dt_stability(self, edgelen):
        """
        This function computes the maximal timestep to ensure computation stability
        of the hillslope processes. This CFL-like condition is computed using diffusion
        coefficients and distances between TIN nodes.
        It is worth noticing that the approach does not rely on the elevation of the nodes
        and therefore the maximal hillslope timestep to ensure stability just needs to be
        computed once for each given TIN grid.

        **Parameters**

        variable : edgelen
            Numpy arrays containing the edges of the TIN surface for the considered partition.
        """

        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Get the tin edges lengths
        maxCD = max(self.CDaerial,self.CDmarine)
        edgedist = edgelen.flatten()
        distIDs = numpy.where(edgedist > 0.)

        # First-order, forward-in-time scheme
        CFL = numpy.zeros(1)
        CFL[0] = 0.05*numpy.amin(edgedist[distIDs]**2)/maxCD

        # Global mimimum value for diffusion stability
        comm.Allreduce(MPI.IN_PLACE,CFL,op=MPI.MIN)
        self.CFL = CFL[0]

    def sedflux(self, diff_flux, sea, elevation, area):
        """
        This function computes the sedimentary fluxes induced by hillslope processes based
        on a linear diffusion approximation.
        The linear diffusion process is implemented through the FV approximation and is based on
        the area of each node voronoi polygon and the sum over all the neighbours of the slope of the
        segment (i.e. height differences divided by the length of the mesh edge) as well as the length
        of the corresponding voronoi edge.

        **Parameters**

        variable : diff_flux
            Numpy arrays representing for each node the sum of the ratio between the height differences
            and  the length of the mesh edge multiply by the length of the corresponding voronoi edge.

        variable : sea
            Real value giving the sea-level height at considered time step.

        variable : elevation
            Numpy arrays containing the elevation of the nodes.

        variable : area
            Numpy arrays containing the area of the voronoi polygon for each TIN nodes.
        """

        ids = numpy.where(area > 0)[0]
        
        areacoeff = numpy.zeros(len(area))
        areacoeff[ids] = 1./area[ids]
        coeff = numpy.where(elevation >= sea, self.CDaerial, self.CDmarine)

        return numpy.nan_to_num(diff_flux * areacoeff * coeff)
