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
based on non-linear diffusion.
"""

import math
import numpy
import warnings
from pyBadlands.libUtils import FLOWalgo
import mpi4py.MPI as MPI

class diffnLinear:
    """
    Class for handling hillslope computation using a non-linear diffusion equation.
    """

    def __init__(self):

        '''Initialization.
        '''
        self.CDaerial = None
        self.CDmarine = None
        self.Sc = None
        self.CFL = None

        return

    def dt_stability(self, diffcfl):
        """
        This function computes the maximal timestep to ensure computation stability
        of the non-linear hillslope processes.

        Parameters
        ----------
        variable : diffcfl
            Numpy arrays containing the cfl-like condition at each nodes of the considered partition.
        """

        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Get the tin edges lengths
        maxCD = max(self.CDaerial,self.CDmarine)

        cflIDs = numpy.where(diffcfl > 0.)

        # First-order, forward-in-time scheme
        CFL = numpy.zeros(1)
        CFL[0] = 0.1*numpy.amin(diffcfl[cflIDs])/maxCD

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

        Parameters
        ----------
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

        flux = numpy.zeros(len(elevation))
        tmpIDs = numpy.where(area > 0)
        flux[tmpIDs] = diff_flux[tmpIDs] / area[tmpIDs]
        iddd = numpy.where(flux==flux.max())

        dryIDs = numpy.where(elevation >= sea)
        flux[dryIDs] = flux[dryIDs] * self.CDaerial

        wetIDs = numpy.where(elevation < sea)
        flux[wetIDs] = flux[wetIDs] * self.CDmarine

        return flux
