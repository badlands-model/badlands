##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to **badlands** hillslope computation based on **diffusion**.


.. image:: img/diff.png
   :scale: 50 %
   :alt: diffusion
   :align: center

Erosion/deposition induced after 130 ka of hillslope diffusion using the *linear* and *non-linear* formulations.
Left: Linear diffusion produces convex upward slopes (κhl = κhn = 0.05). Right: non-linear approach tends to have convex to planar profiles as hillslope processes dominate when slopes approach or exceed the critical slope (Sc = 0.8)


"""

import math
import numpy
import warnings
import os
if 'READTHEDOCS' not in os.environ:
    from badlands import sfd

class diffLinear:
    """
    Class for handling hillslope computation using a linear/non-linear diffusion equation.
    """

    def __init__(self):
        '''Initialization.'''
        self.CDaerial = None
        self.CDmarine = None
        self.CDriver = None
        self.CFL = None
        self.CFLms = None
        self.CFLfail = None
        self.ids = None
        self.Sc = 0.
        self.Sfail = None
        self.Cfail = None
        self.updatedt = 0

        return

    def dt_stability(self, edgelen):
        """
        This function computes the maximal timestep to ensure computation stability
        of the hillslope processes. This CFL-like condition is computed using diffusion
        coefficients and distances between TIN nodes.

        Note:
            It is worth noticing that the approach does not rely on the elevation of the nodes
            and therefore the maximal hillslope timestep to ensure stability just needs to be
            computed once for each given TIN grid.

        Args:
            edgelen: numpy array containing the edges of the TIN surface for the considered partition.
        """

        # Get the tin edges lengths
        maxCD = max(self.CDaerial,self.CDmarine)
        edgedist = edgelen.flatten()
        distIDs = numpy.where(edgedist > 0.)

        # First-order, forward-in-time scheme
        CFL = numpy.zeros(1)
        if maxCD > 0.:
            CFL[0] = 0.05*numpy.amin(edgedist[distIDs]**2)/maxCD
        else:
            CFL[0] = 1.e6

        # Global mimimum value for diffusion stability
        self.CFL = CFL[0]
        self.updatedt = 1

        return

    def dt_stabilityCs(self, elev, neighbours, distances, globalIDs, borders):
        """
        This function computes the maximal timestep to ensure computation stability
        of the non-linear hillslope processes.

        Note:
            This CFL-like condition is computed using non-linear diffusion coefficients and distances between TIN nodes.

        Args:
            edgelen: numpy array containing the edges of the TIN surface for the considered partition.
            elevation: numpy array containing the edges of the TIN surface for the considered partition.
        """

        # Get the tin edges lengths
        maxCD = max(self.CDaerial,self.CDmarine)

        CFL = numpy.zeros(1)
        if maxCD > 0.:
            Sc = numpy.zeros(1)
            Sc[0] = self.Sc
            mCD = numpy.zeros(1)
            mCD[0] = maxCD
            CFL = sfd.diffnlcfl(Sc, mCD, elev, borders, neighbours, distances, globalIDs)
        else:
            CFL[0] = 1.e6

        # Global mimimum value for diffusion stability
        self.CFL = CFL[0]
        self.updatedt = 0

        return

    def dt_stability_ms(self, edgelen):
        """
        This function computes the maximal timestep to ensure computation stability
        of the river deposited marine sediment processes. This CFL-like condition is computed using diffusion
        coefficients and distances between TIN nodes.

        Note:
            It is worth noticing that the approach does not rely on the elevation of the nodes
            and therefore the maximal hillslope timestep to ensure stability just needs to be
            computed once for each given TIN grid.

        Args:
            edgelen: numpy array containing the edges of the TIN surface for the considered partition.
        """

        # Get the tin edges lengths
        edgedist = edgelen.flatten()
        distIDs = numpy.where(edgedist > 0.)

        # First-order, forward-in-time scheme
        CFL = numpy.zeros(1)
        if self.CDriver == 0.:
            CFL[0] = 0.01*numpy.amin(edgedist[distIDs]**2)
        else:
            CFL[0] = 0.01*numpy.amin(edgedist[distIDs]**2)/self.CDriver

        # Global mimimum value for diffusion stability
        self.CFLms = CFL[0]

        return

    def dt_stability_fail(self, edgelen):
        """
        This function computes the maximal timestep to ensure computation stability
        of the slope failure processes. This CFL-like condition is computed using diffusion
        coefficients and distances between TIN nodes.
        
        Args:
            edgelen: numpy array containing the edges of the TIN surface for the considered partition.
        """

        # Get the tin edges lengths
        edgedist = edgelen.flatten()
        distIDs = numpy.where(edgedist > 0.)

        # First-order, forward-in-time scheme
        CFL = numpy.zeros(1)
        if self.Cfail== 0.:
            CFL[0] = 0.01*numpy.amin(edgedist[distIDs]**2)
        else:
            CFL[0] = 0.01*numpy.amin(edgedist[distIDs]**2)/self.Cfail

        # Global mimimum value for diffusion stability
        self.CFLfail = CFL[0]

        return

    def sedflux(self, sea, elevation, area):
        """
        This function computes the sedimentary fluxes induced by hillslope processes based
        on a linear diffusion approximation.

        Important:
            The linear diffusion process is implemented through the FV approximation and is based on
            the area of each node voronoi polygon and the sum over all the neighbours of the slope of the
            segment (i.e. height differences divided by the length of the mesh edge) as well as the length
            of the corresponding voronoi edge.

        Args:
            sea: float value giving the sea-level height at considered time step.
            elevation: numpy array containing the edges of the TIN surface for the considered partition.
            area: numpy array containing the area of the voronoi polygon for each TIN nodes.

        Returns:
            - CA - numpy array containing the diffusion coefficients.
        """

        if self.ids is None:
            self.ids = numpy.where(area > 0)[0]

        areacoeff = numpy.zeros(len(area))
        areacoeff[self.ids] = 1./area[self.ids]
        coeff = numpy.where(elevation >= sea, self.CDaerial, self.CDmarine)

        return numpy.nan_to_num(areacoeff * coeff)

    def sedfluxmarine(self, sea, elevation, area):
        """
        This function computes the diffusion of marine sediments transported by river processes using
        a linear diffusion approximation.

        Args:
            sea: float value giving the sea-level height at considered time step.
            elevation: numpy array containing the edges of the TIN surface for the considered partition.
            area: numpy array containing the area of the voronoi polygon for each TIN nodes.

        Returns:
            - CA - numpy array containing the diffusion coefficients.
        """

        if self.ids is None:
            self.ids = numpy.where(area > 0)[0]

        areacoeff = numpy.zeros(len(area))
        areacoeff[self.ids] = 1./area[self.ids]
        coeff = numpy.where(elevation >= sea, 0., self.CDriver)

        return numpy.nan_to_num(areacoeff * coeff)

    def sedfluxfailure(self, area):
        """
        This function computes the diffusion of slope failure using a linear diffusion approximation.

        Args:
            area: numpy array containing the area of the voronoi polygon for each TIN nodes.

        Returns:
            - CA - numpy array containing the diffusion coefficients.

        """

        if self.ids is None:
            self.ids = numpy.where(area > 0)[0]

        areacoeff = numpy.zeros(len(area))
        areacoeff[self.ids] = 1./area[self.ids]

        return numpy.nan_to_num(areacoeff * self.Cfail)
