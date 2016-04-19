##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module is a wrapper for A. Wickert's gFlex flexure model and defines several
functions used to compute flexural isostasy response based on Badlands cumulative
thickness evolution computed over the TIN.
"""
import os
import numpy
import gflex
import pandas
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator

class isoFlex:
    """
    This class uses the gFlex model from Wickert 2015 to compute flexural isostasy. For
    more information regarding the code and the available features it is recommended to
    have a look at:
        - gFlex github repository: https://github.com/awickert/gFlex
        - Wickert, A. D.: Open-source modular solutions for flexural isostasy:
            gFlex v1.0, Geosci. Model Dev. Discuss., 8, 4245-4292, 2015
    """
    def __init__(self, nx, ny, youngMod, mantleDensity, sedimentDensity,
                    elasticH, Boundaries, xyTIN, cumflex):
        """
        gFlex initialisation function.

        Parameters
        ----------
        variable : nx, ny
            Number of points along the X and Y axis

        variable : youngMod
            Young's modulus

        variable : mantleDensity
            Mantle density

        variable : sedimentDensity
            Sediment density

        variable : elasticH
            The elastic thickness. Can be scalar or an array

        variable : Boundaries
            List of string describing boundary conditions for West, East, South and North
            Choose from '0Slope0Shear', 'Dirichlet0', '0Displacement0Slope',
                        '0Moment0Shear', 'Periodic', 'Mirror'

        variable : cumflex
            Cumulative flexural displacement.
        """
        # Build the flexural grid
        self.nx = nx
        self.ny = ny

        self.xyTIN = xyTIN
        xmin, xmax = min(self.xyTIN[:,0]), max(self.xyTIN[:,0])
        ymin, ymax = min(self.xyTIN[:,1]), max(self.xyTIN[:,1])
        np.linspace(2.0, 3.0, num=5)
        self.xgrid = numpy.linspace(xmin,xmax,num=self.nx)
        self.ygrid = numpy.linspace(ymin,ymax,num=self.ny)
        self.xi, self.yi = numpy.meshgrid(self.xgrid, self.ygrid)
        self.xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]

        # Call gFlex instance
        self.flex = gflex.F2D()
        flex = self.flex

        # Set-up the grid variable
        self.flex.dx = self.xgrid[1]-self.xgrid[0]
        self.flex.dy = self.ygrid[1]-self.ygrid[0]

        tindx = self.xyTIN[:,1] - self.xyTIN[:,0]

        self.searchpts = max(int(self.flex.dx*self.flex.dy/(tindx*tindx)),4)

        # Solution method finite difference
        self.flex.Method = 'FD'
        self.flex.Quiet = True

        # van Wees and Cloetingh (1994)
        self.flex.PlateSolutionType = 'vWC1994'
        self.flex.Solver = 'direct'

        # Acceleration due to gravity
        self.flex.g = 9.8
        # Poisson's Ratio
        self.flex.nu = 0.25
        # Infill Material Density
        self.flex.rho_fill = 0.
        # Young's Modulus
        self.flex.E = YoungMod
        # Mantle Density
        self.flex.rho_m = mantleDensity
        # Sediment Density
        self.rho_s = sedimentDensity
        # Sea Water Density
        self.rho_w = 1029.0
        # Basement depth
        self.basedepth = 50000.
        # Elastic thickness [m]
        if isinstance(elasticH, basestring):
            TeMap = pandas.read_csv(elasticH, sep=r'\s+', engine='c', header=None,
                na_filter=False, dtype=numpy.float, low_memory=False)
            self.flex.Te = numpy.reshape(TeMap.values, (len(self.nx), len(self.ny)), order='F')
        else:
            self.flex.Te = elasticH * np.ones((self.nx, self.ny))

        # Surface load stresses
        self.flex.qs = numpy.zeros((self.nx, self.ny), dtype=float)

        # Boundary conditions can be:
        self.flex.BC_W = str(Boundaries[0])
        self.flex.BC_E = str(Boundaries[1])
        self.flex.BC_S = str(Boundaries[2])
        self.flex.BC_N = str(Boundaries[3])

        # State of the previous flexural grid used for updating current
        # flexural displacements.
        self.previous_flex = numpy.zeros((self.nx, self.ny), dtype=float)

        self.tree = cKDTree(self.xyTIN)
        self.cumflex = cumflex

        return

    def _compute_flexure(self):
        """
        Use gFlex module to compute flexure from surface load.
        """

        self.flex.initialize()

        self.flex.run()

        self.flex.finalize()

        return

    def get_flexure(self, elev, cumdiff, sea, initFlex=False):
        """
        From TIN erosion/deposition values and sea-level compute the
        surface load on the flexural grid.

        Parameters
        ----------
        variable : elev
            Numpy array containing TIN surface elevation

        variable : cumdiff
            Numpy array containing TIN cumulative erosion/deposition

        variable : sea
            Sea-level.

        variable : initFlex
            Initialise simulation flexural values

        Return
        ----------
        variable : flexureTIN
            Numpy array containing flexural deflection values for the TIN.
        """

        # Inverse distance weighting interpolation to estimate elevation
        # and erosion/deposition values on the flexural grid
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)
        elev_val = elev[indices][:,:,0]
        elevFlex = numpy.average(elev_val,weights=(1./distances), axis=1)
        cum_val = cumdiff[indices][:,:,0]
        sedload = numpy.average(cum_val,weights=(1./distances), axis=1)
        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            elevFlex[onIDs] = elev[indices[onIDs,0]]
            sedload[onIDs] = cumdiff[indices[onIDs,0]]

        # Compute surface loads
        # Get surface load associated with sediment thickness
        qs = self.rho_s * self.flex.g * (self.basedepth + sedload )
        # Add surface load associated to sea water load
        marine = numpy.where(elevFlex < sea)[0]
        qs[marine] += self.rho_w * self.flex.g * (sea - elevFlex[marine])
        # Combine the 2
        self.flex.qs = numpy.reshape(qs,(self.nx, self.ny))

        # Compute flexural isostasy with gFlex
        self._compute_flexure()

        # Reinterpolate values on TIN, record new flexural values and compute
        # cumulative flexural values
        if initFlex:
            self.previous_flex = self.flex.w
            flexureTIN = numpy.zeros(len(self.xyTIN[:,0]))
        else:
            flex_diff = self.flex.w - self.previous_flex
            self.previous_flex = self.flex.w
            rgi_flexure = RegularGridInterpolator((self.xgrid, self.ygrid), flex_diff)
            flexureTIN = rgi_flexure((self.xyTIN[:,0],self.xyTIN[:,1]))
            self.cumflex += flexureTIN

        return flexureTIN
