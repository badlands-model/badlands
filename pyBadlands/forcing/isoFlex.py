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
import math
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
    def __init__(self):
        """
        Initialization.
        """
        self.nx = 0
        self.ny = 0
        self.xyTIN = None
        self.xgrid = None
        self.ygrid = None
        self.xi = None
        self.yi = None
        self.xyi = None
        self.flex = None
        self.ball = 0.
        self.rho_s = 2500.
        self.rho_w = 1029.0
        self.previous_flex = None
        self.tree = None
        #self.Area = None
        #self.Acell = None
        self.Te = None
        self.searchpts = None

        return

    def buildGrid(self, nx, ny, youngMod, mantleDensity, sedimentDensity,
                    elasticT, Boundaries, xyTIN):
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

        variable : elasticT
            The elastic thickness. Can be scalar or an array

        variable : Boundaries
            List of string describing boundary conditions for West, East, South and North
            Choose from '0Slope0Shear', 'Dirichlet0', '0Displacement0Slope',
                        '0Moment0Shear', 'Periodic', 'Mirror'

        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m2)
        """
        # Build the flexural grid
        self.nx = nx
        self.ny = ny
        self.xyTIN = xyTIN
        xmin, xmax = min(self.xyTIN[:,0]), max(self.xyTIN[:,0])
        ymin, ymax = min(self.xyTIN[:,1]), max(self.xyTIN[:,1])
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
        self.ball = math.sqrt(0.25*(self.flex.dx*self.flex.dx + self.flex.dy*self.flex.dy))
        tindx = self.xyTIN[1,0] - self.xyTIN[0,0]
        self.searchpts = max(int(self.flex.dx*self.flex.dy/(tindx*tindx)),4)
        '''
        area = numpy.zeros((self.nx, self.ny))
        area.fill(self.flex.dx*self.flex.dy)
        area[0,:] = 0.5 * self.flex.dx*self.flex.dy
        area[self.nx-1,:] = 0.5 * self.flex.dx*self.flex.dy
        area[:,0] = 0.5 * self.flex.dx*self.flex.dy
        area[:,self.ny-1] = 0.5 * self.flex.dx*self.flex.dy
        area[0,0] = 0.25 * self.flex.dx*self.flex.dy
        area[self.nx-1,0] = area[0,0]
        area[0,self.ny-1] = area[0,0]
        area[self.nx-1,self.ny-1] = area[0,0]
        self.Area = area.flatten()
        '''

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
        self.flex.E = youngMod
        # Mantle Density
        self.flex.rho_m = mantleDensity
        # Sediment Density
        self.rho_s = sedimentDensity
        # Sea Water Density
        self.rho_w = 1029.0
        # Elastic thickness [m]
        if isinstance(elasticT, basestring):
            TeMap = pandas.read_csv(elasticT, sep=r'\s+', engine='c', header=None,
                na_filter=False, dtype=numpy.float, low_memory=False)
            self.Te = numpy.reshape(TeMap.values, (self.ny, self.nx))
        else:
            self.Te = elasticT * numpy.ones((self.ny, self.nx))

        # Surface load stresses
        self.flex.qs = numpy.zeros((self.ny, self.nx), dtype=float)

        # Boundary conditions
        self.flex.BC_W = Boundaries[0]
        self.flex.BC_E = Boundaries[1]
        self.flex.BC_S = Boundaries[2]
        self.flex.BC_N = Boundaries[3]

        # State of the previous flexural grid used for updating current
        # flexural displacements.
        self.previous_flex = numpy.zeros((self.ny, self.nx), dtype=float)

        self.tree = cKDTree(self.xyTIN)
        #self.Acell = Acell

        return

    def update_flexure_parameters(self, xyTIN):
        """
        Update TIN variables after 3D displacements.

        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m2)
        """

        self.xyTIN = xyTIN
        self.tree = cKDTree(self.xyTIN)

        return

    def _compute_flexure(self):
        """
        Use gFlex module to compute flexure from surface load.
        """

        self.flex.Te = self.Te

        self.flex.initialize()

        self.flex.run()

        self.flex.finalize()

        return

    def get_flexure(self, elev, cumdiff, sea, boundsPt, initFlex=False):
        """
        From TIN erosion/deposition values and sea-level compute the
        surface load on the flexural grid.

        Parameters
        ----------
        variable : elev
            Numpy array containing TIN surface elevation

        variable : cumdiff
            Numpy array containing TIN cumulative erosion/deposition

        variable : boundsPt
            Number of points along the boundary

        variable : sea
            Sea-level.

        variable : initFlex
            Initialise simulation flexural values

        Return
        ----------
        variable : flexureTIN
            Numpy array containing flexural deflection values for the TIN.
        """

        # Average volume of sediment and water on the flexural grid points
        ballIDs = self.tree.query_ball_point(self.xyi, self.ball)
        sedload = numpy.zeros(len(self.xyi))
        waterload = numpy.zeros(len(self.xyi))
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)
        elev_vals = elev[indices]
        felev = numpy.average(elev_vals,weights=(1./distances), axis=1)
        cum_vals = cumdiff[indices]
        fcum = numpy.average(cum_vals,weights=(1./distances), axis=1)
        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            felev[onIDs] = elev[indices[onIDs,0]]
            fcum[onIDs] = cumdiff[indices[onIDs,0]]
        sedload = fcum
        marine = numpy.where(felev < sea)[0]
        waterload[marine] = sea - felev[marine]
        '''
        for i in range(len(self.xyi)):
            ids = numpy.asarray(ballIDs[i], dtype=numpy.int)
            inIDX = numpy.where(numpy.logical_and(self.xyTIN[ids,0] >= self.xyi[i,0] - 0.5*self.flex.dx,
                                                  self.xyTIN[ids,0] <= self.xyi[i,0] + 0.5*self.flex.dx))[0]
            inIDY = numpy.where(numpy.logical_and(self.xyTIN[ids,1] >= self.xyi[i,1] - 0.5*self.flex.dy,
                                                  self.xyTIN[ids,1] <= self.xyi[i,1] + 0.5*self.flex.dy))[0]
            inIDs = numpy.intersect1d(ids[inIDX], ids[inIDY])
            sedload[i] = numpy.sum(cumdiff[inIDs]*self.Acell[inIDs]) / self.Area[i]

            marine = numpy.where(elev[inIDs] < sea)[0]

            waterload[i] = numpy.sum((sea-elev[inIDs[marine]])*self.Acell[inIDs[marine]]) / self.Area[i]
        '''
        # Compute surface loads
        self.flex.qs = self.rho_w * self.flex.g * numpy.reshape(waterload,(self.ny, self.nx))
        self.flex.qs += self.rho_s * self.flex.g * (self.Te + numpy.reshape(sedload,(self.ny, self.nx)) )
        # Compute flexural isostasy with gFlex
        self._compute_flexure()

        # Reinterpolate values on TIN, record new flexural values and compute
        # cumulative flexural values
        if initFlex:
            self.previous_flex = self.flex.w
            flexureTIN = numpy.zeros(len(self.xyTIN[:,0]))
        else:
            flexureTIN = numpy.zeros(len(self.xyTIN[:,0]))
            flex_diff = self.flex.w - self.previous_flex
            self.previous_flex = self.flex.w
            rgi_flexure = RegularGridInterpolator((self.ygrid, self.xgrid), flex_diff)
            flexureTIN[boundsPt:] = rgi_flexure((self.xyTIN[boundsPt:,1],self.xyTIN[boundsPt:,0]))

        return flexureTIN
