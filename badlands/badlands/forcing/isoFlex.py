##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module is a wrapper for A. Wickert's flexure model (gflex_) and defines several
functions used to compute flexural isostasy response based on **badlands** cumulative
thickness evolution computed over the TIN.

.. image:: img/gflex.png
   :scale: 70 %
   :alt: gflex
   :align: center

Flexural isostasy can be produced in response to a range of geological loads (from Wickert, 2016).


Note:
    Wickert, A. D. (2016), Open-source modular solutions for flexural isostasy: gFlex v1.0,
    Geosci. Model Dev., 9(3), 997–1017, `doi:10.5194/gmd-9-997-2016`_.
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
    This class uses the gFlex model from Wickert to compute flexural isostasy.
    """

    def __init__(self):
        """
        Initialisation.
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
        self.rho_s = 2500.0
        self.rho_w = 1029.0
        self.previous_flex = None
        self.tree = None
        self.Te = None
        self.searchpts = None
        self.Te1 = None
        self.dtime = 0.0
        self.ftime = None

        return

    def buildGrid(
        self,
        nx,
        ny,
        youngMod,
        mantleDensity,
        sedimentDensity,
        elasticT,
        elasticT2,
        Boundaries,
        xyTIN,
        ftime,
    ):
        """
        gFlex initialisation function.

        Args:
            nx: number of points along the X axis
            ny: number of points along the Y axis
            youngMod: Young's modulus
            mantleDensity: mantle density
            sedimentDensity: sediment density
            elasticT: elastic thickness. Can be scalar or an array
            elasticT2: initial elastic thickness.
            Boundaries: list of string describing boundary conditions for West, East, South and North.
            xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN
            ftime: flexure time step
        """
        # Build the flexural grid
        self.nx = nx
        self.ny = ny
        self.xyTIN = xyTIN
        xmin, xmax = min(self.xyTIN[:, 0]), max(self.xyTIN[:, 0])
        ymin, ymax = min(self.xyTIN[:, 1]), max(self.xyTIN[:, 1])
        self.xgrid = numpy.linspace(xmin, xmax, num=self.nx)
        self.ygrid = numpy.linspace(ymin, ymax, num=self.ny)
        self.xi, self.yi = numpy.meshgrid(self.xgrid, self.ygrid)
        self.xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]
        self.ftime = ftime

        # Call gFlex instance
        self.flex = gflex.F2D()
        flex = self.flex

        # Set-up the grid variable
        self.flex.dx = self.xgrid[1] - self.xgrid[0]
        self.flex.dy = self.ygrid[1] - self.ygrid[0]
        self.ball = math.sqrt(
            0.25 * (self.flex.dx * self.flex.dx + self.flex.dy * self.flex.dy)
        )
        tindx = self.xyTIN[1, 0] - self.xyTIN[0, 0]
        self.searchpts = max(int(self.flex.dx * self.flex.dy / (tindx * tindx)), 4)

        # Solution method finite difference
        self.flex.Method = "FD"
        self.flex.Quiet = True

        # van Wees and Cloetingh (1994)
        self.flex.PlateSolutionType = "vWC1994"
        self.flex.Solver = "direct"

        # Acceleration due to gravity
        self.flex.g = 9.8
        # Poisson's Ratio
        self.flex.nu = 0.25
        # Infill Material Density
        self.flex.rho_fill = 0.0
        # Young's Modulus
        self.flex.E = youngMod
        # Mantle Density
        self.flex.rho_m = mantleDensity
        # Sediment Density
        self.rho_s = sedimentDensity
        # Sea Water Density
        self.rho_w = 1029.0
        # Elastic thickness [m]
        if isinstance(elasticT, str):
            TeMap = pandas.read_csv(
                elasticT,
                sep=r"\s+",
                engine="c",
                header=None,
                na_filter=False,
                dtype=numpy.float,
                low_memory=False,
            )
            self.Te = numpy.reshape(TeMap.values, (self.ny, self.nx))
        elif elasticT2 is None:
            self.Te = elasticT * numpy.ones((self.ny, self.nx))
        else:
            self.Te = elasticT2 * numpy.ones((self.ny, self.nx))
            self.Te1 = elasticT

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
        # self.Acell = Acell

        return

    def update_flexure_parameters(self, xyTIN):
        """
        Update TIN variables after 3D displacements.

        Args:
            xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN
        """

        self.xyTIN = xyTIN
        self.tree = cKDTree(self.xyTIN)

        return

    def _compute_flexure(self):
        """
        Use gFlex module to compute flexure from surface load.
        """

        if self.Te1 is None:
            self.flex.Te = self.Te
        else:
            coeff = self.Te1 * numpy.sqrt(self.dtime)
            self.flex.Te = coeff * numpy.ones((self.ny, self.nx)) + self.Te

        self.flex.initialize()

        self.flex.run()

        self.flex.finalize()

        return

    def get_flexure(self, elev, cumdiff, sea, boundsPt, initFlex=False):
        """
        From TIN erosion/deposition values and sea-level compute the
        surface load on the flexural grid.

        Args:
            elev: numpy array containing TIN surface elevation
            cumdiff: numpy array containing TIN cumulative erosion/deposition
            boundsPt: number of points along the boundary
            sea: sea-level.
            initFlex: initialise simulation flexural values

        Returns:
            - flexureTIN - numpy array containing flexural deflection values for the TIN.
        """

        # Average volume of sediment and water on the flexural grid points
        sedload = numpy.zeros(len(self.xyi))
        waterload = numpy.zeros(len(self.xyi))
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)

        if len(elev[indices].shape) == 3:
            elev_vals = elev[indices][:, :, 0]
            cum_vals = cumdiff[indices][:, :, 0]
        else:
            elev_vals = elev[indices]
            cum_vals = cumdiff[indices]

        distances[distances < 0.0001] = 0.0001
        with numpy.errstate(divide="ignore"):
            felev = numpy.average(elev_vals, weights=(1.0 / distances), axis=1)
            fcum = numpy.average(cum_vals, weights=(1.0 / distances), axis=1)

        onIDs = numpy.where(distances[:, 0] <= 0.0001)[0]
        if len(onIDs) > 0:
            felev[onIDs] = elev[indices[onIDs, 0]]
            fcum[onIDs] = cumdiff[indices[onIDs, 0]]
        sedload = fcum
        marine = numpy.where(felev < sea)[0]
        waterload[marine] = sea - felev[marine]

        # Compute surface loads
        self.flex.qs = (
            self.rho_w * self.flex.g * numpy.reshape(waterload, (self.ny, self.nx))
        )
        self.flex.qs += (
            self.rho_s * self.flex.g * numpy.reshape(sedload, (self.ny, self.nx))
        )

        # Compute flexural isostasy with gFlex
        self._compute_flexure()

        # Reinterpolate values on TIN, record new flexural values and compute
        # cumulative flexural values
        if initFlex:
            self.previous_flex = self.flex.w
            flexureTIN = numpy.zeros(len(self.xyTIN[:, 0]))
        else:
            flexureTIN = numpy.zeros(len(self.xyTIN[:, 0]))
            flex_diff = self.flex.w - self.previous_flex
            self.previous_flex = self.flex.w
            rgi_flexure = RegularGridInterpolator((self.ygrid, self.xgrid), flex_diff)
            flexureTIN[boundsPt:] = rgi_flexure(
                (self.xyTIN[boundsPt:, 1], self.xyTIN[boundsPt:, 0])
            )

        self.dtime += self.ftime

        return flexureTIN
